/* sa.gpusort.cu
 * Code for sorting on a GPU using Nvidia Thrust library
 *
 * This module is part of the sortalign package
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 *
 * Created: December 30, 2025
 *
 * Optimisation 1 (April 2026) — precomputed genome BuildRuns:
 *   The run-length encoding (unique seeds, counts, starts) of the genome
 *   index partitions is now computed once inside GPUIndexCreate and stored
 *   in device memory for the lifetime of the index.  Previously BuildRuns
 *   was called on every genome partition for every read block, paying the
 *   cost of a reduce_by_key + exclusive_scan over the full 3 GB genome
 *   index on every block.  The genome index never changes between blocks,
 *   so this computation was entirely redundant.
 *
 * Optimisation 2 (April 2026) — per-call output buffer, no mutex:
 *   saGPUMatchHits now allocates its SEEDMATCH output buffer locally as a
 *   thrust::device_vector, then copies the result into a cudaMallocHost
 *   pinned buffer before returning.  The pinned buffer is returned to the
 *   caller via SEEDMATCH** out_buffer.  Because all writable device state
 *   is local to each call, concurrent agents require no mutex.
 *   The caller frees the pinned buffer with cudaFreeHost (called explicitly
 *   before the standard bigArray destructor so that the CUDA free is
 *   localised to sa.main.c and the bigArray library needs no USEGPU ifdefs).
 *
 * Optimisation 3 (April 2026) — pinned CW upload:
 *   The caller (sa.main.c) calls cudaHostRegister on each cwsN bigArray
 *   base pointer before calling saGPUMatchHits, and cudaHostUnregister
 *   afterwards.  This pins the existing posix_memalign CW buffers in place
 *   so the PCIe upload inside saGPUMatchHits runs at full bandwidth with
 *   no copy and no change to the bigArray allocator.
 */

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/set_operations.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/adjacent_difference.h>
#include <thrust/gather.h>
#include <chrono>
#include <iostream>

#include "sa.gpusort.h"
#include "sa.common.h"


// Comparators for thrust::sort(), reimplemented versions of comparators in
// sa.sort.c so that the compiler can better optimize the code.
struct compare_CW {
    // compare code words, same as cwOrder in sa.sort.c
    __host__ __device__
    bool operator()(const CW& a, const CW& b) const
    {return a.seed <= b.seed ; }  // we do not need a more detailed comparisons

// {return a.seed < b.seed || (a.seed == b.seed && a.nam < b.nam) ||
//    (a.seed == b.seed && a.nam == b.nam && a.pos < b.pos);}
};


struct compare_HIT {
    // compare hits, same as hitOrder in sa.sort.c
    __host__ __device__
    bool operator()(const HIT& a, const HIT& b) const
    {return a.read < b.read || a.read == b.read && a.chrom < b.chrom ||
    a.read == b.read && a.chrom == b.chrom && a.a1 < b.a1 ||
    a.read == b.read && a.chrom == b.chrom && a.a1 == b.a1 && a.x1 < b.x1;}
};


struct compare_HIT_pairs {
    // compare hits for read pairs, same as hitPairOrder in sa.sort.c
    __host__ __device__
    bool operator()(const HIT& a, const HIT& b) const
    {
        if ((a.read >> 1) < (b.read >> 1)) {
            return true;
        }
        if (a.read == b.read) {
            if (a.chrom < b.chrom) {
                return true;
            }
            int n1 = a.a1 + (a.x1 >> NSHIFTEDTARGETREPEATBITS);
            int n2 = b.a1 + (b.x1 >> NSHIFTEDTARGETREPEATBITS);
            if (n1 < n2) {
                return true;
            }
            return n1 == n2 && a.x1 < b.x1;
        }
        return false;
    }
};

// sort on a GPU
template<typename T, typename CMP>
void saGPUSort(T* cp, long int number_of_records)
{
    auto start = std::chrono::high_resolution_clock::now();
    // copy data to a GPU
    thrust::device_vector<T> d_vec(cp, cp + number_of_records);
    auto end = std::chrono::high_resolution_clock::now();
    std::cerr << "Copy data to GPU: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    // sort
    thrust::sort(d_vec.begin(), d_vec.end(), CMP());
    end = std::chrono::high_resolution_clock::now();
    // std::cerr << "Sort: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    // copy sorted data back to the host
    thrust::copy(d_vec.begin(), d_vec.end(), cp);
    end = std::chrono::high_resolution_clock::now();
    std::cerr << "Copy sorted data to back: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
}

// the sort function callable from C
void saGPUSort (char *cp, long int number_of_records, int type)
{

    switch (type) {
        case 1:
            saGPUSort<CW, compare_CW>(reinterpret_cast<CW*>(cp), number_of_records);
            break;

        case 2:
            saGPUSort<HIT, compare_HIT>(reinterpret_cast<HIT*>(cp), number_of_records);
            break;

        case 3:
            saGPUSort<HIT, compare_HIT_pairs>(reinterpret_cast<HIT*>(cp), number_of_records);
            break;
    };

    return ;
}


struct SeedOfCW {
    __host__ __device__
    std::uint32_t operator()(const CW& x) const {
        return x.seed;
    }
};

// Find runs of the same CW::seed value in input and return start position
// and counts.  Input must be sorted by CW::seed.
static void BuildRuns(const thrust::device_vector<CW>& input,
                      thrust::device_vector<std::uint32_t>& unique_seeds,
                      thrust::device_vector<std::uint32_t>& counts,
                      thrust::device_vector<std::uint32_t>& starts)
{
    unique_seeds.resize(input.size());
    counts.resize(input.size());

    auto keys_begin = thrust::make_transform_iterator(input.begin(), SeedOfCW{});
    auto keys_end   = keys_begin + input.size();

    auto end_pair = thrust::reduce_by_key(
        keys_begin, keys_end,
        thrust::make_constant_iterator<std::uint32_t>(1U),
        unique_seeds.begin(),
        counts.begin()
    );

    std::size_t n_groups = static_cast<std::size_t>(end_pair.first - unique_seeds.begin());
    unique_seeds.resize(n_groups);
    counts.resize(n_groups);
    starts.resize(n_groups);

    // find prefix sums
    thrust::exclusive_scan(counts.begin(), counts.end(), starts.begin(), std::uint32_t{0});
}


// -----------------------------------------------------------------------
// Per-partition storage: CW data plus its precomputed run-length encoding.
// The run-length encoding is computed once in GPUIndexCreate and reused
// for every read block that passes through saGPUMatchHits.
// -----------------------------------------------------------------------
struct IndexPartition {
    thrust::device_vector<CW>            cws;           // raw sorted CW records on device
    thrust::device_vector<std::uint32_t> unique_seeds;  // one entry per run
    thrust::device_vector<std::uint32_t> counts;        // length of each run
    thrust::device_vector<std::uint32_t> starts;        // start offset of each run in cws
};

// Data used to search genome index.
// out_pairs has been removed: each saGPUMatchHits call allocates its own
// local output buffer, so concurrent calls from different agents do not
// share any writable device memory and need no mutex.
struct GPUIndexType
{
    std::vector<IndexPartition>  partitions;  // one per NN sub-table, read-only after creation
};


GPUIndex* GPUIndexCreate(CW** index_parts, long int* sizes, unsigned int num_parts)
{
    GPUIndexType* index = new GPUIndexType();
    index->partitions.resize(num_parts);

    for (unsigned int i = 0; i < num_parts; i++) {
        IndexPartition& p = index->partitions[i];

        // --- upload partition to device ---
        auto t0 = std::chrono::high_resolution_clock::now();
        p.cws.assign(index_parts[i], index_parts[i] + sizes[i]);
        auto t1 = std::chrono::high_resolution_clock::now();
        std::cerr << "Copy index partition " << i
                  << " to GPU (" << sizes[i] << "): "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
                  << "ms" << std::endl;

        // --- precompute run-length encoding for this partition ---
        // This replaces the per-block BuildRuns call that the original code
        // performed on every read block.  The genome index is read-only after
        // this point, so paying the cost once here is correct.
        auto t2 = std::chrono::high_resolution_clock::now();
        BuildRuns(p.cws, p.unique_seeds, p.counts, p.starts);
        auto t3 = std::chrono::high_resolution_clock::now();
        std::cerr << "BuildRuns index partition " << i
                  << " (" << p.unique_seeds.size() << " unique seeds): "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()
                  << "ms" << std::endl;
    }

    return index;
}


GPUIndex* GPUIndexFree(GPUIndex* idx)
{
    GPUIndexType* idxobj = static_cast<GPUIndexType*>(idx);
    if (idxobj) {
        delete idxobj;
    }
    return nullptr;
}


__global__
void EmitCartesianProduct(const CW* index,
                          const std::uint32_t* idx_starts,
                          const std::uint32_t* idx_counts,
                          const CW* words,
                          const std::uint32_t* w_starts,
                          const std::uint32_t* w_counts,
                          const std::uint32_t* idx_common,
                          const std::uint32_t* w_common,
                          const std::uint32_t* out_offsets,
                          SEEDMATCH* out_pairs,
                          std::size_t num_common)
{
    std::size_t ind = static_cast<std::size_t>(blockIdx.x);
    if (ind >= num_common) {
        return;
    }

    std::uint32_t ind_idx = idx_common[ind];
    std::uint32_t ind_w = w_common[ind];

    std::uint32_t start_idx = idx_starts[ind_idx];
    std::uint32_t start_w = w_starts[ind_w];
    std::uint32_t num_idx = idx_counts[ind_idx];
    std::uint32_t num_w = w_counts[ind_w];

    std::size_t total = num_idx * num_w;
    std::size_t start = out_offsets[ind];

    for (std::size_t i=static_cast<std::size_t>(threadIdx.x);i < total;
         i += static_cast<std::size_t>(blockDim.x)) {

        std::size_t k_w = i / num_idx;
        std::size_t k_idx = i % num_idx;

        const CW& idx = index[start_idx + k_idx];
        const CW& w = words[start_w + k_w];

        out_pairs[start + i] = SEEDMATCH {w.seed, w.nam, w.pos, w.intron,
                                          idx.seed, idx.nam, idx.pos, idx.intron};
    }
}


struct diff_CW {
    __host__ __device__
    unsigned int operator()(const CW& a, const CW& b) {return b.seed - a.seed;}
};

// ---------------------------------------------------------------------------
// Pinned host memory allocation for SEEDMATCH and CW arrays.
//
// cudaMallocHost allocates page-locked (pinned) memory on the host.  Pinned
// memory has two advantages over ordinary malloc memory in this pipeline:
//
//   1. PCIe transfers to/from pinned buffers bypass the driver's internal
//      bounce buffer and run at full PCIe bandwidth — typically 2× faster
//      than transfers from pageable memory.
//
//   2. The GPU can DMA directly into the pinned SEEDMATCH buffer during the
//      final thrust::copy, so no extra device→host staging is needed.
//
// cudaMallocHost guarantees at least 256-byte alignment, which satisfies
// both SEEDMATCH (__attribute__((aligned(32)))) and CW (__attribute__((aligned(16)))).
//
// IMPORTANT: pinned memory MUST be released with cudaFreeHost, NOT free().
// The bigArray caller must use bigArraySwitchCudaBase so that bigArray's
// destructor calls cudaFreeHost rather than free().
//
// saGPUAllocHostCW and saGPUAllocHostSeedMatches are separate entry points
// so that the caller can apply bigArraySwitchCudaBase with the correct
// element size in each case.
// ---------------------------------------------------------------------------

// saGPUFreeHostBuffer: release a cudaMallocHost buffer.
// Called explicitly by sa.main.c before the standard bigArray destructor,
// so that cudaFreeHost is localised to the magic2 GPU layer and the bigArray
// library needs no USEGPU ifdefs.
void saGPUFreeHostBuffer(void* ptr)
{
    if (ptr)
        cudaFreeHost(ptr);
}

// ---------------------------------------------------------------------------
// saGPUMatchHits
//
// Find all SEEDMATCH records between the current read block and the (resident)
// genome index, sort by read id, and return both the count and a pointer to a
// freshly allocated pinned host buffer containing the results.
//
// The caller owns the returned buffer and must free it with
// saGPUFreeHostBuffer() (via bigArraySwitchCudaBase on destruction).
//
// All writable device memory is allocated locally within this call.
// The genome index partitions (GPUIndexType::partitions) are read-only.
// Concurrent calls from different agents therefore share no writable state
// and require NO MUTEX.
//
// Workflow:
//   1. For each of the NN seed partitions:
//      a. Upload read CW array to device (fast if source is pinned via
//         saGPUAllocHostCW / bigArraySwitchCudaBase).
//      b. Sort by seed.
//      c. BuildRuns on reads (genome runs are precomputed).
//      d. set_intersection to find common seeds.
//      e. gather counts, exclusive_scan for output offsets.
//      f. EmitCartesianProduct kernel writes SEEDMATCHes into local
//         device buffer, growing it as needed.
//   2. Sort all SEEDMATCHes by read id.
//   3. cudaMallocHost a pinned host buffer of exactly N records.
//   4. thrust::copy the sorted device results into the pinned buffer.
//      Because the destination is pinned, this DMA is direct and fast.
//   5. Return N and the pinned pointer via *out_buffer.
// ---------------------------------------------------------------------------
unsigned int saGPUMatchHits(GPUIndex* idx, CW** words, long int* sizes,
                            unsigned int num_parts, SEEDMATCH** out_buffer)
{
    auto start = std::chrono::high_resolution_clock::now();
    GPUIndexType* index = static_cast<GPUIndexType*>(idx);
    *out_buffer = nullptr;

    // Local device output buffer — not shared between concurrent callers.
    thrust::device_vector<SEEDMATCH> out_pairs(1 << 20);
    std::size_t last_pair = 0;

    for (unsigned int i = 0; i < num_parts; i++) {
        // Genome partition: read-only, precomputed in GPUIndexCreate.
        const IndexPartition& gpart = index->partitions[i];

        // Upload read CW array.  If words[i] was allocated via
        // saGPUAllocHostCW (pinned), this transfer runs at full PCIe
        // bandwidth.  If it is ordinary pageable memory it still works,
        // just slightly slower.
        thrust::device_vector<CW> word_vec(words[i], words[i] + sizes[i]);
        thrust::sort(word_vec.begin(), word_vec.end(), compare_CW());

        // Build runs for reads (genome runs already built).
        thrust::device_vector<std::uint32_t> w_unique_seeds;
        thrust::device_vector<std::uint32_t> w_counts;
        thrust::device_vector<std::uint32_t> w_starts;
        BuildRuns(word_vec, w_unique_seeds, w_counts, w_starts);

        std::size_t max_common_words = std::min(gpart.unique_seeds.size(),
                                                w_unique_seeds.size());
        if (max_common_words == 0)
            continue;

        thrust::device_vector<std::uint32_t> common_words(max_common_words);
        auto common_end = thrust::set_intersection(w_unique_seeds.begin(),
                                                   w_unique_seeds.end(),
                                                   gpart.unique_seeds.begin(),
                                                   gpart.unique_seeds.end(),
                                                   common_words.begin());

        std::size_t num_common = static_cast<std::size_t>(
                                     common_end - common_words.begin());
        common_words.resize(num_common);
        if (num_common == 0)
            continue;

        thrust::device_vector<std::uint32_t> idx_common(num_common);
        thrust::device_vector<std::uint32_t> w_common(num_common);

        thrust::lower_bound(gpart.unique_seeds.begin(), gpart.unique_seeds.end(),
                            common_words.begin(), common_words.end(),
                            idx_common.begin());
        thrust::lower_bound(w_unique_seeds.begin(), w_unique_seeds.end(),
                            common_words.begin(), common_words.end(),
                            w_common.begin());

        thrust::device_vector<std::uint32_t> idx_matched_counts(num_common);
        thrust::device_vector<std::uint32_t> w_matched_counts(num_common);
        thrust::gather(idx_common.begin(), idx_common.end(),
                       gpart.counts.begin(), idx_matched_counts.begin());
        thrust::gather(w_common.begin(), w_common.end(),
                       w_counts.begin(), w_matched_counts.begin());

        thrust::device_vector<std::uint32_t> out_counts(num_common);
        thrust::transform(idx_matched_counts.begin(), idx_matched_counts.end(),
                          w_matched_counts.begin(),
                          out_counts.begin(),
                          thrust::multiplies<std::uint32_t>());

        thrust::device_vector<std::uint32_t> out_offsets(num_common);
        thrust::exclusive_scan(out_counts.begin(), out_counts.end(),
                               out_offsets.begin(), std::uint32_t{0});

        std::uint32_t last_offset = 0, last_count = 0;
        thrust::copy_n(out_offsets.begin() + (num_common - 1), 1, &last_offset);
        thrust::copy_n(out_counts.begin()  + (num_common - 1), 1, &last_count);
        std::size_t total_out = last_offset + last_count;

        const CW*           d_idx        = thrust::raw_pointer_cast(gpart.cws.data());
        const CW*           d_w          = thrust::raw_pointer_cast(word_vec.data());
        const std::uint32_t* d_idx_starts = thrust::raw_pointer_cast(gpart.starts.data());
        const std::uint32_t* d_idx_counts = thrust::raw_pointer_cast(gpart.counts.data());
        const std::uint32_t* d_w_starts   = thrust::raw_pointer_cast(w_starts.data());
        const std::uint32_t* d_w_counts   = thrust::raw_pointer_cast(w_counts.data());
        const std::uint32_t* d_idx_common = thrust::raw_pointer_cast(idx_common.data());
        const std::uint32_t* d_w_common   = thrust::raw_pointer_cast(w_common.data());
        const std::uint32_t* d_out_offsets= thrust::raw_pointer_cast(out_offsets.data());

        while (last_pair + total_out > out_pairs.size()) {
            std::size_t new_size = out_pairs.size() * 2;
            out_pairs.resize(new_size);
            std::cerr << "Resize\t" << new_size << std::endl;
        }
        SEEDMATCH* d_pairs = thrust::raw_pointer_cast(out_pairs.data() + last_pair);

        EmitCartesianProduct<<<static_cast<int>(num_common), 256>>>(
            d_idx, d_idx_starts, d_idx_counts,
            d_w,   d_w_starts,   d_w_counts,
            d_idx_common, d_w_common,
            d_out_offsets, d_pairs, num_common);

        last_pair += total_out;
    }

    out_pairs.resize(last_pair);

    // Sort all matches by read id.
    thrust::sort(out_pairs.begin(), out_pairs.end(),
                 [] __device__ (const SEEDMATCH& a, const SEEDMATCH& b)
                 { return a.read < b.read; });

    auto end = std::chrono::high_resolution_clock::now();
    std::size_t N = out_pairs.size();
    std::cerr << "Found " << N << " matches in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << "ms" << std::endl;

    if (N == 0) {
        *out_buffer = nullptr;
        return 0;
    }

    // Allocate a pinned host buffer of exactly N records.
    // cudaMallocHost guarantees >= 256-byte alignment, satisfying
    // SEEDMATCH's __attribute__((aligned(32))).
    // The caller must free this with saGPUFreeHostBuffer() (bigArray
    // destructor will do so automatically via bigArraySwitchCudaBase).
    SEEDMATCH* pinned = nullptr;
    cudaError_t err = cudaMallocHost(reinterpret_cast<void**>(&pinned),
                                     N * sizeof(SEEDMATCH));
    if (err != cudaSuccess) {
        std::cerr << "saGPUMatchHits: cudaMallocHost failed for " << N
                  << " SEEDMATCH records: " << cudaGetErrorString(err) << std::endl;
        // fall back: caller gets nullptr and 0 — treat as empty block
        return 0;
    }

    // Direct DMA from device to pinned host — no staging, no extra copy.
    thrust::copy_n(out_pairs.begin(), N, pinned);

    *out_buffer = pinned;
    // out_pairs goes out of scope here; device memory freed automatically.
    return static_cast<unsigned int>(N);
}
