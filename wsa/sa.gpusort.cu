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

// -----------------------------------------------------------------------
// GPUIndexType: shared genome index plus persistent output buffers.
//
// Under mutex protection a single caller uses this struct at a time, so
// the output buffers (device and pinned-host) can be preallocated once in
// GPUIndexCreate and reused across every block.  They grow by doubling
// when a block produces more matches than the current capacity, but they
// never shrink — this minimises cudaMallocHost / cudaFreeHost calls to at
// most O(log N_max) over the entire run rather than once per block.
// -----------------------------------------------------------------------
struct GPUIndexType
{
    std::vector<IndexPartition>  partitions;   // one per NN sub-table, read-only after creation

    // Persistent device output buffer, reused every block.
    thrust::device_vector<SEEDMATCH>  out_pairs;

    // Persistent pinned host buffer, reused every block.
    // Allocated with cudaMallocHost so device→host DMA is direct.
    SEEDMATCH*   host_buf;
    std::size_t  host_buf_capacity;  // in SEEDMATCH records
};


GPUIndex* GPUIndexCreate(CW** index_parts, long int* sizes, unsigned int num_parts)
{
    GPUIndexType* index = new GPUIndexType();
    index->partitions.resize(num_parts);

    for (unsigned int i = 0; i < num_parts; i++) {
        IndexPartition& p = index->partitions[i];

        auto t0 = std::chrono::high_resolution_clock::now();
        p.cws.assign(index_parts[i], index_parts[i] + sizes[i]);
        auto t1 = std::chrono::high_resolution_clock::now();
        std::cerr << "Copy index partition " << i
                  << " to GPU (" << sizes[i] << "): "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
                  << "ms" << std::endl;

        auto t2 = std::chrono::high_resolution_clock::now();
        BuildRuns(p.cws, p.unique_seeds, p.counts, p.starts);
        auto t3 = std::chrono::high_resolution_clock::now();
        std::cerr << "BuildRuns index partition " << i
                  << " (" << p.unique_seeds.size() << " unique seeds): "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()
                  << "ms" << std::endl;
    }

    // Preallocate output buffers.  1<<20 (~1M) records is a comfortable
    // starting capacity; doubling on demand costs at most a handful of
    // reallocations over an entire multi-gigabase run.
    const std::size_t initial = 1 << 20;
    index->out_pairs.resize(initial);

    index->host_buf_capacity = initial;
    cudaError_t err = cudaMallocHost(reinterpret_cast<void**>(&index->host_buf),
                                     initial * sizeof(SEEDMATCH));
    if (err != cudaSuccess) {
        std::cerr << "GPUIndexCreate: cudaMallocHost failed: "
                  << cudaGetErrorString(err) << std::endl;
        index->host_buf = nullptr;
        index->host_buf_capacity = 0;
    }

    return index;
}


GPUIndex* GPUIndexFree(GPUIndex* idx)
{
    GPUIndexType* idxobj = static_cast<GPUIndexType*>(idx);
    if (idxobj) {
        if (idxobj->host_buf)
            cudaFreeHost(idxobj->host_buf);
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

// saGPUFreeHostBuffer: exposed so sa.main.c can explicitly release the
// pinned SEEDMATCH buffer before the standard bigArray destructor runs,
// keeping cudaFreeHost out of the bigArray library.
void saGPUFreeHostBuffer(void* ptr)
{
    if (ptr)
        cudaFreeHost(ptr);
}

// ---------------------------------------------------------------------------
// saGPUMatchHits
//
// Called under pthread_mutex_lock by sa.main.c.  Because only one agent
// enters at a time the persistent output buffers in GPUIndexType are safe
// to reuse.  Both buffers grow by doubling when needed but never shrink,
// so the expensive cudaMallocHost / cudaFreeHost occur at most O(log N_max)
// times over the entire run.
//
// Returns N (number of SEEDMATCH records found).
// *out_buffer receives a pointer into the persistent pinned host buffer;
// the caller must copy or consume the data before releasing the mutex and
// before the next saGPUMatchHits call overwrites the buffer.
// In practice sa.main.c does:
//   bb.sms = bigArrayHandleCreate(N+1, SEEDMATCH, bb.h);
//   bigArrayMax(bb.sms) = N;
//   memcpy(bigArrayp(bb.sms,0,SEEDMATCH), *out_buffer, N*sizeof(SEEDMATCH));
// and then releases the mutex, so the pinned buffer is never aliased.
// ---------------------------------------------------------------------------
unsigned int saGPUMatchHits(GPUIndex* idx, CW** words, long int* sizes,
                            unsigned int num_parts, SEEDMATCH** out_buffer)
{
    auto start = std::chrono::high_resolution_clock::now();
    GPUIndexType* index = static_cast<GPUIndexType*>(idx);
    *out_buffer = nullptr;

    std::size_t last_pair = 0;

    for (unsigned int i = 0; i < num_parts; i++) {
        const IndexPartition& gpart = index->partitions[i];

        thrust::device_vector<CW> word_vec(words[i], words[i] + sizes[i]);
        thrust::sort(word_vec.begin(), word_vec.end(), compare_CW());

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

        // Grow persistent device buffer if needed (doubling, never shrinks).
        while (last_pair + total_out > index->out_pairs.size()) {
            std::size_t new_size = index->out_pairs.size() * 2;
            index->out_pairs.resize(new_size);
            std::cerr << "Resize device buffer\t" << new_size << std::endl;
        }

        const CW*            d_idx        = thrust::raw_pointer_cast(gpart.cws.data());
        const CW*            d_w          = thrust::raw_pointer_cast(word_vec.data());
        const std::uint32_t* d_idx_starts = thrust::raw_pointer_cast(gpart.starts.data());
        const std::uint32_t* d_idx_counts = thrust::raw_pointer_cast(gpart.counts.data());
        const std::uint32_t* d_w_starts   = thrust::raw_pointer_cast(w_starts.data());
        const std::uint32_t* d_w_counts   = thrust::raw_pointer_cast(w_counts.data());
        const std::uint32_t* d_idx_common = thrust::raw_pointer_cast(idx_common.data());
        const std::uint32_t* d_w_common   = thrust::raw_pointer_cast(w_common.data());
        const std::uint32_t* d_out_offsets= thrust::raw_pointer_cast(out_offsets.data());
        SEEDMATCH*           d_pairs      = thrust::raw_pointer_cast(
                                               index->out_pairs.data() + last_pair);

        EmitCartesianProduct<<<static_cast<int>(num_common), 256>>>(
            d_idx, d_idx_starts, d_idx_counts,
            d_w,   d_w_starts,   d_w_counts,
            d_idx_common, d_w_common,
            d_out_offsets, d_pairs, num_common);

        last_pair += total_out;
    }

    // Trim logical size (does not free memory, just adjusts end iterator).
    index->out_pairs.resize(last_pair);

    // Sort matches by read id.
    thrust::sort(index->out_pairs.begin(), index->out_pairs.end(),
                 [] __device__ (const SEEDMATCH& a, const SEEDMATCH& b)
                 { return a.read < b.read; });

    std::size_t N = last_pair;

    auto end = std::chrono::high_resolution_clock::now();
    std::cerr << "Found " << N << " matches in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << "ms" << std::endl;

    if (N == 0) {
        return 0;
    }

    // Grow persistent pinned host buffer if needed (doubling, never shrinks).
    if (N > index->host_buf_capacity) {
        std::size_t new_cap = index->host_buf_capacity ? index->host_buf_capacity : 1;
        while (new_cap < N) new_cap *= 2;
        if (index->host_buf)
            cudaFreeHost(index->host_buf);
        cudaError_t err = cudaMallocHost(reinterpret_cast<void**>(&index->host_buf),
                                         new_cap * sizeof(SEEDMATCH));
        if (err != cudaSuccess) {
            std::cerr << "saGPUMatchHits: cudaMallocHost failed for " << new_cap
                      << " records: " << cudaGetErrorString(err) << std::endl;
            index->host_buf = nullptr;
            index->host_buf_capacity = 0;
            return 0;
        }
        index->host_buf_capacity = new_cap;
        std::cerr << "Resize pinned host buffer\t" << new_cap << std::endl;
    }

    // Direct DMA from device to pinned host — no staging, no extra copy.
    thrust::copy_n(index->out_pairs.begin(), N, index->host_buf);

    // Restore full capacity for next block (resize to capacity, not to N,
    // so the next block starts with the full pre-allocated device memory).
    index->out_pairs.resize(index->out_pairs.capacity());

    *out_buffer = index->host_buf;
    return static_cast<unsigned int>(N);
}
