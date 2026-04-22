/* sa.gpusort.cu
 * Code for sorting on a GPU using Nvidia Thrust library
 *
 * This module is part of the sortalign package
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 *
 * Created: December 30, 2025
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


// Data used to search genome index
struct GPUIndexType
{
    std::vector< thrust::device_vector<CW> > d_vecs;
};

struct SeedOfCW {
    __host__ __device__
    std::uint32_t operator()(const CW& x) const {
        return x.seed;
    }
};

// Find runs of the the same CW::seed value in input and return start position
// and counts. Input must be sorted by CW::seed.
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
    thrust::exclusive_scan(counts.begin(), counts.end(), starts.begin(), std::uint64_t{0});
}


GPUIndex* GPUIndexCreate(CW** index_parts, long int* sizes, unsigned int num_parts)
{
    GPUIndexType* index = new GPUIndexType();
    index->d_vecs.reserve(num_parts);
    for (unsigned int i=0;i < num_parts;i++) {
        auto start = std::chrono::high_resolution_clock::now();
        index->d_vecs.emplace_back(index_parts[i], index_parts[i] + sizes[i]);
        auto end = std::chrono::high_resolution_clock::now();
        std::cerr << "Copy index partition " << i << " to GPU (" << sizes[i] << "): "  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
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

        out_pairs[start + i] = SEEDMATCH {0, w.nam, w.pos, w.intron,
                                          0, idx.nam, idx.pos, idx.intron};
    }
}


struct diff_CW {
    __host__ __device__
    unsigned int operator()(const CW& a, const CW& b) {return b.seed - a.seed;}
};

void saGPUMatchHits(GPUIndex* idx, CW** words, long int* sizes,
                    unsigned int num_parts)
{
    GPUIndexType* index = static_cast<GPUIndexType*>(idx);
    // maximum size of the output vector
    std::size_t output_size = 0;
    for (unsigned int i=0;i < num_parts;i++) {
        output_size += min(index->d_vecs[i].size(), sizes[i]);
    }
    thrust::device_vector<SEEDMATCH> out_pairs(output_size);
    // points to the last output match saved so far
    std::size_t last_pair = 0;

    for (unsigned int i=0;i < num_parts;i++) {
        auto start = std::chrono::high_resolution_clock::now();
        thrust::device_vector<CW> word_vec(words[i], words[i] + sizes[i]);
        auto end = std::chrono::high_resolution_clock::now();
        std::cerr << "Copied word partition " << i << " to GPU (" << sizes[i] << "): "  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

        start = std::chrono::high_resolution_clock::now();
        thrust::sort(word_vec.begin(), word_vec.end(), compare_CW());
        end = std::chrono::high_resolution_clock::now();
        std::cerr << "Sorted word partition " << i << " in "  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

        // build runs of unique words for index and reads
        thrust::device_vector<std::uint32_t> idx_unique_seeds;
        thrust::device_vector<std::uint32_t> idx_counts;
        thrust::device_vector<std::uint32_t> idx_starts;
        start = std::chrono::high_resolution_clock::now();
        BuildRuns(index->d_vecs[i], idx_unique_seeds, idx_counts, idx_starts);
        end = std::chrono::high_resolution_clock::now();
        std::cerr << "Build runs for index partition " << i << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << ", " << idx_unique_seeds.size() << " unique seeds" << std::endl;

        thrust::device_vector<std::uint32_t> w_unique_seeds;
        thrust::device_vector<std::uint32_t> w_counts;
        thrust::device_vector<std::uint32_t> w_starts;
        start = std::chrono::high_resolution_clock::now();
        BuildRuns(word_vec, w_unique_seeds, w_counts, w_starts);
        end = std::chrono::high_resolution_clock::now();
        std::cerr << "Build runs for read partition " << i << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << ", " << w_unique_seeds.size() << " unique seeds" << std::endl;


        // find common words
        std::size_t max_common_words = std::min(idx_unique_seeds.size(),
                                                w_unique_seeds.size());
        if (max_common_words == 0) {
            continue;
        }


        start = std::chrono::high_resolution_clock::now();
        thrust::device_vector<std::uint32_t> common_words(max_common_words);
        auto common_end = thrust::set_intersection(w_unique_seeds.begin(),
                                                   w_unique_seeds.end(),
                                                   idx_unique_seeds.begin(),
                                                   idx_unique_seeds.end(),
                                                   common_words.begin());

        std::size_t num_common = static_cast<std::size_t>(
                                        common_end - common_words.begin());
        common_words.resize(num_common);
        end = std::chrono::high_resolution_clock::now();
        std::cerr << "Found " << common_words.size() << " matching seeds in partition " << i << " in "  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

        if (num_common == 0) {
            continue;
        }


        // find locations of common words in the index and reads
        thrust::device_vector<std::uint32_t> idx_common(num_common);
        thrust::device_vector<std::uint32_t> w_common(num_common);

        thrust::lower_bound(idx_unique_seeds.begin(), idx_unique_seeds.end(),
                            common_words.begin(), common_words.end(),
                            idx_common.begin());

        thrust::lower_bound(w_unique_seeds.begin(), w_unique_seeds.end(),
                            common_words.begin(), common_words.end(),
                            w_common.begin());


        // collect counts of common words in the index and reads
        thrust::device_vector<std::uint32_t> idx_matched_counts(num_common);
        thrust::device_vector<std::uint32_t> w_matched_counts(num_common);
        thrust::gather(idx_common.begin(), idx_common.end(),
                       idx_counts.begin(), idx_matched_counts.begin());
        thrust::gather(w_common.begin(), w_common.end(),
                       w_counts.begin(), w_matched_counts.begin());


        // find number of pairs for each matched word
        thrust::device_vector<std::uint32_t> out_counts(num_common);
        thrust::transform(idx_matched_counts.begin(), idx_matched_counts.end(),
                          w_matched_counts.begin(),
                          out_counts.begin(),
                          thrust::multiplies<std::uint32_t>());


        // compute prefix sums for matched words
        thrust::device_vector<std::uint32_t> out_offsets(num_common);
        thrust::exclusive_scan(out_counts.begin(), out_counts.end(),
                               out_offsets.begin(), std::uint32_t{0});


        // last_offset = out_offsets[num_common - 1] with explicit copy to host
        std::uint32_t last_offset = 0, last_count = 0;
        thrust::copy_n(out_offsets.begin() + (num_common - 1), 1, &last_offset);
        thrust::copy_n(out_counts.begin() + (num_common - 1), 1, &last_count);

        std::size_t total_out = last_offset + last_count;
        std::cerr << "total " << total_out << std::endl;

        // generate matches
        const CW* d_idx = thrust::raw_pointer_cast(index->d_vecs[i].data());
        const CW* d_w = thrust::raw_pointer_cast(word_vec.data());

        const std::uint32_t* d_idx_starts = thrust::raw_pointer_cast(idx_starts.data());
        const std::uint32_t* d_idx_counts = thrust::raw_pointer_cast(idx_counts.data());

        const std::uint32_t* d_w_starts = thrust::raw_pointer_cast(w_starts.data());
        const std::uint32_t* d_w_counts = thrust::raw_pointer_cast(w_counts.data());

        const std::uint32_t* d_idx_common = thrust::raw_pointer_cast(idx_common.data());
        const std::uint32_t* d_w_common = thrust::raw_pointer_cast(w_common.data());

        const std::uint32_t* d_out_offsets = thrust::raw_pointer_cast(out_offsets.data());

        SEEDMATCH* d_pairs = thrust::raw_pointer_cast(out_pairs.data() + last_pair);

        int threads = 256;
        int blocks = static_cast<int>(num_common);

        EmitCartesianProduct<<<blocks, threads>>>(
                                        d_idx, d_idx_starts, d_idx_counts,
                                        d_w, d_w_starts, d_w_counts,
                                        d_idx_common, d_w_common,
                                        d_out_offsets,
                                        d_pairs,
                                        num_common);

        last_pair += num_common;

    }
    out_pairs.resize(last_pair);
    // sort matches by read id
    thrust::sort(out_pairs.begin(), out_pairs.end(), [] __device__ (const SEEDMATCH& a, const SEEDMATCH& b) {return a.read < b.read;});

    // copy to host
//    thrust::copy(out_pairs.begin(), out_pairs.end(), cp);


//    CUDA_CHECK(cudaGetLastError());

}
