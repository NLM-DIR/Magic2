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


typedef std::vector< thrust::device_vector<CW> > GPUIndexType;

GPUIndex* GPUIndexCreate(CW** index_parts, long int* sizes, unsigned int num_parts)
{
    GPUIndexType* d_vecs = new GPUIndexType();
    d_vecs->reserve(num_parts);
    for (unsigned int i=0;i < num_parts;i++) {
        auto start = std::chrono::high_resolution_clock::now();
        d_vecs->emplace_back(index_parts[i], index_parts[i] + sizes[i]);
        auto end = std::chrono::high_resolution_clock::now();
        std::cerr << "Copy index partition " << i << " to GPU (" << sizes[i] << "): "  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
    }

    return d_vecs;
}


GPUIndex* GPUIndexFree(GPUIndex* idx)
{
    GPUIndexType* idxobj = static_cast<GPUIndexType*>(idx);
    if (idxobj) {
        delete idxobj;
    }
    return nullptr;
}


struct diff_CW {
    __host__ __device__
    unsigned int operator()(const CW& a, const CW& b) {return b.seed - a.seed;}
};

void saGPUMatchHits(GPUIndex* idx, CW** words, long int* sizes,
                    unsigned int num_parts)
{
    GPUIndexType& index_vecs = *static_cast<GPUIndexType*>(idx);

    for (unsigned int i=0;i < num_parts;i++) {
        auto start = std::chrono::high_resolution_clock::now();
        thrust::device_vector<CW> word_vec(words[i], words[i] + sizes[i]);
        auto end = std::chrono::high_resolution_clock::now();
        std::cerr << "Copied word partition " << i << " to GPU (" << sizes[i] << "): "  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

        start = std::chrono::high_resolution_clock::now();
        thrust::sort(word_vec.begin(), word_vec.end(), compare_CW());
        end = std::chrono::high_resolution_clock::now();
        std::cerr << "Sorted word partition " << i << " in "  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

        start = std::chrono::high_resolution_clock::now();
        thrust::device_vector<std::size_t> upper_it(word_vec.size());
        auto it = thrust::upper_bound(index_vecs[i].begin(), index_vecs[i].end(), word_vec.begin(), word_vec.end(), upper_it.begin(), compare_CW());
        end = std::chrono::high_resolution_clock::now();
        std::cerr << "Found matching seeds " << i << " in "  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

        /* for debugging, make sure that upper_bound finds what we need */
        auto word_it = word_vec.begin();
	auto idx_it = upper_it.begin();
	for (int kk=0;kk < 10;kk++,word_it+=20,idx_it+=20) {
            if (idx_it != upper_it.end()) {
		CW w = *word_it;
		//auto ind = index_vecs[i].begin() + *idx_it;
                CW d = index_vecs[i][*idx_it];
		//CW dd = index_vecs[i][*idx_it - 1];
                std::cerr << "seeds:\t" << w.seed << "\t"
                          << d.seed << "\t"
                          << *idx_it
                          << std::endl;
            }
	}



        /* these may come handy */
        start = std::chrono::high_resolution_clock::now();
        thrust::device_vector<CW> diffs(word_vec.size());
        thrust::adjacent_difference(word_vec.begin(), word_vec.end(), diffs.begin(), [] __device__ (const CW& a, const CW& b) {CW r; r.seed = b.seed - a.seed; return r;});
        end = std::chrono::high_resolution_clock::now();
        std::cerr << "Diffs in partition " << i << " in "  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

        thrust::device_vector<std::size_t> positions(diffs.size());
        auto first = thrust::make_counting_iterator(0);
        auto last = first + diffs.size();
        auto last_pos = thrust::copy_if(first, last, diffs.begin(),
                                        positions.begin(),
                        [] __device__ (const CW& x) {return x.seed != 0;});

        thrust::device_vector<CW> unique_hashes(last_pos - positions.begin());
        thrust::gather(positions.begin(), last_pos, word_vec.begin(), unique_hashes.begin());



    }

}
