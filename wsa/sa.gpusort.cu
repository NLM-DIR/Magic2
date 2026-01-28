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
#include <chrono>
#include <iostream>

#include "sa.gpusort.h"
#include "sa.common.h"


// Comparators for thrust::sort(), reimplemented versions of comparators in
// sa.sort.c so that the compiler can better optimize the code.
struct compare_CW {
    // compare code words, same as cwOrder in sa.sort.c
    __host__ __device__
    bool operator()(const CW& a, const CW& b)
    {return a.seed <= b.seed ; }  // we do not need a more detailed comparisons

// {return a.seed < b.seed || (a.seed == b.seed && a.nam < b.nam) ||
//    (a.seed == b.seed && a.nam == b.nam && a.pos < b.pos);}
};


struct compare_HIT {
    // compare hits, same as hitOrder in sa.sort.c
    __host__ __device__
    bool operator()(const HIT& a, const HIT& b)
    {return a.read < b.read || a.read == b.read && a.chrom < b.chrom ||
    a.read == b.read && a.chrom == b.chrom && a.a1 < b.a1 ||
    a.read == b.read && a.chrom == b.chrom && a.a1 == b.a1 && a.x1 < b.x1;}
};


struct compare_HIT_pairs {
    // compare hits for read pairs, same as hitPairOrder in sa.sort.c
    __host__ __device__
    bool operator()(const HIT& a, const HIT& b)
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
