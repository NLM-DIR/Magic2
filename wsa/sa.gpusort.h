/* sa.gpusort.h
 * A function that sorts on a GPU callable from C
 *
 */

#include "sa.common.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef void GPUIndex;

void saGPUSort (char *cp, long int number_of_records, int type);

GPUIndex* GPUIndexCreate(CW** index_parts, long int* sizes, unsigned int num_parts);

GPUIndex* GPUIndexFree(GPUIndex* idx);

void saGPUMatchHits(GPUIndex* idx, CW** words, long int* sizes,
                    unsigned int num_parts);

#ifdef __cplusplus
}
#endif
