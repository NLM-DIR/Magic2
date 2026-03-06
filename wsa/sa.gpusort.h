/* sa.gpusort.h
 * A function that sorts on a GPU callable from C
 *
 */

#include "sa.common.h"

#ifdef __cplusplus
extern "C" {
#endif

void saGPUSort (char *cp, long int number_of_records, int type);

void saGPUMatchHits(CW** index_parts, long int* sizes, unsigned int num_parts);

#ifdef __cplusplus
}
#endif
