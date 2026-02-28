/* sa.gpusort.h
 * A function that sorts on a GPU callable from C
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

void saGPUSort (char *cp, long int number_of_records, int type);

#ifdef __cplusplus
}
#endif
