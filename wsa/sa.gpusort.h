/* sa.gpusort.h
 * GPU-accelerated sort and seed-matching, callable from C.
 *
 * Optimisation (April 2026) — precomputed genome BuildRuns:
 *   Genome index run-length encoding is computed once in GPUIndexCreate
 *   and reused for every read block.  The shared out_pairs device buffer
 *   is protected by pthread_mutex in sa.main.c and grows by doubling but
 *   never shrinks, avoiding per-block cudaMalloc/cudaFree overhead.
 */

#include "sa.common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ */
/* Minimal CUDA runtime declarations for plain C callers.
 * Including cuda_runtime.h in a .c file pulls in C++-only constructs. */
/* ------------------------------------------------------------------ */
typedef int cudaError_t;
int  cudaSetDevice      (int device);
int  cudaGetDeviceCount (int *count);
int  cudaHostRegister   (void *ptr, size_t size, unsigned int flags);
int  cudaHostUnregister (void *ptr);

#define cudaHostRegisterDefault  0x00U

/* ------------------------------------------------------------------ */
/* Opaque genome index handle                                          */
/* ------------------------------------------------------------------ */
typedef void GPUIndex;

/* ------------------------------------------------------------------ */
/* Radix sort (unchanged)                                              */
/* ------------------------------------------------------------------ */
void saGPUSort (char *cp, long int number_of_records, int type);

/* ------------------------------------------------------------------ */
/* Genome index lifecycle                                              */
/* ------------------------------------------------------------------ */

/* Upload all NN genome partitions to device memory and precompute their
 * run-length encodings.  Called once at startup.                      */
GPUIndex* GPUIndexCreate (CW** index_parts, long int* sizes,
                          unsigned int num_parts);

/* Release all device memory associated with the index.                */
GPUIndex* GPUIndexFree (GPUIndex* idx);

/* ------------------------------------------------------------------ */
/* Seed matching — call both under pthread_mutex_lock                  */
/* ------------------------------------------------------------------ */

/* Pass 1: run the full GPU join pipeline, sort results by read id,
 * store in the shared device buffer inside idx, return match count N.
 * Caller allocates bb.sms = bigArrayHandleCreate(N+1, SEEDMATCH, h).  */
unsigned int saGPUMatchHits (GPUIndex* idx, CW** words, long int* sizes,
                             unsigned int num_parts);

/* Pass 2: copy N SEEDMATCH records from the device buffer into the
 * caller-allocated host buffer, then reset the device buffer to full
 * capacity for the next block.                                        */
void saGPUMatchHitsCopyToHost (GPUIndex* idx, SEEDMATCH* out_buffer);

#ifdef __cplusplus
}
#endif
