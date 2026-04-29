/* sa.gpusort.h
 * GPU-accelerated sort and seed-matching, callable from C.
 *
 * Pinned memory strategy (April 2026)
 * ------------------------------------
 * All host buffers exchanged between the GPU pipeline and the C pipeline
 * (SEEDMATCH output, and optionally CW input) are allocated with
 * cudaMallocHost.  Pinned memory:
 *
 *   - transfers to/from the GPU at full PCIe bandwidth (no driver bounce
 *     buffer), typically 2x faster than pageable memory transfers;
 *   - is guaranteed to be aligned to at least 256 bytes, satisfying both
 *     SEEDMATCH (__attribute__((aligned(32)))) and CW (__attribute__((aligned(16)))).
 *
 * Pinned buffers MUST be freed with saGPUFreeHostBuffer(), not free().
 * Use bigArraySwitchCudaBase(aa, ptr, saGPUFreeHostBuffer) so that bigArray's
 * destructor calls saGPUFreeHostBuffer automatically.
 *
 * Concurrency
 * -----------
 * saGPUMatchHits allocates all writable device memory locally per call.
 * The genome index (GPUIndex) is read-only after GPUIndexCreate returns.
 * Concurrent calls from different agents share no writable state and
 * require NO MUTEX.
 */

#include "sa.common.h"

#ifdef __cplusplus
extern "C" {
#endif

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
GPUIndex* GPUIndexCreate (CW** index_parts, long int* sizes, unsigned int num_parts);

/* Release all device memory associated with the index.                */
GPUIndex* GPUIndexFree   (GPUIndex* idx);

/* ------------------------------------------------------------------ */
/* Pinned host memory helpers                                          */
/* ------------------------------------------------------------------ */

/* Release a cudaMallocHost buffer allocated by saGPUMatchHits.
 * Call this explicitly in sa.main.c before the standard bigArray
 * destructor runs, so that cudaFreeHost stays localised to the magic2
 * GPU layer and the bigArray library needs no USEGPU ifdefs.          */
void saGPUFreeHostBuffer (void* ptr);

/* ------------------------------------------------------------------ */
/* Seed matching                                                       */
/* ------------------------------------------------------------------ */

/* Find all SEEDMATCH records between the read block (words/sizes, NN
 * partitions) and the resident genome index, sort by read id, allocate
 * a pinned host buffer of exactly N records, copy the results into it,
 * and return N.
 *
 * *out_buffer receives a pointer to the pinned SEEDMATCH array.
 * The caller must attach it to a bigArray with bigArraySwitchCudaBase
 * so that destruction calls saGPUFreeHostBuffer.
 *
 * If words[i] was allocated via saGPUAllocHostCW (pinned), the PCIe
 * upload of each partition runs at full bandwidth.
 *
 * No mutex is needed: all writable device state is local to this call. */
unsigned int saGPUMatchHits (GPUIndex* idx, CW** words, long int* sizes,
                             unsigned int num_parts, SEEDMATCH** out_buffer);

#ifdef __cplusplus
}
#endif
