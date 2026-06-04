/*
 * sa.torch.h
 *
 * Public C interface to sa_torch.cpp (LibTorch GPU backend).
 *
 * Design rules (mirrors wsra/sra_read.h):
 *   - Pure C: compiled cleanly by both C and C++ translation units.
 *   - No acedb headers, no C++ types, no LibTorch types.
 *   - All C++ internals are hidden behind the opaque SATorchObj pointer.
 *   - Buffers returned by saTorchJoinSort / saTorchCodeJoinSort are
 *     owned by the SATorchObj and remain valid until the next call to
 *     either function or until saTorchFree().  The C caller must NOT
 *     free them.
 *
 * Typical call sequence (C side):
 *
 *   SATorchObj *tor = saTorchNew(device_id, seedLength, nIndex);
 *   for each partition p:
 *       saTorchIndexUpload(tor, p, seeds, K, offsets, cw, M);
 *   ...
 *   // PATH A: CPU seeds already extracted into cws_ptr[]/cws_n[]
 *   long sm_count;
 *   const int32_t *sm = saTorchJoinSort(tor, NN,
 *                                        cws_ptr, cws_n, &sm_count);
 *   if (sm) { memcpy into bb->sms ... }
 *
 *   // PATH B: GPU seed extraction from flat DNA buffer
 *   const int32_t *sm = saTorchCodeJoinSort(tor,
 *                                            dna_buf, total_bases,
 *                                            off_ptr, len_ptr, iaMax,
 *                                            &sm_count);
 *   if (sm) { memcpy into bb->sms ... }
 *   ...
 *   saTorchFree(tor);
 *
 * Authors: Magic2 team / NLM-NIH
 * Created: June 2026
 */

#ifndef SA_TORCH_H
#define SA_TORCH_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ *
 * Opaque handle — hides all C++ and LibTorch internals               *
 * ------------------------------------------------------------------ */
typedef void SATorchObj ;

/* ------------------------------------------------------------------ *
 * Lifecycle                                                           *
 * ------------------------------------------------------------------ */

/*
 * saTorchNew
 *
 * Allocate a new SATorchObj, select the GPU device, and initialise
 * seed-packing LUTs and weight vectors for the given word length.
 *
 *   device_id   GPU index (0 for first GPU; ignored if no GPU found)
 *   wLen        k-mer length (seedLength, typically 18 for human/mouse)
 *   NN          partition count (must equal TORCH_NN = 16)
 *
 * Returns NULL on failure.
 */
SATorchObj *saTorchNew (int device_id, int wLen, int NN) ;

/*
 * saTorchFree
 *
 * Release all device tensors and the SATorchObj itself.
 * Always returns NULL so the caller can write: tor = saTorchFree(tor);
 */
SATorchObj *saTorchFree (SATorchObj *tor) ;

/* ------------------------------------------------------------------ *
 * Genome index upload — call once per partition before any block     *
 * ------------------------------------------------------------------ */

/*
 * saTorchIndexUpload
 *
 * Upload partition p of the genome index to the device.
 * The three arrays are read-only; the C caller retains ownership.
 *
 *   p            partition index in [0, NN-1]
 *   seeds        [K]     sorted unique seeds for this partition (uint32_t)
 *   K            number of unique seeds
 *   offsets      [K+1]   CSR offset array (uint32_t)
 *   cw           [M*4]   genome CW records; 4 x int32_t per record:
 *                        [seed, nam, pos, flags]
 *   M            number of CW records
 *
 * Returns 1 on success, 0 on failure.
 */
int saTorchIndexUpload (SATorchObj     *tor,
                        int             p,
                        const uint32_t *seeds,    long K,
                        const uint32_t *offsets,
                        const uint32_t *cw,       long M) ;

/* ------------------------------------------------------------------ *
 * PATH A — join + sort only (CPU has already extracted seeds)        *
 * ------------------------------------------------------------------ */

/*
 * saTorchJoinSort
 *
 * GPU join of read seeds against the genome index, followed by sort
 * by read nam.  CPU seed extraction (saCodeSequenceSeeds) must have
 * been called before this function.
 *
 *   NN           partition count
 *   cws_ptr[p]   pointer to read CW array for partition p
 *                (4 x int32_t per record: seed, nam, pos, flags)
 *                NULL if partition p has no read seeds
 *   cws_n[p]     number of CW records in partition p
 *   sm_count     OUTPUT: number of SEEDMATCH records produced
 *
 * Returns a pointer to an internal buffer of (sm_count * 8) int32_t
 * values (8 fields per SEEDMATCH).  The buffer is owned by *tor and
 * remains valid until the next call to saTorchJoinSort,
 * saTorchCodeJoinSort, or saTorchFree.
 * Returns NULL on failure; caller must fall back to CPU path.
 */
const uint32_t *saTorchJoinSort (SATorchObj            *tor,
                                int                    NN,
                                const uint32_t *const *cws_ptr,
                                const long            *cws_n,
                                long                  *sm_count) ;

/* ------------------------------------------------------------------ *
 * PATH B — GPU seed extraction + join + sort                         *
 * ------------------------------------------------------------------ */

/*
 * saTorchCodeJoinSort
 *
 * Extract seeds from a flat IUPAC DNA buffer entirely on the GPU,
 * then join against the genome index and sort.
 *
 *   dna_buf      flat buffer: 1 byte per base, IUPAC 4-bit (low nibble)
 *                Must be the globalDna buffer for the current block.
 *   total_bases  number of bytes in dna_buf
 *   off_ptr      [iaMax] byte offset of read ia in dna_buf (uint32_t)
 *                ia=0 is sentinel (offset=0, len=0)
 *   len_ptr      [iaMax] length of read ia in bases (uint32_t)
 *                0 = absent or rejected (too short / bad entropy)
 *   iaMax        length of off_ptr and len_ptr
 *   sm_count     OUTPUT: number of SEEDMATCH records produced
 *
 * Returns pointer to internal buffer (same ownership as saTorchJoinSort).
 * Returns NULL on failure.
 */
const uint32_t *saTorchCodeJoinSort (SATorchObj    *tor,
                                    const uint8_t *dna_buf,
                                    long           total_bases,
                                    const uint32_t *off_ptr,
                                    const uint32_t *len_ptr,
                                    int            iaMax,
                                    long          *sm_count) ;

/* ------------------------------------------------------------------ *
 * Diagnostics                                                         *
 * ------------------------------------------------------------------ */

/*
 * saTorchIsAvailable
 * Returns 1 if a non-CPU device (CUDA, MPS) was selected, 0 otherwise.
 * Call after saTorchNew.
 */
int saTorchIsAvailable (SATorchObj *tor) ;

/*
 * saTestTorch
 * Standalone smoke-test: allocates a tiny tensor on device_id.
 * Returns 0 on success, 1 on failure.
 * Called by magic2_torch --testTorch N; does not require saTorchNew.
 */
int saTestTorch (int device_id) ;

#ifdef __cplusplus
}
#endif

#endif /* SA_TORCH_H */
