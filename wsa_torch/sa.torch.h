/*
 * sa.torch.h
 *
 * Public C interface to sa_torch.cpp (LibTorch GPU backend).
 *
 * Design rules (mirrors wsra/sra_read.h):
 *   - Pure C: compiled cleanly by both C and C++ translation units.
 *   - No acedb headers, no C++ types, no LibTorch types.
 *   - All C++ internals are hidden behind the opaque SATorchObj pointer.
 *   - The buffer returned by saTorchJoinSort / saTorchCodeJoinSort is
 *     owned by the SATorchObj and remains valid until the next call to
 *     either function or until saTorchFree().  The C caller must NOT
 *     free it.
 *
 * Output record: TORCHMATCH  (5 x uint32 = 20 bytes)
 *
 *   [0] read         read CW nam  (ia<<1 | strand_bit, ia 1-based)
 *   [1] x1           read CW pos  (1-based position in read)
 *   [2] target       genome CW nam
 *   [3] a1           genome CW pos (1-based position in genome)
 *   [4] targetFlags  genome CW flags (intron / junction bits)
 *
 * The two seed placeholder fields (readSeed, targetSeed) and readFlags
 * present in SEEDMATCH are dropped: they are always zero and the aligner
 * never uses them.  Under USE_TORCH the aligner reads TORCHMATCH directly.
 *
 * Genome index memory model (v2 — resident g_cw):
 *
 *   g_keys[p]     [K_p]     int32  sorted unique seeds   — resident
 *   g_offsets[p]  [K_p+1]   int32  CSR offset array      — resident
 *   g_cw[p]       [M_p, 3]  int32  genome payload cols   — resident
 *                                   col 0 = nam
 *                                   col 1 = pos
 *                                   col 2 = flags
 *
 * The seed column (CW col 0) is NOT stored in g_cw: the CSR index makes
 * it redundant, and dropping it reduces resident device memory from
 * ~22.8 GB to ~17.1 GB on a human/mouse index (NN=16, 1.43 B records).
 *
 * Typical call sequence (C side):
 *
 *   SATorchObj *tor = saTorchNew(device_id, seedLength, nIndex);
 *   for each partition p:
 *       saTorchIndexUpload(tor, p, seeds, K, offsets, cw, M);
 *       // cw may be freed immediately on return
 *   ...
 *   // PATH A: CPU seeds already extracted into cws_ptr[]/cws_n[]
 *   long tm_count;
 *   const uint32_t *tm = saTorchJoinSort(tor, NN,
 *                                         cws_ptr, cws_n, &tm_count);
 *   if (tm) { memcpy into bb->sms (TORCHMATCH) ... }
 *
 *   // PATH B: GPU seed extraction from flat DNA buffer
 *   const uint32_t *tm = saTorchCodeJoinSort(tor,
 *                                             dna_buf, total_bases,
 *                                             off_ptr, len_ptr, iaMax,
 *                                             &tm_count);
 *   if (tm) { memcpy into bb->sms (TORCHMATCH) ... }
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
 * Output record                                                       *
 * ------------------------------------------------------------------ */

/*
 * TORCHMATCH  — 5 fields, 20 bytes.
 * Sorted by read (field 0) ascending before return, so pairs are
 * grouped and mate-1 precedes mate-2 (same guarantee as SEEDMATCH).
 * Under USE_TORCH, sa.align.c reads this struct directly.
 */
typedef struct torchMatchStruct {
  uint32_t read ;        /* read CW nam  (ia<<1 | strand_bit) */
  uint32_t x1 ;          /* read CW pos  (1-based)            */
  uint32_t target ;      /* genome CW nam                     */
  uint32_t a1 ;          /* genome CW pos (1-based)           */
  uint32_t targetFlags ; /* genome CW flags                   */
} TORCHMATCH ;

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
 * Upload partition p of the genome index to the device as resident
 * tensors.  All three tensors (g_keys, g_offsets, g_cw) remain on the
 * device for the full run; the C caller may free its source buffers
 * immediately after this call returns.
 *
 *   p            partition index in [0, NN-1]
 *   seeds        [K]     sorted unique seeds for this partition (uint32_t)
 *                        Seeds are treated as int32 (signed) for comparison;
 *                        the on-disk index must be sorted the same way
 *                        (use saRadixSort3Signed when building the index).
 *   K            number of unique seeds
 *   offsets      [K+1]   CSR offset array (uint32_t)
 *   cw           [M*4]   genome CW records; 4 x uint32_t per record:
 *                        [seed, nam, pos, flags]
 *                        Only columns 1,2,3 (nam,pos,flags) are uploaded;
 *                        the seed column is redundant given the CSR index.
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
 * GPU join of read seeds against the resident genome index, followed
 * by sort by read nam.  CPU seed extraction (saCodeSequenceSeeds) must
 * have been called before this function.
 *
 *   NN           partition count
 *   cws_ptr[p]   pointer to read CW array for partition p
 *                (4 x uint32_t per record: seed, nam, pos, flags)
 *                NULL if partition p has no read seeds
 *   cws_n[p]     number of CW records in partition p
 *   tm_count     OUTPUT: number of TORCHMATCH records produced
 *
 * Returns a pointer to an internal buffer of (tm_count * 5) uint32_t
 * values (5 fields per TORCHMATCH).  The buffer is owned by *tor and
 * remains valid until the next call to saTorchJoinSort,
 * saTorchCodeJoinSort, or saTorchFree.
 * Returns NULL on failure; caller must fall back to CPU path.
 */
const uint32_t *saTorchJoinSort (SATorchObj            *tor,
                                 int                    NN,
                                 const uint32_t *const *cws_ptr,
                                 const long            *cws_n,
                                 long                  *tm_count) ;

/* ------------------------------------------------------------------ *
 * PATH B — GPU seed extraction + join + sort                         *
 * ------------------------------------------------------------------ */

/*
 * saTorchCodeJoinSort
 *
 * Extract seeds from a flat IUPAC DNA buffer entirely on the GPU,
 * then join against the resident genome index and sort.
 *
 *   dna_buf      flat buffer: 1 byte per base, IUPAC 4-bit (low nibble)
 *   total_bases  number of bytes in dna_buf
 *   off_ptr      [iaMax] byte offset of read ia in dna_buf (uint32_t)
 *   len_ptr      [iaMax] length of read ia in bases (uint32_t)
 *                0 = absent or rejected (too short / bad entropy)
 *   iaMax        length of off_ptr and len_ptr
 *   tm_count     OUTPUT: number of TORCHMATCH records produced
 *
 * Returns pointer to internal buffer (same ownership as saTorchJoinSort).
 * Returns NULL on failure.
 */
const uint32_t *saTorchCodeJoinSort (SATorchObj     *tor,
                                     const uint8_t  *dna_buf,
                                     long            total_bases,
                                     const uint32_t *off_ptr,
                                     const uint32_t *len_ptr,
                                     int             iaMax,
                                     long           *tm_count) ;

/* ------------------------------------------------------------------ *
 * Diagnostics                                                         *
 * ------------------------------------------------------------------ */

/*
 * saTorchIsAvailable
 * Returns 1 if a non-CPU device (CUDA, MPS) was selected, 0 otherwise.
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
