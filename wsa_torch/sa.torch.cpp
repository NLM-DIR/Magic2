/*
 * sa_torch.cpp
 *
 * LibTorch GPU backend for Magic2 sort-and-join pipeline.
 *
 * Design rules (mirrors wsra/sra_read.cpp):
 *   - Includes ONLY LibTorch headers and sa.torch.h.
 *   - Never includes sa.h, array.h, ac.h or any acedb header.
 *   - All data crosses the C/C++ boundary as flat pointers + counts.
 *   - Output buffers are owned by SATorchState and loaned to the C caller;
 *     never freed by the C side.
 *
 * Two public entry points for Agent S (Sorter):
 *
 *   saTorchJoinSort      PATH A: CPU has filled cws_ptr[p] arrays.
 *                                GPU does join + sort only.
 *
 *   saTorchCodeJoinSort  PATH B: GPU extracts seeds from dna_buf,
 *                                then join + sort.
 *
 * Both paths share sTorchJoinSortDo() for the join+sort body.
 *
 * Genome index memory layout:
 *
 *   g_keys[p]     [K_p]     uint32   sorted unique seeds — resident on device
 *   g_offsets[p]  [K_p+1]   uint32   CSR offset array   — resident on device
 *   h_cw[p]                 uint32*  raw pointer into mmap'd cwsN on bbG->h
 *   h_M[p]                  long     number of CW records in partition p
 *
 * g_cw is NOT kept resident on the device: it is uploaded partition by
 * partition inside sTorchJoinSortDo, used for the join, then released.
 * This keeps peak device memory at max(M_p) * 16 bytes rather than
 * sum(M_p) * 16 bytes, making a 32 GB GPU viable for human/mouse.
 *
 * The cwsN mmap buffers contain jumper lines at every row index divisible
 * by 256.  Jumper lines store a look-ahead seed (used by the CPU sort-merge
 * for fast pointer advancement) and must be excluded from the GPU Cartesian
 * product.  They are filtered by row index inside sTorchJoinSortDo — no
 * extra data field is needed because the GPU knows its own row numbers.
 *
 * CW field layout (4 x uint32, matches sa.h CW struct):
 *   [0] seed    sort/lookup key (also used in jumper lines — excluded by row)
 *   [1] nam     genome sequence id
 *   [2] pos     1-based coordinate of first base of seed in genome
 *   [3] flags   auxiliary (passenger — copied verbatim to SEEDMATCH)
 *
 * SEEDMATCH layout (8 x uint32, matches sa.h SEEDMATCH struct):
 *   [0] readSeed   = 0 placeholder
 *   [1] read       = read CW nam    (used for sort; ia 1-based, strand in bit 0)
 *   [2] x1         = read CW pos    (passenger)
 *   [3] readFlags  = read CW flags  (passenger)
 *   [4] targetSeed = 0 placeholder
 *   [5] target     = genome CW nam  (passenger)
 *   [6] a1         = genome CW pos  (passenger)
 *   [7] targetFlags= genome CW flags (passenger)
 *
 * Seed encoding (matches saCodeSequenceSeedsStep1 in sa.seeds.c):
 *   2-bit fwd: A=0 C=1 G=2 T=3
 *   Z = min(fwd_kmer, rc_kmer)  strand-agnostic
 *   strand_bit = 1 if rc < fwd
 *   nam = (ia << 1) | strand_bit   ia is 1-based index in bb->dnas
 *   pos = 1-based position of first base of seed in read
 *   partition  = Z & (NN-1)
 *   stored     = (Z >> nHidden2) & 0xFFFFFFFF
 *   nHidden2   = 2*(wLen-16) for wLen>16, else 0
 *
 * Authors: Magic2 team / NLM-NIH
 * Created: June 2026
 */

/* LibTorch — must come first, before any system headers that might
 * be pulled in with conflicting macros                               */
#include <torch/torch.h>

/* Standard C++ headers                                               */
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <vector>

/* Public C interface — pure C header, no acedb types               */
#include "sa.torch.h"

/* ================================================================== *
 * §1  Internal state object                                          *
 * ================================================================== */

/*
 * SATorchState holds all C++ / LibTorch state for one run.
 * Exposed to C as the opaque SATorchObj pointer.
 */
struct SATorchState {

  /* ---- device ---------------------------------------------------- */
  torch::Device dev ;

  /* ---- genome index: keys and offsets resident on device --------- */
  static constexpr int TORCH_NN = 16 ;
  torch::Tensor g_keys    [TORCH_NN] ;   /* [K_p]   uint32, on device  */
  torch::Tensor g_offsets [TORCH_NN] ;   /* [K_p+1] uint32, on device  */
  bool          g_ready = false ;

  /* ---- genome CW records: mmap'd on CPU, streamed per join call -- */
  /* Pointers into the C caller's mmap'd cwsN arrays (bbG->h).       *
   * The C caller owns the memory and keeps it alive for the run.    *
   * sa_torch.cpp never frees these pointers.                        */
  const uint32_t *h_cw [TORCH_NN] ;     /* raw pointer per partition  */
  long            h_M  [TORCH_NN] ;     /* CW record count            */

  /* ---- seed-packing tables (device tensors) --------------------- */
  torch::Tensor g_fwd_lut ;      /* [16] uint32                       */
  torch::Tensor g_valid_lut ;    /* [16] uint32                       */
  torch::Tensor g_fwd_weights ;  /* [wLen] int64  4^(wLen-1-j)        */
  torch::Tensor g_rev_weights ;  /* [wLen] int64  4^j                 */
  bool          init_done = false ;
  int           wLen = 0 ;
  int           NN   = 0 ;

  /* ---- output buffer (loaned to C caller, never freed by C) ----- */
  std::vector<uint32_t> sm_buf ;

  /* ---- constructor ---------------------------------------------- */
  explicit SATorchState (torch::Device d) : dev(d)
  {
    for (int p = 0 ; p < TORCH_NN ; p++) {
      h_cw[p] = nullptr ;
      h_M [p] = 0 ;
    }
  }
} ;

/* ================================================================== *
 * §2  IUPAC look-up tables (CPU-side constants, uploaded once)      *
 * ================================================================== */

/* 4-bit IUPAC nibble → 2-bit fwd base (A=0 C=1 G=2 T=3)           */
static const uint32_t kFwdLutData[16] = {
  0, 0, 3, 0,
  2, 0, 0, 0,
  1, 0, 0, 0,
  0, 0, 0, 0
} ;

/* 4-bit IUPAC nibble → 1 if unambiguous ACGT, else 0               */
static const uint32_t kValidLutData[16] = {
  0, 1, 1, 0,
  1, 0, 0, 0,
  1, 0, 0, 0,
  0, 0, 0, 0
} ;

/* ================================================================== *
 * §3  TensorOptions helpers                                          *
 * ================================================================== */

static inline torch::TensorOptions kInt32CPU ()
{
  return torch::TensorOptions ()
           .dtype (torch::kInt32)
           .device (torch::kCPU) ;
}
static inline torch::TensorOptions kInt64Dev (const torch::Device &dev)
{
  return torch::TensorOptions ()
           .dtype (torch::kInt64)
           .device (dev) ;
}
static inline torch::TensorOptions kInt32Dev (const torch::Device &dev)
{
  return torch::TensorOptions ()
           .dtype (torch::kInt32)
           .device (dev) ;
}

/* ================================================================== *
 * §4  Device selection                                               *
 * ================================================================== */

static torch::Device sTorchSelectDevice (int device_id)
{
  if (torch::cuda::is_available ()) {
    int n = torch::cuda::device_count () ;
    if (device_id < 0 || device_id >= n) device_id = 0 ;
    std::fprintf (stderr,
                  "[sa_torch] device: CUDA %d\n", device_id) ;
    return torch::Device (torch::kCUDA, device_id) ;
  }
#ifdef TORCH_MPS_AVAILABLE
  if (torch::mps::is_available ()) {
    std::fprintf (stderr, "[sa_torch] device: MPS (Apple)\n") ;
    return torch::Device (torch::kMPS) ;
  }
#endif
  std::fprintf (stderr, "[sa_torch] device: CPU (no GPU found)\n") ;
  return torch::Device (torch::kCPU) ;
}

/* ================================================================== *
 * §5  Join+sort body (shared by PATH A and PATH B)                  *
 * ================================================================== */

/*
 * sTorchJoinSortDo
 *
 * Precondition: s->g_ready == true.
 *
 * For each partition p:
 *   - g_keys[p] and g_offsets[p] are resident on the device.
 *   - h_cw[p] points to the mmap'd CW records on the CPU (bbG->h).
 *     The CW array contains jumper lines at every row divisible by 256;
 *     these are excluded from the Cartesian product by row-index masking.
 *   - cws_ptr[p] / cws_n[p] are the read-side CW arrays for partition p.
 *
 * The genome CW tensor is uploaded, used, and released partition by
 * partition so that peak device memory is bounded by the largest single
 * partition rather than the total index size.
 *
 * On success: fills s->sm_buf with SEEDMATCH records (8 x uint32 each),
 * returns the total record count.
 * On failure: returns -1.
 */
static long sTorchJoinSortDo (SATorchState           *s,
                               int                     NN,
                               const uint32_t *const  *cws_ptr,
                               const long             *cws_n)
{
  std::vector<torch::Tensor> sm_parts ;

  for (int p = 0 ; p < NN ; p++) {

    if (!cws_ptr[p] || cws_n[p] <= 0) continue ;
    if (!s->h_cw[p] || s->h_M[p] <= 0) continue ;

    long N = cws_n[p] ;

    /* 1. Wrap read CWs as CPU tensor, move to device               */
    torch::Tensor R =
      torch::from_blob (const_cast<uint32_t *>(cws_ptr[p]),
                        {N, 4}, kInt32CPU ())
      .to (s->dev) ;                             /* [N, 4] int32    */

    /* 2. Extract seed column (column 0)                             */
    torch::Tensor r_seeds =
      R.select (1, 0).contiguous () ;            /* [N] int32       */

    /* 3. Bucketize: floor position of each read seed in genome keys */
    torch::Tensor buckets =
      torch::bucketize (r_seeds, s->g_keys[p],
                        /*out_int32=*/false,
                        /*right=*/false) ;        /* [N] int64       */

    /* 4. Verify exact match                                         */
    int64_t K_p = s->g_keys[p].size (0) ;
    torch::Tensor b_clamped = buckets.clamp (0, K_p - 1) ;

    torch::Tensor genome_seeds =
      s->g_keys[p].index_select (0, b_clamped) ; /* [N] int32       */

    torch::Tensor exact_match =
      genome_seeds.eq (r_seeds) ;                /* [N] bool        */

    torch::Tensor r_matched =
      exact_match.nonzero ().squeeze (1) ;        /* [n_match] int64 */

    if (r_matched.numel () == 0) continue ;

    /* 5. Gather CSR starts and counts for matched read seeds        */
    torch::Tensor b_matched =
      b_clamped.index_select (0, r_matched) ;    /* [n_match] int64 */

    torch::Tensor starts =
      s->g_offsets[p]
        .index_select (0, b_matched)
        .to (torch::kInt64) ;                    /* [n_match] int64 */

    torch::Tensor counts =
      (s->g_offsets[p]
         .index_select (0, b_matched + 1)
         .to (torch::kInt64)
       - starts)
      .clamp_min (0) ;                           /* [n_match] int64 */

    int64_t total = counts.sum ().item<int64_t> () ;
    if (total == 0) continue ;

    /* 6. Expand read rows by hit count (Cartesian product LHS)      */
    torch::Tensor R_matched =
      R.index_select (0, r_matched) ;            /* [n_match, 4]    */

    torch::Tensor R_exp =
      R_matched.repeat_interleave (counts, 0) ;  /* [total, 4]      */

    /* 7. Build genome row index vector                              */
    torch::Tensor group_start =
      starts.repeat_interleave (counts) ;        /* [total] int64   */

    torch::Tensor cum = counts.cumsum (0) ;
    torch::Tensor group_base =
      torch::cat ({
          torch::zeros ({1}, kInt64Dev (s->dev)),
          cum.slice (0, 0, -1)
      }) ;                                       /* [n_match] int64 */

    torch::Tensor group_base_exp =
      group_base.repeat_interleave (counts) ;    /* [total] int64   */

    torch::Tensor inner_offset =
      torch::arange (total, kInt64Dev (s->dev))
      - group_base_exp ;                         /* [total] int64   */

    torch::Tensor g_rows =
      (group_start + inner_offset)
      .clamp (0, s->h_M[p] - 1) ;               /* [total] int64   */

    /* 8. Filter jumper lines: exclude every row where (row&255)==0  *
     * Jumper lines store look-ahead seeds for the CPU sort-merge    *
     * and must not produce SEEDMATCH records on the GPU path.       *
     * Row 0 of every partition is also a jumper and is excluded.    */
    torch::Tensor jumper_keep =
      (g_rows & 255).ne (0).nonzero ().squeeze (1) ; /* [n_keep] int64 */

    if (jumper_keep.numel () == 0) continue ;

    g_rows = g_rows.index_select (0, jumper_keep) ; /* [n_keep] int64 */
    R_exp  = R_exp .index_select (0, jumper_keep) ; /* [n_keep, 4]    */

    /* 9. Upload genome CW partition, gather rows, release device    *
     * tensor immediately to bound peak GPU memory usage.            */
    torch::Tensor g_cw_p =
      torch::from_blob (const_cast<uint32_t *>(s->h_cw[p]),
                        {s->h_M[p], 4}, kInt32CPU ())
      .to (s->dev).clone () ;                    /* [M_p, 4] int32  */

    torch::Tensor G_exp =
      g_cw_p.index_select (0, g_rows) ;         /* [n_keep, 4]     */

    g_cw_p = torch::Tensor () ;                 /* release device memory */

    /* 10. Assemble SEEDMATCH rows (8 columns)                       */
    int64_t n_keep = G_exp.size (0) ;
    torch::Tensor zeros =
      torch::zeros ({n_keep, 1}, kInt32Dev (s->dev)) ;

    torch::Tensor SM = torch::cat ({
        zeros,                                   /* [0] readSeed  = 0     */
        R_exp.narrow (1, 1, 1),                  /* [1] read nam          */
        R_exp.narrow (1, 2, 1),                  /* [2] x1  (pos)         */
        R_exp.narrow (1, 3, 1),                  /* [3] readFlags         */
        zeros,                                   /* [4] targetSeed = 0    */
        G_exp.narrow (1, 1, 1),                  /* [5] target nam        */
        G_exp.narrow (1, 2, 1),                  /* [6] a1  (pos)         */
        G_exp.narrow (1, 3, 1),                  /* [7] targetFlags       */
    }, /*dim=*/1) ;                              /* [n_keep, 8]           */

    sm_parts.push_back (SM) ;
  }

  if (sm_parts.empty ()) {
    s->sm_buf.clear () ;
    return 0 ;
  }

  /* 11. Concatenate all partitions                                   */
  torch::Tensor SM_all =
    (sm_parts.size () == 1)
    ? sm_parts[0]
    : torch::cat (sm_parts, 0) ;

  sm_parts.clear () ;

  /* 12. Sort by read nam (column 1)                                  */
  torch::Tensor read_col =
    SM_all.select (1, 1).contiguous () ;
  torch::Tensor sort_idx = torch::argsort (read_col) ;
  torch::Tensor SM_sorted =
    SM_all.index_select (0, sort_idx) ;

  /* 13. Copy back to host into sm_buf                               */
  torch::Tensor SM_cpu =
    SM_sorted.to (torch::kCPU).contiguous () ;

  int64_t Total = SM_cpu.size (0) ;
  s->sm_buf.resize ((size_t)(Total * 8)) ;
  std::memcpy (s->sm_buf.data (),
               SM_cpu.data_ptr<int32_t> (),
               (size_t)Total * 8 * sizeof (uint32_t)) ;

  return (long)Total ;
}

/* ================================================================== *
 * §6  GPU seed extraction (PATH B)                                   *
 * ================================================================== */

/*
 * sTorchSeedExtract
 *
 * Extract k-mer seeds from a flat IUPAC DNA buffer on the GPU.
 * Fills cws_out[p] and cws_n_out[p] for each partition p.
 *
 * iStep is fixed at 1 on the GPU (see architecture doc §6):
 * the CPU min-hash rolling phasing has no GPU equivalent; any stride
 * greater than 1 creates blind spots rather than a speed/sensitivity
 * trade-off.  Extra seeds cost nothing in the join step and are
 * filtered naturally by the genome index.
 *
 * ia in nam is 1-based (matches saCodeSequenceSeedsStep1):
 *   nam = (ia << 1) | strand_bit   ia in [1, iaMax-1]
 * ia=0 is the sentinel; it contributes zero windows (len=0) and
 * never appears in the output.
 *
 * Returns true on success.
 */
static bool sTorchSeedExtract (
    SATorchState                      *s,
    const uint8_t                     *dna_buf,
    long                               total_bases,
    const uint32_t                    *off_ptr,
    const uint32_t                    *len_ptr,
    int                                iaMax,
    std::vector<std::vector<uint32_t>> &cws_out,
    std::vector<long>                  &cws_n_out)
{
  const int wLen     = s->wLen ;
  const int NN       = s->NN ;
  const int nHidden2 = (wLen > 16) ? 2 * (wLen - 16) : 0 ;

  if (total_bases <= 0 || iaMax <= 1) return false ;

  /* 1. Upload flat DNA buffer and decode to 2-bit bases             */
  torch::Tensor raw_cpu =
    torch::from_blob (const_cast<uint8_t *>(dna_buf),
                      {total_bases},
                      torch::TensorOptions ()
                        .dtype (torch::kUInt8)
                        .device (torch::kCPU)) ;

  torch::Tensor raw     = raw_cpu.to (s->dev) ;      /* [total] uint8  */
  torch::Tensor nibbles = raw.to (torch::kInt32) & 0x0F ;

  torch::Tensor bases2  =
    s->g_fwd_lut.index_select (0, nibbles) ;         /* [total] int32  */
  torch::Tensor valid1  =
    s->g_valid_lut.index_select (0, nibbles) ;       /* [total] int32  */

  /* 2. Upload offset and length arrays                              */
  torch::Tensor seq_off =
    torch::from_blob (const_cast<uint32_t *>(off_ptr),
                      {(long)iaMax}, kInt32CPU ())
    .to (s->dev) ;                                   /* [iaMax] int32  */

  torch::Tensor seq_len =
    torch::from_blob (const_cast<uint32_t *>(len_ptr),
                      {(long)iaMax}, kInt32CPU ())
    .to (s->dev) ;                                   /* [iaMax] int32  */

  /* 3. Number of windows per read (zero for sentinel ia=0)          */
  torch::Tensor seq_len64 = seq_len.to (torch::kInt64) ;
  torch::Tensor n_win =
    (seq_len64 - wLen + 1).clamp_min (0) ;           /* [iaMax] int64  */

  int64_t N_win_total = n_win.sum ().item<int64_t> () ;
  if (N_win_total == 0) {
    cws_out.assign (NN, {}) ;
    cws_n_out.assign (NN, 0) ;
    return true ;
  }

  /* 4. Build flat window-start indices                              */
  torch::Tensor seq_off64 = seq_off.to (torch::kInt64) ;

  torch::Tensor read_base =
    seq_off64.repeat_interleave (n_win) ;            /* [N_win] int64  */

  /* ia_idx is 1-based: ia=0 is sentinel, contributes zero windows,  *
   * and never appears in ia_exp.  This matches the C encoding:      *
   * nam = (ia << 1) | strand_bit with ia 1-based.                  */
  torch::Tensor ia_idx =
    torch::arange (1, (int64_t)iaMax + 1,
                   kInt64Dev (s->dev)) ;             /* [iaMax] int64  */

  torch::Tensor ia_exp =
    ia_idx.repeat_interleave (n_win) ;               /* [N_win] int64  */

  torch::Tensor cum = n_win.cumsum (0) ;
  torch::Tensor group_base =
    torch::cat ({
        torch::zeros ({1}, kInt64Dev (s->dev)),
        cum.slice (0, 0, -1)
    }) ;
  torch::Tensor group_base_exp =
    group_base.repeat_interleave (n_win) ;
  torch::Tensor inner =
    torch::arange (N_win_total, kInt64Dev (s->dev))
    - group_base_exp ;                               /* [N_win] int64  */

  torch::Tensor win_start = read_base + inner ;      /* [N_win] int64  */

  /* 5. Gather k-base windows                                        */
  torch::Tensor offsets_w =
    torch::arange ((int64_t)wLen, kInt64Dev (s->dev)) ;

  torch::Tensor win_idx =
    win_start.unsqueeze (1) + offsets_w ;            /* [N_win, wLen]  */

  win_idx = win_idx.clamp (0, total_bases - 1) ;

  torch::Tensor flat_idx = win_idx.reshape ({-1}) ;

  torch::Tensor W =
    bases2.index_select (0, flat_idx)
          .reshape ({N_win_total, (int64_t)wLen}) ;  /* [N_win, wLen]  */

  torch::Tensor V =
    valid1.index_select (0, flat_idx)
          .reshape ({N_win_total, (int64_t)wLen}) ;

  /* 6. Validity mask: all wLen bases must be unambiguous ACGT       */
  torch::Tensor all_valid =
    std::get<0> (V.min (1)).eq (1) ;                /* [N_win] bool   */

  /* 7. Pack forward and RC k-mer values                             */
  torch::Tensor W64 = W.to (torch::kInt64) ;

  torch::Tensor S =
    (W64 * s->g_fwd_weights).sum (1) ;              /* [N_win] int64  */

  torch::Tensor Wc = 3 - W64 ;
  torch::Tensor C  =
    (Wc * s->g_rev_weights).sum (1) ;               /* [N_win] int64  */

  int64_t maskSeedLn = (wLen < 32)
    ? ((int64_t)1 << (2 * wLen)) - 1
    : (int64_t)0x7FFFFFFFFFFFFFFFLL ;
  S = S & maskSeedLn ;
  C = C & maskSeedLn ;

  /* 8. Strand-agnostic Z, partition, stored seed                    */
  torch::Tensor minus  = C.lt (S) ;
  torch::Tensor Z      = torch::where (minus, C, S) ;

  int64_t maskNN = (int64_t)(NN - 1) ;
  torch::Tensor part =
    (Z & maskNN).to (torch::kInt32) ;

  torch::Tensor z_stored =
    (torch::bitwise_right_shift (Z, nHidden2) & (int64_t)0xFFFFFFFF)
    .to (torch::kInt32) ;

  /* 9. Assemble CW columns                                          */
  torch::Tensor strand_bit = minus.to (torch::kInt32) ;

  /* ia_exp already holds 1-based ia values (see step 4 comment)    */
  torch::Tensor nam =
    torch::bitwise_left_shift (ia_exp.to (torch::kInt32), 1)
    | strand_bit ;                                   /* [N_win] int32  */

  torch::Tensor pos =
    (inner + 1).to (torch::kInt32) ;                /* 1-based        */

  torch::Tensor flags =
    torch::zeros ({N_win_total}, kInt32Dev (s->dev)) ;

  /* 10. Apply validity filter                                       */
  torch::Tensor keep =
    all_valid.nonzero ().squeeze (1) ;               /* [n_valid] int64 */

  if (keep.numel () == 0) {
    cws_out.assign (NN, {}) ;
    cws_n_out.assign (NN, 0) ;
    return true ;
  }

  torch::Tensor seed_f  = z_stored.index_select (0, keep) ;
  torch::Tensor nam_f   = nam     .index_select (0, keep) ;
  torch::Tensor pos_f   = pos     .index_select (0, keep) ;
  torch::Tensor flags_f = flags   .index_select (0, keep) ;
  torch::Tensor part_f  = part    .index_select (0, keep) ;

  /* 11. Scatter into NN partition vectors                           */
  cws_out.resize (NN) ;
  cws_n_out.resize (NN, 0) ;

  for (int p = 0 ; p < NN ; p++) {

    torch::Tensor pmask =
      part_f.eq (p).nonzero ().squeeze (1) ;

    if (pmask.numel () == 0) {
      cws_out[p].clear () ;
      cws_n_out[p] = 0 ;
      continue ;
    }

    torch::Tensor CW_p = torch::stack ({
        seed_f .index_select (0, pmask),
        nam_f  .index_select (0, pmask),
        pos_f  .index_select (0, pmask),
        flags_f.index_select (0, pmask)
    }, /*dim=*/1) ;                                 /* [N_p, 4] int32 */

    torch::Tensor CW_cpu = CW_p.to (torch::kCPU).contiguous () ;
    int64_t N_p = CW_cpu.size (0) ;

    cws_out[p].resize ((size_t)(N_p * 4)) ;
    std::memcpy (cws_out[p].data (),
                 CW_cpu.data_ptr<int32_t> (),
                 (size_t)N_p * 4 * sizeof (uint32_t)) ;
    cws_n_out[p] = (long)N_p ;
  }

  return true ;
}

/* ================================================================== *
 * §7  extern "C" implementations                                     *
 * ================================================================== */

extern "C" {

/* ------------------------------------------------------------------ *
 * saTorchNew                                                         *
 * ------------------------------------------------------------------ */
SATorchObj *saTorchNew (int device_id, int wLen, int NN)
{
  try {
    torch::Device dev = sTorchSelectDevice (device_id) ;
    SATorchState *s   = new SATorchState (dev) ;

    s->wLen = wLen ;
    s->NN   = NN   ;

    /* LUT tensors                                                   */
    s->g_fwd_lut =
      torch::from_blob (const_cast<uint32_t *>(kFwdLutData),
                        {16}, kInt32CPU ())
      .to (dev).clone () ;

    s->g_valid_lut =
      torch::from_blob (const_cast<uint32_t *>(kValidLutData),
                        {16}, kInt32CPU ())
      .to (dev).clone () ;

    /* Forward weight vector: fwd_weights[j] = 4^(wLen-1-j)        */
    {
      std::vector<int64_t> fw ((size_t)wLen) ;
      int64_t pw = 1 ;
      for (int j = wLen - 1 ; j >= 0 ; j--, pw *= 4) fw[(size_t)j] = pw ;
      s->g_fwd_weights =
        torch::from_blob (fw.data (), {(long)wLen},
                          torch::TensorOptions ()
                            .dtype (torch::kInt64)
                            .device (torch::kCPU))
        .to (dev).clone () ;
    }

    /* Reverse weight vector: rev_weights[j] = 4^j                  */
    {
      std::vector<int64_t> rw ((size_t)wLen) ;
      int64_t pw = 1 ;
      for (int j = 0 ; j < wLen ; j++, pw *= 4) rw[(size_t)j] = pw ;
      s->g_rev_weights =
        torch::from_blob (rw.data (), {(long)wLen},
                          torch::TensorOptions ()
                            .dtype (torch::kInt64)
                            .device (torch::kCPU))
        .to (dev).clone () ;
    }

    s->init_done = true ;
    std::fprintf (stderr,
                  "[sa_torch] saTorchNew: wLen=%d NN=%d ready\n",
                  wLen, NN) ;
    return static_cast<SATorchObj *>(s) ;

  } catch (const std::exception &e) {
    std::fprintf (stderr,
                  "[sa_torch] saTorchNew: exception: %s\n", e.what ()) ;
    return nullptr ;
  }
}

/* ------------------------------------------------------------------ *
 * saTorchFree                                                        *
 * ------------------------------------------------------------------ */
SATorchObj *saTorchFree (SATorchObj *tor)
{
  if (tor) {
    SATorchState *s = static_cast<SATorchState *>(tor) ;
    delete s ;
  }
  return nullptr ;
}

/* ------------------------------------------------------------------ *
 * saTorchIsAvailable                                                 *
 * ------------------------------------------------------------------ */
int saTorchIsAvailable (SATorchObj *tor)
{
  if (!tor) return 0 ;
  SATorchState *s = static_cast<SATorchState *>(tor) ;
  return (s->dev.type () != torch::kCPU) ? 1 : 0 ;
}

/* ------------------------------------------------------------------ *
 * saTestTorch                                                        *
 * ------------------------------------------------------------------ */
int saTestTorch (int device_id)
{
  try {
    torch::Device dev = sTorchSelectDevice (device_id) ;
    if (dev.type () == torch::kCPU) {
      std::fprintf (stderr, "[sa_torch] --testTorch: no GPU found\n") ;
      return 1 ;
    }
    torch::Tensor t =
      torch::zeros ({1}, torch::TensorOptions ()
                           .dtype (torch::kFloat32)
                           .device (dev)) ;
    (void)t ;
    std::fprintf (stderr, "[sa_torch] --testTorch: OK\n") ;
    return 0 ;
  } catch (const std::exception &e) {
    std::fprintf (stderr,
                  "[sa_torch] --testTorch: %s\n", e.what ()) ;
    return 1 ;
  }
}

/* ------------------------------------------------------------------ *
 * saTorchIndexUpload                                                 *
 *                                                                    *
 * Upload strategy:                                                   *
 *   - g_keys and g_offsets are small (K*4 and (K+1)*4 bytes) and    *
 *     are copied to device with .clone() — the C caller may         *
 *     ac_free the source mmap buffers immediately on return.        *
 *   - cw (cwsN) is large (~4 GB for partition 0 on human/mouse)     *
 *     and is NOT copied to device here.  The raw pointer is stored  *
 *     in h_cw[p] and the C caller must keep the mmap buffer alive   *
 *     for the duration of the run (allocate on bbG->h, same as the  *
 *     CPU path).  The tensor is uploaded per-block inside            *
 *     sTorchJoinSortDo and released immediately after use.          *
 * ------------------------------------------------------------------ */
int saTorchIndexUpload (SATorchObj     *tor,
                        int             p,
                        const uint32_t *seeds,    long K,
                        const uint32_t *offsets,
                        const uint32_t *cw,       long M)
{
  if (!tor) return 0 ;
  SATorchState *s = static_cast<SATorchState *>(tor) ;

  if (p < 0 || p >= s->TORCH_NN) {
    std::fprintf (stderr,
                  "[sa_torch] partition %d out of range (NN=%d)\n",
                  p, s->TORCH_NN) ;
    return 0 ;
  }

  try {
    /* Keys and offsets: small — copy to device, source can be freed */
    torch::Tensor cpu_keys =
      torch::from_blob (const_cast<uint32_t *>(seeds),
                        {K}, kInt32CPU ()) ;
    torch::Tensor cpu_offsets =
      torch::from_blob (const_cast<uint32_t *>(offsets),
                        {K + 1}, kInt32CPU ()) ;

    s->g_keys    [p] = cpu_keys   .to (s->dev).clone () ;
    s->g_offsets [p] = cpu_offsets.to (s->dev).clone () ;

    /* CW records: large — store raw pointer only, do NOT copy now   */
    s->h_cw [p] = cw ;
    s->h_M  [p] = M ;

    std::fprintf (stderr,
                  "[sa_torch] partition %02d registered: K=%ld M=%ld\n",
                  p, K, M) ;

    if (p == s->NN - 1) {
      s->g_ready = true ;
      std::fprintf (stderr,
                    "[sa_torch] all partitions registered, index ready\n") ;
    }
    return 1 ;

  } catch (const std::exception &e) {
    std::fprintf (stderr,
                  "[sa_torch] partition %02d registration failed: %s\n",
                  p, e.what ()) ;
    return 0 ;
  }
}

/* ------------------------------------------------------------------ *
 * saTorchJoinSort  (PATH A)                                          *
 * ------------------------------------------------------------------ */
const uint32_t *saTorchJoinSort (SATorchObj            *tor,
                                 int                    NN,
                                 const uint32_t *const *cws_ptr,
                                 const long            *cws_n,
                                 long                  *sm_count)
{
  *sm_count = 0 ;
  if (!tor) return nullptr ;
  SATorchState *s = static_cast<SATorchState *>(tor) ;
  if (!s->g_ready) return nullptr ;

  try {
    long n = sTorchJoinSortDo (s, NN, cws_ptr, cws_n) ;
    if (n < 0) return nullptr ;
    *sm_count = n ;
    return s->sm_buf.empty () ? nullptr : s->sm_buf.data () ;
  } catch (const std::exception &e) {
    std::fprintf (stderr,
                  "[sa_torch] saTorchJoinSort: %s\n", e.what ()) ;
    return nullptr ;
  }
}

/* ------------------------------------------------------------------ *
 * saTorchCodeJoinSort  (PATH B)                                      *
 * ------------------------------------------------------------------ */
const uint32_t *saTorchCodeJoinSort (SATorchObj     *tor,
                                     const uint8_t  *dna_buf,
                                     long            total_bases,
                                     const uint32_t *off_ptr,
                                     const uint32_t *len_ptr,
                                     int             iaMax,
                                     long           *sm_count)
{
  *sm_count = 0 ;
  if (!tor) return nullptr ;
  SATorchState *s = static_cast<SATorchState *>(tor) ;
  if (!s->g_ready || !s->init_done) return nullptr ;

  try {
    std::vector<std::vector<uint32_t>> cws_out ;
    std::vector<long>                  cws_n_out ;

    if (!sTorchSeedExtract (s, dna_buf, total_bases,
                            off_ptr, len_ptr, iaMax,
                            cws_out, cws_n_out))
      return nullptr ;

    int NN = s->NN ;
    std::vector<const uint32_t *> ptrs ((size_t)NN) ;
    for (int p = 0 ; p < NN ; p++)
      ptrs[(size_t)p] = cws_out[(size_t)p].empty ()
                        ? nullptr
                        : cws_out[(size_t)p].data () ;

    long n = sTorchJoinSortDo (s, NN, ptrs.data (), cws_n_out.data ()) ;
    if (n < 0) return nullptr ;
    *sm_count = n ;
    return s->sm_buf.empty () ? nullptr : s->sm_buf.data () ;

  } catch (const std::exception &e) {
    std::fprintf (stderr,
                  "[sa_torch] saTorchCodeJoinSort: %s\n", e.what ()) ;
    return nullptr ;
  }
}

} /* extern "C" */
