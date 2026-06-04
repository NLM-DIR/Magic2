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
 * Genome index tensors (uploaded once per run by saTorchIndexUpload):
 *   g_keys[p]     [K_p]     uint32   sorted unique seeds for partition p
 *   g_offsets[p]  [K_p+1]   uint32   CSR offset array
 *   g_cw[p]       [M_p, 4]  uint32   genome CW records
 *
 * CW field layout (4 x uint32, matches sa.h CW struct):
 *   [0] seed    sort/lookup key
 *   [1] nam     read id: (ia << 1) | strand_bit
 *   [2] pos     1-based coordinate of first base of seed in read
 *   [3] flags   auxiliary
 *
 * SEEDMATCH layout (8 x uint32, matches sa.h SEEDMATCH struct):
 *   [0] readSeed   = 0 placeholder
 *   [1] read       = read CW nam
 *   [2] x1         = read CW pos
 *   [3] readFlags  = read CW flags
 *   [4] targetSeed = 0 placeholder
 *   [5] target     = genome CW nam
 *   [6] a1         = genome CW pos
 *   [7] targetFlags= genome CW flags
 *
 * Seed encoding (matches saCodeSequenceSeedsStep1 in sa.seeds.c):
 *   2-bit fwd: A=0 C=1 G=2 T=3
 *   Z = min(fwd_kmer, rc_kmer)  strand-agnostic
 *   strand_bit = 1 if rc < fwd
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

  /* ---- genome index (device tensors, one per partition) ---------- */
  static constexpr int TORCH_NN = 16 ;
  torch::Tensor g_keys    [TORCH_NN] ;   /* [K_p]     uint32 */
  torch::Tensor g_offsets [TORCH_NN] ;   /* [K_p+1]   uint32 */
  torch::Tensor g_cw      [TORCH_NN] ;   /* [M_p, 4]  uint32 */
  bool          g_ready = false ;

  /* ---- seed-packing tables (device tensors) --------------------- */
  torch::Tensor g_fwd_lut ;      /* [16] uint32 */
  torch::Tensor g_valid_lut ;    /* [16] uint32 */
  torch::Tensor g_fwd_weights ;  /* [wLen] uint64  4^(wLen-1-j) */
  torch::Tensor g_rev_weights ;  /* [wLen] uint64  4^j           */
  bool          init_done = false ;
  int           wLen = 0 ;
  int           NN   = 0 ;

  /* ---- output buffer (loaned to C caller, never freed by C) ----- */
  std::vector<uint32_t> sm_buf ;   /* holds last SEEDMATCH output */

  /* ---- constructor ---------------------------------------------- */
  explicit SATorchState (torch::Device d) : dev(d) {}
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

/* 4-bit IUPAC nibble → 1 if pure ACGT, else 0                      */
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
 * For partition p in [0..NN-1]:
 *   cws_ptr[p]  read CW records as uint32_t*, 4 per CW, or NULL
 *   cws_n[p]    number of CW records
 *
 * On success: fills s->sm_buf with SEEDMATCH records (8 x uint32 each),
 * returns the number of records.
 * On failure: returns -1.
 */
static long sTorchJoinSortDo (SATorchState           *s,
                               int                     NN,
                               const uint32_t *const   *cws_ptr,
                               const long             *cws_n)
{
  std::vector<torch::Tensor> sm_parts ;

  for (int p = 0 ; p < NN ; p++) {

    if (!cws_ptr[p] || cws_n[p] <= 0) continue ;

    long N = cws_n[p] ;

    /* 1. Wrap read CWs as CPU tensor, move to device               */
    torch::Tensor R =
      torch::from_blob (const_cast<uint32_t *>(cws_ptr[p]),
                        {N, 4}, kInt32CPU ())
      .to (s->dev) ;                             /* [N, 4] uint32   */

    /* 2. Extract seed column                                        */
    torch::Tensor r_seeds = R.select (1, 0).contiguous () ; /* [N] */

    /* 3. Bucketize: floor position in sorted genome keys            */
    torch::Tensor buckets =
      torch::bucketize (r_seeds, s->g_keys[p],
                        /*out_uint32=*/false,
                        /*right=*/false) ;         /* [N] int64     */

    /* 4. Verify exact match                                         */
    int64_t K_p = s->g_keys[p].size (0) ;
    torch::Tensor b_clamped = buckets.clamp (0, K_p - 1) ;

    torch::Tensor genome_seeds =
      s->g_keys[p].index_select (0, b_clamped) ;

    torch::Tensor exact_match = genome_seeds.eq (r_seeds) ; /* bool */

    torch::Tensor r_matched =
      exact_match.nonzero ().squeeze (1) ;         /* [n_match] int64 */

    if (r_matched.numel () == 0) continue ;

    /* 5. Gather CSR starts and counts                               */
    torch::Tensor b_matched =
      b_clamped.index_select (0, r_matched) ;

    torch::Tensor starts =
      s->g_offsets[p]
        .index_select (0, b_matched)
        .to (torch::kInt64) ;

    torch::Tensor counts =
      (s->g_offsets[p]
         .index_select (0, b_matched + 1)
         .to (torch::kInt64)
       - starts)
      .clamp_min (0) ;

    int64_t total = counts.sum ().item<int64_t> () ;
    if (total == 0) continue ;

    /* 6. Expand read rows by hit count                              */
    torch::Tensor R_matched =
      R.index_select (0, r_matched) ;

    torch::Tensor R_exp =
      R_matched.repeat_interleave (counts, 0) ;    /* [total, 4]    */

    /* 7. Build genome row index vector                              */
    torch::Tensor group_start =
      starts.repeat_interleave (counts) ;

    torch::Tensor cum = counts.cumsum (0) ;
    torch::Tensor group_base =
      torch::cat ({
          torch::zeros ({1}, kInt64Dev (s->dev)),
          cum.slice (0, 0, -1)
      }) ;

    torch::Tensor group_base_exp =
      group_base.repeat_interleave (counts) ;

    torch::Tensor inner_offset =
      torch::arange (total, kInt64Dev (s->dev)) - group_base_exp ;

    torch::Tensor g_rows =
      (group_start + inner_offset)
      .clamp (0, s->g_cw[p].size (0) - 1) ;

    /* 8. Gather genome CWs and assemble SEEDMATCH                  */
    torch::Tensor G_exp =
      s->g_cw[p].index_select (0, g_rows) ;       /* [total, 4]    */

    torch::Tensor zeros =
      torch::zeros ({total, 1}, kInt32Dev (s->dev)) ;

    /* SEEDMATCH: 8 columns                                         */
    torch::Tensor SM = torch::cat ({
        zeros,                                      /* [0] readSeed  */
        R_exp.narrow (1, 1, 1),                     /* [1] read nam  */
        R_exp.narrow (1, 2, 1),                     /* [2] x1        */
        R_exp.narrow (1, 3, 1),                     /* [3] readFlags */
        zeros,                                      /* [4] tgtSeed   */
        G_exp.narrow (1, 1, 1),                     /* [5] target    */
        G_exp.narrow (1, 2, 1),                     /* [6] a1        */
        G_exp.narrow (1, 3, 1),                     /* [7] tgtFlags  */
    }, /*dim=*/1) ;                                 /* [total, 8]    */

    sm_parts.push_back (SM) ;
  }

  if (sm_parts.empty ()) {
    s->sm_buf.clear () ;
    return 0 ;
  }

  /* 9. Concatenate all partitions                                   */
  torch::Tensor SM_all =
    (sm_parts.size () == 1)
    ? sm_parts[0]
    : torch::cat (sm_parts, 0) ;

  sm_parts.clear () ;

  /* 10. Sort by read nam (column 1)                                 */
  torch::Tensor read_col =
    SM_all.select (1, 1).contiguous () ;
  torch::Tensor sort_idx = torch::argsort (read_col) ;
  torch::Tensor SM_sorted =
    SM_all.index_select (0, sort_idx) ;

  /* 11. Copy back to host into sm_buf                              */
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
 * Fills cws_out[p] and cws_n_out[p] for each partition.
 *
 * On return the caller owns the vectors; they are passed straight into
 * sTorchJoinSortDo via pointer arrays built on the stack.
 *
 * Returns true on success.
 */
static bool sTorchSeedExtract (
    SATorchState                  *s,
    const uint8_t                 *dna_buf,
    long                           total_bases,
    const uint32_t                 *off_ptr,
    const uint32_t                 *len_ptr,
    int                            iaMax,
    std::vector<std::vector<uint32_t>> &cws_out,   /* [NN][N_p * 4] */
    std::vector<long>              &cws_n_out)     /* [NN]          */
{
  const int wLen    = s->wLen ;
  const int NN      = s->NN ;
  const int nHidden = wLen > 16 ? wLen - 16 : 0 ;
  const int nHidden2= nHidden * 2 ;

  /* iStep fixed at 1 on GPU — see architecture doc §6              */
  /* (CPU min-hash phasing has no GPU equivalent; stride >1 would   *
   *  create blind spots)                                           */

  if (total_bases <= 0 || iaMax <= 1) return false ;

  /* 1. Upload flat DNA buffer                                       */
  torch::Tensor raw_cpu =
    torch::from_blob (const_cast<uint8_t *>(dna_buf),
                      {total_bases},
                      torch::TensorOptions ()
                        .dtype (torch::kUInt8)
                        .device (torch::kCPU)) ;

  torch::Tensor raw = raw_cpu.to (s->dev) ;        /* [total_bases] */

  torch::Tensor nibbles =
    raw.to (torch::kInt32) & 0x0F ;

  torch::Tensor bases2 =
    s->g_fwd_lut.index_select (0, nibbles) ;       /* [total] uint32 */

  torch::Tensor valid1 =
    s->g_valid_lut.index_select (0, nibbles) ;     /* [total] uint32 */

  /* 2. Upload pre-built offset and length arrays                   */
  torch::Tensor seq_off =
    torch::from_blob (const_cast<uint32_t *>(off_ptr),
                      {(long)iaMax}, kInt32CPU ())
    .to (s->dev) ;                                 /* [iaMax] uint32 */

  torch::Tensor seq_len =
    torch::from_blob (const_cast<uint32_t *>(len_ptr),
                      {(long)iaMax}, kInt32CPU ())
    .to (s->dev) ;                                 /* [iaMax] uint32 */

  /* 3. Number of windows per read                                   */
  torch::Tensor seq_len64 = seq_len.to (torch::kInt64) ;
  torch::Tensor n_win =
    (seq_len64 - wLen + 1).clamp_min (0) ;         /* [iaMax] int64 */

  int64_t N_win_total = n_win.sum ().item<int64_t> () ;
  if (N_win_total == 0) {
    cws_out.assign (NN, {}) ;
    cws_n_out.assign (NN, 0) ;
    return true ;
  }

  /* 4. Build flat window-start indices                              */
  torch::Tensor seq_off64 = seq_off.to (torch::kInt64) ;

  torch::Tensor read_base =
    seq_off64.repeat_interleave (n_win) ;          /* [N_win] int64 */

  torch::Tensor ia_idx =
    torch::arange (0, (int64_t)iaMax, kInt64Dev (s->dev)) ;
  torch::Tensor ia_exp =
    ia_idx.repeat_interleave (n_win) ;             /* [N_win] int64 */

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
    - group_base_exp ;

  torch::Tensor win_start = read_base + inner ;    /* [N_win] int64 */

  /* 5. Gather k-base windows                                        */
  torch::Tensor offsets_w =
    torch::arange ((int64_t)wLen, kInt64Dev (s->dev)) ;

  torch::Tensor win_idx =
    win_start.unsqueeze (1) + offsets_w ;          /* [N_win, wLen] */

  win_idx = win_idx.clamp (0, total_bases - 1) ;

  torch::Tensor flat_idx = win_idx.reshape ({-1}) ;

  torch::Tensor W =
    bases2.index_select (0, flat_idx)
          .reshape ({N_win_total, (int64_t)wLen}) ;  /* uint32        */

  torch::Tensor V =
    valid1.index_select (0, flat_idx)
          .reshape ({N_win_total, (int64_t)wLen}) ;

  /* 6. Validity mask                                                */
  torch::Tensor all_valid =
    std::get<0> (V.min (1)).eq (1) ;
 
  /* 7. Pack forward and RC seeds                                    */
  torch::Tensor W64 = W.to (torch::kInt64) ;

  torch::Tensor S =
    (W64 * s->g_fwd_weights).sum (1) ;             /* [N_win] int64 */

  torch::Tensor Wc = 3 - W64 ;
  torch::Tensor C =
    (Wc * s->g_rev_weights).sum (1) ;

  int64_t maskSeedLn = (wLen < 32)
    ? ((int64_t)1 << (2 * wLen)) - 1
    : (int64_t)0x7FFFFFFFFFFFFFFFLL ;
  S = S & maskSeedLn ;
  C = C & maskSeedLn ;

  /* 8. Strand-agnostic Z, partition, stored seed                   */
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
  torch::Tensor nam =
    torch::bitwise_left_shift (ia_exp.to (torch::kInt32), 1) | strand_bit ;

  torch::Tensor pos =
    (inner + 1).to (torch::kInt32) ;              /* 1-based         */

  torch::Tensor flags =
    torch::zeros ({N_win_total}, kInt32Dev (s->dev)) ;

  /* 10. Apply validity filter                                       */
  torch::Tensor keep =
    all_valid.nonzero ().squeeze (1) ;

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

  /* 11. Scatter into NN partitions                                  */
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
    }, /*dim=*/1) ;                               /* [N_p, 4] uint32 */

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
      int64_t p = 1 ;
      for (int j = wLen - 1 ; j >= 0 ; j--, p *= 4) fw[(size_t)j] = p ;
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
      int64_t p = 1 ;
      for (int j = 0 ; j < wLen ; j++, p *= 4) rw[(size_t)j] = p ;
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
 * ------------------------------------------------------------------ */
int saTorchIndexUpload (SATorchObj     *tor,
                        int             p,
                        const uint32_t *seeds,    long K,
                        const uint32_t *offsets,
                        const uint32_t  *cw,       long M)
{
  if (!tor) return 0 ;
  SATorchState *s = static_cast<SATorchState *>(tor) ;

  try {
    torch::Tensor cpu_keys =
      torch::from_blob (const_cast<uint32_t *>(seeds),
                        {K}, kInt32CPU ()) ;
    torch::Tensor cpu_offsets =
      torch::from_blob (const_cast<uint32_t *>(offsets),
                        {K + 1}, kInt32CPU ()) ;
    torch::Tensor cpu_cw =
      torch::from_blob (const_cast<uint32_t *>(cw),
                        {M, 4}, kInt32CPU ()) ;

    s->g_keys    [p] = cpu_keys   .to (s->dev) ;
    s->g_offsets [p] = cpu_offsets.to (s->dev) ;
    s->g_cw      [p] = cpu_cw     .to (s->dev) ;

    std::fprintf (stderr,
                  "[sa_torch] partition %02d uploaded: K=%ld M=%ld\n",
                  p, K, M) ;

    if (p == s->NN - 1) {
      s->g_ready = true ;
      std::fprintf (stderr,
                    "[sa_torch] genome index ready on device\n") ;
    }
    return 1 ;

  } catch (const std::exception &e) {
    std::fprintf (stderr,
                  "[sa_torch] partition %02d upload failed: %s\n",
                  p, e.what ()) ;
    return 0 ;
  }
}

/* ------------------------------------------------------------------ *
 * saTorchJoinSort  (PATH A)                                          *
 * ------------------------------------------------------------------ */
const uint32_t *saTorchJoinSort (SATorchObj          *tor,
                                int                  NN,
                                const uint32_t *const *cws_ptr,
                                const long           *cws_n,
                                long                 *sm_count)
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
const uint32_t *saTorchCodeJoinSort (SATorchObj    *tor,
                                    const uint8_t *dna_buf,
                                    long           total_bases,
                                    const uint32_t *off_ptr,
                                    const uint32_t *len_ptr,
                                    int            iaMax,
                                    long          *sm_count)
{
  *sm_count = 0 ;
  if (!tor) return nullptr ;
  SATorchState *s = static_cast<SATorchState *>(tor) ;
  if (!s->g_ready || !s->init_done) return nullptr ;

  try {
    std::vector<std::vector<uint32_t>> cws_out ;
    std::vector<long>                 cws_n_out ;

    if (!sTorchSeedExtract (s, dna_buf, total_bases,
                            off_ptr, len_ptr, iaMax,
                            cws_out, cws_n_out))
      return nullptr ;

    /* Build pointer arrays on the stack for sTorchJoinSortDo       */
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
