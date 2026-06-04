/*
 * sa_torch.cpp
 *
 * LibTorch GPU backend for Magic2 sort-and-join pipeline.
 *
 * Design rules (mirrors wsra/sra_read.cpp):
 *   - Includes ONLY LibTorch headers and sa.torch.h.
 *   - Never includes sa.h, array.h, ac.h or any acedb header.
 *   - All data crosses the C/C++ boundary as flat pointers + counts.
 *   - Output buffer is owned by SATorchState and loaned to the C caller;
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
 * Genome index memory layout (v2   fully resident):
 *
 *   g_keys[p]     [K_p]     int32  sorted unique seeds   resident on device
 *   g_offsets[p]  [K_p+1]   int32  CSR offset array     resident on device
 *   g_cw[p]       [M_p, 3]  int32  genome payload        resident on device
 *                                   col 0 = nam
 *                                   col 1 = pos
 *                                   col 2 = flags
 *
 * The seed column of the raw CW struct is NOT stored in g_cw.  The CSR
 * index (g_keys + g_offsets) encodes which rows belong to each seed;
 * the seed value itself is never needed again after the join.  Dropping
 * it reduces resident device memory from ~22.8 GB to ~17.1 GB for a
 * human/mouse index (NN=16, 1.43 B records), which fits comfortably in
 * a 32 GB GPU alongside working tensors (~1-2 GB peak per block).
 *
 * Seeds are stored and compared as int32 (signed).  The on-disk index
 * must be built with saRadixSort3Signed so that sort order matches.
 * (LibTorch 2.x has no kUInt32 for bucketize; int32 signed comparison
 * gives identical ordering provided both sides use the same convention.)
 *
 * Output record: TORCHMATCH (5 x uint32 = 20 bytes):
 *   [0] read         read CW nam
 *   [1] x1           read CW pos
 *   [2] target       genome CW nam
 *   [3] a1           genome CW pos
 *   [4] targetFlags  genome CW flags
 *
 * The jumper lines present in cwsN (every row index divisible by 256)
 * are excluded from the Cartesian product by row-index masking inside
 * sTorchJoinSortDo   no extra data field is needed.
 *
 * Authors: Magic2 team / NLM-NIH
 * Created: June 2026
 */

/* LibTorch   must come first                                         */
#include <torch/torch.h>

/* Standard C++ headers                                               */
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <vector>

/* Public C interface                                                  */
#include "sa.torch.h"

/* ================================================================== *
 * §1  Internal state object                                          *
 * ================================================================== */

struct SATorchState {

  /* ---- device ---------------------------------------------------- */
  torch::Device dev ;

  /* ---- genome index: all three tensors resident on device --------- */
  static constexpr int TORCH_NN = 16 ;
  torch::Tensor g_keys    [TORCH_NN] ;  /* [K_p]    int32  sorted seeds  */
  torch::Tensor g_offsets [TORCH_NN] ;  /* [K_p+1]  int32  CSR offsets   */
  torch::Tensor g_cw      [TORCH_NN] ;  /* [M_p, 3] int32  nam,pos,flags */
  bool          g_ready = false ;

  /* ---- seed-packing tables (device tensors) ---------------------- */
  torch::Tensor g_fwd_lut ;     /* [16] int32                           */
  torch::Tensor g_valid_lut ;   /* [16] int32                           */
  torch::Tensor g_fwd_weights ; /* [wLen] int64  4^(wLen-1-j)           */
  torch::Tensor g_rev_weights ; /* [wLen] int64  4^j                    */
  bool          init_done = false ;
  int           wLen = 0 ;
  int           NN   = 0 ;

  /* ---- output buffer (loaned to C caller) ------------------------ */
  /* 5 uint32 per TORCHMATCH record                                   */
  std::vector<uint32_t> tm_buf ;

  /* ---- constructor ---------------------------------------------- */
  explicit SATorchState (torch::Device d) : dev(d) {}
} ;

/* ================================================================== *
 * §2  IUPAC look-up tables                                           *
 * ================================================================== */

static const uint32_t kFwdLutData[16] = {
  0, 0, 3, 0,
  2, 0, 0, 0,
  1, 0, 0, 0,
  0, 0, 0, 0
} ;

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
  return torch::TensorOptions ().dtype (torch::kInt32).device (torch::kCPU) ;
}
static inline torch::TensorOptions kInt64CPU ()
{
  return torch::TensorOptions ().dtype (torch::kInt64).device (torch::kCPU) ;
}
static inline torch::TensorOptions kInt64Dev (const torch::Device &dev)
{
  return torch::TensorOptions ().dtype (torch::kInt64).device (dev) ;
}
static inline torch::TensorOptions kInt32Dev (const torch::Device &dev)
{
  return torch::TensorOptions ().dtype (torch::kInt32).device (dev) ;
}

/* ================================================================== *
 * §4  Device selection                                               *
 * ================================================================== */

static torch::Device sTorchSelectDevice (int device_id)
{
  if (torch::cuda::is_available ()) {
    int n = torch::cuda::device_count () ;
    if (device_id < 0 || device_id >= n) device_id = 0 ;
    std::fprintf (stderr, "[sa_torch] device: CUDA %d\n", device_id) ;
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
 * Precondition: s->g_ready == true, all g_keys/g_offsets/g_cw resident.
 *
 * For each partition p:
 *   cws_ptr[p] / cws_n[p]    read-side CW arrays (4 cols, int32).
 *   g_keys[p], g_offsets[p], g_cw[p]    resident genome tensors.
 *
 * Jumper lines (row index % 256 == 0) are excluded by masking.
 *
 * On success: fills s->tm_buf with TORCHMATCH records (5 x uint32 each),
 * sorted by read nam (col 0).  Returns total record count.
 * On failure: returns -1.
 */
static long sTorchJoinSortDo (SATorchState           *s,
                               int                     NN,
                               const uint32_t *const  *cws_ptr,
                               const long             *cws_n)
{
  std::vector<torch::Tensor> tm_parts ;

  for (int p = 0 ; p < NN ; p++) {

    if (!cws_ptr[p] || cws_n[p] <= 0) continue ;
    if (!s->g_cw[p].defined () || s->g_cw[p].size (0) == 0) continue ;

    long N = cws_n[p] ;

    /* 1. Wrap read CWs as CPU tensor (4 cols), move to device       */
    torch::Tensor R =
      torch::from_blob (const_cast<uint32_t *>(cws_ptr[p]),
                        {N, 4}, kInt32CPU ())
      .to (s->dev) ;                             /* [N, 4] int32    */

    /* 2. Extract seed column (col 0) as int32   signed comparison.
     * Seeds on both sides are stored as int32 (saRadixSort3Signed),
     * so bucketize with g_keys (also int32) gives correct order.    */
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
      genome_seeds.eq (r_seeds) ;               /* [N] bool        */

    torch::Tensor r_matched =
      exact_match.nonzero ().squeeze (1) ;       /* [n_match] int64 */

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

    /* 6. Expand read rows by hit count (Cartesian product LHS)     */
    /* Keep only cols 1,2 (nam, pos) from read  flags always zero  */
    torch::Tensor R_matched =
      R.index_select (0, r_matched)
       .narrow (1, 1, 2) ;                       /* [n_match, 2]    */

    torch::Tensor R_exp =
      R_matched.repeat_interleave (counts, 0) ;  /* [total, 2]      */

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

    int64_t M_p = s->g_cw[p].size (0) ;
    torch::Tensor g_rows =
      (group_start + inner_offset)
      .clamp (0, M_p - 1) ;                      /* [total] int64   */

    /* 8. Filter jumper lines: exclude row index % 256 == 0         */
    torch::Tensor jumper_keep =
      (g_rows & 255).ne (0).nonzero ().squeeze (1) ;

    if (jumper_keep.numel () == 0) continue ;

    g_rows = g_rows.index_select (0, jumper_keep) ; /* [n_keep]     */
    R_exp  = R_exp .index_select (0, jumper_keep) ; /* [n_keep, 2]  */

    /* 9. Gather genome payload from resident g_cw (3 cols)         */
    torch::Tensor G_exp =
      s->g_cw[p].index_select (0, g_rows) ;     /* [n_keep, 3]     */
    /* cols: 0=nam, 1=pos, 2=flags                                  */

    /* 10. Assemble TORCHMATCH rows (5 columns)
     *   [0] read  = R_exp col 0  (read nam)
     *   [1] x1    = R_exp col 1  (read pos)
     *   [2] target= G_exp col 0  (genome nam)
     *   [3] a1    = G_exp col 1  (genome pos)
     *   [4] targetFlags = G_exp col 2                              */
    torch::Tensor TM = torch::cat ({
        R_exp.narrow (1, 0, 1),                  /* [0] read        */
        R_exp.narrow (1, 1, 1),                  /* [1] x1          */
        G_exp.narrow (1, 0, 1),                  /* [2] target      */
        G_exp.narrow (1, 1, 1),                  /* [3] a1          */
        G_exp.narrow (1, 2, 1),                  /* [4] targetFlags */
    }, /*dim=*/1) ;                              /* [n_keep, 5]     */

    tm_parts.push_back (TM) ;
  }

  if (tm_parts.empty ()) {
    s->tm_buf.clear () ;
    return 0 ;
  }

  /* 11. Concatenate all partitions                                  */
  torch::Tensor TM_all =
    (tm_parts.size () == 1)
    ? tm_parts[0]
    : torch::cat (tm_parts, 0) ;

  tm_parts.clear () ;

  /* 12. Sort by read nam (column 0)                                 */
  torch::Tensor read_col =
    TM_all.select (1, 0).contiguous () ;
  torch::Tensor sort_idx = torch::argsort (read_col) ;
  torch::Tensor TM_sorted =
    TM_all.index_select (0, sort_idx) ;

  /* 13. Copy back to host into tm_buf                              */
  torch::Tensor TM_cpu =
    TM_sorted.to (torch::kCPU).contiguous () ;

  int64_t Total = TM_cpu.size (0) ;
  s->tm_buf.resize ((size_t)(Total * 5)) ;
  std::memcpy (s->tm_buf.data (),
               TM_cpu.data_ptr<int32_t> (),
               (size_t)Total * 5 * sizeof (uint32_t)) ;

  return (long)Total ;
}

/* ================================================================== *
 * §6  GPU seed extraction (PATH B)                                   *
 * ================================================================== */

/*
 * sTorchSeedExtract   unchanged from v1.
 * Fills cws_out[p] and cws_n_out[p] with read CW records (4 cols).
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

  torch::Tensor raw     = raw_cpu.to (s->dev) ;
  torch::Tensor nibbles = raw.to (torch::kInt32) & 0x0F ;

  torch::Tensor bases2  = s->g_fwd_lut  .index_select (0, nibbles) ;
  torch::Tensor valid1  = s->g_valid_lut.index_select (0, nibbles) ;

  /* 2. Upload offset and length arrays                              */
  torch::Tensor seq_off =
    torch::from_blob (const_cast<uint32_t *>(off_ptr),
                      {(long)iaMax}, kInt32CPU ())
    .to (s->dev) ;

  torch::Tensor seq_len =
    torch::from_blob (const_cast<uint32_t *>(len_ptr),
                      {(long)iaMax}, kInt32CPU ())
    .to (s->dev) ;

  /* 3. Number of windows per read                                   */
  torch::Tensor seq_len64 = seq_len.to (torch::kInt64) ;
  torch::Tensor n_win =
    (seq_len64 - wLen + 1).clamp_min (0) ;

  int64_t N_win_total = n_win.sum ().item<int64_t> () ;
  if (N_win_total == 0) {
    cws_out.assign (NN, {}) ;
    cws_n_out.assign (NN, 0) ;
    return true ;
  }

  /* 4. Build flat window-start indices                              */
  torch::Tensor seq_off64 = seq_off.to (torch::kInt64) ;

  torch::Tensor read_base =
    seq_off64.repeat_interleave (n_win) ;

  torch::Tensor ia_idx =
    torch::arange (0, (int64_t)iaMax, kInt64Dev (s->dev)) ;

  torch::Tensor ia_exp =
    ia_idx.repeat_interleave (n_win) ;

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

  torch::Tensor win_start = read_base + inner ;

  /* 5. Gather k-base windows                                        */
  torch::Tensor offsets_w =
    torch::arange ((int64_t)wLen, kInt64Dev (s->dev)) ;

  torch::Tensor win_idx =
    win_start.unsqueeze (1) + offsets_w ;
  win_idx = win_idx.clamp (0, total_bases - 1) ;

  torch::Tensor flat_idx = win_idx.reshape ({-1}) ;

  torch::Tensor W =
    bases2.index_select (0, flat_idx)
          .reshape ({N_win_total, (int64_t)wLen}) ;

  torch::Tensor V =
    valid1.index_select (0, flat_idx)
          .reshape ({N_win_total, (int64_t)wLen}) ;

  /* 6. Validity mask                                                */
  torch::Tensor all_valid =
    std::get<0> (V.min (1)).eq (1) ;

  /* 7. Pack forward and RC k-mer values                             */
  torch::Tensor W64 = W.to (torch::kInt64) ;

  torch::Tensor S = (W64 * s->g_fwd_weights).sum (1) ;

  torch::Tensor Wc = 3 - W64 ;
  torch::Tensor C  = (Wc * s->g_rev_weights).sum (1) ;

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

  /* Store seed as int32 (signed)   matches saRadixSort3Signed order */
  torch::Tensor z_stored =
    (torch::bitwise_right_shift (Z, nHidden2) & (int64_t)0xFFFFFFFF)
    .to (torch::kInt32) ;

  /* 9. Assemble CW columns                                          */
  torch::Tensor strand_bit = minus.to (torch::kInt32) ;

  /* ia_exp is 0-based; nam = ((ia+1) << 1) | strand_bit            */
  torch::Tensor nam =
    torch::bitwise_left_shift ((ia_exp).to (torch::kInt32), 1)
    | strand_bit ;

  torch::Tensor pos   = (inner + 1).to (torch::kInt32) ;
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
    }, /*dim=*/1) ;                              /* [N_p, 4] int32  */

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
                          kInt64CPU ())
        .to (dev).clone () ;
    }

    /* Reverse weight vector: rev_weights[j] = 4^j                  */
    {
      std::vector<int64_t> rw ((size_t)wLen) ;
      int64_t pw = 1 ;
      for (int j = 0 ; j < wLen ; j++, pw *= 4) rw[(size_t)j] = pw ;
      s->g_rev_weights =
        torch::from_blob (rw.data (), {(long)wLen},
                          kInt64CPU ())
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
  if (tor) delete static_cast<SATorchState *>(tor) ;
  return nullptr ;
}

/* ------------------------------------------------------------------ *
 * saTorchIsAvailable                                                 *
 * ------------------------------------------------------------------ */
int saTorchIsAvailable (SATorchObj *tor)
{
  if (!tor) return 0 ;
  return (static_cast<SATorchState *>(tor)->dev.type () != torch::kCPU) ? 1 : 0 ;
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
    std::fprintf (stderr, "[sa_torch] --testTorch: %s\n", e.what ()) ;
    return 1 ;
  }
}

/* ------------------------------------------------------------------ *
 * saTorchIndexUpload                                                 *
 *                                                                    *
 * Upload strategy (v2):                                              *
 *   g_keys     int32 on device (seeds stored signed, same as index) *
 *   g_offsets   int32 on device                                      *
 *   g_cw       3-column int32 on device (nam, pos, flags)           *
 *              seed column dropped; saves ~25% device memory.       *
 *   All three are cloned to device immediately; the C caller may    *
 *   free its source buffers as soon as this function returns.       *
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
    /* g_keys: store as int32 (signed)   matches saRadixSort3Signed  */
    s->g_keys[p] =
      torch::from_blob (const_cast<uint32_t *>(seeds),
                        {K}, kInt32CPU ())
      .to (s->dev).clone () ;

    /* g_offsets: int32                                              */
    s->g_offsets[p] =
      torch::from_blob (const_cast<uint32_t *>(offsets),
                        {K + 1}, kInt32CPU ())
      .to (s->dev).clone () ;

    /* g_cw: 3 columns only (nam=col1, pos=col2, flags=col3).
     * The raw cw array is [M, 4] uint32; we view it as such,
     * slice off columns 1-3, and clone to device.
     * from_blob sees the full 4-col layout; narrow selects cols 1..3. */
    torch::Tensor cw_cpu =
      torch::from_blob (const_cast<uint32_t *>(cw),
                        {M, 4}, kInt32CPU ()) ;

    s->g_cw[p] =
      cw_cpu.narrow (1, 1, 3)   /* cols 1,2,3: nam,pos,flags       */
            .contiguous ()
            .to (s->dev).clone () ;

    std::fprintf (stderr,
                  "[sa_torch] partition %02d uploaded: K=%ld M=%ld"
                  " (%.2f GB on device)\n",
                  p, K, M,
                  (double)(K * 4 + (K + 1) * 4 + M * 12) / 1e9) ;

    if (p == s->NN - 1) {
      s->g_ready = true ;
      std::fprintf (stderr,
                    "[sa_torch] all partitions ready, index resident\n") ;
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
const uint32_t *saTorchJoinSort (SATorchObj            *tor,
                                 int                    NN,
                                 const uint32_t *const *cws_ptr,
                                 const long            *cws_n,
                                 long                  *tm_count)
{
  *tm_count = 0 ;
  if (!tor) return nullptr ;
  SATorchState *s = static_cast<SATorchState *>(tor) ;
  if (!s->g_ready) return nullptr ;

  try {
    long n = sTorchJoinSortDo (s, NN, cws_ptr, cws_n) ;
    if (n < 0) return nullptr ;
    *tm_count = n ;
    return s->tm_buf.empty () ? nullptr : s->tm_buf.data () ;
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
                                     long           *tm_count)
{
  *tm_count = 0 ;
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
    *tm_count = n ;
    return s->tm_buf.empty () ? nullptr : s->tm_buf.data () ;

  } catch (const std::exception &e) {
    std::fprintf (stderr,
                  "[sa_torch] saTorchCodeJoinSort: %s\n", e.what ()) ;
    return nullptr ;
  }
}

} /* extern "C" */
