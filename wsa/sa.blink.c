/*
 * sa.blink.c

 * This module is part of the magic2 package
 * Created May 22, 2026

 * Authors: Danielle Thierry-Mieg, Jean Thierry-Mieg, Greg Boratyn, NCBI/NLM/NIH

 * This code is public.
 * 
 * This module count words in reads and detects SNPs
 */

#include "sa.h"
#include <stdint.h>


/********************************************************************************************************************/
/*************************** construct combs of alternate parity of length KMAX in a dna ****************************/
/********************************************************************************************************************/
/* ln must be odd and se select the strand based on the central letter */
static kmerCountOne (Array dna, Array kmers, int ln)
{
  int i, iMax = dna ? arrayMax (dna) : 0 ;
  int iMax2 ;
  const unsigned char *cp = iMax ? arrayp (dna, 0, unsigned char) ;
  unsigned long w0 = 0, w1 = 0, w1r = 0 ;
  int dn = 2 * (ln - 1) ;
  unsigned long mask = ((unsigned long) 1) << (2 * ln) - 1 ;
  unsigned long mask2 = ((unsigned long) 1) << (2 * ln) - 3 ; 
  if (cp)
    {
      for (i = 0, iMax2 = iMax - 2 * ln + 1 ; i < iMax2 ; i++, cp++)
	{
	  w0 = ((w0 << 2) | (cp[0] & 0x3)) & mask ;
	  w1 = ((w1 << 2) | (cp[1] & 0x3)) & mask ;
	  w0r = ((w0r >> 2) & mask2) | (((cp[1] ^ 0x3) & 0x3) << dn) ;
	  w1r = ((w1r >> 2) & mask2) | (((cp[1] ^ 0x3) & 0x3) << dn) ;
	}
      if (i >= ln)
	{ /* words are complete */

	}
    }
} /* kmerCount */

  static kmerCount (PP *pp, BB *bb, Array kmers)
/*
 * sa.kmer_bitmap.c
 *
 * K-mer existence bitmap: each agent accumulates up to NMAX k-mer indices,
 * radix-sorts them, then mutex-locks the global 2GB bitmap and writes all
 * bits sequentially (cache-line by cache-line) before releasing the lock.
 *
 * Design rationale
 * ----------------
 *   WL  = word length in bases (k).  Each base = 2 bits.
 *   k   = 17  =>  4^17 = 2^34  possible k-mers.
 *   Bitmap size = 2^34 bits = 2 GB  (1 << 31 bytes).
 *   One cache line = 64 bytes = 512 bits.
 *   Number of cache lines in bitmap = 2^34 / 512 = 2^25 = 33,554,432.
 *
 *   Each agent collects up to NMAX k-mer indices (uint64_t).
 *   After a block, it radix-sorts them (2-pass, 17-bit radix, see below)
 *   then acquires the global mutex and writes sequentially, one cache line
 *   at a time, using 64-bit OR into the bitmap.  Because the list is sorted
 *   the writes are strictly ascending in memory: the hardware prefetcher
 *   loves this and cache-line eviction/reload is minimised.
 *
 * Radix sort
 * ----------
 *   2^34 values fit in 34 bits.  We use a 2-pass radix sort:
 *     pass 1: low  17 bits  (buckets 0..131071)
 *     pass 2: high 17 bits  (buckets 0..131071)
 *   Each pass is O(N).  Count arrays are 128K entries of uint32_t = 512 KB,
 *   fits comfortably in L2.  Scratch buffer of NMAX uint64_t reused each
 *   agent block (allocated once per agent at startup).
 *
 * AVX-512 / AVX2 notes
 * --------------------
 *   The bitmap write loop uses plain 64-bit stores; on x86-64 these are
 *   naturally 8-byte aligned if the bitmap is aligned (see posix_memalign).
 *   For the count/scatter passes we use plain C; the compiler with
 *   -O3 -march=native will auto-vectorise the histogram accumulation.
 *   If you want explicit AVX-512 scatter you can replace the scatter loop
 *   with _mm512_i64scatter_epi64, but in practice the bottleneck is RAM
 *   bandwidth, not the scatter instruction itself.
 *
 * Compile flags expected:
 *   -O3 -march=native -std=c11 -pthread
 *
 * Authors: Magic2 team, NIH/NLM
 */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <assert.h>

/* ------------------------------------------------------------------ */
/*  Parameters                                                          */
/* ------------------------------------------------------------------ */

#define WL        17                      /* k-mer word length in bases */
#define KMER_BITS (2 * WL)               /* 34 bits to represent a k-mer */
#define KMER_N    (1ULL << KMER_BITS)    /* 4^17 = 2^34 possible k-mers  */

#define BITMAP_BYTES  (KMER_N / 8)       /* 2 GB                          */
#define CACHE_LINE    64                  /* bytes; universally true x86-64 & ARM64 */
#define CL_BITS       (CACHE_LINE * 8)   /* 512 bits per cache line        */
#define CL_UINT64     (CACHE_LINE / 8)   /* 8 x uint64_t per cache line    */
#define N_CACHE_LINES (BITMAP_BYTES / CACHE_LINE)  /* 2^25 = 33,554,432   */

/* Radix parameters: 2 passes of 17 bits each cover 34-bit k-mer index */
#define RADIX_BITS  17
#define RADIX_SIZE  (1 << RADIX_BITS)    /* 131072 buckets                 */
#define RADIX_MASK  (RADIX_SIZE - 1)

/* Maximum k-mers accumulated per agent before flush.
   Tune between 3M and 32M.  32M gives ~12.5% cache-line collision
   probability per write => ~87.5% of mutex acquisitions are uncontested. */
#define NMAX  (32 * 1024 * 1024)         /* 32M k-mer slots per agent      */


/* ------------------------------------------------------------------ */
/*  Global bitmap + mutex                                               */
/* ------------------------------------------------------------------ */

static uint8_t  *g_bitmap;              /* 2 GB, 64-byte aligned           */
static pthread_mutex_t g_bitmap_mutex = PTHREAD_MUTEX_INITIALIZER;


/* ------------------------------------------------------------------ */
/*  Per-agent working storage                                           */
/* ------------------------------------------------------------------ */

typedef struct {
    uint64_t *kmers;    /* accumulation buffer, NMAX entries             */
    uint64_t *scratch;  /* scratch buffer for radix sort, NMAX entries   */
    uint32_t *cnt0;     /* radix histogram pass 0, RADIX_SIZE entries    */
    uint32_t *cnt1;     /* radix histogram pass 1, RADIX_SIZE entries    */
    int       n;        /* number of valid entries currently in kmers[]  */
} AgentBuf;


/* ------------------------------------------------------------------ */
/*  Initialisation (call once at program start)                         */
/* ------------------------------------------------------------------ */

/* Allocate the 2 GB global bitmap, 64-byte aligned so every cache line
   is naturally aligned.  Returns 0 on success. */
int kmer_bitmap_init (void)
{
    int rc = posix_memalign ((void **)&g_bitmap, CACHE_LINE, BITMAP_BYTES);
    if (rc) return rc;
    memset (g_bitmap, 0, BITMAP_BYTES);
    return 0;
}

/* Allocate per-agent buffers.  Call once per agent thread at startup.
   Returns 0 on success. */
int agent_buf_init (AgentBuf *ab)
{
    /* k-mer arrays: 8 bytes each */
    ab->kmers   = malloc (NMAX * sizeof(uint64_t));
    ab->scratch = malloc (NMAX * sizeof(uint64_t));
    ab->cnt0    = malloc (RADIX_SIZE * sizeof(uint32_t));
    ab->cnt1    = malloc (RADIX_SIZE * sizeof(uint32_t));
    ab->n       = 0;
    if (!ab->kmers || !ab->scratch || !ab->cnt0 || !ab->cnt1)
        return -1;
    return 0;
}

void agent_buf_destroy (AgentBuf *ab)
{
    free (ab->kmers);   ab->kmers   = NULL;
    free (ab->scratch); ab->scratch = NULL;
    free (ab->cnt0);    ab->cnt0    = NULL;
    free (ab->cnt1);    ab->cnt1    = NULL;
    ab->n = 0;
}


/* ------------------------------------------------------------------ */
/*  K-mer extraction  (inline, called millions of times)               */
/* ------------------------------------------------------------------ */

/* 2-bit encoding: A=0, C=1, G=2, T=3.
   Adjust to match Magic2's existing encoding in sa.seeds.c            */
static const uint8_t base2bits[256] = {
    /* fill everything with 0, override the 4 bases */
    ['A'] = 0, ['a'] = 0,
    ['C'] = 1, ['c'] = 1,
    ['G'] = 2, ['g'] = 2,
    ['T'] = 3, ['t'] = 3,
    ['N'] = 0, ['n'] = 0,   /* N maps to A; k-mer is still recorded */
};

/* Push a k-mer into the agent buffer.  Flush automatically when full. */
/* Declaration forward; flush is defined below. */
static void agent_flush (AgentBuf *ab);

static inline void agent_push_kmer (AgentBuf *ab, uint64_t kmer)
{
    ab->kmers[ab->n++] = kmer;
    if (__builtin_expect (ab->n == NMAX, 0))
        agent_flush (ab);
}

/* Slide a WL-base window along 'seq' of length 'seqlen', pushing each
   valid k-mer into the agent buffer.
   'N' bases reset the window (k-mer spanning an N is skipped).        */
void agent_extract_kmers (AgentBuf *ab,
                          const unsigned char *seq, int seqlen)
{
    uint64_t w    = 0;           /* current k-mer shift register         */
    uint64_t mask = KMER_N - 1; /* low KMER_BITS bits                   */
    int      fill = 0;           /* how many valid bases in w            */

    for (int i = 0; i < seqlen; i++) {
        uint8_t b = base2bits[(unsigned char)seq[i]];
        if (seq[i] == 'N' || seq[i] == 'n') {
            fill = 0; w = 0;    /* restart after N                      */
            continue;
        }
        w = ((w << 2) | b) & mask;
        if (++fill >= WL)
            agent_push_kmer (ab, w);
    }
}


/* ------------------------------------------------------------------ */
/*  Radix sort  (2-pass, 17 + 17 bits, fully in-place via scratch)    */
/* ------------------------------------------------------------------ */

static void radix_sort_kmers (uint64_t *src, uint64_t *dst,
                              uint32_t *cnt0, uint32_t *cnt1,
                              int n)
{
    /* -- histogram -- */
    memset (cnt0, 0, RADIX_SIZE * sizeof(uint32_t));
    memset (cnt1, 0, RADIX_SIZE * sizeof(uint32_t));

    for (int i = 0; i < n; i++) {
        cnt0[  src[i]               & RADIX_MASK ]++;
        cnt1[ (src[i] >> RADIX_BITS) & RADIX_MASK ]++;
    }

    /* -- prefix sums -- */
    uint32_t s0 = 0, s1 = 0;
    for (int b = 0; b < RADIX_SIZE; b++) {
        uint32_t t0 = cnt0[b]; cnt0[b] = s0; s0 += t0;
        uint32_t t1 = cnt1[b]; cnt1[b] = s1; s1 += t1;
    }

    /* -- pass 1: scatter by low 17 bits, src -> dst -- */
    for (int i = 0; i < n; i++) {
        uint32_t bucket = src[i] & RADIX_MASK;
        dst[cnt0[bucket]++] = src[i];
    }

    /* -- pass 2: scatter by high 17 bits, dst -> src -- */
    for (int i = 0; i < n; i++) {
        uint32_t bucket = (dst[i] >> RADIX_BITS) & RADIX_MASK;
        src[cnt1[bucket]++] = dst[i];
    }
    /* result is back in src, fully sorted ascending */
}


/* ------------------------------------------------------------------ */
/*  Bitmap write  (called under mutex, sequential cache-line writes)   */
/* ------------------------------------------------------------------ */

/* Write all k-mers in sorted array 'kmers[0..n-1]' into g_bitmap.
   Because the array is sorted we walk it once, accumulating bits for
   the current 64-byte cache line, then store them with a single 8-way
   OR (8 x uint64_t) before advancing to the next cache line.
   This keeps write traffic sequential and minimises cache misses.     */
static void bitmap_write_sorted (const uint64_t *kmers, int n)
{
    if (n == 0) return;

    uint64_t *bm = (uint64_t *)g_bitmap;  /* treat bitmap as uint64_t array */
    /* Each uint64_t covers 64 bits.  A cache line covers 8 uint64_t.
       Cache-line index of a k-mer: kmer / CL_BITS = kmer >> 9           */

    int i = 0;
    while (i < n) {
        /* cache line that kmer[i] belongs to */
        uint64_t cl_idx = kmers[i] >> 9;   /* divide by 512 bits */
        uint64_t *cl    = bm + cl_idx * CL_UINT64;

        /* accumulate a local copy of this cache line */
        uint64_t tmp[CL_UINT64];
        memcpy (tmp, cl, CACHE_LINE);      /* one 64-byte load             */

        /* OR in all k-mers that fall in this same cache line */
        do {
            uint64_t bit_pos = kmers[i] & (CL_BITS - 1); /* bit within CL */
            tmp[bit_pos >> 6] |= (1ULL << (bit_pos & 63));
            i++;
        } while (i < n && (kmers[i] >> 9) == cl_idx);

        memcpy (cl, tmp, CACHE_LINE);      /* one 64-byte store            */
    }
}


/* ------------------------------------------------------------------ */
/*  Flush: sort + lock + write + unlock                                 */
/* ------------------------------------------------------------------ */

static void agent_flush (AgentBuf *ab)
{
    if (ab->n == 0) return;

    /* --- sort (parallel, no lock needed) --- */
    radix_sort_kmers (ab->kmers, ab->scratch,
                      ab->cnt0, ab->cnt1, ab->n);

    /* --- write under mutex --- */
    pthread_mutex_lock   (&g_bitmap_mutex);
    bitmap_write_sorted  (ab->kmers, ab->n);
    pthread_mutex_unlock (&g_bitmap_mutex);

    ab->n = 0;
}


/* ------------------------------------------------------------------ */
/*  Public API summary                                                  */
/* ------------------------------------------------------------------ */
/*
 *  Startup (main thread):
 *      kmer_bitmap_init();
 *
 *  Each agent thread:
 *      AgentBuf ab;
 *      agent_buf_init(&ab);
 *      ...
 *      agent_extract_kmers(&ab, seq, len);   // called per read
 *      ...
 *      agent_flush(&ab);                     // call at end of block too
 *      agent_buf_destroy(&ab);
 *
 *  Query (after all agents done, no locking needed):
 *      static inline int kmer_exists(uint64_t kmer) {
 *          return (g_bitmap[kmer >> 3] >> (kmer & 7)) & 1;
 *      }
 */
