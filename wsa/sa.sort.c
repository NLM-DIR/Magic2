/*
 * sa.sort.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module implements a sort-merge algorithm
 * and ends on the correct p[arity in an insertion sort to avoid a global copy
 * If the table with n lines is already sorted,
 * the code just performs n-1 comparisons
 * The order functions are inlined to acccelarate the system
 * If the hardware allows it, 128bit copies are performed
 * Notice that all our structures are 128 bits aligned
 */

#include "sa.h"

#ifdef __SSE2__
#define VECTORIZED_MEM_CPY
#include <emmintrin.h> // SSE2
#include <immintrin.h>
#endif /* __SSE2__ */

#ifdef TIME_EVAL      
#include <time.h>
#endif

#ifdef USEGPU
#include "sa.gpusort.h"
#endif /* USEGPU */

/**************************************************************/
/**************************************************************/
/* sort treating the full struct as a 256 bits unsigned int */

/**************************************************************/
/**************************************************************/

static inline int cwOrder (const void *va, const void *vb)
{
  const CW *up = va ;
  const CW *vp = vb ;
  int n ;
  n = (up->seed > vp->seed) - (up->seed < vp->seed) ; if (n) return n ;

  return 0 ;
} /* cwOrder */

/**************************************************************/
/* up->read is the second field in the struct which contains 8 unsigned ints */
static inline int seedMatchOrder (const void *va, const void *vb)
{
  const SEEDMATCH *up = va ;
  const SEEDMATCH *vp = vb ;
  int n ;
  n = (up->read > vp->read) - (up->read < vp->read) ; if (n) return n ;

  return 0 ;
} /* seedMatchOrder */

/**************************************************************/
/* a0 = a1 - x1 is the putative position of base 1 of the read 
 * It also works for the negative strand (a1 < 0, x1 > 0).
 */
static inline int hitReadOrder (const void *va, const void *vb)
{
  const HIT *up = va ;
  const HIT *vp = vb ;
  int n ;

  n = ((up->read > vp->read) - (up->read < vp->read)) ; if (n) return n ;
  n = ((up->chrom > vp->chrom) - (up->chrom < vp->chrom)) ; if (n) return n ; 
  n = ((up->a1 > vp->a1) - (up->a1 < vp->a1)) ;  if (n) return n ;
  n = ((up->x1 > vp->x1) - (up->x1 < vp->x1)) ; 
  return n ;
} /* hitReadOrder */

/**************************************************************/
/* a0 = a1 - x1 is the putative position of base 1 of the read 
 * It also works for the negative strand (a1 < 0, x1 > 0).
 */
static inline int hitPairOrder (const void *va, const void *vb)
{
  const HIT *up = va ;
  const HIT *vp = vb ;
  int n ;
  /*
   * this line will be removed after a check on all test datasets, it is supposed to never occur 
   *  if (up->a1 < 0) messcrash ("negative a1") ;
   */
  n = ((up->read >> 1) > (vp->read >> 1)) -  ((up->read >> 1) < (vp->read >> 1)) ; if (n) return n ; 
  n = ((up->chrom > vp->chrom) - (up->chrom < vp->chrom)) ; if (n) return n ; 
  long unsigned int n1 = (long unsigned int) up->a1 + (up->x1 >> NSHIFTEDTARGETREPEATBITS) ;
  long unsigned int n2 = (long unsigned int) vp->a1 + (vp->x1 >> NSHIFTEDTARGETREPEATBITS) ;
  n = (n1 > n2) - (n1 < n2) ; if (n) return n ;
  n = ((up->x1 > vp->x1) - (up->x1 < vp->x1)) ; 

  return n ;
} /* hitPairOrder */

/**************************************************************/
/**************************************************************/
/* cwOrder specialized version, sort of first uint out of 4 */
static BOOL checkClean1 (const char *b, mysize_t n)
{
  const unsigned int *up  = (const unsigned int *) b ;
  const unsigned int *end = up + 4 * (n - 1) ;

  while (up < end)
    {
      if (*(up + 4) < *up)
        return FALSE ;
      up += 4 ;
    }
  return TRUE ;
} /* checkClean1 */

/**************************************************************/
/* hitReadOrder specialized version */
static BOOL checkClean2 (const char *b, mysize_t n)
{
  const HIT *up  = (const HIT *) b ;
  const HIT *end = up + (n - 1) ;

  while (up < end)
    {
      if (hitReadOrder (up+1, up) < 0)
        return FALSE ;
      up++ ;
    }
  return TRUE ;
} /* checkClean2 */

/**************************************************************/
/* hitPairOrder specialized version */
static BOOL checkClean3 (const char *b, mysize_t n)
{
  const HIT *up  = (const HIT *) b ;
  const HIT *end = up + (n - 1) ;

  while (up < end)
    {
      if (hitPairOrder (up+1, up) < 0)
        return FALSE ;
      up++ ;
    }
  return TRUE ;
} /* checkClean3 */

/**************************************************************/
/* seedMatchOrder specialized version, sort on second uint out of 8 */
static BOOL checkClean4 (const char *b, mysize_t n)
{
  const unsigned int *up = (const unsigned int *) b + 1 ;  /* start at field 1 */
  const unsigned int *end = up + 8 * (n - 1) ;
  while (up < end)
    {
      if (*(up + 8) < *up)
	return FALSE ;
      up += 8 ;
    }
  return TRUE ;
} /* checkClean4 */

/**************************************************************/
/**************************************************************/
/* saSort algorithm minimizing memcpy */

/* Taquin insertion algorithm
 * works en place
 */
static BOOL newInsertionSort (char *b, mysize_t n, int s, int (*cmp)(const void *va, const void *vb))
{
  mysize_t i, j ;
  char buf[s] ;
  BOOL clean = TRUE ;

  for (i = 1 ; i < n ; i++)
    {
      j = i - 1 ;
      if ((*cmp) (b + i*s, b + j*s) >= 0)
	continue ;
      clean = FALSE ;
      memcpy (buf, b + i*s, s) ;
      memcpy (b + i*s, b + j*s, s) ;
      while (j > 0 && (*cmp) (buf, b + (j-1)*s) < 0)
	{
	  memcpy (b + j*s, b + (j-1)*s, s) ;
	  j-- ;
	}
      memcpy (b + j*s, buf, s) ;
    }
  return clean ;
} /* newInsertionSort */

/* #ifndef USEGPU */
/* recursivelly split the table
 * the partially sorted data oscillate between b and buf
 * they end up correctly in b because for small n
 * we switch to the insertion taquin algoright
 * on correct parity, as speed is 2 persent higher with
 * insertion n>0,  relative to n==0
 * but n=8,16,32 are equivalent speeds
 */

/* __attribute__((target("avx2"))) */
static BOOL saSortDo (char *b, long int nn, int s, char *buf, BOOL hitIsTarget, int (*cmp)(const void *va, const void *vb))
{
 char *up, *vp, *wp ;
  long int n1 = nn / 2 ;
  long int n2 = nn - n1 ;
  char *b1 = b ;
  char *b2 = b + n1 * s ;
  char *b01 = buf ;
  char *b02 = buf + n1 * s ;
  BOOL ok = FALSE ;
  BOOL clean1, clean2, clean = TRUE ;
  /* for small n,
   * sort en place using the insertion algorithm (game of taquin)
   */
  if (hitIsTarget && nn <= 16)
    {
      clean = newInsertionSort (b, nn, s, cmp) ;
      return clean ;
    }

  /* otherwise: sort the 2 halves exchanging hit and buf */
  clean1 = saSortDo (b01, n1, s, b1, ! hitIsTarget, cmp) ;
  clean2 = saSortDo (b02, n2, s, b2, ! hitIsTarget, cmp) ;
  
  /* then merge the 2 sorted halves */
  up = b01 ;
  vp = b02 ;
  wp = b1 ;

  if  ((*cmp) (b02 - s, b02) <= 0)
    {
      /* sortmerge is not needed, copy whole blocks */
      /* do we need to copy */
      if (! clean1)
	{ /* copy n1 (or nn idf ! clean2)  records back to b */
	  clean = FALSE ;
	  memcpy(wp, up, (clean2 ? n1 : nn)  * s) ;	  
	}
      else if (! clean2)
	{ /* just copy n2 to the second part of b */
	  clean = FALSE ;
	  memcpy(wp + n1 * s, vp, n2 * s);	  
	}
      /* if clean1 && clean2, no copying is needed */
    return clean  ;
    }
  clean = FALSE ;
  
#ifdef VECTORIZED_MEM_CPY
  if (cmp == cwOrder)
    {
      while (n1 > 0 && n2 > 0)
	{
	  __m128i u = _mm_load_si128((__m128i*)up) ;
	  __m128i v = _mm_load_si128((__m128i*)vp) ;

	  int n = (*(unsigned int*)up <= *(unsigned int*)vp) ? 1 : 0 ;
	  
	  _mm_store_si128((__m128i*)wp, n  ? u : v) ;
	  wp += s ;
	  up = n ? up + s : up ;
	  vp = n ? vp : vp + s ;
	  n1 -= n ;
	  n2 -= ! n ;
	}
      ok = TRUE ;
    }
  
  else if (cmp == seedMatchOrder) 
    {
      while (n1 > 0 && n2 > 0)
	{
	  __m128i u1 = _mm_load_si128((__m128i*)up) ;
	  __m128i v1 = _mm_load_si128((__m128i*)vp) ;
	  __m128i u2 = _mm_load_si128((__m128i*)up+1) ;
	  __m128i v2 = _mm_load_si128((__m128i*)vp+1) ;
	  
	  int n = (*((unsigned int*)up + 1) <= *((unsigned int*)vp + 1)) ? 1 : 0 ;
	  
	  _mm_store_si128((__m128i*)wp, n  ? u1 : v1) ;
	  _mm_store_si128((__m128i*)wp+1, n  ? u2 : v2) ;
	  wp += s ;
	  up = n ? up + s : up ;
	  vp = n ? vp : vp + s ;
	  n1 -= n ;
	  n2 -= ! n ;
	}
      ok = TRUE;
    }

#endif  

  if (! ok) /* either we do not have _mm_store_si128, or size s is not 16 */
    { /* classic code */
      while (n1 > 0 && n2 > 0)
	{
	  int n = ((*cmp) (up, vp) <= 0) ? 1 : 0 ;

	  memcpy (wp, (n ? up : vp), s) ;
	  wp += s ;
	  up = n ? up + s : up ;
	  vp = n ? vp : vp + s ;
	  n1 -= n ;
	  n2 -= ! n ;
	}
    }

  /* I also tried to count all greater cases and bulk copy
   * but this code was more complex and not faster
   */

  /* bulk copy the reminders */
    if (n1 > 0) memcpy(wp, up, n1 * s);
    else if (n2 > 0) memcpy(wp, vp, n2 * s);

    return clean ;
} /* saSortDo */
/* #endif // USEGPU */

/**************************************************************/
static int saRadixSort (void *base, size_t N, size_t stride, int keyIndex) ;
static void *saRadixSort3 (void *base, size_t N, size_t stride, int keyIndex) ;
int saSort (BigArray aa, int type)
{
  long int N = bigArrayMax (aa) ;
  char *cp = N ?  (char *) aa->base : 0 ;
  int s = aa->size ;
  int (*cmp)(const void *va, const void *vb) = NULL ;
  int usedGPU = 0 ;
  int useRadix = 1 ;
  BOOL done = FALSE ;
  
  if (N < 2) return FALSE ;   /* i did not use the GPU */

  switch (type)
    {
    case 1:
      if (s != 16) messcrash ("Wrong call to saSort type 1") ;
      if (checkClean1 (cp, N))
	return FALSE ;  /* i did not use the GPU */
      cmp = cwOrder ;
      break ;
    case 2:
      if (s != sizeof(HIT)) messcrash ("Wrong call to saSort type 2") ;      
      if (checkClean2 (cp, N))
	return FALSE ;  /* i did not use the GPU */
      cmp = hitReadOrder ;
      break ;
    case 3:
      if (s != sizeof(HIT)) messcrash ("Wrong call to saSort type 3") ;      
      if (checkClean3 (cp, N))
	return FALSE ;  /* i did not use the GPU */
      cmp = hitPairOrder ;
      break ;
    case 4:
      if (s != 32)
	messcrash ("Wrong call to saSort type 4") ;
      if (checkClean4 (cp, N))
	return FALSE ;  /* i did not use the GPU */
      cmp = seedMatchOrder ;
      break ;
    case 5:
      if (s != 32)
	messcrash ("Wrong call to saSort type 5") ;
      break ;
    default:
      messcrash ("Wrong call to saSort type = %d > 4", type) ;
    }
  if (N < 128)
    {
      newInsertionSort (cp, N, s, cmp) ;
      done = TRUE ;
    }

#ifdef TIME_EVAL
  struct timespec start, end;
  double ellapsed;
  timespec_get(&start, TIME_UTC);
#endif
  
#ifdef USEGPU
  if (! done && type < 3 && (size_t)N * s > (1 << 20))
    {
      usedGPU = 1 ;
      saGPUSort (cp, N, type) ;
      done = TRUE ;
    }
#endif

  /* GPU threshold not met: fall through to CPU sort below */
  if (! done && useRadix)
    {
      void *dst = 0 ;
      switch (type)
	{
	case 1:
	  dst = saRadixSort3 (bigArrp (aa, 0, CW), bigArrayMax (aa), sizeof (CW), 0) ;
	  bigArraySwitchBase (aa, bigArrayMax (aa), dst) ;
	  done = TRUE ;
	  break ;
	case 4:
	  dst = saRadixSort3 (bigArrp (aa, 0, SEEDMATCH), bigArrayMax (aa), sizeof (SEEDMATCH), 1) ;
	  bigArraySwitchBase (aa, bigArrayMax (aa), dst) ;
	  done = TRUE ;
	  break ;
	}
    }

  if (! done)
    {
      size_t alloc_size = ((size_t)N * s + 15) & ~15 ;
      char *buf = aligned_alloc(16, alloc_size) ;
      if (! buf) messcrash ("\nsa.sort.c alloc failure, consider lowering --bMax\n") ;
      memcpy (buf, cp, N * s) ;
      saSortDo (cp, N, s, buf, TRUE, cmp) ;
      free (buf) ;
    }
  
#ifdef TIME_EVAL      
  timespec_get(&end, TIME_UTC);
  ellapsed = (double)(end.tv_sec - start.tv_sec) +
    (double)(end.tv_nsec - start.tv_nsec) / 1000000000.0;
  fprintf(stderr, "Sorted %ld elements in %f seconds (type: %d)\n", N, ellapsed, type);
#endif      

  return usedGPU ;
}/* saSort */


/**************************************************************/
/**************************************************************/
/**************************************************************/

/* Two-pass LSB radix sort on an array of fixed-size structs.
 *
 * Prototype:
 *   saRadixSort(void *base, size_t N, size_t stride, int keyIndex)
 *
 *   base      : pointer to the array, guaranteed 128-byte aligned
 *   N         : number of records
 *   stride    : size in bytes of each struct (must be a multiple of
 *               sizeof(unsigned int) and >= 2*sizeof(unsigned int))
 *   keyIndex  : 0 => sort on first  unsigned int (k1, offset 0)
 *               1 => sort on second unsigned int (k2, offset 4)
 *
 * The sort is a stable, two-pass LSB radix sort, 16 bits per pass.
 *   Pass 1: bits  0..15 (low  half of key)
 *   Pass 2: bits 16..31 (high half of key)
 *
 * The result is left in base[].
 * Returns 0 on success, -1 on allocation failure.
 *
 * Fast paths (compiled only when VECTORIZED_MEM_CPY is defined):
 *   stride == 16 : one _mm_load/store_si128  per record (1 x 128-bit op)
 *   stride == 32 : two _mm_load/store_si128  per record (2 x 128-bit ops)
 *   Requires SSE2 only (-msse2), available on all x86-64 targets.
 *   Key extraction uses a direct unsigned int* read from the aligned
 *   pointer — no SSE4.1 (_mm_extract_epi32) needed.
 *
 *  N==0 guard made explicit (bigArrayMax can return 0 with null aa).
 *  All VLAs removed; counts/offsets are static-sized on the stack.
 **************************************************************/

#define RADIX_BITS  16
#define RADIX_SIZE  (1u << RADIX_BITS)   /* 65 536 buckets */
#define RADIX_MASK  (RADIX_SIZE - 1u)

#ifdef VECTORIZED_MEM_CPY

/* Read the sort key from an aligned struct pointer.
 * keyOff is 0 (k1) or 4 (k2).
 * Casting through unsigned int* is valid in C: the struct is aligned
 * and unsigned char* -> unsigned int* aliasing is well-defined here
 * because the underlying object really is an unsigned int.            */
#define SA_KEY_FROM_PTR(ptr)  (*(unsigned int *)((ptr) + keyOff))

#endif  /* VECTORIZED_MEM_CPY */

static int saRadixSort (void *base, size_t N, size_t stride, int keyIndex)
{
  unsigned char       *src = (unsigned char *) base ;
  unsigned char       *dst ;
  long int             counts[2][RADIX_SIZE] ;
  long int             offsets[RADIX_SIZE] ;
  long int             i ;
  unsigned int         key ;
  const size_t         keyOff = (size_t) keyIndex * sizeof (unsigned int) ;

  if (!base || N <= 1)
    return 0 ;

  /* scratch buffer: N records of stride bytes, 128-byte aligned */
  dst = (unsigned char *) aligned_alloc (128, N * stride) ;
  if (!dst)
    messcrash ("malloc failure in saRadixSort") ;

  /* ------------------------------------------------------------------ */
  /* Histogram pass: count low-16 and high-16 bits in one linear scan   */
  /* ------------------------------------------------------------------ */
  memset (counts, 0, sizeof (counts)) ;
  for (i = 0 ; i < (long int) N ; i++)
    {
      memcpy (&key, src + (size_t) i * stride + keyOff, sizeof (unsigned int)) ;
      counts[0][ key         & RADIX_MASK]++ ;   /* low  16 bits */
      counts[1][(key >> 16u) & RADIX_MASK]++ ;   /* high 16 bits */
    }

  /* ------------------------------------------------------------------ */
  /* Pass 1: sort on low 16 bits of key,  src -> dst                    */
  /* ------------------------------------------------------------------ */
  offsets[0] = 0 ;
  for (i = 1 ; i < (long int) RADIX_SIZE ; i++)
    offsets[i] = offsets[i - 1] + counts[0][i - 1] ;

#ifdef VECTORIZED_MEM_CPY

  if (stride == 16)
    {
      for (i = 0 ; i < (long int) N ; i++)
        {
          unsigned char  *sp     = src + (size_t) i * stride ;
          unsigned int    bucket = SA_KEY_FROM_PTR (sp) & RADIX_MASK ;
          __m128i         rec    = _mm_load_si128 ((__m128i *) sp) ;
          _mm_store_si128 ((__m128i *)(dst + (size_t) offsets[bucket] * stride), rec) ;
          offsets[bucket]++ ;
        }
    }
  else if (stride == 32)
    {
      for (i = 0 ; i < (long int) N ; i++)
        {
          unsigned char  *sp     = src + (size_t) i * stride ;
          unsigned int    bucket = SA_KEY_FROM_PTR (sp) & RADIX_MASK ;
          __m128i         lo     = _mm_load_si128 ((__m128i *) sp) ;
          __m128i         hi     = _mm_load_si128 ((__m128i *)(sp + 16)) ;
          unsigned char  *dp     = dst + (size_t) offsets[bucket] * stride ;
          _mm_store_si128 ((__m128i *) dp,       lo) ;
          _mm_store_si128 ((__m128i *)(dp + 16), hi) ;
          offsets[bucket]++ ;
        }
    }
  else
#endif  /* VECTORIZED_MEM_CPY */
    {
      for (i = 0 ; i < (long int) N ; i++)
        {
          unsigned int bucket ;
          memcpy (&key, src + (size_t) i * stride + keyOff, sizeof (unsigned int)) ;
          bucket = key & RADIX_MASK ;
          memcpy (dst + (size_t) offsets[bucket] * stride,
                  src + (size_t) i              * stride,
                  stride) ;
          offsets[bucket]++ ;
        }
    }

  /* ------------------------------------------------------------------ */
  /* Pass 2: sort on high 16 bits of key, dst -> src                    */
  /* ------------------------------------------------------------------ */
  offsets[0] = 0 ;
  for (i = 1 ; i < (long int) RADIX_SIZE ; i++)
    offsets[i] = offsets[i - 1] + counts[1][i - 1] ;

#ifdef VECTORIZED_MEM_CPY
  if (stride == 16)
    {
      for (i = 0 ; i < (long int) N ; i++)
        {
          unsigned char  *dp     = dst + (size_t) i * stride ;
          unsigned int    bucket = (SA_KEY_FROM_PTR (dp) >> 16u) & RADIX_MASK ;
          __m128i         rec    = _mm_load_si128 ((__m128i *) dp) ;
          _mm_store_si128 ((__m128i *)(src + (size_t) offsets[bucket] * stride), rec) ;
          offsets[bucket]++ ;
        }
    }
  else if (stride == 32)
    {
      for (i = 0 ; i < (long int) N ; i++)
        {
          unsigned char  *dp     = dst + (size_t) i * stride ;
          unsigned int    bucket = (SA_KEY_FROM_PTR (dp) >> 16u) & RADIX_MASK ;
          __m128i         lo     = _mm_load_si128 ((__m128i *) dp) ;
          __m128i         hi     = _mm_load_si128 ((__m128i *)(dp + 16)) ;
          unsigned char  *sp     = src + (size_t) offsets[bucket] * stride ;
          _mm_store_si128 ((__m128i *) sp,       lo) ;
          _mm_store_si128 ((__m128i *)(sp + 16), hi) ;
          offsets[bucket]++ ;
        }
    }
  else
#endif  /* VECTORIZED_MEM_CPY */
    {
      for (i = 0 ; i < (long int) N ; i++)
        {
          unsigned int bucket ;
          memcpy (&key, dst + (size_t) i * stride + keyOff, sizeof (unsigned int)) ;
          bucket = (key >> 16u) & RADIX_MASK ;
          memcpy (src + (size_t) offsets[bucket] * stride,
                  dst + (size_t) i              * stride,
                  stride) ;
          offsets[bucket]++ ;
        }
    }

  /* result is back in src (== base); dst was scratch only */
  free (dst) ;
  return 0 ;
}  /* saRadixSort */

/**************************************************************/
/**************************************************************/
/**************************************************************/

/**************************************************************/
/**************************************************************/
/**************************************************************/

/* Three-pass LSB radix sort on an array of fixed-size structs.
 *
 * Prototype:
 *   saRadixSort3(void *base, size_t N, size_t stride, int keyIndex)
 *
 *   base      : pointer to the array, guaranteed 128-byte aligned
 *   N         : number of records
 *   stride    : size in bytes of each struct (must be a multiple of
 *               sizeof(unsigned int) and >= 2*sizeof(unsigned int))
 *   keyIndex  : 0 => sort on first  unsigned int (k1, offset 0)
 *               1 => sort on second unsigned int (k2, offset 4)
 *
 * The sort is a stable, three-pass LSB radix sort:
 *   Pass 1: bits  0..10  (11 bits, 2048 buckets)
 *   Pass 2: bits 11..21  (11 bits, 2048 buckets)
 *   Pass 3: bits 22..31  (10 bits, 1024 buckets)
 *
 * Working set per pass: 2048 * sizeof(long int) = 16 kB — fits in L1.
 * This eliminates the LLC cache-miss thrashing that made the 2-pass
 * 65536-bucket version memory-bound.
 *
 * The result is left in base[].
 * After an odd number of passes the final result must be copied back;
 * this is done with the same SIMD path used for scatter, so the extra
 * copy is cheap.
 *
 * Fast paths (compiled only when VECTORIZED_MEM_CPY is defined):
 *   stride == 16 : one _mm_load/store_si128  per record
 *   stride == 32 : two _mm_load/store_si128  per record
 *   Requires SSE2 only (-msse2).
 *
 * Returns 0 on success, -1 on allocation failure.
 **************************************************************/

#define R3_BITS_A   11
#define R3_BITS_B   11
#define R3_BITS_C   10

#define R3_SIZE_A   (1u << R3_BITS_A)   /* 2048 */
#define R3_SIZE_B   (1u << R3_BITS_B)   /* 2048 */
#define R3_SIZE_C   (1u << R3_BITS_C)   /* 1024 */

#define R3_MASK_A   (R3_SIZE_A - 1u)
#define R3_MASK_B   (R3_SIZE_B - 1u)
#define R3_MASK_C   (R3_SIZE_C - 1u)

#define R3_SHIFT_A   0
#define R3_SHIFT_B   R3_BITS_A                     /* 11 */
#define R3_SHIFT_C   (R3_BITS_A + R3_BITS_B)       /* 22 */

#ifdef VECTORIZED_MEM_CPY
#ifndef _EMMINTRIN_H_INCLUDED   /* guard against double-include */
#include <emmintrin.h>
#endif

#define SA3_KEY_FROM_PTR(ptr)  (*(unsigned int *)((ptr) + keyOff))

/* scatter one 16-byte record from sp into dp */
#define SA3_SCATTER_16(dp, sp)                                  \
  do {                                                          \
    __m128i _r = _mm_load_si128  ((__m128i *)(sp)) ;           \
    _mm_store_si128 ((__m128i *)(dp), _r) ;                     \
  } while (0)

/* scatter one 32-byte record from sp into dp */
#define SA3_SCATTER_32(dp, sp)                                  \
  do {                                                          \
    __m128i _lo = _mm_load_si128 ((__m128i *)(sp)) ;           \
    __m128i _hi = _mm_load_si128 ((__m128i *)((sp) + 16)) ;    \
    _mm_store_si128 ((__m128i *)(dp),        _lo) ;             \
    _mm_store_si128 ((__m128i *)((dp) + 16), _hi) ;             \
  } while (0)

#endif  /* VECTORIZED_MEM_CPY */

/* ------------------------------------------------------------------ */
/* Internal helper: one scatter pass                                   */
/*   from[0..N-1] -> to[0..N-1]                                       */
/*   key bits selected by (key >> shift) & mask                        */
/*   counts[] is the already-computed histogram for this pass          */
/* ------------------------------------------------------------------ */
static void saR3Pass (unsigned char *from,
                      unsigned char *to,
                      size_t         N,
                      size_t         stride,
                      size_t         keyOff,
                      unsigned int   shift,
                      unsigned int   mask,
                      long int       counts[],   /* R3_SIZE_A >= all sizes */
                      long int       offsets[])
{
  long int      i ;
  unsigned int  sz = mask + 1u ;   /* number of buckets for this pass */

  /* prefix sum */
  offsets[0] = 0 ;
  for (i = 1 ; i < (long int) sz ; i++)
    offsets[i] = offsets[i - 1] + counts[i - 1] ;

  /* scatter */
#ifdef VECTORIZED_MEM_CPY
  if (stride == 16)
    {
      for (i = 0 ; i < (long int) N ; i++)
        {
          unsigned char *sp     = from + (size_t) i * stride ;
          unsigned int   bucket = (SA3_KEY_FROM_PTR (sp) >> shift) & mask ;
          SA3_SCATTER_16 (to + (size_t) offsets[bucket] * stride, sp) ;
          offsets[bucket]++ ;
        }
      return ;
    }
  if (stride == 32)
    {
      for (i = 0 ; i < (long int) N ; i++)
        {
          unsigned char *sp     = from + (size_t) i * stride ;
          unsigned int   bucket = (SA3_KEY_FROM_PTR (sp) >> shift) & mask ;
          SA3_SCATTER_32 (to + (size_t) offsets[bucket] * stride, sp) ;
          offsets[bucket]++ ;
        }
      return ;
    }
#endif  /* VECTORIZED_MEM_CPY */
  /* generic fallback */
  {
    unsigned int key ;
    for (i = 0 ; i < (long int) N ; i++)
      {
        unsigned int bucket ;
        memcpy (&key, from + (size_t) i * stride + keyOff, sizeof (unsigned int)) ;
        bucket = (key >> shift) & mask ;
        memcpy (to   + (size_t) offsets[bucket] * stride,
                from + (size_t) i              * stride,
                stride) ;
        offsets[bucket]++ ;
      }
  }
}  /* saR3Pass */

/* ------------------------------------------------------------------ */
/* Copy-back helper: dst -> src, full array, same SIMD paths           */
/* ------------------------------------------------------------------ */
static void saR3CopyBack (unsigned char *dst,
                          unsigned char *src,
                          size_t         N,
                          size_t         stride)
{
#ifdef VECTORIZED_MEM_CPY
  size_t i ;
  if (stride == 16)
    {
      for (i = 0 ; i < N ; i++)
        SA3_SCATTER_16 (src + i * stride, dst + i * stride) ;
      return ;
    }
  if (stride == 32)
    {
      for (i = 0 ; i < N ; i++)
        SA3_SCATTER_32 (src + i * stride, dst + i * stride) ;
      return ;
    }
#endif
  memcpy (src, dst, N * stride) ;
}  /* saR3CopyBack */

/* ------------------------------------------------------------------ */
/* Main entry point                                                     */
/* ------------------------------------------------------------------ */
static void *saRadixSort3 (void *base, size_t N, size_t stride, int keyIndex)
{
  unsigned char  *src = (unsigned char *) base ;
  unsigned char  *dst ;
  /* largest pass has R3_SIZE_A == 2048 buckets */
  long int        counts[3][R3_SIZE_A] ;
  long int        offsets[R3_SIZE_A] ;
  long int        i ;
  unsigned int    key ;
  const size_t    keyOff = (size_t) keyIndex * sizeof (unsigned int) ;

  if (!base || N <= 1)
    return 0 ;

  /* scratch buffer: N records, 128-byte aligned */
  dst = (unsigned char *) aligned_alloc (128, N * stride) ;
  if (!dst)
    messcrash ("malloc failure in saRadixSort3") ;

  /* ------------------------------------------------------------------ */
  /* Histogram pass: fill all three count arrays in one linear scan     */
  /* ------------------------------------------------------------------ */
  memset (counts, 0, sizeof (counts)) ;
  for (i = 0 ; i < (long int) N ; i++)
    {
      memcpy (&key, src + (size_t) i * stride + keyOff, sizeof (unsigned int)) ;
      counts[0][(key >> R3_SHIFT_A) & R3_MASK_A]++ ;
      counts[1][(key >> R3_SHIFT_B) & R3_MASK_B]++ ;
      counts[2][(key >> R3_SHIFT_C) & R3_MASK_C]++ ;
    }

  /* ------------------------------------------------------------------ */
  /* Pass 1: bits  0..10,  src -> dst                                   */
  /* Pass 2: bits 11..21,  dst -> src                                   */
  /* Pass 3: bits 22..31,  src -> dst                                   */
  /* After pass 3 the sorted data is in dst; copy back to src (base).   */
  /* ------------------------------------------------------------------ */
  saR3Pass (src, dst, N, stride, keyOff, R3_SHIFT_A, R3_MASK_A, counts[0], offsets) ;
  saR3Pass (dst, src, N, stride, keyOff, R3_SHIFT_B, R3_MASK_B, counts[1], offsets) ;
  saR3Pass (src, dst, N, stride, keyOff, R3_SHIFT_C, R3_MASK_C, counts[2], offsets) ;
  if (0) saR3CopyBack (dst, src, N, stride) ;

  return dst  ;
}  /* saRadixSort3 */

/**************************************************************/
/**************************************************************/
/**************************************************************/



#ifdef JUNK
# AI prompt to analyze the code

This code module is part of the sortalign package (RNA aligner, NCBI/NLM/NIH).
  Key context:
  - C code targeting portability, SSE2 available, AVX2 avoided
  - BigArray, HIT, CW, SEEDMATCH, BOOL, messcrash() are defined in sa.h
  - Arrays are elastic tables, possibly reallocated when accessed with a large i index as bigArrayp(a,i,type)
     but accessed via macro as bigArrp(a,i,type) )(checked in debugged mode, unchecked in production mode)
     - DNA is coded on the 4 low bits of an unsigned char  A_=0x1, T_=0x2, C_=0x4, G_=Ox8
  
  This code runs on very large dataset (Tera bases of dna) and will be distributed
  and recompiled in many places
  
  The purpose of this chat is to analyze this module against any code weakness
  and suggesting any speed optimization

  the function i wish to analyze  loops on an array of dnas.
       for each one we wish to extract a seed of length seedLength (18)
       and store the best canditate among the next 'step' seeds
       also the best between the seed w and its complement wr.

	 this is performed on the reads and on the target genome
	 a joint between these 2 tables will recover the coordinates of the matching seeds

	 Can we optimize this function better. For example using restric pointers, or
	 removing branch points or any other suggestion
#endif
  
/**************************************************************/
