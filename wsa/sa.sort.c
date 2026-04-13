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
static int saRadixSortSeedMatch (BigArray aa) ;
int saSort (BigArray aa, int type)
{
  long int N = bigArrayMax (aa) ;
  char *cp = N ?  (char *) aa->base : 0 ;
  int s = aa->size ;
  int (*cmp)(const void *va, const void *vb) = NULL ;
  int usedGPU = 0 ;
  
  if (N < 2) return FALSE ;

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
    newInsertionSort (cp, N, s, cmp) ;
  else
    {
      
#ifdef TIME_EVAL
      struct timespec start, end;
      double ellapsed;
      timespec_get(&start, TIME_UTC);
#endif
      
#ifdef USEGPU
      if (type < 3 && (size_t)N * s > (1 << 20))
	{
	  usedGPU = 1 ;
	   saGPUSort (cp, N, type) ;
	}
      else /* GPU threshold not met: fall through to CPU sort below */
#endif
	{
	  if (type == 5)  /* 20260412 to be tested for speed against case 4 */
	    saRadixSortSeedMatch (aa) ;
	  else
	    {
	      size_t alloc_size = ((size_t)N * s + 15) & ~15 ;
	      char *buf = aligned_alloc(16, alloc_size) ;
	      if (! buf) messcrash ("\nsa.sort.c alloc failure, consider lowering --bMax\n") ;
	      memcpy (buf, cp, N * s) ;
	      saSortDo (cp, N, s, buf, TRUE, cmp) ;
	      free (buf) ;
	    }
	}
      
      
#ifdef TIME_EVAL      
      timespec_get(&end, TIME_UTC);
      ellapsed = (double)(end.tv_sec - start.tv_sec) +
	(double)(end.tv_nsec - start.tv_nsec) / 1000000000.0;
      fprintf(stderr, "Sorted %ld elements in %f seconds (type: %d)\n", N, ellapsed, type);
#endif      
    }
  return usedGPU ;
}/* saSort */

/**************************************************************/
/**************************************************************/
/**************************************************************/
/* Two-pass LSB radix sort on SEEDMATCH.read (unsigned int key at offset 0)
 * vp  : pointer to SEEDMATCH array, must be 64-byte aligned
 * N   : number of records
 * The sort is stable and out-of-place internally (allocates one scratch buffer).
 * Returns 0 on success, -1 on allocation failure.
 
 
 * Two-pass LSB radix sort, 16 bits per pass, key = first unsigned int (read) 
 * Pass 1: sort on bits  0..15 (low  half of read)
 * Pass 2: sort on bits 16..31 (high half of read)
 */

#define RADIX_BITS  16
#define RADIX_SIZE  (1 << RADIX_BITS)   /* 65536 buckets */
#define RADIX_MASK  (RADIX_SIZE - 1)

static int saRadixSortSeedMatch (BigArray aa)
{
  long int N = aa ? bigArrayMax (aa) : 0 ;
  SEEDMATCH *src = bigArrp (aa, 0, SEEDMATCH) ;
  SEEDMATCH *dst ;
  long int   counts[2][RADIX_SIZE] ;
  long int   offsets[RADIX_SIZE] ;
  long int   i, bucket ;
  unsigned int key ;
  
  if (N <= 1)
    return 0 ;
  
  /* allocate scratch buffer, same alignment guarantee as input */
  dst = (SEEDMATCH *) aligned_alloc (64, (size_t) N * sizeof (SEEDMATCH)) ;
  if (!dst)
    return -1 ;
  
  /* --- histogram pass: count both 16-bit halves in one sequential scan --- */
  memset (counts, 0, sizeof (counts)) ;
  for (i = 0 ; i < N ; i++)
    {
      key = src[i].read ;
      counts[0][ key        & RADIX_MASK]++ ;   /* low  16 bits */
      counts[1][(key >> 16) & RADIX_MASK]++ ;   /* high 16 bits */
    }
  
  /* --- pass 1: sort on low 16 bits, src -> dst --- */
  offsets[0] = 0 ;
  for (i = 1 ; i < RADIX_SIZE ; i++)
    offsets[i] = offsets[i-1] + counts[0][i-1] ;
  
  for (i = 0 ; i < N ; i++)
    {
      bucket = src[i].read & RADIX_MASK ;
      dst[offsets[bucket]++] = src[i] ;
    }
  
  /* --- pass 2: sort on high 16 bits, dst -> src --- */
  offsets[0] = 0 ;
  for (i = 1 ; i < RADIX_SIZE ; i++)
    offsets[i] = offsets[i-1] + counts[1][i-1] ;
  
  for (i = 0 ; i < N ; i++)
    {
      bucket = (dst[i].read >> 16) & RADIX_MASK ;
      src[offsets[bucket]++] = dst[i] ;
    }
  
  /* result is back in src (original vp), dst is scratch only */
  free (dst) ;
  return 0 ;
}

/**************************************************************/
/**************************************************************/
/**************************************************************/
/**************************************************************/
/**************************************************************/

#ifdef JUNK
# AI prompt to analyze the code

This code module is part of the sortalign package (RNA aligner, NCBI/NLM/NIH).
  Key context:
  - C code targeting portability, SSE2 available, AVX2 avoided
  - USEGPU and TIME_EVAL are compile-time flags
  - BigArray, HIT, CW, SEEDMATCH, BOOL, messcrash() are defined in sa.h
  - aligned_alloc(16, ...) is used for SSE2-aligned buffers
  - mysize_t is the array size type
  
  This code runs on very large dataset (Tera bases of dna) and will be distributed
  and recompiled in many places
  
  The pupose of this chat is to analyze this module against any code weakness
  and suggesting any speed optimization
#endif
  
/**************************************************************/
