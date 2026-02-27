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
#endif /* __SSE2__ */

#ifdef USEGPU
#include <time.h>
#include "sa.gpusort.h"
#endif /* USEGPU */

/**************************************************************/

static inline int cwOrder (const void *va, const void *vb)
{
  const CW *up = va ;
  const CW *vp = vb ;
  int n ;
  n = (up->seed > vp->seed) - (up->seed < vp->seed) ; if (n) return n ;
  n = up->nam - vp->nam ; if (n) return n ;
  n = (up->pos > vp->pos) - (up->pos < vp->pos) ; if (n) return n ;
  return 0 ;
} /* cwOrder */

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
  int n, n1, n2 ;
  
  n = ((up->read >> 1) > (vp->read >> 1)) -  ((up->read >> 1) < (vp->read >> 1)) ; if (n) return n ; 
  n = ((up->chrom > vp->chrom) - (up->chrom < vp->chrom)) ; if (n) return n ; 
  n1 = up->a1 + (up->x1 >> NSHIFTEDTARGETREPEATBITS) ;
  n2 = vp->a1 + (vp->x1 >> NSHIFTEDTARGETREPEATBITS) ;
  n = n1 - n2 ; if (n) return n ;
  n = ((up->x1 > vp->x1) - (up->x1 < vp->x1)) ; 

  return n ;
} /* hitPairOrder */

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
} /* insertionSort */

/* #ifndef USEGPU */
/* recursivelly split the table
 * the partially sorted data oscillate between b and buf
 * they end up correctly in b because for small n
 * we switch to the insertion taquin algoright
 * on correct parity, as speed is 2 persent higher with
 * insertion n>0,  relative to n==0
 * but n=8,16,32 are equivalent speeds
 */

static BOOL saSortDo (char *b, long int nn, int s, char *buf, BOOL hitIsTarget, int (*cmp)(const void *va, const void *vb))
{
 char *up, *vp, *wp ;
  long int n1 = nn / 2 ;
  long int n2 = nn - n1 ;
  char *b1 = b ;
  char *b2 = b + n1 * s ;
  char *b01 = buf ;
  char *b02 = buf + n1 * s ;
  int n = 0 ;
  BOOL ok = FALSE ;
  BOOL clean1, clean2, clean = TRUE ;
  /* for small n,
   * sort en place using the insertion algorithm (game of taquin)
   */
  if (hitIsTarget && nn < 8)
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
  /* code generated by Grok, loads and stores 16bytes (128 bits) */
  if (cmp == cwOrder)
    {
      while (n1 > 0 && n2 > 0)
	{
	  __m128i u = _mm_load_si128((__m128i*)up) ;
	  __m128i v = _mm_load_si128((__m128i*)vp) ;

	  int n = (*(unsigned int*)up <= *(unsigned int*)vp) ;
	  
	  _mm_store_si128((__m128i*)wp, n  ? u : v) ;
	  wp += s ;
	  up += n * s ;
	  vp += (1 - n) * s ;
	  n1 -= n ;
	  n2 -= 1 - n ;
	}
      ok = TRUE ;
    }
  else if (s == 16)
    {
      while (n1 > 0 && n2 > 0)
	{
	  __m128i u = _mm_load_si128((__m128i*)up) ;
	  __m128i v = _mm_load_si128((__m128i*)vp) ;
	  
	  int n = ((*cmp)(up, vp) <= 0) ? 1 : 0 ;
	  
	  _mm_store_si128((__m128i*)wp, n  ? u : v) ;
	  wp += s ;
	  up += n * s ;
	  vp += (1 - n) * s ;
	  n1 -= n ;
	  n2 -= 1 - n ;
	}
      ok = TRUE ;
    }
#endif
  if (! ok) /* either we do not have _mm_store_si128, or size s is not 16 */
    { /* classic code */
      while (n1 > 0 && n2 > 0)
	{
	  n = ((*cmp) (up, vp) <= 0) ? 1 : 0 ;

	  memcpy (wp, (n<=0) ? up : vp, s) ;
	  wp += s ;
	  up += n * s ;
	  vp += (1 - n) * s ;
	  n1 -= n ;
	  n2 -= 1 - n ;
	}
    }

  /* I also tried to count all greater cases and bulk copy
   * but this code was more complex and not faster
   */

  /* bulk copy the reminders */
    if (n1 > 0) memcpy(wp, up, n1 * s);
    if (n2 > 0) memcpy(wp, vp, n2 * s);

    return clean ;
} /* saSortDo */
/* #endif // USEGPU */

/**************************************************************/

int saSort (BigArray aa, int type)
{
  long int N = bigArrayMax (aa) ;
  char *cp = N ?  (char *) aa->base : 0 ;
  int s = aa->size ;
  int (*cmp)(const void *va, const void *vb) ;
  int usedGPU = 0 ;
  
  switch (type)
    {
    case 1:
      cmp = cwOrder ;
      break ;
    case 2:
      cmp = hitReadOrder ;
      break ;
    case 3:
      cmp = hitPairOrder ;
      break ;
    default:
      messcrash ("Wrong call to saSort typw = %dd>4", type) ;
    }
  if (N <= 1) return 0 ;
  else if (N < 128)
    newInsertionSort (cp, N, s, cmp) ;
  else
    {

        struct timespec start, end;
        double ellapsed;
        timespec_get(&start, TIME_UTC);


#ifdef USEGPU
	if (type < 3 && N > (1<<12))
	  {
	    usedGPU = 1 ;
	    saGPUSort (cp, N, type) ;
	  }
	else
	  {
	    char *buf = malloc (N * s) ;
	    memcpy (buf, cp, N * s) ;
	    saSortDo (cp, N, s, buf, TRUE, cmp) ;
	    free (buf) ;
	  }
#else
      char *buf = malloc (N * s) ;
      memcpy (buf, cp, N * s) ;
      saSortDo (cp, N, s, buf, TRUE, cmp) ;
      free (buf) ;
#endif

      timespec_get(&end, TIME_UTC);
      ellapsed = (double)(end.tv_sec - start.tv_sec) +
          (double)(end.tv_nsec - start.tv_nsec) / 1000000000.0;
      if (0) fprintf(stderr, "Sorted %ld elements in %f seconds (type: %d)\n", N, ellapsed, type);
    }
  return usedGPU ;
}/* saSsort */

/**************************************************************/
/**************************************************************/
/**************************************************************/
