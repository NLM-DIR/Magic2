 /*  File: bigarraysub.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copuyright (C) J Thierry-Mieg yand R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description:
 *              Extension of Array, indexed by long integers
 * Exported functions:
 *              See Header file: bigarray.h (includes lots of macros)
 * HISTORY:
 * Last edited:
 * Created: Aug 24 2014 (mieg)
 *-------------------------------------------------------------------
 */

/* Warning : we have no provision to store a big array in the database
 * but we can memory-map them from and to disk with a simple user interface
 * bigArrayMapRead / bigArrayMapWrite
 */

#include "regular.h"
#include "bigarray.h"
#include "bitset.h"
#include <sys/mman.h>
#include <stdint.h>

  /* Defines bitField bit array from bitset.h */

extern BOOL finalCleanup ;	/* in messubs.c */


/************ Array : class to implement variable length arrays ************/

static long int bigTotalAllocatedMemory = 0 ;
static mysize_t bigTotalNumberCreated = 0 ;
static mysize_t bigTotalNumberActive = 0 ;
static Array reportBigArray = 0 ;
static void uBigArrayFinalise (void *cp) ;

static char * bigArrayAlloc (long int n, int size, int *shift)
{
  mysize_t nn = (mysize_t) n * size + 64 ;
  char *cp = malloc (nn) ;
  memset(cp,0,nn) ;

    uintptr_t a = (uintptr_t)cp;
    uintptr_t aligned = (a + 63) & ~63ULL;

    *shift = (int)(aligned - a);          // 0..63

    // cp + *shift   is always 64-byte aligned
    // and cp is located inside the allocated area

    return cp;
} /* bigArrayAlloc */

/**************/

#ifndef MEM_DEBUG
  BigArray uBigArrayCreate (long int n, int size, AC_HANDLE handle)
{ mysize_t id = bigTotalNumberCreated++ ;
  BigArray neuf = (BigArray) handleAlloc (uBigArrayFinalise, 
				   handle,
				   sizeof (struct BigArrayStruct)) ;
#else
BigArray   uBigArrayCreate_dbg (long int n, int size, AC_HANDLE handle,
					      const char *hfname,int hlineno) 
{ mysize_t id = bigTotalNumberCreated++ ;  
  BigArray neuf = (BigArray) handleAlloc_dbg (uBigArrayFinalise, 
				   handle,
				   sizeof (struct BigArrayStruct),
				   dbgPos(hfname, hlineno, __FILE__), __LINE__) ;
#endif

  int shift = 0 ;
  if (!reportBigArray)
    { reportBigArray = (Array)1 ; /* prevents looping */
      reportBigArray = arrayCreate (512, BigArray) ;
    }
  if (size <= 0)
    messcrash("negative size %d in uArrayCreate", size) ;
  if (n < 1)
    n = 1 ;
  if (reportBigArray != (Array)2)
    bigTotalAllocatedMemory += n * size ;
  neuf->trueBase = bigArrayAlloc (n , size, &shift) ;
  neuf->base = neuf->trueBase + shift ;
  neuf->dim = n ;
  neuf->max = 0 ;
  neuf->size = size ;
  neuf->id = ++id ;
  neuf->magic = BIG_ARRAY_MAGIC ;
  bigTotalNumberActive++ ;
  if (reportBigArray != (Array)1 &&
      reportBigArray != (Array)2) 
    { if (neuf->id < 20000)
	array (reportBigArray, neuf->id, BigArray) = neuf ;
      else
	{ Array aa = reportBigArray ;
	  reportBigArray = (Array)1 ; /* prevents looping */
	  arrayDestroy (aa) ;
	}
    }
  return neuf ;
}

/**************/

mysize_t bigArrayReportMark (void)
{
  return reportBigArray != (Array)1 && reportBigArray != (Array)2 ?  
      arrayMax(reportBigArray) : 0 ;
}

/**************/
/* call with j==-2 to block all static reports, needed in multithreaded code */
void bigArrayReport (mysize_t j)
{ mysize_t i ;
  BigArray a ;

  if (j == -2)
    {
      if (0) fprintf (stderr, "// bigArrayReport(-2) blocking all static reports in wego code\n") ;
      reportBigArray = (Array)2 ;
    }
  if (reportBigArray == (Array)2) 
    return ;
  
  fprintf(stderr,
	  "\n\n %lu active bigArrays, %lu created, %lu kb allocated\n\n ",   
	  (unsigned long) bigTotalNumberActive,   (unsigned long) bigTotalNumberCreated , (unsigned long) bigTotalAllocatedMemory/1024) ;
  
  if (reportBigArray == (Array)1) 
    return ;
  
  fprintf(stderr,"\n\n") ;
  
  i = arrayMax (reportBigArray) ;
  while (i-- && i > j)
    { a = arr (reportBigArray, i, BigArray) ;
      if (bigArrayExists(a))
	fprintf (stderr, "BigArray %lu  size=%u max=%lu dim=%lu\n",  (unsigned long) i, a->size,  (unsigned long) a->max,  (unsigned long) a->dim) ;
    }
}

/**************/

void bigArrayStatus (mysize_t *nmadep, mysize_t *nusedp,  long int *memAllocp, long int  *memUsedp)
{ 
  mysize_t i ;
  BigArray a, *ap ;

  *nmadep = bigTotalNumberCreated ; 
  *nusedp = bigTotalNumberActive ;
  *memAllocp = bigTotalAllocatedMemory ;
  *memUsedp = 0 ;

  if (reportBigArray == (Array)1) 
    return ;
  if (reportBigArray == (Array)2) 
    return ;

  i = arrayMax(reportBigArray) ;
  ap = arrp(reportBigArray, 0, BigArray) - 1 ;
  while (ap++, i--)
    if (bigArrayExists (*ap))
      { a = *ap ;
	*memUsedp += a->max * a->size ;
      }
}

/**************/

BigArray uBigArrayReCreate (BigArray a, long int n, int size)
{ 
  if (!bigArrayExists(a))
    return  uBigArrayCreate(n, size, 0) ;

  if(a->size != size)
    messcrash("Type  missmatch in uBigArrayRecreate, you should always "
	      "call recreate using the same type") ;

  if (n < 1)
    n = 1 ;
  if (a->dim < n || 
      (a->dim - n)*size > (1 << 19) ) /* free if save > 1/2 meg */
    {
      int shift = 0 ;
      if (reportBigArray != (Array)2)
	bigTotalAllocatedMemory -= a->dim * size ;
      free (a->trueBase) ;
      a->dim = n ;
      if (reportBigArray != (Array)2)
	bigTotalAllocatedMemory += n * size ;
      a->trueBase = bigArrayAlloc (n, size, &shift) ;
      a->base = a->trueBase + shift ;
    }
  else
    memset(a->base,0,(mysize_t)(a->dim*size)) ;

  a->max = 0 ;
  return a ;
} /* uBigArrayReCreate */

/**************/

void uBigArrayDestroy (BigArray a)
/* Note that the finalisation code attached to the memory does the work, 
   see below */
{
  if (!a) return;

  if (a->magic != BIG_ARRAY_MAGIC)
    messcrash ("uBigArrayDestroy received corrupt array->magic");

  if (a->lock)
     messcrash ("bigArrayDestroy called on locked bigArray") ;
  a->magic = 0 ;
  messfree(a);
}

static void uBigArrayFinalise (void *cp)
{
  BigArray a = (BigArray)cp;
  
  if (reportBigArray != (Array)2)
    bigTotalAllocatedMemory -= a->dim * a->size ;

  if (a->map)
    { /* memory mapped read only bigArray */
      if (a->readOnly && a->max != a->readOnly)
	messcrash ("A read only memory mapped bigArray has been modified : %s", a->fName) ;
      if (a->map != (void*)1 && munmap(a->map, a->max * a->size) == -1)
	messcrash ("Failed to unmap bigArray %s", a->fName) ;
      if (a->fd != -1)
	close (a->fd) ;
      a->base = a->trueBase = a->map = 0 ; a->fd = 0 ; a->readOnly = 0 ;
      messfree (a->fName) ;
    }
  else if (!finalCleanup)
    {
      free (a->trueBase) ;
      a->base = a->trueBase = 0 ;
    }
  a->magic = 0 ;
  bigTotalNumberActive-- ;
  if (!finalCleanup && reportBigArray != (Array)1 && reportBigArray != (Array)2) 
    arr(reportBigArray, a->id, BigArray) = 0 ;
}

/******************************/

void bigArrayLock (BigArray a)
{
  if (!bigArrayExists(a))
    messcrash ("bigArrayLock called on non exiting bigArray") ;
  if (a->lock)
    messcrash ("bigArrayLock called on already locked bigArray") ;
  a->lock = TRUE ;
  return ;
} /* bigArrayLock */

void bigArrayUnlockLock (BigArray a)
{
  if (!bigArrayExists(a))
    messcrash ("bigArrayUnlock called on non exiting bigArray") ;
  if (!a->lock)
    messcrash ("bigArrayLock called on unlocked bigArray") ;
  a->lock = FALSE ;
  return ;
} /* bigArrayLock */



/******************************/

#ifndef MEM_DEBUG
void bigArrayExtend (BigArray a, long int n)
#else
  void bigArrayExtend_dbg (BigArray a, long int  n, const char *hfname,int hlineno) 
#endif
{
  char *base, *trueBase ;
  int shift = 0 ;
  
  if (!a || n < a->dim)
    return ;
  if (a->lock)
     messcrash ("bigArrayExtend called on locked bigArray") ;
  if (n > 0x8000000000000000)
    messcrash ("BigArray in acedb cannot exceed max-long-int (2^63) bytes, sorry, n=%ld max=%ld cell size = %d\n"
	       , (long)n, a->max, a->size) ;

  if (reportBigArray != (Array)2)
    bigTotalAllocatedMemory -= a->dim * a->size ;
  if (a->dim*a->size < 1 << 30)  /* 256 megs of keys, or 1G of ram */
    a->dim *= 2 ;
  else
    a->dim += 1024 + ((1 << 30) / a->size) ;
  if (n >= a->dim)
    a->dim = n + 1 ;

  if (reportBigArray != (Array)2)
    bigTotalAllocatedMemory += a->dim * a->size ;
  trueBase = bigArrayAlloc (a->dim, a->size, &shift) ;
  base = trueBase + shift ;
  
  memcpy (base, a->base,a->size*a->max) ;
  free (a->trueBase) ;

  a->trueBase = trueBase ;
  a->base = base ;
}

/***************/

char *uBigArray (BigArray a, long int i)
{
  if (i < 0)
    messcrash ("referencing BigArray element %ld < 0", i) ;
  if (!a)
    messcrash ("uBigArray called with NULL Array struc");

  if (i >= a->max)
    { if (i >= a->dim)
        bigArrayExtend (a,i) ;
      a->max = i+1 ;
    }
  return a->base + i*a->size ;
}

/***************/

 char *uBigArrayCheck (BigArray a, long int i, int size)
{
  if (! a)
    messcrash ("dereferecning a null BigArray") ;
  if (size != a->size)
    messcrash ("BigArray size mismatch accessing array of size %d with pointer if typse-size %d", a->size, size) ;
  if (i < 0)
    messcrash ("referencing BigArray element %ld < 0", i) ;

  return uBigArray(a, i) ;
}

/***************/

 char *uBigArrCheck (BigArray a, long int i, int size)
{
  if (! a)
    messcrash ("dereferecning a null BigArray") ;
  if (size != a->size)
    messcrash ("BigArray size mismatch accessing array of size %d with pointer if typse-size %d", a->size, size) ;
  if (i >= a->max || i < 0)
    messcrash ("BigArray index %ld out of bounds [0,%ld]",
	       i, a->max - 1) ;
  return a->base + i*a->size ;
}

/**************/

#ifndef MEM_DEBUG
  BigArray bigArrayHandleCopy (BigArray a, AC_HANDLE handle)
#else
  BigArray bigArrayHandleCopy_dbg(BigArray a, const char *hfname,int hlineno, AC_HANDLE handle) 
#endif
{ BigArray b ;
  
  if (bigArrayExists (a) && a->size)
    {
#ifndef MEM_DEBUG
      b = uBigArrayCreate (a->max, a->size, handle) ;
#else
	  b = uBigArrayCreate_dbg (a->max, a->size, handle,
				dbgPos(hfname, hlineno, __FILE__), __LINE__) ;
#endif
      memcpy(b->base, a->base, a->max * a->size);
      b->max = a->max ;
      return b;
    }
  else
    return 0 ;
}
#ifndef MEM_DEBUG
BigArray bigArrayCopy (BigArray a)
{
  return bigArrayHandleCopy (a, 0) ;
}
#else
BigArray bigArrayCopy_dbg(BigArray a, const char *hfname,int hlineno) 
{
  return bigArrayHandleCopy_dbg (a, hfname, hlineno, 0) ;
}
#endif 
/**************/

BigArray bigArrayTruncatedCopy (BigArray a, long int x1, long int x2)
{ BigArray b = 0 ;
  
  if (x1 < 0 || x2 < x1 || x2 > a->max)
    messcrash 
      ("Bad coordinates x1 = %ld, x2 = %ld in arrayTruncatedCopy",
       x1, x2) ;
  if (bigArrayExists (a) && a->size)
    { if (x2 - x1)
	{ b = uBigArrayCreate (x2 - x1, a->size, 0) ;
	  b->max = x2 - x1 ;
	  memcpy(b->base, a->base + x1, b->max * b->size);
	}
      else
	b = uBigArrayCreate (10, a->size, 0) ;
    }
  return b;
}

/**************/ 

void bigArrayCompress(BigArray a)
{
  long int i, j , as ;
  char *x, *y, *ab  ;
  
  if (!a || !a->size || bigArrayMax(a) < 2 )
    return ;

  ab = a->base ; 
  as = a->size ;
  for (i = 1, j = 0 ; i < bigArrayMax(a) ; i++)
    { x = ab + i * as ; y = ab + j * as ;
      /*
	for (k = a->size ; k-- ;)		
	if (*x++ != *y++) 
	goto different ;
      */
      if (! memcmp (x, y, as))
	continue ;
      
      /*     different: */
      if (i != ++j)
	{ x = ab + i * as ; y = ab + j * as ;
	  /* 
	     for (k = a->size ; k-- ;)	 
	     *y++ = *x++ ;
	     */
	  memcpy (y , x, as) ;
	}
    }
  bigArrayMax(a) = j + 1 ;
}

/****************/

/* 2025_09_10 , added arraySortSlice */
void bigArraySortSlice (BigArray a, long int pos1, long int pos2, int (*order)(const void*, const void*))
{
  long int n = pos2 - pos1 ;
  
  if (! bigArrayExists (a))
    messcrash("bigArraySortSlice called on non existing bigArray, please edit the source code") ;
  if (pos1 < 0)
    messcrash("pos1 = %ld < 0 in bigArraySortSlice", pos1) ;
  if (pos2 > a->max)
    messcrash("pos2 = %ld > max = %ld in bigArraySortPosSlice", pos2, a->max) ;
  
  if (n > 1) 
    {
      int s = a->size ;
      void *v = a->base + s * pos1 ;
      mSort (v, n, s, order) ;
    }
}

/* 31.7.1995 dok408  added arraySortPos() - restricted sorting to tail of array */
void bigArraySortPos (BigArray a, long int  pos, int (*order)(const void*, const void*))
{ return bigArraySortSlice(a, pos, a ? a->max : 0, order) ;  }

void bigArraySort (BigArray a, int (*order)(const void*, const void*))
{ return bigArraySortSlice(a, 0, a ? a->max : 0, order) ;  }


/***********************************************************/

BOOL bigArrayIsEntry (BigArray a, long int i, void *s)
{
  char *cp = (char *) uBigArray(a,i), *cq = (char *)s ;
  int j = a->size;

  while (j--)
    if (*cp++ != *cq++) 
      return FALSE ;
  return TRUE;
}

/***********************************************/
       /* Finds Entry s from Array  a
        * sorted in ascending order of order()
        * If found, returns TRUE and sets *ip
        * if not, returns FALSE and sets *ip one step left
        */

BOOL bigArrayFind(BigArray a, void *s, long int  *ip, int (* order)(const void*, const void*))
{
  int ord;
  long int i = 0 , j = arrayMax(a), k;

  if(!j || (ord = order(s,uBigArray(a,0)))<0)
    { if (ip)
	*ip = -1; 
      return FALSE;
    }   /* not found */

  if (ord == 0)
    { if (ip)
	*ip = 0;
      return TRUE;
    }

  if ((ord = order(s,uBigArray(a,--j)))>0 )
    { if (ip)
	*ip = j; 
      return FALSE;
    }
  
  if (ord == 0)
    { if (ip)
	*ip = j;
      return TRUE;
    }

  while(TRUE)
    { k = i + ((j-i) >> 1) ; /* midpoint */
      if ((ord = order(s, uBigArray(a,k))) == 0)
	{ if (ip)
	    *ip = k; 
	  return TRUE;
	}
      if (ord > 0) 
	(i = k);
      else
	(j = k) ;
      if (i == (j-1) )
        break;
    }
  if (ip)
    *ip = i ;
  return FALSE;
} /* bigArrayFind(BigArray */

/**************************************************************/
/* Removes Entry s from Array  a
 * sorted in ascending order of order()
 */
BOOL bigArrayRemove (BigArray a, void * s, int (* order)(const void*, const void*))
{
  long int i;

  if (bigArrayFind(a, s, &i,order))
    {
      /* memcpy would be faster but regions overlap
       * and memcpy is said to fail with some compilers
       */
      char *cp = uBigArray(a,i),  *cq = cp + a->size ;
      mysize_t j = (bigArrayMax(a) - i)*(a->size) ;
      while(j--)
	*cp++ = *cq++;

      bigArrayMax(a) --;
      return TRUE;
    }
  else

    return FALSE;
}

/**************************************************************/
       /* Insert Segment s in Array  a
        * in ascending order of s.begin
        */

BOOL bigArrayInsert(BigArray a, void * s, int (*order)(const void*, const void*))
{
  long int i, j;

  if (bigArrayFind(a, s, &i,order))
    return FALSE;  /* no doubles */
  
  j = bigArrayMax(a) + 1;
  uBigArray(a,j-1) ; /* to create space */

	/* avoid memcpy for same reasons as above */
  {
    char* cp = uBigArray(a,j - 1) + a->size - 1,  *cq = cp - a->size ;
    long int k = (j - i - 1)*(a->size);
    while(k--)
      *cp-- = *cq--;
    
    cp = uBigArray(a,i+1); cq = (char *) s; k = a->size;
    while(k--)
      *cp++ = *cq++;
  }
  return TRUE;
}

/**********************************************************************/
/********* BitSet - inherits from BigArray **********/
/* mieg 2001, turn these into function to systematically make room
   allowing to use bit() without errors
   which is ok since anyway bitSetMax is not a defined method
*/
BitSet bitSetCreate (unsigned long int n, AC_HANDLE h)
{
  BigArray bb  = 0 ;
  if (n < 256) n = 256 ;
  bb = bigArrayHandleCreate (1 + (n >> 5), unsigned int, h) ;
  bigArray (bb, (n >> 5), unsigned int) = 0 ;
  if (sizeof(int) != 4)
    messcrash ("sorry, the biset code in w1/bigArraysubs.c assummes sizeof(int) = %d should be 4", sizeof (int)) ;
  return bb ;
}

BitSet bitSetReCreate (BitSet bb, unsigned long int n) 
{
  if (n == 0) 
    n = (bigArrayMax (bb)) >> 5 ;
  if (n < 256) n = 256 ;
  bb = bigArrayReCreate (bb, 1 + (n >> 5), unsigned int) ;
  bigArray (bb, (n >> 5), unsigned int) = 0 ;
  return bb ;
}

void bitExtend (BigArray bb, unsigned long int n)
{
  if (n < 256) n = 256 ;
  bigArrayExtend(bb, n >> 5) ;
  bigArray (bb, (n >> 5), unsigned int) = 0 ;
}

unsigned long int bitSetCount (BitSet a)
{
  register unsigned int cc, *cp ;
  register unsigned long int i1 ;
  register unsigned long int j = bigArrayExists(a) ? bigArrayMax(a) : 0 , n = 0 ;

  if(!j)
    return 0 ;
  cp  = bigArrp(a,0,unsigned int) ;
  while (j--)
    {
      if(!*cp) n+=32 ;
      else if(cc = ~(*cp), cc != 0)
	for (i1 = 0 ; i1 < 32 ; i1++, cc >>= 1) if (cc & 0x1) n++ ;
      cp++;
    }
  return 32*bigArrayMax(a) - n ;
}

/* performs b1 = b1 & b2, returns count (b1) */
unsigned long int bitSetAND (BitSet b1, BitSet b2)
{
  register unsigned int cc, *cp, *cq ;
  register unsigned long int i1, n = 0 ;
  register unsigned long int i = bigArrayExists(b1) ? bigArrayMax(b1) : 0 ;
  register unsigned long int j = bigArrayExists(b2) ? bigArrayMax(b2) : 0 ;

  if (i > j) i = j ;
  bigArrayMax (b1) = i ;
  if (i == 0)
    return 0 ;

  cp  = bigArrp(b1,0,unsigned int) ;
  cq  = bigArrp(b2,0,unsigned int) ;
  while (i--)
    { 
      *cp &= *cq ;
      if(!*cp) n+=32 ;
      else if(cc = ~(*cp), cc != 0)
	for (i1 = 0 ; i1 < 32 ; i1++, cc >>= 1) if (cc & 0x1) n++ ;
      cp++; cq++ ;
    }
  return 32*bigArrayMax(b1) - n ;
} /* bitSetAND */

/* count b1 & b2, but without modifying b1 or b2 */
unsigned long int bitSetANDcount (BitSet b1, BitSet b2)
{
  register unsigned int cc, *cp, *cq ;
  register unsigned long int n = 0 ;
  register unsigned long int i = bigArrayExists(b1) ? bigArrayMax(b1) : 0 ;
  register unsigned long int j = bigArrayExists(b2) ? bigArrayMax(b2) : 0 ;
  
  if (i > j) i = j ;
  if (i == 0)
    return 0 ;
  
  cp  = bigArrp(b1,0,unsigned int) ;
  cq  = bigArrp(b2,0,unsigned int) ;
  while (i--)
    { 
      cc = *cp & *cq ;
      switch ((int) cc)
	{
	case 0: n += 32 ;  break ;
	case 0xffffffff: break ;
	default:
	  j = 32 ; cc = -cc ;
	  while (j--)
	    if (cc & 0x1) n++ ;
	}
      cp++; cq++ ;
    }
  return 32*bigArrayMax(b1) - n ;
} /* bitSetANDcount */

/* performs b1 = b1 | b2, returns count (b1) */
unsigned long int bitSetOR (BitSet b1, BitSet b2)
{
  register unsigned int cc, *cp, *cq ;
  register unsigned long int i1, n = 0 ;
  register unsigned long int i = bigArrayExists(b1) ? bigArrayMax(b1) : 0 ;
  register unsigned long int j = bigArrayExists(b2) ? bigArrayMax(b2) : 0 ;

  if (i < j)
    {
      i = j ;
      bigArray(b1,i - 1,unsigned int) = 0 ;
    }
  if (i > j)
    bigArray(b2, i - 1,unsigned int) = 0 ;


  if (i == 0)
    return 0 ;

  cp  = bigArrp(b1,0,unsigned int) ;
  cq  = bigArrp(b2,0,unsigned int) ;
  while (i--)
    { 
      *cp |= *cq ; /* *cp OR *cq */
      if(!*cp) n+=32 ;
      else if(cc = ~(*cp), cc != 0)
	for (i1 = 0 ; i1 < 32 ; i1++, cc >>= 1) if (cc & 0x1) n++ ;
      cp++; cq++ ;
    }
  return 32*bigArrayMax(b1) - n ;
} /* bitSetOR */

/* performs b1 = b1 ^ b2, returns count (b1) */
unsigned long int bitSetXOR (BitSet b1, BitSet b2)
{
  register unsigned int cc, *cp, *cq ;
  register unsigned long int i1, n = 0 ;
  register unsigned long int i = bigArrayExists(b1) ? bigArrayMax(b1) : 0 ;
  register unsigned long int j = bigArrayExists(b2) ? bigArrayMax(b2) : 0 ;

  if (i < j)
    {
      i = j ;
      bigArray(b1,i - 1,unsigned int) = 0 ;
    }
  if (i > j)
    bigArray(b2, i - 1,unsigned int) = 0 ;


  if (i == 0)
    return 0 ;

  cp  = bigArrp(b1,0,unsigned int) ;
  cq  = bigArrp(b2,0,unsigned int) ;
  while (i--)
    { 
      *cp ^= *cq ; /* *cp XOR *cq */
      if(!*cp) n+=32 ;
      else if(cc = ~(*cp), cc != 0)
	for (i1 = 0 ; i1 < 32 ; i1++, cc >>= 1) if (cc & 0x1) n++ ;
      cp++; cq++ ;
    }
  return 32*bigArrayMax(b1) - n ;
} /* bitSetXOR */

/* performs b1 = b1 \ b2, returns count (b1) */
unsigned long int bitSetMINUS (BitSet b1, BitSet b2)
{
  register unsigned int cc, *cp, *cq ;
  register unsigned long int i1, n = 0 ;
  register unsigned long int i = bigArrayExists(b1) ? bigArrayMax(b1) : 0 ;
  register unsigned long int j = bigArrayExists(b2) ? bigArrayMax(b2) : 0 ;

  if (i < j)
    {
      i = j ;
      bigArray(b1,i - 1,unsigned int) = 0 ;
    }
  if (i > j)
    bigArray(b2, i - 1,unsigned int) = 0 ;


  if (i == 0)
    return 0 ;

  cp  = bigArrp(b1,0,unsigned int) ;
  cq  = bigArrp(b2,0,unsigned int) ;
  while (i--)
    { 
      *cp = *cp ^ (*cp & *cq) ; /* *cp XOR ( *cp AND *cq) */
      if(!*cp) n+=32 ;
      else if(cc = ~(*cp), cc != 0)
	for (i1 = 0 ; i1 < 32 ; i1++, cc >>= 1) if (cc & 0x1) n++ ;
      cp++; cq++ ;
    }
  return 32*bigArrayMax(b1) - n ;
} /* bitSetMINUS */

/*************************************************************/
/* perfmeters */
long int assBounce = 0, assBounceValue = 0, assBucketMax = 0, assFound = 0, assNotFound = 0, assInserted = 0, assRemoved = 0 ;

/* Associator package is for associating pairs of pointers. 
   Assumes that an "in" item is non-zero and non -1.
   Implemented as a hash table of size 2^m.  
   
   Originally grabbed from Steve Om's sather code by Richard Durbin.

   Entirelly rewritten by mieg, using bouncing by relative primes
   and deletion flagging.

   User has access to structure member ->n = # of pairs
*/

#define VSIZE (sizeof(void*))
#define VS5 (VSIZE/5)
#define VS7 (VSIZE/7)
#define moins_un ((void*) (-1))

#define HASH_old(_xin) { register unsigned long z = VS5, x = (char*)(_xin) - (char*)0 ; \
		  for (hash = x, x >>= 5 ; z-- ; x >>= 5) hash ^= x ; \
		  hash &= a->mask ; \
		}

#define DELTA_old(_xin)   { register unsigned long z = VS7, x = (char*)(_xin) - (char*)0 ; \
		       for (delta = x, x >>= 7 ; z-- ; x >>= 7) delta ^= x ; \
		       delta = (delta & a->mask) | 0x01 ; \
		   }  /* delta odd is prime relative to  2^m */

/*
 Multiplication Method (Cormen). Choose m to be a power of 2. Let A be some random-looking real number. Knuth suggests M = 0.5*(sqrt(5) - 1). Then do the following:

     s = k*A
     x = fractional part of s
     h(k) = floor(m*x)

This seems to be the method that the theoreticians like.

To do this quickly with integer arithmetic, let w be the number of bits in a word (e.g. 32) and suppose m is 2^p. Then compute:

     s = floor(A * 2^w)
     x = k*s
     h(k) = x >> (w-p)      // i.e. right shift x by (w-p) bits
                            // i.e. extract the p most significant 
                            // bits from x
*/
static long unsigned int  _A_HASH = 0 ;
static long unsigned int  _D_HASH = 0 ;

#define HASH(_xin) { hash  =   (((long unsigned int)(_xin - 0x0) * _A_HASH) >> ((VSIZE << 3) - a->nbits)) & a->mask ;}
#define DELTA(_xin){ delta =  ((((long unsigned int)(_xin - 0x0) * _D_HASH) >> ((VSIZE << 3) - a->nbits)) & a->mask) | 0x1 ;}

#ifdef JUNK
from Yann
/****************************************************************************/
2	/*                                                                          */
3	/* This file is part of libDDD, a library for manipulation of DDD and SDD.  */
4	/*                                                                          */
5	/*     Copyright (C) 2001-2008 Yann Thierry-Mieg, Jean-Michel Couvreur      */
6	/*                             and Denis Poitrenaud                         */
7	/*    Based on a file written by Alexandre Duret-Lutz for Spot,             */
8	/*                                     Alexandre.Duret-Lutz@lip6.fr         */
9	/*                                                                          */
10	/*     This program is free software; you can redistribute it and/or modify */
11	/*     it under the terms of the GNU Lesser General Public License as       */
12	/*     published by the Free Software Foundation; either version 3 of the   */
13	/*     License, or (at your option) any later version.                      */
14	/*     This program is distributed in the hope that it will be useful,      */
15	/*     but WITHOUT ANY WARRANTY; without even the implied warranty of       */
16	/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
17	/*     GNU LEsserGeneral Public License for more details.                   */
18	/*                                                                          */
19	/* You should have received a copy of the GNU Lesser General Public License */
20	/*     along with this program; if not, write to the Free Software          */
21	/*Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
22	/*                                                                          */
23	/****************************************************************************/
24	#ifndef __DDD_MISC_HASHFUNC_HH
25	#define __DDD_MISC_HASHFUNC_HH
26	
27	/******************************************************************************/
28	
29	#include <stdint.h>
30
31	namespace ddd
32	{
33	  /// \addtogroup hash_funcs Hashing functions
34	  /// \ingroup misc_tools
35	  /// @{
36	
37	  /// \brief Thomas Wang's 32 bit hash function.
38	  ///
39	  /// Hash an integer amongst the integers.
40	  /// http://www.concentric.net/~Ttwang/tech/inthash.htm
41	  inline size_t
42	  wang32_hash(size_t key)
43	  {
44	    // We assume that size_t has at least 32bits.
45	    key += ~(key << 15);
46	    key ^=  (key >> 10);
47	    key +=  (key << 3);
48	    key ^=  (key >> 6);
49	    key += ~(key << 11);
50	    key ^=  (key >> 16);
51	    return key;
52	  }
53	
54	  /// Another of Wang's fast hash with a magic number.
55	  /// good for (sequence of) integers
56	  inline uint32_t int32_hash(uint32_t a) {
57	    a = (a ^ 61) ^ (a >> 16);
58	    a = a + (a << 3);
59	    a = a ^ (a >> 4);
60	    a = a * 0x27d4eb2d;
61	    a = a ^ (a >> 15);
62	    return a;
63	  }
64	
65	
66	  /// \brief Knuth's Multiplicative hash function.
67	  ///
68	  /// This function is suitable for hashing values whose
69	  /// high order bits do not vary much (ex. addresses of
70	  /// memory objects).  Prefer spot::wang32_hash() otherwise.
71	  /// http://www.concentric.net/~Ttwang/tech/addrhash.htm
72	  inline size_t
73	  knuth32_hash(size_t key)
74	  {
75	    // 2654435761 is the golden ratio of 2^32.  The right shift of 3
76	    // bits assumes that all objects are aligned on a 8 byte boundary.
77	    return (key >> 3) * 2654435761U;
78	  }
79	  /// @}
80	}
81

#endif

#define roll(_h) ((_h << 7) | (_h >> (32 - 7)))
#define roll7(_h) ((_h << 7) | (_h >> ((VSIZE<<3) - 7)))
#define roll9(_h) ((_h << 9) | (_h >> ((VSIZE<<3) - 9)))
#define BARKER 0x1111100110101


#ifdef OLEG_GOOD_IDEA_does_not_quite_work_as_written
/* barker code have very low autocorrelation properties 
 * see wikipedia

 ** 2010_08_28
 ** hasing the whole genome works faster with the original acedb HASH function
 */
/* OLEG_RAND should be a random numer, i wrote in hexa some decimals of Pi */
#define OLEG_RAND 0xc47fd0d2ef287c75
/* static unsigned long oleg_start = OLEG_RAND ; */
static unsigned long T1[256], T2[256] ;
static void assHashOlegInit (void)
{
  /* const int *TAB = isDiff ? (const int *)dictHashOleg : (const int *)stackTokeniseTextOn ;  a random set ? */
  static int T= 0 ;

  if (T==0)
    {
      unsigned long k = (unsigned long) BARKER ;
      for (T = 0 ; T < 256 ; T++)
	{ T1[T] = k ^ 0x532814ea ; k = roll(k) ; T2[T]= k ^ 0x47e2132 ; } 
    }
  return ;
}
#define DELTA3(_xin)   { const unsigned long *TAB = T1  ; register int i=sizeof(void*)/sizeof(char) ; \
              register unsigned long h=BARKER, x = (char*)(_xin) - (char*)0 ; \
		       while(i--) {h = roll9(h) + (TAB[(unsigned char)(h ^ (x & 0xff))]^(x & 0xff)) ; x>>=8;} \
                       h ^= h >> 16 ; h ^= (h >> 8) ; delta = (h & a->mask) | 1 ; \
		   }  /* delta odd is prime relative to  2^m */
#define HASH3(_xin)   { const unsigned long *TAB = T2 ; register int i=sizeof(void*)/sizeof(char) ; \
              register unsigned long h=BARKER, x = (char*)(_xin) - (char*)0 ; \
		       while(i--) {h = roll7(h) + (TAB[(unsigned char)(h ^ (x & 0xff))]^(x & 0xff)) ; x>>=8;} \
                       h ^= h >> 16 ; h ^= (h >> 8) ; hash = (h & a->mask) ; \
		   }  /* delta odd is prime relative to  2^m */


static unsigned int assHashOleg (void *s, unsigned mask, BOOL isDiff)
{ /* contributed by Oleg Khovayko, 2007 */
  unsigned long h = (unsigned int)BARKER ;
  unsigned char *cp = (unsigned char *)s ;
  unsigned int i, h1 ;
  unsigned long *TAB = isDiff ? T1 : T2 ;

  for (i = 0 ; i < 8 ; cp++, i++)
    h = roll(h) + (TAB[(unsigned char)(h ^ *cp)] ^ *cp) ;
  h ^= h >> 32 ; /* accumulate the bits on the right */
  h ^= h >> 16 ; /* accumulate the bits on the right */
  h ^= (h >> 8) ;
  /* h &= (1 << n) - 1 ;      compress down to n bits */
  h1 = ((unsigned int)h) & mask ;
  if (isDiff)			/* return odd number */
    h1 |= 1 ;

  return h1 ;
}  /* assHashOleg */

#define HASH2(_xin)  { hash=assHashOleg (_xin,a->mask,0);}
#define DELTA2(_xin) { delta=assHashOleg (_xin,a->mask,1);}

#endif

/**************** Destruction ****************************/

static void assFinalise(void *cp)
{ Associator a = (Associator)cp;

  a->magic = 0 ;
  if (!finalCleanup)
    { 
      messfree(a->in) ;
      messfree(a->out);
      messfree (a->bbh) ;
    }
}

void uAssDestroy (Associator a)
{ if (assExists(a))
    messfree (a) ;  /* calls assFinalise */
}

/**************** Creation *******************************/
/* before we had one bitset for hash, one for delta with protection 3
 * following oleg, i now have a single bitset wit protection 4
 */
#define PROTECTION 4

#ifndef MEM_DEBUG
  static Associator assDoCreate (int nbits, BOOL inOnly, AC_HANDLE handle)
#else
    static Associator assDoCreate_dbg (int nbits, BOOL inOnly, AC_HANDLE handle,
				 const char *hfname, int hlineno)
#endif
{
  static int nAss = 0 ;
  Associator a ;
  unsigned long int size, vsize, size2 ;

  if ( _A_HASH == 0)
    {
      long unsigned int un = 1 ;
      double z = un << (VSIZE << 2) ;
      _A_HASH = ((double)0.5 * (sqrt(5.0) - 1)) * z * z ;
      _D_HASH = ((double)0.5 * (sqrt(13.0) - 1)) * z * z ;
    }


  size = 1 << nbits ;  /* size must be a power of 2 */
  size2 = size << PROTECTION ;  /* size must be a power of 2 */
  vsize = size * VSIZE ; 

  /*   assHashOlegInit () ; */
#ifndef MEM_DEBUG
  a = (Associator) handleAlloc(assFinalise, 
			       handle, 
			       sizeof (struct AssStruct)) ;
  a->in = (const void**) messalloc (vsize) ;
  if (! inOnly)
    {
      a->out = (const void**) messalloc (vsize) ;
    }
#else
  a = (Associator) handleAlloc_dbg(assFinalise, 
				   handle, 
				   sizeof (struct AssStruct),
				   dbgPos(hfname, hlineno, __FILE__), __LINE__) ;
  a->in = (const void**) messalloc_dbg(vsize,dbgPos(hfname, hlineno, __FILE__), __LINE__) ;
  if (! inOnly)
    {
      a->out = (const void**) messalloc_dbg(vsize,dbgPos(hfname, hlineno, __FILE__), __LINE__) ;
    }
#endif
  a->magic = ASS_MAGIC ;
  a->id = ++nAss ;
  a->readOnly = FALSE ;
  a->n = 0 ;
  /* a->i = 0 ; was start up position for recursive calls */
  a->nbits = nbits ;
  a->mask = size2 - 1 ;
  if (nbits > 8)
    {
      a->bbh = bitSetCreate (size2, 0) ;
    }
  a->inOnly = inOnly ;
  a->floatingLine = 2 ; /* do not fill in/out over 1/2, and the bit table over 1/(1<<floatingLine), i.e. the default is 1/16 */
  /* this allows to create a very empty bit-array and recognize very fast non existing
     words, as for example when we scan the genome to look for a hit to a small number of words
  */

  return a ;
}

/*******************/

#ifndef MEM_DEBUG
Associator assBigHandleCreate (long int size, AC_HANDLE h)
#else
  Associator assBigHandleCreate_dbg(long int size, AC_HANDLE h, const char *hfname, int hlineno)
#endif
{
  int n = 2 ; /* make room, be twice as big as needed */

  if (size <= 0) 
    messcrash ("assBigCreate called with size = %d <= 0", size) ;

  --size ;
  while (size >>= 1) n++ ; /* number of left most bit + 1 */
#ifndef MEM_DEBUG
  return assDoCreate(n, FALSE, h) ;
#else
  return assDoCreate_dbg(n, FALSE, h, hfname, hlineno) ;
#endif
}

#ifndef MEM_DEBUG
Associator assInOnlyHandleCreate (long int size, AC_HANDLE h)
#else
  Associator assInOnlyHandleCreate_dbg(long int size, AC_HANDLE h, const char *hfname, int hlineno)
#endif
{
  int n = 2 ; /* make room, be twice as big as needed */

  if (size <= 0) 
    messcrash ("assBigCreate called with size = %d <= 0", size) ;

  --size ;
  while (size >>= 1) n++ ; /* number of left most bit + 1 */
#ifndef MEM_DEBUG
  return assDoCreate(n, TRUE, h) ;
#else
  return assDoCreate_dbg(n, TRUE, h, hfname, hlineno) ;
#endif
}

/*******************/

#ifndef MEM_DEBUG
  Associator assHandleCreate(AC_HANDLE handle) { return assDoCreate(5, FALSE, handle) ;}
#else
  Associator assHandleCreate_dbg(AC_HANDLE handle, const char *hfname,int hlineno)
  { return assDoCreate_dbg(5, FALSE, handle, dbgPos(hfname, hlineno, __FILE__), __LINE__) ; }
#endif

/*******************/

void assClear (Associator a)
{
  long int size, size2, vsize ;
  if (!assExists(a)) 
    return ;

  if (a->readOnly)
    messcrash ("assClear called on a readOnly associator") ;
  
  a->n = 0 ;
  size2 = a->mask + 1 ;
  size = size2 >> PROTECTION ;
  vsize = size * VSIZE ;
  
  memset(a->in, 0, vsize) ;
  if (a->inOnly)
    {
      messfree (a->out) ;
    }
  else
    {    
      memset(a->out, 0, vsize) ;
    }
  if (a->bbh)
    {
      a->bbh = bitSetReCreate (a->bbh, size2) ;
    }
}

/********************/

Associator assReCreate (Associator a)
{ if (!assExists(a))
    return assHandleCreate(0) ;
  else
    {
      if (a->readOnly)
	messcrash ("assReCreate called on a readOnly associator") ;
  
      assClear (a) ;  
      return a ;
    }
}

/********************/
/* forbids assInsert allowing multithreading */
BOOL assSetReadOnly (Associator a, BOOL b) 
{
  if (!assExists(a))
    messcrash ("SetReadOnly received corrupted associator") ;
  a->readOnly = b ;
  return b ;
}

static void assDouble (Associator a)
{ 
  unsigned long int oldsize = 1, newsize = 1, newsize2, newvsize ;
  int nbits ;
  register unsigned long hash, delta, mask2 ;
  const void **old_in = a->in, **old_out = a->out, **xin, **xout ;
  long unsigned int i ;

  nbits = a->nbits ;
  oldsize <<= nbits ;
  nbits++ ;
  newsize <<= nbits ; /* size must be a power of 2 */
  newsize2 = newsize << PROTECTION ;  /* size must be a power of 2 */
  newvsize = newsize * VSIZE ; 

  a->n = 0 ;
  /* a->i = 0 ; was start up position for recursive calls */
  a->nbits = nbits ;
  a->mask = newsize2 - 1 ;

  mask2 = (a->mask >> PROTECTION) ; 
  if (nbits > 8) /* a one bit bitset for presence of the main hash, no bouncing */
    {
      if (a->bbh)
	{
	  a->bbh = bitSetReCreate (a->bbh, newsize2) ;
	}
      else
	{
	  a->bbh = bitSetCreate (newsize2, 0) ;
	}
    }
  
  a->in  = (const void**) messalloc (newvsize) ;
  if (! a->inOnly)
    a->out = (const void**) messalloc (newvsize) ;
  
  for (i = 0, xin = old_in, xout = old_out ; i < oldsize ; i++, xin++, xout++)
    if (*xin && *xin != moins_un)
      { 
	HASH(*xin) ; 
	DELTA(*xin) ;
	if (a->bbh) /* a one bit bitset for presence of the main hash, no bouncing */
	  {
	    /* register unsigned long x = (char*)(*xin) - (char*)0 ;  */
	    bitSet (a->bbh, hash) ;
	    bitSet (a->bbh, delta) ;
	  }
	hash &= mask2 ; delta &= mask2 ;
	while (TRUE)
	  { 
	    if (!a->in[hash])  /* don't need to test moins_un */
	      { 
		a->in[hash] = *xin ;
		if (! a->inOnly)
		  a->out[hash] = *xout ; /* works as is for buckets, since we do not look at the content */
		++a->n ;
		break ;
	      }
	    hash = (hash + delta) & mask2 ;
	  }
      }

  messfree (old_in) ;
  messfree (old_out) ; 
} /* assDouble */

/************************ Searches  ************************************/

static BOOL assDoFind (Associator a, const void* xin, const void** pout, unsigned long int *hashPosp)
/* if found, updates *pout and returns TRUE, else returns FALSE	*/
{ 
  unsigned long hash, delta = 0, mask2 ;
  const void* test ;
  mask2 = (a->mask >> PROTECTION) ; 

  if (!assExists(a))
    messcrash ("assFind received corrupted associator") ;
  if (!xin || xin == moins_un) return FALSE ;

  HASH(xin) ;
  if (a->bbh && ! bit (a->bbh, hash))
    return FALSE ; 
  DELTA(xin) ;
  if (a->bbh && ! bit (a->bbh, delta))
    return FALSE ; 
      
  if (a->inOnly)
    return TRUE ;

  hash &= mask2 ; delta &= mask2 ;
  while (TRUE)
    { test = a->in[hash] ;
      if (test == xin)
	{ if (pout)
	    *pout = a->out[hash] ;
	  assFound++ ;
	  /* a->i = hash ; */
	  if (hashPosp) *hashPosp = hash + 1 ; /* + 1 so we do not pass 0, reserved for initialisation */
	  return TRUE ;
	}
      if (!test)
	break ;
      assBounce++ ;
      hash = (hash + delta) & mask2 ;
    }
  assNotFound++ ;
  return FALSE ;
} /* assDoFind */

/********************/
  /* Usage: Array bucket = 0 ; int iBucket = 0  ; while (assFindNext(..., &bucket, &iBucket)) ;
   * This new method is thread safe, since the place holder is held   * by the client rather than by the array structure
   */

BOOL assFindNext (Associator ass, const void* xin, const void** pout, Array *bucketp, int *iBucketp)
/* if found, updates *pout and returns TRUE, else returns FALSE	*/
{ 
  long unsigned int v ;
  int k = *iBucketp ;
  Array aa = *bucketp ;

  if (k == 0) aa = *bucketp = 0 ;
   if (aa && !arrayExists(aa))
    messcrash ("assFindNext received corrupted array") ;
   if (!assExists(ass))
    messcrash ("assFindNext received corrupted associator") ;
  if (ass->inOnly && pout)
    messcrash ("assFindNext pout != NULL called on inOnly associator") ;
  if (!xin)
    return FALSE ;
  if (xin == moins_un)
    return FALSE ;
     
  if (pout) *pout = 0 ;
  if (! aa && k == 1)  /* previous was not a bucket */
    return FALSE ; 
  if (! aa && k)  /* previous was not a bucket */
    messcrash ("assFindNext called with null *bucketp AND non null iBucket = %d", k) ;
  if (aa && k <= 0)  /* previous was a bucket but the index is wrong  */
    messcrash ("assFindNext called with *iBucketp = 0 AND non null *bucket") ;

  if (aa && arrayMax(aa) && k)
    {
      if (k > 0 && k < arrayMax(aa))
	{
	  *iBucketp = k+1 ;
	  if (pout) *pout = arr (aa, k, void*) ;
	  return TRUE ;
	}
      *iBucketp = 0 ;
      *bucketp = 0 ;
      return FALSE ;
    }
  
  *bucketp = 0 ;
  *iBucketp = 0 ; 
  if (! assDoFind (ass, xin, pout, 0))
    return FALSE ;
  if (!pout)
    { 
      *iBucketp = 1 ; 
      return TRUE ; 
    }
  v = *pout - NULL ;
  if (((v >> 62) & 0x3) == 0x1) /* we are in a bucket */
    {
       long int n = v & 0x3fffffffffffffff ;
       *pout = 0 ;
       if (ass->aa && n < bigArrayMax (ass->aa))
	 {
	   Array a = bigArr (ass->aa, n, Array) ;
	   k = a ? arrayMax (a) : 0 ;
	   if (k)
	     {
	       *bucketp = a ;
	       *iBucketp = 1 ;
	       *pout = arr (a, 0, const void*) ;
	       return TRUE ;
	     }
	   return FALSE ;
	 }
     }
   else
     {
       *bucketp = 0 ;
       *iBucketp = 1 ; 
       return TRUE ;
     }
  return FALSE ;
} /* assFindNext */

/************************/

BOOL uAssFind (Associator ass, const void* xin, const void** pout, unsigned long int *hashPosp)
{
  Array bucket = 0 ;
  int iBucket = 0 ;
  return assFindNext (ass, xin, pout, &bucket, &iBucket) ;
} /* uAssFind */

/********************/

void  bucketAssTest (void)
{
  Associator a = assHandleCreate(0) ;
  Array bucket = 0 ;
  int x , y, iBucket ; 
  const void *vp, *wp ;
  
  x = 17 ; y = 1 ;
  y = 1 ; assMultipleInsert (a, assVoid(x), assVoid (y)) ;
  y = 2 ; assMultipleInsert (a, assVoid(x), assVoid (y)) ;
  y = 3 ; assMultipleInsert (a, assVoid(x), assVoid (y)) ;
  x = 12 ;
  y = 11 ; assMultipleInsert (a, assVoid(x), assVoid (y)) ;
  y = 12 ; assMultipleInsert (a, assVoid(x), assVoid (y)) ;
  y = 13 ; assMultipleInsert (a, assVoid(x), assVoid (y)) ;
  x = 20 ;
  y = 20 ; assMultipleInsert (a, assVoid(x), assVoid (y)) ;
  x = 21 ;
  y = 21 ; assMultipleInsert (a, assVoid(x), assVoid (y)) ;


  printf ("Expected result : 17 11 to 17 13, then 12 1 to 12 3 then 20 20, 21 21  : 8 lines total \n") ;
  
  for (x = 10 ; x < 22 ; x++)
    {
      iBucket = 0 ;  bucket = 0 ;
      while (assFindNext (a, assVoid (x), &vp, &bucket, &iBucket)) 
	printf ("x = %d, y = %d\n", x, assInt (vp)) ;
    }

  vp = wp = 0 ;
  
  printf ("\n\nsecond export iterating on assNext; expoect same reult in a different order\n") ;
  if (1) while (assNext (a, &wp, &vp))
    printf ("x = %d, y = %d\n",  assInt (wp), assInt (vp)) ;
  exit (0) ;
}

/************************ Insertions  ************************************/

     /* if already there returns FALSE, else inserts and returns TRUE */
/* XQX
*
* Inserts a (key,value) pair into an associator.  Returns true if
* the (key,value) was inserted, or false if the insert failed because
* the key is already there.
*
* if noMultiples is false, the insert cannot fail.
*/

static BOOL assDoInsert (Associator a, const void* xin, const void* xout, BOOL noMultiples)
{ 
  unsigned long hash, delta = 0, mask2 ;
  const void* test ;

  if (!assExists(a))
    messcrash ("assInsert received corrupted associator") ;

  if (!xin || xin == moins_un) 
    messcrash ("assInsert received forbidden value xin == 0") ;

  if (a->readOnly)
    messcrash ("assInsert called on a readOnly associator") ;

  while (a->n+1  >= (a->mask >> (PROTECTION + a->floatingLine))) /* reaching floating line */
    assDouble (a) ;

  mask2 = (a->mask >> PROTECTION) ; 
  HASH (xin) ; DELTA (xin) ;

  if (a->bbh) /* a one bit bitset for presence of the main hash, no bouncing */
    {
      bitSet (a->bbh, hash) ;  
      bitSet (a->bbh, delta) ;  
    }
  hash &= mask2 ; delta &= mask2 ;
  while (TRUE)
    { 
      test = a->in[hash] ;
      if (!test || test == moins_un)  /* reuse deleted slots */
	{ 
	  a->in[hash] = xin ;
	  if (! a->inOnly) a->out[hash] = xout ;
	  ++a->n ;
	  assInserted++ ;
	  return TRUE ;
	}
      /*  if (noMultiples && test == xin)		// already there 
	return FALSE ;
      */
      if (test == xin)		/* already there */
	{
	  if (noMultiples ||  a->inOnly)
	    return FALSE ;
	  if (a->out[hash] == xout)
	    return FALSE ;
	  else
	    assBounceValue++ ;
	}
      else
	assBounce++ ;
     hash = (hash + delta) & mask2 ;
    }
} /* assDoInsert */

/*****************/
     /* This one does not allow multiple entries in one key. */
BOOL assInsert (Associator a, const void* xin, const void* xout)
{ 
  return assDoInsert (a, xin, xout, TRUE) ;
} /* assInsert */
 
/*****************/

BOOL assMultipleInsert (Associator ass, const void* xin, const void* xout)
     /* This one allows multiple entries in one key. */
{
  if (ass->inOnly)
    messcrash ("bucketAssMultipleInsert called on inOnly associator") ;
  unsigned long hash, delta = 0, mask2 ;
  const void* test ;

  if (!assExists(ass))
    messcrash ("assInsert received corrupted associator") ;

  if (!xin || xin == moins_un) 
    messcrash ("assInsert received forbidden value xin == 0") ;

  if (ass->readOnly)
    messcrash ("assInsert called on a readOnly associator") ;

  while (((ass->n+1) << (PROTECTION + ass->floatingLine)) >= ass->mask) /* reaching floating line */
    assDouble (ass) ;

  mask2 = (ass->mask >> PROTECTION) ; 
  { Associator a = ass ; HASH (xin) ; DELTA (xin) ; } /* the macro refers to a, not very kosher, sorry */

  if (ass->bbh) /* a one bit bitset for presence of the main hash, no bouncing */
    {
      bitSet (ass->bbh, hash) ;  
      bitSet (ass->bbh, delta) ;  
    }
  hash &= mask2 ; delta &= mask2 ;
  while (TRUE)
    { 
      test = ass->in[hash] ;
      if (!test || test == moins_un)  /* reuse deleted slots */
	{ 
	  long unsigned int v = ass->out[hash] - NULL ;
	  long int k ;
	  ass->in[hash] = xin ; 
	  if (((v >> 62) & 0x3) == 0x1) /* the entry mimics a bucket */
	    {
	      /* create a bucket */
	      Array a ;
	      if (! ass->aa)
		ass->aa = bigArrayHandleCreate (256, Array, ass->h) ;
	      k = bigArrayMax (ass->aa) ;
	      a = bigArray (ass->aa, k, Array) = arrayHandleCreate (8, const void*, ass->h) ;
	      array (a, 0, const void *) = xout ; 
	      ass->out[hash] = (void*)0 + ((((long unsigned int)0x1) << 62) | k) ;	  
	    }
	  else
	    ass->out[hash] = xout ;
	  ++ass->n ;
	  assInserted++ ;
	  return TRUE ;
	}
      if (test == xin)		/* already there */
	{
	  long unsigned int v = ass->out[hash] - NULL ;
	  if (((v >> 62) & 0x3) == 0x1) /* we are in a bucket */
	    {
	      long int n = v & 0x3fffffffffffffff ;
	      Array a = bigArr (ass->aa, n, Array) ;
	      int i, k = a ? arrayMax (a) : 0 ;
	      if (k)
		{
		  for (i = 0 ; i < k && i < 3 ; i++)
		    if (arr (a, i, void*) == xout)
		      return FALSE ;
		  array (a, k, const void*) = xout ;
		  if (k > assBucketMax)  assBucketMax = k ;
		  return TRUE ;
		}
	      else
		messcrash ("bucketAssMultipleInsert  k == 0, internal error") ;
	    } 
	  if (ass->out[hash] == xout)
	    return FALSE ;
	  /* create a bucket */
	  if (! ass->aa)
	    ass->aa = bigArrayHandleCreate (256, Array, ass->h) ;
	  {
	    long int n = bigArrayMax (ass->aa) ;
	    Array a = bigArray (ass->aa, n, Array) = arrayHandleCreate (8, const void*, ass->h) ;
	    array (a, 0, const void *) = ass->out[hash] ;
	    array (a, 1, const void *) = xout ; 
	    if (assBucketMax < 2) assBucketMax = 2 ;
	    ass->out[hash] = (void*)0 + ((((long unsigned int)0x1) << 62) | n) ;	 
	  } 
	  return TRUE ;
	}
      else
	assBounce++ ;
     hash = (hash + delta) & mask2 ;
    }
  return FALSE ;
} /* bucketAssMultipleInsert */

/************************ Removals ************************************/
   /* if found, removes entry and returns TRUE, else returns FALSE	*/
   /* No a->n--, because the entry is blanked but not actually removed */
BOOL assRemove (Associator ass, const void* xin)
{ 
  unsigned long int hashPos = 0 ;
  if (ass->readOnly)
    messcrash ("assRemove called on a readOnly associator") ;

  if (assExists(ass) && assDoFind (ass, xin, 0, &hashPos))
    { 
      long unsigned int v = ass->out[hashPos-1] - NULL ;
      if (((v >> 62) & 0x3) == 0x1) /* we are in a bucket */
	{
	  long int n = v & 0x3fffffffffffffff ;
	  Array a = (ass->aa && n < bigArrayMax (ass->aa) ? bigArr (ass->aa, n, Array) : 0) ;
	  if (a) 
	    {
	      arrayDestroy (a) ;
	      bigArr (ass->aa, n, Array) = 0 ;
	    }
	}
      ass->in[hashPos-1] = moins_un ; ass->out[hashPos-1] = 0 ;
      assRemoved++ ;
      return TRUE ;
    }
  else
    return FALSE ;
} /* assRemove */

/*******************/

/* if found, removes entry and returns TRUE, else returns FALSE	*/
/* Requires both xin and xout to match */

/************************ dumpers ********************************/
     /* lets you step through all members of the table */
BOOL uAssNext (Associator ass, const void* *pin, const void* *pout)
{ 
  long int size = 1 ;
  const void *test ;
  long unsigned int ii ;
  long unsigned int v ;

  if (!assExists(ass))
     messcrash("uAssNext received a non existing associator") ;
  if (ass->inOnly)
    messcrash ("assNext called on inOnly associator") ;
  size <<= ass->nbits ;

  ii = ass->i ;
  if (!*pin)
    { ass->i = ii = 0 ; ass->iBucket = -1 ; }
  else if (*pin != ass->in[ii])
    { messerror ("Non-consecutive call to assNext()") ;
      return FALSE ;
    }
  ass->iBucket++ ;

 lao:
  test = ass->in[ii] ;  
  if (test && test != moins_un) /* not empty or deleted */
    {
      v = ass->out[ii] - NULL ; 
      if (pout && ((v >> 62) & 0x3) == 0x1) /* we are in a bucket */
	{
	  long int n = v & 0x3fffffffffffffff ;
	  Array a = (ass->aa && n < bigArrayMax (ass->aa) ? bigArr (ass->aa, n, Array) : 0) ;
	  if (a && arrayMax (a) >  ass->iBucket)
	    {
	      *pin = ass->in[ii] ;
	      *pout = arr (a, ass->iBucket, void*) ;
	      return TRUE ;
	    }
	}
      else if (ass->iBucket == 0)
	{
	  *pin = ass->in[ii] ;
	  if (pout) *pout = ass->out[ii] ;
	  return TRUE ;
	}
    }
  /* else we must loop till we find another entry */
   
  while (++ii < size)
    { test = ass->in[ii] ;
      if (test && test != moins_un) /* not empty or deleted */
	{ 
	   ass->i = ii ;
	   *pin = ass->in[ii] ;
	   ass->iBucket = 0 ;
           goto lao ;
	}
    }

  return FALSE ;
} /* uAssNext */

/************************ dumpers ********************************/
     /* lets you step through all members of the table */
BOOL uNewAssNext (Associator a, const void* *pin, const void* *pout, unsigned long int *hashPosp)
{ 
  long int size = 1 ;
  const void *test ;
  long int ii ;

  if (!assExists(a))
     messcrash("uAssNext received a non existing associator") ;
  if (a->inOnly)
    messcrash ("assNext called on inOnly associator") ;
  size <<= a->nbits ;

  ii = *hashPosp - 1 ;
  if (!*pin)
    ii = -1 ;
  else if (*pin != a->in[ii])
    { messerror ("Non-consecutive call to assNext()") ;
      return FALSE ;
    }

  while (++ii < size)
    { test = a->in[ii] ;
      if (test && test != moins_un) /* not empty or deleted */
	{ *pin = a->in[ii] ;
	  if (pout)
	    *pout = a->out[ii] ;
	  *hashPosp  = ii + 1 ;
	  return TRUE ;
	}
    } 

  return FALSE ;
} /* uNewAssNext */

/*******************/

void assDump (Associator a)
{ 
  long int i, x ; 
  const void **in, **out ;
  
  if (!assExists(a)) 
    {
      fprintf(stderr,"assDump fails assExists\n");
      return ;
    }
  
  if (a->inOnly)
    messcrash ("assDump called on inOnly associator") ;
  i = a->mask;
  in = a->in; out = a->out;
  /* keep stderr here since it is for debugging */
  fprintf (stderr,"Associator %lx : %ld pairs %ld slots\n",(unsigned long)a,(unsigned long)a->n,(unsigned long)i) ;
  for (x=0; x<=i; x++)
    if (in[x])
      fprintf(stderr,"\t%16p - %16p\n", in[x],  out[x]) ;
} /*assDump */

/**********************************************************************/
/**********************************************************************/
/* File MEMORY MAPPING
 *  This functionalilty is covered inside the BigArray interface
 *  allowing easy array access and automatic handling of deallocation
 *   Notice that in ReadOnly case the array is ReadOnly
 *   It can only be freed via ac)free (a) or ac_free(h)
 *   but it should not be modified
 *   Since array aree accessed directly via pointers
 *   It is the responsability of the client to do this properly
 *   Checks are made if the array is extended, but not generally
 */
/* bigArray can be memory mapped to and from file */

/*
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>
*/
#include <sys/mman.h>

#ifdef JUNK
/* possibilities to study */
#include <sys/statvfs.h>
 struct statvfs stat;
 statvfs(".", &stat);
 uint64_t free_space = stat.f_bavail * stat.f_frsize;
 if (free_space < total_size) {
    fprintf(stderr, "Insufficient disk space\n");
    return -1;
 }
 write_buffers("temp.bin", buf, sizes, N);
 rename("temp.bin", "buffers.bin");char *chunk;
 posix_memalign((void **)&chunk, 4096, chunk_size);
#endif

/********************************************************/
/* write a buffer to a memory mapped file */
static BOOL mmapWrite (const char *fName, char *buf, long int size)
{
  /*  Write a buffer using mmap */
  /* Open the file */
  int fd = open(fName, O_RDWR | O_CREAT | O_TRUNC, 0666);
  if (fd == -1)
    messcrash ("Failed to create a memory map file of size %ld : %s", size, fName) ;

  /* Extend the file to the requested size */
  if (ftruncate(fd, size) == -1)
    messcrash ("Failed to extend a memory map file of size %ld : %s", size, fName) ;

  /* map the file into memory */
  char *map = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  if (map == MAP_FAILED)
    messcrash ("Failed to map a memory map file of size %ld : %s", size, fName) ;

  /* copy the buffer to the mapped region */
  memcpy (map, buf, size);

  /* ensure the data are  written to disk */
  if (msync (map, size, MS_SYNC) == -1) 
    messcrash ("Failed to sync a memory map file of size %ld : %s", size, fName) ;

  /* iunmap and close */
  if (munmap (map, size) == -1)
    messcrash ("Failed to unmap a memory map file of size %ld : %s", size, fName) ;
  
  close(fd);
  return TRUE ;
} /* mmapWrite */

/********************************************************/

static unsigned char *mmapCreate (const char *fName, long int *size, int *fdp, unsigned char **mapp, BOOL readOnly)
{
  unsigned char  *buf = 0 ;
  /* Open the file */
  int fd = open (fName, O_RDONLY);
  if (fd == -1)
    {
      messerror ("Failed to open for reading a memory map file : %s", fName) ;
      return NULL ;
    }

  /*  Get file size */
  struct stat st;
  if (fstat(fd, &st) == -1)
    {
      messerror ("Failed to stat for reading a memory map file of size %ld : %s", size, fName) ;
      close(fd);
      return NULL ;
    }
  *size = st.st_size;

  /* Map the file to memory */
  /* the last arg, here 0, says that mapping should start n pages down whwere pageSize = sysconf(_SC_PAGE_SIZE). */
  /* PROT can be PROT_NONE or   'PROT_READ | PROT_WRITE | PROT_EXE == 777' */
  /* MAP_PRIVATE | MAP_NORESERVE ; do not write do not reserve swap as we do not write */
  /* MAP_PRIVATE | MAP_NORESERVE ; do not write do not reserve swap as we do not write */
#ifdef DARWIN
  /* MAP_POPULATE and MAP_NORESERVE are not available on MacOS */
  unsigned char *map = mmap (NULL, *size, PROT_READ, MAP_PRIVATE, fd, 0);
#else
  unsigned char *map = mmap (NULL, *size, PROT_READ, MAP_PRIVATE | MAP_NORESERVE | MAP_POPULATE , fd, 0);
#endif

  if (map == MAP_FAILED)
    {
      messerror ("Failed to map for reading a memory map file of size %ld : %s", size, fName) ;
      close(fd);
      return NULL ;
    }  /* the last arg, hemre 0, says that mapping should start n pages down whwere pageSize = sysconf(_SC_PAGE_SIZE). */
  
  /* Real size */
  if (0) memcpy (size, map, 8) ;

  /*  Allocate output buffer and copy
   * use map directly if read-only
   */
  if (readOnly)
    buf = map ; /* this is not allocated */
  else
    {
      buf = (unsigned char *) messalloc (*size);
      memcpy (buf, map, *size);
    }

#ifdef JUNK
  /* Verify CRC32, requires linking with -lz  (zlib compression library)  */
  long unsigned int crc, expected_crc = crc32(0, buf, size);
  memcpy(&crc, map + *size, 8) ;
  if (crc != expected_crc)
    messcrash ("Failed to read the CRC memory map file of size %ld : %s", size, fName) ;
#endif
  if (! readOnly)
    { /*  Unmap and close */
      if (munmap(map, *size) == -1)
	messcrash ("Failed to unmap memory map file of size %ld : %s", size, fName) ;
      close(fd);
      map = 0 ; fd = -1 ; /* prevent double deallocation and closure */
    }
  *mapp = map ;
  *fdp = fd ;
  return  buf ;
} /* mmapCreate */

/********************************************************/

BigArray uBigArrayMapRead (const char *fName, int recordSize, BOOL readOnly, AC_HANDLE h)
{
  unsigned char *map = 0 ;
  BigArray aa = uBigArrayCreate (8, recordSize, h) ;
  long int size = 0 ;
  int fd = -1 ;
  if (! fName)
    messcrash ("bigArrayMapRead called with NULL file name") ;
  free (aa->trueBase) ; aa->base = aa->trueBase = 0 ;
  aa->base = aa->trueBase = (char *) mmapCreate (fName, &size, &fd, &map, readOnly) ;
  aa->readOnly = readOnly ;
  aa->size = recordSize ;
  aa->max = aa->dim = size/recordSize ;
  if (readOnly)
    aa->readOnly = aa->max ;
  aa->map = map ? map : (void *)1 ;
  aa->fd = fd ;
  aa->fName = strnew (fName, 0) ;
  return aa ;
} /* uBigArrayMapRead */

/**********************************************************************/

BOOL bigArrayMapWrite (BigArray aa, const char *fName)
{
  if (aa->readOnly)
    messcrash ("cannot MapWrite a readOnly bigArray") ;
  return mmapWrite (fName, aa->base, aa->max * aa->size) ;
} /* bigArrayMapWrite  */

/**********************************************************************/
/**********************************************************************/
/**** BigStack : arbitrary BigStack class - inherits from BigArray ****/
/**********************************************************************/

static void uBigStackFinalise (void *cp) ;

#ifndef MEM_DEBUG
BigStack bigStackHandleCreate (long int n, AC_HANDLE handle)  /* n is initial size */
{
  BigStack s = (BigStack) handleAlloc (uBigStackFinalise, 
				       handle,
				       sizeof (struct BigStackStruct)) ;

#else
BigStack bigStackHandleCreate_dbg (long int n, AC_HANDLE handle,  /* n is initial size */
				   const char *hfname, int hlineno)
{
  BigStack s = (BigStack) handleAlloc_dbg (uBigStackFinalise, 
					   handle,
					   sizeof (struct BigStackStruct),
					   hfname, hlineno) ;
#endif
  s->magic = BIGSTACK_MAGIC ;
  s->a = bigArrayCreate (n,char) ;
  s->pos = s->ptr = s->a->base ;
  s->safe = s->a->base + s->a->dim - 16 ; /* safe to store objs up to size 8 */
  return s ;
}

BigStack bigStackReCreate (BigStack s, long int n)               /* n is initial size */
{
  if (!bigStackExists(s))
    messcrash ("bigStackReCeate called on NULL bigSatck, use bigStackhandleCreate") ;

  s->a = bigArrayReCreate (s->a,n,char) ;
  s->pos = s->ptr = s->a->base ;
  s->safe = s->a->base + s->a->dim - 16 ; /* safe to store objs up to size 8 */
  return s ;
}

BigStack bigStackCopy (BigStack old, AC_HANDLE handle)
{
  BigStack neuf = 0 ;

  if (bigStackExists(old))
    {
      neuf = bigStackHandleCreate (old->a->dim, handle) ;
      memcpy (neuf->a->base, old->a->base, old->a->dim) ;
      neuf->ptr = neuf->a->base + (old->ptr - old->a->base) ;
      neuf->pos = neuf->a->base + (old->pos - old->a->base) ;
    }
  return neuf ;
}

void uBigStackDestroy(BigStack s)
{ if (s && s->magic == BIGSTACK_MAGIC) messfree(s);
} /* the rest is done below as a consequence */

static void uBigStackFinalise (void *cp)
{ BigStack s = (BigStack)cp;
  if (!finalCleanup) bigArrayDestroy (s->a) ;
  s->magic = 0 ;
}

void bigStackExtend (BigStack s, long int n)
{
  long int ptr = s->ptr - s->a->base ;
  long int pos = s->pos - s->a->base ;
  s->a->max = s->a->dim ;	/* since only up to ->max copied over */
  bigArrayExtend (s->a,ptr+n+16) ;	/* relies on arrayExtend mechanism */
  s->ptr = s->a->base + ptr ;
  s->pos = s->a->base + pos ;
  s->safe = s->a->base + s->a->dim - 16 ;
}

long int bigStackMark (BigStack s)
{ return (s->ptr - s->a->base) ;
}

long int bigStackPos (BigStack s)
{ return (s->pos - s->a->base) ;
}

void bigStackCursor (BigStack s, long int mark)
{ s->pos = s->a->base + mark ;
}

void bigStackClear(BigStack s)
{
  if (bigStackExists(s))
    {
      s->pos = s->ptr = s->a->base;
      s->a->max = 0 ;
    }
}

/**********************************************************************/
/************************  end of file ********************************/
/**********************************************************************/
 
 
 
 
