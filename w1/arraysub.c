/*  File: arraysub.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg yand R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description:
 *              Arbitrary length arrays, stacks, associators
 *              line breaking and all sorts of other goodies
 *              These functions are declared in array.h
 *               (part of regular.h - the header for libfree.a)
 * Exported functions:
 *              See Header file: array.h (includes lots of macros)
 * HISTORY:
 * Last edited: Dec  4 11:12 1998 (fw)
 * * Nov  1 16:11 1996 (srk)
 *		-	MEM_DEBUG code clean-up 
 *                      (some loose ends crept in from WIN32)
 *		-	int (*order)(const void*, const void*) prototypes used 
 *                                            uniformly throughout
 * * Jun  5 00:48 1996 (rd)
 * * May  2 22:33 1995 (mieg): killed the arrayReport at 20000
      otherwise it swamps 50 Mbytes of RAM
 * * Jan 21 16:25 1992 (mieg): messcrash in uArrayCreate(size<=0)
 * * Dec 17 11:40 1991 (mieg): stackTokeniseTextOn() tokeniser.
 * * Dec 12 15:45 1991 (mieg): Stack magic and stackNextText
 * Created: Thu Dec 12 15:43:25 1989 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: arraysub.c,v 1.74 2020/05/30 15:03:10 mieg Exp $	 */

   /* Warning : if you modify Array or Stack structures or 
      procedures in this file or array.h, you may need to modify 
      accordingly the persistent array package asubs.c.
   */

#include "regular.h"

/*#include <limits.h>*/

extern BOOL finalCleanup ;	/* in messubs.c */

/******** tells how much system stack used *********/

mysize_t stackused (void)
{ return 0 ;
}

/************ Array : class to implement variable length arrays ************/

static mysize_t totalAllocatedMemory = 0 ;
static mysize_t totalNumberCreated = 0 ;
static mysize_t totalNumberActive = 0 ;
static Array reportArray = 0 ;
static void uArrayFinalise (void *cp) ;

#ifndef MEM_DEBUG
  Array uArrayCreate (mysize_t n, int size, AC_HANDLE handle)
{ mysize_t id = totalNumberCreated++ ;
  Array neuf = (Array) handleAlloc (uArrayFinalise, 
				   handle,
				   sizeof (struct ArrayStruct)) ;
#else
Array   uArrayCreate_dbg (mysize_t n, int size, AC_HANDLE handle,
					      const char *hfname,int hlineno) 
{ mysize_t id = totalNumberCreated++ ;  
  Array neuf = (Array) handleAlloc_dbg (uArrayFinalise, 
				   handle,
				   sizeof (struct ArrayStruct),
				   dbgPos(hfname, hlineno, __FILE__), __LINE__) ;
#endif

  if (!reportArray)
    { reportArray = (Array)1 ; /* prevents looping */
      reportArray = arrayCreate (512, Array) ;
    }
  if (size <= 0)
    messcrash("negative size %d in uArrayCreate", size) ;
  if (n < 1)
    n = 1 ;
  if (reportArray != (Array)2)
    totalAllocatedMemory += n * size ;
#ifndef MEM_DEBUG
  neuf->base = (char *) messalloc (n*size) ;
#else
  neuf->base = (char *) messalloc_dbg (n*size,dbgPos(hfname, hlineno, __FILE__), __LINE__) ;
#endif
  neuf->dim = n ;
  neuf->max = 0 ;
  neuf->size = size ;
  neuf->id = ++id ;
  neuf->magic = ARRAY_MAGIC ;
  totalNumberActive++ ;
  if (reportArray != (Array)1 &&
      reportArray != (Array)2) 
    { if (neuf->id < 20000)
	array (reportArray, neuf->id, Array) = neuf ;
      else
	{ Array aa = reportArray ;
	  reportArray = (Array)1 ; /* prevents looping */
	  arrayDestroy (aa) ;
	}
    }
  return neuf ;
}

/**************/

mysize_t arrayReportMark (void)
{
  return reportArray != (Array)1 && reportArray != (Array)2 ?  
      arrayMax(reportArray) : 0 ;
}

/**************/
/* call with j==-2 to block all static reports, needed in multithreaded code */
void arrayReport (mysize_t j)
{ mysize_t i ;
  Array a ;

  if (j == -2)
    {
      if (0) fprintf (stderr, "// arrayReport(-2) blocking all static reports in wego code\n") ;
      reportArray = (Array)2 ;
    }
  if (reportArray == (Array)2) 
    return ;
  
  fprintf(stderr,
	  "\n\n %lu active arrays, %lu created, %lu kb allocated\n\n ",   
	  (unsigned long) totalNumberActive,   (unsigned long) totalNumberCreated , (unsigned long) totalAllocatedMemory/1024) ;
  
  if (reportArray == (Array)1) 
    return ;
  
  fprintf(stderr,"\n\n") ;
  
  i = arrayMax (reportArray) ;
  while (i-- && i > j)
    { a = arr (reportArray, i, Array) ;
      if (arrayExists(a))
	fprintf (stderr, "Array %lu  size=%u max=%lu dim=%lu\n",  (unsigned long) i, a->size,  (unsigned long) a->max,  (unsigned long) a->dim) ;
    }
}

/**************/

void arrayStatus (mysize_t *nmadep, mysize_t *nusedp, mysize_t *memAllocp, mysize_t *memUsedp)
{ 
  mysize_t i ;
  Array a, *ap ;

  *nmadep = totalNumberCreated ; 
  *nusedp = totalNumberActive ;
  *memAllocp = totalAllocatedMemory ;
  *memUsedp = 0 ;

  if (reportArray == (Array)1) 
    return ;
  if (reportArray == (Array)2) 
    return ;

  i = arrayMax(reportArray) ;
  ap = arrp(reportArray, 0, Array) - 1 ;
  while (ap++, i--)
    if (arrayExists (*ap))
      { a = *ap ;
	*memUsedp += a->max * a->size ;
      }
}

/**************/

Array uArrayReCreate (Array a, int n, int size)
{ if (!arrayExists(a))
    return  uArrayCreate(n, size, 0) ;

  if(a->size != size)
    messcrash("Type  missmatch in uArrayRecreate, you should always "
	      "call recreate using the same type") ;

  if (n < 1)
    n = 1 ;
  if (a->dim < n || 
      (a->dim - n)*size > (1 << 19) ) /* free if save > 1/2 meg */
    { 
      if (reportArray != (Array)2)
	totalAllocatedMemory -= a->dim * size ;
      messfree (a->base) ;
      a->dim = n ;
      if (reportArray != (Array)2)
	totalAllocatedMemory += a->dim * size ;
      a->base = (char *) messalloc (a->dim*size) ;
    }
  memset(a->base,0,(mysize_t)(a->dim*size)) ;

  a->max = 0 ;
  return a ;
}

/**************/

void uArrayDestroy (Array a)
/* Note that the finalisation code attached to the memory does the work, 
   see below */
{
  if (!a) return;

  if (a->magic != ARRAY_MAGIC)
    messcrash ("uArrayDestroy received corrupt array->magic");
  if (a->lock)
     messcrash ("arrayDestroy called on locked array") ;
  a->magic = 0 ;
  messfree(a);
}

static void uArrayFinalise (void *cp)
{
  Array a = (Array)cp;
  
  if (reportArray != (Array)2)
    totalAllocatedMemory -= a->dim * a->size ;
  if (!a->lock && !finalCleanup) messfree (a->base) ;
  a->magic = 0 ; a->base = 0 ;
  totalNumberActive-- ;
  if (!finalCleanup && reportArray != (Array)1 && reportArray != (Array)2) 
    arr(reportArray, a->id, Array) = 0 ;
}

/******************************/

void arrayLock (Array a)
{
  if (!arrayExists(a))
    messcrash ("arrayLock called on non exiting array") ;
  if (a->lock)
    messcrash ("arrayLock called on already locked array") ;
  a->lock = TRUE ;
  return ;
} /* arrayLock */

void arrayUnlock (Array a)
{
  if (!arrayExists(a))
    messcrash ("arrayUnlock called on non exiting array") ;
  if (!a->lock)
    messcrash ("arrayLock called on unlocked array") ;
  a->lock = FALSE ;
  return ;
} /* arrayLock */

/******************************/

#ifndef MEM_DEBUG
  void arrayExtend (Array a, int n)
#else
  void arrayExtend_dbg (Array a, int n, const char *hfname,int hlineno) 
#endif
{
  char *neuf ;

  if (!a || n < a->dim)
    return ;
  
  if (a->lock)
     messcrash ("arrayExtend called on locked array") ;
  
  if (n > 0x80000000)
    messcrash ("Array in acedb cannot exceed 2G bytes, sorry, n=%ld max=%d cell size = %d\n", (long)n, a->max, a->size) ;

  if (reportArray != (Array)2)
    totalAllocatedMemory -= a->dim * a->size ;
  if (a->dim*a->size < 1 << 23)  /* 2 megs of keys, or 8 megs of ram */
    a->dim *= 2 ;
  else
    a->dim += 1024 + ((1 << 23) / a->size) ;
  if (n >= a->dim)
    a->dim = n + 1 ;

  if (reportArray != (Array)2)
    totalAllocatedMemory += a->dim * a->size ;
#ifndef MEM_DEBUG
  neuf = (char *) messalloc (a->dim*a->size) ;
#else
  neuf = (char *) messalloc_dbg (a->dim*a->size,dbgPos(hfname, hlineno, __FILE__), __LINE__) ;
#endif
  memcpy (neuf,a->base,a->size*a->max) ;
  messfree (a->base) ;
  a->base = neuf ;
}

/***************/

char *uArray (Array a, int i)
{
  if (i < 0)
    messcrash ("referencing array element %d < 0", i) ;
  if (!a)
    messcrash ("uArray called with NULL Array struc");

  if (i >= a->max)
    { if (i >= a->dim)
        arrayExtend (a,i) ;
      a->max = i+1 ;
    }
  return a->base + ((long int)i)*a->size ;
}

/***************/

 char *uArrayCheck (Array a, int i, int size)
{
  if (! a)
    messcrash ("dereferecning a null Array") ;
  if (size != a->size)
    messcrash ("Array size mismatch accessing array of size %d with pointer if typse-size %d", a->size, size) ;
  if (i < 0)
    messcrash ("referencing array element %d < 0", i) ;

  return uArray(a, i) ;
}

/***************/

 char *uArrCheck (Array a, int i, int size)
{
  if (! a)
    messcrash ("dereferencing a null Array") ;
  if (size != a->size)
    messcrash ("Array size mismatch accessing array of size %d with pointer if typse-size %d", a->size, size) ;
  if (i >= a->max || i < 0)
    messcrash ("array index %d out of bounds [0,%d]",
	       i, a->max - 1) ;
  return a->base + ((long int)i)*a->size ;
}

/**************/

#ifndef MEM_DEBUG
 Array arrayHandleCopyExtend (Array a, int extendBy, AC_HANDLE handle)
#else
   Array arrayHandleCopyExtend_dbg(Array a,  const char *hfname,int hlineno, int extendBy, AC_HANDLE handle) 
#endif
{ Array b ;
  
  if (arrayExists (a) && a->size)
    {
#ifndef MEM_DEBUG
      b = uArrayCreate (a->max+ extendBy, a->size, handle) ;
#else
	  b = uArrayCreate_dbg (a->max+extendBy, a->size, handle,
				dbgPos(hfname, hlineno, __FILE__), __LINE__) ;
#endif
      memcpy(b->base, a->base, a->max * a->size);
      b->max = a->max ;
      return b;
    }
  else
    return 0 ;
}

/**************/

Array arrayTruncatedCopy (Array a, int x1, int x2)
{ Array b = 0 ;
  
  if (x1 < 0 || x2 < x1 || x2 > a->max)
    messcrash 
      ("Bad coordinates x1 = %d, x2 = %d in arrayTruncatedCopy",
       x1, x2) ;
  if (arrayExists (a) && a->size)
    { if (x2 - x1)
	{ b = uArrayCreate (x2 - x1, a->size, 0) ;
	  b->max = x2 - x1 ;
	  memcpy(b->base, a->base + x1, b->max * b->size);
	}
      else
	b = uArrayCreate (10, a->size, 0) ;
    }
  return b;
}

/**************/

void arrayCompress(Array a)
{
  mysize_t i, j, k , as ;
  char *x, *y, *ab  ;
  
  if (!a || !a->size || arrayMax(a) < 2 )
    return ;

  ab = a->base ; 
  as = a->size ;
  for (i = 1, j = 0 ; i < arrayMax(a) ; i++)
    { x = ab + i * as ; y = ab + j * as ;
      for (k = a->size ; k-- ;)		
	if (*x++ != *y++) 
	  goto different ;
      continue ;
      
    different:
      if (i != ++j)
	{ x = ab + i * as ; y = ab + j * as ;
	  for (k = a->size ; k-- ;)	 
	    *y++ = *x++ ;
	}
    }
  arrayMax(a) = j + 1 ;
}

/****************/

/****************/

/* 2025_09_10 , added arraySortSlice */
void arraySortSlice (Array a, mysize_t pos1, mysize_t pos2, int (*order)(const void*, const void*))
{
  mysize_t n = pos2 - pos1 ;
  
  if (! arrayExists (a))
    messcrash("arraySortSlice called on non existing array, please edit the source code") ;
  if (pos1 < 0)
    messcrash("pos1 = %ld < 0 in arraySortSlice", pos1) ;
  if (pos2 > a->max)
    messcrash("pos2 = %ld > max = %ld in arraySortSlice", pos2, a->max) ;
  
  if (n > 1) 
    {
      int s = a->size ;
      void *v = a->base + s * pos1 ;
      mSort (v, n, s, order) ;
    }
}

/* 31.7.1995 dok408  added arraySortPos() - restricted sorting to tail of array */
void arraySortPos (Array a, mysize_t  pos, int (*order)(const void*, const void*))
{ return arraySortSlice(a, pos, a ? a->max : 0, order) ;  }

void arraySort (Array a, int (*order)(const void*, const void*))
{ return arraySortSlice(a, 0, a ? a->max : 0, order) ;  }


/***********************************************************/

BOOL arrayIsEntry (Array a, int i, void *s)
{
  char *cp = (char *) uArray(a,i), *cq = (char *)s ;
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

BOOL arrayFind(Array a, void *s, int *ip, int (* order)(const void*, const void*))
{
  int ord;
  int i = 0 , j = arrayMax(a), k;

  if(!j || (ord = order(s,uArray(a,0)))<0)
    { if (ip)
	*ip = -1; 
      return FALSE;
    }   /* not found */

  if (ord == 0)
    { if (ip)
	*ip = 0;
      return TRUE;
    }

  if ((ord = order(s,uArray(a,--j)))>0 )
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
      if ((ord = order(s, uArray(a,k))) == 0)
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
}

/**************************************************************/
       /* Removes Entry s from Array  a
        * sorted in ascending order of order()
        */

BOOL arrayRemove (Array a, void * s, int (* order)(const void*, const void*))
{
  int i;

  if (arrayFind(a, s, &i,order))
    {
      /* memcpy would be faster but regions overlap
       * and memcpy is said to fail with some compilers
       */
      char *cp = uArray(a,i),  *cq = cp + a->size ;
      mysize_t j = (arrayMax(a) - i)*(a->size) ;
      while(j--)
	*cp++ = *cq++;

      arrayMax(a) --;
      return TRUE;
    }
  else

    return FALSE;
}

/**************************************************************/
       /* Insert Segment s in Array  a
        * in ascending order of s.begin
        */

BOOL arrayInsert(Array a, void * s, int (*order)(const void*, const void*))
{
  int i, j;

  if (arrayFind(a, s, &i,order))
    return FALSE;  /* no doubles */
  
  j = arrayMax(a) + 1;
  uArray(a,j-1) ; /* to create space */

	/* avoid memcpy for same reasons as above */
  {
    char* cp = uArray(a,j - 1) + a->size - 1,  *cq = cp - a->size ;
    mysize_t k = (j - i - 1)*(a->size);
    while(k--)
      *cp-- = *cq--;
    
    cp = uArray(a,i+1); cq = (char *) s; k = a->size;
    while(k--)
      *cp++ = *cq++;
  }
  return TRUE;
}

/********* Stack : arbitrary Stack class - inherits from Array **********/

static void uStackFinalise (void *cp) ;

#ifndef MEM_DEBUG
Stack stackHandleCreate (int n, AC_HANDLE handle)  /* n is initial size */
{
  Stack s = (Stack) handleAlloc (uStackFinalise, 
				 handle,
				 sizeof (struct StackStruct)) ;
#else
Stack stackHandleCreate_dbg (int n, AC_HANDLE handle,  /* n is initial size */
					      const char *hfname, int hlineno)
{
  Stack s = (Stack) handleAlloc_dbg (uStackFinalise, 
				 handle,
				 sizeof (struct StackStruct),
				 hfname, hlineno) ;
#endif
  s->magic = STACK_MAGIC ;
  s->a = arrayCreate (n,char) ;
  s->pos = s->ptr = s->a->base ;
  s->safe = s->a->base + s->a->dim - 16 ; /* safe to store objs up to size 8 */
  s->pushPop = s->textOnly = FALSE ;  /* clean slate */
  return s ;
}

Stack stackReCreate (Stack s, int n)               /* n is initial size */
{
  if (!stackExists(s))
    return stackCreate(n) ;

  s->a = arrayReCreate (s->a,n,char) ;
  s->pos = s->ptr = s->a->base ;
  s->safe = s->a->base + s->a->dim - 16 ; /* safe to store objs up to size 8 */
  s->pushPop = s->textOnly = FALSE ; /* clean slate */
  return s ;
}

Stack stackCopy (Stack old, AC_HANDLE handle)
{
  Stack neuf = 0 ;

  if (stackExists(old))
    {
      neuf = stackHandleCreate (old->a->dim, handle) ;
      memcpy (neuf->a->base, old->a->base, old->a->dim) ;
      neuf->ptr = neuf->a->base + (old->ptr - old->a->base) ;
      neuf->pos = neuf->a->base + (old->pos - old->a->base) ;
      neuf->textOnly = old->textOnly ;
      neuf->pushPop = old->pushPop ;
    }
  return neuf ;
}

Stack arrayToStack (Array a)
{ 
  Stack s ;
  int n ;

  if (!arrayExists(a) || a->size != 1 )
    return 0 ;
    
  n = arrayMax(a) ;
  s = stackCreate(n  + 32) ;
              
  memcpy(s->a->base, a->base, n) ;
                
  s->pos = s->a->base ;
  s->ptr = s->a->base + n ;
  s->safe = s->a->base + s->a->dim - 16 ; /* safe to store objs up to size 8 */
  while ((long)s->ptr % STACK_ALIGNMENT)
    *(s->ptr)++ = 0 ;   
  return s ;
}

void uStackDestroy(Stack s)
{ if (s && s->magic == STACK_MAGIC) messfree(s);
} /* the rest is done below as a consequence */

static void uStackFinalise (void *cp)
{ Stack s = (Stack)cp;
  if (!finalCleanup) arrayDestroy (s->a) ;
  s->magic = 0 ;
}

void stackExtend (Stack s, int n)
{
  int ptr = s->ptr - s->a->base,
      pos = s->pos - s->a->base ;
  s->a->max = s->a->dim ;	/* since only up to ->max copied over */
  arrayExtend (s->a,ptr+n+16) ;	/* relies on arrayExtend mechanism */
  s->ptr = s->a->base + ptr ;
  s->pos = s->a->base + pos ;
  s->safe = s->a->base + s->a->dim - 16 ;
}

int stackMark (Stack s)
{ return (s->ptr - s->a->base) ;
}

int stackPos (Stack s)
{ return (s->pos - s->a->base) ;
}

void stackCursor (Stack s, int mark)
{ s->pos = s->a->base + mark ;
}


/* access doubles on the stack by steam on systems where we don't
   align the stack strictly enough, this code assumes that ints are 32bits
   and doubles 64 */

union mangle { double d;
	       struct { int i1;
			int i2;
		      } i;
	     };

void ustackDoublePush(Stack stk, double x)
{ union mangle m;

  stk->pushPop = TRUE ;
  if (stk->textOnly)
    messcrash("pushDouble called on textOnly stack\n");

  m.i.i1 = m.i.i2 = 0 ; /* for compiler happiness */
  m.d = x;
  push(stk, m.i.i1, int);
  push(stk, m.i.i2, int);
}

double ustackDoublePop(Stack stk)
{ union mangle m;

  stk->pushPop = TRUE ;
  if (stk->textOnly)
    messcrash("popDouble called on textOnly stack\n");

  m.i.i2 = pop(stk, int);
  m.i.i1 = pop(stk, int);

  return m.d;
}

double ustackDoubleNext(Stack stk)
{ union mangle m;

  stk->pushPop = TRUE ;
  if (stk->textOnly)
    messcrash("nextDouble called on textOnly stack\n");

  m.i.i1 = stackNext(stk, int);
  m.i.i2 = stackNext(stk, int);

  return m.d;
}

int pushText (Stack s, const char* text)
{
  int n = stackMark (s) ;

  s->textOnly = TRUE ;
  if (s->pushPop)
    messcrash("pushText called on pushPop stack\n");

  while (s->ptr + strlen(text)  > s->safe) 
    stackExtend (s,strlen(text)+1) ;
  while ((*(s->ptr)++ = *text++)) ;
#ifdef MEM_DEBUG
  if (!s->textOnly)
    messcrash("pushText on non-text stack\n");
#else
  /* old code preserved in case something depends on it */
  if (!s->textOnly)
    while ((long)s->ptr % STACK_ALIGNMENT)
      *(s->ptr)++ = 0 ;   
#endif
  return n ;
}

char* popText (Stack s)
{
  char *base = s->a->base ;

  s->textOnly = TRUE ;
  if (s->pushPop)
    messcrash("popText called on pushPop stack\n");

  while (s->ptr > base && !*--(s->ptr)) ;
  while (s->ptr >= base && *--(s->ptr)) ;
  return ++(s->ptr) ;
}

void catText (Stack s, const char* text)
{
  s->textOnly = TRUE ;
  if (s->pushPop)
    messcrash("catText on pushPop stack\n");
  while (s->ptr + strlen(text) > s->safe)
    stackExtend (s,strlen(text)+1) ;
  *s->ptr = 0 ;
  while (s->ptr >= s->a->base && *s->ptr == 0)
    s->ptr -- ;
  s->ptr ++ ;
  while ((*(s->ptr)++ = *text++)) ;
  /* old code preserved in case something depends on it */
  if (!s->textOnly)
    while ((long)s->ptr % STACK_ALIGNMENT)
      *(s->ptr)++ = 0 ;   
}

void catBinary (Stack s, const void *data, int size)
{
  int total;
  total = size + 1;

  s->textOnly = TRUE ;
  if (s->pushPop)
    messcrash("catBinary on pushPop stack\n");
  while (s->ptr + total > s->safe)
    stackExtend (s,size+1) ;
  memcpy(s->ptr,data,size);
  s->ptr += size;
}

char* stackNextText (Stack s)
{ 
  char *text = s->pos ;

  s->textOnly = TRUE ;
  if (s->pushPop)
    messcrash("stackNextText called on pushPop stack\n");

  if (text>= s->ptr)
    return 0 ;  /* JTM, so while stackNextText makes sense */
  while (*s->pos++) ;
  /* old code preserved in case something depends on it */
  if (!s->textOnly)
    while ((long)s->pos % STACK_ALIGNMENT)
      ++s->pos ;
  return text ;
}

/*********/
     /* Push text in stack s, after breaking it on delimiters */
     /* You can later access the tokens with command
	while (token = stackNextText(s)) work on your tokens ;
	*/
void  stackTokeniseTextOn(Stack s, char *text, char *delimiters)
{
  char *cp, *cq , *cend, *cd, old, oldend ;
  int i, n ;

  if(!stackExists(s) || !text || !delimiters)
    messcrash("stackTextOn received some null parameter") ;

  s->textOnly = TRUE ;
  if (s->pushPop)
    messcrash("stackTokenizeTextOn called on pushPop stack\n");

  n = strlen(delimiters) ;
  cp = cq  = text ;
  while(TRUE)
    {
      while(*cp == ' ')
	cp++ ;
      cq = cp ;
      old = 0 ;
      while(*cq)
	{ for (cd = delimiters, i = 0 ; i < n ; cd++, i++)
	    if (*cd == *cq)
	      { old = *cq ;
		*cq = 0 ;
		goto found ;
	      }
	  cq++ ;
	}
    found:
      cend = cq ;
      while(cend > cp && *--cend == ' ') ;
      if (*cend != ' ') cend++ ;
      oldend = *cend ; *cend = 0 ;
      if (*cp && cend > cp)
	pushText(s,cp) ;
      *cend = oldend ;
      if(!old)
	{ stackCursor(s, 0) ;
	  return ;
	}
      *cq = old ;
      cp = cq + 1 ;
    }
}

void stackClear(Stack s)
{ if (stackExists(s))
    { s->pos = s->ptr = s->a->base;
      s->a->max = 0 ;
      s->pushPop = s->textOnly = FALSE ; /* clean slate */
    }
}

/************************  end of file ********************************/
/**********************************************************************/
 
 
 
 
