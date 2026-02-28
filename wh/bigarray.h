/*  File: bigarray.h
 *  Author: Richar Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: header for bigarraysub.c
 *              Allow larger Arrays, indexed by long int
 *              I contemplated modifying Array to use long but would imply
 *              editing the code in a zillion places and would slow everything
 *              when bigarrays are only useful in particular situations
 *              NOT to be included by the user, included by regular.h and ac.h
 * Exported functions:
 *              the BigArray type and associated functions
 * HISTORY:
 * Created: Aug 24 2014 (mieg)
 *-------------------------------------------------------------------
 */

#ifndef DEF_BIG_ARRAY_H
#define DEF_BIG_ARRAY_H
#include "regular.h"  /* may set ARRAY_CHECK and MEM_DEBUG*/ 
mysize_t stackused (void) ;
 
/************* Array package ********/

typedef struct BigArrayStruct
  { char* base ;    /* char* since need to do pointer arithmetic in bytes, always aligned on 64bytes boundary */
    unsigned long int dim ;     /* length of alloc'ed space */
    int   size ;            /* length of a cell */
    long int  max ;     /* largest element accessed via array() */
    mysize_t   id ;      /* unique identifier */
    int   magic ;
    BOOL lock ; /* a locked array cannot be destroyed or reallocated */
    long int readOnly ; /* original size, cannot be modified */
    void *map ;         /* address to be unmapped */
    char *trueBase ;    /* true unaligned alloc, can be freed */
    const char *fName ;
    int fd ; /* file descriptor */
  } *BigArray ;
 
    /* NB we need the full definition for arr() for macros to work
       do not use it in user programs - it is private.
    */

#define BIG_ARRAY_MAGIC 2914574
#define BIGSTACK_MAGIC 4318275
#if !defined(MEM_DEBUG)
  BigArray   uBigArrayCreate (long int n, int size, AC_HANDLE handle) ;
  void    bigArrayExtend (BigArray a, long int  n) ;
  BigArray bigArrayCopy (BigArray a) ;
  BigArray bigArrayHandleCopy (BigArray a, AC_HANDLE handle) ;
#else
  BigArray   uBigArrayCreate_dbg (long int  n, int size, AC_HANDLE handle,
			    const char *hfname,int hlineno) ;
  void    bigArrayExtend_dbg (BigArray a, long int n, const char *hfname,  int hlineno) ;
  BigArray	bigArrayCopy_dbg(BigArray a, const char *hfname,int hlineno) ; 
  BigArray	bigArrayHandleCopy_dbg(BigArray a, const char *hfname,int hlineno, AC_HANDLE handle) ; 
#define bigArrayCopy(a) bigArrayCopy_dbg(a, __FILE__, __LINE__)
#define uBigArrayCreate(n, s, h) uBigArrayCreate_dbg(n, s, h, __FILE__, __LINE__)
#define bigArrayExtend(a, n ) bigArrayExtend_dbg(a, n, __FILE__, __LINE__)
#define bigArrayHandleCopy(a,h) bigArrayHandleCopy_dbg(a, __FILE__, __LINE__, h)

#endif

BOOL bigArrayMapWrite (BigArray aa, const char *fName) ;
BigArray uBigArrayMapRead (const char *fName, int size, BOOL readOnly, AC_HANDLE h) ;
#define bigArrayMapRead(fName,type,readOnly,h) uBigArrayMapRead(fName,sizeof(type),readOnly,h) 
BigArray   uBigArrayReCreate (BigArray a,long int n, int size) ;
void    uBigArrayDestroy (BigArray a);
char    *uBigArray (BigArray a, long int index) ;
char    *uBigArrCheck (BigArray a, long int index, int size) ;
char    *uBigArrayCheck (BigArray a, long int index, int size) ;
#define bigArrayCreate(n,type)	uBigArrayCreate(n,sizeof(type), 0)
#define bigArrayHandleCreate(n,type,handle) uBigArrayCreate(n, sizeof(type), handle)
#define bigArrayReCreate(a,n,type)	uBigArrayReCreate(a,n,sizeof(type))
#define bigArrayDestroy(x)		((x) ? uBigArrayDestroy(x), x=0, TRUE : FALSE)

#if (defined(ARRAY_CHECK))
#define bigArrp(ar,i,type)	((type*)uBigArrCheck(ar,i,sizeof(type)))
#define bigArr(ar,i,type)	(*(type*)uBigArrCheck(ar,i,sizeof(type)))
#define bigArrayp(ar,i,type)	((type*)uBigArrayCheck(ar,i,sizeof(type)))
#define bigArray(ar,i,type)	(*(type*)uBigArrayCheck(ar,i,sizeof(type)))
#else
#define bigArr(ar,i,type)	((*(type*)((ar)->base + ((long int)i)*(ar)->size)))
#define bigArrp(ar,i,type)	(((type*)((ar)->base + ((long int)i)*(ar)->size)))
#define bigArrayp(ar,i,type)	((type*)uBigArray(ar,i))
#define bigArray(ar,i,type)	(*(type*)uBigArray(ar,i))
#endif /* ARRAY_CHECK */

            /* only use arr() when there is no danger of needing expansion */
BigArray   bigArrayTruncatedCopy (BigArray a, long int x1, long int x2) ;
void    bigArrayStatus  (mysize_t *nmadep, mysize_t *nusedp, long int *memAllocp, long int *memUsedp) ;
mysize_t    bigArrayReportMark (void) ; /* returns current array number */
void    bigArrayReport (mysize_t j) ;	/* write stderr about all arrays since j */
#define bigArrayMax(ar)            ((ar)->max)
#define bigArrayForceFeed(ar,j) (uBigArray(ar,j), (ar)->max = (j))
#define bigArrayExists(ar)		((ar) && (ar)->magic == BIG_ARRAY_MAGIC ? (ar)->id : 0 ) 
            /* JTM's package to hold sorted arrays of ANY TYPE */
/*
BOOL    arrayInsert(BigArray a, void * s, int (*order)(const void*, const void*));
BOOL    arrayRemove(BigArray a, void * s, int (*order)(const void*, const void*));
*/
void bigArrayLock (BigArray a) ;
void bigArrayUnlockLock (BigArray a) ;
/* sort the whole array */
void    bigArraySort(BigArray a, int (*order)(const void*, const void*)) ;
/* sort from pos to end */
void    bigArraySortPos (BigArray a, long int pos, int (*order)(const void*, const void*));
/* sort from pos1 to pos 2, pos2 excluded [pos1, pos2[ */
void    bigArraySortSlice (BigArray a, long int pos1, long int pos2, int (*order)(const void*, const void*)); 
void    bigArrayCompress(BigArray a) ;
BOOL    bigArrayFind(BigArray a, void *s, long int *ip, int (*order)(const void*, const void*));
BOOL    bigArrayIsEntry(BigArray a, long int i, void *s);
void bigMSort (void *b, long int n, int s, int (*cmp)(const void *va, const void *vb)) ; 

/* ATTENTION: to  optimize the computation of huge lists of number array a is double-barrelled
 * for i < LL  bigArray(a, i, float) contains the number of times i has been seen
 * for i >= LL bigArray(a, i, float) enumerates the values
 * this is a very strong optimisation when measuring huge SNP tables
 */
int bigFloatVariance (BigArray a, int LL, float *medianp, float *averagep, float *sigmap) ;

/************** BigStack package **************/

/*
* BigStack only implements dynamically sized strings.  In the interface,
* this looks like a C string, but the actual content of the string
* is copied into the bigStack.
*
* This would be just like any other data item except for one thing:
* If the thing on the top of the bigStack is a Text, you can concatenate
* additional text onto it.  That is, after
* 	pushText(s,"a");
* 	catText(s,"b");
* 	catText(s,"c");
* 	cp = popText(s)
* 	assert( strcmp(cp,"abc") == 0 )
* 
* Another feature is a cursor.  You can set the cursor to an
* arbitrary position in the bigStack, and then you can peek at what
* is there.  You can also walk the cursor toward the top of the
* bigStack, so you can examine the bigStack content in the order it was
* pushed.
*
*
* Another common use of the bigStack is as a place to store strings
* without explicitly writing many malloc calls.
*	n = bigStackMark(s);
*	pushText(s,"xyzzy");
*	// do something
*	printf("%s\n",bigStackText(s,n));
*
*
* Creating a BigStack:
*	s = bigStackHandleCreate( initial size in bytes, handle )
*		makes an empty bigStack on a handle
*
*	s = bigStackReCreate( bigStack, initial size in bytes )
*		The existing bigStack is emptied.  If the memory allocated to
*		the bigStack is more than 500k, the data area will be freed
*		and re-created.
*
*	s = bigStackCopy( bigStack )
*		make another bigStack just like one that you have
*
* Destroying bigStacks:
*	bigStackDestroy( bigStack )
*
* Misc functions:
*	bigStackExtend( bigStack, size in bytes )
*		makes the bigStack bigger if necessary.  Since the bigStack
*		grows automatically, you don't really need this, but 
*		if you know you have a lot of data coming in, you can
*		grow it all at once.
*
*	bigStackClear( bigStack )
*		Makes the bigStack empty.  This is just like bigStackReCreate,
*		except it will not free memory.
*
*	bigStackEmpty( bigStack )
*		returns true if the bigStack contains 0 bytes of data
*
*	bigStackExists( bigStack )
*		returns true if the bigStack passed in appears to be
*		correctly initialized (i.e. has been created and
*		not destroyed)
*
* Cursor for peeking into the bigStack:
*	There is one cursor in each bigStack.  It can be set to an arbitrary
*	position in the bigStack, and then moved toward the top of the bigStack
*	as you examine data in the bigStack.  None of these functions modify
*	the data in the bigStack.
*
*	bigStackMark( bigStack )
*		Returns an int that describes the current position of the top
*		of the bigStack.  You can pass this value to bigStackCursor later.
*
*	bigStackPos( bigStack )
*		returns an int that describes the current position of the
*		cursor in the bigStack.  You can pass this value to bigStackCursor 
*		later.
*
*	bigStackCursor( bigStack, position )
*		set the cursor to a particular position in the bigStack.  The
*		only values you can really trust are those returned from
*		bigStackMark() or bigStackPos(), or 0 which is the bottom of
*		the bigStack.
*
*	bigStackAtEnd( bigStack )
*		returns true if the cursor has reached the top of the
*		bigStack.
*
*	bigStackText( bigStack, position )
*		returns the text at the indicated position in the bigStack.
*		This does not modify the cursor.
*
*/
 
typedef struct BigStackStruct      /* assumes objects <= 16 bytes long */
  { BigArray a ;
    int magic ;
    char* ptr ;         /* current end pointer */
    char* pos ;         /* potential internal pointer */

    char* safe ;        /* need to extend beyond here */
  } *BigStack ;
 
        /* as with ArrayStruct, the user should NEVER access BigStackStruct
           members directly - only through the subroutines/macros
        */
#if !defined(MEM_DEBUG)
  BigStack   bigStackHandleCreate (long int n, AC_HANDLE handle) ;
#else
  BigStack   bigStackHandleCreate_dbg (long int n, AC_HANDLE handle,
				 const char *hfname,int hlineno) ;
#define bigStackHandleCreate(n, h) bigStackHandleCreate_dbg(n, h, __FILE__, __LINE__)
#endif

BigStack   bigStackReCreate (BigStack s, long int n) ;
BigStack   bigStackCopy (BigStack, AC_HANDLE handle) ;

void    uBigStackDestroy (BigStack s);
#define bigStackDestroy(x)	 ((x) ? uBigStackDestroy(x), (x)=0, TRUE : FALSE)
void    bigStackExtend (BigStack s, long int n) ;
void    bigStackClear (BigStack s) ;
#define bigStackEmpty(stk)  ((stk)->ptr <= (stk)->a->base)
#define bigStackExists(stk) ((stk) && (stk)->magic == BIGSTACK_MAGIC ? bigArrayExists((stk)->a) : 0)

long int     bigStackMark (BigStack s) ;              /* returns a mark of current ptr */
long int     bigStackPos (BigStack s) ;              /* returns a mark of current pos, useful with bigStackNextText */
void    bigStackCursor (BigStack s, long int mark) ;  /* sets ->pos to mark */
 
#define bigStackText(stk,mark) ((char*)((stk)->a->base + (mark)))

#endif   /* BIG_ARRAY_DEF */
