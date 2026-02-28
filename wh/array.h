/*  File: array.h
 *  Author: Richar Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: header for arraysub.c
 *              NOT to be included by the user, included by regular.h
 * Exported functions:
 *              the Array type and associated functions
 *              the Stack type and associated functions
 *              the Associator functions
 * HISTORY:
 * Last edited: Dec  4 11:03 1998 (fw)
 * Created: Fri Dec  4 11:01:35 1998 (fw)
 *-------------------------------------------------------------------
 */

#ifndef DEF_ARRAY_H
#define DEF_ARRAY_H
#include "regular.h"  /* may set ARRAY_CHECK */
 
mysize_t stackused (void) ;
 
/************* Array package ********/

/* #define ARRAY_CHECK either here or in a single file to
   check the bounds on arr() and arrp() calls
*/
typedef struct { KEY k ; float x, y ;} POINT2D ;

typedef struct ArrayStruct
  { char* base ;    /* char* since need to do pointer arithmetic in bytes */
    unsigned long int dim ;     /* length of alloc'ed space */
    int   size ;
    unsigned int  max ;     /* largest element accessed via array() */
    mysize_t   id ;      /* unique identifier */
    int   magic ;
    BOOL lock ;
  } *Array ;
 
    /* NB we need the full definition for arr() for macros to work
       do not use it in user programs - it is private.
    */

#define ARRAY_MAGIC 8918274
#define STACK_MAGIC 8918275
#define   ASS_MAGIC 8918276

#if !defined(MEM_DEBUG)
  Array   uArrayCreate (mysize_t n, int size, AC_HANDLE handle) ;
  void    arrayExtend (Array a, int n) ;
  void    arrayExtend (Array a, int n) ;
Array arrayHandleCopyExtend(Array a, int extendBy, AC_HANDLE h) ;
#define arrayCopy(a) arrayHandleCopyExtend(a,0,0)
#define arrayHandleCopy(a,h) arrayHandleCopyExtend(a, 0,h)
#define dnaCopy(_a) arrayHandleCopyExtend(_a,1,0)
#define dnaHandleCopy(_a,h) arrayHandleCopyExtend(_a,1,h)
#else
  Array   uArrayCreate_dbg (mysize_t n, int size, AC_HANDLE handle,
			    const char *hfname,int hlineno) ;
  void    arrayExtend_dbg (Array a, int n, const char *hfname,int hlineno) ;
  Array	arrayCopy_dbg(Array a, const char *hfname,int hlineno) ; 
Array	arrayHandleCopyExtend_dbg(Array a, const char *hfname,int hlineno, int extendBy, AC_HANDLE handle) ; 
#define uArrayCreate(n, s, h) uArrayCreate_dbg(n, s, h, __FILE__, __LINE__)
#define arrayExtend(a, n ) arrayExtend_dbg(a, n, __FILE__, __LINE__)
#define arrayCopy(a) arrayHandleCopyExtend_dbg(a, __FILE__, __LINE__,0,0)
#define arrayHandleCopy(a,h) arrayHandleCopyExtend_dbg(a, __FILE__, __LINE__, 0, h)
#define dnaCopy(_dna) arrayHandleCopyExtend_dbg(_dna, __FILE__, __LINE__,1,0)
#define dnaHandleCopy(_dna,h) arrayHandleCopyExtend_dbg(_dna, __FILE__, __LINE__, 1,h)

#endif
void arrayLock (Array a) ;
void arrayUnlock (Array a) ;

Array   uArrayReCreate (Array a, int n, int size) ;
void    uArrayDestroy (Array a);
char    *uArray (Array a, int index) ;
char    *uArrCheck (Array a, int index, int size) ;
char    *uArrayCheck (Array a, int index, int size) ;
#define arrayCreate(n,type)	uArrayCreate(n,sizeof(type), 0)
#define arrayHandleCreate(n,type,handle) uArrayCreate(n, sizeof(type), handle)
#define arrayReCreate(a,n,type)	uArrayReCreate(a,n,sizeof(type))
#define arrayDestroy(x)		((x) ? uArrayDestroy(x), x=0, TRUE : FALSE)

#if (defined(ARRAY_CHECK))
#define arrp(ar,i,type)	((type*)uArrCheck(ar,i,sizeof(type)))
#define arr(ar,i,type)	(*(type*)uArrCheck(ar,i,sizeof(type)))
#define arrayp(ar,i,type)	((type*)uArrayCheck(ar,i,sizeof(type)))
#define array(ar,i,type)	(*(type*)uArrayCheck(ar,i,sizeof(type)))
#else
#define arr(ar,i,type)	((*(type*)((ar)->base + ((long int)(i))*(ar)->size)))
#define arrp(ar,i,type)	(((type*)((ar)->base + ((long int)(i))*(ar)->size)))
#define arrayp(ar,i,type)	((type*)uArray(ar,i))
#define array(ar,i,type)	(*(type*)uArray(ar,i))
#endif /* ARRAY_CHECK */

            /* only use arr() when there is no danger of needing expansion */
Array   arrayTruncatedCopy (Array a, int x1, int x2) ;
void    arrayStatus (mysize_t *nmadep, mysize_t *nusedp, mysize_t *memAllocp, mysize_t *memUsedp) ;
mysize_t     arrayReportMark (void) ; /* returns current array number */
void    arrayReport (mysize_t j) ;	/* write stderr about all arrays since j */
#define arrayMax(ar)            ((ar)->max)
#define arrayForceFeed(ar,j) (uArray(ar,j), (ar)->max = (j))
#define arrayExists(ar)		((ar) && (ar)->magic == ARRAY_MAGIC ? (ar)->id : 0 ) 
            /* JTM's package to hold sorted arrays of ANY TYPE */
BOOL    arrayInsert(Array a, void * s, int (*order)(const void*, const void*));
BOOL    arrayRemove(Array a, void * s, int (*order)(const void*, const void*));
/* sort the whole array */
void    arraySort(Array a, int (*order)(const void*, const void*)) ;
/* sort from pos to end */
void    arraySortPos (Array a, mysize_t pos, int (*order)(const void*, const void*));
/* sort from pos1 to pos 2, pos2 excluded [pos1, pos2[ */
void    arraySortSlice (Array a, mysize_t pos1, mysize_t pos2, int (*order)(const void*, const void*));
void    arrayCompress(Array a) ;
BOOL    arrayFind(Array a, void *s, int *ip, int (*order)(const void*, const void*));
BOOL    arrayIsEntry(Array a, int i, void *s);
void mSort (void *b, mysize_t n, int s, int (*cmp)(const void *va, const void *vb)) ;

/* average and variance of an array of float, stable algo by Knuth */
int floatVariance (Array a, float *medianp, float *averagep, float *sigmap) ;
int doubleVariance (Array a, double *medianp, double *averagep, double *sigmap) ;
BOOL diVariance (Array a, float *splitp, float *avp, float *sigmap, float *av1p, float *sigma1p, float *av2p, float *sigma2p, float *sigmaCombinedp) ;
BOOL linearRegression (Array xy, double *ap, double *bp, double *rp, double *wp) ;

/************** Stack package **************/

/*
* Stack really implements three different major features in the same package.
*
* The first feature is a conventional stack.   If you only use push()
* and pop() to manipulate the stack, it behaves as you might expect.
* The caveat is that when you pop a value, you must know what it
* is going to be.  The pop() macro just picks up the right number
* of bytes and hands them to you.
* 
* The second feature is a dynamically sized string.  One of the
* special data types used by the stack is "Text".  In the interface,
* this looks like a C string, but the actual content of the string
* is copied into the stack.
*
* This would be just like any other data item except for one thing:
* If the thing on the top of the stack is a Text, you can concatenate
* additional text onto it.  That is, after
* 	pushText(s,"a");
* 	catText(s,"b");
* 	catText(s,"c");
* 	cp = popText(s)
* 	assert( strcmp(cp,"abc") == 0 )
* 
* In fact, I see this feature used far more often that I ever see
* stack-like functions.
*
* If you use "Text" in the stack, do not also use other types of
* values.
*
* The third feature is a cursor.  You can set the cursor to an
* arbitrary position in the stack, and then you can peek at what
* is there.  You can also walk the cursor toward the top of the
* stack, so you can examine the stack content in the order it was
* pushed.
*
*
* Another common use of the stack is as a place to store strings
* without explicitly writing many malloc calls.
*	n = stackMark(s);
*	pushText(s,"xyzzy");
*	// do something
*	printf("%s\n",stackText(s,n));
*
*
* Creating a Stack:
*	s = stackCreate( initial size in bytes )
*		makes an empty stack
*
*	s = stackHandleCreate( initial size in bytes, handle )
*		makes an empty stack on a handle
*
*	s = stackReCreate( stack, initial size in bytes )
*		The existing stack is emptied.  If the memory allocated to
*		the stack is more than 500k, the data area will be freed
*		and re-created.
*
*	s = stackCopy( stack )
*		make another stack just like one that you have
*
* TextOnly stacks:
*
*	Stacks must be used either for texts, using
*                 pushText, catText, stackText
*              OR as push-pop stacks for fundamental types of fixed size
*                 push
*		The library enforces this restriction.
*
* Destroying stacks:
*	stackDestroy( stack )
*
* Misc functions:
*	stackExtend( stack, size in bytes )
*		makes the stack bigger if necessary.  Since the stack
*		grows automatically, you don't really need this, but 
*		if you know you have a lot of data coming in, you can
*		grow it all at once.
*
*	stackClear( stack )
*		Makes the stack empty.  This is just like stackReCreate,
*		except it will not free memory.
*
*	stackEmpty( stack )
*		returns true if the stack contains 0 bytes of data
*
*	stackExists( stack )
*		returns true if the stack passed in appears to be
*		correctly initialized (i.e. has been created and
*		not destroyed)
*
* Stack functions that implement an actual stack:
*	push( stack, value, type )
*		push a value onto the top of the stack.  this is a macro.
*
*	pop(stack, type)
*		pops a value from the stack and returns it.  this is a macro.
*
*	push and pop really compute a pointer and then use *(type)pointer to
*	access the stack data.  That means you cannot push or pop any data
*	item that cannot be aligned on STACK_ALIGNMENT on the local processor.
*
*	It appears that on at least one machine supported by acedb, the
*	stack alignment used is not suitable for storing doubles.  I guess
*	it is either something about the Alpha or trying to keep doubles
*	aligned on the Pentium external data bus to avoid making two bus
*	transactions.  So, there are:
*
*	pushDouble( stack, double ), pushPtr (stack, void *)
*	popDouble( stack ), popPtr (stack, void *)
*		just like push/pop, but with double values and alignments
*
* Text functions:
*	pushText( stack, char * )
*		copy the actual bytes of the string, including the \0 at
*		the end, onto the top of the stack.  You should only
*		use this with TextOnly stacks.
*
*	popText( stack )
*		pops a Text from the stack; returns char *.  The pointer
*		returned actually points into the stack memory.  The contents
*		become invalid when you next alter the stack in any way.
*		You should only use this with TextOnly stacks.
*
*	catText( stack, char *)
*		appends bytes to a Text that is already on top of the stack.
*			pushText(s, "abcd");
*			catText(s, "efgh");
*			cp = poptext(s);
*			assert( strcmp(cp,"abcdefgh") == 0);
*		This function should work fine on an empty stack, though
*		the code often says
*			pushText(s,"");
*		before using catText()
*
*		This function pops \0 bytes off the stack before appending
*		the data.  It keeps going until it finds something that
*		is not \0 and then appends there.  This means that you
*		should never use catText after a catBinary that may contain
*		\0 at the end of the data.  For example, in this code the 
*		data inserted by catBinary is eliminated completely:
*			pushText(s,"a");
*			catBinary(s,"\0\0",2);
*			catText(s,"b");
*
*	catBinary (Stack s, char* data, mysize_t size) ;
*		Like catText, except that it takes a size for
*		the string to be inserted so the string can contain \0.
*
*		catBinary is only useable with TextOnly stacks, no checks.
*
*		To extract your data from a TextOnly stack that has been
*		manipulated with catBinary, use this code:
*			mysize_t n = stackMark(s);
*			char *s = stackText(s, 0);
*			n = n - 1;
*			memcpy(my_buffer, s, n);
*
*	stackTokeniseTextOn( stack, string, delimiters )
*		The string you pass in is split into words at characters found
*		in the string delimiters.  The resulting words are pushed on
*		the stack in the order they exist in the string.  Multiple 
*		delimiters are allowed between words.  That is, "a\t\tb" 
*		tokenizes into two words, "a", "b", not "a", "", "b".
*
*		The delimiters in the string passed in are written over 
*		with \0.
*
*		The expectation is that you would use something like this to
*		extract words from the resulting stack:
*			stackCursor(s, 0);
*			while ((cp = stackNextText(s)))
*				printf("->%s<-\n",cp);
*
* Cursor for peeking into the stack:
*	There is one cursor in each stack.  It can be set to an arbitrary
*	position in the stack, and then moved toward the top of the stack
*	as you examine data in the stack.  None of these functions modify
*	the data in the stack.
*
*	stackMark( stack )
*		Returns an int that describes the current position of the top
*		of the stack.  You can pass this value to stackCursor later.
*
*	stackPos( stack )
*		returns an int that describes the current position of the
*		cursor in the stack.  You can pass this value to stackCursor 
*		later.
*
*	stackCursor( stack, position )
*		set the cursor to a particular position in the stack.  The
*		only values you can really trust are those returned from
*		stackMark() or stackPos(), or 0 which is the bottom of
*		the stack.
*
*	stackAtEnd( stack )
*		returns true if the cursor has reached the top of the
*		stack.
*
*	stackSkip( stack)
*	stackNext( stack, type )
*	stackNextText( stack )
*	stackDoubleNext( stack )
*		returns the data item at the current cursor position 
*		and moves the cursor up to the next data item.
*
*	stackText( stack, position )
*		returns the text at the indicated position in the stack.
*		This does not modify the cursor.
*
*	stackTextForceFeed( stack , size )
*		NOTE: This is an older function.  A better way to do 
*		the same thing is:
*			s = stackReCreate(s, length);
*			catBinary(s, answer, length );
*		
*		Grow the stack to the indicated size, and make it look like that
*		much data has been place in the stack.  This function does NOT
*		put any data in the stack -- it just moves all the pointers
*		around.  You can see it used like this:
*			s = stackReCreate (s,length) ;
*			stackTextForceFeed (s, length) ;
*			memcpy (stackText(s,0), answer, length) ;
*		That is, the code knows how Stack is implemented and is directly
*		manipulating the stack internals.
*
*
*/
 
typedef struct StackStruct      /* assumes objects <= 16 bytes long */
  { Array a ;
    int magic ;
    char* ptr ;         /* current end pointer */
    char* pos ;         /* potential internal pointer */

    char* safe ;        /* need to extend beyond here */
    BOOL  pushPop ;     /* push/pop stack, incompatible with textStack */
    BOOL  textOnly; /* If this is set, don't align the stack.
		       This (1) save space (esp on ALPHA) and
		       (2) provides stacks which can be stored and got
		       safely between architectures. Once you've set this,
		       using stackTextOnly() only pushText, popText, etc, 
		       no other types. */   
  } *Stack ;
 
        /* as with ArrayStruct, the user should NEVER access StackStruct
           members directly - only through the subroutines/macros
        */
#if !defined(MEM_DEBUG)
  Stack   stackHandleCreate (int n, AC_HANDLE handle) ;
#else
  Stack   stackHandleCreate_dbg (int n, AC_HANDLE handle,
				 const char *hfname,int hlineno) ;
#define stackHandleCreate(n, h) stackHandleCreate_dbg(n, h, __FILE__, __LINE__)
#endif

#define stackCreate(n) stackHandleCreate(n, 0)
Stack   stackReCreate (Stack s, int n) ;
Stack   stackCopy (Stack, AC_HANDLE handle) ;

void    uStackDestroy (Stack s);
#define stackDestroy(x)	 ((x) ? uStackDestroy(x), (x)=0, TRUE : FALSE)
void    stackExtend (Stack s, int n) ;
void    stackClear (Stack s) ;
#define stackEmpty(stk)  ((stk)->ptr <= (stk)->a->base)
#define stackExists(stk) ((stk) && (stk)->magic == STACK_MAGIC ? arrayExists((stk)->a) : 0)


/* Stack alignment: we use two strategies: the smallest type we push is
   a short, so if the required alignment is to 2 byte boundaries, we 
   push each type to its size, and alignments are kept.

   Otherwise, we push each type to STACK_ALIGNMENT, this ensures 
   alignment but can waste space. On machines with 32 bits ints and
   pointers, we make satck alignment 4 bytes, and do the consequent unaligned
   access to doubles by steam.

   Characters and strings are aligned separately to STACK_ALIGNMENT.
*/ 


/* mieg 2104_01_13, in the pop messcrash i change  *((type*)0 to  *((type*)1 to remove compiler warnings
 * it makes no difference since messcrash stops the program before dereferencing 1
 */

#if (STACK_ALIGNMENT<=2)
#define push(stk,x,type) ((stk)->ptr < (stk)->safe ? \
                           ( *(type *)((stk)->ptr) = (x) , (stk)->ptr += sizeof(type)) : \
			    (stackExtend (stk,16), \
			     *(type *)((stk)->ptr) = (x) , (stk)->ptr += sizeof(type)) )
#define pop(stk,type)    (  ((stk)->ptr -= sizeof(type)) >= (stk)->a->base ? \
			    *((type*)((stk)->ptr)) : \
                          (messcrash ("User stack underflow"), *((type*)1)) )
#define stackNext(stk,type) (*((type*)(  (stk)->pos += sizeof(type) )  - 1 )  )
#define stackSkip(stk) ((stk)->pos += sizeof(int) )

#else

#define push(stk,x,type) ((stk)->ptr < (stk)->safe ? \
                           ( *(type *)((stk)->ptr) = (x) , (stk)->ptr += STACK_ALIGNMENT) : \
			    (stackExtend (stk,16), \
			     *(type *)((stk)->ptr) = (x) , (stk)->ptr += STACK_ALIGNMENT) )
#define pop(stk,type)    (  ((stk)->ptr -= STACK_ALIGNMENT) >= (stk)->a->base ? \
			    *((type*)((stk)->ptr)) : \
                          (messcrash ("User stack underflow"), *((type*)1)) )
#define stackNext(stk,type) (*((type*)(  ((stk)->pos += STACK_ALIGNMENT ) - \
                                             STACK_ALIGNMENT ))  )
#define stackSkip(stk) ((stk)->pos += STACK_ALIGNMENT )
#endif

#if STACK_DOUBLE_ALIGNMENT > STACK_ALIGNMENT
void ustackDoublePush(Stack stk, double x);
double ustackDoublePop(Stack stk);
double stackDoubleNext(Stack stk);
#define pushDouble(stk,x) ustackDoublePush(stk, x)
#define popDouble(stk) ustackDoublePop(stk)
#define stackDoubleNext(stk) ustackDoubleNext(stk) 
#else
#define pushDouble(stk,x) push(stk, x, double)
#define popDouble(stk) pop(stk, double)
#define stackDoubleNext(stk) stackNext(stk, double)
#endif  


int     pushText (Stack s, const char *text) ; /* returns n, i.e. stackText (s, n) recovers the pushed text */
char*   popText (Stack s) ;	/* returns last text and moves pointer before it */
void    catText (Stack s, const char *text) ;  /* like strcat */
void	stackTokeniseTextOn(Stack s, char *text, char *delimiters) ; /* tokeniser */

int     stackMark (Stack s) ;              /* returns a mark of current ptr */
int     stackPos (Stack s) ;              /* returns a mark of current pos, useful with stackNextText */
void    stackCursor (Stack s, int mark) ;  /* sets ->pos to mark */
#define stackAtEnd(stk)     ((stk)->pos >= (stk)->ptr)
char*   stackNextText (Stack s) ;
 
#define stackText(stk,mark) ((char*)((stk)->a->base + (mark)))
#define stackTextForceFeed(stk,j) (arrayForceFeed((stk)->a,j) ,\
                 (stk)->ptr = (stk)->pos = (stk)->a->base + (j) ,\
                 (stk)->safe = (stk)->a->base + (stk)->a->dim - 16 )
void    catBinary (Stack s, const void* data, int size) ;
 
/********** Line breaking package **********/
 
int     uLinesText (char *text, int width) ;
char    *uNextLine (char *text) ;
char    *uPopLine (char *text) ;
char    **uBrokenLines (char *text, int width) ; /* array of lines */
char    *uBrokenText (char *text, int width) ; /* \n's intercalated */

/********** Associator package *************/

/*
*
* The associator package is a keyed lookup table.  The keys are
* void pointers - the address is the key, not anything at the address. 
* The values stored in the table are also void pointers, presumably 
* the address of some interesting thing.
*
* Keys may not have a value of ((void *)-1) or ((void *)0)
*
* You add (key,value) pairs to the associator, then you can look
* them up later.  If you choose, you can allow more than one value
* with the same key.
*
* Associator *a = assCreate();
* assDestroy(a);
* 	create/destroy associator.
*
* BOOL uAssFind (Associator a, void* key, void **value) ;
*	key is the key to look up
*	*value is set to the value found
*	returns true if key was found, false otherwise
*
* assFind() is the same as uAssFind() with a cast for the value parameter
*
* assFindNext is thread safe in read mode, but should not be intermingled with assMultipleInsert
* usage:
*   Array bucket = 0 ;
*   int iBucket = 0 ;
*   void *out ;
*   while (assFindnext (ass, assVoid(in), &out, &bucket, &iBucket))
*        { do something with *out ; }
*
* BOOL assInsert (Associator a, void* key, void* value) ;
*	If the key is already in use in the associator, the 
*	(key,value) are not inserted and the function returns FALSE.
*	Otherwise the (key,value) are inserted and the function returns
*	TRUE.
*
* void assMultipleInsert(Associator a, void* key, void* value);
*	Adds the (key,value) pair to the associator.  You can have
*	arbitrarily many values associated with the same key.
* 	This call cannot fail.
* 
* BOOL assRemove (Associator a, void* key) ;
*	Removes the first (key,value) pair associated with the given key.
*	
* BOOL uAssNext (Associator a, void**key, void**value, unsigned long int *hashPosp);
*	Iterates through entire associator table:
*		void *key, *value;
*               unsigned long int hashPos = 0 ;
*		key = 0;
*		while (uAssNext(a,&key,&value, &hashPos))
*			printf("key %x value %x\n",key,value)
*
* assVoid(i)
*	turn an integer into a void * so you can use it as a key or value
*
* assInt(v)
*	turn a void * back into an integer
*	
*
* Associator is implemented as a conventional hash table with
* collision resolution, not buckets.
* It is combined with a Bloom which accelerates negative searches.
*
*/
#include "bitset_.h"

typedef struct AssStruct
  { int magic ;                 /* Ass_MAGIC */
    int id ;                    /* unique identifier */
    long int n ;			/* number of items stored */
    int nbits ;			/* size = 2^nbits */
    long int i ; int iBucket  ;   /* non reentrant counters used by assNext */
    int floatingLine;       /* default 1= fill to 1/2; 2->to 1/4 etc */
    const void **in,**out ;  /* **out leftmost bit indicates an entry in aa */
    unsigned long mask ;		/* 2^(nbits+floatingLine) - 1 */
    BitSet bbh, bbd, bb3 ;       /* bitsets marking probable existence */
    BOOL inOnly ;         /* if TRUE, **out is disabled */
    BOOL readOnly ;       /* while TRUE, all modifs, except destroy, are forbidden */
    AC_HANDLE h ;
    BigArray aa ;     /* array of arrays of values, allocated on h, soring  multiple **out values */
  } *Associator ; 
 
#define assExists(a) ((a) && (a)->magic == ASS_MAGIC ? (a)->id : 0 )

#if !defined(MEM_DEBUG)
  Associator assHandleCreate (AC_HANDLE handle) ;
  Associator assBigHandleCreate (long int size, AC_HANDLE handle) ;
  Associator assInOnlyHandleCreate (long int size, AC_HANDLE handle) ;
#define assBigCreate(s) assBigHandleCreate(s,0)
#define assCreate() assHandleCreate(0)
#else
  Associator assHandleCreate_dbg (AC_HANDLE handle,
				 const char *hfname, int hlineno) ;
  Associator assBigCreate_dbg (long int size, const char *hfname, int hlineno) ;
  Associator assBigHandleCreate_dbg (long int size, AC_HANDLE handle, const char *hfname, int hlineno) ;
  Associator assInOnlyHandleCreate_dbg (long int size, AC_HANDLE handle, const char *hfname, int hlineno) ;
#define assInOnlyHandleCreate(s,h) assInOnlyHandleCreate_dbg(s,h,__FILE__, __LINE__)
#define assCreate() assHandleCreate_dbg(0, __FILE__, __LINE__)
#define assHandleCreate(h) assHandleCreate_dbg(h, __FILE__, __LINE__)
#define assBigCreate(s) assBigHandleCreate_dbg(s,0, __FILE__, __LINE__)
#define assBigHandleCreate(s,h) assBigHandleCreate_dbg(s,h, __FILE__, __LINE__)
#endif

Associator assReCreate (Associator a) ;
void    uAssDestroy (Associator a) ;
BOOL assSetReadOnly (Associator a, BOOL b) ; /* forbids assInsert allowing multithreading */
#define assDestroy(x)  ((x) ? uAssDestroy(x), x = 0, TRUE : FALSE)

BOOL    uAssFind (Associator a, const void* xin, const void** pout, unsigned long int *hashPosp) ; 
#define assFind(ax,xin,pout)    uAssFind((ax),(xin),(const void**)(pout), 0)
          /* if found, updates *pout and returns TRUE, else returns FALSE */
BOOL assFindNext (Associator a, const void* xin, const void** pout, Array *bucketp, int *iBucketp) ;

BOOL    assInsert (Associator a, const void* xin, const void* xout) ;
            /* if already there returns FALSE, else inserts and returns TRUE */
BOOL    assMultipleInsert(Associator a, const void* xin, const void* xout) ;
           /* allow multiple Insertions */
BOOL    assRemove (Associator a, const void* xin) ;
            /* if found, removes entry and returns TRUE, else returns FALSE */
void    assDump (Associator a) ;
           /* for debug - uses printf */
void    assClear (Associator a) ;
BOOL    uAssNext (Associator a, const void* *pin, const void* *pout) ;
BOOL    uNewAssNext (Associator a, const void* *pin, const void* *pout, unsigned long int *hashPosp) ;
#define assNext(_ax,_pin,_pout)	uAssNext((_ax),(const void**)(_pin),(const void**)(_pout))
#define newAssNext(_ax,_pin,_pout,_hp)	uNewAssNext((_ax),(const void**)(_pin),(const void**)(_pout),_hp)
/* convert an integer to a void * without generating a compiler warning */
#define assVoid(i) ((void *)(((char *)0) + (i))) 
#define assInt(v) ((int)(0xffffffff & (((char *)v) - ((char *)0))))
#define assULong(v) ((unsigned long)(((char *)v) - ((char *)0)))

#endif /* defined(DEF_ARRAY_H) */

/**************************** End of File ******************************/

 
