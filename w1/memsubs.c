/*  File: memsubs.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  4 11:24 1998 (fw)
 * Created: Thu Aug 20 16:54:55 1998 (rd)
 *-------------------------------------------------------------------
 */
/* define MALLOC_CHECK here to check mallocs - or preferably in regular.h */

#include "regular.h"

#if defined(NEXT) || defined(HP)
  extern void* malloc (mysize_t size) ;
#elif defined(MAC_X)
  extern void* malloc (size_t size) ;
#elif !defined(WIN32) 
#include <malloc.h>   /* normal machines  */
#endif

static void* myMalloc (mysize_t size, BOOL recursion) ;
static void myFree (void *vp) ;

/********** primary type definition **************/

typedef struct _AC_HANDLE_STRUCT {
  AC_HANDLE next ;	/* for chaining together on Handles */
  AC_HANDLE back ;	/* to unchain */
  void (*final)(void*) ;	/* finalisation routine */
  mysize_t size ;			/* of user memory to follow */
#ifdef MALLOC_CHECK
  int check1 ;			/* set to known value */
#endif
} AC_HANDLE_STRUCT ;

/*********************************************************************/
/********** memory allocation - messalloc() and handles  *************/

static long nMessAlloc = 0 ;
static long numMessAlloc = 0 ;
static long totMessAlloc = 0 ;
static long maxMessAlloc = 0 ;
static long totMessServed = 0 ;


  /* Calculate to size of an AC_HANDLE_STRUCT rounded to the nearest upward
     multiple of sizeof(double). This avoids alignment problems when
     we put an AC_HANDLE_STRUCT at the start of a memory block

 */
#if ! defined(GPLUSPLUS)
/* no idea why this is refused by g++, march 22, 2005 */
#define STORE_OFFSET ((((sizeof(AC_HANDLE_STRUCT)-1)/MALLOC_ALIGNMENT)+1)\
                             * MALLOC_ALIGNMENT)
#else
#define STORE_OFFSET (6*sizeof(void*))
#endif

  /* macros to convert between a void* and the corresponding AC_HANDLE */

#define toAllocUnit(x) (AC_HANDLE) ((char*)(x) - STORE_OFFSET)
#define toMemPtr(unit)((void*)((char*)(unit) + STORE_OFFSET))

#ifdef MALLOC_CHECK
BOOL handlesInitialised = FALSE;
static Array handles = 0 ;

  /* macro to give the terminal check int for an AC_HANDLE_STRUCT */
  /* unit->size must be a multiple of sizeof(int) */
#define check2(unit)  *(int*)(((char*)toMemPtr(unit)) + (unit)->size)

static void checkUnit (AC_HANDLE unit) ;
static int handleOrder (const void *a, const void  *b)
{ return (*((const AC_HANDLE *)a) == *((const AC_HANDLE *)b)) ? 0 : 
	  (*((const AC_HANDLE *)a) > *((const AC_HANDLE *)b)) ? 1 : -1 ;
}
#endif
#if defined(MALLOC_CHECK) || defined(MEM_DEBUG)
AC_HANDLE_STRUCT handle0 ;
#endif

/************** halloc(): key function - messalloc() calls this ****************/
static BOOL isThreadSafe = FALSE ;
void memSetIsMultiThreaded (void )
{
  isThreadSafe = TRUE ;
}

#ifdef MEM_DEBUG
void *halloc_dbg(mysize_t size, AC_HANDLE handle,const char *hfname, int hlineno) 
#else
void *halloc(mysize_t size, AC_HANDLE handle)
#endif
{ 
  AC_HANDLE unit ;
  
#ifdef MALLOC_CHECK
  if (!handlesInitialised)		/* initialise */
    { handlesInitialised = TRUE;
      /* BEWARE, arrayCreate calls handleAlloc, line above must precede
         following line to avoid infinite recursion */
      handles = arrayCreate (16, AC_HANDLE) ;
      array (handles, 0, AC_HANDLE) = &handle0 ;
      handle0.next = 0 ;
    }

  while (size % INT_ALIGNMENT) size++ ; /* so check2 alignment is OK */
  unit = (AC_HANDLE) myMalloc (STORE_OFFSET + size + sizeof(int), FALSE) ;
#else
  unit = (AC_HANDLE) myMalloc (STORE_OFFSET + size, FALSE) ;
#endif

  if (!unit)			/* out of memory -> messcrash */
    {
      messcrash (
 "Memory allocation failure when requesting %u bytes, %ld Mb already allocated (max %ld Mb)", 
  size, totMessAlloc/(1024*1024), maxMessAlloc/(1024*1024)) ;
    }
#if defined(MALLOC_CHECK) || defined(MEM_DEBUG)
  if (! isThreadSafe && !handle)
    handle = &handle0 ;
#endif
  if (handle) 
    { unit->next = handle->next ;
      unit->back = handle ;
      if (handle->next) (handle->next)->back = unit ;
      handle->next = unit ;
    }

  unit->size = size ;
#ifdef MALLOC_CHECK
  unit->check1 = 0x12345678 ;
  check2(unit) = 0x12345678 ;
  if (0 && totMessAlloc > 100000000 && size > 20000000)
    invokeDebugger() ;

#endif

  ++numMessAlloc ;
  ++nMessAlloc ;
  totMessAlloc += size ;
  totMessServed += size ;
  if (totMessAlloc > maxMessAlloc)
    maxMessAlloc = totMessAlloc ;
  return toMemPtr(unit) ;
}

void blockSetFinalise(void *block, void (*final)(void *))
{ AC_HANDLE unit = toAllocUnit(block);
  unit->final = final ;
}  

/***** handleAlloc() - does halloc() + blockSetFinalise() - archaic *****/

#ifdef MEM_DEBUG
void *handleAlloc_dbg (void (*final)(void*), AC_HANDLE handle, mysize_t size,
		   const char *hfname, int hlineno)
{
  void *result = halloc_dbg(size, handle, hfname, hlineno) ;
#else
void *handleAlloc (void (*final)(void*), AC_HANDLE handle, mysize_t size)
{
  void *result = halloc(size, handle);
#endif
  if (final) 
    blockSetFinalise(result, final);

  return result;
}

/****************** useful utility ************/
/* ensure 2 terminal zeroes */
#ifdef MEM_DEBUG
char *strnew_dbg(const char *old, AC_HANDLE handle, const char *hfname, int hlineno)
{ char *result = 0 ;
  if (old)
    { result = (char *)halloc_dbg(2+strlen(old), handle, hfname, hlineno) ;
#else
char *strnew(const char *old, AC_HANDLE handle)
{ 
  char *result = 0 ;
  if (old)    /* 06_05_18 changed 2->PTR_ALIGNMENT to try to please purify 
	      *    please change back to 2 as soon as possible
	      */
    { 
      result = (char *)halloc(PTR_ALIGNMENT+strlen(old), handle);
#endif
      strcpy(result, old);
    }
  return result;
}

/****************** messfree ***************/

void umessfree (void *cp)
{
  AC_HANDLE unit = toAllocUnit(cp) ;

#ifdef MALLOC_CHECK
  checkUnit (unit) ;
  unit->check1 = 0x87654321; /* test for double free */
#endif

  if (unit->final)
    (*unit->final)(cp) ;

  if (unit->back) 
    { (unit->back)->next = unit->next;
      if (unit->next) (unit->next)->back = unit->back;
    }
  
  --numMessAlloc ;
  totMessAlloc -= unit->size ;
#ifdef MALLOC_CHECK
  memset(cp, 0xaa, unit->size);
  memset(unit, 0xaa, sizeof(struct _AC_HANDLE_STRUCT));
#endif
  myFree (unit) ;
}

/************** create and destroy handles **************/

/* NOTE: handleDestroy is #defined in regular.h to be messfree */
/* The actual work is done by handleFinalise, which is the finalisation */
/* routine attached to all AC_HANDLEs. This allows multiple levels */
/* of free-ing by allocating new AC_HANDLES on old ones, using */
/* handleHandleCreate. handleCreate is simply defined as handleHandleCreate(0) */

static void handleFinalise (void *p)
{
  AC_HANDLE handle = (AC_HANDLE)p;
  AC_HANDLE next, unit = handle->next ;

/* do handle finalisation first */  
  if (handle->final)
    (*handle->final)((void *)handle->back);

      while (unit)
    { 
#ifdef MALLOC_CHECK
      checkUnit (unit) ;
      unit->check1 = 0x87654321; /* test for double free */
#endif
      if (unit->final)
	(*unit->final)(toMemPtr(unit)) ;
      next = unit->next ;
      --numMessAlloc ;
      totMessAlloc -= unit->size ;
      myFree (unit) ;
      unit = next ;
    }

#ifdef MALLOC_CHECK
  arrayRemove (handles, &p, handleOrder) ;
#endif

/* This is a finalisation routine, the actual store is freed in messfree,
   or another invokation of itself. */
}
  
void handleSetFinalise(AC_HANDLE handle, void (*final)(void *), void *arg)
{ handle->final = final;
  handle->back = (AC_HANDLE)arg;
}

AC_HANDLE handleHandleCreate(AC_HANDLE handle)
{ 
  AC_HANDLE res = (AC_HANDLE) handleAlloc(handleFinalise, 
						handle,
						sizeof(AC_HANDLE_STRUCT));
#ifdef MALLOC_CHECK
  /* NB call to handleAlloc above ensures that handles is initialised here */
  arrayInsert (handles, &res, handleOrder) ;
#endif
  res->next = res->back = 0 ; /* No blocks on this handle yet. */
  res->final = 0 ; /* No handle finalisation */
  return res ;
}

BOOL finalCleanup = FALSE ;
#ifdef MEM_DEBUG
void handleCleanUp (void) 
{ finalCleanup = TRUE ;
  handleFinalise ((void *)&handle0) ;
}
#endif

AC_HANDLE ac_new_handle (void)
{ 
  /*
   * creates a new AC_HANDLE, 
   * freeing the handle will free recursivelly all memory allocated on the handle
   * The handle library was written for acedb by Simon Kelley
   */
  return (AC_HANDLE) handleCreate () ;
}

AC_HANDLE ac_handle_handle (AC_HANDLE parent_handle)
{    
  /*
   * creates a new AC_HANDLE on a parent handle, 
   * freeing the parent_handle will free the
   * new handle and recursivelly all memory allocated on it 
   * The handle library was written for acedb by Simon Kelley
   */
  return (AC_HANDLE) handleHandleCreate ((AC_HANDLE) parent_handle) ;
}

/* Reentrant utility function, export time in a buffer allocated on the handle */
extern char *timeBufShowNow (char *timeBuf) ;
char *timeHandleShowNow (AC_HANDLE handle)
{
  char *timeBuf = halloc (25, handle) ;
  return timeBufShowNow (timeBuf) ;
}


/* Reentrant utility function, do printf into dynamically allocated memory,
   allocate the memory on a handle, if provided.
   Consider this to replace messprintf everywhere....
*/
char *hprintf (AC_HANDLE h, char *format, ...)
{
  /* ACFORMAT ("")  is not reentrant, we need a version allocating on h */
  int len ;
  char *message ; 
  va_list ap; 

  va_start(ap, format);
  len = utPrintfSizeOfArgList (format, ap) ;
  va_end (ap) ;
 
  message = halloc (len + 1, h) ;

  va_start(ap, format);
  vsprintf(message, format, ap) ;
  va_end (ap) ;

  return message ;
}

/************** checking functions, require MALLOC_CHECK *****/

#ifdef MALLOC_CHECK
static void checkUnit (AC_HANDLE unit)
{
  if (unit->check1 == 0x87654321)
    messerror ("Block at %x freed twice - bad things will happen.",
	       toMemPtr(unit));
  else
    if (unit->check1 != 0x12345678)
      messerror ("Malloc error at %x length %d: "
		 "start overwritten with %x",
		 toMemPtr(unit), unit->size, unit->check1) ;
  
  if (check2(unit) != 0x12345678)
    messerror ("Malloc error at %x length %d: "
	       "end overwritten with %x",
	       toMemPtr(unit), unit->size, check2(unit)) ;
}

void messalloccheck (void)
{
  int i ;
  AC_HANDLE unit ;

  if (!handles) return ;

  for (i = 0 ; i < arrayMax(handles) ; ++i) 
    for (unit = arr(handles,i,AC_HANDLE)->next ; unit ; unit=unit->next)
      checkUnit (unit) ;
}
#else
void messalloccheck (void) {}
#endif

/******************* status monitoring functions ******************/

void handleInfo (AC_HANDLE handle, int *number, mysize_t *size)
{
  AC_HANDLE unit = handle->next;

  *number = 0;
  *size = 0;

  while (unit)
    { ++*number ;
      *size += unit->size ;
      unit = unit->next ;
    }
}

int messAllocStatus2 (int *ser)
{ 
  *ser = (int) (totMessServed/(1024*1024)) ;
  return (int) numMessAlloc ;
}

int messAllocStatus (int *mem)
{ 
  *mem = (int) (totMessAlloc/(1024*1024)) ;
  return (int) numMessAlloc ;
}

int messAllocMaxStatus (int *mem)
{ 
  *mem = (int) (maxMessAlloc/(1024*1024)) ;
  return (int) numMessAlloc ;
}

/****************************************************************/
/****************************************************************/
/**********  New malloc packge , mieg, dec 2000 *****************/
/****************************************************************/
/* 
mieg: dec 2000
i tried in this short package to better the memory usage
but it does not seem very efficient

more testing needed
originaly i thiught i had a pseudo memory leak on solaris
when malloc/free random size buffers, but more recently
i could not reproduce the problem
*/


#define ggZZNEWMALLOC
#if defined(ZZNEWMALLOC)

attention, this draft stupid packeged is not thread safe

#define nwAllocMax 1204 
static BOOL zAlloc = 0 ;
static void** vAlloc[256] ; /* address of the free buffers */
static void* wAlloc[nwAllocMax] ; /* adress of first Gb, to please purify */
static int nAlloc[256] ;   /* number of available buffers */
static int sAlloc[256] ;   /* size of vAlloc list */
static int tAlloc = 0, uAlloc = 0, nwAlloc = 0 ; /* global stats */

static void myMallocInit (void)
{
  int i = 256 ;
  while (i--)
    { vAlloc[i] = 0 ; nAlloc[i] = 0 ; sAlloc[i] = 0 ;}
  zAlloc = TRUE ;
  tAlloc = uAlloc = 0 ;
}

static void* myMalloc (mysize_t size, BOOL recursion) 
{
  mysize_t p = 5, pp, ppw ;
  char *cp ;

  if (size <= 0 || size > (1 << 30))
    messcrash ("Excessive call to mymalloc, requesting %d bytes", size) ;

  pp = (1 << p) ;
  while (pp < size) { p++ ; pp <<= 1 ; }
  ppw = pp + MALLOC_ALIGNMENT ;

  if (!zAlloc)
    myMallocInit () ;

  if (recursion)
    {
      cp = (char*)  malloc (ppw) ;
      if (!cp)
	messcrash ("cannot allocate %d bytes, %d/%d allready used", tAlloc, uAlloc) ;
      tAlloc += pp + MALLOC_ALIGNMENT ; 
      if (nwAlloc < nwAllocMax)
	wAlloc[nwAlloc++] = cp ;
    }
  else
    {
      if (!sAlloc[p])
	{
	  int r = 10, rr, rrw ;

	  /* prepare a buffer that can later be recovered in myfree */
	  rr = (1 << r) ; rrw = rr + MALLOC_ALIGNMENT ;

	  cp = (char*)  malloc (rrw) ;
	  if (!cp)
	    messcrash ("cannot allocate %d bytes, %d/%d allready used", tAlloc, uAlloc) ;
	  tAlloc += rrw ;
	  if (nwAlloc < nwAllocMax)
	    wAlloc[nwAlloc++] = cp ;

	  cp[0] = r ; cp[1] = 247 ; /* magic */ 
	  cp += MALLOC_ALIGNMENT ;
	  vAlloc[p] = (void**) (cp) ;
	  sAlloc[p] = rr/(sizeof(void*)) ;
	  nAlloc[p] = 0 ;
	}
      
      if (!nAlloc[p])
	{
	  int j, pp = 1 << p , ppw = pp + MALLOC_ALIGNMENT ;
	  
	  for (j = 1 ; j < sAlloc[p] && j * pp < (1 << 12) ; j++) ; /* at least 4M of anything */

	  cp = (char *) malloc (j * ppw) ;
	  if (!cp)
	    messcrash ("cannot allocate %d bytes, %d/%d allready used", tAlloc, uAlloc) ;
	  tAlloc += j * ppw ;
	  if (nwAlloc < nwAllocMax)
	    wAlloc[nwAlloc++] = cp ;

	  while (j--)
	    { vAlloc[p][nAlloc[p]++] = (void*) cp ; *cp = 0 ; cp += ppw ; }
	}
      
      cp = vAlloc[p][--nAlloc[p]] ;
    }
  /* printf ("p=%d nAlloc[p]=%d sAlloc[p]=%d\n", p,  nAlloc[p], sAlloc[p]) ; */
  memset (cp, 0, ppw) ;
  cp[0] = p ;
  cp[1] =  247 ; /* magic */ 
  cp += MALLOC_ALIGNMENT ;

  uAlloc += pp ;
  return (void*) cp ;
}

static void myFree (void *vp)
{
   unsigned char *cp = (unsigned char *) vp - MALLOC_ALIGNMENT ;
   int p = cp[0] ;
   int q = cp[1] ;

   if (q != 247)  /* magic */
     messcrash ("Bad %d magic in myfree, you probably overwrote the memory in you application try to use PURIFY", q) ;
   if (p > 255 || p < 0 || !sAlloc[p])
     messcrash ("Non allocated p=%d  in myfree,  you probably overwrote the memory in you application try to use PURIFY", q) ;

   if (nAlloc[p] >= sAlloc[p])
     {
       void *zp, *wp = myMalloc (2 * sAlloc[p] * (sizeof(void*)), TRUE) ; /* recursive call ! */
       memcpy (wp, vAlloc[p], sAlloc[p] * (sizeof(void*))) ;
       zp = vAlloc[p] ;
       vAlloc[p] = wp ;
       sAlloc[p] *= 2 ;
       myFree (vAlloc[p]) ; /* now safe i think becase vAlloc is already doubled */
     }
   vAlloc[p][nAlloc[p]++] = cp ;
   return ;
}

void messAllocShutdown (void) 
{
  while(nwAlloc--)
    free (wAlloc[nwAlloc]) ;
}

#else

#if 0

FILE * malloc_fp = NULL;

static void* myMalloc (mysize_t size, BOOL recursion) 
{
  char *cp = malloc (size) ;
  if (cp)  memset (cp, 0, size) ;
	if (! malloc_fp )
		{
		malloc_fp = fopen("malloc_log", "w");
		if (! malloc_fp )
			messcrash("cannot create malloc_log for debugging");
		}
	putc( 'm', malloc_fp );
	putc( ((int)cp) , malloc_fp );
	putc( ((int)cp) >> 8 , malloc_fp );
	putc( ((int)cp) >> 16 , malloc_fp );
	putc( ((int)cp) >> 24 , malloc_fp );

	putc( size , malloc_fp );
	putc( size >> 8 , malloc_fp );
	putc( size >> 16 , malloc_fp );
	putc( size >> 24 , malloc_fp );

  return cp ;
}

static void myFree (void *cp) 
{
	putc( 'f', malloc_fp );
	putc( ((int)cp) , malloc_fp );
	putc( ((int)cp) >> 8 , malloc_fp );
	putc( ((int)cp) >> 16 , malloc_fp );
	putc( ((int)cp) >> 24 , malloc_fp );
	free (cp) ;
}

void messAllocShutdown (void) { fclose (malloc_fp) ;return ; }

#else

static void *myMalloc (mysize_t size, BOOL recursion) 
{ 
  int ii ;
  void *cp = 0 ;

  /* ii >= 0 look silly, but on 2011_02_12, we looped with ii == -266904305 */
  for (ii = 0 ; ! cp && ii < 30 ; ii++)
    {
      cp = malloc (size) ;
      if (cp)
	{ memset (cp, 0, size) ;  }
      else  /* wishful thinking for 5 minutes */
	{
	  invokeDebugger () ;
	  fprintf (stderr, "myMalloc failed allocating %ld bytes allready alloc %ld, waiting loop %d\n", (long)size, totMessAlloc, ii) ;
	  sleep (10) ;
	}
    }
  return cp ;
}
static void myFree (void *vp) { free (vp) ;}
void messAllocShutdown (void) { return ; }

#endif

#endif /* NEWMALLOC */


/*************************** end of file ************************/
/****************************************************************/


