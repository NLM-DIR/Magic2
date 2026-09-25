/*
 * scratch.c

 * This module is part of the acedb package
 * Authors: Jean Thierry-Mieg, NCBI/NLM/NIH
 * Created Sep 24 2026

 * This is public.

 * This module implements temporary array storage on disk
 
 * Usage:
     SCR *s = scratchCreate (dirName, h1) ;
     Array aa ;
     int nn = scratchArrayPut (scr, aa) ;    : put the array in scr 
Array bb = scratchArrayGet (scr, nn, h2) ;   : recover a new array allocated on h2 
     ac_free (s) ;    or ac_free (h1) ;


SCR scr = scratchCreate (h, "/scratch/myjob") ;
int k = scratchArrayPut (scr, aa) ; arrayMax (aa) = 0 ;


Array bb = scratchArrayGet (scr, k, 0, MYTYPE, h) ;
MYTYPE *up = arrp (bb, 0, MYTYPE) ;

 or reuse a buffer across many reads (better in the merge loop: no realloc) 
Array buf = arrayHandleCreate (100000, MYTYPE, h) ;
for (k = 0 ; k < scratchArrayCount (scr) ; k++)
  {
    scratchArrayGet (scr, k, buf, MYTYPE, 0) ;
    for (i = 0, up = arrp (buf, 0, MYTYPE) ; i < arrayMax (buf) ; i++, up++)
      ...
  }


*/
#define _GNU_SOURCE        /* O_TMPFILE, mkostemp: before ALL #include */
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include "ac.h"
#include "scratch.h"
#define SCRATCHMAGIC 3154992

typedef struct scrBlockStruct { off_t pos ; int aMax ; int aSize ; } SCR_ ;
struct scrStruct { int magic ; int fd ; off_t end ; Array blocks ; } ;
static void scratchCheck (SCR scr, const char *caller) ;

int scratchArrayCount (SCR scr)
{
  scratchCheck (scr, "scratchArrayCount") ;
  return arrayMax (scr->blocks) ;
} /* scratchArrayCount */

/**************************************************************/

static int fdReadWrite (int fd, char *p, size_t n, off_t off, int isWrite)
{
  while (n > 0)
    {
      ssize_t r = isWrite ? pwrite (fd, p, n, off) : pread (fd, p, n, off) ;
      if (r < 0) { if (errno == EINTR) continue ; return -1 ; }
      if (r == 0) { errno = EIO ; return -1 ; }
      p += r ; off += r ; n -= r ;
    }
  return 0 ;
} /* fdReadWrite */

/**************************************************************/

static void scratchCheck (SCR scr, const char *caller)
{
  if (! scr || scr->magic != SCRATCHMAGIC)
    messcrash ("%s called on a bad or freed SCR", caller) ;
} /* scratchCheck */

/**************************************************************/

static void scratchFinalize (void *vp)
{
  SCR scr = (SCR) vp ;
  scratchCheck (scr, "scratchFinalize") ;
  scr->magic = 0 ;
  if (scr->fd >= 0) close (scr->fd) ;   /* disk space released here */
  arrayDestroy (scr->blocks) ;  
} /* scratchFinalize */

/**************************************************************/

SCR scratchCreate (const char *dirName, AC_HANDLE h)
{
  SCR scr = (SCR) handleAlloc (scratchFinalize, h, sizeof (struct scrStruct)) ;
  int fd ;

  scr->magic = SCRATCHMAGIC ; scr->fd = -1 ; scr->end = 0 ;
  scr->blocks = arrayCreate (128, SCR_) ;
  if (! dirName || ! *dirName) dirName = "." ;

  fd = open (dirName, O_TMPFILE | O_RDWR | O_CLOEXEC, 0600) ;  /* no O_CREAT, no name */
  if (fd < 0)            /* filesystem without O_TMPFILE: named file, then unlink */
    {
      char buf[4096] ;
      snprintf (buf, sizeof (buf), "%s/scratch.XXXXXX", dirName) ;
      fd = mkostemp (buf, O_CLOEXEC) ;
      if (fd < 0)
        messcrash ("scratchCreate cannot create a file in %s: %s", dirName, strerror (errno)) ;
      unlink (buf) ;
    }
  scr->fd = fd ;
  return scr ;
} /* scratchCreate */

/**************************************************************/
/* appends aa, returns its block number 0,1,2..., aa is not modified */
int scratchArrayPut (SCR scr, Array aa)
{
  int k ;
  SCR_ *s ;
  size_t n ;

  scratchCheck (scr, "scratchArrayPut") ;
  if (! arrayExists (aa))
    messcrash ("scratchArrayPut called on a bad Array") ;

  k = arrayMax (scr->blocks) ;
  s = arrayp (scr->blocks, k, SCR_) ;
  s->pos = scr->end ;
  s->aMax = arrayMax (aa) ;
  s->aSize = aa->size ;  /* set even if empty */

  n = (size_t) s->aMax * s->aSize ;
  if (n && fdReadWrite (scr->fd, aa->base, n, scr->end, 1))
    messcrash ("scratchArrayPut failed: %s", strerror (errno)) ;
  scr->end += n ;
  return k ;
} /* scratchArrayPut */

/**************************************************************/
/* Reuses aa or create a new Array on handle h
 */
Array scratchArrayDoGet (SCR scr, int k, Array aa, int typeSize, AC_HANDLE h)
{
  SCR_ *s ;

  scratchCheck (scr, "scratchArrayGet") ;
  if (k < 0 || k >= arrayMax (scr->blocks))
    messcrash ("scratchArrayGet: block %d not in [0, %d[", k, arrayMax (scr->blocks)) ;
  s = arrp (scr->blocks, k, SCR_) ;
  if (typeSize != s->aSize)
    messcrash ("scratchArrayGet: block %d has records size %d, not %d"
               , k, s->aSize, typeSize) ;

  if (! aa)
    aa = uArrayCreate (s->aMax, s->aSize, h) ;  
  else if (aa->size != s->aSize)
    messcrash ("scratchArrayGet: block %d has records of size %d, the Array has size %d"
               , k, s->aSize, aa->size) ;

  if (s->aMax > 0)
    {
      uArray (aa, s->aMax - 1) ;                 /* = arrayp (aa, iMax-1, TYPE): make room */
      if (fdReadWrite (scr->fd, aa->base, (size_t) s->aMax * s->aSize, s->pos, 0))
        messcrash ("scratchArrayGet block %d failed: %s", k, strerror (errno)) ;
    }
  arrayMax (aa) = s->aMax ;
  return aa ;
} /* scratchArrayDoGet */

/***************************/

