/*
 * scratch.c

 * This module is part of the acedb package
 * Authors: Jean Thierry-Mieg, NCBI/NLM/NIH
 * Created Sep 24 2026

 * This code is public.

 * This module implements transparent temp disk storage of BigArrays
 * hence diminishing the live RAM used by a long program
 * Usage:
 * The stored arrays are numbered sequencially
 * They must be of homogeneous TYPE (inforced weakly via aa->size)
 * The disk space is released by calling ac_free (scr) or ac_free (h)
 *
 * int k = scratchPut (scr, aa) : returns a chronological integer
 * int kMax = scratchCount (scr) : return the number of Put
 * scratchGet (scr, aa, k, TYPE, h)
 *    k > kMax is forbidden
 *    if (aa == 0) Get allocates a new array of TYPE on h 
 *    if (aa != 0) Get checks aa->size, reuses aa and adjust aa->max

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

int scratchCount (SCR scr)
{
  scratchCheck (scr, "scratchArrayCount") ;
  return arrayMax (scr->blocks) ;
} /* scratchCount */

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
  SCR scr = (SCR ) vp ;
  scratchCheck (scr, "scratchFinalize") ;
  scr->magic = 0 ;
  if (scr->fd >= 0) close (scr->fd) ;   /* disk space released here */
  arrayDestroy (scr->blocks) ;  
} /* scratchFinalize */

/**************************************************************/

SCR scratchCreate (const char *dirName, AC_HANDLE h)
{
  SCR scr = (SCR ) handleAlloc (scratchFinalize, h, sizeof (struct scrStruct)) ;
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
int scratchPut (SCR scr, BigArray aa)
{
  int k ;
  SCR_ *s ;
  size_t n ;

  scratchCheck (scr, "scratchBigArrayPut") ;
  if (! arrayExists (aa))
    messcrash ("scratchBigArrayPut called on a bad BigArray") ;

  k = arrayMax (scr->blocks) ;
  s = arrayp (scr->blocks, k, SCR_) ;
  s->pos = scr->end ;
  s->aMax = arrayMax (aa) ;
  if (k == 0) s->aSize = aa->size ;  /* set even if aa is empty */

  n = (size_t) s->aMax * s->aSize ;
  if (n && fdReadWrite (scr->fd, aa->base, n, scr->end, 1))
    messcrash ("scratchBigArrayPut failed: %s", strerror (errno)) ;
  scr->end += n ;
  return k ;
} /* scratchPut */

/**************************************************************/
/* Reuses aa or create a new BigArray on handle h
 */
BigArray uScratchGet (SCR scr, int k, BigArray aa, int typeSize, AC_HANDLE h)
{
  SCR_ *s ;

  scratchCheck (scr, "scratchBigArrayGet") ;
  if (k < 0 || k >= arrayMax (scr->blocks))
    messcrash ("scratchBigArrayGet: block %d not in [0, %d[", k, arrayMax (scr->blocks)) ;
  s = arrp (scr->blocks, k, SCR_) ;
  if (typeSize != s->aSize)
    messcrash ("scratchBigArrayGet: block %d has records size %d, not %d"
               , k, s->aSize, typeSize) ;

  if (! aa)
    aa = uBigArrayCreate (s->aMax, s->aSize, h, TRUE) ;  
  else if (aa->size != s->aSize)
    messcrash ("scratchBigArrayGet: block %d has records of size %d, the BigArray has size %d"
               , k, s->aSize, aa->size) ;

  if (s->aMax > 0)
    {
      uBigArray (aa, s->aMax - 1) ;                 /* = arrayp (aa, iMax-1, TYPE): make room */
      if (fdReadWrite (scr->fd, aa->base, (size_t) s->aMax * s->aSize, s->pos, 0))
        messcrash ("scratchBigArrayGet block %d failed: %s", k, strerror (errno)) ;
    }
  arrayMax (aa) = s->aMax ;
  return aa ;
} /* uScratchGet */

/***************************/

