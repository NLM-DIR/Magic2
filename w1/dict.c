/*  File: dict.c
 *  Author: Richard Durbin and Jean Thierry-Mieg (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1995
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: 
 * Exported functions:
 * HISTORY:
	   DICT library was rewritten nov 2002 to ensure
	   that the char* returned by dictName are valid
	   untill the whole dict is destroyed.
	   i.e. they are never reallocated even if you
	   keep adding names foreever
	   and dictRemove was added
 * Last edited: Dec  4 14:50 2002 (mieg)
 * Created: Tue Jan 17 17:33:44 1995 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: dict.c,v 1.29 2020/06/12 01:58:42 mieg Exp $ */
#include "regular.h"
                             /* every thing here is private */

static char *DICT_MAGIC = "dict-magic" ;
typedef struct { unsigned curr ; int pos, free, max ; char *base ; } DVOC ;
typedef struct dict_struct {
  BOOL caseSensitive ;          /* default FALSE, public */
  int dim ;         /* dimension of the hash table, has to be a power of 2 */
  int max ;           /* 2 to the  dict->dim */
  int count ;         /* number of active names in the dict */
  int nVoc ;           /* current dVoc */
  int newPos ;          /* moving index in the hash table */
  char *magic ;
  Array table ;			/* hash table */
  Array keys ;		/* ofsets in the set of vocs */
  Array dVocs ;                  /* array of DVOC */
  AC_HANDLE handle ;
  void * lock ;
  void (*channelLock)(void* lock) ;
  void (*channelUnlock)(void* lock) ;
} DICT_STRUCT ;

#include "dict.h"
typedef DICT _DICT ;

/*
  The names are stored incrementally in large buffers, dVoc->buff, which
  are never reallocated, so that the pointers returned by dictName
  remain valid until the dict is destroyed.

  The key is as usual in acedb a composite 1 byte for the dVoc number
  3 bytes for the offsett in dVoc->buff. In reality we divide by 8
  the offset since we align the names on 64 bit boundaries.

  The table is a double hashing hash table, its dimension is a power of 2
  so it is prime with the delta hash value which is always odd. In this
  way the orbit (first hash value modulo delta covers the whole table).

  The first 0 found during the travle indicates that the name is absent
  in the hash table. The first negative value indicates a spot than has 
  been freed and can be reused, but does not stop the bouncing search .

  We retrieve a name by direct dereferencing.
  We recognize a name by hash search.
*/

/* #define DEBUG_MODE */

#ifndef DEBUG_MODE
/* optimised mode */
/* these 2 private complex macros are useful to accelerate the hashing code */
/* whereas the public interface always checks the range of the parameters */

#define dictKey2Name(_dict,_k) ((arrp (_dict->dVocs, (((_k) >> 24) & 0x000000ff) - 1, DVOC))->base + ((0x00ffffff & (_k)) << 3))
#define dictIndex2Name(_dict,_i) (dictKey2Name((_dict),arr((_dict)->keys,(_i),KEY)))

#else
/* debugging mode */
static const char *dictKey2Name (_DICT *dict, KEY k) 
{
  int nVoc = ((k >> 24) & 0x000000ff) - 1 ;
  int pos = (0x00ffffff & k) << 3 ;
  DVOC *dVoc ;
  
  if (nVoc < 0 || nVoc >= arrayMax (dict->dVocs))
    messcrash ("uDictKey2Name bad nVoc = %d", nVoc) ;
  dVoc = arrp (dict->dVocs, nVoc, DVOC) ;
  if (pos >= dVoc->curr)
    messcrash ("uDictKey2Name bad pos = %d >= curr = %d", pos, dVoc->curr) ;
  return dVoc->base + pos ;
}

static const char *dictIndex2Name (_DICT *dict, int ii)
{
  return dictName (dict, ii) ;
}
#endif /* DEBUG_MODE */

/************* standard utility from Jean *************/
#define  SIZEOFINT (8 * sizeof (int))
static unsigned int dictHash (const char *cp, int n, BOOL isDiff)
{
  register int i ;
  register unsigned int j, x = 0 ;
  register int rotate = isDiff ? 21 : 13 ;
  register int leftover =  SIZEOFINT - rotate ;

  while (*cp)
    x = ace_upper (*cp++) ^ ((x << rotate) | ( x >> leftover)) ; 

				/* compress down to n bits */
  for (j = x, i = n ; i < SIZEOFINT ; i += n)
    j ^= (x >> i) ;
  j &= (1 << n) - 1 ;

  if (isDiff)			/* return odd number */
    j |= 1 ;

  return j ;
}  /* dictHash */

/************* Oleg Khovaykov *************/
/* Khovayko, Oleg  435-5885 olegh@ncbi, bulding 45 5 AN12D-40 1 */
/* mieg:2007/09/06, Oleg came by chance and gave me this code */
/* he says, buckets of size 3 to 5 would use less memory */
#define  SIZEOFINT (8 * sizeof (int))

/* roll is a magic construct that gets compiled 
   into a single roll assembly instruction 
   
   to verify: write a  sample code

unsigned int hash(const char *s) {
  unsigned int h = 1;
  unsigned char c;
  while(c = *s++) {
   h = ((h << 7) | (h >> (32 - 7))) + c;
  }
  return h;
}

then compile in in assembly and look
   cc -S x.c
   vi x.s
you will see the instruction  rorl    $25, %eax
*/
#define roll(_h) ((_h << 7) | (_h >> (64 - 7)))
#define BARKER 0x1111100110101
/* barker code have very low autocorrelation properties 
 * see wikipedia
 */
#ifdef JUNK

static unsigned int dictHashOleg (const char *s, int n, BOOL isDiff)
{ /* contributed by Oleg Khovayko, 2007 */
  unsigned long h = (unsigned int)BARKER ;
  unsigned char c;

   /* const int *TAB = isDiff ? (const int *)dictHashOleg : (const int *)stackTokeniseTextOn ;  a random set ? */
  static unsigned long T1[256], T2[256] ;
  static int T= 0 ;
  unsigned long *TAB = isDiff ? T1 : T2 ;

  if (T==0)
    {
      unsigned long k = (unsigned long) BARKER ;
      for (T = 0 ; T < 256 ; T++)
	{ T1[T] = k ^ ((unsigned long) 0x532814ea) ; k = roll(k) ; T2[T]= k ^  ((unsigned long)0x47e2132) ; } 
    }
  while ((c = *s++)) 
    h = roll(h) + (TAB[0xff & ((unsigned char)(h ^ c))] ^ c) ;
  h ^= h >> 32 ; /* accumulate the bits on the right */
  h ^= h >> 16 ; /* accumulate the bits on the right */
  h ^= (h >> 8) ;
  h &= (1 << n) - 1 ;	/* compress down to n bits */
  if (isDiff)			/* return odd number */
    h |= 1 ;

  return (unsigned int) h ;
}  /* dictHashOleg */

#endif

/****************************************************************************/

static void dictReHash (_DICT *dict, int newDim)
{
  int ii ;
  KEY *kp ;

  if (newDim <= dict->dim)
    return ;
  dict->dim = newDim ;
  dict->max = 1 << newDim ;
				/* remake the table */
  arrayDestroy (dict->table) ;
  dict->table = arrayHandleCreate (dict->max, int, dict->handle) ;
  array (dict->table, dict->max-1, int) = 0 ;	/* set arrayMax */

				/* reinsert all the names */
  for (ii = 1, kp = arrp(dict->keys, ii, KEY)  ; ii < arrayMax(dict->keys) ; ii++, kp++)
    { dictFind (dict, dictKey2Name (dict, *kp), 0) ;
				/* will fail, but sets dict->newPos */
      arr(dict->table, dict->newPos, int) = ii ;
    }
} /* dictReHash */

/****************************************************************************/

static DVOC *dictAddVoc (_DICT *dict, int s)
{
  DVOC *dVoc ;
  int n1 = 0 ;

  dVoc = arrayp (dict->dVocs, dict->nVoc, DVOC) ;
  if (s >  0x00ffffff << 3)    /* 134M */
    messcrash ("You cannot enter in a dict a word longer than 130 Million characters, sorry") ;
  if (dict->nVoc > 200)
    n1 = 0x80000000 ;   /* 2G */
  else if (dict->nVoc > 3)
    n1 = 0x40000000 ; /* 1G */
  else if (dict->nVoc > 2)
    n1 = 0x10000000 ; /* 250M */
  else if (dict->nVoc > 1)
    n1 = 0x01000000 ; /* 16M */
  else if (dict->nVoc > 0)
    n1 = 0x0100000 ; /* 1M */
  else
    n1 = 0x1000 ;    /* 4k */
  if (s > n1)
    n1 = s ;

  if (n1 >  0x00ffffff << 3)
    n1 =  0x00ffffff << 3 ; /* 134M = maximal offset which can be stored in a key */

  dVoc->base = (char *) halloc (n1, dict->handle) ; /* never realloc */
  dVoc->curr = 0 ;
  dVoc->max = n1 ;
  dVoc->free = n1 ;
  
  dict->nVoc++ ;
  if (dict->nVoc > 255)
    messcrash ("The dictionary is full, it already contains 256 table of 8Gb, please check  w1/dict.c") ;
  return dVoc ;
}  /* dictAddVoc */


void uDictDestroy (_DICT *dict) 
{ 
  if (dict && dict->magic == DICT_MAGIC) 
    {  /* we need 2 lines, because otherwise messfree sets (freed) dict->handle=0 */
      AC_HANDLE handle = dict->handle ;
      messfree (handle) ;
    }
}

/****************************************************************************/

_DICT *dictHandleCreate (int size, AC_HANDLE handle)
{
  AC_HANDLE h = handleHandleCreate (handle) ;
  _DICT *dict = (_DICT*) halloc (sizeof (_DICT), h) ;

  dict->handle = h ;
  dict->magic = DICT_MAGIC ;
  dict->caseSensitive = FALSE ;
  dict->count = 0 ;
  for (dict->dim = 6, dict->max = 64 ; 
       dict->max < size && dict->max < 0x00ffffff ; 
       ++dict->dim, dict->max *= 2) ;
  dict->table = arrayHandleCreate (dict->max, int, dict->handle) ;
  array (dict->table, dict->max-1, int) = 0 ;	/* set arrayMax */
  dict->dVocs = arrayHandleCreate (8, DVOC, dict->handle) ;
  dictAddVoc (dict, 8 * dict->max) ;  /* wild guess average word length will be 16 bytes */
  dict->keys = arrayHandleCreate (dict->dim/4, KEY, dict->handle) ;
  array (dict->keys, 0, KEY) = 0 ;		/* reserved for empty table entry */
  return dict ;
} /* dictHandleCreate */

_DICT *dictCreate (int size) 
{
  return dictHandleCreate (size, 0) ; 
} /* dictCreate */

_DICT *dictCaseSensitiveHandleCreate (int size, AC_HANDLE handle) 
{
  _DICT *dict = dictHandleCreate (size, handle) ;
  dict->caseSensitive = TRUE ;
  
  return dict ;
} /* dictCaseSensitiveHandleCreate */

/****************************************************************************/

void dictSetLock (DICT *dict, void *lock, void (*cLock)(void*), void (*cUnlock)(void*))
{
  dict->lock = lock ;
  dict->channelLock = cLock ;
  dict->channelUnlock = cUnlock ;
} /* dictSetLock */

/****************************************************************************/

BOOL dictFind (_DICT *dict, const char *s, int *ip)
{
  register int ii ;
  register unsigned int h, dh = 0 ;
  int (*mystrcmp)(const char *, const char *) = dict->caseSensitive ?
    strcmp : strcasecmp ;

  if (!dict || !s || !*s)
    return FALSE ;
  int len = strlen (s) ;
  if (len >  0x00ffffff << 3)    /* 134M */
    messcrash ("You cannot enter in a dict a word longer than 130 Million characters, sorry") ;

  dict->newPos = 0 ; /* will become first reusable spot */
  h = dictHash (s, dict->dim, FALSE) ;
  
  while (TRUE)
    { ii = arr (dict->table, h, KEY) ;
      if (!ii)			/* empty slot, s is unknown */
	{ 
	  if (ip) 
	    *ip = - 1 ;
	  if (dict->newPos == 0)
	    dict->newPos = h ;
	  return FALSE ;
	}
      else if (ii < 0)		/* freed stop */
	{ 
	  if (dict->newPos == 0)
	    dict->newPos = h ;
	  /* continue ; bug: provokes an infinite loop, we do want to increment h */
	}
      else if (!mystrcmp (s, dictIndex2Name(dict,ii)))
	{ if (ip) 
	    *ip = ii ;
	  dict->newPos = h ;
	  return TRUE ;
	}
      if (!dh)
	dh = dictHash  (s, dict->dim, TRUE) ;
      h += dh ;
      if (h >= dict->max)
	h -= dict->max ;
    }
} /* dictFind */

/****************************************************************************/

BOOL dictRemove (_DICT *dict, const char *s) 
{
  int ii = 0 ;

  if (!dict || 
      !s ||
      !dictFind (dict, s, &ii))	/* word unkown */
    return FALSE ;

  ii++ ;
  arr (dict->keys, ii, KEY) = 0 ;  /* will not be rehashed */
  arr (dict->table, dict->newPos, int) = -ii ; /* will be reusable */
  dict->count-- ;
  return TRUE ;
} /* dictRemove */

/****************************************************************************/
/* always fills ip, returns TRUE if added, FALSE if known */
BOOL dictAdd (_DICT *dict, const char *s, int *ip)
{
  int ii = 0, len ;
  DVOC *dVoc  ;
  char *buf ;

  if (! dict || dict->magic != DICT_MAGIC)
    messcrash ("Wrong call to dictAdd, the dict is non existent") ;

  if (!s || ! *s)
    return FALSE ;

  if (dictFind (dict, s, &ii))	/* word already known */
    {
      if (ip)
	*ip = ii ;
      return FALSE ;
    }

  if (dict->lock) 
    {
      dict->channelLock (dict->lock) ;
      if (dictFind (dict, s, &ii))	/* word added by a parallel thread */
	{
	  if (ip)
	    *ip = ii ;
	  dict->channelUnlock (dict->lock) ;
	  return FALSE ;
	}
    }
  ii++ ;
  if (ii < 0)
    ii = -ii ; /* reuse */
  else
    ii = arrayMax(dict->keys) ;
  array (dict->table, dict->newPos, int) = ii ;
  
  buf = strnew (s, 0) ; /* needed if s in users space is a pointer inside dict */
  dVoc = arrayp (dict->dVocs, dict->nVoc - 1, DVOC) ;
  len = strlen (buf) ;
  len++ ;            /* count the terminal zero */
  while (len%8) len++ ; /* adjust on word boundary */
  if (len + 1 >= dVoc->free)
    dVoc = dictAddVoc (dict, len) ;

  array (dict->keys, ii, KEY) = (0xff000000 & ((dict->nVoc)<<24)) | (0x00ffffff & (dVoc->curr)>>3) ;
  strcpy (dVoc->base + dVoc->curr, buf) ;
  dVoc->curr += len ;
  dVoc->free -= len ;
  dict->count++ ;
  if (arrayMax(dict->keys) > 0.4 * dict->max)
    dictReHash (dict, dict->dim+1) ;

  messfree (buf) ;
  if (ip)
    *ip = ii ;

  if (dict->lock) dict->channelUnlock (dict->lock) ;
  return TRUE ;
} /* dictAdd */

/********************** utilities ***********************/
/* dictName returns a pointer that never gets reallocated */

const char *dictName (const _DICT *dict, int ii)
{
  KEY key ;

  if (ii < 0 || ii >= arrayMax(dict->keys))
    {
      messcrash ("Call to dictName() out of bounds: %d not in [1,%d] ", ii, dictMax(dict) - 1) ;
      return "(Dict error, NULL NAME)" ;
    }
  else if (ii == 0)
    {
      messcrash ("Call to dictName() out of bounds: %d not in [1,%d] ", ii, dictMax(dict) - 1) ;
      return "(NULL)" ;
    }

  key = arr (dict->keys, ii, KEY) ;
  return dictKey2Name(dict, key) ;
} /* dictName */

int dictCount (const _DICT *dict)
{
  return dict->count ;   /* number of active names */
}  /* dictMax */

int dictMax (const _DICT *dict)   /* max to be used if looping on all entries [1,dictMax] inclusively */
{
  return arrayMax(dict->keys) - 1 ;   /* 0 == reserved pseudo key */
}  /* dictMax */

void dictStatus (const _DICT *dict)
{
  fprintf(stderr, "DICT status: count = %d max=%d nTables=%d ndVocs=%d nKeys=%d\n"
	  ,  dict->count, dict->max
	  , arrayMax(dict->table)
	  , arrayMax(dict->dVocs)
	  , arrayMax(dict->keys)
	  ) ;
}

/******************** end of file **********************/
/****************************************************************************/
/****************************************************************************/
