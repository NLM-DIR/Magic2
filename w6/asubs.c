/*  File: asubs.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 **  Creates saves and updates permanent Arrays and Stacks
 * Exported functions:
 * HISTORY:
 *   added u format for BSunits, assumes swaps like k
 * Last edited: Dec  2 00:08 1996 (rd)
 * Created: Wed Jan 22 15:38:55 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: asubs.c,v 1.12 2017/10/02 19:48:53 mieg Exp $ */

#define DEF_BP /* avoid typedefing it as void* in disk.h */
typedef struct ablock *ABLOCK ;
typedef ABLOCK BP ;

#include "acedb.h"
#include "a_.h"
#include "disk_.h"  /* for size d */
#include "block.h"
#include "../wh/byteswap.h"
#include "client.h"
#include "acache.h"
#include "bigarray.h"

/**************************************************************/
static int naread = 0, nawritten = 0 ;  /* for the status info */

void aStatus(int *nrp, int *nwp)
{ *nrp = naread ; *nwp = nawritten ;
}

/**************************************************************/

/* Determine how the compiler aligns data types in structs */
/* There are also some of these in mystdlib.h for short, int, etc. */

#define KEY_ALIGNMENT (sizeof(struct { char c; KEY k; }) - sizeof(KEY))
#define DISK_ALIGNMENT (sizeof(struct { char c; DISK d; }) - sizeof(DISK))

void swapABlock(BP bp)
{ bp->hat.h.disk = swapDISK(bp->hat.h.disk);
  bp->hat.h.nextdisk = swapDISK(bp->hat.h.nextdisk);
  bp->hat.h.session = swapInt(bp->hat.h.session);
  bp->hat.h.key = swapKEY(bp->hat.h.key);
  bp->hat.size = swapInt(bp->hat.size);
  bp->hat.max = swapInt(bp->hat.max);
  bp->hat.localmax = swapInt(bp->hat.localmax);
  /*  bp->hat.parent = swapKEY(bp->hat.parent) ; */
}

void swapArray(Array a, char *format)
{ int i, offset = 0;
  char c, *vp;

  if (sizeof(float) != sizeof (KEY) ||
      sizeof(int) != sizeof (KEY) ||
      sizeof(DISK) != sizeof (KEY) ||
      sizeof(mytime_t) != sizeof (KEY))
    messcrash 
      ("%s%s", "On this platform BSdata and KEY differ",
       " in byte size, Rewrite the Swap code ! ") ;
    
  while ((c = *format++))
    switch (c)
      { 
      case 'c': 
	offset += 1;
	break;

/* We require that everything in bsData is 32 bits and swapable in the 
   same way, and that keys, ints and floats are in BSdata */	
      case 'k': 
      case 'a':   /* BSdata */
      case 'i': 
      case 'f': 
	while(offset%INT_ALIGNMENT) offset++;
	vp = a->base + offset;
	for (i= 0;  i<a->max; i++)
	  { *((int *)vp) = swapInt(*((int *)vp));
	    vp += a->size;
	  }
	offset += sizeof(float);
	break;

      case 'v':
	messcrash("swapArray : Cannot byte swap a pointer");
	break;

      case 'd':
	while(offset%DISK_ALIGNMENT) offset++;
	vp = a->base + offset;
	for (i= 0;  i<a->max; i++)
	  { *(DISK *)vp = swapDISK(*(DISK *)vp);
	    vp += a->size;
	  }
	offset += sizeof(DISK);
	break;
      }
}


 /* expected format looks like "fvvkkk" "iii" etc. */
int sizeOfFormattedArray(char *format)
{ int n = 0;
  int mostRestrictiveAlignment = 1;
  int s = 0, a = 0;
  char *cp = format ;
  
  if (!cp)
    messcrash("sizeOfFormattedArray received a null pointer") ;

  cp-- ;
  while (*++cp)
    {
	switch(*cp)
	  {
	  case 'c': s = 1 ; a = 1; break ; /* single char */
	  case 'a':
	  case 'k': s = sizeof(KEY); a = KEY_ALIGNMENT; break ;  
	  case 'i': s = sizeof(int); a = INT_ALIGNMENT; break ;  
	  case 'f': s = sizeof(float); a = FLOAT_ALIGNMENT; break ; 
	  case 'v': s = sizeof(void *); a = PTR_ALIGNMENT; break ;  
	  case 'd': s = sizeof(DISK); a = DISK_ALIGNMENT; break ; 
	  default:
	    messcrash("Unknown format %s at %s in sizeOfFormattedArray", 
		      format, cp) ;
	  }
	if (a > mostRestrictiveAlignment) mostRestrictiveAlignment = a;
	while (n%a) n++;
	n += s;
	
    }
  while (n%mostRestrictiveAlignment) n++;

  return  n;
}
	
/**************************************************************/

BOOL dumpFormattedArray (FILE *f, Array a, char *format)
{ int i ;
  char *cp = format ,  *cq ;
  KEY key ;

  if (!arrayExists(a) || !f)
    messcrash("null a or f in dumpFormattedArray") ;
 
  for (i = 0 ; i < arrayMax(a) ; i++)
    { fputc('\n', f) ;
      cq = a->base + i * a->size ;
      cp = format - 1 ;
      while(*++cp)
	{
	    switch (*cp)
	      {
	      case 'a':
		while (((long)cq)%KEY_ALIGNMENT) cq++;
		key = *(KEY*) cq ;
		fprintf(f, " bsunit ") ;
		cq += sizeof(KEY) ;
		break ;
	      case 'k':
		while (((long)cq)%KEY_ALIGNMENT) cq++;
		key = *(KEY*) cq ;
		fprintf(f, " %s : \"%s\" ",
			className(key),
			name(key)) ;
		cq += sizeof(KEY) ;
		break ;
	      case 'v':
		messcrash("dumpFormattedArray tries to dump a pointer") ;
		break ;
	      case 'i':
		while (((long)cq)%INT_ALIGNMENT) cq++;
		fprintf(f, " %d ", *(int*)cq) ;
		cq += sizeof(int) ;
		break ;
	      case 'f':
		while (((long)cq)%FLOAT_ALIGNMENT) cq++;
		fprintf(f, " %f ", *(float*)cq) ;
		cq += sizeof(float) ;
		break ;
	      case 'c':
		fputc(*cq++, f) ;
		break ;
	      default:
		messcrash("Unknown format %s at %s in dumpFormattedArray", 
			  format, cq ) ;
	      }
	}
    }
  fputc('\n', f) ;
  fputc('\n', f) ;	    
  return TRUE ;
}
	
/**************************************************************/

        /* No check if size = -1 */
Array uArrayHandleGet(KEY key, int size, char *format, AC_HANDLE handle)
{
  Array a=0;  
#ifndef ACEDB5
  BP bp ;
  int i = 0, j;
  BOOL doSwap = swapData ;
  char type = 'A' ;
#endif

  key = lexAliasOf(key) ;
  switch (pickType(key))
    {
    case 'A':
      if (!format || size == -1)
	messcrash("Format strings are no linger optional in calls to arrayGet");
      if (size != sizeOfFormattedArray(format))
	messcrash("Actual size of array  %d != size %d of format %s",
		  size, sizeOfFormattedArray(format), format) ;
      break ;
    default:
      messcrash ("Attempt to arrayGet a non-array key %s:%s",
		 className(key), name(key)) ;
    }
  
#ifdef ACEDB5
  if (aCacheGet (key, 0, &a, format)) naread++ ;
  return a ;
#else /***** ACEDB.1/2/3/4 system *******/
  if (externalServer) externalServer (key, 0, 0, 0) ;
  if(!iskeyold(key))
    return 0;
  naread++ ;

  blockpinn(key,&bp);

  if (bp->hat.h.type != type)
    messcrash("arrayGet(%s) read a block of h.type %c",
	      name(key),bp->hat.h.type);
  if (size != -1 && bp->hat.size != size)
    messcrash("arrayGet(%s) read an array of incorrect size %d != %d",
	      name(key),bp->hat.size,size);
  a =  uArrayCreate(bp->hat.max + (bp->hat.size == 1 ? 1 : 0) ,bp->hat.size, handle);
  arrayMax(a) = bp->hat.max;
  while (TRUE)
    {
      j = bp->hat.localmax ;
      memcpy(a->base + i*a->size, bp->n, j*a->size) ;
      i += j ;
      if (!blockNext(&bp))
	break ;
    }
  blockunpinn(key);
  if(i != arrayMax(a))
    messcrash("Size inconsistency in arrayGet (%s)",
	      name(key));

  
  if (doSwap)
    swapArray(a, format);

  return a;
#endif  
}

BigArray uBigArrayHandleGet(KEY key, int size, char *format, AC_HANDLE h)
{
  BigArray b = 0 ;
  Array a = uArrayHandleGet (key, size, format, 0) ; 
  if (a && arrayMax (a))
    {
      int dx = size == 1 ? 1 : 0 ; /* zero terminate the dna */
      b = uBigArrayCreate (arrayMax (a) + dx, size, h, TRUE) ;
      bigArrayMax (b) = arrayMax (a) ;
      memcpy(b->base, a->base, a->max * a->size) ;
    }
  arrayDestroy (a) ;
  return b ;
}

BigArray uBigArrayGet(KEY key, int size, char *format)
{ return uBigArrayHandleGet (key, size, format, 0) ; }
  
/**********************************************************/

void arrayStore(KEY key, Array a, char* format)
{ 
#ifndef ACEDB5
  BP bp;
  int i = 0 , local, k;
  char type = 'A' ;
#endif
  int j = arrayExists(a) ? a->max : 0 ;

  key = lexAliasOf(key) ;
  
  switch (pickType(key))
    {
    case 'A':
      if (!format)
	messcrash("Format string not provided in call to arrayStore");
      if (a->size != sizeOfFormattedArray(format))
	messcrash("Actual size of array  %d != size %d of format %s",
		  a->size, sizeOfFormattedArray(format), format) ;
      break ;
    default:
      messcrash ("Attempt to arrayStore a non-array key %s:%s",
		 className(key), name(key)) ;
    }
  nawritten++ ;

  sessionAutoSave (-8, 0) ; /* register date of last modif */

#ifdef ACEDB5
  if (j)
    aCacheStore (key, 0, a, format) ;  /* no swap, use second array */ 
  else
    arrayKill(key) ;
#else      /***** ACEDB.1/2/3/4 system *******/

  if (externalSaver)   /* dump old version */
    externalSaver (key, 1) ;

  if (j)
   { 
     if (swapData)
       swapArray(a, format);
     blockpinn (key, &bp);
     
     local = AMAX / a->size ;
     while (TRUE)
       {
	 bp->hat.h.key = key;
	 bp->hat.h.type = type ;
	 bp->hat.size = a->size;
	 bp->hat.max  = a->max;

	 k = bp->hat.localmax =  (j>local) ? local : j ;
	 memcpy(bp->n, a->base + i*a->size, k*a->size) ;
	 if (k < local)
	   memset(bp->n + k*a->size, 0, (local-k)*a->size) ;
	 if (local*a->size < AMAX)  /* mieg, 2002, vanquish uninit memory read warning of purify */
	   memset(bp->n + local*a->size, 0, AMAX - local*a->size) ;

	 j -= k ;
	 i += k ;
	 if(!j) 
	   {
	     blockSetEnd(bp) ;
	     break ;
	   }
	 blockSetNext(&bp) ;
       }
     blockrewrite(key) ;

     if (swapData)
       swapArray(a, format); /* put it back now */  
     
     if (externalSaver)   /* dump new  version and acediff it */
       externalSaver (key, 2) ;

     lexAlphaMark (class (key)) ;
     lexUnsetStatus (key, EMPTYSTATUS) ;
     lexSetStatus (key, TOUCHSTATUS) ;
   }
 else
   arrayKill(key) ;
#endif
}

void bigArrayStore(KEY key, BigArray b, char* format)
{
  Array a = 0 ;
  long int n1 = b->max, n2 = n1 * b->size ;

  if (n1 * b->size > 0x8fffffff)
    messcrash ("bigArrayStore  received an array which is too large ") ;
  a = uArrayCreate (n1, b->size, 0) ;
  a->max = n1 & 0x8fffffff ;
  memcpy (a->base, b->base, n2) ;
  arrayStore (key, a, format) ;
  arrayDestroy (a) ;
}

/**********************************************************/

  /* it is assumed here that array has 
     no entry in the objcache
     */

KillFunc killFunc[256] ;  /* == MAXTABLE of lexsubs.c */

static void arrayDoKill (KEY key)
{ 
#ifdef ACEDB5
  aCacheKill (key) ;
#else
  BP bp ;

  if (killFunc[class(key)])
    killFunc[class(key)] (key) ;
  if (iskeyold(key)) 
    { 
      blockpinn(key, &bp); /* unpinned by blockSetEmpty */
      blockSetEmpty(bp) ;
      lexkill(key) ;      
    }
  lexAlphaMark (class (key)) ;
  lexSetStatus (key, EMPTYSTATUS) ;
  lexSetStatus (key, TOUCHSTATUS) ;  
  if (externalSaver)   /* dump old version */
    externalSaver (key, 41) ;
#endif
}

void arrayKill (KEY key)
{
  key = lexAliasOf(key) ;
  if (pickType(key) != 'A')
    messcrash ("Attempt to arrayKill a non-array key %s:%s",
	       className(key), name(key)) ;
  
  arrayDoKill (key) ;
}

/*********************************************/
/*********************************************/

void stackStore(KEY key, Stack s)
{
 if(!stackExists(s))
   messcrash("stackStore received a non initialised stack");
 if (!s->textOnly)
   messcrash("stackStore got non text only stack");
 arrayMax(s->a) = stackMark(s) ;
 arrayStore(key, s->a,"c") ;
}

/*********************************************/

Stack stackHandleGet(KEY key, AC_HANDLE handle)
{
  Stack s ;
  Array a = arrayHandleGet(key, char,"c", handle) ;
  
  if(!a)
    return 0 ;
  
  if(a->size != sizeof(char))
    messcrash("stackGet read a non char array") ;
  
  s = stackHandleCreate(32, handle);
  arrayDestroy(s->a) ;
  s->a = a ;
  s->pos = s->ptr = s->a->base + arrayMax(s->a) ;
  s->safe = s->a->base + s->a->dim - 8 ; 
  s->textOnly = TRUE;

  return s ;
}

Stack stackGet(KEY key)
{ return stackHandleGet (key, 0) ; }

/**********************************************************/
/**********************************************************/
   /* These new routines are designed for tableGetStore
    * they allow to store under a single key an array and a stack
    */

void arrayStackStore(KEY key, Array a, char *format, Stack s) 
{
  int na, as, ns, bMax ;
  Array b = 0 , hh ;  /* hh is a header array */
  char *cp ;

  key = lexAliasOf(key) ;
  if (pickType(key) != 'A')
    messcrash ("Attempt to arrayStackStore a non-array key %s:%s",
	       className(key), name(key)) ;

  if (!format)
    messcrash("Format string not provided in call to arrayStackStore");
  if (a->size != sizeOfFormattedArray(format))
   messcrash("Actual size of array  %d != size %d of format %s",
	     a->size, sizeOfFormattedArray(format), format) ;
  if (!arrayExists(a) || !stackExists(s))
    messcrash (" arrayStackStore received a non valid array or stack") ;

  na = arrayMax(a) ; as = a->size ; ns = stackMark(s) ;

  hh = arrayCreate (3, KEY) ;
  array (hh, 0, KEY) = na ;   array (hh, 1, KEY) = as ;   array (hh, 2, KEY) = ns ; 
 
  if (na && swapData)
    swapArray(a, format);
  if (swapData) swapArray (hh, "k") ;
  
  bMax = hh->size * hh->max + a->size * a->max + ns ;
  b = arrayCreate (bMax, char) ;

  cp = b->base ;
  memcpy (cp, hh->base, hh->max * hh->size) ; cp +=  3 * hh->size ;
  memcpy (cp,  a->base,  na * as) ; cp +=  na * as ;
  memcpy (cp, s->a->base, ns) ; cp += ns ;
  b->max = cp - b->base ;
  if (b->max != bMax) messcrash ("arrayStackStore internal error") ;
  arrayStore(key, b,"c") ;

  arrayDestroy (b) ;
  arrayDestroy (hh) ;
}  

/**********************************************************/
/**********************************************************/

BOOL  arrayStackHandleGet (KEY key, Array* aap, char *format, 
			   Stack *sp, AC_HANDLE handle)
{
  Array a = 0, b = 0, hh = 0 ; Stack s = 0 ;
  char *cp ;
  int na, as, ns ;

  key = lexAliasOf(key) ;
  if (pickType(key) != 'A')
    messcrash ("Attempt to arrayStackHandleGet a non-array key %s:%s",
	       className(key), name(key)) ;
  if (!format)
    messcrash("Format string not provided in call to arrayStackHandleGet");

  b = arrayGet (key, char, "c") ;
  if (!b) return FALSE ;

  hh = arrayCreate (3, KEY) ;
  arrayMax(hh) = 3 ;
  
  /* recover the header */
  cp = b->base ; 
  memcpy (hh->base, cp, 3 * hh->size) ; cp += 3 * hh->size ;
  if (swapData) swapArray (hh, "k") ;

   na = array (hh, 0, KEY) ; as = array (hh, 1, KEY) ; ns =  array (hh, 2, KEY) ;

  if (as != sizeOfFormattedArray(format))
   messcrash("Actual size of array  %d != size %d of format %s",
	     as, sizeOfFormattedArray(format), format) ;



  /* recover the array */
  a =  uArrayCreate(na, as, handle) ;
  a->max = na ;
  memcpy (a->base, cp, na * as) ; cp += na * as ;
  if (swapData) swapArray (a, format) ;

  /* recover the satck */
  s = stackHandleCreate(32,handle);
  arrayDestroy(s->a) ;
  s->a = uArrayCreate(ns, 1, handle) ;
  s->pos = s->ptr = s->a->base + ns ;
  s->safe = s->a->base + s->a->dim - 8 ; 
  s->textOnly = TRUE;
  s->a->max = na ;
  memcpy (s->a->base, cp, ns) ; cp += ns ;

  if (b->max != cp - b->base) messcrash ("arrayStackGet internal error") ;
  arrayDestroy (b) ;
  arrayDestroy (hh) ;

  *aap = a ; *sp = s ;

  return TRUE ;
}

/**********************************************************/
/**********************************************************/
