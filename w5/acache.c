/*  File: acache.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
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
    For acedb5, we are develloping a new disk handling and cache system
    using adisk.c acache.c bs2block.c
    removing disknew.c blocksubs.c and the lower level b*.c routines

    This cache handles A objects
    The idea is to keep them around for a while before sending them to
    disk so that multiple calls to arrayGet/Store, as happens during 
    parse operations and elsewhere cost nothing

    The B objects have their own B-develloped cache, that B-cache, handled
    in objcache.c stores its own objects as arrays but calls directly
    adisk.c and does not travel through the present cache

    Very large Arrays, larger than 1/2 the whole cache size are not cached.
    Also all accessed Arrays  are considered write enabled and no control is 
    exercised here. 

    We have to copy the array before serving it, because the client
    may be updating it, so it would be invalid for a second client, and may
    never store it back if he did not update it, so it would be lost.
    This is a bit costly, but certainly less than accessing the disk.

    We also have to copy the modified array when storing it, because the client
    may keep using it. This could be fixed by distinguishing arrayStore and
    arrayStoreDestroy in a few places in the code.

    Swapping is handled inside the cache, just when read/writing to disk

    Client updating is handled just when disk storing, assuming that acediff
    does not handle A object in any significant way, 
    I do not know if this is correct
 * Exported functions:
    BOOL aCacheGet
    void aCacheKill
    void aCacheStore
    void aCacheStoreDestroy 
    int aCacheSaveAll         // called from saveAll in adisk,c
 * HISTORY:
 * Last edited: Feb  8 14:29 1999 (edgrif)
 * Created: Wed Feb 25 1997
 *-------------------------------------------------------------------
 */

/* $Id: acache.c,v 1.4 2020/05/30 16:50:32 mieg Exp $ */

      /***************************************************************/
      /**  File acache.c :                                          **/
      /**  Handles the array cache     of the ACeDB program.        **/
      /***************************************************************/
      /***************************************************************/
      /*                                                             */
      /*  ? routines are public :                                    */
      /*         R.Durbin & J.Thierry-Mieg.                          */
      /*                    last modified  7/3/1991 by JTM.          */
      /*                                                             */
      /***************************************************************/

#include "acedb.h"

#include "lex_bl_.h"  /* for q->aCache */
#include "lex.h"
#include "pick.h"
#include "adisk.h"
#include "acache.h"
#include "client.h"

extern BOOL swapData ;
extern void swapArray(Array a, char *format) ;
extern int MAXKNOWNACACHE ;      /* cumulated array size, set in adisk,c  */
static int totalCache = 0, nKnownACache = 0, nAllocatedACache = 0 ;
typedef struct acacheStruct *ACACHE ;
struct acacheStruct { void *magic ; KEY key, parent ; Array a, b ; 
                 ACACHE up , next ; char format[32] ;
                 BOOL  hasBeenModified ; } ;
static ACACHE  knownACacheTop = NULL ,
	       knownACacheEnd  = NULL ;
static int acacheMagic = 0 ;
#define  ACACHEMAGIC  &acacheMagic 
#define SIZE_OF_ACACHE sizeof(struct acacheStruct)
extern KillFunc killFunc[] ;

static int naCache1read = 0, naCache1written = 0 , naCache2served = 0, naCache2saved = 0 ;

static void aCachePop(ACACHE v);
static void aCacheKOadd(ACACHE v);
static void aCacheKOremove(ACACHE v);

static ACACHE aCacheAlloc (void) ;     /* self managed calloc */
static void aCacheFree (ACACHE v);

/*****************************************************************/
/*****************************************************************/
                 /* General utility routines */
/*****************************************************************/
/*****************************************************************/
/********************    self managed calloc  ********************/         
static Stack freeStack = 0 ;

/***********************/

static ACACHE aCacheAlloc (void)   
{
  static int blocSize = 256 ;
  ACACHE p ;
  int i = -12 ;

  if (!freeStack)
    freeStack = stackCreate (4*blocSize) ;
  if (stackEmpty (freeStack))
    { p = (ACACHE) messalloc (blocSize * SIZE_OF_ACACHE) ;
      for (i = blocSize ; i-- ; ++p)
          push (freeStack,p,ACACHE) ;
/* printf ("Adding %d to ACACHE free list\n",blocSize) ; */
      nAllocatedACache += blocSize ;
      blocSize *= 2 ;
    }
  p = pop (freeStack,ACACHE) ;
  p->magic = ACACHEMAGIC ;
  return p ;
}

/***********************/

static void aCacheFree (ACACHE v)
{
  if (v->a)
    totalCache -= v->a->max * v->a->size ;
  arrayDestroy (v->a) ;
  memset (v, 0, SIZE_OF_ACACHE) ; 
  v->magic = 0 ;
  push (freeStack,v,ACACHE) ;
}

/*****************************************************************/
/*****************************************************************/

/************************************************************************/
                          /* Gives the aCache status*/
                          /*    Public;  Calls nothing*/
int aCacheStatus (int *used, int*known, int *modified, 
		 int *nc2rp, int *nc2wp, int *nccrp, int *nccwp)
{
  register ACACHE v;
  register int i, l, m;

  i = l = m = 0; v = knownACacheTop;
  while(v) 
    { i++ ;
      if (v && v->magic != ACACHEMAGIC)
	messcrash ("wrong achache in aCacheStatus") ;

      if(v->hasBeenModified)
	m++ ;
      v = v->next; 
    }
  *used = nAllocatedACache ;
  *known = i ;
  *modified = m ;
  *nc2rp = naCache2served ; *nc2wp = naCache2saved ;
  *nccrp = naCache1read ; *nccwp = naCache1written ;
  
  return TRUE;
}
/**************************************************************/


void blockStatus(int *nc1rp, int* nc1wp, int *ncc1rp, int* ncc1wp)
{
  *nc1rp =  naCache2served ; *nc1wp = totalCache >> 10 ; /* naCache2saved ;  */
  *ncc1rp = naCache1read ; *ncc1wp = naCache1written ;
}


/**************************************************************/
/**************************************************************/
    /* Linked Lists manipulations */

/************************************************************************/
 /*********************************************************************/
                /* to add an aCacheEntry to the Known ACacheEntry List */
                /*  Static; called by aCacheCreate and aCacheSave        */

static void aCacheKOadd(ACACHE v)
{
  ACACHE w ;
  
  nKnownACache++;
  w = knownACacheTop;
  knownACacheTop=v;
  
  if (!knownACacheEnd)
    knownACacheEnd = v;
  
  if (w) 
    w->up = v ;
  
  v->up =  NULL;
  v->next = w ; 
}

/*********************************************************************/
                /* to remove an aCacheEntry to the Known ACacheEntry List */
                /*  Static; called by aCacheUpdate and aCacheNew        */
                /*          Calls nothing.      */
static void aCacheKOremove(ACACHE v)
{
  ACACHE u, w;
  
  nKnownACache--; 
  u=v->up;
  w=v->next;
  
  if (knownACacheTop == v) 
    knownACacheTop = w ;
  
  if(knownACacheEnd == v)
    knownACacheEnd = u ;
  
  if(u)
    u->next = w;
  if (w) 
    w->up = u;
}

/**************************************************************/
                      /*pops a aCache to the top of the knownACache line*/
                       /*  Static; called by aCacheCreate; Calls nothing*/
static void aCachePop(ACACHE v)
{
  aCacheKOremove(v) ;
  aCacheKOadd(v) ;
}

/**************************************************************/

static void  aCacheDiskStore(ACACHE v, KEY key, KEY parent, 
			     Array a, char *format, BOOL destroy)
{
 
  Array b = 0 ;
  if (v)
    {
      v->hasBeenModified = FALSE ;
      a = v->a ;
      key = v->key ;
      format = v->format ;
      parent = v->parent ; 
      if (destroy)
	{ 
	  totalCache -= v->a->max * v->a->size ;
	  v->a = 0 ;
	}
    }

  if (swapData && format && *format)
    {
      b = destroy ? a : arrayCopy (a) ;
      swapArray(b, format);
      aDiskArrayStore (key, parent, 0, b) ;  /* no swap, use second array */	 
      arrayDestroy (b) ; 
      if (destroy) a = 0 ;
    }
  else
    {
      aDiskArrayStore (key, parent, 0, a) ;  /* no swap, use second array */
      if (destroy)  arrayDestroy (a) ; 
    }  
}

static void aCacheKOMakeRoom (int size)
{
  ACACHE v = knownACacheEnd, w ;
  LEXP q ;
  while (v && totalCache + size > MAXKNOWNACACHE)
    { w = v->up ;
      if(v->hasBeenModified) 
	aCacheDiskStore (v,0,0,0,0,TRUE) ;
      q = KEY2LEX(v->key);
      if (q && q->cache == v) q->cache = 0 ;
      aCacheKOremove(v);
      aCacheFree(v);
      v = w ;
    }
}

/**************************************************************/

static BOOL  aCacheKOPut (KEY key, KEY parent, Array a, char *format, BOOL modif)
{
  LEXP q = KEY2LEX(key);
  ACACHE v = q ? (ACACHE)q->cache : 0 ;  

  if (v && v->magic != ACACHEMAGIC)
    messcrash ("wrong achache in aCacheKOPut(%s)", name(key)) ;

  if (q) q->cache = 0 ;
  if (v) { aCacheKOremove(v); aCacheFree(v) ;}

  q = KEY2LEX(key); 
  if (format && strlen(format) > 32)
    return FALSE ;  /* do not store monsters */
  if (a->dim * a->size * 2 > MAXKNOWNACACHE)
    return FALSE ;       /* do not store giants */
  if (!q)
    return FALSE ;       /* do not store anonymous */

  aCacheKOMakeRoom(a->dim * a->size) ;
  v = aCacheAlloc() ; 

  v->key = key ;
  v->parent = parent ;
  v->hasBeenModified = modif ;
  strncpy(v->format, format, 32) ;
  v->a = arrayCopy (a) ;
  totalCache += v->a->max * v->a->size ;
  q = KEY2LEX(key) ;   /* reevealuate, may have moved */
  if (!q) messcrash("acacheKOput cannot register in q") ;
  q->cache = v ;
  aCacheKOadd(v) ;
  return TRUE ;
}

/**************************************************************/

static void aCacheDoStore (KEY key, KEY parent, Array a, char *format, BOOL destroy)
{
  if (externalSaver)   /* dump old version */
    externalSaver (key, 1) ;
  
  if (!aCacheKOPut (key, parent, a, format, TRUE))
    aCacheDiskStore (0, key, parent, a, format, destroy) ;

  if (externalSaver)   /* dump new  version and acediff it */
    externalSaver (key, 2) ;
  
  lexUnsetStatus (key, EMPTYSTATUS) ;
  lexSetStatus (key, TOUCHSTATUS) ;
}

/**************************************************************/
/********************  public functions ***********************/
/**************************************************************/

/* swapping NOT treated */
void aCacheKill (KEY key)
{
  LEXP q = KEY2LEX(key);
  ACACHE v = q ? (ACACHE)q->cache : 0 ;
  
  if (v && v->magic != ACACHEMAGIC)
    messcrash ("wrong achache in aCacheKOKill(%s)", name(key)) ;

  if (killFunc[class(key)])
    killFunc[class(key)] (key) ;

  if (q) q->cache = 0 ;
  if (v)
    { 
      lex2clear(v->key) ;
      aCacheKOremove(v) ; 
      aCacheFree(v);
      lexkill(key) ;     
    }
  aDiskKill (key) ;  
  lexSetStatus (key, EMPTYSTATUS) ;
  lexSetStatus (key, TOUCHSTATUS) ;  
  if (externalSaver)   /* dump old version */
    externalSaver (key, 41) ;
}

/**************************************************************/

BOOL aCacheGet (KEY key, KEY *pp, Array *ap, char *format)
{
  LEXP q = KEY2LEX(key);
  ACACHE v = q ? (ACACHE)q->cache : 0 ;

  if (v && v->magic != ACACHEMAGIC)
    messcrash ("wrong achache in aCacheKOGet(%s)", name(key)) ;
  if (v)
    {
      aCachePop (v) ;
      if (pp) *pp = v->parent ;
      *ap = arrayCopy (v->a) ;
      return TRUE ;
    }

  if (externalServer) externalServer (key, 0, 0, 0) ;
  if(!iskeyold(key))
    return 0;
  if (!aDiskArrayGet (key, pp, 0, ap)) return 0 ;

  if (swapData && format && *format)
    swapArray(*ap, format);
  aCacheKOPut (key, pp ? *pp : 0, *ap, format, FALSE) ;
  
  return TRUE ;
}
  

/**************************************************************/
 
void aCacheStore (KEY key, KEY parent, Array a, char *format)
{
  aCacheDoStore (key, parent, a, format, FALSE) ;
}

void aCacheStoreDestroy (KEY key, KEY parent, Array a, char *format)
{
  aCacheDoStore (key, parent, a, format, TRUE) ;
}

/**************************************************************/

int aCacheSaveAll(void)
{
  ACACHE v = knownACacheEnd, w ;
  int n = 0, total = 0 ;

  while (v)
    { w = v->up ; 
      total += v->a->max * v->a->size ;
      if(v->hasBeenModified) 
	{
	  n++ ;
	  aCacheDiskStore(v,0,0,0,0,FALSE) ;
	}
      v = w ;
    }
  if (total != totalCache)
    messcrash ("bad total count in aCacheSaveAll") ;
  return n ;
}

/**************************************************************/
/**************************************************************/
