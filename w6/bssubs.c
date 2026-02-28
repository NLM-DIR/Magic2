/*  file: bssubs.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
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
 * Description: interface to the 'B' classes - tree objects
 		public declarations in bs.h
		private stuff in bs_.h (restricted to other w6 files)
 * Exported functions:
      The following preceeded by obj - used by objcache.c
        Get, Store, StoreNoCopy, Destroy, Copy
      The following all preceeded by bs
 	Key
	Save
	GetTag, GetKey, GetClass, GetData, Flatten
 	AddTag, AddKey, AddData, AddComment
	Type, TypeCheck
	Mark, Goto
	Remove, Prune
	FuseMode	(treedisp.c)
	MakePaths	(model.c)
        TagsInClass, KeySet
 * HISTORY:
 * Last edited: Feb 21 11:42 2001 (edgrif)
 * * Nov 17 16:40 2000 (edgrif): Fix bug in alias which was introduced by
 *              recent fix for longstanding bug in xreffing.
 * * Sep  2 14:55 1998 (edgrif): Remove incorrect declaration of lexAlphaMark.
 * * Nov 24 22:39 1993 (rd): allow Xreffing from subobjects
 * * Apr 26 11:30 1993 (cgc): reintroduced bsAddKey(bsHere)
 * * Sep 23 15:56 1992 (mieg): replaced get/addObj by pushObj
 * * Sep 19 18:47 1992 (rd): fixed UNIQUE and removed FREE (== ANY REPEAT)
 * * Aug 27 17:30 1992 (mieg): scope of UNIQUE is now just 1 move right
 * * Aug 16 10:57 1992 (mieg): no XREF on models
 * * Jul 17 12:17 1992 (mieg): modified bsKeySet and added TagsInClass
 * * Feb 21 14:05 1992 (mieg): getObj, addObj etc - not complete yet (remove etc)
 * * Jan 13 16:47 1992 (mieg): more asn stuff, rd's bsfuse
 * * Jan  8 19:02 1992 (mieg): messcrah when model is missing (I guess in wspec)
 * * Oct 22 13:30 1991 (mieg): added bsVisible.. to better queryGrep searches
 * * Oct 16 18:04 1991 (rd): fixed bug in xref deleting
 * Created: Wed Oct 16 17:33:54 1991 (rd)
 * CVS info:   $Id: bssubs.c,v 1.35 2020/05/30 16:50:32 mieg Exp $
 *-------------------------------------------------------------------
 */

#include <ctype.h>

#define DEF_BP /* avoid typedefing it as void* in disk.h */
typedef struct bblock *BBLOCK ;
typedef BBLOCK BP ;
#define DEFINE_OBJ
typedef struct sobj *OBJ ;

#include "acedb.h"
#include "cache_.h"
#include "bs_.h"  /* private functions of bs package */
#include "cache.h"
#include "bs.h"   /* public functions of bs package */
#include "dump.h"
#include "parse.h"
#include "check.h"
#include "client.h"
#include "model.h"

BOOL isTimeStamps = TRUE ;		/* GLOBAL for using timeStamps */
BOOL XREF_DISABLED_FOR_UPDATE = FALSE ;

/* if a tag points to more than than MAXDOWN keys, we create an associator */
#define MAXDOWN 12   
int BSTEST = 0 ;
#undef DELETE_DEBUG

static void     OBJfree (OBJ obj) ;
static void     makePaths (Associator ass, BS bs, BOOL isUnique) ;
static BOOL     findTag (OBJ obj, KEY tag, BOOL makeIt) ;
static void	xrefPrune (OBJ obj, BS bs, BS bsm) ;

BOOL     isXrefing = TRUE ;      /* disabled during doXref
				  * Used in cacheUpdate() and xClientSave()
				  */

static Array tabModel = 0;

/* XQX
* Conjecture: tabAssPath is indexed by the object class number.  Look
* up tabAssPath[class number] to get an associator.  Look up a tag
* number in the associator to get another array.  In that array, there
* is some kind of information about the branches of the objects tree
* structure that you follow to get to the tag.
*/
static Array tabAssPath = 0;

#define   CHECK_MODEL(ccc) if (!tabModel ||                      \
			       !array(tabModel,KEYKEY(ccc),BS))  \
                                  bsMakePaths (ccc) ;

#define OBJ_MAGIC 294756

/* NB I assume that the first call to get or add functions in a
   sequence uses a tag, not _bsRight, _bsDown.  This means that I
   only need to check that an object is valid in bsFind/addTag.
   Also, I always leave x->curr pointing to a node that fits the
   model, and x->modCurr pointing to the respective model node.
*/

/**************************************************************/
/**************************************************************/
/* Manage OBJ allocation
   hope they speed things up/prevent fragmentation
   in any case, they will help monitoring memory usage
   NB they must be paired.
*/
/**************************************************************/
/**************************************************************/
static int nCreated = 0, nUpdated = 0, nDestroyed = 0, nSaved = 0, nKilled = 0 ;

void BsSubsStatus (int *used, int *locked)
{ *used = nCreated - nDestroyed ; *locked = nUpdated - nSaved - nKilled ;
}

static AC_HANDLE OBJhandle = 0 ;
static Stack freeOBJstack = 0 ;
static int nOBJused = 0, nOBJalloc = 0 ;        /* useful for debug */

OBJ OBJalloc (void)       /* self managed calloc */
{
  static int blocSize = 2048 ;
  OBJ p ;
  int i ;

  if (!freeOBJstack)
    {
      OBJhandle = handleCreate () ;
      freeOBJstack = stackHandleCreate (4*blocSize, OBJhandle) ;
    }
  if (stackEmpty (freeOBJstack))
    { p = (OBJ) halloc (blocSize * sizeof(struct sobj), OBJhandle) ;
      for (i = blocSize ; i-- ; ++p)
        push (freeOBJstack,p,OBJ) ;
      nOBJalloc += blocSize ;
      blocSize *= 2 ;
    }
  p = pop (freeOBJstack,OBJ) ;
  memset (p, 0, sizeof (struct sobj)) ;
  ++nOBJused ;
  p->magic = OBJ_MAGIC ;
  return p ;
}

/* OBJfree must just release memory, not manage cache etc.
   see bsDoSave()
*/

static void OBJfree (OBJ obj)
{
 if (!obj)
    return ;
  if (obj->magic == OBJ_MAGIC)
    { obj->magic = 0 ;
      stackDestroy(obj->xref) ;
      push (freeOBJstack,obj,OBJ) ;
      --nOBJused ;
    }
  else
    messcrash("Bad obj in OBJfree") ;
}

void OBJstatus (int *used, int *alloc)
{ *used = nOBJused ; *alloc = nOBJalloc ;
}

static KEY localModel (KEY key)
{
#if !defined(NEW_MODELS)
  return class (key) ;
#else
  if (class(key) == _VModel)
    return key ;
  else
    return pickList[class(key)].model ;
#endif
}

void OBJShutDown (void)
{
  messfree (OBJhandle) ;
}

/**************************************************************/
/**************************************************************/

KEY bsKey (OBJ obj)
{
  if (!obj || obj->magic != OBJ_MAGIC)
    return 0 ;
  return obj->key ;
}

char *bsName (OBJ obj)
{
  if (!obj || obj->magic != OBJ_MAGIC)
    return "(NULL-OBJ)" ;
  return name(obj->key) ;
}

/**************************************************************/
/**************************************************************/
/*             Object oriented operations                     */
/********************   Creations *****************************/
/**************************************************************/

                 /* cannot fail */
static OBJ objDecorate (KEY key, CACHE x, BS bs)
{
 OBJ obj = OBJalloc() ;
 KEY model = localModel(key) ;

 obj->key = key ;
 obj->magic = OBJ_MAGIC ;
 obj->x = x ;
 obj->localModel = model ;
 obj->curr = obj->root = obj->localRoot = bs ;

 if (!isModel(key))
   { CHECK_MODEL(model) ;
     obj->modCurr = array(tabModel,KEYKEY(model),BS) ;
   }

 obj->xref = 0 ;
 obj->flag = 0 ;

 sessionAutoSave (-9, 0) ; /* register date of last event */
 /* this is a nice place to register acedb activity,
  * because it belongs to libace and is raised by
  * tace and by xace
  */

 return obj ;
}

/*************************************/

OBJ bsCreate (KEY key)
{
  void *bs ;
  CACHE x ;

  if (!key) return 0 ;
  if (pickType(key) != 'B')
    { chronoReturn(); return 0 ; }

  chrono("bsCreate") ;

  key = lexAliasOf(key) ;

  if (externalServer) externalServer (key, 0, 0, FALSE) ;

  if (iskey(key) != 2)
    { chronoReturn(); return 0 ; }

  nCreated++ ;
  x = cacheCreate(key, &bs, FALSE) ;   /* cannot fail since iskey == 2 */
  chronoReturn();
  return
    objDecorate(key, x, bs) ;
}

/**************************************************************/
          /* You can freely edit this copy, but can t save it */
OBJ bsCreateCopy (KEY key)
{
  void *bs ;
  CACHE x ;

  chrono("bsCreate") ;

  key = lexAliasOf(key) ;
  if (pickType(key) != 'B')
    { chronoReturn(); return 0 ; }

  if (externalServer) externalServer (key, 0, 0, FALSE) ;

  if (iskey(key) != 2)
    { chronoReturn(); return 0 ; }

  nCreated++ ;
  x = cacheCreate(key, &bs, TRUE) ;   /* cannot fail since iskey == 2 */
  chronoReturn();
  return
    objDecorate(key, x, bs) ;
}

/**************************************************************/
#ifndef READONLY

extern BOOL READING_MODELS;

OBJ bsUpdate (KEY key)
{
  void *bs ;
  CACHE x ;

  if (!key)
    return 0 ;

  key = lexAliasOf(key) ;
  if (pickType(key) != 'B')
    messcrash ("bsUpdate called on non B key %s:%s",
	       className(key), name(key)) ;

  if (isModel(key) && !READING_MODELS)
    messcrash ("Erroneous attempt to update a model object!!!!!! Seek help now!");


  if (externalServer)    /* get a fresh copy, unless locally touched */
    {
      if ((lexGetStatus(key) & SERVERSTATUS) &&
	  !(lexGetStatus(key) & TOUCHSTATUS))
	lexUnsetStatus(key ,  SERVERSTATUS) ;
      externalServer (key, 0, 0, FALSE) ;
    }

  x = cacheUpdate(key, &bs) ;

  if (!x)
    { chronoReturn(); return 0 ; }  /* may be locked */

  if (BSTEST)
    cacheCheck(key, x) ;

  nUpdated++ ;

  chrono("bsUpdate") ;
  chronoReturn();

  return
    objDecorate(key, x, bs) ;
}

/* create or update object refered to by key; throw away existing contents
   (if any) and replace with source */

OBJ bsClone (KEY key, OBJ source)
{
  void *bs ;
  CACHE x ;

  if (class(key) != class(source->key))
    messcrash("Attempt to clone into object of different class in bsClone");

  if (!key) return 0 ;
  key = lexAliasOf(key) ;
  if (pickType(key) != 'B')
    messcrash ("bsUpdate called on non B key %s:%s",
	       className(key), name(key)) ;

  x = cacheUpdate(key, &bs) ;

  if (!x)
    { chronoReturn(); return 0 ; }  /* may be locked */

  cacheClone(source->x, x, &bs);
  ((BS)bs)->key = key;  /* give it it's new identity */

  return
    objDecorate(key, x, bs) ;
}


#endif
/**************************************************************/
/* needed jun 2005 becuase we must transfer the texts tot he neObj->x->handle */

static void bsTransferHandle (AC_HANDLE h, BS bs)
{
  while (bs)
    {
      if (bs->right)
	bsTransferHandle (h, bs->right) ;
      if (bs->bt && bs->bt->cp)
	bs->bt->cp = strnew (bs->bt->cp, h) ;

      bs = bs->down ;
    }  
}


BS bsTransferBranch (OBJ newObj, OBJ obj, BOOL whole)
{
  BS bs1 = obj->curr ;

  /* isolate the bs1 branch */
  if (whole && bs1->right) 
    { 
      bs1 = bs1->right ;
      if (bs1->down && !bs1->right) bs1 = bs1->down ;
    }
  else if (bs1->down) /* mieg, dec 98  else needed*/
    {
      bsTreePrune (bs1->down) ; 
      bs1->down = 0 ;
    }
  
  if (bs1->up)
    {
      if (bs1->up->down == bs1) 
	bs1->up->down = 0 ;
      else if (bs1->up->right == bs1) 
	bs1->up->right = 0 ;
    }
  /* transfer owner ship */
  
  bsTransferHandle (cacheStoreHandle (newObj->x), bs1) ;  
  
  return bs1 ;
}


/**************************************************************/
#ifndef JUNK
                 /* cannot fail */
OBJ bsGhost (int classe)
{
  void *bs ;
  CACHE x ;
  KEY key = 0 ;

  if (pickList[classe & 255 ].type != 'B')
    messcrash ("bsGhost called on non B key %s",
	       classe) ;

  x = cacheUpdate(key, &bs) ;

  if (!x)
    return 0 ;  /* may be locked */

  chrono("bsGhost") ;
  chronoReturn();
  return
    objDecorate(key, x, bs) ;
}

#endif
/**************************************************************/
      /* Moves to the model of the subtree, sets objCurr to bsHere  */
BOOL bsPushObj (OBJ obj)
{
  int model ;
  BS bsm = obj->modCurr ;

  if (!bsm || !bsm->right || !bsIsType (bsm->right))
    return FALSE ;

  model = localModel(bsm->right->key) ;

  obj->localModel = model ;
  obj->localRoot = obj->curr ;
  CHECK_MODEL(model) ;
  obj->modCurr = array(tabModel,KEYKEY(model),BS) ;

  return TRUE ;
}

/**************************************************************/
/**************  Destructions   *******************************/
/**************************************************************/

void  bsDoDestroy(OBJ obj)
{
  if (!obj)
    return ;
  if (obj->magic != OBJ_MAGIC)
    messcrash ("bsDoDestroy received  a non magic B obj (%s)",name(obj->key)) ;

  ac_free (obj->x->assHandle) ; obj->x->tagAss = 0 ;
  cacheDestroy (obj->x) ; /*bsTreePrune (obj->root) ; could be used for dummys */

  OBJfree(obj) ;
  nDestroyed++ ;
}

/**************************************************************/

#ifndef READONLY
static void doXref (OBJ obj) ;

void  bsDoKill (OBJ obj)
{
  KEY key ;

  if (!obj)
    return ;
  key = obj->key ;
  if (externalSaver)   /* dump old version */
    externalSaver (key, 4) ;

  bsGoto (obj, 0) ;
  ac_free (obj->x->assHandle) ; obj->x->tagAss = 0 ;

  while (bsGetKeyTags (obj, _bsRight, 0))
    bsRemove (obj) ;		 /* does bsGoto(obj,0) after remove */
  obj->flag |= OBJFLAG_TOUCHED ; /* in case nothing removed */

  bsSave (obj) ;	/* safer this way - e.g. checks subclasses */

  nSaved-- ;
  nKilled++ ;
  if (externalSaver)   /* release kill lock */
    externalSaver (key, 5) ;
}

/**************************************************************/

BOOL bsIsDirty (OBJ obj)
{
  return (obj->flag & OBJFLAG_TOUCHED) ? TRUE : FALSE ;
}

/**************************************************************/

void bsDoSave(OBJ obj)
{
  if (!obj)
    return ;
  nSaved++ ;

  if (!(obj->flag & OBJFLAG_TOUCHED))
    {
#ifdef SAVE_DEBUG
      printf("bsDoSave2 %s: end - untouched obj\n", nameWithClass(obj->key)) ;
#endif

      bsDestroy (obj) ;

      nDestroyed-- ;
      return ;
    }

  if (obj->magic != OBJ_MAGIC)
    messcrash ("bsDoSave received  a non magic B obj (%s)",name(obj->key)) ;
  if (!obj->root)
    messcrash ("bsDoSave says : Jean owes a beer to Richard") ;

  sessionAutoSave (-8, 0) ; /* register date of last modif */

  if (pickList[class(obj->key)].constraints)
    {
    }
  if (obj->key != obj->root->key)
    messcrash ("key mismatch in objStore %s != %s",
	       name(obj->key), name(obj->root->key)) ;

  if (externalSaver)   /* dump old version */
    externalSaver (obj->key, 1) ;

  if (BSTEST)
    cacheCheck(obj->key, obj->x) ;
  if (obj->xref)
    doXref (obj) ;
  if (BSTEST)
    cacheCheck(obj->key, obj->x) ;

				/* set subclass lex bits */
  bIndexObject(obj->key, obj) ;

  ac_free (obj->x->assHandle) ; obj->x->tagAss = 0 ;
  cacheSave (obj->x) ;		/* must come last */

  lexSetStatus (obj->key, TOUCHSTATUS) ;
  if (obj->root->right)
    { if (lexGetStatus (obj->key) & EMPTYSTATUS)
	{ lexUnsetStatus (obj->key, EMPTYSTATUS) ;
	  lexAlphaMark (class (obj->key)) ;
	  lexmark(class(obj->key));
	}
    }
  else				/* setSession might have emptied it */
    { if (!(lexGetStatus (obj->key) & EMPTYSTATUS))
	{ lexSetStatus (obj->key, EMPTYSTATUS) ;
	  lexAlphaMark (class (obj->key)) ;
	  lexmark(class(obj->key));
	}
    }

  /* All cache and lex operations must be complete by this
     point, so a new obj can be made on the same key.  
     We rely on OBJfree() just releasing memory.  
  */
  if (externalSaver)   /* dump new version and acediff it */
    externalSaver (obj->key, 2) ;

  OBJfree(obj) ;

#ifdef SAVE_DEBUG
  printf("bsDoSave2 %s: end - touched obj\n", nameWithClass(obj->key)) ;
#endif
  return ;
}

#endif

/**************************************************************/
/*                  Local surgery                             */
/**************************************************************/
/**************** first the get package ***********************/

/* if tag present, returns TRUE and moves curr, else nothing */

BOOL bsRemoveTag (OBJ obj, KEY target)
{
  if (bsFindTag (obj, target))
    { 
      /* it is impossible to know how many tags get removed
       * right of the tag we wish to remove 
       */
      ac_free (obj->x->assHandle) ; obj->x->tagAss = 0 ; 
      bsRemove (obj) ; 
      return TRUE ;
    }
  return FALSE ;
}

/* if tag present, returns TRUE and moves curr, else nothing */

BOOL bsFindTag (OBJ obj, KEY target)
{
  BS    bs, bsm ;

  if (!obj)
    return FALSE ;
  if (obj->magic != OBJ_MAGIC)
    messcrash ("bsFindTag called with null or bad object");

  switch (target)
    {
    case _bsRight :		/* used in queryexe */
      bs = obj->curr->right ;
      while (bs && bsIsComment(bs))
	bs = bs->down ;
      bsm = bsModelRight(obj->modCurr) ;
      if (!bs || !bsm)		/* off the end of the model */
	return FALSE ;
      bsModelMatch (bs, &bsm) ;
      obj->curr = bs ;
      obj->modCurr = bsm ;
      return TRUE ;
    case _bsDown :
      messcrash ("bsFindTag called on bsDown, the result cannot be guaranteed") ;
    case _bsHere:
      return TRUE ;
    default:
      return findTag (obj, target, FALSE) ;      /* FALSE -> do not grow path  */
    }
}

/*******************************************************************/
/*******************************************************************/
/* simple routines, avoid opening the object in the calling routine */

KEY keyGetKey(KEY key, KEY tag)
{
  OBJ obj = 0 ;
  KEY val = 0 ;

  if (!tag)
    messcrash ("keyGetKey called with NULL tag, may be this tag is not initialised") ;
  if (key && bIndexFind(key, tag) &&
      (obj = bsCreate(key)))
    { bsGetKeyTags (obj, tag, &val) ;
      bsDestroy (obj) ;
    }
  return val ;
}

/*******************************************************************/
/* same thing adapted to tag2 system, in particular SMAP */

KEY keyGetTag2Key(KEY key, KEY tag)
{
  OBJ obj = 0 ;
  KEY val = 0, tag2 ;
  static BSMARK mark ;

  static Array units = 0 ;

  units = arrayReCreate (units, 32, BSunit) ;
  if (key && bIndexFind(key, tag) &&
      (obj = bsCreate(key)))
    {
      if (bsGetKeyTags (obj, tag, &tag2))
	do  /* do while construct */
	  {
	    mark = bsMark (obj, mark) ;
	    if (bsGetKey (obj, _bsRight, &val))
	      break ;
	    bsGoto (obj, mark) ;
	  } while (bsGetKeyTags (obj, _bsDown, &tag2)) ;
      bsDestroy (obj) ;
    }
  return val ;
}

/*******************************************************************/

int keyGetInt (KEY key, KEY tag)
{
  OBJ obj = 0 ;
  int val = 0 ;

  if (!tag)
    messcrash ("keyGetKey called with NULL tag, may be this tag is not initialised") ;
  if (key && bIndexFind(key, tag) &&
      (obj = bsCreate(key)))
    { bsGetData (obj, tag, _Int, &val) ;
      bsDestroy (obj) ;
    }
  return val ;
}

/*******************************************************************/

float keyGetFloat (KEY key, KEY tag)
{
  OBJ obj = 0 ;
  int val = 0 ;

  if (!tag)
    messcrash ("keyGetKey called with NULL tag, may be this tag is not initialised") ;
  if (key && bIndexFind(key, tag) &&
      (obj = bsCreate(key)))
    { bsGetData (obj, tag, _Float, &val) ;
      bsDestroy (obj) ;
    }
  return val ;
}

/*******************************************************************/

const char *keyGetText (KEY key, KEY tag) /* not paralelisable */
{
  OBJ obj = 0 ;
  static char buf[1024] ;
  char *cp ;

  buf[0] = 0 ;
  if (!tag)
    messcrash ("keyGetKey called with NULL tag, may be this tag is not initialised") ;
  if (key && bIndexFind(key, tag) &&
      (obj = bsCreate(key)))
    {
      if (bsGetData (obj, tag, _Text, &cp))
	strncpy (buf, cp, 1023) ;
      buf[1023] = 0 ;
      bsDestroy (obj) ;
    }
  return buf ; /* not paralelisable */}

/*******************************************************************/


KEY keyFindTag (KEY key, KEY tag)
{
  OBJ obj = 0 ;
  BOOL val = FALSE ;

  if (key && tag >= _Date && class(tag) == 0)
    switch (bIndexFind(key, tag))
      {
      case 0:
	return FALSE ;
      case 1:
	if ((obj = bsCreate(key)))
	  {
	    val = bsFindTag (obj, tag) ;
	    bsDestroy (obj) ;
	  }
	break ;
      case 2:
	return TRUE ;
      }
  return val ;
}

/*********************************************************************/

int keySetAddTag (KEYSET ks, KEY tag)
{
  OBJ obj = 0 ;
  int ii, nn = 0 ;
  KEY key ;
  
  for (ii = 0 ; ks && ii < keySetMax (ks) ; ii++)
    {
      key = keySet (ks, ii) ;
      if (bIndexFind(key, tag) < 2 && 
	  (obj = bsUpdate (key)))
	{
	  bsAddTag (obj, tag) ;
	  if (bsFindTag (obj, tag))
	    nn++ ;
	  bsSave (obj) ;
	}
    }
  return nn ;
} /* keySetAddTag */

/*********************************************************************/

int keySetRemoveTag (KEYSET ks, KEY tag)
{
  OBJ obj = 0 ;
  int ii, nn = 0 ;
  KEY key ;
  
  for (ii = 0 ; ks && ii < keySetMax (ks) ; ii++)
    {
      key = keySet (ks, ii) ;
      if (bIndexFind(key, tag) > 0 && 
	  (obj = bsUpdate (key)))
	{
	  if (bsRemoveTag (obj, tag))
	    nn++ ;
	  bsSave (obj) ;
	}
    }
  return nn ;
} /* keySetRemoveTag */

/*******************************************************************/
/*******************************************************************/
        /* basic routine to read from an object
           target = tag or _bsRight or _bsDown or _bsHere (relative to curr)
	   ignores comments
           if there is something there,
             fills found, updates curr and returns TRUE
           else
             returns FALSE
        */

static BOOL bsGetKey2 (OBJ obj, KEY target, KEY *found, BOOL alsoGetTags)
{
  BS    bs, bsm ;

  bs = obj->curr ;
  bsm = obj->modCurr ;

  switch (target)
    {
    case _bsRight :
      bs = bs->right ;
      bsm = bsModelRight (bsm) ;
      break ;
    case _bsDown :
      bs = bs->down ;
      break ;
    case _bsHere:
      break ;
    default:
      if (!findTag (obj, target, FALSE))
        return FALSE ;
      bs = obj->curr->right ;   /* go right from tag */
      bsm = bsModelRight (obj->modCurr) ;
    }

  if (!bsm)		/* off the end of the model */
    return FALSE ;

  while (bs)
    { if (!bsIsComment(bs) && (class(bs->key) || alsoGetTags))
        { if (found)
            *found = bs->key ;
	  bsModelMatch (bs, &bsm) ;
          obj->curr = bs ;
	  obj->modCurr = bsm ;
          return TRUE ;
        }
      bs = bs->down ;     /* iterate downwards, ignoring user comments */
    }

   return FALSE ;
}

BOOL bsGetKey (OBJ obj, KEY target, KEY *found)
{
  return
    bsGetKey2(obj,target,found,FALSE) ;
}

   /* Also get tags, used by query package */
BOOL bsGetKeyTags (OBJ obj, KEY target, KEY *found)
{
  return
    bsGetKey2(obj,target,found,TRUE) ;
}

KEY bsParentKey(OBJ obj)
{ int c ;
  BS bs = obj->curr ;
  while(bs)
    { c = class(bs->key) ;
      if (c && c != _VComment && c != _VUserSession)
	return bs->key ;
      bs = bs->up ;
    }
  return
    obj->key ;
}
				/* for tag2 systems */

BOOL bsFindKey2 (OBJ obj, KEY tag, KEY key, KEY *ptag2)
{
  static BSMARK mark = 0 ;

  if (!bsFindTag (obj, tag))
    return FALSE ;

  if (bsGetKeyTags (obj, _bsRight, ptag2)) do
    { mark = bsMark (obj, mark) ;
      if (bsFindKey (obj, _bsRight, key))
	return TRUE ;
      bsGoto (obj, mark) ;
    } while (bsGetKeyTags (obj, _bsDown, ptag2)) ;

  return FALSE ;
}

/*******************************************************************/

        /* matches to bs->key type and returns TRUE if found,
             else goes to bottom of column and returns FALSE
           can call with x=0 to test if match found and update curr
        */

BOOL bsGetData (OBJ obj, KEY target, KEY type, void *x)
     /* WARNING:
      * if called with type==_Text it is currently not threadsafe
      * as another process could reset the pointer before the
      * calling process could secure the text using strnew or the like.
      * The solution would be to have a new function
      * bsGetTextData, which is given a AC_HANDLE and the pointer
      * to a string (i.e. a char **) so it can do the strnew
      * inside that function */
{
  KEY k ;

  if (!bsGetKey2 (obj,target,&k, TRUE))
    return FALSE ;

  while (k != type)
    if (!bsGetKey2 (obj,_bsDown,&k, TRUE))
      return FALSE ;

  if (!x)
    return TRUE ;

  switch(type)
    {
    case _Int :
      *(int*)x = obj->curr->n.i ;
      break ;
    case _Float :
      *(float*)x = obj->curr->n.f ;
      break ;
    case _DateType :
      *(mytime_t*)x = obj->curr->n.time ;
      break ;
    default:
      if(type<=_LastC)
        *(char**)x = bsText(obj->curr) ;
      else
        messcrash("Wrong type in bsGetData") ;
      break ;
    }

  return TRUE ;
}

/* Note: to get a vector of, say, int's into Array a, do something like:
   arrayMax(a) = 0 ; i = 0 ;
   if (bsFindTag (obj,tag))
     while (bsGetData (obj,_bsRight,_Int,&ix))
       array(a,i++,int) = ix ;
*/

/***********************************************************/
/* here we should use hprintf */
BOOL bsGetText (OBJ obj, KEY target, char **cpp)
{
  static char buffer[64] ;
  KEY k ;
  char timeBuf[25] ;

  if (!bsGetKey2 (obj,target,&k,TRUE))
    return FALSE ;

  if (!cpp)
    return TRUE ;

  switch (k)
    {
    case _Int :
      sprintf (buffer, "%d", obj->curr->n.i) ;
      *cpp = buffer ;
      break ;
    case _Float :
      sprintf (buffer, "%1.7g", obj->curr->n.f) ;
      *cpp = buffer ;
      break ;
    case _DateType :
      *cpp = timeShow (obj->curr->n.time, timeBuf, 25) ;
      break ;
    default:
      if (k <= _LastC)
        *cpp = bsText(obj->curr) ;
      else
	*cpp = name(k) ;
      break ;
    }

  return TRUE ;
}

/***********************************************************/

    /* bsGetArray(), bsFlatten,
       flattens a subtree into an array of bsUnits
       provided by the user. The array contains a multiple of n
       BSunits, each containing an int, float, char*
       or KEY according to the model, which the user is presumed
       to know.  Missing data are represented by 0.
    */

BOOL bsGetArray (OBJ obj, KEY target, Array a, int width)
{
  if (!arrayExists(a))
    messcrash("bsgetArray received a non existing array %d", a) ;
  if (a->size != sizeof(BSunit))
    messcrash("bsgetArray received an array of wrong type (should be BSunit)") ;
  arrayMax(a) = 0 ;

  return
    bsFindTag(obj, target) && bsFlatten (obj, width, a)
      ? TRUE : FALSE ;
}

/************/

BOOL bsFlatten (OBJ obj, int n, Array a)
{
  int i, j, m = 0 ;
  KEY direction = _bsRight ;
  BS curr, bs ;
  BSunit *u ;
  static BSunit zero ;		/* rely on this being 0 */
  static Array currStack = 0, unit = 0 ;
  BOOL gotit = TRUE ;		/* don't want an empty list */
  BOOL added ;

  if (!arrayExists(a))
    messcrash("bsFlatten received a non existing array %d", a) ;
  if (a->size != sizeof(BSunit))
    messcrash("bsFlatten received an array of wrong type (should be BSunit)") ;
  arrayMax(a) = 0 ;

  if (!obj)
    return FALSE ;
  if (obj->magic != OBJ_MAGIC)
    messcrash ("bsFlatten got a bad obj") ;
  curr = obj->curr ;
  if (!curr->right)
    { arrayMax (a) = 0 ;
      return FALSE ;
    }

  currStack = arrayReCreate (currStack, n, BS) ; arrayMax(currStack) = n ;
  unit = arrayReCreate (unit, n, BSunit) ; arrayMax(unit) = n ;

  for (i = 0 ; i >= 0 ;)
    { added = FALSE ;
      bs = (direction == _bsRight) ? curr->right : curr->down ;
      while (bs)
	{ if (!bsIsComment(bs)) /* add this node to the unit list */
	    { curr = bs ;
	      arr(currStack,i,BS) = curr ;
	      if (class (curr->key) || curr->key > _LastN)
		arr(unit,i,BSunit).k = curr->key ;
	      else
		switch (curr->key)
		  {
		  case _Int:
		    arr(unit,i,BSunit).i = curr->n.i ;
		    break ;
		  case _Float:
		    arr(unit,i,BSunit).f = curr->n.f ;
		    break ;
		  case _DateType:
		    arr(unit,i,BSunit).time = curr->n.time ;
		    break ;
		  default:
		    if (curr->key <= _LastC)
		      arr(unit,i,BSunit).s = bsText(curr) ? 
			strnew (bsText(curr), cacheStoreHandle (obj->x)) : 0 ;
		    else
		      messcrash ("Bad type in bsFlatten") ;
		    break ;
		  }
	      ++i ;
	      gotit = FALSE ;
	      direction = _bsRight ;
	      added = TRUE ;
	      break ;
	    }
	  bs = bs->down ;	/* iterate downwards, ignoring user comments */
	}
      if (!added || i >= n)
	{ if (!gotit)
	    { array(a,n*(++m)-1,BSunit) = zero ; /* set max */
	      u = arrp(a,n*(m-1),BSunit) ;
	      for (j = 0 ; j < i ;)
		*u++ = arr(unit,j++,BSunit) ;
	      while (j++ < n)
		(*u++).k = 0 ;
	      --i ;
	      gotit = TRUE ;
	    }
	  if (direction == _bsDown)
	    { --i ;
	      if (i >= 0)
		curr = arr(currStack,i,BS) ;
	    }
	  direction = _bsDown ;
	}
    }

  return arrayMax(a) ? TRUE : FALSE ;
}

/***********************************************************/

static KEYSET viKeySet = 0 ;

static void bsViGet (BS bs)
{
  keySet(viKeySet, keySetMax(viKeySet)) = bs->key ;

  if (bs->right)
    bsViGet (bs->right) ;
  if (bs->down)
    bsViGet (bs->down) ;
}

/*********************************/

     /* Returns a keySet containing all the  Keys
	that can be found in bsCreate(key)
     */

KEYSET bsKeySet (KEY key)
{
  OBJ  obj ;

  if (!key || pickType(key) != 'B' || !(obj  = bsCreate(key)))
    return 0 ;

  viKeySet = keySetCreate() ;
  bsViGet (obj->root) ;
  bsDestroy (obj) ;

  keySetSort (viKeySet) ;
  keySetCompress (viKeySet) ;
  return viKeySet ;
}

/*********************************/

KEYSET bsTagsInClass (int table)
{
  KEY model ;
  int j = 0 ;
  char *vt = 0;
  KEYSET result = keySetCreate() ;
  Associator ass ;

  if (pickList[table].type != 'B')
    return result ;

  model = localModel(KEYMAKE(table,0)) ;
  CHECK_MODEL(model) ;
  ass = array(tabAssPath,KEYKEY(model),Associator) ;
  while (assNext(ass, &vt, 0))
    keySet(result,j++) = assInt (vt) ;

  keySetSort (result) ;
  keySetCompress (result) ;
  return result ;
}

/*********************************/

BOOL bsIsClassInClass (int table, int maClasse)
{
  KEY model = 0 ;
  int i = 0 ;
  KEYSET ks = 0 ;
  KEY *kp ;
  BOOL found = FALSE ;

  if (pickList[table].type != 'B')
    return FALSE ;

  model = localModel(KEYMAKE(table,0)) ;
#if !defined(NEW_MODELS)
  model = KEYMAKE (model, 0) ;
#else
  messcrash (" bsIsClassInClass not coded for new_models, sorry") ;
#endif
  ks =  bsKeySet(model) ;
  if (!ks) return FALSE ;
  kp = arrp(ks, 0, KEY) - 1 ;
  i = keySetMax(ks) ;
  while (kp++, i--)
    if (class(*kp) == maClasse)
      { found = TRUE ; break ; }
  keySetDestroy (ks) ;
  return found ;
}

/*********************************/

BOOL bsIsTagInObj (OBJ obj, KEY key, KEY tag)
{
  char *vt ;
  KEY model ;

  if (obj)
    model = obj->localModel ;
  else if (pickType(key) == 'B')
    model = localModel (key) ;
  else
    return FALSE ;

  if (!model)
    return FALSE ;

  CHECK_MODEL(model) ;

  if (tag == _bsHere ||
      tag == _bsRight ||
      tag == _bsDown)
    return TRUE ;

  vt = assVoid(tag) ;
  return assFind (array(tabAssPath,KEYKEY(model),Associator), vt, 0) ;
}

BOOL bsIsTagInClass (int class, KEY tag)
{
  KEY key = KEYMAKE (class, 1) ;

  if (pickType (key) != 'B')
    return FALSE ;
  return
    bsIsTagInObj (0, key, tag) ;
}

/*********************************/
/* suzi, these next 2 routines are for the query builder */
/*********************************/

static void addKeyTag (KEYSET tagks, KEYSET keyks, BS bs, BS prev)
{
    if (iskey (bs->key) == 2
        && class(bs->key) != _VText
        && prev
        && iskey (prev->key) != 2) {
        keySet (tagks, tagks->max) = prev->key ;
        keySet (keyks, keyks->max) = bs->key ;
    }
    if (bs->right && bsModelRight(bs) != bs)
        addKeyTag (tagks, keyks, bsModelRight(bs), bs) ;
    if (bs->down)
        addKeyTag (tagks, keyks, bs->down, 0) ;
}

/*********************************/

void bsKeyTagsInClass (int classe, KEYSET *tagks, KEYSET *keyks)
{
  KEYSET tagresult = keySetCreate() ;
  KEYSET keyresult = keySetCreate() ;
  KEY model = localModel (KEYMAKE(classe,0)) ;

  if (pickList[classe].type == 'B')
    {
      CHECK_MODEL(model) ;
      addKeyTag (tagresult, keyresult,
		 (array(tabModel,classe,BS))->right, 0) ;
    }
  *tagks = tagresult ;
  *keyks = keyresult ;
}

/***************************************************************/
/***************************************************************/
/************************ Add package **************************/


#ifndef READONLY

/* main function of bsAddTag is to leave curr in the right place.
   if necessary adds tag according to model to get there.
   returns TRUE if it adds stuff.
*/

BOOL bsAddTag (OBJ obj, KEY tag)
{
  if (!obj)
    return FALSE ;
  if (obj->magic != OBJ_MAGIC)
    messcrash ("bsAddTag received  a non magic B obj (%s)",name(obj->key)) ;

  if (!isCacheLocked(obj->x))
    { messout ("bsAddTag fails because %s is not locked\n",
               name (obj->key)) ;
      return FALSE ;
    }

  return findTag (obj,tag,TRUE) ;       /* i.e. new tags made */
}

/***************************************************************/

KEY bsGetTimeStamp (OBJ obj)
{
  if (obj && obj->curr)
    return obj->curr->timeStamp ;
  else
    return 0 ;
}

#ifdef JUNK
mytime_t bsGetNodeTime(OBJ obj)
{ if (obj && obj->curr && obj->curr->timeStamp )
    return sessionUserSession2Time(obj->curr->timeStamp) ;
  else
    return 0 ;
}
#endif /* JUNK */

void bsSetTimeStamp (OBJ obj, KEY stamp)
{
  obj->curr->timeStamp = stamp ;
  obj->flag |= OBJFLAG_TOUCHED ;
  cacheMark (obj->x) ;
}

/******************************************************************/

/* basic routine to add to an object.
   searches down the column for key new.
   if there is a key new there it moves curr to it and returns FALSE
   else it makes a BS with key new and puts curr on it and returns TRUE
   handles user comments properly (now, we hope)
   if target is _bsHere, both this and bsAddData change whatever you
     were on.  Caveat user.
*/

BOOL bsAddKey (OBJ obj, KEY target, KEY new)
{
  BS    bs, bsNew, bsUp = 0,bs1 = 0, bsm ;
  BOOL isNewAnAlias ;
  KEY target0 = 0 ;
  Associator ass = 0 ;
  int nDown = 0 ;
  /* BOOL debug = FALSE ; */

  if (!obj)
    return FALSE ;
  if (obj->magic != OBJ_MAGIC)
    messcrash ("bsAddKey received a non magic B obj (%s)", name(obj->key)) ;

  /* if (class (new) == _VSequence) debug = TRUE ; */

  switch (target)
    {
    case _bsHere:		/* do everything here */
      if (new == obj->curr->key)
	return FALSE ;
      if (!class(new))
	messcrash ("bsAddKey(bsHere) is forbidden for a tag") ;
      if (class(new) != class(obj->curr->key))
	messcrash ("type mismatch in bsAddKey(bsHere)") ;
      if (!isXrefing)
	messcrash ("bsAddKey(bsHere) while X referencing") ;
      if (obj->curr == obj->root)
	messcrash ("bsAddKey(bsHere) at root of %s", name(obj->key)) ;
      bsm = obj->modCurr ;
      if (KEYKEY(bsm->n.key) && class(bsm->n.key) != _VConstraint)
/* used to be KEYKEY(bsm->n.key) &&
   (obj->root == obj->localRoot || bsm->n.key == _Quoted_in))
   to prevent XREF from within subobjects, but why do this?
*/
	{ if (!obj->xref) obj->xref = stackCreate (64) ;
	  push (obj->xref, obj->curr->key, KEY) ;
	  push (obj->xref, KEYKEY(bsm->n.key) | DELETE_BIT, KEY) ;
	  push (obj->xref, new, KEY) ;
	  push (obj->xref, KEYKEY(bsm->n.key), KEY) ;
	}
      obj->curr->key = new ;
      obj->curr->timeStamp = sessionUserKey() ;
      obj->flag |= OBJFLAG_TOUCHED ;
      cacheMark (obj->x) ;
      return TRUE ;
    case _bsDown:
      bsUp = obj->curr ;  /* node on which we will add the hook */
      bs = bsUp->down ;
      bsm = obj->modCurr ;
      break ;
    default:
      target0 = target ;
      if (!bsAddTag (obj,target))
	return FALSE ;
      target = _bsRight ;
    case _bsRight:	/* note fall through from default */
      bsUp = obj->curr ;	/* must come after bsAddTag() */
      if (bsUp->key && ! class (bsUp->key) && obj->root == obj->localRoot)
	target0 = bsUp->key ;
      bs = bsUp->right ;
      bsm = bsModelRight (obj->modCurr) ;
      if (!bsm)
	{ messerror ("Trying to add %s off end of model for %s",
		     name(new), name(obj->key)) ;
	  return FALSE ;
	}
    }

/* from here on target must be _bsDown or _bsRight */
  isNewAnAlias = lexIsAlias (new) ;

  if (bs && target0 && obj->x->tagAss &&
      assFind (obj->x->tagAss, assVoid (target0), &ass) &&
      ass
      )
    {
      bs1 = 0 ;
      if (assFind (ass, assVoid (new), &bs1) || (isNewAnAlias && assFind (ass, assVoid (new), &bs1)))
	{
	  /* if (debug) fprintf (stderr, "## bssubs recovered %s bs=%p bs->up=%p\n", name(bs1->key), bs1, bs1->up) ;  */
	  bsModelMatch (bs, &bsm) ;
	  obj->curr = bs1 ;
	  obj->modCurr = bsm ;
	  return FALSE ;
	}
      else
	if (assFind (ass, assVoid (-2), &bs1))
	  { 
	    bs = bs1 ;             /* we are immediatly positionned correctly */
	    nDown = MAXDOWN + 1 ;  /* we pretend we have to loop to force registering the new value */
	    /* if (debug) fprintf (stderr, "## bssubs located %s bs=%p bs->up=%p\n", name(bs1->key), bs, bs->up) ;  */
	  }
    }
	

  while (bs)			/* go to bottom of column */
    {
      if (bs->key == new || 	/* check that new isn't there already */
	  (isNewAnAlias && lexAliasOf (bs->key) == new )  /* mieg 12-11-2004 */
	  ) 
        { 
	  bsModelMatch (bs, &bsm) ;
	  obj->curr = bs ;
	  obj->modCurr = bsm ;
	  return FALSE ;
	}
      bsUp = bs ;
      bs = bs->down ;
      target = _bsDown ;
      nDown++ ;
    }

  bsNew = BSalloc () ;
  bsNew->timeStamp = sessionUserKey() ;
  bsNew->key = new ;
  if (!bsModelMatch (bsNew, &bsm))
    { messout ("Type of %s does not check in B object \"%s\"",
	       name(new), name(obj->key)) ;
      BSfree (bsNew) ;
      return FALSE ;
    }

  if (bsm->n.key && class(bsm->n.key) == _VConstraint &&
      !checkKey(new, bsm->n.key))
    { messout ("Key %s does not check constraint in B object \"%s\"",
	       name(new), name(obj->key)) ;
      BSfree (bsNew) ;
      return FALSE ;
    }

  if (isXrefing && KEYKEY(bsm->n.key) &&
      (!XREF_DISABLED_FOR_UPDATE || pickList[class(new)].Xref))
    { if (!obj->xref) obj->xref = stackCreate (64) ;
      push (obj->xref, new, KEY) ;
      push (obj->xref, KEYKEY(bsm->n.key), KEY) ;
    }
	/* if unique and something there, delete noncomment entries */
  if ((bsm->n.key & UNIQUE_BIT) && target == _bsDown)
    {				/* start again at the top */
      nDown = 0 ;
      for (bs = bsUp ; bs->up->down == bs ; bs = bs->up) ;
      bsUp = bs->up ; target = _bsRight ;
      while (bs)
	if (bsIsComment (bs))
	  { bsUp = bs ; bs = bs->down ; target = _bsDown ; }
	else
	  { if (bsUp->right == bs)
	      bsUp->right = bs->down ;
	    else
	      bsUp->down = bs->down ;
	    if (bs->down)
	      bs->down->up = bsUp ;
	    bs->down = 0 ;
/* RD 990518 important bug fix: we used to not xref deletions caused by replacing
   a UNIQUE entry during XREF.  This was wrong, leading to database inconsistency.
   Forcing isXrefing to TRUE below fixes this.
*/
	    if (isXrefing)
	      xrefPrune (obj, bs, bsm) ;
	    else
	      { isXrefing = TRUE ;
	        xrefPrune (obj, bs, bsm) ;
		isXrefing = FALSE ;
	      }
	    bs = (target == _bsDown) ? bsUp->down : bsUp->right ;
	  }
    }

  bsNew->up = bsUp ;
  if (target == _bsRight)	/* hook right */
    bsUp->right = bsNew ;
  else				/* hook down */
    bsUp->down = bsNew ;

  obj->curr = bsNew ;
  obj->modCurr = bsm ;
  obj->flag |= OBJFLAG_TOUCHED ;
  cacheMark (obj->x) ;
  
  if (!bs && target0 && obj->x->tagAss)
    assFind (obj->x->tagAss, assVoid (target0), &ass) ;

  if ((ass || nDown > MAXDOWN) &&
      (target0 == _Quoted_in || target0 > _LastSystemTag) && 
      class (obj->key) >  _VMainClasses &&
      obj->root == obj->localRoot
      ) 
    {
      if (! obj->x->assHandle)
	obj->x->assHandle = ac_new_handle () ;
      if (! obj->x->tagAss)
	obj->x->tagAss = assHandleCreate (obj->x->assHandle) ;
      if (1 && ! ass)
	{
	  if (! assFind (obj->x->tagAss, assVoid (target0), &ass))
	    {
	      ass = assHandleCreate (obj->x->assHandle) ;
	      assInsert (obj->x->tagAss, assVoid (target0), ass) ;
	      bs1 = bsNew ;
	      while (bs1->up && bs1->up->down == bs1)
		{
		  bs1 = bs1->up ;
		  assInsert (ass, assVoid (bs1->key), bs1) ;
		}
	    }
	}

      if (ass)
	{
	  assInsert (ass, assVoid (bsNew->key), bsNew) ;
	  /* if (debug) fprintf (stderr, "## bssubs registered %s bsNew=%p bsNew->up=%p\n", name(bsNew->key), bsNew, bsNew->up) ;  */
	  assRemove (ass, assVoid (-2)) ;
	  assInsert (ass, assVoid (-2), bsNew) ;  /* last bs of the column */
	  bs1 = 0 ;
	  assFind (ass, assVoid (-2), &bs1) ;
	  /* fprintf (stderr, "## bssubs registered %s\n", name(bs1->key)) ; */
	  bs1 = 0 ;
	}
    }

  return TRUE ;
}

/***************************************************************/

#if 0
void
bs_set_lock(OBJ obj, int i)
{
  obj->curr->locked = i;
  obj->flag |= OBJFLAG_TOUCHED ;
  obj->modCurr = 0 ;
  cacheMark (obj->x) ;
}
#endif

/*
* Add data to the object.
*	obj - the object to add the data to
*	target - the tag name of where to put the data
*		_bsRight and _bsHere are magic values that mean
*		to place the tag relative to the current cursor
*		position in the object.
*	type - the data type of the data we are placing
*	xp - address of the data
*
*/
BOOL bsAddData (OBJ obj, KEY target, KEY type, const void *xp)
         /* NB must call with the ADDRESS of the numbers */
{
  BS bs ;

  if (!obj)
    return FALSE ;
  if (obj->magic != OBJ_MAGIC)
    messcrash ("bsAddData received  a non magic B obj (%s)",name(obj->key)) ;

  while (!bsAddKey (obj,target,type))   /* got a key match */
    {
      bs = obj->curr ;
      if (bs->key != type)
        return FALSE ; /* happens if bsAddKey failed because typeCheck did */
      switch (type)                     /* first check for data match */
        {
        case _Int:
          if (bs->n.i == *(int*)xp)
            return FALSE ;
          break ;
        case _Float:
	  {
	    double zf = *(float*)xp, zf2 ;

	    if (zf == 0 && bs->n.f == 0)  /* test first, frequent case */
	      return FALSE ;
	    zf2 = (zf - bs->n.f) * (zf - bs->n.f) ;
	    if (zf != 0 && bs->n.f != 0 && zf2 < zf * bs->n.f * ACE_FLT_RESOLUTION)  /* quasi equal */
	      return FALSE ;
	  }
          break ;
        case _DateType:
          if (bs->n.time == *(mytime_t*)xp)
            return FALSE ;
          break ;
        default:
          if (type <= _LastC)
            { 
	      char *cr = bsText(bs) ;
	      if (cr)
		{
		  if (!strcmp (cr,(char*)xp))
		    return FALSE ;
		  if (!strcasecmp (cr,(char*)xp)) /* force capitals */
		    { /* do not break because we do not want to clean to the right */
		      strcpy(cr,(char*)xp);	
		      goto done ;
		    }
		}
            }
          else
            messcrash ("Unknown type in bsAddData") ;
        }
      if (target == _bsHere)	/* overwrite */
	break ;

      if (obj->modCurr->n.key & UNIQUE_BIT)
	{ if (bs->right)
	    { xrefPrune (obj, bs->right, obj->modCurr->right) ;
	      bs->right = 0 ;
	    }
	  if (bs->down)
	    { xrefPrune (obj, bs->down, obj->modCurr) ;
	      bs->down = 0 ;
	    }
	  break ;
	}

      target = _bsDown ;
    }

  bs = obj->curr ;

/* NB problem here - maybe you have made a node that you must delete
  if (obj->modCurr->n.key &&
      class(obj->modCurr->n.key) == _VConstraint &&
      !checkData(xp, type, obj->modCurr->n.key))
    return FALSE ;
*/

  switch(type)
    {
    case _Int :
      bs->n.i = *(int*)xp ;
      break ;
    case _Float :
      bs->n.f = *(float*)xp ;
      break ;
    case _DateType :
      bs->n.time = *(mytime_t*)xp ;
      break ;
    default:
      if (type<=_LastC)
        { 
	  if (!bs->bt)
	    bs->bt = BTalloc() ;
	  else if (bs->bt->cp)
	    bs->bt->cp = 0 ; /* do not messfree, it belongs on 
			      *	the handle and the user may be referencing it
			      */
          /* bs->bt->cp = halloc (strlen ((char*)xp) + 1, cacheStoreHandle (obj->x)) ; */
	  bs->bt->cp = strnew ((char*)xp, cacheStoreHandle (obj->x)) ; 
	  {
	    char *cp = (char*)xp ;
	    while (*cp == ' ') cp++ ; /* clean initial spaces */
	    strcpy (bs->bt->cp, cp) ;
	    cp = bs->bt->cp + strlen(bs->bt->cp) - 1 ;
	    while (*cp == ' ' && cp >= bs->bt->cp)
	      *cp = 0 ;  /* clean terminal spaces */
	  }
        }
      else
        messcrash("Wrong type in bsAddData") ;
      break ;
    }

done:
  obj->flag |= OBJFLAG_TOUCHED ;
  bs->timeStamp = sessionUserKey() ;
  cacheMark (obj->x) ;
  return TRUE ;
}

/********************************************/

/* XQX
* Inserts a comment in the object at the current cursor location.
*
* type is the type of comment to insert:
*	C - regular comment, text is the comment to add
*
* Strange - it seems that "new" should have it's tag set to
* _Vcomment somewhere...
*
*/

BOOL bsAddComment (OBJ obj, const char* text, char type)
{
  BS bs ;
  BS new ;

  if (!obj)
    return FALSE ;
  if (obj->magic != OBJ_MAGIC)
    messcrash ("bsAddComment received a non magic B obj (%s)",
	       name(obj->key)) ;

  if (! text || !strlen(text))
    return FALSE;

	/* add to end of list of comments */
  bs = obj->curr ;
  if (bs->right && bsIsComment(bs->right))
    { if (!strcmp(text,name(bs->right->key)))
	return FALSE ;
      for (bs = bs->right ;
	   bs->down && bsIsComment(bs->down) ;
	   bs = bs->down) 
	if (!strcmp(text,name(bs->down->key)))
	  return FALSE ;
    }

  new = BSalloc() ;
  new->timeStamp = sessionUserKey() ;
  new->up = bs ;

  if (type == 'C')
    lexaddkey (text, &new->key, _VComment) ;
  else
    messcrash ("Attempt to add unknown comment type %c", type) ;

  if (isXrefing)
    { if (!obj->xref) obj->xref = stackCreate (64) ;
      push (obj->xref, new->key, KEY) ;
      push (obj->xref, _Quoted_in, KEY) ;
    }

  if (bs == obj->curr)
    { new->down = bs->right ;
      bs->right = new ;
    }
  else
    { new->down = bs->down ;
      bs->down = new ;
    }
  if (new->down)
    new->down->up = new ;

  obj->curr = new ; obj->modCurr = 0 ; /* needed so can change timeStamp */
  obj->flag |= OBJFLAG_TOUCHED ;
  cacheMark (obj->x) ;
  return TRUE ;
}

/***********************************************************/

static BOOL bsAddCell (OBJ obj, BSunit *u)
{
  KEY type = bsType (obj, _bsRight) ;

  if (!type)
   return TRUE ;
  while (class (type) && KEYKEY(type) == 1)
    {
      if (!bsPushObj (obj))
	return FALSE ;
      else
	type = bsType (obj, _bsRight) ;
    }
  if (type < _LastC)
    return u->s ? (bsAddData (obj, _bsRight, type, u->s), TRUE) : FALSE ;
  else if (type < _LastN)
    return (bsAddData (obj, _bsRight, type, u), TRUE)  ;
  else
    return u->k ? (bsAddKey (obj, _bsRight, u->k), TRUE)  : FALSE ;
}

/***************************************/

    /* bsFlatten() flattens a subtree into an array of bsUnits
       provided by the user. The array contains a multiple of n
       BSunits, each containing an int, float, char*
       or KEY according to the model, which the user is presumed
       to know.  Missing data are represented by 0.
    */

BOOL bsAddTruncatedArray (OBJ obj, KEY target, Array a, int arrWidth, int objWidth)
{
  int i, j ;
  BSunit *u ;

  if (!arrayExists(a))
    messcrash("bsAddArray received a non existing array %d", a) ;
  if (a->size != sizeof(BSunit))
    messcrash("bsAddArray received an array of wrong type (should be BSunit)") ;
  if (!objWidth || arrWidth < objWidth)
    messcrash ("Incompatible width in array and obj") ;
  if (target <= _bsDown)
    messcrash ("bsAddArray called on bad target %s",
	       target ? name(target) : "zero") ;

  if (!bsAddTag (obj, target))
    return FALSE ; /* checks obj, magic, writeaceess etc */
  bsRemoveTag(obj, target) ;
  bsAddTag (obj, target) ; /* second try, cannot fail */

  for (i = 0 ; i < arrayMax(a) ; i += arrWidth)
    {
      bsGoto (obj, 0) ;
      bsAddTag (obj, target) ; /* second try, cannot fail */
      for (j = 0 ; j < arrWidth ; j++)
	{
	  u = arrp(a, i + j, BSunit) ;
	  if (!bsAddCell (obj, u))
	    break ; /* stop this line */
	}
    }
  bsGoto (obj, 0) ; /* pop out of possible constructed types */
  return TRUE ;
}

#endif /* READONLY */

/******* mark package to allow return to the same point ********/

#define MARK_MAGIC 843843
typedef struct MarkStruct  { int magic ;
			     BS curr, modCurr, localRoot ;
			     int localModel ;
			   } *MARK ;

BSMARK bsHandleMark (OBJ obj, BSMARK mark, AC_HANDLE hh)
{
  if (!mark)
    { mark = (BSMARK) halloc (sizeof (struct MarkStruct), hh) ;
      ((MARK)mark)->magic = MARK_MAGIC ;
    }
  else if (((MARK)mark)->magic != MARK_MAGIC)
    messcrash("Bad mark in bsMark") ;

  ((MARK)mark)->curr = obj->curr ;
  ((MARK)mark)->modCurr = obj->modCurr ;
  ((MARK)mark)->localRoot = obj->localRoot ;
  ((MARK)mark)->localModel = obj->localModel ;
  return mark ;
}

BSMARK bsMark (OBJ obj, BSMARK mark)
{ return bsHandleMark (obj, mark, 0) ; }

char *bsGetMarkAsString (BSMARK mark)
{
  return mark ? messprintf("%lx ->curr=%lx", (unsigned long)mark, (unsigned long)((MARK)mark)->curr) : "(NULL mark)" ;
}

void bsGoto (OBJ obj, BSMARK mark)
{
  if (mark)
    { if (((MARK)mark)->magic != MARK_MAGIC)
	messcrash("Bad mark in bsGoto") ;

      obj->curr = ((MARK)mark)->curr ;
      obj->modCurr = ((MARK)mark)->modCurr ;
      obj->localRoot = ((MARK)mark)->localRoot ;
      obj->localModel = ((MARK)mark)->localModel ;
    }
  else
    { obj->curr = obj->localRoot = obj->root ;
      obj->localModel = localModel (obj->key) ;
      obj->modCurr = array(tabModel,KEYKEY(localModel(obj->key)),BS) ;
    }
}

/********** for deletion from BS trees ***************/

static void xrefPrune (OBJ obj, BS bs, BS bsm)
     /* recursive, NB can be called on comments with bsm == 0 */
{
#ifdef DELETE_DEBUG
  static int level = 1 ;
  printf ("%*s%s",2*level,"",name(bs->key)) ;
#endif

	/* check for cross referencing */
  if (isXrefing)
    { if (class(bs->key) == _VComment)
	{ if (!obj->xref) obj->xref = stackCreate (64) ;
	  push (obj->xref, bs->key, KEY) ;
	  push (obj->xref, (_Quoted_in | DELETE_BIT), KEY) ;
#ifdef DELETE_DEBUG
      printf ("\txref") ;
#endif
	}
      else if (bsm && bsModelMatch (bs, &bsm) && KEYKEY(bsm->n.key))
	{ if (!obj->xref) obj->xref = stackCreate (64) ;
	  push (obj->xref, bs->key, KEY) ;
	  push (obj->xref, (KEYKEY(bsm->n.key) | DELETE_BIT), KEY) ;
#ifdef DELETE_DEBUG
	  printf ("\txref") ;
#endif
	}
    }

/* if not Xrefing bsm becomes meaningless, but it doesn't matter! */

#ifdef DELETE_DEBUG
  printf ("\n") ;
#endif

  if (bs->down)
    xrefPrune (obj,bs->down,bsm) ;

#ifdef DELETE_DEBUG
  ++level ;
#endif

  if (bs->right)
    {
      if (bsm)
	xrefPrune (obj,bs->right,bsm->right) ;
      /* not bsModelRight(bsm) so not to jump into subObjects */
      else
	xrefPrune (obj,bs->right,0) ;
    }

#ifdef DELETE_DEBUG
  --level ;
#endif

  BSfree (bs) ;
}

/*********************************/

BOOL bsRemove (OBJ obj)	  /* removes curr and subtree to the right
			     returns curr to the root of the tree */
{
  BS    bs, bsm ;

  if (!obj)
    return FALSE ;
  if (obj->magic != OBJ_MAGIC)
    messcrash ("bsRemove received  a non magic B obj (%s)",name(obj->key)) ;

  bs = obj->curr ;
  bsm = obj->modCurr ;
   /* for this routine only, modCurr can be 0, if bs is a comment */

  if (!bs || bs == obj->root)	/* can't remove the root node */
    return FALSE ;

  /* it is impossible to know how many tags get removed
   * right of the tag we wish to remove 
   */
  ac_free (obj->x->assHandle) ; obj->x->tagAss = 0 ;

  if (bs->up->right == bs)
    bs->up->right = bs->down ;
  else
    bs->up->down = bs->down ;
  if (bs->down)
    bs->down->up = bs->up ;
  bs->down = 0 ;
#ifdef DELETE_DEBUG
  printf ("bsRemove from %s\n", name(obj->key)) ;
#endif
  xrefPrune (obj, bs, bsm) ;

	/* must return to root because of REPEAT, FREE problems */
  bsGoto (obj, 0) ;
  obj->flag |= OBJFLAG_TOUCHED ;
  cacheMark (obj->x) ;
  return TRUE ;
}

/**********************************/

BOOL bsPrune (OBJ obj)
     /* removes curr and everything left of it that is unique.
	NB must go right back to root and come out to keep model
	straight because of possible REPEAT etc.
	Returns curr, modCurr to root at end */
{
  BS bs, bsm ;
  KEY targetKey = 0 ;
  static Stack path = 0 ;
  Associator ass = 0 ;
  /*   BOOL debug = FALSE ; */

  if (!obj)
    return FALSE ;
  if (obj->magic != OBJ_MAGIC)
    messcrash ("bsPrune received  a non magic B obj (%s)",name(obj->key)) ;

  path = stackReCreate (path, 32) ;

  bs = obj->curr ;  
  
  /*   if (bs && class (bs->key) == _VSequence) debug = TRUE ; */
  targetKey = obj->modCurr && obj->root == obj->localRoot ? bs->key : 0 ;
  
  while (!bs->down && bs->up && bs->up->right == bs
	 && bs->up != obj->root)		/* unique non-root */
    bs = bs->up ;

  if (obj->x->tagAss &&
      targetKey
      )
    {
      BS bs1 = bs ; 
      while (bs1->up && bs1->up->down == bs1
	     && bs1->up != obj->root)		/* look for a tag */
	bs1 = bs1->up ;
      if (bs1->up && bs1->up->right == bs1 &&  bs1->up != obj->root)
	{
	  bs1 = bs1->up ;
	  if (! class (bs1->key) && 
	      assFind (obj->x->tagAss, assVoid (bs1->key), &ass) &&
	      ass)
	    {
	      BS bs2 = 0, bs3 = 0 ;
	      if (assFind (ass, assVoid (targetKey), &bs2) && bs2)
		{
		  assRemove (ass, assVoid (targetKey)) ; 
		  /* if (debug) fprintf (stderr, "## bsPruned removed bs1->key = %s   targetKey=%s, bs=%p bs->up=%p\n", name(bs1->key) ,name(targetKey), bs2, bs2->up) ;  */
		}
	      if (assFind (ass, assVoid (-2), &bs3) && bs3 == bs2)
		{
		  assRemove (ass, assVoid (-2)) ; 
		  if (bs3->up)
		    assInsert (ass, assVoid (-2), bs3->up) ;
		}
	    }
	}
    }

  while (bs != obj->root)
    { push (path,bs,BS) ;
      while (bs->up && bs->up->down == bs)
	bs = bs->up ;
      bs = bs->up ;
    }

  bsm = array(tabModel,KEYKEY(localModel(obj->key)),BS) ;
  while (!stackEmpty(path))
    { 
      bs = pop (path,BS) ; ass = 0 ;
 	  
      bsm = bsModelRight(bsm) ;
      if (bsIsComment (bs))
	{ bsm = 0 ;
	  break ;		/* must be the final elt */
	}
      if (!bsModelMatch (bs, &bsm))
	messcrash ("Model screwup in bsPrune: %s unmatched in %s",
		   name(bs->key), name(obj->key)) ;
    }

  obj->curr = bs ;
  obj->modCurr = bsm ;
  return bsRemove (obj) ;
}

/******************************************************/

BOOL bsFuseObjects (KEY old, KEY new)
{
#ifdef READONLY
  return FALSE ;
#else
  static Stack s = 0 , s1 = 0 , s2 = 0 ;
  char *cp ;
  extern void cacheKill(KEY) ;

  s = stackReCreate(s, 300) ;
  s1 = stackReCreate(s1, 300) ;
  s2 = stackReCreate(s2, 300) ;

  /* for use in lexalias() - call this BEFORE aliasing the keys.
     idea is to ace dump old then new to a Stack both with object name new,
     then kill the old object to undo the Xrefs and parse the buffer.
  */

  dumpKey(old, 0, s) ;
  { OBJ obj = bsUpdate(old) ;

    bsKill(obj) ;    /* To unXref and clean the second cache */
    cacheKill (old) ; /* to unreference it in cache2 */
    bsTreeKill(old) ; /* To remove old from first cache, hence from disk */
  }

  dumpKey(new, 0, s2) ;  /* one object result in a single string */

  pushText(s1, messprintf("%s : \"%s\"", className(new), name(new))) ;
  catText(s1, "\n") ;

  cp = stackText(s, 0) ;
  while (*cp++ != '\n') ;    /* skip name(key) */
  catText(s1, cp) ;
  catText(s1, "\n") ;
  catText(s1, "\n") ;

  cp = stackText(s2, 0) ;     /* add back new, to give higher priority */
  catText(s1, cp) ;
  catText(s1, "\n") ;
  catText(s1, "\n") ;

#ifdef DELETE_DEBUG
  stackCursor(s1, 0) ;
  while (cp = stackNextText(s1))
    printf("%s\n",cp) ;
#endif

  return parseBuffer (stackText(s1,0), 0) ;
#endif
}

/***************************************************************************/
/******* Local Routines - all concern the formation and use of paths *******/

        /* make set of paths for model only when needed by bsUpdate
           set of paths held as an associator: tag -> path
           these are stored in a table over models: tabAssPath
           findHook() sets ->curr to the end of the path
        */

static Array workPath ;         /* 2 statics needed for makepaths */
static int depth ;
static int pathTable ;  /* JTM, to help debugging the errors in models.wrm*/

/*****************************************/

void bsMakePaths (int table)
{
  OBJ Model = 0 ;
  KEY key ;
  if (!tabModel)
    { tabModel = arrayCreate(256,BS) ;
      tabAssPath = arrayCreate(256,Associator) ;
    }

#if !defined(NEW_MODELS)
  key = KEYMAKE(table,0) ;
  if (pickList[table].type != 'B')
    messcrash ("Trying to update a non bs voc with bs package");
#else
  key = table ; table = KEYKEY(key) ;
#endif
  if (array(tabAssPath,table,Associator))
    assClear (array(tabAssPath,table,Associator)) ;
  else
    array(tabAssPath,table,Associator) = assCreate () ;

  if (array(tabModel, table, BS))
    bsTreePrune (array(tabModel, table, BS)) ;

#if !defined(NEW_MODELS)
    if (table == _VModel)  /* case Model */
#else
     if (!strcmp(name(table), "#Model"))
#endif
    { BS bs = BSalloc() ;
      array(tabModel, table, BS) = bs ;
      bs->key = _ANY ; /* trivial model TREE, empty tabAssPath */
      return ;
    }

  Model = bsCreate(key) ;
  if (!Model)
    return ;

  if (!Model || ! Model->root || ! Model->root->right)
    messcrash("%s: %s\n%s",
	      "Missing model of class", className(key),
	      "use Read-Models from main-menu and/or check wspec/models.wrm") ;

  array(tabModel,table,BS) = bsTreeCopy(Model->root, 0) ; /* no memeory handle needed here */
  bsDestroy (Model) ;

  if (!workPath)
    workPath = arrayCreate (4,KEY) ;
  depth = 0 ;
  pathTable = table ;

  makePaths (array(tabAssPath,table,Associator), (array(tabModel,table,BS))->right, FALSE) ;
}

/******************************************/

static void bsRFree (BS bs)
{
  if (!bs) return ;
  if (bs->right && bs != bs->right) /* models are special */
    bsRFree (bs->right) ;
  if (bs->down && bs != bs->down) /* models are special */
    bsRFree (bs->down) ;
  bs->right = bs->down = 0 ;
  BSfree (bs) ;
}

void bsMakePathShutDown (void)
{
  Associator *asp ;
  int n = tabAssPath ? arrayMax(tabAssPath) : 0 ;
  BS *bsp ;

  asp = n ?  arrp (tabAssPath, 0, Associator) - 1 : 0 ;
  while (asp++, n--)
    {
      void *vp = 0, *wp = 0;
      Array aa ;
      if (*asp)
	while (assNext (*asp, &vp, &wp))
	{
	  aa = (Array) wp ;
	  arrayDestroy (aa) ;
	}
      assDestroy (*asp) ;
    }
  arrayDestroy (tabAssPath) ;
  arrayDestroy (workPath) ;
  n = tabModel ? arrayMax(tabModel) : 0 ;
  bsp = n ? arrp (tabModel, 0, BS) - 1 : 0 ;
  while (bsp++, n--)
    bsRFree (*bsp) ;
  arrayDestroy (tabModel) ;
}

/******************************************/

static void makePaths (Associator ass, BS bs, BOOL isUnique)
{ BS rehook = 0 ;		/* recursively makes paths */
  char *vp ;

  if (bs->key == _UNIQUE)	/* attach UNIQUE to right-hand nodes */
    { BS bs_u = bs ;
               /* present limitations that we may wish to undo */
      if (bs->down)
	messcrash ("UNIQUE right above something else after tag %s in model %s",
	       name(bs->up->key), pickClass2Word(pathTable)) ;
      if (!bs->up || bs != bs->up->right)
	messcrash ("UNIQUE right below something else after tag %s in model %s",
	       name(bs->up->key), pickClass2Word(pathTable)) ;
      if (!(bs = bs_u->right))
	return ;
      if (bs_u->down)
	rehook = bs_u->down ;
      if (bs_u->up->right == bs_u)
	bs_u->up->right = bs ;
      else
	bs_u->up->down = bs ;
      bs->up = bs_u->up ;
      BSfree (bs_u) ;
      isUnique = TRUE ;
    }

  if (bs->down)
    makePaths (ass, bs->down, isUnique) ;  /* don't include bs in path */

  if (rehook)
    { BS temp = bs ;
      while (temp->down)
	temp = temp->down ;
      rehook->up = temp ;
      temp->down = rehook ;
    }

  /* Note that we need to eat up XREF and UNIQUE all the way right */
  /* so if (!bsIsTag(bs) || bs->key == _ANY) return ; would be wrong */

  array(workPath,depth,KEY) = bs->key ;
  ++depth ;
  array(workPath,depth,KEY) = 0 ;

  vp = assVoid(bs->key) ;
  if (bsIsTag(bs) && bs->key != _ANY &&
      !assInsert (ass,vp,arrayCopy(workPath)))
    messcrash ("Duplicate tag %s in model %s",
	       name(bs->key), pickClass2Word(pathTable)) ;

  bs->n.key = 0 ; /* compress unique and xref keys into ->n field
                   * and swallows the XREF tag nodes
		   */
  if (bs->right && bs->right->key == _XREF)
    { BS xrefnode = bs->right, tagnode = xrefnode->right ;
      if (!tagnode)
	messcrash ("Missing tag after XREF in model of class %s",
		   pickClass2Word(pathTable)) ;
      bs->n.key = tagnode->key ; /* pick up back ref tag */
      bs->right = tagnode->right ;
      if (bs->right)
	bs->right->up = bs ;
      BSfree (tagnode) ;
      BSfree (xrefnode) ;
    }
  else if (pickXref(class(bs->key)))
    bs->n.key = _Quoted_in ;

  if (bs->right && bs->right->key == _COORD)
    { BS temp = bs->right ;
      bs->n.key |= COORD_BIT ;
      bs->right = bs->right->right ;
      if (bs->right)
	bs->right->up = bs ;
      BSfree (temp) ;
    }

  if (bs->right && bs->right->key == _SORTED)
    { BS temp = bs->right ;
      bs->n.key |= SORTED_BIT ;
      bs->right = bs->right->right ;
      if (bs->right)
	bs->right->up = bs ;
      BSfree (temp) ;
    }

  /* RMD How does the _VConstraint key get set up? : in models.c (mieg) */

  if (bs->right && class(bs->right->key) == _VConstraint)
    { BS constraintnode = bs->right ;
      bs->n.key = constraintnode->key ; /* pick up back ref tag */
      bs->right = constraintnode->right ; /* RMD **IMPORTANT** */
      if (bs->right)
	bs->right->up = bs ;
      BSfree (constraintnode) ;
    }

  if (isUnique)
    bs->n.key |= UNIQUE_BIT ;

  if (bs->right)
    {
      if (bs->right->key == _REPEAT)
	bs->right = bs ;
      else
	makePaths (ass, bs->right, FALSE) ; /* UNIQUE does not prop right */
    }

  --depth ;         /* important to restore to original value */
}

/******************************************/

KEYSET bsGetPath (int classe, KEY tag)
{
  Associator ass ;
  KEY model = localModel (KEYMAKE (classe,0)) ;
  int table = KEYKEY(model) ;
  KEYSET path ;
  char * vt ;

  if (!model)
    return 0 ;

  CHECK_MODEL(model) ;

  ass = array(tabAssPath,table,Associator) ;

  if (tag == 0)
    return 0 ;

  if (tag >= lexMax(0))
    messcrash ("bsGetPath receives tag %d < 0 or > tagMax = %d, \n%s", tag, lexMax(0),
	       "Probably, the tag variable used was never initialised") ;
   vt = assVoid(tag) ;
   if (!assFind (ass,vt,&path))
     return 0 ;
   return arrayCopy(path) ;
}

/*****************************************/

        /*  set curr to the tag position in obj.  Also set modCurr.
            construct path to tag if necessary only if makeIt is TRUE.
            return TRUE if found (made), else FALSE
        */

static BOOL findTag (OBJ obj, KEY tag, BOOL makeIt)
{
  Array path ;
  int model = obj->localModel ;
  int table = KEYKEY(model) ;
  Associator ass = array(tabAssPath,table,Associator) ;
  BS    bs = obj->localRoot, bsm, new, bs1, top ;
  int   i = 0 ;
  BOOL  ismodif = FALSE, isDown ;
  KEY   key ;
  char * vt ;

  if (tag == 0)
    return FALSE ;

  if (tag >= lexMax(0))
    messcrash ("findTag receives tag %d < 0 or > tagMax = %d, \n%s", tag, lexMax(0),
	       "Probably, the tag variable used was never initialised, have a look at tags.c") ;

  /* XQX
  * Conjecture: here, ass is a list of path data for each tag in the
  * object.  path data is an array.  it tells you which tags you must
  * pass through to get to a given tag.
  *
  * So, here we fill in the array "path" with the series of tags to walk
  * through to get to tag "tag".
  */

  vt = assVoid (tag) ;
  if (!assFind (ass,vt,&path))
    return FALSE ;

  /* XQX
  * Conjecture:
  * Extract the "bsm" for this object class from the model structure.
  * Presumably, this bsm points at the root.
  */
  bsm = array(tabModel,table,BS) ;

  /* XQX
  * bs a pointer into the actual object data.  It begins at the root
  * of the object.
  */

  /* XQX
  * Conjecture:  All objects and models contain a dummy node at the
  * top.  You have to follow ->right to get to the first top-level
  * tag.  (That is why this loop begins by moving right on the model
  * and right on the object.)
  */

  i = 0;
  while ((key = arr(path,i++,KEY)))
    {
      /* XQX
      * Each successive key we extract is the next one in the path to
      * the tag we are looking for.
      */
      bsm = bsModelRight (bsm) ;

      /* XQX
      * Skip over comment nodes in the tree.  Comment nodes have
      * only ->down, not ->right.  I wonder which node the comment
      * applies to...
      */
      isDown = FALSE ;
      for (top = bs->right ;
	   top && bsIsComment(top) ;
	   top = top->down)
	{ bs = top ;
	  isDown = TRUE ;
	}

      /* XQX
      * Remembering that the root of the tree is to the left and
      * branches are arranged one above the other, we walk down
      * the branches at the current level looking for the one that
      * has the tag we need.
      */
      for (bs1 = top ; bs1 ; bs1 = bs1->down)
	if (bs1->key == key)
	  {
	    /* XQX
	    * we found the row that has the key we want, so
	    * remember that.  Walk the model cursor down to the
	    * same key.
	    */
   	    bs = bs1 ;
	    while (bsm->key != key)
	      bsm = bsm->down ;
	    goto gotIt ;
	  }

      /* if get here key is not there */
      if (!makeIt)
	return FALSE ;

      /* XQX
      * if we get here, we are about to add the tag to the object tree.
      *
      * top is still pointing at the first non-comment row of the tree.
      */

      /* XQX
      * Don't quite know what this is doing.  The ->up field is the node
      * "above" us.  If it's down pointer points at us, then the tag
      * we are looking at is the first one - we want to trash the current
      * pointer and insert a new one instead.
      *
      * but if we get here, we know that the tag we are looking for is
      * NOT PRESENT IN THE OBJECT, and therefore CANNOT be conflicting
      * with a "unique" attribute.
      */
      if (top && bsm->n.key & UNIQUE_BIT) /* UNIQUE */
	{ if (top->up->down == top)
	    top->up->down = 0 ;
	  else
	    top->up->right = 0 ;
	  xrefPrune (obj, top, bsm) ;
	  top = 0 ;
	}

      /* XQX
      * create a new tag node to store in the object tree; label it with
      * the tag we are inserting
      */
      new = BSalloc() ;
      new->timeStamp = sessionUserKey() ;
      new->key = key ;
      ismodif = TRUE ;

      /* XQX
      * Now we proceed through the model tree looking for the key that we
      * are trying to add.  I don't understand the rest, but I think it
      * leaves bsm pointing at the model entry for the tag that we just
      * created in new, and new->up pointing at the tag "above" where we
      * want to be inserted.
      */
      for (bs1 = top ; bsm->key != key ; bsm = bsm->down)
	while (bs1 && (class(bs1->key) == _VUserSession || bs1->key == bsm->key))
	  { new->up = bs1 ;
	    bs1 = bs1->down ;
	  }

      /* XQX
      * now we connect the new node into the tree
      */
      if (new->up)
	{ 
          /* XQX
	  * if there is a node in the tree above us, insert ourselves as
	  * the first row of data below it.
   	  */
          new->down = new->up->down ;	/* previous element is below us */
	  new->up->down = new ;		/* we are new first element below parent */
	  if (new->down)
	    new->down->up = new ;	/* the guy below us points up to us */
	}
      else
	{ 
	  /* XQX
	  * we are inserting at the top of the object.
	  */
	  new->down = top ;		/* down from us is everything in the object */
	  if (top)
	    top->up = new ;		/* the thing below us has us for "up" */

	  /* XQX
	  * sim sala bim... :(
	  */
	  new->up = bs ;
	  if (isDown)
	    bs->down = new ;
	  else
	    bs->right = new ;
	}

      bs = new ;

    /* XQX
    * Here, we either found the tag we were looking for, or we inserted it
    * in the place we were hoping to see it.  But remember we are in a loop
    * walking through the list of tags we must traverse to get to the one
    * we want.  This is where we go on to the next iteration of the loop.
    */
gotIt : ;
    }

  /* XQX
  * having fallen out of the loop, we leave the cursor sitting on the last
  * tag that we looked at.
  */
  obj->modCurr = bsm ;
  obj->curr = bs ;
  if (ismodif)
    { obj->flag |= OBJFLAG_TOUCHED ;
      cacheMark (obj->x) ;
    }
  return TRUE ;
}

/*********************************************/

        /* for testing ahead in the model to see what types allowed */
        /* currently you must go once right then down the column
           to make sure - not really done properly */

KEY bsType (OBJ obj, KEY target)
{ KEY key ;
  BS bsm = obj->modCurr ;
  int table ;

  if (!bsm)
    return FALSE ;

  switch (target)
    {
    case _bsHere:  break ;
    case _bsDown:  bsm = bsm->down ; break ;
    case _bsRight: bsm = bsm->right ; break ;
    default:   messcrash ("bsType called with %s - must be right or down",
                          name(target)) ;
    }

  if (!bsm)
    return 0 ;
  key = bsm->key ;
  if (class(key) == _VModel) /* newsystem */
    { table = pickWord2Class (name(key) + 1) ;
      if (*name(key) == '#')
	return KEYMAKE (table,1) ;
      else
	return KEYMAKE (table,0) ;
    }
  return key ;
}

/********** typematch bsm in column to bs ************/

BOOL bsModelMatch (BS bs, BS *bsmp)
{
  KEY bskey = bs->key ;
  BS  bsm = *bsmp ;
#ifndef NEW_MODELS
  int type = class(bskey) ? KEYMAKE(class(bskey),0) : bskey ;
#else
  int type = class(bskey) ? pickList[class(bskey)].classe : bskey ;
#endif
  if (!bsm || !bsm->up)
    return FALSE ;

  if (bsm->key == type)
    return TRUE ;

  while (bsm->up->down == bsm)    /* go to top of allowed column */
    bsm = bsm->up ;

  while (bsm)                  /* go back down looking for a match */
    { if (bsm->key == _ANY || bsm->key == type)
        { *bsmp = bsm ;
	  return TRUE ;
        }
      bsm = bsm->down ;
    }
  return FALSE ;
}

/*****************************************************/

#ifndef READONLY

#undef XREF_DEBUG

static void doXref (OBJ originalObj)
{
  Stack s  = originalObj->xref ;
  OBJ obj ;
  KEY key, tag , original = originalObj->key ;
  static KEY trueOriginal = 0, level = 0 ;

  isXrefing = FALSE ; /* This flag is also used in cacheUpdate and xclient */

  if (!level++)
    trueOriginal = original ;
  stackCursor(s,0) ;
  while (!stackAtEnd (s))
    {
      key = stackNext(s,KEY) ;
      tag = stackNext(s,KEY) ;

      if (!KEYKEY(key))
	continue ; /* Can't xref to a model */
      if (key == original ||
	  key == trueOriginal)
	continue ; /* Do not self Xref */
      if (BSTEST)
	cacheCheck(original, originalObj->x) ;
      if ((obj = bsUpdate (key)))
        { if (tag & DELETE_BIT) /* delete: obj->curr is in place */
	    { if (bsFindKey (obj, KEYKEY(tag), original)) /* locates obj->curr */
		bsPrune (obj) ;
	    }
	  else
	    bsAddKey (obj, KEYKEY(tag), original) ; /* fails if original does not type check */

	  if (BSTEST) cacheCheck(key, obj->x) ;
	  bsSave (obj) ;
        }
      else
	messout ("Can't crossref %s onto %s because %s is locked",
		 name(original),name(key),name(key)) ;
      /* could put it onto a global stack to try at the session end? */
    }

  level-- ;
  isXrefing = TRUE ;
}

#endif

/****************************************************/
/***** routines for updating in treedisp.c **********/
/* mieg: apr 20-2004
 * justKnownTags is used in tag chooser in table maker
 * to limit the portion of model fused to existing tags
 */
static void bsSubFuse (BS bs, BS bsm, int cl, KEY tag1, int justKnownTags)
     /* recursively fuses model into object */
{    /* assumes bs matches bsm - does the whole column right */
  BS new,bsr,top,old ;

  if (!bsm) return ;
  if ((bsm->n.key & UNIQUE_BIT) && !bsIsTag(bs))
    bs->size |= UNIQUE_FLAG ;
  if (!bs->bt)
    bs->bt = BTalloc() ;
  bs->bt->bsm = bsm ;

  if (!bsModelRight(bsm) ||
      ((bs->size & MODEL_FLAG) && bsModelRight(bsm) == bsm))
				/* 2nd poss. if REPEAT */
    return ;

  	/* check all comments have bt->bsm */
  for (top = bs->right ; top ; top = top->down)
    if (!top->bt)
      top->bt = BTalloc() ;

	/* top will be the top of the column - hook it on later */
  for (top = bs->right ; top && bsIsComment (top) ; top = top->down) ;
  if (top)
    top->up = 0 ;

     /* JTM fuse to types */
  if (top)
    bsm = bsModelRight(bsm) ;	/* keep expanding the model */
  else
    bsm = bsm->right ;		/* don't expand subtypes */

  if (!bsIsTag(bsm))
    justKnownTags = 0 ;
  
    /* take the model keys in the column, one by one */
  for ( ; bsm ; bsm = bsm->down)
    {
      switch (justKnownTags)
	{
	case 0: /* show all tags */
	  break ;
	case 1: /* show only tags used in the database */
	  if (bIndexTagCount (cl, bsm->key) < 1) 
	    /* -1, dont know in bindex
	     *  0, not evaluated
	     */
	    continue ;
	  break ;
	case 2: /* show only already filled tags */
	  if (tag1 == bsm->key)
	    justKnownTags = 5 ; /* == 1 | 4 to restore */
	  else if (top && bsIsTag (top)  && ! (top->size & MODEL_FLAG)) ;
	  else
	    continue ;
	  break ;
	}
      /* first fuse to anything already there that matches */
      for (bsr = top, old = 0 ; bsr ; old = bsr, bsr = bsr->down)
	if (bsr->key == bsm->key ||
	    (class(bsr->key) && class(bsr->key) == class(bsm->key)))
	  { 
	    bsSubFuse (bsr, bsm, cl, tag1, justKnownTags & 0x3) ; 
	    if (bsIsTag(bsm) || bsm->n.key & UNIQUE_BIT)
	      goto skipNewNode ;
	    if (justKnownTags == 5)
	      goto skipNewNode ;
	  }
      if (justKnownTags == 0x2)
	goto skipNewNode ;
      /* now make new node */
      new = BSalloc () ;
      new->size |= MODEL_FLAG ;		/* flag that it is from model */
      new->key = bsm->key ;
      if (bsIsType(bsm))
	new->size |= SUBTYPE_FLAG ;
      
      /* insert below last occurence, or previous tag, so search backwards */
      
      if (old)
	for (bsr = old ; bsr ; bsr = bsr->up)
	  if (bsIsTag (bsr) || bsr->key == bsm->key || bsr->key == bsm->up->key ||
	      (class(bsr->key) && class(bsr->key) == class(bsm->key)))
	    break ;
      
      if (bsr)			/* put new below it */
	{ 
	  new->up = bsr ;
	  new->down = bsr->down ;
	  bsr->down = new ;
	  if (new->down)
	    new->down->up = new ;
	}
      else			/* put new at the top */
	{
	  new->down = top ;
	  if (top)
	    top->up = new ;
	  top = new ;
	}
      
      if (justKnownTags != 5 && !bsIsType(bsm))
	bsSubFuse (new, bsm, cl, tag1, justKnownTags & 0x3) ;
skipNewNode: ;
      if (justKnownTags == 5) justKnownTags = 2 ; /* restore */
    }

	/* now hook our new column back onto bs (+ comments) */

  if (bs->right && bsIsComment(bs->right))
    { 
      for (bs = bs->right ;
	   bs->down && bsIsComment(bs->down) ;
	   bs = bs->down)
	;
      bs->down = top ;
    }
  else
    bs->right = top ;
  if (top)
    top->up = bs ;
}

void bsSubFuseType (OBJ obj, BS bs)
{
  KEY model = localModel(bs->key) ;
  BOOL isUp = FALSE ;
  BS new, bsup ;

  /* unhook*/
  bsup = bs->up ;
  if (bsup->down == bs)
    isUp = TRUE ;
  else if (bsup->right == bs)
    isUp = FALSE ;
  else
    messcrash("hook error in bsSubFuseType") ;
  bs->up = 0 ;
  bsTreePrune(bs) ;

  new = BSalloc () ;

  CHECK_MODEL(model) ;
  bsSubFuse (new, array(tabModel,KEYKEY(model),BS), class (bs->key), 0, 0) ;

  /*  rehook */
  if (isUp)
    bsup->down = new->right ;
  else
    bsup->right = new->right ;
  new->right->up = bsup ;

  BSfree(new) ;
}

/* mieg: april 2004, added fine control of which tags are fused */
void bsFuseModel (OBJ obj, KEY tag1, int justKnownTags)
{ KEY model = localModel(obj->key) ;

  CHECK_MODEL(model) ;
  ac_free (obj->x->assHandle) ; obj->x->tagAss = 0 ;
  bsSubFuse (obj->root, array(tabModel,KEYKEY(model),BS), class (obj->key), tag1, justKnownTags) ;
}

/*********************************/

BS bsModelRight (BS bsm)
{
  if (!bsm)
    return 0 ;

  bsm = bsm->right ;

  if (bsm && bsIsType (bsm))	/* if subtype, jump into it */
    { KEY model = localModel (bsm->key) ;
      BS bs ;
      CHECK_MODEL(model) ;
      bs = array(tabModel,KEYKEY(model),BS) ;
      bsm = bs ? bs->right : 0 ;
    }

  return bsm ;
}

BS bsModelRoot (OBJ obj)
{
  KEY model ;
  if (!obj)
    return 0 ;

  model = localModel(obj->key) ;
  CHECK_MODEL(model) ;
  return array(tabModel,KEYKEY(model),BS) ;
}

/*********************************/

static float fShift, dfShift ;
static int iShift, diShift ;
static BOOL isShiftAfter, tShift ;
static KEY shiftTimeStamp ;

static void bsCoordShift2 (BS bs, BS bsm)
{
		/* do column starting with bs */
  do
    if (bsModelMatch (bs, &bsm))
      { if (bsm->n.key & COORD_BIT)
	  switch (bsm->key)
	    {
	    case _Int:
	      if ((isShiftAfter && bs->n.i >= iShift) ||
		  (!isShiftAfter && bs->n.i < iShift))
		{ bs->n.i += diShift ;
		  tShift = TRUE ;
		  bs->timeStamp = shiftTimeStamp ;
		}
	      break ;
	    case _Float:
	      if ((isShiftAfter && bs->n.f >= fShift) ||
		  (!isShiftAfter && bs->n.f < fShift))
		{
		  bs->n.f += dfShift ;
		  {
		    char zbuf[128] ;
		    float zf =  bs->n.f ;
		    sprintf (zbuf, "%1.7g", zf) ; sscanf (zbuf, "%f", &zf) ;
		    bs->n.f = zf ;
		  }
		  tShift = TRUE ;
		  bs->timeStamp = shiftTimeStamp ;
		}
	      break ;
	    default:
	      messcrash ("bsCoordShift set for illegal type %s",
			 name (bsm->key)) ;
	    }
	if (bs->right)
	  bsCoordShift2 (bs->right, bsModelRight (bsm)) ;
      }
  while ((bs = bs->down)) ;
}

void bsCoordShift (OBJ obj, KEY tag, float x, float dx, BOOL isAfter)
{
  if (!tag)			/* go to root */
    bsGoto (obj, 0) ;
  else if (!bsFindTag (obj, tag))
    return ;

  if (!dx) return ;
  if (obj->curr == obj->root)
    bsFindTag (obj, _bsRight) ;	/* get into real object */

  tShift = FALSE ;
  fShift = x ; dfShift = dx ;
  iShift = (int)(x + ((x>0) ? 0.5 : -0.5)) ;
  diShift = (int)(dx + ((dx>0) ? 0.5 : -0.5)) ;
  shiftTimeStamp = sessionUserKey() ;
  isShiftAfter = isAfter ;

  bsCoordShift2 (obj->curr, obj->modCurr) ;
  if (tShift)
    { obj->flag |= OBJFLAG_TOUCHED ;
      cacheMark (obj->x) ;
    }
}

/******* next set handle an int *index lookup function ***********/

static int *cIndex ;

static void bsCoordIndex2 (BS bs, BS bsm)
{
		/* do column starting with bs */
  do
    if (bsModelMatch (bs, &bsm))
      { if (bsm->n.key & COORD_BIT)
	  switch (bsm->key)
	    {
	    case _Int:
	      if (bs->n.i != cIndex[bs->n.i])
		{ bs->n.i = cIndex[bs->n.i] ;
		  tShift = TRUE ;
		  bs->timeStamp = shiftTimeStamp ;
		}
	      break ;
	    default:
	      messcrash ("bsIndex set for illegal type %s",
			 name (bsm->key)) ;
	    }
	if (bs->right)
	  bsCoordIndex2 (bs->right, bsModelRight (bsm)) ;
      }
  while ((bs = bs->down)) ;
}

void bsCoordIndex (OBJ obj, int *index)
{
  bsGoto (obj, 0) ;
  bsFindTag (obj, _bsRight) ;	/* get into real object */

  tShift = FALSE ;
  cIndex = index ;
  shiftTimeStamp = sessionUserKey() ;
  bsCoordIndex2 (obj->curr, obj->modCurr) ;
  if (tShift)
    { obj->flag |= OBJFLAG_TOUCHED ;
      cacheMark (obj->x) ;
    }
}

/********** end of file **********/
