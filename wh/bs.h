/*  File: bs.h
 *  Author: J Thierry-Mieg and R Durbin
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
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
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Public functions of the bs package
 * Exported functions: 
 * HISTORY:
 * Last edited: Dec  7 09:46 2000 (edgrif)
 * Created: Thu Aug 26 17:54:10 1999 (fw)
 * CVS info:   $Id: bs.h,v 1.10 2011/10/27 22:45:33 mieg Exp $
 *-------------------------------------------------------------------
 */

#ifndef ACEDB_BS_H
#define ACEDB_BS_H

#include "acedb.h"

#ifndef DEFINE_OBJ
#define DEFINE_OBJ
  typedef void* OBJ ;
#endif

#include "bindex.h"

/************************************************************/
/*  classType 'B' : BS Trees, public functions.             */
/************************************************************/

typedef union			/* structure for unknown types */
  { int	  i ;
    float f ;
    KEY   k ;
    char* s ;
    mytime_t time ;
  } BSunit ;
 

KEY bsKey (OBJ obj) ;  /* return the key of the object */
char *bsName (OBJ obj) ; /* returns the name of the obj */

/* 4 functions to grab or release a BS tree in the object cache */
                                
OBJ bsCreate(KEY key) ;         /* read only */
                                
OBJ bsCreateCopy(KEY key) ;     /* temporary read/write without possibility of saving */
OBJ bsUpdate(KEY key) ;         /* read write 
				 * If you write priviledged code using bs__.h  
                                 * please call objMark() before touching
                                 * to prevent read only access by another
                                 * window */
OBJ bsClone(KEY key, OBJ source); /* like update, but replace data with a 
				     copy from source */
    
                            /* disregard recent updates */
void  bsDoDestroy(OBJ obj) ;
#define bsDestroy(obj)   (bsDoDestroy(obj) , (obj) = 0 ) 

BOOL bsIsDirty (OBJ obj) ; /* TRUE if obj needs to be saved */

                                 /* validates the latest update */
                                /* provided a call to objMark was made */
void bsDoSave        (OBJ obj) ;
#define bsSave(obj)    (bsDoSave(obj) , (obj) = 0 ) 

                            /* Zero and object and lexSets EMPTYSTATUS */
                            /* obj must be grabbed bsUpdate() */
void bsDoKill        (OBJ obj) ;
#define bsKill(obj)    (bsDoKill(obj) , (obj) = 0 ) 
void keySetKill (KEYSET ks) ; /* kills all 'A' and 'B' members of ks */

/* Several sorts of dumps */
void bsDump        (OBJ obj) ;             /* dumps on stdout in debug format */

void bsAceDump     (OBJ obj, FILE *fil, Stack s, BOOL isInternal) ;  /* ace dump a B Tree  */
#define showObj(_obj) bsAceDump(_obj,stdout,0,0)

void bsMinilibDumpB(OBJ obj, FILE *fil, Stack s) ;
 
    /* functions for local manipulations of branches and leaves */

/* simple routines, avoid opening the object in the calling routine */
KEY keyGetKey (KEY key, KEY tag) ;      /* get first key 1 to the right of tag */
KEY keyGetTag2Key (KEY key, KEY tag) ;  /* tag2 system, get first key 2 to the right of tag */
int   keyGetInt (KEY key, KEY tag) ;      /* get first key 1 to the right of tag */
float keyGetFloat (KEY key, KEY tag) ;      /* get first key 1 to the right of tag */
const char *keyGetText (KEY key, KEY tag) ;    /* idem for text, not parallelisable */
  
KEY keyFindTag (KEY key, KEY tag) ;
int keySetAddTag (KEYSET ks, KEY tag) ;  /* return number of keys where tag was actually added */
int keySetRemoveTag (KEYSET ks, KEY tag) ;   /* return number of keys where tag was  actually removed */

/* when you start a series of operations on a BS 
   never start with a _bsDown, _bsRight.
*/

BOOL bsFindTag     (OBJ obj, KEY tag) ;                 /* locates curr at tag */
#define bsFindKey(obj,tag,key)  bsGetData(obj, tag, key,0)     /* BOOL, locates curr at key */
BOOL bsFindKey2    (OBJ obj, KEY tag, KEY key, KEY *ptag2) ;
BOOL bsGetKey      (OBJ obj, KEY target, KEY *found) ;  /* Just keys */
BOOL bsGetKeyTags  (OBJ obj, KEY target, KEY *found) ;  /* Also tags */
BOOL bsGetData     (OBJ obj, KEY target, KEY type, void *x) ;
BOOL bsGetText     (OBJ obj, KEY target, char **cpp) ;  /* makes cp to point to a NON-REENTRANT ascii value of the cell */
BOOL bsGetArray    (OBJ obj, KEY target, Array units, int width) ; /* fills array(units, BSunit) */

BOOL bsPushObj     (OBJ obj) ;  /* Gets a #Class sub tree */
int  bsLocalClass  (OBJ obj) ;  /* current class, differs from class(key) after a bsPushObj */
BOOL bsFlatten   (OBJ obj, int n, Array units) ; /* flattens tree n deep */
KEY bsParentKey    (OBJ obj) ;  /* first key/class != 0 up from obj->curr */

KEYSET bsTagsInClass (int classe) ; /* Collects all tags in class */
BOOL bsIsTagInClass (int class, KEY tag) ; /* true if tag belongs to model of this class */
BOOL   bsIsTagInObj (OBJ obj, KEY key, KEY tag) ; /* looks for tag in obj  */
BOOL   bsIsClassInClass (int c1, int maClasse) ; /* looks for maclassein neighbours of c1 */
KEYSET bsKeySet (KEY key) ; /* Collects all keys in object bsCreate(key) */
KEYSET bsGetPath (int classe, KEY tag) ; /* returns path to given tag */

BOOL bsAddTag      (OBJ obj, KEY tag) ;
BOOL bsAddKey      (OBJ obj, KEY target, KEY new) ;
BOOL bsAddData     (OBJ obj, KEY target, KEY type, const void *data) ;
BOOL bsAddComment  (OBJ obj, const char* text, char type) ;
#define bsAddArray(obj,target,units,width) bsAddTruncatedArray (obj, target, units, width, width)
BOOL bsAddTruncatedArray (OBJ obj, KEY target, Array units, 
			  int arrWidth, int objWidth) ; /* reciprocal to GetArray */

void bsCoordShift (OBJ obj, KEY tag, float x, float dx, BOOL isAfter) ;
/* shift all nodes for which COORD is set in model with value greater
   or less than x (according to isAfter) by dx - for coordinate editing */
void bsCoordIndex (OBJ obj, int *index) ;
/* same but apply index[] to remap the coords (for assemblies) */

KEY  bsType        (OBJ obj, KEY target) ;  /* what if we add to curr ? */

BOOL bsRemove      (OBJ obj) ; /* removes curr and subtree to right */
BOOL bsRemoveTag      (OBJ obj, KEY target) ; /* locate and removes tag and subtree to right */
BOOL bsPrune	   (OBJ obj) ; /* prunes back unique branch to curr */

typedef void*      BSMARK ;	       /* private handle */
BSMARK  bsMark     (OBJ obj, BSMARK mark) ; /* preserve's current tree position */
BSMARK bsHandleMark (OBJ obj, BSMARK mark, AC_HANDLE hh) ; /* allocates mark on hh */
char *bsGetMarkAsString (BSMARK mark); /* for debugging */
void    bsGoto     (OBJ obj, BSMARK mark) ; /* returns to mark */
#define bsMarkFree(x)	messfree(x)

				/* the following act on the current node */
KEY  bsGetTimeStamp (OBJ obj) ;
mytime_t bsGetNodeTime(OBJ obj);
void bsSetTimeStamp (OBJ obj, KEY stamp) ;

int cacheSaveAll(void) ;
int cacheStatus (int *used, int *locked, int*known, int *modified, 
		 int *nc2rp, int *nc2wp, int *nccr, int *nccw) ;
 
void BSstatus (int *used, int *alloc, int *mem) ;	/* for status */
void BTstatus (int *used, int *alloc, int *mem) ;

void BSshutdown (void) ;
void BTshutdown (void) ;
void cacheShutDown (void) ;
void OBJShutDown (void) ;

void bsMakePathShutDown (void) ;

/******************************************************************/

/* globals defined in bssubs.c */

extern BOOL XREF_DISABLED_FOR_UPDATE ;
extern BOOL isTimeStamps;

#endif /* !ACEDB_BS_H */

/******************************************************************/
/******************************************************************/

