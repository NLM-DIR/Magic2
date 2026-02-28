/*  File: trace.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
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
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description:
 *              linked in with graphical acembly programs
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 22 10:39 1998 (fw)
 * Created: Mon Apr 18 17:41:24 1994 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: trace.c,v 1.29 2017/10/02 19:49:25 mieg Exp $ */
#include "acedb.h"
#include "acembly.h"
#include "freeout.h"
#include "parse.h"
#include "a.h"
#include "bs.h"

#include "display.h"
#include "../wh/menu.h"
#include "plot.h"
#include "key.h"
#include "sysclass.h"
#include "session.h"
#include "pick.h"
#include "cdna.h"
#include "map.h"		/* opaque decl of LOOK */
#include "basecall.h"		/* completes LookStruct -> TRACELOOK */
#include "fmap.h"		/* public fmap function headers */

/************************************************************/

static int nextError (LOOK look, BOOL down) ;

static LOOK traceLook = 0 ;
static Graph traceGraph = 0 ;
static float lastMag = 0, genex = 0 ;

/************************************************************/

static int BASE_COLOR [] = {A_COLOR, G_COLOR, C_COLOR, T_COLOR, N_COLOR} ;
static void traceDraw (LOOK look) ;
static BOOL fineTune(LOOK look, LANE **lane1p, LANE *lane) ;
static void traceSelect(LOOK look, int box) ;
static void traceMiddleDown (double x, double y) ;
static void traceLeftDown (int box, double xx, double yy) ;
static BOOL isLaneShown (LOOK look, LANE *lane) ;
static void traceDnaLane (LOOK look, MAP map, LANE *lane, 
			  float *offset, int min, int max) ;
static void traceDnaLaneErrors (LOOK look, MAP map, LANE *lane, 
				float *offset, int min, int max, BOOL foundSplice) ;
static int jackPotSearch (LOOK look, MAP map, LANE *lane, 
				    float *offset, int min, int max) ;
static void traceEdLaneBases (LOOK look, LANE *lane, float *offset,
			       int start, int stop, int type, BOOL isAligned) ;
static void laneShowExtrema(LOOK look, LANE *lane, float *offset, 
			    int min, int max) ;
BOOL traceGetLane (LOOK look, LANE *lane) ;
/* static void fixButton (void) ; */
static void traceRealignAll (void) ;
static void laneDoTag (LOOK look, LANE *lane, KEY key, char *buf) ;
static void newFixButton (void) ;
static void newBigFixButton (void) ;
static	BOOL traceLaneCrossing (Array lanes, LANE *lane) ;

#define   CHECKLANE { if (!lane->scf)  traceGetLane (look, lane) ; \
                      if (lane->scf > 1 && !lane->errArray) \
		       laneMakeErrArray (look, lane) ; }
#define   CHECKSEQ { if (lane->scf == 2) { baseCallGetSeq(lane) ; }}
#define   CHECKPOS { if (lane->scf == 3 && !lane->basePos) baseCallFakePos (lane) ; \
                     else if (lane->scf > 3 && !lane->basePos) baseCallGet(lane); }

static BOOL showSelect = TRUE ;
int laneGlobalOrder  (const void *a, const void *b) ;

#define SHOW_ALL 0 
#define HIDE_UP 1
#define HIDE_DOWN 2
#define SHOW_BEST_ERR 4
#define HIDE_NON_BEST 12
#define SHOW_CLIP 16
#define SHOW_SPLICER 32
#define HIDE_NON_FULL 64
#define HIDE_POOR 128

static int showDefaultMode = SHOW_BEST_ERR ;
static int showChosenMode = SHOW_ALL ;
static LANE* lastUnclippedLane = 0 ;
static KEY lastVectorClippedRead = 0 ;

#define NEXT_ERROR 1
#define NEXT_PROBLEM 2
#define NEXT_COMPARE 3
#define NEXT_CLIP 4
#define NEXT_HOLE 5
#define NEXT_TAG 6

static FREEOPT nextMenu[] = 
{ {6, "Next What"},
  {NEXT_PROBLEM, "Problem"},
  {NEXT_ERROR, "Difference"},
  {NEXT_CLIP, "Clip"},
  {NEXT_HOLE, "Hole"},
  {NEXT_TAG, "Tag"}, 
  {NEXT_COMPARE, "Comparison"}
} ;

static BOOL sortNeeded = TRUE ;
static int  dragFast = -1, dragBox = 0 ;

/************************************************************/

static void* TRACELOOK_ASSOC;	/* find the (TRACE)LOOK on the active graph */
static void* TRACELOOK_MAGIC;	/* verify the (TRACE)LOOK pointer */

#define TRACELOOKGET(name) \
  LOOK look ; \
  if (!graphAssFind (&TRACELOOK_ASSOC, &look)) \
    messcrash ("graph not found in %s", name) ; \
  if (look->magic != &TRACELOOK_MAGIC) \
    messcrash ("%s received a wrong pointer", name)

/************************************************************/
/***************** Editor     *******************************/
#ifndef NON_GRAPHIC
/************************************************************/

void cdnaStopHisto (void) ;

static void updateConsensus (LOOK look)
{ 
  if (!look->dnaKey || look->mode != SHOTGUN)
    return ;
  defCptForget (look->link, look->dnaKey);
  dnaStoreDestroy (look->dnaKey, dnaCopy(look->dna)) ;
  if (look->fMapLook)
    fMapPleaseRecompute (look->fMapLook) ;
}

/************************************************************/

static FREEOPT errorChoiceAll[] =
{ 
  {22, "Define"}, 

  {2, "A->T"},
  {3, "A->G"},
  {4, "A->C"},

  {1002, "T->A"},
  {1004, "T->G"},
  {1003, "T->C"},

  {5, "G->A"},
  {6, "G->T"},
  {7, "G->C"},

  {1006, "C->A"},
  {1005, "C->T"},
  {1007, "C->G"},

  {101, "Delete A"},
  {1101, "Delete T"},
  {102, "Delete G"},
  {1102, "Delete C"},

  {201, "Insert A"},
  {1201, "Insert T"},
  {202, "Insert G"},
  {1202, "Insert C"},

  {400, "Other"},
  {1400, "Other"}

} ;

static FREEOPT errorChoiceA[] =
{ 
  {9, "Define"}, 

  {2, "A->T"},
  {3, "A->G"},
  {4, "A->C"},

  {101, "Delete A"},

  {201, "Insert A"},
  {1201, "Insert T"},
  {202, "Insert G"},
  {1202, "Insert C"},
  {400, "Other"}
} ;

static FREEOPT errorChoiceT[] =
{ 
  {9, "Define"}, 

  {1002, "T->A"},
  {1004, "T->G"},
  {1003, "T->C"},

  {1101, "Delete T"},

  {201, "Insert A"},
  {1201, "Insert T"},
  {202, "Insert G"},
  {1202, "Insert C"},
  {400, "Other"}
} ;

static FREEOPT errorChoiceG[] =
{ 
  {9, "Define"}, 

  {5, "G->A"},
  {6, "G->T"},
  {7, "G->C"},

  {102, "Delete G"},

  {201, "Insert A"},
  {1201, "Insert T"},
  {202, "Insert G"},
  {1202, "Insert C"},
  {400, "Other"}
} ;

static FREEOPT errorChoiceC[] =
{ 
  {9, "Define"}, 

  {1006, "C->A"},
  {1005, "C->T"},
  {1007, "C->G"},

  {1102, "Delete C"},

  {201, "Insert A"},
  {1201, "Insert T"},
  {202, "Insert G"},
  {1202, "Insert C"},
  {400, "Other"}
} ;

static void doRegisterSequenceError (OBJ Cosmid, int x, int k, char *cp)
{
  bsAddData (Cosmid, str2tag ("Possible_error_at_base"),_Int,&x) ;
  freekey2text(k, errorChoiceAll) ;
  bsAddData (Cosmid, _bsRight, _Text,  freekey2text(k, errorChoiceAll)) ;
  bsAddData (Cosmid, _bsRight, _Text,  cp) ;
}

static void registerSequenceError (LOOK look, KEY type, int pos, char cBase)
{ 
  KEY gene = look->gene, cosmid = 0 ;
  OBJ Cosmid = 0 ;
  int a1, a2, x, k ; KEY key ;
  FREEOPT *errorChoice ;

  cosmid = keyGetKey (gene,str2tag ("Genomic_sequence")) ;
  Cosmid = bsUpdate (cosmid) ;
  if (!Cosmid || !bsFindKey (Cosmid, _Transcribed_gene, gene) ||
      !bsGetData (Cosmid, _bsRight,_Int, &a1) ||
      !bsGetData (Cosmid, _bsRight,_Int, &a2) )
    goto abort ;
  if (a1 < a2) x = pos + a1 - 1 ;
  else x = a1 - pos + 1 ;

  switch (cBase)
    {
    case A_: errorChoice = errorChoiceA ; break ;
    case T_: errorChoice = errorChoiceT ; break ;
    case G_: errorChoice = errorChoiceG ; break ;
    default: errorChoice = errorChoiceC ; break ;
    }

  if (graphSelect (&key, errorChoice))
    {
      k = key ;
      if (a1 > a2) 
	{ if (k < 1000) k += 1000 ; else k -= 1000 ; }
      if (messPrompt ("Please comment","","t"))
	{
	  char *comment = strnew (freeword (),0) ;
	  if (bsFindTag (Cosmid, str2tag("Genomic")))
	    doRegisterSequenceError (Cosmid, x, k, comment) ;
	  else if (bsFindTag (Cosmid, str2tag("Junction")))
	    {
	      KEYSET parts = 0 ;
	      KEY source = 0 ;

	      bsSave (Cosmid) ;
	      parts = queryKey (cosmid, ">Parts") ;
	      source = keyGetKey (cosmid, _Source) ;
	      if (source && keySetMax (parts))
		{
		  int u1, u2 ;
		  OBJ Source = bsCreate (source) ;
		  
		  if (Source && 
		      bsFindKey (Source, _Subsequence, cosmid) &&
		      bsGetData (Source, _bsRight,_Int, &u1) &&
		      bsGetData (Source, _bsRight,_Int, &u2))
		    {
		      int iPart = keySetMax (parts) ;
		      int xx = u1 < u2 ? u1 + x - 1 : u1 - x + 1 ; /* in source */

		      while (iPart--)
			{
			  KEY pp = keySet(parts, iPart) ;
			  int z1, z2 ;
			  
			  if (bsFindKey (Source, _Subsequence, pp) &&
			      bsGetData (Source, _bsRight,_Int, &z1) &&
			      bsGetData (Source, _bsRight,_Int, &z2))
			    {
			      int xxx = z1 < z2 ? xx - z1 + 1 : z1 - xx + 1 ;
			      int longueur = z1 < z2 ? z2 - z1 + 1 : z1 - z2 + 1 ;

			      if (xxx > 0 && xxx < longueur)
				{
				  OBJ Part = bsUpdate (pp) ;
				  if (Part)
				    doRegisterSequenceError (Part, xxx, k, comment) ;
				  bsSave (Part) ;
				}
			    }
			}
		      keySetDestroy (parts) ;
		    }
		  bsDestroy (Source) ;
		}
	      keySetDestroy (parts) ;
	    }
	  messfree (comment) ;
	}
    }
abort:
  bsSave (Cosmid) ;
}

/************************************************************/

static FREEOPT baseDoEditMenu[] = 
{
  {13,  "edit"},

  {'a', "A"},
  {'t', "T"},
  {'g', "G"},
  {'c', "C"},
  {'n', "n"},
  {'d', "Delete"},
  {'N', "Insert"},
 /*  modifying this to register hand clips
  {'x', "Clip Start"},
  {'y', "Clip End"},
  */
  {'J', "Clip Start"},
  {'j', "Clip End"},
  
  {'X', "Clip Start Insert"},
  {'Y', "Clip End Insert"}, 

  {'U', "Tag this base"},
  {'W', "Accept SNP"}
/*
  {1, "Clip Excellent"},
  {2, "Clip Good"},
  {3, "Clip Fair"},

  {'B', "Neural-net base call"},

  {14,"NN Excellent"},
  {13, "NN Good"},
  {12, "NN Fair"},
  {11, "NN Poor"},
  {10, "NN Bad"},

   {'w', "push down"},
   {'u', "push up"},
   {'s', "Show Trace"},
   {'h', "Hide Trace"},
   {'i', "more clip"},  // used by clipMore button 
   {400, "Close down jack pot"},
   {401, "Close up jack pot"},
   {'I', "Other Hint"},

*/
} ;

static FREEOPT baseDoEditcDNAMenu[] = 
{
  {21,  "edit"},
  
  {'a', "A"},
  {'t', "T"},
  {'g', "G"},
  {'c', "C"},
  {'n', "n"},
  {'d', "Delete"},
  {'P', "Insert"},
  {'W', "Accept SNP"},
  {'X', "Clip Start Insert:"},
  {'Y', "Clip End Insert"},  
  /* {'Z', "Clip End Quality"}, */
  { 300, "PolyA"},
  { 311, "Fake internal priming"},
  { 301, "SL1:cccaagtttgag"},
  { 302, "SL2:cagttactcaag"},
  { 303, "SL3:cagttaaccaag"},
  { 304, "SL0: capped clone"},
  { 305, "gccgtgctc"},
  { 306, "other motif"},
  /* { 320, "Favor normal"}, */
  { 321, "Favor delete"},
  { 322, "Favor insert"},
  /*
  { 312, "Stops matching above"},
  { 313, "Stops matching under"},
  */
  { 330, "Error in genomic"}
} ;
/*
static FREEOPT baseDoSuperEditMenu[] = 
{
  {6,  "Super Edit"},
  {'a', "A"},
  {'t', "T"},
  {'g', "G"},
  {'c', "C"},
  {'n', "n"},
  {'d', "Delete"}
} ;
*/
static FREEOPT baseDoInsertMenu[] = 
{
  {5,  "insert"},
  {'A', "Insert A"},
  {'T', "Insert T"},
  {'G', "Insert G"},
  {'C', "Insert C"},
  {'N', "Insert n"}
} ;

static FREEOPT *baseEditMenu = &baseDoEditMenu[0];
static void baseEdit1 (KEY key, int box) ;
static void (*baseEdit)(KEY key, int box) = baseEdit1 ;

static void saveHandClip (LOOK look, LANE *lane)
{
  OBJ obj = 0 ;
  int x ;
  
  if ((obj = bsUpdate (lane->key)))
    { 
      if (bsIsTagInObj (obj, lane->key, _Hand_Clipping))
	{
	  x = lane->handClipTop + 1 ;
	  bsAddData (obj, _Hand_Clipping, _Int, &x) ;
	  x = lane->handClipEnd ;
	  bsAddData (obj, _bsRight, _Int, &x) ;
	}
      bsSave (obj) ;
    }
}

static void baseDoEdit (LOOK look, LANE *lane, KEY key, 
			int iLong, int iShort, BOOL doRedraw)
{ int i, min, max, move = 0, pos = 0 ;
  char cc, *base, *cBase, *cp ;
  short *basePos ;
  float offSet, mag ;
  MAP map = look->map ;
  char buf [2] ;
  int oldErrMax = arrayMax(look->dna) ;
  BOOL isHand = FALSE ;
  HINT *hh = 0 ;
  static BOOL firstPass = TRUE ;

  if (look && firstPass)
    {
      firstPass = FALSE ;
      sessionAutoSave (60 * 30, 60 * 5) ;   /* save every 30 minutes, or after 5 minutes inactivity */
    }

  look->activeLane = lane ;
  if (key == 'P') 
    { /* solves a tiny interface problem when sequence is very dense */
      key = 'N' ;
      if (!lane->upSequence)
	{ iShort++ ; iLong++ ; }
    }
  if (!lane->base || iShort < 0 || iShort >= arrayMax(lane->base))
    return ;
  base = arrp (lane->base, iShort, char) ;
  if (iLong >= 0 && iLong < arrayMax(look->dna))
    cBase = arrp (look->dna, iLong, char) ;
  else
    return ;
  basePos = arrp (lane->basePos, iShort, short) ;
 lao:
  move = 0 ;
  switch (key)  /* non edit functions */
    {
    case 'J':
      isHand = TRUE ;
      key = 'x' ;
      goto lao ;
      break ;
    case 'j':
      isHand = TRUE ;
      key = 'y' ;
      goto lao ;
      break ;
    case 'h':
      if (lane->hide)
	return ;
      lane->hide = TRUE ;
      traceDraw (look) ;
      return ;
    case 'X': case 'Y':
      isHand = TRUE ;
    case 1: case 2: case 3: case 'x': case 'y': case 'Z': case 'i': 
    case 400: case 401:
    case 'B':
    case 10: case 11: case 12: case 13: case 14:
      switch (key)
	{
	case 'B': 
	  cc = arr (lane->base, pos, char) ;
	  nnBaseCall (lane, iShort) ;
	  if (cc != arr (lane->base, pos, char))
	    arr (lane->base, pos, char) |=  BC_HAND ;
	  break ;
	case 10: case 11: case 12: case 13: case 14:
	  nnQualityTrain (lane, iShort,key - 10) ;
	  return ;
	case 1:
	  baseCallUnclipLane (lane, 'E', &i) ;
	  break ;
	case 2:
	  baseCallUnclipLane (lane, 'G', &i) ;
	  break ;
	case 3:
	  baseCallUnclipLane (lane, 'F', &i) ;
	  break ;
	  
	case 'X':
	  i = iShort - lane->vectorTop ;
	  if (i*i > 36 && iShort > 40 &&
	      !messQuery(
			 messprintf("Do you really want to move the vector clipping by %d bases",
				    i))) return ;	
	  lane->vectorTop = iShort ;
	  if (lane->vectorTop < 0) lane->vectorTop = 0 ;
	  lastVectorClippedRead = lane->key ;
	  /* fall thru */
	case 'x': case 401:
	  i = iShort - lane->clipTop ;
	  lane->clipTop = iShort ;
	  if (lane->clipTop < 0) lane->clipTop = 0 ;
	  lane->x1 += i * (lane->upSequence ? -1 : 1) ;
	  break ;
	case 'y': case 'Y': case 'Z': case 'J':  case 400: case 305:
	  i = iShort - lane->clipEnd ;
	  if ( i + lane->clipEnd < lane->clipTop + 10)
	    i = lane->clipTop + 10 - lane->clipEnd ;
	  if ( i + lane->clipEnd > arrayMax (lane->dna))
	    i = arrayMax (lane->dna) - lane->clipEnd ;
	  if (key != 400 && key != 500 && i*i > 2500 &&
	      !messQuery(
			 messprintf("Do you really want to move the clipping by %d bases",
				    i))) return ;	
	  if (lane->upSequence)
	    {	
	      arrayDestroy (lane->errArray) ;
	      lane->clipEnd += i + 1; 
	      lane->x2 -= (i + 1) ; lane->x2 += -1 ;
	      /* key = 'l' ; */
	      if (look->fMapLook)
		fMapPleaseRecompute (look->fMapLook) ;
	      sortNeeded = TRUE ;
	    }
	  else
	    {
	      lane->clipEnd += i ;
	      lane->x2 += i ;
	    }
	  if (key == 'Y')
	    {
	      lane->vectorEnd = lane->clipEnd ;
	      lastVectorClippedRead = lane->key ;
	    }
	  else if (key == 'Z')
	    lane->qualityEnd = lane->clipEnd ;
	  break ;
	case 'i':
	  i = lane->clipEnd ;
	  if (lane->vectorEnd > 0 && 
	      lane->vectorEnd <= arrayMax(lane->base) &&
	      60 + iShort > lane->vectorEnd)
	    lane->clipEnd =  lane->vectorEnd ;
	  else if (60 + iShort > arrayMax(lane->base))
	    lane->clipEnd = arrayMax(lane->base) ;
	  else
	    lane->clipEnd += 60 ;
	  lane->x2 += (lane->clipEnd - i) * (lane->upSequence ? -1 : 1) ;
	  break ;
	}
      if (isHand)
	switch (key)
	  {
	  case 'X': 
	    lane->handClipTop = lane->clipTop ;
	    lane->handClipEnd = lane->vectorEnd ;
	    saveHandClip (look, lane) ;
	    break ;
	  case 'Y':
	    lane->handClipTop = lane->vectorTop ; 
	    lane->handClipEnd = lane->clipEnd ;
	    saveHandClip (look, lane) ;
	    break ;
	  default:
	    break ;
	  }
      lastUnclippedLane = lane ;
      key = 'l' ;
      if (lane->upSequence)
	{
	  if (look->fMapLook)
	    fMapPleaseRecompute (look->fMapLook) ;
	  sortNeeded = TRUE ;
	  goto lao ;
	}
      findXclipping (lane) ;
      goto saveedits ;
    case 'l':
      lane->x2 += lane->upSequence ? +1 : -1 ;
      dnaAlignRecale(look->dna, &(lane->x1), &(lane->x2),
		     lane->dna, lane->clipTop, lane->clipEnd - 1) ;
      lane->x2 -= lane->upSequence ? +1 : -1 ;
      findXclipping (lane) ;
      goto saveedits ;
    case 'I':   /* other hint */
      lane->getNewSpliceHint = lane->lastSpliceHint ? 2 : 0 ;
      goto redraw ;
    }

  if (lane->scf < 4)  /* trace not visible, only n are acceptables */
    switch (key)
      {
      case 'a': 
      case 't': 
      case 'g': 
      case 'c': 
      case 'W':
      case 'n': key = 'n' ; break ;
      case 'A': 
      case 'T': 
      case 'G': 
      case 'C':
      case 'N': key = 'N' ; break ; 
      }

  if (iLong >= 0)
    if (*cBase == N_ || !*cBase) /* remove N from consesnsus */
      switch (ace_lower(key))
	{
	case 'a': *cBase = A_ ; updateConsensus (look) ; break ;
	case 't': *cBase = T_ ; updateConsensus (look) ; break ;
	case 'g': *cBase = G_ ; updateConsensus (look) ; break ;
	case 'c': *cBase = C_ ; updateConsensus (look) ; break ;
	}

  if (lane->upSequence)
    switch (key)
      {
      case 'a': key = 't' ; break ;
      case 't': key = 'a' ; break ;
      case 'g': key = 'c' ; break ;
      case 'c': key = 'g' ; break ;
      case 'A': key = 'T' ; break ;
      case 'T': key = 'A' ; break ;
      case 'G': key = 'C' ; break ;
      case 'C': key = 'G' ; break ;
      case 'u': key = 'w' ; break ;
      case 'w': key = 'u' ; break ;
      }

  switch (key)
    {
    case 'a': *base &= ~BC_TAG ; *base = BC_LOW | BC_HAND | A_ ; break ;
    case 't': *base &= ~BC_TAG ; *base = BC_LOW | BC_HAND | T_ ; break ;
    case 'g': *base &= ~BC_TAG ; *base = BC_LOW | BC_HAND | G_ ; break ;
    case 'c': *base &= ~BC_TAG ; *base = BC_LOW | BC_HAND | C_ ; break ;
    case 'n': *base &= ~BC_TAG ; *base = BC_LOW | BC_HAND | N_ ; break ;
    case 'W': 
      {
	char mycBase = lane->upSequence ? complementBase(*cBase) : *cBase ;
	char snpBase = (*base | mycBase)  & 0xf ; /* accept wobble */
	*base = BC_LOW | BC_HAND | snpBase ; 
	break ;
      }
    case 'U': *base |= BC_TAG ; 
      { OBJ obj ;
	if ((obj = bsUpdate (lane->key)))
	  { buf [1] = 0 ;
	    buf [0] = dnaDecodeChar [*base & 0x0f] ;
	    cp = buf ; i = iShort + 1 ;
	    if (bsAddData (obj, _Significant_bases, _Int, &i)) 
	      bsAddData (obj, _bsRight, _Text, cp) ;
	    bsSave (obj) ;
	    if (look->fMapLook)
	      fMapPleaseRecompute (look->fMapLook) ;
	  }
      }
      break ;

    case 301: case 302: case 303: case 304: case 305: case 306:      
      {
	OBJ obj = 0 ;
	KEY tsl, oldTsl = 0 ;
	KEY _Manual_transpliced_to = str2tag ("Manual_transpliced_to") ;

	if ((obj = bsUpdate (lane->key)))
	  { 
	    switch (key)
	      {
	      case 301: 
		lexaddkey ("SL1", &tsl, _VMotif) ;
		break ;
	      case 302: 
		lexaddkey ("SL2", &tsl, _VMotif) ;
		break ;
	      case 303: 
		lexaddkey ("SL3", &tsl, _VMotif) ;
		break ;
	      case 304: 
		lexaddkey ("SL0", &tsl , _VMotif) ;
		break ;
	      case 305: 
		lexaddkey ("gccgtgctc", &tsl , _VMotif) ;
		break ;
	      case 306: 
		if (!messPrompt ("Motif de bas en haut","","t"))
		  { bsDestroy (obj) ; return ; }
		cp = freeword() ;
		lexaddkey (cp, &tsl , _VMotif) ;
		break ;
	      }

	    if (bsGetKey (obj,  _Transpliced_to, &oldTsl) &&
		strcmp (name(oldTsl),"SL0") &&
		strcmp (name(oldTsl),name(tsl)))
	      bsRemove (obj) ;
	    iShort++ ;
	    bsAddKey (obj,  _Manual_transpliced_to, tsl) ;
	    if (bsFindKey (obj, _Manual_transpliced_to, tsl))
	      bsAddData (obj,  _bsRight, _Int, &iShort) ; 
	    bsAddKey (obj,  _Transpliced_to, tsl) ;
	    if (bsFindKey (obj,  _Transpliced_to, tsl))
	      bsAddData (obj,  _bsRight, _Int, &iShort) ; 
	    iShort-- ;

	    bsSave (obj) ;
	    if (look->fMapLook)
	      fMapPleaseRecompute (look->fMapLook) ;

	    if (lane->upSequence) { key = 'y' ; iShort-- ; }
	    else { key = 'x' ; iShort++ ; }

	    goto lao ;
	  }
      }
      break ;

    case 300: case 311: case 312: case 313:
      {
	OBJ obj = 0 ;

	if ((obj = bsUpdate (lane->key)))
	  { 
	    switch (key)
	      {
	      case 300: 
		cBase = arrp (look->dna, iLong, char) ;
		while (*cBase == A_)
		  {
		    if (lane->upSequence)	
		      { cBase-- ; iLong-- ; iShort++ ;}
		    else
		      { cBase++ ; iLong++ ; iShort++ ;}
		  }
		if (lane->upSequence)
		  iShort += 0 ;
		else
		  iShort += 0 ;
		bsAddData (obj, _PolyA_after_base, _Int, &iShort) ;  
		bsAddData (obj, str2tag("Manual_PolyA"), _Int, &iShort) ;  
		break ;
	      case 311: 
		if (lane->upSequence)
		  iShort += 2 ;
		bsAddTag (obj, str2tag("Fake_internal_poly_A")) ;
		{ 
		  KEYSET ks2 = queryKey (lane->key, ">cDNA_clone ; >Read") ;
		  int iks = keySetMax(ks2) ;
		  OBJ obj2 = 0 ; KEY key2 ;
		  while (iks--)
		    { 
		      key2 = keySet(ks2, iks) ;
		      if (key2 != lane->key &&
			  (obj2 = bsUpdate (key2)))
			{
			  bsAddTag (obj2, str2tag("Fake_internal_poly_A")) ;
			  bsSave (obj2) ;
			}
		    }
		  keySetDestroy (ks2) ;
		}			  
		{ 
		  KEY clone = keyGetKey (lane->key, _cDNA_clone) ;
		  OBJ obj2 = clone ? bsUpdate (clone) : 0 ;
		  
		  if (obj2)
		    {
		      bsAddTag (obj2, str2tag("Internal_priming")) ;
		      bsAddTag (obj2, str2tag("Internal_priming_manual")) ;
		      bsSave (obj2) ;
		    }
		}
		break ;
	      case 312:
		  { 
		    KEY clone = keyGetKey (lane->key, _cDNA_clone) ;
		    OBJ obj2 = clone ? bsUpdate (clone) : 0 ;
		    
		    if (obj2)
		      {
			bsAddTag (obj2,  str2tag("Incomplete_5p_Match")) ;
			bsSave (obj2) ;
		      }
		  }
		i = -1 ;/* cut above */
		cp = messprintf("stops matching above %d",iShort) ;
		goto cutbelow ;
	      case 313:
		  { 
		    KEY clone = keyGetKey (lane->key, _cDNA_clone) ;
		    OBJ obj2 = clone ? bsUpdate (clone) : 0 ;
		    
		    if (obj2)
		      {
			bsAddTag (obj2,  str2tag("Incomplete_3p_Match")) ;
			bsSave (obj2) ;
		      }
		  }
		i = 1 ; /* cut below */
		cp = messprintf("stops matching below %d",iShort) ;
	      cutbelow:
		bsDestroy (obj) ;
		obj = bsUpdate (lane->clone) ;
		if (!obj)
		  return ;
		if (!messPrompt ("Stops Matching ?",cp,"w"))
		  { bsDestroy (obj) ; return ; }
		cp = freeword() ;
		iShort++ ;
		bsAddData (obj, _Problem, _Text, cp) ;
		iShort-- ;
		if (i == -1) { key = 'X' ; iShort++ ; }
		else { key = 'Y' ; iShort-- ; }
		break ;
	      }
	    bsSave (obj) ;
	    if (look->fMapLook)
	      fMapPleaseRecompute (look->fMapLook) ;
	    switch (key)
	      {
	      case 300:
		if (!lane->upSequence) { key = 'y' ; iShort++ ; } /* was key = 305 */
		else { key = 'x' ; iShort-- ; } /* -- since += 2 above */
		break ;
	      case 311:
		goto saveedits ;
		break ;
	      case 312:
		break ;
	      }
	    goto lao ;
	  }
      }
      break ;
    case 320: /* favor normal */
      lane->favorDelete = 1 ;
      break ;
    case 321: /* favor delete */
      lane->favorDelete = 2 ;
      key = 'd' ;
      goto lao ;
      break ;
    case 322: /* favor insert */
      lane->favorDelete = 3 ;
      key = 'N' ;
       if (!lane->upSequence)
	{ iShort++ ; iLong++ ; }
      goto lao ;
      break ;
    case 330:  /* Error in genomic */
      if (lane->upSequence)
	i = lane->x1 - iShort + lane->clipTop + 1 ;
      else
	i = lane->x1 + iShort - lane->clipTop + 1 ;
      { int j = lane->upSequence ? i - 1 : i - 1 ;
      i -= look->origin ;
      if (j >= 0 && look->dna && j < arrayMax(look->dna)) 
	registerSequenceError (look, key, i, arr(look->dna,j,char)) ;
      }
      return ;
    case 'd':
      *(basePos - 1) += (*basePos - *(basePos - 1)) / 3 ;
      *(basePos + 1) += (*basePos - *(basePos + 1)) / 3 ;
      i = arrayMax (lane->base) - iShort  - 1 ;
      if (i > 0) while (i--)
	{ *base = *(base + 1) ; base++ ;
	  *basePos = *(basePos + 1) ; basePos++ ;
	}
      arrayMax(lane->base)-- ;
      arrayMax(lane->basePos)-- ;
      if (iShort < lane->clipTop) lane->clipTop-- ;
      if (iShort < lane->clipEnd) lane->clipEnd-- ;
      if (iShort < lane->vectorTop) lane->vectorTop-- ;
      if (iShort < lane->vectorEnd) lane->vectorEnd-- ;
      if (iShort < lane->clipExtend) lane->clipExtend-- ;
      if (lane->clipEnd > arrayMax(lane->dna))
	  lane->clipEnd = arrayMax(lane->dna) ;
      pos = iShort ; /* for update in laneEdit Save */
      move = -1 ;
      break ;
    case 'u':
      *basePos -= (*basePos - *(basePos - 1)) / 3 ;
      break ;
    case 'w':
      *basePos -= (*basePos - *(basePos + 1)) / 3 ;
      break ;
    case 'A': case 'T': case 'G': case 'C': case 'N':
         /* extend */
      if (iShort < lane->clipTop) lane->clipTop++ ;
      if (iShort < lane->clipEnd) lane->clipEnd++ ;
      if (iShort < lane->clipExtend) lane->clipExtend++ ;
      if (iShort < lane->vectorTop) lane->vectorTop++ ;
      if (iShort < lane->vectorEnd) lane->vectorEnd++ ;

      base = arrayp(lane->base, arrayMax(lane->base) , char) ;
      basePos = arrayp(lane->basePos, arrayMax(lane->basePos) , short) ;
      i = arrayMax (lane->base) - iShort - 1 ;
      while (i--)
	{ *base = *(base - 1) ; base-- ;
	  *basePos = *(basePos - 1) ; basePos-- ;
	}

      base = arrp (lane->base, iShort, char) ;
      basePos = arrp (lane->basePos, iShort, short) ;

      *(basePos + 1) += ( *(basePos + 2) - *(basePos + 1))/ 3 ;
      *(basePos) -= ( *(basePos) - *(basePos - 1))/ 3 ;

      switch (key)
	{
	case 'A':
	  *base &= ~BC_TAG ; *base = BC_LOW | BC_HAND | A_ ;
	  break ;

	case 'T':
	  *base &= ~BC_TAG ; *base = BC_LOW | BC_HAND | T_ ;
	  break ;

	case 'G':
	  *base &= ~BC_TAG ; *base = BC_LOW | BC_HAND | G_ ;
	  break ;

	case 'C':
	  *base &= ~BC_TAG ; *base = BC_LOW | BC_HAND | C_ ;
	  break ;

	case 'N':
	  *base &= ~BC_TAG ; *base = BC_LOW | BC_HAND | N_ ;
	  break ;
	}
      pos = iShort ;		/* ? update seulement dans les ins/del ? */
      move = 1 ;		/* for update in laneEditSave */
      break ;
    }

 saveedits:

  laneEditSave (look, lane, iShort /* pos */, move) ;
  if (!doRedraw)
    return ;
  oldErrMax = arrayMax(lane->errArray) ;
  if (arrayMax(lane->errArray) > oldErrMax + 6) /* something stupid happenned */
    {
      if (look->mode == CDNA && key == 'x')
	;  /* just open intron */
      else
	{ /* realign */
	  key = 'l' ; goto lao ;
	}
    }
     
  if (!isWriteAccess ())	/* may occur is somebody else grabed it */
    sessionGainWriteAccess() ;	/* try to grab it */

  if (look->fMapLook)
    fMapPleaseRecompute (look->fMapLook) ;

  if (key == 'l')
    { traceDraw(look) ;
      return ;
    }

 redraw:
  i = arrayMax(look->hints) ;
  while (i--)
    {
      hh = arrp(look->hints,i, HINT) ;
      if (lane == hh->lane)
	hh->lane = 0 ;
    } 
  /* for some reason these boxes sometimes get wrong */
  if (lane->laneBaseCallBox < graphLastBox())
    graphBoxClear (lane->laneBaseCallBox) ;
  if (lane->edLaneBaseBox < graphLastBox())
    graphBoxClear (lane->edLaneBaseBox) ;

  map->centre -= lane->dy ;
  mag = map->mag ;
  map->mag *= lane->ddy ;
  lane->laneBaseCallBox = graphBoxStart () ;
  offSet = lane->laneBaseCallBoxOffSet ;
  min = lane->laneBaseCallBoxMin ;
  max = lane->laneBaseCallBoxMax ;

  traceDnaLane (look, look->map, lane, &offSet, min, max) ;
  graphBoxEnd () ;
  map->centre += lane->dy ;
  map->mag = mag ;

  offSet = lane->edLaneBaseBoxOffSet ;

  map->mag *= look->edMag ;
  if (lane->scf >= 3)
    traceEdLaneBases (look, lane, &offSet, look->edMin,
		      look->edMax, lane->scf, FALSE) ;
  map->mag = mag ;

  graphBoxDraw (lane->laneBaseCallBox, -1, -1) ;
  graphBoxDraw (lane->edLaneBaseBox, -1, -1) ;
}

/*****************/

static void baseEdit1 (KEY key, int box)
{ HINT* hh ;
  TRACELOOKGET("baseedit1") ;
  
  if (box < arrayMax(look->baseBoxes) &&
      (hh = arrp(look->baseBoxes,box, HINT)))	  
    {
      if (hh->iLong == 1 && hh->y > 0)
	hh->iLong = GRAPH2MAP (look->map, hh->y) ;	
      if (hh->iLong < 0 ||  hh->iLong >= arrayMax(look->dna))
	hh->iLong = 1 ;
      
      baseDoEdit (look, hh->lane, key, hh->iLong, hh->iShort, TRUE) ;
    }
}

/********************/

static void hintAccept (LOOK look, HINT *hh)
{ if (hh && hh->lane)
    {
      if (hh->lane->hide)
	{ hh->lane->hide = FALSE ;
	  traceDraw (look) ;
	}
      else
	{ 
	  if (hh->isSpliceHint == 2)
	    hh->lane->getNewSpliceHint = 1 ; /* stabilize these hints */
	  baseDoEdit (look, hh->lane, hh->key, hh->iLong, hh->iShort, TRUE) ;
	  hh->lane = 0 ; /* do not reuse */
	}
    }
}

/************************************************************/

static void doAcceptZoneEdits (BOOL ok, BOOL doDraw, BOOL isV, LANE *myLane, float y1, float y2)
{
  int i ;
  HINT *hh ;
  TRACELOOKGET("doAcceptZoneEdits") ;

  i = arrayMax(look->hints) ;
  while (i--)
    {
      hh = arrp(look->hints,i, HINT) ;
      if (!hh || !hh->lane || !hh->iLong || hh->lane->hide || 
	  hh->isSpliceHint ||   /* no zone editing of splice hints */
	  hh->lane->scf < 4)  /* no zone editing if trace not seen  */
	continue ;
      if (myLane && (hh->lane != myLane))
	continue ;
      if (!hh->lane->upSequence && !hh->lane->previousEnd && hh->iShort < hh->lane->clipTop + 6 &&
	  hh->lane->vectorTop < 5)
	{
	  if (!isV)
	    continue ; /* do not edit too close to vector */
	  else
	    {
	      int xv =-1, yv = -1, j ;
	      HINT *h2 ;
	      LANE *lane = hh->lane ;

	      for (j = 0,  h2 = arrp(look->hints,j, HINT) ;
		   j < arrayMax(look->hints) ; j++, h2++)
		{
		  if (h2 && h2->lane &&
		      h2->lane == lane && 
		      h2->iShort < 80 &&
		      h2->iShort < h2->lane->clipTop + 8 &&
		      xv < h2->iShort )
		    { xv = h2->iShort ; yv = h2->iLong ; h2->lane = 0 ; }
		}
	      if (xv > 0) 
		baseDoEdit (look, lane, 'X', yv + 1, xv + 1, FALSE) ;
	      continue ;
	    }
	}
      switch (hh->key)
	{
	case 'A': case 'T': case 'G': case 'C': case 'N': if (!ok) hh->key = 'N' ; break ;
	case 'a': case 't': case 'g': case 'c':  if (!ok) hh->key = 'n' ; break ;
	case 'd': break ;
	default: continue ;
	} 
      if (!y1 && !y2)
	{
	  y1 = look->topMargin + 2 ; 
	  y2 = look->graphHeight - 2 ;
	}
      if (hh->key == 'A' || hh->key == 'T' || hh->key == 'G' || hh->key == 'C' || hh->key == 'N' || 
	  !(BC_HAND & arr (hh->lane->base, hh->iShort, char))) /* do not Re edit in  zone mode */
	{
	  if (hh->y > y1 && hh->y + .9 < y2 
	      /* && (BC_LOW & arr (hh->lane->baseCall, hh->iShort, char)) 25/6/2002, try to avoid high peaks */
	      )
	    baseDoEdit (look, hh->lane, hh->key, hh->iLong, hh->iShort, FALSE) ;
	}
      hh->lane = 0 ; /* no recursion */
    } 
  showSelect = myLane ? FALSE : TRUE ;   /* FALSE ; si on veut stabiliser et voir les editions */
  if (doDraw) traceDraw(look) ;
  showSelect = TRUE ;
}

/*****************/

static void acceptZoneEdits (void)
{ doAcceptZoneEdits (TRUE, TRUE, FALSE, 0, 0, 0) ; }

/*****************/
 
static void nifyZoneEdits (void)
{ doAcceptZoneEdits (FALSE, TRUE, FALSE, 0, 0, 0) ; }

/*****************/
 
static void vectorifyZoneEdits (void)
{ doAcceptZoneEdits (FALSE, TRUE, TRUE, 0, 0, 0) ; }

/*****************/
 
#ifdef APPARENTLY_NOT_NEEDED
static void undoLastVectorClip (void)
{ 
  OBJ obj ;
  BOOL done = FALSE ;
  
  if ((obj = bsUpdate(lastVectorClippedRead)))
    {
      if (bsFindTag (obj, _Vector_Clipping))
	{ 
	  bsRemove (obj) ;
	  done = TRUE ;
	}
      bsSave (obj) ;
    }
  lastVectorClippedRead = 0 ;
  if (done)
    newFixButton () ;
  else
    messout ("No previous vector cliping, sorry") ;
}
#endif /* APPARENTLY_NOT_NEEDED */

/************************************************************/

static void traceDiscard (LOOK look, KEY key)
{
  int ii ; 

  KEY clone = 0, gene = look->gene ;
  LANE *lane ;

  ii = arrayMax (look->lanes) ;
  if (!ii) return ; 
  
  while (ii--)
    { 
      lane = arrp (look->lanes, ii, LANE) ;
      if (key == lane->key) lane->scf = 1 ;
    } 
  showSelect = TRUE ;
      
  clone = keyGetKey (key, _cDNA_clone) ;
  cDNARemoveCloneFromGene (clone, gene) ;

  if (look->fMapLook)
    fMapPleaseRecompute (look->fMapLook) ;
}

/************************************************************/

static void laneTagger (KEY key, int box)
{
  extern void fMapTraceSuppress (KEY key, KEY tag, BOOL keep); /* fmaptrace.c */
  LANE *lane ;
  int i, ii, c1, c2, v1, v2 ;
  char *cp ;
  OBJ obj = 0 ;
  KEY kk ;
  TRACELOOKGET("laneTagger") ;
 
  ii = arrayExists(look->laneShown) ? 
    arrayMax (look->laneShown) : 0 ;
  if (!ii) return ; 
  
  while (ii--)
    { lane = arr (look->laneShown, ii, LANE*) ;
      if (box == lane->boxName)
	{ switch (key)
	    {
	    case 'e':
	      lanePlotDerivee (look, lane) ;
	      return ;
	    case 'p':
	      lanePlotPeriodicite (look, lane) ;
	      return ;
	    case 'B':
	      nnLaneTrain (look, lane) ;
	      return ;
	    case 'q':
	      if (look->mode == CDNA) 
		{ messout ("Open just a single est trace from the Est menu in the sequence window re-base call") ;
		return ;
		}
	      if (!messQuery(messprintf("%s\n%s\n\n%s",
   "I will irreversibly replace the present base call by the acembly basecall",
   "corresponding to the colored horizontal lines",
   "Should I proceed (y/n)"))) return ;
	      CHECKLANE ;
	      if (lane->scf < 3) break ;
	      baseCallRedoLaneBaseCall (lane, 0) ;
	      baseDoEdit(look, lane, 'l', 1,1, TRUE) ;
	      if (arrayMax(look->lanes) == 1) /* single read, then redraw all */
		graphDestroy () ;
	      return ;
	    case 's':
	    case 'S':
	      lane->scf = 1 ; /* will suppress further display */
	      if (look->fMapLook)
		{ fMapPleaseRecompute (look->fMapLook) ;
		  if (graphExists(look->fMapGraph))
		    { Graph old = graphActive() ;
		      graphActivate (look->fMapGraph) ;
		      fMapTraceSuppress (lane->key, _Assembled_into,key == 's' ? TRUE : FALSE) ;
		      graphActivate (old) ;
		    }
		}
	      traceDraw (look) ;
	      return ;
	    case 'h': /* hide */
	      lane->hide = TRUE ;
	      showSelect = FALSE ;
	      traceDraw (look) ;
	      break ;
	    case 'd':
	      traceDiscard (look, lane->key) ;
	      break ;	      
	    case 'z':
	      doAcceptZoneEdits(TRUE, TRUE, FALSE, lane, 0, 0) ;  /* see also lane->boxClose */
	      return ;
	      break ;
	    case 'w':
	      display(lane->key, 0, TREE) ;
	      return ;
	    case 'V':
	      lane->vectorTop = lane->vectorEnd = 0 ;
	      {
		Stack s = stackCreate (1000) ;
		int level = freeOutSetStack (s) ;
		i = trackVector (0, lane->key, TRUE) ;
		freeOutClose (level) ;
		parseBuffer (stackText (s, 0), 0) ;
		stackDestroy (s) ;
	      }
	      if (!i)
		return ;
	      c1 = c2 = v1 = v2 = 0 ;
	      obj = bsCreate (lane->key) ;
	      if (!obj || !bsGetData (obj, _Clipping, _Int, &c1) ||
		  !bsGetData (obj, _bsRight, _Int, &c2))
		i = 0 ;
	      if (!obj || !bsGetData (obj, _Vector_clipping, _Int, &v1) ||
		  !bsGetData (obj, _bsRight, _Int, &v2))
		i = 0 ;
	      bsDestroy (obj) ;
	      if (! i)
		return ; 
	      i = c1 - lane->clipTop ;
	      lane->clipTop += i ;
	      lane->x1 += i * (lane->upSequence ? -1 : 1) ;
	      if (c2 < v2) c2 = v2 ; /* open a max */
	      i = c2 - lane->clipEnd ;
	      if ( i + lane->clipEnd < lane->clipTop + 10)
		i = lane->clipTop + 10 - lane->clipEnd ;
	      if ( i + lane->clipEnd > arrayMax (lane->dna))
		i = arrayMax (lane->dna) - lane->clipEnd ;
	      lane->clipEnd += i ;
	      lane->x2 += i * (lane->upSequence ? -1 : 1) ;
	      lane->vectorTop = v1 - 1 ;
	      lane->vectorEnd = v2 ;
	      lastUnclippedLane = lane ;
	      sortNeeded = TRUE ;
	      if (lane->upSequence && look->fMapLook)
		fMapPleaseRecompute (look->fMapLook) ;
	      lane->x2 += lane->upSequence ? +1 : -1 ;
	      dnaAlignRecale(look->dna, &(lane->x1), &(lane->x2),
			     lane->dna, lane->clipTop, lane->clipEnd - 1) ;
	      lane->x2 -= lane->upSequence ? +1 : -1 ;
	      findXclipping (lane) ;
	      laneEditSave (look, lane, 0, 0) ;

	      if (!isWriteAccess ())	/* may occur is somebody else grabed it */
		sessionGainWriteAccess() ; /* try to grab it */
	      break ;
	    case 1000: /* unedit */
	      {
		OBJ obj = bsUpdate (lane->key) ;
		KEY kk = 0 ;
		
		if (obj && 
		    bsGetKey (obj, _SCF_Position, &kk) &&
		    bsFindTag (obj, _SCF_Position))
		  {
		    bsRemove (obj) ;
		    arrayKill (kk) ;
		  }
		if (obj && 
		    bsGetKey (obj, _BaseCall, &kk) &&
		    bsFindTag (obj, _BaseCall))
		  {
		    bsRemove (obj) ;
		    arrayKill (kk) ;
		  }
		bsSave (obj) ;
		arrayDestroy (lane->errArray) ;
		arrayDestroy (lane->dna) ;
		arrayDestroy (lane->base) ;
		arrayDestroy (lane->basePos) ;
		arrayDestroy (lane->baseCall) ;
		arrayDestroy (lane->tags) ;
		lane->scf = 0 ;
		CHECKLANE ; CHECKSEQ ; CHECKPOS ;
		baseDoEdit (look, lane, 'l', 1,1, TRUE) ;
		newFixButton () ;
	      }
	      break ;
	    case 999: /* basecall */
	      {
		OBJ obj = 0 ;
		KEY kk = 0 ;
		
		if (lane->scf < 3) break ;
		if (!messQuery(messprintf("%s\n%s\n\n%s",
					  "I will irreversibly replace the present base call by the acembly basecall",
					  "corresponding to the colored horizontal lines",
					  "Should I proceed (y/n)"))) break ;
		
		CHECKLANE ;
		obj = bsUpdate (lane->key) ;
		if (obj && 
		    bsGetKey (obj, _SCF_Position, &kk) &&
		    bsFindTag (obj, _SCF_Position))
		  {
		    bsRemove (obj) ;
		    arrayKill (kk) ;
		  }
		if (obj && 
		    bsGetKey (obj, _BaseCall, &kk) &&
		    bsFindTag (obj, _BaseCall))
		  {
		    bsRemove (obj) ;
		    arrayKill (kk) ;
		  }
		if (obj && 
		    bsGetKey (obj, _DNA, &kk) &&
		    bsFindTag (obj, _DNA))
		  {
		    bsRemove (obj) ;
		    arrayKill (kk) ;
		  }
		if (obj && 
		    bsGetKey (obj, str2tag("OriginalBaseCall"), &kk) &&
		    bsFindTag (obj, str2tag("OriginalBaseCall")))
		  {
		    bsRemove (obj) ;
		    arrayKill (kk) ;
		  }
		bsSave (obj) ;

		baseCallRedoBaseCall (lane->key, 0) ;
		baseDoEdit(look, lane, 'l', 1,1, TRUE) ;
 
		newFixButton () ;
	      }
	      break ;
	    case 1001: /* select */
	      if (!keySetExists(look->selectedKs))
		look->selectedKs = keySetCreate() ;
	      kk = lane->clone ? lane->clone : lane->key ;
	      keySetInsert (look->selectedKs, kk) ;
	      break ;
	    case 1002: /* unselect */ 
	      if (keySetExists(look->selectedKs))
		{
		  kk = lane->clone ? lane->clone : lane->key ;
		  keySetRemove (look->selectedKs, kk) ;
		}
	      break ;
	    default:
	      if (!messPrompt ("Comments ?:", "", "t"))
		return ;
	      cp = freeword () ;
	      laneDoTag (look, lane, key, cp) ;
	      break ;
	    }
	  fMapPleaseRecompute (look->fMapLook) ;
	  traceDraw (look) ;
	  break ;
	}
    }
}

/************/

static FREEOPT laneTagMenu [] = {
  {7, "TagMenu"},
  {'w', "Who am I"},
  {'t', "Tag"}, 
  {'s', "Move to new contig"},
  {'S', "Discard"},
/*
  {'h', "Hide Lane"},
  {'e', "energie de la derivee"},
  {'B', "Neural-Net Training"},
*/
  {'p', "periodicite"},
  {'V', "Search Vector"},
  {'q', "Base Call"}  
    } ;

static FREEOPT laneTagcDNAMenu [] = {
  { 7, "TagMenu"},
  {'w', "Who am I"},
  {'V', "Search Vector"},
  {'d', "Discard from this gene"},
  {1000, "Unedit"},
  {1001, "Select"},
  {1002, "Unselect"}, 
  {999, "Base Call"}  
    } ;

static FREEOPT traceTagMenu [] = {
  {4, "TagMenu"},
  {'t', "Comment"},
  {'c', "Compression"},
  {'a', "Ambiguity"},
  {'r', "Resolved"}
   } ;

/************************************************************/
/***************** Destruction ******** *********************/

static void traceDoDestroy (LOOK look)
{ int i ;
  LANE *lane ;

  if (look->lanes &&
      (i = arrayMax(look->lanes)))
    while (i--)
      { lane = arrp(look->lanes, i, LANE) ;
	laneDestroy (lane) ;
      }
  arrayDestroy (look->lanes) ;
  keySetDestroy (look->modifiedTraces) ;
  arrayDestroy(look->dna) ;
  arrayDestroy(look->dnaR) ;
  arrayDestroy(look->hints) ;
  arrayDestroy(look->laneShown) ;
  arrayDestroy(look->baseBoxes) ;
  arrayDestroy(look->hits) ;
  arrayDestroy(look->selectedKs) ;
}

static void traceDestroy (void)
{ LOOK look = traceLook ;

  if (look && graphExists (traceGraph))
    { 
      mapDestroy (look->map) ;
      graphAssRemove (&MAP2LOOK_ASSOC) ; /* detach this graph from
						 a map (or vice versa) */
      
      arrayDestroy (look->box2seg) ;
      
      traceDoDestroy (look) ;
    }
  traceGraph = 0 ;
  traceLook = 0 ;
  if (look)
   { look->magic = 0 ;
     messfree (look) ;
   }
}

/************************************************************/
/***************** Mouse and Kbd ****************************/
/************************************************************/

static void traceSelect(LOOK look, int box)
{
  HINT *hh ;
  
  if (box < arrayMax(look->hints))
    {
      hh = arrp(look->hints,box, HINT) ;
      if (hh)
	{
	  hintAccept (look, arrp(look->hints,box, HINT)) ;
	  return ;
	}
    }

  if (box < arrayMax(look->baseBoxes) &&
      (hh = arrp(look->baseBoxes,box, HINT)))
    look->activeLane = hh->lane ;

}

static void tracePick (int box, double x, double y) 
{ LANE *lane = 0 ;
  int i ;
  TRACELOOKGET("tracePick") ;

  if (!box) return ;
  i = arrayMax(look->lanes) ;
  while (i--)
    { lane = arrp(look->lanes, i, LANE) ;
      if (!isLaneShown (look, lane))
	continue ;
      if (box == lane->edLaneBaseBox)
	{ 
	  lane->hide = !lane->hide ;
	  showSelect = FALSE ;
	  traceDraw (look) ;
	  return ;
	}
      if (!lane->isShown)
	continue ;
      if (box == lane->boxMoreClip)
	{ baseDoEdit (look, lane, 'i', 0, lane->clipEnd - 1, TRUE) ;
	  return ;
	}
      else if (box == lane->boxLessClip)
	{ if (lane->clipEnd > lane->clipTop + 22)
	    baseDoEdit (look, lane, 'y', 0, lane->clipEnd - 10, TRUE) ;
	  return ;
	}
      else if (box == lane->boxLessClipTop)
	{
	  /*lane->showTraceBeforeClipTop = TRUE ;
	    traceDraw (look) ;
	  */
	  if (lane->clipTop > 0)
	    baseDoEdit (look, lane, 'x', 0, 
			lane->vectorTop > 0 && lane->clipTop > lane->vectorTop ? lane->vectorTop : 
			(lane->clipTop > 20 ? lane->clipTop - 20 : 0) , 
			TRUE) ;
	  return ;
	}
      else if (box == lane->boxClose)
	{    
	  laneTagger ('z', lane->boxName) ;
	  /*
	  lane->hide = TRUE ;
	  showSelect = FALSE ;
	  traceDraw (look) ;
	  return ;
	  */
	  return ;

	}
      else if (box == lane->boxName)
	{ display(lane->key, 0, TREE) ;
	  return ;
	}
      else if (box == lane->boxNewSpliceHint)
	{ 
	  baseDoEdit (look, lane, 'I', 0, 0, TRUE) ;
	  return ;
	}
    }

  if (box == look->geneBox)
    {
      look->map->centre = look->geneMin + 
	y * (look->geneMax - look->geneMin) / (.8 * ((float)look->graphHeight - look->topMargin)) ;
      (look->map->draw) () ;
      return ;
    }

  if (box == look->map->thumb.box)
    { 
      int  cc = look->map->centre - look->origin ;
      
      if (look->mode == CDNA)
	fMapDrawFromTrace (look->fMapLook, look->gene, cc, 1) ;
      else
	{
	  i = arrayMax(look->lanes) ;
	  while (i--)
	    { lane = arrp(look->lanes, i, LANE) ;
	    if (cc >= lane->x1 && cc < lane->x2)
	      break ;
	    }
	  /* take smallest  value, because seg in fmap are always downwards */
	  i = lane->x1 < lane->x2 ? lane->x1 : lane->x2 ;
	  
	  if (graphActivate (look->fMapGraph) &&
	      look->fMapLook)
	    fMapDrawFromTrace (look->fMapLook, lane->key, 
			       (int)( look->map->centre - i - look->origin), 1) ;
	  else
	    display (look->link, lane->key, 0) ;
	}
      graphActivate (look->graph) ;
    }
  else
    { if (box > look->minLiveBox) 
	traceSelect (look, box) ;
      look->activeBox = box ;
    }
}

/*********************************************************/
/****************** get sequences ************************/
/*********************************************************/
  /* lane->seq == the ABI file, 
   * lane->dna == the edited-assembled sequence, clipTop...
   */
/*
extern BOOL dnaAlignForceMatch (Array longDna, int x1, int x2, 
				Array shortDna, int y1, int y2, 
				int*topp, int *endp, int *sensp) ;
*/

/*********************************************************/

static int lanePlotOrder  (const void *a, const void *b)
{ const LANE* la = *(const LANE**)a,  *lb = *(const LANE**) b ;

  BOOL
    ra = la->upSequence ,
    rb = lb->upSequence ;
  int 
    xa = la->x1 ,
    xb = lb->x1 ;
  if (ra && !rb)
    return -1 ;
  if (rb && !ra)
    return 1 ;

  return xa - xb ; /* best traces will be on the outside */
    /* else favor closest clipTop
  if (ra)  //thus rb 
    return xa - xb ; 
  else
    return xb - xb ; */
}

/***************************************************************************/

static BOOL traceIsPoorLane (KEY tg, KEY est)
{
  static Array units = 0 ;
  static KEYSET ks = 0 ;
  static KEY myTg = 0, myTg2 = 0 ;
  BSunit *uu ;
  int ii, jj, j2 ;
  OBJ Mrna = 0 ;

  if (!est) /* reinitialise the need to search */
    { myTg = myTg2 = 0 ; return FALSE ; } 

  /*
    if (keyFindTag (est, _Transpliced_to))
    return FALSE ;
  */
  if (tg != myTg)
    {
      myTg = tg ;
      keySetDestroy (ks) ;
      ks = queryKey (tg, ">mrna ; {>Best_available_clone ; >read} $| {>mrna_covered_by} $| {>cds_covered_by}") ;
    }
  
  if (keySetFind (ks, est, 0))
    return FALSE ;
  
  j2 = 0 ;
  if (tg != myTg2)
    {
      KEYSET mrnas = queryKey (tg, ">mrna") ;
      KEY mrna ;

      myTg2 = tg ;
      units = arrayReCreate (units, 120, BSunit) ;
      for (jj = 0 ; jj < keySetMax(mrnas) ; jj++)
	{
	  mrna = keySet (mrnas, jj) ;
	  if ((Mrna = bsCreate (mrna)))
	    {
	      bsGetArray (Mrna, str2tag ("Tiling"), units, 5) ;
	      for (ii = 0 ; ii < arrayMax(units) ; ii += 5)
		{
		  uu = arrp (units, ii, BSunit) ;
		  keySetInsert (ks, uu[2].k) ;
		}
	      j2 += arrayMax(units) ;
	      bsDestroy (Mrna) ;
	    }
	} 
      keySetDestroy (mrnas) ;
    } 
  
  if (j2 && keySetFind (ks, est, 0))
    return FALSE ;
  
  return TRUE ;
}

/***************************************************************************/

static BOOL traceMakeLanes (LOOK look, int sens)
{ Array aa = arrayCreate (1000, BSunit), aa1 = arrayCreate (30, BSunit) ;
  int i, ii = 10000, j, x1, x2 , ct, ce, max = look->map->max ;
  LANE *lane ; 
  BSunit *u ;
  OBJ obj = 0 ;
  KEY key ;

  if ((obj = bsCreate (look->key)))
    { 
      if (bsGetArray (obj, _Assembled_from, aa, 5)) ;
      else   /* a read, artificially create one lane */
	{ array(aa, 0, BSunit).k = look->key ;
	  x1 = 1 ; x2 = arrayMax(look->dna) ;
	  array(aa, 1, BSunit).i = x1 ;
	  array(aa, 2, BSunit).i = x2 ;
	  array(aa, 3, BSunit).i = x1 ;
	  array(aa, 4, BSunit).i = x2 ;
	}
      if (bsGetArray (obj, _Aligned, aa1, 5)) /* pour rester en periodicite 5 */
	{ ii = i = arrayMax(aa) ; 
	  for (j = 0 ; j < arrayMax (aa1) ; i++, j++)
	    array (aa, i, BSunit) = array (aa1, j, BSunit) ;
	}
      bsDestroy (obj) ;
    }

  look->lanes = arrayCreate (12, LANE) ;
  if (sens > 0)
    { for (i = 0 , j = 0 ; i < arrayMax(aa) ; i += 5)
	{ u = arrp (aa, i, BSunit) ;
	  key = u[0].k ;
	  x1 = u[1].i ;
	  x2 = u[2].i ;
	  ct = u[3].i ;
	  ce = u[4].i ;
	  /* x3 = (x1 < x2 ? 1 : -1) * 50 + x2 ; */
	  lane = arrayp (look->lanes, j++, LANE) ;
	  lane->key = key ;
	  lane->hide = look->hide ;
	  if (x1 < x2)
	    { lane->upSequence = FALSE ;
	      lane->x1 = x1 - 1 ;
	      lane->x2 = x2 ;
	      lane->x3 = x2 ; /* later += (lane->clipExtend - lane->clipEnd) */
	    }
	  else
	    { lane->upSequence = TRUE ;
	      lane->x1 = x1 - 1 ;
	      lane->x2 = x2 - 2 ;
	      lane->x3 = x2 - 2 ;
	    }
	  lane->clipTop = ct - 1 ;
	  lane->clipEnd = ce ;
	  lane->scf = 0 ; /* will search it when needed */
	  if (i >= ii)
	    lane->isAligned = TRUE ;
	  if (look->mode ==  CDNA) 
	    {
	      lane->clone = keyGetKey (lane->key,_cDNA_clone) ;
	      lane->isPoor = traceIsPoorLane (look->key, lane->key) ;
	    }
	}
    }
  else /* contig a l'envers dans le link */
    { for (i = 0 , j = 0 ; i < arrayMax(aa) ; i+= 5)
	{ u = arrp (aa, i, BSunit) ;
	  key = u[0].k ;
	  x1 = u[1].i ;
	  x2 = u[2].i ;
	  ct = u[3].i ;
	  ce = u[4].i ;
	  /* x3 = (x1 < x2 ? 1 : -1) * 50 + x2 ; */
	  lane = arrayp (look->lanes,j++, LANE) ;
	  lane->key = key ;
	  lane->hide = look->hide ;
	  if (x1 > x2)
	    { lane->upSequence = FALSE ;
	      lane->x1 = max - x1 ;
	      lane->x2 = max - x2 + 1 ;
	      lane->x3 = max - x2 + 1 ;
	    }
	  else
	    { lane->upSequence = TRUE ;
	      lane->x1 = max - x1 ;
	      lane->x2 = max - x2 - 1 ;
	      lane->x3 = max - x2 - 1 ;
	    }
	  lane->clipTop = ct - 1 ;
	  lane->clipEnd = ce ;
	  lane->scf = 0 ; /* will search it when needed */
	  if (i >= ii)
	    lane->isAligned = TRUE ;
	  if (look->mode ==  CDNA) 
	    {
	      lane->clone = keyGetKey (lane->key,_cDNA_clone) ;
	      lane->isPoor = traceIsPoorLane (look->key, lane->key) ;
	    }
	}
    }
  arrayDestroy (aa) ;
  arrayDestroy (aa1) ;
  arraySort(look->lanes, laneGlobalOrder) ;
  sortNeeded = FALSE ;
  return arrayMax (look->lanes) ;
}

/***************************************************************************/

static BOOL traceMakeHitLanes (LOOK look, int sens)
{ Array aa = arrayCreate (1000, BSunit), aa1 = arrayCreate (30, BSunit) ;
  int j, j1, x1, x2 , c3, ct, ce, max = look->map->max ;
  LANE *lane ; 
  HIT *hh ;

  look->lanes = arrayCreate (12, LANE) ;
  if (sens > 0)
    { 
      for (j = j1 = 0 ; j1 < arrayMax(look->hits) ; j1++)
	{ 
	  hh = arrp (look->hits, j1, HIT) ;
	  x1 = hh->a1 ;
	  x2 = hh->a2 ;
	  ct = hh->x1 ;
	  ce = hh->x2 ;
	  lane = arrayp (look->lanes, j++, LANE) ;
	  lane->key = hh->est ;
	  lane->hide = look->hide ;
	  if (ct < ce)
	    { lane->upSequence = FALSE ;
	      lane->x1 = x1 - 1 ;
	      lane->x2 = x2 ;
	      lane->x3 = x2 ; /* later += (lane->clipExtend - lane->clipEnd) */
	      if (j1 > 0 &&
		  hh->est == (hh-1)->est)
		lane->previousEnd = (hh-1)->a2 ;
	      if (j1 < arrayMax(look->hits) - 1 && 
		  hh->est == (hh+1)->est)
		lane->nextTop = (hh+1)->a1 - 1 ;
	    }
	  else
	    { lane->upSequence = TRUE ;
	      lane->x1 = x2 - 1 ;
	      lane->x2 = x1 - 2 ;
	      lane->x3 = x1 - 2 ;
	      c3 = ce ; ce = ct ; ct = c3 ;
	      if (j1 > 0 &&
		  hh->est == (hh-1)->est)
		lane->nextTop = (hh-1)->a2 - 1 ;
	      if (j1 < arrayMax(look->hits) - 1 && 
		  hh->est == (hh+1)->est)
		lane->previousEnd = (hh+1)->a1 - 2 ;
	    }
	  lane->clipTop = ct - 1 ;
	  lane->clipEnd = ce ;
	  lane->scf = 0 ; /* will search it when needed */
	  lane->clone = keyGetKey (lane->key,_cDNA_clone) ; 
	  lane->isPoor = traceIsPoorLane (look->key, lane->key) ;
	}
    }
  else /* contig a l'envers dans le link */
    { 
      for (j = j1 = 0 ; j1 < arrayMax(look->hits) ; j1++)
	{ 
	  hh = arrp (look->hits, j1, HIT) ;
	  x1 = hh->a1 ;
	  x2 = hh->a2 ;
	  ct = hh->x1 ;
	  ce = hh->x2 ;
	  lane = arrayp (look->lanes,j++, LANE) ;
	  lane->key = hh->est ;
	  lane->hide = look->hide ;
	  if (x1 > x2)
	    { lane->upSequence = FALSE ;
	      lane->x1 = max - x1 ;
	      lane->x2 = max - x2 + 1 ;
	      lane->x3 = max - x2 + 1 ;
	      if (j1 > 0 &&
		  hh->est == (hh-1)->est)
		lane->previousEnd = max - (hh-1)->a2 + 1;
	      if (j1 < arrayMax(look->hits) - 1 && 
		  hh->est == (hh+1)->est)
		lane->nextTop = max - (hh+1)->a1 ;
	    }
	  else
	    { lane->upSequence = TRUE ;
	      lane->x1 = max - x1 ;
	      lane->x2 = max - x2 - 1 ;
	      lane->x3 = max - x2 - 1 ;
	      if (j1 > 0 &&
		  hh->est == (hh-1)->est)
		lane->nextTop =  max -  (hh-1)->a2 - 1;
	      if (j1 < arrayMax(look->hits) - 1 && 
		  hh->est == (hh+1)->est)
		lane->previousEnd = max - (hh+1)->a1 ;
	    }
	  lane->clipTop = ct - 1 ;
	  lane->clipEnd = ce ;
	  lane->scf = 0 ; /* will search it when needed */
	  lane->clone = keyGetKey (lane->key,_cDNA_clone) ;
	  lane->isPoor = traceIsPoorLane (look->key, lane->key) ;
	}
    }
  arrayDestroy (aa) ;
  arrayDestroy (aa1) ;
  arraySort(look->lanes, laneGlobalOrder) ;
  sortNeeded = FALSE ;
  return arrayMax (look->lanes) ;
}

/*********************************************************/
/****************** columns       ************************/
/*********************************************************/

static void traceConsensus (LOOK look, float *offset)
{
  int i, j, box ;
  HINT *hh ;
  LANE *lane ;
  BOOL old = FALSE ;
  char *cp, buf[2] ;
  float x = *offset, oldW, y, yy, yy1, yy2 ;
  int start, stop ;
  MAP map = look->map ;
  
  start = GRAPH2MAP (map, look->topMargin) ;
  stop = GRAPH2MAP (map, look->graphHeight) ;
  buf[1] = 0 ;
  if (start >= arrayMax(look->dna)) 
    start = arrayMax(look->dna) - 1 ;
  if (start < 0) 
    start = 0 ;
  if (stop > arrayMax(look->dna)) 
    stop = arrayMax(look->dna) ;
  i = stop - start ;
  cp = arrp (look->dna, start, char) ;
  j = start ;
  if (i>0) while (i--)
    {
      buf[0] = dnaDecodeChar[(int)*cp++] ;
      y = MAP2GRAPH(map, j) ; /* no zero */
      box = graphBoxStart () ;
      hh = arrayp (look->baseBoxes, box, HINT) ;
      hh->lane = (LANE *)(-1) ;
      hh->iShort = 1 ;
      hh->iLong = j ;
      hh->y = y ;
      hh->key = 0 ;
      graphText(buf, x, y) ;
      
      /*  graphColor (WHITE) ;
	  graphRectangle (x - .3, y, x, y) ;
	  */
      graphBoxEnd() ;
      /*
	graphBoxFreeMenu (box, baseEdit, baseDoSuperEditMenu) ;
	crash le code pour le moment
	*/
      j++ ; yy = MAP2GRAPH(map, j) - .1;
      if (map->mag > 0)
	{ yy1 = y + 1 ; yy2 = yy - .1 ; }
      else
	{ yy1 = yy + 1 ; yy2 = y - .1 ; }
      if (TRUE || yy2 > yy1)
	{
	  box = graphBoxStart () ;
	  hh = arrayp (look->baseBoxes, box, HINT) ;
	  hh->lane = (LANE *)(-1) ;
	  hh->iShort = 1 ;
	  hh->iLong = j ;
	  hh->key = 0 ;
	  hh->y = (yy1 + yy2) / 2 ;
	  graphRectangle (*offset - 3, yy1, *offset + 1, yy2) ;
	  graphBoxEnd() ;
	  graphBoxDraw (box, WHITE, WHITE) ;
	  graphBoxFreeMenu (box, baseEdit, baseDoInsertMenu) ;
	}
    }
  *offset += 1 ; /* after consensus */
  
  
  if (look->consensus) 
    {
      for (i = 0, old = FALSE ; i < arrayMax (look->laneShown) ; i++) 
	{
	  lane = arr (look->laneShown ,i, LANE*) ;
	  if (!lane->isAligned)
	    continue ;
	  CHECKSEQ ; 
	  if (!old)
	    { old = TRUE ;
	    *offset += 2 ; /* betweeen up and down */
	    }
	  traceEdLaneBases (look, lane, offset, look->edMin,
			    look->edMax, 3, TRUE) ;
	  lane->nerr = -1 ; /* don t use as best trace */
	}
    }
  
 /* line after consensus */
  oldW = graphLinewidth(.3) ;
  *offset += 1 ;
  graphLine (*offset, look->topMargin, *offset, look->graphHeight) ;
  *offset += 6 ; /* 1.5 ; */
  graphLinewidth(oldW) ;
}

/***************************************************************************************/
/***************************************************************************************/

static void showDoIt (LOOK look)
{ 
  LANE *lane ;
  int i,j ;
  static KEY _Complete_CDS_of = 0 ;

  if (!_Complete_CDS_of) _Complete_CDS_of = str2tag ("Complete_CDS_of") ;

  switch (look->hide)
    {
    case SHOW_BEST_ERR:
    case HIDE_NON_BEST: break ;
    case SHOW_ALL:
      if (look->lanes && (i = arrayMax(look->lanes)))
	{
	  lane = arrp (look->lanes, 0, LANE) - 1 ;
	  while (lane++, i--)
	    lane->hide = FALSE ;
	}
      break ;
    }
  switch (look->hideUp)
    {
    case HIDE_UP:
    case HIDE_DOWN:
    case HIDE_POOR:
      if (look->lanes && (i = arrayMax(look->lanes)))
	{ 
	  lane = arrp (look->lanes, 0, LANE) - 1 ;
	  j = 0 ;
	  while (lane++, i--)
	    switch (look->hideUp)
	      {
	      case HIDE_UP: lane->hide = lane->upSequence ; break ;
	      case HIDE_DOWN: lane->hide = !lane->upSequence ; break ;
	      case HIDE_POOR: lane->hide = lane->isPoor ; if(lane->isPoor) j++ ; break ;
	      }
	  if (j) printf("Hide %d poor lanes\n",j) ;
	}
      break ;
    case HIDE_NON_FULL:
      if (look->lanes && (i = arrayMax(look->lanes)))
	{
	  lane = arrp (look->lanes, 0, LANE) - 1 ;
	  while (lane++, i--)
	    lane->hide = keyFindTag (lane->clone, _Complete_CDS_of) ? FALSE : TRUE ; 
	}
      break ;
    default:
      break ;
    }
  showSelect = TRUE ; 
  traceDraw (look) ;
}

/***********/

static void showChooser (KEY k, int box)
{ TRACELOOKGET("showChooser") ;

   switch (k)
    {
    case HIDE_UP:
    case HIDE_DOWN:
    case HIDE_NON_FULL:
    case HIDE_POOR:
      look->hideUp = k ;
      break ;
    case SHOW_ALL:
      if (look->hideUp)
	look->hide = SHOW_BEST_ERR ; /* to provoke redraw */
      look->hideUp = 0 ;
    default:
      if (look->hide == k)
	return ;
      look->hide = showChosenMode = k ;
    }
  showDoIt (look) ;
}

/***********/

static void showButton (void) /* toggles SHOW_BEST_ERR */
{ TRACELOOKGET("showButton") ;

  if (look->hide != showChosenMode)
    look->hide = showChosenMode  ;
  else if (showDefaultMode != showChosenMode)
    look->hide = showDefaultMode ;
  else
    { 
      if (look->hide != SHOW_ALL)
	look->hide = SHOW_ALL ;
      else
	look->hide = SHOW_BEST_ERR ;
    }
  if (look->hideUp == HIDE_POOR)
    look->hideUp = 0 ;
  showDoIt (look) ;
}

/****************************************************/
/****************************************************/

static BOOL traceMakeLaneTags (LANE *lane)
{ OBJ obj ;
  Array units = 0 ;

  if (lane->hasTag == 1)
    return FALSE ;
  lane->hasTag = 1 ;
  obj = bsCreate (lane->key) ;
  if (!obj)
    return FALSE ;
  units = arrayReCreate (lane->tags, 12, BSunit) ;
  if (bsFindTag (obj, _Assembly_tags) &&
      bsFlatten (obj, 4, units) &&
      arrayMax (units))
    lane->tags = units ;
  else
    { arrayDestroy (units) ;
      bsDestroy (obj) ;
      return FALSE ;
    }
  bsDestroy (obj) ;

  lane->hasTag = 2 ;
  return TRUE ;
}

static BOOL traceDrawDnaTag (LOOK look, MAP map, LANE *lane,
			 float *offset) 
{ int i , box, col, x, x1, x2, u1, u2, min, max ;
  float y1, y2 ;
  static Array units = 0 ;
  BSunit *u ; char *cp ;
  BOOL found = FALSE ;

  if (!traceMakeLaneTags(lane)) 
    return FALSE ;
  
  min = look->edMin ; max = look->edMax ;
  u1 = lane->x1 ; u2 = lane->x2 ;
  if (u1 > u2) { i = u1 ; u1 = u2 ; u2 = i ; }
  if (min < u1)
    min = u1 ;
  if (max > u2)
    max = u2 ;
  if (min >= max) 
    return FALSE ;
  
  units = lane->tags ;
  for ( i = 0 ; i < arrayMax(units) ; i += 4)
    { u = arrp(units,i,BSunit) ;
      cp = u[0].s ;
      x1 = u[1].i ;
      x2 = u[2].i ;
      if (lane->upSequence)
	{ x1 = u2 - x1 + lane->clipTop ;
	  x2 = u2 - x2 + lane->clipTop ;
	}
      else
	{ x1 = u1 + x1 - lane->clipTop ;
	  x2 = u1 + x2 - lane->clipTop ;
	}
      if (x1 > x2)
	{ x = x1 ; x1 = x2 ; x2 = x ; }
      if ( x1 > max || x2 < min)
	continue ;
      if (x1 < min) x1 = min ;
      if (x2 > max) x2 = max ;
      y1 = MAP2GRAPH (map, x1) ;
      y2 = MAP2GRAPH (map, x2) ;
      box = graphBoxStart () ;
      graphRectangle (*offset, y1, *offset + .8, y2) ;
      graphBoxEnd () ;
      if (!strcmp (cp, "Compression"))
	col = LIGHTRED ;
      else if (!strcmp (cp, "Ambiguity"))
	col = LIGHTRED ;
      else 
	col = LIGHTGREEN ;
      graphBoxDraw (box, BLACK, col) ;
      found = TRUE ;
    }
  if (found) *offset += 1 ;
  return found ;
}

/***********/

static void traceDnaLane (LOOK look, MAP map, LANE *lane, float *offset, int min, int max)
{ int n = 0, top, fin, foundSplice ;
  float y = 0, yy = 0,  yy1 = 0, yy2 = 0 ;
  int  nBase , box , color, bgcolor, vv ;
  char c, *base, *cp ; short *basePos ;
  char buf[2] ;
  /*  BOOL showCoord = look->coord ; */
  BOOL compl = map->mag < 0 ? TRUE : FALSE ;
  HINT *hh ;
  KEY key ;

  CHECKSEQ ;  CHECKPOS ;
  if (!lane->base || !lane->basePos)
    return ;

  *offset += 1.5;  /* shift right (for errors sake) */
  buf[1] = 0 ;
  nBase = arrayMax (lane->basePos) ;
  if (!nBase)
    return ;

  if (min < 0)
    min = 0 ;
  if (max > lane->maxPos)
    max = lane->maxPos ;
  
  if (min < lane->xClipTop )
    { int box = graphBoxStart() ;
      graphRectangle (*offset, MAP2GRAPH(map, min),
		      *offset + 1, MAP2GRAPH(map, lane->xClipTop)) ;
      graphBoxEnd () ;
      color = YELLOW ;
      if (look->mode == CDNA && lane->previousEnd)
	color = PALEMAGENTA ;
      graphBoxDraw (box, BLACK, color) ;
    }

  if (lane->vectorTop > 0 && lane->vectorTop < arrayMax(lane->basePos))
    { vv = (arr(lane->basePos, lane->vectorTop, short) +
	arr(lane->basePos, lane->vectorTop - 1 , short)) / 2 ;
      if (min < vv )
	{ int box = graphBoxStart() ;
	  graphRectangle (*offset, MAP2GRAPH(map, min),
			  *offset + 1, MAP2GRAPH(map, vv)) ;
	  graphBoxEnd () ;
	  graphBoxDraw (box, BLACK, LIGHTBLUE) ;
	}
    }

  if (max > lane->xClipEnd )
    { int box = graphBoxStart() ;
      graphRectangle (*offset, MAP2GRAPH(map, lane->xClipEnd),
		      *offset + 1, MAP2GRAPH(map, max)) ;
      graphBoxEnd () ;
      color = YELLOW ; 
      if (look->mode == CDNA && lane->nextTop)
	color = PALEMAGENTA ;
      graphBoxDraw (box, BLACK, color) ;
    }

  if (lane->vectorEnd > 0 && lane->vectorEnd < arrayMax(lane->basePos))
    { vv = (arr(lane->basePos, lane->vectorEnd, short) +
	arr(lane->basePos, lane->vectorEnd + 1 , short)) / 2 ;
      if (max > vv )
	{ int box = graphBoxStart() ;
	  graphRectangle (*offset, MAP2GRAPH(map, vv),
			  *offset + 1, compl ? look->topMargin : look->graphHeight) ;
	  graphBoxEnd () ;
	  graphBoxDraw (box, BLACK, ORANGE) ;
	}
    }

  base = arrp (lane->base, 0, char) ;
  basePos = arrp (lane->basePos, 0, short) ;

  while(nBase > 1 && *basePos < min)
    { basePos++ ; base++ ; nBase-- ; n++ ;
    }
  top = n + 1 ; /* Plato */
  graphTextFormat (BOLD) ; 
  lane->minBaseBox = graphBoxStart() ;
  graphBoxEnd () ;
  lane->minBase = base - arrp (lane->base, 0, char) ;
  color = WHITE ;
  while(nBase && *basePos <= max)  
    { n++ ;
      box = graphBoxStart() ;
      hh = arrayp (look->baseBoxes, box, HINT) ;
      hh->lane = lane ;
      hh->iShort = base - arrp (lane->base, 0, char) ;
      hh->iLong = 1 ;
      hh->key = 0 ;
      y = MAP2GRAPH(map, *basePos) ;
      hh->y = y ;
      c = *base & 0x0f ;
      *buf = dnaDecodeChar[(int)(compl ? complementBase(c) : c)] ;
      if (*base & BC_LOW) 
	*buf = ace_lower (*buf) ;
      else
	*buf = ace_upper (*buf) ;
      if (c == N_) *buf = '?' ;

      graphText("   ", *offset - 1.8, y) ; /* make the box  bit bigger */
      graphText(buf, *offset, y) ; 
      graphBoxEnd() ;

      color = BLACK ; 
      switch (ace_lower(*buf))
	{
	case 'a': color = A_COLOR ; break ;
	case 't': color = T_COLOR ; break ;
	case 'g': color = G_COLOR ; break ;
	case 'c': color = C_COLOR ; break ;
	case '?': color = MAGENTA ; break ;
	}

      bgcolor = WHITE ;
      if (*base & BC_ACE) bgcolor = LIGHTGREEN ;
      if (*base & BC_HAND) bgcolor = YELLOW ;
      if (*base & BC_TAG)
	{ bgcolor = color ;
	  color = BLACK ;
	} 
      graphBoxDraw(box, color, bgcolor) ;
      graphBoxFreeMenu (box, baseEdit, baseEditMenu) ;
      if (nBase > 1)
	{
	  yy = MAP2GRAPH(map, *(basePos + 1)) ;
	  if (map->mag > 0)
	    { yy1 = y + 1 ; yy2 = yy - .1 ; }
	  else
	    { yy1 = yy + 1 ; yy2 = y - .1 ; }
	  if (hh->iShort >= lane->clipTop - 1 &&
	      hh->iShort < lane->clipEnd &&
	      yy2 > yy1)
	    { box = graphBoxStart () ;
	      hh = arrayp (look->baseBoxes, box, HINT) ;
	      hh->lane = lane ;
	      hh->iShort = base + 1 - arrp (lane->base, 0, char) ;
	      hh->iLong = 1 ;
	      hh->key = 0 ;
	      hh->y = (yy1 + yy2) / 2 ;
	      graphRectangle (*offset - 3, yy1, *offset + .8 , yy2) ;
	      graphBoxEnd() ;
	      graphBoxDraw (box, WHITE, WHITE) ;
	      graphBoxFreeMenu (box, baseEdit, baseDoInsertMenu) ;
	    }
	}
      basePos++ ; base++ ; nBase-- ;
    }
  lane->maxBase = fin = n ;

  if (!lane->hide)
    { int x = *offset - .8 ;
      float oldW, y ;
      lane->xHideButton = x ;

      lane->boxName = box = graphBoxStart () ;
      y = look->topMargin ;
      oldW = graphLinewidth(.3) ;

      if (compl)
	{ graphLine (x, y, x + 2 , y) ;
	  graphLine (x, y, x + 1 , y - 1) ;
	  graphLine (x + 1 , y - 1, x + 2 , y) ;
	}
      else
	{ graphLine (x, y  - 1, x + 2 , y - 1) ;
	  graphLine (x, y - 1, x + 1 , y) ;
	  graphLine (x + 1 , y, x + 2 , y - 1) ;
	}
      graphLinewidth(oldW) ;
      graphBoxEnd () ;
      graphBoxFreeMenu (box, laneTagger, 
	look->mode == CDNA ? laneTagcDNAMenu : laneTagMenu) ;
      color = TRANSPARENT ;
      if (lane->color) color = lane->color ;
      key = lane->clone ? lane->clone : lane->key ;
      if (keySetExists(look->selectedKs) &&
	  keySetFind (look->selectedKs, key ,0))
	color = RED ;
      graphBoxDraw (box, BLACK, color) ;

      box = graphBoxStart () ;
      cp = messprintf("%d", compl ? fin : top) ;
      graphText (cp, x - strlen(cp) + 1.9, y - 2.3) ;
      graphBoxEnd () ;
      if (lane == look->activeLane)
	graphBoxDraw (box, BLACK, LIGHTBLUE) ;
    }
  else
    lane->boxName = 0 ;
 
  graphTextFormat (PLAIN_FORMAT) ; 
  /*  was:  *offset +=  showCoord ? 1.5 + strlen(messprintf("%d",n)) : 1.8 ; */
  foundSplice = look->mode == CDNA ?
    jackPotSearch (look, map, lane, offset, min, max) : 0 ;
  lane->getNewSpliceHint = 0 ;
  if (foundSplice < 2)
    traceDnaLaneErrors (look, map, lane, offset, min, max, foundSplice) ;
  *offset += 1.95 ;  /* was 1.8 */
}     

/***********/
#ifdef JUNK
static FREEOPT jackPotMenu [] = {
  {4, "Jack-Pot"},
  {'a', "Accept"},
  {'c', "Close"},
  {'r', "Reduce match"},
  {'i', "Increase match"},
} ;

/************/

static void jackPotAction (KEY key, int box)
{ LANE *lane, *lane1 ;
  int i, ii ;
  TRACELOOKGET("jackPotAction") ;

 
  ii = arrayExists(look->laneShown) ? 
    arrayMax (look->laneShown) : 0 ;
  if (!ii) return ; 
  
  while (ii--)
    { 
      lane = arr (look->laneShown, ii, LANE*) ;
      if (box == lane->boxNewSpliceHint)
	{ 
	  switch (key)
	    {
	    case 'a': /* accept */
	      lane1 = 0 ;
		{
		  if (lane->jackPotSens == 1 && lane->nextTop)
		    {
		      i = arrayMax (look->lanes) ; 
		      while (i--)
			{
			  lane1 = arrp (look->lanes, i, LANE) ;
			  if (lane != lane1 && lane->key == lane1->key &&
			      lane1->x1 == lane->nextTop)
			    break ;
			}
		      if (i == -1) lane1 = 0 ;
		    }
		  else if (lane->jackPotSens == -1 && lane->previousEnd)
		    {
		      i = arrayMax (look->lanes) ; 
		      while (i--)
			{
			  lane1 = arrp (look->lanes, i, LANE) ;
			  if (lane != lane1 && lane->key == lane1->key &&
			      lane1->x2 == lane->previousEnd) /* was lane1->x1 == lanePreviousEnd + 1 */
			    break ;
			}
		      if (i == -1) lane1 = 0 ;
		    }

		  if (lane->jackPotSens == 1)
		    {
		      if (lane1) 
			{
			  baseDoEdit (look, lane1, 'x', lane->jackPotOriginLong + lane->jackPotJump, lane->jackPotOriginShort + 1, FALSE) ;
			  lane->nextTop = lane->jackPotOriginLong + lane->jackPotJump ;
			  lane1->previousEnd = lane->jackPotOriginLong ;
			}
		      baseDoEdit (look, lane, 'y', lane->jackPotOriginLong, lane->jackPotOriginShort, TRUE) ;
		    }
		  else
		    {
		      if (lane1) 
			{ 
			  baseDoEdit (look, lane1, 'y', lane->jackPotOriginLong - lane->jackPotJump, lane->jackPotOriginShort - 1, FALSE) ;
			  lane->previousEnd = lane->jackPotOriginLong - lane->jackPotJump ;
			  lane1->nextTop = lane->jackPotOriginLong ;
			}

		      baseDoEdit (look, lane, 'x', lane->jackPotOriginLong, lane->jackPotOriginShort, TRUE) ;
		    }
		}
	      break ;
	    case 'c': /* close */ 
	      if (lane->jackPotSens == 1)
		baseDoEdit (look, lane, 400, lane->jackPotOriginLong, lane->jackPotOriginShort,  TRUE) ;
	      else if (lane->jackPotSens == -1)
		baseDoEdit (look, lane, 401, lane->jackPotOriginLong, lane->jackPotOriginShort,  TRUE) ;
	      break ; 
	    case 'r':  /* reduce */
	      if (lane->jackPotMatchLength > 4)
		lane->jackPotMatchLength-- ;  
	      baseDoEdit (look, lane, 'I', 0, 0, TRUE) ;
	    case 'i':  /* increase */
	      if (lane->jackPotMatchLength < 9)
		lane->jackPotMatchLength++ ;  
	      baseDoEdit (look, lane, 'I', 0, 0, TRUE) ;
	    }
	  break ;
	}
    }
}
#endif
/***********/

static BOOL jackPotHint (LOOK look, MAP map, LANE *lane, float *offset, 
				    int min, int max, int iShort, int iLong, int jump, 
				    int sens, BOOL justTheButton)
{ 
  float y ;
  int n = 0, nBase , box , color, sensShort, origin = look->origin ;
  char *base ; short *basePos , bp ;
  char c, cc, *cp, *cq, buf[2] ;
  HINT *hh ;
  BOOL done = FALSE ;
  CHECKSEQ ;
  if (lane->scf < 3)
    return FALSE ;
  CHECKPOS ;
  buf[1] = 0 ;
  nBase = arrayMax (lane->base) ;

  if (min < 0)
    min = 0 ;
  if (max > lane->maxPos)
    max = lane->maxPos ;

  base = arrp (lane->base, 0, char) ;
  basePos = arrp (lane->basePos, 0, short) ;

  while(*basePos < min  && nBase)
    { basePos++ ; base++ ; nBase-- ; n++ ;
    }

  lane->minBase = base - arrp (lane->base, 0, char) ;
  if (lane->upSequence)
    {
    }
  else
    { 
    }
  *offset -= 3.2 ;

  sensShort = lane->upSequence ? -sens : sens ;
  bp = (arr (lane->basePos, iShort - sensShort, short) + 
	arr (lane->basePos, iShort, short))/2 ;
  if (bp > min && bp < max)
    {
      y = MAP2GRAPH(map, bp) + (1.0 - sens) / 4.0 ;
      if (!justTheButton)
	{
	  graphLine (*offset -1, y, *offset+3, y) ;
	  graphLine (*offset -1, y +.3 * sens, *offset+3, y +.3 * sens) ;
	  if (sens == -1) y -= .8 ;
	  if (sens == -1)
	    graphText (messprintf("%d", iLong - origin + (lane->upSequence ? +1 : +1)), *offset-1, y - sens*1.9) ;
	  else
	    graphText (messprintf("%d", iLong - origin + (lane->upSequence ? +1 : +1)), *offset-1, y - sens*1.9) ;
	  if (sens == -1)
	    graphText (messprintf("%d", iLong - origin + jump * sens + (lane->upSequence ? +2 : +2)), *offset-1, y - sens*1.0) ;
	  else
	    graphText (messprintf("%d", iLong - origin + jump * sens + (lane->upSequence ? +0 : +0)), *offset-1, y - sens*1.0) ;
	  if (sens == -1) y += .8 ;
	}
      if ((!lane->upSequence && sens == 1) || (lane->upSequence && sens == -1))
	{
	  lane->jackPotOriginLong = iLong - 1 ;
	  lane->jackPotOriginShort = iShort - 1 ;
	  lane->jackPotSens = 1 ;
	}
      else
	{
	  lane->jackPotOriginLong = iLong + 1 ;
	  lane->jackPotOriginShort = iShort + 1 ;
	  lane->jackPotSens = -1 ;
	}
      lane->jackPotJump = jump ;
      lane->boxNewSpliceHint = box = graphBoxStart() ;
      graphCircle (*offset, y - 3*sens, 1) ;
      graphCircle (*offset, y - 3*sens, .8) ;
      graphBoxEnd () ;
      /*
	graphBoxFreeMenu (box, jackPotAction,jackPotMenu) ;
	ce truc est tout faux, a verifier avant de reactiver 
       */
    }
  if (justTheButton)
    { lane->lastSpliceHint = 0 ; *offset += 3.2 ; return FALSE ; }
  lane->lastSpliceHint = sens * jump ;

  iLong += jump * sens ;
  color = WHITE ;
  iLong -= sens ; iShort -= sensShort ; 
  cp = arrp (look->dna, iLong, char) ; cq = arrp (lane->dna, iShort, char) ;
  while (cp += sens , cq += sensShort, iLong += sens, iShort += sensShort,
	 iShort >= 0 &&
	 iShort < arrayMax(lane->basePos) &&
	 iLong >= 0 &&
	 iLong < arrayMax(look->dna) 
	 )	 
    { bp = arr (lane->basePos, iShort, short) ;
      if (bp > max)
	break ;
      y = MAP2GRAPH(map, bp) ;
      if (y < look->topMargin || y > look->graphHeight) continue ;
      done = TRUE ;
      cc = !lane->upSequence ? *cq : complementBase(*cq) ;
      c = dnaDecodeChar[(int)(*cp)] ;
      c = ace_lower(c) ;
      box = graphBoxStart () ;
      color = WHITE ;
      if (*cp != cc)
	{
	  hh = arrayp (look->hints, box, HINT) ;
	  hh->lane = lane ;
	  hh->iShort = iShort ;
	  hh->iLong = iLong ;
	  hh->y = y ;
	  hh->key = c ;
	  hh->isSpliceHint = 2 ;
	  buf[0] = ' ' ;
	  color = RED ; 
	  switch(c)
	    {
	    case 'a':
	      color = A_COLOR ;
	      break ;
	    case 't':
	      color = T_COLOR ;
	      break ;
	    case 'g':
	      color = G_COLOR ;
	      break ;
	    case 'c':
	      color = C_COLOR ;
	      break ;
	    default:
	      color = ORANGE ;
	      break ;
	    }
	}
      buf[0] = c ;
      graphText(buf, *offset, y) ;
      graphBoxEnd() ;
      graphBoxDraw (box, BLACK, color) ;
    }
 
 *offset += 3.2 ; /* shift back (for errors sake) */
 return done ;
}    

/***********/

static int jackPotSearch (LOOK look, MAP map, LANE *lane,
				    float *offset, int min, int max)
{ 
  int i, maxA, n = 0, n1, n2, j, j1 = 0, oldsens = 0, sens ;
  Array a = lane->errArray ;
  A_ERR * ep = 0 ;
  char *cp, *cq ;
  int isQ = 0 ;  /* useless now, */
  
  lane->boxNewSpliceHint = 0 ;
  if (!a || (maxA = arrayMax(a)) < 5 || look->mode != CDNA)
    return 0 ;
  /* try to propose alternative editions */
  /* going down */
  
  /*   if (lane->upSequence)     return 0 ; */

  if (lane->getNewSpliceHint)
    oldsens = lane->lastSpliceHint < 0 ? -1 : 1 ;
  sens = oldsens ? oldsens : 1 ;
  
 for (;sens > -2 ; sens -= 2)
   {
     if (oldsens && sens != oldsens)
       break ;
     if (sens == 1)
       {
	 for (i = 0, ep = arrp (a, i, A_ERR) ; 
	       i < maxA  && ep->iLong < look->wmin  ; i++, ep++)
	   {} ;
	 n2 = i < maxA ? look->wmax - ep->iLong : 0 ; j = i ;
	 for (n1 = 1 ; i < maxA && ep->iLong < look->wmax  ; i++, n1++, ep++) ;
	 ep = j < maxA ? arrp (a, j, A_ERR) : 0 ;
       }
     else
       {
	 for (i = arrayMax (a) - 1, ep = arrp (a, i, A_ERR) ; 
	      i >= 0 && ep->iLong > look->wmax ; i--, ep--) ;
	 n2 = i >= 0 ? ep->iLong - look->wmin : 0 ; j = i ;
	 for (n1 = 0 ; i >= 0 && ep->iLong > look->wmin ; i--, n1++, ep--) ;
	 ep = j >= 0 ? arrp (a, j, A_ERR) : 0 ;
       }
     
     if (!ep || n1 < 4 || 3*n1 < n2) 
       continue ;
     
     switch (lane->getNewSpliceHint)
       {
       case 0: j1 = 20 ; break ; 
       case 1: j1 = oldsens * lane->lastSpliceHint ; break ;
       case 2: j1 = oldsens * lane->lastSpliceHint + 20 ; break ;
       }
     
     for (; j1 < 12000 ; j1++)
       {
	 j = ep->iLong  + j1 * sens ;
	 if (j < 0 || j >= arrayMax(look->dna))
	   break ;
	 
	 cp = arrp (lane->dna, ep->iShort, char) ;
	 cq = arrp (look->dna, j, char) ;
	 n = lane->jackPotMatchLength ? lane->jackPotMatchLength : 6 ; 
	 lane->jackPotMatchLength = n ;
	 while (j1 > (1 << (2 * n - 4))) n++ ;
	 isQ = 0 ;
	 while (n--) 
	   if (!lane->upSequence)
	     {
	       if ( !(*cp & *cq))
		 break ;
	       if (*cp != *cq)
		 { n++ ; isQ++ ; }
	       if (isQ > 1) break ;
	       if (sens == 1) { cp++ ; cq++ ; }
	       else { cp-- ; cq-- ; }
	     }
	   else
	     {
	       if ( !(*cp & complementBase(*cq)))
		 break ;
	       if (*cp !=  complementBase(*cq))
		 { n++ ; isQ++ ; }
	       if (isQ > 1) break ;
	       if (sens == 1) { cp-- ; cq++ ; }
	       else { cp++ ; cq-- ; }
	     }
	 if (n > -1)
	   continue  ;
	 /* now count the new set of errors in region */
	 cp = arrp (lane->dna, ep->iShort, char) ;
	 cq = arrp (look->dna, j, char) ;
	 if (sens == 1)
	   {
	     n = look->wmax - ep->iLong ; n2 = 0 ;
	     while (n--) 
	       {
		 if (!lane->upSequence)
		   { if (*cp++ != *cq++)  n2++ ;}
		 else
		   { if (*cp-- !=  complementBase(*cq++))  n2++ ;}
	       }
	   }
	 else
	   {
	     n = ep->iLong - look->wmin ; n2 = 0 ;
	     while (n--) 
	       {
		 if (!lane->upSequence)
		   { if (*cp-- != *cq--)  n2++ ;}
		 else
		   { if (*cp++ !=  complementBase(*cq--))  n2++ ;}
	       }
	   }
	 
	 if (n2 >  n1 - 3)
	   continue ;
	 
	 /* success */
	 if (!jackPotHint (look, map, lane, offset, min, max, ep->iShort, ep->iLong, j1, sens, FALSE))
	   return 0 ;
	 return 2 ; /* 2 => do not draw also local hints */
       }
   }
 if (oldsens)
   jackPotHint (look, map, lane, offset, min, max, ep->iShort, ep->iLong, j1, oldsens, TRUE) ;
 return 0 ;
}     

/***********/

static void traceDnaLaneErrors (LOOK look, MAP map, LANE *lane, float *offset, int min, int max, BOOL foundSplice)
{ int n = 0 ;
  float y ;
  int  i, nBase , box , color, sens ;
  char *base ; short *basePos , bp ;
  char c, buf[2] ;
  Array a = 0 ;
  A_ERR * ep ; HINT *hh ;
  
  CHECKSEQ ;
  if (lane->scf < 3)
    return ;
  CHECKPOS ;
  a = lane->errArray ;
  buf[1] = 0 ;
  nBase = arrayMax (lane->base) ;

  if (min < 0)
    min = 0 ;
  if (max > lane->maxPos)
    max = lane->maxPos ;

  i = arrayMax(a) ;
  if (!i)
    return ;

  base = arrp (lane->base, 0, char) ;
  basePos = arrp (lane->basePos, 0, short) ;

  while(*basePos < min  && nBase)
    { basePos++ ; base++ ; nBase-- ; n++ ;
    }

  lane->minBase = base - arrp (lane->base, 0, char) ;
  if (lane->upSequence)
    { sens = -1 ;
      ep = arrp(a, i - 1, A_ERR)  ;
    }
  else
    { sens = 1 ;
      ep = arrp(a, 0, A_ERR)  ;
    }
  while (i && ep->iShort < lane->minBase) { ep += sens ; i-- ;}
  if (i < 0)
    return ;

  *offset -= 1.2 ;

  color = WHITE ;
  i++ ; ep -= sens ;
  while (ep += sens , 
	 --i > 0 && 
	 ep->iShort >= 0 &&
	 ep->iShort < arrayMax(lane->basePos) &&
	 ep->iShort < lane->clipEnd &&
	 ep->iLong >= 0 &&
	 ep->iLong < arrayMax(look->dna) 
	 )	 
    { bp = arr (lane->basePos, ep->iShort, short) ;
      if (bp > max)
	break ;
      y = MAP2GRAPH(map, bp) ;
      c = dnaDecodeChar[(int)arr(look->dna, ep->iLong, char)] ;
      c = ace_lower(c) ;
      box = graphBoxStart () ;
      hh = arrayp (look->hints, box, HINT) ;
      hh->lane = lane ;
      hh->iShort = ep->iShort ;
      hh->iLong = ep->iLong ;
      hh->y = y ;
      if (foundSplice) hh->isSpliceHint = 1 ;
      switch(ep->type)
	{
	case TROU_DOUBLE:
	case TROU_TRIPLE:
	case TROU:
	  hh->key = ace_upper(c) ;
	  if(lane->upSequence)
	    hh->iShort++ ;
	  break ;
	case AMBIGUE: 
	case ERREUR: 
	  hh->key = c ;
	  break ;
	case INSERTION: /*_AVANT: */
	case INSERTION_DOUBLE: /*_AVANT: */
	case INSERTION_TRIPLE: /*_AVANT: */
	  hh->key = 'd' ; 
	  break ;
	case TYPE80:
	  break ;
	}
 
      buf[0] = ' ' ;
      switch(ep->type)
	{
	case TROU:
	  graphText(buf, *offset, y - 1) ;
	case AMBIGUE: 
	case ERREUR: 
	  color = RED ; 
	  switch(c)
	    {
	    case 'a':
	      color = A_COLOR ;
	      break ;
	    case 't':
	      color = T_COLOR ;
	      break ;
	    case 'g':
	      color = G_COLOR ;
	      break ;
	    case 'c':
	      color = C_COLOR ;
	      break ;
	    }
	  break ;
	case INSERTION: /*_AVANT: */
	  color = CYAN ;
	  c = '*' ;
	  break ;
	case INSERTION_DOUBLE: /*_AVANT: */
	case INSERTION_TRIPLE:
	  color = MAGENTA ;
	  c = '#' ;
	  break ;
	case TROU_DOUBLE:
	  graphText(buf, *offset, y - 2) ;
	  graphText(buf, *offset, y - 1) ;
	  color = MAGENTA ;
	  c = 'X' ;
	  break ;
	case TROU_TRIPLE:
	  graphText(buf, *offset, y - 3) ;
	  graphText(buf, *offset, y - 2) ;
	  graphText(buf, *offset, y - 1) ;
	  color = MAGENTA ;
	  c = 'X' ;
	  break ;
	case TYPE80:
	  break ;
	}
      buf[0] = c ;
      graphText(buf, *offset, y) ;
      graphBoxEnd() ;

      
      graphBoxDraw (box, BLACK, color) ;
    }
  *offset += 1.2 ;
}    

/***************************************************************************************/
/****************************** plots **************************************************/

static	BOOL traceLaneCrossing (Array lanes, LANE *lane)
{
  int nl = arrayMax(lanes) ;
  LANE *ll = arrp(lanes, 0, LANE) - 1 ;

  while (ll++, nl--)
    if (ll->key != lane->key && 
	ll->clone &&
	ll->clone == lane->clone &&
	((lane->upSequence && !ll->upSequence && lane->x2 < ll->x2) ||
	 (!lane->upSequence && ll->upSequence && lane->x2 > ll->x2)
	 ))
      return TRUE ;
  return FALSE ;
}

/***************************************************************************************/

static void laneHeight (LANE *lane, TRACE **ap)
{ int k, i ;
  TRACE *sp, hMin = 0, hMax = 0 ; 
  int min = 0, max = lane->maxPos ;

  if (lane->basePos && arrayMax(lane->basePos) > 2)
    min = arr (lane->basePos,1,short) ;
  hMin = hMax = 0 ;
  for (k=0 ; k<4 ; k++)
    { sp = ap[k] ;
      for (i = min, sp += min ; i < max ; i++, sp++)
	{ if (hMax < *sp && *sp  < 1000) hMax = *sp ;
	  if (hMin > *sp) hMin = *sp ;
	}
    }
  lane->laneMax = hMax > hMin ? hMax : -1 ; /* -1 will prevent recursion */
  lane->laneMin = hMin ;
}

/***************************************************************************************/
static Array
  histoTrace0 = 0,
  histoTrace1 = 0, 
  histoTrace2 = 0,
  histoTrace3 = 0,
  histoTrace4 = 0 ;

static void traceDoPlot (LOOK look, 
			 LANE *lane, int min, int max, float *offset,
			 MAP map, 
			 TRACE *ap[4],
			 int laneMin, int laneMax, 
			 int width, int seqmax)
{ int k, i, p, origin = look->origin ;
  float hScale, x0, x1, x2, y1, y2, dy;
  int box ;
  float oldW ;
  BOOL compl ;
  static Array points = 0 ;

  if (laneMax <= 0 || width < 5)
    return ;

  hScale = (width - 3) / (float) (laneMax - laneMin) ;
  x0 = *offset - laneMin * hScale ;
  if (min < 0)
    min = 0 ;
  if (max > seqmax)
    max = seqmax ;

  if (!lane->showTraceBeforeClipTop &&
      min < lane->xClipTop )
    min = lane->xClipTop ; 
  if (max > lane->xClipEnd )
    max = lane->xClipEnd ;
  if (max <min+1) max = min+1 ;
  oldW = graphLinewidth(.3) ;
  for (k=0 ; k<4 ; k++)
    { int
	z0, z1 , z2, z3, z4, /* various predictors */
	u1, u2, u3, u4 ;    /* previous values */
      TRACE
	*sp = ap[k] ;   /* actual trace */

      u1 = u2 = u3 = u4 = 0 ;
      box = graphBoxStart() ;
      
      
      dy = MAP2GRAPH(map,min + 1)  - MAP2GRAPH(map,min) ; 
      points = arrayReCreate (points, 2*(max - min) + 8 , float) ; p = 0 ;
      for (i = min, sp += min, x1 = x0 + *sp * hScale , y1 = y2 = MAP2GRAPH(map,min), 
	  array (points, p++, float) = x1, 
	  array (points, p++, float) = y1 ;
	   i < max ; i++, sp++)
	{ x2 = x0 + *sp * hScale ;
	  y2 += dy ;
	
	  z0 = 0 ; 
	  z1 = u1 ;
	  if (z1 < 0) z1 = 0 ;
	  if (z1 > 255) z1 = 255 ;
	  z2 = 2*u1 - u2 ;
	  if (z2 < 0) z2 = 0 ;
	  if (z2 > 255) z2 = 255 ;
	  z3 = 3*u1 - 3*u2 + u3 ;
	  if (z3 > 255) z3 = 255 ;
	  if (z3 < 0) z3 = 0 ;
	  z4 = 4*u1 - 6*u2 + 4*u3 - u4;
	  if (z4 > 255) z4 = 255 ;
	  if (z4 < 0) z4 = 0 ;
	  if (i < min + 1) z1 = z0 ;
	  if (i < min + 2) z2 = z1 ;
	  if (i < min + 3) z3 = z2 ;
	  if (i < min + 4) z4 = z3 ;
	  if (histoTrace0)
	    array(histoTrace0, (300 + (*sp - z0 )) , int)++ ;
	  if (histoTrace1)
	    array(histoTrace1, (300 + (*sp - z1 )) , int)++ ;
	  if (histoTrace2)
	    array(histoTrace2, (300 + (*sp - z2 )) , int)++ ;
	  if (histoTrace3)
	    array(histoTrace3, (300 + (*sp - z3 )) , int)++ ;
	  if (histoTrace4)
	    array(histoTrace4, (300 + (*sp - z4 )) , int)++ ;
	  u4 = u3 ; u3 = u2 ; u2 = u1 ; u1 = *sp ;

/* 	here is the place to make any transform you like
     on the display 
     example:   x2 = x0 + z3 * hScale ;
*/
	  array (points, p++, float) = x2 ;
	  array (points, p++, float) = y2 ;
	  x1 = x2 ; y1 = y2 ;
	}
      graphLineSegs (points) ;
      graphBoxEnd() ;
      graphBoxDraw(box, 
		   BASE_COLOR[lane->upSequence ? 
			      3 - k : k], TRANSPARENT) ;
    }

  graphLinewidth(oldW) ;
  lane->boxLessClip = lane->boxMoreClip = 0 ;
  compl = lane->x1 > lane->x2 ;
  if ((compl && y2 > look->topMargin + 4 ) ||
      (!compl && y2 < look->graphHeight - 4))
    { 
      float 
	x = x0 + 1, y,
	dy = compl ? -1 : 1 ;
      
      if (lane->nextTop)
	{
	  y = compl ? y2 - 2 : y2 + 2 ;
	  graphText (messprintf ("%d", lane->x2 - origin + 1), x, !compl ? y2 + 1.5 : y2 - 1.5) ;
	  graphText (messprintf ("%d", lane->nextTop - origin + (compl ? 1 : - 1)), x, !compl ? y2 + 2.5 : y2 - 2.5) ;
	  y = compl ? y2 - 4 : y2 + 3 ; 
	}
      else
	{
	  y = compl ? y2 - 2 : y2 + 2 ; 
	  
	  box = graphBoxStart () ;
	  
	  lane->boxLessClip = box ;
	  graphLine (x, y, x + 2 , y) ;
	  graphLine (x, y, x + 1 , y - dy ) ;
	  graphLine (x + 1 , y - dy, x + 2 , y) ;
	  
	  graphBoxEnd () ;
	  if (lane->handClipTop && lane->clipTop == lane->handClipTop)
	    graphBoxDraw(box, BLACK, ORANGE) ;
	}

      y += compl ? - .7 : .7 ;
      box = graphBoxStart () ;
      if (lane->clipEnd < arrayMax(lane->dna))
	{ 
	  lane->boxMoreClip = box ;
	  graphLine (x, y , x + 2 , y ) ;
	  graphLine (x, y, x + 1 , y + dy) ;
	  graphLine (x + 1 , y + dy , x + 2 , y) ;
	}
      else
	{ 
	  lane->boxMoreClip = 0 ;
	  graphRectangle (x, y , x + 2 , y + dy) ;
	}
      graphBoxEnd () ;
      if (lane->handClipEnd && lane->clipEnd == lane->handClipEnd)
	    graphBoxDraw(box, BLACK, ORANGE) ;
      if (traceLaneCrossing (look->lanes, lane)) 
	graphBoxDraw(box, BLACK, ORANGE) ;
    }
  lane->boxLessClipTop = lane->boxMoreClipTop = 0 ;
  y1 = MAP2GRAPH(map,min) ; 
  if (look->mode == CDNA &&
      lane->clipTop > 0   /* lane->vectorTop */ &&
      (
       (!compl && y1 > look->topMargin + 4 ) ||
       (compl && y1 < look->graphHeight - 4))
      )
    { 
      float 
	x = x0 + 1, y = compl ? y1 + 4 : y1 - 4 , 
	dy = compl ? -1 : 1 ;
      
      if (lane->previousEnd)
	{
	  graphText (messprintf ("%d", lane->x1 - origin + (lane->upSequence ? 2 : 0)), x, compl ? y1 + 1.5 : y1 - 1.5) ;
	  graphText (messprintf ("%d", lane->previousEnd - origin + (lane->upSequence ? -1 : 1)), x, compl ? y1 + 2.5 : y1 - 2.5) ;
	}
      box = graphBoxStart () ;
      
      lane->boxLessClipTop = box ;
      graphLine (x, y, x + 2 , y) ;
      graphLine (x, y, x + 1 , y - dy ) ;
      graphLine (x + 1 , y - dy, x + 2 , y) ;
      
      graphBoxEnd () ;
    }
  
  laneShowExtrema (look, lane, offset, min, max) ;
  *offset += width ;
}

/***********/

static void tracePlotLane(LOOK look, LANE *lane, float *offset, int min, int max)
{
  TRACE *bp[4] ;
  
  lane->boxClose = graphBoxStart () ;
  graphCircle ( *offset + 2.345, look->topMargin - .9, 1.1) ;
  graphText(messprintf("e"),  *offset + 2, look->topMargin - 1.45) ;
  graphBoxEnd() ;
  if (lane->favorDelete == 2)
    graphText(messprintf("*"),  *offset + 4, look->topMargin - 1.45) ;
  else if (lane->favorDelete == 3)
    graphText(messprintf("#"),  *offset + 4, look->topMargin - 1.45) ;

  bp[0] = lane->seq->traceA ;
  bp[1] = lane->seq->traceG ;
  bp[2] = lane->seq->traceC ;
  bp[3] = lane->seq->traceT ;
  
  if (!lane->laneMax)
    laneHeight (lane, bp) ;

  traceDoPlot(look, lane, min, max, 
	      offset, 
	      look->map,
	      bp, 
	      lane->laneMin, lane->laneMax,
	      look->traceWidth, lane->maxPos
	      ) ;
}

/***********/

static BOOL isLaneShown (LOOK look, LANE *lane)
{
  int
    x1, x2,
    u1 = look->wmin + 1 ,
    u2 = look->wmax - 1 ;
  
  if (look->hideUp == HIDE_POOR &&
      lane->isPoor)
    return FALSE ;

  if (lane->upSequence)
    { x1 = lane->x2 ; 
      x2 = lane->x1 ;
    }
  else
    { x1 = lane->x1 ; 
      x2 = lane->x2 ;
    }

  if (x1 > u2 || x2 < u1)
    return FALSE ;
  
  CHECKLANE ;
  return lane->scf > 1;
}

/***********/

static BOOL laneNoError (LOOK look, LANE *lane, BOOL isClipTop)
{
  int j ;
  A_ERR *ep ;

  if (!lane->errArray)  /* happens after a fix */
    {
      CHECKLANE ;
    }

  if (!arrayMax(lane->errArray)) 
    return TRUE ;
  if (isClipTop)
    {
      if (!lane->upSequence)
	{
	  ep = arrp (lane->errArray, 0, A_ERR) - 1 ;
	  j = arrayMax (lane->errArray) ;
	  while (ep++, j--)
	    if (ep->iLong > lane->x1)
	      {  /* first error further than splice site + 9 */
		if (ep->iLong > lane->x1 + 7)
		  return TRUE ;
		break ;
	      }
	}
      else
	{
	  ep = arrp (lane->errArray, arrayMax(lane->errArray) -1 , A_ERR) + 1 ;
	  j = arrayMax (lane->errArray) ;
	  while (ep--, j--)
	    if (ep->iLong < lane->x1 - 8)
	      {  /* first error further than splice site + 9 */
		
		break ;
	      } 
	  if (j == arrayMax (lane->errArray) - 1)   return TRUE ;
	}
    }
  else
    {
      if (lane->upSequence)
	{
	  ep = arrp (lane->errArray, 0, A_ERR) - 1 ;
	  j = arrayMax (lane->errArray) ;
	  while (ep++, j--)
	    if (ep->iLong > lane->x2 + 1)
	      {  /* first error further than splice site + 9 */
		if (ep->iLong > lane->x2 + 8)
		  return TRUE ;
		break ;
	      }
	}
      else
	{
	  ep = arrp (lane->errArray, arrayMax(lane->errArray) - 1 , A_ERR) + 1 ;
	  j = arrayMax (lane->errArray) ;
	  while (ep--, j--)
	    if (ep->iLong < lane->x2 - 8)
	      {  /* first error further than splice site + 9 */
		break ;
	      }
	  if (j == arrayMax (lane->errArray) - 1)  return TRUE ;
	}
    } 
  return FALSE ;
}

/***********/

static void selectBestTraces (LOOK look)
{ MAP map = look->map ;
  LANE
    *bestUp, *bestDown,
    *lane, *lane1 = 0 ;
  int 
    i, e, eUp, eDown, 
    wmin = look->wmin, wmax = look->wmax,
    wc1 = wmin + (wmax - wmin) / 6, 
    wc2 = wmax - (wmax - wmin) / 6, 
    nUp = 0, nDown = 0,
    nPerfectUp = 0, nPerfectDown = 0, nClippedNoError = 0,
    j,  min, max ;
  KEY key ;
  float mag ;
  A_ERR *ep ;

     /* check and count */
  i = arrayMax (look->laneShown) ;
  if (!i) return ;
  while (i--)
    { lane = arr (look->laneShown ,i, LANE*) ;
      if (lane->isAligned) continue ;
      lane->hide = 0 ;  /* abuse lane->hide */
      if (look->hideUp == HIDE_POOR &&
	  lane->isPoor)
	continue ;
      key = lane->clone ? lane->clone : lane->key ;
      if (look->select &&
	  keySetExists(look->selectedKs) &&
	  !keySetFind (look->selectedKs, key, 0))
	continue ;
      if (!fineTune (look, &lane1, lane))
	continue ;
      CHECKLANE ;
      if (lane->scf < 2) /* data is lacking */
	continue ;
      
      lane->hide = 1 ;              /* 1: data exists */
      if (lane->upSequence)
	nUp++ ;
      else
	nDown++ ;
      
      if (arrayMax(lane->errArray))
	{ ep = arrp (lane->errArray, 0, A_ERR) - 1 ;
	  j = arrayMax (lane->errArray) ;
	  while (ep++, j--)
	    { if (ep->iLong < wc1)
		continue ;
	      if (ep->iLong < wc2)
		lane->hide |=                     /* 2: ambigue */
		  (ep->type == AMBIGUE ?          /* 4: ambigue editee */
		   ((BC_HAND & arr (lane->base, ep->iShort, char)) ?
		    4 : 2)  : 8) ;                /* 8: error   */
 
	      if (ep->iLong >= wc2)
		break ;
	    }

	}
      if (lane->x1 > look->wmin && lane->x1 < look->wmax)
	{
	  lane->hide |= 16 ;                        /* 16: clip top */
	  if (look->hide == SHOW_SPLICER)
	    {   /* see if there is no error close to the splice site */
	      if (laneNoError (look, lane, TRUE))
		lane->hide |= 128 ; ;
	    }
	  if (lane->hide & 128) /* verify that the other side of the intron is exact */
	    {
	      lane->hide &= ~128 ; /* reset false untill we verify other side */
	      for (j = 0, lane1 = arrp (look->lanes,j, LANE) ; j < arrayMax(look->lanes) ; j++, lane1++)
		if (lane1->key == lane->key && lane1->x2 == lane->previousEnd)
		  { /* check no error */
		    if (laneNoError (look, lane1, FALSE))
		      lane->hide |= 128 ; ;
		  }
	    }
	}
      
      if (lane->x2 > look->wmin && lane->x2 < look->wmax)
	{ 
	  lane->hide |= 32 ;                        /* 32: clip end */ 
	  if (look->hide == SHOW_SPLICER)
	    {   /* see if there is no error close to the splice site */
	       if (laneNoError (look, lane, FALSE))
		lane->hide |= 128 ; ;
	    }
	  if (lane->hide & 128) /* verify that the other side of the intron is exact */
	    {
	      lane->hide &= ~128 ; /* reset false untill we verify other side */
	      for (j = 0, lane1 = arrp (look->lanes,j, LANE) ; j < arrayMax(look->lanes) ; j++, lane1++)
		if (lane1->key == lane->key && lane1->x1 == lane->nextTop)
		  { /* check no error */
		    if (laneNoError (look, lane1, TRUE))
		      lane->hide |= 128 ; ;
		  }
	    }
	}
      
      if (lane->hide == 1)
	{ 
	  if (lane->upSequence)
	    nPerfectUp++ ;
	  else
	    nPerfectDown++ ;
	}
      else if (!(lane->hide & 10))
	nClippedNoError++ ;
      
    }
  switch (look->hide)
    { 
    case SHOW_ALL:
      goto ok ;
    }
     /* look for best among perfect */
  bestUp = bestDown = 0 ; eUp = eDown = 0 ;

  i = arrayMax (look->laneShown) ;
  while (i--)
    { lane = arr (look->laneShown ,i, LANE*) ;
      if (lane->isAligned) continue ;
      if (!(lane->hide & 1)) continue ;
      if (lane->upSequence)
	{
	  switch (nPerfectUp)
	    {
	    case 0:
	      if (nClippedNoError && (lane->hide & 10))
		continue ;
	      break ; /* measure the energy of all bads */
	    case 1:  /* single candidate */
	      if (lane->hide == 1)
		bestUp = lane ;
	      continue ;
	    default:  /* only measure the best ones */
	      if (lane->hide != 1)
		continue ;
	    }
	}
      else
	{
	  switch (nPerfectDown)
	    {
	    case 0:
	      if (nClippedNoError && (lane->hide & 10))
		continue ;
	      break ; /* measure the energy of all bads */
	    case 1:  /* single candidate */
	       if (lane->hide == 1)
		 bestDown = lane ;
	      continue ;
	    default:  /* only measure the best ones */
	      if (lane->hide != 1)
	      continue ;
	    }
	}

      map->centre -= lane->dy ;
      mag = map->mag ;
      map->mag *= lane->ddy ;
      
      min = GRAPH2MAP(map, look->topMargin) ;
      max = GRAPH2MAP(map, look->graphHeight) ;
      if (lane->upSequence)
	{ j = min ; min = max ; max = j ; }
      e = seqEnergyOfDerivee(lane->seq, min, max) ;
      map->centre += lane->dy ;
      map->mag = mag ;
      
      if (lane->upSequence && e > eUp)
	{ eUp = e ;
	  bestUp = lane ;
	}
      if (!lane->upSequence && e > eDown)
	{ eDown = e ;
	  bestDown = lane ;
	}
    }
	
  if (bestUp) bestUp->hide |= 64 ;                 /* 64: best */
  if (bestDown) bestDown->hide |= 64 ;

  /* restore */
 ok:
  i = arrayMax (look->laneShown) ;
  if (look->mode == CDNA && look->hide == SHOW_BEST_ERR)
    {
      int x1, x2 ;
      KEYSET ks = keySetCreate () ;
      while (i--)
	{
	  lane = arr (look->laneShown, i, LANE*) ;
	  /*
	    printf("%s hide = %d x1 = %d x2 = %d\n",
		 name(lane->key), lane->hide, lane->x1 - look->origin , lane->x2 - look->origin) ;
	  */

	  switch (look->hideUp)
	    {
	    case HIDE_UP:
	      if (lane->upSequence)
		{ lane->hide = TRUE ; continue ; }
	      break ;
	    case HIDE_DOWN:
	      if (!lane->upSequence)
	    	{ lane->hide = TRUE ; continue ; }
	      break ;
	    case HIDE_POOR:
	      if (lane->isPoor)
	    	{ lane->hide = TRUE ; continue ; }
	      break ;
	    case HIDE_NON_FULL:
	      if (!lane->clone || !keyFindTag (lane->clone, str2tag("Complete_CDS_of")))
		{ lane->hide = TRUE ; continue ; }
	      break ;
	    }
	  if (!(lane->hide & 1))                       /* can't fineTune */
	    { lane->hide = TRUE ; continue ; }
	  if (lane->hide & 10)                          /* errors */
	    { lane->hide = FALSE ; continue ; }
	  
	  if (
	      (
	       (lane->hide & 16) &&                  /* clip top 5p */
	       !lane->upSequence &&
	       lane->clipTop > 0 &&
	       (lane->clipTop != lane->vectorTop )
	       ) ||
	      (
	       (lane->hide & 16) &&                  /* clip top 3p  */
	       lane->upSequence &&
	       lane->clipTop > 0 &&
	       (lane->clipTop != lane->vectorTop )
	       ) ||
	       (
		(lane->hide & 32) &&                  /* clip end 5p */
		!lane->upSequence &&
		(lane->nextTop || !traceLaneCrossing (look->lanes, lane))
		) ||
	      (
	       (lane->hide & 32) &&                  /* clip end 3p  */
	       lane->upSequence &&
	       (lane->nextTop || !traceLaneCrossing (look->lanes, lane))
	       )
	      )	      
	    {
	      x1 = lane->x1 ; x2 = lane->x2 ;
	      lane->hide = TRUE ;

	      if (x1 > 0 && x1 > look->wmin && x1 < look->wmax)
		{ 
		  if (!keySet( ks, x1))
		    { 
		      keySet (ks, x1) = 1 ;
		      lane->hide = FALSE ;
		    }
		  x1 = lane->previousEnd ;
		  if (x1>= 0 && !keySet( ks, x1))
		    { 
		      keySet (ks, x1) = 1 ;
		      lane->hide = FALSE ;
		    }
		}
	      if (x2 > 0 && x2 > look->wmin && x2 < look->wmax)
		{ 
		  if (!keySet( ks, x2))
		    { 
		      keySet (ks, x2) = 1 ;
		      lane->hide = FALSE ;
		    }
		  x2 = lane->nextTop ;
		  if (x2 >= 0 && !keySet( ks, x2))
		    { 
		      keySet (ks, x2) = 1 ;
		      lane->hide = FALSE ;
		    }
		}
		
	      continue ;
	    }
	  else
	    {
	      if (lane->hide & 64)              /* best */
		{ lane->hide = FALSE ; continue ; }
	      
	      lane->hide = TRUE ;               /* default */
	    }
	}
      keySetDestroy (ks) ;
      return ;
    }
    
  i = arrayMax (look->laneShown) ;
  while (i--)
    { 
      lane = arr (look->laneShown, i, LANE*) ;
      if (!(lane->hide & 1))                       /* can't fineTune */
	{ lane->hide = TRUE ; continue ; }

      switch (look->hideUp)
	{
	case HIDE_UP:
	  if (lane->upSequence)
	    { lane->hide = TRUE ; continue ; }
	  break ;
	case HIDE_DOWN:
	  if (!lane->upSequence)
	    	{ lane->hide = TRUE ; continue ; }
	  break ;
	case HIDE_POOR:
	  if (lane->isPoor)
	    	{ lane->hide = TRUE ; continue ; }
	  break ;
	case HIDE_NON_FULL:
	  if (!lane->clone || !keyFindTag (lane->clone, str2tag("Complete_CDS_of")))
	    { lane->hide = TRUE ; continue ; }
	  break ;
	}

      if ((look->next == NEXT_CLIP || look->next == NEXT_HOLE) &&
	  lane == lastUnclippedLane)
	{ lane->hide = FALSE ; }
      else
	switch (look->hide)
	{
	case SHOW_ALL:
	  lane->hide = FALSE ;
	  break ;
	case SHOW_BEST_ERR:
	  if (lane->hide & (64 + 8 + 4 + 2))
	    lane->hide = FALSE ;
	  else
	    lane->hide = TRUE ;
	  break ;
	case HIDE_NON_BEST:
	  if (lane->hide & 64)
	    lane->hide = FALSE ;
	  else
	    lane->hide = TRUE ;
	  break ;
	case SHOW_CLIP:
	  if (lane->hide & (64 + 32 + 16))
	    lane->hide = FALSE ;
	  else
	    lane->hide = TRUE ;
	  break ;
	case SHOW_SPLICER: /* show lane contributing to splicing */
	  if (lane->hide & 128)
	    lane->hide = FALSE ;
	  else
	    lane->hide = TRUE ;
	  break ;
	}

    }
}

/****************/

static void tracePlotAllLanes (LOOK look, float *offset)
{ LANE *lane, *lane1 = 0 ;
  MAP map = look->map ;
  int i, j1, j2, min, max ;
  /*   BOOL old = TRUE ; */
  float mag, y, x0 = *offset ;

  j1 = j2 = 0 ;
  i = arrayMax (look->laneShown) ;
  if (!i)    return ;
  while (i--)
    { lane = arr (look->laneShown, i, LANE*) ;
      CHECKSEQ ; CHECKPOS ;
      if (!lane->hide && !lane->isAligned && lane->scf >= 3) j1++ ;
      if (!lane->hide && !lane->isAligned && lane->scf >= 4) j2++ ;
    }
  if (!j1)
    return ;
  if (!j2) j2 = 1 ;
  look->traceWidth = (look->graphWidth - *offset - 5 - 3.8 * j1) / j2 ;
  if (look->traceWidth < 12)  /* was 8 */
    look->traceWidth = 12 ;    /* was 6 */
  if (look->traceWidth > 18)
    look->traceWidth = 18 ;

  for (i = 0 ; i < arrayMax (look->laneShown) ; i++)
    { 
      lane = arr (look->laneShown, i, LANE*) ;
      if (lane->hide || lane->isAligned)
	continue ; 
      if (*offset > look->graphWidth - 3)
	{
	  if (lane->greenBox)
	    graphBoxDraw (lane->greenBox, WHITE, TRANSPARENT) ;
	  continue ;
	}
      if (!showSelect &&  /* did fineTune */
	  !fineTune (look, &lane1, lane) )
	continue ;

      map->centre -= lane->dy ;
      mag = map->mag ;
      map->mag *= lane->ddy ;

      
      min = GRAPH2MAP(map, look->topMargin) ;
      max = GRAPH2MAP(map, look->graphHeight) ;

      if (min > max)	  /* mag < 0, reverse to simplify maths */
	{ int
	    tmp = max ; max = min ; min = tmp ;
	}

#ifdef JUNK
      if (old != lane->upSequence)
	{ old = lane->upSequence ;
	  *offset += .8 ;
	  graphLine (*offset, min, *offset, max) ;
	  *offset += .4 ;
	  graphLine (*offset, min, *offset, max) ;    
	  *offset += .8 ;
	}

      fin = look->showClip ?  lane->xClipExtend : lane->xClipEnd ;
#endif

      lane->laneBaseCallBox = graphBoxStart () ;
      lane->laneBaseCallBoxOffSet = *offset ;
      lane->laneBaseCallBoxMin = min ;
      lane->laneBaseCallBoxMax = max ;
      CHECKSEQ ; CHECKPOS ;
      traceDnaLane (look, look->map, lane, offset, min, max) ;
      lane->isShown = TRUE ; /* enables box pick */
      graphBoxEnd () ;
      if (!lane->hide)
	{	  
	  if (lane->scf >= 4)
	    tracePlotLane (look,lane, offset, min, max) ;
	  else
	    { 
	      char *cp = name(lane->key) ;
	      char buf[2] ;
	      int y = 12 ;
	      
	      buf[1] = 0 ;
	      while (*cp && y < look->graphHeight - 4)
		{
		  buf[0] = *cp++ ;
		  graphText (buf, *offset + .2, y++) ;
		}
	      y += 6 ;
	      cp = "no trace" ;
	      while (*cp && y < look->graphHeight - 4)
		{
		  buf[0] = *cp++ ;
		  graphText (buf, *offset + .2, y++) ;
		}
	      
	      graphRectangle (*offset, look->topMargin + 2, 
			      *offset + 1.4, look->graphHeight - 2) ;
	      *offset += 3.0 ;
	    }
	}
      map->centre += lane->dy ;
      map->mag = mag ;
    }
  y = MAP2GRAPH (look->map, look->map->centre) ;
  graphLine (x0, y, *offset, y) ;
}

/************************************/

void traceScale (LOOK look, float *offset)
{
  float y, unit, subunit ;
  int
    x, origin = look->origin, 
    iUnit, iSubunit,
    max = look->wmax + 1 - origin , /* no zero */
    min = look->wmin + 1 - origin ;
 
  unit = subunit = 1.0 ;
  mapFindScaleUnit (look->map, &unit, &subunit) ;
  iUnit = unit + 0.5 ;
  if (iUnit < 1) iUnit = 1 ;
  iSubunit = subunit + 0.5 ;
  if (iSubunit < 1) iSubunit = 1 ;

  x = iUnit * (min/iUnit) ;
  while (x < min) x += iUnit ;
  for ( ; x < max ; x += iUnit)
    { y = MAP2GRAPH(look->map, x + origin - 1) + .5; /* NO ZERO */
      graphLine (*offset+1.0, y, *offset+1.8, y) ;
      graphText (messprintf ("%d", x),
			     *offset+1.95, y - .5) ;
    }
      
  x =  iSubunit * (min/iSubunit) ;
  while (x < min) x += subunit ;
  for ( ; x < max ; x += subunit)
    { y = MAP2GRAPH(look->map, x + origin - 1) + .5;
      graphLine (*offset+1.0, y, *offset+1.5, y) ;
    }
    
  graphLine (*offset+1, look->topMargin+2, *offset+1, look->graphHeight-0.5) ;
  *offset += 8 ;
}

/*****************************************************************/
/****************** left most, gene center    ********************/
/*****************************************************************/

void traceGene (LOOK look, float *offset)
{
   OBJ obj = 0 ;
   float  x, y1, y2 ;
  int i, x1, x2, box ; 
  int min = ACEDB_MAXINT, max = 0 ;
  KEY tt ;
  MAP map = look->map ;
  float mag = map->mag, centre = map->centre;
  BSunit *u ;
  Array geneSplicing = 0 ;

  if (look->mode != CDNA || !look->gene)
    return ;
  
  obj = bsCreate (look->gene) ;
  if (!obj)
    return ;
  
  geneSplicing = arrayCreate (36, BSunit) ;  
  bsGetArray (obj, _Splicing, geneSplicing, 4) ;
  i = arrayMax(geneSplicing);
  for (i = 0 ; i < arrayMax(geneSplicing) ; i += 4)
    { 
      u = arrp (geneSplicing, i, BSunit) ;
      u[0].i  += look->origin ;
      u[1].i  += look->origin ;
      if (u[0].i < min) min = u[0].i ;
      if (u[1].i > max) max = u[1].i ;
    }
  if (max < min) { min = 0 ; max = 1000 ; }
  look->geneMin = min ;
  look->geneMax = max ;
  
 if (!arrayMax(geneSplicing)) 
   goto done ;

  genex = x = *offset + .2 ;
  
  map->centre = (look->geneMax + look->geneMin) / 2 ;
  map->mag = .8 * ((float)look->graphHeight - look->topMargin) / (look->geneMax - look->geneMin) ;

  y1 = MAP2GRAPH (map, look->edMin) ;
  y2 = MAP2GRAPH (map, look->edMax) ;
  box = graphBoxStart () ;
  graphRectangle (x -1.0, y1, x + 1.4, y2) ;
  graphBoxEnd () ;
  graphBoxDraw (box, GREEN, GREEN) ;

  y1 = MAP2GRAPH (map, look->wmin) ;
  y2 = MAP2GRAPH (map, look->wmax) ;
  box = graphBoxStart () ;
  graphRectangle (x +.2, y1, x + 3, y2) ;
  graphBoxEnd () ;
  graphBoxDraw (box, GREEN, GREEN) ;

  look->geneBox = graphBoxStart () ;
  graphColor (MAGENTA) ;  
  for (i = 0 ; i < arrayMax(geneSplicing) ; i += 4)
    {
      u = arrp (geneSplicing, i, BSunit) ;
      x1 = u[0].i ;
      x2 = u[1].i ;
      tt = u[2].k ;

      y1 = MAP2GRAPH (map, x1) ;
      y2 = MAP2GRAPH (map, x2) ;
      if (tt == _Intron)
	{
	  float mid = (y1 + y2) / 2 ;
	  if (!u[3].s || strcmp("gt_ag", u[3].s))
	    graphColor (BLUE) ; 
	  graphLine (x +.2, y2, x + 1.6 , mid) ;
	  graphLine (x +1.6, mid, x + .2, y1 +.1) ;
	}
      else if (tt == _Alternative_intron)
	{
	  float mid = (y1 + y2) / 2 ;
	  if (!u[3].s || strcmp("gt_ag", u[3].s))
	    graphColor (BLUE) ; 
	  graphLine (x +.5, y2, x + 1.9 , mid) ;
	  graphLine (x +1.9, mid, x + .5, y1 +.1) ;
	}
      else if (tt == _Exon || tt == _First_exon || tt == _Last_exon)
	{ 
	  box = graphBoxStart () ;
	  graphRectangle (x +.2, y1, x + 1.4 , y2) ; 
	  graphBoxEnd () ;
	  graphBoxDraw (box,MAGENTA, PALEMAGENTA) ;
	  graphBoxSetPick (box, FALSE) ;
	}
      else if (tt == _Alternative_exon || tt == _Alternative_partial_exon)
	{ float old = graphLinewidth (.2) ;
	  box = graphBoxStart () ;
	  graphRectangle (x -.2, y1, x + 1.9 , y2) ; 
	  graphBoxEnd () ;
	  graphBoxDraw (box,ORANGE, TRANSPARENT) ;
	  graphBoxSetPick (box, FALSE) ;
	  graphLinewidth (old) ;
	}
      else if (tt == _Partial_exon) 
	{ 
	  graphRectangle (x +.2, y1, x + 1.4 , y2) ;
	}
      else if (tt == _Gap)
	{ 
	  graphColor (BLACK) ;
	  graphLine (x +.8, y1, x + .8 , y2) ; 
	}
      else
	{
	  graphColor (RED) ;
	  graphRectangle (x +.1, y1, x + .4 , y2) ;
	}
      graphColor (MAGENTA) ;
    }
  graphColor (BLACK) ;
  graphBoxEnd () ;
  graphBoxDraw (look->geneBox,MAGENTA, TRANSPARENT) ;
 
  map->mag = mag ;
  map->centre = centre ;
  *offset += 2 ;

 done:
  bsDestroy (obj) ;
  return ;
}

/*****************************************************************/
/****************** left half/sequence editor ********************/
/*****************************************************************/
#ifdef JUNK
static void countVisibleTraces (float x)
{
   int i = arrayMax(look->lanes) ;
   LANE *lane = i ? arrp(look->lanes,0,LANE) - 1 : 0 ;

   x += 1.6   /* edconcensus */
     + 1 ;   /* trace eds */

   while(lane++, i--) 
     if (isLaneShown (look, lane))
}
#endif
/*****************************************************************/

void traceEdCoords (LOOK look, float *offset)
{
  float  x, y ;
  int n, i, di, origin = look->origin ;
  char *cp = 0 ;
  int start = look->edMin - origin, stop = look->edMax - origin ;
  MAP map = look->map ;
  float mag = map->mag ;
  
  map->mag *= look->edMag ;

  graphTextHeight (0.75) ;
  n = stop - start  ; n /= 4 ;
  di = 1 ;
 lao:
  if (2*di > n) goto done ;
  di *= 2 ;
  if (5*di > 2*n) goto done ;
  di *= 5 ; di /= 2 ;
  if (2*di > n) goto done ;
  di *= 2 ;
  goto lao ;
 done:
  for (i = ((start + di - 1)/di) * di ; i < stop ; i += di)
    { x = i ;
      y = MAP2GRAPH(map,x + origin) - 1 ;
      graphText (cp = messprintf("%d", i),  *offset+0.9, y) ;
    }
  graphTextHeight (0) ;

  map->mag = mag ;
  *offset += .5 + ( cp ? strlen(cp) : 0 );
  /* countVisibleTraces (*offset) ; */
}

/***************************************************************************************/
  /* start, stop are the first and last consensus base on screen */
static void traceEdConsensus (LOOK look, float *offset)
{
  char *cp, buf[2] ;
  float x = *offset, y ;
  int start = look->edMin, stop = look->edMax ;
  MAP map = look->map ;
  int i , j ;
  float mag = map->mag ;

  map->mag *= look->edMag ;

  buf[1] = 0 ;
  if (start >= arrayMax(look->dna)) 
    start = arrayMax(look->dna) - 1 ;
  if (start < 0) 
    start = 0 ;
  if (stop > arrayMax(look->dna)) 
    stop = arrayMax(look->dna) ;
  i = stop - start ;
  if (i > 0)
    {
      cp = arrp (look->dna, start, char) ;
      j = start ;
      while (i--)
	{
	  buf[0] = dnaDecodeChar[(int)*cp++] ;
	  y = MAP2GRAPH(map, j++) ; /* no zero */
	  graphText(buf, x, y) ;
	}
    }
  map->mag = mag ;
  *offset += 1.6 ; /* see traceEdCoords */
}

/**********************************************************************/

  /* start, stop are the first and last consensus base on screen */
static void traceEdLaneBases (LOOK look, LANE *lane, float *offset,
			       int start, int stop, int type, BOOL isAligned)
{
  char c, buf[2] ;
  int 
    color, i, nn, u1, u2, u3, nerr = 0, 
    v1 = lane->vectorTop, v2 = lane->vectorEnd, w1, w2,
    x1 = lane->x1, x2 = lane->x2 , x3 = lane->x3 ;
  BOOL upSequence = lane->upSequence ;
  float x = *offset, y = look->topMargin, ddx, ddy ;
  Array a ;
  MAP map = look->map ;
  A_ERR * ep ;

  if (look->hideDots & 0x1)   /* 1:auto or 3:hand set to hide */
    return ;
  a = lane->errArray ;
/*   baseCallGet (lane) ; */
  if (lane->scf != type || lane->isAligned != isAligned)
    return ;
  if (lane->scf == 2) *offset += .5 ;
  if (!look->showClip)
    x3 = x2 ;
  if (upSequence)
    { u1 = x3 ; u2 = x2 ; u3 = x1 ; }
  else
    { u1 = x1 ; u2 = x2 ; u3 = x3 ; }
 
  if (lane->hasTag != 1)
    traceDrawDnaTag (look, look->map, lane, offset) ;
  x = *offset ;
  i = GRAPH2MAP (look->map, look->topMargin) ;
  if (start < i) start = i ;
  i = GRAPH2MAP (look->map, look->graphHeight) ;
  if (stop > i) stop = i ;
  if (u3 < start || u1 >= stop)
    return ;

  lane->edLaneBaseBox = graphBoxStart () ;
  lane->edLaneBaseBoxOffSet = *offset ;
  if (u1 < start) u1 = start ;
  if (u2 < start) u2 = start ;
  if (u2 > stop) u2 = stop ;
  if (u3 > stop) u3 = stop ;

/* show clipping */
   if (!lane->hide)
     {
       int ngreen = 0, box = graphBoxStart() ; 
       graphRectangle (*offset - .2, MAP2GRAPH(map, start - .8),
		       *offset + 1.0, MAP2GRAPH(map, stop)) ;
       graphBoxEnd () ;
       if (ngreen++ < 1000) /* look->maxVisibleTraces) */
	 graphBoxDraw (box, PALEGREEN, PALEGREEN) ;
       graphBoxSetPick (box, FALSE) ;
       lane->greenBox = box ;
     }

  if (!upSequence && start < u1) 
    { int box = graphBoxStart() ;
      graphRectangle (*offset - .1, MAP2GRAPH(map, start - .8),
		      *offset + .9, MAP2GRAPH(map, u1 - .15)) ;
      graphBoxEnd() ;
      color = YELLOW ;
      if (look->mode == CDNA && lane->previousEnd)
	color = PALEMAGENTA ;
      graphBoxDraw (box, BLACK, color) ;
      graphBoxSetPick (box, FALSE) ;
    }
  if (upSequence && stop > x1) 
    { int box = graphBoxStart() ;
      graphRectangle (*offset - .1, MAP2GRAPH(map, x1 + .80),
		      *offset + .9, MAP2GRAPH(map, stop)) ;
      graphBoxEnd() ;
      color = YELLOW ;
      if (look->mode == CDNA && lane->previousEnd)
	color = PALEMAGENTA ;
      graphBoxDraw (box, BLACK, color) ;
      graphBoxSetPick (box, FALSE) ;
    }

  if (!upSequence && stop > u2) 
    { int box = graphBoxStart() ;
      graphRectangle (*offset - .1, MAP2GRAPH(map, u2),
		      *offset + .9, MAP2GRAPH(map, stop)) ;
      graphBoxEnd() ;
      color = YELLOW ;
      if (look->mode == CDNA && lane->nextTop)
	color = PALEMAGENTA ;
      graphBoxDraw (box, BLACK, color) ;
      graphBoxSetPick (box, FALSE) ;
    }
  if (upSequence && start < u2) 
    { int box = graphBoxStart() ;
      graphRectangle (*offset - .1, MAP2GRAPH(map, start - .8),
		      *offset + .9, MAP2GRAPH(map, u2 + .85)) ;
      graphBoxEnd() ;
      color = YELLOW ;
      if (look->mode == CDNA && lane->nextTop)
	color = PALEMAGENTA ; 
      graphBoxDraw (box, BLACK, color) ;
      graphBoxSetPick (box, FALSE) ;
      u1++ ;
    }

/* show vector clipping */
  if (lane->dna && lane->clipEnd == arrayMax(lane->dna))
    v2 = lane->clipEnd ;
  if (lane->dna && lane->clipTop == 0)
    v1 = 0 ;
  if (v1 >= 0 || v2 > 0)
    { float oldW = graphLinewidth(.3) ;
      w1 = w2 = -1 ;
      if (!lane->upSequence)
	{ if (v1 > 0) w1 = lane->x1 + v1 - lane->clipTop ;
	if (v2 > 0) w2 = lane->x2 + v2 - lane->clipEnd /* - 1 */ ;
	}
      else
	{ if (v1 >= 0) w1 = lane->x1 - v1 + lane->clipTop ;
	if (v2 > 0) w2 = lane->x2 - v2 + lane->clipEnd + 1 ;
	}
      
      if (!upSequence && start < w1) 
	{ int box = graphBoxStart() ;
	  graphRectangle (*offset - .1, MAP2GRAPH(map, start - .8),
			  *offset + .9, MAP2GRAPH(map, w1 - .15)) ;
	  graphBoxEnd() ;
	  graphBoxDraw (box, BLACK, LIGHTBLUE) ;
	  graphBoxSetPick (box, FALSE) ;
	}
      if (w1 > 0 && upSequence && stop > w1) 
	{ int box = graphBoxStart() ;
	  graphRectangle (*offset - .1, MAP2GRAPH(map, w1 + .85),
			  *offset + .9, MAP2GRAPH(map, stop)) ;
	  graphBoxEnd() ;
	  graphBoxDraw (box, BLACK, LIGHTBLUE) ;
	  graphBoxSetPick (box, FALSE) ;
	}
      
      if (w2 > 0 && !upSequence && stop >= w2) 
	{ int box = graphBoxStart() ;
	  graphRectangle (*offset - .1, MAP2GRAPH(map, w2),
			  *offset + .9, MAP2GRAPH(map, stop)) ;
	  graphBoxEnd() ;
	  graphBoxDraw (box, BLACK, PALEORANGE) ;
	  graphBoxSetPick (box, FALSE) ;
	  /* if (u2 == w2) u3-- ; */
	}
      if (w2 > 0 && upSequence && start <= w2) 
	{ int box = graphBoxStart() ;
	  graphRectangle (*offset - .1, MAP2GRAPH(map, start - .8),
			  *offset + .9, MAP2GRAPH(map, w2 - .15)) ;
	  graphBoxEnd() ;
	  graphBoxDraw (box, BLACK, PALEORANGE) ;
	  graphBoxSetPick (box, FALSE) ;
	  if (u1 < w2) u1++ ;
	}
      graphLinewidth(oldW) ;
    }

/********/

  buf[0] ='.' ; buf[1] = 0 ;

  if (upSequence) u3++ ;
  for (i = u1 ; i < u3 ; i++)
    graphText(buf, x, MAP2GRAPH(map, i - .2)) ;

  if (!arrayExists(a) || !arrayMax(a))
    goto done;
    

  ep = arrp(a, 0, A_ERR)  - 1 ;
  i = arrayMax(a) ; nerr = 0 ;
  while (ep++, i--)
    { c = dnaDecodeChar[upSequence ? 
			complementBase(ep->baseShort) : (int)ep->baseShort] ;
      if (c == 'w') c = 'X' ;
      ddx = 0 ; ddy = 0 ;
      switch(ep->type)
	{
	case ERREUR: 
	  color = RED ; 
	  c = ace_lower(c) ;
	  break ;
	case TROU:
	  color = BLUE ; 
	  c = '*' ;
	  break ;
	case TROU_DOUBLE:
	  color = BLUE ; 
	  c = '#' ;
	  break ;
	case TROU_TRIPLE:
	  color = BLUE ; 
	  c = '#' ;
	  break ;
	case INSERTION: /*_AVANT: */
	  color = BLUE ;
	  c = ace_upper(c) ;
	  ddx = -.5 ; ddy = -.5 ;
	  break ;
/*	case INSERTION_APRES:
	  color = BLUE ;
	  c = ace_upper(c) ;
	  ddx = .5 ; ddy = .5 ;
	  break ;  */
	case INSERTION_DOUBLE: /*_AVANT: */
	  color = ORANGE ;     /* BLUE ; */
	  c = 'Y' ;
	  ddx = -.5 ; ddy = -.5 ;
	  break ;
/*	case INSERTION_DOUBLE_APRES:
	  color = BLUE ;
	  c = 'X' ;
	  ddx = .5 ; ddy = .5 ;
	  break ; */
	case INSERTION_TRIPLE: /*_AVANT: */
	  color = ORANGE ;     /* BLUE ; */
	  c = 'Z' ;
	  ddx = -.5 ; ddy = -.5 ;
	  break ;
/*	case INSERTION_TRIPLE_APRES:
	  color = BLUE ;
	  c = 'z' ;
	  ddx = .5 ; ddy = .5 ;
	  break ; */

	case AMBIGUE: 
	  color = GREEN ; 
	  c = ace_lower(c) ;
	  break ;
	case TYPE80:
	  break ;
	}
      if (upSequence)
	{ ddx *= -1 ; /* ddy *= -1 ; */ }
      nn = ep->iLong ; 
      y = MAP2GRAPH(map, nn) ;
      buf[0] = c ;
      if (nn >= start && nn < stop &&
	  ep->iShort >= lane->clipTop &&
	  ep->iShort < lane->clipEnd)
	{ graphText(buf, x + ddx, y + ddy) ; nerr++ ;
	}
    }

/* this is a duplicate ?
  if (lane->scf == 2)
    { int box ; 
      float oldW = graphLinewidth(.3) ;
      lane->boxName = box = graphBoxStart () ;
      y = look->topMargin ;
      x = *offset - .2 ;

      if (lane->upSequence)
	{ graphLine (x, y, x + 2 , y) ;
	  graphLine (x, y, x + 1 , y - 1) ;
	  graphLine (x + 1 , y - 1, x + 2 , y) ;
	}
      else
	{ graphLine (x, y  - 1, x + 2 , y - 1) ;
	  graphLine (x, y - 1, x + 1 , y) ;
	  graphLine (x + 1 , y, x + 2 , y - 1) ;
	}
      graphLinewidth(oldW) ;
      graphBoxEnd () ;
      graphBoxDraw (box, BLACK, 
		    look->activeLane &&
		    lane == look->activeLane ?
		    LIGHTBLUE : TRANSPARENT) ;

    }
  else
    lane->boxName = 0 ;
*/
 done:
   graphBoxEnd () ; 
/*    
  if (lane->boxName)
    { if (lane->upSequence)
	graphText (messprintf("%d", u3), x +.1, y - 2.3) ; 
      else
	graphText (messprintf("%d", u1 + 1), x +.1, y - 2.3) ; 
    }
*/
  lane->nerr = nerr ;
  *offset += lane->scf >= 3 ? 1.4 : 2.4 ;
}

/*****************************************************************/

void traceEds (LOOK look, float *offset)
{ int k = 0, i = arrayMax(look->lanes) ;
  LANE *lane = i ? arrp(look->lanes,0,LANE) - 1 : 0 ;
  BOOL old = FALSE ;
  float mag = look->map->mag ;
  Array a ;


  a = look->laneShown = arrayReCreate (look->laneShown, 20, LANE*) ;

  look->map->mag *= look->edMag ;
  *offset += 1 ; /* after consensus */

  while(lane++, i--)
    { lane->isShown = FALSE ; lane->greenBox = 0 ;
      if (isLaneShown (look, lane))
	{ CHECKSEQ ;
	  if (lane->isAligned || lane->scf >= 3)
	    array (a, k++, LANE*) = lane ;
	}
    }

  if (showSelect)
    selectBestTraces (look) ; /* will fineTune */

  if (sortNeeded)
    arraySort(look->lanes, laneGlobalOrder) ;
  sortNeeded = FALSE ;
  arraySort (a, lanePlotOrder) ;

  if (look->hideDots < 2) /* auto choice, count and choose */
    { 
      int nn = 0 ;
      for (i = 0 ; i < k ; i++)
	{
	  lane = arr (a, i, LANE*) ;
	  if (!lane->isAligned && lane->scf >= 3)
	    nn++ ;
	}
      look->hideDots = nn > 6 ? 1 : 0 ;	
    }

  if (! (look->hideDots & 0x1))
    for (i = 0 ; i < k ; i++)
      {
	lane = arr (a, i, LANE*) ;
	if (old != lane->upSequence)
	  *offset += 1 ; /* betweeen up and down */
	
	old = lane->upSequence ;
	if (!lane->isAligned && lane->scf >= 3)
	  traceEdLaneBases (look, lane, offset, look->edMin,
			    look->edMax, lane->scf, FALSE) ;
      }
  
  look->map->mag = mag ;
  *offset += 1 ; /* see traceEdCoords */
}

/*****************************************************************/
/*************************** Header ******************************/


/**********************************************************/

static void exportLane0(void)
{ TRACE *bp[4], *sp ;
  int u1, u2, u3, u4, z ;
  LANE *lane ;
  FILE *ff = 0 ;
  int i, j, n ;
  Array a ;  
  char *cp ;
  TRACELOOKGET("export") ;

  if (!look->lanes ||
      !arrayMax(look->lanes))
    return ;

  lane = arrp(look->lanes, 0, LANE) ;
  traceGetLane (look, lane) ;
  CHECKSEQ ;
  if (lane->scf < 4)
    return ;
  bp[0] = lane->seq->traceA ;
  bp[1] = lane->seq->traceG ;
  bp[2] = lane->seq->traceC ;
  bp[3] = lane->seq->traceT ;
  
  if ((ff = filopen (messprintf("%s", name(lane->key)), "trx", "w")))
    { fprintf (ff, "       pos     A     T     G     C\n\n") ;
      for (n=0 ; n < lane->maxPos ; n++)
	{ fprintf(ff, "\n %6d: ", n) ;
	  j = *(bp[0] + n) ;
	  j = *(bp[3] + n) ;
	  j = *(bp[1] + n) ;
	  j = *(bp[2] + n) ;
	  fprintf (ff, "%6d", j) ;
	}
      filclose (ff) ;
    }
  if ((ff = filopen (messprintf("%s", name(lane->key)), "tr0", "w")))
    { i = 4 ;
      while (i--)
	{ cp = (char*) bp[i] ;
	  n = lane->maxPos * sizeof (TRACE) ;
	  fwrite (cp, n, 1, ff) ;
	}
      filclose (ff) ;
    }
  if ((ff = filopen (messprintf("%s", name(lane->key)), "tr1", "w")))
    { i = 4 ;
      while (i--)
	{ a = arrayCreate (lane->maxPos, char) ;
	  array (a, lane->maxPos - 1 , char) = 0 ;
	  cp = arrp (a, 0, char) ;
	  sp = bp[i] ;
	  u1 = u2 = u3 = 0 ;
	  for (j=0 ; j < lane->maxPos ; cp++, sp++, j++)
	    { z = u1 ;
	      *cp =  (char) ((*sp - z) & 255) ;
	      u3 = u2 ; u2 = u1 ; u1 = *sp ;
	    }
	  n = arrayMax (a) ;  cp = arrp (a, 0, char) ;
	  fwrite (cp, n, 1, ff) ;
	  arrayDestroy (a) ;
	}
      filclose (ff) ;
    }
  if ((ff = filopen (messprintf("%s", name(lane->key)), "tr2", "w")))
    { i = 4 ;
      while (i--)
	{ a = arrayCreate (lane->maxPos, char) ;
	  array (a, lane->maxPos - 1 , char) = 0 ;
	  cp = arrp (a, 0, char) ;
	  sp = bp[i] ;
	  u1 = u2 = u3 = 0 ;
	  for (j=0 ; j < lane->maxPos ; cp++, sp++, j++)
	    { z = 2*u1 - u2 ;
	      if (z < 0) z = 0 ;
	      if (z > 255) z = 255 ;
	      *cp =  (char) ((*sp - z) & 255) ;
	      u3 = u2 ; u2 = u1 ; u1 = *sp ;
	    }
	  n = arrayMax (a) ;  cp = arrp (a, 0, char) ;
	  fwrite (cp, n, 1, ff) ;
	  arrayDestroy (a) ;
	}
      filclose (ff) ;
    }
  if ((ff = filopen (messprintf("%s", name(lane->key)), "tr3", "w")))
    { i = 4 ;
      while (i--)
	{ a = arrayCreate (lane->maxPos, char) ;
	  array (a, lane->maxPos - 1 , char) = 0 ;
	  cp = arrp (a, 0, char) ;
	  sp = bp[i] ;
	  u1 = u2 = u3 = 0 ;
	  for (j=0 ; j < lane->maxPos ; cp++, sp++, j++)
	    { z = 3*u1 - 3*u2 + u3 ;
	      if (z > 255) z = 255 ;
	      if (z < 0) z = 0 ;
	      *cp =  (char) ((*sp - z) & 255) ;
	      u3 = u2 ; u2 = u1 ; u1 = *sp ;
	    }
	  n = arrayMax (a) ;  cp = arrp (a, 0, char) ;
	  fwrite (cp, n, 1, ff) ;
	  arrayDestroy (a) ;
	}
      filclose (ff) ;
    }
  if ((ff = filopen (messprintf("%s", name(lane->key)), "tr4", "w")))
    { i = 4 ;
      while (i--)
	{ a = arrayCreate (lane->maxPos, char) ;
	  array (a, lane->maxPos - 1 , char) = 0 ;
	  cp = arrp (a, 0, char) ;
	  sp = bp[i] ;
	  u1 = u2 = u3 = u4 = 0 ;
	  z = 4*u1 - 6*u2 + 4*u3 - u4;
	  for (j=0 ; j < lane->maxPos ; cp++, sp++, j++)
	    { z = 3*u1 - 3*u2 + u3 ;
	      if (z > 255) z = 255 ;
	      if (z < 0) z = 0 ;
	      *cp =  (char) ((*sp - z) & 255) ;
	      u4 = u3 ; u3 = u2 ; u2 = u1 ; u1 = *sp ;
	    }
	  n = arrayMax (a) ;  cp = arrp (a, 0, char) ;
	  fwrite (cp, n, 1, ff) ;
	  arrayDestroy (a) ;
	}
      filclose (ff) ;
    }
} 


/**********************************************************/

static void traceHisto (void);
static void baseCallPatchButton (void) ;

static FREEOPT showMenu[] = 
{ {8, "Show What"},
  {SHOW_ALL, "All Traces"},
  {HIDE_UP, "Hide Up Traces"},
  {HIDE_DOWN, "Hide Down Traces"},
  {HIDE_POOR, "Hide Poor Traces"},
  {SHOW_BEST_ERR, "Best + Errors"},
  {HIDE_NON_BEST, "Best traces"},
  {HIDE_NON_FULL, "Complete CDS clones"},
  {SHOW_SPLICER, "Exact splice sites"}
} ;

static void dummyButton (void) { return ; }

static void nextChooser (KEY k, int box)
{ TRACELOOKGET("nextChooser") ;
  look->next = k ;

  switch (k)
    {
    case NEXT_ERROR: showDefaultMode = SHOW_BEST_ERR ; break ;
    case NEXT_PROBLEM: showDefaultMode = SHOW_BEST_ERR ; break ;
    case NEXT_COMPARE: showDefaultMode = HIDE_NON_BEST ; look->consensus = TRUE ; break ;
    case NEXT_CLIP: showDefaultMode = HIDE_NON_BEST ; break ;
    case NEXT_HOLE: showDefaultMode = HIDE_NON_BEST ; break ;
    case NEXT_TAG: showDefaultMode = SHOW_BEST_ERR ; break ;
    }
}

static void nnTrain (void)
{ TRACELOOKGET("nnTrain") ;

  nnAssemblyTrain (look) ;
}

static MENUOPT emptyMenu[] = {
  {menuSpacer,""},
  {0, 0} } ;

/*
static MENUOPT globalEditMenu[] = {
  {undoLastVectorClip,"Undo last vector clipping"}, 
  {0, 0} } ;
*/

static MENUOPT complexMenu[] = {
  {mapColControl, "Display Control"},
  {traceHisto, "Signal Histograms"},
  {exportLane0, "Export lane 0"},
  {nnTrain, "Train Neural-Net"},
  {cdnaStopHisto, "Stop Histo"},
  {0, 0} } ;

static MENUOPT fixMenu[] = {
  {newFixButton, "Fix Window"}, 
  {newBigFixButton, "Fix Contig"},
  /* {fixButton, "Fix Whole Contig"}, */
  {traceRealignAll, "Realign"},
  {baseCallPatchButton, "Auto Edit BaseCall"},
  {0, 0} } ;

#ifdef APPARENTLY_NOT_NEEDED
static void tagButton(void)
{ messout ("This is a menu. Please, use the right mouse button") ;
  return ;
}
#endif /* APPARENTLY_NOT_NEEDED */


void traceMapZoomIn (void)
{ 
  MAP map = currentMap("traceMapZoomIn");

  lastMag = map->mag *= 1.414 ;
  (map->draw) () ;

  return;
} /* traceMapZoomIn */


void traceMapZoomOut (void)
{
  MAP map = currentMap("traceMapZoomOut");

  lastMag = map->mag /= 1.414 ; 
  (map->draw)() ;

  return;
} /* traceMapZoomOut */

void selectButton (void)
{
   TRACELOOKGET("selectButton") ;

  look->select = TRUE ; 
  traceDraw (look) ;
}

void hideDotsButton (void)
{
   TRACELOOKGET("selectButton") ;

  if (look->hideDots & 0x1)
    look->hideDots = 2 ;
  else
    look->hideDots = 3 ;
  traceDraw (look) ;
}

void selectArrows (void)
{
   TRACELOOKGET("selectButton") ;

  look->select = TRUE ; 
  traceDraw (look) ;
}

void unselectButton (void)
{
   TRACELOOKGET("unselectButton") ;

  look->select = FALSE ; 
  traceDraw (look) ;
}

Array traceAddTagMenu = 0 ;

static MENUOPT buttonOpts[] = {
  {graphDestroy, "Quit"},
  {help, "Help"},
  {traceMapZoomIn, "Zoom In"},
  {traceMapZoomOut, "Zoom Out"},
  {graphPrint, "Print"},
  {newFixButton, "Fix"},
  /*  {tagButton, "Tag..."}, */
  {0, 0} } ;

/************************************************************/
/***************** Registered routines *********************/

static void traceKbd (int k)
{ 
  TRACELOOKGET("traceKbd") ;

  if (look->mode == CDNA)
    {
      switch (k)
	{
	case 'e': case RETURN_KEY: case SPACE_KEY:
	  acceptZoneEdits () ;
	  break ;
	case 'n':
	  nifyZoneEdits () ;
	  break ;
	case 'v':
	  vectorifyZoneEdits () ;
	  break ;
	case 'f':
	  newFixButton () ;
	  graphPop () ;
	  return ;
	case 'g':
	  {
	    Graph old = graphActive () ;
	    newFixButton () ;  
	    if (graphActivate (look->fMapGraph))
	      graphPop () ;
	    graphActivate (old) ;
	  }
	  return ;
	case LEFT_KEY: 
	  traceMapZoomIn() ;
	  break ;
	case RIGHT_KEY:
	  traceMapZoomOut() ;
	  break ;
	case UP_KEY:
	  look->map->centre -= 5 ;
	  look->map->centre = nextError (look, FALSE) ;
	  (look->map->draw) () ;
	  break ;
	case DOWN_KEY:
	  look->map->centre += 5 ;
	  look->map->centre = nextError (look, TRUE) ;
	  (look->map->draw) () ;
	  break ;
	default:
	  break ;
	}
      return ;
    }
}
                
/*************************************************************/

static void traceCentre (LOOK look, float from)
{ int ny ;
  
  graphFitBounds (&look->graphWidth,&look->graphHeight) ;
  halfGraphHeight = 0.5*(look->graphHeight - look->topMargin) ;
 
  /* int nx = look->graphWidth ;  */
  ny = look->graphHeight - look->topMargin - 5 ;

  look->map->centre = from ;
  if (!look->map->mag)
    look->map->mag = 
      ny / (look->mode == CDNA ? 50.0 : 22.0) ;
   
  if (look->map->mag <= 0)  /* safeguard */
    look->map->mag = 0.05 ;
}

/*************************************************************/

static void traceDraw (LOOK look)
{ float mag, x ;
  char *cp ;
  int box, min, max ;


  if (graphActivate(look->graph))
    graphPop() ;
  else
    return ;

  graphClear () ;
  if (!look->dna)
    { graphRedraw() ; graphDestroy () ; return ; }
  arrayDestroy (histoTrace0) ;
  arrayDestroy (histoTrace1) ;
  arrayDestroy (histoTrace2) ;
  arrayDestroy (histoTrace3) ;
  arrayDestroy (histoTrace4) ;

  histoTrace0 = arrayCreate (256, int) ;
  histoTrace1 = arrayCreate (256, int) ;
  histoTrace2 = arrayCreate (256, int) ;
  histoTrace3 = arrayCreate (256, int) ;
  histoTrace4 = arrayCreate (256, int) ;

  genex = 0 ;
  dragFast = -1 ;
  look->topMargin = 5 ; /* see traceGene */

  traceCentre (look, look->map->centre) ;
  look->box2seg = arrayReCreate (look->box2seg, 256, int) ;

  look->hints = arrayReCreate (look->hints, 100, HINT) ;
  look->baseBoxes = arrayReCreate (look->baseBoxes, 100, HINT) ;
  look->wmin = GRAPH2MAP(look->map, look->topMargin) + .8 ;
  look->wmax = GRAPH2MAP(look->map, look->graphHeight) + .9 ;

  mag = look->map->mag ;
  look->map->thumb.box = 0 ;
  look->map->mag = 1.2 ;
  look->edMag = look->map->mag / mag ;

  min = GRAPH2MAP(look->map, look->topMargin) ;
  max = GRAPH2MAP(look->map, look->graphHeight) ;

  if (min < 0) min = 0 ;
  if (max > arrayMax(look->dna)) max = arrayMax(look->dna) ;
  
  if (max < 0) max = 0 ;
  if (min > arrayMax(look->dna)) min = arrayMax(look->dna) ;
  
  look->edMin = look->map->min = min ;
  look->edMax = look->map->max = max ;

  look->map->mag = mag ;
  graphMenu (emptyMenu) ;

  look->summaryBox = graphBoxStart () ;
  graphTextPtr (look->summary, look->graphWidth - 6, 2, 4) ;
  graphBoxEnd() ;
  graphBoxDraw (look->summaryBox, BLACK, LIGHTBLUE) ;

  look->activeBox = 0 ;
  mapDrawColumns (look->map) ;

  box = graphButtons(buttonOpts, 1,1.21,look->graphWidth) ;
  graphBoxMenu (box + 4, complexMenu) ;
  if (look->mode != CDNA)
    graphBoxMenu (box + 5, fixMenu) ;
  /*
    graphBoxFreeMenu (box + 6, (FreeMenuFunction) traceTagger,
		    traceTagMenu) ;
		    */
  x = 56.0 ;
  box = graphButton ("Show...", showButton, x, .2) ;
  graphBoxFreeMenu (box, (FreeMenuFunction) showChooser, showMenu) ;
  cp = freekey2text (look->hide, showMenu) ;
  graphText (cp, x + 8.5, .2) ;
  if (look->mode != CDNA)
    {
      box = graphButton ("Next...", dummyButton, x, 1.30) ;
      graphBoxFreeMenu (box, (FreeMenuFunction) nextChooser, nextMenu) ;
      cp = freekey2text (look->next, nextMenu) ;
      graphText (cp, x + 8.5, 1.32) ;
    }
  else
    {
      if (look->select)
	graphButton ("Unselect", unselectButton, x, 1.30) ;
      else
	graphButton ("Select", selectButton, x, 1.30) ;
    }
  if (look->hideDots & 0x1)
    graphButton ("Show dots", hideDotsButton, x + 10, 1.30) ;
  else
    graphButton ("Hide dots", hideDotsButton, x + 10, 1.30) ;

  switch (look->hideUp)
    {
    case HIDE_UP:
      graphText ("and Hide Up traces", x + 25, .2) ;
      break ;
    case HIDE_DOWN:
      graphText ("and Hide Down traces", x + 25, .2) ;
      break ;
    case HIDE_POOR:
      graphText ("and Hide poor traces", x + 25, .2) ;
      break ;
    case HIDE_NON_FULL:
      graphText ("and Hide Non full CDS", x + 25, .2) ;
      break ;
    }
  graphRedraw() ;
}

static void traceHisto (void)
{
  if (histoTrace0 && arrayMax(histoTrace0))
    plotHisto ("Distribution des traces", histoTrace0) ;

  if (histoTrace1 && arrayMax(histoTrace1))
    plotHisto ("Distribution des derivees premieres", histoTrace1) ;

  if (histoTrace2 && arrayMax(histoTrace2))
    plotHisto ("Distribution des derivees secondes", histoTrace2) ;

  if (histoTrace3 && arrayMax(histoTrace3))
    plotHisto ("Distribution des derivees troisiemes", histoTrace3) ;

  if (histoTrace4 && arrayMax(histoTrace4))
    plotHisto ("Distribution des derivees quatriemes", histoTrace4) ;

   histoTrace0 = 0 ; /* destroys in plot */
   histoTrace1 = 0 ; /* destroys in plot */
   histoTrace2 = 0 ; /* destroys in plot */
   histoTrace3 = 0 ; /* destroys in plot */
   histoTrace4 = 0 ; /* destroys in plot */
}

/************** little function to act as a callback ***********/

static void drawVoid (void)
{
  TRACELOOKGET("drawVoid") ;
  traceDraw (look) ;
  pickRememberDisplaySize (DtMultiTrace) ;
}

/*****************************************************************/
/*********************************************************/
/*********************************************************/
  /* needed after some of the edition operations in fmaptrace */
void traceGraphDestroy (void)
{ Graph old = graphActive() ;
  if (graphActivate (traceGraph))
    graphDestroy () ;
  graphActivate (old) ;
}

BOOL multiTraceDisplay (KEY key, KEY from, BOOL isOldGraph)
{ KEYSET ks ;
  Array linkDna = 0 ;
  LOOK look ;
  KEY dnaKey, linkKey ;
  Graph fMapGraph ;
  int sens, p1, p2 ;
  OBJ Link = 0 ;
  void *fLook = 0 ;

  fMapActive(&linkDna, 0, &linkKey, &fLook) ;
  fMapGraph = fMapActiveGraph () ;
  ks = queryKey(key, "{ >Assembled_from} $| {>Aligned}") ;  /*; SCF_File*/
  if (!keySetMax(ks))
    keySet(ks, 0) = key ;

  if (!traceLook)
    { traceLook = (LOOK) messalloc (sizeof(struct LookStruct)) ;
      traceLook->magic = &TRACELOOK_MAGIC ;
    }
  look = traceLook ;
  look->activeBox = 0 ;
  look->mode = acemblyMode ;
  if (keyFindTag (key, str2tag("Is_Reference")))
     look->mode = MUTANT ;
  baseEditMenu = look->mode == CDNA ? &baseDoEditcDNAMenu[0] : &baseDoEditMenu[0];
  showSelect = TRUE ;
  if (!graphActivate (traceGraph))
    { traceGraph = look->graph = displayCreate (DtMultiTrace) ;
      if (!look->graph)
	return FALSE ;
      
      graphAssociate (&TRACELOOK_ASSOC, look) ;
      graphAssociate (&MAP2LOOK_ASSOC, look) ;
      look->box2seg = arrayCreate (64,int) ;

#ifdef JUNK
      
      look->map = mapCreate2 (colInfo, drawVoid) ;
      static MapColRec2 colInfo[] = {
	{  1.1, TRUE, "Locator", mapShowLocator},
	{  1.5, TRUE, "EditorCoor", traceEdCoords},
	{  2.0, TRUE, "EditorConsensus", traceEdConsensus},
	{  3.0, TRUE, "Editor", traceEds},
	{  4.0, TRUE, "Gene", traceGene},
	{  5.0, TRUE,  "Scale", traceScale},
	{  6.0, TRUE, "Consensus", traceConsensus},
	{  7.0, TRUE, "Traces", tracePlotAllLanes},
	{ 0.0 , 0, 0, 0}			/* terminator */
      } ;
#endif
      look->map = mapCreate (drawVoid) ;
      mapInsertCol (look->map,1.1, TRUE, "Locator", mapShowLocator) ;
      mapInsertCol (look->map,1.5, TRUE, "EditorCoor", traceEdCoords) ;
      mapInsertCol (look->map,2.0, TRUE, "EditorConsensus", traceEdConsensus) ;
      mapInsertCol (look->map,3.0, TRUE, "Editor", traceEds) ;
      mapInsertCol (look->map,4.0, TRUE, "Gene", traceGene) ;
      mapInsertCol (look->map,5.0, TRUE,  "Scale", traceScale) ;
      mapInsertCol (look->map,6.0, TRUE, "Consensus", traceConsensus) ;
      mapInsertCol (look->map,7.0, TRUE, "Traces", tracePlotAllLanes) ;

      look->map->mag = lastMag ;
      look->next = NEXT_PROBLEM ;
      look->hide = SHOW_BEST_ERR ;
      look->hideUp = 0 ;
      showDefaultMode = SHOW_BEST_ERR ;
      graphRegister (RESIZE, drawVoid) ;
      graphRegister (DESTROY, traceDestroy) ;
      graphRegister (PICK, traceLeftDown) ;
      graphRegister (KEYBOARD, traceKbd) ;
      graphRegister (MIDDLE_DOWN, traceMiddleDown) ;
      graphHelp ("Trace_Editor") ;
    }

  look->map->min = 0 ;
  
  if (fLook) /* if i registered from fmap */
    { look->fMapLook = fLook  ;
      look->fMapGraph = fMapGraph ;
      look->link = linkKey ;
    }
  Link = bsCreate (look->link) ;
  if (Link && 
      (bsFindKey (Link, _Subsequence, key) || bsFindKey (Link, _Transcribed_gene, key)) &&
      bsGetData (Link, _bsRight, _Int, &p1) && 
      bsGetData (Link, _bsRight, _Int, &p2))
    sens = p2 - p1 ;
  else 
    sens = 1 ;
  bsDestroy (Link) ;
  traceIsPoorLane (0, 0) ;
  if (fMapcDNAReferenceDna || key != look->key || look->sens * sens < 0)
    { traceDoDestroy (look) ;
      look->key = key ;
      look->sens = sens ;
      look->hide = SHOW_BEST_ERR ; 
      if (key != look->key || look->sens * sens < 0)
	look->hideUp = 0 ;
      look->next = NEXT_PROBLEM ;
      showDefaultMode = SHOW_BEST_ERR ;
      showChosenMode = SHOW_ALL ;
      if (fMapcDNAReferenceDna)
	{ 
	  sens = 1 ;
	  look->dnaKey = 0 ;
	  look->mode = CDNA ;
	  baseEditMenu = &baseDoEditcDNAMenu[0];
	  look->dna = fMapcDNAReferenceDna ;
	  look->hits = fMapcDNAReferenceHits ;
	  look->gene = fMapcDNAReferenceGene ;
	  look->origin = fMapcDNAReferenceOrigin ;
	  fMapcDNAReferenceHits = 0 ; fMapcDNAReferenceDna = 0 ; /* now mine */
	  fMapcDNAReferenceGene = 0 ;
	  fMapcDNAReferenceOrigin = 0 ;
	  look->hide = SHOW_BEST_ERR ;  
	  look->hideUp = HIDE_POOR ;
	  look->next = NEXT_ERROR ;
	  showDefaultMode = SHOW_BEST_ERR ;
	  showChosenMode = SHOW_ALL ;
	}
      else if (!lexReClass (look->key, &dnaKey, _VDNA) ||
	       !(look->dnaKey = dnaKey) ||
	       !(look->dna = dnaGet(dnaKey)))
	{ look->dnaKey = 0 ; /* never edit somebody else dna */
	  if (!(look->dna = dnaGet (look->key))) /* try to construct the dna */
	    { messout ("no refernce dna") ; graphDestroy () ; return FALSE ; }
	}

      look->dnaR = dnaCopy (look->dna) ;
      if (sens > 0)
	reverseComplement (look->dnaR) ;
      else /* contig a l'envers dans le link */
	reverseComplement (look->dna) ;
      look->map->max = arrayMax (look->dna) ;
      look->wmin = 0 ;
      look->wmax = look->map->max ;
      if (look->hits)
	traceMakeHitLanes (look, sens) ;
      else
	traceMakeLanes (look, sens) ;
      look->hideDots = 0 ;
    }
  graphRetitle (messprintf("Traces: %s", name(key))) ;
  look->map->centre = from + look->origin ;
  if (fMapcDNAReferenceEst)
    {
      LANE *lane ;
      int i = arrayMax(look->lanes) ;
      while (i--)
	{
	  lane = arrp(look->lanes, i, LANE) ;
	  if (lane->key == fMapcDNAReferenceEst)
	    { look->activeLane = lane ; break ; }
	}
    }    
  traceDraw (look) ;

  return TRUE ;
}

/*********************************************************/
/*********************************************************/

/* Adapt precisely by sliding an heptamer */

static BOOL fineTuneDx1(Array dna, LANE *lane, int u0, int *zp) 
{  char *cp, *cp0 ;
   Array bc = lane->baseCall ;
   int 
     dy, s, j0, x0 = *zp, 
     max = arrayMax(bc),
     delta = 2, 
     delta1 = 2,   /* align on 5-mers */
     i, n ;
   BOOL up = lane->upSequence ;
   BASECALL *bb, *bb0, *bestbb ;
   char mybase[] = { A_, G_, C_, T_ } ;

   if (u0 < 0 || u0 >= arrayMax(dna))
     return FALSE ;
   cp0 = arrp(dna, u0, char) ;
   s = -10 ;

    /* search for the last base left of x0 */
   j0 = 0 ; 
   bb0 = arrp (bc, j0, BASECALL) ;
   while (j0 < max && bb0->x < x0)
     { j0 += 20 ; bb0 += 20 ; }
   if (j0 > 0)
     { bb0 -= 20 ; j0 -= 20 ; }
   while (j0 < max && bb0->x < x0)
     { j0 ++ ; bb0++ ;}
   if (j0 > 0)
     { bb0-- ; j0-- ;}

   bestbb = bb0 ;
   if (!up)
     { i = bb0 - arrp (bc, 0, BASECALL) ;
       if (i < delta + delta1 || i + delta + delta1 >= max)
	 return FALSE ;
       dy = 0 ; /* hope fineTuneDx did a good job */
       { cp = cp0 - delta1 ;
	 bb = bb0 - delta1 + dy ;
	 i = 2*delta1 + 1 ; n = 0 ;
	 while (i--)
	   if (mybase[(int)((bb++)->t)] == *cp++)
	     n++ ;
	 if (n > s)
	   { s = n ;
	     bestbb = bb0 + dy ;
	   }
       }
       if (s < 2*delta1 + 1)
	 for (dy = -delta ; dy < delta ; dy++)
	   { cp = cp0 - delta1 ;
	     bb = bb0 - delta1 + dy ;
	     i = 2*delta1 + 1 ; n = 0 ;
	     while (i--)
	       if (mybase[(int)((bb++)->t)] == *cp++)
		 n++ ;
	     if (n > s)
	       { s = n ;
		 bestbb = bb0 + dy ;
	       }
	   }
     }
   else
     { i = bb0 - arrp (bc, 0, BASECALL) ;
       if (i < delta + delta1 || i + delta + delta1 >= max)
	 return FALSE ; 
       dy = 0 ;
       { cp = cp0 + delta1 ;
	 bb = bb0 - delta1 + dy ;
	 i = 2*delta1 + 1 ; n = 0 ;
	 while (i--)
	   if (mybase[(int)3 - (bb++)->t] == *cp--)
	     n++ ;
	 if (n > s)
	   { s = n ;
	     bestbb = bb0 + dy ;
	   }
       }
       if (s < 2*delta1 + 1)
	 for (dy = -delta ; dy < delta ; dy++)
	   { cp = cp0 + delta1 ;
	     bb = bb0 - delta1 + dy ;
	     i = 2*delta1 + 1 ; n = 0 ;
	     while (i--)
	       if (mybase[(int)3 - (bb++)->t] == *cp--)
		 n++ ;
	     if (n > s)
	       { s = n ;
		 bestbb = bb0 + dy ;
	       }
	   }
     }

   if (s >= 2*delta1)  /* one error only */
     *zp = bestbb->x ;
   else
     *zp = bb0->x ;
   *zp += up ? +1 : +1 ;
   return TRUE ;
}

/**********/

static int fineTuneDx(LOOK look, LANE *lane, int u0, int *zp) 
{  
  A_ERR *ep ;
  Array dna = look->dna ;
  int x0, di = 0, n ;
  int
    x1 = lane->x1, x2 = lane->x2 , x3 = lane->x3  ;
 
  if (!look->showClip)
    x3 = x2 ;

  if (!lane->upSequence &&
      (x3 < u0 || x1 >= u0))
    return FALSE ;
  if (lane->upSequence &&
      (x1 < u0 || x3 >= u0))
    return FALSE ;

  if (!lane->errArray)  /* happens after a fix */
    laneMakeErrArray (look, lane) ;
 
  n = arrayMax (lane->errArray) ;

#ifdef JUNK
  BASECALL *xp ;
  /* in case of perfect match, verify that the whole base call is not shifted */
  if (!n ||
      arrp(lane->errArray, 0, A_ERR)->iLong > u0 )  /* perfect match */
    {
      if (!lane->upSequence)
	x0 = u0 - lane->x1 + lane->clipTop ;
      else
	x0 = - u0 + lane->x2 + lane->clipEnd ;
      if (x0 < 200 && x0 > 100 && /* verifier l'alignement global des basecall */
	  lane->baseCall && arrayMax(lane->baseCall))
	{
	  int 
	    min = x0 - 20, max = x0 + 20, 
	    max = x0 + 20,
	    delta[] = {0,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,7,-7,8,-8,9,-9,10,-10} ;


	  xp = arrp(lane->baseCall, 0, BASECALL) ; 
	  i = arrayMax(lane->baseCall) ;
	  while(xp->x < min && i)
	    xp++ , i-- ;
	  while(i && xp->x < max)
	    { if (xp->t == k)
	      { y = MAP2GRAPH(map, xp->x + 1) ;
	      if (xp->flag & BC_LOW) ;
	      }
	    xp++ ; i-- ;
	    }

	  k1 = lane->upSequence ? 3 - k : k ;
	}	  
    }
#endif

  /* align using the edited data */
  if (!n ||
      arrp(lane->errArray, 0, A_ERR)->iLong > u0 + 30)  /* perfect match */
    {
      if (!lane->upSequence)
	x0 = u0 - lane->x1 + lane->clipTop ;
      else
	x0 = - u0 + lane->x2 + lane->clipEnd ;
    }
  else
    { ep = arrp(lane->errArray, 0, A_ERR) ;
      if (n-- > 1 && ep->iLong < u0) 
	{ while (ep++, n--)
	    if (ep->iLong > u0) break ;
	  if (ep > arrp(lane->errArray, 0, A_ERR)) 
	    ep-- ; /* last correct error, <= u0 */
	  else
	    {
	      if (!lane->upSequence)
		x0 = u0 - lane->x1 + lane->clipTop ;
	      else
		x0 = - u0 + lane->x2 + lane->clipEnd ;
	      goto ok ;
	    }
	}
      
      if (ep->iLong != u0)
	switch (ep->type)
	  {
	  case TROU: 
	    di = -1 ; break ;
	  case TROU_DOUBLE: 
	    di = -2 ; break ;
	  case INSERTION : /*_AVANT: */
	    di = 1 ; break ;
	  case INSERTION_DOUBLE: /* _AVANT: */
	    di = 2 ; break ;
	  default:
	    break ;
	  }
    
      if (!lane->upSequence)   
	x0 = di + u0 - ep->iLong + ep->iShort ;
      else
	x0 = - di + ep->iLong  - u0 + ep->iShort ;
    }
  
 ok:
  if (lane->scf >= 3)
    { if (x0 < 0 || !lane->basePos || x0 >= arrayMax(lane->basePos))
	return FALSE ;
      
      *zp = arr (lane->basePos, x0, short) ;
      
      /* now  adapt to the unedited sequence */
      if (lane->baseCall && !lane->hide)
	fineTuneDx1 (dna, lane, u0, zp) ; 
    }
  return TRUE ; 
}

/***********/

static void selectBestLane(LOOK look, LANE **lane1p) 
{ int i = arrayMax (look->lanes), nn = arrayMax (look->dna) ;
  LANE *lane = arrp (look->lanes, 0, LANE) - 1 ;

  while (lane++, i--)
    if (lane->scf >=3 && lane->nerr >= 0 && lane->nerr < nn)
      { nn = lane->nerr ; *lane1p = lane ; }
  if (!fineTune (look, lane1p, *lane1p))
    *lane1p = (LANE *)1 ;  
}

static BOOL fineTuneTrace (LOOK look, LANE **lane1p, LANE *lane, 
			   int x0, int u0, int u1, int u2)
{ int dx ;
  LANE *lane1 ;
  BOOL direct ;

  if (!*lane1p) 
    selectBestLane(look, lane1p) ;
  lane1 = *lane1p ;
  if (lane1 < (LANE *)2) return FALSE ;

  direct = (lane->upSequence == (*lane1p)->upSequence ? TRUE : FALSE ) ;
  if (lane1->seq && lane->seq &&
      baseCorrel (lane1->seq, lane1->t1, direct, lane->seq, x0, lane1->t2, 80, 1, &dx))
    { lane->ddy = lane1->ddy * (direct ? 1 : -1) ; 
      lane->dy =  look->map->centre - (x0  + dx) ;
      return TRUE ; /* ulrich */
    }
  else
    return FALSE ;
}

/***********/

static BOOL fineTune (LOOK look, LANE **lane1p, LANE *lane)
{ int 
    x0 = -100000, x1, x2, x01, x02, 
    ww = look->wmax - look->wmin,
    u1 = look->wmin + 1 ,
    u01 = look->wmin + ww/4,
    u2 = look->wmax - 1,
    u02 = look->wmax - ww/4,
    u0 = look->map->centre ;
  BOOL   r0, r1, r2, r01, r02 ;

  if (u1 + 5 > u2)
    return FALSE ;

  CHECKSEQ ; CHECKPOS ;
  if (lane->isAligned || lane->scf < 3)
    return FALSE ;

  r0 = fineTuneDx(look, lane, u0, &x0) ;
  if (r0 &&
      lane->nerr > 3 && /* difficult case */
      *lane1p != (LANE *)1  &&  /* bestlane did not fine tune */
      *lane1p != lane && 
      fineTuneTrace (look, lane1p, lane, x0, u0, u1, u2)) /* align lane against *lane1p */
    return TRUE ;
  x1 = x2 = x0 ;
  r1 = fineTuneDx(look, lane, u1, &x1) ;
  r2 = fineTuneDx(look, lane, u2, &x2) ;
  if (!r0)
    { return FALSE ; /* BUG, j'arrive pas a centre  */
      r01 = r02 = FALSE ;
      if (r1)
	r01 = fineTuneDx(look, lane, u01, &x01) ;
      else if (r2)
	r02 = fineTuneDx(look, lane, u02, &x02) ;
      if (!r01 && !r02)
	return FALSE ;
      else if (r01)
	x0 = x1 + (x01 - x1) * (u0 - u1)/(u01 - u1) ; 
      else if (r02)
	x0 = x2 - (x2 - x02) * (u2 - u0)/(u2 - u02) ;
    }
  lane->t1 = x0 ; lane->t2 = x2 > x1 ? x2 - x1 : x1 - x2 ;
 
  if (r0 && r1 && r2)
    { if (x2 == x1)
	return FALSE ;
      lane->ddy = (u2 - u1) / ((float)(x2 - x1)) ;
    }
  else if (r0 && r1)
    { if (x0 == x1)
	return FALSE ;
      lane->ddy = (u0 - u1) / ((float)(x0 - x1)) ;
      u2 = 2*u0 - u1 ;
    }
  else if (r0 && r2)
    { if (x0 == x2)
	return FALSE ;
      lane->ddy = (u0 - u2) / ((float)(x0 - x2)) ;
      u1 = 2*u0 - u2 ;
    }
  else
    return FALSE ;

  lane->dy = look->map->centre - x0 ;
  return TRUE ;
}

/***********/

static void newFixButton (void)
{ 
  extern BOOL fMapcDNAPickTrace (KEY gene, KEY est, int from);
  int z1, z2, dz, zmax, x ;
  KEY gene = 0, newGene ;
  KEYSET selKs = 0 ;
  void* fMapLook = 0 ;
  Graph fMapGraph = 0, traceGraph = 0 ;
  TRACELOOKGET("newFixButton") ;
  
  traceGraph = look->graph ;
  graphActivate (look->graph) ;
  switch (look->mode)
    { 
      case 0: break ;
    case CDNA:
      selKs = look->selectedKs ;
      look->selectedKs = 0 ;
      newGene = cDNARealignGene (look->gene, 0, 0, 0, FALSE, 0, 0, 0) ; /* locally, no repeats */
      if (look->fMapLook)
	{
	  Graph old = graphActive () ;
	  fMapPleaseRecompute (fMapLook = look->fMapLook) ;
	  if (graphActivate (look->graph))
	    tracePick (look->map->thumb.box, MAP2GRAPH (look->map, look->map->centre), 0) ;
	  gene = look->gene ; fMapGraph = look->fMapGraph ;
	  x = look->map->centre - look->origin ;
	  if (graphActivate (look->graph))
	    {
	      /* graphDestroy () ; */
	      traceDoDestroy (look) ; /* destroys look */
	      graphClear () ;
	    }
	  if (graphActivate (fMapGraph))
	    {
	      if (!fMapcDNAPickTrace (gene, 0, x))
		gene = 0 ; /* trace graph needed to reconstruct fixed gene */
	      graphActivate (old) ;
	    }
	  else if (graphActivate (old))
	    graphDestroy () ; /* trace graph needed to reconstruct fixed gene */
	}
      if (newGene != gene) /* gene changed name, coords are lost */
	graphDestroy () ;
      if (selKs &&
	  graphActivate (traceGraph) &&
	  graphAssFind (&TRACELOOK_ASSOC, &look))
	{ look->selectedKs = selKs ; selKs = 0 ; traceDraw (look) ; }

      keySetDestroy (selKs) ;
      return ;
    case MUTANT:
      break ; 
    default:
      messout("In the current mode, you cannot edit the consensus") ;
     return ;
    }

  if (arrayMax(look->lanes) == 1) /* single read, then redraw all */
    { multiTraceDisplay (look->key, 0, TRUE) ; return ;}

  dz = (look->wmax -  look->wmin)/2 ;
  zmax = arrayMax (look->dna) ;
  z1 = look->wmin - dz ; if (z1 < 0) z1 = 0 ;
  z2 = look->wmax + dz ; if (z2 > zmax ) z2 = zmax ;

  trackContig (look, z1, z2, FALSE) ;
  
  updateConsensus (look) ;
  traceIsPoorLane (0, 0) ;
  traceDraw (look) ;
}

/***********/

static void newBigFixButton (void)
{ int zmax ;
  TRACELOOKGET("newFixButton") ;
       
  if (look->mode)
    { messout("In the current, you cannot edit the consensus") ;
     return ;
    }

  if (arrayMax(look->lanes) == 1) /* single read, then redraw all */
    { multiTraceDisplay (look->key, 0, TRUE) ; return ;}

  zmax = arrayMax (look->dna) ;
 
  trackContig (look, 0, zmax, TRUE) ;
  
  updateConsensus (look) ;
  traceDraw (look) ;
  /*
  traceRealignAll () ; 
  traceRealignAll () ; 
  */
}

/***********/

static BOOL isProblem (LOOK look, LANE *ln1, int z, BOOL down)
{ int i, j, k, y, u1, u2 ;
  LANE *lane ;
  A_ERR *ep ;
  BOOL upClip = FALSE, downClip = FALSE ;
  static Array a, aUp = 0, aDown = 0 ;
  int nLanesUp = 0, nLanesDown = 0 ;
  int nVoisins = 0, nExactUp = 0, nExactDown = 0;
  
     /* comptons les sequences up et down couvrant la region */
  aUp = arrayReCreate (aUp, 12, LANE*) ;
  aDown = arrayReCreate (aDown, 12, LANE*) ;
  i = arrayMax(look->lanes) ;
  lane = arrp (look->lanes, 0, LANE) - 1 ;
  while (lane++, i--)
    { if (lane->upSequence)
	{ u1 = lane->x2 ; u2 = lane->x1 ; }
      else
	{ u1 = lane->x1 ; u2 = lane->x2 ; }

      if (u1 > z || u2 < z) 
	    continue ;

      CHECKLANE ; 
      if (lane->scf < 2)
	continue ;
      if (lane->x2 == z)
	{ if (look->next == NEXT_CLIP)
	    return TRUE ;
	  if (lane->upSequence)
	    upClip = TRUE ;
	  else
	    downClip = TRUE ;
	}
      if (lane->upSequence)
	array (aUp, nLanesUp++, LANE*) = lane ;
      else
	array (aDown, nLanesDown++, LANE*) = lane ;
    }

  if (look->next == NEXT_HOLE)
    { if ( (upClip && nLanesUp < 3) ||
	  (downClip && nLanesDown < 3) )
	return TRUE ;
      else
	return FALSE ;
    }

  if ((ln1->upSequence && nLanesUp < 3) ||
      (!ln1->upSequence && nLanesDown < 3))
    return TRUE ;
  k = nLanesDown ;
  while (k--)
    { lane = arr (aDown, k, LANE*) ;
      a = lane->errArray ;
      j = arrayMax(a) ;
      if (!j)
	continue ;
      ep = arrp(lane->errArray, 0, A_ERR) - 1 ;
      while(ep++, j--)
	{ y = ep->iLong ;
	  if (y > z - 3 && y < z + 3)
	    nVoisins++ ; 
	  if (y > z + 3)
	    break ;
	  if (y == z)
	    nExactDown++ ;
	}
    }
      
  k = nLanesUp ;
  while (k--)
    { lane = arr (aUp, k, LANE*) ;
      a = lane->errArray ;
      j = arrayMax(a) ;
      if (!j)
	continue ;
      ep = arrp(lane->errArray, 0, A_ERR) - 1 ;
      while(ep++, j--)
	{ y = ep->iLong ;
	  if (y > z - 3 && y < z + 3)
	    nVoisins++ ; 
	  if (y > z + 3)
	    break ;
	  if (y == z)
	    nExactUp++ ;
	}
    }

  if (nLanesUp < 2*nExactUp ||
      nLanesDown < 2*nExactDown ||
      (nLanesUp + nLanesDown) < 3*(nExactUp + nExactDown) ||
      2 * nVoisins > 3 *(nLanesUp + nLanesDown))
    { if (nLanesUp > 5 && nLanesUp - nExactUp > 2 &&
	  nLanesDown > 5 && nLanesDown - nExactDown > 2)
	return FALSE ;
      return TRUE ;
    }
  else
    return FALSE ;  
}

/***************************************************************/

static int nextError (LOOK look, BOOL down)
{ int i, j, i1, x = look->map->centre, y, ys, z, z0, 
    top, end, u1, u2, minTop, maxEnd ;
  LANE *lane, *ln1 = 0 ;
  A_ERR *ep ;
  KEY key ;
  Array a ;
  int dlane = down ? 1 : -1 , isEnd = 0 ;
  int nextType = look->next ;

  minTop = 0 ; maxEnd = arrayMax (look->dna) ;
  look->hide = showDefaultMode ;
  x += down ? 3 : -3 ;
 encore:
  z = z0 = x + (down ? 2 * arrayMax(look->dna) : - arrayMax(look->dna)) ;
  i = arrayMax(look->lanes) ;
  lane = down ?
    arrp (look->lanes, 0, LANE) - 1 : 
    arrp (look->lanes, arrayMax(look->lanes) - 1, LANE) + 1 ;
  while (lane += dlane, i--)
    { 
      if ((nextType == NEXT_COMPARE) ^ lane->isAligned)
	continue ;
      switch (look->hideUp)
	{
	case HIDE_UP:
	  if (lane->upSequence)
	    continue ;
	  break ;
	case HIDE_DOWN:
	  if (!lane->upSequence)
	    continue ;
	  break ;
	case HIDE_POOR:
	  if (lane->isPoor)
	    continue ;
	  break ;
	case HIDE_NON_FULL:
	  if (!lane->clone || !keyFindTag (lane->clone, str2tag("Complete_CDS_of")))
	    continue ;
	  break ;
	}
      
      if (lane->upSequence)
	{ u1 = lane->x2 ; u2 = lane->x1 ; }
      else
	{ u1 = lane->x1 ; u2 = lane->x2 ; }
      
      if (u1 < minTop) minTop = u1 ;
      if (u2 > maxEnd) maxEnd = u2 ;

      if ((down && (u1 >= z || u2 <= x)) ||
	  (!down && (u1 >= x || u2 <= z)))
	    continue ;

      if (look->mode == CDNA || nextType == NEXT_CLIP || nextType == NEXT_HOLE)
	{ if ( (down && x < lane->x2  && z > lane->x2  ) ||
	      (!down && x > lane->x2  && z < lane->x2  ))
	    { CHECKLANE ;
	      if (lane->dna && lane->clipEnd < arrayMax(lane->dna) - 5)
		{ 
		  BOOL forget = FALSE ;
		  int nl = arrayMax(look->lanes) ;
		  LANE *ll = arrp(look->lanes, 0, LANE) - 1 ;
		  while (ll++, nl--)
		    if (ll->key != lane->key && 
			ll->clone &&
			ll->clone == lane->clone &&
			((lane->upSequence && !ll->upSequence && lane->x2 < ll->x2) ||
			 (!lane->upSequence && ll->upSequence && lane->x2 > ll->x2)
			 ))
		      { forget = TRUE ; break ; }
		  if (!forget)
		    {
		      ln1 = lane ; z = lane->x2 ;
		      isEnd = lane->upSequence ? -2 : 2 ;
		    }
		}
	    }
	}

      if (look->mode == CDNA)
	{ 
	  if ( !lane->upSequence &&
	       (
		(down && x < lane->x1  && z > lane->x1  ) ||
		(!down && x > lane->x1  && z < lane->x1  ))
	       )
	    { 
	      CHECKLANE ;
	      if (lane->clipTop > lane->vectorTop ||  
		  ( (down && i == arrayMax(look->lanes) - 1 ) ||
		    (!down && i == 0)
		    ))		    
		{ ln1 = lane ; z = lane->x1 ;
		isEnd = lane->upSequence ? +2 : -2 ;
		}
	    }
	  if ( lane->upSequence &&
	       (
		(down && x < lane->x1  && z > lane->x1 ) ||
		(!down && x > lane->x1  && z < lane->x1  ))
	       )
	    { 
	      CHECKLANE ;
	      if (lane->clipTop > lane->vectorTop ||  
		  ( (!down && i == arrayMax(look->lanes) - 1 ) ||
		    (down && i == 0)
		    ))		    
		{ ln1 = lane ; z = lane->x1 ;
		isEnd = lane->upSequence ? +2 : -2 ;
		}
	    }
	}

      if (nextType == NEXT_CLIP || nextType == NEXT_HOLE)
	continue ;

      if (nextType == NEXT_TAG)
	{ if ( (down && x < lane->x2  && z > lane->x2  ) ||
	      (!down && x > lane->x2  && z < lane->x2  ))
	    { CHECKLANE ;
	      if (lane->hasTag != 1 && traceMakeLaneTags(lane))
		{ Array units = lane->tags ;
		  BSunit *u ;
		  for (i1 = 0 ; i1 < arrayMax (units) ; i1 += 4)
		    { u = arrp(units, i1, BSunit) ;
		      y  = (u[1].i + u[2].i) / 2 ;
		      if (lane->upSequence)
			y = lane->x1 - y + lane->clipTop ;
		      else
			y = lane->x1 + y - lane->clipTop ;
		      if ( (down && y > x && y < z) ||
			  (!down && y <  x && y > z))
			{ ln1 = lane ; z = y ; }
		    }
		}
	    }
	  continue ;
	}
      key = lane->clone ? lane->clone : lane->key ;
      if (look->select &&
	  keySetExists(look->selectedKs) &&
	  !keySetFind (look->selectedKs, key,0))
	continue ;


      /* next difference */ 
      CHECKLANE ; CHECKSEQ ; CHECKPOS ;
      if (lane->scf < 3 || !lane->base) /* data is lacking */
	continue ;
      a = lane->errArray ;
      j = arrayMax(a) ;
      if (!j)
	continue ;
      top = lane->clipTop ; end = lane->clipEnd ;
      ep = arrp(lane->errArray, 0, A_ERR) - 1 ;
      while(ep++, j--)
	{ 
	  if (look->mode == CDNA &&
	      ep->type == AMBIGUE &&
	      (BC_HAND & arr (lane->base, ep->iShort, char))     
	      )
	    continue ;
	  y = ep->iLong ;
	  if (down)
	    { if (y > x)
		{ ys = ep->iShort ;
		  if (y < z &&
		      ys >= top &&
		      ys < end)
		    { ln1 = lane ; z = y ; }
		  if (y >= z) break ;
		}
	    }
	  else
	    { if (y >= x)
		break ;
	      if (y > z &&
		  (ys = ep->iShort)  &&
		  ys >= top &&
		  ys < end)
		{ ln1 = lane ; z = y ; }
	    }
	}      
    }

  if (z == z0)
    {
      if (ln1) ln1->hide = FALSE ;
      return down ? maxEnd - 2 : minTop + 2 ;
    }

  x = z ;
  switch (look->next)
    {   
    case NEXT_PROBLEM:   /* ATTENTION general fall thru */
    case NEXT_HOLE: 
      if (!isProblem(look, ln1, z, down))
	goto encore ;
    case NEXT_CLIP:
    case NEXT_ERROR:
      if (!baseCallGetSeq (ln1))
	goto encore ; 
    case NEXT_COMPARE:
    case NEXT_TAG:
      break ;
    }

  if (ln1) ln1->hide = FALSE ;
  return z - 1 * isEnd ;
}

/***************************************************************/
/* dragging */
/*************************************************/
/************* middle button for thumb **********/

static double oldx, oldy, oldDy, xOrigin, yOrigin ;
static double dfx1, dfx2, dfy1, dfy2 ;
/* static int  dragFast = -1, dragBox = 0 ; declared at top of file */

static void traceMiddleDrag (double x, double y) 
{ MAP map ;
  TRACELOOKGET("traceMiddleDrag") ;

  map = look->map ;
  switch (dragFast)
    {
    case 0:
      graphXorLine (map->thumb.x, oldy, look->graphWidth, oldy) ;
      /* register zone edit in mid button
	 if (x > xOrigin + 3 || x < xOrigin - 3)
	dragFast = 2 ;
	*/
      break ;
    case 1:
      graphXorLine (0, oldy - oldDy, map->thumb.x, oldy - oldDy) ;
      graphXorLine (0, oldy + oldDy, map->thumb.x, oldy + oldDy) ;
      break ;
    case 2:
      graphXorLine (dfx1, dfy1, dfx1, dfy2) ;
      graphXorLine (dfx1, dfy1, dfx2, dfy1) ;
      graphXorLine (dfx2, dfy1, dfx2, dfy2) ;
      graphXorLine (dfx1, dfy2, dfx2, dfy2) ;
      break ;
    }

  oldy = y ;

  switch (dragFast)
    {
    case 0:
      graphXorLine (map->thumb.x, y, look->graphWidth, y) ;
      break ;
    case 1:
      oldDy *= exp ((x - oldx) / 25.) ;
      oldx = x ;
      graphXorLine (0, y - oldDy, map->thumb.x, y - oldDy) ;
      graphXorLine (0, y + oldDy, map->thumb.x, y + oldDy) ;
      break ;
    case 2:
      if (x > xOrigin) 
	{ dfx1 = xOrigin ; dfx2 = x ; }
      else
	{ dfx2 = xOrigin ; dfx1 = x ; }
      if (y > yOrigin) 
	{ dfy1 = yOrigin ; dfy2 = y ; }
      else
	{ dfy2 = yOrigin ; dfy1 = y ; }
      graphXorLine (dfx1, dfy1, dfx1, dfy2) ;
      graphXorLine (dfx1, dfy1, dfx2, dfy1) ;
      graphXorLine (dfx2, dfy1, dfx2, dfy2) ;
      graphXorLine (dfx1, dfy2, dfx2, dfy2) ;
      break ;
    }
}

static void traceMiddleUp (double x, double y) 
{ float x1,x2,y1,y2 ;
  MAP map ;
  LANE *lane ;
  int i ;
  BOOL doDraw = FALSE ;
  TRACELOOKGET("traceMiddleUp") ;

  map = look->map ;

  switch (dragFast)
    {
    case 0:
      if (y - yOrigin > 2)
	map->centre = nextError (look, TRUE) ;
      else if (y - yOrigin <  -2)
	map->centre = nextError (look, FALSE) ;
      else
	map->centre = GRAPH2MAP(map,y) ;
      doDraw = TRUE ;
      break ;
    case 1:
      graphBoxDim (map->thumb.box, &x1, &y1, &x2, &y2) ;
      map->mag *= (y2 - y1) / (2. * oldDy) ;
      map->centre = WHOLE2MAP(map, y) ;
      doDraw = TRUE ;
      break ;
    case 2:
      i = arrayMax(look->lanes) ;
      while (i--)
	{ 
	  lane = arrp (look->lanes, i, LANE) ;
	  if (!lane->isShown)
	    continue ;
	  if (dfx1 < lane->laneBaseCallBoxOffSet - .1 &&
	      dfx2 > lane->laneBaseCallBoxOffSet + 1.2)
	    { doDraw = TRUE ; doAcceptZoneEdits (TRUE, FALSE,FALSE, lane, dfy1, dfy2) ; }
	}
      if (doDraw)  traceDraw(look) ;
      break ;
    }
  dragBox = 0 ;
  dragFast = -1 ;
  (map->draw) () ;
}

static void traceMiddleDown (double x, double y) 
{ float x1,x2,y1,y2 ;
  MAP map ;
  TRACELOOKGET("traceMiddleDown") ;

  map = look->map ;
  if (map->thumb.box)
    graphBoxDim (map->thumb.box, &x1, &y1, &x2, &y2) ;
  oldDy = (y2 - y1) / 2. ;
  
  dragFast = (x < map->thumb.x) || x < genex ? 1 : 0  ;

  if (dragFast)
    { 
      if (genex)
	map->thumb.x = genex ;
      graphXorLine (0, y - oldDy, map->thumb.x, y - oldDy) ;
      graphXorLine (0, y + oldDy, map->thumb.x, y + oldDy) ;
    }
  else
    { graphXorLine (map->thumb.x, y, look->graphWidth, y) ;
    }
  oldx = xOrigin = x ;
  oldy = yOrigin = y ;

  graphRegister (MIDDLE_DRAG, traceMiddleDrag) ;	/* must redo */
  graphRegister (MIDDLE_UP, traceMiddleUp) ;
}

/****************************************************/
/************* left button for thumb **********/

static void traceLeftDrag (double x, double y) 
{ 
  TRACELOOKGET("traceLeftDrag") ;

  switch (dragFast)
    {
    case 0:
      break ;
    case 1:
      graphXorLine (dfx1, dfy1, dfx1, dfy2) ;
      graphXorLine (dfx1, dfy1, dfx2, dfy1) ;
      graphXorLine (dfx2, dfy1, dfx2, dfy2) ;
      graphXorLine (dfx1, dfy2, dfx2, dfy2) ;
      break ;
    default:
      return ;
    }

  if (x > xOrigin) 
    { dfx1 = xOrigin ; dfx2 = x ; }
  else
    { dfx2 = xOrigin ; dfx1 = x ; }
  if (y > yOrigin) 
    { dfy1 = yOrigin ; dfy2 = y ; }
  else
    { dfy2 = yOrigin ; dfy1 = y ; }

  if (dfx2 > dfx1 + 1 || dfy2 > dfy1 + 1)
    dragFast = 1 ;
  switch (dragFast)
    {
    case 0:
      break ;
    case 1:
      graphXorLine (dfx1, dfy1, dfx1, dfy2) ;
      graphXorLine (dfx1, dfy1, dfx2, dfy1) ;
      graphXorLine (dfx2, dfy1, dfx2, dfy2) ;
      graphXorLine (dfx1, dfy2, dfx2, dfy2) ;
      break ;
    }
}

static void traceLeftUp (double x, double y) 
{
  float x1,x2,y1,y2 ;
  LANE *lane ;
  int i ; 
  BOOL doDraw = FALSE ;
  TRACELOOKGET("traceLeftUp") ;


  switch (dragFast)
    {
    case 0:
      graphBoxDim (dragBox, &x1, &y1, &x2, &y2) ;
      tracePick (dragBox, xOrigin - x1, yOrigin - y1) ;
      dragBox = 0 ;
      /* doDraw = TRUE ; */
      break ;
    case 1:
      graphXorLine (dfx1, dfy1, dfx1, dfy2) ;
      graphXorLine (dfx1, dfy1, dfx2, dfy1) ;
      graphXorLine (dfx2, dfy1, dfx2, dfy2) ;
      graphXorLine (dfx1, dfy2, dfx2, dfy2) ;
      i = arrayMax(look->lanes) ;
      while (i--)
	{ 
	  lane = arrp (look->lanes, i, LANE) ;
	  if (!lane->isShown)
	    continue ;
	  if (dfx1 < lane->laneBaseCallBoxOffSet - .1 &&
	      dfx2 > lane->laneBaseCallBoxOffSet + 1.2)
	    { doDraw = TRUE ; doAcceptZoneEdits (TRUE, FALSE,FALSE, lane, dfy1, dfy2) ; }
	} 
      if (doDraw)  traceDraw(look) ;
      dragBox = 0 ;
      break ;
    }
  dragFast = -1 ;
}

static void traceLeftDown (int box, double xx, double yy) 
{
  float x1,x2,y1,y2 ;
  double x, y ;
  TRACELOOKGET("traceLeftDown") ;

  dragBox = box ;
  
  graphBoxDim (box, &x1, &y1, &x2, &y2) ;
  x = xx + x1 ;
  y = yy + y1 ;
  
  oldx = xOrigin = x ;
  oldy = yOrigin = y ;
  dragFast = 0 ;
  graphRegister (LEFT_DRAG, traceLeftDrag) ;	/* must redo */
  graphRegister (LEFT_UP, traceLeftUp) ;
}

/****************************************************/

LANE *tagLane = 0 ;

/****************************************************/

static void laneDoTag (LOOK look, LANE *lane, KEY key, char *buf)
{ int x1, x2 ;
  OBJ obj = 0 ;
  char *myTag = freekey2text (key, traceTagMenu) ;

  lane->hasTag = 0 ; /* force reevaluation */

  if (!lane->upSequence)
    x1 = look->map->centre - lane->x1 + lane->clipTop + 1 ;
  else
    x1 = - look->map->centre + lane->x1 + lane->clipTop ;	
  x1 -= 10 ; x2 = x1 + 20 ;
  if ((obj = bsUpdate (lane->key)))
    { bsAddData (obj, _Assembly_tags, _Text, myTag) ;
      bsAddData (obj, _bsRight, _Int, &x1) ;
      bsAddData (obj, _bsRight, _Int, &x2) ;
      if (buf)
	bsAddData (obj, _bsRight, _Text, buf) ;
      bsSave (obj) ;
    }
}

/****************************************************/

#ifdef APPARENTLY_NOT_NEEDED
static void traceTagger (KEY key)
{ LANE *lane ;
  int i ;
  char *cp ;
  TRACELOOKGET("traceTagger") ;
 
  i = arrayExists(look->laneShown) ? 
    arrayMax (look->laneShown) : 0 ;
  if (!i) return ; 
  if (!messPrompt ("Comments ?:", "", "t"))
    return ;

  cp = freeword () ;
  while (i--)
    { lane = arr (look->laneShown, i, LANE*) ;
      if (!lane->isAligned)
	laneDoTag (look, lane, key, cp) ;
    }
  traceDraw (look) ;
}
#endif /* APPARENTLY_NOT_NEEDED */


/****************************************************/
#ifdef OLDJUNK
static int tagBegin = -1 ;

static void traceDoAddTag (KEY key, int box)
{ HINT* hh ;
  MAP map ; LANE *lane ; 
  int min, max, tagEnd ;
  float mag, offSet ;
  OBJ obj ;
  TRACELOOKGET("traceDoAddTag") ;
  
  if (key == 1)
    return ;
  if (box < arrayMax(look->baseBoxes) &&
      (hh = arrp(look->baseBoxes,box, HINT)))
    { tagEnd = hh->iShort ;
      if (tagBegin > tagEnd)
	{ min = tagEnd ; tagEnd = tagBegin ; tagBegin = min ; } 
      if (tagBegin > 0 &&
	  tagLane == hh->lane &&
	  (obj = bsUpdate(hh->lane->key)))
	{ if (bsAddData (obj, key, _Int, &tagBegin))
	    bsAddData (obj, _bsRight, _Int, &tagEnd) ;
	  bsSave (obj) ;
	}
      /* reregister the menu */
      tagLane = 0 ;
      tagBegin = -1 ;
      lane = hh->lane ;

      graphBoxClear (lane->laneBaseCallBox) ;
      
      map = look->map ;
      map->centre -= lane->dy ;
      mag = map->mag ;
      map->mag *= lane->ddy ;
      lane->laneBaseCallBox = graphBoxStart () ;
      offSet = lane->laneBaseCallBoxOffSet ;
      min = lane->laneBaseCallBoxMin ;
      max = lane->laneBaseCallBoxMax ;
      
      traceDnaLane (look, look->map, lane, &offSet, min, max) ;
      graphBoxEnd () ;
      map->centre += lane->dy ;
      map->mag = mag ;
      
      graphBoxDraw (lane->laneBaseCallBox, -1, -1) ;
    }
}

/****************************************************/

static void traceAddTag (LOOK look, LANE *lane, int nn)
{ FREEOPT *fp ;
  MAP map = look->map ; int min, max ;
  float mag, offSet ;
  OBJ obj ; Array aa ;
  KEY model, tag ; int i, j ;
  
  tagBegin = nn ; tagLane = lane ;

#if !defined(NEW_MODELS)
 model =  KEYMAKE (_VSequence, 0) ;
#else
  model = pickList[_VSequence].model ;
#endif

  if (!traceAddTagMenu)
    { traceAddTagMenu = arrayCreate(20, FREEOPT) ;
      obj = bsCreate (model) ;
      j = 0 ;
      aa = arrayCreate(12, BSunit) ;
      if (obj)
	{ if (bsFindTag (obj, _Ace_mbly_tags) &&
	      bsFlatten (obj, 2, aa))
	    for (i = 1 ; i < arrayMax(aa) ; i+= 2)
	      { tag = arr(aa, i, BSunit).k ;
		fp = arrayp(traceAddTagMenu,++j,FREEOPT) ;
		fp->key = tag ; fp->text = name(tag) ;
	      } while (bsGetKeyTags (obj, _bsDown, &tag)) ;
	  bsDestroy (obj) ;
	}
      arrayDestroy (aa) ;
      fp = arrayp(traceAddTagMenu,++j,FREEOPT) ;
      fp->key = 1 ; fp->text = "Cancel" ;

      fp = arrayp(traceAddTagMenu,0,FREEOPT) ;
      fp->key = j ; fp->text = "toto" ;
      if (!j)
	{ arrayDestroy (traceAddTagMenu) ;
	  return ;
	}
    }
        /* reregister the menu */
  baseEditMenu = arrp(traceAddTagMenu, 0, FREEOPT) ;
  baseEdit = traceDoAddTag ;
  graphBoxClear (lane->laneBaseCallBox) ;

  map->centre -= lane->dy ;
  mag = map->mag ;
  map->mag *= lane->ddy ;
  lane->laneBaseCallBox = graphBoxStart () ;
  offSet = lane->laneBaseCallBoxOffSet ;
  min = lane->laneBaseCallBoxMin ;
  max = lane->laneBaseCallBoxMax ;

  traceDnaLane (look, look->map, lane, &offSet, min, max) ;
  graphBoxEnd () ;
  map->centre += lane->dy ;
  map->mag = mag ;

  graphBoxDraw (lane->laneBaseCallBox, -1, -1) ;
  
  baseEditMenu = look->mode == CDNA ? &baseDoEditcDNAMenu[0] : &baseDoEditMenu[0];
  baseEdit = baseEdit1 ;
}
  
/****************************************************/

static FREEOPT *traceMakeAddTagMenu (void)
{ FREEOPT *fp ;
  OBJ obj ; Array aa ;
  KEY model, tag ; int i, j ;
  
#if !defined(NEW_MODELS)
 model =  KEYMAKE (_VSequence, 0) ;
#else
  model = pickList[_VSequence].model ;
#endif

  if (!traceAddTagMenu)
    { traceAddTagMenu = arrayCreate(20, FREEOPT) ;
      obj = bsCreate (model) ;
      j = 0 ;
      aa = arrayCreate(12, BSunit) ;
      if (obj)
	{ if (bsFindTag (obj, _Ace_mbly_tags) &&
	      bsFlatten (obj, 2, aa))
	    for (i = 1 ; i < arrayMax(aa) ; i+= 2)
	      { tag = arr(aa, i, BSunit).k ;
		fp = arrayp(traceAddTagMenu,++j,FREEOPT) ;
		fp->key = tag ; fp->text = name(tag) ;
	      } while (bsGetKeyTags (obj, _bsDown, &tag)) ;
	  bsDestroy (obj) ;
	}
      arrayDestroy (aa) ;
      fp = arrayp(traceAddTagMenu,++j,FREEOPT) ;
      fp->key = 1 ; fp->text = "Cancel" ;

      fp = arrayp(traceAddTagMenu,0,FREEOPT) ;
      fp->key = j ; fp->text = "toto" ;
      if (!j)
	{ arrayDestroy (traceAddTagMenu) ;
	  return ;
	}
    }
        /* reregister the menu */
  return  arrp(traceAddTagMenu, 0, FREEOPT) ;
}
#endif 
/****************************************************/
/****************************************************/
/***************** segs ***********************************/

static void laneShowExtrema(LOOK look, LANE *lane, float *offset, 
			    int min, int max)
{ int i , k, k1 ;
  float y,
    x1 = *offset + look->traceWidth, 
    x2 = *offset,  x3 = (x1 + x2) / 2 ;
  MAP map = look->map ;
  int box ;
  BASECALL *xp ;
    
  CHECKLANE ; CHECKSEQ ; CHECKPOS ;
  if (lane->scf < 3)
    return ;

  if (! arrayExists(lane->baseCall) &&
      ! findBaseCall (lane))
    return ;
 
  if (! arrayMax(lane->baseCall))
    return ;

  for (k=0 ; k<4 ; k++)
    { box = graphBoxStart() ;
      
      xp = arrp(lane->baseCall, 0, BASECALL) ; 
      i = arrayMax(lane->baseCall) ;
      while(xp->x < min && i)
	xp++ , i-- ;
      while(i && xp->x < max)
	{ if (xp->t == k)
	    { y = MAP2GRAPH(map, xp->x + 1) ;
	      if (xp->flag & BC_LOW)
		graphLine(x3, y, x2, y) ;
	      else
		graphLine(x1, y, x2, y) ;
	    }
	  xp++ ; i-- ;
	}
      graphBoxEnd() ;
      k1 = lane->upSequence ? 3 - k : k ;
      
      graphBoxDraw(box, BASE_COLOR[k1], TRANSPARENT) ;
    }
  *offset += 0 ;  /* Redraw en place */
}

/***************************************************/

static void baseCallPatchButton (void)
{ int i ;
  LANE *lane ;
  TRACELOOKGET("patch") ;

  messStatus ("Autoedit") ;
  lane = arrp (look->lanes, 0, LANE) -1 ;
  i = arrayMax (look->lanes) ;
  while (lane++, i--)
    { if (!lane->isShown)
	continue ;
      CHECKLANE ;
      if (lane->scf < 2)
	continue ;
      if (findBaseCall(lane))
	{       
	  findXclipping (lane) ;
	  baseCallPatchLane (look->dna, look->dnaR, lane) ;
	  findXclipping (lane) ;
	  laneEditSave (look, lane, 0, 0) ;
	}
    }
  if (look->fMapLook)
    fMapPleaseRecompute (look->fMapLook) ;
  traceDraw (look) ;
}

/***************************************************/

static void traceRealignAll (void)
{ int i ;
  LANE *lane ;
  TRACELOOKGET("realign") ;

  messStatus ("Realigning") ;
  lane = arrp (look->lanes, 0, LANE) - 1 ;
  i = arrayMax (look->lanes) ;
  while (lane++, i--)
    { if (!lane->isShown && !lane->isAligned)
	continue ;
      CHECKLANE ;
      if (lane->scf < 2)
	continue ;
      lane->x2 += lane->upSequence ? +1 : -1 ;
/*       lane->clipEnd = arrayMax(lane->dna) ; */
      dnaAlignRecale(look->dna, &(lane->x1), &(lane->x2),
		     lane->dna, lane->clipTop, lane->clipEnd - 1) ;
      lane->x2 -= lane->upSequence ? +1 : -1 ;
      findXclipping (lane) ;
      laneEditSave (look, lane, 0, 0) ;
    }
  traceDraw (look) ;
}
#endif /* non graphic */

/***************************************************/
/***************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/
#ifndef NON_GRAPHIC

/*********************************************************************/
/* get all coords of transcribed genes on chrom */
static Array getTgHits (KEY chrom, int type)
{
  KEYSET tgs = 0 ;
  Array hits = 0 ;
  int i, a1, a2, j ;
  KEY tg, map ;
  OBJ Tg = 0 ;
  HIT *hh ;

  char *cp ;
  hits = arrayCreate (10000, HIT) ;
  if (type == 0)
    cp = messprintf ("Find Map \"%s\" ;> Transcribed_gene", name(chrom)) ;
  else if  (type == 1)
    cp = messprintf ("Find Map \"%s\" ;> Sequence", name(chrom)) ;
  else
    return hits ;
  tgs = queryKey (0, cp) ;

  for (i = 0, j= 0 ; i < arrayMax (tgs) ; i++)
    {
      tg = keySet (tgs, i) ;

      if ((Tg = bsCreate (tg)))
	{
	  if (bsGetKey (Tg, _IntMap, &map) &&
	      bsGetData (Tg, _bsRight, _Int, &a1) &&
	      bsGetData (Tg, _bsRight, _Int, &a2))
	    {
	      hh = arrayp (hits, j++, HIT) ;
	      hh->gene = tg ;
	      hh->a1 = a1 ;
	      hh->a2 = a2 ;
	    }
	  bsDestroy (Tg) ;
	}
    }
  
  keySetDestroy (tgs) ;

  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;

  return hits ;
}

/***********************************************************/

static Array getStops (Array dna, int frame)
{
  int max = arrayMax (dna) - 2, i , nn ;
  char *cp ;
  Array stops = arrayCreate (2000000, int) ;
  int j = 0 ;

  for (i = frame, cp = arrp (dna, i, char), nn = 0 ; i < max ; i += 3, cp += 3)
    {
      if (*cp == N_)
	nn++ ;
      if (*cp == T_ && 
	  (
	   (*(cp+1)== A_ && *(cp+2) == A_) ||
	   (*(cp+1)== A_ && *(cp+2) == G_) ||
	   (*(cp+1)== G_ && *(cp+2) == A_) 
	   )
	  )
	{
	  array (stops, j++, int) = i ;
	  array (stops, j++, int) = nn ;
	  nn = 0 ;
	}
    }  
  return stops ;
}

/***********************************************************/

static void getStopHisto (KEY chrom, Array histo, Array histo2, Array stops, Array hits, int type)
{
  int i, *s, s0, dx, *h, j, a1, a2, nn ;
  HIT *hh ;

  a1 = a2 = j = -1 ;

  if (histo && histo2 && hits && stops && arrayMax(stops))
    for (i = 0, j = 0, s0 = 0, s = arrp (stops, 0, int);
	 i < arrayMax(stops) ; i +=2 , s += 2)
      { 
	dx = (*s - s0 - 3)/3 ; 
	nn = *(s+1) ;
	if (type == 0 && nn < 5 && dx > 0)
	  {
	    h = arrayp (histo, dx, int) ;
	    (*h)++ ;
	    if (dx > 20000)
	      printf ("%d\t= %d\t- %d   %s\n", dx , *s, s0, name(chrom)) ;
	  }
	while (a2 < s0 && j < arrayMax(hits) - 1)
	  { 
	    hh = arrp(hits, ++j, HIT) ;
	    a1 = hh->a1 ; a2 = hh->a2 ;
	  }
	if (dx > 0 && nn < 5 && a2 > s0 && a1 < *s)
	  {
	    h = arrayp (histo2, dx, int) ;
	    (*h)++ ;
	  }
	s0 = *s ;
      }
}

/***********************************************************/

static int chromHisto (KEY chrom, Array histo, Array histo2, Array histo3)
{
  Array dna = 0, stops = 0, tgs = 0, pgs = 0 ;
  int nn = 0, frame ;

  dna = dnaGet (chrom) ;
  if (dna && arrayMax(dna))
    {
      tgs = getTgHits (chrom, 0) ;
      pgs = getTgHits (chrom, 1) ;
      for (frame = 0 ; frame < 3 ; frame++)
	{
	  stops = getStops (dna, frame) ;
	  getStopHisto (chrom, histo, histo2, stops, tgs, 0) ;
	  getStopHisto (chrom, histo, histo3, stops, pgs, 1) ;
	  nn += arrayMax(stops)/2 ;
	  arrayDestroy (stops) ;
	}
      reverseComplement (dna) ;
      for (frame = 6 ; frame < 3 ; frame++)
	{
	  stops = getStops (dna, frame) ;
	  getStopHisto (chrom, histo, histo2, stops, tgs, 0) ;
	  getStopHisto (chrom, histo, histo3, stops, pgs, 1) ;
	  nn += arrayMax(stops)/2 ;
	  arrayDestroy (stops) ;
	}
    }
  arrayDestroy (tgs) ;
  arrayDestroy (pgs) ;
  arrayDestroy (dna) ;
  return nn ;
}

/***********************************************************/

void cdnaStopHisto (void)
{
  KEYSET ks = query (0, "Find sequence genomic ; !Overlap_left") ;
  KEY chrom, cosmid ;
  int i, nn, ntot = 0 ;
  Array histo = arrayCreate (30000, int) ;
  Array histo2 = arrayCreate (30000, int) ;
  Array histo3 = arrayCreate (30000, int) ;

  for (i = 0 ; i < keySetMax (ks) ; i++)
    {
      cosmid = keySet (ks, i) ;
      while (( chrom = keyGetKey (cosmid, _Source)))
	{ cosmid = chrom ; }
      nn = chromHisto (cosmid, histo, histo2, histo3) ;
      ntot += nn ;
    }
  for (i = 0 ; i < arrayMax (histo2) ; i++)
    {
      int 
	x = 100 * array (histo2, i, int),
	y = array (histo, i, int) ;
      if (y > 0)
	array (histo2, i, int) = x/y ;
    }
  for (i = 0 ; i < arrayMax (histo3) ; i++)
    {
      int 
	x = 100 * array (histo3, i, int),
	y = array (histo, i, int) ;
      if (y > 0)
	array (histo3, i, int) = x/y ;
    }
  if (histo)
    plotHisto ("ORF lengths", histo) ; /* will destroy histo */
  if (histo2)
    plotHisto ("percent in transcribed genes/ORF length", histo2) ; /* will destroy histo */
  if (histo3)
    plotHisto ("percent in predicted genes/ORF length", histo3) ; /* will destroy histo */
  keySetDestroy (ks) ;
  printf("total i= %d stop = %d\n", i, ntot) ;
}

#endif /* non graphic */

/***********************************************************/
/***********************************************************/
/***********************************************************/
