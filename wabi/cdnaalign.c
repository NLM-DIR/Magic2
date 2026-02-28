/*  File: cdnaalign.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg. 1997
 *-------------------------------------------------------------------
 * This file is part of the ACEMBLY extension to 
        the ACEDB genome database package, written by
 *         Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *        Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
   Align cDNA
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  9 15:18 1998 (fw)
 * Created: Thu Dec  9 00:01:50 1997 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */
#define ARRAY_CHECK
#define MALLOC_CHECK


#include "acembly.h"
#include "freeout.h"
#include "vtxt.h"
#include "bs.h"
#include "a.h"
#include "systags.h"
#include "tags.h"
#include "classes.h"
#include "dna.h"
#include "cdna.h"
#include "dnaalign.h"
#include "bump.h"
#include "query.h"
#include "session.h"
#include "bindex.h"
#include "lex.h"
#include "mytime.h"
#include "parse.h"
#include "pick.h"
#include "dump.h"
#include "command.h"
#include "cdnapath.h"
#include "cdnatr.h"
#include "makemrna.h" 
#include "chrono.h"
#include "percolate.h"
#include "basecall.h"


#define SWAP(_u1,_u2,_type) { _type _dummy = (_u1) ; _u1 = (_u2) ; _u2 = _dummy ; }
static int OLIGO_LENGTH = 15 ;
static int OLIGO_LENGTH_DIFFICULT = 12 ;
int acemblyMode = SHOTGUN ;

static void showSpl (Array allSpl) ;

static BOOL getBestUpPosition (Array cDNA, Array longDna, int x1, int *x2p, int a1, int *a2p, BOOL ignoreN, BOOL slide) ;
static BOOL getBestVpPosition (Array cDNA, Array longDna, int *y1p, int y2, int *b1p, BOOL ignoreN, BOOL slide) ;
static void showOligo (unsigned int oligo) ;
static void countErrors (KEY cosmid, Array hits, Array dna, Array dnaR, int maxX, BOOL countN) ;

static Associator estAssCreate (Array words,
                                KEYSET estSet, BOOL reverse, 
                                KEYSET cDNA_clones, KEYSET clipTops, KEYSET clipEnds, KEYSET trackingClones,
                                KEYSET estMaps, KEYSET estMap1, KEYSET estMap2, 
                                int *n1p, int *n2p, int type, int oligoLength, BOOL isDifficult, 
				KEY chrom, int c1, int c2
				) ;
static void getCosmidHits (KEY cosmid, Array hits, 
                           Associator ass, Associator assR, 
                           Array dna, BOOL isUp, KEYSET clipEnds, int nn, 
			   KEYSET estMaps, KEYSET estMap1, KEYSET estMap2, 
			   KEY cosmidMap, int cosmidC1, int cosmidC2,
			   int direction) ;

int ficheAceRunKeySet(KEYSET ks,void * funccall,void * callback,int params);
Array getSeqDna(KEY key) ;
static Array teamCreateIntrons (Array hits, Array dna, Array dnaR) ;
static void cDNARealignGeneKeySet (KEYSET ks, BOOL doFuse, int locally, int searchRepeats, int splitCloud) ;

/*********************************************************************/
/*********************************************************************/
/************** Strategy, imported from clone MainClone **************/

/*********************************************************************/
/* set == 0 return current value 
   set == -1 set to unknown (will reask from the data)
   set == -2 set to false
   set == 1 set to true
*/
   
int getPleaseDoRepeats (int set)
{
  static int style = -1 ;

  switch (set)
    {
    case 0: break ;               /* return current value */
    case 1: style = 1 ; break ;   /* set to true  */
    case -2: style = 0 ; break ;  /* set to false  */
    case -1: style = -1 ; break ; /* set to unknown */
      
    }
  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && SearchForTandemRepeats") ;
      style = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    style = 0 ; /* default */

  return style ;
}

/*********************************************************************/

static int getPleaseMinIntronSize (void)
{
  static int minIntronSize = -1 ;

  if (minIntronSize == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && MinIntronSize") ;
      minIntronSize = keySetMax (ks) ? 1 : -1 ;
      
      if (minIntronSize == 1)
        {
          OBJ Clone = bsCreate (keySet (ks,0)) ;
          
          if (Clone)
            {
              bsGetData (Clone, str2tag("MinIntronSize"), _Int, &minIntronSize) ;
              bsDestroy (Clone) ;
            }
        }
      keySetDestroy (ks) ;
    }
  if (minIntronSize == -1)
    minIntronSize = 0 ; /* default */

  return minIntronSize ;
} /* getPleaseMinIntronSize */

/*********************************************************************/

static int getIgnoreWalls (void)
{
  static int ignoreWalls = -1 ;
  
  if (ignoreWalls == -1)
    {
      KEYSET ks = query (0, "Find clone  main_Clone && IgnoreWalls") ;
      ignoreWalls = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  return ignoreWalls ;
}

/*********************************************************************/

static int getIgnoreDiscard (void)
{
  static int ignoreDiscard = -1 ;
  
  if (ignoreDiscard == -1)
    {
      KEYSET ks = query (0, "Find clone  main_Clone && IgnoreDiscard") ;
      ignoreDiscard = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  return ignoreDiscard ;
}

/*********************************************************************/

static int getAddPolyA (void)
{
  static int addPolyA = -1 ;
  
  if (addPolyA == -1)
    {
      KEYSET ks = query (0, "Find clone  main_Clone && AddPolyA") ;
      addPolyA = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  return addPolyA ;
}

/*********************************************************************/

static char* getMaskFrequentWords (void)
{
  static char *wName = (char *) 0x1 ;
  
  if (wName == (char *) 0x1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && MaskFrequentWords") ;
      KEY clone = keySetMax (ks) ? keySet(ks, 0) : 0 ;
      OBJ Clone = clone ? bsCreate (clone) : 0 ;
      char *cp = 0 ;

      if (Clone)
        {
          bsGetData (Clone, str2tag("MaskFrequentWords"), _Text, &cp) ;
          bsDestroy (Clone) ;
        }
      keySetDestroy (ks) ;

      if (cp)
        wName = strnew (cp, 0) ;
      else
        wName = 0 ;
    }
  return wName ;
}

/*********************************************************************/

static int  getMaxIntronSize (void)
{
  static int limit = -1 ;
  
  if (limit == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && MaxIntronSize") ;
      KEY key = keySetMax(ks) ? keySet (ks, 0) : 0 ;
      OBJ Clone = 0 ;
      
      if (key && (Clone = bsCreate(key)))
        {
          bsGetData (Clone, str2tag("MaxIntronSize"), _Int, &limit) ;
          bsDestroy (Clone) ;
        }
      if (limit < 0) 
        limit = 0 ;
      keySetDestroy (ks) ;
    }
  return limit ;
} /* getMaxIntronSize */

static void errRemoveAmbigue (Array err)
{
  A_ERR *ep, *ep1 ;
  int i, j ;
  
  if (err && arrayMax(err))
    {
      for (i = j = 0, ep = ep1 = arrp (err, 0, A_ERR) ; i < arrayMax(err) ; ep++, i++)
        {
          if (ep->type == AMBIGUE)
            continue ;
          if (i > j)
            *ep1 = *ep ;
          ep1++ ; j++ ;
        }
      arrayMax (err) = j ;
    }
}

/*********************************************************************/
/* the idea here is to clip the hanging bp if they look
 * remotely like a vector or a very low quality stuff
 *
 * err is com,puted from x1 to *x2p
 * accepted ali before is from x1 to oldx2
 * accepted ali after is from x1 to *x2p
*/
static int vectorReclip (Array shortDna, int x1, int oldx2, int *x2p, 
                          Array longDna, int a1, int olda2, int *a2p, 
                          Array err) 
{
  A_ERR *ep, *ep1 ;
  int i, jj, n, dx, myjj, myx2 ;
  char *cp, *cq ;
   
  if (err && arrayMax(err))
    {
      myx2 = *x2p /* oldx2 */  ; myjj = arrayMax (err) ; ep1 = 0 ;
      for (jj = arrayMax (err) - 1, ep = arrp (err, jj, A_ERR) ; jj >= 0 ; ep--, jj--)
        {
          if (myx2 - ep->iShort < 3) /* this error does not bring 3 bp */
            { myx2 = ep->iShort ; ep1 = ep ; myjj = jj ; }
          else
            break ;
        }
      if (ep1)  /* AL832583 */
        { 
          dx = *x2p - ep1->iShort ;
          *a2p -= dx ;
          *x2p -= dx ;
          return myjj ; /* cut before the stupid errors */
        }

      
      for (jj = 0, ep = arrp (err, 0, A_ERR) ; jj < arrayMax(err) ; ep++, jj++)
        { 
          i = oldx2 - ep->iShort ;
          
          if (i > 2 && i < 30) /* check for suspect sequencer tuning or poly A and bad orientation */
            {
              if (i > 18) i = 18 ; /* strat 4 bp upstream because of TTTTT leads */
              for (n = - i/2, cp = arrp (shortDna, ep->iShort, char) - 4, cq = cp + 1 ; 
                   i > 0 ; cp++, cq++, i--)
                if (*cp & *cq) n++ ;
              if (n > 0) /* 50% repeats */
                { 
                  dx = *x2p - ep->iShort ;
                  *a2p -= dx ;
                  *x2p -= dx ;
                  return jj ; /* cut before the stupid repetition */
                }            
            }
        }
      return arrayMax(err) ;
    }
  return 0 ;
}

/*********************************************************************/
/* makeMrnaGene  only called from makeMrnaGene , so gmrna->estHits is always positive strand a1 < a2 */
BOOL cdnaReExtendHits (S2M *s2m, SC* sc, SMRNA *gmrna, KEYSET clipTops, KEYSET clipEnds)
{
  Array dnaD = 0, dnaR = 0, shortDna, hits = gmrna->estHits ;
  Array err = 0, longDna, longD, longR ;
  HIT *up, *vp, *vpMin, *vpMax ;
  int ii, ga1, a1, a2, olda2, x1, x2, oldx2, nerr, Nn ;
  KEY oldEst = 0 ;
  int maxJump = 4 ;
  BOOL cross, ok = FALSE, isBest ;
  KEY clone ;

  if (sc->a1 < sc->a2)
    { 
      longD = s2m->dnaD ;  
      longR = s2m->dnaR ;  
      ga1 = sc->a1 ;
      /* ga2 = sc->a2 ; */
    }
  else
    { 
      longD = s2m->dnaR ;  
      longR = s2m->dnaD ; 
      ga1 = (long)arrayMax (longD) - sc->a1 + 1 ;
      /*  ga2 = (long)arrayMax (longD) - sc->a2 + 1 ; */
    }
   
  vpMin = arrp(hits, 0, HIT) ;
  vpMax = arrp(hits, arrayMax(hits) - 1, HIT) ;
  for (ii = 0, up = arrp(hits, 0, HIT) ; ii < arrayMax(hits) ; ii++, up++)
    {
      clone = up->cDNA_clone ;
      if (up->est != oldEst)
        {
          dnaD = 0 ;
          arrayDestroy (dnaR) ;
          oldEst = up->est ;
        }

      if (up->type & gLink) /* 3p 5poverlapping pair of reads  */
        continue ;

      if ((up->type & gX) && !(up->type & gD) &&  (ii == 0 || up->est != (up-1)->est))
        {
          if (up->x1 > up->clipTop && up->x1 < up->x2)  /* case 1: begin of down 5' read */
            { 
              if (!dnaD)
                dnaD = getSeqDna (up->est) ;
              if (!dnaD)
                continue ;

              a1 = ga1 + up->a1 - 1 ; x1 = up->x1 ; a2 = ga1 + up->a2 - 1 ; x2 = up->x2 ;
              if ((isBest = getBestVpPosition (dnaD, longD, &x1, x2, &a1, FALSE, FALSE)))
                {
                  x2 = x1 + 1 ; a2 = a1 + 1 ; /* nice matching pos */
                }
              if (isBest && x1 <= up->clipTop) /* adjust directly */
                {
                  int dx = up->clipTop - up->x1 ;

                  up->x1 += dx ;
                  up->a1 += dx ;
                }
              else
                {
                  longDna = longR ;
                  a1 = (long)arrayMax(longDna) - a2 + 1 ; /* nice match */
                  a2 = olda2 = (long)arrayMax(longDna) - (ga1 + up->a1 - 1) + 1 ;
                  
                  if (!dnaR) { dnaR = dnaCopy(dnaD) ; reverseComplement (dnaR) ; }
                  shortDna = dnaR ;
                  x1 = (long)arrayMax (dnaR) - x2 + 1 ; /* nice match */
                  oldx2 = (long)arrayMax (dnaR) - up->x1 + 1 ;
                  x2 = (long)arrayMax (dnaR) - keySet(clipTops, KEYKEY(up->est)) + 1 ;
                  a2 += 12 ;
                  a1-- ; x1-- ; /* we start at the end of the exon and move back up */
                  err = aceDnaTrackErrors(shortDna, x1, &x2, longDna, a1, &a2, 0, err, maxJump, 0, FALSE, 0, FALSE) ;
                  nerr = vectorReclip (shortDna, x1, oldx2, &x2, longDna, a1, olda2, &a2, err) ;
                  if (x2 != oldx2 && a2 > olda2 - 12 && a2 > a1 + 1)         
                    { 
		      up->x1 -= (x2 - oldx2) ; up->a1 -= (a2 - olda2) ;
		      x1 = up->x1 - 1 ; oldx2 = x2 = up->x2 ;  Nn = 0 ;
		      a1 =  (long)arrayMax(longDna) - up->a2 + 1 - 1 ; 
		      a2 = olda2 =  (long)arrayMax(longDna) - up->a1 + 1 ;
		      err = aceDnaTrackErrors(shortDna, x1, &x2, longDna, a1, &a2, &Nn, err, maxJump, 0, FALSE, 0, FALSE) ;
		      up->nerrAll = (long)arrayMax (err) - Nn ;
		      ok = TRUE ; 
		    }
                }
            }
          
          if (up->x1 > up->x2 && /* case 3: end of up 3' read */
              !(up->type & gLink)) /* 3p 5poverlapping pair of reads  */
            {  
              if (!dnaD)
                dnaD = getSeqDna (up->est) ;
              if (!dnaD)
                continue ;

              /*sans doute inutile, a cause de gLink */
               cross = FALSE ; 
              vp = up ;
              while (!cross && --vp >= vpMin && vp->cDNA_clone == clone)
                if (vp->est != oldEst && vp->x1 < vp->x2 && vp->a1 < up->a1 && vp->a2 > up->a1)
                  cross = FALSE ;
              vp = up ;
              while (!cross && ++vp <= vpMax && vp->cDNA_clone == clone)
                if (vp->est != oldEst && vp->x1 < vp->x2 && vp->a1 < up->a1 && vp->a2 > up->a1)
                  cross = FALSE ;
              if (cross)
                continue ;
             
              longDna = longR ;
              a1 = (long)arrayMax(longDna) - (ga1 + up->a2 - 1) + 1 ;
              a2 = olda2 = (long)arrayMax(longDna) - (ga1 + up->a1 - 1) + 1 ;

              shortDna = dnaD ;
              x1 = up->x2 ;
              x2 = up->x1 ;
                        
              if ((isBest = getBestUpPosition (dnaD, longDna, x1, &x2, a1, &a2, FALSE, FALSE)))
                {
                  x1 = x2 - 1 ; a1 = a2 - 1 ; /* nice matching pos */
                }
              if (isBest && x2 >= up->clipEnd) /* adjust directly */
                {
                  int dx = up->clipEnd - up->x1 ;

                  up->x1 += dx ;
                  up->a1 -= dx ;
                }
              else
                {
                  a2 = olda2 + 12 ;
                  x2 = keySet(clipEnds, KEYKEY (up->est)) ; oldx2= up->x1 ;
                  a1-- ; x1-- ;
                  err = aceDnaTrackErrors(shortDna, x1, &x2, longDna, a1, &a2, 0, err, maxJump, 0, FALSE, 0, FALSE) ;
                  nerr = vectorReclip (shortDna, x1, oldx2, &x2, longDna, a1, olda2, &a2, err) ; 
                  if (x2 != oldx2 && a2 > olda2 - 12 && a2 > a1 + 1)         
                    { 
		      up->x1 += (x2 - oldx2) ; up->a1 -= (a2 - olda2) ; up->nerrAll = nerr ;
		      x1 = up->x2 - 1 ; oldx2 = x2 = up->x1 ;  Nn = 0 ;
		      a1 = (long)arrayMax(longDna) - (ga1 + up->a2 - 1) + 1 - 1 ;
		      a2 = olda2 = (long)arrayMax(longDna) - (ga1 + up->a1 - 1) + 1 ;
		      err = aceDnaTrackErrors(shortDna, x1, &x2, longDna, a1, &a2, &Nn, err, maxJump, 0, FALSE, 0, FALSE) ;
		      up->nerrAll = (long)arrayMax (err) - Nn ;
		      ok = TRUE ; 
		    }
                }
            }
        } 
 
      if (!(up->type & (gF | gA)) &&  (ii ==  arrayMax(hits) - 1 || up->est != (up+1)->est))
        {
          if (up->x1 < up->x2 && /* case 2: end of down 5' read */
              !(up->type & gLink)) /* 3p 5poverlapping pair of reads  */
            { 
              /*sans doute inutile, a cause de gLink */
              cross = FALSE ; 
              vp = up ;
              while (!cross && --vp >= vpMin && vp->cDNA_clone == clone)
                if (vp->est != oldEst && vp->x1 > vp->x2 && vp->a1 < up->a2 && vp->a2 > up->a2)
                  cross = TRUE ;
              vp = up ;
              while (!cross && ++vp <= vpMax && vp->cDNA_clone == clone)
                if (vp->est != oldEst && vp->x1 > vp->x2 && vp->a1 < up->a2 && vp->a2 > up->a2)
                  cross = TRUE ;
              if (cross)
                continue ;

              longDna = longD ;
              a1 = up->a1 + ga1 - 1 ;
              a2 = olda2 = up->a2 + ga1 - 1 ;
              
              if (!dnaD)
                dnaD = getSeqDna (up->est) ;
              if (!dnaD)
                continue ;

              shortDna = dnaD ;
              x1 = up->x1 ;
              x2 = up->x2 ;
              if ((isBest = getBestUpPosition (dnaD, longDna, x1, &x2, a1, &a2, FALSE, FALSE)))
                {
                  x1 = x2 - 1 ; a1 = a2 - 1 ; /* nice matching pos */
                }
              if (isBest && x2 >= up->clipEnd) /* adjust directly */
                {
                  int dx = up->clipEnd - up->x2 ;

                  up->x2 += dx ;
                  up->a2 += dx ;
                }
              else
                {
                  a2 = olda2 + 12 ;
                  x2 = keySet(clipEnds, KEYKEY (up->est)) ; oldx2= up->x2 ;
                  
                  a1-- ; x1-- ;
                  err = aceDnaTrackErrors(shortDna, x1, &x2, longDna, a1, &a2, 0, err, maxJump, 0, FALSE, 0, FALSE) ;
                  nerr = vectorReclip (shortDna, x1, oldx2, &x2, longDna, a1, olda2, &a2, err) ;
                  if (x2 != oldx2 && a2 > olda2 - 12 && a2 > a1 + 1)                   
                    {
		      up->x2 += (x2 - oldx2) ; up->a2 += (a2 - olda2) ;  Nn = 0 ;
		      x1 = up->x1 - 1 ; oldx2 = x2 = up->x2 ; a1 = up->a1 + ga1 -1 ; a2 = olda2 = up->a2 + ga1 - 1 ;
		      err = aceDnaTrackErrors(shortDna, x1, &x2, longDna, a1, &a2, &Nn, err, maxJump, 0, FALSE, 0, FALSE) ;
		      up->nerrAll = (long)arrayMax (err) -Nn ;
		      ok = TRUE ;
		    }
                }
            }
          
          if (up->x1 > up->x2)  /* case 4: begin of up 3' read */
            {   
              longDna = longD ;
              a1 = up->a1 + ga1 - 1 ; x1 = up->x1 ; 
              a2 = olda2 = up->a2 + ga1 - 1 ; x2 = up->x2 ;

              if (!dnaD)
                dnaD = getSeqDna (up->est) ;
              if (!dnaD)
                continue ;

              if (!dnaR) { dnaR = dnaCopy(dnaD) ; reverseComplement (dnaR) ; }
              shortDna = dnaR ;
              x1 = (long)arrayMax (dnaR) - up->x1 + 1 ;
              x2 = oldx2 = (long)arrayMax (dnaR) - up->x2 + 1 ;

              if ((isBest = getBestUpPosition (dnaR, longDna, x1, &x2, a1, &a2, FALSE, FALSE)))
                {
                  x1 = x2 - 1 ; a1 = a2 - 1 ; /* nice matching pos */
                }
              if (isBest && (long)arrayMax (dnaR) - x2 + 1 <= up->clipTop) /* adjust directly */
                {
                  int dx = up->clipTop - up->x2 ;

                  up->x2 += dx ;
                  up->a2 -= dx ;
                }
              else
                {
                  a2 = olda2 + 12 ;
                  x2 = (long)arrayMax (dnaR) - keySet(clipTops, KEYKEY(up->est)) + 1 ;
                  
                  a1-- ; x1-- ;
                  err = aceDnaTrackErrors(shortDna, x1, &x2, longDna, a1, &a2, 0, err, maxJump, 0, FALSE, 0, FALSE) ;
                  nerr = vectorReclip (shortDna, x1, oldx2, &x2, longDna, a1, olda2, &a2, err) ; 
                  if (x2 != oldx2 && a2 > olda2 - 12 && a2 > a1 + 1)                   
                    { 
		      up->x2 -= (x2 - oldx2) ; up->a2 += (a2 - olda2) ; 
		      x1 =  (long)arrayMax (dnaR) - up->x2 - 1 ; 
		      oldx2 = x2 =  (long)arrayMax (dnaR) - up->x1 ; 
		      a1 = up->a1 + ga1 - 2 ; a2 = olda2 = up->a2 + ga1 - 1 ; Nn = 0 ;
		      err = aceDnaTrackErrors(shortDna, x1, &x2, longDna, a1, &a2, &Nn, err, maxJump, 0, FALSE, 0, FALSE) ;
		      up->nerrAll = (long)arrayMax (err) - Nn ;
		      ok = TRUE ; }
                }
            }
        }
      oldEst = up->est ;
    } 

  arrayDestroy (dnaR) ;
  arrayDestroy (err) ;

  return ok ;
}  /* cdnaReExtendHits */

/*********************************************************************/
/* makeMrnaGene  only called from makeMrnaGene , so gmrna->estHits is always positive strand a1 < a2 */
void JUNKOLDcdnaReExtendHits (S2M *s2m, SC* sc, SMRNA *gmrna, KEYSET clipTops, KEYSET clipEnds)
{
  /* dec11 ???? this code is absolute junk in the way it seems to use dnad dnar */
  Array dnaD = 0, dnaR = 0, shortDna, hits = gmrna->estHits ;
  Array err = 0, longDna = 0 ;
  HIT *up, *vp, *vpMin, *vpMax ;
  int ii, ga1 = sc->a1, ga2 = sc->a2, a1, a2, olda2, x1, x2, oldx2 ;
  KEY oldEst = 0 ;
  int maxJump = 4 ;
  BOOL cross ;
  KEY clone ;

  vpMin = arrp(hits, 0, HIT) ;
  vpMax = arrp(hits, arrayMax(hits) - 1, HIT) ;
  for (ii = 0, up = arrp(hits, 0, HIT) ; ii < arrayMax(hits) ; ii++, up++)
    {
      clone = up->cDNA_clone ;

      if (up->est != oldEst && !(up->type & gD))
        {
          dnaD = getSeqDna (up->est) ;
          arrayDestroy (dnaR) ;
          if (!dnaD)
            continue ;

          if (ga1 < ga2 && up->x1 < up->x2)  /* case 1: begin of down 5' read */
            { 
              longDna = s2m->dnaR ;
              a1 = (long)arrayMax(longDna) - (ga1 + up->a2 - 1) + 1 ;
              a2 = olda2 = (long)arrayMax(longDna) - (ga1 + up->a1 - 1) + 1 ;
              
              if (!dnaR) { dnaR = dnaCopy(dnaD) ; reverseComplement (dnaR) ; }
              shortDna = dnaR ;
              x1 = (long)arrayMax (dnaR) - up->x2 + 1 ;
              oldx2 = (long)arrayMax (dnaR) - up->x1 + 1 ;
              x2 = (long)arrayMax (dnaR) - keySet(clipTops, KEYKEY(up->est)) + 1 ;
             
              a1-- ; x1-- ; /* we start at the end of the exon and move back up */
              err = aceDnaTrackErrors(shortDna, x1, &x2, longDna, a1, &a2, 0, err, maxJump, 0, TRUE, 0, FALSE) ;
              
              if ((x2 > oldx2 && a2 > olda2) ||
                  vectorReclip (shortDna, x1, oldx2, &x2, longDna, a1, olda2, &a2, err))
                { up->x1 -= (x2 - oldx2) ; up->a1 -= (a2 - olda2) ; }
            }
          
          if (ga1 < ga2 && up->x1 > up->x2) /* case 3: end of up 3' read */
            {  
              cross = FALSE ; 
              vp = up ;
              while (!cross && --vp >= vpMin && vp->cDNA_clone == clone)
                if (vp->est != oldEst && vp->x1 < vp->x2 && vp->a1 < up->a1 && vp->a2 > up->a1)
                  cross = FALSE ;
              vp = up ;
              while (!cross && ++vp <= vpMax && vp->cDNA_clone == clone)
                if (vp->est != oldEst && vp->x1 < vp->x2 && vp->a1 < up->a1 && vp->a2 > up->a1)
                  cross = FALSE ;
              if (cross)
                continue ;

              longDna = s2m->dnaR ; 
              a1 = (long)arrayMax(longDna) - (ga1 + up->a2 - 1) + 1 ;
              a2 = olda2 = (long)arrayMax(longDna) - (ga1 + up->a1 - 1) + 1 ;
                     
              shortDna = dnaD ;
              x1 = up->x1 ;
              x2 = keySet(clipEnds, KEYKEY (up->est)) ; oldx2= up->x2 ;
              
              a1-- ; x1-- ;
              err = aceDnaTrackErrors(shortDna, x1, &x2, longDna, a1, &a2, 0, err, maxJump, 0, TRUE, 0, FALSE) ;
              
              if (x2 > up->x2 && a2 > olda2)
                { up->x2 += (x2 - oldx2) ; up->a1 -= a2 - olda2 ; }
            }

          if (ga1 > ga2 && up->x1 < up->x2)  /* case 5: begin of up 5' read */
            { 
              longDna = s2m->dnaD ;
              a1 = ga1 - up->a2 + 1 ;
              a2 = olda2 = ga1 - up->a1 + 1 ;
              
              if (!dnaR) { dnaR = dnaCopy(dnaD) ; reverseComplement (dnaR) ; }
              shortDna = dnaR ;
              x1 = (long)arrayMax (dnaR) - up->x2 + 1 ;
              oldx2 = (long)arrayMax (dnaR) - up->x1 + 1 ;
              x2 = (long)arrayMax (dnaR) - keySet(clipTops, KEYKEY(up->est)) + 1 ;

              a1-- ; x1-- ; 
              err = aceDnaTrackErrors(shortDna, x1, &x2, longDna, a1, &a2, 0, err, maxJump, 0, TRUE, 0, FALSE) ;
              
              if ((x2 > oldx2 && a2 > olda2) ||
                  vectorReclip (shortDna, x1, oldx2, &x2, longDna, a1, olda2, &a2, err))
                { up->x1 -= (x2 - oldx2) ; up->a1 -= (a2 - olda2) ; }
            }
          
          if (ga1 > ga2 && up->x1 > up->x2) /* case 7: end of down 3' read */
            {   
              cross = FALSE ; 
              vp = up ;
              while (!cross && --vp >= vpMin && vp->cDNA_clone == clone)
                if (vp->est != oldEst && vp->x1 < vp->x2 && vp->a1 < up->a1 && vp->a2 > up->a1)
                  cross = FALSE ;
              vp = up ;
              while (!cross && ++vp <= vpMax && vp->cDNA_clone == clone)
                if (vp->est != oldEst && vp->x1 < vp->x2 && vp->a1 < up->a1 && vp->a2 > up->a1)
                  cross = FALSE ;
              if (cross)
                continue ;

              longDna = s2m->dnaD ;
              a1 = ga1 - up->a2 + 1 ;
              a2 = olda2 = ga1 - up->a1 + 1 ;
              
              shortDna = dnaD ;
              x1 = up->x1 ;
              x2 = keySet(clipEnds, KEYKEY (up->est)) ; oldx2= up->x2 ;
              
              a1-- ; x1-- ;
              err = aceDnaTrackErrors(shortDna, x1, &x2, longDna, a1, &a2, 0, err, maxJump, 0, TRUE, 0, FALSE) ;
              
              if (x2 > up->x2 && a2 > olda2)
                { up->x2 += (x2 - oldx2) ; up->a1 -= (a2 - olda2) ; }
            }
          
        } 
 
      if (dnaD && !(up->type & gF) &&  (ii ==  arrayMax(hits) - 1 || up->est != (up+1)->est))
        {
          if (ga1 < ga2 && up->x1 < up->x2)  /* case 2: end of down 5' read */
            { 
              cross = FALSE ; 
              vp = up ;
              while (!cross && --vp >= vpMin && vp->cDNA_clone == clone)
                if (vp->est != oldEst && vp->x1 > vp->x2 && vp->a1 < up->a2 && vp->a2 > up->a2)
                  cross = TRUE ;
              vp = up ;
              while (!cross && ++vp <= vpMax && vp->cDNA_clone == clone)
                if (vp->est != oldEst && vp->x1 > vp->x2 && vp->a1 < up->a2 && vp->a2 > up->a2)
                  cross = TRUE ;
              if (cross)
                continue ;

              longDna = s2m->dnaD ;
              a1 = up->a1 + ga1 - 1 ;
              a2 = olda2 = up->a2 + ga1 - 1 ;
              
              shortDna = dnaD ;
              x1 = up->x1 ;
              x2 = keySet(clipEnds, KEYKEY (up->est)) ; oldx2= up->x2 ;

              a1-- ; x1-- ;
              err = aceDnaTrackErrors(shortDna, x1, &x2, longDna, a1, &a2, 0, err, maxJump, 0, TRUE, 0, FALSE) ;
              
              if (0 && x2 > oldx2 && a2 > olda2)
                { up->x2 += (x2 - oldx2) ; up->a2 += (a2 - olda2) ; }

              if ((x2 > oldx2 && a2 > olda2) ||
                  vectorReclip (shortDna, x1, oldx2, &x2, longDna, a1, olda2, &a2, err))
                { up->x2 += (x2 - oldx2) ; up->a2 += (a2 - olda2) ; }
            }
          
          if (ga1 < ga2 && up->x1 > up->x2)  /* case 4: begin of up 3' read */
            {   
              longDna = s2m->dnaD ; 
              a1 = up->a1 + ga1 - 1 ;
              a2 = olda2 = up->a2 + ga1 - 1 ;

              if (!dnaR) { dnaR = dnaCopy(dnaD) ; reverseComplement (dnaR) ; }
              shortDna = dnaR ;
              x1 = (long)arrayMax (dnaR) - up->x2 + 1 ;
              oldx2 = (long)arrayMax (dnaR) - up->x1 + 1 ;
              x2 = (long)arrayMax (dnaR) - keySet(clipTops, KEYKEY(up->est)) + 1 ;

              a1-- ; x1-- ;
              err = aceDnaTrackErrors(shortDna, x1, &x2, longDna, a1, &a2, 0, err, maxJump, 0, TRUE, 0, FALSE) ;
              
              if (x2 > oldx2 && a2 > olda2)
                { up->x1 -= (x2 - oldx2) ; up->a2 += (a2 - olda2) ; }
            }
          
          if (ga1 > ga2 && up->x1 < up->x2)  /* case 6: end of up 5' read */
            { 
              cross = FALSE ; 
              vp = up ;
              while (!cross && --vp >= vpMin && vp->cDNA_clone == clone)
                if (vp->est != oldEst && vp->x1 > vp->x2 && vp->a1 < up->a2 && vp->a2 > up->a2)
                  cross = TRUE ;
              vp = up ;
              while (!cross && ++vp <= vpMax && vp->cDNA_clone == clone)
                if (vp->est != oldEst && vp->x1 > vp->x2 && vp->a1 < up->a2 && vp->a2 > up->a2)
                  cross = TRUE ;
              if (cross)
                continue ;

              longDna = s2m->dnaR ; 
              a1 = (long)arrayMax(longDna) - (ga1 - up->a1 + 1) + 1 ;
              a2 = olda2 = (long)arrayMax(longDna) - (ga1 - up->a2 + 1) + 1 ;

              shortDna = dnaD ;
              x1 = up->x1 ;
              x2 = keySet(clipEnds, KEYKEY (up->est)) ; oldx2= up->x2 ;
              
              a1-- ; x1-- ;
              err = aceDnaTrackErrors(shortDna, x1, &x2, longDna, a1, &a2, 0, err, maxJump, 0, TRUE, 0, FALSE) ;
              
              if (x2 > oldx2 && a2 > olda2)
                { up->x2 += (x2 - oldx2) ; up->a2 += (a2 - olda2) ; }
            }
          
          if (ga1 > ga2 && up->x1 > up->x2)  /* case 8: begin of down 3' read */
            {   
              longDna = s2m->dnaR ;
              a1 = (long)arrayMax(longDna) - (ga1 - up->a1 + 1) + 1 ;
              a2 = olda2 = (long)arrayMax(longDna) - (ga1 - up->a2 + 1) + 1 ;

              if (!dnaR) { dnaR = dnaCopy(dnaD) ; reverseComplement (dnaR) ; }
              shortDna = dnaR ;
              x1 = (long)arrayMax (dnaR) - up->x2 + 1 ;
              oldx2 = (long)arrayMax (dnaR) - up->x1 + 1 ;
              x2 = (long)arrayMax (dnaR) - keySet(clipTops, KEYKEY(up->est)) + 1 ;
              
              a1-- ; x1-- ;
              err = aceDnaTrackErrors(shortDna, x1, &x2, longDna, a1, &a2, 0, err, maxJump, 0, TRUE, 0, FALSE) ;
              
              if (x2 > oldx2 && a2 > olda2)
                { up->x1 -= (x2 - oldx2) ; up->a2 += (a2 - olda2) ; }
            }
        }
      oldEst = up->est ;
    }
  arrayDestroy (dnaR) ;
  arrayDestroy (err) ;
}


/*********************************************************************/
/*********************************************************************/
/*   realignement des cDNA sur le genomique */


/*********************************************************************/
static int reportDna = 0, reportBase = 0 ;

void monDnaReport (int *nDnap, int *nBasep)
{ *nDnap = reportDna ; *nBasep = reportBase ;
}

Array getSeqDna(KEY key)
{
  static Associator ass = 0 ;
  void *vp, *vq ;
  Array a = 0 ;
  
  if (key == KEYMAKE (_VCalcul, 12345))
    {
      if (!ass) return 0 ;
      vp = vq = 0 ;
      while (assNext(ass, &vp, &vq))
        {
          a = (Array) vq ;
          if (a == (Array) assVoid(1))
            continue ;
          reportDna-- ; reportBase -= (long)arrayMax(a) ;
          if (arrayExists(a))
            arrayDestroy (a) ;
          else 
            messcrash ("getSeqDna(%s) misused", name(assInt(vp))) ;
        }
      assDestroy (ass) ; 
      return 0 ;
    }
  if (!key)
    return 0 ;
  if (!ass) ass = assCreate() ;
  if (assFind (ass, assVoid(key), &vq))
    {
      if (vq == assVoid(1))
        return 0 ;
      a = (Array) vq ;
      if (!arrayExists (a))
        messcrash ("getSeqDna(%s) misused", name(key)) ;
      return a ;
    }
  a = dnaGet (key) ; 
  if (a)
    { reportDna++ ; reportBase += (long)arrayMax (a) ; }
  vq = a ? (void*)a : assVoid(1) ;
  assInsert (ass, assVoid(key), vq) ;
  return a;
}

/**********************************************************************/

static BOOL checkOligoEntropy (int nn, unsigned int oligo)
{
  int na = 0, nc = 0, nt = 0, ng = 0, i = nn ;
  static int ee15[16],  ee12[13], ee6[7] ;
  static BOOL isFirst = TRUE, ok = TRUE ;

  if (isFirst)
    {
      double s = 0, log4 = log(4.0) ; /* divide by log4 to be in base equivalent */ 
      int j ;
      for (j = 1 ; j <= 15 ; j++)
        { s = j ; s = s/15.0 ; ee15[j] = - (int) (1000.0 * j * log(s)/log4 ) ; }
      ee15[0] = 0 ;
      for (j = 1 ; j <= 12 ; j++)
        { s = j ; s = s/12.0 ; ee12[j] = - (int) (1000.0 * j * log(s)/log4 ) ; }
      ee12[0] = 0 ;
      for (j = 1 ; j <= 6 ; j++)
        { s = j ; s = s/6.0 ; ee6[j] = - (int) (1000.0 * j * log(s)/log4 ) ; }
      ee6[0] = 0 ;
      isFirst = FALSE ;
    }
  i = nn ; 
  while (i--)
    {
      switch ( oligo & 0x03)
        {
        case 0x00: na++ ; break ;
        case 0x01: ng++ ; break ;
        case 0x02: nc++ ; break ;
        case 0x03: nt++ ; break ;
        }
      oligo >>= 2 ;
    }
  switch (nn)
    {
    case 15:
      i = (ee15[na] + ee15[ng] + ee15[nc] + ee15[nt]) / 1000 ;
      ok = i >= 9.0 ? TRUE : FALSE ; break ;
    case 12:
      i = (ee12[na] + ee12[ng] + ee12[nc] + ee12[nt]) / 1000 ;
      return i >= 7.0 ? TRUE : FALSE ; break ;
    case 6:
      i = (ee6[na] + ee6[ng] + ee6[nc] + ee6[nt]) / 1000 ;
      ok = i >= 4.0 ? TRUE : FALSE ; break ;
    default:
      ok = TRUE ; break ;
    }
  return ok ;
}

/**********************************************************************/

static int assSequence (Array words, 
                        KEY key, BOOL isUp, int nn /* OLIGO_LEMGTH */,
                        Associator ass, int step, int minPos,
                        KEYSET cDNA_clones, KEYSET clipTops, KEYSET clipEnds, KEYSET trackingClones, 
			KEYSET estMaps, KEYSET estMap1, KEYSET estMap2, 
			KEY chrom, int c1, int c2, 
			int type
			)
{
  int ii, pos, pos1, max, maxpos, j, n1 = 0, eshift = 0, vtop, vend, vendold, nw, nwmin, nwmax ;
  int ipa = 0, ipb = 0 ;
  int leftShift = 2 * (nn - 1) ;
  unsigned int oligo, oligoreverse, anchor, smask, word, wordpos, 
    mask = (1<<(2 * nn)) -1 , keyMask = (1<<22) -1, wordMask = (1<<30) - 1 ;
  unsigned int  leftMask = mask & ~(3 << leftShift) ;
  unsigned char wordchar ;
  Array aa = 0, ctm = 0, dna = 0 ; /* NOT to be stored in local dna cache ! */
  KEY eMap = 0, eMap1 = 0, eMap2 = 0, cDNA_clone = 0 ;
  char *cp ;
  OBJ obj = 0 ;
  KEY topFlag = 0, endFlag = 0 ;
  BOOL isDifficult = FALSE ;
  BOOL debug = FALSE ;

  if (nn > 15) messcrash ("Oligo_length > 15 in assSequence") ;
  
  obj = bsUpdate (key) ;
  bsGetKey (obj, _cDNA_clone, &cDNA_clone) ;
  
  isDifficult = bsFindTag (obj, str2tag("Difficult_to_align"));

  if (obj && type != 9)
    {
      KEY map= 0 ;
      int e1, e2 ;
      
      if (bsGetKey (obj, _IntMap, &map) &&
	  bsGetData (obj, _bsRight, _Int, &e1) &&
	  bsGetData (obj, _bsRight, _Int, &e2) &&
	  e1 > 0 && e2 > 0
	  )
	{
	  eMap = map ; eMap1 = (KEY)e1 ; eMap2 = (KEY)e2 ; 
	}
    }

  if (eMap && chrom)
    {
      BOOL isForward = TRUE ;
      
      if (bsFindTag (obj, _Reverse))
	isForward = FALSE ;
      if (c1 > c2)
	{ int c0 = c1 ; c1 = c2 ; c2 = c0 ; isForward = ! isForward ; }
      if (
	  (eMap != chrom) ||
	  (
	   ! isUp &&  /* forward */
	   ( (isForward && eMap1 > eMap2) ||  (! isForward && eMap1 < eMap2)) 
	   ) ||
	  (
	   isUp &&  /* reverse */
	   ( (! isForward && eMap1 > eMap2) ||  (isForward && eMap1 < eMap2)) 
	   ) ||
	  ( eMap1 < eMap2 && (eMap1 > c2 || eMap2 < c1)) ||
	  ( eMap1 > eMap2 && (eMap2 > c2 || eMap1 < c1))
	  )
	goto abort ;
    }

  dna =  dnaGet(key) ; 
  if (!dna)
    goto abort ;

  vtop = 1 ;
  vend = max = arrayMax(dna) ;

  if (max < 100)
    isDifficult = TRUE ;

  ii = 0 ;
  if (obj && type != 9 &&
      bsGetData (obj, _Vector_clipping, _Int, &vtop) &&
      bsGetData (obj, _bsRight, _Int, &vend))
    {
      if (vtop > max) vtop = max ;
      if (vtop < 1) vtop = 1 ;
      if (ii && vend > ii) vend = ii ;
      if (vend < 20 || vend > max) vend = max ;
    }

  if (!isUp && obj &&
      bsGetData (obj, str2tag("Initial_polyA"), _Int, &ipa) &&
      bsGetData (obj, _bsRight, _Int, &ipb))
    {
      ipa = ipa + ipb ;
      if (ipa > vtop && ipa < vtop + 200)
        vtop = ipa ;
      if (vtop > max) vtop = max ;
    }
  if (!isUp && obj &&
      bsGetData (obj, str2tag("Initial_polyT"), _Int, &ipa) &&
      bsGetData (obj, _bsRight, _Int, &ipb))
    {
      ipa = ipa + ipb ;
      if (ipa > vtop && ipa < vtop + 200)
        vtop = ipa ;
      if (vtop > max) vtop = max ;
    }

  ii = 0 ;
  if (!isUp && obj &&
      ! keyFindTag (cDNA_clone, _Internal_priming) &&
      bsGetData (obj, _PolyA_after_base, _Int, &ii))
    {
      if (ii < vtop + 50 && ! bsFindTag (obj, _Composite))
        ii = 0 ; /* limit 200 to prevent problem with bad orientations */
      if (ii > 0 && vend > ii) 
        { vend = ii ; endFlag = KEYMAKE (1,0) ; } 
      if (ii > 0 && vtop > ii) 
        { vtop = ii ; }
    }

  ctm = arrayCreate (36, BSunit) ;
  if (!obj || 
      !bsGetArray (obj, _Contamination, ctm, 2) ||
      !arrayMax (ctm))
    arrayDestroy (ctm) ;

  if (ctm)
    {
      int iCtm ;
      BSunit *uc ;

      for (iCtm = 0 ; ctm && iCtm < arrayMax(ctm) ; iCtm += 2)
        {
          uc = arrp(ctm, iCtm, BSunit) ;
          if (uc[0].i <= 10 && uc[1].i > vtop)
            vtop = uc[1].i ;
          if (uc[0].i < vend && uc[1].i > vend - 50)
            vend = uc[1].i ;
        }
    }
  aa = arrayCreate (6, BSunit) ;
  if (obj && bsGetArray (obj, _Transpliced_to, aa, 2))
    {
      for (ii = 0, j = 0 ; j < arrayMax(aa) ; j += 2)
        if (!ii || 
            (!isUp && arrp(aa, j + 1, BSunit)->i > ii) ||
            (isUp && arrp(aa, j + 1, BSunit)->i < ii)
            )
          ii = arrp(aa, j + 1, BSunit)->i ;
      if (ii > 0)
        {
          if (!isUp && vtop < ii + 1) 
            { vtop = ii + 1 ; topFlag = 2 ; }
          if (isUp && vend > ii)
            { vend = ii - 1 ; endFlag = 2 ; }
        }
    }
  arrayDestroy (aa) ;

  if (isUp && obj)
    {
      if (bsFindTag (obj, str2tag("Real_5prime")))
        topFlag = 2 ;
      if (! keyFindTag (cDNA_clone, _Internal_priming) &&
          bsGetData (obj, _PolyA_after_base, _Int, &ii))
        {
          if (ii > 150 && ! bsFindTag (obj, _Manual_polyA) && ! bsFindTag (obj, _Composite))
            ii = 0 ; /* limit 200 to prevent problem with bad orientations */
          if (ii > 0 && vtop < ii)
            { vtop = ii ; topFlag = 2 ; }
        }
      else if (! keyFindTag (cDNA_clone, _Internal_priming) &&
               getAddPolyA ())
        {
          ii = 0 ;
          cp = arrp (dna , 0, char) - 1 ;
          while (*++cp == T_) ii++ ;
          
          if (vtop < ii)
            {
              vtop = ii ; topFlag = 2 ;
              if (ii > 6)
                bsAddData (obj, _PolyA_after_base, _Int, &ii) ;
            }
        }
    }
        
  bsSave (obj) ;
  
  /* vtop /vend = the actual sequence 
   * check for startches of n 
   */
 cutAgain: /* CF y2L761b8.3 vecteur start=543, pure N after 585, length = 980 */
  j = 0 ; vendold = vend ;
  if (2 * arrayMax (dna) > vtop + vend)
    for (pos = pos1 = ii = (vtop+vend)/2 , cp = arrp (dna, ii, char) ; ii < vend ; ii++, cp++)
      {
	switch ((int)*cp)
	  {
	  case A_: case T_: case G_: case C_: j-- ; break ;
	  default: j += 2 ; break ;
	  }
	if (j < 0) { j = 0 ; pos = ii ; }
	if (j > 20) /* about 10 N in a stretch */
	  { vend = pos + 5 ; break ; }
      }
  if (vend < vendold - 20 && pos < pos1 + 20 && vend > vtop + 60) goto cutAgain ;
  j = 0 ;
  if (2 * arrayMax (dna) > vtop + vend)
    for (pos = ii = (vtop+vend)/2 , cp = arrp (dna, ii, char) ; ii > vtop ; ii--, cp--)
      {
	switch ((int)*cp)
	  {
	  case A_: case T_: case G_: case C_: j-- ; break ;
	  default: j += 2 ; break ;
	  }
	if (j < 0) { j = 0 ;pos = ii ; }
	if (j > 20) /* about 10 N in a stretch */
	  { vtop = pos - 5 ; break ; }
      }
  
  if (minPos < vtop) minPos = vtop ;
  if (! eMap && minPos + nn + 10 > vend)
    goto abort ;
  
  topFlag = endFlag = 0 ; /* don t activate this yet, 
                           *  if yes, check all usage of clipTops/Ends 
                           */
  if (type != 9)
    {
      if (clipTops) keySet (clipTops, KEYKEY(key)) = KEYMAKE (topFlag, vtop) ;
      if (clipEnds) keySet (clipEnds, KEYKEY(key)) = KEYMAKE (endFlag, vend) ;
      if (estMaps) keySet (estMaps, KEYKEY(key)) = eMap ;
      if (estMap1) keySet (estMap1, KEYKEY(key)) = eMap1 ;
      if (estMap2) keySet (estMap2, KEYKEY(key)) = eMap2 ;

      if (cDNA_clones) 
        {
          keySet (cDNA_clones, KEYKEY(key)) = keyGetKey (key, _cDNA_clone) ;
        }
      if (cDNA_clone && 
          (
           keyFindTag (cDNA_clone, str2tag("Tracking_clone")) || 
           (
            keyFindTag (cDNA_clone, _Suspected_internal_deletion) &&
            !keyFindTag (cDNA_clone, _Manual_no_internal_deletion) &&
            !keyGetKey (cDNA_clone, _Suspected_internal_deletion)  /* in this subcase, just delete one intron */
            ) ||
           keyFindTag (cDNA_clone, str2tag("Duplicate_clone")) ||
           keyFindTag (cDNA_clone, str2tag("Ignore_this_clone_automatic")) ||
           keyFindTag (cDNA_clone, str2tag("Ignore_this_clone"))
           )
          )
        keySetInsert (trackingClones, cDNA_clone) ;
    }
  else
    {
      vtop = keySet (clipTops, KEYKEY(key)) ;
      vend = keySet (clipEnds, KEYKEY(key)) ;
    }

  /* eshift, smask is used to anchor seqs longer than 1024
   * only 10 bits can be spared for the anchoring
   * for longer seqs, i position modulo smask
   */ 
  smask = vend ;
  smask >>= 10 ; eshift = 0 ;
  while (smask) { smask >>= 1 ; eshift++ ; }
  smask = (1<<eshift) - 1 ;
  if (step < (1<<eshift)) step = 1 << eshift ;
  
  maxpos = vend - nn ;
  if (isDifficult) 
    {
      step = 3 ;
      if (maxpos < 1500 && maxpos > 500) maxpos = 500 ;
    }
  if (isUp) 
    { 
      reverseComplement (dna) ;
      pos = minPos ;
      minPos = max - maxpos ;
      maxpos = max - pos ;
    }

  if (step < (maxpos - minPos) / 20)
    step = (maxpos - minPos) / 20 ;
  if (step > 30)  /* adequate for mrna which are > 2000; size of an exon */
    step = 30 ;

  if (isDifficult) 
    step = 3 ;
  nwmin =  (maxpos - minPos) / step ;
  nwmax = 2 ;
  pos = minPos ;
  ii = -10000000 ;
  while (pos < maxpos)
    {
      if (ctm)
        {
          BOOL isCtm = FALSE ;
          BSunit *uc ;
          int iCtm ;

          for (iCtm = 0 ; !isCtm && iCtm < arrayMax(ctm) ; iCtm += 2)
            {
              uc = arrp(ctm, iCtm, BSunit) ;
              if (uc[0].i <= pos && uc[1].i >= pos - nn)
                isCtm = TRUE ;
            }
          if (isCtm)
            { pos++ ; continue ; }
        }
      j = nn ; oligo = 0 ; oligoreverse = 0 ;
      pos1 = isUp ? max - pos : pos  ;
      pos1 = ((pos1 + smask) >> eshift) << eshift ; /* +smask prevents forever */
      if (pos1 == ii) { pos++ ; continue ; }
      ii = pos1 ;
      pos = isUp ? max - pos1 : pos1  ;
      if (pos >= maxpos) /* may happen because of eshift rounding */ continue ;
      cp = arrp(dna, pos, char) ;
      while (cp++, j--)
        { 
          if (!*cp || *cp == N_)
            goto suite ;            
          oligo <<= 2 ; oligo |= B2[(int)(*cp)] ; oligo &= mask ; 
          oligoreverse >>= 2 ;  oligoreverse &= leftMask ;
          oligoreverse |=  ((B2r[(int)(*cp)]) << leftShift ); 
        }
      cp -= nn ; /* reset to pos */

      anchor = (((pos1 >> eshift) & 1023) << 22) | (key & keyMask) ;
      nw = 1 ; /* default */
      if (nn == 15 && words)
        {
          word = oligo & wordMask ;
          wordpos = word >> 1 ; 
          wordchar = array(words, wordpos, unsigned char) ;
          if (word & 0x1)
            wordchar >>= 4 ;
          nw = wordchar & ((unsigned char) 0x0f) ;

          if (nw < nwmax)
            {
              word = oligoreverse & wordMask ;
              wordpos = word >> 1 ; 
              wordchar = array(words, wordpos, unsigned char) ;
              if (word & 0x1) 
                wordchar >>= 4 ;
              nw += wordchar & ((unsigned char) 0x0f) ;
            }
        }
      if (oligo   && oligoreverse && /* not the poly A */
          nw > 0 && nw < nwmax &&
          checkOligoEntropy(nn, oligo))
        { 
          assMultipleInsert (ass, assVoid(oligo), assVoid(anchor)) ; n1++ ; 
          if (debug) printf("%s %d nwmax=%d\n", name(key), pos1, nwmax) ; 
          pos += step ;  continue ;
        }
      else if (debug)
        printf("NO: %s %d \n", name(key), pos) ; 
          
      if (pos > 4 * step && n1 < nwmin && nwmax < 12)
        { pos = minPos ; nwmax++ ; ii = -10000000 ; }
    suite:
      pos++ ;
    }
  
abort:
  bsSave (obj) ;
  arrayDestroy (ctm) ;
  arrayDestroy (dna) ;
  return n1 ;
}

/*********************************************************************/
/*********************************************************************/
/* we know there is an exact hit, try to extend it
   a1, a2 are cosmid coord, x1, x2, cDNA coord 
   */

static BOOL extendHits2 (KEY cosmid, Array dna, Array dnaR,
                        KEY est, Array dnaEst,
                        int *a1p, int *a2p, int *x1p, int *x2p,
                         int v1, int v2, int *nerrp, int maxJump, BOOL doExtend, BOOL countN, int *maxExactp)
{
  BOOL isUp = FALSE ;
  int a1, a2, x1, x2, u1, u2, nn, nN = 0 ;
  Array dnaLong = 0 ;
  static Array aa = 0 ;
  int maxError = 0 ;

  x1 = *x1p ; x2 = *x2p ; a1 = *a1p ; a2 = *a2p ;
  if (a1 > a2) isUp = TRUE ;

  nn = arrayMax (dna) ;
  if (isUp)
    { dnaLong = dnaR ; u1 = nn - a1 + 1 ; u2 = nn - a2 + 1; }
  else
    { dnaLong = dna ; u1 = a1 ; u2 = a2 ; }
  if (doExtend) x2 = v2 ;
  if (x1 < v1) x1 = v1 ;
  x1-- ; u1-- ; /* plato */
  chrono("extendHits2") ; chronoReturn () ;
  if (keyFindTag (est, _Composite))
    {
      maxError = -2 ; maxJump = 0 ; 
    }
  aa = aceDnaTrackErrors (dnaEst, x1, &x2,
                    dnaLong, u1, &u2, &nN, aa, maxJump, maxError, FALSE, maxExactp, FALSE) ;
  x1++ ; u1++ ; /* plato */
  if (isUp) a2 = nn - u2 + 1 ; 
  else  a2 = u2 ;
  
  *x1p = x1 ; *x2p = x2 ; *a1p = a1 ; *a2p = a2 ;
  *nerrp = (long)arrayMax (aa) - (countN ? 0 : nN) ;
  return TRUE ;
}

/*********************************************************************/

static BOOL extendHits (KEY cosmid, Array dna, Array dnaR, HIT *vp)
{
  BOOL estUp = FALSE, ok, debug = FALSE ;
  int nerr1 = 0, nerr2 = 0, a1, a2, a3, x1, x2, x3, v1, v2, v3, nn, b2, y2 ;
  Array dnaEst = 0 ;
  
  /* placer les coords dans l'ordre naturel */
  a1 = vp->a1 ; a2 = vp->a2 ; x1 = vp->x1 ; x2 = vp->x2 ;
  if (x1 > x2) 
    { 
      estUp = TRUE ;
      x3 = x1 ; x1 = x2 ; x2 = x3 ; 
      a3 = a1 ; a1 = a2 ; a2 = a3 ;
    }
  else
    estUp = FALSE ;

  dnaEst = getSeqDna (vp->est) ;
  if (!dnaEst) return FALSE ;
  nn = arrayMax(dnaEst) ;
  if (debug) printf("extend 1:: a1 = %d a2 = %d x1 = %d x2 = %d\n",
                    a1, a2, x1, x2) ;

  v1 = vp->clipTop ; v2 = vp->clipEnd ;
  vp->maxExact = 0 ;
  extendHits2(cosmid, dna,dnaR,vp->est,dnaEst,&a1,&a2,&x1,&x2,v1,v2, &nerr1, 2, TRUE, FALSE, &(vp->maxExact)) ;
  if (debug) printf("extend 2:: a1 = %d a2 = %d x1 = %d x2 = %d\n",
                    a1, a2, x1, x2) ;

  ok = x1 < 2 ? TRUE : FALSE ;
  /* retablir */
  if (estUp) 
    {  x3 = x1 ; x1 = x2 ; x2 = x3 ; 
       a3 = a1 ; a1 = a2 ; a2 = a3 ; 
    }
  b2 = a2 ; y2 = x2 ;
  if (ok)
    { nerr1 *= 2 ; goto done ; }

  /* extend backwards, starting from the exact oligo  */
  a1 = vp->a1 ; a2 = vp->a2 ; x1 = vp->x1 ; x2 = vp->x2 ;

    /* inverser */
  reverseComplement (dnaEst) ;
  x3 = x1 ; x1 = nn - x2 + 1 ; x2 = nn - x3 + 1 ;
  a3 = a1 ; a1 = a2 ; a2 = a3 ;
  v3 = v1 ; v1 = nn - v2 + 1 ; v2 = nn - v3 + 1 ;

  /* placer les coords dans l'ordre naturel */
  if (x1 > x2) 
    { estUp = TRUE ; x3 = x1 ; x1 = x2 ; x2 = x3 ; 
    a3 = a1 ; a1 = a2 ; a2 = a3 ; }
  else
    estUp = FALSE ;

  if (debug) printf("extend 3:: a1 = %d a2 = %d x1 = %d x2 = %d\n",
                    a1, a2, x1, x2) ;
  extendHits2(cosmid, dna,dnaR,vp->est,dnaEst,&a1,&a2,&x1,&x2,v1, v2, &nerr2, 2, TRUE, FALSE, &(vp->maxExact)) ;
  if (debug) printf("extend 4:: a1 = %d a2 = %d x1 = %d x2 = %d\n",
                    a1, a2, x1, x2) ;
  
  /* retablir */
  if (estUp) 
    {  x3 = x1 ; x1 = x2 ; x2 = x3 ; 
       a3 = a1 ; a1 = a2 ; a2 = a3 ;  
    }
  reverseComplement (dnaEst) ;
  x3 = x1 ; x1 = nn - x2 + 1 ; x2 = nn - x3 + 1 ;
  a3 = a1 ; a1 = a2 ; a2 = a3 ;
 done:
  if (debug) printf("extend 5:: a1 = %d a2 = %d x1 = %d x2 = %d\n",
                    a1, b2, x1, y2) ;
  if (vp->x1 < vp->x2)
    {
      if (y2 > vp->x2)
        {  vp->a2 = b2 ; vp->x2 = y2 ; }
      if (x1 < vp->x1)
        {  vp->a1 = a1 ; vp->x1 = x1 ; }
    }
  else
    {
      if (y2 < vp->x2)
        {  vp->a2 = b2 ; vp->x2 = y2 ; }
      if (x1 > vp->x1)
        {  vp->a1 = a1 ; vp->x1 = x1 ; }
    }
  vp->nerr = (nerr1 + nerr2)/2 ;

  return TRUE ;
}

/*********************************************************************/

static BOOL countOneErrors (KEY cosmid, Array dna, Array dnaR, HIT *vp, BOOL countN)
{
  BOOL estUp = FALSE, ok, debug = FALSE ;
  int nerr1, nerr2, a1, a2, a3, x1, x2, x3, v1, v2, v3, nn, b2, y2, dx, dx1, dx2, maxExact = 0 ;
  Array dnaEst = 0 ;
  
  if (vp->a1 == vp->ea1 && vp->a2 == vp->ea2 && (vp->x1 + (countN ? (1 << 30) : 0))== vp->ex1  && vp->x2 == vp->ex2)
    {
      chrono ("countOneErrorsRepeat") ; chronoReturn () ;
      return TRUE ; /* already computed */
    } 
  chrono ("countOneErrors") ;
  /* placer les coords dans l'ordre naturel */
  a1 = vp->a1 ; a2 = vp->a2 ; x1 = vp->x1 ; x2 = vp->x2 ;
  if (x1 > x2) 
    { 
      estUp = TRUE ;
      x3 = x1 ; x1 = x2 ; x2 = x3 ; 
      a3 = a1 ; a1 = a2 ; a2 = a3 ;
    }
  else
    estUp = FALSE ;

  dnaEst = getSeqDna (vp->est) ;
  if (!dnaEst) return FALSE ;
  nn = arrayMax(dnaEst) ;
  if (debug) printf("extend 1:: a1 = %d a2 = %d x1 = %d x2 = %d\n",
                    a1, a2, x1, x2) ;
  dx = x2 - x1 ;
  v1 = vp->clipTop ; v2 = vp->clipEnd ;
  nerr1 = nerr2 = 0 ;
  extendHits2(cosmid, dna,dnaR,vp->est,dnaEst,&a1,&a2,&x1,&x2,v1,v2,&nerr1,8, FALSE,countN, &maxExact) ;
  if (debug) printf("extend 2:: a1 = %d a2 = %d x1 = %d x2 = %d\n",
                    a1, a2, x1, x2) ;
  dx1 = x2 - x1 ; dx2 = 0 ;
  ok =  (dx1 >= dx) ? TRUE : FALSE ;
  /* we are just counting, we measure twice but only in bad cases where we failed to track */
  /* retablir */
  if (estUp) 
    {  x3 = x1 ; x1 = x2 ; x2 = x3 ; 
       a3 = a1 ; a1 = a2 ; a2 = a3 ; 
    }
  b2 = a2 ; y2 = x2 ;

  if (ok)
    { nerr1 *= 2 ; goto done ; }

  /* extend backwards, starting from the exact oligo  */
  a1 = vp->a1 ; a2 = vp->a2 ; x1 = vp->x1 ; x2 = vp->x2 ;

    /* inverser */
  reverseComplement (dnaEst) ;
  x3 = x1 ; x1 = nn - x2 + 1 ; x2 = nn - x3 + 1 ;
  a3 = a1 ; a1 = a2 ; a2 = a3 ;
  v3 = v1 ; v1 = nn - v2 + 1 ; v2 = nn - v3 + 1 ;

  /* placer les coords dans l'ordre naturel */
  if (x1 > x2) 
    { estUp = TRUE ; x3 = x1 ; x1 = x2 ; x2 = x3 ; 
    a3 = a1 ; a1 = a2 ; a2 = a3 ; }
  else
    estUp = FALSE ;

  if (debug) printf("extend 3:: a1 = %d a2 = %d x1 = %d x2 = %d\n",
                    a1, a2, x1, x2) ;
  extendHits2(cosmid, dna,dnaR,vp->est,dnaEst,&a1,&a2,&x1,&x2,v1, v2, &nerr2, 8, FALSE, countN, &maxExact) ;
  if (debug) printf("extend 4:: a1 = %d a2 = %d x1 = %d x2 = %d\n",
                    a1, a2, x1, x2) ;
  dx2 = x2 - x1 ;
  /* retablir */
  if (estUp) 
    {  x3 = x1 ; x1 = x2 ; x2 = x3 ; 
       a3 = a1 ; a1 = a2 ; a2 = a3 ;  
    }
  reverseComplement (dnaEst) ;
  x3 = x1 ; x1 = nn - x2 + 1 ; x2 = nn - x3 + 1 ;
  a3 = a1 ; a1 = a2 ; a2 = a3 ;
 done:
  if (debug) printf("extend 5:: a1 = %d a2 = %d x1 = %d x2 = %d\n",
                    a1, b2, x1, y2) ;
  vp->maxExact = maxExact ;
  vp->nerr = (nerr1 + nerr2)/2 ;
  nn = dx - dx1 - dx2 ;
  if (nn > 50) nn = 50 ; /* over 50, happens in long reads, my tracking went nuts */
  if (nn > 0) vp->nerr += nn ;
  vp->ea1 = vp->a1 ; vp->ea2 = vp->a2 ; vp->ex1 = vp->x1 ; vp->ex2 = vp->x2 ; 
  if (countN) vp->ex1 += (1 << 30) ;
  chronoReturn () ;
  return TRUE ;
}

/***************************/

static void countErrors (KEY cosmid, Array hits, Array dna, Array dnaR, int maxX, BOOL countN)
{
  int ii = arrayMax(hits), x2, a2, dx ;
  Array dnaEst = 0 ;
  BOOL mm = FALSE ;
  HIT *up ;

  cDNASwapX (hits) ;
  arraySort (hits, cDNAOrderByX1) ;
  
  while (ii--)
    {
      up = arrp (hits, ii, HIT) ;
      if (up->est == 1) /* ghost */
        continue ;
      if (up->x1 < up->clipTop)
        { 
          dx = up->clipTop - up->x1 ;
          if (up->a1 < up->a2) up->a1 += dx ;
          else up->a1 -= dx ;
          up->x1 += dx ;
        }
      if (up->x2 > up->clipEnd)
        { 
          dx = up->x2 - up->clipEnd ;
          if (up->a1 < up->a2) up->a2 -= dx ;
          else up->a2 += dx ;
          up->x2 -= dx ;
        }

      if ((up->type & gI) ||
          up->est == 1 ||
          (up->x1 > up->x2 - 5))
        { up->nerrAll = up->nerr = 0 ; continue ; }

      if (keyFindTag (up->est, _mRNA) ||
          keyFindTag (up->est, str2tag ("Ref_mRNA")) || 
          keyFindTag (up->est, str2tag ("Ref_seq")) || 
          ((dnaEst = getSeqDna (up->est)) && arrayMax(dnaEst) > 1500))
        mm = TRUE ;
      else
        mm = FALSE ;

      x2 = up->x2 ; a2 = up->a2 ;
      countOneErrors (cosmid, dna, dnaR, up, countN) ;
      up->nerrAll = up->nerr ;   /* always count all errors */
      dx = 0 ;
      if (!mm && maxX && up->x2 > maxX && up->nerr) /* count again in restricted mode */
        {
          dx = up->x2 - maxX ;
          up->x2 -= dx ;
          if (up->a1 < up->a2)
            up->a2 -= dx ;
          else
            up->a2 += dx ;
        }
      if (dx)
        {
          if (up->x1 < up->x2 - 3) /* may be wrong with maxX shift */
            countOneErrors (cosmid, dna, dnaR, up, countN) ;
          else
            up->nerr = 0 ;
          /* restore */
          up->x2 = x2 ; up->a2 = a2 ; 
        }

    }
}
    
/*********************************************************************/
/*********************************************************************/

static int cDNAOrderEstByA1 (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;

  if (up->gene != vp->gene)
    return up->gene - vp->gene ;
  if (up->a1 != vp->a1)
    return up->a1 - vp->a1 ;
  else if (up->a2 != vp->a2)
    return up->a2 - vp->a2 ;      
  else  
    return up->x1 - vp->x1 ;
}

/*********************************************************************/

int cDNAOrderGloballyByA1 (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;

  if (up->a1 != vp->a1)
    return up->a1 - vp->a1 ;
  else if (up->a2 != vp->a2)
    return up->a2 - vp->a2 ;      
  else if (up->reverse != vp->reverse)  /* order by strand */
    return up->reverse ? 1 : -1 ;
  else 
    return up->est - vp->est ;
}

/*********************************************************************/

int cDNAOrderGloballyByA1Errors (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;

  if (up->a1 != vp->a1)
    return up->a1 - vp->a1 ;
  else if (up->a2 != vp->a2)
    return up->a2 - vp->a2 ;      
  else if (up->reverse != vp->reverse)  /* order by strand */
    return up->reverse ? 1 : -1 ;
  else if (up->nerr != vp->nerr)
    return up->nerr - vp->nerr ;      
  else 
    return up->est - vp->est ;
}

/*********************************************************************/

static int cDNAOrderGenesGloballyByA1 (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;

  if (up->reverse != vp->reverse)  /* order by strand */
    return up->reverse ? 1 : -1 ;

  if (up->a1 != vp->a1)
    return up->a1 - vp->a1 ;
  else if (up->a2 != vp->a2)
    return up->a2 - vp->a2 ;      
  else  
    return up->type - vp->type ;    
}

/*********************************************************************/

static int cDNAOrderGenesByA1 (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;

  if (up->reverse != vp->reverse)  /* order by strand */
    return up->reverse ? 1 : -1 ;

  if (!up->reverse)
    {
      if (up->a1 != vp->a1)
        return up->a1 - vp->a1 ;
      else if (up->a2 != vp->a2)
        return up->a2 - vp->a2 ;      
      else  
        return up->type - vp->type ;    
    }
  else
    {
      if (up->a2 != vp->a2)
        return up->a2 - vp->a2 ;      
      else if (up->a1 != vp->a1)
        return up->a1 - vp->a1 ;
      else  
        return up->type - vp->type ;    
    }
}

/*********************************************************************/
/* used for debugging */
#ifdef JUNK
static int cDNAOrderBySize (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;
  int da, db ;

  if (up->reverse != vp->reverse)  /* order by strand */
    return up->reverse ? 1 : -1 ;

  da = up->a2 > up->a1 ? up->a2 - up->a1 : up->a1 - up->a2 ;
  db = vp->a2 > vp->a1 ? vp->a2 - vp->a1 : vp->a1 - vp->a2 ;

  if (0) cDNAOrderBySize (0,0) ; /* to please the compiler */
  return db - da ; /* long first */
}
#endif
/*********************************************************************/

int cDNAOrderIntronsTeams (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;

  if (up->nerr != vp->nerr)
    return vp->nerr - up->nerr ; /* priviledge frequently used introns */
  if (up->cDNA_clone != vp->cDNA_clone)
    return up->cDNA_clone - vp->cDNA_clone ;  /* the cDNA_clone */
  else if (up->est != vp->est)
    return up->est - vp->est ;  /* the cDNA_clone */
  else if (up->a1 != vp->a1)
    return up->a1 - vp->a1 ;
  else if (up->a2 != vp->a2)
    return up->a2 - vp->a2 ;      
  else  
    return up->x1 - vp->x1 ;    
}

/*********************************************************************/

int cDNAOrderByA1 (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;

  if (up->cDNA_clone != vp->cDNA_clone)
    return up->cDNA_clone - vp->cDNA_clone ;  /* the cDNA_clone */
  else if (up->reverse != vp->reverse)
    return up->reverse - vp->reverse ;  /* the cDNA_clone */
  else if (up->est != vp->est)
    return up->est - vp->est ;  /* the ESt */
  else if (up->a1 != vp->a1)
    return up->a1 - vp->a1 ;
  else if (up->a2 != vp->a2)
    return up->a2 - vp->a2 ;      
  else  
    return up->x1 - vp->x1 ;    
}

/*********************************************************************/

int cDNAOrderByX1 (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;

  if (up->cDNA_clone != vp->cDNA_clone)
    return up->cDNA_clone - vp->cDNA_clone ;  /* the cDNA_clone */
  else if (up->est != vp->est)
    return up->est - vp->est ;  /* the cDNA_clone */
  else if (up->x1 != vp->x1)
    return up->x1 - vp->x1 ;
  else
    return up->x2 - vp->x2 ;
}

/*********************************************************************/

int cDNAOrderByGene (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;

  if (up->gene != vp->gene)
    return up->gene - vp->gene ;  /* the gene */
  return up->cDNA_clone - vp->cDNA_clone ;  /* the cDNA_clone */
}

/*********************************************************************/
/* in extend multiple hits, i try to contract the gene */
static int cDNAOrderByX1Merge (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;

  if (up->cDNA_clone != vp->cDNA_clone)
    return up->cDNA_clone - vp->cDNA_clone ;  /* the cDNA_clone */
  else if (up->est != vp->est)   /* one read then the other */
    return up->est - vp->est ;
  else if ((up->a2 - up->a1) * (vp->a2 - vp->a1) < 0)
    return (up->a2 - up->a1) > 0 ? -1 : 1 ; /* one strand then the other */
  else if (up->x1 != vp->x1)   /* nearest hit along read */
    return up->x1 - vp->x1 ;
  else if (up->x2 != vp->x2)   /* then longuest, big x2 first */
    return -up->x2 + vp->x2 ;
  else if((up->a2 - up->a1) * (up->x2 - up->x1) > 0) /* parallel */
    return vp->a1 - up->a1 ;
  else
    return up->a1 - vp->a1 ;
}

/*********************************************************************/

/* i need to sort by cDNA_clone */
void cDNASwapX (Array hits)
{
  int tmp, i ;
  HIT *vp ;

  for (i = 0 ; i < arrayMax(hits) ; i++)
    { 
      vp = arrp (hits, i, HIT) ; 
      if (vp->x1 > vp->x2)
        {
          tmp = vp->x1 ; vp->x1 = vp->x2 ; vp->x2 = tmp ;
          tmp = vp->a1 ; vp->a1 = vp->a2 ; vp->a2 = tmp ;
        }
    }
}

/*********************************************************************/

void cDNASwapA (Array hits)
{
  int tmp, i ;
  HIT *vp ;

  for (i = 0 ; i < arrayMax(hits) ; i++)
    { 
      vp = arrp (hits, i, HIT) ; 
      if (vp->a1 > vp->a2)
        {
          tmp = vp->x1 ; vp->x1 = vp->x2 ; vp->x2 = tmp ;
          tmp = vp->a1 ; vp->a1 = vp->a2 ; vp->a2 = tmp ;
        }
    }
}

/*********************************************************************/

static void  getcDNA_clonesAndClips (Array hits, KEYSET cDNA_clones, KEYSET clipTops, KEYSET clipEnds)
{
  int  i ;
  HIT *vp ;

  for (i = 0 ; i < arrayMax(hits) ; i++)
    { 
      vp = arrp (hits, i, HIT) ; 
      vp->cDNA_clone = keySet (cDNA_clones, KEYKEY(vp->est)) ;
      vp->clipTop = keySet (clipTops, KEYKEY(vp->est)) ;
      vp->clipEnd = keySet (clipEnds, KEYKEY(vp->est)) ;
    }
}

/*********************************************************************/

static HIT* cDNA_clone2gene(HIT *ah, Array genes, Array sets, BOOL fromHits)
{
  int j1 = arrayMax(genes) ;
  HIT *gh ;
  KEY cDNA_clone = ah->cDNA_clone ;
  BOOL reverse ;

  if (!fromHits) /* from geneHits */
    reverse = ah->reverse ;
  else
    {
      reverse = (ah->a2 - ah->a1) * (ah->x2 - ah->x1) < 0 ? TRUE : FALSE ;
      if (ah->reverse)   /* 3 prime read */
        reverse = !reverse ;
    }
  while (j1--)
    {
      gh = arrp(genes,j1,HIT) ;
      if (!(gh->type & gGene)) 
        continue ;
      if (
          (
           (!reverse && gh->a1 < gh->a2) ||
           (reverse && gh->a1 > gh->a2) 
           ) && 
          keySetFind (arr(sets, gh->zone,KEYSET), cDNA_clone, 0))
        return gh ;
    }
  return 0 ;
}

/*********************************************************************/

static void  saveHits (KEY cosmid, Array genes, Array hits, Array sets, int direction) 
{
  int j ;
  OBJ Cosmid = 0, Est = 0 ;
  KEY old ; 
  HIT *gh, *up ;

  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;

  Cosmid = bsUpdate(cosmid) ;
  if (direction != -1 && bsFindTag (Cosmid, _Hit))
    bsRemove (Cosmid) ;

  
  for (j = 0 ; j < arrayMax (hits); j ++)
    {
      up = arrp(hits, j, HIT) ;
      if (up->est < 2) /* ghost */ 
        continue ;
      if (genes && sets)
        {
          gh = cDNA_clone2gene (up, genes, sets, TRUE) ;
          if (gh)
            up->gene = gh->gene ;
          else
            continue ;
        }
      bsAddKey (Cosmid, _Hit, up->est) ;
      bsAddData (Cosmid, _bsRight, _Int, &(up->a1));
      bsAddData (Cosmid, _bsRight, _Int, &(up->a2)) ;
      bsAddData (Cosmid, _bsRight, _Int, &(up->x1)) ;
      bsAddData (Cosmid, _bsRight, _Int, &(up->x2)) ;
    }
  bsSave(Cosmid) ;
  
  cDNASwapX (hits) ;
  arraySort (hits, cDNAOrderByX1) ;
  for (j = 0, old = 0; j < arrayMax (hits); j ++)
    {
      up = arrp(hits, j, HIT) ;

      if (up->est != old)
        {
          bsSave (Est) ;
          old = up->est ;
          if (up->est > 1)
            Est = bsUpdate (old) ;
        }
      if (!Est) continue ;        
      bsAddKey (Est, _Hit, cosmid) ;
      bsAddData (Est, _bsRight, _Int, &(up->x1)) ;
      bsAddData (Est, _bsRight, _Int, &(up->x2)) ;
      bsAddData (Est, _bsRight, _Int, &(up->a1));
      bsAddData (Est, _bsRight, _Int, &(up->a2)) ;
    }
  bsSave (Est) ;
}

/*********************************************************************/

static void  saveAlignment (KEY cosmid, Array hits, Array dna, char *nom) 
{
  int i, j ;
  OBJ obj = 0 ;
  KEY alignment = 0, subAlignment = 0 ;
  HIT *up ;

  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;
  cDNASwapX (hits) ;

  i = 0 ;
  if (!nom || !*nom) /* let the system decide, no reuse */
    {
      nom = "alignment" ;
      if (lexword2key (messprintf("%s_%s",name(cosmid), nom), &alignment, _VSequence))
        {
          while (i++, lexword2key (messprintf("%s_%s.%d",name(cosmid), nom, i), &alignment, _VSequence)) ;
          lexaddkey (messprintf("%s_%s.%d",name(cosmid), nom,i), &alignment, _VSequence) ;
        }
      else
        lexaddkey (messprintf("%s_%s",name(cosmid), nom), &alignment, _VSequence) ;
      
      if (lexword2key ("z_1", &subAlignment, _VSequence))
        {
          i = 0 ;
          while (i++, lexword2key (messprintf("z_.%d",i), &subAlignment, _VSequence)) ;
          lexaddkey (messprintf("z_.%d",i), &subAlignment, _VSequence) ;
        }
      else
        lexaddkey ("z_1", &subAlignment, _VSequence) ;
    }
  else  /* reuse */
    {
      lexaddkey (messprintf("%s_%s",name(cosmid), nom), &alignment, _VSequence) ;
      lexaddkey (messprintf("z_%s_%s",name(cosmid), nom), &subAlignment, _VSequence) ;
    }
  
  /*
  lexaddkey(name(subAlignment), &key, _VDNA) ;
  dnaStoreDestroy (key, dnaCopy (dna)) ;  // sets length in subAlignment 
  */
  obj = bsUpdate(subAlignment) ;
  if (bsFindTag (obj, _Assembled_from))
    bsRemove (obj) ;
  bsAddTag (obj, str2tag("Is_Reference")) ;

  for (j = 0 ; j < arrayMax (hits); j ++)
    {
      up = arrp(hits, j, HIT) ;

      bsAddKey (obj, _Assembled_from, up->est) ;
      bsAddData (obj, _bsRight, _Int, &(up->a1));
      bsAddData (obj, _bsRight, _Int, &(up->a2)) ;
      bsAddData (obj, _bsRight, _Int, &(up->x1)) ;
      bsAddData (obj, _bsRight, _Int, &(up->x2)) ;
    }
  bsSave(obj) ;

  obj = bsUpdate(alignment) ;
  if (obj)
    {
      bsAddKey (obj, _Subsequence, cosmid) ;
      i = 1 ;
      bsAddData (obj, _bsRight, _Int, &i) ;
      i = arrayMax(dna) ;
      bsAddData (obj, _bsRight, _Int, &i) ;
      bsAddKey (obj, _Subsequence, subAlignment) ;
      i = 1 ;
      bsAddData (obj, _bsRight, _Int, &i) ;
      i = arrayMax(dna) ;
      bsAddData (obj, _bsRight, _Int, &i) ;
      bsSave (obj) ;
    }
}

/*********************************************************************/
/*********************************************************************/
/* option 1: remove from here, 2 ignore, 3 trash */
void cDNARemoveCloneFromGene (KEY clone, KEY gene)
{
  KEY seq = 0 ;
  OBJ obj = 0 ;
  int x1, x2 ;
  char *cp ;

  if ((obj = bsUpdate(gene)))
    {
      if (bsFindKey (obj, _cDNA_clone, clone))
        bsRemove (obj) ;
      bsGetKey (obj, _Genomic_sequence, &seq) ;
      x1 = x2 = 0 ;
      if (bsGetData (obj, _Covers, _Int, &x1) &&
          bsGetData (obj, _bsRight, _Text, &cp) &&
          bsGetData (obj, _bsRight, _Int, &x1) &&
          bsGetData (obj, _bsRight, _Int, &x2))
	{}  ;
      bsAddKey (obj, str2tag("Homolog_cDNA"), clone) ; 
      bsSave (obj) ;
      cDNARealignGene (gene, 0, 0, 0, 0, 2, 0, 0) ; /*nofuse, locally, does the clean up */

      if (seq && (obj = bsUpdate (seq)))
        {
          bsAddKey (obj, _Discarded_cDNA, clone) ;
          if (x2)
            { 
              if (x1 > x2) { int tmp = x1 ; x1 = x2 ; x2 = tmp ; }
              bsAddData (obj, _bsRight, _Int, &x1) ;
              bsAddData (obj, _bsRight, _Int, &x2) ;
            }
          bsSave (obj) ;
        }
   }
}

/*********************************************************************/
/*********************************************************************/
 
static void fixEstIntronsNew (KEY cosmid, KEY est, 
                          Array dnaD, Array dnaR, 
                          Array oneHits, Array newHits, BOOL doIt)
{
  int i, j, iii, max = arrayMax(dnaD), maxCdna, dx ;
  int x1, x2, y1, y2, a1, a2, b1, b2, aBest, xBest, bBest, yBest, nn, nBest ;
  int maxjump = 10 ;
  HIT *up, *vp ; 
  A_ERR *eUp, *eVp, *eVpBest ;
  Array cDna = 0, cDnaR = 0, errUp = 0, errVp = 0 ;
  BOOL bestUp, bestVp ;
  
  for (iii = 0 ; iii + 1 < arrayMax(oneHits) ; iii++)
    {
      up = arrp(oneHits, iii, HIT) ;
      vp = arrp(oneHits, iii + 1, HIT) ;

      if (up->x1 == up->clipTop)
	{
	  if (!cDna)
	    cDna = getSeqDna (est) ;
	  a1 = up->a1 - 1 ;  x1 = up->x1 - 1 ; x2 = up->x2 ;
	  if (0 && (bestVp = getBestVpPosition (cDna, dnaD, &x1, x2, &a1, TRUE, FALSE)))
	    {
	      nn = x1 - up->clipTop ;
	      if (nn < 5)
		up->a1 = a1 - nn ;
	    }
	}

      /* nothing to do */
      if (up->x2 <= vp->x1 - 1)
        continue ;
      
      if (!cDna)
        cDna = getSeqDna (est) ;
      /* position a2 x2 on last exact base */    
      a1 = up->a1 - 1 ; a2 = up->a2 ; x1 = up->x1 - 1 ; x2 = up->x2 ;
      bestUp = getBestUpPosition (cDna, dnaD, x1, &x2, a1, &a2, TRUE, FALSE) ;
      b1 = vp->a1 - 1 ; b2 = vp->a2 ; y1 = vp->x1 - 1 ; y2 = vp->x2 ;
      bestVp = getBestVpPosition (cDna, dnaD, &y1, y2, &b1, TRUE, FALSE) ;
    
      /* perfect solution exists */
      if (bestUp && bestVp && x2 >= y1 - 1)
        {
          dx = x2 - y1 + 1 ;
          up->x2 = x2 - dx ; up->a2 = a2 - dx ;
          vp->x1 = y1 ; vp->a1 = b1 ;

          continue ;
        }
      if (!doIt)    /* on first pass do not do fancy things */
        continue ;

      /* need to minimize the total number of errors */
      /* if possible use best to minimize errMax */
      if (bestUp) { a1 = a2 - 1 ; x1 = x2 - 1 ;}
      else { a1 = up->a1 ; x1 = up->x1 ; }
      a2 = up->a2 ; x2 = up->x2 ;
      if (x1 > x2) 
        {
          x2 = x1 + 1 ; a2 = a1 + 1 ; /* AK027782 [chromo human 19] */
          if (x1 > y2)
            {
              dx = x1 - y2 + 2 ;
              x1 -= dx ; a1 -= dx ;
            }
        }
      if (bestVp) { b2 = b1 + 1 ; y2 = y1 + 1 ; }
      else { b2 = vp->a2 ; y2 = vp->x2 ; }
      b1 = vp->a1 ; y1 = vp->x1 ;
      if (y1 > y2) 
        { 
          y1 = y2 - 1 ; b1 = b2 - 1 ; /* al117481 [chromo human 13] */
          if (y2 < x1)
            {
              dx = x1 - y2 + 2 ;
              y2 += dx ; b2 += dx ;
            }
        }
      if (x2 < y1 + 5)
        { x2 += 5 ; a2 += 5 ; y1 -= 5 ; b1 -= 5 ; }
      /* compute all errors combing towards the intron
         then restore in est direction coordinates
      */
      {
        x1-- ; a1-- ; /* deny Plato */
        errUp = aceDnaTrackErrors (cDna, x1, &x2, dnaD, a1, &a2, 0, errUp, maxjump, -1, FALSE, 0, FALSE) ;
        x1++ ; a1++ ; /* Plato */
        if (arrayMax (errUp))
          for (i = 0, eUp = arrp (errUp, 0, A_ERR) ; i < arrayMax(errUp) ; eUp++, i++)
            { /* Plato */
              eUp->iShort++ ;
              eUp->iLong++ ;
            }
      }
      if (!cDnaR) 
        {
          cDnaR = dnaCopy (cDna) ;
          reverseComplement (cDnaR) ;
        }

      {
        maxCdna = arrayMax (cDna) ;
        y1 = maxCdna - y1 + 1 ; y2 = maxCdna - y2 + 1 ; /* swap */
        b1 = max - b1 + 1 ; b2 = max - b2 + 1 ;
        y2-- ; b2-- ; /* deny Plato */
        errVp = aceDnaTrackErrors (cDnaR, y2, &y1, dnaR, b2, &b1, 0, errVp, maxjump, -1, FALSE, 0, FALSE) ;
        y2++ ; b2++ ; /* Plato */
        y1 = maxCdna - y1 + 1 ; y2 = maxCdna - y2 + 1 ; /* swap */
        b1 = max - b1 + 1 ; b2 = max - b2 + 1 ;        
        if (arrayMax (errVp))
          for (i = 0, eVp = arrp (errVp, 0, A_ERR) ; i < arrayMax(errVp) ; eVp++, i++)
            { /* Plato + swap */
              eVp->iShort = maxCdna - (eVp->iShort + 1) + 1 ;
              eVp->iLong = max - (eVp->iLong + 1) + 1 ;
            }
      }
      /* global minimisation */
      /* for convenience, add an error at the end of both err sets */
      
      eUp = arrayMax (errUp) ? arrayp (errUp, arrayMax (errUp) - 1, A_ERR) : 0 ;
      if (! eUp || eUp->iShort < x2 + 1)
        {
          eUp = arrayp (errUp, arrayMax (errUp), A_ERR) ;
          eUp->iShort = x2 + 1 ; eUp->iLong = a2 + 1;
        }
      /* insure to stop at y2 */
      for (i = 0, eUp = arrp (errUp, 0, A_ERR) ; i < arrayMax(errUp) ; eUp++, i++)
        {
          if (eUp->iShort >= y2)
            {
              dx = eUp->iShort - y2 ;
              eUp->iLong -= dx ;
              eUp->iShort -= dx ;
              arrayMax (errUp) = i + 1 ;
              break ;
            }
        }

      eVp = arrayMax (errVp) ? arrayp (errVp, 0, A_ERR) : 0 ;
      if (! eVp || eVp->iShort < y2)
        { 
          for (j = arrayMax(errVp), eVp = arrayp (errVp, j, A_ERR) ; j >= 1 ; eVp--, j--)
            *eVp = *(eVp - 1) ;
          eVp->iShort = y2 ; eVp->iLong = b2 ;
        }

      /* initialise at a1 */
      nBest = arrayMax (errVp) ;
      xBest = x1 ;
      aBest = a1 ;
      for (i = 0 ; /* initialise i(eUp) here */
           x1 < y1 && i < arrayMax(errUp) ; i++)
        { 
          eUp = arrp (errUp, i, A_ERR) ; 
          if (eUp->iShort >= y1) break ; 
          nBest-- ;
          xBest = eUp->iShort ;
          aBest = eUp->iLong ;
        }
      /* position b1 */
      eVp = 0 ;
      for (j = arrayMax(errVp) - 1 ; j >= 0 ; j--)
        {
           eVp = arrp (errVp, j, A_ERR) ;
           if (eVp->iShort <= xBest)
            nBest-- ;
          else
            break ;
        }
      yBest = eVp->iShort ; bBest = eVp->iLong ; eVpBest = eVp ;
     
      /* iterate on eUp */
      for (nn = nBest; nBest > 0 && i < arrayMax(errUp) ; i++, nn++)
        { 
          eUp = arrp (errUp, i, A_ERR) ;
          for ( ; j >= 0 ; j--)
            {
              eVp = arrp (errVp, j, A_ERR) ;
              if (eVp->iShort < eUp->iShort)
                nn-- ;
              else
                break ;
            }
          if (nn <= nBest) /* <= because initial a1/x1 was arbitrary */
            { 
              nBest = nn ; 
              xBest = eUp->iShort - 1 ; aBest = eUp->iLong - 1 ;
              yBest = eVp->iShort ; bBest = eVp->iLong ; eVpBest = eVp ;
            }
        }
 
      /* finalize */
      up->a2 = aBest ; /* we know these coords match */
      up->x2 = xBest ;
      
      vp->a1 = bBest ; /* we know these coords match */
      vp->x1 = yBest ;

      /* now we eat up the overlap */
      dx = up->x2 - vp->x1 + 1 ;
      vp->x1 += dx ;
      vp->a1 += dx ;
      if (dx < 0) /* we span the last error , compensate */
        switch (eVpBest->type)
          {                    
          case TROU:  vp->a1 -= 1 ; break ;
          case TROU_DOUBLE:  vp->a1 -= 2 ; break ;
          case INSERTION:   vp->a1 += 1 ; break ;
          case INSERTION_DOUBLE:  vp->a1 += 2 ; break ;
          default: break ;
          }
      /* loop on next intron */
    }
  
  /* drop tiny fuzzy introns */
  for (iii = 0 ; iii < arrayMax(oneHits) - 1 ; iii++)
    {
      up = arrp(oneHits, iii, HIT) ;
      vp = arrp(oneHits, iii + 1, HIT) ;
      /* nothing to do on introns >= 9 bp or if there is a double cover */
      if (vp->a1 > up->a2 + 8 ||
          up->x2 > vp->x1 + 8)
        continue ;
      
      /* between 3 and 5 bp count the errors locally */
      if (vp->a1 > up->a2 + 2)
        {
          char *cp, *cq ;
          int nerr = 0 ;

          if (!cDna)
            cDna = getSeqDna (est) ;
          
          cp = arrp (cDna, up->x2 - 1, char) ;
          cq = arrp (dnaD, up->a2 - 1, char) ;
          
          for (i = 0 ; i < 6 ; cp--, cq--, i++)
            if (*cp != *cq) nerr++ ;
          
          cp = arrp (cDna, vp->x1 - 1, char) ;
          cq = arrp (dnaD, vp->a1 - 1, char) ;
          
          for (i = 0 ; i < 6 ; cp++, cq++, i++)
            if (*cp != *cq) nerr++ ;
          
          if (nerr < 2) /* keep the short introns that match perfectly */
            continue ;
        }
      /* drop the others */
      up->a2 = vp->a2 ;
      up->x2 = vp->x2 ;
      for (i = iii+1 ; i < arrayMax(oneHits) - 1  ; i++)
        {
          up = arrp(oneHits, i, HIT) ;
          vp = arrp(oneHits, i + 1, HIT) ;
          *up = *vp ;
        }
      arrayMax(oneHits) -= 1 ;
    }
  arrayDestroy (cDnaR) ;
  arrayDestroy (errUp) ;
  arrayDestroy (errVp) ;

  for (i = 0, j = arrayMax(newHits) ; i < arrayMax (oneHits) ; j++, i++)
    {
      up = arrayp (oneHits, i, HIT) ;
      vp = arrayp (newHits, j, HIT) ;
      *vp = *up ;
    }
    
  return ;
} /* fixEstIntronsNew */

/*********************************************************************/
 
static void complementHits (Array hits, int max, int i1, int i2)
{
  int i ;
  HIT *vp ;
  
  for (i = i1; i < i2 ; i++)
    { 
      vp = arrp(hits, i, HIT) ;
      vp->a1 = max - vp->a1 + 1 ; vp->a2 = max - vp->a2 + 1 ;
    }
} /* complementHits */

/*********************************************************************/
 
static void complementXHits (Array hits, int max, int i1, int i2)
{
  int i ;
  HIT *vp ;
  
  for (i = i1; i < i2 ; i++)
    { 
      vp = arrp(hits, i, HIT) ;
      vp->x1 = max - vp->x1 + 1 ; vp->x2 = max - vp->x2 + 1 ;
    }
} /* complementHits */

/*********************************************************************/
 
static void fixOneIntrons (KEY cosmid, KEY est, Array oneHits, Array hits, Array dnaD, Array dnaR, BOOL doIt)
{
  BOOL isUp ;
  int i1, max = arrayMax(dnaD) ;
  HIT *vp ; 
  Array myD, myR ;

  if (arrayMax(oneHits))
    {
      vp = arrp(oneHits, 0, HIT) ;
      if (vp->a1 > vp->a2)
        {
          isUp = TRUE ;  myD = dnaR ; myR = dnaD ;
          complementHits (oneHits, max, 0, arrayMax(oneHits)) ; 
        }
      else
        {
          isUp = FALSE ; myD = dnaD ; myR = dnaR ;
        }
      i1 = arrayMax(hits) ;
      fixEstIntronsNew (cosmid, est, myD, myR, oneHits, hits, doIt) ;
      if (isUp)
        complementHits (hits, max, i1, arrayMax(hits)) ; 
    }
  return ;
} /* fixOneIntrons */

/*********************************************************************/

static void fixIntrons (KEY cosmid, Array hits, Array dna, Array dnaR, BOOL doIt)
{
  Array oldHits = arrayCopy (hits), oneHits = 0 ;
  KEY est, old ;
  int j, j1 ;
  HIT *up ;


  cDNASwapX (oldHits) ;
  arraySort (oldHits, cDNAOrderByX1) ;
  arrayMax (hits) = 0 ;

  oneHits = arrayCreate (12, HIT) ; j1 = 0 ; old = 0 ;
  for (old = 0, j = 0 ; j < arrayMax(oldHits) ; j ++)
    {
      up = arrp(oldHits, j, HIT) ;
      est = up->est ;
      if (est != old)
        {
          if (old) fixOneIntrons (cosmid, old, oneHits, hits, dna, dnaR, doIt) ;
          oneHits = arrayReCreate (oneHits, 12, HIT) ;  j1 = 0 ; old = est ;
        }
      array (oneHits, j1++, HIT) = *up ;
    }
  if (old) fixOneIntrons (cosmid, old, oneHits, hits, dna, dnaR, doIt) ;

  arrayDestroy (oldHits) ;
  arrayDestroy (oneHits) ;
}

/*********************************************************************/

static BOOL cdnaDropIsLoner (Array hits, int kk, int type)
{
  HIT *up, *vp ;
  int i, a1, a2 ;
  KEY clone, est ;
  
  up = arrp (hits, kk, HIT) ;
  a1 = up->a1 ; a2 = up->a2 ;
  clone = up->cDNA_clone ;
  est = up->est ;
  
  if (!clone || !est)
    return TRUE ;
  
  /* first check i am not overshooting my own other read */
  vp = arrp (hits, 0, HIT) - 1 ; i = arrayMax(hits) ;
  while (vp++, i--)
    if (vp->est != est && vp->cDNA_clone == clone)
      switch (type)
        {
        case -10:
          if (!vp->reverse && vp->a2 > a2 + 30 &&  /* select only 5p reads */
              vp->a2 - vp->a1 > up->a2 - up->a1)
              return TRUE ; /* overshooting */
          break ;
        case -30:
          if (vp->reverse && vp->a2 > a2 + 30 && /* select only 3p reads */
              vp->a2 - vp->a1 > up->a2 - up->a1)
            return TRUE ; /* overshooting */
          break ;
        case 10:
          if (!vp->reverse && vp->a1 > a1 + 30 &&  /* select only 5p reads */
              vp->a2 - vp->a1 > up->a2 - up->a1)
            return TRUE ; /* overshooting */
          break ;
        case 30:
          if (vp->reverse && vp->a1 < a1 - 30 &&  /* select only 3p reads */
              vp->a2 - vp->a1 > up->a2 - up->a1)
            return TRUE ; /* overshooting */
          break ;
        }

  /* now check if loner */
  vp = arrp (hits, 0, HIT) - 1 ; i = arrayMax(hits) ;
  while (vp++, i--)
    if (vp->est != est)
      switch (type)
        {
        case -1: case -10:
          if (!vp->reverse && vp->a1 < a1)  /* select only 5p reads */
            return FALSE ; /* not a loner */
          if (vp->a2 == a2 && vp->a1 <= a1 && vp->a2 - vp->a1 - 2 * vp->nerr > 30)
            return FALSE ; /* not a loner */
          break ;
        case -2:
          if (vp->a2 > a1 && vp->a1 < a2)  /* somebody intersects me */
            return FALSE ; /* not a loner */
          if (vp->a2 == a2 && vp->a1 <= a1 && vp->a2 - vp->a1 - 2 * vp->nerr > 30)
            return FALSE ; /* not a loner */
          break ;
        case -30:
          if (vp->reverse && vp->a1 < a2)  /* select only 3p reads */
            return FALSE ; /* not a loner */
          if (vp->a2 == a2 && vp->a1 <= a1 && vp->a2 - vp->a1 - 2 * vp->nerr > 30)
            return FALSE ; /* not a loner */
          break ;
        case 1: case 10:
          if (!vp->reverse && vp->a1 > a1)  /* select only 5p reads */
            return FALSE ; /* not a loner */
          if (vp->a1 == a1 && vp->a2 >= a2 && vp->a2 - vp->a1 - 2 * vp->nerr > 30)
            return FALSE ; /* not a loner */
          break ;
        case 2:
          if (vp->a2 > a1 && vp->a1 < a2)  /* somebody intersects me */
            return FALSE ; /* not a loner */
          if (vp->a1 == a1 && vp->a2 >= a2 && vp->a2 - vp->a1 - 2 * vp->nerr > 30)
            return FALSE ; /* not a loner */
          break ;
        case 30:
          if (vp->reverse && vp->a2 > a1)  /* select only 3p reads */
            return FALSE ; /* not a loner */
          if (vp->a1 == a1 && vp->a2 >= a2 && vp->a2 - vp->a1 - 2 * vp->nerr > 30)
            return FALSE ; /* not a loner */
          break ;
        }
  return TRUE ;
}
              
/*********************************************************************/

static void dropTerminalExons (Array hits, int *jp, int j1, int j2)
{
  int jj, iu, iv, da, da1, da2, ba, j01 = j1, j02 = j2 ;
  HIT *up = 0, *vp ;
  BOOL isGeneUp = FALSE ;
  
  up = arrp(hits, j1, HIT) ;
  isGeneUp = up->x2 < up->x1 ? TRUE : FALSE ;
  if (up->reverse)
    isGeneUp = ! isGeneUp ;
  for (jj = j1; jj < j2 - 1; jj++)
    {
      iu = jj ; iv = jj + 1 ;
      up = arrp(hits, iu, HIT) ;
      vp = arrp(hits, iv, HIT) ;
      da = vp->a1 - up->a2 ;
      /* ba = nb of base to be more significant than intron * 64 */
      /* june 4 2002, i increase penalty for errors */
      for (ba = 4 ; (1 << (2*ba-6)) < da ; ba++) ;
      da1 = up->a2 - up->a1 - up->nerr - 0.8 * (up->nerr > 2 ? up->nerr : 0) - 0.8 * (up->nerr > 4 ? up->nerr : 0) - ba - 0.3 * vp->nerr ; if (up->x1 == up->clipEnd || up->x1 == vp->clipTop) da1 += 2 ;
      da2 = vp->a2 - vp->a1 - vp->nerr - 0.8 * (vp->nerr > 2 ? vp->nerr : 0) - 0.8 * (vp->nerr > 4 ? vp->nerr : 0) - ba - 0.3 * up->nerr ; if (vp->x2 == vp->clipEnd || vp->x2 == vp->clipTop) da2 += 2 ;
      if (!up->est || !vp->est) continue ;
      if (iu == j1 || (jj>j1 && !((up-1)->est)))   /* do not trust very long terminal introns in a given group */
        { /* june 4, 2002, i add +10 to calling dropIsLoner 
             and pass da1-z,z as old limit, to kill overshooting */
          if (! isGeneUp && up->reverse && da1 < 8 && /* short tail of a 3p read going up */
              (da1 < 2 || (up->nerr && cdnaDropIsLoner (hits, iu, -10))))  /* nothing above from same clone */
            { up->est = 0 ;  j1++; jj-- ; continue ; } /* drop it */
          if (! isGeneUp && !up->reverse && da1 < 6 && /* first exon of a 5p read going down */
              (da1 < 0 || (up->nerr && cdnaDropIsLoner (hits, iu, -1))))  /* nothing above from same clone */
            { up->est = 0 ;  j1++; jj-- ; continue ; } /* drop it */
          if (isGeneUp && up->reverse && da1 < 10 && /* first exon of a 3p read going down */
              (da1 < 6 || (up->nerr && cdnaDropIsLoner (hits, iu, -2))))  /* nothing above from same clone */
            { up->est = 0 ;  j1++; jj-- ; continue ; } /* drop it */
          if (isGeneUp && !up->reverse && da1 < 10 && /* short tail of a 5p read going up */
              (da1 < 2 || (up->nerr && cdnaDropIsLoner (hits, iu, -30))))  /* nothing above from same clone */
            { up->est = 0 ; j1++; jj-- ; continue ; } /* drop it */
        }
      if (iv == j2 - 1)   /* do not trust very long terminal introns */
        {
          if (! isGeneUp && !vp->reverse && da2 < 10 && /* short tail of a 5p read going down */
              (da2 < 2 || (vp->nerr && cdnaDropIsLoner (hits, iv, 30))))  /* nothing above from same clone */
            { vp->est = 0 ; if (jj > j1 && (vp-1)->est) jj -= 2 ; j2-- ; continue ; } /* drop it */
          if (! isGeneUp && vp->reverse && da2 < 8 && /* first exon of a 3p read going up */
              (da2 < 4 || (vp->nerr && cdnaDropIsLoner (hits, iv, 2))))  /* nothing above from same clone */
            { vp->est = 0 ; if (jj > j1 && (vp-1)->est) jj -= 2 ; j2-- ; continue ; } /* drop it */
          if (isGeneUp && !vp->reverse && da2 < 6 && /* first exon of a 5p read going up */
              (da2 < 0 || (vp->nerr && cdnaDropIsLoner (hits, iv, 1))))  /* nothing above from same clone */
            { vp->est = 0 ; if (jj > j1 && (vp-1)->est) jj -= 2 ; j2-- ; continue ; } /* drop it */
          if (isGeneUp && vp->reverse && da2 < 10 && /* short tail of a 3p read going down */
              (da2 < 2 || (vp->nerr && cdnaDropIsLoner (hits, iv, 10))))  /* nothing above from same clone */
            { vp->est = 0 ; if (jj > j1 && (vp-1)->est) jj -= 2 ; j2-- ; continue ; } /* drop it */
        } 
    }

  for (jj = j01; jj < j02 ; jj++)  /* register happy few */
    { 
      up = arrp(hits, jj, HIT) ;
      if (up->est)
        {
          if (jj != *jp)
            arr (hits, (*jp), HIT) = arr(hits, jj, HIT) ;  
          (*jp)++ ;
        }
    }
}

/*********************************************************************/

static void dropLowEntropyHits (Array dna, Array dnaR, Array hits, int *jp, int j1, int j2)
{
  int i, j, nn, ntotal, x1, x2, na, nt, ng, nc,  nmatch = 0 ;
  Array dnaEst = 0 ;
  HIT *up = 0 ;
  KEY est = 0 ;
  char *cp ;
  double stotal = 0, s = 0, x, y ;

  x1 = x2 = na = nt = ng = nc = nn =  ntotal = 0 ;
  for (i = j1; i < j2 && ! ntotal ; i++)
    {
      up = arrp(hits, i, HIT) ;
      if (up->est)
	{
	  if (keyFindTag (up->est, _IntMap))
	    goto ok ;
	  ntotal = up->clipEnd - up->clipTop + 1 ;
	}
    }
  
  if (ntotal < 50) /* only keep very best single exon quasi perfect hits */
    {
      if (!keyFindTag (up->est, str2tag("Ref_mRNA")))
        return ;
      /* keep only best hits */
      na = j = 0 ;
      for (i = j1; i < j2 ; i++)
        {
          up = arrp(hits, i, HIT) ;
          nt = up->a2 - up->a1 + 1 - up->nerr ;
          if (na < nt) na = nt ;
        }
      if (na < ntotal - 1 || (ntotal < 25 && na < ntotal))
        na = 1000 ; /* discard all */
      for (i = j1; i < j2  ; i++)
        {
          up = arrp(hits, i, HIT) ;
          nt = up->a2 - up->a1 + 1 - up->nerr ;
          if (na > nt) 
            up->est = 0 ;
        }
    }
      
  x1 = x2 = na = nt = ng = nc = nn = 0 ;
  for (i = j1; i < j2 ; i++)
    {
      up = arrp(hits, i, HIT) ;
      up->zone = 0 ;
      est = up->est ; 
      if (! est)
        continue ;

      na = nt = ng = nc = nn = 0 ;
      if (up->x1 < up->x2)
	{ x1 = up->x1 - 1 ; x2 = up->x2 ; }
      else
	{ x1 = up->x2 - 1 ; x2 = up->x1 ; }
      dnaEst = getSeqDna (up->est) ;
      if (!dnaEst) continue ;
      if (x1 < 0) x1 = 0 ;
      if (x2 > arrayMax(dnaEst))
	x2 = arrayMax(dnaEst) ;
      cp = arrp(dnaEst, x1, char) ;
      for (j = x1 ; j < x2 ; j++, cp++)
        switch (*cp)
          {
          case A_: na++ ; break ;
          case T_: nt++ ; break ;
          case G_: ng++ ; break ;
          case C_: nc++ ; break ;
          default: nn++ ; break ;
          }
      y = x2 - x1 ;
      nmatch += y ;
      s = 0 ;  /* casting necessary */
      x = na ; if (na && y > 0) s += x * log(x/y) ;
      x = nt ; if (nt && y > 0) s += x * log(x/y) ;
      x = ng ; if (ng && y > 0) s += x * log(x/y) ;
      x = nc ; if (nc && y > 0) s += x * log(x/y) ;

      s = - s / log (4.0) ; /* so i count in base equivalent */
      if ((na > nt + ng + nc || nt > na + ng + nc) && up->nerr) 
        s -= 4 ; /* do not trust polyA polyT */
      /* 4a779 (worm) -> ajoute la cond nerr >0 pour retrancher 4 , 1 0ct 2003 */
      if (s - 2 * up->nerr < 5.0)
        up->zone = 1 ; /* drop it */
      if ((na > 2*(nt + ng + nc) && nt + ng + nc < 2 * up->nerr) ||
          (nt > 2*(na + ng + nc) && na + ng + nc < 2 * up->nerr))
        up->zone = 1;
      if ((na > 5*(nt + ng + nc - up->nerr)) ||
          (nt > 5*(na + ng + nc - up->nerr)))
        up->zone = 1;
      if (! up->zone)
        stotal += s ;
      if (ntotal < 50) /* only allow single exon hits */
        break ;
    }
  /* now we cut out all the baddies at both ends */
  for (i = j1; i < j2 ; i++)
    {
      up = arrp(hits, i, HIT) ;
      if (up->zone)
        up->est = 0 ;
      else
        break ;
    }
  for (i = j2 - 1; i >= j1 ; i--)
    {
      up = arrp(hits, i, HIT) ;
      if (up->zone)
        up->est = 0 ;
      else
        break ;
    }

  if (nmatch < 10) return ; /* reject */

  /* filter on entropy evaluation */
  /* ntotal = est length, nerr, nb of errors, nmatch, lenght of match */
  if (ntotal < 25)
    {
      if (stotal < 15)
        return ;
    }
  else if (ntotal < 50)
    {
      if (stotal < 20)
        return ;
    }
  else if (stotal < 30) 
    return ;
  /* printf("entropy %s s=%d nn=%d = %d + %d + %d + %d\n", 
   *  name(est), (int) s, nmatch, na, nt, ng, nc) ;
   */
  /* register happy few */
 ok:
  for (i = j1; i < j2 ; i++)
    { 
      up = arrp(hits, i, HIT) ;
      up->zone = 0 ; /* reestablish */
      if (up->est)
        {
          if (*jp != i)
            arr (hits, (*jp), HIT) = arr(hits, i, HIT) ;  
          (*jp)++  ;
        }
    }
}

/*********************************************************************/

static void dropSmallGlobalHits (Array hits, int *jp, int j1, int j2, int mini)
{
  int i, n, nn, nt, gap = 0, oldx1 = 0 ;
  HIT *up = 0 ;
  KEY est = 0 ;
  int v1 = 0 , v2 = 0 ;

  nn = 0 ;
  for (i = j1; i < j2 ; i++)
    {
      up = arrp(hits, i, HIT) ;
      if (!est) { oldx1 = up->x1 ; }
      est = up->est ; 
      v1 = up->clipTop ; v2 = up->clipEnd ;
      n = up->reverse ? up->x1 - up->x2  : up->x2 - up->x1 ;
      if (n < 0) n-- ; else n++ ;
      nn += n ;
    }
  /* note that we do an algebraic sum so zitterbewegung cancels out */
  n = nn < 0 ? -nn : nn ;
  /* counting negative gap penalty should eliminate random
     long full genes that match on repeats 
  */
  if (up->x2 > oldx1) gap = up->x2 - oldx1 - n ;
  else gap = oldx1 - up->x2 - n ;
  n -= gap/2 ;
  nt = v2 - v1 + 1 ;

  if (n >= 150)      
    {  if ( 10 * n < nt)   return ;}  /* 10% coverage at least */
  else if (n >= 50)   
    {  if ( 10 * n < 5 * nt)   return ; }  /* 50% coverage at least */
  else
    {
      if (10 * n < 7 * nt)    /* 70% coverage at least */
        return ;
      if (nt < 50 && n < nt - 1)
        return ;
    }

  nn = nn > 0 ? 1 : -1 ;
  for (i = j1; i < j2 ; i++)
    { 
      up = arrp(hits, i, HIT) ;
      n = up->reverse ? up->x1 - up->x2 + 1 : up->x2 - up->x1 + 1 ;
      if (n * nn < mini) /* drop backwards or mini hits */
        continue ;
      if (*jp != i)
        arr (hits, (*jp), HIT) = arr(hits, i, HIT) ;  
      (*jp)++  ;
    }
}

/*********************************************************************/
/* discard baq quality reads but not if 5' read also matches */
static void dropBadQualityGlobalHits (Array dna, Array dnaR, Array hits, int *jp, int j1, int j2)
{
  int i=0, maxExact = 0, nn, nt, nerr = 0, a1, a2, x1min, x1, x2, nmatch = 0, aa1, aa2, xx1, xx2 ;
  HIT *up = 0 ;
  KEY est = 0, clone = 0 ;
  Array err1 = 0, dnaLong = 0, dnaEst = 0 ;
  int ok = 0, quality = 0 ;
  BOOL isMrna = 0 ;
  KEY  _Transfered = str2tag ("Transfered") ;

  a1 = a2 = x1 = x2 = nn = 0 ;
  x1min = 10000 ;

  up = arrp(hits, j1, HIT) ;
  clone =   up->cDNA_clone ;
  est = up->est ;
  if (keyFindTag (est, _Transfered))
    goto laba ;
  /* just check previous clones */
  nt = 0 ;
  a1 = a2 = up->a1 ; x1 = x2 = up->x2 ; /* good initialisation */
  for (i = 0 ; i < arrayMax(hits) ; i++)
    {
      up = arrp(hits, i, HIT) ;
      if (clone == up->cDNA_clone)
        { 
          if (up->x1 < x1min)
            x1min = up->x1 ;
          if (up->x2 < x1min)
            x1min = up->x2 ;
          if (a1 > up->a1) 
            a1 = up->a1 ;
          if (a1 > up->a2) 
            a1 = up->a2 ;
          if (a2 < up->a1)
            a2 = up->a1 ;
          if (a2 < up->a2)
            a2 = up->a2 ;
          if (est != up->est)
            { ok = 2 ;  }
	  else
	    maxExact += up->maxExact ;
          nt += up->a2 - up->a1 ;
        }
    }
  if (maxExact > 50)
    ok = 2 ;
  else if (ok == 2)  /* verify length and x1min */
    {
      if (nt < 300 || x1min > 200 || a2 - a1 > 10000)
        ok = 0 ;
    }


  for (i = j1; ok < 2 && i < j2 ; i++)
    {
      up = arrp(hits, i, HIT) ;
      est = up->est ; 
      nt = up->clipEnd - up->clipTop + 1 ;
      nmatch = nerr = maxExact = 0 ;

      if (up->x1 <= up->x2) 
        { 
          a1 = aa1 = up->a1 - 1 ; a2 = aa2 = up->a2 ; 
          x1 = xx1 = up->x1 - 1 ; x2 = xx2 = up->x2 ;
          dnaLong = dna ; 
        } 
      else if (up->x1 > up->x2) /* reject length zero */
        { 
          a1 = aa1 = (long)arrayMax(dna) - up->a2 ; 
          a2 = aa2 = (long)arrayMax(dna) - up->a1 + 1 ; 
          x1 = xx1 = up->x2 - 1 ; x2 = xx2 = up->x1 ;
          dnaLong = dnaR ; 
        }
      
      if ((dnaEst = getSeqDna (est)))
        {
          int maxjump = 4, k = 0, mx, nerr1 ;
          A_ERR *ep ;
          int ok1 = 0 ;
	  maxExact = 0 ;
          for (k = 0 ; k < 3 ; k++) /*do not loop forever ! */
            {
	      mx = 0 ;
              x1 = xx1 + 10 * k ; a1 = aa1 + 10 * k ;
              x2 = xx2 ; a2 = aa2 ; 
              err1 = aceDnaTrackErrors (dnaEst, x1, &x2, 
                                        dnaLong, a1, &a2, 0, err1, maxjump, 0, FALSE, &mx, FALSE) ;
              if (x2 > xx2 - maxjump || a2 > aa2 - maxjump)
                break ;
            }
          errRemoveAmbigue (err1) ;
          nmatch += a2 - a1 ;
          nerr1 = arrayMax(err1) ;
          nerr += nerr1 ;
	  maxExact += mx ;
          if (x1 < x1min) x1min = x1 ;

	  if (maxExact > 50) ok = 2 ;
          if (ok < 2 && a2 - a1 > 50 && nerr1 * 20 < a2 - a1) ok = 2 ;
          if (ok < 2 && a2 - a1 > 50 && nerr1 * 10 < a2 - a1) ok = 1 ;
          
          if (ok<2 && nerr1 > 0)
            {
              ep = arrp (err1, 0, A_ERR) ;
              if (ok < 2 && ep->iShort - x1 > 50)
                { if (!ok1) ok++ ; ok1 = 1 ; }
              if (ep->iShort - x1 > 100)
                { ok = 2 ; ok1 = 1 ; }
            }
          
          if (ok < 2 && nerr1 > 1)
            {  
              for (k = 1, ep = arrp (err1, k, A_ERR) ;
                   ok < 2 && k < arrayMax(err1) ; k++, ep++)
                {
                  if (ep->iShort < 350 && ep->iShort > (ep-1)->iShort + 100)
                    { ok = 2 ; ok1 = 1 ; }
                  else if (ep->iShort < 350 && ep->iShort > (ep-1)->iShort + 50)
                    { if (!ok1) ok++ ; ok1 = 1 ; }
                }
            }
          if (ok < 2 && arrayMax(err1) > 2)
            {
               for (k = 2, ep = arrp (err1, k, A_ERR) ;
                    ok < 2 && k < arrayMax(err1) ; k++, ep++)
                 {
                   if (ep->iShort < 350 && ep->iShort > (ep-2)->iShort + 50)
                     { if (!ok1) ok++ ; ok1 = 1 ; }
                   if (ep->iShort < 350 && ep->iShort > (ep-2)->iShort + 100)
                     { ok = 2 ; ok1 = 1 ; }
                 }
            }
          
          /*
            printf("%s a1=%d a2=%d up->a2=%d x1=%d x2=%d up->x2=%d nmatch=%d nerr=%d   ok=%d\n",
            name(est),  a1,   a2,   up->a2,   x1,   x2,   up->x2,   nmatch,   nerr,     ok) ;
            */
        }
      dnaEst = 0 ;  /* do not destroy, got by getSeqDna */
    }
  arrayDestroy (err1) ;

  /* filter on quality evaluation */
  /* nt = est length, nerr, nb of errors, nmatch, lenght of match */
  switch (ok)
    {
    case 0: /* discard */
      if (20 * nerr < nmatch && nmatch > 450)       /* very long match, accept */
        break ;
      if (14 * nerr > nmatch)  /* 14 per 100 error rate everywhere , no ok zone reject */
        return ;
      /* fall thru */
    case 1:   /* either globally or locally better than 10 % */
      if (20 * nerr <  nmatch)  /* globally less than 5% error, accept */
        break ;
      if (nmatch > 450)       /* very long match, accept */
        break ;
      if (14 * nerr > nmatch) /* 14 per 100 global error rate, one ok zone  , reject */
        return ;
      if (x1min < 250 && nmatch > 150 &&
          (keyFindTag (up->est, _CTF_File) || keyFindTag (up->est, _SCF_File)))
        break ; /* could be salvaged by editing */
      else      /* hopeless */
        return ; 
    case 2: /* good read,  locally better than 5  % */
      break ;
    }

  if (ok < 2)
    { 
      OBJ Est = 0 ;
      BSunit *uu ;
      Array units = 0 ;
      int iUnit ;

      up = arrp(hits, j1, HIT) ;
      isMrna = keyFindTag (up->est, _mRNA) ;
      quality = mrnaQualityEvaluate (up->clipEnd - up->clipTop + 1, nmatch, nerr, isMrna, 0, 0) ;
      
      if (quality >= 11) /* forget disgusting hits */
        return ;
      
      if ((Est = bsCreate (est)))
        {
          units = arrayCreate (64, BSunit) ;
          if (bsGetArray (Est, _From_gene, units, 6))
            {
              for (iUnit = 0 ; quality < 11 && iUnit < arrayMax(units) ; iUnit += 6)
                {
                  uu = arrp (units, iUnit, BSunit) ;
                  if (uu[5].i > 0 && uu[5].i < quality - 3) /* 4 units of quality, but only one est */
                    quality = 11;
                }
            }
          bsDestroy (Est) ; 
          arrayDestroy (units) ;
          if (quality >= 11) /* forget disgusting hits */
            return ;
        }
    }  
  /* register happy few */
 laba:
  for (i = j1; i < j2 ; i++)
    { 
      up = arrp(hits, i, HIT) ;
      if (*jp != i)
        arr (hits, (*jp), HIT) = arr(hits, i, HIT) ;  
      (*jp)++ ; 
    }
  return ; 
} /* dropBadQualityGlobalHits */

/*********************************************************************/
/* discard bad quality composite reads: i.e. pseudo-ESTs supported by deep seq copied from the genome */
static void dropBadCompositeGlobalHits (Array dna, Array dnaR, Array hits, int *jp, int j1, int j2)
{
  int i=0, nerr, nmatch = 0 ;
  HIT *up = 0 ;
  BOOL ok = TRUE ;

  up = arrp(hits, j1, HIT) ;
  nmatch = nerr = 0 ;
  if (keyFindTag (up->est, _Composite))
    {
      for (i = j1; ok && i < j2 ; i++)
	{
	  up = arrp(hits, i, HIT) ;
	  if (up->x1 <= up->x2) 
	    nmatch += up->x2 - up->x1 + 1 ;
	  else
	    nmatch += up->x1 - up->x2 + 1 ;
	  nerr += up->nerr ;
	  if (nerr > 5)
	    ok = FALSE ;
	}
      if (100 * nerr > 2 * nmatch)
	ok = FALSE ;
    }
  
  /* register happy few */
  for (i = j1; ok && i < j2 ; i++)
    { 
      up = arrp(hits, i, HIT) ;
      if (*jp != i)
        arr (hits, (*jp), HIT) = arr(hits, i, HIT) ;  
      (*jp)++ ; 
    }
  return ; 
} /* dropBadCompositeGlobalHits  */

/*********************************************************************/

static void dropBadHits (Array dna, Array dnaR, Array hits, int type)
{
  KEY est, old ;
  int j, j0, j1;
  HIT *up ;

  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;

  j0 = j1 = 0 ;
  for (old = 0, j = 0 ; j < arrayMax(hits) ; j ++)
    {
      up = arrp(hits, j, HIT) ;
      est = up->est ;
      if (est != old)
        {
          if (j > j1)
            {
              if (type == 2) dropTerminalExons (hits, &j0, j1, j) ;
              /* else if (type == 3) dropNonColinearHits (hits, &j0, j1, j) ; */
              else if (type == 5) dropSmallGlobalHits (hits, &j0, j1, j, 5) ;
              else if (type == 6) dropLowEntropyHits (dna, dnaR, hits, &j0, j1, j) ;
              else if (type == 7) dropBadQualityGlobalHits (dna, dnaR, hits, &j0, j1, j) ;
              else if (type == 8) dropBadCompositeGlobalHits (dna, dnaR, hits, &j0, j1, j) ;
            }
          j1 = j ;
        }
      old = up->est ;
    }
  if (j > j1)
    {
      if (type == 2) dropTerminalExons (hits, &j0, j1, j) ;
      /* else if (type == 3) dropNonColinearHits (hits, &j0, j1, j) ; */
      else if (type == 5) dropSmallGlobalHits (hits, &j0, j1, j, 5) ; 
      else if (type == 6) dropLowEntropyHits (dna, dnaR, hits, &j0, j1, j) ;
      else if (type == 7) dropBadQualityGlobalHits (dna, dnaR, hits, &j0, j1, j) ; 
      else if (type == 8) dropBadCompositeGlobalHits (dna, dnaR, hits, &j0, j1, j) ;
    }
  arrayMax (hits) = j0 ; 
  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;
  arrayCompress (hits) ; 
}

/*********************************************************************/

static void dropOvernumerousHits (Array hits, int hMax)
{
  int j, j1, dx ;
  BOOL ok = FALSE ;
  KEY est = 0, _is_mrna ;
  HIT *up ;
  
  if (arrayMax(hits) < hMax)
    return ;
    
  _is_mrna = str2tag ("Is_mRNA") ;
  for (j = 0, j1 = 0 ; j < arrayMax(hits) ; j++)
    {
      up = arrp(hits, j, HIT) ;
      if (est != up->est)
        {
          est = up->est ;
          ok = keyFindTag (est, _is_mrna) || keyFindTag (est, _Ref_mRNA) ||  keyFindTag (est, _Ref_Seq) ||  keyFindTag (est, _IntMap) ;
	  if (!strncmp(name(est),"U454",4))
	    ok = TRUE ;
        }
      if (! ok)
	{
	  dx =  up->a2 > up->a1 ? up->a2 - up->a1 : up->a1 - up->a2 ;
	  if ((50 * up->nerr < dx && dx > 28) || (20 * up->nerr < dx && dx > 100))
	    ok = TRUE ;  /* 2013_09_21 28 (was 30) because of XK uses 30bp as foot */
	  if (!ok)
	    continue ;
	}
      if (j != j1)
        arr (hits, j1, HIT) = arr(hits, j, HIT) ;   
      j1++ ;
    } 
  messout ("dropOvernumerousHits dropped %d/%d hits", (long)arrayMax (hits) - j1, arrayMax (hits)) ;
  arrayMax (hits) = j1 ;
}

/*********************************************************************/

static void dropOutOfZoneHits (Array hits, int z1, int z2)
{
  int j, j1 ;
  HIT *up ;

  if (z1 >= z2) return ;

  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;

  for (j = 0, j1 = 0 ; j < arrayMax(hits) ; j++)
    {
      up = arrp(hits, j, HIT) ;
      if (up->a2 < z1 || up->a1 > z2)
        continue ;
      if (j != j1)
        arr (hits, j1, HIT) = arr(hits, j, HIT) ;   
      j1++ ;
    } 
  arrayMax (hits) = j1 ;
} /* dropOutOfZoneHits */

/*********************************************************************/
/* compare premapped tags to the IntMap coordinates of the cosmid */
static void dropOutOfZonePreMappedHits (Array hits, 
					KEYSET estMaps, KEYSET estMap1, KEYSET estMap2, 
					KEY chrom, int c1, int c2,
					int direction
					)
{
  int j, j1, e1 = 0, e2 = 0 ;
  HIT *up ;
  KEY oldEst = 0, eChrom = 0 ;

  if (chrom == 0 || !estMaps) return ;
  
  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;
  
  for (j = 0, j1 = 0 ; j < arrayMax(hits) ; j++)
    {
      up = arrp(hits, j, HIT) ;
      if (oldEst != KEYKEY(up->est))
	{
	  oldEst = KEYKEY(up->est) ;
	  eChrom = keySet (estMaps, oldEst) ;
	  e1 = keySet (estMap1, oldEst) ;
	  e2 = keySet (estMap2, oldEst) ;
	}
      if (eChrom)
	{
	  if (eChrom != chrom)
	    continue ;
	  if (1)
	    {
	      if (direction >= 2) /* forward */
		{
		  if (
		      (! up->reverse && e1 > e2) ||
		      (up->reverse && e1 < e2)
		      )
		    continue ;
		}
	      else /* reverse */
		{
		  if (
		      (up->reverse && e1 > e2) ||
		      (! up->reverse && e1 < e2)
		      )
		    continue ;
		}
	    }
	  if (
	      (e1 < e2 && (c1 + up->a2 < e1 - 1000 || c1 + up->a1 > e2 + 1000)) ||
	      ( e1 > e2 && (c1 + up->a2 > e1 + 1000 || c1 + up->a1 < e2 - 1000))
	      )
	    continue ;
	}
      if (j != j1)
        arr (hits, j1, HIT) = arr(hits, j, HIT) ;   
      j1++ ;
    } 
  arrayMax (hits) = j1 ;
} /* dropOutOfZonePreMappedHits */

/*********************************************************************/

static void dropOneBackToBackPair (Array hits, BOOL isUp, int *jp, int j1, int j2)
{
  int ii, x1, x2, y1, y2, a1, a2, b1, b2 ;
  HIT *up ;
  KEY est = 0 ;
  int  s1, s2, t1, t2, top1, top2 ;
  BOOL u1, u2 ;

  x1 = x2 = y1 = y2 = 0 ;
  a1 = a2 = b1 = b2 = 0 ;
  s1 = s2 = t1 = t2 = 0 ;
  top1 = top2 = 0 ;
  for (ii = j1; ii < j2 ; ii++)
    {
      up = arrp(hits, ii, HIT) ;
      if (est != up->est) /* all exons allready colinear, checkone */
        {
          if (up->reverse)
            { 
              b1 = up->a1 ; b2 = up->a2 ; y1 = up->x1 ; y2 = up->x2 ; s2 = up->x2 - up->x1 ;
              if (up->x1 < up->x2) b1 -= up->x1 ;
              else b2 += up->x2 ;
               /* if (s2 * t2 < 0) { up->est = 0 ; continue ; } */
              t2 = s2 ; top2 = up->clipTop ;
            }
          else
            { 
              a1 = up->a1 ; a2 = up->a2 ; x1 = up->x1 ; x2 = up->x2 ; s1 = up->x2 - up->x1 ;
              if (up->x1 < up->x2) a1 -= up->x1 ;
              else a2 += up->x2 ;
              /* if (s1 * t1 < 0)  { up->est = 0 ; continue ; } */
              t1 = s1 ; top1 = up->clipTop ;
            }
        }
      else
        {
          if (up->reverse)
            { 
              if (up->x1 < up->x2)
                { if (b2 < up->a2) b2 = up->a2 ; }
              else
                { if (b2 < up->a2 + up->x2) b2 = up->a2 + up->x2 ; }
            }
          else
            { 
              if (up->x1 < up->x2)
                { if (a2 < up->a2) a2 = up->a2 ; }
              else
                { if (a2 < up->a2 + up->x2) a2 = up->a2 + up->x2 ;}
            }
          if (up->reverse)
            {
              s2 += up->x2 - up->x1 ;
              if (y1 < y2)
                {
                  if (y1 > up->x1) y1 = up->x1 ;
                  if (y2 < up->x2) y2 = up->x2 ;
                }
              else
                {
                  if (y1 < up->x1) y1 = up->x1 ;
                  if (y2 > up->x2) y2 = up->x2 ;
                }
            }
          else
            {
              s1 += up->x2 - up->x1 ;
              if (x1 < x2)
                {
                  if (x1 > up->x1) x1 = up->x1 ;
                  if (x2 < up->x2) x2 = up->x2 ;
                }
              else
                {
                  if (x1 < up->x1) x1 = up->x1 ;
                  if (x2 > up->x2) x2 = up->x2 ;
                }
            }
        }
      est = up->est ;
    }
  
  u1 = u2 = TRUE ;

#ifdef JUNK
  This would be the global orientation of the cosmid, no good
  if (s1 > 0 && isUp) u1 = FALSE ;
  if (s2 < 0 && isUp) u2 = FALSE ;
  if (s1 < 0 && !isUp) u1 = FALSE ;
  if (s2 > 0 && !isUp) u2 = FALSE ;
#endif

  /*   s1 = x2 - x1 ;  s2 = y2 - y1 ; */
  if ((s1 * s2 > 0)  || /* case parallel */
      (
       (s1 * s2 < 0) &&    /* antiparallel */
       ((s1 > 0 && a1 > b2) ||
        (s1 < 0 && a2 < b1))
       )
      )
    {
      if (s1 * s1 > s2 * s2) 
        u2 = FALSE ;
      else
        u1 = FALSE ;
    }
  
  
  /* register happy few  */
  for (ii = j1; ii < j2 ; ii++)
    {
      up = arrp(hits, ii, HIT) ; 
      if (!up->est) continue ;
      if (up->reverse && !u2)
        continue ;
      if (!up->reverse && !u1)
        continue ;
      if (u1 && u2 && s1 && s2)  /* drop overshooting */
        {
          /* but do not drop pieces longer than the opposite strand */
          if (!up->reverse && s1 > 0 && up->a1 > b2)
            {
              if (x2 - up->x1 > - 1.7 * s2 || x2 - up->x1 < y2 - top2)
                s1 = 0 ; /*do not drop this exon or any further exon of this est */
              else
                continue ;  /* drop this and obviously folling ones */
            }
          if (!up->reverse && s1 < 0 && up->a2 < b1)
            {
              if (x1 - up->x2 > 1.7 * s2 || x1 - up->x2 < y1 - top2)
                s1 = 0 ; /*do not drop this exon or any further exon of this est */
              else
                continue ;  /* drop this and obviously folling ones */
            }
          if (up->reverse && s2 < 0 && up->a2 < a1)
            { 
              if (y1 - up->x2 > 1.7 * s1 || y1 - up->x2 < x1 - top1)
                s2 = 0 ; /*do not drop this exon or any further exon of this est */
              else
                continue ;  /* drop this and obviously folling ones */
            }
          if (up->reverse && s2 > 0 && up->a1 > a2)
            {
              if (y2 - up->x1 > - 1.7 * s1 || y2 - up->x1 < x2 - top1)
                s2 = 0 ; /*do not drop this exon or any further exon of this est */
              else
                continue ;  /* drop this and obviously folling ones */
            }
        }
#ifdef JUNK
      if (u1 * u2 * s1 * s2 &&  /* drop overshooting */
          (  /* but do not drop pieces longer than the opposite strand */
           (!up->reverse && s1 > 0 && up->a1 > b2 && x2 - up->x1 < - 1.7 * s2) ||
           ( up->reverse && s2 < 0 && up->a2 < a1 && y1 - up->x2 < 1.7 * s1) ||
           (!up->reverse && s1 < 0 && up->a2 < b1 && x1 - up->x2 < 1.7 * s2) ||
           ( up->reverse && s2 > 0 && up->a1 > a2 && y2 - up->x1 < - 1.7 * s1))
          )
        continue ;
#endif
      if (ii != *jp)
        {
          arr (hits, *jp, HIT) = arr(hits,ii, HIT) ;  
        }
      (*jp)++ ;
    }  
}

/*********************************************************************/
/* all exons of a given est allready colinear, just add up */
static BOOL dropBackGetOrientation (Array hits)
{
  int ii = arrayMax(hits), dx = 0 ;
  HIT *up = arrp(hits, 0, HIT) - 1 ;
         
  while (up++, ii--)
    dx +=  up->reverse ? up->x1 - up->x2 : up->x2 - up->x1 ;
   
  return dx < 0 ? TRUE : FALSE ;
}
/*********************************************************************/

static void dropBackToBackPairs (Array hits)
{
  KEY old, clone ;
  int j, j0, j1;
  HIT *up ;
  BOOL isUp = FALSE ;

  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;

  isUp = dropBackGetOrientation (hits) ;
  j0 = j1 = 0 ;
  for (old = 0, j = 0 ; j < arrayMax(hits) ; j ++)
    {
      up = arrp(hits, j, HIT) ;
      clone = up->cDNA_clone ;
      if (clone != old)
        {
          if (j > j1)
            dropOneBackToBackPair (hits, isUp, &j0, j1, j) ;
          j1 = j ;
        }
      old = up->cDNA_clone ;
    }
  if (j > j1)
    dropOneBackToBackPair (hits, isUp, &j0, j1, j) ;

  arrayMax (hits) = j0 ;
  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;
  arrayCompress (hits) ;
}

/*********************************************************************/

static BOOL dropOneHugeIntronIsGtag (Array hits, Array dna, Array dnaR, int ii, int x1, int x2)
{
  HIT *up, *vp ;
  char *cp, *cq ;

  up = arrp(hits, ii - 1, HIT) ;
  vp = arrp(hits, ii, HIT) ;

  cp = arrp (dna, up->a2, char) ; /* first base in the intron, PLATO */
  cq = arrp (dna, vp->a1 - 3, char) ; /* last base - 1 of the intron, PLATO */

  if (
      (*cp == G_ && (*(cp+1) == T_ || *(cp+1) == T_) && *cq == A_ && *(cq+1) == G_) ||
      (*cp == C_ && *(cp+1) == T_  && (*cq == A_ || *cq == G_) && *(cq+1) == C_)
      )
    return TRUE ;
  return FALSE ;
}

/*********************************************************************/

static BOOL dropOneHugeIntronQuality (Array hits, Array dna, Array dnaR, int j1, int j2, int ii, int x1, int x2, BOOL gt_ag)
{
  BOOL ok = TRUE ;
  HIT *up ;
  int da, dx1, dx2, ne1, ne2 ;
  int j ;

  dx1 = dx2 = ne1 = ne2 = 0 ;
  for (j = ii - 1; j >= j1 ; j--)
    {
      up = arrp(hits, j, HIT) ;
      dx1 += up->x2 - up->x1 ;
      ne1 += up->nerr ;
    }
  for (j = ii; j < j2 ; j++)
    {
      up = arrp(hits, j, HIT) ;
      dx2 += up->x2 - up->x1 ;
      ne2 += up->nerr ;
    }

  da = arr(hits, ii , HIT).a1 - arr(hits, ii - 1, HIT).a2 ;

  if (dx1 < 0) dx1 = - dx1 ;
  if (dx2 < 0) dx2 = - dx2 ;
  if (da < 0) da = - da ;
  /*
    
  if (!gt_ag) { ne1 -= 5 ; ne2 -= 5 ;} 
  i.e. before intron adjustment, -10 not good see ak057234
  it blocks a go thru getNewExon to look for a better solution 
  feb 11, 2004:
  exactly same problem with AK096206 with -5
  i could try to always run fix introns first
  in fast mode 
  see also ak128063
  */
  
  if (ii == j1 + 1)  { ne1 -= 3 ; ne2 -= 3 ;} /* errors above intron should not count */
  if (ne1 < 0) ne1 = 0 ; 
  if (ne2 < 0) ne2 = 0 ;
  ne1 = (100.0 * ne1)/(1 + dx1) ;
  ne2 = (100.0 * ne2)/(1 + dx2) ; /* beware of zero */

  if ( 
      ( da < 50000 &&
        (
         (dx1 < 20 && ne1 > 3) || 
         (dx2 < 20 && ne2 > 3) ||
         (dx1 < 50 && ne1 > 10) || 
         (dx2 < 50 && ne2 > 10)
         )
        ) ||
      ( da >= 50000 &&
        (
         (dx1 < 50 && ne1 > 2) || 
         (dx2 < 50 && ne2 > 2) ||
         (dx1 < 200 && ne1 > 8) || 
         (dx2 < 200 && ne2 > 8) ||
         (dx1 < 500 && ne1 > 12) || 
         (dx2 < 500 && ne2 > 12)
         )
        ) ||
      ( da >= 100000 && !gt_ag &&
        (
         (dx1 < 20) ||
         (dx2 < 20) ||
         (dx1 < 50 && ne1 > 0) || 
         (dx2 < 50 && ne2 > 0) ||
         (dx1 < 200 && ne1 > 2) || 
         (dx2 < 200 && ne2 > 2) ||
         (dx1 < 500 && ne1 > 5) || 
         (dx2 < 500 && ne2 > 5)
         )
        ) ||
      ( da >= 50000 &&
        gt_ag &&
        !dropOneHugeIntronIsGtag (hits, dna, dnaR, ii, x1, x2) && 
        (
         (dx1 < 20) ||
         (dx2 < 20) ||
         (dx1 < 50 && ne1 > 0) || 
         (dx2 < 50 && ne2 > 0) ||
         (dx1 < 200 && ne1 > 2) || 
         (dx2 < 200 && ne2 > 2) ||
         (dx1 < 500 && ne1 > 8) || 
         (dx2 < 500 && ne2 > 8)
         )
        )
      )
    ok = FALSE ;
  
  return ok ;
}

/*********************************************************************/

static void dropOneHugeIntrons (Array hits, Array dna, Array dnaR, 
                                int *jp, int j1, int j2, int dxMax, BOOL gt_ag)
{
  int i, ii, a2, x1, x2 ;
  HIT *up ;

  a2 = -1 ;
  x1 = x2 = -1 ;
  
  x1 =  arr(hits, j1, HIT).x1 ; 
  x2 =  arr(hits, j2 - 1, HIT).x2 ; 
  a2 =  arr(hits, j1, HIT).a2 ;

  for (ii = j1 + 1 ; ii < j2 ; ii++)
    {
      up = arrp(hits, ii, HIT) ;
      if (!up->est) continue ;
      if (a2 >= 0 && 
          (
           (dxMax && up->a1 > a2 + dxMax) ||
           (up->a1 > a2 + 25000 && !dropOneHugeIntronQuality (hits, dna, dnaR, j1, j2, ii, x1, x2, gt_ag))
           )
          )           
        {
          if (((up - 1)->x2 - x1) * (x2 - x1) < (x2 - up->x1) * (x2 - x1))
            {
              for (i = j1 ; i < ii ; i++)  arr(hits, i, HIT).est = 0 ;
            }
          else
            {
              for (i = ii ; i < j2 ; i++)  arr(hits, i, HIT).est = 0 ;
            }
        }
      a2 = up->a2 ;
    }
  
  /* register happy few  */
  for (ii = j1; ii < j2 ; ii++)
    {
      up = arrp(hits, ii, HIT) ; 
      if (!up->est) continue ;
      if (ii != *jp)
        {
          arr (hits, *jp, HIT) = arr(hits,ii, HIT) ;  
        }
      (*jp)++ ;
    }  
}

/*********************************************************************/

static void dropHugeIntrons (Array hits, Array dna, Array dnaR, int dxMax, BOOL gt_ag)
{
  KEY old, est ;
  int j, j0, j1;
  HIT *up ;

  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;

  j0 = j1 = 0 ;
  for (old = 0, j = 0 ; j < arrayMax(hits) ; j ++)
    {
      up = arrp(hits, j, HIT) ;
      est = up->est ;
      if (est != old)
        {
          if (j > j1)
            dropOneHugeIntrons (hits, dna, dnaR, &j0, j1, j, dxMax, gt_ag) ;
          j1 = j ;
        }
      old = up->est ;
    }
  if (j > j1)
    dropOneHugeIntrons (hits, dna, dnaR, &j0, j1, j, dxMax, gt_ag) ;

  arrayMax (hits) = j0 ;
}

/*********************************************************************/

static void dropOneDoubleRead (Array hits, KEYSET clipEnds, int *jp, int j1, int j2, Array dna, Array dnaR)
{
  Array longDna ;
  Array cDNA = 0 ;
  int i, nA, nB, superia, superib, a1, a2, b1, b2 , milieu, amax, bmin, drop, goodAmax, goodBmin, superAmax, superBmin ;
  HIT *up, *upa, *upb ;
  KEY est = 0 ;
  BOOL firstA = TRUE, firstB = TRUE ;
  BOOL ok = FALSE ;

  nA = nB = a1 = a2 = b1 = b2 = goodAmax = goodBmin = superAmax = superBmin = 0 ;

  /* first a hack
     it is possible that on of the reads has problems
     in that case we rather avoid it
  */

  superia = superib = 0 ;
  for (i = j1; i < j2 ; i++)
    {
      up = arrp(hits, i, HIT) ;
      if (up->x1 > up->x2)
        continue ;
      nA++ ;
      if (firstA && est != up->est)
        {
          firstA = FALSE ;
          a1 = up->a1 ; a2 = up->a2 ; 
          if (up->nerr < 30)
            goodAmax = up->a2 ;
          if (up->nerr == 0 && a2 - a1 > 20)
            { superia = i ; superAmax = up->a2 ; }
          else
            superia = i - 1 ;
        }
      else
        {
          if (a1 > up->a1) a1 = up->a1 ;
          if (a2 < up->a2) a2 = up->a2 ;
          if ( goodAmax < up->a2 && up->nerr < 2)
            goodAmax = up->a2 ; 
          if (up->nerr == 0 && a2 - a1 > 20)
            { superia = i ; superAmax = up->a2 ; }
          else if (i > j1 && est == (up-1)->est &&
                   (up-1)->nerr < 5 && (up-1)->a2 - (up-1)->a1 > 20)
            { superia = i - 1 ; superAmax = (up - 1)->a2 ; }
        }
      est = up->est ;
    }
  
  for (i = j2 - 1 ; i >= j1 ; i--)
    {
      up = arrp(hits, i, HIT) ; 
      if (up->x1 < up->x2)
        continue ;
      nB++ ;
      if (firstB && est != up->est)
        {
          firstB = FALSE ;
          b1 = up->a1 ; b2 = up->a2 ; 
           if (up->nerr < 30)
             goodBmin = up->a1 ;
          if (up->nerr == 0 && a2 - a1 > 20)
            { superib = i ; superBmin = up->a1 ; }
          else
            superib = i + 1 ;
        }
      else
        {
          if (b1 > up->a1) b1 = up->a1 ;
          if (b2 < up->a2) b2 = up->a2 ;
          if (goodBmin > up->a1 && up->nerr < 2)
            goodBmin = up->a1 ; 
          if (up->nerr == 0 && a2 - a1 > 20)
            { superib = i ; superBmin = up->a1 ; }
          else if (i < j2 - 1 && est == (up+1)->est &&
                   (up+1)->nerr < 5 && (up+1)->a2 - (up+1)->a1 > 20)
            { superib = i + 1 ; superBmin = (up + 1)->a1 ; }
        }
      est = up->est ;
    }
  
  if (!nA || !nB) 
    {
      /* register happy few  */
      for (i = j1; i < j2 ; i++)
        {
          if (i != *jp)
            arr (hits, *jp, HIT) = arr(hits, i, HIT) ; 
          (*jp)++ ;
        } 
      return ;
    }

  if (! (superAmax && superBmin && (superAmax > superBmin))) /* not ideal yet */
    { /* to extend the a ideal */
      if (superia + 1 < j2 &&
          (up = arrp(hits, superia + 1, HIT)) &&
          class(up->est) &&
          up->x1 < up->x2)
        {
          /* look for error free zone in this exon */ 
          BOOL reverse = up->a1 > up->a2 ? TRUE : FALSE ;
          int dz, ux1, maxDna = arrayMax(dna) ;
          char *cp, *cq ;
          int aa1 = up->a1 ;

          cDNA = getSeqDna (up->est) ;
          if (reverse)
            { 
              longDna = dnaR ; 
              aa1 = maxDna - aa1 + 1 ;  
            }
          else
            longDna = dna ;
          
          dz = 0 ; ux1 = up->x1 ;
          cp = arrp(cDNA, ux1, char) ;
          cq = arrp (longDna, aa1, char) ;
          while (*cp == *cq && ux1 < up->x2)
            { ux1++ ; aa1++ ; cp++ ; cq++ ; dz++ ; }
          if (dz > 8)
            {
              superia++ ; 
              superAmax = reverse ?  maxDna - aa1 + 1 : aa1 ; 
            }
        }
    }

  if (! (superAmax && superBmin && (superAmax > superBmin))) /* not ideal yet */
    { /* to extend the a ideal */
      if (superib - 1 >= j1 &&
          (up = arrp(hits, superib - 1, HIT)) &&
          class(up->est) &&
          up->x1 > up->x2)
        {
          /* look for error free zone in this exon */
          BOOL reverse = up->a1 < up->a2 ? TRUE : FALSE ;
          int bb2, dz, uy2, maxDna = arrayMax(dna) ;
          char *cp, *cq ;

          cDNA = getSeqDna (up->est) ;
          bb2 = up->a2 ;
          if (reverse)
            { 
              longDna = dnaR ; 
              bb2 = maxDna - bb2 + 1 ; 
            }
          else
            longDna = dna ;
          
          dz = 0 ; uy2 = up->x2 ;
          cp = arrp(cDNA, uy2, char) ;
          cq = arrp (longDna, bb2, char) ;
          while (*cp == *cq && uy2 < up->x1)
            { uy2++ ; bb2++ ; cp++ ; cq++ ; dz++ ; }
          if (dz > 8)
            {
              superib-- ; 
              superBmin = reverse ?  maxDna - bb2 + 1 : bb2 ;
            }
        }

    }

  if (superBmin < goodBmin) /* happens beacuse of above extension */
    goodBmin = superBmin ;
  if (superAmax > goodAmax) /* happens beacuse of above extension */
    goodAmax = superAmax ;
 
 if (superAmax && goodBmin && superAmax > superBmin + 30)
    {
      goodAmax = superAmax ; goodBmin = superBmin ; 
    }
  else if (superBmin && superBmin < goodAmax - 30)
        goodAmax = superBmin + 30 ;
  
  amax = a2 ; bmin = b1 ; milieu = (amax + bmin)/2 ; 
  if (goodBmin && a2 > goodBmin + 200) /* see 4b586, the err of b are all at 3' end */
    { amax = goodBmin + 120 ; if (amax < goodAmax) amax = goodAmax ; }
  else if (goodBmin && a2 > goodBmin + 100) /* see 4b586, the err of b are all at 3' end */
    { amax = goodBmin + 60 ; if (amax < goodAmax) amax = goodAmax ; }
  else if (goodBmin && a2 > goodBmin + 30)
    { amax = goodBmin + 30 ; if (amax < goodAmax) amax = goodAmax ; }
  if (goodAmax && b1 > 0 && b1 < goodAmax - 100)
    { bmin = goodAmax - 120 ; if (goodBmin && bmin > goodBmin) bmin = goodBmin ; }
  else if (goodAmax && b1 > 0 && b1 < goodAmax - 100)
    { bmin = goodAmax - 60 ; if (goodBmin && bmin > goodBmin) bmin = goodBmin ; }
  else if (goodAmax && b1 > 0 && b1 < goodAmax - 30)
    { bmin = goodAmax - 30 ; if (goodBmin && bmin > goodBmin) bmin = goodBmin ; }
  if (bmin > 0 && amax > bmin + 30)
    { milieu = (amax + bmin)/2 ; 
      amax = milieu + 15 ; bmin = milieu - 15 ; 
    }
  /* do not stop inside an intron */
  upa = upb = 0 ;
  for (i = j1; i < j2 ; i++)
    {
      up = arrp(hits, i, HIT) ;
      if (up->x1 > up->x2)
        continue ;
      else
        {
          if (amax < up->a2)
            {
              if (amax < up->a1 + 30)
                {
                  amax = up->a1 + 30 ;
                  if (amax >= up->a2)
                    amax = up->a2 - 1 ;
                }
              upa = up ;
              superia = i ;
              break ;
            }
        }
    }
  for (i = j2 - 1; i >= j1 ; i--)
    {
      up = arrp(hits, i, HIT) ;
      if (up->x1 < up->x2)
        continue ;
      else
        {
          if (bmin > up->a1)
            {
              if (bmin > up->a2 - 30)
                {
                  bmin = up->a2 - 30 ;
                  if (bmin <= up->a1)
                    bmin = up->a1 + 1 ; 
                }
              upb = up ;
              superib = i ;
              break ;
            }
        }
    }
  /* mieg 2006_10_30
   * verify that either we cut in the same exon, 
   * or on both side of a single exon
   * but never over 2 introns (yk194b12, OST_30001C2
   */
  if (1)
    {
      ok = FALSE ;
      if (upa && upb)
        {
          if (amax >= upb->a1 && amax <= upb->a2 &&
              bmin >= upa->a1 && bmin <= upa->a2)
            ok = TRUE ; /* same exon */
          else if (superia > j1 && superib < j2 - 1 &&
                   upa->est == (upa-1)->est && 
                   upb->est == (upb+1)->est && 
                   amax >= (upb+1)->a1 && amax <= (upb+1)->a2 &&
                   bmin >= (upa-1)->a1 && (bmin+1) <= upa->a2)
            ok = TRUE ; /* common jumped intron */
          else  /* on a un probleme on va couper ce qui perd le moins */
            {
              int c1, c2, ga = 0, gb = 0 ; /* total length to be cut in a and b */
              c1 = upa->a1 > upb->a1 ? upa->a1 : upb->a1 ;
              c2 = upa->a2 < upb->a2 ? upa->a2 : upb->a2 ;
              if (c1 < c2 && milieu >= c1 && milieu <= c2) /* they intersect but not totally */
                {
                  if (c1 >= upa->a1)
                    {
                      if (c2 - c1 > upa->a2 - c2) amax = c2 ;
                      else bmin = c2 + 1 ;
                    }
                  else if (c2 < upb->a2)
                    {
                      if (c2 - c1 > c1 - upb->a1) bmin = c1 ;
                      else amax = c1 - 1 ;
                    }
                  else /* both overshoot the other exon */
                    {
                      amax = c2 ; bmin = c1 ;
                    }
                }
              else /* no intersect */
                {
                  for (up = upa, i = superia ; i >= j1 ; i--, up--)
                    if (up->a2 < bmin) 
                      break ;
                    else 
                      ga += ((amax < up->a2) ? amax : up->a2) 
                        - ((up->a1 > bmin) ? up->a1 : bmin)
                        + 1 ;
                  for (up = upb, i = superib ; i < j2 ; i++, up++)
                    if (up->a1 > amax) 
                      break ;
                    else 
                      gb += ((amax < up->a2) ? amax : up->a2) 
                        - ((up->a1 > bmin) ? up->a1 : bmin)
                        + 1 ;
                  if (ga < gb &&
                      superia > j1 &&  
                      upa->est == (upa-1)->est)
                    amax = (upa-1)->a2 ;
                  else if (ga > gb &&
                           superib < j2 - 1 &&  
                           upb->est == (upb+1)->est)
                    bmin = (upb+1)->a1 ;
                  else
                    amax = bmin = milieu ;
                }
            }
        }   
    }
  /* register happy few  */
  for (i = j1; i < j2 ; i++)
    { /* drop overshooting */
      up = arrp(hits, i, HIT) ;
      /* do not trim reads that have a polyA 
         if (clipEnds && KEYKEY(up->est) < keySetMax(clipEnds) &&
         (class(keySet(clipEnds, KEYKEY(up->est))) & 1)) ;
      */
      if (TRUE)
        {          
          if (up->x1 > up->x2)
            {
              if (up->reverse ||  
                  (!up->reverse && (b1 > a1 + 30)))
                {
                  if (up->a2 < bmin)
                    continue ;
                  if (up->a1 < bmin)
                    { 
                      drop = bmin - up->a1 ;
                      up->a1 += drop ;
                      up->x1 -= drop ; 
                      up->type |= gLink ;
                    }
                }
            }
          else
            {
              if (up->reverse || 
                  (!up->reverse && (a2 < b2 - 30))) /* do not cut  reads hitting all the way to the end */
                {
                  
                  if (up->a1 > amax)
                    continue ;
                  if (up->a2 > amax)
                    { 
                      drop = up->a2 - amax ;
                      up->a2 -= drop ;
                      up->x2 -= drop ;
                      up->type |= gLink ;
                    }
                }
            }
        }
      if (i != *jp)
        arr (hits, *jp, HIT) = arr(hits, i, HIT) ; 
      (*jp)++ ;
    } 
  i = ok ; /* to please the compiler */
}

/*********************************************************************/

static void dropDoubleReads (Array hits, KEYSET clipEnds, Array dna, Array dnaR)
{
  KEY old, clone ;
  int j, j0, j1;
  HIT *up ;

  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;

  j0 = j1 = 0 ;
  for (old = 0, j = 0 ; j < arrayMax(hits) ; j ++)
    {
      up = arrp(hits, j, HIT) ;
      clone = up->cDNA_clone ;
      if (clone != old)
        {
          if (j > j1)
            dropOneDoubleRead (hits, clipEnds, &j0, j1, j, dna, dnaR) ;
          j1 = j ;
        }
      old = up->cDNA_clone ;
    }
  if (j > j1)
    dropOneDoubleRead (hits, clipEnds, &j0, j1, j, dna, dnaR) ;

  arrayMax (hits) = j0 ;
}

/*********************************************************************/

static void mergeOverlappingHits (Array hits)
{
  int j1, j2, milieu, dx, da ;
  HIT *up1, *up2 ;
  BOOL isUp1, isUp2 ;

  cDNASwapX (hits) ;
  arraySort (hits, cDNAOrderByX1Merge) ;

  for (j1 = 0 ; j1 < arrayMax(hits) ; j1++)
    {
      up1 = arrp(hits, j1, HIT) ;
      if (!up1->est)
        continue ;
      if (up1->a1 == up1->a2 || up1->x1 == up1->x2) /* drop nils */
        { up1->est = 0 ; continue ; }

    }

  /* merge contiguous hits */
  for (j1 = 0 ; j1 < arrayMax(hits) ; j1++)
    {
      up1 = arrp(hits, j1, HIT) ;
      if (!up1->est)
        continue ;
      isUp1 = up1->a1 <= up1->a2 ? FALSE : TRUE ;

      for (j2 = j1 + 1 ; j2 < arrayMax(hits) ; j2++)
        {
          up2 = arrp(hits, j2, HIT) ;
          if (!up2->est) /* already eaten up */
            continue ;
          if (up2->est != up1->est)  /* goto next j1 */
            break ;

          isUp2 = up2->a1 <= up2->a2 ? FALSE : TRUE ;
          if (isUp1 != isUp2) /* antiparallel */
            continue ;

          dx = up2->x1 - up1->x2 - 1 ; /* should be zero */
          if (dx <= -3 || dx >= 3)
            continue ;
          da = isUp1 ? up1->a2 - up2->a1 - 1 : up2->a1 - up1->a2 - 1 ;
          if (da > -3 && da < 3) /* fuse */
            {
              up2->est = 0 ;
              up1->x2 = up2->x2 ;
              up1->a2 = up2->a2 ;
            }
        }
    }

  /* merge overlapping hits */
  for (j1 = 0 ; j1 < arrayMax(hits) ; j1++)
    {
      up1 = arrp(hits, j1, HIT) ;
      if (!up1->est)
        continue ;
      isUp1 = up1->a1 <= up1->a2 ? FALSE : TRUE ;

      for (j2 = j1 + 1 ; j2 < arrayMax(hits) ; j2++)
        {
          up2 = arrp(hits, j2, HIT) ;
          if (!up2->est) /* already eaten up */
            continue ;
          if (up2->est != up1->est)  /* goto next j1 */
            break ;

          isUp2 = up2->a1 <= up2->a2 ? FALSE : TRUE ;
          if (isUp1 != isUp2) /* antiparallel */
            continue ;
          milieu = (up2->a2 + up2->a1)/2 ;
          if (    /* in case of equality in x better to contract the gene  */
              up1->x1 <= up2->x1 + 5 && /* the +5 is to help absorbing good in genes like yk1924 which has a repeat */
              up1->x2 >= up2->x2 - 5 &&  
              (
               (milieu > up1->a1 &&  milieu < up1->a2) || /* so old aa is inside new */
               (milieu < up1->a1 &&  milieu > up1->a2)
               )
              )
            { up2->est = 0 ; continue ; }    /* absorb new in old, next j2 */         
          milieu = (up1->a2 + up1->a1)/2 ;
          if (up1->x1 >= up2->x1 - 5 &&
              up1->x2 <= up2->x2 + 5 && /* old inside new */
              (
               (milieu > up2->a1 &&  milieu < up2->a2) || /* so old aa is inside new */
               (milieu < up2->a1 &&  milieu > up2->a2)
               )
              )
            { up1->est = 0 ; break ; } /* absorb old in new, next j1 */


          if (up1->x2 >= up2->x1 - 1 &&
              (
               ( !isUp1 && 
                 up1->a2 + 20 >= up2->a1 - 1 + up1->x2 - up2->x1 &&
                 up1->a1 < up2->a1 
                 ) ||
               ( isUp1 && 
                 up1->a2 - 20 <= up2->a1 + 1 - up1->x2 + up2->x1 &&
                 up1->a1 > up2->a1 
                 )
               )
              )
            {
              int dx = up2->x2 - up1->x2, da = up2->a2 - up1->a2 ;
              if (isUp1) da = - da ;
              if (dx - da > -8 && dx - da < 8) /* otherwise we are creating a stupid knot */
                {
                  up1->x2 = up2->x2 ;
                  if (!isUp1 && up1->a2 < up2->a2)
                    up1->a2 = up2->a2 ;
                  else if (isUp1 && up1->a2 > up2->a2)
                    up1->a2 = up2->a2 ;
                  up2->est = 0 ; continue ;
                }
            }
        }
    }
  if (arrayMax(hits))
    {  /* register happy few */
      for (j1 = 0, j2 = 0, up1 = arrp(hits, j1, HIT) ; j1 < arrayMax(hits) ; j1++, up1++) 
        { 
          if (up1->est)
            {
              if (j1 != j2)
                arr (hits, j2, HIT) = arr(hits, j1, HIT) ;  
              j2++ ;
            }
        }
      arrayMax (hits) = j2 ;
    }
}

/*********************************************************************/
/* if isUp, all loops should go backwards to reconcile 3' and 5' reads, CF clone yk1727g10 */
/* do all searchs with & operator to allow for rna editing */
static void slideIntron (KEY cosmid, Array dna, KEY est, Array cDna, Array newHits, Array hits)
{
  int a1, a2, b1, b2, x2, y1, x, z1=0, z2=0, z3, z4, dz, da, db ;
  int test ;
  int i, j1, nerr, ninsert, ndelete, ne, nd, ni ;
  HIT *up, *vp ;
  char *cp, *cq, *cr ;
  BOOL debug = FALSE ;

  for (j1 = 0 ; j1 < arrayMax(hits) - 1 ; j1++)
    {
      up = arrp (hits, j1, HIT) ;
      vp = up+1 ;
      a1 = up->a1 ; a2 = up->a2 ; x2 = up->x2 ;
      b1 = vp->a1 ; b2 = vp->a2 ; y1 = vp->x1 ;
      da = db = 0 ;
      for (test = 0 ; test < 7 ; test++)
        {
          /* test = 0: search gt_ag gc_ag, introduce no error 
           *        1: search gt_ag gc_ag, accept one error upstream
           *        2: search gt_ag gc_ag, accept one deletion in exon 1
           *        3: search gt_ag gc_ag, accept one insert in exon 1
           *        4: search gt_ag gc_ag, accept one deletion in exon 2
           *        5: search gt_ag gc_ag, accept one insert in exon 2
           *        6: search others, introduce no error 
           */
           switch (test)
            {
            case 0: break ;
            case 1: 
              if (b1 - a2 < 60) continue ; /* do not force fix a length polymorphism */
              break ;
            case 6: break ;
              /* june 25 2004, after looking at GOLD Danielle no longer likes indel at intron boundaries */
            case 2: 
            case 3: 
            case 4: 
            case 5: continue ;
            }
           
          /* slide back, slide forth, if slideIsPossible, 
             stop on GT/AG, then AT/AC then GC/AG */
          da = db = nerr = ninsert = ndelete = 0 ;
          switch (test)
            {
            case 0: break ; 
            case 1: nerr = 1 ; break ;
              
            case 2: da = 1 ; db = 0 ; ndelete = 1 ; break ; /* missing last letter in cdna before intron */
            case 3: da =  0 ; db = -1 ; ndelete = 1 ; break ; /* missing first letter in cdna after intron */
            case 4: da = -1 ; db =  0 ; ninsert = 1 ; break ; /* duplication in cdna just before the intron */
            case 5: da =  0 ; db =  1 ; ninsert = 1 ; break ; /* duplication in cdna just after the intron */
              
            case 6: break ;
            }

          /* the search area is z1 z3 on exon 1, z2, z4 on exon 2
             we zip cdna and exon2 to move left, cdna and exon 1 to move right
          */
          z1 = a2 ; z2 = b1 ; x = x2 ; 
          cp = arrp(dna, a2 - 1, char) ; /* last base of top exon */
          cq = arrp(dna, b1 - 2, char) ; /* last base of intron */
          cr = arrp (cDna, x - 1 , char) ;
          ne = nerr ; nd = db ? ndelete : 0 ; ni = db ? ninsert : 0 ;
          while (z1 > a1)
            { 
              if (!(*cq & *cr))
                {
                  if (ndelete && nd)
                    {
                      nd-- ;
                      cq-- ; /* no move on cp cr */
                    }
                  else if (ninsert && ni)
                    {
                      ni-- ;
                      cp-- ; cr-- ; z1-- ; z2-- ; /* no move on cq */
                    }
                  else if (nerr && ne)
                    {
                      if (*cp & *cr) /* otherwise i just trade a point error */
                        ne-- ;
                      cp-- ; cq-- ; cr-- ; z1-- ; z2-- ;
                    }
                  else 
                    break ;
                }
              else
                {
                  cp-- ; cq-- ; cr-- ; z1-- ; z2-- ;
                }
            }
          
          z3 = a2 ; z4 = b1 ; x = y1 ; 
          cp = arrp(dna, a2, char) ; /* first base of intron */
          cq = arrp(dna, b1 - 1, char) ; /* first base of down exon */
          cr = arrp (cDna, x - 1 , char) ;
          ne = nerr ; nd = da ? ndelete : 0 ; ni = da ? ninsert : 0 ;
          while (z4 < b2)
            { 
              if (!(*cp & *cr))
                {
                  if (ndelete && nd)
                    {
                      nd-- ;
                      cp++ ; /* no move on cq cr, z */
                    }
                  else if (ninsert && ni)
                    {
                      ni-- ;
                      cq++ ; cr++ ;  z3++ ; z4++ ; /* no move on cp, z */
                    }
                  else if (nerr && ne)
                    {
                      if (*cq & *cr) /* otherwise i just trade a point error */
                        ne-- ;
                      cp++ ; cq++ ; cr++; z3++ ; z4++ ;
                    }
                  else 
                    break ;
                }
              else
                {
                  cp++ ; cq++ ; cr++; z3++ ; z4++ ;
                }
            }

          /* if (z1 == z3)  no sliding possible but victory is possible ! */

          /* recherche GT en aval de z1 */
          
          cp = arrp(dna,z1+da, char) ; cq = arrp(dna,z2 - 3 +db, char) ;
          for (i=0 ; i <= z3 - z1 ; i++, cp++, cq++) /* GT_AG majeur */
            if ((*cp & G_) && (*(cp+1) & T_) && (*cq & A_) && (*(cq+1) & G_)) 
              goto ok ;
          cp = arrp(dna,z1+da, char) ; cq = arrp(dna,z2 - 3 +db, char) ;
          for (i=0 ; i <= z3 - z1 ; i++, cp++, cq++) /* GT_AG majeur */
            if ((*cp & G_) && (*(cp+1) & T_) && (*cq & A_) && (*(cq+1) &  G_))
              goto ok ;
          cp = arrp(dna,z1+da, char) ; cq = arrp(dna,z2 - 3 +db, char);
          for (i=0 ; i <= z3 - z1 ; i++, cp++, cq++) /* GC_AG mineur */
            if ((*cp & G_) && (*(cp+1) & C_) && (*cq & A_) && (*(cq+1) & G_))
              goto ok ;
          
          if (test == 6)
            {
              cp = arrp(dna,z1, char) ; cq = arrp(dna,z2 - 3, char);
              for (i=0 ; i <= z3 - z1 ; i++, cp++, cq++) /* gt_ag on other strand */
                if ((*cp & C_) && (*(cp+1) & T_) && (*cq & A_) && (*(cq+1) & C_))
                  goto ok ;
              cp = arrp(dna,z1, char) ; cq = arrp(dna,z2 - 3, char);
              for (i=0 ; i <= z3 - z1 ; i++, cp++, cq++) /* gc_ag on other strand */
                if ((*cp & C_) && (*(cp+1) & T_) && (*cq & G_) && (*(cq+1) & C_))
                  goto ok ;
              cp = arrp(dna,z1, char) ; cq = arrp(dna,z2 - 3, char) ;
              for (i=0 ; i <= z3 - z1 ; i++, cp++, cq++)  /* U12 */
                if ((*cp & A_) && (*(cp+1) & T_) && (*cq & A_) && (*(cq+1) & C_))
                  goto ok ;
              
              cp = arrp(dna,z1, char) ; cq = arrp(dna,z2 - 3, char);
              for (i=0 ; i <= z3 - z1 ; i++, cp++, cq++) /* rationalize GT donor */
                if ((*cp & G_) && (*(cp+1) & T_))
                  goto ok ;
              cp = arrp(dna,z1, char) ; cq = arrp(dna,z2 - 3, char);
              for (i=0 ; i <= z3 - z1 ; i++, cp++, cq++) /* rationalize GC donor */
                if ((*cp & G_) && (*(cp+1) & C_))
                  goto ok ;
              cp = arrp(dna,z1, char) ; cq = arrp(dna,z2 - 3, char);
              for (i=0 ; i <= z3 - z1 ; i++, cp++, cq++) /* nemat 2 ? */
                if ((*cq & A_) && (*(cq+1) & G_))      /* rationalize acceptor */
                  goto ok ;
            } /* endif test==6 */
        } /* loop on test */
      
      i = 0 ;  /* failed to slide, pack left */
    ok:      /* succes, we broke out of the test loop with a goto */
      dz = z1 + i - up->a2 ;
      if (dz)
        if (debug) printf("Cosmid %s Moved %s %d %d %d %d by %d bases\n",
                          name(cosmid), name(up->est), up->a2, vp->a1, z1 + i, z2 + i, dz) ;
      up->a2 = z1 + i + da ; vp->a1 = z2 + i + db ; 
      if (up->x1 < up->x2) { up->x2 += dz ; vp->x1 += dz ; }
        else  { up->x2 -= dz ; vp->x1 -= dz ; }
    }
  
  for (j1 = 0 ; j1 < arrayMax(hits) ; j1++)
    array(newHits, arrayMax(newHits), HIT) = arr (hits, j1, HIT) ;
}

/*********************************************************************/
 
static void slideOneIntrons (KEY cosmid, KEY est, Array oneHits, Array hits, Array dna, Array dnaR)
{
  int i1, max = arrayMax(dna) ;
  HIT *vp ; 
  Array dnaLong, cDna, mycDna ;
  BOOL isUp, isXup ;

  if (arrayMax(oneHits))
    {
      vp = arrp(oneHits, 0, HIT) ;
      if (  /* to be always in the direction of the gene */
          (vp->x1 > vp->x2 && !vp->reverse) ||
          (vp->x1 < vp->x2 && vp->reverse) 
          )
        {
          dnaLong = dnaR ; isUp = TRUE ;
          complementHits (oneHits, max, 0, arrayMax(oneHits)) ; 
          cDNASwapA (oneHits) ;
          arraySort (oneHits, cDNAOrderByA1) ;
        }
      else
        {
          dnaLong = dna ; isUp = FALSE ;
        }
      i1 = arrayMax(hits) ;
      vp = arrp(oneHits, 0, HIT) ;
      cDna = getSeqDna (vp->est) ; 
      if (vp->x1 > vp->x2)
        {
          isXup = TRUE ;
          mycDna = dnaCopy (cDna) ;
          reverseComplement (mycDna) ;
          complementXHits (oneHits, arrayMax(cDna), 0, arrayMax(oneHits)) ; 
        }
      else
        {
          isXup = FALSE ;
          mycDna = cDna ;
        }
      slideIntron (cosmid, dnaLong, est, mycDna, hits, oneHits) ;
      if (isXup)
        {
          arrayDestroy (mycDna) ;
          complementXHits (hits, arrayMax(cDna), i1, arrayMax(hits)) ; 
        }
      if (isUp)
        complementHits (hits, max, i1, arrayMax(hits)) ; 
    }
}

/*********************************************************************/

static void slideIntrons (KEY cosmid, Array hits, Array dna, Array dnaR)
{
  Array oldHits = arrayCopy (hits), oneHits = 0 ;
  KEY est, old ;
  int j, j1 ;
  HIT *up ;


  cDNASwapA (oldHits) ;
  arraySort (oldHits, cDNAOrderByA1) ;
  arrayMax (hits) = 0 ;

  oneHits = arrayReCreate (oneHits, 12, HIT) ; j1 = 0 ; old = 0 ;
  for (old = 0, j = 0 ; j < arrayMax(oldHits) ; j ++)
    {
      up = arrp(oldHits, j, HIT) ;
      est = up->est ;
      if (est == 1) continue ; /* remove the ghosts */

      if (class(est) && est != old)
        {
          if (old) slideOneIntrons (cosmid, old, oneHits, hits, dna, dnaR) ;
          oneHits = arrayReCreate (oneHits, 12, HIT) ;  j1 = 0 ; old = est ;
        }
      array (oneHits, j1++, HIT) = *up ;
    }
  if (old) slideOneIntrons (cosmid, old, oneHits, hits, dna, dnaR) ;

  arrayDestroy (oldHits) ;
  arrayDestroy (oneHits) ;
}

/*********************************************************************/
/*********************************************************************/
#ifdef JUNK
static void cDNAReverseOneClone (KEY cosmid, KEY est, Array hits, Array introns)
{
  int j1 ;
  HIT *vp ;
  OBJ obj = 0, obj1;
  KEY clone, key ;
  
  for (j1 = 0 ; j1 < arrayMax(hits) ; j1 ++)
    {
      vp = arrp(hits, j1, HIT) ;
      if (vp->est == est) vp->reverse = vp->reverse ? FALSE : TRUE ;
    }
  for (j1 = 0 ; j1 < arrayMax(introns) ; j1 ++)
    {
      vp = arrp(introns, j1, HIT) ;
      if (vp->est == est) vp->reverse = vp->reverse ? FALSE : TRUE ;
    }
  
  obj = bsUpdate (est) ;
  if (!obj) return ;
  if (bsFindTag (obj, _Reverse))
    bsAddTag (obj, _Direct) ;
  else
    bsAddTag (obj, _Reverse) ;
  bsGetKey (obj, _cDNA_clone, &clone) ;
  bsSave (obj) ;  
  
  if (clone)
    obj = bsUpdate (clone) ;
  if (!clone) return ;
  if (bsGetKey (obj, _Read, &key))
    do
      {
        if (key != est) /* already done */
          {
            obj1 = bsUpdate (est) ;
            if (!obj1) break ;
            if (bsFindTag (obj1, _Reverse))
              bsAddTag (obj1, _Direct) ;
            else
              bsAddTag (obj1, _Reverse) ;
            bsSave (obj1) ;  
          }
      } while (bsGetKey (obj, _bsDown, &key)) ;
  bsSave (obj) ;
}

/*********************************************************************/

static void cdnaReverseDubiousClone (KEY cosmid, Array hits, Array dna, Array dnaR)
{
  Array introns = 0 ;
  int j, j1, n1, n2, n3 ;
  KEY est, old = 0 ;
  HIT *up, *vp ;


  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;

  introns = arrayCreate (40, HIT) ;
  for (old = 0, j = 0, j1 = 0 ; j < arrayMax(hits) ; j ++)
    {
      up = arrp(hits, j, HIT) ;
      est = up->est ;
      if (est == old)
        { 
          vp = arrayp (introns, j1++, HIT) ;
          vp->reverse = up->reverse ;
          vp->a1 = (up -1)-> a2 ; vp->a2 = up->a1 ; 
          
        }
      old = up->est ;
    }
  cDNASwapA (introns) ;
  arraySort (introns, cDNAOrderGloballyByA1) ;

  /* search backwards clones */ 
  for (n3 = 0, j = 0 ; j < arrayMax(introns) ; j ++)
    {
      up = arrp(introns, j, HIT) ;
      n1 = 0 ; n2 = 0 ;
      if (up->reverse) n2++ ; else n1++ ;
      for (j1 = j + 1 ; j1 < arrayMax(introns) ; j1 ++)
        {
          vp = arrp(introns, j1, HIT) ;
          if (vp->a1 > up->a1) break ;
          if (vp->reverse) n2++ ; else n1++ ;
        }
      if (n1 * n2)  /* problem */
        {
          if (n2 == 1 && n1 > 1) /* then n2 is dubious */
            for (j1 = j + 1 ; j1 < arrayMax(introns) ; j1 ++)
              {
                vp = arrp(introns, j1, HIT) ;
                if (vp->a1 > up->a1) break ;
                if (vp->reverse)
                  { n3++ ; cDNAReverseOneClone (cosmid, vp->est, hits, introns) ; }
              }
          else if (n1 == 1 && n2 > 1) /* then n1 is dubious */
            { n3++ ;  cDNAReverseOneClone (cosmid, up->est, hits, introns) ; }
        }
    }
  arrayDestroy (introns) ;
}
#endif
/*********************************************************************/
/*********************************************************************/

static void removeIntronWobbling (KEY cosmid, Array hits, Array dna, Array dnaR)
{
  Array introns = 0 ;
  KEYSET plus = 0, minus = 0 ;
  int j, j1, a1, a2, np, nm ;
  HIT *up, *vp ;
  KEY old = 0 ;
  char *cp ;


  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;

  introns = arrayCreate (40, HIT) ;
  for (old = 0, j = 0, j1 = 0 ; j < arrayMax(hits) ; j ++)
    {
      up = arrp(hits, j, HIT) ;
      if (up->est == old)
        { 
          vp = arrayp (introns, j1++, HIT) ;
          vp->reverse = up->reverse ^ (up->x1 > up->x2) ;
          vp->a1 = (up -1)-> a2 ; vp->a2 = up->a1 ; 
        }
      old = up->est ;
    }
  cDNASwapA (introns) ;
  arraySort (introns,  cDNAOrderGloballyByA1 ) ;
  arrayCompress (introns) ;
  /* search backwards clones */ 
  nm = np = 0 ;
  plus = keySetCreate () ; minus = keySetCreate () ;
  for (j = 0 ; j < arrayMax(introns) ; j ++)
    {
      up = arrp(introns, j, HIT) ;
      a1 = up->a1 ; a2 = up->a2 ;
      for (j1 = j + 1 ; j1 < arrayMax(introns) ; j1++)
        {
          vp = arrp(introns, j1, HIT) ;
          if (vp->a1 > up->a2) break ;
          if (up->reverse != vp->reverse) continue ;
          if (vp->a1 == (a1 + 1))
            {
              if (!up->reverse)
                {
                  cp = arrayp (dna, a1, char) ;
                  if (*cp == G_ && *(cp+1) == T_)
                    keySet (minus, nm++) = a1 + 1 ;
                  if (*(cp + 1) == G_ && *(cp+2) == T_)
                    keySet (plus, np++) = a1 ;
                }
              else
                {
                  cp = arrayp (dna, a1, char) ;
                  if (*cp == C_ && *(cp+1) == T_)
                    keySet (minus, nm++) = a1 + 1 ;
                  if (*(cp + 1) == C_ && *(cp+2) == T_)
                    keySet (plus, np++) = a1 ;
                }
            }
          if (vp->a2 == (a2 + 1))
            {
              if (!up->reverse)
                {
                  cp = arrayp (dna, a2 - 3, char) ;
                  if (*cp == A_ && *(cp+1) == G_)
                    keySet (minus, nm++) = a2 + 1 ;
                  if (*(cp + 1) == A_ && *(cp+2) == G_)
                    keySet (plus, np++) = a2 ;
                }
              else
                {
                  cp = arrayp (dna, a2 - 3, char) ;
                  if (*cp == A_ && *(cp+1) == C_)
                    keySet (minus, nm++) = a2 + 1 ;
                  if (*(cp + 1) == A_ && *(cp+2) == C_)
                    keySet (plus, np++) = a2 ;
                }
            }
        }
    }
  if (keySetMax(plus))
    {
      keySetSort (plus) ; keySetCompress (plus) ;
      for (j = 0 ; j < arrayMax(hits) ; j ++)
        {
          up = arrp(hits, j, HIT) ;
          if (keySetFind (plus, up->a1, 0))
            up->a1++ ; /* creates as desided an error in the est */
          if (keySetFind (plus, up->a2, 0))
            up->a2++ ; /* creates as desided an error in the est */
        }
    }
  if (keySetMax(minus))
    {
      keySetSort (minus) ; keySetCompress (minus) ;
      for (j = 0 ; j < arrayMax(hits) ; j ++)
        {
          up = arrp(hits, j, HIT) ;
          if (keySetFind (minus, up->a1, 0))
            up->a1-- ; /* creates as desided an error in the est */
          if (keySetFind (minus, up->a2, 0))
            up->a2-- ; /* creates as desided an error in the est */
        }
    }
  keySetDestroy (plus) ;
  keySetDestroy (minus) ;
  
  arrayDestroy (introns) ;
}

/*********************************************************************/
/*********************************************************************/
#ifdef JUNK

/**********************/

static void removeOneLongNonClassicIntrons  (Array hits, Array dna, Array dnaR, int *jp, int j1, int j2, int dxMax)
{
  int i, ii, a2, x1, x2 ;
  HIT *up ;

  a2 = -1 ;
  x1 = x2 = -1 ;
  
  x1 =  arr(hits, j1, HIT).x1 ; 
  x2 =  arr(hits, j2 - 1, HIT).x2 ; 
  a2 =  arr(hits, j1, HIT).a2 ;

  for (ii = j1 + 1 ; ii < j2 ; ii++)
    {
      up = arrp(hits, ii, HIT) ;
      if (!class(up->est)) continue ;
      if (a2 >= 0 && 
          dxMax && 
          up->a1 > a2 + dxMax &&
          !removeOneLongIntronIsClassic (hits, dna, dnaR, ii))
        {
          if ((up->x1 - x1) * (x2 - x1) < (x2 - up->x1) * (x2 - x1))
            {
              for (i = j1 ; i < ii ; i++)  arr(hits, i, HIT).est = 0 ;
            }
          else
            {
              for (i = ii ; i < j2 ; i++)  arr(hits, i, HIT).est = 0 ;
            }
        }
      a2 = up->a2 ;
    }
  
  /* register happy few  */
  for (ii = j1; ii < j2 ; ii++)
    {
      up = arrp(hits, ii, HIT) ; 
      if (!up->est) continue ;
      if (ii != *jp)
        {
          arr (hits, *jp, HIT) = arr(hits,ii, HIT) ;  
        }
      (*jp)++ ;
    }  
} /* removeOneLongNonClassicIntrons */

/*********************************************************************/

static void removeLongNonClassicIntrons (KEY cosmid, Array hits, Array dna, Array dnaR, int dxMax)
{
  KEY old, est ;
  int j, j0, j1;
  HIT *up ;

  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;

  j0 = j1 = 0 ;
  for (old = 0, j = 0 ; j < arrayMax(hits) ; j ++)
    {
      up = arrp(hits, j, HIT) ;
      est = up->est ;
      if (est != old)
        {
          if (j > j1)
            removeOneLongNonClassicIntrons (hits, dna, dnaR, &j0, j1, j, dxMax) ;
          j1 = j ;
        }
      old = up->est ;
    }
  if (j > j1)
    removeOneLongNonClassicIntrons (hits, dna, dnaR, &j0, j1, j, dxMax) ;

  arrayMax (hits) = j0 ;
} /* removeLongNonClassicIntrons */
#endif

/*********************************************************************/
/*********************************************************************/

static void cDNAFilterUsedQualities (Array hits, Array usedQualities) 
{
  HIT *up, *vp, *up1 ;
  int ii, jj, ii1, iMax, jMax, a1, a2, nerr, ali ;
  int ok ;
  KEY clone ;

  if (!hits || !usedQualities || !arrayMax(hits))
    return ;
  
  cDNASwapA (usedQualities) ;
  arraySort (usedQualities,  cDNAOrderByA1 ) ;

  jj = 0 ; jMax = arrayMax (usedQualities) ;
  ii = 0 ; up = arrp (hits, 0, HIT) ; iMax = arrayMax (hits) ;
  while (ii < iMax)
    {
      up1 = up ; ii1 = ii ; ok = 0 ;
      nerr = 0 ; ali = 0 ;
      if (up->a1 < up->a2)
        { a1 = up->a1 ; a2 = up->a2 ; }
      else
        { a1 = up->a2 ; a2 = up->a1 ; }
      while (ii < iMax && up1->cDNA_clone == up->cDNA_clone)
        {
          ali += up1->x2 - up1->x1 + 1 ; nerr += up1->nerr ;
          if (a2 < up1->a1) a2 = up1->a1 ;
          if (a2 < up1->a2) a2 = up1->a2 ;
          ii1++ ; up1++ ;
        }
      for (; jj < jMax ; jj++)
        {
          vp = arrp (usedQualities, jj, HIT) ;
          if (vp->cDNA_clone > up->cDNA_clone) break ;
          if (vp->cDNA_clone == up->cDNA_clone)
            {
              ok = 1 ;
              if (vp->x1 - vp->nerr > ali - nerr + 5 ||
                  (vp->x1 - vp->nerr > ali - nerr - 3 && 3*(vp->a2 - vp->a1) < 2*(a2 - a1)))
                ok = 3 ; /* discard , previous qual was better */
              else if (vp->x1 - vp->nerr >= ali - nerr)
                ok = 2 ; /* do nothing */
              break ;
            }
        }
      switch (ok)
        {
        case 0: /* new quality */
          vp = arrayp (usedQualities, arrayMax (usedQualities), HIT) ;
        case 1:
          vp->cDNA_clone = up->cDNA_clone  ;
          vp->x1 = ali ; vp->nerr = nerr ; vp->a1 = a1 ; vp->a2 = a2 ;
          break ;
        case 2:
          break ;
        case 3: /* discard */
          up1 = up ; ii1 = ii ; clone =  up->cDNA_clone ;
          while (ii1 < iMax && up1->cDNA_clone == clone)
            { up1->cDNA_clone = 0 ; ii1++ ; up1++ ; }
          break ;
        }

      ii = ii1 ; up = up1 ;
    }
     
  /* keep happy few */
  for (ii = jj = 0, up = vp = arrp (hits, 0, HIT) ; ii < arrayMax (hits) ; up++, ii++)
    {
      if (! up->cDNA_clone)
        continue ;
      if (vp < up) 
        *vp = *up ;
      vp++; jj++ ; 
    }
  arrayMax (hits) = jj ;

} /* cDNAFilterUsedQualities */

/*********************************************************************/

static void cDNAFilterUsedZones (Array hits, Array usedZones) 
{
  HIT *up, *vp, *up0 ;
  int ii, jj, ii0 ;
  BOOL ok ;

  if (!hits || !usedZones || !arrayMax(hits) || !arrayMax(usedZones))
    return ;
  cDNASwapA (hits) ;
  cDNASwapA (usedZones) ;
  arraySort (hits,  cDNAOrderGloballyByA1 ) ;
  arraySort (usedZones,  cDNAOrderGloballyByA1 ) ;

  for (ii = ii0 = 0, up0 = up = arrp (hits, 0, HIT) ; ii < arrayMax (hits) ; up++, ii++)
    {
      ok = 1 ;
      for (jj = 0, vp = arrp (usedZones, 0, HIT) ; ok && jj < arrayMax (usedZones) ; vp++, jj++)
        if (up->cDNA_clone == vp->cDNA_clone &&
            up->a1 < vp->a2 && up->a2 > vp->a1) /* intersect */
          ok = 0 ;
      if (ok)
        {
          if (up0 != up) *up0 = *up ;
          up0++; ii0++ ;
        }
    }
  arrayMax (hits) = ii0 ;
} /* cDNAFilterUsedZones */

/*********************************************************************/

static BOOL cDNAFilterUsedZonesC1 (int c1, KEY clone, Array usedZones, int reverse, int sens) 
{
  HIT *vp ;
  int jj ;
  static KEY oldClone = 0 ;
  static int oldjj = 0 ;

  if (clone == 0)  /* reinitialise */
    { oldClone = 0 ; return FALSE ; }
  if (!usedZones || !arrayMax(usedZones))
    return FALSE;
  
  if (reverse) c1 = reverse - c1 ; 
  if (clone >= oldClone && oldjj >= 0 && oldjj < arrayMax (usedZones))
    jj = oldjj ;
  oldjj = -1 ;
  for (jj = 0, vp = arrp (usedZones, 0, HIT) ; jj < arrayMax (usedZones) ; vp++, jj++)
    if (clone == vp->cDNA_clone)
      {
        if (oldjj == -1) 
          oldjj = jj ; /* first for this clone */
        if (c1 <= vp->a2 && c1 >= vp->a1) /* intersect */
          {
            if (sens > 0)
              {
                if (! reverse) return vp->a2 - c1 + 1 ;
                else return  c1 - vp->a1 + 1 ; 
              }
            else
              {
                if (! reverse) return c1 - vp->a1 + 1 ;
                else return vp->a2 - c1 + 1 ;
              }
          }
      }
    else if (clone < vp->cDNA_clone)
      break ;

  return FALSE ;
} /* cDNAFilterUsedZonesC1 */

/*********************************************************************/

static void cDNAFilterPreviousLink (int direction, KEY cosmid, KEY link, Array hits)
{
  int j, jj, j1, ii, iks;
  HIT *up ;
  KEYSET ks = 0, ks1 = 0 ;
  OBJ obj = 0 ;
  BSunit *uu ;
  Array units  = 0 ;
  KEY tg ;
  int da ;

  ks1 = keySetCreate () ; iks = 0 ;
  for (jj = 0 ; jj < 2 ; jj++)
    {
      if (jj==0)
        {
          obj = cosmid ? bsCreate (cosmid) : 0 ;
          units = arrayReCreate (units, 60, BSunit) ;
        }
      else
        {
          obj = link ? bsCreate (link) : 0 ;
          units = arrayReCreate (units, 60, BSunit) ;
        }
      if (obj)
        {
          if (bsGetArray (obj, _Transcribed_gene, units, 3))
            {
              for (ii = 0 ; ii < arrayMax(units) ; ii += 3)
                {
                  uu = arrp(units, ii, BSunit) ;
                  tg = uu[0].k ;
                  da = uu[2].i - uu[1].i ; 
                  switch (direction)
                    {
                    case 1: /* reverse strand */
                    case -1:
                      if (da < 0) keySet (ks1, iks++) = tg ;
                      break ;
                    default: /* direct strand */ 
                      if (da > 0) keySet (ks1, iks++) = tg ;
                      break ;
                    }
                }
            }
          bsDestroy (obj) ;
        }
      }
  /* the tracking clones do not belong to any mrna we keep them here */
  ks = query (ks1,  "{ >MRNA ; >cDNA_clone ; >Read} $& { >Read } ; ! Super_long ") ;
  for (j = 0, j1 = 0 ; j < arrayMax(hits) ; j ++)
    {
      up = arrp(hits, j, HIT) ;
      if (keySetFind (ks, up->est, 0))
        continue ; /* eliminate that hit */
      if (j != j1) arr (hits,j1,HIT) = *up ;  /* EN PLACE */
      j1++ ;
    }
  arrayMax(hits) = j1 ;
  keySetDestroy (ks) ;  keySetDestroy (ks1) ; arrayDestroy (units) ;
} /* cDNAFilterPreviousLink */

/*********************************************************************/

static void cDNAFilterDirectionalHits (int direction, Array hits, BOOL verbose)
{
  int nn, j, j1 ;
  HIT *up ;

  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;

  for (j = 0, j1 = 0 ; j < arrayMax(hits) ; j ++)
    {
      up = arrp(hits, j, HIT) ;
      nn = up->x2 - up->x1 ;
      if (up->reverse)
        nn = - nn ;
      /* if (nn == 0) both directions are ok */
      switch (direction)
        {
        case 1: /* reverse strand */
        case -1:
          if (nn > 0) 
            {
              if (verbose)
                messerror ("cDNAFilterDirectionalHits: unexpected drop in %s %d %d, direction=%d",
                           name(up->est), up->x1, up->x2, direction) ;
              continue ;
            }
          break ;
        default: /* direct strand */
          if (nn < 0) 
            {
              if (verbose)
                messerror ("cDNAFilterDirectionalHits: unexpected drop in %s %d %d, direction=%d",
                           name(up->est), up->x1, up->x2, direction) ;
              continue ;
            }
          break ;
        }
      if (j1 < j)
        arr (hits,j1,HIT) = *up ;  /* EN PLACE */
      j1++ ;
    }
  arrayMax(hits) = j1 ;
} /* cDNAFilterDirectionalHits */

/*********************************************************************/

static void cDNAFilterDiscardedHits (KEY cosmid, Array hits)
{
  KEY cDNA_clone ;
  int j, j1, a1, a2 ;
  HIT *up ;
  OBJ Cosmid = 0 ;

  if (getIgnoreDiscard ())
     return ;

  Cosmid = cosmid ? bsCreate (cosmid) : 0 ;
  if (!Cosmid) return ;

  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;

  for (j = 0, j1 = 0 ; j < arrayMax(hits) ; j ++)
    {
      up = arrp(hits, j, HIT) ;
      cDNA_clone = up->cDNA_clone ;
      if (bsFindKey(Cosmid, _Discarded_cDNA, cDNA_clone))
        {
          if (bsGetData (Cosmid, _bsRight, _Int, &a1))
            {
              if (bsGetData (Cosmid, _bsRight, _Int, &a2) &&
                  a1 < up->a1 &&
                  a2 > up->a2)
                continue ; /* eliminate that hit */
            }
          else
            continue ; /* if coordinates of discard are unkown, assume discard */
        }
      if (j != j1)
        arr (hits,j1,HIT) = *up ;  /* EN PLACE */
      j1++ ;
    }
  arrayMax(hits) = j1 ;
  bsDestroy (Cosmid) ;
} /* cDNAFilterDiscardedHits */

/*********************************************************************/

static void cDNAFilterBestAlignments (KEY cosmid, Array hits, int deltaQual)
{
  KEY bestCosmid = 0, oldEst = 0 ;
  int j, j1 ;
  static KEY _Best_cosmid_alignment = 0 ;
  HIT *up ;
  BOOL ok = TRUE ;
  KEYSET ks1, ks2, ks3 ;
  OBJ Est = 0 ;

  if (!_Best_cosmid_alignment)
    lexword2key ("Best_cosmid_alignment", &_Best_cosmid_alignment, _VSystem) ;
  if (!arrayMax(hits) || !_Best_cosmid_alignment)
    return ;
  if (keyFindTag (cosmid, str2tag ("is_gene_tile")))
    return ;
  ks1 = queryKey (cosmid, "{IS *} $| { >Parts }") ;

  for (j = 0, j1 = 0, up = arrp(hits, 0, HIT) ; j < arrayMax(hits) ; j++, up++)
    {
      if (up->est != oldEst)
        {
          ok = TRUE ;
          oldEst = up->est ;

          if (keyFindTag (up->est, str2tag("super_long")))
            ok = TRUE ;
          else if ((Est = bsCreate (up->est)))
            {
              if (bsGetKey (Est, _Best_cosmid_alignment, &bestCosmid))
                do
                  {    /* $| { >subsequence } */
                    ks2 = queryKey (bestCosmid, "{IS *} $| { >Parts }") ;
                    ks3 = keySetAND (ks1, ks2) ;
                    ok = keySetMax (ks3) ? TRUE : FALSE ;
                    keySetDestroy (ks2) ;
                    keySetDestroy (ks3) ;
                  } while (! ok && bsGetKey (Est, _bsDown, &bestCosmid)) ;
              bsDestroy (Est) ;
            }
        }
      if (ok)
        {
          if (j1 < j)
            arr (hits, j1, HIT) = *up ;  /* EN PLACE */
          j1++ ;
        }
    }
  arrayMax(hits) = j1 ; 
  keySetDestroy (ks1) ;
} /* cDNAFilterBestAlignments */

/*********************************************************************/

static int cDNAFilterMultipleHits (KEY cosmid, Array hits)
{
  KEY cDNA_clone, old ;
  int j, j1,  nn2 = 0 ;
  HIT *up ;
  int min = 3 ; /* at least min hits */
 
  cDNASwapX (hits) ;
  arraySort (hits, cDNAOrderByX1) ;

  for (old = 0, j = 0, j1 = 0 ; j < arrayMax(hits) ; j ++)
    {
      up = arrp(hits, j, HIT) ;
      cDNA_clone = up->cDNA_clone ;
      if (cDNA_clone != old &&
	  ! keyFindTag (up->est, _IntMap) &&
          (j + min > arrayMax(hits)  || (up+min-1)->cDNA_clone != cDNA_clone))
        continue ; /* not at least 3 hits to the same est */
      if (cDNA_clone!= old) nn2++ ;
      old = cDNA_clone ;
      if (j != j1)
        arr (hits,j1,HIT) = *up ;  /* EN PLACE */
      j1++ ;
    }
  arrayMax(hits) = j1 ;
  return nn2 ;
} /* cDNAFilterMultipleHits */

/*********************************************************************/

static void extendMultipleHits (KEY cosmid, Array hits, Array dna, Array dnaR)
{
  int j, j1, oldx1 = 0, oldx2 = 0, olda1, olda2, milieu, oldolda1 = 0, oldoldx1 = 0 ;
  KEY oldEst = 0 ;
  HIT *vp ;

  cDNASwapX (hits) ;
  arraySort (hits, cDNAOrderByX1Merge) ;

  for (j = 0, oldEst = 0, olda1 = olda2 = 0, j1 = 0 ; j < arrayMax(hits) ; j++)
    {
      vp = arrp(hits, j, HIT) ;

      if (vp->x1 < vp->clipTop) /* mauvaise amorce */
        continue ;
      if (vp->x2 > vp->clipEnd)
        continue ; 
      if (vp->est == oldEst &&
          vp->x1 >= oldx1 && vp->x2 <= oldx2 &&
          (milieu = (vp->a2 + vp->a1)/2) && 
          (
           (milieu > olda1 &&  milieu < olda2) ||
           (milieu < olda1 &&  milieu > olda2)
           ) &&
          (vp->a2 - vp->a1) * (olda2 - olda1) > 0 && /* parallel */
          ( /* the 2 oligos must be approximatelly in phase AB002292, S_8_14, nov 24 2003 */
           (
            vp->a1 < vp->a2 && 
            (vp->x1 - oldoldx1) - (vp->a1 - oldolda1) < 6 &&
            (vp->x1 - oldoldx1) - (vp->a1 - oldolda1) > -6
            ) ||
           (
            vp->a1 > vp->a2 && 
            (vp->x1 - oldoldx1) + (vp->a1 - oldolda1) < 6 &&
            (vp->x1 - oldoldx1) + (vp->a1 - oldolda1) > -6
            )
           ) 
          )
        continue ;  /* drop overlapping oligos */
      oldEst = vp->est ; 
      oldoldx1 = vp->x1 ; 
      oldolda1 = vp->a1 ; 
      extendHits (cosmid, dna, dnaR, vp) ;
      oldx1 = vp->x1 ; oldx2 = vp->x2 ;
      olda1 = vp->a1 ; olda2 = vp->a2 ;
      if (vp != arrp (hits, j1, HIT))
        array(hits, j1, HIT) = *vp ;
      j1++ ;
    }
  arrayMax (hits) = j1 ; /* drop overlapping oligos */
  cDNASwapX (hits) ;
  arraySort (hits, cDNAOrderByX1Merge) ;
  arrayCompress (hits) ;
} /* extendMultipleHits */

/*********************************************************************/

static void dropEst (Array hits, int ii) 
{
  int i ;
  
  for (i = ii ; i < arrayMax(hits) - 1 ; i++)
    array (hits, i, HIT) = array(hits, i+1, HIT) ;
  arrayMax (hits) -= 1 ;
}

static BOOL getBestUpPosition (Array cDNA, Array longDna,
                               int x1, int *x2p, 
                               int a1, int *a2p, BOOL ignoreN, BOOL slide)
{
  int i, nz, nz1, iz, jz, dz, bestnz1, bestdz ;
  char *cp, *cq, *cpmin, *cpmax, *cqmin, *cqmax ;
  int x2 = *x2p, a2 = *a2p, loop = 0, loop2 = 0, nN ;
  int LIMIT = 12 ; /* was 8 in july 2002 */
  static Array err = 0 ;

  loop2 = 0 ;
  cpmin = arrp (cDNA, 0, char)  ;
  cpmax = arrp (cDNA, arrayMax(cDNA) - 1, char)  ;
  cqmin = arrp (longDna, 0, char)  ;
  cqmax = arrp (longDna, arrayMax(longDna) - 1, char)  ;
  if (0 && a1 > *a2p) invokeDebugger() ;
 toulao:
  nz = x2 - x1 ;  
  if (nz > 18) nz = 18 ;
  loop = 0 ;
 lao:
  bestnz1 = 0 ; bestdz = 0 ;
  for (iz = 0 ; bestnz1 < LIMIT && iz < 6 ; iz++)
    for (jz = -1 ;  bestnz1 < LIMIT && jz < 2 ; jz += 2)
      {
        if (iz == 0 && jz == -1) continue ;
        dz = iz * jz ;
        
        if (x2 - nz < 0 || x2 - nz >= arrayMax(cDNA) ||
            a2 - nz - dz< 0 || a2 - nz - dz >= arrayMax(longDna))
          continue ;
        cp = arrp(cDNA, x2 - nz, char) ;
        cq = arrp(longDna, a2 - nz - dz, char) ;
        
        nz1 = 0 ;
        if (ignoreN)
          while ((*cp & *cq) && cp <= cpmax && cq <= cqmax) { cp++ ; cq++ ; nz1++ ;}
        else
          while (*cp == *cq  && cp <= cpmax && cq <= cqmax) { cp++ ; cq++ ; nz1++ ;}
        if (loop2 &&  cp <= cpmax && cq <= cqmax) /* accept one error and 2 N */
          {
            cp++ ; cq++ ; nN = 2 ;
            while (nN > 0 && (*cp & *cq) && cp <= cpmax && cq <= cqmax)
              {
                if (*cp != *cq) nN-- ;
                cp++ ; cq++ ; nz1++ ;
              }
          }
        if (nz1 > bestnz1)
          { bestnz1 = nz1 ; bestdz = dz ; }
      }
  if (bestnz1 < LIMIT && nz < x2 - x1 - 10 && loop++ < 1) 
    { nz += 15 ; goto lao ;}
  /* decale jusqu'a la derniere base juste */
  nz1 = 0 ;
  if (bestnz1)
    {
      x2 = x2 - nz + bestnz1 ; a2 = a2 - nz - bestdz + bestnz1 ; 
      
      cp = arrp(cDNA, x2 - 1, char) ;
      cq = arrp(longDna, a2 - 1, char) ;
      
      i = x2 - 1;
      while (i >= 0  && nz1 < LIMIT && (*cp & *cq) && cp >= cpmin && cq >= cqmin) { i-- ; cp-- ; cq-- ; nz1++; }
    }
  if (nz1 < LIMIT && !loop2++)
    { LIMIT = 8 ; goto toulao ; }
  if (nz1 >= 8 && slide) /* try to accept isolated errors */
    {
      int ierr, ux1 = x2 - 1, ux2 = *x2p, ua1 = a2 - 1, ua2 = *a2p ;
      A_ERR *ep1, *ep2 ;

      err = aceDnaTrackErrors(cDNA, ux1, &ux2, longDna, ua1, &ua2, 0, err, 2, 1, FALSE, 0, FALSE) ;
      if (arrayMax(err) < 2)
        { if (ux2 - x2 >= 8) { x2 = ux2 ; a2 = ua2 ;} }
      else
        { 
          for (ierr = 0, ep1 = arrp (err, 0, A_ERR), ep2 = ep1 + 1 ; ierr < arrayMax(err) - 1 ; ierr++, ep1++, ep2++)
           {
             if (ep2->iShort - ep1->iShort >= 8)
               { x2 = ep2->iShort ; a2 = ep2->iLong ; }
             else
               break ;
           }
        }        
    }
  *x2p = x2 ; *a2p = a2 ;
  return nz1 >= 8 ? TRUE : FALSE ;
}

static BOOL getBestVpPosition (Array cDNA, Array longDna,
                             int *y1p, int y2, int *b1p, BOOL ignoreN, BOOL slide)
{
  int i, nz, nz1, iz, jz, dz, bestnz1, bestdz, nN ;
  char *cp, *cq, *cpmin, *cpmax, *cqmin, *cqmax  ;
  int y1 = *y1p, b1 = *b1p, loop = 0, loop2 = 0 ; 
  int LIMIT = 12 ; /* was 8 in july 2002 */

  loop2 = 0 ; 
  cpmin = arrp (cDNA, 0, char)  ;
  cpmax = arrp (cDNA, arrayMax(cDNA) - 1, char)  ;
  cqmin = arrp (longDna, 0, char)  ;
  cqmax = arrp (longDna, arrayMax(longDna) - 1, char)  ;

 toulao:
  nz = y2 - y1 ;  
  if (nz > 18) nz = 18 ;
  loop = 0 ;
 lao:
  bestnz1 = 0 ; bestdz = 0 ;
  for (iz = 0 ; bestnz1 < LIMIT && iz < 6 ; iz++)
    for (jz = -1 ;  bestnz1 < LIMIT && jz < 2 ; jz += 2)
      {
        if (iz == 0 && jz == -1) continue ;
        dz = iz * jz ;
        
        if (y1 - 2 + nz < 0 || y1 - 2 + nz >= arrayMax(cDNA) ||
            b1 - 2 + nz + dz < 0 || b1 - 2 + nz + dz >= arrayMax(longDna))
          continue ;
        cp = arrp(cDNA, y1 - 2 + nz, char) ;
        cq = arrp(longDna, b1 - 2 + nz + dz, char) ;
        
        nz1 = 0 ; i = y1 - 2 + nz ;
        if (ignoreN)       
          while (i >= 0 && (*cp & *cq) && cp >= cpmin && cq >= cqmin) { i-- ; cp-- ; cq-- ; nz1++; }
        else
          while (i >= 0 && (*cp == *cq) && cp >= cpmin && cq >= cqmin) { i-- ; cp-- ; cq-- ; nz1++; }
        if (loop2 && cp >= cpmin && cq >= cqmin) /* accept one error and 2 N */
          {
            cp-- ; cq-- ; i-- ; nN = 2 ; 
            while (i >= 0  && nN > 0 && (*cp & *cq) && cp >= cpmin && cq >= cqmin)
              {
                if (*cp != *cq) nN-- ;
                i-- ; cp-- ; cq-- ; nz1++ ;
              }
          }
        if (nz1 > bestnz1)
          { bestnz1 = nz1 ; bestdz = dz ; }
      }
  if (bestnz1 < LIMIT && nz < y2 - y1 - 10 && loop++ < 1) 
    { nz += 15 ; goto lao ;}
  /* decale jusqu'a la premiere base juste */  
  nz1 = 0 ;
  if (bestnz1)
    {
      y1 = y1 + nz - bestnz1 ; b1 = b1 + nz + bestdz - bestnz1 ; 
      
      cp = arrp(cDNA, y1 - 1, char) ;
      cq = arrp(longDna, b1 - 1, char) ;

      while (nz1 < LIMIT && (*cp & *cq) && cp <= cpmax && cq <= cqmax) { i++ ; cp++ ; cq++ ; nz1++; }
    }
  if (nz1 < LIMIT && !loop2++)
    { LIMIT = 8 ; goto toulao ; }
  if (nz1 >= 8 && slide) /* try to accept isolated errors, i cannot use track on reverse strand
                          so it only works for one point error 
                         */
    {
      int uy1 = y1 - 1, ub1 = b1 - 1, ui = 0 ;
      cp = arrp(cDNA, uy1, char) - 1 ;
      cq = arrp(longDna, ub1, char)  - 1 ;
      while (cp >= cpmin && cq >= cqmin && uy1 >= *y1p && *cp == *cq)
        { cp-- ; cq-- ; uy1--; ub1--; ui++ ;}
      if (ui >= 8)
        { y1 = uy1; b1 = ub1 ; }   
    }
  *y1p = y1 ; *b1p = b1 ;
  return nz1 >= 8 ? TRUE : FALSE ;
} /* getBestVpPosition */

/**************************************************************************************/

static BOOL locateNewExon (KEY clone, HIT *up, HIT *vp, int cmin, int cmax, int zmin, int zmax, int nn, int nz, 
                           int *nz1p, int *errminp, int *cbestp, int *c2bestp, int *nc1bestp,
                           Array cDNA, Array longDna, Array usedZones, int reversedLongDna)
{
  int 
    c1, c2, c3, dc3, z1, z2, c1Old, nc1, /* new */
    pos, i, imax,
    errmin, maxError,
    cbest, 
    c2best,
    nc1best,
    nz1, 
    maxDna = arrayMax (longDna) ;

  unsigned int oligo,  mask ;
  char *cp ;
#ifdef USEASS
  Associator ass = 0 ;
  void *voidp ;
#else
  static Array estWords = 0 ;
#endif
  static Array err = 0 ;
  static KEYSET c1OldKs = 0 ; /* i make it static because it is quite big, so irealloc it less often */
  BOOL debug = FALSE, oligoOk ;

  /* cherche un 20 mer */
  nz1 = nz > 20 ? 20 : nz ;
  if (up)
    { z1 = zmin ; z2 = zmin + nz1 ; }
  else
    { z2 = zmax ; z1 = zmax - nz1 ; }
  errmin = nz1/4 ;
  cbest = 1 ; c2best = 0 ;
  c1OldKs = arrayReCreate (c1OldKs, zmax + cmax - cmin, int) ;
  
  oligoOk = FALSE ;
  oligo = 0 ; 
  mask = (1<<(2 * nn)) -1 ;
#ifdef USEASS
  assDestroy (ass) ;
  ass = assBigCreate (8*nz) ; /* room for 10 oligos per est */
#else
  estWords = arrayReCreate (estWords, mask, int) ;
#endif
  
  for (i = 0, cp = arrp(cDNA, z1 - 1, char) ; i < nn - 1 ; i++, cp++)
    {
      char fixedbase = *cp ;
      switch(fixedbase)
        {
        case A_: case T_: case G_: case C_:
          fixedbase = * cp ;
          break ;
        case 0:
          goto abort ;
          break ;
        default:
          fixedbase = A_ ;
          break ;
        }
      oligo <<= 2 ;  
      oligo |= B2[(int)(fixedbase)] ; oligo &= mask ; 
    }
  
  imax = nz1 ;
  if (up) { imax = 3 * nz1 ; if (imax > nz - nn) imax = nz - nn ; }
  for (i = 0 ; i < imax ; cp++, i++)   /* was i < nz1 aug2001: i tried i < nz - nn, which in jan 2002 seems bizare  */
    {
      char fixedbase = *cp ;
      switch(fixedbase)
        {
        case A_: case T_: case G_: case C_:
          fixedbase = * cp ;
          break ;
        case 0:
          goto abort ;
          break ;
        default:
          fixedbase = A_ ;
          break ;
        }
      pos = i ;
      oligo <<= 2 ; oligo |= B2[(int)(fixedbase)] ; oligo &= mask ; 
      
#ifdef USEASS 
      if (oligo && checkOligoEntropy (nn, oligo))  /* not the poly A */
        {
          oligoOk = TRUE ;
          assMultipleInsert (ass, assVoid(oligo), assVoid(pos)) ;
        }
#else
      if (oligo && checkOligoEntropy (nn, oligo))  /* not the poly A */
        {
          oligoOk = TRUE ;
          while (array (estWords, oligo, int))
            oligo++ ; /* this leads to inexact oligos but accepts multiple matches at low cost */
          array (estWords, oligo, int) = pos + 1 ;
        }
#endif
      
      if (0 && debug)
        { printf ("z1=%4d  pos=%4d ", z1, pos) ; showOligo(oligo) ; }
    }
#ifdef USEASS 
#else
#endif
  chrono ("getLocateNewExon") ;
  /* now look in cosmid for a hit and work from there */
  if (up)
    {
      int posN = 0 ;
      char *cpMax =  arrp(longDna, arrayMax(longDna) -1, char) ;
      
      chrono ("getLocateNewExon1") ;
      
      oligo = 0 ;  /* preload oligo */
      for (c3 = cmin,  cp = arrp(longDna, c3 - 1, char), i = 0 ;
           i < nn - 1 && cp < cpMax ; c3++, cp++, i++)
        {
          oligo <<= 2 ;
          posN++ ;
          if (!*cp || *cp == N_)
            posN = 0 ;
          else
            oligo |= B2[(int)(*cp)] ; 
          oligo &= mask ; 
        }
      
      for (c1Old = -1, nc1best = 0 ; c3 < cmax - nz1 + nn ; c3++, cp++)
        {
          oligo <<= 2 ;
          posN++ ;
          if (usedZones && (dc3 = cDNAFilterUsedZonesC1 (c3, clone, usedZones, reversedLongDna, 1)))
            { 
              break ; /* no need to intron over a cpy of myself */
              c3 += dc3 ; cp += dc3 ; continue ;
            }
          if (!*cp || *cp == N_)
            posN = 0 ;
          else
            oligo |= B2[(int)(*cp)] ; 
          oligo &= mask ;
          if (posN < nn)
            continue ;
#ifdef USEASS 
          /* 
	     Array bucket = 0 ;
	     int iBucket = 0 ;

            if (checkOligoEntropy (nn, oligo) &&
            assFind (ass, assVoid(oligo), &voidp))
            while (assFindNext(ass, assVoid(oligo), &voidp,&bucket, &iBucket))
            {
            pos = assInt(voidp) ;
          */
#else
          if (checkOligoEntropy (nn, oligo) &&
              array (estWords, oligo, int))
            while ((pos = array (estWords, oligo++, int))) 
              {
                pos-- ;
#endif          
                oligoOk = TRUE ;
                c1 = c3 - nn - pos + 1 ; z1 = zmin ;
                if (c1 < cmin) /* happens if pos is large and c3 = cmin + nn */
                  continue ;
                c1Old = keySet (c1OldKs, c1 - cmin + zmax - z1) ;
                if (c1 > c1Old - 6 && c1 < c1Old + 6) /* already tested */
                  continue ;
                if (usedZones && cDNAFilterUsedZonesC1 (c1, clone, usedZones, reversedLongDna, 1))
                  continue ;                
                keySet (c1OldKs, c1 - cmin + zmax - z1) = c1 ;
                c1 = c3 - nn + 1 ; z1 = zmin + pos ;
                c2 = c1 + nz1 + 10 ; z2 = z1 + nz1 ;
		isComposite |= keyFindTag (up->est, _Composite) ;
		maxError = isComposite ? -2 :  errmin ;
                err = aceDnaTrackErrors (cDNA, z1 , &z2, longDna, c1 , &c2, 0, err, 2, maxError, TRUE, 0, FALSE) ;
                if (0)
                  {
                    printf ("locatenewExons1 : est %s c1:%d c2:%d z1:%d pos:%d nc1:%d nerr:%u\n",
                            name(up->est), c1, c2, z1, pos, nc1, arrayMax(err)) ;
                    showDna(cDNA,z1 - 1) ;
                    showDna(longDna,c1 - 1) ;
                  }
                
                c1 = c3 - nn - pos + 1 ; z1 = zmin ;
                c2 = c1 + nz1 + 10 ; z2 = zmin + nz1 ;
                /* chrono("locatenewExons1") ; chronoReturn () ; */
		maxError = isComposite ? -2 :  errmin ;
                err = aceDnaTrackErrors (cDNA, z1, &z2, longDna, c1, &c2, 0, err, 2, maxError, FALSE, 0, FALSE) ;
                nc1 = c2 - c1 ;
                if (up && vp && arrayMax(err) && c2 - c1 < 12)
                  continue ;
                if (debug && arrayMax(err) < 5)
                  printf ("est %s c1:%d c2:%d z1:%d pos:%d nc1:%d nerr:%u\n",
                          name(up->est), c1, c2, z1, pos, nc1, arrayMax(err)) ;
                if (arrayMax(err) == errmin && 1000*nc1 - (cbest-cmin) > 1000*nc1best - (c1 - cmin))
                  { errmin = arrayMax(err) ; cbest = c1 ; c2best = c2 ; nc1best = nc1 ; }
                if (arrayMax(err) < errmin && !(nc1 <  nc1best - 20))
                  { errmin = arrayMax(err) ; cbest = c1 ; c2best = c2 ; nc1best = nc1 ; }
                /* nov 1, chromo S_20_1 BC004237, je peche un exon de 7 bp au lieu de prolonger mon exon 2 */
                if (arrayMax(err) == errmin && nc1 == nc1best && (c1 - cmin < 5 || cmax - c2 < 5))
                  { errmin = arrayMax(err) ; cbest = c1 ; c2best = c2 ; nc1best = nc1 ; }
                if (0 && errmin == 0)
                  break ;
              }
          if (0 && errmin == 0)
            break ;        
          if (errmin < nz1/4 && c3 > cbest + 500 && cbest > cmin + 500 && cbest > cmin + (cmax - cmin)/3)
            break ;
        }
      chronoReturn () ;
      if (debug && errmin < 5)
        printf ("BEST  est:%s cbest:%d  c2best:%d nc1best:%d errmin:%d\n",name(up->est), cbest, c2best, nc1best, errmin) ;
    }
  else /* search backwards */
    { 
      chrono ("getLocateNewExon2") ;
      
      for (c1Old = -1, nc1best = 0, c3 = cmax - nz1 - 1 ; 
           c3 >= cmin ; c3--)
        {
          int posN = 0 ;

          if (usedZones && (dc3 = cDNAFilterUsedZonesC1 (c3, clone, usedZones, reversedLongDna, -1)))
            { 
              break ; /* no need to intron over a cpy of myself */
              c3 -= dc3 ;  continue ;
            }
          oligo = 0 ; cp = arrp (longDna, c3 - 1, char) ;
          for (i = 0 ; i < nn && i < maxDna - c3 ; i++, cp++)
            {
              oligo <<= 2 ;
              posN++ ;
              if (!*cp || *cp == N_)
                posN = 0 ;
              else
                oligo |= B2[(int)(*cp)] ; 
              oligo &= mask ;
            }
          if (posN < nn)
            continue ;
#ifdef USEASS 
          /* 
	     Array bucket = 0 ;
	     int iBucket = 0 ;
	     if (checkOligoEntropy (nn, oligo) &&
             assFind (ass, assVoid(oligo), &voidp))
               while (assFindNext(ass, assVoid(oligo), &voidp, bucket, &iBucket))             
                 pos = assInt(voidp) ;
          */
#else
          if ( checkOligoEntropy (nn, oligo) &&
               array (estWords, oligo, int))
            while ((pos = array (estWords, oligo++, int)))
              {
                pos-- ;
#endif          
                oligoOk = TRUE ;
                c1 = c3 - pos ; z1 = zmax - nz1 ;
                if (c1 < cmin) /* happens if pos is large and c3 = cmin + nn */
                  continue ;
                if (usedZones && cDNAFilterUsedZonesC1 (c1, clone, usedZones, reversedLongDna, -1))
                  continue ;                
                c1Old = array (c1OldKs, c1 - cmin + zmax - z1, int) ;
                if (c1 > c1Old - 6 && c1 < c1Old + 6) /* already tested */
                  continue ;
                array (c1OldKs, c1 - cmin + zmax - z1, int) = c1 ;
                
                c1 = c3 - pos ; z1 = zmax - nz1 ;
                c2 = c1 + nz1 + 10 ; z2 = zmax ;
		isComposite |= keyFindTag (vp->est, _Composite) ;
		maxError = isComposite ? -2 :  errmin ;
                err = aceDnaTrackErrors (cDNA, z1 , &z2, longDna, c1 , &c2, 0, err, 2, maxError, TRUE, 0, FALSE) ;
                nc1 = c2 - c1 ;
                while (c1 >= cmin && z1 > 0 && arr(cDNA,--z1, char) == arr (longDna, --c1, char))
                  nc1++ ;
                
                c1 = c3 - pos ; z1 = zmax - nz1 ;
                c2 = c1 + nz1 + 10 ; z2 = zmax ;
                
                /* chrono("locatenewExons2") ; chronoReturn () ; */
		maxError = isComposite ? -2 :  errmin ;
                err = aceDnaTrackErrors (cDNA, z1, &z2, longDna, c1, &c2, 0, err, 2, maxError, FALSE, 0, FALSE) ;
		nc1 = c2 - c1 ;
                if (0)
                  {
                    printf ("locatenewExons1 : est %s c1:%d c2:%d z1:%d pos:%d nc1:%d nerr:%u\n",
                            name(vp->est), c1, c2, z1, pos, nc1, arrayMax(err)) ;
                    showDna(cDNA,z1 - 1) ;
                    showDna(longDna,c1 - 1) ;
                  }
                
                if (debug && arrayMax(err) < 8)
                  printf ("est %s c1:%d c2:%d z1:%d pos:%d nc1:%d nerr:%u\n",
                          name(vp->est), c1, c2, z1, pos, nc1, arrayMax(err)) ;
                if (arrayMax(err) == errmin && nc1 > nc1best && 1000*nc1 - (cmax - cbest) > 1000*nc1best - (cmax - c1))
                  { errmin = arrayMax(err) ; cbest = c1 ; c2best = c2 ; nc1best = nc1 ; }
                if (arrayMax(err) < errmin && !(nc1 <  nc1best - 20))
                  { errmin = arrayMax(err) ; cbest = c1 ; c2best = c2 ; nc1best = nc1 ; }
                if (0 && errmin == 0)
                  break ;
              }
          if (0 && errmin == 0)
            break ;
          if (errmin < nz1/4 &&  c3 < cbest  - 500 && cbest < cmax - 500 && cbest < cmax  - (cmax - cmin)/3)
            break ;
        } 
      chronoReturn() ;
      if (debug && errmin < 5)
        printf ("BEST  est:%s cbest:%d  c2best:%d nc1best:%d errmin:%d\n",name(vp->est), cbest, c2best, nc1best, errmin) ;
    }
  /* arrayDestroy (c1OldKs) ; now static */
#ifdef USEASS 
  assDestroy (ass) ;
#endif        
  *nz1p = nz1 ; *errminp = errmin ; *cbestp = cbest ; *c2bestp = c2best ;
  *nc1bestp = nc1best ;
  chronoReturn() ;
  return oligoOk ;
 abort:
  *c2bestp = -1 ;

  return TRUE ; /* stop looping */
} /* locateNewExon */

/**************************************************************************************/

static BOOL getNewExon (KEY cosmid, KEY est, KEY clone, Array hits, int upIndex, int type, 
                        Array dna, Array dnaR, int maxIntronSize, Array usedZones, int zz1, int zz2)
{
  int 
    v1, v2,
    a1, a2, x1, x2, /* in up */
    b1, b2, y1, y2, /* in vp */
    c1, c2, z1, z2, nc1best, /* new */
    nn, maxOverlap, vpIndex = -1, wpIndex = -1,
    errmin, cbest, c2best, cmin, cmax, zmin, zmax, 
    nz, nz1, maxDna = 0 ;
  Array cDNA = 0, longDna = 0 ;
  BOOL reverse, debug = FALSE ;
  HIT *up = 0, *vp = 0 , *wp = 0 ;
  HIT originalUp, originalVp ;
  int  goodUp = 0, goodVp = 0 ;

  memset (&originalUp, 0, sizeof(HIT)) ;
  memset (&originalVp, 0, sizeof(HIT)) ;
  cmin = cmax = zmin = zmax = 0 ;
  a1 = a2 = b1 = b2 = v1 = v2 = y1 = y2 = x1 = x2 = 0 ;
  switch (type)
    {
    case 1: 
      vpIndex = upIndex ;
      up = 0 ;
      vp = arrp(hits, upIndex, HIT) ;
      v1 = vp->clipTop ; v2 = vp->clipEnd ;
      break ;
    case 2:
      up = arrp(hits, upIndex, HIT) ;
      vp = 0 ;
      v1 = up->clipTop ; v2 = up->clipEnd ;
      break ;
    case 3:
      vpIndex = upIndex + 1 ;
      up = arrp(hits, upIndex, HIT) ;
      vp = arrp(hits, upIndex + 1, HIT) ;
      v1 = vp->clipTop ; v2 = vp->clipEnd ;
      break ;
    }
  if (up)
    { a1 = up->a1 ; a2 = up->a2 ; x1 = up->x1 ; x2 = up->x2 ;  }
  if (vp)
    { b1 = vp->a1 ; b2 = vp->a2 ; y1 = vp->x1 ; y2 = vp->x2 ;  }

  cDNA = getSeqDna (est) ;
  if (!cDNA) return FALSE ;

  chrono ("getNewExon1") ;

  reverse = up ? up->a1 > up->a2 : vp->a1 > vp->a2 ;
  maxDna = arrayMax(dna) ;
  if (reverse)
    { 
      if (zz1 > 0 || zz2 > 0) 
	{ int zz = zz1 ; zz1 = maxDna - zz2 + 1 ; zz2 = maxDna - zz + 1 ; }
      longDna = dnaR ; 
      a1 = maxDna - a1 + 1 ;  a2 = maxDna - a2 + 1 ; b1 = maxDna - b1 + 1 ; b2 = maxDna - b2 + 1 ; 
    }
  else
    longDna = dna ;

  if (up && !vp)
    { 
      originalUp = *up ;
      goodUp = 0 ;
      if (1)
        {
          int ux1 = x2 - 12, ux2 = x2, ua1 = a2 - 12 , ua2 = a2 ;
          if ((goodUp =  getBestUpPosition (cDNA, longDna, ux1, &ux2, ua1, &ua2, TRUE, TRUE)))
            {
              x2 = ux2 ; a2 = ua2 ;
              zmin = ux2 + 1 ; cmin = ua2 + 1 ;
            }
        }
      if (!goodUp)
        {
          goodUp =  getBestUpPosition (cDNA, longDna, x1, &x2, a1, &a2, TRUE, TRUE) ;
          zmin = x2 + 1 ; cmin = a2 + 1 ;
        }

      zmax = zmin + 80 ; 
      if (zmax > v2) zmax = v2 ;
      nz = 3000 ;
      if (zmax - zmin > 10) nz =  10000 ;
      if (zmax - zmin > 20) nz =  20000 ;
      if (zmax - zmin > 50) nz =  maxIntronSize ? maxIntronSize  : 100000 ; /* 2019_06_10, was 100000 */  
      if (maxIntronSize && nz > maxIntronSize)
        nz = maxIntronSize ;
      cmax = arrayMax (longDna) > a2 + nz ? a2 + nz : arrayMax (longDna) ;
    }
  else if (!up && vp)
    { 
      originalVp = *vp ;
      goodVp = 0 ;
      if (1)
        {
          int uy1 = y1, uy2 = y1 + 12, ub1 = b1 ;
          if ((goodVp =  getBestVpPosition (cDNA, longDna, &uy1, uy2, &ub1, TRUE, TRUE)))
            {
              y1 = uy1 ; b1 = ub1 ;
              zmax = uy1 ; cmax = ub1 ;
            }
        }
      if (!goodVp)
        {
          goodVp =  getBestVpPosition (cDNA, longDna, &y1, y2, &b1, TRUE, TRUE) ;
          zmax = y1 ; cmax = b1 ;
        }

      zmin = zmin - 80 ;
      if (zmin < v1) zmin = v1 ; 
      nz = 3000 ;
      if (zmax - zmin > 10) nz =  10000 ;
      if (zmax - zmin > 20) nz =  20000 ;
      if (zmax - zmin > 50) nz =  maxIntronSize ? maxIntronSize  : 100000 ; /* 2019_06_10, was 100000 */
      if (maxIntronSize && nz > maxIntronSize)
        nz = maxIntronSize ;
      cmin = b1 > nz ? b1 - nz : 1 ;
    }
  else if (up && vp)
    { 
      originalUp = *up ;
      originalVp = *vp ;
      goodUp = 0 ;
      if (1) /* added nov 24 2003, AK125716 S_16_10 */
        {
          int ux1 = x2 - 12, ux2 = x2, ua1 = a2 - 12 , ua2 = a2 ;
          if ((goodUp =  getBestUpPosition (cDNA, longDna, ux1, &ux2, ua1, &ua2, TRUE, TRUE)))
            {
              x2 = ux2 ; a2 = ua2 ;
              zmin = ux2 + 1 ; cmin = ua2 + 1 ;
            }
        }
      if (!goodUp)
        {
          goodUp =  getBestUpPosition (cDNA, longDna, x1, &x2, a1, &a2, TRUE, TRUE) ;
          zmin = x2 + 1 ; cmin = a2 + 1 ;
        }
      goodVp = 0 ;
      if (1)
        {
          int uy1 = y1, uy2 = y1 + 12, ub1 = b1 ;
          if ((goodVp =  getBestVpPosition (cDNA, longDna, &uy1, uy2, &ub1, TRUE, TRUE)))
            {
              y1 = uy1 ; b1 = ub1 ;
              zmax = uy1 ; cmax = ub1 ;
            }
        }
      if (!goodVp)
        {
          goodVp =  getBestVpPosition (cDNA, longDna, &y1, y2, &b1, TRUE, TRUE) ;
          zmax = y1 ; cmax = b1 ;
        }

      if (goodUp && goodVp && x2 >= y1 - 1)
        {
          int dx ;

          up->a2 = a2 ; up->x2 = x2 ;
          vp->a1 = b1 ; vp->x1 = y1 ;
          dx = (x2 - y1)/2 ;
          up->x2 -= dx ;
          up->a2 -= dx ;
          vp->a1 += up->x2 + 1 - vp->x1 ;
          vp->x1 = up->x2 + 1 ;
          if (reverse)
            { up->a2 = maxDna - up->a2 + 1 ; vp->a1 = maxDna - vp->a1 + 1 ; } 
          chronoReturn () ;
          return FALSE ; /* i fixed the gene but not created the an exon */

        }
    }
  else if (!up && !vp)
    return FALSE ; /* goto abort ; */

  if (zz1 > 0 && cmin < zz1) cmin = zz1 ; 
  if (zz2 > 0 && cmax > zz2) cmax = zz2 ;
  
  if ( cmin > cmax - 12 ||
      (up && !vp && cmin > arrayMax(longDna) - 12) ||
      (!up && vp && cmax < 12) 
      )
     return FALSE ; /* goto abort ; */
  
  nz = zmax - zmin ;
  chronoReturn() ;
  if (nz < 6)  /* size of hashed word */
    {
      if (up)          /* allonge a droite */
        {
          up->x2 = zmax ; a2 = a2 + zmax - x2 ;
          if (vp && y1 > 10 && b1 > 10) { y1 -= 10 ; b1 -= 10 ; }
          if (reverse)
            { a2 = maxDna - a2 + 1 ; b1 = maxDna - b1 + 1 ; }
          up->a2 = a2 ; if (vp) { vp->x1 = y1 ; vp->a1 = b1 ; }
          return /* zmax > x2 ? TRUE : */ FALSE ; /* FALSE: do not iterate */
        }
      else if (vp && b1 > y1) /* allonge a gauche */
        {
          vp->x1 = zmin ; b1 += zmin - y1 ;  
          if (reverse)  { b1 = maxDna - b1 + 1 ; }
          vp->a1 = b1 ;
          return /* zmin < y1 ? TRUE : */ FALSE ;
        }
      return FALSE ;
    }


  /* hash this 20 mer */
  nn = 6 ; /* mieg jan 16 2002, increased the oligo size from 6 to 8 => looses the 9bp exon of 2G876 */

  /* cherche un 20 mer */
  nz1 = nz > 20 ? 20 : nz ;
  errmin = nz1/4 ;
  cbest = 1 ; c2best = nc1best = 0 ;
  
  while (! locateNewExon (clone, up, vp, cmin, cmax, zmin, zmax, nn, nz, &nz1, &errmin, &cbest, &c2best, &nc1best, cDNA, longDna, usedZones, reverse ? maxDna : 0) && 
         nn < nz - 3 && nn < 12) /* && !usedZones) do not look too hard for repeats (usedZones != 0) */
    { nn += 3 ; }

  chrono ("getNewExon4") ;
  if (c2best < cbest + nz1 - 2)
    goto abort ;

  c1 = cbest ; c2 = c1 + nz1 ;  
  if (up)
    { z1 = zmin ; z2 = zmin + nz1 ; }
  else
    { z2 = zmax ; z1 = zmax - nz1 ; }

  nz = arrayMax(hits) ;
  if (errmin > 0 && !up && vp && errmin > vp->nerr/4 && z2 - nc1best >= vp->x1 - 3)
    goto abort ;
  if (errmin > 0 && !vp && up && errmin > up->nerr/4 && z1 + nc1best <= up->x2 + 3)
    goto abort ;

  if (reverse)
    { /* go back to original coords before registering in hits */
      a1 = maxDna - a1 + 1 ;  a2 = maxDna - a2 + 1 ; b1 = maxDna - b1 + 1 ; b2 = maxDna - b2 + 1 ; 
      c1 = maxDna - c1 + 1 ; c2 = maxDna - c2 + 1 ;
    }

  wp = arrayp(hits, nz, HIT) ; /* make room */
  for (; nz > upIndex ; nz--, wp--) 
    *wp = *(wp - 1) ;
  /* thus, entry upIndex is duplicated */
  /* remember that hits may be relocated */

  if (type != 1) wpIndex = upIndex + 1 ; /* insert after */
  else  wpIndex = upIndex ; /* before */

  wp = arrayp (hits, wpIndex , HIT) ;
  wp->a1 = c1 ; wp->a2 = c2 ; wp->x1 = z1 ; wp->x2 = z2 ;
  if (debug) printf("%s %s: new exon[%d] avant extend  %s %d %d  %d %d\n", 
                name(cosmid), name(est), type, reverse ? "reverse" : "direct",
                wp->a1, wp->a2, wp->x1, wp->x2) ;
  
  up = type == 1 ? 0 : wp - 1 ; 
  vp = type == 2 ? 0 : wp + 1 ;

  /* avoid looping */
  if (
      (up && wp && wp->x2 <= originalUp.x2  &&
        (
         (!reverse && wp->a2 <= originalUp.a2) ||
         (reverse && wp->a2 >= originalUp.a2) 
        ) 
       ) ||
      (vp && wp && wp->x1 >= originalVp.x1 &&
       (
        (!reverse && wp->a1 >= originalVp.a1) ||
        (reverse && wp->a1 <= originalVp.a1) 
        ) 
       )
      )
    goto abort ;
    
  
  if (goodUp && wp->x1 <= x2 + 1)
    { /* set exact matches */
      int dx = x2 + 1 -  wp->x1 ;
      wp->x1 += dx ;
      if (reverse) wp->a1 -= dx ;
      else wp->a1 += dx ;

      up->x2 = x2 ; up->a2 = a2 ;  
      goodUp = 2 ; c1 = wp->a1 ;
    }
  else
    {
      if (!errmin &&
          up && wp && wp->x1 <= up->x2)
        {
          int dx = up->x2 - wp->x1 + 1 ;
          up->x2 -= dx ;
          if (reverse) up->a2 += dx ;
          else up->a2 -= dx ; 
          if (reverse) { up->a1 = maxDna - up->a1 + 1 ;  up->a2 = maxDna - up->a2 + 1 ; }
           getBestUpPosition (cDNA, longDna, x1, &up->x2, up->a1, &up->a2, TRUE, TRUE) ; 
          if (reverse) { up->a1 = maxDna - up->a1 + 1 ;  up->a2 = maxDna - up->a2 + 1 ; }
        }
    }

  if (goodVp && wp->x2 >= y1 - 1)
    { /* set exact matches */
      int dx = wp->x2 - y1 + 1 ;
      wp->x2 -= dx ;
      if (reverse) wp->a2 += dx ;
      else wp->a2 -= dx ;

      vp->x1 = y1 ; vp->a1 = b1 ;
      goodVp = 2 ; c2 = wp->a2 ;
    }
  else
    {
      if (!errmin &&
          vp && wp && wp->x2 >= vp->x1)
        {
          int dx = wp->x2 - vp->x1 + 1 ;
          vp->x1 += dx ;
          if (reverse) vp->a1 -= dx ;
          else vp->a1 += dx ;
          if (reverse) { vp->a1 = maxDna - vp->a1 + 1 ;  vp->a2 = maxDna - vp->a2 + 1 ; }
          getBestVpPosition (cDNA, longDna, &vp->x1, y2, &vp->a1, TRUE, TRUE) ;
          if (reverse) { vp->a1 = maxDna - vp->a1 + 1 ;  vp->a2 = maxDna - vp->a2 + 1 ; }
        }
    }

  extendHits( cosmid, dna, dnaR, wp) ;

  if (goodUp == 2) /* restore */
    { wp->a1 = c1 ; wp->x1 = up->x2 + 1 ; }
  if (goodVp == 2) /* restore */
    { wp->a2 = c2 ; wp->x2 = vp->x1 - 1 ; }

  if (
      (up && wp && wp->x2 <= originalUp.x2  &&
        (
         (!reverse && wp->a2 <= originalUp.a2) ||
         (reverse && wp->a2 >= originalUp.a2) 
        ) 
       ) ||
      (vp && wp && wp->x1 >= originalVp.x1 &&
       (
        (!reverse && wp->a1 >= originalVp.a1) ||
        (reverse && wp->a1 <= originalVp.a1) 
        ) 
       )
      )
    goto abort ;

  chronoReturn() ;
  chrono ("getNewExon5") ;
  
  if (
      (goodUp == 2 && goodVp == 2) ||
      ( !up && goodVp == 2) ||
      ( !vp && goodUp == 2)
      )
    goto gagne ;

  if (up && wp && up->x1 > wp->x1 - 6 && 
      up->x2 > wp->x1 - 50 && /* feb 12 2004 */
      up->x2 < wp->x2 - 18)  /* very rare case, see yd906H10 */
    {
      if (wp->nerr * (originalUp.x2 - originalUp.x1) < originalUp.nerr * (wp->x2 - wp->x1))
        {
        /* eatup up */
          nz = arrayMax(hits) - wpIndex ; 
          while (up++, nz-- > 0)
            *(up -1) = *up ;
          arrayMax(hits) -= 1 ;
          goto gagne ;
      }
      else
        goto abort ;
    }

  if (vp && wp && vp->x2 < wp->x2 + 6 && 
      vp->x1 < wp->x1 + 50 && /* feb 12 2004 ak126465, human chr2 */
      vp->x1 > wp->x1 + 18)  /* very rare case, see ?? */
    {
     
      if (wp->nerr * (originalVp.x2 - originalVp.x1) < originalVp.nerr * (wp->x2 - wp->x1))
        {
          /* eatup up */
          nz = arrayMax(hits) - wpIndex - 2; 
          while (vp++, nz-- > 0)
            *(vp -1) = *vp ;
          arrayMax(hits) -= 1 ;
          goto gagne ;
        }
       else
        goto abort ;
    }
     

  /* a priori the new exon is usefull
     so i clip the overlap favourizing 2;1 the new guy
  */
  if (reverse)
    { 
      cmin = maxDna - cmin + 1 ; cmax = maxDna - cmax + 1 ;
    }

  if (up && wp && !reverse && wp->a1 < cmin - 4)
    {
      int dx = cmin - 4 - wp->a1 ;
      wp->a1 += dx ; wp->x1 += dx ; 
      
      dx = up->a2 - cmin - 2 ;
      up->a2 -= dx ; up->x2 -= dx ;
    }
  if (up && wp && reverse && wp->a1 > cmin + 4)
    {
      int dx = wp->a1 - cmin - 4 ;
      wp->a1 -= dx ; wp->x1 += dx ;
      
      dx = cmin - up->a2 - 2 ;
      up->a2 += dx ; up->x2 -= dx ; 
    }
  if (vp && wp && !reverse && wp->a2 > cmax + 4)
    {
      int ddx = goodVp ? 0 : 2 ;
      int dx = wp->a2 - cmax - (ddx + 2) ;
      wp->a2 -= dx ; wp->x2 -= dx ;
      
      dx = cmax - vp->a1 - ddx ;
      vp->a1 += dx ; vp->x1 += dx ;
    }
  if (vp && wp && reverse && wp->a2 < cmax - 4)
    { 
      int ddx = goodVp ? 0 : 2 ;
      int dx = cmax - (ddx + 2) - wp->a2 ;
      wp->a2 += dx ; wp->x2 -= dx ;
      
      dx = vp->a1 - cmax - ddx;
      vp->a1 -= dx ; vp->x1 += dx ; 
    }

    chronoReturn() ;
  chrono ("getNewExon6") ;

  maxOverlap = 30 ;
 limitOverlap:
  if (vp && wp && wp->x2 >= vp->x1 +  maxOverlap)
    {
      int dx = (wp->x2 - vp->x1 + 1)/3 - maxOverlap/3 ;
      wp->x2 -= dx ;
      vp->x1 += 2 * dx ;
      if (reverse) { vp->a1 -= 2 * dx ; wp->a2 += dx ; }
      else { vp->a1 += 2 * dx ; wp->a2 -= dx ; }

    }
  if (up && wp && wp->x1 <= up->x2 -  maxOverlap)
    {
      int dx = (up->x2 - wp->x1 + 1)/3 -  maxOverlap/3 ;
      up->x2 -= 2 * dx ;
      wp->x1 += dx ;
      if (reverse) { up->a2 += 2 * dx ; wp->a1 -= dx ; }
      else { up->a2 -= 2 * dx ; wp->a1 += dx ; }
    }
  /* now i look if overlap is larger twice what remains of new
     in which case i clip the overlap a second time so that
     the topology will not remove it 
  */
  if (wp && maxOverlap > 6 && wp->x2 - wp->x1 < 2 * maxOverlap)
    {
      maxOverlap -= 6 ;
      goto limitOverlap ;
    }

 gagne:
  if (debug)
    printf("%s %s: new exon[%d]  %s %d %d  %d %d\n", 
           name(cosmid), name(est), type, reverse ? "reverse" : "direct",
           wp->a1, wp->a2, wp->x1, wp->x2) ;

  chronoReturn () ;
  return TRUE ; 
abort: 
  if (wpIndex > -1)
    {
       nz = wpIndex ;  /* deallocate */ 
       wp = arrp (hits, wpIndex, HIT) ; 
       for (; nz + 1 < arrayMax(hits) ; nz++, wp++)
         *wp = *(wp + 1) ;
       arrayMax(hits) -- ;
     }

  if (up) { up = arrp (hits, upIndex, HIT) ; *up = originalUp ; }
  if (vp) { vp = arrp (hits, vpIndex, HIT) ; *vp = originalVp ; }

  /*   arrayDestroy (c1OldKs) ; now static */
  chronoReturn () ;
  return FALSE ;
} /* getNewExon */

/**************************************************************************************/

static void getOtherExons (KEY cosmid, Array hits, Array dna, Array dnaR, int maxIntronSize
			   , Array usedZones, int z01, int z02
			   , KEYSET estMap1, KEYSET estMap2 
			   , int cosmidC1, int cosmidC2
			   )
{
  int j, n, oldx1 = 0, oldx2 = 0, dx, iter=0, iter2 = 0, nHits;
  KEY est, oldEst = 0 ;
  HIT *up, *vp ;
  int e1, e2, z1, z2 ;
  BOOL debug = FALSE ;

  cDNASwapX (hits) ;
  arraySort (hits, cDNAOrderByX1) ;

  n = arrayMax(hits) ;
  if (n) 
    { 
      up = arrayp (hits, 2*n + 64, HIT) ;  /* make room */
      arrayMax(hits) = n ;
    }

  z1 = z01 ; z2 = z02 ;
  for (j = 0, nHits = oldEst = 0 ; j < arrayMax(hits) ; j++)
    {
      up = arrp(hits, j, HIT) ;  /* hits may be extended by getNewExon () */
      vp = j < arrayMax(hits) - 1 ? arrp(hits, j + 1, HIT) : 0 ;
      est = up->est ;
      if (vp && vp->est != est) vp = 0 ;

      if (est != oldEst)
        { 
	  z1 = z01 ; z2 = z02 ;
	  iter = iter2 = oldx1 = oldx2 = 0 ; nHits = arrayMax(hits) ; 
	  e1 = estMap1 ? keySet (estMap1, KEYKEY(est)) : 0 ;
	  e2 = estMap2 ? keySet (estMap2, KEYKEY(est)) : 0 ;
	  if (e1 && e2 && cosmidC1 > 0 && cosmidC1 < cosmidC2)
	    { 
	      z1 = (e1 < e2 ? e1 : e2) ; z2 = (e1 >= e2 ? e1 : e2) ; 
	      z1 = z1 - cosmidC1 ; z2 = z2 - cosmidC1 + 1 ;
	    }
	}      
      if (iter >= 6) /* test case: align him-4 genbank sequence */
        {
          if (iter2 >= 100 && iter2 < 200)
            { iter2 = 200 ; messout ("Looping iter2 at clone %s max(hits)=", name(up->cDNA_clone), arrayMax(hits)) ; }
          if (iter > 100 || iter2 >= 100)
            continue ;
          if (arrayMax(hits) > nHits) /* we are gaining something */
            { iter = 0 ; nHits = arrayMax(hits) ; }
          else  /* forget extending this est */
            iter = 100 ;
          continue ;
        }
      if ((!j || up->est != (up-1)->est) &&
          (up->nerr > 4 || up->x1 + 2 * up->nerr >= up->clipTop + 9))
        {
          if (debug) printf("Call new exon[1]: j=%d %s %s x1=%d x2 = %d ll=%d\n",
                 j, name(cosmid), name(up->est), up->x1,up->x2, up->clipEnd) ;
          if (getNewExon (cosmid, est, up->cDNA_clone, hits, j, 1, dna, dnaR, maxIntronSize, usedZones, z1, z2))
            { j-- ; iter++ ; iter2++ ; oldEst = est ; continue ; }
        }
      if (vp &&  oldEst && vp->est == oldEst && 
          vp->x1 >= oldx1 && vp->x2 <= oldx2)
        { dropEst (hits, j) ; j-- ; iter++ ; iter2++ ; oldEst = up->est ; continue ; }  /* drop overlapping oligos */
      if (vp && 
          (
           (vp->nerr > 4 && up->nerr > 4) ||
           up->x2 != vp->x1 - 1)) /* in these cases, try then round up */
        {
          if (debug) printf("Call new exon[3]: j=%d %s %s x1=%d x2=%d ll=%d\n",
                 j, name(cosmid), name(up->est), up->x1,up->x2, up->clipEnd) ;
          
          if (getNewExon (cosmid, est,  up->cDNA_clone, hits, j, 3, dna, dnaR, maxIntronSize, usedZones, z1, z2))
            { j-- ;  iter++ ; iter2++ ; oldEst = est ; continue ; };
        }
      if (vp && up->x2 >= vp->x1 - 18)  /* was -9 , round up */
        {
          dx = (vp->x1 - up->x2)/2 ; 
          if (dx > 0) /* a big overlap is ok */
            { 
              dx += 3 ; /* so we create a  6 base overlap */
              up->x2 += dx ; vp->x1 -= dx ;
              if (up->a1 < up->a2) { up->a2 += dx ; vp->a1 -= dx ;}
              else { up->a2 -= dx ; vp->a1 += dx ; }
            }
        }
      else if (!vp && (up->nerr > 4 || up->x2 <= up->clipEnd - 9))
        { 
          if (debug) printf("Call new exon[2]: j=%d %s %s x1=%d x2=%d ll=%d\n",
                 j, name(cosmid), name(est), up->x1,up->x2, up->clipEnd) ;
          if (getNewExon (cosmid, est, up->cDNA_clone, hits, j, 2, dna, dnaR, maxIntronSize, usedZones, z1, z2))
              { j-- ;  iter++ ; iter2++ ; oldEst = est ; continue ; }
        }
      if (j >=0 )
        {
          up = arrp(hits, j, HIT) ;
          oldx1 = up->x1 ; oldx2 = up->x2 ;
          oldEst = up->est ;
        }
      else
        { oldEst = 0 ; oldx1 = oldx2 = 0 ; }
    }
} /* getOtherExons */

/*********************************************************************/

static void cdnaAlignCountOneIntron (HIT *ah1, HIT *ah2, 
                                     Array dnaD, Array dnaR,
                                     int *gtp, int *gcp, int *ctp)
{ 
  int dx1 = ah1->x2 - ah1->x1 ;
  int dx2 = ah2->x2 - ah2->x1 ;
  int da = ah2->a1 - ah1->a2 ;
  int a1, a2 ;
  char *cp, feet[6] ;
  Array dna = 0 ;
  BOOL isUp = FALSE ;

  feet[2] = '_' ; feet[5] = 0 ;
  if (dx1 < 0) { isUp = TRUE ; dx1 = - dx1 ; }
  if (ah1->reverse) isUp = !isUp ;

  if (dx2 < 0) dx2 = - dx2 ;
  if (dx1 > 8 && dx2 > 8 && da > 4 && dnaD && dnaR)
    {
      /* a1 a2 = first and last base of Intron, in base 1 mode */
      if (!isUp)
        { dna = dnaD ; a1 = ah1->a2 + 1 ; a2 =  ah2->a1 -1 ; }
      else
        {
          dna = dnaR ; a1 = arrayMax(dna) - ah2->a1 + 2 ; a2 = arrayMax(dna) - ah1->a2 ; 
        }
      cp = arrp (dna, a1 - 1, char) ; /* first base of intron in base 0 mode */
      feet[0] = dnaDecodeChar [(int)(*cp++)] ;
      feet[1] = dnaDecodeChar [(int)(*cp)] ;
      cp = arrp (dna, a2 - 2, char) ; /* prelast base  of intron in base 0 mode */
      feet[3] = dnaDecodeChar [(int)(*cp++)] ;
      feet[4] = dnaDecodeChar [(int)(*cp)] ;

      if (!strcmp (feet,"gt_ag")) (*gtp)++ ;
      if (!strcmp (feet,"gc_ag")) (*gcp)++ ;
      if (!strcmp (feet,"ct_ac")) (*ctp)++ ;
    }  
  return ;
} /* cdnaAlignCountOneIntron */

/*********************************************************************/

static void cdnaAlignCountIntrons (Array hits, Array dnaD, Array dnaR, 
                                   KEY est, KEY cDNA_clone, int j1, 
                                   int *gtp, int *gcp, int *ctp)
{
  HIT *ah1, *ah2 ;
  int j2 ;
  
  *gtp = *gcp = *ctp = 0 ;
  for (j2 = j1 ; j2 + 1 < arrayMax(hits) ; j2++)
    {
      ah1 = arrp(hits, j2, HIT) ;
      ah2 = arrp(hits, j2 + 1, HIT) ;
      if (ah2->est != est)
        break ; 
      if (ah1->reverse == ah2->reverse)
        cdnaAlignCountOneIntron (ah1, ah2, dnaD, dnaR, gtp, gcp, ctp) ;
    }
} /* cdnaAlignCountIntrons */

/*********************************************************************/

static KEYSET saveEst2CosmidQuality (KEY cosmid, Array hits, Array linkPos,
                                     Array dnaD, Array dnaR)
{
  KEY cDNA_clone, est = 0 ;
  int x, j1, j2, nmatch, nerr, nReads, quality, lengthEst, maxExact ;
  int minx, maxx, mina, maxa ;
  int gt, gc, ct ;
  OBJ Est = 0 ;
  HIT *hh, *ah ;
  KEYSET ks = 0 ;
  BOOL isUp = FALSE, isMrna = FALSE ;
  KEYSET genes = keySetCreate() ;
  KEY  _Transfered = str2tag ("Transfered") ;

  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;
  ks = keySetCreate () ;

  for (j1 = 0 ; j1 < arrayMax(hits) ; j1++)
    {
      ah = arrp(hits, j1, HIT) ;
      if (!ah->cDNA_clone || ah->est == 1) 
        continue ;        
      est = ah->est ;
      cDNA_clone = ah->cDNA_clone ;
      if (keySetFind (ks, est, 0))
        continue ;
      keySetInsert (ks, est) ;

      nReads = 0 ;
      nmatch = 0 ;
      nerr = 0 ;
      minx = maxx = mina = maxa = maxExact = 0 ;
      lengthEst = ah->clipEnd - ah->clipTop + 1 ;
      isUp = ah->x1 < ah->x2 ? FALSE : TRUE ;
      for (j2 = j1 ; j2 < arrayMax(hits); j2++)
        {
          ah = arrp(hits, j2, HIT) ;
          if (ah->est != est)
            continue ; 
          if (ah->cDNA_clone != cDNA_clone)
            break ; /* hits are ordered by clone */
          nReads |= ah->reverse ? 0x1 : 0x2 ; 
          if (minx == maxx)
            {  minx = maxx = ah->x1 ; }
          if (ah->x1 < minx) minx = ah->x1 ;
          if (ah->x2 < minx) minx = ah->x2 ;
          if (ah->x1 > maxx) maxx = ah->x1 ;
          if (ah->x2 > maxx) maxx = ah->x2 ;

          if (mina == maxa)
            {  mina = maxa = ah->a1 ; }
          if (ah->a1 < mina) mina = ah->a1 ;
          if (ah->a2 < mina) mina = ah->a2 ;
          if (ah->a1 > maxa) maxa = ah->a1 ;
          if (ah->a2 > maxa) maxa = ah->a2 ;
          
          nerr += ah->nerr ;
          nmatch += ah->a2 - ah->a1 + 1 ;
	  if (maxExact < ah->maxExact) maxExact = ah->maxExact ;
        }
      isMrna = keyFindTag (est, _mRNA) ;


      if (keyFindTag (est, _Transfered))
	nerr = 0 ; /* just evaluate on coverage */

      quality = mrnaQualityEvaluate (lengthEst, nmatch, nerr, isMrna, 0, 0) ;
      if (quality >= 11) /* forget disgusting hits */
        continue ; 
      /* keep only genes which are more left than the begin of cosmid2 */
      if (linkPos && arrayMax(linkPos))
        {
          hh = arrp(linkPos, arrayMax(linkPos) - 1, HIT) ;
          x = hh->x1  < hh->x2 ? hh->x1 : hh->x2 ;
          if (mina >= x)
            continue ; /* drop this gene */
        }
      gt = gc = ct = 0 ;
      cdnaAlignCountIntrons (hits, dnaD, dnaR, est, cDNA_clone, j1, &gt, &gc, &ct) ;
      if ((Est = bsUpdate (est)))
        {
          int a1, a2 ;
          if (isUp) { a1 = maxa ; a2 = mina ;}
          else { a1 = mina ; a2 = maxa ; }
          
          if (gt) bsAddData (Est, _gt_ag, _Int, &gt) ;
          if (gc) bsAddData (Est, _gc_ag, _Int, &gc) ;
          if (ct) bsAddData (Est, _ct_ac, _Int, &ct) ;
          bsAddKey (Est, str2tag("From_cosmid"), cosmid) ;
          bsAddData (Est, _bsRight, _Int, &a1) ;
          bsAddData (Est, _bsRight, _Int, &a2) ;
          bsAddData (Est, _bsRight, _Int, &minx) ;
          bsAddData (Est, _bsRight, _Int, &maxx) ;
          bsAddData (Est, _bsRight, _Int, &lengthEst) ;
          bsAddData (Est, _bsRight, _Int, &nmatch) ;
          bsAddData (Est, _bsRight, _Int, &nerr) ; /* note, here we do not use nerrAll */
          bsAddData (Est, _bsRight, _Int, &quality) ;
	  bsAddData (Est, _bsRight, _Int, &maxExact) ;
          bsSave (Est) ;
          keySetInsert (genes, est) ;
        }
    }
  keySetDestroy (ks) ;
  return genes ;
}

/*********************************************************************/
/*********************************************************************/

static Associator estAssCreate (Array words,
                                KEYSET estSet, BOOL reverse, 
                                KEYSET cDNA_clones, KEYSET clipTops, KEYSET clipEnds, KEYSET trackingClones,
                                KEYSET estMaps, KEYSET estMap1, KEYSET estMap2, 
                                int *n1p, int *n2p, int type, int oligoLength, BOOL isDifficult, 
				KEY chrom, int c1, int c2
				)
{
  KEYSET ks = 0 ;
  int ii, n1 = 0, n2 = 0 ;
  KEY *kp ;
  Associator ass = 0 ;
  int step = 30 ;

  chrom = 0 ; /* avoid filtering on IntMap, becuase the ass is created only once */
  /* create associator and fill it with the est */
  
  chrono ("estAssCreate") ;

  switch (type)
    {
    case 3:
      if (!estSet)
        estSet = keySetCreate () ;
      /* fall thru */
    case 0: case 1: case 2:
      step = 30 ;
      if (estSet)
        { 
          char *cp3 = isDifficult ?  "AND (Difficult_to_align || intmap)" : "AND NOT (Difficult_to_align || intmap)" ;
          if (reverse)
            ks = query (estSet,  messprintf("Reverse AND (NOT Discarded)  AND Is_read %s", cp3)) ;
          else
            ks = query (estSet,  messprintf("(NOT Reverse) AND (NOT Discarded) AND Is_read %s", cp3)) ;
        }
      else
        {
          char *cp1 = "" ; /* "IS yk350e12* AND" */
          char *cp2 = "cDNA_clone  AND Is_read AND (NOT Discarded)" ;
          char *cp3 = isDifficult ?  "AND (Difficult_to_align || intmap)" : "AND NOT (Difficult_to_align || intmap)" ;

         if (reverse)
            ks = query (0, messprintf("FIND Sequence  %s %s %s AND Reverse ", cp1, cp2, cp3)) ;
          else
            ks = query (0, messprintf("FIND Sequence  %s %s %s AND  (NOT Reverse)", cp1, cp2, cp3)) ;
        }
      step = 10 ;
      break ;
    case 9:
      step = 6 ;
      ks = query (estSet, "Is_read && NOT Discarded ") ;
      break ;
    case 1002: case 1001: case 1003:
      step = 10 ;
      if (estSet)
        {
          if (reverse)
            ks = query (estSet, "Reverse AND (NOT Discarded) ") ;
          else
            ks = query (estSet, "(NOT Reverse) AND (NOT Discarded) ") ;
        } 
      step = 10 ;
      break ;
    }
  n1 = ii = ks ? keySetMax(ks) : 0 ;
  if (ii)
    {
      ass = assBigCreate (15 * keySetMax(ks)) ; /* room for 10 oligos per est */
      kp = arrp(ks, 0, KEY) -1  ;
      while (kp++, ii--) 
        {
          if (KEYKEY(*kp) >= (1<<22))
            messcrash ("sorry the alignment subroutine is limited to 4M sequences in the database") ;
          n2 += assSequence (words, *kp, reverse, oligoLength, ass, step, 0, 
                             cDNA_clones, clipTops, clipEnds, trackingClones
			     , estMaps, estMap1, estMap2
			     , chrom, c1, c2
			     , type) ;
        }
      if (n1p) *n1p += n1 ; 
      if (n2p) *n2p += n2 ;
    }

  keySetDestroy (ks) ;
  chronoReturn () ;
  return ass ;
}  /* estAssCreate */

/*********************************************************************/

static void getCosmidHits (KEY cosmid, Array hits, 
                           Associator ass, Associator assR, 
                           Array dna, BOOL isUp, KEYSET clipEnds, int nn, 
			   KEYSET estMaps, KEYSET estMap1, KEYSET estMap2, 
			   KEY cosmidMap, int cosmidC1, int cosmidC2,
			   int direction)
{
  KEY est ;
  const void *vp = 0 ; HIT *up ;
  int i, pos, epos, max, a1, a2, x1, x2, eshift, new = 0 ;
  unsigned int oligo, anchor, smask, mask = (1<<(2 * nn)) - 1 , keyMask = (1<<22) - 1 ;
  char *cp ;
  Array bucket = 0 ;
  int iBucket = 0 ;

  if (! ass && ! assR)
    return ;
  
  i = nn - 1 ;
  oligo = 0 ;
  max = arrayMax (dna) ;
  if (max < 4*nn) return ;
  cp = arrp(dna, 0, char) - 1  ;
  while (cp++, i--)
    { oligo <<= 2 ; oligo |= B2[(int)(*cp)] ; } 
  oligo &= mask ; 
  /* cp-- ;  was added fev 19 */
  pos = -1 ; i = arrayMax(dna) - nn + 1 ; /* was pos = -2 fev 19 */
  if (i > 20) while (pos++, cp++, i--)
    { 
      oligo <<= 2 ; oligo |= B2[(int)(*cp)] ; oligo &= mask ;
      bucket = 0 ; iBucket = 0 ;
      if (oligo && ass && assFind(ass, assVoid(oligo), &vp))
        while (assFindNext(ass, assVoid(oligo), &vp, &bucket, &iBucket))
        {
          anchor = assInt(vp) ;
          est = KEYMAKE (_VSequence, (anchor & keyMask)) ;
          smask = keySet (clipEnds, KEYKEY(est)) ;
          smask >>= 10 ; eshift = 0 ;
          while (smask) { smask >>= 1 ; eshift++ ; }
          epos = ((anchor >> 22) & 1023) << eshift  ;

          if (isUp)
            { a1 = max - pos ; a2 = max - pos - nn + 1 ;}
          else
            { a1 = pos + 1 ; a2 = pos + nn ; }
          x1 = epos + 1 ; x2 = epos + nn ;  
          up = arrayp (hits, arrayMax(hits), HIT) ;
          up->est = est ; up->reverse = FALSE ;
          up->a1 = a1 ; up->a2 = a2 ; up->x1 = x1 ; up->x2 = x2 ;
	  if (estMaps && new++ > 1000000)
	    {
	      dropOutOfZonePreMappedHits (hits, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
	      new = 0 ;
	    }
        }
      bucket = 0 ; iBucket = 0 ;
      if (oligo && assR && assFind(assR, assVoid(oligo), &vp))
        while (assFindNext(assR, assVoid(oligo), &vp, &bucket, &iBucket))
        {
          anchor = assInt(vp) ;
          est = KEYMAKE (_VSequence, (anchor & keyMask)) ;
          smask = keySet (clipEnds, KEYKEY(est)) ;
          smask >>= 10 ; eshift = 0 ;
          while (smask) { smask >>= 1 ; eshift++ ; }
          epos = ((anchor >> 22) & 1023) << eshift ;

          if (isUp)
            { a1 = max - pos ; a2 = max - pos - nn + 1 ;}
          else
            { a1 = pos + 1 ; a2 = pos + nn ; }
          /*
          maxEst = keySet(clipEnds, KEYKEY(est)) ;
          x1 = maxEst - epos ; x2 = maxEst - epos - nn + 1 ;  
          */
          x1 = epos ; x2 = epos -nn + 1 ;
          up = arrayp (hits, arrayMax(hits), HIT) ;
          up->est = est ; up->reverse = TRUE ;
          up->a1 = a1 ; up->a2 = a2 ; up->x1 = x1 ; up->x2 = x2 ;
	  if (estMaps && new++ > 1000000)
	    {
	      dropOutOfZonePreMappedHits (hits, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
	      new = 0 ;
	    }
        }
     }
}

/*********************************************************************/

static int checkIntron (KEY cosmid, Array dna, Array dnaR, Array dnaEst, HIT *up, HIT *vp, BOOL secondPass)
{
  Array dnaLong = 0 ;
  char *cp, *cq ;
  int i, a1, a2, x1, x2, n ;

  a1 = a2 = x1 = x2 = 0 ;
  if (up->x2 + 1 != vp->x1)
    {
      int da = up->a1 < up->a2 ? vp->a1 - up->a2 : -vp->a1 + up->a2 ;
      if (vp->x1 > up->x2 - 6  && vp->x1 < up->x2 + 6  && da > 30 + vp->x1 -  up->x2)
        return 1 ;
      return 0 ;
    }
  /* verify if already checked */
   if ( 1 && secondPass &&
       (up->type & gF) && (vp->type & gD)
        )
     return 2 ; /* excellent intron */
  /* check up  exon */
  if (up->a1 < up->a2) 
    { a1 = up->a1 ; a2 = up->a2 ; dnaLong = dna ; } 
  else if (up->a1 > up->a2) /* reject length zero */
    { a1 = arrayMax(dna) - up->a1 + 1 ; a2 = arrayMax(dna) - up->a2 + 1 ; 
      dnaLong = dnaR ; 
    }
  if (!dnaLong)
    return 0 ;
  cp = arrp (dnaLong, a2 - 1, char) ;
  x1 = up->x1 ; x2 = up->x2 ;
  cq = arrp (dnaEst, x2 - 1, char) ;
  if (a2 < a1 + 8 || x2 < x1 + 8) return 1 ;
  n = 0 ;
  for (i = 0 ; i < 8 && a1 <= a2 - i && x1 <= x2 - i ; i++, cp--, cq--)
    if (*cp != *cq) /* n prevents sliding so do not say : !(*cp & *cq)) */
      n++ ;
  if (n) 
    return 1 ;
  /* check down  exon */
  if (vp->a1 < vp->a2) 
    { a1 = vp->a1 ; a2 = vp->a2 ; dnaLong = dna ; }
  else if (vp->a1 > vp->a2) /* reject length zero */
    { a1 = arrayMax(dna) - vp->a1 + 1 ; a2 = arrayMax(dna) - vp->a2 + 1 ; 
      dnaLong = dnaR ; 
    } 
  if (!dnaLong)
    return 0 ;  /* gap */
  cp = arrp (dnaLong, a1 - 1, char) ;
  x1 = vp->x1 ; x2 = vp->x2 ;
  cq = arrp (dnaEst, x1 - 1, char) ;
  if (a2 < a1 + 8 || x2 < x1 + 8) return 1 ;
  for (i = 0 ; i < 8 && a1 + i < a2 && x1 + i < x2 ; i++, cp++, cq++)
    if (*cq == N_ || !(*cp & *cq))
      n++ ;

  if (n)        /* ill defined intron */
    return 1 ;
  return 2 ;   /* excellent intron */
} /* checkIntron */

/*********************************************************************/

static int checkSL (KEY cosmid, Array dna, Array dnaR, Array dnaEst, HIT *up, BOOL isTop)
{
  Array dnaLong = 0 ;
  char *cp, *cq ;
  int i, a1, a2, x1, x2 ;

  a1 = a2 = x1 = x2 = 0 ;

  /* check up  exon */
  if (up->a1 < up->a2) 
    { a1 = up->a1 ; a2 = up->a2 ; dnaLong = dna ; } 
  else if (up->a1 > up->a2) /* reject length zero */
    { a1 = arrayMax(dna) - up->a1 + 1 ; a2 = arrayMax(dna) - up->a2 + 1 ; 
      dnaLong = dnaR ; 
    }
  if (!dnaLong)
    return 0 ;
  if (isTop)
    {
      cp = arrp (dnaLong, a1 - 1, char) ;
      x1 = up->x1 ; x2 = up->x2 ;
      cq = arrp (dnaEst, x1 - 1, char) ;
      if (a2 < a1 + 4 || x2 < x1 + 4) return 0 ;

      for (i = 0 ; i < 5 ; i++, cp++, cq++)
        if (!(*cp & *cq)) /* allow n */
          return 0 ;
    }
  else
    {
      cp = arrp (dnaLong, a2 - 1, char) ;
      x1 = up->x1 ; x2 = up->x2 ;
      cq = arrp (dnaEst, x2 - 1, char) ;
      if (a1 > a2 + 4 || x1 > x2 + 4) return 0 ;
    
      for (i = 0 ; i < 5 ; i++, cp--, cq--)
        if (!(*cp & *cq)) /* allow n */
          return 0 ;
    }

  return 1 ; 
} /* checkSL */

/*********************************************************************/
/* move loose exon boundaries to exact position */
static void checkExon (KEY cosmid, Array dna, Array dnaR, Array cDNA, HIT *gh, BOOL fuzzyAbove, BOOL fuzzyBelow)
{
  Array longDna = 0;
  int a1, a2, b1, b2, x1, x2, y1, y2, maxDna, dx ;
  BOOL reverse = gh->a1 > gh->a2 ;
  int mygF, mygD ;

  if (gh->nerr < 3)
    return ;
  if (!fuzzyAbove && !fuzzyBelow)
    return ;

  maxDna = arrayMax(dna) ; 
  a1 = gh->a1 ; a2 = gh->a2 ; x1 = gh->x1 ; x2 = gh->x2 ;
  if (reverse)
    { 
      longDna = dnaR ; 
      a1 = maxDna - a1 + 1 ;  a2 = maxDna - a2 + 1 ;
      mygF = gD ; mygD = gF ;
    }
  else
    { 
      longDna = dna ; 
      mygF = gF ; mygD = gD ;
    } 
  if (gh->reverse)
    { int mygg = mygF ; mygF = mygD ; mygD = mygg ; }
  if (gh->reverse ^ reverse) /* exclusive or */
    { BOOL myFuzzy = fuzzyAbove ; fuzzyAbove = fuzzyBelow ; fuzzyBelow = myFuzzy ; }

  y1 = x1 ; y2 = x2 ; b1 = a1 ; b2 = a2 ; /* b y = revised positions */

  dx = gh->nerr > 10 ? 15 : 8 ;
  if (! (gh->type & mygF) && fuzzyBelow)
    {
      if (getBestUpPosition (cDNA, longDna, x1, &y2, b1, &b2, TRUE, TRUE))
        { /* never extend, do not cut too much */
          if (y2 < x2)
            {
              if (y2 > x2 - 20 && y2 > x2 - 8 - 2 * gh->nerr)
                { x2 = y2 ; a2 = b2 ; }
              else
                { x2 -= dx ; a2 -= dx ; }
            }
        }
      else
        { x2 -= dx ;a2 -= dx ; }
    }
  if (! (gh->type & mygD) && fuzzyAbove) 
    {
      if (getBestVpPosition (cDNA, longDna, &y1, x2, &b1, TRUE, TRUE))
        { /* never extend, do not cut too much */
          if (y1 > x1)
            {
              if (y1 < x1 + 20 && y1 < x1 + 8 + 2 * gh->nerr)
                { x1 = y1 ; a1 = b1 ; }
              else
                { x1 += dx ; a1 += dx ;}
            }
        }
      else
        { x1 += dx ; a1 += dx ;}
    }

  if (reverse)
    { 
      a1 = maxDna - a1 + 1 ;  a2 = maxDna - a2 + 1 ;
    }
  gh->a1 = a1 ; gh->a2 = a2 ; gh->x1 = x1 ; gh->x2 = x2 ;
}

/*********************************************************************/

static void cleanGeneReverseHits (Array hits)  /* simply exchange a1, a2, do not touch anything else,
                                              the flags are in the direction of the gene */
{
  int i = arrayMax (hits) ;
  HIT *gh = arrp (hits, 0, HIT) -1 ;

  while (gh++, i--)
    { int dummy = gh->a1 ; gh->a1 = 100000000 - gh->a2 ; gh->a2 = 100000000 - dummy ; }
}


static void cleanGeneHits (KEY cosmid, Array geneHits, KEYSET trackingClones, BOOL upGene)
{
  int limit = 10 ; /* allowed wobble */
  int  j0, j1, j2, newj = 0 ;
  HIT *gh1, *gh2, *ghnew ;
  Array oldHits = 0 ;
  
  if (!arrayMax(geneHits))
    return ;

  if (upGene) cleanGeneReverseHits (geneHits) ;
  for (j2 = 0 ; j2 < arrayMax (geneHits) ; j2++)
    { 
      gh2 = arrayp (geneHits, j2, HIT) ;
      gh2->type &= ~gLink ; /* useful only for clean gene hits */
    }

  cDNASwapA (geneHits) ;

  if (0 && upGene)  /* on n'a pas le bug */
    { 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }

  if (1) /* remove micro introns */
    {
      oldHits = geneHits ; /* no copy needed */
      arraySort (oldHits, cDNAOrderByA1) ; /* respect reads */
      for (j2 = 1 ; j2 < arrayMax (oldHits) - 1 ; j2++)
        { 
          gh2 = arrayp (oldHits, j2, HIT) ;
          if (gh2->zone == 9999)
            continue ;
          if (gh2->nerr > 6 && 4 * gh2->nerr > gh2->a2 - gh2->a1) /* was 7 * nerr, but bad for 1k832 */
            continue ;
          if (gh2->type & gSuspect) 
            continue ;
          if (gh2->reverse != upGene || (gGhost & gh2->type))
            continue ;
          if (!(gh2->type & gMicro))
            continue ;

	  if (1)
	    {
	      (gh2-1)->a2 =  (gh2+1)->a2 ;
	      (gh2-1)->x2 =  (gh2+1)->x2 ;
	      (gh2-1)->type &=  ~gF ;
	      (gh2-1)->type |= ((gh2+1)->type & gF) ;
	      gh2->type |= gSuspect ;
	      (gh2+1)->type |= gSuspect ; /* so gh2+1 will be jumped */
	      (gh2)->type |= gSuspect ;
	      
	      continue ;
	    }
	  /* else more accurate undebuged possibility:  compare and upgrade previous */
          for (j1 = 0 ; j1 <  j2 ; j1++)
            { 
              gh1 = arrayp (geneHits, j1, HIT) ;
              
              if (gh1->reverse != upGene)  /* different strands */
                continue ;
              if ((gh1->type & gSuspect) || !(gX & gh1->type))
                continue ;
              if (gh1->a1 > gh2->a1 - 10 ||  gh1->a2 < gh2->a2 - 10)
                continue ;


              if (gh1->a1 == (gh2-1)->a1 && gh1->a2 == (gh2+1)->a2)
                { 
                  gh1->type |= ((gh2+1)->type & gFin) ;
                  gh1->type |= ((gh2-1)->type & gDebut) ;
                  if ((gFF & gh1->type) && ! (gFF & (gh2+1)->type))
                    gh1->type &= ~gFF ;
                  if ((gDF & gh1->type) && ! (gDF & (gh2-1)->type))
                    gh1->type &= ~gDF ;
                  gh2->type |= gSuspect ;
                  (gh2+1)->type |= gSuspect ; /* so gh2+1 will be jumped */
                  (gh2-1)->type |= gSuspect ;
                }
              else if (gh1->a1  < (gh2-1)->a1 && ((gh2 - 1)->type & gDF))
                {
                  gh2->type |= gSuspect ;
                  (gh2-1)->type |= gSuspect ; 
                  (gh2+1)->type |= gDF ;
                }
              else if (gh1->a2  > (gh2+1)->a2 && ((gh2 + 1)->type & gFF))
                {
                  gh2->type |= gSuspect ;
                  (gh2+1)->type |= gSuspect ;
                  (gh2-1)->type |= gFF ;
                }
              else
                {
                  gh2->type |= gSuspect ; 
                  (gh2+1)->type |= gDF ;
                  (gh2-1)->type |= gFF ;
                }
            }
        }
    }

  if (0 && upGene)  /* on A le bug */
    { 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }


  for (j0 = 0 ; j0 < 2 ; j0++)   /* first one round to match 3p and 5p of one clone 
                                  * second round to remove micro introns
                                  */
    {
      oldHits = arrayCopy (geneHits) ; newj = 0 ;
      geneHits = arrayReCreate (geneHits, arrayMax(oldHits), HIT) ;
      arraySort (oldHits, cDNAOrderGenesGloballyByA1) ; /* Globally since dec 26 99, to cure YK2233 */

      for (j2 = 0 ; j2 < arrayMax (oldHits) ; j2++)
        { 
          gh2 = arrayp (oldHits, j2, HIT) ;
          if (gh2->zone == 9999)
            continue ;

	  if (0 && j2>=0 && upGene)  /* on A le bug already at j0=j2=0 */
	    { 
	      arrayDestroy (geneHits) ;
	      invokeDebugger () ;
	    }

            /* we work with coordinates in the orientation of the gene and flags with their natural meaning */

          if (gh2->nerr > 6 && 4 * gh2->nerr > gh2->a2 - gh2->a1) /* was 7 * nerr, but bad for 1k832 */
            continue ;
          if (gh2->type & gSuspect) 
            continue ;
          if (gh2->reverse != upGene || (gGhost & gh2->type))
            goto done ;
          if (keySetFind (trackingClones, gh2->cDNA_clone, 0))  /* "XM_"  */
            continue ;  /* do not register this gh2 */
          if (gh2->nerr == 0 &&
              (
               gh2->x1 == gh2->clipTop ||
               gh2->x2 == gh2->clipTop ||
               gh2->x1 == gh2->clipEnd ||
               gh2->x2 == gh2->clipEnd
               ) &&
              keyFindTag (gh2->cDNA_clone, _RT_PCR)
              )
            gh2->nerr = 1 ; /* this will make the ends of a pcr product effectivelly fuzzy */

          /* compare and upgrade previous */
          for (j1 = 0 ; j1 < newj ; j1++)
            { 
              gh1 = arrayp (geneHits, j1, HIT) ;
              
              if (gh1->type & gSuspect) 
                continue ;
              if (!j0 && gh1->cDNA_clone != gh2->cDNA_clone)
                continue ;
              
               if (gh1->reverse != upGene)  /* different strands */
                continue ;
              
              if (gh2->type == gh1->type &&
                  gh1->a1 == gh2->a1 && gh1->a2 == gh2->a2)
                goto nextj2 ;    /* do not double */
              
              if ((gI & gh1->type) &&  (gI & gh2->type))
                {
                  if (!(gJ & gh1->type) &&  (gJ & gh2->type)) /* gh1 may win */
                    {
                      if (gh2->a1 > gh1->a1 - 15 && gh2->a1 < gh1->a1 + 15 &&
                          gh2->a2 > gh1->a2 - 15 && gh2->a2 < gh1->a2 + 15)
                        goto nextj2 ; /* forget gh2  */
                    }
                  if ((gJ & gh1->type) &&  (!(gJ & gh2->type)   || (gJ & gh2->type) )) /* gh2 may win */
                      {                                                   /* even if both are weird */
                      if (gh2->a1 > gh1->a1 - 15 && gh2->a1 < gh1->a1 + 15 &&
                          gh2->a2 > gh1->a2 - 15 && gh2->a2 < gh1->a2 + 15)
                        {
                          if ((gJ & gh1->type) &&  (!(gJ & gh2->type)))
                            gh1->type &= ~gJ ;
                          *gh1 = *gh2 ;
                          goto nextj2 ; /* forget gh1 !  */
                        }
                    }

                }
              if (gI & gh1->type)  /* introns */
                continue ;
              if (gI & gh2->type)  /* introns */
                continue ;
              
              if (gh2->a1 > gh1->a2 + 5 || gh2->a2 + 5 < gh1->a1)
                continue ; /* no quasi intersect */
              
              if (gGhost & gh1->type) /* ghosts */
                continue ;

              /********* treat debut of exons ********/
              if ((gD & gh1->type) && (gD & gh2->type)) ;  /* debut debut */
              else if (gD & gh1->type) /* debut/bof */
                {
                  if (gh2->a1 > gh1->a1 -  (gh2->nerr ? limit + gh2->nerr : 5) && gh2->a1 < gh1->a2 &&
                      ! (gh2->a1 >  gh1->a1 + 9 && (gDF & gh2->type)))
                    gh2->a1 = gh1->a1 ;
                }
              else if (gD & gh2->type) /* bof/debut */
                {
                  if (gh1->a1 > gh2->a1 -  (gh1->nerr ? limit + gh1->nerr : 5) && gh1->a1 < gh2->a2 &&
                      ! (gh1->a1 >  gh2->a1 + 9 && (gDF & gh1->type)))
                    gh1->a1 = gh2->a1 ;
                }
              else if ((gDF & gh1->type) && !(gDF & gh2->type))  /* debutFuzzy/bof */
                {
                  if (gh2->a1 > gh1->a1 -  (gh2->nerr ? limit + gh2->nerr : 5) && gh2->a1 < gh1->a2 &&
                      ! (gh2->a1 >  gh1->a1 + 9 && (gDF & gh2->type)))
                    gh2->a1 = gh1->a1 ;
                }
              else if (!(gDF & gh1->type) && (gDF & gh2->type)) /* bof/debutFuzzy */
                {
                  if (gh1->a1 > gh2->a1 -  (gh1->nerr ? limit + gh1->nerr : 5) && gh1->a1 < gh2->a2 &&
                      ! (gh1->a1 >  gh2->a1 + 9 && (gDF & gh1->type)))
                    gh1->a1 = gh2->a1 ;
                }
              else                     /* bof/bof */
                {
                  if (gh1->a1 > gh2->a1 &&
                      ! (gh1->a1 >  gh2->a1 + 9 && (gDF & gh1->type)))
                    gh1->a1 = gh2->a1 ;
                  else if ( gh2->a1 > gh1->a1 &&
                            ! (gh2->a1 >  gh1->a1 + 9 && (gDF & gh2->type)))
                    gh2->a1 = gh1->a1 ;
                }
              if (gh1->a1 == gh2->a1)
                {
                  if ((gDF & gh1->type) && (gDebut & gh2->type) && ! (gDF & gh2->type))
                    gh1->type &= ~gDF ;
                  if ((gDF & gh2->type) && (gDebut & gh1->type) && ! (gDF & gh1->type))
                    gh2->type &= ~gDF ;
                  gh1->type |= (gh2->type & gDebut) ;
                  gh2->type |= (gh1->type & gDebut) ;
                }

              /********* treat fin of exons ********/
              if ((gF & gh1->type) && (gF & gh2->type))  /* fin fin */
                {
                  /* do not rationalise exact intron boundaries , but ok with polyA sites */
                  if ((gA & gh1->type) && (gA & gh2->type))
                    {
                      if (gh1->a2 < gh2->a2 && gh1->a2 > gh2->a2 - 99)
                        gh1->a2 = gh2->a2 ;
                      else if (gh2->a2 < gh1->a2 && gh2->a2 > gh1->a2 - 99)
                        gh2->a2 = gh1->a2 ;
                    }
                }
              else if (gF & gh1->type) /* fin/bof */ /*(gh2->nerr ? da + gh2->nerr : 5) is also in s2mMakeOneShadow */
                {
                  if ((gFF & gh2->type) && gh1->a2 > gh2->a2 + 9) ;
                  else if (gA & gh2->type)  /* do not extend aribitrarily a polyA initiated 3' read */
                    { /* and gA ends are rationalised only to other exact gA ends */
                      if ((gA & gh1->type) &&
                          (gh2->a2 < gh1->a2 + (gh2->nerr ? limit + gh2->nerr : 5) && gh2->a2 > gh1->a1 && 
                           gh2->a2 > gh1->a2 - 99))
                        gh2->a2 = gh1->a2 ;
                    }
                  else
                    {
                      if (gh2->a2 < gh1->a2 + (gh2->nerr ? limit + gh2->nerr : 5) && gh2->a2 > gh1->a1)
                        gh2->a2 = gh1->a2 ;
                    }
                }
              else if (gF & gh2->type) /* bof/fin */
                {
                  if ((gFF & gh1->type) && gh2->a2 > gh1->a2 + 9) ;
                  else if (gA & gh1->type)
                    { 
                      if ((gA & gh2->type) &&
                          (gh1->a2 < gh2->a2 +  (gh1->nerr ? limit + gh1->nerr : 5) && gh1->a2 > gh2->a1 && 
                           gh1->a2 > gh2->a2 - 99))
                        gh1->a2 = gh2->a2 ;
                    }
                  else
                    {
                      if (gh1->a2 < gh2->a2 +  (gh1->nerr ? limit + gh1->nerr : 5) && gh1->a2 > gh2->a1)
                        gh1->a2 = gh2->a2 ;
                    }
                }
              else if ((gFF & gh1->type) && ! (gFF & gh2->type)) /* finFuzzy/bof */ /*(gh2->nerr ? da + gh2->nerr : 5) is also in s2mMakeOneShadow */
                {
                  if (gA & gh2->type)  /* do not extend aribitrarily a polyA initiated 3' read */
                    { /* and gA ends are rationalised only to other exact gA ends */
                      if ((gA & gh1->type) &&
                          (gh2->a2 < gh1->a2 + (gh2->nerr ? limit + gh2->nerr : 5) && gh2->a2 > gh1->a1 && 
                           gh2->a2 > gh1->a2 - 99))
                        gh2->a2 = gh1->a2 ;
                    }
                  else
                    {
                      if (gh2->a2 < gh1->a2 + (gh2->nerr ? limit + gh2->nerr : 5) && gh2->a2 > gh1->a1)
                        gh2->a2 = gh1->a2 ;
                    }
                }
              else if (!(gFF & gh1->type) && (gFF & gh2->type)) /* bof/finFuzzy */
                {
                  if (gA & gh1->type)
                    { 
                      if ((gA & gh2->type) &&
                          (gh1->a2 < gh2->a2 +  (gh1->nerr ? limit + gh1->nerr : 5) && gh1->a2 > gh2->a1 && 
                           gh1->a2 > gh2->a2 - 99))
                        gh1->a2 = gh2->a2 ;
                    }
                  else
                    {
                      if (gh1->a2 < gh2->a2 +  (gh1->nerr ? limit + gh1->nerr : 5) && gh1->a2 > gh2->a1)
                        gh1->a2 = gh2->a2 ;
                    }
                }
              else                     /* bof/bof */
                {
                  if (gh1->a2 < gh2->a2 && gh1->a2 > gh2->a1 &&
                      ! (gh2->a2 >  gh1->a2 + 9 && (gFF & gh1->type)))  
                    {
                      if (gA & gh1->type)
                        { 
                          if (gh1->a2 > gh2->a2 - 99)
                            gh1->a2 = gh2->a2 ;
                        }
                      else
                        {
                          gh1->a2 = gh2->a2 ;
                        }
                    }
                  else if (gh2->a2 < gh1->a2 && gh2->a2 > gh1->a1 &&
                           ! (gh1->a2 >  gh2->a2 + 9 && (gFF & gh2->type)))
                    {
                      if (gA & gh2->type)
                        { 
                          if (gh2->a2 > gh1->a2 - 99)
                            gh2->a2 = gh1->a2 ;
                        }
                      else
                        {
                          gh2->a2 = gh1->a2 ;
                        }
                    }
                }
              /* Intron start remains different from exact polyA */
              if (gh1->a2 == gh2->a2 &&
                  ! ((gF & gh1->type) && (gF & gh2->type) &&
                     ((gA & gh1->type) != (gA & gh2->type)))
                  )
                {
                  if ((gFF & gh1->type) && (gF & gh2->type) && ! (gFF & gh2->type))
                    gh1->type &= ~gFF ;
                  if ((gFF & gh2->type) && (gF & gh1->type) && ! (gFF & gh1->type))
                    gh2->type &= ~gFF ;
                  gh1->type |= (gh2->type & gFin) ;
                  gh2->type |= (gh1->type & gFin) ;
                }

              /********* if now identical, remove gh2 ****************/
              if (gh1->a1 == gh2->a1 && gh1->a2 == gh2->a2 &&
                  gh1->type == gh2->type)
                {
                  if (gh1->nerr < gh2->nerr) 
                    gh1->nerr = gh2->nerr ; /* notice that nerr is only used in bof cases */
                  goto nextj2 ;
                }
              /*******************************************************/
            }  /* end internal j1 loop */
          /* register */
        done:
          ghnew = arrayp (geneHits, newj++, HIT) ;
          *ghnew = *gh2  ;
          continue ;
        nextj2:  /* do not register this gh2 */
	  if (0 && upGene)  /* on A le bug le bug et j2 = 2 */
	    { 
	      arrayDestroy (geneHits) ;
	      invokeDebugger () ;
	    }
          continue ;  
        }  /* end j2 loop */
      arrayMax(geneHits) = newj ;
      arrayDestroy (oldHits) ;
    }  /* end j0 loop */
      if (0 && upGene)  /* on A le bug */
	{ 
	  arrayDestroy (geneHits) ;
	  invokeDebugger () ;
	}
  if (arrayMax(geneHits))
    {
      if (upGene) 
        cleanGeneReverseHits (geneHits) ;  
      cDNASwapA (geneHits) ;
      arraySort (geneHits, cDNAOrderGenesByA1) ;
    }
}

/*********************************************************************/
/*********************************************************************/

static void addOneGhostHits (int dnaMax, Array hits, int j1, int j2, Array walls)
{
  int i ;
  float xf = 0 ; 
  int cloneLength = 0 ;
  HIT *up ;
  KEY upEst = 0, downEst = 0, est = 0, est5p = 0, clone = 0, fuseTo = 0 ;
  int a1, a2, b1, b2, fuseMin, fuseMax ;
  OBJ Clone = 0 ;
  BOOL reverse = FALSE, downGene = TRUE ;
  
  a1 = a2 = b1 = b2 = 0 ;
  
  for (i = j1; i < j2 ; i++)
    {
      up = arrp(hits, i, HIT) ;
      
      if (! (up->type & gX) ||
          (!keyFindTag (up->cDNA_clone,  str2tag("Fuse_to_clone")) &&
           ( keyFindTag (up->cDNA_clone, _Mosaic) ||
             TRUE || /* we have had to many problems with this system */
             keyFindTag (up->cDNA_clone, str2tag("Length_anomaly")))))
        continue ;
      if (up->x1 < up->x2)
        { 
          if (!downEst)
            { a1 = up->a1 ; a2 = up->a2 ; }
          downEst = up->est ;
          if (!up->reverse) { est5p = up->est ; downGene = TRUE ; }
          else downGene = FALSE ;
          if (a1 > up->a1) a1 = up->a1 ;
          if (a2 < up->a2) a2 = up->a2 ;
          cloneLength -= up->a2 - up->a1 ;
          reverse = up->reverse ; 
        }
      else
        {
          if (!upEst)
            { b1 = up->a1 ; b2 = up->a2 ; }
          upEst = up->est ;
          if (up->reverse) { est5p = up->est ; downGene = FALSE ; }
          else  downGene = TRUE ;
          if (b1 > up->a1) b1 = up->a1 ;
          if (b2 < up->a2) b2 = up->a2 ;
          cloneLength -= up->a2 - up->a1 ; /* cumul both in same clone ! */
          reverse = up->reverse ; /* multiply defined but already rationalised */
        }
    }
  est = upEst ? upEst : downEst ;
  clone = keyGetKey (est, _cDNA_clone) ;
  if (!clone)
    return ;
  if (upEst && downEst && !keyGetKey (clone,str2tag("Fuse_to_clone") ) )
    return ;
  if (!upEst && !downEst)
    return ;
  if (est5p && keyFindTag (est5p, _PolyA_after_base) && 
      !keyFindTag (clone, _Internal_priming) &&
      !keyFindTag (clone, _Internal_priming_on_A_rich)
      )
    return ;
  if (!(Clone = bsCreate(clone)))
    return ;
  if (bsGetData (Clone, _PCR_product_size, _Float, &xf))
    cloneLength += 1000 * xf ;
  fuseMin = 10000000 ;
  fuseMax = -1000000 ;
  if (bsGetKey (Clone, str2tag("Fuse_to_clone"), &fuseTo))
    {
      if (downEst)
        for (i = 0; i < arrayMax(hits); i++)
          {
            up = arrp(hits, i, HIT) ;
            
            if (! (up->type & gX) ||
                up->cDNA_clone != fuseTo ||
                up->a1 < 0*a2)
              continue ;
            fuseMin = up->a1 ;
            fuseMax = up->a2 ;
            cloneLength = fuseMax - a2 ;
            if (cloneLength < 0)
              cloneLength = - cloneLength ;
            break ;
          }
      else
        for (i = arrayMax(hits) - 1 ; i >= 0 ; i--)
          {
            up = arrp(hits, i, HIT) ;
            
            if (! (up->type & gX) ||
                up->cDNA_clone != fuseTo ||
                0*up->a2 > b1)
              continue ;
            fuseMin = up->a1 ;
            fuseMax = up->a2 ;
            cloneLength = b1 - fuseMin ;
            if (cloneLength < 0)
              cloneLength = - cloneLength ;
            break ;
          }
    }
  bsDestroy (Clone) ;
  /* cloneLength  -= 350 ; feb 12, do not substract., there is always some introns */
  if (cloneLength < 0 || (fuseTo && fuseMax < fuseMin) )
    return ;
  
  if (walls) /* ghosts should not cut walls */
    {
      int iWall ; HIT *zw ;
      for (iWall = 0 ; iWall < arrayMax(walls) ; iWall++)
        {
          zw = arrayp (walls, iWall, HIT) ;
          if (zw->reverse == downGene)  /* opposite strands */
            continue ;
          if (downEst && a2 < zw->a1 + 50 && 
              a2 + cloneLength >  zw->a1 - 50)
            { cloneLength = -1 ; break ; }
          if (!downEst && b1 > zw->a1 - 50 && 
              b1 - cloneLength <  zw->a1 + 50)
            { cloneLength = -1 ; break ; }
        }
    }
  if (cloneLength < 0 || cloneLength > 100000) /* ignore huge ghosts */
    return ;
  
  /* sept 6, rather create just a long ghost */
  up = arrayp (hits, arrayMax(hits), HIT) ;
  up->reverse = reverse ;
  up->est = 1 ; /* ghost */
  up->cDNA_clone = clone ;
  up->type = gGhost ;

  if (downEst) /* create a ghost upEst */
    { 
      if (fuseTo)
        {
          up->a1 = fuseMin ; up->a2 = fuseMax ;
          up->x1 = 10000 ; up->x2 =  up->x1 +  fuseMax - fuseMin ; 
          up->type |= gFuseToGhost ;
        }
      else
        { 
          up->a1 = a2 + 50 ; up->a2 = up->a1 + cloneLength - 50 ;
          up->x1 =  cloneLength - 50 ; up->x2 = 1 ;
        }
    }
  else
    {
      if (fuseTo)
        {
          up->a1 = fuseMin ; up->a2 = fuseMax ;
          up->x1 = 1 ; up->x2 = fuseMax - fuseMin ;
          up->type |= gFuseToGhost ;
        }
      else
        { 
          up->a2 = b1 - 50 ; up->a1 = up->a2 - cloneLength + 50 ;
          up->x1 = 1 ; up->x2 = cloneLength - 50 ;
        }  
    }
}

/*********************************************************************/

static void flagTrackingCloneHits (Array hits, KEYSET trackingClones)
{
  int ii; 
  HIT *up ;
  
  if (hits && arrayMax(hits) && trackingClones && keySetMax(trackingClones))
    for (ii = 0, up = arrp (hits, 0, HIT) ; ii < arrayMax(hits)  ; ii++, up++)
      if (keySetFind (trackingClones, up->cDNA_clone, 0))
        up->type |= gDroppedGene ;
}

/*********************************************************************/

static void addGhostHits (int dnaMax, Array hits, Array walls)
{
  KEY old, clone ;
  int j, j1, jmax;
  HIT *up ;

  cDNASwapA (hits) ; 
  arraySort (hits, cDNAOrderByA1) ;

  jmax = arrayMax(hits)  ;
  for (old = 0, j = 0, j1 = 0 ; j < jmax ; j++)
    {
      up = arrp(hits, j, HIT) ;
      clone = up->cDNA_clone ;
      if (clone != old)
        {
          if (j > j1)
            addOneGhostHits (dnaMax, hits, j1, j, walls) ;
          j1 = j ;
        }
      up = arrp(hits, j, HIT) ; /* may have moved */
      old = up->cDNA_clone ;
    }
  if (j > j1)
    addOneGhostHits (dnaMax, hits, j1, j, walls) ; 
  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;
}

/*********************************************************************/

static void coalignAllIntrons (KEY cosmid, Array hits)
{
  int ii, jj ;
  HIT *hh, *hh2 ;
  KEY est = 0 ;
  Array dnaEst = 0 ;

  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;

  for (ii = 0, hh = arrp (hits, ii, HIT) ; ii < arrayMax (hits) ; hh++, ii++)
    { 
      if (hh->nerr < 5)
        continue ;
      if (dnaEst && est != hh->est)
        dnaEst = 0 ;
      est = hh->est ;
       if (! (hh->type & gD))
        { /* adapt beginning */
          for (jj = 0, hh2 = arrp (hits, jj, HIT) ; jj < arrayMax (hits) ; hh2++, jj++)
            if (hh2->est != est && (hh2->type & gD) && hh2->a2 == hh->a2)
              { /* trim excess */
                if (hh2->a1 < hh->a1 + 20 && hh2->a1 > hh->a2 - 3 && 
                    ( ii == 0 || (hh-1)->est != est))
                  {
                    int da = hh2->a1 - hh->a1 - 3 ;
                    hh->a1 += da ;
                    if (hh->x1 < hh->x2) 
                      hh->x1 += da ;
                    else
                      hh->x2 -= da ;                    
                  }
              }
        }

      if (! (hh->type & gF))
        { /* adapt end */
          for (jj = 0, hh2 = arrp (hits, jj, HIT) ; jj < arrayMax (hits) ; hh2++, jj++)
            if (hh2->est != est && (hh2->type & gF) && hh2->a1 == hh->a1)
              {
                /* find a model */
                if (hh2->a2 > hh->a2 - 20 && hh2->a2 < hh->a2 - 3 && 
                    ( ii == arrayMax(hits) || (hh+1)->est != est))
                  {
                    int da = hh->a2 - hh2->a2 - 3 ;

                    if (!dnaEst) dnaEst = getSeqDna (est) ;
                    if (dnaEst)
                      {

                      }
                    /* trim excess */
                    da -= 3 ;
                    hh->a2 -= da ;
                    if (hh->x1 < hh->x2) 
                      hh->x2 -= da ;
                    else
                      hh->x1 += da ;                    
                  }
              }
        }    
    }
}

/*********************************************************************/

static void checkAllIntrons (KEY cosmid, Array hits, Array dna, Array dnaR)
{
  int ii ;
  HIT *hh1, *hh2 ;
  Array dnaEst = 0 ;
  KEY est = 0 ;

  if (arrayMax(hits))
    {
      cDNASwapX (hits) ;
      arraySort (hits, cDNAOrderByX1) ;
      
      for (ii = 0, hh1 = arrp (hits, ii, HIT) ; ii < arrayMax (hits); hh1++, ii++)
        hh1->type &= gLink ; /* clean up needed because team jump abuses the same flags */
      for (ii = 0, hh1 = arrp (hits, ii, HIT), hh2 = hh1 + 1 ; ii < arrayMax (hits) - 1 ; hh1++, hh2++, ii++)
        { 
          if (hh1->est != est) dnaEst = 0 ;
          est = hh1->est ;
          if (!dnaEst) dnaEst = getSeqDna (est) ;
          if (!dnaEst)
            continue ;
          if (hh1->est == hh2->est && 
              checkIntron (cosmid, dna, dnaR, dnaEst, hh1, hh2, FALSE) == 2)
            { /* mark in A1 order */
              if (hh1->a1 < hh1->a2)
                { hh1->type |= gF ; hh2->type |= gD ; }
              else
                { hh1->type |= gD ; hh2->type |= gF ; }
            }
        }
    }
} /* checkAllIntrons  */

/*********************************************************************/

static Array getClonesWithBadIntron (KEY cosmid, Array hits, Array dna, Array dnaR)
{
  int ii, i1, nn, a1, a2, intronMax ;
  BOOL isBad ;
  HIT *hh1, *hh2, *vp ;
  Array introns = 0 , baddies = keySetCreate () ;
  KEY badClone = 0 ;

  if (arrayMax(hits))
    {
      cDNASwapA (hits) ; 
      arraySort (hits, cDNAOrderByA1) ;
      introns = teamCreateIntrons (hits, dna, dnaR) ; 
      if (introns &&  arrayMax (introns))
        {
          intronMax = arrayMax (introns) ;
          
          for (ii = nn = 0, hh1 = arrp (hits, ii, HIT), hh2 = hh1 + 1 ; ii < arrayMax (hits) - 1 ; hh1++, hh2++, ii++)
            { 
              if (hh1->cDNA_clone != badClone && hh1->est == hh2->est)
                {
                  a1 = hh1->a2 ; a2 = hh2->a1 ;
                  for (i1 = 0, isBad = TRUE, vp = arrp (introns, i1, HIT) ; 
                       isBad && i1 < intronMax ; vp++, i1++)
                    {
                      if (vp->a1 == a1 && vp->a2 == a2) /* classic intron, ok */
                        isBad = FALSE ;
                    }
                  if (isBad)
                    {
                      badClone = keySet (baddies, nn++) = hh1->cDNA_clone ;
                    }
                }
            }
	}
      arrayDestroy (introns) ;
    }
  keySetSort (baddies) ;
  keySetCompress (baddies) ;

  return baddies ;
} /* getClonesWithBadIntron */

/*********************************************************************/

static void avoidSuspectedInternalDeletions (Array geneHits)
{
  int ii, x ;
  HIT *hh ;
  KEY old = 0 ;
  OBJ Clone = 0 ;
  BOOL cleanUpNeeded = FALSE ; 

  cDNASwapX (geneHits) ;
  arraySort (geneHits, cDNAOrderByX1) ;


  for (ii = 0, hh = arrp (geneHits, 0, HIT) ; ii < arrayMax (geneHits) - 1 ; hh++, ii++)
    { 
      if (hh->zone == 9999)
        continue ;
      if (hh->cDNA_clone != old)
        {
          if (Clone)
            bsDestroy (Clone) ;
          old = hh->cDNA_clone ;
        }
      if (!Clone && keyGetKey (old,  _Suspected_internal_deletion) && !keyFindTag (old, str2tag("Manual_no_internal_deletion")))
        Clone = bsCreate (old) ;
      if (!Clone)
        continue ;

      /* flagged but the coordinates are missing */
      if (
          ( bsFindKey (Clone, _Suspected_internal_deletion, hh->est) &&
            ! bsGetData (Clone, _bsRight, _Int, &x) 
            )
          ||
          ! bsGetKey (Clone, _Suspected_internal_deletion, 0)
          )
        {   /* discard this est from the gene hits */
          int  j ;
          HIT *hh2 ;
          
          for (j = ii, hh2 = hh ; j < arrayMax (geneHits) ; hh2++, j++)
            if (hh2->est == hh->est)
              hh2->zone = 9999 ;
            else
              break ;
          continue ;
        }

      if (!(hh->type & (gI|gJ)))
        continue ;

      if (bsFindKey (Clone, _Suspected_internal_deletion, hh->est) && 
          bsGetData (Clone, _bsRight, _Int, &x) 
          )
        {
          /* gSuspect hits are discarded from cleanGeneHits
           * then in cdnaTr they are considered incompatible with
           * introns of at least 40 bp 
           * gGap are eaten by grignotte and fixed by stealFromNeighbours at the end 
           */
          do
            {
              if (hh->x1 >= x - 10 && hh->x1 <= x + 10)
                {
                  hh->type |= gSuspect ; 
                  hh->type |= gGap ;
                  hh->type |= gJ ; /* any gap is suspect */
                    
                  if (ii > 0 && (hh-1)->est == hh->est) 
                    (hh-1)->type &= ~gF ;
                  if (ii < arrayMax(geneHits) - 1&& (hh+1)->est == hh->est) 
                    (hh+1)->type &= ~gD ;
                }
            } while (bsGetData (Clone, _bsDown, _Int, &x)) ;
          if (FALSE &&
              (hh->type & gJ) &&
              ! (hh->type & gSuspect)
              )  /* we did not find these coordinates
                  * may be we should do something, 
                  * between apr 19 2005 and june 7, i was setting a gap for such introns
                  */
            { hh->type |= gSuspect ;  hh->type |= gGap ; }
        }
    }
  bsDestroy (Clone) ;

  if (cleanUpNeeded) /* keep happy few */
    {
      HIT *up, *vp ;
      int j, jj ;
      
      for (j = jj = 0, up = vp = arrp (geneHits, 0, HIT)  ; j < arrayMax(geneHits) ; j++, vp++)
        {
          if (vp->est)
            {
              if (jj < j) *up = *vp ;
              jj++ ; up++ ;
            }
        } 
    }

  return ;
}

/*********************************************************************/

static Array  getGeneHits (KEY cosmid, Array hits, Array clipTops, Array clipEnds,
                           Array walls, KEYSET trackingClones, Array dna, Array dnaR)
{
  Array geneHits = 0 ;
  KEY dummy, est = 0;
  unsigned int type ;
  int ii, jj, z1, z2, g1, clipTop, clipEnd, pA, xsl ;
  HIT *hh, *hh1, *gh=0, *lastConfirmedhh = 0, *lastConfirmedFuzzyhh = 0 ;
  Array dnaEst = 0 ;
  OBJ obj = 0 ;
  int intronFound ;
  BSMARK mark = 0 ;
  BOOL ignoreSL0 = FALSE ; /* was TRUE before feb 10 2005, see klp-19 */
  int minIntronSize = getPleaseMinIntronSize () ;

  cDNASwapX (hits) ;
  arraySort (hits, cDNAOrderByX1) ;

  for (ii = 0, hh = arrp (hits, ii, HIT) ; ii < arrayMax (hits) ; hh++, ii++)
    if (hh->type && hh->a1 > hh->a2)
      {
        /* do not trust checkAllIntrons which was done a la swapA */
        type = hh->type ;
        hh->type = 0 ;
        if (type & gD) hh->type |= gF ;
        if (type & gF) hh->type |= gD ;
        if (type & gLink) hh->type |= gLink ;
      }

  geneHits = arrayCreate (10, HIT) ; g1 = 0 ;
  for (jj = 0 ; jj < arrayMax (hits) ; jj++)
    { 
      hh = arrayp (hits, jj, HIT) ;
      hh1 = jj < arrayMax (hits) - 1 ? hh + 1 : 0 ;
      if (hh->est != est) dnaEst = 0 ;
      if (hh->a1 < 1)  /* i met this case, sept 01, yk36h8.5  a1=0 a2=115   x1=578 x2=693 */
        {
          if (hh->x1 < hh->x2) hh->x1 += hh->a1 - 1 ;
          else hh->x2 -= hh->a1 - 1 ;
          hh->a1 = 1 ;
        }
      est = hh->est ;
      /* discard only in cleangenehits, because genehits is passed to makemrna
      if (keySetFind (trackingClones, hh->cDNA_clone, 0))  // "XM_" 
        continue ;
        */
      if (!dna || 
          (!dnaR && (hh->a1 > hh->a2 || (hh+1)->a1 > (hh+1)->a2)))
        messcrash ("no dna in getGeneHits") ;
      /*dna = getSeqDna (cosmid) ;
      if (!dna) continue ;
      if (!dnaR && (hh->a1 > hh->a2 || (hh+1)->a1 > (hh+1)->a2))
        { dnaR = dnaCopy (dna) ; reverseComplement (dnaR) ; }
        */
      if (hh->est == 1) /* ghost */
        {
          gh = arrayp (geneHits, g1++, HIT) ;
          *gh = *hh ; /* a1 a2 x1 x2 est, reverse etc */
          gh->nerr = gh->nerrAll = 0 ;
          gh->ex1 = gh->ex2 = 0 ;
          gh->reverse = hh->a1 < hh->a2  ? hh->reverse : ! hh->reverse ;
          gh->type = gGhost ;
          continue ;
        }
      if (!dnaEst) dnaEst = getSeqDna (est) ;
      if (!dnaEst)
        continue ;
      if (!hh1 || hh->est != hh1->est) 
        goto lastExon ;
      intronFound = checkIntron (cosmid, dna, dnaR, dnaEst, hh, hh+1,TRUE);

      if (intronFound)
        {  /* register intron */
          if (hh->a2 > hh->a1 )
            { z1 = hh->a2 + 1 ; z2 =  hh1->a1 - 1 ;}
          else
            { z1 = hh->a2 - 1 ; z2 = hh1->a1 + 1 ; }
          
          gh = arrayp (geneHits, g1++, HIT) ;
          
          *gh = *hh ; /* est, reverse, nerr etc */
          gh->nerr = gh->nerrAll = 0 ;
          gh->ex1 = gh->ex2 = 0 ;
          gh->a1 = z1 ; gh->a2 = z2 ;
          gh->type = gI ;
          if (minIntronSize && 
              (
               (z1 < z2 && z2 - z1 + 1< minIntronSize) ||
               (z1 > z2 && z1 - z2 + 1< minIntronSize)
               )
              )
            gh->type |= gMicro ;
          if (intronFound == 1)
            gh->type |= gJ ;
          gh->x1 = hh->x2 ;  gh->x2 = hh->x2 + 1 ;
          if (hh->a1 < hh->a2)
            {    
              if (gh->a1 > gh->a2) 
                {  gh->a2 = gh->a1 + 1 ;  gh->type |= gJ  ; }
                /* should be a negative intron 
                { int aa = gh->a1 ; gh->a1 = gh->a2 ; gh->a2 = aa ; gh->type |= gN  ; }
                                    */
              gh->reverse = hh->reverse ;
            }
          else
            {
              if (gh->a1 < gh->a2)
                {  gh->a2 = gh->a1 - 1 ;  gh->type |= gJ  ; }
                /* should be a negative intron 
                { int aa = gh->a1 ; gh->a1 = gh->a2 ; gh->a2 = aa ; gh->type |= gN  ; }
                */
              gh->reverse = !hh->reverse ;
            }
        }
      /* register previous exon */
      gh = arrayp (geneHits, g1++, HIT) ;
      *gh = *hh ; /* a1 a2 x1 x2 est, reverse etc */
      gh->nerr = gh->nerrAll = 0 ;
      gh->ex1 = gh->ex2 = 0 ;
      gh->reverse = hh->a1 < hh->a2  ? hh->reverse : ! hh->reverse ;
      gh->type = gX ;
      if (hh->type & gLink) gh->type |= gLink ;
      if (hh->reverse) gh->type |= g3 ;
      else gh->type |= g5 ;
      if (hh == lastConfirmedhh)  /* Intron above Found */ 
        { 
          if (hh->reverse) gh->type |= gF ;
          else gh->type |= gD ;
        }
      else if (hh == lastConfirmedFuzzyhh)  /* Fuzzy Intron above Found */ 
        { 
          if (hh->reverse) gh->type |= gFF ;
          else gh->type |= gDF ;
        }
      if (intronFound == 2)     /* Intron below Found */ 
        {
          if (hh->reverse) gh->type |= gD ;
          else gh->type |= gF ;
        }
      else if (intronFound == 1)     /* Fuzzy Intron below Found */ 
        {
          if (hh->reverse) gh->type |= gDF ;
          else gh->type |= gFF ;
        }

      obj = bsCreate (hh->est) ;
      clipTop =  keySet(clipTops, KEYKEY(hh->est)) ;
      clipEnd =  keySet(clipEnds, KEYKEY(hh->est)) ;
              if (obj && !hh->reverse &&
                  bsGetKey (obj, _Transpliced_to, &dummy))
                do
                  {
                    mark = bsMark (obj, mark) ;
                    if (bsGetData (obj, _bsRight, _Int, &xsl))
                      {
                        if (!strncmp(name(dummy),"SL", 2) &&
                            clipTop == xsl + 1 &&
                            clipTop >= hh->x1-3)
                          {
                            if (clipTop == hh->x1 && strncmp(name(dummy),"SL0", 3) &&
                                checkSL (cosmid, dna, dnaR, dnaEst, hh, TRUE)) /* treat fuzzy SL1 as SL0 */
                              { gh->type |= (gD | gS) ; gh->type &= ~gS0 ; }
                            else if (! (gh->type & gS))
                              {
                                if (ignoreSL0)
                                  gh->type |= (gS0) ;
                                else
                                  gh->type |= (gDF | gS0) ;
                              }
                          }
                      }
                    bsGoto (obj, mark) ;
                  } while (bsGetKey (obj, _bsDown, &dummy)) ;
              
              if (obj && hh->reverse &&
                  bsGetKey (obj, _Transpliced_to, &dummy))
                do
                  {
                    mark = bsMark (obj, mark) ;
                    if (bsGetData (obj, _bsRight, _Int, &xsl))
                      {
                        if (!strncmp(name(dummy),"SL", 2) &&
                            clipEnd == xsl - 1 &&
                            clipEnd <= hh->x2 + 2)
                          {
                            if (clipEnd == hh->x2 && strncmp(name(dummy),"SL0", 3) &&
                                checkSL (cosmid, dna, dnaR, dnaEst, hh, FALSE)) /* treat fuzzy SL1 as SL0 */
                              { gh->type |= (gD | gS) ; gh->type &= ~gS0 ; }
                            else if (! (gh->type & gS))
                              {
                                if (ignoreSL0)
                                  gh->type |= (gS0) ;
                                else
                                  gh->type |= (gDF | gS0) ;
                              }
                          }
                      }
                    bsGoto (obj, mark) ;
                  } while (bsGetKey (obj, _bsDown, &dummy)) ;

              if (obj && 
                  1 &&   /* ignore danielle 2006_10_07 */
                  (
                   (!hh->reverse && clipTop >= hh->x1 - 5) ||
                   (hh->reverse && clipEnd <= hh->x2 + 5)
                  )
                  &&
                  bsFindTag (obj, str2tag("Complete_CDS")) &&
                  !bsFindTag (obj, str2tag("Is_partial")) 
                  )
                gh->type |= (gCompleteCDS) ;

              if (obj && !hh->reverse &&
                  !bsFindTag (obj, _Problem) &&
                  (        /* yk522a6 cas complique */           
                   bsFindTag (obj, str2tag("Real_5prime")) 
                   ) &&
                  clipTop >= hh->x1 - 5)
                gh->type |= (gD | gS | gReal5p) ;
              
              if (obj && hh->reverse &&
                  !bsFindTag (obj, _Problem) &&
                  (                  
                   bsFindTag (obj, str2tag("Real_5prime")) 
                   ) &&
                  clipEnd <= hh->x2 + 5)
                gh->type |= (gD | gS | gReal5p) ;
              
              if (obj && !hh->reverse &&
                  !bsFindTag (obj, _Problem) &&
                  (
                   !keyFindTag (hh->cDNA_clone, _Internal_priming) && /* retabli sept 22 2005, dan */
                   /* !keyFindTag (hh->cDNA_clone, _Internal_priming_on_A_rich) && supprime feb 7 2007 */
                   bsFindTag (obj, str2tag("Real_3prime"))
                   ) &&
                  clipEnd <= hh->x2 + 5)
                gh->type |= (gA | gReal3p) ;
              
              if (obj && hh->reverse &&
                  !bsFindTag (obj, _Problem) &&        
                  (
                   !keyFindTag (hh->cDNA_clone, _Internal_priming) && /* retabli sept 22 2005, dan */
                   /* !keyFindTag (hh->cDNA_clone, _Internal_priming_on_A_rich) && supprime feb 7 2007 */
                   bsFindTag (obj, str2tag("Real_3prime")) 
                   ) &&
                  clipTop >= hh->x1 - 5)
                gh->type |= (gA | gReal3p) ;
      
              if (obj && hh->reverse &&
                  !keyFindTag (hh->cDNA_clone, _Internal_priming) &&
                  !keyFindTag (hh->cDNA_clone, _Internal_priming_on_A_rich) &&
                  !bsFindTag (obj, _Problem) &&
                  bsGetData (obj, _PolyA_after_base, _Int, &pA) &&
                  pA < 100 &&
		  clipTop < pA + 10 &&
                  clipTop >= hh->x1 - 5)
                gh->type |= gA ; /* (gF | gA) ; */
              
              if (obj && !hh->reverse &&
                  !keyFindTag (hh->cDNA_clone, _Internal_priming) &&
                  !keyFindTag (hh->cDNA_clone, _Internal_priming_on_A_rich) &&
                  !bsFindTag (obj, _Problem) && 
                  bsGetData (obj, _PolyA_after_base, _Int, &pA) &&
                  (pA > 200 || 2*pA > hh->clipEnd) &&
		  clipEnd > pA - 10 &&
                  clipEnd <= gh->x2 + 5)  /* use gh->x2, which may now be != hh->x2 */
                gh->type |= gA ; /* (gF | gA) ; */
              
              if (hh->reverse &&  /* 3prime read primed on polyA */
                  hh->x1 < 30 &&  
                  !keyFindTag (hh->cDNA_clone, _Internal_priming) &&
                  !keyFindTag (hh->cDNA_clone, _Internal_priming_on_A_rich) &&
                  !keyFindTag (hh->cDNA_clone, _Mosaic) &&
                  !bsFindTag (obj, _Problem) &&
                  hh->cDNA_clone && keyFindTag (hh->cDNA_clone, _Primed_on_polyA))
                { 
                   int dx = .007 * (hh->x1 - clipTop) ;
                   gh->a2 += dx ; hh->x1 -= dx ;
                   gh->type |= gA ;  
                }
              /* check exon may move exon boundaries */
     checkExon (cosmid, dna, dnaR, dnaEst, gh, 
                !(hh->type & gD) && (hh == lastConfirmedFuzzyhh), 
                !(hh->type & gF) && (intronFound == 1)) ; 
     bsDestroy (obj) ;



      if (intronFound == 2)
        lastConfirmedhh = hh1 ;
      if (intronFound == 1)
        lastConfirmedFuzzyhh = hh1 ;
      continue ;
  
    lastExon:   /* special case */
      if (hh == lastConfirmedhh || hh->x1 <= hh->clipTop + 10 || hh->nerr < 10 ||
          hh->x2 > hh->clipEnd - 10 || hh->x2 > hh->x1 + 30 ||
          jj == 0 || (jj > 0 && hh->est != (hh-1)->est)
          )
        {
          gh = arrayp (geneHits, g1++, HIT) ; 
          *gh = *hh ; /* a1 a2 x1 x2 est, reverse etc */
          gh->nerr = gh->nerrAll = 0 ;
          gh->ex1 = gh->ex2 = 0 ;
          gh->reverse = hh->a1 < hh->a2  ? hh->reverse : ! hh->reverse ;
          gh->type = gX ;
          if (hh->type & gLink) gh->type |= gLink ;
          if (hh->reverse) gh->type |= g3 ;
          else gh->type |= g5 ;
          if (hh == lastConfirmedhh)  /* Intron above Found */ 
            { 
              if (hh->reverse) gh->type |= gF ;
              else gh->type |= gD ;
            }
          if (hh == lastConfirmedhh)  /* Intron above Found */ 
            { 
              if (hh->reverse) gh->type |= gF ;
              else gh->type |= gD ;
            }
          else if (hh == lastConfirmedFuzzyhh)  /* Fuzzy Intron above Found */ 
            { 
              if (hh->reverse) gh->type |= gFF ;
              else gh->type |= gDF ;
            }

            {
              obj = bsCreate (hh->est) ;
              clipTop =  keySet(clipTops, KEYKEY(hh->est)) ;
              clipEnd =  keySet(clipEnds, KEYKEY(hh->est)) ;
              if (obj && !hh->reverse &&
                  bsGetKey (obj, _Transpliced_to, &dummy))
                do
                  {
                    mark = bsMark (obj, mark) ;
                    if (bsGetData (obj, _bsRight, _Int, &xsl))
                      {
                        if (!strncmp(name(dummy),"SL", 2) &&
                            clipTop == xsl+1 &&
                            clipTop >= hh->x1 - 3)
                          {
                            if (clipTop == hh->x1 && strncmp(name(dummy),"SL0", 3) &&
                                checkSL (cosmid, dna, dnaR, dnaEst, hh, TRUE)) /* treat fuzzy SL1 as SL0 */
                              { gh->type |= (gD | gS) ; gh->type &= ~gS0 ; }
                            else if (! (gh->type & gS))
                              {
                                if (ignoreSL0)
                                  gh->type |= (gS0) ;
                                else
                                  gh->type |= (gDF | gS0) ;
                              }
                          }
                      }
                    bsGoto (obj, mark) ;
                  } while (bsGetKey (obj, _bsDown, &dummy)) ;
              
              if (obj && hh->reverse &&
                  bsGetKey (obj, _Transpliced_to, &dummy))
                do
                  {
                    mark = bsMark (obj, mark) ;
                    if (bsGetData (obj, _bsRight, _Int, &xsl))
                      {
                        if (!strncmp(name(dummy),"SL", 2) &&
                            clipEnd == xsl - 1 &&
                            clipEnd <= hh->x2 + 2)
                          {
                            if (clipEnd == hh->x2 && strncmp(name(dummy),"SL0", 3) &&
                                checkSL (cosmid, dna, dnaR, dnaEst, hh, FALSE)) /* treat fuzzy SL1 as SL0 */
                              { gh->type |= (gD | gS) ; gh->type &= ~gS0 ; }
                            else if (! (gh->type & gS))
                              {
                                if (ignoreSL0)
                                  gh->type |= (gS0) ;
                                else
                                  gh->type |= (gDF | gS0) ;
                              }
                          }
                      }
                    bsGoto (obj, mark) ;
                  } while (bsGetKey (obj, _bsDown, &dummy)) ;
              
              if (obj && 
                  1 &&   /* ignore danielle 2006_10_07 */
                  (
                   (!hh->reverse && clipTop >= hh->x1 - 5) ||
                   (hh->reverse && clipEnd <= hh->x2 + 5)
                  )
                  &&
                  bsFindTag (obj, str2tag("Complete_CDS")) &&
                  !bsFindTag (obj, str2tag("Is_partial"))
                  ) 
                gh->type |= (gCompleteCDS) ;

              if (obj && !hh->reverse &&
                  !bsFindTag (obj, _Problem) &&
                  (
                  
                    bsFindTag (obj, str2tag("Real_5prime")) 
                   ) &&
                  clipTop >= hh->x1 - 5)
                gh->type |= (gD | gS | gReal5p) ;
              
              if (obj && hh->reverse &&
                  !bsFindTag (obj, _Problem) &&
                  (
                   
                   bsFindTag (obj, str2tag("Real_5prime")) 
                   ) &&
                  clipEnd <= hh->x2 + 5)
                gh->type |= (gD | gS | gReal5p) ;
              
              if (obj && !hh->reverse &&
                  !bsFindTag (obj, _Problem) &&
                  (
                   !keyFindTag (hh->cDNA_clone, _Internal_priming) && /* retabli sept 22 2005, dan */
                   /* !keyFindTag (hh->cDNA_clone, _Internal_priming_on_A_rich) && supprime feb 7 2007 */
                   bsFindTag (obj, str2tag("Real_3prime"))
                   ) &&
                  clipEnd <= hh->x2 + 5)
                gh->type |= (gA | gReal3p) ;
              
              if (obj && hh->reverse &&
                  !bsFindTag (obj, _Problem) &&        
                  (
                   !keyFindTag (hh->cDNA_clone, _Internal_priming) && /* retabli sept 22 2005, dan */
                   /* !keyFindTag (hh->cDNA_clone, _Internal_priming_on_A_rich) && supprime feb 7 2007 */
                   bsFindTag (obj, str2tag("Real_3prime")) 
                   ) &&
                  clipTop >= hh->x1 - 5)
                gh->type |= (gA | gReal3p) ;

              if (obj && hh->reverse &&
                  !keyFindTag (hh->cDNA_clone, _Internal_priming) &&
                  !keyFindTag (hh->cDNA_clone, _Internal_priming_on_A_rich) &&
                  !bsFindTag (obj, _Problem) &&
                  bsGetData (obj, _PolyA_after_base, _Int, &pA) &&
                  pA < 100 &&
		  clipTop < pA + 10 &&
                  clipTop >= hh->x1 - 5)
                gh->type |= gA ; /* (gF | gA) ; */
              
              if (obj && !hh->reverse &&
                  !keyFindTag (hh->cDNA_clone, _Internal_priming) &&
                  !keyFindTag (hh->cDNA_clone, _Internal_priming_on_A_rich) &&
                  !bsFindTag (obj, _Problem) &&
                  bsGetData (obj, _PolyA_after_base, _Int, &pA) &&
                  (pA > 200 || 2*pA > hh->clipEnd) &&
		  clipEnd > pA - 10 &&
                  clipEnd <= gh->x2 + 5)  /* use gh->x2, which may now be != hh->x2 */
                gh->type |= gA ; /* (gF | gA) ; */
              
              if (hh->reverse &&  /* 3prime read primed on polyA */
                  hh->x1 < 30 &&  /* on this lib there is NO vector, do not use cliptop */
                  !keyFindTag (hh->cDNA_clone, _Internal_priming) &&
                  !keyFindTag (hh->cDNA_clone, _Internal_priming_on_A_rich) &&
                  !keyFindTag (hh->cDNA_clone, _Mosaic) &&
                  !bsFindTag (obj, _Problem) &&
                  hh->cDNA_clone && keyFindTag (hh->cDNA_clone, _Primed_on_polyA))
                { 
                   int dx = .007 * (hh->x1 - clipTop) ;
                   gh->a2 += dx ; hh->x1 -= dx ;
                   gh->type |= gA ;  
                }
              checkExon (cosmid, dna, dnaR, dnaEst, gh,
                         !(hh->type & gD) && (hh == lastConfirmedFuzzyhh), 
                         FALSE) ;  
              bsDestroy (obj) ;
            }
               
        }
      else if (hh == lastConfirmedFuzzyhh)  /* delete that intron */
        {
          *(gh - 1) = *gh ; arrayMax(geneHits)-- ; g1-- ;
        }
      lastConfirmedhh = 0 ;
    }
  avoidSuspectedInternalDeletions (geneHits) ; /* after gI,gJ,gMicro are established */
  cDNASwapA (geneHits) ;
  arraySort (geneHits, cDNAOrderByA1) ; /* needed by diamino etc */

  messfree (mark) ;
  return geneHits ;
} /* getGeneHits */

/*********************************************************************/
/*********************************************************************/
/******* Attribute clones to transcripts *****************************/

/*********************************************************************/
/*********************************************************************/

/* allSpl, allDna :: one per poossible product
   spl:  KEY int int int int 
     tag (intron/gap/exon..)
     begin+end of this exon in (product) dna
     number of this exon in original
     */

/*********************************************************************/
typedef struct { KEY gene, type ; 
                int a1, a2, x1, x2, exonNb, iGap ; 
  BOOL isExon, isIntron, isGap, isOrfGap, isFirstExon, isLastExon, isBack ; } SPL ;
/* a1 a2 are in genomic, x1 x2 in mrna in construction */

/*********************************************************************/

static int splOrder (const void *va, const void *vb)
{
  const SPL *up = (const SPL *)va, *vp = (const SPL *)vb ;

  if (up->a1 != vp->a1)
    return up->a1 - vp->a1 ;
  else if (up->a2 != vp->a2)
    return up->a2 - vp->a2 ;      
  return 0 ;
} /* splOrder */

/********************************************************************************/
/* clean up all data associated to the gene */
void cdnaCleanUp (KEY cosmid, KEY tg, KEYSET tgs)
{
  KEYSET 
    genes = 0, xs = 0, xs1 = 0, xs2 = 0 ;
  KEY geneBox ;
  OBJ Gene = 0 ;
  int i ;

  mrnaCleanUp (cosmid, tg, tgs) ;

  if (cosmid)
    {
      xs = queryKey (cosmid, ">Transcribed_gene") ;
      xs1 = queryKey (cosmid, "{>MRNA} SETOR {>Transcribed_gene;>mRNA}") ;
    }
  else if (tgs)
    {
      xs = keySetCopy (tgs) ;
      xs1 = query (tgs, ">MRNA") ;
    }
  else if (tg)
    {
      xs = keySetCreate () ;
      keySet(xs, 0) = tg ;
      xs1 = queryKey (tg, ">mrna") ;
    }
  xs2 = query (xs1,  "{>DNA} SETOR {>Product} SETOR {>Product;>Peptide}") ;

  genes = query (xs, ">gene") ;
  i = keySetMax (genes) ;
  while (i--)
    {
      geneBox = keySet (genes, i) ;
      if (geneBox && (Gene = bsUpdate (geneBox)))
        {
          if (bsFindTag (Gene, str2tag ("Structure")))
            bsRemove (Gene) ;
          if (bsFindTag (Gene, str2tag ("Expression_profile")))
            bsRemove (Gene) ;
          if (bsFindTag (Gene, str2tag ("Title")))
            bsRemove (Gene) ;
          bsSave (Gene) ;
        }
    }
  
  keySetKill(xs) ;
  keySetKill(xs1) ;
  keySetKill(xs2) ;

  keySetDestroy (genes) ;
  keySetDestroy (xs) ;
  keySetDestroy (xs1) ;
  keySetDestroy (xs2) ;

  return ;
} /* cdnaCleanUp */

/*********************************************************************/
/*********************************************************************/

static BOOL teamInitialJump (KEY cosmid, Array introns, Array oldHits, Array newHits, int ii, int *jjp, Array dna, Array dnaR)
{
  HIT *up, *vp ;
  KEY est ;
  int aa, a1, a2, xx, i1, intronMax = arrayMax(introns) ;
  Array dnaEst = 0 ;
  char *cp, *cq, *cpMin, *cpMax ;

  up = arrp (oldHits, ii, HIT) ;
  if (!up->est) return 0 ;
  est = up->est ; a1 = up->a1 ; a2 = up->a2 ;
  
  for (i1 = 0 ; i1 < intronMax ; i1++)
    {
      vp = arrp (introns, i1, HIT) ;
      if (! vp->x1) /* only regularize towards good introns */
        continue ; 
      if (up->nerr == 0 && (vp->a2 <= up->a1 || vp->a2 > up->a1 + 8)) continue ;
      if (vp->a2 >= a1 && vp->a2 < a1 + 50 &&  vp->a2 < a2) /* intron ends in reasonable place */
        {
          if (!dnaEst)
            dnaEst = getSeqDna (est) ;
          if (!dnaEst)
            return FALSE ;
          cpMin = arrp (dnaEst, up->clipTop - 1, char) ;
          cpMax = arrp (dnaEst, up->clipEnd - 1, char) ;
          cp = arrp (dnaEst, up->x2 - 1, char) ; /* last base of my exon */
          if (up->x1 < up->x2)
            {
              aa = up->a2 ; xx = up->x2 ;
              cq = arrp (dna, aa - 1, char) ;
              while (cp > cpMin && (*cp & *cq)) { cp-- ; cq-- ; aa-- ;xx-- ; }
            }
          else
            {
              aa = up->a2 ; xx = up->x2 ;
              cq = arrp (dna, aa - 1, char) ;
              while (cp < cpMax && complementBase(*cp & 255) & *cq) { cp++ ; cq-- ; aa-- ; xx++ ; }
            }
          if (aa < vp->a2) /* so i was perfect up and including vp->a2 */
            {
              int i = 0, iend = 0 ;
              if (up->x1 < up->x2)
                xx = up->x2 - (up->a2 - vp->a2) - 1 ;
              else
                xx = up->x2 + (up->a2 - vp->a2) + 1 ;
              aa = vp->a1 ;
              cp = arrp(dnaEst, xx - 1, char) ;
              cq = arrp (dna, aa - 1, char) ;
              for (i = 0 ; i < 64  ; i++)
                if (up->x1 < up->x2)
                  {
                    if (cp < cpMin) { if (vp->x1) /* gt_ag */ iend = 10 ; break ; }
                    if (*cp & *cq) { cp-- ; cq-- ; }
                    else break ;
                  }
                else
                  { 
                    if (cp > cpMax) { iend = 10 ;  break ; }
                    if (complementBase(*cp & 255) & *cq) { cp++ ; cq-- ; }
                    else break ;
                  }
              if (i + iend > 6)
                {
                  HIT *newUp ;
                  
                  newUp = arrayp (newHits, (*jjp - 1), HIT) ;
                  *newUp = *up ;
                  if (i > 1)
                    {
                      if (up->x1 < up->x2)
                        {
                          newUp->x1 = up->x2 - (up->a2 - vp->a2) - i ;
                          newUp->x2 = up->x2 - (up->a2 - vp->a2) - 1 ;
                          newUp->a1 = vp->a1 - i + 1 ;
                          newUp->a2 = vp->a1 ;
                        }
                      else
                        {
                          newUp->x1 = up->x2 + (up->a2 - vp->a2) + i ;
                          newUp->x2 = up->x2 + (up->a2 - vp->a2) + 1 ;
                          newUp->a1 = vp->a1 - i + 1 ;
                          newUp->a2 = vp->a1 ;
                        }
                      newUp->type |= gF ;
                      newUp->nerr = 0 ; /* prevent teamCentralJump */
                      
                      newUp = arrayp (newHits, (*jjp)++, HIT) ;
                      *newUp = *up ;
                    }
                  if (up->x1 < up->x2)
                    {
                      newUp->x1 = up->x2 - (up->a2 - vp->a2) ;
                      newUp->a1 = vp->a2 ;
                    }
                  else
                    {
                      newUp->x1 = up->x2 + (up->a2 - vp->a2) ;
                      newUp->a1 = vp->a2 ;
                    } 
                  newUp->nerr = 0 ; /* prevent teamCentralJump */
                  newUp->type |= gD ;
                  return TRUE ;
                }
              
            }
        }
    }
  if (ii < arrayMax (oldHits) - 1 && up->est == (up+1)->est &&
      up->a2 - up->a1 - 2 * up->nerr < 20) /* do not accept a short non classic intron */
    {
      a1 = up->a2 ; a2 = (up+1)->a1 ;
      for (i1 = 0 ; i1 < intronMax ; i1++)
        {
          vp = arrp (introns, i1, HIT) ;
          if (vp->a1 == a1 && vp->a2 == a2) /* classic intron, ok */
            return FALSE ;
        }
      (*jjp)-- ; /* kill this hit */
      return TRUE ;
    }
  /* look for an intron ending close to my a1 */
  /* do not destroy dnaEst */
  return FALSE ;
} /* teamInitialJump */

/*********************************************************************/

static int teamFinalJump (KEY cosmid, Array introns, Array oldHits, Array newHits, int ii, int *jjp, Array dna, Array dnaR)
{
  HIT *up, *vp ;
  KEY est ;
  int aa, a1, a2, xx, i1, intronMax = arrayMax(introns) ;
  Array dnaEst = 0 ;
  char *cp, *cq, *cpMin, *cpMax ;

  up = arrp (oldHits, ii, HIT) ;
  if (!up->est) return 0 ;
  est = up->est ; a1 = up->a1 ; a2 = up->a2 ;
    
  for (i1 = 0 ; i1 < intronMax ; i1++)
    {
      vp = arrp (introns, i1, HIT) ;
      if (up->nerr == 0 && (vp->a1 >= up->a2 || vp->a1 < up->a2 - 8)) continue ;
      if (! vp->x1) /* only regularize towards good introns */
        continue ; 
      if (vp->a1 < a2 && vp->a1 > a2 - 50 && vp->a1 > a1) /* intron ends in reasonable place */
        {
          if (!dnaEst)
            dnaEst = getSeqDna (est) ;
          if (!dnaEst)
            return FALSE ;
          cpMin = arrp (dnaEst, up->clipTop - 1, char) ;
          cpMax = arrp (dnaEst, up->clipEnd - 1, char) ;
          cp = arrp (dnaEst, up->x1 - 1, char) ; /* first base of my exon */
          if (up->x1 < up->x2)
            {
              aa = up->a1 ; xx = up->x1 ;
              cq = arrp (dna, aa - 1, char) ;
              while (*cp & *cq) { cp++ ; cq++ ; aa++ ;xx++ ; }
            }
          else
            {
              aa = up->a1 ; xx = up->x1 ;
              cq = arrp (dna, aa - 1, char) ;
              while (xx > 0 && complementBase(*cp & 255) & *cq) { cp-- ; cq++ ; aa++ ; xx-- ; }
            }
          if (aa > vp->a1) /* so i was perfect up and including vp->a1 */
            {
              int i = 0, iend = 0 ;
              if (up->x1 < up->x2)
                xx = up->x1 + (vp->a1 - up->a1) + 1 ;
              else
                xx = up->x1 - (vp->a1 - up->a1) - 1 ;
              aa = vp->a2 ;
              cp = arrp(dnaEst, xx - 1, char) ;
              cq = arrp (dna, aa - 1, char) ;
              for (i = 0 ; i < 64  ; i++)
                if (up->x1 < up->x2)
                  {
                    if (cp > cpMax) { iend = 10 ; break ; }
                    if (*cp & *cq) { cp++ ; cq++ ; }
                    else break ;
                  }
                else
                  {
                    if (cp < cpMin) { iend = 10 ; break ; }
                    if (complementBase(*cp & 255) & *cq) { cp-- ; cq++ ; }
                    else break ;
                  }
              if (i + iend > 6)
                {
                  HIT *newUp ;
                  int ddx = vp->a1 - up->a1 ;

                  newUp = arrayp (newHits, (*jjp - 1), HIT) ;
                  *newUp = *up ;
                  if (up->x1 < up->x2)
                    {
                      newUp->x2 = up->x1 + ddx ;
                      newUp->a2 = vp->a1 ;
                    }
                  else
                    {
                      newUp->x2 = up->x1 - ddx ;
                      newUp->a2 = vp->a1 ;
                    }         
                  newUp->type |= gF ;

                  if (i > 1)  /* so that orientation is kept */
                    {
                      newUp->nerr = 0 ; /* prevent teamCentralJump */
                      newUp = arrayp (newHits, (*jjp)++, HIT) ;
                      *newUp = *up ;
                      if (up->x1 < up->x2)
                        {
                          newUp->x2 = up->x1 + ddx + i ;
                          newUp->x1 = up->x1 + ddx + 1 ;
                          newUp->a2 = vp->a2 + i - 1 ;
                          newUp->a1 = vp->a2 ;
                        }
                      else
                        {
                          newUp->x2 = up->x1 - ddx - i ;
                          newUp->x1 = up->x1 - ddx - 1 ;
                          newUp->a2 = vp->a2 + i - 1 ;
                          newUp->a1 = vp->a2 ;
                        } 
                      newUp->nerr = 0 ; /* prevent teamCentralJump */
                      newUp->type |= gD ;
                    }
                  return 1 ;
                }
              
            }
        }
    }

  if (ii > 0 && up->est == (up-1)->est &&
      up->a2 - up->a1 - 2 * up->nerr < 20) /* do not accept a short non classic intron */
    {
      a1 = (up-1)->a2 ; a2 = up->a1 ;
      for (i1 = 0 ; i1 < intronMax ; i1++)
        {
          vp = arrp (introns, i1, HIT) ;
          if (vp->a1 == a1 && vp->a2 == a2) /* classic intron, ok */
            return 0 ;
        }
      (*jjp)-- ; /* kill this hit */
      return 2 ; /* please reanalyse previous exon */
    }

  /* look for an intron ending close to my a1 */
  /* do not destroy dnaEst */
  return 0 ;
}   /* teamFinalJump */

/*********************************************************************/

static BOOL teamDoubleJump (KEY cosmid, Array introns, HIT *up, HIT *xp, HIT *wp,
                            Array dna, Array dnaR)
{
  HIT *vp ;
  int a1, a2, b1, b2, c1, c2, bb, da, db, dc, ddb, ddb1, ddb2 ;
  int ii, nn, intronMax = arrayMax (introns) ;
  
  a1 = up->a1 ; a2 = up->a2 ; 
  b1 = xp->a1 ; b2 = xp->a2 ;  /* short central exon */
  c1 = wp->a1 ; c2 = wp->a2 ;

  /* check if both sides are gt_ag || gc_ag */
  for (ii = nn = 0 ; ii < intronMax ; ii++)
    {
      vp = arrp (introns, ii, HIT) ;
      if (vp->a1 == a2 && vp->a2 == b1 && vp->x1)
        nn++ ;
      if (vp->a1 == b2 && vp->a2 == c1 && vp->x1)
        nn++ ;
    }
  if (nn == 2)  /* we are happy with this nice small exon */
    return FALSE ;
  /* may be we could kill this exon
   * we look for an intron that will gobble it
   */
  bb = c2 - c1 + b2 - b1 + a2 - a1 + 3 ; /* total length of the 3 exons */
  for (ii = nn = 0 ; ii < intronMax ; ii++)
    {
      vp = arrp (introns, ii, HIT) ;
      if (! vp->x1) /* only consider good introns */
        continue ; 
      if (vp->a1 < a1 || vp->a2 > c2 || vp->a1 < a2 - 3 || vp->a2 > c1 + 3)
        continue ;
      nn = vp->a1 - a1 + c2 - vp->a2 + 2 ;
      if (nn < bb - 3 || nn > bb + 3)
        continue ;
      db = b2 - b1 + 1 ;
      da = vp->a1 - a2 ;
      dc = c1 - vp->a2 ;
      ddb = db - da - dc ;
      ddb1 = ddb/2 ; ddb2 = ddb - ddb1 ;
      
      if (xp->x1 < xp->x2)
        { up->x2 += da + ddb1 ; wp->x1 -= dc + ddb2 ; }
      else
        { up->x2 -= da + ddb1 ; wp->x1 += dc + ddb2 ; }
      up->a2 += da ; wp->a1 -= dc ;
      xp->est = 0 ;
      
      return TRUE ;
    }

  return FALSE ;    
} /* teamDoubleJump */

/*********************************************************************/

static BOOL teamFirstIntron (KEY cosmid, Array introns,  HIT *up, HIT *wp, Array dna, Array dnaR)
{
  HIT *vp, *xp, *yp ;
  int a2, b1, i1, nerr, newnerr, derr, da, di, intronMax = arrayMax(introns), nn, nnmax ;
  BOOL foundJump = FALSE, gt_ag = FALSE ;
  static Array myHits = 0 ;
  
  if (!myHits)
    myHits = arrayCreate (2, HIT) ;
  
  /* est = up->est ;  a1 = up->a1 ; */ a2 = up->a2 ; b1 = wp->a1 ; /* b2 = wp->a2 ; */
  nerr = up->nerr + wp->nerr ;

  for (i1 = nn = nnmax = 0 ; i1 < intronMax ; i1++) /* notice that i may try several solutions */
    {
      vp = arrp (introns, i1, HIT) ;
      if (vp->a1 == a2 && vp->a2 == b1) /* self */
        { 
          if (!vp->x1) /* not good intron */
	    nerr += 3 ;
	  else
	    gt_ag = TRUE ;
          break ;
        }
    }

  for (i1 = nn = nnmax = 0 ; i1 <= intronMax ; i1++) /* notice that i may try several solutions */
    {
      if (i1 < intronMax)
	{
	  vp = arrp (introns, i1, HIT) ;
	  if (! (vp->x1 && /* only regularize towards good introns */
		 vp->a2 >= b1 && vp->a2 < b1 + 8 && vp->a1 > a2)) /* closer jump */         
	    continue ;
	  di = vp->a1 - a2 ; /* intron shortening */
	  derr = di > 3000 ? 2 : (di > 1000 ? 1 : 0) ;
	  xp = arrayp (myHits, 0, HIT) ;  *xp = *up ;
	  yp = arrayp (myHits, 1, HIT) ;  *yp = *wp ;
	  xp->a2 = vp->a1 ; xp->a1 = vp->a1 - (up->a2 - up->a1) ;
	  da = vp->a2 - b1 ;
	  if (da)
	    {
	      yp->a1 += da ; 
	      xp->a1 -= da ;
	      if (up->x1 < up->x2) { xp->x2 += da ; yp->x1 += da ; }
	      else { xp->x2 -= da ; yp->x1 -= da ; }
	    }
	}
      else /* try locally */
	{
	  xp = arrayp (myHits, 0, HIT) ; *xp = *up ;
	  yp = arrayp (myHits, 1, HIT) ; 
	  xp->a2 = wp->a1 - 1 ;
	  xp->a1 = xp->a2 - (up->a2 - up->a1) ;
	  nerr = up->nerr + (gt_ag ? 0 : 3) ;
	  yp->nerr = 0 ; 
	  derr = gt_ag ? 0 : +3 ; /* i like to kill non gt_ag */
	  arrayMax (myHits) = 1 ;
	}
      countErrors(cosmid, myHits, dna, dnaR, 0, FALSE) ; /* implies cDNASwapX */
      newnerr = xp->nerr + yp->nerr ;
      if (newnerr <= nerr + derr && xp->maxExact > up->maxExact - 10)
	{ 
	  foundJump = TRUE ;  cDNASwapA (myHits) ; arraySort (myHits, cDNAOrderByA1) ;  
	  if (i1 < intronMax)
	    { *up = *xp ; *wp = *yp ; }
	  else
	    { 
	      /* reextend */
	      xp->a2 = wp->a2 ; xp->x2 = wp->x2 ;
	      countErrors(cosmid, myHits, dna, dnaR, 0, FALSE) ; /* implies cDNASwapX */
	      *wp = *xp ; up->est = 0 ;
	    }
	  /* a1 = up->a1 ; */  a2 = up->a2 ; b1 = wp->a1 ; /* b2 = wp->a2 ; */
	}
    }
  return foundJump ;
} /* teamFirstIntron */

/*********************************************************************/

static BOOL teamLastIntron (KEY cosmid, Array introns,  HIT *up, HIT *wp, Array dna, Array dnaR)
{
  HIT *vp, *xp, *yp ;
  int a2, b1, i1, nerr, newnerr, derr, di, da, intronMax = arrayMax(introns), nn, nnmax ;
  BOOL foundJump = FALSE, gt_ag = FALSE ;
  static Array myHits = 0 ;
  
  if (!myHits)
    myHits = arrayCreate (2, HIT) ;
  
  /* est = wp->est ; a1 = up->a1 ; */ a2 = up->a2 ; b1 = wp->a1 ; /* b2 = wp->a2 ; */
  nerr = up->nerr + wp->nerr ;

  for (i1 = nn = nnmax = 0 ; i1 < intronMax ; i1++) /* notice that i may try several solutions */
    {
      vp = arrp (introns, i1, HIT) ;
      if (vp->a1 == a2 && vp->a2 == b1) /* self */
        { 
          if (!vp->x1) /* not good intron */
	    nerr += 3 ;
	  else
	    gt_ag = TRUE ;
          break ;
        }
    }

  for (i1 = nn = nnmax = 0 ; i1 <= intronMax ; i1++) /* notice that i may try several solutions */
    {
      if (i1 < intronMax)
	{
	  vp = arrp (introns, i1, HIT) ;
	  if (! (vp->x1 && /* only regularize towards good introns */
		 vp->a1 <= a2 && vp->a1 > a2 - 8 && vp->a2 < b1)) /* closer jump */         
	    continue ;
	  di = b1 - vp->a2 ; /* intron shortening */
	  derr = di > 3000 ? 2 : (di > 1000 ? 1 : 0) ;
	  xp = arrayp (myHits, 0, HIT) ;  *xp = *up ;
	  yp = arrayp (myHits, 1, HIT) ;  *yp = *wp ;
	  yp->a1 = vp->a2 ; yp->a2 = vp->a2 + (wp->a2 - wp->a1) ;
	  da = a2 - vp->a1 ;
	  if (da)
	    {
	      xp->a2 -= da ; 
	      yp->a2 += da ;
	      if (up->x1 < up->x2) { xp->x2 -= da ; yp->x1 -= da ; }
	      else { xp->x2 += da ; yp->x1 += da ; }
	    }
	}
      else /* try locally */
	{
	  xp = arrayp (myHits, 0, HIT) ;  *xp = *wp ;
	  yp = arrayp (myHits, 1, HIT) ; 
	  xp->a1 = up->a2 + 1 ;
	  xp->a2 = up->a2 + wp->a2 - wp->a1 ;
	  nerr = wp->nerr + (gt_ag ? 0 : 3) ;
	  yp->nerr = 0 ; 
	  derr = gt_ag ? 0 : +3 ; /* i like to kill non gt_ag */
	  arrayMax (myHits) = 1 ;
	}
      countErrors(cosmid, myHits, dna, dnaR, 0, FALSE) ; /* implies cDNASwapX */
      newnerr = xp->nerr + yp->nerr ;
      if (newnerr <= nerr + derr && xp->maxExact > up->maxExact - 10)
	{ 
	  foundJump = TRUE ;  cDNASwapA (myHits) ; arraySort (myHits, cDNAOrderByA1) ; 
	  if (i1 < intronMax)
	    { *up = *xp ; *wp = *yp ; }
	  else
	    { 
	      /* reextend */
	      up->a2 = xp->a2 ; up->x2 = xp->x2 ;
	      *xp = *up ;
	      countErrors(cosmid, myHits, dna, dnaR, 0, FALSE) ; /* implies cDNASwapX */
	      cDNASwapA (myHits) ;
	      *up = *xp ; *wp = *up ; wp->est = 0 ; 
	    }
	  a2 = up->a2 ; b1 = wp->a1 ; 
	}
    }
  return foundJump ;
} /* teamLasttIntron */

/*********************************************************************/

static BOOL teamCentralJump (KEY cosmid, Array introns,  HIT *up, HIT *wp, Array dna, Array dnaR)
{
  int dSmall = 18 ;
  HIT *vp, *xp, *yp ;
  int a2, b1, i1, da, di, da1, nerr, newnerr, intronMax = arrayMax(introns), nn, nnmax, gt_ag ;
  BOOL foundJump = FALSE, reverse ;
  static Array myHits = 0 ;
  
  if (!myHits)
    myHits = arrayCreate (2, HIT) ;
  
  a2 = up->a2 ; b1 = wp->a1 ; 
  da = b1 - a2 ;
  nerr = up->nerr + wp->nerr ;
  reverse = up->x1 > up->x2 ? TRUE : FALSE ;
  gt_ag = 0 ; /* intron is gt_ag */

  for (i1 = nn = nnmax = 0 ; i1 < intronMax ; i1++) /* notice that i may try several solutions */
    {
      vp = arrp (introns, i1, HIT) ;
      if (vp->a1 == a2 && vp->a2 == b1) /* self */
        { 
          nn = vp->nerr ;  /* number of supporting clones */
          gt_ag = vp->x1 ; /* good intron */
          continue ;
        }
      if (vp->nerr > nnmax && 
          vp->x1 && /* only regularize towards good introns */
          ((vp->a1 <= a2 + dSmall && vp->a2 >= a2 - dSmall) || (vp->a1 <= b1 + dSmall && vp->a2 >= b1 - dSmall))
          )  
        nnmax = vp->nerr ; /* max number of clones sustaining a candidate */
    }

  if (!nnmax || nn > 3 * nnmax || (gt_ag && nn > 3 && 3*nn > nnmax && nerr < 3))
    return FALSE ;
  
  if (!gt_ag) nerr += 2 ;
  if (5*nn < nnmax && nerr > 3) nerr += 5 ; /* disadvantage rare introns */
  /* if same length +- 1 and shifted by up to 6bp
   * regularize no questions asked
   * a blue is worse than one more error
   * notice that this is much weaker than imposing gt-ag
  */
  {
    int bestda1 = 999 ;
    HIT *bestVp = 0 ;
    for (i1 = 0 ; i1 < intronMax ; i1++) /* notice that i try just one solution */
      {
        vp = arrp (introns, i1, HIT) ;
        if (vp->a1 == a2 && vp->a2 == b1) /* self */
          continue ;
        if (! vp->x1) /* only regularize towards good introns */
          continue ; 
        da1 = vp->a1 - a2 ;
        if (da1 < -5 || da1 > 5)
          continue ;
        di = vp->a2 - vp->a1 ;
        if (da - di < -1 || da - di > 1) /* size difference of the introns */
          continue ;
        if (da1 * da1 < bestda1)
          { 
            bestda1 = da1 * da1 ; bestVp = vp ;
            foundJump = TRUE ; 
          }
      }        
    if (foundJump)
      {
        up->a2 = bestVp->a1 ;
        wp->a1 = bestVp->a2 ;
        return TRUE ;
      }
  }

  if (nn < 3 || 5*nn < nnmax)  /* try to slide, but just to very close neighbours */
    {
      for (i1 = 0 ; i1 < intronMax ; i1++) /* notice that i may try several solutions */
        {
          vp = arrp (introns, i1, HIT) ;
          if (vp->a1 == a2 && vp->a2 == b1) /* self */
            continue ;
          if (! vp->x1) /* only regularize towards good introns */
            continue ; 
          di = vp->a2 - vp->a1 ;
          if (da - di < 2 && da - di > - 2 &&
              3 * vp->nerr >= nn  &&  /* not a minority intron */
              vp->a1 < a2 + dSmall && vp->a1 > a2 - dSmall &&
              vp->a2 < b1 + dSmall && vp->a2 > b1 - dSmall) /* intron is close enough to be a good candidate */
            {
              /* ddi = (da != di ? 1 : 0) ;  frame shift */
              xp = arrayp (myHits, 0, HIT) ;  *xp = *up ;
              yp = arrayp (myHits, 1, HIT) ;  *yp = *wp ;
              xp->a2 = vp->a1 ; if (reverse) xp->x2 -= vp->a1 - up->a2  ; else xp->x2 += vp->a1 - up->a2 ; 
              yp->a1 = vp->a2 ; if (reverse) yp->x1 -= vp->a2 - wp->a1  ; else yp->x1 += vp->a2 - wp->a1 ; 
              countErrors(cosmid, myHits, dna, dnaR, 0, FALSE) ; /* implies cDNASwapX */
              newnerr = xp->nerr + yp->nerr ;
              if (5*vp->nerr < nnmax && newnerr > 3) newnerr += 5 ; /* disadvantage rare introns */
              if (newnerr - di < nerr) /* to create a frameshift usually creates 2 errors */
                {   /* arp 14 2005, i changed from newnerr + 2 to newnerr -1 */
                  cDNASwapA (myHits) ; arraySort (myHits, cDNAOrderByA1) ; 
                  *up = *xp ; *wp = *yp ; nerr = xp->nerr + yp->nerr ; 
                  up->type |= gF ; wp->type |= gD ;
                  nerr = newnerr ;
                  nn = vp->nerr ;
                  foundJump = TRUE ; 
                  return TRUE ;
                }
            }
        }
    }
  
  if (! foundJump)
    for (i1 = 0 ; nerr > 0 && i1 < intronMax ; i1++) /* notice that i may try several solutions */
      {
        vp = arrp (introns, i1, HIT) ;
        if (vp->a1 == a2 && vp->a2 == b1) /* self */
          continue ;
        di = vp->a2 - vp->a1 ;
        if (da - di < 6 && da - di > - 6 &&
            3 * vp->nerr > nn  &&
            vp->a1 < a2 + dSmall && vp->a1 > a2 - dSmall &&
            vp->a2 < b1 + dSmall && vp->a2 > b1 - dSmall) /* intron is close enough to be a good candidate */
          {
            xp = arrayp (myHits, 0, HIT) ;  *xp = *up ;
            yp = arrayp (myHits, 1, HIT) ;  *yp = *wp ;
            xp->a2 = vp->a1 ; if (reverse) xp->x2 -= vp->a1 - up->a2  ; else xp->x2 += vp->a1 - up->a2 ; 
            yp->a1 = vp->a2 ; if (reverse) yp->x1 -= vp->a2 - wp->a1  ; else yp->x1 += vp->a2 - wp->a1 ; 
            countErrors(cosmid, myHits, dna, dnaR, 0, FALSE) ; /* implies cDNASwapX */
            newnerr = xp->nerr + yp->nerr ;
            if (!vp->x1) newnerr += 2 ;
            if (5*vp->nerr < nnmax && newnerr > 3) newnerr += 5 ; /* disadvantage rare introns */
            if (newnerr < nerr)
              { 
                cDNASwapA (myHits) ; arraySort (myHits, cDNAOrderByA1) ; 
                *up = *xp ; *wp = *yp ; nerr = xp->nerr + yp->nerr ; foundJump = TRUE ; 
                up->type |= gF ; wp->type |= gD ;
                nerr = newnerr ;
                nn = vp->nerr ;
              }
          }
      }

  if (! foundJump && 5*nn < nnmax && nerr >= 6)
    for (i1 = 0 ; nerr > 0 && i1 < intronMax ; i1++) /* notice that i may try several solutions */
      {
        vp = arrp (introns, i1, HIT) ;
        if (vp->a1 == a2 && vp->a2 == b1) /* self */
          continue ;
        if (! vp->x1) /* only regularize towards good introns */
          continue ; 
        if (vp->a1 > a2 + 25)
          break ;
        di = vp->a2 - vp->a1 ;
        if (da - di < 6 && da - di > - 6 &&  
            3 * vp->nerr > nn  &&
            vp->a1 < a2 + dSmall && vp->a1 > a2 - dSmall &&
            vp->a2 < b1 + dSmall && vp->a2 > b1 - dSmall) /* intron is close enough to be a good candidate */
          continue ; /* already tested in previous paragraph */
        if (da - di < 6 && da - di > - 6 &&
            vp->a1 < a2 + 25 && vp->a1 > a2 - 25 &&
            vp->a2 < b1 + 25 && vp->a2 > b1 - 25) 
          /* intron is close enough to be a good candidate */
          {
            xp = arrayp (myHits, 0, HIT) ;  *xp = *up ;
            yp = arrayp (myHits, 1, HIT) ;  *yp = *wp ;
            xp->a2 = vp->a1 ; if (reverse) xp->x2 -= vp->a1 - up->a2  ; else xp->x2 += vp->a1 - up->a2 ; 
            yp->a1 = vp->a2 ; if (reverse) yp->x1 -= vp->a2 - wp->a1  ; else yp->x1 += vp->a2 - wp->a1 ; 
            countErrors(cosmid, myHits, dna, dnaR, 0, FALSE) ; /* implies cDNASwapX */ 
            newnerr = xp->nerr + yp->nerr ;
            if (!vp->x1) newnerr += 2 ;
            if (5*vp->nerr < nnmax && newnerr > 3) newnerr += 5 ; /* disadvantage rare introns */
            if (newnerr < nerr - 3) /* real benefit */
              { 
                cDNASwapA (myHits) ; arraySort (myHits, cDNAOrderByA1) ; 
                *up = *xp ; *wp = *yp ; nerr = xp->nerr + yp->nerr ; foundJump = TRUE ; 
                up->type |= gF ; wp->type |= gD ; 
                nerr = newnerr ;
                nn = vp->nerr ;
              }
          }
      }
  return foundJump ;
} /* teamCentralJump */

/*********************************************************************/

static Array teamCreateIntrons (Array hits, Array dna, Array dnaR)
{
  int ii, jj, j, imax = hits ? arrayMax(hits) : 0 ;
  HIT *up, *vp ;
  Array introns = 0 ;
  char *cp, *cq ;
  
  /* register occurences of each intron */
  if (arrayMax(hits) < 2)
    return 0 ;
  introns = arrayCreate (imax, HIT) ;
  for (ii = 1, jj = 0, up = arrp(hits, ii, HIT); ii < imax ; up++, ii++)
    {      
      if (keyFindTag (up->cDNA_clone, _Suspected_internal_deletion))
        continue ;
      if (up->est == (up-1)->est)
        {
          for (j = 0 ; j < jj ; j++)
            {
              vp = arrp (introns, j, HIT) ;
              if (vp->a1 == (up-1)->a2 &&
                  vp->a2 == up->a1)
                { vp->nerr++ ; break ; } /* number of supporting clones */
            }
          if (j == jj)
            {
              vp = arrayp (introns, jj++, HIT) ;
              vp->a1 = (up-1)->a2 ;
              vp->a2 = up->a1 ;
              vp->reverse = (up->x1 < up->x2 ? up->reverse : !up->reverse) ;
              vp->nerr++ ; /* number of supporting clones */
            }
        }
    }
  
  arraySort (introns, cDNAOrderIntronsTeams) ;
  arrayCompress (introns) ;
  
  /* check if gt_ag */
  for (j = 0 ; j < arrayMax(introns) ; j++)
    {
      vp = arrp (introns, j, HIT) ;
      if (!vp->reverse)
        {
	  if (vp->a2 >= 3)
	    {
	      cp = arrayp (dna, vp->a1, char) ;
	      cq = arrayp (dna, vp->a2 - 3, char) ;
	      if (*cp == G_ && (*(cp+1) == T_ || *(cp+1) == C_) &&
		  *cq == A_ && *(cq+1) == G_)
		vp->x1 = 1 ;   /*   GOOD INTRON */
	    }
	}
      else
        {
	  if (vp->a2 >= 3)
	    {
	      cp = arrayp (dna, vp->a1, char) ;
	      cq = arrayp (dna, vp->a2 - 3, char) ;
	      if (*cp == C_ && *(cp+1) == T_ &&
		  (*cq == A_ || *cq == G_) && *(cq+1) == C_)
		vp->x1 = 1 ;  /*   GOOD INTRON */
	    }
        }
    }

#if 0
/* keep happy few */
  for (j = jj = 0, up = vp = arrp (introns, 0, HIT)  ; j < arrayMax(introns) ; j++, vp++)
    {
      if (vp->x1)
        {
          if (jj < j) *up = *vp ;
          jj++ ; up++ ;
        }
    }
  arrayMax (introns) = jj ;
#endif

  if (0) showHits(introns) ;
  return introns ;
}

/*********************************************************************/

static BOOL teamJump (KEY cosmid, Array hits, Array dna, Array dnaR)
{
  BOOL result = FALSE ;
  int ii, jj, imax = hits ? arrayMax(hits) : 0 ;
  HIT *up, *vp ;
  KEY est ;
  Array introns, oldHits ;
  
  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;
  
  introns = teamCreateIntrons (hits, dna, dnaR) ;

  if (!introns)
    return FALSE ;

  imax = arrayMax(hits) ;
  for (ii = 1 ; ii < imax ; ii++)
    {
      vp = arrayp (hits, ii, HIT) ;
      if (vp->est && vp->est == (vp - 1)->est &&
	  (ii == 1 || vp->est != (vp - 2)->est))  /* inside first intron */
        result |= teamFirstIntron (cosmid, introns, vp - 1, vp, dna, dnaR) ;
    }

  imax = arrayMax(hits) ;
  for (ii = 1 ; ii < imax ; ii++)
    {
      vp = arrayp (hits, ii, HIT) ;
      if (vp->est && vp->est == (vp - 1)->est &&
	  (ii == imax - 1 || vp->est != (vp + 1)->est))  /* inside last intron */
        result |= teamLastIntron (cosmid, introns, vp - 1, vp, dna, dnaR) ;
    }

  imax = arrayMax(hits) ;
  for (ii = 1 ; ii < imax ; ii++)
    {
      vp = arrayp (hits, ii, HIT) ;
      if (1 && vp->est && vp->est == (vp - 1)->est)  /* inside an intron */
        result |= teamCentralJump (cosmid, introns, vp - 1, vp, dna, dnaR) ;
    }

  if (result)
    {
      arrayDestroy (introns) ;
      introns = teamCreateIntrons (hits, dna, dnaR) ;
    }

  oldHits = arrayCopy (hits) ;
  hits = arrayReCreate (hits, imax, HIT) ;
  
  if (0) showHits(oldHits) ;
  for (ii = jj = 0 ; ii < imax ; ii++)
    {
      up = arrp(oldHits, ii, HIT) ; 
      est = up->est ;
      if (!est) continue ;
      vp = arrayp (hits, jj++, HIT) ;
      *vp = *up ;
      
      /*
      da = 0 ;
      if (ii < imax - 1 && up->est == (up+1)->est)
        da = (up+1)->a1 - up->a2 ;
      */
      if (
          (up->nerr >= 0 || (up->x1 < up->x2 && up->x1 > up->clipTop)) &&
          (ii == 0 || up->est != (up - 1)->est)  /* entering a new est */
          )
        result |= teamInitialJump (cosmid, introns, oldHits, hits, ii, &jj, dna, dnaR) ;
      if (
          (up->nerr >= 0 || (up->x1 > up->x2 && up->x2 > up->clipTop)) &&
          (ii == imax - 1 || up->est != (up+1)->est)  /* leaving an est */
          )
        switch (teamFinalJump (cosmid, introns, oldHits, hits, ii, &jj, dna, dnaR))
          {
          case 0: break ;
          case 1: result = TRUE ; break ;
          case 2: /* the last intron was removed we must go back in oldhits */
            up->est = 0 ; /* so previous seems to be last */
            if (ii >= 2 &&
                (up-1)->est == est)
              ii -= 2 ;
            result = TRUE ;
            break ;
          }
    }

  imax = arrayMax(hits) ;
  /* look for short exon between 2 introns */
  for (ii = 1 ; ii < imax - 1 ; ii++)
    {
      vp = arrayp (hits, ii, HIT) ;
      if (vp->a2 < vp->a1 + 15 &&
          vp->nerr &&
          vp->est &&
          vp->est == (vp - 1)->est &&
          vp->est == (vp + 1)->est)
        result |= teamDoubleJump (cosmid, introns, vp - 1, vp, vp+1, dna, dnaR) ;
    }

  arrayDestroy (introns) ;
  /* keep happy few */  
  for (ii = jj = 0, up = vp = arrp (hits, 0, HIT) ; ii < arrayMax (hits) ; up++, ii++)
    {
      if (! up->est)
        continue ;
      if (vp < up) 
        *vp = *up ;
      vp++; jj++ ; 
    }
  arrayMax (hits) = jj ;

  /* to see the printout */
  if (0)
    {
      introns = teamCreateIntrons (hits, dna, dnaR) ;
      arrayDestroy (introns) ;
    }
  arrayDestroy (oldHits) ;
  return result ;
} /* teamJump */

/*********************************************************************/
/*********************************************************************/

static void exportConfirmedGeneIntrons (KEY gene, BOOL only_alter, BOOL nonSliding, int nBase, int isShort)
{
  Array units = 0 ;
  int i, ii, jj, x1, x2, z1, z2, a1, a2 ;
  OBJ Cosmid = 0, Gene = 0 ;
  BSunit *up, *vp ;
  Array dna = 0, dnaR = 0 ;
  char *cp, *cq, buffer[60001] ;
  KEY cosmid  = keyGetKey (gene, _Genomic_sequence) ;
  int inContext = nonSliding ? 1 : 0 ;

  if (!cosmid || !(Cosmid = bsCreate(cosmid)) ||
      !bsFindKey (Cosmid, _Transcribed_gene, gene) ||
      !bsGetData (Cosmid, _bsRight, _Int, &a1) ||
      !bsGetData (Cosmid, _bsRight, _Int, &a2) ||
      !(Gene = bsUpdate (gene)))
    goto abort ;

  units = arrayCreate (20, BSunit) ;
  dna = dnaGet (cosmid) ;
  if (!dna) goto abort ;
  bsGetArray (Gene, _Splicing, units, 4) ;
  for (ii = 0 ; ii < arrayMax(units) ; ii+= 4)
    {
      up = arrp(units, ii, BSunit) ;
      if (up[2].k != _Intron &&
          up[2].k != _Alternative_intron)
        continue ;
      x1 = up[0].i ; x2 = up[1].i ;
      if (only_alter)
        {
          BOOL isAlter = FALSE ;
          for (jj = 0 ; jj < arrayMax(units) ; jj+= 3)
            { 
              vp = arrp(units, jj, BSunit) ;
              if (ii != jj && 
                  vp[0].i <= x2 &&  vp[1].i >= x1)
                isAlter = TRUE ;              
            }
          if (!isAlter)
            continue ;
        }
      switch (isShort)
        {
        case 0:
          if (x1 > x2 - 3)
            continue ;
          break ;
        case 1:  /* short exons */
          if (x1 > x2 - 3 || x1 <  x2 - 64)
            continue ;
          break ;
        case 2:  /* long exons */
          if (x1 > x2 - 65)
            continue ;
          break ;
        }
      if (a1 < a2) { x1 += a1 - 1 ; x2 += a1 - 1 ; }
      else { x1 = a1 - x1 + 1 ; x2 = a1 - x2 + 1 ; }
      if (x1 < 1 || x1 >= arrayMax(dna) ||
          x2 < 1 || x2 >= arrayMax(dna))
        continue ;
      z1 = x1 ; z2 = x2 ;
      if (!dnaR && z1 > z2)
        { dnaR = dnaCopy (dna) ; reverseComplement (dnaR) ; }
      
      if (z1 == z2 ||
          z1 < inContext + 1 ||
          z2 + inContext >= arrayMax(dna))
        continue ;
      cq = buffer ; i = 0 ;
      if (z1 <= z2)
        { 
          if (nonSliding &&
              ( arr (dna, z1 - 1, char) == arr (dna, z2, char) ||
                arr (dna, z1 - 2, char) == arr (dna, z2 - 1, char))
              )
            continue ;
          
          cp = arrp(dna, z1 - 1, char) ;
          while (z1 <= z2 && i < 60000 && z1 < arrayMax(dna) && z1 >= 0)
            { *cq++ = dnaDecodeChar[(int)(*cp++)] ; z1++ ; i++ ; }
          if (z1 != z2 + 1 || !i)
            continue ;
        }
      else
        { 
          z1 =  arrayMax(dna) - z1 ; z2 = arrayMax(dna) - z2 ; 
          if (nonSliding &&
              ( arr (dna, z1 - 1, char) == arr (dna, z2, char) ||
                arr (dna, z1 - 2, char) == arr (dna, z2 - 1, char))
              )
            continue ;
          
          cp = arrp(dnaR, z1, char) ;
          cq = buffer ; i = 0 ;
          while (z1 <= z2 && i < 60000 && z1 < arrayMax(dna) && z1 >= 0)
            { *cq++ = dnaDecodeChar[(int)(*cp++)] ; z1++ ; i++ ; }
          if (z1 != z2 + 1 || !i)
            continue ;
        }
      *cq = 0 ;
      freeOut (buffer) ; /* too big for messprintf! */
      freeOutf (" %s %d %d\n", name(cosmid), x1, x2) ;
    }
  bsSave (Gene) ;
abort:
  bsDestroy (Gene) ;
  bsDestroy (Cosmid) ;
  arrayDestroy (dna) ; arrayDestroy (dnaR) ; 
  arrayDestroy (units) ; 
}

/*********************************************************************/

static void exportIntronsTrCoordinates (KEY trans, int offset)
{
  Array units = 0 ;
  int ii, x1, x2,  a1, a2 ;
  OBJ Trans = 0 ;
  BSunit *up ;
  KEY map = 0 ;
  char *foot = 0 ;

  map = keyGetKey (trans, _IntMap) ;
  if (!map)
    return ;
  
  if (!trans || !map || !(Trans = bsCreate(trans)) ||
      !bsFindKey (Trans, _IntMap, map) ||
      !bsGetData (Trans, _bsRight, _Int, &a1) ||
      !bsGetData (Trans, _bsRight, _Int, &a2))
    goto abort ;

  units = arrayCreate (600, BSunit) ;

  bsGetArray (Trans, _Splicing, units, 4) ;
  for (ii = 0 ; ii < arrayMax(units) ; ii+= 4)
    {
      up = arrp(units, ii, BSunit) ;
      if (up[2].k != _Intron &&
          up[2].k != _Alternative_intron)
        continue ;
      x1 = up[0].i ; x2 = up[1].i ;
      foot = up[3].s ;

      if (a1 < a2) { x1 += a1 - 1 ; x2 += a1 - 1 ; }
      else { x1 = a1 - x1 + 1 ; x2 = a1 - x2 + 1 ; }
      
      x1 += 1 - offset ;
      x2 += 1 - offset ;
      if (foot &&
          (
           !strcmp("gt_ag", foot) ||
           !strcmp("gc_ag", foot))
          )
          freeOutf ("%s\t%s\t%s\t%d\t%d\n", name(map), name (trans), foot, x1, x2) ;
      if (foot &&
          (
           !strcmp("ct_ac", foot) ||  /* change strand x1  <--> x2 */
           !strcmp("ct_gc", foot))
          )
          freeOutf ("%s\t%s\t%s\t%d\t%d\n", name(map), name (trans), "gt_ag", x2, x1) ;
    }

abort:
  bsDestroy (Trans) ;
  arrayDestroy (units) ; 
}

static void exportIntronsCoordinates (KEY gene, int offset)
{
  KEYSET ks = queryKey (gene, ">Transcript") ;
  int ii, iTrans = keySetMax (ks), n, bestn = -1 , besti = -1 ;
  OBJ Trans = 0 ;
  Array units ;
  BSunit *up ;
  KEY trans ;
  char *foot ;

  units = arrayCreate (600, BSunit) ;

  while (iTrans--)
    {
      trans = keySet(ks, iTrans) ;
      Trans = bsCreate (trans) ;
      if (!Trans)
        continue ;
      n = 0 ;
      bsGetArray (Trans, _Splicing, units, 4) ;
      for (ii = 0 ; ii < arrayMax(units) ; ii+= 4)
        {
          up = arrp(units, ii, BSunit) ;
          if (up[2].k != _Intron &&
              up[2].k != _Alternative_intron)
            continue ;
          foot = up[3].s ;
          
          if (foot &&
              (
               !strcmp("gt_ag", foot) ||
               !strcmp("gc_ag", foot))
              )
            n++ ;
        }
      bsDestroy (Trans) ;   
      if (n > bestn)
        { bestn = n ; besti = iTrans ; }
    }

  arrayDestroy (units) ; 
  if (besti > -1)
    exportIntronsTrCoordinates (keySet(ks, besti), offset) ;
  keySetDestroy (ks) ;
}

/*********************************************************************/

static void exportConfirmedIntrons (KEYSET ks, BOOL only_alter, KEY only_coord, BOOL nonSliding, int nBase, int ishort)
{
  KEY *kp, source = 0 ;
  int ii = 0, offset = 0 ;
  OBJ Source = 0 ;
  
  if (!ks || !keySetMax(ks))
    messout("No Transcribed_gene in active keyset, sorry") ;
  else
    {
      if (class(only_coord) == _VSequence &&
          (source = keyGetKey (only_coord, _Source)))
        {
          if ((Source = bsCreate (source)))
            {
              if (bsFindKey (Source, _Subsequence, only_coord))
                bsGetData (Source, _bsRight, _Int, &offset) ;
              bsDestroy (Source) ;
            }
        }
      
      ii = keySetMax(ks) ;
      messout ("// Exporting the%s introns%s of %d genes",
               nonSliding ? " nonSliding" : "",
               only_coord ? " coords" : "",
               keySetMax(ks)) ;
      kp = arrp(ks, 0, KEY) -1  ;
      while (kp++, ii--) 
        if (only_coord)
          exportIntronsCoordinates (*kp, offset) ;
        else
          exportConfirmedGeneIntrons (*kp, only_alter, nonSliding, nBase, ishort) ;
    }
}

/*********************************************************************/

static void exportOne3p (KEY tr, int nnn)
{
  int n, ne, nl, i, j, nstop = -1, utrl, npa = 0 ;
  OBJ SS, Tr = bsCreate (tr) ;
  KEY seq = 0, subSeq = 0, newName = 0, tg = 0 ;
  KEYSET ks = 0 ;
  Array dna = 0 ;
  char *buf = 0 ;
  /*
    les 300 bp
    tag du stop
    nom du gene
    position genetique 
    new name
    corresponding gf
    cosmid
    taille transcript
    */
  if (!Tr)
    return ;
  buf = messalloc(nnn + 1) ;
  bsGetKey (Tr, _Spliced_sequence, &seq) ;
  if (seq) subSeq = keyGetKey (seq, _Subsequence) ;
  bsGetKey (Tr, _From_gene, &tg) ;
  newName = keyGetKey (tg, _NewName) ;
  bsGetData (Tr, _Length_3prime_UTR, _Int, &utrl) ;
  bsGetData (Tr, _PolyA_after_base, _Int, &npa) ;
  bsDestroy (Tr) ;

  ks = queryKey (tr, ">Assembled_from_cDNA_clone ; Library = \"E*\" ") ;
  ne = keySetMax (ks) ; keySetDestroy (ks) ;
  ks = queryKey (tr, ">Assembled_from_cDNA_clone ; !(Library = \"E*\") ") ;
  nl = keySetMax (ks) ; keySetDestroy (ks) ;

  if (subSeq)
    if ((SS = bsCreate (subSeq)))
      {
        Array aa = arrayCreate (64, BSunit) ;

        if (bsGetArray (SS, _Source_Exons, aa, 2))
          {
            i = arrayMax(aa) ;
            if (i > 0)
              nstop = array(aa, i - 1, BSunit).i ;
          }
        bsDestroy (SS) ;
        arrayDestroy (aa) ;
      }
      

  dna = dnaGet (seq) ;
  if (dna)
    {
      n = arrayMax(dna) ;
      freeOutf ("%s\t%s\tlib_E %d\tlib_M %d\tstop %d\tlength %d\t3putrlength %d\tpolya %d\n", 
                name(tr), name(newName), ne, nl, nstop, n, utrl, npa) ;
      i = nnn + 1 ; while (i--) buf[i] = 0 ;
      for (i = n - nnn, j = 0 ; i < 0 ; i++, j++) buf [j] = ' ' ;
      for (; i < n && j <= nnn ; i++, j++) buf[j] = dnaDecodeChar[(int)arr(dna, i, char)] ;
      freeOutf ("%s\n",buf) ;
      arrayDestroy (dna) ;
    }

  freeOut("\n") ;
  messfree (buf) ;
}

/*********************************************************************/

static void export3p (KEYSET ks, int nnn)
{
  KEY *kp ;
  int ii = 0 ;

  if (!keySetMax(ks))
    messout("No Transcribed_gene in active keyset, sorry") ;
  else
    {
      ii = keySetMax(ks) ;
      messout ("// Exporting the %d last bases of %d transcripts",
               nnn, keySetMax(ks)) ;
      kp = arrp(ks, 0, KEY) -1  ;
      while (kp++, ii--) 
        exportOne3p (*kp, nnn) ;
    }
}

/*********************************************************************/
/*********************************************************************/

static void exportOneOligoExons (KEY gene)
{
  KEY clone = keyGetKey (gene, _Longest_cDNA_clone) ;

  freeOutf("cDNA_clone %s\n",name(clone)) ;

}

/*********************************************************************/

static void exportOligoExons (void)
{
  KEY *kp ;
  int ii ;
  KEYSET ks = 0 ;
  
  ks =  query (0, "Find Transcribed_gene") ;

  if (!keySetMax(ks))
    messout("No Transcribed_gene in active keyset, sorry") ;
  else 
    for (ii = 0, kp = arrp(ks, 0, KEY) ; 
         ii <  keySetMax(ks) ; kp++, ii++)
      exportOneOligoExons (*kp) ;

  keySetDestroy (ks) ;
}

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

static void showOligo (unsigned int oligo)
{
  int ii = 16 ;
  unsigned int uu ;
  while (ii--)
    {
      uu = ( oligo >> (2 * ii)) & 0x3 ;
      switch (uu)
        {
        case 0: printf ("A") ; break ;
        case 3: printf ("T") ; break ;
        case 1: printf ("G") ; break ;
        case 2: printf ("C") ; break ;
        }
    }
  printf ("\n") ;
  if (ii > 16) { showOligo(oligo) ; showDna (0,0) ; } /* to please the compiler */
}

/*********************************************************************/

void showKeySet (KEYSET ks)
{
  int i ;
  KEY key ;

  for (i = 0 ; ks && i < keySetMax(ks) ; i++)
    {
      if (0) showKeySet (0) ; /* self reference trick to avoid compiler warning */
      key = keySet (ks, i) ;
      printf ("%03d %s:%s\n", i, className(key), name(key)) ;
    }
}

/*********************************************************************/

void showKeySetKey (KEYSET ks, char *cp)
{
  int i ;
  KEY key ;

  for (i = 0 ; ks && i < keySetMax(ks) ; i++)
    {
      if (0) showKeySet (0) ; /* self reference trick to avoid compiler warning */
      key = keySet (ks, i) ;
      if (cp && *cp && !strcasecmp (name(key), cp))
        printf ("%03d %s:%s\n", i, className(key), name(key)) ;
    }
}

/*********************************************************************/

void showCHits(Array hits, char *cName)
{
  int j, aMax = 0 ;
  HIT *vp ;
  KEY clone = 0 ;

  if (cName)
    {
      if (sscanf (cName, "%d", &j) == 1 && j > 100)
        aMax = j ;
      else
        lexword2key (cName, &clone, _VcDNA_clone) ;
    }
  printf("\n") ; /* because under debugger this will purge stdout buffer */
  if (arrayExists(hits)) 
    for (j = 0  ; j < arrayMax(hits) ; j++)
      {
        vp = arrp(hits, j, HIT) ;
        if (clone && clone != vp->cDNA_clone) continue ;
        if (aMax && vp->a1 > aMax && vp->a2 > aMax)
          continue ;
        printf("%2d:gene: %s e:%s cl:%s  rv=%d zone=%d  a1=%d a2=%d   x1=%d x2=%d  nr=%d nra=%d cTop=%d cEnd=%d", 
               j, vp->gene ? name(vp->gene) : "0", 
               vp->est ? name(vp->est) : "0", 
               vp->cDNA_clone ? name(vp->cDNA_clone) : "0", 
               vp->reverse, vp->zone, vp->a1, vp->a2, vp->x1, vp->x2, vp->nerr, vp->nerrAll, vp->clipTop, vp->clipEnd) ;
        printf (" ex=%d type=", vp->maxExact) ;
        if (gGene & vp->type) printf ("Gene ") ;
        if (gGap & vp->type) printf ("Gap ") ;
        if (gLink & vp->type) printf ("Link ") ;
        if (gGhost & vp->type) printf ("Ghost ") ;
        if (gFuseToGhost & vp->type) printf ("FuseToGhost ") ;
        if (gDroppedGene & vp->type) printf ("DroppedGene ") ;
        if (gSuspect & vp->type) printf ("Suspect ") ;
        if (gMicro & vp->type) printf ("Micro-") ;
        if (gX & vp->type) printf ("Exon ") ;
        if (gI & vp->type) printf ("Intron ") ;
        if (gJ & vp->type) printf ("Fuzzy ") ;
        if (gD & vp->type) printf ("Debut ") ;
        if (gF & vp->type) printf ("Fin ") ;
        if (gDF & vp->type) printf ("DebutFuzzy ") ;
        if (gFF & vp->type) printf ("FinFuzzy ") ;
        if (gReal3p & vp->type) printf ("r3p ") ;
        if (gReal5p & vp->type) printf ("r5p ") ;
        if (gS & vp->type) printf ("SL ") ;
        if (gS0 & vp->type) printf ("SL0 ") ;
        if (gA & vp->type) printf ("polyA ") ;
        if (gB & vp->type) printf ("Alter ") ;
        if (gStolen & vp->type) printf ("Stolen ") ;
        if (gPredicted & vp->type) printf ("Predicted ") ;
        if (gCompleteCDS & vp->type) printf ("CDS ") ;
        
        printf ("\n") ;
      }
}

/***********/

void showHits(Array hits)
{ showCHits (hits, 0) ; }

/*********************************************************************/

static void showSpl (Array allSpl)
{
  int i, j ;
  Array spl ;
  SPL *ss ;
  
  if (arrayExists(allSpl)) 
    for (i = 0  ; i < arrayMax(allSpl) ; i++)
      {
        spl = arr (allSpl, i, Array) ;
        if (arrayExists(spl))
          for (j = 0 ; j < arrayMax(spl) ; j++)
            {
              ss = arrp(spl, j, SPL) ;
              printf("gene: %18s %5d %5d    %5d %5d ",
                     name(ss->gene), 
                     ss->a1, ss->a2, ss->x1, ss->x2) ;
              if (gGene & ss->type) printf ("Gene ") ;
              if (gGap & ss->type) printf ("Gap ") ;
              if (gGhost & ss->type) printf ("Ghost ") ;
              if (gFuseToGhost & ss->type) printf ("FuseToGhost ") ;
              if (gDroppedGene & ss->type) printf ("DroppedGene ") ;
              if (gX & ss->type) printf ("Exon ") ;
              if (gI & ss->type) printf ("Intron ") ;
              if (gD & ss->type) printf ("Debut ") ;
              if (gF & ss->type) printf ("Fin ") ;
              if (gS & ss->type) printf ("SL ") ;
              if (gS0 & ss->type) printf ("SL0 ") ;
              if (gA & ss->type) printf ("polyA ") ;
              printf ("\n") ;
            }
        printf ("\n") ;
  
	showSpl (0) ; /* for compiler happiness */
      }
}

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

static KEY alignChain (KEY cosmid1, Array linkPos)
{
  KEY link = 0, chromMap = 0 ;
  KEY key, cosmid2 = 0 ;
  KEY s1 = 0, s2 = 0 ;
  OBJ Source = 0, Link = 0 ;
  int ii1, ii2, i1, i2, a1, a2, b1, b2, x1, x2 ;
  HIT *h1, *h2, *hh, *hh1, *hh2 ;
  float u1 = 0, u2 = 0 ;

  showSpl (0) ; /* for compiler happiness */
  cosmid2 = cosmid1 ;
  while ((key = keyGetKey (cosmid2, _Overlap_right)))
    cosmid2 = key ;

  /* construct the list of sources for cosmid1 */
  s1 = cosmid1 ;
  a1 = a2 = b1 = b2 = 0 ; hh1 = hh2 = 0 ;
  x1 = x2 = -1 ; i1 = 0 ;  arrayMax (linkPos) = 0 ;
  while (s1 && (s2 = keyGetKey (s1, _Source)))
    {
      if ((Source = bsCreate (s2)) &&
          bsFindKey (Source, _Subsequence, s1) &&
          bsGetData (Source, _bsRight, _Int, &a1) &&
          bsGetData (Source, _bsRight, _Int, &a2))
        {
          hh = arrayp (linkPos, i1++, HIT) ;
          hh->gene = s2 ; hh->est = cosmid1 ;
          if (x1 == -1 && x2 == -1)
            { hh->a1 = x1 = a1 ; hh->a2 = x2 = a2 ; }
          else if (a1 < a2)
            { hh->a1 = a1 + x1 - 1 ; hh->a2 = a1 + x2 - 1 ; x1 = hh->a1 ; x2 = hh->a2 ;}
          else
            { hh->a1 = a1 - x1 + 1 ; hh->a2 = a1 - x2 + 1 ; x1 = hh->a1 ; x2 = hh->a2 ;}
        }
      bsDestroy (Source) ;
      s1 = s2 ;
    }

  /* construct the list of sources for cosmid2 */
  s1 = cosmid2 ;
  x1 = x2 = -1 ; i2 = i1 ;
  while (s1 && (s2 = keyGetKey (s1, _Source)))
    {
      if ((Source = bsCreate (s2)) &&
          bsFindKey (Source, _Subsequence, s1) &&
          bsGetData (Source, _bsRight, _Int, &a1) &&
          bsGetData (Source, _bsRight, _Int, &a2))
        {
          hh = arrayp (linkPos, i2++, HIT) ;
          hh->gene = s2 ; hh->est = cosmid2 ;
          if (x1 == -1 && x2 == -1)
            { hh->a1 = x1 = a1 ; hh->a2 = x2 = a2 ; }
          else if (a1 < a2)
            { hh->a1 = a1 + x1 - 1 ; hh->a2 = a1 + x2 - 1 ; x1 = hh->a1 ; x2 = hh->a2 ;}
          else
            { hh->a1 = a1 - x1 + 1 ; hh->a2 = a1 - x2 + 1 ; x1 = hh->a1 ; x2 = hh->a2 ;}
        }
      bsDestroy (Source) ; 
      s1 = s2 ;
    }

  /* gene->source  est->cosmid a1 a2 coords of cosmid in source */
  ii2 = i2 ; ii1 = i1 ; hh1 = 0 ;
  while (i1-- > 0 && i2-- > ii1)
    {
      h1 = arrayp (linkPos, i1, HIT) ;
      h2 = arrayp (linkPos, i2, HIT) ;
      if (h1->gene == h2->gene)
        { hh1 = h1 ; hh2 = h2 ; }
      else 
        break ;
    }
  if (hh1) /* success, something in common */
    {
      int aa1, aa2, bb1, bb2, v1, v2, vv ;
      
      /* define the junction and its coordinates */
      lexaddkey (messprintf("Chain_%s_%s", name(cosmid1), name(cosmid2)), &link, _VSequence) ;

      /* now topmost coordinates */

      hh1 = arrayp (linkPos, ii1 - 1, HIT) ;
      hh2 = arrayp (linkPos, ii2 - 1, HIT) ;
      if (pickMatch (name(hh1->gene), "CHROMO*"))
        {
          lexaddkey (name(hh1->gene), &chromMap, _VMap) ;
          aa1 = hh1->a1 ;  aa2 = hh1->a2 ; 
          bb1 = hh2->a1 ;  bb2 = hh2->a2 ; 
          
          v1 = aa1 ; 
          if (aa2 < v1) v1 = aa2 ;
          if (bb1 < v1) v1 = bb1 ;
          if (bb2 < v1) v1 = bb2 ;
          v2 = aa1 ; 
          if (aa2 > v2) v2 = aa2 ;
          if (bb1 > v2) v2 = bb1 ;
          if (bb2 > v2) v2 = bb2 ;
          /* reorient the link from cosmid 1 towards cosmid 2 */
          if (v1 == bb1 || v1 == bb2)
            { vv = v1 ; v1 = v2 ; v2 = vv ;}
        }
      if ((Link = bsUpdate (link)))
        {
          if (chromMap)
            {
              bsAddTag (Link, _Is_chain) ;
              vv = v2 - v1 ; 
              bsAddData (Link, _Left, _Int, &vv) ;
              bsAddKey (Link, _IntMap, chromMap) ;
              bsAddData (Link, _bsRight, _Int, &v1) ;
              bsAddData (Link, _bsRight, _Int, &v2) ;
              bsAddKey (Link, _Map, chromMap) ;
              bsPushObj (Link) ;
              u1 = v1 ; u2 = v2 ; /* cast to float */
              bsAddData (Link, _Left, _Float, &u1) ;
              bsAddData (Link, _Right, _Float, &u2) ;
              bsGoto (Link, 0) ;
            }
          bsSave (Link) ;
        } 
    }

  return link ;
}

/*********************************************************************/

static KEY alignLink (KEY cosmid1, KEY cosmid2, Array linkPos)
{
  KEY link = 0, chromMap = 0 ;
  BOOL ok = FALSE ;
  KEY s1 = 0, s2 = 0, dnaKey ;
  OBJ Cosmid = 0, Source = 0, Link = 0 ;
  int ii1, ii2, i1, i2, a1, a2, b1, b2, x1, x2, dna1, dna2 ;
  int aa1, aa2, bb1, bb2, v1, v2, vv ;

  HIT *h1, *h2, *hh, *hh1, *hh2 ;
  float u1 = 0, u2 = 0 ;

  if (!cosmid1 || !cosmid2)
    return 0 ;
  /* construct the list of sources for cosmid1 */
  s1 = cosmid1 ;
  a1 = a2 = b1 = b2 = v1 = v2 = 0 ; hh1 = hh2 = 0 ;
  x1 = x2 = -1 ; i1 = 0 ;  arrayMax (linkPos) = 0 ;
  while (s1 && (s2 = keyGetKey (s1, _Source)))
    {
      if ((Source = bsCreate (s2)) &&
          bsFindKey (Source, _Subsequence, s1) &&
          bsGetData (Source, _bsRight, _Int, &a1) &&
          bsGetData (Source, _bsRight, _Int, &a2))
        {
          hh = arrayp (linkPos, i1++, HIT) ;
          hh->gene = s2 ; hh->est = cosmid1 ;
          if (x1 == -1 && x2 == -1)
            { hh->a1 = x1 = a1 ; hh->a2 = x2 = a2 ; }
          else if (a1 < a2)
            { hh->a1 = a1 + x1 - 1 ; hh->a2 = a1 + x2 - 1 ; x1 = hh->a1 ; x2 = hh->a2 ;}
          else
            { hh->a1 = a1 - x1 + 1 ; hh->a2 = a1 - x2 + 1 ; x1 = hh->a1 ; x2 = hh->a2 ;}
        }
      bsDestroy (Source) ;
      s1 = s2 ;
    }

  /* construct the list of sources for cosmid2 */
  s1 = cosmid2 ;
  x1 = x2 = -1 ; i2 = i1 ; ok = TRUE ;
  while (s1 && (s2 = keyGetKey (s1, _Source)))
    {
      ok = FALSE ;
      if ((Source = bsCreate (s2)) &&
          bsFindKey (Source, _Subsequence, s1) &&
          bsGetData (Source, _bsRight, _Int, &a1) &&
          bsGetData (Source, _bsRight, _Int, &a2))
        {
          hh = arrayp (linkPos, i2++, HIT) ;
          hh->gene = s2 ; hh->est = cosmid2 ;
          if (x1 == -1 && x2 == -1)
            { hh->a1 = x1 = a1 ; hh->a2 = x2 = a2 ; }
          else if (a1 < a2)
            { hh->a1 = a1 + x1 - 1 ; hh->a2 = a1 + x2 - 1 ; x1 = hh->a1 ; x2 = hh->a2 ;}
          else
            { hh->a1 = a1 - x1 + 1 ; hh->a2 = a1 - x2 + 1 ; x1 = hh->a1 ; x2 = hh->a2 ;}
        }
      bsDestroy (Source) ; 
      s1 = s2 ;
    }

  /* gene->source  est->cosmid a1 a2 coords of cosmid in source */
  ii2 = i2 ; ii1 = i1 ; hh1 = 0 ;
  while (i1-- > 0 && i2-- > ii1)
    {
      h1 = arrayp (linkPos, i1, HIT) ;
      h2 = arrayp (linkPos, i2, HIT) ;
      if (h1->gene == h2->gene)
        { hh1 = h1 ; hh2 = h2 ; }
      else 
        break ;
    }
  if (hh1) /* success, something in common */
    {
      Source = bsUpdate (hh1->gene) ;
      a1 = hh1->a1 ;  a2 = hh1->a2 ; 
      b1 = hh2->a1 ;  b2 = hh2->a2 ; 

      /* define the junction and its coordinates */
      lexaddkey (messprintf("%s_%s", name(cosmid1), name(cosmid2)), &link, _VSequence) ;

      x1 = a1 ; 
      if (a2 < x1) x1 = a2 ;
      if (b1 < x1) x1 = b1 ;
      if (b2 < x1) x1 = b2 ;
      x2 = a1 ; 
      if (a2 > x2) x2 = a2 ;
      if (b1 > x2) x2 = b1 ;
      if (b2 > x2) x2 = b2 ;
      /* reorient the link from cosmid 1 towards cosmid 2 */
      if (x1 == b1 || x1 == b2) 
        { int xx = x1 ; x1 = x2 ; x2 = xx ; }
      bsAddKey (Source,_Subsequence, link) ;
      bsAddData (Source, _bsRight, _Int, &x1) ;
      bsAddData (Source, _bsRight, _Int, &x2) ;
      ok = TRUE ;

      bsSave (Source) ;


      /* now topmost coordinates */

      hh1 = arrayp (linkPos, ii1 - 1, HIT) ;
      hh2 = arrayp (linkPos, ii2 - 1, HIT) ;
      if (TRUE || pickMatch (name(hh1->gene), "CHROMO*"))
        {
          lexaddkey (name(hh1->gene), &chromMap, _VMap) ;
          aa1 = hh1->a1 ;  aa2 = hh1->a2 ; 
          bb1 = hh2->a1 ;  bb2 = hh2->a2 ; 
          
          v1 = aa1 ; 
          if (aa2 < v1) v1 = aa2 ;
          if (bb1 < v1) v1 = bb1 ;
          if (bb2 < v1) v1 = bb2 ;
          v2 = aa1 ; 
          if (aa2 > v2) v2 = aa2 ;
          if (bb1 > v2) v2 = bb1 ;
          if (bb2 > v2) v2 = bb2 ;
          /* reorient the link from cosmid 1 towards cosmid 2 */
          if (v1 == bb1 || v1 == bb2)
            { vv = v1 ; v1 = v2 ; v2 = vv ;}
        }

      if ((Link = bsUpdate (link)))
        {
          bsAddKey (Link, str2tag("Parts"), cosmid1) ;
          bsAddKey (Link, str2tag("Parts"), cosmid2) ;
          if (chromMap)
            {
              bsAddKey (Link, _IntMap, chromMap) ;
              bsAddData (Link, _bsRight, _Int, &v1) ;
              bsAddData (Link, _bsRight, _Int, &v2) ;
              bsAddKey (Link, _Map, chromMap) ;
              bsPushObj (Link) ;
              u1 = v1 ; u2 = v2 ; /* cast to float */
              bsAddData (Link, _Left, _Float, &u1) ;
              bsAddData (Link, _Right, _Float, &u2) ;
              bsGoto (Link, 0) ;
            }
          bsSave (Link) ;
        } 
      hh = arrayp (linkPos, arrayMax(linkPos), HIT) ;
      hh->gene = link ;
      /* coord of cosmids 1 and 2 in the link coords */
      if (x1 < x2)
        {
          hh->a1 = a1 - x1 + 1 ;  hh->a2 = a2 - x1 + 1 ;
          hh->x1 = b1 - x1 + 1 ;  hh->x2 = b2 - x1 + 1 ;
        }
      else
        {
          hh->a1 = x1 - a1 + 1 ;  hh->a2 = x1 - a2 + 1 ;
          hh->x1 = x1 - b1 + 1 ;  hh->x2 = x1 - b2 + 1 ;
        }
      aa1 = hh->a1 ; aa2 = hh->a2 ; bb1 = hh->x1 ; bb2 = hh->x2 ;
      if (aa1 > aa2) { int tmp = aa1 ; aa1 = aa2 ; aa2 = tmp ; }
      if (bb1 > bb2) { int tmp = bb1 ; bb1 = bb2 ; bb2 = tmp ; }

      /* recompute real length using available dna */
      dna1 = dna2 = 0 ;
      if ((Cosmid = bsCreate (cosmid1)))
        {
          if (bsGetKey (Cosmid, _DNA, &dnaKey))
            bsGetData (Cosmid, _bsRight, _Int, &dna1) ;
          bsDestroy (Cosmid) ;
        }
      if ((Cosmid = bsCreate (cosmid2)))
        {
          if (bsGetKey (Cosmid, _DNA, &dnaKey))
            bsGetData (Cosmid, _bsRight, _Int, &dna2) ;
          bsDestroy (Cosmid) ;
        }
      aa2 = aa1 + dna1 - 1 ;
      bb1 = bb2 - dna2 + 1 ;
      if (aa2 < bb1 - 1)   /* there is a gap */
        {
          KEY gap = 0 ;
          OBJ Gap = 0 ;

          i1 = aa2 + 1 ; i2 = bb1 - 1 ;
          
          lexaddkey (messprintf("Gap_%s_%s", name(cosmid1), name(cosmid2)), &gap, _VSequence) ;
          if (gap && link && (Link = bsUpdate (link)))
            {
              bsAddKey (Link, _Subsequence, gap) ;
              bsAddData (Link, _bsRight, _Int, &i1) ;
              bsAddData (Link, _bsRight, _Int, &i2) ;
              bsSave (Link) ;
            } 
          if (v1 <= v2)
            { x1 = v1 + i1 - 1 ; x2 = v1 + i2 - 1 ; }
          else
            { x1 = v1 - i1 + 1 ; x2 = v1 - i2 + 1 ; }

          if ((Gap = bsUpdate(gap)))
            {
             bsAddTag (Gap, _Is_gap) ;
             b1 = i2 - i1 + 1 ;
             bsAddData (Gap, _bsRight, _Int, &b1) ;
             if (chromMap)
               {
                 bsAddKey (Gap, _IntMap, chromMap) ;
                 bsAddData (Gap, _bsRight, _Int, &x1) ;
                 bsAddData (Gap, _bsRight, _Int, &x2) ;
                 bsAddKey (Gap, _Map, chromMap) ;
                 bsPushObj (Gap) ;
                 u1 = x1 ; u2 = x2 ; /* cast to float */
                 bsAddData (Gap, _Left, _Float, &u1) ;
                 bsAddData (Gap, _Right, _Float, &u2) ;
                 bsGoto (Gap, 0) ;
               }
             bsSave (Gap) ;
            }
        }
    }

  return ok ? link : 0 ;
}

/*********************************************************************/

static Array alignGetOneCosmidWalls (KEY cosmid, Array walls, int a1, int a2)
{
  OBJ Cosmid = 0 ;
  Array units = 0 ;
  BSunit *up ;
  int ii, iWall, x1, x2, b1 = 0, b2 = 0 ;
  BOOL ok ;
  HIT *wp ;

  if (keyFindTag (cosmid, _Gene_wall))
    {
      if (!(Cosmid = bsCreate (cosmid)))
        return 0 ;

      units = arrayCreate (24, BSunit) ; 
      bsGetArray (Cosmid, _Gene_wall, units, 2) ;
      iWall = walls ? arrayMax(walls) : 0 ;
      for (ii = 0 ; ii < arrayMax(units) ; ii+= 2)
        {
          ok = FALSE ;
          up = arrp (units, ii, BSunit) ;
          x1 = up[0].i ;  x2 = up[1].i ;
          if (a1 < a2 && a1 <= x1 && a1 <= x2 && a2 >= x1 && a2 >= x2)
            {
              if (! getIgnoreWalls() || (x2 >= x1 + 10 || x1 >= x2 + 10))
                {
                  ok = TRUE ;
                  b1 = x1 - a1 + 1 ; b2 = x2 - a1 + 1 ;
                }
            }
          else if (a1 > a2 && a2 <= x1 && a2 <= x2 && a1 >= x1 && a1 >= x2)
            {
              if (! getIgnoreWalls() || (x2 >= x1 + 10 || x1 >= x2 + 10))
                {
                  ok = TRUE ;
                  b1 = a1 - x1 + 1 ; b2 = a1 -x2 + 1 ;
                }
            }
          if (ok)
            {
              if (!walls)
                walls = arrayCreate (12, HIT) ;
              up = arrp(units, ii, BSunit) ;
              wp = arrayp(walls, iWall++, HIT) ; 
              wp->a1 = b1 ; wp->a2 = b2 ;
              wp->reverse = b1 > b2 ? TRUE : FALSE ;
            }
        }
      
      bsDestroy (Cosmid) ;
    }
  arrayDestroy (units) ;
  return walls ;
}

/*********************************************************************/
static KEY alignGetparent (KEY cosmid, int *a1p, int *a2p)
{
  KEY parent = keyGetKey (cosmid, _Source) ;
  OBJ Parent = 0 ;

  if (parent &&
      (Parent = bsCreate (parent)))
    {
      if (!bsFindKey (Parent, _Subsequence, cosmid) ||
          !bsGetData (Parent, _bsRight, _Int, a1p) ||
          !bsGetData (Parent, _bsRight, _Int, a2p))
        parent = 0 ;
      bsDestroy (Parent) ;
    }
  return parent ;
}

static Array alignCosmidWalls (KEY cosmid)
{
  Array units = 0, walls = 0, walls1 = 0, walls2 = 0, links = 0 ;
  HIT *wp, *ww = 0, *lh ;
  int type = 0, ii, iWall ;
  KEY part1 = 0, part2 = 0 ;
  int a1, a2 ;

  if (!cosmid)
    return 0 ;

  if (keyFindTag (cosmid, str2tag("Parts")))
    type = 2 ;
  else
    type = 1 ;

  units = arrayCreate (24, BSunit) ; 
  switch (type)
    {
    case 1:
      
      a1 = 1 ; a2 = ACEDB_MAXINT ;
      while (cosmid)
        {
          walls =  alignGetOneCosmidWalls (cosmid, walls, a1, a2) ;
          cosmid = alignGetparent (cosmid, &a1, &a2) ;
        }
      break ;
    case 2:
      {
        KEYSET pp = queryKey (cosmid, ">Parts") ;
        if (keySetMax(pp) > 1)
          { /* order is not garanteed in acedb */
            part1 = keySet (pp, 0) ;
            part2 = keySet (pp, 1) ;
            if (part1 ==  keyGetKey (part2, _Overlap_right))
              { KEY p3 = part1 ; part1 = part2 ; part2 = p3 ; }
          }
        keySetDestroy (pp) ;
      }
      links = arrayCreate (4, HIT) ;
      alignLink (part1, part2, links) ;
      walls1 = alignCosmidWalls (part1) ;
      walls2 = alignCosmidWalls (part2) ;
      walls = arrayCreate (32, HIT) ;
      iWall = 0 ;
      if (arrayMax(links))
        {
          lh = arrayp (links, arrayMax(links) - 1, HIT) ;
          if (walls1)
            for (ii = 0 ; ii < arrayMax(walls1) ; ii++)
              {
                wp = arrayp(walls, iWall++, HIT) ;
                ww = arrp(walls1, ii, HIT) ;
                if (lh->a1 < lh->a2)
                  {
                    wp->a1 = lh->a1 + ww->a1 - 1 ;
                    wp->a2 = lh->a1 + ww->a2 - 1 ;
                    wp->reverse = ww->reverse ;
                  }
                else
                  {
                    wp->a1 = lh->a1 - ww->a1 + 1 ;
                    wp->a2 = lh->a1 - ww->a2 + 1 ;
                    wp->reverse = !ww->reverse ;
                  }
              }
          if (walls2)
            for (ii = 0 ; ii < arrayMax(walls2) ; ii++)
              {
                wp = arrayp(walls, iWall++, HIT) ;
                ww = arrp(walls2, ii, HIT) ;
                if (lh->x1 < lh->x2)
                  {
                    wp->a1 = lh->x1 + ww->a1 - 1 ;
                    wp->a2 = lh->x1 + ww->a2 - 1 ;
                    wp->reverse = ww->reverse ;
                  }
                else
                  {
                    wp->a1 = lh->x1 - ww->a1 + 1 ;
                    wp->a2 = lh->x1 - ww->a2 + 1 ;
                    wp->reverse = !ww->reverse ;
                  }
              }
        }
        
      arrayDestroy (walls1) ;
      arrayDestroy (walls2) ;
      arrayDestroy (links) ;
      break ;
    default:
      break ;
    }

  if (walls && arrayMax(walls))
    {
      cDNASwapA (walls) ;
      arraySort (walls, cDNAOrderByA1) ;
    }
  else
    arrayDestroy (walls) ;

  arrayDestroy (units) ;
  return walls ;
}

/***********/

static Array alignLinkWalls (KEY link, KEY cosmid, KEY cosmid2, Array linkPos)
{
  Array walls = 0, brick = 0 ;
  HIT *bp, *wp, *hh ;
  int ib, iw ;

  if (!cosmid) /* then create it */
    {
      KEYSET pp = queryKey (link, ">Parts") ;
      if (keySetMax(pp) > 1)
        { /* order is not garanteed in acedb */
          cosmid = keySet (pp, 0) ;
          cosmid2 = keySet (pp, 1) ;
          if (cosmid ==  keyGetKey (cosmid2, _Overlap_right))
            { KEY p3 = cosmid2 ; cosmid2 = cosmid ; cosmid = p3 ; }
        }
      keySetDestroy (pp) ;
      alignLink (cosmid, cosmid2, linkPos) ;
    }
        
  walls = alignCosmidWalls (cosmid) ;
  brick = alignCosmidWalls (cosmid2) ;
  if (!brick)
    return walls ;  /* may be null */

  if (!walls)
    walls = arrayCreate (12, HIT) ;
  
  hh = arrp(linkPos, arrayMax(linkPos) -1, HIT) ;
  iw = arrayMax (walls) ;
  for (ib = 0 ; ib < arrayMax(brick) ; ib++)
    {
      bp = arrp (brick, ib, HIT) ;
      wp = arrayp (walls, iw++, HIT) ;
      *wp = *bp ;
      if (hh->a1 < hh->a2)
        {
          wp->a1 = hh->x1 + bp->a1 - 1 ; 
          wp->a2 = hh->x1 + bp->a2 - 1 ;
        }
      else
        {
          wp->a1 = hh->x2 - bp->a1 - 1 ; 
          wp->a2 = hh->x2 - bp->a2 - 1 ;
        }
    }
  arrayDestroy (brick) ;

  if (!arrayMax(walls))
    arrayDestroy (walls) ;

  return walls ;
}

/*********************************************************************/
 /* direction 2 forward, 1 reverse, 3 forward puis reverse as -1 */
static KEYSET alignEst2Cosmid (KEY cosmid, KEYSET estSet, KEYSET currentGenes,
                               int isAlign, int type, int z1, int z2, char *nom, int direction,
                               int searchRepeats)
{
  static int pass = 0 ;
  int n1 = 0, n2 = 0, nDoRepeats = 0, nDoRepeats0 = 0, rNnewgenes = 0 ;
  KEY cosmidMap = 0 ; int cosmidC1 = 0, cosmidC2 = 0 ; /* IntMap of the cosmid */
  Array 
    dna = 0, dnaR = 0, linkPos = 0, walls = 0, 
    hits = 0, geneHits = 0, plainGeneHits = 0, 
    allSpl = 0, genes = 0, sets = 0 ;
  KEYSET newgenes = 0 ;
  static int maxIntronSize = -1 ;
  static Array words = 0 ;
  static Associator ass = 0 , assR = 0, assDifficult = 0, assDifficultR = 0 ;
  static KEYSET cDNA_clones = 0, clipTops = 0, clipEnds = 0, trackingClones = 0 ;
  static KEYSET estMaps = 0, estMap1 = 0, estMap2 = 0 ;
  KEY cosmid1 = 0, cosmid2 = 0, link = 0, previousCosmid = 0 , previousLink = 0 ;
  static int nDone = 0 ;
  Array rHitsAll = 0, usedZones = 0, usedQualities = 0 ;

  pass++ ;
  if (type == 9999) goto abort ; /* cleanup */
  if (maxIntronSize == -1)
    maxIntronSize = getMaxIntronSize () ;


  chrono ("alignEst2Cosmid") ;

  linkPos = arrayCreate (24, HIT) ;
  switch (type)
    {
    case 2: case 2001: case 2222:
      nDoRepeats = nDoRepeats0 = searchRepeats ? 20 : 0 ; /* 0: suppress , 20:  do ; */
      if (type != 2222)
        cdnaCleanUp (cosmid, 0, 0) ;
      link = cosmid2 = 0 ;
      cosmid2 = keyGetKey (cosmid, _Overlap_right) ;
      if (cosmid2)
        {
          link = alignLink (cosmid, cosmid2, linkPos) ;
          walls = alignLinkWalls (link, cosmid, cosmid2, linkPos) ;
        }
      else
        walls = alignCosmidWalls (cosmid) ;
      if ((previousCosmid = keyGetKey (cosmid, _Overlap_left)))
        lexword2key (messprintf ("%s_%s", name(previousCosmid), name(cosmid)),
                     &previousLink, _VSequence) ;                   
      if (link)
        {
          cosmid1 = cosmid ;
          cosmid = link ;
           if (type != 2222)
             {
               cdnaCleanUp (cosmid2, 0, 0) ;
               cdnaCleanUp (link, 0, 0) ;
             }
        }
      if (type == 2222) type = 2 ;
      break ;
    case 3:
      nDoRepeats = nDoRepeats0 = searchRepeats ? 20 : 0 ; /* 0: suppress , 20:  do ; */
      if (keyFindTag (cosmid, _Genomic))
        walls = alignCosmidWalls (cosmid) ;
      else if (keyFindTag (cosmid, str2tag("Parts"))) 
        walls = alignLinkWalls (cosmid, 0, 0, linkPos) ;
      break ;
    default:
      break ;
    }
  if (0)
    { 
      arrayDestroy(genes) ; 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }
  if (1)   /* get the cosmid global coordinates */
    {
      KEY key = 0, target = link ? link : cosmid ;
      OBJ Target = bsCreate (target) ;
      if (Target)
	{
	  if (bsGetKey (Target, _IntMap, &key) &&
	      bsGetData (Target, _bsRight, _Int, &cosmidC1) &&
	      bsGetData (Target, _bsRight, _Int, &cosmidC2)
	      )
	   cosmidMap = key ; /* IntMap of the cosmid */
	  bsDestroy (Target) ;
	}
    }
  dna = getSeqDna(link ? link : cosmid) ;
  if (!dna || arrayMax(dna) < 50) { direction = 0 ; goto abort ; }
  dnaR = dnaCopy (dna) ;
  reverseComplement (dnaR) ;
  if (0)
    { 
      arrayDestroy(genes) ; 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }

  if (type == 3 || !clipEnds) 
    {
      n1 = n2 = 0 ;
      cDNA_clones = keySetReCreate (cDNA_clones) ;
      clipTops = keySetReCreate(clipTops) ;
      clipEnds = keySetReCreate (clipEnds) ; 
      trackingClones = keySetReCreate (trackingClones) ; 
      estMaps = keySetReCreate(estMaps) ;
      estMap1 = keySetReCreate(estMap1) ;
      estMap2 = keySetReCreate(estMap2) ;
      assDestroy (ass) ; assDestroy (assR) ;

      if (!words) /*   was  && (type != 3) ) mieg jul 9 2002 */
        {
          char *wName = getMaskFrequentWords() ;

          if (wName)   /* this is really quite costly, do only for human genome */ 
            {
              unsigned long wnb = 0, wnw = 0 ;
              if (strcasecmp (wName, "Self"))
                words = dnaGetWordUsage (0, 1, &wnb, &wnw, FALSE, wName, 0) ;
              else /* extract the table from the dna of the section */
                {
                  KEYSET ksGenomic = query (0, "Find sequence genomic ; > source") ;
                  words = dnaGetWordUsage (ksGenomic, 30, &wnb, &wnw, TRUE, 0, 0) ;
                  keySetDestroy (ksGenomic) ;
                }
              if (0) messout ("Words: indexed %lu bases, found %lu different 15-mers\n", wnb, wnw) ;
            }
        }
      n2 = 0 ;
      ass = estAssCreate (words, estSet, FALSE
			  , cDNA_clones, clipTops, clipEnds, trackingClones
			  , estMaps, estMap1, estMap2
			  , &n1, &n2, type, OLIGO_LENGTH, FALSE
			  , cosmidMap, cosmidC1, cosmidC2
			  ) ;
      assR = estAssCreate (words, estSet, TRUE
			   , cDNA_clones, clipTops, clipEnds, trackingClones
			   , estMaps, estMap1, estMap2
			   , &n1, &n2, type, OLIGO_LENGTH, FALSE
			   , cosmidMap, cosmidC1, cosmidC2
			   ) ;
      assDifficult = estAssCreate (words, estSet, FALSE
				   , cDNA_clones, clipTops, clipEnds, trackingClones
				   , estMaps, estMap1, estMap2
				   , &n1, &n2, type, OLIGO_LENGTH_DIFFICULT, TRUE
				   , cosmidMap, cosmidC1, cosmidC2
				   ) ;
      assDifficultR = estAssCreate (words, estSet, TRUE
				    , cDNA_clones, clipTops, clipEnds, trackingClones
				    , estMaps, estMap1, estMap2
				    , &n1, &n2, type, OLIGO_LENGTH_DIFFICULT, TRUE
				    , cosmidMap, cosmidC1, cosmidC2
				    ) ;
      if (type != 3)
        messout ("Created an associator for %d cDNA, %d oligos, %s\n", n1, n2, timeShowNow()) ;
    }
  if (!ass && !assR && ! assDifficult && ! assDifficultR) { direction = 0 ; goto abort ; }
  if (0)
    { 
      arrayDestroy(genes) ; 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }

 ladirection:
  chrono ("getHitsLaDirection") ;
  hits  = arrayReCreate (hits, 4096,HIT) ;
  /* printf ("\n#%3d:: %s\t%s\t%s\n", direction, name(cosmid),  name(cosmid1),  name(cosmid2)) ; */
  if (0)
    { 
      arrayDestroy(genes) ; 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }

  if (direction * direction != 1)
    {
      getCosmidHits (cosmid, hits, ass, 0, dna, FALSE, clipEnds, OLIGO_LENGTH, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
      dropOutOfZonePreMappedHits (hits, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
      getCosmidHits (cosmid, hits, 0, assR, dna, FALSE, clipEnds, OLIGO_LENGTH, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
      dropOutOfZonePreMappedHits (hits, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
      getCosmidHits (cosmid, hits, assDifficult, 0, dna, FALSE, clipEnds, OLIGO_LENGTH_DIFFICULT, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
      getCosmidHits (cosmid, hits, 0, assDifficultR, dna, FALSE, clipEnds, OLIGO_LENGTH_DIFFICULT, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
      dropOutOfZonePreMappedHits (hits, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
    }
  else
    {
      getCosmidHits (cosmid, hits, 0, assR, dnaR, TRUE, clipEnds, OLIGO_LENGTH, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
      dropOutOfZonePreMappedHits (hits, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
      getCosmidHits (cosmid, hits, ass, 0, dnaR, TRUE, clipEnds, OLIGO_LENGTH, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
      dropOutOfZonePreMappedHits (hits, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
      getCosmidHits (cosmid, hits, 0, assDifficultR, dnaR, TRUE, clipEnds, OLIGO_LENGTH_DIFFICULT, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
      getCosmidHits (cosmid, hits, assDifficult, 0, dnaR, TRUE, clipEnds, OLIGO_LENGTH_DIFFICULT, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
      dropOutOfZonePreMappedHits (hits, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
    }

  if (0)
    { 
      arrayDestroy(genes) ; 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }
  cDNAFilterDirectionalHits (direction, hits, FALSE) ;
  dropOutOfZonePreMappedHits (hits, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
  if (!arrayMax(hits))
    goto abort ;

  arrayCompress (hits) ;
  getcDNA_clonesAndClips (hits, cDNA_clones, clipTops, clipEnds) ;
  cDNAFilterPreviousLink (direction, previousCosmid, previousLink, hits) ;

  if (0)
    { 
      arrayDestroy(genes) ; 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }


  if (!arrayMax(hits))
    goto abort ;

  if (type == 2 && !isAlign)
    cDNAFilterBestAlignments (cosmid, hits, 2) ;
  cDNAFilterDiscardedHits (cosmid, hits) ;
  cDNAFilterDiscardedHits (cosmid1, hits) ;
  cDNAFilterDiscardedHits (cosmid2, hits) ;
  cDNAFilterDirectionalHits (direction, hits, FALSE) ;

  if (z1 || z2) dropOutOfZoneHits (hits, z1, z2) ;
  cDNAFilterMultipleHits (cosmid, hits) ;

  if (!arrayMax(hits))
    goto abort ;

  extendMultipleHits (cosmid, hits, dna, dnaR) ;
  /* dropBackwardHits (dna, dnaR, hits) ; this code used to choose the best strand */
  mergeOverlappingHits  (hits) ;
  cDNAFilterDirectionalHits (direction, hits, FALSE) ;
  dropBadHits (dna, dnaR, hits, 6) ;  /* drop low entropy hits, s < 30 bp */

  if (0)
    { 
      arrayDestroy(genes) ; 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }

  dropOvernumerousHits(hits, 50000) ;
  countErrors (cosmid, hits, dna, dnaR, 0, FALSE) ;

  rHitsAll = arrayCopy (hits) ; /* all hits */
  rNnewgenes = currentGenes ? arrayMax (currentGenes) : 0 ;
  if (!arrayMax(hits))
    goto abort ;
  chronoReturn() ;
 doRepeats:

  chrono ("cPathSelectBestPath") ;
  cPathSelectBestPath (hits, FALSE, searchRepeats) ;  /* drop noncolinear hits */
  chronoReturn() ; chrono ("getHitsDoRepeats") ; 

  fixIntrons (cosmid, hits, dna, dnaR, 0) ;
  countErrors (cosmid, hits, dna, dnaR, 0, FALSE) ;
  dropHugeIntrons (hits, dna, dnaR, maxIntronSize, FALSE) ;  /* introns larger than the limit */

  getOtherExons  (cosmid, hits, dna, dnaR, maxIntronSize, usedZones, z1, z2, estMap1, estMap2, cosmidC1, cosmidC2) ;
  cDNAFilterDirectionalHits (direction, hits, FALSE) ; /* useful in some rare cases */
  cDNAFilterUsedZones (hits, usedZones) ;
  if (z1 || z2) dropOutOfZoneHits (hits, z1, z2) ;
  dropOutOfZonePreMappedHits (hits, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;

  mergeOverlappingHits  (hits) ; /* the new may overlap or mixup with the old in rare cases, so we must redo selectPath */ 
  countErrors (cosmid, hits, dna, dnaR, 0, FALSE) ; 

  if (0)
    { 
      arrayDestroy(genes) ; 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }

  chronoReturn() ;
  chrono ("cPathSelectBestPath") ;

  cPathSelectBestPath (hits, TRUE, searchRepeats) ;  /* drop noncolinear hits IMAGE:30524261 CBLAEG01 */

  if (!arrayMax(hits))
    goto abort ;
  dropBackToBackPairs (hits) ;  /* 3' and 5' turning away from each other */
  chronoReturn() ;

  chrono ("fixIntrons") ;
  fixIntrons (cosmid, hits, dna, dnaR, 1) ;

  cDNAFilterDirectionalHits (direction, hits, FALSE) ;  /* in case sliding went crazy */
  slideIntrons (cosmid, hits, dna, dnaR) ;  /* slide towards gt_ag gc_ag ... */
  countErrors (cosmid, hits, dna, dnaR, 0, TRUE) ; 
  dropBadHits (dna, dnaR, hits, 8) ;  /* drop composite reads that have over 1% error rate */
  dropBadHits (dna, dnaR, hits, 6) ;   /* drop low entropy hits, added jan 25 2002 fro human, mieg */
  dropBadHits (dna, dnaR, hits, 5) ;  /* drop small coverage of reads and very small or reverse hits, may be redundant on below */
  dropBadHits (dna, dnaR, hits, 7) ;  /* drop too short or quality coverage of reads and  quality 11 hits*/
  dropHugeIntrons (hits, dna, dnaR, 0, TRUE) ;  /* drop introns large  and large and not gt_ag */

  chronoReturn() ;

  if (0)
    { 
      arrayDestroy(genes) ; 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }

  chrono ("dropBadHits") ;

  countErrors (cosmid, hits, dna, dnaR, 0, TRUE) ;
  dropBadHits (dna, dnaR, hits, 2) ;  /* drop silly terminal exons */
  dropDoubleReads  (hits, clipEnds, dna, dnaR) ; /* clip double reading of 3p 5p */
  /* cdnaReverseDubiousClone (cosmid, hits,dna,dnaR) ;  */

  chronoReturn() ;
  chrono ("IntronWobbling") ;

  removeIntronWobbling (cosmid, hits, dna, dnaR) ; /* create an error in est rather than non gt_ag */

  /* one last clean up to protect mrnaCreate from any bug in this part of the code */
  cDNAFilterDirectionalHits (direction, hits, FALSE) ;
  dropBadHits (dna, dnaR, hits, 6) ;   /* drop again low entropy hits, added 2008_09_20 worm 5I318 */
  countErrors (cosmid, hits, dna, dnaR, 0, FALSE) ;

  /* one last turn to entraine */
  /* watch out, team jump may create one bp exon 
   * so from then on x1 < x2 becomes unreliable to judge orientation
   */
  if (1 && teamJump (cosmid, hits, dna, dnaR)) /* try to reuse the same intron at the bad edges of est */
    {
      cDNAFilterDirectionalHits (direction, hits, FALSE) ;
      dropOutOfZonePreMappedHits (hits, estMaps, estMap1, estMap2, cosmidMap, cosmidC1, cosmidC2, direction) ;
    }
  if (!arrayMax(hits))
    goto abort ;

  if (0)
    { 
      arrayDestroy(genes) ; 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }

  { /* force error recount */
    int ii = arrayMax(hits) ;
    HIT *up ;	
    while (ii--)
      { up = arrp (hits, ii, HIT) ; ; up->ex1 = -2 ; }
  }
  countErrors (cosmid, hits, dna, dnaR, 0, FALSE) ;
  chronoReturn() ;
  chrono ("constructGenes") ;

  if (0) /* here we DO NOT have the bug */
    { 
      arrayDestroy(genes) ; 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }

  if (isAlign == 0) /* cdna case */
    {
      KEYSET clonesWithBadIntron ;

      checkAllIntrons (cosmid, hits, dna, dnaR) ;
      coalignAllIntrons (cosmid, hits) ;
      clonesWithBadIntron = getClonesWithBadIntron (cosmid, hits, dna, dnaR) ;
      geneHits = getGeneHits (cosmid, hits, clipTops, clipEnds, walls, trackingClones, dna, dnaR) ;
      addGhostHits (arrayMax(dna), geneHits, walls) ;
      countErrors (cosmid, geneHits, dna, dnaR, 0, FALSE) ;
      if (0)  /* PAS de bug */
	{ 
	  arrayDestroy(genes) ; 
	  arrayDestroy (geneHits) ;
	  invokeDebugger () ;
	}
      if (nDoRepeats)
        {
          if (! usedQualities)
            usedQualities = arrayCreate (1024, HIT) ;
          cDNAFilterUsedQualities (geneHits, usedQualities) ;
        }
      plainGeneHits = arrayCopy (geneHits) ; 
      if (0)  /* on A pas le bug */
	{ 
	  arrayDestroy(genes) ; 
	  arrayDestroy (geneHits) ;
	  invokeDebugger () ;
	}
      cleanGeneHits (cosmid, geneHits, trackingClones, FALSE) ;
      if (0)  /* ??? on ?? le bug */
	{ 
	  arrayDestroy(genes) ; 
	  arrayDestroy (geneHits) ;
	  invokeDebugger () ;
	}
      cleanGeneHits (cosmid, geneHits, trackingClones, TRUE) ;
      if (0)  /* on A le bug */
	{ 
	  arrayDestroy (geneHits) ;
	  invokeDebugger () ;
	}

      {
        int j2 ;
        HIT *gh2 ;
        
        for (j2 = 0 ; j2 < arrayMax (plainGeneHits) ; j2++)
          { 
            gh2 = arrayp (plainGeneHits, j2, HIT) ;
            if (gh2->zone == 9999)
              gh2->zone = 0 ;
          }
      }

      if (arrayMax(geneHits))
        {
          flagTrackingCloneHits (plainGeneHits, trackingClones) ;
          /*      flagReverseReads (plainGeneHits) ; */
          if (1) 
            { 
              countErrors (cosmid, plainGeneHits, dna, dnaR, 600, FALSE) ; /* implies cDNASwapX */
              cDNASwapA (plainGeneHits) ;
              arraySort (plainGeneHits, cDNAOrderByA1) ; /* needed by mrnaCreate */
              newgenes = mrnaCreate (type, cosmid, cosmid1, dna, dnaR, plainGeneHits, geneHits, linkPos, clipTops, clipEnds, clonesWithBadIntron) ; 
              if (0) /* 2007_01_21: very costly (30000 hits in some tiles, and useless */ 
                saveHits (cosmid, 0, hits, 0, direction) ;
            }
        }
      keySetDestroy (clonesWithBadIntron) ;
    }
  else if (isAlign == 1)    /* align to reference case  g_yk5503_53 T10B10_H03A11 */
    {
      saveAlignment (cosmid, hits, dna, nom) ; 
    }
  else
    {
      countErrors (cosmid, hits, dna, dnaR, 600, FALSE) ; /* in view of global clean up */
      saveEst2CosmidQuality (cosmid, hits, linkPos, dna, dnaR) ;
    }
  chronoReturn() ;

abort:
  
  if (!++nDone % 100)
    tStatus () ;

  /* destroy what is strand specific */

  if (0)  /* here we have the bug */
    { 
      arrayDestroy(genes) ; 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }

  if (allSpl)
    {
      int i ;
      for (i = 0 ; i < arrayMax(allSpl) ; i++)
        {
          Array spl = array(allSpl, i, Array) ;
          arrayDestroy (spl) ;
        }
      arrayDestroy (allSpl) ;
    }

  if (sets)
    {
      int j1 ;
      for (j1 = 0 ; j1 < arrayMax (sets) ; j1++)
        {
          KEYSET ks = array (sets, j1, KEYSET) ;
          if (ks) keySetDestroy (ks) ;
        }
      arrayDestroy (sets) ;
    }

  if (0)
    { 
      arrayDestroy(genes) ; 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }

  if (newgenes) 
    {
      int i, j ;
      
      if (! keySetExists(currentGenes))
        currentGenes = keySetCreate () ;
      for (i = 0, j = keySetMax (currentGenes) ; i < keySetMax (newgenes) ; i++)
        {
          /* printf (" %s ", name( keySet (newgenes, i))) ; */
          keySet (currentGenes, j++) = keySet (newgenes, i) ;
        }
      keySetSort (currentGenes) ;
      keySetCompress (currentGenes) ;
    }

    if (0)
    { 
      arrayDestroy(genes) ; 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }

  if ((usedZones || nDoRepeats > 0) && plainGeneHits && arrayMax (plainGeneHits))
    {
      int i, j, ok =  1 ;
      HIT *up, *vp ;
      
      cDNASwapA (plainGeneHits) ;
      arraySort (plainGeneHits,  cDNAOrderByA1 ) ;


      arraySort (plainGeneHits,  cDNAOrderGloballyByA1 ) ;
      if (!usedZones)
        usedZones = arrayCreate (256, HIT) ;
      for (i = 0, up = arrp (plainGeneHits, 0, HIT) ; i < arrayMax (plainGeneHits) ; i++, up++)
        {
          ok = 1 ; 
          if (arrayMax (usedZones))
            for (j = 0, vp = arrp (usedZones, 0, HIT) ; ok && j < arrayMax (usedZones) ; j++, vp++) 
              {
                if (up->cDNA_clone == vp->cDNA_clone && up->a1 < vp->a2 && up->a2 > vp->a1)
                  { 
                    ok = 0 ; 
                    if (vp->a1 > up->a1)  vp->a1 = up->a1 ;
                    if (vp->a2 < up->a2)  vp->a2 = up->a2 ;
                    break ;
                  }
                if (up->a2 < vp->a1)
                  break ;
              }
          if (ok)
            {
              vp = arrayp (usedZones, arrayMax (usedZones), HIT) ;
              vp->a1 = up->a1 ; vp->a2 = up->a2 ;
              vp->cDNA_clone = up->cDNA_clone ;
            }
        }
    }
  if (0)
    { 
      arrayDestroy(genes) ; 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }
  if (usedZones)
    arraySort (usedZones, cDNAOrderEstByA1) ;

  if (0)
    { 
      arrayDestroy(genes) ; 
      arrayDestroy (geneHits) ;
      invokeDebugger () ;
    }
  cDNAFilterUsedZonesC1 (0, 0, 0, 0, 0) ; /* clean up the static since usedZones array changed */

  arrayDestroy (geneHits) ;
  arrayDestroy (genes) ;
  arrayDestroy (plainGeneHits) ;
  { /* keep for second turn only those reads that have a first round gene */
    HIT *u1p, *u2p ;
    KEYSET rHit = keySetCreate () ;
    int i1, i2 ;

    if (hits && nDoRepeats > 0 && arrayMax (hits))
      for (i1 = 0, u1p = arrp (hits, 0, HIT) ; i1 < arrayMax (hits) ; i1++, u1p++)
        keySetInsert (rHit, u1p->cDNA_clone) ;
    arrayDestroy (hits) ;
    keySetDestroy (newgenes) ;

    if (rHitsAll && arrayMax (rHitsAll) && nDoRepeats > 0)
      { /* keep happy few */
        for (i1 = i2 = 0, u1p = u2p = arrp (rHitsAll, 0, HIT) ; i1 < arrayMax (rHitsAll) ; i1++, u1p++)
          {
            if (!keySetFind (rHit, u1p->cDNA_clone, 0))
              continue ;
            if (i1 > i2) *u2p = *u1p ;
            i2++ ; u2p++ ;
          }
        arrayMax (rHitsAll) = i2 ;
      }
    keySetDestroy (rHit) ;
    hits = arrayCopy (rHitsAll) ;
    cDNAFilterUsedZones (hits, usedZones) ;
  }
 
  if (nDoRepeats-- > 0 && hits && arrayMax (hits) && currentGenes && arrayMax (currentGenes) > rNnewgenes) /* pleaseDoRepeats ()) */
    {
      rNnewgenes = arrayMax (currentGenes) ;
      goto doRepeats ;
    }

  /* run second strand */
  if (direction == 3)
    {
      nDoRepeats = nDoRepeats0 ;
      arrayDestroy (usedZones) ;
      arrayDestroy (rHitsAll) ;
      direction = -1 ;
      goto ladirection ;
    }

  /* destroy what is common to both strands */
  arrayDestroy (linkPos) ;
  arrayDestroy (walls) ;
  arrayDestroy (dnaR) ; /* dna created by getSeqDna */
  arrayDestroy (rHitsAll) ;
  arrayDestroy (usedZones) ;
  arrayDestroy (usedQualities) ;
  arrayDestroy (hits) ;

  /* usually, keep ass, assR, clipTops, estlenghts around as static */
  if (type == 3 || type == 9999)
    {
      keySetDestroy (cDNA_clones) ;
      keySetDestroy (clipTops) ;
      keySetDestroy (clipEnds) ;
      keySetDestroy (trackingClones) ;
      assDestroy (ass) ; assDestroy (assR) ;
      assDestroy (assDifficult) ; assDestroy (assDifficultR) ;
    }
  getSeqDna (KEYMAKE (_VCalcul, 12345)) ; /* dna clean Up */

  return currentGenes ;
} /* alignEst2Cosmid */

/*********************************************************************/

static KEYSET getReads (KEY gene)
{
  KEYSET kClones = queryKey (gene,"Follow read ; Follow from_gene ; { Follow cDNA_clone } $| {Follow Ignored_clone } ; Follow read ; { IS *} $| {>fuse_to_clone} ; ! discarded && Is_read") ;

  if (keySetMax(kClones) > 300)
    printf ("Gene %s %u reads", name(gene), keySetMax(kClones)) ; 
  return kClones ;
}

/*********************************************************************/

void autoAnnotateGene (KEY gene, KEY annot)
{
  KEY key ;
  OBJ Gene, Annot ;
  int n ;

  Gene = bsUpdate(gene) ;
  Annot = bsUpdate (annot) ;
  if (!Gene || !Annot)
    goto done ;
  
  n = 0 ;
  if (bsGetKey (Gene, _cDNA_clone, &key)) 
    do { n++; 
    } while (bsGetKey (Gene, _bsDown, &key)) ;
  if (n == 1)
    bsAddTag (Annot,str2tag("Single_clone")) ;
  else if (bsFindTag (Annot,str2tag("Single_clone")))
    bsRemove (Annot) ;                      
       
  /* 5 prime end */
  if (bsGetKey (Gene, _Transpliced_to, &key)) 
    {
      KEY tag = lexReClass (key, &tag, _VSystem) ;
      if (tag && bsIsTagInObj (Annot, annot, tag))
        bsAddTag (Annot, str2tag(name(key))) ;
      else if (!strcmp(name(key), "gccgtgctc"))
        /* bsAddTag (Annot, str2tag("gccgtgctc")) */ ;
      else 
        if (!bsGetData(Annot, str2tag("motif_5p"), _Text, 0))
          bsAddData (Annot, str2tag("motif_5p"), _Text, name(key)) ;
    }
  /* many clones */
  if (!bsFindTag (Gene, _Total_gap_length))
    {
      if (bsFindTag (Annot, str2tag("Gaps")))
        bsRemove (Annot) ;
      if (bsFindTag (Annot, str2tag("Single_gap")))
        bsRemove (Annot) ;
      bsAddTag (Annot,str2tag("No_Gap")) ;
    }
  else
    {
      if (bsFindTag (Annot, str2tag("No_Gap")))
        bsRemove (Annot) ;
      /* do not add gap, dan does not like it 
         if (!bsFindTag (Annot, str2tag("Single_gap")))
         bsAddTag (Annot,str2tag("Gaps")) ;
         */
    }
  /* 3 prime end */ 
  if (bsFindTag (Gene, _PolyA_after_base))
    bsAddTag (Annot,str2tag("Poly_A_seen_in_5p_read")) ;
  /* translation */
  /* genefinder */
done:
  bsSave(Gene) ;
  bsSave (Annot) ;
}

/*********************************************************************/

static KEYSET alignGene2WalledCosmid (KEY cosmid, KEYSET genes, KEYSET reads,
                                      int z1, int z2, char *nom, BOOL isUp, int doLimit,
                                      int searchRepeats)
{
  Array walls = alignCosmidWalls (cosmid) ;
  HIT *wp ;
  int pass, oldz1, iWall, myz1, myz2 ;
  KEYSET clones = 0, oldr = 0, oldr1 = 0, oldr2 = 0, oldr3 = 0, oldr4 = 0 ;
  Array dna = getSeqDna(cosmid) ;  /* do NOT destroy, do not reuse after calling align */
  int dnaMax = dna ? arrayMax(dna) : 0 ;

  if (!dna)
    goto abort ;

  clones = reads ? query (reads, ">cdna_clone ; ignored_in_gene") : 0 ;
  if (clones)
    {
      int i ;
      OBJ Clone = 0 ;

      for (i = 0 ; i < keySetMax (clones) ; i++)
        if ((Clone = bsUpdate (keySet (clones, i))))
          {
            if (bsFindTag (Clone, str2tag("ignored_in_gene")))
              bsRemove (Clone) ;
            bsSave (Clone) ;
          }
    }

  oldr = keySetCopy (reads) ;
  BOOL composite = TRUE ;
  if (z1 > z2) { int dummy = z1 ; z1 = z2 ; z2 = dummy ; }
  if (z1 < 1) z1 = 1 ;
  if (z2 > dnaMax)
    z2 = dnaMax ;

  for (pass = (doLimit ? 1 : 0) ; pass < 2 ; pass++)
    {
      if (!doLimit) { myz1 = 0 ; myz2 = dnaMax ; }
      else { myz1 = z1 ; myz2 = z2 ; }
      if (pass == 0 && composite)
	continue ;
      oldz1 = myz1 ;
      if (walls)                                 /* forbiden, because it duplicates reads */
        /* hence the walls will still cut the genes but a read can still be aligned just once to the cosmid */
        for (iWall = 0 ; iWall < arrayMax(walls) && keySetMax (oldr) ; iWall++)
          {
            wp = arrayp(walls, iWall, HIT) ;
            if (wp->reverse == isUp &&
                wp->a1 > oldz1 && wp->a1 < myz2)
              {
                if (pass == 1 ||
                    (z1 < wp->a1 && z2 > oldz1))
                  genes = alignEst2Cosmid (cosmid, oldr, genes, 0, 3, oldz1, wp->a1, nom, isUp ? 1 : 2, searchRepeats) ;
                oldz1 = wp->a1 ; 
                oldr1 = oldr ;
                oldr2 = query (genes,"{ Follow cDNA_clone } $| {Follow Ignored_clone } ; FOLLOW Read IS_read") ;
                oldr3 = keySetMINUS (oldr1, oldr2) ;
                oldr4 = query (oldr3, ">cdna_clone ; ! ignored_in_gene ; >read") ;
                oldr = keySetAND (oldr3, oldr4) ;
                keySetDestroy (oldr1) ;
                keySetDestroy (oldr2) ;
                keySetDestroy (oldr3) ;
                keySetDestroy (oldr4) ;
              }
          }
      if (keySetMax (oldr))
        genes = alignEst2Cosmid (cosmid, oldr, genes, 0, 3, oldz1, myz2, nom, isUp ? 1 : 2, searchRepeats) ;
      oldr1 = oldr ;
      oldr2 = query (genes,"{ Follow cDNA_clone } $| {Follow Ignored_clone  } ; FOLLOW Read IS_read") ;
      oldr3 = keySetMINUS (oldr1, oldr2) ;
      oldr4 = query (oldr3, ">cdna_clone ; ! ignored_in_gene ; >read") ;
      oldr = keySetAND (oldr3, oldr4) ;
      keySetDestroy (oldr1) ;
      keySetDestroy (oldr2) ;
      keySetDestroy (oldr3) ;
      keySetDestroy (oldr4) ;
    }
  
 abort:
  alignEst2Cosmid (0, 0, 0, 0, 9999, 0, 0, 0, 0, 0) ; /* clean up */
  arrayDestroy (walls) ;
  keySetDestroy (clones) ;
  keySetDestroy (oldr) ;
  keySetDestroy (oldr1) ;
  keySetDestroy (oldr2) ;
  keySetDestroy (oldr3) ;
  keySetDestroy (oldr4) ;

  return genes ;
} /*alignGene2WalledCosmid */

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

static BOOL cDNAIsGeneUp (KEY cosmid, KEY gene)
{
  OBJ Cosmid = 0 ;
  KEY p1 = 0, p2 = 0 ;
  int a1 = 0, a2 = 0, pass = 0 ;
  BOOL isUp = FALSE, keyFound = FALSE ;

 lao: 
  Cosmid = bsCreate (cosmid) ;
  if (Cosmid && 
      (keyFound = bsFindKey (Cosmid, _Transcribed_gene, gene)) &&
      bsGetData (Cosmid, _bsRight, _Int, &a1) &&
      bsGetData (Cosmid, _bsRight, _Int, &a2) &&
      a1 > a2)
    isUp = TRUE ;
  if (!keyFound && !pass &&
      bsGetKey (Cosmid, str2tag("Parts"), &p1))
    bsGetKey (Cosmid, _bsDown, &p2) ;
  bsDestroy (Cosmid) ;
  pass = 1 ;
  if (!keyFound && p2)
    { cosmid = p2 ; p2 = 0 ; goto lao ; }
  
  if (!keyFound && p1)
    { cosmid = p1 ; p1 = 0 ; goto lao ; }
 
  return isUp ;
}

/*********************************************************************/
/* fuse recursivelly all genes that should be fused
 * keeping as representative the longest one
 */
static int cDNAFuseGeneKeySet (KEYSET tgs) 
{
  KEY map, cosmid, tg, longestTg ;
  KEY _To_be_fused_with = str2tag ("To_be_fused_with") ;
  KEYSET ks0, ks1, ks2, ks3 ;
  KEYSET clones = 0 ;
  int ii, j, dx, dxMax, ntgs ;
  int a1, a2, da1, da2, u1min, u2max, u1, u2 ;
  OBJ Tg = 0, Cosmid = 0 ;
  KEYSET parts = 0 ;

  if (keySetMax (tgs) == 0)
    return 0 ;

  keySetSort (tgs) ;
  ks0 = keySetCopy (tgs) ;
  ntgs = 0 ; 
  for (ii = 0 ; ii < keySetMax(ks0) ; ii++)
    {
      tg = keySet (ks0, ii) ;
      ks2 = queryKeyTransitiveClosure (tg, _To_be_fused_with) ;
      if (!keySetMax (ks2))
        continue ;
      parts = queryKey (tg, ">Genomic_sequence ; >Parts") ;
      if (keySetMax (parts)) /* I am already aligned in a double section */
	ks3 = queryKey (tg, ">Genomic_sequence ; {IS *} SETOR {>parts} ; >Transcribed_gene") ;
      else
	{
	  parts = query (ks2, ">Genomic_sequence ; {!parts} SETOR {>parts}") ;
	  if (keySetMax (parts) > 5) /* m2024_04_11, was : i cannot fuse 3 subsections */
	    ks3 = queryKey (tg, ">Genomic_sequence ; >Transcribed_gene") ;
	  else  /* ok i can fuse everybody */
	    ks3 = keySetCopy (ks2) ;
	}
      keySetDestroy (parts) ;
      ks1 = keySetAND (ks2, ks3) ;
      keySetDestroy (ks2) ;
      keySetDestroy (ks3) ;
      ks2 = ks0 ;
      ks0 = keySetMINUS (ks2, ks1) ;
      keySetDestroy (ks2) ;

      ii = -1 ; /* iterate on the remaining tgs */
      /* find the longest tg in ks1 */
      longestTg = 0 ; a1 = a2 = 0 ; da1 = da2 = 0 ;
      u1min = 1 ; u2max = -1 ;
      for (j = dxMax = 0 ; j < keySetMax (ks1) ; j++)
        {
          tg = keySet (ks1, j) ;
          dx = 0 ; u1 = u2 = 0 ;
          if ((Tg = bsCreate (tg)))
            {
              bsGetData (Tg, _Covers, _Int, &dx) ;
              if (bsGetKey (Tg, _IntMap, &map) &&
                  bsGetData (Tg, _bsRight, _Int, &u1) &&
                  bsGetData (Tg, _bsRight, _Int, &u2)
                  )
                {
                  if (u1min > u2max) u1min = u2max = u1 ;
                  if (u1 < u1min) u1min = u1 ;
                  if (u2 < u1min) u1min = u2 ;
                  if (u1 > u2max) u2max = u1 ;
                  if (u2 > u2max) u2max = u2 ;
                }
              bsDestroy (Tg) ;
            }
          if (j == 0 || dx > dxMax)
            {
              dxMax  = dx ;
              longestTg = tg ;
              a1 = u1 ; a2 = u2 ;
            }
        }
      /* da1, da2 is the needed extention of the gene */
      if (a1 < a2) { da1 = a1 - u1min ; da2 = u2max - a2 ; }
      else if (a1 > a2) { da1 = u2max - a1 ; da2 = a2 - u1min ; }

      keySet (tgs, ntgs++) = longestTg ;
      for (j = 0 ; j < keySetMax (ks1) ; j++)
        {
          if (keySet (ks1, j) == longestTg)
            keySet (ks1, j) = 0 ;
        }
      if (keySetMax (ks1) > 1)
        {
          keySetCompress (ks1) ; /* remove longest tg */
          /* attibute all the clones to the longest tg */
          clones = query (ks1, ">cDNA_clone") ;
          cdnaCleanUp (0, 0, ks1) ;

          cosmid = 0 ;
          Tg = bsUpdate (longestTg) ;
          if (Tg)
            {
              for (j = 0 ; j < keySetMax(clones) ; j++)
                bsAddKey (Tg, _cDNA_clone, keySet (clones, j)) ;
              bsGetKey (Tg, _Genomic_sequence, &cosmid) ;
              if (bsFindTag (Tg, _Assembled_from))
                bsRemove (Tg) ; /* forces recalculation */
              bsSave (Tg) ;
            }
          if (cosmid && (Cosmid = bsUpdate(cosmid)))
            {
              if (bsFindKey (Cosmid, _Transcribed_gene, longestTg) &&
                  bsGetData (Cosmid, _bsRight, _Int, &a1) &&
                  bsGetData (Cosmid, _bsRight, _Int, &a2))
                {
                  if (a1 < a2) { a1 -= da1 ; a2 += da2 ; }
                  else { a1 += da1 ; a2 -= da2 ; }
                  bsAddKey (Cosmid, _Transcribed_gene, longestTg) ;
                  bsAddData (Cosmid, _bsRight, _Int, &a1) ;
                  bsAddData (Cosmid, _bsRight, _Int, &a2) ;
                }              
              bsSave (Cosmid) ;
            }
          keySetDestroy (clones) ;
          for (j = 0 ; j < keySetMax (ks1) ; j++)
            {
              tg = keySet (ks1, j) ;
              if (tg) 
                cdnaCleanUp (0, tg, 0) ; 
	      Tg = bsUpdate (tg) ;
	      if (Tg)
		{
		  if (bsFindTag (Tg, _Assembled_from))
		    bsRemove (Tg) ; 
		  if (bsFindTag (Tg, _cDNA_clone))
		    bsRemove (Tg) ; 
		  bsSave (Tg) ;
		}

            }
        }
      keySetDestroy (ks1) ;
    }
  keySetDestroy (ks0) ;
  keySetMax (tgs) = ntgs ;
  keySetSort (tgs) ;
  keySetCompress (tgs) ;
  return keySetMax (tgs) ;
} /* cDNAFuseGeneKeySet */

/*********************************************************************/

static KEY cDNAFuseGenes (KEY tg) 
{
  KEYSET tgs =  keySetCreate () ;

  keySet (tgs, 0) = tg ;
  if (cDNAFuseGeneKeySet (tgs))
    tg = keySet (tgs, 0) ; /* now i am the longest and i own all the clones */
  keySetDestroy (tgs) ;

  return  tg ;
}

/*********************************************************************/

static int cDNAlimitRepeatedGene (KEY cosmid, KEY gene, int c1, int c2, int *z1p, int *z2p)
{
  int i, z1, z2, a1, a2 = 0, b1, b2 ;
  int foundLimit = 0 ;
  KEY mm, map, gg ;
  KEYSET genes = queryKey (gene, ">cdna_clone ; >from_gene") ;
  OBJ Gene = 0 ;
  
  if (keySetMax (genes) < 2)
    goto done ;

  if ((Gene = bsCreate (gene)))
    {
      if (bsGetKey (Gene, _IntMap, &map) &&
          bsGetData (Gene, _bsRight, _Int, &a1) &&
          bsGetData (Gene, _bsRight, _Int, &a2))
	{} ;
      bsDestroy (Gene) ;
    }
  if (a2 == 0)
    goto done ;

  if (a1 < a2) { z1 = 0 ; z2 = a2 +100000000 ; }
  else { z1 = a1 + 100000 ; z2 = 0 ; }

  /* limit between previous and next repeat in IntMap coords */
  for (i = 0 ; i < keySetMax(genes) ; i++)
    {
      gg = keySet (genes, i) ;
      if (gg == gene)
        continue ;
      if ((Gene = bsCreate (gg)))
        {
          if (bsGetKey (Gene, _IntMap, &mm) &&
              mm == map &&
              bsGetData (Gene, _bsRight, _Int, &b1) &&
              bsGetData (Gene, _bsRight, _Int, &b2))
            {
              if (a1 < a2 && b1 < b2)
                {
                  if (b1 > a2 && b1 < z2) z2 = b1 - 1 ;
                  if (b2 < a1 && b2 > z1) z1 = b2 + 1 ;
                }
              if (a1 > a2 && b1 > b2)
                {
                  if (b1 < a2 && b1 > z2) z2 = b1 + 1 ;
                  if (b2 > a1 && b2 < z1) z1 = b2 - 1 ;
                }
            }
          bsDestroy (Gene) ;
        }
    }
  /* transform z into distance from ends of gene */
  if (a1 < a2) { z1 = a1 - z1 ; z2 = z2 - a2 ; }
  else { z1 = z1 - a1 ; z2 = a2 - z2 ; }
  /* now transform into contig coordinates */
  if (c1 < c2) { *z1p = c1 - z1 ; *z2p = c2 + z2 ; }
  else { *z2p = c1 + z1 ; *z1p = c2 - z2 ; }

  foundLimit = 1 ;

done:
  keySetDestroy (genes) ;
  return foundLimit ;
} /* cDNAlimitRepeatedGene */

/*********************************************************************/

static KEYSET cDNARealignGenePart (KEY cosmid,
                                   KEYSET originalReads,
                                   BOOL isUp, int z1, int z2, 
                                   int doLimit, int searchRepeats)
{
  KEYSET genes = 0, oldr = 0, newr = 0, error = 0 ;
  BOOL debug = FALSE ;

  oldr = arrayCopy (originalReads) ;
  genes = alignGene2WalledCosmid (cosmid, genes, oldr, z1, z2, 0, isUp, doLimit, searchRepeats) ;

  if (genes && debug && keySetMax(genes))
    printf ("Pass 1 %s\n", name(keySet(genes,0))) ;
  BOOL composite = TRUE ;
  if (composite)
    goto done ;
  newr = genes ?
    query (genes,"{ Follow cDNA_clone } $| {Follow Ignored_clone ; IS keepIgnoringThatJunk } ; FOLLOW Read IS_read") :
    keySetCreate () ;
  if (debug) { printf ("succces on first strand:\n") ; showKeySet (newr) ; }
  error = keySetMINUS (oldr, newr) ; keySetDestroy (oldr) ; keySetDestroy (newr) ;
  oldr = error ; error = 0 ; 
  /* if some left over, try without boundaries */
  if (doLimit != 2 && arrayMax(oldr))
    {
      if (debug) showKeySet (oldr) ;
      genes = alignGene2WalledCosmid (cosmid, genes, oldr, 0, 0, 0, isUp, 0, searchRepeats) ; /* no limit */
      if (genes &&  debug && keySetMax(genes))
        printf ("Pass 2 %s\n", name(keySet(genes,0))) ;
      
      newr = genes ?
        query (genes,"ignored Follow cDNA_clone  ; FOLLOW Read IS_read") :
        keySetCreate () ;
      if (debug) { printf ("succces on larger aera same strand:\n") ; showKeySet (newr) ; }
      error = keySetMINUS (oldr, newr) ; keySetDestroy (oldr) ; keySetDestroy (newr) ;
      oldr = error ; error = 0 ; 
    }
  /* if some left over, try the other strand */
  if (doLimit != 2 && arrayMax(oldr))
    {
      if (debug) showKeySet (oldr) ;
      genes = alignGene2WalledCosmid (cosmid, genes, oldr, 0, 0, 0, !isUp, FALSE, searchRepeats) ; 
      if (genes &&  debug && keySetMax(genes))
        printf ("Pass 3 %s\n", name(keySet(genes,0))) ;
      newr = genes ?
        query (genes,"ignored Follow cDNA_clone  ; FOLLOW Read IS_read") :
        keySetCreate () ;
      if (debug) { printf ("succces on larger other strand:\n") ; showKeySet (newr) ; }
      error = keySetMINUS (oldr, newr) ; keySetDestroy (oldr) ; keySetDestroy (newr) ;
      oldr = error ; error = 0 ; 
    }
 done:
  /* look for dead mrna */
  
  keySetDestroy (oldr) ;
  keySetDestroy (error) ;
  keySetDestroy (newr) ;

  return genes ;
}  /* cDNARealignGenePart */

/*********************************************************************/
/*********************************************************************/
BOOL gcpDebug = FALSE ;

typedef struct gcpStruct { KEY mrna ; Array spl ; Array exons ; Array clones ; BOOL isReversed ; } GCP ;

static void showGCP (Array aa)
{
  int ii, i ;
  GCP *gcp ;

  if (arrayExists (aa))
    for (ii = 0 ; ii < arrayMax (aa) ; ii++)
      {
        gcp = arrp (aa, ii, GCP) ;
        if (!gcp->clones)
          { 
            if (0) fprintf (stderr, "%d : Empty\n", ii) ;
            continue ;
          }
        
        fprintf (stderr, "%d : %s\n", ii, name(gcp->mrna)) ;
        fprintf (stderr, "\t%u spl::", gcp->spl ? keySetMax (gcp->spl) : 0) ;
        for (i = 0 ; gcp->spl && i < keySetMax (gcp->spl) ; i++)
          fprintf (stderr, " %u",  keySet (gcp->spl, i)) ;
        fprintf (stderr, "\n") ;
        fprintf (stderr, "\t%u clones::", gcp->clones ? keySetMax (gcp->clones) : 0) ;
        for (i = 0 ; gcp->clones && i < 6 && i < keySetMax (gcp->clones) ; i++)
          fprintf (stderr, " %s",  name(keySet (gcp->clones, i))) ;
        fprintf (stderr, "\n") ;
      }
} /* showGCP */

/*********************************************************************/
/* get the reads which have changed starnd and are currently backwards */
static void cDNARealignGetReversedConnectedPart (KEY gene, Array aa, BOOL isBack)
{
  OBJ Gene = bsCreate (gene) ;
  GCP *gcp = 0 ;
  Array units = arrayCreate (32000, BSunit) ;
  BSunit *uu ;
  BOOL isDown ;
  KEY est, clone ;
  KEYSET clones = 0 ;
  int ii, a1, a2, x1, x2 ;
  
  if (Gene && bsGetArray (Gene, _Assembled_from, units, 5))
    {
      for (ii = 0 ; ii < arrayMax (units) ; ii += 5)
	{
	  uu = arrayp (units, ii, BSunit) ;
	  a1 = uu[0].i ; a2 = uu[1].i ;
	  est = uu[2].k ;
	  x1 = uu[3].i ; x2 = uu[4].i ;
	  isDown = ((a1 < a2 && x1 < x2) || (a1 > a2 && x1 > x2)) ;
	  if (isBack) isDown = !isDown ;
	  /* use a negative query to recover as ok the reads
	   * where the orienation is not specified
	   */
	  if ((!keyFindTag (est, _Reverse) && isDown) ||
	      (!keyFindTag (est, _Forward) && !isDown)
	      )
	    {
	      clone = keyGetKey (est, _cDNA_clone) ;
	      if (bsFindKey (Gene, _cDNA_clone, clone))
		{
		  if (!gcp)
		    gcp = arrayp (aa, arrayMax(aa), GCP) ;
		  if (!gcp->clones)
		    gcp->clones = keySetCreate () ;
		  gcp->isReversed = isBack ;
		  keySet (gcp->clones, keySetMax(gcp->clones)) = clone ;
		}	  
	    }
	}
      if (gcp && gcp->clones)
	{
	  keySetSort (gcp->clones) ;
	  keySetCompress (gcp->clones) ;
	}
      /* add the floating clones */
      clones = queryKey (gene,">cDNA_clone;COUNT{>read from_gene}==0") ;
      if (keySetMax (clones))
	{
	  if (! gcp)  
	    gcp = arrayp (aa, arrayMax(aa), GCP) ;
	  if (!gcp->clones)
	    gcp->clones = clones ;
	  else
	    {
	      KEYSET kstmp = keySetOR (gcp->clones, clones) ;
	      keySetDestroy (clones) ;
	      keySetDestroy (gcp->clones) ;
	      gcp->clones = kstmp ;
	    }
	}
    }
  else if (Gene && !isBack && (clones = queryKey (gene,">cDNA_clone")))
    {
      if (!gcp)
	gcp = arrayp (aa, arrayMax(aa), GCP) ;
      gcp->clones = clones ;
    }

  arrayDestroy (units) ;
  bsDestroy (Gene) ;
} /* cDNARealignGetReversedConnectedPart */

/*********************************************************************/

static BOOL cDNARealignGetOneConnectedPart (KEY gene, KEY mrna, GCP *gcp, BOOL isFirst)
{
  OBJ Gene = bsCreate (gene) ;
  OBJ Mrna = 0 ;
  int m1, m2, i, jj = 0, jjx = 0 ; /* coords of mrna in gene */
  BSunit *uu ;
  KEY tsl = 0 ;
  static Array aa = 0 ;
  
  gcp->mrna = mrna ;
  gcp->spl = keySetCreate () ;
  gcp->exons = keySetCreate () ;
  if (bsFindKey (Gene, _mRNA, mrna) && 
      bsGetData (Gene, _bsRight, _Int, &m1) &&
      bsGetData (Gene, _bsRight, _Int, &m2) &&
      (Mrna = bsCreate (mrna)) &&
      (aa = arrayReCreate (aa, 120, BSunit)) &&
      bsGetArray (Mrna, _Splicing, aa, 5)
      ) 
    {  /* SL1 counts as an intron boundary */
      if (bsGetKey (Mrna, _Transpliced_to, &tsl)
          && lexstrcmp ("sl0", name(tsl)) < 0)
        keySet (gcp->spl, jj++) = m1 ;
      for (i = 0 ; i < arrayMax(aa) ; i += 5) 
      {
        uu = arrp (aa, i, BSunit) ;
        if (strstr(name(uu[4].k), "ntron"))
          {
            keySet (gcp->spl, jj++) = m1 + uu[0].i ;
            keySet (gcp->spl, jj++) = m1 + uu[1].i ;
          }        
        if (strstr(name(uu[4].k), "xon"))
          {
            keySet (gcp->exons, jjx++) = m1 + uu[0].i ;
            keySet (gcp->exons, jjx++) = m1 + uu[1].i ;
          }        
      }
    }
  bsDestroy (Gene) ;
  bsDestroy (Mrna) ;

  if (isFirst)
    gcp->clones = queryKey (mrna
                            , "({>cDNA_clone} $& {>from_gene;>cdna_clone}) $| {>cdna_clone ; >Fuse_to_clone} $| {>from_gene ; >cdna_clone ; ! In_mrna} ") ;
  else
    gcp->clones = queryKey (mrna, "({>cDNA_clone}  $& {>from_gene;>cdna_clone}) $| {>cdna_clone ; >Fuse_to_clone}") ;
  return TRUE ;
} /* cDNARealignGetOneConnectedPart */

/*********************************************************************/
/* 2 mrnas are connected if they share a clone or an exon boundary 
 * or if they touch and have no introns
 * or if they overlap by 2/3 or shortest even with introns 
 */
static BOOL cDNARealignIsConneted (Array aa, int n1, int n2)
{
  GCP *gcp1 , *gcp2 ;
  KEYSET ks1 = 0, ks2 = 0 ;
  BOOL ok = FALSE ;
  int ii, jj, ln1, ln2, ln3 ;
  int x1, x2, y1, y2, z1, z2 ;
  BOOL b1, b2 ;

  gcp1 = arrp (aa, n1, GCP) ;
  gcp2 = arrp (aa, n2, GCP) ;

  ks1 = keySetAND (gcp1->spl, gcp2->spl) ;
  if (keySetMax (ks1))
    ok = TRUE ; /*  the 2 mrnas share a donor or an acceptor or an SL1 acceptor site */
  else
    {
      ks2 = keySetAND (gcp1->clones, gcp2->clones) ;
      if (keySetMax (ks2))
        ok = TRUE ; /* the 2 mrnas share a clone */
    }
  if (!ok && gcp1->spl && keySetMax (gcp1->spl) && gcp2->spl && keySetMax (gcp2->spl))
    {
      /* both sets have introns, check the exon overlap */

      /* get the length of each */
      for (ii = ln1 = 0 ; ii < keySetMax (gcp1->exons) ; ii += 2)
        ln1 += keySet (gcp1->exons, ii + 1) - keySet (gcp1->exons, ii) + 1 ;
      for (ii = ln2 = 0 ; ii < keySetMax (gcp2->exons) ; ii += 2)
        ln2 += keySet (gcp2->exons, ii + 1) - keySet (gcp2->exons, ii) + 1 ; 

      /* get the length of the intersect */
      ln3 = ii = jj = 0 ;
      x1 = x2 = y1 = y2 = 0 ;
      b1 = b2 = FALSE ;
      while (1)
        {
          if (!b1 && ii < keySetMax (gcp1->exons))
            { b1 = TRUE ; x1 = keySet (gcp1->exons, ii++) ; x2 = keySet (gcp1->exons, ii++) ; }
          if (!b2 && jj < keySetMax (gcp2->exons))
            { b2 = TRUE ; y1 = keySet (gcp2->exons, jj++) ; y2 = keySet (gcp2->exons, jj++) ; }
          if (! b1 || !b2)
            break ;
          z1 = x1 > y1 ? x1 : y1 ;
          z2 = x2 < y2 ? x2 : y2 ;
          if (z1 <= z2) ln3 += z2 - z1 + 1 ; /* intersect of the 2 exons */
          if (x2 < y2) b1 = FALSE ;
          else b2 = FALSE ;
        }
      /* compare the intersect to each part */
      if (3 * ln3 > 1 * ln1 || 3 * ln3 > 1 * ln2)
        ok = TRUE ; /* 2 spliced mrnas must share 1/3 of the shortest sequence 
                     * was 2/3 before 2007_02_10 
                     */
    }

  if (!ok && /* added 2007_02_10 : spliced or not, check clear overlap */
      gcp1->exons &&
      gcp2->exons
      )
    {
      /* check the genebox overlap */
      int a1, a2, b1, b2, c1, c2 ;

      /* get the length of each */
      a1 =  keySet (gcp1->exons, 0) ;
      b1 =  keySet (gcp2->exons, 0) ;
      ii = keySetMax (gcp1->exons) - 1;
      a2 =  keySet (gcp1->exons, ii) ;
      ii = keySetMax (gcp2->exons) - 1;
      b2 =  keySet (gcp2->exons, ii) ;

      ln1 = a2 - a1 + 1 ;
      ln2 = b2 - b1 + 1 ;
      c1 = a1 > b1 ? a1 : b1 ;
      c2 = a2 < b2 ? a2 : b2 ;
      ln3 = c2 - c1 ;

      /* compare the intersect to each part */
      if (3 * ln3 > 2 * ln1 || 3 * ln3 > 2 * ln2)
        ok = TRUE ; /* added 2007_02_10 : gene boxes must be clearly overlapping */
    }

  keySetDestroy (ks1) ;
  keySetDestroy (ks2) ;
  
  return ok ;
} /* cDNARealignIsConneted */

/*********************************************************************/

static BOOL cDNARealignFuseConnected (Array aa, int n1, int n2)
{
  GCP *gcp1 , *gcp2 ;
  KEYSET ks1 = 0, ks2 = 0, ks3 = 0 ;
  int ii, jj, kk, x1, x2, y1, y2 ;
  BOOL b1, b2 ;
  gcp1 = arrp (aa, n1, GCP) ;
  gcp2 = arrp (aa, n2, GCP) ;
  
  ks1 = keySetOR (gcp1->spl, gcp2->spl) ;
  ks2 = keySetOR (gcp1->clones, gcp2->clones) ;
  ks3 = keySetCreate () ;
  /* i cannot OR the 2 exons keyset because i must see the start stops in each */
  ii = jj = kk = 0 ;
  x1 = x2 = y1 = y2 = 0 ;
  b1 = b2 = FALSE ;
  while (1)
    {
      if (!b1 && ii < keySetMax (gcp1->exons))
        { b1 = TRUE ; x1 = keySet (gcp1->exons, ii++) ; x2 = keySet (gcp1->exons, ii++) ; }
      if (!b2 && jj < keySetMax (gcp2->exons))
        { b2 = TRUE ; y1 = keySet (gcp2->exons, jj++) ; y2 = keySet (gcp2->exons, jj++) ; }
      if (! b1 && !b2)
        break ;
      if (b1 && (!b2 || x1 < y1))
        {
          if (!b2 || x2 < y1) /* register x1 exon */
            {
              keySet (ks3, kk++) = x1 ; 
              keySet (ks3, kk++) = x2 ; 
              b1 = FALSE ;
            }
          else /* register the union */
            {
              keySet (ks3, kk++) = x1 ; 
              keySet (ks3, kk++) = x2 > y2 ? x2 : y2 ; 
              b1 = b2 = FALSE ;
            }
        }
      else
        {
          if (!b1 || y2 < x1) /* register y1 exon */
            {
              keySet (ks3, kk++) = y1 ; 
              keySet (ks3, kk++) = y2 ; 
              b2 = FALSE ;
            }
          else /* register the union */
            {
              keySet (ks3, kk++) = y1 ; 
              keySet (ks3, kk++) = x2 > y2 ? x2 : y2 ; 
              b1 = b2 = FALSE ;
            }
        }
    }

  keySetDestroy (gcp1->clones) ;
  keySetDestroy (gcp2->clones) ;
  keySetDestroy (gcp1->exons) ;
  keySetDestroy (gcp2->exons) ;
  keySetDestroy (gcp1->spl) ;
  keySetDestroy (gcp2->spl) ;

  gcp1->spl = ks1 ; ks1 = 0 ;
  gcp1->clones = ks2 ; ks2 = 0 ;
  gcp1->exons = ks3 ; ks3 = 0 ;

  return TRUE ;
} /* cDNARealignFuseConnected */

/*********************************************************************/

static Array cDNARealignGetConnectedParts (KEY gene, BOOL splitCloud)
{
  KEYSET mrnas = queryKey (gene, " > mrna") ;
  int n1, n2, n3, nn = keySetMax (mrnas) ;
  Array aa = arrayCreate (nn+1, GCP) ;
  GCP *gcp ;

  if (!splitCloud)
    {
      cDNARealignGetReversedConnectedPart (gene, aa, FALSE) ;      
    }
  else
    {
      for (nn = 0 ; nn < keySetMax (mrnas) ; nn++)
	{ /* for each mrna get the splice junctions and clone list */
	  gcp = arrayp (aa, nn, GCP) ;
	  cDNARealignGetOneConnectedPart (gene, keySet (mrnas, nn), gcp, nn == 0 ? TRUE : FALSE ) ;
	}
      if (0) showGCP (aa) ; 
      for (n2 = 1 ; n2 < arrayMax (aa) ; n2++)
	{
	  gcp = arrp (aa, n2, GCP) ;
	  if (! gcp->clones) continue ;
	  /* search if n2 is connected to n1 above */
	  for (n1 = 0 ; n1 < n2 ; n1++)
	    if (cDNARealignIsConneted (aa, n1, n2))
	      {
		/* other also touching n2 also fuse to n1 */
		for (n3 = n1 + 1 ; n3 < n2 ; n3++)
		  if (cDNARealignIsConneted (aa, n3, n2))
		    cDNARealignFuseConnected (aa, n1, n3) ;
		cDNARealignFuseConnected (aa, n1, n2) ;
	      }
	}
      keySetDestroy (mrnas) ;
      /* keep happy few */
      {
	int ii, jj ;
	GCP *up, *vp ;
	for (ii = jj = 0, up = vp = arrp (aa, 0, GCP) ; ii < arrayMax (aa) ; up++, ii++)
	  {
	    if (! up->clones)
	      continue ;
	    if (vp < up) 
	      *vp = *up ;
	    vp++; jj++ ; 
	  }
	arrayMax (aa) = jj ;
      }
    }
  cDNARealignGetReversedConnectedPart (gene, aa, TRUE) ;
  return aa ;
}  /* cDNARealignGetConnectedParts */

/*********************************************************************/
/* 2023_12_01 :  splitClup == 4: splitOnGeneId
 */

static Array cDNARealignGetGeneIdParts (KEY gene)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *errors = 0 ;
  int nId = 0, nM = 0, nGCP = 0 ;
  AC_DB db = ac_open_db (0, &errors) ; /* local database, cannot fail */
  AC_OBJ Gene = ac_get_obj (db, className (gene), name (gene), h) ;
  AC_TABLE ids = 0, tbl = 0 ;
  Array aa = 0, aaM, aaIds ;
  KEY typeMax = 0 ;
  const char *ccp ;
  GCP *gcp ;
  SPL *up ;
  int tg1, tg2 ;
  ccp = "select tg,id from tg in @, g in tg->gene, pg in g->genefinder, id  in pg->GeneId " ;
  ids = ac_obj_bql_table (Gene, ccp, 0, &errors, h) ; 
  if (! ids || ids->rows < 2)
    goto done ;
   
  aa = arrayHandleCreate (ids->rows, GCP, 0) ; 
  aaIds = arrayHandleCreate (ids->rows, SPL, h) ; 

  /* get the genome coordinate of each geneid */
  ccp = "select id,chrom,a1,a2, tg1, tg2 from tg in @, g in tg->gene, pg in g->genefinder, id in pg->GeneId, chrom in pg->intmap, a1 in chrom[1], a2 in chrom[2], tchrom in tg->intmap, tg1 in tchrom[1], tg2 in tchrom[2]" ;
  tbl = ac_obj_bql_table (Gene, ccp, 0, &errors, h) ; 
  for (int ir = 0 ; ir < tbl->rows ; ir++)
    {
      int a1 = 0, a2 = 0 ;
      KEY id = ac_table_key (tbl, ir, 0, 0) ;
      for (int jr = ir ; jr < tbl->rows ; jr++)
	{
	  int b1, b2 ;
	  KEY id2 = ac_table_key (tbl, jr, 0, 0) ;
	  if (id2 != id)
	    break ;
	  b1 = ac_table_int (tbl, jr, 2, 0) ;
	  b2 = ac_table_int (tbl, jr, 3, 0) ;
	  if (b1 > b2) 
	    { int b0 = b1 ; b1 = b2 ; b2 = b0 ; }
	  tg1 = ac_table_int (tbl, jr, 4, 0) ;
	  tg2 = ac_table_int (tbl, jr, 5, 0) ;
	  if (tg1 > tg2) 
	    { int tg0 = tg1 ; tg1 = tg2 ; tg2 = tg0 ;}

	  if (a1 == 0)
	    { a1 = b1 ; a2 = b2 ; continue ; }
	  if (a1 > b1) a1 = b1 ;
	  if (a2 < b2) a2 = b2 ;
	}
      up = arrayp (aaIds, nId++, SPL) ;
      up->gene = id ; up->a1 = a1 ; up->a2 = a2 ;
    }
  /* merge overlapping Ids and extend them to the end of the genes */
  if (nId > 1)
    {
      BOOL ok = TRUE ;
      arraySort (aaIds, splOrder) ;
      while (ok)   /* merge the geneid that overlap by half size */
	{
	  ok = FALSE ;
	  for (int ii = 0 ; ! ok && ii < nId && ii < sizeof(KEY) ; ii++)
	    {
	      SPL *up = arrp (aaIds, ii, SPL) ;
	      if (up->gene == 0)
		continue ;
	      for (int jj = ii + 1 ; !ok && jj < nId ; jj++)
		{
		  SPL *vp = arrp (aaIds, jj, SPL) ;
		  if (vp->gene == 0)
		    continue ;
		  if (up->a2 > (vp->a1 + vp->a2)/2) /* fuse */
		    { 
		      if (up->a2 < vp->a2) 
			up->a2 = vp->a2 ;
		      vp->gene = 0 ;
		      ok = TRUE ;
		    }
		  if ((up->a1 + up->a2)/2 > vp->a1) /* fuse */
		    { 
		      if (up->a2 < vp->a2) 
			up->a2 = vp->a2 ;
		      vp->gene = 0 ;
		      ok = TRUE ;
		    }
		}
	    }
	}
      ok = TRUE ;
      while (ok)   /* shorten geneid that overlap a little */
	{
	  ok = FALSE ;
	  for (int ii = 0 ; ! ok && ii < nId && ii < sizeof(KEY) ; ii++)
	    {
	      SPL *up = arrp (aaIds, ii, SPL) ;
	      if (up->gene == 0)
		continue ;
	      for (int jj = ii + 1 ; !ok && jj < nId ; jj++)
		{
		  SPL *vp = arrp (aaIds, jj, SPL) ;
		  if (vp->gene == 0)
		    continue ;
		  if (up->a2 > vp->a1)
		    { 
		      int x = vp->a1 ;
		      vp->a1 = up->a2 ;
		      up->a2 = x ;
		      ok = TRUE ;
		    }
		}
	    }
	}
      arraySort (aaIds, splOrder) ;
      if (1)
	{
	  int jj = 0 ;
	  SPL *up = arrp (aaIds, 0, SPL) ;
	  SPL *vp = arrp (aaIds, 0, SPL) ;
	  for (int ii = 0 ; ii < arrayMax (aaIds) ; ii++, up++)
	    {
	      if (up->gene)
		{
		  if (jj < ii) {*vp = *up ; }
		  vp++ ; jj++ ;
		}
	    }
	  nId = arrayMax (aaIds) = jj ;
	}
    }
  if (nId < 2)
    goto done ;
  else
    {  /* expand the first and last geneid to the end of the gene */
      SPL *up = arrp (aaIds, 0, SPL) ;
      if (up->a1 > tg1) up->a1 = tg1 ;
      up = arrp (aaIds, nId - 1, SPL) ;
      if (up->a2 < tg2) up->a2 = tg2 ;
    }
  /* get the genome coordinate of each mrna */
  ccp = "select m, chrom, a1,a2 from tg in @, m in tg->mrna, chrom in m->intmap, a1 in chrom[1], a2 in chrom[2]" ;
  tbl = ac_obj_bql_table (Gene, ccp, 0, &errors, h) ; 
  aaM = arrayHandleCreate (1+tbl->rows, SPL, h) ; 
  for (int ir = 0 ; ir < tbl->rows ; ir++)
    {
      int a1 = 0, a2 = 0 ;
      KEY id = ac_table_key (tbl, ir, 0, 0) ;
      BOOL isBack = FALSE ;
      for (int jr = ir ; jr < tbl->rows ; jr++)
	{
	  int b1, b2 ;
	  KEY id2 = ac_table_key (tbl, jr, 0, 0) ;
	  if (id2 != id)
	    break ;
	  b1 = ac_table_int (tbl, jr, 2, 0) ;
	  b2 = ac_table_int (tbl, jr, 3, 0) ;
	  if (b1 > b2) { int b0 = b1 ; b1 = b2 ; b2 = b0 ; isBack = TRUE ; }
	    {
	      if (a1 == 0)
		{ a1 = b1 ; a2 = b2 ; continue ; }
	      if (a1 > b1) a1 = b1 ;
	      if (a2 < b2) a2 = b2 ;
	    }
	}
      up = arrayp (aaM, nM++, SPL) ;
      up->isBack = isBack ;
      up->gene = id ; up->a1 = a1 ; up->a2 = a2 ;
      /* check if mrna intersects the NM */
      for (int ii = 0 ; ii < nId && ii < sizeof(KEY) ; ii++)
	{
	  SPL *vp = arrp (aaIds, ii, SPL) ;
	  if (up->a1 < vp->a2 - 50 && up->a2 > vp->a1 + 50)
	    up->type |= (1<<ii) ;
	  if (up->type > typeMax)
	    typeMax = up->type ;
	}
    }

  /* now we create the GCP with all mrna of a given type */
  for (KEY type = 0 ; type <= typeMax ; type++)
    {
      for (int ii = 0 ; ii < nM ; ii++)
	{
	  up = arrayp (aaM, ii, SPL) ;
	  if (up->type == type)
	    {
	      AC_OBJ Mrna = ac_get_obj (db, className (up->gene), name (up->gene), h) ;
	      AC_TABLE clones = ac_obj_bql_table (Mrna, ">cdna_clone", 0, &errors, 0) ; 
	      
	      if (clones)
		{
		  gcp = arrayp (aa, type, GCP) ;
		  gcp->isReversed = up->isBack ;
		  if (! gcp->clones)
		    gcp->clones = keySetCreate () ;
		  for (int i = 0 ; i < clones->rows ; i++)
		    keySet (gcp->clones, keySetMax(gcp->clones)) = ac_table_key (clones, i, 0, 0) ;
		}
	      ac_free (Mrna) ;
	      ac_free (clones) ;
	    }
	}
    }

  /* compress the GCP */
  nGCP = arrayMax (aa) ;
  if (nGCP > 0)
    {
      int jj = 0 ;
      GCP *up = arrp (aa, 0, GCP) ; 
      GCP *vp = arrp (aa, 0, GCP) ; 
      for (int ii = 0 ; ii < arrayMax (aa) ; up++, ii++)
	{
	  if (up->clones)
	    {
	      if (ii > jj)
		vp->clones = up->clones ;
	      keySetSort (vp->clones) ;
	      keySetCompress (vp->clones) ;
	      jj++ ; vp++ ;
	    }
	}
      nGCP = arrayMax (aa) = jj ;
    }
 
 done:
  ac_free (db) ;
  ac_free (h) ;
  return aa ;
}  /* cDNARealignGetGeneIdParts */

/*********************************************************************/
/* cut the clones of this genes in subsets corresponding to 
 * connected groups of mRNAs: sharing at least a clone or a boundary
 * 2007_02_10 or touching if no intron, or touching by 1/3 if introns
 * 2017_02_07 splitCloud == 4 split genes with 2 geneids
 */
static KEYSET cDNARealignCutPart (KEY cosmid, KEY gene, 
                               int z1, int z2, int direction,
                               int doLimit, int searchRepeats,
                               int splitCloud)
{
  BOOL isUp = FALSE ;
  int old = OLIGO_LENGTH ;
  KEYSET genes = 0, ksDebug = 0 ;
  BOOL debug = FALSE ;
  static int recursionLevel = 0 ;

  recursionLevel++ ;
  OLIGO_LENGTH = 12 ; 
  switch (direction)
    {
    case 1:
      isUp = TRUE ;
      break ;
    case 2:
      isUp = FALSE ; 
      break ;
    default:
      isUp = cDNAIsGeneUp (cosmid, gene) ;
      break ;
    }


  if (recursionLevel == 1 || splitCloud) /* == 2, always recompute, == 1 only if (!read OR will drop pieces) */
    {
      Array aa = 0 ;
      KEYSET originalReads = 0, genes1 = 0, genes2 = 0 ;
      GCP *gcp ;
      int ii ;
      
      if (splitCloud == 4)
	{
	  aa = cDNARealignGetGeneIdParts (gene) ;
	  if (! aa || arrayMax (aa) < 2)
	    goto done ;
	}
      else
	aa = cDNARealignGetConnectedParts (gene, splitCloud) ;

      if (! keyFindTag (gene, _Assembled_from) ||
          arrayMax(aa)
          )
        {
          cdnaCleanUp (0, gene, 0) ;
          if (gcpDebug) showGCP (aa) ;                         

          for (ii = 0, gcp = arrp (aa, 0, GCP) ; ii < arrayMax (aa) ; gcp++, ii++)
            if (gcp->clones)
              {
                originalReads = query (gcp->clones, ">Read ! discarded && Is_read") ;
                genes1 = cDNARealignGenePart (cosmid, originalReads, gcp->isReversed ? !isUp : isUp, z1, z2, doLimit, searchRepeats) ;
                if (debug && genes1)
                  {
                    ksDebug = query (genes1, ">cdna_clone") ;
                    showKeySet (ksDebug) ;
                    keySetDestroy (ksDebug) ;
                  }
                keySetDestroy (originalReads) ;
                keySetDestroy (gcp->spl) ;
                keySetDestroy (gcp->exons) ;
                keySetDestroy (gcp->clones) ;
                if (genes)
                  {
                    genes2 = genes ;
                    genes = keySetOR (genes1, genes2) ;
                    keySetDestroy (genes1) ;
                    keySetDestroy (genes2) ;
                  }
                else
                  { genes = genes1 ; genes1 = 0 ; }
              }
        }
      else
        gene = 0 ; /* avoid destruction */
      arrayDestroy (aa) ;
    }
  else if (0)
    {
      KEYSET originalReads =  getReads (gene) ;

      if (debug && gene)
        {
          ksDebug = queryKey (gene, ">cdna_clone") ;
          printf ("tg %s has %u clones\n", name(gene), keySetMax (ksDebug)) ;
          showKeySet (ksDebug) ;
          keySetDestroy (ksDebug) ;
        }
      
      cdnaCleanUp (0, gene, 0) ;
      genes = cDNARealignGenePart (cosmid, originalReads, isUp, z1, z2, doLimit, searchRepeats) ;
      keySetDestroy (originalReads) ;
      if (debug && genes)
        {
          ksDebug = query (genes, ">cdna_clone") ;
          showKeySet (ksDebug) ;
          keySetDestroy (ksDebug) ;
        }
    }
  /* this part is to kill gene if that name was not reused */
 
  if (genes)
    {
      int i ;
      
      for (i = 0 ; i < keySetMax(genes) ; i++)
        {
          if (gene == keySet(genes, i)) 
            { gene = 0 ; break ; }
        }
      if (! isComposite)
	{
	  mrnaAnalyseNeighbours (genes) ;
	  mrnaSplitDoubleGenes (genes) ;
	}
    }
  if (gene)
    {
      OBJ Gene ;
      if ((Gene = bsUpdate(gene)))
        bsKill (Gene) ;
    }

 done:  
  OLIGO_LENGTH = old ;
  recursionLevel-- ;
  return genes ;
} /* cDNARealignGeneCutPart */

/*********************************************************************/
/* doLimit = 0 => please fuse and compute inside the current double-tile
 * doLimit = 1 => please compute using limits z1 z2
 * doLimit = 2 => please compute using current (fused) gene size 
 *
 * searchRepeats = 1 => align recursivelly from bottom of cosmid up
 * searchRepeats = 2 => rubber, repeat with increased inron cost
 *  splitCloud == 4 => splitGeneId
 */
KEY cDNARealignGene (KEY gene, int z1, int z2, int direction,
                     BOOL doFuse, int doLimit, int  searchRepeats, int splitCloud)
{
  KEY cosmid = 0, newGene = 0 ;
  int a1, a2 ;
  OBJ Cosmid = 0 ;
  static int nRecursion = 0 ;

  if (nRecursion > 1)
    return gene ; 
  
  if (!sessionGainWriteAccess())
    return gene ;  /* stallmate, nothing changed */
  cDNAAlignInit() ;
  nRecursion++ ;
  if (doFuse &&  
      keyFindTag (gene, str2tag("To_be_fused_with")))
    {
      gene = cDNAFuseGenes (gene) ;
    }        
  doFuse = FALSE ; /* since we just did it */
  cosmid = keyGetKey (gene, _Genomic_sequence) ;
  if ((Cosmid = bsCreate(cosmid)))
    {
      if (bsFindKey (Cosmid, _Transcribed_gene, gene) &&
          bsGetData (Cosmid, _bsRight, _Int, &a1))
        bsGetData (Cosmid, _bsRight, _Int, &a2) ;
      bsDestroy (Cosmid) ;
    }
  /* limit to current size +- 500 bp to allow vector clip extension */
  if (doLimit == 2)  /* locally */
    {
      int dz = a1 < a2 ? 500 : -500 ;
      direction = a1 < a2 ? 2 : 1 ;
      z1 = a1 - dz  ; z2 = a2 + dz ;
    }
  if (!keyFindTag (cosmid, str2tag("Is_gene_tile")) &&
      !keyFindTag (cosmid, str2tag("Junction")))
    {
      /* i prefer to assemble into a junction */
      int rx, lx ;
      KEY leftCosmid = 0, rightCosmid = 0 ;
      
      lx = rx = -1 ;
      /* find the relative position of the left and right cosmid */
      if ((Cosmid = bsCreate(cosmid)))
        {
          if (bsGetKey (Cosmid, str2tag ("Overlap_right"), &rightCosmid) &&
              bsGetData (Cosmid, _bsRight, _Int, &rx)) ;
          else 
            rightCosmid = 0 ; /* reset to zero if rx not set */
          bsGetKey (Cosmid, str2tag ("Overlap_left"), &leftCosmid) ;
          bsDestroy (Cosmid) ;
        }
      if (leftCosmid && (Cosmid = bsCreate (leftCosmid)))
        {
          if (bsFindKey (Cosmid, str2tag ("Overlap_right"), cosmid) &&
              bsGetData (Cosmid, _bsRight, _Int, &lx)) ;
          else 
            leftCosmid = 0 ; /* reset to zero if rx not set */
          bsDestroy (Cosmid) ;
        }
      if ((leftCosmid && lx > 0) || (rightCosmid && rx > 0))
        {
          int choice = 0 ;
          if (!(leftCosmid && lx > 0))
            choice = 2 ;
          else if (!(rightCosmid && rx > 0))
            choice = 1 ;
          else if (a1 < a2)
            {
              if (a1 > rx - a2) /* prefer right junction */
                choice = 2 ;
              else
                choice = 1 ;
            }
          else
            {
              if (a2 > rx - a1) /* prefer right junction */
                choice = 2 ;
              else
                choice = 1 ;
            }

          if (choice == 1)
            {
              KEY junction = 0 ;
              if (lexword2key (messprintf("%s_%s", name(leftCosmid), name(cosmid)), &junction, _VSequence) &&
                  keyFindTag (junction, str2tag("Junction")))
                {
                  cosmid = junction ;
                  if (z1 || z2) { z1 += lx - 1 ; z2 += lx - 1 ; }
                  a1 += lx - 1 ; a2 += lx - 1 ;
                }
            }
          else if (choice == 2)
            {
              KEY junction = 0 ;
              if (lexword2key (messprintf("%s_%s", name(cosmid), name(rightCosmid)), &junction, _VSequence) &&
                  keyFindTag (junction,  str2tag("Junction")))
                {
                  cosmid = junction ;
                }
            }
        }
    }
  {
    KEYSET reads = getReads (gene) ;
    if (! keySetMax (reads))
      {
	cdnaCleanUp (0, gene, 0) ;
	keySetDestroy (reads) ;
	goto done ;
      }
    keySetDestroy (reads) ;
  }
  {
    KEYSET reads = getReads (gene) ;
    KEYSET greens2 = 0, greens = keySetCreate () ;
    OBJ Gene ;
    BOOL fix = FALSE ;
    
    for (int i = 0 ; i < keySetMax (reads) ; i++)
      {
	KEY green = keySet (reads, i) ;
	const char *nam = name (green) ;

	if (! strncmp (nam, "XG_", 3) || ! strncmp (nam, "XG_", 3))
	  {
	    fix |= mrnaDesignCutOnePreMrna (gene, green, greens) ;
	  }
	else
	  keySet (greens, keySetMax (greens)) = green ; 
      }
    keySetSort (greens) ;
    greens2 = query (greens, ">cDNA_clone") ;
    if (fix && (Gene = bsUpdate (gene)))
      {
	KEY tag = _bsRight ;
	bsFindTag (Gene, _cDNA_clone) ;
	bsRemove (Gene) ;
	bsAddTag (Gene, _cDNA_clone) ;
	for (int i = 0 ; i < keySetMax (greens) ; i++)
	  {
	    KEY green = keySet (greens, i) ;
	    bsAddKey (Gene, tag, green) ;
	    tag = _bsDown ;
	  }
	bsSave (Gene) ;
      }
    keySetDestroy (greens) ;
    keySetDestroy (greens2) ;
    keySetDestroy (reads) ;
  }
  if (!doLimit && ! searchRepeats)
    {
      doLimit = cDNAlimitRepeatedGene (cosmid, gene, a1, a2, &z1, &z2) ;
    }
  if (cosmid)
    {
      KEYSET newGenes = 0 ;
      int mySplitCloud  ;

      if (searchRepeats)
        {
          /* recompute without splitting-cloud */
          mySplitCloud = 0 ;
          newGenes = cDNARealignCutPart (cosmid, gene, z1, z2
                                         , direction, 2  /* do not move out of the requested region */
                                         , searchRepeats, 0) ;
          /* then fuse locally and shed */ 
          if (newGenes)
            { 
              KEYSET ks1 = query (newGenes, "to_be_fused_with") ;
              KEYSET ks2 = query (newGenes, "NOT to_be_fused_with") ;

              if (keySetMax (ks1))
                {
                  cDNAFuseGeneKeySet (ks1) ;
                  cDNARealignGeneKeySet (ks1, FALSE, 2, 0, 0) ;
                  if (splitCloud)
                    cDNARealignGeneKeySet (ks1, FALSE, 2, 0, 1) ;
                }
              if (splitCloud) 
                 mySplitCloud = 1 ;
              if (keySetMax (ks2))
                cDNARealignGeneKeySet (ks2, FALSE, 2, 0, mySplitCloud) ;
              keySetDestroy (ks1) ;
              keySetDestroy (ks2) ;
            }
          /* then clean up locally */
        }
       else if (splitCloud ==  2) /* force recalculation, then split */
        {
          mySplitCloud = 0 ; /* force recalculation if several mrna-sets  OR if !read */
          newGenes = cDNARealignCutPart (cosmid, gene, z1, z2
                                         , direction, doLimit
                                         , 0, mySplitCloud) ;  
          mySplitCloud = 1 ; /* recalculate if several mrna-sets */
          if (newGenes)
            cDNARealignGeneKeySet (newGenes, FALSE, 2, 0, mySplitCloud) ;
        }
       else if (splitCloud == 4) /* split geneId recalculation, then split */
        {
          mySplitCloud = 0 ; /* force recalculation if several mrna-sets  OR if !read */
          newGenes = cDNARealignCutPart (cosmid, gene, z1, z2
                                         , direction, doLimit
                                         , 0, 4) ;  
          mySplitCloud = 1 ; /* recalculate if several mrna-sets */
          if (newGenes)
            cDNARealignGeneKeySet (newGenes, FALSE, 2, 0, mySplitCloud) ;
        }
       else if (! isComposite)  /* (splitCloud ==  0, 1) AND !searchRepeats */
        {
          newGenes = cDNARealignCutPart (cosmid, gene, z1, z2
                                         , direction, doLimit
                                         , 0, splitCloud) ;  
        }
      keySetDestroy (newGenes) ;
    }
 done:
  nRecursion-- ;
  return newGene ;
}  /* cDNARealignGene */

/*********************************************************************/

static void cDNARealignGeneKeySet (KEYSET ks, BOOL doFuse, int locally, int searchRepeats, int splitCloud)
{
  KEYSET ksa = 0, ks1 = ks ? query (ks, "CLASS Transcribed_gene") :
    query (0, "FIND Transcribed_gene") ;
  int nn, i, j = 0, k = 0 ;
  char timeBuf[26] ;
  char *searchGene = 0 ; /* set to zero usually, to 1P390 to start the looping on that gene */
  time_t now = time (0), now2 ;  /* in seconds */

  if (!sessionGainWriteAccess())
    goto abort ;
  cDNAAlignInit() ;

  if (doFuse)
    cDNAFuseGeneKeySet (ks1) ;
  doFuse = FALSE ; /* since it is done */
  nn = i = keySetMax (ks1) ;
  ksa = keySetAlphaHeap (ks1, nn) ;
  while (i--)
    {
      if (!j--)
        { j = 9 ; printf (" (%d to go) \nRealign: ", i) ; }
      printf ("\n%s:%s", timeShow(timeNow(), timeBuf, 25), name(keySet(ksa, i))) ;
      if (searchGene && lexstrcmp (name(keySet(ksa,i)), searchGene))
        continue ;
      searchGene = 0 ;
      cDNARealignGene (keySet(ksa, i), 0, 0, 0, doFuse, locally, searchRepeats, splitCloud) ;
      messalloccheck() ;
      if (! ++k%30)
        { alignEst2Cosmid (0, 0, 0, 0, 9999, 0, 0, 0, 0, 0) ; }
      if (k == 200) /* check time every 200 genes */
        { 
          k = 0 ;  
          now2 = time (0) ;
            if (now2 - now > 1800) /* save every 30 minutes elapsed time */
            {
              now = now2 ;
              sessionDoSave (TRUE) ; tStatus () ;
              messout ("%d gene realigned, %d togo", nn -i, i) ;
            }
        }
    }

 abort:
  keySetDestroy (ks1) ;
  keySetDestroy (ksa) ;
  cDNAEliminateDeadMrnas () ;
} /* cDNARealignGeneKeySet */

/*********************************************************************/
#ifdef JUNK
static void cDNACreateVirtualTileOld (void)
{
  KEYSET ksEst = 0, ksGene = 0, ksClo = 0, ksTile = 0, ksSource = 0 ;
  Array aa = 0 ;
  int x1, x2, a1, a2, ii, i ;
  KEY est, vTile, source ;
  OBJ Source = 0 ;
  vTXT txt = 0 ;

  if (!sessionGainWriteAccess())
    goto abort ;

  ksEst = query (0, "Find est super_long && from_gene ; COUNT { >Genomic_sequence ; is_gene_tile } = 0; { is_read} SETOR {>from_gene;>read;NM*}") ;
  ii  = keySetMax (ksEst) ;
  while (ii--)
    {
      est = keySet (ksEst, ii) ;
      ksGene = queryKey (est, ">from_gene ; COUNT { >Genomic_sequence ; is_gene_tile } = 0") ;
      if (keySetMax (ksGene) > 1)
        {
          /* create a virtual tile by assembling all the tiles of this gene */
          ksTile = query (ksGene, ">genomic_sequence") ;
          ksSource = query (ksTile, ">Source") ;
          if (keySetMax (ksSource) != 1)
            goto abort ;
          /* locate extremities of the tiles */
          source = keySet (ksSource, 0) ;
          a1 = a2 = -1 ;
          if ((Source = bsCreate (source)))
            {
               for (i = 0 ; i < keySetMax(ksTile) ; i++)
                 if (bsFindKey (Source, _Subsequence, keySet (ksTile, i)) &&
                     bsGetData (Source, _bsRight, _Int, &x1) &&
                     bsGetData (Source, _bsRight, _Int, &x2))
                   {
                     if (a1 == -1)
                       { if (x1 < x2) { a1 = x1 ; a2 = x2 ; } else { a1 = x2 ; a2 = x1 ; } }
                     if (x1 < a1) a1 = x1 ;
                     if (x2 < a1) a1 = x2 ;
                     if (x1 > a2) a2 = x1 ;
                     if (x2 > a2) a2 = x2 ;
                   }
               else
                 goto abort ;
              bsDestroy (Source) ;
            }
          for (i = 0 ; ; i++)
            if (lexaddkey (messprintf ("%s_%d_v", name(source), i), &vTile, _VSequence))
              break ;
              
          /* destroy pieces */
          ksClo = query (ksGene, ">cdna_clone") ;
          txt = vtxtCreate () ;
          for (i = 0 ; i < keySetMax(ksGene) ; i++)
            vtxtPrintf (txt, "Transcribed_gene %s\n-D To_be_fused_with\n-D cDNA_clone\n\n", name(keySet (ksGene,i)));
          parseBuffer (vtxtPtr (txt), 0) ;
          for (i = 0 ; i < keySetMax(ksGene) ; i++)
            cDNARealignGene (keySet (ksGene,i), 0, 0, 0, 0, 0, 0, 0) ;  /*no repeats  clean up */
        
          /* create a gene on this tile in a minimal way */
          vtxtClear (txt) ;
          vtxtPrintf (txt, "Transcribed_gene %s\n", name (keySet (ksGene, 0))) ;
          for (i = 0 ; i < keySetMax(ksClo) ; i++)
            vtxtPrintf (txt, "cDNA_clone %s\n", name (keySet (ksClo, i))) ;
          vtxtPrintf (txt, "\nSequence %s\n", name (source)) ;
          vtxtPrintf (txt, "Subsequence %s %d %d\n\n", name (vTile), a1, a2) ;
          vtxtPrintf (txt, "\nSequence %s\n", name (vTile)) ;
          vtxtPrintf (txt, "Is_gene_tile\n") ;
          vtxtPrintf (txt, "Transcribed_gene %s 1 10000\n", name (keySet (ksGene, 0))) ;
          vtxtPrintf (txt, "\n") ;
          parseBuffer (vtxtPtr (txt), 0) ;
          /* recompute this gene */
          cDNARealignGene (keySet (ksGene,0), 0, 0, 0, 0, 0, 0, 0) ;  /* no repeats */
        }
    }
  sessionDoSave (TRUE) ;
 abort:
  keySetDestroy (ksEst) ;
  keySetDestroy (ksGene) ;
  keySetDestroy (ksTile) ;
  keySetDestroy (ksClo) ;
  keySetDestroy (ksSource) ;
  arrayDestroy (aa) ;
  bsDestroy (Source) ;
} /* cDNACreateVirtualTile */
#endif
/*********************************************************************/
/*********************************************************************/

static int cDNAVirtualGetAllPg (Array pgs, Array pgsBig)
{
  int TILE_SIZE = 600000 ;
  int ii, jj, nBig = 0, a1, a2, a0 ;
  HIT *up, *vp ;
  KEYSET ks = query (0, "Find predicted_gene") ;
  OBJ Pg = 0 ;
  KEY map, cosmid, pg ;

  for (ii = jj = 0 ; ii < keySetMax (ks) ; ii++)
    {
      pg = keySet (ks, ii) ;
      if ((Pg = bsCreate (pg)))
        {
          if (bsGetKey (Pg, _IntMap, &map) &&
              bsGetData (Pg, _bsRight, _Int, &a1) &&
              bsGetData (Pg, _bsRight, _Int, &a2)
              )
            {
              /* register */
              up = arrayp (pgs, jj++, HIT) ;
              if (a1 > a2)
                { up->reverse = TRUE ; a0 = a1 ; a1 = a2 ; a2 = a0 ; }
              up->gene = pg ;
              up->x1 = up->a1 = a1 ; up->x2 = up->a2 = a2 ; up->cDNA_clone = map ;
              if (bsGetKey (Pg, _Source, &cosmid) &&
                  keyFindTag (cosmid, str2tag ("Is_gene_tile"))) ;
              else 
                {
                  if (a2 - a1 >  TILE_SIZE)
                    {
                      /* do not trust the current limits */
                      up->a1 -= 100000 ; up->a2 += 100000 ;
                      vp = arrayp (pgsBig, nBig++, HIT) ;
                      *vp = *up ;
                    }
                }
            }
          bsDestroy (Pg) ;
        }
    }
  keySetDestroy (ks) ;

  return nBig ;
} /* cDNAVirtualGetAllPg  */

/*********************************************************************/
/*********************************************************************/

static int cDNAVirtualGetAllTg (Array tgs)
{
  int ii, jj, nBig = 0, a0, a1, a2 ;
  HIT *up ;
  KEYSET ks = query (0, "Find transcribed_gene") ;
  OBJ Tg = 0 ;
  KEY map, tg ;

  for (ii = jj = 0 ; ii < keySetMax (ks) ; ii++)
    {
      tg = keySet (ks, ii) ;
      if ((Tg = bsCreate (tg)))
        {
          if (bsGetKey (Tg, _IntMap, &map) &&
              bsGetData (Tg, _bsRight, _Int, &a1) &&
              bsGetData (Tg, _bsRight, _Int, &a2)
              )
            {
              /* register */
              up = arrayp (tgs, jj++, HIT) ;
              if (a1 > a2)
                { up->reverse = TRUE ; a0 = a1 ; a1 = a2 ; a2 = a0 ; }
              up->gene = tg ;
              up->x1 = up->a1 = a1 ; up->x2 = up->a2 = a2 ; up->cDNA_clone = map ;
            }
          bsDestroy (Tg) ;
        }
    }
  keySetDestroy (ks) ;

  return nBig ;
} /* cDNAVirtualGetAllTg  */

/*********************************************************************/
/* extend the zone of each big pg to include the neighbouring g */
static int cDNAVirtualMerge (Array aa, Array bb)
{
  int ii, jj, zone = 0 ;
  HIT *up, *vp ;
  BOOL ok ;

  for (ii = 0 ; ii < arrayMax (aa) ; ii++)
    {
      up = arrp (aa, ii, HIT) ;
      up->zone = 0 ;
    }

  for (ii = 0 ; ii < arrayMax (aa) ; ii++)
    {
      up = arrp (aa, ii, HIT) ;
      if (up->zone)
        continue ;
      ok = TRUE ;
      up->zone = ++zone ;
      while (ok)
        {
          ok = FALSE ; /* iterate if the box becomes bigger */
          for (jj = 0 ; jj < arrayMax (bb) ; jj++)
            {
              vp = arrp (bb, jj, HIT) ;
              if (vp->zone || 
                  vp->cDNA_clone != up->cDNA_clone ||
                  vp->reverse != up->reverse
                  )
                continue ;
              if (vp->a2 > up->a1 && vp->a1 < up->a1)
                { up->a1 = vp->a1 ; ok = TRUE ; }
              if (vp->a1 < up->a2 && vp->a2 > up->a2)
                { up->a2 = vp->a2 ; ok = TRUE ; }
              if (vp->a1 < up->a2 && vp->a2 > up->a1)        
                vp->zone = zone ; /* absorbed */
            }
          for (jj = 0 ; jj < arrayMax (aa) ; jj++)
            {
              vp = arrp (aa, jj, HIT) ;
              if (vp->zone || 
                  vp->cDNA_clone != up->cDNA_clone ||
                  vp->reverse != up->reverse
                  )
                continue ;
              if (vp->a2 > up->a1 && vp->a1 < up->a1)
                { up->a1 = vp->a1 ; ok = TRUE ; }
              if (vp->a1 < up->a2 && vp->a2 > up->a2)
                { up->a2 = vp->a2 ; ok = TRUE ; }
              if (vp->a1 < up->a2 && vp->a2 > up->a1)        
                vp->zone = -1 ; /* absorbed */
            }
        }
    }
  /* keep happy few */
 
  for (ii = jj = 0 ; ii < arrayMax (aa) ; ii++)
    {
      up = arrp (aa, ii, HIT) ;
      if (up->zone == -1)
        continue ;
      if (jj < ii)
        { vp = arrp (aa, jj, HIT) ; *vp = *up ; }
      jj++ ;
    }
  arrayMax (aa) = jj ;
 
  return jj ;  
} /* cDNAVirtualMerge */

/*********************************************************************/
/* gather the clones and desroy the tgs */
static KEYSET cDNAVirtualGetClones (HIT *up, Array tgs, KEY *tgp)
{
  KEYSET clones = keySetCreate () ;
  int ii, jj, kk ;
  HIT *vp ;
  KEYSET ks = 0 ;
  KEY tg ;
  OBJ Tg = 0 ;

  *tgp = 0 ;
  for (ii = kk = 0 ; ii < arrayMax (tgs) ; ii++)
    {
      vp = arrp (tgs, ii, HIT) ;
      if (vp->zone != up->zone)
        continue ;
      tg = vp->gene ;
      if (! *tgp) 
        *tgp = tg ;
      ks = queryKey (tg, ">cDNA_clone") ;
      for (jj = 0 ; jj < keySetMax (ks) ; jj++)
        keySet (clones, kk++) = keySet (ks, jj) ;
      keySetDestroy (ks) ;
      if ((Tg = bsUpdate (tg)))
        {
          if (bsFindTag (Tg, _cDNA_clone))
            bsRemove (Tg) ;
          bsSave (Tg) ;
          cDNARealignGene (tg, 0, 0, 0, 0, 2, 0, 0) ; /*nofuse, locally, does the clean up */
        }
    }
  keySetSort (clones) ;
  keySetCompress (clones) ;

  return clones ;
} /* cDNAVirtualGetClones */

/*********************************************************************/

static KEY cDNAVirtualCreateTile (HIT *up, KEY tg)
{
  KEY source, tile = 0 ;
  OBJ Source = 0, Tile = 0 ;
  int i, x1, x2, dx1, dx2 ;

  source = keyGetKey (up->gene, _Source) ;
  /* locate the pg in the source */
  if ((Source = bsUpdate (source)))
    {
      if (bsFindKey (Source, _Subsequence, up->gene) &&
          bsGetData (Source, _bsRight, _Int, &x1) &&
          bsGetData (Source, _bsRight, _Int, &x2)) 
        {
          dx1 = up->x1 - up->a1 ;
          dx2 = up->a2 - up->x2 ; /* extension,of both ends allways positive */
          
          if (x1 <= x2 && !up->reverse)
            { x1 -= dx1 ; x2 += dx2 ; }
          if (x1 > x2 && !up->reverse)
            { x1 += dx1 ; x2 -= dx2 ; }
          if (x1 <= x2 && up->reverse)
            { x1 -= dx2 ; x2 += dx1 ; }
          if (x1 > x2 && up->reverse)
            { x1 += dx2 ; x2 -= dx1 ; }
        }

      for (i = 0 ; ; i++)
        if (lexaddkey (messprintf ("%s_%d_v", name(source), i), &tile, _VSequence))
          break ;
      bsAddKey (Source, _Subsequence, tile) ;
      bsAddData (Source, _bsRight, _Int, &x1) ;
      bsAddData (Source, _bsRight, _Int, &x2) ; 
      
      bsSave (Source) ;
    }

  if (tile && (Tile = bsUpdate (tile)))
    {
      bsAddTag (Tile, str2tag("Is_gene_tile")) ;
      bsAddKey (Tile, _IntMap, up->cDNA_clone) ;
      if (!up->reverse)
        {
          bsAddData (Tile, _bsRight, _Int, &(up->x1)) ;
          bsAddData (Tile, _bsRight, _Int, &(up->x2)) ;
        }
      else
        {
          bsAddData (Tile, _bsRight, _Int, &(up->x2)) ;
          bsAddData (Tile, _bsRight, _Int, &(up->x1)) ; 
        }
      bsAddKey (Tile, _Transcribed_gene, tg) ;
      x1 = 1 ;
      bsAddData (Source, _bsRight, _Int, &x1) ;
      x2 = up->a2 - up->a1 + 1 ;
      bsAddData (Source, _bsRight, _Int, &x2) ; 
      
      bsSave (Tile) ;
    }
  if (1) /* debug */
    {
      Array dna = dnaGet (tile) ;
      x1 = dna ? arrayMax (dna) : 0 ;
      arrayDestroy (dna) ;
    }
              
  return tile ;
} /* cDNAVirtualCreateTile */

/*********************************************************************/

static KEY cDNAVirtualCreateTg (HIT *up, KEY cosmid, KEY tg, KEYSET clones)
{
  int ii ;
  OBJ Tg = bsUpdate (tg) ;

  bsAddKey (Tg, _IntMap, up->cDNA_clone) ;
  if (! up->reverse)
    {
      bsAddData (Tg, _bsRight, _Int, &(up->a1)) ;
      bsAddData (Tg, _bsRight, _Int, &(up->a2)) ;
    }
  else
    {
      bsAddData (Tg, _bsRight, _Int, &(up->a2)) ;
      bsAddData (Tg, _bsRight, _Int, &(up->a1)) ;
    }
  for (ii = 0 ; ii < keySetMax (clones) ; ii++)
    bsAddKey (Tg, _cDNA_clone, keySet (clones, ii)) ;
  bsAddKey (Tg, _Source, cosmid) ;
              
  bsSave (Tg) ;
  return tg ;
} /*  cDNAVirtualCreateTg */

/*********************************************************************/
/* consider all large predicted models not yet in a virtual tile
 * construct a large enough virtual tiles  
 * transfer all the relevant pg/tg to this tile
 * recompute the gene
 */

static void cDNAVirtual (void)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii,  nBig = 0 ;
  KEYSET clones = 0 ;
  Array pgs = arrayHandleCreate (10000, HIT, h) ;
  Array pgsBig = arrayHandleCreate (10000, HIT, h) ;
  Array tgs = arrayHandleCreate (10000, HIT, h) ;
  HIT *up ;
  KEY cosmid, tg ;

  /* get all large predicted models not yet in a virtual tile */
  nBig = cDNAVirtualGetAllPg (pgs, pgsBig) ;
  if (!nBig)
    goto done ;
  cDNAVirtualGetAllTg (tgs) ;

  /* extend the pgsBig to include neighbouring pgs */
  cDNAVirtualMerge (pgsBig, pgs) ;
  /* extend the pgsBig to include neighbouring tgs */
  cDNAVirtualMerge (pgsBig, tgs) ;
  /* foreach pgsBig create and assemble a merged virtual tile */
  for (ii = 0 ; ii < arrayMax (pgsBig) ; ii++)
    {
      up = arrayp (pgsBig, ii, HIT) ;
      clones = cDNAVirtualGetClones (up, tgs, &tg) ; /* gather the clones and desroy the tgs */
      if (tg)
        {
          cosmid = cDNAVirtualCreateTile (up, tg) ;
          cDNAVirtualCreateTg (up, cosmid, tg, clones) ;
          cDNARealignGene (tg, 1, up->a2 - up->a1 + 1, 0, 1, 1, 1, 1) ;
	  /* searchRepeat not rubber since we hope for some large introns */
        }
      keySetDestroy (clones) ;
    }

 done:
  messfree (h) ;
} /* cDNAVirtual */

/*********************************************************************/
/*********************************************************************/

static void alignEst2CosmidChain (KEYSET ks, int type, KEYSET estKs, BOOL doContinue)
{
  KEY *kp, cosmid1, cosmid2 ;
  int ii, jj, ng, ng1, nng, num ;
  KEYSET genes = 0 ;
  Array linkPos = arrayCreate (24, HIT) ;
  char *searchCosmid = 0 ; /* 5F60 t22_Hs22_11677_31_6  0" ; usually 0, set to "ZK637" to start looping on that cosmid */
  int ncosmid = 0 ;  /* 1 was 0,  if searchCosmid != 0, use this counter */
  int ncosmidMax = 2 ; /* 3 was 2, if searchCosmid != 0, start there, stop after ncosmidMax */
  if (!ks || !keySetMax(ks))
    return ;

  if (!sessionGainWriteAccess())
    return ;

  num = 10 ;
  ii = keySetMax(ks) ; jj = 0 ;
  kp = arrp(ks, 0, KEY) -1  ;
  ng = nng = ng1 = 0 ;

  if (type == 21)
    { /* clean up */
      KEYSET ks = query (0, "Find Sequence Is_chain OR Is_gap") ;
      keySetKill (ks) ;
      keySetDestroy (ks) ;
    }
  while (kp++, ii--) 
    {
      if (0 &&  ng > 200)   /* for debugging */
        break ;
      getSeqDna ( KEYMAKE (_VCalcul, 12345)) ;    /* cleanup */
      cosmid1 = *kp ;

      if (type == 21)
        alignChain (cosmid1, linkPos) ;
      while (cosmid1 && ncosmid < ncosmidMax)
        {
          cosmid2 = keyGetKey (cosmid1, _Overlap_right) ;
	  if (ncosmid % 10 == 9)
	    sessionDoSave (TRUE) ;	    
          switch (type)
            {
            case 21:
              if (cosmid2)
                alignLink (cosmid1, cosmid2, linkPos) ;
              break ;
            case 1: case 2: case 2001: case 2002:
              if (
                  (!doContinue || !keyFindTag (cosmid1, str2tag("Transcribed_gene"))) 
                  &&
                  (!searchCosmid || !strcasecmp(name(cosmid1), searchCosmid))
                  )
                {
                  if (searchCosmid)
                    { searchCosmid = 0 ; ncosmid = 1 ; } /* initialise countDown to limit the looping */
                  if (keyFindTag (cosmid1, str2tag("Genomic")))
                    {
                      if (type == 2)
                        genes = alignEst2Cosmid (cosmid1, estKs, genes, type/1000, 2222, 0, 0, 0, 3, 0) ;
                      else if (type != 2002)
                        genes = alignEst2Cosmid (cosmid1, 0, genes, type/1000, 2, 0, 0, 0, 3, 0) ;
                      else
                        genes = alignEst2Cosmid (cosmid1, estKs, genes, type/1000, 2, 0, 0, 0, 3, 0) ;
                    }
                }
              if (ncosmid) ncosmid++ ;
              ng = genes ? keySetMax(genes) : 0 ;
              break ;
            }
          if (!jj)  printf ("\n//") ; /* this printf must come after alignEst2Cosmid initialisation */
          printf (" %s  ",name(cosmid1)) ;
          if (!((++jj)%num)) 
            { 
              printf ("  ---> %d genes", ng - ng1) ; 
              ng1 = ng ;
              /* printf ("\n//") ;              tStatus() ; */
              printf ("\n//") ;
            }
          if (ng - nng > 20000)
            {
              nng = ng ;
              printf ("\nSaving") ;
              if (cosmid1)
              {
                OBJ Cosmid = bsCreate (cosmid1) ;
                KEY mm = 0 ; 
                float px = 0 ;
                
                if (Cosmid)
                  {
                    if (bsGetKey (Cosmid, _Map, &mm) && 
                        bsPushObj (Cosmid) &&
                        bsGetData (Cosmid, _Left, _Float, &px))
                      printf (" %s %g ", name(mm), px) ;
                    bsDestroy (Cosmid) ;
                  }
              }
              printf (" %d genes %d sequences togo ; %s\n", ng, ii, timeShowNow()) ;
              fprintf (stderr,"Saving, %d sequences togo ; %s\n", ii, timeShowNow()) ;
              sessionDoSave (TRUE) ;
              tStatus () ;
              ng = 0 ;
            }
          cosmid1 = cosmid2 ;
        } 
    }
  if (type == 1)
    printf ("\n// Constructed a total of %d genes %s\n", ng, timeShowNow()) ;
  arrayDestroy (linkPos) ;

  alignEst2Cosmid (0, 0, 0, 0, 9999, 0, 0, 0, 0, 0) ;

  if (genes)
    {
      mrnaAnalyseNeighbours (genes) ;
      keySetDestroy (genes) ;
    }
  sessionDoSave (TRUE) ;
}

/*********************************************************************/

static void alignAllEst2Cosmids (int type, KEYSET estKs)
{
  KEYSET ks = 0 ;
  Stack s = 0 ;
  char *cp ;
  BOOL doContinue = FALSE ;

  switch (type)
    {
    case 1: case 2: case 21:
    case 2001: case 2002:
      s = stackCreate (200) ;
      pushText (s, "Find sequence Genomic AND NOT Overlap_left") ;

      if ((cp = freeword()))
        {
          if (!strcmp (cp, "-continue"))
            {
              doContinue = TRUE ;
              cp = freeword () ;
            }
        }
      if (cp)
        { 
          int nchrom = 0 ;
          
          
            
          catText (s, " AND ( ") ;
          if (!strcmp(cp,"NOT")) 
            catText (s, " NOT IntMap") ;
          else do
            {
              if (nchrom++)
                catText (s,"  OR  ") ;
              catText (s, " IntMap = ") ;
              catText (s, cp) ;
            } while ((cp = freeword())) ;
          catText (s, " )") ;
        }
      ks = query (0, stackText (s,0)) ;
      if (type == 1)
        printf("// Aligning all est on all descendants of %s\n",stackText(s,0)) ;
      else if (type == 2)
        printf("// Aligning active set of %u est on all descendants of %s\n"
                , keySetMax (estKs), stackText(s,0)) ;
      else if (type == 2001)
        printf("// Just aligning all est on all descendants of %s\n",stackText(s,0)) ;
      else if (type == 2002)
        printf("// Just aligning %u active est on all descendants of %s\n", estKs ? keySetMax(estKs) : 0 , stackText(s,0)) ;
      else if (type == 21)
        printf("// Creating Link Gap on all descendants of %s\n",stackText(s,0)) ;
      alignEst2CosmidChain (ks, type, estKs, doContinue) ;
      stackDestroy (s) ;
      break ;
    }
  keySetDestroy (ks) ;
}

/*********************************************************************/
/*********************************************************************/

static Associator sageAssCreate (KEYSET sages, int nn, Stack s, Array nArray, BOOL isSage)
{
  Associator ass = 0 ;
  int i, ii = keySetMax(sages), n1 = 0, gap, dx1, dx2 ;
  unsigned int oligo, mask1 = 0, mask2 = 0, mask ;
  char *cp ;
  BOOL goodName ;

  ass = assBigCreate (2 * ii) ; 

  if (nn < 16) /* single search */
    {
      mask = mask1 = (1 << (2 * nn)) - 1 ;
      mask2 = 0 ; 
      dx1 = nn ; dx2 = 0 ; gap = 0 ;
    }
  else /* split search: 8bp--gap--8bp */
    {
      mask1 = mask2 = 0xffff ;
      mask = 0xffffffff ;
      dx1 = dx2 = 8 ;
      gap = nn - dx1 ;
    }

  if (nn > 0) 
    for (ii = 0; ii < keySetMax(sages) ; ii++)
      {
        cp = stackText (s, array(nArray, ii, int)) ;
        if (!cp)
          continue ;
        /* load the first mask */
        cp = stackText (s, array(nArray, ii, int)) ;
        oligo = 0 ; i = 0 ; goodName = TRUE ;
        
        cp-- ; i = 0 ;
        while (*++cp)
          {
            if (dx2 && i == dx1) { i = gap ; cp += gap - dx1 ; }
            if (i++ >= nn) break ;
            switch (ace_upper(*cp))
              {
              case 'A':
                oligo <<= 2 ; oligo |= B2[A_] ;
                break ;
              case 'T':
                oligo <<= 2 ; oligo |= B2[T_] ;
                break ;
              case 'G':
                oligo <<= 2 ; oligo |= B2[G_] ;
                break ;
              case 'C':
                oligo <<= 2 ; oligo |= B2[C_] ;
                break ;
              default:
                goodName = FALSE ;
                break ;
              } 
          }
        oligo &= mask ;
        if (oligo && oligo != (0xffffffff & mask) && goodName) /* oligo is zero if AAAAAAAAAAAA */
          {
            assMultipleInsert (ass, assVoid(oligo), assVoid(ii)) ;
            n1++ ;
          }
      }

  if (isSage)
    messout ("Created an associator for %d Sage %s\n", n1, timeShowNow()) ;
  else
    messout ("Created an associator for %d Primers %s\n", n1, timeShowNow()) ;
  return ass ;
}

/*********************************************************************/
/* if length > 16 we have to explicitely verify the hit */
static BOOL verifySageHit (Array dna, int anchor, int pos, Stack s, Array nArray
                           , int *nerrp, int *nNp, int errMax)
{
  char *nom = stackText(s, array(nArray, anchor, int)) ;
  int nN = 0, nerr = 0, nn = strlen(nom) ;
  char *cp = nom - 1 ;
  char cc, *cq = arrp (dna, pos, char) - 1 ;

  while (cp++, cq++, nn--)
    {
      cc = dnaEncodeChar[(int)*cp] ;
      if (! (cc & *cq))
        if (++nerr > errMax)
          return FALSE ;
      if (cc != *cq)
        nN++ ;
    }

  if (100 * nN > 10 *  strlen(nom)) /* reject over 10% ambiguous bases */
    return FALSE ;
  *nerrp = nerr ;
  *nNp = nN ;

  return TRUE ;
}

/*********************************************************************/

static void getSageHits (KEY cosmid, Associator ass, int nn, Array dna, BOOL isUp, 
                         Array hits, Stack s, Array nArray, int errMax)
{
  int i, pos, max, a1, a2, ll, gap, dx1, dx2, nerr, nN ;
  char *cp ;
  unsigned int oligo1 = 0,  oligo2 = 0, oligo, anchor, mask1 = 0, mask2 = 0, mask ;
  HIT *up ;
  const void *vp;
  Array bucket = 0 ;
  int iBucket = 0 ;


  if (nn < 16) /* single search */
    {
      mask = mask1 = (1 << (2 * nn)) - 1 ;
      mask2 = gap = 0 ; 
      dx1 = nn ; dx2 = 0 ;
    }
  else if (1) /* split search: 8bp--gap--8bp */
    {
      mask1 = mask2 = 0xffff ;
      mask = 0xffffffff ;
      dx1 = dx2 = 8 ;
      gap = nn - dx1 ;
    }
  else /* split search: gap--16bp--* */
    {
      mask1 = 0 ;
      mask = mask2 = 0xffffffff ;
      dx1 = 0 ; dx2 = 16 ;
      gap = (nn - dx2)/2 ;
    }
  i = nn - 1 ;
  oligo = 0 ;
  max = arrayMax (dna) ;
  if (max < 4*nn) return ;
  /* load the first mask */
  if (dx1)
    {
      cp = arrp(dna, 0, char) - 1  ;
      i = dx1 - 1 ;
      while (cp++, i--)
        { oligo1 <<= 2 ; oligo1 |= B2[(int)(*cp)] ; oligo1 &= mask1 ; } 
    }
  /* load the second mask */
 if (dx2)  
   {
     cp = arrp(dna, gap, char) - 1  ;
     i = dx2 - 1 ;
     while (cp++, i--)
       { oligo2 <<= 2 ; oligo2 |= B2[(int)(*cp)] ; oligo2 &= mask2 ; } 
   }
  pos = -1 ; i = arrayMax(dna) - nn + 1 ; 
  cp = arrp(dna, 0, char) + dx1 - 2  ;
  if (i > 20) while (pos++, cp++, i--)
    { 
      if (dx1)
        {
	  oligo1 <<= 2 ; oligo1 |= B2[(int)(*cp)] ; 
	  oligo1 &= mask1 ;
	}
      if (dx2)
        {
          oligo2 <<= 2 ; oligo2 |= B2[(int)(*(cp+gap))] ; oligo2 &= mask2 ;
          oligo = ((oligo1 << 16) | oligo2) & mask ;
        }
      else
        oligo = oligo1 ;

      bucket = 0 ; iBucket = 0 ;
      if (oligo && ass && assFind(ass, assVoid(oligo), &vp))
        while (assFindNext(ass, assVoid(oligo), &vp, &bucket, &iBucket))
          {
            anchor = assInt(vp) ;
            nerr = nN = 0 ;
            if ((ll = strlen(stackText(s,array(nArray, anchor, int)))) < 17 ||
                verifySageHit (dna, anchor, pos, s, nArray, &nerr, &nN, errMax))
              {
                if (isUp)
                  { a1 = max - pos ; a2 = max - pos - ll + 1 ;}
                else
                  { a1 = pos + 1 ; a2 = pos + ll ; }
                up = arrayp (hits, arrayMax(hits), HIT) ;
                up->gene = cosmid ; up->est = anchor ; up->reverse = isUp ;
                up->a1 = a1 ; up->a2 = a2 ; up->x1 = 1 ; up->x2 = nn ; 
                up->nerr = nerr ; up->nerrAll = nerr + nN ;
              }
          }
    }
}

/*********************************************************************/

static void registerSageHits (KEY cosmid, Array hits, Array sages, BOOL isSage)
{
  OBJ Cosmid = 0, Sage = 0 ;
  KEY sage = 0, sageMethod ;
  HIT *up ;
  int ii = arrayMax(hits) ;
  float zero = 0 ;
  BOOL isReverse = FALSE ;
  KEY sHit0 = isSage ? str2tag ("Sage_hit") : str2tag ("Primer_hit") ;
  KEY sHit1 = isSage ? str2tag ("Sage_quasi_hit") : str2tag ("Primer_quasi_hit") ;
  KEY _QuasiHits = str2tag ("Quasi_hits") ;

  lexaddkey (isSage ? "Sage" : "Primer", &sageMethod, _VMethod) ;
  if (bsIsTagInClass (class(cosmid), sHit0) &&
      (Cosmid = bsUpdate(cosmid)))
    {
      for (ii = 0 ; ii < arrayMax(hits) ; ii++)
        {
          up = arrayp (hits, ii, HIT) ;
          sage = keySet (sages, up->est) ;
          bsAddKey (Cosmid, up->nerr ? sHit1 : sHit0, sage) ;
          bsAddKey (Cosmid, _bsRight, sageMethod) ;
          bsAddData (Cosmid, _bsRight, _Float, &zero) ;
          bsAddData (Cosmid, _bsRight, _Int, &up->a1) ;
          bsAddData (Cosmid, _bsRight, _Int, &up->a2) ;
          bsAddData (Cosmid, _bsRight, _Int, &up->x1) ;
          bsAddData (Cosmid, _bsRight, _Int, &up->x2) ; 
          if (up->nerr)
            bsAddData (Cosmid, _bsRight, _Int, &up->nerr) ; 
        }
      bsSave (Cosmid) ;
    }
  /* add the coords in the sage */
  for (ii = 0 ; ii < arrayMax(hits) ; ii++)
    {
      up = arrayp (hits, ii, HIT) ;
      sage = keySet (sages, up->est) ;

      if ((Sage = bsUpdate(sage)))
        {       
          /* register coordinates in the sage, XREF cannot do that */
          isReverse = bsFindTag (Sage, _Reverse) ;
          if (isReverse)
            { int aa = up->a1 ; up->a1 = up->a2 ; up->a2 = aa ; }
          bsAddKey (Sage, up->nerr ? _QuasiHits : _Hits, up->gene) ;
          bsAddData (Sage, _bsRight, _Int, &up->a1) ;
          bsAddData (Sage, _bsRight, _Int, &up->a2) ; 
          if (up->nerr)
            bsAddData (Sage, _bsRight, _Int, &up->nerr) ;
          bsSave(Sage) ;
        }
    }
}

/*********************************************************************/

static void alignSage2Cosmids (KEYSET sages, KEYSET cosmids, BOOL isSage, int errMax, BOOL oneStrand)
{
  Associator ass = 0 ;
  int ii = 0, nn = 0 ;
  KEY *kp ;
  Array dna = 0, hits = 0,
    nArray = arrayCreate (1000, int) ;
  int sageLengthMin = 0, sageLengthMax = 0 ;
  Stack s = stackCreate (3000) ;
  OBJ Sage = 0 ; 
  char *cp ;

  for (ii = 0; ii < keySetMax(sages) ; ii++)
    {
      kp = arrp (sages, ii, KEY) ;
      cp = 0 ;

      if ((Sage = bsCreate (*kp)))
        {
          bsGetData (Sage, _Motif, _Text, &cp) ;
          if (cp)
            { 
              int n = strlen(cp) ;
              if (! sageLengthMax)
                sageLengthMin = sageLengthMax = n ;
              if (n > sageLengthMax)
                sageLengthMax = n ;
              if (n < sageLengthMin)
                sageLengthMin = n ;
              array (nArray, ii, int) = stackMark(s) ;
              pushText (s, cp) ;
            }
          bsDestroy (Sage) ;
        }
    }
  
  if (sageLengthMin > 7) /* 12 */
    {
      ass = sageAssCreate (sages, sageLengthMin, s, nArray, isSage) ;
      hits = arrayCreate (20000, HIT) ;
      kp = arrp (cosmids, 0, KEY) - 1 ;
      ii = keySetMax(cosmids) ;
      
      while (kp++, ii--)
        {
          dna = dnaGet (*kp) ;
          if (dna)
            {
              getSageHits (*kp, ass, sageLengthMin, dna, FALSE, hits, s, nArray, errMax) ;
              if (! oneStrand)
                {
                  reverseComplement (dna) ;
                  getSageHits (*kp, ass, sageLengthMin, dna, TRUE, hits, s, nArray, errMax) ;
                }
              if (1)
                registerSageHits (*kp, hits, sages, isSage) ;
              nn += arrayMax(hits) ;
              arrayReCreate (hits, 2000, HIT) ;
              arrayDestroy (dna) ; 
            }
        }
    }
  else
    messout ("Cancelled because the shortest sage only has %d <= 7 (was 12) bp", sageLengthMin) ;
  
  if (isSage)
    messout ("Found %d Sage Hits\n", nn) ;
  else
    messout ("Found %d Primer Hits\n", nn) ;
  arrayDestroy (hits) ;
  arrayDestroy (dna) ;
  arrayDestroy (nArray) ;
  assDestroy (ass) ;
  stackDestroy(s) ;
}

/*********************************************************************/

static void alignAllSage (KEYSET taceKs, BOOL onMrna, BOOL isSage, int errMax, BOOL oneStrand)
{
  KEYSET cosmids = 0, sages = 0 ;
  char *cp = onMrna ? "Find mRNA ; DNA" : "Find sequence ; Genomic OR Genomic_canonical" ;

  cosmids = query (0, cp) ;
  sages = query (taceKs, "motif") ;
  if (keySetMax(sages) && keySetMax(cosmids))
    alignSage2Cosmids (sages, cosmids,isSage, errMax, oneStrand) ;

  keySetDestroy (cosmids) ;
  keySetDestroy (sages) ;
}

/*********************************************************************/
/*********************************************************************/

typedef struct darStruct 
{ 
  KEY gene ;
  KEY map ; int m1, m2 ;  /* intMap coord */
  KEY pmap ; int pm1, pm2 ;  /* previous int map */
  KEY cosmid ; int c1, c2 ; /* cosmid coord */
  KEY previousCosmid ; int pc1, pc2 ; /* previous cosmid coord */
  int index ; /* index in mapOrder array */
  int nHits ; /* nb of hits of this gene in whole genome */
  BOOL samePos ;
} DAR ;

/*********************************************************************/

static int  doubleAlignRnai (KEY rnai, Array allDars)
{
  int iAllDars = arrayMax(allDars) ;
  OBJ GF = 0, Rnai = 0, O1 = 0 , O2 = 0 ;
  KEY o1 = 0, o2 = 0, dummy, s1, s2, previousGenefinder = 0, pgfmap = 0, pCosmidMap = 0 ;
  int nali = 0, pgfa1, pgfa2, delta ;
  Array h1 = 0, h2 = 0 ;
  BSunit *u1, *u2 ;
  int i1, i2, x1, x2, y1, y2, a1, a2, b1, b2,
    previousA1, previousA2, bestDelta ;
  BOOL debug = FALSE ;
  OBJ Seq1 = 0, Seq2 = 0 ;
  Array dars = 0 ;
  DAR *hh, *hh2 ;
  int ih = 0, jh, limit ;
  KEY smap1, smap2, previousCosmid, cosmid1, cosmid2 ;
  KEY _Primers = str2tag("Primers") ;
  KEYSET parts = 0 ;

  Rnai = bsUpdate (rnai) ;
  if (!Rnai)
    return 0 ;
  
  if (strncmp(name(rnai),"mv_",3) && class(rnai) == _VRNAi)
    limit = 3000 ;
  else
    limit = 60000 ;

 
  if (bsGetKey (Rnai, _Primers, &o1) &&
      bsGetKey (Rnai, _bsDown, &o2) &&
      (O1 = bsCreate (o1)) &&
      (O2 = bsCreate (o2)))
    {
      if (bsFindTag (O1, _Forward))
        {
          if (bsFindTag (O2, _Reverse)) ;
          else 
            {
              freeOutf ("// ERROR %s forward forward\n", name(rnai)) ;
              goto abort ;
            }
        }
      else
        {
          if (bsFindTag (O2, _Forward)) 
            { OBJ obj = O1 ; O1 = O2 ; O2 = obj ; dummy = o1 ; o1 = o2 ; o2 = dummy ; }
          else 
            {
              freeOutf ("// ERROR %s reverse reverse\n", name(rnai)) ;
              goto abort ;
            }
        }

      h1 = arrayCreate (32, BSunit) ;
      h2 = arrayCreate (32, BSunit) ;
      bsGetArray (O1, str2tag("Hits"), h1, 3) ;
      bsGetArray (O2, str2tag("Hits"), h2, 3) ;

      for (i1 = 0 ; i1 < arrayMax(h1) ; i1 += 3)
        { 
          u1 = arrp (h1, i1, BSunit) ;
          s1 = u1[0].k ; x1 = u1[1].i ; x2 = u1[2].k ;
          for (i2 = 0 ; i2 < arrayMax(h2) ; i2 += 3)
            {
              u2 = arrp (h2, i2, BSunit) ;
              s2 = u2[0].k ; y1 = u2[1].i ; y2 = u2[2].k ;

              if (debug)freeOutf ("// PRR %s  s1:%s x1=%d x2=%d   s2=%s  y1=%d y2=%d\n",
                                name(rnai), name(s1), x1, x2, name(s2), y1, y2) ;
              
              bsDestroy (Seq1) ;                  
              bsDestroy (Seq2) ; 
              
              if ((Seq1 = bsCreate (s1))&&
                  bsGetKey (Seq1, _IntMap, &smap1) &&
                  bsGetData (Seq1, _bsRight, _Int, &a1) &&
                  bsGetData (Seq1, _bsRight, _Int, &a2))
                {
                  if ((Seq2 = bsCreate (s2))&&
                      bsGetKey (Seq2, _IntMap, &smap2) &&
                      bsGetData (Seq2, _bsRight, _Int, &b1) &&
                      bsGetData (Seq2, _bsRight, _Int, &b2))
                    {
                      if (smap1 != smap2)
                        continue ;
                      /* put oligo2 in s1 coordinates */
                      if (a1 < a2 && b1 < b2)
                        { y1 += b1 - a1 ; y2 += b1 - a1 ; }
                      else
                        continue ; /* i am too lazy for that case */
                    }
                }
              
              /* store absolute coordinates */
              if (a1 < a2)
                { int aa = a1 ; a1 = aa + x1 - 1 ; a2 = aa + x2 - 1 ; b1 = aa + y1 - 1 ; b2 = aa + y2 - 1 ; }
              else
                { int aa = a1 ; a1 = aa - x1 + 1 ; a2 = aa - x2 + 1 ; b1 = aa - y1 + 1 ; b2 = aa - y2 + 1 ; }
              if (debug)freeOutf ("// PAIR %s  s1:%s x1=%d x2=%d   s2=%s  y1=%d y2=%d\n",
                                  name(rnai), name(s1), x1, x2, name(s2), y1, y2) ;
              if (x1 < x2)
                {
                  if (y1 < y2)
                    {
                      if (y2 - x1 > 30 && y2 - x1 < limit)
                        {
                          if (!dars) 
                            dars = arrayCreate (32, DAR) ;
                          hh = arrayp (dars, ih++, DAR) ;
                          hh->gene = rnai ;
                          hh->map = smap1 ;
                          hh->cosmid = s1 ;
                          hh->m1 = a1 ;
                          hh->m2 = b2 ;
                          hh->c1 = x1 ; hh->c2 = y2 ;
                        }
                      else
                        if (debug)freeOutf ("// ERROR %s too far\n", name(rnai)) ;
                    }
                }
              else
                {
                  if (y1 > y2)
                    {
                      if (x1 - y2 > 30 && x1 - y2 < limit)
                        {
                          if (!dars) 
                            dars = arrayCreate (32, DAR) ;
                          hh = arrayp (dars, ih++, DAR) ;
                          hh->gene = rnai ;
                          hh->map = smap1 ;
                          hh->cosmid = s1 ;
                          hh->m1 = a1 ;
                          hh->m2 = b2 ;
                          hh->c1 = x1 ; hh->c2 = y2 ;
                        }
                      else
                        if (debug)freeOutf ("// ERROR %s too far\n", name(rnai)) ;
                    }
                  /* else
                     freeOutf ("// ERROR %s parallel\n", name(rnai)) ;
                  */
                }
            }
        }
    }

  if (dars)
    {
      arraySort (dars, cDNAOrderEstByA1) ;
      arrayCompress (dars) ;

      for (nali = ih = 0 ; ih < arrayMax(dars); ih++)
        {
          hh = arrp (dars, ih, DAR) ;
          if (hh->gene)
            nali++ ;
        }

      if (bsGetKey (Rnai, str2tag("Previous_cosmid"), &previousCosmid) &&
          bsGetData (Rnai, _bsRight, _Int, &previousA1) &&
          bsGetData (Rnai, _bsRight, _Int, &previousA2)) ;
      else
        previousCosmid = 0 ;
      
      parts = queryKey (previousCosmid, ">Parts") ;
      cosmid1 = cosmid2 = 0 ;
      if (parts && keySetMax(parts) > 1)
        {
          OBJ Cosmid1 = 0 ;
          cosmid1 = keySet (parts, 0) ;
          cosmid2 = keySet (parts, 1) ;
          if (cosmid1 == keyGetKey (cosmid2, str2tag("Overlap_right")))
            SWAP(cosmid1,cosmid2,KEY) ;
          if ((Cosmid1 = bsCreate (cosmid1)))
            { int dCosmid = 0 ;
              if (bsFindKey (Cosmid1, str2tag("Overlap_right"), cosmid2) &&
                  bsGetData (Cosmid1, _bsRight, _Int, &dCosmid)) ;
              else
                { cosmid1 = cosmid2 = 0 ; }
              bsDestroy (Cosmid1) ;
            }
          previousCosmid = cosmid1 ;
        }

      if (previousCosmid)
        for (ih = 0 ; ih < arrayMax(dars); ih++)
          {
            hh = arrp (dars, ih, DAR) ;
            if (!hh->gene)
              continue ;
            hh->previousCosmid = previousCosmid ;
            hh->pc1 = previousA1 ;
            hh->pc2 = previousA2 ;
          }

      /* search for identical to previous */
      nali = 0 ;
      if (previousCosmid)
        for (ih = 0 ; ih < arrayMax(dars); ih++)
          {
            hh = arrp (dars, ih, DAR) ;
            if (!hh->gene)
              continue ;
            if (
                (  /* single cosmid alignment */
                 hh->cosmid == hh->previousCosmid &&
                 (hh->c1 > hh->pc1 - 10 &&  hh->c1 < hh->pc1 + 10) &&
                 (hh->c2 > hh->pc2 - 10 &&  hh->c2 < hh->pc2 + 10) 
                ) )
              {
                for (jh = 0 ; jh < arrayMax(dars); jh++)
                  {
                    if (ih == jh)
                      arrp (dars, jh, DAR)->samePos = TRUE ;
                    else
                      arrp (dars, jh, DAR)->gene = 0 ;
                  }
                nali = 0 ;
              }
            nali++ ;
          }


      bsGetKey (Rnai, str2tag("Previous_genefinder"), &previousGenefinder) ;
      
      if (previousGenefinder && nali > 1)
        {
          if ((GF = bsCreate (previousGenefinder)))
            {
              if (bsGetKey (GF, _IntMap, &pgfmap) &&
                  bsGetData (GF, _bsRight, _Int, &pgfa1) &&
                  bsGetData (GF, _bsRight, _Int, &pgfa2)) 
                {
                  for (ih = 0 ; nali > 1 && ih < arrayMax(dars) ; ih++)
                    {
                      hh = arrp (dars, ih, DAR) ;
                      if (!hh->gene)
                        continue ;
                      if (  /* single cosmid alignment */
                          hh->map == pgfmap &&
                          (hh->m1 > pgfa1 - 400 &&  hh->m1 < pgfa1 + 400) &&
                          (hh->m2 > pgfa2 - 400 &&  hh->m2 < pgfa2 + 400) 
                          )
                        {
                          for (jh = 0 ; jh < arrayMax(dars); jh++)
                            {
                              if (ih == jh)
                                arrp (dars, jh, DAR)->samePos = TRUE ;
                              else
                                arrp (dars, jh, DAR)->gene = 0 ;
                            }
                          nali = 1 ;
                        }
                    }
                }
              else
                previousGenefinder = 0 ; /* non usable partial data, intmap missing */
              
              bsDestroy (GF) ;
            }
        }

      pgfmap = 0 ;
      if (previousCosmid && nali > 1)
        {
          if ((GF = bsCreate (previousCosmid)))
            {
              if (bsGetKey (GF, _IntMap, &pCosmidMap) &&
                  bsGetData (GF, _bsRight, _Int, &pgfa1) &&
                  bsGetData (GF, _bsRight, _Int, &pgfa2)) 
                {
                  for (ih = 0 ; nali > 1 && ih < arrayMax(dars) ; ih++)
                    {
                      hh = arrp (dars, ih, DAR) ;
                      if (!hh->gene)
                        continue ;
                      hh->pmap = pgfmap ;
                      hh->pm1 = pgfa1 + hh->pc1 - 1 ; 
                      hh->pm2 = pgfa1 + hh->pc2 - 1 ; 
                    }
                }
              bsDestroy (GF) ;
            }
        }

      if (previousCosmid && nali > 1)
        {
          /* look for length accuracy */
          bestDelta = 100000 ;
          for (jh = -1, ih = 0, hh2 = 0 ; ih < arrayMax(dars); ih++)
            {
              hh = arrp (dars, ih, DAR) ;
              if (!hh->gene || hh->map != pCosmidMap)
                continue ;
              delta = hh->m2 - hh->m1 - (hh->pm2 - hh->pm1) ;
              if (delta < 0) delta = - delta ;
              if (delta < bestDelta)
                { hh2 = hh ; bestDelta = delta ; }
            }
          /* eliminate macth less accurate or on wrong chromosome */
          if (bestDelta < 100000)
            for (jh = -1, ih = 0 ; ih < arrayMax(dars); ih++)
              {
                hh = arrp (dars, ih, DAR) ;
                if (!hh->gene)
                  continue ;
                delta = hh->m2 - hh->m1 - (hh->pm2 - hh->pm1) ;
                if (delta < 0) delta = - delta ;
                if (delta > bestDelta + 5 || hh->map != pCosmidMap)
                  { hh->gene = 0 ; nali-- ; }
                else if (delta <= bestDelta + 5 && delta > bestDelta && 
                         hh->m1 > hh2->m1 - 10 && hh->m1 < hh2->m1 + 10 &&
                         hh->map == pCosmidMap)
                  { hh->gene = 0 ; nali-- ; }
              }
          
          /* coordinate accuracy not meaning full because chromo may change */
          if (bestDelta < 30)
            {
              bestDelta = 1000000 ;
              for (jh = -1, ih = 0 ; ih < arrayMax(dars); ih++)
                {
                  hh = arrp (dars, ih, DAR) ;
                  if (!hh->gene || hh->map != pCosmidMap)
                    continue ;
                  delta = hh->m2 - hh->pm2 ;
                  if (delta < 0) delta = - delta ;
                  if (delta < bestDelta)
                    bestDelta = delta ;
                }
              /* eliminate macth less accurate or on wrong chromosome */
              if (bestDelta < 10000)
                for (jh = -1, ih = 0 ; ih < arrayMax(dars); ih++)
                  {
                    hh = arrp (dars, ih, DAR) ;
                    if (!hh->gene)
                      continue ;
                    delta = hh->m2 - hh->pm2 ;
                    if (delta < 0) delta = - delta ;
                    if (delta > bestDelta  || hh->map != pCosmidMap)
                      { hh->gene = 0 ; nali-- ; }
                    else 
                      hh->samePos = TRUE ;
                  }
            }
        }
      
      if (nali > 1)
        {
          /* look for length accuracy in overlapping matches */
          bestDelta = 100000 ;
          for (ih = 0, hh2 = 0 ; ih < arrayMax(dars); ih++)
            {
              hh = arrp (dars, ih, DAR) ;
              if (!hh->gene)
                continue ;
              delta = hh->m2 - hh->m1 - (hh->pc2 - hh->pc1) ;
              if (delta < 0) delta = - delta ;
              if (delta < bestDelta)
                { hh2 = hh ; bestDelta = delta ; }
            }
          /* eliminate overlapping macth less accurate matches */
          if (bestDelta < 100000)
            for (ih = 0 ; ih < arrayMax(dars); ih++)
              {
                hh = arrp (dars, ih, DAR) ;
                if (hh == hh2) continue ;
                if (!hh->gene) continue ;

                delta = hh->m2 - hh->m1 - (hh->pc2 - hh->pc1) ;
                if (delta < 0) delta = - delta ;
                if (delta >= bestDelta && 
                    hh->map == hh2->map &&
                    (
                     (hh->m1 < hh2->m2 && hh->m2 > hh2->m1) ||
                     (hh->m1 > hh2->m2 && hh->m2 < hh2->m1)
                     ))
                  { hh->gene = 0 ; nali-- ; }
              }
        }

      for (ih = 0, hh2 = 0, nali = 0 ; ih < arrayMax(dars); ih++)
        {
          /* remove lame ducks */
          hh = arrp (dars, ih, DAR) ;
          if (!hh->gene)
            continue ;
          
          /* remove duplicates */
          if (hh && hh2 &&
              hh->map == hh2->map &&
              hh->m1 == hh2->m1 && hh->m2 == hh2->m2)
            continue ; 
          hh2 = hh ;
          array (allDars, iAllDars++, DAR) = *hh ;
        }
      keySetDestroy (parts) ;
    }
  
 abort:
  bsDestroy (Seq1) ;                  
  bsDestroy (Seq2) ; 

  bsDestroy (O1) ;
  bsDestroy (O2) ;
  bsSave (Rnai) ;
  arrayDestroy (h1) ;
  arrayDestroy (h2) ;

  arrayDestroy (dars) ;
  return 0 ;
}

/*********************************************************************/

static int doubleReAlignMapOrder (const void *a, const void *b)
{
  const DAR *up = (const DAR *)a, *vp = (const DAR *)b ;
  
  if (up->map != vp->map) 
    return up->map - vp->map ; 
  if (up->m1 != vp->m1)  /* IntMap:2 */
    return up->m1 - vp->m1 ; 
  if (up->m2 != vp->m2)  /* IntMap:3 */
    return up->m2 - vp->m2 ; 

  return up->gene - vp->gene ; 
}

/*********************************************************************/

static int doubleReAlignGeneOrder (const void *a, const void *b)
{
  const DAR *up = (const DAR *)a, *vp = (const DAR *)b ;
  
  if (up->gene != vp->gene) 
    return up->gene - vp->gene ; 
  if (up->map != vp->map) 
    return up->map - vp->map ; 
  if (up->m1 != vp->m1)  /* IntMap:2 */
    return up->m1 - vp->m1 ; 
 
  return up->m2 - vp->m2 ; 
}

/*********************************************************************/

static void doubleReAlignRnai (Array dars, int *n2p, int *nzp)
{
  DAR *up, *vp, *wp, *zp, *hp ;
  KEY gene ;
  int n, j1, j2, j3, j4, nHit, nBetween ;
  char *cName, gName[1000] ; 
  Array dars2 = 0 ;
  Stack s = stackCreate(10000) ;

  /* register happy few */
  if (arrayMax(dars))
    {
      for (j1 = 0, j2 = 0, up = arrp(dars, j1, DAR) ; j1 < arrayMax(dars) ; j1++, up++)
        { 
          if (up->gene)
            {
              if (j1 != j2)
                arr (dars, j2, DAR) = arr(dars, j1, DAR) ;  
              j2++ ;
            }
        }
      arrayMax (dars) = j2 ;
    }

  /* count ambiguous guys */
  for (n = j1 = 0, gene = 0, up = arrp (dars, 0, DAR) ; j1 < arrayMax(dars); j1++, up++)
    {
      if (!up->gene)
        continue ;
      if (gene != up->gene)
        {
          if (n)
            for (j2 = j1 - 1, vp = up - 1 ; j2 >= 0 ; j2--, vp--)
              if (vp->gene == gene) vp->nHits = n ;
          gene = up->gene ; n = 0 ;
        }
      n++ ;
    }

  if (n)
    for (j2 = j1 - 1, vp = up - 1 ; j2 >= 0 ; j2--, vp--)
      if (vp->gene == gene) vp->nHits = n ;
  gene = up->gene ;


  /* publish non amiguous guys */
  for (n = nHit = j1 = 0, gene = 0 ; j1 < arrayMax(dars); j1++, up++)
    {
      up = arrp (dars, j1, DAR) ; 
      if (!up->gene)
        continue ;
      if (up->nHits < 1)
        continue ; 
      if (gene != up->gene) 
        {
          nHit = 0 ;
          (*n2p)++ ; /* nb aligned */ 
          if (up->nHits > 1) 
            (*nzp)++ ;
        }
      nHit++ ;
      gene = up->gene ;
      if (up->nHits > 1)
        continue ;

      cName = class(up->gene) == _VGene ? "Genes" : className(up->gene) ;
      if (nHit > 1)
        sprintf(gName, "%s_%d", name(up->gene), nHit) ;
      else
        sprintf(gName, "%s", name(up->gene)) ;

      catText (s, messprintf 
               ("Sequence \"%s\"\n%s \"%s\" %d %d\n\n",
                name (up->cosmid), cName, gName, up->c1, up->c2)) ; 

      if (up->samePos)
        catText (s, messprintf ("%s %s\nSame_position\n\n",
                                className(up->gene), name(up->gene)));
      catText (s, messprintf ("%s %s\nIntMap %s %d %d\n\n",
                              className(up->gene), gName, name(up->map), up->m1, up->m2)) ;
      if (up->nHits > 1)
        catText (s, messprintf ("%s %s\nAmbiguous\n\n",
                                className(up->gene), gName)) ;
    }
  /* store non amibuous guys in the data base */
  parseBuffer (stackText (s,0), 0) ;
  stackClear (s) ;

  /* deal with amiguous guys */
  arraySort (dars, doubleReAlignMapOrder) ;  /* in the chromosome order */
  for (j1 = 0, up = arrp (dars, 0, DAR) ; j1 < arrayMax(dars); j1++, up++)
    up->index = j1 ; /* link back to self */
  dars2 = arrayCopy (dars) ;
  arraySort (dars, doubleReAlignGeneOrder) ;  /* in gene order */

  for (nHit = j1 = 0, gene = 0 ; j1 < arrayMax(dars); j1++, up++)
    {
      up = arrp (dars, j1, DAR) ;
      if (up->nHits < 2)  /* non amiguous, already treated */
        continue ; 
      if (!up->gene)
        continue ;
      /* locate in chromo order the non abiguous genes above (wp) and below (zp) */
      j2 =  up->index ;
      vp = arrp (dars2, j2, DAR) ;
      if (!vp->gene)
        continue ;
      if (vp->gene != up->gene)
        {
          messerror ("confusion in  doubleReAlignRnai at %s, sorry", name (up->gene)) ;
          continue ;
        }
      for (j3 = j2, wp = vp ; j3 >= 0 ; wp--, j3--)
        if (wp->nHits < 2)
          break ;
      if (j3 < 0) wp = 0 ;
      for (j4 = j2, zp = vp ; j4 < arrayMax(dars2) ; zp--, j4--)
        if (zp->nHits < 2)
          break ;
      if (j4 >= arrayMax(dars2)) zp = 0 ;
      
      if (wp && zp && wp->map == zp->map)
        {
          /* search if at least one hit has same map */
          /* if yes kill others */
          nBetween = 0 ;
          for (vp = up, j2 = j1 ; j2 >= 0 ; vp--, j2--)
            if (vp->gene != gene) break ;
          vp++ ; j2++ ;
          for (j3 = j2, hp = vp ; j3 < arrayMax(dars); j3++, hp++)
            {
              if (!hp->gene) continue ;
              if (hp->gene != up->gene) break ;
              if (hp->map == wp->map &&
                  hp->pm2 - wp->pm2 >= 0 &&
                  hp->pm2 - zp->pm2 <= 0)
                {
                  nBetween++ ;
                }
            }
          if (nBetween)
            for (j3 = j2, hp = vp ; j3 < arrayMax(dars); j3++, hp++)
              {
                if (!hp->gene) continue ;
                if (hp->gene != up->gene) break ;
                if (hp->map == wp->map &&
                    hp->pm2 - wp->pm2 >= 0 &&
                    hp->pm2 - zp->pm2 <= 0) 
                  hp->nHits = nBetween ;
                else
                  hp->gene = 0 ;
              }
        }
      if (gene != up->gene) 
        {
          nHit = 0 ;
        }

      nHit++ ;
      gene = up->gene ;
      cName = class(up->gene) == _VGene ? "Genes" : className(up->gene) ;
      if (nHit > 1)
        sprintf(gName, "%s_%d", name(up->gene), nHit) ;
      else
        sprintf(gName, "%s", name(up->gene)) ;

      catText (s, messprintf 
               ("Sequence \"%s\"\n%s \"%s\" %d %d\n\n",
                name (up->cosmid), cName, gName, up->c1, up->c2)) ; 

      if (up->samePos)
        catText (s, messprintf ("%s %s\nSame_position\n\n",
                                className(up->gene), name(up->gene)));
      catText (s, messprintf ("%s %s\nIntMap %s %d %d\n\n",
                              className(up->gene), gName, name(up->map), up->m1, up->m2)) ;
      if (up->nHits > 1)
        catText (s, messprintf ("%s %s\nAmbiguous\n\n",
                                className(up->gene), gName)) ;
    }
  /* store non amibuous guys in the data base */
  parseBuffer (stackText (s,0), 0) ;
  stackClear (s) ;
  
  stackDestroy (s) ;
  arrayDestroy (dars2) ;
}

/*********************************************************************/

static void doubleAlignRnaiSet (KEYSET ks)
{
  KEY key ;
  int i, n, n1, n2, nz ;
  Array dars = arrayCreate (3 * keySetMax(ks), DAR) ;

  n1 = n2 = nz = 0 ;
  for (i = 0 ; i < keySetMax(ks) ; i++)
    {
      key = keySet (ks, i) ;
      if (!keyFindTag (key, str2tag("Primers")))
        continue ;
      n1++ ;
      n = doubleAlignRnai (key, dars) ;
    }
   
  doubleReAlignRnai (dars, &n2, &nz) ;
  freeOutf ("// Found %d Rnai, Aligned %d/%d, found %d ambiguities doubleAlign=%d\n", n1, n2, n1, nz, n) ;
  arrayDestroy (dars) ;
}

/*********************************************************************/
/*********************************************************************/

static void alignKeysetRead2KeysetRef (KEYSET ks1, KEYSET ks2, int type, char *nom)
{
  KEY *kp ;
  int ii, ng = 0, nng = 0 ;
  KEYSET genes = 0 ;

  if (!ks1 || !keySetMax(ks1) || !ks2 || !keySetMax(ks2))
    return ;

  if (!sessionGainWriteAccess())
    return ;

  genes = alignEst2Cosmid (0, 0, 0, 0, 9999, 0, 0, 0, 0, 0) ; /* clean up */
  keySetDestroy (genes) ;
  ii = keySetMax(ks2) ;
  kp = arrp(ks2, 0, KEY) -1  ; 
  printf ("\n//") ;
  while (kp++, ii--) 
    {
      printf (" %s  ",name(*kp)) ;
      if (!ii%5) printf ("\n// ") ;
      genes = alignEst2Cosmid (*kp, ks1, genes, 1, type, 0, 0, nom, 3, 0) ;
      if (genes) 
        ng = keySetMax(genes) ;
      if (ng - nng > 1000)
        {
          nng = ng ;
          printf ("\nSaving, %d sequences togo ; %s\n", ii, timeShowNow()) ;
          fprintf (stderr,"Saving, %d sequences togo ; %s\n", ii, timeShowNow()) ;
          sessionDoSave (TRUE) ;
          ng = 0 ;
        }
    } 
  printf ("\n// ") ;    
  if (genes)
    {
      mrnaAnalyseNeighbours (genes) ;
      keySetDestroy (genes) ;
    }

  genes = alignEst2Cosmid (0, 0, 0, 0, 9999, 0, 0, 0, 0, 0) ; /* clean up */ 
  keySetDestroy (genes) ;
  sessionDoSave (TRUE) ;
} /* alignKeysetRead2KeysetRef  */

/*********************************************************************/

static void alignReads2Reference (int type, KEYSET ks01, KEYSET ks02, char *nom)
{
  KEYSET ks1 = 0, ks2 = 0 ;

  ks1 = ks01 ? ks01 : query (0, "Find Sequence Is_read") ; 
  ks2 = ks02 ? ks02 : query (0, "Find sequence Is_Reference && !Assembled_from") ; 

  alignKeysetRead2KeysetRef (ks1, ks2, type, nom) ;
  if (!ks01) keySetDestroy (ks1) ;
  if (!ks02) keySetDestroy (ks2) ;
}

/*********************************************************************/
/*********************************************************************/

static int cDNAOrderToAssigncDNA (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;

  if ((up->gene || vp->gene) && ! (up->gene && vp->gene))  /* better to have a gene */
    return up->gene ? -1 : 1 ;
  if (up->x1 != vp->x1)  /* better to have 2 reads */
     return up->x1 > vp->x1 ? -1 : 1 ;
  if (up->nerr != vp->nerr)  /* better to have few errors */
    return up->nerr < vp->nerr ? -1 : 1 ;
  if (up->a1 != vp->a1)  /* better to have long match */
    return up->a1 > vp->a1 ? -1 : 1 ;
  return 
    up->gene - vp->gene ; /* not to be ambiguous */
}

/*********************************************************************/

static void assignDoRemoveOne (OBJ Gene, Array hitsR, int i0, int i1)
{
  int ii ;
  HIT *hh ;
  
  for (ii = i0; ii < i1 ; ii++)
    {
      hh = arrp (hitsR, ii, HIT) ;
      if (bsFindKey (Gene, _cDNA_clone, hh->cDNA_clone))
        bsRemove (Gene) ;
    }
}  /* assignDoRemoveOne */ 

/*********************************************************************/

static void assignDoRemoveTwo (Array hitsR, int i0, int i1)
{
  KEYSET ks = 0 ;
  KEY gene ;
  OBJ Clone = 0, Est = 0 ;
  int ii, i ;
  HIT *hh ;

  for (ii = i0; ii < i1 ; ii++)
    {
      hh = arrp (hitsR, ii, HIT) ;
      gene = hh->gene ;
      /* may be the read and clone were not updated */
      if ((Clone = bsUpdate (hh->cDNA_clone)))
        {
          if (bsFindKey (Clone, _From_gene, gene))
            bsRemove (Clone) ;
          bsSave (Clone) ;
        }
      ks = queryKey (hh->cDNA_clone, ">Read ") ;
      for (i = 0 ; i < keySetMax (ks) ; i++)
        {
          if ((Est = bsUpdate (keySet (ks, i))))
            {
              if (bsFindKey (Est, _From_gene, gene))
                bsRemove (Est) ;
              if (bsFindTag (Est, _Intron_boundaries)) /* no longer trustable */
                bsRemove (Est) ;
              bsSave (Est) ;
            }
        }
      keySetDestroy (ks) ;
    }
}  /* assignDoRemoveTwo */ 

/*********************************************************************/

static KEYSET assignDoRemove (Array hitsR)
{
  OBJ Gene = 0 ;
  KEYSET genes = keySetCreate () ;
  int ii, jj = 0, ng = 0 ;
  KEY gene = 0 ;
  HIT *h1 ;

  arraySort (hitsR, cDNAOrderByGene) ;
  arrayCompress (hitsR) ;

  for (ii = 0 ; ii <= arrayMax(hitsR)  ; ii++)
    { 
      h1 =  ii < arrayMax(hitsR) ? arrayp (hitsR, ii, HIT) : 0 ;
      if (!h1 || gene != h1->gene)
        {  
          if (gene && ii > jj)
            {
              Gene = bsUpdate (gene) ;
              assignDoRemoveOne (Gene, hitsR, jj, ii) ;
              bsSave (Gene) ;
              assignDoRemoveTwo (hitsR, jj, ii) ;
              gene = 0 ; 
            }
          if (h1)
            keySet (genes, ng++) = gene = h1->gene ;
          jj = ii ;
        }
   }

  keySetSort (genes) ;
  keySetCompress (genes) ;
  return genes ;
} /* assignDoRemove */

/*********************************************************************/

static int assignPreRemove (Array hitsR, HIT *hh, KEY cDNA, BOOL debug, char *title)
{
  KEY gene = hh->gene ;
  int jj = arrayMax (hitsR) ;
  HIT *h1 ;

  if (debug)
    printf("Assign %s removed from %s %s\n", name(cDNA), name(gene), title) ;
  h1 = arrayp (hitsR, jj++, HIT) ;
  h1->gene = gene ; h1->cDNA_clone = cDNA ;

  hh->gene = 0 ;
  return 1 ;
} /* assignPreRemove */

/*********************************************************************/

static KEYSET assignOncecDNAKeySet (KEYSET originalKeySet)
{
  int i, ii, jj, icdna, iread, ne, neg, nkept, nrecomputed, ndestroyedbyrecomputing, nHits ;
  KEY cDNA, gene, newGene, est = 0 ;
  OBJ Clone = 0, Est = 0 , Gene = 0 ; 
  HIT *hh, *h2 ;
  BSunit *u ;
  BOOL debug = FALSE ;
  BOOL isTracking = FALSE, discard11 = TRUE ;
  Array aa = 0, aa2 = 0, hits = 0, hitsR = 0 ;
  KEYSET cdnaKs = 0, reads = 0, genes = 0 ;
  KEYSET rReads = keySetCreate () ;
  int deltaQuality = 1 ;
  KEY _Mosaic = str2tag("Mosaic") ;
  KEY _Is_Mrna = str2tag("Is_Mrna") ; 
  KEY _Possible_mosaic = str2tag("Possible_mosaic") ;

  chrono ("assigncDNA") ;

  cdnaKs = originalKeySet ? query(originalKeySet, "CLASS cDNA_Clone") : query (0, "FIND  cDNA_clone") ;
  if (!keySetMax(cdnaKs))
    { messout ("No cDNA clone in list, sorry") ; goto abort ; }

  ne = neg = nrecomputed = ndestroyedbyrecomputing = 0 ;


  chrono ("assign getscores") ;
  hitsR = arrayCreate (1000, HIT) ;
  for (icdna = 0 ; icdna < keySetMax(cdnaKs) ; icdna++)
    {
      int nc2g = 0 ;
      cDNA = keySet (cdnaKs, icdna) ;

      isTracking = keyFindTag (cDNA, str2tag("Tracking_clone")) ;
      {
        /* no point loosing time on clones assigned to single gene */
        KEYSET c2g = queryKey (cDNA, "FOLLOW From_gene") ;
        nc2g = keySetMax(c2g) ;
        keySetDestroy (c2g) ;
        if (nc2g < (discard11 ? 1 : 2))   /* 1 to be able to discard quality 11 hits */
          continue ;
        else
          {
            if (debug)
              printf("Assign %s \n", name(cDNA)) ;
          }
      }
      reads = queryKey (cDNA, "FOLLOW Read") ;
      hits = arrayReCreate (hits, 64, HIT) ;
      /* find all genes hit by all reads of this clone */
      for (iread = 0 ; iread < keySetMax(reads) ; iread++)
        {
          est = keySet(reads, iread) ;
          aa = arrayReCreate (aa, 12, BSunit) ;
          if ((Est = bsCreate (est)))
            {
              BOOL reverse = bsFindTag (Est, _Reverse) ;
              if (bsGetArray(Est, _From_gene, aa, 9))  /* gene, nmatch, nerr */
                for (i = 0 ; i < arrayMax(aa) ; i += 9)
                  {
                    hh = arrayp(hits, arrayMax(hits), HIT) ;
                    u = arrp (aa, i, BSunit) ;
                    hh->gene = u[0].k ;  hh->est = est ; 
                    hh->x2 = u[8].i ; /* len to be aligned = nb of bp in the est between clips */
                    hh->a1 = u[3].i ;  /* nb of bp aligned */
                    if (hh->a1 > hh->x2) hh->a1 = hh->x2 ; /* this occurs if many indels */
                    hh->nerr = u[4].i ;  /* nb of errors */
                    if (nc2g > 1 && keyFindTag (hh->gene, _Other)) /* penalty fro blue in tg */
                      {
                        OBJ Tg ;
                        int j ;
                        BSunit *u2 ;

                        aa2 = arrayReCreate (aa2, 12, BSunit) ;
                        if ((Tg = bsCreate (hh->gene)))
                          {
                            if (bsGetArray(Est, _Other, aa2, 5))  /* gene, nmatch, nerr */
                              for (j = 0 ; j < arrayMax(aa2) ; j += 5)
                                { 
                                  u2 = arrp (aa2, j, BSunit) ;
                                  if (u2[4].k == cDNA)
                                    hh->nerr +=3 ; 
                                }
                            bsDestroy (Tg) ;
                          }
                      }
                    hh->reverse = reverse ;  /* strand */
                    if (reverse)
                      {
                        hh->ex2 = u[6].i ; /* offset = nb of bp between cliptop and first aligne dbase */
                        
                        hh->clipTop = u[5].i ; /* quality */
                        hh->clipEnd = 10000 ;
                      }
                    else
                      {
                        hh->ex1 = u[2].i ; /* offset = nb of bp between cliptop and first aligne dbase */

                        hh->clipTop = 10000 ;
                        hh->clipEnd = u[5].i ; /* quality */
                      }
                  }
              bsDestroy (Est) ;
            }
        }
      keySetDestroy (reads) ;

      /* do NOT   cDNASwapA (hits) ; a1 and a2 mean align length offset */
      arraySort (hits,cDNAOrderToAssigncDNA) ;  /* relies on clone = 0 */

      /* accumulate per gene the number of hitting reads and total match */
      for (ii = 0 ; ii < arrayMax(hits) ; ii++)
        {
          hh = arrp (hits, ii, HIT) ;
          if (!hh->gene)
            continue ;
          if (hh->a1 > 20)
            hh->x1 = hh->reverse ? 2 : 1 ;
          else
            hh->x1 = 0 ;
          if (hh->a1 == 0)
            printf("Assign gene %s clone %s a1 missing\n", name(hh->gene), name(hh->est)) ;
          for (jj = ii + 1; jj < arrayMax(hits) ; jj++)
            {
              h2 = arrp (hits, jj, HIT) ;
              if (hh->gene == h2->gene)
                { 
                  if (hh->est != h2->est && h2->a1 > 20) 
                    hh->x1 |= h2->reverse ? 2 : 1 ;      /* 1 forward, 2 reverse */
                  hh->a1 += h2->a1 ; hh->nerr += h2->nerr ;
                  hh->clipTop = hh->clipTop < h2->clipTop ? hh->clipTop : h2->clipTop ;
                  hh->clipEnd = hh->clipEnd < h2->clipEnd ? hh->clipEnd : h2->clipEnd ;
                  if (h2->ex1) hh->ex1 = h2->ex1 ;
                  if (h2->ex2) hh->ex2 = h2->ex2 ;
                  h2->gene = 0 ; 
                }
            }
        }
      if (isTracking)
        {  /*  we try to match gene: G_t1_Hs1_22215_21_1_711     to intmap:c_Hs1_22215_21 */
          int nSame = 0 ;
          KEY intMap = est ? keyGetKey (est, _IntMap) : 0 ;
          
          if (intMap && strlen(name(intMap)) > 3)
            {
              int i, j=0, nn_ = 0 ;
              char *uu = name(intMap) + 3, buf[1000] ;

              buf[j++]='*' ;
              for (i = 1, nn_ = 0 ; i < 800 ; i++)
                {
                  if (*uu == '_') nn_++ ;
                  if (nn_ >= 2) break ;
                  if (!*uu) break ;
                  buf[j++] = *uu++ ;
                }

              buf[j++]='*' ;
              buf[j++]=0;

              for (ii = 0 ; ii < arrayMax(hits) ; ii++)
                {
                  hh = arrp (hits, ii, HIT) ;
                  if (!hh->gene)
                    continue ;
                  if (pickMatch (name(hh->gene), buf))
                    nSame++ ;
                } 
              if (nSame)
                for (ii = 0 ; ii < arrayMax(hits) ; ii++)
                  {
                    hh = arrp (hits, ii, HIT) ; 
                    if (!hh->gene)
                      continue ;
                    if (nSame && pickMatch (name(hh->gene), buf))
                      { hh->x1 = 3 ; nSame = 0 ; } /* keep only one */
                    else
                      hh->x1 = 1 ; 
                  }
            }
        }
            
      /* again, because now we know x1 = 2 if 3' and 5' do hit */
      arraySort (hits,cDNAOrderToAssigncDNA) ;  /* relies on clone = 0 */
      
      if (discard11)
        for (nkept = 0, ii = 0 ; ii < arrayMax(hits) ; ii++)
          {
            hh = arrp (hits, ii, HIT) ;
            if (hh->gene &&  /* mieg jan 16, 2002 */
                ((!hh->clipTop || hh->clipTop >= 11) && (!hh->clipEnd || hh->clipEnd >= 11)))
              ne += assignPreRemove (hitsR, hh, cDNA, debug, "q >= 11") ;
          }

      for (nHits = 0, ii = 0 ; ii < arrayMax(hits) ; ii++)
        {
          hh = arrp (hits, ii, HIT) ;
          if (hh->gene)
            nHits++ ;
        }

      for (nkept = 0, ii = 0 ; ii < arrayMax(hits) ; ii++)
        {
          hh = arrp (hits, ii, HIT) ;
          if (!hh->gene)
            continue ;

          for (jj = ii + 1 ; hh->gene && jj < arrayMax(hits) ; jj++)
            {
              int discard = 0 ;
              h2 = arrp (hits, jj, HIT) ;
              if (!h2->gene)
                continue ;
              else if (hh->x1 == 3 &&  h2->x1 < 3 &&  /* beware of mosaic */
		       ( /* no real overlap */
			(h2->x1 == 1 && h2->a1 + h2->ex1 < hh->ex1 + 20) ||
			(h2->x1 == 2 && h2->a1 + h2->ex2 < hh->ex2 + 20)
			) && 
		       100 * h2->nerr < h2->a1)
		{
		  if (!keyFindTag (cDNA, _Mosaic) &&
		      keyFindTag (hh->est, _Is_Mrna) &&
		      (Clone = bsUpdate (cDNA)))
		    {
		      bsAddTag (Clone,  _Possible_mosaic) ;
		      bsSave (Clone) ;
		    } 
		}
	      else if (h2->x1 == 3 &&  hh->x1 < 3 &&  /* not both reads match gene 2 */
		       ( /* no real overlap */
			(hh->x1 == 1 && hh->a1 + hh->ex1 < h2->ex1 + 20) ||
                        (hh->x1 == 2 && hh->a1 + hh->ex2 < h2->ex2 + 20) 
			)  &&
		       100 * hh->nerr < hh->a1)  /* beware of mosaic */
		{
		  if (!keyFindTag (cDNA, _Mosaic) &&
		      keyFindTag (hh->est, _Is_Mrna) &&
		      (Clone = bsUpdate (cDNA)))
		    {
		      bsAddTag (Clone,  _Possible_mosaic) ;
		      bsSave (Clone) ;
		    }
		}

	      else if (hh->x1 == 3 &&  h2->x1 == 2 &&  /* not both reads match gene 2 */
		       h2->clipTop >= hh->clipTop - 1
		       )
		discard = 1 ;
	      else if (hh->x1 == 3 &&  h2->x1 == 2 &&  /* not both reads match gene 2 */
		       h2->clipTop < hh->clipTop - 1 && /* single hit way better */
		       h2->clipTop < hh->clipEnd - 1
		       )
		discard = -1 ;
	      else if (hh->x1 == 3 &&  h2->x1 == 1 &&  /* not both reads match gene 2 */
		       h2->clipEnd >= hh->clipEnd - 1
		       )
		discard = 1 ;
	      else if (hh->x1 == 3 &&  h2->x1 == 1 &&  /* not both reads match gene 2 */
		       h2->clipEnd < hh->clipTop - 1 && /* single hit way better */
		       h2->clipEnd < hh->clipEnd - 1
		       )
		discard = -1 ;
	      
	      else if (h2->x1 == 3 &&  hh->x1 == 2 &&  /* not both reads match gene 2 */
		       hh->clipTop >= h2->clipTop - 1
		       )
		discard = -1 ;
	      else if (h2->x1 == 3 &&  hh->x1 == 2 &&  /* not both reads match gene 2 */
		       hh->clipTop < h2->clipTop - 1 && /* single hit way better */
		       hh->clipTop < h2->clipEnd - 1
		       )
		discard = 1 ;
	      else if (h2->x1 == 3 &&  hh->x1 == 1 &&  /* not both reads match gene 2 */
		       hh->clipEnd >= h2->clipEnd - 1
		       )
		discard = -1 ;
	      else if (h2->x1 == 3 &&  hh->x1 == 1 &&  /* not both reads match gene 2 */
		       hh->clipEnd < h2->clipTop - 1 && /* single hit way better */
		       hh->clipEnd < h2->clipEnd - 1
		       )
		discard = 1 ;


              else if (hh->x1 != h2->x1)  /* do not discard 5' in favor of 3' */
                continue ;
              else if (hh->x1 != 3 &&
                       hh->x1 == h2->x1 &&  /* beware of mosaics */
                       ( /* no real overlap */
                        (hh->x1 == 1 && hh->a1 + hh->ex1 < h2->ex1 + 20) ||
                        (hh->x1 == 2 && hh->a1 + hh->ex2 < h2->ex2 + 20)
                        ) && 
                       30 * hh->nerr < hh->a1 && /* both bits of less than 3% error */
                       30 * h2->nerr < h2->a1)
                { 
                  if (!keyFindTag (cDNA, _Mosaic) &&
                      keyFindTag (hh->est, _Is_Mrna) &&
                      (Clone = bsUpdate (cDNA)))
                    {
                      bsAddTag (Clone,  _Possible_mosaic) ;
                      bsSave (Clone) ;
                    }
                  continue ;
                }
              else if (hh->x1 != 3 &&
                       hh->x1 == h2->x1 &&  /* beware of mosaics */
                       ( /* no real overlap */
                        (h2->x1 == 1 && h2->a1 + h2->ex1 < hh->ex1 + 20) ||
                        (h2->x1 == 2 && h2->a1 + h2->ex2 < hh->ex2 + 20) 
                        ) && 
                       100 * hh->nerr < hh->a1 && /* both bits of less than 1% error */
                       100 * h2->nerr < h2->a1)
                { 
                  if (!keyFindTag (cDNA, _Mosaic) &&
                      keyFindTag (hh->est, _Is_Mrna) &&
                      (Clone = bsUpdate (cDNA)))
                    {
                      bsAddTag (Clone,  _Possible_mosaic) ;
                      bsSave (Clone) ;
                    }
                  continue ;
                }
              else if (h2->clipTop >= hh->clipTop && 
                       h2->clipEnd >= hh->clipEnd &&
                       h2->clipTop +  h2->clipEnd >= hh->clipTop + hh->clipEnd + deltaQuality)
                discard = 1 ; 
              else if (hh->clipTop >= h2->clipTop && 
                       hh->clipEnd >= h2->clipEnd &&
                       hh->clipTop +  hh->clipEnd >= h2->clipTop + h2->clipEnd + deltaQuality)
                discard = -1 ; 
              
              else if (hh->a1 - hh->nerr > h2->a1 - h2->nerr + 5 &&
                       hh->clipTop <= h2->clipTop)  /* very strict filter introduced feb 1 2003 */
                { discard = 1 ; goto laba ; }
              else if (h2->a1 - h2->nerr > hh->a1 - hh->nerr + 5 &&
                       h2->clipTop <= hh->clipTop)  /* very strict filter introduced feb 1 2003 */
                { discard = -1 ; goto laba ; }
              
              else if (hh->a1 - hh->nerr > h2->a1 - h2->nerr + 2 && nHits > 4 &&
                       hh->clipTop <= h2->clipTop)  /* extra strict filter introduced feb 17 2007 */
                { discard = 1 ; goto laba ; }
              else if (h2->a1 - h2->nerr > hh->a1 - hh->nerr + 2 && nHits > 4 &&
                       h2->clipTop <= hh->clipTop)  /* extra strict filter introduced feb 17 2007 */
                { discard = -1 ; goto laba ; }
              
              else if (hh->a1 > h2->a1 + 200)  /* the new match is very much shorter */
                { discard = 1 ; goto laba ; }
              else if (hh->a1 > h2->a1 + 100)  /* the new match is much shorter */
                { 
                  if (9 * hh->nerr * h2->a1 < 10 * h2->nerr * hh->a1) /* and new error ratio not relatively worse */
                    { discard = 1 ; goto laba ; }
                  if (hh->nerr < h2->nerr + 5 || 9 * hh->nerr < 10 * h2->nerr) /* or not obsolutelly worse */
                    { discard = 1 ; goto laba ; }
                  if (2 * h2->a1 > hh->a1 &&
                      h2->nerr * hh->a1 < 4 * hh->nerr * h2->a1)
                    { discard = -1 ; goto laba ; }
                  if (hh->a1 > 3 * h2->a1 && hh->a1 > 15 * hh->nerr) /* the new much shorter, old better than 6% errors */
                    { discard = 1 ; goto laba ; }
                }
              else /* if (hh->a1 > 50 && hh->a1 > h2->a1 - 20)   the new match about same */
                {
                  if (hh->nerr < 5)
                    {  
                      if (hh->nerr + 5 < h2->nerr)  { discard = 1 ; goto laba ; }
                      if ( 3 * hh->nerr * h2->a1 < h2->nerr * hh->a1)  { discard = 1 ; goto laba ; }
                    }
                  else if (hh->nerr < 10)        
                    {  
                      if (hh->nerr + 10 < h2->nerr)  { discard = 1 ; goto laba ; }
                      if ( 2 * hh->nerr * h2->a1 < h2->nerr * hh->a1)  { discard = 1 ; goto laba ; }
                    }
                  else
                    {  
                      if (hh->nerr + 20 < h2->nerr)  { discard = 1 ; goto laba ; }
                      if ( 2 * hh->nerr * h2->a1 < h2->nerr * hh->a1)  { discard = 1 ; goto laba ; }
                    }
                  
                  if (h2->nerr < 5)
                    {
                      if (h2->nerr + 5 < hh->nerr)  { discard = -1 ; goto laba ; }
                      if ( 3 * h2->nerr * hh->a1 < hh->nerr * h2->a1)  { discard = -1 ; goto laba ; }
                    }
                  else if (h2->nerr < 10)        
                    {  
                      if (h2->nerr + 10 < hh->nerr)  { discard = -1 ; goto laba ; }
                      if ( 2 * h2->nerr * hh->a1 < hh->nerr * h2->a1)  { discard = -1 ; goto laba ; }
                    }
                  else
                    {  
                      if (h2->nerr + 20 < hh->nerr)  { discard = -1 ; goto laba ; }
                      if ( 2 * h2->nerr * hh->a1 < hh->nerr * h2->a1)  { discard = -1 ; goto laba ; }
                    }
                }

            laba:
              if (discard == 1)
                { ne += assignPreRemove (hitsR, h2, cDNA, debug, name(hh->gene)) ; }
              else if (discard == -1)
                { ne += assignPreRemove (hitsR, hh, cDNA, debug, name (h2->gene)) ; }
            }
        }

      for (ii = 0, nkept = 0 ; ii < arrayMax(hits) ; ii++) 
        {
          hh = arrp (hits, ii, HIT) ;
          if (hh->gene)
            nkept++ ;
        }
      if (nkept >= 1 && debug)
        for (ii = 0; ii < arrayMax(hits) ; ii++) 
          {
            hh = arrp (hits, ii, HIT) ;
            if (hh->gene)
              printf("Assign nkept=%d  ::  ali=%d q=%d  %s in %s \n", nkept, hh->a1, hh->clipEnd, name(cDNA), name(hh->gene)) ;
          }
      
      if (nkept >= 9)  /* repeated element, no best place, discard from everywhere */
        {
          for (ii = 0 ; ii < arrayMax(hits) ; ii++)
            { 
              hh = arrp (hits, ii, HIT) ;
              if (hh->nerr >= 7 && hh->gene)
                {
                  ne += assignPreRemove (hitsR, hh, cDNA, debug, "over 9 and q >= 7") ; 
                }
            }
        }
    }
  genes = assignDoRemove (hitsR) ;
  chronoReturn() ;
  chrono ("assign kill dead fishes") ;

  messout ("// %s: Analysed %d cDNA, eliminated %d clones from %d genes",
           timeShowNow(), keySetMax(cdnaKs), ne, keySetMax(genes)) ;

  /* to avoid naming problems, first kill */
  ii = keySetMax (genes) ; neg = 0 ;
  if (1)  while (ii--)
    {
      gene = keySet (genes, ii) ;
      if (!keyGetKey (gene, _cDNA_clone))
        {
          cdnaCleanUp (0, gene, 0) ;
          if ((Gene = bsUpdate (gene)))
            { bsKill (Gene) ; neg++ ; }
          if (debug)        
            printf("Assign %s killed\n", name(gene)) ;
          keySet (genes, ii) = 0 ;
        }
      else
        {
          if (debug)
            {        
              printf("Assign %s kept ", name(gene)) ; 
              if (!keyFindTag (gene, _Splicing))
                printf (" NO splicing\n") ;
              else
                printf (" YES splicing\n") ;
            }
        }
    }
  chronoReturn() ;
  chrono ("assign recalculate") ;

  /* then recalculate, this will reuse small numbers */
  ii = keySetMax (genes) ;
  if (1) while (ii--)
    {
      if (0 && ii != 234) continue ;
      gene = keySet (genes, ii) ;
      if (gene)
        {
          reads = getReads (gene) ;
          if (keySetMax(reads))
            {
              int i, j ;
              KEYSET reads4, reads3, reads2 = query (reads, ">from_gene;>cdna_clone;>read;from_gene") ;
              nrecomputed++ ;
             
              newGene = cDNARealignGene (gene, 0, 0, 0, 0, 2, 0, 0) ; /* no fuse, locally, no repeats */
              if (nrecomputed%1000 == 0)
                {
                  if (debug)
                    {
                      printf("Assign  nrecomputed=%d\n",  nrecomputed) ;
                      tStatus() ; 
                      printf("Assign  nrecomputed=%d cleaning the cache\n",  nrecomputed) ;
                    }
                  alignEst2Cosmid (0, 0, 0, 0, 9999, 0, 0, 0, 0, 0) ;   /* cleanup */
                  getSeqDna ( KEYMAKE (_VCalcul, 12345)) ;    /* cleanup */ 
                  if (debug)
                    {
                      tStatus() ; 
                      printf("Saving\n") ;
                    }
                  sessionDoSave (TRUE) ;
                }
              if (newGene)
                {
                  if (debug)
                    printf("Assign %s realigned -> %s %s\n", 
                           name(gene), name(newGene), timeShowNow()) ;
                }
              else
                {
                  ndestroyedbyrecomputing++ ;
                  if (debug )
                    printf("Assign %s destroyed by recomputing %s\n", 
                           name(gene), timeShowNow()) ;
                }
              /* see if we included new 3' ends */
              reads3 = query (reads, ">from_gene;>cdna_clone;>read;from_gene") ;
              reads4 = keySetMINUS (reads3, reads2) ;
              i = keySetMax(reads4), j = keySetMax(rReads) ;
              while (i--)
                keySet(rReads, j++) = keySet (reads4, i) ;
              keySetDestroy (reads2) ; keySetDestroy (reads3) ;        keySetDestroy (reads4) ;
            }
          else
            { 
              ndestroyedbyrecomputing++ ; 
              if (debug)
                printf("Assign %s killed2\n", name(gene)) ;
            }
          keySetDestroy (reads) ;
        }
    }

  chronoReturn() ;
  chrono ("assign cleanup") ;
  messout ("// cleaning up ") ;
  if (TRUE)
    {
      Array kg ;
      kg = query (0, "Find Transcribed_gene ! Splicing") ;
      keySetKill (kg) ; 
      keySetDestroy (kg) ;
      messout ("// %s: Analysed %d cDNA, eliminated %d clones from %d genes, eliminated %d genes, recomputed %d genes, destroying %d of them",
                timeShowNow(), keySetMax(cdnaKs), ne, keySetMax(genes), neg, nrecomputed, ndestroyedbyrecomputing++) ;
    }
  chronoReturn() ;
abort:
  keySetDestroy (genes) ;
  keySetDestroy (cdnaKs) ;
  arrayDestroy (aa) ;
  arrayDestroy (aa2) ;
  arrayDestroy (hits) ;

  chronoReturn () ; 
  return rReads ;
}  /* assigncDNA */

static void  assigncDNAKeySet (KEYSET originalKeySet)
{
  /* i need to do it twice because recomputing the genes sometimes recalls the 3' */
  KEYSET r2 = 0, rClones = 0, rReads = 0 ;
  
  chrono ("assigncDNAKeySet1") ;
  rReads = assignOncecDNAKeySet (originalKeySet) ;

  chronoReturn() ;
  chrono ("assigncDNAKeySet2") ;
  if (keySetMax (rReads))
    {
      keySetSort (rReads) ;
      keySetCompress (rReads) ;
      
      rClones = query (rReads, ">cDNA_clone ; COUNT read > 1") ;
      r2 = assignOncecDNAKeySet (rClones) ;
    } 
  chronoReturn() ;
  keySetDestroy (rClones) ;
  keySetDestroy (rReads) ;
  keySetDestroy (r2) ;
}

/*********************************************************************/
/*********************************************************************/

static int  buryReadInMrna (KEY tg, KEY mrna, KEYSET doBury, KEYSET fuzzyReads)
{
  int ii, jj, kk, nr = 0 ;
  Array hits = arrayCreate (1000, HIT) ;
  HIT *hh, *h2 ;
  KEYSET clones = 0, ests = 0 ;
  OBJ Est = 0 ;
  KEY clone , est, bestEst, bb ;
  int match, cMatch, cN, cQual, nerr, qual, len, offset ;
  KEY _Ref_seq = str2tag("Ref_seq") ;
  KEY _Ref_mrna = str2tag("Ref_mrna") ;
  KEYSET bset = 0 ;
  
  clones = queryKey (mrna, ">cdna_clone") ;
  for (ii = jj = 0 ; ii < keySetMax(clones) ; ii++)
    { 
      clone = keySet (clones, ii) ;
      ests = queryKey (clone, ">Read") ;
      cMatch = 0 ; cQual = 0 ; cN = 0 ;
      hh = 0 ;
      for (kk = 0 ; kk < keySetMax(ests) ; kk++)
        {
          est = keySet (ests, kk) ;
          if ((Est = bsCreate (est)))
            {
              if (bsFindKey (Est, _From_gene, tg) &&
                  bsGetData (Est, _bsRight, _Int, &len) &&
                  bsGetData (Est, _bsRight, _Int, &offset) &&
                  bsGetData (Est, _bsRight, _Int, &match) &&
                  bsGetData (Est, _bsRight, _Int, &nerr) &&
                  bsGetData (Est, _bsRight, _Int, &qual))
                {
                  bsGetKey (Est, _cDNA_clone, &clone) ;
                  
                  hh = arrayp (hits, jj++, HIT) ;
                  hh->gene = tg ; 
                  hh->cDNA_clone = clone ;
                  hh->est = est ;
                  {
                    int n = 0 ;
                    KEY dummy ;
                    if (bsGetKey (Est, _From_gene, &dummy))
                      { n++ ;} while (bsGetKey (Est, _bsDown, &dummy));
                    if (n > 1) qual += 4 ;
                  }
                  /* order a1/a2 will pick best qual then best match len */
                  qual += 4 ; /* degrade the ESTs relative to the mrnas */
                  if (keyFindTag (clone, _Anomalous_clone) &&
                      ! keyFindTag (clone, _Internal_priming_on_A_rich))
                    qual += 4 ;        
                  if (bsFindTag (Est, str2tag("mRNA_tiling")))
                    qual -= 1 ;
                  if (bsFindTag (Est, _Ref_mrna)) qual -= 2 ;
                  if (bsFindTag (Est, _Ref_seq)) qual = 1 ; /* only ref_seq may be 1 */
                  cMatch += match ;
                  cQual += qual ; cN++ ;
                  hh->a1 = qual ; hh->a2 = 1000000 - match ; 
                  hh->x1 = offset ; hh->x2 = offset + match ;
                  hh->nerr = nerr ;
                }
              bsDestroy (Est) ;
            }
        }
      /* favor read pairs */
      if (cN > 0)  cQual /= cN ;
      for (kk = jj, h2 = hh ; hh && kk > 0 && h2->cDNA_clone == clone ; h2--, kk--)
        { h2->a2 =  1000000 - cMatch ; h2->a1 = cQual ; }
      keySetDestroy (ests) ;
    }
  keySetDestroy (ests) ;

  bset = query (clones, ">read ; >buries") ;
  if (arrayMax (hits))
    {
      arraySort (hits, cDNAOrderGloballyByA1Errors) ;
      
      if (keySetMax (bset))
        bb = keySet (bset, 0) ;
      else
        { 
          hh = arrayp (hits, 0, HIT) ;
          bestEst = hh->est ;
          for (jj = 1 ; hh->a1 == 1 && jj < arrayMax (hits) - 1 ; hh++, jj++) ;
          if (hh->est && hh->a1 > 1 && hh->a1 <= 5 && !keyFindTag (hh->est, _Ref_seq))
            bestEst = hh->est ;
          lexaddkey (name(bestEst), &bb, _VClone) ;
          if ((Est = bsUpdate (bestEst)))
            {
              bsAddKey (Est, str2tag("Buries"), bb) ;
              bsSave (Est) ;
            }
        }
      ests = query (doBury, ">Read !Is_buried_under && !buries") ;
      for (ii = 0, jj = arrayMax (hits) ; ii < keySetMax (doBury) ; ii++, jj++)
        {
          hh = arrayp (hits, jj, HIT) ;
          hh->est = keySet (ests, ii) ;
        }

      for (ii = 0, jj = arrayMax (hits) ; ii < keySetMax (fuzzyReads) ; ii++, jj++)
        {
          hh = arrayp (hits, jj, HIT) ;
          hh->est = keySet (fuzzyReads, ii) ;
        }
    
      for (jj = 200 ; jj < arrayMax (hits) ; jj++)
        {
          hh = arrayp (hits,jj, HIT) ;
          est = hh->est ;
	  if (!strncmp(name(est), "U454.", 5))
	    continue ;
          if ((Est = bsUpdate (est)))
            { 
              if (!bsFindTag (Est, _Is_buried_under) &&
		  !bsFindTag (Est, _Ref_seq) &&
                  !bsFindTag (Est, str2tag("Buries")))
                {
                  nr++ ;
                  bsAddKey (Est, _Is_buried_under, bb) ;
                }
              bsSave (Est) ;
            }
        }
    }
  keySetDestroy (bset) ;
  arrayDestroy (hits) ;
  keySetDestroy (ests) ;
  return nr ;
} /* buryReadInMrna */

/*********************************************************************/
/* in active Tg set
 * look for tg with > 600  clones
 * look inside for mrna with > 200 clones
 * bury the extra clones in a quality ordered way
 * bury the clones making bad underrepresented mrnas 
 */

static int buryReadInTgKeyset (KEYSET originalKeySet)
{
  /* i need to do it twice because recomputing the genes sometimes recalls the 3' */
  KEYSET tgs, mrnas, ests, doBury = 0, fuzzyReads = 0 ;
  int ntg, itg, nmrna = 0, imrna, nr = 0, nm1 ;
  KEY tg, mrna ;

  tgs = query (originalKeySet, "CLASS Transcribed_gene && COUNT cdna_clone > 600") ;
  ntg = keySetMax (tgs) ;
  for (itg = 0 ; itg < keySetMax (tgs) ; itg++)
    {
      doBury = keySetReCreate (doBury) ; 
      fuzzyReads = keySetReCreate (fuzzyReads) ;
      tg = keySet (tgs, itg) ;
      ests = queryKey (tg, ">read") ;

      cDnaFlagOneSuspectIntrons (tg, TRUE, doBury) ;
      fuzzyReads = queryKey (tg, "{>mrna ; COUNT cdna_clone < 6 ; {fuzzy || other ;  >cdna_clone ; > read; } SETOR { >cdna_clone ; > read ; fuzzy || other ; }} SETOR { >cdna_clone ; Ignore_this_clone_automatic || Ignore_this_clone ; > read;}" ) ; /* problem all fuzzy are not declared everywhere */
      mrnas = queryKey (tg, ">mrna ; COUNT cdna_clone >= 200") ;
      nmrna += keySetMax (mrnas) ;
      nm1 = keySetMax (mrnas) ;
      if (!nm1 && keySetMax (fuzzyReads) + keySetMax (doBury) > 0)
        { 
          keySetDestroy (mrnas) ;
          mrnas = queryKey (tg, ">mrna") ;
          nm1 = keySetMax (mrnas) ;
          if (nm1 > 1) nm1 = 1 ;
        }
      for (imrna = 0 ; imrna < nm1 ; imrna++)
        {
          mrna = keySet (mrnas, imrna) ;
          nr += buryReadInMrna (tg, mrna, doBury, fuzzyReads ) ;
          doBury = keySetReCreate (doBury) ;
          fuzzyReads = keySetReCreate (fuzzyReads) ;
        }
      keySetDestroy (mrnas) ;
      keySetDestroy (ests) ;
    }
  messout ("Buried %d reads in %d mrnas in %d tg", nr, nmrna, ntg) ;
  keySetDestroy (tgs) ;
  keySetDestroy (doBury) ; 
  keySetDestroy (fuzzyReads) ;
  return nr ;
} /* buryReadInTgKeyset */

/*********************************************************************/
/*********************************************************************/

typedef struct geneAliStruct { KEY est, tg, db, bestGene ; int ln, ali, start, nerr, qual, ix, iy, bestQual ; BOOL isMrna ; BOOL reject ;} ALI ;

static int allAliOrder (const void *a, const void *b)
{
  const ALI *va = (const ALI *)a ;
  const ALI *vb = (const ALI *)b ;
  
  int dv = va->db - vb->db ;
  if (dv) return dv ;
  dv = va->est - vb->est ;
  if (dv) return dv ;
  dv = va->qual - vb->qual ;

  return - dv ;
}

static int aliOrder (const void *a, const void *b)
{
  const ALI *va = (const ALI *)a ;
  const ALI *vb = (const ALI *)b ;
  
  int dv = va->qual - vb->qual ;
  return dv ;
}


static void assignOneEstInWholeGenome (Array allAli, int *totalp, BOOL doExport, 
                                       KEYSET matrix1, KEYSET matrix2, KEYSET matrix4, 
                                       KEYSET multi1, KEYSET multi2, KEYSET multi4, 
                                       KEY est)
{
  static Array units = 0 ;
  int ii, jj, bestQual ;
  OBJ obj = 0 ;
  static Array alis = 0 ;
  BSunit *uu ;
  ALI *al ;
  int ix, iy ;
  BOOL isMrna = FALSE ;

  units = arrayReCreate (units, 100, BSunit) ;
  alis = arrayReCreate (alis, 100, ALI) ;

  obj = bsCreate (est) ;

  if (obj)
    {
      isMrna = bsFindTag (obj, str2tag ("ref_mrna")) ;

      if (bsGetArray(obj, _From_gene, units, 7)) 
        {
          for (ii = jj = 0 ; ii < arrayMax(units) ; ii += 7)
            {
              uu = arrp(units, ii, BSunit) ;
              al = arrayp (alis, jj++, ALI) ;
              
              (*totalp)++ ;
              al->est = est ;
              al->tg = uu[0].k ;
              al->ln = uu[1].i ;
              al->start = uu[2].i ;
              al->ali = uu[3].i ;
              al->nerr = uu[4].i ;
              al->qual = uu[5].i ;
              al->isMrna = isMrna ;
              if (uu[6].s && strlen(uu[6].s))
                lexaddkey (uu[5].s, &(al->db), _VCalcul) ;
              else
                al->db = 0 ;
            }
        }
      bsDestroy (obj) ;
    }
  for (ii = 0 ; ii < arrayMax(alis) ; ii++)
    {
      int m = 0 ;
      al = arrayp (alis, ii, ALI) ;
      al->qual = mrnaQualityEvaluate (al->ln, al->ali, al->nerr, al->isMrna, &ix, &iy) ;

      m =  100 * (isMrna ? 1 : 0) + 10 * ix + iy ;
      keySet (matrix1, m)++ ;
      al->ix = ix ; al->iy = iy ; al->isMrna = isMrna ;
    }
    
   
  if (arrayMax(alis))
      {
        ALI* alr ;
        BOOL reject = FALSE ;
        KEY bestGene = 0 ;
        int jj, nn1, nn2, nn4, nkept1 = 0, nkept2 = 0 ;

        arraySort (alis, aliOrder) ;
        bestQual = arr(alis, 0, ALI).qual ;
        bestGene = arr(alis, 0, ALI).tg ;
        nn1 = nn2 = nn4 = 0 ;
        if (bestQual >= 7)
          for (ii = 0 ; ii < arrayMax(alis) ; ii++)
            {
              al = arrp (alis, ii, ALI) ;
              if (al->qual <= bestQual) 
                nkept1++ ; 
              if (al->qual <= bestQual + 1) 
                nkept2++ ; 
            }

        for (ii = 0, jj = arrayMax(allAli) ; ii < arrayMax(alis) ; ii++)
          {
            reject = FALSE ;
            al = arrp (alis, ii, ALI) ;
            nn1++ ;
            if (al->qual >= bestQual + 1)
              reject = TRUE ;
            if (al->qual >= 10)
               reject = TRUE ;
            if (nkept1 >= 9 && al->qual >= 7)
              reject = TRUE ;
            if (!reject)
              {
                int m =  100 * (al->isMrna ? 1 : 0) + 10 * al->ix + al->iy ;
                keySet (matrix2, m)++ ;
                nn2++ ;
              }
            
            reject = FALSE ;
            if (al->qual >= bestQual + 2)
              reject = TRUE ;
            if (al->qual >= 11)
               reject = TRUE ;
            if (nkept2 >= 9 && al->qual >= 7)
              reject = TRUE ;

            if (doExport)
              {
                alr = arrayp (allAli, jj++, ALI) ;
                *alr = *al ; 
                alr->reject = reject ; alr->bestQual = bestQual ; alr->bestGene = bestGene ;
              }
            if (!reject)
              {
                int m =  100 * (al->isMrna ? 1 : 0) + 10 * al->ix + al->iy ;
                keySet (matrix4, m)++ ;
                nn4++ ;
              }            
          }
        if (isMrna)
          {
            switch (nn1)
              {
              case 0: break ;
              case 1: keySet(multi1, 1)++ ; break ;
              case 2: keySet(multi1, 3)++ ; break ;
              case 3: keySet(multi1, 5)++ ; break ;
              default: keySet(multi1, 7)++ ; break ;
              }
            switch (nn2)
              {
              case 0: break ;
              case 1: keySet(multi2, 1)++ ; break ;
              case 2: keySet(multi2, 3)++ ; break ;
              case 3: keySet(multi2, 5)++ ; break ;
              default: keySet(multi2, 7)++ ; break ;
              }
            switch (nn4)
              {
              case 0: break ;
              case 1: keySet(multi4, 1)++ ; break ;
              case 2: keySet(multi4, 3)++ ; break ;
              case 3: keySet(multi4, 5)++ ; break ;
              default: keySet(multi4, 7)++ ; break ;
              }
          }
        else
          {
            switch (nn1)
              {
              case 0: break ;
              case 1: keySet(multi1, 0)++ ; break ;
              case 2: keySet(multi1, 2)++ ; break ;
              case 3: keySet(multi1, 4)++ ; break ;
              default: keySet(multi1, 6)++ ; break ;
              }
            switch (nn2)
              {
              case 0: break ;
              case 1: keySet(multi2, 0)++ ; break ;
              case 2: keySet(multi2, 2)++ ; break ;
              case 3: keySet(multi2, 4)++ ; break ;
              default: keySet(multi2, 6)++ ; break ;
              }
            switch (nn4)
              {
              case 0: break ;
              case 1: keySet(multi4, 0)++ ; break ;
              case 2: keySet(multi4, 2)++ ; break ;
              case 3: keySet(multi4, 4)++ ; break ;
              default: keySet(multi4, 6)++ ; break ;
              }
          }

      }
}

/**********************/

static void assignEstInWholeGenomeExport (char *dir, Array allAli)
{
  FILE *ff = 0 ;
  int ii ;
  KEY oldDb = 0 ;
  ALI *al ;
  char *cp, *cq ;
  if (!(ff = filopen (messprintf("%s/test", dir), 0, "w")))
    goto abort ;
  filclose(ff) ; ff = 0 ;
  
  arraySort (allAli, allAliOrder) ;
  
  for (ii = 0 ; ii < arrayMax(allAli) ; ii++)
    {
      al = arrp (allAli, ii, ALI) ;
      if (!al->db)
        continue  ;
      if (oldDb != al->db)
        {
          oldDb = al->db ;
          if (ff) 
            {
              filclose(ff) ; ff = 0 ;
            }
          cp = name(al->db) ;
          cq = cp ;
          while (*cq) cq++ ; /* gets terminal zero */
          while (cq > cp && *(cq-1) != '/') cq-- ;
          ff =  filopen (messprintf("%s/%s.bestqual.ace", dir, cq), 0, "w") ;
          if (!ff) goto abort ;
        }
      fprintf (ff, "Sequence %s\nBestHit %d %s\n\n", name(al->est), al->bestQual, name(al->bestGene)) ;
      if (al->reject) 
        fprintf (ff, "Transcribed_gene %s\nRejected_read %s %d\n\n", name(al->tg), name(al->est), al->qual ) ;
    }
  
 abort:
  filclose (ff) ; ff = 0 ;
}

static void assignEstInWholeGenome (char *cp, KEYSET ks)
{
  char *dir = strnew(cp, 0) ;
  FILE *ff = 0 ;
  BOOL doExport = FALSE ;

  int ii, total = 0, ix, iy, isMrna, iMatrix;
  char *yCaption[] = {" ", "zero","<1%","<2%","<3%","<4%","<5%",">5%"};
  char *xCaption[] = {" ", ">98%", ">90%",">80%",">50%","hit"} ;
  KEYSET matrix = 0;
  KEYSET matrix1 = 0;
  KEYSET matrix2 = 0;
  KEYSET matrix4 = 0;
  KEYSET multi = 0;
  KEYSET multi1 = 0;
  KEYSET multi2 = 0;
  KEYSET multi4 = 0;

  int nnMrna, nnEst ;

  Array allAli = arrayCreate (1000000, ALI) ;
  int level ;
  if (!dir || !*dir || !(ff = filopen (messprintf("%s.quality", dir), 0, "w"))) 
    {
      messout ("Sorry, i cannot create file %s", dir) ;
      goto abort ;
    }
  level = freeOutSetFile (ff) ;

  matrix1 = keySetCreate() ;
  matrix2 = keySetCreate() ;
  matrix4 = keySetCreate() ;
  multi1 = keySetCreate() ;
  multi2 = keySetCreate() ;
  multi4 = keySetCreate() ;
  for (ii = 0 ; ii < keySetMax(ks)  ; ii++)
    assignOneEstInWholeGenome (allAli, &total, doExport, matrix1, matrix2, matrix4, 
                               multi1, multi2, multi4, keySet (ks, ii)) ;

  freeOutf ("%d est, %d alignments\n\n",
           keySetMax(ks), total) ;

  for (iMatrix = 0 ; iMatrix < 3 ; iMatrix++)
    {
      switch (iMatrix)
        {
        case 0: matrix = matrix1 ; multi = multi1 ; break ;
        case 1: matrix = matrix2 ; multi = multi2 ; break ;
        case 2: matrix = matrix4 ; multi = multi4 ; break ;
        }
      nnEst = 0 ; nnMrna = 0 ;
      for (iy = 7 ; iy >0 ; iy--)
        {
          for (isMrna = 0 ; isMrna < 2 ; isMrna++)
            { 
              for (ix = 1 ; ix < 6 ; ix++)
                {
                  int pp =  keySet (matrix, 100*isMrna + 10 * ix + iy) ;
                  if (isMrna) nnMrna += pp ;
                  else nnEst += pp ;
                }
            }
        }
      
      if (nnEst) 
        {
          freeOutf ("Matrice des qualites %d est alignments\n", nnEst) ;
          if (iMatrix == 0)  freeOutf ("All sequences\n") ;
          if (iMatrix == 1)  freeOutf ("Cleanup at 1 points\n") ;
          if (iMatrix == 2)  freeOutf ("Cleanup at 2 point\n") ;
          freeOutf ("EST:: %d(1)  %d(2) %d(3) %d(>3)\n", 
                    keySet(multi, 0),  keySet(multi, 2),  keySet(multi, 4),  keySet(multi, 6)) ;
          
          for (iy = 7 ; iy >0 ; iy--)
            {
              for (isMrna = 0 ; isMrna < 1 ; isMrna++)
                { 
                  freeOutf (" %8s ", yCaption[iy]) ;
                  for (ix = 1 ; ix < 6 ; ix++)
                    {
                      float pp =  100.0 * keySet (matrix, 100*isMrna + 10 * ix + iy) / (1 + nnEst) ;
                      freeOutf ("%8d(%03.1f)", keySet (matrix, 100*isMrna + 10 * ix + iy), pp) ;
                    }
                  freeOutf("\n") ;
                }  
              freeOutf("\n") ;
            }
          freeOutf("       ") ;
          for (ix = 1 ; ix < 6 ; ix++)  freeOutf ("%13s", xCaption[ix]) ;
          freeOutf("\n") ;
          freeOutf("\n\n") ;        
        }

      if (nnMrna)
        {
          freeOutf ("Matrice des qualites %d mRNA alignemnts\n", nnMrna) ; 
          freeOutf ("mRNA:: %d(1)  %d(2) %d(3) %d(>3)\n", 
                    keySet(multi, 1),  keySet(multi, 3),  keySet(multi, 5), keySet(multi, 7)) ;
          for (iy = 7 ; iy >0 ; iy--)
            {
              for (isMrna = 1 ; isMrna < 2 ; isMrna++)
                { 
                  freeOutf (" %8s ", yCaption[iy]) ;
                  for (ix = 1 ; ix < 6 ; ix++)
                    {
                      float pp =  100.0 * keySet (matrix, 100*isMrna + 10 * ix + iy) / (1 + nnMrna) ;
                      freeOutf ("%8d(%03.1f)", keySet (matrix, 100*isMrna + 10 * ix + iy), pp) ;
                    }
                  freeOutf("\n") ;
                } 
              freeOutf("\n") ;
            }
          freeOutf("       ") ;
          for (ix = 1 ; ix < 6 ; ix++)  freeOutf ("%13s", xCaption[ix]) ;
          freeOutf("\n") ;
          freeOutf("\n\n") ;        
        }
    }  
      
  freeOutClose (level) ;
  if (doExport)
    assignEstInWholeGenomeExport (dir, allAli) ;

 abort:
  arrayDestroy (allAli) ;
  keySetDestroy (matrix1) ;
  keySetDestroy (matrix2) ;
  keySetDestroy (matrix4) ;

  keySetDestroy (multi1) ;
  keySetDestroy (multi2) ;
  keySetDestroy (multi4) ;
  /* DO NOT keySetDestroy (ks) ; not mine */

  filclose (ff) ;
  messfree(dir) ;
}

/*********************************************************************/

static int cdnaPreassignHash  (Array words, KEYSET cosmids, int iCosmid,
                               KEY est, Array dna, Associator ass, int step, 
                               BOOL create, int limit, BOOL justTag, KEY _WordHit)
{ 
  int i, jj, n, s, nh = 0 ;
  int max = 7 ; /* sept 23/2005 was n1 < 5, now i use n1 = n2 < 7 */
  unsigned char *cq, *cp ;
  const void *vp ;
  unsigned int oligo = 0, roligo = 0, myOligo, pos1, rest1, n1, pos2, rest2, n2 ;
  unsigned int     mask = 0x3fffffff, selectMask = 0x4000 ; /* selects G or T as base 8 */
  /*
  unsigned int leftMask = 0x3fff0000, centralMask = 0xC000, rightMask = 0x3fff ; 
  */
  KEYSET ks = 0 ;
  BOOL debug = FALSE ; 
  OBJ Est = 0 ;

  cq = arrp(dna, 0, unsigned char) - 1 ;
  n = 0 ; i = arrayMax(dna) ; s = step ; jj = 0 ;

  while (n++, cq++, s--, jj++, i--)
    {
      switch (*cq)
        {  /* ATTENTION copied from B2[] in cdnainit.h, you must keep the 2 identical */
        case A_: oligo <<= 2 ; oligo |= 0x0 ; roligo >>= 2 ; roligo |= (0x3 << 28) ; break ;
        case G_: oligo <<= 2 ; oligo |= 0x1 ; roligo >>= 2 ; roligo |= (0x2 << 28) ; break ;
        case C_: oligo <<= 2 ; oligo |= 0x2 ; roligo >>= 2 ; roligo |= (0x1 << 28) ; break ;
        case T_: oligo <<= 2 ; oligo |= 0x3 ; roligo >>= 2 ; roligo |= (0x0 << 28) ; break ;
        default: n = 0 ; break ;
        }
      if (n < 15 || s > 0)
        continue ;
      oligo &= mask ;
      roligo &= mask ;
      if (oligo & selectMask)
        myOligo = oligo ;
      else
        myOligo = roligo ;
      /*
	mypos = ((myOligo & centralMask) << 16) | ((myOligo & leftMask) >> 2) | (myOligo & rightMask) ;
      */
      if (1)
        {
          if (n >= 15 && oligo && roligo && checkOligoEntropy(15, oligo)) 
            { 
              if (create) /* scanning the cosmids */
                {
                  pos1 = oligo >> 1 ; /* divide by 2 to be in half char */
                  rest1 = oligo & 0x1 ;
                  n1 = 0 ;
                  if (words)
                    {
                      cp =  arrp(words, pos1, unsigned char) ;
                      n1 = *cp ;
                      if (rest1) n1 >>= 4 ;
                      n1 &= 0x0f ;
                    }

                  pos2 = roligo >> 1 ; /* divide by 2 to be in half char */
                  rest2 = oligo & 0x1 ;
                  n2 = 0 ;
                  if (words)
                    {
                      cp =  arrp(words, pos2, unsigned char) ;
                      n2 = *cp ;
                      if (rest2) n2 >>= 4 ;
                      n2 &= 0x0f ;
                    }
                  else
                    n2 = 1 ;

                  if (n1 + n2 > 0 && n1 + n2 < max) 
                    {
                      BOOL found = FALSE ;
		      Array bucket = 0 ;
		      int iBucket = 0 ;

                      s = step ;

                      if (assFind (ass, assVoid (myOligo), &vp))
                        while (assFindNext (ass, assVoid (myOligo), &vp, &bucket, &iBucket))
                          {
                            int ic = assInt (vp) ;
                            if (ic == iCosmid)
                              { found = TRUE ; break ; }
                          }
                      if (!found)
                        assMultipleInsert (ass, assVoid(myOligo), assVoid (iCosmid)) ;  
                      nh++ ;
                    }
                }
              else /* scanning the est */
                {
		  Array bucket = 0 ;
		  int iBucket = 0 ;

                  vp = 0 ;
                  s = step ;
                  if (assFind (ass, assVoid (myOligo), &vp))
                    while (assFindNext (ass, assVoid (myOligo), &vp, &bucket, &iBucket))
                      {
                        int ic = assInt (vp) ;
                        if (ic >=0 && ic < keySetMax(cosmids))
                          {
                            if (!ks)
                              {
                                ks = arrayCreate (keySetMax(cosmids), KEY) ;
                                keySetMax(ks) = keySetMax(cosmids) ;
                              }
                            (keySet (ks, ic))++ ;
                            if (debug && (keySet(ks, ic) < 116))
                              {
                                char buf[16] ;
                                int j ;
                                for (j = 0 ; j < 15 ; j++)
                                  buf[j] = dnaDecodeChar[*(cq - 14 + j)] ;
                                buf[15] = 0 ;
                                printf ("%s\t%s\t%5d%5d %s\n", 
                                        name(est), name(keySet(cosmids, ic)) , jj, jj+14, buf) ;
                              }
                          }
                      }                  
                }
            }
        }
    }
  
  if (ks)
    {
      int nn = 0 ;

      for (i = 0 ; i < keySetMax(ks) ; i++)
        { 
          if (debug)
            printf ("Total %d hits to cosmid %s\n", keySet(ks, i), name(keySet(cosmids, i) ) );
          if (keySet(ks, i) >= limit ||
              (keySet(ks, i) > 0 && i > 0 && keySet(ks, i) + keySet(ks, i-1) >= limit) ||
              (keySet(ks, i) > 0 && i < keySetMax(ks) - 1 && keySet(ks, i) + keySet(ks, i+1) >= limit) )
            {
              KEY cos2 = keySet (cosmids, i) ;

              if (justTag)
                {
                  if ((Est = bsUpdate (est)))
                    {
                      bsAddTag (Est, _WordHit) ;
                      bsSave (Est) ;
                    }
                }
              else
                {
                  if (!nn++)
                    freeOutf ("Sequence %s\n", name(est)) ;
                  freeOutf ("WordHit %s %d\n", name(cos2), keySet(ks, i)) ;
                }
              nh++ ;
            }
        }
      if (nn && !justTag)
        freeOutf ("\n") ;
      keySetDestroy (ks) ;
    }
  return nh ;
}

/*********************************************************************/

static int cdnaPreassignEsts2Cosmid (Array words, KEYSET cosmids, KEYSET ests, 
                                     int *nep, int limit, BOOL justTag, KEY _WordHit) 
{
  Array dna = 0 ;
  Associator ass = 0 ;
  int i, j = 0, iCosmid, nh = 0, n, ne = 0 ;
  KEY *kp, cosmid ;
 
  ass = assBigCreate (1000000) ;
  iCosmid = keySetMax(cosmids) ;
  while (iCosmid--)
    {
      cosmid = keySet (cosmids, iCosmid) ;
      dna = dnaGetWithErrors(cosmid) ;
      if (dna)
        {
          if (!justTag)
            { 
              if (j++ %10 == 0) printf ("\n//") ;
              printf (" %s", name(cosmid)) ;
            }
          nh += cdnaPreassignHash (words, cosmids, iCosmid, 0, dna, ass, 10, TRUE, limit, justTag, _WordHit) ;
          
          arrayDestroy(dna) ;
        }
    }
  if (!justTag)
    printf("\n") ;
  messout ("cdnaPreassignEsts2Cosmid uses %d words in %d cosmids", 
           nh, keySetMax(cosmids)) ;
  kp = arrp (ests, 0, KEY) - 1 ; i = keySetMax(ests) ;
  while (kp++, i--)
    {
      dna = dnaGetWithErrors(*kp) ;
      if (dna)
        {
          n = cdnaPreassignHash (words, cosmids, 0, *kp, dna, ass, 1, FALSE, limit, justTag, _WordHit) ;
          if (n) { ne++ ; nh += n ; }
          arrayDestroy(dna) ;
        }
    }
  assDestroy (ass) ;
  if (nep) *nep += ne ;
  return nh ;
}

/*********************************************************************/

static void cdnaPreassignAllReads (Array words, Array precosmids, int limit, BOOL justTag)
{
  KEYSET ests = 0, cosmids = query (precosmids, "Genomic") ;
  int n, nee = 0, ne = 0, nc = 0, nh = 0 ;
  KEY _WordHit = str2tag ("WordHit") ;

  if (keySetMax (cosmids))
    {
      ests = query (0, "Find sequence Is_read") ;
      if (justTag && ! bsIsTagInClass(_VSequence, _WordHit))
        messout ("No tag WordHit in class sequence, i cannot proceed with -t option") ;
      else if (keySetMax (ests))
        {
          n = cdnaPreassignEsts2Cosmid (words, cosmids, ests, &ne, limit, justTag, _WordHit) ;
          if (n) { nc++ ; nee += ne ; nh += n ; }
        }
    }

  freeOutf ("// Preassigned %d est to %d/%d cosmids, using %d hits\n", ne, nc,keySetMax(cosmids), nh) ;
  keySetDestroy (cosmids) ;
  keySetDestroy (ests) ;
}

/*********************************************************************/
/*********************************************************************/

static void cDNAMakeOneCloneInfo (KEY key)
{
  OBJ Est = 0 ;
  KEY est = 0, clone = 0 ;
  char *cp, buf[256] ;

  if (class(key) == _VDNA)
    dnaReClass (key, &est) ;
  else if (class(key) == _VSequence)
    est = key ;
  else
    est = 0 ;

  if (est && (Est = bsUpdate (est)))
    {
      if (!bsFindTag (Est, _Genomic) && !bsFindTag (Est, _From_gene))
        {
          bsAddTag (Est, str2tag("isRead")) ;
          if (!bsFindTag(Est, str2tag("Strand"))) 
            bsAddTag (Est, _Forward) ;
          bsAddTag (Est, _cDNA) ;
          strncpy (buf, name(est), 255) ;
          cp = buf + strlen(buf) - 2 ;
          if (cp > buf && *cp == '.')
            *cp = 0 ;
          lexaddkey (buf, &clone, _VcDNA_clone) ;
          bsAddKey (Est, _cDNA_clone, clone) ;
        }
      bsSave (Est) ;
    }
}

static void cDNAMakeCloneInfo (KEYSET ks)
{
  KEYSET ks1 = query (ks, 
   "(CLASS Sequence AND DNA AND NOT Genomic) OR (CLASS DNA)") ;
  int i = keySetMax (ks1) ;

  if (!sessionGainWriteAccess())
    return ;
  while (i--)
    cDNAMakeOneCloneInfo (keySet(ks1, i)) ;
  keySetDestroy (ks1) ;
}

/*********************************************************************/
/*********************************************************************/

static void cDNADuplicateOneCloneInfo (KEY clone, int *n1p, int *n2p)
{
  KEYSET reads = queryKey (clone, ">Read") ;
  int jj, ii = keySetMax (reads) ;
  Stack s = 0 ;
  char *info, buf [1024] ;
  OBJ Clone = 0 ;
  vTXT txt = 0 ;
  KEY rr, cl, parentClone = 0 ;

  if (ii > 1 && (Clone = bsUpdate (clone)))
    {
      arraySort (reads, keySetAlphaOrder) ; /* use alphabetic ordering */
      (*n1p) += 1 ; (*n2p) += ii - 1 ;
      s = stackCreate (1000) ;
      txt = vtxtCreate () ;

      /* find the original parent */
      if (!bsGetKey (Clone, str2tag ("Parent_clone"), &parentClone))
        parentClone = clone ;
      /* remove all reads from original clone name */
      if (_Read && bsFindTag (Clone, _Read))
        bsRemove (Clone) ;
      bsSave (Clone) ;

      /* export all clone date */
      dumpKey(clone, 0, s) ;  
      info = stackText(s, 0) ;
      while (*info && *info++ != '\n') ;    /* skip name(key) */
      
      /* add back read 1 */
      vtxtPrintf (txt, "cDNA_clone \"%s\"\n", name(clone)) ;
      rr = keySet (reads, 0) ;
      vtxtPrintf (txt, "Read \"%s\"\n\n", name(rr)) ;

      /* for each other read, create a new clone and add the info */
      for (ii = 1 ; ii < keySetMax (reads) ; ii++)
        {
          rr = keySet (reads, ii) ;
          for (jj = 1 ; ; jj++)
            {
              sprintf (buf, "%s__%d", name(parentClone), jj) ;
              if (! lexword2key (buf, &cl, _VcDNA_clone))
                break ;
            }
          lexaddkey (buf, &cl, _VcDNA_clone) ;                 
          vtxtPrintf (txt, "\n\ncDNA_clone \"%s\"\n", name(cl)) ;
          vtxtPrintf (txt, "Parent_clone \"%s\"\n", name(parentClone)) ;
          vtxtPrintf (txt, "Read \"%s\"\n\n", name(rr)) ;
          if (info && *info)
            {
              vtxtPrintf (txt, "\ncDNA_clone \"%s\"\n", name(cl)) ;
              vtxtPrint (txt, info) ;
            }
        }        
      /* parse all data */
      parseBuffer (vtxtPtr (txt), 0) ;
    }

  keySetDestroy (reads) ;
  stackDestroy (s) ;
  vtxtDestroy (txt) ;
} /* cDNADuplicateOneCloneInfo */

static void cDNADuplicateCloneInfo (KEYSET ks, int *n0p, int *n1p, int *n2p)
{
  KEYSET ks1 = 0 ;
  int ii ;

  if (! bsIsTagInClass (_VcDNA_clone, str2tag("Parent_clone")))
    {
      messout ("Missing tag Parent_clone in class cDNA_clone, sorry") ;
      return ; 
    }
  if (!sessionGainWriteAccess())
    return ;
  
  ks1 = query (ks, "CLASS cDNA_clone && Read") ;
  ii = *n0p = keySetMax (ks1) ; *n1p = *n2p = 0 ;
  while (ii--)
    cDNADuplicateOneCloneInfo (keySet(ks1, ii), n1p, n2p) ;
  keySetDestroy (ks1) ;
} /* cDNADuplicateCloneInfo */

/*********************************************************************/
/*********************************************************************/


void fMapcDNADoSelect (KEY k, KEYSET ks, KEYSET taceKs)
{
  KEY key ;
  char *cp = 0, *cq = 0 ;
  int type = k/10 ;
  KEYSET genes = 0 ;

  cDNAAlignInit () ;

  switch (k)
    {
    case 1:  case 2: case 21:
    case 2001: case 2002:
      alignAllEst2Cosmids(k, taceKs) ;
      break ;
    case 3: case 2003:
      if (!messPrompt ("Please enter a cosmid name:","ZK637","w"))
        break ;
      cp = freeword() ;
      if (!lexword2key(cp,&key,_VSequence))
        { messout ("Unknown Sequence, sorry") ; break ; }
      cdnaCleanUp (key, 0, 0) ;
      if ((genes = alignEst2Cosmid (key, 0, genes, type/1000, type%1000, 0, 0, 0, 3, 0))) /* no repeats */
        {
          mrnaAnalyseNeighbours (genes) ;
          keySetDestroy (genes) ;
        }
      break ;
    case 301:
      {
        int step = 10 ;
        int comb = 0 ;

        cq = 0 ;
        while ((cp = freeword()))
          {
            if (!strcmp (cp, "-f"))
              {
                cp = freeword() ;
                if (!cp)
                  { messout ("Usage cDNA_p1 -f frequency_file [-s step]  [-c comb] // filename missing, sorry") ; break ; }
                cq = strnew (cp, 0) ;
              }
            if (!strcmp (cp, "-s"))
              {
                freenext () ;
                if (!freeint (&step))
                  { messout ("Usage cDNA_p1 -f frequency_file [-s step]  [-c comb] // step value missing, sorry") ; break ; }
              }
            if (!strcmp (cp, "-c"))
              {
                freenext () ;
                if (!freeint (&comb))
                  { messout ("Usage cDNA_p1 -f frequency_file [-s step] [-c comb] // comb value missing, sorry") ; break ; }
              }
          }
        if (!cq)
          { messout ("Usage cDNA_p1 -f frequency_file  [-s step]  [-c comb] // filename missing, sorry") ; }
        else if (!keySetExists(taceKs))
          { messout ("no active genomic set, sorry") ; break ; }
        else
          {
            Array words = dnaGetWordUsage (taceKs, step, 0, 0, TRUE, cq, comb) ;
            arrayDestroy (words) ;
          }
        messfree (cq) ;
      }
      break ;
    case 302:
      {
        int limit = 3 ; /* min number of hit to consider */
        int outLevel = 0 ;
        BOOL justTag = FALSE ;
        FILE *f = 0 ;

        cq = 0 ;

        while ((cp = freeword()))
          {
            if (!strcmp (cp, "-f"))
              {
                cp = freeword() ;
                if (!cp)
                  { messout ("Usage cDNA_p2 -f frequency_file [-t] -l limit -o outfile // filename missing, sorry") ; break ; }
                cq = strnew (cp, 0) ;
              }
            else if (!strcmp (cp, "-t"))
              {
                justTag = TRUE ;
              }
            else if (!strcmp (cp, "-o") && !outLevel)
              {
                char *cp2 ;
                cp = freeword() ;
                if (!cp)
                  { messout ("Usage cDNA_p2 -f frequency_file  [-t] -l limit -o outfile // outfilename missing, sorry") ; break ; }
                cp2 = strnew (cp, 0) ;
                if (! (f = filopen (cp2, 0, "w")))
                  { messout ("Sorry, cannot open outfile %s", cp2) ; break ; }
                else
                  outLevel = freeOutSetFile (f) ;
                messfree (cp2) ;
              }
            else if (!strcmp (cp, "-l"))
              {
                freenext () ;
                if (!freeint (&limit))
                  { messout ("Usage cDNA_p2 -f frequency_file [-t] -l limit -o outfile // limit value missing, sorry") ; break ; }
              }
          }
        /*
        if (!cq)
          messout ("Usage cDNA_p2 -f frequency_file -l limit // filename missing, sorry") ; 
        else
        */
        if (outLevel || justTag)
          {
            if (!keySetExists(taceKs))
              { messout ("no active genomic set, sorry") ; break ; }
            else
              {
                unsigned long wnb = 0, wnw = 0;
                Array words = cq && strcasecmp (cq, "Self") ? dnaGetWordUsage (0, 0, 0, 0, FALSE, cq, 0) :  
                  dnaGetWordUsage (taceKs, 30, &wnb, &wnw, TRUE, 0, 0) ;
                messout ("calling cdnaPreassignAllReads file_name=%s words=%s wnb=%lu wnw=%lu limit=%d"
                         ,cq ? cq : "NULL", words ? "TRUE" : "NULL"
                         , wnb, wnw, limit) ;
                if (!cq || words)
                  cdnaPreassignAllReads (words, taceKs, limit, justTag) ;
                arrayDestroy (words) ;
              }
          }
        messfree (cq) ;
        if (outLevel)
          freeOutClose (outLevel) ;
        if (f)
          filclose (f) ;
      }
      break ;
    case 31: 
      assigncDNAKeySet(0) ;
      break ;      
    case 32:  
      if (!keySetExists(taceKs))
        { messout ("no active cDNA set, sorry") ; break ; }
      assigncDNAKeySet (taceKs) ;
      break ;
    case 38:  
      if (!keySetExists(taceKs))
        { messout ("no active Tg set, sorry") ; break ; }
      buryReadInTgKeyset (taceKs) ;
      break ;
    case 39:   /* report alignment quality */
      cp = freeword() ;
      if (!cp)
        { messout ("Usage cDNA_39 dir // dir missing, sorry") ; break ; }
      if (!keySetExists(taceKs))
        { messout ("no active cDNA set, sorry") ; break ; }
      assignEstInWholeGenome (cp, taceKs) ;
      break ;       
    case 47:/* translate est */
      mrnaTranslateEst (taceKs) ;
      break ;
    case 48:  /* Transfer_pg Transfer the predicted genes to the mRNA class */
      mrnaTransferPg2PredictedMrna (taceKs) ;
      mrnaTransferRefseqMakerKeySet2Product (taceKs) ;
      break ;      
    case 49:   
      cp = freeword() ;
      mrnaAnalyseClusterAllPg (cp) ;
      break ;      
    case 50: 
      mrnaAnalyseAllNeighbours () ;
      break ;      
    case 51:
      if (!keySetExists(taceKs))
        { messout ("no active EST set, sorry") ; break ; }
      fixVector (taceKs) ; /* SlideVector */
      break ; 
    case 52: 
      if (!keySetExists(taceKs))
        { messout ("no active EST set, sorry") ; break ; }
      fixPolyA (taceKs) ; /* SlidePolyA */
      break ;      
    case 53: 
      if (!keySetExists(taceKs))
        { messout ("no active mRNA set, sorry") ; break ; }
      else
        { 
          int nc, n0, n1 ;
          abiFixLabelPolyA (taceKs, &nc, &n0, &n1) ; /* SlidePolyA */
          messout ("Found %d/%d mRNA with a valid polyA tail, flagged %d clones primed on a-rich", n1, n0, nc) ;
        }
      break ;      
    case 62:
    case 63: 
      {
        BOOL onMrna = FALSE, oneStrand = FALSE ;
        BOOL isSage = k == 62 ? TRUE : FALSE ;
        char *cp ;
        int errMax = 0 ;

        while ((cp = freeword()))
          {
            if (! strcasecmp (cp, "-on_mrna"))
              onMrna = TRUE ;
            else if (!strcmp (cp, "-errMax") &&
                     freeint (&errMax)) ;
            else if (!strcmp (cp, "-oneStrand"))
              oneStrand = TRUE ;
            else
              {
                messout ("Align_primers [-on_mrna] [-oneStrand] [-errMax <int>] // bad parameter %s, sorry"
                         , cp ? cp : "missing parameter") ;
                break ;
              }
          }
        
        alignAllSage (taceKs, onMrna, isSage, errMax, oneStrand) ;
        break ;
      }      
    case 64:
      {
        char *cp ;
        int outLevel = 0 ;
        FILE *f = 0 ;
        BOOL isError = FALSE ;
        
        if  ((cp = freeword()) && *cp == '-')
          {
            if (!strcmp("-f",cp))
              { 
                if (( cp = freepath() ))
                  /* && (cq = filName (cp, 0, "w"))) bad, goes in $ACEDB */
                  {
                    if ((f = filopen (cp, 0 , "w")))
                      outLevel = freeOutSetFile (f) ;
                    else
                      {
                        isError = TRUE ;
                        messout ("// sorry, i cannot open %s", cp) ;
                      }
                  }
              }
            else
              isError = TRUE ;
          }
        if (!keySetExists(taceKs))
          { messout ("no active Rnai set, sorry") ; break ; }
        if (!isError)
          doubleAlignRnaiSet (taceKs) ;
        if (outLevel)
          freeOutClose (outLevel) ;
        if (f)
          filclose (f) ;
      }
      break ;      
    case 5:
      {
        FILE *f = 0 ;
        int outLevel = 0, is3p = 0, ishort = 0 ;
        BOOL only_alter = FALSE ;
        KEY  only_coord = 0 ;
        BOOL nonSliding = FALSE ;
        BOOL isError = FALSE ;

        while ((cp = freeword()) && *cp == '-')
          {
            if (!strcmp("-nonSliding",cp))
              nonSliding = TRUE ;             
            else if (!strcmp("-f",cp))
              { 
                if (( cp = freepath() ))
                  /* && (cq = filName (cp, 0, "w"))) bad, goes in $ACEDB */
                  {
                    if ((f = fopen (cp, "w")))
                      outLevel = freeOutSetFile (f) ;
                    else
                      {
                        isError = TRUE ;
                        messout ("// sorry, i cannot open %s", cp) ;
                      }
                  }
              }
            else if (!strcmp("-alter",cp))
              only_alter = TRUE ;
            else if (!strcmp("-3p",cp) && freeint (&is3p) && is3p > 0) ;
            else if (!strcmp("-short", cp)) 
              {
                if (freeint (&ishort) && (ishort == 1 || ishort == 2)) ;
                else
                  {
                    isError = TRUE ;
                    messout ("// -short param should be 1 or 2") ;
                  }
              }                  
            else if (!strcmp("-coord", cp)) 
              {
                only_coord = 1 ;
                if ((cp = freeword()))
                  {
                    if (*cp == '-') 
                      freeback() ;
                    else
                      {
                        KEY origine = 0 ;
                        if (lexword2key (cp, &origine, _VSequence))
                          only_coord = origine ;
                      }                        
                  }
              }  
            else
                  {
                    isError = TRUE ;
                    messout ("// usage [-f file] [-nonSliding] [-alter] [-short 1|2 ] -coord]") ;
                  }
              
          }
        
        if (!isError)
          {
            if (is3p > 0)
              export3p (taceKs, is3p) ;
            else
              exportConfirmedIntrons (taceKs, only_alter, only_coord, nonSliding, 0, ishort) ;
          }
        if (outLevel)
          freeOutClose (outLevel) ;
        if (f) filclose(f) ;
      }
      break ;
    case 6:
      {
        FILE *f = 0 ;
        int outLevel = 0 ;

        while ((cp = freeword()) && *cp == '-')
          {
            if (!strcmp("-f",cp))
              { 
                if (( cp == freeword() ) &&
                    (cq = filName (cp, 0, "w")) &&
                    (f = filopen (cq, 0, "w")))
                  outLevel = freeOutSetFile (f) ;
              }
          }
        exportOligoExons () ;
        if (outLevel)
          freeOutClose (outLevel) ;
        if (f) filclose(f) ;
      }
      break ;
    case 61: 
      cDNAVirtual () ;
      break ;      
    case 71: 
    case 73: 
      {
        BOOL doFuse = TRUE, doSplitCloud = FALSE ;
	int searchRepeats = 0, locally = 0, doClean_killed_mRNA = 0 ;
        KEYSET myKs = 0 ;

        while ((cp = freeword()))
          {
            if (!strcasecmp (cp, "-locally"))
              { locally = 2 ; doFuse = FALSE ; }
            if (!strcasecmp (cp, "-fuse_locally"))
              { locally = 2 ; doFuse = TRUE ; }
            if (!strcasecmp (cp, "-clean_killed_mRNA"))
              doClean_killed_mRNA = 1 ;
            if (!strcasecmp (cp, "-split_cloud"))
              { doSplitCloud = 2 ; locally = 2 ; doFuse = FALSE ; }
            if (!strcasecmp (cp, "-repeats"))
              searchRepeats = 1 ;
            if (!strcasecmp (cp, "-rubber"))
              searchRepeats = 2 ;
          }
        if (k == 73)
          {
            if (!keySetExists(taceKs))
              { messout ("no active transcribed_gene set, sorry") ; break ; }
            else
              myKs = taceKs ;
          }
        if (doClean_killed_mRNA)
          {
            int ii, jj, i ;
            KEY tg, clo ;
            KEYSET clones = 0 ;
            OBJ Clo ;

            for (ii = jj = 0 ; ii < keySetMax (myKs) ; ii++)
              {
                tg = keySet (myKs, ii) ;
                clones = queryKey (tg, "{>cdna_clone} SETMINUS {>mrna ; >cdna_clone}") ;
                for (i = 0 ; i < keySetMax (clones) ; i++)
                  {
                    clo = keySet (clones, i) ;
                    if ((Clo = bsUpdate (clo)))
                      {
                        if (bsFindKey (Clo, _From_gene, tg))
                          bsRemove (Clo) ;
                        bsSave (Clo) ;
                      }
                  }
                if (keySetMax (clones))
                  keySet (myKs, jj++) = tg ;
                keySetDestroy (clones) ;
              }
            keySetMax (myKs) = jj ;
          }
        cDNARealignGeneKeySet (myKs, doFuse, locally, searchRepeats, doSplitCloud) ;
        break ;      
      }
    case 74: /* Create pseudo-EST precolating through the most favorable set of exons */
      {
        KEYSET tgs = taceKs ? query (taceKs, "CLASS Transcribed_gene") : 0 ;
        if (! taceKs || !keySetMax (tgs))
	  {
	    messout ("no active set of transcribed_gene, sorry") ; 
	  }
	else
	  cdnaPercolate (tgs) ;
      }
      break ;	
    case 75: /* Tag_shedded_genes Tag the principal and shedded active Genes (Genebox, not Tg) */
      if (!keySetExists(taceKs))
        { messout ("no active set of Gene boxes, sorry") ; break ; }
      {
        KEYSET genes = query (taceKs, "CLASS Transcribed_gene") ;
        if (!keySetMax (genes))
          messout ("no active set of transcribed_gene, sorry") ; 
        else
          {
            int nn = cdnaTagSheddedTranscribedGenes (genes) ;
            messout ("Tagged %d/%d transcribed genes as shedded", nn, keySetMax (genes)) ;
          }
        keySetDestroy (genes) ;
      }
    case 76: /* split in the active set of TG all XG_ exons on Xends_ signals */
      if (!keySetExists(taceKs))
	{ messout ("no active transcribed_gene set, sorry") ; break ; }
      mrnaDesignCutAllPreMrna (taceKs) ;
      break ;      
      
    case 77: /* kill premrna echo on other strand */
      mrnaDesignCleanEcho () ;
      break ;      
      
    case 80: 
       if (!keySetExists(taceKs))
         { messout ("no active sequence set, sorry") ; break ; }
      cDNAMakeCloneInfo (taceKs) ;
      break ;      
    case 81: 
      {
        int n0 = 0, n1 = 0, n2 = 0 ;
        
        if (!keySetExists(taceKs))
          { messout ("no active sequence set, sorry") ; break ; }
        cDNADuplicateCloneInfo (taceKs, &n0, &n1, &n2) ;
        messout ("Analysed %d clones, duplicated %d, creating %d new ones", n0, n1, n2) ;
        break ; 
      }     
     case 90: 
       if (!keySetExists(taceKs))
         { messout ("no active transcript set, sorry") ; break ; }
       messout ("not yet coded sorry") ;
       /* ficheAceRunKeySet (taceKs,kantorProduct,0,0);  */
      break ;      
    case 91: 
      {
        int nn = 0, nn0 = 0 ;
        if (!keySetExists(taceKs) || !(nn0 = keySetMax (taceKs)))
          { messout ("no active gene/mrna set, sorry") ; break ; }
        nn = mrnaAddKeysetKantorInfo (taceKs) ;
        messout ("Added Kantor Info to %d/%d genes/mrnas", nn,nn0) ;
      }
      break ;      
    case 92: 
      {
        int nn = 0 ;
        BOOL isCoding = FALSE ;

        cp = freeword() ; 
        if (cp && !strcmp(cp, "-coding"))
          isCoding = TRUE ;
        if (!keySetExists(taceKs))
          { messout ("no active gene set, sorry") ; break ; }
        nn = mcAddKeysetAlterSpliceDetails (taceKs, 0, isCoding) ;
        messout ("Added Splice comparison Details to %d geneboxes", nn) ;
      }
      break ;      
    case 95: 
      {
        int nn = 0 ;
        
        if (!keySetExists(taceKs))
          { messout ("no active mrna set, sorry") ; break ; }
        nn = mrnaAddKeysetTiling (taceKs) ;
        messout ("Added Tiling paths to %d mrna", nn) ;
      }
      break ;      
    case 96: 
      {
        int nn = 0 ;
        if (!keySetExists(taceKs))
          { messout ("no active gene set, sorry") ; break ; }
        nn = giwAddKeysetIntronClass (taceKs, 0) ;
        messout ("Added Intron class to %d geneboxes", nn) ;
      }
      break ;      
    case 97: 
      {
        int nn = 0 ;
        if (!keySetExists(taceKs))
          { messout ("no active gene set, sorry") ; break ; }
        nn = giwAddKeysetProbeWalls (taceKs, 0) ;
        messout ("Added probe walls to %d geneboxes", nn) ;
      }
      break ;      
    case 98: 
      {
        int nn = 0 ;
        nn = giwAddIntronHierarchy () ;/* look for cassette introns */
        messout ("Found %d introns including another one", nn) ;
      }
      break ;      
    case 1001:  
      cp = freeword() ; 
      cq = cp && *cp ? strnew(cp, 0) : 0 ;
      alignReads2Reference(k, ks, 0, cq) ;
      messfree (cq) ;
      break ;
    case 1002:   /* align all taceactive sequences to Reference genomic */
      if (keySetExists(taceKs))
        alignReads2Reference(k, taceKs, 0, cq) ;
      else
        messout ("no tace active sequence set, sorry") ; 
      break ;
    case 1003:   /* align all taceactive sequences to all genomic */
      if (keySetExists(taceKs))
        {
          KEYSET allGenomic = query (0, "Find sequence genomic") ;
          alignReads2Reference (k, taceKs, allGenomic, cq) ;
          keySetDestroy (allGenomic) ;
        }
      else
        messout ("no tace active sequence set, sorry") ; 
      break ;
    case 101: 
      cDnaFlipGeneKeySet (0) ;
      break ;      
    case 102:  
      if (!keySetExists(taceKs))
        { messout ("no active transcribed_gene set, sorry") ; break ; }
      cDnaFlipGeneKeySet (taceKs) ;
      break ;
    case 111:  /* FlagSuspect */
      {
        KEYSET genes = 0 ;
        int ng ;
        char *cp = freeword () ;
        BOOL doIgnore = cp && !strcasecmp (cp, "-ignore") ? TRUE : FALSE ;

        if (!keySetExists(taceKs))
          { messout ("no active transcribed_gene set, sorry") ; break ; }
        genes = query (taceKs, "CLASS transcribed_gene") ;

        if (!keySetMax (genes))
          messout ("no transcribed_gene in the active set, sorry") ; 
        else
          {
            
            ng = cDnaFlagSuspectIntrons (genes, doIgnore) ;
            freeOutf("// cDnaFlagSuspectIntrons realigning the %d genes\n", ng) ;
            cDNARealignGeneKeySet (genes, FALSE, 2, 0, 0) ; /* nofuse, locally, no repeats, no splitCloud */
          }
        
        keySetDestroy (genes) ;
      }
      break ;
    case 131: /* rename using phonems */
      {
        BOOL selectNewName = TRUE, renameGene = FALSE ;
        KEYSET genes = 0 ;
        int nn = 0 ;
        char fName[1000] ;
        char chromName[12] ;

        *fName = 0 ; *chromName = 0 ;
        while ((cp = freeword()))
          {
            if (!strcmp (cp, "-c"))
              { 
                if ((cp = freeword()))
                  strcpy (chromName, cp) ;
              }
            if (!strcmp (cp, "-g"))
              { 
                selectNewName = FALSE ; renameGene = TRUE ; 
              }
            else if (!strcmp (cp, "-f"))
              {
                if ((cp = freeword()))
                  {
                    char *cq = filName(cp, "", "r") ;
                    if (!cq) 
                      messout ("cannot read file %s",cp) ;
                    else
                      strcpy (fName, cp) ;
                  }
                else
                  messout ("missing file name after -f") ;
              }
          }
        if (selectNewName && !*fName)
          { messout ("missing reserved_name file, sorry") ; break ; }
          
        if (!keySetExists(taceKs))
          { messout ("no active Gene set, sorry") ; break ; }
        genes = query (taceKs, "CLASS Gene") ;
        nn = cDnaRenameGenes (taceKs, fName, chromName, selectNewName, renameGene) ;
        messout ("renames (phase %s) %d genes", chromName, nn) ;

        keySetDestroy (genes) ;
      }
      break ;
    case 121: /* rna editing */
      if (keySetExists(taceKs))
        {
          KEYSET tgs = query (taceKs, "CLASS Transcribed_gene") ;
          int i, ntg, nr, nag, nagr, nagtg ;
          KEY tg ;

          messStatus("Flag RNA editing") ;
          ntg = nr = nag = nagr = nagtg = 0  ;
          for (i = 0 ; i < keySetMax(tgs) ; i++)
            { 
              tg = keySet(tgs, i) ;
              ntg++ ;
              if (baseCallFlagRnaEditing (tg, &nr, &nag, &nagr))
                nagtg++ ;
            }
          messout ("Analysed %d genes, containing %d reads with traces, found %d a->g in %d reads in %d genes",
                   ntg, nr,nag,nagr,nagtg) ;

          keySetDestroy (tgs) ;
        }
      else
        messout ("no active transcribed_gene set, sorry") ; 
      break ;
    case 313: /* export pair of primers for geneboxes */
      /* myKs = taceKs ; */
    case 311:
      {
        FILE *f = 0 ;
        int nn, outLevel = 0 ;

        while ((cp = freeword()) && *cp == '-')
          {
            if (!strcmp("-f",cp))
              { 
                if (( cp == freeword() ) &&
                    (cq = filName (cp, 0, "w")) &&
                    (f = filopen (cq, 0, "w")))
                  outLevel = freeOutSetFile (f) ;
              }
          }
        nn = cDnaExportGeneboxPrimers (taceKs) ;
        freeOutf ("// Exported %d gene primer pairs\n", nn) ;
        if (outLevel)
          freeOutClose (outLevel) ;
        if (f) filclose(f) ;
      }
      break ;
#ifdef JUNK
    case 11003:  break ;
      if ((type == 0 && !messPrompt ("Please enter a cosmid name:","ZK637","w")) ||
          (type == 1 && !messPrompt ("Please enter a link name:","K06H7_C14B9","w")))
        break ;
      cp = freeword() ;
      if (!lexword2key(cp,&key,_VSequence))
        { messout ("Unknown Sequence, sorry") ; break ; }
      alignReads2Reference(k, ks, 0, cp) ;
      break ;
#endif
    }
}

/*********************************************************************/
/*********************************************************************/


/* cosmides danielle:  zk637, r05d3 c41g7 c05d11 f22b5  */

/*********************************************************************/
/* this code worked well to give a single trace */
#ifdef JUNK
static Array cDNAGetHitDna (KEY est, Array dna, Array dnaR, Array hits, int *fromp)
{
  int i, j1, j2 = 0, sens ;
  char *cp, *cq ;
  HIT *hh ;
  Array new = arrayCreate (1024, char) ;

  cDNASwapX (hits) ;
  arraySort (hits, cDNAOrderByX1) ;
  for (j1 = 0; j1 < arrayMax (hits); j1 ++)
    {
      hh = arrp(hits, j1, HIT) ;

      if (hh->est != est)
        continue ;
      sens = hh->a1 < hh->a2 ? 1 : -1 ;
      if (!j2)
        {  /* complete the dna with wahtever is there */
          if (sens == 1) hh->a1 -= hh->x1 - 1 ;
          else hh->a1 += hh->x1 - 1 ;
          if (hh->a1 < 1)
            {
              i = 1 - hh->a1 ;
              hh->a1 = 1 ;
              array (new, j2 + i, char) = 0 ; /* open hh */
              cq = arrp(new, j2, char) ; j2 += i ;
              while (i--) *cq++ = N_ ;
            }
          else if (hh->a1  > arrayMax(dna))
            {
              i = hh->a1 - arrayMax(dna) ;
              hh->a1 = arrayMax(dna) ;
              array (new, j2 + i, char) = 0 ; /* open hh */
              cq = arrp(new, j2, char) ; j2 += i ;
              while (i--) *cq++ = N_ ;
            }
        }
      if (fromp && sens == 1 && *fromp >= hh->a1 - 1 && *fromp <= hh->a2)
        *fromp = j2 + *fromp - hh->a1 + 1 ;
      if (fromp && sens == 11 && *fromp >= hh->a2 - 1 && *fromp <= hh->a1)
        *fromp = j2 - *fromp + hh->a1 - 1 ;
      if (sens == 1)
        { 
          i = hh->a2 - hh->a1 + 1 ;
          cp = arrp (dna, hh->a1 -1, char) ;
          array (new, j2 + i, char) = 0 ; /* open hh */
          cq = arrp(new, j2, char) ; j2 += i ;
          while (i--) *cq++ = *cp++ ;
        }
      else if (sens == -1)
        { 
          i = hh->a1 - hh->a2 + 1 ;
          cp = arrp (dna, hh->a1 - 1, char) ;
          array (new, j2 + i, char) = 0 ; /* open hh */
          cq = arrp(new, j2, char) ; j2 += i ;
          while (i--) *cq++ = complementBase(*cp--) ;
        }
    }
  if (!j2) arrayDestroy (new) ;
  return new ;
}
#endif
/*********************************************************************/
/*********************************************************************/
/* now for multitraces */

/*********************************************************************/

Array cDNAGetReferenceHits (KEY gene, int origin)
{
  Array units = arrayCreate (60, BSunit) ;
  Array hits = arrayCreate (12, HIT) ;
  OBJ Gene = 0 ;
  int j1, j2 ;
  BSunit *up ; HIT *hh ;

  cDNAAlignInit () ;
  if (gene &&
      (Gene = bsCreate (gene)) &&
      bsGetArray (Gene, _Assembled_from, units, 5))
    for (j1 = j2 = 0 ; j1 < arrayMax(units) - 4 ; j1 += 5)
      {
        up = arrp(units, j1, BSunit) ;
        hh = arrayp (hits, j2++, HIT) ;
        hh->est = up[2].k ;
        hh->a1 = up[0].i + origin ; hh->a2 = up[1].i + origin ; 
        hh->x1 = up[3].i ; hh->x2 = up[4].i ;
      }
  bsDestroy (Gene) ;
  arrayDestroy (units) ;
  if (!arrayMax(hits)) 
    arrayDestroy (hits) ;
  else
    {
      cDNASwapA (hits) ;
      arraySort (hits, cDNAOrderByA1) ;
    }
  return hits ;
}

/*********************************************************************/
/*********************************************************************/

