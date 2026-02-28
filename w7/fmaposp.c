 /*  File: fmaposp.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Interface to OSP, Oligo selection program
 * HISTORY:
 * Last edited: Dec 17 03:11 1998 (rd)
 * Created: Tue Oct 7 1996 (mieg)
 *-------------------------------------------------------------------
 */

/* @(#)fmaposp.c	1.23 2/24/98 */

#include "fmap_.h"

#include "display.h"
#include "dna.h"
#include "systags.h"
#include "query.h"
#include "dna.h"
#include "bump.h"
#include "call.h"
#include "session.h"
#include "parse.h"
#include "pick.h"  
#include "main.h" 
#include "wac/ac.h"
#include "vtxt.h"

/************************************************************/

#define OSP_S_BIT 0x8000  /* used as a number in fmapselectbox */
#define OSP_O_BIT 0x4000  /* used as a number in fmapselectbox */

static void ospSubmitPairs (void) ;
static void fMapOspControl (void) ;
static void fMapOspControlDraw (void) ;
static void ospRedraw (void) ;

static Graph ospControlGraph = 0 ;
static Array selectedOligos = 0 ;

static int 
  oligoMaxScore = 16,
  oligoMinLength = 18, 
  oligoMaxLength = 20, 
  oligoTmMin = 49, 
  oligoTmMax = 52, 
  productMinLength = 1000, 
  productMaxLength = 2000,
  productMaxScore = 80,
  scoreLimit = 10,
  fillLimit = 75 ;
static int oligoNameBox = 0 ;

/**************************************************************************/
/**************************************************************************/

static Array mySegs = 0 ;	/* Array of SEG - global because used
				   by arraySort function for Array aa
				   of segment scores.
				   This will become look->segs */

/**************************************************************************/

static int ospScoreOrder (const void *a, const void *b)
{ 
  const int x = *(const int*) a, y = *(const int*) b ;
  SEG *seg1 =  arrp (mySegs,x,SEG), *seg2 = arrp (mySegs,y,SEG) ;
  int f1 = seg1->data.i & 0xfff, f2 = seg2->data.i  & 0xfff ;

  if (x <= 0 || x > arrayMax (mySegs) ||
      y <= 0 || y > arrayMax (mySegs)) messcrash("ospScoreOrder") ;
 
  if (f1 * f2 > 0)
    return f1 > f2 ? 1 : -1 ;

  return 0 ;
} /* ospScoreOrder */

/**************************************************************************/
/* return active oligos exactly in correct order to draw them */
static KEYSET getActiveOligos (LOOK look, int from, int to)
{
  int try, i, j, jj, j1, x1, x2,  ln;
  int iSeg, nf, nks = 0, maxMask, beginMask, showUp ;
  SEG *seg ; 
  KEYSET ks = 0 ;
  BOOL up , found ;
  int ss, tm ;
  unsigned char *mask = 0 , *maskup = 0 , *mymask ;
  unsigned char j2 ;
  int ddx = to - from ;
  AC_HANDLE h = ac_new_handle () ;
  Array aa = arrayHandleCreate(200, int, h) ;

  if (ddx < 3*fillLimit) ddx = 3*fillLimit ;
  if (ddx <10) return aa ;
  beginMask = from - ddx/2 ; maxMask = 2 * ddx ;
  maxMask /= 8 ;
  if (maxMask <= 32) maxMask = 32 ;
  mask = (unsigned char*) halloc(maxMask  + 8, h) ;
  maskup = (unsigned char*) halloc(maxMask  + 8, h) ;

  found = FALSE ;
  for (iSeg = 1; iSeg <  arrayMax(look->segs) ; iSeg++)
    {
      seg = arrp(look->segs,iSeg,SEG) ;
      if  ((seg->type | 0x1) != OLIGO_UP) continue ;
      if (!seg->x1 || !seg->x2) continue ;
      x1 = seg->x1 ; x2 = seg->x2 ;
      ln = x2 - x1 + 1 ;
      if ( ln <   oligoMinLength ||
	   ln >   oligoMaxLength)
	continue ;
      tm = seg->data.i >> 16 & 0x3ff ; 
      if (tm && 
	  (tm < 10 * oligoTmMin || tm > 10 * oligoTmMax)
	  )
	continue ;

      found = TRUE ;
      ss =  seg->data.i ;
      /* selected = ss &  OSP_S_BIT ? TRUE : FALSE ; */
      up = seg->type & 0x1 ? TRUE : FALSE ;
      if (x1 > to + ddx/2 || x2 < from - ddx/2)
	continue ;
      ss &= 0xfff ;
      if (ss > scoreLimit - .1) 
	{ array(aa, arrayMax(aa), int) = iSeg ;
	continue ;
	} 
       /* now mask */
      mymask = up ? maskup : mask ;
      for (i = seg->x1 - fillLimit ; i < seg->x2 + fillLimit ; i++)
	{ j = i - beginMask ; j1 = j/8 ; j2 = 1 << (j%8) ;
	 if (j > 0 && j1 >= 0 && j1 < maxMask)
	   mymask[j1] |= j2 ;
	} 
    }
  if (!found) goto abort ;
  /* now sort by score and re add them or destroy them */
  mySegs = look->segs ;
  arraySort (aa, ospScoreOrder) ;
  for (jj = 0 ; jj < arrayMax(aa) ; jj++)
    {
      iSeg = arr(aa, jj, int) ;
      if (iSeg <= 0 || iSeg >= arrayMax(look->segs))
	messcrash("osp 3") ;
      seg = arrp(look->segs,arr(aa, jj, int),SEG) ;
      up = seg->type & 0x1 ? TRUE : FALSE ;
      mymask = up ? maskup : mask ;
      nf = fillLimit * .8 ;
      found = FALSE ;
      for (try = 0 ; try < 2 ; try++)
	{ 
	  
	  for (i = seg->x1 - fillLimit ; i < seg->x2 + fillLimit ; i++)
	    { j = i - beginMask ; j1 = j/8 ;  j2 = 1 << (j%8) ;
	      if (j > 0 && j1 >= 0 && j1 < maxMask)
		{
		  if (try)
		    mymask[j1] |= j2 ;
		  else
		    if (! (mymask[j1] & j2))
		      { if(!nf--)  { found = TRUE  ; break ;}}
		}
	    }
	  if (!found)
	    { arr(aa, jj, int) = 0 ;
	      break ;
	    }
	}
    } 
  arraySort (aa, intOrder) ;
  arrayCompress (aa) ;
  nks = 0 ; ks = keySetCreate () ;
  for (showUp = 0 ; showUp < 2 ; showUp++) 
    {
      for (iSeg = 1; iSeg <  arrayMax(look->segs) ; iSeg++)
	{
	  int dummy ;

	  seg = arrp(look->segs,iSeg,SEG) ;
	  if  ((seg->type | 0x1) != OLIGO_UP) continue ;
	  ss =  seg->data.i ;
	  /* selected = ss &  OSP_S_BIT ? TRUE : FALSE ; */
	  ss &= 0xfff ;
	  if ((ss > scoreLimit - .1) && 
	      !arrayFind (aa, &iSeg, &dummy, intOrder))
	    continue ;
	  up = (seg->type & 0x01) ? TRUE : FALSE ;
	  if ((up && !showUp) || (!up && showUp)) continue ;
	  x1 = seg->x1 ; x2 = seg->x2 ;
	  if (x2 < from || x1 > to) continue ;
	  keySet(ks, nks++) = iSeg ;
	}
    }

abort:
  ac_free (h) ;
  return ks ;
} /* getActiveOligos */


/******************************/

static void showOligos (LOOK look, float *offset)
{ 
  int from, to, box, x1, x2, x, xmax = 0,  x2max = 0, iSeg ;
  int ii, showUp ;
  SEG *seg ;
  BUMP bump = 0 ;
  float xx, y1, y2, yb ;
  BOOL up , selected, old ;
  int ss ;
  Array aa = arrayCreate(200, int) ;
  KEYSET ks = 0 ;

  from = look->zoneMin ;
  to =  look->zoneMax ;
  ks = getActiveOligos (look, from, to) ;
  if (!ks) return ;
  bump = bumpCreate (20, 0) ; showUp = 0 ;
  for (ii = 0 ; ii < keySetMax (ks) ; ii++)
    {
      iSeg = keySet(ks, ii) ;
      seg = arrp(look->segs,iSeg,SEG) ;
      ss =  seg->data.i ;
      selected = ss &  OSP_S_BIT ? TRUE : FALSE ;
      old = ss &  OSP_O_BIT ? TRUE : FALSE ;
      ss &= 0xFFF ;
      up = (seg->type & 0x01) ? TRUE : FALSE ;
      if (up && !showUp)
	{ showUp = TRUE ;
	  *offset += xmax*1.2 ; x2max += xmax*1.2 ; xmax = 0 ;   
	  bumpDestroy(bump) ;  
	  bump = bumpCreate (20, 0) ; 
	}
      x1 = seg->x1 ; x2 = seg->x2 ;
      y1 = MAP2GRAPH(look->map,x1) ;
      y2 = MAP2GRAPH(look->map,x2+1) ;
      if (y2 < topMargin  || y1 > mapGraphHeight)
	continue ;
      
      x = 0 ; yb = seg->x1 ; /* bump in seg coords, since increasing */
      x = ss/4 ;
      if (x2 <= x1 || x < 0 || yb < 0) messcrash("") ;
      bumpItem (bump, 1, x2 - x1 + 2, &x, &yb) ; 
      xx = *offset + 1.2*x ;
      box = graphBoxStart () ;
      array(look->boxIndex,box,int) = iSeg ;
      if (y1 < topMargin ) y1 =  topMargin  ;
      graphLine (xx, y1, xx, y2) ;
      if (up)
	{
	  graphLine (xx + .4, y1 +.7, xx, y1) ;
	  graphLine (xx - .4, y1 +.7, xx, y1) ;
	}
      else
	{
	  graphLine (xx + .4, y2 - .7, xx, y2) ;
	  graphLine (xx - .4, y2 - .7, xx, y2) ;
	}
      graphBoxEnd() ;
      graphBoxDraw  (box, BLACK, selected ? GREEN : old ? ORANGE : WHITE) ;
      if (xmax < x) xmax = x ;
    }
  *offset += xmax*1.2 ; x2max += xmax*1.2 ;
  if (x2max < 12)  *offset += 12 - x2max ;
  *offset +=  2 ;

  arrayDestroy (aa) ;
  bumpDestroy (bump) ;
} /* showOligos */

/***************************************************************************************/
/***************************************************************************************/

BOOL fMapOspPositionOligo (LOOK look, SEG *seg, KEY oligo, int *a1p, int *a2p)

{
  unsigned char *cp, *cq, *buf = 0, *buf2 = 0 ;
  OBJ Oligo = 0, Parent = 0 ;
  int x1, x2, pos, oligoLength ;
  int i, maxError = 0 , maxN = 0, from, len ;
  BOOL pUp = seg->type & 0x1 ? TRUE : FALSE ;
  char *motif ;
  BOOL up ;
  int y1, y2 ;
  AC_HANDLE h = ac_new_handle () ;

  Oligo = bsCreate (oligo) ;
  if (!Oligo || !look->dna) return FALSE ;
  
  if (!bsGetData(Oligo, _Sequence, _Text, &motif))
    goto abort ;
  bsDestroy (Oligo) ;
  if ((Parent = bsCreate (seg->key)))
      { 
	if (bsFindKey (Parent, _Oligo, oligo) &&
	    bsGetData (Parent, _bsRight, _Int, &y1) &&
	    bsGetData (Parent, _bsRight, _Int, &y2))
	  {
	    bsDestroy (Parent) ;
	    goto done ;
	  }
      }
 
  bsDestroy (Parent) ;
  oligoLength = strlen(motif) ;
  buf = (unsigned char *) strnew(motif, h) ;
  cp = buf - 1 ;
  while(*++cp)
    *cp = ace_lower(*cp) ;
  dnaEncodeString((char *)buf) ;
  
  /* try downwards */
  up = FALSE ;
  from = look->zoneMin ; /* seg->x1 - look->start ; */
  if (from < 0) from = 0 ;
  len = look->zoneMax - from ; /* seg->x2 - look->start - from ; */
  if (from + len >= arrayMax(look->dna))
    len =  arrayMax(look->dna) - 2 - from ;
  if (from < 0 || len < 0 || from + len >= arrayMax(look->dna))
    goto abort ;
 
  cp = arrp(look->dna,from,unsigned char) ;
  pos = dnaPickMatch (cp, len, buf, maxError, maxN) ;
  if (pos && pos >= len) invokeDebugger() ;
  if (pos && pos < len)
    goto ok ;
  /* try upwards */
  up = TRUE ;
  cp = buf + oligoLength - 1 ;
  buf2 = halloc (oligoLength + 1, h) ;
  cq = buf2 ; 

  pos = 0 ;
  i = oligoLength ;

  if (i>0)
    while (i--)
      *cq++ = complementBase(*cp--) ;
 
  
  cp = arrp(look->dna,from, unsigned char) ;
  pos = dnaPickMatch (cp, len, buf2, maxError, maxN) ;
  if (!pos || pos >= len)
    goto abort ;
  
ok:
  /* coordinates in look coord */
  x1 = pos - 1 + from ; 
  x2 = pos + oligoLength - 1 + from ; 
  /* coordinates relative to seg */
  if (a1p) /* not in case of new oligo */
    { x1 -= seg->x1 ; x2 -= seg->x1 ; }

  if (!up) { y1 = x1 + 1 ; y2 = x2 ; }
  else     { y1 = x2 ; y2 = x1 + 1 ; }
  if (pUp) 
    { y1 = len - y1 + 1 ; y2 = len - y2 + 1 ; }
  
  /* now store */

  if ((Parent = bsUpdate (seg->key)))
      { 
	bsAddKey(Parent, _Oligo, oligo) ;
	if (bsFindKey(Parent, _Oligo, oligo))
	  { 
	    bsAddData (Parent, _bsRight, _Int, &y1) ;
	    bsAddData (Parent, _bsRight, _Int, &y2) ;
	  }
	bsSave(Parent) ;
      }
done: 
  ac_free (h) ;
  if (a1p) *a1p = y1 ; 
  if (a2p) *a2p = y2 ;
  return TRUE ;

abort: 
  ac_free (h) ;
  bsDestroy (Oligo) ;
  return FALSE ;
} /* fMapOspPositionOligo */

void fMapOspDestroy (LOOK look)
{ return ;

} /* fMapOspDestroy */

/******************************/


void fMapOspShowOligo (LOOK look, float *offset)
{ 
  showOligos (look, offset) ;
} /* fMapOspShowOligo */

/***************************************************************************************/

typedef struct keyseg {
  KEY o1, o2 ;
  int i1, i2, a1, a2, b1, b2 , f1, f2, score ; 
  float Tm1, Tm2, Tm ;
} OLPR ;

/***************************************************************************************/

void fMapOspShowOligoPairs (LOOK look, float *offset)
{
  int iSeg, i, box, x1, x2, x0 , xmax = 0, ss ;
  float xx, y1, y2,yb ;
  SEG *seg ;
  BUMP  bump = 0 ;
  Array mm = look->oligopairs ;
  OLPR* mmp ;

  if (!mm)
    return ; 
  bump = bumpCreate (mapGraphWidth,0) ;

  for (iSeg = 1; iSeg <  arrayMax(look->segs) ; iSeg++)
    {
      seg = arrp(look->segs,iSeg,SEG) ;
      /*
	does not work, because pairs coord do not flip
	if  ((seg->type | 0x1) != OLIGO_PAIR_UP) continue ;
	*/
      if  ((seg->type) != OLIGO_PAIR) continue ;
      if (!seg->x1 || !seg->x2) continue ;
  
      x1 = seg->x1 ; x2 = seg->x2 ;
      if (x2 - x1 <   productMinLength ||
	  x2 - x1 >   productMaxLength)
	continue ;

      y1 = MAP2GRAPH(look->map,x1) ;
      y2 = MAP2GRAPH(look->map,x2) ;
      if (y2 < topMargin || y1 > mapGraphHeight)
	continue ;
      i = seg->data.i & 0xffff ;

      if (i < 0 || i >= arrayMax(look->oligopairs))
	  continue ;
      mmp = arrp(look->oligopairs, i, OLPR) ;
      ss = mmp->score ;  
      if (ss > productMaxScore)
	continue ;

      x0 = 0 ; yb = seg->x1 ; /* bump in seg coords, since increasing */
     
      if (x2 <= x1 || yb < 0) messcrash("") ;
      bumpItem (bump, 1, x2 - x1 + 2.0, &x0, &yb) ; 
      xx = *offset + 1.2*x0 ;
      box = graphBoxStart () ;
      array(look->boxIndex,box,int) = iSeg ;
      if (y1 < topMargin ) y1 =  topMargin  ;
      if (y2 > mapGraphHeight - .6) y2 = mapGraphHeight - .6 ;
      graphLine (xx, y1, xx, y2) ;

      if (seg->type & 0x1)  /* up */
	{
	  y1 = MAP2GRAPH(look->map,mmp->b1) ;
	  y2 = MAP2GRAPH(look->map,mmp->b2) ;
	  if (y1 > topMargin && y1 < mapGraphHeight)
	    { 
	      box = graphBoxStart () ;
	      array(look->boxIndex,box,int) = mmp->i2 ;
	      graphLine (xx, y1, xx, y2) ;
	      graphLine (xx - .4, y1, xx + .4, y1) ;
	      graphLine (xx + .4, y2 -.7, xx, y2) ;
	      graphLine (xx - .4, y2 -.7, xx, y2) ;
	      graphBoxEnd() ;
	      if (keySetFind (selectedOligos, mmp->o2, 0))
		graphBoxDraw (box, BLACK, GREEN) ;
	    }

	  y1 = MAP2GRAPH(look->map,mmp->a1) ;
	  y2 = MAP2GRAPH(look->map,mmp->a2) ;
	  if (y1 > topMargin && y1 < mapGraphHeight)
	    {
	      box = graphBoxStart () ;
	      array(look->boxIndex,box,int) = mmp->i1 ;
	      graphLine (xx, y1, xx, y2) ;
	      graphLine (xx - .4, y2, xx + .4, y2) ;
	      graphLine (xx + .4, y1 +.7, xx, y1) ;
	      graphLine (xx - .4, y1 +.7, xx, y1) ;
	      graphBoxEnd() ;
	      if (keySetFind (selectedOligos, mmp->o1, 0))
		graphBoxDraw (box, BLACK, GREEN) ;
	    }
	}
      else /* down */
	{
	  y1 = MAP2GRAPH(look->map,mmp->a1) ;
	  y2 = MAP2GRAPH(look->map,mmp->a2) ;
	  if (y1 > topMargin && y1 < mapGraphHeight)
	    { 
	      box = graphBoxStart () ;
	      array(look->boxIndex,box,int) = mmp->i1 ;
	      graphLine (xx, y1, xx, y2) ;
	      graphLine (xx - .4, y1, xx + .4, y1) ;
	      graphLine (xx + .4, y2 -.7, xx, y2) ;
	      graphLine (xx - .4, y2 -.7, xx, y2) ;
	      graphBoxEnd() ;
	      if (keySetFind (selectedOligos, mmp->o1, 0))
		graphBoxDraw (box, BLACK, GREEN) ;
	      else if (mmp->f1 &  OSP_O_BIT)
		graphBoxDraw (box, BLACK, ORANGE) ;
	    }
	  else if (y1 <= topMargin)
	    graphCircle(xx,topMargin, .4) ;

	  y1 = MAP2GRAPH(look->map,mmp->b1) ;
	  y2 = MAP2GRAPH(look->map,mmp->b2) ;
	  if (y1 > topMargin && y1 < mapGraphHeight)
	    {
	      box = graphBoxStart () ;
	      array(look->boxIndex,box,int) = mmp->i2 ;
	      graphLine (xx, y1, xx, y2) ;
	      graphLine (xx - .4, y2, xx + .4, y2) ;
	      graphLine (xx + .4, y1 +.7, xx, y1) ;
	      graphLine (xx - .4, y1 +.7, xx, y1) ;
	      graphBoxEnd() ;

	      if (keySetFind (selectedOligos, mmp->o2, 0))
		graphBoxDraw (box, BLACK, GREEN) ;  
	      else if (mmp->f2 &  OSP_O_BIT)
		graphBoxDraw (box, BLACK, ORANGE) ;
	    }
	  else if (y1 > mapGraphHeight)
	    graphCircle(xx,mapGraphHeight - .6, .4) ;
	}

      graphBoxEnd() ;
      if (xmax < x0) xmax = x0 ;
    }
  *offset += xmax*1.2 ; 
  *offset +=  2 ;

  bumpDestroy (bump) ;
  return ;
} /* fMapOspShowOligoPairs */

/***********************************************************************/

void fMapOspSelectOligoPair (LOOK look, SEG *seg) 
{
  OLPR* mpp ;
  int i = seg->data.i & 0xffff ;

  if (seg->type != OLIGO_PAIR &&
      seg->type != OLIGO_PAIR_UP)      
    return ;
  i = seg->data.i & 0xffff ;

  if (!arrayExists(look->oligopairs) || i < 0 ||
      i >= arrayMax(look->oligopairs)) 
    return ;

  mpp = arrp(look->oligopairs, i, OLPR) ;
	  
  strncpy (look->segTextBuf, 
	   messprintf(
  "%d %d (l= %d), %s %s, score %d , Tm_o1 %3.1f, Tm_o2 %3.1f, Tm_p %3.1f",
	     seg->x1 + 1, seg->x2, seg->x2 - seg->x1, 
	     name(mpp->o1), name(mpp->o2), mpp->score, 
	     mpp->Tm1, mpp->Tm2, mpp->Tm),
	   125) ;

  return;
} /* fMapOspSelectOligoPair */

/***************************************************************************************/

static int OLPROrder (const void *a, const void *b)
{ 
  int aa = ((const OLPR*)a)->a1, ab = ((const OLPR*)b)->a1 ;
  int sa = ((const OLPR*)a)->score, sb = ((const OLPR*)b)->score ;
  
  if (sa < sb) return -1 ; 
  if (sa > sb) return 1 ;
  return aa - ab ; /* equal score, sort on position */
} /* OLPROrder */

void fMapOspFindOligoPairs (LOOK look) 
{
  static KEY  _Pairwise_scores = 0 ;
  OBJ obj ;
  KEY o2 ;
  float score = 0, Tm ;
  int i, j , i2 ;
  Array mm = 0 ;
  OLPR *mmp ;
  SEG *seg , *seg2 ;
  void *v ;
  BSMARK mark = 0 ;
  Associator aa = assCreate() ;

  if (!keySetExists(selectedOligos))
    selectedOligos = keySetCreate () ;
  /* accumulate all the oligos */
  for (j = 0, i = 1 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      
      switch (seg->type)
	{
	case OLIGO:
	case OLIGO_UP:
	  j++ ;
	  assInsert (aa, assVoid(seg->key), assVoid(i)) ;
	  if (keySetFind (selectedOligos, seg->key, 0))
	    seg->data.i |=  OSP_S_BIT ;
	  break ;
	default:
	  break ;
	}
    }
  if (!j) { assDestroy (aa) ; return ; }


  if (!_Pairwise_scores)
    lexaddkey("Pairwise_scores",&_Pairwise_scores,0) ;

  mm = arrayCreate (32, OLPR) ;
  for (j = 0, i = 1 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      
      switch (seg->type)
	{
	case OLIGO:
	case OLIGO_UP:
	  if ((obj = bsCreate (seg->key)))
	    { if (bsGetKey (obj, _Pairwise_scores, &o2))
	      do
		{ 
		  mark = bsMark (obj, mark) ;
		  if (bsGetData (obj, _bsRight, _Float, &score) &&
		      assFind (aa, assVoid(o2), &v))
		    { 
		      i2 = assInt (v) ;
		      seg2 = arrp(look->segs,i2,SEG) ;
		      mmp = arrayp (mm, j++, OLPR) ;
		      mmp->i1 = i ;
		      mmp->i2 = i2 ;
		      mmp->o1 = seg->key ;
		      mmp->o2 = o2 ;
		      mmp->score = score ;
		      mmp->Tm1 = ((seg->data.i >> 16) & 0x3ff)/10.0 ;
		      mmp->Tm2 = ((seg2->data.i >> 16) & 0x3ff)/10.0 ;
		      mmp->f1 = seg->data.i ; mmp->f2 = seg2->data.i ;
		      Tm = 0 ;
		      bsGetData (obj, _bsRight, _Float, &Tm) ;
		      mmp->Tm = Tm ;
		      switch (seg->type)
			{
			case OLIGO:
			  mmp->a1 = seg->x1 ; mmp->a2 = seg->x2 ;
			  mmp->b1 = seg2->x1 ; mmp->b2 = seg2->x2 ;
			  break ;
			case OLIGO_UP:
			  mmp->a1 = seg2->x1 ; mmp->a2 = seg2->x2 ;
			  mmp->b1 = seg->x1 ; mmp->b2 = seg->x2 ;
			  break ;
			default:
			  break ;
			}
		    }
		  bsGoto (obj, mark) ;
		} while (bsGetKey (obj, _bsDown, &o2)) ;
	      messfree (mark) ;
	      bsDestroy (obj) ;
	    }
	  break ;
	default:
	  break ;
	}
    }
  assDestroy (aa) ;
  if (arrayMax(mm))
    {
      arraySort (mm, OLPROrder) ;
      arrayCompress (mm) ;

      for (j = arrayMax(look->segs), i = 0 ; i < arrayMax(mm) ; i++)
	{
	  mmp = arrayp (mm, i, OLPR) ;
	  seg = arrayp(look->segs, j++, SEG) ;
	  seg->type = OLIGO_PAIR ;
	  seg->parent = 0 ;
	  seg->key = mmp->o1 ;
	  seg->x1 = mmp->a1 ;
	  seg->x2 = mmp->b2 ;
	  seg->data.i = i ;
	}
    }
  else
    arrayDestroy (mm) ;
  look->oligopairs = mm ;

  return ;
} /* fMapOspFindOligoPairs */

/***********************************************************************/

void fMapOspInit (void)
{
  LOOK look = 0 ;

  if (graphAssFind (&GRAPH2FMAPLOOK_ASSOC, &look))
    {
      /* re-use exisitng fmaplook */
      if (look->magic != &FMAPLOOK_MAGIC)
	messcrash ("fMapOspInit found non magic fMap-look");

      mapColSetByName ("Oligos", TRUE) ;
      mapColSetByName ("Oligo_pairs", TRUE) ;  
      look->pleaseRecompute = TRUE ;
      fMapDraw(look,0) ;
    }

  fMapOspControl() ;

  return;
} /* fMapOspInit */

/***********************************************************************/
/************ Osp control graph, actual calls to osp *******************/
/***********************************************************************/

static PickFunc oldPick = 0 ;
static int selecting = 0 ;
static int selectButtonBox = 0 ;
static int renaming = 0 ;
static int renameButtonBox = 0 ;

static void popOspControl(void)
{
  if (graphActivate (ospControlGraph))
    graphPop() ;
} /* popOspControl */

static LOOK OSPGET (char *title)
{
  LOOK look ; 
  void *vp ;
  
  if (!graphActivate (ospControlGraph))
    { ospControlGraph = 0 ; return 0 ; }
  graphCheckEditors(ospControlGraph, TRUE) ;

  /*  if (selecting == 1) selecting = 0 ; */
  fMapOspControlDraw () ;
  if (!fMapActive (0,0,0,&vp)) 
    { 
      messout("Sorry, no active sequence, I cannot proceed") ;
      return 0 ;
    }
  graphUnMessage () ;
  look = (LOOK) vp ;
  if (graphActivate (look->graph))
    {
      if (oldPick) graphRegister (PICK, oldPick) ;
      if (!mapColSetByName("Oligos", -1) ||
	  !mapColSetByName("Oligo_pairs", -1))
	{ 
	  graphPop() ;
	  mapColSetByName ("Oligos", TRUE) ;
	  mapColSetByName ("Oligo_pairs", TRUE) ;  
	  look->pleaseRecompute = TRUE ;
	  fMapDraw(look,0) ;
	  messout("Sorry, the active sequence has changed,  try again") ;
	  return 0 ;
	}
    }
  return look ;
} /* OSPGET */

/******************************/

static void ospRedraw (void)
{
  LOOK look = OSPGET ("ospRedraw") ;
  if (!look) return ;

  fMapDraw(look,0) ;   
  popOspControl() ;
  fMapOspControlDraw () ;

  return;
} /* ospRedraw */


static void ospLimitUp (void)
{
  scoreLimit += 2 ;
  ospRedraw () ;

  return;
} /* ospLimitUp */

static void ospLimitDown (void)
{
  if (scoreLimit > 1) scoreLimit -= 2 ;
  ospRedraw () ;

  return;
} /* ospLimitDown */

/******************************/

static void ospFillUp (void)
{
  int n = fillLimit/2 ;

  if (n < 10) n = 10 ;
  fillLimit += n ;
  fillLimit -= fillLimit % 20 ;
  ospRedraw () ;

  return;
} /* ospFillUp */

static void ospFillDown (void)
{
  int n = fillLimit/3 ;

  if (n < 1) n = 1 ;
  fillLimit -= n ;
  fillLimit -= fillLimit % n ;
  if (fillLimit < n) fillLimit = n ;
  ospRedraw () ;

  return;
} /* ospFillDown */

/******************************/

static void ospSubmitPairs (void)
{ 
  int ii, iSeg, from = 0, to ;
  static BOOL showUp = FALSE, up ;
  char *tmpFasta = 0, *tmpOligos= 0, *cdName = 0, *codeName = 0 ;
  FILE *fil, *filFasta, *filOligos ;
  Stack sparm = 0 ;
  KEYSET  ks = 0 ;
  SEG *seg ;
  AC_HANDLE h = ac_new_handle () ;
  LOOK look = OSPGET ("fMapOspSubmit") ;

  if (!look) return ;

  from = look->zoneMin ;
  to = look->zoneMax ; 
  if (look->flag & FLAG_COMPLEMENT)
    { 
      messout 
	("Sorry, please un-reverse the dna, i can only work on the down strand") ;
      return ;
    }
  ks = getActiveOligos (look, from, to) ;
  if (!ks || keySetMax(ks) < 2)
    { messout ("First get oligos") ; keySetDestroy (ks) ; return ; }
  if (keySetMax(ks) > 100 && 
      !messQuery("This may take several minutes, do you want to proceed"))
    goto abort ;
  if (!(filFasta = filtmpopen(&tmpFasta,"w")))
    goto abort ;
  if (!(filOligos = filtmpopen(&tmpOligos,"w")))
    goto abort ;

  /* export correctly ? */
  dnaDumpFastA (look->dna, from, to, name(look->seqKey), filFasta, 0) ;
  filclose (filFasta) ;
  fprintf(filOligos, "// top strand\n") ; showUp = FALSE ;
  for (ii = 0 ; ii < keySetMax(ks) ; ii++)
    { 
      iSeg = keySet(ks, ii) ;
      seg = arrp(look->segs,iSeg,SEG) ;
      up = (seg->type & 0x01) ? TRUE : FALSE ;  
      if (up && !showUp)
	{ showUp = TRUE ;
	fprintf(filOligos, "// bottom strand\n") ;
	}
      fprintf(filOligos, 
	      "%s %d %d\n", name(seg->key), seg->x1 - from, seg->x2 - from) ;
    }
       
  filclose (filOligos) ;
  /*   messout ("Exported in %s and %s", tmpFasta, tmpOligos ) ;  */
  messStatus ("OSP: searching pairs") ;
 
  sparm = stackCreate (50) ;
  pushText(sparm," -s  ") ;  
  catText(sparm,tmpFasta) ; 
  catText(sparm," -g ") ;  
  catText(sparm,tmpOligos) ; 
  catText (sparm, messprintf(" -l %d -L %d ", productMinLength, productMaxLength)) ;

  codeName = strnew(filName("wscripts/osp2.pl","","x"), h) ;
  if (!codeName)
    { messout ("Sorry, cannot find wscripts/osp2.pl") ; goto abort ; } 

  codeName = strnew(messprintf("osp2.pl < %s", tmpOligos), h) ;
  cdName = strnew(filName("wscripts","","x"), h) ; 
  if (!cdName)
    { messout ("Sorry, cannot find dir wscripts") ; goto abort ; } 
  messStatus("Oligo  Pairs Scoring") ;

  sessionGainWriteAccess() ;
  if (!(fil = callCdScriptPipe (cdName, "osp2.csh", stackText (sparm, 0))))
    goto abort ;

  ks = keySetCreate () ;
  parseKeepGoing = TRUE ;
  parsePipe (fil, 0, ks) ;
  parseKeepGoing = FALSE ;
    
 /* report in fmap window */
  if (graphActivate (look->graph))
    {
      mapColSetByName ("Oligo_pairs", TRUE) ;
      
      look->pleaseRecompute = TRUE ;
      fMapDraw(look,0) ; 
    }

 abort:
  graphUnMessage() ;
  stackDestroy (sparm);
  if (tmpFasta) filtmpremove (tmpFasta) ;
  if (tmpOligos) filtmpremove (tmpOligos) ; 
  ac_free (h) ;
} /* ospSubmitPairs */

static void fMapOspSubmit (void)
{
  int i, from = 0, to, level ;
  OBJ obj = 0 ;
  char *tmpFasta = 0, *codeName = 0, *fnam = 0 ;
  char *seqMrna = "sequence" ;
  FILE *fil, *filFasta, *ftmp ;
  KEYSET  kSet = 0, ksold = 0, ksnew = 0 ;
  KEY _Temporary ;
  AC_HANDLE h = ac_new_handle () ;
  vTXT txt = vtxtHandleCreate (h) ;
  Graph graph = graphActive () ;
  LOOK look = OSPGET ("fMapOspSubmit") ;
  
  if (!look) return ;
 
  from = look->zoneMin ;
  to = look->zoneMax ; 
 
  if (!strcasecmp (className(look->seqKey), "mRNA"))
    seqMrna = "mRNA" ;
  if (to > from + 5000 && 
      !messQuery(
      messprintf("%s%d%s\n%s\n%s",
		 "Searching oligos on ", to - from,
		 " bases may take several minutes",
		 " You may consider reducing the active zone",
		 "Do you wish to continue ?")))
    goto abort ;
  if (!(filFasta = filtmpopen(&tmpFasta,"w")))
    goto abort ;

  dnaDumpFastA (look->dna, from, to, name(look->seqKey), filFasta, 0) ;
  filclose (filFasta) ;

  messStatus ("OSP: searching Oligos") ;

  codeName = strnew(filName("wscripts/osp1.csh","","x"), h) ;
  if (!codeName)
    { messout ("Sorry, cannot find wscripts/osp1.csh") ; goto abort ; } 
  messStatus("Oligo Searching") ;

  ftmp = filtmpopen (&fnam, "w") ;
  filclose (ftmp) ;
  vtxtPrintf (txt, "%s %s %d %d %d %d %d %s %s > %s",
	      codeName ,
	      name(look->seqKey), /* to get the oligos back in place */
	      oligoMaxScore, oligoMinLength,  oligoMaxLength ,  
	      oligoTmMin,oligoTmMax,
              seqMrna,
	      tmpFasta,
	      fnam
	      ) ; /* the sequence fatsa file */

  sessionGainWriteAccess() ;
  
  graphActivate (graph) ;

  /* 2006_07_01
   * i would prefer to use a pipe but as of now the program crashes
   * with a memory error, possibly a misshandling of a file pointer
   * always sensitive issue on linux, so i tried to use a file rather than
   * a pipe and it changes nothing, so we shoud really revert to a pipe
   * unfortunatelly purify not longer works for me at NCBI
   */

  if (0)
    {
      if (!(fil = callCdScriptPipe (codeName, "osp1.csh", vtxtPtr (txt))))
	goto abort ;
    }
  else
    callSystem (vtxtPtr (txt)) ;

  ksold = query (0, "FIND Oligo") ;
  kSet = keySetCreate () ;

  /* bugged 
     level = freesetfile (fil, "") ;
     freespecial ("\n/\\\"\t@") ;
     parseLevel(level, kSet) ; 
     */

  if (0)
    parsePipe (fil, 0, kSet) ;
  else
    {
      fil = filopen (fnam, 0, "r") ;
      level = freesetfile (fil, "") ;
      freespecial ("\n/\\\"\t@") ;
      parseLevel(level, kSet) ; /* will close fil i think */
    }
  graphActivate (graph) ;
  ksnew = keySetMINUS (kSet, ksold) ;

  /* flag just the new oligos as temporary */
  lexaddkey ("Temporary", &_Temporary, 0) ;
  i = keySetMax (ksnew) ;
  while (i--)
    {
      if ((obj = bsUpdate (keySet (ksnew, i))))
	{
	  bsAddTag (obj, _Temporary) ;
	  bsSave (obj) ;
	}
    }
  keySetDestroy (ksold) ;
  keySetDestroy (ksnew) ;
  keySetDestroy (kSet) ;
  /* report in dnacpt window */

  graphActivate (graph) ;
  look->pleaseRecompute = TRUE ;
  fMapDraw(look,0) ;
  graphActivate (graph) ;

 abort:
  graphUnMessage() ;
  ac_free (h) ;

  if (tmpFasta) filtmpremove (tmpFasta) ;  
  if (fnam) filtmpremove (fnam) ;
} /* fMapOspSubmit */

/********************************************/

static void fMapOspSelectBox (int box)
{
  int index ;
  SEG *seg ;
  int ss ;
  BOOL selected = FALSE ;
  Graph old ;
  FMAPLOOKGET("fMapOspSelectBox") ;

  if (!selecting) 
    { 
      if (oldPick) graphRegister (PICK, oldPick) ;
      return ;
    }
  if (!(index = arr(look->boxIndex, box, int)))
    return ;
  seg = arrp(look->segs, index, SEG) ;
  if ((seg->type | 0x1) != OLIGO_UP) return ;
  if (!(seg->data.i & OSP_O_BIT))
    {
      if (box == look->activeBox)
	seg->data.i ^= OSP_S_BIT ; /* toggle */ 
      ss =  seg->data.i ;
      selected = ss  &  OSP_S_BIT ? TRUE : FALSE ;
      if (selected)
	keySetInsert (selectedOligos, seg->key) ;
      else
	keySetRemove (selectedOligos, seg->key) ;
    }
  old = graphActive () ;
  if (oldPick) oldPick (box, 0, 0) ;
  graphActivate (old) ;
  graphBoxDraw  (box, BLACK, selected ? GREEN : WHITE) ; 
} /* fMapOspSelectBox */

static void fMapOspSelectEnds (void)
{ 
  LOOK look = 0 ; 
  
  selecting = 0 ;
  if (graphActivate (ospControlGraph))
    look = OSPGET ("selectEnds") ;
  if (look)
    look->isOspSelecting  = 0 ;
} /* fMapOspSelectEnds */


static void fMapOspSelectOligos (void)
{ 
  LOOK look = 0 ;
  if (isDisplayBlocked())
    {
      messout ("%s", "Finish what you are doing first") ;
      return ;
    }
  selecting = selecting ? 1 : 2 ;
  if (!selecting || !OSPGET ("fMapOspSelectOligos")) /* if (selecting == 1) resets it  to 0 */
    return ;
  if (isDisplayBlocked())
    {
      messout ("%s", "Finish what you are doing first") ;
      return ;
    }
  if(!selecting)
    { graphUnMessage(); 
      return ;
    }
  look = OSPGET ("fMapOspSelectOligos") ;
  selecting = 1 ;
  if (look) 
  {
    PickFunc tmp  = (PickFunc) graphRegister (PICK,  (GraphFunc)fMapOspSelectBox) ;
    if (tmp != (PickFunc) fMapOspSelectBox)
      oldPick = tmp ; 
    look->isOspSelecting  = 1 ;
  }
  graphActivate (ospControlGraph) ;
  graphMessage  ("Click on the oligos you want to select,\n they will turn green\n\n"
		"Reclick on the select button or remove this box when you are done") ;
  graphRegister (MESSAGE_DESTROY, fMapOspSelectEnds) ; 

  return;
} /* fMapOspSelectOligos */

static void fMapOspDiscardOligos (void)
{ 
  KEYSET ks = query (0,"FIND Oligo Temporary") ;
  int i = keySetMax (ks) ;
  OBJ obj = 0 ;
  LOOK look ;

  if (selecting == 1)
  {
    messout ("%s", "Finish what you are doing first") ;
    return ;
  }
  look = OSPGET ("fMapOspDiscardOligos") ;
  if (!look) return ;

  sessionGainWriteAccess() ;  

  while (i--)
    if ((obj = bsUpdate (keySet (ks, i))))
      bsKill (obj) ;	      
  look->pleaseRecompute = TRUE ;
  ospRedraw () ;

  return;
} /* fMapOspDiscardOligos */

static void *klookOrder = 0 ; /*mhmp 16.10.98 */
static void fMapOspOrderOligos (void)
{
  KEYSET oldSet = 0, nKs = 0 ;
  OBJ Oligo ;
  KEY oligo, _Temporary ;
  int i ;
  void *klook ;

  if (selecting == 1)
    {
      messout ("%s", "Finish what you are doing first") ;
      return ;
    }
  if (!OSPGET ("Order"))
    return ;

  sessionGainWriteAccess() ;

  lexaddkey ("Temporary", &_Temporary, 0) ;
  for (i = 0 ; i < keySetMax (selectedOligos) ; i++)
    {
      oligo = keySet (selectedOligos, i) ;
      Oligo = bsUpdate (oligo) ;
      if (Oligo)
	{ 
	  if (bsFindTag (Oligo, _Temporary))
	    bsRemove(Oligo) ;
	  bsSave (Oligo) ;
	}
    }
  /* On cherchait visiblement a garder la meme fenetre. En vain ... */  
  nKs = keySetCopy (selectedOligos) ;
  if(keySetActive(&oldSet, &klook) && (klook ==  klookOrder)) 
    {
      keySetShow (nKs, klook) ;
    }
  else
    {
      newKeySet ("Selected Oligos To be Ordered") ;
      keySetActive(&oldSet, &klook) ;
      keySetShow (nKs, klook) ;      
      klookOrder = klook ;
    }
} /* fMapOspOrderOligos */


static void fMapOspImportOligos (void)
{
  KEYSET ks = 0 , ks0 = 0, ks1 = 0 ;  
  SEG Seg ;
  int i, nfit = 0, nnofit = 0, nbadl = 0 ;

  LOOK look = OSPGET ("fMapOspImportOligos") ;
  if (!look) return ;

  if (!keySetActive(&ks, 0)) 
    { messout("First select a KeySet containing oligos") ;
    return ;
    }
  /* never destroy ks, it belongs to keySetActive */
  ks0 = query (ks, "(CLASS Oligo) AND  (Sequence = *) AND (NOT Temporary)") ; 
  ks1 = keySetAlphaHeap(ks0, keySetMax(ks0)) ; keySetDestroy (ks0) ;
  i = keySetMax (ks1) ;
  if (!i)
    { messout("First select a KeySet containing oligos") ;
    keySetDestroy (ks1) ;
    return ;
    }

  Seg.type = 0 ;
  Seg.x1 =  look->zoneMin ;
  Seg.x2 =  look->zoneMax ;
  Seg.key = look->seqKey ;
  
  i = keySetMax (ks1) ;
  while (i--)
    { 
      OBJ obj = bsCreate (keySet(ks1, i)) ;
      char *cp ;

      if (!obj) goto jump ;
      if (!bsGetData (obj, _Sequence, _Text, &cp))
	goto jump ;
      bsDestroy (obj) ;
      continue ;
    jump:
      nbadl++ ;
      bsDestroy (obj) ;
      keySet(ks1, i) = 0 ;
    }
  i = keySetMax (ks1) ;
  while (i--)
    if ( keySet(ks1,i))
      {
	if (fMapOspPositionOligo (look, &Seg, keySet(ks1, i), 0, 0)) 
	  nfit++ ; 
	else
	  nnofit++ ;
      }
  i = keySetMax (ks1) ;
  if (nfit)
    { look->pleaseRecompute = TRUE ;
    fMapDraw(look,0) ;
    popOspControl() ;
    }
  messout ("%d Oligos accepted, %d don't fit, %d do", 
	   keySetMax (ks1) - nbadl, nnofit, nfit) ;
  keySetDestroy (ks1) ;

  return;
} /* fMapOspImportOligos */

/*****************************************************************/

static MENUOPT ospControlMenu[] =
{ { graphDestroy, "Quit" },
  { help,         "Help" },
  { ospRedraw,    "Redraw" },
  { 0, 0 }
} ;


/*****************************************************************/
static BOOL checkTC (int tc)
{
  if (tc < 20 || tc > 99)
    return FALSE ;
  return TRUE ;
} /* checkTC */

static void fMapOspShowAllOldOligos (void)
{ 
  KEYSET ks = query (0,"FIND OLigo") ;

  keySetNewDisplay (ks, "All oligos") ;

  return;
} /* fMapOspShowAllOldOligos */

static void fMapOspShowFilteredOldOligos (void)
{
  KEYSET ks = 0 ;

  ks = query (0, messprintf (
  "FIND Oligo ; Tm <= %d AND Tm >= %d AND Length <= %d AND Length >= %d",
  oligoTmMax, oligoTmMin,  oligoMaxLength,  oligoMinLength)) ;
  keySetNewDisplay (ks, "Filtered oligos") ;

  return;
} /* fMapOspShowFilteredOldOligos */

static int fMapOspAliasOligo(KEY key, char* newName)
{
  char *oldName ;
  int  classe = class(key) ;
  KEY  newKey ;
  BOOL isCaseSensitive = pickList[classe & 255 ].isCaseSensitive ;
  int (*lexstrIsCasecmp)() = isCaseSensitive ?
    strcmp : strcasecmp ;

     /********** Protection *************/ 

  if (lexiskeylocked(key))
    { 
       messout("Sorry, alias fails because %s is locked elsewhere",
	      name(key)) ;
      return 2 ;
    }

      /************* Identity *************/

  newName = lexcleanup(newName, 0) ;
  oldName = name(key) ;
  if (!strcmp(oldName, newName))
    { messfree(newName) ; return 1 ; }
  if (!lexstrIsCasecmp(oldName, newName))
    { messfree(newName) ; return 1 ; }

     /************* Unknown name case *************/

  if (!lexword2key (newName, &newKey, classe))
    { messfree(newName) ; return 1 ; }
  else
    { messfree(newName) ; return 0 ; }
} /* fMapOspAliasOligo */


static void fMapOspPickOligo (KEY key)
{
  int i, pos, suffix, iNew ;
  char *cp ;
  char oligoName [35], oligoNameBegin [25] ;
  SEG *seg ;
  LOOK look = OSPGET("fMapOspPickOligo") ;

  displayRepeatBlock () ;  
  if (!keySetFind (selectedOligos, key, 0))
    return ;
  for (i = 1 ; i < arrayMax(look->segs) ; ++i)
    if (arrp(look->segs,i,SEG)->key == key)
      { 
	seg = arrp(look->segs,i,SEG) ;
	if (seg->type != OLIGO && seg->type != OLIGO_UP)
	  return ;
	strcpy(oligoName, look->oligoNameBuffer) ;
	if (seg->type & 0x01)
	  strcat (oligoName, "-") ;
	else
	  strcat (oligoName, "+") ;
	pos = COORD (look, seg->x1) ;
	if (pos < 0)
	  strcat (oligoName, "x") ;
	pos = abs(pos)/100 ;
	strcat (oligoName, messprintf("%d",pos)) ;

	strcpy (oligoNameBegin, oligoName) ;
	suffix = 1 ;
	while (!(iNew =fMapOspAliasOligo (key, oligoName))) 
	  {
	    strcpy (oligoName, oligoNameBegin) ;
	    strcat (oligoName, ".") ;
	    strcat (oligoName, messprintf("%d", suffix)) ;
	    suffix++ ;
	  }
	if (iNew == 1)
	  {
	    /* pour permettre de choisir "son" suffixe ex: z9+x12bam */
	    if (messPrompt 
		("This is the new name.\n"
		 "You can modify this name. Then type OK"
		 , oligoName,"w"))
	      {
		freenext () ;
		cp = freepos() ;
		strcpy (oligoName,cp) ;
		if (!fMapOspAliasOligo (key, cp))
		  messout ("Sorry! This name already exists.") ;
		else
		  {
		    lexAlias (&key, oligoName, TRUE, FALSE) ;
		    if (!(look->flag & FLAG_HIDE_HEADER)) 
		      { fMapReportLine (look, seg, FALSE, 0) ;
		        graphBoxDraw (look->segBox, BLACK, LIGHTBLUE) ;
		      }
		  }
	      }
	  } 
	fMapOspOrderOligos () ;
	return ; 
      }

  return;
} /* fMapOspPickOligo */


static void renameOligoDestroy (void)
{
  renaming = 0 ; 
  OSPGET("renameOligoDestroy") ;
  displayUnBlock() ;

  return;
} /* renameOligoDestroy */

static char oligoName [16] ;

static void fMapOspRenameOligo (void) 
{ 
  LOOK look ;

  if (renaming)
    { renameOligoDestroy () ;
      return ;
    }
  if (selecting == 1)
    {
      messout ("%s", "Finish what you are doing first") ;
      return ;
    }
  look = OSPGET ("fMapOspRenameOligo") ;
  if (!look) return  ;
  if (selecting == 1)
    {
      messout ("%s", "Finish what you are doing first") ;
      return ;
    }
  if (!keySetMax(selectedOligos))
    {
      messout ("%s%s", "Select first the oligos\n",
	       "you want to rename");
      return ;
    }
  if(!checkWriteAccess())
    { messout("Sorry, you do not have Write Access");
    return ;
    }
  renaming = 1 ;
  strcpy (look->oligoNameBuffer, oligoName) ;
  graphRegister (MESSAGE_DESTROY, renameOligoDestroy) ;
  if (!strlen(look->oligoNameBuffer))
     messout("%s%s",
		  "Type first the begin of the new name. ",
		  "You have also to fix the Origin ") ;
  else
    { displayBlock (fMapOspPickOligo,
		  "Pick an oligo to rename it.\n "
		  "Remove this message to cancel.\n\n"
		  "If you haven't fixed the Origin:\n"
		  "   Remove\n"
		  "   Fix the Origin.\n"
		  "   Then come back to Rename Oligo.") ;
      fMapOspOrderOligos () ; 
    }

  return;
} /* fMapOspRenameOligo */

/***********************************/

static void fMapOspControlDraw (void)
{ 
  int line = 3;

  graphFitBounds (&mapGraphWidth, &mapGraphHeight) ;

  graphClear () ;
  graphText ("An interface by D&J T-M to OSP",
	     2, line++) ;
  graphText (" the Oligo Selection Program",
	     2, line++) ;
  graphText (" from Ladeana Hillier",
	     2, line++) ;
  line += 1.2 ;
  graphLine(0,line, mapGraphWidth, line) ;
  line += 1.2 ;

  graphText ("Filters: ", 1, line++) ;
  graphText ("Oligo length", 3, line) ;
  graphIntEditor ("Min", &oligoMinLength, 18, line++, 0) ; 
  graphIntEditor ("Max", &oligoMaxLength, 18, line++, 0) ;line += .3 ;
  graphIntEditor ("Tm min", &oligoTmMin, 3, line, checkTC) ;
  graphIntEditor ("Tm max", &oligoTmMax, 21, line, checkTC) ;

  line += 2.2 ;
  graphLine(0, line, mapGraphWidth, line) ;
  line += 1.2 ;

  graphButton("List all old oligos", fMapOspShowAllOldOligos, 2,line) ; line += 2 ; 
  graphButton("List old oligos satisfying the filters", fMapOspShowFilteredOldOligos, 2,line) ;  line += 2 ;
  graphButton("Get Selected list of Oligos", fMapOspImportOligos, 2,line) ;  line += 2 ;

  line += 1 ;

  graphButton("Find New Oligos satisfying the filters",fMapOspSubmit, 2,line) ; 

  line += 2.2 ;
  graphLine(0, line, mapGraphWidth, line) ;
  line += 1.2 ;

  graphText ("Display Oligos", 2,line) ; line += 1.2 ;
  graphText(messprintf("with score <= %3d", scoreLimit),5 + 3.27, line + .7) ;
  graphButton("<", ospLimitDown,5 , line + .5) ;
  graphButton(">", ospLimitUp , 5 + 22.20, line + .5) ;
  graphText(messprintf("or every %5d bp", fillLimit * 2),5 + 3.27, line + 2.2) ;
  graphButton("<", ospFillDown, 5, line + 2.0) ;
  graphButton(">", ospFillUp , 5 + 22.20, line + 2.0) ;

  line += 4.2 ;
  graphLine(0, line, mapGraphWidth, line) ;
  line += 1.2 ;

  graphText ("Product length", 3, line) ;
  graphIntEditor ("Min", &productMinLength, 18, line++, 0) ;
  graphIntEditor ("Max", &productMaxLength, 18, line++, 0) ;
  graphIntEditor ("Max Score", &productMaxScore, 12, line, 0) ;
  line += 2.2 ;
  graphButton("Score new pairs of displayed Oligos", ospSubmitPairs, 2,line) ; 

  line += 2.2 ;
  graphLine(0,line, mapGraphWidth, line) ;
  line += 1.2 ;

  selectButtonBox = graphButton("Select/Unselect Oligos", fMapOspSelectOligos, 2,line) ; line += 2.1 ;
  if (selecting) graphBoxDraw (selectButtonBox,BLACK, YELLOW) ;
  graphButton("Order Selected Oligos", fMapOspOrderOligos, 2,line) ; line += 2.1 ;
  graphButton("Discard New OLigos, except if ordered", fMapOspDiscardOligos, 2,line) ;

  line += 2.1 ;
  graphLine(0,line, mapGraphWidth, line) ;
  line += 1.2 ;

  renameButtonBox = graphButton ("Rename Oligo:",  fMapOspRenameOligo, 2, line) ;
   if (renaming) graphBoxDraw (renameButtonBox,BLACK, YELLOW) ;
  /* mhmp 15.09.98  je dois mettre 16 / 10 !!! */
 oligoNameBox = graphTextEditor (" ",oligoName, 16, 17, line,0); 

  line += 2.1 ;
  graphButtons (ospControlMenu, 1, 1, 80) ;
  graphRedraw () ;

  return;
} /* fMapOspControlDraw */

static void fMapOspControlDestroy (void)
{
  ospControlGraph = 0 ;
  selecting = 0 ;
  oldPick = 0 ;
  selecting = 0 ;
  selectButtonBox = 0 ;
  renaming = 0 ;
  renameButtonBox = 0 ;
} /* fMapOspControlDestroy */

static void fMapOspControl (void)
{
  Graph old = graphActive () ;
  
  if (graphActivate (ospControlGraph))
    { graphPop () ;
      return ;
    }
  
  ospControlGraph = graphCreate (TEXT_SCROLL, "OSP parameters", 
				 0, 0, 0.4, 0.7) ;
  graphHelp("OSP") ;
  graphRegister (DESTROY, fMapOspControlDestroy) ;
  graphMenu (ospControlMenu) ;
  
  fMapOspControlDraw () ;
  graphActivate (old) ;

  return;
} /* fMapOspControl */

/***************************** eof ************************************/
/**********************************************************************/
