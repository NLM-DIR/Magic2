/*  File: fmapsequence.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: sequence display and manipulation for fmap
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 15 16:16 1998 (fw)
 * * Jul 27 10:34 1998 (edgrif): Add fMapInitialise() calls for externally
 *      called routines.
 * * Jul 16 10:07 1998 (edgrif): Introduce private header fmap_.h
 * Created: Sat Jul 25 22:36:47 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: fmapsequence.c,v 1.13 2017/10/16 13:54:41 mieg Exp $ */

#include "fmap_.h"

#include "chrono.h"
#include "dna.h"
#include "peptide.h"
#include "session.h"
#include "systags.h"
#include "restriction.h"

/*************************************************************/

/* Bit masks for colours in the dna coloring code. These are significant 
   in graphColorSquares. If more than one bit is set, the first colour
   in this list actually gets displayed. 
   These were in graph.h until graphColorSquares was generalised.
   The following must correspond with TINT_* #defines in fmap.h.
*/

static int fMapTints[8] = { PALERED, PALEBLUE, RED, LIGHTGRAY, 
			      MAGENTA, CYAN, LIGHTGREEN, YELLOW } ;

/*****************************************************************/
/***************** show sequences/coordinates ********************/

static void makeStartSkip (LOOK look, float offset, int *pwidth,
			   int *pstart, int *pskip)
{
  int width = *pwidth ;
  int start, skip ;
  float dy ;

 lao:
  start = ((COORD (look, look->min) - 1) / width) * width ;
  if (look->flag & FLAG_REVERSE)
    { start = look->length - look->origin - start - 1 ; 
      for (skip = width ; 
	   -skip * look->map->mag < 1.0 ; skip += width) ;
    }
  else
    { start += look->origin ;
      for (skip = width ; 
	   skip * look->map->mag < 1.0 ; skip += width) ;
    }

  if (start > look->min)
    start -= width ;

  look->dnaWidth = width ;	/* this is correct column width */

  if (skip > width && width > 3)
    width -= 3 ;		/* for the "..." */

  *pwidth = width ;
  *pstart = start ;
  *pskip = skip ;

  if (look->flag & FLAG_AUTO_WIDTH)
    { dy = skip ;
      dy = MAP2GRAPH (look->map, dy) - MAP2GRAPH (look->map, 0) ;
      if (dy * dy > 4.1 && width > 1)
	{ width /= 2 ;
	  if (width == 2) width = 1 ;
	  if (width > 2) width = 3 * (width / 3) ;
	  goto lao ;
	}
      if (skip > width && skip > 2 && (width+=3) && width < 60)
	{ width = (skip < 60 ? skip : 60) ;  /* should be widthmax */
	  if (width > 1) width = 3 * ((width + 2) / 3) ;
	  goto lao ;
	}    
    }

  look->dnaStart = start ;
  look->dnaSkip = skip ;
}

static char *geneTranslation = 0 ;
static int *geneTransPos = 0 ;  /*mhmp 11.06.98 */
static BOOL upTranslation = FALSE ;

static void showSequence (LOOK look, float offset, 
			  int width, BOOL isProtein)
{ 
  BOOL isProt = FALSE, isFirst ;
  int decalCoeff = 0 ; 
  int start, skip, end, box, tframe = 0 ; int oldStart ;
  int ii= 0, oldDecal, decalCoord = -ACEDB_MAXINT, jj = 0, delta ;
  float y ;
  int seqEnd, x ;
  register char *cp ;
  register int i ;
  char letter[2] ; 
  char *translationTable ;
  SEG *seg ;
  BOOL doComplementSurPlace = 
    (look->flag & FLAG_COMPLEMENT_SUR_PLACE)  &&  ! isProtein ;
  
  if (!look->dna)
    return ;
  translationTable = pepGetTranslationTable (look->seqKey, 0) ;
  letter[1] = 0 ;
  if (isProtein && geneTranslation)
    {
      arrayDestroy(look->minDna) ;
      arrayDestroy(look->maxDna) ;
      arrayDestroy(look->decalCoord) ;
      arrayDestroy(look->tframe) ;
      look->minDna = arrayCreate(1000, int) ;
      look->maxDna = arrayCreate(1000, int) ;
      look->decalCoord = arrayCreate(1000, int) ;
      look->tframe = arrayCreate(1000, int) ; 
      if (look->flag & FLAG_REVERSE)
	array(look->maxDna, jj++, int) = ACEDB_MAXINT ;
      else
	array(look->maxDna, jj++, int) = -ACEDB_MAXINT ;
      decalCoeff = -1 ;
      if (upTranslation)
	decalCoeff *= -1 ;
      if (look->flag & FLAG_REVERSE)      
	decalCoeff *= -1 ;    
    }
  makeStartSkip (look, offset, &width, &start, &skip) ;
  
  if (isProtein)
    { seqEnd = look->max - 2 ;
      if (!geneTranslation)
	while (start%3 != look->min%3)
	  ++start ; /*mhmp 08.10.98 ++  --> -- */
    }
  else
    seqEnd = look->max ;
  graphTextFormat (FIXED_WIDTH) ;
  /* mhmp modif 26.10 deb */
  oldStart = start ;
  for (; start < seqEnd ; start += 60)
    {
      end = start+60 ; if (end > seqEnd) end = seqEnd + 2 ;
      if (isProtein && !geneTranslation && start%3 != look->min%3) continue ;
      tframe = 0 ;
      if (isProtein && geneTranslation)
	{
	  for (i = start ; i < end ; i++)
	    if (i>look->min &&  geneTranslation[i-look->min] != '%'  &&
		geneTranslation[i-look->min] != 0  &&
		geneTranslation[i-look->min] != '-')
	      break ;
	  if (i < end) /* success */
	    { tframe = (i - start + look->min) % 3 ; if (tframe < 0 ) tframe += 3 ;}
	}
      isFirst = TRUE ; delta = 1 ;
      if (isProtein && geneTranslation)
	for (i = start ; i < end ; i +=delta)
	  if (i >= look->min)
	    {
	      if (geneTransPos[i-look->min + tframe])
		{ 
		  if (isFirst)
		    {
		      isFirst = FALSE ;
		      delta = 3 ;
		    }
		  oldDecal = decalCoord ;
		  decalCoord = COORD(look, i + tframe) + 3 * decalCoeff * geneTransPos[i-look->min + tframe] ;
		  if (oldDecal != decalCoord)
		    { 
		      array(look->minDna, ii, int) = COORD(look, i + tframe) ;
		      array(look->decalCoord, ii, int) = decalCoord ;
		      array(look->tframe, ii, int) = tframe ;
		      ii++ ;
		      if (ii > jj)  
			{array(look->maxDna, jj++, int) = COORD(look, i - 1 + tframe) ; } 
		      isProt = TRUE ;
		    }
		}
	      else
		{
		  if (isProt)
		    { 
		      array(look->maxDna, jj++, int) = COORD(look, i - 1 + tframe) ;
		      isProt = FALSE ;
		    } 
		}
	    }
    }
  if (isProtein && geneTranslation) 
    if (arrayMax(look->maxDna) == arrayMax(look->minDna))
      array(look->maxDna,arrayMax(look->maxDna) , int) = COORD(look, seqEnd + tframe) ;
  start = oldStart ;
  /* mhmp modif 26.10 fin */ 
  for (; start < seqEnd ; start += skip)
    { int s ;
      y = MAP2GRAPH(look->map, start) - 0.5 ;
      if (y < topMargin)
	y = topMargin ;
      if (y > mapGraphHeight - 1)
	y = mapGraphHeight - 1 ;         /* mhmp 08.10.98 le +2 */
      end = start+width ; if (end > seqEnd) end = seqEnd + 2 ;
      if (isProtein && !geneTranslation && start%3 != look->min%3) continue ;

  /* tframe: translation frame, needed to pick the correct codon */
      tframe = 0 ;
      if (isProtein && geneTranslation)
	{ /* ATTENTION i - look->min peut etre negatif! mhmp 12.06.98 */
	  /* mais on teste i > look->min ?!? que je remets ici */
	  for (i = start ; i < end ; i++)
	    if (i>look->min &&  geneTranslation[i-look->min] != '%'  &&
		geneTranslation[i-look->min] != 0  &&
		geneTranslation[i-look->min] != '-')
	      /*    if (i>look->min) */
	      break ;
	  if (i < end) /* success */
	    { tframe = (300000 + i - start + look->min) % 3 ; if (tframe < 0 ) tframe += 3 ;}
	for (i = start ; i < end ; i++)
	  if (i >= look->min)
	    if (geneTransPos[i-look->min + tframe])
	      {
		graphText (messprintf("%7d",geneTransPos[i-look->min + tframe]), offset, y) ;
		break ;
	      }
    	offset +=8 ; /*mhmp 12.06.98*/
	}

      box = graphBoxStart () ;
      if (isProtein)
	{
	  s = 3 ;
	  if (start + tframe >= 0)	/* fixes nasty colors at ends */
	    graphColorSquares (arrp(look->colors, start + tframe, char),
			       offset, y, (end-start)/s, s, fMapTints) ;
	  else if (end > 0)
	    graphColorSquares (arrp(look->colors, tframe, char),
			       offset-start, y, end/s, s, fMapTints) ;
	  else		/* can happen because width shortened by 3 */
	    { graphBoxEnd () ; continue ; }
	}
      else
	{  /* mieg, in dna case use seqEnd not end */
	  int end2 = end < seqEnd ? end : seqEnd ;
	  if (start + tframe >= 0)	/* fixes nasty colors at ends */
	    graphColorSquares (arrp(look->colors, start + tframe, char),
			       offset, y, (end2-start), 1, fMapTints) ;
	  else if (end2 > 0)
	    graphColorSquares (arrp(look->colors, tframe, char),
			       offset-start, y, end2, 1, fMapTints) ;
	  else		/* can happen because width shortened by 3 */
	    { graphBoxEnd () ; continue ; }
	}


      for (x = 0, i = start ; i < 0 ; i+=3, ++x) ;
      if (i < end)
	{
	  cp = arrp(look->dna, i, char) ; 
	  for ( ; i < end ; ++x)
	    { if (isProtein)
	      { 
		if (geneTranslation)
		  *letter = i < seqEnd ? geneTranslation[i-look->min + tframe] : 0 ;
		else
		  *letter =  i < seqEnd ? e_codon(cp, translationTable) : 0 ;
		cp += 3 ;
		i += 3 ;
	      }
	    else
	      { 
		if (i < seqEnd)  /* end may be seqEnd+2 for protein sake */
		  *letter = 
		    ( doComplementSurPlace ?
		      dnaDecodeChar[(int)complementBase(*cp & 15)] : 
		      dnaDecodeChar[(int)(*cp & 15)]
		      ) ;
		else *letter = 0 ;
		++cp ;
		++i ;
	      }
	    if (*letter && i > look->min) /* mhmp 16.06.98 apres i +=3 ?!? */
	      graphText (letter, offset+x, y) ;
	    }
	  if (skip > width && width > 1 && i < seqEnd)
	    graphText (isProtein ? "." : "...", offset+x, y) ;
	}
      graphBoxEnd () ;
      if (isProtein && geneTranslation) 
	offset -=8 ; /*mhmp 12.06.98*/
      graphBoxDraw (box, BLACK, TRANSPARENT) ;
      array(look->boxIndex, box, int) = arrayMax(look->segs) ;
      seg = arrayp(look->segs, arrayMax(look->segs), SEG) ;
      seg->key = isProtein ? _Peptide : _DNA ; 
      seg->type = isProtein ? 
	(geneTranslation ? 
	 (upTranslation ? TRANS_SEQ_UP : TRANS_SEQ) : PEP_SEQ) : DNA_SEQ ;
      seg->parent = 0 ;
      seg->x1 = start + tframe ; seg->x2 = end-1 + tframe ;
    }

  graphColor (BLACK) ;
  graphTextFormat (PLAIN_FORMAT) ;
}
  
static int dnaWidth = 30 ;

void fMapSetDnaWidth (void)
{ 
  int i ;
  FMAPLOOKGET("fMapSetDnaWidth") ;

  if (look->flag & FLAG_AUTO_WIDTH)
    messout("Warning, you need to turn off auto dna width "
	    "for this change to take effect.\n");

  if (graphPrompt ("Give new width: ", messprintf ("%d", dnaWidth), "iz"))
    { freeint (&i) ;
      if (i >= 3)
	{ dnaWidth = i ; /* mieg: i think we should fix it here to i - i%3 ; */
	  fMapDraw (look, 0) ;
	}
    }
}

void fMapToggleAutoWidth (void)
{
  FMAPLOOKGET("fMapToggleAutoWidth") ;

  look->flag ^= FLAG_AUTO_WIDTH ;
  fMapDraw (look, 0) ;
}

void fMapShowDNA (LOOK look, float *offset)
{ 
  showSequence (look, *offset, dnaWidth, FALSE) ; 
  if (!(look->flag & FLAG_AUTO_WIDTH))
    graphButton ("Set DNA width", fMapSetDnaWidth, *offset, topMargin+0.1) ;
  *offset += look->dnaWidth + 3 ;
}

void fMapShow3FrameTranslation (LOOK look, float *offset)
{ 
  fMapInitialise() ;
  
  dnaWidth -= (dnaWidth % 3) ; /* absolute necessity here */
  showSequence (look, *offset, dnaWidth, TRUE) ; 
  *offset += look->dnaWidth/3 + 1 ;
}
void fMapShowTransCoords (LOOK look, float *offset) /*mhmp 11.06.98*/
{  
  int width = dnaWidth ;
  int start, skip ;
  float y ;
  int start3 ;
  graphColor(BLACK) ;

  makeStartSkip (look, *offset, &width, &start, &skip) ;

  for (; start < look->max ; start += skip)
    { y = MAP2GRAPH(look->map, start) - 0.5 ;
      if (y < topMargin)
	y = topMargin ;
      if (y > mapGraphHeight - 1)
	y = mapGraphHeight - 1 ;
      start3 = COORD(look, start) ;
      start3 = start3 / 3 + 1 ;
      graphText (messprintf("%7d", start3), *offset, y) ;
    }
  *offset += 8 ;
}
void fMapShowCoords (LOOK look, float *offset)
{  
  int width = dnaWidth ;
  int start, skip ;
  float y ;

  graphColor(BLACK) ;

  makeStartSkip (look, *offset, &width, &start, &skip) ;
  for (; start < look->max ; start += skip)
    { y = MAP2GRAPH(look->map, start) - 0.5 ;
      if (y < topMargin)
	y = topMargin ;
      if (y > mapGraphHeight - 1)
	y = mapGraphHeight - 1 ;
      graphText (messprintf("%7d", COORD(look, start)), *offset, y) ;
    }
  *offset += 8 ;
}

/******************************************************************/
/*** to display the DNA as a yelllow rectangle with colored marks */


void fMapShowSummary (LOOK look, float *offset)
{
  float oldoff ;
  int i, i1 ;
  float xmin, xmax, x1, x2 , foundSites = 0 ;
  unsigned char *cp, old, older ;
  BOOL stype1 = FALSE ;
  Array sites = look->sites ; /* mhmp */
  Site *s ;

  if (!fmapView (look->view, str2tag("Fmap_Summary")))
    return ;

  if (look->summaryBox)
    { graphBoxDim (look->summaryBox, offset, &x1, &x1, &x1) ;
      if (arrayExists(sites))
	{
	  s = arrp(sites, 0, Site) - 1 ; i = arrayMax(sites); 
	  while (s++, i--)
	    if (s->type & 1)
	      { stype1 = TRUE ;
		break ;
	      }
	}
      *offset -= 0.4 ;
      if (stype1)
	*offset += 0.6 ;
      graphBoxClear (look->summaryBox) ;
    }

  look->summaryBox = graphBoxStart () ;

  xmin = MAP2GRAPH(look->map, look->min) ;
  xmax = MAP2GRAPH(look->map, look->max) ;

  graphColor(YELLOW) ;
  graphFillRectangle (*offset+0.4, xmin, *offset+1.6, xmax) ;
 
  old = older = WHITE; 
  cp = arrp(look->colors, look->min, unsigned char) ;
  for (i = i1 = look->min ; i < look->max ; i++, cp++)
    { int BMNeeded, color;
      BMNeeded = *cp & ~(TINT_HIGHLIGHT1 | TINT_HIGHLIGHT2);
      /* Ignore highlighting */
      switch(BMNeeded^(BMNeeded-1)) /* trick to find first 1 bit */
	{ 
	default: color = WHITE; break;
	case 0x07: color = RED; break;
	case 0x0f: color = DARKGRAY; break;
	case 0x1f: color = MAGENTA; break;
	case 0x3f: color = CYAN; break;
	case 0x7f: color = LIGHTGREEN; break;
	case 0xff: color = YELLOW; break;
	}
      if (color != old) 
	{ if (old != WHITE)
	    { if (old != older)
		{ older = old ;
		  graphColor(old);
		}
	      x1 = MAP2GRAPH(look->map, i1) ;
	      x2 = MAP2GRAPH(look->map, i + 1)  ;
	      graphFillRectangle (*offset+0.5, x1, *offset+1.5, x2) ;
	    }
	  old = color;
	  i1 = i ;
	}
    }
  if (old != WHITE)
    { if (old != older) graphColor(old);
      x1 = MAP2GRAPH(look->map, i1) ;
      x2 = MAP2GRAPH(look->map, i + 1)  ;
      graphFillRectangle (*offset+0.5, x1, *offset+1.5, x2) ;
    }

  x1 = MAP2GRAPH(look->map, look->zoneMin - 0.5) ;
  x2 = MAP2GRAPH(look->map, look->zoneMax - 0.5) ;
  if (x1 > x2)			/* reverse to simplify maths */
    { float tmp = x2 ; x2 = x1 ; x1 = tmp ;
      tmp = xmax ; xmax = xmin ; xmin = tmp ;
    }
  if (x1 > xmin || x2 < xmax)
    { graphColor (BLUE) ;
      if (x2 < xmin || x1 > xmax)
	graphFillRectangle (*offset+0.4, xmin, *offset+0.9, xmax) ;
      else 
	{ if (x1 > xmin)
	    graphFillRectangle (*offset+0.4, xmin, *offset+0.9, x1) ;
	  if (x2 < xmax)
	    graphFillRectangle (*offset+0.4, x2, *offset+0.9, xmax) ;
	}
    }


/******************************************************************/
/*** to display arrows along the DNA yelllow rectangle */

  if (arrayExists(sites))
    {
      s = arrp(sites, 0, Site) - 1 ; i = arrayMax(sites); 
      while (s++, i--)
	{
	  if (s->i > look->min && s->i < look->max)
	    { 
	      cp = arrp(look->colors, s->i, unsigned char) ;
	      {
		int BMNeeded = *cp & ~(TINT_HIGHLIGHT1 | TINT_HIGHLIGHT2);
		/* Ignore highlighting */
		switch(BMNeeded^(BMNeeded-1)) /* trick to find first 1 bit */
		  { 
		  default: graphColor(WHITE); break;
		  case 0x07: graphColor(RED); break;
		  case 0x0f: graphColor(DARKGRAY); break;
		  case 0x1f: graphColor(MAGENTA); break;
		  case 0x3f: graphColor(CYAN); break;
		  case 0x7f: graphColor(LIGHTGREEN); break;
		  case 0xff: graphColor(YELLOW); break;
		  }
	      }
	      x1 = MAP2GRAPH(look->map, s->i) ;
	      x2 = s->i + 1 ;
	      while (x2 < look->max && *cp == *(cp+1))
		{ x2++ ; cp++ ;}
	      x2 = MAP2GRAPH(look->map, x2) ;
	      if (s->type & 1)
		graphFillRectangle (*offset-0.2, x1, *offset+0.5, x2) ;
	      if (s->type & 2)
		{ graphFillRectangle (*offset+1.5, x1, *offset+2.2, x2) ;
		  foundSites = .5 ;
		}
	    }
	}
    }
  /***************************************/
  graphColor(BLACK) ;
  graphRectangle (*offset+0.4, xmin, *offset+1.6, xmax) ;

  graphBoxEnd () ;
  graphBoxDraw (look->summaryBox, -1, -1) ; /* make sure to redraw */
  graphBubbleInfo (look->summaryBox, 0, 0, 
		   messprintf("SummaryBar %d\t%d\t%d", look->min + 1, look->max, look->max - look->min)) ;

  *offset += 2.4 + foundSites ; 
  oldoff = *offset ;
  graphBoxDim (look->summaryBox, offset,&x1,&x1,&x1) ;
  *offset = oldoff ;
}


/******************************************************************/
/*** to display arrows along the DNA yelllow rectangle */

void fMapShowCptSites (LOOK look, float *offset)
{ Site *s ;
  Array sites = look->sites ;
  int i, maxlen = 0 ;
  char *cp ;

  if(!arrayExists(sites) || ! arrayMax(sites))
    return ;

  s = arrp(sites, 0, Site) - 1 , i = arrayMax(sites); 
  while (s++, i--)
    if (s->i > look->min && s->i < look->max)
      { cp = stackText(look->siteNames, s->mark) ;
	while(*cp == ' ') cp++ ;
	graphText(messprintf("%6.9s",cp),
		  *offset- 1.6, MAP2GRAPH(look->map, s->i + 3 ) - .5) ;
	if (strlen(cp) > maxlen)
	  maxlen = strlen(cp) ;
      }
  *offset += maxlen > 5 ? 7 : 1 + maxlen ;
}

void fMapRegisterSites(void *v, Array sites, Stack sname)
{
  LOOK look = (LOOK) v ; 
  
  fMapInitialise() ;
  
  graphActivate (look->graph) ;
  if (sites != look->sites)
    arrayDestroy(look->sites) ;
   if (sname != look->siteNames)
     stackDestroy(look->siteNames) ;
  look->sites = sites ;
  look->siteNames = sname ;
}

/* RD 970917 - removed fMapShowFinished: it relied on SEQUENCE seg->data.k
   values that have not been filled in for years.  We should generalise/
   use Mike Holman's code, if we want this back (for now we use the gmap
   for this).
*/

/***************************************************************/

void fMapReDrawDNA (void *v)
{ LOOK look = (LOOK) v ;
  int i ;
  float dummy ;
  Graph old = graphActive() ;

  fMapInitialise() ;

  graphActivate (look->graph) ;

  for (i = look->minLiveBox ; i < arrayMax(look->boxIndex) ; ++i)
    if (BOXSEG(i)->type == DNA_SEQ)
      graphBoxDraw (i, -1, -1) ;

  if (look->summaryBox)
    fMapShowSummary (look, &dummy) ;

  graphActivate (old) ;
}

/***************************************************************/
/**************** code to zip together DNA sequences ***********/

int sequenceLength (KEY seq)
{
  OBJ obj = bsCreate (seq) ;
  KEY dnaKey ;
  Array dna = 0 ;
  int i, length = 0 ;
  static Array aa = 0 ;
            
  if (!obj)
    return 0 ;

  if (bsGetKey (obj, _DNA, &dnaKey) && 
      (bsGetData (obj, _bsRight, _Int, &length) || 
       (dna = dnaGet (dnaKey))))
    { if (dna)
	{ length = arrayMax(dna) ;
	  arrayDestroy (dna) ;
	}
    }
  else
    { aa = arrayReCreate (aa, 12, BSunit) ;
      if (bsGetArray (obj,str2tag("S_Child"), aa, 4))
	for (i=0 ; i < arrayMax(aa) ; i += 4)
	  { if (arr(aa,i+2,BSunit).i > length)
	      length = arr(aa,i+2,BSunit).i ;
	    if (arr(aa,i+3,BSunit).i > length)
	      length = arr(aa,i+3,BSunit).i ;
	  }
      if (bsFindTag(obj,_Subsequence) && bsFlatten(obj, 3, aa))
	for (i=0 ; i < arrayMax(aa) ; i += 3)
	  { if (arr(aa,i+1,BSunit).i > length)
	      length = arr(aa,i+1,BSunit).i ;
	    if (arr(aa,i+2,BSunit).i > length)
	      length = arr(aa,i+2,BSunit).i ;
	  }
    }
  
  bsDestroy (obj) ;
  return length ;
}

Array fMapFindDNA (KEY seq, int *start, int *stop)
{
  int length ;
  int seqLen = sequenceLength(seq) ;
  Array assembly ;

  fMapInitialise() ;

  if (!*stop)
    *stop = seqLen ;
  if (!*stop)
    return 0 ;

  if (*start < *stop)
    { length = *stop - *start ;
      if (length < 10000)
	length = 10000 ;

      *start -= length ; 
      if (*start < 1)
	*start = 1 ;
      *stop += length ; 
      if (*stop > seqLen) 
	*stop = seqLen ;

      if (*stop-*start <=0) return 0 ; /* mieg a protection */
      assembly = arrayCreate (*stop-*start+2, char) ;
      array(assembly,*stop-*start,char) = 0 ;
    }
  else
    { length = *start - *stop ;
      if (length < 10000)
	length = 10000 ;

      *stop -= length ; 
      if (*stop < 1)
	*stop = 1 ;
      *start += length ; 
      if (*start > seqLen)
	*start = seqLen ;

      if (*start-*stop <=0) return 0 ; /* mieg a protection */
      assembly = arrayCreate (*start-*stop+2, char) ;
      array(assembly,*start-*stop,char) = 0 ;
    }
    
  dnaAdd (assembly, seq, *start-1, *stop-1, FALSE) ; /* allow mismatch */
  return assembly ;
}  

/********************************************/

void fMapClearDNA (LOOK look)	/*  */
{ 
  register int i ;
  char *cp, cc = WHITE ;
  
  fMapInitialise() ;
  
  if (look->colors)
    { 
      cp = arrp(look->colors, 0, char) ;
      i = arrayMax(look->colors) ;
      memset (cp, cc, i) ;
    }
}

void fMapClear (void)
{ 
  FMAPLOOKGET("fMapClear");
  
  fMapClearDNA (look) ;
  look->chosen = assReCreate(look->chosen) ;
  look->antiChosen = assReCreate(look->antiChosen) ;
  arrayDestroy(look->sites) ;
  stackDestroy(look->siteNames) ;
  fMapDraw (look, 0) ;
}

/*********************************************************/
/********************* translate genes *******************/

BOOL fMapGetCDS (LOOK look, KEY parent, Array *cds, Array *index)
{ 
  int   i, j=0, jCode = 0, cds1=0, cds2=0, tmp, dnaMax ;
  SEG   *seg ;
  BOOL isUp = FALSE ;
  Array dna ;

  if (!iskey(parent) || !look->dna)
    return FALSE ;

  dna = look->dna ;
  dnaMax = look->dna ? arrayMax(look->dna) : 0 ;
  
  if (look->segs)
    for (i = 1 ; i < arrayMax(look->segs) ; ++i)
      {
	seg = arrp(look->segs,i,SEG) ;
	if (seg->parent == parent)
	  switch (seg->type)
	    {
	    case CDS: case CDS_UP:
	      cds1 = seg->x1 ;
	      cds2 = seg->x2 ;
	      break ;
	    case SEQUENCE: 
	      isUp = FALSE ;
	      break ;
	    case SEQUENCE_UP: 
	      isUp = TRUE ;
	      break ;
	    case MPRODUCT: 
	    case MRNA: 
	      if (seg->key == seg->parent)
		{
		  cds1 = seg->x1 ;
		  cds2 = seg->x2 ;
		}
	      isUp = FALSE ;
	      break ;
	    case MPRODUCT_UP: 
	    case MRNA_UP: 
	      if (seg->key == seg->parent)
		{
		  cds1 = seg->x1 ;
		  cds2 = seg->x2 ;
		}
	      isUp = TRUE ;
	      break ;
	    default: break ;
	    }
      }
  
  if (cds1 < 0) 
    cds1 = (cds1 + 9999999)%3 ;
  if (cds1 >= cds2)
    return FALSE ;
  
  if (cds)
    *cds = arrayCreate (1000, char) ; /* wild guess */
  if (index)
    *index = arrayCreate (1000, int) ;
  
  if (look->segs)
    for (i = 1 ; i < arrayMax(look->segs) ; ++i)
      { 
	seg = arrp(look->segs,i,SEG) ;
	if (seg->parent == parent)
	  switch (seg->type)
	    {			/* pick cds */
	    case MRNA: case MRNA_UP:
	    case MPRODUCT: case MPRODUCT_UP:
	      if (seg->key && !strstr(name(seg->key), "Exon") &&  !strstr(name(seg->key), "Gap"))
		break ;
	      /* else fall thru */
	    case EXON: case EXON_UP:
	      for (j = seg->x1 ; j <= seg->x2 ; ++j)
		if (j >= cds1 && j <= cds2)
		  { if (cds)
		      {
			if (j >= 0 && j < dnaMax)
			  array(*cds,jCode,char) = arr(dna,j,char) ;
			else
			  array(*cds,jCode,char) = 0 ;
		      }
		    if (index)
		      array(*index,jCode,int) = j ;
		    ++jCode ;
		  }
	      break ;
	    default: break ;
	    }
      }
  /* zero terminate */
  if (cds && *cds)
    {
      int max = arrayMax (*cds) ;
      array (*cds, max, char) = 0 ;
      arrayMax (*cds) = max ;
    }
  
  if (isUp)
    { 
      if (cds)
	reverseComplement (*cds) ;
      if (index && *index)
	for (j = 0, i = jCode-1 ; j < i ; ++j, --i)
	  { 
	    tmp = arr(*index,j,int) ;
	    arr(*index,j,int) = arr(*index,i,int) ;
	    arr(*index,i,int) = tmp ;
	  }
    }
  return cds && *cds && arrayMax (*cds) ? TRUE : FALSE ; ;
}

/********************/

void fMapColorIntronsExons (int box)
{
  int i, j1, j2, j3, colour, length ;
  Array cds, index, colors ;
  SEG *seg ;
  char *translationTable ;
  FMAPLOOKGET("fMapColorIntronsExons");

  length = look->length ;

  if (box >= arrayMax(look->boxIndex) ||
      !(seg = BOXSEG(box)) ||
      !fMapGetCDS (look, seg->parent, &cds, &index))
    { messerror ("ColorIntronsExons only works on coding genes") ;
      return ;
    }

  if (arrayMax(cds) % 3)
    messout ("Coding length is %d mod 3", arrayMax(cds) % 3) ;

  fMapClearDNA (look) ;
  colors = look->colors ;

  j1 = arr(index,0,int) ;
  j2 = arr(index,arrayMax(index)-1,int) ;
  if (j2 < j1) { j3 = j1 ; j1 = j2 ; j2 = j3 ; }
  if (j1 < 0) j1 = 0 ; 
  if (j2 >= length) j2 = length-1 ;
  /* do not color the introns they print same as exons
    if CYAN and LIGHTGREEN
     for (i = j1 ; i <= j2 ; ++i)
    arr(colors,i,char) |= TINT_LIGHTGREEN ;
    */ 
  translationTable = pepGetTranslationTable (look->seqKey, 0) ;
  for (i = 0 ; i + 2 < arrayMax(index) ; i += 3)
    { j1 = arr(index,i,int) ;
      j2 = arr(index,i+1,int) ;
      j3 = arr(index,i+2,int) ;
      if (j1 < 0 || j3 < 0 || j1 >= length || j3 >= length)
	continue ;
      if (j3 == j1+2 || j1 == j3+2)
	{
	  if (e_codon(arrp(cds,i,char), translationTable) == '*')
	    colour = TINT_RED ;
	  else
	    colour = TINT_YELLOW ; /* was  TINT_CYAN ; */
	}
      else
	colour = TINT_MAGENTA ;
      arr(colors,j1,char) |= colour ; 
      arr(colors,j2,char) |= colour ; 
      arr(colors,j3,char) |= colour ; 
    }

  arrayDestroy (cds) ;
  arrayDestroy (index) ;

  mapColSetByName ("DNA Sequence", TRUE) ;
  look->activeBox = 0 ;
  fMapDraw (look, 0) ;
}

/***************************************/

static void translationTitle (FILE *fil, KEY seq)
{
  OBJ obj = bsCreate (seq) ;
  KEY key ;

  fprintf (fil, ">%s", name(seq)) ;
  if (obj)
    { if (bsGetKey (obj, _Locus, &key))
	fprintf (fil, " %s:", name(key)) ;
      if (bsGetKey (obj, _Brief_identification, &key))
	fprintf (fil, " %s", name(key)) ;
      if (bsFindTag (obj, _Start_not_found))
	fprintf (fil, " (start missing)") ;
      if (bsFindTag (obj, _End_not_found))
	fprintf (fil, " (end missing)") ;
      bsDestroy (obj) ;
    }
  else
    messerror ("Missing object for %s in translationTitle", name(seq)) ;
  fprintf (fil, "\n") ;
}

static char fname[FIL_BUFFER_SIZE], dname[DIR_BUFFER_SIZE] ;

void fMapExportTranslation (int box)
{
  FILE *fil ;
  Array cds ;
  int i ;
  SEG *seg ;
  char *translationTable ;
  FMAPLOOKGET("fMapExportTranslation") ;
    
  translationTable = pepGetTranslationTable (look->seqKey, 0) ;
  if (box >= arrayMax(look->boxIndex) ||
      !(seg = BOXSEG(box)) ||
      !fMapGetCDS (look, seg->parent, &cds, 0))
    { messerror ("ExportTranslation only works on coding genes") ;
      return ;
    }

  if (!(fil = filqueryopen (dname, fname, "pep", "a",
			    "File to add translation to")))
    { arrayDestroy (cds) ;
      return ;
    }
  translationTitle (fil, seg->parent) ;
  for (i = 0 ; i < arrayMax(cds) ; i+=3)
    { fputc (e_codon(arrp(cds,i,char), translationTable), fil) ;
      if ((i+3) % 150 == 0)
	fputc ('\n', fil) ;
    }
  if (i % 150)
    fputc('\n',fil) ;
  filclose (fil) ;
  arrayDestroy (cds) ;
}

void fMapExportcDNA (int box)
{
  FILE *fil ;
  Array cds = 0 ;
  int i ;
  SEG *seg ;
  FMAPLOOKGET ("fMapExportcDNA") ;
    
  if (box >= arrayMax(look->boxIndex) ||
      !(seg = BOXSEG(box)) ||
      !fMapGetCDS (look, seg->parent, &cds, 0))
    { messerror ("ExportcDNA only works on coding genes") ;
      return ;
    }

  if (!(fil = filqueryopen (dname, fname, "cdna", "a",
			    "File to add cDNA to")))
    { arrayDestroy (cds) ;
      return ;
    }
  translationTitle (fil, seg->parent) ;
  for (i = 0 ; i < arrayMax(cds) ; ++i)
    { fputc (dnaDecodeChar[(int)arr(cds,i,char)], fil) ;
      if (((i+1) % 50) == 0)
	fputc ('\n', fil) ;
    }
  if (i % 50)
    fputc('\n',fil) ;
  filclose (fil) ;
  arrayDestroy (cds) ;
}

/*********************/

void fMapExportTranslations (void)
{
  FILE *fil ;
  Array cds = 0 ;
  KEY seq ;
  SEG *seg ;
  int i, j ;
  char *translationTable ;
  FMAPLOOKGET("fMapExportTranslations") ;

  if (!(fil = filqueryopen (dname, fname, "pep", "a",
			    "File to add translation to")))
    { arrayDestroy (cds) ;
      return ;
    }
  translationTable = pepGetTranslationTable (look->seqKey, 0) ;
  for (j = 1 ; j < arrayMax(look->segs) ; ++j)
    { seg = arrp(look->segs, j, SEG) ;
      if ((seg->type == CDS || seg->type == CDS_UP) &&
	  seg->x1 >= look->zoneMin && seg->x2 <= look->zoneMax)
	{ seq = seg->parent ;
	  if (!fMapGetCDS (look, seq, &cds, 0))
	    { messout ("screwup finding CDS %s", name(seq)) ;
	      continue ;
	    }
	  translationTitle (fil, seq) ;
	  for (i = 0 ; i < arrayMax(cds) ; i+=3)
	    { fputc (e_codon(arrp(cds,i,char), translationTable), fil) ;
	      if ((i+3) % 150 == 0)
		fputc ('\n', fil) ;
	    }
	  if (i % 150)
	    fputc('\n',fil) ;
	}
    }
  filclose (fil) ;
  arrayDestroy (cds) ;
}

/* mhmp 18.06.98 */
void fMapShowUpGeneTranslation (LOOK look, float *offset)
{
  Array cds, index ;
  if (fMapGetCDS (look, look->translateGene, &cds, &index))
    fMapShowGeneTranslation ( look,  offset) ;
  else
    { messout ("Select an up gene to translate with the "
	       "pulldown menu on the genes") ;
      mapColSetByName ("Up Gene Translation", FALSE) ; 
      look->pleaseRecompute = TRUE ; /* so the above FALSE takes effect */
    }
}

/* mhmp 18.06.98 */
void fMapShowDownGeneTranslation (LOOK look, float *offset)
{
  Array cds, index ;
  if (fMapGetCDS (look, look->translateGene, &cds, &index))
    fMapShowGeneTranslation ( look,  offset) ;
  else
    { messout ("Select a down gene to translate with the "
	       "pulldown menu on the genes") ;
     mapColSetByName ("Down Gene Translation", FALSE) ; 
     look->pleaseRecompute = TRUE ; /* so the above FALSE takes effect */
    }
}
/*********************/

void fMapShowGeneTranslation (LOOK look, float *offset)
{ 
  int i, j, window, min, max ;
  char c ;
  Array cds, index ;
  int k = 1;
  char *translationTable ;

  fMapInitialise() ;
  translationTable = pepGetTranslationTable (look->seqKey, 0) ;
  if (fMapGetCDS (look, look->translateGene, &cds, &index))
    { window = look->max - look->min ;
      geneTranslation = (char*) messalloc (window) ;
      geneTransPos = (int*) messalloc (window * sizeof(int)) ;
      min = arr(index,0,int) - look->min ; 
      max = arr(index,arrayMax(index)-1,int) - look->min ;
      upTranslation = FALSE ;
      if (min > max) { upTranslation = TRUE ; i = min ; min = max ; max = i ; }
      if (min < 0) min = 0 ;
      if (max > window) max = window ;
      for (j = min ; j < max ; ++j)
	{
	  geneTranslation[j] = '-' ;
	  geneTransPos[j] = 0 ;
	}
      for (i = 0 ; i < arrayMax(index) ; i+=3, k++)
	{ c = e_codon(arrp(cds,i,char), translationTable) ;
	  j = arr(index,i,int) - look->min ;
	  if (j >= 0 && j+2 < window)
	    { 
	      geneTranslation[j] = c ;
	      geneTransPos[j] = k ;
	      geneTransPos[j+1] = k ;
	      geneTransPos[j+2] = k ;
	      geneTranslation[j+1] = '%' ; /* impossible letter value */
	      geneTranslation[j+2] = '%' ;
	    }
	}
      showSequence (look, *offset, dnaWidth, TRUE) ; 
      messfree (geneTranslation) ;
      messfree (geneTransPos) ;
      arrayDestroy (cds) ;
      arrayDestroy (index) ;
      *offset += look->dnaWidth/3 + 1 + 8;/* mhmp 12.06.98 */
    }
}

/********************/

void fMapFindCoding (LOOK look)
{
  int i, j, j1 ;
  Array index ;
  SEG *seg ;
  KEY parent ;

  for (i = 1 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (seg->type == CDS && (parent = seg->parent) &&
	  fMapGetCDS (look, parent, 0, &index))
	{ j1 = 0 ;
	  for (j = 0 ; j < arrayMax(index) ; ++j)
	    if (j == arrayMax(index) - 1 ||
		arr(index,j+1,int) > arr(index,j,int)+1) /* boundary */
	      { seg = arrayp (look->segs, arrayMax(look->segs), SEG) ;
		seg->key = _CDS ;
		seg->type = CODING ;
		seg->x1 = arr(index, j1, int) ; 
		seg->x2 = arr(index, j, int) ;
		seg->parent = parent ;
		seg->data.i = j1 ;
		j1 = j+1 ;
	      }
	  arrayDestroy (index) ;
	}
      else if (seg->type == CDS_UP && (parent = seg->parent) &&
	       fMapGetCDS (look, parent, 0, &index))
	{ j1 = 0 ;
	  for (j = 0 ; j < arrayMax(index) ; ++j)
	    if (j == arrayMax(index) - 1 ||
		arr(index,j+1,int) < arr(index,j,int)-1) /* boundary */
	      { seg = arrayp (look->segs, 
			      arrayMax(look->segs), SEG) ;
		seg->key = _CDS ;
		seg->type = CODING_UP ;
		seg->x1 = arr(index,j,int) ; 
		seg->x2 = arr(index,j1,int) ;
		seg->parent = parent ;
		seg->data.i = j1 ;
		j1 = j+1 ;
	      }
	  arrayDestroy (index) ;
	}
				/* remove duplicate INTRON segs */
      if (seg->type == INTRON || seg->type == INTRON_UP)
	{ SEG *seg2 = seg ;
	  int flags = 0 ;
	  BOOL existsParent = FALSE ;
	  SEQINFO *sinf ;

	  while (i + 1 < arrayMax(look->segs) && seg->type == seg2->type &&
		 seg->x1 == seg2->x1 && seg->x2 == seg2->x2)
	    { sinf = arrp(look->seqInfo, seg->data.i, SEQINFO) ;
	      flags |= sinf->flags ;
	      if (seg->parent)
		existsParent = TRUE ;
	      ++i ; ++seg ;
	    }
	  for (--i, --seg ; seg2 <= seg ; ++seg2)
	    if (seg2->parent || !existsParent)
	      { sinf = arrp(look->seqInfo, seg2->data.i, SEQINFO) ;
		if (sinf->flags != flags)
		  { sinf = arrayp(look->seqInfo, arrayMax(look->seqInfo), SEQINFO) ;
		    *sinf = arr(look->seqInfo, seg2->data.i, SEQINFO) ;
		    sinf->flags = flags ;
		    seg2->data.i = arrayMax(look->seqInfo)-1 ;
		  }
		existsParent = TRUE ;
	      }
	    else
	      { *seg2 = array(look->segs,arrayMax(look->segs)-1,SEG) ;
		--arrayMax(look->segs) ;
	      }
	}
    }
}

/************************************************/ 

void fMapShowCoding (LOOK look, float *offset) /* width 0.5 */
{
  float y1, y2 ;
  int i ;
  SEG *seg ;
  int frame = look->min % 3 ;

  if (!look->map->isFrame) 
    return ;

  graphColor (BLUE) ;

  for (i = 1 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (seg->type != CODING ||
	  (seg->x1 - seg->data.i)%3 != frame ||
	  seg->x1 > look->max || 
	  seg->x2 < look->min)
	continue ;
      y1 = MAP2GRAPH(look->map,seg->x1) ;
      y2 = MAP2GRAPH(look->map,seg->x2+1) ;	/* to cover full base */
      if (y1 < 2 + topMargin) y1 = 2 + topMargin ;
      if (y2 < 2 + topMargin) y2 = 2 + topMargin ;
      
      array (look->boxIndex, graphBoxStart(), int) = i ;
      graphRectangle (*offset, y1, *offset+0.4, y2) ;
      graphBoxEnd () ;
    }

  graphColor (BLACK) ;
  *offset += 0.5 ;
}

/************************************************************/

void fMapShowORF (LOOK look, float *offset) /* width 2 */
{
  char *cp ;
  float y1 = 0, y2 ;
  int i1, i, box ;
  int frame ;
  char *translationTable ;
  SEG *seg ;

  fMapInitialise() ;
  translationTable = pepGetTranslationTable (look->seqKey, 0) ;
  frame = look->min % 3 ;

  cp = arrp (look->dna, look->min, char) ;
  for (i1 = look->min ; i1 >= 0 ; i1-=3, cp-=3)
    if (e_codon(cp, translationTable) == '*')
      break ;
  i1 += 3 ;
  y1 = MAP2GRAPH(look->map,look->min) ;

  cp = arrp (look->dna, look->min+3, char) ;
  for (i = look->min+3 ; i1 < look->max ||
       			 i + 3 < arrayMax(look->dna) ; i+=3, cp+=3)
    if (i + 3 >= arrayMax(look->dna) ||
	e_codon(cp, translationTable) == '*')
      { y2 = MAP2GRAPH(look->map, i) ;
	if (i > i1)
	  { box = graphBoxStart() ;
	    if (i1 > look->min)
	      graphLine (*offset+0.25, y1, *offset+1.75, y1) ;
	    else
	      { graphLine (*offset+0.25, y1, *offset+0.5, y1) ;
		graphLine (*offset+1.5, y1, *offset+1.75, y1) ;
	      }
	    if (i < look->max)
	      graphLine (*offset+0.25, y2, *offset+1.75, y2) ;
	    else
	      { y2 = MAP2GRAPH(look->map, look->max) ;
		graphLine (*offset+0.25, y2, *offset+0.5, y2) ;
		graphLine (*offset+1.5, y2, *offset+1.75, y2) ;
	      }
	    graphBoxEnd () ;
	    array(look->boxIndex, box, int) = arrayMax(look->segs) ;

	    seg = arrayp(look->segs, arrayMax(look->segs), SEG) ;
	    seg->key = _ORF ;
	    seg->type = ORF ;
	    seg->x1 = i1 ; seg->x2 = i-1 ;
	    seg->data.i = frame ;
	    seg->parent = arr(look->boxIndex, box, int) ;

	    if (assFind (look->chosen, SEG_HASH(seg), 0))
	      graphBoxDraw (box, BLACK, GREEN) ;
	    else if (assFind (look->antiChosen, SEG_HASH(seg), 0))
	      graphBoxDraw (box, BLACK, LIGHTGREEN) ;
	    else
	      graphBoxDraw (box, BLACK, WHITE) ;
	    graphBoxFreeMenu (box, fMapChooseMenuFunc, 
			      fMapChooseMenu) ;
	  }
	else
	  graphLine (*offset+0.25, y2, *offset+1.75, y2) ;
	y1 = y2 ;
	i1 = i+3 ;
	if (i1 > look->max)
	  break ;
      }
  *offset += 2 ;
}

/******************************************/


void fMapShowATG (LOOK look, float *offset) /* width 1 */
{
  char *cp ;
  int i, box, iseg ;
  int atg, atgMax ;
  float y, width, height ;
  SEG *seg, *atgSeg ;
  static KEY _ATG_ = 0 ;
  char *translationTable ;

  if (!_ATG_)
    lexaddkey ("ATG", &_ATG_, _VMethod) ;
  methodSet ("ATG", YELLOW, 
	     METHOD_STRAND_SENSITIVE | METHOD_FRAME_SENSITIVE,
	     4.195, 1, 'A', 0, 0) ;

  *offset += 0.5 ;		/* centre line */
 
  fMapFindSegBounds (look, ATG, &atg, &atgMax) ;
  translationTable = pepGetTranslationTable (look->seqKey, 0) ;

  cp = arrp (look->dna, look->min, char) ;
  for (i = look->min ; i < look->max &&
		       i + 3 < arrayMax(look->dna) ; i+=3, cp+=3)
    { 
      while (atg < atgMax && 
	     (atgSeg = arrp(look->segs, atg, SEG))->x1 < i)
	++atg;
      if (atg < atgMax && atgSeg->x1 == i)
	{ iseg = atg ;
	  seg = atgSeg ;
	}
      else if (e_codon(cp, translationTable) == 'M')
	{ iseg = arrayMax(look->segs) ;
	  seg = arrayp(look->segs, arrayMax(look->segs), SEG) ;
	  seg->key = _ATG_ ;
	  seg->type = ATG ;
	  seg->x1 = i ; seg->x2 = i+2 ;
	  seg->data.f = 0.0 ;
	}
      else 
	continue;

      width = 0.5 * fMapSquash (seg->data.f, 1.0) ;
      box = graphBoxStart();
      height = look->map->mag * 3 ;
      if (height < 0.3) 
	height = 0.3; ;
      if (width < 0.2) 
	width = 0.2 ;
      y = MAP2GRAPH(look->map, i) ;
      graphRectangle (*offset-width, y, 
		      *offset+width, y+height) ;
      graphBoxEnd () ;
      if (assFind (look->chosen, SEG_HASH(seg), 0))
	graphBoxDraw (box, BLACK, GREEN) ;
      else if (assFind (look->antiChosen, SEG_HASH(seg), 0))
	graphBoxDraw (box, BLACK, LIGHTGREEN) ;
      else
	graphBoxDraw (box, BLACK, YELLOW) ;
      graphBoxFreeMenu (box, fMapChooseMenuFunc, fMapChooseMenu) ;
      array(look->boxIndex, box, int) = iseg ;
  
    }
  *offset += 0.5 ;
}

/************************************************/
/************************************************/
