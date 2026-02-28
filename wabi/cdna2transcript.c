/*  File: cdna2transcript.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg. 1997
 *-------------------------------------------------------------------
 * This file is part of the ACEMBLY extension to 
        the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
    Construct Transcrips out of the alignments of cdna
 * Exported functions:
 * HISTORY:
 * Created: Aug 2000 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */
/* #define ARRAY_CHECK */
/* #define MALLOC_CHECK */
/* #define CHRONO */


#include "acedb.h"
#include "cdna.h"
#include "chrono.h"

/*********************************************************************/
/*********************************************************************/
/* for each clone, */
static Array c2tSplitGenes (int *j0p, Array genes)
{
  KEY gene = 0 ;
  HIT *gh ;
  int j1, gg = 0;
  Array oneGene = 0 ;

  for (j1 = *j0p ; j1 < arrayMax(genes) ; j1++)
    {
      gh = arrp (genes, j1, HIT) ; 
      if (gGene & gh->type) /* new gene */
	{
	  if (gene)
	    break ;
	  gene = gh->gene ;
	  oneGene = arrayCreate (200, HIT) ;
	  gg = 0 ;
	}
      else if (gDroppedGene & gh->type) /* dropped gene */
	{
	  if (gene)
	    break ;
	  gene = 0 ;
	}
      if (gene && oneGene)
	{
	  array (oneGene, gg++, HIT) = arr (genes, j1, HIT) ;
	}
    }
  *j0p = j1 ;
  if (!(arrayMax(oneGene) < 2))
    arrayDestroy (oneGene) ;
  return oneGene ;
}
  
/*********************************************************************/

static int c2tOrderByA1 (const void *va, const void *vb)
{
  HIT *up = (HIT *)va, *vp = (HIT *)vb ;

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

static Array c2tSegmentize (Array oneGene)
{
  HIT *gh, *sh ; 
  KEY gene = 0 ;
  Array segs = arrayCreate (256, HIT) ;
  int j1, jSeg = 0, a2max = -1 ;

  gh = arrp (oneGene, 0, HIT) ; 
  gene = gh->gene ;

  for (j1 = 1 ; j1 < arrayMax(oneGene) ; j1++)
    {
      gh = arrp (oneGene, j1, HIT) ; 
      sh = arrayp (segs, jSeg++, HIT) ;
      sh->gene = gene ;
      sh->a1 = gh->a1 ;
      sh->a2 = sh->a1 + 1 ;
      if (a2max < gh->a2)
	a2max = gh->a2 ;
    }
  arraySort (segs, c2tOrderByA1) ;
  arrayCompress (segs) ; 
  for (jSeg = 0 ; jSeg < arrayMax(segs) - 1 ; jSeg++)
    {
      sh = arrayp (segs, jSeg, HIT) ;
      sh->a2 = (sh + 1)->a1 - 1 ;
    }
  sh = arrayp (segs, jSeg, HIT) ;
  sh->a2 = a2max ;
  return segs ; 
}

/*********************************************************************/

void cdna2transcripts (KEY cosmid, KEY cosmid1, 
		       Array dna, Array dnaR,
		       Array genes, Array linkPos,
		       Array geneHits)
{
  int jGene = 0 ;
  Array oneGene = 0 ;
  Array segments = 0 ;

  while ((oneGene = c2tSplitGenes (&jGene, genes))) /* run one gene at a time */
    {
      segments = c2tSegmentize (oneGene) ; /* cut all alternatives into a tiling path of segments */

      showHits(segments) ;
      arrayDestroy (segments) ;
      arrayDestroy (oneGene) ;
    }
}

/*********************************************************************/
/*********************************************************************/
