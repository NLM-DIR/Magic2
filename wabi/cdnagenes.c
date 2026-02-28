/*  File: cdnaalign.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg. 1997
 *-------------------------------------------------------------------
 * This file is part of the ACEMBLY extension to 
        the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
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
#define CHRONO


#include "freeout.h"
#include "lex.h"
#include "bs.h"
#include "a.h"
#include "dna.h"
#include "systags.h"
#include "tags.h"
#include "classes.h"
#include "dna.h"
#include "dnaalign.h"
#include "bump.h"
#include "query.h"
#include "session.h"
#include "bindex.h"
#include "cdna.h"
#include "acembly.h"
#include "mytime.h"
#include "cdnainit.h"
#include "parse.h"
#include "chrono.h"

/***********************************************************/
/***********************************************************/
/***********************************************************/
/* get all coords of transcribed genes on chrom */
static Array tg2chrom (KEY chrom)
{
  KEYSET tgs = 0 ;
  Array hits = 0 ;
  int i, a1, a2, j ;
  KEY tg ;
  OBJ Tg = 0 ;
  HIT *hh ;

  tgs = queryKey (chrom, "> Transcribed_gene") ;
  hits = arrayCreate (10000, HIT) ;
  for (i = 0, j= 0 ; i < arrayMax (tgs) ; i++)
    {
      tg = keySet (tgs, i) ;
      if ((Tg = bsCreate (tg)))
	{
	  if (bsGetData (Tg, _IntMap, _Int, &a1) &&
	      bsGetData (Tg, _bsRight, _Int, &a2))
	    {
	      hh = arrayp (hist, j++, HIT) ;
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
  int max = arrayMax (dna) - 2, i ;
  char *cp ;
  Array stops = arrayCreate (2000000, int) ;
  int j ;

  for (i = frame, cp = arrp (dna, i, char) ; i < max ; i++, cp++)
    {
      if (*cp == T_ && 
	  (
	   (*(cp+1)== A_ && *(cp+2) == A_) ||
	   (*(cp+1)== A_ && *(cp+2) == G_) ||
	   (*(cp+1)== G_ && *(cp+2) == A_) 
	   )
	  )
	array (stops, j++, int) = i ;
    }  
}

/***********************************************************/

static Array getStopHisto (Array histo, Array stops)
{
  int i, *s, s0, dx, *h ;

  if (!arrayExists(histo))
    histo = arrayCreate (30000, int) ;

  if (stops && arrayMax(stops))
    for (i = 0, s0 = 0, s = arrp (stops, 0, int) ;
	 i < arrayMax(stops) ; i++, s++)
      { 
	dx = (*s - s0 - 3)/3 ;
	h = arrayp (histo, dx) ;
	*h++ ;
      }
  return histo ;
}

/***********************************************************/

static Array chromHisto (KEY chrom, Array histo)
{
  Array dna, stops ;
  int frame ;

  if (!arrayExists(histo))
    histo = arrayCreate (30000, int) ;

  dna = dnaGet (chrom) ;
  if (dna && arrayMax(dna))
    {
      for (frame = 0 ; frame < 3 ; frame++)
	{
	  stops = getStops (dna, frame) ;
	  getStopHisto (histo, stops) ;
	  arrayDestroy (stops) ;
	}
      reverseComplement (dna) ;
      for (frame = 0 ; frame < 3 ; frame++)
	{
	  stops = getStops (dna, frame) ;
	  getStopHisto (histo, stops) ;
	  arrayDestroy (stops) ;
	}
    }
  arrayDestroy (dna) ;
  return histo ;
}

/***********************************************************/

void cdnaStopHisto (void)
{
  KEYSET chroms = query (0, "Find sequence genomic ; !Overlap_left") ;
  KEY chrom, cosmid ;
  int i ;

  for (i = 0 ; i < keySetMax (ks) ; i++)
    {
      cosmid = keySet (chroms, i) ;
      while (( chrom = keyGetKey (cosmid, _Source)))
	{ cosmid = chrom ; }
      printf("i=%d\t%s\n", i, name(cosmid)) ;
      histo = chromHisto (cosmid, histo) ;
    }
  plotHisto ("ORF lengths", histo) ;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
