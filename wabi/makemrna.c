/*  File: makemrna.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
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
 * Created: Mars 2001 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */
#define CHRONO   amounts to 1 minute/ 2 hours for whole chromo 1 worm

#include "ac.h"
#include "bitset.h"
#include "cdna.h"
#include "cdnatr.h"
#include "makemrna.h"
#include "bql.h"
#include "lex.h"
#include "dna.h"
#include "query.h"
#include "parse.h"
#include "pick.h"
#include "../whooks/systags.h"
#include "../whooks/tags.h"
#include "../whooks/classes.h"
#include "peptide.h"
#include "session.h"
#include "dnaalign.h"
#include "cdnatr.h"
#include "giw.h"
#include "../wfiche/gtitle.h"
#include "chrono.h"
#include "../wfiche/biolog.h"


extern Array getSeqDna (KEY key) ;
static Array mrnaRefseqMakerDna2Pep (KEY cosmid, Array dna, BOOL isMrna5pComplete) ;
static int mrnaAnalysePg2Tg (KEY g1, int a1, int a2, BOOL isUp1, KEY g2, int b1, int b2, BOOL isUp2) ;
static void mrnaAddFirstMetAndUorf (KEY mrna) ;
static BOOL smrnaDnaDestroy (SMRNA *smrna) ;
static void mrnaFixNonBestProductTitle (KEY mrna) ;

static void showOrfs (Array orfs)
{
  int ii = 0 ;
  ORFS *orf ;
  
  if (orfs) 
    {
      printf ("\n") ;

      for (ii = 0 ; ii < arrayMax(orfs) ; ii++)
	{
	  orf = arrp (orfs, ii, ORFS) ;
	  printf ("orf %d:: d=%2d  f=%1d upStop=%3d met=%3d leu=%3d start=%3d max=%3d downStop=%3d nOpen=%3d cds=%3d nI=%d nX=%3d %s %s\n",
		  ii, orf->iDna,orf->frame,orf->upStop,orf->met,orf->leu,orf->start,orf->max,orf->downStop,orf->nOpen,orf->cds,orf->nIntron, orf->nX,orf->isSL?"SL ":"", orf->isUorf ? "uORF" : "") ;
	  
	}
    }
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
/************** Strategy, imported from clone MainClone **************/

static int getPleaseNewNames (void)
{
  static int style = -1 ;

  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && NewNames") ;
      style = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    style = 0 ; /* default */

  return style ;
}

/*********************************************************************/

static int getPleaseStealUpStream (void)
{
  static int style = -1 ;

  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && StealUpStream") ;
      style = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    style = 0 ; /* default */

  return style ;
}

/*********************************************************************/

static int getPleaseIgnorePolyAInProduct (void)
{
  static int style = -1 ;

  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && IgnorePolyAInProduct") ;
      style = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    style = 0 ; /* default */

  return style ;
}

/*********************************************************************/

static int getPleaseMarkDoubleFuzzy (void)
{
  static int style = -1 ;

  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && MarkDoubleFuzzy") ;
      style = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    style = 0 ; /* default */

  return style ;
}

/*********************************************************************/

static int getPleaseNoKantorInfo (void)
{
  static int style = -1 ;

  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && NoKantorInfo") ;
      style = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    style = 0 ; /* default */

  return style ;
}

/*********************************************************************/

int mrnaPleaseNameBy (void)
{
  static int style = -1 ;

  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && Name_by_section") ;
      style = keySetMax (ks) ? 1 : -1 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && Name_by_gene") ;
      style = keySetMax (ks) ? 2 : -1 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    style = 0 ; /* default */

  return style ;
}

/*********************************************************************/

static int mrnaPleaseMinTranscriptSize (void)
{
  static int minProteinLength = -1 ;

  if (minProteinLength == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && MinTranscriptSize") ;
      minProteinLength = keySetMax (ks) ? 1 : -1 ;
      
      if (minProteinLength == 1)
        {
          OBJ Clone = bsCreate (keySet (ks,0)) ;
          
          if (Clone)
            {
              bsGetData (Clone, str2tag("MinTranscriptSize"), _Int, &minProteinLength) ;
              bsDestroy (Clone) ;
            }
          if (minProteinLength == 1)
            minProteinLength = 80 ;
        }
      keySetDestroy (ks) ;
    }
  if (minProteinLength == -1)
    minProteinLength = 0 ; /* default */

  return minProteinLength ;
}  /* mrnaPleaseMinTranscriptSize */

/*********************************************************************/
/* maximal gap size in kb */
static int mrnaPleaseBigGapSize (void)
{
  static int maxGapSize = -1 ;

  if (maxGapSize == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && MaxGapSize") ;
      maxGapSize = keySetMax (ks) ? -2 : -1 ;
      
      if (maxGapSize == -2)
        {
          OBJ Clone = bsCreate (keySet (ks,0)) ;
          
          if (Clone)
            {
              bsGetData (Clone, str2tag("MaxGapSize"), _Int, &maxGapSize) ;
	      maxGapSize *= 1000 ;
              bsDestroy (Clone) ;
            }
        }
      keySetDestroy (ks) ;
      if (maxGapSize < 0)
	maxGapSize = 0 ; /* default */
    }

  return maxGapSize ;
}

/*********************************************************************/

static int mrnaPleaseStealFromPrediction (void)
{
  static int style = -1 ;

  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && StealFromPrediction") ;
      style = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    style = 0 ; /* default */

  return style ;
}

/*********************************************************************/
#ifdef JUNK
static int mrnaPleaseSplitByGeneId (void)
{
  static int style = -1 ;

  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && SplitByGeneId") ;
      style = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    style = 0 ; /* default */

  return style ;
} /* mrnaPleaseSplitByGeneId */
#endif
/*********************************************************************/
#ifdef JUNK
static int mrnaPleaseFuseGaps (void)
{
  static int style = -1 ;

  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && species = worm") ;
      style = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    style = 0 ; /* default */

  return style ;
}
#endif
/*********************************************************************/

int mrnaPleaseNoAbandon (void)
{
  static int style = -1 ;

  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && NoAbandon") ;
      style = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    style = 0 ; /* default */

  return style ;
}

/*********************************************************************/

int mrnaPleaseKillNonBest (void)
{
  static int style = -1 ;

  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && Kill_non_best_mRNA") ;
      style = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    style = 0 ; /* default */

  return style ;
}

/*********************************************************************/
/* notice the inversion function wants single, the tag is multiple */
static int mrnaPleaseSingleORF (void)
{
  static int style = -1 ;

  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && MultipleORF") ;
      style = keySetMax (ks) ? 0 : 1 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    style = 1 ; /* default */

  return style ;
}

/*********************************************************************/
/* notice the inversion function wants single, the tag is multiple */
static int mrnaPleaseUseLeucine (void)
{
  static int style = -1 ;

  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && UseLeucine") ;
      style = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    style = 0 ; /* default */

  return style ;
}

/*********************************************************************/

static void mrnaCountErrors (S2M *s2m, SC *sc, Array hits) 
{
  HIT *up ;
  A_ERR *ep ;
  int i, ii, a1, a2, x1, x2, x, maxExact ;
  Array dnaEst ;
  static Array err = 0 ;

  messcrash ("mrnaCountErrors is utterly wrong on reverse strand and bizare on direct starnd, see 3b333") ;
  err = arrayReCreate (err, 12, A_ERR) ;
  if (arrayMax(hits))
    for (ii = 0, up = arrp (hits, 0, HIT) ; ii < arrayMax(hits) ; ii++, up++)
      {
        if ((up->type & gI) ||
            up->est == 1 ||
            up->a1 > up->a2 - 5 ||
            ! (dnaEst = getSeqDna (up->est)))
          { up->nerr = up->nerrAll = 0 ; continue ; }
        x1 = up->x1 ; x2 = up->x2 ;
        a1 = up->a1 ; a2 = up->a2 ;

        if (sc->a1 < sc->a2) 
          { a1 = sc->a1 + a1 - 1 ; a2 = sc->a1 + a2 - 1 ; }
        else
          { a1 = sc->a1 - a2 + 1 ; a2 = sc->a1 - a1 + 1 ;
            x = x1 ; x1 = x2 ; x2 = x ;
          }
        err = aceDnaDoubleTrackErrors (dnaEst, &x1, &x2, x1 < x2,
                                  s2m->dnaD, s2m->dnaR, &a1, &a2,
                                 0, err, 8, 0, FALSE, &maxExact) ;
	up->maxExact = maxExact ;
        up->nerrAll = arrayMax(err) ;
        up->nerr = 0 ;
        if (arrayMax (err))
          for (i = 0, ep = arrp (err, 0, A_ERR) ; i < arrayMax(err) ; i++, ep++)
            if (ep->iShort < 600) up->nerr++ ;
      }
} /* mrnaCountErrors */

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

static Array  cDnaGetIntMap (KEYSET ks, AC_HANDLE h)
{
  int j = 0, nn = keySetMax (ks), a1, a2 ;
  Array hits = arrayHandleCreate (nn, HIT, h) ;
  HIT *hh ;
  OBJ Obj = 0 ;
  KEY key, map ;

  while (nn--)
    {
      key = keySet (ks, nn) ;
      if ((Obj = bsCreate (key)))
        {
          a1 = a2 = 0 ; map = 0 ;
          if (bsGetKey (Obj, _IntMap, &map) &&
              bsGetData (Obj, _bsRight, _Int, &a1) &&
              bsGetData (Obj, _bsRight, _Int, &a2))
            {
              hh = arrayp (hits, j++, HIT) ;
              hh->gene = key ; hh->est = map ;
              hh->a1 = a1 ; hh->a2 = a2 ;
              hh->reverse = a1 > a2 ? TRUE : FALSE ;
            }
          bsDestroy (Obj) ;
        }      
    }
   
  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ; /* cdnaclone then est then a1 then x1 */
  return (hits) ;
} /* cDnaGetIntMap */

/*********************************************************************/

static Array cDnaGetAllIntMap (KEY tag, AC_HANDLE h)
{
  Array hits = arrayHandleCreate (10000, HIT, h) ;
  KEYSET maps = query (0, "Find map") ;
  Array units = arrayCreate (100000, BSunit) ;
  OBJ Map = 0 ;
  KEY map ;
  HIT *hh ;
  int a1, a2, j = 0, ir, nn = keySetMax (maps) ;
  BSunit *uu ;

  while (nn--)
    {
      map = keySet (maps, nn) ;
      if ((Map = bsCreate (map)))
        {
	  if (bsGetArray (Map, tag, units, 4))
	    for (ir = 0 ; ir < arrayMax (units) ; ir += 4)
	      {
		uu = arrp (units, ir, BSunit) ;
		hh = arrayp (hits, j++, HIT) ;
		hh->gene = uu[0].k ; hh->est = map ;
		a1 = uu[1].i ; a2 = uu[2].i ;
		if (a1 < a2)
		  { hh->a1 = a1 ; hh->a2 = a2 ; hh->reverse = FALSE ; }
		else
		  { hh->a1 = a2 ; hh->a2 = a1 ; hh->reverse = TRUE ; }
	      }
          bsDestroy (Map) ;
        }      
    }

  arrayDestroy (maps) ;
  arrayDestroy (units) ;

  arraySort (hits, cDNAOrderByA1) ; /* cdnaclone then est then a1 then x1 */
  return (hits) ;
} /* cDnaGetAllIntMap */

/*********************************************************************/

static void cDnaGetOneSources (KEYSET sources, KEY gene, int dx)
{
  OBJ Source = 0 ; 
  KEY source = 0, old = 0 ;
  int a1, a2, a0 = 0, b1, b2, x1, x2, i, j = keySetMax(sources) ;
  Array aa = 0 ;
  BSunit *uu ;
  KEY _Is_predicted_gene = str2tag("Is_predicted_gene") ;
  source = keyGetKey (gene, _Genomic_sequence) ;
  if (source && (Source = bsCreate (source)))
    {
      if (bsFindKey (Source, _Transcribed_gene, gene) &&
          bsGetData (Source, _bsRight, _Int, &a1) &&
          bsGetData (Source, _bsRight, _Int, &a2))
        a0 = (a1 + a2)/2 ;
      else
        source = 0 ;
      bsDestroy (Source) ;
    }

  while (source) /* recurse upwards */
    {
      keySet (sources, j++) = source ; /* register */
      if ((Source = bsCreate (source)))
        {
          if (old && 
              bsFindKey (Source, _Subsequence, old) &&
              bsGetData (Source, _bsRight, _Int, &x1) &&
              bsGetData (Source, _bsRight, _Int, &x2))
            { /* gene gene position in source */
              if (x1 < x2) a0 = x1 + a0 - 1 ;
              else a0 = x1 - a0 + 1 ;
            }
          aa = arrayReCreate (aa, 30, BSunit) ;
          bsGetArray (Source, _Subsequence, aa, 3) ;
          for (i = 0 ; i < arrayMax(aa) ; i += 3) 
            {
              uu = arrp (aa, i, BSunit) ;
              b1 = uu[1].i ; b2 = uu[2].i ;
              if (b1 > b2) { int bb = b1 ; b1 = b2 ; b2 = bb ; }
              if (b1 < a0 + dx && b2 > a0 - dx && ! keyFindTag (uu[0].k, _Is_predicted_gene))
                keySet (sources, j++) = uu[0].k ; /* register */
            }
          old = source ;
          source = 0 ;
          bsGetKey (Source, _Source, &source) ;
          bsDestroy (Source) ;
        }
      else
        source = 0 ;
    }
  arrayDestroy (aa) ;
  return ;
} /* getOneSources  */

/*********************************************************************/

static Array  cDnaGetSources (KEYSET genes, int dx)
{
  KEYSET sources = keySetCreate () ;
  int nn = keySetMax (genes) ;
  
  while (nn--)
    cDnaGetOneSources (sources, keySet (genes, nn), dx) ;

  keySetSort (sources) ;
  keySetCompress (sources) ;
  return sources ;
}

/*********************************************************************/
/* select short predicted genes that cannot really be found using ESTs are are often embedded in longer mRNAs */
static BOOL isPredictedGeneMirLike (KEY seq)
{
  static BQL *bql = 0 ;
  static KEYSET ksIn = 0, ksOut = 0 ; /* no need to realloc each time */
  BOOL ok ;

  if (! bql)
    {      /* precompile the complex query */
      const char *ccp = "select s from s in @ where s#miRNA || s#vault_RNA || s#G26RNA || s#u21RNA || s#snoRNA || s#snRNA || s#scRNA || s#guide_RNA || s#vault_RNA || s#vault_RNA || s#Y_RNA || s#stRNA || s#RNase_P_RNA || s#RNase_MRP_RNA || s#telomerase_RNA || ((COUNT s->source_exons == 1) && s->source_exons[1] == 1 &&  s->source_exons[2] < 300) " ;
      bql = bqlCreate (FALSE, 0) ;
      if (! bqlParse (bql, ccp, FALSE)) 
	{ 
	  freeOutf ("// bql parse error in query %s\n%s\n", ccp,  bqlError (bql)) ;
	  messcrash ("// bad query in function isMirLike, please debug and recompile") ;
	}
      ksIn = keySetCreate () ;
      ksOut = keySetCreate () ;
    }

  keySet (ksIn, 0) = seq ; keySetMax (ksIn) = 1 ; 
  
  bqlRun (bql, ksIn, ksOut) ;
  ok = ksOut && keySetMax (ksOut) ? TRUE : FALSE ;

  return ok ;
} /* isPredictedGeneMirLike */

/*********************************************************************/

static KEY mrnaAnalyseExactGenefinder (KEY g1, int a1, int a2, BOOL isUp1, KEY g2, int b1, int b2, BOOL isUp2, BOOL *isExactp
				       , KEYSET exactMrnas, KEYSET matchingMrnas, KEYSET touchingMrnas, KEYSET mirLike)
{  
  int i, j, itr, score, bestScore = 0, exonScore, intronScore, endScore ;
  KEYSET mrnas = queryKey (g1, "Follow Mrna") ;
  KEY mrna = 0, exactMrna = 0, approximateMrna = 0 ;
  KEY _PolyA_found = str2tag ("PolyA_found") ;
  BOOL ok = FALSE ;
  BOOL isMirLike = isPredictedGeneMirLike (g2) ;
  OBJ G1, G2 ;
  BSunit *uu, *vv ;
  static Array aa = 0, pA = 0, bb = 0, map = 0 ;

  for (itr = 0 ; itr < keySetMax (mrnas) ; itr++)
    {
      int lenA = 0, lenB = 0, nIntronsA = -1, nIntronsB = -1 ;
      exonScore = intronScore = endScore = score = 0 ;
      mrna = keySet (mrnas, itr) ;
      G1 = bsCreate (mrna) ; 
      G2 = bsCreate (g2) ; 
      aa = arrayReCreate (aa, 40, BSunit) ;
      pA = arrayReCreate (pA, 40, BSunit) ;
      bb = arrayReCreate (bb, 40, BSunit) ;
      map = arrayReCreate (map, 4, BSunit) ;
      bsGetArray (G1, _Splicing, aa, 5) ;
      bsGetArray (G1, _IntMap, map, 3) ;
      bsGetArray (G2, _Source_Exons, bb, 4) ;
      
      if (arrayMax (map) < 3)
	continue ;
      uu = arrp (map, 0, BSunit) ;
      a1 = uu[1].i ; a2 = uu[2].i ; /* coordinates of the mRNA on the genome */
      
      for (i = 0 ; i < arrayMax(aa) ; i += 5)
	{
          uu = arrp (aa, i, BSunit) ;
	  if (i == arrayMax (aa) - 5)
	    {
	      /* last mRNA exons, reclip on the polyA */
	      int dx = 0, maxA = 0 ;
	      bsGetArray (G1, _PolyA_found, pA, 3) ;
	      for (int j = 0 ; j < arrayMax (pA) ; j+= 3)
		{
		  BSunit *pp = arrp (pA, j, BSunit) ;
		  if (pp[2].i > maxA)
		    { dx = pp[0].i ; maxA = pp[2].i ; }
		}
	      uu[1].i += dx ;
	    }
          if (a1 < a2)
            { uu[0].i = a1 + uu[0].i - 1 ; uu[1].i = a1 + uu[1].i - 1 ; }
	  else
            { int tmp = a1 - uu[0].i + 1 ; uu[0].i = a1 - uu[1].i + 1 ; uu[1].i = tmp ; }
	  if (strstr (name(uu[4].k),"xon")) uu[3].i = 1 ;
          else uu[3].i = 0 ;
        }
      /* keep happy few */
      for (i = j = 0 ; j < arrayMax(aa) ; j += 5)
        {
	  int k ;
	  
          vv = arrp (aa, j, BSunit) ;
	  if (! vv[3].i)
	    continue ;
	  uu = arrp (aa, i, BSunit) ;
	  if (i < j)
	    for (k = 0 ; k < 5 ; k++)
	      uu[k].i = vv[k].i ;
	  i += 5 ;
	}
      arrayMax(aa) = i ;
      

      /* merge contiguous CDS UTR in predicted genes */
      for (j = 0 ; j < arrayMax(bb) ; j += 4)
        {
          vv = arrp (bb, j, BSunit) ;
	  vv[3].i = 1 ;
	}

      for (j = 0 ; j < arrayMax(bb) ; j += 4)
        {
	  int j1, c2 ;
          vv = arrp (bb, j, BSunit) ;
	  if (! vv[3].i)
	    continue ;
	  for (c2 = vv[1].i, j1 = j + 4 ; j1 < arrayMax(bb) ; j1 += 4)
	    {
	      BSunit *xx ;
	      xx =  arrp (bb, j1, BSunit) ;
	      if ( xx[0].i == c2 + 1 && xx[3].i)
		{ xx[3].i = 0 ; vv[1].i = c2 = xx[1].i ; }
	      else
		break ;
	    }
	}
      /* keep happy few */
      for (i = j = 0 ; j < arrayMax(bb) ; j += 4)
        {
	  int k ;
	  
          vv = arrp (bb, j, BSunit) ;
	  if (! vv[3].i)
	    continue ;
	  uu = arrp (bb, i, BSunit) ;
	  if (i < j)
	    for (k = 0 ; k < 4 ; k++)
	      uu[k].i = vv[k].i ;
	  i+= 4 ;
	}
      arrayMax(bb) = i ;
      

      /* move to genome coordinates */
      for (j = 0 ; j < arrayMax(bb) ; j += 4)
	{
          vv = arrp (bb, j, BSunit) ;
          if (! isUp2)
            { vv[0].i = b1 + vv[0].i - 1 ; vv[1].i = b1 + vv[1].i - 1 ; }
          else
            { int tmp = b2 - vv[0].i + 1 ; vv[0].i = b2 - vv[1].i + 1 ; vv[1].i = tmp ; } 
	}
      
      ok = FALSE ;
      if (aa && bb)
	{
	  int iMax = arrayMax (aa) ;
	  int jMax = arrayMax (bb) ;
	  
	  /* length of the mRNA */
	  for (i = 0 ; i < iMax ; i += 5)
	    {
	      uu = arrp (aa, i, BSunit) ; 
	      if (uu[3].i)
		{
		  lenA += uu[1].i - uu[0].i + 1 ;
		  nIntronsA++ ;
		}
	    }
	  /* length of the Genefinder */
	  for (j = 0 ; j < jMax ; j += 4)
	    {
	      vv = arrp (bb, j, BSunit) ; 
	      if (vv[3].i)
		{
		  lenB += vv[1].i - vv[0].i + 1 ;
		  nIntronsB++ ;
		}
	    }
	  ok = TRUE ;
	  
	  for (i = j = 0 ; i < iMax && j < jMax ;)
	    {
	      uu = arrp (aa, i, BSunit) ; 
	      vv = arrp (bb, j, BSunit) ; 
	      if (!uu[3].i) { i += 5 ; continue ; }
	      if (!vv[3].i) { j += 4 ; continue ; }
	      /* take next exon */ 
	      /* look for intersection */
	      {
		int c1 = uu[0].i > vv[0].i ?   uu[0].i : vv[0].i ;
		int c2 = uu[1].i < vv[1].i ?   uu[1].i : vv[1].i ;
		if (c1 < c2)
		  {
		    score += c2 - c1 + 1 ;
		    exonScore += c2 - c1 + 1 ;
		    if (i + j == 0 && 2 * score > uu[3].i - uu[2].i ) /* same first exon */
		      { score += 100 ; endScore++ ; }
		    if (i + 5 == iMax && j + 4 == jMax && (10 * exonScore > 3* lenA || 10 * exonScore > 3 * lenB)) /* same last exon */
		      { score += 100 ; endScore++ ; }
		  }
		if (! isUp2)
		  {
		    if (i >= 5 && j >= 4 && uu[-4].i == vv[-3].i && uu[0].i == vv[0].i) /* same intron */
		      { score += 100 ; intronScore++ ; }
		  }
		else
		  {
		    if (i >= 5 && j >= 4 && uu[-3].i == vv[-2].i && uu[1].i == vv[1].i) /* same intron */
		      { score += 100 ; intronScore++ ; }
		  }
		
	      }

	      if (! isUp2)
		{
		  if (uu[1].i == vv[1].i) { i+= 5 ; j += 4 ; }
		  else if (uu[1].i < vv[1].i) i+= 5 ;
		  else j += 4 ;
		}
	      else
		{
		  if (uu[0].i == vv[0].i) { i+= 5 ; j += 4 ; }
		  else if (uu[0].i > vv[0].i) i+= 5 ;
		  else j += 4 ;
		}
	    }
	    
	  if (score > bestScore)
	    { 
	      bestScore = score ;
	      approximateMrna = mrna ; 
	      keySetMax (matchingMrnas) = 0 ;   /* reset */
	      keySetMax (mirLike) = 0 ;         /* reset */
	      keySetMax (exactMrnas) = 0 ;      /* reset */
	    }
	  
	  if (score == bestScore)
	    {
	      if (isMirLike)
		keySet (mirLike, 0) = mrna ;
	      else
		{
		  ok = TRUE ;
		  /* do we have the very same introns */
		  if (endScore < 2 || nIntronsA != nIntronsB || intronScore != nIntronsA)
		    ok = FALSE ;
		  /* do we have a quasi complete overlap */
		  if (10 * exonScore < 8 * lenA || 10 * exonScore < 8 * lenB)
		    ok = FALSE ;
		  if (ok)
		    {
		      exactMrna = mrna ; 
		      if (! isMirLike)
			keySet (exactMrnas, keySetMax (exactMrnas)) = mrna ;
		    }
		  
		  ok = FALSE ;
		  /* do we have the same introns */
		  if (intronScore && (2 * intronScore >= nIntronsA || 2 * intronScore >= nIntronsA))
		    ok = TRUE ;
		  /* do we have a substantial overlap */
		  if (10 * exonScore > 6 * lenA || 10 * exonScore > 6 * lenB)
		    ok = TRUE ;
		  /* do we have the same first and last exons */
		  if (endScore >= 2)
		    ok = TRUE ;
		  
		  if (ok)
		    keySet (matchingMrnas, keySetMax (matchingMrnas)) = mrna ;
		}
	    }
	  if (!isMirLike && score > 100) 
	    keySet (touchingMrnas, keySetMax (touchingMrnas)) = mrna ;
	}
      
      bsDestroy (G1) ;
      bsDestroy (G2) ;
    }
  keySetDestroy (mrnas) ; 

  if (exactMrna) 
    { *isExactp = TRUE ; return exactMrna ; }
  else 
    *isExactp = 0 ;
  return approximateMrna ; /* they intersect but we did not check the frame */
} /* mrnaAnalyseExactGenefinder */

/*********************************************************************/

static KEY mrnaAnalyseExactGenefinder2Product (KEY mrna, int a1, int a2, BOOL isUp1, KEY g2, int b1, int b2, BOOL isUp2, BOOL *isExactp)
{  
  int i, j, k, p1, p2, score, bestScore = 0 ;
  KEY product = 0, approximateProduct = 0 ;
  BOOL ok = FALSE ;
  OBJ G1, G2 ;
  BSunit *uu, *vv, *ww ;
  static Array aa = 0, bb = 0, pp = 0,  map = 0  ;

  while (1)
    {
      G1 = bsCreate (mrna) ; 
      G2 = bsCreate (g2) ;  
      map = arrayReCreate (map, 4, BSunit) ;
      bb = arrayReCreate (bb, 40, BSunit) ;
      pp = arrayReCreate (pp, 40, BSunit) ;
      bsGetArray (G2, _Source_Exons, bb, 4) ;
      bsGetArray (G1, _Product, pp, 5) ;
      bsGetArray (G1, _IntMap, map, 3) ;
   
      if (arrayMax (map) < 3)
	break ;
      uu = arrp (map, 0, BSunit) ;
      a1 = uu[1].i ; a2 = uu[2].i ; /* coordinates of the mRNA on the genome */
  
     /* merge contiguous CDS UTR in predicted genes */
      for (j = 0 ; j < arrayMax(bb) ; j += 4)
        {
	  int j1, c2 ;
          vv = arrp (bb, j, BSunit) ;

	  /* get the coding part of the predicted geme */
	  if (strstr (name(vv[2].k),"UTR")) vv[3].i = 0 ;
	  else vv[3].i = 1 ; /* we accept no qualifier as meaning CDS */

	  if (! vv[3].i)
	    continue ;
	  for (c2 = vv[1].i, j1 = j + 4 ; j1 < arrayMax(bb) ; j1 += 4)
	    {
	      BSunit *xx =  arrp (bb, j1, BSunit) ;
	      if ( xx[0].i == c2 + 1 && xx[3].i)
		{ xx[3].i = 0 ; vv[1].i = c2 = xx[1].i ; }
	      else
		break ;
	    }
	}
      /* move to genome coordinates */
      for (j = 0 ; j < arrayMax(bb) ; j += 4)
	{
	  vv = arrp (bb, j, BSunit) ;
	  if (! isUp2)
	    { vv[0].i = b1 + vv[0].i - 1 ; vv[1].i = b1 + vv[1].i - 1 ; }
	  else
	    { int tmp = b2 - vv[0].i + 1 ; vv[0].i = b2 - vv[1].i + 1 ; vv[1].i = tmp ; }
	}       
      
      /* loop on all products */
      for (k = 0 ; ! ok && k < arrayMax (pp) ; k+= 5)
	{
	  score = 0 ;
	  ww = arrp (pp, k, BSunit) ;
	  product = ww[0].k ;
	  p1 = ww[3].i ; p2 = ww[4].i ; /* coordomnnes du produit sur le premrna */
	  
	  /* get the exons of the mRNA and clip them according to the product */
	  aa = arrayReCreate (aa, 40, BSunit) ;
	  bsGetArray (G1, _Splicing, aa, 5) ;
	  
	  for (i = 0 ; i < arrayMax(aa) ; i += 5)
	    {
	      int c1, c2 ;
	      uu = arrp (aa, i, BSunit) ;
	      if ( ! strstr (name(uu[4].k),"xon"))
		{ uu[3].i = 0 ; continue ; } 

	      /* are we in the CDS of this product, look in the spliced mRNA coordinates */
	      c1 = uu[0].i > p1 ? uu[0].i : p1 ;
	      c2 = uu[1].i < p2 ? uu[1].i : p2 ;
	      if (c1 > c2)
		{ 
		  uu[3].i = 0  ;
		  continue ; /* this exon is fully UTR */
		}
	      /* now move to coordicnates on the pre-MRNA */

	      if (a1 < a2)
		{ uu[0].i = a1 + c1 - 1 ; uu[1].i = a1 + c2 - 1 ; }
	      else
		{ uu[0].i = a1 - c2 + 1 ; uu[1].i = a1 - c1 + 1 ; }
	      uu[3].i = 1 ;
	    }

	  ok = TRUE ;
	  if (1)
	    {
	      ok = TRUE ; 
	      /* do the non UTR exons match each other */
	      for (i = 0, j = 0 ; i < arrayMax(aa) && j < arrayMax(bb) ; )
		{
		  uu = arrp (aa, i, BSunit) ; 
		  vv = arrp (bb, j, BSunit) ; 
		  if (!uu[3].i) { i += 5 ; continue ; }
		  if (!vv[3].i) { j += 4 ; continue ; }
		  /* take next exon */
		  /* look for intersection */
		  {
		    int c1 = uu[0].i > vv[0].i ?   uu[0].i : vv[0].i ;
		    int c2 = uu[1].i < vv[1].i ?   uu[1].i : vv[1].i ;
		    if (c1 < c2)
		      score += c2 - c1 + 1 ;
		    if (score > bestScore)
		      { 
			bestScore = score ;
			approximateProduct = product ;
		      }
		 }
		  /* look for exact match */
		  if (uu[0].i != vv[0].i ||
		      uu[1].i != vv[1].i)
		    ok = FALSE ; 
		  if (uu[1].i <= vv[1].i) i+= 5 ;
		  if (uu[1].i >= vv[1].i) j += 4 ;
		  /* i += 5 ; j += 4 ;  */
		}
	      /* do we have remaining coding exons */
	      for ( ; ok && i < arrayMax(aa) ; )
		{
		  uu = arrp (aa, i, BSunit) ; 
		  if (!uu[3].i) { i += 5 ; continue ; } /* jumping UTR is ok */
		  ok = FALSE ;
		}
	      for ( ; ok && j < arrayMax(bb) ; )
		{
		  vv = arrp (bb, j, BSunit) ; 
		  if (!vv[3].i) { j += 4 ; continue ; }
		  ok = FALSE ;
		}
	    }
	}
	  
      if (! ok) 
	product = 0 ;
      break ;
    }

  bsDestroy (G1) ;
  bsDestroy (G2) ;
  if (ok) 
    { *isExactp = TRUE ; return product ;}
  else 
    *isExactp = FALSE ;
  return approximateProduct ; /* they intersect but we did not check the frame */
} /* mrnaAnalyseExactGenefinder2Product */

/*********************************************************************/

static void mrnaAnalyseGenefinder (Stack s, KEY gene, Array allPg, Array allTg, int type) 
{
  HIT *hh ;
  int i, a1, a2, b1, b2, dx, da, db ;
  KEY map = 0, newname = 0 ;
  BOOL isUp = FALSE ;

  a1 = a2 = b1 = b2 = 0 ;

  if (arrayMax(allTg))
    for (i = 0, hh = arrp (allTg, 0, HIT) ; i < arrayMax(allTg) ; i++, hh++)
      {
        if (hh->gene == gene)
          { a1 = hh->a1 ; a2 = hh->a2 ; map = hh->est ; isUp = hh->reverse ; break ; }
      }
  if (!map) return ;
  
  if (arrayMax(allPg))
    for (i = 0, hh = arrp (allPg, 0, HIT) ; i < arrayMax(allPg) ; i++, hh++)
      {
        if (! hh->gene) continue ;
        if (hh->gene == gene) continue ;
        if (map != hh->est || isUp != hh->reverse) continue ; 
        b1 = hh->a1 ; b2 = hh->a2 ;
        
        dx = ( b2 < a2 ? b2 : a2) -  (b1 > a1 ? b1 : a1) ;        
        da = a2 - a1 ; db = b2 - b1 ;
        
        if (
            (type == 1 && (4*dx > da || 4*dx > db) && 
             ( mrnaAnalysePg2Tg (hh->gene, b1, b2, hh->reverse, gene, a1, a2, isUp) > 80))
            ||
            (type == 2 && (3*dx > da || 3*dx > db)) ||
            (type == 3 && (4*dx > 3*da && 4*dx > 3*db)) /* see  mrnaNameGeneByPhenotypicGene */
            )
          {  /* register overlaps */
            switch (type)
              {
              case 1:
                catText (s, messprintf("Transcribed_gene \"%s\"\n", name(gene))) ;
                catText (s, messprintf("Matching_genefinder_gene \"%s\" \n", name (hh->gene))) ;
		catText (s, "\n") ;
                if (dx > 0 && TRUE)
                  {
                    KEY mrna = 0, product = 0 ;
		    int imm ;
		    BOOL isExact = FALSE ;
		    KEYSET exactMrnas = keySetCreate () ;
		    KEYSET matchingMrnas = keySetCreate () ;
		    KEYSET touchingMrnas = keySetCreate () ;
		    KEYSET mirLike = keySetCreate () ;
 		    KEYSET matchingMrnasBis = 0 ;
		    KEYSET touchingMrnasBis = 0 ;
		    KEYSET touchingMrnasTer = 0 ;

		    mrna = mrnaAnalyseExactGenefinder (gene, a1, a2, isUp, hh->gene, b1, b2, hh->reverse, &isExact, exactMrnas, matchingMrnas, touchingMrnas, mirLike) ;          
                    if (mrna) 
                      {
			matchingMrnasBis = keySetMINUS (matchingMrnas, exactMrnas) ;
			touchingMrnasBis = keySetMINUS (touchingMrnas, matchingMrnas) ;
			touchingMrnasTer = keySetMINUS (touchingMrnasBis, exactMrnas) ;

			for (imm = 0 ; imm < keySetMax (mirLike) ; imm++)
			  {
			    KEY mrna1 = keySet (mirLike, imm) ;
			    catText (s, messprintf("mRNA \"%s\"\n", name(mrna1))) ;
			    catText (s, messprintf("%s \"%s\"\n\n"
						   , "Matching_short_RNA" 
						   , name (hh->gene))) ;
			  }
			for (imm = 0 ; imm < keySetMax (exactMrnas) ; imm++)
			  {
			    KEY mrna1 = keySet (exactMrnas, imm) ;
			    catText (s, messprintf("mRNA \"%s\"\n", name(mrna1))) ;
			    catText (s, messprintf("%s \"%s\"\n\n"
						   , "Identical_to_genefinder" 
						   , name (hh->gene))) ;
			  }
			for (imm = 0 ; imm < keySetMax (matchingMrnasBis) ; imm++)
			  {
			    KEY mrna1 = keySet (matchingMrnasBis, imm) ;
			    catText (s, messprintf("mRNA \"%s\"\n", name(mrna1))) ;
			    catText (s, messprintf("%s \"%s\"\n\n"
						   , "Matching_genefinder"
						   , name (hh->gene))) ;
			  }
			for (imm = 0 ; imm < keySetMax (touchingMrnasTer) ; imm++)
			  {
			    KEY mrna1 = keySet (touchingMrnasBis, imm) ;
			    catText (s, messprintf("mRNA \"%s\"\n", name(mrna1))) ;
			    catText (s, messprintf("%s \"%s\"\n\n"
						   , "Touching_genefinder"
						   , name (hh->gene))) ;
			  }
			product = mrnaAnalyseExactGenefinder2Product (mrna, a1, a2, isUp, hh->gene, b1, b2, hh->reverse, &isExact) ;            
			if (product)
			  {
			    catText (s, messprintf("Product \"%s\"\n", name(product))) ;
			    if (isExact)
			      catText (s, messprintf("Identical_to_gene \"%s\"\n", name (hh->gene))) ;
			    catText (s, messprintf("Matching_gene \"%s\" \n\n", name (hh->gene))) ;
			  }
		      }

		    keySetDestroy (exactMrnas) ;
		    keySetDestroy (matchingMrnas) ;
		    keySetDestroy (touchingMrnas) ;
		    keySetDestroy (mirLike) ;
		    keySetDestroy (matchingMrnasBis) ;
		    keySetDestroy (touchingMrnasBis) ;
		    keySetDestroy (touchingMrnasTer) ;
                  }
                break ;
              case 2:
                break ;
              case 3:
                catText (s, messprintf("Transcribed_gene \"%s\"\n", name(gene))) ;
                if ((newname = keyGetKey (hh->gene, _NewName)))
                  catText (s, messprintf ("NewName %s\n", name(newname))) ;
                catText (s, messprintf("Gene \"%s\" \n\n", name (hh->gene))) ;
                break ;
              }
          }
        
      }
}

/*********************************************************************/

static int mrnaAnalyseAntisens (KEY g1, int a1, int a2, BOOL isUp1, KEY g2, int b1, int b2, BOOL isUp2)
{  
  int dx = 0, i, j, c1, c2 ;
  OBJ G1, G2 ;
  BSunit *uu, *vv ;
  static Array aa = 0, bb = 0 ;

  G1 = bsCreate (g1) ; 
  G2 = bsCreate (g2) ; 
  aa = arrayReCreate (aa, 40, BSunit) ;
  bb = arrayReCreate (bb, 40, BSunit) ;
  bsGetArray (G1, _Splicing, aa, 4) ;
  bsGetArray (G2, _Splicing, bb, 4) ;

  for (i = 0 ; i < arrayMax(aa) ; i += 4)
    {
      uu = arrp (aa, i, BSunit) ;
      if (! isUp1)
        { uu[0].i = a1 + uu[0].i - 1 ; uu[1].i = a1 + uu[1].i - 1 ; }
      else
        { int tmp = a2 - uu[0].i + 1 ; uu[0].i = a2 - uu[1].i + 1 ; uu[1].i = tmp ; }
      if (strstr (name(uu[2].k),"xon")) uu[3].i = 1 ;
      else uu[3].i = 0 ;
    }

  for (j = 0 ; j < arrayMax(bb) ; j += 4)
    {
      vv = arrp (bb, j, BSunit) ;
      if (! isUp2)
        { vv[0].i = b1 + vv[0].i - 1 ; vv[1].i = b1 + vv[1].i - 1 ; }
      else
        { int tmp = b2 - vv[0].i + 1 ; vv[0].i = b2 - vv[1].i + 1 ; vv[1].i = tmp ; }
      if (strstr (name(vv[2].k),"xon")) vv[3].i = 1 ;
      else vv[3].i = 0 ;
    }

  c1 = a1 < b1 ? a1 : b1 ;
  c2 = a2 > b2 ?  a2 : b2 ;
  if (c2 > c1)
   {
     BitSet bb1 = bitSetCreate (c2 - c1 + 2000, 0) ;
     BitSet bb2 = bitSetCreate (c2 - c1 + 2000, 0) ;
     int u1, u2, v1, v2 ;

     /* paint gene 1 */
     for (i = 0 ; i < arrayMax(aa) ; i += 4)
       {
         uu = arrp (aa, i, BSunit) ;  
         if (!uu[3].i) continue ;
         u1 = uu[0].i ; u2 = uu[1].i ;
         for (j = u1 ; j <= u2 && j - c1 + 1001 < c2 - c1 + 2000  ; j++)
           if (j - c1 + 1001 > 0)
             bitSet (bb1, j - c1 + 1001) ;
       }
     /* conditionally paint gene 2 */
     for (i = 0 ; i < arrayMax(bb) ; i += 4)
       {
         vv = arrp (bb, i, BSunit) ; 
         if (!vv[3].i) continue ;
         v1 = vv[0].i ; v2 = vv[1].i ;
         for (j = v1 ; j <= v2 && j - c1 + 1001 < c2 - c1 + 2000 ; j++)
           if (j - c1 + 1001 > 0 && bit (bb1,  j - c1 + 1001))
             bitSet (bb2, j - c1 + 1001) ;
       }
     dx = bitSetCount (bb2) ;
     bitSetDestroy (bb1) ;
     bitSetDestroy (bb2) ;
   }

  bsDestroy (G1) ;
  bsDestroy (G2) ;
  return dx ;
}

/*********************************************************************/

static int mrnaAnalysePg2Tg (KEY g1, int a1, int a2, BOOL isUp1, KEY g2, int b1, int b2, BOOL isUp2)
{  
  int dx = 0, i, j, k ;
  OBJ G1, G2 ;
  BSunit *uu, *vv ;
  HIT *up, *vp ;
  static Array aa = 0, bb = 0 ;
  static Array h1 = 0, h2 = 0 ;

  G1 = bsCreate (g1) ; 
  G2 = bsCreate (g2) ; 
  aa = arrayReCreate (aa, 40, BSunit) ;
  bb = arrayReCreate (bb, 40, BSunit) ;
  h1 = arrayReCreate (h1, 40, HIT) ;
  h2 = arrayReCreate (h2, 40, HIT) ;
  bsGetArray (G1, _Source_Exons, aa, 4) ;
  bsGetArray (G2, _Splicing, bb, 4) ;

  for (i = k = 0 ; i < arrayMax(aa) ; i += 4)
    {
      uu = arrp (aa, i, BSunit) ;
      if (! isUp1)
        { uu[0].i = a1 + uu[0].i - 1 ; uu[1].i = a1 + uu[1].i - 1 ; }
      else
        { int tmp = a2 - uu[0].i + 1 ; uu[0].i = a2 - uu[1].i + 1 ; uu[1].i = tmp ; }
      uu[2].k = _Exon ; uu[3].i = 1 ;
      up = arrayp (h1, k++, HIT) ;
      up->a1 = a1 ; up->a2 = a2 ;
    }

  for (j = k = 0 ; j < arrayMax(bb) ; j += 4)
    {
      vv = arrp (bb, j, BSunit) ;
      if (strstr (name(vv[2].k),"xon"))
	{
	  if (! isUp2)
	    { vv[0].i = b1 + vv[0].i - 1 ; vv[1].i = b1 + vv[1].i - 1 ; }
	  else
	    { int tmp = b2 - vv[0].i + 1 ; vv[0].i = b2 - vv[1].i + 1 ; vv[1].i = tmp ; }
	  vp = arrayp (h2, k++, HIT) ;
	  vp->a1 = a1 ; vp->a2 = a2 ;
	}
    }
  
  arraySort (h1, cDNAOrderByA1) ;
  arraySort (h2, cDNAOrderByA1) ;
  for (i = 0 ; i < arrayMax(h1) ; i++)
    {
      up = arrp (h1, i, HIT) ;
      for (j = 0 ; j < arrayMax(h2) ; j ++)
        {
          int u1, u2 ;
	  vp = arrp (h2, j, HIT) ;
          u1 = up->a1 > vp->a1 ? up->a1 : vp->a1 ;
	  u2 = up->a2 < vp->a2 ? up->a2 : vp->a2 ;
          if (u2 >= u1)
            dx += u2 - u1 + 1 ;
	  if ( up->a1 == vp->a1 || up->a2 == vp->a2)
	    dx += 100 ;
	  if (dx > 80) break ;
	  if (up->a2 < vp->a1) break ;
        }
      if (dx > 80) break ;
    }

  bsDestroy (G1) ;
  bsDestroy (G2) ;
  return dx ;
}/*  mrnaAnalysePg2Tg */

/*********************************************************************/

static void mrnaAnalyseOverlap (Stack s, KEY gene, Array allTg) 
{
  HIT *hh ;
  int i, a1, a2, b1, b2, dx, da, db ;
  KEY map = 0, gene1 = 0, gene2 ;
  BOOL isUp = FALSE ;
  int   OPERON_MAX = 350 ;

  if (!arrayMax(allTg))
    return ;

  a1 = a2 = b1 = b2 = 0 ;
  for (i = 0, hh = arrp (allTg, 0, HIT) ; i < arrayMax(allTg) ; i++, hh++)
    {
      if (hh->gene == gene)
        { a1 = hh->a1 ; a2 = hh->a2 ; map = hh->est ; isUp = hh->reverse ; break ; }
    }
  if (!map) return ;
  
  /* find gene downstream chromosome wise */
  for (i = 0, gene1 = 0, hh = arrp (allTg, 0, HIT) ; i < arrayMax(allTg) ; i++, hh++)
    {
      if (hh->gene == gene) continue ;
      if (map != hh->est) continue ; 
      if (hh->a1 < a2) continue ;
      if (gene1 && hh->a1 > b1) continue ;
      if (hh->reverse != isUp) continue ;

      b1 = hh->a1 ; b2 = hh->a2 ; gene1 = hh->gene ;
    }

  if (gene1)
    {          
      if (!isUp)
        { gene2 = gene1 ; gene1 = gene ; }
      else
        { gene2 = gene ; }
      dx = b1 - a2 ;

      catText (s, messprintf("Transcribed_gene \"%s\"\n", name(gene1))) ;
      catText (s, messprintf("Next_gene_in_cis \"%s\" %d\n\n", name (gene2), dx)) ;
      catText (s, messprintf("Transcribed_gene \"%s\"\n", name(gene2))) ;
      catText (s, messprintf("Previous_gene_in_cis \"%s\" %d\n\n", name (gene1), -dx)) ;

      if (dx < OPERON_MAX)
        {  /* register Operon, twice to have the value twice */
          catText (s, messprintf("Transcribed_gene \"%s\"\n", name(gene1))) ;
          catText (s, messprintf("Possible_operon \"%s\" %d\n\n", name (gene2), dx)) ;
          catText (s, messprintf("Transcribed_gene \"%s\"\n", name(gene2))) ;
          catText (s, messprintf("Possible_operon \"%s\" %d\n\n", name (gene1), -dx)) ;
        }
    }

  /* find gene upstream chromosome wise */
  for (i = 0, gene1 = 0, hh = arrp (allTg, 0, HIT) ; i < arrayMax(allTg) ; i++, hh++)
    {
      if (hh->gene == gene) continue ;
      if (map != hh->est) continue ; 
      if (hh->a2 > a1) continue ;
      if (gene1 && hh->a2 < b2) continue ;
      if (hh->reverse != isUp) continue ;

      b1 = hh->a1 ; b2 = hh->a2 ; gene1 = hh->gene ;
    }

  if (gene1)
    {          
      if (!isUp)
        { gene2 = gene ; }
      else
        { gene2 = gene1 ; gene1 = gene ; }
      dx = a1 - b2 ;

      catText (s, messprintf("Transcribed_gene \"%s\"\n", name(gene1))) ;
      catText (s, messprintf("Next_gene_in_cis \"%s\" %d\n\n", name (gene2), dx)) ;
      catText (s, messprintf("Transcribed_gene \"%s\"\n", name(gene2))) ;
      catText (s, messprintf("Previous_gene_in_cis \"%s\" %d\n\n", name (gene1), -dx)) ;

      if (dx < OPERON_MAX)
        {  /* register Operon, twice to have the value twice */
          catText (s, messprintf("Transcribed_gene \"%s\"\n", name(gene1))) ;
          catText (s, messprintf("Possible_operon \"%s\" %d\n\n", name (gene2), dx)) ;
          catText (s, messprintf("Transcribed_gene \"%s\"\n", name(gene2))) ;
          catText (s, messprintf("Possible_operon \"%s\" %d\n\n", name (gene1), -dx)) ;
        }
    }

  /* now look for overlaps and antisens, not necessarilly unique */
  for (i = 0, hh = arrp (allTg, 0, HIT) ; i < arrayMax(allTg) ; i++, hh++)
    {
      if (hh->gene == gene) continue ;
      if (map != hh->est) continue ; 
      b1 = hh->a1 ; b2 = hh->a2 ;
          
      dx = ( b2 < a2 ? b2 : a2) -  (b1 > a1 ? b1 : a1) + 1 ;    /* extent of overlap */
      da = a2 - a1 ; db = b2 - b1 ;

      if (3*dx > da || 3*dx > db) /* && isUp != hh->reverse  for both strands, Dnaous genuis innovation, aug 22, 2002 */
        {  /* register overlaps */
          catText (s, messprintf("Transcribed_gene \"%s\"\n", name(gene))) ;
          catText (s, messprintf("Overlaps \"%s\" \n\n", name (hh->gene))) ;
        }
      if (dx > 0)
        { /* register antisens */
          int ddx = 0 ;

          ddx = mrnaAnalyseAntisens (gene, a1, a2, isUp, hh->gene, b1, b2, hh->reverse) ;          
          if (ddx > 20 && isUp != hh->reverse)
            {
              catText (s, messprintf("Transcribed_gene \"%s\"\n", name(gene))) ;
              catText (s, messprintf("Antisens_to \"%s\" %d\n\n", name (hh->gene), ddx)) ;
              catText (s, messprintf("Transcribed_gene \"%s\"\n", name(hh->gene))) ;
              catText (s, messprintf("Antisens_to \"%s\" %d\n\n", name (gene), ddx)) ;
            }
          if (ddx > 20 && isUp == hh->reverse) /* assymetric, so no error in cdnarealignGene */
            {
              if (da > db || (da == db  && gene < hh->gene))
                {
                  catText (s, messprintf("Transcribed_gene \"%s\"\n", name(gene))) ;
                  catText (s, messprintf("To_be_fused_with \"%s\" %d\n\n", name (hh->gene), ddx)) ;
                }
              else
                {
                  catText (s, messprintf("Transcribed_gene \"%s\"\n", name(hh->gene))) ;
                  catText (s, messprintf("To_be_fused_with \"%s\" %d\n\n", name (gene), ddx)) ;
                }
            }
        }
    }
}

/*********************************************************************/

static void mrnaAnalyseNewNames (Stack s, KEY gene, Array allTg) 
{
  HIT *hh ;
  int i ;
  KEY map = 0 ;
  
  if (!getPleaseNewNames ())
    return ;

  for (i = 0, hh = arrp (allTg, 0, HIT) ; i < arrayMax(allTg) ; i++, hh++)
    {
      if (hh->gene == gene)
        {  map = hh->est ; break ; }
    }
  if (!map) 
    return ;
}

/*********************************************************************/
/* look for operon, identity to genefinder, newNames etc */

void mrnaAnalyseNeighbours (KEYSET tgenes) 
{
  AC_HANDLE h = ac_new_handle () ;
  int dx = 1000000, nn  ;
  KEYSET sources = 0, tgs = 0,  pgs = 0, gns = 0, mrnas = 0, prods = 0 ;

  Array allTg = 0 ;
  Array allPg = 0 ; 
  Array allGn = 0 ;

  if (tgenes)
    {
      sources = cDnaGetSources (tgenes, dx) ; /* all cosmids dx around the genes */
      tgs = query (sources, "FOLLOW Transcribed_gene") ;
      pgs = query (sources, "{FOLLOW Subsequence ; CLASS Predicted_gene} SETOR {FOLLOW Genes ;Follow genefinder;}") ;
      gns = query (sources, "FOLLOW Genes ") ;

      mrnas = query (tgenes, "FOLLOW mRNA") ;
      prods = query (mrnas, "FOLLOW Product") ;
      
      allTg = cDnaGetIntMap (tgs, h) ;
      allPg = cDnaGetIntMap (pgs, h) ;
      allGn = cDnaGetIntMap (gns, h) ;
    }
  else
    {
      messcrash (" cDnaGetAllIntMap cannot work because IntMap is not XREF for all classes ") ;
      allTg = cDnaGetAllIntMap (_Transcribed_gene, h) ; /* NOT XREFed, cannot work */
    }

  Stack s = stackCreate (32000) ;

  nn = arrayMax (tgenes) ;
  while (nn--)
    {
      KEY tgene = keySet (tgenes, nn) ;
      catText (s, messprintf("Transcribed_gene %s\n", name(tgene))) ;
      catText (s, "-D Overlaps\n") ;
      catText (s, "-D Possible_operon\n") ;
      catText (s, "-D Antisens_to\n") ;
      catText (s, "-D Matching_genefinder_gene\n") ;  
      catText (s, "\n") ;
    }

  nn = arrayMax (prods) ;
  while (nn--)
    {
      KEY product = keySet (prods, nn) ;
      catText (s, messprintf("Product %s\n", name(product))) ;
      catText (s, "-D Identical_to_gene\n") ;
      catText (s, "-D Matching_gene\n") ;
      catText (s, "\n") ;
    }

  nn = arrayMax (mrnas) ;
  while (nn--)
    {
      KEY mrna = keySet (mrnas, nn) ;
      catText (s, messprintf("mRNA %s\n", name(mrna))) ;
      catText (s, "-D Touching_genefinder\n") ;
      catText (s, "-D Matching_genefinder\n") ;
      catText (s, "-D Identical_to_genefinder\n") ;
      catText (s, "-D Stop_of\n") ;
      catText (s, "\n") ;
    }

  nn = arrayMax (tgenes) ;
  while (nn--)
    {
      KEY tg = keySet (tgenes, nn) ;
      mrnaAnalyseGenefinder (s, tg, allPg, allTg, 1) ;
      mrnaAnalyseGenefinder (s, tg, allGn, allTg, 3) ;
      mrnaAnalyseOverlap (s, tg, allTg) ;
      mrnaAnalyseNewNames (s, tg, allTg) ;
    }

  parseBuffer (stackText(s,0), 0) ;
  stackDestroy (s) ;

  keySetDestroy (sources) ;
  keySetDestroy (tgs) ;
  keySetDestroy (prods) ;
  keySetDestroy (pgs) ;
  keySetDestroy (gns) ;
  keySetDestroy (mrnas) ;
  
  ac_free (h) ;

  return ;
} /* mrnaAnalyseNeighbours */

/*********************************************************************/

void mrnaAnalyseAllNeighbours (void)
{
  KEYSET tgenes = query (0, "Find Transcribed_gene") ;
  mrnaAnalyseNeighbours (tgenes) ;
  keySetDestroy (tgenes) ;
  return ;
} /* analyseAllNeighbours */

/*********************************************************************/
/* if metp != 0, export from met to stop and report *metp
 * else export the whole ORF 
 */
static Array mrnaPg2pep (Array dna, int *p1, int *p2, int *leup, int *metp)
{
  Array pep = 0 ;
  char *cp, *cq, *cr, cc ;
  int frame, max, i, i0, j, best, len;
  BOOL useLeu = mrnaPleaseUseLeucine () ;

  if (!dna || arrayMax (dna) < 3)
    return 0 ;
  
  best = 0 ; max = arrayMax (dna) ;
  for (frame = 0 ; frame < 3 ; frame++)
    {
      for (i = i0 = frame, cp = arrp (dna, i, char) ; i < max - 2 + 3 ; cp += 3, i += 3)
      {
        cc =  i < max - 2 ? codon(cp) : '*' ;
        if (cc == '*')
          {
            len = (i - i0)/3 ; 
            if (len > best)
              { best = len ; *p1 = i0 + 1 /* plato */; *p2 = i  ; }
            i0 = i + 3 ;
          }
      }
    }
 
  if (metp) *metp = 0 ;
  if (leup) *leup = 0 ;
  if (best && leup && metp && *p1 > 0)
    {
      int met, leu ;

      met = leu = 0 ;
      for (i = *p1 - 1, cp = arrp (dna, i, char), cq = cp+1, cr = cp+2 ; 
           i < *p2 - 32 ; /* impose at least 10 AA */
           cp += 3, cq += 3, cr += 3, i += 3)
        {
          if (
              *cp == A_ && *cq == T_ && *cr == G_
              )
            {
              met = i + 1 ; /* plato */
              break ;
            }
        }
      for (i = *p1 - 1, cp = arrp (dna, i, char), cq = cp+1, cr = cp+2 ; 
           useLeu && i < *p2 - MINI_LEU2MET ; /* impose at least 10 AA */
           cp += 3, cq += 3, cr += 3, i += 3)
        {
          if (*(cp+2) == G_ && 
              ((*(cp-3) == G_ && *(cp+3) == G_) || (*(cp-3) == A_)) &&
              (*(cp) == A_ || *(cp+1) == T_)
              )
            {
              leu = i + 1 ; /* plato */
              break ;
            }
        }
      if (leu && met && leu > met - MINI_LEU2MET)
        leu = met ;
      *leup = leu ;
      *metp = met ;
    }
  if (best > 0)
    {
      int pp1 = (leup && *leup) ? *leup : *p1 ;
      
      pep = arrayCreate (best + 2, char) ;
      for (i = pp1 - 1, cp = arrp (dna, i, char), j = 0 ; i < *p2 ; cp += 3, i += 3)
        array (pep, j++, char) = pepEncodeChar[(int)codon (cp)] ;
      if (pp1 != *p1)
        array (pep, 0, char) = pepEncodeChar[(int)'M'] ; /* force the Met */
    }
  return pep ;
}

/*********************************************************************/

void mrnaTranslateEst (KEYSET ks)
{
  int ii, p1, p2, leu, met, pepL, nX ;
  char *cp ;
  Array pep = 0, dna = 0 ;
  OBJ Est = 0, Product = 0 ;
  KEY product, est ;
  KEYSET ests = query (ks, "CLASS est") ;
  
  for (ii = 0 ; ii < keySetMax (ests) ; ii++)
    {
      est = keySet (ests, ii) ;
      if ((dna = dnaGet (est)) && (pep = mrnaPg2pep (dna, &p1, &p2, &leu, &met)))
        {
          pepL = arrayMax (pep) ; nX = 0 ;
          for (cp = arrp (pep, 0, char) ; *cp ; cp++)
            if (*cp == 'X') { nX++ ; pepL-- ; }        

          lexaddkey (messprintf("%s.est", name(est)), &product, _VmProduct) ;
          peptideStore (product, pep) ;
          if ((Product = bsUpdate (product)))
            {
              bsAddTag (Product, _Best_product) ;
              bsAddKey (Product, _From_EST, est) ;
              pepL *= 3 ; nX *= 3 ; /* store bp lengths */
              if (leu)
                {
                  bsAddData (Product, _Coding_length, _Int, &pepL) ;
                  pepL += (leu - p1) ;
                }
              bsAddData (Product, _Open_length, _Int, &pepL) ;
              if (nX)
                bsAddData(Product, _Coding_gap, _Int, &nX) ;
              if (leu > 0 && met >= leu)
                {
                  int aa_met = (met - leu)/3 + 1 ;
                  bsAddData (Product, _First_ATG, _Int, &aa_met) ;
                  if (aa_met == 1)
                    bsAddTag (Product, _At_position_1) ;
                }
              bsSave (Product) ;
            }
          mrnaDesignSavePepAndKantor (product, pep) ;
          if ((Est = bsUpdate (est)))
            {
              bsAddKey (Est, _EST_translation, product) ;
              bsAddData (Est, _bsRight, _Int, &p1) ;
              bsAddData (Est, _bsRight, _Int, &p2) ;
              bsAddData (Est, _bsRight, _Int, &leu) ;
              bsAddData (Est, _bsRight, _Int, &met) ;
              bsSave (Est) ;
            }
        }
      arrayDestroy (pep) ;
      arrayDestroy (dna) ;
    }
} /* rnaTranslateEst */

/*********************************************************************/

void mrnaTransferPg2PredictedMrna (KEYSET ks)
{
  KEYSET pgs = query (ks, "CLASS Predicted_gene") ;
  int ipgs, jj, a1, a2, x0, x1, x2, lasta2, pepL, dnaL, p1, p2, leu, met ; 
  KEY key, pg, mrna, product, source ;
  Array units = 0, dna = 0, pep = 0 ;
  char *cp, buf[300] ;
  OBJ Pg = 0, Source = 0, Mrna = 0, Product = 0 ;
  BSunit *uu ;
  BOOL isCoding ;

  units = arrayCreate (80, BSunit) ;
  for (ipgs = 0 ; ipgs < keySetMax (pgs) ; ipgs++)
    {
      pg = keySet (pgs, ipgs) ;

      { /* dna */ 
        dna = dnaGet (pg) ; 
        if (!dna) 
          continue ;
        dnaL = arrayMax (dna) ;
      }

      { /* pep */ 
        isCoding = FALSE ;
        pep = peptideGet (pg) ; p1 = 1 ; p2 = dnaL ; leu = met = 0 ;
	pep = mrnaPg2pep (dna, &p1, &p2, &leu, &met) ;
        if (pep) 
	  isCoding = TRUE ;
        else
          { arrayDestroy (dna) ; continue ; }
        if (leu) p1 = leu ;
        pepL = arrayMax (pep) ;
        for (cp = arrp (pep, 0, char) ; *cp ; cp++)
          if (*cp == 'X') pepL-- ;        
      } 
      
      Pg = bsCreate (pg) ;

      { /* Mrna part */
        sprintf (buf, "%s.pg", name(pg)) ;
        lexaddkey (buf, &mrna, _VmRNA) ;
        lexaddkey (buf, &product, _VmProduct) ;        
        dnaStoreDestroy (mrna, dna) ; dna = 0 ;
        if ((Mrna = bsUpdate (mrna)))
          {
            if (bsFindTag (Mrna, str2tag("Structure")))
              bsRemove (Mrna) ;

            bsAddKey (Mrna, str2tag("From_prediction"), pg) ;
            bsAddTag (Mrna, str2tag("Complete")) ;
            if (bsGetKey (Pg, str2tag ("Genetic_code"), &key))  /* k = pg->Genetic_code ;
                                                                   key = acTagClass (pg, "Genetic_code") ; */
              bsAddKey (Mrna, str2tag("Genetic_code"), key) ;

            bsAddKey (Mrna, _Product, product) ;

            bsAddData (Mrna, _bsRight, _Int, &p1) ;
            bsAddData (Mrna, _bsRight, _Int, &p2) ;

            a1 = p2 - p1 + 1 ;
            bsAddData (Mrna, str2tag ("Longest_CDS"), _Int, &a1) ;

            if (bsGetKey (Pg, _IntMap, &key) &&
                bsGetData (Pg, _bsRight, _Int, &a1) &&
                bsGetData (Pg, _bsRight, _Int, &a2))
              {
                bsAddKey (Mrna, _IntMap, key) ;
                bsAddData (Mrna, _bsRight, _Int, &a1) ;
                bsAddData (Mrna, _bsRight, _Int, &a2) ;
            
              }
            /* coding and splicing */
            x1 = x2 = 0 ; lasta2 = 1 ;
            if (bsGetArray (Pg, _Source_Exons, units, 2))
              for (jj = 0 ; jj < arrayMax (units) ; jj += 2)
                {
                  uu = arrayp (units, jj, BSunit) ;
                  a1 = uu[0].i ; a2 = uu[1].i ; /* unspliced coords */
                  if (a1 > lasta2)
                    {
                      x0 = lasta2 + 1 ;
                      bsAddData (Mrna, _Splicing, _Int, &x0) ;
                      x0 = a1 - 1 ;
                      bsAddData (Mrna, _bsRight, _Int, &x0) ;
                      x0 = 0 ;
                      bsAddData (Mrna, _bsRight, _Int, &x0) ;
                      bsAddData (Mrna, _bsRight, _Int, &x0) ;
                      bsPushObj (Mrna) ;
                      bsAddTag (Mrna, _Intron) ;
                      bsAddData (Mrna, _bsRight, _Text, "gt_ag") ;
                      bsGoto (Mrna, 0) ;
                    }
                  lasta2 = a2 ;
                  x1 = x2 + 1 ; x2 = x1 + a2 - a1 ; /* spliced coords */
                  bsAddData (Mrna, _Splicing, _Int, &a1) ;
                  bsAddData (Mrna, _bsRight, _Int, &a2) ;
                  bsAddData (Mrna, _bsRight, _Int, &x1) ;
                  bsAddData (Mrna, _bsRight, _Int, &x2) ;
                  bsPushObj (Mrna) ;
                  bsAddTag (Mrna, _Exon) ;
                  bsGoto (Mrna, 0) ;

                  { /* 5'UTR */
                    int v1, v2, b1, b2 ;
                    v1 = (x1 < p1 - 1 ? x1 : p1 - 1) ;
                    v2 = (x2 < p1 - 1 ? x2 : p1 - 1) ;
                    if (x1 <= p1 - 1)   /* intersect of spliced coords and 5' UTR */
                      {
                        b1 = a1 ; /* intersect of unspliced coords and 5' UTR */
                        b2 = a1 + v2 - v1  ;
                        bsAddData (Mrna, _Coding, _Int, &b1) ;
                        bsAddData (Mrna, _bsRight, _Int,&b2) ;        
                        bsAddData (Mrna, _bsRight, _Int, &v1) ;
                        bsAddData (Mrna, _bsRight, _Int, &v2) ;
                        bsAddData (Mrna, _bsRight, _Text, "5A") ;
                      }
                  }
                  { /* CDS */
                    int v1, v2, b1, b2 ;
                    v1 = (x1 > p1 ? x1 : p1) ;
                    v2 = (x2 < p2 ? x2 : p2) ;
                    if (v1 <= v2)   /* intersect of spliced coords and ORF */
                      {
                        b1 = x1 < v1 ? a1 + (v1 - x1) : a1 ; /* intersect of unspliced coords and ORF */
                        b2 = x2 > v2 ? a2 - (x2 - v2) : a2 ;
                        bsAddData (Mrna, _Coding, _Int, &b1) ;
                        bsAddData (Mrna, _bsRight, _Int,&b2) ;        
                        bsAddData (Mrna, _bsRight, _Int, &v1) ;
                        bsAddData (Mrna, _bsRight, _Int, &v2) ;
                        bsAddData (Mrna, _bsRight, _Text, "A") ;
                      }
                  }
                  { /* 3'UTR */
                    int v1, v2, b1, b2 ;
                    v1 = (x1 > p2+1 ? x1 : p2+1) ;
                    v2 = (x2 > p2+1 ? x2 : p2+1) ;
                    if (x2 > p2+1)   /* intersect of spliced coords and ORF */
                      {
                        b1 = x1 < v1 ? a1 + (v1 - x1) : a1 ; /* intersect of unspliced coords and ORF */
                        b2 = x2 > v2 ? a2 - (x2 - v2) : a2 ;
                        bsAddData (Mrna, _Coding, _Int, &b1) ;
                        bsAddData (Mrna, _bsRight, _Int,&b2) ;        
                        bsAddData (Mrna, _bsRight, _Int, &v1) ;
                        bsAddData (Mrna, _bsRight, _Int, &v2) ;
                        bsAddData (Mrna, _bsRight, _Text, "3A") ;
                      }
                  }
                }
            bsSave (Mrna) ;
          }

        /* locate the mrna on the genome */
        if (bsGetKey (Pg, _Source, &source) &&
            (Source = bsUpdate (source)))
          {
            if (bsFindKey (Source, _Subsequence, pg) &&
                bsGetData (Source, _bsRight, _Int, &a1) &&
                bsGetData (Source, _bsRight, _Int, &a2))
              {
                bsAddKey (Source, str2tag("mRNAs"), mrna) ;
                bsAddData  (Source, _bsRight, _Int, &a1) ;
                bsAddData  (Source, _bsRight, _Int, &a2) ;
              }
            bsSave (Source) ;

            if ((Mrna = bsUpdate (mrna)))
              {
                int ln = a2 - a1 + 1 ;
                
                if (ln < 0) ln = - ln ;
                bsAddData (Mrna, _Covers, _Int, &ln) ;
                bsAddData (Mrna, _bsRight, _Text,  "bp from") ;
                bsAddData (Mrna, _bsRight, _Int, &a1) ;
                bsAddData (Mrna, _bsRight, _Text, "to") ;
                bsAddData (Mrna, _bsRight, _Int, &a2) ;
                
                bsSave (Mrna) ;
                }
          }
      }

      { /* Product */
        if ((Product = bsUpdate (product)))
          {
            if (bsGetKey (Pg, str2tag ("Model_of_gene"), &key))
              bsAddKey (Product, str2tag("GeneBox"), key) ;
            bsAddKey (Product, str2tag("From_prediction"), pg) ;
            if (p1 > 3) bsAddTag (Product, _Complete) ;
            a1 = (met - leu)/3 + 1;
            if (leu)
              bsAddData (Product, _First_ATG, _Int, &a1) ;
            if (leu && a1 == 1)
              bsAddTag (Product, _At_position_1) ;
            pepL *= 3 ;
              bsAddTag (Product, _Best_product) ;

            bsAddData (Product
                       , isCoding ? _Coding_length : _Open_length
                       , _Int, &pepL) ;
            if (bsGetKey (Pg, _IntMap, &key) &&
                bsGetData (Pg, _bsRight, _Int, &a1) &&
                bsGetData (Pg, _bsRight, _Int, &a2))
              {
                if (a1 < a2) { a1 += p1 - 1 ; a2 = a1 + p2 - p1 ; }
                else { a1 -= (p1 - 1) ; a2 = a1  - (p2 - p1) ; }
                bsAddKey (Product, _IntMap, key) ;
                bsAddData (Product, _bsRight, _Int, &a1) ;
                bsAddData (Product, _bsRight, _Int, &a2) ;
              }
            /* coding and splicing */
            x1 = x2 = 0 ; lasta2 = 1 ;
            if (bsGetArray (Pg, _Source_Exons, units, 2))
              for (jj = 0 ; jj < arrayMax (units) ; jj += 2)
                {
                  uu = arrayp (units, jj, BSunit) ;
                  a1 = uu[0].i ; a2 = uu[1].i ; /* unspliced coords */
                
                  x1 = x2 + 1 ; x2 = x1 + a2 - a1 ; /* spliced coords */
                  {
                    int v1, v2, y1, y2 ;
                    v1 = x1 > p1 ? x1 : p1 ;
                    v2 = x2 < p2 ? x2 : p2 ;
                    if (v1 <= v2)   /* intersect of spliced coords and ORF */
                      {
                        /* 
			   int b1 = x1 < p1 ? a1 + (p1 - x1) : a1 ; // intersect of unspliced coords and ORF 
			   int b2 = x2 > p2 ? a2 - (x2 - p2) : a2 ;
			*/
                        y1 = v1 - (p1 - 1) ; y2 = v2 - (p1 - 1) ; /* spliced coords, origin = ORF */
                        bsAddData (Product,  _Source_Exons, _Int, &y1) ;
                        bsAddData (Product, _bsRight, _Int, &y2) ;
                      }
                  }
                }
            bsSave (Product) ;
          }
      }

      mrnaDesignSavePepAndKantor (product, pep) ;

      bsDestroy (Pg) ;
      arrayDestroy (dna) ;
      arrayDestroy (pep) ;
    }

  arrayDestroy (units) ;
  keySetDestroy (pgs) ;
}  /* mrnaTransferPg2PredictedMrna */

/*********************************************************************/

BOOL mrnaTransferRefseqMaker2Product (KEY am)
{
  int x, iOrf ;
  KEY cosmid, mrna, product = 0, tg, geneBox ;
  Array orfs = 0, mdna = 0, pep = 0 ;
  ORFS *orf ;
  OBJ Mrna = 0, Product = 0 ;
  KEY geneticCode = 0 ;
  char *translationTable ;
  BOOL isMrna5pComplete = FALSE ;


  tg = geneBox = 0 ;
  if ((mrna = keyGetKey (am, str2tag("RefSeqMaker"))) &&
      (product = keyGetKey (mrna, _Product)) &&
      (tg = keyGetKey (mrna, _From_gene)))
    geneBox = keyGetKey (tg, _Gene) ;
  if (!geneBox)
    return FALSE ;
  cosmid = keyGetKey (tg, str2tag("Genomic_sequence")) ;
  translationTable = pepGetTranslationTable (cosmid, &geneticCode) ;
  
  { /* dna */ 
    mdna = dnaGet (am) ; 
    if (!mdna) 
      return FALSE ;
  }
  
  isMrna5pComplete = keyFindTag (mrna, str2tag ("Found5p")) ;
  { /* orfs */ 
    orfs = mrnaRefseqMakerDna2Pep (cosmid, mdna, isMrna5pComplete) ;
    if (!orfs || !arrayMax(orfs))
      { arrayDestroy (mdna) ; return FALSE ; }        
    orf = arrp (orfs, 0, ORFS) ;
  }

  
  if ((Mrna = bsUpdate (mrna)))
    {
      bsAddTag (Mrna, _From_AM) ;
      bsRemoveTag (Mrna, _Product) ;
      bsSave (Mrna) ;
    }

  for (iOrf = 0 ; iOrf < arrayMax(orfs) ; iOrf++)
    {
      orf = arrp (orfs, iOrf, ORFS) ;
      if (! orf->isBest)
        continue ;
      {
        /* peptide */
        if (orf->cds > 0)
          {
            int j, nn = 0 ;
            char *cp, *cq ;
            
            j = orf->start - 1 ;  /* will start in frame */
            nn = orf->cds ;
            if (orf->start < 0)
              j += 3 ;
            if (0 && orf->downStop > 0) 
              nn += 3 ; /* do not export the stop */
            pep = arrayCreate (nn/3 + 4, char) ;
            array(pep, nn/3, char) = 0 ; arrayMax(pep)-- ;  /* make room, zero terminate */
            for (cp = arrp(pep, 0, char), cq = arrp (mdna, j, char) ; nn > 0 ; cp++, cq += 3, nn -= 3)
              *cp = pepEncodeChar[(int)e_codon (cq, translationTable)] ;
            *cp = 0 ; 
            if (orf->start == orf->leu) /* force a met */
              {
                cq = arrp (mdna, j, char) ; cp = arrp(pep, 0, char) ;
                *cp = pepEncodeChar[(int)'M'] ;
              }
          }        
        
        if (pep)
          {
            mrnaDesignSavePepAndKantor (product, pep) ;
            arrayDestroy (pep) ;
          }
        else
          { 
            arrayDestroy (orfs) ; arrayDestroy (mdna) ; 
            return FALSE ;
          }
      }        
  
      /* mRNA */
      if ((Mrna = bsUpdate (mrna)))
        {
          int xx ;
          
          bsAddKey (Mrna, _Product, product) ;
          xx = orf->start > 0 ?  orf->start :  orf->start + 3 ;
          bsAddData (Mrna, _bsRight, _Int, &xx) ;
          xx = orf->downStop > 0 ? orf->max + 3 : orf->max ; /* include the stop */
          bsAddData (Mrna, _bsRight, _Int, &xx) ;
          
          /* Product */
          if ((Product = bsUpdate (product)))
            {
              if (orf->isBest)
                bsAddTag (Product, _Best_product) ;

              x = orf->frame ;
              bsAddData (Product, _Frame, _Int, &x) ;
              
              bsRemoveTag (Product, str2tag("Completeness")) ;
              if (bsFindTag (Mrna, _Transpliced_to) ||
                  (bsFindTag (Mrna, str2tag ("Aggregated_5p_clones")) && orf->start > 3))
                bsAddTag (Product, str2tag("mRNA_5p_complete")) ;
              
            x = orf->nX ;
            if (x > 0)
              {
                int xx ;
                
                bsAddData (Product, str2tag("Coding_gap"), _Int, &x) ; 
                bsAddData (Product, _bsRight, _Text, "bp") ;
                xx = x/3 ;
                bsAddData (Product, _bsRight, _Int, &xx) ;
                bsAddData (Product, _bsRight, _Text, "X") ;  
              }

            /* up and down stops */
	    bsRemoveTag (Product, str2tag("Up_stop")) ;
            bsRemoveTag (Product, str2tag("Down_stop")) ;
            if (orf->upStop > 0)
              {
                int j, nn ;
                char *cp ;
                
                j = orf->start - 1 ;  /* start in frame */
                
                for (nn = 0, cp = arrp (mdna, j, char) ; j >= 0 ; cp -= 3, j -= 3)
                  if (e_codon (cp, translationTable) == '*')
                    {
                      x = j - (orf->start - 1) ;
                      if (!nn++)
                        bsAddTag (Product, str2tag("Up_stop")) ;
                      bsAddData (Product, _bsRight, _Int, &x) ; 
                      /* break ; to show just one, changed to show all of them june 30, danielle */
                    }
              }
            if (orf->downStop > 0)
              {
                x = orf->max - orf->start + 2 ;
                bsAddTag (Product, str2tag("Down_stop")) ;
                bsAddData (Product, _bsRight, _Int, &x) ; 
                /* show just one */
              }        
            }
          
          {
            Array units = arrayCreate (300, BSunit) ;
            Array vnits = arrayCreate (300, BSunit) ;
            BSunit *uu, *vv ;
            int jj, jk, u1, u2, v1, v2, w1, w2, x1, x2, x3 ;
            
            
            /* coding */
            x1 = orf->start ; x2 = orf->max ; x3 = arrayMax (mdna) ;
            if (bsGetArray (Mrna, _Splicing, units, 5))
              for (jj = jk = 0 ; jj < arrayMax(units) ; jj+= 5)
                {
                  uu = arrp(units, jj, BSunit) ;
                  if (!strstr (name(uu[4].k), "xon"))
                    continue ;
                  u1 = uu[0].i ; u2 = uu[1].i ; v1 = uu[2].i ; v2 = uu[3].i ;
                  /* the AM may me shorter than the original genome dna */
                  if (v1 > x3)
                    break ;
                  if (v2 > x3)
                    { u2 = u2 - v2 + x3 ; v2 = x3 ; }
                  w1 = x1 <= v2 ? x1 - 1 : v2 ;
                  w2 = x2 < v2 ? x2 : v2 ;
                  if (v2 < x1)
                    {
                      vv = arrayp (vnits, jk+4, BSunit) ; /* make room */
                      vv = arrp (vnits, jk, BSunit) ;   jk += 5 ;
                      vv[0].i = u1 ; vv[1].i = u2 ;
                      vv[2].i = v1 ; vv[3].i = v2 ;
                      vv[4].s = "5A" ;
                    }
                  else if (x2 < v1)
                    {
                      vv = arrayp (vnits, jk+4, BSunit) ; /* make room */
                      vv = arrp (vnits, jk, BSunit) ;   jk += 5 ;
                      vv[0].i = u1 ; vv[1].i = u2 ;
                      vv[2].i = v1 ; vv[3].i = v2 ;
                      vv[4].s = "3A" ;
                    }
                  else 
                    {
                      if (x1 > v1)
                        {
                          vv = arrayp (vnits, jk+4, BSunit) ; /* make room */
                          vv = arrp (vnits, jk, BSunit) ;   jk += 5 ;
                          vv[0].i = u1 ; vv[1].i = u2 - (v2 - x1 - 1) ;
                          vv[2].i = v1 ; vv[3].i = x1 - 1 ;
                          vv[4].s = "5A" ;
                        }
                      w1 = x1 < v1 ? v1 : x1 ;
                      w2 = x2 < v2 ? x2 : v2 ;
                      if (w1 < w2)
                        {
                          vv = arrayp (vnits, jk+4, BSunit) ; /* make room */
                          vv = arrp (vnits, jk, BSunit) ;   jk += 5 ;
                          vv[0].i = u1 + w1 - v1 ; vv[1].i = u2 - (v2 - w2) ;
                          vv[2].i = w1 ; vv[3].i = w2 ;
                          vv[4].s = "A" ;
                        }
                      if (w2 < v2)
                        {
                          vv = arrayp (vnits, jk+4, BSunit) ; /* make room */
                          vv = arrp (vnits, jk, BSunit) ;   jk += 5 ;
                          vv[0].i = u2 - (v2 - w2 - 1) ; vv[1].i = u2 ; 
                          vv[2].i = w2 + 1 ; vv[3].i = v2 ;
                          vv[4].s = "3A" ;
                        }
                    }
                }
            bsRemoveTag (Mrna, _Coding) ;
            bsAddArray (Mrna, _Coding, vnits, 5) ;
            
            arrayDestroy (units) ;
            arrayDestroy (vnits) ;
            
          }
          
          if (bsFindTag (Product, str2tag("NH2_Complete")) && bsFindTag (Product, str2tag("COOH_Complete")))
            bsAddTag (Product, str2tag("Complete")) ;
          else
            bsAddTag (Product, str2tag("Partial")) ;
          /* first met */
          bsSave (Product) ;
          bsSave (Mrna) ;
        }
    }
  
  dnaStoreDestroy (mrna, mdna) ; mdna = 0 ; /* replace the dna in the mrna object */
  arrayDestroy (orfs) ;
  /* mrnaAddFirstMetAndUorf (mrna) ; called by addkantorInfo */

  return TRUE ;
}  /* mrnaTransferRefseqMaker2Product */

/*********************************************************************/

int mrnaTransferRefseqMakerKeySet2Product (KEYSET ks)
{
  KEYSET ams = query (ks, "CLASS Sequence && RefSeqMaker") ; 
  int nn = 0 , iams = keySetMax (ams) ;

  for (iams = 0 ; iams < keySetMax (ams) ; iams++) 
    if (mrnaTransferRefseqMaker2Product (keySet (ams, iams)))
      nn++ ;
  keySetDestroy (ams) ;

  return nn ;
} /* mrnaTransferRefseqMakerKeyset2Product */

/*********************************************************************/

void mrnaAnalyseClusterAllPg (char *cp)
{
  KEYSET tgs = query (0, messprintf ("Find Predicted_gene %s", cp ? cp : "")) ;
  Array allTg = cDnaGetIntMap (tgs, 0) ;
  HIT *h1, *h2, *h3 ;
  int i, j, k, a1, a2 ;
  BOOL reverse = FALSE ;
  KEY map = 0, master = 0 ;
  Stack s = stackCreate (32000) ;
  BOOL same ;   
  Array ex2, ex3 ;
  
  ex2 = arrayCreate (256, BSunit) ;
  ex3 = arrayCreate (256, BSunit) ;

  for (i = 0 ; i <  arrayMax(allTg) ; i++) /* HORRIBLE HACK */
    {
      h1 = arrp (allTg, i, HIT) ; 
      master = h1->gene ;
      if (!strncmp (name(master), "Sp_", 3))
        catText (s, messprintf ("Sequence \"%s\"\nIdentifies \"%s\"\n\n",
                                name(master), name(master) + 3)) ;
    }


  for (i = 0 ; i <  arrayMax(allTg) ; i++)
    {
      h1 = arrp (allTg, i, HIT) ; 

      reverse = h1->reverse ;
      if (h1->cDNA_clone  || keyFindTag (h1->gene, str2tag("Bad_quality"))) continue ;
      a1 = h1->a1 ; a2 = h1->a2 ; map = h1->est ;
      master = h1->gene ;

      for (j = i ; j < arrayMax(allTg) ; j++)
        {
          h2 = arrp (allTg, j, HIT) ;
          if (h2->est != map) break ;
          if (a2 < h2->a1) break ;
          if (reverse != h2->reverse || keyFindTag (h2->gene, str2tag("Bad_quality")))
            continue ;
          if (a2 < h2->a2) a2 = h2->a2 ;
          same = FALSE ;
          for (k = i ; 0 && k < j ; k++) /* Ecoli identical coordinates yet count as distinct */
            { 
              h3 = arrp (allTg, k, HIT) ;
              if (h3->cDNA_clone != master)
                continue ;
              if (h2->a1 == h3->a1 &&
                  h2->a2 == h3->a2)
                {
                  OBJ G2 = 0, G3 = 0 ;
                  
                  G2 = bsCreate (h2->gene) ;
                  G3 = bsCreate (h3->gene) ;
                  if (G2 && G3 &&
                      bsGetArray (G2, _Source_Exons, ex2, 2) &&
                      bsGetArray (G3, _Source_Exons, ex3, 2) &&
                      arrayMax(ex3) == arrayMax(ex2))
                    {
                      int i3  = arrayMax (ex3) ;
                      same = TRUE ;
                      while (i3--)
                        if (arr (ex2,i3,BSunit).i != arr (ex3,i3,BSunit).i)
                          { same = FALSE ; break ; }                  
                    }
                  bsDestroy (G2) ;
                  bsDestroy (G3) ;
                }
              if (same)
                {
                  /* do this with bs to be real time, acedump may export to wrong object */
                  OBJ G2 = bsCreate (h2->gene) ;  
                  OBJ G3 = bsUpdate (h3->gene) ;
                  Array aa = arrayCreate (60, BSunit) ;

                  /* printf ("Eliminating %s in favor of %s", name(h2->gene), name (h3->gene)) ; */
                  h2->cDNA_clone = 1 ;
                  if (!strncmp (name(h2->gene), "Sp_", 3))
                    {
                      KEY mrna = keyGetKey (h2->gene, _DNA_homol) ;
                      if (mrna)
                        catText (s, messprintf ("Sequence \"%s\"\nIdentifies \"%s\"\n\n", name (h3->gene), name(mrna))) ;
                    }

                  if (bsGetArray (G2, _Homol, aa, 2))
                    {
                      int iaa ;
                      BSunit *uu ;
                      for (iaa = 0 ; iaa < arrayMax (aa) ; iaa += 2)
                        {
                          uu = arrp (aa, iaa, BSunit) ;
                          bsAddKey (G3, uu[0].k, uu[1].k) ;
                        }
                    }

                  if (bsGetArray (G2, str2tag("Identifies"), aa, 1))
                    {
                      int iaa ;
                      BSunit *uu ;
                      for (iaa = 0 ; iaa < arrayMax (aa) ; iaa += 1)
                        {
                          uu = arrp (aa, iaa, BSunit) ;
                          bsAddKey (G3, str2tag("Identifies"), uu[0].k) ;
                        }
                    }

                  bsDestroy (G2) ;
                  bsSave (G3) ;
                  arrayDestroy (aa) ;
                  break ;  /* the k loop */
                }
            }  /* end of k loop */              
          if (same)
            continue ; /* the j loop */
          h2->cDNA_clone = master ;
          catText (s, messprintf ("Sequence \"%s\"\nModel_of \"X_%s\"\n\n", name(h2->gene), name(master))) ;
          if (keyFindTag (h2->gene, _Locus))
            {
              KEYSET Loci = queryKey (h2->gene, ">Locus") ;
              int iLoci = keySetMax(Loci) ;
              
              for (iLoci = 0 ; iLoci < keySetMax(Loci) ;  iLoci++)
                catText (s, messprintf ("Sequence \"X_%s\"\nLocus \"%s\"\n\n", name(master), name(keySet(Loci, iLoci)))) ;
              
              keySetDestroy (Loci) ;
            }
        }    /* end of j loop */
      /* should make it unique new name */
      if (reverse) { int tmp = a1 ; a1 = a2 ; a2 = tmp ; }
      catText (s, messprintf ("Sequence \"%s\"\nSubsequence \"X_%s\" %d %d \n\n", name (map), name(master), a1, a2)) ;
      catText (s, messprintf ("Sequence \"X_%s\"\nIntMap \"%s\" %d %d \n\n", name (master), name(map), a1, a2)) ;
    }
  
  parseBuffer (stackText (s, 0) , 0) ; 
  stackDestroy (s) ;
  keySetDestroy (tgs) ;
  arrayDestroy (allTg) ;
  arrayDestroy (ex2) ;
  arrayDestroy (ex3) ;
  return ;
} /* analyseClusterAllPg */

/*********************************************************************/
/*********************************************************************/
/* get name of smallest Clone Group or smallest clone */
static void mrnaSelectCosmid (SC* sc)
{
  KEYSET ks = 0 ;
  KEY map, part, bestCosmid = 0, cosmid = sc->s2m->cosmid ;
  OBJ Cosmid = 0 ;
  int ii, a1 = 0, a2 = 0, b1, b2, c1 = 0, c2 = 0, g1 = 0, g2 = 0 ;

  /* construct g1, g2 the IntMap coords of the tgene */
  if (!cosmid)
    return ;

  if ((Cosmid = bsCreate (cosmid)))
    {
      if (bsGetKey (Cosmid, _IntMap, &map) &&
          bsGetData (Cosmid, _bsRight, _Int, &c1) &&
          bsGetData (Cosmid, _bsRight, _Int, &c2)) 
	{} ;
      bsDestroy (Cosmid) ;
    }
  if (c1 && c2)
    {
      if (c1 < c2) { g1 = c1 + sc->a1 - 1 ; g2 = c1 + sc->a2 - 1 ; }
      else { g1 = c1 - sc->a1 + 1 ; g2 = c1 - sc->a2 + 1 ; }
    }

  /* i try to see if my genes fits inside one of the parts */
  bestCosmid = 0 ;
  if (g1 && g2)
    {
      ks = queryKey (cosmid, ">Parts") ;
      for (ii = 0 ; !bestCosmid && ii < keySetMax (ks) ; ii++)
        {
          part = keySet (ks, ii) ;
          if ((Cosmid = bsCreate (part)))
            {
              if (bsGetKey (Cosmid, _IntMap, &map) &&
                  bsGetData (Cosmid, _bsRight, _Int, &b1) &&
                  bsGetData (Cosmid, _bsRight, _Int, &b2))
		{} ;
              bsDestroy (Cosmid) ;
            }
          if (b1 < b2 && b1 <= g1 && b1 <= g2 && b2 >= g1 && b2 >= g2)
            { a1 = g1 - b1 + 1 ; a2 = g2 - b1 + 1 ; bestCosmid = part ; }
          if (b1 > b2 && b2 <= g1 && b2 <= g2 && b1 >= g1 && b1 >= g2)
            { a1 = b1 - g1 + 1 ; a2 = b1 - g2 + 1 ; bestCosmid = part ; }
        }
      keySetDestroy (ks) ;
    }
  sc->d1 = sc->a1 ;
  sc->d2 = sc->a2 ;

  if (bestCosmid)
    { 
      sc->s2m->cosmid = bestCosmid ;
      sc->a1 = a1 ;
      sc->a2 = a2 ; 
    }
} /* mrnaSelectBestCosmid */

/*********************************************************************/
/* get name of smallest Clone Group or smallest clone */
static KEY  mrnaNameGeneByCloneGroup (KEYSET ks)
{
  KEYSET ks2 = 0, ks3 = 0 ;
  KEY gg = 0 ;
  
  ks2 = query (ks, "FOLLOW Clone_group Y*") ;
  if (keySetMax(ks2)) 
    {
      ks3 = keySetAlphaHeap (ks2, 1) ;
      gg = keySet (ks3, 0) ;
    }
  else
    { 
      keySetDestroy (ks2) ;
      ks2 = query(ks,"IS yk*") ;
      if (keySetMax(ks2)) 
        {
          ks3 = keySetAlphaHeap (ks2, 1) ;
          gg = keySet (ks3, 0) ;
        }
      else
        {
          keySetDestroy (ks2) ;
          ks2 = query(ks,"IS y*") ;
          if (keySetMax(ks2)) 
            {
              ks3 = keySetAlphaHeap (ks2, 1) ;
              gg = keySet (ks3, 0) ;
            }
          else
            {
              keySetDestroy (ks2) ;
              ks2 = query(ks,"IS mv*") ;
              if (keySetMax(ks2)) 
                {
                  ks3 = keySetAlphaHeap (ks2, 1) ;
                  gg = keySet (ks3, 0) ;
                }
              else
                {
                  keySetDestroy (ks2) ;
                  ks2 = query(ks,"IS cm*") ;
                  if (keySetMax(ks2)) 
                    {
                      ks3 = keySetAlphaHeap (ks2, 1) ;
                      gg = keySet (ks3, 0) ;
                    }
                  else
                    {
                      keySetDestroy (ks2) ;
                      ks2 = query(ks,"IS GB*") ;
                      if (keySetMax(ks2)) 
                        {
                          ks3 = keySetAlphaHeap (ks2, 1) ;
                          gg = keySet (ks3, 0) ;
                        }
                      else
                        {
                          keySetDestroy (ks2) ;
                          ks2 = query(ks,"IS GB*") ;
                          if (keySetMax(ks2)) 
                            {
                              ks3 = keySetAlphaHeap (ks2, 1) ;
                              gg = keySet (ks3, 0) ;
                            }
                          else
                            {
                              keySetDestroy (ks2) ;
                              ks2 = query(ks,"IS *") ;
                              if (keySetMax(ks2)) 
                                {
                                  ks3 = keySetAlphaHeap (ks2, 1) ;
                                  gg = keySet (ks3, 0) ;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

  keySetDestroy (ks2) ;  
  keySetDestroy (ks3) ;  

  return gg ;
}   /* mrnaNameGeneByCloneGroup */

/*********************************************************************/
/* get name of smallest Clone Group or smallest clone */
static KEY mrnaNameGeneByPhenotypicGene (SC *sc)
{
  KEYSET ks = 0, ks1 = 0, ks2 = 0, ks3 = 0 ;
  KEY cosmid = sc->s2m->cosmid1 ? sc->s2m->cosmid1 : sc->s2m->cosmid ;
  OBJ Cosmid = 0, Gene = 0 ;
  KEY gene, map = 0, goodGene = 0 ;
  int ii, b1, b2, d1, d2, dg, db, dx, c1 = 0, c2 = 0, g1 = 0, g2 = 0, 
    bestii, bestDx = 0, bestMissmatch = 0 ;

  /* construct g1, g2 the IntMap coords of the tgene */
  if (!cosmid)
    return 0 ;

  if ((Cosmid = bsCreate (cosmid)))
    {
      if (bsGetKey (Cosmid, _IntMap, &map) &&
          bsGetData (Cosmid, _bsRight, _Int, &c1) &&
          bsGetData (Cosmid, _bsRight, _Int, &c2))
	{} ;
      bsDestroy (Cosmid) ;
    }
  if (c1 && c2)
    {
      if (c1 < c2) { g1 = c1 + sc->a1 - 1 ; g2 = c1 + sc->a2 - 1 ; }
      else { g1 = c1 - sc->a1 + 1 ; g2 = c1 - sc->a2 + 1 ; }
    }

  if (!g1 && !g2)
    return 0 ;


  ks1 = queryKey (cosmid, "{IS *} $| {>Parts} $| {>Subsequence} ") ; /* 1 cosmid OR 1 junction + 2 cosmids */
  ks2 = query (ks1, "{IS *} $| {>Source} $| {>In_junction} ") ; /* 3 junctions and 2 cosmids OR 2 junctions + 2 cosmids */
  ks3 = query (ks2, "{IS *} $| {>Source} $| {>Parts}") ;  /* 2 junctions and 3 cosmids OR 3 junctions + 4 cosmids */
  ks = query (ks3, ">Genes") ;

  bestii = -1 ;
  for (ii = 0 ; !goodGene && ii < keySetMax (ks) ; ii++)
    {
      gene = keySet (ks, ii) ;
      if ((Gene = bsCreate (gene)))
        {
          if (bsFindKey (Gene, _IntMap, map) &&
              bsGetData (Gene, _bsRight, _Int, &b1) &&
              bsGetData (Gene, _bsRight, _Int, &b2))
            {
              if (g1 < g2 && b1 < b2 && b1 < g2 && b2 > g1 )
                {
                  d1 = g1 > b1 ? g1 : b1 ;
                  d2 = g2 < b2 ? g2 : b2 ;
                  dg = g2 - g1 ;
                  db = b2 - b1 ;
                  dx = d2 - d1 ;
                  if (bestii == -1 ||
                      dx > bestDx + 20 ||
                      (dx > bestDx - 100 &&  bestMissmatch > dg + db - 2 * dx + dx - bestDx)
                      )
                    { bestDx = dx ; bestii = ii ; bestMissmatch = dg + db - 2 * dx ; }
                }
              else if ( g1 > g2 && b1 > b2 && b1 > g2 && b2 < g1 ) 
                {
                  d1 = g1 < b1 ? g1 : b1 ;
                  d2 = g2 > b2 ? g2 : b2 ;
                  dg = g1 - g2 ;
                  db = b1 - b2 ;
                  dx = d1 - d2 ;
                  if (bestii == -1 ||
                      dx > bestDx + 20 ||
                      (dx > bestDx - 100 &&  bestMissmatch > dg + db - 2 * dx + dx - bestDx)
                      )
                    { bestDx = dx ; bestii = ii ; bestMissmatch = dg + db - 2 * dx ; }
                }
            }
          bsDestroy (Gene) ;
        }
    }
  
  if (bestii >= 0)
    {
      gene = keySet (ks, bestii) ;
      if ((Gene = bsCreate (gene)))
        {
          if (bsFindKey (Gene, _IntMap, map) &&
              bsGetData (Gene, _bsRight, _Int, &b1) &&
              bsGetData (Gene, _bsRight, _Int, &b2))
            {
              if (g1 < g2 && b1 < b2 && b1 < g2 && b2 > g1 )
                {
                  d1 = g1 > b1 ? g1 : b1 ;
                  d2 = g2 < b2 ? g2 : b2 ;
                  dg = g2 - g1 ;
                  db = b2 - b1 ;
                  dx = d2 - d1 ;
                  if (3*dx > dg && 3*dx > db)   /* see mrnaAnalyseGenefinder */
                    lexaddkey (name(gene), &goodGene, _VTranscribed_gene) ;
                }
              else if ( g1 > g2 && b1 > b2 && b1 > g2 && b2 < g1 ) 
                {
                  d1 = g1 < b1 ? g1 : b1 ;
                  d2 = g2 > b2 ? g2 : b2 ;
                  dg = g1 - g2 ;
                  db = b1 - b2 ;
                  dx = d1 - d2 ;
                  if (3*dx > dg && 3*dx > db)   /* see mrnaAnalyseGenefinder */
                    lexaddkey (name(gene), &goodGene, _VTranscribed_gene) ;
                }
            }
          bsDestroy (Gene) ;
        }
    }
  
  keySetDestroy (ks) ;
  keySetDestroy (ks1) ;
  keySetDestroy (ks2) ;
  keySetDestroy (ks3) ;

  return goodGene ;
}

/*********************************************************************/
/* get name of smallest Clone Group or smallest clone */
static KEY mrnaNameGeneBySection (SC* sc)
{
  KEY gg, cosmid = sc->s2m->cosmid ;
  int a1, a2 ;
  int x, block = 100 ; /* attention, 1000 was too big producing name conflicts */

  /* construct g1, g2 the IntMap coords of the tgene */

  a1 = sc->a1 ; a2 = sc->a2 ; 
      
  if (a1 < a2)
    {
      x = a2/block ;
      x |= 1 ; /* should be odd */
      if (block * x > a2) x -= 2 ;
      if (block * x < a1 && 
          a2 - block * x > block * (x + 2) - a2)
        x += 2 ;
    }
  else
    {
      x = a2/block ;
      x &= ~1 ; /* should be even */
      while (block * x < a2) x += 2 ;
      if (block * x > a1 && 
          block * x - a2 > a2 - block * (x - 2))
        x -= 2 ;
    }
  
  lexaddkey (messprintf("G_%s_%d", name(cosmid), x), &gg, _VTranscribed_gene) ;
  return gg ;
}

/*********************************************************************/

KEY mrnaGetGeneName (SC *sc)
{
  int n ;
  KEY cGroup = 0, gene = 0 ;
  int type = 0 ;

  /* first find a candidate gene name */
  switch (mrnaPleaseNameBy ())
    { 
    case 1:    /* human genome */
      type = 1 ;
      cGroup = mrnaNameGeneBySection (sc) ;
      break ;
    case 2:    /* worm genome name by newname-gene */
      type = 1 ;
      cGroup = mrnaNameGeneByPhenotypicGene (sc) ;        
      if (cGroup)
        break ;
      /* else fall thru */
    default:   /* worm case */
      type = 0 ;
      gene = mrnaNameGeneByCloneGroup (sc->sh.clones) ;
      if (class(gene) == _VClone_Group)
        cGroup = gene ;
      else 
        lexaddkey (messprintf("G_%s", name(gene)), &cGroup, _VClone_Group) ;
      sc->cGroup = cGroup ;
    }
  
  /* iterate to avoid mixing doubles */
  for (n = 0 ; TRUE ; n++)       /* iterate on gene_name candidate */
    {
      BOOL isNew = TRUE ;
      
      switch (type)
        { 
        case 1:    /* name by section */
          if (n == 0)
            lexaddkey (messprintf("%s", name(cGroup)), &gene, _VTranscribed_gene) ;
          else
            lexaddkey (messprintf("%s_%d", name(cGroup), n), &gene, _VTranscribed_gene) ;
          break ;
        default:   /* worm case */
          if (n == 0)
            {
              if (strncmp("G_",name(cGroup),2))
                lexaddkey (messprintf("G_%s", name(cGroup)), &gene, _VTranscribed_gene) ;
              else
                lexaddkey (messprintf("%s", name(cGroup)), &gene, _VTranscribed_gene) ;
            }
          else
            { 
              if (strncmp("G_",name(cGroup),2))
                lexaddkey (messprintf("G_%s_%d", name(cGroup), n), &gene, _VTranscribed_gene) ;
              else
                lexaddkey (messprintf("%s_%d", name(cGroup), n), &gene, _VTranscribed_gene) ;
            }
          break ;
        }
      if (keyFindTag (gene, _Splicing)) /* name already in use */
        {
          /* printf(" spl:%s(%d/%d) ", name(gene), i, j1) ; */
          isNew = FALSE ;
        } 
      /* printf(" ok:%s(j1=%d) ", name(gene), j1) ; */
      /* success this gene name is good */
      if (isNew)
        break ;
    }
  
  return gene ;
}

/*****************************************/
/* this is used to assign the gene to the cosmid rather than to the link */
static KEY  mrnaRelevantCosmid (KEY cosmid, KEY cosmid1, 
                                Array linkPos, int x1, int x2) 
{
  KEY relevantCosmid = cosmid ;
  HIT *hh ;
  
  if (cosmid1 && linkPos &&
      (hh = arrp(linkPos, arrayMax(linkPos) -1, HIT)) &&
      (hh->a1 < hh->a2 && hh->x1 < hh->x2 && hh->a1 < hh->x2) &&
      (x1 < hh->a2 && x2 < hh->a2))
    relevantCosmid = cosmid1 ; /* i.e. the first cosmid of the link */

  return relevantCosmid ;
}

/*********************************************************************/

int mrnaQualityEvaluate (int ln, int ali, int  nerr, BOOL isMrna, int *ixp, int *iyp)
{
  float x = 0.0, y = 100.0, ln2 = 0.0 ;
  int ix, iy, qual ;

  /* x = percent alignement */
  ln2 = ln + .1 ;
  if (!isMrna && ln < 1500)
    {
      if (ln2 > 600.0) ln2 = 600.0 ;
      if (ln2 < (float)ali) ln2 = (float) ali ;
    }

  if (ln2 > 2)                   /* x = percent aligned */
    x = (100.0 * (float)ali) / ln2 ;
  if (ali > 2)                   /* y = percent error sur zone alignee */
    y =  (100.0 * (float)nerr) / ali ;
   
  if (x > 98.0)
    {
      ix = 1 ;
      if (y < .1)
        {
          iy = 1 ;
          qual = isMrna ? 1 : 1 ;
        }
      else if (y < 1.0)
        {
          iy = 2 ;
          qual = isMrna ? 2 : 2 ;
        }
      else if (y < 2.0)
        {
          iy = 3 ;
          qual = isMrna ? 3 : 3 ;
        }
      else if (y < 3.0)
        {
          iy = 4 ;
          qual = isMrna ? 4 : 4 ;
        }
      else if (y < 4.0)
        {
          iy = 5 ;
          qual = isMrna ? 5 : 5 ;
        }
      else if (y < 5.0)
        {
          iy = 6 ;
          qual = isMrna ? 7 : 6 ;
        }
      else
        {
          iy = 7 ;
          qual = isMrna ? 9 : 7 ;
        }
    }
  else if (x > 90.0)
    { 
      ix = 2 ;
      if (y < .1)
        {
          iy = 1 ;
          qual = isMrna ? 2 : 2 ;
        }
      else if (y < 1.0)
        {
          iy = 2 ;
          qual = isMrna ? 3 : 3 ;
        }
      else if (y < 2.0)
        {
          iy = 3 ;
          qual = isMrna ? 4 : 4 ;
        }
      else if (y < 3.0)
        {
          iy = 4 ;
          qual = isMrna ? 5 : 5 ;
        }
      else if (y < 4.0)
        {
          iy = 5 ;
          qual = isMrna ? 6 : 6 ;
        }
      else if (y < 5.0)
        {
          iy = 6 ;
          qual = isMrna ? 8 : 7 ;
        }
      else
        {
          iy = 7 ;
          qual = isMrna ? 10 : 8 ;
        }
       
    }
  else if (x > 80.0)
    {
      ix = 3 ;
      if (y < .1)
        {
          iy = 1 ;
          qual = isMrna ? 3 : 3 ;
        }
      else if (y < 1.0)
        {
          iy = 2 ;
          qual = isMrna ? 4 : 4 ;
        }
      else if (y < 2.0)
        {
          iy = 3 ;
          qual = isMrna ? 5 : 5 ;
        }
      else if (y < 3.0)
        {
          iy = 4 ;
          qual = isMrna ? 6 : 6 ;
        }
      else if (y < 4.0)
        {
          iy = 5 ;
          qual = isMrna ? 7 : 7 ;
        }
      else if (y < 5.0)
        {
          iy = 6 ;
          qual = isMrna ? 9 : 8 ;
        }
      else
        {
          iy = 7 ;
          qual = isMrna ? 10 : 9 ;
        }
    }
  else if (x > 50.0 || ali > 800)
    {
      ix = 4 ;
      if (y < .1)
        {
          iy = 1 ;
          qual = isMrna ? 4 : 5 ;
        }
      else if (y < 1.0)
        {
          iy = 2 ;
          qual = isMrna ? 5 : 6 ;
        }
      else if (y < 2.0)
        {
          iy = 3 ;
          qual = isMrna ? 6 : 7 ;
        }
      else if (y < 3.0)
        {
          iy = 4 ;
          qual = isMrna ? 7 : 9 ;
        }
      else if (y < 4.0)
        {
          iy = 5 ;
          qual = isMrna ? 9 : 10 ;
        }
      else if (y < 5.0)
        {
          iy = 6 ;
          qual = isMrna ? 10 : 11 ;
        }
      else
        {
          iy = 7 ;
          qual = isMrna ? 11 : 11 ;
        }
       
    }
  else
    {
      ix = 5 ;
      if (y < .1)
        {
          iy = 1 ;
          qual = isMrna ? 7 : 9 ;
        }
      else if (y < 1.0)
        {
          iy = 2 ;
          qual = isMrna ? 8 : 10 ;
        }
      else if (y < 2.0)
        {
          iy = 3 ;
          qual = isMrna ? 10 : 11 ;
        }
      else if (y < 3.0)
        {
          iy = 4 ;
          qual = isMrna ? 11 : 11 ;
        }
      else if (y < 4.0)
        {
          iy = 5 ;
          qual = isMrna ? 11 : 11 ;
        }
      else if (y < 5.0)
        {
          iy = 6 ;
          qual = isMrna ? 11 : 11 ;
        }
      else
        {
          iy = 7 ;
          qual = isMrna ? 11 : 11 ;
        }
    }
  if (ixp) *ixp = ix ;
  if (iyp) *iyp = iy ;
  return qual ;
}

/*****/

static void mrnaSaveMatchQuality (KEY gene, Array hits, KEYSET ignoredClones)
{
  KEY cDNA_clone, est = 0 ;
  int j1, j2, j3, nmatch, nerr, nerrAll, nReads, quality, lengthEst, a1, a2, b1, b2, iGhits = 0 ;
  int dx, minx, maxx, clipTop = 0, offset = 0, maxExact ;
  OBJ Gene = 0,Est = 0, Clone = 0 ; 
  HIT *ah, *gh ;
  KEYSET ks = 0 ;
  BOOL isMrna = FALSE ;
  static KEY _Ref_mRNA = 0 ;
  Array gHits = 0 ;

  if (! _Ref_mRNA)
    _Ref_mRNA = str2tag ("Ref_mRNA") ;
  
  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;
  ks = keySetCreate () ;
  
  Gene = bsUpdate(gene) ;
  bsRemoveTag (Gene, _Read) ;
  bsSave(Gene) ;  /* save before adding,  Est saving in loop */

  iGhits = 0 ;
  gHits = arrayCreate (arrayMax(hits), HIT) ;

  for (j1 = 0 ; j1 < arrayMax(hits) ; j1++)
    {
      ah = arrp(hits, j1, HIT) ;
      if (!ah->cDNA_clone || ah->est == 1 || !(ah->type & gX) ||
          (ignoredClones && keySetFind (ignoredClones, ah->cDNA_clone, 0))) 
        continue ;        
          
      cDNA_clone = ah->cDNA_clone ;
      nReads = 0 ;
      est = 0 ;
      a1 = a2 = b1 = b2 = 0 ; /* extremities of 3 and 5p reads */
      maxExact = 0 ;
      for (j2 = 0, ah = arrp(hits, j2, HIT); j2 < arrayMax(hits); j2++, ah++)
        {
          if ((ah->type & gX) && cDNA_clone == ah->cDNA_clone)        
            { 
              nReads |= ah->reverse ? 0x1 : 0x2 ;
              if (ah->x1 < ah->x2)
                { 
                  if (!a1 && !a2)
                    { a1 = ah->a1 ; a2 = ah->a2 ; }
                  if (a1 > ah->a1) a1 = ah->a1 ;
                  if (a2 < ah->a2) a2 = ah->a2 ;
                }
              if (ah->x1 > ah->x2)
                { 
                  if (!b1 && !b2)
                    { b1 = ah->a1 ; b2 = ah->a2 ; }
                  if (b1 > ah->a1) b1 = ah->a1 ;
                  if (b2 < ah->a2) b2 = ah->a2 ;
                }
                          
            }
        }
      for (j2 = j1 ; j2 < arrayMax(hits); j2++)
        {
          ah = arrp(hits, j2, HIT) ;
          if (!(ah->type & gX) || ah->cDNA_clone != cDNA_clone)
            continue ;
          est = ah->est ;
          if (keySetFind (ks, est, 0))
            continue ;
                  
          keySetInsert (ks, est) ;
          nmatch = maxExact = 0 ;
          nerr = nerrAll = 0 ;
          minx = maxx = 0 ;
          clipTop = ah->clipTop ;
          lengthEst = ah->clipEnd - clipTop + 1 ;
          for (j3 = j2 ; j3 < arrayMax(hits); j3++)
            {
              ah = arrp(hits, j3, HIT) ;
              if (!(ah->type & gX) || est != ah->est)
                continue ;
                          
              if (minx == maxx)
                {  minx = maxx = ah->x1 ; }
              if (ah->x1 < minx) minx = ah->x1 ;
              if (ah->x2 < minx) minx = ah->x2 ;
              if (ah->x1 > maxx) maxx = ah->x1 ;
              if (ah->x2 > maxx) maxx = ah->x2 ;
	      if (maxExact < ah->maxExact) maxExact = ah->maxExact ;
       
              nerr += ah->nerr ;
              nerrAll += ah->nerrAll ;
              dx = ah->x2 > ah->x1 ? ah->x2 - ah->x1 + 1 : ah->x1 - ah->x2 + 1 ;
              nmatch += dx ;
            }
                  
          offset = minx - clipTop ;
          if (b1 < b2 && a1 < a2 && b1 < a2)
            {
              if ((Clone = bsUpdate (cDNA_clone)))
                {
                  bsAddKey (Clone, str2tag("Fully_sequenced"), gene) ;
                  bsSave (Clone) ;
                }
              lengthEst = nmatch + offset ;
            }
                  
          isMrna = keyFindTag (est, _mRNA) || keyFindTag (est, _Ref_mRNA)  ;
          quality = mrnaQualityEvaluate (lengthEst, nmatch, nerr, isMrna, 0, 0) ;
         
          gh = arrayp (gHits, iGhits++, HIT) ;          
          gh->est = est ;
          gh->a1 = minx ;  /* was lengthEst */
          gh->a2 = maxx ;  /* was offset */
          gh->x1 = nmatch ;
          gh->nerrAll = nerrAll ;
          gh->x2 = quality ;
	  gh->maxExact = maxExact ;
          
          if ((Est = bsUpdate (est)))
            {
              bsAddKey (Est, _From_gene, gene) ;
              bsAddData (Est, _bsRight, _Int, &minx) ; /* was lengthEst */
	      bsAddData (Est, _bsRight, _Int, &maxx) ; /* was offset */
              bsAddData (Est, _bsRight, _Int, &nmatch) ;
              bsAddData (Est, _bsRight, _Int, &nerrAll) ;
              bsAddData (Est, _bsRight, _Int, &quality) ;
              bsAddData (Est, _bsRight, _Int, &maxExact) ;
              bsAddData (Est, _bsRight, _Int, &offset) ;
              bsAddData (Est, _bsRight, _Int, &lengthEst) ; /* to align */
              bsSave (Est) ;
            }
        }
    }
  keySetDestroy (ks) ;

  /* i must save at the end to avoid copying the Gene in cache2 several thousand times */
  Gene = bsUpdate(gene) ;
  if (arrayMax(gHits))
    for (j1 = 0, gh = arrp (gHits, 0, HIT) ; j1 < iGhits ; gh++, j1++)
      {
        bsAddKey (Gene, _Read, gh->est) ;
        bsAddData (Gene, _bsRight, _Int, &(gh->a1)) ;
        bsAddData (Gene, _bsRight, _Int, &(gh->a2)) ;
        bsAddData (Gene, _bsRight, _Int, &(gh->x1)) ;
        bsAddData (Gene, _bsRight, _Int, &(gh->nerrAll)) ;
        bsAddData (Gene, _bsRight, _Int, &(gh->x2)) ;
        bsAddData (Gene, _bsRight, _Int, &(gh->maxExact)) ;
      } 
  bsSave(Gene) ;
  arrayDestroy (gHits) ;
} /* mrnaSaveMatchQuality */

/*********************************************************************/
/* save confirmed introns boundaries as gt_ag */
KEY mrnaIntronBoundary (Array dnaD, Array dnaR,  int a1, int a2)
{
  KEY fKey = 0 ;
  char *cp, *cq, feet[6] ;
  Array dna = 0 ;
  int x1, x2 ;

  feet[2] = '_' ; feet[5] = 0 ;

  if (!dnaD  || !dnaR)
    return 0 ;

  if (a1 < a2)
    { x1 = a1 - 1 ; x2 = a2 - 1 ; dna = dnaD ;}
  else
    { x1 = arrayMax(dnaR) - a1 ; x2 = arrayMax(dnaR) - a2 ; dna = dnaR ;}

  if (x1 < 0 || x1 > x2 || x2 >= arrayMax(dna))
    return 0 ;

  cp = arrp (dna, x1, char) ;
  cq = arrp (dna, x2 - 1, char) ;
  if ((*cp & G_) && (*(cp+1) & T_) &&
      (*cq & A_) && (*(cq+1) & G_))
    { lexaddkey ("gt_ag", &fKey, 0) ; }
  else if ((*cp & G_) && (*(cp+1) & C_) &&
           (*cq & A_) && (*(cq+1) & G_))
    { lexaddkey ("gc_ag", &fKey, 0) ; }
  else
    {
      feet[0] = dnaDecodeChar [(int)(*cp)] ;
      feet[1] = dnaDecodeChar [(int)(*(cp+1))] ;

      feet[3] = dnaDecodeChar [(int)(*cq)] ;
      feet[4] = dnaDecodeChar [(int)(*(cq+1))] ;
      
      fKey = _Other ;
      lexaddkey (feet, &fKey, 0) ;
    }
  return fKey ;
}

/********************************************************************************/
/* clean up all data associated to the gene */
void mrnaCleanUp (KEY cosmid, KEY gene, KEYSET genes)
{
  KEYSET 
    ks = 0, ks01 = 0, ks02 = 0, ks1 = 0, ks2 = 0, ks3 = 0, ks4 = 0;
  OBJ Gene = 0 ;
  int i ;
  if (cosmid)
    {
      ks = queryKey (cosmid, ">Transcribed_gene") ;
      ks01 = queryKey (cosmid, ">mRNA") ;
      ks02 = query (ks, ">mRNA") ;
      ks1 = keySetOR (ks01, ks02) ;
    }
  else if (genes)
    {
      ks = keySetCopy (genes) ;
      ks1 = queryKey (gene, ">mRNA") ;
    }
  else if (gene)
    {
      ks = keySetCreate () ;
      keySet(ks, 0) = gene ;
      ks1 = queryKey (gene, ">mRNA") ;
    }
  
  ks2 = query (ks1,  ">DNA") ;
  ks3 = query (ks1,  ">Product") ;
  ks4 = query (ks3,  ">Peptide") ;


  /* keySetKill(ks) ; done in cdnaCleanUp */
  keySetKill(ks1) ;
  keySetKill(ks2) ;
  keySetKill(ks3) ;
  keySetKill(ks4) ;

  for (i = 0 ; ks && i < keySetMax (ks) ; i++)
    {
      gene = keySet (ks, i) ;
      if (gene && (Gene = bsUpdate (gene)))
        {
          /* printf (" kill:%s ", name(gene) ) ; */
          bsRemoveTag (Gene, _Splicing) ;
          bsSave (Gene) ;
        }
    }

  keySetDestroy (ks) ;
  keySetDestroy (ks01) ;
  keySetDestroy (ks02) ;
  keySetDestroy (ks1) ;
  keySetDestroy (ks2) ;
  keySetDestroy (ks3) ;
  keySetDestroy (ks4) ;

  return ;
} /* mrnaCleanUp */

/*********************************************************************/

static int mrnaAddPrediction (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas, Array linkPos)
{
  int ii, iks, j, a1, a2, b1, b2, c1, c2, oldMaxSmrnas = arrayMax(smrnas), ml1, ml2 ;
  SMRNA *mrna ;
  HIT *up ;
  Array predictedHits = 0 ;
  KEY pg = 0, map = 0, mappg = 0 ;
  OBJ PG = 0 ;
  Array units = 0 ;
  BSunit *uu, *uu2 ;
  KEYSET ks = queryKey (s2m->cosmid, "{> Subsequence} $| {>Parts ; >Subsequence} $| {>Source ; >Subsequence} $| {>genes ; >genefinder} $|  {>Parts ;>genes ; >genefinder} $|  {>Source ;>genes ; >genefinder} ; Source_Exons" ) ;

  a1 = sc->a1 ; a2 = sc->a2 ;
  for (iks = 0 ; iks < keySetMax(ks) ; iks++)
    {
      pg = keySet (ks, iks) ;
      if ((PG = bsCreate (pg)))
        {
          units = arrayReCreate (units, 40, BSunit) ;
          predictedHits = arrayCreate (24, HIT) ;
          c1 = c2 = 0 ;
          map = mappg = 0 ;
          if (bsGetKey (PG, _IntMap, &mappg) &&
              bsGetData (PG, _bsRight, _Int, &c1)  &&
              bsGetData (PG, _bsRight, _Int, &c2))
            {
              OBJ Link = bsCreate (s2m->cosmid) ;
              
              if (Link)
                {
                  if (bsGetKey (Link, _IntMap, &map) &&
                      map == mappg &&
                      bsGetData (Link, _bsRight, _Int, &ml1)  &&
                      bsGetData (Link, _bsRight, _Int, &ml2)) 
                    {
                      if (ml1 < ml2)
                        { c1 = c1 - ml1 + 1 ; c2 = c2 - ml1 + 1 ; }
                      if (ml1 > ml2)
                        { c1 = ml1 - c1 + 1  ; c2 = ml1 - c2 + 1 ; }
                    }
                  bsDestroy (Link) ;
                }
            }
          if (map == mappg &&
              (
               (c1 < c2 && a1 < a2  && c1 < a2 && c2 > a1) || (c1 > c2 && a1 > a2  && c1 > a2 && c2 < a1)
               ))
            { 
              { c1 = c1 ; c2 = c2 ; }
              units = arrayReCreate (units, 40, BSunit) ;
              if (bsGetArray (PG, _Source_Exons, units, 2))
                for (ii = j = 0 ; ii < arrayMax(units) ; ii+= 2)
                  {
                    uu = arrp (units, ii, BSunit) ; uu2 = uu + 2 ;
                    b1 = uu[0].i ; b2 = uu[1].i ;
                    if (ii + 3 < arrayMax(units) &&
                         uu2[0].i < b2 + 12) /* do not steal micro introns */
                      {  uu2[0].i = b1 ; continue ; }
                    up = arrayp (predictedHits, j++, HIT) ;
                    up->a1 = b1 ; up->a2 = b2 ;
                    up->type = gX ;
                    if (!ii) up->a1 -= 50 ; /* add a 5' UTR */

                    if (ii + 3 < arrayMax(units)) /* add an intron */
                      {
                        up = arrayp (predictedHits, j++, HIT) ;
                        up->a1 = b2 + 1 ; up->a2 = uu[2].i - 1 ;
                        up->type = gI ;
                      }
                    else
                      up->a2 += 200 ; /* add a 3' UTR !! */
                  }
              
              
              
              
              if (predictedHits && arrayMax(predictedHits))
                {
                  mrna = arrayp (smrnas, arrayMax(smrnas), SMRNA) ;
                  if (c1 < c2)
                    mrna->a1 = c1 - sc->a1 + 1 ;
                  else
                    mrna->a1 = sc->a1 - c1 + 1 ;
                  
                  mrna->hits = predictedHits ;
                  predictedHits = 0 ;
                }
            } 
          arrayDestroy (predictedHits) ;
          bsDestroy (PG) ;
        }
    }
  
  keySetDestroy (ks) ;
  arrayDestroy (units) ;
  return  arrayMax(smrnas) - oldMaxSmrnas ;
}

/**********************************************************************************/

static int smrnaOrderByGap (const void *a, const void *b)
{
  const SMRNA *sa = (const SMRNA*) a, *sb = (const SMRNA*) b ;

  return sa->nGaps - sb->nGaps ?
    sa->nGaps - sb->nGaps :                        /* less gaps first */
    arrayMax(sb->clones) - arrayMax(sa->clones) ;  /* more clones first */
}

/**********************************************************************************/

static int smrnaOrderByLengthandBadQuality (const void *a, const void *b)
{
  const SMRNA *sa = (const SMRNA*) a, *sb = (const SMRNA*) b ;
  int orfa, orfb ;

  if (sa->cGroup != sb->cGroup)
    return sa->cGroup - sb->cGroup ; /* baddies after goodies */
  orfa = sa->orfs && arrayMax(sa->orfs) ? arr (sa->orfs, 0, ORFS).nOpen - arr (sa->orfs, 0, ORFS).nX : 0 ;
  orfb = sb->orfs && arrayMax(sb->orfs) ? arr (sb->orfs, 0, ORFS).nOpen - arr (sb->orfs, 0, ORFS).nX : 0 ;
  
  return orfb - orfa ;  /* long first */
}

/*********************************************************************/
/* pass = 0, 1: steal from other cdna, pass = 2, steal from prediction */
static int mrnaStealFromAlternative (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas, int pass, Array linkPos)
{
  HIT *up, *vp, *wp, *up2 ;
  SMRNA *mrna, *mrna2; 
  int imrna, jmrna, a1, b1, limit = 5 ;
  int i, ii, ii2, ivp2, iwp2, withGap, foundGap = 0, allowGap = 1 ;
  Array oldHits = 0 ;
  int oldMaxSmrnas = arrayMax(smrnas) ;
  BOOL debug = FALSE ;

  arraySort (smrnas, smrnaOrderByGap) ;
  if (pass == 2 && !mrnaAddPrediction (s2m, sc, gmrna, smrnas, linkPos))
    return 0 ;
  if (pass > 0) allowGap = 0 ;

  for (imrna = 0 ; imrna < oldMaxSmrnas ; imrna++)
    {
      mrna = arrp (smrnas, imrna, SMRNA) ;
      a1 = mrna->a1 ;

      if (debug)
        {
          printf ("before Stealing, pass=%d, imrna=%d a1=%d\n", pass, imrna, a1) ;
          showHits(mrna->hits) ;
        }

      for (ii = 0 ; ii < arrayMax (mrna->hits) ; ii++)
        {
          up = arrp (mrna->hits, ii, HIT) ;
          if ((up->type & (gGap | gJ)) && /* try to resolve it using a neighbour */
              ii > 0 &&
              ((up - 1)->type & gX) &&
              ii < arrayMax(mrna->hits) - 1 &&
              ((up + 1)->type & gX))
            {
              foundGap++ ;

              vp = up - 1 ;
              wp = up + 1 ;

              /* search in other mrnas for a matching set */ 
              for (withGap = 0 ; withGap < 1 + allowGap ; withGap++) /* steal from other gaps only in second turn if allowed */
                for (jmrna = 0 ; jmrna < arrayMax (smrnas) ; jmrna++)
                  {
                    if (imrna == jmrna)
                      continue ;
                    mrna2 = arrp (smrnas, jmrna, SMRNA) ;
                    b1 = mrna2->a1 ;
                    
                    ivp2 = iwp2 = -1 ;
                    
                    limit = 5 ; /* i used 12 in feb 2005, but this is much too large
                                 * the hits have been rationalised in cleanGeneHits,
                                 * there is no reason to do it again differently here
                                 */
                    for (ii2 = 0, up2 = arrp (mrna2->hits, 0, HIT) ; ii2 < arrayMax (mrna2->hits) ; up2++, ii2++)
                      {
                        if (ivp2 >= 0 && !withGap && (up2->type & (gGap | gJ))) ivp2 = -1 ;
                        if (up2->type & gX)
                          {
                            if (
                                (up2->a1 + b1 <= vp->a1 + a1 ||  /* exact or large intersect */
                                 (up2->a1 + b1 > vp->a1 + a1  &&
                                  up2->a1 + b1 < vp->a2 + a1 - limit)
                                 ) && 
                                up2->a2 + b1  > vp->a2 + a1 - limit &&
                                ((up->type & gGap) || up2->a2 + b1  < vp->a2 + a1 + limit))
                              ivp2 = ii2 ;
                            if (ivp2 >= 0 && 
                                ( up2->a2 + b1 >= wp->a2 + a1 ||
                                  (up2->a2 + b1 < wp->a2 + a1 &&
                                   up2->a2 + b1 > wp->a1 + a1 + limit)
                                  ) &&
                                up2->a1 + b1 < wp->a1 + a1 + limit  &&
                                ((up->type & gGap) || up2->a1 + b1  > wp->a1 + a1 - limit))
                              { iwp2 = ii2 ; goto ok ; }
                          }
                      }
                  ok:
                    if (ivp2 >= 0 && iwp2 >= 0 &&
                        (
                         ((up->type & gGap) && (!withGap || iwp2 - ivp2 >= 2)) ||
                         ((up->type & gJ) && !((wp - 1)->type & (gJ & gGap)) && iwp2 - ivp2 == 2)) 
                        )
                      {
                        int olda1, olda2 ;
                        unsigned int oldSla, oldSlb ;
                        BOOL isGJ = (up->type & gJ) ? TRUE : FALSE ;

                        smrnaDnaDestroy (mrna) ;

                        oldHits = arrayCopy (mrna->hits) ;
                        /* ii -1 is the old vp, to be replaced */
                        i = ii - 1 ; 
                        up = arrayp (mrna->hits, ii - 1, HIT) ;
                        olda1 = up->a1 ;
                        oldSla = up->type & (gS | gS0 | gReal5p) ;
                        up = arrayp (mrna->hits, ii + 1, HIT) ;
                        olda2 = up->a2 ;
                        oldSlb = up->type & (gReal3p | gA) ;
                        for (ii2 = ivp2 ; ii2 <= iwp2 ; ii2++)
                          {
                            up = arrayp (mrna->hits, i++, HIT) ;
                            up2 = arrp (mrna2->hits, ii2, HIT) ;
                            *up = *up2 ;
                            if (!isGJ || ii2 == ivp2 + 1)
                              {
                                if (jmrna >= oldMaxSmrnas || (up2->type & gPredicted))
                                  up->type |= gPredicted ;
                                else
                                  up->type |= gStolen ;
                              }
                            up->a1 += b1 - a1 ;
                            up->a2 += b1 - a1 ;
                            if (ii2 == ivp2) 
                              {
                                up->a1 = olda1 ;
                                up->type &= ~(gS | gS0 | gReal5p) ;
                                up->type |= oldSla ;
                              }
                            if (ii2 == iwp2) 
                              {
                                up->a2 = olda2 ;
                                up->type &= ~(gReal3p | gA) ;
                                up->type |= oldSlb ;
                              }
                          }
                        /* analysis will proceed off ii = i (+1 induced in loop control) */
                        for (ii2 = ii + 2, ii = i - 1 ; ii2 < arrayMax(oldHits) ; ii2++)
                          {
                            up = arrayp (mrna->hits, i++, HIT) ;
                            up2 = arrp (oldHits, ii2, HIT) ;
                            *up = *up2 ;
                          }
                        arrayMax (mrna->hits) = i ;
                        arrayDestroy (oldHits) ;
                        goto nextTest ;
                      }
                  }
            }
        nextTest:
          continue ;
        }
      if (debug)
        {
          printf ("After Stealing, pass=%d, imrna=%d a1=%d\n", pass, imrna, a1) ;
          showHits(mrna->hits) ;
        }
    }
  if (pass == 2) /* destroy the predictions */
    {
      for (imrna = oldMaxSmrnas ; imrna < arrayMax(smrnas) ; imrna++)
        { 
          mrna = arrayp (smrnas, imrna, SMRNA) ;
          arrayDestroy (mrna->hits) ;
        }
      arrayMax(smrnas) = oldMaxSmrnas ;
    }

  return foundGap ;
}

/*********************************************************************/
/* disregard very small gaps between 2 exons*/
static int mrnaFillLastGap (Array smrnas, int maxGap)
{
  HIT *up, *vp, *wp ;
  SMRNA *mrna ;
  int ii, jj, imrna, foundGap = 0, foundGap2 = 0 ;

  for (imrna = 0 ; imrna < arrayMax (smrnas) ; imrna++)
    {
      mrna = arrp (smrnas, imrna, SMRNA) ;
      foundGap = 0 ;
      for (ii = jj = 0 ; ii < arrayMax (mrna->hits) ; ii++)
        {
          up = arrp (mrna->hits, ii, HIT) ;
          if (up->type & gI)
	    jj = ii ;
	}
      /* start the process down of last intron */
      for (ii = jj ; ii < arrayMax (mrna->hits) ; ii++)
        {
          up = arrp (mrna->hits, ii, HIT) ;
          if ((up->type & gGap) && /* try to resolve it using a neighbour */
              ii > 0 &&
              ((up - 1)->type & gX) &&
              ii < arrayMax(mrna->hits) - 1 &&
              ((up + 1)->type & gX) &&
	      up->a2 - up->a1 < maxGap
	      )
            {
              foundGap++ ;
              vp = up - 1 ;
              wp = up + 1 ;
	      vp->a2 = wp->a2 ;
	      vp->x2 = wp->x2 ;
	      vp->type |= wp->type ;
              (vp+1)->type = 0 ;
	      wp->type = 0 ;
	    }
	}
      if (foundGap)
	{
	  foundGap2 += foundGap ;
	  for (ii = jj = 0, up = vp = arrp (mrna->hits, 0, HIT) ; ii < arrayMax(mrna->hits) ; up++, ii++)
	    {
	      if (up->type)
		{
		  if (jj < ii) *vp = *up ;
		  jj++ ; vp++ ;
		}
	    }
	  arrayMax(mrna->hits) = jj ;
	}
    }
  return foundGap2 ;
}

/*********************************************************************/
/*********************************************************************/
/* 2008_11_28, exemple worm:srv-4
 * If i have a good-best-product witout stop
 * and if my protein matches the predicted protein
 * steal to the stop
 */ 
static int mrnaStealOne3pFromPrediction (KEY product)
{
  int ok = 0 ;
  int p1, p2, ipg, i1 ;
  char *cp ;
  KEY pg, mrna = keyGetKey (product, _mRNA) ;
  OBJ Mrna = bsUpdate (mrna) ;
  Array pep = peptideGet (product) ;
  Array pep1 = 0 ;
  KEYSET pgs = queryKey (product, ">GeneBox;>Genefinder") ;

  if (keySetMax (pgs) && pep && arrayMax (pep) > 20 &&
      bsFindKey (Mrna, _Product, product) &&
      bsGetData (Mrna, _bsRight, _Int, &p1) &&
      bsGetData (Mrna, _bsRight, _Int, &p2)
      ) ;
  else
    goto done ;

  pepDecodeArray(pep) ;
  if (0)
    {
      for (ipg = 0 ; ipg < keySetMax (pgs) ; ipg++)
	{
	  pg = keySet (pgs, ipg) ;
	  arrayDestroy (pep1) ;
	  pep1 = peptideGet (pg) ;
	  if (pep1)
	    {
	      pepDecodeArray(pep1) ;
	      i1 = arrayMax (pep) - 12 ;
	      cp = strstr (arrp (pep1, 0, char), arrp (pep, i1, char)) ;
	      if (!cp) continue ; /* the end of the proteins do not match */
	      if (0) messout (messprintf("%s could be used to elungate %s", name(pg), name(product))) ;
	    }
	}
    }
    
 done:
  keySetDestroy (pgs) ;
  arrayDestroy (pep) ;
  arrayDestroy (pep1) ;
  bsSave (Mrna) ;
  return ok ;
} /* mrnaStealOne3pFromPrediction */

/*********************************************************************/
/* find products in those genes that could be prolungated 3prime */
static int mrnaSteal3pFromPrediction (KEYSET genes)
{
  KEYSET products = query (genes, "genefinder ; >transcribed_gene ; > mrna ; ! valid3p ;> Product good_product && best_product && ! down_stop") ;
  int nn = 0, ii ;
  
  for (ii = 0 ; ii < keySetMax (products) ; ii++)
    nn += mrnaStealOne3pFromPrediction (keySet (products, ii)) ;
  
  return nn ;
} /* mrnaSteal3pFromPrediction  */

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

static void mrnaGrignotte (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas, int thisMrna, int bigGap)
{
  HIT *up, *vp ;
  KEY oldClone ;
  BOOL isUsed, foundBig = FALSE ;
  int i, j, jMax, ii, iii, j1, a1, a2, amax ;
  SMRNA *smrna ;
  BitSet bb = 0 ;

  for (iii = 0 ; iii < arrayMax(smrnas) ; iii++) 
    {
      if (thisMrna >= 0 && iii != thisMrna)
        continue ;
      smrna = arrp (smrnas, iii, SMRNA) ;
      a1 = smrna->a1 ; a2 = smrna->a2 ;

      bb = bitSetReCreate (bb, a2) ;  
      jMax =  arrayMax (smrna->hits) ;
      if (!jMax)
	continue ;
      for (j1 = 0, up = arrp (smrna->hits, 0, HIT) ; j1 < jMax ; up++, j1++)
        if (up->a2 > 0)
          bitUnSet (bb, a1 + up->a2) ;/* otherwise we sometimes go out of range */

      oldClone = 0 ; isUsed = FALSE ;
      for (j1 = 0, up = arrp (gmrna->estHits, 0, HIT) ; j1 < arrayMax (gmrna->estHits) ; up++, j1++)
        {
          if (up->cDNA_clone != oldClone)
            {
              isUsed = keySetFind (smrna->clones, up->cDNA_clone,0) ;
              oldClone = up->cDNA_clone ; 
            }
          if (isUsed && (up->type & gX))
            for (i = up->a1 > 0 ? up->a1 : 0 ; i <= up->a2 ; i++)
              bitSet (bb, i) ;
        }
      jMax =  arrayMax (smrna->hits) ;
      for (j1 = 0, up = arrp (smrna->hits, 0, HIT) ; j1 < jMax ; up++, j1++)
        {
          if (!(up->type & gX))
            continue ;
          if (up->type & gSuspect)
            continue ;

          for (j = 0, i = up->a1 ; !bit (bb, i + a1 - 1) && i <= up->a2 ; i++, j++) ;
          if (j > 0)
            {
              smrnaDnaDestroy (smrna) ;

              up->a1 += j ; up->type &= ~gS ; /* i must loose the sl if i grignotte */
              if (j1 > 0 && !((up-1)->type & gX))
                {
                  if (((up-1)->type & gI) && ((up-1)->type & gJ))
                    {
                      int k = 0 ; /* size of non supported introns */
                      for (k = 0, i = up->a1 - 1 ; 
                           i + a1 - 1 >= 0 && !bit (bb, i + a1 - 1) && i >= (up-1)->a1 - 20 ; 
                           i--, k++) ;
                      k -= ((up-1)->a2 + j + 1 - (up-1)->a1) ;
                      if (k < 0) k = -k ;
                      if (k > 5)
                         (up-1)->type = gGap ;
                       (up-1)->a2 += j ; 
                    }
                  else
                    {
                      (up-1)->a2 += j ; 
                      (up-1)->type = gGap ;
                    }
                }
            }

          for (j = 0, i = up->a2 ; !bit (bb, i + a1 - 1) && i >= up->a1 ; i--, j++) ;
          if (j > 0)
            {
              up->a2 -= j ;
              smrnaDnaDestroy (smrna) ;
        
              if (j1 < jMax -1 && !((up+1)->type & gX))
                {
                  if (((up+1)->type & gI) && ((up+1)->type & gJ))
                    {
                      int k = 0 ; /* size of non supported introns */
                      for (k = 0, i = up->a2 + 1 ; 
                           !bit (bb, i + a1 - 1) && i < a2 - a1 ;
                           i++, k++) ;
                      k -= ((up+1)->a2 - j + 1 - (up+1)->a1) ;
                      if (k < 0) k = -k ;
                      if (k > 5)
                        (up+1)->type = gGap ;
                       (up+1)->a1 -= j ; 
                    }
                  else
                    {
                      (up+1)->a1 -= j ; 
                      (up+1)->type = gGap ; 
                    }
                }
            }
          if (up->a1 >= up->a2)
            up->type = 0 ;
	  if (bigGap)
	    {
	      jMax =  arrayMax (smrna->hits) ;
	      for (j1 = 0, up = arrp (smrna->hits, 0, HIT) ; j1 < jMax ; up++, j1++)
		if ((up->type & gGap) && up->a2 - up->a1 > bigGap)
		  {
		    int k ;
		    
		    foundBig = TRUE ;
		    if (2*(j1+1) < jMax) /* cut above */
		      {
			for (vp = up, k = j1 ; k >= 0 ; vp--, k--)
			  vp->type = 0 ;
		      }
		    else /* cut below */
		      {
			for (vp = up, k = j1 ; k < jMax ; vp++, k++)
			  vp->type = 0 ;
		      }
		  }
	    }
        }
      /* keep happy few */
      for (ii = 0, j1 = 0, up = arrp (smrna->hits, 0, HIT), vp = up ; ii < arrayMax (smrna->hits) ; up++, ii++)
        {
          if (!j1 && !(up->type & gX))
            continue ;
          if (up->type != 0)
            {
              if (ii != j1)
                *vp = *up ;
              if (j1 > 0 && !(vp->type & gX) && !((vp-1)->type & gX))
                {
                  (vp - 1)->type = gGap ; 
                  (vp - 1)->a2 = vp->a2 ;
                }
              else
                {
                  vp++ ; j1++ ;
                }
            }
        }

      while (j1 > 0 && (up = arrp (smrna->hits, j1 - 1, HIT)) && ! (up->type & gX))
        j1-- ;
      arrayMax (smrna->hits) = j1 ;  
      if (!j1)
        continue ;
      /* shift origin of coordinates */
      arraySort (smrna->hits, cDNAOrderGloballyByA1) ;
      j = arrp(smrna->hits, 0, HIT)->a1 - 1 ;
      amax = 0 ;
      smrna->a1 += j ;
      for (ii = 0, j1 = 0, up = arrp (smrna->hits, 0, HIT); 
           ii < arrayMax (smrna->hits) ; up++, ii++)
        { 
          up->a1 -= j ; 
          up->a2 -= j ;
          if (up->a2 > amax) amax = up->a2 ;
        } 
      smrna->a2 = smrna->a1 + amax - 1 ;          
    }
  bitSetDestroy (bb) ;

  /* reset the mrna origins */
  if (thisMrna >= 0)
    return ;
  a1 = 999999 ; 
  for (ii = iii = j1 = 0 ; ii < arrayMax(smrnas) ; ii++) 
    { 
      smrna = arrp (smrnas, ii, SMRNA) ; 
      if (arrayMax(smrna->hits))
        for (j = 0, j1 = 0, up = arrp (smrna->hits, 0, HIT); 
             !j1 && j < arrayMax (smrna->hits) ; up++, j++)
          if (up->type & gX)
            j1 = 1 ;

      if (j1 && arrayMax(smrna->hits))
        {
	  a2 = -999999 ;
	  for (j = 0,  up = arrp (smrna->hits, 0, HIT); 
	       j < arrayMax (smrna->hits) ; up++, j++)
	    if ((up->type & gX) && a2 < up->a2)
	      a2 = up->a2 ;
	  if (a2 > 0)
	    smrna->a2 = smrna->a1 + a2 - 1 ;
	  j = smrna->a1 ;
          if (j < a1)
            a1 = j ;
          if (iii != ii)
            arr (smrnas, iii, SMRNA) = arr (smrnas, ii, SMRNA) ;
          iii++ ;
        }
    } 
  arrayMax(smrnas) = iii ;
  /* reset the gene origin */
  j = a1 - 1 ;
  if (foundBig || j > 0)
    {
      for (iii =  amax = 0 ; iii < arrayMax(smrnas) ; iii++) 
        {
          smrna = arrp (smrnas, iii, SMRNA) ; 
          smrna->a1 -= j ; smrna->a2 -= j ;
          if (smrna->a2 > amax)
            amax = smrna->a2 ;
          /*
          for (ii = 0, up = arrp (smrna->hits, 0, HIT); 
               ii < arrayMax (smrna->hits) ; up++, ii++)
            { 
              up->a1 -= j ; 
              up->a2 -= j ;
            } 
          */
        } 
      for (ii = 0, up = arrp (gmrna->estHits, 0, HIT); 
               ii < arrayMax (gmrna->estHits) ; up++, ii++)
        { 
          up->a1 -= j ; 
          up->a2 -= j ;
        } 
      for (ii = 0, up = arrp (gmrna->hits, 0, HIT); 
               ii < arrayMax (gmrna->hits) ; up++, ii++)
        { 
          up->a1 -= j ; 
          up->a2 -= j ;
        } 
      up = arrp (gmrna->hits, 0, HIT) ; up->a1 = 1 ;
      up = arrp (gmrna->hits, ii - 1, HIT) ; up->a2 = amax ;
      if (sc->a1 < sc->a2)
        { sc->a1 += j ; sc->a2 = sc->a1 + amax - 1 ; }
      else
        { sc->a1 -= j ; sc->a2 = sc->a1 - amax + 1 ; }        
    }
}  /* mrnaGrignotte */

/*********************************************************************/
/* flags about SL1 CDS polyA which were set before should be revised
 *  because the list of clones belonging to the mRNA has changed
 */

static void mrnaSetOneCompletenessFlag (HIT *up, SMRNA *smrna, Array estHits, BOOL isTop)
{
  HIT *vp ;
  KEY oldClone ;
  BOOL isUsed ;
  int j1 ;
  unsigned int topFlags = gS | gS0 | gReal5p | gCompleteCDS ;
  unsigned int bottomFlags = gReal3p | gA ;
  unsigned int cFlags = isTop ? topFlags : bottomFlags ;

  
  if (!(up->type & gX) || !(up->type & cFlags))
    return ;
  up->type &= ~cFlags ;
  
  cFlags &= ~gCompleteCDS ;  /* do not consider this flag any more, danielle, jan 23, 2004 */
  oldClone = 0 ; isUsed = FALSE ;
  for (j1 = 0, vp = arrp (estHits, 0, HIT) ; j1 < arrayMax (estHits) ; vp++, j1++)
    {
      if (vp->cDNA_clone != oldClone)
        {
          isUsed = keySetFind (smrna->clones, vp->cDNA_clone,0) ;
          oldClone = vp->cDNA_clone ; 
        }
      if (!isUsed || ! (vp->type &  cFlags))
        continue ;
      if (isTop && up->x1 < up->x2 && up->x1 > up->clipTop + 10)
        continue ;
      if (isTop && up->x1 > up->x2 && up->x2 < up->clipEnd - 10)
        continue ;
      if (smrna->a1 + up->a1 - 1 < vp->a2 && smrna->a1 + up->a2 - 1 > vp->a1)
        up->type |= vp->type & cFlags ;
    }
}  /* mrnaSetOneCompletenessFlag */

/******************/

static void mrnaSetCompletenessFlags (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas)
{
  HIT *up ;
  int iii, j ;
  SMRNA *smrna ;

  for (iii = 0 ; iii < arrayMax(smrnas) ; iii++) 
    {
      smrna = arrp (smrnas, iii, SMRNA) ;
      if (smrna->composite)
	continue ;
      if (! smrna->hits)
	continue ;
      /* search the flags of the first exon */
      for (j = 0, up = arrp (smrna->hits, 0, HIT) ; j < 1 && j < arrayMax(smrna->hits);  up++, j++)
        mrnaSetOneCompletenessFlag (up, smrna, gmrna->estHits, TRUE) ;

       /* search the polyA flags of the last exon */
      for (j = arrayMax(smrna->hits) - 1 ; j >= 0 && (up = arrp (smrna->hits, j, HIT)) ; j = -1)
        mrnaSetOneCompletenessFlag (up, smrna, gmrna->estHits, FALSE) ;
    }
} /* mrnaSetCompletenessFlags */

/*********************************************************************/
/*********************************************************************/

typedef struct { KEY product, mrna, pepKey ; int ln ; Array pep ; BOOL complete ; } PHIT ;

static int pHitLengthOrder (const void *a, const void *b)
{
  const PHIT *pa = (const PHIT*)a, *pb = (const PHIT *)b ;
  
  if (pa->ln != pb->ln)
    return pb->ln - pa->ln ;
  return pa->product - pb->product ;
}

/*********************************************************************/
/*********************************************************************/
#ifdef JUNK
static int mrnaCountIntronsInCds (KEY product)
{
  Array units = arrayCreate (80, BSunit) ;
  BSunit *uu ;
  int ii, nn = -1 ;
  OBJ Product = bsCreate (product) ;

  if (Product && bsGetArray (Product, _Source_Exons, units, 3))
    {
      for (ii = 0 ; ii < arrayMax(units) ; ii += 3)
        {
          uu = arrp (units, 0, BSunit) ;
          if (strstr(name(uu[2].k), "xon")) 
            nn++ ;
        }
    }

  arrayDestroy (units) ;
  bsDestroy (Product) ;
  return nn > 0 ? nn : 0 ;
}
#endif
/*********************************************************************/
/*********************************************************************/
#ifdef JUNK
static BOOL mrnaCountCdsTilingClones (OBJ Mrna, int *np)
{
  int n = 0 ;
  KEY dummy = 0 ;
  
  if (bsGetKey(Mrna, str2tag("CDS_covered_by"), &dummy))
    do { n++ ; } while (bsGetKey (Mrna, _bsDown, &dummy)) ;
  if (np) *np = n ;

  return n > 0 ;
}

static void mrnaFlagIntrons (KEYSET mrnaSet, Array introns)
{
  int a1, a2, ii, jj ;
  HIT *up, *vp ;
  
  for (ii = 0 ; ii < arrayMax (introns) ; ii++)
    {
      up = arrp (introns, ii, HIT) ;
      if (keySetFind (mrnaSet, up->gene, 0))
        up->x1 = 1 ;
    }
  /* now we flag the other introns representatives */
  for (ii = 0 ; ii < arrayMax (introns) ; ii++)
    {
      up = arrp (introns, ii, HIT) ;
      if (!up->x1)
        { /* search up */
          a1 = up->a1 ; a2 = up->a2 ;
          vp = up ; jj = ii ;
          while (vp--, jj--)
            {
              if (vp->a1 != a1 || vp->a2 != a2)
                break ;
              else if (vp->x1)
                { up->x1 = 1 ; break ; }
            }
        }
      if (!up->x1) /* search down */
        {
          vp = up ; jj = ii ;
          while (vp++, ++jj < arrayMax (introns))
            {
              if (vp->a1 != a1 || vp->a2 != a2)
                break ;
              else if (vp->x1)
                { up->x1 = 1 ; break ; }
            }
        }
    }
}
#endif
/*********************************************************************/
/* on post bricole les scores en fontions de la hierarchie des products */
#ifdef JUNK
static int mrnaHierarchyOrder (const void *va, const void *vb)
{
   const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;
   
   if (vp->a1 - up->a1)
     return vp->a1 - up->a1 ; /* by decreasing score */
   return 
     up->a2 - vp->a2 ; /* increasing number of covering clones */
}
#endif

/*********************************************************************/

static void mrnaSaveMrnaHierarchy (KEY tg)
{
  KEY mrna, _Best_in_gene ;
  OBJ Mrna = 0 ;
  KEYSET mrnaSet = 0 ;
  int ii ;
  
  _Best_in_gene = str2tag ("Best_in_gene") ;

  /* clean up */
  mrnaSet = queryKey (tg, "> mrna ; best_in_gene") ;
  for (ii = 0 ; ii < keySetMax (mrnaSet) ; ii++)
    {
      mrna = keySet (mrnaSet, ii) ;
      if ((Mrna = bsUpdate (mrna)))
        {
          bsRemoveTag (Mrna, _Best_in_gene) ;
          bsSave (Mrna) ;
        }
    }
  keySetDestroy (mrnaSet) ;

  /* mark good ones */
  /* 2024_01_22  in f4 reconstruction the mrna do NOT ontain cdna_clones, so i changed to COUNT () >=0 rather than >0 */ 
  mrnaSet = queryKey (tg, "{COUNT mrna = 1 ; > mrna ;} SETOR {>mrna ; COUNT {>cdna_clone ; ! double_fuzzy } >= 0 ; gt_ag || gc_ac || (!gap && !fuzzy && !other)}") ;
  for (ii = 0 ; ii < keySetMax (mrnaSet) ; ii++)
    {
      mrna = keySet (mrnaSet, ii) ;
      if ((Mrna = bsUpdate (mrna)))
        {
          bsAddTag (Mrna, _Best_in_gene) ;
          bsSave (Mrna) ;
        }
    }
  keySetDestroy (mrnaSet) ;

  /* mark well supported mrnas */
  mrnaSet = queryKey (tg, "> mrna ; COUNT {>product best_product && complete_cds_clone} > 0") ;
  for (ii = 0 ; ii < keySetMax (mrnaSet) ; ii++)
    {
      mrna = keySet (mrnaSet, ii) ;
      if ((Mrna = bsUpdate (mrna)))
        {
          bsAddTag (Mrna, str2tag ("Well_supported")) ;
          bsSave (Mrna) ;
        }
    }
  keySetDestroy (mrnaSet) ;

  /* save most compact form of genes i would otherwise destroy */
  mrnaSet = queryKey (tg, "COUNT { >mrna ; best_in_gene } == 0 ") ;
  if (!keySetMax (mrnaSet))
    {
      BSunit *uu ;
      Array units = arrayCreate (36, BSunit) ;
      int len , j0, ii, m1, m2 ;
      OBJ Tg = 0 ;

      if ((Tg = bsCreate (tg)))
        {
          if (bsGetArray (Tg, str2tag ("mMRNA"), units, 3))
            {
              for (ii = len = 0, j0 = -1 ; ii < arrayMax (units) ; ii+= 3)
                {
                  uu = arrp (units, ii, BSunit) ;
                  m1 = uu[1].i ;
                  m2 = uu[2].i ;
                  if (len == 0 || len > m2 - m1)
                    { len = m2 - m1 ; j0 = ii ; }
                }
              if (j0 >= 0)
                {
                  uu = arrp (units, ii, BSunit) ;
                  mrna = uu[0].k ;
                  if ((Mrna = bsUpdate (mrna)))
                    {
                      bsAddTag (Mrna, _Best_in_gene) ;
                      bsSave (Mrna) ;
                    }
                }
            }
          arrayDestroy (units) ;
          bsDestroy (Tg) ;
        }
    }
 
  keySetDestroy (mrnaSet) ;
} /* mrnaSaveMrnaHierarchy */

/*********************************************************************/
/*********************************************************************/
/* x1, x2 are est coords, p1 p2 are product coords in the mRNA dna */
static void mrnaDoSaveProductCoverage (KEY product, KEY est, int p1, int p2
                                       , int x1, int x2, Array mrnaPep)
{
  OBJ Est = bsUpdate (est), Product ;
  Array dna = 0, dna1 = 0 ;
  int i, y1 = 0, yy1, y2 = 0, nX = 0, nOk = 0, nV, pepL, delta, jj, iDelta ;
  int yDelta [51] ; /*  = { 0, 1, -1, 2, -2 , 3, -3, 4, -4, .... ,  0 } ; */
  char *cp, cc, *cq, buf[256] ;

  for (i = 1 ; i < 26 ; i++)
    {
      yDelta [2*i - 1] = i ;
      yDelta [2*i] = -i ;
    }
  yDelta [49] = 0 ;
  yDelta [0] = 0 ;
    
  memset (buf, 0, sizeof(buf)) ;
  if ((dna = getSeqDna (est)))
    {
      if (x1 > x2) 
        { 
          dna1 = dnaCopy (dna) ;
          reverseComplement (dna1) ;
          y1 = arrayMax (dna) - x1 + 1 ;
          y2 = arrayMax (dna) - x2 + 1 ;
        }
      else
        { dna1 = dna ; y1 = x1 ; y2 = x2 ; }
      
      for (iDelta = 0 ; iDelta < 50 ; iDelta++)
        {
          yy1 = y1 + yDelta [iDelta] ;
          nOk = nX = pepL = nV = jj = 0 ;
          buf[0] = 0 ;
          if (yy1 > 0 && yy1 <= arrayMax(dna1) &&
              y2 > 0 && y2 <= arrayMax(dna1))
            /* we jump the first letter which may be a forced leu->met */
            for (pepL = nOk = jj = 1, i = yy1 + 3, cp = arrp (dna1, i - 1, char)
                   , cq = arrp (mrnaPep, 1, char) ;
                 i < y2 ; i += 3, cp += 3, cq++, jj++)
              {
                cc = codon (cp) ;
                if (!*cq)
                  break ;
                if (cc == *cq)
                  nOk++ ;
                else
                  {
                    nV++ ;
                    if (nV <= 5)
                      strcat (buf, messprintf ("%c%d%c ", *cq, jj + 1, cc)) ;
                    else if (nV == 6)
                      strcat (buf, "...") ;
                  }
                if (cc == '*')
                  break ;
                if (nOk > 5 && nV < 2) /* i am in frame */
                  iDelta = 1000 ; /* will loop out */
                if (cc == 'X')
                  nX++ ;
                pepL++ ;
              }
        }
    }

  delta = pepL -  arrayMax (mrnaPep) ;
  if (pepL)
    { 
      bsAddKey (Est, str2tag ("Covers_product"), product) ;
      bsPushObj (Est) ;
      bsAddData (Est, str2tag ("CDS_length"), _Int, &pepL) ; 
      bsAddData (Est, _bsRight, _Text, "AA") ;
      bsAddData (Est, str2tag ("From"), _Int, &x1) ; 
      bsAddData (Est, str2tag ("To"), _Int, &x2) ; 
      
      if (nX)
        bsAddData (Est, str2tag ("nX"), _Int, &nX) ;
      bsAddData (Est, str2tag ("delta"), _Int, &delta) ;
      if (delta)
        bsAddData (Est, _bsRight, _Text, "aa") ;
      
      bsAddData (Est, str2tag ("nV"), _Int, &nV) ;
      if (nV)
        {
          bsAddData (Est, _bsRight, _Text, "aa") ;
          bsAddData (Est, str2tag ("Type"), _Text, buf) ;
        }
    }
  /* add coverage details */
  bsSave (Est) ;
  
  if (pepL && (Product = bsUpdate (product)))
    {
      bsAddKey (Product, str2tag ("Covered_by"), est) ;
      bsAddData (Product, _bsRight, _Text, "Delta") ;
      bsAddData (Product, _bsRight, _Int, &delta) ;
      bsAddData (Product, _bsRight, _Text, "nV") ;
      bsAddData (Product, _bsRight, _Int, &nV) ;
      if (buf[0])
        bsAddData (Product, _bsRight, _Text, buf) ;
      bsSave (Product) ;
    }
  
  if (dna1 != dna)
    arrayDestroy (dna1) ;
  dna = 0 ; /* do not destroy */
}  /* mrnaDoSaveProductCoverage */

/*************************************/
     
static void mrnaSaveProductCoverage (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas, SMRNA *smrna)
{
  KEY mrna = smrna->tr, product ;
  KEYSET ks = keySetCreate () ;
  Array units = arrayCreate (32, BSunit) ;
  Array vnits = arrayCreate (32, BSunit) ;
  Array estHits = gmrna->estHits ;
  OBJ Mrna = bsCreate (mrna) ;
  BSunit *uu, *vv ;
  HIT *up, *vp ;
  int ii, iip, i2, jj, p1, p2, q1, q2, a1, a2, x1, x2, y1, y2 ;
  int isComplete ;
  Array peptide = 0 ;

  if (Mrna && bsGetArray (Mrna, _Product, units, 3) &&
      bsGetArray (Mrna, _Coding, vnits, 4) )
    for (iip = 0 ; iip < arrayMax (units) ; iip += 3)
      {
        uu = arrp (units, iip, BSunit) ;
        product = uu[0].k ;
        p1 = uu[1].i ;
        p2 = uu[2].i ; /* spliced mrna coordinates */
        q1 = q2 = -1 ; /* unspliced mrna coordinates */

        for (i2 = 0 ; i2 < arrayMax (vnits) ; i2 += 4)
          {
            vv = arrp (vnits, i2, BSunit) ;
            a1 = vv[0].i ;            
            a2 = vv[1].i ;
            x1 = vv[2].i ;            
            x2 = vv[3].i ;
            if (x1 <= p1 && x1 >= p1) 
              q1 = a1 + p1 - x1 ;
            if (x1 <= p2 && x2 >= p2) 
              q2 = a1 + p2 - x1 ;
          }
        arrayDestroy (peptide) ;
        peptide = peptideGet (product) ;
        if (!peptide || arrayMax(peptide) < 3 || q1 < 0 || q2 < 0)
          continue ;
        pepDecodeArray (peptide) ;
        q1 += smrna->a1 - 1 ;
        q2 += smrna->a1 - 1 ;  /* unspliced gene coordinates */

        for (ii = 0 ; ii < arrayMax (estHits) ; ii++)
          {
            up = arrp (estHits, ii, HIT) ;
            if (!class(up->est)) /* ghosts */
              continue ;
            if (keyFindTag (up->est, str2tag("duplicate_of_read")))
              continue ;
            if (! (up->type & gX))
              continue ;

            if (!keySetFind (smrna->clones, up->cDNA_clone, 0))
              continue ;
            
            if (keySetFind (ks, up->est, 0))
              continue ;
            
            isComplete = 0 ; y1 = y2 = 0 ;
            for (jj = ii ; jj < arrayMax (estHits) ; jj++)
              {
                vp = arrp (estHits, jj, HIT) ;        
                if (vp->est != up->est)
                  continue ;
                x1 = vp->x1 ; x2 = vp->x2 ; 
                a1 = vp->a1 ; a2 = vp->a2 ;
                if (a1 <= q1 && a2 >= q1)
                  { isComplete |= 1 ; y1 = x1 < x2 ? x1 + q1 - a1 : x1 - q1 + a1 ; }
                if (a1 <= q2 && a2 >= q2)
                  { isComplete |= 2  ; y2 = x1 < x2 ? x1 + q2 - a1 : x1 - q2 + a1 ; }
                if (a2 >= q2)
                  break ;
              }
            if (isComplete == 3)
              mrnaDoSaveProductCoverage (product, up->est, p1, p2, y1, y2, peptide) ;
          }        
      }
  bsDestroy (Mrna) ;
  if ((Mrna = bsUpdate (mrna)))
    {
      if (bsGetData (Mrna, str2tag("cds_covered_by "), _Int, &y1) && y1==1)
	bsAddTag (Mrna, str2tag ("Well_supported")) ;
      bsSave (Mrna) ;
    }
  arrayDestroy (units) ;
  arrayDestroy (vnits) ;
  keySetDestroy (ks) ;
  arrayDestroy (peptide) ;
  return ;
} /* mrnaSaveProductCoverage */

/*********************************************************************/
/*********************************************************************/

static void mrnaSaveMicroIntrons (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas)
{
  HIT *up ;
  int ii, jj, nn, u1, u2, a1 ;
  KEY foot ;
  OBJ Est = 0, Gene = 0, Mrna = 0 ;
  KEY _Small_deletion = str2tag ("Small_deletion") ;
  SMRNA *smrna ;

  for (ii = 0, up = arrp (gmrna->estHits, 0, HIT) ; ii < arrayMax(gmrna->estHits) ; ii++, up++)
    if ((up->type & gMicro) && (Est = bsUpdate (up->est))) 
      {
        if (sc->isUp)
          { u1 = sc->d1 - up->a1 + 1 ; u2 = sc->d1 - up->a2 + 1 ; }
        else
          { u1 = sc->d1 + up->a1 - 1 ; u2 = sc->d1 + up->a2 - 1 ; }
        if (0) printf ("Micro found in %s\n", name(up->est)) ;
        foot = mrnaIntronBoundary (s2m->dnaD, s2m->dnaR, u1, u2) ;
        bsAddData (Est, _Small_deletion, _Text, name(foot)) ;
        bsAddData (Est, _bsRight, _Text, "Position") ;
        bsAddData (Est, _bsRight, _Int, &(up->x1)) ;
        bsAddData (Est, _bsRight, _Text, "Length") ;
        nn = up->a2 - up->a1 + 1 ;
        bsAddData (Est, _bsRight, _Int, &nn) ;
        bsSave (Est) ;

        if ((Gene = bsUpdate (sc->gene)))
          {
            bsAddData (Gene, _Small_deletion, _Text, name(foot)) ;
            bsAddData (Gene, _bsRight, _Int, &(up->a1)) ;
            bsAddData (Gene, _bsRight, _Int, &(up->a2)) ;
            bsAddKey ( Gene, _bsRight, up->cDNA_clone) ;
         
            for (jj = 0 ; jj < arrayMax (smrnas) ; jj++)
              { 
                smrna = arrp (smrnas, jj, SMRNA) ;
                if (bsFindKey (Gene, _mRNA, smrna->tr) &&
                    bsGetData (Gene, _bsRight, _Int, &a1) &&
                    (Mrna = bsUpdate (smrna->tr)))
                  {
                    /* in principle smrna->a1 is equal to a1 extracted from gene
                     * but we may have stole one base, example 1B33
                     */
                    if (bsFindKey (Mrna, _cDNA_clone, up->cDNA_clone))
                      {
                        bsAddData (Mrna, _Small_deletion, _Text, name(foot)) ;
                        nn = up->a1 - a1 + 1 ;
                        bsAddData (Mrna, _bsRight, _Int, &nn) ;
                        nn = up->a2 - up->a1 + 1 ;
                        bsAddData (Mrna, _bsRight, _Int, &nn) ;
                        bsAddKey (Mrna, _bsRight, up->cDNA_clone) ;
                      }
                    bsSave (Mrna) ;
                  }
              }
            bsSave (Gene) ;
          }
      }
} /* mrnaSaveMicroIntrons */

/*********************************************************************/
/*********************************************************************/

static BOOL mrnaNMD (KEY mrna, KEY bestProduct)
{
  OBJ Mrna = bsUpdate (mrna) ;
  BOOL ok = FALSE ;
  int ii, p1, p2 ;
  Array units = 0 ;
  BSunit *uu ;
  KEY _NMD = str2tag("NMD") ;

  units = arrayReCreate (units, 8, BSunit) ;
  if (Mrna)
    {
      bsRemoveTag (Mrna, _NMD) ;
      if (
          !bsFindTag (Mrna, _Gap) &&
          !bsFindTag (Mrna, _Nb_predicted_exons) &&
          !bsFindTag (Mrna, _Nb_stolen_exons) &&
          bsFindKey (Mrna, _Product, bestProduct) &&
          bsGetData (Mrna, _bsRight, _Int, &p1) &&
          bsGetData (Mrna, _bsRight, _Int, &p2) &&
          keyFindTag ( bestProduct, _Good_product) &&
          bsGetArray (Mrna, _Splicing, units, 6))
        {
          for (ii = 0 ; ii < arrayMax(units) ; ii += 6)
            {
              uu = arrp (units, ii, BSunit) ;
              if (uu[2].i > p2+50 && (uu[1].i - uu[0].i) > 40 && strstr(name(uu[4].k),"tron"))
                bsAddTag (Mrna, _NMD) ;
            }
        }
      bsSave (Mrna) ;
    }
  arrayDestroy (units) ;
  return ok ;
} /* mrnaNMD  */

/*********************************************************************/

static void mrnaSelectBestProduct (KEY mrna)
{
  KEY product, bestProduct = 0 ;
  KEYSET pFams = 0 ;
  OBJ Mrna = 0 ;
  OBJ Product = 0 ;
  BSunit *uu ;
  HIT *up, *vp ;
  int ii, jj, p1, p2, p1Min, max = 0, db, dc, c1, c2, nvg, xx, met, leu, nIn, nOut ;
  char *cp ;
  Array units = 0 ;
  Array hits = 0 ;

  if (!keyFindTag (mrna, _Product))
    return ;

  max = 0 ; nvg = 0 ;
  units = arrayReCreate (units, 128, BSunit) ;
  hits = arrayReCreate (hits, 12, HIT) ;
  Mrna = bsUpdate (mrna) ;
  bsGetArray (Mrna, _Product, units, 3) ;
  /* collect products coordinates */
  p1Min = 9999999 ;
  for (ii = jj = 0 ; ii < arrayMax(units) ; ii += 3)
    {
      uu = arrp (units, ii, BSunit) ;
      if (p1Min > uu[1].i)
        p1Min = uu[1].i ;
    }
      
   for (ii = jj = 0 ; ii < arrayMax(units) ; ii += 3)
    {
      uu = arrp (units, ii, BSunit) ;
      product = uu[0].k ;
      Product = bsUpdate (product) ;
      
      p1 = uu[1].i ; p2 = uu[2].i ;
      up = arrayp (hits, jj++, HIT) ;
      up->x1 = p1 ; up->x2 = p2 ; up->cDNA_clone = product ;
      
      /* count the points */
      up->a1 = 0 ;
      up->a2 = (p2 - p1)/300 ; up->nerr = 0 ;
      if (p2 - p1 < 180)
        up->a2-- ; /* kill very short peptides unless well supported */
      met = leu = -99 ;
      bsGetData (Product, str2tag ("First_ATG"), _Int, &met) ;
      bsGetData (Product, str2tag ("First_Kozak"), _Int, &leu) ;
      if ( bsFindTag (Product, _Complete) &&
           (
            (met == -99 && leu == 1 && p2 - p1 < MINI_LEU2END) ||
            (met > 1 && leu == 1 && p2 - p1 - 3*met < MINI_LEU2END)
            )
          )
        { up->a2 -= 1 ; up->nerr = 150 ; }
      nIn = nOut = 0 ;
      bsGetData (Product, str2tag ("Nb_introns_in_CDS"), _Int, &nIn) ; /* nb classic introns */
      bsGetData (Product, str2tag ("Nb_Introns_outside_CDS"), _Int, &nOut) ;
      if (nIn == 1 && nOut == 1) ;
      else if (nIn == 0 && nOut == 1) up->a2 -= 1 ; /* was forgotten: added 2007_12_17 */
      else if (nOut > 1) up->a2 += nIn - (nOut < 4 ? nOut : 4) ; /* do not believe double NMD */
      else up->a2 += nIn ;

      /* scale up */
      up->a2 = 1000 * up->a2 ;

      {
	/* upper bound of 1.5 points for any kind of conservation */
	int bonusConserved = 0 ;
	int spc=ficheDBCreatureID (0) ; 
	char *speciesName = SpcI[spc].speciesName ;
	
	if (bsGetData (Product, str2tag ("Tax_common_ancestor"), _Text, &cp) &&
	    speciesName && ! strstr (cp, speciesName))
	   bonusConserved++ ;
	if (0 &&  /* superseaded by Tax_common_ancestor */
	    bsFindTag (Product, _Blastp))
	   bonusConserved++ ;
	if (0 && bsFindTag (Product, _Pfam)) /* do not give points for dubious PFam */
	   bonusConserved++ ;
	pFams = queryKey (product, ">PFam ; ! IS rvt_1 && ! IS transposase_* && ! IS ribosomal_* && ! IS  gag_* && ! IS  rve && ! IS  rvp && ! IS  dde && ! IS  gp36") ;
	if (keySetMax (pFams))
	   bonusConserved++ ;
	keySetDestroy (pFams) ;
	if (bsFindTag (Product, str2tag ("AceKogWorm")))
	  ;
	if (bsFindTag (Product, str2tag ("Transmembrane_domain")) ||
	    bsFindTag (Product, str2tag ("Coiled_coil_region")) ||
	    bsFindTag (Product, str2tag ("ER_retention_domain")) ||
	    bsFindTag (Product, str2tag ("Golgi_transport_domain")) ||
          bsFindTag (Product, str2tag ("N_myristoylation_domain")) ||
	    bsFindTag (Product, str2tag ("prenylation_domain")) ||
	    (
	     bsGetData (Product, str2tag ("psort_title"), _Text, &cp) &&
	     !strstr (cp, "cytoplasmic") &&  !strstr (cp, "nuclear") &&
	     !strstr (cp, "membrane") 
	     )
	    )
	  bonusConserved++ ;

	if (bonusConserved == 1) up->a2 += 1000 ;
	else if (bonusConserved > 1) up->a2 += 1500 ;
      }
      up->a2 +=  ((p2 - p1) % 300) + up->nerr ; /* departage les ex aequo */
      if (up->a2 > 1000 && up->x1 == p1Min && bsFindTag (Product, _NH2_complete)) up->a2 += 800 ; 
      if (up->a2 > max) max = up->a2 ;
      bsSave (Product) ;
    }
  
  for (ii = 0, up = arrp (hits, 0, HIT) ; ii < arrayMax(hits) ; up++, ii++)
    up->a2 = max + 1 - up->a2 ;
  arraySort (hits, cDNAOrderGloballyByA1) ;
  for (ii = 0, up = arrp (hits, 0, HIT) ; ii < arrayMax(hits) ; up++, ii++)
    up->a2 = max + 1 - up->a2 ;
  
  for (ii = 0, up = arrp (hits, 0, HIT) ; ii < arrayMax(hits) ; up++, ii++)
    {
      /* remove included products */
      for (jj = ii + 1, vp = up + 1 ; jj < arrayMax(hits) ; vp++, jj++)
        {
          if (vp->a2)
            {
              c1 = up->x1 > vp->x1 ? up->x1 : vp->x1 ;
              c2 = up->x2 < vp->x2 ? up->x2 : vp->x2 ;
              dc = c2 - c1 ; db = vp->x2 - vp->x1 ;
              if ((10 * dc > 5 * db && vp->a2 < 5000) ||
                  10 * dc > 8 * db)
                vp->a2 = 0 ;
            }
        }
      product = up->cDNA_clone ;
      Product = bsUpdate (product) ;
      if (ii == 0)
        {
          bsAddTag (Product, _Best_product) ;
          if (up->a2 > 1000)
            bestProduct = product ;
        }
      else
	bsRemoveTag (Product, _Best_product) ;
      if (up->a2 > 1000)
        bsAddTag (Product, _Good_product) ;
      else
	bsRemoveTag (Product, _Good_product);
      if (up->a2 > 5000)
        {
          bsAddTag (Product, _Very_good_product) ;
          if (!ii) nvg++ ;
        }
      else
	bsRemoveTag (Product, _Very_good_product) ;
      xx = up->a2/1000 ;
      bsAddData (Product, str2tag("Quality"), _Int, &xx) ;
      bsSave (Product) ; 
    } 
  if (nvg)
    bsAddTag (Mrna, str2tag("Ambiguous_best_product")) ;
  else
    bsRemoveTag (Mrna, str2tag("Ambiguous_best_product")) ;
  bsSave (Mrna) ;

/* check for NMD candidates */
  if (bestProduct)
    mrnaNMD (mrna, bestProduct) ;

  arrayDestroy (units) ;
  arrayDestroy (hits) ;
}  /* mrnaSelectBestProduct */

/*********************************************************************/

static void  mrnaSaveProductHierarchy (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas)
{
  KEY gene = sc->gene, product, geneBox = 0 ;
  OBJ Product = 0 ;
  KEYSET productSet = 0 ;
  Array pHits = 0 ;
  PHIT *ph, *qh ;
  int ii, jj, maxP ;
  char *cp, *cq ;

  /* kill products of non coding genes */
  if ((geneBox = keyGetKey (gene, _Gene)) && 
      keyFindTag (geneBox, str2tag("Non_protein_coding"))
      )
    {
      productSet = queryKey (gene, ">Product ; coding_length > 180") ;
      maxP = keySetMax (productSet) ; 
      keySetDestroy (productSet) ;
      
      if (maxP)
        messout ("Gene %s has tag Non_protein_coding, yet it encode a protein of more than 80 AA, please check",
                   name (geneBox)) ;
      if (1)
        {
          OBJ Mrna = 0 ;
          KEYSET mrnas ;
          
          productSet = queryKey (gene, ">mrna;{>product} SETOR {>Product;>Peptide;}") ;
          mrnas = queryKey (gene, ">mrna") ;
          for (ii = 0 ; ii < keySetMax (mrnas) ; ii++)
            {
              if ((Mrna = bsUpdate (keySet (mrnas, ii))))
                {
                  bsRemoveTag (Mrna, str2tag ("Longest_ORF")) ;
                  bsRemoveTag (Mrna, str2tag ("Longest_CDS")) ;
		  bsRemoveTag (Mrna, str2tag ("Length_5prime_UTR")) ;
		  bsRemoveTag (Mrna, str2tag ("Length_3prime_UTR")) ;
                  bsSave (Mrna) ;
                }
            }
          keySetDestroy (mrnas) ;
          keySetKill (productSet) ;
          keySetDestroy (productSet) ;
          return ;
        }
    }

  productSet = queryKey (gene, ">Product") ;
  maxP = keySetMax (productSet) ;
  if (maxP < 2)
    goto abort ;
  pHits = arrayCreate (maxP, PHIT) ;

  for (ii = 0 ; ii < maxP ; ii++)
    {
      product = keySet (productSet, ii) ;
      if ((Product = bsUpdate (product)))
        {
          ph = arrayp (pHits, ii, PHIT) ;
          ph->product = product ;
          bsGetKey (Product, _mRNA, &(ph->mrna)) ;
          if (bsGetKey (Product, _Peptide, &(ph->pepKey)))
            bsGetData (Product, _bsRight, _Int, &(ph->ln)) ;
          if (bsFindTag (Product, str2tag("Complete")))
            ph->complete = TRUE ;
          bsRemoveTag (Product, str2tag ("Unicity")) ;
          bsSave (Product) ;

          if ((ph->pep = peptideGet (ph->pepKey)))
            pepDecodeArray (ph->pep) ;
        }
    }
  arraySort (pHits, pHitLengthOrder) ;

  for (ii = 0 ; ii < arrayMax (pHits) ; ii++)
    {
      ph = arrayp (pHits, ii, PHIT) ;
      if (!ph->pep)
        continue ;
      cp = arrp (ph->pep, 0, char) ;
      for (jj = 0 ; jj < ii ; jj++)
        {
          qh = arrayp (pHits, jj, PHIT) ;
          if (!qh->pep)
            continue ;
          cq = arrp (qh->pep, 0, char) ;
          if (!strcmp (cp, cq))
            {
              if ((Product = bsUpdate (ph->product)))
                {
                  bsAddKey (Product, str2tag("Identical_to"), qh->product) ;
                  bsSave (Product) ;
                }
            }
          else if (strstr (cq, cp))  /* cp is included in cq */
            {
              if ((Product = bsUpdate (ph->product)))
                {
                  if (qh->complete)
                    bsAddKey (Product, str2tag("Included_in_complete"), qh->product) ;
                  else
                    bsAddKey (Product, str2tag("Included_in_partial"), qh->product) ;
                  bsSave (Product) ;
                }
            }
        }
      if ((Product = bsUpdate (ph->product)))
        {
          if (!bsFindTag (Product, str2tag ("Unicity")))
            bsAddTag (Product, str2tag ("Original")) ;
          bsSave (Product) ;
        }
    }

 abort:
  if (pHits)
    for (ii = 0 ; ii < arrayMax (pHits) ; ii++)
      {
        ph = arrayp (pHits, ii, PHIT) ;
        arrayDestroy (ph->pep) ;
      }

  arrayDestroy (pHits) ;
  keySetDestroy (productSet) ;
}  /* mrnaSaveProductHierarchy */

/*********************************************************************/

static void  mrnaSaveAceKog(KEYSET genes)
{
  KEYSET pfks = 0 ;
  KEY gene ;
  int i, ii ;
  char *cp ;
  OBJ Product = 0 ;
  vTXT buf = vtxtCreate () ;

  for (ii = 0 ; ii < keySetMax (genes) ; ii++)
    {
      gene = keySet (genes, ii) ;
      vtxtPrintf (buf, "Gene \"%s\"\n", name(gene)) ;
      vtxtPrintf (buf, "-D Homolog\n") ;
      pfks = queryKey (gene, " >product ; AKG") ;
      for (i = 0 ; i < keySetMax (pfks) ; i++)
        {
          if ((Product = bsCreate (keySet (pfks, i))))
            {
              if (bsGetData (Product, str2tag ("AceKogHuman"), _Text, &cp))
                do
                  {
                    vtxtPrintf (buf, "In_human \"%s\"\n", cp) ;
                  } while (bsGetData (Product, _bsDown, _Text, &cp)) ;
              bsDestroy (Product) ;
            }
        }
      keySetDestroy (pfks) ;
      vtxtPrintf (buf, "\n\n") ;
    }
  parseBuffer (vtxtPtr (buf), 0) ;
  vtxtDestroy (buf) ;
} /* mrnaSaveAceKog */
  
/*********************************************************************/

static void  mrnaSavePfam2Go (KEYSET genes)
{
  KEYSET goks = 0, pfks = 0, prods = 0 ;
  KEY t, go, gene, prod ;
  int i, jj, ii ;
  const char *cq ;
  vTXT buf = vtxtCreate () ;

  for (ii = 0 ; ii < keySetMax (genes) ; ii++)
    {
      gene = keySet (genes, ii) ;
      vtxtPrintf (buf, "Gene \"%s\"\n", name(gene)) ;
      vtxtPrintf (buf, "-D Pfam\n-D go_m_pfam\n-D go_b_pfam\n-D go_c_pfam\n") ;
      pfks = queryKey (gene, " >product ; best_product && good_product || very_good_product ; >pfam ;") ;
      for (i = 0 ; i < keySetMax (pfks) ; i++)
        vtxtPrintf (buf, "Pfam \"%s\"\n", name(keySet (pfks, i))) ;
      keySetDestroy (pfks) ;
      if (i)
	{
	  prods = queryKey (gene, ">product best_product && good_product || very_good_product") ;
	  for (jj = 0 ; jj < keySetMax (prods) ; jj++)
	    {
	      prod = keySet (prods, jj) ;
	      goks = queryKey (prod, ">pfam ; ! IS pox_pol* ; > go") ;
	      for (i = 0 ; i < keySetMax (goks) ; i++)
		{
		  cq = freeprotect (name(prod)) ;
		  go = keySet (goks, i) ;
		  if ((t = keyGetKey (go, str2tag("Molecular_function"))))
		    vtxtPrintf (buf, "go_m_pfam  \"%s\" %s\n", name(t), cq) ;
		  if ((t = keyGetKey (go, str2tag("Biological_process"))))
		    vtxtPrintf (buf, "go_b_pfam  \"%s\" %s\n", name(t), cq) ;
		  if ((t = keyGetKey (go, str2tag("Cellular_component"))))
		    vtxtPrintf (buf, "go_c_pfam  \"%s\" %s\n", name(t), cq) ;
		}
	      keySetDestroy (goks) ;
	    } 
	  keySetDestroy (prods) ;
	}
      vtxtPrintf (buf, "\n\n") ;
    }
  parseBuffer (vtxtPtr (buf), 0) ;
  vtxtDestroy (buf) ;
} /* mrnaSavePfam2Go */
  
/*********************************************************************/

static void  mrnaSavePsort2Go (KEYSET genes)
{
  KEYSET goks = 0 ;
  KEY go, gene ;
  int ii, i ;
  OBJ Prod = 0 ;
  char *cp, *cq ;
  BOOL ism ;
  vTXT buf = vtxtCreate () ;

 for (ii = 0 ; ii < keySetMax (genes) ; ii++)
    {
      gene = keySet (genes, ii) ;
      vtxtPrintf (buf, "Gene \"%s\"\n", name(gene)) ;
      vtxtPrintf (buf, "-D go_c_psort\n") ;
      goks = queryKey (gene, " >product best_product && good_product ; Transmembrane_domain OR Psort_title") ;
      for (i = 0 ; i < keySetMax (goks) ; i++)
        {
          go = keySet (goks, i) ;
          Prod = bsCreate (go) ;
          ism = bsFindTag (Prod, str2tag("Transmembrane_domain")) ;
          cp = 0 ;
	  cq = freeprotect (name(go)) ;
          bsGetText (Prod, str2tag("Psort_title"), &cp) ;
          if (cp)
            {
              if (!strcasecmp (cp, "Cytoplasmic"))
                vtxtPrintf (buf, "go_c_psort \"cytoplasm\" %s\n", cq) ;
              else if (!strcasecmp (cp, "Cytoskeletal"))
                  vtxtPrintf (buf, "go_c_psort \"cytoskeleton\" %s\n", cq) ;
              else if (!strcasecmp (cp, "Endoplasmic reticulum"))
                vtxtPrintf (buf, "go_c_psort \"endoplasmic reticulum%s\" %s\n"
                            , ism ? " membrane" : "",  cq) ;
              else if (!strcasecmp (cp, "secreted or extracellular"))
                 vtxtPrintf (buf, "go_c_psort \"extracellular space\" %s\n", cq) ;
              else if (!strcasecmp (cp, "Golgi"))
                vtxtPrintf (buf, "go_c_psort \"Golgi%s\" %s\n"
                            , ism ? " membrane" : " apparatus", cq) ;
              else if (!strcasecmp (cp, "Membrane"))
                vtxtPrintf (buf, "go_c_psort \"membrane\" %s\n", cq) ;
              else if (!strcasecmp (cp, "Mitochondrial"))
                vtxtPrintf (buf, "go_c_psort \"%s\" %s\n"
                            , ism ? "mitochondrial membrane" : "mitochondrion", cq) ;
              else if (!strcasecmp (cp, "Nuclear"))
                vtxtPrintf (buf, "go_c_psort \"%s\" %s\n"
                            , ism ? "nuclear membrane" : "nucleus", cq) ;
              else if (!strcasecmp (cp, "Peroxisomal"))
                vtxtPrintf (buf, "go_c_psort \"%s\" %s\n"
                            , ism ? "peroxisomal membrane" : "peroxisome", cq) ;
              else if (!strcasecmp (cp, "Plasma membrane"))
                vtxtPrintf (buf, "go_c_psort \"plasma membrane\" %s\n", cq) ;
              else if (!strcasecmp (cp, "Vacuolar"))
                vtxtPrintf (buf, "go_c_psort \"%s\" %s\n"
                            , ism ? "vacuolar membrane" : "vacuole", cq) ;
              else if (!strcasecmp (cp, "Vesicles of secretory system"))
                vtxtPrintf (buf, "go_c_psort \"secretory vesicle%s\" %s\n"
                            , ism ? " membrane" : "", cq) ;
            }
          else if (ism)
            vtxtPrintf (buf, "go_c_psort membrane %s\n", cq) ;
          bsDestroy (Prod) ;
        } 
      keySetDestroy (goks) ;
      vtxtPrintf (buf, "\n\n") ;
    }
  parseBuffer (vtxtPtr (buf), 0) ;
  vtxtDestroy (buf) ;
} /* mrnaSavePsort2Go */
  
/*********************************************************************/

static int  mrnaSaveGeneTitle (KEYSET genes)
{
  KEY gene ;
  int ii, nn = 0 ;
  char *cp ;
  vTXT txt = vtxtCreate () ; 

  extern char *swormGetGeneTitle (vTXT blkp, char *geneName) ;

  for (ii = 0 ; ii < keySetMax(genes) ; ii++)
    {
      gene = keySet (genes, ii) ;
      if ((cp = swormGetGeneTitle (txt, name(gene))))
        {
          nn++ ;
	  parseBuffer (messprintf ("Gene %s\nTitle \"%s\"\n\n", 
                                   name (gene), cp), 0) ;
        }
    }
  vtxtDestroy (txt) ;
  return nn ;
} /* mrnaSaveGeneTitle */
  
/*********************************************************************/

static void  mrnaSaveGenePastilles (KEYSET genes)
{
  KEYSET mrnas = 0, ks = 0 ;
  KEY mrna, gene ;
  int ii, jj, n, bestq ;
  OBJ Mrna = 0, Gene = 0 ;
  
  /* well supported mrna
     gene :  !cloud_gene
     product : coding_length > 300 || (coding_length > 150 && Tax_common_ancestor != *sapiens*)
     mrna : COUNT cds_covered_by = 1 
  */
  /* mark cloud and also well supported and also genecard_id */
  for (ii = 0 ; ii < keySetMax(genes) ; ii++)
    {
      gene = keySet (genes, ii) ;
      mrnas =  queryKey (gene, "{>mrna} SETOR {>mrna_part}") ;
      for (jj = 0 ; jj < keySetMax(mrnas) ; jj++)
        {
          mrna = keySet (mrnas, jj) ;  
          if ((Mrna = bsUpdate (mrna)))
	    {
	      ks = queryKey (mrna, "COUNT cds_covered_by = 1") ;
	      if (keySetMax (ks) > 0)
		bsAddTag (Mrna, str2tag ("Well_supported")) ;
	      keySetDestroy (ks) ;
	      
	      ks = queryKey (mrna, "found3p && found5p") ;
	      if (keySetMax (ks) > 0)
		bsAddTag (Mrna, str2tag ("Complete")) ;
	      keySetDestroy (ks) ;
	      
	      ks = queryKey (mrna, "from_gene && Complete && !(found3p && found5p)") ;
	      if (keySetMax (ks) > 0)
		bsRemoveTag (Mrna, str2tag ("Complete")) ;
	      keySetDestroy (ks) ;
	      
	      bsSave (Mrna) ;
	    }
        }
      keySetDestroy (mrnas) ;
      
      if ((Gene = bsUpdate (gene)))
	{
	  BOOL hasTg = bsFindTag (Gene, _Transcribed_gene) ;
	  BOOL isSpliced = FALSE ;
	  BOOL isCoding = FALSE ;

	  /* pastille_regulation_complex_locus */
	  if (bsFindTag (Gene, str2tag("Complex_locus")) &&  bsFindTag (Gene, _Transcribed_gene))
	    bsAddTag (Gene, str2tag("pastille_regulation_complex_locus")) ;
	  else
	    bsRemoveTag (Gene, str2tag("pastille_regulation_complex_locus")) ;
	  
	  /* pastille_regulation_structure */
	  if (
	      (
	       bsFindTag (Gene, str2tag("nsce")) ||
	       bsFindTag (Gene, str2tag("noce")) ||
	       bsFindTag (Gene, str2tag("nCitroenIntrons")) ||
	       bsFindTag (Gene, str2tag("nOverlappingCentralExons"))
	       )
	      &&  bsFindTag (Gene, _Transcribed_gene)
	      )
	    bsAddTag (Gene, str2tag("pastille_regulation_structure")) ;
	  else
	    bsRemoveTag (Gene, str2tag("pastille_regulation_struture")) ;
	  
	  /* pastille_regulation_RNA_editing */
	  if (
	      (
	       bsFindTag (Gene, str2tag("RNA_editing"))
	       )
	      &&  bsFindTag (Gene, _Transcribed_gene)
	      )
	    bsAddTag (Gene, str2tag("pastille_regulation_RNA_editing")) ;
	  else
	    bsRemoveTag (Gene, str2tag("pastille_regulation_RNA_editing")) ;
	  
	  /* pastille_regulation_NMD */
	  ks = queryKey (gene, " {>mrna} SETOR {>mrna_part} ; NMD") ;
	  if (keySetMax (ks))
	    bsAddTag (Gene, str2tag("pastille_regulation_NMD")) ;
	  else
	    bsRemoveTag (Gene,  str2tag("pastille_regulation_NMD")) ;
	  keySetDestroy (ks) ;
	  
	  /* pastille_regulation_Valid3p */
	  ks = queryKey (gene, "{>mrna} SETOR {>mrna_part} ; COUNT Valid3p > 1") ;
	  if (keySetMax (ks))
	    bsAddTag (Gene, str2tag("pastille_regulation_Valid3p")) ;
	  else
	    bsRemoveTag (Gene,  str2tag("pastille_regulation_Valid3p")) ;
	  keySetDestroy (ks) ;
	  
	  /* pastille_regulation_uORF */
	  ks = queryKey (gene, "{>mrna} SETOR {>mrna_part} ; COUNT {>product  uORF_candidate && !(best_product || good_product)} > 0 &&  COUNT {>product  ! uORF_candidate && (best_product || very_good_product)} > 0") ;
	  if (keySetMax (ks))
	    bsAddTag (Gene, str2tag("pastille_regulation_uORF")) ;
	  else
	    bsRemoveTag (Gene,  str2tag("pastille_regulation_uORF")) ;
	  keySetDestroy (ks) ;
	  
	  
	  /* non_atg_start */
	  ks = queryKey (gene, "{>mrna} SETOR {>mrna_part} ; > product best_product && good_product && at_position_1  && First_Kozak=1 && (First_ATG > 1 || ! First_ATG)") ;
	  if (keySetMax (ks))
	    bsAddTag (Gene, str2tag("non_atg_start")) ;
	  else
	    bsRemoveTag (Gene,  str2tag("non_atg_start")) ;
	  keySetDestroy (ks) ;
	  
	  /* pastille_conserved_pfam */
	  ks = queryKey (gene, "{>mrna} SETOR {>mrna_part} ; > product ; >PFam ; ! IS rvt_1 && ! IS transposase_* && ! IS ribosomal_* && ! IS  gag_* && ! IS  rve && ! IS  rvp && ! IS  dde && ! IS  gp36") ;
	  if (keySetMax (ks))
	    bsAddTag (Gene, str2tag("pastille_conserved_pfam")) ;
	  else
	    bsRemoveTag (Gene,  str2tag("pastille_conserved_pfam")) ;
	  keySetDestroy (ks) ;
	  
	  /* nStandardIntrons */
	  ks = queryKey (gene, ">transcribed_gene ; >Intron gt_ag || gc_ag") ;
	  n = keySetMax (ks) ;
	  if (n > 0)
	    {
	      bsAddData (Gene, str2tag("nStandardIntrons"), _Int, &n) ;
	      bsAddTag (Gene,  str2tag("Spliced_gene")) ;
	      isSpliced = TRUE ;
	    }
	  else
	    {
	      bsRemoveTag (Gene,  str2tag("nStandardIntrons")) ;
	      bsRemoveTag (Gene,  str2tag("Spliced_gene")) ;
	      isSpliced = FALSE ;
	    }
	  keySetDestroy (ks) ;
 
	  /* nStandardIntronsInModel */
	  ks = queryKey (gene, ">Genefinder ; >Intron gt_ag || gc_ag") ;
	  n = keySetMax (ks) ;
	  if (n > 0)
	    {
	      bsAddData (Gene, str2tag("nStandardIntronsInModel"), _Int, &n) ;
	    }
	  else
	    bsRemoveTag (Gene,  str2tag("nStandardIntronsInModel")) ;
	  keySetDestroy (ks) ;
	  
	  /* pastille_disease */
	  ks = queryKey (gene, "transcribed_gene ;  locus_description || locus_phenotype || gene_id == phenotype || gene_id == essential* || COUNT {{>disease} SETOR {>Mesh} SETOR {>Reference ; mesh && COUNT gene < 3 ; >Mesh  meshkey } SETOR  {>Extern;>Disease} SETOR  {>Extern ; KEGG_disease ; >pathway } SETOR  {>Extern  GAD AND NOT AntiGad ; >GAD_TITLE} SETOR  {>Extern ; OMIM_disease; >OMIM_title } ; {!Hidden_alias_of} SETOR  {>Hidden_alias_of} ; {! alias_of} SETOR  {>alias_of} ;  {IS *} SETMINUS {>meshkey ; >child ; >mesh} ; !meshkey || meshkey = C* || meshkey = F*} > 0") ;
	  if (keySetMax (ks))
	    bsAddTag (Gene, str2tag("pastille_disease")) ;
	  else
	    bsRemoveTag (Gene,  str2tag("pastille_disease")) ;
	  keySetDestroy (ks) ;
	  
	  
	  /* Single_exon_gene */
	  if (hasTg && ! isSpliced)
	    bsAddTag (Gene, str2tag("Single_exon_gene")) ;
	  else
	    bsRemoveTag (Gene,  str2tag("Single_exon_gene")) ;
	  
	  /* Best_product_quality */
	  ks = queryKey (gene, ">product ; quality") ;
	  for (jj = bestq = 0 ; jj < keySetMax(ks) ; jj++)
	    {
	      int q ;
	      OBJ Product = 0 ;
	      KEY product = keySet (ks, jj) ;  

	      if ((Product = bsCreate (product)))
		{
		  if (bsGetData (Product, str2tag("Quality"), _Int, &q) && q > bestq)
		    bestq = q ;
		  bsDestroy (Product) ;
		}
	    }
	  if (bestq > 0)
	    bsAddData (Gene, str2tag("Best_product_quality"), _Int, &bestq) ;
	  else
	    bsRemoveTag (Gene,  str2tag("Best_product_quality")) ;
	  keySetDestroy (ks) ;

	  /* Pastille_coding */
	  ks = queryKey (gene, "! Non_protein_coding && transcribed_gene ; COUNT {>product; ! From_prediction && best_product && good_product} > 0 ;  [2* COUNT {>product;! From_prediction && best_product && good_product} - COUNT {{>mrna} SETOR {>mrna_part}}] >= 0 OR  [4* COUNT {>product;! From_prediction &&  best_product && very_good_product} - COUNT {{>mrna} SETOR {>mrna_part}}] >= 0") ;
	  if (keySetMax (ks))
	    {
	      isCoding = TRUE ;
	      bsAddTag (Gene, str2tag("Pastille_coding")) ;
	    }
	  else
	    bsRemoveTag (Gene,  str2tag("Pastille_coding")) ;
	  keySetDestroy (ks) ;
	  
	  /* Pastille_marginally_coding */
	  bsRemoveTag (Gene,  str2tag("Pastille_marginally_coding")) ;
	  if (! isCoding && isSpliced)
	    {
	      ks = queryKey (gene, "! Non_protein_coding && COUNT {{>mrna} SETOR {>mrna_part} ; >product good_product && best_product} > 0") ;
	      if (keySetMax (ks))
		{
		  isCoding = TRUE ;
		  bsAddTag (Gene, str2tag("Pastille_marginally_coding")) ;
		}
	      keySetDestroy (ks) ;
	    }
	  /* Pastille_spliced_non_coding */
	  ks = queryKey (gene, "! Pastille_coding && ! Pastille_marginally_coding && spliced_gene") ;
	  if (hasTg && ! isCoding && isSpliced)
	    bsAddTag (Gene, str2tag("Pastille_spliced_non_coding")) ;
	  else
	    bsRemoveTag (Gene,  str2tag("Pastille_spliced_non_coding")) ;
	  keySetDestroy (ks) ;
	  
	  
	  /* pastille_conserved_taxblast */
	  if (0)
	    {
	      const char *species_tax_common_ancestor = "" ;
	      const char *species_tax_common_ancestor2 = "" ;
	      ks = queryKey (gene, messprintf(">product; ((best_product && good_product) || very_good_product) ; Tax_common_ancestor !=  %s && Tax_common_ancestor !=  %s" , species_tax_common_ancestor, species_tax_common_ancestor2)) ;
	      if (keySetMax (ks))
		bsAddTag (Gene, str2tag("pastille_conserved_taxblast")) ;
	      else
		bsRemoveTag (Gene,  str2tag("pastille_conserved_taxblast")) ;
	      keySetDestroy (ks) ;
	    }

	  /* balise */	
	  ks = queryKey (gene, "Omim || GeneId") ;
	  if (isSpliced || keySetMax (ks))
	    bsAddTag (Gene, str2tag("Balise")) ;
	  else
	    bsRemoveTag (Gene,  str2tag("Balise")) ;
	  keySetDestroy (ks) ;

	  /* cloud */	
	  bsRemoveTag (Gene, str2tag("Cloud_gene")) ;
	  bsRemoveTag (Gene, str2tag("Cloud")) ;
	  if (! isCoding && ! isSpliced)
	    {
	      ks = queryKey (gene, "! Balise && !Omim && !GeneId && ! Non_protein_coding &&  Nb_cDNA_clone < 5") ;
	      /* 2011_10_27: single_exon not saved by  tables/GenesSingleExonCleanUpByDeep.def  become cloud */
	      if (keySetMax (ks))
		{
		  bsRemoveTag (Gene,  str2tag("Single_exon_gene")) ;
		  bsAddTag (Gene, str2tag("Cloud_gene")) ;
		  bsAddTag (Gene, str2tag("Cloud")) ;
		} 
	      keySetDestroy (ks) ;
	    }
	  if (1)
	    {
	      int ntr, ntrcomplete, ntr5p, ntr3p, ntrpartial ;
	      int untr, untrcomplete, untr5p, untr3p, untrpartial ;
	      int gpr, gprcomplete, gprnh2, gprcooh, gprpartial ;
	      int vgpr, vgprcomplete, vgprnh2, vgprcooh, vgprpartial ;

	      ks = queryKey (gene, "{>mrna} SETOR {>mrna_part} ; gt_ag || gc_ag ") ;
	      ntr = keySetMax (ks) ; keySetDestroy (ks) ;
	      ks = queryKey (gene, "{>mrna} SETOR {>mrna_part} ; gt_ag || gc_ag ; complete") ;
	      ntrcomplete = keySetMax (ks) ; keySetDestroy (ks) ;
	      ks = queryKey (gene, "{>mrna} SETOR {>mrna_part} ; gt_ag || gc_ag ; found5p && ! found3p") ;
	      ntr5p = keySetMax (ks) ; keySetDestroy (ks) ;
	      ks = queryKey (gene, "{>mrna} SETOR {>mrna_part} ; gt_ag || gc_ag ; found3p && ! found5p") ;
	      ntr3p = keySetMax (ks) ; keySetDestroy (ks) ;
	      ks = queryKey (gene, "{>mrna} SETOR {>mrna_part} ; gt_ag || gc_ag ; !found3p && ! found5p") ;
	      ntrpartial = keySetMax (ks) ; keySetDestroy (ks) ;

	      ks = queryKey (gene, "{>mrna} SETOR {>mrna_part} ; !(gt_ag || gc_ag) ") ;
	      untr = keySetMax (ks) ; keySetDestroy (ks) ;
	      ks = queryKey (gene, "{>mrna} SETOR {>mrna_part} ;  !(gt_ag || gc_ag) ; complete") ;
	      untrcomplete = keySetMax (ks) ; keySetDestroy (ks) ;
	      ks = queryKey (gene, "{>mrna} SETOR {>mrna_part} ;  !(gt_ag || gc_ag) ; found5p && ! found3p") ;
	      untr5p = keySetMax (ks) ; keySetDestroy (ks) ;
	      ks = queryKey (gene, "{>mrna} SETOR {>mrna_part} ;  !(gt_ag || gc_ag) ; found3p && ! found5p") ;
	      untr3p = keySetMax (ks) ; keySetDestroy (ks) ;
	      ks = queryKey (gene, "{>mrna} SETOR {>mrna_part} ;  !(gt_ag || gc_ag) ; !found3p && ! found5p") ;
	      untrpartial = keySetMax (ks) ; keySetDestroy (ks) ;

	      bsRemoveTag (Gene, str2tag("nSpliced_mRNA")) ;
	      bsRemoveTag (Gene, str2tag("nUnspliced_mRNA")) ;
	      if (ntr)
		{
		  bsAddData (Gene, str2tag("nSpliced_mRNA"), _Int, &ntr) ;
		  bsAddData (Gene, _bsRight, _Int, &ntrcomplete) ;
		  bsAddData (Gene, _bsRight, _Int, &ntr5p) ;
		  bsAddData (Gene, _bsRight, _Int, &ntr3p) ;
		  bsAddData (Gene, _bsRight, _Int, &ntrpartial) ;
		}		

	      if (untr)
		{
		  bsAddData (Gene, str2tag("nUnspliced_mRNA"), _Int, &untr) ;
		  bsAddData (Gene, _bsRight, _Int, &untrcomplete) ;
		  bsAddData (Gene, _bsRight, _Int, &untr5p) ;
		  bsAddData (Gene, _bsRight, _Int, &untr3p) ;
		  bsAddData (Gene, _bsRight, _Int, &untrpartial) ;
		}		

	      ks = queryKey (gene, ">product ; good_product && best_producr && !very_good_product") ;
	      gpr = keySetMax (ks) ; keySetDestroy (ks) ;
	      ks = queryKey (gene, ">product ; good_product && best_producr && !very_good_product ; complete") ;
	      gprcomplete = keySetMax (ks) ; keySetDestroy (ks) ;
	      ks = queryKey (gene, ">product ; good_product && best_producr && !very_good_product ; NH2_complete && !COOH_complete") ;
	      gprnh2 = keySetMax (ks) ; keySetDestroy (ks) ;
	      ks = queryKey (gene, ">product ; good_product && best_producr && !very_good_product ; !NH2_complete && COOH_complete") ;
	      gprcooh = keySetMax (ks) ; keySetDestroy (ks) ;
	      ks = queryKey (gene, ">product ; good_product && best_producr && !very_good_product ; !NH2_complete && !COOH_complete") ;
	      gprpartial = keySetMax (ks) ; keySetDestroy (ks) ;

	      ks = queryKey (gene, ">product ; very_good_product") ;
	      vgpr = keySetMax (ks) ; keySetDestroy (ks) ;
	      ks = queryKey (gene, ">product ; very_good_product ; complete") ;
	      vgprcomplete = keySetMax (ks) ; keySetDestroy (ks) ;
	      ks = queryKey (gene, ">product ; very_good_product ; NH2_complete && !COOH_complete") ;
	      vgprnh2 = keySetMax (ks) ; keySetDestroy (ks) ;
	      ks = queryKey (gene, ">product ; very_good_product ; !NH2_complete && COOH_complete") ;
	      vgprcooh = keySetMax (ks) ; keySetDestroy (ks) ;
	      ks = queryKey (gene, ">product ; very_good_product ; !NH2_complete && !COOH_complete") ;
	      vgprpartial = keySetMax (ks) ; keySetDestroy (ks) ;

	      bsRemoveTag (Gene, str2tag("Good_product")) ;
	      bsRemoveTag (Gene, str2tag("Very_good_product")) ;
	      if (gpr)
		{
		  bsAddData (Gene, str2tag("Good_product"), _Int, &gpr) ;
		  bsAddData (Gene, _bsRight, _Int, &gprcomplete) ;
		  bsAddData (Gene, _bsRight, _Int, &gprnh2) ;
		  bsAddData (Gene, _bsRight, _Int, &gprcooh) ;
		  bsAddData (Gene, _bsRight, _Int, &gprpartial) ;
		}		
	      if (vgpr)
		{
		  bsAddData (Gene, str2tag("Very_good_product"), _Int, &vgpr) ;
		  bsAddData (Gene, _bsRight, _Int, &vgprcomplete) ;
		  bsAddData (Gene, _bsRight, _Int, &vgprnh2) ;
		  bsAddData (Gene, _bsRight, _Int, &vgprcooh) ;
		  bsAddData (Gene, _bsRight, _Int, &vgprpartial) ;
		}		
	    }
	  if (0)
	    {
	      /* mouse : throw to cloud the following genes
	       * single_exon_gene
	       * ! geneid   !wbid !pfam Best_product_quality<3  Best_product_quality + cover <= 4
	       */
	    }
	  bsSave (Gene) ;
	}
    }
} /* mrnaSaveGenePastilles */

/*********************************************************************/

static void  mrnaSaveAntisens (KEYSET genes)
{
  KEY gene,gene2,  tg, tg2 ;
  int ii, dx ;
  OBJ Gene = 0, Gene2 = 0, Tg = 0 ;

  for (ii = 0 ; ii < keySetMax(genes) ; ii++)
    {
      gene = keySet (genes, ii) ;
      if (! (Gene = bsUpdate (gene)))
        continue ;
      if (bsGetKey (Gene, _Transcribed_gene, &tg) &&
          (keyFindTag (tg, _gt_ag) || keyFindTag (tg, _gc_ag)) &&
          keyFindTag (tg, str2tag ("Antisens_to")) &&
          (Tg = bsCreate (tg)))
        {
          if (bsGetKey (Gene, _Transcribed_gene, &tg2) &&
              bsGetData (Gene, _bsRight, _Int, &dx) &&
              (dx >= 40) &&
              (keyFindTag (tg2, _gt_ag) || keyFindTag (tg2, _gc_ag)) &&
              (gene2 = keyGetKey (tg2, _Gene))
              )
            {
              bsAddKey (Gene, str2tag ("Antisens_to"), gene2) ;
              bsAddData (Gene, _bsRight, _Int, &dx) ;
              bsSave (Gene) ; /* needed because of the XREFs */
              if ((Gene2 = bsUpdate (gene2)))
                {
                  bsAddKey (Gene2, str2tag ("Antisens_to"), gene) ;
                  bsAddData (Gene2, _bsRight, _Int, &dx) ;
                  bsSave (Gene2) ; /* needed because of the XREFs */
                }
            }
          bsDestroy (Tg) ;
        }
      bsSave (Gene) ;
    }
} /* mrnaSaveAntisens */
  
/*********************************************************************/

static void mrnaAddAlterSpliceDetails (KEYSET genes)
{
  int ii = genes ? keySetMax (genes): 0 ;
  while (ii--)
    {
      KEY gene = keySet (genes, ii) ; 
      if (keyFindTag (gene, _IntMap))
	mcAddKeysetAlterSpliceDetails (0, gene, FALSE) ;
   }
  return ;
} /* mrnaAddAlterSpliceDetails */

/*********************************************************************/

static int mrnaAddOneMrnas (KEY gene)
{
  int nm = 0 ;

  /*  
  BOOL split = mrnaPleaseSplitByGeneId () ;
  if (1)
    {
      AC_TABLE mrnas = ac_bql_table ("", h) ;
      KEYSET mrnas = queryKey (gene, ">transcribed_gene ; >mRNA") ;
      int nn = keySetMax (mrnas) ;
      KEYSET pgs = queryKey (gene, ">transcribed_gene ; >mRNA") ;
      int nn = keySetMax (mrnas) ;
      
      keySetDestroy (mrnas) ;
      keySetDestroy (pgs) ;
    }
  */
  if (1)
    {
      KEYSET mrnas = queryKey (gene, ">transcribed_gene ; >mRNA") ;
      int nn = keySetMax (mrnas) ;
      OBJ Gene = bsUpdate (gene) ;
      
      if (Gene)
	{
	  if (nn)
	    {
	      int jj ;
	      for (jj = 0 ; jj < keySetMax (mrnas) ; jj++)
		{
		  KEY mrna = keySet (mrnas, jj) ;
		  bsAddKey (Gene, _mRNA, mrna) ;
		  nm++ ;
		}
	    }
	  else
	    bsRemoveTag (Gene, _mRNA) ;
	}
      bsSave (Gene) ;
      keySetDestroy (mrnas) ;
    }
  return nm ;
} /* mrnaAddOneMrnas */

/*********************************************************************/

static int mrnaAddMrnas (KEYSET genes)
{
  int ii = genes ? keySetMax (genes): 0 ;
  while (ii--)
    mrnaAddOneMrnas (keySet (genes, ii)) ;
  return 0 ;
} /*  mrnaAddMrnas */

/*********************************************************************/

static int mrnaCountClones (KEYSET genes, int *nGp, int *nClop)
{
  KEYSET ks = 0 ;
  Array exons = 0 ; 
  int i, ii, iii, ix, nn, nr, bp, shadow = 0, x1, x2, a1, a2, b1, b2, tg1, tg2 ;
  float cover = 0 ;
  KEY gene, nb = str2tag ("Nb_cDNA_clone") ;
  OBJ Gene = 0 ;
  Array aa = 0 ;
  BSunit *up ;
  HIT *vp ;
  if (nClop) *nClop = 0 ;
  if (nGp) *nGp = 0 ;

  if (!bsIsTagInClass (_VGene, nb))
    return 0 ;

  for (iii = 0 ; iii < keySetMax(genes) ; iii++)
    {
      gene = keySet (genes, iii) ;
      bp = nn = nr = 0 ;	 
      exons = arrayCreate (32, HIT) ;
      ix = 0 ;

      if (keyFindTag (gene, _mRNA))
	{   /*  2017_02 : new direct gene->mrna connection */
	  KEYSET mrnas = queryKey ( gene, "{>mrna} SETOR {>mrna_part}") ;
	  KEYSET reads = keySetCreate () ;
	  KEY r ;

	  int imrna =  keySetMax (mrnas) ;
	  Array aa = arrayCreate (256, BSunit) ;

	  while (imrna--)
	    {
	      KEY mrna = keySet (mrnas, imrna) ;
	      OBJ Mrna = bsCreate (mrna) ;
	      if (! Mrna)
		continue ;

	      tg1 = tg2 = 0 ;
	      if (bsGetArray (Mrna, _IntMap, aa, 3))
		for (i = 0 ; i < arrayMax (aa) ; i += 20000)
		  {
		    up = arrp (aa, i, BSunit) ;
		    tg1 = up[1].i ; tg2 = up[2].i ;
		  }
		  
	      if (bsGetArray (Mrna, _Constructed_from, aa, 3))
		for (i = 0 ; i < arrayMax (aa) ; i += 3)
		  {
		    up = arrp (aa, i, BSunit) ;
		    r = up[2].k ;
		    if ( ! keySetFind (reads, r, 0) && keyFindTag (r, _Is_read))
		      {
			keySetInsert (reads, r) ;
			bp += up[1].i - up[0].i + 1 ;
		      }
		  }

	      if (bsGetArray (Mrna, _Splicing, aa, 5))
		for (i = 0 ; i < arrayMax (aa) ; i += 5)
		  {
		    up = arrp (aa, i, BSunit) ;
		    if (strstr (name(up[4].k), "xon") > 0)
		      {
			x1 = up[0].i ; x2 = up[1].i ;
			vp = arrayp (exons, ix++, HIT) ;
			if (tg1 < tg2)
			  {
			    vp->a1 = tg1 + x1 - 1 ; 
			    vp->a2 = tg1 + x2 - 1 ; 
			  }
			else
			  {
			    vp->a1 = tg1 - x2 + 1 ; 
			    vp->a2 = tg1 - x1 + 1 ; 
			  }
		      }
		  }
	      bsDestroy (Mrna) ;
	    }
	  arrayDestroy (aa) ;

	  ks = query ( reads, "> cdna_clone") ;
	  nn = keySetMax (ks) ;
	  keySetDestroy (ks) ;
	  nr = keySetMax (reads) ;
	  keySetDestroy (reads) ;
	}
      else
	{   /*  before 2017_02 : gene->tg->mRNA connection */
	  ks = queryKey ( gene, "> transcribed_gene ; > cdna_clone") ;
	  nn = keySetMax (ks) ;
	  keySetDestroy (ks) ;

	  ks = queryKey ( gene, "> transcribed_gene ; > read") ;
	  nr = keySetMax (ks) ;
	  keySetDestroy (ks) ;
	
	  bp = 0 ;
	  aa = arrayCreate (256, BSunit) ;
	  ks = queryKey (gene, "> Transcribed_gene  Assembled_from") ;
	  for (ii = 0 ; ii < keySetMax (ks) ; ii++)
	    {
	      OBJ Tg = bsCreate (keySet (ks, ii)) ;
	      if (Tg)
		{
		  if (bsGetArray (Tg, _Assembled_from, aa, 3))
		    for (i = 0 ; i < arrayMax (aa) ; i += 3)
		      {
			up = arrp (aa, i, BSunit) ;
			bp += up[1].i - up[0].i + 1 ;
		      }
		  
		  tg1 = tg2 = 0 ;
		  if (bsGetArray (Tg, _IntMap, aa, 3))
		    for (i = 0 ; i < arrayMax (aa) ; i += 20000)
		      {
			up = arrp (aa, i, BSunit) ;
			tg1 = up[1].i ; tg2 = up[2].i ;
		      }
		  
		  if (bsGetArray (Tg, _Splicing, aa, 5))
		    for (i = 0 ; i < arrayMax (aa) ; i += 5)
		      {
			up = arrp (aa, i, BSunit) ;
			if (strstr (name(up[2].k), "xon") > 0)
			  {
			    x1 = up[0].i ; x2 = up[1].i ;
			    vp = arrayp (exons, ix++, HIT) ;
			    if (tg1 < tg2)
			      {
				vp->a1 = tg1 + x1 - 1 ; 
				vp->a2 = tg1 + x2 - 1 ; 
			      }
			    else
			      {
				vp->a1 = tg1 - x2 + 1 ; 
				vp->a2 = tg1 - x1 + 1 ; 
			      }
			  }
		      }
		  
		  bsDestroy (Tg) ;
		}
	    }
	  arrayDestroy (aa) ;
	}
      arraySort (exons,  cDNAOrderByA1) ;
      arrayCompress (exons) ;
      /* compute the area of the projection of all the exons of all tg in the gene onto the genome */
      for (shadow = 0, ix = b2 = 0, b1 = 1 ; ix < arrayMax (exons) ; ix++)
	{
	  vp = arrp (exons, ix, HIT) ;
	  a1 = vp->a1 ; a2 = vp->a2 ;
	  if (a1 <= b2 + 1)
	    { 
	      if (a2 > b2) b2 = a2; 
	    } 
	  else 
	    { 
	      shadow += b2 - b1 + 1 ; 
		  b1 = a1 ; b2 = a2 ;
	    }
	}
      shadow += b2 - b1 + 1 ; 
      arrayDestroy (exons) ;

      if ((Gene = bsUpdate (gene)))
        {
	  if (shadow > 0)
	    bsAddData (Gene, str2tag("Exon_shadow"), _Int, &shadow) ;
	  else
	    bsRemoveTag (Gene, str2tag("Exon_shadow")) ;

	  cover = shadow > 0 ? (int)(100 * bp/((float)shadow))/100.0 : 0 ;
          if (nn)
            {
              bsAddData (Gene, nb, _Int, &nn) ;
              bsAddData (Gene,_bsRight, _Int, &nr) ;
              bsAddData (Gene,_bsRight, _Int, &bp) ;
              bsAddData (Gene,_bsRight, _Float, &cover) ;
              if (nClop) *nClop += nn ;
              if (nGp) *nGp += 1 ;
            }
          else
	    bsRemoveTag (Gene, nb) ;
          bsSave (Gene) ;
        }
    } 
  return ii ;
} /* mrnaCountClones */

/*********************************************************************/
/*** Expression profile ***/
/* divise par (9.4, 13.9, 19.3, 11.9, 12.8, 20.0, 17.0, dauer */
/* divise par (22.7, 13.3, 37.1, 30.7, 13.0, 20.0, 8.6, dauer */
/* divise par (36.0, 37.1, 30.7, 13.0, 20.0, 8.6, dauer  le fuse les embyos des 2 types */
/* divise par (1.49, 1.53, 1.27, 0.54, 0.83, 0.35, 1.0=dauer june 9 */

static BOOL mrnaSaveExpressionProfile (OBJ obj, KEYSET clones, KEYSET nonSpecificClones) 
{
  KEY clone, library, _ep = str2tag ("Expression_profile") ;
  Array profile = 0 ;
  Array units = 0 ;
  OBJ Library = 0 ;
  int ii, jj, nClones = 0 ;
  BOOL ok = FALSE ;
  double total ;
  double timeNorm [] = {1.49, 1.53, 1.27, 0.54, 0.83, 0.35, 1.0, 1.0, 1.0,  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0} ; 

  if (keySetMax(clones) >= 3)
    for (jj = 0 ; jj < keySetMax(clones) ; jj++)
      {
        clone = keySet (clones, jj) ;
        if (nonSpecificClones &&
            keySetFind (nonSpecificClones, clone, 0))
          continue ;
        if (keyFindTag (clone, str2tag("Duplicate_clone")))
          continue ;
        if ((library = keyGetKey (clone, _Library)) &&
            keyFindTag (library, _ep) &&
            (Library = bsCreate (library)))
          {
            /* add up the profiles */
            nClones++ ;
            units = arrayReCreate (units, 8, BSunit) ;
            if (bsGetArray (Library, _ep, units, 7))
              {
                if (!profile)
                  profile = arrayCreate (12, BSunit) ;
                for (ii = 0 ; ii < arrayMax(units) ; ii++)
                  array (profile, ii, BSunit).i += (int)((1000.0 * (double)(arr (units, ii, BSunit).i) / timeNorm[ii]) );
              }            
            bsDestroy (Library) ;
          }
      }

  /* register */     
  bsRemoveTag (obj, _ep) ;
  
  if (nClones >= 3 && profile)
    {
      for (total = ii = 0 ; ii < arrayMax(profile) ; ii++)
        total += arr(profile, ii, BSunit).i ;
      if (total == 0) total = 1 ;
      for (ii = 0 ; ii < arrayMax(profile) ; ii++)
        arr(profile, ii,BSunit).i = (100.0 * (double)arr(profile, ii,BSunit).i + (total/2) - 1)/total ;
      /* +total/2 - 1 makes the rounding closer to average to zero, so the sum remains 100% */
      ii = arrayMax(profile) ;
      if (ii > 0) 
        bsAddArray (obj, _ep, profile, ii) ;
      arrayDestroy (profile) ;
      ok = TRUE ;
    }
  arrayDestroy (units) ;
  return ok ;
} /* mrnaSaveExpressionProfile */

/*********************************************************************/

static BOOL  mrnaGeneRegisterIntronSupport (OBJ Gene, Array feet, HIT *fp, HIT *vp)
{
  OBJ Est = 0 ;   /* success this est is responsible for this stupid intron */
  HIT *f1 ;
  int i ;

  if (vp->type & gMicro)
    return TRUE ;
  if (fp->cDNA_clone == _Fuzzy && (vp->type & gI))
    {
      for (i = 0, f1 = arrp (feet, i, HIT) ; i < arrayMax (feet) ; i++, f1++)
        if (f1->a1 == vp->a1 && f1->a2 == vp->a2 &&
            f1->cDNA_clone != _Fuzzy)
          return FALSE ;
    }
  if (1) /* (fp->cDNA_clone != _gt_ag && fp->cDNA_clone != _gc_ag) */
    {
      /* these are already in Gene, but we must relocalize ourself */
      bsAddTag (Gene, fp->cDNA_clone) ;
      if (fp->cDNA_clone == _Other)
        bsAddData (Gene, _bsRight, _Text, name(fp->est)) ;
      bsAddData (Gene, _bsRight, _Int, &(fp->x1)) ;
      
      bsAddData (Gene, _bsRight, _Int, &(fp->a1)) ;
      bsAddData (Gene, _bsRight, _Int, &(fp->a2)) ;
      
      /* now at last we can add the clone name */
      bsAddKey (Gene, _bsRight, vp->cDNA_clone) ; 
      if (fp->cDNA_clone == _Fuzzy)
        bsAddData (Gene, _bsRight, _Text, name(fp->est)) ;
      
    }
  if ((Est = bsUpdate (vp->est)))
    {
      int nold = 0 ;
      bsAddTag (Est, fp->cDNA_clone) ; 
      if (fp->cDNA_clone == _Other)
        bsAddData (Est, _bsRight, _Text, name(fp->est)) ;
      bsGetData (Est, _bsRight, _Int, &nold) ;
      nold++ ;
      /* reposition in Est */
      bsAddTag (Est, fp->cDNA_clone) ; 
      if (fp->cDNA_clone == _Other)
        bsAddData (Est, _bsRight, _Text, name(fp->est))  ;
      bsAddData (Est, _bsRight, _Int, &nold) ;
      bsSave (Est) ;
    }
  return TRUE ;
} /* mrnaGeneRegisterIntronSupport */

/*********************************************************************/

static void mrnaSaveGene (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas)
{
  KEY gene = sc->gene ;
  OBJ Gene = 0 ;
  KEY map = 0, cosmid = s2m->cosmid1 ? s2m->cosmid1 : s2m->cosmid ;
  OBJ Cosmid = 0 ;
  int ii, jj, c1 = 0, c2 = 0, d1, d2, mMax = arrayMax (smrnas) ;
  int ga1 = sc->a1, ga2 = sc->a2 ;
  KEY  relevantCosmid = mrnaRelevantCosmid (s2m->cosmid, s2m->cosmid1, s2m->linkPos, ga1, ga2) ;
  KEY geneBox = 0 ;
  Array feet = 0 ;
  int  minIntronSize = getPleaseMinIntronSize () ;
    
  /*** clean up ***/ 
  cdnaCleanUp  (0, gene, 0) ; /* calls   mrnaCleanUp  (0, gene, 0) ; */
  
  /*** Comid info: before we bsUpdate gene ***/
  
  { 
    Cosmid = bsUpdate (relevantCosmid) ;
    if (Cosmid)
      {
        bsAddKey (Cosmid, _Transcribed_gene, gene) ;
        bsAddData (Cosmid, _bsRight, _Int, &ga1) ;
        bsAddData (Cosmid, _bsRight, _Int, &ga2) ;
        bsSave (Cosmid) ; 
      }
  }
  
  /****************/

  Gene = bsUpdate (gene) ;
  if (!Gene)
    return ;
  
  /*** clean up ***/
  {
    bsRemoveTag (Gene, _Splicing) ;
    bsRemoveTag (Gene, _Assembled_from) ;
    bsRemoveTag (Gene, _cDNA_clone) ;
    bsRemoveTag (Gene, _Structure) ;
  }
  
  /*** Transcribed_from ***/
  {
    int ll ;
        
    if (ga2 > ga1) ll = ga2 - ga1 + 1 ;
    else ll = ga1 - ga2 + 1 ;
    
    bsAddData (Gene, _Covers, _Int, &ll) ;
    bsAddData (Gene, _bsRight, _Text, "bp from") ;
    bsAddData (Gene, _bsRight, _Int, &ga1) ;
    bsAddData (Gene, _bsRight, _Int, &ga2) ;
        
    
    if ((Cosmid = bsCreate (cosmid)))
      {
        if (bsGetKey (Cosmid, _IntMap, &map) &&
            bsGetData (Cosmid, _bsRight, _Int, &c1)  &&
            bsGetData (Cosmid, _bsRight, _Int, &c2))
	  {}  ;
        bsDestroy (Cosmid) ;
      }
    if (map && c1 + c2)
      {
        if (c1 < c2)
          { d1 = c1 + ga1 - 1 ; d2 = c1 + ga2 - 1 ; }
        else
          { d1 = c1 - ga1 + 1 ; d2 = c1 - ga2 + 1; }
        bsAddKey (Gene, _IntMap, map) ;
        bsAddData (Gene, _bsRight, _Int, &d1) ;
        bsAddData (Gene, _bsRight, _Int, &d2) ;
      }
        
    { /* GeneBox && new_name */
      KEY newName = 0 ;
          
      if (lexReClass (sc->gene, &geneBox, _VGene))
        {
          OBJ GeneBox ;

          if ((GeneBox = bsUpdate (geneBox)))
            {
              bsRemoveTag (GeneBox, _Product) ;
              bsSave (GeneBox) ;
            }
          bsAddKey (Gene, str2tag("Gene"), geneBox) ;
          if ((newName = keyGetKey (geneBox, _NewName)))
            bsAddKey (Gene, _NewName, newName) ;
        }
    }
    
    { /* Locus */
      KEY locus = 0, gLocus = 0 ;
      int i ;
      KEYSET loci = query (sc->sh.clones, "FOLLOW Read ; Follow Locus") ;
          
      for (i = 0 ; i < keySetMax(loci) ; i++)
        {
          locus = keySet (loci, i) ;
          bsAddKey (Gene, _Locus, locus) ;
          if ((gLocus = keyGetKey (locus, _Gene)))
            bsAddKey (Gene, str2tag("Gene"), gLocus) ; 
        }
      keySetDestroy (loci) ;

      if (sc->sh.ignoredClones)
        {
          loci = query (sc->sh.ignoredClones, "FOLLOW Read ; Follow Locus") ;
                  
          for (i = 0 ; i < keySetMax(loci) ; i++)
            {
              locus = keySet (loci, i) ;
              bsAddKey (Gene, _Locus, locus) ;
              /* but do not link the tg to the genebox otherwise tha by geometry */
              if (0 && (gLocus = keyGetKey (locus, _Gene)))
                bsAddKey (Gene, str2tag("Gene"), gLocus) ; 
            }
          keySetDestroy (loci) ;
        }
    }
  
    if (bsIsTagInClass (_VTranscribed_gene, str2tag("IST_as_fish")))
      { /* IST */
        KEY ist = 0 ;
        int i ;
        KEYSET ists = query (sc->sh.clones, "FOLLOW correspondingIST") ;
        
        for (i = 0 ; i < keySetMax(ists) ; i++)
          {
            ist = keySet (ists, i) ;
            bsAddKey (Gene, str2tag("IST_as_fish"), ist) ;
          }
        keySetDestroy (ists) ;
      }
  }
  
  /*** Transcripts ***/
  for (ii = jj = 0 ; ii < arrayMax(smrnas) ; ii++)
    {
      SMRNA *smrna = arrp (smrnas, ii, SMRNA) ;
          
      if (smrna->bestDna < 0 || smrna->cGroup == 2)
        mMax-- ;
    }
  
  for (ii = jj = 0 ; ii < arrayMax(smrnas) ; ii++)
    {
      KEY tr ;
      SMRNA *smrna = arrp (smrnas, ii, SMRNA) ; 
      OBJ Mrna = 0 ; 
      if (smrna->bestDna < 0 || smrna->cGroup == 2)
        continue ;
      if (FALSE && mMax == 1) /* false, apr 7, 2005, otherwise product 2 of unc-3 single mrna becomes unc-32 */
        lexaddkey (name(sc->gene), &tr, _VmRNA) ;
      else
        {   /* Dec03 */
          if (jj < 21) /* 21 letters a - u */
            lexaddkey (messprintf ("%s.%c", name(sc->gene), 'a' + jj), &tr, _VmRNA) ;
          else if (jj < 21 + 26) /* va - vz */
            lexaddkey (messprintf ("%s.v%c", name(sc->gene), 'a' + jj - 21), &tr, _VmRNA) ;
          else if (jj < 21 + 2*26) /* wa - wz */
            lexaddkey (messprintf ("%s.w%c", name(sc->gene), 'a' + jj - 21 - 26), &tr, _VmRNA) ;
          else if (jj < 21 + 3*26) /* xa - xz */
            lexaddkey (messprintf ("%s.x%c", name(sc->gene), 'a' + jj - 21 - 2 * 26), &tr, _VmRNA) ;
          else if (jj < 21 + 4*26) /* ya - yz */
            lexaddkey (messprintf ("%s.y%c", name(sc->gene), 'a' + jj - 21 - 3 * 26), &tr, _VmRNA) ;
          else if (jj < 21 + 5*26) /* za - zz */
            lexaddkey (messprintf ("%s.z%c", name(sc->gene), 'a' + jj - 21 - 4 * 26), &tr, _VmRNA) ;
          else  /* need z behind number because prods have a number */
            lexaddkey (messprintf ("%s.zz_%dz", name(sc->gene), jj + 1), &tr, _VmRNA) ;
        }
      jj++ ;
      bsSave (Gene) ;
      if ((Mrna = bsUpdate (tr)))
        bsKill (Mrna) ;
      Gene = bsUpdate (gene) ;
      bsAddKey (Gene, _mRNA, tr) ;
      bsAddData (Gene, _bsRight, _Int, &(smrna->a1)) ;
      bsAddData (Gene, _bsRight, _Int, &(smrna->a2)) ;
      smrna->gene = gene ;
      smrna->tr = tr ;
    }
  
  
  /*** status ***/
  
  {
    mytime_t  now = timeNow () ;
    bsAddData (Gene, _Date, _DateType, &now) ;
  }
  
  /*** Splicing and Structure ***/

  if (arrayMax (gmrna->hits))
  {
    int nExons, nAexons, nIntrons, nAintrons, gapl, ll, 
      ifeet = 0, ngt_ag = 0, ngc_ag = 0, nv_rep = 0, nfuzzy = 0 ;
    HIT *up, *vp, *wp, *fp ;
    KEY _v_rep = str2tag ("V_repeat"), actualFoot = 0 ;
    KEY _Fuzzy_gt_ag = str2tag ("Fuzzy_gt_ag") ;
    KEY _Fuzzy_gc_ag = str2tag ("Fuzzy_gc_ag") ;
    nExons = nAexons = nIntrons = nAintrons = gapl = 0 ;
    
    for (ii = 0, up = arrp (gmrna->hits, 0, HIT) ; ii < arrayMax (gmrna->hits) ; up++, ii++)
      {
        bsAddData (Gene, _Splicing, _Int, &(up->a1)) ;
	wp = up ;
	for (vp = up + 1, jj = ii + 1 ; jj < arrayMax (gmrna->hits) ; vp++, jj++)
	  if ((up->type & gX) && vp->a1 == up->a2 + 1 && vp->type == up->type)
	    wp = vp ;
	  else
	    break ;
        bsAddData (Gene, _bsRight, _Int, &(wp->a2)) ;
        bsPushObj (Gene) ; 
        if ((up->type & gI) && minIntronSize && wp->a2 >= up->a1 && wp->a2 - up->a1 < minIntronSize)
          up->type |= gMicro ;
        if (up->type & gX)
          {
	    char buf[20] ;
            nExons++ ;
            if (up->type & gB) 
              {
                nAexons++ ;
                bsAddKey (Gene, _bsRight, _Alternative_exon) ;
              }
            else
              {
                bsAddKey (Gene, _bsRight, _Exon) ;
              }
            ll = wp->a2 - up->a1 + 1 ;
	    sprintf (buf, "%d", nExons) ;
            bsAddData (Gene, _bsRight, _Text, buf) ;
            bsAddData (Gene, _bsRight, _Text, "Length") ;
            bsAddData (Gene, _bsRight, _Int, &ll) ;
            bsAddData (Gene, _bsRight, _Text, "bp") ;
	    up = wp ;
          }
        else if (up->type & gI)
          {
            int u1, u2 ;
            KEY foot = 0 ;
                        
            if (sc->isUp)
              { u1 = sc->d1 - up->a1 + 1 ; u2 = sc->d1 - up->a2 + 1 ; }
            else
              { u1 = sc->d1 + up->a1 - 1 ; u2 = sc->d1 + up->a2 - 1 ; }
                        
            ll = up->a2 - up->a1 + 1 ;
            if (up->type & gMicro)
              {
                actualFoot = mrnaIntronBoundary (s2m->dnaD, s2m->dnaR, u1, u2) ;
                foot = 0 ;
              }
            else if (ll < 0)
              { actualFoot = foot = _v_rep ; nv_rep++ ; }
            else if (up->type & gJ)
              {
                actualFoot = mrnaIntronBoundary (s2m->dnaD, s2m->dnaR, u1, u2) ;
                foot = _Fuzzy ; 
                if (0 && actualFoot == _gt_ag) actualFoot = _Fuzzy_gt_ag ;
                if (0 && actualFoot == _gc_ag) actualFoot = _Fuzzy_gc_ag ;
                nfuzzy++ ;
              }
            else if (ll > 0)
              {
                foot = actualFoot = mrnaIntronBoundary (s2m->dnaD, s2m->dnaR, u1, u2) ;
                if (foot == _gt_ag) ngt_ag++ ;
                if (foot == _gc_ag) ngc_ag++ ;
              } 
            if (foot)
              {
                if (!feet) feet = arrayCreate (30, HIT) ;
                fp = arrayp (feet, ifeet++, HIT) ;
                fp->a1 = up->a1 ; fp->a2 = up->a2 ; 
                if (bsIsTagInClass (_VTranscribed_gene, foot))
                  fp->cDNA_clone = foot ; 
                else
                  fp->cDNA_clone = _Other ;
                fp->est = actualFoot ;
              }

            if (up->type & (gMicro))
              {
                bsAddKey (Gene, _bsRight, str2tag ("Small_deletion")) ;
              }
            else if (up->type & (gPredicted | gStolen))
              {
                nIntrons++ ;
                bsAddKey (Gene, _bsRight, _Stolen_intron) ;
              }
            else if (up->type & gB)
              {
                nIntrons++ ;
                nAintrons++ ; 
                bsAddKey (Gene, _bsRight, _Alternative_intron) ;
              }
            else
              {   
                                /* intl += a2 - a1 + 1 ; */
                nIntrons++ ;
                bsAddKey (Gene, _bsRight, _Intron) ;
              }
            ll = up->a2 - up->a1 + 1 ;
            if (foot != _Fuzzy)
              bsAddData (Gene, _bsRight, _Text, name(actualFoot)) ;
            else
              bsAddData (Gene, _bsRight, _Text, messprintf ("Fuzzy_%s", name(actualFoot))) ;
            bsAddData (Gene, _bsRight, _Text, "Length") ;
            bsAddData (Gene, _bsRight, _Int, &ll) ;
            bsAddData (Gene, _bsRight, _Text, "bp") ;
          }
        else if (up->type & gGap)
          {
            ll = (up+1)->x1 - (up-1)->x2 + 1 ; 
            gapl +=  ll ;
            bsAddKey (Gene, _bsRight, _Gap) ;
            
            /* this length is a genomic length, it is utterly wrong
            bsAddData (Gene, _bsRight, _Text, "Unsequenced") ;
            bsAddData (Gene, _bsRight, _Text, "Length") ;
            bsAddData (Gene, _bsRight, _Int, &ll) ;
            bsAddData (Gene, _bsRight, _Text, "bp") ;
            */
          }
        bsGoto (Gene, 0) ;
      }
        
    if (nExons > 0)
      bsAddData (Gene, _Nb_possible_exons, _Int, &nExons) ;
    if (nAexons > 0)
      bsAddData (Gene, _Nb_alternative_exons, _Int, &nAexons) ;
    if (nIntrons > 0)
      bsAddData (Gene, _Nb_confirmed_introns, _Int, &nIntrons) ;
    if (nAintrons > 0)
      bsAddData (Gene, _Nb_confirmed_alternative_introns, _Int, &nAintrons) ;
    if (gapl > 0)
      bsAddData (Gene, _Total_gap_length, _Int, &gapl) ;
        
    bsRemoveTag (Gene, _Intron_boundaries) ;
    if (ngt_ag)
      bsAddData (Gene,  _gt_ag, _Int, &ngt_ag) ;
    if (ngc_ag)
      bsAddData (Gene,  _gc_ag, _Int, &ngc_ag) ;
    if (nfuzzy)
      bsAddData (Gene,  _Fuzzy, _Int, &nfuzzy) ;
    if (nv_rep)
      bsAddData (Gene,  _v_rep, _Int, &nv_rep) ;
    if (feet)
      {
        int i1, i2, i3, iMax = arrayMax(feet) ;
        HIT *fp2 ;
                
        cDNASwapA (feet) ;
        arraySort (feet, cDNAOrderByA1) ;
                
        for (i1 = 0 ; i1 < iMax ; i1++)
          {
            fp = arrp (feet, i1, HIT) ;
            if (fp->cDNA_clone == _Other &&
                ! fp->est) /* this is a Fuzzy */
              continue ;
            i3 = 0 ;
            for (i2 = i1, fp2 = fp ; 
                 fp2->cDNA_clone == fp->cDNA_clone && fp2->est == fp->est && i2 < iMax ;
                 i2++, fp2++) 
              i3++ ;  /* get the number, first thing to report */
            if (1) /* (fp->cDNA_clone != _gt_ag && fp->cDNA_clone != _gc_ag) */
              {
                for (i2 = i1, fp2 = fp ; fp2->cDNA_clone == fp->cDNA_clone && fp2->est == fp->est && i2 < iMax ; i2++, fp2++) 
                  {
                    fp2->x1 = i3 ; /* for later reference */
                    bsAddTag (Gene, fp2->cDNA_clone) ;
                    if (fp2->cDNA_clone == _Other)
                      bsAddData (Gene, _bsRight, _Text, name(fp2->est)) ;
                    bsAddData (Gene, _bsRight, _Int, &i3) ; 
                    if (1) /* fp2->cDNA_clone != _gt_ag && fp2->cDNA_clone != _gc_ag) */
                      {
                        bsAddData (Gene, _bsRight, _Int, &(fp2->a1)) ;
                        bsAddData (Gene, _bsRight, _Int, &(fp2->a2)) ;
                      }
                  }
              }
            else
              { 
                fp->x1 = i3 ; /* for later reference */
                bsAddTag (Gene, fp->cDNA_clone) ;
                bsAddData (Gene, _bsRight, _Int, &i3) ; 
              }
            i1 = i2 - 1 ;
          }
      }
  }

  /*** clone_group ***/
  {
    KEYSET ks = sc->sh.clones ;
    KEY clone, group, _Clone_group = str2tag ("Clone_group") ;
    Array groups = 0 ;
    HIT *gh ;

    bsRemoveTag (Gene, _Clone_group) ;
    if (ks)
      {
        int i, j ;
                
        for (i = 0 ; i < keySetMax(ks) ; i++)
          {
            clone = keySet (ks, i) ;
            group = keyGetKey (clone, _Clone_group) ;
            if (group)
              {
                if (!groups)
                  groups = arrayCreate (keySetMax(ks), HIT) ;
                for (j = 0 ; group && j < arrayMax(groups) ; j++)
                  {
                    gh = arrp (groups, j, HIT) ;
                    if (gh->cDNA_clone == group) 
                      {group = 0 ; (gh->a1)++ ; }
                  }
                if (group)
                  { 
                    gh = arrayp (groups, arrayMax(groups), HIT) ;
                    gh->cDNA_clone = group ;
                    gh->a1 = 1 ;
                  }
              }
          }
            
        if (groups)
          {
            arraySort (groups, cDNAOrderByA1) ;
            for (j = 0 ; j < arrayMax(groups) ; j++)
              {
                gh = arrp (groups, j, HIT) ;
                bsAddKey (Gene, _Clone_group, gh->cDNA_clone) ;
                bsAddData (Gene, _bsRight, _Int, &(gh->a1)) ;
              }
          }
                
        arrayDestroy (groups) ;
      }
  }
  
  /*** cdna_clones ***/
  {
    KEYSET ks = sc->sh.clones ;
    KEY clone ;
        
    if (ks)
      {
        int i ;
        arraySort (ks, keySetAlphaOrder) ;
        for (i = 0 ; i < keySetMax(ks) ; i++)
          {
            clone = keySet (ks, i) ;
            bsAddKey (Gene, _cDNA_clone, clone) ;
            if (keyFindTag (clone, _Image))
              bsAddKey (Gene, _Image, clone) ;
          }
        keySetSort (ks) ; /* reestablish */
      }
    ks = sc->sh.ignoredClones ;
    if (ks)
      {
        int i ;
        arraySort (ks, keySetAlphaOrder) ;
        for (i = 0 ; i < keySetMax(ks) ; i++)
          {
            clone = keySet (ks, i) ;
            bsAddKey (Gene, str2tag("Ignored_clone"), clone) ;
          }
        keySetSort (ks) ; /* reestablish */
      }
  }
  
        
  { 
    OBJ GeneBox = 0 ;
    BOOL expPat = mrnaSaveExpressionProfile (Gene, sc->sh.clones, 0) ;
    
    if (geneBox &&   /* GeneBox should get the exp pattern flag */
        (GeneBox = bsUpdate (geneBox)))
      {
        if (expPat)
          bsAddTag (GeneBox, str2tag ("Expression_profile")) ;
        else
	  bsRemoveTag (GeneBox, str2tag ("Expression_profile")) ;
        bsSave (GeneBox) ;
      }
  }

  /*** clean intron boundaries ***/

  {
    KEYSET ks = query (sc->sh.clones, ">Read Intron_boundaries") ;
    KEY est ;
    OBJ Est = 0 ;
    int i ;
        
    for (i = 0 ; i < keySetMax(ks) ; i++)
      {
        est = keySet (ks, i) ;
        if ((Est = bsUpdate (est)))
          {
            bsRemoveTag (Est, _Intron_boundaries) ;
            bsRemoveTag (Est,  str2tag("Gap_length")) ;
            bsSave (Est) ;
          }
      }
    keySetDestroy (ks) ;
        
    if (sc->sh.ignoredClones)
      { 
        ks = query (sc->sh.ignoredClones, ">Read Intron_boundaries") ;
        for (i = 0 ; i < keySetMax(ks) ; i++)
          {
            est = keySet (ks, i) ;
            if ((Est = bsUpdate (est)))
              {
                bsRemoveTag (Est, _Intron_boundaries) ;
                bsSave (Est) ;
              }
          }
        keySetDestroy (ks) ;
      }
  }
  
  /* gaps in reads */
  {
    HIT *up, *vp ;
    int gap, totalGap = 0 ;
    Array hits = gmrna->estHits ;
    OBJ Est = 0 ;

    cDNASwapX (hits) ;
    arraySort (hits, cDNAOrderByX1) ;

    for (ii = 0, vp = 0, up = arrp (hits, 0, HIT) ; ii < arrayMax(hits) ; ii++, up++)
      {
        if (!class(up->est)) /* ghosts */
          continue ; 
        if (!(up->type & gX))
          continue ;
        if (vp && up->est == vp->est)
          {
            if (vp->x2 < up->x1 - 3 || vp->x2 > up->x1 + 3)
              { 
                gap = up->x1 - vp->x2 - 1 ;
                totalGap += gap < 0 ? 1 : gap ;
              }
          }
        else
          {
            if (vp && totalGap && (Est = bsUpdate (vp->est)))
              {
                bsAddData (Est, str2tag("Gap_length"), _Int, &totalGap) ;
                bsSave (Est) ;
              }
            totalGap = 0 ;
          }
        vp = up ;
      }
  }

  /*** Assembled_from ***/
  {
    HIT *up, *fp, *vp, *wp ;
    int i1, jj ;
    Array hits = gmrna->estHits ;

    cDNASwapA (hits) ;
    arraySort (hits, cDNAOrderGloballyByA1Errors) ;
    for (ii = 0, up = arrp (hits, 0, HIT) ; ii < arrayMax(hits) ; ii++, up++)
      {
        if (!class(up->est)) /* ghosts */
          continue ; 
        if (!(up->type & gX))
          continue ;
        if (sc->sh.ignoredClones && 
            keySetFind (sc->sh.ignoredClones, up->cDNA_clone, 0))
          continue ;
        bsAddData (Gene, _Assembled_from, _Int, &(up->a1)) ;
        bsAddData (Gene, _bsRight, _Int, &(up->a2)) ;
        bsAddKey (Gene, _bsRight, up->est) ;
        bsAddData (Gene, _bsRight, _Int, &(up->x1)) ;
        bsAddData (Gene, _bsRight, _Int, &(up->x2)) ;
        bsAddData (Gene, _bsRight, _Int, &(up->nerrAll)) ;
      }

    for (i1 = 0 ; feet && i1 < arrayMax(feet) ; i1++)
      {
        BOOL ok = FALSE ;

        fp = arrp (feet, i1, HIT) ;
        if (fp->cDNA_clone == _Other &&
            ! fp->est) /* this is a Fuzzy */
          continue ;
        for (ii = 0, vp = arrp (hits, ii, HIT) ; ii < arrayMax(hits) ; ii++, vp++)
          if ((vp->type & gI) && 
              (
               (vp->a1 == fp->a1  && vp->a2 == fp->a2) ||
               (fp->cDNA_clone == _Fuzzy && (vp->type & (gJ|gGap)) &&
                vp->a1 <= fp->a1 + 5 && vp->a1 >= fp->a1 - 5  &&
                vp->a2 <= fp->a2 + 5 && vp->a2 >= fp->a2 - 5)
               )
              )
            ok = mrnaGeneRegisterIntronSupport (Gene, feet, fp, vp) ;          
        if (!ok &&
            fp->cDNA_clone == _Fuzzy &&
            fp->cDNA_clone != _gt_ag &&
            fp->cDNA_clone != _gc_ag
            ) /* we have not found a good support
                  * search for it per blue mrna 
                  */
          {
            /* yk1586a7 */  
            if (!ok)
              for (ii = 0, vp = arrp (hits, ii, HIT) ; ii + 1 < arrayMax(hits) ; ii++, vp++)
                if ((vp->type & gI) && 
                    (
                     fp->cDNA_clone == _Fuzzy && (vp->type & (gJ|gGap) &&
                      vp->a1 <= fp->a1 + 15 && vp->a1 >= fp->a1 - 15  &&
                      vp->a2 <= fp->a2 + 15 && vp->a2 >= fp->a2 - 15)
                      )
                     )
                  ok = mrnaGeneRegisterIntronSupport (Gene, feet, fp, vp) ;          
            if (!ok)
              for (ii = 0, vp = arrp (hits, ii, HIT) ; ii + 1 < arrayMax(hits) ; ii++, vp++)
                if ((vp->type & gX) &&
                    vp->a2 > fp->a1 - 16 && vp->a2 < fp->a1 + 16 
                    )
                  {
                    for (jj = ii+1, wp = vp+1 ; !ok && jj < arrayMax(hits) ; jj++, wp++)
                      if ((wp->type & gX) && vp->est == wp->est &&
                          wp->a1 > fp->a2 - 16 && wp->a1 < fp->a2 + 16
                          )
                        ok = mrnaGeneRegisterIntronSupport (Gene, feet, fp, vp) ;
                      else if ((wp->type & gI) && vp->est == wp->est)
                        break ;
                  }
            if (!ok)
              for (ii = 0, vp = arrp (hits, ii, HIT) ; ii < arrayMax(hits) ; ii++, vp++)
                if ((vp->type & gI) && 
                    (
                     (vp->a1 == fp->a1  || vp->a2 == fp->a2) &&
                     (fp->cDNA_clone == _Fuzzy && (vp->type & (gJ|gGap)) &&
                      vp->a1 <= fp->a1 + 25 && vp->a1 >= fp->a1 - 25  &&
                      vp->a2 <= fp->a2 + 25 && vp->a2 >= fp->a2 - 25)
                     )
                    )
                  ok = mrnaGeneRegisterIntronSupport (Gene, feet, fp, vp) ;
            if (!ok)
              for (ii = 0, vp = arrp (hits, ii, HIT) ; ii < arrayMax(hits) ; ii++, vp++)
                if ((vp->type & gX) && 
                    vp->a2 <= fp->a1 + 5 && vp->a2 >= fp->a1 - 5)
                  {
                    for (jj = ii+1, wp = vp+1 ; !ok && jj < arrayMax(hits) ; jj++, wp++)
                      if ((wp->type & gX) && vp->cDNA_clone == wp->cDNA_clone &&
                          wp->a1 <= fp->a2 + 5 && wp->a1 >= fp->a2 - 5
                          )
                        ok = mrnaGeneRegisterIntronSupport (Gene, feet, fp, vp) ;
                  }
          }
      }
  }
  
  bsSave (Gene) ;
  
  /*** quality: after bsSave(Gene)  ****/
  {
    mrnaSaveMatchQuality (gene, gmrna->estHits, sc->sh.ignoredClones) ;    
  }
  
  /*** clone group info: after bsSave(Gene)  ***/
  {
    OBJ CGroup = 0 ;
    
    if (sc->cGroup && (CGroup = bsUpdate (sc->cGroup)))
      {
        bsAddKey (CGroup, _Transcribed_gene, gene) ;
        bsAddKey (CGroup, _bsRight, relevantCosmid) ;
        bsSave (CGroup) ;
      }
  }
  
  arrayDestroy (feet) ;
} /* mrnaSaveGene */

/*********************************************************************/

static int mSplOrder (const void *a, const void *b)
{
  const HIT *ha = (const HIT *)a, *hb = (const HIT *)b ;
  int n ;

  n = ha->x1 - hb->x1 ; if (n) return n ;
  n = ha->x2 - hb->x2 ; if (n) return n ;
  return ha->a1 - hb->a1 ;
}

/*********************************************************************/
/* revise shifted data */
static void mrnaSaveGeneKillNonBest (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas)
{
  vTXT txt = 0 ;
  KEY gene = sc->gene ;
  KEY map = 0, cosmid = s2m->cosmid1 ? s2m->cosmid1 : s2m->cosmid ;
  int x, i, c1 = 0, c2 = 0, d1, d2 ;
  int ga1 = sc->a1, ga2 = sc->a2, gb1, gb2, gx1, gx2 ;
  KEY  relevantCosmid = 0 ;
  BSunit *uu, *vv ;
  BOOL isGeneDown = TRUE ;
  OBJ Gene = 0, Mrna = 0 ;
  Array units = 0, vnits = 0 ;
  BSMARK mark = 0 ;
  HIT *hh ;
  AC_HANDLE handle = handleCreate () ;

  /*** look for mrnas to be killed ***/ 

  {
    KEYSET ks = queryKey (gene, ">mRNA ; !Best_in_gene") ;
    int nn = keySetMax(ks) ;

    if (nn)
      {
        keySetKill (ks) ;
        keySetDestroy (ks) ;
      }
    else
      {
        keySetDestroy (ks) ;
        goto done ;
      }
  }

  /****************/

  Gene = bsUpdate (gene) ;
  if (!Gene)
    goto done ;
  
  bsGetKey (Gene, str2tag ("Genomic_sequence"), &relevantCosmid) ;
  /*** clean up ***/ 

  isGeneDown = ga1 < ga2 ? TRUE : FALSE ;
  units = arrayCreate (300, BSunit) ;
  vnits = arrayCreate (300, BSunit) ;
  /* get the active length */
  gx1 = gx2 = 0 ;
  /*** Transcripts ***/
  txt = vtxtCreate () ;

  {
    int j, jj, is, dx ;
    Array mSpl = arrayCreate (256, HIT) ;
    DICT *dict = dictCreate (256) ;

    /* find the top of the first mrna and verify that each mrna 
       has splicing starting at position 1 
    */
    if (bsGetArray (Gene, _mRNA, units, 3))
      {
        for (i = is = 0 ; i < arrayMax (units) ; i += 3)
          {
            uu = arrp (units, i, BSunit) ;
            if ((Mrna = bsCreate (uu[0].k)))
              {
                if (bsGetArray (Mrna, _Splicing, vnits, 8))
                  {  /* register support in surviving mrnas */
                    for (j = dx = 0 ; j < arrayMax (vnits) ; j += 8)
                      {
                        vv = arrp (vnits, j, BSunit) ;
                        hh = arrayp (mSpl, is++, HIT) ;
                        if (j == 0) dx = vv[0].i ;
                        else if (dx > vv[0].i) dx = vv[0].i ;
                        hh->x1 = uu[1].i + vv[0].i - 1 ;
                        hh->x2 = uu[1].i + vv[1].i - 1 ;
                        hh->gene =  vv[4].k ;
                        dictAdd (dict, vv[5].s, &jj) ;
                        hh->a1 = jj ;
                        hh->a2 = vv[7].k ;
                      }
                    if (dx != 1) /* shift the splicing inside the mrna */
                      {
                        int m1, m2 ;
                        OBJ Cosmid = bsCreate (relevantCosmid) ;
                        if (Cosmid &&
                            bsFindKey (Cosmid, str2tag("mrnas"), uu[0].k) &&
                            bsGetData (Cosmid, _bsRight, _Int, &m1) &&
                            bsGetData (Cosmid, _bsRight, _Int, &m2))
                          {
                            if (isGeneDown) { m1 += dx - 1 ; m2 += dx - 1; }
                            else  { m1 -= dx - 1 ; m2 -= dx - 1; }
                          }
                        bsDestroy (Cosmid) ;
        
                        vtxtPrintf (txt, "Sequence \"%s\"\nmRNAs \"%s\" %d %d\n\n"
                                    , name (relevantCosmid), name(uu[0].k), m1, m2
                                    ) ;

                        uu[1].i += dx - 1 ;
                        uu[2].i += dx - 1 ;
                    
                        vtxtPrintf (txt, "mRNA \"%s\"\n-D Splicing\n-D Coding\n", name(uu[0].k)) ;
                        for (j = 0 ; j < arrayMax (vnits) ; j += 8)
                          {
                            vv = arrp (vnits, j, BSunit) ;
                            vtxtPrintf (txt, "Splicing %d %d %d %d %s %s Length %d bp\n"
                                        , vv[0].i - dx + 1, vv[1].i - dx + 1
                                        , vv[2].i, vv[3].i
                                        , name (vv[4].k), vv[5].s, vv[7].i) ;
                          }
                        if (bsGetArray (Mrna, _Coding, vnits, 5))
                          for (j = 0 ; j < arrayMax (vnits) ; j += 5)
                            {
                              vv = arrp (vnits, j, BSunit) ;
                              vtxtPrintf (txt, "Coding %d %d %d %d"
                                          , vv[0].i - dx + 1, vv[1].i - dx + 1
                                          , vv[2].i, vv[3].i
                                          ) ;
                              if (vv[4].s && *vv[4].s)
                                vtxtPrintf (txt, " \"%s\"",  vv[4].s) ;
                              vtxtPrint (txt, "\n") ;
                            }

                        vtxtPrint (txt, "\n") ;
                      }
                  }
                bsDestroy (Mrna) ;
              }
          }
      }
    arraySort (mSpl, mSplOrder) ;
    for (i = jj = 0 ; i < arrayMax (mSpl) ; i++)
      {
        hh = arrp (mSpl, i, HIT)  ;
        if (i == 0) { gx1 = hh->x1 ; gx2 = hh->x2 ; }
        if (gx1 > hh->x1) gx1 = hh->x1 ;
        if (gx2 < hh->x2) gx2 = hh->x2 ;
      }

    /* from now on we are only interested by the offset gx1 - 1 */
    gx1-- ; gx2-- ;

    vtxtPrintf (txt, "Transcribed_gene \"%s\"\n", name(gene)) ;
    /* export splicing in of gene */
    for (i = 0 ; i < arrayMax (units) ; i += 3)
      {
        uu = arrp (units, i, BSunit) ;
        vtxtPrintf (txt, "mRNA \"%s\" %d %d\n", name(uu[0].k), uu[1].i - gx1, uu[2].i - gx1) ;
      }
    for (i = jj = 0 ; i < arrayMax (mSpl) ; i++)
      {
        hh = arrp (mSpl, i, HIT)  ;
        vtxtPrintf (txt, "Splicing %d %d %s ", hh->x1 - gx1, hh->x2 - gx1, name (hh->gene)) ;
        if (strstr (name(hh->gene), "xon"))
          vtxtPrintf (txt, "%d ", ++jj) ;
        else
          vtxtPrint (txt, dictName(dict, hh->a1)) ;
        vtxtPrintf (txt, " Length %d bp\n", hh->a2) ;
      }
    arrayDestroy (mSpl) ;
    dictDestroy (dict) ;
  }
  bsRemoveTag (Gene, _mRNA) ;

  bsRemoveTag (Gene,  _Splicing) ;

  bsSave (Gene) ;

 /*** Comid info: before we bsUpdate gene ***/
  
  { 
    OBJ Cosmid = bsUpdate(relevantCosmid) ;
    if (Cosmid)
      {
        bsAddKey (Cosmid, _Transcribed_gene, gene) ;
        if (isGeneDown)
          {
            x = gb1 = ga1 + gx1  ;
            bsAddData (Cosmid, _bsRight, _Int, &x) ;
            x = gb2 = ga1 + gx2  ;
            bsAddData (Cosmid, _bsRight, _Int, &x) ;
          }
        else
          {
            x = gb1 = ga1 - gx1 ;
            bsAddData (Cosmid, _bsRight, _Int, &x) ;
            x = gb2 = ga1 - gx2 ;
            bsAddData (Cosmid, _bsRight, _Int, &x) ;
          }

        bsSave (Cosmid) ; 
      }
  }
  
  /****************/

  Gene = bsUpdate (gene) ;
  if (!Gene)
    goto done ;
  
  /*** Transcribed_from ***/

  x = gx2 - gx1 + 1 ;
  bsAddData (Gene, _Covers, _Int, &x) ;
  bsAddData (Gene, _bsRight, _Text, "bp from") ;
  bsAddData (Gene, _bsRight, _Int, &gb1) ; 
  bsAddData (Gene, _bsRight, _Int, &gb2) ;
    
  {
    OBJ Cosmid = bsCreate (cosmid) ;
    if (Cosmid)
      {
        if (bsGetKey (Cosmid, _IntMap, &map) &&
            bsGetData (Cosmid, _bsRight, _Int, &c1)  &&
            bsGetData (Cosmid, _bsRight, _Int, &c2))
	  {}  ;
        bsDestroy (Cosmid) ;
      }
  }

  
  if (map && c1 + c2)
    {
      if (c1 < c2)
        { d1 = c1 + gb1 - 1 ; d2 = c1 + gb2 - 1 ; }
      else
        { d1 = c1 - gb1 + 1 ; d2 = c1 - gb2 + 1; }
      bsAddKey (Gene, _IntMap, map) ;
      bsAddData (Gene, _bsRight, _Int, &d1) ;
      bsAddData (Gene, _bsRight, _Int, &d2) ;
    }
        
  /*** intron boundaries ***/

  if (bsGetArray (Gene, _Intron_boundaries, units, 12))
    {
      /* shift all coordinates */
      for (i = 0 ; i < arrayMax (units) ; i += 12)
        {
          uu = arrp (units, i, BSunit) ;
          if (uu[0].k == _Other)
            {
              uu[1].s = strnew (uu[1].s, handle) ;
              uu[3].i -= gx1 ;
              uu[4].i -= gx1 ;
            }
          else if (uu[2].i || uu[3].i)
            { 
              uu[2].i -= gx1 ;
              uu[3].i -= gx1 ;
            }
          
        }
      bsAddArray (Gene, _Intron_boundaries, units, 12) ;
    }
  
  /*** Assembled_from ***/
 
  if (bsGetArray (Gene, _Assembled_from, units, 12))
    {
      for (i = 0 ; i < arrayMax (units) ; i += 12)
        {
          uu = arrp (units, i, BSunit) ;
          uu[0].i -= gx1 ;
          uu[1].i -= gx1 ;
        }
      bsAddArray (Gene, _Assembled_from, units, 12) ;
    }
  
  
 done:

  bsSave (Gene) ;
  messfree (handle) ;
  messfree (mark) ;
  if (txt)
    parseBuffer (vtxtPtr (txt), 0) ;
  vtxtDestroy (txt) ;
  arrayDestroy (units) ;
  arrayDestroy (vnits) ;
} /* mrnaSaveGeneKillNonBest */

/*********************************************************************/

static int mrnaGetGapSize (S2M *s2m, SC* sc, Array estHits, SMRNA *smrna, HIT *gp)
{
  int ii, delta = 0, bestDelta = 9999, eMax = arrayMax (estHits), nGapping = 0 ;
  int jj, x2, y1, pass, u1, u2 ;
  HIT *up, *vp, *w1p, *w2p ;
  KEYSET ks = 0 ; 
  KEY clone ;
  float size = 0 ;
  int gpa1 = smrna->a1 + gp->a1 - 1 ;
  int gpa2 = smrna->a1 + gp->a2 - 1 ;
  int suspect ;
  BOOL singleEst = FALSE ;

  /* first we look for an est going above the gap */ 
  ks = keySetCreate () ;
  for (suspect = 0 ; suspect < 2 ; suspect++)
    for (ii = 0, up = arrp(estHits, ii, HIT) ; ! singleEst && ii < eMax - 1 ; up++, ii++)
      {
        if (!class(up->est)) /* ghosts */
          continue ; 
        if (up->type & gDroppedGene)
          continue ;
        if (suspect == 0 &&
            keyFindTag (up->cDNA_clone, _Suspected_internal_deletion) &&
            !keyFindTag (up->cDNA_clone, _Manual_no_internal_deletion))
          continue ;
        if (suspect == 1 &&
            ! (keyFindTag (up->cDNA_clone, _Suspected_internal_deletion) &&
               !keyFindTag (up->cDNA_clone, _Manual_no_internal_deletion)))
          continue ;
        if (!keySetFind (smrna->clones, up->cDNA_clone, 0))
          continue ;
        if (keySetFind (ks, up->est, 0))
          continue ;
        keySetInsert (ks, up->est) ;
        if ((up->type & gX) && up->a1 < gpa1)
          {
            w1p = up ;
            for (vp = up + 1, jj = ii + 1 ; jj < eMax ; vp++, jj++)
              if (vp->a1 < gpa1 && vp->est == up->est && (vp->type & gX))
                w1p = vp ;
            for (vp = up + 1, jj = ii + 1 ; ! (vp->type & gSuspect) && jj < eMax ; vp++, jj++)
              if (vp->a2 > gpa2 && vp->est == up->est && (vp->type & gX))
                {
                  delta = vp->x1 - w1p->x2 - 1 ;
                  if (delta*delta < bestDelta*bestDelta) bestDelta = delta ;
                  nGapping = 1 ;
                  singleEst = TRUE ;
                  break ;
                }
          }
      }
  if ( nGapping) { delta = bestDelta ; goto done ; }
  /* then we look for a clone where 5 and 3 are both aligned and size is known */
  ks = keySetReCreate (ks) ;
  /* in pass 0, we only consider sized clones, in pass 1 we default the sie to 1800 bp */
  for (pass = 0 ; !nGapping && pass < 2 ; keySetMax(ks) = 0, pass++)
    for (ii = 0, up = arrp(estHits, ii, HIT) ; ii < eMax - 1 ; up++, ii++)
      {
        clone = up->cDNA_clone ;
        if (!keySetFind (smrna->clones, clone, 0))
          continue ;
        
        if (!class(up->est)) /* ghosts */
          continue ; 
        if (! (up->type & gX))
          continue ;
        if (up->type & gDroppedGene)
          continue ;
        if (keySetFind (ks, clone, 0))
          continue ;  
        if (keyFindTag (up->cDNA_clone, _Suspected_internal_deletion))
          continue ;
        keySetInsert (ks, clone) ;
        
        size = pass ? 1800 : 0 ;
        if (keyFindTag (clone, str2tag ("PCR_product_size")))
          {
            OBJ Clone = bsCreate (clone), Library = 0 ;
            KEY library = 0 ;
            
            if (Clone)
              {
                bsGetData (Clone, str2tag ("PCR_product_size"), _Float, &size) ;
                size *= 1000 ;
                bsGetKey (Clone, _Library, &library) ;
                bsDestroy (Clone) ;
                if (library &&
                    (Library = bsCreate (library)))
                  {
                    int vl = 0 ;
                    if (bsGetData (Library, str2tag("Vector_length") , _Int, &vl) &&
                        size > vl)
                      size -= vl ; /* yk library */
                    bsDestroy (Library) ;
                  }
              }
          }
        if (!size)
          continue ;
        x2 = 0 ; w2p = 0;
        /* look for last hit of 5p read */
        for (jj = 0 ; jj < arrayMax(estHits) ; jj++)
          {
            w1p = arrp (estHits, jj, HIT) ; 
            if (! (w1p->type & gX))
              continue ;
            if (w1p->cDNA_clone != clone ||
                w1p->x1 > w1p->x2 ||
                w1p->a1 > gpa1)
              continue ;
            x2 = w1p->x2 ; w2p = w1p ;
          }
        if (!x2)
          continue ;
        
        
        y1 = 0 ; w1p = w2p ;
        /* look for first hit of 3p read */
        for (jj = 0 ; jj < arrayMax(estHits) ; jj++)
          {
            w2p = arrp (estHits, jj, HIT) ; 
            if (! (w2p->type & gX))
              continue ;
            if (w2p->cDNA_clone != clone ||
                w2p->x1 < w2p->x2 ||
                w2p->a2 < gpa2)
              continue ;
            y1 = w2p->x1 ;
            break ;
          }
        if (!y1)
          continue ;
        /* deduce from that the known exons between me and the start of the gap */

        for (jj = 0, vp = arrp (smrna->hits, 0, HIT) ; jj < arrayMax(smrna->hits) ; jj++, vp++)
          {
            u1 = w1p->a2; u2 = gpa1 ; /* questionable area */
            if (!(vp->type & gX)) continue ;
            if (vp->a2 < u2) u2 = vp->a2 ;
            if (vp->a1 > u1) u1 = vp->a1 ;
            if (u2 > u1) size -= u2 - u1 ;

            u1 = gpa2 ; u2 = w2p->a1 ; /* questionable area */
            if (!(vp->type & gX)) continue ;
            if (vp->a2 < u2) u2 = vp->a2 ;
            if (vp->a1 > u1) u1 = vp->a1 ;
            if (u2 > u1) size -= u2 - u1 ;
          }
        nGapping++ ;
        delta += size - x2 - y1 ;
      }

 done:
  if (nGapping)
    delta /= nGapping ;
  else
    delta = 300 ;

  if (delta < 0) delta = - delta ;
  if (delta > gpa2 - gpa1)
    delta = gpa2 - gpa1 ;
  if (singleEst)
    { if (delta < 3) delta = 3 ; }
  else
    { if (delta < 9) delta = 9 ; }

  keySetDestroy (ks) ;
  return delta ;
}

/*********************************************************************/
#define NGAPMAX 3
static void mrnaMakeDnas (S2M *s2m, SC* sc, Array estHits, SMRNA *smrna)
{
  int i, ii, n, iIntron, lMax, offset, iDna = 0, ngap, loop, loopMax, gaps[12] ;
  Array dna = 0, longDna, dnas = 0, hits = 0, introns = 0 ;
  char nextBase = 0 ; /* a hack to see the stop in vidal clones */
  char nextBase2 = 0 ; /* a hack to see the stop in vidal clones */
  HIT *up ;
  ORFS *orf ;
  BOOL isSL = FALSE ;

  arrayDestroy (smrna->dnas) ;
  if (!arrayMax(smrna->hits))
    return ;
  lMax = arrayMax(s2m->dnaD) ;
  if (sc->isUp)
    { longDna = s2m->dnaR ; offset = lMax - sc->a1 + smrna->a1 - 1 ; }
  else
    { longDna = s2m->dnaD  ; offset = sc->a1 + smrna->a1  - 2 ; }
  
  dnas = arrayHandleCreate (12, ORFS, s2m->h) ;

  for (ii = 0, up = arrp (smrna->hits, 0, HIT) ; ii < arrayMax(smrna->hits);  up++, ii++)
    if (up->type & gS) { isSL = TRUE ; }
    
  for (ii = 0, ngap = 0, up = arrp (smrna->hits, 0, HIT) ; ii < arrayMax (smrna->hits) ; up++, ii++)
    if (up->type & gGap)
      { 
        i = mrnaGetGapSize(s2m, sc, estHits, smrna, up) ;
        if (i < 9)
          { up->type &= !gGap ; up->type |= (gI | gJ) ; }
      }

  loopMax = 1 ;
  for (ii = 0, ngap = 0, up = arrp (smrna->hits, 0, HIT) ; ngap < NGAPMAX && ii < arrayMax (smrna->hits) ; up++, ii++)
    if (up->type & gGap)
      { gaps[ngap++] = mrnaGetGapSize(s2m, sc, estHits, smrna, up) ; loopMax *= 4 ; }

  for (loop = 0 ; loop < loopMax && loop < (4 << NGAPMAX) ; loop++)
    {
      n = 0 ; orf = arrayp (dnas, iDna++, ORFS) ; ngap = 0 ;
      dna = orf->dna = arrayHandleCreate (20000, char, s2m->h) ;
      hits = orf->hits = arrayHandleCopy (smrna->hits, s2m->h) ;
      introns = orf->introns = arrayHandleCreate (120, int, s2m->h) ;
      nextBase = nextBase2 = 0 ;

      for (ii = iIntron = 0, up = arrp (hits, 0, HIT) ; ii < arrayMax (hits) ; up++, ii++)
        if (up->type & gX)
          {
            if (!n && !isSL) 
              {
                int j ;

                orf->topExon = ii ; 
                if (getPleaseStealUpStream ())
                  j = 18  < sc->a1 - 1 ? 18 : sc->a1 - 1 ;
                else
                  j = 0 ;
                j = j/3 ; j = 3 * j ; /* steal in codons */
                i = offset + up->a1 - 1 - j ;
                while (i < 0 && j > 0) { i += 3 ; j -= 3 ; }
                orf->stolen = j ;
                for ( ;  j > 0 ; i++ , j--)
                  array (dna, n++, char) = arr (longDna, i, char) ;
              }
            up->x1 = n + 1 ;
            array (dna, n +  up->a2 - up->a1 + 1, char) = 0 ;
            for (i = offset + up->a1 - 1 ; i >= 0 && i < lMax && i <=  offset + up->a2 - 1 ; i++)
              arr (dna, n++, char) = arr (longDna, i, char) ;
            nextBase = i < arrayMax (longDna) ? arr (longDna, i, char) : 0 ;
            nextBase2 = i < arrayMax (longDna) - 1 ? arr (longDna, i + 1, char) : 0 ;
            up->x2 = n ;
          }
        else if (up->type & gI)
          {
            array (introns, iIntron++, int) = n ;
          }
        else if (up->type & gGap)
          {
            int delta = loop ;
            int gap = ngap < NGAPMAX ? gaps[ngap] : 300 ;

            /* delta adds 0, 1, 2 phase base as a function of the loop */   
            for (i = 0 ; i < ngap ; i++) delta /= 4 ;
            delta = delta % 4 ;
         
            if (delta == 0 && ngap < NGAPMAX) /* try open gap */
              { 
                up->x1 = n + 1 ;
                array (dna, n +  up->a2 - up->a1 + 1, char) = 0 ;
                for (i = offset + up->a1 - 1 ; i >= 0 && i < lMax && i <= offset + up->a2 - 1 ; i++)
                  arr (dna, n++, char) = arr (longDna, i, char) ;
                up->type |= gOpenGap ; 
                up->x2 = n ;
                up->zone = ngap ;
              }
            else
              { 
                array (dna, n +  gap + delta + 1, char) = 0 ;
                for (i = 0 ; i < gap + delta - 1 ; i++) /* delta - 1 to use 0 for open,  1 2 3 for n, to favor choice of 0 later*/
                  arr (dna, n++, char) = N_ ;
              }
            ngap++ ;
          }
      array (dna, n++, char) = nextBase ; /* a secret hiding for next base */
      array (dna, n++, char) = nextBase2 ; /* a secret hiding for next base */
      array (dna, n++, char) = 0 ; /* zero terminate again */

      arrayMax(dna) -= 3 ; /* hide the stolen bases */
    }
  smrna->dnas = dnas ; /* for now those are just the dna contigs betwee gaps */
  arrayDestroy (smrna->orfs) ; 
}

/*********************************************************************/
   
static int orfMetOrder (const void *a, const void *b)
{
  const ORFS *oa = (const ORFS *)a, *ob = (const ORFS *)b ;
  int n ;
  
  n = oa->iDna - ob->iDna ;
  if (n) return n ;

  n = oa->met - ob->met ;
  if (n) return n ;
  return 0 ;  /* increasing met position */
}
/*********************************************************************/
   
static int orfOrder (const void *a, const void *b)
{
  const ORFS *oa = (const ORFS *)a, *ob = (const ORFS *)b ;
  int aopen = oa->cds - oa->nX ; /* jan 2004, was before -oa->nX/2 but this favors huge gaps */
  int bopen = ob->cds - ob->nX ;
  int aIntron = oa->nIntron - (oa->nX ? 1 : 0) ; /* count a gap as a negative intron */
  int bIntron = ob->nIntron - (ob->nX ? 1 : 0) ; /* count a gap as a negative intron */

  /* prefer longer and upstream */
  if (oa->max < ob->max && aopen > bopen) return -1 ;
  if (ob->max < oa->max && bopen > aopen) return 1 ;
  /* prefer much longer and downstream */
  if (aopen > 0 && 2 * aopen > 3 * bopen + 90) return -1 ;
  if (bopen > 0 && 2 * bopen > 3 * aopen + 90) return 1 ;
  /* quasi equal length, prefer introns */
  if (aIntron > 0 && bIntron <= 0) return -1 ;
  if (bIntron > 0 && aIntron <= 0) return  1 ;
  /* prefer longer */
  return bopen - aopen ;  /* decreasing cds order */
}
/*
orf 0:: iDna= 0  frame=1  topExon= 0  met= -2 upStop= -2  downStop=103  nIntron= 0 nOpen=102  cds=102 metOpen=102 nX=  0
orf 1:: iDna= 0  frame=2  topExon= 0  met= -1 upStop= -1  downStop= 56  nIntron= 0 nOpen= 54  cds= 54 metOpen= 54 nX=  0
orf 2:: iDna= 0  frame=3  topExon= 0  met= 39 upStop=  0  downStop= 87  nIntron= 0 nOpen= 84  cds= 84 metOpen= 48 nX=  0
orf 3:: iDna= 0  frame=3  topExon= 0  met=150 upStop= 96  downStop=198  nIntron= 1 nOpen= 99  cds= 48 metOpen= 48 nX=  0
*/
/*********************************************************************/

static int orfOrderByUpStop (const void *a, const void *b)
{
  const ORFS *oa = (const ORFS *)a, *ob = (const ORFS *)b ;
  return oa->upStop - ob->upStop ;
}

/*********************************************************************/

static int orfOrderByCds (const void *a, const void *b)
{
  const ORFS *oa = (const ORFS *)a, *ob = (const ORFS *)b ;
  
  int n ;

  n= oa->cds - ob->cds ;
  if (n) return -n ; /* long first */

  n = oa->upStop - ob->upStop ;
  return n ; /* topmost first */
}

/*********************************************************************/

static void mrnaLocateOrfs (S2M *s2m, int maxSmrnas, SMRNA *smrna, BOOL isMrna5pComplete)
{
  Array hits = 0, dna = 0, dnas = smrna->dnas, orfs = 0, orfs2 = 0, introns = 0 ;
  int i, ii, iBase, iDna, j, jj, met, leu, lastStop, start, nOpen, cds, stolen,
    bestCds = 0, bestDna, frame, dMax = arrayMax(dnas), oMax = 0, dnaMax ;
  ORFS *contig, *orf, *orf2 ;
  char tt, *cp ;
  int isOpenGap = 0, isGap = 0, lengthOpenGap = 0, nX = 0, nX1 = 0, firstXsinceStop = 0 ;
  HIT *up, *vp ;
  BOOL isPolyA = FALSE, isReal3p = FALSE, isSL = FALSE ;
  BOOL  singleProduct = mrnaPleaseSingleORF () ; /** SINGLE/MULTIPLE ORF **/
  BOOL ignorePolyA = getPleaseIgnorePolyAInProduct() ;
  BOOL useLeu = mrnaPleaseUseLeucine () ;
  BOOL isInFirstGap, upGap ;
  static int nLocateCall = 0 ; /* count call to this function, debugging purpose */
  char *translationTable = pepGetTranslationTable (s2m->cosmid1 ? s2m->cosmid1 : s2m->cosmid, 0) ;

  if (arrayMax(smrna->hits))
    {
      for (j = 0, up = arrp (smrna->hits, 0, HIT) ; j < arrayMax(smrna->hits) ;  up++, j++)
        {
          if (j < 1)
            {
              if (up->type & (gS | gCompleteCDS))
                { isSL = TRUE ; /* isSL1 = TRUE ; */ }
              if (up->type & (gS0))
                { isSL = TRUE ; /* isSL0 = TRUE ; */ }
            }
          if (j ==  arrayMax(smrna->hits) - 1)
            {
              if (!ignorePolyA &&  (up->type & gA))
                isPolyA = TRUE ;
              if (up->type & gReal3p)
                isReal3p = TRUE ; 
            }
        }
    }

  /* for each dna find the stops and Met of each piece in each frame */ 

  nLocateCall++ ;
  if (0) printf(" mrnaLocateOrfs nLocateCall=%d\n",nLocateCall) ;
  if (0 && nLocateCall == 361)
    invokeDebugger() ;
  orfs = arrayCreate (300 * dMax, ORFS) ; oMax = 0 ;

  for (iDna = 0 ; iDna < dMax ; iDna++)
    {
      bestCds = 0 ;
      contig = arrp (dnas, iDna, ORFS) ; 
      introns = contig->introns ;
      isOpenGap  = isGap = lengthOpenGap = 0 ;
      hits = contig->hits ;
      if (arrayMax(hits))
        for (j = 0, up = arrp (hits, 0, HIT) ; j < arrayMax(hits);  up++, j++)
          if (up->type & gOpenGap) { isOpenGap = up->x1 ; lengthOpenGap = up->x2 - up->x1 + 1 ; break ; }
      if (arrayMax(hits))
        for (j = 0, up = arrp (hits, 0, HIT) ; j < arrayMax(hits);  up++, j++)
          if (up->type & gGap) { isGap = up->x1 ; break ; }
      dna = contig->dna ;
      stolen = contig->stolen ;
      dnaMax = arrayMax(dna) ;

      for (frame = 0 ; frame < 3 ; frame++)
        {
          orf = 0 ; nX = nX1 = 0 ; firstXsinceStop = dnaMax ;
          lastStop = met = leu = -3 + frame ;
          isInFirstGap = upGap = FALSE ;

          for (iBase = frame, cp = arrp (dna, iBase, char) ; iBase < dnaMax - 2  ; iBase += 3, cp += 3)
            {
              tt = e_codon (cp, translationTable) ;
              if (tt == '*')
                {
                  if (met < 0 && lastStop >= 0 && lastStop < stolen)
                    lastStop = met ; /* ignore stolen stop */
                  nOpen = iBase - lastStop - 3 ; 
                  /* discard ORFSs inside ORF-gaps */
                  if (isOpenGap)
                    {
                      if (lastStop >= isOpenGap)
                        goto nextFrame ;
                      if (iBase < isOpenGap)
                        nOpen = 0 ;  /* gap is longer than orf on any side, reject */
                      else if (arrayMax(hits))
                        for (j = 0, up = arrp (hits, 0, HIT) ; j < arrayMax(hits); up++, j++)
                          {
                            if (up->type & gGap)
                              {  /* test on XJ173==yk1495d2, him-4 */
                                if (up->zone < NGAPMAX && 
                                    (up->type & gOpenGap) &&
                                    (up->x1 <= lastStop || up->x2 >= iBase)
                                    )
                                  { nOpen = 0 ; break ; } /* gap is longer than orf on any side, reject */
                              }
                          }
                    }
                  cds = 0 ;
                  if (nOpen)
                    { /* cds starts above the met if the gene seems open from above */
                      start = lastStop + 3 ;
                      if (lastStop < 0  && !(isSL || (isMrna5pComplete && met > MINI_LEU2MET)))
                        { cds = nOpen; start = lastStop + 3 ; }
                      else  if (useLeu && leu >= 0 && 
                                (leu + MINI_LEU2MET < met || (iBase - leu > MINI_LEU2END && met < 0))
                                )
                        { start = leu ; cds = iBase - leu ; }
                      else if (met >= 0)
                        { start = met ; cds = iBase - met ; }
                      upGap = FALSE ;
                      if (nX && start > firstXsinceStop)
                        { 
                          start = firstXsinceStop + nX1 - 9 ; lastStop = met = leu = -3 + frame ;
                          nOpen = cds = iBase - start ;
                          nX = nX - nX1 + 9 ;
                          upGap = TRUE ;
                        }
                    }
                  if (cds > 3 && (cds > bestCds || cds > 3))
                    {
                      orf = arrayp (orfs, oMax++, ORFS) ;
                      orf->iDna = iDna ;
                      orf->nOpen = nOpen ;
                      orf->cds = cds;
                      orf->frame = frame + 1 ;
                      orf->start = start + 1 ;     /* if <= 0 no stop, but tells me the frame */
                      orf->met = met + 1 ;         /* if <= 0 no met, but tells me the frame */
                      orf->leu = leu + 1 ;         /* if <= 0 no met, but tells me the frame */
                      orf->upStop = lastStop + 1 ; /* if <= 0 no stop, but tells me the frame */
                      orf->downStop = iBase + 1 ;
                      orf->max = iBase ; /* last base of last codon bio-coord, does not include the stop */
                      orf->nX = nX ; nX = 0 ;
                      orf->upGap = upGap ;
                      orf->stolen = stolen ;
                      orf->isSL = isSL ;
                      if (cds > bestCds) {  bestCds = cds ; }
                    }
                  lastStop = iBase ; leu = met = -3 + frame ;  firstXsinceStop = dnaMax ;
                }
              else if (tt == 'M' && (met < 0 || (met < stolen && (isSL || iBase < stolen + 20))))
                met = iBase  ;
              else if (iBase > stolen && useLeu && /* 2006/04/21 complex rule from Kosak 2002 */
                       *(cp+2) == G_ && 
                       ((*(cp-3) == G_ && *(cp+3) == G_) || (*(cp-3) == A_)) &&
                       (*(cp) == A_ || *(cp+1) == T_) && 
                       leu < 0)
                leu = iBase  ;
              if (tt == 'X')
                {
                  nX += 3 ; 
                  if (nX1 ||
                      (
                       e_codon (cp+3, translationTable) == 'X' && 
                       e_codon (cp+6, translationTable) == 'X' && 
                       e_codon (cp+9, translationTable) == 'X' && 
                       e_codon (cp+12, translationTable) == 'X')
                      )
                    {                       
                      if (iBase < firstXsinceStop) 
                        { firstXsinceStop = iBase ; isInFirstGap = TRUE ; }
                      if (isInFirstGap)
                        nX1 += 3 ;
                    }
                }
              else if (isInFirstGap)
                isInFirstGap = FALSE ;
            }
          /* register final ORF (open to the right */
          if (met < 0 && lastStop >=0 && lastStop < stolen)
            lastStop = met ; /* ignore stolen stop */
          nOpen = iBase - lastStop - 3 ;
         /* discard ORFs inside ORF-gaps */
          if (isOpenGap && arrayMax(hits))
            {
              for (j = 0, up = arrp (hits, 0, HIT) ; j < arrayMax(hits); up++, j++)
                {
                  if (up->type & gGap)
                    {
                      if (up->zone < NGAPMAX && 
                          (up->type & gOpenGap) &&
                          (up->x1 <= lastStop || up->x2 >= iBase)
                          )
                        { nOpen = 0 ; break ; }
                    }
                }
            }
          if (( isReal3p || isPolyA) && nOpen < 600) /* mieg 2008_09_14, mais prevent later detection of supect-polyA */
	    nOpen = 0 ;
          cds = 0 ;
          if (nOpen)
            { /* cds starts above the met if the gene seems open from above */
              start = lastStop + 3 ;
              if (lastStop < 0  && !(isSL || (isMrna5pComplete && met > MINI_LEU2MET)))
                { cds = nOpen; start = lastStop + 3 ; }
              else  if (useLeu && leu >= 0 && (leu + MINI_LEU2MET < met || (iBase - leu > MINI_LEU2END && met < 0)))
                { start = leu ; cds = iBase - leu ; }
              else if (met >= 0)
                { start = met ; cds = iBase - met ; }
              upGap = FALSE ;
              if (nX && start > firstXsinceStop)
                { 
                  start = firstXsinceStop + nX1 - 9 ; lastStop = met = leu = -3 + frame ;
                  nOpen = cds = iBase - start ;
                  nX = nX - nX1 + 9 ;
                  upGap = TRUE ;
                }
            }
          if (cds > 3  && (cds > bestCds || cds > 3))
            {
              orf = arrayp (orfs, oMax++, ORFS) ;
              orf->iDna = iDna ; 
              orf->nOpen = nOpen ;
              orf->cds = cds;
              orf->frame = frame + 1 ;
              orf->met = met + 1 ;                /* if <= 0 no met, but tells me the frame */
              orf->leu = leu + 1 ;                /* if <= 0 no met, but tells me the frame */
              orf->start = start + 1 ;            /* if <= 0 no stop, but tells me the frame */
              orf->upStop = lastStop + 1 ;        /* if <= 0 no stop, but tells me the frame */
              orf->downStop = - 2 + frame ;       /* if <= 0 no stop, but tells me the frame */
              orf->max = iBase  ;     /* last base of last codon bio-coord, does not include the stop */ 
              orf->nX = nX ; nX = 0 ;
              orf->upGap = upGap ;
              orf->stolen = stolen ;
              orf->isSL = isSL ;
              if (cds > bestCds)  {  bestCds = cds ; }
            }
        nextFrame:
          continue ;
        }
    }
  /* do not steal a stop unless the met is also stolen */
  for (ii = 0 ; ii < arrayMax(orfs) ; ii++)
    {
      int dx, x;
      orf = arrp (orfs, ii, ORFS) ;
      if (orf->upStop > 0 && orf->upStop < orf->stolen && (orf->met <= 0 || orf->met > orf->stolen))
        {
          x = orf->upStop + 3 ;
          while (x <= orf->stolen) x += 3 ;
          dx = x - orf->start ;
          orf->start = x ; orf->cds -= dx ; orf->nOpen = orf->cds ;
        }
    }        
  
  /* count introns in each orf */
  for (ii = 0 ; ii < arrayMax(orfs) ; ii++)
    {
      int i, n = 0, x;
      orf = arrp (orfs, ii, ORFS) ;
      contig = arrp (dnas, orf->iDna, ORFS) ; 
      if ((introns = contig->introns))
        for (i = n = 0 ; i < arrayMax(introns) ; i++)
          {
            x = arr (introns, i, int) ;
            if (x > orf->start && x < orf->max) 
              n++ ;
          }
      if (orf->cds > bestCds)  {  bestCds = orf->cds ; } 
      orf->nIntron = n ;
    }

  /* select the uORF in each dna */
  arraySort (orfs, orfMetOrder)  ;
  for (iDna = 0 ; iDna < dMax ; iDna++)
    {
      for (ii = 0 ; ii < arrayMax(orfs) ; ii++)
	{
	  orf = arrp (orfs, ii, ORFS) ;
	  if (orf->iDna > iDna) break ;
	  if (orf->iDna < iDna) continue ;
	  if (orf->met > 0 && orf->met > orf->stolen && orf->cds > 15)
	    orf->isUorf = 1 ; 
	  if (orf->met > 0 && orf->met > orf->stolen) /* only the first ATG can be a uORF */
	    break ;
	}
    }
  /* sort */
  arraySort (orfs, orfOrder)  ;
  /* select a reasonable set of orfs for each dna */
  if (bestCds < 18) bestCds = 18 ; 
  orfs2 = arrayCreate (arrayMax(orfs) + 1, ORFS) ; /* avoid eventual zero */
  for (ii = jj = 0 ; ii < arrayMax(orfs) ; ii++)
    {
      orf = arrp (orfs, ii, ORFS) ;
      if (orf->cds == bestCds || orf->isUorf ||
          (orf->cds + (orf->nIntron > 0 ? 60 : 0) + (isSL ? 50 : 0) > 150))
        {
          orf2 = arrayp (orfs2, jj++, ORFS) ;
          *orf2 = *orf ;
        }
    } 
  
  arrayDestroy (orfs) ;
  orfs = orfs2 ; 
  
  /* select best dna */
  bestDna = 0 ; 
  if (arrayMax(orfs) > 0 && dMax > 1)
    {
      int bestDnaOpen = 0 ;

      for (ii = 1 ; ii < dMax ; ii <<= 2)
        bestDna += ii ; /* so that all gaps by default are in frame 1, not zero */
      for (ii = 0, orf = arrp (orfs, ii, ORFS) ; ii < arrayMax(orfs) ; orf++, ii++)
        if (orf->cds - orf->nX > bestDnaOpen)
          {  bestDna = orf->iDna ;  bestDnaOpen = orf->cds - orf->nX ; }
            
      /* keep the orfs of the best dna */
      orfs2 = arrayHandleCreate (arrayMax(orfs), ORFS, s2m->h) ;
      {
        for (ii = jj = 0 ; ii < arrayMax(orfs) ; ii++)
          {
            orf = arrp (orfs, ii, ORFS) ;
            if (orf->iDna == bestDna) 
              { /* register */
                orf2 = arrayp (orfs2, jj++, ORFS) ;
                *orf2 = *orf ;
              }
          }
      }
      arrayDestroy (orfs) ;
      orfs = orfs2 ; 
    }

  /* keep the best */
  orfs2 = smrna->orfs = arrayHandleCreate (arrayMax(orfs) + 1, ORFS, s2m->h) ;
  {
    int bestStart = 0, bestCds = 0 ;
    for (ii = jj = 0 ; jj < 25 && ii < arrayMax(orfs) ; ii++)
      {
        BOOL ignore = FALSE ;
        
        orf = arrp (orfs, ii, ORFS) ;            
        /* keep a single ORF per mrna  ? */
        /* reject if < 600 bp */
        if (jj && singleProduct && orf->cds < 600)
          continue ;
        if (jj && orf->leu >= 0 && orf->met < 0 && orf->cds < 600)
	  continue ;
        /* reject gaps */
        if (orf->nX > 0 && orf->cds - orf->nX < 300)
          continue ;
	/* the best or any orf not starting 9bp above cannot be considered to be an uOrf */
        if (!jj) { bestStart = orf->start ; bestCds = orf->cds ; orf->isUorf = 0 ; }
	else if (bestStart <= orf->start + 9) orf->isUorf = 0 ; 
	if (orf->isUorf && orf->start < orf->met)
	  {
	    if (orf->cds < bestCds/4 && orf->cds < 150 && 
		orf->cds > 15 + orf->met - orf->start && ! (orf->start == orf->leu)) 
	      { orf->cds -= orf->met - orf->start ; orf->start = orf->met ; }
	    else
	      orf->isUorf = 0 ;
	  }
        /* compare to previous accepted, avoiding overlaps */
        if (! orf->isUorf)
	  for (i = 0 ; !ignore && i < arrayMax(orfs2) ; i++)
	    {
	      int a1, a2, b1, b2, c1, c2 ;
	      
	      orf2 = arrayp (orfs2, i, ORFS) ;
	      a1 =  orf->start ; a2 = orf->max ; /* mieg 2006_10_02, was nOpen */
	      b1 =  orf2->start ; b2 = orf2->max ;
	      c1 = a1 > b1 ? a1 : b1 ; c2 = a2 < b2 ? a2 : b2 ;
	      
	      if (orf2->nIntron >= orf->nIntron &&  /* reject if large intersect */
		  orf->cds - orf->nX < 300 &&
		  (
		   (c2 - c1) * 10 > orf->cds * 8 ||
		   ((c2 - c1) * 10 > orf->cds * 5 && (b2 - b1) > (a2 - a1) * 3)
		   ) &&
		  (c2 - c1) * 10 < orf2->cds * 5)
		ignore = TRUE ;
	      else if (orf->cds - orf->nX + (orf->nIntron > 0 ? 60 : 0) < 150 &&
		       orf->met > orf2->met)
		ignore = TRUE ; /* throw away downstreams who unduly benefited from  isSl */
	    }
        
        if (!ignore) /* register */
          {
            orf2 = arrayp (orfs2, jj++, ORFS) ;
            *orf2 = *orf ;
            if (jj == 1) orf2->isBest = TRUE ;
            /* if (!orf->nX) hasNoGap = TRUE ; */
          }
      }
  }
  
  smrna->bestDna = bestDna ;
  if (dMax == 4 &&  bestDna == 0 && orfs2 && arrayMax (orfs2) &&
      maxSmrnas > 1)
    {  /* 2007_02_07
        * we actually use an orf gap 
        * promote the mrna gap to be an exon 
        * doing this, we allow fusion to other partial exons
        * but we lose the orf-gap annotation
        */
      orf = arrayp (orfs2, 0, ORFS) ;
      contig = arrp (dnas, bestDna, ORFS) ; 
      hits = contig->hits ;
      for (up = 0, i = 0 ; i < arrayMax (hits) ; i++)
        {
          up = arrp (hits, i, HIT) ;
          if (
              (up->type & gOpenGap) &&
              orf->met < up->x1 &&
              orf->downStop > up->x2)
            break ;
          else
            up = 0 ;
        }
      hits = smrna->hits ;
      for (i = 1 ; up && i < arrayMax (hits) - 1 ; i++)
        {
          vp = arrp (hits, i, HIT) ;
          if ((vp->type & gGap) &&
              vp->a1 == up->a1 &&
              vp->a2 == up->a2 &&
              (vp-1)->type & gX &&
              (vp+1)->type & gX)
            {
              (vp-1)->type |= (vp+1)->type ;
              (vp-1)->x2 = (vp+1)->x2 ;
              (vp-1)->a2 = (vp+1)->a2 ;
              for (j = i ; j < arrayMax (hits) - 2 ; j++)
                {
                  vp = arrp (hits, j, HIT) ;
                  *vp = *(vp+2) ;
                }
              arrayMax (hits) -= 2 ;
              up = 0 ; /* to break the loop */
            }
        }
    }
  arrayDestroy (orfs) ;
  /* showOrfs (smrna->orfs) ; */
} /* mrnaLocateOrfs */
  
/*********************************************************************/

static BOOL smrnaDnaDestroy (SMRNA *smrna)
{
  if (smrna && smrna->dnas)
    {
      if (! arrayExists (smrna->dnas))
        messerror ("Bad call to mrnaDnasDestroy : dnas") ;
      else
	arrayDestroy (smrna->dnas) ;
    }
  if (smrna && smrna->orfs)
    {
      if (! arrayExists (smrna->orfs))
        messerror ("Bad call to mrnaDnasDestroy : orfs") ;
      else
	arrayDestroy (smrna->orfs) ;
    }
  return TRUE ;
} /* smrnaDnaDestroy */

/*********************************************************************/

static BOOL mrnaConstructOrfs  (S2M *s2m, SC *sc, SMRNA *gmrna, Array smrnas)
{
  int ii ;
  SMRNA *smrna ;

  for (ii = 0 ; ii < arrayMax(smrnas) ; ii++)
    {
      smrna = arrp (smrnas, ii, SMRNA) ;
      if (smrna->composite)
	continue ;
      if (smrna->hits)
	{
	  /*       printf(" makeMrnaGene nCall = %d\n",nMrnaCall) ; */
	  chrono ("makeMrnaGeneLocate") ;
	  if (!smrna->dnas)
	    mrnaMakeDnas (s2m, sc, gmrna->estHits, smrna) ;
	  if (smrna->dnas && !smrna->orfs)
	    mrnaLocateOrfs (s2m, arrayMax (smrnas), smrna, FALSE) ;
	  /*
	    mrnaMakeOrfs (smrnas) ;
	    on a le dna en p morceaux correspondants aux exons p = #gap + 1
	    il faut choisir les orf de tout ca
	    rabouter les dna
	    prendre le ou les orf complet
	  */
	  chronoReturn () ;
	}
    }
  return TRUE ;
} /* mrnaConstructOrfs */

/*********************************************************************/

static Array mrnaRefseqMakerDna2Pep (KEY cosmid, Array dna, BOOL isMrna5pComplete)
{ 
  AC_HANDLE h = handleCreate () ;
  SMRNA *smrna ;
  S2M *s2m ; 
  ORFS *contig ;
  Array orfs = 0 ;

  smrna = (SMRNA *) halloc (sizeof (SMRNA), h) ; /* initialised to zero by halloc */
  s2m = (S2M *) halloc (sizeof (S2M), h) ; /* initialised to zero by halloc */
  
  s2m->cosmid = cosmid ; 
  s2m->h = h ;

  smrna->dnas = arrayHandleCreate (12, ORFS, s2m->h) ;
  smrna->hits = arrayHandleCreate (12, HIT, s2m->h) ;

  contig = arrayp (smrna->dnas, 0, ORFS) ; 
  contig->hits = arrayHandleCreate (12, HIT, s2m->h) ;
  contig->dna = dna ;
    
  mrnaLocateOrfs (s2m, 1, smrna, isMrna5pComplete) ;
  orfs = arrayCopy (smrna->orfs) ;
  messfree (h) ;

  return orfs ;
}

/*********************************************************************/

static KEYSET mrnaAlternative3primeLengths (S2M *s2m, SC* sc, Array estHits, Array smrnas, SMRNA *smrna)
{
  KEYSET ks = 0 ;


  ks = keySetCreate () ;
  return ks ;
}

/*********************************************************************/
#ifdef JUNK
static KEY mrnaBestAvailableClone (KEYSET goodClones, KEYSET completeClones, Array cloneHits)
{
  KEYSET clones = 0, myClones = 0 ;
  int ii, jj, loop,  bestTop = 0, j, jMax = arrayMax (cloneHits), bestError = 1000 ;
  KEY clone, bestClone = 0 ;
  HIT *wp ;

  /* always loop 1 and 2,  11 and 12 */
  for (loop = 0 ;  (!bestClone || bestTop > 200 || loop == 2 || loop == 12) && loop < 20 ; loop++)
    {
      switch (loop)
        {
        case 0: /*search among complete clones */
          myClones = completeClones ;
          continue ;
        case 10: /*search among all clones */
          myClones = goodClones ;
          continue ;
          
        case 1:
        case 11:      /* best situation search a yk clone above Met */
          clones = query (myClones, "y* && Capped_5prime") ;
          break ;
        case 2:   
        case 12:      /* best situation search a yk clone above Met */
          clones = query (myClones, "y*") ;
          break ;
        case 3:      /* use a cm  */
        case 13:
          clones = query (myClones, "IS cm*") ;
          break ;
        case 4:      /* use a vidal */
        case 14:
          clones = query (myClones, "IS mv* OR IS dv*") ;
          break ;
        default:
          continue ;
        }
      
      for (jj = (loop < 10 ? 0 : 1) ; 
           (jj == 0 || !bestClone || bestTop > 200)  && jj < 2 ; jj++)
        for (ii = 0 ; ii < keySetMax(clones) ; ii++)
          {          
            clone = keySet (clones, ii) ;
            
            /* get coords */ 
            for (j = 0, wp = arrp (cloneHits, 0, HIT) ; j < jMax ; j++, wp++)
              if (wp->cDNA_clone == clone)
                break ;
            if (j == jMax)
              continue ; 
            if (jj == 0 && !keyFindTag (clone, str2tag("Fully_sequenced")))
              continue ;
            if (!bestClone ||
                (loop > 10 && wp->a1 <  bestTop) ||     /* among full, (loop<10)  topMost not so important */
                (loop < 10 && 
                 (
                  (wp->a1 <  bestTop + 50 &&  wp->nerrAll < bestError) ||
                  (wp->a1 <  bestTop  &&  wp->nerrAll == bestError) 
                  )
                 )
                )                
              { bestClone = clone ; bestTop = wp->a1 ; bestError = wp->nerrAll ; }
          }
      keySetDestroy (clones) ;
    } 

  return bestClone ;
}  /* mrnaBestClone */
#endif
/*********************************************************************/
/*********************************************************************/

static BOOL mrnaMrnaScore (KEY mrna)
{
  int i, dna2 = 0, orf = 0, cds = 0, score = 0, x, xx, nerr, a1 = 0, a2 = 0 ;
  char bscore[3] ;
  OBJ Mrna = bsUpdate (mrna), Clone = 0 ;
  KEYSET ks = 0 ;
  KEY  _at_ac, dummy, clone ; 
  Array units = 0 ;
  BSunit *uu ;
  BOOL ok = FALSE ;
  /*  BOOL isFuzzy = FALSE, isOther = FALSE, isGap = FALSE  ; */
  if (!Mrna)
    return FALSE ;
 
  /* keep these condition in synch with choosing as des as lower down */
  if (! bsFindTag (Mrna, str2tag("Tiling_path")) ||
      ( FALSE &&
       bsGetData (Mrna, str2tag("Tiling_gap"), _Int, &xx) &&
       xx > 50
       ))
     goto abort ;

  nerr = dna2 = 0 ;
  if (bsGetKey (Mrna, _DNA, &dummy))
    bsGetData (Mrna, _bsRight, _Int, &dna2) ;
  if (!dna2)
    goto abort ;
  /*
    if (bsFindTag (Mrna, str2tag("Gap_length")))
    isGap = TRUE ;
  */
  bscore[0]='E' ;
  bscore[1] = bscore[2] = 0 ;
  if (bsGetData (Mrna, str2tag("Tiling_error"), _Int, &nerr))
    {
      if (1000 * nerr < 3 * dna2)
        bscore[0] = 'A' ;
      else if (1000 * nerr < 10 * dna2)
        bscore[0] = 'B' ;
      else if (1000 * nerr < 30 * dna2)
        bscore[0] = 'C' ;
      else if (1000 * nerr < 50 * dna2)
        bscore[0] = 'D' ;
      else
        bscore[0] = 'E' ;
    }
  else
    bscore[0] = 'A' ;

  units = arrayCreate (128, BSunit) ;
  /* get the region of active tiling */
  if (bsGetArray (Mrna, str2tag("Tiling_path"), units, 5))
    for (i = 0 ; i < arrayMax (units) ; i += 5)
      {
        uu = arrp (units, i, BSunit) ;
        if (i == 0) a1 = uu[0].i ;
        a2 = uu[1].i ;
      }
  if (!arrayMax(units))
    goto abort ;


  

  xx = 0 ;
  /* bons introns in tiled region, baddies in whole mrna  */
  if (bsGetArray (Mrna, str2tag("Splicing"), units, 6))
    for (i = 6 ; i < arrayMax (units) - 6 ; i += 6)
      {
        uu = arrp (units, i, BSunit) ;
        if (strstr(name(uu[4].k),"tron"))
          {
            if (!strcmp (uu[5].s,"gt_ag") ||
                !strcmp (uu[5].s,"gc_ag") ||
                !strcmp (uu[5].s,"at_ac")                 
                ) 
              { if (uu[-3].i > a1 &&  uu[8].i < a2) xx += 2 ; }
            else if (!strcasecmp (uu[5].s,"Fuzzy")) /* need to know number of support */
              {
                /* isFuzzy = TRUE ; */
                xx -= 3 ;
                if (i < arrayMax (units) - 17 &&
                    !strcasecmp (uu[17].s,"Fuzzy"))
                  {
                    xx -= 24 ; /* hence -30 for a double fuzzy */
                    if (getPleaseMarkDoubleFuzzy () &&
                        bsGetKey (Mrna, str2tag("Specific_clone"), &clone))
                      do
                        {
                          if ((Clone = bsUpdate (clone)))
                            {
                              bsAddTag (Clone, str2tag("Double_fuzzy")) ;
                              bsSave (Clone) ;
                            }
                        } while (bsGetKey (Mrna, _bsDown, &clone)) ;
                  }                                 
              }
          }
      }
  if (xx > 0) xx += 2 ; /* bonus */
  score += xx ;

  /* introns louches */
  _at_ac = str2tag("at_ac") ;
  if (bsGetArray (Mrna, str2tag("Intron_boundaries"), units, 5))
    {
      int olda1, olda2, ns ; /* nb de soutiens de l'intron bizare */

      for (x = ns = olda1 = olda2 = i = 0 ; i < arrayMax (units) ; i += 5)
        {
          uu = arrp (units, i, BSunit) ;
          if (uu[0].k == _gt_ag) continue ;
          if (uu[0].k == _gc_ag) continue ;
          if (uu[0].k == _at_ac) continue ;
          if (uu[0].k == _Fuzzy) continue ;
          /* if (uu[2].i < a1 ||  uu[3].i > a2) continue ; louches comptes sur tout le mrna */
          if (uu[2].i != olda1 || uu[3].i != olda2)
            { 
              /* isOther = TRUE ; */
              if (ns == 1) x -= 6 ;
              if (ns > 1) x -= 3 ;
              olda1 = uu[2].i ; olda2 = uu[3].i ; 
              ns = 1 ;
            }
          else
            ns++ ;
        }
      if (ns == 1) x -= 6 ;
      if (ns > 1) x -= 3 ;
      score += x ;
    }

  /****** blast pfam *******/
  ks = queryKey (mrna, "COUNT {>product ; Pfam} > 0") ;
  if (keySetMax (ks) > 0)
    score += 1 ; /* 3 */
  else
    {
      keySetDestroy (ks) ;
      ks = queryKey (mrna, "COUNT {>product ; Blastp} > 0") ;
      if (keySetMax (ks) > 0) 
        score += 1 ; /* 2 */
      else
        score += 0 ;   /* no pfam, no blast */
    }
  keySetDestroy (ks) ;

  /****** completness *******/
  /****** orf length ******/
  
  cds = orf = 0 ;
  bsGetData (Mrna, str2tag("Longest_ORF"), _Int, &orf) ;
  bsGetData (Mrna, str2tag("Longest_CDS"), _Int, &cds) ;
  
  /****** wack out the 3' incomplets */
  ks = queryKey (mrna, 
                 messprintf(">product ; Open_length== %d && !COOH_complete", orf)) ;
  if (keySetMax (ks) > 0) 
    bscore[0]++ ;  /* long mrna not 3' complete is not exportable */
   keySetDestroy (ks) ;

  /*****  consider 5' complets  ***/
  ks = queryKey (mrna, 
                 messprintf(">product ; Open_length== %d && NH2_complete", orf)) ;
  if (keySetMax (ks) > 0) 
    {
                     /* mrna complete on 5' end */
      if (cds >= 900) score += 4 ;
      else if (cds >= 600) score += 3 ;
      else if (cds >= 300) score += 1 ;
      else if (cds >= 270) score += 0 ;
      else score -= 2 ;
    }
  else
    {
      if (orf >= 960) score += 3 ;
      else if (orf >= 660) score += 2 ;
      else if (orf >= 360) score += 0 ;
      else if (orf >= 240) score -= 2 ;
      else score -= 3 ;
    }
  keySetDestroy (ks) ;
  

  /****** full antisens *******/
  if (bsGetArray (Mrna, str2tag("Antisens_to"), units, 2))
    {
      for (x = i = 0 ; i < arrayMax (units) ; i += 2)
        {
          uu = arrp (units, i, BSunit) ;
          x += uu[1].i ;
        }
      if (100 * x > 80 * dna2)
        score -=  3 ; 
    }

  /****** coverage ********/
  
  ks = queryKey (mrna, ">cds_covered_by ") ;
  xx = 0 ;
  if (keySetMax (ks) == 0)  /* no coverage */
    { 
      keySetDestroy (ks) ;
      ks = queryKey (mrna, ">mrna_covered_by ") ;
      if (keySetMax (ks) == 0)  /* no coverage */
        xx = 0 ;
      else if (keySetMax (ks) <= 2)
        xx = 2 ;
      else if (keySetMax (ks) == 3)
        xx = 0 ;
      else if (keySetMax (ks)== 4)
        xx = -1 ;
      else
        xx = -2 ; 
    }
  else if (keySetMax (ks) == 1)
    {
      if (score > 18)  xx = 5 ;
      else if (score > 14)  xx = 4 ;
      else if (score > 10)  xx = 3 ; 
      else if (score > 8)  xx = 2 ;
      else xx = 0 ;
    }
  else if (keySetMax (ks) == 2)
    xx = 0 ;
  else if (keySetMax (ks) == 3)
    xx = -1 ;
  else
    xx = -2 ;

   keySetDestroy (ks) ;


   score += xx ;

  /****** bonus de l'as des as des ORFs ********/
  {
    KEYSET mrnas = queryKey (mrna, ">from_gene ; > mrna ; longest_orf && tiling_path && !(tiling_error > 50)") ;
    int  bestGeneOrf = 0, myOrf = 0, dna3, nerr3, xx3 ;
    OBJ *obj ;
    
    for (i = 0 ; i < keySetMax (mrnas) ; i++)
      if ((obj = bsCreate (keySet(mrnas,i))))
        {  
          if (! bsFindTag (obj, str2tag("Tiling_path")) ||
              ( 0 &&
               bsGetData (obj, str2tag("Tiling_gap"), _Int, &xx3) &&
               xx3 > 50
               ))
            goto abort2;
          if (bsGetKey (obj, _DNA, &dummy))
            bsGetData (obj, _bsRight, _Int, &dna3) ;
          if (!dna3 ||
              ( 0 &&
               bsGetData (obj, str2tag("Tiling_error"), _Int, &nerr3) &&
               1000 * nerr3 > 3 * dna3
               ))
            goto abort2;

          if (bsGetData (obj, str2tag("Longest_orf"), _Int, &x))
            {
              if (x > bestGeneOrf)
                bestGeneOrf = x ;
              if (mrna == keySet(mrnas,i))
                myOrf = x ;
            }
        abort2:
          bsDestroy (obj) ;
        }
    if (bestGeneOrf)
      score = score * (( myOrf + bestGeneOrf - 1) /bestGeneOrf) ;
    keySetDestroy (mrnas) ;
  }

  /****** report *******/ 
  bscore[1] = 0 ; bscore[2] = 0 ;
  /*
    if (isGap) bscore[1] |= 0x1000000 ;
    if (isOther) bscore[1] |= 0x2000000 ;
    if (isFuzzy) bscore[1] |= 0x4000000 ;
  */
  bsAddData (Mrna, str2tag("Score"), _Int, &score) ;
  bsAddData (Mrna, str2tag("Rating"), _Text, bscore) ;
  ok = TRUE ;
    
  {
    KEY dnaTiledKey ;
    OBJ AM = 0 ;

    if (bsGetKey (Mrna, str2tag("RefSeqMaker"), &dnaTiledKey))
      {
        if ((AM = bsUpdate (dnaTiledKey)))
          {
            bsAddData (AM, str2tag("Score"), _Int, &score) ;
            bsAddData (AM, str2tag("Rating"), _Text, &bscore) ;
            bsAddKey (AM, _Colour, _PALEYELLOW) ;
            bsSave (AM) ;
          }
      }
  }

 abort:
  bsSave (Mrna) ;
  arrayDestroy (units) ;
  keySetDestroy (ks) ;
  return ok ;
} /* mrnaMrnaScore */

/*********************************************************************/

static BOOL mrnaFindClone2Resequence (KEY tr) 
{
  BOOL ok = FALSE ;
  OBJ Tr = bsUpdate (tr) ;
  Array units = 0 ;
  KEYSET clones = 0 ;
  BSunit *uu, *vv ;
  int ii, jj, type = 0, utr5 = 0, utr3 = 0, gap = 0, dnaLength = 0, iClone = 0 ;
  KEY clone, est, dummy ;

  if (!Tr)
    goto abort ;

  bsGetData (Tr, _Length_5prime_UTR, _Int, &utr5) ;
  bsGetData (Tr, str2tag("gap_length"), _Int, &gap) ;
  if (bsGetKey (Tr, _DNA, &dummy) &&
      bsGetData (Tr, _bsRight, _Int, &dnaLength)) ;
  else
    goto abort ;
  bsGetData (Tr, _Length_3prime_UTR, _Int, &utr3) ;

  if (utr5 > 3 ||  /* there is a stop */
      bsFindTag (Tr, str2tag("Found5p"))) type |= 1 ;
  if (!gap) type |= 2 ;
  if (utr3 > 2 ||
      bsFindTag (Tr, str2tag("Found3p"))) type |= 4 ;

  if (type == 7) /* perfect mrna */
    goto abort ;

  clones = keySetCreate () ;  
  units = arrayCreate (256, BSunit) ;
  bsGetArray (Tr, _Constructed_from, units, 5) ;

  switch (type)
    {
    case 5: /* a gap inside a perfect mrna */
      /* search all clones covering the whole CDS */
      for (ii = 0,  uu = arrp (units, ii, BSunit) ; !ok && ii < arrayMax(units) ; uu += 5, ii += 5)
        {
          est = uu[2].k ;
          if (uu[0].i < utr5 && keyFindTag (uu[2].k, _Forward))
            {
              clone = keyGetKey (est, _cDNA_clone) ;
              for (jj = ii + 1, vv = uu + 5 ;  !ok && jj < arrayMax(units) ; vv += 5, jj += 5)
                {
                  if (vv[1].i > dnaLength - utr3 &&
                      keyFindTag (vv[2].k, _Reverse) &&
                      clone == keyGetKey (vv[2].k, _cDNA_clone))
                    { keySet (clones, iClone++) = clone ; break ; }
                }              
            }
        }
      if (keySetMax(clones))
        bsAddKey (Tr, str2tag("Gap_clone"), keySet (clones, 0)) ;
      break ;
    case 1: /* a 5' gene */
      break ;
    case 4: /* a 3' gene */
      for (ii = 0,  uu = arrp (units, ii, BSunit) ; ii < arrayMax(units) ; uu += 5, ii += 5)
        {
          
        }
      
      break ;
          
    }
  ok = TRUE ;
 abort:
  arrayDestroy (units) ;
  keySetDestroy (clones) ;
  bsSave (Tr) ;
  return ok ;
} /* mrnaFindClone2Resequence */

/*********************************************************************/

static void mrnaReportGenomicError (SC* sc, HIT *tp, Array allErrors)
{
  int i, nn, pMax = arrayMax (sc->s2m->plainHits) ;
  HIT *up ;
  OBJ Cosmid = 0 ;
  KEY est = tp->est ;
  int a1 = tp->clipTop, x1 = tp->clipEnd, type = tp->type, baseShort = (tp->nerr < 255 ? tp->nerr : 255) ;

  if (type == AMBIGUE)
    return ;
  /* a1 is in  mrna coordinates */
  nn = 0 ; /* how many read go through a1 */
  for (i = 0, up = arrp (allErrors, 0, HIT) ; i < arrayMax (allErrors) ; i++, up++) 
    if (up->a1 == a1 && up->type == type && up->nerr == baseShort)
      nn++ ;
  if (nn < 3) 
    return ;

  /* go to genomic coordinates */

  for (i = 0, up = arrp (sc->s2m->plainHits, 0, HIT) ; i < pMax ; i++, up++) 
    {
      if (est != up->est)
        continue ;
      if (up->x1 <= up->x2)
        {
          if (up->x1 <= x1 && up->x2 >= x1)
            { a1 = up->a1 + x1 - up->x1 ; goto ok ; }
        }
      else
        {
          if (up->x1 >= x1 && up->x2 <= x1)
            { a1 = up->a2 - x1 + up->x2 ; goto ok ; }
        }
    }
  
  return ;

 ok:
  a1 = a1 - sc->d1 + sc->a1 ; /* move from double cosmid to single cosmid */
  Cosmid = bsUpdate (sc->s2m->cosmid) ;
  if (Cosmid)
    {
      bsAddKey (Cosmid, str2tag("Possible_genomic_error"), sc->gene) ;
      bsAddData (Cosmid, _bsRight, _Int, &a1) ;
      bsSave (Cosmid) ;
    }
}

/*********************************************************************/

static Array mrnaTiledDna (KEY mrna, Array tiling, int *offsetp, int *aliLenp, BOOL *ok)
{
  char *cp, *cq ;
  HIT *tp ;
  int i, jj, j1, offset = 0, delta = 0 ;
  Array dnaEst = 0, dnaTiled = 0 ;

  dnaTiled = arrayCreate (32000, char) ;
  *ok = TRUE ;
  for (jj = j1 = 0, tp = arrp (tiling, 0, HIT) ; jj < arrayMax(tiling) ; jj++, tp++)
    {
      if (!tp->est)
        continue ;
      dnaEst = getSeqDna (tp->est) ;
      if (!dnaEst || tp->a1 < 1 || tp->a1 > tp->a2 || tp->x1 < 1 || tp->x2 < 1 ||
          tp->x1 > arrayMax(dnaEst) || tp->x2 > arrayMax(dnaEst) ||
          (tp->a2 - tp->a1 != tp->x2 - tp->x1 && tp->a2 - tp->a1 != tp->x1 - tp->x2)
          ) 
        { if (0) messout ("tiling error in %s", name(mrna)) ; *ok = FALSE ; break ; }

      if (!j1)
        { j1 = 1 ; offset = tp->a1 - 1 ; }
      array (dnaTiled, tp->a2 - 1 - offset, char) = 0 ; /* make room */
      cp = arrp (dnaTiled, tp->a1 - 1 - offset, char)  ;
      if (!tp->reverse)
        {
          cq = arrp (dnaEst, tp->x1 - 1, char) ;
          for (i = tp->a1; i <= tp->a2 ; i++) *cp++ = *cq++ ;
        }
      else
        {
          cq = arrp (dnaEst, tp->x1 - 1, char) ;
          for (i = tp->a1; i <= tp->a2 ; i++) 
            *cp++ = complementBase(*cq--) ;
        }
    }   
  jj = arrayMax(dnaTiled) ; 
  if (jj > 0)
    for (i = 0, cp = arrp (dnaTiled, 0, char) ; i < jj ; i++, cp++)
      if (!*cp) *cp = N_ ;
  array (dnaTiled, jj, char) = 0 ; arrayMax(dnaTiled) = jj ;
  
  if (jj)
    {
      Array dna = dnaGet (mrna) ;
      int u1 = offset, u2 = arrayMax (dna), v1 = 0, v2 = jj ;
      Array err = aceDnaTrackErrors (dna, u1, &u2, dnaTiled, v1, &v2, 0, 0, 3, 5, 0, 0, FALSE) ;
      A_ERR *ep ;
      if (err && arrayMax(err))
	for (i = delta = 0, ep = arrp (err, 0, A_ERR) ; i < arrayMax (err) ; ep++, i++)
	  switch(ep->type)
	    {
	    case INSERTION: delta++ ; break ;
	    case INSERTION_DOUBLE: delta += 2 ; break ;
	    case TROU: delta-- ; break ;
	    case TROU_DOUBLE: delta -= 2 ; break ;
	    default: break ;
	    }
      arrayDestroy (dna) ;
      arrayDestroy (err) ;
    }
  *offsetp = offset ;
  *aliLenp = jj + delta ;
  return dnaTiled ;
}

/*********************************************************************/

static BOOL mrnaLocateBestProduct (OBJ Mrna, Array units, int *x1p, int *x2p) 
{
  BOOL ok = FALSE ;
  BSunit *uu ;
  int ii ;

  if (bsGetArray (Mrna, _Product, units, 3))
    for (ii = 0 ; ii < arrayMax(units) ; ii += 3)
      {
        uu = arrp (units, 0, BSunit) ;
        if (keyFindTag (uu[0].k, str2tag ("Best_product")))
          {
            /* xx = uu[2].i - uu[1].i ; */
            *x1p = uu[1].i ; *x2p = uu[2].i ; 
            ok = TRUE ; 
          }
      }
  return ok ;
} /* mrnaLocateBestProduct */

/*********************************************************************/
/* mRNA tiling_path */
BOOL mrnaTiling (SC *sc, KEY mrna, BOOL useErrors)
{
  Array tiling = arrayCreate (60, HIT) ;
  Array cfrom = arrayCreate (60, HIT) ;
  Array cleanPieces = arrayCreate (60, HIT) ;
  Array allErrors = arrayCreate (60, HIT), errors = 0 ;
  Array units = arrayCreate (180, BSunit) ;
  Array dnaMrna = dnaGet (mrna) , dnaMrnaR = 0, dnaFixed = 0, dnaTiled = 0 ;
  Array dnaEst = 0 ;
  KEYSET cdsCoveringReads = 0 , mrnaCoveringReads = 0 ;
  int i, j, x, p1, p2, iTiling, jj, iAllError = 0, previousa2 = 0, dx, delta, 
    tilingGap = 0, tilingError = 0, tilingBaseError = 0, amin, amax ;
  HIT *tp = 0, *tpA, *tq, *wp, *zp ;
  KEY clone, previousEst = 0 ;
  BSunit *uu ;
  A_ERR *ep ;
  OBJ Transcript = bsUpdate (mrna), Est = 0 ;
  BOOL done = FALSE, isProduct;
  char *cp, *cq ;

  if (!Transcript || !dnaMrna)
    goto abort ;

  
  if (useErrors)
    {
      bsRemoveTag (Transcript, str2tag("Quality")) ;
    }

  dnaMrnaR = dnaCopy (dnaMrna) ;
  reverseComplement (dnaMrnaR) ;

  if (useErrors)
    dnaFixed = arrayCreate (arrayMax(dnaMrna) + 50, char) ;

  /* recover the hits from the tag Constructed_from */
  units = arrayReCreate (units, 300, BSunit) ;
  bsGetArray (Transcript, _Constructed_from, units, 5) ;
  
  for (i = j = 0 ; i < arrayMax(units) ; i += 5)
    {
      uu = arrp (units, i, BSunit) ;
      if (!strncmp(name(uu[2].k),"NM_",3))
        continue ;
      clone = keyGetKey (uu[2].k, _cDNA_clone) ;
      if (keyFindTag (uu[2].k, str2tag("IS_AM"))) /* the AM cannot be used in the tiling path ! */
        continue ;
      if (0 && /* 2006_11_15   : we no longer care since we remove non informative clones elsewhere */
          keyFindTag (clone, _Suspected_internal_deletion) &&
          !keyFindTag (clone, _Manual_no_internal_deletion))
        continue ; 
      if (!getSeqDna (uu[2].k)) /* in local cache, do not destroy */
        continue ; /* happens if someone destroyed an est during life of database */
      if (uu[1].i < 1) /* happens in case of very low quality 5' read */
        continue ;
      wp = arrayp (cfrom, j++, HIT) ;
      wp->a1 = uu[0].i ;
      wp->a2 = uu[1].i ;
      wp->est = uu[2].k ;
      wp->cDNA_clone = clone ;
      wp->x1 = uu[3].i ;
      wp->x2 = uu[4].i ;
      if (wp->a1 < 1) /* trim debordements */
        {
          if (wp->x1 < wp->x2)
            wp->x1 += 1 - wp->a1 ;
          else
            wp->x1 -= 1 - wp->a1 ;
          wp->a1 = 1 ;
        }
      if (useErrors &&
          keyFindTag (wp->est, _PolyA_after_base) &&
          (Est = bsCreate(wp->est)))
        {
          if (bsGetData (Est, _PolyA_after_base, _Int, &x) &&
              x ==  wp->x2)
            wp->ex2 = x ;
          bsDestroy (Est) ;
        }
    }

  if (!arrayMax(cfrom))
    goto abort ;

  cDNASwapA (cfrom) ;
  arraySort (cfrom, cDNAOrderGloballyByA1) ;
  

  /* count the errors then split into errorless pieces */
  for (j = jj = 0, wp = arrp (cfrom, 0, HIT) ; j < arrayMax(cfrom) ; j++, wp++)
    {
      BOOL isDown = wp->x1 < wp->x2 ? TRUE : FALSE ;

      dnaEst = getSeqDna (wp->est) ; /* in local cache, do not destroy */
      errors = arrayReCreate (errors, 12, A_ERR) ;
      if (useErrors)
        errors = aceDnaDoubleTrackErrors (dnaEst, &(wp->x1), &(wp->x2), isDown, /* may retrofix the coords */
                                          dnaMrna, dnaMrnaR, &(wp->a1), &(wp->a2), 0, errors, 8, 0, FALSE, 0) ;
      /* register the zero error case */
      tp = arrayp (cleanPieces, jj++, HIT) ;
      tp->est = wp->est ;
      tp->cDNA_clone = wp->cDNA_clone ;
      tp->a1 = wp->a1 ; tp->a2 = wp->a2 ;
      tp->x1 = wp->x1 ; tp->x2 = wp->x2 ;
      tp->zone = wp->a2 ;
      tp->nerr = errors ? arrayMax (errors) : 0 ;
      tp->type = 9999 ;
      tp->reverse = isDown ? FALSE : TRUE ;
      tp->ex2 = wp->ex2 ; /* polyA found */

      if (errors && arrayMax(errors))
        for (i = 0, ep = arrp (errors, i, A_ERR) ; i < arrayMax(errors) ; ep++, i++)
          {
            if (ep->iLong > wp->a1 + 8 && /* do not report genomic error seen from extremities */
                ep->iLong < wp->a2 - 8)
              {
                zp = arrayp (allErrors, iAllError++, HIT) ;
                zp->est = wp->est ; 
                zp->a1 = ep->iLong ;
                zp->x1 = ep->iShort ;
                zp->type = ep->type ;
                zp->nerr = ep->baseShort ; 
                zp->clipTop = ep->iLong ;
                zp->clipTop = ep->baseShort ;
              }
            
            tp->x2 = ep->iShort + (isDown ? -1 : 1) ;
            tp->a2 = ep->iLong - 1 ;
            tp->ex2 = 0 ; /* polyA now on next piece, tp->ex1 != 0 && tp->ex2 == 0 */

            tp = arrayp (cleanPieces, jj++, HIT) ;
            tq = tp - 1 ; 
            tp->est = wp->est ;
            tp->ex2 = wp->ex2 ;
            tp->cDNA_clone = wp->cDNA_clone ;
            tp->a1 = ep->iLong + 1 ; tp->a2 = wp->a2 ;  
            tp->x1 = ep->iShort + (isDown ? 1 : -1) ; tp->x2 = wp->x2 ;
            tp->zone = wp->a2 ;
            tp->nerr = ep->baseShort ; 
            tp->type = ep->type ;
            tp->ex1 = 0 ; /* polyA */
            tp->clipTop = ep->iLong ;
            tp->clipEnd = ep->iShort ;
            tp->reverse = isDown ? FALSE : TRUE ;
            
            if (1 || isDown)
              {
                switch (ep->type)
                  {
                  case INSERTION: tp->a1 -= 1 ; break ;
                  case TROU: tp->x1 -= (isDown ? 1 : -1) ; break ;
                  case TROU_DOUBLE: tp->a1 += 1 ; tp->x1 -= (isDown ? 1 : -1) ; break ;
                  case INSERTION_DOUBLE: tp->a1 -= 1 ; tp->x1 +=  (isDown ? 1 : -1) ; break ;
                  default: break ;
                  }
              }
            else
              {
                switch (ep->type)
                  {
                  case INSERTION: tq->a2 += 1 ; break ;
                  case TROU:  tq->x2 -= 1 ; break ;
                  case TROU_DOUBLE: tq->a2 -= 1 ; tq->x2 -= 1 ; break ;
                  case INSERTION_DOUBLE: tq->a2 += 1 ; tq->x2  += 1 ;  break ;
                  default: break ;
                  }
              }
            /* eliminate dense error stetches */
            if (tp->a1 > tp->a2)
              { arrayMax(cleanPieces) -= 1 ; jj-- ; break ; }
            if (tq->a1 > tq->a2 || 
                (isDown && tq->x1 > tq->x2) ||
                (!isDown && tq->x1 < tq->x2))
              {
                *tq = *tp ; tp = tq ;
                arrayMax(cleanPieces) -= 1 ; jj-- ;
              }
          }
    }
  
  cDNASwapA (cleanPieces) ;
  arraySort (cleanPieces, cDNAOrderGloballyByA1) ;
  
  if (useErrors && arrayMax(allErrors))
    {
      arraySort (allErrors, cDNAOrderGloballyByA1) ;
      /* remove non repeated errors */
      for (iTiling = jj = j = 0, zp = arrp (allErrors, 0, HIT) ; j < arrayMax(allErrors) ; )
        {
          wp = zp+1 ; i = 1 ;
          while (i+j < arrayMax(allErrors) && wp->a1 == zp->a1) { i++ ; wp++ ; }
          if (i >= 3)
            {
              iTiling++ ;
              wp = zp ;
              while (wp->a1 == zp->a1)
                {
                  tp = arrp (allErrors, jj++, HIT) ;
                  if (tp != wp)
                    *tp = *wp ;
                  wp++ ;
                }
            }
          j += i ;
          zp += i ;
        }
      if (jj && iTiling < 3) /* more than 3 grouped error means either misplaced gene or junk sequence */
        arrayMax (allErrors) = jj ;
      else
        arrayDestroy (allErrors) ;
    }  
  else
    arrayDestroy (allErrors) ;

  /* select best polyA */ 
  amax = -999999 ;
  for (j = 0, wp = arrp (cleanPieces, 0, HIT) ; j < arrayMax(cleanPieces) ; j++, wp++)
    if (amax < wp->a2) amax = wp->a2 ;

  for (j = 0, tpA = 0, wp = arrp (cleanPieces, 0, HIT) ; j < arrayMax(cleanPieces) ; j++, wp++)
    if (
        wp->ex2 && wp->ex2 > amax - 100 &&
        ( 
         !tpA ||
         (wp->ex2 > tpA->ex2  && wp->a2 - wp->a1 > 8) ||
         (wp->ex2 > tpA->ex2 - 30 && wp->a2 - wp->a1 > tpA->a2 - tpA->a1 + 5)
         )
        )
      tpA = wp ;

  if (tpA)
    for (j = 0, wp = arrp (cleanPieces, 0, HIT) ; j < arrayMax(cleanPieces) ; j++, wp++)
      if (wp->ex2 && wp != tpA)
        wp->ex2 = 0 ;
  
  /* extract a minimum tiling path */
  amin = 999999; amax = -999999 ;
  for (iTiling = j = 0, wp = arrp (cleanPieces, 0, HIT) ; j < arrayMax(cleanPieces) ; j++, wp++)
    {
      if (!iTiling)
        {
          tp = arrayp (tiling, iTiling++, HIT) ;
          previousa2 = wp->a1 ;
          previousEst = 0 ;
          *tp = *wp ;
          continue ;
        }
      if (
          (!tpA && wp->a2 > tp->a2)
          ||
          (
           tpA &&   /* favor promissed polyA */
           (
            (wp->a2 > tp->a2 && (wp->ex2 || !tp->ex2)) || 
            (wp->ex2 && !tp->ex2))
           )
          )
        {
          if (wp->a1 > previousa2)
            {
              if (wp->a1 > tp->a1 ||
                  (wp->a1 == tp->a1 && tp->est == previousEst))
                {
                  previousa2 = tp->a2 ;
                  previousEst = tp->est ;
                  if (wp->a1 < tp->a2)
                    tp->type = 9999 ; /* initial error has been eaten up */
                  tp = arrayp (tiling, iTiling++, HIT) ;
                }
            }
          *tp = *wp ;
          if (tp->a1 < amin) amin = tp->a1 ;
          if (tp->a2 > amax) amax = tp->a2 ;
        }
    }
  
  /* keep happy few */
  for (j = jj = 0, wp = tp = arrp (tiling, 0, HIT) ; j < arrayMax(tiling) ; j++, tp++)
    {
      if (!tp->est)
        continue ;
      if (jj < j) *wp = *tp ;
      wp++ ; jj++ ;
    }
  arrayMax(tiling) = jj ;
  

  /**** find the minimal set covering the cds, avant de couper a bouts francs  ****/
  /* locate best orf */
  p1 = p2 = 0 ;
  isProduct = mrnaLocateBestProduct (Transcript, units, &p1, &p2) ;

  if  (!useErrors)
    {
      mrnaCoveringReads = keySetCreate () ;
      for (i = j = 0, tp = arrp (tiling, 0, HIT) ; i < arrayMax(tiling) ; i++, tp++)
        keySet (mrnaCoveringReads, j++) = tp->est ;
      keySetSort (mrnaCoveringReads) ;
      keySetCompress (mrnaCoveringReads) ;
    }

  if  (isProduct && !useErrors)
    {
      int i1, i2 ;
      cdsCoveringReads = keySetCreate () ;

      i1 = i2 = -1 ;
      for (i = 0, tp = arrp (tiling, 0, HIT) ; i < arrayMax(tiling) ; i++, tp++)
        {
          if (tp->a1 <= p1) i1 = i ; /* last read covering p1 */
          if (tp->a2 >= p2 && i2 == -1) i2 = i ; /* first read covering p2 */
        }

      if (i1 >= 0 && i2 >= 0) /* a covering exists */
        for (i = i1, j = 0,  tp = arrp (tiling, i, HIT) ; i <= i2 ;i++, tp++)
          keySet (cdsCoveringReads, j++) = tp->est ;

      keySetSort (cdsCoveringReads) ;
      keySetCompress (cdsCoveringReads) ;
    }

  /* prepare des bouts francs */
  if (arrayMax(tiling) > 1)
    for (j = 1, tp = arrp (tiling, 1, HIT) ; j < arrayMax(tiling) ; j++, tp++)
      {
        if (tp->a1 > amax - 8)
          { arrayMax(tiling) = j ; break ; }
        wp = tp - 1 ;
        jj = (wp->a2 + tp->a1)/2 ; /* desired cut point */
        if (jj > tp->a2) /* may happen because of polyA */
          jj = (tp->a1 + tp->a2)/2 ;
        if (jj < wp->a2)
          {
            dx = wp->a2 - jj ;
            if (!wp->reverse) wp->x2 -= dx ; 
            else wp->x2 += dx ;
            wp->a2 -= dx ;
          }
        if (jj + 1 > tp->a1)
          {
            dx = jj + 1 - tp->a1 ;
            if (!tp->reverse) tp->x1 += dx ; 
            else tp->x1 -= dx ;
            tp->a1 += dx ;
            tp->type = 9999 ; /* this error has been eaten up */
          }            
      }
  
  /* remove initial wobble defore decale les erreurs */
  if (arrayMax(tiling) > 0)
    for ( i = j = 0, tp = arrp (tiling, 0, HIT) ; j < arrayMax(tiling) ; j++, tp++) /* second non empty tp */
      {
        if (i < 8 && tp->a2 < tp->a1 + 8)
          tp->est = 0 ;
        if (tp->a1 > tp->a2 - 8 &&
            tp->a2 >= arrayMax (dnaMrna))
          tp->est = 0 ;
        if (!tp->est)
          continue ;
        i += tp->a2 - tp->a1 ;
      }
   
  /* keep happy few */ 
  if (arrayMax(tiling) > 0)
    {
      for (j = jj = 0, wp = tp = arrp (tiling, 0, HIT) ; j < arrayMax(tiling) ; j++, tp++)
        {
          if (!tp->est)
            continue ;
          if (jj < j) *wp = *tp ;
          wp++ ; jj++ ;
        }
      arrayMax(tiling) = jj ;
    }

  if (!arrayMax (tiling))
    goto abort ;

  /* tient compte des erreurs en decalant les coordonnees genomiques */
  tilingError = 0 ; 
  if (arrayMax(tiling) > 0)
    for (delta = 0, j = 0, tp = arrp (tiling, 0, HIT) ; j < arrayMax(tiling) ; j++, tp++)
      {        
        if (!tp->est)
          continue ;
        if (sc && useErrors && allErrors)
          mrnaReportGenomicError (sc, tp, allErrors) ;
        
        if (j > 0)
          {
            wp = tp - 1 ;
            
            if (wp->est == tp->est)
              {
                int xx, xa ;
                
                xx = tp->x1 - wp->x2 ;
                if (tp->reverse) xx = -xx ;
                xx -= 1 ;
                xa = tp->a1 - wp->a2 + delta - 1 ;
                if (xx > xa) tilingError += xx ;
                else if (xx < xa) tilingError += xa ;
                else if (xx > 1) tilingError += xx ;
                else if (xx == 1 && tp->nerr != N_) tilingError++ ;
                x = xx - xa ;
              }
            else
              {
                switch (tp->type)
                  {
                  case INSERTION: x = 1 ; tilingError++ ; break ;
                  case TROU: x = -1 ; tilingError++ ; break ;
                  case TROU_DOUBLE: x = -2 ; tilingError += 2 ; break ;
                  case INSERTION_DOUBLE: x = 2 ; tilingError += 2 ; break ;
                  case ERREUR: x = 0 ; tilingError++ ; break ;
                  case AMBIGUE: x = 0 ; break ;
                  default: x = 0 ; break ;
                  }
              }  
            delta += x ;
            
            if (dnaFixed && x >= 0 && tp->type != 9999 && tp->a1 + delta - 2 >= 0)
              {
                array (dnaFixed, tp->a1 + delta - 1, char) = 0 ; /* make room */
                if (x == 0 || x == 1)
                  arr(dnaFixed, tp->a1 + delta - 2, char) = tp->nerr ;
                else if (x > 1)
                  for (i = tp->a1 + delta - 2 ; i >= 0 && i >= tp->a1 + delta - 1 - x ; i--)
                    arr(dnaFixed, i, char) =  N_ ;
              }
          }
        
        if (dnaFixed && tp->a2 + delta - 2 >= 0 && tp->a1 < tp->a2)
          {
            array (dnaFixed, tp->a2 + delta - 1, char) = 0 ; /* make room again */
            for (i = tp->a1 ; i < 1 && i + delta < 1 ; i++) ;
            for (cp = arrp(dnaMrna, i - 1, char), 
                   cq = arrp(dnaFixed, i + delta - 1, char) ;
                 i <= tp->a2 ; i++, cp++, cq++) *cq = *cp ;
          }
        if (isProduct && tp->a1 <= p1 && tp->a2 >= p1) p1 += delta ;
        if (isProduct && tp->a1 <= p2 && tp->a2 >= p2) p2 += delta ;
        tp->a1 += delta ; tp->a2 += delta ; /* coordonnees in fixedDna */
      }


  /* fuse if same est */
 
  if (arrayMax(tiling) > 1)
    for ( j = 1, tp = arrp (tiling, 1, HIT) ; j < arrayMax(tiling) ; j++, tp++) /* second non empty tp */
      {
        wp = tp - 1 ; i = j - 1;
        while (i > 0 && !wp->est) { i--; wp-- ; } /* in case i eat several tp in succession */
        if (i >= 0 && tp->est == wp->est)
          {
            wp->a2 = tp->a2 ; wp->x2 = tp->x2 ; tp->est = 0 ;
          }
      }
  
  
  /* remove initial wobble */
  if (arrayMax(tiling) > 0)
    for ( i = j = 0, tp = arrp (tiling, 0, HIT) ; j < arrayMax(tiling) ; j++, tp++) /* second non empty tp */
      {
        if (i < 8 && tp->a2 < tp->a1 + 8)
          tp->est = 0 ;
        if (!tp->est)
          continue ;
        i += tp->a2 - tp->a1 ;
      }

  /* keep happy few */
  for (j = jj = 0, tp = tq = arrp (tiling, 0, HIT) ; j < arrayMax(tiling) ; j++, tp++)  
    if (tp->est)
      {
        if (tq != tp) *tq = *tp ;
        tq++ ; jj++ ;
      }
  arrayMax (tiling) = jj ;

  if (!arrayMax (tiling))
    goto abort ;

  done = TRUE ;
  /* destroy tile gaps relative to product */
  if (isProduct)
    {
      int i1, i2, j1, j2 ;

      if (! isProduct)
        {
          tp = arrp (tiling, arrayMax(tiling)/2, HIT) ;
          p1 = tp->a1 ; p2 = tp->a2 ;
        }

      i1 = i2 = -1 ;
      for (i = 0, tp = arrp (tiling, 0, HIT) ; i < arrayMax(tiling) ; i++, tp++)
        {
          if (tp->a1 <= p1) i1 = i ; /* last read covering p1 */
          if (tp->a2 >= p2 && i2 == -1) i2 = i ; /* first read covering p2 */
        }

      if (i2 == -1) i2 =  arrayMax(tiling) - 1 ; /* no covering exists */
      if (i1 == -1) i1 = 0 ;

      j1 = i1 ; j2 = i2 ;
      if (i1 < arrayMax(tiling) - 1)
        for (i = i1 + 1, tp = arrp (tiling, i, HIT) ; i <= i2 ;i++, tp++)
          {
            if (tp->a1 > (tp-1)->a2 + 3) break ; /* j1 = last left connex tile */
            else  j1 = i ;
          }
      
      if (i2 > 0)
        for (i = i2 - 1,  tp = arrp (tiling, i, HIT) ; i >= i1 ;i--, tp--)
          {
            if (tp->a2 < (tp+1)->a1 - 3) break ; /* j2 = first right connex tile */
            else  j2 = i ;
          }
      
      if (j1 < j2) /* gap in cds */
        {
          if (FALSE && /* always prefer to tile the 3' end */
              arr(tiling, j1, HIT).a2 - arr(tiling, i1, HIT).a1 >
              arr(tiling, i2, HIT).a2 - arr(tiling, j2, HIT).a1)
            { j2 = j1 ; j1 = i1 ; }
          else
            { j1 = j2 ; j2 = i2 ; }
        }
      else
        { j1 = i1 ; j2 = i2 ; }
      
      {
        int q1 = p1 > arr(tiling, j1, HIT).a1 ? p1 : arr(tiling, j1, HIT).a1 ;
        int q2 = p2 < arr(tiling, j2, HIT).a2  ? p2 : arr(tiling, j2, HIT).a2 ;

        if ( 100 * (q2 - q1) < 80 * (p2 - p1))
          { done = FALSE ; goto abortAM ;} /* equivalent to j1 = j2 = 0 => destroy the whole tiling */
      }
      /* search and destroy gaps outside of j1, j2 */
      
      if (0) /* keep all tilings */
        {
          if (j2+1 < arrayMax(tiling))
            {
              for (i = j2+1,  tp = arrp (tiling, i, HIT) ; i <= arrayMax(tiling) ;i++, tp++)
                if (tp->a1 > (tp-1)->a2 + 3) break ; /* j1 = last left connex tile */
              for (;i <= arrayMax(tiling) ;i++, tp++) 
                tp->est = 0 ; /* destroy in utr region */
            }
          if (j1 - 1 >= 0)
            {
              for (i = j1-1,  tp = arrp (tiling, i, HIT) ; i >= 0 ; i--, tp--)
                if (tp->a2 < (tp+1)->a1 - 3) break ; /* j2 = first right connex tile */
              for (; i >= 0 ; i--, tp--)
                tp->est = 0 ; /* destroy in utr region */
            }
        }
    }


  /* create the tiled dna */
  tilingGap = 0 ; /* used as a signal in the verification phase */
  tilingBaseError = 0 ;
  if (useErrors)
    {
      OBJ obj ;
      KEY dnaTiledKey, dnaTiledDnaKey ;
      int offset = 0,  aliLen = 0 ;
      BOOL ok = FALSE ;

      dnaTiled = mrnaTiledDna (mrna, tiling, &offset, &aliLen, &ok) ;
      if (!ok) done = FALSE ;

      jj = dnaTiled && dnaFixed ?
        (offset + arrayMax(dnaTiled) < arrayMax(dnaFixed) ?
         offset + arrayMax(dnaTiled) : arrayMax(dnaFixed))
        : 0 ;

      lexaddkey (messprintf("%s.am",name(mrna)), &dnaTiledKey, _VSequence) ;
      lexaddkey (messprintf("%s.am",name(mrna)), &dnaTiledDnaKey, _VDNA) ; 
      
      if ((obj = bsUpdate (dnaTiledKey)))
        {
          bsAddTag (obj, str2tag("RefSeqMaker")) ;
          bsAddTag (obj, _Is_AM) ;
          bsAddTag (obj, _Forward) ;
          bsAddTag (obj, _Colour) ;
          bsPushObj (obj) ;
          bsAddTag (obj, _PALEMAGENTA) ;
          bsSave (obj) ;
        }
      
      bsAddKey (Transcript, str2tag("RefSeqMaker"), dnaTiledKey) ;
      if (jj)
        {
          x = 1 + offset ; bsAddData (Transcript, _Constructed_from, _Int, &x) ;
          x = aliLen ; bsAddData (Transcript, _bsRight, _Int, &x) ;
          bsAddKey (Transcript, _bsRight, dnaTiledKey) ;
          x = 1 ; bsAddData (Transcript, _bsRight, _Int, &x) ;
          x = jj - offset ; bsAddData (Transcript, _bsRight, _Int, &x) ;

          cp = arrp (dnaTiled, 0, char) ;
          cq = arrp (dnaFixed, offset, char) ;
          for (i = offset ; i < jj ; i++, cp++, cq++)
            if (*cp && *cq && !(*cp & *cq))
              {
                if (0) messout ("Error at base %d while verifying tiling path %s",
                         i, name(mrna)) ;
                tilingBaseError++ ;
              }
        }
      dnaStoreDestroy (dnaTiledDnaKey, dnaTiled) ; dnaTiled = 0 ;
    }
 abortAM:
  /* export the results */
  if (useErrors)
    {
      for (j = 0, tp = arrp (tiling, 0, HIT) ; j < arrayMax(tiling) ; j++, tp++)
        {
          if (tp->a1 > tp->a2)
            continue ;
          
          x = tp->a1 ;
          bsAddData (Transcript, str2tag ("Tiling_path"), _Int, &x) ;
          x = tp->a2 ;
          bsAddData (Transcript, _bsRight, _Int, &x) ;
          bsAddKey (Transcript, _bsRight, tp->est) ;
          if (0)
            {
              OBJ Est = bsUpdate (tp->est) ;
              
              if (Est)
                {
                  bsAddKey (Est, _Colour, _PALEGREEN) ;
                  bsSave (Est) ;
                }
            }
          x = tp->x1 ;
          bsAddData (Transcript, _bsRight, _Int, &x) ;
          x = tp->x2 ;
          bsAddData (Transcript, _bsRight, _Int, &x) ;
          if (j > 0)
            tilingGap += tp->a1 - 1 - (tp-1)->a2 ;
        }      
      if (tilingGap > 0)
        bsAddData (Transcript, str2tag("Tiling_gap"), _Int, &tilingGap) ;
      if (tilingError > 0)
        bsAddData (Transcript, str2tag("Tiling_error"), _Int, &tilingError) ;
      if (tilingBaseError > 0)
        bsAddData (Transcript, str2tag("Tiling_base_error"), _Int, &tilingBaseError) ;
    }
  else
    {
      if (cdsCoveringReads && keySetMax (cdsCoveringReads) > 0)
        for (i = 0 ; i < keySetMax(cdsCoveringReads) ; i++)
          bsAddKey (Transcript, str2tag("CDS_covered_by"), keySet (cdsCoveringReads, i)) ;
      if (cdsCoveringReads && keySetMax (cdsCoveringReads) == 1)
        bsAddTag (Transcript, str2tag ("Well_supported")) ;
      else 
	bsRemoveTag (Transcript, str2tag ("Well_supported")) ;
      if (mrnaCoveringReads && keySetMax (mrnaCoveringReads) > 0)
        for (i = 0 ; i < keySetMax(mrnaCoveringReads) ; i++)
          bsAddKey (Transcript, str2tag("Mrna_covered_by"), keySet (mrnaCoveringReads, i)) ;
    }

 abort:

  bsSave (Transcript) ;
  arrayDestroy (tiling) ;
  arrayDestroy (cfrom) ;
  arrayDestroy (cleanPieces) ;
  arrayDestroy (allErrors) ;
  arrayDestroy (errors) ;
  arrayDestroy (units) ;
  arrayDestroy (dnaMrna) ;
  arrayDestroy (dnaMrnaR) ;
  arrayDestroy (dnaFixed) ;
  arrayDestroy (dnaTiled) ;
  keySetDestroy (cdsCoveringReads) ;
  keySetDestroy (mrnaCoveringReads) ;
  return done ;
}  /* mrnaTiling  */

/*********************************************************************/

static BOOL mrnaDoAddKantorInfo (KEY mrna, KEY product)
{
  KEY kantor = keyGetKey (product, str2tag("Kantor")) ;
  OBJ Product = 0, Kantor = 0, Mrna = 0 ;
  BSunit *uu ;
  Array units = 0, vnits = 0 ;
  BOOL done = FALSE ;
  float pepW ;
  static KEYSET psortBadTags = 0 ;

  if (getPleaseNoKantorInfo () ||
      !(Kantor = bsUpdate (kantor)))
    goto abort ;
  if (0) printf ("//  calling kantorGetTitles\t%s\t%s\n", name (mrna), name(kantor)) ;
  
  /*** Title ***/
  {
    bsRemoveTag (Kantor, str2tag("ID")) ; /* no lazyness because the gettitle cod eis unstable */
    bsRemoveTag (Kantor, str2tag("Title_hints")) ; /* no lazyness beacuse the gettitle code is unstable */
    if (1)
      {
        char *cp ;
        vTXT bfr = vtxtCreate () ; 
        
        bsSave (Kantor) ;
        cp = kantorGetTitles (bfr, name(kantor), keyFindTag(product,str2tag("Complete"))) ; 
        if (cp && *cp)
          {
            vTXT bfr1 = vtxtCreate () ; 
            
            vtxtPrintf (bfr1, "Kantor %s\n%s\n\n", name(kantor), cp) ;
            parseBuffer (vtxtPtr (bfr1), 0) ;
            vtxtDestroy (bfr1) ;
          }
        Kantor = bsUpdate (kantor) ;
        
        vtxtDestroy (bfr) ;
      }
  }

  /* 2006_11: we now (at long last) compute pI ourselves
  */
  /*** weight and pI ***/
  if (1)
    {
      float weight = 0, pI = 0 ;
      KEY pepKey = 0 ;
      Array pep = 0 ;
          
      if (bsGetKey (Kantor, _Peptide, &pepKey))
        pep = peptideGet (pepKey) ;
          
      if (pep)
	{
	  int w = pepWeight (pep) ;
	  w = (w + 49.0)/100.0 ; 
	  w = 100 * w ;
	  weight = w/1000 ; /* weight = w/10.0 ; 2024_11_28, why should we divide by 10 ? */
	  pI = pepPI (pep) ;
	  if (bsIsTagInClass (class(kantor), str2tag("Michel")))
	    {
	      bsAddData (Kantor, str2tag("Michel"), _Float, &weight) ;
	      bsAddData (Kantor, _bsRight, _Float, &pI) ;
	    }          
	  arrayDestroy (pep) ;
	}
    }

  /*** Import homol ***/

  if ((!mrna || (Mrna = bsUpdate (mrna))) && (Product = bsUpdate (product)))
    {
      KEY _Tax_tree = str2tag("Tax_tree") ;
      KEY _Tax_count = str2tag("Tax_count") ;
      char *commonAncestor = 0 ;
     
      if (bsFindTag (Product, str2tag("ID")))
        bsRemove (Product) ;
      if (bsFindTag (Product, _Title))
        bsRemove (Product) ;
      if (Mrna && bsFindTag (Mrna, _Title))
        bsRemove (Mrna) ;
      if (bsFindTag (Product, str2tag("AKG")))
        bsRemove (Product) ;
      if (bsFindTag (Product, str2tag("Homol")))
        bsRemove (Product) ;
      if (bsFindTag (Product, str2tag("Tax_common_ancestor")))
        bsRemove (Product) ;
      if (bsFindTag (Product, str2tag("Tax_count")))
        bsRemove (Product) ;
      if (bsFindTag (Product, _Pfam))
        bsRemove (Product) ;
      if (bsFindTag (Product, str2tag("Pfam_title")))
        bsRemove (Product) ;
      if (bsFindTag (Product, _Blastp))
        bsRemove (Product) ;
      if (bsFindTag (Product, _Blastp_title))
        bsRemove (Product) ;
      if (bsFindTag (Product, str2tag("Domain")))
        bsRemove (Product) ;
      if (bsFindTag (Product, str2tag("Motif_domain")))
        bsRemove (Product) ;
      if (bsFindTag (Product, str2tag("Localisation")))
        bsRemove (Product) ;
      if (bsFindTag (Product, str2tag("Psort")))
        bsRemove (Product) ;
      if (bsFindTag (Product, str2tag("Psort_title")))
        bsRemove (Product) ;

      if (bsFindTag (Product, str2tag("Expasy")))
        bsRemove (Product) ;
      if (bsFindTag (Product, str2tag("Molecular_weight")))
        bsRemove (Product) ;
      if (bsFindTag (Product, _Tax_tree))
        bsRemove (Product) ; 
      if (bsFindTag (Product, _Tax_count))
        bsRemove (Product) ;
      if (1)
        { /* common ancestor */
          if ((commonAncestor = fogicArguableCommonAncestors (name(kantor)) ))
            {
              bsAddData (Product, str2tag("Tax_common_ancestor"), _Text, commonAncestor) ; 
              messfree (commonAncestor) ;
            } 
        }
      if (1)
        { /* ufo and co */
          char *cp ;
          
          if (bsGetData(Kantor, str2tag("Suggested_class"), _Text, &cp))
            bsAddData (Product, str2tag("Suggested_class"), _Text, cp);
        }
      
      if (1)
        { /* AceKog */
          units = arrayReCreate (units, 800, BSunit) ;
          bsRemoveTag (Product, str2tag("AKG")) ;
	  if (bsGetArray (Kantor, str2tag ("AKG"), units, 4)) /* AceKogHuman Worm Mouse... */
	     {
	       int i ;
	       KEY ka, kb ;
	       ka = str2tag("AceKog") ;
	       kb = str2tag("AceKog_Date") ;
	       for (i = 0 ; i < arrayMax (units) ; i += 4)
		 {
                  uu = arrp (units, i, BSunit) ;
                  if (uu[0].k != ka && uu[0].k != kb &&
		      bsAddTag (Product, uu[0].k)
		      )
		    if (1 && uu[1].s)
		      {
			bsAddData (Product, _bsRight, _Text, uu[1].s) ;
			if (uu[2].s)
			  {
			    bsAddData (Product, _bsRight, _Text, uu[2].s) ;
			    if (uu[3].s)
			      bsAddData (Product, _bsRight, _Text, uu[3].s) ;
			  }
		      }
                }
	     }
          if (bsGetArray (Kantor, str2tag ("AceKog"), units, 8)) /* the homols */
            {
              int i ;
              for (i = 0 ; i < arrayMax (units) ; i += 8)
                {
                  uu = arrp (units, i, BSunit) ;
                  uu[3].i = 3 * uu[3].i - 2 ; 
                  uu[4].i = 3 * uu[4].i ; 
                }
              bsAddArray (Product, str2tag("AceKog"), units, 8) ;
            }
        }
      if (1)
        { /* homol */
          vnits = arrayReCreate (vnits, 8, BSunit) ;
          units = arrayReCreate (units, 800, BSunit) ;
          if (bsGetArray (Kantor, str2tag ("Pfam"), units, 9))
            {
              int i ;
              KEY _Pox_polyA_pol ;

              lexaddkey ("Pox_polyA_pol", &_Pox_polyA_pol, _VPfam) ;
              for (i = 0 ; i < arrayMax (units) ; i += 9)
                {
                  uu = arrp (units, i, BSunit) ;
                  if (_Pox_polyA_pol == uu[0].k)
                    continue ;
                  uu[3].i = 3 * uu[3].i - 2 ; 
                  uu[4].i = 3 * uu[4].i ; 
                  if (bsIsTagInClass (class(product), str2tag("Pfam_title")))
                    {
                      if (i == 0) 
                        bsAddKey (Product, str2tag("Pfam_title"), uu[0].k) ;
                      else
                        bsAddKey (Product, _bsRight, uu[0].k) ;
                    }
                }
              bsAddArray (Product, str2tag("Pfam"), units, 9) ;
              if (bsFindKey (Product, str2tag("Pfam"), _Pox_polyA_pol))
                bsRemove (Product) ;
            } 
          if (bsGetArray (Kantor, str2tag("Blastp"), units, 8))
            {
              int i ;
              for (i = 0 ; i < arrayMax (units) ; i += 8)
                {
                  uu = arrp (units, i, BSunit) ;
                  uu[3].i = 3 * uu[3].i - 2 ; 
                  uu[4].i = 3 * uu[4].i ; 
                }
              bsAddArray (Product, str2tag("Blastp"), units, 8) ;
            } 
        }

      if (1)
        {  /* domain */
          if (bsGetArray (Kantor, str2tag("Domain"), units, 8))
            {
              int i ;
              float j ;
              KEY tag ;
              
              for (i = 0, j= 1 ; i < arrayMax (units) ; j++, i += 8)
                {
                  uu = arrp (units, i, BSunit) ;
                  uu[3].i = 3 * uu[3].i - 2 ; 
                  uu[4].i = 3 * uu[4].i ; 
                  if (lexword2key (name(uu[0].k), &tag, _VSystem) &&
                      bsIsTagInClass(class (product), tag))
                    {
                      bsAddTag (Product, tag) ; 
                      bsAddKey (Product, _bsRight, uu[1].k) ; 
                      bsAddData(Product, _bsRight, _Float, &j) ; 
                      bsAddData(Product, _bsRight, _Int, &(uu[3].i)) ; 
                      bsAddData(Product, _bsRight, _Int, &(uu[4].i)) ; 
                      bsAddData(Product, _bsRight, _Int, &(uu[5].i)) ; 
                      bsAddData(Product, _bsRight, _Int, &(uu[6].i)) ; 
		      if (uu[7].s)
			bsAddData(Product, _bsRight, _Text, (uu[7].s)) ; 
                    }
                  else
                    {
                      if (!psortBadTags || !keySetFind (psortBadTags, uu[0].k, 0))
                        {
                          if (!psortBadTags)
                            psortBadTags = keySetCreate () ;
                          keySetInsert (psortBadTags, uu[0].k) ;
                          messerror ("Unkown tag %s in kantor %s please edit class Product->Psort_domain", 
                                     name(uu[0].k), name(kantor)) ;
                        }
                    }
                }
            } 
        }
      
      if (1)
        { /* domain title */
          Stack s = stackCreate (200) ;
          char *cp ;
          int nn = 0, n ;

          if (bsGetArray (Product,  str2tag("Transmembrane_domain"), units, 7))
            {
              n = arrayMax(units)/7 ;
              if (n == 1)
                catText (s, messprintf("%sa transmembrane domain", 
                                     nn > 0 ? ", " : "")) ;
              else
                catText (s, messprintf("%sat least %d transmembrane domains", 
                                       nn > 0 ? ", " : "", n)) ;
              nn += n ;
            } 
          if (bsGetArray (Product,  str2tag("Coiled_coil_region"), units, 7))
            {
              n = arrayMax(units)/7 ;  
              if (n == 1)
                catText (s, messprintf("%sa coiled coil domain", 
                                       nn > 0 ? ", " : "")) ;
              else
                catText (s, messprintf("%s%d coiled coil domains", 
                                       nn > 0 ? "," : "",n)) ;

              nn += n ;
            } 
                 
          if (bsGetArray (Product,  str2tag("ER_retention_domain"), units, 7))
            {
              n = arrayMax(units)/7 ;
              catText (s, messprintf("%sER retention signal",
                                     nn > 0 ? ", " : "")) ;
              nn += n ;
            } 
          if (0 && bsGetArray (Product,  str2tag("N_myristoylation_domain"), units, 7))
            {
              n = arrayMax(units)/7 ;
              catText (s, messprintf("%sN myristoylation domain", 
                                     nn > 0 ? ", " : "")) ;
              nn += n ;
            } 
          if (0 && bsGetArray (Product,  str2tag("prenylation_domain"), units, 7))
            {
              n = arrayMax(units)/7 ;
              catText (s, messprintf("%sprenylation domain", 
                                     nn > 0 ? ", " : "")) ;
              nn += n ;
            } 
         
          if (0 && bsGetArray (Product,  str2tag("Nucleic_acid_binding"), units, 7))
            {
              n = arrayMax(units)/7 ;
              if (n == 1)
                catText (s, messprintf("%sa nucleic acid binding domain", 
                                       nn > 0 ? ", " : "")) ;
              else
                catText (s, messprintf("%s nucleic acid binding domains", 
                                       nn > 0 ? ", " : "")) ;
              nn += n ;
            } 
          
          bsRemoveTag (Product,  str2tag("Motif_title")) ;
          cp = stackText(s, 0) ;
          if (*cp)
            bsAddData (Product, str2tag("Motif_title"), _Text, cp) ;
          stackDestroy (s) ;
        }

      pepW = 100000 ;

      if (1)
        {  /* expasy */
          bsRemoveTag (Product, _Expasy) ;
          if (bsGetArray (Kantor, str2tag("Michel"), units, 2))
            {
              uu = arrp(units, 0, BSunit) ;
              bsAddData (Product, _Expasy, _Float, & (uu[0].f)) ;
              bsAddData (Product, _bsRight, _Float, & (uu[1].f)) ;
              if (bsFindTag (Product, str2tag("Complete")))
                {
                  bsAddData (Product, str2tag("Molecular_weight"), _Float, & (uu[0].f)) ;
                  pepW = uu[0].f ;
                  bsAddData (Product, _bsRight, _Float, & (uu[1].f)) ;
                }
            }
        }

      if (1)
        { /* psort */
          int cds = 0 ;
          
          bsGetData (Product, _Coding_length, _Int, &cds) ;
          if (bsGetArray (Kantor, str2tag("Localization"), units, 2))
            {
              int i ;
              KEY tag, bestTag = 0 ;
              float membraneScore = 0 ; ;
              
              for (i = 0 ; i < arrayMax(units) ; i += 2)
                {
                  uu = arrp (units, i, BSunit) ;
                  if (lexword2key (name(uu[0].k), &tag, _VSystem) &&
                      bsIsTagInClass (class (product), tag))
                    {        
                      bsAddTag (Product, tag) ;
                      bsAddData (Product, tag, _Float, &(uu[1].f)) ;
                      if (
                          (!strcmp(name(tag),"Nuclear") &&
                           uu[1].f > 59.5 &&
                           pepW < 80) || 
                          (uu[1].f > 40.0 && !strcmp(name(tag),"Golgi")) ||
                          (strcmp(name(tag),"Nuclear") && 
                           uu[1].f > 50.0 )
                          )
                        bestTag = tag ;
                      if (!strcmp(name(tag), "Endoplasmic_reticulum") ||
                           !strcmp(name(tag), "Golgi") ||
                          !strcmp(name(tag), "Plasma_membrane"))
                        membraneScore += uu[1].f ;                      
                    }
                  else
                    {
                      if (!psortBadTags || !keySetFind (psortBadTags, uu[0].k, 0))
                        {
                          if (!psortBadTags)
                            psortBadTags = keySetCreate () ;
                          keySetInsert (psortBadTags, uu[0].k) ;
                          messerror ("Unkown tag %s in class Product->Psort", name(uu[0].k)) ;
                        }
                    }
                }
              if (!bestTag && membraneScore > 59)
                {
                  lexaddkey ("Membrane", &bestTag, 0) ; /* note en dynamique */ 
                  /* pas la peine de l'exporter en dur 8/2/05 danielle
                     bsAddTag (Product, bestTag) ;
                     bsAddData (Product, bestTag, _Float, &membraneScore) ;
                  */
                }
              if (!bsFindTag (Product, _NH2_complete))
                {
                  bsRemoveTag (Product, str2tag ("Cleavable_signal_peptide")) ;
                }
                             
              if (!bestTag && 
                  bsFindTag (Product, str2tag ("Complete")) &&
                  bsFindTag (Product, str2tag ("Cleavable_signal_peptide"))
                  )
                {
                  lexaddkey ("Secreted", &bestTag, 0) ; /* note en dynamique */ 
                }
              if (bestTag && 
                  cds > 210 &&
                  bsFindTag (Product,  str2tag("Complete")) &&
                  (!mrna || !bsFindTag (Mrna, str2tag("Gap"))))
                {
                  char * cp = name(bestTag) ;
                  if (!strcasecmp(name(bestTag),"Cytoplasmic"))
                    cp = "cytoplasmic" ;
                  else if (!strcasecmp(name(bestTag),"Cytoskeletal"))
                    cp = "cytoskeletal" ;
                  else if (!strcasecmp(name(bestTag),"Endoplasmic_reticulum"))
                    cp = "endoplasmic reticulum" ;
                  else if (!strcasecmp(name(bestTag),"Secreted"))
                    cp = "secreted or extracellular" ;
                  else if (!strcasecmp(name(bestTag),"Golgi"))
                    cp = "Golgi" ;
                  else if (!strcasecmp(name(bestTag),"Mitochondrial"))
                    cp = "mitochondrial" ;
                  else if (!strcasecmp(name(bestTag),"Nuclear"))
                    cp = "nuclear" ;
                  else if (!strcasecmp(name(bestTag),"Peroxisomal"))
                    cp = "peroxisomal" ;
                  else if (!strcasecmp(name(bestTag),"Plasma_membrane"))
                    cp = "plasma membrane" ;
                  else if (!strcasecmp(name(bestTag),"Membrane"))
                    cp = "membrane" ;
                  else if (!strcasecmp(name(bestTag),"Vacuolar"))
                    cp = "lysosomal" ;
                  else if (!strcasecmp(name(bestTag),"Vesicles_of_secretory_system"))
                    cp = "secretory vesicle" ;
                  
                  bsAddData (Product, str2tag ("Psort_title"), _Text, cp) ;
                }
            }
        }
      
      if (1)
        { /* tax_tree */            
          bsRemoveTag (Product, _Tax_tree) ;
          if (0 && /* no need to copy that i will rm it for the schema in september 2003 */
              bsGetArray (Kantor, _Tax_tree, units, 30))
            bsAddArray (Product, _Tax_tree, units, 30) ;
          bsRemoveTag (Product, _Tax_count) ;
          if (bsGetArray (Kantor, _Tax_count, units, 5))
            bsAddArray (Product, _Tax_count, units, 5) ;
        }

      if (mrna)
        {
          /* add product in gene */
          int i ;
          KEYSET genes = queryKey (mrna, "{>from_gene; >Gene SMAP} $/ {>From_prediction;>model_of_gene}") ;
          
          for (i = 0 ; i < keySetMax (genes) ; i++)
            bsAddKey (Product, str2tag("GeneBox"), keySet (genes, i)) ;
          keySetDestroy (genes) ; 
        }
    }

  if (1)
    { /* construct full title */
      BOOL unknown = FALSE ;
      char 
        *cv, *famTitle = 0, *motifTitle = 0, *taxblastTitle = 0,
        *psortTitle = 0, *blastpTitle = 0, *bidTitle = 0 ;
      KEY kantorTitle ;
      vTXT buf = vtxtCreate () ;
      
      kantorTitle = 0 ;
      
      bsGetData (Kantor, str2tag("Blastp_title"), _Text, &blastpTitle) ;
      bsGetData (Kantor, str2tag("Family_title"), _Text, &famTitle) ;
      bsGetData (Kantor, str2tag("Taxblast_title"), _Text, &taxblastTitle) ;

      bsGetData (Product, str2tag("psort_title"), _Text, &psortTitle) ;
      bsGetData (Product, str2tag("Motif_title"), _Text, &motifTitle) ;
      
      { /* supersede blastp title by the hand edited GB->Brief_identifcation */
        KEYSET reads = queryKey (product, ">mrna ; >cdna_clone ; >read Brief_identification") ;
        
        if (keySetMax (reads) > 0)
          {
            KEY bid = keyGetKey (keySet(reads,0), _Brief_identification) ;
            
            if (bid)
              bidTitle = name(bid) ;
          }
        keySetDestroy (reads) ;
      }

      if (famTitle)
        bsAddData (Product, str2tag("Family_title"), _Text , famTitle) ;
      if (blastpTitle)
        bsAddData (Product, str2tag("Blastp_title"), _Text, blastpTitle) ;
      if (bidTitle)
        blastpTitle = bidTitle ;
      if (!blastpTitle ||
          (blastpTitle && 
           pickMatch (blastpTitle, "*unknown function*") &&
           !pickMatch (blastpTitle, "*and*")))
        unknown = TRUE ;

      /* create a kantor_short_title */
      if (unknown)
        {
          vtxtPrintf (buf, "putative") ; 
          if (bsFindTag (Product,  str2tag("N_myristoylation_domain")))
            vtxtPrintf (buf, " N-myristoylated") ;
          if (bsFindTag (Product,  str2tag("prenylation_domain")))
            vtxtPrintf (buf, " prenylated") ;
          if (psortTitle &&
              (
               !strstr (psortTitle, "membrane") ||
               !motifTitle ||
               !strstr (motifTitle , "membrane"))
              )
            vtxtPrintf (buf, " %s", psortTitle) ;
          if (bsFindTag (Product,  str2tag("Transmembrane_domain")) && 
              !strstr(vtxtPtr(buf), "membrane") &&
              (!motifTitle || !strstr (motifTitle , "membrane"))
              )
            vtxtPrintf (buf, " membrane") ;
          vtxtPrintf (buf, " protein") ;
          if (famTitle && !strstr(vtxtPtr(buf), "family"))
            vtxtPrintf (buf, " family member") ; 
          if (psortTitle && 
              strstr (psortTitle, "secreted or extracellular"))
            vtxtPrintf (buf, " precursor") ;
          if (motifTitle) 
            vtxtPrintf (buf, ", with %s,", motifTitle) ; /* er coli-coil coil-coil 4  */
          if (taxblastTitle)
            vtxtPrintf (buf, "%s %s", !motifTitle && *taxblastTitle == 'n' ? "," : "", taxblastTitle) ;  
          /* do not create a kantor_title */          
        }
      else
        {  
          vtxtPrintf (buf, "%s", blastpTitle) ;
          if (psortTitle && 
              strstr (psortTitle, "secreted or extracellular"))
            vtxtPrintf (buf, " precursor") ;
          if (famTitle && !strstr(vtxtPtr(buf), "family"))
            vtxtPrintf (buf, " family member") ;

          if (bsFindTag (Product,  str2tag("N_myristoylation_domain")))
            vtxtPrintf (buf, ", possibly N-myristoylated,") ;
          if (bsFindTag (Product,  str2tag("prenylation_domain,")))
            vtxtPrintf (buf, ", possibly prenylated") ;
        }
      

      if (bsFindTag (Product, str2tag ("Short_kantor_title")))
        bsRemove (Product) ;
    
      if (bsFindTag (Product, str2tag ("Kantor_title")))
        bsRemove (Product) ;
      if (bsFindTag (Kantor, str2tag ("Short_kantor_title")))
        bsRemove (Kantor) ;
      if (bsFindTag (Kantor, str2tag ("Kantor_title")))
        bsRemove (Kantor) ;
      cv = vtxtPtr(buf) ;
      if (cv)
        {
          cv = gtCleanUp (cv) ; /* returns a static buffer */
          gtSetUpper  (cv) ;
          bsAddData (Product, str2tag ("Short_kantor_title"), _Text, cv) ;
          bsAddData (Kantor, str2tag ("Short_kantor_title"), _Text, cv) ;
        }

      /* Hand annotations win */
      if (bsGetData (Kantor, str2tag ("Title"), _Text, &cv)) 
        {
          lexaddkey (cv, &kantorTitle, _VText) ;
          bsAddKey (Product, str2tag ("Kantor_title"), kantorTitle) ;
        }
      if (bsGetData (Kantor, str2tag ("Short_title"), _Text, &cv))
        {
          bsAddData (Product, str2tag ("Short_kantor_title"), _Text, cv) ;
        }
      vtxtDestroy (buf) ;
    }
  done = TRUE ;
  
  bsSave (Mrna) ;
  bsSave (Product) ;
  bsSave (Kantor) ;

  { /* export gtMrna or Product base title in Mrna->Title and Product->title */
    /* the essential effect is to add the gene-class descriptors */
    char *cp = kantorReportTitles (name(product)) ;
    if (cp)
      {
        parseBuffer (cp, 0) ;
        messfree (cp) ;
      }
  }

 abort: 
  arrayDestroy (units) ;
  arrayDestroy (vnits) ;
  bsSave (Mrna) ;
  bsSave (Product) ;
  bsSave (Kantor) ;

  return done ;
} /* mrnaDoAddKantorInfo */

/*********************************************************************/

BOOL mrnaAddKantorInfo (KEY mrna)
{
  /* we can only analysed correctly hooked predictions, otherwise DoAddKantorInfo is unsafe */
  KEYSET products = 0 ;
  int ii ;
  BOOL done = FALSE ;

  if (1) /* 0 if you ust want to recompute the NTG/ backfix for the 35g distrib */
    {
      products = queryKey (mrna, "CLASS mrna ; from_gene || from_prediction ; >Product ; Kantor") ;
      ii = keySetMax (products) ;
      
      if (ii > 0) while (ii--)
        /* backwards, to set mrna title to same as top product */
        done |=  mrnaDoAddKantorInfo (mrna, keySet (products, ii)) ;
      
      keySetDestroy (products) ;
      
      /* now add directly to the Product object */
      products = queryKey (mrna, "CLASS Product") ;
      ii = keySetMax (products) ;
      if (ii > 0) while (ii--)
        /* backwards, to set mrna title to same as top product */
        done |=  mrnaDoAddKantorInfo (0, keySet (products, ii)) ;
      keySetDestroy (products) ;
    }

  mrnaSelectBestProduct (mrna) ; /* after kantor */
  mrnaFixNonBestProductTitle (mrna) ;
  mrnaAddFirstMetAndUorf (mrna) ;   /* after select best */

  return done ;
} /* mrnaAddKantorInfo */

/*********************************************************************/

static void mrnaKillDeadProducts(void)
{
  KEYSET ks = query (0, "Find product !From_EST && !mrna") ;
  keySetKill (ks) ;
  keySetDestroy (ks) ;
}

int mrnaAddKeysetKantorInfo (KEYSET ks) 
{
  int nn = 0, ii = keySetMax(ks) ;
  KEY *kp ;
  KEYSET genes = 0, mrnas = 0 ;
  KEYSET oneGenes = keySetCreate () ;

  if (getPleaseNoKantorInfo () || ii <= 0)
    return 0 ;

  mrnas = query (ks, "CLASS mrna") ;
  if ((ii = keySetMax (mrnas)))
    {
      kp = arrp (mrnas, 0, KEY) - 1 ;
      while (kp++, ii--)
        {
          if (mrnaAddKantorInfo (*kp))
            {
              mrnaMrnaScore (*kp) ;               /* after kantor */
              nn++ ;
              if (!(nn % 100))
                {
                  printf ("//  calling mrnaAddKantorInfo n=%5d\t%s\n", nn, name (*kp)) ; 
                }
              if (!(nn % 3000))
                sessionDoSave (TRUE) ;
            }
        } 
    }
  
  genes = query (ks, "{CLASS gene } SETOR {CLASS mrna ; >product ; > geneBox}") ;
  if (keySetMax (genes))
    {
      printf ("//  calling mrnaAddKantorInfo->add pfam to the genes on %d genes\n", keySetMax(genes)) ;
      if ((ii = keySetMax (genes)))
	{
	  kp = arrp (genes, 0, KEY) - 1 ;
	  while (kp++, ii--)
	    {
	      keySet (oneGenes, 0) = *kp ;

	      mcAddKeysetAlterSpliceDetails (0, *kp, FALSE) ;
	      mrnaCountClones (oneGenes, 0, 0) ;

	      giwAddKeysetIntronClass (0, *kp) ;
	      giwAddKeysetProbeWalls (0, *kp) ;
	      
	      mrnaSaveAceKog (oneGenes) ;
	      mrnaSavePfam2Go (oneGenes) ;
	      mrnaSavePsort2Go (oneGenes) ;
	      
	      mrnaSaveGenePastilles (oneGenes) ;
	      nn+=mrnaSaveGeneTitle (oneGenes) ;
	    }
	}
    }
  keySetDestroy (mrnas) ;
  keySetDestroy (genes) ;
  keySetDestroy (oneGenes) ;

  mrnaKillDeadProducts () ;
  return nn ;  
} /* mrnaAddKeysetKantorInfo */

/*********************************************************************/
/*********************************************************************/
/* add a posteriori ORF, CDS, first_atg, first_ntg */
static BOOL mrnaProductGenomeLocate  (KEY mrna, KEY product
                                     , int *a1p, int *a2p
                                     , int x1, int x2
                                     )
{
  OBJ Mrna = bsCreate (mrna) ;
  Array units = arrayCreate (240, BSunit) ; 
  BSunit *uu ;
  int ii, n = 0 ;

  if (bsGetArray (Mrna, _Coding, units, 4))
    for (ii = 0 ; ii < arrayMax(units) ; ii += 4)
      {
        uu = arrp (units, ii, BSunit) ;
        if (uu[2].i == x1 || (uu[2].i < 0 && uu[2].i > -4)) { n++ ; *a1p = uu[0].i ; } 
        if (uu[3].i == x2) { n++ ; *a2p = uu[1].i ; }
      }
  bsDestroy (Mrna) ;
  arrayDestroy (units) ;
  return n == 2 ? TRUE : FALSE ;
} /* mrnaProductGenomeLocate */

/*********************************************************************/
/*********************************************************************/
/* add a posteriori ORF, CDS, first_atg, first_ntg */
static int mrnaDoAddFirstMet (KEY mrna, KEY product, int x1, int x2
                              , int *cdsp, int *u5p, int *u3p, BOOL *bp)
{
  Array dna = dnaGet (mrna) ;
  OBJ Product = 0 ;
  char *cp, ntgType [8] ;
  int i, first_atg = 0, first_ntg = 0, hasStop = 0, hasUpStop = 0, openLength, ntgLength ;
  BOOL upGap = FALSE ;

  strcpy (ntgType, "...NNN.") ;
  *cdsp = *u3p = *u5p = 0 ; *bp = FALSE ;
  if (!dna)
    goto abort ;
  if (arrayMax(dna) < x2)
    goto abort ;
  if ((x2 - x1 + 1) % 3)
    {
      messerror ("mrnaDoAddFirstAtg : x2-x1+1 = %d !== 0 mod(3) , mrna=%s product = %s"
                 , x2 - x1 + 1, name(mrna), name(product)) ;
      goto abort ;
    }

  /* find first atg downstream of base 1 */
  /* count first xtg in bp origin 0 */
  first_atg = first_ntg = -1 ;
  if (x1 > 0)
    for (i = x1 - 1, cp = arrp (dna, i, char); i < x2 - 2; i += 3, cp += 3)
      if (*cp == A_ && *(cp+1) == T_ && *(cp+2) == G_)
        { first_atg = i - x1 + 1 ; break ; }
  /* look for gap and last stop */
  
  if (x1 >= 4)
    {
      for (i = x1 - 4, cp = arrp (dna, i, char) ; i >= 0 ; i -= 3, cp -= 3)
        { 
          if (codon(cp) == '*')
            { hasUpStop = i - x1 - 1 ; break ; }
          if (*cp == N_ && *(cp+1) == N_ && *(cp+2) == N_)
            { upGap = TRUE ; break ; }
        }
      /* if first atg == 1, look for first_ntg upstream, down of last stop */
      if (! upGap && first_atg == 0 && x1 >= 4)
        {
          if (i < 0) { i += 3 ; cp += 3 ; }
          for ( ; i < x1 - 1 ; i += 3, cp += 3)
            { 
              if (*(cp+2) == G_ && 
                  ((*(cp-3) == G_ && *(cp+3) == G_) || (*(cp-3) == A_)) &&
                  (*(cp) == A_ || *(cp+1) == T_)
                  )
                {
                  first_ntg = i - x1 + 1 ;
                  ntgType[0] = ace_lower (dnaDecodeChar[(int) *(cp-3)]) ;
                  ntgType[1] = '.' ;
                  ntgType[2] = '.' ;
                  ntgType[3] = ace_upper (dnaDecodeChar[(int) *cp]) ;
                  ntgType[4] = ace_upper (dnaDecodeChar[(int) *(cp+1)]) ;
                  ntgType[5] = ace_upper (dnaDecodeChar[(int) *(cp+2)]) ;
                  ntgType[6] = *(cp+3) == G_ ? 'g' : '.' ;
                  ntgType[7] = 0 ;
                  break ;
                }
            }
        }
    }
  /* else look for first_ntg downstream */
  if (!upGap && first_atg != 0 && x1 >= 1)
    for (i = x1 - 1, cp = arrp (dna, i, char); i < x2 - 2; i += 3, cp += 3)
      if (*(cp+2) == G_ && 
          ((*(cp-3) == G_ && *(cp+3) == G_) || (*(cp-3) == A_)) &&
          (*(cp) == A_ || *(cp+1) == T_)
          )
        {
          first_ntg = i - x1 + 1 ;
          ntgType[0] = ace_lower (dnaDecodeChar[(int) *(cp-3)]) ;
          ntgType[1] = '.' ;
          ntgType[2] = '.' ;
          ntgType[3] = ace_upper (dnaDecodeChar[(int) *cp]) ;
          ntgType[4] = ace_upper (dnaDecodeChar[(int) *(cp+1)]) ;
          ntgType[5] = ace_upper (dnaDecodeChar[(int) *(cp+2)]) ;
          ntgType[6] = *(cp+3) == G_ ? 'g' : '.' ;
          ntgType[7] = 0 ;
          break ;
        }

  if (x2 >= 3)
    {
      cp = arrp (dna, x2 - 3, char) ;
      if (codon(cp) == '*')
        hasStop = 3 ;
    }

  openLength = x2 - x1 + 1 - hasStop ;
  ntgLength = openLength - first_ntg ;
  /* count first xtg in aa, origin 1 */
  if (first_atg >= 0)
    first_atg = 1 + first_atg/3 ;
  else 
    first_atg = 0 ;
  if (first_ntg < -1)
    first_ntg -= 3 ; /* Plato */
  if (first_ntg != -1)
    first_ntg = 1 + first_ntg/3 ;
  else 
    first_ntg = 0 ;
  if (ntgLength < MINI_LEU2END)
    first_ntg = ntgLength = 0 ;
  if (first_atg > 0 && first_atg - first_ntg < MINI_LEU2MET/3)
    first_ntg = ntgLength = 0 ;
  if (upGap)
    first_atg = first_ntg = 0 ;
  

  /* report results */
  if ((Product = bsUpdate (product)))
    {
      *bp = bsFindTag (Product, _Best_product) ;
      if (first_atg)
        bsAddData (Product, _First_ATG, _Int, &first_atg) ;
      else if (bsFindTag (Product, _First_ATG))
        bsRemove (Product) ;

      if (bsFindTag (Product, _First_NTG))
        bsRemove (Product) ;
      if (first_ntg)
        {
          bsAddData (Product, _First_Kozak, _Int, &first_ntg) ;
          bsAddData (Product, _bsRight, _Text, ntgType) ;
        }
      else if (bsFindTag (Product, _First_Kozak))
        bsRemove (Product) ;

      if (first_atg == 1 || (first_ntg == 1 && bsFindTag (Product, _NH2_complete)))
        bsAddTag (Product, _At_position_1) ;
      else if (bsFindTag (Product, _At_position_1))
        bsRemove (Product) ;

      if (bsFindTag (Product, _NH2_complete)) /* i must use ntg or atg */
        {
	  if (first_ntg > 0 && 
              (
               (first_atg == 0 &&  ntgLength > MINI_LEU2END) ||
               first_ntg <= first_atg -  MINI_LEU2MET/3)
              )
            *cdsp = x2 - x1 + 1 - hasStop - 3 * (first_ntg - 1) ;
          else if (first_atg > 0)
            *cdsp = x2 - x1 + 1 - hasStop - 3 * (first_atg - 1) ;
          else
            *cdsp = 0 ;
        }
      else /* export from first bp */
        {
          *cdsp = openLength ;
          if (openLength > 3)
            {
              i = openLength ;
              bsAddData (Product, _Open_length, _Int, &i) ;
              cp = "bp" ; bsAddData (Product, _bsRight, _Text, cp) ;
              i /= 3 ; bsAddData (Product, _bsRight, _Int, &i) ;
              cp = "aa" ; bsAddData (Product, _bsRight, _Text, cp) ;
            }
        }
      
      if (*cdsp > 3)
        {
          i = *cdsp ;
          bsAddData (Product, _Coding_length, _Int, &i) ;
          cp = "bp" ; bsAddData (Product, _bsRight, _Text, cp) ;
          i /= 3 ; bsAddData (Product, _bsRight, _Int, &i) ;
          cp = "aa" ; bsAddData (Product, _bsRight, _Text, cp) ;
          
          *u5p = x1 - 1 ;
          *u3p = arrayMax (dna) - x2 ;        
          if (*u5p > 2)
            bsAddData (Product, _Length_5prime_UTR, _Int, u5p) ;
          else if (bsFindTag (Product, _Length_5prime_UTR))
            bsRemove (Product) ;
          if (*u3p > 2)
            bsAddData (Product, _Length_3prime_UTR, _Int, u3p) ;
          else if (bsFindTag (Product, _Length_3prime_UTR))
            bsRemove (Product) ;
        }
      else
	{
	  if (bsFindTag (Product, _Coding_length))
            bsRemove (Product) ;
	  if (bsFindTag (Product, _Length_5prime_UTR))
            bsRemove (Product) ;
	  if (bsFindTag (Product, _Length_3prime_UTR))
            bsRemove (Product) ;
	}
      if (hasStop)
        {
          int x = x2 - x1 - 1 ;
          bsAddData (Product, str2tag("Down_stop"), _Int, &x) ;
        }
      else if (bsFindTag (Product, str2tag("COOH_complete")))
        bsRemove (Product) ;
        
      {
        int x = 0 ;
        KEY dummy ;

        bsGetData (Product, str2tag("Up_stop"), _Int, &x) ;
        if (x && ! hasUpStop && x + x1 < 0 &&  x > -9) /* allow previously performed stealing */
          {
            if (bsFindTag (Product, str2tag("Up_stop")))
              bsRemove (Product) ;
            if (bsFindTag (Product, str2tag("NH2_Complete")) &&
                ! bsGetKeyTags (Product, _bsRight, &dummy))
              {
                if (bsFindTag (Product, str2tag("NH2_Complete")))
                  bsRemove (Product) ;
              }
          }
        if (!x && hasUpStop )
          bsAddData (Product, str2tag("Up_stop"), _Int, &hasUpStop) ;
      }
        
      if (upGap)
        {
          if (bsFindTag (Product, str2tag("NH2_Complete")))
            bsRemove (Product) ;
        }

      if (bsFindTag (Product, str2tag("NH2_Complete")) &&
          bsFindTag (Product, str2tag("COOH_complete")))
        {
          if (bsFindTag (Product,  str2tag("Partial")))
            bsRemove (Product) ;
          bsAddTag  (Product, str2tag("Complete")) ;
        }
      else 
        {
          if (bsFindTag (Product, str2tag("Complete")))
            bsRemove (Product) ;
          bsAddTag (Product, str2tag("Partial")) ;
        }
        
      bsSave (Product) ;
    }

 abort:
  arrayDestroy (dna) ;

  return first_atg ;
} /* mrnaDoAddFirstMet */

/*************************/

static void mrnaAddFirstMetAndUorf (KEY mrna)
{
  OBJ Mrna = 0 ;
  Array units = arrayCreate (50, BSunit) ;
  BSunit *uu, *vv ;
  char *cp ;
  int i, ii, met, cds, a1, a2, x1, x2, mrnaScore, besti, qual, bestQual,upstop ;
  int u3p, u5p ;
  KEY product ;
  BOOL bestProduct ;

  if ((Mrna = bsUpdate (mrna)))
    { 
      if (bsFindTag (Mrna, _Length_3prime_UTR))
        bsRemove (Mrna) ;
      if (bsFindTag (Mrna, _Length_5prime_UTR))
        bsRemove (Mrna) ;
      if (bsFindTag (Mrna, _Longest_CDS))
        bsRemove (Mrna) ; 
      units = arrayReCreate (units, 50, BSunit) ;
      if (bsGetArray (Mrna, _Product, units, 5))
        for (ii = 0 ; ii < arrayMax(units) ; ii += 5)
          {
            uu = arrp (units, ii, BSunit) ;
            product = uu[0].k ;
            x1 = uu[1].i ; x2 = uu[2].i ; a1 = a2 = 0 ;
            met = mrnaDoAddFirstMet (mrna, product, x1, x2, &cds, &u5p, &u3p, &bestProduct) ;
            /* unspliced coords of the product in the mRNA */
            if (mrnaProductGenomeLocate (mrna, product, &a1, &a2, x1, x2))
              { uu[3].i = a1 ; uu[4].i = a2 ; }
            else
              uu[3].i = uu[4].i = 0 ;
            if (bestProduct)
              {
                bsAddData (Mrna, _Longest_CDS, _Int, &cds) ;
                cp = "bp" ; bsAddData (Mrna, _bsRight, _Text, cp) ;
                cds /= 3 ; bsAddData (Mrna, _bsRight, _Int, &cds) ;
                cp = "aa" ; bsAddData (Mrna, _bsRight, _Text, cp) ;
                if (u5p > 2)
                  bsAddData (Mrna, _Length_5prime_UTR, _Int, &u5p) ;
                if (u3p > 2)
                  bsAddData (Mrna, _Length_3prime_UTR, _Int, &u3p) ;
              }
            if (!met)
              continue ;
            
            /* search uORF */
            if (!bsGetData (Mrna, str2tag("Score"), _Int, &mrnaScore) ||
                mrnaScore < 10 || !bsFindTag (Mrna, str2tag("Tiling_path"))) 
              continue ;
            if (!keyFindTag (product, str2tag("COOH_complete")))
              continue ;
            cds = (x2 - x1 + 1)/3 - met + 1 ;
            if (cds > 8 && cds < 200) /* expected number of aa of the uORF */
              {
                BOOL okUp = TRUE, okDown = FALSE ;
                /* am i the first product in my own frame */
                for (i = ii - 3, vv = uu - 3 ; i >= 0 ; i -= 3, vv -= 3)
                  if (!(x1 - vv[1].i % 3))
                    okUp = FALSE ;
                /* who is the best product + upstop */
                bestQual = besti = upstop = 0 ;
                for (i = 0, vv = arrp (units, i, BSunit) ; i < arrayMax(units) ; i += 3, vv += 3)
                  {
                    OBJ Product = bsCreate (vv[0].k) ;
                    
                    if (Product)
                      {
                        if (bsGetData (Product, str2tag("Quality"), _Int, &qual) &&
                            qual > bestQual && qual > 2)
                          { 
                            bestQual = qual ; besti = i ;
                            upstop = 0 ;
                            bsGetData (Product, str2tag("Up_stop"), _Int, &upstop) ;
                          }
                        bsDestroy (Product) ;
                      }
                  }
                /* is there a candidate long gene down stream */
                for (i = ii + 3, vv = uu + 3 ; 
                     upstop && bestQual && okUp && !okDown && i < arrayMax(units) && i <= besti ;
                     i += 3, vv += 3)
                  if (i == besti && vv[1].i /* + upstop */ > x2 && (vv[2].i - vv[1].i > 2 *(x2 - x1)))
                    okDown = TRUE ;
                if (bestQual && okUp && okDown)
                  {
                    OBJ Product = bsUpdate(product) ;
                    
                    if (Product)
                      {
                        bsAddTag (Product, str2tag("uORF")) ;
                        bsSave (Product) ;
                      }
                  }
              }
          }
      if (bsFindTag (Mrna, _Product))
        {
          bsRemove (Mrna) ;
          if (arrayMax (units))
            bsAddArray (Mrna, _Product, units, 5) ;
        }
    }
  bsSave (Mrna) ;
  arrayDestroy (units) ;
  return ;
} /* mrnaAddFirstMetAndUorf */

 /*********************************************************************/
/* non best product should not inherit the gene title */
static void mrnaFixNonBestProductTitle (KEY mrna)
{
  char *title ;
  KEY tt, product ;
  int ii ;
  OBJ Product ;
  KEYSET products = queryKey (mrna, ">Product ; ! Best_product") ;

  for (ii = 0 ; ii < keySetMax (products) ; ii++)
    {
      product = keySet (products, ii) ;
      if ((Product = bsUpdate (product)))
        {
          if (bsGetData (Product, str2tag("Short_kantor_title"), _Text, &title))
            {
              lexaddkey (title, &tt, _VText) ;
              bsAddKey (Product, _Title, tt) ;
            }
          bsSave (Product) ;
        }
    }
  keySetDestroy (products) ;
} /* mrnaFixNonBestProductTitle */

/*********************************************************************/

static BOOL mrnaAddOneTiling (SC *sc, KEY tr)
{
  BOOL ok1, ok2 ;

  ok1 = mrnaTiling (sc, tr, TRUE) ;
  ok2 = mrnaTiling (sc, tr, FALSE) ;
  mrnaFindClone2Resequence (tr) ;

  return ok1 && ok2 ;
}

/*********************************************************************/

int mrnaAddKeysetTiling (KEYSET ks) 
{
  int nn = 0, ii = keySetMax(ks) ;
  KEY *kp ;

  if (ii <= 0)
    return 0 ;

  kp = arrp (ks, 0, KEY) - 1 ;
  while (kp++, ii--)
    {
      if (mrnaAddOneTiling (0, *kp))
        {
          nn++ ;
          if (!(nn % 100))
            {
              printf ("//  calling mrnaTiling n=%5d\t%s\n", nn, name (*kp)) ; 
            }
          if (!(nn % 1000))
            {
              getSeqDna ( KEYMAKE (_VCalcul, 12345)) ;    /* cleanup */
            }
        }      
    } 

  if (1)
    {
      KEYSET genes = query (ks, ">from_gene") ;

      if (keySetMax (genes))
        for (ii = 0, kp = arrp (genes, 0, KEY) ; ii < keySetMax(genes) ; ii++, kp++)
          mrnaSaveMrnaHierarchy (*kp) ;
      keySetDestroy (genes) ;
    }
  return nn ;  
}

/*********************************************************************/
/* a1, a2 is the intersection region */
static int mrnaCountCloneError (int a1, int a2, Array errDown, Array errUp)
{
  int nerr = 0, n = 0, nbest = 0, i, j ;
  A_ERR *e1, *e2 ;
  int n1 = arrayMax (errDown) ;
  int n2 = arrayMax (errUp) ;

  i = j = 0 ; e1 = e2 = 0 ;
  if (n1) /* only use errDown err up to a1 */
    for (i = 0, e1 = arrp (errDown, 0, A_ERR) ; i < n1 && e1->iLong <= a1 - 2 ; i++, e1++)
      switch  (e1->type)
        {
        case AMBIGUE : break ;
        default: nerr++ ; break ;
        }
  if (n2) /* only use errUp  err up to a2 */
    for (j = n2 - 1, e2 = arrp (errUp, j, A_ERR) ; j >= 0 && e2->iLong >= a1 - 1 ; e2--, j--)
      switch  (e2->type)
        {
        case AMBIGUE : break ;
        default: nbest++ ; break ;
        }
  if (n2)
    for (j = n2 - 1, e2 = arrp (errUp, j, A_ERR) ; j >= 0 && e2->iLong >= a2  ; e2--, j--)
      switch  (e2->type)
        {
        case AMBIGUE : break ;
        default: nerr++ ; nbest-- ; break ;
        }
  
  if (n1 && n2 && e1->iLong <= e2->iLong) /* add the errors in the intersect */
    {
      if (j < 0) j++ ;
      e2 = arrp (errUp, j, A_ERR) ;
      while (j >= 0 && i < n1 && e2->iLong <= a2 - 1 && e1->iLong <= a2 - 1)
        {
          /* N do not count as errors, 
           * n++ before n-- so that we do not forget a snip seen on both strands 
           */
          if (e1->type == AMBIGUE) { i++ ; e1++ ; }
          else if (e2->type == AMBIGUE) { j-- ; e2-- ; }
          else if (e2->iLong <  e1->iLong) { i++ ; n++ ; e1++ ;}
          else { j-- ; n-- ; e2-- ; }

          if (n < nbest) nbest = n ;
        }
      nerr += nbest ; 
    }                     

  return nerr ;
} /* mrnaCountCloneError */

/*********************************************************************/

void mrnaSaveMrnaCoding (SMRNA *smrna)
{
  KEY tr = smrna->tr ;
  OBJ Transcript = bsUpdate (tr) ;
  KEY map = 0 ;
  int dmap1 = 0 ;
  int dmap2 = 0 ;
  int givenBack = smrna->givenBack ;
  int stolen = smrna->stolen ;

  if (bsGetKey (Transcript, _IntMap, &map) &&
      bsGetData (Transcript, _bsRight, _Int, &dmap1))
    bsGetData (Transcript, _bsRight, _Int, &dmap2) ;
  /*** Coding, comes after exporting the products ! ***/
  if (smrna->orfs && arrayMax(smrna->orfs))
  {
    int ii, nn, iProd, a1, a2, x1, x2 ;
    HIT *up, *vp, *prod = 0 ;
    ORFS *orf ;
    BSunit *uu ;
    OBJ Product = 0 ;
    Array prods = arrayCreate (12, HIT) ;  /* coords of products */
    Array units = 0, cHits = 0 , hits = 0 ;
    Stack s = stackCreate (64) ;

    /* collect the coords of all products */
    orf = arrp (smrna->dnas, smrna->bestDna, ORFS) ;  
    hits = orf->hits ;

    cHits =  arrayCopy (hits) ;
    units = arrayReCreate (units, 12, BSunit) ;
    bsGetArray (Transcript, _Product, units, 3) ;
    for (ii = 0, iProd = 0 ; ii < arrayMax (units) ; ii += 3)
      {
        uu = arrp (units, ii, BSunit) ;
        prod = arrayp (prods, iProd++, HIT) ;
        prod->est = uu[0].k ;
        prod->x1 = uu[1].i + givenBack ; prod->x2 = uu[2].i  + givenBack ;
      }

    /* subdivide the exons on product boundaries */
    for (iProd = 0 ; iProd < arrayMax(prods) ; iProd++)
      {
        prod = arrp (prods, iProd, HIT) ;
        for (ii = 0 ; ii < arrayMax (cHits) ; ii++)
          {
            up = arrayp (cHits, ii, HIT) ;
            if (!ii && stolen && (up->type & gX)) up->x1 -= stolen ;
            if (up->x1 < prod->x1 && up->x2 >= prod->x1)
              {
                vp = arrayp (cHits, arrayMax(cHits) , HIT) ;
                up = arrp (cHits, ii, HIT) ; /* possible relocation of cHits */
                *vp = *up ;
                nn = prod->x1 - up->x1 ;
                up->a2 = up->a1 + nn - 1 ; vp->a1 = up->a1 + nn ;
                up->x2 = prod->x1 - 1 ; vp->x1 = prod->x1 ;
              }
            if (up->x1 <= prod->x2 && up->x2 > prod->x2)
              {
                vp = arrayp (cHits, arrayMax(cHits) , HIT) ;
                up = arrp (cHits, ii, HIT) ; /* possible relocation of cHits */
                *vp = *up ;
                nn = prod->x2 - up->x1 + 1 ;
                if (!ii && stolen && (up->type & gX)) nn -= stolen ;
                up->a2 = up->a1 + nn - 1 ; vp->a1 = up->a1 + nn ;
                up->x2 = prod->x2 ; vp->x1 = prod->x2 + 1 ;
              }
          }
      }

    arraySort (cHits, cDNAOrderGloballyByA1) ;

    /* export the results */
    for (ii = 0 ; ii < arrayMax (cHits) ; ii++)
      {
        up = arrp (cHits, ii, HIT) ;
         
        if (! (up->type & gX) && /* pick also the ORF_GAP */
            !(prod && prod->x1 <= up->x1 && prod->x2 >= up->x2)) 
          continue ; 
          
        a1 = up->a1 + stolen ; a2 = up->a2 + stolen ;
        if (!ii) a1 -= stolen ;
        x1 = up->x1 - givenBack ; x2 = up->x2 - givenBack ;
        if (a2 - a1 >= 0)
          {
            bsAddData (Transcript, _Coding, _Int, &a1) ;
            bsAddData (Transcript, _bsRight, _Int, &a2) ;
            bsAddData (Transcript, _bsRight, _Int, &x1) ;
            bsAddData (Transcript, _bsRight, _Int, &x2) ;
            bsPushObj (Transcript) ;

            stackClear (s) ;
            for (iProd = 0 ; iProd < 25 && iProd < arrayMax(prods) ; iProd++)
              {
                prod = arrp (prods, iProd, HIT) ;
                if (up->x2 < prod->x1)
                  catText (s, messprintf ("5%c ",'A' + iProd)) ;
                else if (up->x1 > prod->x2)
                  catText (s, messprintf ("3%c ",'A' + iProd)) ;
                else
                  {
                    if (!prod->a1 || a1 < prod->a1) prod->a1 = a1 ;
                    if (!prod->a2 || a2 > prod->a2) prod->a2 = a2 ;
                    catText (s, messprintf ("%c ",'A' + iProd)) ;
                  }
              }
            bsAddData (Transcript, _bsRight, _Text, stackText (s, 0)) ;
            bsGoto (Transcript, 0) ;
          }
      }
 
    for (iProd = 0 ; iProd < arrayMax(prods) ; iProd++)
      if (map && dmap1 && dmap2 && (Product = bsUpdate (prod->est)))
        /* this is bugged, because i do not include the lenghts of the introns
           it should be done inside the coding section 
           but corrected for stolen length
        */
        {
          int m1, m2 ;
          bsAddKey (Product, _IntMap, map) ;
          if (dmap1 < dmap2)
            { m1 = prod->a1 + dmap1 - 1 - stolen ; m2 = prod->a2 + dmap1 - 1 ; }
          else
            { m1 = dmap1 - prod->a1 + 1 + stolen ; m2 = dmap1 - prod->a2 + 1 ; }
          bsAddData (Product, _bsRight, _Int, &m1) ;
          bsAddData (Product, _bsRight, _Int, &m2) ;
          bsSave (Product) ;
        }
          
    stackDestroy (s) ;
    arrayDestroy (prods) ;
    arrayDestroy (cHits) ;
  }
  bsSave (Transcript) ;

  return ;
} /* mrnaSaveMrnaCoding */

/*********************************************************************/

void mrnaSaveMrna (S2M *s2m, SC* sc, Array estHits, Array smrnas, SMRNA *smrna, BOOL composite)
{
  KEY tr = smrna->tr ;
  OBJ Transcript = 0 ;
  KEY map = 0, cosmid = s2m->cosmid1 ? s2m->cosmid1 : s2m->cosmid ;
  OBJ Cosmid = 0 ;
  int ii, c1 = 0, c2 = 0, dmap1 = 0, dmap2 = 0, givenBack = 0, stolen = 0 ;
  int ga1 = sc->a1, ga2 = sc->a2 ;
  KEY  relevantCosmid = mrnaRelevantCosmid (s2m->cosmid, s2m->cosmid1, s2m->linkPos, ga1, ga2) ;
  BOOL isUp = sc->isUp ;
  int ta1, ta2 ; 
  Array feet = 0, units = 0 ;
  Array cloneHits = 0 ;
  int mrnaLength = 0 ;
  KEYSET nonSpecificClones = 0 ;
  HIT *prod_p ;
  Array prod_array = 0 ;
  int minIntronSize = getPleaseMinIntronSize () ;

  cDNASwapA (estHits) ;
  arraySort (estHits, cDNAOrderGloballyByA1) ;
  if (isUp)
    { ta1 = ga1 - smrna->a1 + 1 ; ta2 = ga1 - smrna->a2 + 1 ; }
  else
    { ta1 = smrna->a1 + ga1 - 1 ; ta2 = smrna->a2 + ga1 - 1 ; }


  Transcript = bsUpdate (tr) ;
  if (!Transcript)
    return ;

  /*** Bad quality ***/

  if (smrna->cGroup)
    bsAddTag(Transcript, str2tag("Bad_Quality")) ;


  /*** Transcribed_from ***/
  {
    int ll ;    

    if (ta2 > ta1) ll = ta2 - ta1 + 1 ;
    else ll = ta1 - ta2 + 1 ;

    bsAddData (Transcript, _Covers, _Int, &ll) ;
    bsAddData (Transcript, _bsRight, _Text, "bp from") ;
    bsAddData (Transcript, _bsRight, _Int, &ta1) ;
    bsAddData (Transcript, _bsRight, _Text, "to") ;
    bsAddData (Transcript, _bsRight, _Int, &ta2) ;
  
    
    if ((Cosmid = bsCreate (cosmid)))
      {
        if (bsGetKey (Cosmid, _IntMap, &map) &&
            bsGetData (Cosmid, _bsRight, _Int, &c1)  &&
            bsGetData (Cosmid, _bsRight, _Int, &c2))
	  {}  ;
        bsDestroy (Cosmid) ;
      }
    if (map && c1 + c2)
      {
        if (c1 < c2)
          { dmap1 = c1 + ta1 - 1 ; dmap2 = c1 + ta2 - 1 ; }
        else
          { dmap1 = c1 - ta1 + 1 ; dmap2 = c1 - ta2 + 1; }
        bsAddKey (Transcript, _IntMap, map) ;
        bsAddData (Transcript, _bsRight, _Int, &dmap1) ;
        bsAddData (Transcript, _bsRight, _Int, &dmap2) ;
      }
  }

  /*** SMAP: Splicing ***/
  {
    int nExons, nIntrons, gapl, ll, ifeet = 0, ngt_ag = 0, ngc_ag = 0,  nv_rep = 0, nfuzzy = 0 ;
    HIT *up, *vp, *wp, *fp; ORFS *orf ;
    Array hits = 0 ;
    KEY _v_rep = str2tag ("V_repeat"), actualFoot ;
    KEY _Fuzzy_gt_ag = str2tag ("Fuzzy_gt_ag") ;
    KEY _Fuzzy_gc_ag = str2tag ("Fuzzy_gc_ag") ;

    units = arrayReCreate  (units, 60, BSunit) ;
    nExons = nIntrons = mrnaLength = gapl = 0 ;
     
    orf = smrna->dnas && smrna->bestDna >=0 && smrna->bestDna < arrayMax (smrna->dnas) ? arrp (smrna->dnas, smrna->bestDna, ORFS) : 0 ;  
    hits = orf ? orf->hits : 0 ;

    for (ii = 0, up = hits ? arrp (hits, 0, HIT) : 0  ; hits && ii < arrayMax (hits) ; up++, ii++)
      {
        int jj ;
	bsAddData (Transcript, _Splicing, _Int, &(up->a1)) ;
	wp = up ;
	for (vp = up + 1, jj = ii + 1 ; jj < arrayMax (hits) ; vp++, jj++)
	  if ((up->type & gX) && vp->a1 == up->a2 + 1 && vp->type == up->type)
	    wp = vp ;
	  else
	    break ;
	bsAddData (Transcript, _bsRight, _Int, &(wp->a2)) ;
        if ((up->type & gI) && minIntronSize && up->a2 >= up->a1  && up->a2 - up->a1 < minIntronSize)
          up->type |= gMicro ;
        if (up->type & gX)
          {
            bsAddData (Transcript, _bsRight, _Int, &(up->x1)) ; 
            bsAddData (Transcript, _bsRight, _Int, &(wp->x2)) ;
            bsPushObj (Transcript) ;
            nExons++ ; 
            if (up->type & gStolen)
              bsAddKey (Transcript, _bsRight, str2tag("Stolen_Exon")) ;
            else if (up->type & gPredicted)
              bsAddKey (Transcript, _bsRight, str2tag("Predicted_Exon")) ;
            else if (up->type & gB)
              bsAddKey (Transcript, _bsRight, str2tag("Alternative_Exon")) ;
            else
              bsAddKey (Transcript, _bsRight, _Exon) ;

            ll = wp->a2 - up->a1 + 1 ;
            mrnaLength += ll ;
            bsAddData (Transcript, _bsRight, _Text, messprintf("%d", nExons)) ;
            bsAddData (Transcript, _bsRight, _Text, "Length") ;
            bsAddData (Transcript, _bsRight, _Int, &ll) ;
            bsAddData (Transcript, _bsRight, _Text, "bp") ; 
	    up = wp ;
          }
        else if (up->type & gI)
          {
            int u1, u2 ;
            KEY foot ;
             
            {
              int x1 = up->x1 ;
              int x2 = up->x2 ;
             
              if (ii > 0 && ii < arrayMax (smrna->hits) -1 &&
                  up->a1 == (up-1)->a2 + 1 &&
                  up->a2 == (up+1)->a1 - 1)
                { x1 = (up-1)->x2 ; x2 = (up+1)->x1 ; }
              bsAddData (Transcript, _bsRight, _Int, &x1) ; 
              bsAddData (Transcript, _bsRight, _Int, &x2) ;
            }
            bsPushObj (Transcript) ;
            if (sc->isUp)
              { u1 = ta1 + sc->d1 - sc->a1 - up->a1 + 1 ; u2 = ta1  + sc->d1 - sc->a1 - up->a2 + 1 ; }
            else
              { u1 = ta1  + sc->d1 - sc->a1 + up->a1 - 1 ; u2 = ta1  + sc->d1 - sc->a1 + up->a2 - 1 ; }
            ll = up->a2 - up->a1 + 1 ;
            if (ll < 0)
              { actualFoot = foot = _v_rep ; nv_rep++ ; }
            else if (up->type & gMicro)
              {
                actualFoot = mrnaIntronBoundary (s2m->dnaD, s2m->dnaR, u1, u2) ;
                foot = 0 ;
              }
            else if (up->type & gJ)  /* intron with error */
              {
                actualFoot = mrnaIntronBoundary (s2m->dnaD, s2m->dnaR, u1, u2) ;
                foot = _Fuzzy ; 
                if (0 && actualFoot == _gt_ag) actualFoot = _Fuzzy_gt_ag ;
                if (0 && actualFoot == _gc_ag) actualFoot = _Fuzzy_gc_ag ;
              }
            else
              {
                foot = actualFoot = mrnaIntronBoundary (s2m->dnaD, s2m->dnaR, u1, u2) ;
                if (foot == _gt_ag) ngt_ag++ ;
                if (foot == _gc_ag) ngc_ag++ ;
              }
            if (foot)
              {
                if (!feet) feet = arrayCreate (30, HIT) ;
                fp = arrayp (feet, ifeet++, HIT) ;
                fp->a1 = up->a1 ; fp->a2 = up->a2 ; 
                if (bsIsTagInClass (_VTranscribed_gene, foot))
                  fp->cDNA_clone = foot ; 
                else
                  fp->cDNA_clone = _Other ;
                fp->est = actualFoot ;
              }

            nIntrons++ ;
            /* intl += a2 - a1 + 1 ; */
            if (up->type & gStolen)
              bsAddKey (Transcript, _bsRight, str2tag("Stolen_Intron")) ;
            else if (up->type & gPredicted)
              bsAddKey (Transcript, _bsRight, str2tag("Predicted_Intron")) ;  
            else if (up->type & gB)
              bsAddKey (Transcript, _bsRight, str2tag("Alternative_Intron")) ;
            else
              bsAddKey (Transcript, _bsRight, _Intron) ;

             if (foot != _Fuzzy)
              bsAddData (Transcript, _bsRight, _Text, name(actualFoot)) ;
            else
              bsAddData (Transcript, _bsRight, _Text, messprintf ("Fuzzy_%s", name(actualFoot))) ;
           
            /*
              frame = (mrnaLength % 3) + 1 ;
              bsAddData (Transcript, _bsRight, _Int, &frame) ;
            */
            ll = up->a2 - up->a1 + 1 ;
            bsAddData (Transcript, _bsRight, _Text, "Length") ;
            bsAddData (Transcript, _bsRight, _Int, &ll) ;
            bsAddData (Transcript, _bsRight, _Text, "bp") ;
          }
        else if (up->type & gGap)
          {
            int xx ;
            xx =  (up-1)->x2 + 1 ;
            bsAddData (Transcript, _bsRight, _Int, &xx) ; 
            xx =  (up+1)->x1 - 1 ;
            bsAddData (Transcript, _bsRight, _Int, &xx) ;
            bsPushObj (Transcript) ;
            ll = (up+1)->x1 - (up-1)->x2 - 1 ; 
            if (ll < 3)
              {
                if (up->a2 - up->a1 + 1 < 12)
                  ll = up->a2 - up->a1 + 1 ;
                else
                  ll = 3 ;
              }
            gapl += ll ;
            mrnaLength += ll ;
            if (up->type & gOpenGap)
              bsAddKey (Transcript, _bsRight, str2tag("ORF_Gap")) ;
            else
              bsAddKey (Transcript, _bsRight, _Gap) ;
            bsAddData (Transcript, _bsRight, _Text, "Unsequenced") ;
            bsAddData (Transcript, _bsRight, _Text, "Length") ;
            bsAddData (Transcript, _bsRight, _Int, &ll) ;
            bsAddData (Transcript, _bsRight, _Text, "bp") ;
            {  /* totosp */
              OBJ Tg = bsUpdate (sc->gene) ;
              int u1, u2 ;
              BSMARK tgmark = 0 ;
              if (Tg)
                {
                  if (bsGetData (Tg, _Splicing, _Int, &u1))
                    do
                      {
                        tgmark = bsMark (Tg, tgmark) ;
                        if (u1 == up->a1 + smrna->a1 - 1)
                          {
                            if (bsGetData (Tg, _bsRight, _Int, &u2))
                              do
                                {
                                  if (u2 == up->a2 + smrna->a1 - 1)
                                    {
                                      if (bsPushObj (Tg) &&
                                          bsFindTag (Tg, _Gap))
                                        {
                                          bsAddData (Tg, _bsRight, _Text, "Unsequenced") ;
                                          bsAddData (Tg, _bsRight, _Text, "Length") ;
                                          bsAddData (Tg, _bsRight, _Int, &ll) ;
                                          bsAddData (Tg, _bsRight, _Text, "bp") ;
                                          goto tgDone ;
                                        }
                                    }
                                } while (bsGetData (Tg, _bsDown, _Int, &u2)) ;
                          }
                        bsGoto (Tg, tgmark) ;
                      } while (Tg && bsGetData (Tg, _bsDown, _Int, &u1)) ;
                tgDone:
                  bsSave (Tg) ;
               }
              bsMarkFree (tgmark) ;
            }
          }
        bsGoto (Transcript, 0) ;
      }

    if (mrnaLength > 0)
      bsAddData (Transcript, str2tag ("Total_length"), _Int, &mrnaLength) ;
    if (gapl > 0)
      bsAddData (Transcript, str2tag ("Gap_length"), _Int, &gapl) ;
    if (ngt_ag)
      bsAddData (Transcript,  _gt_ag, _Int, &ngt_ag) ;
    if (ngc_ag)
      bsAddData (Transcript,  _gc_ag, _Int, &ngc_ag) ;
    if (nfuzzy)
      bsAddData (Transcript,  _Fuzzy, _Int, &nfuzzy) ;
    if (nv_rep)
      bsAddData (Transcript,  _v_rep, _Int, &nv_rep) ;
    if (feet)
      {
        int i1, i2, i3, iMax = arrayMax(feet) ;
        HIT *fp2 ;

        cDNASwapA (feet) ;
        arraySort (feet, cDNAOrderByA1) ;
         
        for (i1 = 0 ; i1 < iMax ; i1++)
          {
            fp = arrp (feet, i1, HIT) ;
            if (fp->cDNA_clone == _Other &&
                ! fp->est) /* this is a Fuzzy */
              continue ;
            i3 = 0 ;
            for (i2 = i1, fp2 = fp ;  fp2->cDNA_clone == fp->cDNA_clone && fp2->est == fp->est && i2 < iMax ; i2++, fp2++) 
              i3++ ;  /* get the number, first thing to report */
            for (i2 = i1, fp2 = fp ; fp2->cDNA_clone == fp->cDNA_clone && fp2->est == fp->est && i2 < iMax ; i2++, fp2++) 
              {
                fp2->x1 = i3 ; /* for later reference */
                bsAddTag (Transcript, fp2->cDNA_clone) ;
                if (fp2->cDNA_clone == _Other)
                  bsAddData (Transcript, _bsRight, _Text, name(fp2->est)) ;
                bsAddData (Transcript, _bsRight, _Int, &i3) ;
                if (fp2->cDNA_clone != _gt_ag && fp2->cDNA_clone != _gc_ag)
                  {
                    bsAddData (Transcript, _bsRight, _Int, &(fp2->a1)) ;
                    bsAddData (Transcript, _bsRight, _Int, &(fp2->a2)) ;

                    /* search the justifying clone in the gene */
                    {
                      OBJ Tg = 0 ;
                      BSunit *uu ;
                      int iu ;
                      KEY direction = _bsRight ;
                      
                      if ((Tg = bsCreate (sc->gene)))
                        {
                          BSMARK trmark = 0 ;
                          if (bsGetArray (Tg, str2tag ("Intron_boundaries"),  units, 6))
                            for (iu = 0 ; iu < arrayMax (units) ; iu+= 6)
                              {
                                uu = arrp (units, iu, BSunit) ;
                                if (uu[0].k == fp2->cDNA_clone)
                                  {
                                    if (fp2->cDNA_clone == _Other &&
                                        (uu[1].s && !strcmp(name(fp2->est), uu[1].s)) &&
                                        uu[3].k == fp2->a1 + smrna->a1 - 1 &&
                                        uu[4].k == fp2->a2 + smrna->a1 - 1 &&
                                        uu[5].k)
                                      {
                                        bsAddKey (Transcript,  direction, uu[5].k) ; 
                                        direction = _bsDown ;
                                      }
                                    else if (fp2->cDNA_clone != _Other &&
                                             uu[2].k == fp2->a1 + smrna->a1 - 1 &&
                                             uu[3].k == fp2->a2 + smrna->a1 - 1 &&
                                             uu[4].k)
                                      {
                                        bsAddKey (Transcript, direction, uu[4].k) ; 
                                        if (fp->cDNA_clone == _Fuzzy)
                                          {
                                            trmark = bsMark (Transcript, trmark) ;
                                            bsAddData (Transcript, _bsRight, _Text, name(fp->est)) ;
                                            bsGoto (Transcript, trmark) ;
                                          }
                                        direction = _bsDown ;
                                      }
                                  }                            
                              }
                          bsMarkFree (trmark) ;
                          bsDestroy (Tg) ;
                        }
                    }
                  }
              }
            i1 = i2 - 1 ;
          }
      }
  }
  
  /*** Constructed_from ***/
  if (smrna->dnas && smrna->bestDna < arrayMax (smrna->dnas))
    {
    Array sHits ;
    ORFS *orf ;
    int jj, x1, x2, a1, a2, b1, b2, sMax, eMax = arrayMax (estHits) ;
    HIT *up, *vp, *wp ;
    KEYSET ks = keySetCreate () ;
    BOOL estUp, complete3 = FALSE, complete5 = FALSE ;
    int naggregated = 0 ;
    KEY _Composite = str2tag ("Composite") ;
    KEY est ;
    KEY _Not_real_5prime = str2tag ("Not_real_5prime") ;
    OBJ Est = 0, Clone = 0 ;

    orf = arrp (smrna->dnas, smrna->bestDna, ORFS) ; /* ici ... */
    sHits = orf->hits ;
    sMax = arrayMax (sHits) ;
    cloneHits = arrayCreate (200, HIT) ;

    if (bsFindTag (Transcript, _Constructed_from))
      bsRemove (Transcript) ;

    bsAddTag (Transcript, str2tag("Completeness")) ;
    for (ii = jj = 0 ; ii < eMax; ii++)
      {
        up = arrp(estHits, ii, HIT) ;
        est = up->est ;
        if (!class(up->est)) /* ghosts */
          continue ; 
        if (! (up->type & gX))
          continue ;
        if (keySetFind (ks, est, 0))
          continue ;
        if (up->a1 < smrna->a1 - 12 ||
            up->a2 > smrna->a2 + 12)
          continue ;
        keySetInsert (ks, est) ;
        
        /* was here to count aggregated in single mrna 
        if (!keySetFind (smrna->clones, up->cDNA_clone, 0))
          continue ;
        */        
        x1 = up->x1 ; x2 = up->x2 ; a1 = up->a1 ; a2 = up->a2 ;
        estUp = x1 < x2 ? FALSE : TRUE ;
                
        if (!estUp &&            /* 5' read */
	    !keyFindTag (est, _Composite) &&
            !keySetFind (ks, up->cDNA_clone, 0) &&
            !keyFindTag (up->cDNA_clone, _Not_real_5prime) &&
            !keyFindTag (up->cDNA_clone, _RT_PCR) &&
            /*   smrna->a1 < 200 && suppressed 2007_02_07 */
            up->a1 > smrna->a1 - 20 && up->a1 < smrna->a1 + 20 &&        /* cluster at top of mrna */
            up->x1 == up->clipTop /* avoid genomic gap */
            )
          {
            naggregated++ ;
            keySetInsert (ks, up->cDNA_clone) ;
          }

        if (0 && /* on peut pas lire les flags dans les composites */
	    !estUp &&            /* composite 5' read */
	    keyFindTag (est, _Composite) &&
            !keySetFind (ks, up->cDNA_clone, 0) &&
            !keyFindTag (up->cDNA_clone, _Not_real_5prime) &&
            !keyFindTag (up->cDNA_clone, _RT_PCR) &&
            /*   smrna->a1 < 200 && suppressed 2007_02_07 */
            up->a1 > smrna->a1 - 20 && up->a1 < smrna->a1 + 20 &&        /* cluster at top of mrna */
            up->x1 == up->clipTop /* avoid genomic gap */ 
            )
          {
	    OBJ Est = bsCreate (est) ;
	    int n = 0 ;
	    bsGetData (Est, _Composite, _Int, &n) ;
	    bsDestroy (Est) ;
            naggregated += n ; 
            keySetInsert (ks, up->cDNA_clone) ;
          }

        /* now here to count aggregated in whole gene
         *  but limit constructed_from to single mrna
         */
        if (0 && !keySetFind (smrna->clones, up->cDNA_clone, 0))
          continue ;

        for (wp = 0, jj = 0 ; jj < arrayMax(cloneHits) ; jj++)
          {
            wp = arrp (cloneHits, jj, HIT) ;
            if (wp->cDNA_clone == up->cDNA_clone)
              break ;
            else
              wp = 0 ;
          }
        if (!wp)
          wp = arrayp (cloneHits, arrayMax(cloneHits), HIT) ;
        wp->nerr += up->nerr ;
        wp->nerrAll += up->nerrAll ;
        wp->cDNA_clone = up->cDNA_clone ;
        wp->clipTop |= (estUp ? 2 : 1) ;
        for (vp = up + 1, jj = ii + 1 ; jj < eMax ; jj++, vp++)
          {
            if (! (vp->type & gX))
              continue ; 
            if (wp->cDNA_clone != vp->cDNA_clone)
              continue ;
            if (est != vp->est)
              continue ;
            if (vp->a2 > smrna->a2 + 12)
              continue ;
            wp->nerr += vp->nerr ;
            wp->nerrAll += vp->nerrAll ;
            if (estUp)
              {
                if (vp->x1 > x1) x1 = vp->x1 ;
                if (vp->x2 < x2) x2 = vp->x2 ;
              }
            else
              {
                if (vp->x1 < x1) x1 = vp->x1 ;
                if (vp->x2 > x2) x2 = vp->x2 ;
              }
            a2 = vp->a2 ;
          }
        /* we are in unspliced gmrna coordinates */
        a1 = a1 - smrna->a1 + 1 ; a2 = a2 - smrna->a1 + 1 ;
        /* we are in unspliced smrna coordinates */
        b1 = a1 ; b2 = a2 ;
        for (jj = 0, vp = arrp (sHits, 0, HIT) ; jj < sMax ; vp++, jj++)
          {
            if ((vp->type & gX) &&
                (jj == sMax - 1 || vp->a2 >= a1))
              {
                /* mieg jan 2004, see AL526305 in ZH18 human */
                /* because the last exon of this 3' read was merged in the middle of a big stolen exon */
                if (up->a2 - smrna->a1 + 1 < vp->a2 - 50)
                  /* try to reach next exon on the est */
                  {
                    int jj1 ;
                    HIT *vp1 ;

                    for (vp1 = up + 1, jj1 = ii + 1 ; jj1 < eMax ; jj1++, vp1++)
                      if ((vp1->type & gX) &&
                          vp1->est == est &&
                          vp1->a1 - smrna->a1 + 1 < vp->a2)
                        {
                          b1 = vp1->a1 - smrna->a1 + 1 ;
                          x1 = vp1->x1 ;
                        }
                  }                
                b1 = b1 - vp->a1 + vp->x1 ; 
                break ;
              }
          }
        for (/* same place */ ; jj < sMax  ; vp++, jj++)
          { 
            if (jj == sMax - 1 ||  vp->a2 >= a2)
              { 
                if (vp->type & gX)  /* exon */
                  b2 = b2 - vp->a1 + vp->x1 ; 
                else if (jj > 0) /* gap or intron, take exon above it */
                  {
                    while (jj > 0 && !(vp->type & gX)) { jj-- ; vp-- ; }
                    b2 = b2 - vp->a2 + vp->x2 ;
                  }
                if (b2 > vp->x2 + 300) 
                  b2 = b1 + (x1 < x2 ? x2 - x1 : x1 - x2) ; /* workaround, we have a bug here see W17197, S_20_5 */
                break ; 
              }
          }
        wp->est = up->est ;

        if (!(wp->x1 + wp->x2)) /* initialise */
          {  wp->a1 = b1 ; wp->a2 = b2 ; }
        if (x1 < x2) wp->x1++ ; /* nb 5p read */
        else   wp->x2++ ;       /* nb 3p read */
        if (b1 < wp->a1) wp->a1 = b1  ;
        if (b2 > wp->a2) wp->a2 = b2  ;

        bsAddData (Transcript, _Constructed_from, _Int, &b1) ;
        bsAddData (Transcript, _bsRight, _Int, &b2) ;
        bsAddKey (Transcript, _bsRight, up->est) ;
        bsAddData (Transcript, _bsRight, _Int, &x1) ;
        bsAddData (Transcript, _bsRight, _Int, &x2) ;

        Est = bsCreate(up->est) ;
        Clone = bsCreate(up->cDNA_clone) ;
        if (Est && Clone)
          {
            int x = 0 ;
	    int clip1 = 1, clip2 = 0 ;
	    if (!keyFindTag (up->est, _DNA)) /* clip2 = dna length */
	      { 
                BSunit *uu ;
                        
                units = arrayReCreate (units, 8, BSunit) ;
		if (bsGetArray (Est,  _DNA, units, 2))
		  {
		    uu = arrp (units, 0, BSunit) ;
		    clip2 = uu[1].i ; 
		  }
	      }

	    if (!keyFindTag (up->est, _Vector_Clipping)) /* clip1,clip2 = vector clipping */
	      { 
                BSunit *uu ;
                        
                units = arrayReCreate (units, 8, BSunit) ;
		if (bsGetArray (Est,  _Vector_Clipping, units, 2))
		  {
		    uu = arrp (units, 0, BSunit) ;
		    clip1 = uu[0].i ; 
		    clip2 = uu[1].i ; 
		  }
	      }
	    else if (!keyFindTag (up->est, _Clipping)) /* clip1,clip2 = vector clipping */
	      { 
                BSunit *uu ;
                        
                units = arrayReCreate (units, 8, BSunit) ;
		if (bsGetArray (Est, _Clipping, units, 2))
		  {
		    uu = arrp (units, 0, BSunit) ;
		    clip1 = uu[0].i ; 
		    clip2 = uu[1].i ; 
		  }
	      }

            if (!keyFindTag (up->cDNA_clone, str2tag("Internal_capping")))
              {
                int i, pure = 0 ;
                BSunit *uu ;
		
                units = arrayReCreate (units, 800, BSunit) ;
                if (bsGetArray (Est,  str2tag("Transpliced_to"), units, 2))
                  {
                    int dx = isComposite ? 1 : 60 ;
                    for (i = 0 ; i < arrayMax(units) ; i += 2)
                      {
                        uu = arrp (units, i, BSunit) ;
                        if ((uu[1].i > x1 - 3 && uu[1].i <= x1 + 3) && /* wobble helps in graphic interface */
                            up->a1 > smrna->a1 - dx && up->a1 < smrna->a1 + dx &&        /* cluster at top of mrna */
                            !strncasecmp ("SL", name(uu[0].k), 2))
                          {
                            int s = 0 ;
			    bsGetData (Est, str2tag("Composite"), _Int, &s) ;
			    bsAddKey (Transcript, str2tag("Transpliced_to"), uu[0].k) ;
                            bsAddKey (Transcript, _bsRight, up->cDNA_clone) ;
			    if (s)
			      bsAddData (Transcript, _bsRight, _Int, &s) ;
			    
                            if (!strcasecmp ("SL0", name(uu[0].k))) /* av=61.7 sd=19 susuki EMBO report 2001 */
                              pure |=  1 ;
                            else if (!strcasecmp ("SL1", name(uu[0].k)))
                              pure |=  2 ;
                            else if (!strncasecmp ("SL", name(uu[0].k), 2))
                              pure |=  4 ;
                          }
                      }
                    if (pure)
                      complete5 = TRUE ;
                    if (pure >= 0)
                      {
                        KEY library = 0, stage = 0, tag1 = 0, tag2 = 0, tag4 = 0;
                        OBJ Library = 0 ;
                        
                        library = keyGetKey (up->cDNA_clone, _Library) ;
                        if (library && (Library = bsCreate (library)))
                          {
                            bsGetKeyTags (Library, str2tag("Stage"), &stage) ;
                            bsDestroy (Library) ;
                          }
                        
                        if (pure & 1)
                          {
                            if (!strcasecmp ("Embryo", name(stage)))
                              tag1 = str2tag ("SL0_e") ;
                            else if (!strcasecmp ("Dauer", name(stage)))
                              tag1 = str2tag ("SL0_d") ;
                            else if (!strcasecmp ("L1", name(stage)))
                              tag1 = str2tag ("SL0_L1") ;
                            else if (!strcasecmp ("L2", name(stage)))
                              tag1 = str2tag ("SL0_L2") ;
                            else if (!strcasecmp ("L3", name(stage)))
                              tag1 = str2tag ("SL0_L3") ;
                            else if (!strcasecmp ("L4", name(stage)))
                              tag1 = str2tag ("SL0_L4") ;
                            else 
                              tag1 = str2tag ("SL0_m") ;
                           }
                        if (pure & 2)
                          {
                            if (!strcasecmp ("Embryo", name(stage)))
                              tag2 = str2tag ("SL1_e") ;
                            else if (!strcasecmp ("Dauer", name(stage)))
                              tag2 = str2tag ("SL1_d") ;
                            else if (!strcasecmp ("L1", name(stage)))
                              tag2 = str2tag ("SL1_L1") ;
                            else if (!strcasecmp ("L2", name(stage)))
                              tag2 = str2tag ("SL1_L2") ;
                            else if (!strcasecmp ("L3", name(stage)))
                              tag2 = str2tag ("SL1_L3") ;
                            else if (!strcasecmp ("L4", name(stage)))
                              tag2 = str2tag ("SL1_L4") ;
                            else 
                              tag2 = str2tag ("SL1_m") ;
                          }
                        if (pure & 4)
                          {
                            if (!strcasecmp ("Embryo", name(stage)))
                              tag4 = str2tag ("SL2s_e") ;
                            else if (!strcasecmp ("Dauer", name(stage)))
                              tag4 = str2tag ("SL2_d") ;
                            else if (!strcasecmp ("L1", name(stage)))
                              tag4 = str2tag ("SL2s_L1") ;
                            else if (!strcasecmp ("L2", name(stage)))
                              tag4 = str2tag ("SL2s_L2") ;
                            else if (!strcasecmp ("L3", name(stage)))
                              tag4 = str2tag ("SL2s_L3") ;
                            else if (!strcasecmp ("L4", name(stage)))
                              tag4 = str2tag ("SL2s_L4") ;
                            else 
                              tag4 = str2tag ("SL2s_m") ;
                          }
                        if (tag1 && bsIsTagInObj (Transcript, tr, tag1))
                          {
                            int n = 0 ;
                            bsAddTag (Transcript, tag1) ;
                            bsGetData (Transcript, tag1, _Int, &n) ;
                            n++ ;
                            bsAddData (Transcript, tag1, _Int, &n) ;
                          }
                        if (tag2 && bsIsTagInObj (Transcript, tr, tag2))
                          {
                            int n = 0 ;
                            bsAddTag (Transcript, tag2) ;
                            bsGetData (Transcript, tag2, _Int, &n) ;
                            n++ ;
                            bsAddData (Transcript, tag2, _Int, &n) ;
                          }
                        if (tag4 && bsIsTagInObj (Transcript, tr, tag4))
                          {
                            int n = 0 ;
                            bsAddTag (Transcript, tag4) ;
                            bsGetData (Transcript, tag4, _Int, &n) ;
                            n++ ;
                            bsAddData (Transcript, tag4, _Int, &n) ;
                          }
                      }
                  }
              }
                
            if (bsFindTag (Est, str2tag ("Complete_mRNA")))
              {
		int u1 = x1 < x2 ? x1 : x2 ;
		int u2 = x1 < x2 ? x2 : x1 ;
                if (a1 < 10 && u1 < clip1 + 10)
                  {
                    bsAddKey (Transcript, str2tag ("Submitted_as_5p_complete"), up->cDNA_clone) ;
                    complete5 = TRUE ;
                  }
                if (a2 > mrnaLength - 10 && u2 > clip2 - 10)
                  {
                    bsAddKey (Transcript, str2tag ("Submitted_as_3p_complete"), up->cDNA_clone) ;
                    complete3 = TRUE ;
                  }
              }
            if (bsFindTag (Est, str2tag ("Real_5prime")))
              {
		int u1 = x1 < x2 ? x1 : x2 ;
                if (a1 < 10 && u1 < clip1 + 10)
                  {
		    if (bsFindTag (Est, str2tag ("Composite")))
		      {
			int k = 0 ;
			if (bsGetData (Est, str2tag ("Composite"), _Int, &k))
			  if (k > 0)
			    {
			      bsAddData (Transcript, str2tag ("Aggregated_5p_clones"), _Int, &k) ;
			    }
		      }
		    else
		      bsAddKey (Transcript, str2tag ("Submitted_as_5p_complete"), up->cDNA_clone) ;
                    complete5 = TRUE ;
                  }
              }
            if (bsFindTag (Est, str2tag ("Real_3prime")))
              {
		int u2 = x1 < x2 ? x2 : x1 ;
                if (a2 > mrnaLength - 10 && u2 > clip2 - 10)
                  {
		    if (bsFindTag (Est, str2tag ("Composite")))
		      {
			int k = 0 ;
			if (bsGetData (Est, str2tag ("Composite"), _Int, &k))
			  if (k > 0)
			    {
			      bsAddData (Transcript, str2tag ("Aggregated_3p_clones"), _Int, &k) ;
			    }
		      }
		    else
		      bsAddKey (Transcript, str2tag ("Submitted_as_3p_complete"), up->cDNA_clone) ;
                    complete3 = TRUE ;
                  }
              }

            if ( (vp->type & gX) && vp->a1 + smrna->a1 - 1 <= up->a2 && vp->a2 + smrna->a1 - 1 >= up->a2 && 
                !bsFindTag (Clone, str2tag("Internal_priming")) &&
                bsGetData (Est, _PolyA_after_base, _Int, &x))
              {
                if (/* b2 > mrnaLength + stolen - 20 && stupid, preventslater flagging of internal_priming */
                    x > x2 - 5 && x < x2 + 5)
                  {
		    int s = 0 ;
		    bsGetData (Est, str2tag("Composite"), _Int, &s) ;
                    bsAddData (Transcript, str2tag ("PolyA_found"), _Int, &a2) ;
                    bsAddKey (Transcript, _bsRight, up->cDNA_clone) ;
		    if (s)
		      bsAddData (Transcript, _bsRight, _Int, &s) ;
                    complete3 = TRUE ;
                  }
		/* int aaa2 = 1000000 + a2 ;  debugging */
                   if (a2 == 99995342) invokeDebugger () ;
                bsAddData (Transcript, str2tag ("PolyA_seen"), _Int, &a2) ;
                bsAddKey (Transcript, _bsRight, up->cDNA_clone) ;
              }
            if ( (vp->type & gX) && vp->a1 + smrna->a1 - 1 <= up->a2 && vp->a2 + smrna->a1 - 1 >= up->a2 && 
		bsFindTag (Est, _Reverse) &&
                /* b2 > mrnaLength + stolen - 20 && stupid, prevents later flagging of internal_priming 
                 * because the ORF will not reach the polyA if we have imposed completion
                 */
                x2 < 30 &&
                !bsFindTag (Clone, str2tag("Internal_priming")) &&
                !bsFindTag (Clone, _Mosaic) &&
                (
                 bsFindTag (Est, _PolyA_after_base) || bsFindTag (Clone, str2tag("Primed_on_polyA"))
                 ))
              {
                complete3 = TRUE ;
                bsAddTag (Transcript, str2tag ("PolyA_primed_clone")) ;
                if (bsFindTag (Est, _PolyA_after_base))
                  bsAddData (Transcript, str2tag ("PolyA_found"), _Int, &a2) ;
                else
                  bsAddData (Transcript, str2tag ("PolyA_primed_clone_seen"), _Int, &a2) ;
                bsAddKey (Transcript, _bsRight, up->cDNA_clone) ;
              }
            bsDestroy (Est) ;
            bsDestroy (Clone) ;
          }
      }

    keySetDestroy (ks) ;
    if (naggregated > 2 && 
        100*naggregated > 5*keySetMax(sc->sh.clones) /* 2007_02_07, min 4 percent of all gene clones */
        )
      { bsAddTag (Transcript, str2tag ("Aggregated_5p_clones")) ; complete5 = TRUE ; }
    if (complete5) bsAddTag (Transcript, str2tag ("Found5p")) ;
    if (complete3) bsAddTag (Transcript, str2tag ("Found3p")) ;
    if (complete3 && complete5) bsAddTag (Transcript, str2tag ("Complete")) ;
  }
  
  /*** Child->Products ***/
  if (smrna->orfs && arrayMax(smrna->orfs))
    {
      int oMax = arrayMax(smrna->orfs) ;
      int iProd = 0 ;
      
      /* we need to reorder by upstop to correctly evaluate the amount of stolen and givenBack bases
	 when adjusting from genomic to mrna coordinates in saveMrna 
	 but this must be done AFTER ordering the mRNAs
      */
      
      /* mars 25 2005, i change strategy */
      /* first round to measure givenBack */
      arraySort (smrna->orfs, orfOrderByUpStop) ;
      for (ii = 0 ; ii < 1 ; ii++)
	{
	  BOOL isNh2Complete = FALSE ;
	  ORFS * orf = arrp (smrna->orfs, ii, ORFS) ;
	  int jj ;
	  HIT *up ;
	  KEY lib ;
	  char *cp ;
	  int beDisHonest = 0; /* how much do we steal */
	  
	  if (!orf->upGap &&
	      (
	       orf->isSL || 
	       bsFindTag (Transcript, str2tag("Transpliced_to")) ||
	       (
		(bsFindTag (Transcript, str2tag("Submitted_as_5p_complete")) ||
		 bsFindTag (Transcript, str2tag("Aggregated_5p_clones"))
		 ) &&
		orf->met < orf->start + MINI_LEU2MET/2 &&
		(orf->met <= 0 || (orf->met > 0 && orf->start < 150 + orf->stolen))
		)
	       )
	      )
	    isNh2Complete = TRUE ;
	  
	  if (isNh2Complete && orf->start != orf->met && orf->start != orf->leu)
	    {
	      int s = -1, dx = 0 ;
	      if (orf->met > 0 && (orf->leu <= 0 || orf->leu > orf->met - MINI_LEU2MET))
		s = orf->met ;
	      else if (orf->leu > 0)
		s = orf->leu ;
	      if (s > 0)
		dx = s - orf->start ;
	      if (dx > 0)
		{        orf->start = s ; orf->cds -= dx ;}
	    }
	  
	  if (orf->met > 0 && orf->met <= orf->stolen && orf->isBest)
	    {
	      int dStolen = orf->stolen - orf->met + 1 ;
	      
	      stolen = orf->stolen ;
	      givenBack = 0 ;
	      
	      if (getPleaseStealUpStream () )
		{
		  for (jj = 0, up = arrp (estHits, 0, HIT) ; !beDisHonest && jj < arrayMax(estHits) ; jj++, up++)
		    {
		      if (smrna->a1 + dStolen < up->a1 || up->x1 > up->x2 || 
			  up->x1 > up->clipTop + 3 ||
			  ! keySetFind (smrna->clones, up->cDNA_clone, 0))
			continue ;
		      lib = keyGetKey (up->cDNA_clone, _Library) ;
		      if (lib)
			{
			  cp = name(lib) ;
			  if (*cp++ == 'y' && *cp++ == 'k' && 
			      (*cp == 'e' || *cp == 'm'))
			    beDisHonest = dStolen ;
			}
		      if (! beDisHonest && dStolen == 1 && smrna->a1 == up->a1 &&
			  keyFindTag (up->cDNA_clone, str2tag ("WORF")))
			beDisHonest = dStolen ;
		    }
		}
	      if (! beDisHonest)
		{
		  orf->nOpen += orf->start - orf->frame ;
		  orf->cds = orf->nOpen ;
		  orf->upStop = orf->met = orf->frame - 3  ; /* so negative in frame, and i will give back */
		  orf->start = orf->frame ;
		}
	    }
	  
	  if (orf->met > 0 && orf->met <= orf->stolen && beDisHonest)
	    {
	      stolen = stolen - orf->met + 1 ;
	      if (stolen < 0) stolen = 0 ;
	      givenBack = orf->stolen - stolen ;  
	    }
	  else
	    {
	      stolen = 0 ;
	      givenBack = orf->stolen ; 
	    }
	}
      
      /* true export and product naming */
      arraySort (smrna->orfs, orfOrderByCds) ;
      prod_array = arrayCreate (12, HIT) ;
      for (ii = iProd = 0 ; ii < oMax ; ii++)
	{
	  ORFS * orf = arrp (smrna->orfs, ii, ORFS) ;
	  KEY product = 0 ;
	  OBJ Product = 0 ;
	  int jj, x, x1, x2, u1, u2 ;
	  HIT *up ;
	  BOOL isNh2Complete = FALSE ;
	  Array pep = 0 ;
	  
	  if (!orf->upGap &&
	      (
	       orf->isSL || 
	       bsFindTag (Transcript, str2tag("Transpliced_to")) ||
	       bsFindTag (Transcript, str2tag("Found5p")) ||
	       (
		(bsFindTag (Transcript, str2tag("Submitted_as_5p_complete")) ||
		 bsFindTag (Transcript, str2tag("Aggregated_5p_clones"))
		 ) && 
		orf->met < orf->start + MINI_LEU2MET/2 &&
		orf->met > 0 && orf->start < 150 + orf->stolen))
	      )
	    isNh2Complete = TRUE ;
	  
	  if (0 && /* 2005_09_15, idiot on veut annoter, CF 5E914 */ 
	      isNh2Complete && iProd > 0 &&
	      orf->met <= orf->stolen)
	    continue ; 
	  
	  /* mars 25 2005, i change the naming strategy */
	  if (iProd > 0)
	    lexaddkey (messprintf ("%s%d", name(tr), iProd), &product, _VmProduct) ;
	  else
	    lexaddkey (messprintf ("%s", name(tr)), &product, _VmProduct) ;
	  iProd++ ;
	  
	  prod_p = arrayp (prod_array, ii, HIT) ;
	  prod_p->gene = product ;
	  
	  if (orf->met > 0 && orf->met <= orf->stolen)
	    {
	      /* isNh2Complete = TRUE ; since i stole up to the Met, commented out nov 16 2004 */
	      int dx = orf->met - orf->start ;
	      if (dx > 0)
		{ orf->start += dx ; orf->nOpen -= dx ; orf->cds -= dx ; }
	    }
	  else if (orf->start < givenBack)
	    {
	      int dx =  givenBack ;
	      orf->start += dx ; 
	      if (orf->nOpen == orf->cds) orf->cds -= dx ;
	      orf->nOpen -= dx ;
	    }
	  
	  if (isNh2Complete && orf->start < orf->met)
	    orf->start = orf->met ;
	  
	  x1 = orf->start ;              
	  x2 = orf->max ; 
	  
	  if (orf->downStop > 0) x2 += 3 ; /* export the stop */
	  bsAddKey (Transcript, _Product, product) ; 
	  prod_p->x1 = x1 ; /* needed to flag  complete_cds_clones */
	  prod_p->x2 = x2 ; /* needed to flag  complete_cds_clones */
	  x1 -= givenBack ; x2 -= givenBack ;
	  bsAddData (Transcript, _bsRight, _Int, &x1) ;  
	  bsAddData (Transcript, _bsRight, _Int, &x2) ;
	  if ((Product = bsUpdate(product)))
	    {
	      KEY M_mProduct ;
	      
	      Array mdna = arr (smrna->dnas, orf->iDna, ORFS).dna ;
	      KEY geneticCode = 0 ;
	      char *translationTable = pepGetTranslationTable (s2m->cosmid1 ? s2m->cosmid1 : s2m->cosmid, &geneticCode) ;
	      
	      if (geneticCode)
		bsAddKey (Transcript, _Genetic_code, geneticCode) ;
	      
	      if (orf->isUorf)
		bsAddTag (Product, str2tag("uORF_candidate")) ;
	      bsAddTag (Product, _CDS) ;
	      lexaddkey ("Product", &M_mProduct, _VMethod) ;
	      bsAddKey (Product, _Method, M_mProduct) ;
	      
	      x = (180 + orf->frame - 1 - givenBack) % 3 + 1 ;
	      bsAddData (Product, _Frame, _Int, &x) ;
	      
	      bsAddTag (Product, str2tag("Completeness")) ;
	      if (isNh2Complete) 
		{
		  bsAddTag (Product, str2tag("NH2_Complete")) ;
		  if (bsFindTag (Transcript, str2tag("Found5p")))
		    bsAddTag (Product, str2tag("mRNA_5p_complete")) ;
		}
	      if (orf->met > 0 && isNh2Complete)
		bsAddTag (Product, str2tag("at_position_1 ")) ; 
	      x = orf->nOpen - orf->nX ;
	      if (x > 0)
		{
		  int xx ;
		  
		  if (orf->isBest)
		    bsAddTag (Product, str2tag ("Best_product")) ;
		  
		  bsAddData (Product, str2tag("Open_length"), _Int, &x) ; 
		  bsAddData (Product, _bsRight, _Text, "bp") ;
		  xx = x/3 ;
		  bsAddData (Product, _bsRight, _Int, &xx) ;
		  bsAddData (Product, _bsRight, _Text, "aa") ;  
		  
		  if (!bsGetData (Transcript, str2tag ("Longest_ORF"), _Int, &xx) ||
		      xx < x)
		    {
		      bsAddData (Transcript, str2tag("Longest_ORF"), _Int, &x) ;
		      bsAddData (Transcript, _bsRight, _Text, "bp") ;
		      x /= 3 ;
		      bsAddData (Transcript, _bsRight, _Int, &x) ;
		      bsAddData (Transcript, _bsRight, _Text, "aa") ;
		    }                
		}
	      x = orf->nX ;
	      if (x > 0)
		{
		  int xx ;
		  
		  bsAddData (Product, str2tag("Coding_gap"), _Int, &x) ; 
		  bsAddData (Product, _bsRight, _Text, "bp") ;
		  xx = x/3 ;
		  bsAddData (Product, _bsRight, _Int, &xx) ;
		  bsAddData (Product, _bsRight, _Text, "X") ;  
		}
	      
	      if (orf->cds > 0)
		{
		  x = orf->cds - orf->nX ;
		  if (x > 0 || isNh2Complete)
		    {
		      bsAddData (Product, str2tag("Coding_length"), _Int, &x) ;
		      bsAddData (Product, _bsRight, _Text, "bp") ;
		      x /= 3 ;
		      bsAddData (Product, _bsRight, _Int, &x) ;
		      bsAddData (Product, _bsRight, _Text, "aa") ;
		      x = x1 - 1 ;
		      if (x > 0 && isNh2Complete)
			bsAddData (Product, _Length_5prime_UTR, _Int, &x) ;
		      if (orf->isBest)
			{        
			  x = orf->cds  - orf->nX ;
			  bsAddData (Transcript, str2tag("Longest_CDS"), _Int, &x) ;
			  bsAddData (Transcript, _bsRight, _Text, "bp") ;
			  x /= 3 ;
			  bsAddData (Transcript, _bsRight, _Int, &x) ;
			  bsAddData (Transcript, _bsRight, _Text, "aa") ;
			  x = x1 - 1 ; 
			  if (x > 0 && isNh2Complete)
			    bsAddData (Transcript, _Length_5prime_UTR, _Int, &x) ;
			  bsRemoveTag (Transcript, _Length_3prime_UTR) ;
			}                      
		    }
		}
	      
	      x = arrayMax(mdna) - givenBack - x2 ;
	      bsRemoveTag (Product, _Length_3prime_UTR) ;
	      if (orf->downStop > 0)
		{ 
		  KEYSET alt3 =  mrnaAlternative3primeLengths (s2m, sc, estHits, smrnas, smrna) ;
		  int ia3, xx ;
		  
		  x = arrayMax(mdna) - orf->max - 3 ;
		  if (x > 2)
		    bsAddData (Product, _Length_3prime_UTR, _Int, &x) ;
		  if (orf->isBest)
		    bsAddData (Transcript, _Length_3prime_UTR, _Int, &x) ;
		  
		  if (alt3)
		    for (ia3 = 0 ; ia3 < keySetMax(alt3) ; ia3++)
		      {
			xx = x - keySet (alt3, ia3) ;
			bsAddData (Product, _Length_3prime_UTR, _Int, &xx) ;
		      }
		  keySetDestroy (alt3) ;
		}
	      
	      bsAddKey (Product, _From_gene, sc->gene) ;
	      
	      /* introns in and out */
	      {
		int i ;
		int nIntronsIn = 0, nIntronsOut = 0 ;
		BSunit *uu ;
		
		units = arrayReCreate (units, 300, BSunit) ;
		
		bsGetArray (Transcript, _Splicing, units, 6) ;
		for (i = 0 ; i < arrayMax (units) ; i += 6)
		  {
		    uu = arrp (units, i, BSunit) ;
		    if (strstr (name (uu[4].k), "tron") && uu[5].s &&
			(strstr (uu[5].s, "gt_ag") || strstr (uu[5].s, "gc_ag"))
			)
		      {
			if (prod_p->x1 < uu[2].i && prod_p->x2 > uu[3].i)
			  nIntronsIn++ ;
			else
			  nIntronsOut++ ;
		      }
		  }
		
		if (nIntronsIn)
		  bsAddData (Product, str2tag ("Nb_introns_in_cds"), _Int, &nIntronsIn) ;
		if (nIntronsOut)
		  bsAddData (Product, str2tag ("Nb_introns_outside_cds"), _Int, &nIntronsOut) ;
	      }
	      
	      /* mRNA is filled in the mrna object */
	      /* peptide after the save Product */
	      {
		int lastu2 = 0 ;
		ORFS *orf ;
		Array hits ;
		
		orf = arrp (smrna->dnas, smrna->bestDna, ORFS) ;  
		hits = orf->hits ;
		
		bsRemoveTag (Product, _Source_Exons) ;
		for (jj = 0, up = arrp (hits, 0, HIT) ; jj < arrayMax (hits) ; up++, jj++)
		  {
		    
		    if (up->type & (gX | gGap))
		      {
			u1 = up->x1 - givenBack > x1 ? up->x1 - givenBack : x1 ;
			u2 = up->x2 - givenBack < x2 ? up->x2 - givenBack : x2 ;
			if (u2 >= u1)
			  {
			    u2 = u2 - x1 + 1 ; u1 = u1 - x1 + 1 ;
			    if (u1 > lastu2 + 1)
			      {
				lastu2++ ;
				bsAddData (Product, _Source_Exons, _Int, &lastu2) ;
				lastu2 = u1 - 1 ;
				bsAddData (Product, _bsRight, _Int, &lastu2) ;
				bsPushObj (Product) ;
				bsAddTag (Product, _Gap) ;
				bsGoto (Product, 0) ;
			      }
			    lastu2 = u2 ;
			    bsAddData (Product, _Source_Exons, _Int, &u1) ;
			    bsAddData (Product, _bsRight, _Int, &u2) ;
			    bsPushObj (Product) ;
			    if (up->type & gX)
			      bsAddTag (Product, _Exon) ;
			    else if (up->type & gOpenGap)
			      bsAddTag (Product, str2tag ("ORF_Gap")) ;
			    else if (up->type & gGap) 
			      bsAddTag (Product, _Gap) ;
			    bsGoto (Product, 0) ;
			  }
		      }
		  }
	      }
	      /* up and down stops */
	      if (orf->start > 0)
		{
		  int j, nn ;
		  char *cp ;
		  Array mdna = 0 ;
		  
		  mdna = arr (smrna->dnas, orf->iDna, ORFS).dna ;
		  /* jmax = arrayMax(mdna) ; */
		  
		  j = orf->start  - 1 ;
		  for (nn = 0, cp = arrp (mdna, j, char) ; j >= givenBack ; cp -= 3, j -= 3)
		    if (e_codon (cp, translationTable) == '*')
		      {
			x = j - (orf->start - 1) ;
			if (!nn++)
			  bsAddTag (Product, str2tag("Up_stop")) ;
			bsAddData (Product, _bsRight, _Int, &x) ; 
			/* break ; to show just one, changed to show all of them june 30, danielle */
		      }
		  
		  if (orf->downStop > 0)
		    {
		      x = orf->max - orf->start + 2 ;
		      bsAddTag (Product, str2tag("Down_stop")) ;
		      bsAddData (Product, _bsRight, _Int, &x) ; 
		      /* show just one */
		    }
		  else
		    {
		      bsRemoveTag (Product, str2tag("Down_stop")) ;
		    }
		}        
	      if (bsFindTag (Product, str2tag("NH2_Complete")) && bsFindTag (Product, str2tag("COOH_Complete")))
		bsAddTag (Product, str2tag("Complete")) ;
	      else
		bsAddTag (Product, str2tag("Partial")) ;
	      
	      bsSave (Product) ;
	      /* peptide */
	      if (orf->cds > 0)
		{
		  int j, k, nn = 0 ;
		  char *cp, *cq ;
		  Array mdna = arr (smrna->dnas, orf->iDna, ORFS).dna ;
		  
		  j = orf->start - 1 ;  /* will start in frame */
		  nn = orf->cds ;
		  if (orf->start < 0)
		    j += 3 ;
		  if (1 && orf->downStop > 0) 
		    nn += 3 ;  /* 2023_11_15 DO export, was DO NOT export the stop */
		  pep = arrayCreate (nn/3 + 4, char) ;
		  array(pep, nn/3, char) = 0 ; arrayMax(pep)-- ;  /* make room, zero terminate */
                  if (!strcmp (name(product), "G_tchr22_chr22|NC_060946.1_86_768.b3"))
		    invokeDebugger () ;
		  for (cp = arrp(pep, 0, char), cq = arrp (mdna, j, char), k = 0 ; nn > 0 ; cp++, cq += 3, nn -= 3, k++)
		    {
		      char cc = e_codon (cq, translationTable) ;
		      *cp = pepEncodeChar[(int)cc] ;
		      if (cc == '*')
			nn = -1 ; /* will break */
		    }
		  *cp = 0 ; arrayMax (pep) = k ;
		  if (orf->start == orf->leu) /* force a met */
		    {
		      cq = arrp (mdna, j, char) ; cp = arrp(pep, 0, char) ;
		      *cp = pepEncodeChar[(int)'M'] ;
		    }
		}
	      if (pep)
		{
		  mrnaDesignSavePepAndKantor (product, pep) ;
		  arrayDestroy (pep) ;
		}
	    }
	}
    }
  
  /*** edit the splicing tables ***/
  if (givenBack > 0 || stolen > 0)
    {
      char *dummy ;
      KEY dummyKey ;
      int x, ii ;
      BSMARK mark = 0 ;
#ifdef JUNK
      this does not work because i cannot do
       getArray AddArray in the presence of _Text
       because u.s is a pointer to the bs cell
         which gets destroyed in addarray by bsRemove()

       if (bsGetArray (Transcript, _Splicing, units, 9))
        {
          int i ;
          BSunit *uu ; 
          
          for (i = 0 ; i < arrayMax(units) ; i += 9)
            {
              uu = arrp (units, i, BSunit) ; 
              if (i > 0) uu[0].i += stolen ;
              uu[1].i += stolen ;
              if (uu[2].i > 0) uu[2].i -= givenBack ; 
              if (i == 0) uu[2].i -= stolen ;
              if (uu[3].i > 0) uu[3].i -= givenBack ; 
              if (i == 0)  /* fix length of exon */
                  uu[7].i += stolen ;
            }
          bsAddArray (Transcript, _Splicing, units, 9) ;
        }
#endif
       if (ii = 0, bsGetData (Transcript, _Splicing, _Int, &x))
        do 
          { /* totosp */
            mark = bsMark (Transcript, mark) ;
            
            if (ii++) x += stolen ;
            bsAddData (Transcript, _bsHere, _Int, &x) ;
            bsGetData (Transcript, _bsRight, _Int, &x) ;
            x += stolen ;
            bsAddData (Transcript, _bsHere, _Int, &x) ;
            bsGetData (Transcript, _bsRight, _Int, &x) ;
            if (x > 0) x -= givenBack ; 
            if (ii == 1) x -= stolen ;
            bsAddData (Transcript, _bsHere, _Int, &x) ;
            bsGetData (Transcript, _bsRight, _Int, &x) ;
            if (x > 0) x -= givenBack ;
            bsAddData (Transcript, _bsHere, _Int, &x) ;
            if (ii == 1)
              { /* fix length of exon */
                bsGetKeyTags (Transcript, _bsRight, &dummyKey) ;
                bsGetData (Transcript, _bsRight, _Text, &dummy) ; 
                bsGetData (Transcript, _bsRight, _Text, &dummy) ;
                bsGetData (Transcript, _bsRight, _Int, &x) ;
                x += stolen ; 
                bsAddData (Transcript, _bsHere, _Int, &x) ;
              }
            bsGoto (Transcript, mark) ;
          } while ( bsGetData (Transcript, _bsDown, _Int, &x)) ;

#ifdef JUNK
      if (ii = 0, bsGetData (Transcript, _Constructed_from, _Int, &x))
        do 
          {
            mark = bsMark (Transcript, mark) ;
            
            x -=  givenBack ;
            bsAddData (Transcript, _bsHere, _Int, &x) ;
            bsGetData (Transcript, _bsRight, _Int, &x) ;
            x -=  givenBack ;
            bsAddData (Transcript, _bsHere, _Int, &x) ;
            bsGoto (Transcript, mark) ;
          } while ( bsGetData (Transcript, _bsDown, _Int, &x)) ;
#endif
      if (bsGetArray (Transcript, _Constructed_from, units, 5))
        {
          int i ;
          BSunit *uu ;
          for (i = 0 ; i < arrayMax(units) ; i += 5)
            {
              uu = arrp (units, i, BSunit) ;
              uu[0].i -=  givenBack ;
              uu[1].i -=  givenBack ;
            } 
          bsAddArray (Transcript, _Constructed_from, units, 5) ;
        }

      if (bsGetData (Transcript, _Covers, _Int, &x))
        {        
          x += stolen ;
          bsAddData (Transcript, _bsHere, _Int, &x) ;
          bsGetData (Transcript, _bsRight, _Text, &dummy) ;
          bsGetData (Transcript, _bsRight, _Int, &x) ;
          if (sc->a1 < sc->a2)
            x -= stolen ;
          else
            x += stolen ;
          bsAddData (Transcript, _bsHere, _Int, &x) ;
        } 
      if (bsFindKey (Transcript, _IntMap, map) &&
          bsGetData (Transcript, _bsRight, _Int, &x))
        {
          if (dmap1 < dmap2)
            {
              x -= stolen ; 
              bsAddData (Transcript, _bsHere, _Int, &x)  ;
            }
          else
            {  
              x += stolen ; 
              bsAddData (Transcript, _bsHere, _Int, &x)  ;
            }
        }
      bsMarkFree (mark) ;
    }

  /*** number of various kinds of exons ****/
  {
    int nExons, nAlternativeExons, nStolenExons, nPredictedExons ;
    KEY type ;

    nExons = nAlternativeExons = nStolenExons = nPredictedExons = 0 ;

    if (bsGetArray (Transcript, _Splicing, units, 5))
      for (ii = 0; ii < arrayMax(units) ; ii += 5)
        {
          type = arr (units, ii + 4, BSunit).k ;
          if (type == _Exon) nExons++ ;
          if (type == _Stolen_exon) nStolenExons++ ;
          if (type == _Predicted_exon) nPredictedExons++ ;
          if (type == _Alternative_exon) nAlternativeExons++ ;
        }
    if (nExons)
      bsAddData (Transcript, str2tag("Nb_exons"), _Int, &nExons) ;
    if (nAlternativeExons)
      bsAddData (Transcript, str2tag("Nb_alternative_exons"), _Int, &nAlternativeExons) ;
    if (nStolenExons)
      bsAddData (Transcript, str2tag("Nb_stolen_exons"), _Int, &nStolenExons) ;
    if (nPredictedExons)
      bsAddData (Transcript, str2tag("Nb_predicted_exons"), _Int, &nPredictedExons) ;

  }

  bsSave (Transcript) ;
  smrna->stolen = stolen ;
  smrna->givenBack = givenBack ;
  mrnaSaveMrnaCoding (smrna) ;
  Transcript = bsUpdate (tr) ;
  
  /*** cdna_clones: Best, Specific, Complete_CDS ***/
  if (smrna->orfs && arrayMax(smrna->orfs) && smrna->clones)
  {
    KEYSET ks = smrna->clones ;

    char *condition = "(!Suspected_internal_deletion || Manual_no_internal_deletion ) && !Internal_priming && !Mosaic && !Partly_inverted" ;
    KEYSET goodClones = query (ks, condition) ;
    KEYSET specificClones = keySetCreate () ;
    KEYSET allClones = keySetCreate () ;
    int i, j, iSpecific = 0, x ;
    HIT *wp ;
    KEY clone ;
 
    if (arrayMax(smrnas) > 1)
      {  /* collect clones in other mrnas of same gene */
        int im, j, jj = 0 ;
        SMRNA *sm ;
        
        nonSpecificClones = keySetCreate () ;
        for (im = 0, sm = arrp (smrnas, im, SMRNA) ; im < arrayMax(smrnas) ; sm++, im++)
          if (sm->clones && sm != smrna)
            for (j = 0 ; j < arrayMax (sm->clones) ; j++)
              keySet (nonSpecificClones, jj++) = keySet (sm->clones, j) ;
        
        keySetSort (nonSpecificClones) ;
        keySetCompress (nonSpecificClones) ;
      }

    /* collect complete clones in a1 order */
    
    if (arrayMax (cloneHits))
      {
	cDNASwapA (cloneHits) ;
	arraySort (cloneHits, cDNAOrderGloballyByA1) ;
	
	for (j = 0, wp = arrp (cloneHits, 0, HIT) ; j < arrayMax(cloneHits) ; j++, wp++)
	  {
	    clone = wp->cDNA_clone ;
	    
	    /* register complete clones */
	    if (prod_array && arrayMax (prod_array))
	      {
		int i_prod ;
		OBJ Product = 0 ;
		
		for (i_prod = 0 ; i_prod < arrayMax (prod_array) ; i_prod++)
		  {
		    prod_p = arrp (prod_array, i_prod, HIT) ;
		    if (wp->a1 <= prod_p->x1+1 &&
			wp->a2 >= prod_p->x2 &&
			keySetFind (goodClones, clone, 0) &&
			(Product = bsUpdate (prod_p->gene)))
		      {
			bsAddKey (Product, str2tag ("Complete_CDS_clone"), clone) ;
			bsSave (Product) ;
		      }                
		  }
	      }
	    /* register specific clones */
	    if (nonSpecificClones &&
		!keySetFind (nonSpecificClones, clone, 0)) 
	      keySet (specificClones, iSpecific++) = clone ;
	    
	    /* register all clones and their size */
	    bsAddKey (Transcript, _cDNA_clone, clone) ;
	    
	    x = 0 ;
	    switch (wp->clipTop)
	      {
	      case 1: /* just 5' read */
		x = mrnaLength - wp->a1 + 1 ;
		break ;
	      case 2: /* just 3' read */
		x = wp->a2 ;
		break ;
	      case 3: /* both */
		if (wp->a1 > -20 && wp->a2 > 0)
		  x = wp->a2 - wp->a1 + 1 ;
		break ;
	      }
	    
	    bsAddData (Transcript, _bsRight, _Int, &x) ; 
	    bsAddData (Transcript, _bsRight, _Text, "bp,") ;
	    bsAddData (Transcript, _bsRight, _Int, &(wp->nerrAll)) ;
	    if (wp->nerrAll > 1)
	      bsAddData (Transcript, _bsRight, _Text, "errors") ;
	    else
	      bsAddData (Transcript, _bsRight, _Text, "error") ;
	    if (keyFindTag (clone, str2tag ("PCR_product_size")))
	      {
		OBJ Clone = bsCreate (clone) ;
		float y = 0 ;
		
		if (Clone)
		  {
		    if (bsGetData (Clone, str2tag ("PCR_product_size"), _Float, &y))
		      {
			KEY library = 0 ;
			OBJ Library = 0 ;
			
			/* deduce vector from pcr size to get insert size */
			if (bsGetKey (Clone, _Library, &library) &&
			    (Library = bsCreate (library)))
			  {
			    int dy = 0 ;
			    if (bsGetData (Library, str2tag ("Vector_length"),_Int, &dy ))
			      y -= dy/1000.0 ;
			    bsDestroy (Library) ;
			  }
			
			bsAddData (Transcript, _bsRight, _Float, &y) ;
			bsAddData (Transcript, _bsRight, _Text, "kb") ;
			y *= 1000 ;
			if (x > 100 && y > 100 && y > 1.1 * x)
			  {
			    int dx = (100 * (y - x))/x ;
			    bsAddData (Transcript, _bsRight, _Int, &dx) ;
			    bsAddData (Transcript, _bsRight, _Text, "%") ;
			  }
			else if (x > 100 && y > 100 && y < .9 * x)
			  {
			    int dx = (100 * (y - x))/x ;
			    bsAddData (Transcript, _bsRight, _Int, &dx) ;
			    bsAddData (Transcript, _bsRight, _Text, "%") ;
			  }
			else /* write something for the benefit of @ttach in the model */
			  {
			    int dx = 0 ;
			    bsAddData (Transcript, _bsRight, _Int, &dx) ;
			    bsAddData (Transcript, _bsRight, _Text, "%") ;
			  }
			
			bsDestroy (Clone) ;
		      }
		  }
	      }
	  }
      }

    /* register best clone */

    if (keySetMax (specificClones))
      {
        KEYSET ks = keySetAlphaHeap (specificClones, keySetMax (specificClones)) ;
        
        for (i = 0 ; i < keySetMax(ks) ; i++)
          bsAddKey (Transcript, str2tag ("Specific_clone"), keySet (ks, i)) ;
        keySetDestroy (ks) ;
      }

    if (keySetMax (allClones))
      {
        KEYSET ks = keySetAlphaHeap (allClones, keySetMax (allClones)) ;
        
        for (i = 0 ; i < keySetMax(ks) ; i++)
          bsAddKey (Transcript, _cDNA_clone, keySet (ks, i)) ;
        keySetDestroy (ks) ;
      }

    if (0) 
      {
        /* tout faux parce que completeClones ne sait pas qui est le best product 
           removed feb 12, 2005
        clone = mrnaBestAvailableClone (goodClones, completeClones, cloneHits) ;
        if (clone)
          bsAddKey (Transcript, str2tag ("Best_available_clone"), clone) ;
        */
      }
 
    keySetDestroy (goodClones) ;
    keySetDestroy (specificClones) ;
    keySetDestroy (allClones) ;
  }
    
  /*** Expression profile ***/
  if (smrna->clones)
    mrnaSaveExpressionProfile (Transcript, smrna->clones, nonSpecificClones) ;
  
  /*** PolyA signal ****/
  
  if (smrna->orfs && arrayMax (smrna->orfs))
    {
      KEY product = 0 ;
      
      if (bsFindTag (Transcript, str2tag("Found3p")) ||
	  (
	   bsGetKey (Transcript, _Product, &product) &&
	   keyFindTag (product, str2tag("COOH_Complete"))
	   ))
	{
	  ORFS * orf = arrp (smrna->orfs, 0, ORFS) ;
	  Array dna = arr (smrna->dnas, orf->iDna, ORFS).dna ;
	  unsigned char polyASignal[7] = { A_,A_,T_,A_,A_,A_, 0} ;
	  unsigned char polyASignal2[7] = { A_,T_,T_,A_,A_,A_, 0} ;
	  int nn = 0, n1 = 0, n2, debut = 0, dx ;
	  
	  debut = arrayMax(dna) - 40 ; /* was 40 before 2006_10_16 */
	  if (debut < 0) debut = 0 ;
	  if (debut > 10)
	    {
	      nn = dx = 0 ; /* look for signal closest to end of mRNA */
	      while ((n2 = dnaPickMatch (arrp (dna, debut + dx, unsigned char), 39 - dx, polyASignal, 0, 0)))
		{ nn = n2 + dx ; dx = nn ; }
	      if (!nn) /* favor ATTAAA */
		while ((n2 = dnaPickMatch (arrp (dna, debut + dx, unsigned char), 39 - dx, polyASignal2, 0, 0)))
		  { n1 = n2 + dx ; dx = n1 ; }
	      if (!nn && !n1)
		{
		  dx = 0 ;
		  while ((n2 = dnaPickMatch (arrp (dna, debut + dx, unsigned char), 39 - dx, polyASignal, 1, 0)))
		    { n1 = n2 + dx ; dx = n1 ; }
		}
	    }
	  if (nn > 0)
	    {
	      nn += debut ; nn = arrayMax(dna) - nn + 1 ;
	      bsAddData (Transcript, str2tag ("AATAAA"), _Int, &nn) ;
	    }
	  else if (n1 > 0 && n1+debut+5 < arrayMax(dna)) 
	    { 
	      int i ;
	      n1 += debut ;
	      
	      for (i = 0 ; i < 6 ; i++)
		polyASignal[i] = dnaDecodeChar [(int)arr (dna, n1 + i - 1, char)] ;
	      
	      n1 = arrayMax(dna) - n1 + 1 ;
	      bsAddData (Transcript, str2tag ("Variant"), _Text, polyASignal) ;
	      bsAddData (Transcript, _bsRight, _Int, &n1) ;
	    }
	}
    }
    
  bsSave (Transcript) ; /* so that the query(tr,..) work */
  Transcript = bsUpdate (tr) ; 

  /*** clone error against the dna of the mrna ***/
  if (smrna->orfs && arrayMax(smrna->orfs))
  {
    int a1, a2, x1, x2, a21, a12, n, nerr, x, ii, jj1, jj, jj2 ; 
    Array dna1 = arr (smrna->dnas, smrna->bestDna, ORFS).dna ;
    KEYSET clones = queryKey (tr, ">cDNA_clone") ;
    KEY est, clone ;
    Array errUp = 0, errDown = 0 ;
    Array dnaEst, dna2 = 0 ;
    BSunit *uu ;

    bsGetArray (Transcript, _Constructed_from, units, 5) ;
    for (ii = 0 ; ii < keySetMax (clones) ; ii++)
      {
        clone = keySet (clones, ii) ;
        jj1 = jj2 = -1 ; a21 = a12 = 0 ;
        for (n = jj = 0 ; jj < arrayMax (units) ; jj += 5)
          {
            uu = arrp (units, jj, BSunit) ;
            est = uu[2].k ;
            if (clone == keyGetKey (est, _cDNA_clone))
              {
                n++ ;
                a1 = uu[0].i ; a2 = uu[1].i ; 
                x1 = uu[3].i ; x2 = uu[4].i ; 
                if (x1 < x2) { jj1 = jj ; a21 = a2 ;} 
                if (x1 > x2) { jj2 = jj ; a12 = a1 ; }
              }
          }
        if (n != 2 || jj1 < 0 || jj2 < 0 || a21 < a12)
          continue ;

        if (!dna2)
          {
            dna2 = dnaCopy (dna1) ;
            reverseComplement (dna2) ;
          }
        uu = arrp (units, jj1, BSunit) ;
        est = uu[2].k ;
        a1 = uu[0].i + givenBack ; a2 = uu[1].i + givenBack ; 
        x1 = uu[3].i ; x2 = uu[4].i ; 
        dnaEst = dnaGet (est) ;
        if (!dnaEst) continue ;
        errDown = arrayReCreate (errDown, 64, A_ERR) ; 
        errDown = aceDnaDoubleTrackErrors (dnaEst, &x1, &x2, x1 < x2,
                                  dna1, dna2, &a1, &a2,
                                 0, errDown, 8, 0, FALSE, 0) ;
        uu = arrp (units, jj2, BSunit) ;
        est = uu[2].k ;
        a1 = uu[0].i + givenBack ; a2 = uu[1].i + givenBack ; 
        x1 = uu[3].i ; x2 = uu[4].i ; 
        dnaEst = dnaGet (est) ;
        if (!dnaEst) continue ;
        
        errUp = arrayReCreate (errUp, 64, A_ERR) ; 
        errUp = aceDnaDoubleTrackErrors (dnaEst, &x1, &x2, x1 < x2,
                                  dna1, dna2, &a1, &a2,
                                 0, errUp, 8, 0, FALSE, 0) ;

        /* find the best total intersection */
        nerr = mrnaCountCloneError (a12 + givenBack, a21 + givenBack, errDown, errUp) ;
        bsAddKey (Transcript, _cDNA_clone, clone) ;
        bsGetData (Transcript, _bsRight, _Int, &x) ; 
        bsAddData (Transcript, _bsRight, _Text, "bp,") ;
        bsAddData (Transcript, _bsRight, _Int, &nerr) ;        
        if (nerr > 1)
          bsAddData (Transcript, _bsRight, _Text, "errors") ;
        else
          bsAddData (Transcript, _bsRight, _Text, "error") ;
      }
    arrayDestroy (dna2) ;
    keySetDestroy (clones) ;
    arrayDestroy (errUp) ; arrayDestroy (errDown) ;
  }

  bsSave (Transcript) ;


  /*** edit coordinates of mrna in transcript ***/
  if (stolen)
    {
      int x ;
      OBJ Gene = bsUpdate (sc->gene) ;

      if (Gene)
        {
          if (bsFindKey (Gene, _mRNA, tr) &&
              bsGetData (Gene, _bsRight, _Int, &x))
            {
              x -= stolen ;
              bsAddData (Gene, _bsHere,  _Int, &x) ;
            }
          bsSave (Gene) ;
        }

    }

  /*** SMAP : S_Parent Comid info: after we bsUpdate gene ***/
  { 
    Cosmid = bsUpdate(relevantCosmid) ;
    if (Cosmid)
      {
        bsAddKey (Cosmid, str2tag("mRNAs"), tr) ;
        if (ta1 < ta2) ta1 -= stolen ;  
        else  ta1 += stolen ;  
        bsAddData (Cosmid, _bsRight, _Int, &ta1) ;
        bsAddData (Cosmid, _bsRight, _Int, &ta2) ;
        bsSave (Cosmid) ; 
      }
  }
  
  /*** SMAP: DNA  after we bssave the transcript ***/
  if (smrna->orfs && arrayMax(smrna->orfs))  
  {
    Array dna1, dna2 = 0 ;

    dna1 = arr (smrna->dnas, smrna->bestDna, ORFS).dna ;
    if (dna1 && arrayMax(dna1))
      {
        dna2 = dnaCopy (dna1) ;
        if (givenBack > 0)
          {
            int nn = arrayMax(dna1) - givenBack ;
            memcpy (arrp(dna2, 0, char), arrp(dna1, givenBack, char), nn) ;
            arr (dna2, nn, char) = 0 ;
            arrayMax(dna2) = nn ;
          }
        dnaStoreDestroy  (tr, dna2) ;
      }
  }

  /*** Tiling, can also be called later, + genomic erros, needs sc */

  keySetDestroy (nonSpecificClones) ;
  arrayDestroy (prod_array) ;
  arrayDestroy (feet) ;
  arrayDestroy (units) ;
  arrayDestroy (cloneHits) ;
} /* mrnaSaveMrna */
  
/**********************************************************************************/

static int mrnaCountGaps (Array smrnas)
{
  int i, ii, jj, mx, maxgap = 0 ;
  SMRNA *smrna ;
  HIT *vp ;

  if (arrayMax(smrnas))
    for (ii = jj = maxgap = 0, smrna = arrp (smrnas, 0, SMRNA)  ; ii < arrayMax(smrnas) ;  smrna++, ii++)
      { 
        smrna->nGaps = 0 ; 
        mx = arrayMax(smrna->hits) ;
        if (mx)
          for (i = 0, vp = arrp (smrna->hits, 0, HIT) ; i < mx ; i++, vp++)
            if (vp->type & gGap)
              smrna->nGaps++ ;
        if (maxgap < smrna->nGaps) 
          maxgap = smrna->nGaps ;
      }
  return maxgap ;
}

/**********************************************************************************/

static KEYSET mrnaCloneNoGap (Array smrnas)
{
  int i, ii, mx ;
  SMRNA *smrna ;
  KEY *kp ;
  KEYSET ks = keySetCreate () ;
  
  for (ii = 0 ; ii < arrayMax(smrnas) ; ii++)
    {
      smrna = arrp (smrnas, ii, SMRNA) ;
      mx = keySetMax(smrna->clones) ;
      if (mx && !smrna->nGaps)
        for (i = 0, kp = arrp (smrna->clones, 0, KEY); i < mx ; kp++, i++)
          keySetInsert (ks, *kp) ; 
    }
  return ks ;
}

/**********************************************************************************/
/**********************************************************************************/
/* after stealing mrnas may become identical, if so fuse their clones */

static int smrnaA1Order (const void *a, const void *b)
{
  const SMRNA *sa = (const SMRNA *)a, *sb = (const SMRNA *)b ;
  int da = sa->a1 - sb->a1 ;
  
  return da ;
}

/*************/

static BOOL mrnaFuseRedundant (S2M *s2m, SMRNA *gmrna, Array smrnas)
{
  int i0, ii, jj, i, j, a1, b1, iMax, jMax ;
  HIT *up, *vp ;
  SMRNA *smrna1, *smrna2 ;
  BOOL ok = FALSE, fused = FALSE ;
  int fuseGaps = 1 ; /* = mrnaPleaseFuseGaps () ;
		      * with 2: ca marche pas du tout ce truc
		      * parce que ce code ne sait pas changer le nombre 
		      * d'exons interne , seulement en rajouter dans les bouts
		      */
  if (arrayMax (smrnas) < 2) return FALSE;
  arraySort (smrnas, smrnaA1Order) ;
  /* search mrnas with same coordinates */
  for (ii = 0; ii < arrayMax(smrnas) ; ii++)
    {
      smrna1 = arrp (smrnas, ii, SMRNA) ;
      iMax = arrayMax (smrna1->hits) ; 
      if (!iMax)
        continue ;
      a1 = smrna1->a1 ;
      for (jj = ii + 1; jj < arrayMax(smrnas) ; jj++)
        {
          smrna2 = arrp (smrnas, jj, SMRNA) ;
          b1 = smrna2->a1 ; 
	  if (b1 > smrna1->a2) 
	    continue ;
          jMax = arrayMax (smrna2->hits) ; 
          if (!jMax)
            continue ;
          ok = TRUE ; i0 = -1 ;
          for (i = j = 0 ; ok && i < iMax && j < jMax ; )
            {
              up = arrp (smrna1->hits, i, HIT) ;
              vp = arrp (smrna2->hits, j, HIT) ;
              if ( ( (up->type & gI) && (vp->type & gX) ) ||    /* intron exon */
                   ( (up->type & gX) && (vp->type & gI) ) 
                   )
                {
                  if (
                      up->a2 + a1 > vp->a1 + b1 && /* they intersect */
                      vp->a2 + b1 > up->a1 + a1
                      )
                    {  ok = FALSE ; continue ; }
                }
              else if ((up->type & gI) && (vp->type & gI)) /* intron intron */
                {
                  if (up->a1 + a1 != vp->a1 + b1 ||
                      up->a2 + a1 != vp->a2 + b1)
                    {  ok = FALSE ; continue ; }
                }
              else if (fuseGaps != 2 && (up->type & gGap) && (vp->type & gGap))
                {
                  if (up->a1 + a1 != vp->a1 + b1 ||
                      up->a2 + a1 != vp->a2 + b1)
                    {  ok = FALSE ; continue ; }
                }
              else if (fuseGaps == 2 && (up->type & gGap) && (vp->type & gGap))
                {
                  ok = ok ;
                }
	      else if (
                       (fuseGaps && (up->type & gI) && (vp->type & gGap)) || 
                       (fuseGaps && (up->type & gGap) && (vp->type & gI))
                       )
                {
                  ok = ok ;
                }
              else if (
                       ((up->type & gX) && (vp->type & gX)) ||  /* exon exon */
                       (fuseGaps && (up->type & gX) && (vp->type & gGap)) || 
                       (fuseGaps && (up->type & gGap) && (vp->type & gX))
                       )
                {
                  if (i == 0 && (up->type & gS) && up->a1 + a1 > vp->a1 + b1)
		    ok = FALSE ;
                  if (j == 0 && (vp->type & gS) && up->a1 + a1 < vp->a1 + b1)
		    ok = FALSE ;
		  if (ok && i == 0 && (up->type & gX) && (vp->type & gGap) &&
		      up->a1  + a1 > vp->a1 + b1) /* 1st vp exon does not contact up exon */
		    ok = FALSE ;
		  if (ok && j == 0 && (vp->type & gX) && (up->type & gGap) &&
		      vp->a1  + b1 > up->a1 + a1) /* 1st up exon does not contact vp exon */
		    ok = FALSE ;
		  if (ok &&up->a2 + a1 < vp->a1 + b1)
                    { 
                      if (i == iMax - 1 || /* this mrna is already over */
                          (vp->type & (gS | gS0 | gReal5p))
                          )
                        /* i cannot accept extra 5' exons */
                        ok = FALSE ; /* up intron totally above vp true start */
                      i++ ; continue ;
                    }
                  if (ok && i > 0 &&  up->a1  + a1 > vp->a1 + b1 &&
                      (((up-1)->type & gI) && vp->type & gX)
                      )
                    ok = FALSE ;     /* begin of vp is inside up intron */
                  if (i0 == -1) i0 = i ; /* first relevant i exon */
                  if (ok && up->a1 + a1 != vp->a1 + b1)
                    { 
                      if (j == 0 &&
                          up->a2 + a1 > vp->a1 + b1 && /* they intersect */
                          vp->a2 + b1 > up->a1 + a1 &&
                          ! (
                             ((up->type & gS) && up->a1 + a1 > vp->a1 + b1) ||
                             ((vp->type & gS) && up->a1 + a1 < vp->a1 + b1) 
                             )
                          ) ;
                      else if (up->type & gGap &&
                               ! (
                                  ((vp->type & gS) && up->a1 + a1 < vp->a1 + b1) 
                                  )
                               )  ;
                      else if (vp->type & gGap &&
                               ! (
                                  ((up->type & gS) && up->a1 + a1 > vp->a1 + b1) 
                                  )
                               ) ;
                      else if (
			       (i > 0 && ((up-1)->type & gGap)) ||
			       (j > 0 && ((vp-1)->type & gGap))
			       ) ;
                      else ok = FALSE ; /* incompatible true starts in intersecting introns */
                    }
                  if (ok && i == iMax -1 && j < jMax - 1)
                    {
                      if ( up->a2 + a1 > vp->a1 + b1 && /* they intersect */
                           vp->a2 + b1 > up->a1 + a1 &&
                           ! (up->a2 + a1 > vp->a2 + b1) && /* up exon leaks in vp intron */
                           ! (up->type & (gReal3p | gA))
                           ) ; /* up is one exon  shorter and is a real3p */
                      else if (vp->type & gGap &&
                               ! (up->type & (gReal3p | gA))
                               )  ;
                      else ok = FALSE ; /* up is a true end */
                    }
                  else if (ok && j == jMax -1 && i < iMax - 1)
                    {
                      if ( up->a2 + a1 > vp->a1 + b1 && /* they intersect */
                           vp->a2 + b1 > up->a1 + a1 &&
                           ! (((up+1)->type & gI) && up->a2 + a1 < vp->a2 + b1) && /* vp exon leaks in up intron */
                           ! (vp->type & (gReal3p | gA))
                           ) ; /* vp is one exon shorter and not a real3p */
                      else if (up->type & gGap &&
                               ! (vp->type & (gReal3p | gA))
                               )  ;
                      else ok = FALSE ;  /* vp is a true end */
                    }
                  else if (ok && up->a2 + a1 != vp->a2 + b1)
                    {  /* fuse last exons of i and j */
                      if (j == jMax - 1 &&
                          up->a2 + a1 > vp->a1 + b1 && /* they intersect */
                          vp->a2 + b1 > up->a1 + a1 
                          ) ;
                      else if (up->type & gGap &&
                               ! (vp->type & (gReal3p | gA))
                               )  ;
                      else if (vp->type & gGap &&
                               ! (up->type & (gReal3p | gA))
                               )  ;
                      else ok = FALSE ;
                    }
                }
              else /* if (i==0 && j==0)   anything else: gaps etc */
                {  ok = FALSE ; continue ; }

              if (ok && i == iMax -1 && j < jMax - 1 &&
                  (vp->type & gGap) &&
                  ( up->type & (gReal3p | gA))
                  ) 
                ok = FALSE ; /* up is a true end */
              if (ok &&  j == jMax -1 && i < iMax - 1 &&
                  (up->type & gGap) &&
                  ( vp->type & (gReal3p | gA))
                  ) 
                ok = FALSE ; /* up is a true end */
              if (up->a2 + a1 < vp->a2 + b1)
                i++ ; /* do NOT increase both */
              else
                j++ ;
              continue ;
            }
          if (i0 >= 0 && ok) /* fuse jj in ii */
            {
              int j = keySetMax (smrna1->clones) ;

              smrnaDnaDestroy (smrna1) ;
              smrnaDnaDestroy (smrna2) ;

              for (i = 0 ; i < keySetMax (smrna2->clones) ; i++)
                keySet (smrna1->clones, j++) = keySet (smrna2->clones, i) ;
              keySetSort (smrna1->clones) ;
              keySetCompress (smrna1->clones) ;
              for (i = i0, j = 0 ; i  < iMax && j < jMax ; i++, j++)
                {
                  up = arrp (smrna1->hits, i, HIT) ;
                  vp = arrp (smrna2->hits, j, HIT) ;
		  if ((up->type & gX) && i < iMax - 1 && ((up+1)->type & gGap) &&
		      (vp->type & gX) &&  /* vp exon extends up-1 exon into  up gap */
		      vp->a1 + b1 < up->a2 + a1 &&
		      vp->a2 + b1 > up->a2 + a1
		      ) 
		    {  /* extend the up exon */
		      int x = (vp->a2 + b1 < (up+1)->a2 + a1) ? vp->a2 + b1 - a1 : (up+1)->a2 - 1 ;
			
		      up->a2 = x ; (up+1)->a1 = x+1 ;
		    }
                  if (up->type != vp->type)
                    {
                      unsigned int zz ;
                      zz = (up->type | vp->type) & (gReal3p | gReal5p | gCompleteCDS | gA | gS | gS0) ;
                      up->type |= zz ;
                      if ((up->type & gPredicted) && !(vp->type & gPredicted))
                        up->type &= ~gPredicted ;
                      if ((up->type & gStolen) && !(vp->type & gStolen))
                        up->type &= ~gStolen ;
                      if ((up->type & gJ) && !(vp->type & gJ))
                        up->type &= ~gJ ;
                      if ((up->type & gGap) && !(vp->type & gGap))
                        up->type = vp->type ;
                    }
                  if (i == 0 &&  up->a1  + a1 > vp->a1 + b1)
                    up->a1 = vp->a1 + b1 - a1 ;
                  if (i == iMax - 1 && up->a2  + a1 < vp->a2 + b1)
                    up->a2 = vp->a2 + b1 - a1 ;
                }
	      if (iMax - i0 >= jMax)
		{
		  up = arrayp (smrna1->hits, iMax - 1, HIT) ;
                  vp = arrp (smrna2->hits, jMax - 1, HIT) ;
                  if (up->a2 + a1 > vp->a1 + b1 &&
		      up->a2 + a1 < vp->a2 + b1)
		    up->a2 = vp->a2 + b1 - a1 ;
		}
              for (i = iMax, j = i - i0 ; j < jMax ; i++, j++)
                {
                  up = arrayp (smrna1->hits, i, HIT) ;
                  vp = arrp (smrna2->hits, j, HIT) ;
                  *up = *vp ;
                  up->a1 += b1 - a1 ;                  
                  up->a2 += b1 - a1 ;
                }
              if (smrna2->a2 > smrna1->a2)
                smrna1->a2 = smrna2->a2 ;
              arrayMax (smrna2->hits) = 0 ;
              iMax = arrayMax (smrna1->hits) ;
            }
        }
    } 
  for (ii = jj = 0; ii < arrayMax(smrnas) ; ii++)
    {
      smrna1 = arrp (smrnas, ii, SMRNA) ;
      if (! (arrayMax (smrna1->hits)))
        { fused = TRUE ; continue ; }
      if (jj < ii)
        {
          smrna2 = arrp (smrnas, jj, SMRNA) ;
          *smrna2 = *smrna1 ;
        }
      jj++ ;
    }
  arrayMax(smrnas) = jj ; 
  return fused ;
} /* mrnaFuseRedundant */

/**********************************************************************************/
/**********************************************************************************/
#ifdef JUNK_SEE_BELOW
/* oct-6 2004
 * it happens that there is a long (genomic contaminant ?) est ABCDE
 * now we have 2 read (intron-A exon-B) and (exon-D;intron-E)
 * makechain has considered BCD as a bonafide exon because there is
 * no intron boundary inside, and put the 2 reads in a single mrna
 * and no gap is reported since (BCD) exon is touched
 * and yes if stealing is allowed, this is all perfect
 * but is stealing is disallowed, i wish to split on C
 * so here i make the gap explicit, it will be split 
 * in remove redundant
 */

static BOOL mrnaCutIfGapInMrna (S2M *s2m, SC *sc, SMRNA *gmrna, Array smrnas)
{
  int nn = 0 ;
  int ii, iiMax, iUp, j, jVp, kWp, a1, a2 ;
  KEY oldclo ;
  SMRNA *smrna, *smrnaNew ;
  HIT *up, *vp, *wp ;
  Array cHits = 0, dHits = 0 ;

  /* cut mrna with gaps in pieces 
   * accumulate the rejected part at the end of the smrnas array
   * so they get searched later
   */
  iiMax = arrayMax (smrnas) ; /* array may extend but loop should not */
  for (ii = 0 ; ii < iiMax ; ii++)
    {
      smrna = arrp (smrnas, ii, SMRNA) ; /* may get relocated */

      /* we know all clones are compatible with this mrna
       * so we are only interested in the clone full extent measured on the genome
       * and we merge all connected segments
       */
      cHits = arrayReCreate (cHits, 256, HIT) ; jVp = a1 = a2 = 0 ;
      for (iUp = 0, up = arrp (gmrna->estHits, 0, HIT), oldclo = 0 ; iUp < arrayMax (gmrna->estHits) ; iUp++, up++)
        {
          if (!(up->type & gX))
            continue ;
          if (up->cDNA_clone == oldclo)
            {
              if (a1 > up->a1) a1 = up->a1 ;
              if (a2 < up->a2) a2 = up->a2 ; /* extend */
            }
          else if (! keySetFind (smrna->clones, up->cDNA_clone, 0))
            continue ;
          else /* register */
            {
              if (oldclo)
                {
                  vp = arrayp (cHits, jVp++, HIT) ;
                  vp->a1 = a1 ; vp->a2 = a2 ; vp->cDNA_clone = oldclo ; 
                }
              oldclo = up->cDNA_clone ; a1 = up->a1 ; a2 = up->a2 ;
            }
        }
      if (oldclo)
        {
          vp = arrayp (cHits, jVp++, HIT) ;
          vp->a1 = a1 ; vp->a2 = a2 ; vp->cDNA_clone = oldclo ;
        }
      if (jVp < 2)  /* single clone */
        continue ; /* loop on next ii:smrna */
      /* sort merge cHits */
      kWp = 0 ; wp = 0 ;
      arraySort (cHits, cDNAOrderGloballyByA1) ;
      for (j = 0, vp = arrayp (cHits, 0, HIT), a1 = vp->a1, a2 = vp->a2 ;
           j < jVp ; j++, vp++)
        {
          if (vp->a1 < a2 && /* intersect */
              vp->a2 > a1)
            { if (a2 < vp->a2) a2 = vp->a2 ;    /* extend */ }
          else /* start a new contig */
            {
              if (!kWp) dHits = arrayReCreate (dHits, 6, HIT) ;
              wp = arrayp (dHits, kWp++, HIT) ;
              wp->a1 = a1 ; wp->a2 = a2 ;
              a1 = vp->a1 ; a2 = vp->a2 ;
            }
        }
      if (kWp) /* we have started a second contig, finish it */
        { 
          wp = arrayp (dHits, kWp++, HIT) ;
          wp->a1 = a1 ; wp->a2 = a2 ;
        }
      else  /* there is a single contig, this mrna is ok */
        continue ; /* loop on next ii:smrna */
      
      /* we now cut out the first contig and copy it 
       * at the end of the smrna list
       * note there will not need to reexamine it
       */
      for (kWp = 0 ; kWp < arrayMax(dHits) - 1 ; kWp++)
        {
          int k1 = 0, k2 = 0 ;

          smrnaNew = arrayp (smrnas, arrayMax (smrnas), SMRNA) ; /* may get relocated */
          smrna = arrp (smrnas, ii, SMRNA) ; /* may just have been relocated */
          smrnaDnaDestroy (smrna) ;
          *smrnaNew = *smrna ; 
          smrnaNew->clones = arrayHandleCreate (128, KEY, s2m->h) ;
          keySetMax (smrna->clones) = 0 ;
          smrnaNew->hits = arrayHandleCopy (smrna->hits, s2m->h)  ;

          /* move all relevant clones from smrna to smrnaNew */
          for (j = k1 = k2 = 0, vp = arrayp (cHits, 0, HIT) ; j < jVp ; j++, vp++)
            if (vp->a1 <= wp->a2 && wp->a1 <= vp->a2) /* intersect */
              keySet (smrnaNew->clones, k2++) = vp->cDNA_clone ;
            else
              keySet (smrna->clones, k1++) = vp->cDNA_clone ;
          keySetSort (smrna->clones) ;
          keySetSort (smrnaNew->clones) ;
          /* do not change the hits, they will be grignotted later */
          nn++ ;
        }
    }
  arrayDestroy (cHits) ;
  arrayDestroy (dHits) ;

  return nn ;
} /* mrnaCutIfGapInMrna */
#endif
/**********************************************************************************/
/* feb 7, 2007
 * new algo
 * we prefer to keep a maximum of 1 gap per mrna
 * we separate the mrna on all the gaps and call
 * k1, k2,    kn the clones supporting the n contigs
 * for each pair ki kj
 *   if kij = ki inter kj is not empty
 *   we do  ki <- kij  kj <- kj - kij kinew = ki - kij
 *   and we iterate
 * at the end, each non empty set ki defines s separate mrna
 */

static int mrnaNewCutGapInMrnaGetClones (KEYSET ks, int a1, int a2, SMRNA *gmrna, SMRNA *smrna)
{
  int iUp, i ;
  HIT *up ;
  KEY clone ;
  KEYSET ks2 ;
  for (iUp = 0, up = arrp (gmrna->estHits, 0, HIT) ; iUp < arrayMax (gmrna->estHits) ; iUp++, up++)
    {
      if (!(up->type & gX) || up->a2 < a1 + 5 || up->a1 > a2 - 5)
        continue ;
      clone = up->cDNA_clone ;
      if (! keySetFind (smrna->clones, clone, 0))
        continue ;
      if ( keySetFind (ks, clone, 0))
        continue ;
      keySetInsert (ks, clone) ;
    }

  ks2 = query (ks, ">fuse_to_clone") ;
  for (i = 0 ; i < keySetMax (ks2) ; i++)
    keySetInsert (ks, keySet(ks2,i)) ;
  
  return keySetMax (ks) ;
} /* mrnaNewCutGapInMrnaGetClones */

/****************/

static int mrnaNewCutIfGapInMrna (S2M *s2m, SC *sc, SMRNA *gmrna, Array smrnas, int bigGap)
{
  int a1, iks, iMrna, iMrnaMax, oldMrnaMax, ii, jj ;
  SMRNA *smrna, *smrna1 ;
  HIT *up ;
  Array kss = 0 ;
  KEYSET ks = 0 ;
  BOOL hasIntron ;

  /* cut mrna with gaps in pieces 
   * accumulate the rejected part at the end of the smrnas array
   * so they do not get searched twice 
   */
  iMrnaMax = arrayMax (smrnas) ; /* array may extend but loop should not */
  for (a1 = iMrna = 0 ; iMrna < iMrnaMax ; iMrna++)
    {
      smrna = arrp (smrnas, iMrna, SMRNA) ; /* may get relocated */
      /* do not cut intronless mrnas */
      hasIntron = FALSE ;
      for (jj = 0 ; jj < arrayMax (smrna->hits) ; jj++)
        {
          up = arrp (smrna->hits, jj, HIT) ;
          if (gI & up->type)
	    { hasIntron = TRUE ; break ; }
          if (bigGap && (gGap & up->type) && up->a2 > up->a1 + bigGap)
	    { hasIntron = TRUE ; break ; }
	}
      if (! hasIntron)
	continue ;
      iks = 0 ; 
      kss = arrayReCreate (kss, 12, KEYSET) ;
      for (jj = 0 ; jj < arrayMax (smrna->hits) ; jj++)
        {
          up = arrp (smrna->hits, jj, HIT) ;
          if ((gGap & up->type) && up->a2 - up->a1 > 100)
            {
              ks = array (kss, iks++, KEYSET) = keySetCreate () ;
              mrnaNewCutGapInMrnaGetClones (ks, smrna->a1 + a1 - 1, smrna->a1 + up->a1 - 2
                                            , gmrna, smrna) ;
              a1 = up->a2 + 1 ;
            }
        }
      if (iks == 0) /* this mrna has no gap */
        continue ; 
      ks = array (kss, iks++, KEYSET) = keySetCreate () ; /* get the contig downstream of last gap */
      mrnaNewCutGapInMrnaGetClones (ks, smrna->a1 + a1 - 1, smrna->a2
                                    , gmrna, smrna) ;
      iks = keySetPartition (kss) ; /* create independant clone subsets */
      /* create one mrna for each set */
      if (iks >= 2)
        {      
          oldMrnaMax = arrayMax (smrnas) ;
          ii = oldMrnaMax + iks - 2 ;
          smrna1 = arrayp (smrnas, ii, SMRNA) ; /* make room */
          smrna = arrp (smrnas, iMrna, SMRNA) ; /* may get relocated */
          keySetDestroy (smrna->clones) ;
          ks = arr (kss, 0, KEYSET) ;
          smrna->clones =  arrayHandleCopy (ks, s2m->h) ;
          smrnaDnaDestroy (smrna) ;
          for (ii = 1 ; ii < iks ; ii++)
            {
              smrna1 = arrayp (smrnas, oldMrnaMax + ii - 1, SMRNA) ;
              *smrna1 = *smrna ;
              smrna1->hits = arrayHandleCopy (smrna->hits, s2m->h)  ;
              ks = arr (kss, ii, KEYSET) ;
              smrna1->clones = arrayHandleCopy (ks, s2m->h) ;
              /* do not change the hits, they will be grignotted later */
            }
        }
      for (ii = 0 ; ii < iks ; ii++)
        {
          if ((ks = arr (kss, ii, KEYSET)))
            keySetDestroy (ks) ;
        }
    }
  arrayDestroy (kss) ;
  return arrayMax (smrnas) - iMrnaMax ; /* number of new mrnas */
} /* mrnaNewCutIfGapInMrna */

/**********************************************************************************/
/* 
 * do not allow the same clone to belong to a full and to a gapped mrna
 */
static BOOL mrnaPreferGaplessMrna (S2M *s2m, SC *sc, SMRNA *gmrna, Array smrnas, BOOL notJustGaps)
{
  int i, ii, jj, nn = 0, maxgap ;
  SMRNA *smrna, *smrna2 ;
  KEYSET kss = 0, ks1 = 0, ks2 = 0 ;

  if (! notJustGaps)
    {
      maxgap = mrnaCountGaps (smrnas) ;
      if (!maxgap)
        return 0 ;
      arraySort (smrnas, smrnaOrderByGap) ;
    }
  else
    { 
      mrnaConstructOrfs (s2m, sc, gmrna, smrnas) ;
      arraySort (smrnas, smrnaOrderByLengthandBadQuality) ;
    }
  /* if mrna has gap and clone is  also in another previous mrna
     that's an oeuf (enough)
  */
  for (ii = 0 ; ii < arrayMax(smrnas) ; ii++)
    { 
      smrna = arrp (smrnas, ii, SMRNA) ; /* may be relocalized */
      if (! notJustGaps && ! smrna->nGaps)
        continue ;
      ks1 = smrna->clones ;
      for (jj = 0 ; jj < ii ; jj++)
        { 
          smrna2 = arrp (smrnas, jj, SMRNA) ; /* may be relocalized */
          smrnaDnaDestroy (smrna2) ;

          ks2 = smrna2->clones ;
          kss = keySetMINUS (ks1, ks2) ;
          if (keySetMax (kss) < keySetMax (ks1))
            {
              nn++ ;
              for (i = 0 ; i < keySetMax (kss) ; i++)
                keySet (ks1, i) = keySet (kss, i) ;
              keySetMax (ks1) = i ;
              mrnaGrignotte (s2m, sc, gmrna, smrnas, ii, FALSE) ;  /* grignotte all, needed a second time  */ 
            }
          keySetDestroy (kss) ;
        }
    }
  return nn ;
} /* mrnaPreferGaplessMrna */

/**********************************************************************************/

static BOOL mrnaRemoveRedundant (S2M *s2m, SC *sc, SMRNA *gmrna, Array smrnas)
{
  int i, j, ii, jj, iGap, iUp, xGap, mx, maxgap ;
  SMRNA *smrna ;
  KEYSET ks = 0 ;
  BOOL redundant, ok = FALSE ;
  KEY *kp, fused = 0 ;
  KEYSET kss = 0, ks1 = 0, ks2 = 0 ;
  HIT *up, *vp, *wp ;

  if (arrayMax(smrnas) < 2)
    goto done ;

  maxgap = mrnaCountGaps (smrnas) ;
  if (!maxgap)
    goto done ;

  arraySort (smrnas, smrnaOrderByGap) ;


  /* remove mrnas totally included in previous mrna */
  ks = keySetCreate () ;
  for (ii = 0, jj = 0 ; ii < arrayMax(smrnas) ; ii++)
    {
      redundant = TRUE ;

      smrna = arrp (smrnas, ii, SMRNA) ;
      mx = keySetMax(smrna->clones) ;
      if (!mx)
        continue ;

      redundant = TRUE ;
      for (i = 0, kp = arrp (smrna->clones, 0, KEY); i < mx ; kp++, i++)
        if (! keySetFind (ks, *kp, 0))
          { redundant = FALSE ;  keySetInsert (ks, *kp) ; }
        
      if (!redundant)
        {
          if (jj != ii)
            {
              SMRNA *smrna2 = arrp (smrnas, jj, SMRNA) ;
              *smrna2 = *smrna ;
            }
          jj++ ;
        }
    }
  arrayMax (smrnas) = jj ;
  if (arrayMax(smrnas) < 2)
    goto done ;

 allgapAgain:
  maxgap = mrnaCountGaps (smrnas) ; 
  arraySort (smrnas, smrnaOrderByGap) ;
  
  if (!maxgap)
    goto done ;

  mrnaCloneNoGap (smrnas) ;

  /* cut mrna with gaps in pieces and remove redundant pieces */
  for (ii = 0 ; ii < arrayMax(smrnas) ; ii++)
    {
      redundant = TRUE ;

      smrna = arrp (smrnas, ii, SMRNA) ;
      if (!smrna->nGaps)
        continue ;
      /* locate 1st gap */
      mx = arrayMax(smrna->hits) ;
      iGap = -1 ;
    xGapAgain: 
      xGap = 0; 
      ks1 = keySetCreate () ;
      ks2 = keySetCreate () ;
      kss = ks1 ;
      if (mx)
        for (i = 0, vp = arrp (smrna->hits, 0, HIT) ; i < mx ; i++, vp++)
          {
            if (vp->type & gX)
              {
                int z1 = smrna->a1 + vp->a1 - 1 , z2 = smrna->a1 + vp->a2 - 1 ;
                for (iUp = 0, up = arrp (gmrna->estHits, 0, HIT) ; iUp < arrayMax (gmrna->estHits) ; iUp++, up++)
                  if ((up->type & gX) && up->a1 < z2 && up->a2 > z1 &&
                      keySetFind (smrna->clones, up->cDNA_clone, 0) )
                   {
                     keySetInsert (kss, up->cDNA_clone) ;
                     if ((fused = keyGetKey (up->cDNA_clone, str2tag ("Fuse_to_clone"))))
                       keySetInsert (kss, fused) ;
                   }
              }
            if ((vp->type & gGap) && !xGap && i > iGap)
              { xGap = (vp->a1 + vp->a2)/2 ; kss = ks2 ; iGap = i ; }
          }
      kss = 0 ; /* we do not want to destroy kss */
      if (keySetMax(ks1) && keySetMax(ks2)) /* interesting case */
        {
          KEYSET ks0 = 0 ;
          
          ks0 = keySetAND (ks1, ks2) ;
          if (keySetMax (ks0)) /* this is a gap bridged by a specific clone */
            {
              keySetDestroy (ks0) ;
              keySetDestroy (ks1) ;
              keySetDestroy (ks2) ;
              goto xGapAgain ;
            } 
          keySetDestroy (ks0) ;
          if (1)
            {
              /* kill the first segment if it has no specific clone */
              if (!keySetMax(ks1))
                {
                  /* do it this way because smrna->clones is allocated on handle */
                  ks0 = keySetMINUS (smrna->clones, ks1) ;
                  for (i = 0 ; i < keySetMax (ks0) ; i++)
                    keySet (smrna->clones, i) = keySet (ks0, i) ;
                  keySetMax (smrna->clones) = i ;
                  keySetDestroy (ks0) ;
                  
                  for (i = 0, vp = arrp (smrna->hits, 0, HIT) ; i <= iGap ; i++, vp++)
                    vp->type = 0 ;
                  ok = TRUE ;
                }
              /* kill the second segment if it has no specific clone */
              else if (!keySetMax(ks2))
                {
                  /* do it this way because smrna->clones is allocated on handle */
                  ks0 = keySetMINUS (smrna->clones, ks2) ;
                  for (i = 0 ; i < keySetMax (ks0) ; i++)
                    keySet (smrna->clones, i) = keySet (ks0, i) ;
                  keySetMax (smrna->clones) = i ;
                  keySetDestroy (ks0) ;
                  
                  for (i = iGap, vp = arrp (smrna->hits, i, HIT) ; i < mx ; i++, vp++)
                    vp->type = 0 ;
                  ok = TRUE ;
                } 
              else /* interesting case, create a new mrna */
                {
                  int jj2 = arrayMax (smrnas) ;
                  { /* first restrict current mrna to ks1 */
                    for (i = 0 ; i < keySetMax (ks1) ; i++)
                      keySet (smrna->clones, i) = keySet (ks1, i) ;
                    keySetMax (smrna->clones) = i ; 
                  }
                  { /* now duplicate smrna */
                    SMRNA *smrna2 = arrayp (smrnas, jj2, SMRNA) ;

		    smrna = arrp (smrnas, ii, SMRNA) ;
                    *smrna2 = *smrna ; /* a1 a2 etc */
                    if (smrna->clones) smrna2->clones = arrayHandleCopy (ks2, s2m->h) ;
                    if (smrna->hits) smrna2->hits = arrayHandleCopy (smrna->hits, s2m->h) ; 
		    if (smrna->estHits) smrna2->estHits = arrayHandleCopy (smrna->estHits, s2m->h) ; 
		    if (smrna->dnas) smrna2->dnas = arrayHandleCopy (smrna->dnas, s2m->h) ; 
		    if (smrna->orfs) smrna2->orfs = arrayHandleCopy (smrna->orfs, s2m->h) ; 
                  }
                  mrnaGrignotte (s2m, sc, gmrna, smrnas, ii, FALSE) ;
                  mrnaGrignotte (s2m, sc, gmrna, smrnas, jj2, FALSE) ;

                  ok = TRUE ;
                  keySetDestroy (kss) ;  
                  keySetDestroy (ks1) ;
                  keySetDestroy (ks2) ;
                  goto allgapAgain ;
                }
            }
          keySetDestroy (ks1) ;
          keySetDestroy (ks2) ;
          goto xGapAgain ;
        }
      keySetDestroy (ks1) ;
      keySetDestroy (ks2) ;

      /* cleanup */
      for (i = j = 0, vp = wp =arrp (smrna->hits, 0, HIT) ; i < mx ; i++, vp++)
        if (vp->type)
          {
            if (j < i) { *wp = *vp ; }
            j++ ; wp++ ;
          }
      arrayMax (smrna->hits) = j ;          
    } 

  /* remove dead mrnas */
  for (ii = 0, jj = 0 ; ii < arrayMax(smrnas) ; ii++)
    {
      smrna = arrp (smrnas, ii, SMRNA) ;
      mx = keySetMax(smrna->clones) ;
      if (!mx)
        continue ;
      if (jj != ii)
        {
          SMRNA *smrna2 = arrp (smrnas, jj, SMRNA) ;
          *smrna2 = *smrna ;
        }
      jj++ ;
    }
  arrayMax (smrnas) = jj ;

  /* identify newly created gaps in the mrnas 
     this would be useful in ank.7 (ZH18)
  */
 done:
  keySetDestroy (ks) ;
  return ok ;
} /* mrnaRemoveRedundant */

/**********************************************************************************/
/**********************************************************************************/
/* flag all exons introns not reused inside the conserved mrnas */
static void makeMrnaFilterGeneHits (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas)
{
  int ii, jj, k = 0 ;
  HIT *up, *vp ;
  Array hits = arrayHandleCreate (300, HIT, s2m->h) ;
  unsigned int badType, myBadType ;
  SMRNA *smrna ;
  
  /* accumulate all non discarded mrna hits in single list */
  for (ii = 0 ; ii < arrayMax(smrnas) ; ii++)
    {
      smrna = arrp (smrnas, ii, SMRNA) ;
      if (smrna->hits && smrna->bestDna >= 0 && smrna->cGroup < 2)
        for (jj = 0 ; jj < arrayMax (smrna->hits) ; jj++)
          {
            vp = arrp (smrna->hits, jj, HIT) ;
            if (!vp->type)
              continue ;
            up = arrayp (hits, k++, HIT) ;
            up->type = vp->type & ~gB ;
            up->a1 = vp->a1 + smrna->a1 - 1 ;
            up->a2 = vp->a2 + smrna->a1 - 1 ;
          }
    }
  
  cDNASwapA (hits) ; 
  arraySort (hits, cDNAOrderByA1) ;
  arrayCompress (hits) ; 
  
  /* simplify the stolen flags */
  if (arrayMax (hits))
    for (ii = 0, up = arrp (hits, 0, HIT) ; ii < arrayMax (hits) ; up++, ii++)
      {
        if (!up->type)
          continue ;
        /* find the most exquisite flag */
        badType = gOpenGap | gStolen | gPredicted ;
        for (vp = up + 1, jj = ii + 1 ; jj < arrayMax (hits) && 
               vp->a1 == up->a1 && vp->a2 == up->a2; vp++, jj++)
          {
            if (!vp->type)
              continue ;
            if ((up->type & badType) || (vp->type & badType))
              {
                if (!(up->type & badType) || !(vp->type & badType))
                  { up->type &= ~badType ;  vp->type &= ~badType ; }
                else if ((up->type & gOpenGap) || (vp->type & gOpenGap))
                  { 
                    up->type |= gOpenGap ; 
                    up->type &= ~gStolen ; 
                    up->type &= ~gPredicted ;
                    vp->type |= gOpenGap ; 
                    vp->type &= ~gStolen ; 
                    vp->type &= ~gPredicted ;
                  }
                else if ((up->type & gStolen) || (vp->type & gStolen))
                  { 
                    up->type |= gStolen ; 
                    up->type &= ~gPredicted ;
                    vp->type |= gStolen ; 
                    vp->type &= ~gPredicted ;
                  }
              }
          }
        myBadType = badType & up->type ;          
        for (vp = up + 1, jj = ii + 1 ; jj < arrayMax (hits) && 
               vp->a1 == up->a1 && vp->a2 == up->a2; vp++, jj++)
          {
            if (!vp->type)
              continue ;
            vp->type &= ~badType ;
            vp->type |= myBadType ;
            if (up->type == vp->type)
              vp->type = 0 ;
          }
      }

  /* register happy few */ 
  if (arrayMax (hits))
    {
      for (ii = jj = 0, up = vp = arrp (hits, 0, HIT) ; ii < arrayMax (hits) ; up++, ii++)
      {
        if (up->type)
          {
            if (vp != up) { *vp = *up ; }                  
            vp++ ; jj++ ;
          }
      }
      arrayMax(hits) = jj ;
    }

      /* mark the inside alternatives */
  /* introns are alternative only if they overlap a different intron
   *  but cassette exons count as alternative exons 
   */
   if (arrayMax (hits))
     {
       for (ii = 0, up = arrp (hits, 0, HIT) ; ii < arrayMax (hits) ; up++, ii++)
         {
           for (vp = up + 1, jj = ii + 1 ; jj < arrayMax (hits) && vp->a1 <= up->a2 ; vp++, jj++)
             {
	       if (up->type & gI)
		 {
		   if (vp->type & gI)
		     up->type |= gB ; 
		 }
	       else
		 up->type |= gB ; 
	       if (vp->type & gI)
		 {
		   if (up->type & gI)
		     vp->type |= gB ; 
		 }
	       else
		 vp->type |= gB ; 
             }
         }
     }                  

  arrayDestroy (gmrna->hits) ;
  gmrna->hits = hits ;
}

/**********************************************************************************/
/**********************************************************************************/
/* keep only a resonable number of best models */
/* set smrna->bestDna = -1 to remove it */
/*        query  find mrna  NOT Bad_Quality && ct_ac > 1 */
/*         edit Bad_quality  */
/*         query  find mrna   */
/*         query NOT Bad_Quality && Longest_ORF < 300 && (!gt_ag || Gap_length) && !(COUNT cdna_clone > 5)  */
/*         edit Bad_quality */
/**********************************************************************************/

static BOOL neverDiscard (KEY clone)
{
  KEYSET ks ;
  BOOL ok = FALSE ;

  ks = queryKey (clone, ">Read ; Ref_mrna || Ref_seq") ;
  if (keySetMax(ks))
    ok = TRUE ;
  keySetDestroy (ks) ;
  
  return ok ;  
}   /* neverDiscard  */

/**********************************************************************************/

static void makeMrnaFilterGene (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas)
{
  int ln, i, ii, iGoodClones = 0, cds ;
  int minLength = mrnaPleaseMinTranscriptSize ();
  BOOL doDiscard = FALSE ;
  SMRNA *smrna ;
  KEYSET goodClones = 0, ks = 0 ;
  ORFS * orf ;
  
  if (minLength < 0) 
    { minLength *= -1 ; doDiscard = TRUE ; }
  if (doDiscard)
    goodClones = arrayHandleCreate (keySetMax(sc->sh.clones), KEY, s2m->h) ;
  /* filter on ORF length */
  
  for (ii = 0 ; ii < arrayMax(smrnas) ; ii++)
    {
      smrna = arrp (smrnas, ii, SMRNA) ;
      if (smrna->bestDna < 0  || !smrna->orfs || !arrayMax(smrna->orfs))
        continue ;
      for (cds = i = 0, orf = arrp (smrna->orfs, 0, ORFS) ; 
           i < arrayMax (smrna->orfs) ; orf++, i++)
        if (cds < orf->cds)
          cds = orf->cds ;

      if (cds < 3*minLength &&          /* short orf */
          keySetMax (smrna->clones) < 6)   /* not many clones */
        {
          int ngt_ag = 0, u1, u2 ;
          HIT *up = arrayMax (smrna->hits) ? arrp (smrna->hits, 0, HIT) : 0 ;
          KEY foot ;
          
          if (up)
            for (i = 0 ; i < arrayMax(smrna->hits) ; i++, up++)
              {
                if (sc->isUp)
                  { u1 = sc->a1 - up->a1 + 1 - smrna->a1 + 1 ; u2 = sc->a1 - up->a2 + 1  - smrna->a1 + 1 ; }
                else
                  { u1 = sc->a1 + up->a1 - 1 + smrna->a1 - 1 ; u2 = sc->a1 + up->a2 - 1 + smrna->a1 - 1 ; }
                
                 ln = up->a2 - up->a1 + 1 ;
                 if (ln > 0 && (up->type & gI)  && !(up->type & gJ))
                  {
                    foot =  mrnaIntronBoundary (s2m->dnaD, s2m->dnaR, u1, u2) ;
                    if (foot == _gt_ag || foot == _gc_ag)
                      ngt_ag++ ;
                    else if (foot == _ct_ac)
                      ngt_ag -= 2 ;
                  }
                if ((up->type & gJ) && i > 1 && ((up-2)->type & gJ)) ngt_ag -= 6 ;
                /*  else if ((up->type & gGap) || (up->type & gJ))
                 *   nIntrons -= 2 ; 
                 */
              }
          if (ngt_ag <= 0 || 
              ( 0 && ngt_ag == 1 &&  keySetMax (smrna->clones) < 3 && 5*cds < 4*3*minLength))
            smrna->cGroup = (doDiscard ? 2 : 1) ;
        }
      if (doDiscard && smrna->cGroup < 2)
        for (i = 0 ; i < keySetMax(smrna->clones) ; i++)
          keySet (goodClones, iGoodClones++) = keySet (smrna->clones, i) ;
    }
  if (doDiscard)
    {
      KEY clone ;
      for (i = 0 ; i < keySetMax(sc->sh.clones) ; i++)
        {
          clone = keySet (sc->sh.clones, i) ;
          if (!strncmp(name(clone), "NM_",3) || neverDiscard (clone))
            keySet (goodClones, iGoodClones++) = clone ;
        }
      keySetSort (goodClones) ;
      keySetCompress (goodClones) ;
      ks = keySetMINUS (sc->sh.clones, goodClones) ;
      keySetDestroy (sc->sh.clones) ;
      sc->sh.clones = goodClones ;
      if (sc->sh.ignoredClones)
        {
          KEYSET ks1 = ks ; 
          ks = keySetOR (sc->sh.ignoredClones, ks1) ;
          keySetDestroy (sc->sh.ignoredClones) ;
          
          keySetDestroy (ks1) ;
        }
      sc->sh.ignoredClones = arrayHandleCopy (ks, s2m->h) ;
      keySetDestroy (ks) ;
    }
}  /* makeMrnaFilterGene */

/**********************************************************************************/
/**********************************************************************************/

void showSmrnas (S2M *s2m, Array smrnas, char *mm)
{
  int ii, j, j1 ;
  SMRNA *smrna ;
  
  printf("smrnas:: %s %s\n", mm ? "searching for cDNA_clone" : "", mm && *mm ? mm : "") ;
  if (smrnas && arrayMax(smrnas))
    for (ii = j1 = 0 ; ii < arrayMax(smrnas) ; ii++)
      {
        smrna = arrp (smrnas, ii, SMRNA) ;
        if (smrna->clones)
          for (j = j1 = 0 ; j < keySetMax(smrna->clones) ; j++)
            {
              if (mm && mm != (char*)1 && *mm && !strstr (name(keySet(smrna->clones,j)), mm))
                continue ;
              if (j1 % 10 == 9)
                printf( "\n\t") ;
              if (mm == (char*)1 && j > 5) break ;
              printf("  %s", name(keySet(smrna->clones,j))) ;
              j1++ ;
            }
        if (!smrna->hits && !j1) continue ;
        printf ("\nsmrna %03d\ta1=%d\ta2=%d:", ii, smrna->a1, smrna->a2) ;
        printf("\n") ;
        if (smrna->orfs && arrayMax(smrna->orfs))
          printf ("Longest ORF %d %d %d\n"
                  , smrna->orfs && arrayMax(smrna->orfs) > 0 ?  arr (smrna->orfs, 0, ORFS).nOpen : 0
                  , smrna->orfs && arrayMax(smrna->orfs) > 1 ?  arr (smrna->orfs, 1, ORFS).nOpen : 0
                  , smrna->orfs && arrayMax(smrna->orfs) > 2 ?  arr (smrna->orfs, 2, ORFS).nOpen : 0
                  ) ;
        if (0 && mm != (char*)1) showHits(smrna->hits) ;        
        printf("\n") ;
      }        
  printf("\nmaxgap = %d\n", mrnaCountGaps (smrnas)) ;
}

/**********************************************************************************/

KEY makeMrnaGene (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas, 
                  KEYSET clipTops, KEYSET clipEnds, Array linkPos)
{
  int ii, nn, nLoop, ii1, ii2 ;
  KEY mycosmid = 0, mycosmid1 = 0, geneBox = 0 ;
  SMRNA *smrna ;
  BOOL useAm, debug = FALSE ;
  KEYSET genes = 0 ;
  int doFilter = mrnaPleaseMinTranscriptSize() ; 
  int bigGap = mrnaPleaseBigGapSize() ; 
  static int nnG = 0 ;  
  if (0) messalloccheck () ;
  chrono ("makeMrnaGene") ;

  if (! ((++nnG) % 200))
    sessionClose (TRUE) ;
  cdnaReExtendHits (s2m, sc, gmrna, clipTops, clipEnds) ;
  mrnaGrignotte (s2m, sc, gmrna, smrnas, -1, FALSE) ;  /* grignotte all */
  chronoReturn () ;

  chrono ("makeMrnaGeneRedundant") ; 

  for (nn = nLoop = 0 ; (nLoop < 2 || nn) && nLoop < 5 ; nLoop++)
    {
      nn = 0 ;
      if (0) showSmrnas (s2m, smrnas, 0) ;
      if (1 || mrnaPleaseStealFromPrediction())
        nn += mrnaNewCutIfGapInMrna (s2m, sc, gmrna, smrnas, bigGap) ;
      if (nn) mrnaGrignotte (s2m, sc, gmrna, smrnas, -1, FALSE) ;  /* grignotte all, needed a second time  */         
      if (0) showSmrnas (s2m, smrnas, 0) ;
      nn += mrnaPreferGaplessMrna (s2m, sc, gmrna, smrnas, FALSE) ;
      nn += mrnaRemoveRedundant (s2m, sc, gmrna, smrnas) ; /* cut out non specific contigs */
      if (1)
        {
	  nn += mrnaPreferGaplessMrna (s2m, sc, gmrna, smrnas, TRUE) ;
          nn += mrnaRemoveRedundant (s2m, sc, gmrna, smrnas) ; /* cut out non specific contigs */
        }
      if (0) showSmrnas (s2m, smrnas, 0) ;
      if ((!nLoop && bigGap) || nn) /* need at least 1 pass to cut the BIG clone pairs */ 
	mrnaGrignotte (s2m, sc, gmrna, smrnas, -1, bigGap) ;  /* grignotte all, needed a second time  */ 
      if (0) showSmrnas (s2m, smrnas, 0) ;
    }
  chronoReturn () ;
  chrono ("makeMrnaGeneSteal") ;
  if (mrnaPleaseStealFromPrediction())
    { ii1 = 0 ; ii2 = 3 ; } /* 0,1: from variant, 2 from predictions */
  else
    { ii1 = 1 ; ii2 = 2 ; } /* we use to say  ii1 = 0 ; ii2 = 2 ; to allow steal from neighbours */
  for (ii =ii1 ; ii < ii2 ; ii++)
    mrnaStealFromAlternative (s2m, sc, gmrna, smrnas, ii, linkPos) ;
  mrnaFillLastGap (smrnas, 100) ; /* 2009_04_06,  smmal gaps are frequent in deep sequencing experiments */

  /* loop because until now gaps between 5p 3p read pairs 
   * have remained incompatible with introns
   */
  while (mrnaFuseRedundant (s2m, gmrna, smrnas)) ; /* fuse mrna identical because of stealing*/

  if (0) showSmrnas (s2m, smrnas, 0) ;

  s2m->composite = isComposite = TRUE ;
  if (s2m->composite)
    mrnaDesignUsingCompositeStrategy (s2m, sc, gmrna, smrnas) ;

  /* keep only a reasonable number of best models */
  /* set smrna->bestDna = -1 to remove it */
  if (! s2m->composite)
    {
      mrnaSetCompletenessFlags (s2m, sc, gmrna, smrnas) ;
      
  /* construct the ORFs */
      mrnaConstructOrfs (s2m, sc, gmrna, smrnas) ;
      arraySort (smrnas, smrnaOrderByLengthandBadQuality) ;
      if (0) showSmrnas (s2m, smrnas, 0) ;
      if (0) showOrfs (0) ; /* to please the compiler */

      if (doFilter)
	makeMrnaFilterGene (s2m, sc, gmrna, smrnas) ;
      makeMrnaFilterGeneHits (s2m, sc, gmrna, smrnas) ;
      /* save the remaining mrnas */
      if (0) showSmrnas (s2m, smrnas, 0) ;
      arraySort (smrnas, smrnaOrderByLengthandBadQuality) ;
      if (0) showSmrnas (s2m, smrnas, 0) ;

      /* unfortunate hack, selectCosmid alters sc->s2m because it was conceived for a single sc contig */
      mycosmid =  sc->s2m->cosmid ;
      mycosmid1 =  sc->s2m->cosmid1 ;
      
      if (1) mrnaSelectCosmid (sc) ;  /* attributes gene to best adequate cosmid */
      sc->gene = mrnaGetGeneName (sc) ;

      
      if (doFilter == 0 || arrayMax(smrnas))   /* could be zero because of filtergene */
	{
	  if (0) /* that function is not yet reliable */
	    mrnaCountErrors (s2m, sc, gmrna->estHits) ;
	  mrnaSaveGene (s2m, sc, gmrna, smrnas) ;
	  if ((geneBox = keyGetKey (sc->gene, str2tag("Gene"))) &&
	      keyFindTag (geneBox, str2tag("Use_AM")))
	    useAm = TRUE ;
	  else
	    useAm = FALSE ;
	  for (ii = 0 ; ii < arrayMax(smrnas) ; ii++)
	    {
	      BOOL addOneSuccess ;
	      
	      smrna = arrp (smrnas, ii, SMRNA) ;
	      if (! smrna->hits)
		continue ;
	      if (!(smrna->bestDna >= 0 && smrna->cGroup < 2))
		continue ;

	      mrnaSaveMrna (s2m, sc, gmrna->estHits, smrnas, smrna, FALSE) ;
	      mrnaSaveProductCoverage (s2m, sc, gmrna, smrnas, smrna) ;
	      
	      mrnaAddKantorInfo (smrna->tr) ;
	      addOneSuccess = mrnaAddOneTiling (sc, smrna->tr) ;
	      
	      if (useAm)
		{
		  if (addOneSuccess)
		    {
		      mrnaTransferRefseqMaker2Product (keyGetKey(smrna->tr, str2tag("RefSeqMaker"))) ;
		      mrnaSaveProductCoverage (s2m, sc, gmrna, smrnas, smrna) ;
		      mrnaAddKantorInfo (smrna->tr) ; 
		      mrnaAddOneTiling (sc, smrna->tr) ; /* after kantor info */
		    }
		  else
		    if (debug) messout ("Cannot construct mrna %s from the cDNA sequence because of a tiling error, sorry",
					name (smrna->tr)) ;
		}
	    }
	  mrnaMrnaScore (smrna->tr) ;	       /* after tiling */
	}
      mrnaSaveMrnaHierarchy (sc->gene) ;
      if (mrnaPleaseKillNonBest ())
	mrnaSaveGeneKillNonBest (s2m, sc, gmrna, smrnas) ;
      mrnaSaveMicroIntrons (s2m, sc, gmrna, smrnas) ;
      mrnaSaveProductHierarchy (s2m, sc, gmrna, smrnas) ;
      if (1 || ! s2m->compositeDesignCovering)
	abiFixLabelPolyATg (sc->gene, 0, 0, 0) ;
      /* this code should work on the genen for tg and/or pg */
      genes = queryKey (sc->gene, ">Gene") ;
      if (mrnaPleaseStealFromPrediction())
        mrnaSteal3pFromPrediction (genes) ;
      mrnaSaveAntisens (genes) ;
      mrnaAddMrnas (genes) ;
      mrnaAddAlterSpliceDetails (genes) ;
      mrnaCountClones (genes, 0, 0) ;

      giwAddKeysetIntronClass (0, sc->gene) ;
      giwAddKeysetProbeWalls (0, sc->gene) ;
      if (!getPleaseNoKantorInfo ())
	{
	  mrnaSaveAceKog (genes) ;
	  mrnaSavePfam2Go (genes) ;
	  mrnaSavePsort2Go (genes) ;
	  mrnaSaveGenePastilles (genes) ;
	  mrnaSaveGeneTitle (genes) ;
	}
      keySetDestroy (genes) ;

      sc->s2m->cosmid = mycosmid ;
      sc->s2m->cosmid1 = mycosmid1 ;
    }
  
  chronoReturn () ; 

  return sc->gene ;
} /* makeMrnaGene */

/**********************************************************************************/
/**********************************************************************************/

static KEYSET mrnaSplitOneDoubleGene (KEY tg, KEYSET pgs, KEYSET gids) 
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j,  g ;
  int gMax = arrayMax (gids) ;
  KEY _Gene_part = str2tag ("Gene_part") ;
  KEYSET genes = keySetCreate () ;
  KEYSET mrnas = queryKey (tg, ">mRNA") ;
  int mMax = arrayMax (mrnas) ;
  Array aa = arrayHandleCreate (1+mMax, HIT, h) ;
  vTXT txt = vtxtHandleCreate (h) ;
  OBJ Mrna  = 0 ;

  if (gMax < 8)
    for (i = 0 ; i < mMax ; i++)
      {  /* associate each mRNA with a best score to its gids via a flag in up->zone */
	KEY gene ;
	KEY mrna = keySet (mrnas, i) ;
	HIT *up = arrayp (aa, i, HIT) ;
	up->gene = mrna ;
	KEYSET ks = queryKey (mrna, "{>Identical_to_genefinder} SETOR {>Matching_genefinder} SETOR {>Touching_genefinder} ; >geneid_pg") ;
	int n = 0 ;
	
	Mrna = bsUpdate (mrna) ;
	vtxtClear (txt) ;
	for (j = 0 ; j < keySetMax (ks) ; j++)
	  {
	    KEY k1 = keySet (ks, j) ;
	    for (g = 0 ; g < gMax ; g++)
	      if (k1 == keySet (gids, g))
		{
		  up->zone |= (1 << g) ;
		  if(n++)
		    vtxtPrint (txt, "--") ;
		  vtxtPrint (txt, name (k1)) ;
		  bsAddKey (Mrna, _GeneId, k1) ;
		}
	  }
	keySetDestroy (ks) ;
	
	if (n)
	  {
	    if (bsFindTag (Mrna, _Gene_part))
	      bsRemove (Mrna) ;
	    lexaddkey (vtxtPtr (txt), &gene, _VGene) ;
	    bsAddKey (Mrna, _Gene_part, gene) ;
	    keySetInsert (genes, gene) ; 
	  }
	bsSave (Mrna) ;
	
	if (n)
	  {
	    int g ;
	    OBJ Gene = bsUpdate (gene) ;
	    
	    for (g = 0 ; g < gMax ; g++)
	      if (up->zone & (1 << g))
		{
		  KEY gid = keySet (gids, g) ;
		  bsAddKey (Gene, _GeneId, gid) ;
		}
	    bsSave (Gene) ;
	  }
      }
  
  ac_free (h) ;
  return genes ;
} /*  mrnaSplitOneDoubleGene */

/**********************************************************************************/
/* merge weak genes in better ones */
static BOOL  mrnaSplitOneMerge(KEYSET genes)
{
  BOOL ok = FALSE ;
  int ii, jj ;
  KEY _GeneId = str2tag ("GeneId") ;
  KEY _mRNA_part = str2tag ("mRNA_part") ;
  KEYSET gids1 = 0 ;
  KEYSET gids2 = 0 ;
  KEYSET gids3 = 0 ;

  for (ii = 0 ; ii < keySetMax (genes) ; ii++)
    {
      KEY gene1 = keySet (genes, ii) ;
      gids1 = queryKey (gene1, ">geneId") ;
      if (keySetMax (gids1))
	for (jj = 0 ; jj < keySetMax (genes) ; jj++)
	  if (ii != jj)
	    {
	      KEY gene2 = keySet (genes, jj) ;
	      gids2 = queryKey (gene2, ">geneId") ;
	      
	      gids3 = keySetMINUS (gids1, gids2) ;
	      if (! keySetMax (gids3)) /* we could merge gene1 into gene2 */
		{  /* look for good matches only */
		  KEYSET w1 = queryKey (gene1, ">mRNA_part ; {>Identical_to_genefinder} SETOR {>Matching_genefinder} ; >geneid_pg") ;
		  /*
		    KEYSET w2 = queryKey (gene2, ">mRNA_part ; {>Identical_to_genefinder} SETOR {>Matching_genefinder} ; >geneid_pg") ;
		    KEYSET w3 = keySetXOR (w1, w2) ;
		  */
		  if (! keySetMax (w1) /* || ! keySetMax (w3) */)  /* i do not have a strong case, or we have the same ones */ 
		    {
		      int i ;
		      KEYSET mrnas = queryKey (gene1, ">mRNA_part") ;
		      OBJ Gene = bsUpdate (gene1) ;
		      bsKill(Gene) ;	

		      ok = TRUE ;

		      Gene = bsUpdate (gene2) ;
		      for (i = 0 ; i < keySetMax (gids1) ; i++)
			bsAddKey (Gene, _GeneId, keySet (gids1, i)) ;
		      for (i = 0 ; i < keySetMax (mrnas) ; i++)
			bsAddKey (Gene, _mRNA_part, keySet (mrnas, i)) ;
		      bsSave (Gene) ;

		      keySetDestroy (mrnas) ;
		    }
		  keySetDestroy (w1) ;
		  /*
		    keySetDestroy (w2) ;
		    keySetDestroy (w3) ;
		  */
		  if (ok)
		    goto done ;
		}
	      
	      keySetDestroy (gids2) ;
	      keySetDestroy (gids3) ;
	    }
      keySetDestroy (gids1) ;
    }

 done:
    keySetDestroy (gids1) ;
    keySetDestroy (gids2) ;
    keySetDestroy (gids3) ;
  return ok ;
} /* mrnaSplitOneMerge */

/**********************************************************************************/
/* set the geometry */
static void mrnaSplitOneCoords (KEYSET genes)
{
  int ii ;
  KEY _mRNAs =  str2tag("mRNAs") ;

  for (ii = 0 ; ii < keySetMax (genes) ; ii++)
    {
      KEY mrna = 0, gSeq = 0, gChrom = 0, gene = keySet (genes, ii) ;
      int g1, g2 = -1 ;
      int a1, a2 ;
      OBJ Gene = bsUpdate (gene) ;
      KEYSET mrnas = queryKey (gene, ">mRNA_part") ;
      KEYSET rids = queryKey (gene, ">mRNA_part ; >matching_short_RNA ; >geneid_pg") ;
      int jj ; 

      for (jj = 0 ; jj < keySetMax (rids) ; jj++)
	bsAddKey (Gene, _GeneId, keySet (rids, jj)) ;

      for (jj = 0 ; jj < keySetMax (mrnas) ; jj++)
	{
	  KEY chrom = 0 ; 
	  OBJ Mrna = 0 ;

	  mrna = keySet (mrnas, jj) ;
	  Mrna = bsCreate (mrna) ;

	  bsGetKey (Mrna, _Genomic_sequence, &gSeq) ;
	  if (bsGetKey (Mrna, _IntMap, &chrom) &&
	    bsGetData (Mrna, _bsRight, _Int, &a1)
	    )
	  bsGetData (Mrna, _bsRight, _Int, &a2) ;
	  if (g2 == -1)
	    { g1 = a1 ; g2 = a2 ; gChrom = chrom ; }
	  if (a1 < a2)
	    { if (g1 > a1) g1 = a1 ; if (g2 < a2) g2 = a2 ; }
	  else
	    { if (g1 < a1) g1 = a1 ; if (g2 > a2) g2 = a2 ; }
	}
      if (gChrom) 
	{
	  bsAddKey (Gene, _IntMap, gChrom) ;
	  bsAddData (Gene, _bsRight, _Int, &g1) ;
	  bsAddData (Gene, _bsRight, _Int, &g2) ;
	}
      
      if (1)
	{  /* connect the GenePart (i.e. all 'genes') via the Tg, to the main GeneBox */
	  KEY tg = mrna ? keyGetKey (mrna, _From_gene) : 0 ;
	  KEY gBox = tg ? keyGetKey (tg, _Gene) : 0 ;
	  if (gBox)
	    bsAddKey (Gene, str2tag("GeneBox"), gBox) ;
	}
      bsSave (Gene) ;
      
      /* deal with the geometry */
      if (gSeq)
	{
	  OBJ Gseq = bsUpdate (gSeq) ;
	  int u1, u2 ;
	  int da1, da2 ;
	  if (bsFindKey (Gseq, _mRNAs, mrna) &&
	      bsGetData (Gseq, _bsRight, _Int, &u1) &&
	      bsGetData (Gseq, _bsRight, _Int, &u2)
	      )
	    {
	      if (a1 < a2)
		{ da1 = g1 - a1 ; da2 = g2 - a2 ; }
	      else
		{ da1 = a1 - g1 ; da2 = a2 - g2 ; }
	      if (u1 < u2)
		{ u1 += da1 ; u2 += da2 ; }
	      else
		{ u1 -= da1 ; u2 -= da2 ; }
	      bsAddKey (Gseq, _Genes, gene) ;
	      bsAddData (Gseq, _bsRight, _Int, &u1) ;
	      bsAddData (Gseq, _bsRight, _Int, &u2) ;
	    }
	  bsSave (Gseq) ;
	}

      keySetDestroy (mrnas) ;
      keySetDestroy (rids) ;
    }
} /* mrnaSplitOneCoords */

/**********************************************************************************/
/* set the products */
static void mrnaSplitOneProducts (KEYSET genes)
{
  AC_HANDLE h = ac_new_handle () ;  
  KEY _GeneBox =  str2tag("GeneBox") ;
  int ii ;
  
  for (ii = 0 ; ii < keySetMax (genes) ; ii++)
    {
      OBJ Gene = 0 ;
      KEY gene = keySet (genes, ii) ;
      KEYSET products = queryKey (gene, ">mRNA_part ; >Product") ;
      KEYSET gids = queryKey (gene, ">GeneId") ;
      int jj ;

      for (jj = 0 ; jj < keySetMax (products) ; jj++)
	{
	  KEY product = keySet (products, jj) ;
	  OBJ Product = bsUpdate (product) ;
	  if (bsFindTag (Product, _GeneBox))
	    bsRemove (Product) ;
	  bsAddKey (Product, _GeneBox, gene) ;
	  bsSave (Product) ;
	}
      keySetDestroy (products) ;  
      
      Gene = bsUpdate (gene) ;
      for (jj = 0 ; jj < keySetMax (gids) ; jj++)
	{
	  KEY loc, gid = keySet (gids, jj) ;
	  bsAddKey (Gene, _GeneId, gid) ;
	  if ((loc = keyGetKey (gid, _LocusLink)))
	    bsAddKey (Gene, _LocusLink, loc) ;
	}
      bsSave (Gene) ;
      keySetDestroy (gids) ;  
    }

  ac_free (h) ;
} /*  mrnaSplitOneProducts */

/**********************************************************************************/
/* set the products */
static void mrnaSplitOnePg (KEYSET genes)
{
  AC_HANDLE h = ac_new_handle () ;  
  KEY _Genefinder =  str2tag("Genefinder") ;
  int ii ;
  
  for (ii = 0 ; ii < keySetMax (genes) ; ii++)
    {
      OBJ Gene = 0 ;
      KEY gene = keySet (genes, ii) ;
      KEYSET pgs = queryKey (gene, ">mRNA_part ; {>Matching_short_RNA} SETOR {>identical_to_genefinder} SETOR {>Matching_genefinder} SETOR {>Touching_genefinder}") ;
      int jj ;

      Gene = bsUpdate (gene) ;
      for (jj = 0 ; jj < keySetMax (pgs) ; jj++)
	{
	  KEY pg = keySet (pgs, jj) ;
	  bsAddKey (Gene, _Genefinder, pg) ;
	}
      bsSave (Gene) ;
      keySetDestroy (pgs) ;
    }

  ac_free (h) ;
} /*  mrnaSplitOnePg */

/**********************************************************************************/

void mrnaSplitDoubleGenes (KEYSET tgs)
{
  int i, iMax = keySetMax (tgs) ;

  for (i = 0 ; i < iMax ; i++)
    {
      KEY tg = keySet (tgs, i) ;
      
      KEYSET pgs = queryKey (tg, ">Matching_genefinder_gene") ;
      KEYSET gids = pgs ? query (pgs, ">GeneId_pg") : 0 ;

      if (gids && keySetMax (gids) > 1) 
	{
	  KEYSET subGenes = mrnaSplitOneDoubleGene (tg, pgs, gids) ;
	  while (mrnaSplitOneMerge (subGenes))
	    {} 
	  mrnaSplitOneCoords (subGenes) ;
	  mrnaSplitOneProducts (subGenes) ;
	  mrnaSplitOnePg (subGenes) ;
	  mrnaAddAlterSpliceDetails (subGenes) ; 
	  mrnaCountClones (subGenes, 0, 0) ;

	  giwAddKeysetIntronClass (subGenes, 0) ;
	  giwAddKeysetProbeWalls  (subGenes, 0) ;
	  
	  mrnaSaveAceKog (subGenes) ;
	  mrnaSavePfam2Go (subGenes) ;
	  mrnaSavePsort2Go (subGenes) ;
	  
	  mrnaSaveGenePastilles (subGenes) ;

	  keySetDestroy (subGenes) ;
	}
      else
	{
	  KEY gene = keyGetKey (tg, _Gene) ;
	  if (gene)
	    {
	      KEYSET rids = queryKey (tg, ">mRNA;>matching_short_RNA;>geneid_pg") ;
	      if (keySetMax (rids))
		{
		  int i ;
		  OBJ Gene = bsUpdate (gene) ;
		  for (i = 0 ; i < keySetMax (rids) ; i++)
		    {
		      KEY gid = keySet (rids, i) ;
		      bsAddKey (Gene, _GeneId, gid) ;
		    }
		  bsSave (Gene) ;
		}
	      keySetDestroy (rids) ;
	    }
	}
      
      keySetDestroy (pgs) ;
      keySetDestroy (gids) ;
    }

  return ;
} /* mrnaSplitDoubleGenes */

/**********************************************************************************/
/**********************************************************************************/

