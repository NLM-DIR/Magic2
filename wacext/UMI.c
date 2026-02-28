/*  File: UMI.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2005
 * -------------------------------------------------------------------
 * AceView is free software; you can redistribute it and/or
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
 * This file is part of the AceView project developped by
 *	 D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *
 *  This should be linked with the acedb libraries
 *  available from http://www.acedb.org
 *
 *  Please direct all question about this code to
 *  mieg@ncbi.nlm.nih.gov

 * 2017_09_07
 * strategic idea

 * Analyze the UMI alignments, the bar code is in column 20
 * Rationalize away the multiplets: same a1 coordinates and same barcode
 *    while creating the consensus of the match errors and overlaps
 * Count expression per gene and transcripts in all target classes
 *    add the hierarchic counts to detect the missing transcripts of the top hierarchy
 * Construct count and export all variations 
 * Locate all gene fusion mrna1, mrna2, coords, delta-x
 *   export the delta-x histogram
 * Construct the consesnus (bridge) for all SNPs and gene fusions
 * Rextend the partial reads and count the 'validations'

 ** copied from variant-caller.c (to be modified)
 * 1: Read the hits file, select a chromosomal zone of say 10 Mb
 *    Consider only the extremal coords of each read ali
      Construct the coverage wiggle
      Construct the wiggle of all brks: snp, overhang, deletion

 * 2: Delineate the segments a--------X==X==X---------b
 *    such that X is a brk, x==X is shorter than DELTA a---X and X----b are longer
 *    and X is in a place covered >= 10 and #brk/#cover > 10%

 * 3: Read the fastc file of the reads touching one of the segments

 * 4: Loop on all segments, get the relevant ali, sort them
      Anchor at a, and extend base per base
        each time there is minor letter > 5%, biffurcate

 * 5: Merge the left and right consensus if both exist
      Export the consensus either as a del-ins or as an overhang

 *    DONE

 */

#define VERSION "1.1"

#include "ac.h"
#include "channel.h"
#include "bitset.h"

/* snp chaining distance */
#define DELTA 12
#define ABS(_x) ((_x) < 0 ? -(_x) : (_x))

static char *types[] = {"Wild_type", "Substitution", "Duplication", "Insertion", "Diminution", "Deletion",  "DelIns", "Intron gt_ag forward", "Intron gc_ag forward", "Intron at_ac forward", "Intron gt_ag reverse", "Intron gc_ag reverse", "Intron at_ac reverse", "Unkown", 0} ;

typedef struct UMIStruct {
  AC_HANDLE h ;
  const char *outFileName ;
  BOOL gzi, gzo, debug, is64 ;
  ACEIN ai ; 
  ACEOUT ao ;

  const char *run ;
  const char *lane ;
  const char *method ;
  const char *laneListFileName ;

  const char *fastcDir ;
  const char *hitDir ;

  const char *target_class ;
  const char *target ;
  const char *targetFileName ;
  int offset, t1, t2 ;

  BOOL intron_only ;
  BOOL SAM ;
  BOOL geneFusion ;
  const char *SAM_file ;
  const char *hitFile ;
  const char *geneMapFileName ;
  const char *fastcFile ;
  const char *externalSnpListFileName ;
  const char *geneFusionFileName ;
  const char *geneFusionPositionsFileName ;

  int microDelta, Delta ;
  int minSnpCover, minSnpCount, minSnpFrequency ;
  int minOverhang ;

  int minIntron ;
  int maxIntron ;

  int max_threads ;
  int maxLanes ;

  int nErrTrackings ;
  Array lanes ;
  DICT *laneDict ;
  DICT *geneDict ;

  const char *makeTest ;
  const char *runTest ;

  Array brkIndels ;
  Array geneFusions ;
  KEYSET covers, brks, externalBrks, aliLn, fusionHisto, geneSupport ;
  BitSet externalBitSet ;
  Array externalBridges ;
  Array hits, segs, genes ;
  Array dna, dnaR ; /* target dna */

  int snpExported ;

  CHAN *getLaneInChan, *getLaneOutChan, *analyzeChan ;
  CHAN *doneChan ;
} UMI ;

typedef struct laneStruct {
  AC_HANDLE h ;
  const char *run ;
  const char *target ;
  int lane ;   /* offset in laneDict */
  Array brkIndels ;
  Array geneFusions ;
  KEYSET covers, brks, aliLn, fusionHisto, geneSupport ;
  Array clones,  hits ;
  DICT *cloneDict ;
} LANE ;

typedef struct geneFusionStruct {
  LANE *lane ;
  int clone ;  /* offset in lane->cloneDict */
  int g1, g2 ; /* genes in geneDict */
  int t1, t2 ; /* targets (i.e. mRNA or chrom) in lane->cloneDict */
  int seq, tag ;
  char wo[3] ; /* ++ +- -+ -- */ 
  BOOL isPair ;
  int score1, score2 ;
  int c1, c2, d1, d2 ; /* chain coords */
  int x1, x2, y1, y2 ; /* read coodrs */
  int a1, a2, b1, b2 ; /* target coords */
  int chromDistance ;
} GF ;

typedef struct cloneStruct {
  const char *run ;
  int clone ;  /* offset in cloneDict */
  int lane ;   /* offset in laneDict */
  Array dna1, dna2 ;
} CLONE ;

typedef struct geneStruct {
  const char *run ;
  int gene ;  /* offset in umi->geneDict */
  int chrom ;   /* offset in  umi->geneDict */
  int a1, a2 ; /* coords in gene */
  int b1, b2 ; /* non stranded coords */
} GENE ;

typedef struct hitStruct {
  int lane ;   /* offset in laneDict */
  int clone ;  /* offset in lane->cloneDict */
  int gene ;   /* offset in umi->geneDict */
  int target, target_class ; /* offset in lane->cloneDict */
  int score, mult, uu, a1, a2, c1,c2, x1, x2, z1, z2 ;  /* absolute chrom coordinates */
  Array dna, cigarettes ;
} HIT ;

typedef struct localHitStruct {
  int lane ;   /* offset in laneDict */
  int clone ;  /* offset in cloneDict */
  int mult, uu ;
  int a1, a2, x1, x2 ;   /* local coordinates (umi->t1 is substracted) */
  int b1, b2, y1, y2 ;   /* match b1 b2 of the seg */
  int c1, c2, z1, z2 ;   /* match b1 b2 in other part of same read */
  BOOL isDown, isRead1, isB ;
  Array dna, err ;
  int cc, pos ;
  int Lcc, Rcc ;
} LHIT ;

typedef struct segStruct {
  int n, a1, a2, b1, b2 ;  /* local coordinatess (umi->t1 is substracted) */
  int c1, c2, z1, z2 ;   /* match b1 b2 in other part of same read */
  int protect ;
  BOOL isTooShort ;
  BOOL isSubseg ;
  BOOL isExternal ;
  Array Fcss ;
  Array Rcss ;
  Array bridges ;
  Array subsegs ;
} SEG ;

typedef struct baseCountStruct {
  int n ;
  int base ;
} BASECOUNT ;

typedef struct brkIndelStruct {
  int pos, mult, daa, dxx, c1, c2, protect ;
} BKID ;

typedef struct consensusStruct {
  int cumulsU, cumulsNu ;
   Array dna ;
} CSS ;


typedef struct bridgeStruct {
  int ii, jj, ln, dx, ddx, ddy, ds, dt, nDiff ;
  int type ;
  int a1, a2, b1, b2 ;

  int cover, covernu ;
  int coverp, coverm ;
  int coverpnu, covermnu ;
  int mp, mm, wp, wm ;
  int mpnu, mmnu, wpnu, wmnu ;

  int zcover, zcovernu ;
  int zcoverp, zcoverm ;
  int zcoverpnu, zcovermnu ;
  int zmp, zmm, zwp, zwm ;
  int zmpnu, zmmnu, zwpnu, zwmnu ;
  BOOL isR, isExternal ;
  Array targetSnippet, tsCoded ;
  Array variantSnippet  ;
  Array reverseVariantSnippet  ;
  char buf1[201], buf2[201] ;
} BRIDGE ;

static void umiAnalyzeSeg (const UMI *umi, SEG *seg) ;

/*************************************************************************************/
 
static int umiGeneFusionOrder (const void *va, const void *vb)
{
  const GF *up = (const GF *)va, *vp = (const GF*)vb ;
  int n ;
  
  n = up->lane - vp->lane ; if (n) return n ;
  n = up->isPair - vp->isPair ; if (n) return n ;
  n = up->g1 - vp->g1 ; if (n) return n ;
  n = up->g2 - vp->g2 ; if (n) return n ;
  n = up->t1 - vp->t1 ; if (n) return n ;
  n = up->t2 - vp->t2 ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  n = up->b1 - vp->b1 ; if (n) return n ;
  n = up->b2 - vp->b2 ; if (n) return n ;
  n = up->x1 - vp->x1 ; if (n) return n ;
  n = up->x2 - vp->x2 ; if (n) return n ;
  n = up->y1 - vp->y1 ; if (n) return n ;
  n = up->y2 - vp->y2 ; if (n) return n ;
  
  n = up->wo[0] - vp->wo[0] ; if (n) return n ;
  n = up->wo[1] - vp->wo[1] ; if (n) return n ;

  return 0 ;
} /* umiCssOrder */

/*************************************************************************************/
 
static int umiCssOrder (const void *va, const void *vb)
{
  const CSS *up = (const CSS *)va, *vp = (const CSS*)vb ;
  return -up->cumulsU -up->cumulsNu + vp->cumulsU  + vp->cumulsNu ;
} /* umiCssOrder */

/*************************************************************************************/
 
static int umiBrkidOrder (const void *va, const void *vb)
{
  const BKID *up = (const BKID *)va, *vp = (const BKID*)vb ;
  int n ;
  n = up->pos -vp->pos ; if (n) return n ;
  n = up->daa - vp->daa ; if (n) return -n ;
  return 0 ;
} /* umiBrkidOrder */

/*************************************************************************************/
 
static int umiBaseCountOrder (const void *va, const void *vb)
{
  const BASECOUNT *up = (const BASECOUNT*)va, *vp = (const BASECOUNT*)vb ;
  return -up->n + vp->n ;
} /* baseOrder */

/*************************************************************************************/
 
static int umiA1Order (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;
  int n ;
  
  n = up->a1 - vp->a1 ; if (n)  return n ;
  n = up->a2 - vp->a2 ; if (n)  return n ;
  return 0 ;
} /* a1Order */

/*************************************************************************************/
 
static int umiSegOrder (const void *va, const void *vb)
{
  const SEG *up = (const SEG *)va, *vp = (const SEG *)vb ;
  int n ;
  
  n = up->a1 - vp->a1 ; if (n)  return n ;
  n = up->a2 - vp->a2 ; if (n)  return n ;
  return 0 ;
} /* segOrder */

/*************************************************************************************/
 
static int umiBridgeOrder (const void *va, const void *vb)
{
  const BRIDGE *up = (const BRIDGE *)va, *vp = (const BRIDGE *)vb ;
  int n ;
  
  if (up->type == 0) return -1 ;
  if (vp->type == 0) return 1 ;
  if (up->type < 0) return 1 ;
  if (vp->type < 0) return -1 ;

  n = up->zmm + up->zmp - vp->zmm -vp->zmp ; if (n)  return -n ; /* big cover first */
  n = up->type - vp->type ; if (n)  return n ;

  return 0 ;
} /* bridgeOrder */

/*************************************************************************************/
 
static int a1LOrder (const void *va, const void *vb)
{
  const LHIT *up = (const LHIT *)va, *vp = (const LHIT *)vb ;
  int n ;
  
  n = up->a1 - vp->a1 ; if (n)  return n ;
  n = up->a2 - vp->a2 ; if (n)  return n ;
  return 0 ;
} /* a1Order */

/*************************************************************************************/
#ifdef JUNK
static int a2Order (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;
  int n ;
  
  n = up->a2 - vp->a2 ; if (n)  return n ;
  n = up->a1 - vp->a1 ; if (n)  return n ;
  return 0 ;
} /* a2Order */
#endif
/*************************************************************************************/

static void umiShowBrkIndels (Array brkIndels)
{
  int ii, iMax = arrayMax (brkIndels) ; 
      
  for (ii = 0  ;  ii < iMax ; ii ++)
    {
      BKID *b = arrp (brkIndels, ii, BKID) ;
      fprintf (stderr, "++++ ii=%d pos=%d daa/dxx==%d/%d c1/c2=%d/%d\n",ii,b->pos, b->daa, b->dxx, b->c1,b->c2) ;

    }
  return ;
}


/*************************************************************************************/
 
static Array  umiSeg2Hits (const UMI *umi, SEG *seg, AC_HANDLE h0) 
{
  int ii, jj, iMax, di ;
  HIT *up ;
  int t1 = umi->t1 ? umi->t1 - 1 : umi->offset ;
  int a1 = seg->a1 - umi->Delta + t1 ;
  int a2 = seg->a2 + umi->Delta + t1 ;
  Array aa = arrayHandleCreate (1024, LHIT, h0) ; 

  /* position by dichotomy */
  iMax = arrayMax (umi->hits) ; 
  for (ii = iMax / 2, di = iMax / 4 ; di > 1 ; di /= 2)
    {
      up = arrp (umi->hits, ii, HIT) ;
      if (up->a2 > a1)
	ii -= di ;
      else
	ii += di ;
    }
  for (up = arrp (umi->hits, ii, HIT) ; ii > 0 ; ii--, up--)
    {
      if (up->a2 < a1)
	break ;
     }

  /* collate the relevant hits */
  for (ii = jj = 0, up = arrp (umi->hits, ii, HIT) ; ii < iMax ; ii++, up++) 
    {
      if (up->a1 > a2)
	break ;
      if (up->a1 < a2 && up->a2 > a1)
	{
	  LHIT *vp = arrayp (aa, jj++, LHIT) ;
	  vp->lane = up->lane ;
	  vp->clone = up->clone ;
	  vp->mult = up->mult ;
	  vp->dna = up->dna ;
	  vp->uu = up->uu ;
	  vp->a1 = up->a1  - t1  ;  /* local coords */
	  vp->a2 = up->a2  - t1  ;  /* local coords */
	  vp->x1 = up->x1 ;
	  vp->x2 = up->x2 ;
	  vp->c1 = up->c1 ;
	  vp->c2 = up->c2 ;
	  vp->z1 = up->z1 ;
	  vp->z2 = up->z2 ;
	  if (vp->x1 < 0)
	    { vp->x1 = - up->x1 ; vp->isDown = TRUE ; }
	  if (vp->x2 < 0)
	    { vp->x2 = - up->x2 ; vp->isRead1 = TRUE ; }
	}
    }

  return aa ;
} /* umiSeg2Hits */

/*************************************************************************************/
/* Position segs in the iDown direction */
static void  umiCombSegPosition (const UMI *umi, SEG *seg, Array hits, BOOL isDown, AC_HANDLE h0)
{
  int kk = 0, ii, iMax = arrayMax (hits) ;
  LHIT *up ;
  int b1 = seg->b1  ;
  int b2 = seg->b2  ;
  int c1 = seg->c1  ;
  int c2 = seg->c2  ;
  int Ln = arrayMax (umi->dna) ;
  Array dnaLong = 0 ;
  int delta = umi->Delta < 9 ? umi->Delta - 1 : 8 ;
  if (! isDown)
    {
      int dummy ;
      b1 = Ln - b1 + 1 ;
      b2 = Ln - b2 + 1 ;
      dummy = b1 ; b1 = b2 ; b2 = dummy ;  
       
      c1 = Ln - c1 + 1 ;
      c2 = Ln - c2 + 1 ;
      dummy = c1 ; c1 = c2 ; c2 = dummy ;  
       
      dnaLong = umi->dnaR ;
    }
  else
    dnaLong = umi->dna ;

  for (ii = 0 ; ii < iMax ; ii++)
    {
      int jMax ;
      
      up = arrp (hits, ii, LHIT) ;
      if (isDown != up->isDown)
	continue ;
      
      /* detect the position on the read corresponding to b1 */
      /* 
	 if (! up->err || up->a1 > b1 || up->a2 < b1)
	 continue ; 
      */
      up->b1 = b1 ;
      jMax = arrayMax (up->err) ;
      if (! jMax)
	up->y1 = up->x1 + b1 - up->a1 ; 
      else
	{
	  int j ;
	  A_ERR *ep ;
	  for (j = 0, ep = arrp (up->err, 0, A_ERR) ; j < jMax && ep->iLong < b1 ; j++, ep++)
	    { ; }
	  if (j >= jMax)
	    up->y1 = up->x2 + b1 - up->a2 ; 
	  else
	    up->y1 = ep->iShort + b1 - ep->iLong ;
	}

      up->b2 = b2 ;
      jMax = arrayMax (up->err) ;
      if (! jMax)
	up->y2 = up->x1 + b2 - up->a1 ; 
      else
	{
	  int j ;
	  A_ERR *ep ;
	  for (j = 0, ep = arrp (up->err, 0, A_ERR) ; j < jMax && ep->iLong < b2 ; j++, ep++)
	    { ; }
	  if (j >= jMax)
	    up->y2 = up->x2 + b2 - up->a2 ; 
	  else
	    up->y2 = ep->iShort + b2 - ep->iLong ;
	}

      if (up->b1 > up->a2 - delta|| up->b1 < up->a1)
	up->b1 = up->y1 = 0 ;
      if (up->b2 < up->a1 + delta || up->b2 > up->a2)
	up->b2 = up->y2 = 0 ;
      if (umi->debug)
	{
	  char *strand = isDown ? "DOWN" : "UP" ;
	  if (! kk++) { fprintf (stderr, "TARGET_1    \t%s\t", strand) ; showDna (dnaLong, up->b1) ; }
	  fprintf (stderr, "READ_1 %03d >\t%s%s\t", up->mult, strand, up->isB ? "B" : "") ; showDna (up->dna, up->y1) ;
	}
    }
  return ;
} /* umiCombSegPosition */

/*************************************************************************************/
/* Comb segs in the iDown direction */
static BOOL umiCombSeg (const UMI *umi, SEG *seg, Array hits, BOOL isDown, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  int iMax = arrayMax (hits) ;
  LHIT *up ;
  int pos = 0 ;
  int b1 = seg->b1 ;
  int b2 = seg->b2 ;
  int c1 = seg->c1 ;
  int c2 = seg->c2 ;
  Array dna, css ;
  BOOL ok = TRUE ;
  int *cumulp ;
  Array cumuls ;
  Array nn = arrayHandleCreate (4, BASECOUNT, h0) ;
  int b2b[256] ;
  int Ln = arrayMax (umi->dna) ;
  Array dnaLong = 0 ;
  
  if (! isDown)
    {
      int dummy ;
      b1 = Ln - b1 + 1 ;
      b2 = Ln - b2 + 1 ;
      dummy = b1 ; b1 = b2 ; b2 = dummy ;  
      
      c1 = Ln - c1 + 1 ;
      c2 = Ln - c2 + 1 ;
      dummy = c1 ; c1 = c2 ; c2 = dummy ;  
       
      dnaLong = umi->dnaR ;
    }
  else
    dnaLong = umi->dna ;

  memset (b2b, 0, sizeof(b2b)) ;
  b2b[A_] = 1 ;
  b2b[T_] = 2 ;
  b2b[G_] = 3 ;
  b2b[C_] = 4 ;

  array (nn, 3, BASECOUNT).n = 0 ; /* make room */
  /* array of consensus DNA */
  css = arrayHandleCreate (12, CSS, h0) ;
    if (isDown)
    { 
      seg->Fcss = css ;
    }
  else
    { 
      seg->Rcss = css ;
    }

  dna = array (css, 0, CSS).dna = arrayHandleCreate (b2 - b1 + 128, char, h0) ; 
  cumuls = arrayHandleCreate (12, int, h) ;

  for (pos = 0 ; ok ; pos++) /* threading extension */
    {
      int cc, cMax = arrayMax (css) ;
      int cMaxOld = cMax ;
      
      ok = FALSE ;
      for (cc = 0 ; cc < cMaxOld ; cc++) /* try to extend each existing consensus */
	{
	  BASECOUNT *bc ;
	  int ii, base, base2 ;
	  int nnn = 0 ;
	  dna = array (css, cc, CSS).dna ;
	  cumulp = arrayp (cumuls, cc, int) ;
	  if (*cumulp < 0) 
	    continue ;
	  if (pos && *cumulp == 0) 
	    continue ;
	  for (base = 0 ; base < 4 ; base++)
	    arr (nn, base, BASECOUNT).n = 0 ;
	  for (ii = 0 ; ii < iMax ; ii++)
	    {
	      up = arrp (hits, ii, LHIT) ;
	      if (! up->dna)
		continue ;
	      up->pos = 0 ;
	      if (pos && up->cc != cc + 1)
		continue ;
	      if (isDown != up->isDown)
		continue ;
	      if (! up->b1)
		continue ;
	      
	      if (! up->dna || up->y1 + pos >= arrayMax (up->dna))
		continue ;
	      if (1 && up->y2 &&  up->y1 + pos > up->y2)
		continue ;
	      if (1 && up->z2 &&  up->y1 + pos > up->z2)
		continue ;

	      base = arr (up->dna, up->y1 + pos - 1, char) ;
	      base2 = b2b[base & 0xff] ;
	      if (base2 > 0)
		{
		  arr (nn, base2 - 1, BASECOUNT).base = base ;
		  arr (nn, base2 - 1, BASECOUNT).n += up->mult ;
		  nnn += up->mult ;
		  up->pos = pos ;
		  up->cc = cc + 1 ;
		}
	    }
	  if (! nnn)
	    {
	      *cumulp = 0 ;
	      continue ;
	    }
	  ok = TRUE ;
	  arraySort (nn, umiBaseCountOrder) ;
	  bc = arrp (nn, 0,  BASECOUNT) ;
	  array (dna, pos, char) = bc->base ; 
	  (*cumulp) += bc->n ;
	  
	  /* bifurcate */
	  for (base2 = 1 ; base2 < 4 ; base2++)
	    {
	      BASECOUNT *bc2 = arrp (nn, base2,  BASECOUNT) ;
	      if (bc2->n)
		{  
		  Array dna2 ;
		  int *cumul2p ;
		  
		  dna2 = dnaHandleCopy (dna, h0) ;
		  
		  /* register and restore current pointers */
		  array (css, cMax, CSS).dna = dna2 ;
		  dna = array (css, cc, CSS).dna ;
		  cumul2p = arrayp (cumuls, cMax, int) ;
		  cumulp = arrayp (cumuls, cc, int) ;
		  
		  /* extend */
		  *cumul2p = *cumulp -bc->n + bc2->n ;
		  array (dna2, pos, char) = bc2->base ;
		  cMax++ ;
		  
		  /* assign the reads to the new consensus */
		  {
		    int ii2 ; 
		    int delta = umi->Delta < 9 ? umi->Delta - 1 : 8 ;

		    for (ii2 = 0 ; ii2 < iMax ; ii2++)
		      {
			LHIT *up2 ;
			
			up2 = arrp (hits, ii2, LHIT) ;
			if (! up2->pos)
			  continue ;
			if (up2->cc != cc + 1)
			  continue ;
			if (! up2->dna || up2->y1 + pos >= arrayMax (up2->dna))
			  continue ;
			if (arr (up2->dna, up2->y1 + pos - 1, char) == bc2->base)
			  up2->cc = cMax ;
		      }
		    if (pos < delta) /* only accept a single (wild type) consensus on first delta bases */
		      *cumul2p = - (*cumul2p) ;
		  }
		}
	    }
	}
    }
 
  if (1)
    {
      int ii ;
      
      for (ii = 0 ; ii < iMax ; ii++)
	{
	  up = arrp (hits, ii, LHIT) ;
	  if (isDown) up->Lcc = up->cc ;
	  else up->Rcc = up->cc ;
	  if (up->cc)
	    {
	      if (up->uu == 1)
		array (css, up->cc - 1, CSS).cumulsU += up->mult ;
	      else
		array (css, up->cc - 1, CSS).cumulsNu += up->mult ;
	    }
	  up->cc =0 ;
	}
      for (ii = 0 ; ii < arrayMax (css) ; ii++)
	{
	  CSS *up = arrp (css, ii, CSS) ;
	  if (1 && array (cumuls, ii, int) < 0)
	    {
	      up->cumulsU *= -1 ;
	      up->cumulsNu *= -1 ;
	    } 
	}
    }

  arraySort (css, umiCssOrder) ;
  if (0) 
    fprintf (stderr, "seg->a1/a2=%d/%d b1/b2=%d/%d\n"
	     ,seg->a1,seg->a2
	     ,seg->b1,seg->b2
	     ) ;
  if (0) /* 2018_02_02 */
    {
      int cc, cMax = arrayMax (css) ;
      int u1, u2 ;
      
      u1 = seg->b1 + umi->t1 - 1 ;
      u2 = seg->b2 + umi->t1 - 1 ;
      
      fprintf (stdout, "\nTARGET %s %s b1=%d b2=%d\n"
	       , umi->target
	       , isDown ? "Forward reads" : "Reverse reads"
	       , u1, u2
	       ) ;
      
      fprintf (stdout, "TARGET    \t\t\t") ; showDna (dnaLong, b1 - 1) ; fprintf (stdout, "\n") ;
      for (cc = 0 ; cc < cMax ; cc++)
	{
	  int i ;
	  char *cp, *cq ;
	  dna = array (css, cc, CSS).dna ;
	  if (arrayMax (dna) < 5)
	    continue ;
	  fprintf (stdout, "CONSENSUS %d\t%d\t%d\t"
		   , cc + 1
		   , array (css, cc, CSS).cumulsU
		   , array (css, cc, CSS).cumulsNu
		   ) ; 
	  
	  for (i = 0, cp = arrp (dna, 0, char), cq = arrp (dnaLong, b1 - 1, char) ; i < arrayMax (dna) ; i++, cp++, cq++) 
	    {
	      char cc = dnaDecodeChar [(int) *cp] ;
	      if (*cp != *cq)
		cc = ace_upper (cc) ;
	      fprintf (stdout, "%c", cc) ;
	    }
	  fprintf (stdout, "\n") ;
	}
    }

  ac_free (h) ;  
  return ok ;
} /* umiCombSeg */

/*************************************************************************************/
/* stitches left aligned segs with right aligned segs */
static void  umiMakeBridge (const UMI *umi, SEG *seg, Array hits, AC_HANDLE h0)
{
  BRIDGE *bb ;
  int ii, jj, iMax, jMax ;       
  int t1 = umi->t1 ? umi->t1 : umi->offset ;
  Array dna ;
  Array buf = arrayHandleCreate (8, char, h0) ;
  Array bufT = arrayHandleCreate (8, char, h0) ;
  /*
  KEYSET i2j = keySetHandleCreate (h0) ; 
  KEYSET j2i = keySetHandleCreate (h0) ; 
  */
  seg->bridges = arrayHandleCreate (6, BRIDGE, h0) ;

  /* reverse the Rconsensus */
  iMax = seg->Fcss  ? arrayMax (seg->Fcss)  : 0 ;
  jMax = seg->Rcss ? arrayMax (seg->Rcss)  : 0 ;
  for (ii = 0 ; ii < iMax ; ii++)
    {
      int i ;
      
      dna = arr (seg->Fcss, ii, CSS).dna ;
      /* zero terminate so we can use strstr */
      i = arrayMax (dna) ;
      array (dna, i, char) = 0 ;
      arrayMax (dna) = i ;
      if (0 && i)
	{
	  fprintf (stderr, "seg->a1=%d/%d,ii=%d\t%s\n",seg->a1,seg->a2,ii,arrp(dna,0,char)) ;
	  showDna(dna,0) ;
	}
    }
  for (jj = 0 ; jj < jMax ; jj++)
    {
      int i ;
      
      dna = arr (seg->Rcss, jj, CSS).dna ;
      /* zero terminate so we can use strstr */
      i = arrayMax (dna) ;
      array (dna, i, char) = 0 ;
      arrayMax (dna) = i ;
      reverseComplement (dna) ;
      if (0 && i)
	{
	  fprintf (stderr, "seg->a1=%d/%d,jj=%d\t%s\n",seg->a1,seg->a2,jj,arrp(dna,0,char)) ;
	  showDna(dna,0) ;
	}
    }

  /* create the wild type sequence */
  if (umi->dna && arrayMax (umi->dna) > seg->b2 - seg->b1 + 1 )
    {
      int i, ln = seg->b2 - seg->b1 + 1 ;
      
      for  (i = 0 ; i < ln ; i++)
	array (bufT, i, char) = arr (umi->dna, seg->b1 + i - 1 , char) ;
      array (bufT, i, char)  = 0 ;
      arrayMax (bufT) = i ;
      dnaDecodeArray (bufT) ;
  	 
      /* create a wild type bridge */     	 
      bb = arrayp (seg->bridges, arrayMax (seg->bridges), BRIDGE) ;
      bb->nDiff = 0 ;
      bb->type = 0 ;
      bb->ln = ln ;
      bb->dx = ln ;
      bb->ddx = 0 ;
      bb->ddy = 0 ;
      bb->targetSnippet = bufT ;
      bb->variantSnippet = bufT ;
    }

  /* for each left right pair try to match the L and R consensus */
  for (ii = 0 ; ii < iMax ; ii++)
    {
      BOOL ok = FALSE ;
      Array Ldna = arr (seg->Fcss, ii, CSS).dna ;
      int lnF =  arrayMax (Ldna) ;
      int Dx = 2 * umi->Delta ;

      if (Dx > 16) Dx = 16 ;

      if (arr (seg->Fcss, ii, CSS).cumulsU < 0 || 
	  arr (seg->Fcss, ii, CSS).cumulsNu < 0
	  )
	continue ;
      if (lnF < Dx)
	continue ;
      for (jj = 0 ; ! ok && jj < jMax ; jj++)
	{
	  int ln, dx, nDiff = 0 ;
	  char cc = 0, *cp, *cq, *cr, *cs , *ct = 0 ;
	  Array Rdna = arr (seg->Rcss, jj, CSS).dna ;
	  int lnR =  arrayMax (Rdna) ;
	  if (lnR < Dx)
	    continue ;
	  if (arr (seg->Rcss, jj, CSS).cumulsU < 0 || 
	      arr (seg->Rcss, jj, CSS).cumulsNu < 0
	      )
	    continue ;

	  /* NO: we need all possibilities to find correctly duplications 
	  if (keySet (j2i, jj)) 
	    continue ;
	  */
	  cc = 0 ; cp = arrp (Ldna, 0, char) ; 
	  if (lnF > Dx) 
	    { cc = cp[Dx] ; cp[Dx] = 0 ; }
	  cq = arrp (Rdna, 0, char) ;
	  cr = strstr (cq, cp) ;
	  if (cc)
	    cp[Dx] = cc ;
	  cc = 0 ;
                     
	  cs = 0 ;
	  if (cq && lnR > Dx && lnR > Dx)
	    {
	      cs = strstr (cp, cq + lnR - Dx) ;
	    }
	  if (cr)  /* the consensus is R left clipped */
	    {
	      /* check that cp and cr match over their overlap */
	      if (strstr (cp, cr))
		;    /* ok if cp longer than cr */
	      else if (strstr (cr, cp) )
		;   /* ok if cp longer than cr */
	      else
		continue ;

	      ct = cr ; /* the doubly clipped part */
	    }
	  else if (cs)  /* the consensus is L right-clipped */
	    {
	      ct = strstr (cp, cq) ;
	      if (ct)
		{
		  int dx = ct - cp + strlen (cq) ;
		  cr = cp + dx ; cc = *cr ; *cr = 0 ; /* will be reset to cc */
		  ct = cp ;
		}
	      else
		continue ;
	    }
	  else /* look for matches of the tip of the 2 overhangs */
	    {
	      int i, j ;
	      cs = 0 ;
	      if (cp && lnF >= Dx && lnR >= Dx)
		{
		  cr = strstr (cq, cp + lnF - Dx) ;
		  cc = cq[Dx] ; cq[Dx] = 0;
		  cs = strstr (cp, cq) ; 
		  cq[Dx] = cc ;
		}
	      if (!cs)
		continue ;
	      cr = strstr (cq, cs) ;
	      if (!cr)
		continue ;
	      /* ok we splice */
	      
	      for  (i = 0 ; i <  cs - cp ; i++)
		array (buf, i, char) = cp [i] ;
	      for ( j = 0 ; j < lnR ; i++,  j++)
		array (buf, i, char) = cq [j] ;
	      array (buf, i, char)  = 0 ;
	      arrayMax (buf) = i ;
	      dnaDecodeArray (buf) ;
	      ct = 0 ;
	    }
	  if (ct)
	    {
	      int i, ln = strlen (ct) ;

	      for  (i = 0 ; i <  ln ; i++)
		array (buf, i, char) = ct [i] ;
	      array (buf, i, char)  = 0 ;
	      arrayMax (buf) = i ;
	      dnaDecodeArray (buf) ;
	      if (cc && cr)
		*cr = cc ;
	    }

	  /* count the differences */
	  ln = arrayMax (bufT) ;
	  dx = arrayMax (buf) ;
	  for (cp = arrp (buf, 0, char), cq = arrp (bufT, 0, char), nDiff = 0 ;
	       *cp && *cq ;
	       cp++, cq++)
	    if (*cp != *cq)
	      nDiff++ ;

	  /*
	  keySet (j2i, jj) = ii + 1 ;
	  keySet (i2j, ii) = jj + 1 ;
	  */
	  bb = arrayp (seg->bridges, arrayMax (seg->bridges), BRIDGE) ;
	  bb->ii = ii ;
	  bb->jj = jj ;
	  bb->nDiff = nDiff ;
	  bb->type = nDiff ? 1 : 0 ;
	  bb->ln = ln ;
	  bb->dx = dx ;
	  bb->targetSnippet = bufT ;
	  bb->variantSnippet = dnaHandleCopy (buf, h0)  ;

	  if (1)
	    { /* counts are establised later */
	      if (nDiff == 0)
		{
		  bb->wp = array (seg->Fcss, ii, CSS).cumulsU ;
		  bb->wm = array (seg->Rcss, jj, CSS).cumulsU ;
		  bb->wpnu = array (seg->Fcss, ii, CSS).cumulsNu ;
		  bb->wmnu = array (seg->Rcss, jj, CSS).cumulsNu ;
		}
	      else
		{
		  bb->mp = array (seg->Fcss, ii, CSS).cumulsU ;
		  bb->mm = array (seg->Rcss, jj, CSS).cumulsU ;
		  bb->mpnu = array (seg->Fcss, ii, CSS).cumulsNu ;
		  bb->mmnu = array (seg->Rcss, jj, CSS).cumulsNu ;
		}
	      bb->coverp = bb->mp + bb->wp ;
	      bb->coverm = bb->mm + bb->wm ;
	      bb->cover = bb->mp + bb->mm + bb->wp + bb->wm ;

	      if (seg->a1 == 827 || seg->a2 == 827 ||  umi->debug)
		      fprintf (stdout, "SEG %s:%d:%d\tBRIDGE F %d // R %d dx=%d ln=%d ddx=%d ddy=%d\t%s%s\t%d:%d\t%s>%s\tm=%d\tw=%d\n"
		      , umi->target, seg->a1, seg->a2
		      , bb->ii + 1,  bb->jj + 1
			       , bb->dx, bb->ln, 0, 0
		      , bb->nDiff > 1 ? "MULTI_" : ""
		      , types[bb->type]
		      , t1 + seg->b1 , t1 + seg->b2 
		      , arrp (bb->targetSnippet, 0, char) 
		      , arrp (bb->variantSnippet, 0, char) 
		      , bb->mp + bb->mm
		      , bb->wp + bb->wm
		      ) ;	

	    }
	}
    }

  return ;
} /* umiMakeBridge */

/*************************************************************************************/
/* fuse identical bridges, add the wild type */
static void  umiFuseBridges (const UMI *umi, SEG *seg, Array hits)
{
  int ii, jj, iMax = seg->bridges ? arrayMax (seg->bridges) : 0 ;
  BRIDGE *bb, *bb0 ;
  int coverp, coverm ;
  int coverpnu, covermnu ;
  /*   Array bufT = 0 ; */
  /*   int ddd = 1 ; */
 
  if (! iMax) return ;

  /* remove duplicates */
  for (ii = 0 ; ii < iMax ; ii++)
    {
      char *cp ;
      bb0 = arrp (seg->bridges, ii, BRIDGE) ;
      cp = arrp (bb0->variantSnippet , 0, char) ;
      for (jj = ii + 1 ; jj < iMax ; jj++)
	{
	  char *cq ;
	  bb = arrp (seg->bridges, jj, BRIDGE) ;
	  cq = arrp (bb->variantSnippet , 0, char) ;
	  
	  if (strcmp (cp, cq))
	    {
	      if (cp == strstr (cp, cq))
		{
		  /* cq is inside but shorter than cp, prefer cp */
		}
	      else if (cp == strstr (cp, cq))
		{
		  /* cp is inside but shorter than cq, prefer cq */
		  bb0->variantSnippet = bb->variantSnippet ;
		  bb->type = -1 ;		  
		}
	      else   /* keep both cp and cq, they are different */
		continue ;
	    }

	  if (bb0->type < bb->type)
	    bb0->type = bb->type ;
	  bb->type = -1 ;
	  
	  bb0->wp += bb->wp ;
	  bb0->wm += bb->wm ;
	  bb0->wpnu += bb->wpnu ;
	  bb0->wmnu += bb->wmnu ;
	  
	  bb0->mp += bb->mp ;
	  bb0->mm += bb->mm ;
	  bb0->mpnu += bb->mpnu ;
	  bb0->mmnu += bb->mmnu ;
	  
	  bb0->coverp += bb->coverp ;
	  bb0->coverm += bb->coverm ;
	  bb0->cover += bb->cover ;
	}
    }	
  /* locate the wild type */
  for (ii = 0 ; ii < iMax ; ii++)
    {
      bb0 = arrp (seg->bridges, ii, BRIDGE) ;
      if (! bb0->type)
	break ;
    }
  /*
    bufT = bb0->targetSnippet ;
    ln = arrayMax (bufT) ;
  */

  /* add the wild counts */
    for (ii = 0 ; ii < iMax ; ii++)
    {
      bb = arrp (seg->bridges, ii, BRIDGE) ;
      if (bb->type)
	{
	  bb->wp = bb0->wp ;
	  bb->wm = bb0->wm ;
	  bb->wpnu = bb0->wpnu ;
	  bb->wmnu = bb0->wmnu ;
	}
    }
    /* compute the cover counts */
    coverp = bb0->wp ;
    coverm = bb0->wm ;
    coverpnu = bb0->wpnu ;
    covermnu = bb0->wmnu ;
    for (ii = 0 ; ii < iMax ; ii++)
      {
	bb = arrp (seg->bridges, ii, BRIDGE) ;
	if (bb->type)
	  {
	    coverp += bb->mp ;
	    coverm += bb->mm ;
	    coverpnu += bb->mpnu ;
	    covermnu += bb->mmnu ;
	  }
      }
    /* retrofit and filter */
    for (ii = 0 ; ii < iMax ; ii++)
      {
	bb = arrp (seg->bridges, ii, BRIDGE) ;
	if (bb->type > 0)
	  {
	    bb->coverp = coverp ;
	    bb->coverm = coverm ;
	    bb->cover = coverp + coverm ;
	    bb->coverpnu = coverpnu ;
	    bb->covermnu = covermnu ;
	    bb->covernu = coverpnu + covermnu ;
	  }
      }

   /* select the type */
   for (ii = 0, bb = arrp (seg->bridges, 0, BRIDGE) ;  ii < iMax ; bb++, ii++)
     {  
       if (bb->type >= 0)
	 {
	   if (bb->dx == bb->ln)
	     bb->type = bb->nDiff ? 1 : 0 ;
	   else if (bb->dx < bb->ln)
	     bb->type = 5 ;
	   else if (bb->dx > bb->ln)
	     bb->type = 3 ;
	 }
     }

} /* umiFuseBridges  */

/*************************************************************************************/
/* count reads supporting all briges */
static void  umiBridgeSupport (const UMI *umi, SEG *seg, Array hits, AC_HANDLE h0)
{
  AC_HANDLE h = 0 ;
  int ii, iMax = seg->bridges ? arrayMax (seg->bridges) : 0 ;
  int jj, jMax = hits ? arrayMax (hits) : 0 ;
  int zcoverp = 0, zcoverm = 0 ;
  int zcoverpnu = 0, zcovermnu = 0 ;
  BRIDGE *bb, *bb0 ;
  BOOL isExternal = FALSE ;
  LHIT *up ;
  int isDown ;
  DICT *dict = 0 ;
  char buf[64] ;
  if (! iMax || ! jMax) return ;

  
  if (seg->a1 >800000 && seg->a1 < 1000000)
    fprintf (stderr, "bridgesupport seg    %d %d  %d bridges\n", seg->a1, seg->a2, iMax) ;
  h = ac_new_handle () ;
  dict = dictHandleCreate (jMax, h) ;
  for (isDown = 1 ; isDown >= 0 ; isDown--)
    for (jj = 0 , up = arrp (hits, 0, LHIT) ; jj < jMax ; jj++, up++)
      {
	int ok = 0 ;
	if (
	    (up->isDown &&  !isDown) ||
	    (! up->isDown &&  isDown)
	    )
	  continue ;
	if (! up->b1 && ! up->b2)
	  continue ;
	sprintf (buf, "%d%d%d", up->lane,up->clone,up->isRead1) ; 
	if (dictFind (dict,  buf, 0))
	  continue ;
	dictAdd (dict,  buf, 0) ;
	for (ii= 0, bb = arrp (seg->bridges, ii, BRIDGE) ; ok >= 0 && ii < iMax ; ii++, bb++)
	  {
	    int ln ;
	    char *cp, *cq ;
	    if (bb->type < 0)
	      continue ;
	    cp = arrp (up->dna, up->y1 > 0 ? up->y1 - 1 : 0, char) ;
	    if (! bb->tsCoded)
	      {
		bb->tsCoded = dnaHandleCopy (bb->variantSnippet, h0) ;
		dnaEncodeArray (bb->tsCoded) ;
	      }
	    if (! isDown && ! bb->isR)
	      {
		bb->isR = TRUE ;
		reverseComplement (bb->tsCoded) ;
	      }
	    cq = arrp (bb->tsCoded, 0, char) ;
	    ln = strlen (cq) ;
	    if (! strncasecmp (cp, cq, ln) || strstr (cq, cp))
	      {
		if (ok == 0)
		  ok = ii + 1 ;
		else
		  ok = -1 ; /* i suppport 2 bridges */
	      }
	  }
	if (ok > 0)
	  {
	    bb = arrp (seg->bridges, ok - 1, BRIDGE)  ;
	    if (up->uu == 1)
	      {
		if (isDown ^ up->isB)
		  bb->zmp += up->mult ;
		else
		  bb->zmm += up->mult ;
	      }
	    else
	      {
		if (isDown ^ up->isB)
		  bb->zmpnu += up->mult ;
		else
		  bb->zmmnu += up->mult ;
	      }
	  }	
      }
  /* locate the wild type */
  for (ii= 0, bb0 = 0, bb = arrp (seg->bridges, ii, BRIDGE) ;  !bb0 &&  ii < iMax ; ii++, bb++)
    if (! bb->type) bb0 = bb ;
  /* add the wild type counts */
  for (ii= 0, bb = arrp (seg->bridges, ii, BRIDGE) ;  bb0 &&  ii < iMax ; ii++, bb++)
    {
      if (bb->type < 0)
	continue ;
      bb->zwp = bb0->zmp ;
      bb->zwm = bb0->zmm ;
      bb->zwpnu = bb0->zmpnu ;
      bb->zwmnu = bb0->zmmnu ;
    }
  /* count the coverage */
   for (ii= 0, bb = arrp (seg->bridges, ii, BRIDGE) ;  ii < iMax ; ii++, bb++)
    {
      if (bb->type < 0)
	continue ;
      zcoverp +=  bb->zmp ;
      zcoverm +=  bb->zmm ;
      zcoverpnu +=  bb->zmpnu ;
      zcovermnu +=  bb->zmmnu ;
    }
  /* add the cover */
   for (ii= 0, bb = arrp (seg->bridges, ii, BRIDGE) ; ii < iMax ; ii++, bb++)
     {
       if (bb->type < 0)
	 continue ;
       bb->zcoverp = zcoverp ;
       bb->zcoverm = zcoverm ;
       bb->zcoverpnu = zcoverpnu ;
       bb->zcovermnu = zcovermnu ;
     }

   arraySort (seg->bridges, umiBridgeOrder) ;
   for (ii= 0, bb = arrp (seg->bridges, ii, BRIDGE) ;  ii < iMax ; ii++, bb++)
     {
       int mp = bb->zmp ;
       int mm = bb->zmm ;
       int coverp = bb->zcoverp ;
       int coverm = bb->zcoverm ;
       int cover = coverp + coverm ;
       int m = mm + mp ;
       int mpnu = bb->zmpnu ;
       int mmnu = bb->zmmnu ;
       int coverpnu = bb->zcoverpnu ;
       int covermnu = bb->zcovermnu ;
       int covernu = coverpnu + covermnu ;
       int mnu = mmnu + mpnu ;
       float ff ;
       float ffnu ;
       int ddd = 1 ;

       if (ii == 0 && umi->externalBitSet &&
	   seg->a2 - seg->a1 < 20)
	 {
	   BitSet ebs = umi->externalBitSet ;
	   int x ;
	   int xMax = bitSetCount (ebs) ;
	   for (x = seg->a1 ; x <= seg->a2 && x < xMax ; x++)
	     if (bit (ebs, x))
	       { isExternal = seg->isExternal = TRUE ; break ; }
	 }
       if (bb->type <= 0)
	 continue ;

       if (isExternal)
	 { isExternal = FALSE ; bb->isExternal = TRUE ; continue ; } /* export at least once */
       if (! seg->isSubseg && seg->a2 > seg->a1 + 4 * umi->Delta)
	 {
	   ddd = 2 ;
	   if (ii < 2)
	     ddd = 10 ;
	 }
       
       if (cover + covernu < umi->minSnpCover/ddd) 
	 bb->type = -2 ;
       ff = 100.0 * (m) / cover ;
       
       ffnu = covernu ? 100.0 * (mnu) / covernu : 0 ;
       
       if (1 && ff < umi->minSnpFrequency/ddd  && ffnu  < umi->minSnpFrequency/ddd)
	  bb->type = -3 ;
       if (m + mnu < umi->minSnpCount/ddd)
	  bb->type = -4 ;
     }

   ac_free (h) ;
} /* umiBridgeSupport */

/*************************************************************************************/
/* set common letters to lower */
static BOOL  umiBridgeSetUpper (const UMI *umi, SEG *seg, Array hits, AC_HANDLE h0)
{
  int i, ii ;
  int iMax = seg->bridges ? arrayMax (seg->bridges) : 0 ;
  BRIDGE *bb, *bb0 ;
  int ln, ddx, ddy ;
  char cc, *cp ;
  Array bufT ;
  BOOL ok = FALSE ;

  if (! iMax) return ok ;
 
  for (ii = 0, bb0 = 0, bb = arrp (seg->bridges, 0, BRIDGE) ;  ii < iMax ; bb++, ii++)
    if (bb->type == 0)
      { bb0 = bb ; break ; }
  if (! bb0)
    return ok ;
  bufT = bb0->targetSnippet ;
  ln = arrayMax (bufT) ;

  for (ii = 0, bb = arrp (seg->bridges, 0, BRIDGE) ;  ii < iMax ; bb++, ii++)
    if (bb->type > 0)
      ok = TRUE ;

  /* setup all letters to upper case */
  for (ii = 0, bb = arrp (seg->bridges, 0, BRIDGE) ;  ii < iMax ; bb++, ii++)
    if (bb->type >= 0)
      for (cp = arrp (bb->variantSnippet, 0, char) ; *cp ; cp++)
	*cp = ace_upper (*cp) ;
  for (cp = arrp (bufT, 0, char) ; *cp ; cp++)
     *cp = ace_upper (*cp) ;

   /* eatup the right bases common to all bridges */
   for (cc = 1, ddy = 0 ; cc && ddy < ln - 1 ; ddy++)
     {
       cc = arr (bufT, ln - ddy - 1, char) ;
       for (ii = 0, bb = arrp (seg->bridges, 0, BRIDGE) ;  cc && ii < iMax ; bb++, ii++)
	   if (bb->type >= 0)
	     {
	       if (bb->dx < ddy + 1)
		 cc = 0 ;
	       else 
		 {
		   if (cc != arr (bb->variantSnippet, bb->dx - ddy - 1, char))
		     cc = 0 ;
		 }
	     }
       if (!cc) break ; /* do not increase ddy */
     }
   /* lower the right common letters */
   for (i = 0, cp = arrp (bufT, ln - 1, char) ; *cp && i < ddy ; cp--, i++)
     *cp = ace_lower (*cp) ;
   for (ii = 0, bb = arrp (seg->bridges, 0, BRIDGE) ;  ii < iMax ; bb++, ii++)
     {
       bb->ddy = ddy ;
         if (bb->type >= 0)
	   for (i = 0, cp = arrp (bb->variantSnippet, bb->dx - 1, char)  ; i < ddy ; i++, cp--)
	     *cp = ace_lower (*cp) ;
     }

   /* eatup the left bases common to all bridges */
   for (cc = 1, ddx = 0 ; cc && ddx < ln ; ddx++)
     {
	 cc = arr (bufT, ddx, char) ;
	 for (ii = 0, bb = arrp (seg->bridges, 0, BRIDGE) ;  cc && ii < iMax ; bb++, ii++)
	     if (bb->type >= 0)
	       if (cc != arr (bb->variantSnippet, ddx, char))
		 cc = 0 ;
	 if (!cc) break ; /* do not indrease ddx */
     }
   /* lower the left common letters */
   for (i = 0, cp = arrp (bufT, 0, char) ; *cp && i < ddx ; cp++, i++)
     *cp = ace_lower (*cp) ;
   for (ii = 0, bb = arrp (seg->bridges, 0, BRIDGE) ;  ii < iMax ; bb++, ii++)
     {
       bb->ddx = ddx ;
       if (bb->type >= 0)
	 for (cp = arrp (bb->variantSnippet, 0, char), i = 0 ; *cp && i < ddx ; cp++, i++)
	   *cp = ace_lower (*cp) ;
     }

   /* lower the central  common letters */
   {
     int L1 = strlen (arrp (bufT, 0, char) ) ;
     for (ii = 0, bb = arrp (seg->bridges, 0, BRIDGE) ;  L1 &&  ii < iMax ; bb++, ii++)
       {
	 if (bb->type < 0)
	   continue ;
	 cp = arrp (bb->variantSnippet, 0, char) ;
	 if (strlen (cp) != L1)
	   L1 = 0 ;
       }
     if (L1)
       {
	 char cc, *cq ;
	 for (i = 0, cp = arrp (bufT, 0, char) ; *cp ;  cp++, i++)
	   {
	     cc = *cp ;
	     for (ii = 0, bb = arrp (seg->bridges, 0, BRIDGE) ; cc && ii < iMax ; bb++, ii++)
	       {
		 if (bb->type < 0)
		   continue ;
		 cq = arrp (bb->variantSnippet, i, char) ;
		 if (cc != *cq)
		   cc = 0 ;
	       }
	     if (cc) /* common letter, go lower case */
	       {
		 *cp = ace_lower (*cp) ;
		 for (ii = 0, bb = arrp (seg->bridges, 0, BRIDGE) ; cc && ii < iMax ; bb++, ii++)
		   {
		     if (bb->type < 0)
		       continue ;
		     cq = arrp (bb->variantSnippet, i, char) ;
		     *cq = ace_lower (*cq) ;
		   }
	       }
	     else /* someone differs, go upper case */
	       {
		 *cp = ace_upper (*cp) ;
		 for (ii = 0, bb = arrp (seg->bridges, 0, BRIDGE) ; cc && ii < iMax ; bb++, ii++)
		   {
		     if (bb->type < 0)
		       continue ;
		     cq = arrp (bb->variantSnippet, i, char) ;
		     *cq = ace_upper (*cq) ;
		   }
	       }
	   }
       }
   }
   return ok ;
} /* umiBridgeSetUpper */

/*************************************************************************************/
/* duplicate the reads so that they come x2->x1 so that we can extend them in prime */
static void  umiAddReverseSeg (const UMI *umi, SEG *seg, Array hits,AC_HANDLE h0)
{
  int dummy, jj, ii, iMax = arrayMax (hits) ;
  LHIT *up, *vp ;
  Array dnaR ;
  int ln, Ln = arrayMax (umi->dna) ;

  for (ii = 0, jj = iMax ; ii < iMax ; ii++)
    {
      up = arrp (hits, ii, LHIT) ;
      if (! up->dna)
	continue ;
      vp = arrayp (hits, jj++, LHIT) ;
      up = arrp (hits, ii, LHIT) ; /* relocate */

      *vp = *up ;
      vp->isDown = ! up->isDown ;
      dnaR = dnaHandleCopy (up->dna, h0) ;
      reverseComplement (dnaR) ;
      vp->dna = dnaR ;
      vp->isB = TRUE ;
      ln = arrayMax (up->dna) ;
      if (vp->x1)     
	vp->x1 = ln - vp->x1 + 1 ;
      if (vp->x2)     
	vp->x2 = ln - vp->x2 + 1 ;
      if (vp->y1)     
	vp->y1 = ln - vp->y1 + 1 ;
      if (vp->y2)     
	vp->y2 = ln - vp->y2 + 1 ;
      if (vp->z1)     
	vp->z1 = ln - vp->z1 + 1 ;
      if (vp->z2)     
	vp->z2 = ln - vp->z2 + 1 ;

      /* complement the genome coordinates */
      if (vp->a1)     
	vp->a1 = Ln - vp->a1 + 1 ;
      if (vp->a2)     
	vp->a2 = Ln - vp->a2 + 1 ;
      if (vp->b1)
	vp->b1 = Ln - vp->b1 + 1 ;
      if (vp->b2)
	vp->b2 = Ln - vp->b2 + 1 ;
      if (vp->c1)
	vp->c1 = Ln - vp->c1 + 1 ;
      if (vp->c2)
	vp->c2 = Ln - vp->c2 + 1 ;

      /* swap the coordinates */
      dummy = vp->x1 ; vp->x1 = vp->x2 ; vp->x2 = dummy ;
      dummy = vp->y1 ; vp->y1 = vp->y2 ; vp->y2 = dummy ;
      dummy = vp->z1 ; vp->z1 = vp->z2 ; vp->z2 = dummy ;

      dummy = vp->a1 ; vp->a1 = vp->a2 ; vp->a2 = dummy ;
      dummy = vp->b1 ; vp->b1 = vp->b2 ; vp->b2 = dummy ;
      dummy = vp->c1 ; vp->c1 = vp->c2 ; vp->c2 = dummy ;
    }

   arraySort (hits, a1LOrder) ; 
       
   return ;
}  /* umiAddReverseSeg */

/*************************************************************************************/

static void  umiComplementSeg (const UMI *umi, SEG *seg, Array hits,AC_HANDLE h0)
{
  int Ln, ln, dummy, ii, iMax = arrayMax (hits) ;
  LHIT *up ;
  CLONE *clone ;
  LANE *lane ;
  Array dnaR ;

   for (ii = 0 ; ii < iMax ; ii++)
    {
      up = arrp (hits, ii, LHIT) ;

      if (! umi->SAM)
	{
	  lane = arrp (umi->lanes, up->lane, LANE) ;
	  if (! lane || ! lane->clones)
	    continue ;
	  if (up->clone > arrayMax (lane->clones))
	    continue ;
	  clone = arrp (lane->clones, up->clone, CLONE) ;
	  if (! clone)
	    continue ;
	  up->dna = up->isRead1 ? clone->dna1 : clone->dna2 ;
	}
      if (! up->dna)
	continue ;
      if (up->isDown)
	continue ;
      if ( umi->SAM)
	{ 
	  dnaR = dnaHandleCopy (up->dna, h0) ;
	  reverseComplement (dnaR) ;
	  up->dna = dnaR ;
	}

      /* complement the DNA */
      if (0)
	{
	  dnaR = dnaHandleCopy (up->dna, h0) ;
	  reverseComplement (dnaR) ;
	  up->dna = dnaR ;

	  ln = arrayMax (up->dna) ;
	  if (up->x1)     
	    up->x1 = ln - up->x1 + 1 ;
	  if (up->x2)     
	    up->x2 = ln - up->x2 + 1 ;
	  if (up->y1)     
	    up->y1 = ln - up->y1 + 1 ;
	  if (up->y2)     
	    up->y2 = ln - up->y2 + 1 ;
	  if (up->z1)     
	    up->z1 = ln - up->z1 + 1 ;
	  if (up->z2)     
	    up->z2 = ln - up->z2 + 1 ;
	}

      /* complement the genome coordinates */
      Ln = arrayMax (umi->dna) ;
      if (up->a1)     
	up->a1 = Ln - up->a1 + 1 ;
      if (up->a2)     
	up->a2 = Ln - up->a2 + 1 ;
      if (up->b1)
	up->b1 = Ln - up->b1 + 1 ;
      if (up->b2)
	up->b2 = Ln - up->b2 + 1 ;
      if (up->c1)
	up->c1 = Ln - up->c1 + 1 ;
      if (up->c2)
	up->c2 = Ln - up->c2 + 1 ;

      /* swap the coordinates */
      dummy = up->x1 ; up->x1 = up->x2 ; up->x2 = dummy ;
      dummy = up->y1 ; up->y1 = up->y2 ; up->y2 = dummy ;
      dummy = up->z1 ; up->z1 = up->z2 ; up->z2 = dummy ;
      
      dummy = up->a1 ; up->a1 = up->a2 ; up->a2 = dummy ;
      dummy = up->b1 ; up->b1 = up->b2 ; up->b2 = dummy ;
      dummy = up->c1 ; up->c1 = up->c2 ; up->c2 = dummy ;
     }

   return ;
}  /* umiComplementSeg */

/*************************************************************************************/

static void  umiAlignSeg (const UMI *umi, SEG *seg, Array hits, BOOL isDown, AC_HANDLE h0)
{
  int ii, iMax = arrayMax (hits) ;
  LHIT *up ;
  Array dna = isDown ? umi->dna : umi->dnaR ;

   for (ii = 0 ; ii < iMax ; ii++)
    {
      up = arrp (hits, ii, LHIT) ;
      if (isDown != up->isDown)
	continue ;
      if (! up->dna)
	continue ;

      up->err = arrayHandleCreate (12, A_ERR, h0) ;
      aceDnaTrackErrors (up->dna, up->x1, &(up->x2)
			 , dna, up->a1, &(up->a2)
			 , 0, up->err, 
			 3, 8, FALSE, 0
			 , TRUE) ; /* err was just created and is already set t zero */
      ((UMI*)umi)->nErrTrackings++ ; /* approximative, not thread safe */
    }
  return ;
} /* umiAlignSeg */

/*************************************************************************************/
 
static void umiSetSubsegs (const UMI *umi, SEG *seg)
{
  AC_HANDLE h = ac_new_handle () ;
  Array subsegs = 0 ;
  char *cp ;
  int b, bMax = arrayMax (seg->bridges) ;
  int dx = 0 ;
  
  dx = seg->b2 ; /* infinity */
  for (b = 0 ; b < bMax ; b++)
    {
      int i ;
      BRIDGE *bb = arrp (seg->bridges, b, BRIDGE) ;
      if (bb->type < 0)
	continue ;
      cp = arrp (bb->variantSnippet, 0, char) ;
      for (i = 0 ; *cp && *cp == ace_lower(*cp) ; cp++, i++)
	;
      if (i < dx) dx = i ;
    }
  if (dx > umi->Delta)
    {
      if (! subsegs)
	subsegs = seg->subsegs = arrayHandleCreate (4, SEG, 0) ;
      SEG *subseg = arrayp (subsegs, 0, SEG) ;
      subseg->a1 = seg->b1 + dx ;
      subseg->a2 = seg->a2 ;
      subseg->b1 = subseg->a1 - umi->Delta ;
      subseg->b2 = seg->b2 ;
      subseg->c1 = seg->c1 ;
      subseg->c2 = seg->c2 ;
      subseg->isSubseg = TRUE ;
    }

  if (subsegs && arrayMax (subsegs))
    {
      int j, jMax = arrayMax (subsegs) ;
      SEG *subseg = arrp (subsegs, 0, SEG) ;
      
      for (j = 0 ; j < jMax ; subseg++, j++)
	umiAnalyzeSeg (umi, subseg) ;
    }

  ac_free (h) ;
  return ;
} /*  umiSetSubsegs */

/*************************************************************************************/
 
static void umiAnalyzeSeg (const UMI *umi, SEG *seg)
{
  AC_HANDLE h = ac_new_handle () ;

  Array hits = umiSeg2Hits (umi, seg, h) ;
  int x = 0x3 ; /* both strands, 0x1 forward, 0x2 reverse*/

  if (arrayMax (hits))
    {
      if (x & 0x2) umiComplementSeg (umi, seg, hits, h) ;
      if (x & 0x2) umiAddReverseSeg (umi, seg, hits, h) ;
      if (x & 0x1) umiAlignSeg (umi, seg, hits, TRUE, h) ; /* fills seg->err isDown */
      if (x & 0x2) umiAlignSeg (umi, seg, hits, FALSE, h) ; /* fills seg->err NOT isDown*/
      if (x & 0x1) umiCombSegPosition (umi, seg, hits, TRUE, h) ; /* extends left-right */
      if (x & 0x2) umiCombSegPosition (umi, seg, hits, FALSE, h) ; /* extends right-left */
      if (x & 0x1) umiCombSeg (umi, seg, hits, TRUE, h) ; /* extends left to right */
      if (x & 0x2) umiCombSeg (umi, seg, hits, FALSE, h) ; /* extends right to left */
   
      umiMakeBridge (umi, seg, hits, umi->h) ; /* stich */
      umiFuseBridges (umi, seg, hits) ; /* stich */
      umiBridgeSupport (umi, seg, hits, h) ;
      if (umiBridgeSetUpper (umi, seg, hits, h) && 
	  ! seg->isSubseg &&
	  ! seg->isTooShort && 
	  seg->a2 > seg->a1 + 4 * umi->Delta)
	umiSetSubsegs (umi, seg) ;
    }
  ac_free (h) ;
  return ;
} /* umiAnalyzeSeg */

/*************************************************************************************/

static void umiAnalyzer (const void *vp)
{
  const UMI *umi = (const UMI *)vp ;
  int ii ;

  while (channelGet (umi->analyzeChan, &ii, int))
    {
      SEG *seg = arrayp (umi->segs, ii, SEG) ;
      umiAnalyzeSeg (umi, seg) ;
    }
  ii = 100 ; /* success */
  channelPut (umi->doneChan, &ii, int) ;
  return ;
} /* umiAnalyzer */

/*************************************************************************************/

static void umiAnalyze (UMI *umi)
{
  int n = 0, nn = 4 ; /* number of analysers */
  int ii ;
  int iMax = arrayMax (umi->segs) ;

  umi->analyzeChan = channelCreate (30, int, umi->h) ;
  for (ii = 0 ; ii < nn ; ii++)
    wego_go (umiAnalyzer, umi, UMI) ;  /* will channelPut (umi->doneChan) */ 

  for (ii = 0 ; ii < iMax ; ii++)
    channelPut (umi->analyzeChan, &ii, int) ;
  channelClose (umi->analyzeChan) ;

  while (nn--)
    if (! channelGet (umi->doneChan, &n, int) || n != 100)
    {
      fprintf (stderr, "umiAnalyze received a bad value %d\n", n) ;
      exit (1) ;
    }
  return ;
} /* umiAnalyze */

/*************************************************************************************/

static void umiMergeSubSegs (UMI *umi)
{
  Array segs = umi->segs ;
  int ii, jj, iMax = arrayMax (segs) ;

  /* merge subsegs */
  for (ii = 0 ; ii < iMax ; ii++)
    {
      SEG *seg = arrayp (segs, ii, SEG) ;
      Array subsegs = seg->subsegs ; /* seg may be relocated */
      
      if (subsegs)
	{
	  int j, jMax = arrayMax (subsegs) ;
	  SEG *seg1, *subseg = arrp (subsegs, 0, SEG) ;
	  
	  for (j = 0 ; j < jMax ; subseg++, j++)
	    {
	      seg1 = arrayp (segs, iMax++, SEG) ;
	      *seg1 = *subseg ;
	    }
	}
    }
  /* order bridges by support */
  for (ii = jj = 0 ; ii < iMax ; ii++)
    {
      SEG *seg = arrayp (segs, ii, SEG) ;
      Array bridges = seg->bridges ;
      BOOL ok = FALSE ;
      int b, bMax = bridges ? arrayMax (bridges) : 0 ;

      for (b = 0 ; b < bMax ; b++)
	{
	  BRIDGE *bb = arrp (bridges, b, BRIDGE) ;
	  if (bb->type >= 0)
	    ok = TRUE ;
	}
      
      if (ok)
	{
	  arraySort (bridges, umiBridgeOrder) ;
	  if (jj < ii)
	    {
	      SEG *segj = arrayp (segs, jj, SEG) ;
	      *segj = *seg ;
	    }
	  jj++ ;
	}
    }
  iMax = arrayMax (segs) = jj ;
  arraySort (segs, umiSegOrder) ;
} /* umiMergeSubSegs */

/*************************************************************************************/

static int umiRationalizeOne (UMI *umi, SEG *seg, BRIDGE *bb, int *previousp)
{
  int ddx = 0, ds, dt, ok ;
  char *cp, *cq, *cs, *ct, *cs0, *ct0 ;
  char *buf1 = bb->buf1 ;
  char *buf2 = bb->buf2 ;

  cs0 = arrp (bb->targetSnippet, 0, char) ;
  ct0 = arrp (bb->variantSnippet, 0, char) ;

  cs = cs0 ; ct = ct0 ;
  ds = strlen (cs) ;
  dt = strlen (ct) ;
  cp = cs + ds - 1 ;
  cq = ct + dt - 1 ;

  /* set to lower common bases */
  /* we may have gone too far. for example if beginning of intron is similar to next exon */
  /* we should also look for gt_at */
  
  while (cp >= cs && cq > ct && ace_upper (*cp) == ace_upper(*cq))
    { cp-- ; cq-- ; }
  while (*cs && cs < cp && ct < cq && ace_upper (*cs) == ace_upper(*ct))
    { cs++ ; ct++ ; ddx++ ;}

  ds = cp - cs + 1 ;
  dt = cq - ct + 1 ;



  ok = 0 ;
  if (ds >= umi->minIntron  && ds <= umi->maxIntron  && dt == 1) /* possible intron */
    {
      char *cs1, *cp1, *cs11, *cp11 ;   /* look downstream for gt_ag OR gc_ag OR ct_ac OR ct_ac OR ct_gc OR gt_ag */
      char *cp1Max = cp + strlen (cp) - 3 ;

      if (1 && seg->c1)
	{
	  /* do not slide intron outside of matching zone */
	  int ddx = (cs - cs0) - (seg->c1 - seg->b1) ;
	  if (ddx < 0) ddx = 0 ;
	  cs -= ddx ; cp -= ddx ;
	  ddx = seg->c2 - seg->b1 - (cp - cs0) ;
	  if (ddx < 0) ddx = 0 ;
	  cs += ddx ; cp += ddx ;
	}
      cs1 = cs ;  cp1 = cp ;         /* gt_ag foward */
      while (! ok && cp1 < cp1Max)
	if (ace_lower(cs1[0]) == 'g' && ace_lower(cs1[1]) == 't' && ace_lower(cp1[-2]) == 'a' && ace_lower(cp1[-1]) == 'g' )
	  { ok = 7 ; /* intron gt_ag */ ; cs11 = cs1 ; cp11 = cp1 ; }
	else
	  { cs1++ ; cp1++ ; }

      cs1 = cs ;  cp1 = cp ;         /* gt_ag reverse */
      while (ok < 10 && cp1 < cp1Max)
	if (ace_lower(cs1[0]) == 'c' && ace_lower(cs1[1]) == 't' && ace_lower(cp1[-2]) == 'a' && ace_lower(cp1[-1]) == 'c' )
	  { ok += 10 ;  /* intron ct_ac */ ; }
	else
	  { cs1++ ; cp1++ ; }
      
      if (ok == 17)
	{  /* 2 possible solutions: toss a coin */
	  /*
	    if (seg->b1 + cs1 - cs0  > seg->c1)
	    ok = 7 ;
	  else if (seg->b1 + cs11 - cs0  > seg->c1)
	  ok = 10 ;
	  
	  if (randint () % 2 == 0) 
	  ok = 7 ; 
	  else
	    ok = 10 ;
	  */
	  ok = *previousp ;
	}
      if (ok == 7) /* restore the correct coordinates */
	{ cs1 = cs11 ; cp1 = cp11 ; }
      if (ok)
	{ *previousp = ok ; goto foundIntron ; }

      cs1 = cs ;  cp1 = cp ;         /* gc_ag foward */
      while (! ok && cp1 < cp1Max)
	if (ace_lower(cs1[0]) == 'g' && ace_lower(cs1[1]) == 'c' && ace_lower(cp1[-2]) == 'a' && ace_lower(cp1[-1]) == 'g' )
	   { ok = 8 ; goto foundIntron ; }
	else
	  { cs1++ ; cp1++ ; }

      cs1 = cs ;  cp1 = cp ;         /* gc_ag reverse */
      while (! ok && cp1 < cp1Max)
	if (ace_lower(cs1[0]) == 'c' && ace_lower(cs1[1]) == 't' && ace_lower(cp1[-2]) == 'g' && ace_lower(cp1[-1]) == 'c' )
	   { ok = 11 ; goto foundIntron ; }
	else
	  { cs1++ ; cp1++ ; }

 
     cs1 = cs ;  cp1 = cp ;         /* at_ac foward */
     while (! ok && cp1 < cp1Max)
	if (ace_lower(cs1[0]) == 'a' && ace_lower(cs1[1]) == 't' && ace_lower(cp1[-2]) == 'a' && ace_lower(cp1[-1]) == 'c' )
	   { ok = 9 ; goto foundIntron ; }
	else
	  { cs1++ ; cp1++ ; }

      cs1 = cs ;  cp1 = cp ;        /* at_ac reverse */
      while (! ok && cp1 < cp1Max)
	if (ace_lower(cs1[0]) == 'g' && ace_lower(cs1[1]) == 't' && ace_lower(cp1[-2]) == 'a' && ace_lower(cp1[-1]) == 't' )
	  { ok = 12 ; goto foundIntron ; }
	else
	  { cs1++ ; cp1++ ; }

    foundIntron:
      if (ok)
	{ 
	  int dx = cs1 - cs ;
	  cs += dx ; cp += dx ; ct += dx ; cq += dx ; bb->type = ok ;
	  ddx += dx ; 
	  if (ok < 10)
	    {
	      bb->a1 = cs1 - cs0 ;     /* first base of intron */
	      bb->a2 = cp1 - cs0 - 1 ; /* last base of intron */
	    }
	  else
	    {
	      bb->a2 = cs1 - cs0 ;     /* first base of intron */
	      bb->a1 = cp1 - cs0 - 1 ; /* last base of intron */
	    }
	}
      cs -= ddx ; cp -= ddx ; /* restore */
    }

  if (! ok && 
      ddx > 0 &&
      (
       (ds > 0 && dt == 0) ||
       (dt > 0 && ds == 0)
       )
      )      
    { ds++ ; dt++ ; ddx-- ; cs-- ; ct-- ; }

  buf1[0] = buf2[0] = 0 ;
  if (ds < 200 && dt < 200)
    {
      if (ds > 0 && ds < 200)
	{ strncpy (buf1, cs, ds) ; buf1[ds] = 0 ; }
      if (dt > 0 && dt < 200)
	{ strncpy (buf2, ct, dt) ; buf2[dt] = 0 ; }
      
      cp = buf1 - 1 ; while (*++cp) *cp = ace_upper (*cp) ;
      cq = buf2 - 1 ; while (*++cq) *cq = ace_upper (*cq) ;
      
      cp = buf1 - 1 ; cq = buf2 - 1 ;
      while (*++cp && *++cq) 
	if (ace_upper (*cp) == ace_upper(*cq))
	  { *cp = ace_lower (*cp) ; *cq = ace_lower (*cq) ; }
      cp = buf1 + strlen (buf1) ;
      cq = buf2 + strlen (buf2) ;
      while (--cp >= buf1 && --cq >= buf2)
	if (ace_upper (*cp) == ace_upper(*cq))
	  { *cp = *cq = 0 ; ds-- ; dt-- ; }
	else
	  break ;
      cp = buf1 + strlen (buf1) ;
      cq = buf2 + strlen (buf2) ;
      while (--cp >= buf1 && --cq >= buf2)
	if (ace_upper (*cp) == ace_upper(*cq))
	  { *cp = ace_lower (*cp) ; *cq = ace_lower (*cq) ; }
     }
  else
    {
      cp = cs0 + strlen (cs) ;
      cq = ct0 + strlen (ct) ;
      while (--cp >= cs && --cq >= ct)
	if (ace_upper (*cp) == ace_upper(*cq))
	  { cp-- ; cq-- ; ds-- ; dt-- ; }
	else
	  break ;
      cp = cs0 ; 
      cq = ct0 ;
      while (*cs && *cp && ace_upper (*cp) == ace_upper(*cq))
	{ cp++ ; cq++ ; ds-- ; dt-- ; }
    }
  bb->ds = ds ;
  bb->dt = dt ;
  bb->ddx = ddx ;

  return ddx ;
} /*  umiRationalizeOne */

/*************************************************************************************/

static int umiRationalize (UMI *umi)
{ 
  int ii ;
  Array segs = umi->segs ;
  int iMax = arrayMax (segs) ;
  int previous = 7 ;

  arraySort (segs, umiSegOrder) ;
  for (ii = 0 ; ii < iMax ; ii++)
    {
      int b, bMax ;
      SEG *seg = arrayp (segs, ii, SEG) ;

      if (seg->subsegs)
	continue ;

      while (1)
	{
	  BOOL isTooShort = FALSE ;
	  bMax = arrayMax (seg->bridges) ;
	  for (b = 0 ; b < bMax ; b++) 
	    {
	      BRIDGE *bb = arrp (seg->bridges,  b, BRIDGE) ;
	      
	      if (bb->type > 0)
		{
		  umiRationalizeOne (umi, seg, bb, &previous) ;
		  if (bb->ddx == 0) /* region is too short */
		    isTooShort = TRUE ;
		}
	    }
	  if (seg->isTooShort && isTooShort)
	    {
	      if (seg->b1 > 5)
		{
		  seg->b1 -= 12 ;
		  if (seg->b1 < 1)
		    seg->b1 = 1 ;
		  seg->isTooShort = TRUE ;
		  umiAnalyzeSeg (umi, seg) ;
		}
	    }
	  else
	    break ;
	}
    }

  arraySort (segs, umiSegOrder) ;
  for (ii = 0 ; ii < iMax ; ii++)
    {
      int bMax ;
      SEG *seg = arrayp (segs, ii, SEG) ;
      
      if (seg->subsegs)
	continue ;
      
      bMax = seg->bridges ? arrayMax (seg->bridges) : 0 ;
      if  (bMax > 1)
	arraySort (seg->bridges, umiBridgeOrder) ;
    }

  for (ii = 0 ; ii < iMax ; ii++)
    {
      int jj, b1, b2, bMax, b2Max ;
      SEG *seg = arrayp (segs, ii, SEG) ;
      
      if (seg->subsegs)
	continue ;
      bMax = seg->bridges ? arrayMax (seg->bridges) : 0 ;

      for (jj = ii ; jj < iMax ; jj++)
	{
	  SEG *seg2 = arrayp (segs, jj, SEG) ; 
	  if (seg2->subsegs)
	    continue ;
	  
	  b2Max = seg2->bridges ? arrayMax (seg2->bridges) : 0 ;
	  for (b1 = 0 ; b1 < bMax ; b1++) 
	    {
	      BRIDGE *bb2, *bb1 = arrp (seg->bridges,  b1, BRIDGE) ;
	      if (bb1->type < 0)
		continue ;
	      
	      for (b2 = 0; b2 < b2Max ; b2++) 
		{
		  if (ii == jj && b2 <= b1)
		    continue ;

		   bb2 = arrp (seg2->bridges,  b2, BRIDGE) ;
		   if (bb2->type < 0)
		     continue ;
		   
		   if (bb1->ln == bb2->ln &&
		       bb1->dx == bb2->dx &&
		       bb1->ddx + seg->b1 == bb2->ddx + seg2->b1 &&
		       ! strcmp (bb1->buf1, bb2->buf1) &&
		       ! strcmp (bb1->buf2, bb2->buf2)
		       )
		     {
		       if (bb1->zmp < bb2->zmp)
			 bb1->zmp = bb2->zmp ;
		       if (bb1->zmm < bb2->zmm)
			 bb1->zmm = bb2->zmm ;
		       if (bb1->zwp < bb2->zwp)
			 bb1->zwp = bb2->zwp ;
		       if (bb1->zwm < bb2->zwm)
			 bb1->zwm = bb2->zwm ;
		       if (bb1->zcoverp < bb2->zcoverp)
			 bb1->zcoverp = bb2->zcoverp ;
		       if (bb1->zcoverm < bb2->zcoverm)
			 bb1->zcoverm = bb2->zcoverm ;
		       bb2->type = -7 ;
		     }
		}
	    }
	}
    }
  return 0 ;
} /* umiRationalize */

/*************************************************************************************/

static void umiAliLn (UMI *umi)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (umi->outFileName, ".aliLn", umi->gzo, h) ;
  KEYSET a, aliLn = umi->aliLn = arrayHandleCreate (10000, KEY, umi->h) ;
  int ii, jj, iiMax = arrayMax (umi->lanes) ;
  LANE *lane ;

  for (ii = 1 ; ii < iiMax ; ii++)
    {
      lane = arrp (umi->lanes, ii, LANE) ;
      a = lane->aliLn ;
      if (a)
	for (jj = 0 ; jj < keySetMax (a) ; jj++)
	  keySet (aliLn, jj) += keySet (a, jj) ;
    }
  for (jj = 0 ; jj < keySetMax (aliLn) ; jj++)
    aceOutf (ao, "%d\t%d\n", jj, keySet (aliLn, jj)) ;

  ac_free (h) ;
} /* umiAliLn */

 /*************************************************************************************/

static void umiExportGeneFusions (UMI *umi)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT aoGF = 0 ;

  if (1 && arrayMax (umi->geneFusions) > 1)
    {
      int i, j, iMax ;
      GF *up, *vp ;
      
      arraySort (umi->geneFusions, umiGeneFusionOrder) ;
      iMax = arrayMax (umi->geneFusions) ;
      up = vp =  arrp (umi->geneFusions, 0, GF) ;
      for (i = j = 0 ; i < iMax ; up++, i++)
	{
	  if (up->g1 > 0 &&
	      up->lane == vp->lane && 
	      up->isPair == vp->isPair && 
	      
	      up->g1 == vp->g1 && up->g2 == vp->g2 && 
	      up->t1 == vp->t1 && up->t2 == vp->t2 && 
	      
	      up->x1 == vp->x1 && up->x2 == vp->x2 && 
	      up->y1 == vp->y1 && up->y2 == vp->y2 && 
	      up->a1 == vp->a1 && up->a2 == vp->a2 && 
	      up->b1 == vp->b1 && up->b2 == vp->b2 && 
	      
	      up->wo[0] == vp->wo[0] && up->wo[1] == vp->wo[1]
	      )
	    {
	      if (vp < up) 
		{ vp->tag += up->tag ; }
	      continue ;
	    }
	  if (up->g1 > 0)
	    {
	      vp++ ; j++ ;
	      if (j < i)
		*vp = *up ;
	    }
	}
      arrayMax (umi->geneFusions) = j + 1 ;
    }
  if (umi->geneFusions)
    {
      int i ;
      int iMax =  arrayMax (umi->geneFusions) ;
      GF *gf ;
      const char *lane = umi->lane ? umi->lane : umi->run ;
      aoGF = aceOutCreate (umi->outFileName, ".geneFusion.txt", umi->gzo, h) ;
      aceOutDate (aoGF, "####", "Gene fusion candidates") ; 
      aceOutf (aoGF, "#Run\tRead\tFusion\tGene_A\tGene_B\tChrom_A\tChrom_B\tDistance\tType\tx1 A\tx2 A\tmRNA_A\tfrom\tto\tx1 B\tx2 B\tmRNA_B\tfrom\tto\tDistinct_supports\tSupport\tGene_A supports\tGene_B supports\tscore A\tscore B\tAli A\tAli B\tc1 A\tc2 A\tc1 B\tc2 B\n") ;

      for (i = 0, gf = arrayp (umi->geneFusions, 0, GF) ; i < iMax ; gf++, i++)
	if (gf->g1 > -1)
	  {
	    const char *cp1 = dictName (umi->geneDict, gf->g1) ;
	    const char *cp2 = dictName (umi->geneDict, gf->g2) ;
	    GENE *ga = umi->genes ? arrayp (umi->genes, gf->g1, GENE) : 0 ;
	    GENE *gb = umi->genes ? arrayp (umi->genes, gf->g2, GENE) : 0 ;
	    
	    if (!strncmp (cp1, "X__", 3)) cp1 += 3 ;
	    if (!strncmp (cp2, "X__", 3)) cp2 += 3 ;
	    
	    aceOutf (aoGF, "%s\t%s\t%s__%s%s\t%s\t%s"
		     , lane
		     , dictName (gf->lane->cloneDict, gf->clone) 
		     , cp1, cp2, gf->wo
		     , cp1, cp2
		     ) ;
	    aceOutf (aoGF, "\t%s\t%s\t%d\t%s"
		     , ga ? dictName (umi->geneDict, ga->chrom) : ""
		     , gb ? dictName (umi->geneDict, gb->chrom) : ""
		     , gf->chromDistance
		     , gf->isPair ? "PAIR" : "READ"
		     ) ;
	    aceOutf (aoGF, "\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d"
		     , gf->x1, gf->x2, dictName (gf->lane->cloneDict, gf->t1), gf->a1, gf->a2
		     , gf->y1, gf->y2, dictName (gf->lane->cloneDict, gf->t2), gf->b1, gf->b2
		     ) ;
	    aceOutf (aoGF, "\t%d\t%d\t%d\t%d"
		     , gf->seq, gf->tag
		     , keySet (umi->geneSupport, gf->g1)
		     , keySet (umi->geneSupport, gf->g2)
		     ) ;
	    aceOutf (aoGF, "\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
		     , gf->score1, gf->score2
		     , gf->c2 - gf->c1 + 1, gf->d2 - gf->d1 + 1
		     , gf->c1, gf->c2, gf->d1, gf->d2
		     ) ;
	  }
    }
  if (umi->fusionHisto)
    {
      int i ;
      const char *lane = umi->lane ? umi->lane : umi->run ;
      aoGF = aceOutCreate (umi->outFileName, ".geneFusionHisto.txt", FALSE, h) ;
      aceOutDate (aoGF, "###", "Gene fusion overlap histo") ;
      aceOutf (aoGF, "#Overlap\tNumber of reads\n") ;

      for (i = -100 ; i <= 100 ; i++)
	aceOutf (aoGF, "%s\t%s\t%05d\t%d\n", umi->run, lane, i, keySet (umi->fusionHisto, i + 100)) ;
    }

  ac_free (h) ;
} /* umiExportGeneFusions */

 /*************************************************************************************/

static void umiExport (UMI *umi)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii ;
  int t1 = umi->t1 ? umi->t1 : umi->offset ;
  ACEOUT ao = aceOutCreate (umi->outFileName, ".bridges.ace", umi->gzo, h) ;
  const char *method =  umi->method ? umi->method : "UMI" ;
  const char *run = umi->run ? umi->run : "UMI" ;
  Array segs = umi->segs ;
  int iMax = arrayMax (segs) ;

  arraySort (segs, umiSegOrder) ;
  for (ii = 0 ; ii < iMax ; ii++)
    {
      int b, bMax ;
      SEG *seg = arrayp (segs, ii, SEG) ;

      if (seg->subsegs)
	continue ;
      bMax = seg->bridges ? arrayMax (seg->bridges) : 0 ;
      if  (bMax > 1)
	arraySort (seg->bridges, umiBridgeOrder) ;

      for (b = 0 ; b < bMax ; b++) 
	{
	  BRIDGE *bb = arrp (seg->bridges,  b, BRIDGE) ;
	  int ddx ;
	  int mp = bb->zmp ;
	  int mm = bb->zmm ;
	  int wp = bb->zwp ;
	  int wm = bb->zwm ;
	  int coverp = bb->zcoverp ;
	  int coverm = bb->zcoverm ;
	  int cover = coverp + coverm ;
	  int m = mm + mp ;
	  int w = wm + wp ;
	  int mpnu = bb->zmpnu ;
	  int mmnu = bb->zmmnu ;
	  int wpnu = bb->zwpnu ;
	  int wmnu = bb->zwmnu ;
	  int coverpnu = bb->zcoverpnu ;
	  int covermnu = bb->zcovermnu ;
	  int covernu = coverpnu + covermnu ;
	  int mnu = mmnu + mpnu ;
	  int wnu = wmnu + wpnu ;
	  float ff , fp, fm;
	  float ffnu , fpnu, fmnu;
	  int da, ds, dt ;
	  char *cp, *cq, cc1, cc2, *cs, *ct ;

	  if (bb->isExternal)
	    {
	      int b1 ;
	      for (b1 = 1 ; b + b1 < bMax ; b1++)
		if (bb[b1].type > 0 &&
		    bb[b1].coverp + bb[b1].coverm >=  umi->minSnpCover
		    )
		  { bb->isExternal = FALSE ; break ; }
	    }

	  if (! bb->isExternal &&
	      (bb->type <= 0 || cover + covernu < umi->minSnpCover) 
	      )
	    continue ;
	  ff = 100.0 * (m) / cover ;
	  fp = coverp ? 100.0 * mp / coverp  : 0 ;
	  fm = coverm ? 100.0 * mm / coverm : 0 ;

	  ffnu = covernu ? 100.0 * (mnu) / covernu : 0 ;
	  fpnu = coverpnu ? 100.0 * mpnu / coverpnu  : 0 ;
	  fmnu = covermnu ? 100.0 * mmnu / covermnu : 0 ;
	  
	  cp = arrayp (bb->targetSnippet, bb->ln - bb->ddy, char) ;
	  cc1 = *cp ; *cp = 0 ;
	  cq = arrayp (bb->variantSnippet, bb->dx - bb->ddy, char) ;
	  cc2 = *cq ; *cq = 0 ;

	  ddx = bb->ddx ;
	  cs = bb->buf1 ; ct = bb->buf2 ;
	  ds = bb->ds ; dt = bb->dt ;
	  da = dt - ds ;

	  if (ds + dt > 200)
	    {
	      int n = strlen(cs) ;
	      if (n > 30)
		{
		  cs[30] = 0 ;
		  cs = hprintf (h, "%s...%d_bases", cs, n) ;
		}
	      n = strlen(ct) ;
	      if (n > 30)
		{
		  ct[30] = 0 ;
		  ct = hprintf (h, "%s...%d_bases", ct, n) ;
		}
	    }
	    
	  umi->snpExported++ ;
	  aceOutf (ao, "Variant \"%s:%d:%s:%s\" // seg->a1/a2=%d/%d, b1/b2=%d/%d\n"
		   , umi->target, seg->b1 + t1 + ddx 
		   , cs
		   , ct
		   , seg->a1, seg->a2
		   , seg->b1, seg->b2
		   )  ;
	  aceOutf (ao, "VCF %d %s %s\n"
		   , seg->b1 + t1 + ddx
		   , cs
		   , ct
		   )  ;

	  if (bb->type >=7 && bb->type <= 12) /* intron */
	    {
	      aceOutf (ao, "%s %d %d\n"
		       , types[bb->type]
		       , seg->b1 + bb->a1 + t1
		       , seg->b1 + bb->a2 + t1 
		       ) ;
	    }
	  else if (ds == 1 && dt == 1)
	    aceOutf (ao, "Substitution\n%s2%s\n", cs, ct) ;
	  else if (da == 0)
	    aceOutf (ao, "Multi_substitution %d %s %s\n"
		     , ds
		     , ds < 61 ? cs : ""
		     , dt < 61 ? ct : ""
		     ) ;
	  else if (ds == 1 && dt == 2)
	    aceOutf (ao, "%s%c\n", bb->type == 2 ? "Dup" : "Ins", ct[1]) ;
	  else if (dt == 1 && ds == 2)
	    aceOutf (ao, "%s%c\n", bb->type == 4 ? "Dim" : "Del", cs[1]) ;
	  else if (ds == 1 && dt > 1)
	    aceOutf (ao, "Multi_insertion %d %s\n", da, dt < 61 ? ct + 1 : "") ; 
	  else if (dt == 1 && ds > 1)
	    aceOutf (ao, "Multi_deletion %d %s\n", -da, ds < 61 ? cs + 1 : "") ; 
	  else
	    aceOutf (ao, "DelIns %d %d %s %s\n"
		     , ds, dt
		     , ds + dt < 120 ? cs + 1 : ""
		     , ds + dt < 120 ? ct + 1 : ""
		     ) ;

	  *cp = cc1 ; *cq = cc2 ;
	  if (arrayMax (bb->targetSnippet) && arrayMax (bb->targetSnippet) < 120)
	    aceOutf (ao, "Reference_genomic_sequence %s\n", arrp (bb->targetSnippet, 0, char) );
	  if (arrayMax (bb->variantSnippet) && arrayMax (bb->variantSnippet) < 120)
	    aceOutf (ao, "Observed_genomic_sequence %s\n", arrp (bb->variantSnippet, 0, char)) ;

	  aceOutf (ao, "Method %s\n", method) ;
	  aceOutf (ao, "IntMap \"%s\" %d 1\n"
		   ,  umi->target, seg->b1 + t1 + bb->ddx 
		   ) ;
	  aceOutf (ao, "Found_in_genome\n") ;
	  if (coverp) 
	    aceOutf (ao, "fCounts %s %d %d %d Frequency %.1f\n", run, mp, wp, coverp, fp) ; 
	  if (coverm)
	    aceOutf (ao, "rCounts %s %d %d  %d Frequency %.1f\n", run, mm, wm, coverm, fm) ;
	  if (cover)
	    aceOutf (ao, "nsCounts %s %d %d  %d Frequency %.1f\n", run, m, w, cover, ff) ;
	  if (coverpnu) 
	    aceOutf (ao, "fCounts %sNU %d %d %d Frequency %.1f\n", run, mpnu, wpnu, coverpnu, fpnu) ; 
	  if (covermnu)
	    aceOutf (ao, "rCounts %sNU %d %d  %d Frequency %.1f\n", run, mmnu, wmnu, covermnu, fmnu) ;
	  if (covernu)
	    aceOutf (ao, "nuCounts %s %d %d  %d Frequency %.1f\n", run, mnu, wnu, covernu, ffnu) ;

	  aceOutf (ao, "\n") ;
	}
    }
  ac_free (ao) ;
  ac_free (h) ;
} /*  umiExport */

/*************************************************************************************/
/*************************************************************************************/

static void umiGetTargetFasta (UMI *umi)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *fNam = umi->targetFileName ;
  ACEIN ai = 0 ;
  int state = 0 ;
  int n ;
  int a1, a2, x, dx ;
  int t1 = umi->t1 ;
  int t2 = umi->t2 ;
  const char *target = umi->target ;
  Array dna = umi->dna ;
  Array dnaR = umi->dnaR ;

  if (! fNam)
    goto done ;
  ai = aceInCreate (fNam, 0, h) ;
  if (! ai)
    messcrash ("Cannot open target fasta file -targetFileName %s\n", fNam) ;
 
  if (umi->runTest) 
    {
      t1 -= umi->offset ;
      if (t2)
	t2 -= umi->offset ;
    }
  if (t1 < 1) 
    { umi->t1 += (1-t1) ; t1 = 1 ; }
  a1 = a2 = -1 ; x = 0 ;
  t1-- ; if (t2) t2-- ; umi->t1-- ; /* Plato */
  while (state < 2 && aceInCard (ai))
    {
      char *cp = aceInWord (ai) ;
      if (! cp || ! *cp)
	continue ;
      switch (state)
	{
	case 0:
	  if (*cp == '>')
	    {
	      if (umi->runTest || ! strcmp (target, cp + 1))
		state = 1 ; 
	    }
	  break ;
	case 1:
	  if (*cp == '>') /* start of next target, we are  not interested */
	    { state = 2 ; break ; }
	  dx = strlen (cp) ;
	  a1 = a2 + 1 ; a2 = a1 + dx - 1 ; /* coordinates of this segment */
	  if (a2 < t1)
	    break ;
	  if (t1 > a1)
	    { cp += t1 - a1 ; dx -= (t1 - a1) ;}
	  array (dna, x + dx, char) = 0 ; /* add a double terminal zero */
	  memcpy (arrp (dna, x, char), cp, dx) ;
	  if (t2 && x + dx > t2 - t1 + 1)
	    { array (dna, t2 - t1 + 1, char) = 0 ; arrayMax (dna) = t2 - t1 + 1 ; state = 2 ; }
	  x += dx ;
	  break ;
	}
    }  
  dnaEncodeArray (dna) ;

  n = arrayMax (dna) ;
  if (n)
    { /* we must do it in this way because the umi struct is copied by value
       * so dnaR is preallocated 
       */
      array (dnaR, n, char) = 0 ;
      memcpy (arrp (dnaR, 0, char), arrp (dna, 0, char), n) ;
      arrayMax (dna) = arrayMax (dnaR) = n - 1 ;
      reverseComplement (dnaR) ;
    }

  if (umi->makeTest && umi->t1)
    {
      ACEOUT ao = aceOutCreate (umi->makeTest, ".t1", FALSE, h) ;
      if (! ao)
	messcrash ("Cannot create test file %s.t1\n", umi->makeTest) ;
      aceOutf (ao, "##OFFSET %d\t%d\n", umi->t1, umi->t2) ;
      ac_free (ao) ;
    }

 if (umi->makeTest && dna && arrayMax (dna))
    {
      int i, iMax = arrayMax (dna) ;
      ACEOUT ao = aceOutCreate (umi->makeTest, ".target.fasta", FALSE, h) ;
      Array dna2 = dnaCopy (dna) ;
      dnaDecodeArray (dna2) ;

      if (! ao)
	messcrash ("Cannot create test file %s.target.fasta\n", umi->makeTest) ;
      aceOutf (ao, ">%s\n", umi->target) ;
      for (i = 0 ; 60*i + 60 < iMax ; i++) 
	{
	  char cc, *cp = arrp (dna2, 60*i, char) ;
	  cc = cp[60] ; cp[60] = 0 ;
	  aceOutf (ao, "%s\n",cp) ;
	  cp[60] = cc ;
	}
      aceOutf (ao, "%s\n", arrp (dna2, 60*i, char) ) ;
      ac_free (ao) ;
    }
  wego_log (hprintf (umi->h, "umiGetTargetFasta parsed %d bases in %s\n"
		     , arrayMax (dna), fNam)) ;
  
  if (! arrayMax (dna))
    {
      fprintf (stderr, "FATAL ERROR: no target sequence %s in file %s\n"
	       , umi->target
	       , fNam
	       ) ;
      ac_free (ai) ;
      exit (1) ;
    }
 done:
  ac_free (ai) ;
  ac_free (h) ;
  n = 0 ; /* success */
  t1 = umi->t1 ;
  channelPut (umi->doneChan, &(t1), int) ;

  return ;
} /* umiGetTargetFasta */

/*************************************************************************************/
/*************************************************************************************/

static BOOL umiGetOneFastc (const UMI *umi, LANE *lane)
{
  AC_HANDLE h = ac_new_handle () ;
  char *fNam1 = umi->fastcFile || ! umi->fastcDir ? 0 : 
    hprintf (h, "%s/%s.fastc%s"
	     , umi->fastcDir
	     , dictName (umi->laneDict, lane->lane)
	     , umi->runTest ? "" : ".gz"
	     ) ;
  const char *fNam =  fNam1 ? fNam1 : umi->fastcFile ;
  ACEIN ai = aceInCreate (fNam, 0, h) ;
  ACEOUT ao = 0 ;
  BOOL state = 0 ;
  DICT *dict = lane->cloneDict ;
  Array clones = lane->clones ;
  int n, nn = 0 ;
  CLONE *clone ;

  if (! fNam || ! ai)
    messcrash ("Cannot open lane fastc file -targetFileName %s\n", fNam ? fNam : "unspecified") ;
  if (umi->makeTest)
    ao = aceOutCreate (umi->makeTest, hprintf (h, ".%d.fastc", lane->lane), FALSE, h) ; 
  while (aceInCard (ai))
    {
      char *cq, *cp = aceInWord (ai) ;
      if (! cp || ! *cp || *cp == '#')
	continue ;
      switch (state)
	{
	case 0:
	  if (*cp == '>')
	    {
	      if (dictFind (dict, cp+1, &n))
		{
		  state = 1 ; 
		  clone = arrayp (clones, n, CLONE) ;
		  if (ao)
		    aceOutf (ao, "%s %s\n", cp,  aceInPos (ai)) ;
		}
	    }
	  break ;
	case 1:
	  if (ao)
	    aceOutf (ao, "%s %s\n", cp, aceInPos (ai)) ;
	  cq = strstr (cp, "><") ;
	  if (cq) *cq = 0 ;
	  n = strlen (cp) ;
	  if (n > 2)
	    {
	      nn++ ;
	      clone->dna1 = arrayHandleCreate (n, char, lane->h) ;
	      arrayForceFeed (clone->dna1, n)  ;
	      memcpy (arrp (clone->dna1, 0, char), cp, n) ; 
	      dnaEncodeArray (clone->dna1) ;
	    }
	  if (cq && (cp = cq + 2) && (n = strlen (cp)) && n > 2)
	    {
	      nn++ ;
	      clone->dna2 = arrayHandleCreate (n, char, lane->h) ;
	      arrayForceFeed (clone->dna2, n)  ;
	      memcpy (arrp (clone->dna2, 0, char), cp, n) ; 
	      dnaEncodeArray (clone->dna2) ;
	    }
	  state = 0 ;
	  break ;
	}
    }  

  if (! umi->runTest)
    wego_log (hprintf (umi->h, "umiGetOneFastc found %d reads in %s\n", nn, fNam)) ;
  ac_free (ai) ;
  ac_free (ao) ;
  ac_free (h) ;
  return TRUE ;
} /* umiGetOneFastc */

/*************************************************************************************/

static BOOL umiGetOneHit (const UMI *umi, LANE *lane)
{ 
  AC_HANDLE h = ac_new_handle () ;
  char *fNam1 = umi->hitFile ? 0 : hprintf (h, "%s/%s.hits%s"
					    , umi->hitDir
					    , dictName (umi->laneDict, lane->lane)
					    , umi->runTest ? "" : ".gz"
					   ) ;
  const char *fNam = umi->hitFile ? umi->hitFile : fNam1 ;
  ACEIN ai = aceInCreate (fNam, 0, h) ;
  ACEOUT ao = 0 ;
  int nn = 0 ;
  int score, bestScore ;
  char buf[256] , oldbuf[256], bigBuf[21024], skipBuf[21024] ;
  char bestTarget_class[64] ;
  const char *tc = umi->target_class ;
  const char *target = umi->target ;
  int t1 = umi->t1 ? umi->t1 : umi->offset ;
  int t2 = umi->t2 ;
  DICT *dict = lane->cloneDict ;
  Array hits = lane->hits ;
  KEYSET brks = lane->brks ;
  Array brkIndels = lane->brkIndels ;
  KEYSET covers = lane->covers ;
  Array geneFusions = 0 ;
  KEYSET fusionHisto = 0 ;
  KEYSET geneSupport = 0 ;
  KEYSET geneFusionKs1, geneFusionKs2 ; /* genes touched by read 1 and 2 of a pair */
  HIT *hit ;  
  int nHit = 0, lastGoodGF = 0 ;
  long int nBp = 0, nBrks = 0, nBrkIndels = 0, nOverhangs = 0 ;
  int minOverhang = umi->minOverhang ;

  if (! ai)
    messcrash ("Cannot open lane hits file -targetFileName %s\n", fNam) ;
  
  if (umi->makeTest)
    ao = aceOutCreate (umi->makeTest, hprintf (h, ".%d.hits", lane->lane), FALSE, h) ; 

  if (umi->geneFusion)
    {
      geneFusions = lane->geneFusions = arrayHandleCreate (1000, GF, lane->h) ;
      fusionHisto = lane->fusionHisto = keySetHandleCreate (lane->h) ;
      geneSupport = lane->geneSupport = keySetHandleCreate (lane->h) ;
    }
  memset (bestTarget_class, 0, sizeof (bestTarget_class)) ;
  oldbuf[0] = 0 ;
  bigBuf[1023] = 0 ;
  skipBuf[1023] = 0 ;
  if (t2 == 0)
    t2 = 1 << 30 ;
  while (aceInCard (ai))
    {
      BOOL isRead1 = TRUE ;
      int i, a1, a2, b1, b2, x1, x2, x11, x22, mult = 1, uu = 1, ali, toali, gene= 0, trgt , target_class ;
      char *cp = aceInPos (ai) ; 
      
      if (!cp || ! *cp || *cp == '#' || strlen (cp) > 20001)
	continue ;
      
      strncpy (bigBuf, aceInPos (ai), 20002) ;
      if (! strcmp (bigBuf, skipBuf))
	continue ;
      cp = aceInWord (ai) ;
      strncpy (buf, cp, 255) ;
      
      cp = buf + strlen (buf) - 1 ;
      if (*cp == '<') { isRead1 = FALSE ; *cp = 0 ; }
      else if (*cp == '>') { isRead1 = TRUE ; *cp = 0 ; }
      else  isRead1 = TRUE ;
      if (strcmp (buf, oldbuf))
	{
	  if (umi->geneFusion)
	    {
	      if (! geneFusionKs1)
		{
		  geneFusionKs1 = keySetHandleCreate (h) ;
		  geneFusionKs2 = keySetHandleCreate (h) ;
		}
	      if (keySetMax (geneFusionKs1) + keySetMax (geneFusionKs2) > 2 &&  arrayMax (geneFusions) > lastGoodGF)
		{ /* eliminate pair-fusions g1-g2 where read1 or read 2 touches both genes */
		  int i ;
		  GF *gf ;
		  for (i = lastGoodGF, gf = arrp (geneFusions, i, GF)  ; i < arrayMax (geneFusions) ; gf++, i++)
		    if (1 && gf->isPair)
		      {
			if (
			    (keySetFind (geneFusionKs1, gf->g1, 0) && keySetFind (geneFusionKs1, gf->g2, 0)) ||
			    (keySetFind (geneFusionKs2, gf->g1, 0) && keySetFind (geneFusionKs2, gf->g2, 0))
			    )
			  gf->g1 = gf->g2 = -1 ; /* kill */
		      }
		}
	      lastGoodGF =  arrayMax (geneFusions) ;
	      geneFusionKs1 = keySetReCreate (geneFusionKs1) ;
	      geneFusionKs2 = keySetReCreate (geneFusionKs2) ;
	    }
	  bestTarget_class[0] = 0 ; bestScore = 0 ;
	}
      
      if (0 && 
	  strcmp (buf, "seq.9184011") &&   /* HSG19t1r1 gest3 */
	  strcmp (buf, "seq.1451819") &&   /* HSG19t1r1 gest3 */
	  strcmp (buf, "seq.2300998") &&   /* PFALt1r1 gtest9 */
 	  strcmp (buf, "seq.4763229") &&   /* PFALt1r1 gtest8 */
	  strcmp (buf, "seq.696365") &&
	  strcmp (buf, "seq.1859179") && 
	  strcmp (buf, "seq.178616")
	  )
	continue ;
      
      aceInStep (ai, '\t') ; aceInInt (ai, &score) ;
      
      aceInStep (ai, '\t') ; aceInInt (ai, &mult) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &toali) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &ali) ;
      x1 = x2 = x11 = x22 = 0 ;
      aceInStep (ai, '\t') ; aceInInt (ai, &x1) ; x11 = x1 ;
      aceInStep (ai, '\t') ; aceInInt (ai, &x2) ; x22 = x2 ;
      
      for(i = 8 ; i <= 8 ; i++) /* locate target class in column 8 */
	{ aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; }
      if (! cp)
	continue ;
      if (*bestTarget_class && strcasecmp (cp, bestTarget_class))
	continue ;
      strncpy (bestTarget_class, cp, 63) ; 
      if (tc && strcasecmp (cp, tc))
	continue ;
      if (score < bestScore)
	bestScore = score ;
      strcpy (oldbuf, buf) ;

      dictAdd (dict, cp, &target_class) ;
      aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* gene */
      gene = 0 ;
      if (umi->geneDict)
	dictFind (umi->geneDict, cp, &gene) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &uu) ; /* unicity */
      aceInStep (ai, '\t') ; cp =  aceInWord (ai) ;
      if (! cp)
	continue ;
      if (target && ! umi->runTest && strcasecmp (cp, target))
	continue ;
      dictAdd (dict, cp, &trgt) ;
      a1 = a2 = 0 ;
      aceInStep (ai, '\t') ; aceInInt (ai, &a1) ; b1 = a1 ; 
      aceInStep (ai, '\t') ; aceInInt (ai, &a2) ; b2 = a2 ;
      if (tc && umi->geneDict && !strcmp (tc, "Z_genome"))
	{
	  GENE *ga ;
	  dictAdd (umi->geneDict, messprintf("%s__%d", dictName (dict, trgt), a1/100000), &gene) ;
	  ga = arrayp (umi->genes, gene, GENE) ;	
	  dictAdd (umi->geneDict, dictName (dict, trgt), &(ga->chrom)) ;
	  ga->a1 = ga->b1 = (a1/100000) ;
	  ga->a2 = ga->b2 = ga->a1 + 1 ;
	  ga->a1 = ga->b1 = 100000 * ga->a1 ;
	  ga->a2 = ga->b2 = 100000 * ga->a2 ;
	}
      if (a1 > a2) { int a0 = a1 ; a1 = a2 ; a2 = a0 ; a0 = x1 ; x1 = x2 ; x2 = a0 ; }
      if (t1 + t2 && (a1 < t1 || a2 > t2))
	continue ;
      
      if (ao)
	aceOutf (ao, "%s\n", bigBuf) ;
      /* we are in the correct zone, we register this hit */
      dictAdd (dict, buf, &nn) ;
      hit = arrayp (hits, nHit++, HIT) ;
      hit->lane = lane->lane ;
      hit->clone = nn ;
      hit->gene = gene ;
      if (gene > 0) keySet (geneSupport, gene)++ ;
      if (isRead1)
	keySetInsert (geneFusionKs1, gene) ;
      else
	keySetInsert (geneFusionKs2, gene) ;
      hit->target_class = target_class ;
      hit->target = trgt ;
      hit->mult = mult ;
      hit->uu = uu ;
      hit->a1 = a1 ; hit->a2 = a2 ; /* local coords */
      hit->x1 = x1 ; hit->x2 = x2 ;
      if (x1 < x2) hit->x1 = -x1 ; /* hack isDown flag */
      if (isRead1)  hit->x2 = -x2 ; /* hack isRead1 flag */
      
       /* register the cover wiggle */ 
       if (a1 < t1) a1 = t1 ;
       if (a2 > t2) a2 = t2 ;
       for (i = a1 ; i <= a2 ; i++)
	 {
	   int k = keySet (covers, i - t1) ;
	   k += mult ; nBp += mult ;
	   keySet (covers, i - t1) = k ;
	 }
       /* register the error wiggle */  
       for(i = 14 ; i <= 17 ; i++) /* locate errors in column 17 */
	{ aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; }
       if (umi->intron_only && cp && cp[1]) continue ;
       if (0 || umi->intron_only) cp = 0 ; /* mask the local errors */
       while (cp && *cp && *cp != '-')
	 {
	   int x ;
	   char *cq = strchr (cp, ',') ;
	   if (cq) *cq = 0 ;

	   x = -1 ;
	   if (sscanf (cp, "%d:", &x) == 1 && x >= t1 && x <=  t2)
	     {
	       int k ;
	       cp = strchr (cp, ':') ;
	       if (cp) cp++ ;
	       if (cp && *cp == '*') cp++ ;
	       if (cp)
		 {
		   int dx = 0 ;
		   if (cp && cp[0] == '+')
		     {
		       dx = 1 ;
		       while (cp[dx] == '+')
			 dx++ ;
		     }
		   else if (cp && cp[0] == '-')
		     {
		       dx = 1 ;
		       while (cp[dx] == '-')
			 dx++ ;
		       dx = -dx ;
		     }
		   if (dx)
		     {
		       BKID *b =  arrayp (brkIndels, nBrkIndels++, BKID) ; 
		       x-- ;
		       b->pos = x - t1 ;
		       b->daa = 3 ;
		       b->dxx = - dx  ;
		       b->mult = mult ;
		     }	   
		 }
	       k = keySet (brks, x - t1) ;
	       k += mult ; nBrks += mult ;
	       keySet (brks, x - t1) = k ;
	     }
	   cp = cq ? cq + 1 : 0 ;
	 } 
       /* register the significant overhangs */
       aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* 5p overhang */
       if (
	   (cp && *cp && b1 > t1 && b1 < t2 && *cp == ace_lower (*cp) && strlen (cp) > minOverhang) 
	   ||
	   (x11 > 20 && (! cp || ! cp[1] || *cp == ace_lower (*cp)))
	   ) 
	 { 
	   int k = keySet (brks, b1 - t1) ;
	   k += mult ; nOverhangs += mult ;
	   keySet (brks, b1 - t1) = k ;
	   if (1)
	     {
	       BKID *b =  arrayp (brkIndels, nBrkIndels++, BKID) ; 
	       b->pos = b1 - t1 ;
	       b->daa = 1 ;
	       b->dxx =   (x1 == x11) ? - 1 : 1 ;
	       b->mult = mult ;
	     }	   
	 }
       aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* 3p overhang */
       if (
	   (cp && *cp && b2 > t1 && b2 < t2 && *cp == ace_lower (*cp) && strlen (cp) > minOverhang)
	   ||
	   ( x22 < toali - 20 && (! cp || ! cp[1] || *cp == ace_lower (*cp)))
	   )
	{ 
	  int k = keySet (brks, b2 - t1) ;
	  k += mult ; nOverhangs += mult ;
	  keySet (brks, b2 - t1) = k ;
	  if (1)
	    {
	      BKID *b =  arrayp (brkIndels, nBrkIndels++, BKID) ; 
	      b->pos = b2 - t1 ;
	      b->daa = 1 ;
	      b->dxx =   (x1 == x11) ? 1 : -1 ;
	      b->mult = mult ;
	    }	   
	}   
       aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* 5p target prefix */
       aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* 3p target prefix */
       aceInStep (ai, '\t') ; cp =  aceInWordCut (ai, "\t", 0) ; /* col 22: */
       aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* col 23: */
       aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* col 24: */
       aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* col 25: */
       aceInStep (ai, '\t') ; cp =  aceInWordCut (ai, "\t", 0) ; /* col 26: chain */
       aceInStep (ai, '\t') ; aceInInt (ai, &(hit->c1)) ; /* col 27: chain1 */
       aceInStep (ai, '\t') ; aceInInt (ai, &(hit->c2)) ; /* col 28: chain2 */

       /* register de-uno discovered gene fusions */
       if (umi->geneFusion && nHit > 1)
	 {
	   HIT *old = hit - 1 ;
	   int dc, dcOld, dcc, cc1, cc2, sens = 0, chromDistance = 0 ;
	   int wab ;
	   char *wo ;
	   int minChromDistance =100001 ; /* 100 k , the +1 at the end is needed in Z_genome case */
	   int minAli = 20 ;

	   dc = nHit -1 ;
	   while (dc-- > 0 && 
		  old->clone == hit->clone &&
		  old->target_class != hit->target_class)
	     old-- ;
	   if (old->clone == hit->clone &&
	       old->gene > 0 && hit->gene > 0 &&  
	       hit->gene != old->gene &&
	       /* hit->x1 * old->x1 > 0 &&   same x strand */
	       hit->x2 * old->x2 > 0     /* same read of the pair */
	       )
	     { /* same read */
	       int a1, a2, b1, b2, x1, x2, y1, y2 ; 
	       HIT *up = old, *vp = hit ;
	       int overlap = 999999 ;

	       a1 = old->a1 ; a2 = old->a2 ;
	       x1 = old->x1 ; x2 = old->x2 ;
	       b1 = hit->a1 ; b2 = hit->a2 ;
	       y1 = hit->x1 ; y2 = hit->x2 ;

	       /* unhack */
	       if (x1 < 0) x1 = -x1 ; 
	       if (x2 < 0) x2 = -x2 ; 
	       if (y1 < 0) y1 = -y1 ;
	       if (y2 < 0) y2 = -y2 ;

	       cc1 = hit->c1 ; if (old->c1 > cc1) cc1 = old->c1 ;   /*  overlap start */
	       cc2 = hit->c2 ; if (old->c2 < cc2) cc2 = old->c2 ;   /*  overlap end */
	       dc = hit->c2 - hit->c1 ;      /* chain length */
	       dcOld = old->c2 - old->c1 ;      /* chain length */
	       
	       dcc = 3 * (cc2 - cc1) ;  /* 3 * length of overlap, 
					 * must be smaller than the 2 chain lengths 
					 */

	       sens = (a2 - a1) * (b2 - b1) ; 
	       wab = 0 ;
	       /* are the genes oriented identically relative to the read */
	       if (sens > 0)
		 { 
		   wo =  "++" ;
		   if ((b2 - b1) * (y2 - y1) * (y1 + y2 - x1 - x2) > 0 )
		     { overlap = y1 - x2 - 1 ; wab = 0 ; } /* oldGene ---> ---> gene */
		   else
		     { overlap = x1 - y2 - 1 ; up = hit ; vp = old ; } /* gene  ---> ---> oldGene */
		 }
	       else
		 { 
		   wab = 1 ; /* switchable */
		   if ((a2 - a1) * (x2 - x1) > 0) /* read in sens in old gene */
		     { 
		       if (x1 + x2 < y1 + y2)
			 wo = "+-" ;   /* oldGene ---> <--- gene */
		       else
			 wo = "-+" ; /* gene  ---> <--> oldGene */
		     }
		   else /* read in antisenes in old gene */
		     {
		       if (x1 + x2 < y1 + y2)
			 wo = "-+" ;    /* oldGene ---> <--- gene */
		       else
			 wo = "+-" ; /* gene  ---> <--> oldGene */
		     }
		 }
	       
	       if (wab == 1) /* genes can be switched */
		 {
		   wab = 0 ;
		   if (umi->genes)
		     {
		       GENE *ga = arrayp (umi->genes, old->gene, GENE) ;
		       GENE *gb = arrayp (umi->genes, hit->gene, GENE) ;
		       
		       if (ga && gb && (ga->chrom > gb->chrom || (ga->chrom == gb->chrom && ga->b1 > gb->b1)))
			 wab = 1 ;
		     }
		   else
		     if (old->gene > hit->gene)
		       wab = 1 ;
		 }
	       if (wab == 1) /* switch */
		 { up = hit ; vp = old ; }

	       /* check the minimal distance */
	       chromDistance = 0 ;
	       if (umi->genes)
		 {
		   GENE *ga = arrayp (umi->genes, old->gene, GENE) ;
		   GENE *gb = arrayp (umi->genes, hit->gene, GENE) ;
		   
		   if (ga && gb && ga->chrom == gb->chrom)
		     {
		       int dc2 = (gb->b1 + gb->b2)/2 - (ga->b1 + ga->b2)/2 ; /* distance of mid position a is left of b */
		       int ln = (gb->b2 - gb->b1 + ga->b2 - ga->b1)/2 ;
		       if (dc2 < 0) dc2 = - dc2 ;
		       chromDistance = dc2 - ln ;
		       if (chromDistance < minChromDistance)
			 dc = 0 ; /* eliminate this case */
		     }
		 }

	       if (dc > minAli && dcOld > minAli && dc > dcc && dcOld > dcc && old->gene > 0 && hit->gene > 0 && old->gene != hit->gene)
		 {
		   GF *gf = arrayp (geneFusions, arrayMax (geneFusions), GF) ;
		   
		   gf->lane = lane ;
		   dictAdd (dict, buf, &(gf->clone)) ;
		   gf->g1 = up->gene ;
		   gf->g2 = vp->gene ;
		   gf->t1 = up->target ;
		   gf->t2 = vp->target ;
		   gf->wo[0] = wo[0] ;
		   gf->wo[1] = wo[1] ;
		   gf->wo[2] = 0 ;
		   gf->isPair = FALSE ;
		   gf->score1 = up->score ;
		   gf->score2 = vp->score ;

		   gf->seq = 1 ;
		   gf->tag = up->mult ;
		   gf->c1 = up->c1 ;
		   gf->c2 = up->c2 ;
		   gf->d1 = vp->c1 ;
		   gf->d2 = vp->c2 ;

		   gf->x1 = ABS(up->x1) ;
		   gf->x2 = ABS(up->x2) ;
		   gf->y1 = ABS(vp->x1) ;
		   gf->y2 = ABS(vp->x2) ;

		   gf->a1 = up->a1 ;
		   gf->a2 = up->a2 ;
		   gf->b1 = vp->a1 ;
		   gf->b2 = vp->a2 ;
		   
		   gf->chromDistance = chromDistance ;
		   
		   if (overlap < 999999 && overlap > 100) overlap = 100 ;
		   if (overlap < -100) overlap = -100 ;
		   keySet (fusionHisto, overlap + 100) += 1 ;
		 }
	     }

	   /* recognize gene fusions from pairs */	   dc = nHit -1 ;
	   dc = nHit -1 ;
	   while (dc-- > 0 && 
		  old->clone == hit->clone &&
		  old->target_class != hit->target_class)
	     old-- ;
	   
	   /* eliminate read with a previous hit on hit->gene */
	   if (old->clone == hit->clone && 
	       old->gene > 0 && hit->gene > 0 && 
	       /* old->uu == 1 && hit->uu == 1 &&     the 2 reads align in a single gene */
	       hit->gene != old->gene &&
	       /* hit->x1 * old->x1 > 0 &&  same x strand */
	       hit->x2 * old->x2 < 0     /* opposite reads from same pair */
	       )
	     { /* opposite reads */
	       int a1, a2, b1, b2, x1, x2, y1, y2 ; 
	       HIT *up = old, *vp = hit ;
	       
	       a1 = old->a1 ; a2 = old->a2 ;
	       x1 = old->x1 ; x2 = old->x2 ;
	       b1 = hit->a1 ; b2 = hit->a2 ;
	       y1 = hit->x1 ; y2 = hit->x2 ;
	       
	       /* unhack */
	       if (x1 < 0) x1 = -x1 ; 
	       if (x2 < 0) x2 = -x2 ; 
	       if (y1 < 0) y1 = -y1 ;
	       if (y2 < 0) y2 = -y2 ;
	       
	       cc1 = hit->c1 ; if (old->c1 > cc1) cc1 = old->c1 ;   /*  overlap start */
	       cc2 = hit->c2 ; if (old->c2 < cc2) cc2 = old->c2 ;   /*  overlap end */
	       dc = hit->c2 - hit->c1 ;      /* chain length */
	       dcOld = old->c2 - old->c1 ;      /* chain length */
	       
	       wab = 0 ;
	       if ((x2 - x1) * (a2 - a1) > 0 &&
		   (y2 - y1) * (b2 - b1) > 0
		   ) /* both read are stranded in the gene */
		 { wo = "+-" ;  wab = 1 ; }
	       else if ((x2 - x1) * (a2 - a1) < 0 &&
			(y2 - y1) * (b2 - b1) < 0
			) /* both read are anti-stranded in the gene */
		 { wo = "-+" ; wab = 1 ; }
	       else if ((x2 - x1) * (a2 - a1) > 0 &&
			(y2 - y1) * (b2 - b1) < 0
			) /* both read are anti-stranded in the gene */
		 wo = "++" ;
	       else /* --, switch the genes */
		 { wo = "++" ; up = hit ; vp = old ; }
	       
	       if (wab == 1) /* genes ca be switched */
		 {
		   wab = 0 ;
		   if (umi->genes)
		     {
		       GENE *ga = arrayp (umi->genes, old->gene, GENE) ;
		       GENE *gb = arrayp (umi->genes, hit->gene, GENE) ;
		       
		       if (ga && gb && (ga->chrom > gb->chrom || (ga->chrom == gb->chrom && ga->b1 > gb->b1)))
			 wab = 1 ;
		     }
		   else
		     if (old->gene > hit->gene)
		       wab = 1 ;
		 }

	       /* check the minimal distance */
	        chromDistance = 0 ;
		if (umi->genes)
		 {
		   GENE *ga = arrayp (umi->genes, old->gene, GENE) ;
		   GENE *gb = arrayp (umi->genes, hit->gene, GENE) ;
		   
		   if (ga && gb && ga->chrom == gb->chrom)
		     {
		       if (ga->b1 > gb->b1)
			 chromDistance = gb->b1 - ga->b2 ;
		       else
			 chromDistance = ga->b1 - gb->b2 ;
		       if (chromDistance < minChromDistance)
			 dc = 0 ; /* eliminate this case */
		     }
		 }
	       if (wab == 1) /* switch */
		 { up = hit ; vp = old ; }
	       
	       if (dc > minAli && dcOld > minAli) /* do not check dcc overlap, not the same read */
		 {
		   GF *gf = arrayp (geneFusions, arrayMax (geneFusions), GF) ;
		   
		   gf->lane = lane ;
		   dictAdd (dict, buf, &(gf->clone)) ;
		   gf->g1 = up->gene ;
		   gf->g2 = vp->gene ;
		   gf->t1 = up->target ;
		   gf->t2 = vp->target ;
		   gf->wo[0] = wo[0] ;
		   gf->wo[1] = wo[1] ;
		   gf->wo[2] = 0 ;
		   gf->isPair = TRUE ;
		   gf->score1 = up->score ;
		   gf->score2 = vp->score ;
		   gf->seq = 1 ;
		   gf->tag = up->mult < vp->mult ? up->mult : vp->mult ;

		   gf->c1 = up->c1 ;
		   gf->c2 = up->c2 ;
		   gf->d1 = vp->c1 ;
		   gf->d2 = vp->c2 ;

		   gf->x1 = ABS(up->x1) ;
		   gf->x2 = ABS(up->x2) ;
		   gf->y1 = ABS(vp->x1) ;
		   gf->y2 = ABS(vp->x2) ;

		   gf->a1 = up->a1 ;
		   gf->a2 = up->a2 ;
		   gf->b1 = vp->a1 ;
		   gf->b2 = vp->a2 ;
		   
		   gf->chromDistance = chromDistance ;
		 }
	     }
	}
    
       /* register de-uno discovered deletions */
       if (! umi->geneFusion && nHit > 1)
	 {
	   HIT *old = hit - 1, *secondHit ;
	   int da, dx, dummy, protect = 0, c1, c2 ;

	   if (old->clone == hit->clone &&
	       hit->x1 * old->x1 > 0 &&
	       hit->x2 * old->x2 > 0
	       )
	     { /* same read */
	       int a1, a2, b1, b2, x1, x2, y1, y2 ;
	       if (old->a1 < hit->a1)
		 { 
		   secondHit = hit ;
		   a1 = old->a1 ; a2 = old->a2 ;
		   x1 = old->x1 ; x2 = old->x2 ;
		   b1 = hit->a1 ; b2 = hit->a2 ;
		   y1 = hit->x1 ; y2 = hit->x2 ;
		 }
	       else
		 { 
		   secondHit = old ;
		   b1 = old->a1 ; b2 = old->a2 ;
		   y1 = old->x1 ; y2 = old->x2 ;
		   a1 = hit->a1 ; a2 = hit->a2 ;
		   x1 = hit->x1 ; x2 = hit->x2 ;
		 }

	       c1 = a2 ; c2 = b1 ;
	       if (x1 < 0) x1 = -x1 ; 
	       if (x2 < 0) x2 = -x2 ; 
	       if (y1 < 0) y1 = -y1 ;
	       if (y2 < 0) y2 = -y2 ;
	       dx = y1 - x2 ; 
	       if (y1 < y2) dx = dx - 1 ; else  dx = -dx - 1 ; 

	       /* clean up x duplications */
	       
	       if (b1 > a1)
		 {
		   if (x1 < x2)
		     {
		       if (y1 < x2 + 1)
			 {
			   int dx = (x2 + 1) - y1 ;
			   protect = dx ;
			   y1 += dx ; b1 += dx ; dx = 0 ;
			   secondHit->a1 = b1 ;
			   if (secondHit->x1 < 0) secondHit->x1 = - y1 ; else secondHit->x1 = y1 ;
			 }
		     }
		   else
		     {
		       if (y1 > x2 - 1)
			 {
			   int dx = y1 - (x2 - 1) ;
			   protect = dx ;
			   y1 -= dx ; b1 += dx ; dx = 0 ;
			   secondHit->a1 = b1 ;
			   if (secondHit->x1 < 0) secondHit->x1 = - y1 ; else secondHit->x1 = y1 ;
			 }
		     }
		 }

	       da = b1 - a2 ;
	       if (b1 > a2)
		 {
		   da-- ;
		   if (y1 > y2 && (y1 >= x1 - umi->Delta || x2 < y2 + umi->Delta))
		     da = dx = 0 ;
		   if (y1 < y2 && (x2 > y2 - umi->Delta || y1 < x1 + umi->Delta))
		     da = dx = 0 ;
		 }
	       if (b1 < a1 + umi->Delta || a2 >= b2 -umi->Delta)
		 da = dx = 0 ;

	       if (a2 >= b1)
		 { da = a2 - b1 ; dummy = a2 ; a2 = b1 ; b1 = dummy ;  a2-- ; b1 += da + 1 ; }

	       if (dx < 0 && da > -dx)
		 { a2 += dx ; x2 += dx ; dx = 0 ; da -= dx ; } 

	       if (da || dx)
		  {
		    int k ;
		    BKID *b =  arrayp (brkIndels, nBrkIndels++, BKID) ; 
		    b->pos = a2 - t1 ;
		    b->mult = mult ;
		    b->daa = da ;
		    b->dxx = dx ;  /* da = 0 dx < 0 represents a duplication */
		    b->protect = protect ;
		    b->c1 = c1 - t1 ;
		    b->c2 = c2 - t1 ;
		    k = keySet (brks, b1 - t1) ;
		    k += mult ; nOverhangs += mult ;
		    keySet (brks, b1 - t1) = k ;

		    if (0)
		      fprintf (stderr, "gethits found c1/c2 = %d %d\n", c1, c2);
		    /* if we have a duplication, we must be out of it on both sides */
		    if (dx < 0 && -dx > umi->Delta)
		      b1 += (-dx  - umi->Delta) ;

		    if (dx < 0)
		      a2 += dx - 1 ; 
		    k = keySet (brks, a2 - t1) ;
		    k += mult ; nOverhangs += mult ;
		    keySet (brks, a2 - t1) = k ;

		    if (a2 < b1 && a2 > b1 - 60)
		      for (++a2 ; a2 < b1 ; a2 += umi->Delta/2)
			{
			  k = keySet (brks, a2 - t1) ;
			  k += mult ; nOverhangs += mult ;
			  keySet (brks, a2 - t1) = k ;
			}
		  }	   
	     }
	 }
     }

  if (! umi->runTest)
    wego_log (hprintf (umi->h, "umiGetOneHit found %d hits, covering %ld bp and %ld brks %ld overhangs in %s\n", nHit, nBp, nBrks, nOverhangs, fNam)) ;
  ac_free (ai) ;
  ac_free (ao) ;
  ac_free (h) ;
  return TRUE ;
} /* umiGetOneHit */

/*************************************************************************************/

static BOOL umiGetOneSamHit (const UMI *umi, LANE *lane)
{ 
  AC_HANDLE h = ac_new_handle () ;
  char *fNam1 =  umi->SAM_file  ? 0 : hprintf (h, "%s/%s.genome.sam%s"
					       , umi->hitDir ? umi->hitDir : ""
					       , dictName (umi->laneDict, lane->lane)
					       , umi->runTest ? "" : ".gz"
					       ) ;
  const char *fNam = umi->SAM_file ? umi->SAM_file : fNam1 ;
  ACEIN ai = aceInCreate (fNam, 0, h) ;
  ACEOUT ao = 0 ;
  int nn = 0 ;
  static int nnn = 0 ;
  int iCigar, dummy = 0 ;
  char buf[256] , bigBuf[1024], *cigar ;
  const char *target = umi->target ;
  int t1 = umi->t1 ;
  int t2 = umi->t2 ;
  DICT *dict = lane->cloneDict ;
  Array hits = lane->hits ;
  KEYSET brks = lane->brks ;
  KEYSET covers = lane->covers ;
  HIT *hit ;
  int nHit = 0 ;
  long int nBp = 0, nBrks = 0, nBrkIndels = 0, nOverhangs = 0 ; 
  Array clones = lane->clones ;
  Array cigarettes = arrayHandleCreate (128, SAMCIGAR, h) ; 
  Array brkIndels = lane->brkIndels ;
  Array myDna = arrayHandleCreate (128, char, h) ; 
  SAMCIGAR *cgr ;
  CLONE *clone ;

  if (! ai)
    messcrash ("Cannot open lane hits file -targetFileName %s\n", fNam) ;
  
  if (umi->makeTest)
    ao = aceOutCreate (umi->makeTest, hprintf (h, ".%d.hits", lane->lane), FALSE, h) ; 

  bigBuf[1023] = 0 ;
  if (t2 == 0)
    t2 = 1 << 30 ;
  while (aceInCard (ai))
    {
      BOOL isRead1 = TRUE ;
      int i, flag, a1, a2, x1, x2, mult = 1, uu = 1, score = 0, nerr = 0, ali = 0 ;
      char *seq, *cp = aceInPos (ai) ; 

      if (!cp || ! *cp || *cp == '@' || *cp == '#')
	continue ;

      cp = aceInWord (ai) ;
      strncpy (bigBuf, aceInPos (ai), 1023) ; /* fragment name */
      aceInStep (ai, '\t') ; cp = aceInWord (ai) ;
      strncpy (buf, cp, 255) ;
      cp = strchr (buf, '#') ;
      if (cp && sscanf (cp+1, "%d", &i) == 1 && i > 0)
	mult = i ;

      aceInStep (ai, '\t') ; aceInInt (ai, &flag) ; /* flag hides the read 1/2 distinction */
      isRead1 = (flag & 128 ? FALSE : TRUE) ;

      aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* target */
      if (target && ! umi->runTest && strcasecmp (cp, target))
	continue ;

      aceInStep (ai, '\t') ; aceInInt (ai, &a1) ;   /* first aligned base on top strand of the genome */
      a2 = a1 + 100 ; /* just a guess until we scan the cigar */
      if (t1 + t2 && (a1 < t1 || a2 > t2))
	 continue ;

      aceInStep (ai, '\t') ; aceInInt (ai, &dummy) ; /* a quality coefficient */

      aceInStep (ai, '\t') ; cigar =  aceInWord (ai) ; /* cigar */
      samParseCigar (cigar, cigarettes, a1, &a2, &x1, &x2, 0, 0, 0, 0) ;

      aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* not used */
      aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* not used */
      aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* not used */

      aceInStep (ai, '\t') ; seq =  aceInWord (ai) ; /* read sequence */
      dictAdd (dict, buf, &nn) ;
   
      myDna = arrayReCreate (myDna, strlen (seq) + 1, char) ;
      arrayMax (myDna) = strlen (seq) ;
      memcpy (arrp (myDna, 0, char), seq, strlen (seq)) ;
      dnaEncodeArray (myDna) ;

      aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* not used */

      aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; if (cp && !strncasecmp (cp, "NH:i:",5)) sscanf (cp + 5, "%d", &uu) ;  /* NH:i:%d target multiplicity */
      aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; if (cp && !strncasecmp (cp, "AS:i:",5)) sscanf (cp + 5, "%d", &score) ;  /* AS:i:%d alignment score */
      aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; if (cp && !strncasecmp (cp, "NM:i:",5)) sscanf (cp + 5, "%d", &nerr) ;  /* NM:i:%d edit distance to the reference */
      aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; if (cp && strlen (cp) > 5)  cp += 5 ;     /* MD:Z:%s genome cigarette */

      
      /* 4 more unclear columns */
      
      /* we are in the correct zone, we register this hit */
      if (arrayMax (cigarettes))
	for (iCigar = 0, cgr = arrp (cigarettes, 0, SAMCIGAR) ; iCigar < arrayMax (cigarettes) ; iCigar++, cgr++)
	  {  
	    switch (cgr->type)
	      {
	      case 'M':
		ali += cgr->x2 - cgr->x1 ;
		break ;
	      case 'D':
	      case 'I':
	      case 'N':/* intron */
		{
		  int k ;
		  BKID *b =  arrayp (brkIndels, nBrkIndels++, BKID) ; 
		  int dx = cgr->x2 - cgr->x1 ;
		  if (dx < 0)
		    { cgr->x2 -= dx ; cgr->a2 -= dx ; cgr->a1 += dx ; cgr->x1 += dx ; dx = 0 ; }
		  /* cgr->a1 a2 = first last base of introm, g._.g of gt_ag */
		  b->pos = cgr->a1 - 1 - t1 ;
		  b->mult = mult ;
		  b->daa = cgr->a2 - cgr->a1 + 1 ;
		  b->dxx = 0 ;  /* da = 0 dx < 0 represents a duplication */
		  k = keySet (brks, cgr->a1 -1 - t1) ;
		  k += mult ; nOverhangs += mult ;
		  keySet (brks, cgr->a1 -1 - t1) = k ;
		  
		  k = keySet (brks, cgr->a2 - t1) ;
		  k += mult ; nOverhangs += mult ;
		  keySet (brks, cgr->a2 - t1) = k ;
		  
		  k = cgr->a2 - t1 ;
		}
		continue ;
	      case 'S':
	      default:
		continue ;
	      }
	    keySet (lane->aliLn, ali) ++ ;

	    hit = arrayp (hits, nHit++, HIT) ;
	    
	    hit->lane = lane->lane ;
	    hit->clone = nn ;
	    clone = arrayp (clones, nn, CLONE) ;
	    if (arrayMax (myDna))
	      {
		if (isRead1)
		  {
		    if (! clone->dna1) clone->dna1 = dnaHandleCopy (myDna, lane->h) ;
		    hit->dna = clone->dna1 ;
		  }
		else
		  {
		    if (! clone->dna2) clone->dna2 = dnaHandleCopy (myDna, lane->h) ;
		    hit->dna = clone->dna2 ;
		  }
	      }
	    hit->mult = mult ;
	    hit->uu = uu ;
	    a1 = hit->a1 = cgr->a1 ;
	    a2 = hit->a2 = cgr->a2 ;
	    hit->x1 = cgr->x1 ;
	    hit->x2 = cgr->x2 ;

	    /* register the cover wiggle */ 
	    if (a1 < t1) a1 = t1 ;
	    if (a2 > t2) a2 = t2 ;
	    for (i = a1 ; i <= a2 ; i++)
	      {
		int k = keySet (covers, i - t1) ;
		k += mult ; nBp += mult ;
		keySet (covers, i - t1) = k ;
	      }
	    /* register the error wiggle */  

	    if (hit->dna && umi->dna &&
		cgr->x1 > 0 && cgr->a1 > 0 &&
		cgr->x1 < arrayMax (hit->dna) &&
		cgr->a1 < arrayMax (umi->dna)
		)
	      {
		char *cp = arrp (hit->dna, cgr->x1 - 1, char) ;
		char *cq = arrp (umi->dna, cgr->a1 - 1, char) ;
		int x ;
		int xMax = cgr->x2 - cgr->x1 + 1 ;
		int aMax = cgr->a2 - cgr->a1 + 1 ;
		if (aMax < xMax) xMax = aMax  ;
		for (x = 0 ; x < xMax ; x++, cp++, cq++, nnn++)
		  if (! (*cp & *cq))
		    {
		      int k = keySet (brks, cgr->a1 + x - t1) ;
		      k += mult ; nBrks += mult ;
		      keySet (brks, cgr->a1 + x - t1) = k ;
		    }
	      }
	    if (x1 < x2) hit->x1 = -hit->x1 ; /* hack isDown flag */
	    if (isRead1) hit->x2 = -hit->x2 ; /* hack isRead1 flag */
	  }
    }

  if (! umi->runTest)
    wego_log (hprintf (umi->h, "umiGetOneHit found %d hits, covering %ld bp and %ld brks %ld overhangs in %s\n", nHit, nBp, nBrks, nOverhangs, fNam)) ;
  ac_free (ai) ;
  ac_free (ao) ;
  ac_free (h) ;

  return TRUE ;
} /* umiGetOneSamHit */

/*************************************************************************************/

static void umiGetLaneHits (const void *vp)
{
  AC_HANDLE h = ac_new_handle () ;
  const UMI *umi = (const UMI *)vp ;
  int ii ;
  BOOL ok = TRUE ;
  
  while (channelGet (umi->getLaneInChan, &ii, int))
    {
      LANE *lane = arrp (umi->lanes, ii, LANE) ;

      if (ok) 
	{
	  if (umi->SAM)
	    ok = umiGetOneSamHit (umi, lane) ;
	  else
	    ok = umiGetOneHit (umi, lane) ;
	}
	  
      if (umi->targetFileName && ! umi->SAM)
	ok = umiGetOneFastc (umi, lane) ;

      channelPut (umi->getLaneOutChan, &ii, int) ;
    }

  ac_free (h) ;
  return ;
} /* umiGetLaneHits */

/*************************************************************************************/
 
static void umiFuseLaneHits (UMI *umi)
{
  AC_HANDLE h = ac_new_handle () ;

  int ii, jj, iiMax = arrayMax (umi->lanes) ;
  LANE *lane ;
  KEYSET brks, covers ;
  Array brkIndels = 0 ;
  Array hits ;
  HIT *up, *vp ;

  brks = umi->brks = arrayHandleCreate (10000, KEY, umi->h) ;
  brkIndels = umi->brkIndels = arrayHandleCreate (10000, BKID, umi->h) ;
  covers = umi->covers = arrayHandleCreate (10000, KEY, umi->h) ;
 
  /* fuse the hits sort them and recreate the dict in the same order */
  jj = 1 ; 
  for (ii = 1 ; ii < iiMax ; ii++)
    {
      lane = arrp (umi->lanes, ii, LANE) ;
      jj +=  arrayMax (lane->hits) ;
    }
  hits = umi->hits = arrayHandleCreate (jj, HIT, umi->h) ; 
  vp = arrayp (umi->hits, jj - 1, HIT) ; /* make room */
  for (jj = 0, ii = 1 ; ii < iiMax ; ii++)
    {
      int iMax ;
      lane = arrp (umi->lanes, ii, LANE) ;
      iMax = arrayMax (lane->hits) ;
      if (iMax)
	{
	  up = arrp (lane->hits, 0, HIT) ;
	  vp = arrp (hits, jj, HIT) ;
	  if (iMax) memcpy (vp, up, iMax * sizeof(HIT)) ;
	  jj += iMax ; 
	}
      arrayDestroy (lane->hits) ;
    }
  arraySort (hits, umiA1Order) ;
    
  /* fuse the coverage plots */
  for (ii = 1 ; ii < iiMax ; ii++)
    {
      KEY *kp ;
      int i, iMax ;
      lane = arrp (umi->lanes, ii, LANE) ;
      iMax = arrayMax (lane->covers) ;
      if (iMax)
	{
	  keySet (covers, iMax) += 0 ;
	  for (i = 0, kp = arrayp (lane->covers, i, KEY) ; i < iMax ; i++, kp++)
	    keySet (covers, i) += *kp ;
	}
    }
  /* fuse the brks plots */
  for (ii = 1 ; ii < iiMax ; ii++)
    {
      KEY *kp ;
      int i, iMax ;
      lane = arrp (umi->lanes, ii, LANE) ;
      iMax = arrayMax (lane->brks) ;
       if (iMax)
	{
	  keySet (brks, iMax) += 0 ;
	  for (i = 0, kp = arrayp (lane->brks, i, KEY) ; i < iMax ; i++, kp++)
	    keySet (brks, i) += *kp ;
	}
    } 

  /* fuse the brkIndels plots */
  for (ii = 1 ; ii < iiMax ; ii++)
    {
      BKID *b ;
      int i, j, iMax ;
      lane = arrp (umi->lanes, ii, LANE) ;
      iMax = arrayMax (lane->brkIndels) ;
      j =  arrayMax (brkIndels) ;
       if (iMax)
	{
	  array (brkIndels, j + iMax - 1, BKID).pos = 0 ;
	  for (i = 0, b = arrayp (lane->brkIndels, i, BKID) ; i < iMax ; i++, b++)
	    {
	      array (brkIndels, j++, BKID) = *b ;
	      if (0 && b->c1)
		fprintf (stderr, "=== fuse lane %d, iMax=%d b[%d].da=%d    brkIndels=%p max=%d c1/c2 = %d/%d\n", ii, iMax, i, b->daa, brkIndels, arrayMax (brkIndels), b->c1, b->c2) ;
	    }
	}
    }
  
  /* fuse the geneFusions */
  if (umi->geneFusion)
    {
      KEYSET fusionHisto = umi->fusionHisto = keySetHandleCreate (umi->h) ;
      Array geneFusions = umi->geneFusions = arrayHandleCreate (10000, GF, umi->h) ;
      
      for (ii = 1 ; ii < iiMax ; ii++)
	{
	  GF *b ;
	  int i, j, iMax ;
	  lane = arrp (umi->lanes, ii, LANE) ;
	  iMax = arrayMax (lane->geneFusions) ;
	  j =  arrayMax (geneFusions) ;
	  if (iMax)
	    {
	      array (geneFusions, j + iMax - 1, GF).g1 = 0 ;
	      for (i = 0, b = arrayp (lane->geneFusions, i, GF) ; i < iMax ; i++, b++)
		array (geneFusions, j++, GF) = *b ;
	    }
	  j = keySetMax (lane->fusionHisto) ;
	  for (i = 0 ; i < j ; i++)
	    keySet (fusionHisto, i) += keySet (lane->fusionHisto, i) ;
	}
    }

  /* fuse the geneSupport */
  if (umi->geneFusion)
    {
      KEYSET geneSupport = umi->geneSupport = keySetHandleCreate (umi->h) ;
      
      for (ii = 1 ; ii < iiMax ; ii++)
	{	  
	  int i, iMax ;
	  iMax = arrayMax (lane->geneSupport) ;
	  for (i = iMax - 1 ; i>- 0 ; i--)
	    keySet (geneSupport, i) += keySet (lane->geneSupport, i) ;
	}
    }

  if (1)
    {
      KEY *kp ;
      int i, n, iMax = keySetMax (brks) ;
      
      for (i = n = 0, kp = arrayp (brks, i, KEY) ; i < iMax ; i++, kp++)
	n += *kp ;
     
      wego_log (hprintf (umi->h, "Selected %d brks\n", n)) ;
    }
  ac_free (h) ;
} /* umiFuseLaneHits */

/*************************************************************************************/

static void umiGetHits (UMI *umi)
{
  AC_HANDLE h = ac_new_handle () ;

  int ii, iMax = arrayMax (umi->lanes) ;

  if (umi->SAM_file || umi->hitFile)
    iMax = 2 ;
  if (umi->maxLanes > 0 && umi->maxLanes < iMax) iMax = umi->maxLanes + 1 ;

  umi->getLaneInChan = channelCreate (iMax, int, h) ;
  umi->getLaneOutChan = channelCreate (iMax, int, h) ;

  umi->geneDict = dictHandleCreate (128, umi->h) ;
  umi->genes = arrayHandleCreate (128, GF, umi->h) ;

  for (ii = 1 ; ii < iMax ; ii++)
    channelPut (umi->getLaneInChan, &ii, int) ;

  channelClose (umi->getLaneInChan) ;

  for (ii = 0 ; ii < 16 && ii < iMax - 1 ; ii++)
    wego_go (umiGetLaneHits, umi, UMI) ;      // launch several threads
  for (ii = 1 ; ii < iMax ; ii++)
    {
      int x = 0 ;
      channelGet (umi->getLaneOutChan, &x, int) ;
    }
   umiFuseLaneHits (umi) ;
   if (0) exit(0) ;

   ac_free (h) ;  
} /* umiGetHits  */

/*************************************************************************************/
/* we expect to read the file 
 * tmp/METADATA/$target.mrna_map_ln_gc_gene_geneid.txt
 * col 2: chrom:123-456
 * col 5: gene  one line per exon, we must cluster the coordinates 
 */

static void umiGetGeneMap (UMI *umi)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (umi->geneMapFileName, 0, h) ;

  umi->geneDict = dictHandleCreate (128, umi->h) ;
  umi->genes = arrayHandleCreate (20000, GENE, umi->h) ;
  if (! ai)
    {
      fprintf (stderr, "// FATAL ERROR : cannot open geneMapFile file %s\n// Try UMI -h\n", umi->geneMapFileName) ;
      exit (1) ;
    }
    aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      char *cp, *cq ;
      int a1 = 0, a2 = 0, gene = 0, chrom = 0 ;
      GENE *up ;
     
      cp = aceInWord (ai) ;
      if (! cp || ! *cp || *cp == '#' || *cp == '/')
	continue ;
      aceInStep (ai, '\t') ; cp = aceInWord (ai) ; /* chrom:123-456 */
      if (! cp || ! *cp)
	continue ;
      cq = strchr (cp, ':') ;
      if (! cq || ! *cq)
	continue ;
      *cq = 0 ;
      dictAdd (umi->geneDict, cp, &chrom) ;
      cp = cq + 1 ;
      cq = strchr (cp, '-') ;
      if (! cq || ! *cq)
	continue ;
       *cq = 0 ;
       a1 = atoi (cp) ;
       cp = cq + 1 ;
       if (! *cp)
	continue ;
       a2 = atoi (cp) ;
       aceInStep (ai, '\t') ; cp = aceInWord (ai) ; /* col3 */
       aceInStep (ai, '\t') ; cp = aceInWord (ai) ; /* col 4 */
       aceInStep (ai, '\t') ; cp = aceInWord (ai) ; /* col 5 : gene */
 
      if (! cp || ! *cp)
	continue ;

      /* register */
      dictAdd (umi->geneDict, cp, &gene) ;
      up = arrayp (umi->genes, gene, GENE) ;
      up->gene = gene ;
      up->chrom = chrom ;
      if (! up->a1)
	{ 
	  up->a1 = a1 ; up->a2 = a2 ; 
	}
      if (a1 < a2)
	{
	  if (up->a1 > a1)
	    up->a1 = a1 ;
	  if (up->a2 < a2)
	    up->a2 = a2 ;
	   up->b1 = up->a1 ; up->b2 = up->a2 ;
	}
      else
	{
	  if (up->a1 < a1)
	    up->a1 = a1 ;
	  if (up->a2 > a2)
	    up->a2 = a2 ;
	   up->b1 = up->a2 ; up->b2 = up->a1 ;
	}
    }

  ac_free (h) ;
  return ; 
} /*  umiGetGenemap  */

/*************************************************************************************/

static void umiGetLaneList (UMI *umi)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (umi->laneListFileName, 0, h) ;
  ACEOUT ao = 0 ;
  int nn = 0 ;

  if (umi->laneListFileName)
    ai = aceInCreate (umi->laneListFileName, 0, h) ;
  else if (umi->SAM)
    ai = aceInCreateFromText ("sam", 0, h) ;
  else if (umi->hitFile)
    ai = aceInCreateFromText ("hits", 0, h) ;

  if (! ai)
    {
      fprintf (stderr, "// FATAL ERROR : cannot open LaneList file %s\n// Try UMI -h\n", umi->laneListFileName) ;
      if (0) exit (1) ;
    }
  umi->laneDict = dictHandleCreate (128, umi->h) ;
  umi->lanes = arrayHandleCreate (128, LANE, umi->h) ; 
  if (umi->makeTest)
     ao = aceOutCreate (umi->makeTest, ".LaneList", FALSE, h) ;
 
  while (aceInCard (ai))
    {
      char * cp ;
      LANE *lane ;

      cp = aceInWord (ai) ;
      if (cp && *cp)
	{
	  dictAdd (umi->laneDict, cp, &nn) ;
	  lane = arrayp (umi->lanes, nn, LANE) ;
	  lane->lane = nn ;
	  if (ao)
	    aceOutf (ao, "%s.%d\n", umi->makeTest, nn) ;

	  lane->run = umi->run ;
	  lane->h = ac_new_handle () ;
	  lane->cloneDict = dictHandleCreate (10000, lane->h) ;
	  lane->clones = arrayHandleCreate (10000, CLONE, lane->h) ;
	  lane->brks = arrayHandleCreate (10000, KEY, lane->h) ;
	  lane->brkIndels = arrayHandleCreate (10000, BKID, lane->h) ;
	  lane->covers = arrayHandleCreate (10000, KEY, lane->h) ;
	  lane->aliLn = arrayHandleCreate (10000, KEY, lane->h) ;
	  lane->hits = arrayHandleCreate (10000, HIT, lane->h) ;
	}
    }
  
  if (! umi->debug && ai)
    fprintf (stderr, "// Found %d lanes in file %s\n"
	     , nn
	     , aceInFileName (ai)
	     ) ;

  ac_free (ao) ;
  ac_free (ai) ;
  ac_free (h) ;
  return ;
} /* umiGetLaneList */

/*************************************************************************************/

static void umiGetBrks (UMI *umi)
{
  AC_HANDLE h = ac_new_handle () ;
  KEY *kp, *wp ;
  KEYSET brks, brks2 ;
  KEYSET covers = umi->covers ;
  int minSnpFrequency = umi->minSnpFrequency ;
  int minSnpCount = umi->minSnpCount ;
  int minSnpCover = umi->minSnpCover ;
  int ii, jj, iMax ;

  brks = umi->brks ;
  brks2 = umi->brks = arrayHandleCreate (10000, KEY, umi->h) ;
  iMax = arrayMax (brks) ;

  /* select the important positions */
  if (iMax >  arrayMax (covers))
    iMax = arrayMax (covers) ;
  for (ii = 1, jj = 0, kp = arrayp (brks, ii, KEY), wp = arrayp (covers, ii, KEY)  ;  ii < iMax - 1 ; ii++, kp++, wp++)
    {
      if (
	  *kp &&
	  (*wp >= minSnpCover || *(wp - 1) >= minSnpCover || *(wp + 1) >= minSnpCover))
	{
	  register int i, m = 0 ;
	  for (i = - umi->microDelta ; i <= umi->microDelta ; i++)
	    m += *(kp + i) ;
	  if (m >= minSnpCount &&
	      200 * m >= minSnpFrequency * (*wp) /* 200 because these are not the final counts */
	      )
	    {
	      keySet (brks2, jj++) = ii ;
	      if (200 * m < minSnpFrequency * (*wp) )
		invokeDebugger() ;
	    }
	} 
    }
  if (umi->externalBrks)
    {
      KEYSET ks = umi->externalBrks ;
      for (ii = 0 ; ii < keySetMax (ks) ; ii++)
	keySet (brks2, jj++) = keySet (ks, ii) ;
    }
  keySetSort (brks2) ;

  ac_free (brks) ;
  wego_log (hprintf (umi->h, "Selected %d well covered brks\n", arrayMax (umi->brks))) ;
  ac_free (h) ;
} /* umiGetBrks  */

/*************************************************************************************/
/*     several formats are recognized
 *	 chr7:1234     
 *       chr7:1234:
 *       chr7  1234
 *	   one position per line\n"
 */	
static void umiGetExternalSnpList (UMI *umi)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (umi->externalSnpListFileName, 0, h) ;
  KEYSET ks ;
  BitSet ebs ;
  int jj = 0 ;
  int t1 = umi->t1 ? umi->t1 : umi->offset ;

  if (!ai)
    {
      fprintf (stderr, "// FATAL ERROR: cannot open -externalSnpList %s\n"
	       , umi->externalSnpListFileName) ;
      exit (1) ;
    }

  ks = umi->externalBrks = keySetHandleCreate (umi->h) ;
  ebs = umi->externalBitSet = bitSetCreate (12000000, umi->h) ;
  while (aceInCard (ai))
    {
      char cutter, *cq, *cp = aceInWord (ai) ;
      int x = 0 ;
      
      cp = aceInWordCut (ai, "\t:", &cutter) ;
      if (! cp || ! *cp || *cp == '#')
	continue ;
      
      
 
     cq = strchr (cp, ':') ;
      if (cq)
	{
	  *cq = 0 ;
	  if (umi->target && strcmp (cp, umi->target))
	    continue ;
	  cp = cq + 1 ;
	  if (! *cp || sscanf (cp, "%d:", &x) != 1)
	    continue ;
	}
      else
	{
	  if (umi->target && strcmp (cp, umi->target))
	    continue ;
	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;
	  if (!cp || ! *cp || sscanf (cp, "%d:", &x) != 1)
	    continue ;
	}
      if (! umi->t2 || (umi->t1 <= x && umi->t2 >= x))
	{
	  keySet (ks, jj++) = x - t1  + 1 ;
	  if (umi->debug) fprintf (stderr, "jj = %d x = %d\n", jj, x) ;
	  bitSet (ebs, x - umi->t1  + 1) ;
	}
      else
	continue ;
      
 

    }

  ac_free (ai) ;
  ac_free (h) ;
} /* umiGetExternalSnpList */

/*************************************************************************************/
/* using the genome, snail trail the indels */ 
static int umiSnailTrailSlippingIndels (UMI *umi)
{
  Array brkIndels = umi->brkIndels  ;
  KEYSET brks = umi->brks  ;
  int ii, jj, kk ;
  int iMax = arrayMax (brkIndels) ; 
  int kMax = arrayMax (brks) ; 
  int kNew = arrayMax (brks) ; 
  int Delta = umi->Delta ;
  int pos, dx ;
  Array dna = umi->dna ;
  char *cp, buf[256] ;
  
  arraySort (brkIndels, umiBrkidOrder) ;
  for (ii = jj = 0  ;  ii < iMax ; ii ++)
    {
      BKID *b = arrp (brkIndels, ii, BKID) ;
      BOOL found = FALSE ;
      pos = b->pos ;
      if (ii > 0 && pos == (b-1)->pos && b->daa == (b-1)->daa && b->dxx == (b-1)->dxx)
	continue ;
      if (1 && b->c1) 
	found =TRUE ;
      for (kk = 0 ; kk < kMax && ! found ; kk++)
	if (keySet (brks,  kk) == pos)
	  found = TRUE ;
      if ( !found)
	continue ;
      /* this indel has adequate coverage */
      for (kk = ii ; kk < iMax ; kk++)
	{
	  BKID *b1 = arrp (brkIndels, jj, BKID) ;
	  BKID *b2 = arrp (brkIndels, kk, BKID) ;
	  if (b2->pos != b->pos ||
	      b2->daa != b->daa ||
	      b2->dxx != b->dxx
	      )
	    break ;
	  if (jj < kk)
	    *b1 = *b2 ;
	  jj++;
	}
    }
  iMax = arrayMax (brkIndels) = jj ;
  if (0)
    umiShowBrkIndels (umi->brkIndels) ;
      
  for (ii = 0  ;  ii < iMax ; ii ++)
    {
      BKID *b = arrp (brkIndels, ii, BKID) ;
      pos = b->pos ;
      for (jj = ii + 1 ; jj < iMax ; jj++)
	{
	  BKID *b1 = arrp (brkIndels, jj, BKID) ;

	  if (b1->pos != b->pos)
	    break ;
	  if (b1->daa == b->daa && b1->dxx == b->dxx)
	    {
	      b->mult += b1->mult ;
	      b1->mult = 0 ;
	      if (!b->c1 && b1->c1)
		{ b->c1 = b1->c1 ; b->c2 = b1->c2 ; }
	    }	  
	}
    }
  for (ii = jj = 0  ;  ii < iMax ; ii ++)
    {
      BKID *b = arrp (brkIndels, ii, BKID) ;
      if (b->mult)
	{
	  BKID *b1 = arrp (brkIndels, jj, BKID) ;
	  if (jj < ii)
	    *b1 = *b ;
	  jj++ ;
	}
    }
  iMax = arrayMax (brkIndels) = jj ;

  for (ii = jj = 0  ;  ii < iMax ; ii ++)
    {
      BKID *b = arrp (brkIndels, ii, BKID) ;
      pos = b->pos ;
      dx = b->dxx > 0 ? b->dxx : - b->dxx ;
      
      /* even for an insert we look for homopolymer in the genome
       * because it does not change much and it is as dubious
       * and much easier to analyze 
       */
      if (dna && dx < 255 && pos > dx  && dx > 0 && pos < arrayMax (dna))
	{
	  int i ;
	  if (b->dxx > 0) keySet (brks, kNew++) = pos + dx ;
	  cp = arrp (dna, pos  -dx - 1, char) ;
	  for (i = 0 ; i < dx ; i++, cp++)
	    buf[i] = *cp ;
	  cp = arrp (dna, pos -2, char) ;
	  for (i = -1 ; *cp == buf[(1000 * dx + i) % dx] ; i--, cp--)
	    if ((-i) % (Delta/2) == 0)
	      keySet (brks, kNew++) = pos + i ;
	  keySet (brks, kNew++) = pos + i ;
	}
    }
    
  keySetSort (brks) ;
  keySetCompress (brks) ;
  
  return arrayMax (brks) ;
} /* umiSnailTrailSlippingIndels */
  
/*************************************************************************************/

static void umiGetSegs (UMI *umi)
{
  AC_HANDLE h = ac_new_handle () ;
  int Delta = umi->Delta ;
  KEY *kp ;
  SEG *seg ;
  KEYSET brks = umi->brks  ;
  KEYSET nnn = keySetHandleCreate (h) ;
  KEYSET lns = keySetHandleCreate (h) ;
  Array segs ;
  Array brkIndels = umi->brkIndels  ;
  int ii, jj = 0, iMax = arrayMax (brks) ;
  int n, a, a1 = 0, a2 = 0 ;
  int dx = umi->Delta ;

  segs = umi->segs = arrayHandleCreate (10000, SEG, umi->h) ;

  
  if (dx > 8)
    dx = 8 ;

   /* create quasi continuous segs of brks */
  if (iMax)
    for (ii = jj = n = 0, kp = arrayp (brks, ii, KEY)  ;  ii < iMax + 1 ; ii++, kp++)
      {
	a = *kp ; /* do not add umi->t1, the seg is in local coords */
	if (!a && ii < iMax) continue ;
	if (ii == iMax || a - a2 > Delta)
	  {
	    if (a1)
	      {
		seg = arrayp (segs, jj++, SEG) ;
		seg->a1 = a1  ;
		seg->a2 = a2  ;
		seg->n = n ;
		
  if (seg->a1 >800000 && seg->a1 < 1000000)
    fprintf (stderr, "getsegs seg    %d %d\t*%d %d\n", seg->a1, seg->a2, seg->b1, seg->b2) ;

		seg->b1 = seg->a1 - dx ; 
		if (seg->b1 < 1)
		  seg->b1 = 1 ;
		seg->b2 = seg->a2 + dx ;
		
		keySet (nnn, n)++ ;
		keySet (lns, a2 - a1 + 1)++ ;
	      }
	    a1 = a ; n = 0 ;
	}
	if (! n++) a1 = a ;
	a2 = a ;
      }
  
  
  if (0)
    for  (ii = 0 ; ii < arrayMax (segs) ; ii++)
      {
	seg = arrp (segs, ii, SEG) ;
	fprintf (stderr, "getsegA seg(%d)  %d %d  %d %d brkIndels=%p max=%d\n", ii, seg->b1, seg->a1, seg->a2, seg->b2, brkIndels, arrayMax (brkIndels)) ;
      }
 
  if (brkIndels && arrayMax (brkIndels))
    {
      int ib, ibMax = arrayMax (brkIndels) ;

      for (ib = 0 ; ib < ibMax ; ib++)
	{
	  BKID *b = arrp (brkIndels, ib, BKID) ;  

	  if (b->c1)
	    fprintf (stderr, "=== getSegs brkIndels ibMax=%d b[%d].da=%d\n", ibMax, ib, b->daa) ;
	   
	  if (b->daa > Delta)
	    { /* create a new seg */
	      a1 = b->pos + b->protect ;
	      a2 = a1 + b->daa ;
	      n = 2 ;
	      /* if needed extend */
	      if (iMax)
		{
		  for (ii = 0, kp = arrayp (brks, ii, KEY)  ;  ii < iMax + 1 ; ii++, kp++)
		    if (*kp > a2)
		      {
			if (*kp < a2 + Delta) 
			  a2 = *kp ;
			else
			  break ;
		      }
		  for ( ; ii >= 0 ; ii--, kp--)
		    if (*kp < a1)
		      {
			if (*kp > a1 - Delta) 
			  a1 = *kp ;
			else
			  break ;
		      }
		}
	      seg = arrayp (segs, jj++, SEG) ;
	      seg->a1 = a1  ;
	      seg->a2 = a2  ;
	      seg->n = n ;
	      seg->b1 = seg->a1 - dx ; 
	      if (seg->b1 < 1)
		seg->b1 = 1 ;
	      seg->b2 = seg->a2 + dx ;
	      if (b->dxx > 0) 
		seg->b2 += b->dxx ;
	      seg->protect = b->protect ;
	      seg->c1 = b->c1 ;
	      seg->c2 = b->c2 ;

 	      keySet (nnn, n)++ ;
	      if (a2 < a1 + 120)
		keySet (lns, a2 - a1 + 1)++ ;
	    }

	  else if (b->dxx > 0) 
	    {
	      int j ;
	      for (j = 0 ; j < jj ; j++)
		{
		  seg = arrayp (segs, j, SEG) ;
		  if (seg->a1 <= b->pos && seg->a2 >= b->pos)
		    {
		      if (seg->b1 > seg->a1 - umi->Delta - b->dxx)
			seg->b1 = seg->a1 - umi->Delta - b->dxx ;
		    }
		}
	    }
	}
    }
 

  if (0)
    for  (ii = 0 ; ii < arrayMax (segs) ; ii++)
      {
	seg = arrp (segs, ii, SEG) ;
	fprintf (stderr, "getsegB seg(%d)  %d %d  %d %d\n", ii, seg->b1, seg->a1, seg->a2, seg->b2) ;
      }
  

 /* elongate segs that endup in sliding mess */
  iMax = arrayMax (segs) ;
  for  (ii = 0 ; ii < iMax ; ii++)
    {
      const register char *cp ;
      int b1, b2, i, dx ;
      Array dna = umi->dna ;
      int dMax = dna ? arrayMax (dna) : 0 ;

      if (! dna)
	continue ;
      seg = arrp (segs, ii, SEG) ;
      b1 = seg->b1 ; b2 = seg->b2 ;

      for (dx = 30 ; dx >= 1 ; dx--)
	{
	  int k ;
	  for (k = 0, i = b1, cp = arrp (dna, i - 1 , char) ; i - 1 - k - dx > 0 ; k++, cp--)
	    {
	      if (cp[0] != cp[-dx])
		break ;
	    }
	  if (k < dx)
	    continue ;
	  if (b1 > seg->b1 - k)
	    b1 = seg->b1 - k ;
	}
      seg->b1 = b1 ;
      if (seg->b1 < 1)
	seg->b1 = 1 ;

      for (dx = 30 ; dx >= 1 ; dx--)
	{
	  int k ;

	  if (b2 > arrayMax (dna))
	    b2 = arrayMax (dna) ;
	  for (k = 0, i = b2, cp = arrp (dna, i - 1 , char) ; i - 1 + k + dx < dMax ; k++, cp++)
	    {
	      if (cp[0] != cp[dx])
		break ;
	    }
	  if (k < dx)
	    continue ;
	  if (b2 < seg->b1 + k)
	    b2 = seg->b2 + k ;
	}
      seg->b2 = b2 ;

      if (0)
	fprintf (stderr, "getsegsC seg(%d)  %d %d  %d %d\n", ii, seg->b1, seg->a1, seg->a2, seg->b2) ;
    }	  

  iMax = arrayMax (segs) ;
  if (0)
    {
      for  (ii = jj = 0 ; ii < iMax ; ii++)
	{
	  seg = arrp (segs, ii, SEG) ;
	  a1 = seg->a1 ; a2 = seg->a2 ;
	  
	  if (a1 >= 823625 && a2 > a1 + 50 && a2 <= 825466)
	    arr (segs, jj++, SEG) = *seg ;
	}
      iMax = arrayMax (segs) = jj ;
    }

  for (jj = 0 ; jj < 2 ; jj++)
    {
      ACEOUT ao ;
      char *suffix = jj ? ".segN_histo.txt" : ".segLength_histo.txt" ;
      KEYSET ks = jj  ? nnn : lns ;
      int scumul = 0, cumul = 0 ;

      ao = aceOutCreate (umi->outFileName, suffix, FALSE, h) ;
      aceOutDate (ao, "###", hprintf (h, "Histogram of variable segment %s Delta = %d", jj ? "number of SNPs" :  "lengths in bp", umi->Delta)) ; 
      if (! jj)
	aceOutf (ao, "# Length Delta = %d\tNumber of variable segments\n", umi->Delta) ;
      else
	aceOutf (ao, "# Number of SNPs per fragment Delta = %d\tNumber of variable segments\tNumber of SNPs in these fragments\tCumul\n", umi->Delta) ;
      for (ii = 0 ; ii < keySetMax (ks) ; ii++)
	{
	  int k = keySet (ks, ii) ;
	  cumul += k ;
	  scumul += k * ii ;
	  if (jj) 
	    aceOutf (ao, "%d\t%d\t%d\t%d\n", ii, k, k * ii, scumul) ;
	  else
	    aceOutf (ao, "%d\t%d\n", ii, k) ;
	}
      aceOutf (ao, "Cumul\t%d\n", cumul) ;
      ac_free (ao) ;
    }

  for  (ii = 0 ; ii < arrayMax (segs) ; ii++)
    {
      seg = arrp (segs, ii, SEG) ;
      if (seg->c1 && seg->b1 > seg->c1 - dx)
	seg->b1 = seg->c1 - dx ; 
      if (seg->b1 < 1)
	seg->b1 = 1 ;
      if (seg->c2 && seg->b2 < seg->c2 + dx)
	seg->b2 = seg->c2 + dx ; 
    }

  iMax = arrayMax (segs) ;
  wego_log (hprintf (umi->h, "Selected %d variable segments\n", iMax)) ;
  if (0)
    for  (ii = 0 ; ii < iMax ; ii++)
      {
	seg = arrp (segs, ii, SEG) ;
	fprintf (stderr, "getsegsD seg(%d)  %d %d  %d %d  c1/c2 = %d %d\n", ii, seg->b1, seg->a1, seg->a2, seg->b2, seg->c1, seg->c2) ;
      }
 
  ac_free (h) ;
  return ;
} /* umiGetSegs  */

/*************************************************************************************/

static void umiGetTestOffset (UMI *umi)
{
if (umi->runTest) 
     {
       AC_HANDLE h = ac_new_handle () ;
       ACEIN ai = aceInCreate (hprintf (h, "%s.t1", umi->runTest), 0, h) ;

       while (aceInCard (ai))
	 {
	   char *cp = aceInWord (ai) ;
	   if (cp && !strcmp (cp,  "##OFFSET"))
	     {
	       aceInStep (ai, '\t') ;
	       aceInInt (ai, &(umi->offset)) ;
	       aceInStep (ai, '\t') ;
	       if (! umi->t2) aceInInt (ai, &(umi->t2)) ;
	       break ;
	     }
	   ac_free (ai) ;
	   ac_free (h) ;
	 }
     }
} /* umiGetTestT1 */

/*************************************************************************************/

static void umiInit (UMI *umi)
{
  int nn = umi->t2 - umi->t1 + 1 ;
  umi->doneChan = channelCreate (12, int, umi->h) ;
  
  umiGetTestOffset (umi) ;

  if (umi->externalSnpListFileName)
    umiGetExternalSnpList (umi) ;

  if (! nn) nn = 10000000 ;
  umi->dna = arrayHandleCreate (nn, char, umi->h) ;
  umi->dnaR = arrayHandleCreate (nn, char, umi->h) ;
  
  if (umi->geneMapFileName) 
    umiGetGeneMap (umi) ;
  umiGetLaneList (umi) ;
  wego_max_threads (umi->max_threads) ;
  /* wego_go (umiGetTargetFasta, umi, UMI) ;   will channelPut (umi->doneChan) */
  umiGetTargetFasta (umi) ;  /* will channelPut (umi->doneChan) */
 
  if (! umi->SAM) /* in hit mode we read the errors before we read the genome */
    {
      umiGetHits (umi) ;
      umiGetBrks (umi) ;  
    }
  if  (! channelGet (umi->doneChan, &nn, int) || nn < 0)
   {
     fprintf (stderr, "umiInit received a bad value %d\n", nn) ;
     exit (1) ;
   }
  else
    umi->t1 = nn ;

  if (umi->SAM) /* in SAM mode the errors only exist if we know the genome */
    {
      messcrash ("I do not yet know wnere the UMI-BAR-code is in the BAM file") ;
      umiGetHits (umi) ;
      umiGetBrks (umi) ;  
    }
  /* using the genome, snail trail the indels */ 
  if (1) umiSnailTrailSlippingIndels (umi) ;
  if (umi->dna && arrayMax (umi->dna)) umiGetSegs (umi) ;  
  
  return ;
} /* umiInit */

/*************************************************************************************/

static void umiCleanUp (UMI *umi)
{
  int ii, iMax = arrayMax (umi->lanes) ;
  LANE *lane ;

  for (ii = 0, lane = arrp (umi->lanes, 0, LANE) ; ii < iMax ; ii++, lane++)
    ac_free (lane->h) ;
  return ;
} /* umiCleanUp */
  
/*************************************************************************************/

static void importSra (const char *sra)
{
  AC_HANDLE h = ac_new_handle () ;

  char *s = hprintf (h, "fastq-dump --split-files --fasta 0 --gzip %s &", sra) ;
  fprintf (stderr, "Launched SRA download %s\n", s) ;
  system (s) ;

  ac_free (h) ;
  return ;
}

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/
/* run -run NA12878MOD -laneList tmp/TSNP/NA12878MOD/LaneList -t 20 -target_fasta TARGET/CHROMS/hs.chrom_20.fasta.gz -t1 25000001 -t2 30010000 -minSnpFrequency 18 -minSnpCover 10 -minSnpCount 4 -target_class Z_genome -o tata -maxLanes 4
 */

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: UMI -run runName ..... \n"
	   "//      try: -h --help \n"
	   "// Example\n"
	   "//     UMI -run SRR123456 -t chr3 -f chr3.fasta\n"
	   "// Several sets of parameters must or may be specified\n"
	   "// READS TO BE ANALYZED: \n"
	   "//   -run runName : [mandatory]\n"
	   "//   -laneList fileName :  [optional] list of lanes to be analyzed\n"
	   "//       [defaults to $fastcDir/$run/LaneList]\n"
	   "//    example:\n"
	   "//        -run SRR123456 -laneList Fastc/SRR123456/LaneList\n"
	   "//   -SAM : SAM format\n"
	   "//     the sam files are $hitDir/$lane.sam.gz\n"
	   "//   -SAM_file file_name : the single input alignment files are in SAM format\n"
	   "//   -hitFile file_name : a single hit file exported by clipalign\n"
	   "//    example:\n"
	   "//        -run SRR123456 -hitFile tmp/COUNT/SRR123456/f2.3.hits.gz\n"
	   "// TARGET ZONE:\n"
	   "//   -target_fasta fileName : the name of the reference fasta file, possibly gzipped \n"
	   "//   -target_class : target class (A_mito ... Z_genome) [default Z_genome]\n"
	   "//   -t targetName :  [optional] often a chromosome name\n"
	   "//   -t1 int -t2 int :  [optional] analyze only section [t1,t2] of the target\n"
	   "//    example:\n"
	   "//        -target_class  Z_genome -t chr7 -t1 2000000 -t2 3000000\n"
	   "//       Typically t2 - t1 = 10 Mbases, or if t1, t2 are not specified\n"
	   "//       the whole chromosome is analyzed .\n"
	   "// SEQUENCE FILES: \n"
	   "//   -fastcFile file_name : a single fasta/fastc file for the short reads\n"
	   "//   -fastcDir dirName : the directory containing the .fastc.gz files \n"
	   "//      [optional, default Fastc/$run]\n"
	   "//   -hitDir dirName : the directory containing the .hits.gz files \n"
	   "//       [optional, default tmp/COUNT/$run]\n"
	   "//       The program expects that for each lane listed in the laneList file\n"
	   "//       there exists a file called\n"
	   "//          $hitDir/$lane.hits.gz      containingg the magic hits\n"
	   "//          $fastcDir/$lane.fastc.gz   containing the read sequences]n"
	   "// SNP FILTERS\n"
	   "//   -method arbitrary_name: the method is echoed in the output file\n" 
	   "//   -externalSnpList filename: [optional] additional list of positions to be analysed\n"
	   "//      several formats are recognized\n"
	   "//          chr7:1234     chr7:1234:....   chr7\\t1234\n"
	   "//          one position per line\n"
	   "//   -ddx integer : default 1, must be in [0, 3] interval\n"
	   "//      It does not matter if the positions are off by +- 1 base\n"
  	   "//   -minSnpCover integer : min coverage [default 10] \n"
	   "//   -minSnpCount integer: [default 4]\n"
	   "//   -minSnpFrequency integer: [default 18] minimal MAF precentage\n"
	   "//   -minOverhang integer: [default 12] collect positions of incomplete alignments\n"
	   "//   -intron_only : just detect and report introns\n"
	   "//   -min_intron <int> : [default 30]  just detect and report introns\n"
	   "//   -max_intron <int> : [default 100000]  just detect and report introns\n"
	   "// SEGMENT CLUSTERING\n"
	   "//   -dx int ; [default 5] minimal acceptable value 5\n"
	   "//      A new phased segment starts when the next SNV is furher than dx\n"
	   "// TESTING\n"
	   "//    -make_test prefix : create a test over a small region\n" 
	   "//        prefix is a file prefix, for example test1 or TEST/test2\n"
	   "//        in the second case, the directory TEST must be created in adavance\n"
	   "//        After collecting the data from the usual parameters -t -t1 -t2 -run -fasctDir ...\n"
	   "//        the code will export prefix.fastc prefix.target corresponding to the region\n"
	   "//        You may then run the test quickly while adjusting -t1 -t2 -dx ... as wished\n"
	   "//    -run_test same_prefix : run a test prepared using -make_test\n"
	   "// OUTPUT\n"
	   "//     The program exports a vcf file and several histograms\n"
	   "//   -o fileNamePrefix : output file name, equivalent to redirecting stdout\n"
  	   "//   -gzo : the output files should be gzipped\n"
	   "//   -debug : more diagnostics\n"
	   ) ;

  
  if (argc > 1)
    {
      fprintf (stderr,
	       "//\n// You said: %s\n", commandBuf) ;
      fprintf (stderr,
	       "// ########## ERROR: I do not understand the argument%s ", argc > 2 ? "s" : "") ;
      for (i = 1 ; i < argc ; i++)
	fprintf (stderr, "%s ", argv[i]) ;
      fprintf (stderr, "\n") ;
    }
  exit (1) ;
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  UMI umi ;
  AC_HANDLE h = 0 ;
  char commandBuf [4000] ;

  freeinit () ; 
  messErrorInit (argv[0]) ;

  if (sizeof(void*) == 8)
    umi.is64 = TRUE ;

  h = ac_new_handle () ;
  memset (&umi, 0, sizeof (UMI)) ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
  umi.h = h ;

  {
    const char *sra = 0 ;
    if (getCmdLineOption (&argc, argv, "-SRA", &sra))
      {
	importSra (sra) ;
	exit (0) ;
      }
  }

  if (argc < 2)
    usage (commandBuf, argc, argv) ;
  if (getCmdLineBool (&argc, argv, "-help") ||
      getCmdLineBool (&argc, argv, "--help")||
      getCmdLineBool (&argc, argv, "-h")
      )
    usage (commandBuf, 1, argv) ;
  if (getCmdLineBool (&argc, argv, "-v") ||
      getCmdLineBool (&argc, argv, "--version")
      )
    {
      fprintf (stderr, "UMI: %s\n", VERSION) ;
      exit (0) ;
    }
 
  {
    int ix ;
    char *cp ;
    for (ix = 0, cp = commandBuf ;  ix < argc && cp + strlen (argv[ix]) + 1 < commandBuf + 3900 ; cp += strlen (cp), ix++)
      sprintf(cp, "%s ", argv[ix]) ;
  }

  /* consume optional args */
  umi.gzi = getCmdLineBool (&argc, argv, "-gzi") ;
  umi.gzo = getCmdLineBool (&argc, argv, "-gzo") ;
  umi.debug = getCmdLineBool (&argc, argv, "-debug") ;

  umi.intron_only = getCmdLineBool (&argc, argv, "-intron_only") ;

  getCmdLineOption (&argc, argv, "-o", &(umi.outFileName)) ;
 
  umi.SAM = getCmdLineBool (&argc, argv, "-SAM") ;
  if (getCmdLineOption (&argc, argv, "-SAM_file", &(umi.SAM_file)))
    umi.SAM = TRUE ;
  if (getCmdLineOption (&argc, argv, "-samFile", &(umi.SAM_file)))
    umi.SAM = TRUE ;
  getCmdLineOption (&argc, argv, "-hitFile", &(umi.hitFile)) ;
  getCmdLineOption (&argc, argv, "-fastcFile", &(umi.fastcFile)) ;
  getCmdLineOption (&argc, argv, "-geneMap", &(umi.geneMapFileName)) ;
  umi.geneFusion = getCmdLineBool (&argc, argv, "-geneFusion") ;

  /* TARGET ZONE */
  umi.target_class = 0 ; /* "Z_genome" ; */
  getCmdLineOption (&argc, argv, "-target_class", &(umi.target_class)) ;
  getCmdLineOption (&argc, argv, "-t", &(umi.target)) ;
  getCmdLineOption (&argc, argv, "-make_test", &(umi.makeTest)) ;

  getCmdLineOption (&argc, argv, "-method", &(umi.method)) ;
  
  if (getCmdLineOption (&argc, argv, "-run_test", &(umi.runTest)))
    {
      umi.targetFileName = hprintf (h, "%s.target.fasta", umi.runTest) ;
      umi.laneListFileName = hprintf (h, "%s.LaneList", umi.runTest) ;
      umi.run = umi.runTest ;
      umi.target = umi.runTest ;
      umi.hitDir = umi.fastcDir = "." ;
      umi.target_class = 0 ;
    }

  getCmdLineOption (&argc, argv, "-lane", &(umi.lane)) ;
  getCmdLineOption (&argc, argv, "-run", &(umi.run)) ;
  if (! umi.run)
   {
     fprintf (stderr, "FATAL ERROR: Missing argument -run, please try\n\tUMI -help\n") ;
     exit (1) ;
   }
  

  if (! umi.target && ! umi.geneFusion)
    {
      fprintf (stderr, "FATAL ERROR: Missing argument -t target_name, please try\n\tUMI -help\n") ;
      exit (1) ;
    }

  umi.minIntron = 30 ; /* default */
  umi.maxIntron = 100000 ; /* default */

  getCmdLineInt (&argc, argv, "-min_intron", &(umi.minIntron)) ;
  getCmdLineInt (&argc, argv, "-max_intron", &(umi.maxIntron)) ;

  getCmdLineOption (&argc, argv, "-target_fasta", &(umi.targetFileName)) ;
  getCmdLineInt (&argc, argv, "-t1", &(umi.t1)) ;
  getCmdLineInt (&argc, argv, "-t2", &(umi.t2)) ;
  if (! umi.targetFileName && umi.makeTest)
    {
      fprintf (stderr, "FATAL ERROR: Missing argument -target_fasta while using make_test, please try\n\tUMI -help\n") ;
      exit (1) ;
    }
   
  if (umi.t2 < umi.t1)
    {
      fprintf (stderr, "FATAL ERROR: zone limits expeted t1 < t2, got %d > %d, please try\n\tUMI -help\n", umi.t1, umi.t2) ;
      exit (1) ;
    }
  if (umi.t1 < 0)
    {
      fprintf (stderr, "FATAL ERROR: zone limits expeted t1 > 0 got %d  please try\n\tUMI -help\n", umi.t1) ;
      exit (1) ;
    }
  /* READ SET */
  if (! umi.SAM_file && ! umi.hitFile && ! umi.runTest && ! umi.lane)
    {
      if  (umi.run)
	{ /* usual magic configuration defaults */
	  umi.laneListFileName = hprintf (h, "Fastc/%s/LaneList", umi.run) ;
	  umi.fastcDir = "Fastc" ;
	  umi.hitDir = "tmp/COUNT" ;
	}
      
      getCmdLineOption (&argc, argv, "-laneList", &(umi.laneListFileName)) ;
      getCmdLineOption (&argc, argv, "-hitDir", &(umi.hitDir)) ;
      getCmdLineOption (&argc, argv, "-fastcDir", &(umi.fastcDir)) ;
    }

  getCmdLineOption (&argc, argv, "-externalSnpList", &(umi.externalSnpListFileName)) ;

  if (! umi.laneListFileName && ! umi.hitFile && ! umi.SAM_file)
    fprintf (stderr, "FATAL ERROR: Missing argument -LaneList, please try\n\tUMI -help\n") ;
	
  
  /* SNP FILTERS */
  umi.minSnpCover = 10 ;
  umi.minSnpCount = 4 ;
  umi.minSnpFrequency = 10 ;
  umi.minOverhang = 12 ;
  umi.Delta = 5 ;
  getCmdLineInt (&argc, argv, "-maxLanes", &(umi.maxLanes)) ;
  getCmdLineInt (&argc, argv, "-dx", &(umi.Delta)) ;
  if (umi.Delta < 5) 
    {
      fprintf (stderr, "-dx %d,  should be at least 5 bp\n", umi.Delta) ;
      exit (1) ;
    }
  
  umi.microDelta = 1 ;
  getCmdLineInt (&argc, argv, "-ddx", &(umi.microDelta)) ;
  if (umi.microDelta > 4) 
    {
      fprintf (stderr, "-ddx %d,  should be in [0,3] bp  interval\n", umi.microDelta) ;
      usage (commandBuf, argc, argv) ;
    }
 
  getCmdLineInt (&argc, argv, "-minOverhang", &(umi.minOverhang)) ;
  getCmdLineInt (&argc, argv, "-minSnpCover", &(umi.minSnpCover)) ;
  getCmdLineInt (&argc, argv, "-minSnpCount", &(umi.minSnpCount)) ;
  getCmdLineInt (&argc, argv, "-minSnpFrequency", &(umi.minSnpFrequency)) ;
  
  umi.max_threads = 8 ;
  getCmdLineInt (&argc, argv, "-max_threads", &umi.max_threads) ;
  
  if (umi.minSnpCount > umi.minSnpCover)
    umi.minSnpCount =  umi.minSnpCover ;
  if (argc != 1)
    {
      fprintf (stderr, "unknown argument, sorry\n") ;
      usage (commandBuf, argc, argv) ;
    }
  /* check file existence */
  if (umi.SAM_file)
    {
      ACEIN ai = aceInCreate (umi.SAM_file, 0, h) ; 
      if (! ai)
	messcrash ("Cannot open SAM_file %s\n", umi.SAM_file) ;
      ac_free (ai) ;
    }
  if (umi.hitFile)
    {
      ACEIN ai = aceInCreate (umi.hitFile, 0, h) ; 
      if (! ai)
	messcrash ("Cannot open hitFile %s\n", umi.hitFile) ;
      ac_free (ai) ;
    }
  if (umi.fastcFile)
    {
      ACEIN ai = aceInCreate (umi.fastcFile, 0, h) ; 
      if (! ai)
	messcrash ("Cannot open hitFile %s\n", umi.hitFile) ;
      ac_free (ai) ;
    }
  if (umi.laneListFileName)
    {
      ACEIN ai = aceInCreate (umi.laneListFileName, 0, h) ; 
      if (! ai)
	messcrash ("Cannot open laneList %s\n", umi.laneListFileName) ;
      ac_free (ai) ;
    }
 if (umi.externalSnpListFileName)
    {
      ACEIN ai = aceInCreate (umi.externalSnpListFileName, 0, h) ; 
      if (! ai)
	messcrash ("Cannot open externalSnpList %s\n", umi.externalSnpListFileName) ;
      ac_free (ai) ;
    }
  /* parallelization */
   if (umi.max_threads < 4)
     umi.max_threads = 4 ;
 
  /* check the absolute args */

  umiInit (&umi) ;
  if (umi.geneFusion && arrayMax (umi.geneFusions))
    umiExportGeneFusions (&umi) ;

  if (umi.segs)
    {
      umiAnalyze (&umi) ;
      umiMergeSubSegs (&umi) ;
      umiRationalize (&umi) ; 
      umiAliLn (&umi) ;
      umiExport (&umi) ;
      umiCleanUp (&umi) ;
    }
  wego_flush () ;

  ac_free (umi.ao) ;
  if (1)
    { 
      int mx ;
      messAllocMaxStatus (&mx) ; 
      fprintf (stderr, "// minSnpCount %d minSnpCover %d minSnpFrequency %d SNP\n"
	       , umi.minSnpCount
	       , umi.minSnpCover
	       , umi.minSnpFrequency
	       ) ;
      fprintf (stderr, "// SNP detected %d SEG analysed %d SNP reported %d\n"
	       , umi.brks ? arrayMax (umi.brks) : 0
	       , umi.segs ? arrayMax (umi.segs) : 0
	       , umi.snpExported
	       ) ;
      fprintf (stderr, "// %s done, %d errTrackings max memory %d Mb\n", timeShowNow(),  umi.nErrTrackings, mx) ;
     }
  ac_free (umi.h) ;

  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderrn and to ensure all pipes are closed*/
  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
