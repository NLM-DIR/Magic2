 /*  File: tricoteur.c
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

 * SNP and InDEl DETECTION

 * Analyze the TCT alignments, the bar code is in column 20
 * Rationalize away the multiplets: same a1 coordinates and same barcode
 *    while creating the consensus of the match errors and overlaps
 * Count expression per gene and transcripts in all target classes
 *    add the hierarchic counts to detect the missing transcripts of the top hierarchy
 * Construct count and export all variations 
 * Locate all gene fusion mrna1, mrna2, coords, delta-x
 *   export the delta-x histogram
 * Construct the consesnus (bridge) for all SNPs and gene fusions
 * Rextend the partial reads and count the 'validations'

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

 * SNP VALIDATION

 * export 31mers, count them in all the original fasta/fastq/fastc sequencing files

 * SNP REPORT

 * using the database, export nice reports

 */

#define VERSION "1.1"

#define MALLOC_CHECK   
#define ARRAY_CHECK   

#include "ac.h"
#include "channel.h"
#include "bitset.h"

/* snp chaining distance */
#define DELTA 12
#define ABS(_x) ((_x) < 0 ? -(_x) : (_x))

static char *types[] = {"Wild_type", "Substitution", "Duplication", "Insertion", "Diminution", "Deletion",  "DelIns", "Intron gt_ag forward", "Intron gc_ag forward", "Intron at_ac forward", "Intron gt_ag reverse", "Intron gc_ag reverse", "Intron at_ac reverse", "Unkown", 0} ;
typedef enum { FASTA=0, FASTC, FASTQ, RAW} DNAFORMAT ;

typedef struct variant_callerStruct {
  AC_HANDLE h ;
  const char *outFileName ;
  BOOL gzi, gzo, debug ;
  DNAFORMAT inFormat ;
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
  int targetKey ; /* offset of tct->target in tct->targetDict */
  int offset, t1, t2 ;

  BOOL intron_only ;
  BOOL SAM ;
  BOOL tabular ;
  BOOL qc ;
  BOOL geneFusion ;
  BOOL wantDeUno ;
  BOOL wantBrks ;

  const char *SAM_file ;
  const char *hitFile ;
  const char *geneMapFileName ;
  const char *mrnaRemapFileName ;
  const char *splitMrnaFileName ;
  const char *rawFile ;
  const char *fastaFile ;
  const char *fastcFile ;
  const char *fastqFile ;
  const char *externalSnpListFileName ;
  const char *geneFusionFileName ;
  const char *geneFusionPositionsFileName ;
  Array mrnas, splitMrnas ;
  KEYSET mrna2aa, splitMrna2aa ; /* mrna offSet in array mrnas */
  BOOL snpCount, snpExtend ;
  int stranded ;
  int valWordLength, snpExtendLength, snpTails ;
  Associator valAss ;
  Array valSnps ;
  DICT *valSnpZones ;
  DICT *valSnpNames ;
  DICT *valSnpWords ;
  const char *valSnpListFileName ;
  const char *valSnpList ;

  int microDelta, Delta ;
  int minSnpCover, minSnpCount ;
  float minSnpFrequency ;
  int minOverhang ;
  int subSampling ;
  int minBridge ;
  int minIntron ;
  int maxIntron ;

  int minAli, minAliPerCent ;

  int max_threads ;
  int maxLanes ;
  int minGF ;

  int nAna ; /* default 4, number of analyzers launched in parallel wego */
  int nErrTrackings ;
  Array lanes ;
  DICT *laneDict ;
  DICT *geneDict ;
  DICT *targetDict ;
  DICT *deUnoDict ;
  DICT *overhangDict ;
  const char *makeTest ;
  const char *runTest ;

  BigArray deUno ;
  BigArray brkIndels ;
  Array geneFusions ;
  Array wig ;
  KEYSET covers, brks, externalBrks, aliLn, fusionHisto, geneSupport ;
  KEYSET targetStart, targetStop ;
  BitSet externalBitSet ;
  Array externalBridges ;
  BigArray hits ;
  Array segs, genes ;
  Array dna, dnaR ; /* target dna */

  int snpExported ;

  CHAN *getLaneInChan, *getLaneOutChan, *analyzeChan ;
  CHAN *doneChan ;
} TCT ;

typedef struct laneStruct {
  AC_HANDLE h ;
  const char *run ;
  int lane ;   /* offset in laneDict */
  Array tails[8] ;
  Array brkIndels ;
  Array geneFusions ;
  KEYSET covers, brks, aliLn, fusionHisto, geneSupport ;
  Array genes, clones,  hits, deUno ;
  Array wig, qc ;
  Array valSnpCountsF, valSnpCountsR ;    /* count the word and its complement (systematic sequencing error detection) */
  Array valSnpCountsFF, valSnpCountsRR ;  /* if -(anti)strand flip read (1)2 before counting (intron strand detection)  */
  DICT *cloneDict, *deUnoDict, *qcDict ;
} LANE ;

typedef struct geneFusionStruct {
  LANE *lane ;
  int run ;
  int clone ;  /* offset in lane->cloneDict */
  int g1, g2, g12 ; /* genes in tct->geneDict */
  int m1, m2 ; /* mrna in geneDict */ 
  int chrom1, chrom2 ; /* chrom in tct->geneDict */
  int t1, t2 ; /* targets (i.e. mRNA or chrom) in tct->geneDict lane->cloneDict */
  int seq, tag ;
  char wo[3] ; /* ++ +- -+ -- */ 
  BOOL isPair ;
  int score1, score2 ;
  int c1, c2, d1, d2 ; /* chain coords */
  int x1, x2, y1, y2 ; /* read coodrs */
  int a1, a2, b1, b2 ; /* target coords */
  int chromDistance, da, dx ;
  int overhang1, overhang2 ;
  int n0 ; /* support */
  int gap ;  /* positive if unaliagned gap in the read, negative if overlap (sliding is possible) */ 
  int d2donor ;     /* distance from end of ali in transcript 1 to nearest downstream  donor (or acceptor if anti-mrna) */
  int d2acceptor ;  /* distance from end of ali in transcript 2 to nearest upstream acceptor (or donor if anti-mrna) */
} GF ;

typedef struct cloneStruct {
  const char *run ;
  int clone ;  /* offset in lane->cloneDict */
  int lane ;   /* offset in tct->laneDict */
  Array dna1, dna2 ;
} CLONE ;

typedef struct geneStruct {
  const char *run ;
  int gene ;  /* offset in tct->geneDict */
  int chrom ;   /* offset in  tct->geneDict */
  int a1, a2 ; /* coords in gene */
  int b1, b2 ; /* non stranded coords */
} GENE ;

typedef struct deUnoStruct {
  int target ;
  int gene ;   /* offset in tct->geneDict */
  int type ;   /* full name of the variant */
  int tag ;   /* the acedb tag: DelA, A2G ...Multi_deletion 31 */
  int nt, ns ; /* number of supporting tags, distintct sequences  */
  int a1, a2, da, ddx ; /* absolute chrom coordinates, da is negative for an insert */
  int clone, x1, x2 ; /* coordinates of a long insert in its read */
  int insert ; /* offset in lane->deUnoDict: sequence of the insert i.e. atgc or intron feet */
  int dup ;
  BOOL isRead1, isDown, isExternal ;
} UNO ;

#define QCMAX 6
typedef struct qcStruct {
  int type, fmt ;
  long int x[QCMAX] ;
  int t[QCMAX] ;
  float z[QCMAX] ;
} QC ;

typedef struct hitStruct {
  int lane ;   /* offset in tct->laneDict */
  int clone ;  /* offset in lane->cloneDict */
  int gene ;   /* offset in tct->geneDict */
  int target  ; /* offset in tct->targetDict */
  int target_class ; /* offset in tct->targetDict */
  int score, mult, uu, a1, a2, b2, c1,c2, x1, x2 ;  /* absolute chrom coordinates */
  int overhang5p, overhang3p ;
  BOOL isRead1, isDown ;
  Array dna, cigarettes ;
} HIT ;

typedef struct localHitStruct {
  int target ;
  int lane ;   /* offset in laneDict */
  int clone ;  /* offset in lane->cloneDict */
  int mult, uu ;
  int a1, a2, x1, x2 ;   /* local coordinates (tct->t1 is substracted) */
  int b1, b2, y1, y2 ;   /* match b1 b2 of the seg */
  int c1, c2, z1, z2 ;   /* match b1 b2 in other part of same read */
  BOOL isDown, isRead1, isB ;
  Array dna, err ;
  int cc, pos ;
  int Lcc, Rcc ;
} LHIT ;

typedef struct segStruct {
  int target ;
  int a1, a2, b1, b2 ;  /* local coordinatess (tct->t1 is substracted) */
  int c1, c2, z1, z2 ;     /* match b1 b2 in other part of same read */
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

typedef struct valWordStruct {
  int zone, name, motif ;
} VW ;

typedef struct brkIndelStruct {
  int target, pos, mult, daa, dxx, c1, c2 ;
} BKID ;

typedef struct consensusStruct {
  int cumulsU, cumulsNu ;
  double proba ;
  Array dna ;
} CSS ;

#define INSMAX 3001
#define DELMAX 201

typedef struct bridgeStruct {
  int ii, jj, ln, dx, ddx, ddy, ds, dt, nDiff ;
  int type ;
  int a1, a2, b1, b2, c1, c2 ;

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
  char buf1[DELMAX], buf2[INSMAX] ;
} BRIDGE ;

static void tctAnalyzeSeg (const TCT *tct, SEG *seg) ;
static void tctDeUnoFuseOneLane (const TCT *tct, int myLane) ;

/*************************************************************************************/
 
static int tctGeneFusionOrder (const void *va, const void *vb)
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
} /* tctCssOrder */


/*************************************************************************************/
 
static int tctCssOrder (const void *va, const void *vb)
{
  const CSS *up = (const CSS *)va, *vp = (const CSS*)vb ;

  if (0)
    return up->proba >= vp->proba ? -1 : 1 ;
  else if (1)
    {
      int n = (up->cumulsU + up->cumulsNu) - (vp->cumulsU  + vp->cumulsNu) ;
      if (n >0) return -1 ;
      if (n < 0) return 1 ;
      return 0 ;
    }
  else
   return -up->cumulsU -up->cumulsNu + vp->cumulsU  + vp->cumulsNu ; 
} /* tctCssOrder */

/*************************************************************************************/
 
static int tctBrkidOrder (const void *va, const void *vb)
{
  const BKID *up = (const BKID *)va, *vp = (const BKID*)vb ;
  int n ;
  n = up->target - vp->target ; if (n)  return n ;
  n = up->pos -vp->pos ; if (n) return n ;
  n = up->daa - vp->daa ; if (n) return -n ;
  return 0 ;
} /* tctBrkidOrder */

/*************************************************************************************/
 
static int tctBaseCountOrder (const void *va, const void *vb)
{
  const BASECOUNT *up = (const BASECOUNT*)va, *vp = (const BASECOUNT*)vb ;
  return -up->n + vp->n ;
} /* baseOrder */

/*************************************************************************************/
 
static int tctA1Order (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;
  int n ;
  
  n = up->target - vp->target ; if (n)  return n ;
  n = up->a1 - vp->a1 ; if (n)  return n ;
  n = up->a2 - vp->a2 ; if (n)  return n ;
  return 0 ;
} /* a1Order */

/*************************************************************************************/
 
static int tctX1Order (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;
  int n ;
  
  n = up->clone - vp->clone ;  if (n)  return n ;
  n = up->isRead1 - vp->isRead1 ; if (n)  return -n ;
  n = (up->isDown ? up->x1 : up->x2) - (vp->isDown ? vp->x1 : vp->x2) ; if (n)  return n ;
  n = (up->isDown ? up->x2 : up->x1) - (vp->isDown ? vp->x2 : vp->x1) ; if (n)  return n ;
  n = up->target - vp->target ; if (n)  return n ;
  return 0 ;
} /* a1Order */

/*************************************************************************************/
 
static int tctSegOrder (const void *va, const void *vb)
{
  const SEG *up = (const SEG *)va, *vp = (const SEG *)vb ;
  int n ;
  
  n = up->target - vp->target ; if (n)  return n ;
  n = up->a1 - vp->a1 ; if (n)  return n ;
  n = up->a2 - vp->a2 ; if (n)  return n ;
  return 0 ;
} /* segOrder */

/*************************************************************************************/
 
static int tctUnoOrder (const void *va, const void *vb)
{
  const UNO *up = (const UNO *)va, *vp = (const UNO *)vb ;
  int n ;
  
  n = up->target - vp->target ; if (n)  return n ;
  n = up->a1 - vp->a1 ; if (n)  return n ;
  n = up->a2 - vp->a2 ; if (n)  return n ;
  n = up->type - vp->type ; if (n)  return n ;

  return 0 ;
} /* unoOrder */


/*************************************************************************************/
 
static int tctBridgeOrder (const void *va, const void *vb)
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
  
  n = up->target - vp->target ; if (n)  return n ;
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
  
  n = up->target - vp->target ; if (n)  return n ;
  n = up->a2 - vp->a2 ; if (n)  return n ;
  n = up->a1 - vp->a1 ; if (n)  return n ;
  return 0 ;
} /* a2Order */
#endif
/*************************************************************************************/

static void tctShowSegs (Array segs, const char *title)
{
  SEG *seg ;
  int ii, iMax = arrayMax (segs) ;

  if (! title) title = "===" ;
  for  (ii = 0 ; ii < iMax ; ii++)
    {
      seg = arrp (segs, ii, SEG) ;
      fprintf (stderr, "%s seg(%d)  %d %d  %d %d  c1/c2 = %d %d\n", title, ii, seg->b1, seg->a1, seg->a2, seg->b2, seg->c1, seg->c2) ;
    }
}

/*************************************************************************************/

static void tctShowLhits (Array hits)
{
  int ii, iMax = hits ? arrayMax (hits) : 0 ;
      
  for (ii = 0  ;  ii < iMax ; ii ++)
    {
      LHIT *up = arrp (hits, ii, LHIT) ;
      fprintf (stderr, "++++ ii=%d a1/a2=%d/%d\tx1/x2=%d/%d\tb1/b2=%d/%d\ty1/y2=%d/%d\tc1/c2=%d/%d  dnaLn=%d\n"
	       ,ii
	       , up->a1, up->a2
	       , up->x1, up->x2
	       , up->b1, up->b2
	       , up->y1, up->y2
	       , up->c1, up->c2
	       , up->dna ? arrayMax (up->dna) : 0
	       ) ;
    }
  return ;
}

/*************************************************************************************/

static void tctShowBridges (Array bridges)
{
  int ii, iMax = bridges ? arrayMax (bridges) : 0 ; 
      
  for (ii = 0  ;  ii < iMax ; ii ++)
    {
      BRIDGE *b = arrp (bridges, ii, BRIDGE) ;

      fprintf (stderr, "++++ ii=%d t=%d =%d/%d coverp=%d coverm=%d, mp/mm=%d/%d, wp/wm=%d/%d   ln/dx=%d/%d\n"
	       , ii
	       , b->type
	       , b->a1, b->a2
	       , b->coverp, b->coverm
	       , b->mp, b->mm
	       , b->wp, b->wm 
	       , b->ln, b->dx
	       ) ;
    }
  return ;
}

/*************************************************************************************/

static void tctShowBrkIndels (BigArray brkIndels)
{
  long int ii, iMax = bigArrayMax (brkIndels) ; 
      
  for (ii = 0  ;  ii < iMax ; ii ++)
    {
      BKID *b = bigArrp (brkIndels, ii, BKID) ;
      fprintf (stderr, "++++ ii=%ld pos=%d daa/dxx==%d/%d c1/c2=%d/%d mult=%d\n"
	       ,ii,b->pos, b->daa, b->dxx, b->c1,b->c2,b->mult) ;
    }
  tctShowBridges (0) ;
  tctShowLhits (0) ;
  return ;
}

/*************************************************************************************/
 
static Array  tctSeg2Hits (const TCT *tct, SEG *seg, AC_HANDLE h0) 
{
  int jj ;
  long int kk, kMax, dk ;
  HIT *up ;
  int t1 = tct->t1 ? tct->t1 - 1 : tct->offset ;
  int a1 = seg->a1 - tct->Delta + t1 ;
  int a2 = seg->a2 + tct->Delta + t1 ;
  Array aa = arrayHandleCreate (1024, LHIT, h0) ; 

  /* position by dichotomy */
  kMax = bigArrayMax (tct->hits) ; 
  if (! kMax)
    return aa ;
  kk = kMax - 1 ;
  up = bigArrp (tct->hits, kk, HIT) ;
  /* look by dichotomy in the a1 sorted array for the rightmost hit with non-decreasing b2 < searched a1 */
  if (up->b2 < a1)       /* impossible */
    return aa ;
  if (up->b2 < a1) /* if equal, we win, else we need dichotomy search */
    {
      kk = kMax/2 ; kk-- ; /* iMax = 7 or 8, position at 3 */
      for (dk = kk ; dk > 1 ; dk /= 2)
	{
	  up = bigArrp (tct->hits, kk, HIT) ;
	  if (kk > dk && up->b2 > a1)
	    kk -= dk ;
	  else
	    kk += dk ;
	}
    }
  for (up = bigArrp (tct->hits, kk, HIT) ; kk > 0 ; kk--, up--)
    {
      if (up->b2 < a1)
	break ;
     }

  /* collate the relevant hits */
  for (kk = jj = 0, up = bigArrp (tct->hits, kk, HIT) ; kk < kMax ; kk++, up++) 
    {
      int da = up->a2 - up->a1 ;
      
      if (da < 0) da = -da ;

      if (up->a1 > a2)
	break ;
      if (up->clone && up->a1 < a2 && up->a2 > a1 && (up->a1 < a1 + 2000 || up->a2 > a2 - 2000) )
	{
	  LHIT *vp = arrayp (aa, jj++, LHIT) ;
	  LANE *lane = arrp (tct->lanes, up->lane, LANE) ;
	  CLONE *clone = arrp (lane->clones, up->clone, CLONE) ;

	  vp->lane = up->lane ;
	  vp->clone = up->clone ;
	  vp->mult = up->mult ;
	  vp->uu = up->uu ;
	  vp->a1 = up->a1  - t1  ;  /* local coords */
	  vp->a2 = up->a2  - t1  ;  /* local coords */
	  vp->x1 = up->x1 ;
	  vp->x2 = up->x2 ;
	  vp->c1 = up->c1 ;
	  vp->c2 = up->c2 ;
	  if (vp->x1 < 0)
	    { vp->x1 = - up->x1 ; vp->isDown = TRUE ; }
	  if (vp->x2 < 0)
	    { vp->x2 = - up->x2 ; vp->isRead1 = TRUE ; }
	  vp->dna = vp->isRead1 ? clone->dna1 : clone->dna2 ;	  
	}
    }

  return aa ;
} /* tctSeg2Hits */

/*************************************************************************************/
/* Position segs in the iDown direction */
static void  tctCombSegPosition (const TCT *tct, SEG *seg, Array hits, BOOL isDown, AC_HANDLE h0)
{
  int kk = 0, ii, iMax = arrayMax (hits) ;
  LHIT *up ;
  int b1 = seg->b1  ;
  int b2 = seg->b2  ;
  int c1 = seg->c1  ;
  int c2 = seg->c2  ;
  int Ln = arrayMax (tct->dna) ;
  Array dnaLong = 0 ;
  int delta = tct->Delta < 9 ? tct->Delta - 1 : 8 ;
  BOOL debug = FALSE ;

  if (! isDown)
    {
      int dummy ;
      b1 = Ln - b1 + 1 ;
      b2 = Ln - b2 + 1 ;
      dummy = b1 ; b1 = b2 ; b2 = dummy ;  
       
      c1 = Ln - c1 + 1 ;
      c2 = Ln - c2 + 1 ;
      dummy = c1 ; c1 = c2 ; c2 = dummy ;  
       
      dnaLong = tct->dnaR ;
    }
  else
    dnaLong = tct->dna ;

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
      if (debug)
	{
	  char *strand = isDown ? "DOWN" : "UP" ;
	  if (! kk++) { fprintf (stderr, "TARGET_1    \t%s\t", strand) ; showDna (dnaLong, up->b1) ; }
	  fprintf (stderr, "READ_1 %03d >\t%s%s\t", up->mult, strand, up->isB ? "B" : "") ; showDna (up->dna, up->y1) ;
	}
    }
  return ;
} /* tctCombSegPosition */

/*************************************************************************************/
/* Comb segs in the iDown direction */
static BOOL tctCombSeg (const TCT *tct, SEG *seg, Array hits, BOOL isDown, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  int iMax = arrayMax (hits) ;
  int iiStep = (iMax)/2048 + 1 ;
  LHIT *up ;
  int pos = 0 ;
  int b1 = seg->b1 ;
  int b2 = seg->b2 ;
  int c1 = seg->c1 ;
  int c2 = seg->c2 ;
  double proba ;
  CSS *cssp, *cssp2 ;
  Array dna, css ;
  BOOL ok = TRUE ;
  int *cumulp ;
  Array cumuls ;
  Array nn = arrayHandleCreate (4, BASECOUNT, h0) ;
  int b2b[256] ;
  int minSnpCount = tct->minSnpCount ;
  int Ln = arrayMax (tct->dna) ;
  Array dnaLong = 0 ;
  BOOL debug = FALSE ;
  
  if (! isDown)
    {
      int dummy ;
      b1 = Ln - b1 + 1 ;
      b2 = Ln - b2 + 1 ;
      dummy = b1 ; b1 = b2 ; b2 = dummy ;  
      
      c1 = Ln - c1 + 1 ;
      c2 = Ln - c2 + 1 ;
      dummy = c1 ; c1 = c2 ; c2 = dummy ;  
       
      dnaLong = tct->dnaR ;
    }
  else
    dnaLong = tct->dna ;
  
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
  
  cssp = arrayp (css, 0, CSS) ;
  if (b2 < b1 - 120)
    goto failed ;
  dna = cssp->dna = arrayHandleCreate (b2 - b1 + 128, char, h0) ; 
  proba = cssp->proba = 1.0 ;
  cumuls = arrayHandleCreate (12, int, h) ;

  for (pos = 0 ; ok ; pos++) /* threading extension */
    {
      int cc, cMax = arrayMax (css) ;
      int cMaxOld = cMax ;
      
      arraySort (css, tctCssOrder) ;
      if (pos < 60 && cMax > 96)
	cMaxOld = cMax = 96 ;
      ok = FALSE ;
      for (cc = 0 ; cc < cMaxOld ; cc++) /* try to extend each existing consensus */
	{
	  BASECOUNT *bc ;
	  int ii, base, base2 ;
	  int nnn = 0 ;
	  float proba1 = 1 ;       

	  cssp = arrayp (css, cc, CSS) ;
	  dna = cssp->dna ;
	  cumulp = arrayp (cumuls, cc, int) ;
	  if (*cumulp < 0) 
	    continue ;
	  if (pos && *cumulp == 0) 
	    continue ;
	  for (base = 0 ; base < 4 ; base++)
	    arr (nn, base, BASECOUNT).n = 0 ;
	  for (ii = 0 ; ii < iMax ; ii += iiStep)
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
	  arraySort (nn, tctBaseCountOrder) ;
	  bc = arrp (nn, 0,  BASECOUNT) ;
	  array (dna, pos, char) = bc->base ; 
	  (*cumulp) += bc->n ;
	  
	  proba = cssp->proba ; proba1 = 0 ;
	  /* bifurcate */
	  for (base2 = 1 ; cMax < 128 && base2 < 4 ; base2++)
	    {
	      float proba2 ;
	      BASECOUNT *bc2 = arrp (nn, base2,  BASECOUNT) ;
	      if (bc2->n >= minSnpCount  && bc2->n * proba / (1.0 * nnn) > .0001)  
		{  
		  Array dna2 ;
		  int *cumul2p ;
		  dna2 = dnaHandleCopy (dna, h0) ;
		  
		  /* register and restore current pointers */
		  cssp2 = arrayp (css, cMax, CSS) ;
		  cssp = arrp (css, cc, CSS) ; /* may have been relocated */
		  cssp2->dna = dna2 ;
		  dna = cssp->dna ;
		  cumul2p = arrayp (cumuls, cMax, int) ;
		  cumulp = arrayp (cumuls, cc, int) ;
		  if (1)
		    {
		      proba2 = bc2->n /(1.0 * nnn) ;
		      proba1 += proba2 ;
		      cssp2->proba = proba * proba2 ;
		    }
		  else
		    {
		      cssp2->proba = proba * bc2->n /(1.0 * nnn) ;
		      cssp->proba *= bc->n /(1.0 * nnn) ;
		    }
		  
		  /* extend */
		  *cumul2p = *cumulp -bc->n + bc2->n ;
		  array (dna2, pos, char) = bc2->base ;
		  cMax++ ;
		  
		  /* assign the reads to the new consensus */
		  {
		    int ii2 ; 
		    int delta = tct->Delta < 9 ? tct->Delta - 1 : 8 ;

		    for (ii2 = 0 ; ii2 < iMax ; ii2 += iiStep)
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
	  cssp->proba *= (1.0 - proba1) ;
	}
    }
 
  if (1)
    {
      int ii ;
      
      for (ii = 0 ; ii < arrayMax (css) ; ii++)
	{
	  CSS *up = arrp (css, ii, CSS) ;
	  if (1 && array (cumuls, ii, int) < 0)
	    {
	      up->cumulsU *= -1 ;
	      up->cumulsNu *= -1 ;
	    } 
	}
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
    }

  arraySort (css, tctCssOrder) ;
  if (debug) 
    fprintf (stderr, "seg->a1/a2=%d/%d b1/b2=%d/%d\n"
	     ,seg->a1,seg->a2
	     ,seg->b1,seg->b2
	     ) ;
  if (debug) /* 2018_02_02 */
    {
      int cc, cMax = arrayMax (css) ;
      int u1, u2 ;
      
      u1 = seg->b1 + tct->t1 - 1 ;
      u2 = seg->b2 + tct->t1 - 1 ;
      
      fprintf (stdout, "\nTARGET %s %s b1=%d b2=%d\n"
	       , tct->target
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
	  fprintf (stdout, "CONSENSUS %s\t%d\t%d\t%d\t"
		   , isDown ? "Down" : "Up"
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
 failed:
  ac_free (h) ;  
  return 0 ;
} /* tctCombSeg */

/*************************************************************************************/
/* stitches left aligned segs with right aligned segs */
static void  tctMakeBridge (const TCT *tct, SEG *seg, Array hits, AC_HANDLE h0)
{
  BRIDGE *bb ;
  int ii, jj, iMax, jMax ;       
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
      if (! dna)
	return ;
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
  if (tct->dna && arrayMax (tct->dna) > seg->b2 - seg->b1 + 1 )
    {
      int i, ln = seg->b2 - seg->b1 + 1 ;
      
      for  (i = 0 ; i < ln ; i++)
	array (bufT, i, char) = arr (tct->dna, seg->b1 + i - 1 , char) ;
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
      int Dx = 2 * tct->Delta ;

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
	    }
	}
    }

  return ;
} /* tctMakeBridge */

/*************************************************************************************/
/* fuse identical bridges, add the wild type */
static void  tctFuseBridges (const TCT *tct, SEG *seg, Array hits)
{
  int ii, jj, iMax = seg->bridges ? arrayMax (seg->bridges) : 0 ;
  BRIDGE *bb, *bb0 ;
  int coverp, coverm ;
  int coverpnu, covermnu ;
  /*   Array bufT = 0 ; */
  /*   int ddd = 1 ; */
 
  if (! iMax) return ;

  /* remove duplicates */
  for (ii = 0, bb0 = arrp (seg->bridges, 0, BRIDGE) ; ii < iMax ; bb0++, ii++)
    {
      char *cp ;
      if (bb0->type < 0)
	continue ;
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
  for (ii = 0, bb0 = arrp (seg->bridges, 0, BRIDGE) ; ii < iMax && bb0->type ; bb0++, ii++)
    ; /* empty loop */

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

} /* tctFuseBridges  */

/*************************************************************************************/
/* count reads supporting all briges */
static void  tctBridgeSupport (const TCT *tct, SEG *seg, Array hits, AC_HANDLE h0)
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
	    ! up->dna ||
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

   arraySort (seg->bridges, tctBridgeOrder) ;
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

       if (ii == 0 && tct->externalBitSet &&
	   seg->a2 - seg->a1 < 20)
	 {
	   BitSet ebs = tct->externalBitSet ;
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
       if (! seg->isSubseg && seg->a2 > seg->a1 + 4 * tct->Delta)
	 {
	   ddd = 2 ;
	   if (ii < 2)
	     ddd = 10 ;
	 }
       
       if (ddd * (cover + covernu) < tct->minSnpCover) 
	 bb->type = -2 ;
       ff = 100.0 * (m) / cover ;
       
       ffnu = covernu ? 100.0 * (mnu) / covernu : 0 ;
       
       if (!((ddd*ff >= tct->minSnpFrequency && 2 * m > tct->minSnpCount) || (ddd*ffnu  >= tct->minSnpFrequency &&  2 * mnu > tct->minSnpCount)))
	  bb->type = -3 ;
     }

   arraySort (seg->bridges, tctBridgeOrder) ;
   for (ii= 0, bb = arrp (seg->bridges, ii, BRIDGE) ;  ii < iMax ; ii++, bb++)
       if (bb->type < 0)
	 { iMax = arrayMax (seg->bridges) = ii ; break ; } 

   ac_free (h) ;
} /* tctBridgeSupport */

/*************************************************************************************/
/* set common letters to lower */
static BOOL  tctBridgeSetUpper (const TCT *tct, SEG *seg, Array hits, AC_HANDLE h0)
{
  int i, ii ;
  int iMax = seg->bridges ? arrayMax (seg->bridges) : 0 ;
  BRIDGE *bb, *bb0 ;
  int ln, ddx, ddy ;
  char cc, *cp ;
  Array bufT ;
  BOOL ok = FALSE ;
  BOOL debug = FALSE ;

  if (iMax < 2) return ok ;
 
  for (ii = 0, bb0 = 0, bb = arrp (seg->bridges, 0, BRIDGE) ;  ii < iMax ; bb++, ii++)
    if (bb->type == 0)
      { bb0 = bb ; break ; }
  if (! bb0)
    return ok ;
  bufT = bb0->targetSnippet ;
  ln = arrayMax (bufT) ;

  for (i = ii = 1, bb = arrp (seg->bridges, 0, BRIDGE) ;  i < iMax ; bb++, i++)
    if (bb->type >= 0)
      { ok = TRUE ; ii = i + 1 ; }
  iMax = ii ;
  if (! ok) return ok ;

  if (debug)
    {
      fprintf (stderr,"SetUpper seg a1/a2 %d/%d\n", seg->a1, seg->a2) ;
      tctShowBridges (seg->bridges) ;
    }
  /* setup all letters to lower case */
  for (ii = 0, bb = arrp (seg->bridges, 0, BRIDGE) ;  ii < iMax ; bb++, ii++)
    if (bb->type >= 0)
      for (cp = arrp (bb->variantSnippet, 0, char) ; *cp ; cp++)
	*cp = ace_lower (*cp) ;
  for (cp = arrp (bufT, 0, char) ; *cp ; cp++)
     *cp = ace_lower (*cp) ;

   /* eatup the right bases common to all bridges */
   for (cc = 1, ddy = 0 ; cc && ddy < ln ; ddy++)
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
   /* upper the right common letters */
   for (i = 0, cp = arrp (bufT, ln - 1, char) ; *cp && i < ddy ; cp--, i++)
     *cp = ace_upper (*cp) ;
   for (ii = 0, bb = arrp (seg->bridges, 0, BRIDGE) ;  ii < iMax ; bb++, ii++)
     {
       bb->ddy = ddy ;
       if (bb->type >= 0)
	 for (i = 0, cp = arrp (bb->variantSnippet, bb->dx - 1, char)  ; i < ddy ; i++, cp--)
	   *cp = ace_upper (*cp) ;
       if (ii && ddy == bb->dx)
	 bb->type = -1 ;
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
   /* upper the left common letters */
   if (ddx + ddy > bb0->dx)
     ddx = bb0->dx - ddy ;
   if (ddx < 0) ddx = 0 ;/* should not occur */
   for (i = 0, cp = arrp (bufT, 0, char) ; *cp && i < ddx ; cp++, i++)
     *cp = ace_upper (*cp) ;
   for (ii = 0, bb = arrp (seg->bridges, 0, BRIDGE) ;  ii < iMax ; bb++, ii++)
     {
       bb->ddx = ddx ;
       if (bb->type >= 0)
	 for (cp = arrp (bb->variantSnippet, 0, char), i = 0 ; *cp && i < ddx ; cp++, i++)
	   *cp = ace_upper (*cp) ;
     }
   
   /* upper the central  common letters in each variant */
   for (i = ddx ; i < bb0->dx - ddy ; i++)
     {
       char cc1 ;
       cc = arr (bufT, i, char) ;
       cc1 = ace_upper (cc) ;
       for (ii = 1, bb = arrp (seg->bridges, ii, BRIDGE) ;  cc && ii < iMax ; bb++, ii++)
	 if (bb->type >= 0)
	   {
	     cc = arr (bufT, i, char) ;
	     if (i < bb->dx && cc == arr (bb->variantSnippet, i, char))
	       arr (bb->variantSnippet, i, char) = ace_upper (cc) ;
	     else 
	       cc1 = 0 ;
	   }       
       if (cc1)
	 arr (bufT, i, char) = cc1 ;
     }
   return ok ;
} /* tctBridgeSetUpper */

/*************************************************************************************/
/* duplicate the reads so that they come x2->x1 so that we can extend them in prime */
static void  tctAddReverseSeg (const TCT *tct, SEG *seg, Array hits,AC_HANDLE h0)
{
  int dummy, jj, ii, iMax = arrayMax (hits) ;
  LHIT *up, *vp ;
  Array dnaR ;
  int ln, Ln = arrayMax (tct->dna) ;

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
}  /* tctAddReverseSeg */

/*************************************************************************************/

static void  tctComplementSeg (const TCT *tct, SEG *seg, Array hits,AC_HANDLE h0)
{
  int Ln, ln, dummy, ii, iMax = arrayMax (hits) ;
  LHIT *up ;
  CLONE *clone ;
  LANE *lane ;
  Array dnaR ;

   for (ii = 0 ; ii < iMax ; ii++)
    {
      up = arrp (hits, ii, LHIT) ;

      if (! tct->SAM)
	{
	  lane = arrp (tct->lanes, up->lane, LANE) ;
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
      if ( tct->SAM)
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
      Ln = arrayMax (tct->dna) ;
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
}  /* tctComplementSeg */

/*************************************************************************************/

static void  tctAlignSeg (const TCT *tct, SEG *seg, Array hits, BOOL isDown, AC_HANDLE h0)
{
  int ii, iMax = arrayMax (hits) ;
  LHIT *up ;
  Array dna = isDown ? tct->dna : tct->dnaR ;

   for (ii = 0 ; ii < iMax ; ii++)
    {
      up = arrp (hits, ii, LHIT) ;
      if (isDown != up->isDown)
	continue ;
      up->err = arrayHandleCreate (12, A_ERR, h0) ;
      if (! up->dna)
	continue ;

      aceDnaTrackErrors (up->dna, up->x1, &(up->x2)
			 , dna, up->a1, &(up->a2)
			 , 0, up->err
			 , 3, 8, FALSE, 0
			 , TRUE) ; /* err was just created and is already set t zero */
      /* tct->nErrTrackings++ ; not thread safe */
    }
  return ;
} /* tctAlignSeg */

/*************************************************************************************/
 
static void tctSetSubsegs (const TCT *tct, SEG *seg)
{
  AC_HANDLE h = ac_new_handle () ;
  Array subsegs = 0 ;
  int b, bMax = arrayMax (seg->bridges) ;
  BRIDGE *bb0 = bMax ? arrp (seg->bridges, 0, BRIDGE) : 0 ;
  BOOL debug = FALSE ;

  if (debug) tctShowBridges(seg->bridges) ;

  for (b = 1 ; b < bMax ; b++)
    {
      int aa ;
      SEG *subseg = 0 ;
      BRIDGE *bb = arrp (seg->bridges, b, BRIDGE) ;
      if (bb->type <= 0)
	continue ;
      if (bb0->ln == bb->ddx + bb->ddy) /* perfect insertion */
	continue ;
      if (bb0->ln < bb->ddx + bb->ddy + 2 * tct->Delta) /* short imperfection */
	continue ;
      if (bb->dx > bb0->dx - 5 && bb->dx < bb0->dx + 5)
	{ /* split midway */
	  if (! subsegs)
	    subsegs = seg->subsegs = arrayHandleCreate (4, SEG, 0) ;
	  subseg = arrayp (subsegs, arrayMax (subsegs), SEG) ;
	  subseg->isSubseg = TRUE ;
	  subseg->a1 = seg->a1 ;
	  subseg->b1 = seg->b1 ;
	  subseg->c1 = seg->c1 ;
	  aa = seg->a1 + bb->ddx + (bb0->ln - bb->ddx - bb->ddy)/2 ;
	  subseg->a2 = aa ;  
	  subseg->b2 = aa ;
	  subseg->c2 = aa ;

	  subseg = arrayp (subsegs, arrayMax (subsegs), SEG) ;
	  subseg->isSubseg = TRUE ;
	  subseg->a1 = aa ;
	  subseg->b1 = aa ;
	  subseg->c1 = aa ;
	  subseg->a2 = seg->a2 ;
	  subseg->b2 = seg->b2 ;
	  subseg->c2 = seg->c2 ;
	}
    }
  if (subsegs && arrayMax (subsegs))
    {
      int j, jMax = arrayMax (subsegs) ;
      SEG *subseg = arrp (subsegs, 0, SEG) ;
      
      for (j = 0 ; j < jMax ; subseg++, j++)
	tctAnalyzeSeg (tct, subseg) ;
    }

  ac_free (h) ;
  return ;
} /*  tctSetSubsegs */

/*************************************************************************************/
 
static void tctAnalyzeSeg (const TCT *tct, SEG *seg)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL debug = FALSE ;

  Array hits = tctSeg2Hits (tct, seg, h) ;
  int x = 0x3 ; /* both strands, 0x1 forward, 0x2 reverse*/

  if (seg->a2 - seg->a1 >= tct->minBridge)
    if (arrayMax (hits))
      {
	if (debug)
	  fprintf (stderr, "... analyzeSeg %d %d nhits=%d\n", seg->a1, seg->a2, arrayMax (hits)) ;

	if (x & 0x2) tctComplementSeg (tct, seg, hits, h) ;
	if (x & 0x2) tctAddReverseSeg (tct, seg, hits, h) ;
	if (x & 0x1) tctAlignSeg (tct, seg, hits, TRUE, h) ; /* fills seg->err isDown */
	if (x & 0x2) tctAlignSeg (tct, seg, hits, FALSE, h) ; /* fills seg->err NOT isDown*/
	if (x & 0x1) tctCombSegPosition (tct, seg, hits, TRUE, h) ; /* extends left-right */
	if (x & 0x2) tctCombSegPosition (tct, seg, hits, FALSE, h) ; /* extends right-left */
	if (x & 0x1) tctCombSeg (tct, seg, hits, TRUE, h) ; /* extends left to right */
	if (x & 0x2) tctCombSeg (tct, seg, hits, FALSE, h) ; /* extends right to left */
	
	tctMakeBridge (tct, seg, hits, tct->h) ; /* stitch */
	tctFuseBridges (tct, seg, hits) ; /* stitch */
	tctBridgeSupport (tct, seg, hits, h) ;
	if (tctBridgeSetUpper (tct, seg, hits, h) && 
	    ! seg->isSubseg &&
	    ! seg->isTooShort && 
	    seg->a2 > seg->a1 + 4 * tct->Delta
	    )
	  tctSetSubsegs (tct, seg) ;
      }
  ac_free (h) ;

  return ;
} /* tctAnalyzeSeg */

/*************************************************************************************/

static void tctAnalyzer (const void *vp)
{
  const TCT *tct = (const TCT *)vp ;
  int ii ;
  BOOL debug = FALSE ;

  while (channelGet (tct->analyzeChan, &ii, int))
    {
      SEG *seg = arrayp (tct->segs, ii, SEG) ;
      if (debug)
	fprintf (stderr, "=== tctAnalyzeSeg %d b1=%d a1=%d a2=%d b2=%d %s\n",ii, seg->b1, seg->a1, seg->a2, seg->b2, timeShowNow()) ;
      tctAnalyzeSeg (tct, seg) ;
    }
  ii = 100 ; /* success */
  channelPut (tct->doneChan, &ii, int) ;
  return ;
} /* tctAnalyzer */

/*************************************************************************************/

static int tctCompressSegs (const TCT *tct)
{
  Array segs = tct->segs ;
  int iMax = arrayMax (tct->segs) ;
  int delta = tct->Delta ;

  if (iMax > 1)
    {
      SEG *up, *vp, *wp ;
      int j, ii, jj ;
  
      arraySort (segs, tctSegOrder) ;
      for (ii = jj = 0, up = vp = arrayp (tct->segs, ii, SEG) ; ii < iMax ; up++, ii++)
	{
	  BOOL ok = TRUE ;

	  for (wp = vp - 1, j = jj - 1 ; ok && j >= 0 ; wp--, j--)
	    {
	      if (wp->a1 == vp->a1 && wp->a2 == up->a2)
		ok = FALSE ;
	      else if (up->a1 < wp->a1 + delta && up->a2 > wp->a2 - delta && up->a2 <= wp->a2)
		{ ok = FALSE ; }
	      else if (up->a1 < wp->a1 + delta && up->a2 > wp->a2 && up->a2 < wp->a2 + delta)
		{ wp->a2 = up->a2 ; wp->b2 = up->b2 ; ok = FALSE ; }
	      else if (wp->a1 <= vp->a1 && wp->a2 > up->a1 + .9 * (up->a2 - up->a1) && (up->b2 - up->a1 > 100 * delta) &&
		  (up->a2 - up->a1) + 30 > (wp->a2 - wp->a1) && (up->a2 - up->a1) < (wp->a2 - wp->a1) + 30 &&
		  up->a1 < wp->a1 + 30 && up->a2 < wp->a2 + 30 && wp->a2 < vp->a2 + 30
		  )
		ok = FALSE ;
	    }
	  if (ok)
	    {
	      if (ii > jj)
		*vp = *up ;
	      jj++ ; vp++ ;
	    }
	}
      arrayMax (segs) = iMax = jj ;
    }
  return iMax ;
} /*  tctCompressSegs */

/*************************************************************************************/

static void tctAnalyzeSegs (TCT *tct)
{
  int n = 0, nn = tct->nAna ; /* number of analysers */
  int ii ;
  int iMax = arrayMax (tct->segs) ;

  tct->analyzeChan = channelCreate (30, int, tct->h) ;
  for (ii = 0 ; ii < nn ; ii++)
    wego_go (tctAnalyzer, tct, TCT) ;  /* will channelPut (tct->doneChan) */ 

  for (ii = 0 ; ii < iMax ; ii++)
    channelPut (tct->analyzeChan, &ii, int) ;
  channelClose (tct->analyzeChan) ;

  while (nn--)
    if (! channelGet (tct->doneChan, &n, int) || n != 100)
    {
      fprintf (stderr, "tctAnalyze received a bad value %d\n", n) ;
      exit (1) ;
    }
  return ;
} /* tctAnalyze */

/*************************************************************************************/
/*************************************************************************************/

static void tctMergeSubSegs (const TCT *tct)
{
  Array segs = tct->segs ;
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
	  arraySort (bridges, tctBridgeOrder) ;
	  if (jj < ii)
	    {
	      SEG *segj = arrayp (segs, jj, SEG) ;
	      *segj = *seg ;
	    }
	  jj++ ;
	}
    }
  iMax = arrayMax (segs) = jj ;
  arraySort (segs, tctSegOrder) ;
} /* tctMergeSubSegs */

/*************************************************************************************/

static int tctRationalizeOne (const TCT *tct, SEG *seg, BRIDGE *bb, int *previousp)
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
  if (ds >= tct->minIntron  && ds <= tct->maxIntron  && dt == 1) /* possible intron */
    {
      char *cs1, *cp1, *cs11 = 0, *cp11 ;   /* look downstream for gt_ag OR gc_ag OR ct_ac OR ct_ac OR ct_gc OR gt_ag */
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
  if (ds < DELMAX && dt < INSMAX)
    {
      if (ds > 0 && ds < DELMAX)
	{ strncpy (buf1, cs, ds) ; buf1[ds] = 0 ; }
      if (dt > 0 && dt < INSMAX)
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
      while (ds > 1 && dt > 1 && --cp >= cs && --cq >= ct)
	if (ace_upper (*cp) == ace_upper(*cq))
	  { cp-- ; cq-- ; ds-- ; dt-- ; }
	else
	  break ;
      cp = cs0 ; 
      cq = ct0 ;
      while (ds > 1 && dt > 1 && *cs && *cp && ace_upper (*cp) == ace_upper(*cq))
	{ cp++ ; cq++ ; ds-- ; dt-- ; }
    }
  bb->ds = ds ;
  bb->dt = dt ;
  bb->ddx = ddx ;

  return ddx ;
} /*  tctRationalizeOne */

/*************************************************************************************/

static int tctRationalize (const TCT *tct)
{ 
  int ii ;
  Array segs = tct->segs ;
  int iMax = arrayMax (segs) ;
  int previous = 7 ;

  arraySort (segs, tctSegOrder) ;
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
		  tctRationalizeOne (tct, seg, bb, &previous) ;
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
		  tctAnalyzeSeg (tct, seg) ;
		}
	    }
	  else
	    break ;
	}
    }

  arraySort (segs, tctSegOrder) ;
  for (ii = 0 ; ii < iMax ; ii++)
    {
      int bMax ;
      SEG *seg = arrayp (segs, ii, SEG) ;
      
      if (seg->subsegs)
	continue ;
      
      bMax = seg->bridges ? arrayMax (seg->bridges) : 0 ;
      if  (bMax > 1)
	arraySort (seg->bridges, tctBridgeOrder) ;
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
		   
		   if (bb1->ds == bb2->ds &&
		       bb1->dt == bb2->dt &&
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
} /* tctRationalize */

/*************************************************************************************/

static void tctAliLn (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (tct->outFileName, ".aliLn", tct->gzo, h) ;
  KEYSET a, aliLn = tct->aliLn = arrayHandleCreate (10000, KEY, tct->h) ;
  int ii, jj, iiMax = arrayMax (tct->lanes) ;
  LANE *lane ;

  for (ii = 1 ; ii < iiMax ; ii++)
    {
      lane = arrp (tct->lanes, ii, LANE) ;
      a = lane->aliLn ;
      if (a)
	for (jj = 0 ; jj < keySetMax (a) ; jj++)
	  keySet (aliLn, jj) += keySet (a, jj) ;
    }
  for (jj = 0 ; jj < keySetMax (aliLn) ; jj++)
    aceOutf (ao, "%d\t%d\n", jj, keySet (aliLn, jj)) ;

  ac_free (h) ;
} /* tctAliLn */

 /*************************************************************************************/

static void tctExportGeneFusions (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT aoGF = 0 ;
  Array geneFusions = tct->geneFusions ;

  arraySort (geneFusions, tctGeneFusionOrder) ;
      
#ifdef JUNK
  /* we want to see the full list names of the supporting reads */
  if (arrayMax (geneFusions) > 1)
    {
      int i, j, iMax ;
      GF *up, *vp ;
      
      iMax = arrayMax (geneFusions) ;
      up = vp =  arrp (geneFusions, 0, GF) ;
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
		{ vp->tag += up->tag ; vp->seq += up->seq ; }
	      continue ;
	    }
	  if (up->g1 > 0)
	    {
	      vp++ ; j++ ;
	      if (j < i)
		*vp = *up ;
	    }
	}
      arrayMax (geneFusions) = j + 1 ;
    }
#endif
  if (geneFusions)
    {
      int i ;
      int iMax =  arrayMax (geneFusions) ;
      GF *gf ;
      const char *lane = tct->lane ? tct->lane : tct->run ;
      aoGF = aceOutCreate (tct->outFileName, ".geneFusion.txt", tct->gzo, h) ;
      aceOutDate (aoGF, "####", "Gene fusion candidates") ; 
      aceOutf (aoGF, "#Run\tRead\tFusion\tGene_A\tGene_B\tChrom_A\tChrom_B\tDistance\tType\tx1 A\tx2 A\tmRNA_A\tfrom\tto\tx1 B\tx2 B\tmRNA_B\tfrom\tto\tDistinct_supports\tSupport\tGene_A supports\tGene_B supports\tscore A\tscore B\tAli A\tAli B\tc1 A\tc2 A\tc1 B\tc2 B\tOverlapOrGap\td2donor\td2acceptor\toverhang1\toverhang2\n") ;

      for (i = 0, gf = arrayp (tct->geneFusions, 0, GF) ; i < iMax ; gf++, i++)
	if (gf->g1 > -1)
	  {
	    const char *cp1 = dictName (tct->geneDict, gf->g1) ;
	    const char *cp2 = dictName (tct->geneDict, gf->g2) ;
	    GENE *ga = tct->genes ? arrayp (tct->genes, gf->g1, GENE) : 0 ;
	    GENE *gb = tct->genes ? arrayp (tct->genes, gf->g2, GENE) : 0 ;
	    ga = tct->genes ? arrayp (tct->genes, gf->g1, GENE) : 0 ; /* redo in case vp->gene increased the array */

	    if (!strncmp (cp1, "X__", 3)) cp1 += 3 ;
	    if (!strncmp (cp2, "X__", 3)) cp2 += 3 ;
	    
	    aceOutf (aoGF, "%s\t%s\t%s__%s%s\t%s\t%s"
		     , lane
		     , dictName (gf->lane->cloneDict, gf->clone) 
		     , cp1, cp2, gf->wo
		     , cp1, cp2
		     ) ;
	    aceOutf (aoGF, "\t%s\t%s\t%d\t%s"
		     , ga && ga->chrom ? dictName (tct->geneDict, ga->chrom) : ""
		     , gb && gb->chrom ? dictName (tct->geneDict, gb->chrom) : ""
		     , gf->chromDistance
		     , gf->isPair ? "PAIR" : "READ"
		     ) ;
	    aceOutf (aoGF, "\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d"
		     , gf->x1, gf->x2, dictName (tct->geneDict, gf->t1), gf->a1, gf->a2
		     , gf->y1, gf->y2, dictName (tct->geneDict, gf->t2), gf->b1, gf->b2
		     ) ;
	    aceOutf (aoGF, "\t%d\t%d\t%d\t%d"
		     , gf->seq, gf->tag
		     , keySet (tct->geneSupport, gf->g1)
		     , keySet (tct->geneSupport, gf->g2)
		     ) ;
	    aceOutf (aoGF, "\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n"
		     , gf->score1, gf->score2
		     , gf->c2 - gf->c1 + 1, gf->d2 - gf->d1 + 1
		     , gf->c1, gf->c2, gf->d1, gf->d2
		     , gf->gap, gf->d2donor, gf->d2acceptor
		     , gf->overhang1 ? dictName (tct->overhangDict, gf->overhang1) : "-"
		     , gf->overhang2 ? dictName (tct->overhangDict, gf->overhang2) : "-"
		     ) ;
	  }
    }
  if (tct->outFileName && tct->fusionHisto)
    {
      int i ;
      const char *lane = tct->lane ? tct->lane : tct->run ;
      aoGF = aceOutCreate (tct->outFileName, ".geneFusionHisto.txt", FALSE, h) ;
      aceOutDate (aoGF, "###", "Gene fusion overlap histo") ;
      aceOutf (aoGF, "#Overlap\tNumber of reads\n") ;

      for (i = -100 ; i <= 100 ; i++)
	aceOutf (aoGF, "%s\t%s\t%05d\t%d\n", tct->run, lane, i, keySet (tct->fusionHisto, i + 100)) ;
    }

  ac_free (h) ;
} /* tctExportGeneFusions */

/*************************************************************************************/
#ifdef JUNK
 
static int tctQcOrder (const void *va, const void *vb)
{
  const QC *up = (const QC *)va, *vp = (const QC *)vb ;
  int n ;
  
  n = up->type - vp->type ; if (n)  return n ;

  return 0 ;
} /* qcOrder */


/*************************************************************************************/

static void tctQcExport (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  Array lanes = tct->lanes ;
  int jj ;

  ACEOUT ao = aceOutCreate (tct->outFileName, ".qc.tsf", tct->gzo, h) ;
  aceOutDate (ao, "###", "Magic-BLAST QC") ;
  aceOutf (ao, "# Run %s %s\n", tct->run, "This .tsf file has repeated entries per data block, please use tsf --merge, then parse it in acedb MetaDB") ;

  for (jj = 1 ; jj < arrayMax (lanes) ; jj++)
    {
      LANE *vp = arrp (lanes, jj, LANE) ;
      DICT *qcDict = vp->qcDict ;
      QC *qc ;
      int ii, iiMax = arrayMax (vp->qc) ;

      for (ii = 0, qc = arrp (vp->qc, 0, QC) ; ii < iiMax ; ii++, qc++)
	if (qc->type)
	  {
	    int ki, kf, kt ;
	    const char *fmt =  dictName (qcDict, qc->fmt) ;
	    arraySort (vp->qc, tctQcOrder) ;
	    aceOutf (ao, "%s\t%s\t%s"
		     , dictName (qcDict, qc->type)
		     , tct->run
		     , fmt
		     ) ;
	    for (ki = kf = kt = 0 ; *fmt ; fmt++)
	      switch ((int)*fmt)
		{
		case 'i':
		  aceOutf (ao, "\t%ld", qc->x[ki++]) ;
		  break ;
		case 'f':
		  aceOutf (ao, "\t%.3f", qc->z[kf++]) ;
		  break ;
		case 't':
		  aceOutf (ao, "\t%s", dictName (qcDict, qc->x[kt++])) ;
		  break ;
		}
	    aceOutf (ao, "\n") ;
	  }
    }
  ac_free (h) ;
  return ;
} /* tctQcExport */
#endif

/*************************************************************************************/

static void tctDeUnoFuseOneLane (const TCT *tct, int myLane)
{
  UNO *uno ;
  Array lanes = tct->lanes ;
  BigArray deUno = tct->deUno ;
  DICT *deUnoDict = tct->deUnoDict ;
  long int ii, iiMax = 0 ;
  Array wig = tct->wig ; 
  LANE *lane = arrp (lanes, myLane, LANE) ;
  
  /* fuse the wiggles */
  if (lane->wig)
    {
      
      int iMax = arrayMax (lane->wig) ;
      if (iMax)
	{
	  int i, *ip, *jp ;
	  array (wig, iMax-1, int) += 0 ;
	  for (i = 0, ip = arrayp (lane->wig, 0, int), jp =  arrayp (wig, 0, int) ; i < iMax ; i++)
	    jp[i] += ip[i] ;
	}
    }


  /* fuse the deUno */
  if (lane->deUno)
    {
      UNO *vno ;
      int j = 0, jMax = lane->deUno ? arrayMax (lane->deUno) : 0 ;
      ii = bigArrayMax (deUno) ;
      if (jMax)
	{
	  for (j = 0,  vno = arrp (lane->deUno, 0, UNO) ; j < jMax ; j++, vno++)
	    {
	      uno = bigArrayp (deUno, ii++, UNO) ;
	      *uno = *vno ;
	      if (vno->type)
		dictAdd (deUnoDict, dictName (lane->deUnoDict, vno->type), &(uno->type)) ;
	      if (vno->tag)
		dictAdd (deUnoDict, dictName (lane->deUnoDict, vno->tag), &(uno->tag)) ;
	      if (vno->insert)
		dictAdd (deUnoDict, dictName (lane->deUnoDict, vno->insert), &(uno->insert)) ;
	    }
	  bigArraySort (deUno, tctUnoOrder) ;
	  /* keep happy few */
	  ii = bigArrayMax (deUno) ;
	  if (ii)
	    {
	      long int i, j ;
	      for (i = j = 0,  uno = vno = bigArrp (deUno, 0, UNO) ; i < ii ; i++, uno++)
		{
		  if (uno > vno && 
		      vno->a1 == uno->a1 && vno->a2 == uno->a2 &&
		      vno->target == uno->target &&
		      vno->type == uno->type &&
		      vno->tag == uno->tag &&
		      vno->insert == uno->insert
		      )
		    { vno->nt += uno->nt ; vno->ns += uno->ns ; }
		  else if (vno < uno)
		    { j++, vno++ ; if(vno < uno) *vno = *uno ; } 
		}
	      ii = iiMax = bigArrayMax (deUno) = j + 1 ;
	    }
	}
    }
}  /* tctDeUnoFuseOneLane */

/*************************************************************************************/

static void tctDeUnoExport (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  UNO *uno ;
  BigArray deUno = tct->deUno ;
  DICT *geneDict = tct->geneDict ;
  DICT *deUnoDict = tct->deUnoDict ;
  const char *method =  tct->SAM ? "MagicBlast" : "deUno" ;
  long int ii, iiMax = bigArrayMax (deUno) ;

  if (iiMax)
    {
      Array wig = tct->wig ;
      ACEOUT ao = aceOutCreate (tct->outFileName, ".deUno.tsf", tct->gzo, h) ;
      aceOutDate (ao, "###", "De uno discovery of substitutions, insertions, deletions, introns, gene-fusions") ;
      aceOutf (ao, "# %s", " Target\tRun\tFormat\tLast matching base before\tFirst matching base after\tLength\tNumber of supports\tNumber of distinct supports\tWiggle before\tWiggle after\tMethod\tTag\tLength\tInsert\n") ;

      /* integrate the wiggle, so far, it represents the differential of the wiggle */
      if (wig)
	{
	  int n = 0, i, *ip, iMax = arrayMax (wig) ;
	  for (i = 0, ip = arrayp (wig, 0, int) ; i < iMax ; i++)
	    { n += ip[i] ; ip[i] = n ; }
	}
      if (0) fprintf (stderr, "#2 wig[40]=%d\n", array (wig, 40, int)) ;

      for (ii = 0, uno = bigArrp (deUno, 0, UNO) ; ii < iiMax ; ii++, uno++)
	if (uno->isExternal || (uno->ns && uno->type))
	  {
	    int a1 = uno->a1, a2 = uno->a2 ;
	    int nt = uno->nt ;
	    int ns = uno->ns ;
	    int da = uno->da ;
	    int nd, na, nad ;
	    int minSnpCover = tct->minSnpCover ;
	    int minSnpCount = tct->minSnpCount ;
	    float minSnpFrequency = tct->minSnpFrequency ;
	    int targetStart = keySet (tct->targetStart , uno->target) ; 

	    if (da == 0) /* sub */
	      {
		nd = array (wig, targetStart + a1 + 1, int) ;
		na = array (wig, targetStart + a2 - 1, int) ;
	      }
	    else  /* ins or del */
	      { 
		nd = array (wig, targetStart + a1, int) ;
		na = array (wig, targetStart + a2, int) ;
	      }
	    if (nd < nt) nd = nt ;
	    if (na < nt) na = nt ;
	    nad = na < nd ? na : nd ;
	    if (uno->isExternal ||
		(
		 nt >= minSnpCount &&
		 nad >= minSnpCover &&
		 nt * 100 >= minSnpFrequency * nad 
		 )
		)	    
	      {
		const char *insert = 0 ;
		const char *ccp = dictName (deUnoDict, uno->type) ;
		char *cp, *cq, bufType[1024] ;

		if (uno->insert)
		  insert = dictName (deUnoDict, uno->insert) ;
		strncpy (bufType, ccp, 1022) ;
		bufType[1023] = 0 ;

		cp = strstr (bufType, "Del_") ;
		if (cp)
		  {
		    cp = strstr (cp, "::") ;
		    if (cp && insert && strlen (insert) + strlen (bufType) < 1010)
		      {
			int k = strlen(insert) ; 
			cp[1] = cp[2] ; 
			if (k < 16)
			  strcpy (cp+2, insert) ;
			else
			  {
			    strncpy (cp+2, insert,8) ;
			    cp[2+8] = '_' ;
			    strncpy (cp+2+9, insert + k - 8 , 8) ;
			  }
			cq = bufType + strlen (bufType) ; 
			cq[0] = ':' ;
			cq[1] = cp[1] ;
			cq[2] = 0 ;
		      }
		  }
		cp = strstr (bufType, "Ins_") ;
		if (cp)
		  {
		    cp = strstr (cp, ":") ;
		    if (cp && insert && strlen (insert) + strlen (bufType) < 1010)
		      {
			int k = strlen(insert) ; 
			cp[3] = cp[1] ; 
			if (k < 16)
			  strcpy (cp+4, insert) ;
			else
			  {
			    strncpy (cp+4, insert,8) ;
			    cp[4+8] = '_' ;
			    strncpy (cp+4+9, insert + k - 8 , 8) ;
			  }
		      }
		  }
		aceOutf (ao, "%s:%s\t%s\ttttiiiitttt\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s%s%s\n"
			 , dictName (geneDict, uno->target)
			 , bufType
			 , tct->run
			 , uno->a1, uno->a2, uno->ddx
			 , nt, ns
			 , nd
			 , na 
			 , method
			 , uno->tag ? dictName (deUnoDict, uno->tag) : "X"
			 , insert ? "\t" : ""
			 , insert ? insert : ""
			 ) ;
	      }
	  }
      if (0) fprintf (stderr, "#3 wig[40]=%d\n", array (wig, 40, int)) ;
    }
  ac_free (h) ;
  return ;
} /* tctDeUnoExport */

 /*************************************************************************************/

static void tctExport (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  DICT *targetDict = tct->targetDict ;
  int iSeg ;
  int t1 = tct->t1 ? tct->t1 : tct->offset ;
  ACEOUT ao = aceOutCreate (tct->outFileName, ".bridges.ace", tct->gzo, h) ;
  const char *method =  tct->method ? tct->method : "variant_caller" ;
  const char *run = tct->run ? tct->run : "variant_caller" ;
  Array segs = tct->segs ;
  int iMax = arrayMax (segs) ;
  vTXT vNam = vtxtHandleCreate (h) ;

  arraySort (segs, tctSegOrder) ;
  for (iSeg = 0 ; iSeg < iMax ; iSeg++)
    {
      int b, bMax ;
      SEG *seg = arrayp (segs, iSeg, SEG) ;

      if (seg->subsegs)
	continue ;
      bMax = seg->bridges ? arrayMax (seg->bridges) : 0 ;
      if  (bMax > 1)
	arraySort (seg->bridges, tctBridgeOrder) ;

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
	  BOOL ok = FALSE ;

	  if (bb->isExternal)
	    {
	      int b1 ;
	      for (b1 = 1 ; b + b1 < bMax ; b1++)
		if (bb[b1].type >= 0 &&
		    bb[b1].coverp + bb[b1].coverm >=  tct->minSnpCover
		    )
		  { bb->isExternal = FALSE ; break ; }
	    }

	  if (bb->type <= 0) continue ;
	  if (cover + covernu >= tct->minSnpCover && 100.0 * (m + mnu) >= (cover + covernu) * tct->minSnpFrequency) ok = TRUE ;
	  if (cover >= tct->minSnpCover && 100.0 * (m) >= (cover) * tct->minSnpFrequency) ok = TRUE ;
	  if (covernu >= tct->minSnpCover && 100.0 * (mnu) >= (covernu) * tct->minSnpFrequency) ok = TRUE ;
	  if (!ok) continue ;

	  bb->c1 =  seg->b1 + t1 + bb->ddx ;

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

	  tct->snpExported++ ;

	  if (ds>30) cs = "" ;
	  if (dt>30) ct = "" ;

	  vtxtClear (vNam) ;
	  switch (bb->type) 
	    {
	    case 1:
	      if (ds == 1)
		vtxtPrintf (vNam, "Sub:%s:%s", cs, ct) ;
	      else if (ds == dt)
		vtxtPrintf (vNam, "Sub_%d:%s:%s", ds, cs, ct) ;
	      else 
		vtxtPrintf (vNam, "DelIns_%d_%d:%s:%s", ds, dt, cs, ct) ;
	      break ;
	    case 2: 
	      if (ds == 1 && dt == 2)
		vtxtPrintf (vNam, "Dup:%s:%s", cs, ct) ;
	      else if (ds == 1 && dt > 2)
		vtxtPrintf (vNam, "Dup_%d:%s:%s", dt - ds, cs, ct) ;
	      else 
		vtxtPrintf (vNam, "DelIns_%d_%d:%s:%s", ds, dt, cs, ct) ;
	      break ;
	    case 3: 
	      if (ds == 1 && dt == 2)
		vtxtPrintf (vNam, "Ins:%s:%s", cs,ct) ;
	      else if (ds == 1 && dt > 2)
		vtxtPrintf (vNam, "Ins_%d:%s:%s", dt - ds, cs,ct) ;
	      else 
		vtxtPrintf (vNam, "DelIns_%d_%d:%s:%s", ds, dt, cs, ct) ;
	      break ;
	    case 4: 
	      if (ds == 2 && dt == 1)
		vtxtPrintf (vNam, "Dim:%s:%s", cs,ct) ;
	      else if (ds > 2 && dt == 1)
		vtxtPrintf (vNam, "Dim_%d:%s:%s", ds - dt, cs,ct) ;
	      else 
		vtxtPrintf (vNam, "DelIns_%d_%d:%s:%s", ds, dt, cs, ct) ;
	      break ;
	    case 5: 
	      if (ds == 2 && dt == 1)
		vtxtPrintf (vNam, "Del:%s:%s", cs,ct) ;
	      else if (ds > 2 && dt == 1)
		vtxtPrintf (vNam, "Del_%d:%s:%s", ds - dt, cs,ct) ;
	      else 
		vtxtPrintf (vNam, "DelIns_%d_%d:%s:%s", ds, dt, cs, ct) ;
	      break ;
	    case 6:
	      vtxtPrintf (vNam, "DelIns_%d_%d:%s:%s", ds, dt, cs, ct) ;
	      break ;
	    default:  /* 7-12 : intron */
	      da = bb->a2 > bb->a1 ? bb->a2 - bb->a1 + 1 : bb->a1 - bb->a2 + 1 ;
	      vtxtPrintf (vNam, "%s_%d:%d:%d"
			  , types[bb->type]
			  , da
			  , seg->b1 + bb->a1 + t1
			  , seg->b1 + bb->a2 + t1 
			  ) ;
	      break ;
	    }
	  
	  /* G1103T(ORF2a:L3606F) */
	  aceOutf (ao, "Variant \"%s:%d:%s\" // seg->a1/a2=%d/%d, b1/b2=%d/%d iSeg=%d\n"
		   , dictName (targetDict, seg->target)
		   , seg->b1 + t1 + ddx 
		   , vtxtPtr (vNam)
		   , seg->a1
		   , seg->a2
		   , seg->b1
		   , seg->b2
		   , iSeg
		   ) ;

	  cs = bb->buf1 ; ct = bb->buf2 ;
	  aceOutf (ao, "VCF %d %s %s\n"
		   , seg->b1 + t1 + ddx
		   , cs
		   , ct
		   )  ;

	  if (bb->type >=7 && bb->type <= 12) /* intron */
	    {
	      da = bb->a2 > bb->a1 ? bb->a2 - bb->a1 + 1 : bb->a1 - bb->a2 + 1 ;
	      aceOutf (ao, "%s %d %d %d\n"
		       , types[bb->type]
		       , seg->b1 + bb->a1 + t1
		       , seg->b1 + bb->a2 + t1 
		       , da
		       ) ;
	    }
	  else if (ds == 1 && dt == 1)
	    aceOutf (ao, "Substitution\n%s2%s\n", cs, ct) ;
	  else if (da == 0)
	    aceOutf (ao, "Multi_substitution %d %s %s\n"
		     , ds
		     , ds < 300 ? cs : ""
		     , dt < 300 ? ct : ""
		     ) ;
	  else if (ds == 1 && dt == 2)
	    aceOutf (ao, "%s%c\n", bb->type == 2 ? "Dup" : "Ins", ct[1]) ;
	  else if (dt == 1 && ds == 2)
	    aceOutf (ao, "%s%c\n", bb->type == 4 ? "Dim" : "Del", cs[1]) ;
	  else if (ds <= 1 && dt > 1)
	    aceOutf (ao, "Multi_insertion\t%d\t%s\n", da, dt < 120 ? ct + 1 : "") ; 
	  else if (dt <= 1 && ds > 1)
	    aceOutf (ao, "Multi_deletio\tn%d\t%s\n", -da, ds < 120 ? cs + 1 : "") ; 
	  else
	    aceOutf (ao, "DelIns %d %d %s %s\n"
		     , ds, dt
		     , ds + dt < 120 ? cs : ""
		     , ds + dt < 120 ? ct : ""
		     ) ;

	  if (0)
	    {  /* now established by tsnp.c -db_translate */
	      *cp = cc1 ; *cq = cc2 ;
	      if (arrayMax (bb->targetSnippet) && arrayMax (bb->targetSnippet) < 120)
		aceOutf (ao, "Reference_genomic_sequence %s\n", arrp (bb->targetSnippet, 0, char) );
	      if (arrayMax (bb->variantSnippet) && arrayMax (bb->variantSnippet) < 120)
		aceOutf (ao, "Observed__genomic_sequence %s\n", arrp (bb->variantSnippet, 0, char)) ;
	    }
	  aceOutf (ao, "Method %s\n", method) ;
	  cp = "" ;
	  if (dt == 1 && ds ==1)
	    cp = hprintf (h, "\"Base %d:%s>%s is modified\"", seg->b1 + t1 + bb->ddx, cs, ct) ;
	  else if (dt == 1 && ds == 2)
	    cp = hprintf (h, "\"Base %d:%c is deleted\"", seg->b1 + t1 + bb->ddx + 1, cc1) ;
	  else if (dt == 1 && ds > 2)
	    cp = hprintf (h, "\"Bases %d to %d (%d bases) are deleted\"", seg->b1 + t1 + bb->ddx + 1, seg->b1 + t1 + bb->ddx + ds - 1, ds -1) ;
	  else if (dt == 2 && ds == 1)
	    cp = hprintf (h, "\"Base %c is inserted between base %d and %d\"", cc2, seg->b1 + t1 + bb->ddx, seg->b1 + t1 + bb->ddx + 1) ;
	  else if (dt > 2 && ds == 1)
	    cp = hprintf (h, "\"%d bases are inserted between base %d and %d\"", dt - 1, seg->b1 + t1 + bb->ddx, seg->b1 + t1 + bb->ddx + 1) ;
	  else if (dt < 1 && ds > 1)
	    cp = hprintf (h, "\"Bases %d to %d %s replaced by %d bases\"", seg->b1 + t1 + bb->ddx + 1, seg->b1 + t1 + bb->ddx + ds - 1, ds > 2 ? "are" : "is", dt - 1) ;

	  aceOutf (ao, "IntMap \"%s\" %d %d %s \n"
		   , dictName (targetDict, seg->target)
		   , seg->b1 + t1 + bb->ddx 
		   , seg->b1 + t1 + bb->ddx + ds
		   , cp
		   ) ;
	  aceOutf (ao, "Found_in_genome\n") ;	  
	  aceOutf (ao, "Parent_sequence \"%s\"\n", dictName (targetDict, seg->target)) ;
	  if (coverp) 
	    aceOutf (ao, "fCounts %s %d %d %d Frequency %.2f\n", run, mp, wp, coverp, fp) ; 
	  if (coverm)
	    aceOutf (ao, "rCounts %s %d %d  %d Frequency %.2f\n", run, mm, wm, coverm, fm) ;
	  if (cover)
	    aceOutf (ao, "nsCounts %s %d %d  %d Frequency %.2f\n", run, m, w, cover, ff) ;
	  if (coverpnu) 
	    aceOutf (ao, "fCounts %sNU %d %d %d Frequency %.2f\n", run, mpnu, wpnu, coverpnu, fpnu) ; 
	  if (covermnu)
	    aceOutf (ao, "rCounts %sNU %d %d  %d Frequency %.2f\n", run, mmnu, wmnu, covermnu, fmnu) ;
	  if (covernu)
	    aceOutf (ao, "nuCounts %s %d %d  %d Frequency %.2f\n", run, mnu, wnu, covernu, ffnu) ;

	  aceOutf (ao, "\n") ;
	}
    }
  ac_free (ao) ;
  ac_free (h) ;
} /*  tctExport */

/*************************************************************************************/
/*************************************************************************************/

static void tctGetTargetFasta (const void *vp)
{
  AC_HANDLE h = ac_new_handle () ;
  const TCT *tct = (const TCT *)vp ;
  const char *fNam = tct->targetFileName ;
  ACEIN ai = 0 ;
  int state = 0 ;
  int tt1 = tct->t1 ;
  int n ;
  int a1, a2, x, dx ;
  int t1 = tct->t1 ;
  int t2 = tct->t2 ;
  const char *target = tct->target ;
  Array dna = tct->dna ;
  Array dnaR = tct->dnaR ;
  int trgt = 0 ;

  if (! fNam)
    goto done ;
  ai = aceInCreate (fNam, 0, h) ;
  if (! ai)
    messcrash ("Cannot open target fasta file -targetFileName %s\n", fNam) ;
 
  if (tct->runTest) 
    {
      t1 -= tct->offset ;
      if (t2)
	t2 -= tct->offset ;
    }
  if (t1 < 1) 
    { tt1 += (1-t1) ; t1 = 1 ; }
  a1 = a2 = -1 ; x = 0 ;
  t1-- ; if (t2) t2-- ; tt1-- ; /* Plato */
  while (state < 2 && aceInCard (ai))
    {
      char *cp = aceInWord (ai), *cr ;
      if (! cp || ! *cp)
	continue ;
      switch (state)
	{
	case 0:
	case 1:
	  if (*cp == '>') /* start of next target */	
	    {
	      if (trgt && tct->t1)
		messcrash ("Parameter -t1 is imcompatible with -target_fasta containing more than 1 sequence") ;
	      keySet (tct->targetStop, trgt) = x ;
	      if (tct->runTest || (!target) || ! strcmp (target, cp + 1))
		state = 1 ; 
	      cr = strchr (cp, '|') ; if (cr) *cr = 0 ;
	      dictAdd (tct->geneDict, cp + 1, &trgt) ;
	      keySet (tct->targetStart, trgt) = x ;
	      break ;
	    }
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
    { /* we must do it in this way because the tct struct is copied by value
       * so dnaR is preallocated 
       */
      array (dnaR, n, char) = 0 ;
      memcpy (arrp (dnaR, 0, char), arrp (dna, 0, char), n) ;
      arrayMax (dna) = arrayMax (dnaR) = n - 1 ;
      reverseComplement (dnaR) ;
    }

  if (tct->makeTest && tct->t1)
    {
      ACEOUT ao = aceOutCreate (tct->makeTest, ".t1", FALSE, h) ;
      if (! ao)
	messcrash ("Cannot create test file %s.t1\n", tct->makeTest) ;
      aceOutf (ao, "##OFFSET %d\t%d\n", tct->t1, tct->t2) ;
      ac_free (ao) ;
    }

 if (tct->makeTest && dna && arrayMax (dna))
    {
      int i, iMax = arrayMax (dna) ;
      ACEOUT ao = aceOutCreate (tct->makeTest, ".target.fasta", FALSE, h) ;
      Array dna2 = dnaCopy (dna) ;
      dnaDecodeArray (dna2) ;

      if (! ao)
	messcrash ("Cannot create test file %s.target.fasta\n", tct->makeTest) ;
      aceOutf (ao, ">%s\n", tct->target) ;
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
  wego_log (hprintf (tct->h, "tctGetTargetFasta parsed %d bases in %s\n"
		     , arrayMax (dna), fNam)) ;
  
  if (! arrayMax (dna))
    {
      fprintf (stderr, "FATAL ERROR: no target sequence %s in file %s\n"
	       , tct->target
	       , fNam
	       ) ;
      ac_free (ai) ;
      exit (1) ;
    }
 done:
  ac_free (ai) ;
  ac_free (h) ;
  n = 0 ; /* success */
  channelPut (tct->doneChan, &(tt1), int) ;

  return ;
} /* tctGetTargetFasta */

/*************************************************************************************/
/*************************************************************************************/
typedef struct mrnaStruct {int gene, mrna, chrom, a1, a2, x1, x2, gene2 ; BOOL isRename ;} MRNA ;
/* we expect to read the file 
 * tmp/METADATA/mrnaRemap.gz
 * and obtain the coordinates of the mRNAs and their introns
 */

/*************************************************************************************/
 
static int tctMrnaOrder (const void *va, const void *vb)
{
  const MRNA *up = (const MRNA*)va, *vp = (const MRNA*)vb ;
  int n ;
  
  n = up->mrna - vp->mrna ; if (n) return n ;
  n = up->gene - vp->gene ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  n = up->x1 - vp->x1 ; if (n) return n ;
  n = up->x2 - vp->x2 ; if (n) return n ;

  return 0 ;
} /* tctMrnaOrder */

/*************************************************************************************/

static void tctGetMrnaRemap (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (tct->mrnaRemapFileName, 0, h) ;
  const char *tc = tct->target_class ;
  DICT *dict = tct->geneDict ;
  int nn = 0 ;
  Array aa  ;

  if (! dict)
    dict = tct->geneDict = dictHandleCreate (10000, tct->h) ;
  if (! tc)
    {
      fprintf (stderr, "// FATAL ERROR: argument -mrnaRemap requires aregument -target_class\n") ;
      exit (1) ;
    }

  if (! ai)
    {
      fprintf (stderr, "// FATAL ERROR : cannot open mrnaREmap file %s\n// Try TCT -h\n", tct->mrnaRemapFileName) ;
      exit (1) ;
    }
  aceInSpecial (ai, "\n") ;

  if (! tct->mrnas) tct->mrnas = arrayHandleCreate (100000, MRNA, tct->h) ;
  aa = tct->mrnas ;

  while (aceInCard (ai))
    {
      char *cp ;
      int a1 = 0, a2 = 0, x1, x2, gene = 0, chrom = 0, mrna = 0 ;
      MRNA *up ;

      cp = aceInWord (ai) ;
      if (! cp || strcmp (cp, tc))
	continue ;
      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;
      if (! cp)
	continue ;
      dictAdd (dict, cp, &mrna) ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &x1))
	continue ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &x2))
	continue ;
      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;
      if (! cp)
	continue ;
      dictAdd (dict, cp, &chrom) ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &a1))
	continue ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &a2))
	continue ;
      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;
      if (! cp)
	continue ;
      dictAdd (dict, cp, &gene) ;
      up = arrayp (aa, nn++, MRNA) ;
      up->mrna = mrna ;
      up->gene = gene ;
      up->chrom = chrom ;
      up->a1 = a1 ;
      up->a2 = a2 ;
      up->x1 = x1 ;
      up->x2 = x2 ;
      up->isRename = FALSE ;
    }
  ac_free (h) ;
  arraySort (aa, tctMrnaOrder) ;
  if (nn)
    {
      int i, m = 0, g = 0 ;
      MRNA *up ;
      KEYSET ks = tct->mrna2aa ;
      if (! ks)
	ks = tct->mrna2aa = keySetHandleCreate (tct->h) ;

      for (i = 0, up = arrp (aa, 0, MRNA) ; i < nn ; i++, up++)
	{
	  if (up->mrna != m)
	    { m = up->mrna ; keySet (ks, m) = i + 1 ; }
	  if (up->gene != g)
	    { g = up->gene ; keySet (ks, g) = i + 1 ; }
	}
    }
  return ;
} /* tctGetMrnaRemap */

/*************************************************************************************/

static void tctGetSplitMrnas (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (tct->splitMrnaFileName, 0, h) ;
  const char *tc = tct->target_class ;
  DICT *dict = tct->geneDict ;
  int nn = 0 ;
  Array aa  ;

  if (! dict)
    dict = tct->geneDict = dictHandleCreate (10000, tct->h) ;
  if (! tc)
    {
      fprintf (stderr, "// FATAL ERROR: argument -mrnaRemap requires aregument -target_class\n") ;
      exit (1) ;
    }

  if (! ai)
    {
      fprintf (stderr, "// FATAL ERROR : cannot open mrnaREmap file %s\n// Try TCT -h\n", tct->mrnaRemapFileName) ;
      exit (1) ;
    }
  aceInSpecial (ai, "\n") ;

  if (! tct->splitMrnas) tct->splitMrnas = arrayHandleCreate (100000, MRNA, tct->h) ;
  aa = tct->splitMrnas ;

  while (aceInCard (ai))
    {
      char *cp ;
      int x1, x2, mrna, gene2, gene = 0 ;
      MRNA *up ;
      GENE *ga, *gb ;

      cp = aceInWord (ai) ;
      if (! cp)
	continue ;
      dictAdd (dict, cp, &gene) ; 

      aceInStep(ai,'\t') ;
      cp = aceInWord (ai) ;
      if (! cp)
	continue ;
      mrna = gene ;
      if (cp[0] != '*') /* if * remap all mrnas of this gene */
	dictAdd (dict, cp, &mrna) ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &x1))
	continue ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &x2))
	continue ;
      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;
      if (! cp)
	continue ;
      dictAdd (dict, cp, &gene2) ;
      up = arrayp (aa, nn++, MRNA) ;
      up->mrna = mrna ;
      up->gene = gene ;
      up->gene2 = gene2 ;
      up->x1 = x1 ;
      up->x2 = x2 ;
      up->isRename = TRUE ;

      gb = tct->genes ? arrayp (tct->genes, gene2, GENE) : 0 ;
      ga = tct->genes ? arrayp (tct->genes, gene, GENE) : 0 ;
      gb = tct->genes ? arrayp (tct->genes, gene2, GENE) : 0 ; /* repeat if ga increased the array */
      if (ga && gb && ! gb->chrom)
	{ *gb = *ga ; gb->gene = gene2 ; }
    }
  ac_free (h) ;
  arraySort (aa, tctMrnaOrder) ;
  if (nn)
    {
      int i, m = 0, g = 0 ;
      MRNA *up ;
      KEYSET ks = tct->splitMrna2aa ;
      if (! ks)
	ks = tct->splitMrna2aa = keySetHandleCreate (tct->h) ;

      for (i = 0, up = arrp (aa, 0, MRNA) ; i < nn ; i++, up++)
	{
	  if (up->mrna != m)
	    { m = up->mrna ; keySet (ks, m) = i + 1 ; }
	  if (up->gene != g)
	    { g = up->gene ; keySet (ks, g) = i + 1 ; }
	}
    }
  return ;
} /* tctGetSplitMrnas */

/*************************************************************************************/
/*************************************************************************************/

static BOOL tctGetOneFastc (const TCT *tct, LANE *lane)
{
  AC_HANDLE h = ac_new_handle () ;
  char *fNam1 = tct->fastcFile || ! tct->fastcDir ? 0 : 
    hprintf (h, "%s/%s.fastc%s"
	     , tct->fastcDir
	     , dictName (tct->laneDict, lane->lane)
	     , tct->runTest ? "" : ".gz"
	     ) ;
  const char *fNam =  fNam1 ? fNam1 : tct->fastcFile ;
  ACEIN ai = aceInCreate (fNam, 0, h) ;
  ACEOUT ao = 0 ;
  int state = 0 ;
  DICT *dict = lane->cloneDict ;
  Array clones = lane->clones ;
  int i, k, n = 0, nn = 0, line = 0 ;
  CLONE *clone ;

  if (! fNam || ! ai)
    messcrash ("Cannot open lane fastc file -targetFileName %s\n", fNam ? fNam : "unspecified") ;
  if (tct->makeTest)
    ao = aceOutCreate (tct->makeTest, hprintf (h, ".%d.fastc", lane->lane), FALSE, h) ; 
  while (aceInCard (ai))
    {
      char *cq, *cp = aceInWord (ai) ;
      line++ ;
      if (! cp || ! *cp || *cp == '#')
	continue ;
      cq = 0 ;
      if (*cp == '>')
	{
	  int nCl = 0 ;
	  state = n = 0 ;
	  if (dictFind (dict, cp+1, &nCl))
	    {
	      state = 1 ; 
	      k = 0 ; /* current length */
	      clone = arrayp (clones, nCl, CLONE) ;
	    }
	}
      switch (state)
	{
	case 0:
	  break ;
	case 1: /* outside */
	  k = 0 ;
	  nn++ ;
	  clone->dna1 = arrayHandleCreate (n+2, char, lane->h) ;
	  state = 2 ;
	  break ;
	case 2: /* parsing read1 */
	  cq = strstr (cp, "><") ;
	  if (cq) *cq = 0 ;
	  n = strlen (cp) ;
	  if (n > 0)
	    {
	      arrayForceFeed (clone->dna1, k+n+2)  ; /* double zero terminate the dna */
	      arrayForceFeed (clone->dna1, k+n)  ;
	      memcpy (arrp (clone->dna1, k, char), cp, n) ; 
	      k += n ;
	    }
	  if (!cq)
	    break ;
	  cp = cq + 2 ;
	  state = 3 ;
	case 3: /* initiate read 2 */
	  k = 0 ;
	  n = strlen (cp) ;
	  nn++ ;
	  clone->dna2 = arrayHandleCreate (n+2, char, lane->h) ;
	  state = 4 ;
	  if (!n)
	    break ;
	  /* fall thru */
	case 4: /* parsing read3 */
	  n = strlen (cp) ;
	  if (n > 0)
	    {
	      arrayForceFeed (clone->dna2, k+n+2)  ;  /* double zero terminate the dna */
	      arrayForceFeed (clone->dna2, k+n)  ;
	      memcpy (arrp (clone->dna2, k, char), cp, n) ; 
	      k += n ;
	    }
	  state = 0 ; /* fastc filesgive the dna on a single line */
	  break ;
	}
    }  

  nn = arrayMax (clones) ;
  if (nn)
    for (i = 0, clone = arrp (clones, 0, CLONE) ; i < nn ; i++)
      {
	if (ao && clone[i].dna1)
	  {
	    aceOutf (ao, ">%s\n", dictName (dict, i)) ;
	    aceOutf (ao, "%s", arrp (clone[i].dna1, 0, char)) ;
	    if (clone[i].dna2)
	      aceOutf (ao, "><%s", arrp (clone[i].dna2, 0, char)) ;
	    aceOut (ao, "\n") ;
	  }
	if (clone[i].dna1)  dnaEncodeArray (clone[i].dna1) ;
	if (clone[i].dna2)  dnaEncodeArray (clone[i].dna2) ;
      }

  if (! tct->runTest)
    wego_log (hprintf (lane->h, "tctGetOneFastc found %d reads in %s\n", nn, fNam)) ;
  ac_free (ai) ;
  ac_free (ao) ;
  ac_free (h) ;
  return TRUE ;
} /* tctGetOneFastc */

/*************************************************************************************/
/*************************************************************************************/
/* shift sliding fusions towrds introns and report distance to intron on donor and acceptor sides */
static void tctShiftFusionTowardsIntrons (const TCT *tct, LANE *lane, GF *gf)
{
  int mrna1 = gf->t1 ; 
  int x1 = gf->x1 ;      
  int x2 = gf->x2 ;
  int a1 = gf->a1 ;      
  int a2 = gf->a2 ;

  int mrna2 = gf->t2 ; 
  int y1 = gf->y1 ;
  int y2 = gf->y2 ;
  int b1 = gf->b1 ;      
  int b2 = gf->b2 ;

  int dx = (y1 + y2 > x1 + x2 ? y1 - x2 - 1 : x2 - y1 - 1) ; /* number of sliding bases */
  if (y1 + y2 > x1 + x2)
    dx = (y1 < y2 ? y1 : y2) - (x2 > x1 ? x2 : x1) + 1 ;
  else
    dx = (x1 < x2 ? x1 : x2) - (y2 > y1 ? y2 : y1) + 1 ;
  gf->gap = dx ;

  if (tct->mrna2aa)
    {
      int j1 = keySet (tct->mrna2aa, mrna1) ;
      int j2 = keySet (tct->mrna2aa, mrna2) ;
      Array mmm = tct->mrnas ;

      if (j1 > 0 && mmm && j1 < arrayMax (mmm))
	{
	  MRNA *mp = arrp (mmm, j1 - 1, MRNA) - 1 ;
	  gf->d2donor = 999999 ;
	  while (mp++, mrna1 == mp->mrna) /* locate intron boundaries */
	    {
	      if (mp->isRename) continue ;
	      if (a1 < a2)  /* we are in read order and orientation */
		{
		  int da = 9999999, u2 = mp->x2 ;
		  if (a1 < u2 && a2 >= u2) 
		    {
		      da = u2 - a2 ;
		      if (dx < 0 && da <= 0) /* slide back */
			{
			  int dda = da - dx  < 0 ? dx : da ;
			  if (da - dda <= 10) /* i nearly reached the intron boundary */
			    {
			      if (gf->x1 < gf->x2) 
				gf->x2 += dda ; 
			      else
				gf->x2 -= dda ; 
			      gf->a2 += dda ;
			      da -= dda ;
			      dx -= dda ; /* this quantity of overlap has been used */
			      gf->d2donor = -da ; break ;
			    }
			}
		    }
		  else /*  if (dx > 0 && a1 < u2 && a2 + dx  > u2)  */
		    da = u2 - a2 ;		    
		  if (ABS(da) < ABS(gf->d2donor))
		    gf->d2donor = -da ;
		  if (ABS(da) < 10)
		    break ;
		}
	      else
		{
		  int da = 9999999, u1 = mp->x1 ;
		  if (a2 <= u1 && a1 > u1) 
		    {
		      da = a2 - u1 ;
		      if (dx < 0 && da <= 0) /* slide up */
			{
			  int dda = da - dx  < 0 ? dx : da ;
			  if (da - dda < 10) /* i nearly reached the intron boundary */
			    {
			      if (gf->x1 < gf->x2) 
				gf->x2 += dda ; 
			      else
				gf->x2 -= dda ; 
			      gf->a2 -= dda ;
			      da -= dda ;
			      dx -= dda ; /* this quantity of overlap has been used */
			      gf->d2donor = -da ; break ;
			    }
			}
		    }
		  else /* if (dx > 0 && a2 - dx < u1 && a1 > u1)  */
		    da = a2 - u1 ;
		  if (ABS(da) < ABS(gf->d2donor))
		    gf->d2donor = -da ;
		  if (ABS(da) < 10)
		    break ;
		}
	    }
	}
      if (j2 > 0 && mmm && j2 < arrayMax (mmm))
	{
	  MRNA *mp = arrp (mmm, j2 - 1, MRNA) - 1 ;
	  gf->d2acceptor = 999999 ;
	  while (mp++, mrna2 == mp->mrna) /* locate intron boundaries */
	    {
	      if (mp->isRename) continue ;
	      if (b1 < b2)  /* we are in read order and orientation */
		{
		  int db = 9999999, u1 = mp->x1 ;
		  if (b1 <= u1 && b2 > u1) 
		    {
		      db = u1 - b1 ;
		      if (dx < 0 && db >= 0) /* slide forward */
			{
			  int ddb = db + dx  > 0 ? -dx : db ;
			  if (db - ddb <= 10) /* i nearly reached the intron boundary */
			    {
			      if (gf->y1 < gf->y2) 
				gf->y1 += ddb ; 
			      else
				gf->y1 -= ddb ; 
			      gf->b1 += ddb ;
			      db -= ddb ;
			      dx -= ddb ;
			      gf->d2acceptor = -db ; break ;
			    }
			}
		    }
		  else /* if (dx > 0 && b1 - dx <= u1 && b2 > u1)  */
		    db = u1 - b1 ;
		  if (ABS(db) < ABS(gf->d2acceptor))
		    gf->d2acceptor = db ;
		  if (ABS(db) < 10)
		    break ;
		}
	      else
		{
		  int db = 9999999, u2 = mp->x2 ;
		  if (b2 < u2 && b1 > u2) 
		    {
		      db = b1 - u2 ;
		      if (dx < 0 && db >= 0) /* slide back */
			{
			  int ddb = db + dx  > 0 ? -dx : db ;
			  if (db - ddb < 10) /* i nearly reached the intron boundary */
			    {
			      if (gf->y1 < gf->y2) 
				gf->y1 += ddb ; 
			      else
				gf->y2 -= ddb ; 
			      gf->b1 -= ddb ;
			      db -= ddb ;
			      dx -= ddb ;
			      gf->d2acceptor = -db ; break ;
			    }
			}
		    }
		  else /* if (dx > 0 && b1 + dx >= u2 && b2 <  u2)  */
		    db = u2 - b1 ;
		  if (ABS(db) < ABS(gf->d2acceptor))
		    gf->d2acceptor = db ;
		  if (ABS(db) < 10)
		    break ;
		}
	    }
	}
    }
} /* tctShiftFusionTowardsIntrons */

/*************************************************************************************/
  
static void tctGetOneFusions (const TCT *tct, LANE *lane)
{
  HIT *hit ;
  Array hits = lane->hits ;
  int ii, iMax = arrayMax (hits) ;
  
  arraySort (hits, tctX1Order) ;
  for (ii = 1, hit = arrp (hits, 1, HIT) ; ii < iMax ; ii++, hit++)
    {
      HIT *old = hit - 1 ; /* recognize de-uno discovered gene fusions */
      if (old->clone == hit->clone && old->target != hit->target)
	{
	  int dc, dcOld, cc1, cc2, sens = 0, chromDistance = 0 ;
	  int wab ;
	  char * wo ;
	  int minChromDistance =100001 ; /* 100 k , the +1 at the end is needed in Z_genome case */
	  int minAli = tct->minAli, myMinAli ;
	  
	  if (minAli < 5)
	    minAli = 35 ;
	  dc = ii ; old = hit - 1 ;
	  while (--dc > 0)
	    {
	      if (old->clone != hit->clone ||
		  hit->isRead1 != old->isRead1  ||    /* not same read of the pair */
		  (old->target != hit->target && old->gene != hit->gene && old->x2 < hit->x1 + 20 && old->x1 < hit->x1 - 50)
		  )
		break ;
	      old-- ;
	    }
	  if (old->clone == hit->clone &&
	      old->gene > 0 && hit->gene > 0 &&  
	      hit->gene != old->gene &&
	      /* hit->x1 * old->x1 > 0 &&   same x strand */
	      hit->isRead1 == old->isRead1     /* same read of the pair */
	      )
	    { /* same read */
	      int a1, a2, b1, b2 ; 

	      HIT hit1, hit2 ;
	      int da, dc, gap, dc1, dc2 ;
	      BOOL sameChrom = FALSE ;

	      hit1 = *old ;
	      hit2 = *hit ;

	      if (hit1.x1 > hit1.x2)
		{
		  int x = hit1.x1 ; hit1.x1 = hit1.x2 ; hit1.x2 = x ;
		  int a = hit1.a1 ; hit1.a1 = hit1.a2 ; hit1.a2 = a ;
		}
	      if (hit2.x1 > hit2.x2)
		{
		  int x = hit2.x1 ; hit2.x1 = hit2.x2 ; hit2.x2 = x ;
		  int a = hit2.a1 ; hit2.a1 = hit2.a2 ; hit2.a2 = a ;
		}
	      if (hit1.c1 > hit1.c2)
		{
		  int c = hit1.c1 ; hit1.c1 = hit1.c2 ; hit1.c2 = c ;
		}
	      if (hit2.c1 > hit2.c2)
		{
		  int c = hit2.c1 ; hit2.c1 = hit2.c2 ; hit2.c2 = c ;
		}

	      if (hit1.x1 + hit1.x2 > hit2.x1 + hit2.x2)
		{ HIT h ; h = hit1 ; hit1 = hit2 ; hit2 = h ; }
	      /* we are in the read orientation and hit1 < hit2 */ 

	      a1 = hit1.a1 ; a2 = hit1.a2 ;
	      b1 = hit2.a1 ; b2 = hit2.a2 ;
	      
	      cc1 = (hit1.c1 > hit2.c1 ? hit1.c1 : hit2.c1 ) ;   /*  overlap start */
	      cc2 = (hit1.c2 < hit2.c2 ? hit1.c2 : hit2.c2 ) ;   /*  overlap end */

	      dc1 = hit1.c2 - hit1.c1 + 1 ;      /* chain length */
	      dc2 = hit2.c2 - hit2.c1 + 1 ;      /* chain length */
	      gap = hit2.x1 - hit1.x2 - 1 ;      /* gap if > 0   / overlap if < 0 */
	      dc = hit2.x2 - hit1.x1 + 1 ;           /* total aligned length */

	      if (dc1 + gap < 2 * dc1 / 3 || dc2 + gap < 2 * dc2 / 3)
		dc = 0 ; /* discard if overlap > 1/3 of shorter ali */

	      sens = (a1<a2 ? 0 : 2) + (b1 < b2 ? 0 : 1) ;  /* ++:0 +-:1  -+:2 --:3 */
	      
	      /* check the minimal distance, and relative orientation */
	      da = 0 ;
	      if (dc && tct->genes)
		{
		  GENE *ga = arrayp (tct->genes, hit1.gene, GENE) ;
		  GENE *gb = arrayp (tct->genes, hit2.gene, GENE) ;
		  ga = arrayp (tct->genes, hit1.gene, GENE) ; /* redo in case hit2.gene increased the array */
		  
		  sameChrom = TRUE ;
		  
		  if (ga && gb && ga->chrom == gb->chrom)
		    {
		      int u2, v1 ;
		      int aa2 = a2 + (a1 < a2 ?  gap : - gap) ; /* chrom distance of fused bases */

		      u2 = ga->a1 < ga->a2 ? ga->a1 + aa2 : ga->a1 - aa2 ; /* chrom position */
		      v1 = gb->a1 < gb->a2 ? gb->a1 + b1 : gb->a1 - b1 ; /* chrom position */


		      da = v1 - u2 ;
		      
		      switch (sens)
			{
			case 0: /* ++ */
			case 1: /* +- */
			  if (da > -30 && da < minChromDistance)  /* new intron or short deletion */
			    dc = 0 ;
			  break ;
			case 2: /* -+ */
			case 3: /* -- */
			  if (da < 30 && da > - minChromDistance)  /* new intron or short deletion */
			    dc = 0 ;
			  break ;
			}
		    }
		}

	      if (! tct->stranded)
		{
		  if (! sameChrom)
		    {
		      switch (sens)
			{
			case 2: /* -+ */
			  { HIT hh = hit1 ; hit1 = hit2 ; hit2 = hh ; sens = 1 ;  }
			  break ;
			}
		    }
		  
		  if (sens == 3) /* -- */
		    {
		      HIT hh = hit1 ; 
		      int x ;
		      
		      hit1 = hit2 ; hit2 = hh ; 
		      sens = 0 ;
		      x = hit1.x1 ; hit1.x1 = hit1.x2 ; hit1.x2 = x ;
		      x = hit1.a1 ; hit1.a1 = hit1.a2 ; hit1.a2 = x ;
		      x = hit2.x1 ; hit2.x1 = hit2.x2 ; hit2.x2 = x ;
		      x = hit2.a1 ; hit2.a1 = hit2.a2 ; hit2.a2 = x ;
		    }
		}

	      switch (sens)
		{
		case 0:
		  wo = "++" ;
		  break ;
		case 1:
		  wo = "+-" ;
		  break ;
		case 2:
		  wo = "-+" ;
		  break ;
		case 3:
		  wo = "--" ;
		  break ;
		}

	      myMinAli =  (dc > 500 && 100 > minAli) ?  100 : minAli ;
	      if (dc > 0 && dc1 > myMinAli && dc2 > myMinAli && hit1.gene && hit2.gene)
		{
		  GF *gf = arrayp (lane->geneFusions, arrayMax (lane->geneFusions), GF) ;
		  
		  gf->lane = lane ;
		  gf->clone = hit1.clone ;
		  gf->g1 = hit1.gene ;
		  gf->g2 = hit2.gene ;
		  gf->t1 = hit1.target ;
		  gf->t2 = hit2.target ;
		  gf->wo[0] = wo[0] ;
		  gf->wo[1] = wo[1] ;
		  gf->wo[2] = 0 ;
		  gf->isPair = FALSE ;
		  gf->score1 = hit1.score ;
		  gf->score2 = hit2.score ;
		  
		  gf->seq = 1 ;
		  gf->tag = hit1.mult ;
		  gf->n0 = gf->tag ;
		  gf->c1 = hit1.c1 ;
		  gf->c2 = hit1.c2 ;
		  gf->d1 = hit2.c1 ;
		  gf->d2 = hit2.c2 ;
		  
		  gf->x1 = ABS(hit1.x1) ;
		  gf->x2 = ABS(hit1.x2) ;
		  gf->y1 = ABS(hit2.x1) ;
		  gf->y2 = ABS(hit2.x2) ;
		  
		  gf->a1 = hit1.a1 ;
		  gf->a2 = hit1.a2 ;
		  gf->b1 = hit2.a1 ;
		  gf->b2 = hit2.a2 ;
		
		  gf->gap = gap ;
		  gf->chromDistance = chromDistance ;
		  gf->overhang1 = hit1.x1 < hit1.x2 ? hit1.overhang3p : hit1.overhang5p ;
		  gf->overhang2 = hit2.x1 < hit2.x2 ? hit2.overhang5p : hit2.overhang5p ;

		  tctShiftFusionTowardsIntrons (tct, lane, gf) ;
		}
	    }
	  
	  /* recognize gene fusions from pairs */
	  if (hit->isRead1)
	    continue ;
	  old = hit + 1 ;
	  if (ii < iMax - 1 &&  /* only analyze the highest x exon of read2, and try to connectit ot read 1 */
	      old->clone == hit->clone)
	    continue ;
	  /* look for read1 */
	  dc = ii - 1 ; old = hit - 1 ;
	  while (dc > 0 &&   /* move to highest x coordinates on second read */
		 old->clone == hit->clone &&
		 old->isRead1
		 )
	    { dc-- ; old-- ; }
	  /* eliminate read with a previous hit on hit->gene */
	  if (old->clone == hit->clone && 
	      old->gene > 0 && hit->gene > 0 && 
	      /* old->uu == 1 && hit->uu == 1 &&     the 2 reads align in a single gene */
	      hit->gene != old->gene &&
	      /* hit->x1 * old->x1 > 0 &&  same x strand */
	      hit->isRead1 != old->isRead1     /* opposite reads from same pair */

	      )
	    { /* opposite reads */
	      int a1, a2, b1, b2, x1, x2, y1, y2 ; 
	      HIT *up = old, *vp = hit ;
	      
	      a1 = old->a1 ; a2 = old->a2 ;
	      x1 = old->x1 ; x2 = old->x2 ;
	      b1 = hit->a1 ; b2 = hit->a2 ;
	      y1 = hit->x1 ; y2 = hit->x2 ;
	      
	      cc1 = hit->c1 ; if (old->c1 > cc1) cc1 = old->c1 ;   /*  overlap start */
	      cc2 = hit->c2 ; if (old->c2 < cc2) cc2 = old->c2 ;   /*  overlap end */
	      dc = hit->c2 >  hit->c1  ? hit->c2 - hit->c1 : hit->c1 - hit->c2 ;      /* chain length */
	      dcOld = old->c2 > old->c1 ? old->c2 - old->c1 : old->c1 - old->c2 ;      /* chain length */
	      
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
		  if (tct->genes)
		    {
		      GENE *ga = arrayp (tct->genes, up->gene, GENE) ;
		      GENE *gb = arrayp (tct->genes, vp->gene, GENE) ;
		      ga = arrayp (tct->genes, up->gene, GENE) ;  /* redo in case vp->gene increased the array */
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
	      if (tct->genes)
		{
		  GENE *ga = arrayp (tct->genes, up->gene, GENE) ;
		  GENE *gb = arrayp (tct->genes, vp->gene, GENE) ;
		  ga = arrayp (tct->genes, up->gene, GENE) ;  /* redo in case vp->gene increased the array */
		  if (ga && gb && ga->chrom == gb->chrom)
		    {
		      int dc2 = (gb->b1 + gb->b2)/2 - (ga->b1 + ga->b2)/2 ; /* distance of mid position a is left of b */
		      int ln = (gb->b2 - gb->b1 + ga->b2 - ga->b1)/2 ;
		      if (wab == 0 && gb->a2 < gb->a1) /* genes are on bottom strand */
			dc2 = - dc2 ;
		      if (dc2 > 0) /* tandem */
			{
			  chromDistance = dc2 - ln ;
			  if (chromDistance < minChromDistance)
			    dc = 0 ; /* eliminate this case, it should be detected as a new intron */
			}
		      else /* duplication */
			{
			  chromDistance = -dc2 - ln ;
			  if (chromDistance < 0)
			    dc = 0 ; /* eliminate this case */
			}
		      if (ga->b1 > gb->b1)
			chromDistance = gb->b1 - ga->b2 ;
		      else
			chromDistance = ga->b1 - gb->b2 ;
		      if (chromDistance < minChromDistance)
			dc = 0 ; /* eliminate this case */
		    }
		}
	      
	      myMinAli =  (dc + dcOld > 500 && 100 > minAli) ?  100 : minAli ;
	      if (dc > myMinAli && dcOld > myMinAli) /* do not check dcc overlap, not the same read */
		{
		  GF *gf = arrayp (lane->geneFusions, arrayMax (lane->geneFusions), GF) ;
		  
		  gf->lane = lane ;
		  gf->clone = hit->clone ;
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
		  gf->n0 = gf->tag ;
		  
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
    }
} /* tctGetOneFusions */

/*************************************************************************************/

static BOOL tctGetOneHit (const TCT *tct, LANE *lane)
{ 
  AC_HANDLE h = ac_new_handle () ;
  char *fNam1 = tct->hitFile ? 0 : hprintf (h, "%s/%s.hits%s"
					    , tct->hitDir
					    , dictName (tct->laneDict, lane->lane)
					    , tct->runTest ? "" : ".gz"
					   ) ;
  const char *fNam = tct->hitFile ? tct->hitFile : fNam1 ;
  ACEIN ai = aceInCreate (fNam, 0, h) ;
  ACEOUT ao = 0 ;
  int nn = 0, line = 0 ;
  int score, bestScore ;
  char buf[256] , buf2[256], oldBuf[256], bigBuf[21024], skipBuf[21024] ;
  char bestTarget_class[64] ;
  const char *tc = tct->target_class ;
  const char *target = tct->target ;
  int t1 = tct->t1 ? tct->t1 : tct->offset ;
  int t2 = tct->t2 ;
  DICT *dict = lane->cloneDict ;
  DICT *geneDict = tct->geneDict ;
  DICT *targetDict = tct->targetDict ;
  Array genes = 0 ;
  Array wig = lane->wig ;
  Array hits = lane->hits ;
  KEYSET brks = lane->brks ;
  Array brkIndels = lane->brkIndels ;
  KEYSET covers = lane->covers ;
  Array geneFusions = 0 ;
  KEYSET geneSupport = 0 ;
  KEYSET geneFusionKs1 = 0, geneFusionKs2 = 0 ; /* genes touched by read 1 and 2 of a pair */
  HIT *hit ;  
  int nSeq = 0, nHit = 0, lastGoodGF = 0 ;
  long int nBp = 0, nBrks = 0, nBrkIndels = 0, nOverhangs = 0 ;
  int minOverhang = tct->minOverhang ;
  int minAli = tct->minAli ;
  int minAliPerCent = tct->minAliPerCent ;
  BOOL subJump = FALSE, subSampling = tct->subSampling ? TRUE : FALSE ;
  char uno3p[50], uno3pOld[50], uno5p[50] ;
  BOOL needFastc = FALSE ;


  if (! ai)
    messcrash ("Cannot open lane hits file -targetFileName %s\n", fNam) ;
  
  if (tct->makeTest)
    ao = aceOutCreate (tct->makeTest, hprintf (h, ".%d.hits", lane->lane), FALSE, h) ; 

  geneFusions = lane->geneFusions = arrayHandleCreate (1000, GF, lane->h) ;
  lane->fusionHisto = keySetHandleCreate (lane->h) ;
  geneSupport = lane->geneSupport = keySetHandleCreate (lane->h) ;
  genes = lane->genes = arrayHandleCreate (128, GENE, tct->h) ;

  memset (bestTarget_class, 0, sizeof (bestTarget_class)) ;
  memset (buf, 0, sizeof (buf)) ;
  memset (buf2, 0, sizeof (buf2)) ;
  memset (oldBuf, 0, sizeof (oldBuf)) ;
  memset (bigBuf, 0, sizeof (bigBuf)) ;
  memset (skipBuf, 0, sizeof (skipBuf)) ;

  memset (uno3p, 0, sizeof (uno3p)) ;
  memset (uno5p, 0, sizeof (uno5p)) ;
  memset (uno3pOld, 0, sizeof (uno3pOld)) ;

  if (t2 == 0)
    t2 = 1 << 30 ;
  while (aceInCard (ai))
    {
      char tagBuf[64] ;
      BOOL isRead1 = TRUE ;
      int i, a1, a2, b1, b2, x1, x2, x11, x22, mult = 1, uu = 1, nErr = 0, ali, toBeAligned, gene= 0, trgt , target_class, trgtStart = 0 ;
      char *cp = aceInPos (ai) ; 
      
      line++ ;
      if (!cp || ! *cp || *cp == '#' || strlen (cp) > 20001)
	continue ;

      if (0 && ! strstr (aceInPos (ai), "2189:c>t"))
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

      if (subSampling)
	{
	  if (strcmp (buf, buf2))
	    {
	      strncpy (buf2, buf, 256) ; /* clone name name */
	      subJump = (nSeq++) % tct->subSampling ? TRUE : FALSE ;
	    }
	  if (subJump)
	    continue ;
	}

      if (strcmp (buf, oldBuf))
	{
	  if (! geneFusionKs1)
	    {
	      geneFusionKs1 = keySetHandleCreate (h) ;
	      geneFusionKs2 = keySetHandleCreate (h) ;
	    }
	  if (tct->geneFusion)
	    {
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
	      /********* fusions ******************/
	      if (tct->geneFusion && nHit > 1)
		tctGetOneFusions (tct, lane) ;

	      geneFusionKs1 = keySetReCreate (geneFusionKs1) ;
	      geneFusionKs2 = keySetReCreate (geneFusionKs2) ;

	      nHit = 0 ; 
	      arrayMax(hits) = 0 ; /* we are only interested by single fragments */	      
	    }
	  bestTarget_class[0] = 0 ; bestScore = 0 ;

	}

      aceInStep (ai, '\t') ; aceInInt (ai, &score) ;
      if (0 && score < bestScore)
	continue ;
      bestScore = score ;

      aceInStep (ai, '\t') ; aceInInt (ai, &mult) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &toBeAligned) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &ali) ;
      if (ali < minAli)
	continue ; 
      if (ali < 140 &&                    /* 100 per read translates into 200 for the pair which seems secure */
	  100 * ali < minAliPerCent * toBeAligned)
	continue ;
      x1 = x2 = x11 = x22 = 0 ;
      aceInStep (ai, '\t') ; aceInInt (ai, &x1) ; x11 = x1 ;
      aceInStep (ai, '\t') ; aceInInt (ai, &x2) ; x22 = x2 ;
      
      for(i = 8 ; i <= 8 ; i++) /* locate target class in column 8 */
	{ aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; }
      if (! cp)
	continue ;
      if (tct->target_class && strcasecmp (cp, tct->target_class))
	continue ;
      if (*bestTarget_class && strcasecmp (cp, bestTarget_class))
	continue ;
      strncpy (bestTarget_class, cp, 63) ; 
      if (tc && strcasecmp (cp, tc))
	continue ;
      if (score < bestScore)
	bestScore = score ;

      dictAdd (targetDict, cp, &target_class) ;
      aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* gene */
      gene = 0 ;
      dictAdd (geneDict, cp, &gene) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &uu) ; /* unicity */
      aceInStep (ai, '\t') ; cp =  aceInWord (ai) ;
      if (! cp)
	continue ;
      if (target && ! dictFind (targetDict, cp, &trgt))
	continue ;
      dictAdd (geneDict, cp, &trgt) ;
      trgtStart = keySet (tct->targetStart, trgt) ;
      a1 = a2 = 0 ;
      aceInStep (ai, '\t') ; aceInInt (ai, &a1) ; b1 = a1 ; 
      aceInStep (ai, '\t') ; aceInInt (ai, &a2) ; b2 = a2 ;
      if (tc && !strcmp (tc, "Z_genome"))
	{
	  GENE *ga = 0 ;
	  dictAdd (geneDict, messprintf("%s__%d", dictName (geneDict, trgt), a1/100000), &gene) ;
	  ga = arrayp (genes, gene, GENE) ;	
	  ga->chrom = trgt ;
	  ga->a1 = ga->b1 = (a1/100000) ;
	  ga->a2 = ga->b2 = ga->a1 + 1 ;
	  ga->a1 = ga->b1 = 100000 * ga->a1 ;
	  ga->a2 = ga->b2 = 100000 * ga->a2 ;
	}
      if (t1 + t2 && (a1 < t1 || a2 < t1 || a1 > t2 || a2 > t2))
	continue ;
      
      if (ao)
	aceOutf (ao, "%s\n", bigBuf) ;

      if (gene && tct->splitMrna2aa)
	{
	  Array mmm = tct->splitMrnas ;
	  int  mrna3 = 0, i3 = 0, ok = 0 ;
	  
	  i3 = keySet (tct->splitMrna2aa, gene) ;
	  
	  if (i3 > 0 && mmm && i3 < arrayMax (mmm))
	    {
	      MRNA *mp = arrp (mmm, i3 - 1, MRNA) - 1 ;
	      while (mp++, mrna3 == mp->mrna) /* scan the coordinates to locate the correct gene or the fusion or intron boundaries */
		{
		  if (! mp->isRename) continue ;
		  if (gene == mp->gene && gene == mp->mrna)
		    {
		      ok = 1 ;
		      gene = mp->gene2 ; /* global rename all transcripts of this gene */
		    }
		}
	    }
	  if (! ok && trgt)
	    {
	      mrna3 = trgt ;
	      
	      i3 = keySet (tct->splitMrna2aa, mrna3) ;
	      if (i3 > 0 && mmm && i3 < arrayMax (mmm))
		{
		  MRNA *mp = arrp (mmm, i3 - 1, MRNA) - 1 ;
		  int y1 = x1 < x2 ? x1 : x2 ;
		  int y2 = x1 < x2 ? x2 : x1 ;
		  int g1 = 0, g2 = 0 ;		  
		  
		  while (mp++, mrna3 == mp->mrna) /* scan the coordinates to locate the correct gene or the fusion or intron boundaries */
		    {
		      int u1 = y1 < mp->x1 ? mp->x1 : y1 ;
		      int u2 = y2 > mp->x2 ? mp->x2 : y2 ;

		      if (! mp->isRename) continue ;
		      if (mp->x1 == 0 && mp->x2 ==0)
			{ g2 = mp->gene2 ; break ; }
		      if (u1 < u2) /* some intersect */
			{ if (g2 != g1) g1 = g2 ; g2 = mp->gene ; }
		    }
		  if (g2 && ! g1)
		    gene = g2 ;
		  if (g2 && g1 && g1 != g2)
		    {
		      char buf[2048] ;
		      if (x1 > x2) 
			{ int g0 = g1 ; g1 = g2 ; g2 = g0 ;}
		      sprintf (buf, "%s__%s", dictName (geneDict, g1), dictName (geneDict, g2)) ;
		      dictAdd (geneDict, buf, &gene) ;  /* conflit possible en multi threaded */
		    }
		}
	    }
	}

      if (a1 > a2) { int a0 = a1 ; a1 = a2 ; a2 = a0 ; a0 = x1 ; x1 = x2 ; x2 = a0 ; }

      /* we are in the correct zone, we register this hit */
      strcpy (oldBuf, buf) ;
      dictAdd (dict, buf, &nn) ;
      if (gene > 0) keySet (geneSupport, gene)++ ;

      hit = arrayp (hits, nHit++, HIT) ;
      hit->lane = lane->lane ;
      hit->clone = nn ;
      hit->isRead1 = isRead1 ;
      hit->gene = gene ;
      hit->target_class = target_class ;
      hit->target = trgt ;
      hit->mult = mult ;
      hit->uu = uu ;
      hit->a1 = a1 ; hit->a2 = a2 ; /* local coords */
      hit->x1 = x1 ; hit->x2 = x2 ;
      if (x1 < x2) hit->isDown = TRUE ;
      hit->isRead1 = isRead1 ;

      if (tct->deUno && a2 > a1 + 8)
	{
	  array (wig, a1 + 8 + trgtStart, int) += mult ;
	  array (wig, a2 - 8 + 1 + trgtStart, int) -= mult ;
	  
	  if (0) fprintf (stderr, "#0 wig[40]=%d a1=%d a2=%d mult=%d\n", array (wig, 40, int), a1, a2, mult) ;
	}
      
       /* register the cover wiggle */ 
      if (! tct->geneFusion)
	{
	  if (a1 < t1) a1 = t1 ;
	  if (a2 > t2) a2 = t2 ;
	  for (i = a1 ; i <= a2 ; i++)
	    {
	      int k = keySet (covers, i - t1 + trgtStart) ;
	      k += mult ; nBp += mult ;
	      keySet (covers, i - t1 + trgtStart) = k ;
	    }
	}
      aceInStep (ai, '\t') ; aceInInt (ai, &nErr) ; 
       /* register the error wiggle */  
       for(i = 15 ; i <= 17 ; i++) /* locate target errors in column 17 */
	{ aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; }
       if (tct->intron_only && cp && cp[1]) continue ;
       if (0 || tct->intron_only) cp = 0 ; /* mask the local errors */
       if (tct->geneFusion) 
	 {
	   if (isRead1)
	     keySetInsert (geneFusionKs1, gene) ;
	   else
	     keySetInsert (geneFusionKs2, gene) ;
	 }
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
		       x-- ; /* position on the hook base */
		       if (tct->wantBrks)
			 {
			   BKID *b =  arrayp (brkIndels, nBrkIndels++, BKID) ; 
			   
			   b->target = hit->target ;
			   b->pos = x - t1 ;
			   b->daa = 3 ;
			   b->dxx = - dx  ;
			   b->mult = mult ;
			 }	   
		     }
		   tagBuf[0] = 0 ;
		   if (tct->deUno)
		     {
		       if (dx == 0)
			 {
			   int k ;
			   char buf[64] ;
			   UNO *uno ;
			   
			   sprintf(buf, "%d:Sub:%c:%c", x,ace_upper(cp[0]), ace_upper(cp[2])) ;
			   sprintf(tagBuf, "%c2%c\t0\t-", ace_upper(cp[0]), ace_upper(cp[2])) ;
			   dictAdd (lane->deUnoDict, buf, &k) ;
			   uno = arrayp (lane->deUno, arrayMax (lane->deUno), UNO) ;
			   uno->target = trgt ;
			   uno->a1 = x - 1 ;
			   uno->a2 = x + 1 ;
			   uno->da = 0 ;
			   uno->nt += mult ;
			   uno->ns++ ;
			   uno->type = k ;
			   dictAdd (lane->deUnoDict, tagBuf, &(uno->tag)) ;
			 }
		       else if (dx > 0)
			 {
			   int k ;
			   char buf[64], buf2[dx+1] ; 
			   UNO *uno ;
			   int da = dx ;
			   int xx = x - 1 - t1 + keySet (tct->targetStart, trgt) ;
			   char *cr = arrp (tct->dna, xx, char) ; /* last matching base */
			   cp += dx ; 
			   k = da - 1 ;
			   while (dnaDecodeChar[(int)cr[0]] == cp[k]) { x-- ; cr-- ; k-- ; if(k<0) k = da - 1 ; }
			   for (i = 0 ; i < da && i < INSMAX ; i++)
			     { buf2[i] = ace_upper(cp[k]) ; k = (k+1) % da ; }
			   buf2[i] = 0 ;
			   if (dx ==1)
			     {
			       char cc = ace_upper(buf2[0]) ;
			       sprintf(tagBuf, "Ins%c\t1\t%c", cc, cc) ; 
			       sprintf(buf, "%d:Ins:%c:%c%c", x, dnaDecodeChar[(int)cr[0]], dnaDecodeChar[(int)cr[0]], cc) ;
			     }
			   else if (dx < 30)
			     {
			       sprintf (tagBuf, "Multi_insertion\t%d\t%s", dx, buf2) ;
			       sprintf (buf, "%d:Ins_%d:%c:%c%s", x, da, dnaDecodeChar[(int)cr[0]], dnaDecodeChar[(int)cr[0]], buf2) ;
			     }
			   else
			     {
			       sprintf (tagBuf, "Multi_insertion\t%d\t%s", dx, buf2) ;
			       sprintf(buf, "%d:Ins_%d:%c:", x, da, dnaDecodeChar[(int)cr[0]]) ;
			     }
			   dictAdd (lane->deUnoDict, buf, &k) ;
			   uno = arrayp (lane->deUno, arrayMax (lane->deUno), UNO) ;
			   uno->target = trgt ;
			   uno->a1 = x  ;
			   uno->a2 = x + 1 ;
			   uno->da = da ;
			   uno->ddx = da ;
			   uno->nt += mult ;
			   uno->ns++ ;
			   uno->type = k ;
			   dictAdd (lane->deUnoDict, tagBuf, &k) ;
			   uno->tag = k ;
			 }
		       else if (dx < 0)
			 {
			   int k, ddw = 0 ;
			   char buf[64] ;
			   UNO *uno ;
			   int da = -dx ;
			   int xx = x - 1 - t1 + trgtStart ;
			   char *cr = arrp (tct->dna, xx, char) ; /* last matching base */
			   while (xx > 2 && cr[0] == cr[da]) { x-- ; xx-- ; cr-- ; ddw-- ;}
			   if (ddw ) /* adjust the wiggle */
			     {
			       array (wig, xx, int) -= mult ;
			       array (wig, xx + da, int) += mult ;
			       array (wig, xx - ddw, int) += mult ;
			       array (wig, xx + da -ddw, int) -= mult ;
			     }
			   
			   if (da == 1)
			     {
			       char cc = ace_upper(dnaDecodeChar[(int)cr[1]]) ;
			       sprintf(tagBuf, "Del%c\t1\t%c", cc, cc) ;
			       sprintf (buf, "%d:Del:%c%c:%c", x, dnaDecodeChar[(int)cr[0]], cc, dnaDecodeChar[(int)cr[0]]) ;
			     }
			   else if (da < DELMAX)
			     {
			       int k0 ;
			       sprintf (buf, "%d:Del_%d:%c", x, da, dnaDecodeChar[(int)cr[0]]) ;
			       k0 = k = strlen (buf) ; 
			       for (i = 1 ; i <= da ; i++)
				 buf[k++] = ace_upper (dnaDecodeChar[(int)cr[i]]) ;
			       buf[k] = 0 ;
			       sprintf (tagBuf, "Multi_deletion\t%d\t%s", da, buf+k0) ;
			       buf[k++] = ':' ;
			       buf[k++] = dnaDecodeChar[(int)cr[0]] ;
			       buf[k++] = 0 ;
			     }
			   else
			     {
			       sprintf (tagBuf, "Multi_deletion\t%d", da) ;
			       sprintf (buf, "%d:Del_%d:%c:", x, da, dnaDecodeChar[(int)cr[0]]) ;
			     }
			   dictAdd (lane->deUnoDict, buf, &k) ;
			   uno = arrayp (lane->deUno, arrayMax (lane->deUno), UNO) ;
			   uno->target = trgt ;
			   uno->a1 = x ;
			   uno->a2 = x  + da + 1 ;
			   uno->da = da ;
			   uno->ddx = -da ;
			   uno->nt += mult ;
			   uno->ns++ ;
			   uno->type = k ;
			   dictAdd (lane->deUnoDict, tagBuf, &k) ;
			   uno->tag = k ;
			 }
		     }
		 }
	       if (tct->wantBrks)
		 {
		   k = keySet (brks, x - t1 + trgtStart) ;
		   k += mult ; nBrks += mult ;
		   keySet (brks, x - t1 + trgtStart) = k ;
		 }
	     }
	   cp = cq ? cq + 1 : 0 ;
	 } 
       /* register the significant overhangs */
       aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* 5p overhang */
       if (cp && *cp && tct->geneFusion && hit)
	 dictAdd (tct->overhangDict, cp, &(hit->overhang5p)) ;
       uno5p[0] = 0 ; if (cp) memcpy (uno5p, cp, 49) ; uno5p[49] = 0 ;
       if (! tct->geneFusion && tct->wantBrks &&
	   (
	    (cp && *cp && b1 > t1 && b1 < t2 && *cp == ace_lower (*cp) && strlen (cp) > minOverhang) 
	    ||
	    (x11 > 20 && (! cp || ! cp[1] || *cp == ace_lower (*cp)))
	    ) 
	   )
	 { 
	   int k = keySet (brks, b1 - t1 + trgtStart) ;
	   k += mult ; nOverhangs += mult ;
	   keySet (brks, b1 - t1 + trgtStart) = k ;
	   if (1)
	     {
	       BKID *b =  arrayp (brkIndels, nBrkIndels++, BKID) ; 
	       b->target = trgt ;
	       b->pos = b1 - t1 ;
	       b->daa = 1 ;
	       b->dxx =   (x1 == x11) ? - 1 : 1 ;
	       b->mult = mult ;
	     }	   
	 }
       aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* 3p overhang */
       if (cp && *cp && tct->geneFusion && hit)
	 dictAdd (tct->overhangDict, cp, &(hit->overhang3p)) ;
       strcpy (uno3pOld, uno3p) ;
       uno3p[0] = 0 ; if (cp) memcpy (uno3p, cp, 49) ; uno3p[49] = 0 ;
       if (! tct->geneFusion && tct->wantBrks &&
	   (
	    (cp && *cp && b2 > t1 && b2 < t2 && *cp == ace_lower (*cp) && strlen (cp) > minOverhang)
	    ||
	    ( x22 < toBeAligned- 20 && (! cp || ! cp[1] || *cp == ace_lower (*cp)))
	    )
	   )
	 { 
	   int k = keySet (brks, b2 - t1 + trgtStart) ;
	   k += mult ; nOverhangs += mult ;
	   keySet (brks, b2 - t1 + trgtStart) = k ;
	   if (1)
	     {
	       BKID *b =  arrayp (brkIndels, nBrkIndels++, BKID) ; 
	       b->target = trgt ;
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
       aceInStep (ai, '\t') ; hit->c1 = x1 ; aceInInt (ai, &(hit->c1)) ; /* col 27: chain1 */
       aceInStep (ai, '\t') ; hit->c2 = x2 ; aceInInt (ai, &(hit->c2)) ; /* col 28: chain2 */



       /* register de-uno discovered deletions */
       if (! tct->geneFusion && nHit > 1)
	 {
	   HIT *old = hit - 1 ;
	   int da, dx, protect = 0, c1, c2 ;
	   
	   if (old->clone == hit->clone &&
	       hit->target == old->target &&
	       hit->isDown == old->isDown &&
	       hit->isRead1 == old->isRead1 &&
	       (
		1 ||
		(old->a1 < hit->a1 && hit->x1 > old->x1) ||
		(old->a1 > hit->a1 && hit->x1 < old->x1)
		)
	       )
	     { /* same read */
	       BOOL isDown = TRUE ;
	       int a1, a2, b1, b2, x1, x2, y1, y2 ;
	       int dup= 0 ;

	       if (old->a1 < hit->a1)
		 { 
		   a1 = old->a1 ; a2 = old->a2 ;
		   x1 = old->x1 ; x2 = old->x2 ;
		   b1 = hit->a1 ; b2 = hit->a2 ;
		   y1 = hit->x1 ; y2 = hit->x2 ;
		 }
	       else
		 { 
		   b1 = old->a1 ; b2 = old->a2 ;
		   y1 = old->x1 ; y2 = old->x2 ;
		   a1 = hit->a1 ; a2 = hit->a2 ;
		   x1 = hit->x1 ; x2 = hit->x2 ;
		 }

	       if (x1 > x2)
		 isDown = FALSE ;
	       c1 = a2 ; c2 = b1 ;
	       if (isDown)
		 {
		   dx =  y1 - x2 ;
		   if (y1 < y2) dx = dx - 1 ; else  dx = -dx - 1 ; 
		 }
	       else
		 {
		   dx =  -y1 + x2 ;
		   if (y1 > y2) dx = dx - 1 ; else  dx = -dx - 1 ; 
		 }
	       /* clean up x duplications */
	       
	       if (b1 > a1)
		 {
		   if (x1 < x2)
		     {
		       if (y1 < x2 + 1)
			 {
			   dx = (x2 + 1) - y1 ;
			   protect = dx ;
			   if (! nErr)
			     {
			       x2 -= dx ; a2 -= dx ; dx = 0 ;
			     }
			   else
			     {
			       y1 += dx ; b1 += dx ; dx = 0 ;
			     }
			 }
		     }
		   else
		     {
		       if (y1 > x2 - 1)
			 {
			   dx = y1 - (x2 - 1) ;
			   protect = dx ;
			   x2 += dx ; a2 -= dx ; dx = 0 ;
			   dx = 0 ;
			 }
		     }
		 }

	       da = b1 - a2 ;
	       if (b1 > a2)
		 {
		   da-- ;
		   if (isDown)
		     {
		       if (y1 > y2 && (y1 >= x1 - tct->Delta || x2 < y2 + tct->Delta))
			 da = dx = 0 ;
		       if (y1 < y2 && (x2 > y2 - tct->Delta || y1 < x1 + tct->Delta))
			 da = dx = 0 ;
		     }
		   else
		     {
		       if (y1 < y2 && (y1 < x1 + tct->Delta || x2 > y2 - tct->Delta))
			 da = dx = 0 ;
		       if (y1 > y2 && (x2 < y2 + tct->Delta || y1 > x1 - tct->Delta))
			 da = dx = 0 ;
		     }
		 }
	       if (b1 < a1 + tct->Delta || a2 >= b2 -tct->Delta)
		 da = dx = 0 ;

	       if (dx == 0 && a2 >= b1)
		 { dup = a2 - b1 + 1 ; a2 -= dup  ; x2 += (x1 < x2) ? - dup : dup ; da = dx = 0 ; /* b1 - a2 - 1 ; dx = y1 - x2 + (x1 < x2 ? -1 : + 1 ) - dup ; */ }

	       else if (dx > 0 && a2 >= b1)
		 { 
		   if (isDown)
		     { dx = dx + a2 - b1 + 1 ; b1 = a2 + 1 ; y1 -= da ; }
		   else
		     { dx = dx + a2 - b1 + 1 ; b1 = a2 + 1 ; y1 += da ;}
		   dup = 0 ; da = 0 ; 
		 }
	       if (dx < 0 && da > -dx)
		 { a2 += dx ; x2 += isDown ? dx : -dx ; dx = 0 ; da -= dx ; } 

	       if (tct->wantBrks)
		 if (da || dx || dup)
		   {
		     int k ;
		     BKID *b =  arrayp (brkIndels, nBrkIndels++, BKID) ; 
		     b->target = hit->target ;
		     b->pos = a2 + protect - t1 ;
		     b->mult = mult ;
		     b->daa = da ;
		     b->dxx = dx ;  /* da = 0 dx = 0 dup>0 represents a duplication */
		     b->c1 = c1 - t1 ;
		     b->c2 = c2 - t1 ;
		     k = keySet (brks, b1 - t1) ;
		     k += mult ; nOverhangs += mult ;
		     keySet (brks, b1 - t1) = k ;
		     
		     if (0)
		       fprintf (stderr, "gethits found c1/c2 = %d %d\n", c1, c2);
		     /* if we have a duplication, we must be out of it on both sides */
		     if (dx - dup < 0 && -dx + dup > tct->Delta)
		       b1 += (-dx  - dup - tct->Delta) ;
		     
		     if (dx -dup < 0)
		       a2 += dx - dup - 1 ; 
		     k = keySet (brks, a2 - t1 + trgtStart) ;
		     k += mult ; nOverhangs += mult ;
		     keySet (brks, a2 - t1 + trgtStart) = k ;
		     
		     if (a2 < b1 && a2 > b1 - 60)
		       {
			 int a22 ;
			 for (a22 = a2 + 1 ; a22 < b1 ; a22 += tct->Delta/2)
			   {
			     k = keySet (brks, a22 - t1 + trgtStart) ;
			     k += mult ; nOverhangs += mult ;
			     keySet (brks, a22 - t1 + trgtStart) = k ;
			   }
		       }
		   }	   

	       if (tct->deUno)
		 if (da || dx || dup) /* deletion , insertion duplication */
		   {
		     int k ;
		     UNO *uno ;
		     int xx = a2 - 1 - t1 + keySet (tct->targetStart, trgt) ;
		     char *cr = arrp (tct->dna, xx, char) ;
		     int ddx = dx + dup ;
		     if (dup < 0) dup = -dup ;
		     if (ddx < 0) ddx = -ddx ;
		     if (da > 100 && a2 < 75 && a2 + protect >= 75)
		       {
			 int dd = 75 - a2 ;
			 a2 += dd ; x2 += isDown ? dd : -dd ; b1 += dd ; y1 += isDown ? dd : -dd ;
		       }
		     
		     if (da)
		       {
			 sprintf (tagBuf, "Multi_deletion\t%d", da) ;
			 sprintf (buf, "%d:Del_%d::%c", a2, da, dnaDecodeChar[(int)cr[0]]) ;
		       }
		     else if (ddx)  /* da = 0 dx < 0 represents a duplication */
		       {
			 sprintf (tagBuf, "Multi_insertion\t%d", ddx) ;
			 sprintf (buf, "%d:Ins_%d:%c:", a2 + (isDown && dx > 0 ? dup : 0) + (isDown && dx==0 ? 0 : 0), ddx, dnaDecodeChar[(int)cr[0]]) ;
		       }
		     dictAdd (lane->deUnoDict, buf, &k) ;
		     uno = arrayp (lane->deUno, arrayMax (lane->deUno), UNO) ;
		     uno->target = hit->target ;
		     uno->a1 = a2 + (isDown && dx > 0 ? dup : 0) + (isDown && dx==0 ? 0 : 0) ;
		     uno->a2 = b1 + (isDown && dx > 0 ? dup : 0) + (isDown && dx==0 ? 0 : 0) ;
		     uno->ddx = ddx - da ;
		     uno->da = da ? da : -dx ;
		     uno->dup = dup ;
		     uno->nt += mult ;
		     uno->ns++ ;
		     uno->type = k ;
		     dictAdd (lane->deUnoDict, tagBuf, &k) ;
		     uno->tag = k ;
		     if (da > 0)
		       {
			 if (da < 35)
			   {
			     int i = da ;
			     char *cp = arrp (tct->dna, a2, char) ;
			     char buf[da+1] ;
			     memcpy (buf, cp, da) ;
			     buf[da] = 0 ;
			     cp = buf ;
			     while (i--)
			       {*cp = ace_upper (dnaDecodeChar[(int)*cp]) ; cp++; }
			     dictAdd (lane->deUnoDict, buf, &k) ;
			     uno->insert = k ;
			   }
			 else
			   {
			     int i = da ;
			     char *cp = arrp (tct->dna, a2, char) ;
			     char *cq, buf[18] ;
			     
			     cq = buf ;
			     for (i = 0 ; i < 8 ; i++)
			       cq[i] = ace_upper (dnaDecodeChar[(int)cp[i]]) ;
			     cq[8] = '_' ;
			     for (i = 0 ; i < 8 ; i++)
			       cq[9+i] = ace_upper (dnaDecodeChar[(int)cp[da-8+i]]) ;
			     cq[17] = 0 ;
			     dictAdd (lane->deUnoDict, buf, &k) ;
			     uno->insert = k ;
			   }

		       }
		     else if (isDown && dup == 0 && dx > 0 && dx < 5 && uno3pOld[0] && strlen (uno3pOld) >= dx)
		       {
			 uno3pOld[dx] = 0 ; 
			 char *cp = uno3pOld ;
			 while (1 && *cp)
			   {*cp = ace_upper (*cp) ; cp++; }
			 dictAdd (lane->deUnoDict, uno3pOld, &k) ;
			 uno->insert = k ;
		       }
		     else if (!isDown && dup == 0 && dx > 0 && dx < 5 && uno5p[0] && strlen (uno5p) >= dx)
		       {
			 uno5p[dx] = 0 ; 
			 char *cp = uno5p ;
			 while (1 && *cp)
			   {*cp = ace_upper (*cp) ; cp++; }
			 dictAdd (lane->deUnoDict, uno5p, &k) ;
			 uno->insert = k ;
		       }
		     else if (dx == 0 && dup)
		       {
			 int i = dup ;
			 char *cp = arrp (tct->dna, a2, char) ;
			 char buf[dup+1] ;
			 memcpy (buf, cp, dup) ;
			 buf[dup] = 0 ;
			 cp = buf ;
			 while (i--)
			   {*cp = ace_upper(dnaDecodeChar[(int)*cp]) ; cp++; }
			 dictAdd (lane->deUnoDict, buf, &k) ;
			 uno->insert = k ;
		       }
		     else if (dup && dx > 0 && dx < 500 && uno3pOld[0] && strlen (uno3pOld) >= dx)
		       {
			 int i = dup ;
			 char *cp ;
			 char buf[dx + dup+1] ;
			 uno3pOld[dx] = 0 ; 
			 memcpy (buf, uno3pOld,dx) ;
			 cp = arrp (tct->dna, a2 , char) ;
			 memcpy (buf + dx, cp, dup) ;
			 cp = buf + dx ;
			 i = dup ;
			 while (i--)
			   {*cp = dnaDecodeChar[(int)*cp] ; cp++; }
			 buf[dup+dx] = 0 ;
			 cp = buf ;
			 while (*cp)
			   { *cp = ace_upper(*cp) ; cp++ ; }
			 dictAdd (lane->deUnoDict, buf, &k) ;
			 uno->insert = k ;
		       }
		     else if (dx  && !uno->insert && !uno->clone)
		       {
			 needFastc = TRUE ;
			 uno->isRead1 = hit->isRead1 ;
			 uno->isDown = isDown ;
			 if (isDown)
			   {
			     uno->x1 = x2 + 1 ;
			     uno->x2 = y1 - 1 ;
			   }
			 else
			   {
			     uno->x1 = y1 + 1 ;
			     uno->x2 = x2 - 1 ;
			   }
			 dictAdd (lane->cloneDict, dictName (dict, hit->clone), &uno->clone) ;
		       }
		   }	   
	     }
	 }
     }
  /********* last fusion ******************/
  if (1 && tct->geneFusion && nHit > 1)
    tctGetOneFusions (tct, lane) ;

  if (needFastc && (tct->fastcFile || tct->fastcDir))
    {
      UNO *uno ;
      int ii, iiMax = arrayMax (lane->deUno) ;
      tctGetOneFastc (tct, lane) ;

      for (ii = 0, uno = arrp (lane->deUno, 0, UNO) ; ii < iiMax ; ii++, uno++)
	if (uno->ddx && uno->clone && ! uno->insert)
	  {
	    CLONE *clone = arrayp (lane->clones, uno->clone, CLONE) ;
	    Array dna = uno->isRead1 ? clone->dna1 : clone->dna2 ;

	    if (dna && arrayMax (dna) > uno->x2)
	      {
		int dup = uno->dup ;
		char *cp, *cpg = arrp (tct->dna, uno->a2 , char) ;
		int x1 = uno->x1 ;
		int x2 = uno->x2 ;
		int dx = x2 > x1 ? x2 - x1 + 1 : x1 - x2 + 1 ;
		int i ;
		char buf[dup + dx + 1] ;
		
                if (dx) dup = 0 ;
		dna = arrayCopy (dna) ;
		if (1 && ! uno->isDown)
		  {
		    int x0 = x1 ;
		    reverseComplement (dna)  ;
		    x1 = arrayMax (dna) - x2 + 1 ;
		    x2 = arrayMax (dna) - x0 + 1 ;
		  }
		dx = x2 - x1 + 1 ;


		cp = arrp (dna, x1 - 1, char) ;
		memcpy (buf, cp, dx) ;
		if (dup) 
		  memcpy (buf+dx, cpg, dup) ;
		cp = buf ;
		i = dx + dup ;
		buf[i] = 0 ;
		while (i--)
		  {*cp = dnaDecodeChar[(int)*cp] ; cp++; }
		cp = buf ;
		while (*cp)
		  { *cp = ace_upper(*cp) ; cp++ ; }
		dictAdd (lane->deUnoDict, buf, &uno->insert) ;
		arrayDestroy (dna) ;
	      }
	  }
    }


  if (! tct->runTest)
    wego_log (hprintf (lane->h, "tctGetOneHit found %d hits, covering %ld bp and %ld brks (max %d) %ld overhangs in %s\n", nHit, nBp, nBrks, brks ? arrayMax(brks) : 0, nOverhangs, fNam)) ;
  ac_free (ai) ;
  ac_free (ao) ;
  ac_free (h) ;
  return TRUE ;
} /* tctGetOneHit */

/*************************************************************************************/
/* output: brks, deUno report and all the python  QC */
static BOOL tctGetOneSamHit (const TCT *tct, LANE *lane)
{ 
  AC_HANDLE h = ac_new_handle () ;
  const char *tcl = tct->target_class ? tct->target_class :  "genome"  ;
  const char *tcl2 = tcl ? strchr (tcl, '_') + 1 : 0 ;      /* split av out of ET_av */
  const char *tcl3 = tcl2 ? tcl2 : tcl ;
  char *fNam1 =  tct->SAM_file  ? 0 : hprintf (h, "%s/%s.%s.mbhits%s"
					       , tct->hitDir ? tct->hitDir : ""
					       , dictName (tct->laneDict, lane->lane)
					       , tcl3
					       , tct->runTest ? "" : ".gz"
					       ) ;
  const char *fNam = tct->SAM_file ? tct->SAM_file : fNam1 ;
  ACEIN ai = aceInCreate (fNam, tct->gzi, h) ;
  ACEOUT ao = 0 ;
  int nn = 0 ;
  static int nnn = 0 ;
  int iCigar, dummy = 0 ;
  char buf[DELMAX+256], buf2[DELMAX+256] , bigBuf[1024], *cigar ;
  char tagBuf[DELMAX+256] ;
  int t1 = tct->t1 ;
  int t2 = tct->t2 ;
  DICT *dict = lane->cloneDict ;
  DICT *targetDict = tct->targetDict ;
  Array hits = lane->hits ;
  Array wig = lane->wig ;
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
  BOOL subJump = FALSE, subSampling = tct->subSampling ? TRUE : FALSE ;
  int nSeq = 0 ;

  if (! ai)
    messcrash ("Cannot open lane hits file -targetFileName %s\n", fNam) ;
  
  if (tct->makeTest)
    ao = aceOutCreate (tct->makeTest, hprintf (h, ".%d.hits", lane->lane), FALSE, h) ; 

  memset (buf, 0, sizeof (buf)) ;
  memset (buf2, 0, sizeof (buf2)) ;
  memset (bigBuf, 0, sizeof (bigBuf)) ;

  if (t2 == 0)
    t2 = 1 << 30 ;
  while (aceInCard (ai))
    {
      BOOL isRead1 = TRUE ;
      int i, flag, a1, a2, x1, x2, mult = 1, uu = 1, score = 0, nerr = 0, ali = 0, trgt = 0 ;
      char *seq, *cp = aceInPos (ai) ; 

      if (!cp || ! *cp || *cp == '@' || *cp == '#')
	continue ;

      cp = aceInWord (ai) ;
      strncpy (bigBuf, aceInPos (ai), 1023) ; /* fragment name */
      strncpy (buf, cp, 255) ;
      if (subSampling)
	{
	  if (strcmp (buf, buf2))
	    {
	      strncpy (buf2, buf, 256) ; /* clone name name */
	      subJump = (nSeq++) % tct->subSampling ? TRUE : FALSE ;
	    }
	  if (subJump)
	    continue ;
	}
      cp = strchr (buf, '#') ;
      if (cp && sscanf (cp+1, "%d", &i) == 1 && i > 0)
	mult = i ;
      


      aceInStep (ai, '\t') ; aceInInt (ai, &flag) ; /* flag hides the read 1/2 distinction */
      isRead1 = (flag & 128 ? FALSE : TRUE) ;

      aceInStep (ai, '\t') ; cp =  aceInWord (ai) ; /* target */
      if (tct->target && ! dictFind (targetDict, cp, &trgt))
	continue ;
      dictAdd (tct->geneDict, cp, &trgt) ;

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
		if (cgr->a2)
		  {
		    int a1 = cgr->a1 ;
		    int a2 = cgr->a2 ;
		    if (a1 > a2) 
		      { int a0 = a1 ; a1 = a2 ; a2 = a0 ; }
		    if (wig && a2 > a1 + 8)
		      {
			array (wig, a1 + 8, int) += mult ;
			array (wig, a2 - 8 + 1, int) -= mult ;
		      }
		  }
		if (0) fprintf (stderr, "#0 wig[40]=%d a1=%d a2=%d mult=%d\n", array (wig, 40, int), cgr->a1, cgr->a2, mult) ;
		break ;
	      case 'D':
	      case 'I':
	      case 'N':/* intron */
		{
		  int k, a1, a2, da = 0 ;
		  UNO *uno ;
		  char *cq, buf[DELMAX+64] ;
		  BKID *b =  brkIndels ? arrayp (brkIndels, nBrkIndels++, BKID) : 0 ;
		  int ddw, dda, dx = cgr->x2 - cgr->x1 ;

		  tagBuf[0] = 0 ;
		  if (dx < 0)
		    { cgr->x2 -= dx ; cgr->a2 -= dx ; cgr->a1 += dx ; cgr->x1 += dx ; dx = 0 ; }
		  /* cgr->a1 a2 = first last base of intron, g._.g of gt_ag */
		  b->target = trgt ;
		  b->pos = cgr->a1 - 1 - t1 ;
		  b->mult = mult ;
		  b->daa = cgr->a2 - cgr->a1 + 1 ;
		  b->dxx = 0 ;  /* da = 0 dx < 0 represents a duplication */

		  if (tct->wantBrks)
		    {
		      k = keySet (brks, cgr->a1 -1 - t1) ;
		      k += mult ; nOverhangs += mult ;
		      keySet (brks, cgr->a1 -1 - t1) = k ;
		      
		      k = keySet (brks, cgr->a2 - t1) ;
		      k += mult ; nOverhangs += mult ;
		      keySet (brks, cgr->a2 - t1) = k ;
		    }

		  k = cgr->a2 - t1 ;
		  a1 = cgr->a1 ; a2 = cgr->a2 ;
		  da = a1 < a2 ? a2 - a1 + 1 : a1 - a2 + 1 ;
		  switch (cgr->type)
		    {
		    case 'N':/* intron */
		      a1 = cgr->a1 ; a2 = cgr->a2 ;
		      da = a2 - a1 - 1 ;
		      a1++ ; a2-- ;
		      cp = arrp (tct->dna,0 + keySet (tct->targetStart, trgt) ,char) ; /* last deleted base */
		      if (a1 - t1 > 0 && a2 - t1 < arrayMax (tct->dna) && a1 > 0 && a2 > 0)
			sprintf (buf, "Intron_%d:%d:%d:%c%c_%c%c", da, a1, a2,cp[a1-1],cp[a1],cp[a2-2],cp[a2-1]) ;
		      else
			sprintf (buf, "Intron_%d:%d:%d", da, a1, a2) ;
		      /* fall thru */
		    case 'D':
		      a1 = cgr->a1 ; a2 = cgr->a2 ;
		      x1 = cgr->x1 ; x2 = cgr->x2 ;
		      da = a2 - a1 - 1 ;
		      cp = arrp (tct->dna, cgr->a1-1+da - t1 + keySet (tct->targetStart, trgt) ,char) ; /* last deleted base */
		      cq = arrp (tct->dna, cgr->a1-1 - t1  + keySet (tct->targetStart, trgt) ,char) ;    /* base immediately left of the deletion */
		      ddw = 0 ;
		      /* slide right */
		      while (cp[1] && cp[1] == cq[1]) { a1++ ; a2++ ; x1++ ; x2++ ; cp++ ; cq++ ; ddw++ ; } 
		      /* slide left */
		      dda = 0 ;
		      while (a1 > 1 && *cp == *cq) { a1-- ; a2-- ; x1-- ; x2-- ; cp-- ; cq-- ; dda++ ; ddw-- ; } 
		      /* favor 75 */
		      if (a1 < 75 && a1 + dda >= 75)
			{ dda = 75 - a1 ; a1 += dda ; a2 += dda ; x1 += dda ; x2 += dda ; ddw += dda ; }
		      if (wig && ddw ) /* adjust the wiggle */
			{
			  array (wig, a1 + 1, int) -= mult ;
			  array (wig, a1 + 1 + da, int) += mult ;
			  array (wig, a1 + 1 - ddw, int) += mult ;
			  array (wig, a1 + 1 + da -ddw, int) -= mult ;
			}
			
		      if (da == 1)
			{
			  sprintf(tagBuf, "Del%c", ace_upper(dnaDecodeChar[(int)cq[1]])) ;
			  sprintf (buf, "%d:Del:%c%c:%c", a1, dnaDecodeChar[(int)cq[0]], ace_upper(dnaDecodeChar[(int)cq[1]]), dnaDecodeChar[(int)cq[0]]) ;
			}
		      else if (da <= 3) 
			{
			  sprintf (tagBuf, "Multi_deletion %d", da) ;
			  sprintf (buf, "%d:Del_%d:%c", a1, da, dnaDecodeChar[(int)cq[0]]) ;
			  k = strlen (buf) ; 
			  for (i = 1 ; i <= da ; i++)
			    buf[k++] = ace_upper (dnaDecodeChar[(int)cq[i]]) ;
			  buf[k++] = ':' ;
			  buf[k++] = dnaDecodeChar[(int)cq[0]] ;
			  buf[k++] = 0 ;
			}
		      else
			{
			  sprintf (tagBuf, "Multi_deletion %d", da) ;
			  sprintf (buf, "%d:Del_%d::", a1, da) ;
			}
		      break ;
		    case 'I':
		      a1 = cgr->a1 ; a2 = cgr->a2 ;
		      da = a1 ; a1 = a2 ; a2 = da ;
		      da = cgr->x2 - cgr->x1 + 1 ;
		      x1 = cgr->x1 ; x2 = cgr->x2 ;
		      cp = arrp (myDna, cgr->x1-2+da - t1,char) ;    /* last inserted base */
		      cq = arrp (tct->dna, cgr->a1-2 - t1  + keySet (tct->targetStart, trgt) ,char) ;    /* last matching base */
		      /* slide left */
		      while (x1 > da - 2 && *cp == *cq) { a1-- ; a2-- ; x1-- ; x2-- ; cp-- ; cq-- ; } 
		      if (da == 1)
			{
			  sprintf(tagBuf, "Ins%c", ace_upper (dnaDecodeChar[(int)cp[1]])) ;
			  sprintf (buf, "%d:Ins:%c:%c%c", a1, dnaDecodeChar[(int)cq[0]], dnaDecodeChar[(int)cq[0]], ace_upper (dnaDecodeChar[(int)cp[1]])) ;
			}
		      else  if (da < 30)
			{
			  sprintf (tagBuf, "Multi_insertion\t%d ", da) ;
			  sprintf (buf, "%d:Ins_%d:%c:%c%c", a1, da, dnaDecodeChar[(int)cq[0]], dnaDecodeChar[(int)cq[0]], ace_upper (dnaDecodeChar[(int)cp[1]])) ;
			  k = strlen (buf) ; cp -= da ;
			  for (i = 2 ; i <= da ; i++)
			    buf[k++] = ace_upper (dnaDecodeChar[(int)cp[i]]) ;
			  buf[k++] = 0 ;
			  k = strlen (tagBuf) ;
			  for (i = 1 ; i <= da ; i++)
			    tagBuf[k++] = ace_upper (dnaDecodeChar[(int)cp[i]]) ;
			  tagBuf[k++] = 0 ;
			}
		      else 
			{
			  sprintf (tagBuf, "Multi_insertion\t%d ", da) ;
			  sprintf (buf, "%d:Ins_%d:%c:    ", a1, da, dnaDecodeChar[(int)cq[0]]) ;
			}
		      break ;
		    }
		  uno = arrayp (lane->deUno, arrayMax (lane->deUno), UNO) ;
		  dictAdd (lane->deUnoDict, buf, &(uno->type)) ;
		  dictAdd (lane->deUnoDict, tagBuf, &(uno->tag)) ;
		  uno->target = trgt ;
		  uno->a1 = a1 ;
		  uno->a2 = a2 ;
		  uno->da = da ;
		  uno->nt += mult ;
		  uno->ns++ ;
		  if (cgr->type == 'I') /* insertion */
		    {
		      int x2 = cgr->x2 ; 
		      int x1 = cgr->x1 ;
		      int dx = x2 - x1 + 1 ;
		      if (myDna && dx > 0 && x1 > 0 && x2 < arrayMax (myDna))
			{
			  char buf[dx+1] ;
			  int i ;
			  for (i = 0 ; i < dx && i < INSMAX; i++) 
			    buf [i] = ace_upper(dnaDecodeChar[(int)arr(myDna, x1+i-1,char)]) ;
			  buf[i] = 0 ;
			  dictAdd (lane->deUnoDict, buf, &k) ;
			  uno->insert = k ;
			}
		    }
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
	    if (hit->dna && tct->dna &&
		cgr->x1 > 0 && cgr->a1 > 0 &&
		cgr->x1 < arrayMax (hit->dna) &&
		cgr->a1 - t1 > 0 &&
		cgr->a1 - t1  + keySet (tct->targetStart, trgt) < arrayMax (tct->dna)
		)
	      {
		char *cp = arrp (hit->dna, cgr->x1 - 1, char) ;
		char *cq = arrp (tct->dna, cgr->a1 - 1 - t1  + keySet (tct->targetStart, trgt) , char) ;
		int x ;
		int xMax = cgr->x2 - cgr->x1 + 1 ;
		int aMax = cgr->a2 - cgr->a1 + 1 ;
		if (aMax < xMax) xMax = aMax  ;
		for (x = 0 ; x < xMax ; x++, cp++, cq++, nnn++)
		  if (! (*cp & *cq))
		    {
		      UNO *uno ;
		      int k ;
		      if (tct->wantBrks)
			{
			  int k = keySet (brks, cgr->a1 + x - t1) ;
			  k += mult ; nBrks += mult ;
			  keySet (brks, cgr->a1 + x - t1) = k ;
			}
		      sprintf(buf, "%d:Sub:%c:%c", cgr->a1 + x, ace_upper(dnaDecodeChar[(int)*cq]), ace_upper(dnaDecodeChar[(int)*cp])) ;
		      dictAdd (lane->deUnoDict, buf, &k) ;
		      uno = arrayp (lane->deUno, arrayMax (lane->deUno), UNO) ;
		      uno->target = trgt ;
		      uno->a1 = cgr->a1 + x - 1 ;
		      uno->a2 = cgr->a1 + x + 1 ;
		      uno->da = 0 ;
		      uno->nt += mult ;
		      uno->ns++ ;
		      uno->type = k ;
		      sprintf(tagBuf, "%c2%c", ace_upper(dnaDecodeChar[(int)*cq]), ace_upper(dnaDecodeChar[(int)*cp])) ;
		      dictAdd (lane->deUnoDict, tagBuf, &(uno->tag)) ;
		    }
	      }
	    if (x1 < x2) hit->x1 = -hit->x1 ; /* hack isDown flag */
	    if (isRead1) hit->x2 = -hit->x2 ; /* hack isRead1 flag */
	  }
    }

  /********* fusions ******************/
  if (tct->geneFusion && nHit > 1)
    tctGetOneFusions (tct, lane) ;
  
  if (! tct->runTest)
    wego_log (hprintf (lane->h, "tctGetOneSamHit found %d hits, covering %ld bp and %ld brks %ld overhangs in %s\n", nHit, nBp, nBrks, nOverhangs, fNam)) ;
  ac_free (ai) ;
  ac_free (ao) ;
  ac_free (h) ;

  return TRUE ;
} /* tctGetOneSamHit */

/*************************************************************************************/

static void tctGetLaneHits (const void *vp)
{
  AC_HANDLE h = ac_new_handle () ;
  const TCT *tct = (const TCT *)vp ;
  int ii ;
  BOOL ok = TRUE ;
  
  while (channelGet (tct->getLaneInChan, &ii, int))
    {
      LANE *lane = arrp (tct->lanes, ii, LANE) ;

      if (ok) 
	{
	  if (tct->SAM)
	    ok = tctGetOneSamHit (tct, lane) ;
	  /*
	  else  if (tct->tabular)
	    pmbGetOneTabularHit (tct, lane) ;
	  */
	  else
	    ok = tctGetOneHit (tct, lane) ;
	}
	  
      channelPut (tct->getLaneOutChan, &ii, int) ;
    }

  ac_free (h) ;
  return ;
} /* tctGetLaneHits */

/*************************************************************************************/
 
static void tctFuseLaneHits (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;

  int ii, iiMax = arrayMax (tct->lanes) ;
  long int kk ;
  LANE *lane = 0 ;
  KEYSET brks = 0, covers ;
  BigArray brkIndels = 0 ;
  BigArray hits ;
  HIT *up, *vp ;

  brks = tct->brks = arrayHandleCreate (10000, KEY, tct->h) ;
  brkIndels = tct->brkIndels = bigArrayHandleCreate (10000, BKID, tct->h) ;
  covers = tct->covers = arrayHandleCreate (10000, KEY, tct->h) ;
 
  /* fuse the hits sort them and recreate the dict in the same order */
  kk = 0 ; 
  for (ii = 1 ; ii < iiMax && ii < 8 ; ii++)
    {
      lane = arrp (tct->lanes, ii, LANE) ;
      kk +=  arrayMax (lane->hits) ;
    }
  hits = tct->hits = bigArrayHandleCreate (kk, HIT, tct->h) ; 
  if (kk < 256) kk = 256 ;
  vp = bigArrayp (tct->hits, kk - 1, HIT) ; /* make room */

  for (kk = 0, ii = 1 ; ii < iiMax ; ii++)
    {
      int iMax ;
      lane = arrp (tct->lanes, ii, LANE) ;
      iMax = arrayMax (lane->hits) ;
      if (iMax)
	{
	  up = arrp (lane->hits, 0, HIT) ;
	  vp = bigArrayp (tct->hits, kk + iMax - 1, HIT) ; /* make room */
	  vp = bigArrp (hits, kk, HIT) ; 
	  memcpy (vp, up, iMax * sizeof(HIT)) ;
	  kk += iMax ; 
	}
      arrayDestroy (lane->hits) ;
    }
  bigArrayMax (hits) = kk ;
  bigArraySort (hits, tctA1Order) ;
  /* add in the hit a marker of the largest a2 seen so far , allowing later a dicothomy search */
  if (bigArrayMax (hits))
    {
      long int i, iMax = bigArrayMax (hits) ;
      int b2 = 0 ;

      for (i = 0, up = bigArrp (hits, 0, HIT) ; i < iMax ; i++, up++)
	{
	  if(up->a2 > b2) b2 = up->a2 ;
	  if(up->a1 > b2) b2 = up->a1 ;
	  up->b2 = b2 ;
	}
    }    

  /* fuse the coverage plots */
  for (ii = 1 ; ii < iiMax ; ii++)
    {
      KEY *kp ;
      int i, iMax ;
      lane = arrp (tct->lanes, ii, LANE) ;
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
      lane = arrp (tct->lanes, ii, LANE) ;
      if (lane->brks)
	{
	  iMax = arrayMax (lane->brks) ;
	  if (iMax)
	    {
	      keySet (brks, iMax) += 0 ;
	      for (i = 0, kp = arrayp (lane->brks, i, KEY) ; i < iMax ; i++, kp++)
	    keySet (brks, i) += *kp ;
	    }
	}
    } 

  /* fuse the brkIndels plots */
  for (ii = 1 ; ii < iiMax ; ii++)
    {
      BKID *b ;
      long int i, j, iMax ;
      lane = arrp (tct->lanes, ii, LANE) ;
      if (lane->brkIndels)
	{
	  iMax = arrayMax (lane->brkIndels) ;
	  j =  bigArrayMax (brkIndels) ;
	  if (iMax)
	    {
	      bigArray (brkIndels, j + iMax - 1, BKID).pos = 0 ;
	      for (i = 0, b = arrayp (lane->brkIndels, i, BKID) ; i < iMax ; i++, b++)
		{
		  bigArray (brkIndels, j++, BKID) = *b ;
		  if (0 && b->c1)
		    fprintf (stderr, "=== fuse lane %d, iMax=%ld b[%ld].da=%d    brkIndels=%p max=%ld c1/c2 = %d/%d\n", ii, iMax, i, b->daa, brkIndels, bigArrayMax (brkIndels), b->c1, b->c2) ;
		}
	    }
	}
    }
  
  /* fuse the geneFusions */
  if (tct->geneFusion)
    {
      KEYSET fusionHisto = tct->fusionHisto = keySetHandleCreate (tct->h) ;
      Array geneFusions = tct->geneFusions = arrayHandleCreate (10000, GF, tct->h) ;
      KEYSET geneSupport = tct->geneSupport = keySetHandleCreate (tct->h) ;
      
      for (ii = 1 ; ii < iiMax ; ii++)
	{
	  GF *a, *b ;
	  int i, j, iMax ;
	  lane = arrp (tct->lanes, ii, LANE) ;
	  if (lane->geneFusions)
	    {
	      KEYSET ks = lane->geneSupport ;
	      iMax = arrayMax (lane->geneFusions) ;
	      j =  arrayMax (geneFusions) ;
	      if (iMax)
		{
		  array (geneFusions, j + iMax - 1, GF).g1 = 0 ;
		  for (i = 0, b = arrayp (lane->geneFusions, i, GF) ; i < iMax ; i++, b++)
		    {
		      a = arrayp (geneFusions, j++, GF) ;
		      *a = *b ;
		    }
		}
	      j = keySetMax (lane->fusionHisto) ;
	      for (i = 0 ; i < j ; i++)
		keySet (fusionHisto, i) += keySet (lane->fusionHisto, i) ;
	      if (ks)
		for (j = 0 ; j < keySetMax (ks) ; j++)
		  keySet (geneSupport, j) += keySet (ks, j) ;
	    }
	}
    }

  /* fuse the geneSupport */
  if (tct->geneFusion)
    {
      KEYSET geneSupport = tct->geneSupport = keySetHandleCreate (tct->h) ;
      
      for (ii = 1 ; ii < iiMax ; ii++)
	{	  
	  int i, iMax ;
	  lane = arrp (tct->lanes, ii, LANE) ;
	  if (lane->geneSupport)
	    {
	      iMax = arrayMax (lane->geneSupport) ;
	      for (i = iMax - 1 ; i>- 0 ; i--)
		keySet (geneSupport, i) += keySet (lane->geneSupport, i) ;
	    }
	}
    }

  /* fuse the geneFusion coordinates, now because we cannot reallocate an array in parallel mode */
  if (lane->genes)
    {
      int gene, gMax = arrayMax (lane->genes) ;
      GENE *ga, *gb ;
      
      if (gMax)
	{
	  for (gene = 1, ga = arrp (lane->genes, gene, GENE); gene < gMax ; ga++, gene++)
	    if (ga->chrom)
	      {
		gb = arrayp (tct->genes, gene, GENE) ;
		*gb = *ga ;
	      }
	}
    }

  if (1)
    {
      KEY *kp ;
      int i, n, iMax = keySetMax (brks) ;
      
      for (i = n = 0, kp = arrayp (brks, i, KEY) ; i < iMax ; i++, kp++)
	n += *kp ;
     
      wego_log (hprintf (lane->h, "Selected %d brks\n", n)) ;
    }
  ac_free (h) ;
} /* tctFuseLaneHits */

/*************************************************************************************/

static void tctGetHits (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;

  int ii, iMax = arrayMax (tct->lanes) ;

  if (tct->SAM_file || tct->hitFile)
    iMax = 2 ;
  if (tct->maxLanes > 0 && tct->maxLanes < iMax) iMax = tct->maxLanes + 1 ;

  tct->getLaneInChan = channelCreate (iMax, int, h) ;
  tct->getLaneOutChan = channelCreate (iMax, int, h) ;

  for (ii = 1 ; ii < iMax ; ii++)
    channelPut (tct->getLaneInChan, &ii, int) ;

  channelClose (tct->getLaneInChan) ;

  for (ii = 0 ; ii < tct->nAna && ii < iMax - 1 ; ii++)
    wego_go (tctGetLaneHits, tct, TCT) ;      // launch several threads
  for (ii = 1 ; ii < iMax ; ii++)
    {
      int x = 0 ;
      channelGet (tct->getLaneOutChan, &x, int) ;
      if (tct->deUno)
	tctDeUnoFuseOneLane (tct, x) ;
    }
  if (1) /* tct->wantBrks) */
    tctFuseLaneHits (tct) ;
  if (0) exit(0) ;

   ac_free (h) ;  
} /* tctGetHits  */

/*************************************************************************************/
/* we expect to read the file 
 * tmp/METADATA/$target.mrna_map_ln_gc_gene_geneid.txt
 * col 2: chrom:123-456
 * col 5: gene  one line per exon, we must cluster the coordinates 
 */

static void tctGetGeneMap (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (tct->geneMapFileName, 0, h) ;

  if (! tct->geneDict)
    tct->geneDict = dictHandleCreate (128, tct->h) ;
  tct->genes = arrayHandleCreate (20000, GENE, tct->h) ;
  if (! ai)
    {
      fprintf (stderr, "// FATAL ERROR : cannot open geneMapFile file %s\n// Try TCT -h\n", tct->geneMapFileName) ;
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
      dictAdd (tct->geneDict, cp, &chrom) ;
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
      dictAdd (tct->geneDict, cp, &gene) ;
      up = arrayp (tct->genes, gene, GENE) ;
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
} /*  tctGetGenemap  */

/*************************************************************************************/

static void tctGetLaneList (TCT *tct)
{
  AC_HANDLE h = 0 ;
  ACEIN ai = aceInCreate (tct->laneListFileName, 0, h) ;
  ACEOUT ao = 0 ;
  int nn = 0 ;
  char *fNam ;

  tct->laneDict = dictHandleCreate (128, tct->h) ;
  tct->lanes = arrayHandleCreate (128, LANE, tct->h) ; 

  if (tct->fastcFile)
    {
      LANE *lane ;
      const char *ccp = "fastcFile" ;

      dictAdd (tct->laneDict, ccp, &nn) ;
      lane = arrayp (tct->lanes, nn, LANE) ;
      lane->lane = nn ;
      lane->run = tct->run ;
      lane->h = ac_new_handle () ;
      	      
      lane->cloneDict = dictHandleCreate (10000, lane->h) ; 
      lane->qcDict = dictHandleCreate (1000, lane->h) ; 
      lane->clones = arrayHandleCreate (10000, CLONE, lane->h) ;
      if (tct->wantBrks)
	{
	  lane->brks = arrayHandleCreate (10000, KEY, lane->h) ;
	  lane->brkIndels = arrayHandleCreate (10000, BKID, lane->h) ;
	}
      lane->covers = arrayHandleCreate (10000, KEY, lane->h) ;
      lane->aliLn = arrayHandleCreate (10000, KEY, lane->h) ;
      lane->hits = arrayHandleCreate (10000, HIT, lane->h) ;
      if (tct->wantDeUno)
	{
	  lane->deUno = arrayHandleCreate (10000, UNO, lane->h) ;
	  lane->wig =arrayHandleCreate (1000, int, lane->h) ;
	  lane->deUnoDict = dictHandleCreate (1000, lane->h) ;
	}
      return ;
    }
  h = ac_new_handle () ;
  if (tct->laneListFileName)
    ai = aceInCreate (tct->laneListFileName, 0, h) ;
  else if (tct->SAM)
    ai = aceInCreateFromText ("sam", 0, h) ;
  else if (tct->hitFile)
    ai = aceInCreateFromText ("hits", 0, h) ;

  if (! ai)
    {
      fprintf (stderr, "// FATAL ERROR : cannot open LaneList file %s\n// Try TCT -h\n", tct->laneListFileName) ;
      if (0) exit (1) ;
    }
  if (tct->makeTest)
     ao = aceOutCreate (tct->makeTest, ".LaneList", FALSE, h) ;
 
  fNam =  strnew (aceInFileName (ai), h) ;
  while (aceInCard (ai))
    {
      char * cp ;
      LANE *lane ;

      cp = aceInWord (ai) ;
      if (cp && *cp)
	{
	  dictAdd (tct->laneDict, cp, &nn) ;
	  lane = arrayp (tct->lanes, nn, LANE) ;
	  lane->lane = nn ;
	  if (ao)
	    aceOutf (ao, "%s.%d\n", tct->makeTest, nn) ;

	  lane->run = tct->run ;
	  lane->h = ac_new_handle () ;
	  
	  if (!tct->snpCount)
	    {
	      lane->cloneDict = dictHandleCreate (10000, lane->h) ;
	      lane->qcDict = dictHandleCreate (1000, lane->h) ;
	      lane->clones = arrayHandleCreate (10000, CLONE, lane->h) ;
	      if (tct->wantBrks)
		{
		  lane->brks = arrayHandleCreate (10000, KEY, lane->h) ;
		}
	      lane->brkIndels = arrayHandleCreate (10000, BKID, lane->h) ;
	      lane->covers = arrayHandleCreate (10000, KEY, lane->h) ;
	      lane->aliLn = arrayHandleCreate (10000, KEY, lane->h) ;
	      lane->hits = arrayHandleCreate (10000, HIT, lane->h) ;
	      if (tct->wantDeUno)
		{
		  lane->deUnoDict = dictHandleCreate (1000, lane->h) ;
		  lane->deUno = arrayHandleCreate (10000, UNO, lane->h) ;
		  lane->wig = arrayHandleCreate (1000, int, lane->h) ;
		}
	    }
	}
      if (tct->maxLanes && nn >= tct->maxLanes)
	break ;
    }
  
  if (! tct->debug && ai)
    fprintf (stderr, "// Found %d lanes in file %s\n"
	     , nn
	     , fNam ? fNam : "0"
	     ) ;

  ac_free (ao) ;
  ac_free (ai) ;
  ac_free (h) ;
  return ;
} /* tctGetLaneList */

/*************************************************************************************/

static void tctGetBrks (TCT *tct)
{
  AC_HANDLE h = 0 ;
  KEY *kp, *wp ;
  KEYSET brks, brks2 ;
  KEYSET covers = tct->covers ;
  float minSnpFrequency = tct->minSnpFrequency ;
  int minSnpCount = tct->minSnpCount ;
  int minSnpCover = tct->minSnpCover ;
  int microDelta = tct->microDelta ;
  int dm = 2 * microDelta + 1;
  int ii, jj, iMax ;
  BOOL debug = FALSE ;

  if (!tct->wantBrks)
    return ;
  h = ac_new_handle () ;
  brks = tct->brks ;
  brks2 = tct->brks = arrayHandleCreate (10000, KEY, tct->h) ;
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
	  for (i = - tct->microDelta ; i <= tct->microDelta ; i++)
	    m += *(kp + i) ;
	  if (m >= minSnpCount &&
	      200 * m >= dm * minSnpFrequency * (*wp) /* 200 because these are not the final counts */
	      )
	    keySet (brks2, jj++) = ii ;
	} 
    }
  if (tct->externalBrks)
    {
      KEYSET ks = tct->externalBrks ;
      for (ii = 0 ; ii < keySetMax (ks) ; ii++)
	keySet (brks2, jj++) = keySet (ks, ii) ;
    }
  keySetSort (brks2) ;
  keySetCompress (brks2) ;

  if (debug)
    {
      fprintf (stderr, "...tctGetBrks") ; 
      for (ii = 0 ; ii < keySetMax (brks2) ; ii++)
        fprintf (stderr, " %d", keySet (brks2, ii)) ;
      fprintf (stderr, "\n") ;
    }

  ac_free (brks) ;
  wego_log (hprintf (tct->h, "Selected %d well covered brks\n", arrayMax (tct->brks))) ;
  ac_free (h) ;
} /* tctGetBrks  */

/*************************************************************************************/
/*     format: tab delimited file (filtered on minSupport
 *	 mrna1 a1 a2 mrna2 b1 b2 support
 *	   one fusion per line, jump for position a2 of mrna1 to position b1 of mrna2
 */	

static int tctGFCumul (Array aa)
{
  GF *up, *vp ;
  int ii, jj, iiMax = arrayMax (aa) ;
  
  for (ii = 1, jj = 0, vp = arrp (aa, 0, GF), up = vp + 1 ; ii < iiMax ; up++, ii++)
    {
      if (up->m1 != vp->m1 || 
	  up->m2 != vp->m2 ||
	  up->a2 != vp->a2 ||
	  up->b1 != vp->b1
	  )
	{ 
	  jj++ ;
	  vp++ ;
	  if (vp < up)
	    *vp = *up ;
	}
      else
	{
	  vp->n0 += up->n0 ;
	  vp->seq += up->seq ;
	  vp->tag += up->tag ;
	}

    }
  arrayMax (aa) = jj ;
  return 0 ;
}

static int tctGFOrder (const void *va, const void *vb)
{
  const GF* up = (const GF*) va ;
  const GF* vp = (const GF*) vb ;
  int n ;

  n = up->g1 - vp->g1 ; if (n) return n ;
  n = up->m1 - vp->m1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  n = up->g2 - vp->g2 ; if (n) return n ;
  n = up->m2 - vp->m2 ; if (n) return n ;
  n = up->b1 - vp->b1; if (n) return n ;

  return 0 ;
}

static void tctGetGeneFusionList (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (tct->geneFusionFileName, 0, h) ;
  Array geneFusions = tct->geneFusions ;
  DICT *targetDict = tct->targetDict ;
  /*   DICT *gfDict = dictHandleCreate (10000, h) ; */
  GF *up ;
  int jj = 0 ;
  int g12, da, chrom1, chrom2, run, a1, a2, b1, b2, g1, g2, m1, m2, x1, x2, y1, y2, seq, tag ;
  int minGF = tct->minGF ;
  int dx = tct->minOverhang ;

  if (dx < 15)
    dx = tct->minOverhang = 15 ;
  if (!ai)
    {
      fprintf (stderr, "// FATAL ERROR: cannot open -geneFusion %s\n"
	       , tct->geneFusionFileName) ;
      exit (1) ;
    }

  if (! geneFusions)
    geneFusions = tct->geneFusions = arrayHandleCreate (10000, GF, tct->h) ;
  while (aceInCard (ai))
    {
      char *cp = aceInWord (ai) ;
      if (! cp || *cp == '#') continue ;
      dictAdd (targetDict, cp, &run) ;
      aceInStep (ai, '\t') ; aceInWord (ai) ; /* drop read name */
      aceInStep (ai, '\t') ; cp = aceInWord (ai) ;
      if (!cp)	continue ;
      dictAdd (targetDict, cp, &g12) ;
      aceInStep (ai, '\t') ; cp = aceInWord (ai) ;
      if (!cp)	continue ;
      dictAdd (targetDict, cp, &g1) ;
      aceInStep (ai, '\t') ; cp = aceInWord (ai) ;
      if (!cp)	continue ;
      dictAdd (targetDict, cp, &g2) ;
      aceInStep (ai, '\t') ; cp = aceInWord (ai) ; /* chrom1 */
      dictAdd (targetDict, cp, &chrom1) ;
      aceInStep (ai, '\t') ; cp = aceInWord (ai) ; /* chrom2 */
      if (!cp)	continue ;
      dictAdd (targetDict, cp, &chrom2) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &da) ; /* distance */
      aceInStep (ai, '\t') ; aceInWord (ai) ; /* drop READ */

      aceInStep (ai, '\t') ; if (! aceInInt (ai, &x1)) continue ;
      aceInStep (ai, '\t') ; if (! aceInInt (ai, &x2)) continue ;
      aceInStep (ai, '\t') ; cp = aceInWord (ai) ; /* mRNA1 */
      if (!cp)	continue ;
      dictAdd (targetDict, cp, &m1) ;
      aceInStep (ai, '\t') ; if (! aceInInt (ai, &a1)) continue ;
      aceInStep (ai, '\t') ; if (! aceInInt (ai, &a2)) continue ;

      aceInStep (ai, '\t') ; if (! aceInInt (ai, &y1)) continue ;
      aceInStep (ai, '\t') ; if (! aceInInt (ai, &y2)) continue ;
      aceInStep (ai, '\t') ; cp = aceInWord (ai) ; /* mRNA2 */
      if (!cp)	continue ;
      dictAdd (targetDict, cp, &m2) ;
      aceInStep (ai, '\t') ; if (! aceInInt (ai, &b1)) continue ;
      aceInStep (ai, '\t') ; if (! aceInInt (ai, &b2)) continue ;

      aceInStep (ai, '\t') ; if (! aceInInt (ai, &seq)) continue ;
      aceInStep (ai, '\t') ; if (! aceInInt (ai, &tag)) continue ;
      aceInStep (ai, '\t') ; if (tag < minGF) continue ;

      up = arrayp (geneFusions, jj++, GF) ;
      up->run = tct->run ? 0 : run ;
      up->chrom1 = chrom1 ; up->chrom2 = chrom2 ; up->da = da ;
      up->g12 = g12 ;
      up->g1 = g1 ; up->m1 = m1 ; up->a1 = a1 ; up->a2 = a2 ; up->da = da ; 
      up->dx = (x1 < x2) ? y1 - x2 - 1 : x2 - y1 - 1 ;
      up->g2 = g2 ; up->m2 = m2 ; up->b1 = b1 ; up->b2 = b2 ;
      up->n0 = tag ; /* number of supports in external file */
    }


  if (jj > 1)
    {
      arraySort (geneFusions, tctGFOrder) ;
      tctGFCumul (geneFusions) ;
    }
  ac_free (h) ;
}  /* tctGetGeneFusionList */

/*************************************************************************************/

static void tctExportGeneFusion2tsf (TCT *tct)
{ 
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (tct->outFileName, ".tsf", tct->gzo, h) ;
  Array geneFusions = tct->geneFusions ;
  DICT *targetDict =tct->targetDict ;
  int ii, iiMax = arrayMax (geneFusions) ;
  GF *gf ;

  aceOutf (ao, "# Gene\tRun\ti3ti3t3i\tSupport\tChrom1\tg1\tm1\ta2\tChrom2\tg2\tm2\tb1\tda\tdx\n") ;
  for (ii = 0, gf = arrp (geneFusions, 0, GF) ; ii < iiMax ; ii++, gf++)
    {
      aceOutf (ao, "%s\t%s\ti3ti3t3i\t%d",  dictName (targetDict, gf->g12), gf->run ? dictName (targetDict, gf->run) : tct->run, gf->n0) ;
      aceOutf (ao, "\tch_%s\t%s\t%s\t%d", dictName (targetDict, gf->chrom1), dictName (targetDict, gf->g1), dictName (targetDict, gf->m1), gf->a2) ;
      aceOutf (ao, "\tch_%s\t%s\t%s\t%d", dictName (targetDict, gf->chrom2), dictName (targetDict, gf->g2), dictName (targetDict, gf->m2), gf->b1) ;
      aceOutf (ao, "\t%d\t%d", gf->da, gf->dx) ;
      aceOut (ao, "\n") ;
    }
  ac_free (h) ;
} /* tctExportGeneFusion2tsf */

/*************************************************************************************/
/*     several formats are recognized
 *	 chr7:1234     
 *       chr7:1234:
 *       chr7  1234
 *	   one position per line\n"
 */	
static void tctGetExternalSnpList (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (tct->externalSnpListFileName, 0, h) ;
  KEYSET ks ;
  BitSet ebs ;
  int jj = 0 ;
  int t1 = tct->t1 ? tct->t1 : tct->offset ; 
  DICT *deUnoDict = tct->deUnoDict ;

  if (!ai)
    {
      fprintf (stderr, "// FATAL ERROR: cannot open -externalSnpList %s\n"
	       , tct->externalSnpListFileName) ;
      exit (1) ;
    }

  ks = tct->externalBrks = keySetHandleCreate (tct->h) ;
  ebs = tct->externalBitSet = bitSetCreate (12000000, tct->h) ;
  while (aceInCard (ai))
    {
      char buf[512] ;
      char tagBuf[512] ;

      char cutter, *cq, *cp ;
      int x = 0, trgt = 0, tag = 0, dx = 0, ddx = 0, type = 0 ;
      
      cp = aceInWordCut (ai, "\t:", &cutter) ;
      if (! cp || ! *cp || *cp == '#')
	continue ;
 
     cq = strchr (cp, ':') ;
     if (cq)
       *cq++ = 0 ;
     dictAdd (tct->targetDict, cp, &trgt) ;
     if (cq)
       cp = cq + 1 ;
     else
       {
	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;
       }       
     strncpy (buf, cp, 510) ;
     buf[510] = 0 ;
     if (! *cp || sscanf (cp, "%d:", &x) != 1)
       continue ;

     if (! tct->t2 || (tct->t1 <= x && tct->t2 >= x))
       {
	 keySet (ks, jj++) = x - t1  + 1 ;
	 if (tct->debug) fprintf (stderr, "jj = %d x = %d\n", jj, x) ;
	 bitSet (ebs, x - tct->t1  + 1) ;
       }
      else
	continue ;
     dictAdd (deUnoDict, buf, &tag) ;
     
     cp = strchr (cp, ':') ;
     if (! cp)
       continue ;
     cp++ ;
     ddx = dx = 1 ;
     if (cp[3]=='_')
       {
	 if (sscanf (cp+4, "%d:", &dx) != 1)
	   continue ;
       }
     if (! strncmp (cp, "Sub_", 4))
       {
	 sprintf (tagBuf, "Multi_sub\t%d", dx) ;
	 ddx = 0 ;
       }
     else if (! strncmp (cp, "Sub", 3))
       {
	 sprintf (tagBuf, "%c2%c\t0\t-", cp[5],cp[7]) ;
	 ddx = 0 ;
       }
     else if (! strncmp (cp, "Del_", 4))
       {
	 sprintf (tagBuf, "Multi_deletion \t%d", dx) ;
	 ddx = -dx ;
       }
     else if (! strncmp (cp, "Del", 3))
       {
	 sprintf (tagBuf, "Del%c\t1\t%c", cp[6],cp[6]) ;
	 ddx = -1 ;
       }
    else if (! strncmp (cp, "Ins_", 4))
       {
	 sprintf (tagBuf, "Multi_insertion\t%d", dx) ;
	 ddx = dx ; dx = 0 ;
       }
     else if (! strncmp (cp, "Ins", 3))
       {
	 sprintf (tagBuf, "Ins%c\t1\t%c", cp[8],cp[8]) ;
	 ddx = 1 ;
	 dx = 0 ;
       }
     else
       continue ;
     
     dictAdd (deUnoDict, buf, &type) ;
     dictAdd (deUnoDict, tagBuf, &tag) ;

     UNO *uno = bigArrayp (tct->deUno, bigArrayMax (tct->deUno), UNO) ;
     uno->isExternal = TRUE ;
     uno->target = trgt ;
     uno->type = type ;
     uno->a1 = x ;
     uno->a2 = x + dx + 1 ;
     uno->ddx = ddx ;
     uno->tag = tag ;
    }

  ac_free (ai) ;
  ac_free (h) ;
} /* tctGetExternalSnpList */

/*************************************************************************************/
/* using the genome, snail trail the indels */ 
static int tctSnailTrailSlippingIndels (const TCT *tct)
{
  BigArray brkIndels = tct->brkIndels  ;
  KEYSET brks = tct->brks ;
  KEYSET covers = tct->covers ;
  KEY *wp ;
  int ii, jj, kk ;
  float minSnpFrequency = tct->minSnpFrequency ;
  int minSnpCount = tct->minSnpCount ;
  int minSnpCover = tct->minSnpCover ;
  long int iMax = bigArrayMax (brkIndels) ; 
  int kMax = arrayMax (brks) ; 
  int kNew = arrayMax (brks) ; 
  int Delta = tct->Delta ;
  int pos, dx ;
  Array dna = tct->dna ;
  char *cp, buf[256] ;
  BOOL debug = FALSE ;

  bigArraySort (brkIndels, tctBrkidOrder) ;
  for (ii = jj = kk = 0  ;  ii < iMax ; ii ++)
    {
      BKID *b2, *b1 = bigArrp (brkIndels, ii, BKID) ;
      if (! b1->mult) 
	continue ;
      for (jj = ii+1  ;  jj < iMax ; jj ++)
	{
	  b2 = bigArrp (brkIndels, jj, BKID) ;
	  if (b1->pos == b2->pos &&
	      b1->daa == b2->daa &&
	      b1->dxx == b2->dxx
	      )
	    {
	      b1->mult += b2->mult ;
	      b2->mult = 0 ;
	    }
	}
      if (kk < ii)
	{
	  b2 = bigArrp (brkIndels, kk, BKID) ;
	  *b2 = *b1 ; 
	  b1->mult = 0 ;
	}
      kk++ ;      
    }
  iMax = bigArrayMax (brkIndels) = kk ;
  for (ii = jj = 0  ;  ii < iMax ; ii ++)
    {
      BKID *b = bigArrp (brkIndels, ii, BKID) ;
      BOOL found = FALSE ;
      int mult = b->mult ;
      pos = b->pos ;
  
      wp = arrayp (covers, pos, KEY) ;
      if ( *wp < minSnpCover ||
	   mult < minSnpCount ||
	  500 * mult <  minSnpFrequency * (*wp) /* 300 because these are not the final counts */
	  )
	continue ;
	
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
	  BKID *b1 = bigArrp (brkIndels, jj, BKID) ;
	  BKID *b2 = bigArrp (brkIndels, kk, BKID) ;
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
  iMax = bigArrayMax (brkIndels) = jj ;
  if (0)
    tctShowBrkIndels (tct->brkIndels) ;
      
  for (ii = 0  ;  ii < iMax ; ii ++)
    {
      BKID *b = bigArrp (brkIndels, ii, BKID) ;
      pos = b->pos ;
      for (jj = ii + 1 ; jj < iMax ; jj++)
	{
	  BKID *b1 = bigArrp (brkIndels, jj, BKID) ;

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
      BKID *b = bigArrp (brkIndels, ii, BKID) ;
      if (b->mult)
	{
	  BKID *b1 = bigArrp (brkIndels, jj, BKID) ;
	  if (jj < ii)
	    *b1 = *b ;
	  jj++ ;
	}
    }
  iMax = bigArrayMax (brkIndels) = jj ;

  for (ii = jj = 0  ;  ii < iMax ; ii ++)
    {
      BKID *b = bigArrp (brkIndels, ii, BKID) ;
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
  
  if (debug)
    tctShowBrkIndels (tct->brkIndels) ;

  return arrayMax (brks) ;
} /* tctSnailTrailSlippingIndels */
  
/*************************************************************************************/

static void tctGetSegs (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  int Delta = tct->Delta ;
  KEY *kp ;
  SEG *seg ;
  KEYSET brks = tct->brks  ;
  KEYSET nnn = keySetHandleCreate (h) ;
  KEYSET lns = keySetHandleCreate (h) ;
  Array segs ;
  BigArray brkIndels = tct->brkIndels  ;
  int ii, jj = 0, iMax = arrayMax (brks) ;
  int n, a, a1 = 0, a2 = 0 ;
  int dx = tct->Delta ;
  BOOL debug = FALSE ;

  segs = tct->segs = arrayHandleCreate (10000, SEG, tct->h) ;
  if (1) tctSnailTrailSlippingIndels (tct) ;
  
  if (dx > 8)
    dx = 8 ;

   /* create a SEG where the brks are quasi-continuous IDEA we could limit to positions where the ratio events/coverage is high enough */
  if (iMax)
    for (ii = jj = n = a1 = a2 = 0, kp = arrayp (brks, ii, KEY) ;  ii < iMax + 1 ; ii++, kp++)
      {
	a = *kp ; /* do not add tct->t1, the seg is in local coords */
	if (!a && ii < iMax) continue ;
	if (ii == iMax || a - a2 > Delta)
	  {
	    if (a1)
	      {
		seg = arrayp (segs, jj++, SEG) ;
		seg->a1 = a1  ;
		seg->a2 = a2  ;
		seg->b1 = seg->a1 - dx ; 
		if (seg->b1 < 1)
		  seg->b1 = 1 ;
		seg->b2 = seg->a2 + dx ;

		keySet (nnn, n)++ ;
		keySet (lns, a2 - a1 + 1)++ ;
	      }
	    n = 0 ;
	}
	if (! n++) a1 = a ;
	a2 = a ;
      }
  
  if (debug)    tctShowSegs (segs, "==A") ;

  if (brkIndels && bigArrayMax (brkIndels))
    {
      long int ib, ibMax = bigArrayMax (brkIndels) ;

      for (ib = 0 ; ib < ibMax ; ib++)
	{
	  BKID *b = bigArrp (brkIndels, ib, BKID) ;  

	  if (debug && b->c1)
	    fprintf (stderr, "=== getSegs brkIndels ibMax=%ld b[%ld].da=%d\n", ibMax, ib, b->daa) ;
	   
	  if (b->daa > Delta)
	    { /* create a new seg */
	      a1 = b->pos ;
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
	      seg->b1 = seg->a1 - dx ; 
	      if (seg->b1 < 1)
		seg->b1 = 1 ;
	      seg->b2 = seg->a2 + dx ;
	      if (b->dxx > 0) 
		seg->b2 += b->dxx ;
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
		      if (seg->b1 > seg->a1 - tct->Delta - b->dxx)
			{
			  seg->b1 = seg->a1 - tct->Delta - b->dxx ;
			  seg->b2 = seg->a2 + tct->Delta  ;
			}
		    }
		}
	    }
	}
    }
 

  if (debug)    tctShowSegs (segs, "==B") ;

  iMax = arrayMax (segs) ;
  if (1)
    {
      int dx = 8 ; 
      int t2 = arrayMax (tct->dna) ;
      for  (ii = 0 ; ii < arrayMax (segs) ; ii++)
	{
	  seg = arrp (segs, ii, SEG) ;
	  seg->b1 = seg->c1 = seg->a1 - dx ; if (seg->b1 < 1) seg->b1 = seg->c1 = 1 ; 
	  seg->b2 = seg->c2 = seg->a2 + dx ; if (seg->b2 >= t2) seg->b2 = seg->c2 = t2 - 1 ;
	}
      arraySort (segs, tctSegOrder) ;
      arrayCompress (segs) ;
    }

 /* elongate segs that endup in sliding mess */
  iMax = arrayMax (segs) ;
  for  (ii = 0 ; ii < iMax ; ii++)
    {
      const register char *cp ;
      int b1, b2, i, dx ;
      Array dna = tct->dna ;
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
	    b2 = seg->b1 + k ;
	}
      if (b2 > seg->b2)
	seg->b2 = b2 ;
    }	  

  if (debug)    tctShowSegs (segs, "**A") ;


  if (debug)    tctShowSegs (segs, "**B") ;
  if (0)
    {
      tctCompressSegs (tct) ;
      tctShowSegs (segs, "**C") ;
    }
  if (0)
    {
      for  (ii = 0 ; ii < arrayMax (segs) ; ii++)
	{
	  dx = 8 ;
	  seg = arrp (segs, ii, SEG) ;
	  if (seg->c1 && seg->b1 > seg->c1 - dx)
	    seg->b1 = seg->c1 - dx ; 
	  if (seg->b1 < 1)
	    seg->b1 = 1 ;
	  if (seg->c2 && seg->b2 < seg->c2 + dx)
	    seg->b2 = seg->c2 + dx ; 
	}
      tctShowSegs (segs, "**D") ;
    }

  iMax = arrayMax (segs) ;
  wego_log (hprintf (tct->h, "Selected %d variable segments\n", iMax)) ;
 
  ac_free (h) ;
  return ;
} /* tctGetSegs  */

/*************************************************************************************/
/*********************** SNP validation **********************************************/
/*************************************************************************************/

static void tctValGetSnpList (TCT *tct, BOOL extend)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = 0 ;
  const char *ccp ;
  int i, j, dj, ii = 0, nok, wLn = tct->valWordLength ;
  int extendLn = tct->snpExtend ? tct->snpExtendLength : 0 ;
  VW *up ;
  int shiftR = 2 * (wLn - 1) ;
  DICT *words, *names, *zones ;
  unsigned long long int oligo = 0, oligoR = 0, mask = 1, gtMask = 1 ;
  Associator ass ;

  if (tct->valSnpList) 
    {
      char *cp = strnew (tct->valSnpList, h), *cq = cp ;

      while ((cq = strchr (cq, ',')))
	*cq = '\n' ;
      ai = aceInCreateFromText (cp, 0, h) ;
    }
  else if (tct->valSnpListFileName)
    {
      ai = aceInCreate (tct->valSnpListFileName, 0, h) ;
      if (! ai)
	{
	  fprintf (stderr, "Could not open -wordFile %s\n", tct->valSnpListFileName) ;
	  exit (1) ;
	}
    }
  if (! ai)
    {
      fprintf (stderr, "in -count mode, please specify -words or -wordFile (try -help)\n") ;
      exit (1) ;
    }

  mask = (mask << 2*wLn) ; mask-- ;
  { int i = wLn/2 ; gtMask = (gtMask << 2*i) ;} 

  tct->valSnps = arrayHandleCreate (1000, VW, tct->h) ;
  tct->valSnpZones = zones = dictHandleCreate (1000, tct->h) ;
  tct->valSnpNames = names = dictHandleCreate (1000, tct->h) ;
  tct->valSnpWords = words = dictHandleCreate (1000, tct->h) ;
  ass = tct->valAss = assBigCreate (1000) ;

  dictAdd (tct->valSnpZones , "-", 0) ; /* default zone */
  ii = -1 ;
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      char cutter ;

      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#')
	continue ;
      i = strlen (ccp) ;
      if (i < wLn)
	continue ;

      ii++ ;
      up = arrayp (tct->valSnps, ii, VW) ;
      dictAdd (words, ccp, &(up->motif)) ;

      aceInStep (ai, '\t') ;
      ccp = aceInWordCut (ai, "\t", &cutter) ;
      if (! ccp)
	ccp = hprintf (h, "w%d\n", ii + 1) ;
      dictAdd (names, ccp, &(up->name)) ;

      aceInStep (ai, '\t') ;
      ccp = aceInWordCut (ai, "\t", &cutter) ;
      if (ccp)
	dictAdd (zones, ccp, &(up->zone)) ;
      else
	up->zone = 1 ;

      oligo = oligoR = 0 ;
      ccp = dictName (words, up->motif) ;
      i = strlen (ccp) ;
      dj = (i - wLn)/2 ;
      nok = 0 ; 
      if (extend)
	{
	  nok = extendLn ;
	  dj = i - (wLn -extendLn) ; /* terminal extendable word */
	}
      ccp += dj ;  
      dj = (i - wLn)/10 ; if (dj < 1) dj = 1 ;
      j = 0 ; ccp-- ;
      while (*++ccp)
	{
	  unsigned long long int cc ;
	  switch ((int)*ccp)
	    {
	    case 'a': case 'A': cc = 0 ; break ;
	    case 't': case 'T': cc = 3 ; break ;
	    case 'g': case 'G': cc = 1 ; break ;
	    case 'c': case 'C': cc = 2 ; break ;
	    default: nok = 0 ; oligo = oligoR = 0 ; continue ;
	    }
	  oligo = ((oligo << 2) | cc ) & mask ;
	  oligoR = ((oligoR >> 2) | ((3-cc) << shiftR)  ) & mask ;
	  if (++nok >= wLn && j < 1 && oligo)
	    {
	      if ((j % dj) == 0)
		{
		  int kk = 0 ;
		  unsigned long long int myOligo  ;

		  if (oligo & gtMask)
		    {
		      myOligo = oligo ;
		      kk = 0 ;
		    }
		  else
		    {
		      myOligo = oligoR ;
		      kk = 1 ;
		    }
		    
		  if (! extend)
		    assMultipleInsert (ass, assVoid (myOligo), assVoid (2 * ii + kk)) ; 
		  else
		    {
		      unsigned long long int xOligo = 0 , yOligo, eOligo = (oligo << (2 * extendLn)) ;
		      int i ;
		      int iMax = (1 << (2 * extendLn)) ;
		      for (i = 0 ; i < iMax ; i++)
			{
			  xOligo = i ;
			  yOligo = (eOligo | xOligo) & mask ;
			  assMultipleInsert (ass, assVoid (yOligo), assVoid (iMax * ii + i)) ;
			  }
		    }		    
		  j++ ;
		}
	    }
	}
    }
  ac_free (h) ;
  return ;
} /* tctValGetSnpList  */

/*************************************************************************************/
/* scan the fastc file
 * construct the words stranded or not
 * count them in the hash table
 */
static void tctValCountLaneOne (const TCT *tct, LANE *lane)
{
  Associator ass = tct->valAss ; 
  AC_HANDLE h = ac_new_handle () ;
  char *fNam1 = tct->fastcFile || ! tct->fastcDir ? 0 : 
    hprintf (h, "%s/%s.fastc%s"
	     , tct->fastcDir
	     , dictName (tct->laneDict, lane->lane)
	     , tct->runTest ? "" : ".gz"
	     ) ;
  const char *fNam =  fNam1 ? fNam1 : tct->fastcFile ;
  ACEIN ai = aceInCreate (fNam, FALSE, h) ;
  Array countsF = lane->valSnpCountsF ;
  Array countsR = lane->valSnpCountsR ;
  Array countsFF = lane->valSnpCountsFF ;
  Array countsRR = lane->valSnpCountsRR ;
  unsigned long long int oligo = 0, oligoR = 0, mask = 1, gtMask = 1 ;
  int wLn = tct->valWordLength ;
  int shiftR = 2 * (wLn - 1) ;
  int nok = 0 ;
  int line = 0 ;
  int mult = 1 ;
  int tails = tct->snpTails ;
  int decal = tct->snpExtend ? 1 : 2 ;
  int stranded = tct->stranded ; /* 0: non stranded, 1 :read1 is forward, 2: read1 is reverse */

  mask = (mask << 2*wLn) ; mask-- ;
  { int i = wLn/2 ; gtMask = (gtMask << 2*i) ;} 
  
  /* load the oligo, count only if we have wLn good letters */
  while (aceInCard (ai))
    {
      const char *cp0, *cp = aceInWord (ai) ;
      int strand = stranded ;
      line++ ;
      if (! cp)
	cp = "#" ;
      cp0 = cp ;
      switch ((int)*cp)
	{
	case 'a':
	case 't':
	case 'g':
	case 'c':
	case 'A':
	case 'T':
	case 'G':
	case 'C':
	  break ;
	case '>':
	  {
	    int k = 0 ;
	    char cc ;

	    mult = 1 ;
	    cp = strchr (cp, '#') ;
	    if (cp && sscanf (cp+1, "%d%c", &k, &cc) == 1 && k > 0)
	      mult = k ;
	  }
	  /* fall thru */
	default:
	  oligo = 0 ;
	  nok = 0 ; 
	  continue ; 
	}
      if (tct->inFormat == FASTQ && (line % 4) != 2)
	continue ;
      cp-- ;
      while (*++cp)
	{
	  unsigned long long int cc ;
	  switch ((int)*cp)
	    {
	    case 'a': case 'A': cc = 0 ; break ;
	    case 't': case 'T': cc = 3 ; break ;
	    case 'g': case 'G': cc = 1 ; break ;
	    case 'c': case 'C': cc = 2 ; break ;
	    case '>': if (cp[1] == '<') strand = - strand ;
	      /* fall through */
	    default: nok = 0 ; oligo = oligoR = 0 ; continue ;
	    }
	  oligo = ((oligo << 2) | cc ) & mask ;
	  oligoR = ((oligoR >> 2) | ((3-cc) << shiftR)  ) & mask ;
	  if (++nok >= wLn)
	    {
	      Array bucket = 0 ;
	      int iBucket = 0 ;
	      unsigned long long int myOligo ;
	      int kk = 0 ;
	      const void *op ;
	      const void *vp ;

	      if (oligo & gtMask)
		{
		  myOligo = oligo ;
		  kk = 0 ;
		}
	      else
		{
		  myOligo = oligoR ;
		  kk = 1 ;
		}
	      
	      op = assVoid (myOligo) ;

	      if (assFind (ass, op, 0))   /* (myass) */
		while (assFindNext(ass, op, &vp, &bucket, &iBucket))
		  {
		    int w = assInt(vp) ;
		    int kkk = ((w & 0x1) == kk) ? 0 : 4 ;

		    array (kkk ? countsR : countsF, w/decal, int) += mult ;
		    if (strand) 
		      {  
			int kkk1 = (strand == -1 ? 4 - kkk : kkk) ;
			array (kkk1 ? countsRR : countsFF, w/decal, int) += mult ;
		      }
		    if (tails)
		      {
			int i, k ;
			char bb[] = "atgc" ;
			Array a ;
			const char *cr ;

			w /= 2 ;
			for (k = 0 ; k < 4 ; k++)
			  {  /* suffix */
			    a = array (lane->tails[kkk+k], w, Array) ;
			    if (! a)
			      a = array (lane->tails[kkk+k], w, Array) = arrayHandleCreate (2*tails + wLn, int, lane->h) ;
			    for (i = tails , cr = cp - wLn + 1 ; *cr && i < 2*tails + wLn ; i++, cr++)
			      if (*cr == bb[k])
				array (a, i, int) += mult ;
			    for (i = tails-1 , cr = cp - wLn ; cr >= cp0 && i>=0 ; i--, cr--)
			      if (*cr == bb[k])
				array (a, i, int) += mult ;
			  }
		      }
		  }
	    }
	}
    }
  ac_free (h) ;	      
  return ;
} /* tctValCountLaneOne */

/*************************************************************************************/

static void tctValCountLane (const void *vp)
{
  const TCT *tct = (const TCT *)vp ;
  int ii ;
  BOOL debug = FALSE ;
  int mask2 = 1 << (2 * (tct->snpExtend ? tct->snpExtendLength : 0)) ;
  int nn = mask2 * dictMax (tct->valSnpNames) + 1 ;
  int tails = tct->snpTails ;

  while (channelGet (tct->analyzeChan, &ii, int))
    {
      LANE *lane = arrayp (tct->lanes, ii, LANE) ;
      if (debug)
	fprintf (stderr, "=== tctValCountLanes %d %s\n",ii, timeShowNow()) ;
      lane->valSnpCountsF = arrayHandleCreate (nn, int, lane->h) ;
      lane->valSnpCountsR = arrayHandleCreate (nn, int, lane->h) ;
      array (lane->valSnpCountsF , nn - 1, int) = 0 ;
      array (lane->valSnpCountsR , nn - 1, int) = 0 ;
      if (tct->stranded)
	{
	  lane->valSnpCountsFF = arrayHandleCreate (nn, int, lane->h) ;
	  lane->valSnpCountsRR = arrayHandleCreate (nn, int, lane->h) ;
	  array (lane->valSnpCountsFF , nn - 1, int) = 0 ;
	  array (lane->valSnpCountsRR , nn - 1, int) = 0 ;
	}
      if (tails)
	{
	  int k ;
	  for (k = 0 ; k < 8 ; k++)
	    {
	      lane->tails[k] = arrayHandleCreate (nn, Array, lane->h) ;
	      array (lane->tails[k], tails - 1 , Array) = 0 ;
	    }
	}
      tctValCountLaneOne (tct, lane) ;
    }
  ii = 100 ; /* success */
  channelPut (tct->doneChan, &ii, int) ;
  return ;
} /* tctValCountLane */

/*************************************************************************************/
/* parse a SNP or indel table
 * Construct a list of words of length snpValN using the genome 
 * Hash those words and number thep 0 to K
 * In parallel, read the fastc files and count per lane: array int length K
 * Cumul per run
 * Report
 */

static void tctValCount (TCT *tct)
{
  int n = 0, nn = 4 ; /* number of analysers */
  int ii ;
  int iMax = arrayMax (tct->lanes) ;

  tct->analyzeChan = channelCreate (30, int, tct->h) ;
  for (ii = 0 ; ii < nn ; ii++)
    wego_go (tctValCountLane, tct, TCT) ;  /* will channelPut (tct->doneChan) */ 

  for (ii = 1 ; ii < iMax ; ii++)
    channelPut (tct->analyzeChan, &ii, int) ;
  channelClose (tct->analyzeChan) ;

  while (nn--)
    if (! channelGet (tct->doneChan, &n, int) || n != 100)
    {
      fprintf (stderr, "tctValCount received a bad value %d\n", n) ;
      exit (1) ;
    }
  return ;
} /* tctValCount */

/*************************************************************************************/

static void tctValInit (TCT *tct, BOOL extend)
{
  int nn = tct->t2 - tct->t1 + 1 ;

  if (tct->valSnpListFileName)
    {
      if (! filCheckName (tct->valSnpListFileName, 0, "r"))
	{
	  fprintf (stderr, "Could not open -wordFile %s\n", tct->valSnpListFileName) ;
	  exit (1) ;
	}
  }

  tctGetLaneList (tct) ;

  wego_go (tctGetTargetFasta, tct, TCT) ;  /* will channelPut (tct->doneChan) */
  if  (! channelGet (tct->doneChan, &nn, int) || nn < 0)
    {
      fprintf (stderr, "tctValInit received a bad value %d\n", nn) ;
      exit (1) ;
    }
  else
    tct->t1 = nn ;

  tctValGetSnpList (tct, extend) ; /* requires targetFasta to construct the words */  
  return ;
} /* tctValInit */

/*************************************************************************************/
/*
gcagtaaggatggctagtgtaactagcaagaataccacgaaagcaagaaaaagaagtacgctattaactattaacgtacctgtctcttccGAAACGAatgagtacataagttcgtttagagaacagatctacaaga
*/
static void tctValReport (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (tct->outFileName, (tct->snpTails ? ".tails.txt" : (tct->snpExtend ? ".extend.txt" : ".val.tsf")) , tct->gzo, h) ;
  Array lanes = tct->lanes ;
  Array words = tct->valSnps ;
  int ii, w, wMax = arrayMax (words) ;
  Array countsFF = 0, countsF = arrayHandleCreate (wMax, int, h) ;
  Array countsRR = 0, countsR = arrayHandleCreate (wMax, int, h) ;
  int tails = tct->snpTails ;
  Array runTails[8] ;
  const int stranded = tct->stranded ;
  char suffix[1 << (2*tct->snpExtendLength)][tct->snpExtendLength + 1] ;
  memset (suffix, 0, sizeof(suffix)) ;

  for (ii = 0 ; ii < (1 << (2*tct->snpExtendLength)) ; ii++)
    {
      int i ;
      char *cc = "agct" ;
      
      for (i = 0 ; i < tct->snpExtendLength ; i++)
	suffix[ii][i] = cc[(ii >> 2*(tct->snpExtendLength-i-1) & 0x3)] ;
      suffix[ii][i] = 0 ;
    }

  if (tails)
    {
      int k ;
      for (k = 0 ; k < 8 ; k++)
	{
	  runTails[k] = arrayHandleCreate (tails, Array, h) ;
	  array (runTails[k], wMax - 1, Array) = 0 ;
	}
    }

  /* cumulate the lane count into run counts */
  if (wMax) 
    {
      array (countsF, wMax - 1, int) = 0 ;
      array (countsR, wMax - 1, int) = 0 ;
      if (stranded)
	{
	  countsFF = arrayHandleCreate (wMax, int, h) ;
	  countsRR = arrayHandleCreate (wMax, int, h) ;
	  array (countsFF, wMax - 1, int) = 0 ;
	  array (countsRR, wMax - 1, int) = 0 ;
	}
    }
  for (ii = 1 ; ii < arrayMax (lanes) ; ii++)
    {
      register int *ipF, *jpF, *ipR, *jpR ;
      register int *ipFF, *jpFF, *ipRR, *jpRR ;
      LANE *vp = arrp (lanes, ii, LANE) ;
      Array laneCountsF =  vp->valSnpCountsF ;
      Array laneCountsR =  vp->valSnpCountsR ;
      Array laneCountsFF =  vp->valSnpCountsFF ;
      Array laneCountsRR =  vp->valSnpCountsRR ;
      int iMax = arrayMax (laneCountsF) ;
      if (iMax < arrayMax (laneCountsR))
	iMax = arrayMax (laneCountsR) ;
      w = array (laneCountsF, iMax - 1, int) ; /* make room */
      w = array (countsF, iMax - 1, int) ; /* make room */
      w = array (laneCountsR, iMax - 1, int) ; /* make room */
      w = array (countsR, iMax - 1, int) ; /* make room */
      ipF = arrp (countsF, 0, int) ; jpF = arrp (laneCountsF, 0, int) ; ipR = arrp (countsR, 0, int) ; jpR = arrp (laneCountsR, 0, int) ;
      if (stranded)
	{
	  w = array (laneCountsFF, iMax - 1, int) ; /* make room */
	  w = array (countsFF, iMax - 1, int) ; /* make room */
	  w = array (laneCountsRR, iMax - 1, int) ; /* make room */
	  w = array (countsRR, iMax - 1, int) ; /* make room */
	  ipFF = arrp (countsFF, 0, int) ; jpFF = arrp (laneCountsFF, 0, int) ; ipRR = arrp (countsRR, 0, int) ; jpRR = arrp (laneCountsRR, 0, int) ;
	}
      for (w = 0 ; w < iMax ; w++)
	{
	  ipF[w] += jpF[w] ;
	  ipR[w] += jpR[w] ;
	  if (stranded)
	    {
	      ipFF[w] += jpFF[w] ;
	      ipRR[w] += jpRR[w] ;
	    }
	  
	  if (tct->snpTails && (ipF[w] + ipR[w] > 0))
	    {
	      int k ;
	      for (k = 0 ; k < 8 ; k++)
		{
		  Array a = array (runTails[k], w, Array) ;
		  Array b = array (vp->tails[k], w, Array) ;
		  int i, iMax = b ? arrayMax (b) : 0 ;
		  
		  if (!a)
		    a = array (runTails[k], w, Array) = arrayHandleCreate (tails, int, h) ;
		  for (i = 0 ; i < iMax ; i++)
		    array (a, i, int) += arr (b, i, int) ;
		}
	    }
	}
    }
  /* export in txt format */
  if (1)
    {
      const char *runName = tct->run ? tct->run : "XXX" ;
      VW *wp ;
      int *ipF, *ipR ;
      int *ipFF = 0, *ipRR = 0 ;
      DICT *names = tct->valSnpNames ;
      DICT *zones = tct->valSnpZones ;
      DICT *wordDict = tct->valSnpWords ;
      int extendLn = tct->snpExtendLength ;
      int mask1 = extendLn ? (1 << 2 * extendLn) - 1 : 0 ;
      int mask2 = mask1 + 1 ;
      Array hits = arrayHandleCreate (mask2, HIT, h) ;

      ipF = arrp (countsF, 0, int) ; ipR = arrp (countsR, 0, int) ; wp = arrp (words, 0, VW) ;
      if (stranded)
	{ ipFF = arrp (countsFF, 0, int) ; ipRR = arrp (countsRR, 0, int) ; }
      for (w = 0 ; w < wMax ; wp++, ipF += mask2, ipR += mask2, ipFF += mask2, ipRR += mask2, w++)
	if (wp->motif)
	  {
	    if (extendLn)
	      {
		int i, kF = 0, kR = 0 ;
		for (i = 0 ; i < mask2 ; i++)
		  {
		    HIT *up = arrayp (hits, i, HIT) ;
		    up->a1 = -ipF[i] -ipR[i] ; up->a2 = i ;
		    up->x1 = ipF[i] ; up->x2 = ipR[i] ;
		    kF += ipF[i] ; kR += ipR[i] ;
		  }
		arraySort (hits, tctA1Order) ;
		aceOutf (ao, "%s\t%s\t%s\t%d\t%d\t%s", dictName (zones, wp->zone), dictName (names, wp->name), runName, kF, kR, dictName (wordDict, wp->motif)) ;      
		kF = arrp (hits, 0, HIT)->a1 ;
		for (i = 0 ; i < mask2 ; i++)
		  {
		    HIT *up = arrp (hits, i, HIT) ;
		    if(32*up->a1 < kF) 
		      aceOutf (ao, " %d:%d:%s", up->a2, -up->a1,suffix[up->a2]) ;
		  }
		aceOut (ao, "\n") ;
	      }
	    else if (tails)
	      {
		int i, k, kk = 0 ;
		int wLn = tct->valWordLength ;		
		char *ss[] = {"3prime Forward", "5prime Complement"} ;
		char atgc[] = "atgc" ;
		
		if (*ipF + *ipR <= 0)
		  {
		    aceOutf (ao, "%s\t%s\t%s\t%d\t%d\t%s\n", dictName (zones, wp->zone), dictName (names, wp->name), runName, *ipF, *ipR, dictName (wordDict, wp->motif)) ;      
		    continue ;
		  }
		for (kk = 0 ; kk < 8 ; kk += 4) 
		  {
		    char buf[wLn + 2*tails + 2] ;
		    int last = wLn  ;
		    
		    for (i = 0 ; i < 2*tails + wLn ; i++)
		      {
			int cc = 0, s = 0, bestK = 0, nK = 0 ;
			for (k = 0 ; k < 4 ; k++)
			  {
			    Array a = array (runTails[kk+k], w, Array) ;
			    int iMax = a ? arrayMax (a) : 0 ;
			    int m = (i < iMax ? arr (a, i, int) : 0) ;
			    if (m > nK)
			      { bestK = k ; nK = m ; }
			    s += m ;
			  }
			cc = '-' ;
			if (100 * nK > 50 * s)
			  { cc = atgc[bestK] ; last = i+ 1 ;}
			if (100 * nK > 95 * s)
			  cc = ace_upper(cc) ;
			buf[i] = cc ;
		      }
		    buf[last] = 0 ;
		    aceOutf (ao, "%s\t%s\t%s\t%d\t%d\t%s", dictName (zones, wp->zone), dictName (names, wp->name), runName, *ipF, *ipR, dictName (wordDict, wp->motif)) ;      
		    aceOutf (ao, "%s\tExtended_word %s\t%s\n", kk == 4 ? "\t\t\t" : "" ,ss[kk/4], buf) ;
			
		    for (k = 0 ; k < 4 ; k++)
		      {
			int nn = kk==0 ? *ipF : *ipR ;
			aceOutf (ao, "%s\t%s\t%s\t%s\t%d\t%s\t%s\t%c", dictName (zones, wp->zone), dictName (names, wp->name), nn, dictName (wordDict, wp->motif), runName, ss[kk/4], atgc[k]) ;      
			if (nn)
			  {
			    Array a = array (runTails[kk+k], w, Array) ;
			    int iMax = a ? arrayMax (a) : 0 ;
			    for (i = 0 ; i < iMax ; i++)
			      aceOutf (ao, "\t%d", arr (a, i, int)) ;
			    aceOut (ao, "\n") ;
			  }
		      }
		  }
	      }
	    else if (stranded)
	      aceOutf (ao, "%s\t%s\tiiii\t%d\t%d\t%d\t%d\n", dictName (names, wp->name), runName, *ipF, *ipR, *ipFF, *ipRR) ;
	    else
	      aceOutf (ao, "%s\t%s\tii\t%d\t%d\n", dictName (names, wp->name), runName, *ipF, *ipR) ;
	  }
    }
  ac_free (h) ;
  return ;
} /*tctValReport */
/*
echo "ref\tacacgtaaccctgcttggagaaaagctgtct" > toto2.motif
echo "var\tacacgtaaccTtgcttggagaaaagctgtct" >> toto2.motif
ttacacgtaaccTtgcttggag
*/

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

static void tctGetTestOffset (TCT *tct)
{
  if (tct->runTest) 
    {
       AC_HANDLE h = ac_new_handle () ;
       ACEIN ai = aceInCreate (hprintf (h, "%s.t1", tct->runTest), 0, h) ;

       while (aceInCard (ai))
	 {
	   char *cp = aceInWord (ai) ;
	   if (cp && !strcmp (cp,  "##OFFSET"))
	     {
	       aceInStep (ai, '\t') ;
	       aceInInt (ai, &(tct->offset)) ;
	       aceInStep (ai, '\t') ;
	       if (! tct->t2) aceInInt (ai, &(tct->t2)) ;
	       break ;
	     }
	   ac_free (ai) ;
	   ac_free (h) ;
	 }
     }
} /* tctGetTestT1 */

/*************************************************************************************/

static void tctInit (TCT *tct)
{
  int nn = tct->t2 - tct->t1 + 1 ;

  if (! nn) nn = 10000000 ;
  tct->dna = arrayHandleCreate (nn, char, tct->h) ;
  tct->dnaR = arrayHandleCreate (nn, char, tct->h) ;
  if (tct->wantDeUno)
    {
      tct->wig = arrayHandleCreate (tct->dna ? arrayMax (tct->dna) : 1000, int, tct->h) ;
      tct->deUno = bigArrayHandleCreate (1000, UNO, tct->h) ;
      tct->deUnoDict = dictHandleCreate (1000, tct->h) ;
    }

  if (tct->geneFusion)
    tct->overhangDict = dictHandleCreate (1000, tct->h) ;
  tct->targetDict = dictHandleCreate (1000, tct->h) ;

  if (tct->externalSnpListFileName)
    tctGetExternalSnpList (tct) ;

  if (tct->nAna > 1)
    dictSetLock (tct->targetDict, channelLockCreate (tct->h), channelLock, channelUnlock) ;

  if (! tct->geneDict)
    tct->geneDict = dictHandleCreate (128, tct->h) ;
  if (tct->nAna > 1)
    dictSetLock (tct->geneDict, channelLockCreate (tct->h), channelLock, channelUnlock) ;
  tct->genes = arrayHandleCreate (20000, GENE, tct->h) ;


  tct->targetStart = keySetHandleCreate (tct->h) ;
  tct->targetStop = keySetHandleCreate (tct->h) ;
  if (tct->geneFusion)
    tct->geneFusions = arrayHandleCreate (1024, GF, tct->h) ;

  if (tct->target)
    dictAdd (tct->targetDict, tct->target, &tct->targetKey) ;
  if (tct->geneMapFileName) 
    tctGetGeneMap (tct) ;
  if (tct->mrnaRemapFileName)
    tctGetMrnaRemap (tct) ;
  if (tct->splitMrnaFileName)
    tctGetSplitMrnas (tct) ;

  tctGetLaneList (tct) ;

  wego_go (tctGetTargetFasta, tct, TCT) ;  /* will channelPut (tct->doneChan) */
  if  (! channelGet (tct->doneChan, &nn, int) || nn < 0)
    {
      fprintf (stderr, "tctInit received a bad value %d\n", nn) ;
      exit (1) ;
    }
  else
    tct->t1 = nn ;

  /* the errors only exist if we know the genome */
    {
      tctGetHits (tct) ;
      tctGetBrks (tct) ;  
    }

  return ;
} /* tctInit */

/*************************************************************************************/

static void tctCleanUp (const TCT *tct)
{
  int ii, iMax = tct->lane ? arrayMax (tct->lanes) : 0 ;
  LANE *lane ;

  if (iMax)
    for (ii = 0, lane = arrp (tct->lanes, 0, LANE) ; ii < iMax ; ii++, lane++)
      ac_free (lane->h) ;
  return ;
} /* tctCleanUp */
  
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
	   "// Usage: variant_caller -run runName ..... \n"
	   "//      try: -h --help \n"
	   "// Example\n"
	   "//     variant_caller -run SRR123456 -t chr3 -f chr3.fasta\n"
	   "//     variant_caller -run SRR123456 -count -words f -wLn 21\n"
	   "// Several sets of parameters must or may be specified\n"
	   "// READS TO BE ANALYZED: \n"
	   "//   -run runName : [mandatory]\n"
	   "//   -laneList fileName :  [optional] list of lanes to be analyzed\n"
	   "//       [defaults to $fastcDir/$run/LaneList]\n"
	   "//   -maxLanes int: analyze at most n lanes per run\n"
	   "//    example:\n"
	   "//        -run SRR123456 -laneList Fastc/SRR123456/LaneList\n"
	   "//   -fastcFile file_name : a single fasta/fastc file for the short reads\n"
	   "//   -fastcDir dirName : the directory containing the .fastc.gz files \n"
	   "//      [optional, default Fastc/$run]\n"
	   "// ALIGNMENTS TO BE ANALYZED\n"
	   "//   -SAM : SAM format\n"
	   "//     the sam files are $hitDir/$lane.sam.gz\n"
	   "//   -SAM_file file_name : the single input alignment files are in SAM format\n"
	   "//   -hitFile file_name : a single hit file exported by clipalign\n"
	   "//    example:\n"
	   "//        -run SRR123456 -hitFile tmp/COUNT/SRR123456/f2.3.hits.gz\n"
	   "//   -hitDir dirName : the directory containing the .hits.gz files \n"
	   "//       [optional, default tmp/COUNT/$run]\n"
	   "//       The program expects that for each lane listed in the laneList file\n"
	   "//       there exists a file called\n"
	   "//          $hitDir/$lane.hits.gz      containingg the magic hits\n"
	   "//          $fastcDir/$lane.fastc.gz   containing the read sequences]n"
	   "// TARGET ZONE:\n"
	   "//   -target_fasta fileName : the name of the reference fasta file, possibly gzipped \n"
	   "//   -target_class : target class (A_mito ... Z_genome) [default Z_genome]\n"
	   "//   -t targetName :  [optional] often a chromosome name\n"
	   "//   -t1 int -t2 int :  [optional] analyze only section [t1,t2] of the target\n"
	   "//    example:\n"
	   "//        -target_class  Z_genome -t chr7 -t1 2000000 -t2 3000000\n"
	   "//       Typically t2 - t1 = 10 Mbases, or if t1, t2 are not specified\n"
	   "//       the whole chromosome is analyzed .\n"
	   "// SNP FILTERS, DISCOVERY PHASE\n"
	   "//   -method arbitrary_name: the method is echoed in the output file\n" 
	   "//   -externalSnpList filename: [optional] additional list of positions to be analysed\n"
	   "//      several formats are recognized\n"
	   "//          chr7:1234     chr7:1234:....   chr7\\t1234\n"
	   "//          one position per line\n"
	   "//   -ddx integer : default 1, must be in [0, 3] interval, use only when debugging\n"
	   "//      It does not matter if the positions are off by +- ddx bases\n"
	   "//   -minAli : [default 25] discard hits shorter than the specified limit in bp\n"
	   "//   -minAliPercent :  discard reads not aligning on a n%% of their clipped lenght and ali < 140\n"
 	   "//   -minSnpCover integer : min coverage [default 10] \n"
	   "//   -minSnpCount integer: [default 4]\n"
	   "//   -minSnpFrequency float: [default 18] minimal MAF precentage\n"
	   "//   -minOverhang integer: [default 12] collect positions of incomplete alignments\n"
	   "//   -intron : [default off] diifferentiate introns from deletions\n"
	   "//   -intron_only : just detect and report introns\n"
	   "//   -min_intron <int> : [default 30] min reported intron length\n"
	   "//   -max_intron <int> : [default 0]  max reported intron length\n"
	   "//   -min_bridge <int> : [default 0]  only report lon ger introns or deletions\n"
	   "//   -dx int ; [default 5] minimal acceptable value 5\n"
	   "//      A new phased segment starts when the next SNV is furher than dx\n"
	   "// COUNT\n"
	   "//   -count : count in the raw fastc/fastq files the occurence of some words\n"
	   "//       See the first few lines of this help to see how to specify\n"	
	   "//       the fasta/fasc/fastq files that should be searched\n"
	   "//         The words and their complements are counted separatelly\n"
	   "//       This may show that a SNP is seen as mutant on one starnd, wild-type on the other\n"
	   "//       possibly indicating a systematic sequencing error, But if one specifies\n"
	   "//   -strand : optional, flip read-2 but not read-1\n"
	   "//   -antistrand : optional, flip read-1 but not read-2\n"
	   "//       Then one counts the stranded support of each word. For example, in RNA seq,\n"
	   "//       we can know the measure of an intron, or a deletion.\n"
	   "//   -wLn <int> : [default 17] word length, must be odd\n"
	   "//   -words <a,b> : coma delimited list of words, example: atgcg,tctgg\n"
	   "//   -wordFile fileName : sequence to be searched in the fastc files\n"
	   "//      One or two columns, tab delimited:\n"
	   "//              atgctgac[mandatory]   identifier[optional]\n"
	   "//      Example:\n"
	   "//        -count -wLn 23 -wordFile fNam -run r -o toto -gzo\n"
	   "//      For each identifier in file fNam, export in toto.val.txt.gz\n"
	   "//      the number of occurence of the words listed in the second column.\n"
	   "//      Only the central subword of length wLn is considered.\n"
	   "//   -tails <int n>\n"
	   "//      Search the fast[acq] files for exactly the given words and report\n"
	   "//      the concensus and base profiles extended by n bases in both directions\n"
	   "// GENE FUSION\n"
	   "//   -uno : detect SNPs and short indels detected by the aligner\n"
	   "//   -bridge : search indels, introns, fusion rerrangments\n"
	   "//   -target_class : target class (KT_RefSeq ET_av...) [default ET_av]\n"
	   "//   -min_GF integer : [default 5]  filter geneFusionFile on min support \n" 
	   "//   -minOverhang integer : [default 15] minimal number of bases\n"
	   "//   -geneFusion : scan th hits, report a geneFusionFile per read\n"
	   "//   -geneFusionFile fileName: file of genefusions to be analysed\n"
	   "//      mrna1 a1 a2 mrna2 b1 b2 n (n supports for a jump from mrna1[position a2] to m2[b1]\n"
	   "//      Scan the hit file(s) report for each donor/acceptor read count that support\n"
	   "//         the proposed donor and goes to the acceptor\n"
	   "//         OR align locally OR jump locally OR jump elsewhere\n"
	   "//    -geneMap fileName : tmp/METADATA/$target.mrna_map_ln_gc_gene_geneid.txt\n"
	   "//    -mrnaRemap fileName : tmp/METADATA/mrnaRemap.gz\n"
	   "//    -splitMrnas fileName : tmp/METADATA/splitMrnas.gz\n"
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
  TCT tct ;
  AC_HANDLE h = 0 ;
  char commandBuf [4000] ;

  freeinit () ; 
  messErrorInit (argv[0]) ;

  if (sizeof(void*) != 8)
    messcrash ("This program only works on 64bits machines, sorry") ;

  h = ac_new_handle () ;
  memset (&tct, 0, sizeof (TCT)) ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
  tct.h = h ;

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
      fprintf (stderr, "variant_caller: %s\n", VERSION) ;
      exit (0) ;
    }
 
  {
    int ix ;
    char *cp ;
    for (ix = 0, cp = commandBuf ;  ix < argc && cp + strlen (argv[ix]) + 1 < commandBuf + 3900 ; cp += strlen (cp), ix++)
      sprintf(cp, "%s ", argv[ix]) ;
  }

  /* consume optional args */
  tct.gzi = getCmdLineBool (&argc, argv, "-gzi") ;
  tct.gzo = getCmdLineBool (&argc, argv, "-gzo") ;
  tct.debug = getCmdLineBool (&argc, argv, "-debug") ;


  getCmdLineOption (&argc, argv, "-o", &(tct.outFileName)) ;
 
  tct.qc = getCmdLineBool (&argc, argv, "-qc") ;
  tct.wantDeUno = getCmdLineBool (&argc, argv, "-uno") ;
  tct.wantBrks = getCmdLineBool (&argc, argv, "-bridge") ;
  tct.SAM = getCmdLineBool (&argc, argv, "-SAM") ;
  if (getCmdLineOption (&argc, argv, "-SAM_file", &(tct.SAM_file)))
  tct.tabular = getCmdLineBool (&argc, argv, "-tabular") ;

  if (getCmdLineOption (&argc, argv, "-samFile", &(tct.SAM_file)))
    tct.SAM = TRUE ;
  getCmdLineOption (&argc, argv, "-hitFile", &(tct.hitFile)) ;
  getCmdLineOption (&argc, argv, "-fastcFile", &(tct.fastcFile)) ;
  if (getCmdLineOption (&argc, argv, "-fastaFile", &(tct.fastcFile)))
    tct.inFormat = RAW ;
  if (getCmdLineOption (&argc, argv, "-fastaFile", &(tct.fastcFile)))
    tct.inFormat = FASTA ;
  if (getCmdLineOption (&argc, argv, "-fastqFile", &(tct.fastcFile)))
    tct.inFormat = FASTQ ;
  getCmdLineOption (&argc, argv, "-geneMap", &(tct.geneMapFileName)) ;
  getCmdLineOption (&argc, argv, "-mrnaRemap", &(tct.mrnaRemapFileName)) ;
  getCmdLineOption (&argc, argv, "-splitMrnas", &(tct.splitMrnaFileName)) ;
  tct.geneFusion = getCmdLineBool (&argc, argv, "-geneFusion") ;
  getCmdLineOption (&argc, argv, "-geneFusionFile", &(tct.geneFusionFileName)) ;
  getCmdLineOption (&argc, argv, "-geneFusionPositions", &(tct.geneFusionPositionsFileName)) ;

  /* TARGET ZONE */
  tct.target_class = 0 ; /* "Z_genome" ; */
  getCmdLineOption (&argc, argv, "-target_class", &(tct.target_class)) ;
  getCmdLineOption (&argc, argv, "-t", &(tct.target)) ;
  getCmdLineOption (&argc, argv, "-make_test", &(tct.makeTest)) ;

  getCmdLineOption (&argc, argv, "-method", &(tct.method)) ;
  
  if (getCmdLineOption (&argc, argv, "-run_test", &(tct.runTest)))
    {
      tct.targetFileName = hprintf (h, "%s.target.fasta", tct.runTest) ;
      tct.laneListFileName = hprintf (h, "%s.LaneList", tct.runTest) ;
      tct.run = tct.runTest ;
      tct.target = tct.runTest ;
      tct.hitDir = tct.fastcDir = "." ;
      tct.target_class = 0 ;
      tct.wantBrks = TRUE ;
    }

  getCmdLineOption (&argc, argv, "-lane", &(tct.lane)) ;
  getCmdLineOption (&argc, argv, "-run", &(tct.run)) ;
  if (! tct.run)
   {
     fprintf (stderr, "FATAL ERROR: Missing argument -run, please try\n\tTCT -help\n") ;
     exit (1) ;
   }
  

  getCmdLineInt (&argc, argv, "-min_GF", &(tct.minGF)) ;
  tct.intron_only = getCmdLineBool (&argc, argv, "-intron_only") ;

  tct.minIntron = 30 ; /* default */
  tct.maxIntron = 0 ; /* default : do not search introns */
  if (getCmdLineBool (&argc, argv, "-intron"))
    tct.maxIntron = 100000 ; /* default if we search introns */
  getCmdLineInt (&argc, argv, "-min_bridge", &(tct.minBridge)) ;
  getCmdLineInt (&argc, argv, "-min_intron", &(tct.minIntron)) ;
  getCmdLineInt (&argc, argv, "-max_intron", &(tct.maxIntron)) ;

  tct.minAli = 25 ;
  getCmdLineInt (&argc, argv, "-minAli", &tct.minAli) ;
  getCmdLineInt (&argc, argv, "-minAliPerCent", &tct.minAliPerCent) ;


  getCmdLineOption (&argc, argv, "-target_fasta", &(tct.targetFileName)) ;
  getCmdLineInt (&argc, argv, "-t1", &(tct.t1)) ;
  getCmdLineInt (&argc, argv, "-t2", &(tct.t2)) ;

  tct.snpCount = getCmdLineBool (&argc, argv, "-count") ;
  if (getCmdLineInt (&argc, argv, "-extend", &(tct.snpExtendLength)))
    {
      tct.snpExtend = TRUE ;
      tct.snpCount = TRUE ;
      if (tct.snpExtendLength < 1 || tct.snpExtendLength> 6)
	{
	  fprintf (stderr, "FATAL ERROR: argument -extend must be in [1,6], not %d\n" , tct.snpExtendLength) ;
	  exit (1) ;
	}      
    }
  if (getCmdLineInt (&argc, argv, "-tails", &(tct.snpTails)))
    {
      tct.snpCount = TRUE ;
      if (tct.snpTails < 1 || tct.snpTails > 256)
	{
	  fprintf (stderr, "FATAL ERROR: argument -tails must be in [1,256], not %d\n" , tct.snpTails) ;
	  exit (1) ;
	}      
    }
  
  getCmdLineOption (&argc, argv, "-wordFile", &(tct.valSnpListFileName)) ;
  getCmdLineOption (&argc, argv, "-words", &(tct.valSnpList)) ;
  tct.valWordLength = 17 ; /* default */
  if (getCmdLineBool (&argc, argv, "-strand")) tct.stranded = 1 ;
  if (getCmdLineBool (&argc, argv, "-antistrand")) tct.stranded = -1 ;
  if (getCmdLineInt (&argc, argv, "-wLn", &(tct.valWordLength)))
    tct.valWordLength |= 0x1 ; /* make it odd */
  if (tct.valWordLength > 31)
    messcrash ("Validation word length %d cannot exceed 31", tct.valWordLength) ;
  if (! tct.target && ! tct.geneFusion && ! tct.snpCount && ! tct.qc)
    {
      fprintf (stderr, "FATAL ERROR: Missing argument -t target_name, please try\n\tTCT -help\n") ;
      exit (1) ;
    }

  if (! tct.targetFileName && tct.makeTest)
    { 
      fprintf (stderr, "FATAL ERROR: Missing argument -target_fasta while using make_test, please try\n\tvariant_caller -help\n") ;
      exit (1) ;
    }

   
  if (tct.t2 < tct.t1)
    {
      fprintf (stderr, "FATAL ERROR: zone limits expeted t1 < t2, got %d > %d, please try\n\tvariant_caller -help\n", tct.t1, tct.t2) ;
      exit (1) ;
    }
  if (tct.t1 < 0)
    {
      fprintf (stderr, "FATAL ERROR: zone limits expeted t1 > 0 got %d  please try\n\tvariant_caller -help\n", tct.t1) ;
      exit (1) ;
    }

  /* READ SET */
  if (! tct.fastcFile && ! tct.SAM_file && ! tct.hitFile && ! tct.runTest && !tct.geneFusionFileName)
    { /* usual magic configuration defaults */
      tct.laneListFileName = hprintf (h, "Fastc/%s/LaneList", tct.run) ;
      tct.fastcDir = "Fastc" ;
      tct.hitDir = tct.SAM ? "tmp/MAGICBLAST" : "tmp/COUNT" ;
    }
  
  getCmdLineOption (&argc, argv, "-laneList", &(tct.laneListFileName)) ;
  getCmdLineOption (&argc, argv, "-hitDir", &(tct.hitDir)) ;
  getCmdLineOption (&argc, argv, "-fastcDir", &(tct.fastcDir)) ;

  getCmdLineOption (&argc, argv, "-externalSnpList", &(tct.externalSnpListFileName)) ;

  if (! tct.laneListFileName && ! tct.hitFile && ! tct.SAM_file && ! tct.snpCount && !tct.geneFusionFileName)
    fprintf (stderr, "FATAL ERROR: Missing argument -LaneList, please try\n\tvariant_caller -help\n") ;
	
  
  /* SNP FILTERS */
  tct.minSnpCover = 10 ;
  tct.minSnpCount = 4 ;
  tct.minSnpFrequency = 18 ;
  tct.minOverhang = 12 ;
  tct.Delta = 5 ;
  getCmdLineInt (&argc, argv, "-maxLanes", &(tct.maxLanes)) ;
  getCmdLineInt (&argc, argv, "-dx", &(tct.Delta)) ;
  if (tct.Delta < 5) 
    {
      fprintf (stderr, "-dx %d,  should be at least 5 bp\n", tct.Delta) ;
      exit (1) ;
    }
  
  tct.microDelta = 1 ;
  getCmdLineInt (&argc, argv, "-ddx", &(tct.microDelta)) ;
  if (tct.microDelta > 4) 
    {
      fprintf (stderr, "-ddx %d,  should be in [0,3] bp  interval\n", tct.microDelta) ;
      usage (commandBuf, argc, argv) ;
    }
 
  getCmdLineInt (&argc, argv, "-minOverhang", &(tct.minOverhang)) ;
  getCmdLineInt (&argc, argv, "-minSnpCover", &(tct.minSnpCover)) ;
  getCmdLineInt (&argc, argv, "-minSnpCount", &(tct.minSnpCount)) ;
  getCmdLineFloat (&argc, argv, "-minSnpFrequency", &(tct.minSnpFrequency)) ;
  getCmdLineInt (&argc, argv, "-subSampling", &(tct.subSampling)) ;
  tct.nAna = 4 ;
  getCmdLineInt (&argc, argv, "-nAna", &(tct.nAna)) ;
  
  
  tct.max_threads = 4 ;
  getCmdLineInt (&argc, argv, "-max_threads", &tct.max_threads) ;
  
  if (tct.minSnpCount > tct.minSnpCover)
    tct.minSnpCount =  tct.minSnpCover ;
  if (argc != 1)
    {
      fprintf (stderr, "unknown argument, sorry\n") ;
      usage (commandBuf, argc, argv) ;
    }
  /* check file existence */
  if (tct.SAM_file)
    {
      ACEIN ai = aceInCreate (tct.SAM_file, 0, h) ; 
      if (! ai)
	messcrash ("Cannot open SAM_file %s\n", tct.SAM_file) ;
      ac_free (ai) ;
    }
  if (tct.hitFile)
    {
      ACEIN ai = aceInCreate (tct.hitFile, 0, h) ; 
      if (! ai)
	messcrash ("Cannot open hitFile %s\n", tct.hitFile) ;
      ac_free (ai) ;
    }
  if (tct.fastcFile)
    {
      ACEIN ai = aceInCreate (tct.fastcFile, 0, h) ; 
      if (! ai)
	messcrash ("Cannot open hitFile %s\n", tct.hitFile) ;
      ac_free (ai) ;
    }
  if (tct.laneListFileName)
    {
      ACEIN ai = aceInCreate (tct.laneListFileName, 0, h) ; 
      if (! ai)
	messcrash ("Cannot open laneList %s\n", tct.laneListFileName) ;
      ac_free (ai) ;
    }
 if (tct.externalSnpListFileName)
    {
      ACEIN ai = aceInCreate (tct.externalSnpListFileName, 0, h) ; 
      if (! ai)
	messcrash ("Cannot open externalSnpList %s\n", tct.externalSnpListFileName) ;
      ac_free (ai) ;
    }
 if (0 && tct.geneFusionFileName && ! tct.geneFusionPositionsFileName)
   {
     messcrash (" -geneFusion fileName requires -geneFusionPositions fileName") ;
   }
 if (tct.geneFusionFileName)
    {
      ACEIN ai = aceInCreate (tct.geneFusionFileName, 0, h) ; 
      if (! ai)
	messcrash ("Cannot open -geneFusion file: %s\n", tct.geneFusionFileName) ;
      ac_free (ai) ;
      ai = aceInCreate (tct.geneFusionPositionsFileName, 0, h) ; 
      if (! ai)
	messcrash ("Cannot open -geneFusionPositions file: %s\n", tct.geneFusionPositionsFileName) ;
      ac_free (ai) ;
    }
  /* parallelization */
   if (tct.max_threads < 4)
     tct.max_threads = 4 ;
   if (tct.max_threads < tct.nAna)
     tct.max_threads = tct.nAna ;
 
  /* check the absolute args */

  wego_max_threads (tct.max_threads) ;
  tct.doneChan = channelCreate (12, int, tct.h) ;
  tctGetTestOffset (&tct) ;
  if (!tct.geneFusionFileName)
    tctGetLaneList (&tct) ;

  if (tct.wantDeUno)
    {
      tctInit (&tct) ;
      tctDeUnoExport (&tct) ;
    }
  else if (tct.geneFusionFileName)
    {
      tctGetGeneFusionList (&tct) ;
      if (arrayMax (tct.geneFusions))
	tctExportGeneFusion2tsf (&tct) ;
    }
  else if (tct.geneFusion)
    {
      tctInit (&tct) ;
      if (1 || arrayMax (tct.geneFusions))
	tctExportGeneFusions (&tct) ;
    }
  else if (tct.snpExtend || tct.snpTails)
    {
      tctValInit (&tct, tct.snpExtend) ;
      tctValCount (&tct) ;  /* hybridize in paralle the fastc files to the words */
      tctValReport (&tct) ; /* report */
    }
  else if (tct.snpCount)
    {
      tctValInit (&tct, FALSE) ;
      tctValCount (&tct) ;  /* hybridize in paralle the fastc files to the words */
      tctValReport (&tct) ; /* report */
    }
  else
    {
      tctInit (&tct) ;
      if (tct.wantDeUno)
	tctDeUnoExport (&tct) ;

      tctGetSegs (&tct) ;  
      tctAnalyzeSegs (&tct) ;
      tctMergeSubSegs (&tct) ;
      tctRationalize (&tct) ; 
      if (tct.outFileName) tctAliLn (&tct) ;
      tctExport (&tct) ;
    }
  tctCleanUp (&tct) ;
  wego_flush () ;

  ac_free (tct.ao) ;
  if (1)
    { 
      int mx ;
      messAllocMaxStatus (&mx) ; 
      fprintf (stderr, "// minSnpCount %d minSnpCover %d minSnpFrequency %.1f%% SNP wLn %d\n"
	       , tct.minSnpCount
	       , tct.minSnpCover
	       , tct.minSnpFrequency
	       , tct.valWordLength	     
	       ) ;
      if (tct.brks)
	fprintf (stderr, "// SNP detected %d BRKS analysed %d SEG reported %d\n"
		 , tct.brks ? arrayMax (tct.brks) : 0
		 , tct.segs ? arrayMax (tct.segs) : 0
		 , tct.snpExported
		 ) ;
      fprintf (stderr, "// %s done, %d errTrackings max memory %d Mb\n", timeShowNow(),  tct.nErrTrackings, mx) ;
     }
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderrn and to ensure all pipes are closed*/
  ac_free (tct.h) ;

  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
