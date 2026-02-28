/*  File: clipalign.c
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
 */

#define VERSION "2017_12"
#define MAXJUMP 8

#define LATESORT 0
#include "ac.h"
#include "bitset.h"
#include "bigarray.h"


#ifdef RESTRICT 
#define Restrict restrict
#else
#define Restrict 
#endif

/* max number of bases in prefix/suffix */
#define OVLN 30
#define GENE_STRAND_MALUS 3

typedef struct clipAlignStruct {
  AC_HANDLE h ;
  const char *probeFileName ;
  const char *probeFileName2 ;  /* pair case, one can give the sequences in 2 independant files */
  const char *targetFileName ;
  const char *multiTargetFileNames ;
  const char *selectedProbeFileName ;
  const char *rejectedProbeFileName ;
  const char *previousScoreFileName ;
  const char *targetMaskFileName ;
  const char *outFileName ;
  const char *target2geneFileName ;
  const char* wordFrequencyConstructTable ;
  const char* wordFrequencyTable ;
  const char *target_class ;
  const char *runName ;
  FILE *probeFile, *targetFile, *selectedProbeFile , *rejectedProbeFile ;
  Stack probeStack, qualityStack ;
  unsigned long int mask, mask1, mask2 ;
  int dx1, dx2, gap ;
  int target, targetGene, previousTarget, previousTargetGene, nTargets ;
  BitSet oligoGene ;
  BOOL oligoGeneTouched ;
  int wffSize ;
  unsigned char *wff ; /* word frequency table of words of length wffSize */
  int UMI ; /* length of cell bar code, clip behind, report in col 20, do not map read 1 */
  int jump5_1, jump5_2 ; /* jump that many base at the 5' or 3' end of the probes, for read 1 or read 2 */
  int forceRightClip_1, forceRightClip_2 ; /* force rightClip read 1 and/or read 2 */
  int clipAt ; /* clip all probes at that maximal length */
  int targetStart, targetStop; /* this range on the target, must be included, used to validate snps */
  int exactTargetStart, exactTargetStop; /* no error in this range on the target, used to search exon juntions */
  int overhangLength ; /* used in -splice mode to see that the overhanging prefix exists in the genome */
  int intronMaxLength ; /*  used in -splice mode when veryfying the overhang matches the genome */
  int intronMinLength ; /*  used in -splice mode to classify as intron or short deletion */
  int maxDiscardableAli ;
  BOOL showTargetPrefix ;
  BOOL probeSetName ; /* if provided, export the column */
  BOOL showProbeSequence, showOverhang ;
  BOOL clipPolyA, clipPolyA2 ; /* clip terminal polyA */
  BOOL clipPolyT, clipPolyT2 ; /* clip initial polyT */
  BOOL hasGenes ;
  BOOL MRNAH ; /* hierarchic mapping, drop other tarnscripts in alphabetic order */
  BOOL solid ;  /* we are aligning in color space, apply the error correcting procedure */
  BOOL roche ; /* favor deletions */
  BOOL pacbio ; /* favor deletions */
  BOOL nanopore ; /* favor deletions */
  BOOL ultima; /* favor deletions */
  BOOL fastQ ; /* probe file is in fastQ format */
  BOOL A2G, G2A, T2C, C2T, A2C, C2A, T2G, G2T, G2C, C2G, A2T, T2A, X2X ;  /* search for intense RNA-editing by aligning in 3 letter code */
  int probeQquality ; /* ignore, if fastq format, right clip until a letter of given quality is reached */
  int fastQbase ; /* base of the quality factors, may be 33 (NCBI) or 64 (ILM/solexa), if set, report the quality of the SNP */
  int nbpQclipped ;   /* how many bp where clipped on the right */
  int nseqQclipped ;  /* in how many sequence tags */
  int nTagEmptyInsert ;
  int nTagEmptyOneInsert[64] ;
  int nSeqEmptyInsert ;
  long int nBpEmptyInsert ;
  int stranded ; /* 1 for reads, 2 for pairs */
  BOOL selectPreviousScore ;
  BOOL hasPairs ;
  int pairTampon ;
  int strandBonus ; /* malus if the the tag is antisense relative to the target */
  int nReadPlus, nReadMinus ; /* number of reads observed on each strand of the transcriptome */
  int nIntronPlus, nIntronMinus ; /* number of stranded introns */
  BOOL strand, antistrand, antiprobe, aceOutGenomic, aceOutMrna, aceOutNM, silent ;
  BOOL vectorize ; /* export the alignments -a la Jim Kent- with coma separated muti coordinates */ 
  BOOL bestHit, splice, getTm, golay, multi, preScanGenome, preScanProbes, strandedTarget ;
  BOOL sam ; /* sam output requested */
  BOOL decoy ;
  BOOL avoidPseudoGenes ;   /* kill hits to unspliced gene on the same strand as hits to spliced genes */
  enum { STRATEGY_ZERO=0, STRATEGY_RNA_SEQ, STRATEGY_GENOME} strategy ;
  float minTM ;
  unsigned char **sl, **slR, **entryAdaptor, **entryAdaptorDecoded, **entryAdaptor2, **entryAdaptor2Decoded , **exitAdaptor, **exitAdaptor2, **exitAdaptorDecoded, **exitAdaptor2Decoded ;
  unsigned int *entryAdaptorFound, *entryAdaptor2Found,  *exitAdaptorFound,  *exitAdaptor2Found, *entryAdaptorHasN, *entryAdaptor2HasN,  *exitAdaptorHasN ,  *exitAdaptor2HasN ;
  int targetMutantMalus, mutantMalus, geneStrandMalus ;
  int aliBonusForThisAli, leftAdaptorClipForThisAli, rightAdaptorClipForThisAli ; /* belong to mm, but it saves memory to have them as globals */   
  int errMin, errMax, errRateMax, errCost, endBonus, bonus, slBonus, clipN,
    maxHit, maxHitR, minEntropy, nExportedHits,
    probeMinMultiplicity,
    probeMinLength, probeMaxLength,  /* requested by the user */
    probeLengthMin, probeLengthMax,  /* observed in the prore file */
    errPosMin, errPosMax, minAli, minAliPerCent, minChainAli, 
    seedOffset, seedShift, seedOffsetUsed, seedLength, seedLengthRequired ;
  int decimate ;
  DICT *previousScoreDict, *rejectedProbeDict, *selectedProbeDict, *exportDict, *targetDict, *target2geneDict, *probeDict ;
  BitSet isSplicedTarget, scannedGenes ;
  Associator adaptorAss ;
  Array probePairs ;
  Array originalScores ;
  Array target2geneArray ;
  Array probes, letters, dna0, previousScores ;
  BigArray exportGeneHits ;
  BigArray exportHits, exportDonors, exportAcceptors ;
  BigArray exportIntrons ;
  BigArray exportDoubleIntrons ;
  BigArray geneDonors, geneAcceptors ;
  BigArray geneIntrons ;
  Associator ass, ass1, ass2, assB, assB1, assB2, assMulti[32] ;
  Array scores[1024] ;
  int hashPhase, nOligo, score ;
  Array targetMask ;
  Array oligoFrequency ;
  Array oligoFrequencyB ;
  BOOL doubleSeed ;
  BOOL isShortTarget ;
  BOOL gzi, gzo ;
  BitSet seedOver25 ;
  BitSet doubleGenes ;
  ACEIN ai ; ACEOUT ao, outNoInsert, aoRead2Intron ;

  const char *method ;
  long int nRawReads, nRawBases, nReads, nErrors ;
  long int nAlignedPairs, nAlignedReads, nUnaligned, nMultiAligned[11], nAlignedBases ;
} CLIPALIGN ;

#define NBESTSCORES 25

typedef struct mmStruct {
  int probeLength ;
  int probeLeftClip ;
  int probeRightClip ;
  int probeName ;     /* offsett in CLIPALIGN->probeStack */
  int probeSetName ;  /* offsett in CLIPALIGN->probeStack */
  int probeSequence ; /* offsett in CLIPALIGN->probeStack */
  int probeQuality ;  /* offsett in CLIPALIGN->qualityStack */
  int nScore ;        /* offset in scores array for that loop */
  int lastTarget, lastTargetGene, lastA1, lastA2, lastX1, lastX2, lastDelta ;
  int ipx ;           /* offset in exportHits array */
  int ipxg ;          /* offset in exportGeneHits array */
  int pair ;          /* +- offset + 2 of other read of same fragment */
  int fragment ;      /* fragment number */
  int mult ;          /* multiplicity in fastc format */
  int entropy ;       /* complexity in bp equivalent */
  BOOL lastIsUp ;
  float tm ;          /* TM Maniatis formula */
  short nGeneHits ;   /* number of hits reported */
  short nPerfectHit ; /* number of perfect hits reported */
  short bestHit ;     /* minimal N +err+unaligned so far */
  short bestScore[NBESTSCORES], bestX1[NBESTSCORES], bestX2[NBESTSCORES] ;     /* best score so far */
  short pairScore ;
  short bestPairScore ;
  short same ;          /* number of letters in common with (mm+1) */ 
  short stranded ;      /* 0, 1 or -1, depends on run and if paired run on probeName prefix */
  unsigned char wordCount[16] ; /* count how many time sequence section modulo 8 got a word hit */
  int gt_ag, ct_ac ;
} MM ;


static char B2[256] ;
static int nVerif = 0 ;
static int nPerfect = 0 ;
static int nHashRejected = 0 ;

static BOOL is64 = FALSE ;
#define INTRONBONUS 3
static BOOL slideCheckExonEntropy16 (CLIPALIGN *pp, unsigned const char *exon) ;
static BOOL slideCheckHookEntropy (const unsigned char *hook) ;
static int clipAlignOptimizePairs (CLIPALIGN *pp, BOOL singleGene) ;

static const char *COSMID = 0 ;
static void clipAlignSquish (CLIPALIGN *pp, int nn) ;

/*************************************************************************************/
/*********************************************************************/
typedef struct pExportStruct { 
  int target, targetGene
    , ln, ali, score, pairScore, chain, chainAli
    , a1, a2, x1, x2, s1, s2, uu
    , presuffix, prefix, suffix
    , targetPresuffix 
    , nN, nErr
    , errLeftVal, errRightVal
    , ipxold
    /*     , geneStrandMalus */
    , fragment
    , intron, doubleIntron
    , donor, acceptor
    ;
  MM *mm ;
  char type, isDown, foot[6]
    ; } PEXPORT ;

/*************************************************************************************/

static void showPexport (CLIPALIGN *pp, int nn)
{
  int ii ;
  BigArray bb = 0 ;
  int iMax ;
  PEXPORT *px ;
  const char *title = "" ;

  switch (nn)
    {
    case 1: bb = pp->exportGeneHits ; title = "geneHits" ; break ;
    case 2: bb = pp->geneDonors ; title = "gene donors" ; break ;
    case 3: bb = pp->geneAcceptors ; title = "gene acceptors" ; break ;
    case 4: bb = pp->exportDonors ; title = "export donors" ; break ;
    case 5: bb = pp->exportAcceptors ; title = "export acceptors" ; break ;
    case 6: bb = pp->geneIntrons ; title = "gene introns" ; break ;
    case 7: bb = pp->exportIntrons ; title = "export introns" ; break ;
    case 9: bb = pp->exportDoubleIntrons ; title = "export double introns" ; break ;
    default: bb = pp->exportHits ; title = "hits" ; break ;
    }

  if (bb)
    {
      iMax = bigArrayMax (bb) ;
      fprintf (stderr, "---  showPexport %d lines, %s\n", iMax, title) ;
      if (iMax)
	for (ii = 0, px = bigArrp (bb, 0, PEXPORT) ; ii < iMax ; ii++, px++)
	  if (px->mm && (nn > 5 ||  (px->target && px->score))) 
	    fprintf (stderr, "%d: \ts:%d\tps=%d\tali:%d\ta:%d-%d \tx:%d-%d\t%s\tfrag=%d\tchain=%d:%d:%d chAli=%d\tu=%d\t%s\tpfx:%d\tsfx:%d\ti:%d di:%d d:%d a:%d\n"
		     , ii
		     , px->score,  px->pairScore, px->ali
		     , px->a1, px->a2
		     , px->x1, px->x2
		     , stackText (pp->probeStack, px->mm->probeName)
		     , px->fragment
		     , px->chain, px->s1, px->s2, px->chainAli, px->uu
		     , px->target ? dictName (pp->targetDict, px->target) : "---"
		     , px->prefix, px->suffix
		     , px->intron, px->doubleIntron, px->donor, px->acceptor
		     ) ;
    }
} /* showPexport */

/*************************************************************************************/

static void clipAlignSquish (CLIPALIGN *pp, int nn)
{
  int ii, jj ;
  BigArray aa = 0 ;
  int iMax ;
  PEXPORT *px, *py, *pz ;

  switch (nn)
    {
    case 1: aa = pp->exportGeneHits ; break ;
    case 2: aa = pp->geneDonors ; break ;
    case 3: aa = pp->geneAcceptors ; break ;
    case 4: aa = pp->exportDonors ; break ;
    case 5: aa = pp->exportAcceptors ; break ;
    default: aa = pp->exportHits ; break ;
    }
  iMax = aa ? arrayMax (aa) : 0 ;
  if (iMax)
    {
      for (ii = jj = 0, px = py = pz = bigArrp (aa, 0, PEXPORT) ; ii < iMax ; ii++, px++)
	{	
	  if (px->score)
	    {
	      if (jj < ii && ! memcmp (px, pz, sizeof (PEXPORT)))
		continue ;
	      if (jj < ii) *py = *px ;
	      pz = py ;
	      py++; jj++ ;
	    }
	}
      bigArrayMax (aa) = jj ;
    }
}  /* clipAlignSquish */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/*****************************  actual work ******************************************/
/*************************************************************************************/
/* Compute the product of 2 Solid transitions 
 * It gives at the same time a product transition and the  value of a transformed base
 */
static unsigned char clipAlignSolidProduct (unsigned char a1, unsigned char a2)
{
  unsigned char x1 = 0 ;
  switch ((int)a1)
    {
    case 0: 
    case A_: x1 = a2 ; break ;
    case T_:
      switch ((int)a2)
	{
	case 0: 
	case A_: x1 = T_ ; break ;
	case T_: x1 = A_ ; break ;
	case G_: x1 = C_ ; break ;
	case C_: x1 = G_ ; break ;
	}
      break ;
    case G_:
      switch ((int)a2)
	{
	case 0: 
	case A_: x1 = G_ ; break ;
	case T_: x1 = C_ ; break ;
	case G_: x1 = A_ ; break ;
	case C_: x1 = T_ ; break ;
	}
      break ;
    case C_:
      switch ((int)a2)
	{
	case 0: 
	case A_: x1 = C_ ; break ;
	case T_: x1 = G_ ; break ;
	case G_: x1 = T_ ; break ;
	case C_: x1 = A_ ; break ;
	}
      break ;
    }
  return x1 ;
} /* clipAlignSolidProduct */

/*************************************************************************************/
/* verify that the transition products a1.a2 and b1.b2 are equal
 * so the letter before the double error does transforms into the letter after the double error
 * if yes, we have in reality a single error in the shortDna relative to the longDna
 */
static BOOL clipAlignSolidCheckError (unsigned char a1, unsigned char a2, unsigned char b1, unsigned char b2)
{
  unsigned char x1 = clipAlignSolidProduct (a1, a2) ;
  unsigned char x2 = clipAlignSolidProduct (b1, b2) ;

  return (x1 != 0 && x1 == x2) ? TRUE : FALSE ;
} /* clipAlignSolidCheckError */

/*************************************************************************************/

static int clipAlignSolidErrorCorrect (MM *mm, Array dnaShort, Array dnaLong, Array err, int x1, int x2, BOOL isUp, int *nCorrectedp)
{
  int KK = 1 ;
  int oldx, ii, jj, nn = 0 ;
  unsigned char cc1, cc2, cc3, cc4, cca, ccb, ccc ;
  A_ERR *ep1, *ep2 ;
  BOOL found ;

  *nCorrectedp = 0 ;
  if (arrayMax (err) == 0)
    return 0 ;

  /* eliminate isolated errors if the next KK bp are correct */
  for (ii = oldx = 0, found = FALSE ; ii < arrayMax(err) ; ii++)
    {
      ep1 = arrp (err, ii, A_ERR) ;
      if (ep1->type == AMBIGUE) { ep1->type = TYPE80 ; continue ; }
      if (ep1->type == ERREUR && ep1->iShort > x1 &&  ep1->iShort < x2)
	{
	  if (ep1->iShort > oldx + KK && (ii == 0 || (ep1-1)->type == ERREUR ||  (ep1-1)->type == AMBIGUE) &&
	      (
	       (ii == arrayMax(err) - 1 && ep1->iShort < arrayMax(dnaShort) - KK) ||
	       (ii < arrayMax(err) - 1 && (ep1+1)->type == ERREUR && ep1->iShort < (ep1+1)->iShort - KK)
	       )
	      )
	    { 
	      if (0) { found = TRUE ; ep1->type = 0 ; } 
	      else { ep1->type = TYPE80 ; nn++ ;}
	    }
	}
      oldx = ep1->iShort ;
    }	  
  /* compress */
  if (found)
    {
      for (ii = jj = 0, ep1 = ep2 = arrp (err, 0, A_ERR) ; ii < arrayMax(err) ; ep1++, ii++)
	{	  
	  if (ep1->type)
	    {
	      if (ii > jj) *ep2 = *ep1 ;
	      jj++ ; ep2++ ;
	    }
	  else
	    {
	      if (ep1->iShort + mm->probeLeftClip > 1)
		nn++ ; /* this error has been error corrected (but at position 2 it does not count */
	    }
	}
      arrayMax (err) = jj ;
    }
  *nCorrectedp = nn ;
  for (ii = 0 ; ii < arrayMax(err) ; ii++)
    {
      ep1 = arrp (err, ii, A_ERR) ;
      if (ep1->type && ep1->type != TYPE80 && ep1->iShort > x1 &&  ep1->iShort < x2)
	{
	  BOOL ok = FALSE ;

	  for (ep2 = 0, jj = 1 ; ! ep2 && ii + jj < arrayMax(err) ; jj++)
	    if ((ep1+jj)->type != TYPE80) ep2 = ep1 + jj ;
	  
	  if (ep2 && ep2->iShort == ep1->iShort)
	    ok = TRUE ;
	  else if (ep2 && ep2->iShort == ep1->iShort + 1)
	    ok = TRUE ;
	  else if (ep2 && ep2->iShort == ep1->iShort + 2 && ep1->type == INSERTION_DOUBLE)
	    ok = TRUE ;
	  else if (ep2 && ep2->iShort == ep1->iShort + 3 && ep1->type == INSERTION_TRIPLE)
	    ok = TRUE ;
	  else if (ep2 && ep1->type == ERREUR && ep2->type == INSERTION)
	    {
	      int i ;
	      unsigned char cc = ep2->baseShort ;
	      
	      ok = TRUE ;
	      for (i = ep1->iLong + 1 ; ok && i < ep2->iLong ; i++)
		if (arr (dnaLong, i, unsigned char) != cc) ok = FALSE ;
	    }
	  else if (ep1->type == INSERTION && (!ep2 || ep2->type != ERREUR || ep2->iShort != ep1->iShort + 1))
	    {
	      ok =  TRUE ; ep2 = 0 ;
	    }
	  else if (ep1->type == INSERTION_DOUBLE && (!ep2 || ep2->type != ERREUR || ep2->iShort != ep1->iShort + 2))
	    {
	      ok =  TRUE ; ep2 = 0 ;
	    }
	  else if (ep1->type == INSERTION_TRIPLE && (!ep2 || ep2->type != ERREUR || ep2->iShort != ep1->iShort + 3))
	    {
	      ok =  TRUE ; ep2 = 0 ;
	    }
	  else if (ep2 && ep1->type == TROU && ep2->type == ERREUR && ep1->iLong + 1 == ep2->iLong)
	    {
	      ok = TRUE ;
	    }
	  else if (ep1->type == TROU && (!ep2 || ep1->iLong + 1 != ep2->iLong))
	    {
	      ok = TRUE ; ep2 = 0 ;
	    }
	  else if (ep2 && ep1->type == TROU_DOUBLE && ep2->type == ERREUR && ep1->iLong + 2 == ep2->iLong)
	    {
	      ok = TRUE ;
	    }
	  else if (ep1->type == TROU_DOUBLE && (!ep2 || ep1->iLong + 2 != ep2->iLong))
	    {
	      ok = TRUE ; ep2 = 0 ;
	    }
	  else if (ep2 && ep1->type == TROU_TRIPLE && ep2->type == ERREUR && ep1->iLong + 3 == ep2->iLong)
	    {
	      ok = TRUE ;
	    }
	  else if (ep1->type == TROU_TRIPLE && (!ep2 || ep1->iLong + 3 != ep2->iLong))
	    {
	      ok = TRUE ; ep2 = 0 ;
	    }
	  else if (ep2 && ep1->type == ERREUR && ep2->type == TROU)
	    {
	      int i ;
	      unsigned char cc = arr (dnaLong, ep2->iLong, unsigned char) ;
	      
	      ok = TRUE ;
	      for (i = ep1->iLong + 1 ; ok && i < ep2->iLong ; i++)
		if (arr (dnaLong, i, unsigned char) != cc) ok = FALSE ;
	      for (i = ep1->iShort + 1 ; ok && i < ep2->iShort ; i++)
		if (arr (dnaShort, i, unsigned char) != cc) ok = FALSE ;
	    }
	  if (ok) /* 2 consecutive errors may count as 1 */
	    {
	      cc1 = cc2 = cc3 = cc4 = 0 ;
	      /* verify that the double transition on the genome == double transition on the cdna */
	      if (ep1->type == ERREUR) 
		{
		  cc1 = arr(dnaShort, ep1->iShort, unsigned char) ;
		  cc3 = arr(dnaLong, ep1->iLong, unsigned char) ;
		}
	      else if (ep1->type == INSERTION)
		{
		  cc1 = arr(dnaShort, ep1->iShort, unsigned char) ;
		  cc2 = cc3 = cc4 = A_ ;
		}
	      else if (ep1->type == INSERTION_DOUBLE)
		{
		  cca = arr(dnaShort, ep1->iShort, unsigned char) ;
		  ccb = arr(dnaShort, ep1->iShort+1, unsigned char) ;
		  cc1 = clipAlignSolidProduct (cca, ccb) ;
		  cc2 = cc3 = cc4 = A_ ;
		}
	      else if (ep1->type == INSERTION_TRIPLE)
		{
		  cca = arr(dnaShort, ep1->iShort, unsigned char) ;
		  ccb = arr(dnaShort, ep1->iShort+1, unsigned char) ;
		  ccc = arr(dnaShort, ep1->iShort+2, unsigned char) ;
		  cc1 = clipAlignSolidProduct (cca, ccb) ;
		  cc1 = clipAlignSolidProduct (cc1, ccc) ;
		  cc2 = cc3 = cc4 = A_ ;
		}
	      else if (ep1->type == TROU)
		{
		  cc1 = cc2 = cc4 = A_ ;
		  cc3 = arr(dnaLong, ep1->iLong, unsigned char) ;
		}
	      else if (ep1->type == TROU_DOUBLE)
		{
		  cca = arr(dnaLong, ep1->iLong, unsigned char) ;
		  ccb = arr(dnaLong, ep1->iLong + 1, unsigned char) ;
		  cc2 = clipAlignSolidProduct (cca, ccb) ;
		  cc1 = cc3 = cc4 = A_ ;
		}
	      else if (ep1->type == TROU_TRIPLE)
		{
		  cc1 = A_ ;
		  cca = arr(dnaLong, ep1->iLong, unsigned char) ;
		  ccb = arr(dnaLong, ep1->iLong + 1, unsigned char) ;
		  ccc = arr(dnaLong, ep1->iLong + 2, unsigned char) ;
		  cc3 = clipAlignSolidProduct (cca, ccb) ;
		  cc3 = clipAlignSolidProduct (cc3, ccc) ;
		}
	      else 
		ok = FALSE ;

	      if (ep2 && ep2->type == ERREUR)
		{
		  cc2 = arr(dnaShort, ep2->iShort, unsigned char) ;
		  cc4 = arr(dnaLong, ep2->iLong, unsigned char) ;
		}
	      else if (ep2 && ep2->type == INSERTION)
		{
		  cc2 = arr(dnaShort, ep2->iShort, unsigned char) ;
		  cc4 = A_ ;
		}
	      else if (ep2 && ep2->type == TROU)
		{
		  cc2 = A_ ;
		  cc4 = arr(dnaLong, ep2->iLong, unsigned char) ;
		}
	    }
	  
	  if (ok && clipAlignSolidCheckError (cc1, cc2, cc3, cc4))
	    {
	      if (ep1->type == ERREUR) 
		ep1->type = ep2->type ; /* the idea is to still recognize indels */
	      if (ep1->type == TROU && ! isUp) 
		ep1->iLong-- ;
	      if (ep2) 
		{
		  ep2->type = 0 ;
		  ep1++ ; ii++ ; /* the error at ep2 is already treated */
		}
	    }
	  else  /* non acceptable error : blank them out */
	    { 
	      if (1)  /* clip on this error */
		{ 
		  ep1->type = TYPE80 ; 
		  (*nCorrectedp)++ ;
		  if (0 && ep2) 
		    {
		      ep2->type = TYPE80 ;
		      (*nCorrectedp)++ ;
		      ep1++ ; ii++ ; /* the error at ep2 is already treated */
		    }
		}
	    }     
	}
    }
  /* eliminate an eventual isolated last error */
  if (ii == arrayMax(err) - 1)
    {
      ep1 = arrp (err, ii, A_ERR) ;
      ep1->type = TYPE80 ;   
      (*nCorrectedp) += 1 ;
    }
  /* compress and count */ 
  nn = 0 ;
  for (ii = jj = 0, ep1 = ep2 = arrp (err, 0, A_ERR) ; ii < arrayMax(err) ; ep1++, ii++)
    {
      if (ep1->type == TROU && ep1->sens == 1)
	ep1->iLong-- ; /* isolated repeat */
      if (ep1->type && ep1->iShort > 0)
	{
	  if (ii > jj) *ep2 = *ep1 ;
	  if (ep2->type == TYPE80) nn++ ;
	  jj++ ; ep2++ ;
	}
      /* else  this error has been error corrected but does not count, it was part of a double error */
    }
  arrayMax (err) = jj ;
  (*nCorrectedp) = nn ;

  return nn ; /* number of corrected errors */
} /* clipAlignSolidErrorCorrect */

/*************************************************************************************/

static void clipAlignErrorDoClip (CLIPALIGN *pp, MM *mm, BOOL transitionSpace, Array err
				  , int *a1p, int *a2p,  int *x1p, int *x2p
				  , int errCost, int transitionCost, int transitionFavor, int endBonus
				  , Array dna, Array probeDnaArray
				  )
{
  int ii, jj = 0, i0, i1, i2, oldx = 0, dx, eMax ;
  int score, oldErrCost ;
  int ambErrCost = 2 ; 
  A_ERR *ep1, *ep2 ;

  /* never count the first base */
  /* locate the longuest exact strecth */
  i1 = 0; i0 = -1; oldx = *x1p-1; dx = 0 ; eMax = arrayMax (err) ; 
  if (eMax)
    for (ep1 = arrp (err, i1, A_ERR) ; i1 < eMax && ep1->iShort < *x2p+1 ; ep1++, i1++)
    {
      if (ep1->iShort - oldx > dx)
	{ dx = ep1->iShort - oldx ; i0 = i1 ; }
      oldx = ep1->iShort ;
    }
  if (*x2p  - oldx > dx)
    { dx = *x2p - oldx ; i0 = i1 ; }
  
  /* eliminate error walls  to the right by computing the running score */
  if (i0 < arrayMax (err))
    {
      oldErrCost = errCost ;
      for (i2 = ii = i0, score = oldErrCost, ep2 = arrp (err, ii, A_ERR), oldx = ep2->iShort, eMax = arrayMax (err) ;
	   ii < eMax ; ep2++, ii++)
	{
	  score += ep2->iShort - oldx - oldErrCost ; /* start at zero when ii = i2 = i0 */ 
	  if (score + transitionFavor > 0) /* if the score is positive i2 is the error to be clipped */
	    {i2 = ii ; score = 0 ;}
	  oldx = ep2->iShort ;
	  oldErrCost = errCost ;
	  switch (ep2->type) /* faux ? , on veut sans doute oldErrCost */
	    {
	    case ERREUR: 
	      if (ep2->baseLong & ep2->baseShort)
		oldErrCost = 1 ;  
	      if (0)
		switch ((int)ep2->baseLong)
		  {
		  case A_:
		  case T_:
		  case G_:
		  case C_:
		    break ;
		  default:
		    oldErrCost = 1 ;  
		    break ;
		  }
	      break ;
	    case TROU_DOUBLE: oldErrCost += 2 ; break ; 
	    case TROU_TRIPLE: oldErrCost += 4 ; break ; 
	    case INSERTION_DOUBLE: oldErrCost += 2 ; break ; 
	    case INSERTION_TRIPLE: oldErrCost += 4 ; break ; 
	    case TYPE80: oldErrCost = transitionCost ; break ;
	    case AMBIGUE: oldErrCost = ambErrCost ; break ;
	    default: break ;
	    }
	  if (1 && ii == eMax - 1) /* if the last error is an indel, it cannot must be followed by a homopolymer */ 
	    switch (ep2->type) 
	      {
	      case TROU :
	      case INSERTION:
	      case TROU_DOUBLE: 
	      case TROU_TRIPLE: 
	      case INSERTION_DOUBLE: 
	      case INSERTION_TRIPLE: 
		{
		  int i5 = ep2->iShort + 1 ;
		  int i4 = *x2p - i5 ;
		  unsigned char b,  *cp = arrp (probeDnaArray, i5, unsigned char) ;
		  b = *cp ;
		  while (i4 && (*cp++ == b)) i4-- ;
		  if (i4 == 0)
		    oldErrCost = 10000 ;
		  break ;
		}
	      default:
		break ;
	      }
	}
      score += *x2p - oldx - oldErrCost ;
      if (*x2p == mm->probeLength) score += endBonus ;
      if (score + transitionFavor > 0) /* if the score is positive do not clip any error to the right */
	i2 = ii ;
    }
  else
    i2 = i0 ;

  /* eliminate error walls  to the left by computing the running score */
  if (i0 > 0)
    {
      oldErrCost = errCost ;
      for (i1 = ii = i0 - 1, score = oldErrCost, ep1 = arrp (err, ii, A_ERR), oldx = ep1->iShort ; ii >= 0 ; ep1--, ii--)
	{
	  score += - ep1->iShort + oldx - oldErrCost ; /* start at zero when ii = i1 = i0 - 1 */ 
	  if (score + transitionFavor > 0) /* if the score is positive i2 is the error to be clipped */
	    { 
	      i1 = ii ; score = 0 ; 
	    }
	  oldx = ep1->iShort ;	      oldErrCost = errCost ;
	  switch (ep1->type) /* faux, on veut sans doute oldErrCost */
	    {
	    case TROU_DOUBLE: oldErrCost += 2 ; break ; 
	    case TROU_TRIPLE: oldErrCost += 4 ; break ; 
	    case INSERTION_DOUBLE: oldErrCost += 2 ; break ; 
	    case INSERTION_TRIPLE: oldErrCost += 4 ; break ; 
	    case TYPE80: oldErrCost = transitionCost ; break ;
	    case AMBIGUE: oldErrCost = ambErrCost ; break ;
	    default: break ;
	    }
	  if (1 && ii == 0) /* if the last error is an indel, it cannot must be followed by a homopolymer */ 
	    switch (ep1->type) 
	      {
	      case TROU :
	      case INSERTION:
	      case TROU_DOUBLE: 
	      case TROU_TRIPLE: 
	      case INSERTION_DOUBLE: 
	      case INSERTION_TRIPLE: 
		{
		  int i5 = ep1->iShort ;
		  int i4 = i5 - *x1p + 1 ;
		  if (i5 >= 1 && i4  > 0)
		    {
		      unsigned char b,  *cp = arrp (probeDnaArray, i5 - 1, unsigned char) ;
		      b = *cp ;
		      while (i4 && (*cp--  == b)) i4-- ;
		      if (i4 == 0)
			oldErrCost = 10000 ;
		    }
		  break ;
		}
	      default:
		break ;
	      }

	}
      score +=  oldx + 2 - *x1p - oldErrCost ;
      if (*x1p == 1) score += endBonus ;
      if (score + transitionFavor > 0) /* if the score is positive do not clip any error to the left */
	i1 = ii ;
    }
  else
    i1 = -1 ;

  if (i1 >= 0 && i1 < arrayMax (err))
    {
      ep1 = arrp (err, i1, A_ERR) ;
      if (transitionSpace && ep1->type == ERREUR && ep1->iShort == 0) ;
      else
	{
	  *x1p = ep1->iShort + 2 ;  /* not  ep1->iShort +1,  Plato */ 
	  *a1p = ep1->iLong  + 2 ;
	  
	  if (transitionSpace)
	    {
	      if (ep1->sens == 10 && ep1->type == ERREUR) /* 2010_05_23: because or error correction */
		{ (*x1p)++ ; (*a1p)++ ; }
	      if (ep1->type == TROU)                      /* 2010_05_23: because or error correction */
		{ (*x1p)+=2 ; (*a1p)+=2 ; }
	      if (ep1->sens == 10 && ep1->type == INSERTION) /* 2010_05_23: because or error correction */
		{ (*x1p)+=2 ; (*a1p)+=2 ; }
	    }

	  switch (ep1->type)
	    {
	    case TROU: (*x1p)-- ; break ; 
	    case TROU_DOUBLE: (*x1p) -=1 ; (*a1p)++ ; break ; 
	    case TROU_TRIPLE: (*x1p) -=1 ; (*a1p) += 2 ; break ; 
	    case INSERTION: (*a1p)-- ; break ; 
	    case INSERTION_DOUBLE:(*a1p)-- ; (*x1p)++ ; break ; 
	    case INSERTION_TRIPLE:(*a1p)-- ; (*x1p)+=2 ; break ; 
	    case AMBIGUE:
	    case ERREUR:
	    default: break ;
	    }
	}
    }
  if (i2 >= 0 && i2 < arrayMax (err))
    {
      ep2 = arrp (err, i2, A_ERR) ;
      *x2p = ep2->iShort ; /* not  ep2->iShort -1,  Plato */ 
      *a2p = ep2->iLong  ;
    }

  /* compress */ 

  eMax = arrayMax (err) ; jj = 0 ;
  if (eMax)
    for (ii = jj = 0, ep1 = ep2 = arrp (err, 0, A_ERR) ; ii < eMax ; ep1++, ii++)
      {
	if (0)
	  {
	    if (ep1->sens > 1) ep1->sens = 1 ;
	    if (ep1->sens < -1) ep1->sens = -1 ;
	  }
	if (ii > i1 && ii < i2)
	  {
	  if (ii > jj) *ep2 = *ep1 ;
	  jj++ ; ep2++ ;
	  }
      }
  /* reclip the A2G type errors close to the edges
   *  but NOT recursivelly
   */

 if (0 && pp->X2X)
    {
      char *cp, *cq ;
      int i, j, k1 = 0, k2 = 0 ;
      for (i = 1, j = 0, cp = arrp (dna, *a1p - 1, char), cq = arrp (probeDnaArray, *x1p - 1, char) ; j < errCost && j <= *x2p - *x1p ; i++, j++, cp++, cq++) 
	{
	  if (*cp != *cq)
	    {
	      if (*cp & *cq)
		k1 = i ;
	    }
	}
      for (i = 0, j = 0, cp = arrp (dna, *a2p -errCost, char), cq = arrp (probeDnaArray, *x2p - errCost, char) ; j < errCost - 1 ; i++, j++, cp++, cq++) 
	{
	  if (*cp != *cq)
	    {
	      if (*cp & *cq)
		{
		  k2 = errCost - i ;
		  break ;
		}
	    }
	}
      if (k1)
	{ *x1p += k1 ; *a1p += k1 ; }
      if (k2)
	{ *x2p -= k2 ; *a2p -= k2 ; }
    }
  
 if (0)
   fprintf (stderr ,"max = %d  jj=%d\n", arrayMax(err), jj);

 if (jj >=0 && jj < arrayMax (err))
   arrayMax (err) = jj ;
} /* clipAlignErrorDoClip */

/*************************************************************************************/
/* eliminate the error walls at the 2 ends */
static int clipAlignErrorClip (CLIPALIGN *pp, MM *mm, 
			       Array err, Array probeDnaArray, Array dna, Array dnaR,
			       int *nNp, int *nCorrectedp, int pos, int xpos, int *a1p, int *a2p, int *x1p, int *x2p,
			       BOOL isUp, float *tmp)
{
  int ii, jj, dx, eMax ;
  int score, score2, errCost, oldErrCost ;
  int nCorrected = 0 ;
  int pairTampon = pp->pairTampon ;
  int endBonus = pp->endBonus ;
  A_ERR *ep1 ;
  char *cp ;
  int ooCost ;           /* cost of a corrected base reverted to fasta space */
  int transitionCost ;   /* cost of correcting isolated errors in color space */
  int transitionFavor ;  /* if errCost = 4 we accept when the running score reaches 1
			  * if 8, if it reaches 0, this is implemented by setting favor = 1
			  */
  /* elimininates stretches of n */
  cp = arrp (dna, *a1p - 1, char) ; dx = *a2p - *a1p + 1 ;
  while (*cp++ == N_ && dx > 0) 
    { dx-- ; (*a1p)++ ; (*x1p)++ ; }

  /* Notice we run this condition BEFORE solid error correction */
  if (pp->exactTargetStart > 0 && pp->exactTargetStop > 0)
    {
      int exact1 = 0, exact2 = 0 ;
      /* these numbers are given in bio coordinates on the forward strand of the target */
      if (!isUp)
	{
	  exact1 = pp->exactTargetStart ;
	  exact2 =  pp->exactTargetStop ;
	}
      else
	{
	  exact1 = arrayMax(dna) - pp->exactTargetStop + 1 ;
	  exact2 = arrayMax(dna) - pp->exactTargetStart + 1 ;
	}
      if (*a1p > exact1 || *a2p < exact2)
	return 0 ;
      if ((eMax = arrayMax (err)) > 0)
	for (ii = 0, ep1 = arrp (err, 0, A_ERR) ; ii < eMax ; ii++, ep1++)
	  if (ep1->iLong + 1 >= exact1 && ep1->iLong + 1 <= exact2)
	    return 0 ;
    }

  if (pp->minAli > *x2p - *x1p + 1)
    return 0 ;
    
  switch (pp->strategy)
    {
    case STRATEGY_RNA_SEQ:
      ooCost = 8 ;
      transitionCost = 3 ;
      errCost = pp->errCost ;
      transitionFavor = 1 ;
      break ;
    case STRATEGY_GENOME:
      ooCost = 3 ;
      transitionCost = 2 ;
      errCost = pp->errCost ;
      transitionFavor = 1 ;  
      /* FDA 2017_08  was 4 , see also STRATEGY_ in snp.c */
       break ;
    default:
      ooCost = 0 ;
      errCost = pp->errCost ;
      transitionCost = pp->errCost ;
      transitionFavor = 1 ;
      break ;
    }
  

  /* Solid error correcting code */
  if (arrayMax (err) && pp->solid) 
    {
      /* first hard clip left and right walls */
      clipAlignErrorDoClip (pp, mm, TRUE, err, a1p, a2p, x1p, x2p, transitionCost, transitionCost, transitionFavor, endBonus, dna, probeDnaArray) ; /* transition space */
      /* do not abandon non correctable errors before clipping on the walls */
      clipAlignSolidErrorCorrect (mm, probeDnaArray, dna, err, *x1p, *x2p, isUp, &nCorrected) ;
      *nNp = 0 ;
    }

  if (pp->minAli && (pp->X2X || arrayMax (err) > 0))
    {
      /* eliminate the error walls at the 2 ends */
      clipAlignErrorDoClip (pp, mm, FALSE, err, a1p, a2p, x1p, x2p, errCost, transitionCost, 0, endBonus, dna, probeDnaArray) ; /* fasta space costs */
    }

  if (pp->minAli > *x2p - *x1p + 1)
    return 0 ;

  /* recompute the errors from left to rigth in the correct zone to rationalize
   * the previous left-right error tracking
   */
  

  if (pp->X2X || arrayMax (err))
    {
      if (*x1p < 1) *x1p = 1 ;
      if (*x2p > mm->probeLength)
	*x2p = mm->probeLength ;
      
      if (isUp == FALSE)
	{    
	  int pass = 5 ;
	  while (pass-- > 0) 
	    { /* beware of losing the zone of the anchor */
	      aceDnaDoubleTrackErrors (probeDnaArray, x1p, x2p, TRUE /* bio coordinates  */
				       , dna, dnaR, a1p, a2p
				       , nNp, err, MAXJUMP, pp->errMax > 0 ? 999 : 0, FALSE, 0) ;
	      if (xpos >= *x1p && xpos <= *x2p)
		break ;
	      /* rescan */
	      *x1p = xpos ; *a1p = pos ; 
	      *x2p = xpos + pp->seedLength ; *a2p = pos + pp->seedLength ;
	      aceDnaDoubleTrackErrors (probeDnaArray, x1p, x2p, TRUE /* bio coordinates  */
				       , dna, dnaR, a1p, a2p
				       , nNp, err, 2, pp->errMax > 0 ? 999 : 0, TRUE, 0) ;
	    }

	  if (pp->solid) 
	    {
	      /* we do not need to 	 clipAlignErrorDoClip(), in transition space, it is done already */ 
	      if (clipAlignSolidErrorCorrect (mm, probeDnaArray, dna, err, *x1p, *x2p, isUp, &nCorrected) == -1)
		return 0 ; /* abandon this sequence, it cannot be corrected */
	      *nNp = 0 ;
	      if (*x1p == 0)
		{ (*x1p)++; (*a1p)++ ; }
	    }
	  if (pp->minAli && (pp->X2X ||arrayMax (err) > 0))
	    {
	      /* to stay on the safe side eliminate the error walls at the 2 ends */
	      clipAlignErrorDoClip (pp, mm, FALSE, err, a1p, a2p, x1p, x2p, errCost, transitionCost, 0, endBonus, dna, probeDnaArray) ; /* fasta space costs */
	    }
	}
      else
	{
	  int y1 = arrayMax (probeDnaArray) - *x2p + 1 ;
	  int y2 = arrayMax (probeDnaArray) - *x1p + 1 ;
	  int b1 = arrayMax (dna) - *a2p + 1 ;
	  int b2 = arrayMax (dna) - *a1p + 1 ;
	  Array probeR = dnaHandleCopy (probeDnaArray, 0) ;

	  if (pp->solid) /* reverse but do not complement */
	    {
	      char  c,  *cp = arrp(probeR,0,char), *cq = cp + arrayMax(probeR) - 1 ;
	      
	      while (cp < cq)
		{ 
		  c = *cp ; 
		  *cp++ = *cq ;
		  *cq-- = c  ;
		}
	      b1++ ; b2++ ;
	    }
	  else
	    reverseComplement (probeR) ;

	  aceDnaDoubleTrackErrors (probeR, &y1, &y2, TRUE /* bio coordinates  */
				   , dnaR, dna, &b1, &b2
				   , nNp, err, MAXJUMP, pp->errMax > 0 ? 999 : 0, FALSE, 0) ;
	  
	  if (arrayMax (err) && pp->solid) 
	    {
	      /* first hard clip left and right walls in transition space */
	      clipAlignErrorDoClip (pp, mm, TRUE, err, &b1, &b2, &y1, &y2, transitionCost, transitionCost, transitionFavor, endBonus, dnaR, probeR) ; /* transition space */
	      clipAlignSolidErrorCorrect (mm, probeR, dnaR, err, y1, y2, FALSE, &nCorrected) ;
	      *nNp = 0 ;
	    }
	  
	  clipAlignErrorDoClip (pp, mm, FALSE, err, &b1, &b2, &y1, &y2, errCost, transitionCost, 0, endBonus, dnaR, probeR) ; /* transition space */

	  *x1p = arrayMax (probeDnaArray) - y2 + 1 ;
	  *a1p = arrayMax (dna) - b2 + 1 ;
	  *x2p = arrayMax (probeDnaArray) - y1 + 1 ;
	  *a2p = arrayMax (dna) - b1 + 1 ;

	  arrayDestroy (probeR) ;
	}
    }

   if (pp->targetStart > 0 && pp->targetStop > 0)
     {
       int exact1 = 0, exact2 = 0 ;
       /* these numbers are given in bio coordinates on the forward strand of the target */
       if (!isUp)
	 {
	   exact1 = pp->targetStart ;
	   exact2 =  pp->targetStop ;
	 }
       else
	 {
	   exact1 = arrayMax(dna) - pp->targetStop + 1 ;
	   exact2 = arrayMax(dna) - pp->targetStart + 1 ;
	 }
       if (*a1p > exact1 || *a2p < exact2)
	 return 0 ;
     }

  /* compute the score */
  if (*x1p < 1) *x1p = 1 ;
  if (*x2p > mm->probeLength)
    *x2p = mm->probeLength ;
  dx = *x2p - *x1p + 1 ;  /* aligned length */
  if (*x1p > 1 && pp->solid) /* left error is transition error but good base */
    dx++ ;

  score = dx ;
  if (! arrayExists(err)) 
     messcrash ("err does not exist") ;
  if (arrayMax(err))
    for (ii = 0, ep1 = arrp (err, 0, A_ERR) ; ii < arrayMax(err) ; ep1++, ii++)
      {
	if (! ep1->type) continue ;
	oldErrCost = errCost ;
	switch (ep1->type) /* faux, on veut sans doute oldErrCost */
	  {
	  case TROU_DOUBLE: oldErrCost += 2 ; break ; 
	  case TROU_TRIPLE: oldErrCost += 4 ; break ; 
	  case INSERTION_DOUBLE: oldErrCost += 2 ; break ; 
	  case INSERTION_TRIPLE: oldErrCost += 4 ; break ; 
	  case TYPE80: oldErrCost = ooCost ; break ;
	  case AMBIGUE: oldErrCost = 1 ; break ; /* will count as 2 since we substract nN at the end */
	  default: break ;
	  }
	score -= oldErrCost ;      
      }
  /* remove the in range N present in the tag (am I counting the penalty twice ?) */
  if (! pp->solid && *x1p < arrayMax (probeDnaArray))
    for (*nNp = 0, ii = *x1p, cp = arrp (probeDnaArray, ii, char) ; ii < *x2p ; cp++, ii++)
      if (*cp == N_) { (*nNp)++ ; score-- ; }
  /* report the auto corrected positions */ 
  (*nCorrectedp) += nCorrected ;

  /* penalise wrong strand */
  if (mm->stranded > 0 && isUp)
    score -= pp->strandBonus ;  
  if (mm->stranded < 0 && ! isUp)
    score -= pp->strandBonus ;       
  score += pp->bonus ;

  /* reject the poor alignments */
  if (arrayMax (err) < pp->errMin)
    return 0 ;

  if (pp->minAliPerCent)
    {
      if (100 * dx < pp->minAliPerCent * mm->probeLength)
	return 0 ;
    }
  if (pp->minAli)
    {
      if (dx < pp->minAli)
	return 0 ;
    }
  else
    {
      jj = arrayMax (err) ;
      if (jj && ! pp->errMax && !pp->errMax) 
	return 0 ;
      if (pp->errMax && jj > pp->errMax) 
	return 0 ;

      jj +=  mm->probeLength - dx ;
      if (pp->errRateMax && 100 * jj > dx * pp->errRateMax)
	return 0 ;
    }
  *tmp = 0 ;
  if (pp->minTM > 0)
    {
      float tm = oligoTm (probeDnaArray, *x1p, *x2p, 0) ;
      if (tm < pp->minTM)
	return 0 ;
      *tmp = tm ;
    }

  if (pp->minEntropy)
    {
      unsigned char *cp1, *cp2, cc ;
      BOOL sOk ;

      cp1 = arrp (probeDnaArray, *x1p - 1, unsigned char) ;
      cp2 = arrp (probeDnaArray, *x2p - 1, unsigned char) + 1 ;
      cc = *cp2 ; *cp2 = 0 ;
      sOk = slideCheckExonEntropy16 (pp, cp1) ;
      *cp2 = cc ;
      if (! sOk)
	return 0 ;
    }

  if (*x1p == 1) score += endBonus ;
  if (*x2p == mm->probeLength) score += endBonus ;

  /* compare to previous alignments */
  /* we may gain points by aligning the terminal A/T or 8 by aligning to the SL */
  score2 = score + pp->slBonus + mm->probeLeftClip + mm->probeRightClip +
    ((pp->splice || pp->clipPolyT2) ? INTRONBONUS : 0) +
    ((pp->splice || pp->clipPolyA2) ? INTRONBONUS : 0) +
    pairTampon ;

  if (mm->bestX1[0])
    {
      int i ;
      for (i = 0 ; i < NBESTSCORES ; i++)
	if (mm->bestX1[i])
	  {
	    int c1 = *x1p > mm->bestX1[i] ? *x1p : mm->bestX1[i] ;
	    int c2 = *x2p < mm->bestX2[i] ? *x2p : mm->bestX2[i] ;
	    dx = 3 *(c2 - c1 + 1) ;
	    if (pp->bestHit && !pp->splice &&  mm->bestScore[i] > score2 && dx > *x2p - *x1p)
	      return 0 ;
	  }
    }
  if (pp->bestHit && bigArrayMax (pp->exportGeneHits) > 1 && mm->ipxg)
    {
      int x1 = *x1p, x2 = *x2p, a1 = *a1p, a2 = *a2p ;
      int ipxg, ddx, dx = x2 - x1 + 1 ;
      int ipxgMax = bigArrayMax (pp->exportGeneHits) ;
      PEXPORT *px ;

      for (ipxg = mm->ipxg ; ipxg > 0 && ipxg < ipxgMax ; ipxg = px ? px->ipxold : 0)
	{ 
	  px = bigArrp (pp->exportGeneHits, ipxg, PEXPORT) ;
	  if (pp->target == px->target && px->mm == mm)
	    {
	      if (x1 == px->x1 && x2 == px->x2 && 
		  (
		   (a1 == px->a1 && a2 == px->a2) ||
		   (a1 == px->a2 && a2 == px->a1)
		       )
		  )
		return 0 ;
	      if (px->score > score2)
		{
		  ddx = (px->x2 < x2 ? px->x2 : x2) - (px->x1 > x1 ? px->x1 : x1) ;
		  if (2 * ddx > dx)
			return 0 ;
		}
	    }
	  else
	    break ;
	}
    }

  return score ;
} /* clipAlignErrorClip */

/*************************************************************************************/

static void slInit (CLIPALIGN *pp)
{
  /* if you edit the SL table, also edit the corresponding table in solexa.c */
  char *slu[] = { "GGTTTAATTACCCAAGTTTGAG" ,    /* SL1 */
	 	  "GGTTTTAACCCAGTTACTCAAG" ,    /* SL2 */
	 	  "GGTTTTAACCCAGTTAACCAAG" ,    /* SL3 */
		  "GGTTTTAACCCAGTTTAACCAAG" ,    /* SL4 */
     		  "GGTTTTAACCCAGTTACCAAG" ,    /* SL5 */
 	 	  "GGTTTAAAACCCAGTTACCAAG" ,    /* SL6 */
 	 	  "GGTTTTAACCCAGTTAATTGAG" ,    /* SL7 */
 		  "GGTTTTTACCCAGTTAACCAAG" ,    /* SL8 */
		  "GGTTTATACCCAGTTAACCAAG" ,    /* SL9 */
		  "GGTTTTAACCCAAGTTAACCAAG" ,    /* SL10 */
 		  "GGTTTTAACCAGTTAACTAAG" ,    /* SL11 */
		  "GGTTTTAACCCATATAACCAAG" ,    /* SL12 */
		  /* existe pas
		   * "GTTTTTAACCCAGTTACTCAAG" ,   SL13 
 	 	   * "GGTTTTTAACCCAGTTACTCAAG" ,   SL14 
		   */
		  0 } ;
  /* char *solexaOutAdaptor = "GTAGGCAC" ;  to be extended */
  static unsigned char *sl[15] ;
  static unsigned char *slR[15] ;
  static int firstPass = 0 ;
  unsigned char *cq, *cr ;
  char **cpp, *cp ;
  int i ;
  
  if (!firstPass)
    {
      memset (sl, 0, sizeof (sl)) ;
      memset (slR, 0, sizeof (slR)) ;
      firstPass = 1 ;
      for (i = 0, cpp = slu ; *cpp ; i++, cpp++)
	{
	  sl[i] = messalloc(32) ;
	  slR[i] = messalloc(32) ;
	  cp = *cpp - 1  ;
	  cq = sl[i] ; cr = slR[i] + strlen (cp+1) - 1 ;
	  while (*++cp)
	    {
	      *cq = dnaEncodeChar [(int)*cp] ;
	      *cr-- = complementBase(*cq++) ;
	    }
	}      
      pp->sl = sl ; pp->slR = slR ;
    }
  return ;
} /* slInit */

/*************************************************************************************/
/* enter all vector words and their complement in the adaptorDict
 * to avoid using them as seeds
 */
static void adaptorMask (CLIPALIGN *pp)
{
  char *buf, *cp, *cq, *cr, cc ;
  int n , ns, i, pass ;
  unsigned int oligo ;

  if (pp->entryAdaptor || pp->exitAdaptor || pp->exitAdaptor2)
    {
      if (pp->adaptorAss) 
	ac_free (pp->adaptorAss) ;
      pp->adaptorAss = assBigCreate (10000) ;
      for (pass = 0 ; pass < 3 ; pass++)
	{
	  switch (pass)
	    {
	    case 0:
	      if (! (pp->entryAdaptor && pp->entryAdaptor[0])) continue ;
	      break ;
	    case 1:
	      if (!  (pp->exitAdaptor && pp->exitAdaptor[0])) continue ;
	      break ;
	    case 2:
	      if (!  (pp->exitAdaptor2 && pp->exitAdaptor2[0])) continue ;
	      break ;
	    }
	  
	  for (ns = 0 ; ; ns++)
	    {
	      buf = 0 ;
	      switch (pass)
		{
		case 0:
		  buf = pp->entryAdaptor[ns] ? strnew ((char *)pp->entryAdaptor[ns], 0) : 0 ;
		  break ;
		case 1:
		  buf = pp->exitAdaptor[ns] ? strnew ((char *)pp->exitAdaptor[ns], 0) : 0 ;
		  break ;
		case 2:
		  buf = pp->exitAdaptor2[ns] ? strnew ((char *)pp->exitAdaptor2[ns], 0) : 0 ;
		  break ;
	      }
	      if (! buf) 
		break ;
	      n = strlen (buf) - pp->seedLength + 1 ;
	      for (i = 0, cp = buf, cq = cp + pp->seedLength ; i < n ; cp++, cq++, i++)
		{
		  cc = *cq ; *cq = 0 ;
		  oligo = 0 ;
		  for (cr = cp ; *cr ; cr++)
		    { oligo <<= 2 ; oligo |=  B2[(int)(*cr)] ; }
		  if (oligo)
		    assInsert (pp->adaptorAss, assVoid(oligo), assVoid(1)) ;
		  *cq = cc ;
		}
	      ac_free (buf) ;
	    }	    
	}
    }
} /* adaptorMask */

/**************************/

static void entryAdaptorInit (CLIPALIGN *pp, const char *eV0)
{
  static unsigned char *eV[256], *eVD[256] ;
  static int firstPass = 1 ;
  unsigned char *cq = 0, *cqD = 0 ;
  unsigned char cc ;
  const unsigned char *cp ;
  int i, j, ns ;

  if (firstPass)
    {
      memset (eV, 0, sizeof (eV)) ;
      firstPass = 0 ;
      for (i = j = 0, cp = (const unsigned char *)eV0 ; j < 255 && *cp ; cp++)
	{
	  if (i+j == 0 || *cp == ',')
	    {
	      if (j) *cq++ = 0 ;
	      cq = eV[j] = halloc(256, pp->h) ;
	      cqD = eVD[j] = halloc(256, pp->h) ;
	      j++ ;
	      if (*cp == ',') continue ;
	    }
	  cc = dnaEncodeChar [(int)*cp] ;
	  if (!cc) 
	    messcrash ("Unrecognized char in -entryAdaptor %s", pp->entryAdaptor[0]) ;
	  if (i < 255)
	    {
	      i++ ;
	      *cq++ = cc ;
	      *cqD++ = *cp ;
	    }
	}
      if (j && cq) { *++cq = 0 ; *++cqD = 0 ;}
      for (i = 0 ; i < j ; i++)
	if (strlen ((char*)eV[i])> 15)  /* limit the entry adaptor to 15 bases */
	  { eV[i][15] = 0 ; eVD[i][15] = 0 ; }
      pp->entryAdaptor = &eV[0] ;
      pp->entryAdaptorDecoded = &eVD[0] ;
      pp->entryAdaptorFound = halloc((j+1) * sizeof(int), pp->h) ; 
      pp->entryAdaptorHasN = halloc((j+1) * sizeof(int), pp->h) ;
    } 
  
  if (pp->entryAdaptor)
    for (ns = 0 ; pp->entryAdaptor[ns] ; ns++)
      {
	unsigned char *ccp = pp->entryAdaptor[ns] - 1 ;
	int nN = 0 ;
	while (*++ccp) if (*ccp == N_) nN++ ;
	pp->entryAdaptorHasN[ns] = nN ;
      }

  return ;
} /* entryAdaptorInit */

/**************************/

static void entryAdaptor2Init (CLIPALIGN *pp, const char *eV0)
{
  static unsigned char *eV[256], *eVD[256] ;
  static int firstPass = 1 ;
  unsigned char *cq = 0, *cqD = 0 ;
  unsigned char cc ;
  const unsigned char *cp ;
  int i, j, ns ;

  if (firstPass)
    {
      memset (eV, 0, sizeof (eV)) ;
      firstPass = 0 ;
      for (i = j = 0, cp = (const unsigned char *)eV0 ; j < 255 && *cp ; cp++)
	{
	  if (i+j == 0 || *cp == ',')
	    {
	      if (j) *cq++ = 0 ;
	      cq = eV[j] = halloc(256, pp->h) ;
	      cqD = eVD[j] = halloc(256, pp->h) ;
	      j++ ;
	      if (*cp == ',') continue ;
	    }
	  cc = dnaEncodeChar [(int)*cp] ;
	  if (!cc) 
	    messcrash ("Unrecognized char in -entryAdaptor2 %s", pp->entryAdaptor2[0]) ;
	  if (i < 255)
	    {
	      i++ ;
	      *cq++ = cc ;
	      *cqD++ = *cp ;
	    }
	}
      if (j && cq) { *++cq = 0 ; *++cqD = 0 ;}

      pp->entryAdaptor2 = &eV[0] ;
      pp->entryAdaptor2Decoded = &eVD[0] ;
      pp->entryAdaptor2Found = halloc((j+1) * sizeof(int), pp->h) ; 
      pp->entryAdaptor2HasN = halloc((j+1) * sizeof(int), pp->h) ;
    } 
  
  if (pp->entryAdaptor2)
    for (ns = 0 ; pp->entryAdaptor2[ns] ; ns++)
      {
	unsigned char *ccp = pp->entryAdaptor2[ns] - 1 ;
	int nN = 0 ;
	while (*++ccp) if (*ccp == N_) nN++ ;
	pp->entryAdaptor2HasN[ns] = nN ;
      }

  return ;
} /* entryAdaptor2Init */

/**************************/

static void exitAdaptorInit (CLIPALIGN *pp, const char *eV0)
{
  static unsigned char *eV[256], *eVD[256] ;
  static int firstPass = 1 ;
  unsigned char *cq = 0, *cqD = 0 ;
  unsigned char cc ;
  const unsigned char *cp ;
  int i, j, ns ;

  if (firstPass)
    {
      memset (eV, 0, sizeof (eV)) ;
      firstPass = 0 ;
      for (i = j = 0, cp = (const unsigned char *)eV0 ; j < 255 && *cp ; cp++)
	{
	  if (i+j == 0 || *cp == ',')
	    {
	      if (j) *cq++ = 0 ;
	      cq = eV[j] = halloc(256, pp->h) ;
	      cqD = eVD[j] = halloc(256, pp->h) ;
	      j++ ;
	      if (*cp == ',') continue ;
	    }
	  cc = dnaEncodeChar [(int)*cp] ;
	  if (!cc) 
	    messcrash ("Unrecognized char in -exitAdaptor %s", pp->exitAdaptor[0]) ;
	  if (i < 253)
	    {
	      i++ ;
	      *cq++ = cc ;
	      *cqD++ = *cp ;
	    }
	}
      if (j && cq) { *++cq = 0 ; *++cqD = 0 ;}

      pp->exitAdaptor = &eV[0] ;
      pp->exitAdaptorDecoded = &eVD[0] ;
      pp->exitAdaptorFound = halloc((j+1) * sizeof(int), pp->h) ;
      pp->exitAdaptorHasN = halloc((j+1) * sizeof(int), pp->h) ;
    }  

  if (pp->exitAdaptor)
    for (ns = 0 ; pp->exitAdaptor[ns] ; ns++)
      {
	unsigned char *ccp = pp->exitAdaptor[ns] - 1 ;
	int nN = 0 ;
	while (*++ccp) if (*ccp == N_) nN++ ;
	pp->exitAdaptorHasN[ns] = nN ;
      }

  return ;
} /* exitAdaptorInit */

/*************************************************************************************/

static void exitAdaptor2Init (CLIPALIGN *pp, const char *eV0)
{
  static unsigned char *eV[256], *eVD[256] ;
  static int firstPass = 1 ;
  unsigned char *cq = 0, *cqD = 0 ;
  unsigned char cc ;
  const unsigned char *cp ;
  int i, j, ns ;

  if (firstPass)
    {
      memset (eV, 0, sizeof (eV)) ;
      firstPass = 0 ;
      for (i = j = 0, cp = (const unsigned char *)eV0 ; j < 255 && *cp ; cp++)
	{
	  if (i+j == 0 || *cp == ',')
	    {
	      if (j) *cq++ = 0 ;
	      cq = eV[j] = halloc(256, pp->h) ;
	      cqD = eVD[j] = halloc(256, pp->h) ;
	      j++ ;
	      if (*cp == ',') continue ;
	    }
	  cc = dnaEncodeChar [(int)*cp] ;
	  if (!cc) 
	    messcrash ("Unrecognized char in -exitAdaptor2 %s", pp->exitAdaptor2 ? pp->exitAdaptor2[0] : (unsigned char *)"NULL") ;
	  if (i < 253)
	    {
	      i++ ;
	      *cq++ = cc ;
	      *cqD++ = *cp ;
	    }
	}
      if (j && cq) { *++cq = 0 ; *++cqD = 0 ;}

      pp->exitAdaptor2 = &eV[0] ;
      pp->exitAdaptor2Decoded = &eVD[0] ;
      pp->exitAdaptor2Found = halloc((j+1) * sizeof(int), pp->h) ;
      pp->exitAdaptor2HasN = halloc((j+1) * sizeof(int), pp->h) ;
    }  

  if (pp->exitAdaptor2)
    for (ns = 0 ; pp->exitAdaptor2[ns] ; ns++)
      {
	unsigned char *ccp = pp->exitAdaptor2[ns] - 1 ;
	int nN = 0 ;
	while (*++ccp) if (*ccp == N_) nN++ ;
	pp->exitAdaptor2HasN[ns] = nN ;
      }

  return ;
} /* exitAdaptor2Init */

/*************************************************************************************/

static BOOL chenillette (const unsigned char *dna, int x, int width, int *dxLp, int *dxRp)
{
  const unsigned char *ccp ;
  unsigned char chenille[width], *cp ;
  int i, jL, jR ;

  memcpy (chenille, dna + x, width) ; 
  /* move left */
  for (ccp = dna + x, i = 0, jL = -1, cp = chenille ; *cp == *ccp && i + x >= 0 ; 
       ccp--, i--, cp = chenille + ((x*width+i) % width))
    jL++ ;
  /* move right */
  for (ccp = dna + x + width, i = jR = 0, cp = chenille ; *cp == *ccp ; 
       ccp++, i++, cp = chenille + (i % width))
    jR++ ;

  *dxLp = jL ; *dxRp = jR ;
  return jL + jR >= width ? TRUE : FALSE ;
} /* chenillette */

/*************************************************************************************/
/* if length > 16 we have to explicitely verify the hit */
BOOL clipAlignVerifyProbeHit (CLIPALIGN *pp, MM *mm, Array dna, Array dnaR, int pos, int xpos
			      , const unsigned char *probeDna, int probeLength
			      , BOOL isUp
			      , int *nerrp, int *nNp
			      , int *a1p, int *a2p
			      , int *x1p, int *x2p
			      , char *errLeftVal
			      , char *errRightVal
			      , unsigned char *prefix, unsigned char *suffix
			      , unsigned char *targetPrefix			      
			      )
{
  int i, j, ii, iLeft = 0, iRight = 0, lastShort, dxl, dxr ;
  int a1, a2, x1, x2, nS, x01, x02, score, nCorrected = 0 ;
  static Array err = 0, probeDnaArray = 0 ;
  const unsigned char *ccp ;
  unsigned char *cp ;
  int aliBonus = 0, leftAdaptorClip = 0, rightAdaptorClip = 0 ;
  A_ERR *ep = 0 ;
  float temperature = 0 ;
  int pairTampon = pp->pairTampon ;
  BOOL pAok = ! isUp || ! pp->strandedTarget ; /* if not TRUE, do not report polyA */
  BOOL pTok =   isUp || ! pp->strandedTarget ; /* if not TRUE, do not report polyT */

  nVerif++ ; 
  if (0) fprintf(stderr, "%d\t%s\t%d\t%s\n"
		 , nVerif, COSMID, pos
		 , stackText (pp->probeStack, mm->probeName)
		 ) ;
  if (!err) 
    {
      err = arrayCreate (256, A_ERR) ;
      probeDnaArray = arrayCreate (256, unsigned char) ;
    }
  /* transform the unsigned char buffer into a dna Array needed by aceDnaDoubleTrackErrors
   * notice that it is crucial in solid case to do this dynamically, 
   * since we freely mofify base 0 or probeDnaArray 
   */
  array (probeDnaArray, probeLength, unsigned char) = 0 ; /* make room */
  memcpy (arrp(probeDnaArray, 0, unsigned char), probeDna, probeLength) ;
  arrayMax (probeDnaArray) = probeLength ;

  if (pp->dx2)
    {
      x1 = xpos ;
      x2 = x1 + pp->seedLength - 1 ;
      a1 = pos + 1 ;
      a2 = a1 + pp->seedLength - 1 ;
    }
  else
    {
      x1 = xpos ; x2 = probeLength ;
      a1 = pos + 1 ; a2 = a1 + x2 - x1  ;
    }
  dxl = dxr = 9999 ;
  /* if the tag is inside a hit studied previously, drop it */
  if (pp->bestHit && bigArrayMax (pp->exportGeneHits) > 1 && mm->ipxg)
    {
      int b1, b2, ipxg ;
      PEXPORT *px ;
      int ipxgMax = bigArrayMax (pp->exportGeneHits) ;

      if (isUp)
	{
	  int max = arrayMax (dna) ;
	  b1 = max - a1 + 1 ; b2 = max - a2 + 1 ; 

	  for (ipxg = mm->ipxg ; ipxg> 0 && ipxg < ipxgMax  ; ipxg = px ? px->ipxold : 0)
	    {	     
	      px = bigArrp (pp->exportGeneHits, ipxg, PEXPORT) ;
	      if (pp->target == px->target && px->mm == mm)
		{
		  if (x1 >= px->x1 && x2 <= px->x2 && b1 <= px->a1 && b2 >= px->a2)
		    { nVerif-- ; return FALSE ; }
		}
	      else
		break ;
	    }
	}
      else
	{
	  b1 = a1, b2 = a2 ;

	  for (ipxg = mm->ipxg ; ipxg > 0 && ipxg < ipxgMax  ; ipxg = px ? px->ipxold : 0)  
	    { 
	      px = bigArrayp (pp->exportGeneHits, ipxg, PEXPORT) ;
	      if (pp->target == px->target && px->mm == mm)
		{
		  if (x1 >= px->x1 && x2 <= px->x2 && b1 >= px->a1 && b2 <= px->a2)
		    { nVerif-- ; return FALSE ; }
		} 
	      else
		    break ;
	    }
	}
    }
  

  /* Array aceDnaDoubleTrackErrors (Array  dna1, int *x1p, int *x2p, BOOL isDown,
                                      Array dna2, Array dna2R, int *a1p, int *a2p, 
                                      int *NNp, Array err, int maxJump, int maxError, BOOL doExtend
				      , int *maxExact)
  */
  aceDnaDoubleTrackErrors (probeDnaArray, &x1, &x2, TRUE /* bio coordinates  */
			   , dna, dnaR, &a1, &a2
			   , nNp, err, MAXJUMP, pp->errMax > 0 &&  pp->errMax < 999 ? 999 :  pp->errMax, TRUE, 0) ;
  
  /* clip out terminal error walls and compute the score */
  if (! arrayExists(err)) messcrash ("err does not exist") ;
  score = clipAlignErrorClip (pp, mm, err, probeDnaArray, dna, dnaR, nNp, &nCorrected, pos, xpos, &a1, &a2, &x1, &x2, isUp, &temperature) ;

  if (0) fprintf (stderr, "### score=%d  %d %d   %d %d\n", score, x1, x2, a1, a2) ;

  mm->lastTarget = pp->target ;
  mm->lastTargetGene = pp->targetGene ;
  mm->lastIsUp = isUp ;
  mm->lastA1 = a1 ; mm->lastA2 = a2 ; mm->lastX1 = x1 ; mm->lastX2 = x2 ;
  mm->lastDelta = pos - xpos ;
  if (! score) 
     return FALSE ;

  *a1p = a1 ; *a2p = a2 ; 
  *x1p = x1 ; *x2p = x2 ;
  x01 = x1 ; x02 = x2 ;

  strcpy (errLeftVal, "-") ; 
  strcpy (errRightVal, "-") ; 

  iLeft = iRight = -1 ;

  if (pp->X2X)
    { 
      char *cp, *cq ;
      int b1, b2, y1, y2, i, j, k = 0, ii, iiMax ; 
      int ln = arrayMax(dna) ;
      int lnX = arrayMax(probeDnaArray) ;
      Array err2 = arrayCopy (err) ;
      BOOL A2G = isUp ? pp->T2C : pp->A2G ;
      BOOL T2C = isUp ? pp->A2G : pp->T2C ;
      BOOL G2A = isUp ? pp->C2T : pp->G2A ;
      BOOL C2T = isUp ? pp->G2A : pp->C2T ;
      BOOL A2C = isUp ? pp->T2G : pp->A2C ;
      BOOL C2A = isUp ? pp->G2T : pp->C2A ;
      BOOL T2G = isUp ? pp->A2C : pp->T2G ;
      BOOL G2T = isUp ? pp->C2A : pp->G2T ;
      BOOL C2G = isUp ? pp->G2C : pp->C2G ;
      BOOL G2C = isUp ? pp->C2G : pp->G2C ;
      BOOL A2T = isUp ? pp->T2A : pp->A2T ;
      BOOL T2A = isUp ? pp->A2T : pp->T2A ;


      iiMax = arrayMax (err2) ;
      for (ii = 0 ; ii <= iiMax ; ii++)  
	{ /* for each error free segment */
	  int iiR = isUp ? iiMax - ii - 1 : ii ;
	  if (ii < iiMax)
	    {
	      A_ERR *eq = arrayp (err2, iiR, A_ERR) ;
	      y2 = isUp ? lnX - eq->iShort - 1 : eq->iShort ;
	      b2 = isUp ? ln - eq->iLong - 1 : eq->iLong ;
	      if (0)
		switch (eq->type)
		  {
		  case TROU:
		    if (isUp) b2 -= 0 ;
		    break ;
		  case TROU_DOUBLE:
		    if (isUp) b2 -= 0 ;
		    break ;
		  case TROU_TRIPLE:
		    if (isUp) b2 -= 0 ;
		    break ;
		  case INSERTION:
		    if (0)y2++ ;
		    break ;
		  case INSERTION_DOUBLE:
		    if (0) y2 += 2 ;
		    break ;
		  case INSERTION_TRIPLE:
		    if (0) y2 += 3 ;
		    break ;
		  default:
		    break ;
		  }
	    }
	  else
	    {
	      y2 = x2 ;
	      b2 = a2 ;
	    }
	  if (ii > 0)
	    {
	      A_ERR *er, *eq = arrayp (err2, iiR + (isUp ? 1 : - 1), A_ERR) ;
	      er = arrayp (err, k++, A_ERR) ;
	      *er = *eq ;
	      y1 = isUp ? lnX - eq->iShort - 1 : eq->iShort ;
	      b1 = isUp ? ln - eq->iLong - 1 : eq->iLong ;
	      switch (eq->type)
		{
		case TROU:
		  b1++ ;
		  if (isUp) { y1++ ; }
		  break ;
		case TROU_DOUBLE:
		  b1 += 2 ;
		  if (isUp) { b1-- ; y1++ ; }
		  break ;
		case TROU_TRIPLE:
		  b1 += 3 ;
		  if (isUp) { b1-=2 ; y1++ ; }
		  break ;
		case INSERTION:
		  y1++ ;
		  if (isUp) b1++ ; 
		  break ;
		case INSERTION_DOUBLE:
		  y1 += 2 ;
		  if (isUp) { b1++ ; y1-- ;} 
		  break ;
		case INSERTION_TRIPLE:
		  y1 += 3 ;
		  if (isUp) { b1++ ; y1-=2 ; }
		  break ;
		default:
		  y1++ ;
		  b1++ ;
		  break ;
		}
	    }
	  else
	    {
	      y1 = x1 - 1 ;
	      b1 = a1 - 1 ;
	    }
	  for (i = b1, j = y1, cp = arrp (dna, i, char), cq = arrp (probeDnaArray, j, char) ; i < b2 && j < y2 ; i++, j++, cp++, cq++) 
	    {
	      if (*cp != *cq)
		{
		  switch (*cp)
		    {
		    case R_:
		      if (A2G && *cq == A_)
			continue ;
		      if (G2A && *cq == G_)
			continue ;
		      break ;
		    case Y_:
		      if (T2C && *cq == T_)
			continue ;
		      if (C2T && *cq == C_)
			continue ;
		      break ;
		    case M_:
		      if (A2C && *cq == A_)
			continue ;
		      if (C2A && *cq == C_)
			continue ;
		      break ;
		    case K_:
		      if (T2G && *cq == T_)
			continue ;
		      if (G2T && *cq == G_)
			continue ;
		      break ;
		    case S_:
		      if (G2C && *cq == G_)
			continue ;
		      if (C2G && *cq == C_)
			continue ;
		      break ;
		    case W_:
		      if (A2T && *cq == A_)
			continue ;
		      if (T2A && *cq == T_)
			continue ;
		      break ;
		    }
		  {
		    A_ERR *eq = arrayp (err, k++, A_ERR) ;
		    eq->iLong = (isUp ? ln - i - 1 : i) ;
		    eq->iShort = (isUp ? lnX - j - 1 :  j) ;
		    eq->type = ERREUR ;
		    score-- ;
		  }
		}
	    }
	}
      arrayDestroy (err2) ;
    }

  *nerrp = arrayMax (err) ; 
  if (*nerrp == 0)
     nS = a2 - a1 + 1 ; 
  else
    {
      BOOL A2G = pp->A2G ;
      BOOL T2C = pp->T2C ;
      BOOL G2A = pp->G2A ;
      BOOL C2T = pp->C2T ;
      BOOL A2C = pp->A2C ;
      BOOL C2A = pp->C2A ;
      BOOL T2G = pp->T2G ;
      BOOL G2T = pp->G2T ;
      BOOL C2G = pp->C2G ;
      BOOL G2C = pp->G2C ;
      BOOL A2T = pp->A2T ;
      BOOL T2A = pp->T2A ;
      BOOL X2X = pp->X2X ;

      /* count the  secondary errors */
      for (ii = 0, ep = arrp (err, 0, A_ERR) ; ii < arrayMax (err) ; ii++, ep++)
	{
	  if (ep->type == TYPE80 || ep->type == AMBIGUE)
	    (*nerrp)-- ; 
	}
      /* locate the center of the best segment */
      nS = 0 ; j = (x01 + x02)/2 ; 
      for (lastShort = x1 - 1, ii = 0, ep = arrp (err, 0, A_ERR) ; ii < arrayMax (err) ; ii++, ep++)
	{
	  if (ep->type == TYPE80 || ! ep->type) continue ;
	  if (ep->iShort < x1)
	    continue ;
	  if (ep->iShort > x2)
	    continue ;
	  if (nS < ep->iShort - lastShort)
	    { 
	      iLeft =  ii - 1 ; iRight = ii ;
	      nS = ep->iShort - lastShort - 1 ; 
	    }
	  lastShort =  ep->iShort ;
	}
      if (nS < x02 - lastShort)
	{ iLeft =  ii - 1 ; iRight = -1 ; nS = x02 - lastShort ; }

      if (arrayMax (err))  /* in errLeftVal, report all the errors in the aligned segment */
	{
	  int xShort, xLong  ;
	  char *cp = errLeftVal ;
	  char *cr = errRightVal ;
	  
	  dxl = 0 ;
	  for (ii = 0 ; ii < arrayMax (err) ; ii++)
	    {
	      ep = arrp (err, ii, A_ERR) ;
	      if (ii) 
		{
		  cp += strlen(cp) ;
		  if (cp - errLeftVal > 10000)
		    break ;
		  *cp++ = ',' ; 
		  cr += strlen(cr) ;
		  if (cr - errRightVal > 10000)
		    break ;
		  *cr++ = ',' ; 
		}
	      xShort = isUp ? arrayMax (probeDnaArray) - ep->iShort : ep->iShort + 1 ;
	      xLong = ep->iLong + 1 ;
#ifdef JUNK
	      { 
		char *cqual, cshort ;
		cshort = arr(probeDnaArray, xShort - 1, unsigned char) ;
		if (mm->probeQuality)
		  { 
		    int q =  *(stackText (pp->qualityStack, mm->probeQuality + xShort - 1)) - 64 ;
		    if (q < 0) q = 0 ;
		    cqual = messprintf ("#%d", q) ;
		  }
		else
		  cqual = "" ;
	      }
#endif
	      switch (ep->type)
		{
		case TYPE80:
		  {
		    char cc1a, cc2a ;
		    int xLongR ;
		    
		    xLongR = arrayMax (pp->dna0) - xLong + 1 ;

		    cc1a = arr (pp->dna0, (isUp ? xLongR - 0 : xLong - 2), unsigned char) ;
		    cc2a = arr (pp->dna0, (isUp ? xLongR - 1 : xLong - 1), unsigned char) ;

		    sprintf(cp,"%d:%c%c>oo"
			    , isUp ? xLong - 1 : xLong - 1
			    , isUp ? dnaDecodeChar[(int)complementBase(cc1a)] : dnaDecodeChar[(int)cc1a]
			    , isUp ? dnaDecodeChar[(int)complementBase(cc2a)] : dnaDecodeChar[(int)cc2a]
			    ) ;
		    
		    sprintf(cr,"%d:%c%c>oo"
			    , xShort - 1
			    , isUp ? dnaDecodeChar[(int)cc2a] : dnaDecodeChar[(int)cc1a]
			    , isUp ? dnaDecodeChar[(int)cc1a] : dnaDecodeChar[(int)cc2a]
	 		    ) ;
		  }
		  break ;
		case AMBIGUE:
		case ERREUR:
		  {
		    char cc1o, cc2o, cc1oc, cc2oc, cc1a, cc2a, cc1ac, cc2ac ;
		    int xLongR ;
		    
		    xLongR = arrayMax (pp->dna0) - xLong + 1 ;

		    if (pp->solid)
		      {
	 		cc1a = clipAlignSolidProduct (arr(pp->dna0, (isUp ? xLongR : xLong - 2),unsigned char)
						          , probeDna[xShort - 1]) ;
		      }
		    else
		      {
	 		cc1a = probeDna[xShort - 1] ;
		      }
		    cc1o = arr(pp->dna0, isUp ? xLongR - 1 : xLong - 1,unsigned char) ;

		    if (X2X)
		      {
			if (A2G)
			  {
			    if (cc1o == R_) cc1o = A_ ;
			    if (cc1o == Y_) cc1o = T_ ;
			  }
			if (T2C)
			  {
			    if (cc1o == R_) cc1o = A_ ;
			    if (cc1o == Y_) cc1o = T_ ;
			  }
			if (G2A)
			  {
			    if (cc1o == R_) cc1o = G_ ;
			    if (cc1o == Y_) cc1o = C_ ;
			  }
			if (T2C)
			  {
			    if (cc1o == Y_) cc1o = T_ ;
			    if (cc1o == R_) cc1o = A_ ;
			  }
			if (C2T)
			  {
			    if (cc1o == Y_) cc1o = C_ ;
			    if (cc1o == R_) cc1o = G_ ;
			  }
			if (A2C)
			  {
			    if (cc1o == M_) cc1o = A_ ;
			    if (cc1o == K_) cc1o = T_ ;
			  }
			if (C2A)
			  {
			    if (cc1o == M_) cc1o = C_ ;
			    if (cc1o == K_) cc1o = G_ ;
			  }
			if (T2G)
			  {
			    if (cc1o == K_) cc1o = T_ ;
			    if (cc1o == M_) cc1o = A_ ;
			  }
			if (G2T)
			  {
			    if (cc1o == K_) cc1o = G_ ;
			    if (cc1o == M_) cc1o = C_ ;
			  }
			if (G2C)
			  {
			    if (cc1o == S_) cc1o = isUp ? C_ : G_ ;
			  }
			if (C2G)
			  {
			    if (cc1o == S_) cc1o = isUp ? G_ : C_ ;
			  }
			if (A2T)
			  {
			    if (cc1o == W_) cc1o = isUp ? T_ : A_ ;
			  }
			if (T2A)
			  {
			    if (cc1o == W_) cc1o = isUp ? A_ : T_ ;
			  }
		      }

		    cc2o = isUp ? complementBase(cc1o) : cc1o ; 
		    cc2a = isUp ? complementBase(cc1a) : cc1a ; 
		    cc1oc = dnaDecodeChar[(int)cc1o] ;
		    cc2oc = dnaDecodeChar[(int)cc2o] ;
		    cc1ac = dnaDecodeChar[(int)cc1a] ;
		    cc2ac = dnaDecodeChar[(int)cc2a] ;

		    sprintf(cp,"%d:%c>%c"
			    , xLong 
			    , cc2oc, cc2ac
			    ) ;
		    
		    sprintf(cr,"%d:%c>%c"
			    , xShort + (isUp ? mm->probeLeftClip + (pp->solid ? -1 : 0) : mm->probeLeftClip )
			    , cc1oc, cc1ac
			    ) ;
		  }

		  break ;
		case TROU: 
		  {
		    char *ss = "-", cc1a, cc2a, cc1ac, cc2ac ;
		    int dxL = 0, dxR = 0, xLongR ;
		    
		    if (pp->solid) xLong++ ;
		    xLongR = arrayMax (pp->dna0) - xLong + 1 + (pp->solid ? 0 : 0) ;
		    if (ii == 0)
		      {
			lastShort-- ;
		      }
		    
		    if (1)  /* solid is not an issue for deletions */
		      {
			if (!pp->solid)
			  cc1a = arr(pp->dna0, (isUp ? xLongR - 1 + 0 : xLong - 1 + 0), unsigned char) ;
			else
			  cc1a = arr(pp->dna0, (isUp ? xLongR - 1 - 1 : xLong - 1 + 1), unsigned char) ;
		      }
		    
		    if (pp->X2X)
		      {
			if (cc1a == R_ || cc1a == M_) cc1a = A_ ;
			if (cc1a == Y_ || cc1a == K_) cc1a = T_ ;
		      }
		    cc2a = isUp ? complementBase(cc1a) : cc1a ; 
		    
		    cc1ac = dnaDecodeChar[(int)cc1a] ;
		    cc2ac = dnaDecodeChar[(int)cc2a] ;

		    if (chenillette (arrp(pp->dna0, 0, unsigned char), isUp ? xLongR - 1 - (pp->solid ? 1 : 0) : xLong - 1 + (pp->solid ? 1 : 0), 1, &dxL, &dxR))
		      {
			ss = "*-" ;
			if (0 && pp->solid && isUp) xLong++ ;
			if (0 && pp->solid && isUp) xShort-- ;
		      }
		    
		    sprintf(cp,"%d:%s%c"
			    , xLong + (pp->solid ? 1 : 0)
			    , ss 
			    , cc2ac
			    ) ;
		    
		    sprintf(cr,"%d:%s%c"
			    , xShort + (isUp ? mm->probeLeftClip + (pp->solid ? 0 : 1) : mm->probeLeftClip )
			    , ss
			    , cc1ac
			    ) ;
		  }
		  break ;
		case TROU_DOUBLE:
		  {
		    char *ss = "--", cc1a, cc1b, cc2a, cc2b, cc1ac, cc1bc, cc2ac, cc2bc ;
		    int dxL = 0, dxR = 0, xLongR ;
		    
		    xLongR = arrayMax (pp->dna0) - xLong + 1 ;
		    if (ii == 0)
		      {
			lastShort-- ;
			lastShort-- ;
		      }

		    if (chenillette (arrp(pp->dna0, 0, unsigned char), isUp ? xLongR - 2 : xLong - 1, 2, &dxL, &dxR))
		      {
			ss = "*--" ; 
			xShort += (pp->solid ? 0 : 0) ;
		      }

		    if (1)  /* solid is not an issue for deletions */
		      {
			cc1a = arr(pp->dna0, (isUp ? xLongR - 2 + 0 : xLong -1 + 0), unsigned char) ;
			cc1b = arr(pp->dna0, (isUp ? xLongR - 2 + 1 : xLong -1 + 1), unsigned char) ;
		      }

		    if (pp->X2X)
		      {
			if (cc1a == R_ || cc1a == M_) cc1a = A_ ;
			if (cc1a == Y_ || cc1a == K_) cc1a = T_ ;
			if (cc1b == R_ || cc1b == M_) cc1b = A_ ;
			if (cc1b == Y_ || cc1b == K_) cc1b = T_ ;
		      }

		    cc2a = isUp ? complementBase(cc1b) : cc1a ; 
		    cc2b = isUp ? complementBase(cc1a) : cc1b ; 

		    cc1ac = dnaDecodeChar[(int)cc1a] ;
		    cc1bc = dnaDecodeChar[(int)cc1b] ;
		    cc2ac = dnaDecodeChar[(int)cc2a] ;
		    cc2bc = dnaDecodeChar[(int)cc2b] ;

		    sprintf(cp,"%d:%s%c%c"
			    , xLong 
			    , ss 
			    , cc2ac, cc2bc
			    ) ;
		    
		    sprintf(cr,"%d:%s%c%c"
			    , xShort + (isUp ? mm->probeLeftClip + (pp->solid ? 0 : 1)  : mm->probeLeftClip )
			    , ss
			    , cc1ac, cc1bc
			    ) ;
		  }
		  break ;
		case TROU_TRIPLE:
		  {
		    char *ss = "---", cc1a, cc1b, cc1c, cc2a, cc2b, cc2c, cc1ac, cc1bc, cc1cc, cc2ac, cc2bc, cc2cc ;
		    int dxL = 0, dxR = 0, xLongR ;
		    
		    xLongR = arrayMax (pp->dna0) - xLong + 1 ;
		    if (ii == 0)
		      {
			lastShort-- ;
			lastShort-- ;
			lastShort-- ;
		      }

		    if (chenillette (arrp(pp->dna0, 0, unsigned char), isUp ? xLongR - 3 : xLong - 1, 3, &dxL, &dxR))
		      {
			ss = "*---" ;
			xShort += (pp->solid ? 0 : 0 ) ;
		      }

		    if (1)  /* solid is not an issue for deletions */
		      {
			cc1a = arr(pp->dna0, (isUp ? xLongR - 3 + 0 : xLong -1 + 0), unsigned char) ;
			cc1b = arr(pp->dna0, (isUp ? xLongR - 3 + 1 : xLong -1 + 1), unsigned char) ;
			cc1c = arr(pp->dna0, (isUp ? xLongR - 3 + 2 : xLong -1 + 2), unsigned char) ;
		      }

		    if (pp->X2X)
		      {
			if (cc1a == R_ || cc1a == M_) cc1a = A_ ;
			if (cc1a == Y_ || cc1a == K_) cc1a = T_ ;
			if (cc1b == R_ || cc1b == M_) cc1b = A_ ;
			if (cc1b == Y_ || cc1b == K_) cc1b = T_ ;
			if (cc1c == R_ || cc1c == M_) cc1c = A_ ;
			if (cc1c == Y_ || cc1c == K_) cc1c = T_ ;
		      }

		    cc2a = isUp ? complementBase(cc1c) : cc1a ; 
		    cc2b = isUp ? complementBase(cc1b) : cc1b ; 
		    cc2c = isUp ? complementBase(cc1a) : cc1c ; 

		    cc1ac = dnaDecodeChar[(int)cc1a] ;
		    cc1bc = dnaDecodeChar[(int)cc1b] ;
		    cc1cc = dnaDecodeChar[(int)cc1c] ;
		    cc2ac = dnaDecodeChar[(int)cc2a] ;
		    cc2bc = dnaDecodeChar[(int)cc2b] ;
		    cc2cc = dnaDecodeChar[(int)cc2c] ;

		    sprintf(cp,"%d:%s%c%c%c"
			    , xLong 
			    , ss 
			    , cc2ac, cc2bc, cc2cc
			    ) ;
		    
		    sprintf(cr,"%d:%s%c%c%c"
			    , xShort + (isUp ? mm->probeLeftClip  + (pp->solid ? 0 : 1) : mm->probeLeftClip )
			    , ss
			    , cc1ac, cc1bc, cc1cc
			    ) ;
		  }
		  break ;
		case INSERTION: 
		  {
		    char *ss = "+", cc1a, cc2a, cc1ac, cc2ac ;
		    int dxL = 0, dxR = 0, xLongR ;
		    
		    xLongR = arrayMax (pp->dna0) - xLong + 1 ;

		    if ((ii == 0 || ep->iShort > (ep-1)->iShort + 1) && chenillette (probeDna, xShort - 1, 1, &dxL, &dxR))
		      ss = "*+" ;
		    if (pp->solid)
		      {
			cc1a = clipAlignSolidProduct (arr(pp->dna0, (isUp ? xLongR : xLong - 2),unsigned char)
						          , probeDna[xShort - 1 + (isUp ? -0 : 0)]) ;
		      }
		    else
		      {
			cc1a = probeDna[xShort - 1 + (isUp ? -0 : 0)] ;
		      }

		    if (pp->X2X)
		      {
			if (cc1a == R_ || cc1a == M_) cc1a = A_ ;
			if (cc1a == Y_ || cc1a == K_) cc1a = T_ ;
		      }


		    cc2a = isUp ? complementBase(cc1a) : cc1a ; 
		    cc1ac = dnaDecodeChar[(int)cc1a] ;
		    cc2ac = dnaDecodeChar[(int)cc2a] ;

		    sprintf(cp,"%d:%s%c"
			    , xLong 
			    , ss 
			    , cc2ac
			    ) ;
		    
		    sprintf(cr,"%d:%s%c"
			    , xShort + (isUp ? mm->probeLeftClip - 0 : mm->probeLeftClip )
			    , ss
			    , cc1ac
			    ) ;
		  }
		  break ;
		case INSERTION_DOUBLE:
		  {
		    char *ss = "++", cc1a, cc1b, cc2a, cc2b, cc1ac, cc1bc, cc2ac, cc2bc ;
		    int dxL = 0, dxR = 0, xLongR ;
		    
		    xLongR = arrayMax (pp->dna0) - xLong + 1 ;

		    if ((ii == 0 || ep->iShort > (ep-1)->iShort + 1) && chenillette (probeDna, xShort - 1, 2, &dxL, &dxR))
		      ss = "*++" ;
		    if (pp->solid)
		      {
			cc1a = clipAlignSolidProduct (arr(pp->dna0, (isUp ? xLongR : xLong - 2),unsigned char)
						          , probeDna[xShort - 1 + (isUp ? -0 : 0)]) ;
			cc1b = clipAlignSolidProduct (cc1a, probeDna[xShort - 1 + (isUp ? -1 : 1)]) ;
			if (isUp)
			  { int cc = cc1a ; cc1a = cc1b ; cc1b = cc ; }
		      }
		    else
		      {
			cc1a = probeDna[xShort - 1 + (isUp ? -1 : 0)] ;
			cc1b = probeDna[xShort - 1 + (isUp ? -0 : 1)] ;
		      }


		    if (pp->X2X)
		      {
			if (cc1a == R_ || cc1a == M_) cc1a = A_ ;
			if (cc1a == Y_ || cc1a == K_) cc1a = T_ ;
			if (cc1b == R_ || cc1b == M_) cc1b = A_ ;
			if (cc1b == Y_ || cc1b == K_) cc1b = T_ ;
		      }


		    cc2a = isUp ? complementBase(cc1b) : cc1a ; 
		    cc2b = isUp ? complementBase(cc1a) : cc1b ; 
		    cc1ac = dnaDecodeChar[(int)cc1a] ;
		    cc1bc = dnaDecodeChar[(int)cc1b] ;
		    cc2ac = dnaDecodeChar[(int)cc2a] ;
		    cc2bc = dnaDecodeChar[(int)cc2b] ;

		    sprintf(cp,"%d:%s%c%c"
			    , xLong 
			    , ss 
			    , cc2ac, cc2bc
			    ) ;
		    
		    sprintf(cr,"%d:%s%c%c"
			    , xShort + (isUp ? mm->probeLeftClip - 1 : mm->probeLeftClip )
			    , ss
			    , cc1ac, cc1bc
			    ) ;
		  }
		  break ;
		case INSERTION_TRIPLE:
		  {
		    char *ss = "+++", cc1a, cc1b, cc1c, cc2a, cc2b, cc2c, cc1ac, cc1bc, cc1cc, cc2ac, cc2bc, cc2cc ;
		    int dxL = 0, dxR = 0, xLongR ;
		    
		    xLongR = arrayMax (pp->dna0) - xLong + 1 ;
		    
		    if ((ii == 0 || ep->iShort > (ep-1)->iShort + 1) && chenillette (probeDna, xShort - 1, 3, &dxL, &dxR))
		      ss = "*+++" ;
		    if (pp->solid)
		      {
			cc1a = clipAlignSolidProduct (arr(pp->dna0, (isUp ? xLongR : xLong - 2),unsigned char)
						          , probeDna[xShort - 1 + (isUp ? -0 : 0)]) ;
			cc1b = clipAlignSolidProduct (cc1a, probeDna[xShort - 1 + (isUp ? -1 : 1)]) ;
			cc1c = clipAlignSolidProduct (cc1b, probeDna[xShort - 1 + (isUp ? -2 : 2)]) ;
			if (isUp)
			  { int cc = cc1a ; cc1a = cc1c ; cc1c = cc ; }
		      }
		    else
		      {
			cc1a = probeDna[xShort - 1 + (isUp ? -2 : 0)] ;
			cc1b = probeDna[xShort - 1 + (isUp ? -1 : 1)] ;
			cc1c = probeDna[xShort - 1 + (isUp ? -0 : 2)] ;
		      }


		    if (pp->X2X)
		      {
			if (cc1a == R_ || cc1a == M_) cc1a = A_ ;
			if (cc1a == Y_ || cc1a == K_) cc1a = T_ ;
			if (cc1b == R_ || cc1b == M_) cc1b = A_ ;
			if (cc1b == Y_ || cc1b == K_) cc1b = T_ ;
			if (cc1c == R_ || cc1c == M_) cc1c = A_ ;
			if (cc1c == Y_ || cc1c == K_) cc1c = T_ ;
		      }


		    cc2a = isUp ? complementBase(cc1c) : cc1a ; 
		    cc2b = isUp ? complementBase(cc1b) : cc1b ; 
		    cc2c = isUp ? complementBase(cc1a) : cc1c ; 
		    cc1ac = dnaDecodeChar[(int)cc1a] ;
		    cc1bc = dnaDecodeChar[(int)cc1b] ;
		    cc1cc = dnaDecodeChar[(int)cc1c] ;
		    cc2ac = dnaDecodeChar[(int)cc2a] ;
		    cc2bc = dnaDecodeChar[(int)cc2b] ;
		    cc2cc = dnaDecodeChar[(int)cc2c] ;

		    sprintf(cp,"%d:%s%c%c%c"
			    , xLong 
			    , ss 
			    , cc2ac, cc2bc, cc2cc
			    ) ;
		    
		    sprintf(cr,"%d:%s%c%c%c"
			    , xShort + (isUp ? mm->probeLeftClip - 2 : mm->probeLeftClip )
			    , ss
			    , cc1ac, cc1bc, cc1cc
			    ) ;

		  }
		  break ;
		}
	    }
	}
    }
  /* N in the reference sequence around 'pos' count as  genuine errors */
  ccp = arrp (dna, pos, unsigned char) - 1 ; i = *a2p - pos + 1 ;
  while (i-- > 0  && *++ccp == N_) { (*nNp)-- ; (*nerrp)++ ; score -= pp->errCost ; }
  ccp = arrp (dna, pos, unsigned char) - 2 ; i = pos - *a1p + 0 ;
  while (i-- > 0 && ccp >= (unsigned char*)dna->base && *--ccp == N_) { (*nNp)-- ; (*nerrp)++ ; }
  if (*nNp<0) *nNp = 0 ;
  if (10 * (*nNp) > nS && nS < 1000 && 2*nS < arrayMax (dna) ) /* reject over 10% ambiguous bases */
    return FALSE ;
  *nNp += nCorrected ;

  /* in case of success, report the prefix/suffix */
  memset (prefix, 0, 64+OVLN) ;
  memset (suffix, 0, 64+OVLN) ;
  memset (targetPrefix, 0, 64 + 4*OVLN) ;
  prefix[0] = 0 ; prefix[1] = '-' ; prefix[2] = 0 ;
  suffix[0] = 0 ; suffix[1] = '-' ; suffix[2] = 0 ;
  targetPrefix[0] = '-' ; targetPrefix[1] = 0 ;

  /**************************/
  /* Export the prefix      */
  /**************************/

  if ((pp->splice || pp->clipPolyT || pp->showOverhang || pp->sl || pp->entryAdaptor) && x1 > 1)
    {
      int goodPrefix = 0, genomeLn ;
      const unsigned char *genome ;
      unsigned char buf[x1+OVLN+1] ;
      int n = 0, pN = 0 ;

      memset (buf, 0, sizeof(buf)) ;
      /* construct the prefix, in the orientation opposite to  the probe */
      if (!pp->solid)
	{
	  for (cp = buf, i = x1 - 2, ccp = probeDna + i ;
	       i > 0 ; i--, cp++, ccp--)
	    {
	      *cp = complementBase(*ccp) ;
	    }
	  *cp = 0 ;
	}
      else /* decipher the prefix */
	{
	  unsigned char bb ;
	  
	  /* we start on the last correct  base */
	  ccp = arrp(pp->dna0, a1 - 1, unsigned char) ;
	  bb = *ccp ;
	  /* export the unaligned probe bases */
	  /* start at prefix+0 because on acceptor side the 
	   * first error is an erroneous transition, but a correct base 
	   * so this base cannot be found upstream in the genome 
	   */
	  j = OVLN ; if (j > x1 - 2) j = x1 - 1 ; buf[j] = 0 ; 
	  for (cp = buf + j - 1, i = x1, ccp = probeDna + i - 1 ;
	       i > 1 && j > 0 ; j--, i--, cp--, ccp--)
	    {
	      bb = clipAlignSolidProduct (*ccp, bb) ; /* true value of the erroneous letter */
	      *cp = bb ; /* do not dnaDecode */
	    }
	  cp = buf + strlen((char*)buf) - 1 ; /* kill the last letter */
	}
      *cp = 0 ;
      strncpy ((char *)prefix+1, (char*)buf, OVLN) ;
      
      for (cp = buf + strlen((char*)prefix+1) - 1, pN = 0, i = 0 ; cp >= buf && *cp && i < 8 ; i++, cp--)
	if (*cp == N_) pN++ ;

      if (0 && !pN && pp->clipPolyT) /* look for 9 pure T */
	{
	  for (n = 0, i = strlen ((char*)prefix), ccp = prefix ; 
	       i > 0 && *ccp == T_ ; i++, ccp++)
	    n++ ;
	  if (n >= 9)
	    goodPrefix = 2 ;
	}

      /* locate the upstream genome */
      i = a1 - 1 ;
      i = i - pp->intronMaxLength ; if (i < 0) i = 0 ;
      genome = arrp (pp->dna0, i, unsigned char) ;
      genomeLn = a1 - i - 12 ;

      /* search for transpliced leaders in the prefix of the probe */
      if (! pN && pp->sl && !goodPrefix && 
	  pAok &&
	  strlen((char*)buf) < 24 && 
	  strlen((char*)buf) > 8) /* search for SL1 */
	{
	  int ns, n1 ;
	  unsigned char buf2[OVLN+3] ;
	  
	  for (ns = 0 ; pp->sl[ns] && !goodPrefix ; ns++)
	    {
	      ccp = pp->sl[ns] ;
	      i =  strlen((char *)ccp)  ;
	      n1 = dnaPickMatch (ccp, i, buf, 0, 0) ;
	      if (n1 && ! pp->solid &&  mm->probeLeftClip)
		{  /* left TTTT should also be present in the SL motif */
		  for (i = 0 ; i < mm->probeLeftClip ; i++)
		    if (n1 < 2+i || ccp[n1-2-i] != T_) 
		      n1 = 0 ;
		}
	      if (n1) /* the prefix matches in the SL, */
		{
		  /* locate the tail of the SL present in the tag */
		  ccp = pp->sl[ns] + n1 - 1 ; 
		  /* verify that the ag acceptor consensus was in the aligned part */
		  i = strlen((char *)ccp) ;
		  if (i + 1 < x1 + (pp->solid ? 0 : 0)) continue ;
		  /* verify that the whole tail of the SL matches the ttag */
		  if (! pp->solid && ! dnaPickMatch (probeDna, mm->probeLength, ccp, 0, 0))
		    continue ;
		  /* check that the prefix including an extra donor 'gt' 
		   * is not present in the genome
		   */
		  int kkk = ccp ? strlen ((const char*)ccp) : 0 ;
                  BOOL inGenome = FALSE ;
		  if (kkk < OVLN)
		    {
		      sprintf ((char *)buf2, "%s%c%c", ccp, G_, T_ | C_) ;
		      inGenome = dnaPickMatch (genome, genomeLn, buf2, 0, 0) ;
		    }
		  if (! inGenome)
		    {
		      goodPrefix = ns + 11 ;
		      score += pp->slBonus ;
		      i = 1 + strlen ((char *)ccp) - x1 ;
		      if (i > 0)
			{
			  a1 = (*a1p) += i ;
			  x1 = (*x1p) += i ;
			  strncpy ((char *)buf, (char *)ccp, OVLN) ;
			}
		      aliBonus += 0 ; /*  strlen((char *)buf) + mm->probeLeftClip ; */
		      leftAdaptorClip = -x1 + 1 ; /* substract from the to be aligned length */
		    }
		}
	    }
	}
      /* search for entryAdaptor (possibly a short bar-code) */
      int bufLn = strlen((char*)prefix+1) ;
      unsigned char **entryAdaptor  = (mm->pair < 0 ? pp->entryAdaptor2 : pp->entryAdaptor) ;
      unsigned  int *entryAdaptorHasN = (mm->pair <0 ? pp->entryAdaptor2HasN : pp->entryAdaptorHasN ) ;
      if (! pN && entryAdaptor && !goodPrefix && 
	  bufLn < 200 && 
	  bufLn >= 5) /* search for the bar-code or entry adaptor */
	{
	  int ns, n1, petitMalus = 0 ;
	  unsigned char buf2[OVLN+3] ;
	  
	  for (ns = 0 ; entryAdaptor[ns] && !goodPrefix ; ns++)
	    {
	      int aAT = 0, sAT = 0 ;
	      int adapLn = 0 ;
	      petitMalus = 0 ;
	      ccp = entryAdaptor[ns] ;
	      adapLn = strlen((char *)ccp)  ;
	      if (pp->clipPolyT)
		{
		  for (aAT = 0 ; ccp[aAT] == A_  ; aAT++)
		    ;
		  for (sAT = 0 ; prefix[1+sAT] == A_ ; sAT++)
		    ;
		}
	      if (aAT > sAT)
		continue ;
	      i =  adapLn - aAT ;
	      n1 = dnaPickMatch (ccp + aAT, i, prefix + 1 + sAT, i > 5 ? 1 : 0 , 0) ;
	      if (n1 && ! pp->solid &&  mm->probeLeftClip)
		{  /* left TTTT should also be present in the SL motif */
		  for (i = 0 ; i < mm->probeLeftClip  - sAT ;  i++)
		    if (n1 < 2+i || ccp[n1-2-i] != T_) 
		      n1 = 0 ;
		}
	      if (n1 <= 0)
		continue ;
	      if (adapLn - n1 < 6 && bufLn > 30)
		continue ;
	      if (adapLn - n1 < 10 && bufLn > 100)
		continue ;
	      if (adapLn - n1 < 20 && bufLn > 1000)
		continue ;
	      if (n1 >= 0) /* the prefix matches in the adaptor */
		{
		  /* locate the tail of the adaptor present in the tag */
		  ccp = entryAdaptor[ns] + n1 - 1 + aAT ; 
		  i = strlen((char *)ccp) ;
		  petitMalus = i + 1 - x1 ;
		  /* verify that the whole tail of the adaptor matches the ttag */
		  if (! pp->solid && n1 > 1)
		    {
		      for (i = 0 ; i < n1 - 1 ; i++)
			if (probeDna[ x1 - 2 - i] != complementBase( ccp[ - i-1]))
			  n1 = 0 ;
		    }
		  /* check that the prefix 
		   * is not present in the genome
		   */
                  int kkk = ccp ? strlen ((const char *)ccp) : 0 ;
                  BOOL inGenome = FALSE ;
		  if (kkk && kkk < OVLN)
		    {
		      if (0) fprintf (stderr , "zzzzz   kkk=%d ccp=%s\n", kkk, ccp) ;
		      sprintf ((char *)buf2, "%s%c%c", ccp, G_, T_ | C_) ;
		      inGenome = dnaPickMatch (genome, genomeLn, buf2, 0, entryAdaptorHasN[ns]) ;
		    }
		  if (! inGenome)
		    {
		      goodPrefix = 3 ;  /* report Adaptor */
		      if (goodPrefix > 250) goodPrefix = 250 ;
		      score -= petitMalus ; /* do not count the aligned letters which belong to the adaptor */
		      i = 1 + strlen ((char *)ccp) - x1 ;
		      if (i > 0)
			{
			  a1 = (*a1p) += i ;
			  x1 = (*x1p) += i ;
			}
		      leftAdaptorClip = -x1 + 1 ; /* substract from the to be aligned length */
		    }
		}
	    }
	}
      /* search for a polyT */
      if (! pN && !goodPrefix && 
	  pp->clipPolyT &&
	  bufLn < 200 && 
	  bufLn >= 8) /* search for the bar-code or entry adaptor */
	{
	  int sAT = 0 ;
	  for (sAT = 0 ; prefix[1+sAT] == A_ ; sAT++)
		    ;
	  if (sAT >= 8)
	    goodPrefix = 3 ;
	}

      /* search for an upstream donor exon */
      if (! pN && (pp->splice || pp->showOverhang || pp->clipPolyT2) && !goodPrefix && pp->overhangLength > 3 && strlen ((char *)buf) >= pp->overhangLength)
	{
	  i = pp->overhangLength < x1 - 1 ? pp->overhangLength : x1 - 1 ;
	  goodPrefix = 1 ;  /* so it is available to construct introns */
	  n = 1 ; /* dnaPickMatch (genome, genomeLn, buf + strlen((char*)buf) - i, 0, 0) ; */
	  if (pp->splice) /* 2013_08_13, was if(1) */
	    {
	      goodPrefix = 63 ; 
	      score += INTRONBONUS ; /* needed to find best score will be removed at export time 
				      * was: score += strlen((char*)buf) - i 
				      */
	    }
	}
      
      prefix[0] = goodPrefix ; /* flag it is an SL a poly-A or an intron candidate */
      
      /* fix the prefix */
      for (cp = prefix+1 ; *cp ; cp++)
	{
	  *cp = dnaDecodeChar[(int)*cp] ;
	  if (goodPrefix > 1 && goodPrefix != 63)
	    *cp = ace_upper (*cp) ;
	}
    }
	
  /**************************/
  /* Export the suffix      */
  /**************************/
  
  if ((pp->splice || pp->clipPolyA || pp->clipPolyT || pp->exitAdaptor  || pp->exitAdaptor2 || pp->showOverhang || (pp->sl && !( mm->stranded > 0))) && x2 < mm->probeLength)
    { 
      int pN = 0, goodSuffix = 0, genomeLn, adaptorSuffix = 0 ;
      unsigned char *genome ;
      unsigned char cc, bb = arr (pp->dna0, a2 - 1, unsigned char) ;
      
      i = x2 + 1 ;
      cp = suffix + 1 ; ccp = probeDna + i - 1 ;

      /* construct the suffix */
      for (j = 0 ; *ccp && j < OVLN ; cp++, ccp++, i++, j++)
	{
	  if (pp->solid)
	    {
	      bb = clipAlignSolidProduct (*ccp, bb) ; /* true value of the probe */
	      *cp = bb ; /* do not decodeUnsigned Char */
	    }
	  else
	    *cp = *ccp ;
	}
      *cp = 0 ;
      
     for (cp = suffix + 1, pN = 0, i = 0 ; *cp && i < 8 ; i++, cp++)
       if (*cp == N_) pN++ ;
     
     /* locate the downstream genome */
     {
       int a22 = a2 + 12 ;
       if (a22 > arrayMax (pp->dna0) - 1)
	 a22 = arrayMax (pp->dna0) - 1 ;	 
       
       genome = arrayp (pp->dna0, a22, unsigned char) ;
       genomeLn = arrayMax(dna) - a22 ;
     }
     
     if (genomeLn > pp->intronMaxLength)
       genomeLn = pp->intronMaxLength ;
     
     /* search for transpliced leaders in the suffix of the probe */
     if (! pN && pp->slR && !goodSuffix && !( mm->stranded > 0) &&
	 pTok &&
	 strlen((char *)suffix + 1) < 24 && 
	 strlen((char *)suffix + 1) > 8) /* search for SL1 */
       {
	 int ns, n1 ;
	 unsigned char buf2[OVLN+3] ;
	 
	 for (ns = 0 ; pp->slR[ns] && !goodSuffix ; ns++)
	   {
	     ccp = pp->slR[ns] ;
	     i =  strlen((char *)ccp)  ;  
	     n1 = dnaPickMatch (ccp, i, suffix+1, 0, 0) ;
	     /* the SL must also match the trailing A */
	     if (n1 > 0)
	       {
		 for (i = 0 ; n1 && i < mm->probeRightClip ; i++)
		   if (ccp [n1 + strlen((char *)suffix+1) - 1 + i] != A_)
		     n1 = 0 ;
	       }
	     if (n1 > 0) /* the suffix matches in the SL, 
			   * and the ag acceptor consensus was in the aligned part 
			   */
		{
		  /* verify that the bp of the SL upstream of n1 are in the exon */
		  for (i = 1, j = 1 ; j && i < n1 ; i++)
		    if (ccp[n1-i-1] != arr(pp->dna0, a2 - i, unsigned char))
		      j = 0 ;
		  if (j) /* ok the whole begining of the reverse SL is in the probe */
		    {
		      /* check that the accepted part SL including an extra donor 'gt' 
		       * is not present in the genome
		       */
		      int kkk = ccp ? strlen ((const char *)ccp) : 0 ;
		      BOOL inGenome = FALSE ;
		      if (kkk < OVLN)
			{
			  sprintf ((char *)buf2, "%c%c%s", A_ | G_,  C_, ccp) ;
			  inGenome = dnaPickMatch (genome, genomeLn, buf2, 0, 0) ;
			}
		      if (! inGenome)
			{
			  goodSuffix = ns + 11 ;
			  score += pp->slBonus ;
			  i = n1 - 1 ;
			  if (i > 0)
			    {
			      x2 = (*x2p) -= i ;
			      a2 = (*a2p) -= i ; 
			      strncpy ((char *)(suffix + 1), (char *)ccp, n1) ;
			    }
			  aliBonus += 0 ; /*  strlen((char *)suffix+1)  */
			  rightAdaptorClip = x2 ;
			}
		    }
		}
	    }
	}
      /* search for the exit adaptor */
      unsigned char **exitAdaptor  = (mm->pair < 0 ? pp->exitAdaptor2 : pp->exitAdaptor) ;
      unsigned int *exitAdaptorFound  = (mm->pair < 0 ? pp->exitAdaptor2Found : pp->exitAdaptorFound) ;
      unsigned  int *exitAdaptorHasN = (mm->pair < 0 ? pp->exitAdaptor2HasN : pp->exitAdaptorHasN ) ;

      if (exitAdaptor && exitAdaptor[0] && !adaptorSuffix &&
	  strlen((char *)suffix + 1) >= 6) /* search for adaptor */
	{
	  int ns, n1 = 0, i1, j1, nE, nN ;


	  for (ns = 0 ; exitAdaptor[ns] && !adaptorSuffix ; ns++)
	    {
	      int aAT = 0, sAT = 0 ;
	      ccp = exitAdaptor[ns] ;
	      int adapLn = strlen ((char *)ccp) ;	
	      if (pp->clipPolyA)
		{
		  for (aAT = 0 ; ccp[aAT] == A_  ; aAT++)
		    ;
		  for (sAT = 0 ; suffix[1+sAT] == A_ ; sAT++)
		    ;
		}
	      if (aAT > sAT)
		continue ;

	      if (adapLn <= aAT)
		continue ;
	      i =  adapLn - aAT  ; if (i>6) i = 6 ;
	      for (j1 = -1, n1 = 0 ; !n1 && (++j1) < 16 ; )
		n1 = dnaPickMatch (ccp + j1+aAT, i, suffix+1+sAT, 1, 2) ;
	      if (n1 > 0)
		{
		  i =  strlen((char *)ccp)  - aAT ;
		  i1 = strlen ((char *)suffix+1+sAT-n1+1) ;
		  n1 = nE = nN = 0 ;
		  if (i >= i1)
		    {
		      if (i1< 6) nE = -1 ;
		      else if (i1 < 8 && exitAdaptorFound[ns] < 1000) { nE = -1 ;nN = 2 ; }
		      else if (i1 < 12 && exitAdaptorFound[ns] < 200) { nE = 0 ; nN = 1 ; }
		      else if (i1 > 20 && exitAdaptorFound[ns] > 1000) { nE = 2 ; nN = 2 ; }
		      else { nE = 1 ; nN = 1 ; }
		    }
		  else
		    {
		      if (
			  i < 16 ||
			  (i < 20 && i < i1 + 5)
			  )
			nE = -1 ;
		      else if (i >= 20)
			{ nE = 1 ; nN = 1 ; }
		      /* else keep nE == 0 */
		    }
		  if (i > i1) i = i1 ;
		  if (i > 20) { i = i1 = 20 ; }
		  nN += exitAdaptorHasN[ns] ;
		  if (nE >= 0)
		    n1 = dnaPickMatch (ccp + aAT, i, probeDna + x2 - n1 - j1 + sAT, nE, nN) ;
		  if (i > 999) i = 999 ; /* to allow 100*ns + i multiplex construction */
		}
	      if (n1 > 0 && n1 < 3)
		{
		  (exitAdaptorFound[ns]) += mm->mult ;
		  adaptorSuffix = 1000 * (1+ns) + i ;
		  /* flush the matching letters of the vector back into the probe suffix */
		  a2 = (*a2p) -= j1 ;
		  x2 = (*x2p) -= j1 ;
		  score -= j1 ;
		  for (cp = suffix+1 +sAT; *cp ; cp++) ;
		  for (--cp ; cp > suffix + j1+sAT ; cp--)
		    *cp = *(cp - j1 ) ;
		  for (i = 0 ; i < j1 ; i++)
		    *(suffix + 1 + i) = *(ccp + i) ;
		  rightAdaptorClip = x2 ;
		  goodSuffix = 3 ;   /* report Adaptor */
		}
	    }
	}
   
      int bufLn = strlen ((char *)suffix+1) ;
      if (! pN && !goodSuffix && 
	  pp->clipPolyA &&
	  bufLn < 200 && 
	  bufLn >= 8) /* search for the bar-code or entry adaptor */
	{
	  int sAT = 0 ;
	  for (sAT = 0 ; suffix[1+sAT] == A_ ; sAT++)
		    ;
	  if (sAT >= 8)
	    goodSuffix = 3 ;
	}

      /* search for a downstream acceptor exon */
      if (!pN && ! goodSuffix && (pp->splice || pp->showOverhang) && !goodSuffix && pp->overhangLength > 3 && strlen ((char *)suffix + 1) >= pp->overhangLength)
	{
	  cp = suffix + 1 + pp->overhangLength ; 
	  cc = *cp ; *cp = 0 ;

	  goodSuffix = 1 ;  /* so it is available to construct introns */
	  /* n =  dnaPickMatch (genome, genomeLn, suffix + 1, 0, 0) ;  */
	  if (pp->splice) /* 2013_08_13, was if(1) */
	    {
	      goodSuffix = 63 ;
	      score += INTRONBONUS ; /* needed to find best score will be removed at export time */
	    }
	  *cp = cc ;
	}

      suffix[0] = goodSuffix ;  /* flag it is an SL a poly-A or an intron candidate */
      
     /* decipher the suffix */
      for (cp = suffix + 1, i = 1, j = adaptorSuffix % 1000 ; *cp ; cp++, i++)
	{
	  *cp = dnaDecodeChar[(int)*cp] ;
	  /* if (i <= j) */	
	  if (goodSuffix > 1 && goodSuffix != 63)
	    *cp = ace_upper(*cp) ;
	}
      *cp = 0 ;
      /*      if (! strcasecmp( suffix + 1, "TGAGCTGC")) invokeDebugger() ; */
      adaptorSuffix /= 1000 ; /* demultiplex */

      /* re append the polyA tail */
      if (mm->probeRightClip)
	{	  
	  j = OVLN - strlen ((char *)suffix + 1) ;
	  bb = 'a' ;
	  if (pp->solid)
	    bb = suffix[strlen((char *)suffix+1)] ; /* true value of the homopolymer */
	  for (i = 0 ; i < mm->probeRightClip && i < j ; i++, cp++)
	    *cp = bb ;  
	  *cp = 0 ;
	}
    }

  /****************************/
  /* Rextend the polyT prefix */
  /****************************/

  if (x1 <= 3 && *a1p > 1 && mm->probeLeftClip && ! pp->solid) /* so many T at the lead of the probe */
    {
      int ok = -1 ;

      i = *a1p - 2  ; /* last non aligned base, could be a T, if yes, extend */
      j = mm->probeLeftClip ;


     /* accept one non T over the last 2 bp  because of systematic error */
      ccp = arrp (dna, i, unsigned char) ;
      if (x1 == 1) ok = 0 ;
      else if (x1 == 2 && mm->probeLeftClip >= 6)
	{ 
	  if (mm->probeLeftClip >= pp->errCost) ok = 1 ; 
	  i-- ; j-- ; 
	}
      else if (x1 == 3 && mm->probeLeftClip >= 6)
	{ 
	  if (mm->probeLeftClip >= pp->errCost && *(ccp-1) == arr (probeDnaArray, x1-3, unsigned char)) 
	    ok = 2 ; 
	  i -= 2 ; j -= 2 ; 
	}
      if (ok >= 0)
	{
	  /* even for solid, we truly know the type of the first bases of the probe */
	  if (i >= 0)
	    for (ccp = arrp (pp->dna0, i, unsigned char) ; j > 0 && i >= 0 && *ccp == T_ ; i--, j--, ccp--) ;
	  
	  /* reextend */
	  i = mm->probeLeftClip - j ;
	  a1 = (*a1p) -= i ;
	  nS += i ;
	  x1 = (*x1p) -= i ;
	  if (0 && ! pTok) score += i ; /* 2014_12_04  adding points for the extension seems to create an unbalance 
					 * if the mRNA does not have its AAAA tail but the genome has it
					 * a contrario it always necessary to reextend in the wrong direction
					 * as the clipping i that case was illicit
					 *  2014_12_05, we cannot do that it favorizes the wrong strand in exact test data
					 */
	  leftAdaptorClip = i ; /* add back to the clipped probe length */
	  if (pTok && j >= 8)
	    {
	      char aA = 'a' ;
	      int nA = j, na = 0 ;
	      unsigned char *cq ;
	      
	      if (a1 > 1)
		{
		  cq =  arrp (pp->dna0, a1 - 2, unsigned char) ;
		  for (i = 0 ; i <= nA ; i++, cq--)
		    {
		      if (*cq == T_) na++ ;
		    }
		}
	      if (10 * na <= 5 * nA)
		{
		  *prefix = 2 ;  /* annotate a polyA */
		  aA = 'A' ;
		}
	      for (i = 0, cp = prefix + 1 ; i < j && i < 58 ; i++, cp++) 
		*cp = aA ; /* export the complement */
	      *cp = 0 ;
	    }
	}
    }

  /**************************/
  /* Rextend the polyA tail */
  /**************************/

  if ( mm->probeRightClip && *a2p < arrayMax (dna) && x2 >= mm->probeLength - 2) /* so many A at the tail of the probe */
    {
      unsigned char bb = A_ ;
      int ok = -1 ;
      
      i = *a2p  ; /* first non aligned base, could be an A, if yes, extend */
      j = mm->probeRightClip ;
      
      /* accept one non A over the last 2 bp  because of systematic error */
      ccp = arrp (dna, i, unsigned char) ;
      if (x2 == mm->probeLength) ok = 0 ;
      else if (x2 == mm->probeLength - 1 && mm->probeRightClip >= 6)
	{ if (mm->probeRightClip >= pp->errCost) ok = 1 ; i++ ; j-- ; }
      else if (x2 == mm->probeLength - 2 && mm->probeRightClip >= 6)
	{ 
	  if (mm->probeRightClip >= pp->errCost && *(ccp+1) == arr (probeDnaArray, x2+1, unsigned char)) 
	    ok = 2 ; 
	  i+=2 ; j-=2 ;
	}

      if (ok >= 0)
	{
	  if (!pp->solid)
	    {
	      if ( i < arrayMax (dna))
		for (ccp = arrp (dna, i, unsigned char) ; j > 0 && i < arrayMax (dna) && *ccp == A_ ; i++, j--, ccp++) ;
	    }
	  else
	    {
	      unsigned char t2 ;
	      
	      bb = arr(pp->dna0,a2 - 2,unsigned char) ;  /* base before the homopolymer, read on the genome */
	      t2 = arr (probeDnaArray, x2 - 1, unsigned char) ; /* last transition on the colored probe */
	      bb = clipAlignSolidProduct (bb, t2) ;  /* true value of the homopolymer */
	      if (bb == A_) /* I truly have a A homopolymer on the probe, not a poly TTT or CCC or GGG */
		{ 
		  if (arr(pp->dna0, a2 - 1, unsigned char) == A_)
		    for (ccp = arrp (pp->dna0, i, unsigned char) ; j > 0 && i < arrayMax (dna) && *ccp == A_ ; i++, j--, ccp++) ;
		  else /* my polyA tail on the probe actually started one base earlier */
		    j = mm->probeRightClip + 1 ;
		}
	      else /* the homolpolymer was of the wrong type, all non t5 bases are errors */
		{
		  j = mm->probeRightClip ;
		}
	    }
	  /* reextend */
	  i = mm->probeRightClip - j ;
	  if (pp->solid) i++ ; /* because the base before is also an A */
	  a2 = (*a2p) += i ;
	  nS += i ;
	  x2 = (*x2p) += i ;
	  if (0 && ! pAok) score += i ; /* 2014_12_04  adding points for the extension seems to create an unbalance 
					 * if the mRNA does not have its AAAA tail but the genome has it
					 * a contrario it always necessary to reextend in the wrong direction
					 * as the clipping in that case was illicit
					 *  2014_12_05, we cannot do that it favorizes the wrong strand in exact test data
					 */
	  rightAdaptorClip = x2 ;

	  if (bb == A_ && j >= 8 && pAok && a2 < arrayMax(pp->dna0))
	    {	 
	      char aA  = 'a' ;
	      int nA = j, na = 0 ;
	      unsigned char *cq =  arrp (pp->dna0, a2, unsigned char) ;

	      for (i = 0 ; i < nA ; i++, cq++)
		{
		  if (*cq == A_) na++ ;
		}
	       if (10 * na <= 5 * nA)
		{
		  *suffix = 2 ;  /* annotate a polyA */
		  aA = 'A' ;
		}
	      for (i = 0, cp = suffix + 1 ; i < j && i < 58 ; i++, cp++) 
		*cp = aA ;
	      *cp = 0 ;
	    }
	}
    }

  if (pp->solid && x1 >= 2)  /* 2011_05_15 the first correct transition indicates the letter before is OK
			      * however we never want to include the leading T in the count
			     */
    {
      *a1p = a1 -= 1 ; 
      *x1p = x1 -= 1 ;
      if (x1 == 0)score++ ;
      nS++ ;
    }

  /***********************************************/
  /* Export the exon ends, plus 2 bordering bp   */
  /***********************************************/

  if (1)
    {
      /* give first the genome sequence of the accepting exon in the direction of the probe */
      /* note that in this function the probe and the target are always forward */
      
      cp = targetPrefix ;

      if (pp->splice || pp->showOverhang || (prefix[0] >= 11 && prefix[0]< 63))
	{
	  /* 2bp on target (probably AG) + 28bp in the exon */
	  i = a1 - 2 - (pp->solid ? 1 : 0) ;  j = 0 ;
	  while (j < OVLN && i <= 0)
	    { *cp++ = '-' ; i++ ; j++ ;}
	  ccp = arrp (pp->dna0, i - 1, unsigned char) ;
	  while (j < OVLN && i <= arrayMax(pp->dna0))
	    { *cp++ = dnaDecodeChar[(int)*ccp++] ; i++ ; j++ ; }
	  if (j < OVLN && ! pp->solid) /* continue on the tag itself */
	    {
	      i = x2 + 1 ;
	      ccp = probeDna + i - 1 ;
	      while (j < OVLN && *ccp)
		{ *cp++ = dnaDecodeChar[(int)complementBase(*ccp++)] ; j++ ; }
	    }
	}

      /* now give complement of the genome sequence of the donor exon */  
      /* 2bp on target (probably GT complemented = AC) + 28bp in the exon */

      *cp++ = '*' ; /* separator */ 
      if (pp->splice || pp->showOverhang || (suffix[0] >= 11 && suffix[0]< 63))
	{
	  i = a2 + 2 ; j = 0 ;
	  while (j < OVLN && i > arrayMax(pp->dna0))
	    { *cp++ = '-' ; i-- ; j++ ; }
	  ccp = arrp (pp->dna0, i - 1, unsigned char) ;
	  while (j < OVLN && i >= 1)
	    { *cp++ = dnaDecodeChar[(int)complementBase(*ccp--)] ; i-- ; j++ ; }
	  if (j < OVLN && ! pp->solid) /* continue on the tag itself */
	    {
	      i = x1 - 1 ;
	      ccp = probeDna + i - 1 ;
	      while (j < OVLN && i >= 1)
		{ *cp++ = dnaDecodeChar[(int)complementBase(*ccp--)] ; i-- ; j++ ; }
	    }
	}

      /****************************/
      /* Export the genome prefix */
      /****************************/
      
      if (pp->showTargetPrefix)
	{
	  int i0 = 0 ;
	  /* prefix */
	  *cp++ = '#' ; /* separator */
	  i = a1 - 29  ;  j = 0 ;
	  while (j < 30 && i <= 0)
	    { *cp++ = '-' ; i++ ; j++ ;}
	  ccp = arrp (pp->dna0, i - 1, unsigned char) ;
	  while (j < 30 && i < a1)
	    { *cp = dnaDecodeChar[(int)*ccp++] ; i++ ; if(!i0 && *cp == 'n') *cp = '-' ; else i0 = 1; cp++ ; j++ ; }
	  /* then the suffix */
	  *cp++ = '\t' ; /* separator */
	  i = a2 + 1 ; j = 0 ;
	  if (i-1 < arrayMax (pp->dna0))
	    {
	      ccp = arrp (pp->dna0, i - 1, unsigned char) ;
	      while (j < 30 && *ccp && i <= arrayMax (pp->dna0))
		{ *cp++ = dnaDecodeChar[(int)*ccp++] ; i++ ; j++ ; }
	    }
	  while (j>0 && *cp == 'n') { j-- ; *cp-- = 0 ;}
	 if (! j) /* ensure something if we are at the end of the sequence */
	    { *cp++ = '-' ;}
	}
      *cp++ = 0 ;
    }

  /***************************/
  /* verify the score, after */
  /* adding pA/pT/SL bonuses */
  /* Rextend the polyA tail  */
  /***************************/

  /* here use bestScore, since pA and INTRON bonuses have already been integrated */
   if (mm->bestX1[0])
    {
      int i ;
      for (i = 0 ; i < NBESTSCORES ; i++)
	if (mm->bestX1[i])
	  {
	    int c1 = x1 > mm->bestX1[i] ? x1 : mm->bestX1[i] ;
	    int c2 = x2 < mm->bestX2[i] ? x2 : mm->bestX2[i] ;
	    int dc = 3 *(c2 - c1 + 1) ;
	    int dx = x2 - x1 + 1 ;
	    if (pp->bestHit && !pp->splice && score < mm->bestScore[i] -pairTampon && dc > dx)
	      return FALSE ;
	  }
    }
 
  if (pp->minAli > x2 - x1 + 1) /* 2012_06_07 needed because of reclipping in case of recognized vector or SL */
    return FALSE ; 
  if (pp->targetMutantMalus && score <  array (pp->originalScores, mm - arrayp (pp->probes, 0, MM), short) + pp->mutantMalus -pairTampon)
    return FALSE ; /* idiot: on ne verra pas les cas d'egalite wild type/mutant */
  score -= pp->targetMutantMalus ;

  if (score < 0)
    return FALSE ;
  if (1)
    {
      int i ;
      int done = -1 ;
      for (i = 0 ; done < 0 && i < NBESTSCORES ; i++)
	if (mm->bestX1[i])
	  {
	    int c1 = x1 > mm->bestX1[i] ? x1 : mm->bestX1[i] ;
	    int c2 = x2 < mm->bestX2[i] ? x2 : mm->bestX2[i] ;
	    int dc = 3 *(c2 - c1 + 1) ;
	    int dx = x2 - x1 + 1 ;
	    int dx2 = mm->bestX2[i] - mm->bestX1[i] + 1 ;
	    if (dc > dx || i == NBESTSCORES - 1) /* good overlap position */
	      {
		if (score >= mm->bestScore[i] && 2*dx > dx2)
		  {
		    mm->bestScore[i] = score ;
		    mm->bestX1[i] = x1 ; mm->bestX2[i] = x2 ;
		    mm->nGeneHits = 0 ;
		  }
		done = i ;
	      }
	  }
	else
	  {
	    mm->bestScore[i] = score ;
	    mm->bestX1[i] = x1 ; mm->bestX2[i] = x2 ;
	    mm->nGeneHits = 0 ;
	    done = i ;
	  }

      if (done >= 0) /* eat up the region below me that I overlap */
	{
	  for (i = done + 1 ; i < NBESTSCORES ; i++)
	    if (mm->bestX1[i])
	      {
		int c1 = x1 > mm->bestX1[i] ? x1 : mm->bestX1[i] ;
		int c2 = x2 < mm->bestX2[i] ? x2 : mm->bestX2[i] ;
		int dc = 3 *(c2 - c1 + 1) ;
		int dx = mm->bestX2[i] - mm->bestX1[i] + 1 ;
		if (dc > dx)       /* good overlap position kill it */
		  {
		    int j ; 
		    for (j = i + 1 ; j < NBESTSCORES ; j++)
		      {
			mm->bestScore[j-1] = mm->bestScore[j] ;
			mm->bestX1[j-1] = mm->bestX1[j] ;
			mm->bestX2[j-1] = mm->bestX2[j] ;
		      }
		    mm->bestScore[j-1] = 0 ;
		    mm->bestX1[j-1] = 0 ;
		    mm->bestX2[j-1] = 0 ;
		  }
	      }
	}
      if (done > 0) /* reorder by shiftting it up */
	{
	  for (i = done - 1 ; i >= 0 ; i--)
	    if (mm->bestScore[i] < score) /* swap */
	      {
		mm->bestScore[done] = mm->bestScore[i] ;
		mm->bestX1[done] = mm->bestX1[i] ;
		mm->bestX2[done] = mm->bestX2[i] ;

		mm->bestScore[i] = score ;
		mm->bestX1[i] = x1 ; mm->bestX2[i] = x2 ;
		
		done = i ; /* new position of this record */
	      }
	}
    }
  if (1 && /* this is a second line garantee, we are already protected by mm->wordHits */
      pp->maxHit && mm->nGeneHits * 30 > pp->maxHit * mm->probeLength) /* assume 30bp per exon */
     return FALSE ;
  if (0 && /* this is a second line garantee, we are already protected by mm->wordHits */
      pp->maxHit && mm->nGeneHits > pp->maxHit ) /* assume 30bp per exon */
     return FALSE ;

  pp->score = score ;
  /* x1  >= mm->probeLength because we added back the trailing A */
  if (x1 + mm->probeLeftClip == 1 && x2 >= mm->probeLength && arrayMax (err) == 0)
    { if (! mm->nPerfectHit) nPerfect++ ;  mm->nPerfectHit++ ; }

  mm->tm = temperature ;
  pp->aliBonusForThisAli = aliBonus ;
  pp->leftAdaptorClipForThisAli = leftAdaptorClip ;
  pp->rightAdaptorClipForThisAli = rightAdaptorClip ;

  ii = iLeft + dxl ; /* for compiler happiness, we keep these variable for debugging purpose */
  return TRUE ;
} /* clipAlignVerifyProbeHit */

/* 
ctgggccctggcgatctggtttgaaaatcatggttctatctcaatccaatcaatggaaaactgaatttcaagtgttatctcccctcaaacaaaaaacaaa
*/
/*********************************************************************/

static void debugOligo (char *title,  unsigned int oligo)
{
  unsigned int x ;
  int nn = 16 ;
  unsigned char cc ;
  fprintf (stderr, "// %s = ", title) ;

  while (nn--) 
    {
      x = oligo >> 2*nn ; x &= 0x3 ; 
      switch (x)
	{
	case 0: cc = 'A' ; break ;
	case 1: cc = 'G' ;break ;
	case 2: cc = 'C' ;break ;
	case 3: cc = 'T' ;break ;
	default: cc = 'N' ;break ;
	}
      fprintf  (stderr, "%c", cc) ;
    }
  fprintf (stderr, "\n") ;
} /* debugOligo */


/*********************************************************************/
/* accept an 8bp polyA polyT or 4 g/c in first 12bp or a sequence whose first 12bp have entropy >= 6 */
static BOOL slideCheckHookEntropy (const unsigned char *hook)
{
  int ss = 0 ;
  int na = 0, nt = 0, ng = 0, nc = 0, nn ;
  const unsigned char *ccp ;
  static BOOL oldNn = -1 ;
  static Array ee = 0 ;

  if (! hook)
    return FALSE ;

  /*  count all letters
   */
  if (1)
    {
      ccp = hook - 1 ; nn = 0 ;
      while (nn++ < 8 && *++ccp)
	switch ((int)(*ccp))
	  {
	    case 'a': na++ ; break ;
	    case 't': nt++ ; break ;
	    case 'g': ng++ ; break ;
	    case 'c': nc++ ; break ;
	    default:  break ;
	  }
    }

  nn = na + nt + ng + nc ;
  if (na == 8 || nt == 8) 
    return TRUE ;
  if (nn < 8)
    return FALSE ;

  /* read now up to 12 */
  if (1)
    {
      while (nn++ < 12 && *++ccp)
	switch ((int)(*ccp))
	  {
	    case 'a': na++ ; break ;
	    case 't': nt++ ; break ;
	    case 'g': ng++ ; break ;
	    case 'c': nc++ ; break ;
	    default:  break ;
	  }
    }

  if (ng+nc > 4) return TRUE ;

  /* lazy calculation of all the logs */
  if (! ee || (nn > 1 && nn != oldNn))
    {
      double s = 0, log4 = log(4.0) ; /* divide by log4 to be in base equivalent */ 
      int j ;

      ee = arrayReCreate (ee, nn, int) ;
      for (j = 1 ; j <= nn ; j++)
        { s = j ; s /= nn ; array (ee, j, int) = - (int) (1000.0 * j * log(s)/log4 ) ; }
      array (ee, 0, int) = 0 ;
    }

  ss = arr (ee, na, int) + arr (ee, nt, int) + arr (ee, ng, int) + arr (ee, nc, int) ;
  ss = (ss + 499)/1000 ;
  return ss >= 6 ? TRUE : FALSE ;
} /* hookEntropy */

/*********************************************************************/
/* accept an 8bp polyA polyT or a sequence whose first 8bp have entropy > 4 */
static BOOL slideCheckExonEntropy16 (CLIPALIGN *pp, unsigned const char *exon)
{
  int ss = 0 ;
  int nn ;

  nn = pp->minEntropy ;
  if (! nn || nn > pp->minAli)
    nn = pp->minAli ;
  ss = exon ? oligoEntropy (exon, strlen ((char*)exon), nn) : 0 ;

  return ss > 0 ? TRUE : FALSE ;
} /* exonEntropy */

/*********************************************************************/
/*********************************************************************/

static int  pExportPairAOrder (const void *a, const void *b)
{
  const PEXPORT *px = (const PEXPORT *)a, *py = (const PEXPORT *)b ;
  int n ;

  if (! px->mm) return  1 ; /* zero probe last */
  if (! py->mm) return -1 ;
  n = px->fragment - py->fragment ; if (n) return -n ; /* singletons last */
  n = px->targetGene - py->targetGene ; if (n) return n ;
  n = px->target - py->target ; if (n) return n ;
  n = px->mm - py->mm ; if (n) return n ;
  n = px->isDown - py->isDown ; if (n) return n ;
  n = px->a1 - py->a1 ; if (n) return n ;
  n = px->score - py->score ; if (n) return - n ; /* bad score last */
  n = px->chain - py->chain ; if (n) return n ;
  n = px->x1 - py->x1 ; if (n) return n ;
  return 0 ;  
} /* pExportPairAOrder */

/*********************************************************************/

static int  pExportPairXOrder (const void *a, const void *b)
{
  const PEXPORT *px = (const PEXPORT *)a, *py = (const PEXPORT *)b ;
  int n ;

  if (! px->mm) return  1 ; /* zero probe last */
  if (! py->mm) return -1 ;
  n = px->fragment - py->fragment ; if (n) return -n ; /* singletons last */
  n = px->targetGene - py->targetGene ; if (n) return n ;
  n = px->target - py->target ; if (n) return n ;
  n = px->mm - py->mm ; if (n) return n ;
  n = px->isDown - py->isDown ; if (n) return n ;
  n = px->x1 - py->x1 ; if (n) return n ;
  n = px->score - py->score ; if (n) return - n ; /* bad score last */
  n = px->chain - py->chain ; if (n) return n ;
  n = px->a1 - py->a1 ; if (n) return n ;
  return 0 ;  
} /* pExportPairXOrder */

/*********************************************************************/

static int pExportHitOrder (const void *a, const void *b)
{
  const PEXPORT *px = (const PEXPORT *)a, *py = (const PEXPORT *)b ;
  int n ;

  if (! px->mm) return  1 ; /* zero probe last */
  if (! py->mm) return -1 ;
  n = px->fragment - py->fragment ; if (n) return -n ; /* singletons last */
  n = px->mm - py->mm ; if (n) return n ;
  n = px->score - py->score ; if (n) return - n ; /* bad score last */
  n = px->targetGene - py->targetGene ; if (n) return n ;
  n = px->target - py->target ; if (n) return n ;
  n = px->isDown - py->isDown ; if (n) return n ;
  n = px->chain - py->chain ; if (n) return n ;
  n = px->x1 - py->x1 ; if (n) return n ;
  n = px->a1 - py->a1 ; if (n) return n ;
  return 0 ;  
} /* pExportHitOrder */

/*********************************************************************/

static int pExportChainOrder (const void *a, const void *b)
{
  const PEXPORT *px = (const PEXPORT *)a, *py = (const PEXPORT *)b ;
  int n ;

  if (! px->mm) return  1 ; /* zero probe last */
  if (! py->mm) return -1 ;
  n = px->fragment - py->fragment ; if (n) return -n ; /* singletons last, used in  clipAlignOptimizeIntronCleanUp  */
  n = px->mm - py->mm ; if (n) return n ;
  n = px->x1 - py->x1 ; if (n) return n ;
  n = px->targetGene - py->targetGene ; if (n) return n ;
  n = px->target - py->target ; if (n) return n ;
  n = px->isDown - py->isDown ; if (n) return n ;
  n = px->a1 - py->a1 ; if (n) return n ;
  n = px->score - py->score ; if (n) return - n ; /* bad score last */
 return 0 ;  
} /* pExportChainOrder */

/*********************************************************************/

static int pExportSamOrder (const void *a, const void *b)
{
  const PEXPORT *px = (const PEXPORT *)a, *py = (const PEXPORT *)b ;
  int n ;

  if (! px->mm) return  1 ; /* zero probe last */
  if (! py->mm) return -1 ;
  n = px->fragment - py->fragment ; if (n) return -n ; /* singletons last, used in  clipAlignOptimizeIntronCleanUp  */
  n = px->mm - py->mm ; if (n) return n ;
  n = px->score - py->score ; if (n) return - n ; /* bad score last */
  n = px->targetGene - py->targetGene ; if (n) return n ;
  n = px->target - py->target ; if (n) return n ;
  n = px->isDown - py->isDown ; if (n) return n ;
  n = px->chain - py->chain ;  if (n) return n ;
  n = px->a1 - py->a1 ; if (n) return n ;
  n = px->x1 - py->x1 ; if (n) return n ;
 return 0 ;  
} /* pExportChainOrder */

/*********************************************************************/

static int pExportDonorOrder (const void *a, const void *b)
{
  const PEXPORT *px = (const PEXPORT *)a, *py = (const PEXPORT *)b ;
  int n ;

  n = px->fragment - py->fragment ; if (n) return -n ; /* singletons last */
  n = px->mm - py->mm ; if (n) return n ;
  n = px->target - py->target ; if (n) return n ;
  n = px->x2 - py->x2 ; if (n) return n ;
  n = px->a2 - py->a2 ; if (n) return n ;
  n = px->score - py->score ; if (n) return - n ; /* bad score last */

  return 0 ;  
} /* pExportDonorOrder */

/*********************************************************************/

static int pExportAcceptorOrder (const void *a, const void *b)
{
  const PEXPORT *px = (const PEXPORT *)a, *py = (const PEXPORT *)b ;
  int n ;

  n = px->fragment - py->fragment ; if (n) return -n ; /* singletons last */
  n = px->mm - py->mm ; if (n) return n ; /* zero probe last */
  n = px->target - py->target ; if (n) return n ;
  n = px->x1 - py->x1 ; if (n) return n ;
  n = px->a1 - py->a1 ; if (n) return n ;
  n = px->score - py->score ; if (n) return - n ; /* bad score last */

  return 0 ;  
} /* pExportAcceptorOrder */

/*********************************************************************/

static int pExportIntronOrder (const void *a, const void *b)
{
  const PEXPORT *px = (const PEXPORT *)a, *py = (const PEXPORT *)b ;
  int n ;

  n = px->a2 - py->a2 ; if (n) return n ;
  n = px->x1 - py->x1 ; if (n) return n ;
  n = px->x2 - py->x2 ; if (n) return n ;
  n = px->s1 - py->s1 ; if (n) return n ;

  n = px->a1 - py->a1 ; if (n) return n ;
  n = px->s2 - py->s2 ; if (n) return n ;
  n = px->score - py->score ; if (n) return - n ; /* bad score last */

  return 0 ;  
} /* pExportIntronOrder */

/*********************************************************************/

static int pExportDoubleIntronOrder (const void *a, const void *b)
{
  const PEXPORT *px = (const PEXPORT *)a, *py = (const PEXPORT *)b ;
  int n ;

  n = px->a2 - py->a2 ; if (n) return n ;
  n = px->x1 - py->x1 ; if (n) return n ;
  n = px->x2 - py->x2 ; if (n) return n ;
  n = px->s1 - py->s1 ; if (n) return n ;

  n = px->a1 - py->a1 ; if (n) return n ;
  n = px->s2 - py->s2 ; if (n) return n ;
  n = px->score - py->score ; if (n) return - n ; /* bad score last */

  return 0 ;  
} /* pExportDoubleIntronOrder */

/*********************************************************************/
static Stack probePairOrderStack = 0 ;

static int probePairOrder (const void *a, const void *b)
{
  const MM *mx = (const MM *)a, *my = (const MM *)b ;
  const char *cp, *cq ;

  cp =  stackText (probePairOrderStack, mx->probeSequence) ;
  cq =  stackText (probePairOrderStack, my->probeSequence) ;

  return lexstrcmp(cp, cq) ;
} /* probePairOrder */

/*********************************************************************/

static void clipAlignDoExportPx (CLIPALIGN *pp, PEXPORT *px, int type)
{
  MM *mm = px->mm ;
  const char *ct = "-\t-" ;
  char buf[256] ;

  buf[0] = 0 ;
  if (pp->showTargetPrefix &&  px->targetPresuffix &&
      ! (px->prefix == -1 &&  px->suffix == -1)
      )
    { /* third subfield in * separated targetPrefix */
      const char *ct1 = dictName (pp->exportDict, px->targetPresuffix) ;
      if (ct1) ct1 = strchr (ct1, '#') ;
      if (ct1 && *(ct1+1)) 
	{
	  ct = ct1 + 1 ;
	  if (px->prefix == -1)
	    {
	      ct1 = strchr (ct, '\t') ;
	      if (ct1)     /* kill the target preffix */
		{
		  buf[0] = '-' ; buf[1] = '\t' ;
		  strncpy (buf + 2, ct1 + 1, 90) ;
		  ct = buf ;
		}
	    }
 	  if (px->suffix == -1)
	    {
	      char *ct2 ; /* writable */
	      strncpy (buf, ct, 90) ;
	      ct2 = strchr (buf, '\t') ;
	      if (ct2) { ct2[1] = '-' ; ct2[2] = 0 ; }
	      ct = buf ;  /* kill the target suffix */
	    }
	}
    }
  else if (pp->UMI && pp->UMI < 250) /* so it can be exported inside buf */
    {
      int i = 0 ;
      int n1 = array (pp->probePairs, mm->fragment, int) ;
      const char *UMIbarCode = 0 ;

      n1 = - n1 ;
      if (n1 >= 2)
	UMIbarCode = stackText (pp->probeStack, n1 - 1) ;
      if (UMIbarCode && UMIbarCode[0])
	for (i = 0 ; i < pp->UMI ; i++)
	  buf[i] = dnaDecodeChar[(int)UMIbarCode[i]] ;      
      else
	buf[i++] = '-' ;
      buf[i++] = '\t' ;
      buf[i++] = '-' ;
      buf[i++] = 0 ;
    }
  /*
    if (0 && mm)  // mm1 points to the matching other half of a paired fragment 
    {
      int n = mm->pair ;
      MM *mm1 ;

      if (n < 0) n = -n ;
      if (n > 1) mm1 = arrp (pp->probes, n - 2, MM) ;
    }
  */
  int   jump = (mm->fragment ? pp->jump5_1 : pp->jump5_2) ;
  switch (type)
    {
    case 0: /* hit */
      aceOutf (pp->ao, "%s\t%05d\t%d\t%d\t%d\t%d\t%d"
	       , stackText (pp->probeStack, mm->probeName)
	       , px->score /* + px->geneStrandMalus */
	       , mm->mult
	       , px->ln  /* tag length_to_align */
	       , px->chainAli ?  px->chainAli : px->ali    /* length aligned */
	       , px->x1 + jump, px->x2 + jump
	       ) ;
      
      aceOutf (pp->ao, "\t%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s"
	       , pp->target_class ? pp->target_class : "-"
	       , px->targetGene ? dictName (pp->targetDict, px->targetGene) : "-"
	       , px->uu ? px->uu : 1  /* unicity */
	       , dictName (pp->targetDict, px->target)
	       , px->a1, px->a2
	       , px->nN, px->nErr

	       , px->errRightVal ? dictName (pp->exportDict, px->errRightVal) : "-"  /* error list in tag coordinates and orientation */
	       , px->errLeftVal ? dictName (pp->exportDict, px->errLeftVal) : "-"   /* error list in target coordinates and orientation */

	       , px->prefix > 0 ? dictName (pp->exportDict, px->prefix) : "-"        /* 5' overhang read in the tag */
	       , px->suffix > 0 ? dictName (pp->exportDict, px->suffix) : "-"            /* 3' overhang read in the tag */

	       , ct
	       ) ;  
      aceOutf (pp->ao, "\t-\t-\t-\t-\tchain %d\t%d\t%d", px->chain,px->s1, px->s2) ; /* first column: distance to mate (later added by bestali */
      aceOutf (pp->ao, "\t%d\t%d", px->mm->gt_ag, px->mm->ct_ac) ;
      if (pp->showProbeSequence)
	{
	  char *cp, *cq, buf[mm->probeLength + 2] ;
	  int i = mm->probeLength ;
	  cp = buf ; cq = stackText (pp->probeStack,  mm->probeSequence) - 1 ;
	  while (i--)
	    *cp++ = dnaDecodeChar[(int)*cq++] ;
	  *cp = 0 ;
	  if (pp->solid)
	    {
	      buf[0] = ace_upper(buf[0]) ;
	      for (cp = buf+1 ; *cp ; cp++)
		switch ((int)*cp)
		  {
		  case 'a': *cp = '0' ; break ;
		  case 'c': *cp = '1' ; break ;
		  case 'g': *cp = '2' ; break ;
		  case 't': *cp = '3' ; break ;
		  default: *cp = '.' ; break ;
		  }
	    }
	  aceOutf (pp->ao, "\t%s", buf) ;
	}
      aceOutf (pp->ao, "\n") ;
      break ;
    case 1: /* suffix : for a forward read, probably pA or donor */
      {
	int a1 = px->a1, a2 = px->a2 ;
	char *type = 0, *strand ;
	char slBuf[5] ;
	strcpy (slBuf, "SL1") ;
	
	strand = a1 < a2 ? "Forward" : "Reverse" ;
	switch ((int)px->type)
	  {
	  case 1:
	    if (0 && pp->splice)
	      return ;
	  case 63:
	    if (a1 > a2 && mm->stranded > 0)
	      { type = "DONOR" ; strand = "Reverse" ; }
	    else if (a1 > a2)
	      { type = "ACCEPTOR" ; strand = "Forward" ; }
	    else if (a1 < a2 && mm->stranded < 0)
	      { type = "ACCEPTOR" ; strand = "Reverse" ; }
	    else
	      { type = "DONOR" ; strand = "Forward" ; }
	    break ;
	  case 2:
	    type = "pA" ;
	    break ;
	  case 3:
	    type = "ExitAdaptor" ;
	    break ;
	  case 11:
	  case 12:
	  case 13:
	  case 14:
	  case 15:
	  case 16:
	  case 17:
	  case 18:
	  case 19:
	    type = slBuf ; slBuf[2] = px->type - 10 + '0' ; slBuf[3] = 0 ;
	    strand = a1 > a2 ? "Forward" : "Reverse" ;
	    break ;
	  case 20:
	  case 21:
	  case 22:
	  case 23:
	  case 24:
	    type = slBuf ; slBuf[2] = '1' ; slBuf[3] = px->type - 20 + '0' ; slBuf[4] = 0 ;
	    strand = a1 > a2 ? "Forward" : "Reverse" ;
	    break ;
	  default:
	    break ;
	  }	   
	aceOutf (pp->ao,  "%s\t%05d\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%s\t%d\n"
		 , stackText (pp->probeStack, px->mm->probeName)
		 , px->score, mm->mult
		 , px->x2+jump, px->x1+jump
		 , type
		 , dictName (pp->targetDict, px->target)
		 , strand
		 , a2
		 , px->targetPresuffix ? dictName (pp->exportDict, px->targetPresuffix) : "-"
		 , px->presuffix ? dictName (pp->exportDict, px->presuffix) : "-"
		 , px->score
		 ) ; 
      }
      break ;
    case 2: /* prefix : for a forward read, probably SL or acceptor */
      {
	int a1 = px->a1, a2 = px->a2 ;
	char *type = 0, *strand ;
	char slBuf[5] ;
	strcpy (slBuf, "SL1") ;
	
	strand = a1 < a2 ? "Forward" : "Reverse" ;
	switch ((int)px->type)
	  {
	  default:
	    break ;
	  case 1:
	    return ; 
	  case 63:
	    if (a1 > a2 && mm->stranded > 0)
	      { type = "ACCEPTOR" ; strand = "Reverse" ; }
	    else if (a1 > a2)
	      {  type = "DONOR" ; strand = "Forward" ; } 
	    else if (a1 < a2 && mm->stranded < 0)
	      { type = "DONOR" ; strand = "Reverse" ; }
	    else
	      { type = "ACCEPTOR" ; strand = "Forward" ; }
	    break ;
	  case 2:
	    type = "pT" ; /* a1<a2, i see TTT, => polyA on Reverse strand */
	    break ;
	  case 3:
	    type = "EntryAdaptor" ;
	    break ;
	  case 11:
	  case 12:
	  case 13:
	  case 14:
	  case 15:
	  case 16:
	  case 17:
	  case 18:
	  case 19:
	    type = slBuf ; slBuf[2] = px->type - 10 + '0' ; slBuf[3] = 0 ;
	    break ;
	  case 20:
	  case 21:
	  case 22:
	  case 23:
	    type = slBuf ; slBuf[2] = '1' ; slBuf[3] = px->type - 20 + '0' ; slBuf[4] = 0 ;
	    break ;
	  }

	if (px->score)
	  aceOutf (pp->ao,  "%s\t%05d\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%s\t%d\n"
		   , stackText (pp->probeStack, mm->probeName)
		   , px->score, mm->mult
		   , px->x1+jump, px->x2+jump
		   , type
		   , dictName (pp->targetDict, px->target)
		   , strand
		   , a1
		   , px->targetPresuffix ? dictName (pp->exportDict, px->targetPresuffix) : "-"
		   , px->presuffix ? dictName (pp->exportDict, px->presuffix) : "-"
		   , px->score
		   ) ; 
      }
      break ;

    case 3: /*INTRONS */
      if (px->score)
	{
	  aceOutf (pp->ao, "%s%s\t%s\t%s\t%05d\t%05d\t%s\t%05d\t%05d\t%s\t%d\t%d\t%d"
		   , px->nN & 0x2 ? "NS" : ""
		   , px->nN < 0x4 ? "INTRON" : (px->nN < 0x8 ? "BIGINTRON" : (px->nN < 0x10 ? "HUGEINTRON" : (px->nN < 0x20 ? "DELETION" : (px->nN < 0x40 ? "DUPLICATION" : "TRANSLOCATION"))))
		   , px->nN >= 0x40 ? "-" : (px->nN & 0x1 ? "Reverse" : "Forward")
		   , dictName (pp->targetDict, px->s1)
		   , px->x1+jump, px->a1
		   , dictName (pp->targetDict, px->s2)
		   , px->a2, px->x2+jump
		   , px->foot
		   , px->ali
		   , px->score 
		   , px->pairScore
		   ) ;
#ifdef JUNK
	  if (px->nN >= 0x10 && px->nErr)
	    {  /* export VCF name */
	      if (px->nN < 0x20) /* deletion */
		aceOutf (pp->ao, "\t%s:%d:%s:%c"
			 , dictName (pp->targetDict, px->s1)
			 , px->a1 - px->nErr + 1
			 , px->suffix ? dictName (pp->targetDict, px->suffix) : ""
			 , px->suffix ? (dictName (pp->targetDict, px->suffix))[0] : ""
			 ) ; 
	      else if (px->nN < 0x40) /* duplication */
		aceOutf (pp->ao, "\t%s:%d:%c:%s"
			 , dictName (pp->targetDict, px->s1)
			 , px->a2 - px->nErr + 1
			 , px->suffix ? (dictName (pp->targetDict, px->suffix))[0] : ""
			 , px->suffix ? dictName (pp->targetDict, px->suffix) : ""
			 ) ; 
	  }
#endif
	  aceOutf (pp->ao, "\n") ;
	}
      break ;
    case 4: /* DOUBLE INTRONS */
      if (px->score)
	{
	  aceOutf (pp->ao, "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n"
		   , "DOUBLEINTRON"
		   , dictName (pp->targetDict, px->target)
		   , px->a1, px->a2
		   , px->x1+jump, px->x2+jump
		   , px->s1+jump, px->s2+jump
		   , px->score
		   , px->foot 
		   , (px->mm &&  px->mm->probeName ? stackText (pp->probeStack, px->mm->probeName) : "")
		   ) ; 	
	}
      break ;
    }
} /* clipAlignDoExportPx */

/*********************************************************************/
/* export in vectorized mode
 * add the scores and stack the coordicates of the pieces a la Jim Kent 
 */
static void clipAlignDoExportVectorizedPx (CLIPALIGN *pp, PEXPORT *px, int type)
{
  clipAlignDoExportPx (pp, px, type) ;
}

/*********************************************************************/

static int clipAlignExportPx (CLIPALIGN *pp, int type)
{
  PEXPORT *px ;
  int i, nn = 0 ;
  BigArray hits = 0 ;
  MM *oldMm = 0 ;
  
  switch (type)
    {
    case 0:
      fprintf (stderr, "Got %ld hits %d readPlus %d readMinus %.2f%%\n"
	       , bigArrayMax (pp->exportHits) - 1
	       , pp->nReadPlus, pp->nReadMinus
	       , 100.0 * pp->nReadPlus/(1 + pp->nReadPlus + pp->nReadMinus)
	       ) ;
      hits = pp->exportHits ;
      break ;
    case 1:
      fprintf (stderr, "Got %ld donors\n",  bigArrayMax (pp->exportDonors)) ;
      hits = pp->exportDonors ;
      break ;
    case 2:
      fprintf (stderr, "Got %ld acceptors\n",  bigArrayMax (pp->exportAcceptors)) ;
      hits = pp->exportAcceptors ; 
      break ;
    case 3: 
      fprintf (stderr, "Got %ld introns\n",  bigArrayMax (pp->exportIntrons)) ;
      hits = pp->exportIntrons ;
      break ;
    case 4: 
      fprintf (stderr, "Got %ld double-introns\n",  bigArrayMax (pp->exportDoubleIntrons)) ;
      hits = pp->exportDoubleIntrons ;
      break ;
    }

  if (hits)
    for (i = 0, px = bigArrayp (hits, 0, PEXPORT) ; i < bigArrayMax (hits) ; i++, px++)
      {
	/* vectorize ? */
	if (px->score > 0)
	  { 
	    nn++ ;
	    if (px->mm != oldMm)
	      pp->nAlignedReads++ ;
	    pp->nErrors += px->nErr ;
	    pp->nAlignedBases += px->x2 - px->x1 + 1 ;
	    oldMm = px->mm ;
	    if (type == 0 && pp->vectorize )
	      clipAlignDoExportVectorizedPx (pp, px, type) ; 
	    else
	      clipAlignDoExportPx (pp, px, type) ; 
	  }
      }
  return nn ;
} /* clipAlignExportPx */

/*********************************************************************/
/* hits is  in pExportSamOrder
 * so organized by  read/target/chain in the orientation of the target 
 * as befits a good Cuban CIGAR
 */
static char *clipAlignExportOneSamCigar (CLIPALIGN *pp,PEXPORT *px, int di, vTXT cigar, Array cigarettes, Array dnaShort, Array err, int *nErrp, AC_HANDLE h)
{
  int ii, nErr = 0 ;
  PEXPORT *py ;
  int a1, a2, b2 = 0, x1, x2, y2 ;
  const char *ccp, *errors ;
  int isDown = 1 ;
  int ln = 0 ;

  vtxtClear (cigar) ;

  for (py = px, ii = 0 ; ii < di ; ii++, py++)
    {
      BOOL isGtAg = FALSE ;
      BOOL isGcAg = FALSE ;
      int oldA, ddN = 0, ddA = 0 ;
      int da1 = 0, da2 = 0, intron = 0, intronLn = 0 ;
      a1 = py->a1 ;
      a2 = py->a2 ;
      x1 = py->x1 ;
      x2 = py->x2 ;

      if (a1 < a2) 
	{
	  intron = py->intron ;
	  if (intron)
	    {
	      PEXPORT *pi = bigArrp (pp->exportIntrons, intron - 1, PEXPORT) ;
	      isGtAg = !strcmp (pi->foot, "gt_ag") ;
	      isGcAg = !strcmp (pi->foot, "gc_ag") ;
	      if (pi->a2 > pi->a1)
		{
		  da1 = pi->a2 - a1 ;
		  x1 += da1 ;
		  intronLn = pi->a2 - pi->a1 - 1 ;
		}
	      else
		{
		  da1 = pi->a1 - a1 ;
		  x1 += da1 ;
		  intronLn = pi->a1 - pi->a2 - 1 ;
		}
	    }
	  intron = (ii < di - 1) ? py[1].intron : 0 ;
	  if (intron)
	    {
	      PEXPORT *pi = bigArrp (pp->exportIntrons, intron - 1, PEXPORT) ;
	      if (pi->a2 > pi->a1)
		{
		  da2 = pi->a1 - a2 ;
		  x2 += da2 ;
		}
	      else
		{
		  da2 = pi->a2 - a2 ;
		  x2 += da2 ;
		}
	    }
	}
      else
	{
	  intron = (ii > 0) ? py[-1].intron : 0 ;
	  if (intron)
	    {
	      PEXPORT *pi = bigArrp (pp->exportIntrons, intron - 1, PEXPORT) ;
	      isGtAg = !strcmp (pi->foot, "gt_ag") ;
	      isGcAg = !strcmp (pi->foot, "gc_ag") ;
	      if (pi->a2 > pi->a1)
		{
		  da2 = pi->a2 - a2 ;
		  x2 -= da2 ;
		  intronLn = pi->a2 - pi->a1 - 1 ;
		}
	      else
		{
		  da2 = pi->a1 - a2 ;
		  x2 -= da2 ;
		  intronLn = pi->a1 - pi->a2 - 1 ;
		}
	    }
	  intron = py->intron ;
	  if (intron)
	    {
	      PEXPORT *pi = bigArrp (pp->exportIntrons, intron - 1, PEXPORT) ;
	      if (pi->a2 > pi->a1)
		{
		  da1 = pi->a1 - a1 ;
		  x1 -= da1 ;
		}
	      else
		{
		  da1 = pi->a2 - a1 ;
		  x1 -= da1 ;
		}
	    }
	}
      a1 += da1 ;
      a2 += da2 ;

      errors = py->errLeftVal ? dictName (pp->exportDict, py->errLeftVal) : 0 ;   /* error list in target coordinates and orientation */
      
      if (a1 > a2)
	{
	  int a0 ;
	  a0 = a1 ; a1 = a2 ; a2 = a0 ; a0 = x1 ; x1 = x2 ; x2 = a0 ;
	  isDown = -1 ;
	  if (ii > 0) /* look for an overlap on the read */
	    {
	      ddN = y2 - 1 - x1 ;
	    }
	}
      else if (ii > 0) /* look for an overlap on the read */
	{
	  ddN = x1 - y2 - 1 ;
	}
      b2 += intronLn ;
      ddA = a1 - b2 - 1 ; 


      if (ii == 0) /* gap en tete de sequence */
	{
	  ddN = 0 ; ddA = 0 ;
	  if (isDown == 1 && x1 > 1)
	    { vtxtPrintf (cigar, "%dS", x1 - 1) ;  ln += x1 - 1 ; }
	  else if (isDown == -1 && x1 < arrayMax (dnaShort))
	    { vtxtPrintf (cigar, "%dS", arrayMax (dnaShort) - x1) ; ;  ln += arrayMax (dnaShort) - x1 ; }
	}
      if (intronLn)
	{
	  vtxtPrintf (cigar, "%d%c", intronLn, isGtAg ? 'N' :(isGcAg ? 'n' : 'D')) ;
	}

      if (ddN < 0) /* overlap in x , on detruit cet overlap */ 
	{
	  if (1) { a1 -= ddN ; x1 -= ddN ; ddA -= ddN ; }
	  ddN = 0 ;
	}
      if (ddN > 0 && ddA >= ddN)
	{  /* on alonge l'exon precedent avec des erreurs */ 
	  vtxtPrintf (cigar, "%dX", ddN) ; ln += ddN ;
	  ddA -= ddN ; 
	  b2 += ddN ;
	  y2 += ddN * isDown ;
	  ddN = 0 ;
	}
      if (ddN > 0 && ddA < 0)
	{
	  ddN -= ddA ;
	  a1 -= ddA ;
	  x1 -= ddA * isDown ;
	  ddA = 0 ;
	}
      if (ddN > 0 && ddA > 0 && ddN >= ddA)
	{  /* on alonge l'exon precedent avec des erreurs */ 
	  vtxtPrintf (cigar, "%dX", ddA) ; ln += ddA ;
	  b2 += ddA ;
	  y2 += ddA * isDown ;
	  ddN -= ddA ;
	  ddA = 0 ;
	}
      if (ddN > 0 && ddA == 0)
	{  /* insertion */
	  vtxtPrintf (cigar, "%dI", ddN) ; ln += ddN ;
	  y2 += ddN * isDown ;
	  ddN = 0 ;
	}

      if (ddA < 0) /* duplication */
	{
	  vtxtPrintf (cigar, "%dI", -ddA) ; ln += -ddA ;
	  a1 -= ddA ; x1 -= ddA * isDown ; ddA = 0 ;
	}
      else if (0 && ddA > 30 && pp->strategy == STRATEGY_RNA_SEQ) /* intron */
	{
	  vtxtPrintf (cigar, "%dN", ddA) ; ddA = 0 ;
	}
      else if (ddA > 0)  /* deletion */
	{
	  vtxtPrintf (cigar, "%dD", ddA) ; ddA = 0 ;
	}
 

      oldA = a1 - 1 ;

      if (errors)
	{
	  int a, n, nSub = 0 ;
	  ccp = errors - 1 ;
	  
	  while (ccp[1] && ccp[2])
	    {
	      a = 0 ; while (*++ccp != ':') { a = 10 * a + (*ccp - '0') ; } /* coord of the error */ 
	      if (a < oldA) /* clip this error */
		{
		  while (ccp[1] && *ccp != ',')
		    ccp++ ;
		  continue ;		  
		}
	      
	      n = 0 ; ccp++ ; /* jump the : */
	      if (*ccp == '*') ccp++ ;
	      if (a - oldA  + (*ccp == '+' ? 1 : 0) > 1)
		{
		  int k ;
		  if (nSub)
		    { vtxtPrintf (cigar, "%dX", nSub) ;  ln += nSub ; }
		  nSub = 0 ;
		  k =  a - oldA - 1 + (*ccp == '+' ? 0 : 0) ;
		  if (k)
		    {
		      vtxtPrintf (cigar, "%d=",k) ;
		      ln +=  k ;
		    }
		  oldA  = a ;
		}
	      if (n == 0)
		{
		  while (*ccp == '+') { n++ ; ccp++; }
		  if (n) 
		    { 
		      if (nSub)
			{ vtxtPrintf (cigar, "%dX", nSub) ;  ln += nSub ; }
		      nSub = 0 ; 	
		      oldA = a ;
		      vtxtPrintf (cigar, "%dI", n) ; oldA += 0*n - 1 ; ccp += n ;
		      ln += n ;
		    }
		}
	      if (n == 0)
		{
		  while (*ccp == '-') { n++ ; ccp++; }
		  if (n) 
		    {
		      if (nSub)
			{ vtxtPrintf (cigar, "%dX", nSub) ; ln += nSub ; }
		      nSub = 0 ;
		      oldA = a ;
		      vtxtPrintf (cigar, "%dD", n) ; oldA += n - 1 ; ccp += n ; 
		    }
		}
	      if (n == 0) /* substitution */
		{
		  if (a == oldA  + 1)
		    ;
		  else
		    {
		      if (nSub)
			{ vtxtPrintf (cigar, "%dX", nSub) ; ln +=  nSub ; } 
		      nSub = 0 ;
		      if (a - oldA) 
			{ vtxtPrintf (cigar, "%d=", a - oldA - 1) ; ln += a - oldA - 1 ; }
		    }
		  oldA = a ; 
		  nSub++ ;
		  ccp += 3 ; n = 1 ;
		}
	      nErr += n ; if (! *ccp) break ; 
	    }
	  if (nSub)
	    { vtxtPrintf (cigar, "%dX", nSub) ; ln +=  nSub ; }
	  nSub = 0 ;
	}
      if (a2 > oldA)
	{ vtxtPrintf (cigar, "%d=", a2 - oldA) ; ln +=  a2 - oldA ; }

      y2 = x2 ;
      b2 = a2 ;
    }  

  /* gap en queue */
  if (isDown == 1 && x2 <  arrayMax (dnaShort))
    { vtxtPrintf (cigar, "%dS", arrayMax (dnaShort) - x2) ;  ln += arrayMax (dnaShort) - x2 ; }
  else if (isDown == -1 && x2 > 1)
    { 
      int k = arrayMax (dnaShort) - (ln + x2 - 1) ;
      x2 += k ;
      if (x2 > 1)
	vtxtPrintf (cigar, "%dS", x2 - 1) ;  
      ln += x2 - 1 ; 
    }

  samCheckCigar(0, vtxtPtr(cigar), cigarettes, px->a1, arrayMax (dnaShort)) ;

  *nErrp = nErr ;
  return vtxtPtr (cigar) ;
} /* clipAlignExportOneSamCigar */

/*********************************************************************/
/* count the hit lines summarized in this sam line */ 
static int clipAlignExportOneSam (ACEOUT ao, CLIPALIGN *pp
				  , vTXT cigar, Array cigarettes
				  , Array dnaShort, Array err
				  , BigArray hits
				  , long int ii0, long int iMax)
{
  AC_HANDLE h = 0 ;
  int flag, di = 1 ;					
  long int j ;
  PEXPORT *pz = 0, *py, *px = bigArrayp (hits, ii0, PEXPORT) ;
  int nerr = 0 ;
  int a1 = px->a1, a2 = px->a2 ;
  int b1 = 0, b2 = 0 ;
  int target = px->target ;
  int chain = px->chain ;
  BOOL isDown = TRUE ;
  BOOL mateIsDown = TRUE ;
  BOOL isPrimary = TRUE ;
  BOOL isSecondary = FALSE ;
  BOOL isSupplementary = FALSE ;
  BOOL isOrphan = FALSE ; 
  char *seq ;
  MM *mm = px->mm ;
  int fragment = mm ? mm->fragment : 0 ;
  MM *mate = 0 ;

  if (!px->target || ! px->score || ! px->pairScore)
    return 1 ;
  h = ac_new_handle () ;

  if ((a1 < a2 && px->x1 < px->x2) || (a1 > a2 && px->x1 > px->x2))
    isDown = TRUE ;
  else
    isDown = FALSE ;

  if (a1 > a2)
    { int a0 = a1 ; a1 = a2 ; a2 = a0 ; }

  /* the chain starting at px including di lines will be reported as a single CIGAR */
  for (di = 0, j = ii0, py = px ; j < iMax && py->mm == mm && py->target == target && py->chain == chain ; j++, py++)
    {
      if (py->a1 < a1) a1 = py->a1 ;
      if (py->a2 < a1) a1 = py->a2 ;
      di++ ;
    }

  /* check f I am a primary alignment, i.e. first alignment of the current read */
  for (j = ii0, py = px ; j >= 0 && py->mm == mm && py->target == target && py->chain == chain ; j--, py--)
    ;
  if (j >= 0 && py->mm == mm) /* i changed chain or target */
    {
      isPrimary = FALSE ;
      if (px->uu == 1)
	isSupplementary = TRUE ;
      else
	isSecondary = TRUE ;
    }

  /* proceed to the next read of the same fragment */
  if (! strcmp (stackText (pp->probeStack, mm->probeName), "seq.144399>"))
    invokeDebugger () ;
  if (isPrimary && mm->pair)
    {
      /* try upwards */
      mate = 0 ; pz = 0 ;
      for (j = ii0, py = px ; j >= 0 && py->mm->fragment == fragment ; j--, py--)
	;
      j++ ; py++ ;
      if (j >= 0 && py->mm->fragment == fragment && py->mm != mm)
	{ pz = py ; mate = pz->mm ; } /* first hit of other read */
      
      if (mate == 0)
	{ /* try downwards */
	  for (j = ii0, py = px ; j < iMax && py->mm == mm ; j++, py++)
	    ;
	  if (j < iMax && py->mm && py->mm->fragment == fragment && py->mm != mm)
	    { pz = py ; mate = pz->mm ; } /* first hit of other read */
	}
      if (mate) /* select mateTarget and matePosition and fragment length */
	{
	  b1 = b2 = pz->a1 ;
	  if (pz->x1 < pz->x2 && pz->a1 < pz->a2)
	    mateIsDown = TRUE ;
	  else
	    mateIsDown = FALSE ;
	  for (py = pz ; j < iMax && py->mm == pz->mm && py->target == pz->target && py->chain == pz->chain ; j++, py++)
	    {
	      if (py->a1 < b1) b1 = py->a1 ;
	      if (py->a2 < b1) b1 = py->a2 ;
	      if (py->a1 > b2) b2 = py->a1 ;
	      if (py->a2 > b2) b2 = py->a2 ;
	    }
	}
    }
 
   if (isPrimary && mm->pair && ! mate)  
     isOrphan = TRUE ;

  /* there are 11 mandatory columns */
  aceOutf (ao, "%s", stackText (pp->probeStack, mm->probeName)) ;
  flag = 0 ;
  if (1) 
    {   /* SAM FLAGS */
      /* this read belongs to a pair */
      if (pp->hasPairs && mm->pair) flag |= 0x1 ;
      /* all segments are well and compatibly aligned */
      if (pp->hasPairs && mm->pair && mate) flag |= 0x2 ;
      /* unmapped read */
      if (! px->target) flag |= 0x4 ;  
      /* orphan pair-mate is unmapped */
      if (pp->hasPairs && mm->pair && isOrphan) flag |= 0x8 ;
      /* the seq is given in the orientation of the genome
       * the flags say if the shown read is the complement of 
       * read found on the machine fastq file 16 == 0x10,  32 == 0x20
       */
      if (! isDown) flag |= 0x10 ;
      if (pp->hasPairs && mm->pair && mate && ! mateIsDown) flag |= 0x20 ;
      /* i am read 1 of a pair or first or central part of a composite */
      if (pp->hasPairs && mm->pair > 1) flag |= 0x40 ;
      /* i am read 2 of a pair, or last or central part of a composite 128 == 0x80*/
      /* notice that central read parts are 0x4 | 0x8 == 0xb0 */ 
      if (pp->hasPairs && mm->pair < 0) flag |= 0x80 ;
      /* multi-mapping : all have this flag except first one 256 == 0x100 */
      if (isSecondary) flag |= 0x100 ;
      /* read not passing quality filters 512 = 0x200 */
      if (px->score < 10) flag |= 0x200 ;
      /* i am a PCR duplicate: never used by magic 1024 == 0x400 */
      if (0) flag |= 0x400 ;
      /*  supplementary alignment  2048 == 0x800
       *  divers portion du meme read sont incompatibles
       *  all have this flag except the first one
       */
      if (isSupplementary) flag |= 0x800 ;
    }

  aceOutf (ao, "\t%d", flag) ;
  aceOutf (ao, "\t%s", dictName (pp->targetDict, px->target)) ;
  aceOutf (ao, "\t%d", a1) ;
  { /* do we trust the ali */
    int z = 0 ;
    if (
	(! mm->pair || px->pairScore > px->score) &&
	10*px->ali >= 9*px->ln && 100 * px->nErr < px->ali
	)
      z = 33 ;
    else if (px->nErr < 4)
      z = 20 ;
    else
      z = 10 ;
    z = px->pairScore ;
    aceOutf (ao, "\t%d", z) ;
  }
  
  seq =  stackText (pp->probeStack, mm->probeSequence) ;
  array (dnaShort, strlen (seq) + mm->probeLeftClip + mm->probeRightClip, char) = 0 ; /* zero terminate */
  arrayMax (dnaShort) = strlen (seq) + mm->probeLeftClip + mm->probeRightClip  ;
  if (mm->probeLeftClip)
    memset (arrp (dnaShort, 0, char), T_,  mm->probeLeftClip) ;
  memcpy (arrp (dnaShort, mm->probeLeftClip, char), seq, strlen (seq)) ;
  if (mm->probeRightClip)
    memset (arrp (dnaShort, mm->probeLeftClip + strlen (seq), char), A_,  mm->probeRightClip) ;

  if (! isDown)
    reverseComplement (dnaShort) ;
  
/* CIGAR */
  aceOutf (ao, "\t%s", clipAlignExportOneSamCigar (pp, px, di, cigar, cigarettes, dnaShort, err, &nerr, h)) ;
  if (1)
    aceOutf (ao, "\t%s", isPrimary && pz && pz->target ? (px->target == pz->target ? "=" :dictName (pp->targetDict, pz->target)) : "0") ; /* RNEXT name of mate or next chain of same read */
  else
    aceOutf (ao, "\t%s", isPrimary && mate ? (mm == mate ? "=" : stackText (pp->probeStack, mate->probeName)) : "0") ; /* RNEXT name of mate or next chain of same read */
  aceOutf (ao, "\t%d", b1) ; /* PNEXT pmate position or next string if on same target */
  aceOutf (ao, "\t%d", b2 ? (a1 < b2 ? b2 - a1 + 1 : b1 - a2 - 1) : 0) ; /* TLEN observed template length */
  
  if (0 && ! isPrimary)
    aceOutf (ao, "\t*") ;
  else
    {
      dnaDecodeArray (dnaShort) ;
      aceOutf (ao, "\t%s", arrp (dnaShort, 0, char)) ;
    }
  
  if (! isPrimary)
    aceOutf (ao, "\t*") ;
  else
    {
      if (isDown)
	aceOutf (ao, "\t*") ;  /* QUAL ascii de phread scale quaity + 33 ou * */
      else
	aceOutf (ao, "\t*") ; /* we must reverse the order of the QUALs */
    }
  /* end of the 11 mandatory columns */
  
  aceOutf (ao, "\tNH:i:%d", px->uu) ;
  aceOutf (ao, "\tAS:i:%d", px->score) ;
  aceOutf (ao, "\tNM:i:%d", nerr) ;
  aceOutf (ao, "\tmt:i:%d", mm->mult) ;
  /* aceOutf (ao, "\tMD:Z:0") ; */
  
  aceOut (ao, "\n") ;
  ac_free (h) ;
  return di ;
}  /* clipAlignExportOneSam */

/*********************************************************************/

static int clipAlignExportSam (CLIPALIGN *pp, const char *commandBuf)
{ 
  long int nn = 0 ;
  BigArray hits = pp->exportHits ;
  
  if (hits)
    {
      AC_HANDLE h = ac_new_handle () ;
      Array dnaShort = arrayHandleCreate (1024, char, h) ;
      Array err = arrayHandleCreate (1024, A_ERR, h) ;
      Array cigarettes = arrayHandleCreate (1024, SAMCIGAR, h) ;
      ACEOUT ao = aceOutCreate (pp->outFileName, ".sam", pp->gzo, h) ;
      long int i, iMax = bigArrayMax (hits) ;
      int di = 1 ;
      vTXT cigar = vtxtHandleCreate (h) ;
      aceOutf (ao, "@HD VN:1.5\tSO:queryname\n") ;
      aceOutf (ao, "@PG ID:1\tPN:Magic\tVN:%s", VERSION) ;
      aceOutf (ao, "\tCL:%s", commandBuf) ;
      aceOut (ao, "\n") ;

      bigArraySort (hits, pExportSamOrder) ;
      for (i = 0 ; i < iMax ; i += di)
	di = clipAlignExportOneSam (ao, pp, cigar, cigarettes, dnaShort, err, hits, i, iMax) ;

      fprintf (stderr, "// %s Exported %ld ali in sam file %s\n", timeShowNow(), nn, aceOutFileName (ao)) ;
      ac_free (ao) ;
      ac_free (h) ;
    }

  return nn ;
} /*  clipAlignExportSam */

/*********************************************************************/

static int clipAlignRestoreAcceptorScores (CLIPALIGN *pp)
{
  PEXPORT *up, *vp ;
  int ii, iiMax = bigArrayMax (pp->exportDonors) ;
  int jj, jjMax = bigArrayMax (pp->exportAcceptors) ;

  /* restore the score of the donors/acceptors */
  if (iiMax)
    for (ii = 0, up =  bigArrp (pp->exportDonors, 0, PEXPORT) ; ii < iiMax ; up++, ii++)
      if (up->score < 0) up->score = - up->score ;
  if (jjMax)
    for (jj = 0, vp =  bigArrp (pp->exportAcceptors, 0, PEXPORT) ; jj < jjMax ; vp++, jj++)
      if (vp->score < 0) vp->score = - vp->score ;
  return 0 ;
} /*clipAlignRestoreAcceptorScores */

/*********************************************************************/

static int clipAlignConstructIntrons (CLIPALIGN *pp, int isDown, BOOL singleTarget)
{
  PEXPORT *up, *vp, *wp, *px = 0 ;
  int i, ii, jj, kk, jjA = 0, nIntrons = 0, fragment = 0 ;
  int da, dx, bestda, bestdx, bestjj, bestk, k1, k2 ;
  int a2, b1, x2, y1 ;
  int nnTest = 0 ;
  BigArray exportIntrons       = (singleTarget ? pp->geneIntrons : pp->exportIntrons) ;
  BigArray exportDonors        = (singleTarget ? pp->geneDonors : pp->exportDonors) ;
  BigArray exportAcceptors     = (singleTarget ? pp->geneAcceptors : pp->exportAcceptors) ;
  int iiMax = arrayMax (exportDonors) ;
  int jjMax = arrayMax (exportAcceptors) ;
  BigArray aa ;
  MM *mm ;
  const char *exon1, *exon2, *hook1, *hook2, *cp, *ccp, *ccq ;
  int dx1, dx2, overlap ;
  char *cr, cc ;
  BOOL reverse ;
  char f1[6], f2[6], foot[6] ;
  int shorty = pp->seedLength + pp->seedShift + 1 ;
  BOOL debug = FALSE ;
  int stranded = 0, globallyStranded = 0 ;
  if (debug) showPexport (pp, 2) ;
  if (debug) showPexport (pp, 3) ;

  if (! bigArrayMax (exportDonors) || ! bigArrayMax(exportAcceptors))
    return 0 ;
  if (!singleTarget)
    return 1 ;

  if (1) globallyStranded = 0 ; /* at this stage we ignore the strand of the read, because we have not counted yet */
  if (0)
    {	
      if (pp->nIntronPlus  > 1000 && pp->nIntronPlus  > 20 * pp->nIntronMinus)
	globallyStranded = 1 ;
      else if (pp->nIntronMinus > 1000 && pp->nIntronMinus > 20 * pp->nIntronPlus)
	globallyStranded = -1 ;
    }
  /* 20250418: set k2=0 below so that minority antisens introns are reported  */
  if (pp->stranded > 0)   globallyStranded = 1 ;
  if (pp->stranded < 0)   globallyStranded = -1 ; 
  
  aa = exportIntrons ;
  ii = jj = 0 ; nIntrons = arrayMax (aa) ;
  up = bigArrp (exportDonors, 0, PEXPORT) ;

  for (ii = 0,  jjA = 0 ; up && ii < iiMax && jjMax ; up++, ii++)
    {
      BOOL foundIntron = FALSE ;
      if (up->type != 1 && up->type != 63) continue ;
      if (up->score <= 0 || up->pairScore < 0) continue ;
      mm = up->mm ;
      fragment = up->fragment ;
      x2 = up->x2 ;
      a2 = up->a2 ;
      if ((isDown==1 &&  up->a1 > up->a2) || (isDown == -1 &&  up->a1 < up->a2))
	continue ;
      if (globallyStranded)
	stranded = globallyStranded ;
      else
	stranded = pp->stranded ;
      if (mm->pair < 0) stranded = -stranded ;
      if (debug)
	fprintf (stderr, "ii=%d f=%d s=%d a1=%d a2=%d x1=%d x2=%d\n", ii, mm->fragment, up->score, up->a1, up->a2, up->x1, up->x2) ;

      /* locate on closest corresponding acceptor in x space */
      for (jj = jjA, vp = bigArrp (exportAcceptors, jj, PEXPORT) ;  jj < jjMax && vp->fragment >= fragment && ! foundIntron ; jj++, vp++)
	;
      bestda = 1<<30 ; bestjj = -1 ;
      if (jj == jjMax) { jj-- ; vp-- ; }
      for ( ; jj >= 0 && vp->fragment <= fragment /* && ! foundIntron */ ; jj--, vp--)
	if ((vp->type == 1 || vp->type == 63) && vp->fragment == fragment && 
	    vp->mm == mm && vp->score > 0 && vp->pairScore >= 0)
	  { 
	    if ((isDown==1 &&  vp->a1 > vp->a2) || (isDown == -1 &&  vp->a1 < vp->a2))
	      continue ;
	    if ((isDown==1 &&  vp->a2 < up->a1) || (isDown == -1 &&  vp->a2 > up->a1))
	      continue ;
	    y1 = vp->x1 ; 
	    b1 = vp->a1 ;
	    dx = x2 + 1 - y1 ; /* number of slipage base */
	    if (dx < - shorty)
	      continue ;
	    if (dx > 24)
	      continue ;
	    da = (b1 > a2 ? b1 - a2 : a2 - b1) ;
	    overlap = 0 ;
	    
	    if (up->target != vp->target)
	      continue ;

	    bestda = da ; bestjj = jj ;

	    if (debug) fprintf (stderr, "... jj=%d f=%d s=%d a1=%d a2=%d x1=%d x2=%d\n", jj, vp->mm->fragment, vp->score, vp->a1, vp->a2, vp->x1, vp->x2) ;
	    /* verify that reciprocally there is no better pair */
	    if (0)
	      { /* 2018-11-05 ce truc est idiot si kk < ii on a deja essaye kk */
		y1 = vp->x1 ; 
		b1 = vp->a1 ; 
		for (kk = ii, wp = up ; kk > 0 && wp->fragment >= fragment ; kk--, wp--) ;
		for ( ; kk < iiMax && bestjj >= 0 && wp->fragment <= fragment ; kk++, wp++)
		  {
		    if (wp->score == 0 || wp-> mm != mm || (isDown==1 &&  wp->a1 > wp->a2) || (isDown == -1 &&  wp->a1 < wp->a2))
		      continue ;
		    
		    dx = wp->x2 + 1 - y1 ; /* number of slipage base */
		    if (dx > 24)
		      break ;
		    if (dx < - shorty)
		      continue ;
		    da = (b1 > wp->a2 ? b1 - wp->a2 : wp->a2 - b1) ;
		    if (up->target != vp->target)
		      break ;
		    if (bestda > da) /* ii,jj is not a best reciprocal match */
		      { if (kk != ii) bestjj = -1 ; }
		  }
	      }
	    if (bestjj < 0)
	      continue ;
	    else
	      {
		/* this block documents de uno discoveries and chaining
		 * it has no effect on reporting donors and acceptors or hits geometry
		 */
		int dxMin, dxMax ;
		nnTest++ ;
		jj = bestjj ;
		vp = bigArrp (exportAcceptors, jj, PEXPORT) ;
		
		x2 = up->x2 ;
		a2 = up->a2 ;
		y1 = vp->x1 ; 
		b1 = vp->a1 ; 
		
		da = (b1 > a2 ? b1 - a2 : a2 - b1) ;
		overlap = 0 ;
		dx = x2 + 1 - y1 ; /* number of slipage base */
		
		if (1)
		  { /*  2015_11_02 chaining on very short deletions is a good idea */
		    if (dx < - shorty || dx > 24)
		      continue ;
		    if ((0 && isDown == 1 && a2 - dx > b1 - pp->intronMinLength) || 
			(0 && isDown  == -1 &&  a2 < b1 - dx + pp->intronMinLength) ||
			(isDown  == 1 && a2 < b1 - pp->intronMaxLength) || 
			(isDown  == -1 &&  a2 > b1 + pp->intronMaxLength)
			)
		      continue ;
		  }
		
		if (up->target == vp->target)
		  {
		    if (da < 30)
		      overlap = 12 ; /* small standard intron */ 
		    else if (da< 20000) 
		      overlap = 8 ; /* standard intron */ 
		    else if (da < 100000) 
		      overlap = 12 ; /* big standard intron */ 
		    else if (da < 1000000) 
		      overlap = 16 ; /* big standard intron */ 
		  }
		if (! overlap || isDown == 0)  /* translocations ? */
		  overlap = 20 ;
		
		if (0)
		  fprintf(stderr,  "%d %d   %d %d   ln=%d\n", x2, y1, a2, b1, da) ;
		
		exon1 = dictName (pp->exportDict, up->targetPresuffix) ;
		hook1 = dictName (pp->exportDict, up->presuffix) ;
		
		exon2 = dictName (pp->exportDict, vp->targetPresuffix) ;
		hook2 = dictName (pp->exportDict, vp->presuffix) ;
		
		/* do not report de_uno discoveries if hooks are too short */
		/* check if hooks are long enough */
		if (strlen (hook1) < overlap ||
		    strlen (hook2) < overlap
		    )
		  continue ;
		/* in missed exon case, check if hooks are long enough */
		if (dx < 0 &&
		    ( strlen (hook1) < overlap - dx ||
		      strlen (hook2) < overlap -dx 
		      )
		    )
		  continue ;
		/* check if hook2 hooks into exon1 */
		cr = (char *)hook2 + 2*overlap + (dx < 0 ? -dx : 0) ; cc = *cr ; *cr = 0 ; /* abuse the dictionary privacy */
		cp = strstr (exon1 + 2, hook2 + (dx < 0 ? -dx : 0)) ; *cr = cc ;
		if (!cp)
		  {
		    cr = (char *)hook2 + overlap + (dx < 0 ? -dx : 0) ; cc = *cr ; *cr = 0 ; /* abuse the dictionary privacy */
		    cp = strstr (exon1 + 2, hook2 + (dx < 0 ? -dx : 0)) ; *cr = cc ;
		  }
		if (!cp)
		  continue ;
		dx1 = cp - exon1 ;
		if (dx1 < 2)
		  continue ;
		
		/* check if hook1 hooks into exon2 */
		cr = (char *)hook1 + 2*overlap + (dx < 0 ? -dx : 0) ; cc = *cr ; *cr = 0 ; /* abuse the dictionary privacy */
		cp = strstr (exon2 + dx1, hook1 + (dx < 0 ? -dx : 0)) ; *cr = cc ;
		if (!cp)
		  {
		    cr = (char *)hook1 + overlap + (dx < 0 ? -dx : 0) ; cc = *cr ; *cr = 0 ; /* abuse the dictionary privacy */
		    cp = strstr (exon2 + dx1, hook1 + (dx < 0 ? -dx : 0)) ; *cr = cc ;
		  }
		if (!cp)
		  continue ;
		dx2 = cp - exon2 ;
		if (dx2 < 2)
		  continue ;
		
		if (dx1 != dx2)
		  continue ;
		
		/* verify there is room for the gt-ag feet */
		if (dx1 + dx2  < 4)
		  continue ;
		/* verify that the sliding parts are reverse-complementary */
		for (i = 2 ; i < dx1 ; i++)
		  {
		    ccp = exon1 + i  ; ccq = exon2 + dx1 + 1 - i ;
		    if ( 
			(*ccp == 'a' && *ccq == 't') ||
			(*ccp == 't' && *ccq == 'a') ||
			(*ccp == 'g' && *ccq == 'c') ||
			(*ccp == 'c' && *ccq == 'g')
			 ) ;
		    else
		      { dx1 = dx2 = 0 ; break ; }
		  }
		if (dx1 + dx2 < 4)
		  continue ;
		
		if (dx < 0) /* missing exon, no yet fully analysed */
		  {	
		    if (0 && !isDown)
		      fprintf(stderr,  "....GAP %s\t%s\t%s\t%d\t%d\t%d\t%d\tbestk=%d\tln=%d\n"
			      , isDown ?  "Reverse" : "Forward"
			      , dictName (pp->targetDict, up->target) 
			      , stackText (pp->probeStack, mm->probeName)
			      , up->x2 - bestdx
			      , vp->x1 + dx1 - bestdx - 2
			      , up->a2 + (up->a1 < up->a2 ? - bestdx + 1 : bestdx - 1)
			      , vp->a1 + (up->a1 < up->a2 ?  dx1 - bestdx - 3 : - dx1 + bestdx + 3) 
			      , bestk
			      , da
			      ) ;
		    continue ;
		  }
		/* success, now we adjust the boundary and hope to slide to a gt-ag  */
		
		bestk = -1 ;
		bestdx = 0 ;
		reverse = FALSE ;
		foot[0] = '-' ; foot[1] = 0 ;
		f1[2] = '_' ; f1[5] = 0 ;
		f2[2] = '_' ; f2[5] = 0 ; 
		if (pp->strategy == STRATEGY_RNA_SEQ &&
		    (
		     (isDown == 1 && a2 - dx1 < b1 - pp->intronMinLength) || 
		     (isDown  == -1 &&  a2 > b1 - dx1 + pp->intronMinLength)
		     )
		    )
		  { dxMin = 0 ; dxMax = dx1 - 2 ; }
		else /* always push right */
		  {
		    if (
			(up->a1 < up->a2 && isDown) ||
			(up->a1 > up->a2 && ! isDown) 
			)
		      dxMin = dxMax = dx1 - 2 ;
		    else
		      dxMin = dxMax = 0 ;
		  }	      
		
		/* the valued introns are 
		 * (2019_08_07: added the decoy case, fixed a probable nasty gt_ag strand call error)
		 * s>=0   k1      s<=0   k2 donor_score_acceptor_score
		 * gt_ag 5_4      ct_ac 4_5
		 * gc_ag 4_4      ct_gc 4_4
		 * at_ac 3_3      gt_at 3_3
		 * ct_ac 3_3      gt_ag 3_3
		 * ct_gc 3_2      gc_ag 2_3
		 */
		/* if stranded != 0, we slide but we do not flip */
		for (bestk = -1, i = dxMin ; i <= dxMax ; i++)
		  {
		    char c1, c2 ;
		    k1 = k2 = 0 ; 
		    ccp = exon1 + i  ; ccq = ccp+1 ;
		    f2[3] = *ccp ; f2[4] = *ccq ; 
		    c2 = f1[1] = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)*ccp])] ;
		    c1 = f1[0] = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)*ccq])] ;

		    if (c1 == 'g' && c2 == 't')        /*    d1 == "gt"  */
		      { k1 += 5 ; k2 += 3 ; } 
		    else if (c1 == 'g' && c2 == 'c')   /*    d1 == "gc"  */
		      { k1 += 4 ; k2 += 2 ; } 
		    else if (c1 == 'a' && c2 == 't')   /*    d1 == "at"  */
		      { k1 += 3 ; k2 += 0 ; } 
		    else if (c1 == 'c' && c2 == 't')   /*    d1 == "ct"  */
		      { k1 += 3 ; k2 += 4 ; } 

		    ccp = exon2 + dx1 - i - 2 ; ccq = ccp+1 ;
		    c1 = f1[3] = *ccp ; c2 = f1[4] = *ccq ;

		    f2[1] = complementLetter(*ccp) ;
		    f2[0] = complementLetter(*ccq) ;

		    if (c1 == 'a' && c2 == 'g')        /*    d1 == "ag"  */
		      { k1 += 4 ; k2 += 3 ; } 
		    else if (c1 == 'a' && c2 == 'c')   /*    d1 == "ac"  */
		      { k1 += 3 ; k2 += 5 ; } 
		    else if (c1 == 'g' && c2 == 'c')   /*    d1 == "gc"  */
		      { k1 += 2 ; k2 += 4 ; } 
		    else if (c1 == 'a' && c2 == 't')   /*    d1 == "ac"  */
		      { k1 += 0 ; k2 += 3 ; } 

		    if (stranded < 0) { int kx = k1 ; k1 = k2 ; k2 = kx ; } /* test 2025 */
		    if (globallyStranded) k2 = 0 ;

		    if (k1 > bestk)
		      { bestk = k1 ; bestdx = i ; strcpy (foot, f1) ; reverse = FALSE ; }
		    if ( k2 > bestk)
		      { bestk = k2 ; bestdx = i ; strcpy (foot, f2) ; reverse = TRUE ; }
		    
		    if (0)
		      fprintf (stderr, "XXXXXXXXXX\ta2=%d b1=%d i=%d dxMin=%d dxMax=%d %s stranded=%d %s plus=%d, minus=%d k1=%d  k2=%d\n"
			       ,a2,b1,i,dxMin, dxMax, reverse ? "reverse": "forward", stranded			      
			       , stackText (pp->probeStack, mm->probeName)
			       , pp->nIntronPlus, pp->nIntronMinus
			       , k1, k2
			    ) ;
		  }
		if (bestk < 8 && overlap < 2)
		  continue ;
		
		/* acceptable intron boundaries */
#ifdef JUNK
		{ /* we may want to edit the hits: but this is not inside up, vp but inside the hits table */
		  /* actually, 2015_11_02 we DO NOT want to edit the hits we are happy with double cover of
		   * the repeated region, this is useful in the de_duo analysis and to recognize
		   * sliding introns
		   */
		  up->x2 -= bestdx ; vp->x1 += dx1 - bestdx ; FAUX: up is a donor not  a hit
								if (up->a1 < up->a2)
								  { up->a2 -= bestdx ; vp->a1 += dx1 - bestdx ; }
								else
								  { up->a2 += bestdx ; vp->a1 -= dx1 - bestdx ; }
		}
#endif
		if (debug)
		  fprintf(stderr,  "ZINTRON %s\t%s\t%s\t%d\t%d\t%d\t%d\tbestk=%d\tln=%d\n"
			  , reverse ? "Reverse" : "Forward"
			  , dictName (pp->targetDict, up->target) 
			  , stackText (pp->probeStack, mm->probeName)
			  , up->x2 - bestdx
			  , vp->x1 + dx1 - bestdx - 2
			  , up->a2 + (up->a1 < up->a2 ? - bestdx + 1 : bestdx - 1)
			  , vp->a1 + (up->a1 < up->a2 ?  dx1 - bestdx - 3 : - dx1 + bestdx + 3) 
			  , bestk
			  , da
			  ) ;
		
		foundIntron = TRUE ;
		px = bigArrayp (aa, nIntrons++, PEXPORT) ;
		px->mm = mm ;
		px->intron = nIntrons ;
		px->target = up->target ;
		if (0) {up->donor = 0 ; vp->acceptor = 0 ; }
		px->s1 = up->target ;
		px->s2 = vp->target ;
		px->x1 = up->a1 ;   
		px->a1 = up->a2 + (up->a1 < up->a2 ? - bestdx + 0 : bestdx - 0) ;
		px->x2 = vp->a2 ;
		px->a2 = vp->a1 + (up->a1 < up->a2 ?  dx1 - bestdx - 2 : - dx1 + bestdx + 2) ;
		px->prefix = up->targetPresuffix ;
		px->suffix = vp->targetPresuffix ;
		
		if (pp->strategy != STRATEGY_RNA_SEQ) /* shift right */
		  {
		    if (up->a1 < up->a2)
		      {
			px->a1 += (bestdx - dxMin) ;
			px->a2 += (bestdx - dxMin) ;
		      }
		    else
		      {
			px->a1 += (dxMax - bestdx) ;
			px->a2 += (dxMax - bestdx) ;
		      }
		  }
		
		if (px->mm->bestScore[0]  > 90 || 
		    (px->mm->bestScore[0]  > 60 &&  px->mm->bestScore[0] > .8 * (px->mm->probeLength - px->mm->probeLeftClip - px->mm->probeRightClip))
		    )
		  px->pairScore = mm->mult ;
		px->score = mm->mult ;
		px->fragment = mm->fragment ;
		px->ali = px->a1 < px->a2 ? px->a2 - px->a1 - 1 : px->a1 - px->a2 - 1 ;
		strcpy (px->foot, foot) ;
		px->nN = (px->a1 < px->a2 ? 0x0 : 0x1) ;   /* forward reverse */
		if (!stranded && bestk < 8) px->nN |= 0x2 ;             /* non gt-ag ct-ac */
		if (pp->strategy == STRATEGY_RNA_SEQ &&
		    bestda >=  pp->intronMinLength
		    )
		  {
                    if (bestk >= 9)
		      {
			if (px->a1 < px->a2)
			  mm->gt_ag++ ;
			else
			  mm->ct_ac++ ;
		      }
		    if (bestda > pp->intronMaxLength) px->nN |= 0x4 ;  /* big intron */
		    if (bestda > 200000) px->nN |= 0x8 ;              /* huge intron */
		  }
		else
		  px->nN  |= 0x10 ; /* cannot be called intron */
		
		if ((px->x1 < px->a1 && px->a1 > px->a2) || (px->x1 > px->a1 && px->a1 < px->a2))
		  { px->ali += 2 ; px->nN |=  0x20 ; px->nN ^= 0x1 ; }   /* duplication */
		
		if (px->s1 != px->s2) px->nN |= 0x40 ;       /* translocation between 2 targets */	  
		if ((px->x1 < px->a1 && px->a2 > px->x2) || (px->x1 > px->a1 && px->a2 < px->x2))
		  px->nN |= 0x40 ;                   /* mauvaise topologie */
		
		if ((px->a2 - px->a1) * mm->pair > 0)   px->errLeftVal = mm->mult ; /* my mm is a read > */
		else px->errRightVal = mm->mult ;                /* my mm is a read < */
		if (reverse) 
		  { 
		    int i ;
		    i = px->x1 ; px->x1 = px->x2 ; px->x2 = i ;
		    i = px->a1 ; px->a1 = px->a2 ; px->a2 = i ;
		    i = px->prefix ; px->prefix = px->suffix ; px->suffix = i ;
		    i = up->errLeftVal ; up->errLeftVal = up->errRightVal ; up->errRightVal = i ;
		    px->nN ^= 0x1 ;
		  }
#ifdef JUNK
		too complex to get the dna
		  /* add vcf names */
		  if (px->nN < 0x20 && (px->nN & 0x10))
		    {  /* deletion */
		      char buf[128] ;	
		      int i = 0, dx = px->a2 - px->a1 - 1, dy = 0 ; /* del length */
		      if (dx < strlen(hook1) && dx < 120)
			{
			  const char *cp = hook1 + dx ;
			  while (px->a1 >  dy && cp[0] == cp[dx])
			    { cp-- ; dy++ ; } /* shift left vcf style */
			  memcpy (buf, hook1, dx) ;
			  buf[dx] = 0 ;
			  px->suffix = dictAdd (pp->targetDict,buf, &pp->suffix) ;
			}
		      px->nErr = dy + 1 ; /* shift size */
		    } 
#endif
		/* export the fragment name and intron support */
		if (pp->aoRead2Intron)
		  aceOutf (pp->aoRead2Intron, "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%s\n"
			   , stackText (pp->probeStack, mm->probeName)
			   , mm->bestPairScore
			   , mm->mult
			   , dictName (pp->targetDict, up->target)
			   , px->a1
			   , dictName (pp->targetDict, vp->target)
			   , px->a2
			   , px->foot
			   ) ;
	      }
	  }
    }
  if (debug)
    fprintf (stderr, "// %s  clipAlignConstructIntrons isDown = %d nn=%d introns=%ld\n", timeShowNow(), isDown, nnTest, bigArrayMax (aa)) ;

  return arrayMax (aa) ;
} /* clipAlignConstructIntrons */

/*********************************************************************/
/* merge introns must be done last, because eportHits has ->intron offsets inside exportIntrons */
static int clipAlignMergeIntrons (CLIPALIGN *pp)
{
  PEXPORT *up, *vp ;
  int ii, jj ;
  BigArray aa = pp->exportIntrons ;
  int iMax = arrayMax (aa) ;
  int orientation = 0 ;
  int pass = 0 ;

  for (pass = 0 ; iMax && pass < 2 ; pass += 2)
    {
      bigArraySort (aa, pExportIntronOrder) ;
      bigArrayCompress (aa) ;

      for (ii = 0, up = bigArrp (aa, 0, PEXPORT)  ; ii < iMax ; ii++, up++)
	{
	  if (! up->score) continue ;
	  for (vp = up + 1, jj = ii + 1 ; jj < bigArrayMax (aa) ; jj++, vp++)
	    {
	      if (! vp->score) continue ;
	      if (vp->a1 != up->a1 || vp->a2 != up->a2 || up->s1 != vp->s1 ||  up->s2 != vp->s2)
		break ;
	      if (1)
		{
		  up->score  += vp->score ; up->pairScore += vp->pairScore ;
		  up->errLeftVal += vp->errLeftVal  ; /* supporting fragment orientation */
		  up->errRightVal += vp->errRightVal  ;
		  vp->score = vp->pairScore = 0 ; 

		  if (up->a1 < up->a2 && vp->x1 < up->x1) up->x1 = vp->x1 ;
		  if (up->a1 > up->x2 && vp->x1 > up->x1) up->x1 = vp->x1 ;
		  if (up->a1 < up->a2 && vp->x2 < up->x2) up->x2 = vp->x2 ;
		  if (up->a1 > up->x2 && vp->x2 > up->x2) up->x2 = vp->x2 ;
		}
	    }
	}
      /* keep happy few */
      for (ii = jj = 0, up = vp = bigArrp (aa, 0, PEXPORT)  ; ii < bigArrayMax (aa) ; ii++, up++)
	  if (up->score)
	    {
	      if (ii > jj) *vp = *up ;
	      jj++ ; vp++ ;
	    }
      arrayMax (aa) = jj ;
  
      if (pass) 
	break ;

      /* flip introns */
      if (pp->nIntronPlus + pp->nIntronMinus > 1000)
	{
	  if (pp->nIntronPlus > 90 * pp->nIntronMinus) orientation = 1 ;
	  if (pp->nIntronMinus > 90 * pp->nIntronPlus) orientation = - 1 ;
	}
      for (ii = jj = 0, up = vp = bigArrp (aa, 0, PEXPORT)  ; ii < bigArrayMax (aa) ; ii++, up++)
	{
	  if (
	      ( orientation && /* reverse a NS intron according to dominant stranding of the experiment */
		(up->nN & 0x2) &&  /* non gt-ag, orientation unknown */
		orientation * (up->errLeftVal - up->errRightVal) < 0
		)
	      ||
	      ( ! orientation &&
		!strcmp (up->foot, "ct_ac")
	       ))
	    {
	      char foot[6] ;
	      int i = up->x1 ; up->x1 = up->x2 ; up->x2 = i ;
	      i = up->a1 ; up->a1 = up->a2 ; up->a2 = i ;
	      i = up->prefix ; up->prefix = up->suffix ; up->suffix = i ;
	      i = up->errLeftVal ; up->errLeftVal = up->errRightVal ; up->errRightVal = i ;
	      up->nN ^= 0x1 ; /* flip the orientation flag */
	      strncpy (foot, up->foot, 6) ;
	      up->foot[0] = complementLetter (foot[4]) ;
	      up->foot[1] = complementLetter (foot[3]) ;
	      up->foot[3] = complementLetter (foot[1]) ;
	      up->foot[4] = complementLetter (foot[0]) ;
	      pass = -1 ; /* iterate */
	    }
	}
    }

   return arrayMax (aa) ;
} /* clipAlignMergeIntrons */

/*********************************************************************/

static int clipAlignMergeDoubleIntrons (CLIPALIGN *pp)
{
  PEXPORT *up, *vp ;
  int ii, jj ;
  BigArray aa = pp->exportDoubleIntrons ;
  int iMax = arrayMax (aa) ;

  if (iMax) 
    {
      bigArraySort (aa, pExportDoubleIntronOrder) ;
      bigArrayCompress (aa) ;

      for (ii = 0, up = bigArrp (aa, 0, PEXPORT)  ; ii < iMax ; ii++, up++)
	{
	  if (! up->score) continue ;
	  for (vp = up + 1, jj = ii + 1 ; jj < bigArrayMax (aa) ; jj++, vp++)
	    {
	      if (! vp->score) continue ;
	      if (vp->x1 != up->x1 || vp->x2 != up->x2 ||
		  vp->a2 != up->a2 || vp->s1 != up->s1 ||
		  up->target != vp->target
		  )
		break ;
	      if (1)
		{
		  up->score  += vp->score ;
		  vp->score = 0 ; 
		  if (up->x1 < up->x2 && vp->a1 < up->a1) up->a1 = vp->a1 ;
		  if (up->x1 > up->x2 && vp->a1 > up->a1) up->a1 = vp->a1 ;
		  if (up->x1 < up->x2 && vp->s2 < up->s2) up->s2 = vp->s2 ;
		  if (up->x1 > up->x2 && vp->s2 > up->s2) up->s2 = vp->s2 ;
		}
	    }
	}
      
      /* keep happy few */
      for (ii = jj = 0, up = vp = bigArrp (aa, 0, PEXPORT)  ; ii < bigArrayMax (aa) ; ii++, up++)
	{
	  if (up->score)
	    {
	      if (ii > jj) *vp = *up ;
	      jj++ ; vp++ ;
	    }
	}
      arrayMax (aa) = jj ;
    }

  return arrayMax (aa) ;
} /* clipAlignMergeDoubleIntrons */

/*********************************************************************/
/*********************************************************************/

static void clipAlignExportAll (CLIPALIGN *pp, const char *commandBuf)
{
  int i ;
  BOOL debug = FALSE ;

  fprintf (stderr, "// %s .....clipAlignExportAll\n"
	   , timeShowNow () 
	   ) ;
   
  bigArraySort (pp->exportHits, pExportChainOrder) ;
  bigArrayCompress (pp->exportHits) ;


  if (0 && pp->splice)
    {
      if (! pp->aoRead2Intron)
	pp->aoRead2Intron = aceOutCreate (pp->outFileName, ".frag2intron", pp->gzo, pp->h) ;
      if (pp->strategy == STRATEGY_RNA_SEQ) 
	clipAlignConstructIntrons (pp, 1, FALSE) ;  /* introns on plus strand */
      fprintf (stderr, "// %s after construct introns plus\n", timeShowNow()) ; 
      if (pp->strategy == STRATEGY_RNA_SEQ) 
	clipAlignConstructIntrons (pp, -1, FALSE) ;   /* introns on minus strand */
      fprintf (stderr, "// %s after construct introns minus\n", timeShowNow()) ; 
      if (1) clipAlignConstructIntrons (pp, 0, FALSE) ;  /* rearrangements on remaining donors/acceptors */
      ac_free (pp->aoRead2Intron) ;
      
      fprintf (stderr, "// %s after construct introns 0\n", timeShowNow()) ; 
    }
  
  fprintf (stderr, "// %s .....clipAlignExportAll optimize pair start %ld hits\n"
	   , timeShowNow () 
	   , bigArrayMax (pp->exportHits)
	   ) ;
  clipAlignOptimizePairs (pp, FALSE) ; /* global clean up on best pair score */
  fprintf (stderr, "// %s .....clipAlignExportAll optimize pair done\n"
	   , timeShowNow () 
	   ) ;
  
  if (debug) showPexport (pp, 0) ;
  
  bigArraySort (pp->exportHits, pExportHitOrder) ;
  
  if (debug) showPexport (pp, 0) ;
  
  /* export the hits */
  clipAlignExportPx (pp, 0) ;
  if (pp->sam)
    clipAlignExportSam (pp, commandBuf) ;
  if (!pp->silent) aceOutf (pp->ao, "# %s done\n", timeShowNow ()) ;
  ac_free (pp->ao) ;
  
  
  /* export the list of probes for which each seed is seen over 25 times */
  if (pp->seedOver25 && bitSetCount (pp->seedOver25))
    {
      MM *mm ;
      pp->ao = aceOutCreate (pp->outFileName, ".overRepresentedSeeds", pp->gzo, pp->h) ;
      for (i = 0, mm = arrp (pp->probes, i, MM) ; i < arrayMax (pp->probes) ; mm++, i++)
	if (bit (pp->seedOver25, i))
	  aceOutf (pp->ao, "%s\n", stackText (pp->probeStack, mm->probeName)) ;
      ac_free (pp->ao) ;
    }

  if (pp->splice)   /* flip introns if non stranded run, merge introns, and blank out the donor acceptors */
    {
      clipAlignMergeIntrons (pp) ;
      clipAlignMergeDoubleIntrons (pp) ;
    }

  /* export the donors and acceptors */
  bigArraySort (pp->exportDonors, pExportDonorOrder) ;
  bigArraySort (pp->exportAcceptors, pExportAcceptorOrder) ;
  clipAlignSquish (pp, 4) ;
  clipAlignSquish (pp, 5) ;
  if (0) clipAlignRestoreAcceptorScores (pp) ;  /* restore the score of the donors/acceptors */

  if (bigArrayMax (pp->exportDonors) + bigArrayMax (pp->exportAcceptors))
    {
      pp->ao = aceOutCreate (pp->outFileName, ".overhangs", pp->gzo, pp->h) ;
      
      aceOutf (pp->ao, "#Tag\tScore\tMultiplicity\tPosition\tfrom\tType\tTarget\tStrand\tPosition\tExon (anti if donor)\tTarget overhang\tRead overhang (anti if acceptor)\tNumber of supporting reads\n") ;
      
      clipAlignExportPx (pp, 1) ;
      clipAlignExportPx (pp, 2) ;   
      ac_free (pp->ao) ;
    }

  fprintf (stderr, "// %s before construct introns\n", timeShowNow()) ;
  /* export the introns */
  ac_free (pp->exportHits) ; /* make room */
  
  if (pp->splice)
    {
      fprintf (stderr, "// %s after merge introns\n", timeShowNow()) ;

      if (bigArrayMax (pp->exportIntrons))
	{
	  pp->ao = aceOutCreate (pp->outFileName, ".introns", pp->gzo, pp->h) ;
	  aceOutf (pp->ao, "#Type\tStrand\tFrom target\ta1\ta2\tTo target\tb1\tb2\tBoundary\tLength\tNumber of supporting reads\tNumber of supporting quasi perfect reads\n") ;
	  clipAlignExportPx (pp, 3) ;
	  ac_free (pp->ao) ;
	}
      if (bigArrayMax (pp->exportDoubleIntrons))
	{
	  pp->ao = aceOutCreate (pp->outFileName, ".doubleintrons", pp->gzo, pp->h) ;
	  aceOutf (pp->ao, "#Type\tStrand\tFrom target\ta1\ta2\tx1\tx2\ts1\ts2\tNumber of supporting reads\n") ;
	  clipAlignExportPx (pp, 4) ;
	  ac_free (pp->ao) ;
	}
    }
  if (! pp->sam)
    { /* samStats, useful to compare to the output of sam2gold */
      AC_HANDLE h = ac_new_handle () ;
      const char *METHOD = pp->method ? pp->method : "02_ClipAlign" ;
      ACEOUT ao = aceOutCreate (pp->outFileName, ".s2g.samStats", 0, h) ;
      const char *run = pp->runName ? pp->runName : "-" ;
      
      aceOutf (ao, "%s\t%s\tnRawReads\t%ld\n", run, METHOD, pp->nRawReads) ;
      aceOutf (ao, "%s\t%s\tnRawBases\t%ld\n", run, METHOD, pp->nRawBases) ;
      aceOutf (ao, "%s\t%s\tnReads\t%ld\n", run, METHOD, pp->nReads) ;  
      
      if (0)   aceOutf (ao, "%s\t%s\tnAlignedPairs\t%ld\n", run, METHOD, pp->nAlignedPairs) ;
      aceOutf (ao, "%s\t%s\tnAlignedReads\t%ld\t%.2f%%\n", run, METHOD, pp->nAlignedReads, (100.0 * pp->nAlignedReads)/(pp->nRawReads + .000001)) ;
      
      aceOutf (ao, "%s\t%s\tnAlignedBases\t%ld\t%.2f%%\n", run, METHOD, pp->nAlignedBases, (100.0 * pp->nAlignedBases) / (pp->nRawBases + .000001)) ;
      aceOutf (ao, "%s\t%s\tnMultiAligned %d times\t%ld\t%.2f%%\n", run, METHOD, 0
	       , pp->nUnaligned
	       , 100.0 * pp->nUnaligned / (pp->nRawReads + .000001)
	       ) ;
      for (int j = 1 ; j < 11 ; j++)
	if (pp->nMultiAligned[j])
	  {
	    aceOutf (ao, "%s\t%s\tnMultiAligned %d times\t%ld\t%.2f%%\n", run, METHOD
		     , j, pp->nMultiAligned[j]
		     , 100.0 * pp->nMultiAligned[j] / (pp->nRawReads + .000001)
		     ) ;
	  }
      long int verif = 0 ;
      for (int j = 1 ; j < 11 ; j++)
	verif += pp->nMultiAligned[j] ;
      aceOutf (ao, "nReads = %ld , sum of multiAli = %ld, verif = %ld\n", pp->nRawReads, verif, pp->nRawReads - verif - pp->nUnaligned) ; 
       
      aceOutf (ao, "%s\t%s\tnErrors\t%ld\t%.4f%%\n", run, METHOD, pp->nErrors, (100.0 * pp->nErrors)/(pp->nAlignedBases + 0.00001)) ;
      ac_free(h) ;
    }
  ac_free (pp->ao) ;
} /* clipAlignExportAll */

/********************************************************************/

static void clipAlignTableCaption (CLIPALIGN *pp)
{
  if (!pp->silent) 
    {
      aceOutf (pp->ao, "# %s Magic aligner\n", timeShowNow()) ;
      aceOutf (pp->ao, "# Author: Danielle et Jean Thierry-Mieg, NCBI, mieg@ncbi.nlm.nih.gov\n") ;
      if (0)
	aceOutf (pp->ao, "# Latest version and documentation are  available from https://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/Software\n") ;
      aceOutf (pp->ao, "# Sequences are 1-based\n") ;
      if(0)
	{
	  aceOutf (pp->ao, "# 5\'All DNA sequences are read on the strand going away from the edges of the alignment\n") ;
	  aceOutf (pp->ao, "# So for donors and acceptors the match part (identical on read and target) and the overhang part are read on the opposite strands,\n"
		   "# The advantage of this unusual convention is that the beginning of all these words is independent of the length of the read, \n"
		   "# independent of the strand of the read or of the target, and the overhang at the donor site is identical to the \n"
		   "# match at the acceptor site and vice versa, so the donor acceptor pairs are easy to recognize.\n"
		   ) ;
	}
      aceOutf (pp->ao, "#Read\tAlignment score\tRead multiplicity") ;
      aceOutf (pp->ao, "\tRead length to be aligned, i.e. removing eventual adaptors, leading Ts or trailing As not matching the target") ;
      aceOutf (pp->ao, "\tRead length aligned") ;
      aceOutf (pp->ao, "\tfrom base") ;
      aceOutf (pp->ao, "\tto base, coordinates of the alignment in the read") ;
      
      
      aceOutf (pp->ao, "\tTarget_class (A_mito ... Z_genome) defines a Bayesian hierarchy of most desirable targets, used by some of the post-processing steps. Theses names are prioritized alphabetically") ;
      aceOutf (pp->ao, "\tGene, as given in the target fasta file as : >target_name|Gene|gene_name, to avoid confusions, RefSeq gene names are prefixed by X__") ;
      aceOutf (pp->ao, "\tUniqueness, number of genes or genomic sites on which the read aligns at best score, if the strand of the best alignements varies, the number is negative, hence -2 likely indicates a single alignment on 2 genes antisense to each other") ;
      aceOutf (pp->ao, "\tTarget, the identifier given in the -t fasta file, for example a chromosome or a trancscritpt") ;
      
      aceOutf (pp->ao, "\tfrom base") ;
      aceOutf (pp->ao, "\tto base, coordinates of the alignment in the target") ;
      
      if (pp->solid)
	aceOutf (pp->ao, "\tNumber of SOLiD corrected bases") ;
      else
	aceOutf (pp->ao, "\tNumber of N bases") ;
      aceOutf (pp->ao, "\tNumber of mismatches") ;
      
      
      
      aceOutf (pp->ao, "\tMissmatches: read position:type,... single base variation  (base in target > base in read); single base insertion (+ base inserted in read); single base deletion (- target base missing from read), given in the coordinates and the orientation of the read") ;
      aceOutf (pp->ao, "\tMissmatches: target position:type,... single base variation  (base in target > base in read); single base insertion (+ base inserted in read); single base deletion (- target base missing from read), given in the coordinates and the orientation of the target") ;
      
      aceOutf (pp->ao, "\t5' overhang read backwards, i.e complementary to the unaligned 5' part of the read, max %d bp", OVLN) ;
      aceOutf (pp->ao, "\t3' overhang read forward, i.e. identical to the unaligned 3' part of the read, max %d bp", OVLN) ;
      
      aceOutf (pp->ao, "\tTarget sequence immediately upstream of the alignment") ;
      aceOutf (pp->ao, "\tTarget sequence immediately downstream of the alignment") ;
      aceOutf (pp->ao, "\tFragment length as measured on target in case of paired end sequencing or 0 if not available at this stage of the program") ;
      aceOutf (pp->ao, "\n") ;
      
      aceOutf (pp->ao, "#Read\tAlignment score\tRead multiplicity\tLength to align\tAligned length\tfrom base\tto base in read") ;
      aceOutf (pp->ao, "\tTarget class\tGene (X__name for RefSeq)\tTarget multiplicity\tTarget name (_name for RefSeq)\tfrom base\tto base on target") ;
      
      aceOutf (pp->ao, "\tNumber of %s bases\tNumber of mismatches\tPosition:type in read\tPosition:type in target\t5'overhang (reverse-complemented)\t3'Overhang\tTarget prefix\tTarget suffix\tPair length", pp->solid ? "corrected" : "N") ;
      
      aceOutf (pp->ao, "\n") ;
    }
} /* clipAlignTableCaption */

/********************************************************************/

void clipAlignExportHit (CLIPALIGN *pp, MM *mm, int xpos
			 , int a1, int a2, int x1, int x2
			 , int nErr, int nN
			 , char *errLeftVal, char *errRightVal
			 , unsigned char *prefix, unsigned char *suffix, unsigned char *targetPrefix
			 )
{
  static int nLine = 0 ;
  int intronBonus = 0 ;
  int pairTampon = pp->pairTampon ;
  BOOL debug = FALSE ;
  
  if (*prefix == 63) intronBonus += INTRONBONUS ;
  if (*suffix == 63) intronBonus += INTRONBONUS ;
  
  if (debug) fprintf (stderr, "---------- proposing %s  %d %d    %s %d %d\tscore %d\n"
		      , stackText (pp->probeStack, mm->probeName)
		      , x1, x2
		      , dictName (pp->targetDict, pp->target)
		      , a1, a2
		      , pp->score) ;
  
  if (pp->antiprobe)
    { 
      int aa, ln = mm->probeLength ;
      aa = a1 ; a1 = a2 ; a2 = aa ; 
      x1 = ln - x1 + 1 ; x2 = ln - x2 + 1 ;
      aa = x1 ; x1 = x2 ; x2 = aa ; 
    }
  if (mm->probeLeftClip)
    {
      x1 += mm->probeLeftClip ;
      x2 += mm->probeLeftClip ;
    }
  if (pp->aceOutGenomic)
    {
      const char *cosmidName = dictName (pp->targetDict, pp->target) ;
      
      if (!strncmp(cosmidName, "GM_",3))
	cosmidName += 3 ;
      if (!strncasecmp(cosmidName, "mrna:",5))
	cosmidName += 5 ;
      aceOutf (pp->ao, "Probe %s\n"
	       , freeprotect (stackText (pp->probeStack, mm->probeName))
	       ) ;
      
      if (!nErr)
	aceOutf (pp->ao, "Genome_exact_hit %s %d %d %d\n"
		 , freeprotect (cosmidName)
		 , a1, a2, nErr
		 ) ;
      aceOutf (pp->ao, "Genome_approximate_hit %s %d %d %d\n\n"
	       , freeprotect (cosmidName)
	       , a1, a2, nErr
	       ) ;
      if (mm->probeSetName)
	aceOutf (pp->ao, "Probe %s\nProbe_set %s\n\n"
		 , freeprotect (stackText (pp->probeStack, mm->probeName))
		 , freeprotect (stackText (pp->probeStack, mm->probeSetName))
		 ) ;
    }
  else if (pp->aceOutMrna)
    {
      const char *cosmidName = dictName (pp->targetDict, pp->target) ;
      
      if (!strncmp(cosmidName, "GM_",3))
	cosmidName += 3 ;
      if (!strncasecmp(cosmidName, "mrna:",5))
	cosmidName += 5 ;
      aceOutf (pp->ao, "Probe %s\n"
	       , freeprotect (stackText (pp->probeStack, mm->probeName))
	       ) ;
      if (!nErr && a1 < a2)
	aceOutf (pp->ao, "mRNA_exact_hit %s %d %d %d\n"
		 , freeprotect (cosmidName)
		 , a1, a2, nErr
		 ) ;
      if (a1 < a2)
	aceOutf (pp->ao, "mRNA_hit %s %d %d %d\n\n"
		 , freeprotect (cosmidName)
		 , a1, a2, nErr
		 ) ;
      else
	aceOutf (pp->ao, "mRNA_anti_hit %s %d %d %d\n\n"
		 , freeprotect (cosmidName)
		 , a1, a2, nErr
		 ) ;	
      
      aceOutf (pp->ao, "mRNA %s\n" , freeprotect (cosmidName)) ;
      if (a1 < a2)
	aceOutf (pp->ao, "Probe_hit %s %d  %d\n\n"
		 , freeprotect (stackText (pp->probeStack, mm->probeName))
		 , a1, a2, nErr
		 ) ;
      else
	aceOutf (pp->ao, "Probe_anti_hit %s %d  %d %d\n\n"
		 , freeprotect (stackText (pp->probeStack, mm->probeName))
		 , a1, a2, nErr
		 ) ;
      if (mm->probeSetName)
	aceOutf (pp->ao, "Probe %s\nProbe_set %s\n\n"
		 , freeprotect (stackText (pp->probeStack, mm->probeName))
		 , freeprotect (stackText (pp->probeStack, mm->probeSetName))
		 ) ;
    }
  else if (pp->aceOutNM)
    {
      const char *cosmidName = dictName (pp->targetDict, pp->target) ;
      
      aceOutf (pp->ao, "Probe %s\n"
	       , freeprotect (stackText (pp->probeStack, mm->probeName))
	       ) ;
      if (!nErr && a1 < a2)
	aceOutf (pp->ao, "NM_exact_hit %s %d %d %d\n"
		 , freeprotect (cosmidName)
		 , a1, a2, nErr
		 ) ;
      
      if (a1 < a2)
	aceOutf (pp->ao, "NM_hit %s %d %d %d\n\n"
		 , freeprotect (cosmidName)
		 , a1, a2, nErr
		 ) ;
      else
	aceOutf (pp->ao, "NM_anti_hit %s %d %d %d\n\n"
		 , freeprotect (cosmidName)
		 , a1, a2, nErr
		 ) ;
      aceOutf (pp->ao, "Sequence %s\n" , freeprotect (cosmidName)) ;
      if (a1 < a2)
	aceOutf (pp->ao, "Probe_NM_hit %s %d %d %d\n\n"
		 , freeprotect (stackText (pp->probeStack, mm->probeName))
		 , a1, a2, nErr
		 ) ;
      else
	aceOutf (pp->ao, "Probe_NM_anti_hit %s %d %d %d\n\n"
		 , freeprotect (stackText (pp->probeStack, mm->probeName))
		 , a1, a2, nErr
		 ) ;
      if (mm->probeSetName)
	aceOutf (pp->ao, "Probe %s\nProbe_set %s\n\n"
		 , freeprotect (stackText (pp->probeStack, mm->probeName))
		 , freeprotect (stackText (pp->probeStack, mm->probeSetName))
		 ) ;
    }
  else 
    {
      char *cp ;
      
      if (! nLine++)
	clipAlignTableCaption (pp) ;
      
      if (1)
	{
	  BOOL ok = TRUE ;
	  PEXPORT *px = 0 ;
	  int dx, ddx, dda, ipxg, score ;
	  int ipxgMax = bigArrayMax (pp->exportGeneHits) ;

 	  score = pp->score - intronBonus ;
	  dx = x2 - x1 ;
	  
	  for (ipxg = mm->ipxg ; ok && ipxg > 0 && ipxg < ipxgMax ; ipxg = px ? px->ipxold : 0)
	    { 
	      px = bigArrayp (pp->exportGeneHits, ipxg, PEXPORT) ;
	      if (! px->score)
		continue ;
	      if (pp->target == px->target && x1 == px->x1 && x2 == px->x2 && a1 == px->a1 && a2 == px->a2)
		{ ok = FALSE ; continue ; }
	      if ( ! pp->hasPairs &&   /* we need to wait until the pair is mapped */      
		   pp->MRNAH && score == px->score && pp->targetGene == px->targetGene && pp->target > px->target)
		{ 
		  ok = FALSE ; 
		  continue ;
		}
	      if (0 && score < px->score - pairTampon)
		{
		  if (pp->bestHit)
		    { ok = FALSE ; continue ; }
		  
		  ddx = (px->x2 < x2 ? px->x2 : x2) - (px->x1 > x1 ? px->x1 : x1) ;
		  if (2 * ddx < dx)
		    continue ;
		  if ( px->nErr <= nErr)
		    { ok = FALSE ; continue ; }
		  if (px->nErr == nErr  && ( score < 40 || ddx > dx))
		    { ok = FALSE ; continue ; }
		  if (px->nErr == nErr  && 2 * ddx > dx)
		    { 
		      if (pp->target == px->target)
			{
			  int aa1 = a1 < a2 ? a1 : a2 ;
			  int aa2 = a1 < a2 ? a2 : a1 ;
			  int ba1 = px->a1 < px->a2 ? px->a1 : px->a2 ;
			  int ba2 = px->a1 < px->a2 ? px->a2 : px->a1 ;
			  
			  dda = (aa2 < ba2 ? aa2 : ba2) - (aa1 < ba1 ? ba1 : aa1) ;
			  if (2 * dda > aa2 - aa1)
			    { ok = FALSE ; continue ; }
			}
		    }
		}
	      
	    }
	  if (ok)
	    {
	      if (debug) fprintf (stderr, "----- accepting %s  %d %d    %s %d %d\tscore %d\tmax=%ld\n"
				  , stackText (pp->probeStack, mm->probeName)
				  , x1, x2
				  , dictName (pp->targetDict, pp->target)
				  , a1, a2
				  , score
				  , bigArrayMax (pp->exportGeneHits)
				  ) ;

	      if (1)
		{  /* sanitize because px is not reset to zero */
		  if (mm->ipxg && mm->ipxg >= bigArrayMax (pp->exportGeneHits) - 1) 
		    mm->ipxg = 0 ;
		  if (mm->ipxg)
		    {
		      px = bigArrp (pp->exportGeneHits, mm->ipxg, PEXPORT) ;
		      if (px->mm != mm)
			mm->ipxg = 0 ;
		    }
		}
				
	      px = bigArrayp (pp->exportGeneHits, bigArrayMax (pp->exportGeneHits), PEXPORT) ;
	      /* ATTENTION:  px was not zeroed */
	      memset (px, 0, sizeof (PEXPORT)) ;
	      
	      pp->nExportedHits++ ;
	      px->ipxold = mm->ipxg ;
	      mm->ipxg = bigArrayMax (pp->exportGeneHits) - 1 ;
	      px->mm = mm ; px->score = score ; px->fragment = mm->fragment ;
	      px->ali = x2 - x1 + 1 + pp->aliBonusForThisAli ;
	      
	      px->ln = 
		( pp->rightAdaptorClipForThisAli ? pp->rightAdaptorClipForThisAli : mm->probeLength) 
		+ ( pp->leftAdaptorClipForThisAli ? pp->leftAdaptorClipForThisAli : mm->probeLeftClip) 
		;
	      px->a1 = a1 ; px->a2 = a2 ; 
	      px->x1 = x1 ; px->x2 = x2 ;
	      px->nN = nN ; px->nErr = nErr ;
	      px->target = pp->target ;
	      px->targetGene = pp->targetGene ;
	      cp = (prefix[1] != 0 ? (char*)prefix + 1 : "-") ; dictAdd (pp->exportDict, cp, &(px->prefix)) ;
	      cp = (suffix[1] != 0 ? (char*)suffix + 1 : "-") ; dictAdd (pp->exportDict, cp, &(px->suffix)) ;
	      px->presuffix = ((prefix[0] & 0xff) << 8) | ((suffix[0] & 0xff) ) ;
	      dictAdd (pp->exportDict, (char *)targetPrefix, &(px->targetPresuffix)) ;
	      dictAdd (pp->exportDict, errLeftVal, &(px->errLeftVal)) ;
	      dictAdd (pp->exportDict, errRightVal, &(px->errRightVal)) ;
	    }
	}  
    }
} /* clipAlignExportHit  */

/*********************************************************************/
/*********************************************************************/
typedef struct llStruct LL ;
struct llStruct { int na, nt, ng, nc, p ; LL *next ; }  ;
static int lettersOrder (const void *va, const void *vb)
{
  const LL *a = (const LL *)va, *b = (const LL *)vb ;
  int n ;
  n = a->na - b->na ; if (n) return n ;
  n = a->nt - b->nt ; if (n) return n ;
  n = a->ng - b->ng ; if (n) return n ;
  n = a->nc - b->nc ; if (n) return n ;
  n = a->p - b->p ; 
  return n ;
} /*lettersOrder */

/*********************************************************************/

static void clipAlignCountLetters (CLIPALIGN *pp)
{
  LL *ll, *ll0 ;
  MM *mm ;
  int na, nt, nc, ng, jj, nMax = arrayMax (pp->probes) ;
  const unsigned char *ccp, *probeSeq ;

  pp->letters = arrayHandleCreate (nMax, LL, pp->h) ;
  for (jj = 0 ; jj < arrayMax(pp->probes) ; jj++)
    {
      mm = arrp (pp->probes, jj, MM) ;
      na = nt = nc = ng = 0 ;
      probeSeq = (unsigned char*) stackText (pp->probeStack, mm->probeSequence) ;
      for (ccp = probeSeq ; *ccp ; ccp++)
	switch (*ccp)
	  {
	  case A_: na++ ; break ;
	  case T_: nt++ ; break ;
	  case G_: ng++ ; break ;
	  case C_: nc++ ; break ;
	  }
      ll = arrayp (pp->letters, jj, LL) ;
      ll->p = jj ; ll->na = na ; ll->nt = nt ; ll->ng = ng ;ll->nc = nc ; 
    }
  /* sort by na/nt/ng/nc */
  arraySort (pp->letters, lettersOrder) ;
  /* establish short cuts to next numbers */
  for (jj = 0, ll = ll0 = arrayp (pp->letters, 0, LL) ; jj < arrayMax (pp->letters) ; ll++, jj++)
    {
      if (ll->na != ll0->na ||
	  ll->nt != ll0->nt || 
	  ll->ng != ll0->ng || 
	  ll->nc != ll0->nc
	  )
	{ ll0->next = ll ; ll0 = ll ; }
    }  
} /* clipAlignCountLetters */

/*************************************************************************************/

static void clipAlignInitB2 (CLIPALIGN *pp, BOOL isUp)
{
  int i ; 
 
 
  i = 256 ;  while (i--) B2[i] = 0 ;
  B2[A_] = 0x0 ; B2[T_] = 0x3 ;     /* ATTENTION copied in dnasubs.c:dnagetWordUsage */
  B2[G_] = 0x1 ; B2[C_] = 0x2 ;     /* you must keep the 2 identical */

  if (0) showPexport(0, 0) ; /* for compiler happiness */

  if (pp->A2G) { B2[R_] = B2[A_] ; B2[Y_] = B2[T_] ; }
  if (pp->G2A) { B2[R_] = B2[G_] ; B2[Y_] = B2[C_] ; }
  if (pp->T2C) { B2[Y_] = B2[T_] ; B2[R_] = B2[A_] ; }
  if (pp->C2T) { B2[Y_] = B2[C_] ; B2[R_] = B2[G_] ; }
  if (pp->A2C) { B2[M_] = B2[A_] ; B2[K_] = B2[T_] ; }
  if (pp->C2A) { B2[M_] = B2[C_] ; B2[K_] = B2[G_] ; }
  if (pp->T2G) { B2[K_] = B2[T_] ; B2[M_] = B2[A_] ; }
  if (pp->G2T) { B2[K_] = B2[G_] ; B2[M_] = B2[C_] ; }
  if (pp->G2C) { B2[S_] = isUp ? B2[C_] : B2[G_] ; }
  if (pp->C2G) { B2[S_] = isUp ? B2[G_] : B2[C_] ; }
  if (pp->A2T) { B2[W_] = isUp ? B2[T_] : B2[A_] ; }
  if (pp->T2A) { B2[W_] = isUp ? B2[A_] : B2[T_] ; }

  return ;
} /* clipAlignInitB2 */

/*********************************************************************/
static int clipAlignGetProbeHits (CLIPALIGN *pp, Array dna, Array dnaR, BOOL isUp)
{
  int i, ii, j, jx, jxB, pos, max, dx2, nErr, nN, nHits = 0, same ;
  int a1, a2, xpos, x1, x2 ;
  MM *mm, *mm0 ;
  char *cp=0, *cpB, errLeftVal[10024], errRightVal[10024] ;
  unsigned char prefix[64+OVLN], targetPrefix[64 + 4 * OVLN], suffix[64+OVLN] ;
  int nf = 0 ;
  static int NNN = 0 ;
  int iBucket ;
  int nTargets = pp->nTargets ;
  BOOL oligoGeneTouched = FALSE ;
  BOOL doubleSeed = pp->doubleSeed ;
  BitSet oligoGene = pp->nTargets > 0 ? pp->oligoGene : 0 ;  /* useless on first target */
  unsigned long int oligo = 0,  oligoB = 0, myOligo, anchor, mask1 = 0, mask, maskB ;
  const void *vp;
  int hashPhase = pp->hashPhase ;
  int nOligo ;
  Associator myass ;
  Associator ass = pp->ass ;
  Associator ass1 = pp->ass1 ;
  Associator ass2 = pp->ass2 ;
  Associator assB = pp->assB ;
  Associator assB1 = pp->assB1 ;
  Array oligoFrequency = pp->oligoFrequency ;
  Array oligoFrequencyB = pp->oligoFrequencyB ;
  const char *cosmidName = dictName (pp->targetDict, pp->target) ;
  int seedLength = pp->seedLength, nok ;
  Array bucket = 0 ; 

  COSMID=cosmidName ;

  clipAlignInitB2 (pp, isUp) ;
  if (0)  fprintf (stderr, "get probe hits cosmid=%s isUp=%d  %s\n", cosmidName, isUp, timeShowNow());
  dx2 = pp->dx2 ;
  mask = pp->mask ;
  mask1 = pp->mask1 ;
  maskB = is64 ? 0xffffffffffffffff  : 0xffffffff ;
  if (mask1)
    messcrash ("clipAlignGetProbeHits received obsolete option mask1 != 0") ;
  ii = pp->probeLengthMin - 1 ;
  max = arrayMax (dna) ;
  if (max < 100 && max < pp->probeLengthMin) return 0 ;

  /* load the mask */
  jx = dx2 ;
  jxB = is64 ? 32 : 16 ;
  oligo = 0 ;  oligoB = 0 ; nok = seedLength ;
  cp = arrp(dna, 0, char) - 1 ;
  if (dx2)  
    {
      for (i = 0 ; i < jx - 1 ; i++)
	{
	  cp++ ;
	  oligo <<= 2 ; oligo |= B2[(int)(*cp)] ;    
	  if (*cp == N_) nok = seedLength ;
	  else nok-- ;
	}
      oligo &= mask ;
    }
  if (jxB)
    {
      cpB = arrp(dna, 0, char) - 1 ;
      for (i = 0 ; i < jxB - 1 ; i++)
	{
	  cpB++ ;
	  oligoB <<= 2 ; oligoB |= B2[(int)(*cpB)] ; 
	}
      oligoB &= maskB ;
    }
  j =  pp->probeLengthMin ;
  if (pp->minAli && j > pp->minAli) j = pp->minAli ;
  pos = -1 ; ii = arrayMax(dna) - dx2 ;
  if (ii > 0) while (pos++, cp++, cpB++, ii--)
		{ 
		  if (0 && pp->hashPhase && pos == 4910877  && pp->target == 17 && isUp == FALSE)
		    invokeDebugger() ;
		  if (cp > arrp (dna, arrayMax(dna)-1, char))
		    {
		      char *cq0 = arrp (dna, 0, char) ;
		      char *cq1 = arrp (dna, arrayMax(dna) -1, char) ;
		      char *zero = 0 ;

		      messcrash ("bad pointer too far by %ld bytes cp=%ld, cq0=%ld cq1=%ld, ii=%d maxDna=%d dx2=%d pos=%d hashPhase=%d", cp - cq1, cp - zero, cq0 - zero, cq1 - zero, ii, arrayMax(dna), dx2, pos, hashPhase) ;
		    }
		  oligo = ((oligo << 2) | B2[(int)*cp]) & mask ;
		  if (*cp == N_) nok = seedLength ;
		  else nok-- ;  
		  if (cpB >= arrp (dna, arrayMax(dna)-1, char))
		    cpB-- ;
		  oligoB = ((oligoB << 2) | B2[(int)*cpB]) & maskB ;
		  if (! oligo || nok > 0)
		    continue ; 
      
		  if (hashPhase == 0) 
		    {
		      if (pp->isShortTarget) /* fill ass1 on the fly */
			{ 
			  if (1) 
			    { 
			      if (! assFind (ass1, assVoid(oligo), &vp))
				{
				  pp->nOligo++ ;
				  assInsert (ass1, assVoid(oligo), assVoid(pp->nOligo)) ;
				  nOligo = pp->nOligo ;
				}
			      else
				nOligo = assInt(vp) ;
			      nf = array (oligoFrequency, nOligo, unsigned char) ;
			      if (oligoGene)
				{
				  if (nf < 255 && ! bitt(oligoGene, oligo))
				    {
				      bitSet (oligoGene, oligo) ; oligoGeneTouched = TRUE ;
				      array (oligoFrequency, nOligo, unsigned char) = nf + 1 ;
				    }
				}
			      else
				{
				  if (nf < 5 || (nTargets < 5 && nf < 255)) 
				    array (oligoFrequency, nOligo, unsigned char) = nf + 1 ;
				}
			    }
			  if (doubleSeed) /* big oligos */
			    {
			      if (! assFind (assB1, assVoid(oligoB), &vp))
				{
				  pp->nOligo++ ;
				  assInsert (assB1, assVoid(oligoB), assVoid(pp->nOligo)) ;
				  nOligo = pp->nOligo ;
				}
			      else
				nOligo = assInt(vp) ;
			      nf = array (oligoFrequencyB, nOligo, unsigned char) ;
			      if (nf < 5 || (nTargets < 5 && nf < 255)) 
				array (oligoFrequencyB, nOligo, unsigned char) = nf + 1 ;
			    }
			}
		      else
			{ 
			  if (1)
			    {
			      /*
				if (pos < 20)
				fprintf (stderr, ".... testing pos %d\n", pos) ;
			      */ 
			      if (assFind(ass1, assVoid(oligo), &vp)) 
				{
				  /*
				    if (pos < 20)
				    fprintf (stderr, ".... testing pos %d  success\n", pos) ;
				  */
				  nOligo = assInt(vp) ;
				  nf = array (oligoFrequency, nOligo, unsigned char) ;
				  if (oligoGene)
				    {
				      if (nf < 255 && ! bitt(oligoGene, oligo))
					{
					  bitSet (oligoGene, oligo) ; oligoGeneTouched = TRUE ;
					  array (oligoFrequency, nOligo, unsigned char) = nf + 1 ;
					}
				    }
				  else
				    {
				      if (nf < 5 || (nTargets < 5 && nf < 255)) 
					array (oligoFrequency, nOligo, unsigned char) = nf + 1 ;
				    }
				  if (doubleSeed) /* The short oligo is included in the long one */
				    {
				      if (assFind(assB1, assVoid(oligoB), &vp)) 
					{
					  nOligo = assInt(vp) ;
					  nf = array (oligoFrequencyB, nOligo, unsigned char) ;
					  if (nf < 5 || (nTargets < 5 && nf < 255)) 
					    array (oligoFrequencyB, nOligo, unsigned char) = nf + 1 ;
					}
				    }
				}
			    }
			}
		      continue ;
		    }
		  myass = 0 ;
		  /*
		    if (assFind (ass2, assVoid(oligo), 0) && assFind(ass, assVoid(oligo), &vp))
		    { myass = ass ; myOligo = oligo ; }
		    else if (doubleSeed && assFind (assB2, assVoid(oligoB), 0) && assFind(assB, assVoid(oligoB), &vp))
		    { myass = assB ; myOligo = oligoB ; }
		  */
		  if (doubleSeed)
		    { myass = assB ; myOligo = oligoB ; }
		  else
		    { myass = ass ; myOligo = oligo ; }
		  iBucket = 0 ; bucket = 0 ;
		  if (assFind (ass2, assVoid(oligo), 0))   /* (myass) */
		    while (assFindNext(myass, assVoid(myOligo), &vp, &bucket, &iBucket))
		      {
			int anchorMax = arrayMax (pp->probes) ;
			NNN++ ;
			anchor = assULong(vp)-1 ;
			if (is64)
			  {
			    xpos = (anchor >> 32) & 0xefffffff ; /* avoid the sign bit ! */
			    anchor = anchor & 0xffffffff ; 
			  }
			else
			  {
			    xpos = (anchor >> 22) & 0xefffffff ; /* avoid the sign bit ! */
			    anchor = anchor & 0x3fffff ; 
			  }
			mm0 = arrp (pp->probes, anchor, MM) ;
			for (mm = mm0, same = mm->probeLength ;  mm && anchor < anchorMax ; anchor++, mm = (same  >= xpos  + seedLength - 1 ? mm+1 : 0) )
			  {
			    if (mm->same < same) same = mm->same ; /* part in future mm common with mm0 */
			    if (pp->maxHit && mm->nGeneHits * 30 > pp->maxHit * mm->probeLength) /* assume 30bp per exon */
			      continue ; 
			    if (mm->probeLength < xpos  + seedLength - 1)
			      continue ;
			    if (mm->lastTarget == pp->target)  
			      {   /* same hit as before  ? */
				if (xpos >= mm->lastX1  && xpos <= mm->lastX2 && 
				    pos >= mm->lastA1 && pos <= mm->lastA2  &&
				    pos - xpos >= mm->lastDelta - 20 &&
				    pos - xpos <= mm->lastDelta + 20 &&
				    mm->lastIsUp == isUp
				    )
				  {
				    mm->lastDelta = pos - xpos ;
				    continue ; 
				  }
			      }
			    if (pp->MRNAH && 
				!pp->hasPairs &&  /* 2015_02_02  is hasPairs we cannot do that because if r1 goes just in b and run2 goes in a,b we need to keep b */
				mm->lastTarget &&
				mm->lastTargetGene == pp->targetGene &&
				mm->lastIsUp == isUp &&
				mm->lastTarget < pp->target &&     /* assuming the targets fasta files are sorted in alphabetic order */
				mm->nPerfectHit > 0
				)  
			      continue ;
		
			    /* do not retry the very same word too often throughout the whole run
			     * since we cannot count every position, we count modulo a hash
			     * we do not zero this number when we get a better hit
			     */
			    if (mm->probeLength < 128) /* catastrophic on long reads because we fold back on 8 counters
							* we have 8 * seedStep = 40 so it only work for seq length up to 128
							*/
			      {
				int limit = 64 ; /* must be smaller than 255 = max char value */
				unsigned int xk = xpos ;
				int w ;
				xk /= pp->seedShift ;
				xk = (xk) ^ (xk * 13) ^ (xk * 37) ; xk = xk ^ (xk > 3) ^ (xk >> 5) ;
				xk = xk % 8 ;
				xk <<= 1 ; if (isUp) xk |= 0x1 ;
				w = mm->wordCount[xk] ;
				if (w >= limit)
				  continue ;
				mm->wordCount[xk]++ ;
			      }
			    nErr = nN = 0 ;
			    NNN++ ;
			    if (0 && xpos < 230 && xpos > 130 && pp->target == 17 && isUp == FALSE)
			      invokeDebugger() ;
			    if (1 && !strcmp(stackText (pp->probeStack, mm->probeName), "cYOK-024J17<"))
			      invokeDebugger() ;
			    if (clipAlignVerifyProbeHit (pp, mm, dna, dnaR, pos, xpos
							 , (unsigned char*) stackText (pp->probeStack, mm->probeSequence)
							 , mm->probeLength, isUp
							 , &nErr, &nN
							 , &a1, &a2, &x1, &x2
							 , errLeftVal, errRightVal
							 , prefix, suffix, targetPrefix
							 )
				)
			      {
				nHits++ ;
				if (isUp)
				  {
				    a1 = max - a1 + 1 ; a2 = max - a2 + 1 ; 
				  }
				clipAlignExportHit (pp, mm, xpos, a1, a2, x1, x2, nErr, nN
						    , errLeftVal, errRightVal
						    , prefix, suffix, targetPrefix
						    ) ;
		    
			      }
			  }
		      }
		}
  pp->oligoGeneTouched = oligoGeneTouched ;
  return nHits ;
} /* clipAlignGetProbeHits */

/*************************************************************************************/

static int clipAlignDoSearch (CLIPALIGN *pp, Array dna0)
{
  int nHits = 0 ;
  Array dna = dna0, dnaR = 0, dnaS = 0 ;
  AC_HANDLE h = ac_new_handle () ;

  /* the prefix and walls are read from dna0 since dna may be color coded */
  pp->dna0 = dna0 ;

  if (pp->targetMask)
    {
      int ii, i, iMax = arrayMax (dna0) ;
      register char *cp ;
      PEXPORT *up ;

      for (ii = 0, up = arrayp (pp->targetMask, 0, PEXPORT) ; ii < arrayMax (pp->targetMask) ; ii++, up++)
	if (up->target == pp->target)
	  {
	    for (i = up->a1 - 1, cp = arrp (dna0, i, char) ; i < up->a2 ; cp++, i++)
	      if (i >= 0 && i < iMax) *cp = N_ ; 
	  }
    }

  if (1) /* align first the antistrand  */
    {    
      if (pp->antistrand || ! pp->strand)
	{
	  dnaR = dnaHandleCopy (dna0, h) ;
	  reverseComplement (dnaR) ;
	  if (pp->solid)
	    {
	      pp->dna0 = dnaHandleCopy (dnaR, h) ;
	      dnaSolidEncodeArray (dnaR, TRUE) ;
	      dnaS = dnaHandleCopy (dna0, h) ;
	      dnaSolidEncodeArray (dnaS, TRUE) ;
	    }
	  else
	    { pp->dna0 = dnaR ; dnaS = dna ; }
	  arrayLock (dnaR) ; arrayLock (dnaS) ;
	  nHits += clipAlignGetProbeHits (pp, dnaR, dnaS, TRUE) ;
	  arrayUnlock (dnaR) ; arrayUnlock (dnaS) ;
	}
      if (pp->strand || ! pp->antistrand)
	{
	  if (pp->solid)
	    {
	      if (!dnaS)
		{
		  dnaS = dnaHandleCopy (dna0, h) ;
		  dnaSolidEncodeArray (dnaS, TRUE) ;
		}
	    }
	  else
	    dnaS = dna0 ;
	  pp->dna0 = dna0 ;
	  arrayLock (dnaS) ;
	  nHits += clipAlignGetProbeHits (pp, dnaS, 0, FALSE) ; 
	  arrayUnlock (dnaS) ;
	}
    }
  pp->dna0 = 0 ;
  ac_free (h) ;
  return nHits ;
} /* clipAlignDoSearch */

/*************************************************************************************/
/* select best transcript target for each pair, keep a single one if MRNAH 
 * simply reset score = 0 to discard unwanted its
 */
static void clipAlignOptimizeCleanUp (CLIPALIGN *pp, BigArray aa, int ii1, int ii2)
{
  int ii ;
  PEXPORT *px ;

  for (ii = ii1, px =bigArrp (aa, ii1, PEXPORT) ; ii < ii2 ; px++, ii++)
    if (px->mm)
      {
	px->mm->pairScore = 0 ;
	if (px->mm->pair > 1)
	  {
	    MM *mm1 = arrp (pp->probes, px->mm->pair - 2, MM) ;
	    mm1->pairScore = 0 ;
	  }
	if (px->mm->pair < -1)
	  {
	    MM *mm1 = arrp (pp->probes, - px->mm->pair - 2, MM) ;
	    mm1->pairScore = 0 ;
	  } 
      }
} /* clipAlignOptimizeCleanUp */
/*************************************************************************************/
/* aa has the final list of candidate exons, remove the extraneous donors/acceptors/introns */
static void  clipAlignOptimizeIntronCleanUp (CLIPALIGN *pp, BigArray aa, int ii1, int ii2, BOOL singleTarget,
					     int *t0p, int *w0p, int *d0p, int *c0p)
{  
  int fragment = 0, target = 0 ;
  int ii, jj ;
  int a2, b1, b2 ;
  int nValidatedIntrons = 0 ;
  BigArray introns = singleTarget ? pp->geneIntrons : pp->exportIntrons ;
  BigArray donors =  singleTarget ? pp->geneDonors : pp->exportDonors ;
  BigArray acceptors =  singleTarget ? pp->geneAcceptors : pp->exportAcceptors ;
  int tMax = introns ? bigArrayMax (introns) : 0 ;
  PEXPORT *px, *py, *pt ;
  int t, t1 = tMax + 2, t2 = 0 ;

  return ;
  if (! singleTarget) /* reserve negative score is pairscore to use them when searching rearrangmnets */
    for (ii = ii1, px = bigArrp (aa, ii1, PEXPORT) ; ii < ii2 ; px++, ii++)
      {
	px->pairScore = px->score ;
      }

  if (1)
    for (ii = ii1, px = bigArrp (aa, ii1, PEXPORT) ; ii < ii2 ; px++, ii++)
      {
	if (px->intron)
	  {
	    py = bigArrp (introns, px->intron - 1, PEXPORT) ;
	    py->score = 0 ;
	  }
	if (0 && px->donor)
	  {
	    py = bigArrp (donors, px->donor - 1, PEXPORT) ;
	    py->score = 0 ;
	  }
 	if (0 && px->acceptor)
	  {
	    py = bigArrp (acceptors, px->acceptor - 1, PEXPORT) ;
	    py->score = 0 ;
	  }
	if (px->intron < t1)
	  t1 = px->intron ;
	if (px->intron > t2)
	  t2 = px->intron ;
      }
  for (ii = ii1, px = bigArrp (aa, ii1, PEXPORT) ; ii < ii2 ; px++, ii++)
    {
      if (px->fragment)
	fragment = px->fragment ;
    }

  if (! fragment)
    {
      if (singleTarget)
	{
	  if (donors)    bigArrayMax (donors) = 0 ;
	  if (acceptors) bigArrayMax (acceptors) = 0 ;
	  if (introns)   bigArrayMax (introns) = 0 ;
	}
      return ;
    }
  
  /* scan the accepted exons and reestablish a positive score for the introns/donors/acceptors of those exons */
  for (ii = ii1, px = bigArrp (aa, ii, PEXPORT) ; ii < ii2 ; px++, ii++)
    {
      BOOL intronFound = FALSE ;
      if (! px->score || ! px->mm)
	continue ;   
      target = px->target ;
      a2 = px->a2 ;
      if (1 && introns && bigArrayMax (introns) && t1 < t2)
	for (t = t1, pt = bigArrp (introns, t, PEXPORT) ; t < t2 ; t++, pt++)
	  {
	    int targ = 0 ;
	    int du = px->a1 < px->a2 ? 10 : 2 ;
	    int dv = px->a1 < px->a2 ? 2 : 10 ;
	    
	    if (pt->s1 == target || pt->s2 == target)
	      {
		b1 = pt->a1 ; b2 = pt->a2 ; targ = pt->s1 ; }
	    if (b1 > a2 - du && b1 < a2 + dv && targ == target)
	      { /* found an intron look for next exon */
		int e1 ;
		    
		for (jj = ii + 1, py = px + 1 ; jj < ii2 ; py++, jj++)
		  {
		    e1 = py->a1 ; targ = py->target ; 
		    if (e1 > b2 - du && e1 < b2 + dv &&
			targ == target &&
			py->score && py->mm
			)
		      {
			int mult =  py->mm ?  py->mm->mult : 0 ;
			pt->score += px->mm->mult ; /* success */
			nValidatedIntrons++ ;
			/* count the orientation of the reads supporting gt_ag to influence at export time
			 * the orientation of the NS introns
			 */
			if (! strcmp (pt->foot, "gt_ag"))
			  {
			    if (pt->errLeftVal > pt->errRightVal)
			      pp->nIntronPlus += mult ;
			    else
			      {
				pp->nIntronMinus  += mult ;
			      }
			  }
			else if (! strcmp (pt->foot, "ct_ac"))
			  {
			    if (pt->errLeftVal > pt->errRightVal)
			      {
				pp->nIntronMinus += mult ;
			      }
			    else
			      pp->nIntronPlus  += mult ;
			  }
			py->acceptor = 0 ;
			intronFound = TRUE ;
			break ;  /* break the py loop */
		      }
		  }
		break ; /* break the pt loop */
	      }
	  }
      if (! intronFound && px->donor)
	{	
	  py = bigArrp (donors, px->donor - 1, PEXPORT) ;
	  py->score = px->mm->mult ;
	}
      if (px->acceptor)
	{	
	  py = bigArrp (acceptors, px->acceptor - 1, PEXPORT) ;
	  py->score = px->mm->mult ;
	}
    }
} /* clipAlignOptimizeIntronCleanUp */

/*************************************************************************************/
/* set px-> to the cumulated score of the exon intron chain going through each exon */
static int clipAlignOptimizeHalfChain (CLIPALIGN *pp, BigArray aa, int ii1, int ii2, int chain)
{
  int ii, jj, t, t1 = 0 ;
  int a2, b1, b2, db, c1 ;
  int minChainAli = pp->minChainAli ;
  BigArray introns = pp->geneIntrons ;
  int tMax ;
  PEXPORT *px, *py, *pt ; 
  BOOL debug = FALSE ;

  if (0 && ! introns && pp->exportIntrons)
    introns = pp->exportIntrons ;
  tMax = introns ? bigArrayMax (introns) : 0 ;
  
  /* stitch overlapping exons */
  for (ii = ii1, px = bigArrp (aa, ii, PEXPORT) ; ii < ii2 ; px++, ii++)
    {
      if (! px->score)
	continue ;
      for (jj = ii + 1, py = px + 1 ; jj < ii2 ; py++, jj++)
	if (px->score && py->score && px->mm == py->mm && px->target == py->target)
	  {
	    PEXPORT *pa, *pb ;
	    const char *cpa, *cpb, *cta, *ctb ;
	    char *cp, *ct, bufp[12], buft[12] ;
	    if (px->x1 < py->x1)
	      { pa = px ; pb = py ;}
	    else
	      { pa = py ; pb = px ;}
	    
	    if (pa->x2 < py->x1 || 
		! pa->errLeftVal ||
		! pa->errRightVal ||
		! pb->errLeftVal ||
		! pb->errRightVal 
		)
	      continue ;
	    /* look for matching error */
	    cpa = dictName (pp->exportDict, pa->errRightVal) ;
	    cta = dictName (pp->exportDict, pa->errLeftVal) ;
	    cpb = dictName (pp->exportDict, pb->errRightVal) ;
	    ctb = dictName (pp->exportDict, pb->errLeftVal) ;
	    
	    strncpy (bufp, cpb, 10) ; bufp[11] = 0 ;
	    cp = strchr (bufp, ',') ; if (cp) *++cp = 0 ; else continue ;
	    cp = strstr (cpa, bufp) ;
	    if (!cp)
	      continue ;
	    
	    strncpy (buft, ctb, 10) ; buft[11] = 0 ;
	    ct = strchr (buft, ',') ; if (ct) *++ct = 0 ; else continue ;
	    ct = strstr (cta, buft) ;
	    if (!ct)
	      continue ;
	    else if (0) /* stitch */
	      {
		char buf[20000] ;
		int n ;
		/* we found an exact match, stitch pb behind pa */
		pa->x2 = pb->x2 ; pa->a2 = pb->a2 ; 
		pb->score = pb->pairScore = 0 ; /* kill pb */
		n = cp - cpa ;
		memcpy (buf, cpa, n) ;
		memcpy (buf + n, cpb, strlen(cpb)) ;
		buf[n + strlen(cpb)] = 0 ;
		dictAdd (pp->exportDict, buf, &(pa->errRightVal)) ;

		n = 1 ; cp = buf ; while (*++cp) if (*cp == ',') n++ ;
		pa->nErr = n ;
		pa->ali = pa->x2 - pa->x1 + 1 ;
		pa->score = pa->ali - n * pp->errCost ;
		
		n = ct - cta ;
		memcpy (buf, cta, n) ;
		memcpy (buf + n, ctb, strlen(ctb)) ;
		buf[n + strlen(ctb)] = 0 ;
		dictAdd (pp->exportDict, buf, &(pa->errLeftVal)) ;
	      }
	  }
	  
      
    }

  /* construct an nN chain for all exons */
  if (0)
    for (ii = ii1, px = bigArrp (aa, ii1, PEXPORT) ; ii < ii2 ; px++, ii++)
      px->chain = px->pairScore = 0 ;
  for (ii = ii1, px = bigArrp (aa, ii, PEXPORT) ; ii < ii2 ; px++, ii++)
    { /* tmp space to store the exon chaining and previous score */
      px->s1 = px->s2 = 0 ;
    }
  for (ii = ii1, px = bigArrp (aa, ii, PEXPORT) ; ii < ii2 ; px++, ii++)
    {
      MM *mm = px->mm ;	 
      int isDown = (px->a1 < px->a2 ? 1 : -1) ;
      if (! px->score || ! mm)
	continue ;   
      if (! px->chain)
	{
	  MM *mm = px->mm ;
	  int d0 = px->donor ;
	  int a0 = px->acceptor ;
	  int jj, iMax = bigArrayMax (aa) ;
	  for (jj = ii + 1,  py = px + 1 ; jj <  iMax && py->mm == mm && py->x1 == px->x1 ; jj++, py++)
	    {
	      if (! py->score)
		continue ;
	      else if (py->a1 == px->a1)
		{
		  int d = py->donor ;
		  int a = py->acceptor ;

		  py->donor = d0 ; 
		  py->acceptor = a0 ;
		  if (! memcmp(px, py, sizeof (PEXPORT)))
		    py->score = 0 ;
		  else
		    {
		      py->donor = d ; 
		      py->acceptor = a ;
		    }
		}
	    }

	  px->chain = ++chain ; /* new chain */
	}
 
      /* look for a quasi continuous alignment */
      if (1)
	{
	  int x2 = px->x2 ; 
	  int a2 = px->a2 ; 
	  int target = px->target ;

	  for (jj = ii + 1, py = px + 1 ; jj < ii2 ; py++, jj++)
	    if (py->score && py->mm == mm && py->target == target &&
		py->pairScore <= px->score + px->pairScore)
	      {
		int dx = py->x1 - x2 ;
		int da = isDown * (py->a1 - a2) ;
		if (isDown * (py->a2 - py->a1) > 0 && dx - da > -3 && dx - da < 3 && dx >= -3 && dx < 20)
		  {
		    int ddx = px->x2 - py->x1 + 1 ;
		    if (ddx < -2)
		      ddx = 0 ;
		    py->s1 = ii + 1 ;
		    py->s2 = px->score + px->s2 - ddx ;;
		    break ;  /* break the py loop */
		  }
	      }
	}

      if (! pp->splice)
	continue ;
      /* reclip 8bp at each end and try to find an intron corresponding to the exon */
      a2 = px->a2 ; 

      if (t1 < tMax && introns && bigArrayMax (introns))
	for (t = t1, pt = introns ? bigArrp (introns, t, PEXPORT) : 0 ; t < tMax ; t++, pt++)
	  {
	    b1 = pt->a1 ; b2 = pt->a2 ; 
	    db = b2 - b1 ;

	    if (isDown * db < 0)
	      { int b0 = b1 ; b1 = b2 ; b2 = b0 ; db = -db ; }

	    if (
		( isDown == 1 && b1 < a2 + 2 && b1 > a2 - 20) || 
		( isDown == -1 && b1 > a2 - 2 && b1 < a2 + 20) 
		)
	      { /* found an intron look for next exon */
		for (jj = ii + 1, py = px + 1 ; jj < ii2 ; py++, jj++)
		  {
		    int dx = py->x1 - px->x2 ;
		    int da ;
		    

		    c1 = py->a1 ; da = py->a1 - px->a2 ; 

		    if (isDown * (py->a2 - py->a1) > 0 &&  isDown * (db - da) + dx == 1 &&
			py->pairScore <= px->score + px->pairScore)
		      {	
			if (! py->score || isDown * (py->a2 - py->a1) < 0)
			  continue ;
			
			if (
			    (
			     ( isDown == -1 && c1 > b2 - 2 && c1 < b2 + 20) ||
			     ( isDown == 1 && c1 < b2 + 2 && c1 > b2 - 20) 
			     ) &&
			    py->score && py->mm
			    )
			  {
			    int ddx ;
			    
			    ddx = px->x2 - py->x1 + 1 ;
			    /* do not hook on exon with better predecessor */
			    if (ddx <= -2)
			      ddx = 0 ;
			    if (py->s2 <= px->score + px->s2 - ddx)
			      {
				py->s2 = px->score + px->s2 - ddx;
				py->s1 = ii + 1 ;
				py->intron = t + 1 ;
				break ;  /* break the py loop */
			      }
			  }
		      }
		  }
		/* break ; do not break the pt loop */
	      }
	  }
    }

  /* find recursivelly the best exon chain and set the chain */
  while (1)
    {
      int pairScore, chain, ln = 0, bestS2 = 0 ; ;
      t1 = 0 ;
      for (ii = ii1, px = bigArrp (aa, ii1, PEXPORT) ; ii < ii2 ; px++, ii++)
	if (px->s1)
	  {
	    py = bigArrp (aa, px->s1 - 1, PEXPORT) ;
	    if (py->pairScore)
	      px->s1 = px->s2 = 0 ;
	  }
	
      for (ii = ii1, px = bigArrp (aa, ii1, PEXPORT) ; ii < ii2 ; px++, ii++)    
	if (px->score > 0 && px->pairScore == 0 && px->score + px->s2 > bestS2)
	  {
	    bestS2 = px->score + px->s2 ;
	    t1 = 1 ;
	    t = ii ; 
	  }
      if (! t1)
	break ;
      /* t is the best exon, set its pairScore and those of its predecessor */
      px = bigArrp (aa, t, PEXPORT) ;
      pairScore = px->score ; /* + px->s2 ; */
      /* px->score = px->pairScore = pairScore ; */
      chain = px->chain ;
      ln =  px->x2 - px->x1 + 1 ;
      while (px->s1)
	{
	  int ddx ; 
	  py = px ;
	  px = bigArrp (aa, px->s1 - 1, PEXPORT) ;
	  if (px->pairScore) /* already hook in a strongre chain */
	    break ;
	  px->chain = chain ;
	  px->suffix = py->prefix = -1 ;  /* rm the target presuffix at export time */

	  ddx = px->x2 - py->x1 + 1 ;
	  /* do not hook on exon with better predecessor */
	  if (ddx <= -2)
	    ddx = 0 ;
	  ln += px->x2 - px->x1 + 1 - ddx ;
	  pairScore += px->score - ddx ;
	}
      /* set the global ali */
      px = bigArrp (aa, t, PEXPORT) ;
      if (! minChainAli || ln >=  minChainAli)
	{
	  px->ali = ln ;
	  px->score = px->pairScore = pairScore ;
	  while (px->s1)
	    {
	      px = bigArrp (aa, px->s1 - 1, PEXPORT) ;
	      if (px->pairScore) /* already hook in a stronger chain */
		break ;
	      px->ali = ln ;
	      px->score = px->pairScore = pairScore ;
	    }
	}
      else /* kill */
	{
	  px->score = px->pairScore = px->chain = 0 ;
	  while (px->s1)
	    {
	      px = bigArrp (aa, px->s1 - 1, PEXPORT) ;
	      if (px->pairScore) /* already hook in a strongre chain */
		break ;
	      px->score = px->pairScore = px->chain = px->ali = 0 ;
	    }
	}
 
    }



  /* add a bonus if a chain touches a chain of the other read of the same fragment */

  if (debug) showPexport (pp, 1) ;

  for (ii = ii1, px = bigArrp (aa, ii, PEXPORT) ; ii < ii2 ; px++, ii++)
    { /* restore tmp space to store the exon chaining and previous score */
      px->s1 = px->s2 = 0 ;
    }

  return chain ;
}  /* clipAlignOptimizeHalfChain */

/*************************************************************************************/
/* Tie together with raphia topologically compatible chains
 * starting from thr longest ones
 */

typedef struct raphiaStruct { int nn, a1, a2, x1, x2, isDown, chain, nr, previous, score, rScore, ali ; } RF ;

static void showRaphia (Array aa)
{
  int i, iMax = aa ? arrayMax (aa) : 0 ;
  RF *rf ;

  for (i = 0 ; i < iMax ; i++)
    {
      rf = arrp (aa, i, RF) ;
      fprintf (stderr, "Raphia %3d nn=%d nr=%d <-%d score=%d rScore=%d ali=%d chain=%d x1/x2 %d/%d a1/a2 %d/%d\n"
	       , i
	       , rf->nn, rf->nr, rf->previous
	       , rf->score, rf->rScore, rf->ali, rf->chain
	       , rf->x1, rf->x2, rf->a1, rf->a2
	       ) ;
    }
}

static int raphiaXOrder (const void *va, const void *vb)
{
  const RF *a = (const RF *)va, *b = (const RF *)vb ;
  int n ;
  n = a->isDown - b->isDown ; if (n) return n ;
  n = a->x1 - b->x1 ; if (n) return n ;
  n = a->x2 - b->x2 ; if (n) return -n ;  /* long dx first */
  n = a->nn - b->nn ; if (n) return -n ;  /* long first */

  return n ;
} /*raphiasXOrder */

static int raphiaChainOrder (const void *va, const void *vb)
{
  const RF *a = (const RF *)va, *b = (const RF *)vb ;
  int n ;
  n = a->chain - b->chain ; if (n) return n ;
  n = a->isDown - b->isDown ; if (n) return n ;
  n = a->x1 - b->x1 ; if (n) return n ;
  n = a->x2 - b->x2 ; if (n) return -n ;  /* long dx first */
  n = a->nn - b->nn ; if (n) return -n ;  /* long first */

  return n ;
} /* raphiaChainOrder */

static int raphiaNrOrder (const void *va, const void *vb)
{
  const RF *a = (const RF *)va, *b = (const RF *)vb ;
  int n ;
  n = a->nr - b->nr ; if (n) return n ;
  n = a->x1 - b->x1 ; if (n) return n ;
  n = a->x2 - b->x2 ; if (n) return -n ;  /* long dx first */
  n = a->nn - b->nn ; if (n) return -n ;  /* long first */

  return n ;
} /* raphiaNrOrder */

/*************************************************************************************/
static int clipAlignOptimizeHalfRaphia (CLIPALIGN *pp, BigArray aa, int ii1, int ii2)
{
  AC_HANDLE h = ac_new_handle () ;
  Array raphia = arrayHandleCreate (128, RF, h) ;
  KEYSET scores, chainScores, chain2nr, nr2chain ;
  RF *rf, *rg ;
  int nRF = 0 ;
  int ii, jj, iMax, nr ;
  PEXPORT *px ;
  int isDown = 1 ;
  int daMax = 1000000 ;
  int malus = 2 * pp->errCost ;
  int bestScore, bestii ;
  /* select the orientation */
  /* count the chains and the number of exons in each chain */
  for (ii = ii1, px = bigArrp (aa, ii, PEXPORT) ; ii < ii2 ; px++, ii++)
    {
      if (! px->score)
	continue ;
      rf = arrayp (raphia, px->chain, RF) ;
      rf->chain = px->chain ;
      if (!rf->nn)
	nRF++ ;
      rf->nn++ ;
      rf->score = px->score ;
      rf->ali = px->ali ;
      if (! rf->x1 || px->x1 < rf->x1)
	rf->x1 = px->x1 ;
      if (px->x2 > rf->x2)
	rf->x2 = px->x2 ;
      isDown = px->isDown ? 1 : -1 ;
      if (! rf->a1)
	{ 
	  rf->isDown = isDown ;
	  rf->a1 = px->a1 ;
	  rf->a2 = px->a2 ;
	}
      if (isDown == 1)
	{
	  if (rf->a1 > px->a1)
	    rf->a1 = px->a1 ;
	  if (rf->a2 < px->a2)
	    rf->a2 = px->a2 ;
	}
      else
	{
	  if (rf->a1 < px->a1)
	    rf->a1 = px->a1 ;
	  if (rf->a2 > px->a2)
	    rf->a2 = px->a2 ;
	}
    }
  /*
    mapped_read= 
    6S    20=                 2X17=             1I25=15S
    GTTCACAACGCTTTCATTGCAGGAAATAATGCTATATGGCCATTTGCTGACTTGTCGTAACAAAGTTTCCAATACCATTAGTGATAG   ln= 87
    SSSSSSAACGCTTTCATTGCAGGAAAATATGCTATATGGCCATTGICTGACTTGTCGTAACAAAGTTTCCASSSSSSSSSSSSSSS    ln= 86

    mapped_read= 
    GTCGGCCAACCATTCTCTGAATGTCTGAATGACAGCGACTATTGTTTAAGGTTTT    GTAAATAAGAGACCTTAAGCAATGCATTGACTGCATGGTCTCAACCTCATCTGTGCTGTGAGAA   ln= 119
    SSSGGCCAACCATTCTCTGAATGTCTGAATGACAGCGACTATIGTTTAAGGTTTTTTTTGTAAATAAGAGACCTTAAGCAATGCATTGACTGCATGGTCTCAACCTCATCTGTGCTGTGAGAA    ln= 123

    ttctcacagcacagatgaggttgagaccatgcagtcaatgcattgcttaaggtctcttatttacaaaaA    ccttaaacatagtcgctgtcattcagacattcagagaatggttggcc
    ttctcacagcacagatgaggttgagaccatgcagtcaatgcattgcttaaggtctcttatttacaaaaaaaaaccttaaacatagtcgctgtcattcagacattcagagaatggttggcc

  */
  if (nRF < 2)
    goto done ;
  arraySort (raphia, raphiaChainOrder) ;
  for (ii = jj = 0, rf = rg = arrp (raphia, 0, RF) ; ii < arrayMax (raphia) ; ii++, rf++)
    {
      if (! rf->nn)
	continue ;
      if (jj < ii)
	*rg = *rf ;
      jj++ ; rg++ ;
    }
  arrayMax (raphia) = jj ;
  if (jj < 2)
    goto done ;

  /* merge raphia from same chain */
  iMax = arrayMax (raphia) ;
  for (ii = jj = 0, rf = rg = arrp (raphia, 0, RF) ; ii < iMax ; ii++, rf++)
    {
      if (! rf->nn)
	continue ;
      if (jj < ii && rg->isDown == rf->isDown && rg->chain == rf->chain)
	{
	  if (rg->x1 > rf->x1)
	    rg->x1 = rf->x1 ;
	  if (rg->x2 < rf->x2)
	    rg->x2 = rf->x2 ;
	  if (rf->isDown * (rg->a1 - rf->a1) > 0)
	    rg->a1 = rf->a1 ;
	  if (rf->isDown * (rg->a2 - rf->a2) < 0)
	    rg->a2 = rf->a2 ;
	  rg->nn += rf->nn ;
	  continue ;
	}
      if (jj < ii)
	*rg = *rf ;
      jj++ ; rg++ ;
    }
  iMax = arrayMax (raphia) = jj ;

  if (jj < 2)
    goto done ;

  arraySort (raphia, raphiaXOrder) ;

  for (nr = 1 ; nr <= iMax ; nr++)
    { /* iterate untill there is no rf not assigned to a raphia */
      BOOL found = FALSE ;
      
      for (ii = 0 ; ii < iMax ; ii++)
	{ /* reset */
	  rf = arrp (raphia, ii, RF) ;
	  if (! rf->nn)
	    continue ; 
	  if (rf->nr)
	    continue ;
	  found = TRUE ;
	  rf->previous = 0 ;
	  rf->rScore = 0 ;
	}
      if (! found)
	break ;

      /* find the heaviest path */
      bestScore = bestii = 0 ;
      for (ii = 0 ; ii < iMax ; ii++)
	{
	  int a1 = 0, a2 = 0, y2 = 0 ;
	  int rScore ;
	  
	  rf = arrp (raphia, ii, RF) ;
	  if (! rf->nn)
	    continue ; 
	  if (rf->nr)
	    continue ;
	  
	  rScore = rf->rScore + rf->score ;
	  if (rScore > bestScore)
	    {
	      bestScore = rScore ;
	      bestii = ii ;
	    }
	  for (jj = ii + 1, rg = rf + 1 ; jj <  iMax ; jj++, rg++)
	    {
	      if (! rg->nn)
		continue ;
	      if (rg->nr)
		continue ;

	      if (rg->isDown != rf->isDown)
		break ;
	      
	      /* stop if the y1-y2 exon would fit */
	      if (y2 && rg->x1 > y2 - 4)
		break ;
	      
	      /* set a target wall at the repeats of rf */
	      if (rg->x1 < rf->x1 + 10 && rg->x2 > rf->x2 - 10)
		{
		  if (rg->isDown * (rg->a1 - rf->a1) > 0 && 
		      (!a2 || isDown *( a2  - rg->a1) > 0)
		      )
		    a2 = rg->a1 ;
		  if (rg->isDown * (rg->a1 - rf->a1) < 0 && 
		      (!a1 || isDown *( a1  - rg->a2) < 0)
		      )
		    a1 = rg->a2 ;
		}
	      
	      /* reject if overlap */
	      if (rg->x1 < rf->x2 - 8)
		continue ;
	      if (rg->isDown * (rg->a1 - rf->a2) > 1.1 * (rg->x1 - rf->x2) &&
		  rg->isDown * (rg->a1 - rf->a2) > rg->x1 - rf->x2 + 8)
		continue ;
	      
	      if (rg->rScore > rScore)
		break ;
	      
	      /* reject rg behind walls */
	      if (a2 && (rg->isDown * (rg->a1 - a2) > 0))
		continue ;
	      if (a1 && (rg->isDown * (rg->a2 - a1) < 0))
		continue ;
	      
	      /* reject too distant rg */
	      if (rg->isDown * (rg->a1 - rf->a2) > daMax)
		continue ;
	      
	      /* check the topology */
	      if (rg->isDown * (rg->a1 - rf->a1) * (rg->x1 - rf->x1) < 0)
		continue ;
	      if (rg->isDown * (rf->a2 - rg->a1) - (rf->x2 - rg->x1) > 10 &&
		  rf->x2 > rg->x1)
		continue ;
	      y2 = rg->x2 ;
	      rg->rScore = rScore ; rg->previous = ii + 1 ;
	      if (rg->score + rScore > bestScore)
		{
		  bestScore = rg->score + rScore ;
		  bestii = jj ;
		}
	    }
	}
      /* start from bestii, and assign the nr backwards */
      jj = bestii + 1 ;
      while (jj)
	{
	  rg = arrp (raphia, jj - 1, RF) ;
	  rg->nr = nr ;
	  jj = rg->previous ;
	}
    }
 
  /* compute the score of all the raphias */
  scores = keySetHandleCreate (h) ;
  chainScores = keySetHandleCreate (h) ;
  chain2nr = keySetHandleCreate (h) ;
  nr2chain = keySetHandleCreate (h) ;
 
  arraySort (raphia, raphiaNrOrder) ;
  for (ii = 0, rf = arrp (raphia, ii, RF); ii < iMax ; ii++, rf++)
    {
      if (! rf->nn || ! rf->nr)
	continue ;
      if (keySet (scores, rf->nr))
	{
	  int dx = rf->x1 - rf[-1].x2 - 1 ;
	  if (dx < 0)
	    {
	      dx = -dx ;
	    }
	  else 
	    dx = 0 ;
	  keySet (scores, rf->nr) += rf->score - dx - malus ;
	}
      else
	{
	  keySet (scores, rf->nr) = rf->score ;
	}
    }
  for (ii = 0, rf = arrp (raphia, ii, RF); ii < iMax ; ii++, rf++)
    {
      if (! rf->nn || ! rf->nr)
	continue ;
      keySet (chainScores, rf->chain) =  keySet (scores, rf->nr) ;
      keySet (chain2nr, rf->chain) = rf->nr ;
      jj = keySet (nr2chain, rf->nr) ;
      if (! jj)
	keySet (nr2chain, rf->nr) = rf->chain ;
    }
  /* attribute it to all the chains which are part of the raphia */
  for (ii = ii1, px = bigArrp (aa, ii, PEXPORT) ; ii < ii2 ; px++, ii++)
    {
      if (! px->score)
	continue ;
     
      nr = keySet (chain2nr, px->chain) ;
      px->chain =  keySet (nr2chain, nr) ;
      px->score = px->pairScore = keySet (chainScores, px->chain) ;
    }
  if (0)
    showRaphia (raphia) ;
 done:   
  ac_free (h) ;
  return nRF ;
}  /* clipAlignOptimizeHalfRaphia */

/*************************************************************************************/

static void clipAlignOptimizeSetChainAli  (CLIPALIGN *pp, BigArray aa, int ii1, int ii2)
{  
  AC_HANDLE h = ac_new_handle () ;
  KEYSET chainAlis, x1Alis, x2Alis ;
  PEXPORT *px ;
  int ii, dx, ii3 = ii1, old = -1 ;
  
  chainAlis = keySetHandleCreate (h) ;
  x1Alis = keySetHandleCreate (h) ;
  x2Alis = keySetHandleCreate (h) ;
  
  for (ii = ii1, px = bigArrp (aa, ii, PEXPORT) ; ii < ii2 ; px++, ii++)
    {
      if (! px->score)
	continue ;
      if (ii > ii1 && old >= 0 && px[0].mm != px[old-ii].mm)
	{
	  int i ;
	  PEXPORT *px2 ;
	  
	  ii3 = ii ;
	  for (i = ii1, px2 = bigArrp (aa, i, PEXPORT) ; i < ii ; px2++, i++)
	    {
	      if (! px2->score)
		continue ;
	      px2->chainAli = keySet (chainAlis, px2->chain) ;
	      px2->s1 = keySet (x1Alis, px2->chain) ;
	      px2->s2 = keySet (x2Alis, px2->chain) ;
	    }
	  chainAlis = keySetHandleCreate (h) ;
	  x1Alis = keySetHandleCreate (h) ;
	  x2Alis = keySetHandleCreate (h) ;
	  dx = 0 ;
	}
      else
	dx = px->x1 - keySet (x2Alis, px->chain) - 1 ;
      if (dx < 0)
	dx = -dx ;
      else 
	dx = 0 ;
      old = ii ;
      if (! keySet (x2Alis, px->chain))
	{
	  keySet (x1Alis, px->chain) = px->x1 ;
	  keySet (x2Alis, px->chain) = px->x2 ;
	}
      else 
	{
	  if (keySet (x1Alis, px->chain) > px->x1)
	    keySet (x1Alis, px->chain) = px->x1 ;
	  if (keySet (x2Alis, px->chain) < px->x2)
	    keySet (x2Alis, px->chain) = px->x2 ;
	}
      keySet (chainAlis, px->chain) +=  px->x2 - px->x1 + 1 - dx ;
    }
  
  for (ii = ii3, px = bigArrp (aa, ii, PEXPORT) ; ii < ii2 ; px++, ii++)
    {
      if (! px->score)
	continue ;
      px->chainAli = keySet (chainAlis, px->chain) ;
      px->s1 = keySet (x1Alis, px->chain) ;
      px->s2 = keySet (x2Alis, px->chain) ;
    }
  
  ac_free (h) ;
  return ; 
} /* clipAlignOptimizeSetChainAli */

/*************************************************************************************/

/* separate the 2 reads of the fragment */
static void clipAlignOptimizeChain  (CLIPALIGN *pp, BigArray aa, int ii1, int ii2, BOOL isRaphia)
{
  MM *mm = 0 ;
  int ii = 0, jj1 = ii1 ;
  PEXPORT *px ;
  BOOL isDown ;
  BOOL debug = FALSE ;
  int chain = 0 ;
  int target = 0 ;
  
  for (ii = ii1, px = bigArrp (aa, ii1, PEXPORT) ; ii < ii2 ; px++, ii++)
    {
      px->isDown = px->a1 < px->a2 ? 1 : 0 ;
      if (! target) { isDown = px->isDown ; target = px->target ; }
      if (! mm) mm = px->mm ;
      if (mm != px->mm || target != px->target || isDown != px->isDown)
	{
	  if (! isRaphia)
	    chain = clipAlignOptimizeHalfChain  (pp, aa, jj1, ii, chain) ;
	  else
	    clipAlignOptimizeHalfRaphia (pp, aa, jj1, ii) ;
	  jj1 = ii ; mm = px->mm ; target = px->target ;
	  isDown = px->isDown ;
	}
    }
  
  if (! isRaphia)
    clipAlignOptimizeHalfChain  (pp, aa, jj1, ii2, chain) ;
  else
    clipAlignOptimizeHalfRaphia (pp, aa, jj1, ii2) ; 
  clipAlignOptimizeSetChainAli (pp,  aa, ii1, ii2) ;
  if (debug) 
    showPexport (pp, 1) ;
  return ;
}  /* clipAlignOptimizeChain */

/*************************************************************************************/

static void clipAlignOptimizeOnePair (CLIPALIGN *pp, BOOL MRNAH, BOOL singleTarget, BitSet bb, BigArray aa, int ii1, int ii2
				      , int *t0p, int *w0p, int *d0p, int *c0p)
{
  int ii, jj ;
  PEXPORT * __restrict__ px, * __restrict__ py, * __restrict__ px1, * __restrict__ px2 ;
  BOOL isTranscript = ! pp->splice ;
  int maxDiscardableAli = pp->maxDiscardableAli ;
  BOOL debug = FALSE ;
  static int pass = 0 ;

  if (0 && ! singleTarget && (pass++ % 1000) == 0)
    fprintf (stderr, "-") ;
 
  if (singleTarget)
    clipAlignOptimizeCleanUp (pp, aa, ii1, ii2) ;

  px = bigArrp (aa, ii1, PEXPORT) ;
  if (0 && ! strcmp (stackText (pp->probeStack, px->mm->probeName), "seq.1149586>"))
    debug = TRUE ;
  if (0 && singleTarget)
    clipAlignSquish (pp,1) ;
  
  if (debug)
    {
      fprintf (stderr, "-------0-------clipAlignOptimizeOnePair starts: %s\n"
	       ,  pp->previousTargetGene ? dictName(pp->targetDict, pp->previousTargetGene) : ""
	       ) ;
      showPexport (pp, singleTarget ? 1 : 0) ;
    } 
 
  if (singleTarget) /* genome spliced case, score the exon/intron chains of each read */
    {
      clipAlignOptimizeChain (pp, aa, ii1, ii2, FALSE) ;  /* chain */
      if(1) clipAlignOptimizeChain (pp, aa, ii1, ii2, TRUE)  ;  /* raphia */
    }

  if (debug)
    {
      fprintf (stderr, "-------A-------clipAlignOptimizeOnePair starts: %s\n"
	       ,  pp->previousTargetGene ? dictName(pp->targetDict, pp->previousTargetGene) : ""
	       ) ;
      showPexport (pp, singleTarget ? 1 : 0) ;
    } 

  /* in transcript case, for each read in a given gene keep only the best score in that gene
   * discard suboptimal or multiple hits in the same transcript
   * or, in mRNAH case, keep only hits to the first transcript
   * this module does not look at pairs 
   * it does not matter later  if the best pair uses different transcripts of the same gene
   */
  if (singleTarget && isTranscript && pp->targetGene)  /* non spliced case : MRNA or genome and no splicing */
    {
      MM *mm1 = 0 ;
      bb = bitSetReCreate (bb, ii2) ;
      px = bigArrp (aa, ii1, PEXPORT) ;

      for (;;) /* forever */
	{
	  int score = 0, ali = 0, target, targetGene, chain, x1, x2 ;
	  int  bestii = -1 ;
	  
	  /* find the best score of the read */
	  for (ii = ii1,  py = px ; ii < ii2 ; py++, ii++)
	    {
	      if (! py->score || bit(bb, ii)) /* already done */
		continue ;
	      if (score < py->score || (score == py->score && ali < py->ali))
		{
		  score = py->score ;
		  mm1 = py->mm ;
		  ali = py->ali ;
		  x1 = px->x1 ;
		  x2 = px->x2 ;
		  chain = py->chain ;
		  target = py->target ;
		  targetGene = py->targetGene ;
		  bestii = ii ;
		}
	    }
	  if (bestii < 0)
	    break ;
	  
	  bitSet (bb, bestii) ;
	  /* kill the others */
	  for (ii = ii1,  py = px ; ii < ii2 ; py++, ii++)
	    if (py->score && py->mm == mm1 && ! bit (bb, ii))
	      {
		if (
		    score > py->score ||
		    (score == py->score && ali > py->ali) ||
		    (
		     score == py->score && ali == py->ali && pp->MRNAH && 
		     targetGene == py->targetGene && 
		     ( 
		      target < py->target || (target == py->target && chain != py->chain)
		       )
		     )
		    )
		  {
		    int c1 = x1 > py->x1 ? x1 : py->x1 ;
		    int c2 = x2 < py->x2 ? x2 : py->x2 ;
		    int dc = 3 *(c2 - c1 + 1) ;
		    int dy = py->x2 - py->x1 ;
		    if (dc > dy)       /* good overlap position kill it */
		      {
			py->score = 0 ;
		      }
		  }
	      }
	}
    }

  if (debug)
    {
      fprintf (stderr, "-------1-------clipAlignOptimizeOnePair after selection of the best scores for each read: %s\n"
	       ,  pp->previousTargetGene ? dictName(pp->targetDict, pp->previousTargetGene) : ""
	       ) ;
      showPexport (pp, singleTarget ? 1 : 0) ;
    } 

  /* construct the true pair score 
   * up to now px->pairSCore is the chain score
   */
  if (pp->hasPairs && singleTarget)
    {
      int bestPairScore = 0 ;
      int pass ;
      BOOL foundPair = FALSE ;

      /* first transfer the chain score from pairSCore to score */
      for (ii = ii1, px = bigArrp (aa, ii1, PEXPORT) ; ii < ii2 ; px++, ii++) 
	if (px->score)
	  {
	    px->score = px->pairScore ;
	    if (px->pairScore >  bestPairScore)
	      bestPairScore = px->pairScore ;
	  }

      /* then compute the true pair score in case of compatible pairs */
      bb = bitSetReCreate (bb, ii2) ;
      for (pass = 0 ;  pass < (singleTarget ? 2 : 1) ; pass++)
	{
	  for (;;) /* forever */
	    { 
	      int pairScore = 0 ;
	      int minD, iChain, dx, bestiChain, bestjChain, bestTarget, target, s1 ;
	      PEXPORT *pxa ;
	      MM *bestMm1 = 0, *bestMm2 = 0, *mm1 = bigArrp (aa, ii1, PEXPORT)->mm ;
	      
	      /* search for the shortest distance between 2 reads of a pair */
	      bestiChain = bestjChain = bestTarget = 0 ; iChain = -1 ;
	      for (ii = ii1, px = bigArrp (aa, ii, PEXPORT)  ; ii < ii2 ; ii++, px++)
		{
		  if (bit(bb, ii) || ! px->score || mm1 != px->mm) /* already analyzed */
		    continue ;
		  iChain = px->chain ;
		  pxa = px ; dx = px->a2 - px->a1 ;  minD = 999999999 ; 
		  target = (pass == 1 &&  singleTarget ? 0 : px->target) ;
		  for (jj = ii + 1, py = px + 1 ; jj < ii2 ; jj++, py++)
		    {
		      if (bit(bb, jj) || ! py->score ) /* already analyzed */
			continue ;
		      if (target && py->target != target)
			continue ;
		      if (py->mm == mm1)
			{
			  if (py->chain != iChain)
			    continue ;
			  else
			    pxa = py ; /* new exon of same chain */
			}
		      if (py->mm != mm1)
			{
			  long int dy, ddx1, ddx2 ;
			  /* is this pair compatible */
			  dy = py->a2 - py->a1 ;
			  if (target)
			    {
			      if (dx * dy > 0)
				continue ;
			      ddx1 = py->a1 - px->a1 ; 
			      ddx2 = py->a1 - pxa->a1 ; 
			      if (dx * ddx1 < 0)
				continue ;
			      if (ddx1 * ddx2 < 0) ddx1 = 0 ; /* since the y chain ends inside the x chain */
			      
			      if (ddx1 < 0) ddx1 = -ddx1 ; 
			      if (ddx2 < 0) ddx2 = -ddx2 ;
			      if (ddx2 < ddx1) ddx1 = ddx2 ;
			      if (ddx1 >= minD)
				continue ;
			      if (ddx1 > 1000000)  /* limit compatible pairs to 1 Mb */
				continue ;
			    }
			  s1 = px->score + py->score ; /* - (ddx1/100 < 12 ? ddx1/100 : 12) ; */
			  if (s1 >=  px->pairScore - px->s1 && s1 >=  py->pairScore - py->s1 )
			    {
			      minD = ddx1 ; bestiChain = px->chain ; bestjChain = py->chain ; /* best pair so far */
			      pairScore = px->score + py->score ;
			      bestMm1 = px->mm ; bestMm2 = py->mm ;
			      foundPair = TRUE ;
			      px->pairScore = py->pairScore = pairScore ;
			      /*     px->s1 = py->s1 =  (ddx1/100 < 12 ? ddx1/100 : 12) ; */
			    }
			}
		    }
		}
	      /*
		for (ii = ii1, px = bigArrp (aa, ii, PEXPORT)  ; ii < ii2 ; ii++, px++)
		px->s1 = 0 ;
	      */
	      if (! bestiChain) /* no new chain, this phase is over */
		break ;
	    
	      if (bestMm1->bestPairScore < pairScore)
		{
		  bestMm1->bestPairScore = pairScore ;	     
		  bestMm1->nGeneHits = 0 ;
		}
	      else if (bestMm1->bestPairScore > pairScore)
		pairScore = 0 ;
	      if (bestMm2->bestPairScore < pairScore)
		{
		  bestMm2->bestPairScore = pairScore ; 
		  bestMm2->nGeneHits = 0 ;
		}
	      else if (bestMm2->bestPairScore > pairScore)
		pairScore = 0 ;
	      
	      /* register this new chain */
	      for (ii = ii1, px = bigArrp (aa, ii, PEXPORT)  ; ii < ii2 ; ii++, px++)
		{
		  if (
		      ! bit (bb, ii) &&
		      (
		       (px->mm == bestMm1 && px->chain == bestiChain) ||
		       (px->mm == bestMm2 && px->chain == bestjChain)
		       )
		      )
		    {
		      px->chain = bestiChain ;
		      px->pairScore = pairScore ;
		      bitSet (bb, ii) ;
		    }
		}
	      if (pairScore > bestPairScore)
		bestPairScore = pairScore ;  /* this works for pairs and for orphan cases  */
	    } /* forever */
	  /* scan again and reject suboptimal pair scores */
	  if (foundPair)
	    {
	      px = bigArrp (aa, ii1, PEXPORT) ;
	      bestPairScore = px->mm->bestPairScore ;
	      for (px1 = px , jj = ii1 ; jj < ii2 ; px1++, jj++)
		if (px1->pairScore < bestPairScore)
		  px1->pairScore = px1->score = 0 ;
		else
		  px1->score =  px1->pairScore ;
	      break ; /* do not search for inter transcript pairs if we foun a pair inside a single transcript */
	    }	  
	}  /* pass */
    }
   
   
  if (debug) /*  && pp->hasPairs && isTranscript) */
    {
      fprintf (stderr, "-------2-------clipAlignOptimizeOnePair after slection of the best pair scores: %s\n"
	       ,  pp->previousTargetGene ? dictName(pp->targetDict, pp->previousTargetGene) : ""
	       ) ;
      showPexport (pp, singleTarget ? 1 : 0) ;
    } 

  if (! singleTarget)
    {
      /* scan again and reject suboptimal pair scores */
      MM *mm1 = 0 ;
      bb = bitSetReCreate (bb, ii2) ;
      px = bigArrp (aa, ii1, PEXPORT) ;

      for (;;) /* forever */
	{
	  int score = 0, ali = 0, target, targetGene, chain, x1, x2 ;
	  int  bestii = -1 ;
	  
	  /* find the best score of the read */
	  for (ii = ii1,  py = px ; ii < ii2 ; py++, ii++)
	    {
	      if (! py->score || bit(bb, ii)) /* already done */
		continue ;
	      if (score < py->pairScore || (score == py->pairScore && ali < py->ali))
		{
		  score = py->pairScore ;
		  mm1 = py->mm ;
		  ali = py->ali ;
		  x1 = py->s1 ;
		  x2 = py->s2 ;
		  chain = py->chain ;
		  target = py->target ;
		  targetGene = py->targetGene ;
		  bestii = ii ;
		}
	    }
	  if (bestii < 0)
	    break ;
	  if (px->mm->bestPairScore < score)
	    px->mm->bestPairScore = score ;
	  bitSet (bb, bestii) ;
	  /* kill the others */
	  for (ii = ii1,  py = px ; ii < ii2 ; py++, ii++)
	    if (py->pairScore && py->mm == mm1 && ! bit (bb, ii) &&
		(target != py->target || chain != py->chain)
		)
	      {
		if (
		    score > py->pairScore ||
		    (
		     score == py->pairScore && 
		     pp->MRNAH && 
		     targetGene == py->targetGene && 
		     target < py->target
		     )
		    )
		  {
		    int c1 = x1 > py->s1 ? x1 : py->s1 ;
		    int c2 = x2 < py->s2 ? x2 : py->s2 ;
		    int dc = 3 *(c2 - c1 + 1) ;
		    int dy = py->s2 - py->s1 ;
		    if (dc > dy)       /* good overlap position kill it */
		      {
			py->score = py->pairScore = 0 ;
		      }
		  }
	      }
	}
    }


  if (debug)
    {
      fprintf (stderr, "-------4-------clipAlignOptimizeOnePair: %s\n"
	       ,  pp->previousTargetGene ? dictName(pp->targetDict, pp->previousTargetGene) : ""
	       ) ;
      showPexport (pp, singleTarget ? 1 : 0) ;
    } 


  /* MRNAH : prefer multi hits to same transcript, else keep first transcrip for each read */
  if (singleTarget && isTranscript && pp->hasPairs && pp->MRNAH)  
    {
      int ok = 0 ;
      for (ii = ii1, px = bigArrp (aa, ii1, PEXPORT) ; ok < 3 && ii < ii2 ; px++, ii++)
	{
	  if (px->score)
	    ok |=  (px->mm->pair < 0 ? 0x1 : 0x2) ;
	}
      if (ok < 3) /* keep top transcript */
	{ 
	  if (pp->MRNAH)
	    {
	      for (ii = ii1, px = bigArrp (aa, ii1, PEXPORT) ; ii < ii2 ; px++, ii++)
		if (px->score) break ;
	      for (++ii, ++px ;  ii < ii2 ; px++, ii++)
		px->score = px->pairScore = 0 ;
	    }
	}
      else
	{ /* locate a transcript with 2 hits */
	  int tr, bestTr = 0 ;
	  for (ii = ii1, px = bigArrp (aa, ii1, PEXPORT) ; ! bestTr && ii < ii2 ; px++, ii++)
	    {
	      ok = 0 ;
	      if (! px->score) 
		continue ;
	      ok |=  (px->mm->pair < 0 ? 0x1 : 0x2) ;
	      tr = px->target ;
	      for (jj = ii + 1, py = px + 1 ;  ! bestTr && jj < ii2 ; jj++, py++)  
		{
		  if (! px->score ||
		      tr != py->target ||
		      ok == (px->mm->pair < 0 ? 0x1 : 0x2)
		      )
		    continue ;
		  bestTr = tr ;
		}
	    }
	  ok = 0 ; /* locate first good guy */
	  for (ii = ii1, px = bigArrp (aa, ii1, PEXPORT) ; ii < ii2 ; px++, ii++)
	    {
	      if (! px->score) 
		continue ;
	      if (bestTr && bestTr != px->target)
		{ px->score = px->pairScore = 0 ; continue ; }
	      ok |= (px->mm->pair < 0 ? 0x1 : 0x2) ;
	      break ;
	    }
	  /* kill all others */
	  if (ok)
	    for (++ii, ++px ; ii < ii2 ; px++, ii++)
	      {
		if (! px->score) 
		  continue ;
		if (bestTr && bestTr != px->target)
		  { px->score = px->pairScore = 0 ; continue ; }
		if ( ok & (px->mm->pair < 0 ? 0x1 : 0x2))
		  { px->score = px->pairScore = 0 ; continue ; }
		ok = 0x3 ;  /* will kill all others */
	      }
	}
    }
 
  /* strand counts */
  if (singleTarget && isTranscript)  
    {
      for (ii = ii1, px = bigArrp (aa, ii1, PEXPORT) ; ii < ii2 ; px++, ii++)
	{
	  if (px->score > 80 || 100 * px->score > 90 * px->ln)
	    { 
	      if (px->a1 < px->a2)
		{ if (px->mm->pair < 0) pp->nReadMinus++ ; else pp->nReadPlus++ ; }
	      else
		{ if (px->mm->pair < 0) pp->nReadPlus++ ; else pp->nReadMinus++ ; }
	    }
	}
    }
 
  /* compare different genes and discard lower scores */

  /* Disfavor genes A-and-B also denoted A--B */
  if (! singleTarget || isTranscript)
    {
      BitSet bdg = pp->doubleGenes ;
      
      if (bdg)
	for (ii = ii1, px = bigArrp (aa, ii, PEXPORT)  ; ii < ii2 ; ii++, px++)
	  {
	    MM *mm1 = 0 ;
	    int g1, g2 ;
	    if (px->score)
	      {
		mm1 = px->mm ;
		g1 = px->targetGene ;
		for (jj = ii + 1, py = px + 1 ; jj < ii2 ; jj++, py++)  /* protect against double counting mm */
		  if (py->mm == mm1 && py->score)
		    {
		      g2 = py->targetGene ;
		      if (bit (bdg, g1) && ! bit (bdg,g2))
			{ px->score = 0 ; break ; }
		      if (bit (bdg, g2) && ! bit (bdg,g1))
			py->score = 0 ;
		    }
	      }
	  }
    }

  /*  target multiplicity */
  if (singleTarget || ! isTranscript)
    {
      bb = bitSetReCreate (bb, ii2) ;
      for (ii = ii1, px = bigArrp (aa, ii, PEXPORT)  ; ii < ii2 ; ii++, px++)
	{
	  MM *mm1 = 0 ;
	  if (px->score && ! bit(bb, ii))  /* count each mm only once */
	    {
	      mm1 = px->mm ;
	      for (jj = ii, py = px ; jj < ii2 ; jj++, py++)  /* protect against double counting mm */
		if (py->mm == mm1 && py->score)
		  bitSet (bb, jj) ;
	      mm1->nGeneHits++ ;  /* used to prune the search in future targets */
	    }
	}
    }
  else 
    {   /* ! singleTarget && isTranscript
	   Transcript case count the target multiplicity, if possible select the strand 
	*/
      bb = bitSetReCreate (bb, ii2) ;
      for (ii = ii1, px = bigArrp (aa, ii1, PEXPORT) ; ii < ii2 ; px++, ii++)
	if (px->score && px->mm && ! bit (bb, ii))
	  {
	    int s = 1, nGenePlus = 0,  nGeneMinus = 0, nGenePlusU = 0,  nGeneMinusU = 0, targetGene = 0 ;
	    MM *mm = px->mm ;
	    BOOL isGeneSpliced = FALSE ;
	    int x1 = px->x1 ;
	    int x2 = px->x2 ;
	    int dx = x2 - x1 + 1 ;
	    bitSet (bb, ii) ;
	    
	    for (px1 = px , jj = ii ; jj < ii2 ; px1++, jj++)
	      if (px1->score && px1->mm == mm)
		{
		  int c1 = px1->x1 > x1 ? px1->x1 : x1 ;
		  int c2 = px1->x2 < x2 ? px1->x2 : x2 ;
		  int dc = 3 *(c2 - c1 + 1) ;
		  int dx2 = px1->x2 - px1->x1 + 1 ;
		  if (dc < dx && dc < dx2)
		    continue ;
		
		  bitSet (bb, jj) ;
		  if (targetGene != px1->targetGene)
		    {
		      targetGene = px1->targetGene ;
		      isGeneSpliced = bit (pp->isSplicedTarget, px1->targetGene) ;
		      if (px1->a1 < px1->a2)
			{ 
			  if (mm->pair < 0) 
			    { nGeneMinus++ ; if (! isGeneSpliced) nGeneMinusU++ ; }
			  else 
			    { nGenePlus++ ; if (! isGeneSpliced) nGenePlusU++ ; }
			}
		      else
			{
			  if (mm->pair < 0) 
			    { nGenePlus++ ; if (! isGeneSpliced) nGenePlusU++ ; }
			  else  
			    { nGeneMinus++ ; if(! isGeneSpliced) nGeneMinusU++ ; }
			}
		    }
		}
	    if (nGenePlus && nGeneMinus)
	      {
		s = 0 ;
		if (pp->nReadPlus  > 4000 && pp->nReadPlus  > 4 * pp->nReadMinus)
		  { s = 1 ;  nGeneMinus = 0 ; nGeneMinusU = 0 ; } /* only keep plus genes */
		if (pp->nReadMinus > 4000 && pp->nReadMinus > 4 * pp->nReadPlus)
		  { s = -1 ; nGenePlus = 0 ; nGenePlusU = 0 ; } /* only keep minus genes */
		if (! strcmp ("n.6109037#1>", stackText (pp->probeStack, mm->probeName)))
		  fprintf (stderr, ".... n.6109037#1 s=%d rPlus=%d rMinus=%d\n",s, pp->nReadPlus, pp->nReadMinus) ;


		if (s)
		  for (px1 = px , jj = ii ; jj < ii2 ; px1++, jj++)
		    if (px1->score && px1->mm == mm &&
			(
			 ( mm->pair >= 0 && s * (px1->a2 - px1->a1) < 0) ||
			 ( mm->pair < 0 && s * (px1->a2 - px1->a1) > 0)
			 )
			)
		      px1->score = px1->pairScore = 0 ;
	      }
	    /* disfavor unspliced genes relative to spliced genes, but only on the same strand */
	    if (pp->avoidPseudoGenes && nGenePlusU && nGenePlus > nGenePlusU)
	      {
		nGenePlus -= nGenePlusU ;
		for (px1 = px , jj = ii ; jj < ii2 ; px1++, jj++)
		  if (px1->score && px1->mm == mm &&
		      mm->pair * (px->a2 -  px->a1) > 0 &&
		      ! bit(pp->isSplicedTarget, px1->targetGene)
		      )
		    if (1 || px1->pairScore < maxDiscardableAli)
		      px1->score = px1->pairScore = 0 ;
	      }
	    if (pp->avoidPseudoGenes && nGeneMinusU && nGeneMinus > nGeneMinusU)
	      {
		nGeneMinus -= nGeneMinusU ;
		for (px1 = px , jj = ii ; jj < ii2 ; px1++, jj++)
		  if (px1->score && px1->mm == mm &&
		      mm->pair * (px->a2 -  px->a1) < 0 &&
		      ! bit(pp->isSplicedTarget, px1->targetGene)
		      )
		    if (1 || px1->pairScore < maxDiscardableAli)
		      px1->score = px1->pairScore = 0 ;
	      }

	    /* register */
	    s = nGenePlus + nGeneMinus ;
	    px->mm->nGeneHits = s ;        /* nGeneHits is used in the final report */
	    if (nGenePlus && nGeneMinus) s = - s ;
	    for (px1 = px , jj = ii ; jj < ii2 ; px1++, jj++)
	      if (px1->score)
		px1->uu = s ;
	  }
    }

  if (debug) /*  && ! singleTarget && pp->hasPairs && pp->splice */
    {
      fprintf (stderr, "-------5-------clipAlignOptimizeOnePair after global genomic selection: %s\n"
	       ,  pp->previousTargetGene ? dictName(pp->targetDict, pp->previousTargetGene) : ""
	       ) ;
      showPexport (pp, singleTarget ? 1 : 0) ;
    } 
  
  /* genome  case count the target multiplicity, by counting  the non overlapping best segments*/
  if (! singleTarget &&  ! (isTranscript && pp->strategy == STRATEGY_RNA_SEQ))
    {
      PEXPORT *px3 ;
      int j2, j3 ;
      int nGenePlus = 0,  nGeneMinus = 0 ;

      bb = bitSetReCreate (bb, ii2) ;
      for (ii = ii1, px = bigArrp (aa, ii1, PEXPORT) ; ii < ii2 ; px++, ii++)
	if (px->score && px->mm && ! bit (bb, ii))
	  {
	    MM *mm = px->mm ;
	    int uu , bestScore = 0 ;
	    int bestx1 = 999999, bestx2 = 0, besty1 = 999999, besty2 = 0 ;

	    /* notice that pairScore in genomic case is the chain score */
	    uu = 0 ; 
	    if (1)
	      {
		for (px1 = px, jj = ii ; jj < ii2 ; px1++, jj++)
		  if (px1->score && px1->mm == mm)
		    {
		      if (pp->hasPairs)
			{
			  if (bestScore <= px1->pairScore)
			    { 
			      bestScore = px1->pairScore ; 
			      if (px->mm->pair < -1)
				{
				  if (bestx1 > px1->s1)  bestx1 = px1->s1 ; 
				  if (bestx2 < px1->s2)  bestx2 = px1->s2 ; 
				}
			      else
				{
				  if (besty1 > px1->s1)  besty1 = px1->s1 ; 
				  if (besty2 < px1->s2)  besty2 = px1->s2 ; 
				}
			    }
			}
		      else
			{
			  if (bestScore < px1->score)
			    {
			      bestScore = px1->score ; 
			      if (bestx1 > px1->s1)  bestx1 = px1->s1 ; 
			      if (bestx2 < px1->s2)  bestx2 = px1->s2 ; 
			    }
			}
		    }
	      }

	    if ((1 || ! pp->bestHit) && bestScore > maxDiscardableAli)   /* mieg 2022-12-10   allow muti hits also in transcript to allow fusions */
	      bestScore = maxDiscardableAli ;
	    if (bestScore)
	      for (px1 = px , jj = ii  ; jj < ii2 ; px1++, jj++) 
		{
		  int score ;
		  if (! px1->score || px1->mm != mm)
		    continue ;
		  score = pp->hasPairs ? px1->pairScore : px1->score ;
		  if (score < bestScore)
		    score = px1->score = px1->pairScore = 0 ;
		  if (score && score < bestScore)
		    {
		      int u1, u2 ;
		      if (! pp->hasPairs || px->mm->pair < -1)
			{ u1 = bestx1 ; u2 = bestx2 ; }
		      else
			{ u1 = besty1 ; u2 = besty2 ; }
		      if (u1 < px1->s1) u1 = px1->s1 ;
		      if (u2 >  px1->s2) u2 = px1->s2 ;
		      if (2*(u2 - u1) >  px1->s2 - px1->s1) 
			px1->score = px1->pairScore = 0 ;
		    }
		  if (! px1->score)
		    continue ;
		  if (! bit (bb, jj))
		    {	
		      for (px2 = px1 , j2 = jj ; j2 < ii2 ; px2++, j2++)
			if (px2->score == px1->score && px2->mm == mm && ! bit (bb, j2))  /* last term was ! bit (bb, jj which seems dubious 2018+03_21 */
			  {
			    for (px3 = px2 , j3 = j2 ; j3 < ii2 ; px3++, j3++)
			      if (px3->mm == mm && px3->score == px1->score)
				{
				  bitSet (bb, j3) ;
				  if (
				      (px3->target == px2->target   && px3->x1 < px->x2 - 10 && px3->x2 > px->x1 + 10) /* overlap on same target */
				      || (px3->targetGene != px2->targetGene)                        /* 2 genes */
				      || (px3->targetGene == 0 && px3->target != px2->target )       /* 2 targets (say 2 chromosomes), excluding 2 mrna of the same gene */
				      )
				    {
				      uu++ ;
				      if (px1->a1 < px1->a2)
					{ 
					  if (mm->pair < 0) 
					    nGeneMinus++ ;
					  else 
					    nGenePlus++ ;
					}
				      else
					{
					  if (mm->pair >= 0) 
					    nGeneMinus++ ;
					  else 
					    nGenePlus++ ;			
					}
				    }
				}
			  }
		    }
		}
	    bitSet (bb, ii) ;
	    if (nGenePlus && nGeneMinus) uu = -uu ;
	    for (px1 = px , jj = ii ; jj < ii2 ; px1++, jj++)
	      if (px1->mm == mm && px1->score)
		{ px1->uu = uu ; bitSet (bb, jj) ; }
	  }
    }

  if (debug)
    {
      fprintf (stderr, "-------6-------clipAlignOptimizeOnePair after global gene selection: %s\n"
	       ,  pp->previousTargetGene ? dictName(pp->targetDict, pp->previousTargetGene) : ""
	       ) ;
      showPexport (pp, singleTarget ? 1 : 0) ;
    } 
 

  /* clean up the donor/acceptor/introns lists */
  if (singleTarget)
    clipAlignOptimizeIntronCleanUp (pp, aa, ii1, ii2, singleTarget, t0p, w0p, d0p, c0p) ;


  /* reset the chainAli value */
  if (singleTarget)
    clipAlignOptimizeSetChainAli (pp,  aa, ii1, ii2) ;
  if ( debug)
    {
      fprintf (stderr, "-------7-------after global gene strand selection: %s\n"
	       ,  pp->previousTargetGene ? dictName(pp->targetDict, pp->previousTargetGene) : ""
	       ) ;
      showPexport (pp, singleTarget ? 1 : 0) ;
    } 
}  /* clipAlignOptimizeOnePair */

/*************************************************************************************/

static int clipAlignOptimizePairs (CLIPALIGN *pp, BOOL singleTarget)
{
  PEXPORT * Restrict px ;
  BigArray aa = singleTarget ? pp->exportGeneHits : pp->exportHits ;
  int gene = 0, fragment = 0 ;
  int iOld = 0, ii, iMax = bigArrayMax (aa) ;
  int t0 = 0, w0 = 0, d0 = 0, c0 = 0 ;
  BOOL debug = FALSE ;
  AC_HANDLE h = ac_new_handle () ;
  BitSet bb = bitSetCreate (iMax, h) ;
  BitSet bdg = pp->doubleGenes ;
  
  if (bdg)
    {
      int n2 = dictMax (pp->targetDict) ;
      if (pp->doubleGenes && n2 && bitt (pp->doubleGenes, n2)) {} ; /* adjust size of bit array */
    }

  if (iMax)
    for (ii = 0, px = bigArrp (aa, ii, PEXPORT) ; ii < iMax ; px++, ii++)
      px->isDown = px->a1 < px->a2 ? 1 : 0 ;
  bigArraySort (aa,  pExportPairXOrder) ;
  bigArrayCompress (aa) ; 
  
  bitUnSet (bb, iMax) ; /* make room */
  iMax = bigArrayMax (aa) ;
  for (ii = 0, px = bigArrp (aa, 0, PEXPORT) ; ii < iMax ; ii++, px++)
    {
      if (px->fragment != fragment || (singleTarget && px->targetGene != gene))
	{
	  if (ii > iOld)
	    clipAlignOptimizeOnePair (pp, pp->MRNAH, singleTarget, bb, aa, iOld, ii, &t0, &w0, &d0, &c0) ; 
	  fragment = px->fragment ;
	  gene = px->targetGene ;
	  iOld = ii ;
	}
    }
  

  if (! singleTarget)
    {
      BigArray dd = pp->exportDoubleIntrons ;
      BigArray bb = pp->exportIntrons ;
      BigArray donors = pp->exportDonors ;
      BigArray ac = pp->exportAcceptors ;
      int iMax = bigArrayMax (bb) ;

      if  (iMax)
	for (ii = 0, px = bigArrp (bb, 0, PEXPORT) ; ii < iMax ; ii++)
	  px[ii].score *= -1 ;
      
      iMax =  bigArrayMax (donors) ;
      if  (iMax)
	for (ii = 0, px = bigArrp (donors, 0, PEXPORT) ; ii < iMax ; ii++)
	  px[ii].score *= -1 ;

      iMax =  bigArrayMax (ac) ;
      if  (iMax)
	for (ii = 0, px = bigArrp (ac, 0, PEXPORT) ; ii < iMax ; ii++)
	  px[ii].score *= -1 ;

      iMax =  bigArrayMax (aa) ;
      if  (iMax)
	for (ii = 0, px = bigArrp (aa, 0, PEXPORT) ; ii < iMax ; ii++)
	  {
	    if (px[ii].score && px[ii].intron)
	      {
		PEXPORT *py = bigArrp (bb, px[ii].intron - 1, PEXPORT) ;
		if (py->score < 0)  py->score *= -1 ;
	      }     
 	    if (px[ii].score && px[ii].donor)
	      {
		PEXPORT *py = bigArrp (donors, px[ii].donor - 1, PEXPORT) ;
		if (py->score < 0)  py->score *= -1 ;
	      }     
 	    if (px[ii].score && px[ii].acceptor)
	      {
		PEXPORT *py = bigArrp (ac, px[ii].acceptor - 1, PEXPORT) ;
		if (py->score < 0)  py->score *= -1 ;
	      }  
	  }   
      iMax =  bigArrayMax (ac) ;
      if  (iMax)
	for (ii = 0, px = bigArrp (ac, 0, PEXPORT) ; ii < iMax ; ii++)
	  if (px[ii].score < 0)
	    px[ii].score = 0 ;
      
      iMax = bigArrayMax (dd) ;
      if  (iMax)
	for (ii = 0, px = bigArrp (dd, 0, PEXPORT) ; ii < iMax ; ii++)
	  {
	    if (px[ii].intron && px[ii].doubleIntron)
	      {
		PEXPORT *py = bigArrp (bb, px[ii].intron - 1, PEXPORT) ;
		PEXPORT *pz = bigArrp (bb, px[ii].doubleIntron - 1, PEXPORT) ;
		px[ii].score = (py->score * pz->score  > 0) ? px[ii].score  : 0 ;
	      }
	    else
	      px[ii].score = 0 ;
	  }
      return 1 ;
    }
 
  if (ii > iOld) /* last line is always zero, but in genome case nothing happens before here */
    clipAlignOptimizeOnePair (pp, pp->MRNAH, singleTarget, bb, aa, iOld, ii, &t0, &w0, &d0, &c0) ; 

  clipAlignSquish (pp, singleTarget) ;
  
  if (debug)
    {
      fprintf (stderr, ".....clipAlignOptimizePairs: %s\n"
	       ,  pp->previousTargetGene ? dictName(pp->targetDict, pp->previousTargetGene) : ""
	       ) ;
      showPexport (pp, singleTarget) ;
    } 
  ac_free (h) ;
  return bigArrayMax (pp->exportGeneHits) ;
} /* clipAlignOptimizePairs */

/*************************************************************************************/

static void clipAlignGetDonors (CLIPALIGN *pp)
{						
  int iMax = bigArrayMax (pp->exportGeneHits) ; /* number of objects to copy */
  int ii ;
  PEXPORT *px, *pxg ;

  for (ii = 0,  pxg = bigArrp (pp->exportGeneHits, ii, PEXPORT) ;
       ii <  iMax ; ii++, pxg++)
    {
      if (pxg->score)
	{
	  char *cp, buf[64 + OVLN] ;
	  char type =  pxg->presuffix & 0xff ;

	  strncpy (buf, dictName(pp->exportDict, pxg->suffix), 63 + OVLN) ;
	  if (buf[0] && type &&
	      (type != 1 || type != 63 || slideCheckHookEntropy((unsigned char *)buf))
	      )
	    { 
	      int jj = bigArrayMax (pp->geneDonors) ;
	      px = bigArrayp (pp->geneDonors, jj, PEXPORT) ;
	      *px = *pxg ; px->score =  pxg->score ; 
	      px->type = type ;
	      px->presuffix =  pxg->suffix ;
	      px->donor = px->acceptor = 0 ;
	      if (jj && ! memcmp(px, px-1, sizeof (PEXPORT)))
		bigArrayMax (pp->geneDonors) = jj-- ;
	      if (pxg->targetPresuffix)
		{ 
		  cp = strstr (dictName(pp->exportDict, pxg->targetPresuffix), "*") ;
		  if (cp && cp[1]) 
		    {
		      strncpy (buf,  cp+1, 63 + OVLN) ; 
		      cp = strstr (buf, "#") ; if (cp) *cp = 0 ;
		      dictAdd (pp->exportDict, buf, &(px->targetPresuffix)) ;
		    }
		}
	      else
		px->targetPresuffix = 0 ;
	      if (px->mm->pair < -1)
		{
		  int kk = bigArrayMax (pp->geneAcceptors) ;
		  PEXPORT *py = bigArrayp (pp->geneAcceptors, kk, PEXPORT) ;
		  int x ;
		  MM *mm1 = arrp (pp->probes, -px->mm->pair - 2, MM) ;
		  x = px->a1 ; px->a1 = px->a2 ; px->a2 = x ;
		  x = px->x1 ; px->x1 = px->mm->probeLength  + mm1->probeLength - px->x2 + 1 ; px->x2 =  px->mm->probeLength  + mm1->probeLength - x + 1 ;
		  x = px->prefix ; px->prefix = px->suffix ; px->suffix = x ;	
		  px->mm =  mm1 ;
		  *py = *px ; py->donor = py->acceptor = 0 ;
		  if (kk && ! memcmp(py, py-1, sizeof (PEXPORT)))
		    bigArrayMax (pp->geneAcceptors) = kk-- ;
		  pxg->acceptor = kk + 1 ;
		  bigArrayMax (pp->geneDonors)-- ;
		}
	      else
		pxg->donor = jj + 1 ;
	    }
	}
    }

  return ;
} /* clipAlignGetDonors */

/*************************************************************************************/

static void clipAlignGetAcceptors (CLIPALIGN *pp)
{
  int iMax = bigArrayMax (pp->exportGeneHits) ; /* number of objects to copy */
  int ii ;
  PEXPORT *px, *pxg ;

  for (ii = 0,  pxg = bigArrp (pp->exportGeneHits, ii, PEXPORT) ;
       ii <  iMax ; ii++, pxg++)
    {
      if (pxg->score)
	{
	  char *cp, buf[64 + OVLN] ;
	  char type =   (pxg->presuffix >> 8) & 0xff ;
	  
	  strncpy (buf, dictName(pp->exportDict, pxg->prefix), 63 + OVLN) ;
	  if (buf[0] && type &&
	      (type != 1 || type != 63 ||  slideCheckHookEntropy((unsigned char *)buf))
	      )
	    { 
	      int jj = bigArrayMax (pp->geneAcceptors) ;
	      px = bigArrayp (pp->geneAcceptors, jj, PEXPORT) ;
	      *px = *pxg ; px->score =  pxg->score ;
	      px->type = type ;
	      px->presuffix =  pxg->prefix ;
	      px->donor = px->acceptor = 0 ;
	      if (jj && ! memcmp(px, px-1, sizeof (PEXPORT)))
		bigArrayMax (pp->geneAcceptors) = jj-- ;
	      if (pxg->targetPresuffix)
		{ 
		  strncpy (buf,  dictName(pp->exportDict, pxg->targetPresuffix), 63 + OVLN) ;
		  cp = strstr (buf, "*") ;
		  if (cp) *cp = 0 ;
		  dictAdd (pp->exportDict, buf, &(px->targetPresuffix)) ;
		}
	      else
		px->targetPresuffix = 0 ;
	      if (px->mm->pair < -1)
		{
		  int kk = bigArrayMax (pp->geneDonors) ;
		  PEXPORT *py = bigArrayp (pp->geneDonors, kk, PEXPORT) ;
		  int x ;
		  MM *mm1 = arrp (pp->probes, -px->mm->pair - 2, MM) ;
		  x = px->a1 ; px->a1 = px->a2 ; px->a2 = x ;
		  x = px->x1 ; px->x1 = px->mm->probeLength  + mm1->probeLength - px->x2 + 1 ; px->x2 =  px->mm->probeLength  + mm1->probeLength - x + 1 ;
		  x = px->prefix ; px->prefix = px->suffix ; px->suffix = x ;	
		  px->mm =  mm1 ;
		  *py = *px ; py->donor = py->acceptor = 0 ;
		  pxg->donor = kk + 1 ;
		  bigArrayMax (pp->geneAcceptors)-- ;
		}
	      else
		pxg->acceptor = jj + 1 ;
	    }
	}
    } 

  return ;
} /* clipAlignGetAcceptors */

/*************************************************************************************/

static int clipAlignOptimizeMergeDonors (CLIPALIGN *pp)
{
  PEXPORT *px, *py, *pz ;
  int ii, iiMax = bigArrayMax (pp->exportGeneHits) ;
  int jj = bigArrayMax (pp->exportDonors) ;

  if (iiMax)
    for (ii = 0, px = bigArrp (pp->exportGeneHits, 0, PEXPORT) ; ii < iiMax ; ii++, px++)
      if (px->score && px->donor)
	{ 
	  py = bigArrp (pp->geneDonors, px->donor - 1, PEXPORT) ;
	  pz = bigArrayp (pp->exportDonors, jj++, PEXPORT) ;
	  px->donor = jj ;
	  *pz = *py ; pz->score = 0 ;
	}  
  return jj ;
} /* clipAlignOptimizeMergeDonors */

/*************************************************************************************/

static int clipAlignOptimizeMergeAcceptors (CLIPALIGN *pp)
{
  PEXPORT *px, *py, *pz ;
  int ii, iiMax = bigArrayMax (pp->exportGeneHits) ;
  int jj = bigArrayMax (pp->exportAcceptors) ;

  if (iiMax)
    for (ii = 0, px = bigArrp (pp->exportGeneHits, 0, PEXPORT) ; ii < iiMax ; ii++, px++)
      if (px->score && px->acceptor)
	{ 
	  py = bigArrp (pp->geneAcceptors, px->acceptor - 1, PEXPORT) ;
	  pz = bigArrayp (pp->exportAcceptors, jj++, PEXPORT) ;
	  px->acceptor = jj ;
	  *pz = *py ; pz->score = 0 ;
	}  
  return jj ;
} /* clipAlignOptimizeMergeAcceptors */

/*************************************************************************************/

static int clipAlignOptimizeMergeIntrons (CLIPALIGN *pp)
{
  const int ipxMax = 32 ;
  PEXPORT *px, *py, *pz, *pw, *lastPx = 0, *lastPz = 0, *pxs[ipxMax] ;
  int ii, iiMax = bigArrayMax (pp->geneIntrons), ipx ;
  int jj = bigArrayMax (pp->exportIntrons) ;
  int kk = bigArrayMax (pp->exportDoubleIntrons) ;

  memset (pxs, 0, sizeof(pxs)) ;
  if (iiMax)
    for (ii = 0, px = bigArrp (pp->geneIntrons, 0, PEXPORT) ; ii < iiMax ; ii++, px++)
      px->intron = 0 ;
  
  iiMax = bigArrayMax (pp->exportGeneHits) ;
  if (iiMax)
    for (ii = 0, px = bigArrp (pp->exportGeneHits, 0, PEXPORT) ; ii < iiMax ; ii++, px++)
      {
	ipx = 0 ;
	if (px->score)
	  {
	    if (px->intron)
	      { 
		py = bigArrp (pp->geneIntrons, px->intron - 1, PEXPORT) ;
		pz = bigArrayp (pp->exportIntrons, jj++, PEXPORT) ;
		lastPz = ii && jj > 1 ? bigArrayp (pp->exportIntrons, jj - 2, PEXPORT) : 0 ;
		px->intron = py->intron = jj ;
		*pz = *py ; pz->chain = px->chain ;
		if (1) /* alreadydone  above */
		  {
		    if (! strcmp (py->foot, "gt_ag"))
		      {
			if (py->errLeftVal > py->errRightVal)
			  pp->nIntronPlus += px->mm->mult ;
			else
			  {
			    pp->nIntronMinus  += px->mm->mult ;
			    if (pp->nIntronMinus == 1 || pp->nIntronMinus == 11)
			      invokeDebugger () ;
			  }
		      }
		    if (! strcmp (py->foot, "ct_ac"))
		      {
			if (py->errLeftVal > py->errRightVal)
			  {
			    pp->nIntronMinus += px->mm->mult ;
			    if (pp->nIntronMinus == 1 || pp->nIntronMinus == 11)
			      invokeDebugger () ;
			  }
			else
			  pp->nIntronPlus  += px->mm->mult ;
		      }
		  }
		for (ipx = 0 ; ipx < ipxMax ; ipx++)
		  {
		    lastPx = pxs[ipx] ;
		    if (lastPx && (px->fragment !=  lastPx->fragment || px->chain !=  lastPx->chain))
		      break ;
		    if (lastPx && px->mm == lastPx->mm && ipx)
		      continue ;
		    if (lastPx && lastPx->chain == px->chain && lastPx->mm && px->mm && lastPx->mm->fragment == px->mm->fragment)
		      {
			if (lastPx->intron && px->intron)
			  {
			    lastPx->donor = px->acceptor = 0 ;
			    lastPz = bigArrayp (pp->exportIntrons, lastPx->intron - 1, PEXPORT) ;
			    if (lastPz &&
				lastPz && lastPz->mm && pz->mm && lastPz->mm->fragment == pz->mm->fragment &&
				lastPz->chain == pz->chain &&
				(
				 (pz->a1 < pz->a2 && lastPz->a1 < lastPz->a2 && lastPz->a2 < pz->a1 && lastPz->x2 > pz->x1 - 20) ||			 
				 (pz->a1 < pz->a2 && lastPz->a1 < lastPz->a2 && pz->a2 < lastPz->a1 && pz->x2 > lastPz->x1 - 20) ||
				 (pz->a1 > pz->a2 && lastPz->a1 > lastPz->a2 && lastPz->a2 > pz->a1 && lastPz->x2 < pz->x1 + 20) ||			 
				 (pz->a1 > pz->a2 && lastPz->a1 > lastPz->a2 && pz->a2 > lastPz->a1 && pz->x2 < lastPz->x1 + 20)
				 ) 
				)
			      {
				pw = bigArrayp (pp->exportDoubleIntrons, kk++, PEXPORT) ;
				pw->target = pz->target ;
				pw->score = pz->mm->mult ;
				pw->intron = lastPz->intron ;
				pw->doubleIntron = pz->intron ;
				pw->fragment = pz->fragment ;
				pw->nN = pz->nN ;          /* strand */
				strncpy(pw->foot, pz->foot, 6) ;
				pw->mm = pz->mm ;
				if (
				    (lastPz->x1 < lastPz->a1 && lastPz->a2 < pz->a1) ||
				    (lastPz->x1 > lastPz->a1 && lastPz->a2 > pz->a1)
				    )
				  {
				    pw->a1 = lastPz->x1 ;      /* exon 1: a1-a2 */
				    pw->a2 = lastPz->a1 ;
				    pw->x1 = lastPz->a2 ;      /* exon 2: x1-x2 */
				    pw->x2 = pz->a1 ;     
				    pw->s1 = pz->a2 ;           /* exon 3: s1-s2 */
				    pw->s2 = pz->x2 ;
				  }
				else
				  {
				    pw->a1 = pz->x1 ;      /* exon 1: a1-a2 */
				    pw->a2 = pz->a1 ;
				    pw->x1 = pz->a2 ;      /* exon 2: x1-x2 */
				    pw->x2 = lastPz->a1 ;     
				    pw->s1 = lastPz->a2 ;           /* exon 3: s1-s2 */
				    pw->s2 = lastPz->x2 ;
				  } 
				if (lastPx && px->mm != lastPx->mm)
				  break ; /* at most a single read1 to read2 double intron */
			      }
			  }
		      }
		    else
		      break ;
		  }
	      }
	  }
	ipx++ ;
	if (ipx >= ipxMax)
	  ipx = ipxMax - 1;
	memmove (pxs+1,pxs, ipx * sizeof(PEXPORT*)) ;
	pxs[0] = px ;
      }

  if (iiMax)
    for (ii = 0, px = bigArrp (pp->exportGeneHits, 0, PEXPORT) ; ii < iiMax ; ii++, px++)
      {
	if (px->score && px->donor)
	  {
	    py = bigArrp (pp->exportDonors, px->donor - 1, PEXPORT) ;
	    py->score = px->mm->mult ;
	  }
	if (px->score && px->acceptor)
	  {	
	    py = bigArrp (pp->exportAcceptors, px->acceptor - 1, PEXPORT) ;
	    py->score = px->mm->mult ;
	  }
      }

  return jj ;
} /* clipAlignOptimizeMergeIntrons fragment */

/*************************************************************************************/
/* merge the geneHits into the hits, while resetting the donors */
static int clipAlignMergeHits (CLIPALIGN *pp)
{
  int nn = 0, ii, jj ;
  PEXPORT *px, *py, *pxg ;
  int iMax = bigArrayMax (pp->exportGeneHits) ; /* number of objects to copy */
  int jMax = bigArrayMax (pp->exportHits) ;
  BOOL debug = FALSE ;

  /* merge back in the global tables and destroy pp->geneDonors/Acceptors */
  clipAlignOptimizeMergeDonors (pp) ;
  clipAlignOptimizeMergeAcceptors (pp) ;
  bigArrayDestroy (pp->geneDonors) ;
  bigArrayDestroy (pp->geneAcceptors) ;

  if (pp->splice)
    {
      clipAlignOptimizeMergeIntrons (pp) ;
      bigArrayDestroy (pp->geneIntrons) ;
    }


  for (ii = 0,  jj = jMax, pxg = iMax ? bigArrp (pp->exportGeneHits, ii, PEXPORT) : 0 ;
       ii <  iMax ; ii++, pxg++)
    {
      if (pxg->score)
	{
	  px = bigArrayp (pp->exportHits, jj, PEXPORT) ;
	  *px = *pxg ; 
	  
	  /* remap so that the mm chain points inside exportGeneHits */
	  if (jj &&  px->mm->ipx == jj) 
	    invokeDebugger () ;
	  if (px->mm->ipx >= jj)
	    px->mm->ipx = 0 ;
	  if (px->mm->ipx)
	    {  /* double check we do not create a wrong chain */
	      py = bigArrayp (pp->exportHits, px->mm->ipx, PEXPORT) ;
	      if (py->mm != px->mm)
		px->mm->ipx = 0 ;
	    }
	  px->ipxold = px->mm->ipx ; 
	  px->mm->ipx = jj ; 
	  jj++ ;
	}
      if (pxg->mm)
	pxg->mm->ipxg = 0 ;
    }

  if (debug)
    {
      fprintf (stderr, "clipAlignMergeHits\n") ;
      showPexport (pp, 0) ;
      fprintf (stderr, ".....clipAlignMergeHits\n") ;
      showPexport (pp, 1) ;
    }
  if (0)
    {
      fprintf (stderr, ".....clipAlignMergeHits hits:%ld  geneHits:%ld donors:%ld acc:%ld %d plus %d minus\n"
	       , bigArrayMax (pp->exportHits)
	       , bigArrayMax (pp->exportGeneHits)
	       , bigArrayMax (pp->exportDonors)
	       , bigArrayMax (pp->exportAcceptors) 
	       , pp->nReadPlus
	       , pp->nReadMinus
	       ) ;
    }
  
  bigArrayMax (pp->exportGeneHits) = 1 ; /* ATTENTION : we do not zero it, this may be dangerous */
  px = bigArrp (pp->exportGeneHits, 0, PEXPORT) ;
  memset (px, 0, sizeof (PEXPORT)) ;

  return nn ;
} /* clipAlignMergeHits */

/*************************************************************************************/
/* compute the introns one fragment at a time and optimize that fragemnt chain
 * then clean out the sub optimal donor/acceptor/exons/introns
 * before registering them
 */

static int clipAlignOptimizeGetIntrons (CLIPALIGN *pp)
{ 
  int nn = 0 ;

  pp->geneIntrons = bigArrayCreate (1000, PEXPORT) ;
  
  nn = clipAlignConstructIntrons (pp,  1, TRUE) ;  /* F read on F strand or R-read on R-strand */
  nn += clipAlignConstructIntrons (pp, -1, TRUE) ; /* F read on R strand or R-read on F-strand */
  
  return nn ;
} /* clipAlignGetIntrons */

/*************************************************************************************/
/* get donors/acceptors: both functions fill both reads in case of pairs 
 * the results are available in pp->geneDonors geneAcceptors
 * which are soreted by fragment/read/position
 * and used in construct introns
 */ 
static int  clipAlignOptimizeFragmentHits (CLIPALIGN *pp)
{
  int nn = 0 ; 
    
  pp->geneDonors = bigArrayCreate (100, PEXPORT) ;
  pp->geneAcceptors = bigArrayCreate (100, PEXPORT) ;

  clipAlignGetDonors (pp) ;
  clipAlignGetAcceptors (pp) ;
  
  if (pp->splice)
    nn = clipAlignOptimizeGetIntrons (pp) ;
  clipAlignOptimizePairs (pp, TRUE) ; /* optimize one gene at a time */
  clipAlignMergeHits (pp) ;

  return nn ;
} /* clipAlignOptimizeFragmentHits */

/*************************************************************************************/
/* get donors/acceptors: both functions fill both reads in case of pairs 
 * the results are available in pp->geneDonors geneAcceptors
 * which are sorted by fragment/read/position
 * and used in construct introns
 */ 
static int clipAlignOptimizeGeneHits (CLIPALIGN *pp)
{
  int ii, jj = 0, nn = 0, fragment = 0, iMax ;
  MM *mm = 0 ;
  BigArray aa = pp->exportGeneHits ;
  PEXPORT *px, *py ;
  BOOL debug = FALSE ;

  clipAlignSquish (pp, 1) ;

  if (debug)
    fprintf (stderr, ".....clipAlignOptimizeGeneHits %s  received %ld hits :%s\n"
	     , pp->previousTargetGene ? dictName(pp->targetDict, pp->previousTargetGene) : dictName(pp->targetDict, pp->previousTarget) 
	     , bigArrayMax (pp->exportGeneHits)
	     , timeShowNow () 
	     ) ;
  /* process by fragment or by read if there are no pairs */ 
  aa = pp->exportGeneHits ; 
  bigArraySort (aa,  pExportPairXOrder) ;
  iMax = bigArrayMax (aa) ;
  for (ii = 0,  px = bigArrp (aa, ii, PEXPORT) ;
       ii <  iMax ; ii++, px++)
    {
      if (px->score)
	{
	  int jj ;
	  MM *mm = px->mm ;

	  int d0 = px->donor ;
	  int a0 = px->acceptor ;
	  
	  for (jj = ii + 1,  py = px + 1 ; jj <  iMax && py->mm == mm && py->x1 == px->x1 ; jj++, py++)
	    {
	      if (! py->score)
		continue ;
	      else if (py->a1 == px->a1)
		{
		  int d = py->donor ;
		  int a = py->acceptor ;

		  py->donor = d0 ; 
		  py->acceptor = a0 ;
		  if (! memcmp(px, py, sizeof (PEXPORT)))
		    py->score = 0 ;
		  else
		    {
		      py->donor = d ; 
		      py->acceptor = a ;
		    }
		}
	    }
	}
    }
  for (ii = 0,  px = bigArrp (aa, ii, PEXPORT) ;
       ii <  iMax ; ii++, px++)
    {
      if (px->score > 60)
	{
	  int jj ;
	  MM *mm = px->mm ;
	  
	  for (jj = ii + 1,  py = px + 1 ; jj <  iMax && py->mm == mm && py->x1 <= px->x2 ; jj++, py++)
	    {
	      if (py->score &&
		  py->x2 <= px->x2 && py->score < 30)
		py->score = 0 ;
	    }
	}
    }
  clipAlignSquish (pp, 1) ;

  aa = pp->exportGeneHits ; pp->exportGeneHits = 0 ;
  bigArraySort (aa,  pExportPairAOrder) ;
  iMax = arrayMax (aa) ;
  for (ii = jj = 0, px = bigArrp (aa, 0, PEXPORT) ; ii < iMax ; ii++, px++)
    {
      if (! px->score || ! px->mm)
	continue ;
      if (
	  (px->mm->fragment && fragment != px->mm->fragment) ||
	  ( ! px->mm->fragment && mm != px->mm)
	  )
	{
	  if (jj) 
	    nn += clipAlignOptimizeFragmentHits (pp) ;
	  jj = 0 ;
		    
	  bigArrayDestroy (pp->exportGeneHits) ;
	  pp->exportGeneHits = bigArrayCreate (100, PEXPORT) ;
	}
      py = bigArrayp (pp->exportGeneHits, jj++,  PEXPORT) ;
      *py = *px ;
      fragment = px->mm->fragment ; mm = px->mm ;
    }
  if (jj)
    nn += clipAlignOptimizeFragmentHits (pp) ;
  bigArrayDestroy (pp->exportGeneHits) ;
  bigArrayDestroy (aa) ;
  pp->exportGeneHits = bigArrayCreate (100, PEXPORT) ;

  if (debug)
    fprintf (stderr, ".....clipAlignOptimizeGeneHits %s  done kept %d hits:%s\n"
	     , pp->previousTargetGene ? dictName(pp->targetDict, pp->previousTargetGene) : dictName(pp->targetDict, pp->previousTarget) 
	     , nn 
	     , timeShowNow () 
	     ) ;
  return nn ;
} /* clipAlignOptimizeGeneHits */

/*************************************************************************************/
/* target by target, i.e. gene by gene or chromosome by chromosoe
 * align, optimize the exon/intron chains and register
 */
static int clipAlignSearch (CLIPALIGN *pp, Array dna)
{
  if (0)
    { /* a hack to only keep some chromosome in phase 1, all in phase 0 */
      if ( 
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000007.14") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000008.11") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000019.10") 
	   )
	;
      else if (1 && pp->hashPhase)
	return 0 ;
    }
 
  if (0)
    { /* a hack to only keep some chromosome in phase 1, all in phase 0 */
      if ( 
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000003.12") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000004.12") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000005.10") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000006.12") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000007.14") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000008.11") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000009.12") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000010.11") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000011.10") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000012.12") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000013.11") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000014.9") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000015.10") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000016.10") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000017.11") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000018.10") ||
	  ! strcmp (dictName (pp->targetDict, pp->target), "NC_000019.10") 
	   )
	;
      else if (1 && pp->hashPhase)
	return 0 ;
    }
 

  /* static Array previousDna = 0 ; */
  if (!dna ||
      (! pp->targetGene && pp->target != pp->previousTarget)  ||
      pp->targetGene != pp->previousTargetGene
      )
    {
      if (pp->oligoGene && pp->oligoGeneTouched && pp->nTargets > 1)
	{
	  bitSetReCreate (pp->oligoGene, 0) ; /* recreate at same size */
	  pp->oligoGeneTouched = FALSE ;
	}
      pp->nTargets = 0 ;

      if (bigArrayMax (pp->exportGeneHits) > 1 ||
	  (bigArrayMax (pp->exportGeneHits) == 1 && bigArr (pp->exportGeneHits, 0, PEXPORT).score)
	  )
	clipAlignOptimizeGeneHits (pp) ;
    }
  if (dna)
    clipAlignDoSearch (pp, dna) ;
  /* previousDna = dna ; */
  pp->previousTarget = pp->target ;
  pp->previousTargetGene = pp->targetGene ;
  pp->nTargets++ ;
  if (0 && pp->nTargets > 5 && !pp->oligoGene)
    pp->oligoGene = bitSetCreate (128000, pp->h) ;

  return bigArrayMax (pp->exportGeneHits) ;
} /* clipAlignSearch */

/*************************************************************************************/

static void clipAlignGetTm (CLIPALIGN *pp, int probeSequence, int probeLength, float *entropyP, float *tmP)
{
  unsigned char *probeSeq ;
  static Array dna = 0 ;

  if (! dna) dna = arrayCreate (100000, unsigned char) ;

  probeSeq = (unsigned char*) stackText (pp->probeStack, probeSequence) ;
  array (dna, probeLength, unsigned char) = 0 ; /* make room */
  memcpy (arrp(dna, 0, unsigned char), probeSeq, probeLength) ;
  arrayMax (dna) = probeLength ;
  
  *tmP = oligoTm (dna, 1, probeLength, 0) ;
  *entropyP = oligoEntropy (arrp(dna, 0, unsigned char), probeLength, 0) ;

  return ;
} /* clipAlignGetTm  */

/*************************************************************************************/

static void clipAlignExportTm (CLIPALIGN *pp)
{
  int jj ;
  float gc ;
  MM *mm ;
  Array probeDnaArray = 0 ;
  unsigned char *probeDna ; 

  if (!pp->aceOutGenomic)
    aceOutf (pp->ao, "#Probe\tLength bp\tEntropy  bp-equivalent\tTM\tGC\n") ;
  for (jj = arrayMax(pp->probes), mm = arrp (pp->probes, 0, MM) ; jj >0 ; jj--, mm++)
    {
      /* transform the unsigned char buffer into a dna Array */
      probeDna = (unsigned char*) stackText (pp->probeStack, mm->probeSequence) ;
      probeDnaArray = arrayReCreate (probeDnaArray, 256, unsigned char) ;
      array (probeDnaArray, mm->probeLength, unsigned char) = 0 ; /* make room */
      memcpy (arrp(probeDnaArray, 0, unsigned char), probeDna, mm->probeLength) ;
      arrayMax (probeDnaArray) = mm->probeLength ;
      mm->tm = oligoTm (probeDnaArray, 1, arrayMax(probeDnaArray), &gc) ;

      if (pp->aceOutGenomic)
	aceOutf (pp->ao, "Probe %s\nLength %d\nEntropy %d bp\nTm %.2f\tGC %.1f\n\n"
		 , freeprotect (stackText (pp->probeStack, mm->probeName))
		 , mm->probeLength, mm->entropy, mm->tm, gc
		 ) ;
      else
	aceOutf (pp->ao, "%s\t%d\t%d\t%g\t%.1f\n"
		 , stackText (pp->probeStack, mm->probeName)
		 , mm->probeLength, mm->entropy, mm->tm, gc
		 ) ;
    }

  arrayDestroy (probeDnaArray) ;
  return ;
} /* clipAlignExportTm  */

/*************************************************************************************/
/* Are there more words to prehash in the target or in the probes ? */
static BOOL clipAlignIsShortTarget (CLIPALIGN *pp)
{
  AC_HANDLE h ;
  BOOL ok = TRUE ;
  int nbp = 0, nseq = 0, line = 0 ;
  int isBig = 5000000 ; /* 10 Million creates a hash table of 500 Megabytes */
  BOOL two = TRUE ;
  char *cp ;
  ACEIN ai ;

  if (pp->multiTargetFileNames)
    return FALSE ;
  if (!pp->targetFileName)
    return TRUE ;
  if (!pp->probes)
    return TRUE ;

  if (pp->strand || pp->antistrand)
    two = FALSE ;
 
  if (arrayMax (pp->probes) < isBig)
    isBig = arrayMax (pp->probes) ;
  isBig = stackMark (pp->probeStack)/pp->seedShift ;

  h = ac_new_handle () ;

  ai = aceInCreate (pp->targetFileName, FALSE, h) ;
  aceInSpecial (ai, "\n") ;

  while (ok && aceInCard (ai)) /* will close pp->targetFile */
    {
      line++ ;
      cp = aceInWord (ai) ;
      if (cp) 
	switch (*cp)
	  {
	  case '/': break ;
	  case '#': break ;
	  case '>': nseq++ ; break ;
	  default: nbp += two ? (strlen (cp) << 1) : strlen (cp) ; break ;
	  }
      if (nbp > isBig) 
	ok = FALSE ; 
    }
  ac_free (h) ;
  fprintf (stderr, "//  clipAlignIsShortTarget = %s\n", ok ? "TRUE" : "FALSE") ;
  return ok ;
} /* clipAlignIsShortTarget */

/*************************************************************************************/

static int clipAlignRunOne (CLIPALIGN *pp)
{
  AC_HANDLE h ;
  int n, state, line = 0, nn = 0, seq = 0, nHits = 0 ;
  int lastTarget = 0, lastGene = 0 ;
  char *cp, *cq ;
  Stack s ;
  Array dna ;
  ACEIN ai ;

  if (!pp->targetFileName)
    return FALSE ;

  h = ac_new_handle () ;
  s = stackHandleCreate (300, h) ;
  dna = arrayHandleCreate (100000, unsigned char, h) ;

  ai = aceInCreate (pp->targetFileName, 0, h) ;
  if (!ai)
    messcrash ("cannot read target file %s", pp->targetFileName) ;
  aceInSpecial (ai, "\n") ;

  pp->scannedGenes = bitSetReCreate (pp->scannedGenes, 80000) ;
  state = nn = 0 ;
  while (aceInCard (ai)) /* will close pp->targetFile */
    {
      line++ ;
      cp = aceInWord (ai) ;  /* do NOT accept spaces in target names, this creates havock downstream everywhere */
      if (cp) switch (state)
		{
		case 1:
		case 2:
		  if (*cp != '>')
		    {
		      if (state == 1)
			{
			  state = 2 ;
			  pushText (s, cp) ;
			}
		      else
			catText (s, cp) ;
		      break ;
		    }
		  else
		    {
		      cq = stackText (s, seq) ;
		      n = strlen (cq) ;
		      array (dna , n, unsigned char) = 0 ; /* make room */
		      arrayMax (dna) = n ;
		      memcpy (arrp (dna, 0, unsigned char), cq, n) ;
		      dnaEncodeArray (dna) ;
		      if (pp->A2G)
			{
			  register int i = n ;
			  register char *cr = arrp (dna, 0, char) - 1 ;
			  while (cr++, i--)
			    if (*cr == A_) *cr = R_ ;
			}
		      if (pp->G2A)
			{
			  register int i = n ;
			  register char *cr = arrp (dna, 0, char) - 1 ;
			  while (cr++, i--)
			    if (*cr == G_) *cr = R_ ;
			}
		      if (pp->T2C)
			{
			  register int i = n ;
			  register char *cr = arrp (dna, 0, char) - 1 ;
			  while (cr++, i--)
			    if (*cr == T_) *cr = Y_ ;
			}
		      if (pp->C2T)
			{
			  register int i = n ;
			  register char *cr = arrp (dna, 0, char) - 1 ;
			  while (cr++, i--)
			    if (*cr == C_) *cr = Y_ ;
			}
		      if (pp->A2C)
			{
			  register int i = n ;
			  register char *cr = arrp (dna, 0, char) - 1 ;
			  while (cr++, i--)
			    if (*cr == A_) *cr = M_ ;
			}
		      if (pp->C2A)
			{
			  register int i = n ;
			  register char *cr = arrp (dna, 0, char) - 1 ;
			  while (cr++, i--)
			    if (*cr == C_) *cr = M_ ;
			}
		      if (pp->T2G)
			{
			  register int i = n ;
			  register char *cr = arrp (dna, 0, char) - 1 ;
			  while (cr++, i--)
			    if (*cr == T_) *cr = K_ ;
			}
		      if (pp->G2T)
			{
			  register int i = n ;
			  register char *cr = arrp (dna, 0, char) - 1 ;
			  while (cr++, i--)
			    if (*cr == G_) *cr = K_ ;
			}
		      if (pp->G2C)
			{
			  register int i = n ;
			  register char *cr = arrp (dna, 0, char) - 1 ;
			  while (cr++, i--)
			    if (*cr == G_) *cr = S_ ;
			}
		      if (pp->C2G)
			{
			  register int i = n ;
			  register char *cr = arrp (dna, 0, char) - 1 ;
			  while (cr++, i--)
			    if (*cr == C_) *cr = S_ ;
			}
		      if (pp->A2T)
			{
			  register int i = n ;
			  register char *cr = arrp (dna, 0, char) - 1 ;
			  while (cr++, i--)
			    if (*cr == A_) *cr = W_ ;
			}
		      if (pp->T2A)
			{
			  register int i = n ;
			  register char *cr = arrp (dna, 0, char) - 1 ;
			  while (cr++, i--)
			    if (*cr == T_) *cr = W_ ;
			}

		      if (pp->mutantMalus)
			pp->targetMutantMalus =  strncmp (stackText (s, 0), "M:", 2) ? 0 : pp->mutantMalus ;
		      pp->targetGene = 0 ;
		      cq = (char *)strcasestr(stackText (s, 0), "|Gene|") ;
		      if (cq) 
			{
			  *cq=0 ;
			  cq += strlen( "|GENE|") ;
			  if (*cq)
			    {
			      char *cr ;
			      cr = strstr (cq, "|") ;
			      if (cr) *cr = 0 ;
			      dictAdd (pp->targetDict, cq, &(pp->targetGene)) ;
			      if (! strstr (stackText (s, 0), "unspliced"))
				bitSet (pp->isSplicedTarget, pp->targetGene) ;
			      if (strstr (cq, "--") || strstr (cq, "-AND-") || strstr (cq, "-and-"))
				bitSet (pp->doubleGenes, pp->targetGene) ;
			      if (cr) *cr = '|' ;
			    }
			}
		      dictAdd (pp->targetDict,  stackText (s, 0), &(pp->target)) ;

		      if (! pp->targetGene && pp->target2geneArray &&
			  pp->target < arrayMax (pp->target2geneArray)
			  )
			{
			  pp->targetGene = arr (pp->target2geneArray, pp->target, PEXPORT).targetGene ;
			}
		      if (! pp->targetGene && pp->target2geneArray)
			{
			  cq = stackText (s, 0) ;
			  if (!strncmp (cq, "MRNA:", 4))
			    {
			      int n1 = 0 ;
			      if (dictFind (pp->targetDict, cq + 4, &n1) &&
				  n1 < arrayMax (pp->target2geneArray)
				  )
				pp->targetGene = arr (pp->target2geneArray, n1, PEXPORT).targetGene ;
			    }
			}

		      if (! pp->hashPhase && pp->targetGene && bitt (pp->scannedGenes, pp->targetGene) && ! pp->splice && pp->MRNAH && lastGene && 
			  pp->targetGene != lastGene)
			messcrash ("Target fasta file sorting error\nWhen the target fasta files refers to alternative transcripts denoted >name|GENE|geneName, the fasta file must be sorted by gene, then by transcript so that the aligner can correctly assign reads and pairs to the best transcript of a given gene\n--- Line %d of file %s\n--- In gene %s,\n--- transcript %s|GENE|%s\n--- wrongly appeared below %s|GENE|%s\nPlease sort the target fasta file in alphabetic order of the transcripts names."
				   , aceInStreamLine (ai), pp->targetFileName  /* aceInFileName is now null, since ai is closed, but the line number is ok */
				   , dictName (pp->targetDict, pp->targetGene)
				   , dictName (pp->targetDict,  pp->targetGene), dictName (pp->targetDict,  pp->target)
				   , lastGene ? dictName (pp->targetDict,  lastGene) : "-"
				   , lastTarget ? dictName (pp->targetDict,  lastTarget) : "-"
				   ) ;
		      bitSet (pp->scannedGenes, pp->targetGene) ;
		      if (! pp->hashPhase && pp->MRNAH && lastTarget && lastGene && 
			  pp->targetGene   && ! pp->splice && pp->MRNAH  &&
			  pp->targetGene == lastGene &&
			  lexstrcmp (dictName (pp->targetDict,  lastTarget), dictName (pp->targetDict,  pp->target)) > 0
			  )
			messcrash ("Target fasta file sorting error\nOption -MRNAH requires that, for a given gene, the target fasta file be sorted in alphabetic order of the transcripts, since this affects the hierarchic mapping of the reads to the alternative transcripts.\n--- Line %d of file %s\n--- In gene: %s,\n--- transcript %s|GENE|%s\n--- wrongly appeared below %s|GENE|%s\nPlease sort the target fasta file in alphabetic order of the transcripts names."
				   , aceInStreamLine (ai), pp->targetFileName  /* aceInFileName is now null, since ai is closed, but the line number is ok */
				   , dictName (pp->targetDict, pp->targetGene)
				   , dictName (pp->targetDict,  pp->target)
				   , dictName (pp->targetDict,  pp->targetGene)
				   , lastTarget ? dictName (pp->targetDict,  lastTarget) : "-"
				   , lastGene ? dictName (pp->targetDict,  lastGene) : "-"
				   ) ;
		      lastGene = pp->targetGene ;
		      lastTarget = pp->target ;

		      nHits += clipAlignSearch (pp, dna) ;
		      state = 0 ; /* and fall thru */
		      stackClear (s) ;
		      dna = arrayReCreate (dna,100000, unsigned char) ; 		
		    }
		  /* fall thru to new sequence */
		case 0: /* expecting   >target_name */
		  if (*cp == '#')
		    break ;
		  if (*cp != '>' || ! *(cp+1))
		    {
		      fprintf (stderr, "// bad character in file %s line %d, expecting a '>', got %s"
			       , pp->targetFileName, line, cp) ;
		      return FALSE ;
		    }
		  pushText (s, cp + 1) ;
		  state = 1 ;
		  seq = stackMark (s) ;
		  break ;
		}
    }
  if (state > 0)
    {
      cq = stackText (s, seq) ;
      n = strlen (cq) ;
      array (dna , n, unsigned char) = 0 ; /* make room */
      arrayMax (dna) = n ;
      memcpy (arrp (dna, 0, unsigned char), cq, n) ;
      dnaEncodeArray (dna) ;
      if (pp->A2G)
	{
	  register int i = n ;
	  register char *cr = arrp (dna, 0, char) - 1 ;
	  while (cr++, i--)
	    if (*cr == A_) *cr = R_ ;
	}
      if (pp->G2A)
	{
	  register int i = n ;
	  register char *cr = arrp (dna, 0, char) - 1 ;
	  while (cr++, i--)
	    if (*cr == G_) *cr = R_ ;
	}
      if (pp->T2C)
	{
	  register int i = n ;
	  register char *cr = arrp (dna, 0, char) - 1 ;
	  while (cr++, i--)
	    if (*cr == T_) *cr = Y_ ;
	}
      if (pp->C2T)
	{
	  register int i = n ;
	  register char *cr = arrp (dna, 0, char) - 1 ;
	  while (cr++, i--)
	    if (*cr == C_) *cr = Y_ ;
	}
      if (pp->A2C)
	{
	  register int i = n ;
	  register char *cr = arrp (dna, 0, char) - 1 ;
	  B2[M_] = B2[A_] ;
	  B2[K_] = B2[T_] ;
	  while (cr++, i--)
	    if (*cr == A_) *cr = M_ ;
	}
      if (pp->C2A)
	{
	  register int i = n ;
	  register char *cr = arrp (dna, 0, char) - 1 ;
	  B2[M_] = B2[C_] ;
	  while (cr++, i--)
	    if (*cr == C_) *cr = M_ ;
	}
      if (pp->T2G)
	{
	  register int i = n ;
	  register char *cr = arrp (dna, 0, char) - 1 ;
	  while (cr++, i--)
	    if (*cr == T_) *cr = K_ ;
	}
      if (pp->G2T)
	{
	  register int i = n ;
	  register char *cr = arrp (dna, 0, char) - 1 ;
	  B2[K_] = B2[G_] ;
	  while (cr++, i--)
	    if (*cr == G_) *cr = K_ ;
	}
      if (pp->G2C)
	{
	  register int i = n ;
	  register char *cr = arrp (dna, 0, char) - 1 ;
	  B2[S_] = B2[G_] ;
	  while (cr++, i--)
	    if (*cr == G_) *cr = S_ ;
	}
      if (pp->C2G)
	{
	  register int i = n ;
	  register char *cr = arrp (dna, 0, char) - 1 ;
	  B2[S_] = B2[C_] ;
	  while (cr++, i--)
	    if (*cr == C_) *cr = S_ ;
	}
      if (pp->A2T)
	{
	  register int i = n ;
	  register char *cr = arrp (dna, 0, char) - 1 ;
	  B2[W_] = B2[A_] ;
	  while (cr++, i--)
	    if (*cr == A_) *cr = W_ ;
	}
      if (pp->T2A)
	{
	  register int i = n ;
	  register char *cr = arrp (dna, 0, char) - 1 ;
	  B2[W_] = B2[T_] ;
	  while (cr++, i--)
	    if (*cr == T_) *cr = W_ ;
	}
      if (! cq) cq = strstr(stackText (s, 0), "|GENE|") ;
      if (pp->mutantMalus)
	pp->targetMutantMalus =  strncmp (stackText (s, 0), "M:", 2) ? 0 : pp->mutantMalus ;
      pp->targetGene = 0 ;
      cq = (char *)strcasestr(stackText (s, 0), "|Gene|") ;
      if (cq) 
	{
	  *cq=0 ;
	  cq += strlen( "|GENE|") ;
	  if (*cq)
	    {
	      char *cr ;
	      cr = strstr (cq, "|") ;
	      if (cr) *cr = 0 ;
	      pp->hasGenes = TRUE ;
	      dictAdd (pp->targetDict, cq, &(pp->targetGene)) ;
	      if (cr) *cr = '|' ;
	    }
	}
      dictAdd (pp->targetDict,  stackText (s, 0), &(pp->target)) ;
      if (! pp->targetGene && pp->target2geneArray &&
	  pp->target < arrayMax (pp->target2geneArray)
	  )
	{
	  pp->targetGene = arr (pp->target2geneArray, pp->target, PEXPORT).targetGene ;
	}
      if (! pp->targetGene && pp->target2geneArray)
	{
	  cq = stackText (s, 0) ;
	  if (!strncmp (cq, "MRNA:", 4))
	    {
	      int n1 = 0 ;
	      if (dictFind (pp->targetDict, cq + 4, &n1) &&
		  n1 < arrayMax (pp->target2geneArray)
		  )
		pp->targetGene = arr (pp->target2geneArray, n1, PEXPORT).targetGene ;
	    }
	}
      if (! pp->hashPhase && pp->targetGene && bitt (pp->scannedGenes, pp->targetGene) && ! pp->splice && pp->MRNAH && 
	  pp->targetGene != lastGene)
	messcrash ("Target fasta file sorting error\nWhen the target fasta files refers to alternative transcripts denoted >name|GENE|geneName, the fasta file must be sorted by gene, then by transcript so that the aligner can correctly assign reads and pairs to the best transcript of a given gene\n--- Line %d of file %s\n--- In gene %s,\n--- transcript %s|GENE|%s\n--- wrongly appeared below %s|GENE|%s\nPlease sort the target fasta file in alphabetic order of the transcripts names."
		   , aceInStreamLine (ai), pp->targetFileName  /* aceInFileName is now null, since ai is closed, but the line number is ok */
		   , dictName (pp->targetDict, pp->targetGene)
		   , dictName (pp->targetDict,  pp->targetGene), dictName (pp->targetDict,  pp->target)
		   , dictName (pp->targetDict,  lastGene), dictName (pp->targetDict,  lastTarget)
		   ) ;
      bitSet (pp->scannedGenes, pp->targetGene) ;
      if (! pp->hashPhase && pp->MRNAH && lastTarget  && 
	  pp->targetGene  && ! pp->splice && pp->MRNAH &&
	  pp->targetGene == lastGene &&
	  lexstrcmp (dictName (pp->targetDict,  lastTarget), dictName (pp->targetDict,  pp->target)) > 0
	  )
	messcrash ("Target fasta file sorting error\nOption -MRNAH requires that, for a given gene, the target fasta file be sorted in alphabetic order of the transcripts, since this affects the hierarchic mapping of the reads to the alternative transcripts.\n--- Line %d of file %s\n--- In gene: %s,\n--- transcript %s|GENE|%s\n--- wrongly appeared below %s|GENE|%s\nPlease sort the target fasta file in alphabetic order of the transcripts names."
		   , aceInStreamLine (ai), pp->targetFileName  /* aceInFileName is now null, since ai is closed, but the line number is ok */
		   , dictName (pp->targetDict, pp->targetGene)
		   , dictName (pp->targetDict,  pp->targetGene), dictName (pp->targetDict,  pp->target)
		   , dictName (pp->targetDict,  lastGene), dictName (pp->targetDict,  lastTarget)
		   ) ;
      lastGene = pp->targetGene ;
      lastTarget = pp->target ;

      bitUnSet (pp->isSplicedTarget, dictMax (pp->targetDict) + 1) ; /* make room */

      nHits += clipAlignSearch (pp, dna) ;
      state = 0 ; /* and fall thru */
      stackClear (s) ;
    }


  bitUnSet (pp->isSplicedTarget, dictMax (pp->targetDict) + 1) ;             /* make room */
  if (bigArrayMax (pp->exportGeneHits) > 1 ||
      ( bigArrayMax (pp->exportGeneHits) == 1 && bigArr (pp->exportGeneHits, 0, PEXPORT).score)
      )
    nHits += clipAlignSearch (pp, 0) ;            /* merge the results of the last target */
  ac_free (h) ;
  return nHits ;
} /* clipAlignRunOne */

/*************************************************************************************/
/* loop on multiTargets */
static int clipAlignRun (CLIPALIGN *pp)
{
  int nHits = 0 ;

  if (! pp->multiTargetFileNames)
    return clipAlignRunOne (pp) ;
  else
    {
      AC_HANDLE h = ac_new_handle () ;
      char *cp = strnew (pp->multiTargetFileNames, h) ;
      char *cq, *cr ;
      
      while (cp)
	{
	  int n = 0 ;
	  cq = strchr (cp, ',') ;
	  if (cq) 
	    *cq++ = 0 ;
	  else
	    cq = 0 ;
	  cr = strchr (cp, ':') ;
	  if (cr)
	    { 
	      n++ ; *cr++ = 0 ; pp->target_class = cp ; 
	    }
	  if (n == 1)
	    {
	      char *cp1 = cr ;
	      cr = strchr (cr, ':') ;
	      *cr++ = 0 ; 
	      if (sscanf (cp1, "%d", &pp->bonus) == 1 &&
		  pp->bonus >= 0 &&
		  pp->bonus < 40
		  )
		n++ ;
	      else
		messcrash ("cannot scan the bonus in %s", pp->multiTargetFileNames) ;
	    }
	  if (n == 2 && cr && *cr)
	    {
	      n++ ;
	      pp->targetFileName = cr ;
	    }
	  if (n == 3)
	    {
	      if (0)
		fprintf (stderr, "// ... %s : %s\n", timeShowNow (), pp->targetFileName) ;
	      nHits += clipAlignRunOne (pp) ;
	    }
	  else
	    break ;
	  cp = cq ;
	}
      ac_free (h) ;
    }
  pp->target_class = "multi" ;
  pp->targetFileName = pp->multiTargetFileNames ;

  return nHits ;
} /* clipAlignRun */

/*************************************************************************************/
/*************************************************************************************/

static int clipAlignPreScanGenome (CLIPALIGN *pp)
{
  int nn = 0 ;
  AC_HANDLE h = ac_new_handle () ; 
  
  ac_free (h) ;
  return nn ;
} /* clipAlignPreScanGenome */

/*************************************************************************************/

static int clipAlignPreScanProbes (CLIPALIGN *pp)
{
  int nn = 0 ;
  AC_HANDLE h = ac_new_handle () ; 
  
  ac_free (h) ;
  return nn ;
} /* clipAlignPreScanProbes */

/*************************************************************************************/
/*************************************************************************************/

static int clipAlignCountSequenceWithBadLetters (CLIPALIGN *pp)
{
  int ii, j ;
  char *cp ;
  MM *mm ;
  BOOL ok ;
  int nSeqOK = 0, nSeqBad = 0 ;

  if (arrayMax (pp->probes))
    for (ii = 0 ; ii < arrayMax (pp->probes) ; ii++)
      {
	mm = arrp (pp->probes, ii, MM) ;
	/* load the first mask */
	cp = stackText (pp->probeStack, mm->probeSequence) ;
	if (!cp)
	  continue ;
	cp-- ; ok = TRUE ; j = 0 ;
	while (j++, *++cp)
	  switch ((int)*cp)
	    {
	    case A_: case T_: case G_: case C_: break ;
	    default: 
	      ok = FALSE ;
	      printf ("%c at %5d in %s\n"
		      , dnaDecodeChar[(int)*cp], j 
		      , stackText (pp->probeStack, mm->probeName)
		      ) ;
	      break ;
	    }
	if (ok) nSeqOK++ ;
	else  nSeqBad++ ;
      }
  printf ("Found %d sequences with just ATGC, %d sequence with some other letters\n",
	  nSeqOK, nSeqBad) ;

  if (0) clipAlignCountSequenceWithBadLetters (0) ; /* for compiler happiness */
  return  nSeqBad ;
}

/*************************************************************************************/

static BOOL oligoCheckEntropy (int limit, int na, int nt, int ng, int nc)
{
  int ss = 0 ;
  int nn = na + nt + ng + nc ;
  static BOOL oldNn = -1 ;
  static Array ee = 0 ;

  /* lazy calculation of all the logs */
  if (! ee || (nn > 1 && nn != oldNn))
    {
      double s = 0, log4 = log(4.0) ; /* divide by log4 to be in base equivalent */ 
      int j ;

      ee = arrayReCreate (ee, nn, int) ;
      for (j = 1 ; j <= nn ; j++)
        { s = j ; s /= nn ; array (ee, j, int) = - (int) (1000.0 * j * log(s)/log4 ) ; }
      array (ee, 0, int) = 0 ;
    }

  ss = arr (ee, na, int) + arr (ee, nt, int) + arr (ee, ng, int) + arr (ee, nc, int) ;
  ss = (ss + 499)/1000 ;
  return ss > limit ? TRUE : FALSE ;
} /* oligoTm */


/*************************************************************************************/

static void clipAlignHashTarget (CLIPALIGN *pp)
{
  /* with size >= 1000 we are sure to create the accelerating bitsets in the associators */
  ac_free (pp->ass) ;
  pp->ass = assBigCreate (1000 + arrayMax (pp->probes)) ; 
  ac_free (pp->ass1) ;
  pp->ass1 = assBigCreate (1000 + arrayMax (pp->probes)) ; 
  pp->nOligo = 0 ;
  pp->oligoFrequency = arrayCreate (1000 + arrayMax (pp->probes), unsigned char) ;
  if (pp->doubleSeed)
    {
      pp->assB1 = assBigCreate (1000 + arrayMax (pp->probes)) ; 
      pp->oligoFrequencyB = arrayCreate (1000 + arrayMax (pp->probes), unsigned char) ;
    }
  nHashRejected = 0 ;
  
  clipAlignRun (pp) ;
} /* clipAlignHashTarget */

/*************************************************************************************/

static int clipAlignHashProbes (CLIPALIGN *pp)
{
  int N25 = 25 ;
  int nf, i, jx, jxB, ii, n1 = 0, gap, ok, ok2, nw, phase ;
  int na, nt, ng, nc ;
  int nFoundIn1 = 0 ;
  int delta = pp->seedShift ? pp->seedShift :  pp->probeLengthMax ;
  unsigned long int oligo, oligoB, mask, maskB ;
  char cc, *cp, *cp0 ;
  BOOL doubleSeed = pp->doubleSeed ;
  MM *mm ;
  BOOL goodName, goodNameB ;
  BOOL is25 = FALSE, debug = FALSE ;
  static int shiftLoop = 0 ;
  int hashPhase = pp->hashPhase ;
  int nOligo, nLoop = 0, nSkipped = 0, nProbes = 0, nWords = 0 ;
  int oldFirstGoodProbe = 0, firstGoodProbe = 0  ;
  Array oligoFrequency ;
  Array oligoFrequencyB ;
  Associator ass, ass1, ass2, assB, assB1, assB2 ;

  adaptorMask (pp) ; 
  if (hashPhase == 0)
    {
      /* with size >= 1000 we are sure to create the accelerating bitsets in the associators */
      ac_free (pp->ass) ;
      pp->ass = assBigHandleCreate (1000 + 4 * arrayMax (pp->probes), pp->h) ; 
      ac_free (pp->ass1) ;
      pp->ass1 = assBigHandleCreate (1000 + 4 * arrayMax (pp->probes), pp->h) ; 

      if (pp->doubleSeed)
	{
	  ac_free (pp->assB) ;
	  pp->assB = assBigHandleCreate (1000 + arrayMax (pp->probes), pp->h) ; 
	  ac_free (pp->assB1) ;
	  pp->assB1 = assBigHandleCreate (1000 + arrayMax (pp->probes), pp->h) ; 
	  pp->oligoFrequencyB = arrayHandleCreate (1000 + arrayMax (pp->probes), unsigned char, pp->h) ;
	}
      pp->nOligo = 0 ;
      pp->oligoFrequency = arrayHandleCreate (1000 + arrayMax (pp->probes), unsigned char, pp->h) ;
      nHashRejected = 0 ;
    } 
  if (hashPhase == 1)
    {
      ac_free (pp->ass2) ;
      pp->ass2 = assInOnlyHandleCreate (1000 + pp->ass1->n, pp->h) ; 
      if (pp->doubleSeed)
	pp->assB2 = assInOnlyHandleCreate (1000 + pp->assB1->n, pp->h) ; 
    }
  shiftLoop++ ;

  ass = pp->ass ;
  ass1 = pp->ass1 ;
  ass2 = pp->ass2 ;
  assB = pp->assB ;
  assB1 = pp->assB1 ;
  assB2 = pp->assB2 ;
  oligoFrequency = pp->oligoFrequency ;
  oligoFrequencyB = pp->oligoFrequencyB ;
  
  gap = pp->gap ;
  jx = pp->dx2 ;
  mask = pp->mask ;
  /*
    mask1 = pp->mask1 ;
    mask2 = pp->mask2 ;
  */
  pp->seedOffsetUsed = gap + 1 ;
  pp->dx1 = 0 ;
  jxB = is64 ? 32 : 16 ;
  maskB = is64 ? 0xffffffffffffffff  : 0xffffffff ;

  /* if nOk == 0 no probe is long enough */
  for (ii = 0 ; ii < arrayMax (pp->probes) ; ii++)
    {
      int x00 = 0, deltaHack, deltaGap = 0 ; 
      mm = arrp (pp->probes, ii, MM) ;
      oldFirstGoodProbe = firstGoodProbe ;
      firstGoodProbe = 0 ;
      if (pp->gap + jx > mm->probeLength)
	continue ;
      nProbes++ ; nw = 0 ;
      is25 = FALSE ;
      ok = 0 ;
      if (1 && delta) deltaGap = ii % delta ; 
      for (nLoop = 0, gap = pp->gap + deltaGap, x00 = deltaHack = 0 ;  ; gap += delta + deltaHack, nLoop++)
	{
	  nw = 0 ;
	  if (nLoop == 0)  /* last base of first probe, we hope it is OUT of one of the other probes */
	    x00 = gap + jx - 1 ;
	  if (!pp->seedShift && nLoop > 0)
	    break ;
	  if (1 && gap + jx > mm->probeLength)
	    {
	      if (pp->seedShift && x00 > 0)  /* if needed add an extra probe just behind the firts probe choice */
		{ gap = x00 + 2 ; deltaHack = 2 * pp->probeLengthMax ; x00 = 0 ; } /* try an extra probe, cut the next loop */ 
	    }
	  if (gap + jx > mm->probeLength)
	    break ;

	  if (ii && oldFirstGoodProbe && gap + jx >= oldFirstGoodProbe) 
	    { ok = 1 ; firstGoodProbe = oldFirstGoodProbe ; }
	  if (ii && (mm-1)->same >= gap + jx)
	    { nSkipped++ ; nw++ ; continue ; }
	      
	      
	  /* 2010_09_09
	   * if we have not found a word, try to shift by just 2 bp 
	   * but this does not really works, because words shifted by 2 bp
	   * have not been entered in the ass1, so we gain nothing
	   * we disable that code (! phase in the loop control)
	   */
	  for (phase = 0 ; !phase  && ! nw && phase < delta ; phase++)
	    {
	      if (phase > 0 && (!pp->seedShift || hashPhase == 0))
		break ;
	      goodName = TRUE ;  goodNameB = TRUE ; 
	      if (gap + phase + jx > mm->probeLength)
		break ;
	      if (gap + phase + jxB > mm->probeLength)
		goodNameB = FALSE ;
	      
	      /* load the first mask */
	      cp0 = cp = stackText (pp->probeStack, mm->probeSequence) + gap + phase ;
	      oligo = 0 ; oligoB = 0 ; i = 0 ;
	      na = nt = ng = nc = 0 ;
	      for (i = 0 ;  i < (doubleSeed ? jxB : jx) ; i++)
		{
		  cp = cp0 + i ;
		  cc = *cp & 0xf ;
		  if (i<jx)
		    switch (cc)
		      {
		      case A_: na++ ; break ; 
		      case T_: nt++ ; break ; 
		      case G_: ng++ ; break ; 
		      case C_: nc++ ; break ; 
		      }
		  switch (cc)
		    {
		    case A_:
		    case T_:
		    case G_:
		    case C_:
		      if (i<jx) { oligo <<= 2 ; oligo |= B2[(int)cc] ; }
		      oligoB <<= 2 ; oligoB |= B2[(int)cc] ;
		      break ;
		    default:
		      if (i < jx) goodName = FALSE ;
		      goodNameB = FALSE ;
		      break ;
		    }
		}
	      oligo &= mask ;
	      oligoB &= maskB ;
	      if (debug)
		debugOligo (messprintf ("...                                          probe %d, loop=%d phase=%d gap=%d",ii,nLoop,phase,gap), oligo) ;
	      
	      if (goodName && oligo && oligo != (0xffffffffffffffff & mask) 
		  && oligoCheckEntropy (jx/2, na,nt,ng,nc)
		  && (! pp->adaptorAss || ! assFind (pp->adaptorAss, assVoid(oligo), 0))
		  ) /* oligo is zero if AAAAAAAAAAAA */
		{
		  void *vp ;

		  if (debug)
		    debugOligo (messprintf ("...                               entropy ok probe %d, loop=%d phase=%d gap=%d",ii,nLoop,phase,gap), oligo) ;
		  switch (hashPhase)
		    {
		    case 0:
		      /* a first run just to count the frequency of the actual oligo in the genome */
		      nw++ ;
		      if (assInsert (ass1, assVoid(oligo), assVoid(pp->nOligo)))
			{
			  nWords++ ;
			  pp->nOligo++ ;
			  if (debug) fprintf (stderr, " ass1 ZZ inserted nOligo=%d gap=%d  oligo=%lx\n", pp->nOligo, gap, oligo) ; 
			}
		      if (doubleSeed && goodNameB && assInsert (assB1, assVoid(oligoB), assVoid(pp->nOligo)))
			{
			  nWords++ ;
			  pp->nOligo++ ;
			  if (debug) fprintf (stderr, " assB1 inserted nOligo=%d gap=%d  oligo=%lx\n", pp->nOligo, gap, oligoB) ; 
			}
		      break ;
		    case 1:
		      ok2 = 0 ;
		      if (assFind (ass1, assVoid(oligo), &vp))
			{
			  nFoundIn1++ ;
			  nOligo = assInt (vp) ;
			  if (debug )
			    fprintf(stderr, "....found in 1 nOligo=%d ii=%d\n", nOligo, ii) ;
			  if (nOligo < arrayMax (oligoFrequency) &&
			      (nf = arr (oligoFrequency, nOligo, unsigned char)) &&
			      nf >= 0)
			    {
			      if (nf < N25 || nf <= 2 * pp->maxHit || pp->maxHit > N25)
				{
				  nWords++ ; nw++ ;
				  assMultipleInsert (ass, assVoid(oligo), assVoid(ii+1+ (((unsigned long int)gap+phase+1) << (is64 ? 32 : 22)))) ;
				  n1++ ;
				  assInsert (ass2, assVoid(oligo), assVoid(1)) ;
				  ok = ok2 = 1 ;
				  if (! firstGoodProbe)
				    firstGoodProbe = gap ;
				}
			      else 
				is25 = TRUE ;
			    } 
			}
		      if (ok2 == 0 && doubleSeed && assFind (assB1, assVoid(oligoB), &vp))
			{
			  nOligo = assInt (vp) ;
			  if (nOligo < arrayMax (oligoFrequencyB) &&
			      (nf = arr (oligoFrequencyB, nOligo, unsigned char)) &&
			      nf >= 0 && (nf < 2 *  pp->maxHit || pp->maxHit > N25))
			    {
			      nWords++ ; nw++ ;
			      assMultipleInsert (assB, assVoid(oligoB), assVoid(ii+1+ (((unsigned long int)gap+phase+1) << (is64 ? 32 : 22)))) ;
			      n1++ ;
			      assInsert (assB2, assVoid(oligoB), assVoid(1)) ;
			      ok = ok2 = 1 ;
			    }
			}
		      if (ok2 == 0)
			nHashRejected++ ;
		      break ;
		    }
		}
	    }
	}
      if (! ok && is25)
	{
	  if (! pp->seedOver25)
	    pp->seedOver25 = bitSetCreate (arrayMax (pp->probes), pp->h) ;
	  bitSet (pp->seedOver25, ii) ;
	}
    }
  /* aceOutf (pp->ao, "#################// %s Created an associator for %d Probe\n", timeShowNow(), n1) ; */
  fprintf (stderr, "// hashPhase=%d nProbes = %d nWords=%d nFoundIn1=%d nSkipped = %d nRejected=%d\n"
	   , hashPhase, nProbes, nWords, nFoundIn1, nSkipped, nHashRejected) ;

  return nProbes ;
} /* clipAlignHashProbes */

/*************************************************************************************/
/*************************************************************************************/
/* method of Oleg, circular hash of the error  */
#ifdef JUNKNOTUSED
static int clipAlignHashOLeg (CLIPALIGN *pp)
{
  /*
    sept 18, 2007
    Oleg explained a method for hashing all the words with a circular rotation of the hasher
    if we analyse 1M probes of length 50, this means
    50M * (4 bytes for h2 + 4bytes for address) * 2 = 800Mb 
    last factor 2 is to get some empty space in the hasher

    solexa = 5M times 35bp
    we could move by 4 bp at a time
    we could identify each pairs and code on 1 bit and look for the 3 hashings of 16bp

    h = (h<<<7)^Tab[c]
    h = (h<<<7)^T1[c1]^T2[c2] 
    c1 base going out c2 going in
    t1[a T G C] each is a 32 bit word initialised to sqrt(2 3 5 10)
    t2[x] = t1[x] << ((7*seq length) modulo 32)

    then hash in a simple hash, no bucket no double hashing
    drop the probes that bounce and treat them again later
    insert in the hash not the name of the sequence but a second hash value
    so it is more compact and sufficient to be sure of quasi identity

    scan genome,compute hash value h1, go in hash, read h2, compare to hash2(genome)
    if correct export coordinates
    suppose hash1 is 10 bits, corresponds to size 2^10 table but h2 has 32 bits
    so now the reconnaissance pattern is 10+32=42 bit so is certain and
    we only start careful analysys on certitude of correct seed

  */
  /*
    alternative idea
    store 1M solexa = 35M + 100Mb genome 
    use the search repeats system (uses 4 bytes per base) and get all results
    this uses 500M of memory
    stop at words of minimal length 15 = exact seed
    or get best match immediately
    but we may also recognize series of shorter words with same target ?
    parallelize 30 times on the farm

    Dear Oleg
    i like very much the system of yesterday
    but we may have overlooked on point
    sequences are very repeated, so many
    will have identical  streches
    so ithink we cannot drop subsequences that hit an
    occupied spot of the hash table and wait, there will 
    be too many

    but we can do the following
    first i hash the sequences using h1 and i count the occurence
    so
    counter[hash1(subseq)] is the number of times h1 is used
    which means probably this sequence is not informative

    then i run again using a double hasher, but while i scan the
    subsequences, i first look at counter[] and if this number is
    over a threhold, i do NOT insert that subsequence in my
    double hasher
    this way i will have few bounces

    this is equivalent to your proposal if threshold is counter[]<=1
    but if i use counter[] <= 5, i may be more efficient

    other idea is to say that we scan a solexa sequence of length 45bases
    this creates 30 subword of length 16
    we use your system and count how many subword can be inserted in the
    hash,
    probes for which 16 words or more are used are considered succesful
    the others are reprocessed later 

    yours
    jean

  */
  return 1 ;
} /* clipAlignHashOLeg */
#endif

/*************************************************************************************/
/****************************** utlities *********************************************/

static BOOL clipAlignOpenFiles (CLIPALIGN *pp)
{
  if (pp->probeFileName)
    {
    }
  if (pp->targetFileName)
    {
    }
  return TRUE ;
} /* clipAlignOpenFiles */

/*************************************************************************************/

static void target2GeneParse (CLIPALIGN *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *ccp ;
  ACEIN ai = 0 ;
  int nn = 0, nd = 0, n1, n2, x1;
  PEXPORT *px ;

  ai = aceInCreate (pp->target2geneFileName, 0, h) ;
  if (!ai)
    messcrash ("cannot read target2gene file %s", pp->target2geneFileName) ;
  aceInSpecial (ai, "\n") ;
  if (!  pp->targetDict)
    pp->targetDict = dictHandleCreate (10000, pp->h) ;
  pp->target2geneArray = arrayHandleCreate (10000, PEXPORT, pp->h) ;

  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (! ccp || ! *ccp) continue ;
      dictAdd (pp->targetDict, ccp, &n1) ;
      aceInStep (ai,'\t') ;

      ccp = aceInWord (ai) ;
      if (! ccp || ! *ccp) continue ;
      dictAdd (pp->targetDict, ccp, &n2) ;

      if (strstr (ccp, "--") || strstr (ccp, "-AND-") || strstr (ccp, "-and-"))
	bitSet (pp->doubleGenes, n2) ;

      px = arrayp (pp->target2geneArray, n1, PEXPORT) ;
      if (px->targetGene && px->targetGene != n2)
	{
	  nd++ ; 
	  if (0) messcrash ("Target 2 gene file contains contradictory info %s -> %s OR %s"
			    , dictName (pp->targetDict, n1)
			    , dictName (pp->targetDict, px->targetGene)
			    , dictName (pp->targetDict, n2)
			    ) ;
	}
      if (! px->targetGene)
	{
	  px->targetGene = n2 ;
	  nn++ ;
	}
      aceInStep (ai,'\t') ;
      if (aceInInt (ai, &x1))
	messcrash ("Sorry 4 columns target2gene file is not yet programmed") ;
    }
  fprintf (stderr, "Parsed %d target2gene relations, %d duplicates\n", nn, nd) ;
  ac_free (h) ;
  return ;
} /* target2GeneParse */

/*************************************************************************************/

static BOOL clipAlignGetProbes (CLIPALIGN *pp, int pass)
{
  AC_HANDLE h ; 
  int state, line = 0, mult = 1, nn = 0, n ;
  long int nTagFound = 0, nSeqBpFound = 0, nTagBpFound = 0 ;
  int probeName, probeNextName, sequence, quality, probeLength=0, ignoredProbe = 0, nTooShort = 0 ;
  int probeSetName = 0, probeSetNextName, previousScoreNext, fragment = 0, probeNextPair = 0 ;
  int probeLeftClip, probeRightClip, previousScore =  - 9999, previousPairScore = - 9999 ;
  int nDecimate = pp->decimate > 0  ? pp->decimate : 0 ; /* if 100 keep on read every 100 */
  int stranded = 0, strandedNext = 0 ;
  char cc, *cp, *cq ;
  MM *mm = 0 ;
  BOOL isFastc = FALSE, fastQ = pp->fastQ, found2 = FALSE ;
  char prefix = fastQ ? '@' : '>' ;
  ACEIN ai = 0 ;
  ACEOUT outNoInsert = pp->outNoInsert ;
  int lengthOtherRead = 0 ;

  if (! pp->probeFileName2 && ! pp->hasPairs && pass == 1)
    goto doSortProbes ;

  if (!pp->probeFileName)	
    {
      fprintf (stderr, "// No probe file provided\n") ;
      exit (1) ;
    }

  h = ac_new_handle () ;
  
  if (! pass)
    {
      pp->probes = arrayHandleCreate (1000000, MM, pp->h) ;
      pp->probeStack = stackHandleCreate (10000000, pp->h) ;
      pushText (pp->probeStack, "___toto") ; /* avoid zero */
      pp->probeDict = dictHandleCreate (5000000, pp->h) ; /* tmp storage of probe names to establish the pairs */
      if (LATESORT == 0) pp->probePairs = arrayHandleCreate (10000000, int, pp->h) ;
      pp->probeLengthMin = -1 ;
      
      /* it is imperative that we do not ceclare
	 because we want to divide the offsets in pp->probeDir by STACL_ALIGNMENT
	 this is a horrible hack, we should find a better way to index the pairs 
      */

      if (pp->fastQ && (pp->fastQbase || pp->probeQquality > -1000))
	{
	  pp->qualityStack = stackHandleCreate (10000000, pp->h) ;
	  pushText (pp->qualityStack, "toto") ; /* avoid 0 */
	}
    }

  if (strstr ( pp->probeFileName, "fastc"))
    isFastc = TRUE ;
  if (pp->probeFileName2 && pass == 1)
    ai = aceInCreate (pp->probeFileName2, 0, h) ;
  else 
    ai = aceInCreate (pp->probeFileName, 0, h) ;
  if (!ai)
    messcrash ("cannot read probe file %s", pp->probeFileName) ;
  aceInSpecial (ai, "\n") ;
  
  state = probeNextName = probeSetNextName = previousScoreNext = 0 ; 
  /* 0: Outside, look for 'prefix' 'sequence identifier'
   * 1: identifier found, accept first line of DNA 
   * 2: search for more DNA lines until prefix (except case fastQ)
   * 3: fastQ: identifier again+quality
   * 4: ok: process the sequence and loop back
   *
   * If a bad line is found, set state=0;continue; 
   * If the next identifier is found set state=1;continue;
   */
  while (state >= 0)
    {
      /* look for the probe indentifier */
      if (probeNextName)
	{
	  probeName = probeNextName ;
	  probeSetName = probeSetNextName ;
	  previousScore = previousScoreNext ;
	  stranded = strandedNext ;
	  probeNextName = probeSetNextName = previousScoreNext = strandedNext = 0 ;
	  fragment = probeNextPair ;
	  state = 1 ;
	}
      else
	{
	  state = -1 ; cp = 0 ;
	  while (state < nDecimate && aceInCard (ai)) /* will close pp->probeFile */
	    {
	      line++ ;
	      cp = aceInWord (ai) ;
	      if (cp && *cp == prefix)
		state++ ;
	    }
	  if (! cp) { state = -1 ; continue ; }
	  else state = 1 ;
	  probeName = stackMark (pp->probeStack) ;
	  pushText (pp->probeStack, cp + 1) ;	
	  dictAdd (pp->probeDict, cp+1, &fragment) ;

          stranded = pass ? - pp->stranded : pp->stranded ;
	  catText (pp->probeStack, pass ? "<" : ">") ;

	  previousScore = - 9999 ;
	  if (pp->previousScoreDict)
	    {
	      int n1 ;
	      
	      cp = stackText (pp->probeStack, probeName) ;

	      if (dictFind (pp->previousScoreDict, cp, &n1) &&
		  n1 < arrayMax (pp->previousScores))
		previousScore = arr (pp->previousScores, n1, int) ;
	    }
	  
	  if ((cp = aceInWord(ai)) && (cq = strstr (cp, "ProbeSetID=")))
	    {
	      cq += strlen ("ProbeSetID=") ;
	      if (*cq)
		{
		  cp = cq + strlen((char *)cq) - 1 ;
		  if (*cp == ';') *cp = 0 ;
		  probeSetName = stackMark (pp->probeStack) ;
		  pushText (pp->probeStack, cq) ;
		}
	    }
	}

      /* look for the probe sequence */
      sequence = quality = 0 ;
      while (state < 3 && aceInCard (ai)) /* will close pp->probeFile */
	{
	  line++ ;
	  cp = aceInWord (ai) ;
	  if (!cp)
	    state = 4 ;
	  else if (state == 2)
	    {
	      if (*cp == prefix && ! (*cp == '>' && *(cp+1)=='<')) /* next identifier has been found */
		{
		  if (nDecimate > 1)
		    {
		      int i = 2 * nDecimate - 1 ;
		      while (aceInCard (ai) && i--) ;
		      if (i < 0 || ! (cp = aceInWord (ai)))
			{ 
			  state = 4 ;
			  continue ;
			}
		    }
		  probeNextName = stackMark (pp->probeStack) ;
		  pushText (pp->probeStack, cp + 1) ;
		  dictAdd (pp->probeDict, cp+1, &probeNextPair) ;
		  catText (pp->probeStack, pass ? "<" : ">") ;
		  state = 4 ; /* identifier is available */
		  
		  strandedNext = pass ? - pp->stranded : pp->stranded ;
		  
		  previousScoreNext = - 9999 ;
		  if (pp->previousScoreDict)
		    {
		      int n1 ;
		      cp = stackText (pp->probeStack, probeNextName) ;

		      if (dictFind (pp->previousScoreDict, cp, &n1) &&
			  n1 < arrayMax (pp->previousScores))
			previousScoreNext = arr (pp->previousScores, n1, int) ; 
		    }

		  if ((cp = aceInWord(ai)) && (cq = strstr (cp, "ProbeSetID=")))
		    {
		      cq += strlen ("ProbeSetID=") ;
		      if (*cq)
			{
			  cp = cq + strlen((char *)cq) - 1 ;
			  if (*cp == ';') *cp = 0 ;
			  probeSetNextName = stackMark (pp->probeStack) ;
			  pushText (pp->probeStack, cq) ;
			}
		    }
		}
	      else
		{
		  if (! found2)
		    {
		      if (! pass)
			{      
			  cq = strchr (cp, '<') ;  /* detect and clip the reverse half of the pair */
			  if (cq)
			    {
			      found2 = TRUE ;
			      pp->hasPairs = TRUE ;     /* trigger a second scanning of the whole file */
			      lengthOtherRead = strlen (cq+1) ;
			      *cq-- = 0 ;
			      if (*cq == '>')
				*cq = 0 ;
			    }
			}
		      else       /* register the second half of the read */
			{
			  cq = strchr (cp, '<') ; 
			  if (cq)
			    {
			      lengthOtherRead = cq - cp - 1 ;
			      if (lengthOtherRead < 0) lengthOtherRead = 1 ;
			      cp = cq + 1 ;
			    }
			  else
			    {
			      lengthOtherRead = strlen (cp) ;
			      cp[1] = 0 ;		  
			    }
			}
		    }
		  if (!pass || found2)
		    catText (pp->probeStack, cp) ;
		}
	    }
	  else 
	    {
	      sequence = stackMark (pp->probeStack) ;
	      found2 = FALSE ;
	      if (! pass)
		{      
		  cq = strchr (cp, '<') ;  /* detect and clip the reverse half of the pair */
		  if (cq)
		    {
		      found2 = TRUE ;
		      pp->hasPairs = TRUE ;     /* trigger a second scanning of the whole file */   
		      lengthOtherRead = strlen (cq+1) ;
		      *cq-- = 0 ;
		      if (*cq == '>')
			*cq = 0 ;
		    }
		}
	      else       /* register the second half of the read */
		{
		  cq = strchr (cp, '<') ; 
		  if (cq)
		    {
		      lengthOtherRead = cq - cp - 1 ;
		      if (lengthOtherRead < 0) lengthOtherRead = 1 ;
		      cp = cq + 1 ;
		    }
		  else
		    {
		      lengthOtherRead = strlen (cp) ;
		      if (! pp->probeFileName2)
			cp[1] = 0 ;		  
		    }
		}
	      pushText (pp->probeStack, cp) ;
	      state = fastQ ? 3 : 2 ;
	    }
	}

      /* process the sequence */
      if (sequence)
	{
	  probeLength = strlen (stackText (pp->probeStack, sequence)) ;
	  cp = stackText (pp->probeStack, sequence) ;
	  cq = cp - 1 ;
	  while (*++cq)
	    {
	      cc = 0 ;
	      if (pp->solid) /* convert to the pseudo-DNA ccfa format */
		{
		  switch (*cq)
		    {
		    case '0': cc = A_ ; break ;
		    case '1': cc = C_ ; break ;
		    case '2': cc = G_ ; break ;
		    case '3': cc = T_ ; break ;
		    case '.': cc = N_ ; break ;
		    }
		}
	      if (!cc)
		cc = dnaEncodeChar[(int)*cq] ;
	      if (pp->decoy) /* convert the read to decoy mode */
		{
		  switch (cc)
		    {
		    case A_: cc = T_ ; break ;
		    case T_: cc = A_ ; break ;
		    case G_: cc = C_ ; break ;
		    case C_: cc = G_ ; break ;
		    default: ; break ;
		    }
		}
	      if (!cc && pp->solid)
		cc = A_ ;
	      if (!cc && *cq == '-')
		cc = N_ ;
	      if (cc == RS_)
		pp->hasPairs = TRUE ;
	      if (cc)
		*cq = cc ;
	      else
		{
		  fprintf (stderr, "// bad character in file %s line %d, expecting atgc...n, got %s\n"
			   , pp->probeFileName, line, cq) ;
		  exit(1) ;
		}		     
	    }
	}
      /* look for the fastQ quality identifier */
      if (fastQ && state == 3)
	{
	  if (aceInCard (ai))
	    {
	      line++ ;
	      cp = aceInWord (ai) ;
	      if (cp  && (!strcmp (cp,"+") || !strcmp (cp+1, stackText (pp->probeStack, probeName))))
		{
		  state = 3 ;
		}
	      else 
		state = 0 ; /* all Hell is loose */
	    }
	  else
	    state = 0 ; /* bad record */
	}
      /* look for the fastQ quality values */
      if (fastQ && state == 3)
	{
	  if (! aceInCard (ai))
	    state = 0 ; /* bad record */
	  else
	    {
	      line++ ;
	      n = 0 ;
	      if (pp->probeQquality > -1000)
		{
		  char qBuf[4096] ;

		  cp = aceInWord (ai) ;
		  /* quality is -10 log10 (p/1-p) : range is -5 to 40,
		   * -5 is 1/4 == N all letters equiprobable log10(1/3)=.5
		   * In char coding @=64 is 0, ';'=59 is -5
		   * q=30 corresponds to p=1/1000
		   */ 
		  if (*cp >= '0' && *cp <= '9') /* number coding */
		    {
		      int ee = 0, m = 1, iq ;
		      int eemin = pp->probeQquality ;
		      n = 1 ;
		      while (aceInInt(ai, &ee)) 
			{
			  if (ee < -5) ee = -5 ; 
			  if (ee < 40) ee = 40 ;
			  iq = '@' + ee ;
			  qBuf[m-1] = iq ;
			  m++ ; 
			  if (ee>=eemin) 
			    n=m ;
			}
		      qBuf [m - 1] = 0 ;
		    }
		  else /* condensed coding   qv=char -'@' */
		    {
		      char cc0 = pp->fastQbase ? pp->fastQbase : 64 ;
		      strncpy (qBuf, cp, 4095) ; qBuf[4095] = 0 ;
		      cq = cp + strlen((char *)cp) - 1 ;
		      while (cq > cp && ((*cq - cc0) < pp->probeQquality))
			*cq-- = 0 ;
		      if (cc0 != '@')
			{
			  cq = cp + strlen((char *)cp) - 1 ;
			  while (cq > cp) { *cq += '@' - cc0 ; cq-- ; }
			}
		      n = cq - cp + 1 ;
		    }
		  quality = stackMark (pp->qualityStack) ;
		  pushText (pp->qualityStack, qBuf) ;

		  if (n < probeLength)
		    {
		      if (n<0) n = 0 ;
		      pp->nseqQclipped++ ;
		      pp->nbpQclipped += probeLength - n ;
		      cq = stackText (pp->probeStack, sequence) ;
		      *(cq+n) = 0 ;
		      probeLength = n ;
		    }
		}
	      state = 4 ;
	    }
	}
      
      /* check for to-be-ignored sequences after parsing the entire record */
      if (
	  (pp->selectPreviousScore && previousScore == - 9999)
	  || 
	  (
	   pp->rejectedProbeDict && 
	   dictFind (pp->rejectedProbeDict, stackText (pp->probeStack, probeName), 0)
	   )
	  )
	{ 
	  state = 0 ; ignoredProbe++ ;
	}
      if (pp->selectedProbeDict && 
	  ! dictFind (pp->selectedProbeDict, stackText (pp->probeStack, probeName), 0)
	  )
	{ 
	  state = 0 ; ignoredProbe++ ;
	}
      if (state < 4)
	continue ;

      /* process special clippings */
      probeLeftClip = probeRightClip = 0 ;
      n = probeLength ;
      if (pp->clipAt > 0)
	{
	  char *ce = stackText (pp->probeStack, sequence) ;
	  if (n > pp->clipAt)
	    *(ce + pp->clipAt) = 0 ;
	  n = strlen (ce) ;
	}
      if (pass == 0 && pp->jump5_1 > 0)
	{
	  int dx = n -1 < pp->jump5_1 ? n -1 : pp->jump5_1 ;
	  sequence += dx ; n -= dx ;
	  if (quality) quality += dx ;
	}
      if (pass == 1 && pp->jump5_2 > 0)
	{
	  int dx = n -1 < pp->jump5_2 ? n -1 : pp->jump5_2 ;
	  sequence += dx ; n -= dx ;
	  if (quality) quality += dx ;
	}
      if (pass == 0 && pp->forceRightClip_1 > 0)
	{
	  char *ce = stackText (pp->probeStack, sequence) ;
	  if (n > pp->forceRightClip_1)
	    *(ce + pp->forceRightClip_1) = 0 ;
	  n = strlen (ce) ;
	}
      if (pass == 1 && pp->forceRightClip_2 > 0)
	{
	  char *ce = stackText (pp->probeStack, sequence) ;
	  if (n > pp->forceRightClip_2)
	    *(ce + pp->forceRightClip_2) = 0 ;
	  n = strlen (ce) ;
	}
      if (pass == 0 && pp->UMI > 0)
	{
	  char *ce = stackText (pp->probeStack, sequence) ;
	  if (n > pp->UMI)
	    *(ce + pp->UMI) = 0 ;
	  n = strlen (ce) ;
	}
      if (n>2 && pp->solid)  /* always discard the first letter in solid, probably a T */
	{                    /* do not count it as a clipping */
	  sequence++ ; n-- ;
	  if (quality) quality++ ;
	}
      if (n > 5 && stranded <= 0 && (pp->clipPolyT))
	{
	  char *ce = stackText (pp->probeStack, sequence) ;
	  int seq1, qual1 ;
	  int n1, pLeftClip ;

	  /* accept exact T */
	  while (n > 5 && *ce == (pp->solid ? A_ : T_))   /* continuation T */
	    { 
	      ce++ ; sequence++ ; n-- ; probeLeftClip++ ;
	      if (quality) quality++ ;
	    }
          seq1 = sequence ; qual1 = quality ; n1 = n ; pLeftClip = probeLeftClip ;
	  while (1)
	    {
	      /* accept one non-T every 8 letters */
	      ce++ ; seq1++ ; qual1++ ; n1-- ; pLeftClip++ ;
	      /* scan again */
	      while (n1>5 && *ce == (pp->solid ? A_ : T_))   /* continuation T */
		{ 
		  ce++ ; seq1++ ; n1-- ; pLeftClip++ ;
		  if (quality) qual1++ ;
		}		
	      if (seq1 > sequence + 8) /* leftClip no problem it is external to the pairs */
		{
		  sequence = seq1 ; 
		  if (quality) quality = qual1 ;
		  n = n1 ;  probeLeftClip = pLeftClip ;
		}
	      else
		break ;
	    }
	}

      if (pp->clipN && ! pp->solid) /* clip n on the right */
	{
	  int i, j, i2 ;
	  char *ce = stackText (pp->probeStack, sequence) ;

	  n = strlen (ce) ;
	  /* start at end of sequence and clip the Ns not followed by clipN + 2 consecutive good letters */
	  for (ce += n - 1, i = i2 = n, j = 0, n = 0 ; i >= 0 ; i--, ce--)
	    {
	      if (*ce == N_) 
		{
		  j = 0 ; i2 = i ? i - 1 : 0 ; 
		} /* candidate clip position */
	      else j++ ;
	      if (i2 == 0 || j >= pp->clipN)
		{ 
		  ce = stackText (pp->probeStack, sequence + i2) ;
		  *ce = 0 ; n = i2 ; 
		  break ;
		}
	    }
	}
      if (1 && pp->exitAdaptor) /* very strict before alignment */
	{
	  int iv, xv, xv2, ln ;
	  unsigned char *ce = (unsigned char *) stackText (pp->probeStack, sequence) ;

	  /* loop backwards to favor gag... over Agag... and obtain a cleaner profile */
	  for (iv = 63 ; iv >= 0 ; iv--)
	    {
	      if (!pp->exitAdaptor[iv] || ! *pp->exitAdaptor[iv])
		continue ;
	      xv = 0 ;  ln = strlen((const char*)pp->exitAdaptor[iv]) ;
	      if (pp->solid && pp->hasPairs)
		{ 
		  xv2 = dnaPickMatch (ce + xv, n - xv, pp->exitAdaptor[iv], 1, 1) ;
		  if (xv2 < 3) { xv = xv2 ; }		  
		}
	      else if (ln >= 24)
		{
		  while ((xv2 = dnaPickMatch (ce + xv, n - xv, pp->exitAdaptor[iv], 1, 1)))
		    { xv += xv2 ;  if (n - xv <= ln) break ; }
		}
	      else if (ln >= 24 && pp->exitAdaptorFound[iv] > 1000)
		{
		  while ((xv2 = dnaPickMatch (ce + xv, n - xv, pp->exitAdaptor[iv], 2, 1)))
		    { xv += xv2 ;  if (n - xv <= ln) break ; }
		}
	      else
		{
		  while ((xv2 = dnaPickMatch (ce + xv, n - xv, pp->exitAdaptor[iv], 0, 0)))
		    { xv += xv2 ; if (n - xv <= ln) break ; }
		}
	      
	      if ((xv  && n - xv >= 18 && n - xv + 1 <= ln) || (xv == 1 && n - xv >= 18) || (pp->solid && pp->hasPairs && xv && xv < 3))
		{
		  if (xv == 1 || xv == 2) 
		    {
		      int i99 ;
		      unsigned char sBuf[256], *ce99 ;

		      mult = isFastc ? fastcMultiplicity (stackText (pp->probeStack, probeName), 0, 0) : 1 ;
		      pp->nSeqEmptyInsert++ ;
		      pp->nTagEmptyInsert += mult ;
		      pp->nBpEmptyInsert  += mult * strlen ((char*)ce) ;
		      pp->exitAdaptorFound[iv] += mult ;
		      if (! outNoInsert)
			{
			  outNoInsert = pp->outNoInsert = aceOutCreate (pp->outFileName, ".noInsert", pp->gzo, pp->h)  ;
			  aceOutf (outNoInsert, "#Read\tMultiplicity\tStart\tEnd\tAdaptor\n") ;
			}
		      i99 = 0 ;
		      if (xv == 2) sBuf[i99++] = 'n' ;
		      for (ce99 = ce ; *ce99 && i99 < 254 ; i99++, ce99++)
			sBuf[i99] = dnaDecodeChar [(int)*ce99] ;
		      sBuf[i99] = 0 ;
		      aceOutf (outNoInsert, "%s\t%d\t%d\t%d\t%s\n", stackText (pp->probeStack, probeName), mult, xv, strlen ((char*)ce), sBuf) ;
		      if (1) /* reject the reverse read if the forward read is rejected */
			{
			  char *cp4 =  stackText (pp->probeStack, probeName), *cp5 ;
			  cp5 = cp4 + strlen (cp4) - 1 ;
			  if (*cp5 == '>') 
			    {
			      *cp5 = '<' ;
			      aceOutf (outNoInsert, "%s\t%d\n", cp4, mult) ;
			      if (! pp->rejectedProbeDict) 
				pp->rejectedProbeDict = dictHandleCreate (100000, pp->h) ;
			      dictAdd (pp->rejectedProbeDict, cp4, 0) ;
			      *cp5 = '>' ;
			    }
			}
		      ce[xv - 1] = 0 ; n = xv - 1 ;
		    }
		  break ;
		}
	      n = strlen ((const char *)ce) ;
	      if (xv && ln >= 15 && xv + ln < n) /* hard clip just downstrem of the adaptor */
		{
		  pp->exitAdaptorFound[iv] += mult ;
		  if (pp->exitAdaptorFound[iv] > 1000)
		    ce[xv + ln] = 0 ;
		  probeLength = strlen ((const char *)ce) ;
		}		
	    }
	}

      if (1 && pp->exitAdaptor2) /* very strict before alignment */
	{
	  int iv, xv, xv2, ln ;
	  unsigned char *ce = (unsigned char *) stackText (pp->probeStack, sequence) ;

	  /* loop backwards to favor gag... over Agag... and obtain a cleaner profile */
	  for (iv = 63 ; iv >= 0 ; iv--)
	    {
	      if (!pp->exitAdaptor2[iv] || ! *pp->exitAdaptor2[iv])
		continue ;
	      xv = 0 ;  ln = strlen((const char*)pp->exitAdaptor2[iv]) ;
	      if (pp->solid && pp->hasPairs)
		{ 
		  xv2 = dnaPickMatch (ce + xv, n - xv, pp->exitAdaptor2[iv], 1, 1) ;
		  if (xv2 < 3) { xv = xv2 ; }		  
		}
	      else if (ln >= 24)
		{
		  while ((xv2 = dnaPickMatch (ce + xv, n - xv, pp->exitAdaptor2[iv], 1, 1)))
		    { xv += xv2 ;  if (n - xv <= ln) break ; }
		}
	      else if (ln >= 24 && pp->exitAdaptor2Found[iv] > 1000)
		{
		  while ((xv2 = dnaPickMatch (ce + xv, n - xv, pp->exitAdaptor2[iv], 2, 1)))
		    { xv += xv2 ;  if (n - xv <= ln) break ; }
		}
	      else
		{
		  while ((xv2 = dnaPickMatch (ce + xv, n - xv, pp->exitAdaptor2[iv], 0, 0)))
		    { xv += xv2 ; if (n - xv <= ln) break ; }
		}
	      
	      if ((xv  && n - xv >= 18 && n - xv + 1 <= ln) || (xv == 1 && n - xv >= 18) || (pp->solid && pp->hasPairs && xv && xv < 3))
		{
		  if (xv == 1 || xv == 2) 
		    {
		      int i99 ;
		      unsigned char sBuf[256], *ce99 ;

		      mult = isFastc ? fastcMultiplicity (stackText (pp->probeStack, probeName), 0, 0) : 1 ;
		      pp->nSeqEmptyInsert++ ;
		      pp->nTagEmptyInsert += mult ;
		      pp->nBpEmptyInsert  += mult * strlen ((char*)ce) ;
		      pp->exitAdaptor2Found[iv] += mult ;
		      if (! outNoInsert)
			{
			  outNoInsert = pp->outNoInsert = aceOutCreate (pp->outFileName, ".noInsert", pp->gzo, pp->h)  ;
			  aceOutf (outNoInsert, "#Read\tMultiplicity\tStart\tEnd\tAdaptor\n") ;
			}
		      i99 = 0 ;
		      if (xv == 2) sBuf[i99++] = 'n' ;
		      for (ce99 = ce ; *ce99 && i99 < 254 ; i99++, ce99++)
			sBuf[i99] = dnaDecodeChar [(int)*ce99] ;
		      sBuf[i99] = 0 ;
		      aceOutf (outNoInsert, "%s\t%d\t%d\t%d\t%s\n", stackText (pp->probeStack, probeName), mult, xv, strlen ((char*)ce), sBuf) ;
		      if (1) /* reject the reverse read if the forward read is rejected */
			{
			  char *cp4 =  stackText (pp->probeStack, probeName), *cp5 ;
			  cp5 = cp4 + strlen (cp4) - 1 ;
			  if (*cp5 == '>') 
			    {
			      *cp5 = '<' ;
			      aceOutf (outNoInsert, "%s\t%d\n", cp4, mult) ;
			      if (! pp->rejectedProbeDict) 
				pp->rejectedProbeDict = dictHandleCreate (100000, pp->h) ;
			      dictAdd (pp->rejectedProbeDict, cp4, 0) ;
			      *cp5 = '>' ;
			    }
			}
		      ce[xv - 1] = 0 ; n = xv - 1 ;
		    }
		  break ;
		}
	      n = strlen ((const char *)ce) ;
	      if (xv && ln >= 15 && xv + ln < n) /* hard clip just downstream of the adaptor */
		{
		  pp->exitAdaptor2Found[iv] += mult ;
		  if (pp->exitAdaptor2Found[iv] > 1000)
		    ce[xv + ln] = 0 ;
		  probeLength = strlen ((const char *)ce) ;
		}		
	    }
	}

      if (pp->clipPolyA /* 2019/08/07 reactivated polya if pair and nA > 12 */ 
	  /* && ! pp->hasPairs  2018_01_19 very unlikely a good idea in pair sequencing */
	  )
	{ 
	  int i, j, i2 ;
	  char *ce = stackText (pp->probeStack, sequence) ;
	  if (stranded >= 0) /* in solid i do not know which base this is */
	    {             /* but i tentatively clip any terminal homopolymer */
	      int n0, n1 = 0 ;
	      int bad = 0 ;
	      int nBad = 0 ;
	      /* accept exact A */
	      n0 = n1 = n = strlen (ce) ;
	      ce += n0 - 1 ;
	      while (n>5 && (*ce & A_))
		{ ce-- ; n-- ;  }
	      n1 = n ; 
	      while (bad<3)
		{
		  /* accept one non-A every 8 letters */
		  ce-- ; n1-- ; nBad++ ;
		  /* scan again */
		  if (nBad < 3 && n1 > 5 && !  (*ce & A_))
		    { ce-- ; n1-- ; nBad++ ;}
		  while (n1 > 5 && (*ce & A_))
		    { 
		      ce-- ; n1-- ; 
		    }
		  
		  if (n1 + 8 * nBad < n)
		    {
		      n = n1 ; nBad = 0 ;
		    }
		  else
		    bad++ ;
		}
	      ce = stackText (pp->probeStack, sequence) ;
	      if (n0 > n + (pp->hasPairs ? 12 : 0))
		{
		  ce[n] = 0 ;
		  probeRightClip = n0 - n ;
		}
	    }
	  if (probeRightClip && ! pp->solid) /* continue clipping while we have 5A/6 */
	    {
	      ce = stackText (pp->probeStack, sequence) ;
	      
	      n = strlen (ce) ;
	      /* start at end of sequence and clip the Ns not followed by clipN + 2 consecutive good letters */
	      for (ce += n - 1, i = i2 = n, j = 0 ; i >= 10 ; i--, ce--)
		{
		  if (*ce & A_) 
		    {
		      if (j<0) j++ ; 
		      if (j>=0) {i2 = i - 1 ; }
		    } /* candidate clip position */
		  else 
		    j-= 6 ;
		}
	      if (i2 < n)
		{ 
		  ce = stackText (pp->probeStack, sequence + i2) ;
		  *ce = 0 ; probeRightClip+= n - i2 ; n = i2 ;
		}
	    }
	}
      previousPairScore = previousScore ;
      previousScore -= lengthOtherRead ;
      lengthOtherRead = 0 ;
      if (previousPairScore < 0) 
	previousPairScore = 0 ;
      if (previousScore < 0) 
	previousScore = 0 ;

      mult = isFastc ? fastcMultiplicity (stackText (pp->probeStack, probeName), 0, 0) : 1 ;
      if (mult < pp->probeMinMultiplicity)
	{
	  state = 0 ;
	  nTooShort++ ;
	  continue ;
	}
      probeLength = n ;
      if (n < 2 || n < pp->probeMinLength || (pp->probeMaxLength > 0 && n > pp->probeMaxLength))  /* drop the probe */
	{
	  state = 0 ;
	  nTooShort++ ;
	  continue ;
	}
      if (previousScore >  (pp->hasPairs ? 2 : 1) * (probeLength  +  pp->bonus + 2 * pp->endBonus)) /* may happen in breakpoint search */
	continue ;
      
      if (state >= 4)
	{
	  float entropy = 0, tm = 0 ;
	  
	  clipAlignGetTm (pp,sequence, probeLength, &entropy, &tm) ; 
	  if (probeLength < 10000 && entropy < pp->minEntropy) 
	    continue ;
	  
	  nn = arrayMax (pp->probes) ;
	  if (pp->mutantMalus)
	    array (pp->originalScores, nn, short) = previousScore ;
	  if (0) fprintf(stderr, "nn=%d Nam=%d name=%s\n", nn, probeName, stackText(pp->probeStack, probeName)) ;
	  mm = arrayp (pp->probes, nn++, MM) ;
	  mm->probeName = probeName ;
	  mm->probeSetName = probeSetName ;
	  mm->probeSequence = sequence ;
	  mm->probeLength = probeLength ; 
	  mm->probeLeftClip = probeLeftClip ;
	  mm->probeRightClip = probeRightClip ;
	  mm->bestScore[0] = previousScore ;
	  mm->bestPairScore = previousPairScore ;
	  mm->probeQuality = quality ;
	  mm->entropy = entropy ;
	  mm->stranded = stranded ;
	  mm->mult = mult ;
	  mm->fragment = fragment ;

	  if (LATESORT == 1)
	    {
	      mm->pair = pass ; 
	    }
	  else
	    {
	      if (pass == 0)
		array (pp->probePairs, fragment, int) = nn + 1 ;
	      else
		{
		  int n = array (pp->probePairs, fragment, int) ;
		  MM *mm1 ;
		  
		  if (n >= 2 && n <= nn )
		    {
		      mm->pair = -n ;
		      mm1 = arrp (pp->probes, n - 2, MM) ;
		      mm1->pair = nn + 1 ;
		    }
		  else
		    mm->pair = -1 ; /* so the sign of pair tells us we are backwards */
		}
	    }

	  nTagFound += mm->mult ;
	  nSeqBpFound += probeLength ;
	  nTagBpFound += mm->mult * probeLength ;;
	  if (probeLength > pp->probeLengthMax)
	    pp->probeLengthMax = probeLength ;
	  if (nn == 1 || probeLength < pp->probeLengthMin)
	    pp->probeLengthMin = probeLength ;
	  previousScore = - 9999 ;
	}
    }

  ac_free (h) ;
  if (0) clipAlignCountSequenceWithBadLetters (pp) ; /* will exit */

  /* accelerate in case of conflicting requests */
  if (pp->minTM > 0)
    {
      /* Maniatis : tm = 62.3 + 0.41 P - 500/L */
      /*            tm = 62.3 + (41.0 * nGC - 500)/N  */

      int L = 0 ;
      float x = 62.3 + 41 - pp->minTM ;
      if (x > 0) L = 500/x ;
      if (pp->probeLengthMin < L -1) 
	pp->probeLengthMin = L - 1 ;
    } 

  if (pp->errMax == 0 && pp->errRateMax == 0)
    { 
      if (! pp->seedLengthRequired) /* not specified by user */
	pp->seedLength = pp->probeLengthMin ; 
      if (0 && pp->seedLength < pp->probeLengthMin) /* user was pessimistic ! */
	pp->seedLength = pp->probeLengthMin ; 
      if (is64 && pp->seedLength > 32) pp->seedLength = 32 ;
      if (!is64 && pp->seedLength > 16) pp->seedLength = 16 ;
    }
  if (!pp->silent) fprintf (stderr, "// clipAlignGetProbes found %u seqs %ld Mb, %ld tags %ld Mb\n"
			    , arrayMax (pp->probes), nSeqBpFound/1000000
			    , nTagFound, nTagBpFound/1000000
			    ) ;
  if (!pp->silent && (pp->rejectedProbeDict || pp->selectedProbeDict))
    fprintf (stderr, "// clipAlignGetProbes ignored %d probes found in %s or absent from %s \n"
	     , ignoredProbe
	     , pp->rejectedProbeFileName ? pp->rejectedProbeFileName : "NA"
	     , pp->selectedProbeFileName ? pp->selectedProbeFileName : "NA"
	     ) ;
  if (!pp->silent && fastQ && pp->probeQquality > -1000) 
    fprintf (stderr, "// Quality clipped in %d sequences  %d terminal bp below quality \'%c\'\n"
	     , pp->nseqQclipped
	     , pp->nbpQclipped
	     , (char)pp->probeQquality
	     ) ;

  if (!pp->silent && nTooShort)
    fprintf (stderr, "// Rejected %d probes shorter than %d bp %s %s\n"
	     , nTooShort
	     , pp->probeMinLength > 2 ? pp->probeMinLength : 2
	     , pp->clipPolyA || pp->clipPolyT ? ", probably because you requested PolyA/T clipping" : ""
	     , pp->probeMaxLength > 0 ? ", or because they were longer than requested probeMaxLength "  : ""
	     ) ;
		
  if (! arrayMax (pp->probes))
    {
      fprintf (stderr
	       , "// Cancelled because no sequence tags were found, or they were all rejected"
	       ) ;
      return 0 ; /* exit (1) ; */
    } 
  
  if (!pp->probeLengthMin)
    {
      fprintf (stderr
	       , "// Cancelled because the shortest probe only has %d <= 7 bp\n"
	       , pp->probeLengthMin) ; 
      return 0 ; /* exit (1) ; */
    }


  if (pp->probeLengthMin < 7)
    {
      fprintf (stderr
	       , "// Cancelled because the shortest probe only has %d <= 7 bp\n"
	       , pp->probeLengthMin) ;
      return 0 ; /* exit (1) ; */
    } 

 doSortProbes:

  if (! pp->targetFileName ) /* no need to create a hash in these cases */
    {
      if (pp->probeLengthMax == pp->probeLengthMin)
	clipAlignCountLetters (pp) ;
    }
  else if (pass && LATESORT == 0)
    {
      /* the idea is to sort alphabetically by DNA order the union of the 2 reads */
      int i, iMax = arrayMax (pp->probes), n , n1 ;
      MM *mm1 ;
      BigArray big = bigArrayCreate (stackMark (pp->probeStack), int) ;

      /* replace pair offset by probeName before sorting aalphabetically by sequence */
      for (i = 0, mm = arrp (pp->probes, i, MM) ; i < iMax ; i++, mm++)
	{ 
	  n = 1 ; n1 = mm->pair ; if (n1 < 0) { n = -1 ; n1 = - n1 ; }
	  if (n1 > 1)
	    {
	      mm1 = arrp (pp->probes, n1 - 2, MM) ;
	      mm->pair = n * mm1->probeName ;
	    }
	}
      /* sort */
      probePairOrderStack = pp->probeStack ;
      arraySort (pp->probes, probePairOrder) ;

      /* remap the mm in their new order */
      for (i = 0, mm = arrp (pp->probes, i, MM) ; i < iMax ; i++, mm++)
	bigArray (big, mm->probeName, int) = i + 2 ;
      /* reestablish a direct mapping */
      for (i = 0, mm = arrp (pp->probes, i, MM) ; i < iMax ; i++, mm++)
	{
	  n = 1 ; n1 = mm->pair ; if (n1 < 0) { n = -1 ; n1 = - n1 ; }
	  if (n1 > 1)
	    {
	      mm->pair = n * bigArray (big, n1, int) ;
	    }
	}
      bigArrayDestroy (big) ;
      if (iMax > 1)
	for (i = 1, mm = arrp (pp->probes, 1, MM) ; i < iMax ; i++, mm++)
	  if (nn > 1)
	    {
	      register int same = 0 ;
	      register char *cp, *cq ;
	      
	      cp = stackText (pp->probeStack, (mm-1)->probeSequence) ;
	      cq = stackText (pp->probeStack, mm->probeSequence) ;
	      while (*cp && *cp == *cq) {cp++ ; cq++ ; same++ ; }
	      (mm-1)->same = same ;
	    }
    }  

  if (LATESORT == 0 && pp->hasPairs && pass == 1)
    {
      MM *mm ;
      int i, n, n1, iMax = arrayMax (pp->probes) ;

      for (i = 0, mm = arrp (pp->probes, i, MM) ; i < iMax ; i++, mm++)
	{
	  n = 1 ; n1 = mm->pair ; if (n1 < 0) { n = -1 ; n1 = - n1 ; }
	  n = n1 - 2 ; n1 = i ;
	  if (n >= 0 && n < i) n1 = n ;
	  mm->fragment = n1 + 1 ;
	}
    }
 
  if (pp->hasPairs)
    pp->pairTampon = 10 ;
  return pp->probeLengthMin ; 
} /* clipAlignGetProbes */

/*************************************************************************************/
/* sort all probes from > and < reads and espablish 'same' and 'pair correspondance' */

static Stack clipAlignProbeAlphaOrderStack = 0 ;
static int clipAlignProbeAlphaOrder (const void *va, const void* vb)
{
  const MM *ma = (MM *)va, *mb = (MM *)vb ;
  const char *ccpa =  stackText (clipAlignProbeAlphaOrderStack, ma->probeSequence) ;
  const char *ccpb =  stackText (clipAlignProbeAlphaOrderStack, mb->probeSequence) ;
  return  strcmp (ccpa, ccpb) ;
}

/******/

static void clipAlignSortProbes (CLIPALIGN *pp, int np0, int np1)
{
  register int i, iMax = arrayMax (pp->probes) ; 
  register int same = 0, *ip, *iq ;
  register char *cp, *cq ;
  register int dx = sizeof(int) ;
  register MM *mm ;

  /* sort the probes, mixing the > and the < reads */
  clipAlignProbeAlphaOrderStack = pp->probeStack ;
  arraySort (pp->probes, clipAlignProbeAlphaOrder) ;
  
  /* using the new order, establish the pairs, the sign gives the pass  */
  if (np0 > 0 && np1 > 0)
    {
      Array pairs = arrayCreate (np0 > np1 ? np0 + 1 : np1 + 1, int) ;

      for (i = 1, mm = arrp (pp->probes, i, MM) ; i < iMax ; i++, mm++)
	{
	  same = array (pairs, mm->fragment, int) ;
	  if (same)
	    {
	      register MM *mm1 = arrp (pp->probes, same, MM) ;
	      mm1->pair = mm1->pair > 0 ? i + 2 : -(i + 2) ;
	      mm->pair = mm->pair > 0 ?  same + 2 : - (same + 2) ;
	    }
	  else
	    {
	      array (pairs, mm->fragment, int)  = i ;
	    }
	}
      arrayDestroy (pairs) ;
    }
  /* else we stay with mm->pair == pass ? 1 : -1 */

  /* probes are sorted alphabetically, mm->same tell how many base in the next probe */
  for (i = 1, mm = arrp (pp->probes, i, MM) ; i < iMax ; i++, mm++)
    {
      same = 0 ;
      cp = stackText (pp->probeStack, (mm-1)->probeSequence) ;
      cq = stackText (pp->probeStack, mm->probeSequence) ;
      ip = (int *) cp ; iq = (int *) cq ;
      if (! ip || ! iq) continue ;
      while (*ip == *iq) {ip++ ; iq++ ; same += dx ; }
      cp += same ; cq += same ;
      while (*cp && *cp == *cq) {cp++ ; cq++ ; same++ ; }
      (mm-1)->same = same ;
    }
} /* clipAlignSortProbes */

/*************************************************************************************/
/*************************************************************************************/
/* construct a table of frequency on the sequences given on stdin */
static long wffInserted = 0 ;
static void wordFrequencyDoRegister (CLIPALIGN *pp, Array dna)
{
  int i, n, size = pp->wffSize, ok, step = 0 ;
  long nn = ((long)1) << (2*pp->wffSize) ;
  unsigned char *cp, *bb ;
  unsigned int cc ;
  long  mask = nn - 1, oligo = 0 ;

  if (size == 16) mask = 0xffffffff ;
  if (! pp->wff)
    {
      pp->wff = malloc (nn) ; /* initialize as ZERO */
      memset (pp->wff, 0, nn) ;
    }
  bb = pp->wff ;
  for (i = arrayMax (dna), ok = 0,  cp = arrp (dna, 0, unsigned char) ; i > 0 ; i--, cp++)
    switch ((cc = *cp))
      {
      case A_:
      case T_:
      case G_:
      case C_:
	oligo <<= 2 ; oligo |= B2[cc] ; oligo &= mask ;
	ok++ ;
	if (ok >= size && ! (step++ % 1))
	  {
	    n = bb[oligo] ;
	    if (n < 255) n++ ;
	    bb[oligo] = n & 0xff ;
	    wffInserted++ ;
	  }
	break ;
      default:
	ok = 0 ;
	break ;
      }
  return ;
} /* wordFrequencyDoRegister */

/*************************************************************************************/

static void wordFrequencyRegister (CLIPALIGN *pp, Array dna)
{
  dnaEncodeArray (dna) ;
  wordFrequencyDoRegister (pp, dna) ;
  if (! pp->strand)
    {
      reverseComplement (dna) ;
      wordFrequencyDoRegister (pp, dna) ;
    }
} /* wordFrequencyRegister */

/*************************************************************************************/

static void wordFrequencyReport (CLIPALIGN *pp)
{
  long nn[256] ;
  unsigned char *cp = pp->wff ;
  int i ;
  long n = ((long)1) << (2*pp->wffSize) ;

  memset (nn, 0, sizeof(nn)) ;
  while (n--) nn[(int)*cp++]++ ;

  fprintf (stderr, "Number of occurence\tNumber of words\tCumul\tMissed\n") ;
  for (i = n = 0 ; i < 255 ; i++)
    {
      n += i*nn[i] ;
      fprintf (stderr, "%d\t%ld\t%ld\t%ld\n", i, nn[i], n, wffInserted - n) ;
    }
  fprintf (stderr, ">254\t%ld\t%ld\n", nn[255],  n+255*nn[255]) ;
  fprintf (stderr, "At least %ld words out of %ld\n", n+255*nn[255], wffInserted) ;

  return ;
} /* wordFrequencyReport */

/*************************************************************************************/

static void wordFrequencyExport (CLIPALIGN *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreateToFile (pp->wordFrequencyConstructTable, "wb", h) ;
  char buf[64] ; /* title line */
  long  nn ;

  if (! ao)
    messcrash ("Sorry I cannot open for writing the file wordFrequencyConstructTable %s\n"
	       , pp->wordFrequencyConstructTable
	       ) ;

  memset (buf, 0, sizeof(buf)) ;
  buf[0] = 0xff & pp->wffSize ;
  strcpy (buf+4, "wordFrequencyTable version 1") ; 
  aceOutBinary (ao, buf, 64) ;

  nn = 1 << (2*pp->wffSize) ;
  if (aceOutBinary (ao, (char *)pp->wff, nn))
    messcrash ("Cannot write 4^%d bytes in file %s", pp->wffSize, pp->wordFrequencyConstructTable) ;
  else
    fprintf (stderr, "# Wrote 4^%d = %ld bytes in file %s\n",  pp->wffSize , nn, aceOutFileName(ao)) ;
  ac_free (h) ;
  return ;  
} /* wordFrequencyExport */
 
/*************************************************************************************/

static void wordFrequencyParseTable (CLIPALIGN *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreateFromFile (pp->wordFrequencyTable, "rb", 0, h) ;
  char buf[64] ; /* title line */
  int nn ;

  aceInBinary (ai, buf, 64) ;
  pp->wffSize = buf[0] ;
  if (strcmp (buf+4, "wordFrequencyTable version 1"))
    messcrash ("Cannot read the correct header in file %s", pp->wordFrequencyTable) ;
  nn = 1 << (2*pp->wffSize) ;
  pp->wff = malloc (nn) ; /* avoid halloc (nn, pp->h) so we do not init to zero which is slow */
  if (! aceInBinary (ai, (char *)pp->wff, nn))
    messcrash ("Cannot read 4^%d bytes in file %s", pp->wffSize, pp->wordFrequencyTable) ;
  
  ac_free (h) ;
  return ;  
} /* wordFrequencyParseTable */

/*************************************************************************************/

static void wordFrequencyConstructTable (CLIPALIGN *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  int n = 0, n1 ;
  Array dna = arrayHandleCreate (1<<30, char, h) ;
  ACEIN ai = aceInCreate (pp->targetFileName, pp->gzi, h) ;
  char *cp ;

  if (pp->wffSize > 25)
    messcrash ("Sorry you did set the length of the words to a value %d > 25", pp->wffSize) ;
  if (!pp->wffSize)
    messcrash ("Sorry you did not set the length of the words to be tabulated by wordFrequencyConstructTable, i suggest 14") ;

  while (aceInCard (ai))
    {
      cp = aceInWord (ai) ;
      if (cp && *cp != '>')
	{
	  n1 = strlen (cp) ;
	  array (dna , n+n1, unsigned char) = 0 ; /* make room */
	  arrayMax (dna) = n+n1 ;
	  memcpy (arrp (dna, n, unsigned char), cp, n1) ;
	  n+= n1 ;
	}
      else    
	{
	  if (n) wordFrequencyRegister (pp, dna) ;
	  arrayMax (dna) = n = 0 ;
	}
    }
  if (n) wordFrequencyRegister (pp, dna) ;
  wordFrequencyReport (pp) ;
  wordFrequencyExport (pp) ;
  
  ac_free (h) ;
  return ;
} /* wordFrequencyConstructTable */

/**********************/

static int parseTargetMaskFile (CLIPALIGN *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0, a1, a2, chrom, ok ;
  const char *ccp ;
  PEXPORT *up ;
  Array aa = 0 ;
  ACEIN ai = aceInCreate (pp->targetMaskFileName, FALSE, h) ;

  if (ai)
    {
      aa = arrayHandleCreate (1000, PEXPORT, pp->h) ;
      while (aceInCard (ai))
	{
	  ccp = aceInWord (ai) ;
	  if (! ccp || *ccp == '#' || *ccp == '/')
	    continue ;
	  dictAdd (pp->targetDict, ccp, &chrom) ;
	  ok = 0 ;
	  aceInStep (ai, '\t') ; 
	  if (aceInInt (ai, &a1)) ok++ ;
	  aceInStep (ai, '\t') ; 
	  if (aceInInt (ai, &a2)) ok++ ;
	  if (ok == 2)
	    {
	      up = arrayp (aa, nn++, PEXPORT) ;
	      up->target = chrom ;
	      if (a1 < a2) { up->a1 = a1 ; up->a2 = a2 ; }
	      else  { up->a1 = a2 ; up->a2 = a1 ; }
	    }
	}
      if (nn) arraySort ( aa, pExportHitOrder) ;
      else arrayDestroy (aa) ;
    }
  ac_free (h) ;
  pp->targetMask = aa ;
  return nn ;
}

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: clipalign  -t target_fasta_file -i probe_fasta_file [-errMax] -... \n"
	   "//      try: -h --help --hits_file_caption \n"
	   "// Example:  clipalign -i tags.fasta -t chromX.fasta -errMax 2\n"
	   "// -t fileName : the name of a target fasta file on which we align the tags\n"
	   "//      this file is read as needed, so there is no limit to its size\n"
	   "//      all files named *.gz are automatically decompressed\n"
	   "// -targets : multi targets\n"
	   "//     example A_mito:1:f.mito.fasta.gz,B_rrna:2:f.rrna.fasta.gz\n"
	   "//     each file is processed, the number is the bonus for that class\n"
	   "// -sam : export in SAM format\n"
	   "// -i : the name of a fasta/fastc/fastq file, possibly .gz, containing the tags to be aligned \n"
	   "//      (tested with 10 Million tags, upper limit depends on hardware)\n"
	   "// -fastq33 : the input is in fastq, the quality of the mismatches start at 33=!(NCBI)\n"
	   "// -fastq64 : the input is in fastq, the quality of the mismatches start at 64=@(ILM native)\n"
	   "// -minQuality int : [default ignore the qualities] requires -fastq33 or -fastq64 \n"
	   "//      right clip until a letter of given quality is reached, replace letters up to this quality by n\n"
	   "//      and report the quality of all missmatches\n"
	   
	   "// -solid : tags are coded a la SOLID, as base transitions\n"
	   "//          the target is still given as FASTA (optionally .gz)\n"
	   "//          the SOLID error detection procedure is applied, compatible systems of\n"
	   "//          transition mismatches are decoded and exported as sequence mismatches\n"
	   "//          uncompatible or isolated transitions are counted as ambiguities and listed as oo\n"
	   "//          In the -solid option, YOU MUST specify a strategy, because the error detection\n"
	   "//          fine tunes the cost of the mismatches and of the oo in a rather complex way\n" 
	   "// -pacbio :  expects long reads with a rather high error rate \n"
	   "// -nanopore :  expects very long reads with a rather very high error rate \n"
	   "// -ultima :  expects short reads with a very low error rate \n"
	   "// -hit fileName : hits file, possibly .gz, output of previous search, to be postprocessed\n"
	   
	   "// -o fileName : output file name, equivalent to redirecting stdout\n"
  
	   "// -gzo : the output files will be gziped\n"
  
	   "// -strand : if this option is specified, the tags are only aligned on the forward strand of the target.\n"
	   "//            This option may be useful if the target is mRNA and the tags are stranded\n"
	   
	   
	   "// -antistrand : if this option is specified, the tags are only aligned on the reverse strand of the target,\n"
	   "//            this option may be useful if the target is mRNA and the tags are antistranded\n"
	   
	   
	   "// -antiprobe : if this option is specified the tag is reverse complemented complemented on parsing,\n"
	   "//            this option seems useful when mapping probes from the AFFY Exon Array\n"
	   
	   "// -A2G | -T2C : analysis of intense RNA editing \n"
	   "//           These 2 options are exclusive of each other, only A2G is reasonable for a transcript target\n"
	   "//           In all cases, since the 4 bases are mapped to only 3, the seedLength is scaled up by log(4)/log(3)\n"
	   "//         A2G : edit on the top strand of the genome every A to R == A|G == puRine: analyze top strand RNA-edited genes\n"
	   "//         T2C : edit on the top strand of the genome every T to Y == T|C == pYrimidine : analyze bottom strand RNA-edited genes\n"
	   "//  -G2A, T2C, C2T, A2C, C2A, T2G, G2T, C2G, G2C, -A2T, -T2A : idem, used as negative controls to verify they are less abundant than A2G\n"
	   "// -probeMinMultiplicity : if this option is specified, the tags under this limit are ignored, only useful if the reads are in tc or fastc input format,\n"
	   "// -probeMinLength : if this option is specified, the tags shorter than this limit are ignored,\n"
	   "// -probeMaxLength : if this option is specified, the tags longer than this limit are ignored,\n"
	   
	   
	   "// -aceOutGenomic : export the alignments in .ace format, using the Genomic_hit schema\n"
	   
	   
	   "// -aceOutMrna : export the alignments in .ace format, using the mRNA_hit schema\n"
	   
	   
	   "// -aceOutNM : export the alignmens in .ace format, using the RefSeq NM_hit schema\n"
	   
	   
	   "// -getTm : Report the Tm (melting temperature) and S (entropy counted in bp) of all tags\n"
	   
	   
	   "// -errMax number : if non zero, the program will accept up to that number of mismatches\n"
	   "//           i.e. the max of errMax and length_tag * errRateMax/100\n	   "
	   "//           but the program becomes heuristic and may have false negatives\n"
	   
	   
	   "// -errRateMax number : if non zero, the program will accept up to that %% of mismatches\n"
	   "//           i.e. the max of errMax and length_tag * errRateMax/100\n"
	   "//           but the program becomes heuristic and may have false negatives\n"
	   
	   
	   "// -minAli number [-minAliPerCent number]: Minimal alignment mode\n"
	   "//           This is the prefered mode to find discontinuous alignments, e.g. RNA against DNA\n"
	   "//           In this mode, alignments around the exact seed are extended optimizing\n"
	   "//              score = (length aligned) - (number of errors) * (error cost)\n"
	   "//           and alignmnents with score > minAli are kept\n"
	   "//           The parameters errMax and errRateMax are ignored\n"
	   "//           minAliPerCent is used to impose a minimal %% of aligned bases in the read\n"
	   
	   
	   "// -strategy [ Exome | Genome | RNA_seq ]\n"
	   "//           Implies a default seed scanning -seedOffset 1 -seedShift 5\n"
	   "//           Optimizes the treatment of introns and of alignemnts on\n"
	   "//           multiple transcripts of the same gene\n" 

	   "// -errCost number : used in -minAli mode: cost of an error or mismatch in a local alignment\n"
	   "//           The cost is the same for an insertion/deletion and a substitution.\n"
	   "//           The default value is 8bp. \n"
	   "//           Please do not modify this parameter\n"
	   
	   "// -targetBonus number : systematic bonus to be added to the score\n"
	   "//           A negative number is effectivelly a malus\n"
	   "//           This parameter is useful in recursive alignments to favor the more likely targets\n"
	   
#ifdef JUNK	   
	   "// -mutantMalus : systematic malus subtracted for the score of targtes called M:*\n"
	   "//           This is a very specific malus used in snp validation\n"	   
#endif	   
	   "// -errMin number : if non zero, only report hits with at least this number of errors\n"
	   "//                  This may be useful to run the program recursivelly, in conjunction with the -i option\n"
	   
	   
	   "// -exactTargetStart number : if non zero, no missmatches between this position and exactTargetStop,\n"
           "//                   including no SOLiD auto-correction, in target coordinates\n"
	   
	   
	   "// -exactTargetStop number : these options are useful to find tags exactly matching an exon junction\n"
	   
	   
	   "// -targetStart number : if non zero, this position up to targetStop must be part of the alignment,\n"
	   
	   "// -targetStop number : these options are useful to find tags exactly matching an exon junction\n"
	   
	   
	   "// -seedLength number : must be >7, specifies the seed length, the default is 16 or the length of the shortest tag\n"
	   "//                  However, if errMax is zero, the seed length is allways the length of the shortest tag.\n"
	   
	   
	   "// -seedOffset number : position of the first base of the seed\n"
	   "//                tgas shorter than seedOffset + seedLength are ignored\n"
	   "// -seedShift number : optional , iterate shifting the seedOffset untill te end of the longest tag\n"
	   
	   
	   "// -select fileName : name of a file, possibly .gz,  of selected tags, one name per line\n"
	   "//           any tag not listed in this file will be ignored when reading the tags file\n"
	   "//           This is very useful to perform a recursive mapping\n"
	   
  
	   "// -reject fileName : name of a file, possibly .gz,  of rejected tags, one name per line\n"
	   "//           any tag listed in this file will be ignored when reading the tags file\n"
	   "//           This is also very useful to perform a recursive mapping\n"
	   
	   "// -decimate number : align only one read every n, useful to check the quality of a large set of reads\n"
	   "//               this parameter may be set to say 7 using \"setenv MAGIC_DECIMATE 7\"\n"
  
           "// -maxDiscardableAli number : (default 90)\n"
	   "//         useful to characterize rearrangments, do not discard a secondary chain longer than this limit\n" 
	   "// -best    : only report hits with the least number of errors\n"
	   "// -nonBest : report any hit up to the requested maximal number of errors\n"
	   "//            The default is to report all hits exceeding maxDiscardableAli\n"
	   
	   "// -vectorize :  export the alignments -a la Jim Kent- with coma separated muti coordinates\n"

	   "// -maxHit number : default 10 maximal number of hits reported for any given tag\n"
	   "//                  The exact alignments are always reported, but if maxHit is too low it may happen\n"
	   "//                  that maxHit alignments with 2 errors are found, stopping the search, before the \n"
	   "//                  first alignment with 1 error is found.\n" 
	    
  
	   "// -splice [-overhangLength int]: starting from an exact seed, slide left and right and the position where the tag and the genome diverge\n"
	   "//         To facilitate searching pairs, the exon and overhang are reported, Janus like, on the strand going away from the intron\n"
	   "//         This option is useful when searching new exon junctions:\n"
	   "//           the reported sequence of the exon of the donor matches the reported overhang of the acceptor and vice-versa\n"
	   "//           overhangLength: [default 8] is the minimal nuber of base after the error that should be found in the next 20kb\n"
	   "//    -intronMinLength number -intronMaxLength <number> : defaults 30 and 100000 fit C.elegans, prefer min 60 in human\n"
	   "//           used in -splice mode to in the .overhangs table to report candicate introns and smll deletions\n"
	   "// -stranded n : we recommend 10 or -10, applies to sequencing protocols preserving the strand of the RNA\n"
	   "//         If n > 0, the score for alignemnts to genes antisense relative to the read is diminished by n\n"
	   "//         If n < 0, the score for alignemnts to genes sense relative to the read is diminished by n\n"
	   "//         This applies in particular to cases where there are genes annotated on both strands of a region of the genome\n"
	   "//         If n>0 the alignment to the gene which is sense relative to the read will be preferred\n"
	   "//         If n<0 the number is negative, the alignment to the antisense gene will be preferred\n"
	   "//         In fact the alignment to the genome itself will win over a gene on the wrong strand\n"
	   "//         except if the reads jumps an intron, this occurs in the rather exceptional case of RNAs copied\n"
	   "//         in vivo on the standard messenger RNA by an RNA dependent RNA-polymerase\n"

	   "// -minTM number : if non zero, minimal TM of the longest exact segment (TM=62.3 + 0.41 Percent(G+C) - 500/L)\n"
	   "//              This option is useful to study eventual cross hybridization of micro array tags\n"
	   
  
	   "// -minEntropy number : default 7, tags or alignemnts with lower entropy are skipped (Shannon's sum on ATGC of -p log4(p))\n"
	   
	   "// -previousScore filename : only hits equal to or exceeding the previous score are considered\n"
	   "//         Format is 2 column: tag_name  previous_score\n"
	   "// -selectPreviousScore filename : in addition any tag not listed in this file will be ignored\n"
	   "//         Format is 2 column: tag_name  previous_score\n"
	   "// -clipN <int> : clip on N not preceded or followed by at least <int> good letters\n\n"

	   "// -clipPolyA : discard any number of 3'-terminal A, report poly-A tails not matching the target\n\n"
	   
  
	   "// -SclipPolyA : idem, but only if hitting the plus strand of the target\n\n"
	   
  
	   "// -clipPolyT : discard any number of 5'-leading T, report poly-T beginning not matching the target\n\n"
	   
  
	   "// -SclipPolyT : idem, but only if hitting the minus strand of the target\n"
	   "//            This option is TRUE in -solid mode, to compensate an artefact of that technique\n\n"
	   
	   "// -avoidPseudoGenes : high level filter specific to the AceView annotation\n"
	   "//       if all transcripts of a gene are named *unspliced*, hits to this gene are discarded\n"
	   "//       in favor of hits to other genes on the same strand, exluding the case of antisense pairs of genes\n"
  
	   "// -slBonus number : C.elegans specific transpliced leader bonus\n"
	   "//           Triggers SL search and, if successful, contributes to the score\n\n"
	   
  
	   "// -jump5   number : discard so many bases from the beginning (5' end) of each read 1\n"
	   "// -jump5_2 number : discard so many bases from the beginning (5' end) of each read 2\n" 
	   "// -forceRightClip   number : clip read 1 right of this position (short reads followed by an adaptor\n"
	   "// -forceRightClip_2 number : clip read 2 right of this position (short reads followed by an adaptor\n"

	   "// -clipAt number : clip each tag at this length (estimated before jump5)\n\n"
	   
	   "// -entryAdaptor string : comma delimited list of entry adaptors\n"
	   "//           Example:  -entryAdaptor  atgcggtatagg,ggttatmaycc\n"
	   "//           Each adaptor must be >= 8bp, ambiguous letters are accepted\n"
	   "//           Triggers adaptor search upstream of the aligned part of the read\n"
	   "//           If successful the adaptor is trimmed in this particular alignment\n"
	   "//           The adpator must be present backwards relative to the read\n"
	   "//           so in general exitAdaptor == entryAdaptor2 and exitAdaptor2 == entryAdaptor"
	   "// -entryAdaptor2 string : comma delimited list of entry adaptors\n"
	   "//           Same for the second read in paired sequencing\n"
	   "// -exitAdaptor string : comma delimited list of exit adaptors\n"
	   "//           Example: -exitAdaptor atgcggtatagg,ggttatmaycc\n"
	   "//           Each adaptor must be >= 8bp, ambiguous letters are accepted\n"
	   "//           Triggers adaptor search downstream of the aligned part\n"
	   "//           and at the same time downstream of read 2, i.e. downstream of the fragment\n"
	   "//           If successful the adaptor is trimmed in this particular alignment\n"
	   "//           he adaptor is presented in the orientation of the read\n "
	   "// -exitAdaptor2 string: comma delimited list of exit adaptors\n"
	   "//           Each adaptor must be >= 8bp, ambiguous letters are accepted\n"
	   "//           Triggers adaptor search downstream of the aligned part of read 2 in paired sequencing\n"
	   "//           and at the same time upstream of read 1, i.e. upstream of the fragment\n"
	   "//           If successful the effective length of the tag is trimmed\n"
	   
	   "// -UMI int: length of cell UMI bar code in read 1, export it as column 20\n//             Do not map read 1\n"

  	   "// -showTargetPrefix : export 30bp immediatly upstream of the alignment, to study eventual cloning biais\n"
	   "// -showOverhang : export the sequence bits overhanging out of the aligned segment (implied by -splice)\n"
	   "// -showSequence : export the tag sequence in a additional column\n"
	   "// -target_class : target class (A_mito ... Z_genome) reported in the .hits file, useful in post-processing\n"
	   "// -target2gene file_name : associate the -t target fasta files to gene names\n"
	   "//     The file is tab delimited. It may have 2 columns: sequence gene, i.e. mrna gene\n"
	   "//     or 4 columns sequence gene begin end, i.e. chromosome gene begin end\n"
	   "//     Alternativelly, the gene may be given in the -t identifiers as  >sequence|Gene|Gene_name\n"
	   "//\n"
	   "// -targetMask maskFile : mask part of the target\n"
	   "//     The mask file should have 3 columns : identifier coordinat coordinate (in bp)\n"
	   "//     where the identifiers match the identifiers given in the target file\n"
	   "//     The corrsponding DNA regions are blanked, preventing mapping to those regions\n"
	   "//     We recommend to mask the genomic echos of the mitochondria but not the ALU regions\n"
	   "// -decoy : decoy the target\n"
	   "//     Complement the target without reversing it, this creates \"Martian\" DNA with the same statistics\n"
	   "//     Nothing is expected to align on a decoyed target, this provides a control of the specificity of the alignements\n"
	   "// -wordFrequencyConstructTable filename : construct a binary table of frequency of short words\n"
	   "//                     found in the '-t targetfile'\n"
	   "//                     both strands of the starget are used unless you specify '-strand'\n"
	   "// -wfSize int : size of the words to be tabulated, I suggest 14\n"
	   "//                     This table is then given as a parameter to the aligner using the option\n"
	   "// -wf filename : binary table of frequency of short word, used to select seeds\n"
	   "//          This table is needed by the classic AceView 'tacembly' aligner, but so far not used by this program\n"
	   "//\n"
	   "// Optional parameters used in the samStats output file\n"
	   "//   -run <runName>\n"
	   "//   -method <methodName>\n"
	   "//   -nRawReads <long int>\n"
	   "//   -nRawBases <long int>\n"
	   "//\n"
	   "// -silent : suppress title lines and status reports from the output, just report the hits\n"
	   
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

static void clipAlignMemTest (void)
{
  mysize_t giga = 1<<30 ;
  int ii, mem, mx ;
  void * v ;
  int imax = 5 ;
  AC_HANDLE h1 = ac_new_handle () ; 
  
  fprintf (stderr, "imax=%d\n", imax) ;
  
  for (ii = 0 ; ii < imax ; ii++)
    {
      messAllocStatus (&mem) ;
      messAllocMaxStatus (&mx) ;
      if ((v = halloc (giga, h1)))
	fprintf (stderr, "Memory test allocated %d Gb, total %d Mb, max %d Mb\n", ii+1, mem, mx) ;
      else
	messcrash ("Memory test failed to allocate %d Gb, %d Mb allready allocated, max %d Mb, sorry\n", ii+1, mem, mx) ;
    }
  ac_free (h1) ;
} /* clipAlignMemTest */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  extern long int assBounce , assBucketMax, assFound , assNotFound , assInserted , assRemoved  ;
  char *cp ;
  int nHits = 0, ix, np0, np1 ;
  CLIPALIGN p ;
  AC_HANDLE h = 0 ;
  char commandBuf [32000] ;
  const char *ccp = 0 ;

  freeinit () ; 
  messErrorInit (argv[0]) ;

  if (sizeof(void*) == 8)
    is64 = TRUE ;

  h = ac_new_handle () ;
  memset (&p, 0, sizeof (CLIPALIGN)) ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
  p.h = h ;

  if (argc < 2)
    usage (commandBuf, argc, argv) ;
  if (getCmdLineBool (&argc, argv, "-help") ||
      getCmdLineBool (&argc, argv, "-h")
      )
    usage (commandBuf, 1, argv) ;

  if (getCmdLineBool (&argc, argv, "-version"))
    {
      fprintf (stderr, "clipalign: 1.3.0\n") ;
      exit (0) ;
    }

  if (getCmdLineBool (&argc, argv, "-rna") ||
      getCmdLineBool (&argc, argv, "-RNA")
      )
    p.clipPolyA = p.clipPolyT = TRUE ;
  
  if (1)
    {
      int ii, jj, mx = 0, mega = 1 << 20 ;
      AC_HANDLE h1 = ac_new_handle () ;

      if (0 || getCmdLineInt (&argc, argv, "-memTest", &mx))
	{
	  /*	    bucketAssTest () ; */ 

	  if (0)
	    {
	      clipAlignMemTest () ;
	    }

	  if (0) 
	    {   /* pushText crashes at 4gb */ 
	      Associator ass = 0 ;
	      void *v  ;

	      h1 = ac_new_handle () ; 
	      ass = assBigCreate (2000000) ;
	      for (ii = 0 ; ii < mx ; ii++)
		{
		  for (jj = 0 ; jj < 100*mega ; jj++)
		    {
		      v = assVoid (1 + jj + ii * 100 * mega) ;
		      assInsert (ass, v, v) ;
		    }
		  fprintf (stderr, "Memory test assInserted %u M\n", (unsigned int) ((ass->n) >> 20)) ;
		}
	      ac_free (h1) ;
	    }
	  if (1)  exit (0) ;
	} 
      ac_free (h1) ;
    }

  for (ix = 0, cp = commandBuf ;  ix < argc && cp + strlen (argv[ix]) + 1 < commandBuf + 31000 ; cp += strlen (cp), ix++)
    sprintf(cp, "%s ", argv[ix]) ;
  /* consume optional args */
  getCmdLineOption (&argc, argv, "-o", &(p.outFileName)) ;
  getCmdLineOption (&argc, argv, "-i", &(p.probeFileName)) ;
  getCmdLineOption (&argc, argv, "-i2", &(p.probeFileName2)) ;
  getCmdLineOption (&argc, argv, "-query", &(p.probeFileName)) ;
  getCmdLineOption (&argc, argv, "-targetMask", &p.targetMaskFileName) ;

  if (getCmdLineBool (&argc, argv, "-fastq33"))
    { p.fastQbase = 33 ; p.fastQ = TRUE ; }
  else if (getCmdLineBool (&argc, argv, "-fastq64"))
    { p.fastQbase = 64 ; p.fastQ = TRUE ; }

  p.probeQquality = -1000 ; /* do not look at the qualities (beware they may start at '@' or at '!' */
  if (getCmdLineInt (&argc, argv, "-minQuality", &p.probeQquality) && ! p.fastQbase) 
    usage ("-minQuality requires -fastq33 or -fastq64 because there are several conflicting conventions to measure fastqQ qualities", argc, argv) ;

  getCmdLineOption (&argc, argv, "-t", &(p.targetFileName)) ;
  getCmdLineOption (&argc, argv, "-subject", &(p.targetFileName)) ;
  getCmdLineOption (&argc, argv, "-targets", &(p.multiTargetFileNames)) ;

  if (getCmdLineOption (&argc, argv, "-target2gene", &(p.target2geneFileName)))
    target2GeneParse (&p) ;

  getCmdLineOption (&argc, argv, "-wordFrequencyConstructTable", &(p.wordFrequencyConstructTable)) ;
  getCmdLineOption (&argc, argv, "-wf", &(p.wordFrequencyTable)) ;
  getCmdLineInt (&argc, argv, "-wfSize", &(p.wffSize)) ;

  getCmdLineOption (&argc, argv, "-target_class", &(p.target_class)) ;
  if (getCmdLineOption (&argc, argv, "-selectPreviousScore", &(p.previousScoreFileName)))
    p.selectPreviousScore = TRUE ;
  else
    getCmdLineOption (&argc, argv, "-previousScore", &(p.previousScoreFileName)) ;
  getCmdLineOption (&argc, argv, "-select", &(p.selectedProbeFileName)) ;
  getCmdLineOption (&argc, argv, "-reject", &(p.rejectedProbeFileName)) ;

  p.nReadPlus = p.nReadMinus = 0 ;
  p.vectorize =  getCmdLineBool (&argc, argv, "-vectorize") ;
  p.silent = getCmdLineBool (&argc, argv, "-silent") ;
  if (!p.silent) fprintf (stderr, "// %s start\n", timeShowNow()) ;
  p.strand = getCmdLineBool (&argc, argv, "-strand") ;
  p.antistrand = getCmdLineBool (&argc, argv, "-antistrand") ;
  p.antiprobe = getCmdLineBool (&argc, argv, "-antiprobe") ;
  p.getTm  = getCmdLineBool (&argc, argv, "-getTm") ;
  p.solid  = getCmdLineBool (&argc, argv, "-solid") ;
  p.roche = getCmdLineBool (&argc, argv, "-roche") ;
  p.pacbio  = getCmdLineBool (&argc, argv, "-pacbio") ;
  p.ultima  = getCmdLineBool (&argc, argv, "-ultima") ;
  p.nanopore = getCmdLineBool (&argc, argv, "-nanopore") ;

  p.solid |= getCmdLineBool (&argc, argv, "-solidC") ;
  /* it make no sense not to use the correction code */
  if (getCmdLineBool (&argc, argv, "-SclipPolyA"))
    p.clipPolyA = p.strandedTarget = TRUE ;
  else
    p.clipPolyA  = getCmdLineBool (&argc, argv, "-clipPolyA") ;
  if (getCmdLineBool (&argc, argv, "-SclipPolyT"))
    p.clipPolyT = p.strandedTarget = TRUE ;
  else
    p.clipPolyT  = getCmdLineBool (&argc, argv, "-clipPolyT") ;

  /* intense RNA-editing */
  p.A2G  = getCmdLineBool (&argc, argv, "-A2G") ;
  p.G2A  = getCmdLineBool (&argc, argv, "-G2A") ;
  p.T2C  = getCmdLineBool (&argc, argv, "-T2C") ;
  p.C2T  = getCmdLineBool (&argc, argv, "-C2T") ;
  p.A2C  = getCmdLineBool (&argc, argv, "-A2C") ;
  p.C2A  = getCmdLineBool (&argc, argv, "-C2A") ;
  p.T2G  = getCmdLineBool (&argc, argv, "-T2G") ;
  p.G2T  = getCmdLineBool (&argc, argv, "-G2T") ;
  p.C2G  = getCmdLineBool (&argc, argv, "-C2G") ;
  p.G2C  = getCmdLineBool (&argc, argv, "-G2C") ;
  p.A2T  = getCmdLineBool (&argc, argv, "-A2T") ;
  p.T2A  = getCmdLineBool (&argc, argv, "-T2A") ;


  p.splice  = getCmdLineBool (&argc, argv, "-splice") ;
  p.doubleSeed  = getCmdLineBool (&argc, argv, "-doubleSeed") ;
  if (getCmdLineInt (&argc, argv, "-stranded", &p.strandBonus))
    p.stranded = 1 ;
  if (p.strandBonus < 0)
    { p.strandBonus *= -1 ; p.stranded *= -1 ; }
  p.aceOutGenomic = getCmdLineBool (&argc, argv, "-aceOutGenomic") ;
  p.aceOutMrna = getCmdLineBool (&argc, argv, "-aceOutMrna") ;
  p.aceOutNM = getCmdLineBool (&argc, argv, "-aceOutNM") ;
  p.showTargetPrefix  = getCmdLineBool (&argc, argv, "-showTargetPrefix") ;
  p.showProbeSequence  = getCmdLineBool (&argc, argv, "-showSequence") ;
  p.showOverhang  = TRUE ;
  p.showOverhang  = getCmdLineBool (&argc, argv, "-showOverhang") ;
  p.MRNAH  = getCmdLineBool (&argc, argv, "-MRNAH") ;
  p.maxHit = 10 ;
  getCmdLineInt (&argc, argv, "-maxHit", &p.maxHit) ;
  getCmdLineInt (&argc, argv, "-errMin", &p.errMin) ;
  getCmdLineInt (&argc, argv, "-errMax", &p.errMax) ;
  getCmdLineInt (&argc, argv, "-errRateMax", &p.errRateMax) ;

  p.errCost = 8 ;
  p.seedShift = 5 ; p.seedOffset = 1 ; p.seedLength = 16 ;

  if (getCmdLineOption (&argc, argv, "-strategy", &ccp))
    {
      if (! strcasecmp (ccp, "Exome") ||
	  ! strcasecmp (ccp, "Genome")
	  )
	{
	  p.strategy = STRATEGY_GENOME ;  
	}
      else if (! strcasecmp (ccp, "RNA_seq"))
	{
	  p.strategy = STRATEGY_RNA_SEQ ;
	}
      else
	{
	  fprintf (stderr, "-strategy %s , should be Exome, Genome or RNA_seq, sorry\n", ccp) ;
	  exit (1) ;
	}
    }
  else if (p.solid)
    {
      fprintf (stderr, "In -solid case the -strategy option must be specified as Exome, Genome or RNA_seq\n") ;
      exit (1) ;
    }

  if (p.nanopore)
    { p.errCost =  4; p.seedShift = 2 ;  p.seedLength = 14 ; }
  if (p.pacbio)
    p.errCost =  4;
  if (0 && p.ultima)
    p.errCost =  8;
  getCmdLineInt (&argc, argv, "-errCost", &p.errCost) ;

  getCmdLineInt (&argc, argv, "-seedShift", &p.seedShift) ;
  if (! getCmdLineInt (&argc, argv, "-seedOffset", &p.seedOffset))
    getCmdLineInt (&argc, argv, "-seedOffSet", &p.seedOffset) ; /* frequent typo */
  if (p.seedOffset < 0)
    {
      fprintf (stderr, "-seedOffset should be positive, not %d\n", p.seedOffset) ;
      usage (commandBuf, argc, argv) ;
    }
  if (p.solid && p.seedOffset == 1)
    p.seedOffset = 3 ;

  getCmdLineInt (&argc, argv, "-targetBonus", &p.bonus) ;
  p.decoy = getCmdLineBool (&argc, argv, "-decoy");
  p.sam = getCmdLineBool (&argc, argv, "-sam") ;
  if (getCmdLineInt (&argc, argv, "-slBonus", &p.slBonus))
    slInit (&p) ;
  {
    const char *enV, *enV2, *exV, *exV2 ;
    if (getCmdLineOption (&argc, argv, "-entryAdaptor", &enV))
      entryAdaptorInit (&p, enV) ;
    if (getCmdLineOption (&argc, argv, "-entryAdaptor2", &enV2))
      entryAdaptor2Init (&p, enV2) ;
    if (getCmdLineOption (&argc, argv, "-exitAdaptor", &exV))
      exitAdaptorInit (&p, exV) ;
    if (getCmdLineOption (&argc, argv, "-exitAdaptor2", &exV2))
      exitAdaptor2Init (&p, exV2) ;
  }
  p.minAli = 24 ;
  getCmdLineInt (&argc, argv, "-minAli", &p.minAli) ;
  getCmdLineInt (&argc, argv, "-minAliPerCent", &p.minAliPerCent) ;
  getCmdLineInt (&argc, argv, "-clipN", &p.clipN) ;

  p.gzi = getCmdLineBool (&argc, argv, "-gzi");
  p.gzo = getCmdLineBool (&argc, argv, "-gzo") ;

  getCmdLineInt (&argc, argv, "-overhangLength", &p.overhangLength) ;
  p.intronMaxLength = 120000 ; p.intronMinLength = 30 ;
  if (p.strategy != STRATEGY_RNA_SEQ)
    p.intronMinLength = 1000000000 ;

  getCmdLineInt (&argc, argv, "-intronMaxLength", &p.intronMaxLength) ;
  getCmdLineInt (&argc, argv, "-intronMinLength", &p.intronMinLength) ;

  p.maxDiscardableAli = 40 ;  /* do not discard above 70 2015_10_30, was 130 2018_03_15 */
  getCmdLineInt (&argc, argv, "-maxDiscardableAli", &p.maxDiscardableAli) ;

  getCmdLineInt (&argc, argv, "-mutantMalus", &p.mutantMalus) ;
  p.mutantMalus = 0 ; /* this idea is stupid, i realised that while programming it, we lose wild-type == mutant cases */
  p.avoidPseudoGenes = getCmdLineBool (&argc, argv, "-avoidPseudoGenes") ;

  getCmdLineInt (&argc, argv, "-UMI", &p.UMI) ;
  if (p.UMI > 0) 
    p.showTargetPrefix = FALSE ;
  else
    p.UMI = 0 ;
  getCmdLineInt (&argc, argv, "-jump5", &p.jump5_1) ;
  getCmdLineInt (&argc, argv, "-jump5_2", &p.jump5_2) ; 
  getCmdLineInt (&argc, argv, "-forceRightClip", &p.forceRightClip_1) ;
  getCmdLineInt (&argc, argv, "-forceRightClip_2", &p.forceRightClip_2) ; 

  getCmdLineInt (&argc, argv, "-clipAt", &p.clipAt) ;
  getCmdLineInt (&argc, argv, "-targetStart", &p.targetStart) ;
  getCmdLineInt (&argc, argv, "-targetStop", &p.targetStop) ;
  getCmdLineInt (&argc, argv, "-exactTargetStart", &p.exactTargetStart) ;
  getCmdLineInt (&argc, argv, "-exactTargetStop", &p.exactTargetStop) ;
  getCmdLineInt (&argc, argv, "-minEntropy", &p.minEntropy) ;
  if (p.minEntropy < 1)  p.minEntropy = 1 ;
  getCmdLineFloat (&argc, argv, "-minTM", &p.minTM) ;
  getCmdLineInt (&argc, argv, "-probeMinMultiplicity", &p.probeMinMultiplicity) ;
  getCmdLineInt (&argc, argv, "-probeMinLength", &p.probeMinLength) ;
  getCmdLineInt (&argc, argv, "-probeMaxLength", &p.probeMaxLength) ;

  getCmdLineOption (&argc, argv, "-run", &p.runName) ;
  getCmdLineOption (&argc, argv, "-method", &p.method) ;
  getCmdLineLong (&argc, argv, "-nRawReads", &p.nRawReads) ;
  getCmdLineLong (&argc, argv, "-nRawBases", &p.nRawBases) ;

  getCmdLineInt (&argc, argv, "-decimate", &p.decimate) ;
  if (! p.decimate && (ccp =getenv ("DECIMATE")))
    p.decimate = atoi (ccp) ;
  if (p.decimate < 0)
     {
       fprintf (stderr, "##### ERROR  -decimate must be positive :  %d is not allowed, soory\n", p.decimate) ;
       exit (1) ;
     }

  p.doubleGenes = bitSetCreate (50000, p.h) ;

  if (p.mutantMalus)
    p.originalScores = arrayHandleCreate (1000000, short, h) ;
  /* all this no longer applies now that we have read-pairs */ 
  if (0 && p.antiprobe && p.clipPolyA)
    {
      fprintf (stderr, "##### ERROR -antiprobe and -clipPolyA are incompatible, sorry\n") ;
      exit (1) ;
    }
  if (0 && p.stranded > 0 &&  p.clipPolyT)
    {
      fprintf (stderr, "##### ERROR -clipPolyT -stranded \"n>0\"  are incompatible options, sorry\n") ;
      exit (1) ;
    }
  if (0 && p.stranded < 0 &&  p.clipPolyA)
    {
      fprintf (stderr, "##### ERROR -clipPolyA -stranded \"n<0\" are incompatible options, sorry\n") ;
      exit (1) ;
    }
  if (getCmdLineInt (&argc, argv, "-seedLength", &p.seedLengthRequired) &&
      (p.seedLengthRequired < 7 || p.seedLengthRequired > 200)
      )
    {
      fprintf (stderr, "-seedLength should be be in the range [7,200], not %d\n", p.seedLengthRequired) ;
      exit (1) ;
    } 
  else
    p.seedLength = p.seedLengthRequired ;


  p.bestHit = TRUE ;
  p.bestHit =  getCmdLineBool (&argc, argv, "-best") ; /* efault, but we need to recognize -best on command line */
  p.bestHit = ! getCmdLineBool (&argc, argv, "-nonBest") ;

  if (p.probeMinLength < p.seedLength + p.seedOffset - 1)
    p.probeMinLength = p.seedLength + p.seedOffset - 1 ;
  if (p.probeMinLength < p.minAli)
    p.probeMinLength = p.minAli ;
  if (! p.minEntropy)
    p.minEntropy = p.seedLength/2 ;
  if (! p.minEntropy)
    p.minEntropy = 1 ;
  if (p.splice)
    p.bestHit = FALSE ;
  if (p.splice && ! p.minAli)
    p.minAli = 1 ;
  if (p.errCost <= 0)
    p.errCost = 8 ;
  p.endBonus = 1 && p.errCost > 2 ? p.errCost/2 - 1 : 0 ;
  p.endBonus = 0 ;

  if (p.minAli) 
    { 
      p.errMax = 1000 ; 
      if(p.minAli == 1) 
	p.minAli = p.seedLength; 
      if(p.minAli == 0) 
	p.minAli = p.seedLength = 16 ;
      if (!p.overhangLength)
	p.overhangLength = 8 ;  
      p.minChainAli = p.minAli ;
    }
  if (p.solid) 
    aceDnaSetSolidJumper (TRUE) ;
  else if (p.roche)  
    aceDnaSetRocheJumper (TRUE) ;
  else if (p.pacbio)  
    aceDnaSetPacBioJumper (TRUE) ;
  else if (p.nanopore) 
    aceDnaSetNanoporeJumper (TRUE) ;
  else
    aceDnaSetIlmJumper (TRUE) ;

  if (getCmdLineBool (&argc, argv, "--hits_file_caption"))
    {
      p.silent = FALSE ;
      p.ao = aceOutCreate (p.outFileName, ".hits_table_caption", 0, h) ;
      clipAlignTableCaption (&p) ;
      p.silent = TRUE ;
      ac_free (h) ;
      exit (0) ;
    }

  if (argc != 1)
    {
      fprintf (stderr, "unknown argument, sorry\n") ;
      usage (commandBuf, argc, argv) ;
    }
  if (p.targetFileName)
    fprintf (stderr, "// target %s\n", p.targetFileName) ;
  /* check the absolute args */
  if (! p.errMax && p.errRateMax > 0)
    p.errMax = 1000 ;
  clipAlignInitB2 (&p, FALSE) ;
  if (p.wordFrequencyConstructTable)
    { wordFrequencyConstructTable (&p) ; return 0 ; }

  if (!p.probeFileName || 
      (!p.targetFileName && ! p.multiTargetFileNames && !p.getTm)
      )
    usage (commandBuf, argc, argv) ;
  if (!clipAlignOpenFiles (&p))
    usage (commandBuf, argc, argv) ;
  
  if (! p.targetDict) /* may have been created by target2Geneparse */
    p.targetDict = dictHandleCreate (100, h) ;
  if (! p.isSplicedTarget)
    p.isSplicedTarget = bitSetCreate (8000, h) ;
  p.scannedGenes = bitSetCreate (8000, h) ;

  if (p.targetMaskFileName)
    parseTargetMaskFile (&p) ;

  if (p.wordFrequencyTable)
    wordFrequencyParseTable (&p) ;

  if (p.rejectedProbeFileName)
    {
      ACEIN ai = aceInCreate (p.rejectedProbeFileName, FALSE, h) ;
      if (ai)
	{
	  p.rejectedProbeDict = dictHandleCreate (10000, h) ;
	  while (aceInCard (ai))
	    if ((cp = aceInWord (ai)))
	      {
		if (!strcasecmp (cp, "probe"))
		  cp = aceInWord (ai) ;
		if (cp)
		  dictAdd (p.rejectedProbeDict, cp, 0) ;
	      }
	  ac_free (ai) ;
	}
      else
	exit (1) ;
    }
  if (p.selectedProbeFileName)
    {
      ACEIN ai = aceInCreate (p.selectedProbeFileName, FALSE, h) ;
      if (ai)
	{
	  p.selectedProbeDict = dictHandleCreate (10000, h) ;
	  while (aceInCard (ai))
	    if ((cp = aceInWord (ai)))
	      {
		if (!strcasecmp (cp, "probe"))
		  cp = aceInWord (ai) ;
		if (cp)
		  dictAdd (p.selectedProbeDict, cp, 0) ;
	      }
	  ac_free (ai) ;
	}
      else
	exit (1) ;
    }

  if (p.previousScoreFileName)
    {
      AC_HANDLE h1 = ac_new_handle () ;
      ACEIN ai= aceInCreate (p.previousScoreFileName, FALSE, h1) ;
      int n, x ;

      if (ai)
	{
	  p.previousScoreDict = dictHandleCreate (10000, h) ;
	  p.previousScores = arrayHandleCreate (10000, int, h) ;
	  while (aceInCard (ai))
	    if ((cp = aceInWord (ai)))
	      {
		if (dictAdd (p.previousScoreDict, cp, &n))
		  array(p.previousScores, n, int) = -99999 ;
		aceInStep (ai, '\t') ;
		if (aceInInt (ai, &x) && x > array(p.previousScores, n, int))
		  array(p.previousScores, n, int) = x ;
	      }
	  ac_free (ai) ;
	}
      else
	messcrash ("cannot open -previousScore file", p.previousScoreFileName) ;
      ac_free (h1) ;
    }
  
  np0 = clipAlignGetProbes (&p, 0) ;
  np1 = clipAlignGetProbes (&p, 1) ;
  p.nReads = arrayMax (p.probes) ;
  ac_free (p.probeDict) ; /* no longer needed */
  ac_free (p.probePairs) ; /* no longer needed */

  if (LATESORT && np0 + np1)
    clipAlignSortProbes (&p, np0, np1) ; /* sort alphabetically and establish 'same' and 'pair correspondance' */

  if (!p.preScanGenome && ! (np0+np1 > 0))
    {
      if (p.bonus < 0)
	{
	  fprintf (stderr, 
		   "// Given the negative bonus %d, no sequence was selected\n"
		   , p.bonus
		   ) ;
	  exit (0) ;
	}
	fprintf (stderr, 
		 "// Sorry, I do not understand the read fasta file,\n"
		 "// please check the name and format\n"
		 "// I expect something like\n"
		 ">n.1\n"
		 "atgtgtccagtagctatattagttc\n"
		 ">n.2\n"
		 "CCTGTGCNNCAT\n"
		 "AATCGCTAAGTGTGCA\n"
		 "//or for paired end sequences\n"
		 ">n.3\n"
		 "ATGTGCGTAGTC><TTTGCCAG\n"
		 "// upper lower case can be mixed. the 16 letter dna alphabet is accepted\n"
		 "// We use the standard UPAC coding for representing DNA with a single character\n"
		 "// per base:\n"
		 "// \n"
		 "// Exactly known bases\n"
		 "// \n"
		 "// A\n"
		 "// T  (or U for RNA)\n"
		 "// G\n"
		 "// C\n"
		 "// \n"
		 "// Double ambiguity\n"
		 "// \n"
		 "// R	AG	~Y		puRine\n"
		 "// Y	CT	~R		pYrimidine\n"
		 "// M	AC	~K		aMino\n"
		 "// K	GT	~M		Keto\n"
		 "// S	CG	~S		Strong\n"
		 "// W	AT	~w		Weak\n"
		 "// \n"
		 "// Triple ambiguity\n"
		 "// \n"
		 "// H	AGT	~D		not G	\n"
		 "// B	CGT	~V		not A\n"
		 "// V	ACG	~B		not T\n"
		 "// D	AGT	~H		not C\n"
		 "// \n"
		 "// Total ambiguity\n"
		 "// \n"
		 "// N	ACGT	~N		unkNown\n"
		 "// \n"
		 "// Run without parameter to get help\n"
		 ) ;
      exit (1) ;
    }
  ac_free (p.previousScores) ;
  ac_free (p.previousScoreDict) ;

  p.ao = aceOutCreate (p.outFileName, ".hits", p.gzo, h) ;

  if (p.antiprobe)
    {
      BOOL x = p.strand ;
      p.strand = p.antistrand ;
      p.antistrand = x ;
    }

  if (p.preScanGenome)
    clipAlignPreScanGenome (&p) ;
  else if (p.preScanProbes)
    clipAlignPreScanProbes (&p) ;
#ifdef JUNK
  else if (p.golay)
    golayHashProbes (&p) ;
#endif
  else if (p.getTm)
    {
      clipAlignExportTm (&p) ;
    }
  else
    {
      p.isShortTarget = clipAlignIsShortTarget (&p) ;

      p.exportHits = bigArrayHandleCreate (1000000, PEXPORT, h) ;
      p.exportDonors = bigArrayHandleCreate (1000000, PEXPORT, h) ;
      p.exportAcceptors = bigArrayHandleCreate (1000000, PEXPORT, h) ;
      p.exportIntrons = bigArrayHandleCreate (1000, PEXPORT, h) ;
      p.exportDoubleIntrons = bigArrayHandleCreate (1000, PEXPORT, h) ;

      p.exportGeneHits = bigArrayHandleCreate (1000, PEXPORT, h) ;

      bigArrayMax (p.exportHits) = 1 ; bigArrayMax (p.exportGeneHits) = 1 ;
      p.exportDict = dictCaseSensitiveHandleCreate (1000, h) ;

      if (1) /* optimize the seed length */
	{
	  int nn = p.seedLength ;
	  if (!nn) 
	    {
	      nn = p.errMax ? 16 : 32 ;
	      if (is64 && nn > 32) nn = 32 ;
	      if (!is64 && nn > 16) nn = 16 ;
	      
	      if (nn > p.probeLengthMin)
		nn = p.probeLengthMin ;
	    }	  
	  p.seedLength = nn ;
	}
      if (p.seedLength < 7)
	messcrash ("Please specify a seedLength or a probeMinLength >= 7") ;

      if (p.solid) /* sorry, we cannot do this in solid */
	p.A2G = p.G2A = p.T2C = p.C2T = p.A2C = p.C2A = p.T2G = p.G2T = p.G2C = p.C2G = p.A2T = p.T2A = p.X2X = FALSE ;
      
      if (p.A2G || p.G2A || 
	  p.T2C || p.C2T ||
	  p.A2C || p.C2A ||
	  p.T2G || p.G2T ||
	  p.G2C || p.C2G ||
	  p.A2T || p.T2A
	  )
	{
	  p.X2X = TRUE ; if (1) p.seedLength =  p.seedLength * log(4)/log(3) + .8 ; 
	  if (p.solid)
	    messcrash ("Sorry, A2G is not handled for the SOLiD case") ;
	}

      if (1 && p.minAli >  p.seedLength + 2 && p.splice)
	p.minAli = p.seedLength + 2 ; 	

      /* run with zero errors, unless errMin>0  */
      if (0 || (!p.minAli && !p.errMin && p.seedLength <= p.probeLengthMin))
	{
	  CLIPALIGN p2 ;
	  int np, nq, nq1 ;
	  
	  p2 = p ; /* save the required parameters */
	  /*
	    p2.errMax = 0 ;
	    p2.errRateMax = 0 ;
	  */

	  /* p2.seedLength = p2.probeLengthMin ; */
	  if (is64 && p2.seedLength > 32) p2.seedLength = 32 ;
	  if (!is64 && p2.seedLength > 16) p2.seedLength = 16 ;
	  if (0) p2.seedLength = is64 ? 32 : 16 ;

	  for (p2.hashPhase = 0 ; p2.hashPhase < 2 ; p2.hashPhase++)
	    {
	      if (1) 
		{
		  int nn, jx ;

		  nn = p.seedLength ;
		  p2.mask1 = 0 ;
		  if (nn < (is64 ? 32 : 16))
		    { p2.mask = p2.mask2 = (((unsigned long int)1) << (2 * nn)) - 1 ; jx = nn ; }
		  else if (is64)
		    { p2.mask = p2.mask2 = 0xffffffffffffffff ; jx = 32 ; }
		  else 
		    { p2.mask = p2.mask2 = 0xffffffff ; jx = 16 ; }
		  p2.gap = (p2.probeLengthMin - jx)/2 ;
		  if (p2.seedOffset > 0)
		    p2.gap = p2.seedOffset - 1 ;
		  
		  p2.dx1 = 0 ;
		  p2.dx2 = jx ;
		  
		  if (! p2.minEntropy)
		    p2.minEntropy = (jx)/2 ;
		  if (! p2.minEntropy)
		    p2.minEntropy = 1 ;
		}


	      fprintf (stderr, "// %s, hashing phase %d seedLenght=%d\n", timeShowNow(), p2.hashPhase, p2.seedLength) ;
	      fprintf (stderr, " // %ld assInserted, %ld assRemoved, %ld assFound, %ld assNotFound, %ld assBounce, %ld assBucketMax\n",
		       assInserted, assRemoved, assFound, assNotFound, assBounce, assBucketMax) ;

	      
	      np = clipAlignHashProbes (&p2) ;
	      nq1 =  p2.ass1 ? (unsigned int)p2.ass1->n + (p2.assB1 ? (unsigned int)p2.assB1->n : 0) : 0 ;
	      nq =  p2.ass ? (unsigned int)p2.ass->n + ( p2.assB ? (unsigned int)p2.assB->n : 0) : 0 ;
	      fprintf (stderr, "// %s, hashing done %d tags %d/%d words %d rejected nbits %d/%d\n"
		       , timeShowNow(), np, nq1,nq,nHashRejected,p2.ass->nbits,p2.ass1->nbits) ;
	      fprintf (stderr, " // %ld assInserted, %ld assRemoved, %ld assFound, %ld assNotFound, %ld assBounce, %ld assBucketMax\n",
		       assInserted, assRemoved, assFound, assNotFound, assBounce, assBucketMax) ;
	      if (p2.hashPhase)
		{
		  messfree (p2.ass1) ;
		  messfree (p2.oligoFrequency) ;
		}
	      nHits = clipAlignRun (&p2) ;	
	      if (p2.hashPhase == 0)
		{
		  unsigned char c , *cp = arrp (p2.oligoFrequency, 1, unsigned char) ;
		  int i = arrayMax (p2.oligoFrequency) - 1, nw = i, n0 = 0, n12 = 0, n25 = 0 ;
		  
		  if (i > 0) while (c = *cp++, i--)
		    {
		      if (c == 0) n0++ ;
		      else if (c >= 25) n25++ ;
		      else n12++ ;
		    }
		  fprintf (stderr, "// %s: Evaluated %d distinct words, %d ok, %d not in target, %d seen >= 25 times %.1f%% ok\n"
			   , timeShowNow (), nw, n12, n0, n25, (100.0*(nw-n0-n25))/(nw ? nw : 1)
			   ) ;
		}
	      fprintf (stderr, "// %s: search done, %d tags %d hits %d exported verified=%d perfect=%d\n", timeShowNow(), np, nHits,  p2.nExportedHits, nVerif, nPerfect) ;
	      fprintf (stderr, " // %ld assInserted, %ld assRemoved, %ld assFound, %ld assNotFound, %ld assBounce, %ld assBucketMax\n",
		       assInserted, assRemoved, assFound, assNotFound, assBounce, assBucketMax) ;
	    }
	  p.seedOffsetUsed = p2.seedOffsetUsed ;
	  p.nExportedHits = p2.nExportedHits ;
	  /* kill the exact-probe hash and reopen the target file */
	  assDestroy (p2.ass) ;
	  assDestroy (p2.ass1) ;
	  assDestroy (p2.ass2) ;
	  messfree (p2.oligoFrequency) ;
	}
      
     /* run again with the required number of errors */
      if (1 && (p.errMax || p.errRateMax || p.minAli))
	{
	  int jx, nn, np = 1, nq = 0, nq1 = 0 ;
	  
	  if (! is64)
	    {
	      if (p.probeLengthMax >= (1<<10))
		messcrash ("Sorry, on a 32 bit machine i cannot align a sequence longer than 1000 bp, please clip it") ;
	      if (arrayMax (p.probes) >= (1<<22))
		messcrash ("Sorry, on a 32 bit machine i cannot align more than 4Million sequences at once. Please split in smaller blocks") ;
	    }
      
	  nn = p.seedLength ;
	  if (1) 
	    {
	      p.mask1 = 0 ;
	      if (nn < (is64 ? 32 : 16))
		{ p.mask = p.mask2 = (((unsigned long int)1) << (2 * nn)) - 1 ; jx = nn ; }
	      else if (is64)
		{ p.mask = p.mask2 = 0xffffffffffffffff ; jx = 32 ; }
	      else 
		{ p.mask = p.mask2 = 0xffffffff ; jx = 16 ; }
	      p.gap = (p.probeLengthMin - jx)/2 ;
	      if (p.seedOffset > 0)
		p.gap = p.seedOffset - 1 ;

	      p.dx1 = 0 ;
	      p.dx2 = jx ;

	      if (! p.minEntropy)
		p.minEntropy = (jx)/2 ;
	      if (! p.minEntropy)
		p.minEntropy = 1 ;
	    }
	  
	  for (p.hashPhase = 0 ; p.hashPhase < 2 ; p.hashPhase++)
	    {
	      fprintf (stderr, "// %s, hashing phase %d seedLenght=%d\n", timeShowNow(), p.hashPhase, p.seedLength) ;
	      fprintf (stderr, " // %ld assInserted, %ld assRemoved, %ld assFound, %ld assNotFound, %ld assBounce, %ld assBucketMax\n",
		       assInserted, assRemoved, assFound, assNotFound, assBounce, assBucketMax) ;
	      if (p.isShortTarget && p.hashPhase == 0)
		{
		  clipAlignHashTarget (&p) ;
		}
	      else
		{
		  np = clipAlignHashProbes (&p) ; /* number of long enough probes */
		  nq1 =  p.ass1 ? (unsigned int)p.ass1->n + (p.assB1 ? (unsigned int)p.assB1->n : 0) : 0 ;
		  nq =  p.ass ? (unsigned int)p.ass->n  + (p.assB ? (unsigned int)p.assB->n : 0)  : 0 ;
		}
	      fprintf (stderr, "// %s, hashing done %d tags %d/%d words %d rejected ...nbits %d/%d\n", timeShowNow(), np, nq1,nq,nHashRejected,p.ass->nbits,p.ass1->nbits) ;
	      fprintf (stderr, " // %ld assInserted, %ld assRemoved, %ld assFound, %ld assNotFound, %ld assBounce, %ld assBucketMax\n",
		       assInserted, assRemoved, assFound, assNotFound, assBounce, assBucketMax) ;
	      if (p.hashPhase)
		{
		  messfree (p.ass1) ;
		  messfree (p.oligoFrequency) ;
		}
	      nHits = clipAlignRun (&p) ; 
	      if (p.hashPhase == 0)
		{
		  unsigned char c , *cp = 0 ;
		  int i = arrayMax (p.oligoFrequency) - 1, nw = p.nOligo, n0 = 0, n12 = 0, n25 = 0 ;
		  
		  if (arrayMax (p.oligoFrequency) > 1)
		    cp = arrp (p.oligoFrequency, 1, unsigned char) ;
		  if (arrayMax (p.oligoFrequency) > 1)
		    while (c = *cp++, i--)
		      {
			if (c == 0) n0++ ;
			else if (c >= 25 && c > p.maxHit) n25++ ;
			else n12++ ;
		      }
		  fprintf (stderr, "// Evaluated %d distinct words, %d ok, %d not in target, %d seen >= %d times %.1f%% ok\n"
			   , nw, n12, n0, n25,  p.maxHit < 25 ? 25 :  p.maxHit, (100.0*(nw-n0-n25))/(nw ? nw : 1)
			   ) ;
		}
	      fprintf (stderr, "// %s search done, %d tags %d hits %d exported verified=%d perfect=%d\n"
		       , timeShowNow(), np, nHits, p.nExportedHits, nVerif, nPerfect) ;
	      fprintf (stderr, " // %ld assInserted, %ld assRemoved, %ld assFound, %ld assNotFound, %ld assBounce, %ld assBucketMax\n",
		       assInserted, assRemoved, assFound, assNotFound, assBounce, assBucketMax) ;
	    }
	}
      clipAlignExportAll (&p, commandBuf) ;
    }

  if (1 && p.antiprobe)
    {
      BOOL x = p.strand ;
      p.strand = p.antistrand ;
      p.antistrand = x ;
    }


  ac_free (p.ao) ;
  ac_free (p.h) ;
  if (!p.silent)
    { 
      int mx ;
      messAllocMaxStatus (&mx) ; 
      fprintf (stderr, "// %s done: found %d Hits, nVerif=%d, perfect=%d, max memory %d Mb\n", timeShowNow(), nHits, nVerif, nPerfect, mx) ;
      fprintf (stderr, "// errCost = %d, errMax = %d, errRateMax = %d, minTM = %.2f, minEntropy = %d, maxHit = %d, minAli = %d,  minAliPerCent = %d, overhangLength=%d\n// seedLength required = %d/%d, seedOffset %d/%d, seedShift=%d\n// single-strand = %s, antistrand = %s, antiprobe=%s\n// Required: probeMinLength = %d, probeMaxLength = %d, probeMinMultiplicity = %d,\n// observed: probeLengthMin = %d, probeLengthMax = %d\n"
	       , p.errCost, p.errMax, p.errRateMax
	       , p.minTM, p.minEntropy, p.maxHit, p.minAli, p.minAliPerCent, p.overhangLength
	       , p.seedLengthRequired, p.seedLength
	       , p.seedOffset, p.seedOffsetUsed, p.seedShift
	       , p.strand ? "TRUE" : "FALSE"
	       , p.antistrand ? "TRUE" : "FALSE"
	       , p.antiprobe ? "TRUE" : "FALSE"
	       , p.probeMinLength , p.probeMaxLength, p.probeMinMultiplicity
	       , p.probeLengthMin
	       , p.probeLengthMax
	       ) ;
    }
  fprintf(stderr, "// nReadPlus=%d nReadMinus=%d\n", p.nReadPlus, p.nReadMinus) ;
  fprintf(stderr, "// nIntronPlus=%d nIntronMinus=%d\n", p.nIntronPlus, p.nIntronMinus) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr and to ensure all pipes are closed*/
  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

/*
 * 2008_02_10
 * strategic idea

 * consider a probe of length at least 32
 * i have 7 start positions for a 26mer: 1 2 3 4 5 6 7
 * there are 7*6/2=21 paits of doublets, 21 - 6 = 15 non overlapping
 * but there are only 10 different pairs of letters
 * so in every probe at least one knind of pair is present twice
 * we sort the probes in bins, putting them in the bin where they have
 * most of the rarer doublet

 * we preprocess the genome and hash the various pairs in different
 * preporcessing, this CPU time does not count

 * present the probes to the preprocessed genome of the correct pair
 * if it has 3 doublets we are garanteed to see it, but only if it is exact
 * since the 24 mers overlap !

 * we present the probe then the complemented probe, so we explore the 2 strands at once

 * we could do the same trick with 3 letters, we hash 1 pos of the genome in 64
 * so 150M words hashing uses 300M * 8bytes = 6 GB, too large

 * we could do the same trick with 4 letters, we hash 1 pos of the genome in 256
 * so 40M words hashing uses 80M * 8bytes = 640 Mbytes, ok

 * phase 1:
    count all 4 mers in the genome, sort them by decreasing frequency
    prepare 256 preprocessed genomes with a known 4bp start
    we hash the next 22bp and prepare a binary file representing the hashing table
      exact 22bp/position in words 0f 64 bits
      position has 32 bit which has to be remapped to chrom+pos, bad luck !
      better idea, 22bp use 44 bits, use last 16 bits for chrom number !
      then over 64 bits we give the next 32bp on the genome
 
 * phase 2:
    in each probe find the most significant 4 mer with 22 bp behind
    put the probe in that bin
    probe_number = 64 bits +  pos in probe = 8 bits + 22bp sequence = 44 bits
    then over 64 bits we give the next 32bp on the probe

 * phase 3:
    run the alignment code, get all exact hits, if none export 1 err hits
*/

  
