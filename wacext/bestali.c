/*  File: bestali.c
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
 *  available from http://www.aceview.org
 *
 *  Please direct all question about this code to
 *  mieg@ncbi.nlm.nih.gov
 */

#include "ac.h"
#include "acedb.h"
#include "bitset.h"

typedef struct baStruct { 
  AC_HANDLE h ; 
  AC_DB db ;           /* ACEDB database containing the semantics of the experiment */
  ACEOUT maxHitOut ;
  BitSet maxHitBitSet ;
  const char *outFileName ;
  const char *hitFileName ;
  const char *inFileList ;
  const char *run, *title ;
  const char *mrnaRemapFileName ;
  const char *uniqueGeneSupportFileName ;
  Array hits, hits2, hits3 ;
  Array uniqueGeneSupport ;
  Array gene2map ;
  DICT *tagDict ;
  DICT *cloneDict ;
  DICT *errDict ;
  DICT *targetDict ;
  DICT *runDict ;
  DICT *geneDict ;
  DICT *intronDict ;
  DICT *exonDict ;
  DICT *target_classDict ;
  DICT *dnaDict ;
  DICT *inFileDict ; 
  const char *selected8kbFileName ;
  const char *selected5kbFileName ;
  const char *splitMrnaFileName ;
  KEYSET selected8kbWiggle ;
  KEYSET selected5kbWiggle ;
  DICT *selected8kbDict ;
  DICT *selected5kbDict ;
  KEYSET sigTargetKs ;
  BOOL gzi, gzo ;
  BOOL checkGroupHierearchy ;
  BOOL countBest, exportBest ; /* export the hierarchic count */
  BOOL exportVenn ; /* complete venn diagram of number of reads in all targets combinations */
  BOOL exportSuffix ;
  BOOL keepTargetPrefix ;
  BOOL exportMito ;
  BOOL aliProfile ;
  BOOL seqc ;
  BOOL errorProfile ;
  BOOL intronSupport ;
  BOOL mrnaSupport, hierarchic ;
  BOOL geneSupport ;
  BOOL intergenicSupport ;
  BOOL geneSupport2ace ;
  BOOL mrnaSupport2ace ;
  BOOL intronSupport2ace ;
  BOOL geneIndex2Table ;
  BOOL groupLetterProfile ;
  BOOL autoH ;
  BOOL addQualityFactors ;
  BOOL RNA_seq ;
  BOOL hasPair ; /* pair info found in input hits file */ 
  BOOL SAM ;
  int pair ;
  int lastTranscriptClass ;
  /* filtering */
  int Remove_inserts_shorter_than ;
  const char *filter ; /* apply a named quality filter */
  int  pureNsStrand ; /* what to do with non stranded hits */
  int  stranded ;     /* stranded or antistranded manip = +-1*/
  BOOL strand ;     /* just accept alignments on the forward strand of the target */
  BOOL antistrand ; /* just accept alignments on the reverse strand of the target */
  BOOL unique ;     /* only consider tags hitting a single target */
  int target_class ;   /* entry in targetDict, A_mito ... Z_genome */
  int intronClass ;
  int maxErr, maxErrRate ;
  int subsampling ;
  BOOL solid ;      /* auto-recognized on existence of oo errors, affects the filter */
  BOOL greg2ace ;

  Array genome2geneAtlas ;
  Array target2intronAtlas ;
  Array target2exonAtlas ;
  Array splitMrnaArray ;
  DICT *splitMrnaDict ;
  Associator splitMrnaAss ;
  int maxHit ;      /* only consider alignments with less than this number of targets */
  long int nS, nT ; /* total number of sequences and tags */
  int minAli, errCost ; /* global filtering on aligned length */
  int minIntronOverlap ; /* default 8: how many overlap bases to confirm an intron seen in a read aligned to a transcript */
  int ventilate ;   /* number of group of original id files */ 
  char *prefixId ;  /* prefix common to all the sequence identifiers */
  const char *project ;
  const char *id_filename ;
  const char *fastqAnalysis ;
  const char *lane ;
  const char* geneRemapFile ;
  const char *sigTargetFile ;
  ACEOUT seqco ;
} BA ;

/* the structure should be very small
   we have 5M reads, with may have easilly 20M alignmnets, with 80 bytes per struct this gives 80*25M = 2G
   and the array will overflow
 */

typedef struct hitStruct3 { 
  int errTag          /* entry in errDict, errors in tag coordinates and orientation */
    , errTarget       /* entry in errDict, errors in target coordinates and orientation */
    , prefix, suffix  /* overhangs, read in the tag, on the strand extending the alignment */
    , targetPrefix, targetSuffix  /* 30bp, up and downstream, read on the target in the orientation of the tag */
    , nN, nErr
    , c1, c2 /* chain start/stop in x coordinates */
    , gt_ag, ct_ac /* justified introns */
    , hardClip1, hardClip2 
    ;
} HIT3 ;

typedef struct hitStruct2 { 
  int  a1, a2
    , dPair           /* distance to pair */
    , ln              /* length to be aligned (minus the vectors and pA) */
    , ali             /* length of the alignment (measured on the tag) */
    , target_class    /* entry in targetDict, A_mito ... Z_genome */
    , mult            /* multiplicity of the tag */
    , unicity         /* number of targets in that class, negative if we change strand */
    , badPair ; /* raised only in paired end case, discard from gene-counts: bad topology:1 or length: 2 */
    ;
} HIT2 ;

typedef struct hitStruct { 
  int nn             /* links HIT to HIT2 and HIT3 */
    , tag             /* identifier of the short sequence tag */
    , gene            /* entry in geneDict */
    , target          /* entry in targetDict */
    , clone           /* common identifiers of all reads of a given molecule */
    , score         /* score of the ali */
    , x1, x2
    , class           /* the first letter of the target class, defining the hierarchy */
    , chain
    ;
} HIT ;

typedef struct mhitStruct { 
  int exon, gene             /* identifier of the short sequence tag */
    , a1, a2, x1, x2
    , class           /* the first letter of the target class, defining the hierarchy */
    ;
} MHIT ;


typedef struct splitMrnaStruct { 
  int gene, mrna, gXX, gNewOld
    , x1, x2
    ;
} SPLITMRNA ;


typedef struct autoStruct { 
  const char *run ;
  int target        /* entry in targetDict */
    , tag           /* identifier of the short sequence tag */
    , score
    ;
  long long int best          /* Best mapping reads, reattributed only to the target with most hits */
    , tags          /* Best mapping reads, including those mapping equally well to several targets */
    , seqs          /* distinct sequences */
    , ali           /* mapped bases */
    , ln            /* cumulated read length */ 
    ;
} AUTOHIT ;

static int badMaxHit, badLnS, badLnT, isLongS, isLongT, badScoreLS, badScoreLT, isPartialS, isPartialT, badScorePS, badScorePT ;
static int baHit2Exon (BA *ba, HIT *vp, int *exonp, int *genep) ;
static int baHit2Gene (BA *ba, HIT *vp, int *genep) ;
extern int atoi(const char *cp) ;
static BOOL baParseOneSamTranscriptHit (ACEIN ai, BA *ba, HIT *up, int nn) ;
static BOOL baKeepOneHit (BA *ba, HIT *up, const int nn, const int Z_genome) ;
static int baSplitMrnaRemap (BA *ba, HIT *up, HIT2 *up2) ;

/*************************************************************************************/

int intronHitOrder (const void *a, const void *b)
{
  const HIT *up = (const HIT *)a, *vp = (const HIT *)b ;
  int n ;

  n = up->clone - vp->clone ; if (n) return n ;
  n = up->tag - vp->tag ; if (n) return n ;
  n = up->x1 - vp->x1 ; if (n) return n ;
  n = up->x2 - vp->x2 ; if (n) return n ;
  n = up->class - vp->class ;  if (n) return n ;  
  n = up->gene - vp->gene ; if (n) return n ;  
  n = up->target - vp->target ; if (n) return n ;

  return 0 ;
} /* intronHitOrder */

/*************************************************************************************/

int mhitOrder (const void *a, const void *b)
{
  const MHIT *up = (const MHIT *)a, *vp = (const MHIT *)b ;
  int n ;

  n = up->x1 - vp->x1 ; if (n) return n ;
  n = up->x2 - vp->x2 ; if (n) return n ;
  n = up->exon - vp->exon ; if (n) return n ;
  n = up->gene - vp->gene ; if (n) return n ;
  return 0 ;
} /* mhitOrder */

/*************************************************************************************/

int mhitGeneOrder (const void *a, const void *b)
{
  const MHIT *up = (const MHIT *)a, *vp = (const MHIT *)b ;
  int n ;

  n = up->gene - vp->gene ; if (n) return n ;
  n = up->x1 - vp->x1 ; if (n) return n ;
  n = up->x2 - vp->x2 ; if (n) return n ;
  n = up->exon - vp->exon ; if (n) return n ;

  return 0 ;
} /* mhitGeneOrder */

/*************************************************************************************/
  /* 2013_08_01 
   * Ordering problem, it is necesary to sort alphabetically the 'target'
   * because the hierarchic mapping of the transcripts relies on their naming
   * so before sorting the aa array, we must remap the 'targets' ina sorted dictionary
   */
static DICT *baTargetDict = 0 ;
int baHitOrder (const void *a, const void *b)
{
  const HIT *up = (const HIT *)a, *vp = (const HIT *)b ;
  int n ;

  n = up->clone - vp->clone ; if (n) return n ;
  n = up->tag - vp->tag ; if (n) return n ;
  n = up->class - vp->class ; if (n) return n ;
  n = up->score - vp->score ; if (n) return -n ;
  n = up->gene - vp->gene ; if (n) return n ;
  if (up->target && vp->target)
    n = lexstrcmp(dictName(baTargetDict,up->target),dictName(baTargetDict, vp->target)) ; 
  if (n) return n ;
  n = up->target - vp->target ; if (n) return n ;
  n = up->x1 - vp->x1 ; if (n) return n ;
  n = up->x2 - vp->x2 ; if (n) return n ;

  return 0 ;
} /* baHitOrder */

/*************************************************************************************/

int baCloneOrder (const void *a, const void *b)
{
  const HIT *up = (const HIT *)a, *vp = (const HIT *)b ;
  int n ;

  n = up->clone - vp->clone ; if (n) return n ;
  n = up->tag - vp->tag ; if (n) return n ;
  n = up->class - vp->class ; if (n) return n ;
  n = up->score - vp->score ; if (n) return -n ;
  n = up->gene - vp->gene ; if (n) return n ;
  n = up->target - vp->target ; if (n) return n ;
  n = up->x1 - vp->x1 ; if (n) return n ;
  n = up->x2 - vp->x2 ; if (n) return n ;

  return 0 ;
} /* baCloneOrder */

/*************************************************************************************/
/*************************************************************************************/
static void usage (char *message) ;

/*************************************************************************************/
/*****************************    Utilities     **************************************/
/*************************************************************************************/
/*************************************************************************************/
/*****************************    Actual work   **************************************/
/*************************************************************************************/

static void baHitCompress (BA *ba, Array aa, int  Z_genome)
{
  int ii, jj, kk, aaMax = arrayMax (aa), tag ;
  HIT *up, *vp ;
  HIT2 *up2, *vp2 ;

  arrayCompress (aa) ;
  for (ii = kk = 0, up = arrp (aa, 0, HIT) ; ii < aaMax ; ii++, up++)
    {
      if (! baKeepOneHit (ba,up,ii, Z_genome)) /* will possibly zero the score */
	continue ;
      tag = up->tag ;
      up2 = arrp (ba->hits2,up->nn, HIT2) ;
      for (jj = ii + 1, vp = up + 1 ;
	   tag && jj < aaMax && vp->tag == tag && up->x1 == vp->x1 && up->x2 == vp->x2 && up->target == vp->target ;
	   jj++, vp++)
	{
	  vp2 = arrp (ba->hits2, vp->nn, HIT2) ;
   	  if (up2->a1 == vp2->a1 && up2->a2 == vp2->a2 && up2->target_class == vp2->target_class)
	    tag = up->tag = 0 ; 
	}
      if (tag)
	{
	  if (kk < ii) 
	    {
	      vp = arrp (aa, kk, HIT) ;
	      *vp = *up ;
	    }
	  kk++ ;
	}
    }
  arrayMax (aa) = kk ;
} /* baHitCompress */

/*************************************************************************************/

static void baExportOneHit (ACEOUT ao, BA *ba, HIT *up)
{ 
  HIT2 *up2 = arrp (ba->hits2,up->nn, HIT2) ;
  HIT3 *up3 = arrp (ba->hits3,up->nn, HIT3) ;

  aceOutf (ao, "%s\t%05d\t%d\t%d\t%d\t%d\t%d"	       
	   , up->tag ? dictName (ba->tagDict, up->tag) : "-"
	   , up->score
	   , up2->mult
	   , up2->ln
	   , up2->ali
	   , up->x1, up->x2
	   ) ;
  
  aceOutf (ao, "\t%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s"
	   , up2->target_class ? dictName (ba->target_classDict, up2->target_class) : "-"
	   , up->gene ? dictName (ba->targetDict, up->gene) : "-"
	   , up2->unicity /* unicity */
	   , up->target ? dictName (ba->targetDict, up->target) : "-"
	   , up2->a1, up2->a2
	   , up3->nN, up3->nErr
	   
	   , up3->errTag ? dictName (ba->errDict, up3->errTag)  : "-"    /* error list in tag coordinates and orientation */
	   , up3->errTarget ? dictName (ba->errDict, up3->errTarget)  : "-"     /* error list in target coordinates and orientation */
	   
	   , up3->prefix ? dictName (ba->dnaDict, up3->prefix) : "-"       /* 5' overhang read in the tag */
	   , up3->suffix ? dictName (ba->dnaDict, up3->suffix) : "-"       /* 3' overhang read in the tag */
	   
	   , ba->keepTargetPrefix && up3->targetPrefix ? dictName (ba->dnaDict, up3->targetPrefix) : "-"
	   , ba->keepTargetPrefix && up3->targetSuffix ? dictName (ba->dnaDict, up3->targetSuffix) : "-"
	   ) ;

  if (ba->pair)
    {
      if (class(up2->dPair) == 6)
	aceOutf (ao, "\t-6 links_to:%s", dictName (ba->targetDict, KEYKEY(up2->dPair))) ;
      else if (class(up2->dPair) == 10)
	aceOutf (ao, "\t-10 links_to:%s", dictName (ba->targetDict, KEYKEY(up2->dPair))) ;
      else if (class(up2->dPair) == 2)
	aceOutf (ao, "\t-2 close_links_to:%s", dictName (ba->targetDict, KEYKEY(up2->dPair))) ;
      else
	aceOutf (ao, "\t%d", up2->dPair) ;
    }
  else
    aceOutf (ao, "\t0") ;

  aceOutf (ao, "\t-\t-\t-") ;

  if (up->chain)
    aceOutf (ao, "\tchain %d\t%d\t%d\t%d\t%d"
	     , up->chain
	     , up3->c1, up3->c2
	     , up3->gt_ag, up3->ct_ac
	     ) ;
  else
    aceOut (ao, "\t") ;
  aceOut (ao, "\n") ;
} /* baExportOneHit */

/*************************************************************************************/

static BOOL baParseOneHit (ACEIN ai, BA *ba, HIT *up, int nn)
{ 
  char cutter, *ccp ;
  char buf[1000] ;
  static int nTP = 0, nTS = 0, oldTag = 0 ;
  int dummy, dummy2, Z_genome = 0 ;
  HIT2 *up2 ;
  HIT3 *up3 ;
  int isFirstFragment = 1 ;
  int deltaPair ;
  int chain = 0 ;

  aceInSpecial (ai, "\n") ;

  if (ba->SAM)
    return  baParseOneSamTranscriptHit (ai, ba, up, nn) ;

  memset (up, 0, sizeof (HIT)) ;
  if (ba->target_classDict) 
    dictFind (ba->target_classDict, "Z_genome", &Z_genome) ;
  strncpy (buf, aceInPos (ai), 999) ;
  ccp = aceInWordCut (ai, "\t", &cutter) ;
  if (! ccp || ! *ccp || *ccp == '#')
    return FALSE ;
  if (ba->subsampling)
    {
      static char oldClone[256] ;
      static int keep = 0 ;
      char *cr ;
      char newClone[256] ;

      strncpy (newClone, ccp, 255) ; /* remove <> in my own private buffer */
      newClone[255] = 0 ;
      for (cr = newClone ; *cr ; cr++) 
	;
      cr-- ;
      if (*cr == '>' || *cr == '<')
	*cr = 0 ;
      if (strcmp (newClone, oldClone))
	{ /* new clone */
	  float z = randfloat() ;
	  keep = ba->subsampling * z > 1 ? TRUE : FALSE ;
	  
	  strcpy (oldClone, newClone) ;
	}
      if (! keep)
	return FALSE ;
    }

  dictAdd (ba->tagDict,ccp, &(up->tag)) ;
  up->nn = nn ;
  up2 = arrayp (ba->hits2, up->nn, HIT2) ;
  up3 = arrayp (ba->hits3, up->nn, HIT3) ;
  memset (up2, 0, sizeof (HIT2)) ;
  memset (up3, 0, sizeof (HIT3)) ;
  isFirstFragment = 1 ;
  if (1 || ba->pair)
    {
      char cc, *ccq = ccp + strlen (ccp) - 1  ;
      BOOL ok = FALSE ;

      if (*ccq == '>')
	{ ok = TRUE ; isFirstFragment = 1 ; }
      else if (*ccq == '<')
	{ ok = TRUE ; isFirstFragment = -1 ;}
      else
	{
	  int n = 3 ;
	  while (n && *ccq != '/') { n-- ; ccq -- ; }
	  if (*ccq == '/')
	    { ok = TRUE ; }
	}
      if (ba->pair && ok) /* if no clone, do not enter it in the cloneDict */
	{
	  cc = *ccq ; *ccq = 0 ;
	  dictAdd (ba->cloneDict,ccp, &(up->clone)) ;
	  *ccq = cc ;
	}
    }

  {
    int score = 0, mult = 0, ln = 0, ali = 0, x1 = 0, x2 = 0 ;
    
    if (
	! aceInInt (ai, &(score)) ||
	! aceInStep (ai, '\t') || ! aceInInt (ai, &(mult)) ||
	! aceInStep (ai, '\t') || ! aceInInt (ai, &(ln)) || ln < 5 ||
	! aceInStep (ai, '\t') || ! aceInInt (ai, &(ali)) ||
	! aceInStep (ai, '\t') || ! aceInInt (ai, &(x1)) ||
	! aceInStep (ai, '\t') || ! aceInInt (ai, &(x2)) 
	)
      {
	fprintf(stderr, "cannot read x2 or ln=%d<0 in %s\n",  up2->ln, dictName (ba->tagDict,up->tag)) ;
	return FALSE ;
      }
    up->score = score ; up2->mult = mult ; up2->ln = ln ; up2->ali = ali ; up->x1 = x1 ; up->x2 = x2 ;
  }
  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ;
  if (! ccp || ! *ccp)
    {
      fprintf(stderr, "cannot read target_class\n") ;
      return FALSE ;
    }
  if (strcmp (ccp, "-"))
    {
      dictAdd (ba->target_classDict,ccp, &dummy) ;
      up2->target_class = dummy ; 
    }
  if (
      (! ba->intergenicSupport && ba->target_class && up2->target_class != ba->target_class) ||
      (
       ba->intergenicSupport && ba->target_class && up2->target_class != ba->target_class && 
       ! (Z_genome &&  up2->target_class == Z_genome && up->tag != oldTag)
       )
      )
    return FALSE ;
  if (up2->target_class == ba->target_class)
    oldTag = up->tag ;
  up->class = *ccp ;
  
   ccp = aceInWordCut (ai, "\t", &cutter) ;
  if (! ccp || ! *ccp)
    {
      fprintf(stderr, "cannot read gene\n") ;
      return FALSE ;
    }
  if (strcmp (ccp, "-"))
    dictAdd (ba->targetDict,ccp, &(up->gene)) ;

  if (! aceInInt (ai, &dummy))
    {
      fprintf(stderr, "cannot read uniticity\n") ;
      return FALSE ;
    }
  {
    int ngenes ;
    up2->unicity = ngenes = dummy ;
    if ((ba->stranded || ba->pureNsStrand == 2) && ngenes == -2) ngenes = 1 ;
    if (ngenes < 0 && ba->pureNsStrand) return FALSE ;
    if ((ngenes > 1 || ngenes < -1 ) && ba->unique) return FALSE ;
  }
  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ;
  if (1 && ccp && !strncmp (ccp,"MRNA:",5))ccp += 5 ;
  if (! ccp || ! *ccp)
    return FALSE ;
  dictAdd (ba->targetDict,ccp, &(up->target)) ;

  if ( ! aceInInt (ai, &(up2->a1)) ||
      ! aceInStep (ai, '\t') || ! aceInInt (ai, &(up2->a2)) ||
      ! aceInStep (ai, '\t') || ! aceInInt (ai, &dummy) ||
      ! aceInStep (ai, '\t') || ! aceInInt (ai, &dummy2)
      )
    {
      fprintf(stderr, "cannot read nErr\n") ;
      return FALSE ;
    }

  if (! (ba->geneSupport || ba->mrnaSupport) && ba->stranded * isFirstFragment * (up2->a2 - up2->a1) < 0) 
    return FALSE ;

  up3->nN = dummy ;
  up3->nErr = dummy2 ;	   
  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ;  /* col 16: errTag */
  if (! ccp || ! *ccp) 
    {
       fprintf(stderr, "cannot read errTag %s\n", dictName (ba->tagDict,up->tag)) ;
       return FALSE ;
    }
  if (strcmp (ccp, "-"))
    {
      dictAdd (ba->errDict,ccp, &(up3->errTag)) ;
      if (! ba->solid && strstr (ccp, "oo"))
	ba->solid = TRUE ;
    }
  
  ccp = aceInWordCut (ai, "\t", &cutter) ;  /* col 17: errTarget*/
  if (! ccp || ! *ccp)
    {
      fprintf(stderr, "cannot read errTarget\n") ;
      return FALSE ;
    }
  if (strcmp (ccp, "-"))
    dictAdd (ba->errDict,ccp, &(up3->errTarget)) ;
  
  ccp = aceInWordCut (ai, "\t", &cutter) ; /* col 18: prefix */
  if (! ccp || ! *ccp)
    {
      fprintf(stderr, "cannot read prefix\n") ;
      return FALSE ;
    }
  if (strcmp (ccp, "-"))
    {
      dictAdd (ba->dnaDict,ccp, &(up3->prefix)) ;
      if (*ccp == ace_upper(*ccp))
	up3->hardClip1 = up->x1 ;
    }
  
  ccp = aceInWordCut (ai, "\t", &cutter) ; /* col 19: suffix */
  if (! ccp || ! *ccp)
    {
      fprintf(stderr, "cannot read suffix\n") ;
      return FALSE ;
    }
  if (strcmp (ccp, "-"))
    {
      dictAdd (ba->dnaDict,ccp, &(up3->suffix)) ;
      if (*ccp == ace_upper(*ccp))
	up3->hardClip2 = up->x2 ;
    }
  
  ccp = aceInWordCut (ai, "\t", &cutter) ; /* col 20: targetPrefix */
  if (! ccp || ! *ccp)
    {
      if (! nTP++) fprintf(stderr, "cannot read targetPrefix\n") ;
      return FALSE ;
    }
  if (strcmp (ccp, "-"))
    dictAdd (ba->dnaDict,ccp, &(up3->targetPrefix)) ;
  
  ccp = aceInWordCut (ai, "\t", &cutter) ; /* col 21: targetSuffix */
  if (! ccp || ! *ccp)
    {
       if (! nTS++) fprintf(stderr, "cannot read targetSuffix\n") ;
      return FALSE ;
    }
  if (strcmp (ccp, "-") && (ba->keepTargetPrefix || ba->exportSuffix))
    dictAdd (ba->dnaDict,ccp, &(up3->targetSuffix)) ;
  
  /* if pair == false, do not compute dPair, but if known reexport it */
   deltaPair = 0 ;
   /* if no int, the file says -\t, we must step over the minus to see the t */ 
   aceInInt (ai, & deltaPair)  ; aceInStep (ai, '-') ; aceInStep (ai, '\t') ; /* col 22: pairLength */
   if (deltaPair == -14) deltaPair = 0 ;
   if (deltaPair)
    {
      ba->hasPair = TRUE ;
      
      up2->dPair = deltaPair ;
      if ( deltaPair < 0 && deltaPair >=  NON_COMPATIBLE_PAIR && deltaPair != -2  && deltaPair != -10 &&  deltaPair != -5)
	up2->badPair = 1 ;
    }
   
   up->gene = baSplitMrnaRemap (ba, up, up2) ;


   aceInWordCut (ai, "\t", &cutter) ; /* col 23 */
   aceInWordCut (ai, "\t", &cutter) ; /* col 24 */
   aceInWordCut (ai, "\t", &cutter) ; /* col 25 */
   ccp = aceInWord (ai) ; /* col 26 chain */
   if (ccp && ! strcasecmp (ccp, "chain"))
     {
       aceInInt (ai, &chain) ;
       up->chain = chain ;  
       aceInStep (ai, '\t') ; aceInInt (ai, &(up3->c1)) ;
       aceInStep (ai, '\t') ; aceInInt (ai, &(up3->c2)) ;
       aceInStep (ai, '\t') ; aceInInt (ai, &(up3->gt_ag)) ;
       aceInStep (ai, '\t') ; aceInInt (ai, &(up3->ct_ac)) ;
     }
     
  return TRUE ;
} /* baParseOneHit */

/*************************************************************************************/

static BOOL baParseOneSamTranscriptHit (ACEIN ai, BA *ba, HIT *up, int nn)
{ 
  AC_HANDLE h = ac_new_handle () ;
  char *ccp ;
  char buf[1000] ;
  int Z_genome = 0 ;
  HIT2 *up2 ;
  HIT3 *up3 ;
  int flag ;
  char *cigar ;
  Array cigarettes = arrayHandleCreate (128, SAMCIGAR, h) ; 

  memset (up, 0, sizeof (HIT)) ;
  if (ba->target_classDict) 
    dictFind (ba->target_classDict, "Z_genome", &Z_genome) ;
  strncpy (buf, aceInPos (ai), 999) ;
  ccp = aceInWord (ai) ;
  if (! ccp || ! *ccp || *ccp == '#' || *ccp == '@')
    goto done  ;
  dictAdd (ba->tagDict,ccp, &(up->tag)) ;
  up->nn = nn ;
  up2 = arrayp (ba->hits2, up->nn, HIT2) ;
  up3 = arrayp (ba->hits3, up->nn, HIT3) ;
  memset (up2, 0, sizeof (HIT2)) ;
  memset (up3, 0, sizeof (HIT3)) ;
  if (ba->pair)
    {
      char cc, *ccq = ccp + strlen (ccp) - 1  ;
      BOOL ok = TRUE ;

      if (ok) /* if no clone, do not enter it in the cloneDict */
	{
	  cc = *ccq ; *ccq = 0 ;
	  dictAdd (ba->cloneDict,ccp, &(up->clone)) ;
	  *ccq = cc ;
	}
    }
  aceInStep (ai, '\t') ;   aceInInt (ai, &flag) ;
 
  /*
    int isFirstFragment = 1 ;
    if (flag & 128) isFirstFragment = -1 ;
  */

  aceInStep (ai, '\t') ; ccp = aceInWord (ai) ;
  if (ccp && !strncmp (ccp,"MRNA:",5))ccp += 5 ;
  if (! ccp || ! *ccp)
    goto done  ;
  dictAdd (ba->targetDict,ccp, &(up->target)) ;
  
  up2->a1 = 0 ;
  aceInStep (ai, '\t') ; aceInInt (ai, &(up2->a1)) ;
  if (up2->a1)
    goto done  ;
 
  aceInStep (ai, '\t') ; aceInInt (ai, &(up2->a1)) ; /* quality, discard */
  aceInStep (ai, '\t') ;  cigar = aceInWord (ai) ;
  samParseCigar (cigar, cigarettes, up2->a1, &(up2->a2), &(up->x1), &(up->x2), 0, 0, 0, 0) ;

  up2->target_class = ba->target_class ;
  up2->dPair = 0 ;
  dictAdd (ba->target_classDict,ccp, &(up2->target_class)) ;
  aceInStep (ai, '\t') ;  aceInWord (ai) ; /* drop */
  aceInStep (ai, '\t') ;  aceInWord (ai) ; /* drop */
  aceInStep (ai, '\t') ;  aceInInt (ai, &(up2->dPair)) ;
  aceInStep (ai, '\t') ;  ccp = aceInWord (ai) ;  /* the actual sequence */
  if (!ccp)
    goto done ;
  up2->ln = strlen (ccp) ;

  up->score = 0 ; up2->mult = 1 ;
  while (ccp)
    {
      aceInStep (ai, '\t') ;  ccp = aceInWord (ai) ;  /* the actual sequence */
      if (!strncmp (ccp, "AS:", 3))
	sscanf (ccp+5, "%d", &(up->score)) ;
    }

  ac_free (h) ;  
  return TRUE ;

 done:
  ac_free (h) ;  
  return FALSE ;
} /* baParseOneSamTranscriptHit */

/*************************************************************************************/

static BOOL baFilterOne (BA *ba, HIT *up)
{
  int ngenes, score, ali ;
  int tZ = 0 ; /* dict entry for Z_genome */
  int tGS = 0 ; /* dict entry for s_snp */
  int tTS = 0 ; /* dict entry for t_snp */
  int tBACTERIA = 0 ; /* dict entry for b_bacteria */
  int tVIRUS = 0 ; /* dict entry for v_virus */
  int tTRANSPOSON = 0 ; /* dict entry for D_transposon */
  HIT2 *up2 ;
  HIT3 *up3 ;

  const char *t ;
  int deltaErrCost = 5,  deltaNCost = 2
    , partialAliPerCent = 20

    , absoluteMinAli = 16
    , minAli, minAliPerCent 
    , minTAli, minTAliPerCent
    , minTMAli, minTMAliPerCent 
    , minGAli, minGAliPerCent 
    , minGMAli, minGMAliPerCent 

    , minScore
    , minTScore, minGScore
    , minTMScore, minGMScore

    , partialScore
    , partialTScore, partialTMScore
    , partialGScore, partialGMScore
    ;


  if (ba->filter && (! strcmp (ba->filter, "worm") || ! strcmp (ba->filter, "ce") || ! strcmp (ba->filter, "at") || ! strcmp (ba->filter, "sc") || ! strcmp (ba->filter, "polyo") || ! strcmp (ba->filter, "fly") || ! strcmp (ba->filter, "Dmelanogaster")))
    {
      minTAli = 100 ;
      minTMAli = 130 ;
      minGAli = 120 ;
      minGMAli = 250 ;
      
      absoluteMinAli = 18 ;

      minTAliPerCent = 90 ;
      minTMAliPerCent  = 95 ;
      minGAliPerCent  = 90 ;
      minGMAliPerCent  = 95 ;

      minTScore = 18 ;
      minTMScore = 30 ;
      minGScore = 20 ;
      minGMScore = 32 ;

      partialTScore = 23 ;
      partialTMScore = 32 ;
      partialGScore = 27 ;
      partialGMScore = 50 ;
    }

  else if (ba->filter && *ba->filter >= '1' && *ba->filter <= '9')
    {
      int f ;
      const char *ccp ;

      for (f = 0, ccp = ba->filter ; *ccp ; ccp++)
	f = 10 * f + (*ccp - '0') ;

      minTAli = 100 ;
      minTMAli = 130 ;
      minGAli = 120 ;
      minGMAli = 250 ;
      
      minTAliPerCent = 50 ;
      minTMAliPerCent = 55 ;
      minGAliPerCent  = 50 ;
      minGMAliPerCent  = 55 ;

      absoluteMinAli = f ;

      minTScore = f ;
      minTMScore = f + 3 ;
      minGScore = f ;
      minGMScore = f + 3 ;

      partialTScore = f ;
      partialTMScore = f + 3 ;
      partialGScore = f ;
      partialGMScore = f + 3 ;
    }
  else
    {
      minTAli = 100 ;
      minTMAli = 130 ;
      minGAli = 120 ;
      minGMAli = 250 ;
      
      minTAliPerCent = 70 ;
      minTMAliPerCent  = 95 ;
      minGAliPerCent  = 70 ;
      minGMAliPerCent  = 95 ;
	  
      absoluteMinAli = 24 ;

      minTScore = 24 ;
      minTMScore = 32 ;
      minGScore = 26 ;
      minGMScore = 35 ;
      
      partialTScore = 27 ;
      partialTMScore = 50 ;
      partialGScore = 29 ;
      partialGMScore = 50 ;
    }

  /* 
     unicity < 0 means that a non stranded tag maps to genes on both strands
     if we want to remap on the genome via the transcripts -2 counts as non ambiguous
  */

  up2 = arrp (ba->hits2, up->nn, HIT2) ;
  up3 = arrp (ba->hits3, up->nn, HIT3) ;
  ngenes = up2->unicity ; 
  if (ba->pureNsStrand == 2 && ngenes == -2) ngenes = 1 ;
  if (ngenes < 0 && ba->pureNsStrand) return FALSE ;
  if (ngenes < 0) ngenes = - ngenes ; 

  /* the constitutive errCost is 8, alter to errCost (actually = 3) */

  if (ba->errCost) deltaErrCost = ba->errCost - 3 ;
  if (ba->solid)  deltaNCost = 2 ;

  score = up->score + deltaErrCost * up3->nErr + deltaNCost * up3->nN ;
  ali = up2->ali ;

  if (up2->dPair > 0 || (up2->dPair <  NON_COMPATIBLE_PAIR && ! up2->badPair))
    {
      score *= 1.5 ;
      ali *= 1.5 ;
    }
  up2 = arrp (ba->hits2, up->nn, HIT2) ;
  up3 = arrp (ba->hits3, up->nn, HIT3) ;
  /* alter the aligned length by adding in the AAA */
  {
    const char *ccp ;
    int n = 0 ;
    
    ccp =  up3->suffix ? dictName (ba->dnaDict, up3->suffix) : 0 ;
    if (ccp)
      {
	while (ccp && *ccp++ == 'A') n++ ;
	if (*ccp) n = 0 ;
	ali += n ;
      } 
    if (n == 0)
      {
	ccp =  up3->suffix ? dictName (ba->dnaDict, up3->suffix) : 0 ;
	if (ccp)
	  {
	    while (*ccp++ == 'A') n++ ;
	    if (*ccp) n = 0 ;
	    ali += n ; 
	  }
      }
  }

  dictFind (ba->target_classDict, "b_bacteria", &tBACTERIA) ;
  dictFind (ba->target_classDict, "v_virus", &tVIRUS) ;
  dictFind (ba->target_classDict, "D_transposon", &tTRANSPOSON) ;
  if (
      (ba->maxHit && ngenes >= ba->maxHit) && /* in this case, we may have an incomplete list */
      up2->target_class != tVIRUS &&  up2->target_class != tBACTERIA && up2->target_class != tTRANSPOSON
      )
    {
      if (! bitt (ba->maxHitBitSet, up->tag))
	{
	  bitSet (ba->maxHitBitSet, up->tag) ;
	  aceOutf(ba->maxHitOut, "%s\tMulti\t%d\t%d\n", dictName (ba->tagDict,up->tag), score, up2->mult); 
	  badMaxHit += up2->mult ;
	}
      return FALSE ;
    }  
  dictFind (ba->target_classDict, "s_snp", &tGS) ;
  dictFind (ba->target_classDict, "t_snp", &tTS) ;
  dictFind (ba->target_classDict, "Z_genome", &tZ) ;
   minAli = minGAli ; minAliPerCent = minGAliPerCent ; partialScore =  partialGScore ; minScore = minGScore ;
  if (ngenes == 1 && up2->target_class != tZ)
    {
      t = "T" ; 
      minAli = minTAli ; minAliPerCent = minTAliPerCent ; minScore = minTScore ;
      partialScore = partialTScore ;
    } 
  if ( up2->target_class == tVIRUS ||  up2->target_class == tBACTERIA ||  up2->target_class == tTRANSPOSON)
    {
      if (up2->unicity < 0) up2->unicity *= -1 ;
      if (up2->unicity > 9) up2->unicity = 9 ;
    }
  if (up2->unicity != -2 && ngenes > 1 && up2->target_class != tZ)
    {
      t = "TM" ; 
      minAli = minTMAli ; minAliPerCent = minTMAliPerCent ; minScore = minTMScore ;
      partialScore = partialTMScore ;
    }
  if (ngenes == 1 && up2->target_class == tZ)
    {
      t = "G" ; 
      minAli = minGAli ; minAliPerCent = minGAliPerCent ; minScore = minGScore ;
      if (up2->ln >= minTScore && up2->ln <= minGScore)
	minScore = minTScore ;
      partialScore = partialGScore ;
    }    
  if (ngenes > 1 && up2->target_class == tZ)
    {
      t = "GM" ; 
      minAli = minGMAli ; minAliPerCent = minGMAliPerCent ; minScore = minGMScore ;
      partialScore = partialGMScore ;
    }

  if (ali < absoluteMinAli || (ali <= minAli && 100 * ali < partialAliPerCent * up2->ln) )
    {
      badLnS++ ; badLnT += up2->mult ; 
      if (! bitt (ba->maxHitBitSet, up->tag))
	{
	  bitSet (ba->maxHitBitSet, up->tag) ;
	  aceOutf(ba->maxHitOut, "%s\tPartial\t%d\n", dictName (ba->tagDict,up->tag), score); 
	  badMaxHit += up2->mult ;
	}
      return FALSE ;
    }
 
   if (ali >= minAli || 100 * ali >= minAliPerCent * up2->ln)      
    {
      if (score >= minScore)
	{
	  isLongS++ ; isLongT += up2->mult ;
	}
      else
	{
	  badScoreLS++ ; badScoreLT += up2->mult ;
	  if (! bitt (ba->maxHitBitSet, up->tag))
	    {
	      bitSet (ba->maxHitBitSet, up->tag) ;
	      aceOutf(ba->maxHitOut, "%s\tBadScore\t%d\n", dictName (ba->tagDict,up->tag), score) ;
	      badMaxHit += up2->mult ;
	    }
	  return FALSE ; 
	}
    }
  else
    {
      if (score >= partialScore)
	{
	  isPartialS++ ; isPartialT += up2->mult ;
	}
      else
	{
	  badScorePS++ ; badScorePT += up2->mult ;
	  if (! bitt (ba->maxHitBitSet, up->tag))
	    {
	      bitSet (ba->maxHitBitSet, up->tag) ;
	      aceOutf(ba->maxHitOut, "%s\tPartial\t%d\n", dictName (ba->tagDict,up->tag), score); 
	      badMaxHit += up2->mult ;
	    }
	  return FALSE ; 
	}
    }
  if (ngenes > 1)
    {
      if (! bitt (ba->maxHitBitSet, up->tag))
	{
	  bitSet (ba->maxHitBitSet, up->tag) ;
	  aceOutf(ba->maxHitOut, "%s\tOK_Multi\t%d\t%d\n", dictName (ba->tagDict,up->tag), score,ngenes) ; 
	  badMaxHit += up2->mult ;
	}
    }
      
  if (0) fprintf(stderr, "t=%s\n",t); /* for compiler happiness */
  return TRUE ;
} /* baFilterOne */

/*************************************************************************************/
/* open the nth file on the inFileList */
static ACEIN baNextInFile (BA *ba, int nn, AC_HANDLE h)
{
  ACEIN ai = 0 ;

  if (! ba->inFileDict)
    {
      AC_HANDLE h1 = ac_new_handle () ;
      
      ba->inFileDict = dictHandleCreate (128, ba->h) ;
      ai = aceInCreate (ba->inFileList, FALSE, h1) ;
      while (aceInCard (ai))
	{
	  const char *ccp = aceInWord (ai) ;
	  if (ccp)
	    dictAdd (ba->inFileDict, ccp, 0) ;
	}
      ac_free (h1) ;
    }
  if ( ba->inFileDict && nn < dictMax (ba->inFileDict))
    {
      ai = aceInCreate (dictName (ba->inFileDict, nn+1), FALSE, h) ;
    }
  return ai ;
} /* baNextInFile */

/*************************************************************************************/
/*************************************************************************************/

static BOOL baKeepOneHit (BA *ba, HIT *up, const int nn, const int Z_genome)
{ 
  Array aa2, aa3 ;
  HIT *vp ;
  HIT2 *up2, *vp2 ;
  HIT3 *up3, *vp3 ; 
  int iv ;
  int tag = up->tag ;
  int score = up->score ;

  aa2 = ba->hits2 ;
  aa3 = ba->hits3 ;

  up2 = arrayp (aa2, up->nn, HIT2) ;
  up3 = arrayp (aa3, up->nn, HIT3) ;
  if (ba->minAli && up2->ali < ba->minAli) /* reject */
    return FALSE ;
  
  /* compare to previous hits and discard if overlap and lower score */
  for (iv = nn - 1, vp = up - 1 ; iv >= 0 && vp->tag == tag ; iv--, vp--)
    {
      if (vp->score != score) /* check overlap */
	{
	  vp2 = arrayp (aa2, vp->nn, HIT2) ;
	  vp3 = arrayp (aa3, vp->nn, HIT3) ;
	  
	  /* do not kill big transcript in favor of genome */
	  if (score > 200 &&
	      Z_genome &&
	      vp2->target_class == Z_genome &&
	      up2->target_class < Z_genome
	      )
	    continue ;
	  if (score > vp->score && vp->score > 200 &&
	      Z_genome &&
	      up2->target_class == Z_genome &&
	      vp2->target_class < Z_genome
	      )
	    continue ;
	  
	  /* check if real intersect, eliminate loser */
	  {
	    int z1 = up3->c1 > vp3->c1 ?  up3->c1 : vp3->c1 ;
	    int z2 = up3->c2 < vp3->c2 ?  up3->c2 : vp3->c2 ;
	    int dz = z2 - z1 + 1 ;  /* length of intersect */
	    int lnu = up3->c2 - up3->c1 + 1 ;
	    int lnv = vp3->c2 - vp3->c1 + 1 ;
	    
	    if (dz > .5 * lnu || dz > .5 * lnv)  /* real intersect */
	      {
		if (score < vp->score)
		  {
		    score = up->score = 0 ;
		    break ;
		  }
		else if (score > vp->score)
		  vp->score = 0 ;
	      }
	    if (
		(vp3->hardClip1 && up->x2 < vp3->hardClip1 + 30) ||
		(vp3->hardClip2 && up->x1 < vp3->hardClip2 - 30)
		)
	      {
		score = up->score = 0 ;
		break ;
	      }
	  }
	}
    }
  return score ? TRUE : FALSE ;
}

static int baParseHits (BA *ba) 
{
  Array aa ;
  HIT *up ;
  int nn = 0 ;
  int nInFile = 0 ;
  int Z_genome = 0 ;
  ACEIN ai = 0 ;

  if (! ba->maxHitOut)
    {
      ba->maxHitOut = aceOutCreate (ba->outFileName, ".too_many_hits", ba->gzo, ba->h) ;
      ba->maxHitBitSet = bitSetCreate (10000, ba->h) ;
    }
  if (ba->target_classDict) 
    dictFind (ba->target_classDict, "Z_genome", &Z_genome) ;

  ba->hits = aa = arrayHandleCreate (1000000, HIT, ba->h) ;
  ba->hits2 = arrayHandleCreate (1000000, HIT2, ba->h) ;
  ba->hits3 = arrayHandleCreate (1000000, HIT3, ba->h) ;
  ba->tagDict = dictHandleCreate (10000, ba->h) ;
  if (ba->pair) ba->cloneDict = dictHandleCreate (10000, ba->h) ;
  if (! ba->targetDict)
    ba->targetDict = dictHandleCreate (10000, ba->h) ;
  ba->errDict = dictHandleCreate (10000, ba->h) ;
  ba->dnaDict = dictCaseSensitiveHandleCreate (10000, ba->h) ; /* distinguish aaa from AAA */

  up = arrayp (aa, 0, HIT) ;

  while (1) /* read all input files */
    {
      AC_HANDLE h = ac_new_handle () ;

      if (! ba->inFileList)
	ai = aceInCreate (ba->hitFileName, ba->gzi, h) ;
      else
	ai = baNextInFile (ba, nInFile++, h) ;
      if (ai)
	while (aceInCard (ai))
	  if (baParseOneHit (ai, ba, up, nn) &&
	      baKeepOneHit (ba,up,nn, Z_genome)
	      )     /* register new chain */
	    up = arrayp (aa, ++nn, HIT) ;
    
      ac_free (h) ; 
      if (!ai || ! ba->inFileList)
	break ;
    }
  arrayMax (aa) = nn ;
  fprintf (stderr, "Parsed %d lines %s\n", nn, timeShowNow ()) ;

  baTargetDict = ba->targetDict ;
  arraySort (aa, baHitOrder) ;
  baHitCompress (ba, aa, Z_genome) ;
  nn = arrayMax (aa) ;

  fprintf (stderr, "Sorted %d distinct lines %s\n", nn, timeShowNow ()) ;

  return nn ;
} /* baParseHits */

/*************************************************************************************/
/*************************************************************************************/
  /* filter on quality */
static int baFilter (BA *ba)
{
  int ii, jj, nn = 0 ;
  HIT *up, *vp ;
  HIT2 *up2 ; 
  Array aa = ba->hits ; 
  Array aa2 = ba->hits2 ; 
  int tBACTERIA = 0 ;
  int tVIRUS = 0 ;
  int tTRANSPOSON = 0 ;

  dictFind (ba->target_classDict, "b_bacteria", &tBACTERIA) ;
  dictFind (ba->target_classDict, "v_virus", &tVIRUS) ;
  dictFind (ba->target_classDict, "D_transposon", &tTRANSPOSON) ;

  if (arrayMax (aa))
    {
      for (ii = jj = 0, up = arrp (aa, 0, HIT), vp = up ; ii < arrayMax (aa) ; up++, ii++)
	{ 
	  if (!up->score)
	    continue ;
	  if (baFilterOne (ba, up))
	    {
	      if (jj < ii) *vp = *up ;
	      vp++ ; jj++ ;
	    }
	  else 
	    nn++ ;
	}
      arrayMax (aa) = jj ;
      
      for (ii = jj = 0, up = arrp (aa, 0, HIT), vp = up ; ii < arrayMax (aa) ; up++, ii++)
	{ 
	  up2 = arrayp (aa2, up->nn, HIT2) ;
	  if (ba->maxHit && 
	      (up2->unicity >= ba->maxHit || -up2->unicity >= ba->maxHit) &&
	       up2->target_class != tVIRUS &&  up2->target_class != tBACTERIA && up2->target_class != tTRANSPOSON
	      )
	    {
	      if (! bitt (ba->maxHitBitSet, up->tag))
		{
		  bitSet (ba->maxHitBitSet, up->tag) ;
		  aceOutf(ba->maxHitOut, "%s\tMulti\t%d\t%d\n", dictName (ba->tagDict,up->tag), up->score, up2->mult); 
		  badMaxHit += up2->mult ;
		}
	    }
	}
    }
  fprintf (stderr, "Accepted %d lines %s\n", jj, timeShowNow ()) ;
  
  return nn ;
} /* baFilter */

/*************************************************************************************/
/*************************************************************************************/

static int baPairFilterOne (BA *ba, int iiMin, int iiMax, Array geneLinks)
{
  AC_HANDLE h = ac_new_handle () ;
  Array aa = ba->hits ;
  Array aa2 = ba->hits2 ;
  Array bb = arrayHandleCreate (32, HIT, h) ;
  Array bb2 = arrayHandleCreate (32, HIT2, h) ;
  HIT *up, *vp, *wp ;
  HIT2 *up2, *vp2, *wp2, *xp2 ;
  int tBACTERIA = 0 ; /* dict entry for b_bacteria */
  int tVIRUS = 0 ; /* dict entry for v_virus */
  int tTRANSPOSON = 0 ; /* dict entry for D_transposon */
  
  int cl, i1, i2, ii, jj, bestjj, nn = 0,gene1,  tag1, target1, nHits, nHits1 ;
  int a1, a2, b1, b2, x1, bestDa, score, tZ, tMito, tRrna, tChloro, tSpikein, tDNASpikein, tVirus, tBacteria ;
  KEYSET genes = keySetHandleCreate (h) ;
  int ng, geneStrand ;
  BOOL hasNegativeUnicity = FALSE, inGene ;
  int hasGoodClone = 0 ;
  int hasBestClone = 0 ;
  int maxDa = 1000 ;

  dictFind (ba->target_classDict, "b_bacteria", &tBACTERIA) ;
  dictFind (ba->target_classDict, "v_virus", &tVIRUS) ;
  dictFind (ba->target_classDict, "D_transposon", &tTRANSPOSON) ;
  
  if (ba->pair > 1) maxDa = 3 * ba->pair ;
  if (ba->RNA_seq && maxDa < 10000) maxDa = 10000 ; /* allow for introns */

  dictFind (ba->target_classDict, "Z_genome", &tZ) ;
  dictFind (ba->target_classDict, "A_mito", &tMito) ;
  dictFind (ba->target_classDict, "B_rrna", &tRrna) ;
  dictFind (ba->target_classDict, "C_chloro", &tChloro) ;
  dictFind (ba->target_classDict, "0_SpikeIn", &tSpikein) ;
  dictFind (ba->target_classDict, "1_DNASpikeIn", &tDNASpikein) ;
  dictFind (ba->target_classDict, "v_virus", &tVirus) ;
  dictFind (ba->target_classDict, "b_bacteria", &tBacteria) ;

  for (cl = hasBestClone = hasGoodClone = 0 ; cl <= dictMax (ba->target_classDict) ; cl++)
    {
      hasNegativeUnicity = FALSE ; hasGoodClone = 0 ; 
      for (ii = iiMin, up = arrp (aa, ii, HIT), tag1 = up->tag, nHits1 = nHits = 0 ; ii < iiMax ; up++, ii++)
	{ 
	  HIT2 *bestUp2 = 0, *bestVp2 = 0 ;

	  up2 = arrp (aa2, up->nn, HIT2) ;
	  if (up2->target_class != cl)
	    continue ;
	  if (up2->unicity < 0) 
	    hasNegativeUnicity = TRUE ;
	  nHits++ ;
	  if (up->tag != tag1) 
	    continue ;
	  nHits1++ ;
	  inGene = FALSE ;
	  bestDa = (cl == tZ && ba->RNA_seq ? 1000000 : maxDa) ; /* at most 1 kb except for introns */
	  /* ATTENTION 2014_10_12, i should debug this it false if RNA-seq == 0 which is stupid, i then get zero type -4 nToofarOnGenonme */
	  a1 = up2->a1 ; a2 = up2->a2 ; x1 = up->x1 ;
	  /* search best partner of ii */
	  for (jj = ii + 1, vp = up + 1, bestjj = 0, gene1 = up->gene, target1 = up->target ; jj < iiMax ; vp++, jj++)
	    {
	      vp2 = arrp (aa2, vp->nn, HIT2) ;
	      if (vp2->target_class != cl || vp->tag == tag1)
		continue ;

	      if (vp->target == target1)
		{
		  b1 = vp2->a1 ; b2 = vp2->a2 ;
		  if (a1 < a2 && b1 > b2 && b1 + 10 > a1 && a2 - 10 <= b1 && a1 - 10 <= b2 && b1 - a1 + 1 + x1 + vp->x1 - 2 < bestDa)
		    { bestDa = b1 - a1  + 1 + x1 + vp->x1 - 2 ; bestjj = jj ; bestUp2 = up2 ; bestVp2 = vp2 ;}
		  if (a1 > a2 && b1 < b2 && b1 - 10 < a1 && b2 - 10 <= a1 && b1 - 10 <= a2 && a1 - b1  + 1 + x1 + vp->x1 - 2 < bestDa)
		    { bestDa = a1 - b1  + 1 + x1 + vp->x1 - 2 ; bestjj = jj ; bestUp2 = vp2 ;  bestVp2 = up2 ;}
		}
	      else
		{
		  if (
		      gene1 && vp->gene == gene1 && ! bestjj)
		    { inGene = TRUE ; bestjj = jj ; }
		  continue ;
		} 
	    } 
	  if (bestjj) 
	    {
	      wp = arrayp (bb, nn++, HIT) ;
	      if (bestDa < 0 && bestDa >  NON_COMPATIBLE_PAIR) bestDa =  NON_COMPATIBLE_PAIR  -1 ; /* NON_COMPATIBLE_PAIR - 1 is first compatible value, synchronize to hack, values up to -14 are types of pairs */
	      wp->x1 = ii ; wp->score = 1 ; wp->x2 = bestjj ; 
	      wp->nn = (bestDa < maxDa ? bestDa : 1000000) ;

	      if (bestUp2) bestUp2->dPair = wp->nn ;
	      if (bestVp2) bestVp2->dPair = -wp->nn ;
	      if (bestDa < maxDa && hasGoodClone <= 0) 
		{
		  if (cl == tMito && (hasGoodClone <= 0 || hasGoodClone > 4)) hasGoodClone = 4 ;
		  else if (cl == tChloro && (hasGoodClone <= 0 || hasGoodClone > 4)) hasGoodClone = 4 ;
		  else if (cl == tSpikein && (hasGoodClone <= 0 || hasGoodClone > 3)) hasGoodClone = 3 ;
		  else if (cl == tZ && (hasGoodClone <= 0 || hasGoodClone > 2))  hasGoodClone = 2 ;
		  else if (up->gene && (hasGoodClone <= 0 || hasGoodClone > 1)) hasGoodClone = 1 ;
		  else if (hasGoodClone <= 0 || hasGoodClone > 5) hasGoodClone = 5 ;
		}
	      else if (inGene && !hasGoodClone)
		hasGoodClone = -5 ;
	    }
	}
      /* do we have several ii going to the same jj */
      if (arrayMax (bb) && arrayMax (bb2))
	{
	  HIT *xp ;

	  for (i1 = 0, wp = arrp (bb, 0, HIT), wp2 = arrp (bb2, 0, HIT2) ; i1 < arrayMax (bb) && i1 < arrayMax (bb2) ; i1++, wp++, wp2++)
	    {
	      if (! wp->score || ! wp->x2) continue ;
	      for (i2 = i1 + 1, xp = wp + 1, xp2 = wp2 + 1 ; i2 < arrayMax (bb) ; i2++, xp++, xp2++)
		{
		  if (! xp->score) continue ;
		  if (! xp->gene && ! wp->gene && xp->x2 == wp->x2)
		    { 
		      if (wp->nn < xp->nn) xp->score = 0 ;
		      if (xp->nn < wp->nn) wp->score = 0 ;
		    }
		}
	    }
	 
	  if (1) /* measure pair distances in good pairs */
	    {
	      for (ii = iiMin, up = arrp (aa, ii, HIT) ; ii < iiMax ; up++, ii++)
		{ 
		  up2 = arrp (aa2, up->nn, HIT2) ;
		  if (up2->target_class != cl)
		    continue ;
		  for (i1 = score = 0, wp = arrp (bb, 0, HIT) ; i1 < arrayMax (bb) ; i1++, wp++)
		    {
		      if (! wp->score) continue ;
		      if (wp->x1 == ii)
			{ score = up->score ; up2->dPair = (wp->nn == 1000000 ? -5 : wp->nn) ; break ; }
		      if (wp->x2 == ii)
			{ score = up->score ; up2->dPair = (wp->nn == 1000000 ? -5 : - wp->nn)  ; break ; }
		    }
		  if (0) /* 2019_10_10 : this kills the gene fusions in case of single read and anyway we now better trust the clipalign clean up */
		    if (! ba->geneSupport)  /* 2015_05, delay destruction till after counting all anomalies, a good idea in geneSupport, but not destroying may be bad in case of exportBest */
		      up->score = score ? score : - up->score ;  /* if I have at least one pair destroy all other hits */
		}
	      /* restore the negative scores if we are in a chain, or kill it */
	      for (ii = iiMin, up = arrp (aa, ii, HIT) ; ii < iiMax ; up++, ii++)
		{ 
		  if (up->score >= 0)
		    continue ;
		  if (! up->chain)
		    continue ;
		  for (jj = iiMin,vp = arrp (aa, jj, HIT) ; jj < iiMax ; vp++, jj++)
		    {
		      if (vp->score > 0 && up->tag == vp->tag && up->target == vp->target)
			{
			  if (vp->chain == up->chain)
			    { 
			      up->score = vp->score ;
			      up2 = arrp (aa2, up->nn, HIT2) ;
			      vp2 = arrp (aa2, vp->nn, HIT2) ;
			      up2->dPair = vp2->dPair ;
			      break ; 
			    }
			}
		    }
		  if (up->score < 0)
		    up->score = 0 ;
		}
	    }
	  /* reevaluate the multiplicities */
	  for (ii = iiMin, up = arrp (aa, ii, HIT), tag1 = 0 ; ii < iiMax ; up++, ii++)
	    { 
	      up2 = arrp (aa2, up->nn, HIT2) ;
	      keySetMax (genes) = 0 ;
	      geneStrand = ng = 0 ;
	      if (up2->target_class != cl)
		continue ;
	      if (up->tag == tag1)
		continue ;
	      tag1 = up->tag ;
	      for (jj = ii, vp = up, nHits1 = nHits = 0 ; jj < iiMax && vp->tag == tag1 ; vp++, jj++)
		{
		  vp2 = arrp (aa2, vp->nn, HIT2) ;
		  if (vp->score && vp2->target_class == cl)
		    { 
		      if (! vp->gene)
			nHits++ ;
		      else
			{
			  keySet(genes, ng++) = vp->gene ;
			  if (vp2->a1 < vp2->a2) geneStrand |= 1 ; else  geneStrand |= 2 ;
			}
		    }
		}
	      if (geneStrand)
		{
		  keySetSort (genes) ; 
		  keySetCompress (genes) ;
		  nHits = keySetMax (genes) ;
		  if (geneStrand == 0x3) nHits = -nHits ;
		}
	      else if (hasNegativeUnicity && nHits1 > 0 && nHits1 < nHits) /* hits on both strands */
		nHits = -nHits ;
	      for (jj = ii, vp = up ; jj < iiMax && vp->tag == tag1 ; vp++, jj++)
		{
		  vp2 = arrp (aa2, vp->nn, HIT2) ;
		  if (vp->score && vp2->target_class == cl)
		    {
		      /* because of dicontinuous ali, evaluated by the aligner itself, we can reduce but not augment multiplicity by filstering */
		      int old =  vp2->unicity, i, j ;
		      i = old ; if (i < 0) i = -i ;
		      j =  nHits ; if (j < 0) j = -j ;
		      if (i > j)
			vp2->unicity = nHits ;
		    }
		}
	    }
	}

      if (hasGoodClone > 0 &&  (hasBestClone <= 0 || hasBestClone > hasGoodClone))
	hasBestClone = hasGoodClone ; 
      else if (hasBestClone < 0 && hasGoodClone < 0 && hasBestClone < hasGoodClone)
	hasBestClone = hasGoodClone ; 
      else if (! hasBestClone  && hasGoodClone< 0)
	hasBestClone = hasGoodClone ; 
    }
  hasGoodClone = hasBestClone ;

  /* up to here a good pair is compact, with correct topology, in a single transcript or genome region */  
  if (hasGoodClone)  /* destroy all hits in all targets that do not come as a good pair
		      * but do not destroy transcriptome hits in favor of genome hits 
		      */
    if (0)
      {
	for (ii = iiMin, up = arrp (aa, ii, HIT) ; ii < iiMax ; up++, ii++)
	  { 
	    up2 = arrp (aa2, up->nn, HIT2) ;
	    if (! up2->dPair && (hasGoodClone != 2 || !up->gene))
	      up->score = 0 ;
	  }
      }
  /* label pairs linking two distant genes (not seen as case 6: genome compatible) */
  if (! hasGoodClone || hasGoodClone == 2)
    {
      /* look for gene pairs */
      for (ii = iiMin, up = arrp (aa, ii, HIT) ;ii < iiMax ; up++, ii++)
	{
	  if (! up->gene) continue ;
	  tag1 = up->tag ;
	  up2 = arrp (aa2, up->nn, HIT2) ;
	  for (jj = ii, vp = up, nHits1 = nHits = 0 ; jj < iiMax ; vp++, jj++)
	    { 
	      BOOL ok= FALSE ;
	      BOOL closeCisGenes = FALSE ;
	      if (vp->tag == tag1)
		continue ;

	      vp2 = arrp (aa2, vp->nn, HIT2) ;
	      if (vp->tag != tag1 && vp->gene && vp->gene != up->gene && vp2->target_class == up2->target_class)
		{
		  if (1)
		    {
		      HIT *zp = arrayp (geneLinks, arrayMax(geneLinks), HIT) ;
		      zp->score = up2->mult ;
		      if (lexstrcmp(dictName(ba->targetDict,up->gene),  dictName(ba->targetDict,vp->gene)) < 0)
			{ zp->clone = up->clone ; zp->tag = up->gene ; zp->gene = vp->gene ; }
		      else
			{ zp->clone = up->clone ; zp->tag = vp->gene ; zp->gene = up->gene ; }
		    }
		  if (hasGoodClone == 2)
		    {
		      ok = TRUE ;
		      up2->dPair = KEYMAKE(2,vp->gene) ;
		      vp2->dPair = KEYMAKE(2,up->gene) ;
		    }
		  else
		    {
		      if (ba->gene2map &&
			  up->gene < arrayMax (ba->gene2map) &&
			  vp->gene < arrayMax (ba->gene2map) 
			  )
			{
			  HIT *g1 = arrp (ba->gene2map, up->gene, HIT) ;
			  HIT *g2 = arrp (ba->gene2map, vp->gene, HIT) ;
			  
			  if (g1->target == g2->target && 
			      (  /* genes must be in cis and the reads must face each other: 4 orientations */
			       (
				g1->x1 < g1->x2 && g2->x1 < g2->x2 &&
				g1->x1 < g2->x1 &&
				up2->a1 < up2->a2 && vp2->a1 > vp2->a2
				) ||
			       (
				g1->x1 < g1->x2 && g2->x1 < g2->x2 &&
				g1->x1 > g2->x1 &&
				up2->a1 > up2->a2 && vp2->a1 < vp2->a2
				) ||
			       (
				g1->x1 > g1->x2 && g2->x1 > g2->x2 &&
				g1->x1 > g2->x1 &&
				up2->a1 < up2->a2 && vp2->a1 > vp2->a2
				) ||
			       (
				g1->x1 > g1->x2 && g2->x1 > g2->x2 &&
				g1->x1 < g2->x1 &&
				up2->a1 > up2->a2 && vp2->a1 < vp2->a2
				) 
			       )
			      )
			    closeCisGenes = TRUE ;
			}
		      if (closeCisGenes)
			{
			  ok = TRUE ;
			  up2->dPair = KEYMAKE(10,vp->gene) ;
			  vp2->dPair = KEYMAKE(10,up->gene) ;
			}
		      else
			{
			  ok = TRUE ;
			  up2->dPair = KEYMAKE(6,vp->gene) ;
			  vp2->dPair = KEYMAKE(6,up->gene) ;
			}
		    }
		}
	      if (ok)
		{
		  if (hasGoodClone == 2) /* gene to close-by genome */
		    hasGoodClone = -2 ; 
		  else if (closeCisGenes)
		    hasGoodClone = -10 ; 
		  else
		    hasGoodClone = -6 ; 
		}
	    }
	}
    }
  if (hasGoodClone == 2)  /* good genome genome, is it anchored in a gene ? */
    {
      for (ii = iiMin, up = arrp (aa, ii, HIT) ; ii < iiMax ; up++, ii++)  
	{
	  up2 = arrp (aa2, up->nn, HIT2) ;
	  if (up->score && up->gene && ! up2->dPair)
	    { hasGoodClone = 6 ; up2->dPair = -2 ; }
	}
    }
  /* label incompatible pairs */
  if (! hasGoodClone)
    {
      /* look for distant pairs */
      for (ii = iiMin, up = arrp (aa, ii, HIT) ;ii < iiMax ; up++, ii++)
	{
	  up2 = arrp (aa2, up->nn, HIT2) ;
	  if (up->score &&  up2->dPair > 0)
	    { 
	      if (up->gene) hasGoodClone = -3 ;
	      else if (! hasGoodClone) hasGoodClone = -4 ;
	    }
	}
    }

  if (! hasGoodClone)
    {
      /* look for gene rRNA */
      for (ii = iiMin, up = arrp (aa, ii, HIT) ; ! hasGoodClone && ii < iiMax ; up++, ii++)
	{
	  up2 = arrp (aa2, up->nn, HIT2) ;
	  if (! (up2->target_class == tMito || up2->target_class == tRrna || up2->target_class == tSpikein ||  up->gene))
	    continue ;
	  tag1 = up->tag ; 
	  for (jj = iiMin, vp = arrp (aa, jj, HIT), nHits1 = nHits = 0 ; jj < iiMax ; vp++, jj++)  /* jj = ii, vp = up, nHits1 = nHits = 0 ; jj < iiMax ; vp++, jj++ */
	    { 
	      if (vp->tag == tag1)
		continue ;
	      vp2 = arrp (aa2, vp->nn, HIT2) ;
	      if ((up2->target_class == tSpikein || vp2->target_class == tSpikein) &&   up2->target_class !=  vp2->target_class)
		{ up2->dPair = vp2->dPair = hasGoodClone = -11 ; } 
	      if ((up2->target_class == tRrna || vp2->target_class == tRrna) &&   up2->target_class !=  vp2->target_class)
		{ up2->dPair = vp2->dPair = hasGoodClone = -7 ; } 
	      else if ((up2->target_class == tMito  || vp2->target_class == tMito) &&   up2->target_class !=  vp2->target_class)
		{ up2->dPair = vp2->dPair = hasGoodClone = -8 ;  } 
	    }
	}
    }
  
  /* label possible gene extensions from gene to genome */
  if (! hasGoodClone)
    {
      /* look for gene pairs */
      for (ii = iiMin, up = arrp (aa, ii, HIT) ;ii < iiMax ; up++, ii++)
	{
	  if (! up->gene) continue ;
	  tag1 = up->tag ; 
	  up2 = arrp (aa2, up->nn, HIT2) ;
	  for (jj = ii, vp = up, nHits1 = nHits = 0 ; jj < iiMax ; vp++, jj++)
	    { 
	      if (vp->tag == tag1)
		continue ;
	      vp2 = arrp (aa2, vp->nn, HIT2) ;
	      if (ii != jj && vp->tag == tag1 && vp2->target_class == tZ)
		break ;
	    }
	  if (jj < iiMax) continue ; /* up is not on an intron since it maps the genome */
	  tag1 = up->tag ;
	  for (jj = iiMin, vp = up, nHits1 = nHits = 0 ; jj < iiMax ; vp++, jj++)
	    { 
	      if (vp->tag == tag1)
		continue ;
	      vp2 = arrp (aa2, vp->nn, HIT2) ;
	      if (vp->tag != tag1 && vp2->target_class == tZ)
		{ up2->dPair = vp2->dPair = hasGoodClone = -9 ; }
	    }
	}
    }

  /* label all other as incompatible */
  if (! hasGoodClone)
    hasGoodClone = -14 ;
  for (ii = iiMin, up = arrp (aa, ii, HIT) ;ii < iiMax ; up++, ii++)
    {
      up2 = arrp (aa2, up->nn, HIT2) ;
      if (up->score &&  ! up2->dPair)
	up2->dPair = -14 ;
    }

  if (ba->maxHit)
    for (ii = iiMin, up = arrp (aa, ii, HIT) ; ii < iiMax ; up++, ii++)
      { 
	up2 = arrp (aa2, up->nn, HIT2) ;
	if (
	    (up2->unicity >= ba->maxHit || -up2->unicity >= ba->maxHit) &&
	    up2->target_class != tVIRUS &&  up2->target_class != tBACTERIA && up2->target_class != tTRANSPOSON
	    )
	  {
	    if (! bitt (ba->maxHitBitSet, up->tag))
	      {
		bitSet (ba->maxHitBitSet, up->tag) ;
		aceOutf(ba->maxHitOut, "%s\tMulti\t%d\t%d\n", dictName (ba->tagDict,up->tag), up->score, up2->mult); 
		badMaxHit += up2->mult ;
	      }
	    up->score = 0 ;
	  }
      }
  if (hasGoodClone < -4) /* check if we have 2 hits to the same gene */
    {
      inGene = FALSE ;
      for (ii = iiMin, up = arrp (aa, ii, HIT), tag1 = up->tag ; ii < iiMax ; up++, ii++)
	{
	  if (! up->gene || ! up->score || tag1 != up->tag) continue ;
	  up2 = arrp (aa2, up->nn, HIT2) ;
	  a1 = up2->a1 ; a2 = up2->a2 ;
	  /* search best partner of ii */
	  for (jj = ii + 1, vp = arrp (aa, jj, HIT), bestjj = 0, gene1 = up->gene ; jj < iiMax ; vp++, jj++)
	    {
	      vp2 = arrp (aa2, vp->nn, HIT2) ;
	      if (ii == jj || !vp->score || vp2->target_class != up2->target_class || vp->tag == tag1)
		continue ;
	      if (vp->gene == gene1)
		{ 
		  b1 = vp2->a1 ; b2 = vp2->a2 ;
		  if (up->target == vp->target && 
		      ((a1 < a2 && b1 < b2) || (a1 > a2 && b1 > b2))
		      )
		    { up2->dPair = vp2->dPair = -14 ; } /* incompatible topology */
		  else if (
			   up->target != vp->target || 
			   (a1 < a2 && b1 > b2 && b1 + 10 > a1 && a2 - 10 <= b1 && a1 - 10 <= b2) ||
			   (a1 > a2 && b1 < b2 && b1 - 10 < a1 && b2 - 10 <= a1 && b1 - 10 <= a2)
			   )
		    { inGene = TRUE ; up2->dPair = vp2->dPair = -5 ; }
		}
	    }
	}
      if (inGene) 
	hasGoodClone = -5 ;
    }
  if (ba->Remove_inserts_shorter_than)
     {
       BOOL ok = FALSE ;
       int deltaPair ;
       int min = ba->Remove_inserts_shorter_than ;
       for (ii = iiMin, up = arrp (aa, ii, HIT) ; ii < iiMax ; up++, ii++)  
	 {
	   up2 = arrp (aa2, up->nn, HIT2) ;
	   deltaPair = up2->dPair ;
	   if (deltaPair >= min || deltaPair <= -min)
	     ok = TRUE ;
	   else
	     up->score = 0 ;
	 }
       if (! ok)
	 hasGoodClone = -14 ;
     }
  
  ac_free (h) ;
  return hasGoodClone ;
} /* baPairFilterOne */

/*************************************************************************************/
/*
cat tmp/COUNT2/Rhs230/g.1.pairStats | head -28

Aligned_fragments   1075617
Compatible_pairs        856686
Non_compatible_pairs    170453
Orphans Any     48478
Orphans 1       44159
Orphans 2       2581
Orphans 3       733
Orphans 4       359
Orphans 5       243
Orphans 6       174
Orphans 7       108
Orphans 8       72
Orphans 9       49
Compatible_pairs_inside_gene    368159
Compatible_gene_extension       225174
Compatible_pairs_in_genome      156323
Compatible_pairs_in_spikeIn     5568
Compatible_pairs_in_mito        68827
Compatible_pairs_in_other_target        32635
Too_distant_inside_a_gene       0
Too_distant_on_genome   48317
Links_2_transcripts_of_a_gene   0
Links_gene_genome       11657
Links_2_genes   51924
Incompatible_topology   58555


Links_2_transcripts_of_a_gene   63010
Links_gene_genome       0
Links_2_genes   17793
Incompatible_topology   41333
*/
/*************************************************************************************/
/* filter by paired end coherence */
static int baPairFilter (BA *ba) 
{  
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (ba->outFileName, ".pairStats", FALSE, h) ; 
  Array aa = ba->hits ;
  Array aa2 = ba->hits2 ;
  HIT *up, *vp, *wp ;
  HIT2 *up2 = 0, *wp2 ;
  int ii, jj, j1, j2, n, nPairs = 0, nReads = 0, nR, aMax1 = arrayMax (aa), clone, tag ;
  int nGoodPairs = 0, nBadPairs = 0, nOrphans = 0 ;
  int nGeneGene = 0, nGenomeGenome = 0, nSpikeinSpikein = 0, nMitoMito = 0, nOtherOther = 0, nGeneGenome = 0 ;
  int nTooFarGenes = 0, nTooFarGenome = 0,  nCloseCisGenes = 0, nBadGeneGenome = 0, nTwoTranscripts = 0, nTwoGenes = 0,  nTwoCloseGenes = 0, nOtherBad = 0 ;
  int nGeneRrna = 0, nGeneMito = 0, nToSpikein = 0 ;
  KEYSET histo = keySetHandleCreate (h) ;
  KEYSET histoOrphan = keySetHandleCreate (h) ;
  KEYSET geneLinks = arrayHandleCreate (1000, HIT, h) ;

  arraySort (aa, baCloneOrder) ;
  if (NON_COMPATIBLE_PAIR != -15)
    messcrash ("adjust NON_COMPATIBLE_PAIR in wh/acedna.h to match the -14 used here for bad topology\n") ;
  for (ii = jj = 0, up = arrp (aa, 0, HIT), vp = up ; ii < arrayMax (aa) ; up++, ii++)
    { 
      int mult = 1 ;
      tag = up->tag ;
      up2 = arrp (ba->hits2, up->nn, HIT2) ;
      up2->badPair = 0 ;
      clone = up->clone ;
      if (! clone) 
	continue ;

      mult = up2->mult ; 
      nPairs += mult ; nReads += mult ; 
      /* select the range of hits for a given clone */
      for (nR = 1, jj = ii, vp = up ; jj < arrayMax (aa) && vp->clone == clone ; vp++, jj++) 
	if (vp->tag != tag)  nR++ ;
      if (nR == 1)	
	up2->unicity = 1 ; 
      if (nR > 1)
	{
	  nReads += mult ;
	  switch (baPairFilterOne (ba, ii, jj, geneLinks)) /* filter that clone */
	    {
	    case 1:
	      nGoodPairs += mult ;
	      nGeneGene += mult ;
	      break ;
	    case 2:
	      nGoodPairs += mult ;
	      nGenomeGenome += mult ;
	      break ;
	    case 3:
	      nGoodPairs += mult ;
	      nSpikeinSpikein += mult ;
	      break ;
	    case 4:
	      nGoodPairs += mult ;
	      nMitoMito += mult ;
	      break ;
	    case 5:
	      nGoodPairs += mult ;
	      nOtherOther += mult ;
	      break ;
	    case 6:         /* pair compatible on genome, anchored in a gene => dPair=-2 in gene*/
	      nGoodPairs += mult ;
	      nGeneGenome += mult ;
	      break ;
	    case -1:   /* gene to close-by gene */
	      break ;
	    case -2:   /* gene to close-by genome */
	      nGoodPairs += mult ;
	      nTwoCloseGenes += mult ;
	      break ;
	    case -3:   /* 2 reads in same gene but separated > 3 * ba->pair */
	      nGoodPairs += mult ;
	      nTooFarGenes += mult ;
	      break ;
	    case -10:
	      nGoodPairs += mult ;
	      nCloseCisGenes += mult ;
	      break ;
	    case -4:
	      nBadPairs += mult ;
	      nTooFarGenome += mult ;
	      up2->badPair = 1 ;
	      break ;
	    case -5:  /* synchronize to hack with snp.c and wiggle.c and sv.c : ventilate which accepts -2 and -5 as good snps */
	      nGoodPairs += mult ;
	      nTwoTranscripts += mult ;
	      break ;
	    case -6:   /* distant genes, otherwise they were labelled as +6 */
	      nBadPairs += mult ;
	      nTwoGenes += mult ;
	      up2->badPair = 1 ;
	      break ;
	    case -7:
	      nBadPairs += mult ;
	      nGeneRrna += mult ;
	      up2->badPair = 1 ;
	      break ;
	    case -8:
	      nBadPairs += mult ;
	      nGeneMito += mult ;
	      up2->badPair = 1 ;
	      break ;
	    case -11:
	      nBadPairs += mult ;
	      nToSpikein += mult ;
	      up2->badPair = 1 ;
	      break ;
	    case -9:     
	      nBadPairs += mult ;
	      nBadGeneGenome += mult ;
              up2->badPair = 1 ;
	      break ;
	    case -14:  /* synchronize to hack these reserved values with bestDa = -15, around line 971 and 1244 3195 and with wiggle.c and sv.c */
	    default:
	      nBadPairs += mult ;  /* bad topology */
	      nOtherBad += mult ;
	      up2->badPair = 1 ;
	      break ;
	    }	      
	}
      else
	{
	  int uu = 0 ;
	  nOrphans += mult ;
	  up2->badPair = 1 ;

	  for (j1 = ii, wp = arrp (aa, j1, HIT) ; j1 < jj ; j1++, wp++)
	    {
	      wp2 = arrp (aa2, wp->nn, HIT2) ;
	      if (! uu) uu = wp2->unicity ;
	      wp2->dPair = -1 ; /* ORPHAN */
	    }
	  if (uu < 0) uu = - uu ;
	  keySet (histoOrphan, uu)  += mult ;
	}
      ii = jj - 1 ; up = vp - 1 ; /* reposition at end of clone */
    }


  /* keep happy few */
  for (ii = jj = 0, up = arrp (aa, 0, HIT), vp = up ; ii < arrayMax (aa) ; up++, ii++)
    { 
      if (up->score)
	{
	  if (jj < ii) *vp = *up ;
	  vp++ ; jj++ ;
	}
    }
  arrayMax (aa) = jj ;

  if (0)   goto done ;
  /* export histogram of pair lengths */
  keySet (histo, 10000) = 0 ;
  for (n = ii = j1 = 0, up = arrp (aa, 0, HIT), tag = 0 ; ii < arrayMax (aa) ; up++, ii++)
    { 
      up2 = arrp (ba->hits2, up->nn, HIT2) ;
      if (up->tag == tag || up2->dPair <= 0 || class(up2->dPair) > 0) continue ;
      tag = up->tag ;
      jj = up2->dPair < 10000 ? up2->dPair : 10000 ;
      keySet (histo, jj)++ ;
      if (keySet (histo, jj) > n) { n = keySet (histo, jj) ; j1 = jj ; }
    }
  if (0) goto done ;
  /* cut off the histo at a position long enough to include 99.9% of all pairs */
  for (ii = keySetMax (histo) - 1 , j2 = 0 ; ii > j1 ; ii--)
    if (keySet (histo, ii) > n/1000) break ;
    else j2 += keySet (histo, ii) ;
  keySet (histo, ii+1) = j2 ;
  keySetMax (histo) = ii + 2 ;

  if (0)   goto done ;

  aceOutf (ao, "\nAligned_fragments\t%d", nPairs) ;
  if (0) aceOutf (ao, "\t//Good+Bad+Orphans\t%d", nGoodPairs + nBadPairs +  nOrphans) ;
  aceOutf (ao, "\nCompatible_pairs\t%d", nGoodPairs) ;
  aceOutf (ao, "\nNon_compatible_pairs\t%d", nBadPairs) ;

  if (ba->seqco)
    aceOutf (ba->seqco, "Both_Reads_Mapped\t%s\t2\t%ld\t%ld\n", ba->lane, 2*(nGoodPairs + nBadPairs), 2*nGoodPairs) ;

  aceOutf (ao, "\nOrphans\tAny\t%d", nOrphans) ;
  for (ii = 1 ; ii < keySetMax (histoOrphan) ; ii++)
    aceOutf (ao, "\nOrphans\t%d\t%d", ii, keySet (histoOrphan, ii)) ;

  if (0)   goto done ;

  aceOutf (ao, "\nCompatible_pairs_inside_gene\t%d", nGeneGene) ;
  aceOutf (ao, "\nCompatible_gene_extension\t%d", nGeneGenome + nTwoCloseGenes) ;
  aceOutf (ao, "\nCompatible_pairs_in_genome\t%d", nGenomeGenome) ;
  aceOutf (ao, "\nCompatible_pairs_in_spikeIn\t%d", nSpikeinSpikein) ;
  aceOutf (ao, "\nCompatible_pairs_in_mito\t%d", nMitoMito) ;
  aceOutf (ao, "\nCompatible_pairs_in_other_target\t%d", nOtherOther) ;
  aceOutf (ao, "\nCompatible_but_links_2_transcripts_of_a_gene\t%d", nTwoTranscripts) ;
  aceOutf (ao, "\nCompatible_but_too_distant_inside_a_gene\t%d", nTooFarGenes) ;
  aceOutf (ao, "\nCompatible_but_links_2_close_genes_in_cis\t%d", nCloseCisGenes) ;

  aceOutf (ao, "\nToo_distant_on_genome\t%d", nTooFarGenome) ;
  aceOutf (ao, "\nLinks_gene_to_distant_genome\t%d",  nBadGeneGenome) ;
  aceOutf (ao, "\nLinks_2_genes\t%d", nTwoGenes) ;
  aceOutf (ao, "\nLinks_to_rRNA\t%d", nGeneRrna) ; 
  aceOutf (ao, "\nLinks_to_mito\t%d", nGeneMito) ;
  aceOutf (ao, "\nLinks_to_spike_in\t%d", nToSpikein) ;
  aceOutf (ao, "\nIncompatible_topology\t%d", nOtherBad) ;
  

  aceOutf (ao, "\n\n#Histogram of pair lengths: %d (%.2f%%) compatible pairs, %d (%.2f%%) incompatible pairs, %d (%.2f%%) orphan reads\nLength"
	   , nGoodPairs, 100.0 * nGoodPairs/(nPairs ? nPairs : 1)
	   , nBadPairs, 100.0 * nBadPairs/(nPairs ? nPairs : 1)
	   , nOrphans, 100.0 * nOrphans/(nPairs ? nPairs : 1)
	   ) ;
  aceOutf (ao, "\n#Length") ;
  for (ii = 0 ; ii < keySetMax(histo) ; ii++)
    aceOutf(ao, "\t%d", ii) ;
  aceOutf (ao, "\n#Number of pairs") ;
  for (ii = 0 ; ii < keySetMax(histo) ; ii++)
    aceOutf(ao, "\t%d", keySet (histo, ii)) ;
  aceOutf (ao, "\n") ;

  if (arrayMax (geneLinks))
    {
      ACEOUT ao2 = aceOutCreate (ba->outFileName, ".geneLinks", ba->gzo, h) ;
      
      baTargetDict = ba->targetDict ;
      arraySort (geneLinks, baHitOrder) ;
      arrayCompress (geneLinks) ;
      for (ii = 0, up = arrayp (geneLinks, 0, HIT) ; ii < arrayMax (geneLinks) ; ii++, up++)
	{
	  if (up->score)
	    {
	      for (jj = ii + 1, vp = up + 1 ; jj < arrayMax (geneLinks) && vp->tag == up->tag && vp->gene == up->gene ; vp++, jj++)
		{ up->score += vp->score ; vp->score = 0 ; }
	      aceOutf (ao2, "%s\t%s\t%d\n", dictName (ba->targetDict, up->tag), dictName (ba->targetDict, up->gene), up->score) ;
	    }
	}
      ac_free (ao2) ;
    }

 done:
  /* report */
  fprintf (stderr, "Kept %d/%d hits = %.2f%%, found %d pairs, %d reads, %d (%.2f%%) compatible pairs, %d (%.2f%%) incompatible pairs, %d (%.2f%%) orphans\n"
	   , arrayMax (aa), aMax1, 100.0 * arrayMax (aa)/ (aMax1 > 0 ? aMax1 : 1)
	   , nPairs, nReads
	   , nGoodPairs, 100.0 * nGoodPairs/(nPairs ? nPairs : 1)
	   , nBadPairs, 100.0 * nBadPairs/(nPairs ? nPairs : 1)
	   , nOrphans, 100.0 * nOrphans/(nPairs ? nPairs : 1)
	   ) ;

  ac_free (h) ;

  baTargetDict = ba->targetDict ;
  arraySort (aa, baHitOrder) ;
  arrayCompress (aa) ;

  ac_free (ba->cloneDict) ;
  return nGoodPairs ;
} /* baPairFilter */

/*************************************************************************************/
/*************************************************************************************/

int bestAutoHitOrder (const void *a, const void *b)
{
  const AUTOHIT *up = (const AUTOHIT *)a, *vp = (const AUTOHIT *)b ;
  long long int n ;

  n = up->best - vp->best ;  if (n) return -n ;  
  n = up->target - vp->target ; if (n) return n ;

  return 0 ;
} /* bestAutoHitOrder */

/*************************************************************************************/
 /* Auto-hierarchy sort target by number of hits, keep only hits to top target */
static Array autoHierarchyOrderCounts = 0 ;

static int autoHierarchyOrder (const void *va, const void *vb)
{
  AUTOHIT *up = (AUTOHIT *)va ;
  AUTOHIT *vp = (AUTOHIT *)vb ;
  int nn = 0 ;

  nn = up->tag - vp->tag ; if (nn) return nn ;
  nn = up->score - vp->score ; if (nn) return -nn ; /* high score first */
  if (autoHierarchyOrderCounts)
    nn = 
      arrp (autoHierarchyOrderCounts, up->target, AUTOHIT)->tags - 
      arrp (autoHierarchyOrderCounts, vp->target, AUTOHIT)->tags ; 
  else
    nn = up->target - vp->target ; 
  
  return nn ;
} /* autoHierarchyOrder */

/*************************************************************************************/
/* HACK, may 24 2020, the table seem to eat up a bad value while sorting */

static BOOL baCheckHits (Array hits, int badTarget)
{
  AUTOHIT *up ;
  int jj, jjMax = arrayMax (hits) ;
  BOOL ok = FALSE ;

  for (jj = 0, up = arrp (hits, jj, AUTOHIT) ; jj < jjMax ; jj++, up++)
    if (up->target == badTarget)
      { ok = TRUE ; fprintf (stderr, "...... jj=%d up->target=%d\n",jj,up->target) ; }

  return ok ;
}
  
/*****************************/

static int baAutoHierarchyParseOne (BA *ba, ACEIN ai, Array counts, BOOL hierarchic)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0, jj = 0, jjMax, oldTag, oldScore, oldTarget ;
  AUTOHIT *up = 0 , *vp = 0 ;
  DICT *dict = dictHandleCreate (100000, h) ;
  Array hits = arrayHandleCreate (100000, AUTOHIT, h) ;

  while (aceInCard (ai))
    {  
      int tag = 0, score = 0, mult = 0, ln = 0, ali = 0, x1 = 0, x2 = 0, target, dummy = 0 ;
      char *ccp = aceInWord (ai) ;
      if (! ccp || ! *ccp || *ccp == '#')
	continue ;
      dictAdd (dict, ccp, &tag) ;
      if (
	  ! aceInStep (ai, '\t') || ! aceInInt (ai, &(score)) ||
	  ! aceInStep (ai, '\t') || ! aceInInt (ai, &(mult)) ||
	  ! aceInStep (ai, '\t') || ! aceInInt (ai, &(ln)) || ln < 5 ||
	  ! aceInStep (ai, '\t') || ! aceInInt (ai, &(ali)) ||
	  ! aceInStep (ai, '\t') || ! aceInInt (ai, &(x1)) ||
	  ! aceInStep (ai, '\t') || ! aceInInt (ai, &(x2)) 
	  )
	{
	  fprintf(stderr, "cannot read x2 or ln=%d<0 in %s\n",  ln, dictName (dict, tag)) ;
	  continue ;
	}
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ;
      if (! ccp || ! *ccp)
	{
	  fprintf(stderr, "cannot read target_class\n") ;
	  continue ;
	}
      if (ba->target_class)
	{
	  int cl = 0 ;
	  if (! dictFind (ba->target_classDict, ccp, &cl) || cl != ba->target_class)
	    continue ;
	}
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ;
      if (! ccp || ! *ccp)
	{
	  fprintf(stderr, "cannot read gene\n") ;
	  continue ;
	}
      if (! aceInStep (ai, '\t') || ! aceInInt (ai, &dummy))
	{
	  fprintf(stderr, "cannot read unicity\n") ;
	  continue ;
	}
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ;
      if (1 && ccp && !strncmp (ccp,"MRNA:",5))ccp += 5 ;
      if (! ccp || ! *ccp)
	return FALSE ;
      dictAdd (ba->targetDict, ccp, &(target)) ;

      if (up && tag == up->tag)
	{
	  if (score < up->score || target == up->target)
	    continue ; 
	}
      up = arrayp (hits, jj++, AUTOHIT) ;
      up->target = target ;
      vp = arrayp (counts, up->target, AUTOHIT) ;
      vp->target = up->target ;
      up->tag = tag ;
      up->score = score ;
      up->best = 0 ;
      up->tags = mult ;
      up->seqs = 1  ;
      up->ali = mult * ali ;
      up->ln = mult * ln ;
    }

  if (0 && baCheckHits (hits, 201))
    invokeDebugger () ;
  if (hierarchic)
    {
      autoHierarchyOrderCounts = counts ;
      arraySort (hits, autoHierarchyOrder) ;
      autoHierarchyOrderCounts = 0 ;
    }
  else
    arraySort (hits, autoHierarchyOrder) ;
  if (0 && baCheckHits (hits, 201))
    invokeDebugger () ;

  oldTag = oldScore = 0 ; jjMax = arrayMax (hits) ;
  for (jj = 0, up = arrp (hits, jj, AUTOHIT) ; jj < jjMax ; jj++, up++)
    {
      AUTOHIT *vp ;
      
      vp = arrayp (counts, up->target, AUTOHIT) ;
      vp->target = up->target ;
      if (up->tag != oldTag)
	oldScore = 0 ;
      if (up->score < oldScore)
	continue ;
      if (up->tag == oldTag &&
	  (
	   hierarchic ||
	   up->target == oldTarget
	   )
	  )
	continue ;
      oldTag = up->tag ; oldTarget = up->target ; oldScore = up->score ;
      if (hierarchic)
	{
	  vp->best += up->tags ;
	  vp->seqs += up->seqs ;
	  vp->ali += up->ali ;
	  vp->ln += up->ln ;
	}
      else
	{
	  vp->tags += up->tags ;
	}
    }
  ac_free (h) ;
  return nn ;
} /* baAutoHierarchyParseOne */

/*************************************************************************************/

static int baAutoHierarchyParse (BA *ba, Array counts, BOOL hierarchic)
{
  int nn = 0 ;
  int nInFile = 0 ;

  while (1) /* read all input files */
    {   /* do not use the baParse function, because we only want the read counts */
      ACEIN ai = 0 ;
      AC_HANDLE h = ac_new_handle () ;
      
      if (! ba->inFileList)
	ai = aceInCreate (ba->hitFileName, ba->gzi, h) ;
      else
	ai = baNextInFile (ba, nInFile++, h) ;
      if (ai) 
	nn += baAutoHierarchyParseOne (ba, ai, counts,  hierarchic) ;
      ac_free (h) ; 
      if (!ai || ! ba->inFileList)
	break ;
    }
 return nn ;
} /* baAutoHierarchyParse */

/*************************************************************************************/

static int baAutoHierarchy (BA *ba) 
{  
  AC_HANDLE h = ac_new_handle () ;

  ACEOUT ao = aceOutCreate (ba->outFileName, ".autoHierarchy", FALSE, h) ;
  AUTOHIT *vp ;
  int ii, iiMax ;
  Array counts = arrayHandleCreate (1000, AUTOHIT, h) ;
  
  if (! ba->targetDict)
    ba->targetDict = dictHandleCreate (1000, ba->h) ;
  baAutoHierarchyParse (ba, counts, FALSE) ; /* global hierarchic target counts */
  baAutoHierarchyParse (ba, counts, TRUE) ; /* global hierarchic target counts */
  
  arraySort (counts, bestAutoHitOrder) ;
   
  aceOutDate (ao, "###", "Winners take all: Bayesian read counts across similar targets") ;
  aceOut (ao, "## Column 1: Run name, as declared in the project\n") ;
  aceOut (ao, "## Column 2: Best mapping reads, reattributed only to the target with most hits\n") ;
  aceOut (ao, "## Column 3: Best mapping reads, including those mapping equally well to several targets\n") ;
  aceOut (ao, "## Column 4: Newness of the target, i.e. the percentage of winning hits in the target\n") ;
  aceOut (ao, "## Column 5: Number of distinct sequences hitting the target\n") ;
  aceOut (ao, "## Column 6: Average number of bases aligned on the target\n") ;
  aceOut (ao, "## Column 7: Average length of the reads hitting the target\n") ;
  aceOut (ao, "##   Each read is counted as 1 hit, a fragment pair count as 2 hits\n") ;
  aceOut (ao, "##   For each read, only the hits with highest score are kept\n") ;
  aceOut (ao, "##   In the first pass, the hits to each targets are counted, then the targets are sorted\n") ;
  aceOut (ao, "##   In the second pass, each read is attributted only to its most likely target,\n") ;
  aceOut (ao, "##   i.e. only to the target winner of the first pass, determined independently for each run\n") ;
  aceOut (ao, "##   In all columns excpet 4, the numbers are evealueted after reattributing\n") ;
	  
  aceOutf (ao, "# Run\tTarget\tWinning hits\tBest hits\tNewness\tDistinct sequences\tAverage aligment\tRead length\n") ;

  iiMax = arrayMax (counts) ;
  for (ii = 0, vp = arrp (counts, 0, AUTOHIT) ; ii < iiMax ; vp++, ii++)
    {
      if (vp->target)
	aceOutf (ao, "%s\t%s\t%ld\t%ld\t%.2f\t%ld\t%ld\t%ld\n"
		 , ba->run ? ba->run : "xxx"
		 , dictName (ba->targetDict, vp->target)
		 , vp->best
		 , vp->tags
		 , 100.0 * vp->best / (vp->tags ? vp->tags : 1)
		 , vp->seqs
		 , vp->ali/(vp->best ? vp->best : 1)
		 , vp->ln/(vp->best ? vp->best : 1)
		 ) ;
    }

  ac_free (h) ;
  return 0 ;
} /* baAutoHierarchy */

/*************************************************************************************/
/*************************************************************************************/
/* count or each read the sum of the scores and ali of non overlapping chains */
static void baCountBestGetSens (BA *ba, int ii, int cl
			       , int mult, int isFirstFragment
			       , int *sensp, int *sensFp, int *sensRp
				, int *nAp, int *nAFp, int *nARp, int *uup
			       ) 
{
  int jj, oldChain = 0, oldTarget = 0 ;
  int sens = 0, sensF = 0, sensR = 0 ; 
  int nA = 0, nAF = 0, nAR = 0, uu = 0, ln = 0 ;
  HIT *up, *vp ;
  HIT2 *vp2 ;
 
  up = arrp (ba->hits, ii, HIT) ;
  for (jj = ii, vp = up ; jj < arrayMax (ba->hits) && vp->tag == up->tag ; vp++, jj++)
    {
      int ali = 0 ;
      vp2 = arrp (ba->hits2, vp->nn, HIT2) ; 
      if (vp2->target_class != cl)
	continue ;
      
      if (! oldChain)
	{
	  uu = vp2->unicity ; 
	  ln = vp2->ln ;
	  if (uu == -2) 
	    uu = 11 ; 
	  else 
	    { 
	      if (uu < 0) 
		uu = -uu ;
	      if (uu > 10) 
		uu = 10 ; 
	    }
	}
      if (vp->chain != oldChain || oldTarget != vp->target)
	ali = mult * vp2->ali ;
      oldChain = vp->chain ;
      oldTarget = vp->target ; 
      nA += ali ;
      if (isFirstFragment > 0)
	{
	  nAF += ali ;
	  if (isFirstFragment * (vp2->a2 - vp2->a1) > 0)
	    {
	      sens |= 1 ;
	      sensF |= 1 ;
	    }
	  else
	    {
	      sens |= 2 ;
	      sensF |= 2 ;
	    }
	}
      else
	{
	  nAR += ali ;
	  if (isFirstFragment * (vp2->a2 - vp2->a1) > 0)
	    {
	      sens |= 1 ;
	      sensR |= 1 ;
	    }
	  else
	    {
	      sens |= 2 ;
	      sensR |= 2 ;
	    }
	}
    }
  if (uu != 1)  /* limit the calculatin of stranding to unique alignments */
    sens = sensF = sensR = 0 ;
  *sensp = sens ;
  *sensFp = sensF ;
  *sensRp = sensR ;

  if (nA > ln) nA = ln ;
  if (nAF > ln) nAF = ln ;
  if (nAR > ln) nAR = ln ;

  *nAp =  mult * nA ;
  *nAFp =  mult * nAF ;
  *nARp =  mult * nAR ;
  
  *uup = uu ;

  return ;
} /* baCountBestGetSens */

static int baCountBest (BA *ba) 
{  
  AC_HANDLE h = ac_new_handle () ;

  ACEOUT ao = aceOutCreate (ba->outFileName, ".count", FALSE, h) ;
  ACEOUT seqco = ba->seqco ;
  HIT *up, *vp ;
  HIT2 *up2, *vp2 ;
  int isFirstFragment ;
  int ii, jj, score = 0, cl, ln, uu ;
  long int mult = 0 ; /* to cast all products to long int */
  int MX = dictMax (ba->target_classDict) + 5 ;      /* target_class tables */
  BitSet hasClass = bitSetCreate (MX, h) ;
  BOOL hasTrClass ;
  long int nnS = 0, nnT = 0, nnA = 0, nnL = 0 , nnTp = 0, nnTm = 0, nnTamb = 0 ;
  long int nnSF = 0, nnTF = 0, nnAF = 0, nnLF = 0, nnTpF = 0, nnTmF = 0, nnTambF = 0 ;
  long int nnSR = 0, nnTR = 0, nnAR = 0, nnLR = 0, nnTpR = 0, nnTmR = 0, nnTambR = 0 ;
  long int nS[MX], nT[MX], nA[MX], nL[MX], nTp[MX], nTm[MX], nTamb[MX] ;
  long int ntrS = 0, ntrT = 0, ntrA = 0, ntrL = 0, ntrTp = 0, ntrTm = 0, ntrTamb = 0 ;
  long int hnS[MX], hnT[MX], hnA[MX], hnL[MX], hnTp[MX], hnTm[MX], hnTamb[MX] ;

  long int multiTargets[15*MX] ;
  memset (multiTargets, 0, sizeof(multiTargets)) ;

  memset (nS, 0, sizeof(nS)) ;
  memset (nT, 0, sizeof(nT)) ;
  memset (nA, 0, sizeof(nA)) ;
  memset (nL, 0, sizeof(nL)) ;
  memset (nTp, 0, sizeof(nTp)) ;
  memset (nTm, 0, sizeof(nTm)) ;
  memset (nTamb, 0, sizeof(nTamb)) ;

  memset (hnS, 0, sizeof(hnS)) ;
  memset (hnT, 0, sizeof(hnT)) ;
  memset (hnA, 0, sizeof(hnA)) ;
  memset (hnL, 0, sizeof(hnL)) ;
  memset (hnTp, 0, sizeof(hnTp)) ;
  memset (hnTm, 0, sizeof(hnTm)) ;
  memset (hnTamb, 0, sizeof(hnTamb)) ;
  
  for (ii = score = 0, up = arrp (ba->hits, 0, HIT) ; ii < arrayMax (ba->hits) ; up++, ii++)
    {
      int uu, bestCl ;
      up2 = arrayp (ba->hits2, up->nn, HIT2) ;
      /* up3 = arrayp (ba->hits3, up->nn, HIT3) ; */
      score = up->score ; 
      if (!score)
	continue ;
      ln = up2->ln ;
      cl = bestCl = -1 ; 
      /* global count */
      mult = fastcMultiplicity (dictName (ba->tagDict,up->tag), 0, 0) ; 
      isFirstFragment = 1 ;
      { 
	const char *ccp = dictName (ba->tagDict,up->tag) ;
	int i = strlen (ccp) ;
	ccp += i - 1 ;
	if (*ccp == '<')
	  isFirstFragment = -1 ;
      }
      
      /* jj loop on given tag, detect best score, best class, best unicity */
      hasClass = bitSetReCreate (hasClass, MX) ;
      hasTrClass = FALSE ;
      for (score = uu = 0, cl = -1, jj = ii, vp = up ; jj < arrayMax (ba->hits) && vp->tag == up->tag ; vp++, jj++)
	{
	  vp2 = arrp (ba->hits2, vp->nn, HIT2) ; 
	  bitSet (hasClass, vp2->target_class) ;
	  if (vp->score > score)
	    {
	      score = vp->score ;
	      vp2 = arrp (ba->hits2, vp->nn, HIT2) ;   
	      bestCl = vp2->target_class ; 
	      uu = vp2->unicity ; 
	      if (uu == -2) uu = 11 ; else { if (uu < 0) uu = -uu ; if (uu > 10) uu = 10 ; }
	    }
	}
      
      if (score == 0) /* this read is not aligned */
	{
	  ii = jj - 1 ; up = vp - 1 ;
	  continue ;
	}

      /* check the stranding in the best class */
      for (cl = 0 ; cl < MX ; cl++)
	{
	  int unA, unAF, unAR, sens, sensF, sensR ;
	  if (! bitt (hasClass, cl))
	    continue ;

	  baCountBestGetSens (ba, ii, cl, mult, isFirstFragment, &sens, &sensF, &sensR, &unA, &unAF, &unAR, &uu) ;
	  if (cl == bestCl)  /* count once per ii read */
	    {     /* add the data to the global count and to the hierarchic count */
	      nnS++ ; hnS[cl]++ ;
	      nnT += mult ; hnT[cl] += mult ;
	      nnL += mult * ln ;  hnL[cl] += mult * ln ;
	      nnA +=  unA ; hnA[cl] +=  unA ;
	      nnAF +=  unAF ; 
	      nnAR +=  unAR ; 
	      multiTargets [15 * (0) + uu] += mult ;

	      if (isFirstFragment > 0)
		{
		  nnSF++ ;
		  nnTF += mult ;
		  nnLF += mult * ln ;
		  multiTargets [15 * (1) + uu] += mult ;
		}
	      else
		{
		  nnSR++ ;
		  nnTR += mult ;
		  nnLR += mult * ln ;
		  multiTargets [15 * (2) + uu] += mult ;
		}
	      
	      switch (sens)
		{
		case 1:
		  nnTp += mult ; hnTp[cl] += mult ;
		  break ;
		case 2:
		  nnTm += mult ; hnTm[cl] += mult ;
		  break ;
		case 3:
		  nnTamb += mult ; hnTamb[cl] += mult ;
		  break ;
		}
	      
	      switch (sensF)
		{
		case 1:
		  nnTpF += mult ;
		  break ;
		case 2:
		  nnTmF += mult ;
		  break ;
		case 3:
		  nnTambF += mult ;
		  break ;
		}
	      
	      switch (sensR)
		{
		case 1:
		  nnTpR += mult ;
		  break ;
		case 2:
		  nnTmR += mult ;
		  break ;
		case 3:
		  nnTambR += mult ;
		  break ;
		}
	    }
	  if (1) /* for all classes, increase the non hierachic counts */
	    {
	      nS[cl]++ ;
	      nT[cl] += mult ;
	      nL[cl] += mult * ln ;
	      nA[cl] += unA ;
	      multiTargets [15 * (cl+2) + uu] += mult ;
	      switch (sens)
		{
		case 1:
		  nTp[cl] += mult ;
		  break ;
		case 2:
		  nTm[cl] += mult ;
		  break ;
		case 3:
		  nTamb[cl] += mult ;
		  break ;
		}	      
	    }
	  if (!hasTrClass && cl <= ba->lastTranscriptClass) /* for the T classes, increase the non hierachic anyTranscript */
	    { 
	      hasTrClass = TRUE ;
	      ntrS++ ;
	      ntrT += mult ;
	      ntrL += mult * ln ;
	      ntrA += unA ;
	      multiTargets [15 * (cl+2) + uu] += mult ;
	      switch (sens)
		{
		case 1:
		  ntrTp += mult ;
		  break ;
		case 2:
		  ntrTm += mult ;
		  break ;
		case 3:
		  ntrTamb += mult ;
		  break ;
		}	      
	    }
	}
      ii = jj - 1 ; up = vp - 1 ;
    }  
  
  aceOutf (ao, "\n#Date %s", timeShowNow ()) ; 
  aceOut (ao, "\n#HITS\tTarget\tSample\tSeq\tTag+\tTag-\tTag ambiguous\tTag\tAli\tLn\thSeq\thTag+\thTag-\thTag ambiguous\thTag\thAli\thLn") ;
  for (cl = 0 ; cl < MX ; cl++)
    {
      if (!nS[cl]) continue ;
      if (nS[cl])
	{
	  aceOutf (ao, "\nHITS\t%s", cl ? dictName (ba->target_classDict, cl) : "any") ;
	  aceOutf (ao, "\tany") ;
	  aceOutf (ao, "\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld"
		   , nS[cl], nTp[cl], nTm[cl], nTamb[cl], nT[cl], nA[cl], nL[cl] 
		   , hnS[cl], hnTp[cl], hnTm[cl], hnTamb[cl], hnT[cl], hnA[cl], hnL[cl] 
		   ) ;
	}
    } 
  aceOut (ao, "\nHITS\tAnnotated\tany") ;
  aceOutf (ao, "\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld"
	   , ntrS, ntrTp, ntrTm,  ntrTamb, ntrT, ntrA, ntrL 
	   , ntrS, ntrTp, ntrTm,  ntrTamb, ntrT, ntrA, ntrL 
	   ) ;
  aceOut (ao, "\nHITS\tany\tany") ;
  aceOutf (ao, "\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld"
	   , nnS, nnTp, nnTm,  nnTamb, nnT, nnA, nnL 
	   , nnS, nnTp, nnTm,  nnTamb, nnT, nnA, nnL 
	   ) ;
  aceOut (ao, "\nHITS\tany1\tany") ;
  aceOutf (ao, "\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld"
	   , nnSF, nnTpF, nnTmF,  nnTambF, nnTF, nnAF, nnLF 
	   , nnSF, nnTpF, nnTmF,  nnTambF, nnTF, nnAF, nnLF 
	   ) ;
  aceOut (ao, "\nHITS\tany2\tany") ;
  aceOutf (ao, "\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld"
	   , nnSR, nnTpR, nnTmR,  nnTambR, nnTR, nnAR, nnLR 
	   , nnSR, nnTpR, nnTmR,  nnTambR, nnTR, nnAR, nnLR 
	   ) ;
  aceOut (ao, "\n") ;  

  aceOutf (ao, "\n#Date %s, Multiplicity per target", timeShowNow ()) ; 
  aceOut (ao, "\n#MULT\tTarget") ;
  for (uu = 1 ; uu <= 11 ; uu++)
    aceOutf (ao, "\t%d", uu) ;

  for (cl = -2 ; cl < MX ; cl++)
    {
      long int n ; 
      if (cl > 0 && !nS[cl]) continue ;
      for (n = uu = 0 ; uu <= 11 ; uu++)
	n += multiTargets [ 15 * (cl+2) + uu] ;
      if (n)
	{
	  const char *ccp ;
	  switch (cl)
	    {
	    case -2: ccp = "any" ; break ;
	    case -1: ccp = "any1" ; break ;
	    case 0: ccp = "any2" ; break ;
	    default:
	      ccp = dictName(ba->target_classDict,cl) ;
	      break ;
	    }
	  aceOutf (ao, "\nMULT\t%s", ccp) ;
	  for (uu = 1 ; uu <= 11 ; uu++)
	    aceOutf (ao, "\t%ld", multiTargets [ 15 * (cl+2) + uu]) ;
	}
    }
  aceOutf (ao, "\n") ;
  

  if (seqco)
    {
      int u2 ;
      long int nnna = 0, nnnaF = 0, nnnaR = 0 ;
      for (uu = 1 ; uu <= 11 ; uu++)
	{
	  u2 =  uu < 11 ? uu : -2 ;
	  aceOutf (seqco, "MultiAli:%03d\t%s\t3\t%ld\t%ld\t%ld\n"
		   , u2 
		   , ba->lane
		   , multiTargets [ 15 * (0) + uu]
		   , multiTargets [ 15 * (1) + uu]
		   , multiTargets [ 15 * (2) + uu]
		   ) ;
	  if (u2 == -2) u2 = 2 ;
	  nnna  += u2 * multiTargets [ 15 * (0) + uu] ;
	  nnnaF += u2 * multiTargets [ 15 * (1) + uu] ;
	  nnnaR += u2 * multiTargets [ 15 * (2) + uu] ;
	}
	  
      aceOutf (seqco, "Reads_Mapped_per_strand\t%s\t6\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n"
	       , ba->lane
	       , nnTp + nnTmR - nnTpR , nnTm - nnTmR + nnTpR   /* in this format we do not complement the R2 read */
	       , nnTpF, nnTmF
	       , nnTmR, nnTpR
	       ) ;
      aceOutf (seqco, "Reads_Mapped\t%s\t3\t%ld\t%ld\t%ld\n"
	       , ba->lane
	       , nnT
	       , nnTF
	       , nnTR
	       ) ;
      aceOutf (seqco, "Reads_in_file\t%s\t3\t%ld\t%ld\t%ld\n"
	       , ba->lane
	       , nnT
	       , nnTF
	       , nnTR
	       ) ;
      aceOutf (seqco, "Aligned_bases\t%s\t3\t%ld\t%ld\t%ld\n"
	       , ba->lane
	       , nnA
	       , nnAF
	       , nnAR
	       ) ;
      aceOutf (seqco, "Cumulated_read_length\t%s\t3\t%ld\t%ld\t%ld\n"
	       , ba->lane
	       , nnL
	       , nnLF
	       , nnLR
	       ) ;
      aceOutf (seqco, "Alignments\t%s\t3\t%ld\t%ld\t%ld\n"
	       , ba->lane
	       , nnna
	       , nnnaF
	       , nnnaR
	       ) ;
    }

  ac_free (h) ;
  return 0 ;
} /* baCountBest */

/*************************************************************************************/
/* check that the mult=-2 reads really correspond to genes in antisense */
static int baCheckUnicity (BA *ba) 
{  
  int nnFound = 0, nnConfirmed = 0 ;
  int ii, jj, kk, score, tag, tcl, gene1, gene2, ok ;
  HIT *up, *vp, *wp ;
  HIT2 *up2, *vp2, *wp2 ;
  int ok2, ok3 ;

  for (ii = score = 0, up = arrp (ba->hits, 0, HIT) ; ii < arrayMax (ba->hits) ; up++, ii++)
    {
      if (! up->score)
	continue ;
      up2 = arrp (ba->hits2, up->nn, HIT2) ;
      if (up2->unicity != -2)
	continue ;
      nnFound++ ; ok2 = FALSE ;

      tag = up->tag ; tcl = up2->target_class ;
      for (tcl = 1 ; tcl <= dictMax(ba->target_classDict) ; tcl++)
	{
	  if (strncmp (dictName (ba->target_classDict, tcl) + 1, "T_",2))
	    continue ; /* only check the ET_av KT_RefSeq etc transcript classes */
	  
	  ok3 = 0 ;
	  for (jj = ii, vp = up, ok = 2 ; ok == 2 && jj < arrayMax (ba->hits) && vp->tag == tag ; vp++, jj++)
	    {
	      vp2 = arrp (ba->hits2, vp->nn, HIT2) ;
	      if (vp2->target_class != tcl)
		continue ;
	      if (vp2->unicity != -2)
		continue ;
	      gene1 = vp->gene ; ok3 = 1 ;
	      for (kk = jj + 1, wp = vp + 1 ; kk < arrayMax (ba->hits) && wp->tag == tag ; wp++, kk++)
		{
		  wp2 = arrp (ba->hits2, wp->nn, HIT2) ;
		  if (wp2->target_class != tcl)
		    continue ;
		  if (wp2->unicity != -2)
		    continue ;
		  gene2 = wp->gene ;
		  if (gene1 == gene2) /* jump the alternative transcripts cases, wait for the second gene */
		    continue ;
		  ok3++ ;
		  if (ba->gene2map &&
		      gene1 < arrayMax (ba->gene2map) &&
		      gene2 < arrayMax (ba->gene2map) 
		      )
		    {
		      HIT *g1 = arrp (ba->gene2map, gene1, HIT) ;
		      HIT *g2 = arrp (ba->gene2map, gene2, HIT) ;
			  
		      if (g1->target == g2->target && 
			  (  /* genes must be antisense */
			   (
			    g1->x1 < g1->x2 && g2->x1 > g2->x2 &&
			    g1->x1 < g2->x1 && g1->x2 > g2->x2
			    ) ||
			   (
			    g2->x1 < g2->x2 && g1->x1 > g1->x2 &&
			    g2->x1 < g1->x1 && g2->x2 > g1->x2
			    )
			     )
			  )
			{ ok = -2 ; ok2 |= 1 ; }
		    }
		}

	      /* reset all cases for this target class */
	      if (ok != -2)
		for (jj = ii, vp = up ; jj < arrayMax (ba->hits) && vp->tag == tag ; vp++, jj++) 
		  {
		    vp2 = arrp (ba->hits2, vp->nn, HIT2) ;
		    if (vp2->target_class == tcl && vp2->unicity == -2)
		      vp2->unicity = ok3 ; /* may be 2 (failed to find antisense genes or -2 found antisense genes */ 
		  }
	    }
	}
      if (ok2)
	nnConfirmed++ ;
      /* this tag has been verified, reposition on last line with same tag */
      for (jj = ii, vp = up ; jj < arrayMax (ba->hits) && vp->tag == tag ; vp++, jj++) ;
      ii = jj - 1 ; up = vp - 1 ; 
    }
  fprintf (stderr, "--- baCheckUnicity: Found %d u==-2 reads, confirmed %d\n", nnFound, nnConfirmed) ;

  return nnConfirmed ;
} /* baCheckUnicity */

/*************************************************************************************/
/* export only the best hits */
static int baExportBest (BA *ba) 
{  
  AC_HANDLE h = ac_new_handle () ;

  HIT *up ;
  int ii, nn = 0 ;
  ACEOUT ao = aceOutCreate (ba->outFileName, ".hits", ba->gzo, h) ;

  for (ii = 0, up = arrp (ba->hits, 0, HIT) ; ii < arrayMax (ba->hits) ; up++, ii++)
    {
      if (up->score)
	{
	  nn++ ;
	  baExportOneHit (ao, ba, up) ;
	}  
    }
  ac_free (h) ;
  fprintf (stderr, "Export done %s\n", timeShowNow ()) ;
  return nn ;
} /* baExportBest */

/*************************************************************************************/
/* export only the best hits to the transcripts listed in the sigList */
static int baExportSigHits (BA *ba) 
{  
  AC_HANDLE h = ac_new_handle () ;

  DICT *dict = ba->targetDict ;
  HIT *up ; HIT2 *up2 ;
  int ii, score = 0, nn = 0 , nne = 0, av ;
  char *cp ;
  Array hits = ba->hits ;
  Array hits2 = ba->hits2 ;
  int iiMax = arrayMax (hits) ;
  ACEOUT ao = aceOutCreate (ba->outFileName, ".sig_hits", ba->gzo, h) ;
  ACEIN ai = aceInCreate (ba->sigTargetFile, FALSE, h) ;
  KEYSET ks = keySetHandleCreate (h) ;

  dictFind (ba->target_classDict, "ET_av", &av) ;
  while (aceInCard (ai))
    {
      cp = aceInWord (ai) ;
      if (! cp || *cp == '#') continue ;
      if (dictFind (dict, cp, &ii))
	keySet (ks, nn++) = ii ;
    }
  ac_free (ai) ;

  keySetSort (ks) ; keySetCompress (ks) ; nn = keySetMax (ks) ;
  for (ii = score = 0, up = arrp (hits, 0, HIT) ; ii < iiMax ; up++, ii++)
    {
      up2 = arrp (hits2, up->nn, HIT2) ; 
      if (up->score && up2->target_class == av &&
	  keySetFind (ks, up->target, 0) && 
	  100 * up2->ali > up2->ln
	  /* &&	  (! ba->pair || up->badPair) */
	  )
	{
	  nne++ ;
	  baExportOneHit (ao, ba, up) ;
	}  
    }
  ac_free (h) ;
  fprintf (stderr, "Sig_export done %s  found %d targets exported %d hits\n", timeShowNow (), nn, nne) ;
  return nn ;
} /*  baExportSigHits */

/*************************************************************************************/

static int baExportVenn (BA *ba)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (ba->outFileName, ".Venn", 0, h) ;
  DICT *dict = dictHandleCreate (10000, h) ;
  KEYSET ks = keySetHandleCreate (h) ;
  KEYSET ks1 = keySetHandleCreate (h) ;
  KEYSET ks2 = keySetHandleCreate (h) ;
  KEYSET ks3 = keySetHandleCreate (h) ;
  int k, t, t1, t2, ii, n, nn = 0, tag, MX ;
  long unsigned int z ;
  const char *ccp ;
  HIT *up ;
  HIT2 *up2 ;

  for (ii = 0, up = arrp (ba->hits, 0, HIT) ; ii < arrayMax (ba->hits) ; up++, ii++)
    {
      if (! up->score)
	continue ;
      up2 = arrayp (ba->hits2, up->nn, HIT2) ;
      n = up2->mult ;
      for (z = 0, tag = up->tag ; ii < arrayMax (ba->hits) && tag == up->tag && up->score ; ii++, up++)
	{
	  up2 = arrayp (ba->hits2, up->nn, HIT2) ;
	  if (up2->target_class < 32)
	    z = z | (1 << up2->target_class) ;
	}
      up-- ; ii-- ; /* restore */
      dictAdd (dict, messprintf("%ld",z), &k) ;
      keySet (ks, k) += n ;
      nn++ ;
    }

  MX = dictMax(ba->target_classDict) + 1 ;
  for (k = 1 ; k <= dictMax (dict) ; k++)
    {
      int j ;

      ccp = dictName(dict,k) ;
      n = keySet (ks , k) ;
      if (n < 1) continue ;
      aceOutf (ao, "%d\tDetails", n) ;
      sscanf (ccp, "%ld", &z) ;
      for (j = 0, t = 1 ; t<= dictMax (ba->target_classDict) ; t++)
	if (z & (((long unsigned int) 1) << t))
	  {
	    aceOutf (ao, "%s%s", j++ ? ";" : "\t", dictName (ba->target_classDict, t)) ;
	    keySet(ks1,t) += n ;

	    for (t1 = t+1 ; t1 <= dictMax (ba->target_classDict) ; t1++)
	      if (z & (((long unsigned int) 1) << t1))
		{
		  keySet(ks2,MX*t+t1) += n ;
		  for (t2 = t1+1 ; t2 <= dictMax (ba->target_classDict) ; t2++)
		    if (z & (((long unsigned int) 1) << t2))
		      keySet(ks3,MX*MX*t+MX*t1+t2) += n ;
		}

	  }
      aceOutf (ao, "\n") ;
    }

  for (t = 1 ; t<= dictMax (ba->target_classDict) ; t++)
    {
      n = keySet(ks1,t) ;
      if (n)
	aceOutf (ao, "%d\tSinglet\t%s\n", keySet(ks1,t)
		 , dictName (ba->target_classDict, t)
		 ) ;
    }
  for (t = 1 ; t<= dictMax (ba->target_classDict) ; t++)
    for (t1 = t+1 ; t1 <= dictMax (ba->target_classDict) ; t1++)
      {
	n = keySet(ks2,MX*t+t1) ;
	if (n)	  
	  aceOutf (ao, "%d\tDoublet\t%s;%s\n", n
		   , dictName (ba->target_classDict, t)
		   , dictName (ba->target_classDict, t1)
		   ) ;
      }
  for (t = 1 ; t<= dictMax (ba->target_classDict) ; t++)
    for (t1 = t+1 ; t1 <= dictMax (ba->target_classDict) ; t1++)
      for (t2 = t1+1 ; t2 <= dictMax (ba->target_classDict) ; t2++)
	{
	  n = keySet(ks3,MX*MX*t+MX*t1+t2) ;
	  if (n)
	    aceOutf (ao, "%d\tTriplet\t%s;%s;%s\n", n
		     , dictName (ba->target_classDict, t)
		     , dictName (ba->target_classDict, t1)
		     , dictName (ba->target_classDict, t2)
		     ) ;
	}
  
  ac_free (h) ;
  fprintf (stderr, "Venn-Export done %s\n", timeShowNow ()) ;
  return nn ;
} /* baExportVenn */

/*************************************************************************************/
/* export the error types in the best target
 * to restrict to a single target, use the -target_class param and filter
 */
static int baExportErrorProfile (BA *ba) 
{  
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT seqco = ba->seqco ;
  DICT *tagDict = ba->tagDict ;
  DICT *errDict = dictHandleCreate (1000, h) ;
  Array errors1 = arrayHandleCreate (1024, long int, h) ;
  Array errors2 = arrayHandleCreate (1024, long int, h) ;
  Array errPos1 = arrayHandleCreate (1024, long int, h) ;
  Array errPos2 = arrayHandleCreate (1024, long int, h) ;
  Array errPosType1 = arrayHandleCreate (60024, long int, h) ;
  Array errPosType2 = arrayHandleCreate (60024, long int, h) ;
  Array FB1 = arrayHandleCreate (1024, long int, h) ;
  Array FB2 = arrayHandleCreate (1024, long int, h) ;
  Array LB1 = arrayHandleCreate (1024, long int, h) ;
  Array LB2 = arrayHandleCreate (1024, long int, h) ;
  Array LN1 = arrayHandleCreate (1024, long int, h) ;
  Array LN2 = arrayHandleCreate (1024, long int, h) ;
  Array NNPos1 = arrayHandleCreate (1024, long int, h) ;
  Array NNPos2 = arrayHandleCreate (1024, long int, h) ;
  HIT *up ;
  int ii, jj, kk, score = 0, oldTag = 0, mult, nn = 0, err, pos, u, prefix, type, typeMax  ;
  long int ali1 = 0, ali2 = 0, nn1, nn2, nnN1, nnN2, n1 ;
  ACEOUT ao = aceOutCreate (ba->outFileName, ".errorProfile", ba->gzo, h) ;
  char buf[4096], *cp, *cq, *cr ;
  char *atgc = "atgcrymksw" ;
  const char *ccp ;
  int jTr, tr1 = 0, tr2 = 0, genome = 0 ; 
    
  dictFind (ba->target_classDict, "DT_magic", &tr1) ;      /* first transciptome class */
  dictFind (ba->target_classDict, "PT_tRNA", &tr2) ;    /* last */
  dictFind (ba->target_classDict, "Z_genome", &genome) ; /* genome class */

  
  aceOutf (ao, "# %s : error profile start\n", timeShowNow ()) ;
  for (ii = 0 ; ii < 10 ; ii++)
    for (jj = 0 ; jj < 10 ; jj++)
      {
	if (ii == jj) continue ;
	buf[0] = atgc[ii] ; 
	buf[1] = '>' ;
	buf[2] = atgc[jj] ; 
	buf[3] = 0 ;
	dictAdd (errDict, buf, 0) ;
      }

  for (ii = 0 ; ii < 10 ; ii++)
    {
      buf[0] = '+' ;
      buf[1] = atgc[ii] ; 
      buf[2] = 0 ;
      dictAdd (errDict, buf, 0) ;
    }
  for (ii = 0 ; ii < 10 ; ii++)
    {
      buf[0] = '*' ;
      buf[1] = '+' ;
      buf[2] = atgc[ii] ; 
      buf[3] = 0 ;
      dictAdd (errDict, buf, 0) ;
    }
  for (ii = 0 ; ii < 10 ; ii++)
    {
      buf[0] = '-' ;
      buf[1] = atgc[ii] ; 
      buf[2] = 0 ;
      dictAdd (errDict, buf, 0) ;
    }
  typeMax = dictMax (errDict) ;
  for (ii = 0 ; ii < 10 ; ii++)
    {
      buf[0] = '*' ;
      buf[1] = '-' ;
      buf[2] = atgc[ii] ; 
      buf[3] = 0 ;
      dictAdd (errDict, buf, 0) ;
    }


  for (ii = 0 ; ii < 10 ; ii++)
    for (jj = 0 ; jj < 10 ; jj++)
      {
	buf[0] = '+' ;
	buf[1] = '+' ;
	buf[2] = atgc[ii] ; 
	buf[3] = atgc[jj] ; 
	buf[4] = 0 ;
	dictAdd (errDict, buf, 0) ;
      }
  for (ii = 0 ; ii < 10 ; ii++)
    for (jj = 0 ; jj < 10 ; jj++)
      {
	buf[0] = '*' ;
	buf[1] = '+' ;
	buf[2] = '+' ;
	buf[3] = atgc[ii] ; 
	buf[4] = atgc[jj] ; 
	buf[5] = 0 ;
	dictAdd (errDict, buf, 0) ;
      }
  for (ii = 0 ; ii < 10 ; ii++)
    for (jj = 0 ; jj < 10 ; jj++)
      {
	buf[0] = '-' ;
	buf[1] = '-' ;
	buf[2] = atgc[ii] ; 
	buf[3] = atgc[jj] ; 
	buf[4] = 0 ;
	dictAdd (errDict, buf, 0) ;
      }
  for (ii = 0 ; ii < 10 ; ii++)
    for (jj = 0 ; jj < 10 ; jj++)
      {
	buf[0] = '*' ;
	buf[1] = '-' ;
	buf[2] = '-' ;
	buf[3] = atgc[ii] ; 
	buf[4] = atgc[jj] ; 
	buf[5] = 0 ;
	dictAdd (errDict, buf, 0) ;
      }


  for (ii = 0 ; ii < 4 ; ii++)
    for (jj = 0 ; jj < 4 ; jj++)
     for (kk = 0 ; kk < 4 ; kk++)
      {
	buf[0] = '+' ;
	buf[1] = '+' ;
	buf[2] = '+' ;
	buf[3] = atgc[ii] ; 
	buf[4] = atgc[jj] ; 
	buf[5] = atgc[kk] ; 
	buf[6] = 0 ;
	dictAdd (errDict, buf, 0) ;
      }
  for (ii = 0 ; ii < 4 ; ii++)
    for (jj = 0 ; jj < 4 ; jj++)
     for (kk = 0 ; kk < 4 ; kk++)
      {
	buf[0] = '*' ;
	buf[1] = '+' ;
	buf[2] = '+' ;
	buf[3] = '+' ;
	buf[4] = atgc[ii] ; 
	buf[5] = atgc[jj] ; 
	buf[6] = atgc[kk] ; 
	buf[7] = 0 ;
	dictAdd (errDict, buf, 0) ;
      }
  for (ii = 0 ; ii < 4 ; ii++)
    for (jj = 0 ; jj < 4 ; jj++)
     for (kk = 0 ; kk < 4 ; kk++)
      {
	buf[0] = '-' ;
	buf[1] = '-' ;
	buf[2] = '-' ;
	buf[3] = atgc[ii] ; 
	buf[4] = atgc[jj] ; 
	buf[5] = atgc[kk] ; 
	buf[6] = 0 ;
	dictAdd (errDict, buf, 0) ;
      }
  for (ii = 0 ; ii < 4 ; ii++)
    for (jj = 0 ; jj < 4 ; jj++)
     for (kk = 0 ; kk < 4 ; kk++)
      {
	buf[0] = '*' ;
	buf[1] = '-' ;
	buf[2] = '-' ;
	buf[3] = '-' ;
	buf[4] = atgc[ii] ; 
	buf[5] = atgc[jj] ; 
	buf[6] = atgc[kk] ; 
	buf[7] = 0 ;
	dictAdd (errDict, buf, 0) ;
      }

  for (ii = score = 0, up = arrp (ba->hits, 0, HIT) ; ii < arrayMax (ba->hits) ; up++, ii++)
    {
      int x1, x2, ali ; /* first last base aligned */
      HIT *vp ;
      HIT2 *up2 ;

      score = up->score ;
      if (! score)
	continue ;
      if (up->tag == oldTag)   /* count each tag only once, but explore all exons */
	continue ;
      up2 = arrp (ba->hits2, up->nn, HIT2) ;

      if (up2->target_class < tr1) 
	jTr = 2 ;
      else if (up2->target_class <= tr2) 
	jTr = 4 ;
      else if (up2->target_class == genome)
	{
	  if (ba->RNA_seq) continue ;
	}
      else                  /* reject non-transcriptome/non-genome cases */
	continue ;

      u = up2->unicity ; 
      if (u * u > 1 && ! (u == -2 && jTr == 6))           /* unique, non-unique, accepting the case -2 = single site two genes */
	continue ;

      oldTag = up->tag ;
      x1 = up->x1 ; x2 = up->x2 ;
      ali = up2->ali ; /* already cumulated across all discontinuous exons */

      ccp = dictName (tagDict,up->tag) ;
      mult = fastcMultiplicity (ccp, 0, 0) ;

      for (jj = ii, vp = up ; vp->tag == oldTag && vp->score == score && vp->target == up->target ; vp++, jj++) 
	{
	  if (vp->x1 < x1) x1 = vp->x1 ; /* first base */
	  if (vp->x2 > x2) x2 = vp->x2 ; /* last base */
	}

      /* find the orientation */
      ccp += strlen (ccp) - 1 ;
      prefix = 0 ;
      if (*ccp == '<') prefix = 1 ; 

      /* register the ali, fisrt and last base */
      if (prefix == 0)
	{
	  ali1 += mult * ali ;
	  array (FB1, x1, long int) += mult ;
	  array (LB1, x2, long int) += mult ;
	  array (LN1, ali, long int) += mult ;
	}
      else
	{
	  ali2 += mult * ali ;
	  array (FB2, x1, long int) += mult ;
	  array (LB2, x2, long int) += mult ;
	  array (LN2, ali, long int) += mult ;
	}


      /* look for errors */
      for (jj = ii, vp = up ; vp->tag == oldTag && vp->score == score && vp->target == up->target ; vp++, jj++) 
	{
	  /* look for errors */
	  HIT3 *vp3 = arrp (ba->hits3, vp->nn, HIT3) ;

	  nn++ ;
	  
	  if (vp3->errTag && vp3->nErr + vp3->nN >= 0)
	    {
	      int xE = vp3->nErr , xN = vp3->nN ; /* 2012_08_19 rarely the errTag is longer than nErr or nN, fix that */
	      strncpy (buf, dictName (ba->errDict, vp3->errTag), 4095) ;
	      cp = buf ;
	      while (cp)
		{
		  cq = strstr (cp, ",") ;
		  if (cq) *cq = 0 ;
		  cr = strstr (cp, ":") ;
		  if (cr)
		    {
		      *cr = 0 ;
		      pos = atoi (cp) ;
		      type = 0 ;
		      cp = cr + 1 ;
		      if (*cp && ! strstr (cp, "n") && ! strstr (cp, "o"))
			{
			  dictAdd (errDict, cp, &err) ;
			  if (*cp == '*')
			    dictFind (errDict, cp+1, &type) ;
			  else
			    type = err ;
			  if (type > typeMax)
			    {
			      if (cp[0] == '+' || (cp[0] == '*' && cp[1] == '+')) type = typeMax + 1 ; 
			      if (cp[0] == '-' || (cp[0] == '*' && cp[1] == '-')) type = typeMax + 2 ; 
			    }
			  if (prefix == 0)
			    {
			      if (1 || xE-- > 0)
				{
				  array (errors1, err, long int) += mult ;
				  array (errPos1, pos, long int) += mult ;
				  array (errPosType1, type + (typeMax + 3) * pos, long int) += mult ;
				}
			    }
			  else
			    {
			      if (1 || xE-- > 0)
				{
				  array (errors2, err, long int) += mult ;
				  array (errPos2, pos, long int) += mult ;
				  array (errPosType2, type +  (typeMax + 3) * pos, long int) += mult ;
				}
			    }
			} 
		      else
			{
			   if (1 || xN-- > 0)
			     {
			       if (prefix == 0)
				 array (NNPos1, pos, long int) += mult ;
			       else
				 array (NNPos2, pos, long int) += mult ;
			     }
			} 
		    }
		  cp = (cq ? cq + 1 : 0) ;
		}
	    }
	}  
    }
  aceOutf (ao, "# %s : error profile done\n", timeShowNow ()) ;   
  if (ali1 > 0) aceOutf (ao, "ERRALI\tf1\t%ld\tkb\n",ali1/1000) ;
  if (ali2 > 0) aceOutf (ao, "ERRALI\tf2\t%ld\tkb\n",ali2/1000) ;
  

  if (1)
    {
      long int nMis1, nMis2, nSub1, nSub2, nDel1, nDel2, nIns1, nIns2, nTrv1, nTrv2, nTrs1, nTrs2 ;

      nMis1 = nMis2 = nSub1 = nSub2 = nTrs1 = nTrs2 = nTrv1 = nTrv2 = nDel1 = nDel2 = nIns1 = nIns2 = 0 ;

      for (err = 1 ; err <= dictMax (errDict) ; err++)
	{
	  long int n1 = 0, n2 = 0 ;
	  
	  if (array (errors1, err, long int))
	    aceOutf (ao, "ERRS\tf1\t%d\t%s\t%ld\n", err, dictName (errDict, err), n1 = array (errors1, err, long int)) ; 
	  if (array (errors2, err, long int))
	    aceOutf (ao, "ERRS\tf2\t%d\t%s\t%ld\n", err, dictName (errDict, err), n2 = array (errors2, err, long int)) ; 
	  if (seqco && n1 + n2 > 0)
	    {
	      char buf[32] ;
	      int i ;
	      const char *ccp = dictName (errDict, err) ;
	      if (*ccp == '*') ccp++ ;
	      
	      if (*ccp == '+')
		{ 
		  strcpy (buf,"Ins:") ; i = 4 ; while (*ccp == '+') ccp++ ; ccp-- ; while (*++ccp) buf[i++] = ace_upper(*ccp);
		  buf[i++] = 0 ;
		  nIns1 += n1 ; nMis1 += n1 ;
		  nIns2 += n2 ; nMis2 += n2 ;
		}
	      else 
		{
		  if (*ccp == '-') 
		    { 
		      strcpy (buf,"Del:") ; i = 4 ; while (*ccp == '-')ccp++ ; ccp-- ; while (*++ccp) buf[i++] = ace_upper(*ccp); 
		      buf[i++] = 0 ;
		      nDel1 += n1 ; nDel2 += n2 ;
		      nMis1 += n1 ; nMis2 += n2 ;
		    }
		  else 
		    { 
		      strcpy (buf,"Sub:") ; i = 4 ; buf[i++] = ace_upper(ccp[0]) ; buf[i++] = ace_upper(ccp[2]) ;
		      buf[i++] = 0 ;
		      if (
			  (ccp[0] == 'a' && ccp[2] == 'g') ||
			  (ccp[0] == 'g' && ccp[2] == 'a') ||
			  (ccp[0] == 't' && ccp[2] == 'c') ||
			  (ccp[0] == 'c' && ccp[2] == 't')
			  )
			{ nTrs1 += n1 ; nTrs2 += n2 ; }
		      else
			{ nTrv1 += n1 ; nTrv2 += n2 ; }
		      nSub1 += n1 ; nSub2 += n2 ;
		      nMis1 += n1 ; nMis2 += n2 ;
		    }
		}
	      
	      aceOutf (seqco, "%s\t%s\t3\t%ld\t%ld\t%ld\n", buf, ba->lane, n1 + n2, n1, n2) ;
	    }
	}
      if (seqco)
	{
	  aceOutf (seqco, "Mismatches\t%s\t3\t%ld\t%ld\t%ld\n",    ba->lane, nMis1 + nMis2, nMis1, nMis2) ;
	  aceOutf (seqco, "Substitutions\t%s\t3\t%ld\t%ld\t%ld\n", ba->lane, nSub1 + nSub2, nSub1, nSub2) ;
	  aceOutf (seqco, "Transitions\t%s\t3\t%ld\t%ld\t%ld\n",   ba->lane, nTrs1 + nTrs2, nTrs1, nTrs2) ;
	  aceOutf (seqco, "Transversions\t%s\t3\t%ld\t%ld\t%ld\n",  ba->lane, nTrv1 + nTrv2, nTrv1, nTrv2) ;
	  aceOutf (seqco, "Insertions\t%s\t3\t%ld\t%ld\t%ld\n",    ba->lane, nIns1 + nIns2, nIns1, nIns2) ;
	  aceOutf (seqco, "Deletions\t%s\t3\t%ld\t%ld\t%ld\n",     ba->lane, nDel1 + nDel2, nDel1, nDel2) ;
	}
    }

  if (1)
    {
      int mx1 =  arrayMax (errPos1) ;
      int mx2 =  arrayMax (errPos2) ;
      int mx = mx1 > mx2 ? mx1 : mx2 ;
      long int n1, n2 ;

      for (pos = 1, nn1 = nn2 = 0 ; pos < mx ; pos++)
	{
	  n1 = pos < mx1 ? array (errPos1, pos, long int) : 0 ;
	  n2 = pos < mx2 ? array (errPos2, pos, long int) : 0 ;
	  nn1 += n1 ;
	  nn2 += n2 ;
	  if (n1 > 0)
	    aceOutf (ao, "POS\tf1\t%d\t%ld\n", pos, n1) ; 
	  if (n2 > 0)
	    aceOutf (ao, "POS\tf2\t%d\t%ld\n", pos, n2) ; 
	  if (seqco && n1 + n2 > 0)
	    aceOutf (seqco, "Err_pos:%04d\t%s\t3\t%ld\t%ld\t%ld\n", pos - 1, ba->lane, n1 + n2, n1, n2) ; /* zero-based err_pos */
	}
    }
  if (1)
    {
      int mx1 =  arrayMax (NNPos1) ;
      int mx2 =  arrayMax (NNPos2) ;
      int mx = mx1 > mx2 ? mx1 : mx2 ;
      long int n1, n2 ;

      for (pos = 1, nnN1 = nnN2 = 0 ; pos < mx ; pos++)
	{
	  n1 = pos < mx1 ? array (NNPos1, pos, long int) : 0 ;
	  n2 = pos < mx2 ? array (NNPos2, pos, long int) : 0 ;
	  nnN1 += n1 ;
	  nnN2 += n2 ;
	  if (n1 > 0)
	    aceOutf (ao, "NNPOS\tf1\t%d\t%ld\n", pos, n1) ; 
	  if (n2 > 0)
	    aceOutf (ao, "NNPOS\tf1\t%d\t%ld\n", pos, n1) ; 
	  if (seqco && n1 + n2 > 0)
	    aceOutf (seqco, "Ambiguous_pos:%04d\t%s\t3\t%ld\t%ld\t%ld\n", pos, ba->lane, n1 + n2, n1, n2) ;
	}
    }


  for (pos = 1 ; pos < arrayMax (FB1) ; pos++)
    if (array (FB1, pos, long int) > 0)
      aceOutf (ao, "FB\tf1\t%d\t%ld\n", pos, array (FB1, pos, long int)) ; 
  for (pos = 1 ; pos < arrayMax (LB1) ; pos++)
    if (array (LB1, pos, long int) > 0)
      aceOutf (ao, "LB\tf1\t%d\t%ld\n", pos, array (LB1, pos, long int)) ; 
  for (pos = 1 ; pos < arrayMax (LN1) ; pos++)
    if (array (LN1, pos, long int) > 0)
      aceOutf (ao, "LN\tf1\t%d\t%ld\n", pos, array (LN1, pos, long int)) ; 

  for (pos = 1 ; pos < arrayMax (FB2) ; pos++)
    if (array (FB2, pos, long int) > 0)
      aceOutf (ao, "FB\tf2\t%d\t%ld\n", pos, array (FB2, pos, long int)) ; 
  for (pos = 1 ; pos < arrayMax (LB2) ; pos++)
    if (array (LB2, pos, long int) > 0)
      aceOutf (ao, "LB\tf2\t%d\t%ld\n", pos, array (LB2, pos, long int)) ; 
  for (pos = 1 ; pos < arrayMax (LN2) ; pos++)
    if (array (LN2, pos, long int) > 0)
      aceOutf (ao, "LN\tf2\t%d\t%ld\n", pos, array (LN2, pos, long int)) ; 

  if (seqco)
    {
      int mx1 = arrayMax(LN1) ;
      int mx2 = arrayMax(LN2) ;
      int mx = mx1 > mx2 ? mx1 : mx2 ;
      long int n, n1, n2 ;

      for (pos = 1 ; pos < mx ; pos++)
	{
	  n1 = pos < mx1 ? array (LN1, pos, long int) : 0 ;
	  n2 = pos < mx2 ? array (LN2, pos, long int) : 0 ;
	  n = n1 + n2 ;
	  if (n > 0)
	    aceOutf (seqco, "AliLn:%04d\t%s\t3\t%ld\t%ld\t%ld\n"
		     , pos
		     , ba->lane
		     , n, n1, n2
		     ) ;
	}
    }

  for (pos = 1 ; pos < arrayMax (errPos1) ; pos++)
    {
      long int np1 ;
      n1 = array (errPos1, pos, long int) ;
      if (n1 < 4 * nn1/arrayMax (errPos1)) continue ; 
      for (type = 1 ; type <= typeMax ; type++)
	{
	  np1 = array (errPosType1,  type +  (typeMax + 3) * pos, long int) ;
	  if (20 * np1 <= n1) continue  ;
	  aceOutf (ao, "RPOS\tf1\t%d\tany\t%d\t%s\t%ld\n", pos, n1, dictName (errDict, type), np1) ;
	}
      np1 = array (errPosType1,  typeMax+1 +  (typeMax + 3) * pos, long int) ;  
      if (20 * np1 >  n1)
	aceOutf (ao, "RPOS\tf1\t%d\tany\t%d\t%s\t%ld\n", pos, n1, "MultiInsertion",  np1) ;
      np1 = array (errPosType1,  typeMax+2 +  (typeMax + 3) * pos, long int) ;  
      if (20 * np1 >  n1)
	aceOutf (ao, "RPOS\tf1\t%d\tany\t%d\t%s\t%ld\n", pos, n1, "MultiDeletion", np1) ;
    }
  for (pos = 1 ; pos < arrayMax (errPos2) ; pos++)
    {
      long int n2, np2 ; 
      n2 = array (errPos2, pos, long int) ;
      if (n2 < 4 * nn2/arrayMax (errPos2)) continue ; 
      for (type = 1 ; type <= typeMax ; type++)
	{
	  np2 = array (errPosType2,  type +  (typeMax + 3) * pos, long int) ;
	  if (20 * np2 <= n2) continue  ;
	  aceOutf (ao, "RPOS\tf2\t%d\tany\t%d\t%s\t%ld\n", pos, n2, dictName (errDict, type), np2) ;
	}
      np2 = array (errPosType2,  typeMax+1 +  (typeMax + 3) * pos, long int) ;  
      if (20 * np2 >  n2)
	aceOutf (ao, "RPOS\tf2\t%d\tany\t%d\t%s\t%ld\n", pos, n2, "MultiInsertion", np2) ;
      np2 = array (errPosType2,  typeMax+2 +  (typeMax + 3) * pos, long int) ;  
      if (20 * np2 > n2)
	aceOutf (ao, "RPOS\tf2\t%d\tany\t%d\t%s\t%ld\n", pos, n2, "MultiDeletion", np2) ;
    }

  ac_free (h) ;
  fprintf (stderr, "Error profile done %s\n", timeShowNow ()) ;
  return nn ;
} /* baExportErrorProfile */

/*************************************************************************************/
/* Export the quality reports on tag alignments, per transcriptome/genome, unique/multiple
 *
 * We want to compute 5 tables
 *   Table 1: raw length of the mapped tags                   not available
 *   Table 2: length to be aligned (same minus pA, vector etc) == up2->ln
 *   Table 3: length actually aligned                          == up2->ali
 * 
 *   Table 4: percent: 100 * ln ali/ln to be ali            
 *   Table 5: percent: 100 * #errors/ln to be ali 
 *
 * Each table has 16 columns
 *   $1: length or percent, in case %, x means ]x-.5,x+.5], and is given from 0 to max
 *   $2,$3: sum of number of tags of the line, and if relevant same*ln=total bp of the line
 *   $4-7: hierarchic number of tags mapped to transcriptome unique/mult,then genome unique/mult
 *   $8-11: same but divided by col 2, i.e. a repartition by category Tu,Tm,Gu,Gm
 *   $12-16: same but divided by sum of column, i.e. a % histogram per length, Tu,Tm,Gu,Gm,Any
 *
 * Last line: a cumul computed exactly (not by addition of the lines)
 *
 *
 * For each lane we compute in  baExportAliProfile the raw numbers of tags (col 1, and 4 to 7)
 * Later in merge phase we add the sublibraries and groups and export them as ace files
 * Finally, we parse the ace files and produce the actual documents, tab delimited with legend
 */



static int baExportAliProfile (BA *ba) 
{  
  AC_HANDLE h = ac_new_handle () ;
  HIT *up ;
  HIT2 *up2 ;
  HIT3 *up3 ;
  int ii, n, u, tag = 0, mult, nn = 0, iTr, jTr, kTr, cumulE[8], cumulNN[8], oldTag[8] ;
  BOOL is1 = TRUE ;
  ACEOUT ao = aceOutCreate (ba->outFileName, ".aliProfile", ba->gzo, h) ;
  Array aa, aa1, aa2 ;
  Array tagLns = arrayHandleCreate (1024, KEYSET, h) ;
  Array alis = arrayHandleCreate (1024, KEYSET, h) ;
  Array aliPercents = arrayHandleCreate (1024, KEYSET, h) ;
  Array errPercents = arrayHandleCreate (1024, KEYSET, h) ;
  Array errCounts = arrayHandleCreate (1024, KEYSET, h) ;
  Array errCounts1 = arrayHandleCreate (1024, KEYSET, h) ;
  Array errCounts2 = arrayHandleCreate (1024, KEYSET, h) ;
  Array errNCounts = arrayHandleCreate (1024, KEYSET, h) ;
  Array nNPercents = arrayHandleCreate (1024, KEYSET, h) ;
  int tr1 = 0, tr2 = 0, genome = 0 ;
  int lnMax = 0, aliMax = 0, alipMax = 0, errpMax = 0, errcMax = 0, errcNMax = 0, nNMax = 0 ;
  int nPerfect1 = 0, nPerfect2 = 0, oldPerfectTag = 0 ;
  memset (cumulE, 0, sizeof(cumulE)) ;
  memset (cumulNN, 0, sizeof(cumulNN)) ;
  dictFind (ba->target_classDict, "DT_magic", &tr1) ;      /* first transciptome class */
  dictFind (ba->target_classDict, "PT_tRNA", &tr2) ;    /* last */
  dictFind (ba->target_classDict, "Z_genome", &genome) ; /* genome class */

  aceOutf (ao, "# %s : ali profile start", timeShowNow ()) ;
  for (ii = tag = 0, up = arrp (ba->hits, 0, HIT) ; ii < arrayMax (ba->hits) ; up++, ii++)
    {
      if (up->tag != tag)  /* imply stratified hierarchic mapping */
	memset (oldTag, 0, sizeof (oldTag)) ;
      is1 = 1 ;
      { 
	const char *ccp = dictName (ba->tagDict,up->tag) ;
	int i = strlen (ccp) ;
	ccp += i - 1 ;
	if (*ccp == '<')
	  is1 = 0 ;
      }
    
      up2 = arrp (ba->hits2, up->nn, HIT2) ;
      up3 = arrp (ba->hits3, up->nn, HIT3) ;
      tag = up->tag ;
      /*
	if (! up2->ln || up2->target_class < tr1) // reject non-transcriptome cases 
	continue ; 
      */

      if (up2->target_class < tr1) 
	jTr = 2 ;
      else if (up2->target_class <= tr2) 
	jTr = 4 ;
      else if (up2->target_class == genome)
	jTr = 6 ;
      else                  /* reject non-transcriptome/non-genome cases */
	continue ;

      nn++ ;
      mult = fastcMultiplicity (dictName (ba->tagDict,up->tag), 0, 0) ; 
      
      for (kTr = 0 ; kTr < 2 ; kTr++) /* kTr == 0 : compute the union, kTr == 1 compute each category */
	{
	  int ali = up2->ali ;
	  int nErr = up3->nErr ;
	  int nN = up3->nN ;
	  iTr = kTr ? jTr : 0 ; 
	  BOOL isFirst = TRUE ;

	  if (oldTag[iTr]) continue ; /* in union mode the tag should only be evaluated once */
	  if (1)
	    {
	      HIT *vp ; HIT2 *vp2 ; HIT3 *vp3 ;
	      int jj, k = 0 ;

	      for (jj = ii-1, vp = up-1 ; jj >= 0 && vp->tag == up->tag && vp->target == up->target ; vp--, jj--) 
		{
		  vp2 = arrayp (ba->hits2, vp->nn, HIT2) ;
		  vp3 = arrayp (ba->hits3, vp->nn, HIT3) ;
		 if ( up->chain == vp->chain  && vp2->target_class == up2->target_class)
		   {
		      k++ ; nErr += vp3->nErr ; nN += vp3->nN ;
		    }
		}
	      if (k)  /* this case is already counted */
		isFirst = FALSE ;
	      for (jj = ii+1, vp = up+1 ; jj < arrayMax (ba->hits) && vp->tag == up->tag && vp->target == up->target ; vp++, jj++) 
		{
		  vp2 = arrayp (ba->hits2, vp->nn, HIT2) ;
		  vp3 = arrayp (ba->hits3, vp->nn, HIT3) ;
		  if ( up->chain == vp->chain  && vp2->target_class == up2->target_class)
		    {
		      nErr += vp3->nErr ; nN += vp3->nN ;
		    }
		}
	    }
	  oldTag[iTr] = up->tag ;

	  if (isFirst)
	    {
	      u = up2->unicity ; 
	      if (u*u > 1 && ! (u == -2 && jTr == 6))           /* unique, non-unique, accepting the case -2 = single site two genes  */
		iTr++ ;
	      aa = array (tagLns,  iTr, KEYSET) ;
	      if (! aa)
		aa = array (tagLns,  iTr, KEYSET) = keySetHandleCreate (h) ;
	      keySet (aa, up2->ln) += mult ;
	      if (lnMax < up2->ln) lnMax = up2->ln ;
	      
	      aa = array (alis,  iTr, KEYSET) ;
	      if (! aa)
		aa = array (alis,  iTr, KEYSET) = keySetHandleCreate (h) ;
	      keySet (aa, ali) += mult ;
	      if (aliMax < ali) aliMax = ali ;
	      
	      aa = array (aliPercents, iTr, KEYSET) ;
	      n = 100 * (ali/(double)up2->ln) + .5 ; if(n > 100) n = 100 ;
	      if (! aa)
		aa = array (aliPercents,  iTr, KEYSET) = keySetHandleCreate (h) ;
	      keySet (aa, n) += mult ;
	      if (alipMax < n) alipMax = n ;

	      aa = array (errCounts,  iTr, KEYSET) ;
	      aa1 = array (errCounts1,  iTr, KEYSET) ;
	      aa2 = array (errCounts2,  iTr, KEYSET) ;
	      if (! aa)
		{
		  aa = array (errCounts,  iTr, KEYSET) = keySetHandleCreate (h) ;
		  aa1 = array (errCounts1,  iTr, KEYSET) = keySetHandleCreate (h) ;
		  aa2 = array (errCounts2,  iTr, KEYSET) = keySetHandleCreate (h) ;
		}

	      n = nErr ;
	      cumulE[iTr] += mult * nErr ;

	      keySet (aa, nErr) += mult ;
	      if (is1)
		keySet (aa1, nErr) += mult ;
	      else
		keySet (aa2, nErr) += mult ;
	      if (errcMax < nErr) errcMax = nErr ;
	      
	      if (iTr <= 1 && up2->ln == ali && nErr == 0 && up->tag != oldPerfectTag)
		{
		  if (is1) nPerfect1 += mult ;
		  else nPerfect2 += mult ;
		  oldPerfectTag = up->tag ;
		}
	      else if (0) printf ("%s ln=%d ali=%d nErr=%d\n"
				  , dictName (ba->tagDict,up->tag)
				  , up2->ln, ali, nErr
				  ) ;
	      aa = array (errPercents,  iTr, KEYSET) ;
	      n = (ali > 0 ? 100 * (nErr/(double)ali) + .5 : 0) ;
	      if (! aa)
		aa = array (errPercents,  iTr, KEYSET) = keySetHandleCreate (h) ;
	      keySet (aa, n) += mult ;
	      if (errpMax < n) errpMax = n ;
	  
	      aa = array (errNCounts,  iTr, KEYSET) ;
	      n = nN ;
	      cumulNN[iTr] += mult * nN ;
	      if (! aa)
		aa = array (errNCounts,  iTr, KEYSET) = keySetHandleCreate (h) ;
	      keySet (aa, nN) += mult ;
	      if (errcNMax < n) errcNMax = nN ;
	      
	      aa = array (nNPercents,  iTr, KEYSET) ;
	      n = 100 * (nN/(double)ali) + .5 ;
	      if (! aa)
		aa = array (nNPercents,  iTr, KEYSET) = keySetHandleCreate (h) ;
	      keySet (aa, n) += mult ;
	      if (nNMax < n) nNMax = n ;
	    }
	}
    }

  /* export */
  aceOutf (ao, "\n# Reads with best alignment to a single site, in transcripts or genome		Reads aligned in 2 to 9 sites, in transcripts or genome	Reads aligned in mitochondria, rRNA, ERCC at a single site	Reads aligned in mitochondria, rRNA, ERCC in 2 to 9 sites (rare)	Reads aligned to other transcripts, at a single site	Reads aligned to other transcripts, at 2 to 9 sites	Genomic reads aligned at a single site	Genomic reads aligned at 2 to 9 sites") ;
  for (ii = 0 ; ii <= lnMax ; ii++)
    {
      for (n = iTr = 0 ; iTr < 8 ; iTr++)
	{
	  aa = array (tagLns,  iTr, KEYSET) ;
	  n +=  (aa ? keySet (aa, ii) : 0) ;
	}
      if (! n) continue ; 
      aceOutf (ao, "\nClipped_length\t%d", ii) ;
      for (iTr = 0 ; iTr < 8 ; iTr++)
	{
	  aa = array (tagLns,  iTr, KEYSET) ;
	  aceOutf (ao, "\t%d", aa ? keySet (aa, ii) : 0) ;
	}
    }
  for (ii = 0 ; ii <= aliMax ; ii++)
    {
      for (n = iTr = 0 ; iTr < 8 ; iTr++)
	{
	  aa = array (alis,  iTr, KEYSET) ;
	  n +=  (aa ? keySet (aa, ii) : 0) ;
	}
      if (! n) continue ;
      aceOutf (ao, "\nAligned_length\t%d", ii) ;
      for (iTr = 0 ; iTr < 8 ; iTr++)
	{
	  aa = array (alis,  iTr, KEYSET) ;
	  aceOutf (ao, "\t%d", aa ? keySet (aa, ii) : 0) ;
	}
    }
  for (ii = 0 ; ii <= alipMax ; ii+=2)
    {
      for (n = iTr = 0 ; iTr < 8 ; iTr++)
	{
	  aa = array (aliPercents,  iTr, KEYSET) ;
	  n +=  (aa ? keySet (aa, ii) + keySet (aa, ii+1) : 0) ;
	}
      if (! n) continue ;
      aceOutf (ao, "\nPer_cent_aligned\t%d", ii < 100 ? ii + 1 : ii) ;
      for (iTr = 0 ; iTr < 8 ; iTr++)
	{
	  aa = array (aliPercents,  iTr, KEYSET) ;
	  aceOutf (ao, "\t%d",aa ?  keySet (aa, ii)  + keySet (aa, ii+1) : 0) ;
	}
    }
  for (ii = 0 ; ii <= errpMax ; ii++)
    {
      for (n = iTr = 0 ; iTr < 8 ; iTr++)
	{
	  aa = array (errPercents,  iTr, KEYSET) ;
	  n +=  (aa ? keySet (aa, ii) : 0) ;
	}
      if (! n) continue ;
      aceOutf (ao, "\nPer_cent_mismatch\t%d", ii) ;
      for (iTr = 0 ; iTr < 8 ; iTr++)
	{
	  aa = array (errPercents,  iTr, KEYSET) ;
	  aceOutf (ao, "\t%d", aa ? keySet (aa, ii) : 0) ;
	}
     }
  for (ii = 0 ; ii <= errcMax ; ii++)
    {
      for (n = iTr = 0 ; iTr < 8 ; iTr++)
	{
	  aa = array (errCounts,  iTr, KEYSET) ;
	  n +=  (aa ? keySet (aa, ii) : 0) ;
	}
      if (! n) continue ;
      aceOutf (ao, "\nCount_mismatch\t%d", ii) ;
      for (iTr = 0 ; iTr < 8 ; iTr++)
 	{
	  aa = array (errCounts,  iTr, KEYSET) ;
	  aceOutf (ao, "\t%d", aa ? keySet (aa, ii) : 0) ;
	}
      if (ba->seqco)
	{
	  int n1, n2 ;
	  aa = array (errCounts1,  0, KEYSET) ;
	  n1 = aa ?  keySet (aa, ii) : 0 ;
	  aa = array (errCounts1,  1, KEYSET) ;
	  if (aa) n1 += keySet (aa, ii) ;
	  aa = array (errCounts2,  0, KEYSET) ;
	  n2 = aa ?  keySet (aa, ii) : 0 ;
	  aa = array (errCounts2,  1, KEYSET) ;
	  if (aa) n2 += keySet (aa, ii) ;
	  aceOutf (ba->seqco, "\nMismatchPerRead:%03d\t%s\t3\t%d\t%d\t%d\n"
		   , ii, ba->lane, n1+n2, n1, n2
		   ) ;
	}
    }
  if (ba->seqco)
    aceOutf (ba->seqco, "Perfect_reads\t%s\t3\t%d\t%d\t%d\n"
		   , ba->lane, nPerfect1 + nPerfect2, nPerfect1, nPerfect2
		   ) ;
  aceOutf (ao, "\nPerfect_reads\t%s\t3\t%d\t%d\t%d"
	   , ba->lane, nPerfect1 + nPerfect2, nPerfect1, nPerfect2
	   ) ;
  for (ii = 0 ; ii <= errcNMax ; ii++)
    {
      for (n = iTr = 0 ; iTr < 8 ; iTr++)
	{
	  aa = array (errNCounts,  iTr, KEYSET) ;
	  n +=  (aa ? keySet (aa, ii) : 0) ;
	}
      if (! n) continue ;
      aceOutf (ao, "\nCount_ambiguous\t%d", ii) ;
      for (iTr = 0 ; iTr < 8 ; iTr++)
	{
	  aa = array (errNCounts,  iTr, KEYSET) ;
	  aceOutf (ao, "\t%d", aa ? keySet (aa, ii) : 0) ;
	}
     }
  for (ii = 0 ; ii <= nNMax ; ii++)
    {
      for (n = iTr = 0 ; iTr < 8 ; iTr++)
	{
	  aa = array (nNPercents,  iTr, KEYSET) ;
	  n +=  (aa ? keySet (aa, ii) : 0) ;
	}
      if (! n) continue ;
      aceOutf (ao, "\nPer_cent_Ambiguous\t%d", ii) ;
      for (iTr = 0 ; iTr < 8 ; iTr++)
	{
	  aa = array (nNPercents,  iTr, KEYSET) ;
	  aceOutf (ao, "\t%d", aa ? keySet (aa, ii) : 0) ;
	}
     }
  aceOutf (ao, "\nCumulated_mismatches\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d", cumulE[0], cumulE[1], cumulE[2], cumulE[3], cumulE[4], cumulE[5], cumulE[6], cumulE[7]) ;
  aceOutf (ao, "\nCumulated_ambiguous\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d", cumulNN[0], cumulNN[1], cumulNN[2], cumulNN[3], cumulNN[4], cumulNN[5], cumulNN[6], cumulNN[7]) ;
  aceOutf (ao, "\n# %s : ali profile done\n", timeShowNow ()) ;

  ac_free (h) ;
  fprintf (stderr, "Ali profile done %s\n", timeShowNow ()) ;
  return nn ;
} /* baExportAliProfile */

/************************************************************************************/
/* parse an ace file, Gene->Run_U to later ventilate proportionally the quasi unique reads */
static void baUniqueGeneSupportParse (BA *ba)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (ba->uniqueGeneSupportFileName, FALSE, h) ;
  const char *ccp ;
  Array aa ;
  int nLine = 0, gene = 0, nn = 0 ;
  BOOL isGene = FALSE ;

  if (!ai)
    messcrash ("Sorry, I cannot open the -unique_gene_support file : %s", ba->uniqueGeneSupportFileName) ;
  if (! ba->geneDict)
    ba->geneDict = dictHandleCreate (100000, h) ;

  aa = ba->uniqueGeneSupport = arrayHandleCreate (100000, float, ba->h) ;
  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (! (++nLine % 1000000))
	fprintf (stderr, "%d\n", nLine) ;
      if (! ccp || ! *ccp)
	{
	  isGene = FALSE ; gene = 0 ;
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Gene"))
	{
	  ccp = aceInWord (ai) ;
	  if (ccp && ! strncmp (ccp, "G_Any_", 6))
	    continue ;
	  dictAdd (ba->geneDict, ccp, &gene) ;
	  isGene = TRUE ;
	  continue ;
	}  
      if (ccp && isGene && !strcasecmp (ccp, "Run_U"))
	{
	  float f, *fp ;
	
	  ccp = aceInWord (ai) ; /* run name */
	  if (! ccp || 
	      (ba->run && strcmp (ccp, ba->run))
	      )
	    continue ;
	  if (! (ccp = aceInWord (ai))) /* index */
	    continue ;
	  if (! (ccp = aceInWord (ai))) /* seqs */
	    continue ;
	  if (! (ccp = aceInWord (ai))) /* "seqs" */
	    continue ;
	  aceInNext (ai) ;
	  if (aceInFloat (ai, &f))
	    {
	      fp = arrayp (aa, gene, float) ;
	      *fp += f ;
	      nn++ ;
	    }
	  continue ;
	}
    }

  fprintf (stderr, "baUniqueGeneSupportParse found counts for %d genes\n", nn) ;
  ac_free (h) ;
  return ;
} /* baUniqueGeneSupportParse */

/*************************************************************************************/
/* export the tags and aligned bp per gene */
#define MRNA_WIGGLE_STEP 10
static int baExportGeneSupport (BA *ba, BOOL unique) 
{  
  AC_HANDLE h = ac_new_handle () ;
  HIT *up, *vp, *wp ;
  HIT2 *up2, *vp2, *wp2 ;
  HIT3 *up3, *vp3 ;
  int ii, jj, kk, gene, fromGene = 0, run, score = 0, tag = 0, lastTag = 0, ali, mult, nn = 0, chloroClass = 0, mitoClass = 0, spikeInClass = 0 ;
  double uu, uuu, uuCumul ;
  BOOL isFirstMrna = FALSE ;
  ACEOUT ao = aceOutCreate (ba->outFileName, messprintf(".%sSupport.%s", ba->mrnaSupport ? "mrna" : "gene", unique ? "u" : "nu"), ba->gzo, h) ;	
  ACEOUT wao8 = 0 ;
  ACEOUT wao5 = 0 ;
  int wKeptTranscripts = 0 ;
  int MX = 256 ;      /* target_class tables */
  int nTi[MX] ;
  KEYSET wiggle, wiggle3p8kb = 0, wiggle3p5kb = 0 ;
  Array aa = ba->uniqueGeneSupport ;
  Array gsF = arrayHandleCreate (100000, double,h) ;
  Array gsR = arrayHandleCreate (100000, double,h) ;
  Array aliF = arrayHandleCreate (100000, double,h) ;
  Array aliR = arrayHandleCreate (100000, double,h) ;
  Array tagF = arrayHandleCreate (100000, double,h) ;
  Array tagR = arrayHandleCreate (100000, double,h) ;

  Array nReadsF = arrayHandleCreate (100000, double,h) ;
  Array nReadsR = arrayHandleCreate (100000, double,h) ;
  Array nReadsOkF = arrayHandleCreate (100000, double,h) ;
  Array nReadsOkR = arrayHandleCreate (100000, double,h) ;

  Array nerrF = arrayHandleCreate (100000, double,h) ;
  Array nerrR = arrayHandleCreate (100000, double,h) ;
  Array nA2gF = arrayHandleCreate (100000, double,h) ;
  Array nA2gR = arrayHandleCreate (100000, double,h) ;
  Array nPartialF = arrayHandleCreate (100000, double,h) ;
  Array nPartialR = arrayHandleCreate (100000, double,h) ;
  Array nOrphanF = arrayHandleCreate (100000, double,h) ;
  Array nOrphanR = arrayHandleCreate (100000, double,h) ;
  Array nBadtopoF = arrayHandleCreate (100000, double,h) ;
  Array nBadtopoR = arrayHandleCreate (100000, double,h) ;

  Array nMultiF = arrayHandleCreate (100000, double,h) ;
  Array nMultiR = arrayHandleCreate (100000, double,h) ;
  Array nMulti2F = arrayHandleCreate (100000, double,h) ;
  Array nMulti2R = arrayHandleCreate (100000, double,h) ;

  Array wigglesF = arrayHandleCreate (100000, KEYSET,h) ;
  Array wigglesR = arrayHandleCreate (100000, KEYSET,h) ;

  dictAdd (ba->target_classDict,  "0_SpikeIn", &spikeInClass) ;
  dictAdd (ba->target_classDict,  "A_mito", &mitoClass) ;
  dictAdd (ba->target_classDict,  "C_chloro", &chloroClass) ;

  if (! ba->run)
    messcrash ("-geneSupport requires option run run_name") ;
  if (! ba->target_class)
    messcrash ("-geneSupport requires option -target_class to correctly select the genes") ;

  dictAdd (ba->runDict, ba->run, &(run)) ;

  if (ba->selected8kbDict)
    { 
      int i ;
      wao8 = aceOutCreate (ba->outFileName, ".3pHisto.8kb.txt", FALSE, h) ;
      wiggle3p8kb = keySetHandleCreate (h) ;

      aceOutDate (wao8, "###", ba->run) ;
      aceOutf (wao8, "## Number of contributing transcripts\t%d\tCumul\t%ld\n", dictMax (ba->selected8kbDict)) ;
      aceOutf (wao8,
	       "## The transcripts are selected such that the maximum of their coverage plot, usually representing the major polyA addition site, is further than 8kb from the 5' end of the annotated transcript.\n"
	       "## This implies that the histogram below 8kb represents the 3' biais, and above 8kb a commbination of the 3'bias and the prevalence of very long transcripts\n"
	       "## The coverage plots were piled up, aligning the maxima at position zero. The x coordinates run 3' to 5' and represent the distance to the 3' end\n"
	       "# Run\tDistance from 3' end"
	       ) ;
      for (i = 0 ; i <= 16000 ; i += MRNA_WIGGLE_STEP)
	aceOutf (wao8, "\t%d", i) ;
    }
  if (ba->selected5kbDict)
    { 
      int i ;
      wao5 = aceOutCreate (ba->outFileName, ".3pHisto.5kb.txt", FALSE, h) ;
      wiggle3p5kb = keySetHandleCreate (h) ;

      aceOutDate (wao5, "###", ba->run) ;
      aceOutf (wao5, "## Number of contributing transcripts\t%d\tCumul\t%ld\n", dictMax (ba->selected5kbDict)) ;
      aceOutf (wao5,
	       "## The transcripts are selected such that the maximum of their coverage plot, usually representing the major polyA addition site, is further than 8kb from the 5' end of the annotated transcript.\n"
	       "## This implies that the histogram below 8kb represents the 3' biais, and above 8kb a commbination of the 3'bias and the prevalence of very long transcripts\n"
	       "## The coverage plots were piled up, aligning the maxima at position zero. The x coordinates run 3' to 5' and represent the distance to the 3' end\n"
	       "# Run\tDistance from 3' end"
	       ) ;
      for (i = 0 ; i <= 8000 ; i += MRNA_WIGGLE_STEP)
	aceOutf (wao5, "\t%d", i) ;
    }


  /* compute on a per read basis the quality counters for each gene */
  if (ba->pair) /* otherwise ->clone is not initialized */
    {
      for (ii = 0, up = arrp (ba->hits, 0, HIT) ; ii < arrayMax (ba->hits) ; up++, ii++)
	{
	  int ok = 0, ok1 ;
	  int tag = up->tag ;
	  int clone = up->clone ;

	  up2 = arrp (ba->hits2, up->nn, HIT2) ;
	  if (! ba->target_class || up2->target_class == ba->target_class)
	    for (jj = ii, vp = up ; ok != 3 && jj < arrayMax (ba->hits) && vp->clone == clone ; vp++, jj++)
	      {
		vp2 = arrp (ba->hits2, vp->nn, HIT2) ;
		if (! ba->target_class || vp2->target_class == ba->target_class)
		  {
		    ok1 = (tag == vp->tag) ? 0x1 : 0x2 ;
		    if (ok ^ ok1)
		      {
			vp2 = arrp (ba->hits2, vp->nn, HIT2) ;
			if (! vp2->badPair)
			  {
			    ok |= ok1 ;
			    if (100 * vp2->ali < 80 * vp2->ln && vp2->ali < vp2->ln - 8) 
			      vp2->badPair = 2 ; /* correct topology but not sufficiently long alignments */
			  }
		      }
		  }
	      }
	  
	  if (ok != 0x3)
	    {
	      for (jj = ii, vp = up ; jj < arrayMax (ba->hits) && vp->clone == clone ; vp++, jj++)
		{
		  vp2 = arrp (ba->hits2, vp->nn, HIT2) ;
		  vp2->badPair = 1 ; /* orphan */
		}
	    }
	  for (jj = ii, vp = up ; jj < arrayMax (ba->hits) && vp->clone == clone ; vp++, jj++) ;
	  up = vp - 1 ; ii = jj - 1 ;
	}
    }

  for (ii = score = lastTag = 0, up = arrp (ba->hits, 0, HIT) ; ii < arrayMax (ba->hits) ; up++, ii++)
    {
      int isFirstFragment = 1 ;
	{ 
	  const char *ccp = dictName (ba->tagDict,up->tag) ;
	  int i = strlen (ccp) ;
	  ccp += i - 1 ;
	  if (*ccp == '<')
	    isFirstFragment = -1 ;
	}
      up2 = arrp (ba->hits2, up->nn, HIT2) ;
      up3 = arrp (ba->hits3, up->nn, HIT3) ;

      /* measure the quality metrics even if the pair is not correctly aligned */
      if (ba->mrnaSupport || up2->target_class == spikeInClass || up2->target_class == mitoClass || up2->target_class == chloroClass)
	gene = up->target ;
      else
	gene = up->gene ;

      uu = up2->unicity ;
      if (uu < 0) uu = -uu ;

      mult = up2->mult ;
      if (up->score && tag != up->tag)
	{
	  tag = up->tag ; /* so we export only once per tag */
	  isFirstFragment = 1 ;  fromGene = 0 ;
	  { 
	    const char *ccp = dictName (ba->tagDict,up->tag) ;
	    int i = strlen (ccp) ;
	    ccp += i - 1 ;
	    if (*ccp == '<')
	      isFirstFragment = -1 ;
	  }
	  gene = 0 ; 	
	  uuCumul = 0 ;
	  for (vp = up, jj = ii ;  vp->tag == tag && jj < arrayMax (ba->hits) ; vp++, jj++)
	    {
	      vp2 = arrp (ba->hits2, vp->nn, HIT2) ;
	      vp3 = arrp (ba->hits3, vp->nn, HIT3) ;
	      if (ba->target_class && vp2->target_class != ba->target_class)
		continue ;

	      if (ba->mrnaSupport || vp2->target_class == spikeInClass || vp2->target_class == mitoClass || vp2->target_class == chloroClass)
		{
		  if (! vp->target || vp->target == gene)
		    continue ;
		  if (fromGene != vp->gene)
		    isFirstMrna = TRUE ;
		  else
		    isFirstMrna = FALSE ;
		  fromGene = vp->gene ;
		  gene = vp->target ;
		  if (ba->mrnaSupport && ba->hierarchic && ! isFirstMrna)
		    continue ;
		}
	      else
		{
		  if (! vp->gene || vp->gene == gene)
		    continue ;
		  gene = vp->gene ;
		}
	      memset (nTi, 0, sizeof (nTi)) ;
	      /* mult = fastcMultiplicity (dictName (ba->tagDict,vp->tag), nTi, MX) ;  */
	      mult = vp2->mult ;
	      ali = vp2->ali ;

	      uuu = vp2->unicity ;
	      if (uuu < 0) uuu = -uuu ;
	      if (uuu == 0) uuu = 1 ;
	      uu = 1 ;


	      if (uuu > 1) /* find the correct proportions */
		{
		  int gene2 = 0 ;
		  int trueGene = vp->gene ;

		  if (uuCumul == 0) /* count the read cumul of all genes hit by this tag */
		    {
		      if (aa)
			{
			  for (wp = up, kk = ii ;  wp->tag == tag && kk < arrayMax (ba->hits) ; wp++, kk++)    
			    {
			      if (! wp->gene || wp->gene == gene2)
				continue ;
			      wp2 = arrp (ba->hits2, wp->nn, HIT2) ;
			      if (wp2->target_class != vp2->target_class)
				continue ;
			      gene2 = wp->gene ;
			      uuCumul += 20 + array (aa, gene2, float) ; /* damp the proportion in small numbers */
			    }
			}
		    }
		     
		  if (uuCumul == 0)
		    uu = uuu ;
		  else
		    {
		      uu = uuCumul/(20 + array (aa, trueGene, float)) ;
		      array (aa, trueGene, float) += mult/uu ; /* accumulate, so if we add lots of quasi unique, they become more even */
		    }
		}

	      if (! ba->stranded || ba->stranded * isFirstFragment * (vp2->a2 - vp2->a1) > 0)
		{
		  int k ;
		  int a1 = vp2->a2 > vp2->a1 ? vp2->a1 : vp2->a2 ;
		  int a2 = vp2->a2 > vp2->a1 ? vp2->a2 : vp2->a1 ;

		  array (nReadsF, gene, double) += mult/uu ;
		  
		  if (vp2->dPair == -14)                      /* synchronize to hack */
		    array (nBadtopoF, gene, double) += (vp2->mult/uu) ;
		  else if (ba->pair && vp2->dPair && vp2->badPair == 1)                        /* synchronize to hack */
		    array (nOrphanF, gene, double) += (vp2->mult/uu) ;
		  else if (vp2->badPair == 2)
		    array (nPartialF, gene, double) += (vp2->mult/uu) ;
		  else
		    array (nReadsOkF, gene, double) += mult/uu ;
		  
		  if (vp2->unicity == -2)
		    array (nMulti2F, gene, double) += mult/uu ;
		  else if (uu > 1)
		    array (nMultiF, gene, double) += mult/uu ;
		  
		  if (unique && vp2->unicity != 1)
		    continue ;
		  if (100 * vp2->ali < 80 * vp2->ln)
		    continue ;   
		  if (ba->maxErr && vp3->nErr > ba->maxErr && 100 * vp3->nErr > ba->maxErrRate * vp2->ali)
		    continue ;

		  if (ba->pair && vp2->badPair)  /* non correct topology or not sufficiently long alignments */
		    continue ;
   
		  array (nerrF, gene, double) += vp3->nErr*mult/uu ;
		  if (vp3->errTarget)
		    {
		      const char *cs, *cr = dictName(ba->errDict, vp3->errTarget) ;
		      if (cr &&
			  (cs = strstr (cr, "a>g")) &&
			  (cs = strstr (cs+1, "a>g"))
			  )
			array (nA2gF, gene, double) += mult/uu ;
		    }
		       
		  array (gsF, gene, double) += 1/uu ;
		  array (tagF, gene, double) += mult/uu ;
		  array (aliF, gene, double) += ali*mult/uu ; 
		  if (ba->mrnaSupport || vp2->target_class == spikeInClass || vp2->target_class == mitoClass || vp2->target_class == chloroClass)
		    {
		      wiggle = array (wigglesF, gene, KEYSET) ;
		      if (! wiggle)
			wiggle = array (wigglesF, gene, KEYSET) = keySetHandleCreate (h) ;
		      for (k = a1 ; k < a2 ; k += MRNA_WIGGLE_STEP)
			keySet (wiggle, k/MRNA_WIGGLE_STEP) += 100 * mult/uu ;
		    }
		}
	      else
		{
		  int k ;
		  int a1 = vp2->a2 > vp2->a1 ? vp2->a1 : vp2->a2 ;
		  int a2 = vp2->a2 > vp2->a1 ? vp2->a2 : vp2->a1 ;

		  array (nReadsR, gene, double) += mult/uu ;
		  
		  if (vp2->dPair == -14)                      /* synchronize to hack */
		    array (nBadtopoR, gene, double) += (vp2->mult/uu) ;
		  else if (ba->pair && vp2->dPair && vp2->badPair == 1)                        /* synchronize to hack */
		    array (nOrphanR, gene, double) += (vp2->mult/uu) ;
		  else if (vp2->badPair == 2)
		    array (nPartialR, gene, double) += (vp2->mult/uu) ;
		  else
		    array (nReadsOkR, gene, double) += mult/uu ;
		  
		  if (vp2->unicity == -2)
		    array (nMulti2R, gene, double) += mult/uu ;
		  else if (uu > 1)
		    array (nMultiR, gene, double) += mult/uu ;
		  
		       
		  if (unique && vp2->unicity != 1)
		    continue ;
		  if (100 * vp2->ali < 80 * vp2->ln)
		    continue ;
		  if (ba->maxErr && vp3->nErr > ba->maxErr && 100 * vp3->nErr > ba->maxErrRate * vp2->ali)
		    continue ;

		  if (ba->pair && vp2->badPair)  /* non correct topology or not sufficiently long alignments */
		    continue ;
   
		  array (nerrR, gene, double) += up3->nErr*mult/uu ;
		  if (vp3->errTarget)
		    {
		      const char *cs, *cr = dictName(ba->errDict, vp3->errTarget) ;
		      if (cr &&
			  (cs = strstr (cr, "a>g")) &&
			  (cs = strstr (cs+1, "a>g"))
			  )
			array (nA2gR, gene, double) += mult/uu ;
		    }
		  array (gsF, gene, double) += 0 ;
		  array (gsR, gene, double) += 1/uu ;
		  array (tagR, gene, double) += mult/uu ;
		  array (aliR, gene, double) += ali*mult/uu ;
		  if (ba->mrnaSupport || vp2->target_class == spikeInClass || vp2->target_class == mitoClass || vp2->target_class == chloroClass)
		    {
		      wiggle = array (wigglesR, gene, KEYSET) ;
		      if (! wiggle)
			wiggle = array (wigglesR, gene, KEYSET) = keySetHandleCreate (h) ;
		      for (k = a1 ; k < a2 ; k += MRNA_WIGGLE_STEP)
			keySet (wiggle, k/MRNA_WIGGLE_STEP) += 100 * mult/uu ;
		    }
		}
	    }
	  nn ++ ;
	}   
    }

  /* export */
  for (ii = 0 ; ii < arrayMax (gsF) || ii < arrayMax (nReadsF) || ii < arrayMax (nReadsR) ; ii++)
    {
      if (array (gsF, ii, double) > 0 || array (nPartialF, ii, double) > 0 || array (nOrphanF, ii, double) > 0)
      {
	aceOutf (ao, "%s", dictName (ba->target_classDict, ba->target_class)) ;
	aceOutf (ao, "\t%s", dictName (ba->targetDict, ii)) ;
	aceOutf (ao, "\t%s", unique ? "u" : "nu") ;
	aceOutf (ao, "\t%s", ba->run) ;
	aceOutf (ao, "\t%.1f\t%.2f\t%.0f", array (gsF, ii, double), array (tagF, ii, double), array (aliF, ii, double)) ;

	aceOutf (ao, "\t%.1f", array (nReadsF, ii, double)) ;
	aceOutf (ao, "\t%.1f", array (nReadsOkF, ii, double)) ;

	aceOutf (ao, "\t%.1f", array (nerrF, ii, double)) ;
	aceOutf (ao, "\t%.1f", array (nA2gF, ii, double)) ;
	aceOutf (ao, "\t%.1f", array (nPartialF, ii, double)) ;
	aceOutf (ao, "\t%.1f", array (nOrphanF, ii, double)) ;
	aceOutf (ao, "\t%.1f", array (nBadtopoF, ii, double)) ;
	aceOutf (ao, "\t%.1f", array (nMultiF, ii, double)) ;
	aceOutf (ao, "\t%.1f", array (nMulti2F, ii, double)) ;
	/* example
                                seq     tags      kb    F       Fok     err     a>g     Part   orphan  bad-topo  Multi  Multi2
ET_av	mari	u	xx	2.0	2.00	199	3.5	2.0	4.0	0.0	0.0	0.0	0.0	0.5	0.0
ET_av	mari	u	xx	12.0	16.00	1539	16.0	16.0	2.0	0.0	0.0	0.0	0.0	0.0	0.0
ET_av	mari	u	xx	7.0	11.00	1078	49.2	37.8	11.1	0.0	0.3	0.0	0.0	37.2	0.0
ET_av	mari	u	xx	64.0	86.00	8475	97.3	86.0	24.3	0.0	0.3	0.0	0.0	7.3	0.0
ET_av	POF1B	u	Rhs5347	14.0	16.00	1607	30.0	0.0	18.0	0.0	0.0	0.0	19.0	0.0	0.0
ET_av	SMS	u	Rhs5347	443.0	579.00	58000	3250.0	5.0	1371.0	5.0	64.0	0.0	3220.0	0.0	0.0
ET_av	PGRMC1	u	Rhs5347	206.0	268.00	26817	592.0	0.0	182.0	0.0	22.0	0.0	590.0	0.0	0.0
	*/
	wiggle = array (wigglesF, ii, KEYSET) ;
	if (wiggle && keySetMax (wiggle)) 
	  {
	    int k, kMax = keySetMax (wiggle) ;
	    KEY *kp = arrp (wiggle, 0, KEY) ;

	    aceOutf (ao, "\tWIGGLE_START") ;
	    for (k = 0 ; k < kMax && *kp == 0 ; k++, kp++) ;
	    aceOutf (ao, "\t%d", MRNA_WIGGLE_STEP * k) ;
	    for ( ; k < kMax ; k++, kp++)
	      if (*kp == 100 * ((*kp) / 100))
		aceOutf (ao, "\t%d", *kp/100) ;
	      else
		aceOutf (ao, "\t%.1f", *kp/100.0) ;
	  }
	aceOutf (ao, "\n") ;
	if (wiggle &&
	    ba->selected8kbDict && 
	    dictFind (ba->selected8kbDict, dictName (ba->targetDict, ii), 0)
	    )
	  {
	    int kMax = keySetMax (wiggle) ;
	    int i, nn, wMaxPos ;
	    KEY wMaxValue, *kp, wCumul = 0 ;
	    BOOL wDebug = TRUE ;
	    
	    /* find the max of the histo under the condition that we eat at most 10% of the total */
	    wMaxPos = kMax - 1 ; wMaxValue = 0 ;
	    for (wCumul = 0, i =  wMaxPos, kp = arrp (wiggle, i, KEY) ; i >= 0 ; kp--, i--)
	      wCumul += *kp ;
	    for (nn = 0, i =  wMaxPos, kp = arrp (wiggle, i, KEY) ; i >= 0 ; kp--, i--)
	      {
		nn += *kp ;
		if (10 * nn < wCumul && *kp >  wMaxValue)
		  { wMaxPos = i ; wMaxValue = *kp ;  }
	      }
	    if (wMaxPos >= 8000/MRNA_WIGGLE_STEP)
	      {
		wKeptTranscripts++ ;
		   if (wDebug) /* debugging */
		     {
		       aceOutf (wao8, "\n%s\t%s", ba->run, dictName (ba->targetDict, ii)) ;
		       for (i = 0 ; i <= wMaxPos && MRNA_WIGGLE_STEP * i <= 16000 ; i++)
			 aceOutf (wao8, "\t%d",  keySet (wiggle, wMaxPos - i)) ;
		     }
		   
		   /* register the shifted smoothed histo */
		   for (i = 0 ;  i <= wMaxPos && MRNA_WIGGLE_STEP * i <= 16000 ; i++)
		     keySet (wiggle3p8kb, i) += keySet (wiggle, wMaxPos - i) ; /* count backwards, starting at the max, which is probably close to the 3' end of the transcript */
		   
	      }
	  }
	if (wiggle &&
	    ba->selected5kbDict && 
	    dictFind (ba->selected5kbDict, dictName (ba->targetDict, ii), 0)
	    )
	  {
	    int kMax = keySetMax (wiggle) ;
	    int i, nn, wMaxPos ;
	    KEY wMaxValue, *kp, wCumul = 0 ;
	    BOOL wDebug = TRUE ;
	    
	    /* find the max of the histo under the condition that we eat at most 10% of the total */
	    wMaxPos = kMax - 1 ; wMaxValue = 0 ;
	    for (wCumul = 0, i =  wMaxPos, kp = arrp (wiggle, i, KEY) ; i >= 0 ; kp--, i--)
	      wCumul += *kp ;
	    for (nn = 0, i =  wMaxPos, kp = arrp (wiggle, i, KEY) ; i >= 0 ; kp--, i--)
	      {
		nn += *kp ;
		if (10 * nn < wCumul && *kp >  wMaxValue)
		  { wMaxPos = i ; wMaxValue = *kp ;  }
	      }
	    if (wMaxPos >= 5000/MRNA_WIGGLE_STEP)
	      {
		wKeptTranscripts++ ;
		   if (wDebug) /* debugging */
		     {
		       aceOutf (wao5, "\n%s\t%s", ba->run, dictName (ba->targetDict, ii)) ;
		       for (i = 0 ; i <= wMaxPos && MRNA_WIGGLE_STEP * i <= 16000 ; i++)
			 aceOutf (wao5, "\t%d",  keySet (wiggle, wMaxPos - i)) ;
		     }
		   
		   /* register the shifted smoothed histo */
		   for (i = 0 ;  i <= wMaxPos && MRNA_WIGGLE_STEP * i <= 16000 ; i++)
		     keySet (wiggle3p5kb, i) += keySet (wiggle, wMaxPos - i) ; /* count backwards, starting at the max, which is probably close to the 3' end of the transcript */
		   
	      }
	  }
      }
      
      if (array (gsR, ii, double) > 0 || array (nPartialR, ii, double) > 0 || array (nOrphanR, ii, double) > 0)
      {
	aceOutf (ao, "%s", dictName (ba->target_classDict, ba->target_class)) ;
	aceOutf (ao, "\t%s", dictName (ba->targetDict, ii)) ;
	aceOutf (ao, "\t%s", unique ? "u" : "nu") ;
	aceOutf (ao, "\tAnti_%s", ba->run) ;
	aceOutf (ao, "\t%.1f\t%.2f\t%.0f", array (gsR, ii, double), array (tagR, ii, double), array (aliR, ii, double)) ;

	aceOutf (ao, "\t%.1f", array (nReadsR, ii, double)) ;
	aceOutf (ao, "\t%.1f", array (nReadsOkR, ii, double)) ;

	aceOutf (ao, "\t%.1f", array (nerrR, ii, double)) ;
	aceOutf (ao, "\t%.1f", array (nA2gR, ii, double)) ;
	aceOutf (ao, "\t%.1f", array (nPartialR, ii, double)) ;
	aceOutf (ao, "\t%.1f", array (nOrphanR, ii, double)) ;
	aceOutf (ao, "\t%.1f", array (nBadtopoR, ii, double)) ;
	aceOutf (ao, "\t%.1f", array (nMultiR, ii, double)) ;
	aceOutf (ao, "\t%.1f", array (nMulti2R, ii, double)) ;

	wiggle = array (wigglesR, ii, KEYSET) ;
	if (wiggle && keySetMax (wiggle)) 
	  {
	    int k, kMax = keySetMax (wiggle) ;
	    KEY *kp = arrp (wiggle, 0, KEY) ;

	    aceOutf (ao, "\tWIGGLE_START") ;
	    for (k = 0 ; k < kMax && *kp == 0 ; k++, kp++) ;
	    aceOutf (ao, "\t%d", MRNA_WIGGLE_STEP * k) ;
	    for ( ; k < kMax ; k++, kp++)
	      if (*kp == 100 * ((*kp) / 100))
		aceOutf (ao, "\t%d", *kp/100) ;
	      else
		aceOutf (ao, "\t%.1f", *kp/100.0) ;
	  }
	aceOutf (ao, "\n") ;
      }
    } 
 
  /* export the global histo */
  if (wiggle3p8kb)
    {
      long int wCumul = 0 ;
      int i ;

      for (i = 0 ; i < keySetMax (wiggle3p8kb) && MRNA_WIGGLE_STEP *  i < 16000 ; i++)
	wCumul += keySet (wiggle3p8kb, i) ;
      aceOutf (wao8, "\n%s\tKept %d/%d transcripts", ba->run, wKeptTranscripts, dictMax (ba->selected8kbDict)) ;
      for (i = 0 ; i < keySetMax (wiggle3p8kb) && MRNA_WIGGLE_STEP *  i < 16000 ; i++)
	aceOutf (wao8, "\t%d", keySet (wiggle3p8kb, i)) ;

      aceOutf (wao8, "\n") ;
    }
  if (wiggle3p5kb)
    {
      long int wCumul = 0 ;
      int i ;

      for (i = 0 ; i < keySetMax (wiggle3p5kb) && MRNA_WIGGLE_STEP *  i < 16000 ; i++)
	wCumul += keySet (wiggle3p5kb, i) ;
      aceOutf (wao5, "\n%s\tKept %d/%d transcripts", ba->run, wKeptTranscripts, dictMax (ba->selected5kbDict)) ;
      for (i = 0 ; i < keySetMax (wiggle3p5kb) && MRNA_WIGGLE_STEP *  i < 16000 ; i++)
	aceOutf (wao5, "\t%d", keySet (wiggle3p5kb, i)) ;

      aceOutf (wao5, "\n") ;
    }


  ac_free (h) ;
  return nn ;
} /* baExportGeneSupport */

/*************************************************************************************/

static void baExportIntergenicSupportRegister (BA *ba, HIT *vp, int isFirstFragment, Array genes, Array exons, int gene, int exon
					       , Array gsF, Array gsR, Array gsFg, Array gsRg
					       , Array aliF, Array aliR, Array aliFg, Array aliRg
					       , Array tagF, Array tagR, Array tagFg, Array tagRg
					       , Array errF, Array errR, Array errFg, Array errRg
					       )
{
  int mult, ali ;
  double uu ;
  HIT2 *vp2 = arrp (ba->hits2, vp->nn, HIT2) ;

  mult = fastcMultiplicity (dictName (ba->tagDict,vp->tag), 0, 0) ;
  uu = vp2->unicity ;
  ali = vp2->ali ;
  if (uu < 0) uu = -uu ;
  if (uu == 0) uu = 1 ;
  if (! ba->stranded || ba->stranded * isFirstFragment * (vp2->a2 - vp2->a1) > 0)
    {
      if (exon && keySetInsert (exons, exon)) /* so we export only once per tag */
	{
	  array (gsF, exon, double) += 1/uu ;
	  array (tagF, exon, double) += mult/uu ;
	  array (aliF, exon, double) += ali*mult/uu ;
	}
      if (gene && keySetInsert (genes, gene))
	{
	  array (gsFg, gene, double) += 1/uu ;
	  array (tagFg, gene, double) += mult/uu ;
	  array (aliFg, gene, double) += ali*mult/uu ;
	}
    }
  else
    {
      if (exon && keySetInsert (exons, exon)) /* so we export only once per tag */
	{
	  array (gsF, exon, double) += 0 ;
	  array (gsR, exon, double) += 1/uu ;
	  array (tagR, exon, double) += mult/uu ;
	  array (aliR, exon, double) += ali*mult/uu ;
	}
      if (gene && keySetInsert (genes, gene))
	{
	  array (gsFg, gene, double) += 0 ;
	  array (gsRg, gene, double) += 1/uu ;
	  array (tagRg, gene, double) += mult/uu ;
	  array (aliRg, gene, double) += ali*mult/uu ;
	}
    }
} /* baExportIntergenicSupportRegister */

/*************************************************************************************/
/* foreach tag first best alignment
 * if in KT_RefSeq or ET_av or Z_intron
 * check if we overlap an exon (not caring for local errors since we clip) by 8bp using the remap coords
 * we count sequences, tags and bp (adding the aligned bp, counting each tag touching by at least 8bp CF baHit2Exon)
 */
static int baExportIntergenicSupport (BA *ba) 
{  
  AC_HANDLE h = ac_new_handle () ;
  HIT *up, *vp ;
  HIT2 *up2, *vp2 ;
  int ii, jj, exon, gene, run, score = 0, tag = 0, nn = 0, isDown ;
  int isFirstFragment ;
  int unique = ba->unique ;
  int Z_genome = 0 ;
  KEYSET exons = keySetHandleCreate (h) ;
  KEYSET genes = keySetHandleCreate (h) ;

  ACEOUT ao = aceOutCreate (ba->outFileName, unique ? ".intergenicSupport.u" : ".intergenicSupport.nu", ba->gzo, h) ;

  Array gsF = arrayHandleCreate (100000, double,h) ;
  Array gsR = arrayHandleCreate (100000, double,h) ;

  Array aliF = arrayHandleCreate (100000, double,h) ;
  Array aliR = arrayHandleCreate (100000, double,h) ;

  Array tagF = arrayHandleCreate (100000, double,h) ;
  Array tagR = arrayHandleCreate (100000, double,h) ;

  Array errF = arrayHandleCreate (100000, double,h) ;
  Array errR = arrayHandleCreate (100000, double,h) ;
  
  Array gsFg = arrayHandleCreate (100000, double,h) ;
  Array gsRg = arrayHandleCreate (100000, double,h) ;

  Array aliFg = arrayHandleCreate (100000, double,h) ;
  Array aliRg = arrayHandleCreate (100000, double,h) ;

  Array tagFg = arrayHandleCreate (100000, double,h) ;
  Array tagRg = arrayHandleCreate (100000, double,h) ;
  
  Array errFg = arrayHandleCreate (100000, double,h) ;
  Array errRg = arrayHandleCreate (100000, double,h) ;
  
  if (! ba->run)
    messcrash ("-intergenicSupport requires option run run_name") ;

  if (ba->target_classDict) 
    dictFind (ba->target_classDict, "Z_genome", &Z_genome) ;
  dictAdd (ba->runDict, ba->run, &(run)) ;

  for (ii = score = tag = 0, up = arrp (ba->hits, 0, HIT) ; ii < arrayMax (ba->hits) ; up++, ii++)
    {
      up2 = arrp (ba->hits2, up->nn, HIT2) ;
      if (up->score && tag != up->tag)
	{
	  BOOL foundTranscript = FALSE ;
	  tag = up->tag ; keySetMax (exons) = 0 ;  keySetMax (genes) = 0 ;
	  isFirstFragment = 1 ;
	  
	  { 
	    const char *ccp = dictName (ba->tagDict,up->tag) ;
	    int i = strlen (ccp) ;
	    ccp += i - 1 ;
	    if (*ccp == '<')
	      isFirstFragment = -1 ;
	  }
	  for (vp = up, jj = ii ;  vp->tag == tag && jj < arrayMax (ba->hits) ; vp++, jj++)
	    {
	      vp2 = arrp (ba->hits2, vp->nn, HIT2) ;
	      
	      if (unique && vp2->unicity != 1)
		continue ;
	      if (vp2->target_class == Z_genome)
		{
		  if (!foundTranscript)
		    {
		      int a1 = vp2->a1, a2 = vp2->a2 ;
		      exon = 0 ; gene = 0 ;
		      while ((isDown = baHit2Gene (ba, vp, &gene))) /* the first time we report the whole gene */
			{
			  baExportIntergenicSupportRegister (ba, vp, isDown * isFirstFragment, genes, exons, gene, exon
							     , gsF, gsR, gsFg, gsRg
							     , aliF, aliR, aliFg, aliRg
							     , tagF, tagR, tagFg, tagRg	
							     , errF, errR, errFg, errRg	
							     ) ;
			  nn++ ;
			}
		       vp2->a1 = a1 ; vp2->a2 = a2 ;
		    }
		}
	      else
		{
		  int a1 = vp2->a1, a2 = vp2->a2 ;
		  foundTranscript = TRUE ;
		  exon = 0 ; gene = 0 ;
		  while (baHit2Exon (ba, vp, &exon, &gene)) /* the first time we report the whole gene */
		    {
		      baExportIntergenicSupportRegister (ba, vp, isFirstFragment, genes, exons, gene, exon
							 , gsF, gsR, gsFg, gsRg
							 , aliF, aliR, aliFg, aliRg
							 , tagF, tagR, tagFg, tagRg
							 , errF, errR, errFg, errRg	
							 ) ;
		      nn++ ;
		    }
		  vp2->a1 = a1 ; vp2->a2 = a2 ;
		}
	    }
	}   
    }

  /* export */
  for (ii = 0 ; 0 && ii < arrayMax (gsF) ; ii++) /* this code is not yet debugged 2013_04_02 */
    {
      if (array (gsF, ii, double) > 0)
	{
	  aceOutf (ao, "Exon\t%s", dictName (ba->exonDict, ii)) ;
	  aceOutf (ao, "\t%s", unique ? "u" : "nu") ;
	  aceOutf (ao, "\t%s", ba->run) ;
	  aceOutf (ao, "\t%.1f\t%.2f\t%.0f\t%.0f", array (gsF, ii, double), array (tagF, ii, double), array (aliF, ii, double), array (errF, ii, double)) ;
	  aceOutf (ao, "\n") ;
	}
      
      if (array (gsR, ii, double) > 0)
	{
	  aceOutf (ao, "Exon\t%s", dictName (ba->exonDict, ii)) ;
	  aceOutf (ao, "\t%s", unique ? "u" : "nu") ;
	  aceOutf (ao, "\tAnti_%s", ba->run) ;
	  aceOutf (ao, "\t%.1f\t%.2f\t%.0f\t%.0f", array (gsR, ii, double), array (tagR, ii, double), array (aliR, ii, double), array (errR, ii, double)) ;
	  aceOutf (ao, "\n") ;
	}
    } 
  for (ii = 0 ; ii < arrayMax (gsFg) ; ii++)
    {
      if (array (gsFg, ii, double) > 0)
	{
	  aceOutf (ao, "Gene\t%s", dictName (ba->targetDict, ii)) ;
	  aceOutf (ao, "\t%s", unique ? "u" : "nu") ;
	  aceOutf (ao, "\t%s", ba->run) ;
	  aceOutf (ao, "\t%.1f\t%.2f\t%.0f\t%.0f", array (gsFg, ii, double), array (tagFg, ii, double), array (aliFg, ii, double), array (errFg, ii, double)) ;
	  aceOutf (ao, "\n") ;
	}
      
      if (array (gsRg, ii, double) > 0)
	{
	  aceOutf (ao, "Gene\t%s", dictName (ba->targetDict, ii)) ;
	  aceOutf (ao, "\t%s", unique ? "u" : "nu") ;
	  aceOutf (ao, "\tAnti_%s", ba->run) ;
	  aceOutf (ao, "\t%.1f\t%.2f\t%.0f\t%.0f", array (gsRg, ii, double), array (tagRg, ii, double), array (aliRg, ii, double), array (errRg, ii, double)) ;
	  aceOutf (ao, "\n") ;
	}
    } 
  if (ba->target_class)  /* remove the Z_genome hits */
    {
      for (ii = jj = 0, up = vp = arrp (ba->hits, 0, HIT) ; ii < arrayMax (ba->hits) ; up++, ii++)
	{
	  up2 = arrp (ba->hits2, up->nn, HIT2) ;
	  if (ba->target_class == up2->target_class)
	    {
	      if (jj < ii) *vp = *up ;
	      vp++ ; jj++ ;
	    }
	}
      arrayMax (ba->hits) = jj ;
    }
  ac_free (h) ;
  return nn ;
} /* baExportIntergenicSupport */

/*************************************************************************************/
/*************************************************************************************/

static int baCreateAtlas (BA *ba, BOOL getIntrons)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai ;
  Array atlas, map, targetAtlas, genome2geneAtlas = 0 ;
  BOOL isIntron ;
  char buf[1000], tcL ;
  char chromBuf[1000] ;
  char geneBuf[1000] ;
  char intronBuf[1064] ;
  char exonBuf[1064] ;
  const char *ccp ;
  MHIT *mp ;
  int ii, nn = 0, x1, x2, y1, y2, a1, a2, b1, b2, mrna, gene, target, target_class ;

  if (! ba->targetDict)
    ba->targetDict = dictHandleCreate (10000, ba->h) ;

  x1 = x2 = y1 = y2 = a1 = a2 = b1 = b2 = 0 ;
  if (getIntrons)
    {
      targetAtlas = ba->target2intronAtlas = arrayHandleCreate (4, Array, ba->h) ;
      ba->intronDict = dictHandleCreate (1000000, ba->h) ;
    }
  else
    {
      genome2geneAtlas = ba->genome2geneAtlas = arrayHandleCreate (256, Array, ba->h) ;
      targetAtlas = ba->target2exonAtlas = arrayHandleCreate (4, Array, ba->h) ;
      ba->exonDict = dictHandleCreate (1000000, ba->h) ;
    }
  ai = ba->mrnaRemapFileName ? aceInCreate (ba->mrnaRemapFileName, FALSE, h) : 0 ;
  if (ai) 
    aceInSpecial (ai, "\t\n") ;

  while (ai && aceInCard (ai))
    {
      ccp = aceInWord (ai) ; /* target_class */
      if (! ccp || *ccp == '#')
	continue ;
      dictAdd (ba->target_classDict, ccp, &target_class) ;
      tcL = *ccp ;
      if (target_class == ba->intronClass)
	continue ;
      atlas = array (targetAtlas, target_class, Array) ;
      if (! atlas)
	{
	  atlas = arrayHandleCreate (100000, Array, ba->h) ;
	  array (targetAtlas, target_class, Array) = atlas ;
	}
      ccp = aceInWord (ai) ; /* target */
      if (! ccp || *ccp == '#')
	continue ;
      isIntron =  strcmp (buf, ccp) ? FALSE : TRUE ;
      strncpy (buf, ccp, 999) ;
      if (! aceInInt (ai, &x1) || ! aceInInt (ai, &x2))  /* mRNA exon coordinates */
	continue ;
      ccp = aceInWord (ai) ; /* chrom */
      if (! ccp || *ccp == '#')
	continue ;  
      strncpy (chromBuf, ccp, 899) ; /* leave 100 bytes for the coordinates in intronBuf */
      if (! aceInInt (ai, &a1) || ! aceInInt (ai, &a2))  /* chrom exon coordinates */
	continue ;
      ccp = aceInWord (ai) ; /* gene */
      if (! ccp || *ccp == '#')
	continue ;     
      strncpy (geneBuf, ccp, 899) ;
      if (! getIntrons)
	{
	  dictAdd (ba->targetDict, buf, &mrna) ;
	  if (*geneBuf) dictAdd (ba->targetDict, geneBuf, &gene) ;
	  dictAdd (ba->targetDict, chromBuf, &target) ;

	  map = array (atlas, mrna, Array) ;
	  if (! map)
	    map = array (atlas, mrna, Array) = arrayHandleCreate (8, MHIT, ba->h) ;
	  nn++ ;

	  mp = arrayp (map, arrayMax (map), MHIT) ;
	  mp->gene = gene ;
	  mp->a1 = a1 ; mp->a2 = a2 ;

	  mp->class = tcL ;
	  mp->x1 = x1 ; mp->x2 = x2 ;

	  /* name the pieces, in order to ventilate the alignments by dictionary address */
	  if (a1 < a2)
	    {
	      sprintf (exonBuf, "%s__%d_%d", chromBuf, a1, a2) ;
	    }
	  else
	    {
	      sprintf (exonBuf, "%s__%d_%d", chromBuf, a2, a1) ;
	    }
	  dictAdd (ba->exonDict, exonBuf, &mp->exon) ;

	  map = array (genome2geneAtlas, target, Array) ;
	  if (! map)
	    map = array (genome2geneAtlas, target, Array) = arrayHandleCreate (8, MHIT, ba->h) ;

	  mp = arrayp (map, arrayMax(map), MHIT) ;
	  mp->gene = gene ;
	  if (a1 < a2)
	    { 
	      mp->x1 = a1 ; mp->x2 = a2 ;
	      mp->a1 = 1 ;
	    }
	  else
	    { 
	      mp->x1 = a2 ; mp->x2 = a1 ;
	      mp->a1 = -1 ;
	    }
	}
      if (getIntrons && isIntron && (a1 > b2 + 10 || a1 < b2 - 10))
	{
	  dictAdd (ba->targetDict, buf, &target) ;
	  if (*geneBuf) dictAdd (ba->targetDict, geneBuf, &gene) ;
	  map = array (atlas, target, Array) ;
	  if (! map)
	    map = array (atlas, target, Array) = arrayHandleCreate (8, MHIT, ba->h) ;
	  nn++ ;
	  mp = arrayp (map, arrayMax (map), MHIT) ;
	  mp->gene = gene ;

	  mp->class = 'U' ;
	  mp->x1 = y2 ; mp->x2 = x1 ;
	  mp->a1 = b2 ; mp->a2 = a1 ;
	  if (b2 < a1)
	    {
	      sprintf (intronBuf, "%s__%d_%d", chromBuf, b2 + 1, a1 - 1) ;
	    }
	  else
	    {
	      sprintf (intronBuf, "%s__%d_%d", chromBuf, b2 - 1 , a1 + 1) ;
	    }
	  dictAdd (ba->intronDict, intronBuf, &mp->exon) ;
	}
      b1 = a1 ; b2 = a2 ;      
      y1 = x1 ; y2 = x2 ;      
    }

  fprintf (stderr, "baCreateAtlas found %d introns in file %s\n"
	   , nn
	   , ba->mrnaRemapFileName ? ba->mrnaRemapFileName : "Missing remapFileName"
	   ) ;


  /* sort and fuse the transcriptboxes */
  if (genome2geneAtlas)
    for (ii = 0 ; ii < arrayMax (genome2geneAtlas) ; ii++)
    {
      MHIT *mp2 ;

      map = array (genome2geneAtlas, ii, Array) ;
      if (map)
	{
	  int i, j, k ;
	  arraySort (map, mhitGeneOrder) ;
	  for (i = j = 0, mp = arrp(map,i,MHIT) ; i < arrayMax (map) ; mp++, i++) 
	    {
	      for (k = i + 1, mp2 = mp + 1 ; mp->gene && k < arrayMax (map) && mp2->gene == mp->gene ; mp2++, k++) 
		{
		  if (mp2->x2 > mp->x2) 
		    {
		      if (mp->x2 > mp2->x1) mp2->x1 = mp->x2 ;
		      mp->x2 = mp2->x2 ; 
		      if (mp2->x1 > mp2->x2) { mp2->x1 = mp2->x1 ; }
		    }
		}
	     if (i > j) 
	       {
		 mp2 = arrp(map,j,MHIT) ;
		 *mp2 = *mp ;
	       }
	     j++;
	     i = k - 1 ;
	     mp = arrp(map,i,MHIT) ;
	    }
	  arrayMax (map) = j ;	      
	  arraySort (map, mhitOrder) ;
	  
	  if (0)
	    for (i = j = 0, mp = arrp(map,i,MHIT) ; i < arrayMax (map) ; mp++, i++) 
	      fprintf (stderr, "%s\t%s\t%d\t%d\t%d\t%d\n"
		       , dictName (ba->targetDict, ii)
		       , mp->gene ? dictName (ba->targetDict, mp->gene) : "-"
		       , mp->a1, mp->a2, mp->x1, mp->x2
		       ) ;
	}
    }

  ac_free (h) ;
  return nn ;
} /* baCreateAtlas */

/*************************************************************************************/
/* For each gene we have a map
 * For each target class these are gathered in an atlas
 */
static int baHit2Intron (BA *ba, HIT *vp, int *intronp, int *z1p, int *z2p, int *zx1p, int *zx2p, int *gDownp)
{
  int ii, intron = *intronp ;
  int dx = ba->minIntronOverlap - 1 ;
  HIT2 *vp2 ;
  MHIT *mp ;
  Array atlas, map ;
  
  if (dx < 0) dx = 0 ;
  vp2 = arrp (ba->hits2, vp->nn, HIT2) ;
  if (vp2->target_class == ba->intronClass)
    {
      const char *ccp = dictName(ba->targetDict, vp->target) ;
      if (! strncmp(ccp,"ZI_",3)) ccp += 3 ;
      dictAdd (ba->intronDict, ccp, &intron) ;
      *gDownp = 1 ;
      return intron ;
    }
  else
    {
      *intronp = 0 ;

      if (ba->target2intronAtlas && 
	  vp2->target_class < arrayMax (ba->target2intronAtlas) &&
	  (atlas =  arr (ba->target2intronAtlas, vp2->target_class, Array)) &&
	  vp->target < arrayMax(atlas) &&
	  (map = array (atlas, vp->target, Array))
	  )
	{
	  int a1 = vp2->a1, a2 = vp2->a2, x1 = *zx1p, x2 = *zx2p ;
	  *gDownp = 1 ;

	  if (a1 > a2) {  ii = a1 ; a1 = a2 ; a2 = ii ; *gDownp = -1 ; }
	  for (ii = 0, mp = arrp (map, 0, MHIT) ; ii < arrayMax (map)  ; ii++, mp++)
	    {
	      if (intron)
		{
		  if (x1 != mp->x1)
		    continue ;
		  else if (x2 != mp->x2)
		    return 0 ;
		  ii++ ;  mp++ ;
		}
	      if (ii < arrayMax (map) && a1 < mp->x1 - dx && a2 > mp->x2 + dx)
		{
		  *intronp = mp->exon ; 
		  *z1p = mp->a1 ; *z2p = mp->a2 ; *zx1p = mp->x1 ; *zx2p = mp->x2 ; 
		  return *intronp ;
		}
	      if (intron)
		break ;
	    }
	}
    }
  return 0 ;
} /* baHit2Intron */

/*************************************************************************************/
/* For each gene we have a map, and vp is already on the genome
 */
static int baHit2Gene (BA *ba, HIT *vp, int *genep)
{
  int ii, gene = 0, b1, b2 ;
  HIT2 *vp2 ;
  MHIT *mp ;
  Array map ;
  static int iiOld = -1 ;

  if (*genep == 0) iiOld = -1 ;
  vp2 = arrp (ba->hits2, vp->nn, HIT2) ;
  if (ba->genome2geneAtlas && 
      vp->target < arrayMax (ba->genome2geneAtlas) &&
      (map = array (ba->genome2geneAtlas, vp->target, Array))
      )
    {
      int a1 = vp2->a1, a2 = vp2->a2 ;
      BOOL isUp = FALSE ;
      
      if (a1 > a2) { isUp = TRUE ; ii = a1 ; a1 = a2 ; a2 = ii ; }
      if (1)
	{
	  for (ii = iiOld + 1, mp = arrp (map, 0, MHIT) ; ii < arrayMax (map) ; ii++, mp++)
	    {
	      b1 = a1 > mp->x1 ? a1 : mp->x1 ;
	      b2 = a2 < mp->x2 ? a2 : mp->x2 ;
	      if (b2 - b1 >= 8)
		{
		  gene = mp->gene ;
		  if (!isUp) /* consume the aligned part, to avoid overcounting */
		    vp2->a1 = b2 ; 
		  else
		    vp2->a2 = b2 ;
		  break ; 
		}
	      else if (a2 < mp->x1 + 7)  
		break ;
	    }
	}
    }
 if (*genep >= gene)
    gene = 0 ;
  *genep = gene ;
 
  return gene ? mp->a1 : 0 ; /* isDown */
} /* baHit2Gene */

/*************************************************************************************/
/* For each gene we have a map
 * For each target class these are gathered in an atlas
 */
static int baHit2Exon (BA *ba, HIT *vp, int *exonp, int *genep)
{
  int ii, gene = *genep, exon = *exonp, b1, b2 ;
  HIT2 *vp2 ;
  MHIT *mp ;
  Array atlas, map ;
  static int iiOld = -1 ;

  return 0 ; /* this code is not debugged */

  if (*exonp == 0) iiOld = -1 ;
  vp2 = arrp (ba->hits2, vp->nn, HIT2) ;
  if (ba->target2exonAtlas && 
      vp2->target_class < arrayMax (ba->target2exonAtlas) &&
      (atlas =  arr (ba->target2exonAtlas, vp2->target_class, Array)) &&
      vp->target < arrayMax(atlas) &&
      (map = array (atlas, vp->target, Array))
      )
    {
      int a1 = vp2->a1, a2 = vp2->a2 ;
      BOOL isUp = FALSE ;
      
      if (a1 > a2) { isUp = TRUE ; ii = a1 ; a1 = a2 ; a2 = ii ; }
      if (1)
	{
	  for (ii = iiOld + 1, mp = arrp (map, 0, MHIT) ; ii < arrayMax (map) ; ii++, mp++)
	    {
	      b1 = a1 > mp->x1 ? a1 : mp->x1 ;
	      b2 = a2 < mp->x2 ? a2 : mp->x2 ;
	      if (b2 - b1 >= 8 && mp->exon > exon)
		{
		  exon = mp->exon ; gene = mp->gene ; iiOld = ii ; break ; 
		  if (!isUp) /* consume the aligned part, to avoid overcounting */
		    vp2->a1 = b2 ; 
		  else
		    vp2->a2 = b2 ;
		}
	      else if (a2 < mp->x1 + 7)  
		break ;
	    }
	}
    }
  if (*genep >= gene)
    gene = 0 ;
  *genep = gene ;
  if (*exonp >= exon)
    exon = 0 ;
  *exonp = exon ;
  return exon ;
} /* baHit2Exon */

/*************************************************************************************/
/* foreach tag first best alignment
 * if in KT_RefSeq or ET_av or U_intron or any class with introns in the remap file
 * check if we overlap an intron (not caring for local errors since we clip) by 8bp using the remap coords
 * we count sequences, tags and bp (adding the aligned bp on both side)
 * the search is limited to cases where the 2 reads of the pair align in the same gene
 * separately export the double_introns
 */

static int baIntronSupport (BA *ba, BOOL unique) 
{  
  AC_HANDLE h = ac_new_handle () ;
  HIT *up, *vp, *oldVp = 0, *oldCloneVp = 0 ;
  HIT2 *up2, *vp2 ;
  int ii, jj, kk = 0, intron, oldF = 0, oldF1 = 0, oldCloneF = 0, dIntron, run, score = 0, tag = 0, ali, mult, nn = 0 ;
  int oldA1 = 0, oldA2 = 0, oldX1 = 0, oldX2 = 0, oldCloneA1 = 0, oldCloneA2 = 0, oldCloneX1 = 0, oldCloneX2 = 0, oldCloneZX2 = 0 ;
  int isFirstFragment ;
  Associator introns = assHandleCreate (h) ;
  double uu ;
  int str = ba->stranded ?  ba->stranded : 0 ;
  int z1 = 0, z2 = 0, zx1 = 0, zx2 = 0 ;
  BOOL hasPair = ba->pair > 0 || ba->hasPair ;
  int gDown = 1, mDown = 1 ;
  ACEOUT ao = aceOutCreate (ba->outFileName, unique ? ".intronSupport.u" : ".intronSupport.nu", ba->gzo, h) ;
  BOOL debug = FALSE ;
  DICT *doubleIntronDict = dictHandleCreate (100000, h) ;

  int nIntron = 0, ndIntron = 0 ;
  KEYSET iiF = arrayHandleCreate (100000, KEY, h) ;
  KEYSET iiTag = arrayHandleCreate (100000, KEY, h) ;
  Array gsF = arrayHandleCreate (100000, double,h) ;
  Array aliF = arrayHandleCreate (100000, double,h) ;
  Array tagF = arrayHandleCreate (100000, double,h) ;
  
  KEYSET iidF = arrayHandleCreate (100000, KEY, h) ;
  KEYSET iidTag = arrayHandleCreate (100000, KEY, h) ;
  Array dgsF = arrayHandleCreate (100000, double,h) ;
  Array dtagF = arrayHandleCreate (100000, double,h) ;
  Array daliF = arrayHandleCreate (100000, double,h) ;
  
  char rNam[2][1000] ;
  memset (rNam, 0, sizeof(rNam)) ; 

  if (! ba->run)
    messcrash ("-intronSupport requires option run run_name") ;

  dictAdd (ba->runDict, ba->run, &(run)) ;
  
  for (ii = score = tag = 0 ; ii < arrayMax (ba->hits) ; ii++)
    {
      up = arrp (ba->hits, ii, HIT)  ; /* do not use up++ since ii is incremented to jj-1 after the jj loop */
      up2 = arrp (ba->hits2, up->nn, HIT2) ;
      mDown = (up2->a1 < up2->a2 ? 1 : -1) ; gDown = 1 ; mDown = 1 ;
      if (hasPair &&  up2->dPair > NON_COMPATIBLE_PAIR &&  up2->dPair < 0 &&  up2->dPair != -2 &&  up2->dPair != -5)  /* synchronize to hack these reserved values with bestali.c */
	continue ;
      if (! up->score || tag == up->tag)
	continue ;

      if (1)
	{
	  tag = up->tag ; assClear (introns) ;
	  isFirstFragment = 1 ;
	  { 
	    const char *ccp = dictName (ba->tagDict,up->tag) ;
	    int i = strlen (ccp) ;
	    ccp += i - 1 ;
	    if (*ccp == '<')
	      { isFirstFragment = -1 ; i-- ; }
	     if (*ccp == '>')
	      { isFirstFragment = 1 ; i-- ; }
	     kk = (kk + 1) % 2 ;
	     if (i > 1000)
	       oldCloneF = 0 ;
	     else
	       {
		 strncpy (rNam[kk],  dictName (ba->tagDict,up->tag), i) ;
		 if (strcmp (rNam[0], rNam[1]))
		   oldCloneF = 0 ;
	       }
	  }
	  if (oldCloneVp && oldCloneVp->target != up->target) 
	    oldCloneF = 0 ;
	  intron = 0 ; oldF = 0 ;
	  for (vp = up, jj = ii ;  vp->tag == tag && jj < arrayMax (ba->hits) ; vp++, jj++)
	    {
	      vp2 = arrp (ba->hits2, vp->nn, HIT2) ;
	      oldF1 = 0 ;
	      if (unique && vp2->unicity != 1)
		continue ;
	      if (! unique && vp2->unicity == 1)
		continue ;
	      while (z1 = z2 = 0, baHit2Intron (ba, vp, &intron, &z1, &z2, &zx1, &zx2, &gDown))
		if (assInsert (introns, dictName (ba->intronDict, intron), 0)) /* so we export only once per tag */
		  {
		    mult = fastcMultiplicity (dictName (ba->tagDict,vp->tag), 0, 0) ;
		    uu = vp2->unicity ;
		    ali = vp2->ali ;
		    if (uu < 0) uu = -uu ;
		    if (uu == 0) uu = 1 ;
		    if (uu > 1) uu = 1 ; /* tant pis, un read map dans 2 genes differents on compte 1 les 2 fois ce qui resoud le fait que les deux genes partagent le meme intron si l'annotation de distingue pas genes et transcripts */
		    if (str * isFirstFragment * (vp2->a2 - vp2->a1) >= 0)
		      {
			if (debug) printf ("###INTRON %s\t%s\t%s\n",  dictName (ba->intronDict, intron), rNam[kk], dictName(ba->tagDict, up->tag)) ;
			keySet (iiF, nIntron) = intron ;
			keySet (iiTag, nIntron) = up->tag ;
			array (gsF, nIntron, double) += 1/uu ;
			array (tagF, nIntron, double) += mult/uu ;
			array (aliF, nIntron, double) += ali*mult/uu ;
			nIntron++ ;
			if (oldF && oldF != intron && 
			    oldVp && 
			    oldVp == vp  && vp->target == oldVp->target && 
			    vp->tag == oldVp->tag &&
			    vp->class == oldVp->class &&
			    vp->chain == oldVp->chain &&
			    (
			     (oldA1 < oldA2 && z1 < z2) ||
			     (oldA2 < oldA1 && z2 < z1) 
			    ) &&
			    (
			     (oldVp->x1 < vp->x2 &&  oldVp->x2 > vp->x1 - 3) ||
			     (oldVp->x2 > vp->x1 &&  oldVp->x1 < vp->x2 - 3) 
			     )
			    )
			  {
			    char *cp1 = strstr(dictName (ba->intronDict, intron),"__") ;  
			    char *cp2 = strstr(dictName (ba->intronDict, oldF),"__") ;
			    if (cp1)
			      {
				if ((oldA1 < z1 && z1 < z2) || (oldA1 > z1 && z1 > z2))
				  dictAdd (doubleIntronDict, messprintf("%s_%s", dictName (ba->intronDict,oldF), cp1), &dIntron) ;
				else
				  dictAdd (doubleIntronDict, messprintf("%s_%s", dictName (ba->intronDict,intron), cp2), &dIntron) ;

				if (debug) printf ("###DOUBLE_INTRON_uno %s\t%s\n",  dictName (doubleIntronDict, dIntron), rNam[kk]) ;
				keySet (iidF, ndIntron) = dIntron ;
				keySet (iidTag, ndIntron) = up->tag ;
				array (dgsF, ndIntron, double) += 1/uu ;
				array (dtagF, ndIntron, double) += mult/uu ;
				array (daliF, ndIntron, double) += ali*mult/uu ;
				ndIntron++ ;
			      }
			  }
			if (! oldF1) oldF1 = intron ;
			oldF = intron ; oldVp = vp ;
			if (debug) fprintf(stderr, "AAA oldF %d :: %s\t",  oldF, oldF ? dictName (ba->intronDict, oldF) : 0) ;
			if (debug) fprintf(stderr, "AAA oldCloneF %d :: %s oldCloneF=%d\n",  oldCloneF,  oldCloneF ? dictName (ba->intronDict, oldCloneF) : "", oldCloneF) ;
			if (vp2->a1 < vp2->a2) { oldX1 = vp2->a1 ; oldX2 = vp2->a2 ; }
			else  { oldX1 = vp2->a2 ; oldX2 = vp2->a1 ; }
			oldA1 = z1 ; oldA2 = z2 ;
		
			if (ba->pair && oldCloneVp &&  oldCloneVp->tag != vp->tag &&
			    oldCloneF && oldCloneF != intron && oldCloneVp && vp->target == oldCloneVp->target &&
			    ((oldCloneX1 < oldX2 && oldCloneX2 + 30 > oldX1) || (oldCloneX2 > oldX1 && oldCloneX1 - 30 < oldX2))
			    )
			  {
			    int dddx = oldCloneZX2 - zx2 ; if (dddx < 0) dddx = - dddx ;
			    int ddda = 0 ;
			    if (oldCloneA1 < oldCloneA2 && z1 < z2)
			      {
				if (oldCloneA2 < z1) ddda = z1 - oldCloneA2 ;
				if (z2 < oldCloneA1) ddda = oldCloneA1 - z2 ;
			      }
			    else if (oldCloneA1 > oldCloneA2 && z1 > z2)
			      {
				if (oldCloneA2 > z1) ddda = oldCloneA2 - z1 ;
				if (z2 > oldCloneA1) ddda = z2 - oldCloneA1 ;
			      }
			    if (ddda > 0 && ddda < dddx + 30)
			      {
				char *cp1 = strstr(dictName (ba->intronDict, intron),"__") ;
				char *cp2 = strstr(dictName (ba->intronDict, oldCloneF),"__") ;
				if (cp1)
				  {
				    if ((oldCloneA1 < z1 && z1 < z2) || (oldCloneA1 > z1 && z1 > z2))
				      dictAdd (doubleIntronDict, messprintf("%s_%s", dictName (ba->intronDict,oldCloneF), cp1), &dIntron) ;
				    else
				      dictAdd (doubleIntronDict, messprintf("%s_%s", dictName (ba->intronDict,intron), cp2), &dIntron) ;
				    
				    if (debug) printf ("###DOUBLE_INTRON_duo %s\t%s\n",  dictName (doubleIntronDict, dIntron), rNam[kk]) ;
				    keySet (iidF, ndIntron) = dIntron ;
				    keySet (iidTag, ndIntron) = up->tag ;
				    array (dgsF, ndIntron, double) += 1/uu ;
				    array (dtagF, ndIntron, double) += mult/uu ;
				    array (daliF, ndIntron, double) += ali*mult/uu ;
				    ndIntron++ ;
				  }
				if (0) oldCloneF = 0 ;
			      }
			  }
		      }		  
		    if (! oldCloneF  || ! oldCloneVp ||  (oldCloneVp->tag == tag && zx2 - oldCloneZX2) * mDown * gDown > 0)
		      {
			oldCloneVp = oldVp ;
			oldCloneA1 = oldA1 ; oldCloneA2 = oldA2 ; oldCloneX1 = oldX1 ; oldCloneX2 = oldX2 ; oldCloneZX2 = zx2 ;
			oldCloneF = oldF ;
		      }
		    if (debug) fprintf(stderr, "BBB oldClone %d :: %s oldCloneF = %d\n",  oldCloneF,  oldCloneF ? dictName (ba->intronDict, oldCloneF) : "-", oldCloneF) ;
		  }
		else
		  oldVp = 0 ;
	      if (debug) fprintf(stderr, "DDD1 oldClone %d :: %s oldCloneF = %d\n",  oldCloneF,  oldCloneF ? dictName (ba->intronDict, oldCloneF) : "-", oldCloneF) ;
	      if (! oldCloneF  || ! oldCloneVp ||  (oldCloneVp->tag == tag && zx2 - oldCloneZX2) * mDown * gDown > 0)
		{
		  oldCloneVp = oldVp ;
		  oldCloneA1 = oldA1 ; oldCloneA2 = oldA2 ; oldCloneX1 = oldX1 ; oldCloneX2 = oldX2 ; oldCloneZX2 = zx2 ;
		  oldCloneF = oldF ;
		}  
	      if (debug) fprintf(stderr, "DDD2 oldClone %d :: %s oldCloneF = %d\n",  oldCloneF,  oldCloneF ? dictName (ba->intronDict, oldCloneF) : "-", oldCloneF) ;
	    }
	  ii = jj - 1 ;
	  nn ++ ;
	  if (! oldCloneF || ! oldCloneVp  || (oldCloneVp && oldCloneVp->tag == tag && zx2 - oldCloneZX2) * mDown * gDown > 0)
	    {
	      oldCloneVp = oldVp ;
	      oldCloneA1 = oldA1 ; oldCloneA2 = oldA2 ; oldCloneX1 = oldX1 ; oldCloneX2 = oldX2 ; oldCloneZX2 = zx2 ;
	      oldCloneF = oldF ;
	    }
	  if (debug) fprintf(stderr, "CCC oldClone %d :: %s  oldCloneF = %d\n",  oldCloneF,  oldCloneF ? dictName (ba->intronDict, oldCloneF) : "-", oldCloneF) ;
	}  
    }

  /* export */
  for (ii = 0 ; ii < arrayMax (gsF) ; ii++)
    {
      if (array (gsF, ii, double) > 0)
	{
	  int intron = keySet (iiF, ii) ;
	  int tag = keySet (iiTag, ii) ;
	  aceOutf (ao, "Intron\t%s", dictName (ba->intronDict, intron)) ;
	  aceOutf (ao, "\t%s", unique ? "u" : "nu") ;
	  aceOutf (ao, "\t%s", ba->run) ;
	  aceOutf (ao, "\t%.1f\t%.2f\t%.0f", array (gsF, ii, double), array (tagF, ii, double), array (aliF, ii, double)) ;
	  aceOutf (ao, "\t%s", dictName (ba->tagDict,tag)) ;
	  aceOutf (ao, "\n") ;
	}
    } 
  
  for (ii = 0 ; ii < arrayMax (dgsF) ; ii++)
    {
      if (array (dgsF, ii, double) > 0)
	{
	  int dIntron = keySet (iidF, ii) ;
	  int tag = keySet (iidTag, ii) ;
	  if (1) aceOutf (ao, "DoubleIntron\t%s", dictName (doubleIntronDict,dIntron)) ;
	  aceOutf (ao, "\t%s", unique ? "u" : "nu") ;
	  aceOutf (ao, "\t%s", ba->run) ;
	  aceOutf (ao, "\t%.1f\t%.2f\t%.0f", array (dgsF, ii, double), array (dtagF, ii, double), array (daliF, ii, double)) ;
	  if (1) aceOutf (ao, "\t%s", dictName (ba->tagDict,tag)) ;
	  aceOutf (ao, "\n") ;
	}
    }
  ac_free (h) ;
  return nn ;
} /* baIntronSupport */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/* given a 8 column table produced by the scripts elements.hits2gene.tcsh 
 *  KT_RefSeq   X__NAP1L1 u  LIF_S      lung    820.0   6651.00 491792
 * produce an ace file per gene with the index of expression per tissue
 */
typedef struct sg2iStruct { int target_class, gene, isU, run ; double seq, tag, bp, nReads, nReadsOk, err, a2g, partial, orphan, badTopo, multi, multi2, index ; KEYSET wiggle ;} GS2I ;

static void showGS2I (Array aa) 
{
  int ii ;
  GS2I *up ;

  if (aa) 
    {
      for (ii = 0 ; ii < arrayMax (aa) ; ii++)
	{
	  up = arrp (aa, ii, GS2I) ;
	  printf ("%d\t%d\t%d\t%d\t%d\ts=%.2f\tt=%.2f\tbp=%.2f\tz=%.2f\terr=%.2f\tpartial=%.2f\torphan=%.2f\n"
		  , ii
		  , up->target_class, up->gene, up->isU, up->run
		  , up->seq, up->tag, up->bp, up->index
		  , up->err, up->partial, up->orphan
		  ) ;
	}
    }
} /* showGS2I */

/*************************************************************************************/

static int sg2iOrder (const void *va, const void *vb)
{
  const GS2I *a = (const GS2I *) va, *b = (const GS2I *) vb ;
  int n ;
  
  n = a->target_class - b->target_class ; if (n) return n ;
  n = a->gene - b->gene ; if (n) return n ;
  n = a->run - b->run ; if (n) return n ;
  n = a->isU - b->isU ; if (n) return n ;

  return 0 ;
} /* sg2iOrder */

/*************************************************************************************/
/* to be set before calling arraySort */
static DICT *myGeneDict = 0 ; 
static DICT *myRunDict = 0 ; 

static int baElementSupport2aceAlphaOrder (const void *va, const void *vb)
{
  const GS2I *a = (const GS2I *) va, *b = (const GS2I *) vb ;
  int n ;
  
  n = a->target_class - b->target_class ; if (n) return n ;
  n = a->gene - b->gene ; if (n) return lexstrcmp (dictName(myGeneDict,a->gene), dictName(myGeneDict,b->gene)) ;
  n = a->run - b->run ; if (n) return lexstrcmp (dictName(myRunDict,a->run), dictName(myRunDict,b->run)) ;
  n = a->isU - b->isU ; if (n) return n ;

  return 0 ;
} /* baElementSupport2aceAlphaOrder */

/*************************************************************************************/
/* cumulate identical counts */
static int baElementSupport2aceCoalesce (Array aa, AC_HANDLE h)
{
  int ii, jj, nn ;
  GS2I *up, *vp ;
  
  arraySort (aa, sg2iOrder) ;
  for (ii = 0, up = arrp (aa, ii, GS2I), nn = arrayMax (aa) ; ii < nn ; up++, ii++)
    if (up->seq)
      {
	for (jj = ii + 1, vp = up + 1 ; 
	     jj < nn && up->run == vp->run &&
	       up->isU == vp->isU &&
	       up->gene == vp->gene && up->target_class == vp->target_class ; 
	     vp++, jj++)
	  if (vp->seq)
	    {
	      up->seq += vp->seq ;
	      up->tag += vp->tag ;
	      up->bp += vp->bp ;

	      up->nReads += vp->nReads ;
	      up->nReadsOk += vp->nReadsOk ;
	      up->err += vp->err ;
	      up->a2g += vp->a2g ;
	      up->partial += vp->partial ;
	      up->orphan += vp->orphan ;
	      up->badTopo += vp->badTopo ;
	      up->multi += vp->multi ;
	      up->multi2 += vp->multi2 ;

	      vp->seq = 0 ;
	      if (vp->wiggle)
		{
		  if (! up->wiggle)  /* transfer */
		    {
		      up->wiggle = vp->wiggle ;
		      vp->wiggle = 0 ;
		    }
		  else    /* accumulate */
		    {
		      int i = keySetMax (vp->wiggle) ;
		      KEY *ip, *jp ;
		      KEYSET ksu = up->wiggle, ksv = vp->wiggle ;

		      keySet (ksu, i - 1) += 0 ;
		      ip = arrp (ksu, 0, KEY) ;
		      jp = arrp (ksv, 0, KEY) ;
		      while (i--)
			*ip++  += *jp++ ;
		    }
		}
	    }
      }
  
  /* keep happy few */
  for (ii = jj = 0, up = vp = arrp (aa, ii, GS2I), nn = arrayMax (aa) ; ii < nn ; up++, ii++)
    {
      if (up->seq > 0)
	{
	  if (vp < up) *vp = *up ;
	  vp++ ; jj++ ; 
	}
      else
	{
	  if (up->wiggle)
	    keySetDestroy (up->wiggle) ;
	}
    }
  fprintf (stderr, "// baElementSupport2aceCoalesce  %d gene support lines coalesced to %d\n", arrayMax (aa), jj) ;
  arrayMax (aa) = jj ;
  return jj ;
} /* baElementSupport2aceCoalesce */

/*************************************************************************************/

static int baElementSupport2ace (BA *ba, BOOL isGene, BOOL isMrna, BOOL isIntron)
{
  AC_HANDLE h = ac_new_handle () ;
  GS2I *up, *vp, *upAny ;
  const char *ccp ;
  int ii, jj, k0, nn = 0, nn1, nn2 = 0, anyGene[100] ;
  int nInFile = 0 ;
  float zf ;
  KEYSET badGenes = keySetHandleCreate (h) ;
  ACEIN ai = 0 ;
  ACEOUT ao = 0, ao2 = 0 ;
  DICT *target_classDict = dictHandleCreate (10, h) ;
  DICT *geneDict ;
  DICT *runDict = dictHandleCreate (100, h) ;
  Array aa = arrayHandleCreate (1000000, GS2I, h) ;
  char *antiRun ;

  if (ba->run)
    antiRun = hprintf (h, "Anti_%s", ba->run) ;
  geneDict = ba->geneDict ;
  if (! geneDict)
    geneDict = ba->geneDict = dictHandleCreate (100000, ba->h) ;
  if (isGene)
    ao = aceOutCreate (ba->outFileName, ".geneSupport.ace", ba->gzo, h) ;
  else if (isMrna)
    ao = aceOutCreate (ba->outFileName, ".mrnaSupport.ace", ba->gzo, h) ;
  else if (isIntron)
    {
      ao = aceOutCreate (ba->outFileName, ".intronSupport.ace", ba->gzo, h) ;
      ao2 = aceOutCreate (ba->outFileName, ".doubleIntronSupport.ace", ba->gzo, h) ;
    }
  else
    messcrash ("baElementSupport2ace is neither isGene nor isMrna nor isIntron") ;

  for (ii = 0 ; ii < 100 ; ii++)
    if (isGene)
      dictAdd (geneDict, messprintf("000anyGene_%d",ii), &anyGene[ii]) ;
    else if (isMrna)
      dictAdd (geneDict, messprintf("000anyTranscript_%d",ii), &anyGene[ii]) ;
    else if (isIntron)
      dictAdd (geneDict, messprintf("000anyIntron_%d",ii), &anyGene[ii]) ;

  /* parse the 10 columns file */
  while (1)
    {
      if (! ba->inFileList)
	ai = aceInCreate (ba->hitFileName, ba->gzi, h) ;
      else
	ai = baNextInFile (ba, nInFile++, h) ;
      if (ai)
	while (aceInCard (ai))
	  {
	    up = arrayp (aa, nn, GS2I) ;
	    nn1 = nn ; /* remember this offset */
	    
	    if (! (ccp = aceInWord (ai)))
	      continue ;
	    dictAdd (target_classDict, ccp, &(up->target_class)) ; 
	    aceInStep (ai, '\t') ;
	    
	    if (! (ccp = aceInWord (ai)))
	      continue ;
	    dictAdd (geneDict, ccp, &(up->gene)) ;
	    aceInStep (ai, '\t') ;
	    if (! (ccp = aceInWord (ai)))
	      continue ;
	    if (! strcmp (ccp, "u"))
	      up->isU = TRUE ;
	    else if (! strcmp (ccp, "nu"))
	      up->isU = FALSE ;
	    else
	      continue ;
	    aceInStep (ai, '\t') ;
	    if (! (ccp = aceInWord (ai)))
	      continue ;
	    if (ba->run) 
	      {
		if (! strncasecmp (ccp, "Anti_", 5))
		  ccp = antiRun ; /* force the support to count for that run */
		else
		  ccp = ba->run ; /* force the support to count for that run */
	      }
	    dictAdd (runDict, ccp, &(up->run)) ;
	    aceInStep (ai, '\t') ;
	    if (! aceInFloat (ai, &zf))
	      continue ;
	    up->seq = zf ;  /* zf is float, up->.. is double so we can cumul */
	    aceInStep (ai, '\t') ;
	    if (! aceInFloat (ai, &zf))
	      continue ; \
	    up->tag = zf ;
	    aceInStep (ai, '\t') ;
	    if (! aceInFloat (ai, &zf))
	      continue ; 
	    up->bp = zf ;

	    /* to be back compatible, we must not 'continue' if the 3 next columns are missing 2014-05-19 */
	    aceInStep (ai, '\t') ;
	    if (! aceInFloat (ai, &zf)) 
	      zf = 0 ;
	    up->nReads = zf ;

	    aceInStep (ai, '\t') ;
	    if (! aceInFloat (ai, &zf)) 
	      zf = 0 ;
	    up->nReadsOk = zf ;

	    aceInStep (ai, '\t') ;
	    if (! aceInFloat (ai, &zf)) 
	      zf = 0 ;
	    up->err = zf ;

	    aceInStep (ai, '\t') ;
	    if (! aceInFloat (ai, &zf)) 
	      zf = 0 ;
	    up->a2g = zf ;

	    aceInStep (ai, '\t') ;
	    if (! aceInFloat (ai, &zf)) 
	      zf = 0 ;
	    up->partial = zf ;
	    aceInStep (ai, '\t') ;
	    if (! aceInFloat (ai, &zf)) 
	      zf = 0 ;
	    up->orphan = zf ;
	    aceInStep (ai, '\t') ;
	    if (! aceInFloat (ai, &zf)) 
	      zf = 0 ;
	    up->badTopo = zf ;

	    aceInStep (ai, '\t') ;
	    if (! aceInFloat (ai, &zf)) 
	      zf = 0 ;
	    up->multi = zf ;
	    aceInStep (ai, '\t') ;
	    if (! aceInFloat (ai, &zf)) 
	      zf = 0 ;
	    up->multi2 = zf ;
 
	    aceInStep (ai, '\t') ;
	    if (isMrna &&
		(ccp = aceInWord (ai)) &&
		! strcmp (ccp, "WIGGLE_START") &&
		aceInStep (ai, '\t') &&
		aceInInt (ai, &k0)
		)
	      {
		float zk ;
		KEYSET wiggle ;

		k0 /= MRNA_WIGGLE_STEP ;
		wiggle = up->wiggle = keySetHandleCreate (h) ;
		while (aceInStep (ai, '\t') && aceInFloat (ai, &zk))
		  keySet (wiggle, k0++) = 100 * zk ;
	      }
	    /* accept and duplicate this line */     
	    vp = arrayp (aa, ++nn, GS2I) ; ++nn ;
	    up = arrayp (aa, nn1, GS2I) ;  /* aa may be relocated */
	    *vp = *up ;
	    if (up->wiggle)
	      vp->wiggle = keySetHandleCopy (vp->wiggle, h) ;
	    if (up->target_class > 99)
	      messcrash ("up2->target_class = %d > 99, %s please modify the source code"
			 , up->target_class
			 , dictName (target_classDict, up->target_class)
			 ) ;
	    vp->gene = anyGene[up->target_class] ;
	    
	    nn2 += 2 ;
	    if (nn2 > 10000000)
	      {
		/* cumulate identical counts */
		nn = baElementSupport2aceCoalesce (aa, h) ;
		nn2 = 0 ;
	      }
	  }
      if (! ai || ! ba->inFileList)
	break ;
    }

  nn = baElementSupport2aceCoalesce (aa, h) ;

  /*   ntc = 1 + dictMax(target_classDict) ; */

  myGeneDict = geneDict ;
  myRunDict = runDict ;
  arraySort (aa, baElementSupport2aceAlphaOrder) ;

  /* eliminate the genes contributing over 2 per cent (was 1 per cent) from the cumul */
  if (1 && dictMax(geneDict) > 10000)
    for (ii = 0, up = arrp (aa, ii, GS2I), nn = arrayMax (aa) ; ii < nn ; up++, ii++)
      {
	/* locate the cumul gene */
	for (jj = 0, upAny = arrp (aa, jj, GS2I) ; jj < nn ; upAny++, jj++)
	  if (upAny->gene == anyGene[up->target_class] && upAny->run == up->run)
	    break ;
	if (jj == nn || up == upAny)
	  continue ;
	if (50 * up->bp > upAny->bp && upAny->bp > 100000000
	    )
	  {
	    upAny->bp -= up->bp ;
	    upAny->tag -= up->tag ;
	    upAny->seq -= up->seq ;
	    upAny->nReads -= up->nReads ;
	    upAny->nReadsOk -= up->nReadsOk ;
	    upAny->err -= up->err ;
	    upAny->a2g -= up->a2g ;
	    upAny->partial -= up->partial ;
	    upAny->orphan -= up->orphan ;
	    upAny->badTopo -= up->badTopo ;
	    upAny->multi -= up->multi ;
	    upAny->multi2 -= up->multi2 ;
	    keySet (badGenes, keySetMax(badGenes)) = up->gene ;
	  }
      }

  if (keySetMax (badGenes))
    {
      keySetSort (badGenes) ;
      keySetCompress (badGenes) ;
      fprintf (stderr, "Number_over_expressed_gene\t%d\n", keySetMax (badGenes)) ;
      for (ii = 0 ; ii < keySetMax (badGenes) ; ii++)
	fprintf (stderr, "Over_expressed_gene\t%s\n", dictName (geneDict, keySet (badGenes, ii))) ;
    }
	

  /* export */
  if (arrayMax(aa))
    for (ii = 0, up = arrp (aa, ii, GS2I), nn = arrayMax (aa) ; ii < nn ; up++, ii++)
      {
	ACEOUT myao = ao ;
	
	if (! up->seq)
	  continue ;
	
	if (isGene)
	  {
	    if (up->gene > 100)
	      aceOutf (ao, "Gene %s\n", dictName (geneDict, up->gene)) ;
	    else
	      aceOutf (ao, "Transcript G_Any_%s_gene // %s\n", dictName (target_classDict, up->target_class), dictName (geneDict, up->gene)) ;
	    ccp = strstr (dictName (target_classDict, up->target_class), "_") ;
	    if (ccp)
	      aceOutf (ao, "%s\n", ccp + 1) ;
	  }
	else if (isMrna)
	  {
	    if (up->gene > 100)
	      aceOutf (ao, "Transcript %s\n", dictName (geneDict, up->gene)) ;
	    else
	      aceOutf (ao, "Transcript G_Any_%s_transcript // %s\n", dictName (target_classDict, up->target_class), dictName (geneDict, up->gene)) ;
	    ccp = strstr (dictName (target_classDict, up->target_class), "_") ;
	    if (ccp)
	      aceOutf (ao, "%s\n", ccp + 1) ;
	  }
	else if (isIntron)
	  {
	    if (! strcasecmp (dictName (target_classDict, up->target_class), "DoubleIntron"))
	      myao = ao2 ;
	    
	    if (up->gene < 100)
	      aceOutf (myao, "%s G_Any_%s // %s\n", dictName (target_classDict, up->target_class)
		       , dictName (target_classDict, up->target_class)
		       , dictName (geneDict, up->gene)) ;
	    else
	      aceOutf (myao, "%s %s\n", dictName (target_classDict, up->target_class), dictName (geneDict, up->gene)) ;
	  }
	for (jj = ii, vp = up ; 
	     jj < nn && up->gene == vp->gene && up->target_class == vp->target_class ; 
	     vp++, jj++)
	  if (vp->seq)
	    {
	      const char *ccp, *ccq ;
	      ii = jj ; up= vp ;
	      
	      ccp = dictName (runDict, vp->run) ;
	      if (!strncmp (ccp, "Anti_", 5))
		{
		  if (! vp->isU)
		    continue ;
		  ccp += 5 ;
		  ccq = "Anti_run" ;
		}
	      else
		{
		  ccq = vp->isU ? "Run_U" : "Run_nU" ;
		}
	      
	      aceOutf (myao, "%s \"%s\" %.2f %.2f seqs %.2f tags %.2f kb %.2f reads %.2f compRead  %.2f err  %.2f a2g %.2f partial %.2f orphan %.2f badTopo %.2f Multi %.2f AmbStrand"
		       , ccq, ccp
		       , vp->index, vp->seq, vp->tag, vp->bp/1000.0
		       , vp->nReads, vp->nReadsOk
		       , vp->err, vp->a2g		       
		       , vp->partial, vp->orphan, vp->badTopo
		       , vp->multi, vp->multi2
		       ) ;
	      if (isMrna && up->wiggle && keySetMax (up->wiggle))
		{
		  int k, kMax = keySetMax (up->wiggle) ;
		  KEY *kp = arrp (up->wiggle, 0, KEY) ;
		  
		  aceOutf (ao, "\tWIGGLE_START") ;
		  for (k = 0 ; k < kMax && *kp == 0 ; k++, kp++) ;
		  aceOutf (ao, "\t%d", MRNA_WIGGLE_STEP * k) ;
		  for ( ; k < kMax ; k++, kp++)
		    if (*kp == 100 * ((*kp) / 100))
		      aceOutf (ao, "\t%d", *kp/100) ;
		    else
		      aceOutf (ao, "\t%.1f", *kp/100.0) ;
		}
	      aceOutf (myao, "\n") ;
	    }
	
	aceOutf (myao, "\n") ;
      }   
  
  ac_free (h) ;

  return 0 ;
} /* baElementSupport2ace */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

static int xOrder (const void *a, const void *b)
{
  const POINT2D *up = (const POINT2D *)a, *vp = (const POINT2D *)b ;
  float dx ;

  dx = up->x - vp->x ;
  if (dx < 0) return 1 ; /* big first */
  if (dx > 0) return -1 ; /* big first */
  
  return 0 ;
} /* xOrder */

/*************************************************************************************/
extern AC_OBJ ac_key2obj (AC_DB db, KEY key, BOOL fillIt, AC_HANDLE handle) ;

static void baErccOptimize (BA *ba)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ; 
  AC_DB db = ba->db ;
  AC_ITER iter = ac_query_iter (ba->db, TRUE, "find gene  Deep AND ERCC", 0, h) ;
  AC_OBJ gene = 0 ;
  AC_TABLE dr = 0 ;
  KEY run2 ;
  int kk = 0, ir, mix, ng = 0, pp ;
  const char *error = 0 ;
  char *cp ;
  vTXT txt = vtxtHandleCreate (h) ;
  vTXT txt2 = vtxtHandleCreate (h) ;
  const char *runNam[] = { "Run_nU", "Run_U", "Group_nU", "Group_U", 0} ;
  Array kbA = arrayHandleCreate (200, float, h) ;
  DICT *kbDict = dictHandleCreate (1000, h) ;
  KEYSET runs = keySetHandleCreate (h) ;	  
  KEYSET mixes = keySetHandleCreate (h) ;	  
  KEYSET genes = keySetHandleCreate (h) ;	  
  POINT2D *xyp ;
  float kb ;
  Array xy = arrayHandleCreate (300, POINT2D, h) ;
  ACEOUT ao = aceOutCreate (ba->outFileName, ".ERCC.correl.txt", FALSE, h) ;
  aceOutf (ao, "#%s\n", timeShowNow ()) ;
  aceOutf (ao, "#Correlation between the run anf the best fitting ERCC mix\n") ;
  aceOutf (ao, "#Run\tERCC\tR\tR^2\tshift\n") ;
  
  if (iter)
    while (ac_free (gene), ac_free (h1), gene = ac_iter_obj (iter))
      {
	h1 = ac_new_handle () ;
	ng++ ;
	keySet (genes, ng) = ac_obj_key (gene) ;
	dr =  ac_tag_table (gene, "ERCC", h1) ;
	for (ir = 0 ; dr && ir < dr->rows ; ir++)
	  {
	    run2 = ac_table_key (dr, ir, 1, 0) ;
	    kb = ac_table_float (dr, ir, 3, 0) ;
	    keySetInsert (mixes, run2) ;
	    cp = messprintf("Mix_%d_%d", ng, run2) ;
	    dictAdd (kbDict, cp , &kk) ;
	    array (kbA, kk, float) += kb ;
	  }
	for (pp = 0 ; runNam[pp] ; pp++)
	  {
	    dr = ac_tag_table (gene, runNam[pp], h1) ;
	    for (ir = 0 ; dr && ir < dr->rows ; ir++)
	      {
		run2 = ac_table_key (dr, ir, 0, 0) ;
		keySetInsert (runs, run2) ;
		kb = ac_table_float (dr, ir, 1, 0) ;
		
		cp = messprintf("Idx_%d_%d_%d", ng, pp, run2) ;
		dictAdd (kbDict, cp, &kk) ;
		array (kbA, kk, float) += kb ;
	      }
	  }
      }

  /* for each run, opimize the correlation with the known mix */
  if (ng < 10) 
    goto done ;
  
  for (ir = 0 ; ir < keySetMax (runs) ; ir++) /* each run is independent */
    for (pp = 0 ; runNam[pp] ; pp++) /* U and nU are optimized independently */
      {
	    
	/* compute the correlation coef between all the index  */
	int bestMix = -1, nxy = 0 ;
	double z1, z2, besta = 0, bestR = -99999, bestDy = 0, a = 0, b = 0, r = 0, w = 0 ;

	run2 = keySet (runs, ir) ;
	for (mix = 0 ; mix < keySetMax (mixes) ; mix++)
	  {
	    nxy = keySetMax (xy) = 0 ;
	    for (ng = 1 ; ng < keySetMax (genes) ; ng++)
	      {
		z1 = z2 = -999 ;
		cp = messprintf("Mix_%d_%d", ng, keySet(mixes,mix)) ;
		if (dictFind (kbDict, cp, &kk))
		  z1 = array (kbA, kk, float) ;
		cp = messprintf("Idx_%d_%d_%d", ng, pp, run2) ;
		if (dictFind (kbDict, cp, &kk))
		  z2 = array (kbA, kk, float) ;
		if (!strstr (ac_name(gene),"00154") &&  /* bad ercc to avoid in the fit */
		    !strstr (ac_name(gene),"00022") &&  /* bad ercc to avoid in the fit */
		    z2 > - 800 && z1 > 3 && z2 > 3)   /* measurable according to standard AceView Index */
		  {
		    xyp = arrayp (xy, nxy++, POINT2D) ;
		    xyp->x = z1 ; xyp->y = z2 ;
		    if (0) aceOutf(ao, "%s\t%s\t%.3f\t%.2f\t%s\n", name(run2), name( keySet(mixes,mix)), z1, z2, name(keySet(genes,ng))) ;
		  }
	      }	  
	    /* recover the top half of the ERCC */
	    arraySort (xy, xOrder) ;
	    arrayMax (xy) = nxy = nxy/1 ;
	    if (nxy > 0)
	      {
		/* Solves, y = a x + b , in a and b   a defaults to *ap  if too few data given */
		/* BOOL linearRegression (Array xy, double *ap, double *bp, double *rp, double *wp) */
		a = 1 ; r = 0 ;
		linearRegression (xy, &a, &b, &r, &w) ;
		if (r > bestR && r > .5)
		  { bestMix = mix + 1 ; bestR = r ; bestDy = -b ; besta = a ;}
	      }
	  }
	if (bestMix >= 0) 
	  {
	    for (ng = 0 ; ng < keySetMax (genes) ; ng++)
	      if (dictFind (kbDict, messprintf("Idx_%d_%d_%d", ng, pp, run2), &kk) && 
		  array (kbA, kk, float) > -800)
		array (kbA, kk, float) += bestDy ;
	
	    aceOutf (ao, "%s\t%s\t%.3f\t%.3f\t%.2f x + %.2f nxy=%d\n"
		     , name(run2), name (keySet(mixes, bestMix - 1))
		     , bestR
		     , bestR * bestR
		     , besta, -bestDy, nxy
		     ) ;
	    vtxtPrintf (txt2, "Ali \"%s\"\n-D ERCC_mix \nERCC_mix \"%s\"  %.3f %.3f \"R,R^2\"\n\n"
			, name(run2), name (keySet(mixes, bestMix - 1))
			, bestR
			, bestR * bestR
			) ;
	  }
      }
   /* parse the results */
   for (ng = 1 ; ng < keySetMax (genes) ; ng++)
     {
       ac_free (h1) ;
       h1 = ac_new_handle () ;
       vtxtClear(txt) ;
       /* be very careful not to destroy the read counts to the right of the index */
       gene = ac_key2obj (ba->db, keySet (genes, ng), TRUE, h1) ;
       vtxtPrintf (txt, "Gene \"%s\"\n", ac_name(gene)) ;
       for (pp = 0 ; runNam[pp] ; pp++)
	 {
	    dr = ac_tag_table (gene, runNam[pp], h1) ;
	    for (ir = 0 ; dr && ir < dr->rows ; ir++)
	      {
		int ix ;

		run2 = ac_table_key (dr, ir, 0, 0) ;
		vtxtPrintf (txt, "%s %s", runNam[pp], name (run2)) ;
		if (dictFind (kbDict, messprintf("Idx_%d_%d_%d", ng, pp, run2), &kk))
		  kb = array (kbA, kk, float) ;
		else
		  kb = ac_table_float (dr, ir, 1, 0) ;
		vtxtPrintf (txt, " %.2f ", kb) ;

		kb = ac_table_float (dr, ir, 2, 0) ;
		ix = 100 * kb ;
		if (ix % 100 == 0) vtxtPrintf (txt, " %d seqs", ix/100) ;
		else vtxtPrintf (txt, " %.2f seqs", kb) ;

		kb = ac_table_float (dr, ir, 4, 0) ;
		ix = 100 * kb ;
		if (ix % 100 == 0) vtxtPrintf (txt, " %d tags", ix/100) ;
		else vtxtPrintf (txt, " %.2f tags", kb) ;

		kb = ac_table_float (dr, ir, 6, 0) ;
		ix = 100 * kb ;
		if (ix % 100 == 0) vtxtPrintf (txt, " %d kb", ix/100) ;
		else vtxtPrintf (txt, " %.2f kb", kb) ;

		vtxtPrintf (txt, "\n") ;
	      }
	 }
       vtxtPrintf (txt, "\n") ;
       if (vtxtPtr (txt2)) vtxtPrint (txt, vtxtPtr (txt2)) ;
       if (vtxtPtr (txt))
	 {
	   if (0)
	     fprintf (stderr, "%s", vtxtPtr (txt)) ;
	   else
	     {
	       ac_parse (db, vtxtPtr (txt), &error, 0, h) ; 
	       if (error && *error)
		 fprintf (stderr, "baErccOptimize error=%s\n", error) ;
	     }
	 }
    }
   ac_free (h1) ;
 done:
   ac_free (ao) ;
   ac_free (h) ;
} /* baErccOptimize */

/*************************************************************************************/
/*************************************************************************************/
/* 2016_06_22
 * merge inside the Run objects
 */
static int baMergeObservedStrandedness (vTXT txt, AC_OBJ Group, int groupLevel)
{
  AC_HANDLE h1 = 0 , h = ac_new_handle () ; 
  AC_TABLE runs ;
  AC_TABLE tt = 0 ;
  int irun, target = 0 ;
  DICT *dict = 0 ;
  float fp[16], fm[16], fa[16] ;
 
  if (groupLevel == 0)
    runs = ac_tag_table (Group, "sublibraries", h) ;
  else
    runs = ac_tag_table (Group, "Union_of", h) ;
  
  if (! runs || ! runs->rows)
    goto done ;
  
  dict = dictHandleCreate (8, h) ;
  memset (fp, 0, sizeof(fp)) ;
  memset (fm, 0, sizeof(fm)) ;
  memset (fa, 0, sizeof(fa)) ;

  for (irun = 0 ; irun < runs->rows ;  ac_free (h1), irun++)
    {
      h1 = ac_new_handle () ;
      AC_OBJ run = ac_table_obj (runs, irun, 0, h1) ;
      tt = ac_tag_table (run, "Observed_strandedness_in_ns_mapping", h1) ;
      if (tt &&  tt->rows && tt->cols >= 4)
	{
	  int ir = 0 ;
	  const char *ccp ;
	  
	  for (ir = 0 ; tt && ir < tt->rows ; ir++)
	    {
	      ccp = ac_table_printable (tt, ir, 0, 0) ;
	      if (ccp)
		{
		  dictAdd (dict, ccp, &target) ;
		  if (target < 15)  /* since we declared ft[16] */
		    {
		      fp[target] += ac_table_float (tt, ir, 2, 0) ;
		      fm[target] += ac_table_float (tt, ir, 4, 0) ;
		      fa[target] += ac_table_float (tt, ir, 6, 0) ;
		    }
		}
	    }
	}
    }

  vtxtPrintf (txt, "\nRun %s\n", ac_name (Group)) ;
  for (target = 1 ; target <= dictMax (dict) ; target++)
    {
      if (fp[target] + fm[target] > 0)
	{
	  vtxtPrintf (txt, "Observed_strandedness_in_ns_mapping %s ", dictName (dict, target)) ;
	  vtxtPrintf (txt, "%.2f %f plus %f minus %f ambiguous\n", 100.0 * fp[target]/(fp[target] + fm[target]), fp[target], fm[target], fa[target]) ;
	}
    }
  vtxtPrintf (txt, "\n") ;
  
 done:
  ac_free (h) ;
  return target ;
} /* baMergeObservedStrandedness */

/*************************************************************************************/
/* the strand is no longer used in the mapping but is needed for Wiggle and GeneExpression
 * set it automatically if not already done
 */
static int baSetStrand (BA *ba)
{
  AC_HANDLE h1 = 0,  h = ac_new_handle () ; 
  AC_DB db = ba->db ;
  AC_ITER iter = 0 ;
  AC_OBJ run = 0, ali = 0 ;
  AC_TABLE tbl = 0 ;
  int ir, nn = 0, nF = 0, nR = 0, nNS = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  const char *errors = 0, *qq ;
 
  qq = hprintf (h, "find run Is_run && ! strand %s %s && COUNT {>ali ; stranding} > 0"
		, ba->project ? " && project = " : ""
		, ba->project ? ba->project : ""
		) ;
  
  iter = ac_query_iter (ba->db, TRUE, qq, 0, h) ;
   if (iter)
     while (ac_free (run), ac_free (h1), run = ac_iter_obj (iter))
      {
	h1 = ac_new_handle () ;
	ali = ac_tag_obj (run, "Ali", h1) ;
	tbl = ac_tag_table (ali, "Stranding", h) ;
	for (ir = 0; ir < tbl->rows ; ir++)
	  {

	  }
      }
  /* parse the results */
   if (nn)
     {
       fprintf (stderr, " baSetStrand set %d Forward, %d Reverse, %d NS previously undecided strands:  %s\n"
		, nF, nR, nNS
		, timeShowNow()) ;
       fflush (stderr) ;
       if (vtxtPtr (txt))
	 ac_parse (db, vtxtPtr (txt), &errors, 0, h) ; 
       fprintf (stderr, " baSetStrand exit=%s %s\n", errors && *errors ? errors : "success", timeShowNow()) ;
     }
  
  ac_free (h1) ;
  ac_free (h) ;
  return nn ;
} /* baSetStrand */

/*************************************************************************************/

static int baGroupLetterProfileByLevel (BA *ba, int groupLevel)
{
  BOOL ok = FALSE ;
  AC_HANDLE h1 = 0,  h = ac_new_handle () ; 
  AC_DB db = ba->db ;
  int i, iu, jr, kr, n ;
  int badGroup, errPosMaxCol ;
  const char *errors = 0 ;
  AC_TABLE let ;
  AC_ITER iter1 = 0, iter2 = 0 ;
  AC_OBJ Run, Ali, Gr ;
  vTXT txt = vtxtHandleCreate (h) ;
  int Number_of_lanes, Maximal_read_multiplicity, iPerfect ;
  float raw_seq, raw_tags, raw_kb
    , acc_seq, acc_tags, acc_kb
    , rej_tmh_seq, rej_tmh_tags
    , rej_dfm_seq, rej_dfm_tags
    , rej_ni_seq, rej_ni_tags, rej_ni_kb
    , rej_unali_seq, rej_unali_tags, rej_unali_kb
    , rej_lqm_seq, rej_lqm_tags
    , rej_seq, rej_tags, rej_kb
    , adap1_tags, adap1_kb 
    , adap2_tags, adap2_kb
    , Aligned_fragments
    , Fragment_length_1
    , Fragment_length_5
    , Fragment_length_mode
    , Fragment_length_median
    , Fragment_length_average
    , Fragment_length_95
    , Fragment_length_99
    , intergenic, intergenic_density
    , accLnX, accLnK, accLnD, accLnC
    , partial3_15, partial5_15
    , partial3_25, partial5_25
    , perfect
    ;
  long int spots1, spots2, spots3 ;
  Array ks = 0, firstBase = 0, lastBase = 0 ;
  DICT *bloomDict = dictHandleCreate (100, h) ;
  DICT *spongeDict = dictHandleCreate (100, h) ;
  DICT *prefixDict = dictHandleCreate (100, h) ;
  DICT *targetDict = dictHandleCreate (100, h) ;
  DICT *errTypeDict = dictHandleCreate (100, h) ;
  DICT *compatibleDict = dictHandleCreate (100, h) ;
  DICT *orphanDict = dictHandleCreate (100, h) ;

  Array stranding_plus = arrayHandleCreate (100, float,h) ;
  Array stranding_minus = arrayHandleCreate (100, float,h) ;
  Array stranding_amb = arrayHandleCreate (100, float,h) ;
  Array hali_seq = arrayHandleCreate (100, float,h) ;
  Array hali_tag = arrayHandleCreate (100, float,h) ;
  Array hali_kb = arrayHandleCreate (100, float,h) ;
  Array hali_kb2 = arrayHandleCreate (100, float,h) ;
  
  Array nhali_seq = arrayHandleCreate (100, float,h) ;
  Array nhali_tag = arrayHandleCreate (100, float,h) ;
  Array nhali_kb = arrayHandleCreate (100, float,h) ;
  Array nhali_kb2 = arrayHandleCreate (100, float,h) ;

  Array cpu = arrayHandleCreate (100, float,h) ;
  Array elapsed = arrayHandleCreate (100, float,h) ;
  Array maxmem = arrayHandleCreate (100, float,h) ;

  Array lnDistrib = arrayHandleCreate (12, long int, h) ;

  Array errPos = arrayHandleCreate (100, float,h) ;
  Array NNPos = arrayHandleCreate (100, float,h) ;
  Array ErSpike = arrayHandleCreate (100, int,h) ;
  Array cumulMiss = arrayHandleCreate (8, float, h) ;
  Array cumulNN = arrayHandleCreate (8, float, h) ;
  Array errProf1 = arrayHandleCreate (1000, float, h) ;
  Array errProf2 = arrayHandleCreate (1000, float, h) ;
  Array clippedLn1 = arrayHandleCreate (1000, float, h) ;
  Array clippedLn2 = arrayHandleCreate (1000, float, h) ;
  Array clippedLn3 = arrayHandleCreate (1000, float, h) ;
  Array clippedLn4 = arrayHandleCreate (1000, float, h) ;
  Array clippedLn5 = arrayHandleCreate (1000, float, h) ;
  Array clippedLn6 = arrayHandleCreate (1000, float, h) ;
  Array clippedLn7 = arrayHandleCreate (1000, float, h) ;
  Array clippedLn8 = arrayHandleCreate (1000, float, h) ;
  Array aliLn1 = arrayHandleCreate (1000, float, h) ;
  Array aliLn2 = arrayHandleCreate (1000, float, h) ;
  Array aliLn3 = arrayHandleCreate (1000, float, h) ;
  Array aliLn4 = arrayHandleCreate (1000, float, h) ;
  Array aliLn5 = arrayHandleCreate (1000, float, h) ;
  Array aliLn6 = arrayHandleCreate (1000, float, h) ;
  Array aliLn7 = arrayHandleCreate (1000, float, h) ;
  Array aliLn8 = arrayHandleCreate (1000, float, h) ;
  Array pcaliLn1 = arrayHandleCreate (1000, float, h) ;
  Array pcaliLn2 = arrayHandleCreate (1000, float, h) ;
  Array pcaliLn3 = arrayHandleCreate (1000, float, h) ;
  Array pcaliLn4 = arrayHandleCreate (1000, float, h) ;
  Array pcaliLn5 = arrayHandleCreate (1000, float, h) ;
  Array pcaliLn6 = arrayHandleCreate (1000, float, h) ;
  Array pcaliLn7 = arrayHandleCreate (1000, float, h) ;
  Array pcaliLn8 = arrayHandleCreate (1000, float, h) ;
  Array pcmm1 = arrayHandleCreate (1000, float, h) ;
  Array pcmm2 = arrayHandleCreate (1000, float, h) ;
  Array pcmm3 = arrayHandleCreate (1000, float, h) ;
  Array pcmm4 = arrayHandleCreate (1000, float, h) ;
  Array pcmm5 = arrayHandleCreate (1000, float, h) ;
  Array pcmm6 = arrayHandleCreate (1000, float, h) ;
  Array pcmm7 = arrayHandleCreate (1000, float, h) ;
  Array pcmm8 = arrayHandleCreate (1000, float, h) ;
  Array pcNN1 = arrayHandleCreate (1000, float, h) ;
  Array pcNN2 = arrayHandleCreate (1000, float, h) ;
  Array pcNN3 = arrayHandleCreate (1000, float, h) ;
  Array pcNN4 = arrayHandleCreate (1000, float, h) ;
  Array pcNN5 = arrayHandleCreate (1000, float, h) ;
  Array pcNN6 = arrayHandleCreate (1000, float, h) ;
  Array pcNN7 = arrayHandleCreate (1000, float, h) ;
  Array pcNN8 = arrayHandleCreate (1000, float, h) ;
  Array countmm1 = arrayHandleCreate (1000, float, h) ;
  Array countmm2 = arrayHandleCreate (1000, float, h) ;
  Array countmm3 = arrayHandleCreate (1000, float, h) ;
  Array countmm4 = arrayHandleCreate (1000, float, h) ;
  Array countmm5 = arrayHandleCreate (1000, float, h) ;
  Array countmm6 = arrayHandleCreate (1000, float, h) ;
  Array countmm7 = arrayHandleCreate (1000, float, h) ;
  Array countmm8 = arrayHandleCreate (1000, float, h) ;
  Array countNN1 = arrayHandleCreate (1000, float, h) ;
  Array countNN2 = arrayHandleCreate (1000, float, h) ;
  Array countNN3 = arrayHandleCreate (1000, float, h) ;
  Array countNN4 = arrayHandleCreate (1000, float, h) ;
  Array countNN5 = arrayHandleCreate (1000, float, h) ;
  Array countNN6 = arrayHandleCreate (1000, float, h) ;
  Array countNN7 = arrayHandleCreate (1000, float, h) ;
  Array countNN8 = arrayHandleCreate (1000, float, h) ;
  Array amv_seq = arrayHandleCreate (1000, Array, h) ;
  Array amv = 0 ;
  Array amv_tag = arrayHandleCreate (1000, Array, h) ;
  Array amv_kb = arrayHandleCreate (1000, Array, h) ;
  Array amv_bp = arrayHandleCreate (1000, Array, h) ;
  Array unic[13] ; 
  AC_HANDLE amv_h = 0 ;
  Array compatibleArray =  arrayHandleCreate (50, float, h) ;
  Array orphanArray =  arrayHandleCreate (50, float, h) ;
  char *myquery1, *myquery2 ;
  Array bloomArray =  arrayHandleCreate (50, float, h) ;
  Array spongeArray =  arrayHandleCreate (50, float, h) ;

  dictAdd (errTypeDict, "Any", 0) ;
  for (iu = 0 ; iu < 13 ; iu++)
    unic[iu] = arrayHandleCreate (20, float, h) ;
  
  ks = arrayHandleCreate (1000, float, h) ;
  firstBase = arrayHandleCreate (1000, float, h) ;
  lastBase = arrayHandleCreate (1000, float, h) ;

  /* if level == 0 look for runs containing sublibraries
  *  if level > 0 look recursivelly for each group level
  */
  if (groupLevel == 0)
    myquery1 = messprintf ("find run Sublibraries  %s%s%s"
			   , ba->project ? " AND Project = \"" : ""
			   , ba->project ? ba->project : ""
			   , ba->project ? "\"" : ""
			   ) ;
  else
    myquery1 = messprintf ("find run group_level == %d %s%s%s"
			   , groupLevel
			   , ba->project ? " AND Project == \"" : ""
			   , ba->project ? ba->project : ""
			   , ba->project ? "\"" : ""
			   ) ;

  ac_db_refresh (db) ;
  iter1 = ac_dbquery_iter (db, myquery1, h) ;
  Gr = 0 ;
  while (ac_free (amv_h), ac_free (Gr), (Gr = ac_next_obj (iter1)))
    {
      ok = TRUE ; /* found at least one group at desired group_level */
      amv_h = ac_new_handle () ;

      perfect = iPerfect = 0 ;
      spots1 = spots2 = spots3 = 0 ;
      Aligned_fragments = 0 ;
      Fragment_length_1 = 0 ;
      Fragment_length_5 = 0 ;
      Fragment_length_mode = 0 ;
      Fragment_length_median = 0 ;
      Fragment_length_average = 0 ;
      Fragment_length_95 = 0 ;
      Fragment_length_99 = 0 ;
      Maximal_read_multiplicity = 0 ;

      accLnX = accLnK = accLnD = accLnC = 0 ;

      badGroup = 0 ;
      errPosMaxCol = 0 ;
      dictDestroy (prefixDict) ;
      prefixDict = dictHandleCreate (100, h) ;


      ks = arrayReCreate (ks, 1000, float) ;
      firstBase = arrayReCreate (firstBase, 1000, float) ;
      lastBase = arrayReCreate (lastBase, 1000, float) ;
      raw_seq = raw_tags = raw_kb = acc_seq = acc_tags = acc_kb = rej_seq = rej_tags = rej_kb = rej_tmh_seq = rej_tmh_tags = rej_lqm_seq = rej_lqm_tags = rej_dfm_seq = rej_dfm_tags = rej_ni_seq = rej_ni_tags = rej_ni_kb = rej_unali_seq = rej_unali_tags = rej_unali_kb = adap1_tags = adap1_kb  = adap2_tags = adap2_kb = Number_of_lanes= 0 ;
      intergenic = intergenic_density = 0 ;
      partial3_15 = partial5_15 = 0 ;
      partial3_25 = partial5_25 = 0 ;
      stranding_plus = arrayReCreate (stranding_plus,100,float) ;
      stranding_minus = arrayReCreate (stranding_minus,100,float) ;
      stranding_amb = arrayReCreate (stranding_amb,100,float) ;

      hali_seq = arrayReCreate (hali_seq,100,float) ;
      hali_tag = arrayReCreate (hali_tag,100,float) ;
      hali_kb = arrayReCreate (hali_kb,100,float) ;
      hali_kb2 = arrayReCreate (hali_kb2,100,float) ;

      nhali_seq = arrayReCreate (nhali_seq,100,float) ;
      nhali_tag = arrayReCreate (nhali_tag,100,float) ;
      nhali_kb = arrayReCreate (nhali_kb,100,float) ;
      nhali_kb2 = arrayReCreate (nhali_kb2,100,float) ;

      cpu = arrayReCreate (cpu,100,float) ;
      elapsed = arrayReCreate (elapsed,100,float) ;
      maxmem = arrayReCreate (maxmem,100,float) ;

      errPos = arrayReCreate (errPos,100,float) ;
      NNPos = arrayReCreate (NNPos,100,float) ;
      ErSpike = arrayReCreate (ErSpike,100,int) ;
      cumulMiss = arrayReCreate (cumulMiss, 6, float) ;
      cumulNN = arrayReCreate (cumulNN, 6, float) ;

      lnDistrib = arrayReCreate (lnDistrib, 12, long int) ;

      errProf1 = arrayReCreate (errProf1, 1000, float) ;
      errProf2 = arrayReCreate (errProf2, 1000, float) ;
      clippedLn1 = arrayReCreate (clippedLn1, 1000, float) ;
      clippedLn2 = arrayReCreate (clippedLn2, 1000, float) ;
      clippedLn3 = arrayReCreate (clippedLn3, 1000, float) ;
      clippedLn4 = arrayReCreate (clippedLn4, 1000, float) ;
      clippedLn5 = arrayReCreate (clippedLn5, 1000, float) ;
      clippedLn6 = arrayReCreate (clippedLn6, 1000, float) ;
      clippedLn7 = arrayReCreate (clippedLn7, 1000, float) ;
      clippedLn8 = arrayReCreate (clippedLn8, 1000, float) ;

      aliLn1 = arrayReCreate (aliLn1, 1000, float) ;
      aliLn2 = arrayReCreate (aliLn2, 1000, float) ;
      aliLn3 = arrayReCreate (aliLn3, 1000, float) ;
      aliLn4 = arrayReCreate (aliLn4, 1000, float) ;
      aliLn5 = arrayReCreate (aliLn5, 1000, float) ;
      aliLn6 = arrayReCreate (aliLn6, 1000, float) ;
      aliLn7 = arrayReCreate (aliLn7, 1000, float) ;
      aliLn8 = arrayReCreate (aliLn8, 1000, float) ;
      pcaliLn1 = arrayReCreate (pcaliLn1, 1000, float) ;
      pcaliLn2 = arrayReCreate (pcaliLn2, 1000, float) ;
      pcaliLn3 = arrayReCreate (pcaliLn3, 1000, float) ;
      pcaliLn4 = arrayReCreate (pcaliLn4, 1000, float) ;
      pcaliLn5 = arrayReCreate (pcaliLn5, 1000, float) ;
      pcaliLn6 = arrayReCreate (pcaliLn6, 1000, float) ;
      pcaliLn7 = arrayReCreate (pcaliLn7, 1000, float) ;
      pcaliLn8 = arrayReCreate (pcaliLn8, 1000, float) ;

      pcmm1 = arrayReCreate (pcmm1, 1000, float) ;
      pcmm2 = arrayReCreate (pcmm2, 1000, float) ;
      pcmm3 = arrayReCreate (pcmm3, 1000, float) ;
      pcmm4 = arrayReCreate (pcmm4, 1000, float) ;
      pcmm5 = arrayReCreate (pcmm5, 1000, float) ;
      pcmm6 = arrayReCreate (pcmm6, 1000, float) ;
      pcmm7 = arrayReCreate (pcmm7, 1000, float) ;
      pcmm8 = arrayReCreate (pcmm8, 1000, float) ;
      pcNN1 = arrayReCreate (pcNN1, 1000, float) ;
      pcNN2 = arrayReCreate (pcNN2, 1000, float) ;
      pcNN3 = arrayReCreate (pcNN3, 1000, float) ;
      pcNN4 = arrayReCreate (pcNN4, 1000, float) ;
      pcNN5 = arrayReCreate (pcNN5, 1000, float) ;
      pcNN6 = arrayReCreate (pcNN6, 1000, float) ;
      pcNN7 = arrayReCreate (pcNN7, 1000, float) ;
      pcNN8 = arrayReCreate (pcNN8, 1000, float) ;

      countmm1 = arrayReCreate (countmm1, 1000, float) ;
      countmm2 = arrayReCreate (countmm2, 1000, float) ;
      countmm3 = arrayReCreate (countmm3, 1000, float) ;
      countmm4 = arrayReCreate (countmm4, 1000, float) ;
      countmm5 = arrayReCreate (countmm5, 1000, float) ;
      countmm6 = arrayReCreate (countmm6, 1000, float) ;
      countmm7 = arrayReCreate (countmm7, 1000, float) ;
      countmm8 = arrayReCreate (countmm8, 1000, float) ;

      countNN1 = arrayReCreate (countNN1, 1000, float) ;
      countNN2 = arrayReCreate (countNN2, 1000, float) ;
      countNN3 = arrayReCreate (countNN3, 1000, float) ;
      countNN4 = arrayReCreate (countNN4, 1000, float) ;
      countNN5 = arrayReCreate (countNN5, 1000, float) ;
      countNN6 = arrayReCreate (countNN6, 1000, float) ;
      countNN7 = arrayReCreate (countNN7, 1000, float) ;
      countNN8 = arrayReCreate (countNN8, 1000, float) ;

      amv_seq = arrayReCreate (amv_seq, 1000, Array) ;
      amv_tag = arrayReCreate (amv_tag, 1000, Array) ;
      amv_kb = arrayReCreate (amv_kb, 1000, Array) ;
      amv_bp = arrayReCreate (amv_bp, 1000, Array) ;

      bloomArray =  arrayReCreate (bloomArray, 30, float) ;
      spongeArray =  arrayReCreate (spongeArray, 30, float) ;
 
      compatibleArray = arrayReCreate (compatibleArray, 30, float) ;
      orphanArray = arrayReCreate (orphanArray, 30, float) ;

      for (iu = 0 ; iu < 13 ; iu++)
	unic[iu] = arrayReCreate (unic[iu], 20, float) ;
      /* this is not good enough, in case of intersection and minus we must navigate back to the runs */
      myquery2 = hprintf (h, "{> union_of} SETOR {> sublibraries} ; project == \"%s\" ", ba->project ? ba->project : "*") ;
      iter2 = ac_objquery_iter (Gr, myquery2, h) ;
      Run = 0 ; h1 = 0 ;
      while (ac_free (h1), ac_free (Run), (Run = ac_next_obj (iter2)))
	{
	  h1 = ac_new_handle () ; 
	  Ali = ac_tag_obj (Run, "Ali", h1) ;
	  kr = ac_tag_int (Ali, "Maximal_read_multiplicity", 0) ;
	  if (kr > Maximal_read_multiplicity)
	    Maximal_read_multiplicity = kr ;
	  let = ac_tag_table (Ali, "Letter_profile", h1) ;
	  for (kr = 0 ; let && kr < let->rows ; kr++)
	    {
	      int iPrefix ;

	      dictAdd (prefixDict, ac_table_printable (let, kr, 0, "toto"), &iPrefix) ;
	      n = ac_table_int (let, kr, 1, 0) ;
	      if (n)
		{
		  array (ks, 50*n + 10*iPrefix,float) = n ;
		  for (i = 0 ; i < 7 ; i++)
		    array (ks, 50*n+10*iPrefix+1+i,float) += ac_table_float (let, kr, i+7, 0) ;
		}
	    }
	  let = ac_tag_table (Ali, "First_base_aligned", h1) ;
	  for (kr = 0 ; let && kr < let->rows ; kr++)
	    {
	      n = ac_table_int (let, kr, 0, 0) ;
	      if (n)
		{
		  array (firstBase, 4*n,float) = n ;
		  for (i = 1 ; i <= 2 ; i++)
		    array (firstBase, 4*n+i,float) += ac_table_float (let, kr, i, 0) ;
		}
	    }
	  let = ac_tag_table (Ali, "Accessible_length", h1) ;
	  if (let && let->cols >= 4)
	    {
	      float c  = ac_table_int (let, 0, 3, 0) ;
	      accLnX += c * ac_table_int (let, 0, 0, 0) ;
	      accLnK += c * ac_table_int (let, 0, 1, 0) ;
	      accLnD += c * ac_table_int (let, 0, 2, 0) ;
	      accLnC += c ;
	    }

	  let = ac_tag_table (Ali, "Last_base_aligned", h1) ;
	  for (kr = 0 ; let && kr < let->rows ; kr++)
	    {
	      n = ac_table_int (let, kr, 0, 0) ;
	      if (n)
		{
		  array (lastBase, 4*n,float) = n ;
		  for (i = 1 ; i <= 2 ; i++)
		    array (lastBase, 4*n+i,float) += ac_table_float (let, kr, i, 0) ;
		}
	    }
	  Number_of_lanes += ac_tag_int (Ali, "Number_of_lanes", 1) ; /* any ali implies at least one lane */
	  if ((let = ac_tag_table (Ali, "Raw_data", h1)) || (let = ac_tag_table (Ali, "Accepted", h1)) )
	    { /* if raw_data is absent, we approximate it by the always existing accepted counts, otherwise the group ratio may become very wrong */
	      raw_seq += ac_table_float (let, 0, 0, 0) ;
	      raw_tags += ac_table_float (let, 0, 2, 0) ;
	      raw_kb += ac_table_float (let, 0, 4, 0) ;
	    }

	  n = 0 ;
	  if ((let = ac_tag_table (Ali, "Accepted", h1)))
	    {
	      acc_seq += ac_table_float (let, 0, 0, 0) ;
	      n = ac_table_float (let, 0, 2, 0) ;
	      acc_tags += n ;
	      acc_kb += ac_table_float (let, 0, 4, 0) ;
	    }
	  else
	    if (! ac_has_tag (Ali, "Raw_data") && ! ac_has_tag (Ali, "h_Ali"))
	      badGroup = 2 ;

	  if (n && (let = ac_tag_table (Ali, "Length_distribution_1_5_50_95_99_mode_av", h1)))   /* read length distribution */
	    {
	      int j ;	  
	      array (lnDistrib, 0, long int) += n ;
	      for (j = 0 ; let && j < let->cols ; j++)
		{
		  int n1 = ac_table_int (let, 0, j, 0) ;
		  array (lnDistrib, j+1, long int) += n * n1 ;
		}
	    }

	  if (1) /* perfect */
	    {
	       if ((let = ac_tag_table (Ali, "Perfect_reads", h1)))
		 {
		   perfect += ac_table_float (let, 0, 0, 0) ;
		   iPerfect += ac_table_int (let, 0, 1, 0) ;
		 }
	    }

	  if (1) /* spots */
	    {
	       if ((let = ac_tag_table (Run, "Spots", h1)))
		 {
		   long int z = 0 ;
		   const char *ccp ;

		   z = 0 ;
		   ccp = ac_table_printable (let, 0, 0, "0") ;
		   sscanf (ccp, "%ld", &z) ;
		   spots1 += z ;

		   z = 0 ;
		   ccp = ac_table_printable (let, 0, 2, "0") ;
		   sscanf (ccp, "%ld", &z) ;
		   spots2 += z ;

		   z = 0 ;
		   ccp = ac_table_printable (let, 0, 8, "0") ;
		   sscanf (ccp, "%ld", &z) ;
		   spots3 += z ;
		 }
	    }

	  if (1) /* pair fate */
	    {
	      float z ;
	      z = ac_tag_float (Ali, "Aligned_fragments", 0) ;
	      z = ac_tag_float (Ali, "Aligned_fragments", 0) ;
	      if (z == 0)
		z = ac_tag_float (Ali, "Aligned_pairs", 0) ; /* old tag name prior to 2014_08_21 */
	      if (z > 0)
		{
		  float z1 ;
		  Aligned_fragments += z ;
		  z1 = ac_tag_int (Ali, "Fragment_length_1", 0) ;
		  Fragment_length_1 += z * z1 ;
		  z1 = ac_tag_int (Ali, "Fragment_length_5", 0) ;
		  Fragment_length_5 += z * z1 ;
		  z1 = ac_tag_int (Ali, "Fragment_length_mode", 0) ;
		  Fragment_length_mode += z * z1 ;
		  z1 = ac_tag_int (Ali, "Fragment_length_median", 0) ;
		  Fragment_length_median += z * z1 ;
		  z1 = ac_tag_int (Ali, "Fragment_length_average", 0) ;
		  Fragment_length_average += z * z1 ;
		  z1 = ac_tag_int (Ali, "Fragment_length_95", 0) ;
		  Fragment_length_95 += z * z1 ;
		  z1 = ac_tag_int (Ali, "Fragment_length_99", 0) ;
		  Fragment_length_99 += z * z1 ;
		}
	      if ((let = ac_tag_table (Ali, "Compatible", h1)))
		{
		  for (kr = 0 ; let && kr < let->rows ; kr++)
		    {
		      int iPrefix ;
		      
		      dictAdd (compatibleDict, ac_table_printable (let, kr, 0, "toto"), &iPrefix) ;
		      z = ac_table_float (let, kr, 1, 0) ;
		      if (z)
			{
			  array (compatibleArray, iPrefix,float) += z ;
			}
		    }
		}
	      if ((let = ac_tag_table (Ali, "Bloom", h1)))
		{
		  for (jr = 0 ; jr < let->rows ; jr++)
		    {
		      dictAdd (bloomDict, ac_table_printable (let, jr, 0, "toto"), &n) ;
		      for (kr = 1 ; let && kr < let->cols ; kr++)
			{			  
			  z = ac_table_float (let, jr, kr, 0) ;
			  if (z)
			    {
			      array (bloomArray, 100 * n + kr,float) += z ;
			    }
			}
		      array (bloomArray, 100 * n + kr,float) += z ;
		    }
		}

	      if ((let = ac_tag_table (Ali, "Sponge", h1)))
		{
		  for (jr = 0 ; jr < let->rows ; jr++)
		    {
		      dictAdd (spongeDict, ac_table_printable (let, jr, 0, "toto"), &n) ;
		      z = ac_table_float (let, jr, 1, 0) ;
		      if (z)
			  array (spongeArray, n,float) += z ;
		    }
		}

	      if ((let = ac_tag_table (Ali, "Incompatible", h1)))  /* use same dict since the schema uses tags */
		{
		  for (kr = 0 ; let && kr < let->rows ; kr++)
		    {
		      int iPrefix ;
		      
		      dictAdd (compatibleDict, ac_table_printable (let, kr, 0, "toto"), &iPrefix) ;
		      z = ac_table_float (let, kr, 1, 0) ;
		      if (z)
			{
			  array (compatibleArray, iPrefix,float) += z ;
			}
		    }
		}
	      if ((let = ac_tag_table (Ali, "Orphans", h1)))  /* use s specialized dict since the schema uses text */
		{
		  for (kr = 0 ; let && kr < let->rows ; kr++)
		    {
		      int iPrefix ;
		      
		      dictAdd (orphanDict, ac_table_printable (let, kr, 0, "toto"), &iPrefix) ;
		      z = ac_table_float (let, kr, 1, 0) ;
		      if (z)
			{
			  array (orphanArray, iPrefix,float) += z ;
			}
		    }
		}
	    }
	      
	  if ((let = ac_tag_table (Ali, "Intergenic", h1)))
	    {
	      intergenic += ac_table_float (let, 0, 0, 0) ;
	    }
	  if ((let = ac_tag_table (Ali, "Intergenic_density", h1)))
	    {
	      intergenic_density += ac_table_float (let, 0, 0, 0) ;
	    }

	  if ((let = ac_tag_table (Ali, "Aligned_fragments", h1)))
	    {
	      rej_seq += ac_table_float (let, 0, 0, 0) ;
	      rej_tags += ac_table_float (let, 0, 2, 0) ;
	      rej_kb += ac_table_float (let, 0, 4, 0)/1000.0 ;
	    }
	  else if ((let = ac_tag_table (Ali, "Aligned_pairs", h1))) /* old tag name before 2014_08_21 */
	    {
	      rej_seq += ac_table_float (let, 0, 0, 0) ;
	      rej_tags += ac_table_float (let, 0, 2, 0) ;
	      rej_kb += ac_table_float (let, 0, 4, 0)/1000.0 ;
	    }
	  if ((let = ac_tag_table (Ali, "Rejected", h1)))
	    {
	      rej_seq += ac_table_float (let, 0, 0, 0) ;
	      rej_tags += ac_table_float (let, 0, 2, 0) ;
	      rej_kb += ac_table_float (let, 0, 4, 0)/1000.0 ;
	    }
	  if ((let = ac_tag_table (Ali, "At_least_10_sites", h1)))
	    {
	      rej_tmh_seq += ac_table_float (let, 0, 0, 0) ;
	      rej_tmh_tags += ac_table_float (let, 0, 2, 0) ;
	    }
	  if ((let = ac_tag_table (Ali, "Partial_5p", h1)))
	    {
	      partial5_15 += ac_table_float (let, 0, 0, 0) ;
	      partial5_25 += ac_table_float (let, 0, 1, 0) ;
	    }
	  if ((let = ac_tag_table (Ali, "Partial_3p", h1)))
	    {
	      partial3_15 += ac_table_float (let, 0, 0, 0) ;
	      partial3_25 += ac_table_float (let, 0, 1, 0) ;
	    }
	  if ((let = ac_tag_table (Ali, "Diffuse_mapping", h1)))
	    {
	      rej_dfm_seq += ac_table_float (let, 0, 0, 0) ;
	      rej_dfm_tags += ac_table_float (let, 0, 2, 0) ;
	    }
	  if ((let = ac_tag_table (Ali, "No_insert", h1)))
	    {
	      rej_ni_seq += ac_table_float (let, 0, 0, 0) ;
	      rej_ni_tags += ac_table_float (let, 0, 2, 0) ;
	      rej_ni_kb += ac_table_float (let, 0, 4, 0) ;
	    }
	  if ((let = ac_tag_table (Ali, "Unaligned", h1)))
	    {
	      rej_unali_seq += ac_table_float (let, 0, 0, 0) ;
	      rej_unali_tags += ac_table_float (let, 0, 2, 0) ;
	      rej_unali_kb += ac_table_float (let, 0, 4, 0) ;
	    }
	  if ((let = ac_tag_table (Ali, "Low_quality_mapping", h1)))
	    {
	      rej_lqm_seq += ac_table_float (let, 0, 0, 0) ;
	      rej_lqm_tags += ac_table_float (let, 0, 2, 0) ;
	    }
	  if ((let = ac_tag_table (Ali, "stranding", h1)))
	    {
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  dictAdd (targetDict, ac_table_printable (let, kr, 0, "toto"), &n) ;
		  array(stranding_plus,n,float) += ac_table_float (let, kr, 2, 0) ;
		  array(stranding_minus,n,float) += ac_table_float (let, kr, 4, 0) ;
		  array(stranding_amb,n,float) += ac_table_float (let, kr, 6, 0) ;
		}
	    }
	  if ((let = ac_tag_table (Ali, "h_Ali", h1)))
	    {
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  dictAdd (targetDict, ac_table_printable (let, kr, 0, "toto"), &n) ;
		  array(hali_seq,n,float) += ac_table_float (let, kr, 1, 0) ;
		  array(hali_tag,n,float) += ac_table_float (let, kr, 3, 0) ;
		  array(hali_kb,n,float) += ac_table_float (let, kr, 5, 0) ;
		  array(hali_kb2,n,float) += ac_table_float (let, kr, 9, 0) ;
		}
	    }
	  if ((let = ac_tag_table (Ali, "nh_Ali", h1)))
	    {
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  dictAdd (targetDict, ac_table_printable (let, kr, 0, "toto"), &n) ;
		  array(nhali_seq,n,float) += ac_table_float (let, kr, 1, 0) ;
		  array(nhali_tag,n,float) += ac_table_float (let, kr, 3, 0) ;
		  array(nhali_kb,n,float) += ac_table_float (let, kr, 5, 0) ;
		  array(nhali_kb2,n,float) += ac_table_float (let, kr, 9, 0) ;
		  amv = array(amv_seq,n,Array) ; if (! amv) amv = array(amv_seq,n,Array) = arrayHandleCreate (200, float, amv_h) ; array (amv, arrayMax(amv), float) = ac_table_float (let, kr, 1, 0) ;
		  amv = array(amv_tag,n,Array) ;  if (! amv) amv = array(amv_tag,n,Array) = arrayHandleCreate (200, float, amv_h) ; array (amv, arrayMax(amv), float) = ac_table_float (let, kr, 3, 0) ;
		  amv = array(amv_kb,n,Array) ;  if (! amv) amv = array(amv_kb,n,Array) = arrayHandleCreate (200, float, amv_h) ; array (amv, arrayMax(amv), float) = ac_table_float (let, kr, 5, 0) ;
		  amv = array(amv_bp,n,Array) ;  if (! amv) amv = array(amv_bp,n,Array) = arrayHandleCreate (200, float, amv_h) ; array (amv, arrayMax(amv), float) = ac_table_float (let, kr, 7, 0) ;
		}
	    }
	  if ((let = ac_tag_table (Ali, "Error_position", h1)))
	    {
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  n = ac_table_int (let, kr, 0, 0) ;
		  if (n)
		    {
		      for (i = 1 ; i < let->cols && i < 10 ; i++)
			array(errPos, 10*n+i,float) += ac_table_float (let, kr, i, 0) ;
		      if (let->cols > errPosMaxCol) 
			errPosMaxCol = let->cols ;
		    }
		}
	    }
	  if ((let = ac_tag_table (Ali, "Ambiguous_position", h1)))
	    {
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  n = ac_table_int (let, kr, 0, 0) ;
		  if (n)  
		    {
		      for (i = 1 ; i < let->cols && i < 10 ; i++)
			array(NNPos,10*n+i,float) += ac_table_float (let, kr, i, 0) ;
		    }
		}
	    }
	  if ((let = ac_tag_table (Ali, "Clipped_length", h1)))
	    {
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  n = ac_table_int (let, kr, 0, 0) ;
		  if (n)
		    {
		      array (clippedLn1, n, float) += ac_table_float (let, kr, 1, 0) ;
		      array (clippedLn2, n, float) += ac_table_float (let, kr, 2, 0) ;
		      array (clippedLn3, n, float) += ac_table_float (let, kr, 3, 0) ;
		      array (clippedLn4, n, float) += ac_table_float (let, kr, 4, 0) ;
		      array (clippedLn5, n, float) += ac_table_float (let, kr, 5, 0) ;
		      array (clippedLn6, n, float) += ac_table_float (let, kr, 6, 0) ;
		      array (clippedLn7, n, float) += ac_table_float (let, kr, 7, 0) ;
		      array (clippedLn8, n, float) += ac_table_float (let, kr, 8, 0) ;
		    }
		}
	    }
	  if ((let = ac_tag_table (Ali, "Aligned_length", h1)))
	    {
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  n = ac_table_int (let, kr, 0, 0) ;
		  if (n)
		    {
		      array (aliLn1, n, float) += ac_table_float (let, kr, 1, 0) ;
		      array (aliLn2, n, float) += ac_table_float (let, kr, 2, 0) ;
		      array (aliLn3, n, float) += ac_table_float (let, kr, 3, 0) ;
		      array (aliLn4, n, float) += ac_table_float (let, kr, 4, 0) ;
		      array (aliLn5, n, float) += ac_table_float (let, kr, 5, 0) ;
		      array (aliLn6, n, float) += ac_table_float (let, kr, 6, 0) ;
		      array (aliLn7, n, float) += ac_table_float (let, kr, 7, 0) ;
		      array (aliLn8, n, float) += ac_table_float (let, kr, 8, 0) ;
		    }
		}
	    }
	  if ((let = ac_tag_table (Ali, "Per_cent_aligned", h1)))
	    {
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  n = ac_table_int (let, kr, 0, -1) ;
		  if (n>=0)
		    {
		      array (pcaliLn1, n, float) += ac_table_float (let, kr, 1, 0) ;
		      array (pcaliLn2, n, float) += ac_table_float (let, kr, 2, 0) ;
		      array (pcaliLn3, n, float) += ac_table_float (let, kr, 3, 0) ;
		      array (pcaliLn4, n, float) += ac_table_float (let, kr, 4, 0) ;
		      array (pcaliLn5, n, float) += ac_table_float (let, kr, 5, 0) ;
		      array (pcaliLn6, n, float) += ac_table_float (let, kr, 6, 0) ;
		      array (pcaliLn7, n, float) += ac_table_float (let, kr, 7, 0) ;
		      array (pcaliLn8, n, float) += ac_table_float (let, kr, 8, 0) ;
		    }
		}
	    }
	  if ((let = ac_tag_table (Ali, "Per_cent_mismatch", h1)))
	    {
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  n = ac_table_int (let, kr, 0, 0) ;
		  if (n)
		    {
		      array (pcmm1, n, float) += ac_table_float (let, kr, 1, 0) ;
		      array (pcmm2, n, float) += ac_table_float (let, kr, 2, 0) ;
		      array (pcmm3, n, float) += ac_table_float (let, kr, 3, 0) ;
		      array (pcmm4, n, float) += ac_table_float (let, kr, 4, 0) ;
		      array (pcmm5, n, float) += ac_table_float (let, kr, 5, 0) ;
		      array (pcmm6, n, float) += ac_table_float (let, kr, 6, 0) ;
		      array (pcmm7, n, float) += ac_table_float (let, kr, 7, 0) ;
		      array (pcmm8, n, float) += ac_table_float (let, kr, 8, 0) ;
		    }
		}
	    }
	  if ((let = ac_tag_table (Ali, "Per_cent_ambiguous", h1)))
	    {
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  n = ac_table_int (let, kr, 0, -1) ;
		  if (n>=0)
		    {
		      array (pcNN1, n, float) += ac_table_float (let, kr, 1, 0) ;
		      array (pcNN2, n, float) += ac_table_float (let, kr, 2, 0) ;
		      array (pcNN3, n, float) += ac_table_float (let, kr, 3, 0) ;
		      array (pcNN4, n, float) += ac_table_float (let, kr, 4, 0) ;
		      array (pcNN5, n, float) += ac_table_float (let, kr, 5, 0) ;
		      array (pcNN6, n, float) += ac_table_float (let, kr, 6, 0) ;
		      array (pcNN7, n, float) += ac_table_float (let, kr, 7, 0) ;
		      array (pcNN8, n, float) += ac_table_float (let, kr, 8, 0) ;
		    }
		}
	    }
	  if ((let = ac_tag_table (Ali, "Count_mismatch", h1)))
	    {
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  n = ac_table_int (let, kr, 0, -1) ;
		  if (n>=0)
		    {
		      array (countmm1, n, float) += ac_table_float (let, kr, 1, 0) ;
		      array (countmm2, n, float) += ac_table_float (let, kr, 2, 0) ;
		      array (countmm3, n, float) += ac_table_float (let, kr, 3, 0) ;
		      array (countmm4, n, float) += ac_table_float (let, kr, 4, 0) ;
		      array (countmm5, n, float) += ac_table_float (let, kr, 5, 0) ;
		      array (countmm6, n, float) += ac_table_float (let, kr, 6, 0) ;
		      array (countmm7, n, float) += ac_table_float (let, kr, 7, 0) ;
		      array (countmm8, n, float) += ac_table_float (let, kr, 8, 0) ;
		    }
		}
	    }
	  if ((let = ac_tag_table (Ali, "Count_ambiguous", h1)))
	    {
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  n = ac_table_int (let, kr, 0, -1) ;
		  if (n>=0)
		    {
		      array (countNN1, n, float) += ac_table_float (let, kr, 1, 0) ;
		      array (countNN2, n, float) += ac_table_float (let, kr, 2, 0) ;
		      array (countNN3, n, float) += ac_table_float (let, kr, 3, 0) ;
		      array (countNN4, n, float) += ac_table_float (let, kr, 4, 0) ;
		      array (countNN5, n, float) += ac_table_float (let, kr, 5, 0) ;
		      array (countNN6, n, float) += ac_table_float (let, kr, 6, 0) ;
		      array (countNN7, n, float) += ac_table_float (let, kr, 7, 0) ;
		      array (countNN8, n, float) += ac_table_float (let, kr, 8, 0) ;
		    }
		}
	    }
	  if ((let = ac_tag_table (Ali, "Cumulated_mismatches", h1)))
	    {
	      for (kr = 0 ; let && kr < let->cols ; kr++)
		{
		  array (cumulMiss, kr, float) += ac_table_float (let, 0, kr, 0) ;
		}
	    }
	  if ((let = ac_tag_table (Ali, "Cumulated_ambiguous", h1)))
	    {
	      for (kr = 0 ; let && kr < let->cols ; kr++)
		{
		  array (cumulNN, kr, float) += ac_table_float (let, 0, kr, 0) ;
		}
	    }
	  if ((let = ac_tag_table (Ali, "Error_profile", h1)))
	    {
	      int iPrefix ;
	      
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  dictAdd (prefixDict, ac_table_printable (let, kr, 0, "toto"), &iPrefix) ;
		  dictAdd (errTypeDict, ac_table_printable (let, kr, 1, "toto"), &n) ;
		  array (errProf1, iPrefix + 10 * n, float) += ac_table_float (let, kr, 2, 0) ;
		  array (errProf2, iPrefix + 10 * n, float) += ac_table_float (let, kr, 3, 0) ;
		}
	    }
	  if ((let = ac_tag_table (Ali, "Errror_spike", h1)))
	    {
	      int x, iPrefix ;

	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  dictAdd (prefixDict, ac_table_printable (let, kr, 0, "toto"), &iPrefix) ;
		  dictAdd (errTypeDict, ac_table_printable (let, kr, 2, "toto"), &n) ; 
		  x = ac_table_int (let, kr, 1, 0) ;
		  if (n < 1000 && x > 0)
		    array (ErSpike, iPrefix + 10 * n + 1000 *x, int) += ac_table_int (let, kr, 3, 0) ;
		}
	    }
	  if ((let = ac_tag_table (Ali, "CPU", h1)))
	    {
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  dictAdd (targetDict, ac_table_printable (let, kr, 0, "toto"), &n) ;
		  array (cpu, n,float) += ac_table_int (let, kr, 1, 0) ;
		}
	    }
	  if ((let = ac_tag_table (Ali, "Elapsed", h1)))
	    {
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  dictAdd (targetDict, ac_table_printable (let, kr, 0, "toto"), &n) ;
		  array (elapsed, n,float) += ac_table_int (let, kr, 1, 0) ;
		}
	    }
	  if ((let = ac_tag_table (Ali, "Max_memory", h1)))
	    {
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  dictAdd (targetDict, ac_table_printable (let, kr, 0, "toto"), &n) ;
		  if (array (maxmem, n,float) < ac_table_int (let, kr, 1, 0))
		    array (maxmem, n,float) = ac_table_int (let, kr, 1, 0) ;
		}
	    }
	  if ((let = ac_tag_table (Ali, "Unicity", h1)))
	    {
	      const char *cpp2 ;
	      for (kr = 0 ; let && kr < let->rows ; kr++)
		{
		  cpp2 = ac_table_printable (let, kr, 0, "Number_of_targets") ;
		  if (! strcmp (cpp2, "Number_of_targets")) continue ;
		  dictAdd (targetDict, cpp2, &n) ;
		  for (iu = 0 ; iu < 13 ; iu++)
		    array (unic[iu], n, float) += ac_table_float (let, kr, iu, 0) ;
		}
	    }
	}
      vtxtPrintf (txt, "Run %s\n-D Ali\n", ac_name(Gr)) ;
      if (spots1 > 0 && badGroup < 2)
	{
	  vtxtPrintf (txt, "Spots %ld bases_in_SRA %ld Average_length %ld Insert_size 0 spots_with_mate %ld\n", spots1, spots2, spots2/spots1, spots3) ;
	}
      vtxtPrintf (txt, "\nAli %s\nRun %s\n-D Rejected_alignments\n-D Rejected_too_many_hits\n-D Counts\n-D Pair_fate\n-D Strandedness\n-D Ali\n-D Unicity\n-D Alignments\n-D Letter_profile\n-D ATGC_kb\n-D Errors\n-D Computer_ressource"
		  , ac_name(Gr)
		  , ac_name(Gr)
		  ) ;
      
      if (badGroup < 2)
	{
	  if (intergenic > 0)
	    vtxtPrintf (txt, "\nIntergenic %.1f kb"
			, intergenic) ;

	  if (intergenic_density > 0)
	    vtxtPrintf (txt, "\nIntergenic_density %.4f"
			, intergenic_density) ;


	  vtxtPrintf (txt, "\nPerfect_reads %f %d", perfect, iPerfect) ;

	  vtxtPrintf (txt, "\nNumber_of_lanes %d", Number_of_lanes) ;
	  vtxtPrintf (txt, "\nRaw_data %f Id %f Accepted %f kb  %d bp_per_tag"
		      , raw_seq, raw_tags, raw_kb, raw_kb > 1 ? (int)(.5 + (1000.0 * raw_kb)/raw_tags) : 0);
	  vtxtPrintf (txt, "\nAccepted %f Seq %f Tags %f kb %d bp_per_tag"
		      , acc_seq, acc_tags, acc_kb, acc_kb > 1 ? (int)(.5 + (1000.0 * acc_kb)/acc_tags): 0) ;
	  vtxtPrintf (txt, "\nRejected %f NA %f Tags %f \"kb clipped or rejected\""
		      , 0*rej_seq, rej_tags, rej_kb
		      /* , rej_kb > 1 ? (int)(.5 + (1000.0 * rej_kb)/rej_tags) : 0 */
		      ) ;
	  vtxtPrintf (txt, "\nAt_least_10_sites %f Seq %f Tags"
		      , rej_tmh_seq, rej_tmh_tags) ;
	  vtxtPrintf (txt, "\nDiffuse_mapping %f Seq %f Tags"
		      , rej_dfm_seq, rej_dfm_tags) ;
	  vtxtPrintf (txt, "\nNo_insert %f Seq %f Tags %f kb"
		      , rej_ni_seq, rej_ni_tags, rej_ni_kb) ;
	  vtxtPrintf (txt, "\nUnaligned %f Seq %f Tags %f kb"
		      , rej_unali_seq, rej_unali_tags, rej_unali_kb) ;
	  vtxtPrintf (txt, "\nLow_quality_mapping %f Seq %f Tags"
		      , rej_lqm_seq, rej_lqm_tags) ;
	  vtxtPrintf (txt, "\nPartial_5p %f %f", partial5_15, partial5_25) ;
	  vtxtPrintf (txt, "\nPartial_3p %f %f", partial3_15, partial3_25) ;

	  for (n = 1 ; n <= arrayMax (errPos)/10 ; n++)
	    {
	      vtxtPrintf (txt, "\nError_position %d", n) ;
	      for (i = 1 ; i < errPosMaxCol && i < 10 ; i++)
		vtxtPrintf (txt, " %g", array (errPos, 10*n + i,float) ) ;
	    }
	  for (n = 1 ; n <= arrayMax (NNPos)/10 ; n++)
	    {
	      vtxtPrintf (txt, "\nAmbiguous_position %d", n) ;
	      for (i = 1 ; i < errPosMaxCol && i < 10 ; i++)
		vtxtPrintf (txt, " %g", array (NNPos, 10*n + i,float) ) ;
	    }

	  vtxtPrintf (txt, "\nCumulated_mismatches ") ;
	  for (n = 0 ; n < 8 ; n++)
	    vtxtPrintf (txt, " %g", array(cumulMiss, n, float)) ;
	  vtxtPrintf (txt, "\nCumulated_ambiguous ") ;
	  for (n = 0 ; n < 8 ; n++)
	    vtxtPrintf (txt, " %g", array(cumulNN, n, float)) ;
	  for (n = 0 ; n < arrayMax (stranding_plus) ; n++)
	    {
	      double x  = array (stranding_plus, n,float) + array (stranding_minus, n,float) ; /* avoid the ambiguous */
	      if (x >= 1)
		vtxtPrintf (txt, "\nstranding %s %g %g plus %g minus %g ambiguous"
			    , dictName (targetDict, n)
			    ,   (100.0 * array (stranding_plus, n,float)/x)
			    , (double) array (stranding_plus, n,float)
			    , (double) array (stranding_minus, n,float)
			    , (double) array (stranding_amb, n,float)
			    ) ;
	    }
	  if (dictFind (errTypeDict, "Any", &n))
	    {
	      int n0 = n, iPrefix ;
	     
	      for (iPrefix = 1 ; iPrefix <= dictMax (prefixDict) ; iPrefix++)
		for (n = 1 ; n <= dictMax (errTypeDict) && iPrefix + 10 * n < arrayMax (errProf1) ; n++)
		  {
		    float x, y, z ;
		    x  = arr (errProf1, iPrefix + 10 * n, float) ;
		    y  = arr (errProf2, iPrefix + 10 * n0, float) ;  /* alway use the total count */
		    z = y > 1 ? x/y : 0 ;
		    if (x >= 1 || n == n0)   /* alway use the total count */
		      vtxtPrintf (txt, "\nError_profile %s  %s %g %g %f"
				  , dictName (prefixDict, iPrefix)
				  , dictName (errTypeDict, n)
				  , x, y, z
				  ) ;
		  }
	      /* ErSpike goes here */
	      for (iPrefix = 1 ; 0 && iPrefix <= dictMax (prefixDict) ; iPrefix++)
		for (n = 1 ; n <= dictMax (errTypeDict) && iPrefix + 10 * n < arrayMax (errProf1) ; n++)
		  {
		    float x, y, z ;
		    x  = arr (errProf1, iPrefix + 10 * n, float) ;
		    y  = arr (errProf2, iPrefix + 10 * n0, float) ;  /* alway use the total count */
		    z = y > 1 ? x/y : 0 ;
		    if (x >= 1 || n == n0)   /* alway use the total count */
		      vtxtPrintf (txt, "\nError_profile %s  %s %g %g %f"
				  , dictName (prefixDict, iPrefix)
				  , dictName (errTypeDict, n)
				  , x, y, z
				  ) ;
		  }
	    }
	  for (n = 0 ; n < arrayMax (aliLn1) ; n++)
	    {
	      float x1, x2, x3, x4, x5, x6, x7, x8 ;
	      x1  = arr (aliLn1, n, float) ;
	      x2  = arr (aliLn2, n, float) ;
	      x3  = arr (aliLn3, n, float) ;
	      x4  = arr (aliLn4, n, float) ;
	      x5  = arr (aliLn5, n, float) ;
	      x6  = arr (aliLn6, n, float) ;
	      x7  = arr (aliLn7, n, float) ;
	      x8  = arr (aliLn8, n, float) ;
	      if (x1 >= 1)
		vtxtPrintf (txt, "\nAligned_length %d %f %f %f %f %f %f %f %f"
			    , n
			    , x1, x2, x3, x4, x5, x6, x7, x8
			    ) ;
	    }
	  for (n = 0 ; n < arrayMax (pcaliLn1) ; n++)
	    {
	      float x1, x2, x3, x4, x5, x6, x7, x8 ;
	      x1  = arr (pcaliLn1, n, float) ;
	      x2  = arr (pcaliLn2, n, float) ;
	      x3  = arr (pcaliLn3, n, float) ;
	      x4  = arr (pcaliLn4, n, float) ;
	      x5  = arr (pcaliLn5, n, float) ;
	      x6  = arr (pcaliLn6, n, float) ;
	      x7  = arr (pcaliLn7, n, float) ;
	      x8  = arr (pcaliLn8, n, float) ;
	      if (x1 >= 1)
		vtxtPrintf (txt, "\nPer_cent_aligned %d %f %f %f %f %f %f %f %f"
			    , n
			    , x1, x2, x3, x4, x5, x6, x7, x8
			    ) ;
	    }
	  for (n = 0 ; n < arrayMax (pcmm1) ; n++)
	    {
	      float x1, x2, x3, x4, x5, x6, x7, x8 ;
	      x1  = arr (pcmm1, n, float) ;
	      x2  = arr (pcmm2, n, float) ;
	      x3  = arr (pcmm3, n, float) ;
	      x4  = arr (pcmm4, n, float) ;
	      x5  = arr (pcmm5, n, float) ;
	      x6  = arr (pcmm6, n, float) ;
	      x7  = arr (pcmm7, n, float) ;
	      x8  = arr (pcmm8, n, float) ;
	      if (x1 >= 0)
		vtxtPrintf (txt, "\nPer_cent_mismatch %d %f %f %f %f %f %f %f %f"
			    , n
			    , x1, x2, x3, x4, x5, x6, x7, x8
			    ) ;
	    }
	  for (n = 0 ; n < arrayMax (pcNN1) ; n++)
	    {
	      float x1, x2, x3, x4, x5, x6, x7, x8 ;
	      x1  = arr (pcNN1, n, float) ;
	      x2  = arr (pcNN2, n, float) ;
	      x3  = arr (pcNN3, n, float) ;
	      x4  = arr (pcNN4, n, float) ;
	      x5  = arr (pcNN5, n, float) ;
	      x6  = arr (pcNN6, n, float) ;
	      x7  = arr (pcNN7, n, float) ;
	      x8  = arr (pcNN8, n, float) ;
	      if (x1 >= 0)
		vtxtPrintf (txt, "\nPer_cent_ambiguous %d %f %f %f %f %f %f %f %f"
			    , n
			    , x1, x2, x3, x4, x5, x6, x7, x8
			    ) ;
	    }
	  for (n = 0 ; n < arrayMax (countmm1) ; n++)
	    {
	      float x1, x2, x3, x4, x5, x6, x7, x8 ;
	      x1  = arr (countmm1, n, float) ;
	      x2  = arr (countmm2, n, float) ;
	      x3  = arr (countmm3, n, float) ;
	      x4  = arr (countmm4, n, float) ;
	      x5  = arr (countmm5, n, float) ;
	      x6  = arr (countmm6, n, float) ;
	      x7  = arr (countmm7, n, float) ;
	      x8  = arr (countmm8, n, float) ;
	      if (x1 >= 0)
		vtxtPrintf (txt, "\nCount_mismatch %d %f %f %f %f %f %f %f %f"
			    , n
			    , x1, x2, x3, x4, x5, x6, x7, x8
			    ) ;
	    }
	  for (n = 0 ; n < arrayMax (countNN1) ; n++)
	    {
	      float x1, x2, x3, x4, x5, x6, x7, x8 ;
	      x1  = arr (countNN1, n, float) ;
	      x2  = arr (countNN2, n, float) ;
	      x3  = arr (countNN3, n, float) ;
	      x4  = arr (countNN4, n, float) ;
	      x5  = arr (countNN5, n, float) ;
	      x6  = arr (countNN6, n, float) ;
	      x7  = arr (countNN7, n, float) ;
	      x8  = arr (countNN8, n, float) ;
	      if (x1 >= 0)
		vtxtPrintf (txt, "\nCount_ambiguous %d %f %f %f %f %f %f %f %f"
			    , n
			    , x1, x2, x3, x4, x5, x6, x7, x8
			    ) ;
	    }

	  if (arrayMax (lnDistrib))    /* read length distribution */
	    {
	      long int nn = array (lnDistrib, 0, long int) ;
	      if (nn)
		{
		  int i ;
		  vtxtPrint (txt, "\nLength_distribution_1_5_50_95_99_mode_av ") ;
		  for (i = 1 ; i < arrayMax (lnDistrib) ; i++)
		    vtxtPrintf (txt, " %ld ", array (lnDistrib, i, long int)/nn) ;
		}
	    }

	  for (n = 0 ; n < arrayMax (clippedLn1) ; n++)
	    {
	      float x1, x2, x3, x4, x5, x6, x7, x8 ;
	      x1  = arr (clippedLn1, n, float) ;
	      x2  = arr (clippedLn2, n, float) ;
	      x3  = arr (clippedLn3, n, float) ;
	      x4  = arr (clippedLn4, n, float) ;
	      x5  = arr (clippedLn5, n, float) ;
	      x6  = arr (clippedLn6, n, float) ;
	      x7  = arr (clippedLn7, n, float) ;
	      x8  = arr (clippedLn8, n, float) ;
	      if (x1 >= 1)
		vtxtPrintf (txt, "\nClipped_length %d %f %f %f %f %f %f %f %f"
			    , n
			    , x1, x2, x3, x4, x5, x6, x7, x8
			    ) ;
	    }
	  
	  for (n = 0 ; n < arrayMax (hali_tag) ; n++)
	    {
	      double x  =  array (hali_tag, n,float) ;
	      if (x >= 1)
		vtxtPrintf (txt, "\nh_Ali %s %g seq %g tag %g kb_aligned %.2f bp %g kb_to_be_aligned %.2f bp_to_be_aligned"
			    , dictName (targetDict, n)
			    , (double) array (hali_seq, n,float)
			    , (double) array (hali_tag, n,float)
			    , (double) array (hali_kb, n,float)
			    , (1000.0 * array (hali_kb, n,float)/array (hali_tag, n,float))
			    , (double) array (hali_kb2, n,float)
			    , (1000.0 * array (hali_kb2, n,float)/array (hali_tag, n,float))
			    ) ;
	    } 
	  for (n = 0 ; n < arrayMax (amv_seq) ; n++)
	    {
	      float z0, z1 ;
	      
	      amv = array(amv_seq,n,Array) ;
	      
	      if (amv)
		{
		  vtxtPrintf (txt, "\nnh_median  %s ", dictName (targetDict, n)) ;
		  
		  amv = array(amv_seq,n,Array) ;
		  floatVariance (amv, &z0, &z1, 0) ;
		  vtxtPrintf (txt, " %f seq ", z0) ;
		  
		  amv = array(amv_tag,n,Array) ;
		  floatVariance (amv, &z0, &z1, 0) ;
		  vtxtPrintf (txt, " %f tag ", z0) ;
		  
		  amv = array(amv_kb,n,Array) ;
		  floatVariance (amv, &z0, &z1, 0) ;
		  vtxtPrintf (txt, " %f kb ", z0) ;
		  
		  amv = array(amv_bp,n,Array) ;
		  floatVariance (amv, &z0, &z1, 0) ;
		  vtxtPrintf (txt, " %f bp ", z0) ;
		} 
	    }
	  for (n = 0 ; n < arrayMax (amv_seq) ; n++)
	    {
	      float z1 ;
	      
	      amv = array(amv_seq,n,Array) ;
	      if (amv)
		{
		  vtxtPrintf (txt, "\nnh_average  %s ", dictName (targetDict, n)) ;
		  
		  amv = array(amv_seq,n,Array) ;
		  floatVariance (amv, 0, &z1, 0) ;
		  vtxtPrintf (txt, " %f seq ", z1) ;
		  
		  amv = array(amv_tag,n,Array) ;
		  floatVariance (amv, 0, &z1, 0) ;
		  vtxtPrintf (txt, " %f tag ", z1) ;
		  
		  amv = array(amv_kb,n,Array) ;
		  floatVariance (amv, 0, &z1, 0) ;
		  vtxtPrintf (txt, " %f kb ", z1) ;
		  
		  amv = array(amv_bp,n,Array) ;
		  floatVariance (amv, 0, &z1, 0) ;
		  vtxtPrintf (txt, " %f bp ", z1) ;
		} 
	    }
	  
	  for (n = 0 ; n < arrayMax (amv_seq) ; n++)
	    {
	      float z1, z2 ;
	      
	      amv = array(amv_seq,n,Array) ;
	      if (amv)
		{
		  vtxtPrintf (txt, "\nnh_sigma  %s ", dictName (targetDict, n)) ;
		  
		  amv = array(amv_seq,n,Array) ;
		  floatVariance (amv, 0, &z1, &z2) ;
		  vtxtPrintf (txt, " %f seq ", z2) ;
		  
		  amv = array(amv_tag,n,Array) ;
		  floatVariance (amv, 0, &z1, &z2) ;
		  vtxtPrintf (txt, " %f tag ", z2) ;
		  
		  amv = array(amv_kb,n,Array) ;
		  floatVariance (amv, 0, &z1, &z2) ;
		  vtxtPrintf (txt, " %f kb ", z2) ;
		  
		  amv = array(amv_bp,n,Array) ;
		  floatVariance (amv, 0, &z1, &z2) ;
		  vtxtPrintf (txt, " %f bp ", z2) ;
		} 
	    }
	  for (n = 0 ; n < arrayMax (cpu) ; n++)
	    {
	      int x  =  array (cpu, n,float) ;
	      if (x >= 1)
		vtxtPrintf (txt, "\nCPU %s %d seconds", dictName (targetDict, n), (int)array (cpu, n,float)) ; 
	    } 
	  for (n = 0 ; n < arrayMax (elapsed) ; n++)
	    {
	      int x  =  array (elapsed, n,float) ;
	      if (x >= 1)
		vtxtPrintf (txt, "\nElapsed %s %d seconds", dictName (targetDict, n), (int)array (elapsed, n,float)) ; 
	    } 
	  for (n = 0 ; n < arrayMax (maxmem) ; n++)
	    {
	      int x  =  array (maxmem, n,float) ;
	      if (x >= 1)
		vtxtPrintf (txt, "\nMax_memory %s %d Mb", dictName (targetDict, n), (int)array (maxmem, n,float)) ;
	    } 
	  vtxtPrintf (txt, "\nUnicity Number_of_targets 1 2 3 4 5 6 7 8 9 10 -2 -3") ;
	  for (n = 0 ; n < arrayMax (unic[1]) ; n++)
	    {
	      double x  =  array (unic[1], n, float) ;
	      if (x >= 1)
		{
		  if (strcmp (dictName (targetDict, n), "Target")) /* bug in previous code layer, dismiss */
		    {
		      vtxtPrintf (txt, "\nUnicity %s", dictName (targetDict, n)) ;
		      for (iu = 1 ; iu < 13 ; iu++)
			vtxtPrintf (txt, "  %g", array (unic[iu], n, float)) ;
		    }
		}
	    }
	  for (n = 0 ; n < arrayMax (nhali_tag) ; n++)
	    {
	      double x  =  array (nhali_tag, n,float) ;
	      if (x >= 1)
		vtxtPrintf (txt, "\nnh_Ali %s %g seq %g tag %g kb %.2f bp %g kb_to_be_aligned %.2f bp_to_be_aligned"
			    , dictName (targetDict, n)
			    , (double) array (nhali_seq, n,float)
			    , (double) array (nhali_tag, n,float)
			    , (double) array (nhali_kb, n,float)
			    , ( 1000.0 * array (nhali_kb, n,float))/array (nhali_tag, n,float)
			    , (double) array (nhali_kb2, n,float)
			    , (1000.0 * array (nhali_kb2, n,float))/array (nhali_tag, n,float)
			    ) ;
	    }
	  if (arrayMax (ks))
	    {
	      int iPrefix ;
	      long int nnn = 0, nn, x[8], xx ;
	      for (i = 0 ; i < 8 ; i++)
		x[i] = 0 ;
	      
	      for (iPrefix = 1 ; iPrefix <= dictMax (prefixDict) ; iPrefix++)
		{
		  for (kr = 1 ; 50*kr+10*iPrefix+1+5 < arrayMax (ks) ; kr++)
		    {
		      vtxtPrintf (txt, "\nLetter_profile %s %d", dictName (prefixDict, iPrefix), kr) ; /* coordinate */
		      for (nn = 0, i = 1 ; i < 6 ; i++)
			{ nn += array (ks, 50*kr+10*iPrefix+1+i,float) ; x[i] +=  array (ks, 50*kr+10*iPrefix+1+i,float) ; }
		      /* export the proportions */
		      if (nn == 0)
			vtxtPrintf (txt, " 0 0 0 0 0 ") ;
		      else
			{
			  for (i = 1 ; i < 6 ; i++)
			    {
			      xx = 100.0 * array (ks, 50*kr+10*iPrefix+1+i,float) ;
			      xx /= nn ;
			      vtxtPrintf (txt, " %ld ", xx) ;
			    }
			}
		      /* export the counts */
		      for (i = 0 ; i < 7 ; i++)
			vtxtPrintf (txt, " %g", array (ks, 50*kr+10*iPrefix+1+i,float)) ;
		    }
		  for (nnn = 0, i = 1 ; i < 6 ; i++)
		    nnn += x[i] ;
		}	      
	
	      if (nnn)
		{
		  vtxtPrintf (txt, "\nATGC_kb ") ;
		  for (i = 1 ; i < 6 ; i++)
		    {
		      xx = 1000.0 * x[i] ;
		      xx /= nnn ;
		      vtxtPrintf (txt, " %ld ", xx) ;
		    }
		  /* vtxtPrintf (txt, " %ld", nnn/1000) ;  faux, on veut pas la somme 2012_01_11 */
		  for (i = 1 ; i < 6 ; i++)
		    vtxtPrintf (txt, " %ld", x[i]/1000) ;
		}
	    }
	
	  if (arrayMax (firstBase))
	    {
	      BOOL hasN2 = FALSE ;

	      for (kr = 1 ; ! hasN2 && kr <= arrayMax (firstBase)/4 ; kr++)
		if (array (firstBase, 4*kr + 2, float) > 0) hasN2 = TRUE ;
	      for (kr = 1 ; kr <= arrayMax (firstBase)/4 ; kr++)
		{
		  if (array (firstBase, 4*kr, float) < 1) continue ;
		  
		  vtxtPrintf (txt, "\nFirst_base_aligned %d", kr)  ;
		  vtxtPrintf (txt, " %g ", array (firstBase, 4*kr + 1, float)) ;
		  if (hasN2) 
		    vtxtPrintf (txt, " %g ", array (firstBase, 4*kr + 2, float)) ;
		}
	    }
	  if (arrayMax (lastBase))
	    {
	      BOOL hasN2 = FALSE ;

	      for (kr = 1 ; ! hasN2 && kr <= arrayMax (lastBase)/4 ; kr++)
		if (array (lastBase, 4*kr + 2, float) > 0) hasN2 = TRUE ;
	      for (kr = 1 ; kr <= arrayMax (lastBase)/4 ; kr++)
		{
		  if (array (lastBase, 4*kr, float) < 1) continue ;
		  
		  vtxtPrintf (txt, "\nLast_base_aligned %d", kr)  ;
		  vtxtPrintf (txt, " %g ", array (lastBase, 4*kr + 1, float)) ;
		  if (hasN2) 
		    vtxtPrintf (txt, " %g ", array (lastBase, 4*kr + 2, float)) ;
		}
	    }
	}


      if (Aligned_fragments > 0)
	{
	  vtxtPrintf (txt, "\nAligned_fragments %f", Aligned_fragments) ;
	  vtxtPrintf (txt, "\nFragment_length_1 %d", (int) (Fragment_length_1/Aligned_fragments)) ;
	  vtxtPrintf (txt, "\nFragment_length_5 %d", (int) (Fragment_length_5/Aligned_fragments)) ;
	  vtxtPrintf (txt, "\nFragment_length_mode %d", (int) (Fragment_length_mode/Aligned_fragments)) ;
	  vtxtPrintf (txt, "\nFragment_length_median %d", (int) (Fragment_length_median/Aligned_fragments)) ;
	  vtxtPrintf (txt, "\nFragment_length_average %d", (int) (Fragment_length_average/Aligned_fragments)) ;
	  vtxtPrintf (txt, "\nFragment_length_95 %d", (int) (Fragment_length_95/Aligned_fragments)) ;
	  vtxtPrintf (txt, "\nFragment_length_99 %d", (int) (Fragment_length_99/Aligned_fragments)) ;
	}
      for (kr = 0 ; kr < arrayMax (compatibleArray) ; kr++)
	{
	  float z = array (compatibleArray, kr,float) ;
	  if (z > 0)
	    vtxtPrintf (txt, "\n%s %f", dictName (compatibleDict, kr), z) ;
	}
      for (kr = 0 ; kr < arrayMax (orphanArray) ; kr++)
	{
	  float z = array (orphanArray, kr,float) ;
	  if (z > 0)
	    vtxtPrintf (txt, "\nOrphans %s %f", dictName (orphanDict, kr), z) ;
	}

      vtxtPrintf (txt, "\n-D Bloom") ;
      if (arrayMax (bloomArray))
	{
	  for (n = 1 ; n <= dictMax (bloomDict) ; n++)
	    {
	      /* look for highest point to export */
	      for (kr = 99 ; kr > 0 ; kr--)
		if (array (bloomArray, 100 * n + kr,float) > 0)
		  {
		    vtxtPrintf (txt, "\nBloom %s ",  dictName (bloomDict, n)) ;
		    for (jr = 1 ; jr <= kr ; jr++)
		      vtxtPrintf (txt, " %.0f", array (bloomArray, 100 * n + jr,float)) ;
		  }
	    }
	}

      if (arrayMax (spongeArray))
	{
	  for (n = 1 ; n <= dictMax (spongeDict) ; n++)
	    if (array (spongeArray, n,float) > 0)
	      vtxtPrintf (txt, "\nSponge %s %f \"Mb aligned\""
			  ,  dictName (spongeDict, n)
			  , array (spongeArray, n,float)
			  ) ;
	}

      if (accLnC > 0)
	{
	  int x = accLnX / accLnC ;
	  int K = accLnK / accLnC ;
	  int D = accLnD / accLnC ;
	  int C = accLnC ;

	  vtxtPrintf (txt, "\nAccessible_length %d %d %d %d", x, K, D, C) ; 
	}
      
      if (Maximal_read_multiplicity > 0) 
	vtxtPrintf (txt, "\nMaximal_read_multiplicity %d", Maximal_read_multiplicity) ;

      vtxtPrintf (txt, "\n\n") ;
      baMergeObservedStrandedness (txt, Gr, groupLevel) ;  /* the merging works for the runs */

    }
  
  /* parse the results */
  fprintf (stderr, " baGroupLetterProfileByLevel %d  parsing the results %s\n", groupLevel, timeShowNow()) ;
  fflush (stderr) ;
  if (vtxtPtr (txt))
    ac_parse (db, vtxtPtr (txt), &errors, 0, h) ; 
  fprintf (stderr, " baGroupLetterProfileByLevel %d exit=%s %s\n", groupLevel, errors && *errors ? errors : "success", timeShowNow()) ;
  
  ac_free (h1) ;
  ac_free (h) ;
  return ok ;
} /* baGroupLetterProfileByLevel */

/*************************************************************************************/
/*************************************************************************************/

static int  baAddQualityFactors (BA *ba)
{
  AC_HANDLE h = 0,  h1 = ac_new_handle () ; 
  AC_DB db = ba->db ;
  AC_ITER iter = 0 ;
  AC_TABLE LP, EP, FB, LB ;
  int ir, jr, x, xmax, nOk = 0 ;
  double a, z, q, L10 ;
  const char* ccp ;
  const char *errors = 0 ;
  char prefix[8] ;
  AC_OBJ Ali = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  Array ks1 = arrayHandleCreate (200, double, h) ;
  Array ks2 = arrayHandleCreate (200, double, h) ;

  L10 = log(10.0) ;

  fprintf (stderr, "Add quality_profile start: %s\n", timeShowNow()) ;
  vtxtPrint (txt, "Find run  ") ;
  if (ba->project)
    vtxtPrintf (txt, " project = \"%s\"", ba->project) ;
  vtxtPrint (txt, " ; FOLLOW ali ; First_base_aligned && Last_base_aligned && Letter_profile && Error_position") ;
  
  iter = ac_query_iter (db , TRUE, vtxtPtr (txt) , 0, h1) ;

  vtxtClear (txt) ; 
  while (ac_free (h), ac_free (Ali), Ali = ac_iter_obj (iter))
    {
      vtxtPrintf (txt, "\n\nAli \"%s\"\n-D Quality_profile", ac_name(Ali)) ;
      nOk++ ;

      FB = ac_tag_table (Ali, "First_base_aligned", h) ;
      LB = ac_tag_table (Ali, "Last_base_aligned", h) ;
      EP = ac_tag_table (Ali, "Error_position", h) ;
      LP = ac_tag_table (Ali, "Letter_profile", h) ;
      memset (prefix, 0, sizeof (memset)) ;
      
      /* scan the various runs separatelly */
      for (jr = 1 ; jr < EP->cols ; jr++)
	{
	  ks1 = arrayReCreate (ks1, 200, double) ;
	  ks2 = arrayReCreate (ks2, 200, double) ;
	  
	  /* find the prefix and the length of the reads */
	  for (xmax = 0, ir = 0 ; ir < LP->rows ; ir++)
	    {
	      ccp = ac_table_printable (LP, ir, 0, "xx") ;
	      if (jr == 1 && ccp[1] != 0 && ccp[1] != '1')
		continue ;
	      if (jr > 1 && ccp[1] != '0' + jr)
		continue ;
	      strcpy (prefix, ccp) ;

	      x = ac_table_int (LP, ir, 1, 0) ;
	      if (x > xmax) xmax = x ;
	    }
	  /* Count the reads aligned at each position */
	  /* start by the integral of the FB, do not assume that the order is correct */
	  for (ir = 0 ; ir < FB->rows ; ir++)
	    {
	      x = ac_table_int (FB, ir, 0, 0) ;
	      a = ac_table_float (FB, ir, jr, 0) ;
	      array (ks1, x, double) = a ;
	    }
	  /* integrate */
	  for (x = a = 0 ; x <= xmax ; x++)
	    {
	      a += array (ks1, x, double) ;
	      array (ks1, x, double) = a ;
	    }
	  /* analyse the LB */
	  for (ir = 0 ; ir < LB->rows ; ir++)
	    {
	      x = ac_table_int (LB, ir, 0, 0) ;
	      a = ac_table_float (LB, ir, jr, 0) ;
	      array (ks2, x, double) = a ;
	    }
	  /* integrate */
	  for (x = a = 0 ; x <= xmax ; x++)
	    {
	      a += array (ks2, x, double) ;
	      array (ks2, x, double) = a ;
	    }
	  
	  /* substract the integral of the LB */
	  for (x = 1 ; x <= xmax ; x++)
	    {
	      a = array (ks1, x, double) ;
	      array (ks1, x, double) = a - array (ks2, x - 1, double) ;
	    }
	  
	  /* now count the errors at each position
	   * the quality per position is
	   -10 * log10 (nb_error/bases_aligned) 
	  */
	  for (ir = 0 ; ir < EP->rows ; ir++)
	    {
	      x = ac_table_int (EP, ir, 0, 0) ;
	      z = ac_table_float (EP, ir, jr, 0) ;
	      a = array (ks1, x, double) ;
	      q = 0 ;
	      if (a < 1000) ;
	      else if (z < 1000 && x < 10) ;
	      else if (z < 1000 && x > xmax - 10) ;
	      else if (z < a/1.0E10) q = 100 ;
	      else q = -10 * log (z/a)/L10 ;
	      
	      vtxtPrintf (txt, "\nQuality_profile %s %d %.0f  %.0f %.0f ", prefix, x, q, z, a) ;
	    }
	  /* complete to xmax, otherwise the graphic displays have a problem */
	  for (++x; x <= xmax ; x++)
	    vtxtPrintf (txt, "\nQuality_profile %s %d 0 0 0 ", prefix, x) ;
	}
    }
  
  if (nOk)
    {
      vtxtPrintf (txt, "\n") ;
      if ( ! ac_parse (db, vtxtPtr (txt), &errors, 0, h))
	fprintf (stderr, "baAddQualityFactors error while writing to the database: %s\n", errors ? errors : "unknown error, sorry") ;
    }
  
  ac_free (h) ;
  ac_free (Ali) ;
  ac_free (h1) ;
  fprintf (stderr, "Added quality_profile to %d runs\n", nOk) ;
  return nOk  ;
} /* baAddQualityFactors */

/*************************************************************************************/
/*************************************************************************************/

static int baGeneIndex2TableDo (BA *ba, BOOL hasGeneId, int iERCC)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ; 
  ACEOUT ao ;
  AC_DB db = ba->db ;
  AC_ITER iter = 0 ;
  AC_KEYSET ks = 0 ;
  AC_OBJ g_any = 0, gene = 0, Run = 0 ;
  AC_TABLE groups = 0, runs = 0, rgs = 0, dr = 0, imap, tbl = 0, erccs = 0 ;
  const char *error = 0 ;
  KEY run ;
  int ir, jr, typ, uu, pass, ok = 0, nnnn = 0 ;
  const char *deepNam[] = { "Group_U", "Group_nU", "Run_U", "Run_nU", 0 } ;
  const char *typNam[] = { "Index", "Seq", "Reads", "kb", "Variance" , 0 } ;
  const char *titleNam[] = { "Index", "Sequences aligned", "Reads aligned", "kb aligned", "Variance" , 0 } ;
  vTXT txt[5] ;

  switch (iERCC)
    {
    case 0:
      ks = ac_dbquery_keyset (db, hprintf (h, "Find Run Is_group && project == \"%s\"", ba->project ? ba->project : "*"), h) ;
      groups = ac_keyset_table (ks, 0, 0, TRUE, h) ;

      ks = ac_dbquery_keyset (db, hprintf (h, "Find Run Is_run && project == \"%s\"", ba->project ? ba->project : "*"), h) ;
      runs = ac_keyset_table (ks, 0, 0, TRUE, h) ;

      iter = ac_query_iter (db, TRUE, hprintf (h, "find Gene %s AND %s"
					       , dictName (ba->target_classDict, ba->target_class)
					   , hasGeneId ? "GeneId" : "NOT cloud AND NOT GeneId"
					       )
			    , 0, h) ;
      fprintf(stderr, "find Gene %s AND %s\n"
	      , dictName (ba->target_classDict, ba->target_class)
	      , hasGeneId ? "GeneId" : "NOT cloud AND NOT GeneId"
	      ) ;
      break ;
    case 1: 
      baErccOptimize (ba) ;
      ks = ac_dbquery_keyset (db, "Find Run Is_group && project == \"ERCC1\"", h) ;
      groups = ac_keyset_table (ks, 0, 0, TRUE, h) ;

      ks = ac_dbquery_keyset (db, "Find Run Is_run && project == \"ERCC1\"", h) ;
      runs = ac_keyset_table (ks, 0, 0, TRUE, h) ;
      erccs = ac_bql_table (db, "select g,m,m2,mm2 from g in class \"gene\", m in g->Mix1 where m>0, m2 in m[1], mm2 = 1000 * m2 order by -4 ",   0, 0, &error, h) ;
      break ;
    case 2:
      ks = ac_dbquery_keyset (db, "Find Run Is_group && project == \"ERCC2\"", h) ;
      groups = ac_keyset_table (ks, 0, 0, TRUE, h) ;

      ks = ac_dbquery_keyset (db, "Find Run Is_run && project == \"ERCC2\"", h) ;
      runs = ac_keyset_table (ks, 0, 0, TRUE, h) ;

      erccs = ac_bql_table (db, "select g,m,m2,mm2 from g in class \"gene\", m in g->Mix1 where m>0, m2 in m[1], mm2 = 1000 * m2 order by -4 ",   0, 0, &error, h) ;
      break ;
    }
  g_any = ac_get_obj (db, "Gene", hprintf(h,"G_Any_*_%s_gene", dictName (ba->target_classDict, ba->target_class)), h) ;
  memset (txt, 0, sizeof(txt)) ;

  for (typ = 0 ; typNam[typ] ; typ++) /* index, seqs, tags, kb, var */
    {
      txt[typ] = vtxtHandleCreate (h) ;

      if (ba->target_class  && ! strcasecmp (dictName (ba->target_classDict, ba->target_class) , "ERCC"))
	vtxtPrintf (txt[typ], "Gene\tLength\t%%GC\t%%Non unique mapping\tChromosome\tStart\tEnd\tTitle\tERCC Mix 1:log2(concentration in attomole/ul)\tERCC Mix2:log2(concentration in attomole/ul)") ;
      else
	vtxtPrintf (txt[typ], "Gene\tGeneId\tNM\tLength\t%%GC\t%%Unicity\tChromosome\tStart\tEnd\tTitle") ;

      for (pass = 0 ; pass < 4 ; pass+=2)
	{
	  rgs = pass < 2 ? groups : runs ;
	  if (typ == 4 && rgs == runs) continue ;
	  /* name the run/group */
	  for (ir = 0 ; rgs && ir < rgs->rows ; ir++)
	    {
	      run = ac_table_key (rgs, ir, 0, 0) ;
	      for (uu = ba->unique ? 0 : 1 ; uu < 2 ; uu+= 2)
		{
		  vtxtPrintf (txt[typ], "\t%s %s %s"
			      , ac_key_name (run)
			      , typNam[typ]
			      , uu == 0 ? "Unique" : "Non unique"
			      ) ;
		}
	    }
	}
      vtxtPrintf (txt[typ], "\nMachine\t\t\t\t\t\t\t\t\t") ;
      for (pass = 0 ; pass < 4 ; pass+=2)
	{
	  rgs = pass < 2 ? groups : runs ;
	  if (typ == 4 && rgs == runs) continue ;
	  
	  /* name the run/group */
	  for (ir = 0 ; rgs && ir < rgs->rows ; ir++)
	    for (uu = ba->unique ? 0 : 1 ; uu < 2 ; uu+= 2)
	      {
		Run = ac_table_obj (rgs, ir, 0, h) ;
		vtxtPrintf (txt[typ], "\t%s %s", Run ? ac_tag_printable (Run, "Machine", "-") : "-", Run ? ac_tag_printable (Run, "Machine_type", "-") : "-") ;
		ac_free (Run) ;
	      }
	}
      vtxtPrintf (txt[typ], "\nSystem\t\t\t\t\t\t\t\t\t") ;
      for (pass = 0 ; pass < 4 ; pass+=2)
	{
	  AC_TABLE sst  = 0 ;
	  AC_OBJ Sample = 0 ;
	  int iSys, iSys1 ;

	  rgs = pass < 2 ? groups : runs ;
	  if (typ == 4 && rgs == runs) continue ;
	  
	  /* name the run/group */
	  for (ir = 0 ; rgs && ir < rgs->rows ; ir++)
	    for (uu = ba->unique ? 0 : 1 ; uu < 2 ; uu+= 2)
	      {
		Run = ac_table_obj (rgs, ir, 0, h) ;
		Sample = Run ? ac_tag_obj (Run, "Sample", h) : 0 ;
		sst = Sample ? ac_tag_table (Sample, "Systm", h) : 0 ;
		for (iSys = iSys1 = 0 ; sst && iSys < sst->rows ; iSys++)
		  {
		    vtxtPrintf (txt[typ], "%s%s%s%s"
				, iSys1++ ? ", " : "\t"
				, ac_table_printable (sst, iSys, 0, "")
				, ac_table_printable (sst, iSys, 1, 0) ? ", " : ""
				, ac_table_printable (sst, iSys, 1, "") 
				) ;
		  } 
		if (!iSys1)  vtxtPrintf (txt[typ], "\t") ;
		ac_free (Run) ;
		ac_free (Sample) ;
		ac_free (sst) ;
	      }
	}
      vtxtPrintf (txt[typ], "\nTissue\t\t\t\t\t\t\t\t\t") ;
      for (pass = 0 ; pass < 4 ; pass+=2)
	{
	  AC_TABLE sst  = 0 ;
	  AC_OBJ Sample = 0 ;
	  int iSys, iSys1 ;

	  rgs = pass < 2 ? groups : runs ;
	  if (typ == 4 && rgs == runs) continue ;
	  
	  /* name the run/group */
	  for (ir = 0 ; rgs && ir < rgs->rows ; ir++)
	    for (uu = ba->unique ? 0 : 1 ; uu < 2 ; uu+= 2)
	      {
		Run = ac_table_obj (rgs, ir, 0, h) ;
		Sample = Run ? ac_tag_obj (Run, "Sample", h) : 0 ;
		sst = Sample ? ac_tag_table (Sample, "Tissue", h) : 0 ;
		for (iSys = iSys1 = 0 ; sst && iSys < sst->rows ; iSys++)
		  {
		    vtxtPrintf (txt[typ], "%s%s%s%s"
				, iSys1++ ? ", " : "\t"
				, ac_table_printable (sst, iSys, 0, "")
				, ac_table_printable (sst, iSys, 1, 0) ? ", " : ""
				, ac_table_printable (sst, iSys, 1, "") 
				) ;
		  } 
		if (!iSys1)  vtxtPrintf (txt[typ], "\t") ;
		ac_free (Run) ;
		ac_free (Sample) ;
		ac_free (sst) ;
	      }
	}
      vtxtPrintf (txt[typ], "\nSample\t\t\t\t\t\t\t\t\t") ;
      for (pass = 0 ; pass < 4 ; pass+=2)
	{
	  rgs = pass < 2 ? groups : runs ;
	  if (typ == 4 && rgs == runs) continue ;
	  
	  /* name the run/group */
	  for (ir = 0 ; rgs && ir < rgs->rows ; ir++)
	    for (uu = ba->unique ? 0 : 1 ; uu < 2 ; uu+= 2)
	      {
		Run = ac_table_obj (rgs, ir, 0, h) ;
		vtxtPrintf (txt[typ], "\t%s", Run ? ac_tag_printable (Run, "Sample", "-") : "-") ;
		ac_free (Run) ;
	      }
	}
      vtxtPrintf (txt[typ], "\n%%GC in run\t\t\t\t\t\t\t\t\t") ;
      for (pass = 0 ; pass < 4 ; pass+=2)
	{
	  AC_OBJ Ali = 0 ;
	  AC_TABLE atgc = 0 ;
	  float gc ;
	  rgs = pass < 2 ? groups : runs ;
	  if (typ == 4 && rgs == runs) continue ;
	  
	  /* name the run/group */
	  for (ir = 0 ; rgs && ir < rgs->rows ; ir++)
	    for (uu = ba->unique ? 0 : 1 ; uu < 2 ; uu+= 2)
	      {	
		Run = ac_table_obj (rgs, ir, 0, h) ;
		Ali = ac_tag_obj (Run, "Ali", h) ;
		atgc = Ali ? ac_tag_table (Ali, "ATGC_kb", h) : 0 ;
		gc = atgc ? ac_table_float (atgc, 0, 2, 0) + ac_table_float (atgc, 0, 3, 0) : 0 ;
		vtxtPrintf (txt[typ], "\t%.1f", gc/10) ;
		ac_free (atgc) ;
		ac_free (Ali) ;
		ac_free (Run) ;
	      }
	}
      vtxtPrintf (txt[typ], "\nRun description\t\t\t\t\t\t\t\t\t") ;
      for (pass = 0 ; pass < 4 ; pass+=2)
	{
	  rgs = pass < 2 ? groups : runs ;
	  if (typ == 4 && rgs == runs) continue ;
	  
	  /* name the run/group */
	  for (ir = 0 ; rgs && ir < rgs->rows ; ir++)
	    for (uu = ba->unique ? 0 : 1 ; uu < 2 ; uu+= 2)
	      {
		Run = ac_table_obj (rgs, ir, 0, h) ;
		vtxtPrintf (txt[typ], "\t%s", Run ? ac_tag_printable (Run, "Title", ac_name(Run)) : "-") ;
		ac_free (Run) ;
	      }
	}
      
      vtxtPrintf (txt[typ], "\nLower measurable index\t\t\t\t\t\t\t\t\t") ;
      for (pass = 0 ; pass < 4 ; pass+=2)
	{
	  AC_OBJ Ali = 0 ;
	  rgs = pass < 2 ? groups : runs ;
	  if (typ == 4 && rgs == runs) continue ;
	  
	  /* name the run/group */
	  rgs = pass < 2 ? groups : runs ;
	  for (ir = 0 ; rgs && ir < rgs->rows ; ir++)
	    {
	      Run = ac_table_obj (rgs, ir, 0, 0) ;
	      Ali = Run ? ac_tag_obj (Run, "Ali", h1) : 0 ;

	      vtxtPrintf (txt[typ], "\t%d", Ali ? (int)(ac_tag_float (Ali, "Low_index", 0)) : 0) ;
	      ac_free (Run) ;
	      ac_free (Ali) ;
	    }
	}
      vtxtPrint (txt[typ], "\n") ;
      
      if (g_any)
	{
	  switch (typ)
	    {
	    case 0: vtxtPrint (txt[typ], "kb aligned in ") ; break ;
	    case 1: vtxtPrint (txt[typ], "Seq aligned in ") ; break ;
	    case 2: vtxtPrint (txt[typ], "Reads aligned in ") ; break ;
	    case 3: vtxtPrint (txt[typ], "kb aligned in ") ; break ;
	    case 4: vtxtPrint (txt[typ], "variance in ") ; break ;
	    }
	  {
	    vtxtPrintf (txt[typ], "%s, excluding genes capturing over 1/1000 of the total\t\t\t\t\t\t\t\t\t", dictName (ba->target_classDict, ba->target_class)) ;
	  }
	  h1 = ac_new_handle () ;
	  for (pass = 0 ; pass < 4 ; pass+=2)
	    {
	      rgs = pass < 2 ? groups : runs ;
	      if (typ == 4 && rgs == runs) continue ;
	      
	      /* name the run/group */
	      for (ir = 0 ; rgs && ir < rgs->rows ; ir++)
		for (uu = ba->unique ? 0 : 1 ; uu < 2 ; uu+= 2)
		  {
		    dr = g_any ? ac_tag_table (g_any, deepNam[pass + uu], h1) : 0 ;
		    ok = FALSE ;
		    run = ac_table_key (rgs, ir, 0, 0) ;
		    for (jr = 0 ; dr && run && jr < dr->rows ; jr++)
		      if (ac_table_key (dr, jr, 0, 0) == run)
			{
			  ok = TRUE ;
			  switch (typ)
			    {
			    case 0: vtxtPrintf (txt[typ], "\t%.2f", ac_table_float (dr, jr, 6, 0)) ; break ;
			    case 1: vtxtPrintf (txt[typ], "\t%.2f", ac_table_float (dr, jr, 2, 0)) ; break ;
			    case 2: vtxtPrintf (txt[typ], "\t%.2f", ac_table_float (dr, jr, 4, 0)) ; break ;
			    case 3: vtxtPrintf (txt[typ], "\t%.2f", ac_table_float (dr, jr, 6, 0)) ; break ;
			    case 4: vtxtPrintf (txt[typ], "\t%.2f", ac_table_float (dr, jr, 8, 0)) ; break ;
			    }
			  break ;
			}
		    if (! ok)
		      {
			vtxtPrintf (txt[typ], "\t") ;
		      }
		  }
	    }
	  ac_free (h1) ;
	  vtxtPrint (txt[typ], "\n") ;
	}
      switch (typ)
	{
	case 0: vtxtPrint (txt[typ], "kb aligned in ") ; break ;
	case 1: vtxtPrint (txt[typ], "Seq aligned in ") ; break ;
	case 2: vtxtPrint (txt[typ], "Reads aligned in ") ; break ;
	case 3: vtxtPrint (txt[typ], "kb aligned in ") ; break ;
	case 4: vtxtPrint (txt[typ], "variance in ") ; break ;
	}
      {
	vtxtPrintf (txt[typ], "any target\t\t\t\t\t\t\t\t\t") ;
      }
      h1 = ac_new_handle () ;
      for (pass = 0 ; pass < 4 ; pass+=2)
	{
	  AC_TABLE tbl1 ;
	  AC_OBJ Ali = 0 ;
	  rgs = pass < 2 ? groups : runs ;
	  if (typ == 4 && rgs == runs) continue ;
	  
	  /* name the run/group */
	  for (ir = 0 ; rgs && ir < rgs->rows ; ir++)
	    {
	      Run = ac_table_obj (rgs, ir, 0, 0) ;
	      Ali = Run ? ac_tag_obj (Run, "Ali", h1) : 0 ;
	      tbl1 = Ali ? ac_tag_table (Ali, "nh_Ali", h1) : 0 ;
	      for (jr = 0 ; tbl1 && jr < tbl1->rows ; jr++)
		if (! strcasecmp ("any", ac_table_printable (tbl1, jr, 0, "toto")))
		  {
		    ok = TRUE ;
		    switch (typ)
		      {
		      case 0: vtxtPrintf (txt[typ], "\t%.2f", ac_table_float (tbl1, jr, 5, 0)) ; break ;
		      case 1: vtxtPrintf (txt[typ], "\t%.2f", ac_table_float (tbl1, jr, 1, 0)) ; break ;
		      case 2: vtxtPrintf (txt[typ], "\t%.2f", ac_table_float (tbl1, jr, 3, 0)) ; break ;
		      case 3: vtxtPrintf (txt[typ], "\t%.2f", ac_table_float (tbl1, jr, 5, 0)) ; break ;
		      case 4: vtxtPrintf (txt[typ], "\t%.2f", ac_table_float (tbl1, jr, 7, 0)) ; break ;
		      }
		    break ;
		  }
	      if (! ok)
		{
		  vtxtPrintf (txt[typ], "\t0") ;
		}
	      ac_free (Run) ;
	      ac_free (Ali) ;
	      ac_free (tbl1) ;
	    }
	}
      ac_free (h1) ;
      vtxtPrintf (txt[typ], "\n") ;

      if (ba->target_class  && ! strcasecmp (dictName (ba->target_classDict, ba->target_class) , "ERCC"))
	vtxtPrintf (txt[typ], "Gene\tLength\t%%GC\t%%Non unique mapping\tChromosome\tStart\tEnd\tTitle\tERCC Mix 1:log2(concentration in attomole/ul)\tERCC Mix2:log2(concentration in attomole/ul)") ;
      else
	vtxtPrintf (txt[typ], "Gene\tGeneId\tNM\tLength\t%%GC\t%%Unicity\tChromosome\tStart\tEnd\tTitle") ;

      for (pass = 0 ; pass < 4 ; pass+=2)
	{
	  rgs = pass < 2 ? groups : runs ;
	  if (typ == 4 && rgs == runs) continue ;
	  
	  /* name the run/group */
	  for (ir = 0 ; rgs && ir < rgs->rows ; ir++)
	    {
	      run = ac_table_key (rgs, ir, 0, 0) ;
	      for (uu = ba->unique ? 0 : 1 ; uu < 2 ; uu+= 2)
		{
		  vtxtPrintf (txt[typ], "\t%s %s %s"
			      , ac_key_name (run)
			      , typNam[typ]
			      , uu == 0 ? "Unique" : "Non unique"
			      ) ;
		}
	    }
	}

      vtxtPrintf (txt[typ], "\nRun description\t\t\t\t\t\t\t\t\t") ;
      for (pass = 0 ; pass < 4 ; pass+=2)
	{
	  rgs = pass < 2 ? groups : runs ;
	  if (typ == 4 && rgs == runs) continue ;
	  
	  /* name the run/group */
	  for (ir = 0 ; rgs && ir < rgs->rows ; ir++)
	    for (uu = ba->unique ? 0 : 1 ; uu < 2 ; uu+= 2)
	      {
		Run = ac_table_obj (rgs, ir, 0, h) ;
		vtxtPrintf (txt[typ], "\t%s", Run ? ac_tag_printable (Run, "Title", ac_name(Run)) : "-") ;
		ac_free (Run) ;
	      }
	}
    }	  
  nnnn = 0 ;
  while (TRUE)
    {
      ac_free (gene), ac_free (h1) ;
      tbl = 0 ;
      h1 = ac_new_handle () ;
      switch (iERCC)
	{
	case 0:
	  gene = ac_iter_obj (iter) ;
	  break ;
	default:
	  gene = erccs && nnnn < erccs->rows ? ac_table_obj (erccs, nnnn, 0, h1) : 0 ; 
	  break ;
	}
      if (! gene)
	break ;
      nnnn++ ;

      imap = ac_tag_table (gene, "IntMap", h1) ;   
      if (! imap)
	{
	  AC_OBJ pg2 = ac_tag_obj (gene, "Gene_model", h1) ;
	  if (pg2)
	    imap = ac_tag_table (pg2, "IntMap", h1) ;   
	}
      for (typ = 0 ; typNam[typ] ; typ++) /* index, seqs, tags, kb, var */
	{
	  const char *ccp = ac_name(gene) ;
	  if (!strncmp (ccp, "X__", 3)) ccp += 3 ;
	  vtxtPrintf (txt[typ], "\n%s", ccp) ;

	  if (ba->target_class  && ! strcasecmp (dictName (ba->target_classDict, ba->target_class) , "ERCC")) ;
	  else
	    { 
	      vtxtPrintf (txt[typ], "\t") ;
	      tbl = ac_tag_table (gene, "GeneId", h1) ;
	      if (tbl && tbl->rows)
		{
		  for (ir = 0 ; ir < tbl->rows ; ir++)
		    vtxtPrintf (txt[typ],"%s%s", ir ? ";" : "", ac_table_printable (tbl,ir,0,"-")) ;
		}
	      else
		vtxtPrint (txt[typ],"-") ;
	      vtxtPrint (txt[typ],"\t") ;
	      ac_free (tbl) ;
	      tbl = ac_tag_table (gene, "NM_Id", h1) ;
	      if (tbl && tbl->rows)
		{
		  for (ir = 0 ; ir < tbl->rows ; ir++)
		    vtxtPrintf (txt[typ],"%s%s", ir ? ";" : "", ac_table_printable (tbl,ir,0,"-")) ;
		}
	      else
		vtxtPrint (txt[typ],"-") ;	
	    }

	  vtxtPrintf (txt[typ],"\t%d", ac_tag_int (gene, "Length", 0)) ;
	  vtxtPrintf (txt[typ],"\t%d", ac_tag_int (gene, "GC_percent", 0)) ;
	  vtxtPrintf (txt[typ],"\t%.1f", ac_tag_float (gene, " Unicity_index", 0)) ;

	  if (imap)
	    vtxtPrintf (txt[typ], "\t%s\t%d\t%d"
			, ac_table_printable (imap, 0, 0, "-")
			, ac_table_int (imap, 0, 1, 0)
			, ac_table_int (imap, 0, 2, 0)
			) ;
	  else
	    vtxtPrintf (txt[typ], "\t-\t0\t0") ;

	  vtxtPrintf (txt[typ], "\t%s"
		      , ac_tag_printable (gene, "Title", "-")
		      ) ;

	  if (ba->target_class  && ! strcasecmp (dictName (ba->target_classDict, ba->target_class) , "ERCC"))
	    {
	      if (ac_has_tag (gene, "ERCC"))
		{
		  tbl = ac_tag_table (gene, "Mix1", h1) ;
		  if (tbl && tbl->rows && tbl->cols >= 2)
		    vtxtPrintf (txt[typ],"\t%s", ac_table_printable (tbl,0,1,"-")) ;
		  else
		    vtxtPrint (txt[typ],"\t-") ;
		  ac_free (tbl) ;
		  
		  tbl = ac_tag_table (gene, "Mix2", h1) ;
		  if (tbl && tbl->rows && tbl->cols >= 2)
		    vtxtPrintf (txt[typ],"\t%s", ac_table_printable (tbl,0,1,"-")) ;
		  else
		    vtxtPrint (txt[typ],"\t-") ;
		  ac_free (tbl) ;
		}
	      else
		vtxtPrint (txt[typ],"\t-\t-") ;
	    }	      
	}

      for (pass = 0 ; pass < 4 ; pass+=2)
	{
	  rgs = pass < 2 ? groups : runs ;
	  if (typ == 4 && rgs == runs) continue ;

	  /* scan the gene to locate the run/group */
	  for (ir = 0 ; rgs && ir < rgs->rows ; ir++)
	    {
	      /* for each run, give the values unique/non-unique together */
	      for (uu = ba->unique ? 0 : 1 ; uu < 2 ; uu+= 2)
		{
		  dr = ac_tag_table (gene, deepNam[pass + uu], h1) ;
		  ok = FALSE ;
		  run = ac_table_key (rgs, ir, 0, 0) ;
		  for (jr = 0 ; dr && run && jr < dr->rows ; jr++)
		    if (ac_table_key (dr, jr, 0, 0) == run)
		      {
			ok = TRUE ;
			if (txt[0]) 
			  {
			    float z = ac_table_float (dr, jr, 1, -999) ;
			    vtxtPrintf (txt[0], "\t") ;
			    if (z > 0)
			      vtxtPrintf (txt[0], "%.2f", z) ;
			  }
			if (txt[1]) vtxtPrintf (txt[1], "\t%.2f", ac_table_float (dr, jr, 2, 0)) ;
			if (txt[2]) vtxtPrintf (txt[2], "\t%.2f", ac_table_float (dr, jr, 4, 0)) ;
			if (txt[3]) vtxtPrintf (txt[3], "\t%.2f", ac_table_float (dr, jr, 6, 0)) ;
			if (txt[4]) 
			  {
			    float z = ac_table_float (dr, jr, 8, -999) ;
			    vtxtPrint (txt[4], "\t") ;
			    if (z > 0)
			      vtxtPrintf (txt[4], "%.2f", z) ;
			  }
			break ;
		      }
		  if (! ok)
		    {
		      for (typ = 0 ; typNam[typ] ; typ++) 
			if (txt[typ]) 
			  vtxtPrint (txt[typ], typ == 0 || typ == 4 ? "\t" : "\t0") ;
		    }
		}
	    }
	}
    }
  for (typ = 0 ; typNam[typ]  ; typ++) /* index, seqs, tags, kb, variance */
    {
      ao = aceOutCreate (ba->outFileName, messprintf (".%s.%s.txt"
						      , iERCC ? (iERCC == 1 ? "ERCC1" : "ERCC2") : (hasGeneId ? "geneId" : "noGeneId")
						      , typNam[typ]
						      )
			 , FALSE, h) ;

      aceOutf (ao, "%s %s%s\n" 
		  , titleNam[typ]
		  , ba->title ? ba->title : ""
		  , hasGeneId ? "with geneId" : ""	     
	       ) ;
      if (typ == 0)
	aceOutf (ao, "Below a minimal number of reads, the index is skipped since it cannot be evaluated accurately,\ni.e. the low indexes can only be measured in larger runs or in groups of runs\n")  ;
      if (typ == 4)
	aceOutf (ao, "Below a minimal number of runs with the minimal number of reads, the variance in the group cannot be evaluated accurately, hence the value is skipped\n") ;

      aceOutf (ao, "\n%s Genes %s %s\n",   dictName (ba->target_classDict, ba->target_class), typNam[typ], timeShowNow()) ;
      aceOut (ao, vtxtPtr(txt[typ])) ;
      aceOut (ao, "\n") ;
    }
  ac_free (h1) ;
  ac_free (h) ;

  return 0 ;
} /* baGeneIndex2TablesDo */

static int baGeneIndex2Table (BA *ba, BOOL hasGeneId)
{
  if (! strcasecmp (dictName (ba->target_classDict, ba->target_class), "ERCC"))
    {
      if (! hasGeneId)
	{
	  baGeneIndex2TableDo (ba, FALSE, 0) ;
	  baGeneIndex2TableDo (ba, FALSE, 1) ;
	  baGeneIndex2TableDo (ba, FALSE, 2) ;
	}
    }
  else
    {
      baGeneIndex2TableDo (ba, hasGeneId, 0) ;
    }
  return 0 ;
} /* baGeneIndex2Tables */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/* export only the best hits */
static int baExportReadSuffix (BA *ba, int pass, int *p5p15, int *p3p15, int *p5p25, int *p3p25) 
{  
  AC_HANDLE h = ac_new_handle () ;
  DICT *dnaDict = ba->dnaDict ;
  HIT *up ;  HIT2 *up2 ;  HIT3 *up3 ;
  int i, ii, score = 0, tag = 0, nn = 0 ;
  ACEOUT aos, aop = aceOutCreate (ba->outFileName, hprintf (h, ".prefix.%d", pass), ba->gzo, h) ;
  ACEOUT aosnew = aceOutCreate (ba->outFileName,   hprintf (h, ".suffix.%d", pass), ba->gzo, h) ;
  ACEOUT aosknown = aceOutCreate (ba->outFileName, hprintf (h, ".clipped_adaptor.%d", pass), ba->gzo, h) ;
  ACEOUT aotp = aceOutCreate (ba->outFileName,     hprintf (h, ".target_prefix.%d", pass), ba->gzo, h) ;
  ACEOUT aots = aceOutCreate (ba->outFileName,     hprintf (h, ".target_suffix.%d", pass), ba->gzo, h) ;
  
  for (ii = score = 0, up = arrp (ba->hits, 0, HIT) ; ii < arrayMax (ba->hits) ; up++, ii++)
    { 
      const char *ccp = dictName (ba->tagDict, up->tag) ;
      const char *ccq = ccp + strlen (ccp) - 1  ;

      if ((pass == 1 && *ccq == '<') || (pass == 2 && *ccq == '>'))
	continue ;
      up2 = arrp (ba->hits2, up->nn, HIT2) ;
      up3 = arrp (ba->hits3, up->nn, HIT3) ;
      if (up->score && tag != up->tag)
	{
	  tag = up->tag ; /* so we export only once per tag */
	  nn ++ ;
	  ccp = up3->suffix ? dictName (dnaDict, up3->suffix) : 0 ;
	  if (ccp && (!ccp[1] || !ccp[2] || !ccp[3] || !strncmp (ccp, "AAAA", 4)))ccp = 0 ;
	  aos = (ccp && (*ccp == ace_upper(*ccp))) ? aosknown : aosnew ;
	  for (i = 0 ; i < up2->mult ; i++)
	    {
	      if (up3->prefix) 
		{ 
		  ccq =  dictName (dnaDict, up3->prefix) ; 
		  aceOutf (aop, "%s\n", ccq) ; 
		  if (ccq && strlen (ccq) >= 15) *(p5p15)+=up2->mult; 
		  if (ccq && strlen (ccq) >= 25) *(p5p25)+=up2->mult; 
		}
	      if (ccp) 
		{ 
		  aceOutf (aos, "%s\n", ccp) ; 
		  if (strlen (ccp) >= 15) (*p3p15)+=up2->mult; 
		  if (strlen (ccp) >= 25) (*p3p25)+=up2->mult; 
		}
	      if (up3->targetPrefix) aceOutf (aotp, "%s\n", dictName (dnaDict, up3->targetPrefix)) ;
	      if (up3->targetSuffix) aceOutf (aots, "%s\n", dictName (dnaDict, up3->targetSuffix)) ;
	    }
	}  
    }
  ac_free (h) ;

  return nn ;
} /* baExportReadSuffix */

/*************************************************************************************/
/* export only the suffix, separating the 2 reads in case of paired end sequencing */
static int baExportSuffix (BA *ba)
{
  int pass, nn = 0 ;
  int p5p15 = 0 ;
  int p5p25 = 0 ;
  int p3p15 = 0 ;
  int p3p25 = 0 ;
  AC_HANDLE h = ac_new_handle () ; 
  ACEOUT ao = aceOutCreate (ba->outFileName,  ".partial.p3p5", FALSE, h) ;

  for (pass = 1 ; pass < (ba->pair ? 3 : 2) ; pass++)
    nn += baExportReadSuffix (ba,  ba->pair ? pass : 0, &p5p15, &p3p15, &p5p25, &p3p25) ;
  aceOutf (ao, "Partial_5p %d %d\tPartial_3p %d %d\n"
	   , p5p15, p3p15, p5p25, p3p25
	   ) ;
 
  fprintf (stderr, "Suffix done %s\n", timeShowNow ()) ;
  ac_free (h) ;
  return nn ;
} /* baExportSuffix */

/*************************************************************************************/
/* export only the best hits */
static int baExportMito (BA *ba) 
{  
  AC_HANDLE h = ac_new_handle () ;
  HIT *up ;
  HIT2 *up2 ;
  int mito, chloro, spikeIn, virus, bacteria, ii, score = 0, nn = 0 ;
  ACEOUT ao = aceOutCreate (ba->outFileName, ".mito", ba->gzo, h) ;

  dictAdd (ba->target_classDict, "A_mito", &mito) ;
  dictAdd (ba->target_classDict,  "C_chloro", &chloro) ;
  dictAdd (ba->target_classDict,  "0_SpikeIn", &spikeIn) ;
  dictAdd (ba->target_classDict,  "v_virus", &virus) ;
  dictAdd (ba->target_classDict,  "b_bacteria", &bacteria) ;


  for (ii = score = 0, up = arrp (ba->hits, 0, HIT) ; ii < arrayMax (ba->hits) ; up++, ii++)
    {
      up2 = arrp (ba->hits2, up->nn, HIT2) ;
      if (up2->target_class != mito && up2->target_class != chloro && up2->target_class != spikeIn  && up2->target_class != virus && up2->target_class != bacteria) continue ;
      if (up->score)
	{
	  nn ++ ;
	  baExportOneHit (ao, ba, up) ;
	}  
    }
  
  ac_free (h) ;
  fprintf (stderr, "Mito done %s\n", timeShowNow ()) ;
  return nn ;
} /* baExportMito */

/*************************************************************************************/
/*************************************************************************************/
/* Parse the coordinates of the target genes */
static int baParseGeneRemapFile (BA *ba)
{ 
  AC_HANDLE h = ac_new_handle () ;
  register char *cp ;
  HIT* up ;
  int nn = 0, gene, a1, a2, x1, x2, target_class, chrom ;
  ACEIN ai = aceInCreate (ba->geneRemapFile, FALSE, h) ;

  if (! ba->targetDict)
    ba->targetDict = dictHandleCreate (10000, ba->h) ;
  if (! ba->geneDict)
    ba->geneDict = dictHandleCreate (100000, ba->h) ;
  ba->gene2map = arrayHandleCreate (100000, HIT, ba->h) ;
  
  if (ai) 
    { 
      aceInSpecial (ai, "\n") ;
      while (aceInCard (ai))
	{
	  cp = aceInWord (ai) ;
	  if (! cp || *cp == '#')
	    continue ;
	  dictAdd (ba->target_classDict, cp, &target_class) ;
	  
	  if (! aceInStep (ai, '\t') || ! (cp = aceInWord (ai)))
	    continue ;
	  /* drop the transcript name */
	  
	  if (! aceInStep (ai, '\t') || ! aceInInt (ai, &(x1)) ||
	      ! aceInStep (ai, '\t') || ! aceInInt (ai, &(x2)) 
	      )
	    continue ;
	  
	  if (! aceInStep (ai, '\t') || ! (cp = aceInWord (ai)))
	    continue ;
	  dictAdd (ba->targetDict, cp, &chrom) ;
	  
	  if (! aceInStep (ai, '\t') || ! aceInInt (ai, &(a1)) ||
	      ! aceInStep (ai, '\t') || ! aceInInt (ai, &(a2)) 
	      )
	    continue ;

	  if (! aceInStep (ai, '\t') || ! (cp = aceInWord (ai)))
	    continue ;
	  if (! dictFind (ba->targetDict, cp, &gene))
	    continue ;

	  up = arrayp (ba->gene2map, gene, HIT) ;
	  if (! up->gene)
	    { 
	      nn++ ;
	      up->gene = gene ;
	      up->target = chrom ;
	      up->class = target_class ;
	      up->x1 = a1 ; up->x2 = a2 ;
	    }
	  if (a1 < a2 && a1 < up->x1) up->x1 = a1 ; 
	  if (a1 < a2 && a2 > up->x2) up->x2 = a2 ; 
	  if (a1 > a2 && a1 > up->x1) up->x1 = a1 ; 
	  if (a1 > a2 && a2 < up->x2) up->x2 = a2 ; 
	}
    }

  if (! arrayMax (ba->gene2map))
    arrayDestroy (ba->gene2map) ;
  ac_free (h) ;
  fprintf (stderr, "Parsed the coordinates of %d genes %s\n", nn, timeShowNow ()) ;
  return nn ;
} /* aParseGeneSpongeFile */

/*************************************************************************************/
/*************************************************************************************/

static int baParseIdFile (BA *ba, DICT *dict, KEYSET tags, KEYSET tag2names, KEYSET tag2stacks, Stack ss[64], AC_HANDLE h0)
{ 
  AC_HANDLE h ;
  register const char *ccp, *ccq ;
  register char *cp ;
  char buf[10000], prefix[100] ;
  int i, n1 = 0, n2 = 0, nid, mult, pass, pos, tag, ventilate = ba->ventilate, nStack = 0 ;
  int prefixLn = 0 ;
  ACEIN ai ;
  Stack s = ss[0] ;

  memset (prefix, 0, sizeof(prefix)) ;
  for (pass = 0 ; pass < 2 ; pass++)
    {
      h = ac_new_handle () ;
      ai = aceInCreate (ba->id_filename, ba->gzi, h) ;
      n1 = n2 = nid = 0 ;
      if (pass)
	{
	  s = ss[0] = stackHandleCreate (1000000, h0) ;
	}
      memset (buf, 0, sizeof(buf)) ;

      aceInSpecial (ai, "\n") ;
      while (aceInCard (ai))
	{
	  ccp = aceInWord (ai) ;
	  if (!ccp || ! *ccp || *ccp == '#') continue ;
	  n1++ ;

	  if (pass)
	    {
	      dictAdd (dict, ccp, &tag) ;
	      keySet (tags, tag) = n2 ; 
	    }
	  mult = fastcMultiplicity (ccp, 0, 0) ;
	  for (i = 0 ; i < mult ; i++)
	    if (aceInStep (ai, '\t') && (ccp = aceInWord (ai)))
	      {
		if (1) /* hack, the probe was named /1/1 rather than /1 */
		  {
		    char *cqh ;
		    cqh = strstr (ccp, "/1/1") ;
		    if (cqh && cqh[4]) cqh[2] = 0 ;
		    cqh = strstr (ccp, "/2/2") ;
		    if (cqh && cqh[4]) cqh[2] = 0 ;
		  }
		nid++ ;
		if (pass == 0) /* look for longest constant prefix */
		  {
		    if (prefix[0])
		      {
			for (cp = prefix, ccq = ccp ; *cp && *cp == *ccq ; cp++, ccq++) ;
			*cp = 0 ;
		      }
		    else
		      strncpy (prefix, ccp, 99) ;
		  }
		else   /* store 1+hash as the first char of the name, so we hash only once per name */
		  {
		    if (strlen (ccp) < prefixLn)
		      messcrash ("In %s name %s shorter than prefix %s", dictAdd (dict, ccp, &tag), ccp, prefix) ;
		    ccp += prefixLn ;
		    {
		      register unsigned char x = 0 ;
		      register int rotate = 3 ;
		      register int leftover =  5 ;
		      /* drop the last char which is the /1 /2 of the orientation of the read */
		      ccq = ccp + strlen (ccp) - 1 ; x = 3 ;
		      while (x>0 && ccq > ccp && *ccq != '/') { ccq-- ; x-- ; }
		      for (x = 0 ; ccq >= ccp ; ccq--)
			x = (*(unsigned char *)ccq) ^ ((x << rotate) | ( x >> leftover)) ; 
		      buf[0] = 1 + (x % ventilate) ;
		    }
		    
		    pos = stackMark (s) ;
		    keySet (tag2names, n2) = pos ;
		    keySet (tag2stacks, n2) = nStack ;
		    n2++ ;
		    pushText (s, buf) ;
		    catText (s, ccp) ;
		    if (pos > 1000000 - 100)
		      {
			s = ss[++nStack] = stackHandleCreate (1000000, h0) ;
			if (nStack > (1 << 12) ) messcrash ("nStack = %d > 4096 in baParseIdFile", nStack) ;
		      }
		  }
	      }
	}
      if (pass == 0) 
	{
	  prefixLn = strlen (prefix) ;
	  fprintf (stderr, "Prefix length %d :: %s common to %d identifiers\n", prefixLn, prefix, nid) ;
          ac_free (h) ;
	}
    }
  nStack++ ;
  fprintf (stderr, "%s Parsed file %s, %d stacks found %d tag names corresponding to %d original names\n",  timeShowNow (), ba->id_filename, nStack, n1, n2) ;
  {
    int mx ;
    messAllocMaxStatus (&mx) ;   fprintf (stderr, "// max memory %d Mb\n", mx) ;
  }

  if (prefixLn)
    {
      ACEOUT ao = aceOutCreate (ba->outFileName, ".prefix", ba->gzo, h) ; 
      
      aceOutf (ao, "%s\n", prefix) ;
      ba->prefixId = strnew (prefix, ba->h) ;
    }

  ac_free (h) ;
  return nStack ;
} /* baParseIdFile */

/*************************************************************************************/

static void baVentilateHitFilePerOriginalIdentifier (BA *ba)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, n1 = 0, n2 = 0, n3 = 0, t1, tag, score, mult, hh, ventilate = ba->ventilate, nStack ;
  int t2, t3 ;
  BOOL isFirst ;
  ACEIN ai = aceInCreate (ba->hitFileName, ba->gzi, h) ;
  ACEOUT ao[ventilate] ;
  const char *ccp, *ccr, *prefix ;
  char *cp, *cr ;
  DICT *dict = dictHandleCreate (5000000, h) ;
  Stack ss[5000], s ;
  KEYSET tags = keySetHandleCreate (h) ;
  KEYSET tag2names = keySetHandleCreate (h) ;
  KEYSET tag2stacks = keySetHandleCreate (h) ;

  memset (ss, 0, sizeof(ss)) ;
  nStack = baParseIdFile (ba, dict, tags, tag2names, tag2stacks, ss, h) ;
  prefix = ba->prefixId ;
  for (i = 0 ; i < ba->ventilate ; i++)
    ao[i] = aceOutCreate (messprintf ("%s.%d",ba->outFileName, i+1), ".hits", ba->gzo, h) ; 

  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      cp = aceInWord (ai) ;
      if (!cp || ! *cp || *cp == '#') continue ;
      n1++ ;
      i = strlen (cp) ;
      cr = cp + i - 1 ;
      isFirst = TRUE ; 
      if (*cr == '>') 
	*cr = 0 ;
      if (*cr == '<') 
	{ isFirst = FALSE ;  *cr = 0 ; }

      if (cp && dictFind (dict, cp, &tag) &&
	  aceInStep (ai, '\t') &&
	  aceInInt (ai, &score) &&
	  aceInStep (ai, '\t') &&
      	  aceInInt (ai, &mult) &&
	  aceInStep (ai, '\t') &&
	  (ccr = aceInPos (ai))
	  )
	{
	  n2++ ;
	  t1 = keySet (tags, tag) ;	
	  for (j = 0 ; j < mult ; j++)
	    {
	      t2 = keySet (tag2names, t1 + j) ;
	      t3 = keySet (tag2stacks, t1 + j) ;
	      s = ss [t3] ;
	      if (!stackExists(s)) messcrash ("messup in %s (alias j=%d) ss[%d > %d] in baVentilateHitFilePerOriginalIdentifier", dictName (dict,tag), j, t3, nStack) ;
	      ccp = stackText (s, t2) ;
	      hh = (*ccp) - 1 ;
	      n3++ ;
	      aceOutf (ao[hh], "%s%s%c\t%d\t1\t%s\n", prefix ? prefix : "", ccp + 1, isFirst ? '>' : '<', score, ccr) ;
	    }
	}
    }

  for (i = 0 ; i < ba->ventilate ; i++)
    ac_free (ao[i]) ;

  fprintf (stderr, "Parsed %d lines, recognized %d tags, exported %d ventilated lines in %s.[1-%d]\n", n1, n2, n3, ba->outFileName ?  ba->outFileName : "stdout", ba->ventilate) ;
  if (n2 < n1)
    fprintf (stderr, "ERROR: %d tag names where not recognized in the id_file\n", n2 - n1) ;

  ac_free (h) ;

  return ;
} /* baVentilateHitFilePerOriginalIdentifier */

/*************************************************************************************/
/*************************************************************************************/
/* histogram of number or exact or false base per position and quality factor */
static void baFastqAnalysis (BA *ba)
{
  typedef struct fahStruct { char x1, x2 ; char errors[14] ; } FAH ; 
  AC_HANDLE h = ac_new_handle () ;
  Array aa = arrayHandleCreate (10000000, FAH, h) ;
  Array hh = arrayHandleCreate (10000, int, h) ;
  FAH *up ;
  ACEIN ai = 0 ;
  int probe, n = 0, n1, pos ;
  int q, baseQual = 64, qmax = 0 ;
  int i, good, err, x, x1, x2, xmax = 0, reclip = 8 ;
  const char *ccp, *ccq ;
  char *cp, *cq, *cr ;
  DICT *dict = dictHandleCreate (1000000, h) ;
  ACEOUT ao = 0 ; 
  double l10 = log(10.0) ;
  /* BOOL isFirst ; */

  ai = aceInCreate (ba->hitFileName, ba->gzi, h) ;
  aceInSpecial (ai, "\n\t") ;
  while (aceInCard (ai))
    {
      cp = aceInWord (ai) ;
      if (! cp || ! *cp || *cp == '#') continue ;

      i = strlen (cp) ;
      cr = cp + i - 1 ;
      /* isFirst = TRUE ;  */
      if (*cr == '>') 
	*cr = 0 ;
      if (*cr == '<') 
	{ /*  isFirst = FALSE ; */  *cr = 0 ; }

      if (dictFind (dict, cp, &probe))
	continue ;
      dictAdd (dict, cp, &probe) ;
      aceInWord (ai) ; /* score */
      aceInWord (ai) ; /* mult */
      aceInWord (ai) ; /* toali */
      aceInWord (ai) ; /* ali */
      aceInInt (ai, &x1) ; /* */
      aceInInt (ai, &x2) ; /* */
      aceInWord (ai) ; /* class */
      aceInWord (ai) ; /* gene */
      aceInWord (ai) ; /* target_mult */
      aceInWord (ai) ; /* target */
      aceInWord (ai) ; /* a1 */
      aceInWord (ai) ; /* a2 */
      aceInWord (ai) ; /* nN */
      aceInWord (ai) ; /* nerr */

      up = arrayp (aa, probe, FAH) ;
      if (x2 - reclip > 255)
	x2 = reclip + 255 ;
      if (x2 - reclip > 255)
	messcrash ("sorry to gains space we limit to [1,255] this read %s is aligned up to %d", dictName(dict,probe), x2) ;
      up->x1 = x1 + reclip ; up->x2 = x2 - reclip ;

      cp = aceInWord (ai) ; /* err types and pos */
      i = 0 ;
      while (cp)
	{
	  cq = strstr (cp, ",") ;
	  if (cq) *cq = 0 ;
	  cr = strstr (cp, ":") ;
	  if (cr)
	    {
	      *cr = 0 ;
	      pos = atoi (cp) ;
	      if (i < 13 && pos < 256)
		(up->errors)[i++] = pos ;
	      cp = cr + 1 ;
	      if (*cp && ! strstr (cp, "n") && ! strstr (cp, "o"))
		{
		  /*
		    dictAdd (errDict, cp, &err) ;
		    array (errors, err, long int) += mult ;
		    array (errPos, pos, long int) += mult ;
		  */
		} 
	      else
		{
		  /*
		    array (NNPos, pos, long int) += mult ;
		  */
		} 
	    }
	  cp = (cq ? cq + 1 : 0) ;
	}
    }  
  ac_free (ai) ;

  probe = 0 ;
  ai = aceInCreate (ba->fastqAnalysis, FALSE, h) ;
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (! ccp) continue ;
      n1 = ++n % 4 ;
      switch (n1) 
	{
	case 1: 
	  if (*ccp != '@') 
	    messcrash ("No @ line %d got %s\n", n, ccp) ;
	  probe = 0 ;
	  dictFind (dict, ccp+1, &probe) ;
	  continue ;
	case 2: 
	  continue ;
	case 3:
	  continue ;
	case 0:
	  break ;
	}
      if (probe <= 0)
	continue ;
      /* perform the actual analysis */
      up = arrayp (aa, probe, FAH) ;
      probe = 0 ;
      x1 = up->x1 ; x2 = up->x2 ;
      if (x1 > x2) continue ;
      if (x2 > xmax) xmax = x2 ;
      i = 0 ; err = *up->errors ;
      for (x=1, ccq = ccp ; *ccq ; x++, ccq++) /* bio coordinates */
	{
	  /* x must be inside the re-clipped aligned segment */
	  if (x < x1) continue ;
	  if (x > x2) break ;
	  /* do we have or do we not have an error at x */
	  good = 1 ;
	  if (x == err)
	    {
	      good = 0 ;
	      if (err)
		err = up->errors[++i] ;
	    }
	  q = *ccq - baseQual ;
	  if (q > qmax) qmax = q ;
	  array (hh, 1000 * x + 2*q + good, int)++ ;
	}
    }

  /* export */
  ao = aceOutCreate (ba->outFileName, ".fastqHisto.txt", FALSE, h) ;
  aceOutf (ao, "# %s", timeShowNow () ) ;
  aceOutf (ao, "\n\tCounts at the given position of the exact and mismatching bases of a given quality") ;
  for (q = 1 ; q <= qmax ; q++)
    aceOutf (ao, "\t") ;
  aceOutf (ao, "\tRun position\tPercentage at the given position of the mismatching bases of a given quality") ;
  for (q = 1 ; q <= qmax ; q++)
    aceOutf (ao, "\t\t") ;
  aceOutf (ao, "\tRun position\tObserved quality = -10 log10(F/F+T) among the bases of a given fastq quality") ;
  for (q = 1 ; q <= qmax ; q++)
    aceOutf (ao, "\t\t") ;
  aceOutf (ao, "\nRun position") ;
  for (q = 1 ; q <= qmax ; q++)
    aceOutf (ao, "\tF%d\tT%d", q, q) ;
  aceOutf (ao, "\t\t") ;
  for (q = 1 ; q <= qmax ; q++)
    aceOutf (ao, "\t%%F%d", q) ;
  aceOutf (ao, "\t\t") ;
  for (q = 1 ; q <= qmax ; q++)
    aceOutf (ao, "\tq%d", q) ;

  aceOutf (ao, "\nExpected") ;
  for (q = 1 ; q <= qmax ; q++)
    aceOutf (ao, "\t\t") ;
  aceOutf (ao, "\t\tExpected") ;
  for (q = 1 ; q <= qmax ; q++)
    aceOutf (ao, "\t.3f", exp((double)q * l10/10.0)) ;
  aceOutf (ao, "\t\tExpected") ;
  for (q = 1 ; q <= qmax ; q++)
    aceOutf (ao, "\td",  q) ;
  
  for (x = 1 ; x <= xmax ; x++)
    { 
      aceOutf (ao, "\n%s %d", ba->run, x) ;
      for (q = 1 ; q <= qmax ; q++)
	{
	  x1 =  array (hh, 1000 * x + 2*q + 0, int) ;
	  x2 =  array (hh, 1000 * x + 2*q + 1, int) ;
	  aceOutf (ao, "\t%d\t%d", x1, x2) ;
	}
      aceOutf (ao, "\t\t%s %d", ba->run, x) ;
      for (q = 1 ; q <= qmax ; q++)
	{
	  x1 =  array (hh, 1000 * x + 2*q + 0, int) ;
	  x2 =  array (hh, 1000 * x + 2*q + 1, int) ;
	  if (x1 + x2)
	    aceOutf (ao, "\t%.3f", x1 ? 100.0*x1/(x1+x2) : 0) ; 
	  else
	    aceOutf (ao, "\t") ;
	}
      aceOutf (ao, "\t\t%s %d", ba->run, x) ;
      for (q = 1 ; q <= qmax ; q++)
	{
	  x1 =  array (hh, 1000 * x + 2*q + 0, int) ;
	  x2 =  array (hh, 1000 * x + 2*q + 1, int) ;
	  if (x1 > 0)
	    aceOutf (ao, "\t%.1f", -10 * log((double)x1/(x1+x2))/l10) ;
	  else
	    aceOutf (ao, "\t") ;
	}
    }
  aceOutf (ao, "\n") ;
	   
	
 ac_free (h) ;

  return ;
} /* baFastqAnalysis */


/*************************************************************************************/
/* remap the read to the correct gene, in case splitMrna are defined */

static int baSplitMrnaRemap (BA *ba, HIT *up, HIT2 *up2)
{
  Associator ass = ba->splitMrnaAss ;
  int mrna = up->target ;
  Array bucket ;
  int iBucket ;
  int x1 = up2->a1, x2 = up2->a2 ;
  int da = up2->badPair ? 0 : up2->dPair ;

  const void *vp ;

  if (! ass || ! ba->splitMrnaArray)
    return up->gene ;

  if (da * (x2 - x1) < 0)
    da = 0 ; 
  if (x1 < x2)
    {
      if (x1 + da - 1 > x2) /* use the pair */
	x2 = x1 + da ;
    }
  else
    { 
      int x0 ;
      if (x1 + da < x2) /* use the pair */
	x2 = x1 + da ;
      x0 = x1 ; x1 = x2 ; x2 = x0 ; 
    }

  iBucket = 0 ; bucket = 0 ;
  if (! assFind (ass, assVoid (mrna), 0))
    {
      mrna = up->gene ;
    }
  if (mrna && assFind (ass, assVoid (mrna), 0))
    while (assFindNext (ass, assVoid (mrna), &vp, &bucket, &iBucket))
      {
	int nn = assInt (vp) ;
	SPLITMRNA *wp = arrayp (ba->splitMrnaArray, nn, SPLITMRNA) ;                
	
	if (wp->x1 == 0)
	  return wp->gNewOld ; 
	else if (wp->x1 - 8 <= x1 && wp->x2 + 8 >= x2)
	  return wp->gXX ;
	else if (wp->x2 - 8 > x1 && wp->x2 + 8 < x2)
	  return wp->gNewOld ;
      }
  return up->gene ;
} /* baSplitMrnaRemap */

/*************************************************************************************/

static int baSplitMrnaParse (BA *ba)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
  int gene, mrna, x1, x2, gXX, gNewOld ;
  const char *ccp ;
  DICT *dict = ba->targetDict ;
  Array aa  = 0 ;
  ACEIN ai = aceInCreate (ba->splitMrnaFileName, 0, h) ;
  SPLITMRNA *up ;
  Associator ass ;

  if (! dict)
    dict = ba->targetDict = dictHandleCreate (100000, ba->h) ; 
  if (! ai)
    messcrash ("Sorry, i cannot find the -split_mRNAs file : %s", ba->splitMrnaFileName) ;
  aa = ba->splitMrnaArray = arrayHandleCreate (10000, SPLITMRNA, ba->h) ;
  ass = ba->splitMrnaAss = assHandleCreate (ba->h) ;

  while (aceInCard (ai))
    {
      BOOL anyMrna = FALSE ;

      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#')
	continue ;
      dictAdd (dict, ccp, &gene) ;

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp)
	continue ;
      if (! strcmp (ccp, "*"))
	{
	  anyMrna = TRUE ;
	  mrna = gene ;
	}
      else
	dictAdd (dict, ccp, &mrna) ;
      
      x1 = x2 = -1 ;
      aceInStep (ai, '\t') ;
      aceInInt (ai, &x1) ;
      aceInStep (ai, '\t') ;
      aceInInt (ai, &x2) ;

      if (x1 < 0 || x2 < 0)
	messcrash ("FATAL ERROR: In -split_mRNAs, missing coordinates in mRNA %s, line %d, file %s"
		   , dictName (dict, mrna) 
		   , aceInStreamLine (ai) 
		   , ba->splitMrnaFileName
		   ) ;
      
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp || ! strcmp (ccp, "NULL"))
	continue ;
      gXX = 0 ;
      if (! anyMrna && x1)
	dictAdd (dict, hprintf (h, "%s(%s)", ccp, dictName (dict, gene)), &gXX) ; /* a new name for that fraction of the read */

      aceInStep (ai, '\t') ; /* jump col 5 */
      ccp = aceInWord (ai) ;
      if (! ccp || strstr (ccp, "NULL__") == ccp)
	continue ;
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp)
	continue ;
      dictAdd (dict, ccp, &gNewOld) ; /* name(oldNam) */
      
      assMultipleInsert (ass, assVoid(mrna), assVoid (nn)) ;
      up = arrayp (aa, nn++, SPLITMRNA) ;
      up->gene = gene ;
      up->gNewOld = gNewOld ;
      if (gXX) 
	{
	  up->gXX = gXX ;
	  up->x1 = x1 ;
	  up->x2 = x2 ;
	}
    }

  ac_free (h) ;
  return nn ;
} /* baSplitMrnaParse */

/*************************************************************************************/

static int baParseSelected8kbTargetFile (BA *ba)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
  const char *ccp ;
  ACEIN ai = aceInCreate (ba->selected8kbFileName, 0, h) ;

  if (! ai)
    messcrash ("Sorry, i cannot find the -selected8kbList file : %s", ba->selected8kbFileName) ;

  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (ccp && *ccp != '#')
	dictAdd (ba->selected8kbDict, ccp, &nn) ;
    }

  ac_free (h) ;
  return nn ;
} /* baParseSelected8kbTargetFile */

/*************************************************************************************/

static int baParseSelected5kbTargetFile (BA *ba)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
  const char *ccp ;
  ACEIN ai = aceInCreate (ba->selected5kbFileName, 0, h) ;

  if (! ai)
    messcrash ("Sorry, i cannot find the -selected5kbList file : %s", ba->selected5kbFileName) ;

  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (ccp && *ccp != '#')
	dictAdd (ba->selected5kbDict, ccp, &nn) ;
    }

  ac_free (h) ;
  return nn ;
} /* baParseSelected5kbTargetFile */

/*************************************************************************************/
/*************************************************************************************/
/* Special utility : check correct group hierarchy */
static void baCheckGroupHierearchyError (const char *ccp)
{
  fprintf (stderr, "%s\n%s\n"
	   , "Sorry, the group hierachy of the MetaDB is not acceptable, please modify\n"
	   "The expected structure is as follows\n"
	   , ccp
	   ) ;
  exit (1) ;
} /* baCheckGroupHierearchyError */

/*************************************************************************************/

static void baSetGroupLevel (AC_DB db, int level,  AC_KEYSET ks, AC_HANDLE h)
{
  char *command ;

  if (level <= 0)
    command = "edit -D Group_level" ;
  else
    command = hprintf (h, "edit Group_level %d", level) ;
  ac_command_keyset (db, command, ks, h) ;
  return ;
} /* baSetGroupLevel */

/*************************************************************************************/

static BOOL baCheckGroupHierearchy (BA *ba)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_KEYSET runs = 0, ks1, ks2 ;
  AC_DB db = ba->db ;
  int n, ng, ii ;

  runs = ac_dbquery_keyset (db, hprintf (h, "Find project IS \"%s\" ; follow run ", ba->project), h) ;
	  
  n = ac_keyset_count (ac_ksquery_keyset (runs, "Sublibrary_of && In_group", h)) ;
  if (1)
    {
      if (n > 0)
	baCheckGroupHierearchyError ("The MetaDB should not contain runs \"Sublibrary_of && In_group\"") ;
      n = ac_keyset_count (ac_ksquery_keyset (runs, "COUNT Sublibrary_of > 1", h)) ;
      if (n > 0)
	baCheckGroupHierearchyError ("The MetaDB should not contain runs \"COUNT Sublibrary_of > 1\"") ;
    }
  n = ac_keyset_count (ac_ksquery_keyset (runs, "Add_counts ; COUNT {>Union_of; ! Is_run && ! sublibraries && !add_counts} > 0", h)) ;
  {
    if (n > 0)
      baCheckGroupHierearchyError ("Some add_count groups depend on non-add-count subgroups, try and fix:\nQuery Find Run Add_counts ; COUNT {>Union_of ; ! Is_run && ! sublibraries && !add_counts} > 0)") ;
  }

  ng = ac_keyset_count (ac_ksquery_keyset (runs, "Is_group", h)) ; /* total number of groups */
  
  ks2 = ac_ksquery_keyset (runs, "Sublibrary_of", h) ; /* level -1 */
  baSetGroupLevel (db, -1, ks2, h) ;
  ks2 = ac_ksquery_keyset (runs, "(NOT Sublibrary_of) && (NOT Is_group)", h) ; /* level 0 */
  baSetGroupLevel (db, 0, ks2, h) ;
  for (ii = 0 ; ii < ng + 2 ; ii++)
    {
      ks1 = ks2 ;
      ks2 = ac_ksquery_keyset (ks1, hprintf(h, ">group ; project == %s",ba->project) , h) ;
      if (ii == ng+1)
	{
	  int i ;
	  AC_TABLE tbl = ac_keyset_table (ks2,0, -1, 0, h) ;
	  fprintf (stderr, "ERROR: Looping on group Level %d\n", ii) ;
	  for (i = 0 ; i < tbl->rows ; i++)
	    fprintf (stderr, "... %s\n", ac_table_printable (tbl, i, 0, "")) ;
	  baCheckGroupHierearchyError ("The group hierarchy should be acyclic: >In_group;>In_group... should converge") ;
	}
      baSetGroupLevel (db, ii+1, ks2, h) ;
      if (!ac_keyset_count (ks2))
	break ;
    }
  if (ii >= ng + 1)
    baCheckGroupHierearchyError ("The group hierarchy should b acyclic: >In_group;>In_group... should converge") ;
  return TRUE ;
} /* baCheckGroupHierearchy */

/*************************************************************************************/

static int baGroupLetterProfile (BA *ba)
{
  int groupLevel = 0 ; 

  baCheckGroupHierearchy (ba) ;
  for (groupLevel = 0 ;  ; groupLevel++)
    {        
      BOOL ok = baGroupLetterProfileByLevel (ba, groupLevel) ;
      if (groupLevel > 0 && ! ok)
      break ;
    }
  return 0 ;
} /* baGroupLetterProfile */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  fprintf  (stderr,
	    "// bestali.c:\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, October 2009, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "//   to analyse the output of the aligner:\n"
	    "//   Select the best alignments and collate the data\n"
	    "// Database\n"
	    "//   -db ACEDB : acedb database holding the semantics of the experiments\n"
	    "// Input: the program expects a .hits file exported by clipalign or by the -exportBest option\n"
	    "//   -i file_name [-gzi] : default stdin\n"
	    "//      if the file is called .gz or if the option -gzi is specified, invokes gunzip\n"
	    "//   -o out_file_name [-gzo] : default stdout\n"
	    "//      redirect the output, each option adds its own suffix to the file name\n"
	    "//      if -gzo is specified, invokes gzip and add a .gz suffix to the file name\n"
	    "//   -inFileList file_name\n"
	    "//      The program reads all file listed in the given file_name\n"
	    "// Actions\n"
	    "//   -ventilate n -id_file f.id :\n"
	    "//      ventilate the input (best)ali into a set of n files keyed on the original identifiers\n"
	    "//      The id_file f.id is exported in phase a0 by dna2dna -keepName -O fastc -o f \n"
	    "//   -countBest [-gzo] : \n"
	    "//      count for each sample the aligned bases flat and hierarchic\n"
	    "//      .count suffix is added to out_file_name, gz is never applied\n"
	    "//   -seqc \n"
	    "//      export in addition all required reports in .seqc.tsv format\n"
	    "//   -exportBest : export the best alignments as a .hits file\n"
	    "//      .hits[.gz] suffix is added automatically if missing\n"
	    "//   -exportSuffix : export 4 dna files\n"
	    "//      .prefix[.gz] 5' overhangs, read in the tag, on the antistrand of the tag extending the alignment\n"
	    "//      .suffix[.gz] 3' overhangs, read in the tag, on the strand of the tag extending the alignment\n"
	    "//      .target_prefix[.gz] 30bp on the target upstream of the alignment, in the orientation of the tag\n"
	    "//      .target_suffix[.gz] 30bp on the target downstream of the alignment, in the orientation of the tag\n"
	    "//   -keepTargetPrefix : preserve the target prefix/suffix columns\n"
	    "//       The default is to clean them out (and optionally use -exportSuffix)\n"
	    "//   -exportMito : export .mito file, with all mito position with a mismatch\n"
	    "//   -aliProfile : export .aliProfile file, with ali length and percent and error rates histograms\n"
	    "//   -errorProfile : export .errorProfile file, with the count of all error types and positions in the tags\n"
	    "//   -intergenicSupport -run run_name -mrnaRemap mrnaRemap.av_pg.txt\n"
	    "//       export .intergenicSupport, kb supporting exons, intron_regions, 2kb_UTR, intergenic, non-stranded\n"  
	    "//       remap file format: class, transcript, x1, x2, chrom, a1, a2, gene\n"
	    "//       example : KT_RefSeq   _LOC100128190   322     449     Y       13537431        13537304       X__LOC100128190\n"
	    "//   [-geneSupport | -mrnaSupport [-hierarchic]] -run run_name [-stranded | -antistranded]:\n"
	    "//       export .geneSupport file, gene run supporting reads and aligned bp\n"
            "//      -unique_gene_support fileName: expects a .ace file previouly exported by option -geneSupport2ace\n"
	    "//       This file contains the count per gene of the uniquely mapped reads\n"
	    "//         the quasi-unique mapped reads will be attributed in damped proportions to the same genes\n"
	    "//       If this file is given all reads are counted, otherwise only the unique reads are counted\n"
	    "//    -maxErr n -maxErrRate x : reaed with nerr > n AND 100*nerr > r * ali are not counted as gene/mrnaSupportn"
	    "//   -geneSupport2ace \n"
	    "//   -mrnaSupport2ace [-selected8kbList f8kb -selected5kbList f5kb]\n"
	    "//       Read a 10 column tables target_class,gene|mrna,u|nu,run,seq,tag,bp support,err,partial,orphan per tissue\n"
	    "//       exports a .ace file containing a dummy gene index\n"
	    "//       for all mRNAs the ace file also contains the mRNA wiggle\n"
	    "//       for mRNAS in f8Kb/f5kb 1 column file, export the 3'->5' cumulated wiggle\n"
	    "//   -geneRemap fileName\n"
	    "//       METADATA/mrnaRemap.gz file, here used to characterize fragments bridging genes close in cis\n"
	    "//   -geneIndex2Table -db db -target_class name -project project_name\n"
	    "//       Export a table reporting for each gene/group/run in project the index/seq/tags/kb\n"
	    "//   -intronSupport -mrnaRemap mrnaRemap.av_pg.txt\n"
	    "//       Export for all introns in mrnaRemapFile the number of reads in the hits file spanning 8bp on each side\n"
	    "//   -intronSupport2ace \n"
	    "//       Read a 6 column tables intron_class,intron,u|nu,seq,tag,bp support per run\n"
	    "//       exports an an ace containing a dummy intron index\n"
	    "//   -autoH\n"
	    "//       Auto-hierarchy sort target by number of hits, then keep only hits to top target\n"
	    "//   -groupLetterProfile -db db\n"
	    "//       Add up the run letter profiles into the groups\n"
	    "//   -addQualityFactors -db db\n"
	    "//       Add up the quality factors in runs and groups\n"
	    "//   -fastqAnalysis fastq_file : requires -run name -hits hits_file\n"
	    "//       Compare a fastq file to a hist file\n"
	    "//       For the identifiers in common, count the exact/errors per positions and qualities\n"
	    "//       Export the corresponding histograms\n"
	    "//   -checkGroupHierarchy -db MetaDB -project $MAGIC\n"
	    "//       Check that the meta-database defines an acceptable group hierarhy\n"
	    "//       runs with tag sublibrary_of should all be assigne to a single run, not to a group\n"
	    "//       runs and groups can belong to groups, but the graph must be acyclic\n"
	    "//       starting from (non-group) level 0 runs, one move to up to larger groups, but never down\n"
	    "//\n"      
	    "// Filters\n"
	    "//   -pair  int : int is the average length of the pair\n"
	    "//       Use paired end coherence to select the best ali, export pair statistics\n"
	    "//   -minAli int : exclude shorter alignments\n"
	    "//   -maxHit int : exclude sequences with more than this number of best exon targets\n"
	    "//             (alternative transcripts of the same gene do not matter)\n"
	    "//   -target_class name : only consider hits to this target class\n"
	    "//   -unique : only consider unique hits (within a target class)\n"
	    "//   -pureNsStrand int\n"
	    "//        If 1 or 2, tags hitting targets on opposite strands are discarded\n"
	    "//        except that if 2, tags hitting exactly 2 targets on opposite strands are kept\n"
	    "//   -subsampling int\n"
	    "//        optional, if s>0 keep randomly 1 fragment out of every s fragments seen in the hits file\n" 
	    "//   -filter [filter_name | minAli_number]\n"
	    "//        apply a named quality filter, (worm | ce | ara | at | sc | any other name treated as human)\n"
	    "//        or give a number usually equal to minAli if you use minAli < 24\n"
            "// -strategy [ Exome | Genome | RNA_seq ]\n"
	    "//        optional, should match the choice in clipalign\n"
	    "// -lane lane_name\n"
	    "//        optional, used in the -seqc .seqc.tsv format\n"
	    "// Help\n"
	    "//   -help : export this on line help\n"
	    "// Caveat:\n"
	    "//   Lines starting with '#' in the input files are considered as comments and dropped out\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/
/*************************************************************************************/
/* greg is a cigar string for blast_mapper looking like 20AC12GT20-A-T-T100G-T-A-30
 * transform it into magic type x and a lists of variations
 * The coordinates of the alignemnt of greg x1,x2 -> a1,a2 follows the magic conventions
 */
static int greg2aceCigarette (int x1, int x2, int a1, int a2, const char *greg, vTXT xSnps, vTXT aSnps) 
{
  const char *cp ;
  char cc, cc1 ;
  int strand , k = 0 ;
  int step ;
  int aa = a1 ; 
  int xx = x1 ;
  int nerr = 0 ;

  vtxtClear (xSnps) ;
  vtxtClear (aSnps) ;

  if (x1 > x2) messcrash ("greg2aceCigars expects x1=%d < x2=%d", x1, x2) ;
  strand = a1 <= a2 ? 1 : -1 ;

  cp = greg - 1 ;
  while (*++cp)
    { 
      step = 1 ;
      if (*cp >= '0' && *cp <= '9')
	{
	  step = 0 ;
	  while (*cp >= '0' && *cp <= '9')
	    step = 10 * step + (*cp++ - '0') ;
	  cp-- ; /* retrograde to last parsed digit */
	  xx += step ; aa += step * strand ;
	  continue ;
	}
       if (*cp == '-')
	{
	  char aBuf[8], xBuf[8], zBuf[8] ;
	  int j = 0 ;
	  
	  cc = ace_lower (*++cp) ;
	  xBuf[j] = cc ; zBuf[j] = '-' ; nerr ++ ;
	  if (strand < 0) cc = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)cc])] ; 
	  aBuf[j++] = cc ; 
	  
	  while (j < 3 && *(cp+1) == '-')
	    {
	      cp += 2 ;  cc = ace_lower (*cp) ; xBuf[j] = cc ; zBuf[j] = '-' ;  
	      if (strand < 0) cc = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)cc])] ;
	      aBuf[j++] = cc ;
	    }
	  xBuf[j] = 0 ; aBuf[j] = 0 ; zBuf[j] = 0 ; aa += j * strand ;
	  if (strand == -1)
	    {
	      int i ;
	      for (i = 0 ; i < j/2 ; i++)
		{ cc = aBuf[i] ; aBuf[i] = aBuf[j - 1 - i] ; aBuf[j - 1 - i] = cc ; }
	    }
	  vtxtPrintf (xSnps, "%s%d:", k ? "," : "", xx) ;
	  vtxtPrintf (aSnps, "%s%d:", k ? "," : "", aa + (strand == 1 ? -j : 1)) ;
	  k++ ;
	  vtxtPrintf (xSnps, "%s%s", zBuf, xBuf) ;
	  vtxtPrintf (aSnps, "%s%s", zBuf, aBuf) ;
	  
	}
      else
	{
	  cc = ace_lower (*cp++) ; cc1 = ace_lower (*cp) ;
	  if (cc1 == '-')
	    {
	      char aBuf[8], xBuf[8], zBuf[8] ;
	      int j = 0 ;
	  
	       xBuf[j] = cc ; zBuf[j] = '+' ; nerr++ ;
	       if (strand < 0) cc = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)cc])] ;
	       aBuf[j++] = cc ;
	      while (*(cp+2) == '-')
		{
		  cp++ ;  cc = ace_lower (*cp++) ; xBuf[j] = cc ; zBuf[j] = '+' ; 
		  if (strand < 0) cc = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)cc])] ;
		  aBuf[j++] = cc ;
		}
	      xBuf[j] = 0 ; aBuf[j] = 0 ; zBuf[j] = 0 ; xx += j ;
	      if (strand == -1)
		{
		  int i ;
		  for (i = 0 ; i < j/2 ; i++)
		    { cc = aBuf[i] ; aBuf[i] = aBuf[j - 1 - i] ; aBuf[j - 1 - i] = cc ; }
		}
	      vtxtPrintf (xSnps, "%s%d:", k ? "," : "", xx + (strand == 1 ? -j : -2)) ;
	      vtxtPrintf (aSnps, "%s%d:", k ? "," : "", aa + (strand == 1 ? 0 : 1)) ;
	      k++ ;
	      vtxtPrintf (xSnps, "%s%s", zBuf, xBuf) ;
	      vtxtPrintf (aSnps, "%s%s", zBuf, aBuf) ;
	    } 
	  else
	    {
	      vtxtPrintf (xSnps, "%s%d:", k ? "," : "", xx) ;
	      vtxtPrintf (aSnps, "%s%d:", k ? "," : "", aa) ;
	      k++ ;
	      vtxtPrintf (xSnps, "%c>%c", cc, cc1) ; 
	      if (strand < 0) cc = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)cc])] ;
	      if (strand < 0) cc1 = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)cc1])] ;
	      vtxtPrintf (aSnps, "%c>%c", cc, cc1) ; xx++ ; aa += strand ; nerr++ ;
	    }
	}
    }
  if (0)
    {
      printf ("xSnps :: %s\n", vtxtPtr (xSnps) ? vtxtPtr (xSnps) : "") ;
      printf ("aSnps :: %s\n", vtxtPtr (aSnps) ? vtxtPtr (aSnps) : "") ;
    }
  xx-- ; aa -= strand ;
  if (xx != x2 || aa != a2)
    fprintf (stderr, "ERROR greg2aceCigarette x1=%d x=%d x2=%d a1=%d a=%d a2=%d strand =%d :: %s\n", x1, xx, x2, a1, aa, a2,strand, greg) ;
  return nerr ;
}  /* greg2aceCigarette */

/*************************************************************************************/

static void baGreg2ace (BA *ba) 
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (ba->hitFileName, ba->gzi, h) ;
  ACEOUT ao = ba->outFileName ? aceOutCreate (ba->outFileName, ".hits", FALSE, h) : aceOutCreateToStdout (h) ;
  vTXT xSnps = vtxtHandleCreate (h) ;
  vTXT aSnps = vtxtHandleCreate (h) ;
  char *cp ;
  int mult, score, a1, a2, a3, da, x1, x2, dx, ln, targetMult, compartment, nerr ;
  char readName[101], targetName[101], over5p[101], over3p[101] ;
  BOOL isReadDown ;

  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      cp = aceInWord (ai) ; if (! cp || *cp == '#') continue ;
      strncpy (readName, cp, 99) ; readName[100]= 0 ;

      mult = fastcMultiplicity (cp, 0, 0) ;

      aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (! cp || *cp == '#') continue ;
      if (cp) { strncpy (targetName, cp, 99) ; targetName[100] = 0 ; } else targetName[0] = 0 ;

      aceInStep (ai, '\t') ; aceInWord (ai) ; /* per cent identity, drop */
      aceInStep (ai, '\t') ; aceInWord (ai) ; /* place holder, drop */
      aceInStep (ai, '\t') ; aceInWord (ai) ; /* place holder, drop */
      aceInStep (ai, '\t') ; aceInWord (ai) ; /* place holder, drop */

      aceInStep (ai, '\t') ; aceInInt (ai, &x1) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &x2) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &a1) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &a2) ;

      aceInStep (ai, '\t') ; aceInWord (ai) ; /* place holder, drop */
      aceInStep (ai, '\t') ; aceInWord (ai) ; /* place holder, drop */
      aceInStep (ai, '\t') ; aceInInt (ai, &score) ;
      
      aceInStep (ai, '\t') ; cp = aceInWord (ai) ; isReadDown = ! strcmp (cp, "minus") ? FALSE : TRUE ; 
      aceInStep (ai, '\t') ; cp = aceInWord (ai) ; /* isTargetDown = ! strcmp (cp, "minus") ? FALSE : TRUE ; */

      aceInStep (ai, '\t') ; aceInInt (ai, &ln) ;
      aceInStep (ai, '\t') ; cp = aceInWord (ai) ;  nerr = greg2aceCigarette (x1, x2, a1, a2, cp, xSnps, aSnps) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &targetMult) ; targetMult = 1 ; /* ERROR this col is broken */
      aceInStep (ai, '\t') ; aceInWord (ai) ; /* polyA exon.. */

      aceInStep (ai, '\t') ; aceInInt (ai, &compartment) ;

      aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (cp) { strncpy (over5p, cp, 99) ; over5p[100] = 0 ; }   else over5p[0] = 0 ;
      aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (cp) { strncpy (over3p, cp, 99) ; over3p[100] = 0 ; }   else over3p[0] = 0 ;

      aceInStep (ai, '\t') ; aceInWord (ai) ; /* place holder, drop */

      a3 = 0 ; aceInStep (ai, '\t') ; aceInInt (ai, &a3) ; /* pair pos */ if (a3 && a1 < a2) da = a3 - a1 ;  if (a3 && a1 > a2) da = a1 - a3 ; 

      /* success, export that line */
      if (! isReadDown) { int x = x1 ; x1 = x2 ; x2 = x ; x = a1 ; a1 = a2 ; a2 = x ; }
      dx = x2 - x1 + 1 ; 
      if (ba->pair)
	{
	  char *cp = readName + strlen (readName) - 2 ;
	  if (cp[0] == '.' && cp[1] == '1') { cp[0] = '>' ; cp[1] = 0 ; }
	  if (cp[0] == '.' && cp[1] == '2') { cp[0] = '<' ; cp[1] = 0 ; }
	}
      /* score = 200 + compartment ; // do this systematically untill we use the compartment correctly */
      aceOutf (ao, "%s\t%d\t%d"
		 , readName
		 , score, mult
	       ) ;
      aceOutf (ao, "\t%d\t%d\t%d\t%d\t%s\t%s"
		 , ln, dx, x1, x2
	       , ba->target_class ? dictName (ba->target_classDict, ba->target_class) : "Z_genome"
	       , "-"
	       ) ;
      aceOutf (ao, "\t%d\t%s\t%d\t%d\t%d\t%d"
	       , targetMult
	       , targetName, a1, a2
	       , 0, nerr
	       ) ;
      aceOutf (ao, "\t%s\t%s\t%s\t%s"
	       , vtxtPtr (xSnps) ? vtxtPtr (xSnps) : "-"
	       , vtxtPtr (aSnps) ? vtxtPtr (aSnps) : "-" 
	       , over5p
	       , over3p
	       ) ;
       aceOutf (ao, "\t-\t-\t%d"
		, da
		) ;      
      aceOutf (ao, "\n") ; 
    }

  ac_free (h) ;
} /* greg2ace */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  BA ba ;
  const char *errors = 0 ;
      
  AC_HANDLE h = ac_new_handle () ;
  const char *ccp, *dbName = 0 ;

  if (1)   freeinit () ; 
  messErrorInit (argv[0]) ;

  memset (&ba, 0, sizeof (BA)) ;
  ba.h = h ;

  if (argc == 1)
    usage (0) ;

  if (!strcmp (argv[1],"-test"))
    {
      char *txt = "select g,r from g in class \"run\", r in g->union_of" ;
      AC_TABLE tbl ;
      AC_DB db = ac_open_db ("../MagicMiniTest/MetaDB", 0) ; 
      vTXT v = vtxtHandleCreate (h) ;
      const char *errors = 0 ;

      tbl = ac_bql_table (db, txt, 0, 0, 0, h) ;
      
      if (1)
	{
	  txt = "select a,b,x,z from a in  class ali where a == \"subfr\", b in a->gene_expression , x in b[1], z in x[1]" ;
	  tbl = ac_bql_table (db, txt, 0, 0, &errors, h) ;
	  if (tbl)
	    ac_table_display (v, tbl, 0, 0, 0, 0, 0, 0, 0) ;
	  else if (*errors)
	    fprintf (stderr, "%s\n", errors) ;
	  if (vtxtPtr (v))
	    printf ("%s\n", vtxtPtr (v)) ;
	  printf ("nothing to do \n") ;
	}

      if (0)
	{
	  txt = "select a,b,x,z from a in  class ali, b in a->gene_expression , x in b[1],z=x" ;
	  tbl = ac_bql_table (db, txt, 0, 0, &errors, h) ;
	  if (tbl)
	    ac_table_display (v, tbl, 0, 0, 0, 0, 0, 0, 0) ;
	  else if (*errors)
	    fprintf (stderr, "%s\n", errors) ;
	  if (vtxtPtr (v))
	    printf ("%s\n", vtxtPtr (v)) ;
	  printf ("nothing to do \n") ;
	}
      exit (0) ;
    }


  ba.runDict = dictHandleCreate (100, h) ;

  /* optional arguments */

  getCmdLineOption (&argc, argv, "-db", &dbName) ;
  getCmdLineOption (&argc, argv, "-i", &(ba.hitFileName)) ;
  getCmdLineOption (&argc, argv, "-o", &(ba.outFileName)) ;
  getCmdLineOption (&argc, argv, "-mrnaRemap", &(ba.mrnaRemapFileName)) ;
  getCmdLineOption (&argc, argv, "-unique_gene_support", &(ba.uniqueGeneSupportFileName)) ; 
  getCmdLineOption (&argc, argv, "-run", &(ba.run)) ;
  getCmdLineOption (&argc, argv, "-title", &(ba.title)) ;
  getCmdLineInt (&argc, argv, "-maxHit", &(ba.maxHit)) ;
  getCmdLineInt (&argc, argv, "-minAli", &(ba.minAli)) ;
  ba.minIntronOverlap = 8 ;
  getCmdLineInt (&argc, argv, "-minIntronOverlap", &(ba.minIntronOverlap)) ;
  getCmdLineInt (&argc, argv, "-ventilate", &(ba.ventilate)) ;
  getCmdLineOption (&argc, argv, "-id_file", &(ba.id_filename)) ;
  getCmdLineOption (&argc, argv, "-inFileList", &(ba.inFileList)) ;
  getCmdLineOption (&argc, argv, "-geneRemap", &(ba.geneRemapFile)) ;
  getCmdLineOption (&argc, argv, "-sigTargets", &(ba.sigTargetFile)) ;
  getCmdLineOption (&argc, argv, "-split_mRNAs", &(ba.splitMrnaFileName)) ;
  ba.gzi = getCmdLineOption (&argc, argv, "-gzi", 0) ;
  ba.gzo = getCmdLineOption (&argc, argv, "-gzo", 0) ;
  ba.unique = getCmdLineOption (&argc, argv, "-unique", 0) ;
  ba.greg2ace = getCmdLineOption (&argc, argv, "-greg2ace", 0) ;
  ba.autoH = getCmdLineOption (&argc, argv, "-autoH", 0) ;
  getCmdLineInt (&argc, argv, "-pair", &(ba.pair)) ;
  getCmdLineInt (&argc, argv, "-Remove_inserts_shorter_than", &(ba.Remove_inserts_shorter_than)) ;
  getCmdLineOption (&argc, argv, "-filter", &ba.filter) ; 
  getCmdLineOption (&argc, argv, "-fastqAnalysis", &ba.fastqAnalysis) ;
  ba.lane = "xxx" ;
  getCmdLineOption (&argc, argv, "-lane", &ba.lane) ;
  ba.countBest = getCmdLineOption (&argc, argv, "-countBest", 0) ;
  ba.exportBest = getCmdLineOption (&argc, argv, "-exportBest", 0) ;
  ba.exportVenn = getCmdLineOption (&argc, argv, "-exportVenn", 0) ;
  ba.seqc = getCmdLineOption (&argc, argv, "-seqc", 0) ;
  if (ba.seqc) ba.countBest = TRUE ;
  ba.exportSuffix = getCmdLineOption (&argc, argv, "-exportSuffix", 0) ;
  ba.keepTargetPrefix = getCmdLineOption (&argc, argv, "-keepTargetPrefix", 0) ;
  ba.exportMito = getCmdLineOption (&argc, argv, "-exportMito", 0) ;
  ba.errorProfile = getCmdLineOption (&argc, argv, "-errorProfile", 0) ;
  ba.aliProfile = getCmdLineOption (&argc, argv, "-aliProfile", 0) ;
  ba.intronSupport = getCmdLineOption (&argc, argv, "-intronSupport", 0) ;
  ba.intronSupport2ace = getCmdLineOption (&argc, argv, "-intronSupport2ace", 0) ;
  ba.intergenicSupport = getCmdLineOption (&argc, argv, "-intergenicSupport", 0) ;
  ba.geneSupport = getCmdLineOption (&argc, argv, "-geneSupport", 0) ;
  ba.hierarchic = getCmdLineOption (&argc, argv, "-hierarchic", 0) ;
  ba.mrnaSupport = getCmdLineOption (&argc, argv, "-mrnaSupport", 0) ;
  if (ba.mrnaSupport) ba.geneSupport = TRUE ;
  ba.maxErr = ba.maxErrRate = 0 ; /* defaults */
  getCmdLineInt(&argc, argv, "-maxErr", &(ba.maxErr)) ;
  getCmdLineInt(&argc, argv, "-maxErrRate", &(ba.maxErrRate)) ;
  getCmdLineInt(&argc, argv, "-subsampling", &(ba.subsampling)) ;

  ba.geneSupport2ace = getCmdLineOption (&argc, argv, "-geneSupport2ace", 0) ;
  ba.mrnaSupport2ace = getCmdLineOption (&argc, argv, "-mrnaSupport2ace", 0) ;
  ba.geneIndex2Table = getCmdLineOption (&argc, argv, "-geneIndex2Table", 0) ;
  ba.groupLetterProfile = getCmdLineOption (&argc, argv, "-groupLetterProfile", 0) ;
  ba.addQualityFactors = getCmdLineOption (&argc, argv, "-addQualityFactors", 0) ;
  ba.stranded = 0 ;
  if (getCmdLineOption (&argc, argv, "-stranded", 0))
    ba.stranded = 1 ;
  if (getCmdLineOption (&argc, argv, "-antistranded", 0))
    ba.stranded = -1 ;
  
  getCmdLineOption (&argc, argv, "-selected8kbList", &ba.selected8kbFileName) ;
  getCmdLineOption (&argc, argv, "-selected5kbList", &ba.selected5kbFileName) ;

  getCmdLineInt (&argc, argv, "-pureNsStrand", &ba.pureNsStrand) ;
  getCmdLineOption (&argc, argv, "-project" , &ba.project) ;
  if (getCmdLineOption (&argc, argv, "-strategy", &ccp))
    {
      if (! strcasecmp (ccp, "Exome") ||
	  ! strcasecmp (ccp, "Genome")
	  )
	{
	  ba.errCost = 4 ;
	}
      else if (! strcasecmp (ccp, "RNA_seq"))
	{
	  ba.errCost = 8 ; ba.RNA_seq = TRUE ;
	}
    }
  ba.target_classDict = dictHandleCreate (10000, ba.h) ;
    {
        /* synchronize with bin/target2target_class.txt with c5.h_Ali.awk and with baExportAliProfile() */
      const char **cl ;
      const char *classes[] = { "0_SpikeIn", "1_DNASpikeIn", "A_mito", "B_rrna", "C_chloro", "D_transposon"
				, "DT_magic", "ET_av", "FT_av2008", "FT_extra", "FT_simul", "FT_cloud", "FT_Eco", "FT_Tjaden", "KT_RefSeq", "LT_RefSeqCurrent",  "LT_seqc", "LT_magic", "LT_UCSC", "MT_EBI", "MT_Gaj", "NT_miRNA", "NT_MiT",  "NT_HINV", "NT_FBK", "NT_FlyBase", "OT_rnaGene", "PT_tRNA", "QT_smallRNA"
				, "S_est", "U_introns", "W_Line", "X_Bami", "Y_Pfluo", "Z_genome", "b_bacteria", "v_virus", "z_gdecoy"
			       , 0 } ;
      for (cl = classes ; *cl ; cl++)
	{
	  int k = 0 ;
	  dictAdd (ba.target_classDict, *cl, &k) ;
	  if ((*cl)[1] == 'T' && (*cl)[2] == '_')
	    ba.lastTranscriptClass = k ;      
	}
    }
  if (getCmdLineOption (&argc, argv, "-target_class" , &ccp))
    {
      dictAdd (ba.target_classDict, ccp, &(ba.target_class)) ;
    }

  if (getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  ba.checkGroupHierearchy = getCmdLineBool (&argc, argv, "-checkGroupHierarchy") ;
  if (argc > 1)
    {
      fprintf (stderr, "Unknown argument %s, try -help\n", argv[argc-1]) ;
      exit (1) ;
    }

  /* init global counters */
  badMaxHit = badLnS = badLnT = isLongS = isLongT = badScoreLS = badScoreLT = isPartialS = isPartialT = badScorePS = badScorePT = 0 ;

  if (ba.splitMrnaFileName)
    baSplitMrnaParse (&ba) ;

  if (ba.selected8kbFileName)
    {
      ba.selected8kbDict = dictHandleCreate (100, ba.h) ;
      ba.selected8kbWiggle = keySetHandleCreate (ba.h) ;
      baParseSelected8kbTargetFile (&ba) ;
    }
  if (ba.selected5kbFileName)
    {
      ba.selected5kbDict = dictHandleCreate (100, ba.h) ;
      ba.selected5kbWiggle = keySetHandleCreate (ba.h) ;
      baParseSelected5kbTargetFile (&ba) ;
    }

  /* actions */

  fprintf (stderr, "Start %s\n", timeShowNow ()) ;

  if (dbName)
    {
      ba.db = ac_open_db (dbName, &errors);
      if (! ba.db)
	messcrash ("Failed to open db %s, error %s", dbName, errors) ;
    }

  if (ba.seqc)
    {
      ba.seqco = aceOutCreate (ba.outFileName, ".seqc.tsv", FALSE, ba.h) ;
    }

  if (ba.geneSupport)
    { 
      if (ba.exportBest)
	messcrash ("-exportBest and -geneSupport option are incompatible") ;
      if (! ba.run)
	messcrash ("-geneSupport requires option -run run_name") ;
      if (! ba.target_class)
	messcrash ("-geneSupport requires option -target_class to correctly select the genes") ;
      if (ba.uniqueGeneSupportFileName)  /* parse the unique support to allow bayesian ventilation of the multi hits */
	baUniqueGeneSupportParse (&ba) ;
    }
  if (ba.intronSupport)
    {
      if (! ba.run)
	messcrash ("-intronSupport requires option -run run_name") ;
      /*
	if (! baCreateAtlas ())
	messcrash ("-intronSupport requires the construction of a remapping atlas") ;
      */
    }
  if (ba.intergenicSupport && ! baCreateAtlas (&ba, 0))
    ba.intergenicSupport = FALSE ;

  if (ba.greg2ace)
    {
      baGreg2ace (&ba) ;
    }
  else if (ba.checkGroupHierearchy)
    {
      if (! ba.db)
	messcrash ("Option -checkGroupHierarchy requires option -db pointing to the MetaDB database") ;
      if (! ba.project)
	messcrash ("Option -checkGroupHierarchy requires option -project project") ;
      baCheckGroupHierearchy (&ba) ;
    }
  else if (ba.geneSupport2ace)
    {
      baElementSupport2ace (&ba, TRUE, FALSE, FALSE) ;
    }
  else if (ba.mrnaSupport2ace)
    {
      baElementSupport2ace (&ba, FALSE, TRUE, FALSE) ;
    }
  else if (ba.intronSupport2ace)
    {
      baElementSupport2ace (&ba, FALSE, FALSE, TRUE) ;
    }
  else if (ba.geneIndex2Table)
    {
      if (! dbName)
	messcrash ("-geneIndex2Table requires -db dbName") ;
      if (! ba.db)
	messcrash ("-geneIndex2Table requires -db dbName, cannot open the acedb database, sorry: %s", errors ? errors : "") ;
      if (! ba.project)
	messcrash ("-geneIndex2Table requires -project project") ;
      if (! ba.target_class)
	messcrash ("-geneIndex2Table requires -target_class") ;
      baGeneIndex2Table (&ba, TRUE) ;
      baGeneIndex2Table (&ba, FALSE) ;
    }
  else if (ba.groupLetterProfile)
    {
      if (! dbName)
	messcrash ("-groupLetterProfile requires -db dbName") ;
      if (! ba.db)
	messcrash ("-groupLetterProfile requires -db dbName, cannot open the acedb database, sorry: %s", errors ? errors : "") ;
      baGroupLetterProfile (&ba) ;  /* 2016_03_24: recurse through group levels; 2014_08_24, use a sigle pass{ >runs} SETOR {>subgroup;>runs} */
      if (0) baSetStrand (&ba) ; /* obsolete */
      if (ba.addQualityFactors)
	baAddQualityFactors (&ba) ;
   }
  else if (ba.addQualityFactors)
    {
      if (! dbName)
	messcrash ("-addQualityFactors requires -db dbName") ;
      if (! ba.db)
	messcrash ("-addQualityFactors requires -db dbName, cannot open the acedb database, sorry: %s", errors ? errors : "") ;
      baAddQualityFactors (&ba) ;
    }
  else if (ba.ventilate)
    {
      if (! ba.id_filename)
	messcrash ("-ventilate requires -id_file filename") ;
      if (ba.ventilate <= 0)
	messcrash ("-ventilate n requires a positive value, you specified %d < 0", ba.ventilate) ;
      baVentilateHitFilePerOriginalIdentifier (&ba) ;
    }
  else if (ba.fastqAnalysis)
    {
      if (! ba.run)
	messcrash ("-fastqAnalysis requires option -run run_name") ;
      baFastqAnalysis (&ba) ;
    }
  else  if (ba.autoH) /* Auto-hierarchy sort target by number of hits, keep only hits to top target */
    baAutoHierarchy (&ba) ; 
  else
    {
      int filtered = 0 ;

      baParseHits (&ba) ;  /* rejects on the fly overlapping fragments with lower score */
      if (1 &&
	  ba.pair && arrayMax (ba.hits))
	filtered += baPairFilter (&ba) ;
      if (1 || ba.filter)
	filtered += baFilter (&ba) ;  /* filter the chains independently of each other */
      if (ba.geneRemapFile)
	baParseGeneRemapFile (&ba) ;
      if (filtered && ba.gene2map)
	baCheckUnicity (&ba) ;   /* double check that the mult=-2 cases correspond to pairs of genes in antisense */
      if (ba.intergenicSupport)
	baExportIntergenicSupport (&ba) ; /* must come first after parse to clean up the Z_genome */
      if (ba.countBest)
	baCountBest (&ba) ;
      if (ba.exportBest)
	baExportBest (&ba) ;
      if (ba.sigTargetFile)
	baExportSigHits (&ba) ; /* export a separate file with just the hits to this set */
      if (ba.exportVenn)
	baExportVenn (&ba) ;
      if (ba.exportSuffix)
	baExportSuffix (&ba) ;
      if (ba.exportMito)
	baExportMito (&ba) ;
      if (ba.errorProfile)
	baExportErrorProfile (&ba) ;
      if (ba.aliProfile)
	baExportAliProfile (&ba) ;
      if (ba.geneSupport)
	{
	  if (ba.uniqueGeneSupportFileName)
	    baExportGeneSupport (&ba, FALSE) ;
	  else
	    baExportGeneSupport (&ba, TRUE) ;
	}
       if (ba.intronSupport)
	{
	  dictAdd (ba.target_classDict,  "U_introns", &(ba.intronClass)) ;
	  arraySort (ba.hits, intronHitOrder) ;
	  if (baCreateAtlas (&ba, 1) || 1)
	    {
	      baIntronSupport (&ba, FALSE) ;
	      baIntronSupport (&ba, TRUE) ;
	    }
	}
    }
     
  /* clean up */
  fprintf (stderr, "// done: %s\n", timeShowNow()) ;
  {
    int mx ;
    messAllocMaxStatus (&mx) ;   fprintf (stderr, "// max memory %d Mb\n", mx) ;
  }
  
  /* report the global counters */
  if (ba.filter)
    fprintf(stderr,"filter=%s\nbadMaxHit=%d\nbadLnS=%d badLnT=%d isLongS=%d isLongT=%d\nbadScoreLS=%d badScoreLT=%d\nisPartialS=%d isPartialT=%d badScorePS=%d badScorePT=%d\n", ba.filter, badMaxHit,badLnS, badLnT, isLongS, isLongT, badScoreLS, badScoreLT, isPartialS, isPartialT, badScorePS, badScorePT) ;

  showGS2I (0) ; /* for comiler happiness */

  ac_free (ba.db) ;
  fprintf (stderr, "Done %s\n", timeShowNow ()) ;
  ac_free (h) ;

  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

