/*
 * Analyse the solexa alignments generated from Bowtie
 * Analyse the letter distribution in the solexa fasta file
 * Analyse the pattern of errors
 * Search for the vector
 */


#include "../wac/ac.h"
#include "keyset.h"
#include "dna.h"
#include <errno.h>
#include <math.h>
#include "keyset.h"
#include "freeout.h"
#include "vtxt.h"
#include "dict.h"
#include "aceio.h"
#include "wiggle.h"

typedef enum {ZERO, INTERGENIC, INTRON, EXON, CODINGEXON } zTYPE ;
char *zType[] = {"ZERO", "Intergenic", "Intronic", "Exonic non-coding", "Exonic coding" }  ;

typedef struct zoneStruct { 
  AC_HANDLE h ; const char *title ; 
  const char *manip, *tissue, *lane, *run, *group, *target, *nmWiggle, *support, *cosmidWiggle,
    *trackName, *hitFile, *wiggleDir, *wiggleFeature, *ignoredProbeFileName
    , *wiggleFileName, *wiggleFileName_ns, *wiggleFileName_f, *wiggleFileName_r, *mrnaRemapFileName
    , *sxxNewIntronsFileName, *sxxDeUnoIntronsFileName, *sxxKnownIntronsFileName, *sxxChromosome, *spongeFileName, *bestHitFileName
    , *spongeWindowExcludeFileName, *genomeMaskFileName, *sdfFileName
    ;
  int prefix, prefixq, delta, period, minX, maxX, saturation, sponge, intronSponge, transsplicing ;  
  int minimalSupport, minimalIntronSupport, seedLength, overhangLength, intronMinLength, intronMaxLength, exonLevel ; /* thresholds when constructing new introns */
  int non_gt_ag, width ;
  int maxChromNameLength ;
  KEYSET mrna2intron ;
  KEYSET splitMrnas ;
  DICT *ignoredProbeDict, *mapDict, *remapDict ;
  Array chromDnaD, chromDnaR, introns ;
  const char *outFileName ;
  const char *splitMrnasFileName ;
  const char *strategy ;
  BOOL RNA_seq ;
  BOOL plot, wordProfile, venn, intron, count, exp, cmp, nPlatform, hammer, stranded, mrnaWiggle, hits2composite
    , strand, antistrand, intronSpongeMerge,gs2index, unique, exomeDosage, geneDosage, variantDosage, spongeWindow
    , ventilate, newIntrons, newDoubleIntrons, newExons, flagTranscriptEnds, gzi, gzo ;} SX ;

typedef struct pfStruct { int probe, nerr, nN, nhit, targetHit[6], targetErr[5], isUp ; } PF ;
typedef struct pf2Struct { int prefix, nn, ali, ne[12] ; } PF2 ;

static  void sxParseHitsFiles (SX *sx, int lane, char **targets, Array aa, DICT *dict, BOOL isU, BOOL justBestLayer) ;
static BOOL sxTrackLocate (DICT *dict
			   , const char *nam1, char type1
			   , const char *nam2, char type2
			   , const char *nam3, char type3
			   , int *t1p, int *t2p, int *t3p) ;

/*************************************************************************************/
/*************************************************************************************/
/* plot a distribution of (position, words) i.e. sequencing errors
   stdin example : 
     12  A>T
     28  -T
*/
static  void sxWordProfile (SX *sx)
{
  DICT *dict = dictCreate (1000) ;
  Array aa = arrayCreate (100000, int) ;
  KEYSET histo = keySetCreate () ;
  KEYSET types = keySetCreate () ;
  ACEIN ai = aceInCreate (0, 0, sx->h) ; /* aceInCreate allows tab delimited, aceInCreateFrom masks \t into spaces */
  int i, n = 0, x, xmin, *xp, type, typeMin = 0 ;
  const char *ccp ;
  char *cp, buf[6] ;
  unsigned int imax = 0 ;

  dictAdd (dict, "--", &typeMin) ;

  dictAdd (dict, "A>T", &typeMin) ;
  dictAdd (dict, "A>G", &typeMin) ;
  dictAdd (dict, "A>C", &typeMin) ;

  dictAdd (dict, "T>A", &typeMin) ;
  dictAdd (dict, "T>G", &typeMin) ;
  dictAdd (dict, "T>C", &typeMin) ;

  dictAdd (dict, "G>A", &typeMin) ;
  dictAdd (dict, "G>T", &typeMin) ;
  dictAdd (dict, "G>C", &typeMin) ;

  dictAdd (dict, "C>A", &typeMin) ;
  dictAdd (dict, "C>T", &typeMin) ;
  dictAdd (dict, "C>G", &typeMin) ;


  dictAdd (dict, "Insert A", &typeMin) ;
  dictAdd (dict, "Insert T", &typeMin) ;
  dictAdd (dict, "Insert G", &typeMin) ;
  dictAdd (dict, "Insert C", &typeMin) ;

  dictAdd (dict, "Delete A", &typeMin) ;
  dictAdd (dict, "Delete T", &typeMin) ;
  dictAdd (dict, "Delete G", &typeMin) ;
  dictAdd (dict, "Delete C", &typeMin) ;

  if (0)
    {
      dictAdd (dict, "Double deletion", &typeMin) ;
      dictAdd (dict, "Double insertion", &typeMin) ;
    }

  while (n++, aceInCard (ai))
    {
      if (aceInInt (ai, &x))
	{
	  aceInStep(ai, '\t') ;	    
	  if ((ccp = aceInWord (ai)))
	    {
	      strncpy (buf, ccp, 5) ;
	      for (cp = buf ; *cp ; cp++)
		*cp = ace_upper(*cp) ;
	      if (buf[0] == '-' && buf[1] == '-' && buf[2] == 0)
		continue ;
	      if (buf[0] == '+' && buf[1] == '+' && buf[2] == 0)
		continue ;
	      else if (buf[0] == '-' && buf[1] != 0 && buf[2] == 0)
		ccp = messprintf ("Delete %c", buf[1]) ;
	      else if (buf[0] == '+' && buf[1] != 0 && buf[2] == 0)
		ccp = messprintf ("Insert %c", buf[1]) ;
	      else if (buf[0] != 0 && buf[1] == '>' && buf[2] != 0 && buf[3] == 0)
		ccp = buf ; 
	      dictAdd (dict, ccp, &type) ;
	      keySet (types, type)++ ;
	      array (aa, imax++, int) = x ;
	    }
	}
    }
  if (!imax)
    messcrash ("There was no int number on stdin, sorry") ;
  arraySort (aa, intOrder) ;

  for (i = 0, xp = arrp (aa, 0, int), xmin = *xp ; i < imax ; i++, xp++)
    keySet (histo, *xp - xmin)++ ;

  freeOutf ("Type\tNumber of occurences\n") ;
  imax = keySetMax (types) ;
  for (i = 2 ; i < imax ; i++)
    {
      if (i <= typeMin || keySet(types,i) > 0)
	freeOutf ("%s\t%d", dictName(dict, i), keySet(types, i)) ;
      freeOutf ("\n") ;
    }

  freeOutf ("\nPosition\tNumber of reads\n") ;
  imax = keySetMax (histo) ;
  for (i = 0 ; i < imax ; i++)
    {
      freeOutf ("%d\t%d", xmin + i, keySet (histo, i)) ;
      freeOutf ("\n") ;
    }

  ac_free (dict) ;
  keySetDestroy (histo) ;
  keySetDestroy (types) ;
} /* sxErrorProfile */

/*************************************************************************************/
/*************************************************************************************/

/*************************************************************************************/
/*************************************************************************************/
typedef struct h2cStruct {int target, a1, a2, probe, seq ;} H2C ;
static int H2Corder (const void *a, const void *b)
{
  const H2C *up = (const H2C *) a, *vp = (const H2C *) b ;
  if (up->target < vp->target) return -1 ;
  if (up->target > vp->target) return 1 ;
  
  return up->a1 - vp->a1 ;
}

static void showH2c (DICT *dict, Array hits)
{
  H2C *up ;
  unsigned int i = hits ? arrayMax (hits) : 0 ;

  if (i)
    for (i = 0 ; i < arrayMax (hits) ; i++)
      {
	up = arrp (hits, i, H2C) ;
	printf ("%s %d %d %s\n"
		, up->probe ? dictName(dict, up->target) : "null"
		, up->a1, up->a2
		, up->probe ? dictName(dict, up->probe) : "null"
		) ; 
      }
}
/*************************************************************************************/
/* Create a composite EST from  */
static void sxHits2composite (SX *sx)
{
  int dx = -16, kmin = 5, kk ;
  ACEIN ai = 0 ;
  char oldProbe[1000] ;
  char *fNam ;
  const char *ccp ;
  int target, nn, a1, a2, olda1 = 0, olda2 = 0, ii ;
  H2C *up, *vp, *pf ;
  Array hits = arrayCreate (10000, H2C) ;
  Array targets = arrayCreate (10000, H2C) ;
  DICT *dict = dictHandleCreate (10000, sx->h) ;
  Stack s = stackHandleCreate (10000, sx->h) ;

  /* parse the fasta */
  if (!sx->target)
    messcrash ("argument -t target_fasta_file is missing") ;
  fNam = filName (sx->target, 0, "r") ;
  if (!fNam)
    messcrash ("Cannot find the target fasta file to be grignotted") ;
  ai = aceInCreate (sx->target, 0, sx->h) ;
  if (!ai)
    messcrash ("Cannot open the target fasta file to be grignotted") ;

  ii = 0 ;
  while (ai && aceInCard(ai))
    {
      /* this is fasta, get the tag name */
      ccp = aceInWord(ai) ;
      if (!ccp || !*ccp)
	{ ac_free (ai) ; continue ; }
    lao:
      if (*ccp != '>')
	{ continue ; }
      nn = -1 ;
      dictAdd (dict, ccp+1, &nn) ;
      pf = arrayp (targets, nn, H2C) ;
      pf->probe = nn ;  
      pf->seq = stackMark (s) ;  
      /* get the sequence */
      while (aceInCard(ai))
	{
	  ccp = aceInWord(ai) ;
	  if (!ccp || !*ccp)
	    break ;
	  if (*ccp == '>')
	    { goto lao ;}
	  catText (s, ccp) ;
	  ii += strlen (ccp) ;
	}
    }
  fprintf (stderr, "parsed %d bp in the target fasta file", ii) ;
  ac_free (ai) ;

  /* read all the hits */
  ai = aceInCreate (sx->hitFile, 0, sx->h) ;
  ii = 0 ;
  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ; /* probe name */
      if (!ccp || *ccp == '#' || !strcmp (oldProbe, ccp))
	continue ;
      strncpy (oldProbe, ccp, 999) ;
      aceInWord (ai) ;
      aceInWord (ai) ;

      aceInStep(ai,'\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp || ! *ccp)
	continue ;
      dictAdd (dict, ccp, &target) ;
      aceInStep(ai,'\t') ;
      if (! aceInInt (ai, &a1)) continue ;
      aceInStep(ai,'\t') ;
      if (! aceInInt (ai, &a2)) continue ;
      
      up = arrayp (hits, ii++, H2C) ;
      up->target = target ;
      up->a1 = a1 ;
      up->a2 = a2 ;
    }
  if (0) showH2c (dict, hits) ;
  arraySort (hits, H2Corder) ;
  if (0) showH2c (dict, hits) ;
  dx = 0 ; kmin = 1 ;
  for (ii = kk = 0, vp = 0 ; ii < arrayMax (hits) ; vp = up, ii++)
    {
      up = arrayp (hits, ii, H2C) ;
      if (ii == 0 || up->target != vp->target || up->a1 > olda2 + dx)
	{
	  if (ii > 0 && kk > kmin)
	    freeOutf ("%s\t%d\t%d\n", dictName (dict, vp->target), olda1,olda2);
	  olda1 = up->a1 ;
	  olda2 = 0 ;
	  kk = 0 ; 
	}
      kk++ ;
      if (up->a2 > olda2)
	olda2 = up->a2 ;
    }
} /* sxHits2composite */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

/* count for each 6mer prefix the number of sequence and error rate */
static  void sxParseHitsFiles (SX *sx, int lane, char **targets, Array aa, DICT *dict, BOOL isU, BOOL justBestLayer)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai ;
  int i, nTarget, nn, nN, nError, offset, cumul, a1, a2, x1, x2, dummy ;
  char *fNam ;
  const char *ccp ;
  char **target, buf[1024] ;
  PF *pf ;

  messcrash ("obsolete, this function was written in 2008, prior ot the current format, if you need it rewrite it ") ;
  /* the ai file shoudl contain the union of all probealign alignment for a given track of probes */
  /* read the number of errors for the best alignment of each probe */

  for (nTarget = 0, target = targets ; *target ; nTarget++, target++)
    for (offset = 1 ; offset < (isU ? 2 : 128) ; offset++)
      {
	if (nTarget >= 5) messcrash ("too many targets > 5 in sxParseHitsFiles") ;
	if (isU)
	  fNam = hprintf (h, "PHITS_%s.u/%s.%s.%d.%s.u.hits", *target, sx->manip, sx->tissue, lane, *target) ;
	else
	  fNam = hprintf (h, "PHITS_%s/%s.%s.%d.%s.off%d.hits", *target, sx->manip, sx->tissue, lane, *target, offset) ;
	if (filName (fNam, 0, "r"))
	  {
	    cumul = 0 ;
	    ai = aceInCreate (fNam, 0, h) ;
	    if (ai)
	      {
		/* aceInSpecial (ai, "\n\t") ; */
		messcrash ("you need to reimplement aceInStep('\t') everywhere") ;
		while (aceInCard (ai))
		  {
		    if (isU && (!strcmp (*target, "main_aceview") || !strcmp (*target, "RefSeq") || !strcmp (*target, "HINV") || !strcmp (*target, "EBI")))
		      {
			aceInStep(ai,'\t') ;
			ccp = aceInWord(ai) ; 
			if (!ccp || !*ccp || *ccp == '#')
			  goto nextCard ;
			aceInStep(ai,'\t') ;
                        aceInInt (ai, &dummy) ; /* number of gene hits */
			aceInStep(ai,'\t') ;
			aceInWord(ai) ; /* probe name */
			aceInStep(ai,'\t') ;
                        aceInInt (ai, &dummy) ; /* score */
			aceInStep(ai,'\t') ;
			aceInWord(ai) ; /* gene Name */
		      }
		    /* target name is column 4 */
		    aceInStep(ai,'\t') ;
		    ccp = aceInWord(ai) ;
		    if (!ccp || !*ccp || *ccp == '#')
		      goto nextCard ;

		    aceInStep(ai,'\t') ;
		    if (!aceInInt (ai, &a1))
		      goto nextCard ;
		    aceInStep(ai,'\t') ;
		    if (!aceInInt (ai, &a2))
		      goto nextCard ;

		    /* probe name is next column */
		    aceInStep(ai,'\t') ;
		    ccp = aceInWord(ai) ;
		    if (!ccp || !*ccp || *ccp == '#')
		      goto nextCard ;
		    strncpy (buf, ccp, 1000) ;
		    /* coordinates on probe */

		    aceInStep(ai,'\t') ;
		    if (!aceInInt (ai, &x1))
		      goto nextCard ;
		    aceInStep(ai,'\t') ;
		    if (!aceInInt (ai, &x2))
		      goto nextCard ;

		    /* NN, nError */
		    aceInStep(ai,'\t') ;
		    if (!aceInInt (ai, &nN))
		      goto nextCard ;
		    aceInStep(ai,'\t') ;
		    if (!aceInInt (ai, &nError))
		      goto nextCard ;

		    dictAdd (dict, buf, &nn) ;
		    pf = arrayp (aa, nn, PF) ;
		    pf->probe = nn ;
		    cumul++ ;
		    /* we register the best hits, possibly only in best layer */
		    if (!pf->nerr)
		      {
			pf->targetErr[nTarget] = pf->nerr = nError + 1 ;
			pf->targetHit[nTarget] = pf->nhit = 1 ;
			pf->nN = nN ;
			if (a1 > a2) pf->isUp = 1 ;
			else pf->isUp = 0 ;
		      }
		    else if (justBestLayer)
		      {
			if (pf->targetHit[nTarget] > 0) /* we are in best layer */
			  {
			    if (nError + 1 < pf->nerr || (nError + 1 == pf->nerr && pf->nN > nN))		
			      {
				pf->targetErr[nTarget] = pf->nerr = nError + 1 ;
				pf->targetHit[nTarget] = pf->nhit = 1 ;
				pf->nN = nN ;
				if (a1 > a2) pf->isUp = 1 ;
				else pf->isUp = 0 ;
			      }
			    else if (nError + 1 == pf->nerr)
			      {
				pf->targetHit[nTarget]++ ;
				pf->nN = nN ;
				pf->nhit++ ;
				if (a1 > a2) pf->isUp = 1 ;
				else pf->isUp = 0 ;
			      }
			  }
		      }
		    else
		      {
			if (nError + 1 < pf->nerr || (nError + 1 == pf->nerr && pf->nN > nN))
			  {
			    for (i = 0 ; i < nTarget ; i++)
			      pf->targetErr[i] = pf->targetHit[i] = 0 ; /* destroy previous hits */
			    pf->targetErr[nTarget] = pf->nerr = nError + 1 ;
			    pf->targetHit[nTarget] = pf->nhit = 1 ;
			    pf->nN = nN ;
			    if (a1 > a2) pf->isUp = 1 ;
			    else pf->isUp = 0 ;
			  }
			else if (nError + 1 == pf->nerr)
			  {
			    pf->targetHit[nTarget]++ ;
			    pf->nhit++ ;
			    pf->nN = nN ;
			    if (a1 > a2) pf->isUp = 1 ;
			    else pf->isUp = 0 ;
			  }
		      }
		  nextCard:
		    continue ;
		  }
		ac_free (ai) ;
		fprintf (stderr, "Found %d alignments in %s\n", cumul, fNam) ;
	      } 
	  }
      }
} /* sxParseHitsFiles */

/*************************************************************************************/
/* count for each tag where it belongs, how many hits, how many errors */
static  void sxTargetCount (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  Array aa = 0 ;
  DICT *dict = dictHandleCreate (10000000, h) ;
  int i, ii, nTarget, NN=10 ;
  char **target, *targets[] = {"mito", "main_aceview", "intron", 0} ;
  PF *pf ;
  int isUp[NN] ;
  int targetHit[NN] ;
  int targetMultiHit[NN] ;
  int targetErr[NN*NN] ;

  memset (isUp, 0, sizeof(isUp)) ;
  memset (targetHit, 0, sizeof(targetHit)) ;
  memset (targetMultiHit, 0, sizeof(targetMultiHit)) ;
  memset (targetErr, 0, sizeof(targetErr)) ;

  dictAdd (dict, "toto", 0) ;
  /* obtain best hits from any layer */
  aa = arrayHandleCreate (10000000, PF, h) ;
  sxParseHitsFiles (sx, 1, targets, aa, dict, FALSE, TRUE) ;
  for (ii = 0 ; ii < arrayMax (aa) ; ii++)
    {
      pf = arrp (aa, ii, PF) ;
      for (i = 0 ; i < 5 ; i++)
	if (pf->targetHit[i] > 0) 
	  {
	    targetHit[i]++ ;
	    if (pf->targetHit[i] > 1) 
	      targetMultiHit[i]++ ;
	    targetErr[i*NN+pf->nerr-1]++ ;
	    if (pf->isUp)
	      isUp[i]++ ;
	  }
    }
  /* report */
  freeOutf ("#Target\tUniqueTag\tMultiTag\tForward\tReverse\tExact\t1\t2\t3\t4\t5\n") ;
  for (nTarget = 0, target = targets ; *target ; nTarget++, target++)
    {
      freeOutf ("%s", *target) ;
      freeOutf ("\t%d", targetHit[nTarget]) ;
      freeOutf ("\t%d", targetMultiHit[nTarget]) ;
      freeOutf ("\t%d", targetHit[nTarget] - isUp[nTarget]) ;
      freeOutf ("\t%d", isUp[nTarget]) ;
      for (i = 0 ; i < 6 ; i++)
	freeOutf ("\t%d", targetErr[nTarget*NN+i]++) ;
      freeOutf ("\n") ;
    }
  ac_free (h) ;
} /* sxTargetCount */

/*************************************************************************************/
/* count for each 6mer prefix the number of sequence and error rate */
static  void sxPrefixDistrib (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  Array aa = 0, aa2 = 0 ;
  DICT *dict = dictHandleCreate (10000000, h) ;
  DICT *dictP = dictHandleCreate (1000, h) ;
  ACEIN ai ;
  int prefixLn = sx->prefix ?  sx->prefix :  sx->prefixq ;
  int i, nn, prefix, cumul, cumulali, cumule[12] ;
  char *fNam ;
  char **target, *targets[] = {"main_aceview", 0} ; /* {"mito", 0, "main_aceview", "genome", "intron", 0} ; */
  const char *ccp ;
  char buf[256] ;
  PF *pf ;
  PF2 *pf2 ;
  BOOL ok = TRUE ;

  sxTargetCount (sx) ; return ;
  memset (cumule, 0, sizeof(cumule)) ;

  dictAdd (dict, "toto", 0) ;
  dictAdd (dictP, "toto", 0) ;
  /* the ai file shoudl contain the union of all probealign alignment for a given track of probes */
  /* read the number of errors for the best alignment of each probe */
  aa = arrayHandleCreate (10000000, PF, h) ;
  aa2 = arrayHandleCreate (256, PF2, h) ;

  /* obtain best hits from any layer */
  sxParseHitsFiles (sx, sx->lane ? atoi(sx->lane) : 1, targets, aa, dict, FALSE, FALSE) ;

  /* we now analyse the fasta file, given on stdin */
  fNam = hprintf (h, "Fasta/%s/%s.%s.fa", sx->manip, sx->tissue, sx->lane, target) ;
  if (!filName (fNam, 0, "r"))
    fNam = hprintf (h, "Fasta/%s/%s.%s.ccfa", sx->manip, sx->tissue, sx->lane, target) ;
  if (filName (fNam, 0, "r"))
    {
      ai = aceInCreate (fNam, 0, h) ;
      if (ai)
	{
	  while (ai && aceInCard(ai))
	    {
	      /* this is fasta, get the tag name */
	      ccp = aceInWord(ai) ;
	      if (!ccp || !*ccp)
		{ ac_free (ai) ; continue ; }
	    lao:
	      if (*ccp != '>')
		{ continue ; }
	      nn = -1 ;
	      if (sx->ignoredProbeDict && dictFind (sx->ignoredProbeDict, ccp+1, 0))
		ok = FALSE ;
	      dictAdd (dict, ccp+1, &nn) ;
	      pf = arrayp (aa, nn, PF) ;
	      pf->probe = nn ;  
	      
	      /* get the sequence */
	      aceInCard(ai) ;
	      ccp = aceInWord(ai) ;
	      if (!ccp || !*ccp)
		{ ac_free (ai) ; continue ; }
	      if (*ccp == '>')
		{ goto lao ; }
	      
	      /* we have a probe and its number of errors , we construct its primer */
	      if (ok)
		{
		  strncpy (buf, ccp, prefixLn) ; buf[prefixLn] = 0 ;
		  dictAdd (dictP, buf, &prefix) ;
		  pf2 = arrayp (aa2, prefix, PF2) ;
		  pf2->prefix = prefix ;
		  pf2->nn++ ;
		  if (pf->nerr > 11) pf->nerr = 11 ;
		  (pf2->ne[pf->nerr])++ ;
		  if(pf->nerr) pf2->ali++ ;
		}
	      ok = TRUE ;
	    }
	  ac_free (ai) ;
	}
    }
  
  /* report */
  freeOutf ("#Prefix\tTag\t%% aligned\tNot aligned\tExact\t1 error\t2\t3\t4\t5\n") ;
  for (prefix = cumul = cumulali = 0 ; prefix < arrayMax (aa2) ; prefix++)
    {
      pf2 = arrayp (aa2, prefix, PF2) ;
      if (pf2->prefix)
	{
	  cumul += pf2->nn ;
	  cumulali += pf2->ali ;
	  freeOutf ("%s\t%d", dictName (dictP, pf2->prefix), pf2->nn) ;
	  for (i = 0 ; i < 7 ; i++)
	    {
	      freeOutf ("\t%d", pf2->ne[i]) ;
	      cumule[i] += pf2->ne[i] ;
	    }
	  freeOutf ("\n") ;
	}
    }
  freeOutf ("Cumul\t%d", cumul) ;
  for (i = 0 ; i < 7 ; i++)
    freeOutf ("\t%d", cumule[i]) ;
  freeOutf ("\n") ;
  ac_free (h) ;

  return ;
} /* sxPrefixDistrib */

/*************************************************************************************/
/* hard coded for SEQ-C 2008_12_14 */
static  void sxCompareCoverage (SX *sx)
{
  unsigned int flag, flagNx = 0, flagNx2 = 0 ;
  int i, ii, jj, n, nn, nmax, nUnion ;
  char **type, *types[] = { "1.", "2.", "3.", "4.", "6.", "7.", 0} ;
  unsigned int *ip ;
  char *cp ;
  const char **tissue, *tissues[] = {"uhr", 0, "brain", 0} ;
  const char **target, **target2, *targets[] = {"genome", "genome", "aceview", "av_main", "av_cloud", "nm_nr", "xm_xr", "nx", "i454", "avi", 0, 0} ;
  const char **title, *titles[] = {"Genome unique", "Genome", "AceView", "Main genes", "Cloud genes", "NM/NR RefSeq", "XM/XR RefSeq", "NM/NR/XM/XR", "454 introns", "AceView introns", "Super transcriptome", 0} ;
  char buf[1024] ;
  ACEIN ai = 0 ;
  DICT *dict = 0 ;
  Array flags = arrayCreate (50000000, unsigned int) ;

  if (sx->tissue) /* override the predefined list */
    { tissues[0] = sx->tissue ;   tissues[1] = 0 ; }
  fprintf(stderr, "t0=%s len=%d\n", tissues[0], (int)strlen(tissues[0])) ;
  if (strlen(tissues[0]) == 1) 
    { types[0] = "" ;types[1] = 0 ; }
  /* read all the lists and create an array of flags, one bit per target */
  for (nn = 0, tissue = tissues ; *tissue ; tissue++)
    for (type = types ; *type ; type++)
      {
	dict = dictHandleCreate (50000000, 0) ;
	nmax = 0 ;
	for (ii=1, target = targets ; *target ; ii<<=1, target++)
	  {
	    if (!strcmp (*target, "nx")) { flagNx = ii ; continue ; }
	    sprintf (buf, "HITS_%s/%s.%s%s.error.list0", ii<0?"u":"error",*tissue, *type, *target) ; 
	    fprintf (stderr, "// Opening %s  type=#%s#\n", buf, *type) ;
	    ai = aceInCreate (buf, 0, 0) ;
	    while (aceInCard (ai))
	      if ((cp = aceInWord (ai)))
		{
		  if (sx->ignoredProbeDict && dictFind (sx->ignoredProbeDict, cp, 0))
		    continue ;
		  dictAdd (dict, cp, &n) ;
		  if (n > nmax) nmax = n ;
		  ip = arrayp (flags, n + nn, unsigned int) ;
		  *ip |= ii ;
		}
	    ac_free (ai) ;
	   }
	nn += nmax ;
	dictDestroy (dict) ;
      }
  /* create a global Refseq column */ 
  if (flagNx)
    flagNx2 = (flagNx >> 1) | (flagNx >> 2) ;
  /* create a global exonic column */ 
  for (flag = 0, ii=4, target = targets+2 ; *target ; ii<<=1, target++)
    flag |= ii ;
  jj = ii ; *target = "Super transcriptome" ;
  for (nUnion = i = 0, ip = arrp (flags, i, unsigned int) ; i < arrayMax(flags) ; ip++, i++)
    {
      if (*ip) nUnion++ ;
      if (*ip & flag) *ip |= jj ;
      if (*ip & flagNx2) *ip |= flagNx ;
    }
  
  /* Count and report the categories */
  freeOutf("%s", tissues[0]) ;
  for (title = titles ; *title ; title++)
    freeOutf("\t%s", *title) ;
  
  for (ii=1, target = targets, title = titles ; *target ; ii<<=1, target++, title++)
    {
      freeOutf ("\n%s", *title) ;
      for (jj=1, target2 = targets ; *target2 ; jj<<=1, target2++)
	{
	  flag = ii | jj ;
	  for (n = i = 0, ip = arrp (flags, i, unsigned int) ; i < arrayMax(flags) ; ip++, i++)
	    if ((*ip & flag) == flag) n++ ;
	  freeOutf("\t%d", n) ;
	  
	}
    }
  freeOutf ("\n%d tags are mapped in the Union of all these targets\n", nUnion) ;
  return ;
}

/*************************************************************************************/
/*************************************************************************************/
#define GGMAX 400
#define DAMPER (double)30.0

typedef struct ggStruct { int gene, ln ; float s[GGMAX] ; } GG ;
static Array sxParseTagCounts (DICT *dict, AC_HANDLE h)
{
  Array aa = arrayCreate (100000, GG) ;
  ACEIN ai = aceInCreate (0, 0, h) ;
  double logDeux = log(2.0) ;
  char *cp, buf[256] ;
  float x,y,z = 0 ;
  double total = 0, totalb = 0 ;
  GG *gg = 0 ;
  int t, t2 ;
  int nn = 0, i, ln ;
  BOOL isTaq, isX ;

  aceInSpecial (ai, "\n\t\"\\") ;
  while (aceInCard (ai))
    {
      cp = aceInWord (ai) ; if (!cp) continue ;
      if (!strcasecmp (cp, "Element") || !strcasecmp (cp, "Intron"))
	{
	  cp = aceInWord (ai) ; if (!cp) goto done ;	  
	  dictAdd (dict, cp, &nn) ;
	  gg = arrayp (aa, nn, GG) ;
	  gg->gene = nn ;
	  if (0) for (i = 0 ; i < GGMAX ; i++) gg->s[i] = -9999 ;
	  continue ;
	}
      if (strstr(cp, "Length"))
	{
	  if (aceInFloat (ai, &z)) 
	    gg->ln = z ;
	  continue ;
	}
       if (! strstr(cp, "_A") & ! strstr(cp, "_B"))
	continue ;
      if (!strncmp(cp, "X_",2))
	isX = TRUE ;
      else
	isX = FALSE ;
      if (strstr(cp, "Anti"))
	continue ;
      if (strstr(cp, "TAQ_") || strstr(cp, "MAQC_"))
	isTaq = TRUE ;
      else
	isTaq = FALSE ;
      /*
	BOOL isDg ;
	if (strstr(cp, "HELdg"))
	isDg = TRUE ;
	else
	isDg = FALSE ;
      */
      strcpy (buf, cp) ;
      if (dictAdd (dict, cp, &t))
	messcrash ("Tag %s was not declared in sxExp ()", cp) ;
      if (t >= GGMAX) messcrash ("Too many tracks (%d > %d), edit the definition of the GG structure", t, GGMAX) ;
      sxTrackLocate (dict, cp, 't', 0,0,0,0, &t, 0, 0) ;
      if (aceInFloat (ai, &z)) 
	{
	  if (isTaq)
	    {
	      gg->s[t]= exp(z*logDeux) - DAMPER ;
	      ln = gg->ln ;
	      if (ln < 200) ln = 200 ;
	      gg->s[t] *= ln/200.0 ;
	    }
	  else
	    {
	      gg->s[t]= z ; /* number of tags */
	      if (! isX) total += z ;
	    }
	}
      else 
	continue ; 

      if (1)
	{
	  cp = strstr(buf, "_B") ; if (cp && *(cp+2)=='e') cp =  strstr(cp+1, "_B") ; 
	  if (cp) 
	    {
	      *(cp+1) = 'A' ;
	      if (sxTrackLocate (dict, buf, 't', 0,0,0,0, &t2, 0, 0))
		{
		  x = gg->s[t] ; y = gg->s[t2] ;
		  if (x > 0 || y > 0)
		    gg->s[t2+3] = 1000000 + log(DAMPER + x)/logDeux - log(DAMPER + y)/logDeux ; 
		}
	    }
	}
      if (isTaq)
	continue ;
      aceInNext (ai) ;
      if ((i= aceInFloat (ai, &z)) )
	{
	  gg->s[t+1]= z ; /* number of base aligned */
	  totalb += z ;
	}
      aceInNext (ai) ;
      if (aceInFloat (ai, &z)) 
	gg->s[t+2]=z ; else continue ; 
    }
 done:
  printf ("Counted %d Million tags, %d Mb\n", (int)(total/1000000), (int)(totalb/1000000)) ;
  return aa ;
} /* sxParseTagCounts */

/*************************************************************************************/

static BOOL sxTrackLocate (DICT *dict
			   , const char *nam1, char type1
			   , const char *nam2, char type2
			   , const char *nam3, char type3
			   , int *t1p, int *t2p, int *t3p)
{
  int t1, t2, t3 ;

  /* locate the lines */
  t1 = t2 = t3 = -1 ;

  if (nam1)
    {
      if (!dictFind (dict, nam1, &t1)) messcrash ("Bad name %s in sxExpExp", nam1) ;
      switch (type1)
	{
	case 't': t1 = 5 * t1 + 0 ; break ; /* tags rationalise: i.e divided by number of target genes */
	case 'b': t1 = 5 * t1 + 1 ; break ; /* bp rationaliased: i.e divided by number of target genes */
	case 'c': t1 = 5 * t1 + 2 ; break ; /* corrected in case of Helicos HH */
	case 'L': t1 = 5 * t1 + 3 ; break ; /* log ratio */
	default:
	  messcrash ("Bad type1 %c in sxExpExp", type1) ;
	}
      *t1p = t1 ;
    }
  if (nam2)
    {
      if (!dictFind (dict, nam2, &t2)) messcrash ("Bad name %s in sxExpExp", nam2) ;
      switch (type2)
	{
	case 't': t2 = 5 * t2 + 0 ; break ; /* tags rationalise: i.e divided by number of target genes */
	case 'b': t2 = 5 * t2 + 1 ; break ; /* bp rationaliased: i.e divided by number of target genes */
	case 'c': t2 = 5 * t2 + 2 ; break ; /* corrected */
	case 'L': t2 = 5 * t2 + 3 ; break ; /* log ratio */
	default:
	  messcrash ("Bad type2 %c in sxExpExp", type2) ;
	}
      *t2p = t2 ;
    }
  if (nam3)
    {
      if (!dictFind (dict, nam3, &t3)) messcrash ("Bad name %s in sxExpExp", nam3) ;
      switch (type3)
	{
	case 't': t3 = 5 * t3 + 0 ; break ; /* tags */
	case 'b': t3 = 5 * t3 + 1 ; break ; /* bp */
	case 'c': t3 = 5 * t3 + 2 ; break ; /* corrected */
	case 'L': t3 = 5 * t3 + 3 ; break ; /* log ratio */
	default:
	  messcrash ("Bad type3 %c in sxExpExp", type3) ;
	}
      *t3p = t3 ;
    }
  return TRUE ;
} /* sxExpExp */

/*************************************************************/
/* find the value that leaves out trim per cent of the histogram */
static float sxExpTrimMax (Array aa, int t1, int trim)
{
  int ii = arrayMax (aa) ;
  Array zz = arrayCreate (ii, float) ;
  float mx = 0, *zp = arrayp (zz, ii-1, float) + 1 ;
  GG *gg = arrp (aa, 0, GG) - 1 ;
  
  while (zp--, gg++, ii--)
    if (gg->gene) *zp = gg->s[t1] ;
  arraySort (zz, floatOrder) ;
  ii = arrayMax (aa) - (int)(trim * arrayMax(aa)/100.0) ;
  if (ii == arrayMax(aa)) ii-- ;
  if (ii<0) ii = 0 ;
  mx = array (zz, ii, float) ;
  
  arrayDestroy (zz) ;
  
  return mx ;
} /* sxExpTrimMax */

/*************************************************************************************/
typedef struct e2Struct { float x, y ; } E2 ;
static Array sxExp2d (Array aa, int t1, int t2, float mx1, float mx2
		      , char type1, char type2, char *nam1, char *nam2
		      , AC_HANDLE h)
{
  int ii, jj ;
  Array hh = arrayHandleCreate (arrayMax(aa), E2, h) ;
  E2 *up ;
  float x, y ;
  BOOL isTaq1 = FALSE, isTaq2 = FALSE ;
  double logDeux = log(2.0) ;
  GG *gg ;

  if (strstr (nam1, "TAQ") || strstr (nam1, "MAQC"))
    isTaq1 = TRUE ;
  if (strstr (nam2, "TAQ") || strstr (nam2, "MAQC"))
    isTaq2 = TRUE ;
  for (ii = jj = 0, gg = arrp (aa, ii, GG) ; ii < arrayMax (aa) ; gg++, ii++)
    {
      if (! gg->gene) 
	continue ;
      x = gg->s[t1] ; y = gg->s[t2] ;
      if (x < DAMPER && y < DAMPER)
	continue ;
      if (type1 == 'L' && 
	  (
	   (x == 0 || y == 0) || (x < 10 &&  y < 10)
	   ))
	continue ;
      if ((isTaq1 && x == 0) || (isTaq2  && y == 0))
	continue ;
      if ((!mx1 || x <= mx1) && (!mx2 || y <= mx2))
	{
	  up = arrayp (hh, jj++, E2) ;
	  if (type1 != 'L') 
	    x = log(DAMPER + x)/logDeux ;
	  else
	    x -= 1000000 ;
	  if (type2 != 'L') 
	    y = log(DAMPER + y)/logDeux ;
	  else
	    y -= 1000000 ;
	  up->x = x ;
	  up->y = y ;
	  if (0) fprintf(stderr, "j=%d t1:t2 = %d:%d x:y=%f:%f\n",jj,t1,t2,x,y);
	}
    }
  return hh ;
} /* sxExp2d */

/*************************************************************************************/

static float sxExpCorrel (SX *sx, Array hh)
{
  double r, x, y, xx, x2, yy, y2, xy, X, Y, XY ;
  int n = 0, ii ;
  E2 *up ;

  xx = x2 = yy = y2 = xy = 0 ;
  for (ii = 0, up = arrp (hh, 0, E2) ; ii < arrayMax(hh) ; up++, ii++)
    {
      x = up->x ; y = up->y ;
      n++ ; xx += x ; yy += y ; x2 += x*x ;y2 += y*y ; xy += x*y ;
    }
  if (!n) 
    n = 1 ;
  xx /= n ; yy /= n ;
  X = x2-n*xx*xx ; Y = y2 - n*yy*yy ; 
  XY = xy - n*xx*yy ;
  
  r=0 ; if(X*Y > 0) r = XY / sqrt (X*Y) ;
  if (sx->count) r = n ; /* count */
  return r ;
} /* sxExpCorrel */

/*************************************************************************************/

static int sxExpPlot (ACEOUT ao, SX *sx, Array hh)
{
  double x, y ;
  int n = 0, ii ;
  E2 *up ;

  for (ii = 0, up = arrp (hh, 0, E2) ; ii < arrayMax(hh) ; up++, ii++)
    {
      x = up->x ; y = up->y ;
      n++ ; 
      aceOutf (ao, "%g\t%g\n", x, y) ;
    }
  return n ;
} /* sxExpPlot */

/*************************************************************************************/

static int sxExpVennPlot (SX *sx, Array aa, DICT *dict, 
			  const char *nam1, const char *nam2, const char *nam3, int trim)
{
  double x, y, z, limit = 1 ;
  int sign, ii, t1=0, t2=0, t3=0 ;
  int nABC, nAB, nBC, nCA, nA, nB, nC ;
  GG *gg ;
  char *nam0 = "-" ;

  
  if (!sxTrackLocate (dict, nam1, 'L', nam2, 'L', nam3, 'L', &t1, &t2, &t3)) 
    goto done ;
  for (sign = -1 ; sign < 2 ; sign += 2)
    {
      nABC = nAB = nBC = nCA = nA = nB = nC = 0 ; 
      for (ii = 0, gg = arrp (aa, 0, GG) ; ii < arrayMax(aa) ; gg++, ii++)
      {
	if (! gg->gene) 
	  continue ;
	x = gg->s[t1] ; y = gg->s[t2] ; z = gg->s[t3] ;
	x = sign * gg->s[t1] ; y = sign * gg->s[t2] ; z = sign * gg->s[t3] ;
	if(x > limit && y > limit && z > limit) nABC++ ;
	else if(x <= limit && y > limit && z > limit) nBC++ ;
	else if(x > limit && y <= limit && z > limit) nCA++ ;
	else if(x > limit && y > limit && z <= limit) nAB++ ;
	else if(x > limit && y <= limit && z <= limit) nA++ ;
	else if(x <= limit && y > limit && z <= limit) nB++ ;
	else if(x <= limit && y <= limit && z > limit) nC++ ;
      }
      printf ("Venn diagram, Over expressed in %s sample at least 2-fold\n"
	      , sign == -1 ? "A" : "B") ;
      printf ("%d\t%s\t%s\t%s\n", nABC, nam1,nam2,nam3) ;
      printf ("%d\t%s\t%s\t%s\n", nBC, nam0,nam2,nam3) ;
      printf ("%d\t%s\t%s\t%s\n", nCA, nam1,nam0,nam3) ;
      printf ("%d\t%s\t%s\t%s\n", nAB, nam1,nam2,nam0) ;
      printf ("%d\t%s\t%s\t%s\n", nA, nam1,nam0,nam0) ;
      printf ("%d\t%s\t%s\t%s\n", nB, nam0,nam2,nam0) ;
      printf ("%d\t%s\t%s\t%s\n", nC, nam0,nam0,nam3) ;
    }
 done:
  return 1 ;
} /* sxExpVennPlot */

/*************************************************************************************/
/* parse the file sx->support which should look like
Element "1_random__36978_35538"
Intron
OtherSupport    "ensemblGeneFeb09" "ENSG00000215719"
OtherSupport    "geneidPredictionsSpainMay06" "chr1_random_2"
*/

static int sxExpIntronSupport (SX *sx, Array aa, KEYSET i2s, DICT *dict, DICT *eDict, DICT *sDict, KEYSET sTotal, KEYSET estCount)
{
  int iss ;
  int nE = 0, ne, line ;
  const char *ccp ;
  char *cp ;
  GG *gg ;
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (sx->support, 0, h) ;
  aceInSpecial (ai, "\t\n\\\"") ;
  
  dictAdd (sDict, "SwissInstBioinfoGeneSept07", &ne) ;
  dictAdd (sDict, "AV", &ne) ;
  dictAdd (sDict, "EST", &ne) ;
  dictAdd (sDict, "SS", &ne) ;
  dictAdd (sDict, "Hinv", &ne) ;
  dictAdd (sDict, "UCSCknownGeneSep08", &ne) ;
  dictAdd (sDict, "JJ", &ne) ;
  dictAdd (sDict, "ensemblGeneFeb09", &ne) ;
  dictAdd (sDict, "RefGene", &ne) ;
  dictAdd (sDict, "NM", &ne) ;
  dictAdd (sDict, "RefSeq", &ne) ;
  dictAdd (sDict, "NscanPasaGeneOct07", &ne) ;
  dictAdd (sDict, "NscanGeneApr06", &ne) ;
  dictAdd (sDict, "augustusHintsJuly07", &ne) ;
  dictAdd (sDict, "augustusXRAJuly07", &ne) ;
  dictAdd (sDict, "augustusAbinitioJuly2007", &ne) ;
  dictAdd (sDict, "vegaGeneJune08", &ne) ;
  dictAdd (sDict, "HINV", &ne) ;
  dictAdd (sDict, "vegaPseudoGeneJune2008", &ne) ;
  dictAdd (sDict, "genscanBurgeApr06", &ne) ;
  dictAdd (sDict, "sgpGeneOct07", &ne) ;
  dictAdd (sDict, "geneidPredictionsSpainMay06", &ne) ;
  dictAdd (sDict, "XM", &ne) ;
  dictAdd (sDict, "exoniphy", &ne) ;


  line = 0 ;
  while (aceInCard (ai))
    {
      line++ ;
      ccp = aceInWord (ai) ; 
      if (ccp && ! strcmp (ccp, "Element"))
	{
	  ne = 0 ;
	  ccp = aceInWord (ai) ; 
	  if (ccp)
	    {
	      cp = ac_unprotect (ccp, h) ;
	      if (cp && *cp)
		{
		  KEY x ;

		  dictAdd (dict, cp, &ne) ;
		  gg = arrayp (aa, ne, GG) ;
		  gg->gene = ne ;
		  dictAdd (eDict, cp, &ne) ;
		  nE++ ; 
		  x = keySet (i2s, ne) ; /* make room */
		  if (x) {} ; /* for compile happiness */
		}
	    }
	  continue ;
	}
      if (ccp && 
	  (
	   !strcmp (ccp, "AV") ||
	   !strcmp (ccp, "JJ") ||
	   !strcmp (ccp, "NM") ||
	   !strcmp (ccp, "XM") ||
	   !strcmp (ccp, "SS") ||
	   !strcmp (ccp, "RefSeq")
	   )
	  )

	{
	  if (ne && ccp)
	    {
	      dictAdd (sDict, ccp, &iss) ;
	      keySet (i2s, ne) |= (0x1 << iss) ;
	      keySet (sTotal, iss)++ ;
	    }
	  continue ;
	}
      if (ccp && !strcmp (ccp, "EST"))
	{ /* add this count as if it was found in the SEQC file */
	  int x = 0 ;
	  if (aceInInt(ai, &x))
	    {
	      keySet (estCount, ne) = x ;
	      dictAdd (sDict, "EST", &iss) ;
		  if (! (keySet (i2s, ne) & (0x1 << iss)))
		    {
		      keySet (sTotal, iss)++ ;
		      keySet (i2s, ne) |= (0x1 << iss) ;
		    }
	    }
	  continue ;
	}
      if (ccp && ! strcmp (ccp, "OtherSupport"))
	{
	  ccp = aceInWord (ai) ; 
	  if (ne && ccp)
	    {
	      cp = ac_unprotect (ccp, h) ;
	      if (cp && *cp)
		{
		  dictAdd (sDict, cp, &iss) ;
		  if (! (keySet (i2s, ne) & (0x1 << iss)))
		    {
		      keySet (sTotal, iss)++ ;
		      keySet (i2s, ne) |= (0x1 << iss) ;
		    }
		}
	    }
	  continue ;
	}
    }

  ac_free (h) ;
  return nE ;
} /* sxExpIntronSupport */

/*************************************************************************************/

static int sxExpIntron (SX *sx, Array aa, DICT *dict)
{
  GG *gg ;
  int nI[100], nJ[100], nn, nEst, total ;
  int iss, nS, issNM = 0 ;
  int i, ii, jj, igg, t, tt[128], nIS[1024], nJS[1024] ;
  char **s, *sample[] = {"A","B",0}, buf[1000] ;
  char **t1, *tag[] = {"X", "ILMnoSt", "ILM", "ILM_tot", "ILMnoSt_tot","LIF", "LIF_tot", "HEL", "HH", 0} ;
  double x ;
  KEY inNM = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  KEYSET i2n = keySetHandleCreate (h) ; /* histo of support of all introns */
  KEYSET j2n = keySetHandleCreate (h) ; /* histo of support of non NM introns */
  KEYSET i2s = keySetHandleCreate (h) ; /* actual support of all each intron */
  KEYSET sTotal = keySetHandleCreate (h) ; /* actual support of all each intron */
  KEYSET estCount = keySetHandleCreate (h) ; /* EST suppprt of this intron */
  DICT *eDict = dictHandleCreate (300000, h) ; /* to hold the intron names to support association */
  DICT *sDict =dictHandleCreate (30, h) ; /* all the methods name */

  if (sx->support)
    sxExpIntronSupport (sx, aa, i2s, dict, eDict, sDict, sTotal, estCount) ;
  dictFind (sDict, "NM", &issNM) ;
  /* locate all tracks */
  memset (tt, 0, sizeof(tt)) ;
  memset (nIS, 0, sizeof(nIS)) ;
  memset (nJS, 0, sizeof(nJS)) ;
  memset (nI, 0, sizeof(nI)) ;
  memset (nJ, 0, sizeof(nJ)) ;
  for (ii = 0, t1 = tag ; *t1 ; ii++, t1++)
    for (jj = 0, s = sample ; *s ; jj++, s++)
      {
	sprintf (buf, "%s_%s", *t1, *s) ;
	tt[3*ii+jj] = -1 ;
	sxTrackLocate (dict, buf, 't', 0, 0, 0, 0, &tt[3*ii+jj], 0, 0) ;
      }
  /* scan all the introns, get the counts, attribute them to the various prediction and their genes */
  for (igg = 0, gg = arrp (aa, 0, GG) ; igg < arrayMax(aa) ; gg++, igg++)
    {
      if (! gg->gene) 
	continue ;
      nn = 0 ; nS = 0 ; nEst = 0 ;
      for (ii = 0, t1 = tag ; *t1 ; ii++, t1++)
	for (jj = 0, s = sample ; *s ; jj++, s++)
	  {
	    t = tt[3*ii+jj] ;
	    if (t>=0)
	      {
		x = gg->s[t] ;
		nn += x ;
	      }
	  }
      if (nn > 0 && dictFind (eDict, dictName(dict,gg->gene), &iss))
	{
	  nS = keySet (i2s, iss) ;
	  nEst = keySet (estCount, iss) ;
          nn += nEst ;
	}
      inNM = issNM ? ((nS>>issNM) & 0x1) : 0 ;
      if (nn)
	{
	  int j = 1 ;
	  for (i=0;i<10;i++, j++)
	    if (nn >= j) 
	      {
		nI[i]++ ;
		if (!inNM) nJ[i]++ ;

		for (iss = 0 ; nS && iss <= dictMax(sDict) ; iss++)
		  if ((nS>>iss) & 0x1)
		    {
		      nIS[i+10*iss]++ ;
		      if (!inNM)
			nJS[i+10*iss]++ ;
		    }
	      }
	}
      jj = 0 ;
      while (nn>0)
	{
	  if (nn < 10)
	    {
	      keySet (i2n, 10*jj+nn) += 1 ;
	      if (!inNM)
		keySet (j2n, 10*jj+nn) += 1 ;
	    }
	  nn/=10 ; jj++ ;
	}
    }
  printf ("Method\tAll introns");
  for (i=1 ; i<= 10 ; i++)
    printf ("\t>=%d tag in NM\tNot in NM\tNot validated\t%% Validated NM in this annotation (sensitivity)\t%% Novel relative to NM (discovery)\t%% Annotated supported (Specificity)", i) ;
  printf ("\n") ;

  total = 0 ;
  /*
    for (iss = 1 ; iss <= dictMax(sDict) && iss<31 ; iss++)
    total += keySet (sTotal, iss) ;
  */
 
  total = dictMax (eDict) ;
  printf ("Any\t%d",total) ;

  for (i = 0 ; i < 10 ; i++)
    {
      printf ("\t%d\t%d\t%d", nI[i]-nJ[i], nJ[i], total - nI[i]) ;
      printf ("\t%.2f", 100.0*(nI[i]-nJ[i])/nIS[i+10*issNM]) ;
      printf ("\t%.2f", 100.0*(nJ[i])/nI[i]) ;
      printf ("\t%.2f", 100.0*(nI[i])/total) ;
    }
  printf ("\n") ;
  if (dictMax(sDict) > 30)
    printf("More than 30 methods, only 30 can be reported, sorry\n") ;
  for (iss = 1 ; iss <= dictMax(sDict) && iss<31 ; iss++)
    {
      printf ("%s\t%d", dictName (sDict, iss), keySet (sTotal, iss)) ;
      for (i = 0 ; i < 10 ; i++)
	{
	  printf ("\t%d\t%d\t%d", nIS[i+10*iss] -  nJS[i+10*iss], nJS[i+10*iss], keySet (sTotal, iss)- nIS[i+10*iss] ) ;
	  printf ("\t%.2f", 100.0*(nIS[i+10*iss]-nJS[i+10*iss])/nIS[i+10*issNM]) ;
	  printf ("\t%.2f", 100.0*(nJS[i+10*iss])/nIS[i+10*iss]) ;
	  printf ("\t%.2f", 100.0*(nIS[i+10*iss])/keySet (sTotal, iss)) ;
	}
      printf ("\n") ;
    }
	
  printf ("\nIntron supported by at least so many tags\tAny\tNot in NM\n") ;
  for (ii = 0 ; ii < keySetMax (i2n) ; ii++)
    {
      int k = 1 ;
      jj = ii ;
      while (jj >= 10) { k = 10*k ; jj -= 10 ;}
      if (jj) k *= jj ;
      printf ("%d\t%d\t%d\n", k, keySet(i2n,ii), keySet(j2n,ii)) ;
    }
  
  ac_free (h) ;
  return 1 ;
} /* sxExpIntron */

/*************************************************************************************/

static float sxExpExp (SX *sx, Array aa, DICT *dict, char *nam1, char type1, char *nam2, char type2, int trim)
{
  AC_HANDLE h = ac_new_handle () ;
  int t1, t2 ;
  float r = 0, mx1 = 0, mx2 = 0 ;
  Array hh ;
  trim=0;
  /* locate the lines */
  if (!sxTrackLocate (dict, nam1, type1, nam2, type2, 0, 0, &t1, &t2, 0)) 
    goto done ;

  if (trim > 0 && trim< 10)
    {
      mx1 = sxExpTrimMax (aa, t1, trim) ;
      mx2 = sxExpTrimMax (aa, t2, trim) ;
    }
  else if (trim)
    {
      mx1 = trim ;
      mx2 = trim ;
    }
  /* create the 2d histo */
  hh = sxExp2d (aa, t1, t2, mx1, mx2, type1, type2, nam1, nam2, h) ;
  /* export, compute the correl etc */
  if (sx->plot)
    {
      ACEOUT ao = aceOutCreateToFile (hprintf(h,"%s_%s.%c%c.txt",nam1, nam2, type1, type2),"w",h) ;
      sxExpPlot (ao, sx, hh) ;
      fprintf(stderr, "%s_%s.%c%ctxt",nam1, nam2, type1, type2) ;
    }
  else
    r = sxExpCorrel (sx, hh) ;

 done:
  ac_free (h) ;
  return r ;
} /* sxExpExp */

/*************************************************************************************/

static int sxExpOne (SX *sx, Array aa, DICT *dict, int limit)
{
  double x, y ;
  int i, j, ii ;
  int nA, nB, nAB, nU ;
  int npA[20], npB[20], npAB[20], npU[20] ;
  int npAX[20], npBX[20], npABX[20], npUX[20] ;
  int nnA, nnB, nnAB, nnU, any ;
  int nnA2, nnB2, nnAB2, nnU2 ;
  int nnA3, nnB3, nnAB3, nnU3 ;
  int nnA4, nnB4, nnAB4, nnU4 ;
  int nnA5, nnB5, nnAB5, nnU5 ;
  GG *gg ;

  char * tag[] = {"X", 0, "R454", 0, "HEL", 0,  "ILMnoSt", "ILM", "ILM_tot", "ILMnoSt_tot", 0, "LIF", "LIF_tot", 0, (char*)1, 0} ;
  char * tag2[] = {"X", "454", "HEL", "ILM", "LIF", 0, 0} ;
  char buf1[100], buf2[100] ;
  int tA[20], tB[20] ;
  KEYSET anyHisto = keySetCreate () ;
  
  memset (tA, 0, sizeof(tA)) ;
  memset (tB, 0, sizeof(tB)) ;
  
  /* locate everybody */
  for (i = 0 ; tag[i] != (char*)1 ; i++)
    {
      if (tag[i] == 0)
	continue ;
      sprintf (buf1, "%s_%s", tag[i], "A") ;
      sprintf (buf2, "%s_%s", tag[i], "B") ;
      
      sxTrackLocate (dict, buf1, 't', buf2, 't', 0, 'L', tA+i, tB+i, 0) ;
    }
  for (i = 0 ; i < 20 ; i++)
    {
      npA[i] = npB[i] = npAB[i] = npU[i] = 0 ;
      npAX[i] = npBX[i] = npABX[i] = npUX[i] = 0 ;
    }
  nnA = nnB = nnAB = nnU = 0 ;
  nnA2 = nnB2 = nnAB2 = nnU2 = 0 ;
  nnA3 = nnB3 = nnAB3 = nnU3 = 0 ;
  nnA4 = nnB4 = nnAB4 = nnU4 = 0 ;
  nnA5 = nnB5 = nnAB5 = nnU5 = 0 ;
  for (ii = 0, gg = arrp (aa, 0, GG) ; ii < arrayMax(aa) ; gg++, ii++)
    {
      if (! gg->gene) 
	continue ;
      nA = nB = nAB = nU = 0 ;
      any = 0 ;
      for (i = 0 ; i < 20 ; i++)
	{
	  x = 0 ; y = 0 ;
	  for (j = 0 ; tA[i+j] ; j++)
	    {
	      x += gg->s[tA[i+j]] ;
	      y += gg->s[tB[i+j]] ;
	    }
	  any += x + y ;
	  if(any < 0) messcrash ("any < 0, x=%g y=%g g=%s", x, y, dictName (dict, gg->gene)) ;
	  if (x > limit) { npA[i]++ ; nA++ ; }
	  if (y > limit) { npB[i]++ ; nB++ ; }
	  if (x > limit && y > limit)  { npAB[i]++ ; nAB++ ; }
	  if (x > limit ||  y > limit)  { npU[i]++ ; nU++ ; }
	  if (j>1) i += j - 1 ;
	}
     for (i = 0 ; i < 20 ; i++)
	{
	  x = 0 ; y = 0 ;
	  for (j = 0 ; tA[i+j] ; j++)
	    {
	      x += gg->s[tA[i+j]] ;
	      y += gg->s[tB[i+j]] ;
	    }
	  any += x + y ;
	  if(any < 0) messcrash ("any < 0") ;
	  if (nA == 1 && x > limit) { npAX[i]++ ; }
	  if (nB == 1 && y > limit) { npBX[i]++ ; }
	  if (nAB == 1 && x > limit && y > limit)  { npABX[i]++ ; }
	  if (nU == 1 && (x > limit ||  y > limit))  { npUX[i]++ ; } 
	  if (j>1) i += j - 1 ;
	}

      if (!limit)
	keySet (anyHisto, any) += 1 ;
      if (nA) nnA++ ;
      if (nB) nnB++ ;
      if (nAB) nnAB++ ;
      if (nU) nnU++ ;

      if (nA > 1) nnA2++ ;
      if (nB > 1) nnB2++ ;
      if (nAB > 1) nnAB2++ ;
      if (nU > 1) nnU2++ ;

      if (nA > 2) nnA3++ ;
      if (nB > 2) nnB3++ ;
      if (nAB > 2) nnAB3++ ;
      if (nU > 2) nnU3++ ;

      if (nA > 3) nnA4++ ;
      if (nB > 3) nnB4++ ;
      if (nAB > 3) nnAB4++ ;
      if (nU > 3) nnU4++ ;

      if (nA > 4) nnA5++ ;
      if (nB > 4) nnB5++ ;
      if (nAB > 4) nnAB5++ ;
      if (nU > 4) nnU5++ ;
    }

  if (limit > 0)
    {
      printf ("Limit %d\tA or B\tA\tB\tA and B\tA only\tB only\n", limit) ;
      for (i = j = 0 ; i < 20 ; i++, j++)
	{
	  if (tA[i])
	    {
	      printf ("%s\t%d\t%d\t%d\t%d\t%d\t%d\n", tag2[j], npU[i], npA[i], npB[i], npAB[i],  npA[i] - npAB[i], npB[i] - npAB[i]) ;
	      while (tA[i]) i++ ;
	    }
	}
      printf ("%s\t%d\t%d\t%d\t%d\t%d\t%d\n", "Any", nnU, nnA, nnB, nnAB, nnA - nnAB, nnB - nnAB) ;
      printf ("%s\t%d\t%d\t%d\t%d\t%d\t%d\n", "At least 2", nnU2, nnA2, nnB2, nnAB2, nnA2 - nnAB2, nnB2 - nnAB2) ;
      printf ("%s\t%d\t%d\t%d\t%d\t%d\t%d\n", "At least 3", nnU3, nnA3, nnB3, nnAB3, nnA3 - nnAB3, nnB3 - nnAB3) ;
      printf ("%s\t%d\t%d\t%d\t%d\t%d\t%d\n", "At least 4", nnU4, nnA4, nnB4, nnAB4, nnA4 - nnAB4, nnB4 - nnAB4) ;
      printf ("%s\t%d\t%d\t%d\t%d\t%d\t%d\n", "At least 5", nnU5, nnA5, nnB5, nnAB5, nnA5 - nnAB5, nnB5 - nnAB5) ;
      for (i = j = 0 ; i < 20 ; i++, j++)
	{
	  if (tA[i])
	    {
	      printf ("Just %s\t%d\t%d\t%d\t%d\t%d\t%d\n", tag2[j], npUX[i], npAX[i], npBX[i], npABX[i], npAX[i] - npABX[i], npBX[i] - npABX[i]) ;
	      while (tA[i]) i++ ;
	    }
	}
      printf ("\n") ;      
    }
  else if(0) /* pas beau, on a trouve un meilleur pour les introns */
    {
      int i, n ;
      printf ("n\tGenes with n tags\t[10n,10n+9]\t[100n,100n+99]\t[1000n,1000n+9999]\t[10000n,10000n+999]\n") ;
      for (ii = 0 ; ii < keySetMax (anyHisto) && ii < 100 ; ii++)
	{
	  printf ("%d\t%d", ii, keySet (anyHisto, ii)) ;
	  for (n = 0, i = 10*ii ;  i < keySetMax (anyHisto) && i < 10 * (ii+1) ; i++)
	    n += keySet (anyHisto, i) ; 
	  printf ("\t%d", n) ;
	  for (n = 0, i = 100*ii ;  i < keySetMax (anyHisto) && i < 100 * (ii+1) ; i++)
	    n += keySet (anyHisto, i) ; 
	  printf ("\t%d", n) ;
	  for (n = 0, i = 1000*ii ;  i < keySetMax (anyHisto) && i < 1000 * (ii+1) ; i++)
	    n += keySet (anyHisto, i) ; 
	  printf ("\t%d", n) ;
	  for (n = 0, i = 10000*ii ;  i < keySetMax (anyHisto) && i < 10000 * (ii+1) ; i++)
	    n += keySet (anyHisto, i) ; 
	  printf ("\t%d\n", n) ;
	}
      printf("Remainder") ;
      for (n = 0, i = ii ;  i < keySetMax (anyHisto) ; i++)
	n += keySet (anyHisto, i) ; 
      printf ("\t%d", n) ;
      for (n = 0, i = 10*ii ;  i < keySetMax (anyHisto) ; i++)
	n += keySet (anyHisto, i) ; 
      printf ("\t%d", n) ;
      for (n = 0, i = 100*ii ;  i < keySetMax (anyHisto) ; i++)
	n += keySet (anyHisto, i) ; 
      printf ("\t%d", n) ;
      for (n = 0, i = 1000*ii ;  i < keySetMax (anyHisto) ; i++)
	n += keySet (anyHisto, i) ; 
      printf ("\t%d", n) ;
      for (n = 0, i = 10000*ii ;  i < keySetMax (anyHisto) ; i++)
	n += keySet (anyHisto, i) ; 
      printf ("\t%d\n", n) ;
    }

  keySetDestroy (anyHisto) ;
  return 0 ;
} /* sxExpOne */

/*************************************************************************************/

static int sxHammer (SX *sx, Array aa, DICT *dict, int limit)
{
  double x, y ;
  int i, j, ii ;
  int npA[20], npB[20], npAB[20], npU[20] ;
  int any ;
  GG *gg ;

  char * tag[] = {"X", 0, 0 } ;
  char buf1[100], buf2[100] ;
  int tA[20], tB[20] ;
  
  memset (tA, 0, sizeof(tA)) ;
  memset (tB, 0, sizeof(tB)) ;
  
  /* locate everybody */
  for (i = 0 ; tag[i] != 0 ; i++)
    {
      if (tag[i] == 0)
	continue ;
      sprintf (buf1, "%s_%s", tag[i], "A") ;
      sprintf (buf2, "%s_%s", tag[i], "B") ;
      
      sxTrackLocate (dict, buf1, 'b', buf2, 'b', 0, 'L', tA+i, tB+i, 0) ;
    }
  for (i = 0 ; i < 20 ; i++)
    {
      npA[i] = npB[i] = npAB[i] = npU[i] = 0 ;
    }

  for (ii = 0, gg = arrp (aa, 0, GG) ; ii < arrayMax(aa) ; gg++, ii++)
    {
      if (! gg->gene || ! gg->ln) 
	continue ;

      any = 0 ;
      for (i = 0 ; i < 20 ; i++)
	{
	  x = 0 ; y = 0 ;
	  for (j = 0 ; tA[i+j] ; j++)
	    {
	      x += gg->s[tA[i+j]]/gg->ln ;
	      y += gg->s[tB[i+j]]/gg->ln ;
	    }
	  any += x + y ;
	  if(any < 0) messcrash ("any < 0, x=%g y=%g g=%s", x, y, dictName (dict, gg->gene)) ;
	  for (limit = 1 ; limit <= 5 ; limit++)
	    {
	      if (x > limit) { npA[limit]++ ; }
	      if (y > limit) { npB[limit]++ ; }
	      if (x > limit && y > limit)  { npAB[limit]++ ; }
	      if (x > limit ||  y > limit)  {  npU[limit]++ ; }
	    }
	  if (j>1) i += j - 1 ;
	}
    }

  printf ("Limit %d\tA or B\tA\tB\tA and B\tA only\tB only\n", limit) ;
  for (i = 1 ; i <= 5 ; i++)
    {
      printf ("%d\t%d\t%d\t%d\t%d\t%d\t%d\n", limit, npU[i], npA[i], npB[i], npAB[i],  npA[i] - npAB[i], npB[i] - npAB[i]) ;
    }

  return 0 ;
} /* sxHammer */

/*************************************************************************************/

static void sxExp (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  DICT *dict = dictHandleCreate (100000, h) ;
  Array aa = 0 ;
  char **s, *sample[] = {"A","B",0} ;
  char u1, u2, **t1, **t2 ;
  char * tag[] = {"X", "X_st", "ILMnoSt", "ILM", "LIF", "HEL", "HELdg", "R454", "MAQC", "TAQ", "ILMnoSt_tot", "ILM_tot", "LIF_tot", 0, 0} ;
  char buf1[100], buf2[100] ;
  int trim = 1 ;

  dictAdd (dict, "TAQ_A", 0) ;
  dictAdd (dict, "TAQ_B", 0) ;
  dictAdd (dict, "MAQC_A", 0) ;
  dictAdd (dict, "MAQC_B", 0) ;
  dictAdd (dict, "R454_A", 0) ;
  dictAdd (dict, "R454_B", 0) ;
  dictAdd (dict, "ILM_A", 0) ;
  dictAdd (dict, "ILM_B", 0) ;
  dictAdd (dict, "ILMnoSt_A", 0) ;
  dictAdd (dict, "ILMnoSt_B", 0) ;
  dictAdd (dict, "ILM_tot_A", 0) ;
  dictAdd (dict, "ILM_tot_B", 0) ;
  dictAdd (dict, "ILMnoSt_tot_A", 0) ;
  dictAdd (dict, "ILMnoSt_tot_B", 0) ;
  dictAdd (dict, "HEL_A", 0) ;
  dictAdd (dict, "HEL_B", 0) ;
  dictAdd (dict, "HELdg_A", 0) ;
  dictAdd (dict, "HELdg_B", 0) ;
  dictAdd (dict, "HH_A", 0) ;
  dictAdd (dict, "HH_B", 0) ;
  dictAdd (dict, "LIF_A", 0) ;
  dictAdd (dict, "LIF_B", 0) ;
  dictAdd (dict, "LIF_tot_A", 0) ;
  dictAdd (dict, "LIF_tot_B", 0) ;
  dictAdd (dict, "S_cap_A", 0) ;
  dictAdd (dict, "S_cap_B", 0) ;
  dictAdd (dict, "X_A", 0) ;
  dictAdd (dict, "X_B", 0) ;
  dictAdd (dict, "X_st_A", 0) ;
  dictAdd (dict, "X_st_B", 0) ;
  dictAdd (dict, "X_pA_A", 0) ;
  dictAdd (dict, "X_pA_B", 0) ;
  dictAdd (dict, "X_tot_A", 0) ;
  dictAdd (dict, "X_tot_B", 0) ;
  dictAdd (dict, "X_tot_ns_A", 0) ;
  dictAdd (dict, "X_tot_ns_B", 0) ;
  dictAdd (dict, "X_cap_A", 0) ;
  dictAdd (dict, "X_cap_B", 0) ;
  aa = sxParseTagCounts (dict, h) ; /* must come after adding TAQ in the dictionary */

  if (sx->intron)
    {
      sxExpIntron (sx, aa, dict) ;
      goto done ;
    }
  if (sx->plot)
    {
      if (0)
	{
	  sxExpExp (sx, aa, dict, "ILM_A", 'L', "MAQC_A", 'L', trim) ;
	  sxExpExp (sx, aa, dict, "LIF_A", 'L', "MAQC_A", 'L', trim) ;
	  sxExpExp (sx, aa, dict, "LIF_A", 'L', "R454_A", 'L', trim) ;
	  
	  sxExpExp (sx, aa, dict, "ILM_A", 't', "LIF_A", 't', trim) ;
	  sxExpExp (sx, aa, dict, "ILM_B", 't', "LIF_B", 't', trim) ;
	  sxExpExp (sx, aa, dict, "ILM_A", 'L', "LIF_A", 'L', trim) ;
	  
	  sxExpExp (sx, aa, dict, "ILM_A", 'b', "LIF_A", 'b', trim) ;
	  sxExpExp (sx, aa, dict, "ILM_B", 'b', "LIF_B", 'b', trim) ;
	  sxExpExp (sx, aa, dict, "ILM_A", 'b', "R454_A", 'b', trim) ;
	}
      else
	{
	  sxExpExp (sx, aa, dict, "TAQ_A", 'L', "X_A", 'L', trim) ;
	  sxExpExp (sx, aa, dict, "TAQ_A", 'L', "ILMnoSt_A", 'L', trim) ;
	  sxExpExp (sx, aa, dict, "TAQ_A", 'L', "ILM_A", 'L', trim) ;
	  sxExpExp (sx, aa, dict, "TAQ_A", 'L', "LIF_A", 'L', trim) ;
	  sxExpExp (sx, aa, dict, "TAQ_A", 'L', "HEL_A", 'L', trim) ;
	  sxExpExp (sx, aa, dict, "TAQ_A", 'L', "HELdg_A", 'L', trim) ;
	  sxExpExp (sx, aa, dict, "ILMnoSt_A", 'L', "LIF_A", 'L', trim) ;
	  sxExpExp (sx, aa, dict, "ILM_A", 'L', "LIF_A", 'L', trim) ;
	  sxExpExp (sx, aa, dict, "ILMnoSt_A", 'L', "ILM_A", 'L', trim) ;
	  sxExpExp (sx, aa, dict, "ILM_A", 'L', "HEL_A", 'L', trim) ;
	  sxExpExp (sx, aa, dict, "ILMnoSt_A", 'L', "HELdg_A", 'L', trim) ;
	  sxExpExp (sx, aa, dict, "LIF_A", 'L', "HEL_A", 'L', trim) ;
	  sxExpExp (sx, aa, dict, "LIF_A", 'L', "HELdg_A", 'L', trim) ;
	  sxExpExp (sx, aa, dict, "X_A", 't', "X_B", 't', trim) ;
	  sxExpExp (sx, aa, dict, "X_st_A", 't', "X_st_B", 't', trim) ;
	  sxExpExp (sx, aa, dict, "ILM_A", 't', "ILM_B", 't', trim) ;
	  sxExpExp (sx, aa, dict, "LIF_A", 't', "LIF_B", 't', trim) ;
	  sxExpExp (sx, aa, dict, "ILM_A", 't', "LIF_A", 't', trim) ;
	  sxExpExp (sx, aa, dict, "LIF_B", 't', "LIF_tot_B", 't', trim) ;
	  sxExpExp (sx, aa, dict, "ILM_tot_B", 't', "LIF_tot_B", 't', trim) ;
	}
      goto done ;
    }
  
  if (sx->venn)
    {
      sxExpVennPlot (sx, aa, dict, "HELdg_A", "ILM_A", "LIF_A",  trim) ;
      sxExpVennPlot (sx, aa, dict, "HEL_A", "ILM_A", "LIF_A",  trim) ;
      sxExpVennPlot (sx, aa, dict, "HEL_A", "ILMnoSt_A", "LIF_A",  trim) ;
      sxExpVennPlot (sx, aa, dict, "HEL_B", "ILM_B", "LIF_B",  trim) ;
      sxExpVennPlot (sx, aa, dict, "ILMnoSt_A", "ILM_A", "LIF_A",  trim) ;      
      goto done ;
    }

  if (sx->hammer) /* count genes supported by 1-n tags for 1-p platforms */
    {
      sxHammer (sx, aa, dict, 1) ;
      goto done ;
    }

  if (sx->nPlatform) /* count genes supported by 1-n tags for 1-p platforms */
    {
      sxExpOne (sx, aa, dict, 1) ;
      sxExpOne (sx, aa, dict, 2) ;
      sxExpOne (sx, aa, dict, 5) ;
      sxExpOne (sx, aa, dict, 10) ;
      sxExpOne (sx, aa, dict, 20) ;
      sxExpOne (sx, aa, dict, 50) ;
      sxExpOne (sx, aa, dict, 100) ;
      sxExpOne (sx, aa, dict, 200) ;
      sxExpOne (sx, aa, dict, 500) ;
      sxExpOne (sx, aa, dict, 1000) ;
      sxExpOne (sx, aa, dict, 2000) ;
      sxExpOne (sx, aa, dict, 5000) ;
      sxExpOne (sx, aa, dict, 10000) ;
      sxExpOne (sx, aa, dict, 20000) ;
      sxExpOne (sx, aa, dict, 50000) ;
      sxExpOne (sx, aa, dict, 100000) ;
      sxExpOne (sx, aa, dict, 0) ; /* histogram of all tags */
      goto done ;
    }
  if (1)
    for (s = sample ; *s ; s++)
      {
	printf ("\nCorrel") ;    
	for (t2 = tag ; *t2 ; t2++)
	  {
	    sprintf (buf2, "%s_%s", *t2, *s) ;
	    printf ("\t%s", buf2) ;
	  }
	for (t1= tag ; *t1 ; t1++)
	  {
	    sprintf (buf1, "%s_%s", *t1, *s) ;
	    printf ("\n%s", buf1) ;
	    u1 = 't' ; if (strstr(buf1,"TAQ")) u1 = 'c' ;
	    for (t2 = tag ; *t2 ; t2++)
	      {
		sprintf (buf2, "%s_%s", *t2, *s) ;
		u2 = 't' ; if (strstr(buf2,"TAQ")) u2 = 'c' ;
		printf ("\t%.3f", sxExpExp (sx, aa, dict, buf1, u1, buf2, u2, trim)) ;
	      }
	  }
	printf ("\n") ;
      }

  if (0)
    {
      printf ("\n\nHEL A\tAll\tUnique\tCorrected\n") ;
      printf ("All\t%.3f\t%.3f\t%.3f\n"
	      , sxExpExp (sx, aa, dict, "HH_A", 't', "HH_A", 't', trim)
	      , sxExpExp (sx, aa, dict, "HH_A", 't', "HH_A", 't', trim)
	      , sxExpExp (sx, aa, dict, "HH_A", 't', "HH_A", 'c', trim)
	      ) ;
      printf ("Unique\t%.3f\t%.3f\t%.3f\n"
	      , sxExpExp (sx, aa, dict, "HH_A", 't', "HH_A", 't', trim)
	      , sxExpExp (sx, aa, dict, "HH_A", 't', "HH_A", 't', trim)
	      , sxExpExp (sx, aa, dict, "HH_A", 't', "HH_A", 'c', trim)
	      ) ;
      printf ("Corrected\t%.3f\t%.3f\t%.3f\n"
	      , sxExpExp (sx, aa, dict, "HH_A", 'c', "HH_A", 't', trim)
	      , sxExpExp (sx, aa, dict, "HH_A", 'c', "HH_A", 't', trim)
	      , sxExpExp (sx, aa, dict, "HH_A", 'c', "HH_A", 'c', trim)
	      ) ;
      
      printf ("\n\nHEL B\tAll\tUnique\tCorrected\n") ;
      printf ("All\t%.3f\t%.3f\t%.3f\n"
	      , sxExpExp (sx, aa, dict, "HH_B", 't', "HH_B", 't', trim)
	      , sxExpExp (sx, aa, dict, "HH_B", 't', "HH_B", 't', trim)
	      , sxExpExp (sx, aa, dict, "HH_B", 't', "HH_B", 'c', trim)
	      ) ;
      printf ("Unique\t%.3f\t%.3f\t%.3f\n"
	      , sxExpExp (sx, aa, dict, "HH_B", 't', "HH_B", 't', trim)
	      , sxExpExp (sx, aa, dict, "HH_B", 't', "HH_B", 't', trim)
	      , sxExpExp (sx, aa, dict, "HH_B", 't', "HH_B", 'c', trim)
	      ) ;
      printf ("Corrected\t%.3f\t%.3f\t%.3f\n"
	      , sxExpExp (sx, aa, dict, "HH_B", 'c', "HH_B", 't', trim)
	      , sxExpExp (sx, aa, dict, "HH_B", 'c', "HH_B", 't', trim)
	      , sxExpExp (sx, aa, dict, "HH_B", 'c', "HH_B", 'c', trim)
	      ) ;
    }

  if (1)
    {
      printf ("\n\nLog(B/A)") ;
      for (t1= tag ; *t1 ; t1++)
	printf ("\t%s", *t1) ;
      for (t1= tag ; *t1 ; t1++)
	{
	   printf ("\n%s", *t1) ;
	   sprintf (buf1, "%s_%s", *t1, "A") ;
	   for (t2 = tag ; *t2 ; t2++)
	     {
	       sprintf (buf2, "%s_%s", *t2, "A") ;
	       printf ("\t%.3f", sxExpExp (sx, aa, dict, buf1, 'L', buf2, 'L', trim)) ;
	     }
	}
      printf ("\n\n") ;
    }

  if (0)
    {
      /* obsolete */
      printf ("\n\nLog(B/A)\tTAQ\tMAQC\tILMnoSt\tILM\tLIF\tHEL\tHELdg\tX\n") ;
      printf ("TAQ\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n"
	      , sxExpExp (sx, aa, dict, "TAQ_A", 'L', "TAQ_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "TAQ_A", 'L', "MAQC_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "TAQ_A", 'L', "ILMnoSt_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "TAQ_A", 'L', "ILM_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "TAQ_A", 'L', "LIF_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "TAQ_A", 'L', "HEL_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "TAQ_A", 'L', "HELdg_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "TAQ_A", 'L', "X_A", 'L', trim)
	      ) ;
      printf ("MAQC\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n"
	      , sxExpExp (sx, aa, dict, "MAQC_A", 'L', "TAQ_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "MAQC_A", 'L', "MAQC_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "MAQC_A", 'L', "ILMnoSt_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "MAQC_A", 'L', "ILM_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "MAQC_A", 'L', "LIF_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "MAQC_A", 'L', "HEL_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "MAQC_A", 'L', "HELdg_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "MAQC_A", 'L', "X_A", 'L', trim)
	      ) ;
      printf ("ILMnoSt\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n"
	      , sxExpExp (sx, aa, dict, "ILMnoSt_A", 'L', "TAQ_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "ILMnoSt_A", 'L', "MAQC_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "ILMnoSt_A", 'L', "ILMnoSt_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "ILMnoSt_A", 'L', "ILM_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "ILMnoSt_A", 'L', "LIF_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "ILMnoSt_A", 'L', "HEL_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "ILMnoSt_A", 'L', "HELdg_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "ILMnoSt_A", 'L', "X_A", 'L', trim)
	      ) ;
      printf ("ILM\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n"
	      , sxExpExp (sx, aa, dict, "ILM_A", 'L', "TAQ_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "ILM_A", 'L', "MAQC_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "ILM_A", 'L', "ILMnoSt_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "ILM_A", 'L', "ILM_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "ILM_A", 'L', "LIF_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "ILM_A", 'L', "HEL_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "ILM_A", 'L', "HELdg_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "ILM_A", 'L', "X_A", 'L', trim)
	      ) ;
      printf ("LIF\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n"
	      , sxExpExp (sx, aa, dict, "LIF_A", 'L', "TAQ_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "LIF_A", 'L', "MAQC_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "LIF_A", 'L', "ILMnoSt_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "LIF_A", 'L', "ILM_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "LIF_A", 'L', "LIF_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "LIF_A", 'L', "HEL_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "LIF_A", 'L', "HELdg_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "LIF_A", 'L', "X_A", 'L', trim)
	      ) ;
      printf ("HEL\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n"
	      , sxExpExp (sx, aa, dict, "HEL_A", 'L', "TAQ_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "HEL_A", 'L', "MAQC_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "HEL_A", 'L', "ILMnoSt_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "HEL_A", 'L', "ILM_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "HEL_A", 'L', "LIF_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "HEL_A", 'L', "HEL_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "HEL_A", 'L', "HELdg_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "HEL_A", 'L', "X_A", 'L', trim)
	      ) ;
      printf ("HELdg\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n"
	      , sxExpExp (sx, aa, dict, "HELdg_A", 'L', "TAQ_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "HELdg_A", 'L', "MAQC_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "HELdg_A", 'L', "ILMnoSt_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "HELdg_A", 'L', "ILM_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "HELdg_A", 'L', "LIF_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "HELdg_A", 'L', "HEL_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "HELdg_A", 'L', "HELdg_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "HELdg_A", 'L', "X_A", 'L', trim)
	      ) ;
      printf ("X\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n"
	      , sxExpExp (sx, aa, dict, "X_A", 'L', "TAQ_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "X_A", 'L', "MAQC_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "X_A", 'L', "ILMnoSt_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "X_A", 'L', "ILM_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "X_A", 'L', "LIF_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "X_A", 'L', "HEL_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "X_A", 'L', "HELdg_A", 'L', trim)
	      , sxExpExp (sx, aa, dict, "X_A", 'L', "X_A", 'L', trim)
	      ) ;
    }

 done:
  ac_free (h) ;
  return ;
} /* seExp */

/*************************************************************************************/
/******************************* Saturation ******************************************/
/*************************************************************************************/
/* count how many targets are seen with 0 to 100% of the tags being considered */
static void sxSaturation (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  DICT *dict = dictHandleCreate (100000, h) ;
  int ii, n, nn = 0, n0 = sx->saturation, t ;
  long int line = 0 ;
  char *cp ;
  ACEIN ai = 0 ;
  KEYSET ks = keySetHandleCreate (h) ;
  Array aa, a ;
  int *s, ss[] = {1,2,5,10,20,50,100,200,500,1000,0} ;

  aa = arrayHandleCreate (300000, Array, h) ;
  ai = aceInCreate (sx->hitFile, 0, h) ;

  while (ai && aceInCard (ai))
    {
      cp = aceInWord (ai) ;
      if (cp)
	{
	  line++ ;
	  dictAdd (dict, cp, &t) ;
	  a = array (aa, t, Array) ;
	  if (! a)
	    a = array (aa, t, Array) = arrayHandleCreate (n0, int, h) ; 
	  nn = n0 * randfloat () * 51.916/50.0 ;
	  nn = (nn % n0) ;
	  /* count the hits to target t in quantile nn */
	  array (a, nn, int)++ ; 
	}
    }

  /* count the cumul of the hist to each target up to each bin */
  for (t = dictMax (dict) ; t >= 1 ; t--)
    {
      a = array (aa, t, Array) ;
      for (n = nn = 0 ; nn < n0 ; nn++)
	{
	  n += array (a, nn, int) ; 
	  array (a, nn, int) = n ; 
	}
    }

  /* for each threshold, count the number of target seen in that bin */
  for (t = dictMax (dict) ; t >= 1 ; t--)
    {
      a = array (aa, t, Array) ;
      for (nn = 0 ; nn < n0 ; nn++)
	{
	  n = array (a, nn, int) ; 
	  for (ii = 0, s = ss ; *s ; ii++, s++)
	    if (n >= *s)
	      keySet (ks, ii*n0 + nn)++ ;
	}
    }
  printf ("bin") ;
  for (ii = 0, s = ss ; *s ; ii++, s++)
    printf ("\t%d", *s) ;
  printf ("\n0") ;
  for (ii = 0, s = ss ; *s ; ii++, s++)
    printf ("\t0") ;
  for (nn = 0 ; nn < n0 ; nn++)
    {
      printf ("\n%d", nn+1) ;
      for (ii = 0, s = ss ; *s ; ii++, s++)
	printf ("\t%d", keySet (ks, ii*n0 + nn)) ;
    }
  printf ("\nTotal\t%u targets\t%ld hits\n", dictMax (dict), line) ;
  ac_free (h) ;
  return ;
} /* sxSaturation */

/*************************************************************************************/
/******************************* New Introns******************************************/
/*************************************************************************************/

typedef struct daStruct { BOOL isDown ; int chrom, a1, exon, hook, n, nT ; KEY SL ; } DA ;

static void showDA (char *title, DICT *dict, Array hits)
{
  DA *up ;
  unsigned int i = hits ? arrayMax (hits) : 0 ;
  char *cp, *cq = 0 ;
  AC_HANDLE h = ac_new_handle () ;

  if (title) 
    printf ("%s\n", title) ;
  if (i)
    for (i = 0 ; i < arrayMax (hits) ; i++)
      {
	up = arrp (hits, i, DA) ;
	cp = hprintf (h, "%s %s %d %d %s %s\n"
		      , up->isDown ? "down" : "  up"
		      , dictName (dict, up->chrom)
		      , up->n, up->a1
		      , dictName (dict, up->exon)
		      , dictName (dict, up->hook)
		      ) ; 
	if (! cq || strcmp (cp, cq))
	  printf ("%s", cp) ;
	cq = cp ;
      }
  ac_free (h) ;
} /* showDA */

/*************************************************************************************/

static int DAorderHook (const void *a, const void *b)
{
  const DA *up = (const DA *) a, *vp = (const DA *) b ;
  if (up->isDown < vp->isDown) return -1 ;
  if (up->isDown > vp->isDown) return 1 ;
  if (up->chrom < vp->chrom) return -1 ;
  if (up->chrom > vp->chrom) return 1 ;
  if (up->a1 < vp->a1) return -1 ;
  if (up->a1 > vp->a1) return 1 ;
  if (up->hook < vp->hook) return -1 ;  /* by hook */
  if (up->hook > vp->hook) return 1 ;
  if (up->n < vp->n) return 1 ;  /* most support first */
  if (up->n > vp->n) return -1 ;
  
  return 0 ;
} /* DAorder */

/*************************************************************************************/

static int DAorderN (const void *a, const void *b)
{
  const DA *up = (const DA *) a, *vp = (const DA *) b ;
  if (up->isDown < vp->isDown) return -1 ;
  if (up->isDown > vp->isDown) return 1 ;
  if (up->chrom < vp->chrom) return -1 ;
  if (up->chrom > vp->chrom) return 1 ;
  if (up->a1 < vp->a1) return -1 ;
  if (up->a1 > vp->a1) return 1 ;
  if (up->n < vp->n) return 1 ;  /* most support first */
  if (up->n > vp->n) return -1 ;
  if (up->hook < vp->hook) return -1 ;  /* by hook */
  if (up->hook > vp->hook) return 1 ;
  
  return 0 ;
} /* DAorder */

/*************************************************************************************/
/* collate the number of tags supporting the same jump */
static int sxNewIntronsCompress (Array donors, int limit)
{
  int ii, jj = 0 ;
  DA *up, *vp ;

  for (ii = 0, up = arrp (donors, 0, DA) ; ii < arrayMax (donors) ; ii++, up++)
    {
      if (! up->n) continue ;
      for (jj = ii + 1, vp = up + 1 ; jj < arrayMax (donors) ; jj++, vp++)
	if (up->isDown == vp->isDown &&
	    up->chrom == vp->chrom &&
	    up->a1 == vp->a1 &&
	    up->hook == vp->hook &&
	    up->SL == vp->SL
	    )
	  { up->n += vp->n ; up->nT += vp->nT ; vp->n = 0 ; vp->nT = 0 ; }
	else
	  break ;
    }

  /* keep happy few */
  for (ii = jj = 0, vp = up = arrp (donors, 0, DA); ii < arrayMax (donors) ; up++, ii++)
    {
      if (up->n >= limit)
	{
	  if (jj < ii)
	    *vp = *up ;
	  vp++ ; jj++ ;
	}
    }
  arrayMax (donors) = jj ;
  return jj ;
}  /* sxNewIntronsCompress */

/********************************************************************/
/* bp equivalent of the oligo complexity
 *
 * input is dna, x1/x2 in bio-coords
 */
int hookEntropy (const char *hook)
{
  int ss = 0 ;
  int na = 0, nt = 0, ng = 0, nc = 0, nn ;
  const char *ccp ;
  static BOOL oldNn = -1 ;
  static Array ee = 0 ;
  /*  count all letters
   */
  if (1)
    {
      ccp = hook - 1 ;
      while (*++ccp)
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
  return ss ;
} /* hookEntropy */

/*************************************************************************************/
/*************************************************************************************/

typedef struct dbliStruct { int typea, typeb, typec, a1, a2, b1, b2, c1, c2, lna, lnc, beta1, beta2, chrom, donorSupport, acceptorSupport, support ; BOOL isDown ;} DBLI ;

/*************************************************************************************/

static void showDBLI (BigArray introns)
{
  DBLI *up ;
  int i ;

  for (i = 0 ; i < bigArrayMax(introns) ; i++)
    {
      up = bigArrp (introns, i, DBLI) ;
      fprintf (stderr, "%d\ttype %d:%d:%d\ta:%d:%d\tb:%d:%d\tc:%d:%d\tsup:%d:%d:%d\n"
	       , i
	       , up->typea, up->typeb, up->typec
	       , up->a1, up->a2, up->b1, up->b2, up->c1, up->c2
	       , up->donorSupport, up->acceptorSupport, up->support
	       ) ;
    }
}

static int dbliA1Order (const void *a, const void *b)
{
  const DBLI *up = (const DBLI *) a, *vp = (const DBLI *) b ;
  int n ;

  n = up->chrom - vp->chrom ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  n = up->isDown - vp->isDown ; if (n) return n ;
  n = up->b1 - vp->b1 ; if (n) return n ;
  n = up->beta1 - vp->beta1 ; if (n) return n ;
  n = up->b2 - vp->b2 ; if (n) return n ;
  n = up->beta2 - vp->beta2 ; if (n) return n ;
  n = up->c1 - vp->c1 ; if (n) return n ;
  n = up->c2 - vp->c2 ; if (n) return n ;

  return 0 ;
}  /* dbliA1Order */

/*************************************************************************************/
/* do NOT add type in the order function because of deUnoCompress */
static int dbliB1Order (const void *a, const void *b)
{
  const DBLI *up = (const DBLI *) a, *vp = (const DBLI *) b ;
  int n ;

  n = up->chrom - vp->chrom ; if (n) return n ;
  n = up->isDown - vp->isDown ; if (n) return n ; /* mieg: 2016_08_20: moved isDown to here to please deUnoIntronsCompress, we may need to verify this order funct when used for rearaengments */
  n = up->b1 - vp->b1 ; if (n) return n ;
  n = up->beta1 - vp->beta1 ; if (n) return n ;
  n = up->b2 - vp->b2 ; if (n) return n ;
  n = up->beta2 - vp->beta2 ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  n = up->c1 - vp->c1 ; if (n) return n ;
  n = up->c2 - vp->c2 ; if (n) return n ;

  return 0 ;
}  /* dbliB1Order */

/*************************************************************************************/

static int dbliB2Order (const void *a, const void *b)
{
  const DBLI *up = (const DBLI *) a, *vp = (const DBLI *) b ;
  int n ;

  n = up->chrom - vp->chrom ; if (n) return n ;
  n = up->b2 - vp->b2 ; if (n) return n ;
  n = up->beta2 - vp->beta2 ; if (n) return n ;
  n = up->b1 - vp->b1 ; if (n) return -n ; /* shortest intron first */
  n = up->beta1 - vp->beta1 ; if (n) return n ;
  n = up->c1 - vp->c1 ; if (n) return n ;
  n = up->c2 - vp->c2 ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  n = up->isDown - vp->isDown ; if (n) return n ;

  return 0 ;
}  /* dbliB2Order */

/*************************************************************************************/

static Array sxParseFasta (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (sx->target, FALSE, h) ;
  int i = 0, state = 0 ;
  char *cp, c, c1 ;
  Array dna = 0 ;
  
  if (ai)
    {
      dna = arrayHandleCreate (1000000, unsigned char, sx->h) ;
      while (aceInCard (ai))
	{
	  cp = aceInWord (ai) ;
	  if (! cp || *cp == '#')
	    continue ;
	  if (*cp == '>')
	    {
	      state = 1 ; 
	      continue ;
	    }
	  if (!state)
	    continue ;
	  
	  while ((c = *cp++))
	    {
	      if ((c1 = dnaEncodeChar[((int)c) & 0x7f]))        /* accepted base codes */
		array(dna,i++,char) = c1 ;
	      else
		switch(c)
		  {			/* accepted tabulations */
		  case '-': case 'x':  case 'X': /* x is vector masking */
		    array(dna,i++,char) = 0 ;
		    break ;
		  case '*':  /* phil green padding character, ignore */
		  case '\n': case '\t': case ' ':
		  case '0': case '1': case '2': case '3': case '4': 
		  case '5': case '6': case '7': case '8': case '9': 
		    break;
		  default:		/* error in the input file */
		    c1 = 0xff ; goto abort ;
		  }
	    }
	}
    }
 abort:

  ac_free (h) ;
  return dna ;
} /* sxParseFasta  */

/*************************************************************************************/

static BOOL sxGetTargetFasta (SX *sx)
{
  Array dnaD = 0, dnaR = 0 ;

  if (sx->target)
    {
      dnaD = sxParseFasta (sx) ;
      if (dnaD)
	{
	  dnaR = dnaHandleCopy (dnaD, sx->h) ;
	  reverseComplement (dnaR) ;
	  dnaDecodeArray (dnaD) ;
	  dnaDecodeArray (dnaR) ;
	  sx->chromDnaD = dnaD ;
	  sx->chromDnaR = dnaR ;
	  return TRUE ;
	}
    }
  ac_free (dnaD) ;
  ac_free (dnaR) ;

  return FALSE ;
} /* sxGetTargetFasta */

/*************************************************************************************/

static void sxNewDoubleIntronsExport (SX *sx, BigArray dblis, DICT *dict, BOOL debug)
{
  int ii, ns, u1, u2 ;
  int dnaMax ; 
  DBLI *dbli ;
  Array dna, dnaD = 0, dnaR = 0 ;
  char cc, *cp, aaaBuf[1000], tttBuf[1000], slbuf[40] ;
  char *slu[] = { "ggtttaattacccaagtttgag",     /* SL1 */
		  "ggttttaacccagttactcaag" ,    /* SL2 */
		  "ggttttaacccagttaaccaag" ,    /* SL3 */
		  "ggttttaacccagtttaaccaag" ,    /* SL4 */
		  "ggttttaacccagttaccaag" ,    /* SL5 */
		  "ggtttaaaacccagttaccaag" ,    /* SL6 */
		  "ggttttaacccagttaattgag" ,    /* SL7 */
		  "ggtttttacccagttaaccaag" ,    /* SL8 */
		  "ggtttatacccagttaaccaag" ,    /* SL9 */
		  "ggttttaacccaagttaaccaag" ,    /* SL10 */
		  "ggttttaaccagttaactaag" ,    /* SL11 */
		  "ggttttaacccatataaccaag" ,    /* SL12 */
		  "gtttttaacccagttactcaag" ,    /* SL13 */
		  "ggtttttaacccagttactcaag" ,    /* SL14 */
		  0 } ;


  bigArraySort (dblis, dbliA1Order) ;
  memset (aaaBuf, 'A', sizeof (aaaBuf)) ;
  memset (tttBuf, 'T', sizeof (tttBuf)) ;

  sxGetTargetFasta (sx) ;
  dnaD = sx->chromDnaD ; dnaR = sx->chromDnaR ;
  dnaMax = dnaD ? arrayMax (dnaD) : 0 ;
  freeOutf ("#Chromosome\tStrand\tTypea\tTypeb\tTypec\ta1\ta2\tb1\tb2\tc1\tc2\tDNAa\tDNAb\tDNAc\tSupport\n") ;
  for (ii = 0, dbli = bigArrp (dblis, 0, DBLI) ; ii < bigArrayMax (dblis) ; dbli++, ii++)
    {
      int isIntron = 0 ;
      dna = dbli->isDown ? dnaD : dnaR ;

      if (dbli->a1 > dbli->a2 || dbli->a2 > dbli->b1 || dbli->b1 > dbli->b2 || dbli->b2 > dbli->c1 || dbli->c1 > dbli->c2)
	{
	  fprintf( stderr, "INVERSION %s\t%s\typea=%d\tb=%d\tc=%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
		   , dbli->chrom ? dictName (dict, dbli->chrom) : "-"
		   , dbli->isDown ? "Forward" : "Reverse"
		   , dbli->typea, dbli->typeb , dbli->typec
		   , dbli->a1, dbli->a2, dbli->b1, dbli->b2, dbli->c1, dbli->c2
		   ) ;
	  continue ;
	}
      if (! dbli->isDown)
	{
	  int x ;
	  x = dbli->a1 ; dbli->a1 = dbli->c2 ;  dbli->c2 = x ;
	  x = dbli->a2 ; dbli->a2 = dbli->c1 ;  dbli->c1 = x ;
	  x = dbli->b2 ; dbli->b2 = dbli->b1 ;  dbli->b1 = x ;
	  x = dbli->typea ; dbli->typea = dbli->typec ; dbli->typec = x ;
	  x = dbli->lna ; dbli->lna = dbli->lnc ; dbli->lnc = x ;
	}

      freeOutf ("%s", dbli->chrom ? dictName (dict, dbli->chrom) : "-") ;
      freeOutf ("\t%s", dbli->isDown ? "Forward" : "Reverse") ;
      switch (dbli->typea)
	{
	case 0: 
	  freeOutf ("\tOpen") ;
	  break ;
	case 1:
	  freeOutf ("\tExon") ;
	  isIntron++ ;
	  break ;
	case 2:
	  freeOutf ("\tIntron") ;
	  break ;
	case 3:
	  freeOutf ("\tpA") ;
	  break ;
	default:
	  freeOutf ("\tSL%d", dbli->typea - 3) ;
	  break ;
	}
      switch (dbli->typeb)
	{
	case 0: 
	  freeOutf ("\tOpen") ;
	  break ;
	case 1:
	  freeOutf ("\tExon") ;
	  break ;
	case 2:
	  freeOutf ("\tIntron") ;
	  isIntron++ ;
	  break ;
	case 3:
	  freeOutf ("\tpA") ;
	  break ;
	default:
	  freeOutf ("\tSL%d", dbli->typeb - 3) ;
	  break ;
	}

      switch (dbli->typec)
	{
	case 0: 
	  freeOutf ("\tOpen") ;
	  break ;
	case 1:
	  freeOutf ("\tExon") ;
	  isIntron++ ;
	  break ;
	case 2:
	  freeOutf ("\tIntron") ;
	  break ;
	case 3:
	  freeOutf ("\tpA") ;
	  break ;
	default:
	  freeOutf ("\tSL%d", dbli->typec - 3) ;
	  break ;
	}


      freeOutf ("\t%d\t%d\t%d\t%d\t%d\t%d", dbli->a1, dbli->a2, dbli->b1, dbli->b2, dbli->c1, dbli->c2) ;

      if (dbli->a1 < dnaMax && dbli->c2 < dnaMax)
	{
	  switch (dbli->typea)
	    {
	    case 0: /* open ended */
	      u1 = dbli->isDown ? dbli->a1 : dnaMax - dbli->a1 + 1 ;
	      u2 = dbli->isDown ? dbli->a2 : dnaMax - dbli->a2 + 1 ;
	      if (u2 > u1 + 4)
		{
		  cp = arrp (dna, u2, char) ;
		  cc = *cp ; *cp = 0 ;
		  if (debug)freeOutf ("\tDNA_A0 %d_%d=%d", u1, u2, u2-u1+1) ;
		  else freeOutf ("\t%s", arrp (dna, u1 - 1, char)) ;
		  *cp = cc ;
		}
	      else
		freeOutf ("\t") ; 
	      break ;
	    case 1:
	    case 2:
	      u1 = dbli->isDown ? dbli->a1 : dnaMax - dbli->a1 + 1 ;
	      u2 = dbli->isDown ? dbli->a2 : dnaMax - dbli->a2 + 1 ;
	      cp = arrp (dna, u2, char) ;
	      cc = *cp ; *cp = 0 ;
	      if (debug)freeOutf ("\tDNA_A2 %d_%d=%d", u1, u2, u2-u1+1) ; 
	      else freeOutf ("\t%s", arrp (dna, u1 - 1, char)) ;
	      *cp = cc ;
	      break ;
	    case 3:
	      cp = tttBuf + dbli->lna ;
	      *cp = 0 ;
	      freeOutf ("\t%s", tttBuf) ;
	      *cp = 'T' ;
	      break ;
	    default:
	      ns = dbli->typea - 4 ;
	      cp = slu[ns] + strlen (slu[ns]) - dbli->lna ;
	      strcpy (slbuf, cp) ;
	      cp = slbuf - 1 ; while (*++cp) { *cp = ace_upper(*cp) ; } 
	      freeOutf ("\t%s", slbuf) ;
	      break ;
	    }
	  
	  switch (dbli->typeb)
	    {
	    case 1:
	      u1 = dbli->isDown ? dbli->b1 : dnaMax - dbli->b1 + 1 ;
	      u2 = dbli->isDown ? dbli->b2 : dnaMax - dbli->b2 + 1 ;
	      cp = arrp (dna, u2, char) ;
	      cc = *cp ; *cp = 0 ;
	      if (debug)freeOutf ("\tDNA_B1 %d_%d=%d", u1, u2, u2-u1+1) ; 
	      else freeOutf ("\t%s", arrp (dna, u1 - 1, char)) ;
	      *cp = cc ;
	      break ;
	    case 2:
	      freeOutf ("\t") ; 
	      break ;
	    case 3:
	      cp = aaaBuf + dbli->lnc ;
	      *cp = 0 ;
	      freeOutf ("\t%s", aaaBuf) ;
	      *cp = 'A' ;
	      break ;
	    default:
	      ns = dbli->typeb - 4 ;
	      cp = slu[ns] + strlen (slu[ns]) - dbli->lna ;
	      strcpy (slbuf, cp) ;
	      cp = slbuf - 1 ; while (*++cp) { *cp = ace_upper(*cp) ; } 
	      freeOutf ("\t%s", slbuf) ;
	      break ;
	    }

	  switch (dbli->typec)
	    {
	    case 0: /* open ended */
	      u1 = dbli->isDown ? dbli->c1 : dnaMax - dbli->c1 + 1 ;
	      u2 = dbli->isDown ? dbli->c2 : dnaMax - dbli->c2 + 1 ;
	      if (u2 > u1 + 4)
		{
		  cp = arrp (dna, u2, char) ;
		  cc = *cp ; *cp = 0 ;
		  if (debug)freeOutf ("\tDNA_C0 %d_%d=%d", u1, u2, u2-u1+1) ; 
		  else freeOutf ("\t%s", arrp (dna, u1 - 1, char)) ;
		  *cp = cc ;
		}
	      else
		freeOutf ("\t") ; 
	      break ;
	    case 1:
	    case 2:
	      u1 = dbli->isDown ? dbli->c1 : dnaMax - dbli->c1 + 1 ;
	      u2 = dbli->isDown ? dbli->c2 : dnaMax - dbli->c2 + 1 ;
	      cp = arrp (dna, u2, char) ;
	      cc = *cp ; *cp = 0 ;
	      if (debug)freeOutf ("\tDNA_C2 %d_%d=%d", u1, u2, u2-u1+1) ;
	      else freeOutf ("\t%s", arrp (dna, u1 - 1, char)) ;
	      *cp = cc ;
	      break ;
	    case 3:
	      cp = aaaBuf + dbli->lnc ;
	      *cp = 0 ;
	      freeOutf ("\t%s", aaaBuf) ;
	      *cp = 'A' ;
	      break ;
	    default:
	      cp = aaaBuf + dbli->lnc ;
	      *cp = 0 ;
	      freeOutf ("\tSLXXX") ;
	      *cp = 'A' ;
	      break ;
	    }
	}
      else
	freeOutf ("\t-\t-\t-") ;
      freeOutf ("\t%d", dbli->support) ;
      if (isIntron == 3)
	{
	  u1 = dbli->isDown ? dbli->a2 + 1 : dnaMax - dbli->a2 + 2 ;
	  u2 = u1 + 1 ;
	  cp = arrp (dna, u2, char) ;
	  cc = *cp ; *cp = 0 ;
	  freeOutf ("\t%s", arrp (dna, u1 - 1, char)) ;
	  *cp = cc ;
	  u1 = dbli->isDown ? dbli->c1 - 2 : dnaMax - dbli->c1 - 1 ;
	  u2 = u1 + 1 ;
	  cp = arrp (dna, u2, char) ;
	  cc = *cp ; *cp = 0 ;
	  freeOutf ("_%s", arrp (dna, u1 - 1, char)) ;
	  *cp = cc ;
	}
      freeOutf ("\n") ;
    }
} /* sxNewDoubleIntronsExport */

/*************************************************************************************/
/* compress and merge support */
static void sxNewDoubleIntronsDbliCompress (BigArray aa)
{
  long int ii, jj, iMax = bigArrayMax (aa), support, s2,  dSupport, ds2, aSupport, as2 ;
  DBLI *up, *vp ;
  BOOL ok ;

  for (ii = 0, up = bigArrp (aa, 0, DBLI) ; ii < iMax ; up++,  ii++) 
    {
      if (up->support == 0)
	continue ;
      support = up->support ;
      dSupport = up->donorSupport ;
      aSupport = up->acceptorSupport ;
      up->support = up->donorSupport = up->acceptorSupport = 0 ;
      for (jj = ii + 1, vp = up + 1 ; jj < iMax ; vp++,  jj++) 
	{
	  s2 = vp->support ; ds2 =  vp->donorSupport ; as2 = vp->acceptorSupport ;
	  vp->support = vp->donorSupport = vp->acceptorSupport = 0 ;
	  ok = TRUE ;
	  {
	    register int k ;
	    register const char *x, *y ;
	    for (x = (char *)up, y = (char *)vp, k = sizeof (DBLI) ; k > 0 ; x++, y++, k--)
	      if (*x != *y)
		ok = FALSE ;
	  }
	  if (ok)
  	    { support += s2 ; dSupport += ds2 ; aSupport += as2 ; }
	  else
	    {
	      vp->support = s2 ; vp->donorSupport = ds2 ; vp->acceptorSupport = as2 ;
	      break ;
	    }
	}
      up->support = support ; up->donorSupport = dSupport ; up->acceptorSupport = aSupport ;
    }

  for (ii = jj = 0, up = vp = bigArrp (aa, 0, DBLI) ; ii < iMax ; up++,  ii++) 
    {
      if (up->support)
	{
	  if (jj < ii) 
	    *vp = *up ;
	  jj++ ; vp++ ;
	}
    }
  bigArrayMax (aa) = jj ;
}  /* sxNewDoubleIntronsDbliCompress */

/*************************************************************************************/
/* using the intron-exon-intron segments, contruct the exon-intron-exon segments 
 * intronsB1/2 contain the open edge exon-intron-exon candidates 
 */
static void sxNewDoubleIntronsBridge (SX *sx, BigArray dblis, BigArray intronsB1, BigArray intronsB2)
{
  long int ii, jj ;
  DBLI *up, *vp ;
  BigArray aa = 0 ; 

  /* treat the open ended introns */
  aa = bigArrayCopy (intronsB1) ;
  for (ii = bigArrayMax (intronsB1), jj = 0 ; jj < bigArrayMax (intronsB2) ; ii++, jj++)
    {
      up = bigArrayp (aa, ii, DBLI) ;
      vp = bigArrayp (intronsB2, jj, DBLI) ;
      *up = *vp ;
    }
  bigArraySort (aa, dbliB1Order) ;
  /* merge the 2 left and right ends */
  for (ii = 0, up = bigArrp (aa, 0, DBLI) ; ii < bigArrayMax (aa) ; up++, ii++)
    {
      if (! up->typeb) continue ;
      for (jj = ii + 1, vp = up + 1 ; jj < bigArrayMax (aa) ; vp++, jj++)
	{
	  if (up->chrom != vp->chrom || 
	      up->beta1 != vp->beta1 || up->b1 != vp->b1 || 
	      up->beta2 != vp->beta2 || up->b2 != vp->b2
	      )
	    break ;
	  if (vp->typea && ! up->typea)
	    { up->typea = vp->typea ; up->a1 = vp->a1 ; up->a2 = vp->a2 ; }
	  if (vp->typec && ! up->typec)
	    { up->typec = vp->typec ; up->c1 = vp->c1 ; up->c2 = vp->c2 ; }
	  /*   up->support += vp->support ; stupid: intronsB2 is a copy of intronsB1, same supports ! */
	  up->donorSupport += vp->donorSupport ;
	  up->acceptorSupport += vp->acceptorSupport ;
	  if (vp->a1 < up->a1) up->a1 = vp->a1 ;
	  if (vp->c2 > up->c2) up->c2 = vp->c2 ;
	  vp->typeb = 0 ;
	}
    }
  /* compress and merge in dblis */
  for (ii = 0, up = bigArrayp (aa, 0, DBLI), jj = bigArrayMax (dblis) ; ii < bigArrayMax (aa) ; up++, ii++)
    {
      if (! up->typeb) continue ;
      if (! up->typea && !up->typec) continue ;
      vp = bigArrayp (dblis, jj++, DBLI) ;
      if (vp != up) *vp = *up ;
    }
  /* fix the coordinates */
  for (ii = 0, up = bigArrp (dblis, 0, DBLI) ;  ii < bigArrayMax (dblis) ; up++, ii++)
    {
      if (up->isDown && up->typea >= 4 && up->typeb == 1)    /* SL  exon ... */
	{ up->a2 = up->b1 ; up->a1 = up->b1 - 1 ; }
      if (! up->isDown && up->typec >= 4 && up->typeb == 1)  /* ... exon SL */
	{ up->c1 = up->b2 ; up->c2 = up->b2 + 1 ; }
      if (up->isDown && up->typeb == 1 && up->typec == 3)    /* ... exon pA */
	up->c1 = up->b2 ;
      if (! up->isDown && up->typeb == 1 && up->typea == 3)  /* pA  exon ...*/
	up->a2 = up->b1 ;

      if (up->isDown && up->typeb == 3 && up->typea == 1)    /* exon pA ... */
	up->a2 = up->b1 ;
      if (! up->isDown && up->typeb == 3 && up->typec == 1)  /* ...  pA exon */
	up->c1 = up->b2 ;
      
      if (up->isDown && up->typeb >= 4)
	up->b2++ ;
      if (! up->isDown && up->typeb >= 4)
	up->b1-- ;

      if (up->isDown && up->typeb >= 4 && up->typea == 0)    /* open SL ... */
	up->a1 = up->a2 = up->b1 ;
      if (! up->isDown && up->typeb >= 4 && up->typec == 0)  /* ...  SL open  */
	up->c1 = up->c2 = up->b2 ;
      if (! up->isDown && up->typeb >= 4 && up->typea == 1)  /* exon SL ... */
	up->a2 = up->b1 ;
      if (up->isDown && up->typeb >= 4 && up->typec == 1)    /* ...  SL exon */
	up->c1 = up->b2 ; 
    }

  bigArrayDestroy (aa) ;
  return ;
} /* sxNewDoubleIntronsBridge */

/*************************************************************************************/

static long int sxNewDoubleIntronsGetDeUnoIntronsCompress (BigArray introns)
{
  DBLI *up, *vp, *wp ;
  long int ii, iiMax = bigArrayMax (introns), jj, kk ;

  bigArraySort (introns, dbliB1Order) ;
  
  /* merge identical introns, cumulating the support */
  for (up = vp = bigArrayp (introns, 0, DBLI), ii = jj = 0 ; ii < iiMax ; ii++, up++)
    {
      if (! up->support) continue ;
      for (wp = up + 1, kk = ii + 1 ; kk < iiMax ; kk++, wp++)
	{
	  if (wp->chrom != up->chrom ||
	      wp->b1 > up->b1 ||
	      wp->b2 > up->b2 ||
	      wp->isDown != up->isDown
	      ) break ;
	  if (! wp->support || up->typea != wp->typea || up->typeb != wp->typeb || up->typec != wp->typec) continue ;
	  if (wp->a1 < up->a1) up->a1 = wp->a1 ;
	  if (wp->c2 > up->c2) up->c2 = wp->c2 ;
	  if (up->typeb < wp->typeb) { up->typeb = wp->typeb ; up->isDown = wp->isDown ; }
	  up->support += wp->support ;
	  up->donorSupport += wp->donorSupport ;
	  up->acceptorSupport += wp->acceptorSupport ;
	  wp->support = 0 ;
	}
      if (vp < up) *vp = *up ;

      jj++ ; vp++ ;
    }
  bigArrayMax (introns) = jj ;
  return jj ;
} /*sxNewDoubleIntronsGetDeUnoIntronsCompress */

/*************************************************************************************/
/* parse the de duo intron  file
   INTRON Reverse CHROMOSOME_II 000011647 000011602 DONOR:tctctctgtgctttgttaagcaa 1 ACCEPTOR:agttttcgaacagctgtgcactt 1 dx1=2  bestdx=0
   INTRON Forward CHROMOSOME_II 000013061 000013103 DONOR:cttccatcgaaatagtagtttg 5 ACCEPTOR:gggtcaatctggaaatcggtg 5 dx1=5  bestdx=1
*/
static BigArray sxNewDoubleIntronsGetDeDuoIntrons (SX *sx, DICT *dict, AC_HANDLE h)
{
  ACEIN ai = 0 ;
  DBLI *dbli ;
  BigArray introns = 0 ;
  long int nn = 0 ;
  int b1, b2, chrom, beta1, beta2, lna = 0, lnc = 0 ;
  int minA1 = sx->minX ;
  int maxA1 = sx->maxX ;
  BOOL isDown = TRUE, isDown2 = TRUE ;
  char *cp, *cq ;
  int type ;
  int donorSupport, acceptorSupport ;

  introns = bigArrayHandleCreate (10000, DBLI, h) ;
  ai = 0 ;
  if (sx->sxxNewIntronsFileName)
    ai = aceInCreate (sx->sxxNewIntronsFileName, 0,  h) ;
  while (ai && aceInCard (ai))
    {
      type = 0 ;
      if (! (cp = aceInWord (ai)))
	continue ;
      if (! strcmp (cp, "INTRON"))
	type = 2 ;
      else if (! strcmp (cp, "pA"))
	type = 3 ;
      else if (! strcmp (cp, "pT"))
	type = 3 ;
      else if (! strncmp (cp, "SL", 2))
	{
	  cq = cp + 2 ;
	  while (*cq >= '0' && *cq <= '9')
	    type = 10 * type + (*cq++ - '0') ;
	  type += 3 ; /* so SL1 is represented as type == 4 */
	}
      if (! type)
	continue ;

      donorSupport =  acceptorSupport = 0 ;
      aceInStep(ai,'\t') ;
      if (! (cp = aceInWord (ai)))
	continue ;
      if (! strcmp (cp, "Forward"))
	isDown = TRUE ;
      else if (! strcmp (cp, "Reverse"))
	isDown = FALSE ;
      else
	continue ;

      aceInStep(ai,'\t') ;
      cp = aceInWord (ai) ;    /* chromosome */
      if (! cp)  continue ;
      if (sx->sxxChromosome && strcmp (sx->sxxChromosome, cp))
	continue ;
      dictAdd (dict, cp , &chrom) ;

      aceInStep(ai,'\t') ;
      if (!aceInInt (ai, &b1))     /* position of first matching base */
	continue ;

      aceInStep(ai,'\t') ;
      cp = aceInWord (ai) ;    /* chromosome */
      if (! cp)  continue ;
      if (strcmp (cp, dictName(dict, chrom)))
	continue ;

      aceInStep(ai,'\t') ;
      if (! (cp = aceInWord (ai)))
	continue ;
      if (! strcmp (cp, "Forward"))
	isDown2 = TRUE ;
      else if (! strcmp (cp, "Reverse"))
	isDown2 = FALSE ;
      else
	continue ;
      if (isDown != isDown2)
	continue ;

      aceInStep(ai,'\t') ;
      if (!aceInInt (ai, &b2))     /* position of last matching base */
	continue ;
      
      if (! isDown) 
	{ int b0 = b1 ; b1 = b2 ; b2 = b0 ; }

      if ((minA1 && b1 < minA1) || (maxA1 && b1 > maxA1))
	continue ;
      if ((minA1 && b2 < minA1) || (maxA1 && b2 > maxA1))
	continue ;

      aceInStep(ai,'\t') ;
      cp = aceInWord (ai) ;    /* length */
      aceInStep(ai,'\t') ;
      cp = aceInWord (ai) ;    /* gt_ag */
      aceInStep(ai,'\t') ;
      cp = aceInWord (ai) ;    /* SUPPORT */
     aceInStep(ai,'\t') ;
      cp = aceInWord (ai) ;    /* - */

      aceInStep(ai,'\t') ;
      if (!aceInInt (ai, &donorSupport))     /* Number of supports of donor */
	continue ;
      aceInStep(ai,'\t') ;
      if (!aceInInt (ai, &acceptorSupport))     /* Number of supports of donor */
	continue ;
      lna = lnc = 30 ;
      beta1 = beta2 = 0 ;
      if (0)
	{
	  if (type == 2 || type >= 4) /* intron */
	    {
	      aceInStep(ai,'\t') ;
	      cp = aceInWord (ai) ;    /* donor of previous exon or SL, going into the SL*/
	      if (! cp)  continue ;
	      dictAdd (dict, cp+6, &beta1) ;
	      lna = strlen (cp+6) ;
	      
	      aceInStep(ai,'\t') ;
	      if (!aceInInt (ai, &donorSupport))     /* Number of supports of donor */
		continue ;
	      aceInStep(ai,'\t') ;
	      cp = aceInWord (ai) ;    /* acceptor of following exon */
	      if (! cp)  continue ;
	      dictAdd (dict, cp+9, &beta2) ;
	      lnc = strlen (cp+9) ;
	      
	      if (type >= 4)
		{
		  if (isDown) b2 = b1 ; /* so b2+1 is inside the exon */
		  else b1 = b2 ;
		} 
	      aceInStep(ai,'\t') ;
	      if (!aceInInt (ai, &acceptorSupport))     /* Number of supports of acceptor */
		continue ;
	    }
	  else if (type == 3) /* pA */
	    {
	      aceInStep(ai,'\t') ;
	      cp = aceInWord (ai) ;    /*p[AT]:AAAA the sequence of the polyA */
	      if (! cp)  continue ;
	      dictAdd (dict, cp+3, &beta2) ; 
	      lna = 0 ; lnc = strlen (cp+3) ;
	    }
	  if (! isDown)
	    {
	      int x = lna ; lna = lnc ; lnc = x ;
	      x = beta1 ; beta1 = beta2 ; beta2 = x ; 
	    }
	}

      dbli = bigArrayp (introns, nn++, DBLI) ;
      dbli->chrom = chrom ;
      dbli->isDown = isDown ;
      dbli->b1 = b1 ;
      dbli->b2 = b2 ;
      dbli->beta1 = beta1 ;
      dbli->beta2 = beta2 ;
      dbli->lna = lna ;
      dbli->lnc = lnc ;
      dbli->typea = 1 ;
      dbli->typeb = type ;
      dbli->typec = 1 ;
      dbli->donorSupport = donorSupport ;
      dbli->acceptorSupport = acceptorSupport ;
      dbli->support = donorSupport < acceptorSupport ?  donorSupport : acceptorSupport ;

      dbli->a1 = dbli->a2 = dbli->b1 ; /* default values */
      dbli->c1 = dbli->c2 = dbli->b2 ;
      if (type == 2)
	{
	  dbli->a1 = dbli->b1 - lna ; /* left intron foot */
	  dbli->a2 = dbli->b1 - 1 ; 
	  dbli->c1 = dbli->b2 + 1 ;
	  dbli->c2 = dbli->b2 + lnc ; /* right intron foot */
	}
     }
  bigArraySort (introns, dbliB1Order) ;
  sxNewDoubleIntronsGetDeUnoIntronsCompress (introns) ;

  if (bigArrayMax (introns) && sx->minimalIntronSupport > 1)
    {
      DBLI *up, *vp ;
      long int i, j, min = sx->minimalIntronSupport, iMax = bigArrayMax (introns) ;

      for (i = j = 0 , up = vp = bigArrp (introns, 0, DBLI) ; i < iMax ; up++, i++)
	{
	  if (up->support < min)
	    continue ;
	  if (j < i) *vp = *up ;
	  j++ ; vp++ ;
	}
      bigArrayMax (introns) = j ;
    }

  return introns ;
} /* sxNewDoubleIntronsGetIntrons */

/*************************************************************************************/
/* parse the deUno intron  file exported by clipalign as xxx.introns.gz 
   BIGINTRON       Forward 20      9049701 9049853 20      9198037 9198099 gt_ag   148183  1
*/
static BigArray sxNewDoubleIntronsGetDeUnoIntrons (SX *sx, DICT *dict, AC_HANDLE h)
{
  ACEIN ai = 0 ;
  DBLI *dbli ;
  BigArray introns = 0 ;
  long int nn = 0 ;
  int support, nn1, a1, a2, b1, b2, c1, c2, chrom, foot, ln ;
  int minA1 = sx->minX ;
  int maxA1 = sx->maxX ;
  BOOL isDown = TRUE ;
  BOOL stranded = sx->stranded ;
  const char *ccp ;
  int type ;
  
  introns = bigArrayHandleCreate (10000, DBLI, h) ;
  ai = 0 ;
  if (sx->sxxDeUnoIntronsFileName)
    ai = aceInCreate (sx->sxxDeUnoIntronsFileName, 0,  h) ;
  while (ai && aceInCard (ai))
    {
      type = 0 ;
      if (! (ccp = aceInWord (ai)))
	continue ;
      if (! strcmp (ccp, "INTRON"))
	type = 2 ;
      else if (! strcmp (ccp, "BIGINTRON"))
	type = 2 ;
      else if (! strcmp (ccp, "HUGEINTRON"))
	type = 2 ;
      else if (! strcmp (ccp, "NSINTRON"))
	type = 2 ;
      else if (! strcmp (ccp, "NSBIGINTRON"))
	type = 2 ;

      if (! type)
	continue ;

      aceInStep(ai,'\t') ;
      if (! (ccp = aceInWord (ai)))
	continue ;
      if (! strcmp (ccp, "Forward"))
	isDown = TRUE ;
      else if (! strcmp (ccp, "Reverse"))
	isDown = FALSE ;
      else
	continue ;

      aceInStep(ai,'\t') ;
      ccp = aceInWord (ai) ;    /* chromosome of first exon */
      if (! ccp)  continue ;
      if (sx->sxxChromosome && strcmp (sx->sxxChromosome, ccp))
	continue ;
      dictAdd (dict, ccp , &chrom) ;

      aceInStep(ai,'\t') ;
      if (!aceInInt (ai, &a1))     /* position of first matching base of upstream exon */
	continue ;
      aceInStep(ai,'\t') ;
      if (!aceInInt (ai, &a2))     /* position of last matching base */
	continue ;

      aceInStep(ai,'\t') ;
      ccp = aceInWord (ai) ;    /* chromosome of second exon */
      if (! ccp)  continue ;
      
      aceInStep(ai,'\t') ;
      if (!aceInInt (ai, &c1))     /* position of first matching base of upstream exon */
	continue ;
      aceInStep(ai,'\t') ;
      if (!aceInInt (ai, &c2))     /* position of last matching base */
	continue ;

     if (isDown)
	{
	  b1 = a2 + 1 ; b2 = c1 - 1 ;
	}
      else
	{ 
	  int x ;
	  b1 = a2 - 1 ; b2 = c1 + 1 ;
	  x = a1 ; a1 = c2 ; c2 = x ; x = a2 ; a2 = c1 ; c1 = x ; x = b2 ; b2 = b1 ; b1 = x ;
	}

      if ((minA1 && b1 < minA1) || (maxA1 && b1 > maxA1))
	continue ;
      if ((minA1 && b2 < minA1) || (maxA1 && b2 > maxA1))
	continue ;

      aceInStep(ai,'\t') ;
      ccp = aceInWord (ai) ;    /* intron foot gt_ag ... */
      if (! ccp)  
	ccp = "-" ;
      dictAdd (dict, ccp, &foot) ;
      if (! stranded && ! strcmp (ccp, "ct_ac"))
	isDown = ! isDown ;

      aceInStep(ai,'\t') ;
      aceInInt (ai, &ln) ;
      support = 1 ;
      aceInStep(ai,'\t') ;
      aceInInt (ai, &support) ;

      dbli = bigArrayp (introns, nn++, DBLI) ;
      dbli->chrom = chrom ;
      dbli->isDown = isDown ; /* stranded ? isDown : TRUE ; */
      dbli->a1 = a1 ;           /* left intron foot */
      dbli->a2 = a2 ;
      dbli->b1 = b1 ;
      dbli->b2 = b2 ;
      dbli->c1 = c1 ;              /* right intron foot */
      dbli->c2 = c2 ;
      dbli->beta1 = 0 ;
      dbli->beta2 = 0 ;
      dbli->lna = a2 - a1 + 1 ;
      dbli->lnc = c2 - c1 + 1 ;
      dbli->typea = 1 ;
      dbli->typeb = type ; 
      dbli->typec = 1 ;
      dbli->support = support;
      dbli->donorSupport = support ;
      dbli->acceptorSupport = support ;

      if (nn1 > 1000000)  /* 10M */
	{ nn1 = 0 ; nn = sxNewDoubleIntronsGetDeUnoIntronsCompress (introns) ; }
     }

  sxNewDoubleIntronsGetDeUnoIntronsCompress (introns) ;


  if (bigArrayMax (introns) && sx->minimalIntronSupport > 1)
    {
      DBLI *up, *vp ;
      long int i, j, min = sx->minimalIntronSupport ;

      for (i = j = 0 , up = vp = bigArrp (introns, 0, DBLI) ; i < bigArrayMax (introns) ; up++, i++)
	{
	  if (up->support < min)
	    continue ;
	  if (j < i) *vp = *up ;
	  j++ ; vp++ ;
	}
      bigArrayMax (introns) = j ;
    }

  if (! bigArrayMax (introns))
    bigArrayDestroy (introns) ;
  return introns ;
} /* sxNewDoubleIntronsGetDeUnoIntrons */

/*************************************************************************************/
/* search for intron-exon-intron segments */
static long int sxNewDoubleIntrons (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  long int ndbli = 0 ;
  DBLI *up = 0,  *dbli ;
  BigArray dblis = bigArrayHandleCreate (1000, DBLI, h) ;
  BigArray intronsB1, intronsB2, intronsDeUno ;
  DICT *dict = dictHandleCreate (10000, h) ;
  BOOL debug = FALSE ;

  fprintf (stderr, "// sxNewDoubleIntrons start: %s\n", timeShowNow ()) ;

  intronsDeUno = sxNewDoubleIntronsGetDeUnoIntrons (sx, dict, h) ;

  dblis = sxNewDoubleIntronsGetDeDuoIntrons (sx, dict, h) ;
  intronsB1 = bigArrayHandleCopy (dblis, h) ;
  intronsB2 = bigArrayHandleCopy (intronsB1, h) ;
  bigArraySort (intronsB1, dbliB1Order) ;
  bigArraySort (intronsB2, dbliB2Order) ;

  
  /* add the de uno introns */
  if (intronsDeUno)
    {
      long int i, iMax = bigArrayMax (intronsDeUno) ;
      ndbli = bigArrayMax (dblis) ;
      for (i = 0 ; i < iMax ; i++)
	{
	  dbli = bigArrayp (dblis, ndbli++, DBLI) ;
	  up = bigArrp (intronsDeUno, i, DBLI) ;
	  *dbli = *up ;
	}
      sxNewDoubleIntronsGetDeUnoIntronsCompress (dblis) ;
    }
  if (0) 
    {
      sxNewDoubleIntronsBridge (sx, dblis, intronsB1, intronsB2) ;
      /* merge */
      bigArraySort (dblis, dbliB1Order) ;
      sxNewDoubleIntronsDbliCompress (dblis) ;
    }

  sxNewDoubleIntronsExport (sx, dblis, dict, debug) ;
  ac_free (h) ;
  return ndbli ;
} /* sxNewDoubleIntrons */

/*************************************************************************************/
/*************************************************************************************/

typedef struct ixpStruct { BOOL reverse, isDown, isDownA ; int type, chrom, chromA, dan, ddn, bestdx, dx1, exon1, exon2, u1, u2, du ; char foot[6] ; } IXP ;

static int ixpFilterOrder (const void *a, const void *b)
{
  const IXP *up = (const IXP *) a, *vp = (const IXP *) b ;
  int n ;

  n = up->exon1 - vp->exon1 ; if (n) return n ;
  n = up->exon2 - vp->exon2 ; if (n) return n ;
  n = up->du - vp->du ; if (n) return n ;
  n = up->chrom - vp->chrom ; if (n) return n ;
 
  return 0 ;
} /* ixpFilterOrder  */

/**********************/

static int ixpExportOrder (const void *a, const void *b)
{
  const IXP *up = (const IXP *) a, *vp = (const IXP *) b ;
  int n ;

  n = up->chrom - vp->chrom ; if (n) return n ;
  n = up->u1 - vp->u1 ; if (n) return n ;
  n = up->u2 - vp->u2 ; if (n) return n ;

  return 0 ;
} /* ixpExportOrder  */

/**********************/

static Array sxParseGenomeMask (SX *sx, DICT *dict, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0, a1, a2, chrom ;
  const char *ccp ;
  H2C *up ;
  Array aa = 0 ;
  ACEIN ai = aceInCreate (sx->genomeMaskFileName, FALSE, h) ;

  aceInSpecial (ai, "\t\n\\\"") ;
  if (ai)
    {
      aa = arrayHandleCreate (1000, H2C, h0) ;
      while (aceInCard (ai))
	{
	  ccp = aceInWord (ai) ;
	  if (! ccp || *ccp == '#' || *ccp == '/')
	    continue ;
	  dictAdd (dict, ccp, &chrom) ;
	  if (aceInInt (ai, &a1) && aceInInt (ai, &a2))
	    {
	      up = arrayp (aa, nn++, H2C) ;
	      up->target = chrom ;
	      if (a1 < a2) { up->a1 = a1 ; up->a2 = a2 ; }
	      else  { up->a1 = a2 ; up->a2 = a1 ; }
	    }
	}
      if (nn) arraySort ( aa, H2Corder) ;
      else arrayDestroy (aa) ;
    }
  ac_free (h) ;
  return aa ;
}

/**********************/

static int sxNewIntrons (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  IXP *up, *vp ;
  Array exportIntrons = 0 ; 
  int nOk0 = 0, nOk1 = 0, nOk2 = 0, nOkRejected = 0, i, j, ia = 0, id = 0, ipA = 0, iSL = 0
    , ja, k1, k2, bestk, bestk1, bestk2, bestdx, bestdx1, bestdx2, a1, x1, x2,  dx1, dx2, score, chrom, exon, hook, maxA, maxD, maxSL ;
  const char *exon1, *exon2, *hook1, *hook2 ;
  char *cp, buf[1001] ;
  const char *ccp, *ccq ;
  int overhangLength = sx->overhangLength ? sx->overhangLength : 8 ;
  int seedLength = sx->seedLength > 0 ? sx->seedLength : 25 ;
  int minimalSupport = sx->minimalSupport > 0 ? sx->minimalSupport : 1 ;
  int minA1 = sx->minX ;
  int maxA1 = sx->maxX ;
  int entropy, mult, multT ;
  int intronMinLength = sx->intronMinLength > 0 ? sx->intronMinLength : 30 ;
  int intronMaxLength = sx->intronMaxLength > 0 ? sx->intronMaxLength : 10000 ;
  char polyABuf[overhangLength + 1] ;
  int polyA, polyT, SL ;
  BOOL isDonor, isAcceptor, isHit, isPa, isPt, isDown, stranded = sx->stranded, reverse ;
  ACEIN ai = 0 ;
  DA *dd = 0, *da = 0 ;
  Array polyAs = arrayHandleCreate (100000, DA, h) ;
  Array SLs = arrayHandleCreate (100000, DA, h) ;
  Array donors = arrayHandleCreate (100000, DA, h) ;
  Array acceptors = arrayHandleCreate (100000, DA, h) ;
  Array genomeMask = 0 ;
  DICT *dict = dictHandleCreate (10000, h) ;
  char aaaBuf[101] ;
  char f1[6], f2[6], foot[6] ;

  if (sx->transsplicing &&  minimalSupport  < sx->transsplicing)
    minimalSupport = sx->transsplicing ;

  memset (aaaBuf, 'A', 100) ;
  memset (polyABuf,'a', overhangLength) ;
  polyABuf [overhangLength] = 0 ;
  dictAdd (dict, polyABuf, &polyA) ;
  memset (polyABuf,'t', overhangLength) ;
  polyABuf [overhangLength] = 0 ;
  dictAdd (dict, polyABuf, &polyT) ;

  if (sx->genomeMaskFileName)
    genomeMask = sxParseGenomeMask (sx, dict, h) ;
  /* Parse the output of probealign -slide
   * Input should look like:
   * CHROMOSOME_I    2918821      probe  21      taggcttaggcttaggcttaggc caagcctaagcccaa DONOR   Forward
   */
  fprintf (stderr, "// sxNewIntrons start: %s\n", timeShowNow ()) ;
  if (1) ai = aceInCreate (sx->hitFile, FALSE,  h) ;

  while (ai && aceInCard (ai))
    {
      cp = aceInWord (ai) ;    /* probe id, count the multiplicity */
      if (! cp)  continue ;

      aceInStep(ai,'\t') ; aceInInt (ai, &score) ;     /* score, drop it */
      multT = 0 ;
      aceInStep(ai,'\t') ; aceInInt (ai, &mult) ;     /* multiplicity */

      aceInStep(ai,'\t') ; aceInInt (ai, &x1) ;     /* probe coord, drop it */
      aceInStep(ai,'\t') ; aceInInt (ai, &x2) ;     /* probe coord, drop it */


      /* type */
      aceInStep(ai,'\t') ; cp = aceInWord (ai) ;    /* HIT DONOR ACCEPTOR */
      if (! cp)  continue ;
      isHit = ! strcmp (cp, "HIT") ;

      isDonor = ! strcmp (cp, "DONOR") ;
      isAcceptor = ! strcmp (cp, "ACCEPTOR") ;
      
      isPa = ! strcmp (cp, "pA") ;
      isPt = ! strcmp (cp, "pT") ;

      if (! strncmp (cp, "SL", 2))
	dictAdd (dict, cp, &SL) ;
      else
	SL = 0 ;
      
      if (isHit)
	continue ;
      
      aceInStep(ai,'\t') ; cp = aceInWord (ai) ;    /* chromosome */
      if (! cp)  continue ;
      if (sx->sxxChromosome && strcmp (sx->sxxChromosome, cp))
	continue ;
      strncpy (buf, cp, 1000) ;
      
      aceInStep(ai,'\t') ; cp = aceInWord (ai) ;    /* strand */
      if (! cp)  continue ;
      
      /* since we drop the probe coord, we change the meaning */
      if (! strcmp (cp, "Forward"))
	isDown = ((x1 > x2  && isDonor) ||  (x1 < x2  && isAcceptor)) ? TRUE : FALSE ;
      else
	isDown = ((x1 > x2  && isDonor) ||  (x1 < x2  && isAcceptor)) ? FALSE : TRUE ;

      if (! strcmp (cp, "Forward"))
	isDown = TRUE ;
      else
	isDown = FALSE ;

      aceInStep(ai,'\t') ; if (!aceInInt (ai, &a1))     /* position of last matching base */
	continue ;
      
      if (isPt)
	isDown = !isDown ;

      if ((minA1 && a1 < minA1) || (maxA1 && a1 > maxA1))
	continue ;
      
      if (genomeMask)
	{
	  H2C *gm ;
	  int igm, igmMax = arrayMax(genomeMask) ;
	  BOOL bad = FALSE ;

	  for (gm = arrp (genomeMask, 0, H2C), igm = 0 ; igm < igmMax ; igm++, gm++)
	    {
	      if (chrom > gm->target || a1 > gm->a2) break ;
	      if (chrom == gm->target && a1 >= gm->a1 && a1 <= gm->a2) { bad = TRUE ; break ; }
	    }
	  if (bad)
	    continue ;
	}
      aceInStep(ai,'\t') ; cp = aceInWord (ai) ;    /* sequence matching the exon, oriented into the exon */
      if (! cp)  continue ;
      if (strlen (cp) < seedLength) continue ;
      cp[seedLength] = 0 ;             /* trim at say 15, so they look identical */
      if (cp) 
	dictAdd (dict, cp, &exon) ;
      
      aceInStep(ai,'\t') ; cp = aceInWord (ai) ;    /* sticky overhanging sequence, oriented into the intron, i.e. into the previous exon or the SL */
      if (! cp)  continue ;
      if (isPa || isPt)
	multT = strlen (cp) ;
      if (! SL)
	{
	  if (strlen (cp) < overhangLength) continue ;
	  cp[overhangLength] = 0 ;             /* trim at say 8, so they look identical */
	}
      if (cp) 
	dictAdd (dict, cp, &hook) ;
      entropy = hookEntropy (cp) ;
  
      dictAdd (dict, buf, &chrom) ;
      
      da = 0 ;
      if (isPa || isPt)
	{ da = arrayp (polyAs, ipA++, DA) ; da->isDown = isDown ; }
      else if (SL)
	{ da = arrayp (SLs, iSL++, DA) ; da->isDown = isDown ; }
      else if (0 && 2*entropy < overhangLength) /* entropy is verified in probealign -splice */
	{ nOkRejected++ ; continue ; }
      /* forget the strands in the input file to be able to cumulate stranded and non stranded data */
      else if ((isDonor && isDown) || (isAcceptor && ! isDown))
	{ da = arrayp (donors, id++, DA) ; da->isDown = TRUE ; } /* forget the strand */
      else if ((isDonor && ! isDown) || (isAcceptor && isDown))
	{ da = arrayp (acceptors, ia++, DA) ; da->isDown = TRUE ; } 

      if (!da)
	continue ;
      da->chrom = chrom ;
      da->a1 = a1 ;
      da->exon = exon ;
      da->hook = hook ; 
      da->SL = SL ;
      da->n = mult ; /* if fastc : use the multiplicity */
      da->nT = mult * multT ; /* if fastc : use the multiplicity */
    }

  fprintf (stderr, "// sxNewIntrons parsed %u donors %d acceptors %d SL %d polyA,  rejected %d overhangs with low entropy: %s\n"
	   , arrayMax (donors), arrayMax (acceptors), arrayMax (SLs), arrayMax (polyAs), nOkRejected, timeShowNow ()) ;

  /* sort */
  if (0) showDA ("Donors raw", dict, donors) ;
  arraySort (donors, DAorderHook) ;
  if (0) showDA ("Donors sorted", dict, donors) ;
  fprintf (stderr, "// sxNewIntrons sorted the donors %u donors %u acceptors : %s\n"
	   , arrayMax (donors), arrayMax (acceptors), timeShowNow ()) ;

  if (0) showDA ("Acceptors raw", dict, acceptors) ;
  arraySort (acceptors, DAorderHook) ;
  if (0) showDA ("Acceptors sorted", dict, acceptors) ;

  fprintf (stderr, "// sxNewIntrons sorted the acceptors %u donors %u acceptors : %s\n"
	   , arrayMax (donors), arrayMax (acceptors), timeShowNow ()) ;

  arraySort (SLs, DAorderHook) ;
  arraySort (polyAs, DAorderHook) ;

  fprintf (stderr, "// sxNewIntrons sorted %d Sls, %d polyA : %s\n"
	   , arrayMax (SLs), arrayMax (polyAs), timeShowNow ()) ;

  /* merge identical boundaries, and keep those with enough support */
  maxD = sxNewIntronsCompress (donors, minimalSupport) ;
  maxA = sxNewIntronsCompress (acceptors, minimalSupport) ;
  maxSL = sxNewIntronsCompress (SLs, minimalSupport) ;
  /* int maxpA = sxNewIntronsCompress (polyAs, minimalSupport) ; */

  arraySort (donors, DAorderN) ;
  arraySort (acceptors, DAorderN) ;
  arraySort (polyAs, DAorderN) ;
  arraySort (SLs, DAorderN) ;

  if (0) showDA ("Donors compressed", dict, donors) ;
  if (0) showDA ("Acceptors compressed", dict, acceptors) ;
  fprintf (stderr, "// sxNewIntrons consolidated %u donors %u acceptors %u SL %u pA : %s\n"
	   , arrayMax (donors), arrayMax (acceptors)
	   , arrayMax (SLs), arrayMax (polyAs)	   
	   , timeShowNow ()) ;

  {
    /* match the pairs and report the candidate introns */
    /* we do 3 pass 
     * pass 0: match donors to acceptors 
     * pass 1: match donors to donors to detect inversion and transloc
     * pass 1: match acceptors to acceptors to detect inversion and transloc
     */
    
    char *typeName[] = {"ZERO","INTRON","MICRODELETION","DELETION","MICRODUPLICATION","DUPLICATION","INVERSION","PALINDROME","TRANSLOCATION", "CIRCLE", 0} ;
  
    Array trueDonors = donors ;
    Array trueAcceptors = acceptors ;
    int pass, type, bestType, bestDu, iiStrands ; 
 
    arrayDestroy (exportIntrons) ;
    exportIntrons = arrayHandleCreate (10000, IXP, h) ;
    for (pass = 0 ; pass < (sx->transsplicing ? 4 : 1) ; pass++)
      {
	int ddChrom, oldDdChrom ;
	BOOL ddDown, oldDdDown ;
	int ddA1 ;

	switch (pass)
	  {
	  case 0: 
	  case 3:
	    donors = trueDonors ;
	    acceptors = trueAcceptors ;
	    break ;
	  case 1: 
	    donors = trueDonors ;
	    acceptors = trueDonors ;
	    break ;
	  case 2: 
	    donors = trueAcceptors ;
	    acceptors = trueAcceptors ;
	    break ;
	  }
	maxD = arrayMax (donors) ;
	maxA = arrayMax (acceptors) ;

	ddChrom = oldDdChrom = -1 ; ja = 0 ;
	ddDown = oldDdDown = FALSE ;
	for (id = ja = 0 ; id < maxD ; id++)
	  {
	    dd = arrp (donors, id, DA) ; 
	    if (dd->n < minimalSupport) continue ;

	    ddDown = dd->isDown ;
	    switch (pass)
	      {
	      case 1:
	      case 2:
		if (! ddDown) 
		  continue ;
		break ;
	      }

	    ddChrom = dd->chrom ;
	    ddA1 = dd->a1 ;

	    if (dd->chrom != oldDdChrom || ddDown != oldDdDown)
	      ja = 0 ;
	    oldDdChrom = ddChrom  ;
	    oldDdDown = ddDown ;

	    bestType = 999 ; bestDu = 0 ;

	    if (pass > 0 && ja <= id) ja = id + 1 ; /* same true-set. avoid double counting */
	    /* printf ("Pass = %d\n", pass) ; */
	    for (ia = ja, da = arrp (acceptors, ia, DA) ; ia < maxA ; da++, ia++)
	      {  
		if (da->n < minimalSupport) continue ;
		if (! sx->transsplicing)
		  {  /* pass == 0 */
		    if  (da->chrom > ddChrom || da->isDown > ddDown || da->a1 > ddA1 + intronMaxLength)
		      break ;
		    if  (da->chrom < ddChrom || da->isDown <  ddDown || da->a1 < ddA1 - 20)
		      { ja = ia ; continue ; } /* start next acceptor on same chrm same strand donor at this position ; */
		  }
		else if (pass < 3) /* same strand, same chrom */
		  {
		    if (da->chrom < ddChrom || da->isDown < ddDown)   /* impossible topology the begin of a read must hook to the end of a read */
		      { ja = ia ; continue ; } /* start next acceptor on correct strand */
		    if (da->chrom > ddChrom || da->isDown > ddDown)  /* there is no further compatible read-topology */
		      break ;
		  }
		else /* transposlocation inter chromosomal */
		{
		  if  (da->chrom == ddChrom )
		    continue ;
		}
		  
		type = 0 ; /* classify the single break points */
		if (da->chrom != ddChrom)
		  type = 8 ; /* TRANSLOCATION */
		else if (pass > 0)  /* same chrom same strand */
		  {
		      {
			if (da->a1 > ddA1 - 12 && da->a1 < ddA1 + 12)
			  type = 7 ;   /* PALINDROME */
			else
			  type = 6 ; /* INVERSION */
		      }
		  }
		else  if (da->a1 > ddA1)  /* [pass == 0 and same chrom same strand */
		  {
		    type = 3 ; /* DELETION */
		    if (da->a1 < ddA1 + intronMinLength)
		       type = 2 ; /* MICRODELETION */
		  }
		else if (da->a1 < ddA1)
		  {
		    type = 5 ; /* DUPLICATION */
		    if (da->a1 > ddA1 - intronMinLength)
		       type = 4 ; /* MICRODUPLICATION */
		  }

		if (!type || type > bestType)
		  continue ;
		switch (type)
		  {
		  case 0: 
		    continue ;
		  case 1:
		  case 2:
		  case 3:
		    if (bestType <= 3 && da->a1 - dd->a1 > bestDu)
		      continue ;
		    break ;
		  case 4:
		  case 5:
		    if (bestType <= 3 && dd->a1 - da->a1 > bestDu)
		      continue ;
		    break ;
		  }
		nOk0++ ;

		/* check if hook2 hooks into exon1 */
		exon1 = dictName (dict, dd->exon) ;
		hook2 = dictName (dict, da->hook) ;
		cp = strstr (exon1 + 2, hook2) ;
		if (!cp)
		  continue ;
		dx1 = cp - exon1 ;
		if (dx1 < 2)
		  continue ;
		
		/* check if hook1 hooks into exon2 */
		hook1 = dictName (dict, dd->hook) ;
		exon2 = dictName (dict, da->exon) ;
		cp = strstr (exon2 + dx1, hook1) ;
		if (!cp)
		  continue ;
		dx2 = cp - exon2 ;

		if (dx1 != dx2)
		  continue ;
		/* verify there is room for the gt-ag feet */
		if (dx1 + dx2 < 4)
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
		
		/* success, now we adjust the boundary and hope to slide to a gt-ag  */
		nOk1++ ;
		/*		dx = dx1 + dx2 - 4 ; */
		bestk = bestk1 = bestk2 = -1 ;
		bestdx = bestdx1 = bestdx2 = 0 ;
		reverse = FALSE ;
		foot[0] = '-' ; foot[1] = 0 ;
		f1[2] = '_' ; f1[5] = 0 ;
		f2[2] = '_' ; f2[5] = 0 ;
                if (type == 3 && sx->RNA_seq)
		  {
		    for (i = 0 ; i <= dx1 - 2 ; i++)
		      {
			k1 = k2 = 0 ; 
			ccp = exon1 + i  ; ccq = ccp+1 ;
			f2[3] = *ccp ; f2[4] = *ccq ;
			f1[1] = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)*ccp])] ;
			f1[0] = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)*ccq])] ;
			if (
			    (*ccp == 'a' && *ccq == 'c')   /*    d1 == "gt"  */
			    )
			  k1+=2 ;
			if (
			    (*ccp == 'g' && *ccq == 'c')     /* || d1 == "gc"  */
			    )
			  k1++ ;
			else if (!stranded)
			  {
			    if (*ccp == 'a' && *ccq == 'g')   /*    d1 == "ct"  */
			      k2+=2 ;  /*    d1 == "ag"  */
			  }
			
			ccp = exon2 + dx1 - i - 2 ; ccq = ccp+1 ;
			f1[3] = *ccp ; f1[4] = *ccq ;
			f2[1] = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)*ccp])] ;
			f2[0] = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)*ccq])] ;
			if (*ccp == 'a' && *ccq == 'g')   /*    d1 == "ag"  */
			  k1+=2 ;
			else if (!stranded)
			  {
			    if (
				(*ccp == 'a' && *ccq == 'c')     /*    d1 == "ac"  */
				)
			      k2+=2 ;
			    if (
				(*ccp == 'g' && *ccq == 'c')     /* || d1 == "gc"  */
				)
			      k2++ ;
			  }
			if ( k1 > bestk)
			  { bestk = k1 ; bestdx = i ; reverse = FALSE ; strcpy (foot, f1) ; }
			if ( k2 > bestk)
			  { bestk = k2 ; bestdx = i ; reverse = TRUE ; strcpy (foot, f2) ; }
			if ( k1 >= bestk)
			  { bestk1 = k1 ; bestdx1 = i ; }
			if ( k2 >= bestk)
			  { bestk2 = k2 ; bestdx2 = i ; }
		      }
		    if (bestk < 4 && ! (sx->non_gt_ag && dd->n >= sx->non_gt_ag && da->n >= sx->non_gt_ag ))
		      continue ;
		    if (bestk == 4) type = 1 ;
		    if (bestk1 < bestk) bestk1 = 0 ;
		    if (bestk2 < bestk) bestk2 = 0 ;
		    if (bestk1 && bestk2 && bestk < 4) bestk2 = 0 ; /* no ambiguous sliding non gt_ag */
		  }
		/* acceptable intron boundaries */
		for (iiStrands = 0 ; iiStrands < 2 ; iiStrands++)
		  {
		    DA *mydd = dd, *myda = da ;
		    
		    if (iiStrands == 0)
		      {
			if (bestk > bestk1) continue ; 
			reverse = FALSE ; bestdx = bestdx1 ;
		      }
		    else
		      {
			if (bestk > bestk2 || bestk1 == bestk2) continue ; 
			reverse = TRUE ; bestdx = bestdx2 ;
		      }
		    if (dd->chrom > da->chrom || (pass > 0 && dd->a1 > da->a1))
		      { mydd = da ; myda = dd ; }
		    up = arrayp (exportIntrons, nOk2++, IXP) ;
		    up->type = type ;
		    up->isDown = (pass == 2 ? !mydd->isDown : mydd->isDown) ;
		    up->isDownA = (pass == 1 ? !myda->isDown : myda->isDown) ;
		    up->reverse = reverse ;
		    up->chrom = mydd->chrom ;
		    up->chromA = myda->chrom ;
		    up->ddn = mydd->n ;
		    up->dan = myda->n ;
		    up->bestdx = bestdx ;
		    up->dx1 = dx1 ;
		    strcpy (up->foot, foot) ;
		    
		    if (! reverse)
		      {
			up->exon1 = mydd->exon ;
			up->exon2 = myda->exon ;
		      }
		    else
		      {
			up->exon1 = myda->exon ;
			up->exon2 = mydd->exon ;
		      }
		    
		    up->u1 = mydd->a1 + (mydd->isDown ? - bestdx + 1 : bestdx - 1) ;
		    up->u2 = myda->a1 + (myda->isDown ?  dx1 - bestdx - 3 : - dx1 + bestdx + 3) ;
		    if (myda->chrom == mydd->chrom)
		      {
			if (pass == 0 && mydd->isDown == myda->isDown)
			  up->du = (mydd->isDown ? up->u2 - up->u1 + 1 : up->u1 - up->u2 + 1) ;
			else
			  up->du = (up->u1 < up->u2 ?  up->u2 - up->u1 : up->u1 - up->u2 ) ; 
		      }
		    else
		      up->du = 0 ;
		    
		    bestType = type ;
		    bestDu = up->du ;
		  }
		if (pass == 0 && ! sx->transsplicing)
		  break  ; /* do not try another acceptor further away */
	      }
	  }
      }
    /* keep only the shortest intron with given feet, 
     * this avoids contruction mosaics in repeated genes or in
     * a case like 3D156 where exons 2 and 5 (length 202bp) are identical
     */
    arraySort (exportIntrons, ixpFilterOrder) ;
    for (i = 0 ; i < nOk2 - 1 ; i++)
      {
	up = arrp (exportIntrons, i, IXP) ;
	if (! up->chrom) continue ;
	for (j = i+1, vp = up + 1 ; up->chrom && j <  nOk2 && vp->exon1 == up->exon1 && vp->exon2 == up->exon2 ; vp++, j++)
	  {
	    if (vp->chrom ==  up->chrom && ((vp->du > 2*up->du) || (vp->du > up->du + 10 && up->u1 >= vp->u1 && up->u2 <= vp->u2)))
	      vp->chrom = 0 ;
	    if (vp->chrom  ==  up->chrom && ((vp->du > 2*up->du) || (up->du > vp->du + 10 && vp->u1 >= up->u1 && vp->u2 <= up->u2)))
	      up->chrom = 0 ;
	  }
      }
    
    arraySort (exportIntrons, ixpExportOrder) ;
    arrayCompress (exportIntrons) ;
    /* adjust the coordinates */
    if (1)
      for (i = 0, vp = 0 ; i < arrayMax (exportIntrons) ; i++)
	{
	  up = arrp (exportIntrons, i, IXP) ;
	  if (! up->chrom) continue ; 
	  if (up->u1 < up->u2)
	    up->du = up->u2 - up->u1 + 1 ;
	  else
	    up->du = up->u1 - up->u2 + 1 ;
	  switch (up->type)
	    {
	    case 2: /* MICRODELETION */
	    case 3: /* DELETION */

	      if (up->u1 < up->u2)
		{ up->u1++ ; up->u2-- ; }
	      else
		{ up->u1-- ; up->u2++ ; }
	      break ;
	    case 4: /* MICRODUPLICATION */
	    case 5: /* DUPLICATION */
	      if (up->u1 < up->u2)
		{ 
		  up->u1++ ; up->u2-- ; 
		  if (up->u1 == 1) up->type = 9 ; /* CIRCLE */
		}
	      else
		{ 
		  int u ;
		  up->u1-- ; up->u2++ ; 
		  if (up->u2 == 1) up->type = 9 ; /* CIRCLE */
		  u = up->u1 ; up->u1 = up->u2 ; up->u2 = u ;
		}
	      break ;
	    case 6: /* INVERSION */
	    case 7: /* PALINDROME */
	      if (up->u1 < up->u2)
		{ 
		  up->u1-- ; up->u2++ ; 
		  if (up->u2 == 1) up->type = 9 ; /* CIRCLE */
		}
	      else
		{ 
		  up->u1++ ; up->u2-- ;
		  if (up->u1 == 1) up->type = 9 ; /* CIRCLE */
		}
	      break ;
	    case 8: /* TRANSLOCATION */
	      break ;
	    }
	}
    /* export */
    for (i = 0, vp = 0 ; i < arrayMax (exportIntrons) ; i++)
      {
	up = arrp (exportIntrons, i, IXP) ;
	if (! up->chrom) continue ;
	
	if (vp 
	    && up->type == vp->type
	    && up->isDown == vp->isDown 
	    && up->isDownA == vp->isDownA 
	    && up->chrom == vp->chrom
	    && up->chromA == vp->chromA
	    && up->u1 == vp->u1
	    && up->u2 == vp->u2
	    && up->exon1 == vp->exon1
	    && up->exon2 == vp->exon2
	    && up->du == vp->du
	    && up->dx1 == vp->dx1
	    && up->bestdx == vp->bestdx
	    && up->ddn == vp->ddn
	    && up->dan == vp->dan
	    )
	  continue ;
	
	exon1 = dictName (dict, up->exon1) ;
	exon2 = dictName (dict, up->exon2) ;
	
	if (0) freeOutf ("exon1=%d exon2=%d\n", up->exon1,up->exon2) ;
	
	if (pass == 0 && up->type == 0)
	  {
	    if (! up->reverse)
	      freeOutf ("INTRON\t%s\t%s\t%09d\t%09d\tDONOR:%s\t%d\tACCEPTOR:%s\t%d\t%s\tdx1=%d  bestdx=%d  ln=%d\n"
			, up->isDown ? "Forward" : "Reverse"
			, dictName (dict, up->chrom)	
			, up->u1
			, up->u2
			, exon1 + up->bestdx + 2
			, up->ddn
			, exon2 + up->dx1 - up->bestdx  
			, up->dan
			, up->foot
			, up->dx1
			, up->bestdx
			, up->du
			) ;
	    else
	      freeOutf ("INTRON\t%s\t%s\t%09d\t%09d\tDONOR:%s\t%d\tACCEPTOR:%s\t%d\t%s\tdx1=%d  bestdx=%d  ln=%d \n"
			, up->isDown ? "Reverse" :  "Forward"
			, dictName (dict, up->chrom)
			, up->u2
			, up->u1
			, exon1 + up->dx1 - up->bestdx 
			, up->dan
			, exon2 + up->bestdx + 2
			, up->ddn
			, up->foot
			, up->dx1
			, up->bestdx
			, up->du
			) ;
	  }
	else
	  {
	    if (! (up->type == 1 && up->reverse))
	      freeOutf ("%s\t%s\t%s\t%09d\t%s\t%s\t%09d\t%d\t%s\tSUPPORT\t%s\t%d\t%d\t%d\t%d\tSequence\t%s\t%s\n"
			, typeName[up->type]
			, up->isDown ? "Forward" : "Reverse"
			, dictName (dict, up->chrom)	
			, up->u1
			, dictName (dict, up->chromA)	
			, up->isDownA ? "Forward" : "Reverse"
			, up->u2
			, up->du
			, up->foot
			, sx->run ? sx->run : "-"
			, up->ddn
			, up->dan
			, 0
			, 0
			, exon1 + up->bestdx + 2
			, exon2 + up->dx1 - up->bestdx  
			) ;
	    else
	      freeOutf ("%s\t%s\t%s\t%09d\t%s\t%s\t%09d\t%d\t%s\tSUPPORT\t%s\t%d\t%d\t%d\t%d\tSequence\t%s\t%s\n"
			, typeName[up->type]
			, ! up->isDown ? "Forward" : "Reverse"
			, dictName (dict, up->chromA)	
			, up->u2
			, dictName (dict, up->chrom)	
			, ! up->isDownA ? "Forward" : "Reverse"
			, up->u1
			, up->du
			, up->foot
			, sx->run ? sx->run : "-"
			, up->ddn
			, up->dan
			, 0
			, 0
			, exon1 + up->dx1 - up->bestdx 
			, exon2 + up->bestdx + 2
			) ;
	      }
	vp = up ;
      }
  }

  fprintf (stderr, "// sxNewIntrons Verified %d introns, tested %d, attempted %d, minimalSupport %d start: %s\n"
	   , nOk2, nOk1, nOk0, minimalSupport, timeShowNow ()) ;


  /* verify and report the candidate SLs */
  {
    int ns, ln, dx ;
    char *slu[] = { "ggtttaattacccaagtttgag",     /* SL1 */
		    "ggttttaacccagttactcaag" ,    /* SL2 */
		    "ggttttaacccagttaaccaag" ,    /* SL3 */
		   "ggttttaacccagtttaaccaag" ,    /* SL4 */
		     "ggttttaacccagttaccaag" ,    /* SL5 */
		    "ggtttaaaacccagttaccaag" ,    /* SL6 */
		    "ggttttaacccagttaattgag" ,    /* SL7 */
		    "ggtttttacccagttaaccaag" ,    /* SL8 */
		    "ggtttatacccagttaaccaag" ,    /* SL9 */
		   "ggttttaacccaagttaaccaag" ,    /* SL10 */
		     "ggttttaaccagttaactaag" ,    /* SL11 */
		    "ggttttaacccatataaccaag" ,    /* SL12 */
		    /* do not exist 
		     * "gtttttaacccagttactcaag" ,    SL13 
		     * "ggtttttaacccagttactcaag" ,     SL14 
		     */
		    0 } ;

    for (id = 0 ; id < maxSL ; id++)
    {
      char buf[40] ;
      dd = arrp (SLs, id, DA) ; 
      if (! dd->n) continue ;
      
      exon1 = dictName (dict, dd->exon) ;
      hook1 = dictName (dict, dd->hook) ;
      
      ns = (dictName (dict,dd->SL))[2] - '0' ; /* identify the SL */
      if ((dictName (dict,dd->SL))[3])
	ns = 10 * ns + (dictName (dict,dd->SL))[3] - '0' ;
      /* count the sliding bases */
      ccp = slu[ns-1] ; 
      ln = strlen (ccp) ;
      for (i = 12, dx = 0 ; i >= 2 ; i--)
	if (!strncmp (exon1+2, ccp + ln - i, i))
	  { dx = i ; break ; }
      /* at least ag should slide to define a correct acceptor */
      if (dx < 2)
	continue ;
      for (i = 0 ; i < dx ; i++)
	{
	  switch (*(ccp + ln - 1 - i))
	    {
	    case 'a': buf[i] = 't' ; break ;
	    case 't': buf[i] = 'a' ; break ;
	    case 'g': buf[i] = 'c' ; break ;
	    case 'c': buf[i] = 'g' ; break ;
	    }
	}
      strcpy (buf + dx, hook1) ;
      /* else, ok, export the correct coordinate */
      freeOutf ("SL%d\t%s\t%s\t%09d\t%09d\tDONOR:%s\t%d\tACCEPTOR:%s\t%d\tdx1=%d\n"
		, ns
		, dd->isDown ? "Forward" : "Reverse"
		, dictName (dict, dd->chrom)
		, dd->a1 + (dd->isDown ? + dx - 1 : - dx + 1)
		, dd->a1 + (dd->isDown ? + dx - 0 : - dx + 0)
		, buf
		, dd->n
		, exon1+2+dx
		, dd->n
		, dx
		) ;
    }
  }

  /* verify and report the candidate pA */
  {
    for (id = 0 ; id < arrayMax (polyAs) ; id++)
      {
	dd = arrp (polyAs, id, DA) ; 
	if (! dd->n) continue ;
	
	/* count the A */
	if (dd->n)
	  {
	    if (dd->nT/dd->n < 100)
	      aaaBuf[dd->nT/dd->n] = 0 ;
	    freeOutf ("pA\t%s\t%s\t%09d\t%09d\tpA:%s\t%d\n"
		      , dd->isDown ? "Forward" : "Reverse"
		      , dictName (dict, dd->chrom)
		      , dd->a1
		      , dd->a1 + (dd->isDown ?  1 : - 1)
		      , aaaBuf
		      , dd->n
		      ) ;
	    aaaBuf[dd->nT/dd->n] = 'A' ;
	  }
      }
  }


  ac_free (h) ;
  return nOk2 ;
} /* sxNewIntrons */
/*
INTRON  Forward 3J90    000044099       000095932       DONOR:ttcctcca  8       ACCEPTOR:acaggcaa       23      dx=3  bestdx=1
INTRON  Forward 3J90    000047794       000054722       DONOR:ctgtaatt  4       ACCEPTOR:tgcgattt       3       dx=1  bestdx=0
INTRON  Forward 3J90    000047794       000097378       DONOR:ctgtaatt  4       ACCEPTOR:tagatttg       22      dx=0  bestdx=0
INTRON  Forward 3J90   * 000095566       000095648       DONOR:acatcacg  34      ACCEPTOR:agcttaat       24      dx=0  bestdx=0
INTRON  Forward 3J90   * 000095824       000095932       DONOR:aacctcca  28      ACCEPTOR:acaggcaa       23      dx=3  bestdx=1
INTRON  Forward 3J90   * 000096060       000096310       DONOR:tcacctct  18      ACCEPTOR:acaggttg       17      dx=4  bestdx=2
INTRON  Forward 3J90   * 000096060       000096753       DONOR:cacctctt  8       ACCEPTOR:acaggccg       3       dx=3  bestdx=1
INTRON  Forward 3J90   * 000096060       000097097       DONOR:acctcttc  12      ACCEPTOR:gcagcacg       8       dx=2  bestdx=0
INTRON  Forward 3J90   * 000096464       000097379       DONOR:acccaaga  16      ACCEPTOR:tagatttg       11      dx=1  bestdx=0
INTRON  Forward 3J90   * 000097218       000097379       DONOR:accgtaat  18      ACCEPTOR:tagatttg       22      dx=1  bestdx=0
INTRON  Forward 3J90   * 000097588       000097643       DONOR:acccttca  23      ACCEPTOR:cagattcc       19      dx=1  bestdx=0
INTRON  Forward 3J90   -- 000098897       000099867       DONOR:tacctgtc  13      ACCEPTOR:tccaggca       14      dx=4  bestdx=1
INTRON  Forward 3J90    000098898       000099960       DONOR:tacctgtc  13      ACCEPTOR:aagcaaca       3       dx=1  bestdx=0
INTRON  Forward 3J90   *  000099141       000100120       DONOR:accaattg  3       ACCEPTOR:cagagtgt       3       dx=1  bestdx=0
INTRON  Forward 3J90    000100007       000100120       DONOR:acctgatg  17      ACCEPTOR:attcagag       26      dx=4  bestdx=0
INTRON  Forward 3J90   *  000100339       000100504       DONOR:acgagcat  41      ACCEPTOR:agagctct       37      dx=0  bestdx=0
INTRON  Forward 3J90   *  000100675       000100895       DONOR:aaccaatg  33      ACCEPTOR:cagggtcg       28      dx=2  bestdx=1


INTRON  Forward 3J90    000095566       000095648       DONOR:atcacgaaactga     34      ACCEPTOR:cttaatcccgatg  24      dx=0  bestdx=0
INTRON  Forward 3J90    000095824       000095932       DONOR:ctccaaatcaat      28      ACCEPTOR:gcaacattcga    23      dx=3  bestdx=1
INTRON  Forward 3J90    000096060       000097097       DONOR:ctcttcgaagaaa     12      ACCEPTOR:cacgaagacat    8       dx=2  bestdx=0
INTRON  Forward 3J90    000096060       000096753       DONOR:ctcttcgaagaa      8       ACCEPTOR:gccgggactgg    3       dx=3  bestdx=1
INTRON  Forward 3J90    000096060       000096310       DONOR:ctcttcgaaga       18      ACCEPTOR:gttgatcatga    17      dx=4  bestdx=2
INTRON  Forward 3J90    000096464       000097379       DONOR:ccaagacgtagca     16      ACCEPTOR:atttgttgctgg   11      dx=1  bestdx=0
INTRON  Forward 3J90    000097218       000097379       DONOR:cgtaattgcagtt     18      ACCEPTOR:atttgttgctgg   22      dx=1  bestdx=0
INTRON  Forward 3J90    000097588       000097643       DONOR:ccttcacaaattt     23      ACCEPTOR:attccgcgcaac   19      dx=1  bestdx=0
INTRON  Forward 3J90    000098897       000099867       DONOR:ctgtccaggata      13      ACCEPTOR:gcaacaatcg     14      dx=4  bestdx=1
INTRON  Forward 3J90    000098897       000099959       DONOR:ctgtccaggata      13      ACCEPTOR:gcaacaatcgaga  3       dx=1  bestdx=1
INTRON  Forward 3J90    000099141       000100120       DONOR:caattgtctatgg     3       ACCEPTOR:agtgttcgtgct   3       dx=1  bestdx=0
INTRON  Forward 3J90    000100007       000100120       DONOR:ctgatgttgattt     17      ACCEPTOR:agtgttcgt      26      dx=4  bestdx=0
INTRON  Forward 3J90    000100339       000100504       DONOR:gagcatgagccaa     41      ACCEPTOR:agctctctgatgt  37      dx=0  bestdx=0
INTRON  Forward 3J90    000100675       000100895       DONOR:caatgaagacga      33      ACCEPTOR:ggtcgagttcca   28      dx=2  bestdx=1
/
*/
/*************************************************************************************/
/*************************************************************************************/
/* open and scans the overhangs file fNam, ventilate it out into an array of 'outs' files */
static int sxDoVentilate (SX *sx, const char *fNam, DICT *dict, Array outs, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai ;
  ACEOUT ao ;
  const char *ccp ;
  char buf[4000], outNam[1000], type[32], chrom[128], *cp ;
  int i, ok = 0, outId ;

  ai = aceInCreate (fNam, FALSE,  h) ;
  if (! ai)
    messcrash ("Cannot open the overhangs file %s", fNam) ;
  aceInSpecial (ai, "\n") ;
  while (ai && aceInCard (ai))
    {
      ccp = aceInPos (ai) ;  /* full line, that we will ventilate in sub directories */
      if (*ccp == '#')
	continue ;
      strncpy (buf, ccp, 3999) ;
      
      /* get the type: OR, pA, SL in column 6 */
      for (i = 0 ; i < 6 ; i++)
	{ aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; }
      if (! ccp || strlen (ccp) > 31) continue ;  
      /* merge the donor/acceptor cases */
      if (! strcmp (ccp, "DONOR") || ! strcmp (ccp, "ACCEPTOR"))
	ccp = "OR" ;
      strcpy (type, ccp) ;
      
      /* get the chromosome from column 7 */
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; 
      if (! ccp || strlen (ccp) > 127) continue ;
      if (sx->maxChromNameLength &&  strlen (ccp) > sx->maxChromNameLength) continue ;
      if (0 && ! strncmp (ccp, "CHROMOSOME_", 11))
	ccp += 11 ;
      if (! *ccp) continue ;
      strcpy (chrom, ccp) ;
      if (1) /* discard the floating contigs */
	{
	  for (cp = chrom ; *cp ; cp++)
	    if (*cp == '|') { break ;}
	  if (*cp == '|') continue ;
	}
      else
	{
	  for (cp = chrom ; *cp ; cp++)
	    if (*cp == '|') { *cp = 0 ; break ;}
	}

      /* compose the out file name */
      sprintf (outNam, "tmp/%s/%s/%s/%s.clean", type, sx->run, chrom, sx->run) ;
      dictAdd (dict, outNam, &outId) ;

      /* create or recover the ACEOUT handle */
      ao = array (outs, outId, ACEOUT) ;
      if (! ao) /* open the file just once */
	{
	  int dummy ;
	  for (cp = outNam ; *cp ; cp++)
	    if (*cp == '/')
	      {
		*cp = 0 ;
		dummy = system (messprintf ("mkdir -p %s", outNam)) ;
		*cp = '/' ;
		if (dummy) {} ; /* for compiler happiness */
	      }
	  ao = aceOutCreate (outNam, 0, TRUE, h0) ; /*  aceOutCreateToPipe (hprintf (h, "gzip  > %s ", outNam), h0) ;*/
	  array (outs, outId, ACEOUT) = ao ; 
	}

      /* export */
      ok++ ;
      aceOutf (ao, "%s", buf) ; /* copy the whole line */
      aceOutf (ao, "\n") ; /* copy the whole line */
    }
  
  ac_free (h) ;
  return ok ;
} /* sxDoVentilate */

/*************************************************************************************/

static int sxVentilate (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  int nFiles = 0, nLines = 0, nOuts = 0 ;
  ACEIN fList ;
  const char *fNam ;
  Array outs = arrayHandleCreate (1000, ACEOUT, h) ;
  DICT *dict = dictHandleCreate (1000, h) ;

  if (! sx->run)
    messcrash ("-ventilate requires -run run_name, sorry") ;
  fList = aceInCreateFromPipe (messprintf ("ls tmp/PHITS_genome/%s/*.overhangs.gz", sx->run), "r", 0, h) ;
  aceInSpecial (fList, "\t\n") ;

  while (fList && aceInCard (fList))
    {
      while ((fNam = aceInWord (fList)))
	{
	  nFiles++ ;
	  nLines += sxDoVentilate (sx, fNam, dict, outs, h) ;
	}
    }
  nOuts =  dictMax (dict) ;
  ac_free (h) ; /* will close all aceout handles */

  fprintf (stderr, "sxVentilate run %s, exported %d lines from %d files in %d out files\n"
	   , sx->run, nLines, nFiles, nOuts) ;
  return nLines ;
} /* sxVentilate */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/* given a global wiggle, a set of known and candidate introns
 * scan the genome in back and forth and create segments and compute their median
 * a segment stops at an intron boundary and if the local median wiggle drops below a lowwater mark
 * or the local median drops to 1/8 of its previous value
 *
 * the report segments over a high-water mark as exons (or new exons)
 * and report candidate introns bordered by exons as new exons-exons boundaries
 */

/*************************************************************************************/
typedef struct sxxStruct { int a1, a2, b1, b2, m, n, type, eClass, donor, acc, aLevel, nDonorSupport, iLevel, nAcceptorSupport, bLevel ; BOOL isDown ; } SXX ;
typedef struct sxwStruct { int a1, m, n ; } SXW ;

static void showSXX (Array aa)
{
  int i ;
  SXX *s ;

  for (i = 0 ; aa && i < arrayMax (aa) ; i++)
    {
      s = arrp (aa, i, SXX) ;
      fprintf (stderr, "%d: %d %d %d %d m=%d n=%d type=%d e=%d nds=%d nas=%d\n"
	       , i
	       , s->a1, s->a2, s->b1, s->b2
	       , s->m, s->n, s->type, s->eClass
	       , s->nDonorSupport, s->nAcceptorSupport
	       ) ;
    }
  if (0) showSXX (0) ; /* for compiler happiness */
} /*showSXX */



static int sxxOrder (const void *a, const void *b)
{
  const SXX *x = (const SXX *) a, *y = (const SXX *) b ;
  if (x->a1 < y->a1) return -1 ;
  if (x->a1 > y->a1) return 1 ;
  if (x->a2 < y->a2) return -1 ;
  if (x->a2 > y->a2) return 1 ;
  if (x->b1 < y->b1) return -1 ;
  if (x->b1 > y->b1) return 1 ;
  if (x->b2 < y->b2) return -1 ;
  if (x->b2 > y->b2) return 1 ;
  if (x->eClass < y->eClass) return -1 ;
  if (x->eClass > y->eClass) return 1 ;
  if (x->n < y->n) return -1 ;
  if (x->n > y->n) return 1 ;

  return 0 ;
} /* sxxOrder */

static int sxxPaOrder (const void *a, const void *b)
{
  const SXX *x = (const SXX *) a, *y = (const SXX *) b ;
  if (x->b1 < y->b1) return -1 ;
  if (x->b1 > y->b1) return 1 ;
  if (x->a1 < y->a1) return -1 ;
  if (x->a1 > y->a1) return 1 ;
  if (x->b2 < y->b2) return -1 ;
  if (x->b2 > y->b2) return 1 ;
  if (x->a2 < y->a2) return -1 ;
  if (x->a2 > y->a2) return 1 ;
  if (x->eClass < y->eClass) return -1 ;
  if (x->eClass > y->eClass) return 1 ;

  return 0 ;
} /* sxxPaOrder */

/*************************************************************************************/

static Array sxNewExonsGetWiggle (SX *sx, int *stepp, const char *fileName, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  Array aa = 0, wiggle = arrayHandleCreate (100000, SXW, h0) ;
  SXW *sxw ;
  WIGGLEPOINT *wp ;
  int ii ;
  int minX, maxX ;
  BOOL noFilter = TRUE ;
  char *type = "BF" ;

  minX = sx->minX ;
  maxX = sx->maxX ;
  if (minX || maxX) noFilter = FALSE ;

  if (fileName)
    {
      char *cp, *cq, *buf ;
      cp = buf = strnew (fileName, h) ;
      
      if (! strstr (fileName, ".BF.") && ! strstr (fileName, ".bf.") &&
	  (strstr (fileName, ".BV.") || ! strstr (fileName, ".bv."))
	  )
	type = "BV" ;
    
      /* parse a collection of sponge files */
      *stepp = 0 ; /* start agnostic , but the first file imposes the step of the other ones */
      while (cp)
	{
	  cq = strstr (cp, ",") ;
	  if (cq)
	    *cq++ = 0 ;
	  aa = sxGetWiggleZone (aa, cp, type, stepp, 0, 0, 0, h) ;
	  cp = cq ;
	}
    }
  if (aa && arrayMax (aa))
    {
      if (noFilter)
	for (ii = 0, wp = arrp (aa, 0, WIGGLEPOINT) ; ii < arrayMax (aa) ; ii++, wp++) 
	  {
	    sxw = arrayp (wiggle, ii, SXW) ;
	    sxw->a1 = wp->x ;
	    sxw->n = wp->y ;
	  }
      else
	for (ii = 0, wp = arrp (aa, 0, WIGGLEPOINT) ; ii < arrayMax (aa) ; ii++, wp++) 
	  if (wp->x > minX && wp->x < maxX)
	    {
	      sxw = arrayp (wiggle, ii, SXW) ;
	      sxw->a1 = wp->x ;
	      sxw->n = wp->y ;
	    }
    }

  ac_free (h) ;
  return wiggle ;
} /* sxNewExonsGetWiggle */

/*************************************************************************************/

static void sxIntronsAddCounts (Array introns)
{
  SXX *sx, *sy ;
  int ii, jj, nD, nA = 0, nD1, nA1 = 0 ;
    
  if (0)
    {
      arrayCompress (introns) ;
      return ;
    }

  for (ii = 0, sx = arrp (introns, 0, SXX) ; ii < arrayMax (introns) ; sx++, ii++)
    {
      if (sx->m == -1)
	continue ;
      nD = sx->nDonorSupport ; sx->nDonorSupport = 0 ;
      if (sx->eClass != 3) { nA = sx->nAcceptorSupport ; sx->nAcceptorSupport  = 0 ; }
      /* add up the counts */
      for (jj = ii + 1, sy = sx + 1 ;  jj < arrayMax (introns) ; sy++, jj++)
	{  
	  if (sx->m == -1)
	    continue ;
	  nD1 = sy->nDonorSupport ; sy->nDonorSupport = 0 ;
	  if (sy->eClass != 3) { nA1 = sy->nAcceptorSupport ; sy->nAcceptorSupport  = 0 ; }
	  if (! memcmp (sx, sy, sizeof (SXX)))
	    { nD += nD1 ; nA += nA1 ; sy->m = -1 ; }
	  else /* restore */
	    { sy->nDonorSupport = nD1 ; if (sy->eClass != 3) sy->nAcceptorSupport = nA1 ; break ; }
	}
       sx->nDonorSupport = nD ; if (sx->eClass != 3) sx->nAcceptorSupport = nA ; 
    }
  /* keep happy few */
  for (ii = jj = 0, sx = sy = arrp (introns, 0, SXX) ; ii < arrayMax (introns) ; sx++, ii++)
    {
      if (sx->m == -1)
	continue ;
      if (jj < ii) *sy = *sx ;
      jj++ ; sy++ ;
    }
  arrayMax (introns) = jj ;
  return ;
} /* sxIntronsAddCounts */

/*************************************************************************************/

static void sxIntronCleanUp (SX *sx0, Array introns, DICT *dict)
{
  SXX *sx, *sy ;
  int ii, jj, n1, n2, pass, dx, dy ;
  int minS = sx0->minimalIntronSupport ;
  const char *ccp, *ccq ;

  /* in first round, regularize on shortest jump the donors with several very similar acceptors
   * then flip the coordinates and run again
   * this regularizes the acceptors with several donors
   * the reflip the coordinates to be back in the original frame
   */
  for (pass = 0 ; pass < 2 ; pass++)
    {
      for (ii = 0, sx = arrp (introns, 0, SXX) ; ii < arrayMax (introns) ; sx++, ii++)
	{
	  if (sx->m == -1)
	    continue ;
	  if (sx->eClass != 1)
	    {
	      /* add up the counts */
	      for (jj = ii + 1, sy = sx + 1 ;  jj < arrayMax (introns) && sy->a1 == sx->a1 && sy->eClass == sx->eClass; sy++, jj++)
		{  
		  if (sx->eClass == 3 && sx->nAcceptorSupport != sy->nAcceptorSupport )
		    continue ;
		  sy->m = -1 ;
		  sx->nDonorSupport += sy->nDonorSupport ;
		  if (sx->eClass != 3) sx->nAcceptorSupport += sy->nAcceptorSupport ;
		}
	      continue ;
	    }
	  if ((sx->type & 0x3) != 0x3)
	    continue ;
	  dx = sx->a2 - sx->a1 ;
	  for (jj = ii + 1, sy = sx + 1 ;  jj < arrayMax (introns) && sy->a1 == sx->a1 ; sy++, jj++)
	    {
	      if (sy->eClass != 1)
		continue ;
	      dy = sy->a2 - sy->a1 ;
	      if (dx == dy)
		{ /* eliminate doublets */
		  sy->m = -1 ;
		  sx->nDonorSupport += sy->nDonorSupport ;
		  sx->nAcceptorSupport += sy->nAcceptorSupport ;
		}
	      if (dy < 3 * dx)
		continue ; /* the idea is to separate repeated genes */
	      n1 = n2 = 0 ;
	      if ((pass == 0 && ( sx->type & 0x4)) || (pass == 1 && (sx->type & 0x8)))
		{
		  ccp = dictName (dict, sx->acc) ;
		  ccq = dictName (dict, sy->acc) ;
		}
	      else if ((pass == 0 && (sx->type & 0x8)) || (pass == 1 && (sx->type & 0x4)))
		{
		  ccp = dictName (dict, sx->donor) ;
		  ccq = dictName (dict, sy->donor) ;
		}
	      else
		continue ;
	      while (*ccp && *ccq)
		{
		  if (*ccp == *ccq) n1++ ;
		  else n2++ ;
		  ccp++ ; ccq++ ;
		}
	      if (5*n2 < n1+n2 && n2 < 4) /* too close to call */
		{ /* eliminate longest */
		  sy->m = -1 ;
		}
	    }
	}

      /* keep happy few  and inverse the coordinates */
      for (ii = jj = 0, sx = sy = arrp (introns, 0, SXX) ;  ii < arrayMax (introns) ; sx++, ii++)
	{
	  if (sx->m != -1 &&  
	      (
	       (sx->eClass == 1 && sx->nDonorSupport >= minS && sx->nAcceptorSupport >= minS) ||
	       (sx->eClass != 1 && sx->nDonorSupport + sx->nAcceptorSupport >= minS)
	       )
	      )
	    {
	      if (ii != jj)
		*sy = *sx ;
	      jj++ ; sy++ ;	  
	    }
	}
      arrayMax (introns) = jj ;
      for (ii = 0, sx = arrp (introns, 0, SXX) ;  ii < arrayMax (introns) ; sx++, ii++)
	{
	  n1 = sx->a1 ;
	  sx->a1 = - sx->a2 ;
	  sx->a2 = - n1 ;
	}
      arraySort (introns, sxxOrder) ;
      sxIntronsAddCounts (introns) ;
    }

  return ;
} /* sxIntronCleanUp */

/*************************************************************************************/

static Array sxNewExonsGetIntrons (SX *sx, BOOL isNew, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  Array introns = arrayHandleCreate (100000, SXX, h0) ;
  ACEIN ai = 0 ;
  int jj = 0, a1, a2, donor, acc, nds, nas, nsl = 0, eClasse ; 
  const char *ccp ;
  BOOL isDown, ispA, isExon, isSL ;
  SXX *sxx ;
  DICT *dict = dictHandleCreate (10000, h) ;

  if (isNew && sx->sxxNewIntronsFileName)
    ai = aceInCreate (sx->sxxNewIntronsFileName, 0, h) ;
  else if (!isNew && sx->sxxKnownIntronsFileName)
    ai = aceInCreate (sx->sxxKnownIntronsFileName, 0, h) ;
  if (ai)
    {
      /* aceInSpecial (ai, "\n\t") ; */
      while (aceInCard (ai))
	{
	  eClasse = 0 ; ispA = isSL = isExon = FALSE ;
	  ccp = aceInWord(ai) ;
	  if (!ccp || !*ccp)
	    continue ;
	  if (!strcmp (ccp, "INTRON"))
	    eClasse = 1 ;
	  else if (!strcmp (ccp, "pA"))
	    { eClasse = 2 ; ispA = TRUE ; }
	  else if (!strcmp (ccp, "pa")) /* weak poly A */
	    { eClasse = 22 ; ispA = TRUE ; }
	  else if (!strcmp (ccp, "paw")) /* very weak poly A */
	    { eClasse = 23 ; ispA = TRUE ; }
	  else if (!strncmp (ccp, "SL", 2))
	    { eClasse = 3 ; isSL = TRUE ; nsl = 0 ; if (ccp[2]) nsl = ccp[2] - '0' ;  if (ccp[3]){nsl = 10 * nsl +  (ccp[3] - '0') ; }}
	  else if (!strcmp (ccp, "EXON"))
	    { eClasse = 4 ; isExon = TRUE ; }
	  else if (!strcmp (ccp, "CDS"))
	    { eClasse = 5 ; isExon = TRUE ; }
	  else if (!strcmp (ccp, "Primer"))
	    { eClasse = 6 ; isExon = TRUE ; }
	  else if (!strcmp (ccp, "Stop"))
	    { eClasse = 7 ; isExon = TRUE ; }
	  else
	    continue ;

	  aceInStep(ai,'\t') ; ccp = aceInWord(ai) ;
	  if (!ccp || !*ccp)
	    continue ;	
	  if (!strcmp (ccp, "Forward"))
	    isDown = TRUE ;
	  else if (!strcmp (ccp, "Reverse"))
	    isDown = FALSE ;
	  else
	    continue ;
	  
	  aceInStep(ai,'\t') ; ccp = aceInWord(ai) ;
	  if (!ccp || !*ccp || (sx->sxxChromosome && strcmp (ccp, sx->sxxChromosome)))
	    continue ;	
	  aceInStep(ai,'\t') ; if (! aceInInt (ai, &a1))
	    continue ; 
	  if (1)
	    {
	      if (eClasse == 1)
		{ /* chromosome_name forrward/reverse is repeated because it may be different in a rearr */
		  aceInStep(ai,'\t') ; ccp = aceInWord(ai) ;
		  if (!ccp || !*ccp || (sx->sxxChromosome && strcmp (ccp, sx->sxxChromosome)))
		    continue ;	
		  aceInStep(ai,'\t') ; ccp = aceInWord(ai) ;
		  if (!ccp || !*ccp)
		    continue ;	
		}
	      aceInStep(ai, '\t') ;
	      if (! aceInInt (ai, &a2))
		continue ;
	    }
	  aceInStep(ai, '\t') ;
	  donor = acc = 0 ;
	  nds = nas = 0 ;
	  if (ispA || isSL)
	    {
	      aceInStep(ai,'\t') ; ccp = aceInWord(ai) ; /* method */
	      if (ccp)
		{
		  dictAdd (dict, ccp+1, &donor) ;
		}
	      aceInStep(ai,'\t') ; aceInInt (ai, &nds) ; /* number of pA  support */
	    }
	  else if (isNew)
	    {
	      aceInStep(ai,'\t') ;   ccp = aceInWord(ai) ; /* length ignore */
	      aceInStep(ai,'\t') ;   ccp = aceInWord(ai) ; /* gt_ag ignore */
	      aceInStep(ai,'\t') ;   ccp = aceInWord(ai) ; /* SUPPORT ignore */
	      aceInStep(ai,'\t') ;   ccp = aceInWord(ai) ; /* -  ignore */
	      aceInStep(ai,'\t') ; aceInInt (ai, &nds) ; /* number of donor support */
	      aceInStep(ai,'\t') ; aceInInt (ai, &nas) ; /* number of acceptor support */
	      aceInStep(ai,'\t') ;   ccp = aceInWord(ai) ; /* 0  ignore */
	      aceInStep(ai,'\t') ;   ccp = aceInWord(ai) ; /* 0  ignore */
	      aceInStep(ai,'\t') ;   ccp = aceInWord(ai) ; /* Sequence  ignore */
	      if (! ccp || strcmp (ccp, "Sequence"))
		  continue ;
	      aceInStep(ai,'\t') ; ccp = aceInWord(ai) ; /* donor */
	      if (! ccp) continue ;
	      dictAdd (dict, ccp, &donor) ;

	      aceInStep(ai,'\t') ; ccp = aceInWord(ai) ; /* acceptor */
	      if (! ccp) continue ;
	      dictAdd (dict, ccp, &acc) ;
	    }
	  if (
	      (a1 < a2 && (sx->minX > a2 || (sx->maxX && sx->maxX < a1))) ||
	      (a1 > a2 && (sx->minX > a1 || (sx->maxX && sx->maxX < a2)))
	      )	      
	    continue ;

	  sxx = arrayp (introns, jj++, SXX) ;
	  sxx->isDown = isDown ;

	  if (isExon) /* oriented  along chromosome */
	    {
	      sxx->a1 = a1 <= a2 ? a1 : a2 ;
	      sxx->a2 = a1 > a2 ? a1 : a2 ;
	    }
	  else
	    {  /* oriented in their own direction */
	      sxx->a1 = isDown ? a1 : a2 ;
	      sxx->a2 = isDown ? a2 : a1 ;
	    }

	  /* ATTENTION EXONS and STOP are neither donors nor acceptors */
	  if ((ispA && isDown) || (isSL && !isDown))
	    sxx->type |= 0x1 ;
	  else if ((ispA && !isDown) || (isSL && isDown))
	    sxx->type |= 0x2 ;
	  else if (eClasse == 1)
	    sxx->type |= 0x3 ;
	  if (isDown)
	    sxx->type |= 0x4 ;
	  else
	    sxx->type |= 0x8 ;
	  
	  sxx->eClass = eClasse ;
	  
	  sxx->donor = donor ; /* donor sequence or method name */
	  sxx->acc = acc ;
	  sxx->nDonorSupport = nds ;
	  sxx->nAcceptorSupport = isSL ? nsl : nas ;
	}
    }
  arraySort (introns, sxxOrder) ;
  if (isNew)
    sxIntronCleanUp (sx, introns, dict) ;
  else
    sxIntronsAddCounts (introns) ;
  ac_free (h) ;
  return introns ;
} /* sxNewExonsGetIntrons */

/*************************************************************************************/

static int sxNewExonsRemoveCDSpA (Array candidates, Array knownElements, Array newElements, BOOL isDown)
{
  int nn = 0, ii, jj,
    iMax = candidates ? arrayMax (candidates) : 0,
    jMax = knownElements ? arrayMax (knownElements) : 0 ;
  int iNew = arrayMax (newElements) ;
  SXX *sx, *sy, *snew ;
  
  if (! iMax || !jMax)
    return 0 ;

  arraySort (candidates, sxxOrder) ;
  arraySort (knownElements, sxxOrder) ;

  /* 2009_09_06: i do not understand why i have to run once on each strand, but
   *  there seem otherwise to be a problem while adavancing sx++ and sy++, i jump some slots 
   */
  
  /* promote the reliability of the pA when possible */
  for (ii = 0, sx = arrp (candidates, 0, SXX) ; ii < iMax ; ii++, sx++)
    {
      for (sy = sx + 1, jj = ii + 1 ; sx->eClass == 2 && sx->isDown == sy->isDown && sy->a1 == sx->a1 && jj < iMax ; sy++, jj++)
	{
	  if (sy->eClass == 22 || sy->eClass == 23) sy->eClass = sx->eClass ;
	}
    }
  for (ii = 0, sx = arrp (candidates, 0, SXX) ; ii < iMax ; ii++, sx++)
    {
      for (sy = sx + 1, jj = ii + 1 ; sx->eClass == 22 && sx->isDown == sy->isDown && sy->a1 == sx->a1 && jj < iMax ; sy++, jj++)
	{
	  if (sy->eClass == 23) sy->eClass = sx->eClass ;
	}
    }
  
  sx = arrp (candidates, 0, SXX) ; ii = 0 ;
  sy = arrp (knownElements, 0, SXX) ; jj = 0 ;
  while (ii <= iMax && jj <= jMax)
    {
      while (ii < iMax && ((sx->eClass != 2 && sx->eClass != 22 && sx->eClass != 23) ||  sx->isDown != isDown))
	{ ii++ ; sx++ ;}
      while (jj < jMax && sy->eClass != 5 && sy->eClass != 6) { jj++ ; sy++ ;} /* CDS or Primer */
      if (sx->a1 > sy->a1 - 5 && sx->a1 < sy->a2 + 5)
	/* note that on purpose we remove the test on both strands */
	{ 
	  if (ii <= iMax && sx->eClass < 100)
	    {
	      sx->eClass = 101 ;
	      snew = arrayp (newElements, iNew++, SXX) ;
	      nn++ ;
	      *snew = *sx ;
	    }
	  ii++ ; sx++ ;
	}
      else if (sx->a1 > sy->a1)
	{
	  jj++ ; sy++ ;
	}
      else
	{
	  ii++ ; sx++ ;
	}
    }

  /* now remove if pA allready inside an exon even if not CDS, but on the same strand */
  sx = arrp (candidates, 0, SXX) ; ii = 0 ;
  sy = arrp (knownElements, 0, SXX) ; jj = 0 ;
  while (ii <= iMax && jj <= jMax)
    {
      while (ii < iMax && ((sx->eClass != 2 && sx->eClass != 22 && sx->eClass != 23) || sx->isDown != isDown))
	{ ii++ ; sx++ ;}
      while (jj < jMax && sy->eClass != 4) { jj++ ; sy++ ;}
      if (
	  (sx->type & 0xc) == (sy->type & 0xc) &&
	  sx->a1 >= sy->a1 -5 && sx->a1 <= sy->a2 + 5
	  )
	/* note that on purpose we remove the test on both strands */
	{ 
	  if (ii <= iMax && sx->eClass < 100)
	    {
	      sx->eClass = 102 ;
	      snew = arrayp (newElements, iNew++, SXX) ;
	      nn++ ;
	      *snew = *sx ;
	    }
	  ii++ ; sx++ ;
	}
      else if (sx->a1 > sy->a1)
	{
	  jj++ ; sy++ ;
	}
      else
	{
	  ii++ ; sx++ ;
	}
    }
  /* now remove vewr weak pA allready inside an exon on the oppositestrand */
  sx = arrp (candidates, 0, SXX) ; ii = 0 ;
  sy = arrp (knownElements, 0, SXX) ; jj = 0 ;
  while (ii <= iMax && jj <= jMax)
    {
      while (ii < iMax && (sx->eClass != 23 || sx->isDown != isDown))
	{ ii++ ; sx++ ;}
      while (jj < jMax && sy->eClass != 4) { jj++ ; sy++ ;}
      if (sx->a1 > sy->a1 - 5 && sx->a1 < sy->a2 + 5)
	/* note that on purpose we remove the test on both strands */
	{ 
	  if (ii <= iMax && sx->eClass < 100)
	    {
	      sx->eClass = 103 ;
	      snew = arrayp (newElements, iNew++, SXX) ;
	      nn++ ;
	      *snew = *sx ;
	    }
	  ii++ ; sx++ ;
	}
      else if (sx->a1 > sy->a1)
	{
	  jj++ ; sy++ ;
	}
      else
	{
	  ii++ ; sx++ ;
	}
    }
  return nn ;
} /* sxNewExonsRemoveCDSpA */

/*************************************************************************************/
/* add in each wiggle entry the local median */
static void sxNewExonsRollingMedian (int width, Array wiggle)
{
  SXW *sxw, *sxw2 ;
  int old, new, i, j, ii, jj, median = 0, x ;

  int NN = (width > 0 ? 2*(width/2) + 1 : 5) ; /* width for the median window */
  int nn[NN] ;  /* historic of the last NN elements */
  int s[NN] ;   /* sorting  of the last NN elements */

  fprintf (stderr, "%s %s\n", "median start ", timeShowNow()) ;
   for (i = 0 ; i < NN ; i++)
     s[i] = ACEDB_MAXINT ;

  /* compute the rolling median following a neat idea by Vahan */
  for (ii = 0, sxw = arrp (wiggle, 0, SXW), sxw2 = sxw - NN/2 ; ii < arrayMax (wiggle) ; ii++, sxw++, sxw2++)
    {
      if (0)
	{
	  if (ii>NN && s[0]!=nn[0] && s[0] != nn[1] && s[0]!=nn[2])
	    invokeDebugger () ;
	  if (ii>NN && s[1]!=nn[0] && s[1] != nn[1] && s[1]!=nn[2])
	    invokeDebugger () ;
	  if (ii>NN && s[2]!=nn[0] && s[2] != nn[1] && s[2]!=nn[2])
	    invokeDebugger () ;
	}
      new = sxw->n ;
      jj = ii % NN ;
      old = nn[jj] ;
      nn[jj] = new ;
      if (ii >= NN && old == new) 
	{ sxw2->m = median ; }
      else
	{
	  /* delete old */
	  if (ii >= NN)
	    for (i = 0 ; i < NN ; i++)
	      {
		x = s[i] ;
		if (x == old)
		  {
		    for (j = i ; j < NN - 1 ; j++)
		      s[j] = s[j+1] ;
		    break ;
		  }
	      }
	  /* insert */
	  if (new >= s[NN-2])
	    s[NN-1] = new ;
	  else
	    for (i = 0 ; i < NN - 1 ; i++)
	      {
		x = s[i] ;
		if (x > new)
		  {
		    for (j = NN - 1 ; j > i ; j--)
		      s[j] = s[j-1] ;
		    s[j] = new ;
		    break ;
		  }
	      }
	  if (ii >= NN)
	    sxw2->m = median  = s[NN/2] ;
	}
    }

  fprintf (stderr, "%s %d   %s\n", "median done ", ii, timeShowNow()) ;
  return ;
} /* sxNewExonsRollingMedian */

/*************************************************************************************/
/* export a set of walls centered at 5,15..... 5+10x where the median of the wiggle goes up or down
 * in revise == 1 mode, we subdivide if there is a factor 4 between third neighbours 
 */
static Array sxNewExonsGetWalls (SX *sx, Array wallsOld, Array wiggle, int revise, BOOL debug, AC_HANDLE h0)
{
  Array walls = 0 ;
  int new, old = 0, old3 = 0, old2 = 0, old1 = 0, ii, j, jj, aOld ;
  int exonReviseLevel = 3 * sx->exonLevel ;
  SXW *sxw ;
  SXX *sxx, *syy ;
  BOOL debug2 = FALSE ;

  if (revise == 0) exonReviseLevel = sx->exonLevel ;
  if (revise == 1) exonReviseLevel = sx->exonLevel * 3 ;

  walls = arrayHandleCreate (100000, SXX, h0) ;
  /* create a wall if the median signal goes abruptly up or down or there is an intron or a pA */
  jj = 0 ;
  sxx = arrayp (walls, jj++, SXX) ; /* create a wall at position zero */
  sxx->n = 1 ;
  sxx->a1 = sxx->a2 = aOld = new = 0 ; /* create a wall at position zero */
  for (ii = 0, sxw = arrp (wiggle, 0, SXW), old = sxw->m ; ii < arrayMax (wiggle) ; aOld = sxw->a1, ii++, sxw++)
    {
      new = sxw->m ; 
      if (debug2) fprintf (stderr, ".. %d %d\t", sxw->a1, new) ;
      if (new < exonReviseLevel) new = 0 ;
      if (new > 4*old && (revise == 0 || old > exonReviseLevel))
	{ 
	  if (debug2) fprintf (stderr, "up-wall") ;
	  sxx = arrayp (walls, jj++, SXX) ; 
	  sxx->a1 = sxw->a1 - (revise ? 15 : 5) ;
	}
      if (4*new < old && (revise == 0 || new > exonReviseLevel))
	{
	  sxx = arrayp (walls, jj++, SXX) ; 
	  sxx->a1 = sxw->a1 - (revise ? 15 : 5) ;
	}
      if (revise)
	{ 
	  if (debug2) fprintf (stderr, "down-wall") ;
	  old3 = old2 ; old2 = old1 ; old1 = new ;
	  old = old3 ;
	}
      else
	old = new ; 
      if (debug2) fprintf (stderr, "\n") ;
    }

  if (wallsOld) /* only add walls not too close to previous walls */
    {
      if (! arrayMax (wallsOld))
	{ /* just copy */
	  for (ii = jj - 1 
		 , sxx = arrp (walls, ii, SXX)
		 , syy = arrayp (wallsOld, ii, SXX) 
	     ;   ii >= 0 ; sxx--, syy--, ii--)
	    *syy = *sxx ;
	}
      else if (arrayMax (walls))
	{
	  int a1 ;

	  /* elminate redundant new walls */
	  for (ii = 0, sxx = arrayp (walls, ii, SXX)  ; ii < arrayMax (walls) ; ii++, sxx++)
	    {
	      a1 = sxx->a1 ;
	      if (a1 <= 0) continue ;
	      for (jj = ii, syy = arrp (walls, jj, SXX) ; jj < arrayMax (walls) && syy->a1 < a1+25 && syy->a1 > 0 ; jj++, syy++)
		a1 = syy->a1 ;
	      for (j = ii + (jj - ii)/2 + 1 ; j < jj ; j++)
		arr (walls, j, SXX).a1 = -1 ;
	      for (j = ii + (jj - ii)/2 - 1 ; j >= ii ; j--)
		arr (walls, j, SXX).a1 = -1 ;
	    }

	  for (ii = jj = 0, sxx = arrayp (walls, ii, SXX)  ; ii < arrayMax (walls) ; ii++, sxx++)
	    {
	      if (sxx->a1 <= 0)
		continue ;
	      /* search for a close by old walls */
	      a1 = sxx->a1 - 20 ;
	      for (j = jj, syy = arrp (wallsOld, j, SXX) ; j >= 0 && syy->a1 > a1 ; syy--, j--) ;
	      if (syy->a1 >= a1) 
		sxx->a1 = -1 ; /* kill */
	      else 
		{
		  a1 = sxx->a1 + 20 ;
		  for ( ; j < arrayMax (wallsOld) && syy->a1 < a1 ; syy++, j++) ;
		  if (syy->a1 <= a1) 
		    sxx->a1 = -1 ; /* kill */
		}
	    }
	  /* add the surviving new walls in wallsOld */
	  for (ii = 0, jj = arrayMax (wallsOld), sxx = arrp (walls, ii, SXX) 
		 ; ii < arrayMax (walls) ; ii++, sxx++)
	    if (sxx->a1 > 0)
	      {
		syy = arrayp (wallsOld, jj++, SXX)  ;
		*syy = *sxx ;
	      }
	  arraySort (wallsOld, sxxOrder) ;
	  arrayCompress (wallsOld) ;
	}
      arrayDestroy (walls) ;
    }
  else
    {
      sxx = arrayp (walls, jj++, SXX) ;
      sxx->a1 = aOld + 5 ; /* create a wall at last position */
      
      wallsOld = walls ;
    }

  if (debug)
    {
      fprintf (stderr, "Walls\n") ;
      for (ii =  0 ; ii < arrayMax (wallsOld) ; ii++)
	{
	  sxx = arrp (wallsOld, ii, SXX) ;
	  fprintf (stderr, "%d: %d\n", ii, sxx->a1) ;
	}
    }
  
  return wallsOld ;
} /* sxNewExonsGetWalls */

/*************************************************************************************/
/* export a set of walls where there is a donor/acceptor/SL/polyA 
 * and push the coord by 1bp inside the exon if the coord is divisible by 10
 */
static Array sxNewExonsGetIntronsWalls (SX *sx, Array walls, Array knownElements, Array candidates, AC_HANDLE h0)
{
  int a1, ii, jj, pass ;
  SXX *sxw ;
  SXX *sxx ;
  Array introns ;

  /* create a wall at each intron boundary */
  jj = arrayMax (walls) ;
  for (pass = 0 ; pass < 2 ; pass++)
    {
      introns = pass ? knownElements : candidates ;
      if (introns && arrayMax(introns))
	for (ii = 0, sxx = arrp (introns, 0, SXX) ;  ii < arrayMax (introns) ; ii++, sxx++)
	  {
	    if (sxx->eClass == 22 || sxx->eClass == 23) /* do not flag a reflecting wall for a weak polyA */
	      { sxx->eClass = 2 ; continue ; }

	    if (sxx->eClass == 1)  /* an intron */
	      {
		a1 = sxx->a1 + 0 ; /* so the whole exon is before the boundary */ 
		a1 = a1/10 ; a1 = 10 * a1 ;
		a1 += 5 ;  
		
		sxw = arrayp (walls, jj++, SXX) ;
		sxw->type = 0x1 | (sxx->type & 0xc) ; /* 0xc = 0x4 | 0x8 == strand mask */
		sxw->eClass = 1 ;
		sxw->a1 = a1 ;      /* rounded position */
		sxw->a2 = sxx->a1 ; /* original position */
		
		a1 = sxx->a2 - 1 ; /* so the whole exon is after the boundary */ 
		a1 = a1/10 ; a1 = 10 * a1 ;
		a1 += 5 ; 
		
		sxw = arrayp (walls, jj++, SXX) ;
		sxw->type = 0x2 | (sxx->type & 0xc) ; /* 0xc = 0x4 | 0x8 == strand mask */
		sxw->eClass = 1 ;
		sxw->a1 = a1 ;      /* rounded position */
		sxw->a2 = sxx->a2 ; /* original position */
	      }
	    else if ((sxx->eClass == 2 || sxx->eClass == 3 || sxx->eClass == 7) && (sxx->type & 0x1))
	      {  /* up SL or down pA or down stop */
		a1 = sxx->a1 + 0 ; /* so the whole exon is before the boundary */ 
		a1 = a1/10 ; a1 = 10 * a1 ;
		a1 += 5 ;  
		
		sxw = arrayp (walls, jj++, SXX) ;
		sxw->type = sxx->type ;
		sxw->eClass = sxx->eClass ;
		sxw->a1 = a1 ;      /* rounded position */
		sxw->a2 = sxx->a1 ; /* original position */
	      }
	    else if ((sxx->eClass == 2 || sxx->eClass == 3 || sxx->eClass == 7) && (sxx->type & 0x2))
	      {   /* down SL or up pA or up Stop */
		a1 = sxx->a2 - 1 ; /* so the whole exon is after the boundary */ 
		a1 = a1/10 ; a1 = 10 * a1 ;
		a1 += 5 ; 
		
		sxw = arrayp (walls, jj++, SXX) ;
		sxw->type = sxx->type ;
		sxw->eClass = sxx->eClass ;
		sxw->a1 = a1 ;      /* rounded position */
		sxw->a2 = sxx->a2 ; /* original position */
	      }
		
	  }
    }

  sxw = arrayp (walls, jj++, SXX) ;
  sxw->a1 = sxw->a2 = ACEDB_MAXINT ;

  arraySort (walls, sxxOrder) ;
  arrayCompress (walls) ;
 return walls ;
  /* destroy the segment immediatly inside the introns */
  for (ii = 0, sxw = arrp (walls, 0, SXX) ;  ii < arrayMax (walls) ; ii++, sxx++)
    {
      if (sxw->type & 0x2) /* accepteur */
	{
	  /* destroy wiggle walls very close upstream */
	  for (jj = ii - 1, sxx = sxw - 1 ; sxx->a1 > sxw->a1 - 15 && jj > 0 ; jj--, sxx--)
	    if (! sxx->type)
	      sxx->eClass = 999 ;
	}
      if (sxw->type & 0x1) /* donor */
	{
	  /* destroy wiggle walls very close downstream */
	  for (jj = ii + 1, sxx = sxw + 1 ; sxx->a1 < sxw->a1 + 15 && jj < arrayMax (walls) ; jj++, sxx++)
	    if (! sxx->type)
	      sxx->eClass = 999 ;
	}
    }
  /* keep happy few */
  for (ii = jj = 0, sxx = sxw = arrp (walls, 0, SXX) ;  ii < arrayMax (walls) ; ii++, sxx++)
    {
      if (sxx->eClass != 999)
	{
	  if (jj < ii)
	    *sxw = *sxx ;
	  sxw++ ; jj++ ;
	}
    }
  arrayMax (walls) = jj ;
 
  return walls ;
} /* sxNewExonsGetIntronsWalls */

/*************************************************************************************/

static void sxNewExonsMergeWeakPa (Array candidates)
{
  SXX *sxx, *sxw, *sxMax ;
  int ii, jj = 0, a1 ;
  int found = 0, merged = 0, kept = 0, killed = 0 ;
  int nx, nMax, minClasse ;

  for (ii = 0, sxw = arrp (candidates, 0, SXX) ; ii < arrayMax (candidates); ii++, sxw++)
    if (sxw->eClass == 2 || sxw->eClass == 22 || sxw->eClass == 23) /* pA */
      {
	a1 = sxw->a1 ; nMax = nx = 0 ; sxMax = 0 ; minClasse = sxw->eClass ;
	for (sxx = sxw, jj = ii ; jj < arrayMax (candidates) && sxx->a1 < a1 + 20 ; jj++, sxx++)
	  if ((sxx->eClass == 2 || sxx->eClass == 22 || sxx->eClass == 23))
	    {
	      found++ ; nx++ ;
	      if (sxx->eClass <= minClasse && sxx->nDonorSupport >= nMax)
		{ nMax = sxx->nDonorSupport ; sxMax = sxx ; }
	    }
	if (nx > 1)
	  for (sxx = sxw, jj = ii ; jj < arrayMax (candidates) && sxx->a1 < a1 + 20 ; jj++, sxx++)
	    if ((sxx->eClass == 2 || sxx->eClass == 22 || sxx->eClass == 23) &&  sxMax != sxx)
	      { sxMax->nDonorSupport += sxx->nDonorSupport ; sxx->eClass = 999 ; merged++ ; }
	if (sxMax)
	  {
	    if (sxMax->nDonorSupport < 5)
	      { sxMax->eClass = 999 ; killed++ ; }
	    else 
	      { sxMax->eClass = -sxMax->eClass ; kept++ ; }
	  }
      }

  /* keep happy few */
  for (ii = jj = 0, sxx = sxw = arrp (candidates, 0, SXX) ; ii < arrayMax (candidates); ii++, sxw++)
    {
      if (sxw->eClass < 0)
	sxw->eClass = - sxw->eClass ;
      if (sxw->eClass != 999)
	  {
	    if (jj < ii) 
	      *sxx = *sxw ;
	    jj++ ; sxx++ ;
	  }
    }
  arrayMax (candidates) = jj ;
  if (jj)
    {     
      arraySort (candidates, sxxOrder) ;
      arrayCompress (candidates) ;
      
      fprintf (stderr, "sxNewExonsMergeWeakPa found %d pA sites, merged %d cases, killed %d sites, kept %d merged sites\n", found, merged, killed, kept) ;
    }
  return ;
} /* sxNewExonsMergeWeakPa */

/*************************************************************************************/
/* create a segments in between all the walls, falling at 10, 20, ... modulo 10 */
static Array  sxNewExonsSegments (Array walls, AC_HANDLE h0)
{
  int ii, jj, old ;
  SXX *sxx ;
  SXX *sxw ;
  int iMax = arrayMax (walls) ;
  Array segs = arrayHandleCreate (1000, SXX, h0) ;
  
  old = -10 ;
  for (ii = jj = 0, sxw = arrp (walls, 0, SXX) ; ii < iMax ; sxw++, ii++)
    {
      if (sxw->eClass > 0) /* SL or pA */
	continue ;
      if (old >=  sxw->a1)
	continue ;
      sxx = arrayp (segs, jj++, SXX) ;
      sxx->a1 = old ;
      old = sxw->a1/10 ; old *= 10 ; /* the walls are at 5,15,25 ... */
      sxx->a2 = old ;
      old += 10 ;
    }

  /* add a length 0 seg at each intron/SL/pA/Stop wall position */
  for (ii = 0, sxw = arrp (walls, 0, SXX) ; ii < iMax ; sxw++, ii++)
    {
      if (sxw->eClass == 0) /* wiggle wall */
	continue ;
      sxx = arrayp (segs, jj++, SXX) ;
      sxx->type = sxw->type ;
      sxx->eClass = sxw->eClass ;
      sxx->a1 = sxx->a2 = sxw->a2 ; /* original position */
    }

  arraySort (segs, sxxOrder) ;
  arrayCompress (segs) ;

  return segs ;
} /* sxNewExonsSegments */

/*************************************************************************************/
/* for each segment between 2 walls, count the hits */
static void sxNewExonsSegmentsSignal (Array wiggle, Array segs, BOOL debug)
{
  int ii, jj, kk ;
  SXX *sxx ;
  SXW *sxw ;
  int jMax = arrayMax (wiggle) ;
  Array aa = arrayCreate (1000, int) ;

  sxw = arrp (wiggle, 0, SXW) ;
  for (ii = jj = 0, sxx = arrp (segs, 0, SXX) ; ii < arrayMax(segs) ; sxx++, ii++)
    {
      if (sxx->type)
	continue ;
      kk = 0 ;
      while (jj < jMax && sxw->a1 < sxx->a1) 
	{ jj++ ; sxw++ ; }
      while (jj < jMax && sxw->a1 <= sxx->a2) 
	{ array (aa, kk++, int) = sxw->n ; jj++ ; sxw++ ; }
      if (kk)
	{
	  arrayMax (aa) = kk ;
	  if (kk > 1) arraySort (aa, intOrder) ;
	  sxx->m = arr (aa, kk/2, int) ; /* the median */
	}
    }
  arrayDestroy (aa) ;
  if (debug) showSXX(segs) ;
  return ;
} /* sxNewExonsSegmentsSignal */

/*************************************************************************************/
/* add in new elements the probable new exons 
 * i.e. a segment which is
 *    higher then newExonLevel
 *    its 2 neighbours are either higher but no more than 20*myself or lower than myself/20
 *    add it on both strand unless the whole program is stranded
 */

static int sxNewExonsAddNewExons (SX *sx, Array segs, Array candidates, Array knownElements, Array newElements)
{
  SXX *sxx, *sxa, *sxb, *sxe, *sxw1 ;
  int ok, a1, i, ii, jj, m, ma, mb ;
  int k1 = 0, kMax1 = knownElements ? arrayMax (knownElements) : 0 ;
  int k2 = 0, kMax2 = newElements ? arrayMax (newElements) : 0 ;
  int nn = 0, newExonLevel = sx->exonLevel ;
  int bLevel =  newExonLevel ;

  jj = arrayMax (candidates) ;
  for (ii = 1, sxx = arrp (segs, 1, SXX) ; ii < arrayMax (segs) - 1 ; ii++, sxx++)
    {
      if (sxx->eClass)
	continue ;
      m = sxx->m ;
      if (m < newExonLevel)
	continue ;

      bLevel =  newExonLevel ; 
      if (0) /* 2015_06_20 je ne comprend rien a ce code, qui perd les exoncs les plus haut */
	{
	  if (0 && m > 20 * bLevel)
	    bLevel = m / 20 ;
	  for (i = ii -1, sxa = sxx - 1 ; i > 1 &&  (sxa->a2 - sxa->a1 < 0 || sxa->type) ; i--, sxa--) ;
	  ma = sxa->m ;
	  if (0 * ma > bLevel && ma < m) continue ;
	  mb = sxx->m ; sxe = sxx ;
	  for (i = ii + 1, sxb = sxx + 1 ; i < arrayMax (segs) - 1 &&  (sxa->a2 - sxa->a1 < 12 || sxb->type || sxb->m > mb) ; i++, sxb++) 
	    if (! sxb->type && sxa->a2 - sxa->a1 > 12 && sxb->m > mb) 
	      { sxe = sxb ; mb = sxb->m ; }
	  if (sxe > sxx) 
	    continue ;
	  if (mb > bLevel && mb < m) continue ;
	}
      a1 = (sxx->a1 + sxx->a2)/2 ;

      /* check if we are at least at 5% of the neighbouring genes */
      if (0)
	{
	  bLevel =  newExonLevel ;
	  for (i = ii,  sxw1 = sxx ; i >= 0 ; sxw1--, i--)
	    {
	      if (sxw1->a1 < a1 - 300)
		break ;
	      if (sxw1->m > 20 * bLevel)
		bLevel = sxw1->m / 20 ;
	    }
	  for (i = ii , sxw1 = sxx ; i <  arrayMax (segs) ; sxw1++, i++)
	    {
	      if (sxw1->a1 > a1 + 300)
		break ;
	      if (sxw1->m > 20 * bLevel)
		bLevel = sxw1->m / 20 ;
	    }
	  if (m < bLevel)
	    continue ;
	}

      /* check if there is already an exon on this strand */
      ok = 0 ;
      if (kMax1 > 1)
	{ 
	  if (k1 < 1) k1 = 1 ;
	  for (i = k1,  sxa = arrp (knownElements, k1, SXX) ; i > 0 && sxa->a1 > a1 ; i--, sxa--) ;
	  k1 -= 100 ; if (k1 < 0) k1 = 0 ;

	  for (i = k1,  sxa = arrp (knownElements, k1, SXX) ; i < kMax1 ; i++, sxa++)
	    {
	      if ((sxa->eClass == 4) && 
		  sxa->a1 <= a1+5 && sxa->a2 >= a1-5
		  )
		ok |= (sxa->type & 0xc) ;
	      if (sxa->a1 > a1)
		break ;
	    }
	}


      if (kMax2 > 1 && ok != 0xc)
	{ 
	  if (k2 < 1) k2 = 1 ;
	  for (i = k2,  sxa = arrp (newElements, k2, SXX) ; i > 0 && (sxa->b1 > a1 || sxa->a2 > a1) ; i--, sxa--)
	    k2 -= 100 ; 
	  if (k2 < 0) k2 = 0 ;

	  for (i = k2,  sxa = arrp (newElements, k2, SXX) ; i < kMax2 ; i++, sxa++)
	    {
	      if (sxa->b1 <= a1+5 && sxa->a1 >= a1-5)
		ok |= (sxa->type & 0xc) ;
	      if (sxa->a2 <= a1+5 && sxa->b2 >= a1-5)
		ok |= (sxa->type & 0xc) ;
	      if (sxa->a1 > a1)
		break ;
	    }
	}


      if (ok == 0xc) ; /* known on both strands, drop the new exon */
      else
	{ /* report on that strand */
	  sxe = arrayp (candidates, jj++, SXX) ; nn++ ;
	  sxe->eClass = 4 ; /* new exon */
	  sxe->type = 0x3 | ok | 0x10 ; /* this strand new exons */
	  sxe->a1 = a1 ;
	  sxe->a2 = a1 ;

	  if (ok == 0) /* report on the 2 strands even if we require -strand
			* so that assigning strand by inertia works
			* otherwise -strand then -antistrand would not be
			* equivalent to no strand option
			*/
	    {
	      sxe = arrayp (candidates, jj++, SXX) ; nn++ ;
	      *sxe = *(sxe-1) ;
	      sxe->type |= 0x8 ;
	      (sxe - 1)->type |= 0x4 ;
	    }
	}

      /* we could compute the local extension now, as in the else clause
       * but it is better to do it in the general code to forbid unwanted extensions into introns pA
       */
    }
      
    arraySort (newElements, sxxOrder) ;
  return nn ;
} /* sxNewExonsAddNewExons */

/*************************************************************************************/
/* create new objects and register them into knowIntrons */
static int sxNewExonsCreate (SX *sx, Array segs, Array introns, Array newElements, BOOL type, int eClass)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0, ii, iKnown, jj, kk = 0, a1, a2, b1, b2, b01, b02, jMax = arrayMax (segs) ;
  int aLevel, bLevel, iLevel ; /* extend only as long as the expression inceases */
  SXX *sxw, *sxw1, *sxi, *sxk, *szwa, *szwb ;
  int exonLevel = 0 ;
  int MAXXA = 1000 ; /* do not extend a pA upstream more than 1000bp, there should be a new exon */
  BOOL aadd ;

  switch (eClass)
    {
    case 1: /* intron */
      exonLevel = sx->exonLevel ;
      break ;
    case 2: /* polyA */
      exonLevel = sx->exonLevel/4 ;
      break ;
    case 3: /* SL */
      exonLevel = sx->exonLevel/4 ;
      break ;
    case 4: /* new exon */
      exonLevel = sx->exonLevel ;
      break ;
    case 5: /* known primer */
      return 0 ;
    case 6: /* Primer */
      return 0 ;
    case 7: /* new exon */
      exonLevel = 0 ;
      break ;
    }

  sxw = arrp (segs, 0, SXX) ; jj = 0 ;
  for (ii = jj = iKnown = 0, sxi = arrp (introns, 0, SXX) ; ii < arrayMax(introns) ; sxi++, ii++)
    {
      if (sx->minX > sxi->a1 || (sx->maxX && sx->maxX < sxi->a1)) 
	continue ;
      if (sxi->eClass != eClass || (sxi->type & 0xc) != type) /* treat successively each category */
	continue ;
      aLevel = bLevel = iLevel = a1 = a2 = b1 = b2 = b01 = b02 = 0 ; 
      szwa = szwb = 0 ;
      if (sxi->type & 0x1) /* try to extend upwards */
	{
	  kk = 0 ;
	  b1 = a1 = sxi->a1 ; /* donor or polyA */
	  /* locate the segment just right of the donor, in the intron */
	  while (jj < jMax && (sxw->a1 < a1 || sxw->type)) /* ignore artificial segments */ 
	    { jj++ ; sxw++ ; } 
	  if (eClass == 1) iLevel = sxw->m ;
	  else iLevel = -1 ;
	  for (kk = jj, sxw1 = sxw ; kk >= 0  && (sxw1->a1 > a1 || sxw1->type) ; sxw1--, kk--) ;
	  if (sxw1->a1 <= a1 && sxw1->a2 >= a1)
	    b1 = sxw1->a1 - 5 ;
	      /* extend upstream from the donor into the donor exon */

	  /* level of the intron or of the exon supporting the pA or the SL */
	  if (sxw1->m < exonLevel && 
	      sxw1->a1 < a1 - 20 && sxw1->a2 > a1 + 20
	      )
	    goto failed ;

	  aLevel = exonLevel ;
	  if (0 && eClass == 2) /* pA, verify we are not in the gresil of an intron */
	    {
	      for (kk = jj, sxw1 = sxw ; kk >= 0 ; sxw1--, kk--)
		{
		  if (sxw1->a1 < a1 - 1000)
		    break ;
		  if (sxw1->eClass == 1 && (sxw1->type  & 0x2)) /* intron accepteur */
		    break ;
		  if (sxw1->m > 20 * aLevel)
		    { aLevel = sxw1->m / 20 ; szwa = sxw1 ; break ; }
		}
	      for (kk = jj , sxw1 = sxw ; kk >= 0 ; sxw1++, kk++)
		{
		  if (sxw1->a1 > a1 + 200)
		    break ;
		  if (sxw1->m == 0 && ! sxw1->type) /* found a gap */
		    break ;
		  if (sxw1->m > 20 * aLevel)
		    { aLevel = sxw1->m / 20 ; szwa = sxw1 ; break ; }
		}
	    }
	  for (kk = jj - 1, sxw1 = sxw - 1 ; kk >= 0 ; sxw1--, kk--)
	    {
	      if (sxw1->type)
		continue ;
	      else if (sxw1->a1 > a1 - 12)
		continue ;
	      else if (sxw1->a2 < sxw1->a1  + 12)
	        continue ; 
	      else if (sxw1->a1 < a1  - 22 && sxw1->m < 3 * exonLevel && sxw1->m < aLevel) /* 20bp known from the construction stage */
		{ break ; }
	      else if (sxw1->a2 < a1 - 32) /* 20bp known from the construction stage */
		break ;
	      else if (sxw1->m > aLevel) 
		{ aLevel = sxw1->m ; b1 = sxw1->a1 ; szwa = sxw1 ; break ; }
	    }

	  /* check for a SL polyA inside the b1-a1 segment */
	  kk = jj ; sxw1 = sxw ;
	  for (; kk < jMax ; kk++, sxw1++)
	    if (sxw1->a1 > a1)
	      break ;
	  b01 = b1 + 1 ;
	  aadd = FALSE ;
	  while ( sxw1--, kk-- , kk > 0  && sxw1->a1 > b1)
	    {
	      if (!sxw1->type || sxw1->eClass > 20)
		continue ;
	      if (eClass != 4 && sxw1->a1 >= a1 - 10)
		continue ;
	      if (1) /* 2010_06_20 let the exon extend as they wish */
		{
		  /* AA DD rule, stop on Acceptor while Ascending, on Donor while Descending */
		  if ((type & 0x4) == (sxw1->type & 0x4) && (sxw1->type & 0x2)) /* acceptor or SL on my strand */
		    { aadd = TRUE ; b1 = sxw1->a1 ; break ; }
		  if ((type & 0x4) != (sxw1->type & 0x4) && (sxw1->type & 0x1)) /* donor or pA on opposite strand */
		    { aadd = TRUE ; b1 = sxw1->a1 ; break ; }
		  /* NMD prevention rule while descending from an intron, do not fuse into a same strand pA */
		  if (eClass == 1 && (type & 0x4) == (sxw1->type & 0x4) && 
		      sxw1->eClass == 2)
		    { aadd = TRUE ; b1 = sxw1->a1 ; break ; }
		}
	      /* greedy stops, stop 50bp upstream */
	      if (eClass == 2 && (type & 0x4) == (sxw1->type & 0x4) && 
		  sxw1->eClass == 7)
		{ if (b1 < sxw1->a1 - 50) b1 = sxw1->a1 - 50 ; break ; }
	    }

	  b1 += 1 ;  /* eat up 9 bp so we do NOT leak in next intron */
	  if (eClass != 4 && ! aadd) b1 += 8 ;
	  if (eClass != 4 && b1 >= a1 - 10)      /* no donor exon */
	    goto failed ;
	  if (eClass == 2 && b1 < a1 - MAXXA)
	    b1 = a1 - MAXXA ;
	}
      if (sxi->type & 0x2) /* extend downwards */
	{
	  b2 = a2 = sxi->a2 ; /* acceptor or SL */
	  /* locate the segment adjacent to the acceptor */
	  kk = jj ; sxw1 = sxw ;
	  while (kk > 0 && (sxw1->a1 > a2 || sxw1->type))
	    { kk-- ; sxw1-- ;}
	  /* sxw1 is the seg which includes the acceptor exon */
	  if (kk >= jMax) /* no acceptor exon */
	    goto failed ;
	  b2 = sxw1->a2 + 5 ;
	  /* level of the exon supporting the pA or the SL */
	  if (eClass != 1) iLevel = -1 ;
	  if (sxw1->m < exonLevel && 
	      sxw1->a1 < a2 - 20 && sxw1->a2 > a2 + 20
	      )
	    goto failed ;

	  bLevel = exonLevel ; szwb = 0 ;
	  if (0 && eClass == 2) /* pA, verify we are not in the gresil of an intron */
	    {
	      SXX *sxw2 ;

	      for (kk = jj , sxw2 = sxw ; kk >= 0 ; sxw2++, kk++)
		{
		  if (sxw2->a1 > a1 + 1000)
		    break ;
		  if (sxw2->eClass == 1 && (sxw2->type  & 0x1)) /* intron donor */
		    break ;
		  if (sxw2->m > 20 * bLevel)
		    { bLevel = sxw2->m / 20 ;  szwb = sxw2 ; break ; }
		}
	      for (kk = jj, sxw2 = sxw ; kk >= 0 ; sxw2--, kk--)
		{
		  if (sxw2->a1 < a1 - 200)
		    break ;
		  if (sxw2->m == 0 && ! sxw2->type) /* found a gap */
		    break ;
		  if (sxw2->m > 20 * bLevel)
		    { bLevel = sxw2->m / 20 ; szwb = sxw2 ; break ; }
		}
	    }

	  /* extend downstream into the acceptor exon */
	  kk = jj ; sxw1 = sxw ;
	  for (; kk > 0 && ( sxw1->a2 > a2 || sxw1->type) ; sxw1--, kk--) ;
	  for ( ; kk < jMax ; sxw1++, kk++)
	    {
	      if (sxw1->type)
		continue ;
	      else if (sxw1->a2 < a2 + 12)
		continue ;
	      else if (sxw1->a2 < sxw1->a1 + 12)
		continue ;
	      else if (sxw1->a2 > a2 + 22 && sxw1->m < 3 * exonLevel && sxw1->m < bLevel) /* 20bp known from the construction stage */
		break ;
	      else if (sxw1->a1 > a2 + 32) /* 20bp known from the construction stage */
		break ;
	      else if (sxw1->m > bLevel) 
		{ bLevel = sxw1->m ; b2 = sxw1->a2 ; szwb = sxw1 ; break ; }
	    }
	  /* check for a SL polyA inside the a2-b2 segment */
	  kk = jj ; sxw1 = sxw ;
	  for (; kk > 0 ; kk--, sxw1--)
	    if (sxw1->a2 < a2 && ! sxw1->type)
	      break ;
	  b02 = b2 - 1 ;
	  aadd = FALSE ;
	  while (sxw1++, kk++, kk < jMax && sxw1->a1 < b2)
	    {
	      if (!sxw1->type || sxw1->eClass > 20)
		continue ;
	      if (sxw1->a2 < a2+2)
		continue ;
	      if (eClass != 4 && sxw1->a1 <= a2 + 10)
		continue ;
	      if (1) /* 2010_06_20 let the exon extend as they wish */
		{
		  /* AA DD rule, stop on Acceptor while Ascending, on Donor while Descending */
		  if ((type & 0x4) != (sxw1->type & 0x4) && (sxw1->type & 0x2)) /* acceptor or SL opposite strand */
		    { aadd = TRUE ; b2 = sxw1->a1 + 0 ; break ; }
		  if ((type & 0x4) == (sxw1->type & 0x4) && (sxw1->type & 0x1)) /* donor or pA on my strand */
		    { aadd = TRUE ; b2 = sxw1->a1 + 0 ; break ; }
		  /* NMD prevention rule while descending from an intron, do not fuse into a same strand pA */
		  if (eClass == 1 && (type & 0x4) == (sxw1->type & 0x4) && 
		      sxw1->eClass == 2 && (sxw1->type & 0x2)) 
		    { aadd = TRUE ; b2 = sxw1->a1 + 0 ; break ; }
		}
	      /* greedy stops, stop 50bp upstream */
	      if (eClass == 2 && (type & 0x4) == (sxw1->type & 0x4) && 
		  sxw1->eClass == 7)
		{ if (b2 > sxw1->a1 + 50) b2 = sxw1->a1 + 50 ; break ; }
	    }
	  b2 -= 1 ;  /* eat up 9 bp so we do NOT leak in next intron */
	  if (! aadd && eClass != 4) b2 -= 8 ;
	  if (eClass != 4 && b2 < a2 + 12)      /* no acceptor exon */
	    goto failed ;
	  if (eClass == 2 && b2 > a2 + MAXXA)
	    b2 = a2 + MAXXA ;
	}
      if (eClass == 7 && sxi->type & 0x4) /* down peche a la ligne a partir d'un stop */
	{
	  /* we do not care about the wiggle */
	  a1 = b1 = sxi->a1 ; a2 = sxi->a2 ; b2 = a2 + 500 ;
	  /* find the position of the stop */
	  while (jj > 0 && (sxw->a1 > a2 || sxw->type))
	    { jj-- ; sxw-- ;}
	  while (jj < jMax && sxw->a1 < b2 + 20)
	    { 
	      if (sxw->a1 > a2 &&
		  ( sxw->eClass == 1 || /* intron */
		    sxw->eClass == 3 || /* SL */
		    sxw->eClass == 4  /* exon (G_U454 are not imported) */
		    )
		  )
		{ b2 = sxw->a1 - 20 ; break ; }
	      jj++ ; sxw++ ;
	    }  
	  while (jj > 0 && sxw->a1 > b2)
	    { jj-- ; sxw-- ;}
	  b2 = a2 ;
	  /* look back for a pA */
	  while (jj > 0 && sxw->a1 > a2)
	    { 
	      if (sxw->eClass == 2 && (sxw->type & 0x4 )) /* pA */
		{ b2 = sxw->a1 ; break ; } /* success */
	      jj-- ; sxw-- ;
	    }
	  if (b2 < a2 + 30)
	    goto failed ;
	}
      if (eClass == 7 && sxi->type & 0x8) /* up peche a la ligne a partir d'un stop */
	{
	  /* we do not care about the wiggle */
	  a2 = b2 = sxi->a2 ; a1 = sxi->a1 ; b1 = a1 - 500 ;
	  /* find the position of the stop */
	  while (jj < jMax && (sxw->a1 < a1 || sxw->type))
	    { jj++ ; sxw++ ;}
	  while (jj > 0 && sxw->a1 > b1)
	    { 
	      if (sxw->a1 < b1 + 20 &&
		  (sxw->eClass == 1 || /* intron */
		   sxw->eClass == 3 || /* SL */
		   sxw->eClass == 4  /* exon (G_U454 are not imported) */
		   )
		  )
		{ b1 = sxw->a1 + 20 ; break ; }
	      jj-- ; sxw-- ;
	    }  
	  while (jj < jMax && sxw->a1 < b1)
	    { jj++ ; sxw++ ;}
	  b1 = a1 ;
	  /* look back for a pA */
	  while (jj < jMax && sxw->a1 < a1)
	    { 
	      if (sxw->eClass == 2 && (sxw->type & 0x8 )) /* pA */
		{ b1 = sxw->a1 ; break ; } /* success */
	      jj++ ; sxw++ ;
	    }
	  if (b1 > a1 - 30)
	    goto failed ;
	}
      
      if (eClass == 1 && szwa && szwb && szwa < szwb)
	for (iLevel = aLevel ; szwa <= szwb ; szwa++)
	  if (szwa->m < iLevel && ! szwa->type)
	    iLevel = szwa->m ;

      if (eClass != 4 || b2 - b1 > 22)
	{
	  /* success: the boundaries are exon(a2,a1) intron exon(b1,b2) */
	  nn++ ;
	  
	  /* register the exonic parts as new exon */
	  sxk = arrayp (newElements, arrayMax (newElements), SXX) ;
	  sxk->type = sxi->type ;
	  sxk->eClass = eClass ;
	  sxk->a1 = a1 ; sxk->a2 = a2 ;
	  sxk->b1 = b1 ; sxk->b2 = b2 ;
	  sxk->aLevel = aLevel ; sxk->iLevel = iLevel ; sxk->bLevel = bLevel ;
	  sxk->nDonorSupport = sxi->nDonorSupport ;
	  sxk->nAcceptorSupport = sxi->nAcceptorSupport ;
	}
      if (eClass == 4 && b01 > 0 && b02 > 0 && (b01 < b1 - 10 || b02 > b2 + 10))
	{
	  if (b01 >= b1 - 10) b01 = b1 ; 
	  if (b02 <= b2 + 10) b02 = b2 ; 
	  /* register also the longest piece */
	  sxk = arrayp (newElements, arrayMax (newElements), SXX) ;
	  sxk->type = sxi->type ;
	  sxk->eClass = eClass ;
	  sxk->a1 = a1 ; sxk->a2 = a2 ;
	  sxk->b1 = b01 ; sxk->b2 = b02 ;
	  sxk->aLevel = aLevel ; sxk->iLevel = iLevel ; sxk->bLevel = bLevel ;
	  sxk->nDonorSupport = sxi->nDonorSupport ;
	  sxk->nAcceptorSupport = sxi->nAcceptorSupport ;
 	}
      continue ;
    failed:
      if (eClass == 2 || eClass == 22)
	{
	  /* register the failure of the pA */
	  sxk = arrayp (newElements, arrayMax (newElements), SXX) ;
	  *sxk = *sxi ;
	  sxk->eClass = 999 ; /* paFailed */
	}
    }
  
  ac_free (h) ;
  return nn ;
} /* sxNewExonsCreate */

/*************************************************************************************/
/* export a set of walls when we change strand with hysteresis
 * if we enter a strand we stay there untill contradicted
 */
static Array sxNewExonsGetStrandedWalls (SX *sx, Array w_f, Array w_r, AC_HANDLE h)
{
  Array walls ;
  int new, state1 = 0, state8 = 0, state, ii, jj, olda1 = 0 ;
  int level = 8 ;
  SXW *sxw ;
  SXX *sxx, *syy ;

  walls = arrayHandleCreate (1000, SXX, h) ;

  /* locate on plus strand the zone above or below theshold */
  state1 = state8 = 0 ;
  for (ii = jj = 0, sxw = arrp (w_f, 0, SXW), state = 0 ; ii < arrayMax (w_f) ; ii++, sxw++)
    {
      new = sxw->n ;
      if (! state8 && new >= level)
	{
	  sxx = arrayp (walls, jj++, SXX) ; /* create a wall at position zero */
	  sxx->n = 10 ; /* entering a plus region */
	  sxx->a1 = sxx->a2 = sxw->a1 - 1 ;
	  state8 = 1 ;
	}
      else if (state8 && new < level)
	{
	  sxx = arrayp (walls, jj++, SXX) ; /* create a wall at position zero */
	  sxx->n = -10 ; /* exiting a plus region */
	  sxx->a1 = sxx->a2 = sxw->a1 - 1 ;
	  state8 = 0 ;
	}
      if (! state1 && new >= 1)
	{
	  sxx = arrayp (walls, jj++, SXX) ; /* create a wall at position zero */
	  sxx->n = 1 ; /* entering a plus region */
	  sxx->a1 = sxx->a2 = sxw->a1 - 1 ;
	  state1 = 1 ;
	}
      else if (state1 && new < 1)
	{
	  sxx = arrayp (walls, jj++, SXX) ; /* create a wall at position zero */
	  sxx->n = -1 ; /* exiting a plus region */
	  sxx->a1 = sxx->a2 = sxw->a1 - 1 ;
	  state1 = 0 ;
	}
      olda1 = sxw->a1 ;
    }
  if (state8)
    {
      sxx = arrayp (walls, jj++, SXX) ; /* create a wall at position zero */
      sxx->n = -10 ; /* exiting a plus region */
      sxx->a1 = sxx->a2 = olda1 + 9 ;
      state1 = 0 ;
    }
  if (state1)
    {
      sxx = arrayp (walls, jj++, SXX) ; /* create a wall at position zero */
      sxx->n = -1 ; /* exiting a plus region */
      sxx->a1 = sxx->a2 = olda1 + 9 ;
      state1 = 0 ;
    }
    
  /* locate on minus strand the zone above or below theshold */
  state1 = state8 = 0 ;
  for (ii = 0, sxw = arrp (w_r, 0, SXW), state = 0 ; ii < arrayMax (w_r) ; ii++, sxw++)
    {
      new = sxw->n ;
      if (! state8 && new >= level)
	{
	  sxx = arrayp (walls, jj++, SXX) ; /* create a wall at position zero */
	  sxx->n = 20 ; /* entering a minus region */
	  sxx->a1 = sxx->a2 = sxw->a1 - 1 ;
	  state8 = 1 ;
	}
      else if (state8 && new < level)
	{
	  sxx = arrayp (walls, jj++, SXX) ; /* create a wall at position zero */
	  sxx->n = -20 ; /* exiting a minus region */
	  sxx->a1 = sxx->a2 = sxw->a1 - 1 ;
	  state8 = 0 ;
	}
      if (! state1 && new >= 1)
	{
	  sxx = arrayp (walls, jj++, SXX) ; /* create a wall at position zero */
	  sxx->n = 2 ; /* entering a minus region */
	  sxx->a1 = sxx->a2 = sxw->a1 - 1 ;
	  state1 = 1 ;
	}
      else if (state1 && new < 1)
	{
	  sxx = arrayp (walls, jj++, SXX) ; /* create a wall at position zero */
	  sxx->n = -2 ; /* exiting a minus region */
	  sxx->a1 = sxx->a2 = sxw->a1 - 1 ;
	  state1 = 0 ;
	}
      olda1 = sxw->a1 ;
    }
  if (state8)
    {
      sxx = arrayp (walls, jj++, SXX) ; /* create a wall at position zero */
      sxx->n = -20 ; /* exiting a plus region */
      sxx->a1 = sxx->a2 = olda1 + 9 ;
      state1 = 0 ;
    }
  if (state1)
    {
      sxx = arrayp (walls, jj++, SXX) ; /* create a wall at position zero */
      sxx->n = -2 ; /* exiting a plus region */
      sxx->a1 = sxx->a2 = olda1 + 9 ;
      state1 = 0 ;
    }
  arraySort (walls, sxxOrder) ;

  /* when we exit a +8 region, we extend it as long as we do not hit a negative or we are positive > 0 */
  for (ii = 0, sxx = arrp (walls, 1, SXX) ; ii < arrayMax (walls)  ; ii++, sxx++)
    {
      if (sxx->n == 1) state |= 1 ;
      if (sxx->n == 2) state |= 2 ;
      if (sxx->n == -1) state &= 0x2 ;
      if (sxx->n == -2) state &= 0x1 ;
      
      sxx->a2 = 0 ;
      /* extend to the right the positive region we are exiting until contradicted */
      if (sxx->n == -10)
	for (jj = ii, state1 = state & 0x2, syy = arrp (walls, jj, SXX) ; jj < arrayMax (walls)  ; jj++, syy++)
	  {
	    if (syy->n == -1) /* no more support for plus strand */
	      { sxx->a2 = syy->a1 ; if (state1)  break ; }
	    if (syy->n == 1)  /* again supported */
	      sxx->a2 = 0 ; 
	    if (syy->n == 2) /* entering a negative zone */
	      {
		state1 = 2 ;
		if (sxx->a2 > 0) { sxx->a2 = syy->a1 ; break ; }
	      }
	    if (syy->n == -2) /* exiting a negative zone */
	      state1 = 0 ;
	  }
      /* extend to the right the positive region we are exiting until contradicted */
      if (sxx->n == -20)
	for (jj = ii, state1 = state & 0x1, syy = arrp (walls, jj, SXX) ; jj < arrayMax (walls)  ; jj++, syy++)
	  {
	    if (syy->n == -2) /* no more support for minus strand */
	      { sxx->a2 = syy->a1 ; if (state1)  break ; }
	    if (syy->n == 2)  /* again supported */ 
	      sxx->a2 = 0 ; 
	    if (syy->n == 1) /* entering a positive zone */
	      {
		state1 = 2 ;
		if (sxx->a2 > 0) { sxx->a2 = syy->a1 ; break ; }
	      }
	    if (syy->n == -1) /* exiting a positive zone */
	      state1 = 0 ;
	  }

      /* extend left, notice that entering from the right is like exiting from the left */
      /* extend to the left the positive region we are entering until contradicted */
      if (sxx->n == 10)
	for (jj = ii, state1 = state & 0x2, syy = arrp (walls, jj, SXX) ; jj >= 0 ; jj--, syy--)
	  {
	    if (syy->n == 1) /* no more support for plus strand */
	      { sxx->a2 = syy->a1 ; if (state1)  break ; }
	    if (syy->n == -1)  /* again supported */
	      sxx->a2 = 0 ; 
	    if (syy->n == -2) /* entering a negative zone */
	      {
		state1 = 2 ;
		if (sxx->a2 > 0) { sxx->a2 = syy->a1 ; break ; }
	      }
	    if (syy->n == 2) /* exiting a negative zone */
	      state1 = 0 ;
	  }
      /* extend to the left the positive region we are entering until contradicted */
      if (sxx->n == 20)
	for (jj = ii, state1 = state & 0x1, syy = arrp (walls, jj, SXX) ; jj >= 0 ; jj--, syy--)
	  {
	    if (syy->n == 2) /* no more support for minus strand */
	      { sxx->a2 = syy->a1 ; if (state1)  break ; }
	    if (syy->n == -2)  /* again supported */ 
	      sxx->a2 = 0 ; 
	    if (syy->n == -1) /* entering a positive zone */
	      {
		state1 = 2 ;
		if (sxx->a2 > 0) { sxx->a2 = syy->a1 ; break ; }
	      }
	    if (syy->n == 1) /* exiting a positive zone */
	      state1 = 0 ;
	  }
    }
                
  /* remove the micro walls */   
   for (ii = jj = 0, sxx = syy = arrp (walls, 0, SXX) ; ii < arrayMax (walls)  ; ii++, sxx++)
    {
      if (sxx->n > 3 || sxx->n < -3)
	{
	  if (jj < ii) *syy = *sxx ;
	  jj++, syy++ ;
	}
    }
  arrayMax (walls) = jj ;

  /* after removing the +1 walls
   * if in between 2 +8 region on plus strand, there is no +8 on opposite
   * make a single region 
   */
  for (ii = 0, sxx = arrp (walls, 1, SXX) ; ii < arrayMax (walls)  ; ii++, sxx++)
    {
      /* extend to the right the positive region we are exiting until contradicted */
      if (sxx->n == -10)
	{
	  if ((sxx+1)->n == 10)
	    sxx->n = (sxx+1)->n = 0 ;
	}
      if (sxx->n == -20)
	{
	  if ((sxx+1)->n == 20)
	    sxx->n = (sxx+1)->n = 0 ;
	}
    }

  /* keep happy few */
   for (ii = jj = 0, sxx = syy = arrp (walls, 0, SXX) ; ii < arrayMax (walls)  ; ii++, sxx++)
    {
      if (sxx->n)
	{
	  if (jj < ii) 
	    *syy = *sxx ;
	  jj++; syy++ ;
	}
    }
  arrayMax (walls) = jj ;


  /* we are now left with setting a wall at the end of the tails of the remaining walls */
  /* the first wall and the last wall do not count */ 
  if (arrayMax (walls))
    {
      for (ii = 0, sxx = arrp (walls, 0, SXX) ; ii < 1 && ii < arrayMax (walls)  ; ii++, sxx++)
	sxx->n = 0 ;
      for (ii = arrayMax (walls) - 1, sxx = arrp (walls, ii, SXX) ; ii < arrayMax (walls)  ; ii++, sxx++)
	sxx->n = 0 ;
    }
  /* keep happy few */
  for (ii = jj = 0, sxx = syy = arrp (walls, 0, SXX) ; ii < arrayMax (walls)  ; ii++, sxx++)
    {
      if (sxx->n)
	{
	  if (jj < ii) *syy = *sxx ;
	  jj++ ; syy++ ;
	}
    }
  arrayMax (walls) = jj ;
  
  /* exchange the end walls */
  for (ii = 0, sxx = arrp (walls, 0, SXX), syy = sxx+1 ; ii + 1 < arrayMax (walls)  ; ii++, sxx++, syy++)
    {
      switch (sxx->n)
	{
	case -10:
	case -20:
	  if (syy->a2 > sxx->a2)
	   { int a0 = sxx->a2 ; sxx->a2 = syy->a2 ; syy->a2 = a0 ; }
	  break ;
	}
    }
  /* indicate the strand */
  for (ii = 0, sxx = arrp (walls, 0, SXX) ; ii < arrayMax (walls)  ; ii++, sxx++)
    {
      switch (sxx->n)
	{
	case 10:
	case -10:
	  sxx->isDown = TRUE ;
	  break ;
	default:
	  sxx->isDown = FALSE ;
	  break ;
	}
    }
  /* the remaining walls are significant, we need to break the exons on sxx->a2 of those */
  if (0)
    for (ii = 0, sxx = arrp (walls, 0, SXX) ; ii < arrayMax (walls)  ; ii++, sxx++)
      freeOutf ("WALL\t%s\t%s\t%09d\t%09d\t%d\n"
		, sxx->isDown ? "Forward" : "Reverse"
		, sx->sxxChromosome, sxx->a1, sxx->a2 , sxx->n
		) ;

  /* flip the coordinates before sorting */
  for (ii = 0, sxx = arrp (walls, 0, SXX) ; ii < arrayMax (walls)  ; ii++, sxx++)
    {
      int a1 = sxx->a1 ; sxx->a1 = sxx->a2 ; sxx->a2 = a1 ;
    }
  arraySort (walls, sxxOrder) ;
  arrayCompress (walls) ;


  return walls ;
} /* sxNewExonsGetStrandedWalls */

/*************************************************************************************/
/* use the stranded wiggles to select the strand of the new exons 
 *
 */
static int sxNewExonsSelectStrand (SX *sx, Array knowElements, Array newElements, Array w_f, Array w_r)
{
  AC_HANDLE h = ac_new_handle () ;
  int a1, a2, b1 = 0, b2, ii, jj, jj1, iir = 0, iif = 0, iifMax = 0, iirMax = 0, wMax, nr, nf, oldnf ;
  int xnf = 0, xnr = 0, xnf1 = 0, xnr1 = 0 ;
  SXX *sxx, *syy, *sxx1 ;
  SXW *sxw ;
  Array walls = 0 ;

  iifMax = w_f ? arrayMax (w_f) : 0 ;
  iirMax = w_r ? arrayMax (w_r) : 0 ;
  if (iifMax && iirMax && arrayMax (newElements))
    {
      walls = sxNewExonsGetStrandedWalls (sx, w_f, w_r, h) ;
      wMax = arrayMax (walls) ;
      arraySort (newElements, sxxOrder) ;

      /* cut the new exons in several new exons */
      jj = jj1 = arrayMax (newElements) ;
      sxx = arrayp (newElements, 2*jj, SXX) ; /* make room */
      for (ii = iif = 0  ; wMax && ii < jj ; ii++)
	{
	  sxx = arrayp (newElements, ii, SXX) ;
	  if (sxx->eClass != 4) /* exons */
	    continue ;
	  if (sxx->type & 0x8) sxx->isDown = FALSE ;
	  if (sxx->type & 0x4) sxx->isDown = TRUE ;
	  b1 = sxx->b1 ; b2 = sxx->b2 ;
	  /* search for a wall in the b1 b2 segment */
	  /* bactrack */
	  if (iif == iifMax) iif-- ;
	  for (syy= arrp (walls, iif, SXX) ; iif > 0 && syy->a1 > b1 ; syy--, iif--) ; 
	  /* locate */
	  for ( ; iif < wMax && syy->a1 < b1 ; syy++, iif++) ; 
	  for ( ; iif < wMax && syy->a1 < b2 ; syy++, iif++) 
	    {
	      if (sxx->isDown != syy->isDown) continue ;
	      sxx1 = arrayp (newElements, jj1++, SXX) ; /* may relocate newElements */
	      sxx = arrayp (newElements, ii, SXX) ;  /* restore sxx */
	      *sxx1 = *sxx ;
	      sxx->b1 = syy->a1 ; 
	      sxx1->b2 = syy->a1 ;
	    }
	}
      arrayMax (newElements) = jj1 ;     

      /* cut the other types (introns, SL pA and drop the left overs */
      for (ii = iif = 0  ; wMax && ii < jj ; ii++)
	{
	  sxx = arrayp (newElements, ii, SXX) ;
	  if (sxx->eClass == 4) /* exons */
	    continue ;
	  if (sxx->type & 0x8) sxx->isDown = FALSE ;
	  if (sxx->type & 0x4) sxx->isDown = TRUE ;

	  /* cut the forward a1 side */
	  if (sxx->b1 + 10 < sxx->a1)
	    { 
	      a1 = sxx->b1 ; a2 = sxx->a1 ;
	      /* search for a wall in the b1 b2 segment */
	      if (a2 > a1 + 10)
		{
		  /* bactrack */
		  if (iif == iifMax) iif-- ;
		  for (syy= arrp (walls, iif, SXX) ; iif > 0 && syy->a1 > a2 ; syy--, iif--) ; 
		  /* locate */
		  for ( ; iif < wMax && syy->a1 < a2 ; syy++, iif++) ; 
		  if (iif) { iif--; syy-- ; }
		  if (syy->a1 > a1 && syy->a1 + 10 < a2)
		    sxx->b1 = syy->a1 ;
		}
	    }
	  /* cut the reverse a1 side */
	  if (sxx->b1 - 10 > sxx->a1)
	    { 
	      a1 = sxx->a1 ; a2 = sxx->b1 ;
	      /* search for a wall in the b1 b2 segment */
	      if (a2 > a1 + 10)
		{
		  /* bactrack */
		  if (iif == iifMax) iif-- ;
		  for (syy= arrp (walls, iif, SXX) ; iif > 0 && syy->a1 > a1 ; syy--, iif--) ; 
		  /* locate */
		  for ( ; iif < wMax && syy->a1 < a1 ; syy++, iif++) ; 
		  if (syy->a1 > a1 + 10 && syy->a1 < a2)
		    sxx->b1 = syy->a1 ;
		}
	    }

	  /* cut forward the b1 side */  
	  if (sxx->a2 + 10 <  sxx->b2)
	    b1 = sxx->a2 ; 
	  b2 = sxx->b2 ;
	  /* search for a wall in the b1 b2 segment */
	  if (b2 > b1 + 10)
	    {
	      /* bactrack */
	      if (iif == iifMax) iif-- ;
	      for (syy= arrp (walls, iif, SXX) ; iif > 0 && syy->a1 > b1 ; syy--, iif--) ; 
	      /* locate */
	      for ( ; iif < wMax && syy->a1 < b1 ; syy++, iif++) ; 
	      if (syy->a1 > b1 + 10 && syy->a1 < b2)
		sxx->b2 = syy->a1 ;
	    }
	  /* cut reverse the b2 side */  
	  if (sxx->b2 - 10 >  sxx->a2)
	    b1 = sxx->b2 ; 
	  b2 = sxx->a2 ;
	  /* search for a wall in the b1 b2 segment */
	  if (b2 > b1 + 10)
	    {
	      /* bactrack */
	      if (iif == iifMax) iif-- ;
	      for (syy= arrp (walls, iif, SXX) ; iif > 0 && syy->a1 > b1 ; syy--, iif--) ; 
	      /* locate */
	      for ( ; iif < wMax && syy->a1 < b1 ; syy++, iif++) ; 
	      if (iif) { iif-- ; syy-- ; }
	      if (syy->a1 > b1 + 10 && syy->a1 < b2)
		sxx->b2 = syy->a1 ;
	    }
	}

      /* select the strand */

      for (ii = 0, sxx = arrayp (newElements, 0, SXX), oldnf = 0  ; ii < arrayMax (newElements) ; sxx++, ii++)
	{
	  if (sxx->type & 0x4)
	    sxx->isDown = TRUE ;
	  if (sxx->eClass != 4) /* exons */
	    continue ;
	  b1 = sxx->b1 ; b2 = sxx->b2 ;
	  /* count hits in the forward track */
	  /* bactrack */
	  if (iif == iifMax) iif-- ;
	  for (sxw = arrp (w_f, iif, SXW) ; iif > 0 && sxw->a1 > b1 ; sxw--, iif--) ; 
	  /* count */
	  for ( ; iif < iifMax && sxw->a1 < b1 ; sxw++, iif++) ; 
	  for (nf = 0 ; iif < iifMax && sxw->a1 < b2 ; sxw++, iif++) 
	    nf += sxw->n ; 

	  /* count hits in the reverse track */
	  /* bactrack */
	  if (iir == iirMax) iir-- ;
	  for (sxw = arrp (w_r, iir, SXW) ; iir > 0 && sxw->a1 > b1 ; sxw--, iir--) ; 
	  /* count */
	  for ( ; iir < iirMax && sxw->a1 < b1 ; sxw++, iir++) ; 
	  for (nr = 0 ; iir < iirMax && sxw->a1 < b2 ; sxw++, iir++) 
	    nr += sxw->n ;
	  if (0 && b2 < 8500) fprintf (stderr, "// %d:%d    nf:nr = %d:%d\n", b1, b2, nf, nr) ;
	  /* select best, by evidence or by inertia */
	  {
	    if (nf+5 > 2 * (nr + 5) || (nf > nr && oldnf == 1) || (!nf  && !nr && oldnf == 1))
	      { 
		if (sxx->type & 0x8) xnf1++ ;
		sxx->type &= ~0x8 ; sxx->type |= 0x4 ; xnf++ ; sxx->isDown = TRUE ;
		if (sx->antistrand) sxx->type = 0 ;
		oldnf = 1 ;
	      }
	    if (nr+5 > 2 * (nf + 5) || (nr > nf && oldnf == 2) || (!nf  && !nr && oldnf == 2))
	      {
		if (sxx->type & 0x4) xnr1++ ;  
		sxx->type &= ~0x4 ; sxx->type |= 0x8 ; xnr++ ; sxx->isDown = FALSE ;
		if (sx->strand) sxx->type = 0 ;	
		oldnf = 2 ;
	      }
	  }
	}
      fprintf (stderr, "// oriented %d:%d f:r exons, modifed %d:%d\n", xnf, xnr, xnf1, xnr1) ;
      arraySort (newElements, sxxOrder) ;
      arrayCompress (newElements) ;
    }
  ac_free (h) ;
  return xnf + xnr ;
} /* sxNewExonsSelectStrand */
      
/*************************************************************************************/
/* remove pA amplified doubles and pA in CDS
 *
 */
static int sxNewExonsCleanUp (SX *sx, Array knowElements, Array newElements)
{
  int nn = 0 ;
  int ii, jj, iMax = arrayMax (newElements) ;
  SXX *sxx, *sxy ;

  /* remove the new pA included in another one
   */
  arraySort (newElements, sxxPaOrder) ;
  for (ii = 0, sxx = arrayp (newElements, 0, SXX)  ; ii < iMax ; sxx++, ii++)
    {
      if (sxx->eClass != 2 || !(sxx->type & 0x4))
	continue ;
      
      for (jj = ii + 1, sxy = sxx + 1 ; sxx->eClass && jj < iMax ; jj++, sxy++)
	{
	  if (sxy->eClass != sxx->eClass || sxy->type != sxx->type)
	    continue ;
	  if (sxy->b1 > sxx->b1)
	    break ;
	  if (sxy->a1 > sxx->a1)
	    { nn++ ; sxx->eClass = 100 ; } /* pAincluded */
	  }
    }

  for (ii = 0, sxx = arrayp (newElements, 0, SXX)  ; ii < iMax ; sxx++, ii++)
    {
      if (sxx->eClass != 2 || !(sxx->type & 0x8))
	continue ;
      
      for (jj = ii + 1, sxy = sxx + 1 ; sxx->eClass && jj < iMax ; jj++, sxy++)
	{
	  if (sxy->eClass != sxx->eClass || sxy->type != sxx->type)
	    continue ;
	  if (sxy->a2 > sxx->b2)
	    break ;
	  if (sxy->b2 == sxx->b2)
	    { nn++ ; sxy->eClass = 100 ; }
	}
    }
  return nn ;
} /* sxNewExonsCleanUp */

/*************************************************************************************/
/* export the median of each segment making up the exon, comma delimited */
static int sxNewExonsExportExonSupport (SX *sx, Array segs, int a1, int a2)
{
  int nn = 0, u1, u2, jMax = arrayMax (segs) ;
  SXX *sxw ;
  char sep = '\t' ;
  static int jj = 0 ;

  if (a1 < a2) 
    { 
      sxw = arrp (segs, jj, SXX) ;
      while (jj > 0 && (sxw->a2 > a1 || sxw->type)) /* ignore artificial segments */ 
	{ jj-- ; sxw-- ; } 
      while (jj < jMax && (sxw->a2 < a1 || sxw->type)) /* ignore artificial segments */ 
	{ jj++ ; sxw++ ; } 
      while (jj < jMax && (sxw->a1 < a2 || sxw->type)) /* ignore artificial segments */ 
	{ 
	 if (!sxw->type)
	   {
	     u1 = sxw->a1 > a1 ? sxw->a1 : a1 ; 
	     u2 = sxw->a2 < a2 ? sxw->a2 : a2 ;
	     freeOutf ("%c%d:%d:%d", sep, u1 - a1 + 1, u2 - a1 + 1, sxw->m) ;
	     sep = ',' ;
	   }
	  jj++ ; sxw++ ;
	} 
    } 
  if (a1 > a2) 
    { 
      sxw = arrp (segs, jj, SXX) ;
      while (jj < jMax && (sxw->a1 < a1 || sxw->type)) /* ignore artificial segments */ 
	{ jj++ ; sxw++ ; } 
      while (jj > 0 && (sxw->a1 > a1 || sxw->type)) /* ignore artificial segments */ 
	{ jj-- ; sxw-- ; } 
      while (jj < jMax && (sxw->a2 > a2 || sxw->type)) /* ignore artificial segments */ 
	{  
	  if (!sxw->type)
	    {
	      u1 = sxw->a1 > a2 ? sxw->a1 : a2 ; 
	      u2 = sxw->a2 < a1 ? sxw->a2 : a1 ;
	      freeOutf ("%c%d:%d:%d", sep, a1 - u2 + 1, a1 - u1 + 1, sxw->m) ;
	      sep = ',' ;
	    }
	  jj-- ; sxw-- ; 
	}
    } 

  return nn ;
} /* sxNewExonsExportExonSupport */

 /*************************************************************************************/
/* export the results */
static int sxNewExonsExportDna (SX *sx, int a1, int a2)
{
  int nn = 0, u1, u2, dnaMax ;
  Array dna ;

 if (a1 < a2)
    {
      dna = sx->chromDnaD ; 
      dnaMax = dna ? arrayMax (dna) : 0 ;
      u1 = a1 ; u2 = a2 ; 
    }
  else
    {
      dna = sx->chromDnaR ;
      dnaMax = dna ? arrayMax (dna) : 0 ;
      u1 = dnaMax - a1 + 1 ;
      u2 = dnaMax - a2 + 1 ;
    }
  if (dna && u1 >= 1 && u2 <= dnaMax)
    {
      char *cp, cc ;
      
      nn = u2 - u1 + 1 ;
      cp = arrp (dna, u2, char) ;
      cc = *cp ; *cp = 0 ;
      freeOutf ("\t%s", arrp (dna, u1 - 1, char)) ;
      *cp = cc ;
    }
  return nn ;
} /* sxNewExonsExportDna */

/*************************************************************************************/
/* export the results */
static int sxNewExonsExport (SX *sx, Array segs, Array newElements)
{
  int ii ;
  SXX *sxk ;

  sxGetTargetFasta (sx) ;
  if (! sx->sxxChromosome)
    sx->sxxChromosome = "ANY" ;
  /* The strand is the strand of the gene
   * the coordinates are the coordinates of the DST (deep composite tag)
   * the DST for SL and introns are to be interpreted as Forward reads, the pA as Reverse reads
   */
  for (ii = 0 ; ii < arrayMax (newElements) ; ii++)
    {
      sxk = arrayp (newElements, ii, SXX) ;
      if (! sxk->eClass)
	continue ;
      /* export */
      if (sxk->eClass == 2 && (sxk->type & 0x4))  /* down pA, with up stream extension */
	{
	  freeOutf ("pA\tForward\t%s\t%09d\t%09d\t%d\t%d\t%d"
		    , sx->sxxChromosome, sxk->a1, sxk->b1, sxk->nDonorSupport, sxk->aLevel, sxk->bLevel
		    ) ;
	  sxNewExonsExportDna (sx, sxk->a1, sxk->b1) ;
	  freeOutf ("\n") ;
	}
      if (sxk->eClass == 2 && (sxk->type & 0x8))  /* up pA, with down stream extension */
	{
	  freeOutf ("pA\tReverse\t%s\t%09d\t%09d\t%d\t%d\t%d"
		    , sx->sxxChromosome, sxk->a2, sxk->b2, sxk->nAcceptorSupport, sxk->bLevel, sxk->aLevel
		    ) ;
	  sxNewExonsExportDna (sx, sxk->a2, sxk->b2) ;
	  freeOutf ("\n") ;
	}
      
      if (sxk->eClass == 100 && (sxk->type & 0x4))  /* down pA, with up stream extension */
	freeOutf ("paIncluded\tForward\t%s\t%09d\t%09d\t%d\n"
		  , sx->sxxChromosome, sxk->a1, sxk->b1, sxk->nDonorSupport
		  ) ;
      if (sxk->eClass == 100 && (sxk->type & 0x8))  /* up pA, with down stream extension */
	freeOutf ("paIncluded\tReverse\t%s\t%09d\t%09d\t%d\n"
		  , sx->sxxChromosome, sxk->a2, sxk->b2, sxk->nAcceptorSupport
		  ) ;
      
      if (sxk->eClass == 101 && (sxk->type & 0x4))  /* down pA, with up stream extension */
	freeOutf ("paInCDS\tForward\t%s\t%09d\t%09d\t%d\n"
		  , sx->sxxChromosome, sxk->a1, sxk->b1, sxk->nDonorSupport
		  ) ;
      if (sxk->eClass == 101 && (sxk->type & 0x8))  /* up pA, with down stream extension */
	freeOutf ("paInCDS\tReverse\t%s\t%09d\t%09d\t%d\n"
		  , sx->sxxChromosome, sxk->a2, sxk->b2, sxk->nAcceptorSupport
		  ) ;
            
      if (sxk->eClass == 102 && (sxk->type & 0x4))  /* down pA, with up stream extension */
	freeOutf ("paInExonSameStrand\tForward\t%s\t%09d\t%09d\t%d\n"
		  , sx->sxxChromosome, sxk->a1, sxk->b1, sxk->nDonorSupport
		  ) ;
      if (sxk->eClass == 102 && (sxk->type & 0x8))  /* up pA, with down stream extension */
	freeOutf ("paInExonSameStrand\tReverse\t%s\t%09d\t%09d\t%d\n"
		  , sx->sxxChromosome, sxk->a2, sxk->b2, sxk->nAcceptorSupport
		  ) ;
            
      if (sxk->eClass == 103 && (sxk->type & 0x4))  /* down pA, with up stream extension */
	freeOutf ("paInExonOppositeStrand\tForward\t%s\t%09d\t%09d\t%d\n"
		  , sx->sxxChromosome, sxk->a1, sxk->b1, sxk->nDonorSupport
		  ) ;
      if (sxk->eClass == 103 && (sxk->type & 0x8))  /* up pA, with down stream extension */
	freeOutf ("paInExonOppositeStrand\tReverse\t%s\t%09d\t%09d\t%d\n"
		  , sx->sxxChromosome, sxk->a2, sxk->b2, sxk->nAcceptorSupport
		  ) ;
            
      if (0 && sxk->eClass == 999 && (sxk->type & 0x4))  /* down pA, with up stream extension */
	freeOutf ("paFailed\tForward\t%s\t%09d\t%09d\t%d\n"
		  , sx->sxxChromosome, sxk->a1, sxk->b1, sxk->nDonorSupport
		  ) ;
      if (0 && sxk->eClass == 999 && (sxk->type & 0x8))  /* up pA, with down stream extension */
	freeOutf ("paFailed\tReverse\t%s\t%09d\t%09d\t%d\n"
		  , sx->sxxChromosome, sxk->a2, sxk->b2, sxk->nAcceptorSupport
		  ) ;
            
      if (sxk->eClass == 3 && (sxk->type & 0x4)) /* down SL, with down stream extension */
	{
	  freeOutf ("SL%d\tForward\t%s\t%09d\t%09d\t%d"
		    , sxk->nAcceptorSupport
		    , sx->sxxChromosome, sxk->a2 - 1, sxk->b2, sxk->nDonorSupport
		    ) ;
	  sxNewExonsExportDna (sx, sxk->a2 - 1, sxk->b2) ;
	  freeOutf ("\n") ;
	}
      if (sxk->eClass == 3 && (sxk->type & 0x8)) /* up SL, with up stream extension */
	{
	  freeOutf ("SL%d\tReverse\t%s\t%09d\t%09d\t%d"
		    , sxk->nAcceptorSupport
		    , sx->sxxChromosome, sxk->a1 + 1, sxk->b1, sxk->nDonorSupport
		    ) ;
	  sxNewExonsExportDna (sx, sxk->a1 + 1, sxk->b1) ;
	  freeOutf ("\n") ;
	}
      
      if (sxk->eClass == 1 && (sxk->type & 0x4) &&
	  1000 * (sxk->nDonorSupport + sxk->nAcceptorSupport) >  sxk->aLevel + sxk->bLevel
	  ) /* down introns */
	{
	  freeOutf ("INTRON\tForward\t%s\t%09d\t%09d\t%09d\t%09d\t%d\t%d\t%d\t%d\t%d"
		    , sx->sxxChromosome, sxk->b1, sxk->a1-1, sxk->a2+1, sxk->b2
		    , sxk->aLevel, sxk->nDonorSupport, sxk->iLevel, sxk->nAcceptorSupport, sxk->bLevel
		    ) ;
	  sxNewExonsExportDna (sx, sxk->b1, sxk->a1 - 1) ;
	  sxNewExonsExportDna (sx, sxk->a2 + 1, sxk->b2) ;
	  freeOutf ("\n") ;
	}
      if (sxk->eClass == 1 && (sxk->type & 0x8) &&
	  1000 * (sxk->nDonorSupport + sxk->nAcceptorSupport) >  sxk->aLevel + sxk->bLevel	  
	  ) /* up introns */
	{
	  freeOutf ("INTRON\tReverse\t%s\t%09d\t%09d\t%09d\t%09d\t%d\t%d\t%d\t%d\t%d"
		    , sx->sxxChromosome, sxk->b2, sxk->a2+1, sxk->a1-1, sxk->b1
		    , sxk->bLevel, sxk->nAcceptorSupport, sxk->iLevel, sxk->nDonorSupport, sxk->aLevel
		    ) ;
	  sxNewExonsExportDna (sx, sxk->b2, sxk->a2 + 1) ;
	  sxNewExonsExportDna (sx, sxk->a1 - 1, sxk->b1) ;
	  freeOutf ("\n") ;
	}
      
      if (sxk->eClass == 4 && (sxk->type & 0x4)) /* down new exons */
	{
	  freeOutf ("EXON\tForward\t%s\t%09d\t%09d"
		    , sx->sxxChromosome, sxk->b1, sxk->b2
		    ) ;
	  sxNewExonsExportExonSupport (sx, segs, sxk->b1, sxk->b2) ;
	  sxNewExonsExportDna (sx, sxk->b1, sxk->b2) ;
	  freeOutf ("\n") ;
	}
      
      if (sxk->eClass == 4 && (sxk->type & 0x8)) /* up new exons */
	{
	  freeOutf ("EXON\tReverse\t%s\t%09d\t%09d"
		    , sx->sxxChromosome, sxk->b2, sxk->b1
		    ) ;
	  sxNewExonsExportExonSupport (sx, segs, sxk->b2, sxk->b1) ;
	  sxNewExonsExportDna (sx, sxk->b2, sxk->b1) ;
	  freeOutf ("\n") ;
	}
      if (sxk->eClass == 7 && (sxk->type & 0x4)) /* down STOP */
	{
	  freeOutf ("STOP\tForward\t%s\t%09d\t%09d"
		    , sx->sxxChromosome, sxk->a1, sxk->b2
		    ) ;
	  freeOutf ("\n") ;
	}
      if (sxk->eClass == 7 && (sxk->type & 0x8)) /* up STOP */
	{
	  freeOutf ("STOP\tReverse\t%s\t%09d\t%09d"
		    , sx->sxxChromosome, sxk->a2, sxk->b1
		    ) ;
	  freeOutf ("\n") ;
	}
    }
  return 0 ;
} /* sxNewExonsExport */

/*************************************************************************************/
/* at each position, assume the stranded wiggle leaks up to a fraction 1/echo, and remove it from the other strand */
static void sxNewExonsStrandedWiggleRemoveEchoes (SX *sx, Array w_f, Array w_r)
{  
  float echo = 1.0/20.0 ; /* remove up to 5% echo */
  SXW *s1, *s2 ;
  int i, j ;
  int u, v, u1, v1 ;
  int iMax = arrayMax (w_f) ;
  int jMax = arrayMax (w_r) ;

  s1 = arrp (w_f, 0, SXW) ;
  s2 = arrp (w_r, 0, SXW) ; 
  i = j = 0 ;

  while (i < iMax && j < jMax)
    {
      if (s1->a1 == s2->a1)
	{ 
	  /* remove echo */
	  u = s1->m ; v = s2->n ;
	  u1 = u - v * echo ;
	  v1 = v - u * echo ;  
	  if (u1 < 0) {  s1->n = 0 ; s2->n += u ; } /* do not loose any signal, redistribute it */
	  if (v1 < 0) {  s2->n = 0 ; s1->n += v ; } 
	  if (u1 > 0 && v1 > 0)
	    {
	      s1->n = u1 + (u1/(u1+v1)) * echo * (u + v) ;
	      s2->n = v1 + (v1/(u1+v1)) * echo * (u + v) ;
	    }
	  s1++ ; s2++ ; i++ ; j++ ;
	}
      else if (s1->a1 < s2->a1)
	{
	  s1++ ; i++ ; 
	}
      else
	{
	  s2++ ; j++ ; 
	}
    }

  return ;
} /* sxNewExonsStrandedWiggleRemoveEchoes */

/*************************************************************************************/
int SOLEXA_STEP = 0 ; /* agnostic */
static int sxNewExons (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  int n, nn = 0 ;
  Array walls, segs ;
  Array wiggle = 0 ;
  Array knownElements = 0 ;
  Array candidates = 0 ;
  Array newElements = arrayHandleCreate (10000, SXX, h) ;
  BOOL debug = FALSE ;
  Array w_f = 0, w_r = 0  ;

  /* get the stranded wiggles */
   if (sx->wiggleFileName_f && sx->wiggleFileName_r)
    {
      w_f = sxNewExonsGetWiggle (sx, &SOLEXA_STEP, sx->wiggleFileName_f, h) ;
      w_r = sxNewExonsGetWiggle (sx, &SOLEXA_STEP, sx->wiggleFileName_r, h) ;
      if (w_f && w_r)
	sxNewExonsStrandedWiggleRemoveEchoes (sx, w_f, w_r) ; /* remove echoes */
    }


  if ((wiggle = sxNewExonsGetWiggle (sx, &SOLEXA_STEP, sx->wiggleFileName_ns, h)))        /* parse the wiggle file */
    {
      if (debug) fprintf (stderr, "get wiggle done %d\n", arrayMax (wiggle)) ;
      sxNewExonsRollingMedian (sx->delta ? sx->delta : 3, wiggle) ;         /* add the rolling median */
      if (w_f) sxNewExonsRollingMedian (5, w_f) ;  
      if (w_r) sxNewExonsRollingMedian (5, w_r) ;
      if (debug) fprintf (stderr, "get median done %d\n", arrayMax (wiggle)) ;
      
      candidates = sxNewExonsGetIntrons (sx, TRUE, h) ; /* get the new introns and avoid repeats */
      sxNewExonsMergeWeakPa (candidates) ;
	      
      if (debug) fprintf (stderr, "get candidates done %d\n", arrayMax (candidates)) ;
      knownElements = sxNewExonsGetIntrons (sx, FALSE, h) ; /* get the known introns */
      if (debug) fprintf (stderr, "get known elements done %d\n", arrayMax (knownElements)) ;
      if (1)
	{
	  n = sxNewExonsRemoveCDSpA (candidates, knownElements, newElements, FALSE) ;
	  n += sxNewExonsRemoveCDSpA (candidates, knownElements, newElements, TRUE) ;
	  if (1) fprintf (stderr, "eliminated %d pA in CDS\n", n) ;
	}
      walls = sxNewExonsGetWalls (sx, 0, wiggle, 0, debug, h) ;
      if (1) fprintf (stderr, "get walls done %d\n", arrayMax (walls)) ;
      sxNewExonsGetIntronsWalls (sx, walls, knownElements, candidates, h) ;
      if (1) fprintf (stderr, "introns walls done %d\n", arrayMax (walls)) ;
      walls = sxNewExonsGetWalls (sx, walls, wiggle, 1, debug, h) ;
      segs = sxNewExonsSegments (walls, h) ; 
      if (1) fprintf (stderr, "new exons segs done %d\n", arrayMax (segs)) ;
      sxNewExonsSegmentsSignal (wiggle, segs, debug) ; /* compute the median of each segment */
      if (debug) fprintf (stderr, "signal done %d\n", arrayMax (segs)) ;
      if (candidates)
	{
	  if (1) /* search for Introns */
	    {
	      n = 0 ;
	      if (! sx->antistrand)
		n += sxNewExonsCreate (sx, segs, candidates, newElements, 0x4, 1) ; /* down intron */
	      if (! sx->strand)
		n += sxNewExonsCreate (sx, segs, candidates, newElements, 0x8, 1) ; /* up intron */
	      nn += n ;
	      fprintf (stderr, "Constructed %d exon-exon junctions\n", n) ;
	      if (debug) 
		{ 
		  showSXX (newElements) ;
		}
	    }
	  if (1) /* search for pA */
	    {
	      n = 0 ;
	      if (! sx->antistrand)
		n += sxNewExonsCreate (sx, segs, candidates, newElements, 0x4, 2) ; /* down pA */
	      if (! sx->strand)
		n += sxNewExonsCreate (sx, segs, candidates, newElements, 0x8, 2) ; /* up pA */
	      nn += n ;
	      fprintf (stderr, "Constructed %d polyA extensions\n", n) ;
	      if (debug) 
		{
		  showSXX (newElements) ;
		}
	    }
	  if (1) /* search for SL */
	    {
	      n = 0 ;
	      if (! sx->antistrand)
		n += sxNewExonsCreate (sx, segs, candidates, newElements, 0x4, 3) ; /* down SL */
	      if (! sx->strand)
		n += sxNewExonsCreate (sx, segs, candidates, newElements, 0x8, 3) ; /* up SL */
	      nn += n ;
	      if (1) fprintf (stderr, "Constructed %d SL extensions\n", n) ;
	    }
	  if (1) /* search for new exons */
	    {
	      n = sxNewExonsAddNewExons (sx, segs, candidates, knownElements, newElements) ;
	      if (debug) fprintf (stderr, "Proposed %d exon candidates\n", n) ;
	      
	      n = 0 ;
	      if (! sx->antistrand)
		n += sxNewExonsCreate (sx, segs, candidates, newElements, 0x4, 4) ; /* down newExon */
	      if (! sx->strand)
		n += sxNewExonsCreate (sx, segs, candidates, newElements, 0x8, 4) ; /* up newExon */
	      nn += n ;
	      fprintf (stderr, "Constructed %d exons\n", n) ;
	      if (debug) 
		{
		  showSXX (newElements) ;
		}
	    }
	  if (0) /* search for Stop, MUST come last */
	    {
	      n = 0 ;
	      if (! sx->antistrand)
		n += sxNewExonsCreate (sx, segs, candidates, newElements, 0x4, 7) ; /* down Stop */
	      if (! sx->strand)
		n += sxNewExonsCreate (sx, segs, candidates, newElements, 0x8, 7) ; /* up Stop */
	      nn += n ;
	      if (1) fprintf (stderr, "Constructed %d Stop extensions\n", n) ;
	    }
	}

      sxNewExonsCleanUp (sx, knownElements, newElements) ;
      if (! sx->sxxChromosome)
	sx->sxxChromosome = "sxxChromX" ;
      sxNewExonsSelectStrand (sx, knownElements, newElements, w_f, w_r) ;
      sxNewExonsExport (sx, segs, newElements) ;
    }

  ac_free (h) ;
  return nn ;
} /* sxNewExons */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

typedef struct spongeStruct { int gene, exon, chrom, a1, a2, level ; long int ns, nbp ; BOOL isDown ; } SPONGE ;

static int sdfOrder (const void *a, const void *b)
{
  int n ;
  const SPONGE *x = (const SPONGE *) a, *y = (const SPONGE *) b ;

  n = x->level - y->level ; if (n) return n ;
  n = x->chrom - y->chrom  ; if (n) return n ;
  n = x->a1 - y->a1 ; if (n) return n ;
  n = x->a2 - y->a2 ; if (n) return n ;
  n = x->gene - y->gene ; if (n) return n ;
  n = x->exon - y->exon ; if (n) return n ;

  return 0 ;
} /* sdfOrder */

static int spongeOrder (const void *a, const void *b)
{
  int n ;
  const SPONGE *x = (const SPONGE *) a, *y = (const SPONGE *) b ;

  n = x->level - y->level ; if (n) return n ;
  n = x->chrom - y->chrom  ; if (n) return n ;
  n = x->gene - y->gene ; if (n) return n ;
  n = x->a1 - y->a1 ; if (n) return n ;
  n = x->a2 - y->a2 ; if (n) return n ;
  n = x->exon - y->exon ; if (n) return n ;

  return 0 ;
} /* spongeOrder */

static void showSponge (Array aa)
{
  int i ;
  SPONGE *s ;
  
  for (i = 0 ; i < arrayMax (aa) ; i++)
    {
      s = arrp (aa, i, SPONGE) ;
      fprintf (stderr, "%d:level_%d, %d\t%d\tc=%d:%d:%d\t%d\t%ld\t%ld\n"
	       , i
	       , s->level
	       , s->gene, s->exon, s->chrom
	       , s->a1, s->a2, s->a2 - s->a1
	       , s->ns, s->nbp
	       ) ;
    }
  if (0) showSponge (aa) ; /* for compiler happiness */
} /* showSponge */

/*************************************************************************************/

static void sxSpongeParseSplitMrnas (SX *sx, DICT *dict)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = 0 ;
  const char *ccp ;
  KEYSET ks = sx->splitMrnas = keySetHandleCreate (sx->h) ;
  const char *fNam = sx->splitMrnasFileName ;

  ai = aceInCreate (fNam, 0, h) ;
  if (! ai)
    messcrash ("Cannot find the split_mrnas File %s, sorry", fNam ? fNam : "NULL") ;
  else
    {
      char buf [1024] ;
      int a1 = 0, a2 = 0 ;

      aceInSpecial (ai, "\n\t") ;
      while (aceInCard (ai))
	{
	  int mrna = 0, gene = 0 ;

	  ccp = aceInWord(ai) ; /* old gene */
	  if (!ccp || !*ccp)
	    continue ;
	  strncpy (buf, ccp, 1023) ;
	  aceInStep(ai,'\t') ; ccp = aceInWord(ai) ; /* mrna */
	  if (!ccp || !*ccp)
	    continue ;	
	  dictAdd (dict, ccp, &mrna) ;

	  aceInStep(ai,'\t') ; 
	  if (!aceInInt (ai, &a1)) /* a1 */
	    continue ;	
	  aceInStep(ai,'\t') ; 
	  if (!aceInInt (ai, &a2)) /* a2 */
	    continue ;	

	  aceInStep(ai,'\t') ; ccp = aceInWord(ai) ; /* gene */
	  if (!ccp || !*ccp)
	    continue ;	
	  aceInStep(ai,'\t') ; ccp = aceInWord(ai) ; /* gene */
	  if (!ccp || !*ccp || !strcasecmp (ccp, buf))
	    continue ;	
	  aceInStep(ai,'\t') ; ccp = aceInWord(ai) ; /* gene(gene) */
	  if (!ccp || !*ccp)
	    continue ;	
	  
	  dictAdd (dict, ccp, &gene) ;
	  keySet (ks, mrna) = a1 > 0 ? (1<<30) : gene ; /* remove the split mrnas from the mapping game */
	}
    }
  ac_free (h) ;
  return ;
} /* sxSpongeParseSplitMrnas */

/*************************************************************************************/

static void sxSpongeParseOneFile (SX *sx, const char *fNam, DICT *dict, Array segs, KEYSET ks, int level, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = 0 ;
  int jj = 0, kk, a1, a2, gene, exon, chrom ;
  const char *ccp ;
  SPONGE *sxx ;

  jj = arrayMax (segs) ;
  kk = keySetMax (ks) ;

  ai = aceInCreate (fNam, 0, h) ;
  if (! ai)
    messerror ("Cannot find the spongeFile %s, sorry", fNam ? fNam : "NULL") ;
  else
    {
      aceInSpecial (ai, "\n\t") ;
      while (aceInCard (ai))
	{
	  int mrna, geneSplit = 0 ;

	  ccp = aceInWord(ai) ; /* gene or mrna*/
	  if (!ccp || !*ccp)
	    continue ;	
	  if (sx->splitMrnas)    /* else we do not care about the mrna */
	    {
	      dictAdd (dict, ccp, &mrna) ;
	      geneSplit = keySet (sx->splitMrnas, mrna) ;
	      if (geneSplit == (1 << 30))
		continue ;
	    }
	  aceInStep(ai,'\t') ; ccp = aceInWord(ai) ; /* exon mumber */
	  if (!ccp || !*ccp)
	    continue ;	
	  dictAdd (dict, ccp, &exon) ;
	  
	  aceInStep(ai,'\t') ; ccp = aceInWord(ai) ; /* chromosome */
	  if (!ccp || !*ccp)
	    continue ;
	  if (sx->sxxChromosome && strcmp (ccp, sx->sxxChromosome))	
	    continue ;
	  dictAdd (dict, ccp, &chrom) ;
	  
	  /* a1,a2 : chromosome coordinates */
	  aceInStep(ai,'\t') ; if (! aceInInt (ai, &a1) || a1 <1)
	    continue ; 
	  aceInStep(ai, '\t') ;
	  if (! aceInInt (ai, &a2) || a2<1)
	    continue ; 
	  if (sx->strand && a1 > a2) continue ;
	  if (sx->antistrand && a1 < a2) continue ;

          /* gene, either identical to column 1, or the gene of the mrna */
	  aceInStep(ai,'\t') ; ccp = aceInWord(ai) ; /* gene */
	  if (!ccp || !*ccp)
	    continue ;	
	  dictAdd (dict, ccp, &gene) ;
	  if (geneSplit)
	    gene = geneSplit ;
	  
	  sxx = arrayp (segs, jj++, SPONGE) ;
	  sxx->isDown = a1 < a2 ? TRUE : FALSE ;
	  if (1) /* oriented  along chromosome */
	    {
	      sxx->a1 = a1 <= a2 ? a1 : a2 ;
	      sxx->a2 = a1 > a2 ? a1 : a2 ;
	      keySet (ks, kk++) = a1 ;
	      keySet (ks, kk++) = a2 ;
	    }
	  sxx->gene = gene ;
	  sxx->exon = exon ;
	  sxx->chrom = chrom ;
	  sxx->level = level ;
	}
    }
  
  ac_free (h) ;

  return ;
} /* sxSpongeParseOneFile */

/*************************************************************************************/
/* remove duplicated areas */
static void sxPressSponge (Array segs)
{
  int ii, jj, iMax = arrayMax (segs) ;
  SPONGE *up, *vp ;

  for (ii = 0 ; ii < iMax ; ii++)
    {
      up = arrp (segs, ii, SPONGE) ;
      int level = up->level ;
      int chrom = up->chrom ;
      int gene = up->gene ;

      if (gene == -1)
	continue ;
      for (vp = up +1, jj = ii + 1 ;  jj < iMax && 
	     vp->gene == gene &&
	     vp->level == level &&
	     vp->chrom == chrom &&
	     vp->a1 <= up->a2 + 1 
	     ; jj++, vp++)
	{ 
	  if (up->a2 < vp->a2)
	    up->a2 = vp->a2 ;
	      vp->gene = -1 ;
	}
    }
  for (ii = jj = 0, up = vp = arrp (segs, ii, SPONGE) ; ii < iMax ; up++, ii++)
    {
      if (up->gene == -1)
	continue ;
      if (jj < ii)
	*vp = *up ;
      jj++ ; vp++ ;
    }
  arrayMax (segs) = jj ;
} /* sxPressSponge */

/*************************************************************************************/
/* return the true segments
 * report in segsNR the non redundant segments
 * temporarilly store in ks ALL the boundaries 
 */ 

static int sxSpongeParseFile (SX *sx, DICT *dict, DICT *levelNames, Array segs, Array segsNR, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  int jj = 0, kk = 0, a1, a2, level = 0 ;
  SPONGE *sxx ;
  KEYSET ks = keySetHandleCreate (h) ;

  if (! sx->spongeFileName)
    {
      keySet (ks, 0) = 0x7fffffff ;
    }
  else
    {   /* parse a hierarchic collection of sponge files */
      char *cp, *cq, *cr, *buf ;
      cp = buf = strnew (sx->spongeFileName, h) ;
      
      while (cp)
	{
	  while (*cp == ',') cp++ ;
	  cq = strstr (cp, ",") ;
	  if (cq) *cq = 0 ;
	  if (*cp)
	    {
	      level++ ;
	      cr = strstr (cp, ":") ;
	      if (cr)
		{
		  *cr = 0 ;
		  dictAdd (levelNames, cp, 0) ;
		  cp = cr + 1 ;
		}
	      else
		dictAdd (levelNames, hprintf (h, "level_%d", level), 0) ;
	      sxSpongeParseOneFile (sx, cp, dict, segs, ks, level, h0) ;
	    }
	  cp = cq ? cq + 1 : 0 ;
	}

      arraySort (segs, spongeOrder) ;
      arrayCompress (segs) ;
      sxPressSponge (segs) ;
      keySetSort (ks) ; 
      keySetCompress (ks) ;
    }

  /* NR is a non overlapping set of segments covering the whole zone */
  if (segsNR)
    for (jj = kk = a1 = 0 ; kk < keySetMax (ks) ; kk++)
      {
	a2 = keySet (ks, kk) ;
	if (a2 > a1)
	  {
	    sxx = arrayp (segsNR, jj++, SPONGE) ; 
	    sxx->a1 = a1 ; sxx->a2 = a2 - 1 ;
	  }
	sxx = arrayp (segsNR, jj++, SPONGE) ; 
	sxx->a1 = a2 ; sxx->a2 = a2  ; /* each boundary point is isolated */
	a1 = a2 + 1 ;
      }

  ac_free (h) ;
  return level ;
} /* sxSpongeParseFile */

/*************************************************************************************/

static long int sxSponge (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  int level, levelMax = 0, wiggleExtent = 0, chromExtent = 0 ;
  int a1, a2, ii, jj, jjMax, sponge = sx->sponge ;
  Array segs = arrayHandleCreate (100000, SPONGE, h) ;
  Array wiggle = 0 ;
  BOOL debug = FALSE ;
  SPONGE *seg, *seg2 ;
  SXW *sxw ;
  DICT *dict = dictHandleCreate (10000, h) ;
  DICT *levelNames = dictHandleCreate (10000, h) ;
  Array segsNR = arrayHandleCreate (100000, SPONGE, h) ;
  Array nss, nbps, nextent ;
  long int Nbp = 0 , Ns = 0 ;

  if (sx->splitMrnasFileName)
    sxSpongeParseSplitMrnas (sx, dict) ; 
  levelMax = sxSpongeParseFile (sx, dict, levelNames, segs, segsNR, h) ;
  jjMax = arrayMax (segsNR) ; 
  if (jjMax)
    {
      seg = arrp (segsNR, jjMax - 1, SPONGE) ;
      chromExtent = seg->a2 ;
    }
  if ((wiggle = sxNewExonsGetWiggle (sx, &SOLEXA_STEP, sx->wiggleFileName, h)) &&
      arrayMax (wiggle)
      )        /* parse the wiggle file and set SOLEXA_STEP */
    {
      /* add a segment for the whole wiggle */
      wiggleExtent = arr (wiggle, arrayMax(wiggle) - 1, SXW).a1 - arr (wiggle, 0, SXW).a1 ;

      if (debug) fprintf (stderr, "get wiggle done %u\n", arrayMax (wiggle)) ;
      for (ii = jj = 0, sxw = arrp (wiggle, 0, SXW), seg = arrp (segsNR, 0, SPONGE) ; ii < arrayMax (wiggle) ; ii++, sxw++)
	{
	  if (sxw->n < sponge) 
	    continue ;
	  
	  Ns += SOLEXA_STEP ;
	  Nbp += ((long int)sxw->n) * SOLEXA_STEP ;
	  a1 = sxw->a1 ;
	  while (jj < jjMax && seg->a2 < a1) 
	    { seg++ ; jj++ ; }
	  if (seg->a1 <= a1 && seg->a2 >= a1) 
	    {
	      seg->nbp += ((long int)sxw->n) * SOLEXA_STEP   ;
	      seg->ns += SOLEXA_STEP ; /* number of supported positions in this segement */
	    }	    
	}
    }
  else
    {
      if (sx->wiggleFileName)
	messcrash ("could not open wiggle file : %s", sx->wiggleFileName) ;
    }
  nss = arrayHandleCreate (2 + levelMax, long int, h) ;
  nbps = arrayHandleCreate (2 + levelMax, long int, h) ;
  nextent = arrayHandleCreate (2 + levelMax, long int, h) ;

  if (wiggleExtent > chromExtent)
    chromExtent = wiggleExtent ;
  /* fan out the NR seg in the original segs */
  for (ii = jj = 0, seg = arrp (segs, 0, SPONGE) ; ii < arrayMax (segs) ; seg++, ii++)
    {
      if (! seg->level) continue ;
      a1 = seg->a1 ; a2 = seg->a2 ;
      if (jj == jjMax) jj-- ;
      for (seg2 = arrp (segsNR, jj, SPONGE) ; jj > 0 && seg2->a1 > a1 ; seg2--, jj--) ;
      for (seg2 = arrp (segsNR, jj, SPONGE) ; jj < jjMax && seg2->a1 <= a2 ; seg2++, jj++)
	if (seg2->a1 >= a1 && seg2->a2 <= a2)
	  {
	    seg->nbp += seg2->nbp ; seg->ns += seg2->ns ; 
	    
	    if (! seg2->level) /* label the subsegment */
	      {
		seg2->level = seg->level ;
	      }
	  } /* number of supported positions in this segment */
    }


  /* first we impose a strict hiearchy
   * In this way we guarantee that CDS are included in exons, exons in genes1
   */
  for (jj = 0, seg2 = arrp (segsNR, 0, SPONGE) ; jj < jjMax ; seg2++, jj++)
    { 
      if (! seg2->level) continue ;
      for (level = levelMax ; level >= seg2->level ; level--)
	{
	  array (nextent, level, long int) +=  seg2->a2 - seg2->a1 + 1 ;
	  array (nbps, level, long int) += seg2->nbp ;
	  array (nss, level, long int) += seg2->ns ;
	}
    }
  /* fill the outside level */
  level = levelMax + 1 ;
  array (nextent, level, long int) = chromExtent ;
  array (nbps, level, long int) = Nbp ;
  array (nss, level, long int) = Ns ;

  if (debug)
    for (level = 1 ; level <= levelMax + 1 ; level++)
      fprintf (stderr, "%d: %ld\n", level, array (nextent, level, long int)) ;
  
  /* Then we substract the cumuls from each other */
  for (level = levelMax + 1 ; level > 0 ; level--)
    {
      arr (nextent, level, long int) -= arr (nextent, level - 1, long int) ;
      arr (nbps, level, long int) -= arr (nbps, level - 1, long int) ;
      arr (nss, level, long int) -= arr (nss, level - 1, long int) ;
    }
  if (debug)
    for (level = 1 ; level <= levelMax + 1 ; level++)
      fprintf (stderr, "%d: %ld\n", level, array (nextent, level, long int)) ;

  /* report */
  if (arrayMax (segs) == 0)
    {
      seg = arrayp (segsNR, 0, SPONGE) ;
      freeOutf ("# No spongeFile\n#Wiggle\tThreshold\tPositions\tbp\n%s\t%d\t%ld\t%ld\n"
		, sx->wiggleFileName, sponge, seg->ns, seg->nbp
		) ;
    }
  else
    {
      double tot = 0 ;
      for (level = levelMax ; level >= 2 ; level--)
        if (tot < arr (nbps, level, long int))  
	  tot = arr (nbps, level, long int) ;
      if (tot < 1) tot = 1 ;

      freeOutf("#Total\tLevel\tThreshold\tCategory\tExtent bp\tAbove threshold bp\tRNA bp\t%% extent supported\tFold coverage of supported part\tFold coverage of extent\t%% base in area\n", sponge) ;
      for (level = levelMax ; level >= 2 ; level--)
	freeOutf ("Total\t%d\t%d\t%s_not_%s\t%ld\t%ld\t%ld\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n"
		  , levelMax+1 - level
		  , sponge
		  , dictName (levelNames, level), dictName (levelNames, level - 1)
		  , arr (nextent, level, long int)
		  , arr (nss, level, long int)
		  , arr (nbps, level, long int)
		  , 100.0 * arr (nss, level, long int)/((double)1.0+arr (nextent, level, long int))
		  , (float)arr (nbps, level, long int)/((double)1.0+arr (nss    , level, long int))
		  , (float)arr (nbps, level, long int)/((double)1.0+arr (nextent, level, long int))
		  , 100.0 * ((double)arr (nbps, level, long int))/tot
		  ) ;
      /* cumulate again, except for the outside region */
      if (debug)
	for (level = 1 ; level <= levelMax + 1 ; level++)
	  fprintf (stderr, "%d: %ld\n", level, array (nextent, level, long int)) ;

      for (level = 2 ; level <= levelMax ; level++)
	{
	  arr (nextent, level, long int) += arr (nextent, level - 1, long int) ;
	  arr (nbps, level, long int) += arr (nbps, level - 1, long int) ;
	  arr (nss, level, long int) += arr (nss, level - 1, long int) ;
	}

      if (debug)
	for (level = 1 ; level <= levelMax + 1 ; level++)
	  fprintf (stderr, "%d: %ld\n", level, array (nextent, level, long int)) ;

      for (level = 1 ; level <= levelMax + 1  ; level++)
	freeOutf ("Total\t%d\t%d\t%s\t%ld\t%ld\t%ld\t%.1lf\t%.1lf\t%.1lf\t%.1lf\n"
		  , 10+level
		  , sponge
		  , level <= levelMax ? dictName (levelNames, level) : "Outside"
		  , arr (nextent, level, long int)
		  , arr (nss, level, long int)
		  , arr (nbps, level, long int)
		  , 100.0 * arr (nss, level, long int)/((double)1.0+arr (nextent, level, long int))
		  , (double)arr (nbps, level, long int)/((double)1.0+arr (nss, level, long int))
		  , (double)arr (nbps, level, long int)/((double)1.0+arr (nextent, level, long int)) 
		  , (double)arr (nbps, level, long int)/tot
		  ) ;

      freeOutf("#Level\tRun\tGene\tPseudo-exon\tChrom\tfrom\tto\tLength\tThreshold\tLength covered\tAligned bp\tFold coverage of supported part\n") ;   /* \tFold coverage of pseudo-exon */
      for (level = 1 ; level <= levelMax ; level++)
	for (jj = 0, seg = arrp (segs, 0, SPONGE) ; jj < arrayMax(segs) ; jj++, seg++)
	  {
	    if (seg->level != level || ! seg->ns) 
	      continue ;
	    
	    freeOutf ("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d"
		      , dictName (levelNames, level)
		      , sx->run ? sx->run : "-"
		      , dictName (dict, seg->gene)
		      , dictName (dict, seg->exon)
		      , dictName (dict, seg->chrom)
		      , seg->isDown ? seg->a1 : seg->a2
		      , seg->isDown ? seg->a2 : seg->a1
		      , seg->a2 - seg->a1 + 1 	
		      , sponge	      
		      ) ;
	    freeOutf ("\t%ld\t%ld\t%.1lf\n"
		      , seg->ns, seg->nbp
		      , (double)seg->nbp/(seg->ns ? seg->ns : 1)
		      /*  , (float)seg->nbp/(ln ? ln : 1) */
		    ) ;
	}
    }
  ac_free (h) ;

  return levelMax ;
} /* sxSponge */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/* position a set of SDF relative to a hierarchical geometry
 * the geometry is specified by a hierarchical set of wiggleremap 6 column files
 * each SDF is attributed to the first successful interval
 */

static int sxMapSdf (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (sx->outFileName, ".sdfPosition.txt", FALSE, h) ;
  int ok = 0, levelMax = 0, strand ;
  int a1, ii, jj, jjMax, chrom, level ;
  Array segs = arrayHandleCreate (100000, SPONGE, h) ;
  SPONGE *seg, *sdf ;
  DICT *spongeDict = dictHandleCreate (10000, h) ;
  DICT *levelNames = dictHandleCreate (10000, h) ;
  Array sdfs = arrayHandleCreate (10000, SPONGE, h) ;
  KEYSET sdfKs = keySetHandleCreate (h) ;
  BOOL stranded = sx->stranded ;
  BOOL isDown ;

  levelMax = sxSpongeParseFile (sx, spongeDict, levelNames, segs, 0, h) ;

  /* parse the objects to be positioned
   * we expect 5 columns: gene_nam exon_nam chrom a1 a2
   *           or may be: intron donor      chrom a1 a2
   * the mapping can be stranded or not stranded 
   */

  keySet (sdfKs, 0) = 0x7fffffff ;
  sxSpongeParseOneFile(sx, sx->sdfFileName, spongeDict, sdfs, sdfKs, 1, h) ;
  arraySort (segs, sdfOrder) ;
  arraySort (sdfs, sdfOrder) ;
  arrayCompress (sdfs) ;
  keySetSort (sdfKs) ; 
  keySetCompress (sdfKs) ;

  jjMax = arrayMax (segs) ;
  for (strand = 0 ; strand < (stranded ? 2 : 1) ; strand++)
    for (level = 1, jj = 0 ; level <= levelMax ; level++)
      for (ii = 0, sdf = arrp (sdfs, 0, SPONGE), seg = arrp (segs, jj, SPONGE) ; ii < arrayMax(sdfs) ; ii++, sdf++)
	{
	  if (sdf->level == 99)
	    continue ;
	  if (seg->level > level)
	    break ;
	  isDown = sdf->isDown ;
	  if (stranded && (isDown != strand ? TRUE : FALSE))
	    continue ;
	  a1 = sdf->a1 ;
	  chrom = sdf->chrom ; 
	  while (jj < jjMax && seg->level < level)
	    { seg++ ; jj++ ; }
	  while (jj < jjMax && seg->chrom < chrom && seg->level == level) 
	    { seg++ ; jj++ ; }
	  while (jj < jjMax && ((stranded && seg->isDown != isDown) || seg->a2 < a1)  && seg->chrom == chrom && seg->level == level) 
	    { seg++ ; jj++ ; }
	  if ( seg->chrom == chrom && seg->a1 <= a1 && seg->a2 >= a1 && seg->level == level && (!stranded || seg->isDown == isDown)) 
	    {
	      ok++ ;
	      aceOutf (ao, "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\n"
		       , dictName (spongeDict, sdf->gene)
		       , dictName (spongeDict, sdf->exon)
		       , dictName (spongeDict, sdf->chrom)
		       , isDown ? sdf->a1 : sdf->a2
		       , isDown ? sdf->a2 : sdf->a1
		       , dictName (levelNames, seg->level)
		       , dictName (spongeDict, seg->gene)
		       , dictName (spongeDict, seg->exon)   
		       , dictName (spongeDict, sdf->chrom)
		       , seg->isDown ? seg->a1 : seg->a2
		       , seg->isDown ? seg->a2 : seg->a1
		       ) ;
	      sdf->level = 99 ;
	    }	    
      }
  fprintf (stderr, "// sxMapSDF located %d SDF\n", ok) ;
  return ok ;
} /* sxMapSdf */

/*************************************************************************************/
/******************************** intron support *************************************/
/*************************************************************************************/

typedef struct sxIntronStruct { int mrna, chrom, target, x1, x2, a1, a2 ; double nd1, nd5, na1, na5, n1, n5, n8 ; } SXINTRON ;

/*************************************************************************************/
/* expect 6 columns, 
 * target x1 x2 chrom a1 a2
 * target could be a cosmid or an NM
 *   target a1 a2 is as found in the hits file
 *   is is remapped to chrom a1 a2
 * we extract the intron positions
 */
static int sxIntronSpongeMrnaRemap (SX *sx, AC_HANDLE h)
{
  int nn = 0 ;
  int oldMrna = 0,  mrna = 0, oldChrom = 0, chrom = 0, oldx2 = 0, target, x1, x2 = 0, olda2 = 0, a1, a2 = 0 ;
  char *cp ;
  ACEIN ai = 0 ;
  SXINTRON *iip ;
  AC_HANDLE h1 = ac_new_handle () ;

  if (! sx->mrnaRemapFileName)
    return 0 ;
  sx->mapDict = dictHandleCreate (1024, h) ;
  sx->remapDict = dictHandleCreate (1024, h) ;
  dictAdd (sx->mapDict,   "______toto", 0) ;
  dictAdd (sx->remapDict, "______toto", 0) ;

  ai = aceInCreate (sx->mrnaRemapFileName, 0,  h1) ;
  if (! ai)
     messcrash ("Cannot open the intron map file") ;

  sx->introns = arrayHandleCreate (10000, SXINTRON, h) ;
  sx->mrna2intron = keySetHandleCreate (h) ;

  while (ai && aceInCard (ai))
    {
      mrna = chrom = x1 = x2 = a1 = a2 = 0 ;
      /* target av HINV ... */
      cp = aceInWord (ai) ;
      if (!cp) 
	continue ;
      dictAdd (sx->remapDict, cp, &target) ;

      /* mrna name */
      aceInStep(ai,'\t') ; cp = aceInWord (ai) ;
      if (!cp) 
	continue ;
      dictAdd (sx->mapDict, messprintf ("%s:%s",dictName (sx->remapDict,target),cp), &mrna) ;

      aceInStep(ai,'\t') ; if (aceInInt (ai, &x1))
	{
	  aceInStep(ai,'\t') ; 
	  if (! aceInInt (ai, &x2))
	    messcrash ("missing coordinates in column 3 in %s line %d", 
		       sx->mrnaRemapFileName, aceInStreamLine (ai)) ;

	  aceInStep(ai,'\t') ; cp = aceInWord (ai) ;
	  if (cp) 
	    {
	      dictAdd (sx->remapDict, cp, &chrom) ;
	      aceInStep(ai,'\t') ; 
	      if (! aceInInt (ai, &a1) || ! aceInInt (ai, &a2))
		messcrash ("missing coordinates in column 5 or 6 in %s line %d", 
			   sx->mrnaRemapFileName, aceInStreamLine (ai)) ;
	    }
	  else
	    { /* duplicate column 1,2,3 */
	      messcrash ("missing chromo name in column 4 in %s line %d", 
			 sx->mrnaRemapFileName, aceInStreamLine (ai)) ;
	    }
	}
      else
	{
	   messcrash ("missing coordinates in column 2 in %s line %d", 
		       sx->mrnaRemapFileName, aceInStreamLine (ai)) ;
	}
	  
      if (x1 > x2)
	{ 
	  int x ;

	  x = x1 ; x1 = x2 ; x2 = x ; 
	  x = a1 ; a1 = a2 ; a2 = x ;
	}

      if (
	  (a1 < a2 && a2 - a1 != x2 - x1) || 
	  (a1 > a2 && a1 - a2 != x2 - x1)
	  )
	messcrash ("In %s x2-x1=(%d - %d) = %d  !=  a2-a1 = (%d-%d) = %d", sx->mrnaRemapFileName, x2, x1, x2 - x1, a2, a1, a2 - a1) ;

      if (mrna == oldMrna && chrom == oldChrom && x1 == oldx2 + 1)
	{
	  iip = arrayp (sx->introns, nn++, SXINTRON) ;
	  iip->mrna = mrna ;
	  iip->chrom = chrom ;
	  iip->target = target ;

	  if (keySet (sx->mrna2intron, mrna) == 0)
	    keySet (sx->mrna2intron, mrna) = nn ; /* ATTENTION, nn is 1 too big so never zero */
	  
	  iip->x1 = oldx2 ;
	  iip->x2 = x1 ;
	  /* the a1/a2 coordinates are shifted to give the standard name of the intron */
	  if (olda2 < a1)
	    {
	      iip->a1 = olda2 + 1 ;
	      iip->a2 = a1 - 1 ;
	    }
	  else
	    {
	      iip->a1 = olda2 - 1 ;
	      iip->a2 = a1 + 1;
	    }
	}
      oldMrna = mrna ; oldChrom = chrom ; oldx2 = x2 ; olda2 = a2 ;
    }

  fprintf (stderr, "Found %d intron in file %s\n", nn, sx->mrnaRemapFileName) ;
  ac_free (h1) ;

  return nn ;
} /* sxIntronSpongeMrnaRemap */

/*************************************************************************************/
/* scan the bestg .hits file
 * count the support of the last base of the donor exon,
 *                      the first base of the acceptor exon
 *                      the read aligned from -5 to 5
 *                      the read aligned from -8 to 8
 *                      the read aligned from -5 to 5 without error
 *                      the read aligned from -8 to 8 without error

 * 2013_09_02: Possibly this code is superseded by MAGIC phase ii2a,ii2b
 *  where we also count more directly the read supporting the introns
 *  in a recursive way
 */

static int sxIntronSpongeCountSupport (SX *sx, AC_HANDLE h)
{
  int nn = 0 ;
  int score, target, oldTarget = 0, nT, x1, x2, mrna, ii, ln, ali, mult ;
  double y ;
  const char *ccp ;
  char tagBuf[1024], tagBuf2[1024], geneBuf[1024] ;
  AC_HANDLE h1 = ac_new_handle () ;
  ACEIN ai = 0 ;
  SXINTRON *iip ;

  ai = aceInCreate (sx->bestHitFileName, 0,  h1) ;
  if (! ai)
    messcrash ("Cannot open the intron map file") ;

  memset (tagBuf, 0, sizeof(tagBuf)) ;


  if (! ai) 
    goto done ;
  while (aceInCard (ai))
    {
      /* tag name */
      ccp = aceInWord (ai) ;
      if (!ccp || *ccp == '#')
	continue ;
      if (strcmp (tagBuf, ccp)) 
	oldTarget = 0 ;
      strncpy (tagBuf, ccp, 1023) ;
      
	  /* score */
      aceInStep(ai,'\t') ; aceInNext (ai) ;
      if (! aceInInt (ai, &score))
	{
	  messcrash ("Missing score column 2 in AHIT input file %s, line %d:\n"
		     , sx->bestHitFileName, aceInStreamLine (ai)) ;
	}

      /* target class */
      aceInStep(ai,'\t') ; aceInNext (ai) ;
      aceInStep(ai,'\t') ; if (!(ccp = aceInWord (ai)))
	{
	  messcrash ("Missing target class column 3 in AHIT input file %s, line %d:\n"
		     , sx->bestHitFileName, aceInStreamLine (ai)) ;
	}
      if (! dictFind (sx->remapDict, ccp+2, &target))
	continue ; /* we cannot remap this target class */

      /* gene */
      aceInStep(ai,'\t') ; aceInNext (ai) ;
      aceInStep(ai,'\t') ; if (!(ccp = aceInWord (ai)))
	{
	  messcrash ("Missing gene column 4 in AHIT input file %s, line %d:\n"
		     , sx->bestHitFileName, aceInStreamLine (ai)) ;
	}
      if (target == oldTarget && ! strcmp (geneBuf, ccp)) continue ;
      strncpy (geneBuf, ccp, 1023) ;

      if (oldTarget && target != oldTarget && ! strcmp (tagBuf, tagBuf2)) continue ;
      strcpy (tagBuf2, tagBuf) ;
      oldTarget = target ;

      /* read the mrna so it is used as target */
      aceInStep(ai,'\t') ; if (!(ccp = aceInWord (ai)))
	{
	  messcrash ("Missing mRNA column 5 in AHIT input file %s, line %d:\n"
		     , sx->bestHitFileName, aceInStreamLine (ai)) ;
	}
      if (! dictFind (sx->mapDict, messprintf ("%s:%s",dictName (sx->remapDict,target),ccp), &mrna))
	continue ;
      ii = keySet (sx->mrna2intron, mrna) ;
      if (! ii)
	continue ;
      aceInStep(ai,'\t') ; aceInNext (ai) ;
      if (! aceInInt (ai, &ln))
	{
	  messcrash ("Missing length column 6 in AHIT input file %s, line %d:\n"
		     , sx->bestHitFileName, aceInStreamLine (ai)) ;
	}
      
      /* ali */
      aceInStep(ai,'\t') ; aceInNext (ai) ;
      if (! aceInInt (ai, &ali))
	{
	  messcrash ("Missing alignment-length 7 in AHIT input file %s, line %d:\n"
		     , sx->bestHitFileName, aceInStreamLine (ai)) ;
	}
      
      /* multiplicity */
      aceInStep(ai,'\t') ; aceInNext (ai) ;
      if (! aceInInt (ai, &mult))
	{
	  messcrash ("Missing multiplicity column 8 in AHIT input file %s, line %d:\n"
		     , sx->bestHitFileName, aceInStreamLine (ai)) ;
	}
      if (mult == 0) mult = 1;
      /* x1 */
      aceInStep(ai,'\t') ; aceInNext (ai) ;
      if (! aceInInt (ai, &x1))
	{
	  messcrash ("Missing x1 position column 9 in AHIT input file %s, line %d:\n"
		     , sx->bestHitFileName, aceInStreamLine (ai)) ;
	}
      
      /* x2 */
      aceInStep(ai,'\t') ; aceInNext (ai) ;
      if (! aceInInt (ai, &x2))
	{
	  messcrash ("Missing x2 position column 9 in AHIT input file %s, line %d:\n"
		     , sx->bestHitFileName, aceInStreamLine (ai)) ;
	}
      if (x1 > x2) { int dummy = x1 ; x1 = x2 ; x2 = dummy ; }
      /* snip */
      aceInStep(ai,'\t') ; if (!(ccp = aceInWord (ai)))
	{
	  messcrash ("Missing missmatch column 11 in AHIT input file %s, line %d:\n"
		     , sx->bestHitFileName, aceInStreamLine (ai)) ;
	}
      /* strncpy (snpBuf, ccp, 999) ; drop the snip */
      
      /* number of targets in this class */
      aceInStep(ai,'\t') ; aceInNext (ai) ;
      if (! aceInInt (ai, &nT))
	nT = 1  ;
      if (nT < 0) nT = -nT ;

      y = ((double)mult) / nT ;

      /* if we are not covering an intron, continue */

      for (ii--, iip = arrp (sx->introns, ii, SXINTRON) ; ii < arrayMax (sx->introns) ; iip++, ii++)
	{
	  if (iip->mrna != mrna || iip->x1 > x2)
	    break ;

	  if (x1 <= iip->x1 && x2 >= iip->x1)
	    iip->nd1 += y ;
	  if (x1 <= iip->x2 && x2 >= iip->x2)
	    iip->na1 += y ; 
	  if (x1 <= iip->x1 - 4 && x2 >= iip->x1)
	    iip->nd5 += y ;
	  if (x1 <= iip->x2 && x2 >= iip->x2 + 4)
	    iip->na5 += y ; 
	  if (x1 <= iip->x1 - 0 && x2 >= iip->x2 + 0)
	    iip->n1 += y ; 
	  if (x1 <= iip->x1 - 4 && x2 >= iip->x2 + 4)
	    iip->n5 += y ; 
	  if (x1 <= iip->x1 - 7 && x2 >= iip->x2 + 7)
	    iip->n8 += y ; 
	}
    }

 done:
  ac_free (h1) ;
  return nn ;
} /* sxIntronSpongeCountSupport */

/*************************************************************************************/

static int sxIntronSpongeReportOrder (const void *va, const void *vb)
{
  const SXINTRON *a = (const SXINTRON *) va, *b = (const SXINTRON *) vb ;
  int x ;
  x = a->chrom - b->chrom ; if (x) return x ;
  x = a->a1 - b->a1 ; if (x) return x ;
  x = a->a2 - b->a2 ; if (x) return x ;
  x = a->mrna - b->mrna ; if (x) return x ;
  return 0 ;
} /* sxIntronSpongeReportOrder */

/*************************************************************************************/

static void sxIntronSpongeDoMerge (SX *sx, AC_HANDLE h)
{
  AC_HANDLE h1 = ac_new_handle () ;
  ACEIN ai = aceInCreate (0, 0, h1) ;
  int ii = 0, chrom, a1, a2 ;
  float  na1, nd1, na5, nd5, n1, n5, n8 ;
  char *cp ;
  SXINTRON *up ;

  sx->introns = arrayHandleCreate (10000, SXINTRON, h) ;
  sx->remapDict = dictHandleCreate (1024, h) ;
  dictAdd (sx->remapDict, "______toto", 0) ;

  aceInSpecial (ai, "\t\n\\\"") ; /* so we do nothave to  aceInStep(ai, '\t') each time */

  while (aceInCard (ai))
    {
      cp = aceInWord (ai) ;
      if (!cp || *cp == '#' ) 
	continue ;
      dictAdd (sx->remapDict, cp, &chrom) ;    
      if (aceInInt (ai, &a1) && aceInInt (ai, &a2) &&
	  aceInFloat (ai, &nd1) && aceInFloat (ai, &na1) && 
	  aceInFloat (ai, &nd5) && aceInFloat (ai, &na5) && 
	  aceInFloat (ai, &n1) && aceInFloat (ai, &n5) && aceInFloat (ai, &n8)
	  )
	{
	  up = arrayp (sx->introns, ii++, SXINTRON) ;
	  up->chrom = chrom ;
	  up->mrna = 1 ;
	  up->a1 = a1 ;
	  up->a2 = a2 ;
	  up->nd1 = nd1 ;
	  up->na1 = na1 ;
	  up->nd5 = nd5 ;
	  up->na5 = na5 ;
	  up->n1 = n1 ;
	  up->n5 = n5 ;
	  up->n8 = n8 ;
	}
    }

  fprintf (stderr, "sxIntronSpongeMerge parsed %d lines\n", ii) ;
  ac_free (h1) ;
} /*sxIntronSpongeDoMerge */

/*************************************************************************************/

static int sxIntronSpongeReport (SX *sx)
{
  int nn = 0, ii, jj, iiMax = arrayMax (sx->introns) ;
  SXINTRON *up, *vp ;

  /* keep happy few */
  if (0) 
    {
      for (ii = jj = 0, up = vp = arrp (sx->introns, 0, SXINTRON) ; ii < iiMax ; up++, ii++)
	if (up->nd1 + up->na1 > 0)
	  {
	    if (ii > jj) *vp = *up ;
	    vp++ ; jj++ ;
	  }
      arrayMax (sx->introns) = iiMax = jj ;
    }
  /* sort and fuse, this destroys the mrna2intron correspondance */
  arraySort (sx->introns, sxIntronSpongeReportOrder) ;
  /* fuse double counts */
  for (ii = 0, up = arrp (sx->introns, 0, SXINTRON) ; ii < iiMax ; up++, ii++)
    for (jj = ii + 1, vp = up + 1 ; 
	 up->mrna && vp->mrna && jj < iiMax && vp->chrom == up->chrom && vp->a1 == up->a1 && vp->a2 == up->a2 ; 
	 vp++, jj++)
      {
	up->nd1 += vp->nd1 ;
	up->na1 += vp->na1 ;
	up->nd5 += vp->nd5 ;
	up->na5 += vp->na5 ;
	up->n1 += vp->n1 ;
	up->n5 += vp->n5 ;
	up->n8 += vp->n8 ;
	vp->mrna = 0 ;
      }
  /* report */
  freeOutf ("#Chromosome\tIntron start\tend\tTags supporting the last base of the donor exon\tTags supporting the first base of the acceptor exon\tTags supporting 5bp on the donor side\tTags supporting 5bp on the acceptor side\tTags supporting 1bp on each side\tTags supporting 5bp on each side\t8bp on each side (not necessarily without error)\n") ;
  for (ii = nn = 0, up = arrp (sx->introns, 0, SXINTRON) ; ii < iiMax ; up++, ii++)
    if (up->mrna)
      {
	nn++ ;
	freeOutf ("%s\t%d\t%d\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n"
		  , dictName (sx->remapDict, up->chrom)
		  , up->a1, up->a2
		  , up->nd1, up->na1
		  , up->nd5, up->na5
		  , up->n1, up->n5, up->n8
		  ) ;
      }
  fprintf (stderr, "sxIntronSpongeReport exported %d lines\n", nn) ;
  return nn ;
} /* sxIntronSpongeReport */

/*************************************************************************************/

static int sxIntronSponge (SX *sx)
{
  int nn = 0 ;
  AC_HANDLE h = ac_new_handle () ;

  sxIntronSpongeMrnaRemap (sx, h) ;
  sxIntronSpongeCountSupport (sx, h) ;
  nn = sxIntronSpongeReport (sx) ;   /* this destroys the mrna2intron correspondance */

  ac_free (h) ;
  return nn ;
} /* sxIntronSponge */

/*************************************************************************************/

static int sxIntronSpongeMerge (SX *sx)
{
  int nn = 0 ;
  AC_HANDLE h = ac_new_handle () ;

  sxIntronSpongeDoMerge (sx, h) ;
  nn = sxIntronSpongeReport (sx) ;   /* this destroys the mrna2intron correspondance */

  ac_free (h) ;
  return nn ;
} /* sxIntronSpongeMerge */

/*************************************************************************************/
/*************************************************************************************/
/* export useful statistics about power laws */
static int sxPowerTable (SX *sx)
{
  double a, n0, n[101], z, z0 ;
  int i, t, nt ;
  printf ("Power law parameters\nalpha\texp(alpha)\tT\tN") ;
  for (i = 0; i <= 14 ; i++)
    printf ("\tN(%d)",i) ;
  for (z0 = 8 ; z0 <= 12 ; z0 += 2)
    for (n0=10000000;n0<20000000;n0+=2000000)
      {
	a = log(z0) ;
	z = 1/z0 ; 
	n[0] = n0 ;
        for (nt = n0, t = 0, i = 1 ; i <= 100 ; i++)
	  {
	    n[i] = z * n[i-1] ;
	    t += i * n[i] ;
	    nt += n[i] ;
	  }
	printf ("\n%g\t%g\t%d\t%d", a, 1/z,(int)t, (int)nt) ;
	for (i = 0; i <= 14 ; i++)
	  freeOutf ("\t%d",(int)n[i]) ;
      }
  return 0 ;
} /* sxPowerTable */

/*************************************************************************************/

/* Histogramme des nombres de tags dans des fenetres de sx->width bp */
typedef struct siteStruct { int site, map, a1, a2 ; } SITE ;
static int sxSpongeWindowHisto (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  char *cp ;
  int NN = sx->width > 0 ? sx->width : 200, old1 = -1, old2 = -1, nn1 = 0, nn2 = 0, x,x1,x2, n1 = 0 ;
  Array excludedMaps = 0 ;
  Array *excludep ;
  int ix, ns = 0, nx = 0, a1, a2, map, oldmap = -1 ;
  Array hh = arrayHandleCreate (1000, int, h) ;
  ACEIN ai ; 
  SITE *s ;
  DICT *dict = dictHandleCreate (64,h) ;
  ACEOUT ao = aceOutCreate (sx->outFileName, "spongeHisto.txt", FALSE, h) ;

  if (sx->spongeWindowExcludeFileName)
    {
      if(0) printf ("sx->spongeWindowExclude=%s\n", sx->spongeWindowExcludeFileName);
      ai = aceInCreate (sx->spongeWindowExcludeFileName, 0, h) ;
      if (ai)
	{
	  excludedMaps = arrayHandleCreate (200, Array, h) ;
	  while (aceInCard (ai))
	    {
	      cp = aceInWord (ai) ;
	      if (!cp) continue ;
	      dictAdd (dict, cp, &map) ;
	      if (!aceInInt (ai, &a1) || !aceInInt (ai, &a2))
		continue ;
	      excludep = arrayp (excludedMaps, map, Array) ;
	      if (! *excludep)
		*excludep = arrayHandleCreate (2000, SITE, h) ;
	      nx = arrayMax (*excludep) ;
	      s = arrayp (*excludep, nx, SITE) ;
	      s->map = map ;
	      s->a1 = a1 ; s->a2 = a2 ;
	    }
	  ac_free (ai) ;
	}
    }

  sxNewExonsGetWiggle (sx, &SOLEXA_STEP, sx->wiggleFileName_f, h) ;
  
  messcrash ("2013_01_29: This code is incomplete, the direct parsing of stdin must be replaced by analysing the w_f array") ;
  ai = aceInCreate (0, 0, h) ;
  while (aceInCard (ai))
    {
      if (!aceInInt (ai, &x)   ||  /* position in cosmid */
	  ! aceInInt (ai, &n1) ||  /* number of tags */
	  !(cp = aceInWord (ai))   /* chrom name */
	  )
	continue ;
      dictAdd (dict, cp, &map) ;
      if (!aceInInt (ai, &a1))    /* position in chrom */
	continue ;
      if (excludedMaps &&
	  map < arrayMax(excludedMaps) &&
	  (excludep = arrp (excludedMaps, map, Array))
	  )
	{
	  BOOL found = FALSE ;
	 
	  nx = arrayMax (*excludep) ;
	  for (ix = 0 ; !found && ix < nx ; ix++)
	    {
	      s = arrayp (*excludep, ix, SITE) ;
	      if (s->a1 <= a1 && s->a2 >= a1)
		found = TRUE ; 
	      if (s->a1 > a1)
		break ;
	    }
	  if (found) /* excude the tags belonging to the known peaks */
	    { ns += n1 ; continue ;}
	}
      if (1)	    
	{
	  /* plain tiles */
	  x1 = (a1 + randfloat() - .5)/NN ; /* avoid sharp boundaries */
	  if (x1 != old1)
	    {
	      if (old1 >=0 )
		array (hh, nn1, int)++ ;
	      if(nn1 > 20)
		aceOutf (ao, "#\t%s\t%d\t%d\n", dictName (dict, map), NN*old1, nn1); 
	      nn1 = 0 ;
	      if (map == oldmap && x1 > old1+1) /* all the sections in between are empty */
		array (hh, 0, int) += x1 - old1 - 1 ;
	    }
	  old1 = x1 ;
	  oldmap = map ;
	  nn1 += n1 ;
	  /* shifted tiles */
	  if (0)
	    {
	      x2 = (a1+NN/2)/NN ;
	      if (x2 != old2)
		{
		  if (old2 >=0 )
		    array (hh, nn2, int)++ ;
		  nn2 = 0 ;
		  if (x2 > old2+1) /* all the sections in between are empty */
		    array (hh, 0, int) += x2 - old2 - 1 ;
		  nn2 = 0 ;
		}
	      old2 = x2 ;
	      nn2 += n1 ;
	    }
	}
    }

  aceOutf (ao, "#tags\t#windows\t Number of windows of size %dbp containing s given number of Solexa tags, %d tags excluded\n", NN, ns) ;
  for (n1 = 0 ; n1 < arrayMax (hh) ; n1++)
    aceOutf (ao, "%d\t%d\n", n1, arr (hh, n1, int)) ;
  ac_free (h) ;
  /* if (0)  gnuPlot (sx->spongeWindowFileName, sx->title) ; */
  return 0 ;
} /* sxSpongeWindowHisto */

/*************************************************************************************/
/********************************** Wiggles ******************************************/
/*************************************************************************************/

/* Accumulate result from 1 type
 * Hard coded 
 * Dict indexs the names of the Nms
 * aa is an Array of Array, each beeing a wiggle, with a tag count for every bp along the NM
 */
static void sxNmOneWiggle (SX *sx, char **manip, char **tissue
			   , int iTissue, BOOL stranded, DICT *dict, Array aa
			   , AC_HANDLE h)
{
  ACEIN ai = 0 ;
  char *fNam = 0 ;
  Array a = 0 ;
  int a1, a2, *ip, i, jj, nm, lane, line ;
  const char *ccp, *ccq, *prefix = sx->nmWiggle ;
  int lenPrefix = strlen (prefix) ;
  AC_HANDLE h1 = ac_new_handle () ;

  lenPrefix = strlen (prefix) ;
  for (lane = 0 ; lane < 60 ; lane++)
    {
      fNam = hprintf(h1, "PHITS_RefSeq.u/%s.%s.%d.RefSeq.u.hits",*manip,*tissue,lane) ;
      fNam = hprintf(h1, "PHITS_Sasha.u/%s.%s.%d.Sasha.u.hits",*manip,*tissue,lane) ;
      if (filName (fNam, 0, "r"))
	{
	  ai = aceInCreate (fNam, 0, h1) ;
	  line = 0 ;
	  while (aceInCard (ai))
	    {
	      line++ ;
	      for (i=0;i<4;i++)
		{ aceInStep(ai,'\t') ; ccp = aceInWord (ai) ; } /* NM name */
	      if (!ccp) continue ;
	      ccq = ccp + strlen(ccp) - lenPrefix ;
	      if (strcmp (ccq, prefix))
		continue ;
	      dictAdd (dict, ccp, &nm) ;
	      jj = 6*nm + 3*iTissue + stranded ;
	      a = array (aa, jj, Array) ;
	      if (a==0)
		array (aa, jj, Array) = a = arrayHandleCreate (1024, int, h) ;
	      aceInStep(ai,'\t') ; 
	      if (! aceInInt (ai, &a1)) continue ;
	      aceInStep(ai,'\t') ; 
	      if (aceInInt (ai, &a2))
		{
		  if (a1 > a2)
		    {
		      int a0 = a1 ; a1 = a2 ; a2 = a0 ;
		      if (stranded) jj++ ;
		    }
		  for (i = a1 ; i <= a2 ; i++)
		    {
		      ip = arrayp (a, i, int) ;
		      (*ip)++ ;
		    }
		}
	    }
	  ac_free (ai) ;
	}
    }
  ac_free (h1) ;
} /* sxNmOneWiggle */

/*************************************************************************************/
/* Export an ace file for the 4 wiggle of every registed NM */
static void sxNmWiggleExport (SX *sx, DICT *dict, Array aa, char **tissues)
{
  const char *nm ;
  char **tissue ;
  Array a ;
  char tag[800], *sName[] = {"not_stranded", "stranded", "anti_stranded", 0} ;
  int i, iNm, iTissue, jj, stranded, v, oldv ;

  for (iNm = 1 ; iNm <= dictMax(dict) ; iNm++)
    {
      nm = dictName (dict, iNm) ;
      freeOutf ("Sequence \"%s\"\n", nm) ; 
      for (stranded = 0 ; stranded < 3 ; stranded++)
	for (iTissue = 0, tissue = tissues ; *tissue ; iTissue++, tissue++)
	  {
	    sprintf(tag, "%s_%s", *tissue, sName[stranded]) ;
	    jj = 6*iNm + 3*iTissue + stranded ;
	    a = array (aa, jj, Array) ;
	    if (a)
	      {
		oldv = -1 ;
		for (i = 1; i < arrayMax(a) ; i++)
		  {
		    v = arr (a, i, int) ;
		    if (v!=oldv)
		      printf ("Solexa %s %d %d\n", tag, i, v) ;
		    oldv = v ;
		  }
	      }
	  }
      freeOut ("\n") ; 
    }
} /* sxNmWiggleExport */

/*************************************************************************************/
/* create a standed + a non stranded wiggle for all the NMs with correct suffix
 * Hard coded 
 */
static void sxNmWiggle (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  int iTissue ;
  DICT *dict = dictHandleCreate (100000, h) ;
  Array aa = arrayHandleCreate (30000, Array, h) ;
  char **manip, **tissue, *tissues[] = {"Brain", "UHR", 0} ;
  char *stranded[] = {"ILM", "LIF", "HEL", 0} ;
  char *notStranded[] = {"ILMnoSt", "HH", 0} ;

  /* ventilate all data per NM */
  for (manip = stranded ; *manip ; manip++)
    for (iTissue = 0, tissue = tissues ; *tissue ; iTissue++, tissue++)
      sxNmOneWiggle (sx, manip, tissue, iTissue, TRUE, dict, aa, h) ;
  for (manip = notStranded ; *manip ; manip++)
    for (iTissue = 0, tissue = tissues ; *tissue ; iTissue++, tissue++)
      sxNmOneWiggle (sx, manip, tissue, iTissue, FALSE, dict, aa, h) ;

  /* export */
  sxNmWiggleExport (sx, dict, aa, tissues) ;

  ac_free (h) ;
  return ;
} /* seNmWiggle */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/* create a binary wiggle
  fprintf (stderr, 
	   "// -cosmidWiggle target : Read stdin or -hit, export a wiggle per target,\n"
	   "//       Target has 4 columns: target map a1 a2, target extends on map from  a1 to a2,\n"
	   "//       map a1 a2 is expected in col 1,2,3 of input, each tag is trimmed by 8bp at both end\n"
	   "//       the wiggle is exported a variable step with 5bp precision\n"
	   );
 * Hard coded 
 */

typedef struct sxwqStruct { int t, map, a1, a2 ; BOOL isDown ; KEYSET hits, upHits, downHits ;} SXWQ ;

/*************************************************************************************/

int sxwtOrder (const void *va, const void *vb)
{
  const SXWQ *a = (const SXWQ *) va, *b = (const SXWQ *) vb ;
  int x ;
  x = a->map - b->map ; if (x) return x ;
  x = a->a1 - b->a1 ; if (x) return x ;
  x = a->t - b->t ;
  return x ;
}


/*************************************************************************************/
/* expect 4 columns, 
 * target map a1 a2
 * target could be a cosmid or an NM, map a1 a2 is as found col 1/2/3 of the stdin hits file
 */
static int sxTargetWiggleGetTargets (SX *sx, Array targets, DICT *dict)
{
  int nn = 0, t, map, a1, a2 ;
  char *cp ;
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (sx->cosmidWiggle, 0,  h) ;
  SXWQ *sxwt, *sxwt2 ;

  while (ai && aceInCard (ai))
    {
      t = map = a1 = a2 = 0 ;
      cp = aceInWord (ai) ;
      if (cp) 
	dictAdd (dict, cp, &t) ;
      aceInStep(ai,'\t') ; cp = aceInWord (ai) ;
      if (cp) 
	dictAdd (dict, cp, &map) ;
      aceInStep(ai,'\t') ; 
      if (! aceInInt (ai, &a1))
	continue ;
      aceInStep(ai,'\t') ; 
      if (t && map && aceInInt (ai, &a2))
	{
	  sxwt = arrayp (targets, nn++, SXWQ) ;
	  sxwt->t = t ;
	  sxwt->map = map ;
	  if (a1 < a2)
	    {
	      sxwt->a1 = a1 ;
	      sxwt->a2 = a2 ;
	      sxwt->isDown = TRUE ;
	    }
	  else
	    {
	      sxwt->a1 = a2 ;
	      sxwt->a2 = a1 ;
	      sxwt->isDown = FALSE ;
	    }
	}
    }
  arraySort (targets, sxwtOrder) ;
  if (sx->wiggleFeature && arrayMax(targets) > 1)
    {
      int i ;
      for (i = 1, sxwt = arrp (targets, 0, SXWQ), sxwt2 = sxwt - 1 ; targets && i < nn ; sxwt++, sxwt2++, i++)
	{
	  if (sxwt2->map == sxwt->map && sxwt2->a2 >= sxwt->a1)
	    sxwt2->a2 = sxwt->a1 - 1 ;
	}
    }
  ac_free (h) ;
  return nn ;
} /* sxTargetWiggleGetTargets */

/*************************************************************************************/

static int sxTargetWiggleGetHits(SX *sx, Array targets, DICT *dict, BOOL stranded, AC_HANDLE h0)
{
  int nn = 0, map, t, a1, a2, ii, u1, u2, x, x1, x2 ;
  char *cp ;
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = 0 ;
  SXWQ *sxwt ;
  BOOL isDown ;
  Array hits = 0 ;
  int delta = sx->delta ;
  int mult ;

  ai = aceInCreate (sx->hitFile, 0, h) ;
 
  aceInSpecial (ai, "\t\n\\\"") ;

  while (ai && aceInCard (ai))
    {
      map = a1 = a2 = 0 ;
      cp = aceInWord (ai) ;
      if (cp)
	{
	  if (sx->mrnaWiggle)
	    dictAdd (dict, cp, &t) ;
	  if (dictFind (dict, cp, &map) &&
	      aceInInt (ai, &a1) && aceInInt (ai, &a2)
	      )
	    {
	      mult = 1 ;
	      if ((cp = aceInWord (ai)))
		mult = fastcMultiplicity (cp, 0, 0) ;
	      if (a1 <= a2)
		{ 
		  isDown = TRUE ;
		  if (a2 == a1 + 1) a2 = a1 ; /* assume a width 2 tag is just an oriented width 1 tag */
		}
	      else
		{
		  int a0 ;
		  if (a2 == a1 - 1) a2 = a1 ; /* assume a width 2 tag is just an oriented width 1 tag */
		  a0 = a1 ; a1 = a2 ; a2 = a0 ;
		  isDown = FALSE ;
		}
	      
	      if (a1 < a2 - 2*delta)
		{ a1 += delta ; a2 -= delta ; }
	      else
		a1 = a2 = (a1 + a2)/2 ;
	      if (sx->mrnaWiggle) /* export the continuous shape of the tag, cutting delta at each end */
		{
		  sxwt = arrayp (targets, t, SXWQ) ;
		  sxwt->t = t ;
		  sxwt->isDown = TRUE ;
		  if (!stranded || isDown)
		    {
		      x1 = a1 ;
		      x2 = a2 ;
		      if (x1 <= x2)
			{
			  hits = sxwt->hits ;
			  if (!hits) 
			    hits = sxwt->hits = keySetHandleCreate (h0) ;	
			  x = keySet (hits, x2+1) ;
			  if (x1>1)
			    x = keySet (hits, x1-1) ;
			  for (x = x2 ; x >= x1 ; x--)
			    keySet (hits, x)+=mult ;
			}
		    }
		}
	      else
		{
		  for (ii = 0, sxwt = arrayp (targets, 0, SXWQ) ; ii < arrayMax(targets) ; ii++, sxwt++)
		    {
		      if (sxwt->map == map)
			{
			  u1 = sxwt->a1 ; if (u1 < a1) u1 = a1 ;
			  u2 = sxwt->a2 ; if (u2 > a2) u2 = a2 ;
			  if (u1 <= u2)
			    {
			      if (! sx->stranded)
				{
				  hits = sxwt->hits ;
				  if (!hits) 
				    hits = sxwt->hits = keySetHandleCreate (h0) ;
				}
			      else if (isDown)
				{
				  hits = sxwt->downHits ;
				  if (!hits) 
				    hits = sxwt->downHits = keySetHandleCreate (h0) ;
				}
			      else
				{
				  hits = sxwt->upHits ;
				  if (!hits) 
				    hits = sxwt->upHits = keySetHandleCreate (h0) ;
				}
			      if (sxwt->isDown)
				{ x1 = u1 - sxwt->a1 + 1 ; x2 = u2 - sxwt->a1 + 1 ; }
			      else
				{ x1 = sxwt->a2 - u2 + 1 ; x2 = sxwt->a2 - u1 + 1 ; }
			      x = keySet (hits, x2+1) ;
			      if (x1>1)
				x = keySet (hits, x1-1) ;
			      for (x = x2 ; x >= x1 ; x--)
				keySet (hits, x) += mult ;
			    }
			}
		    }
		}
	    }
	}
    }
  
  ac_free (h) ;
  return nn ;
} /* sxTargetWiggleGetHits */

/*************************************************************************************/

static void sxTargetWiggleExport (SX *sx, Array targets, DICT *dict)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii, i, j, n, oldn ;
  SXWQ *sxwt ;
  BOOL isBin = FALSE ;
  int isDownTrack = 1 ;
  char *bName ;
  unsigned char cc ;
  Array hits ;

  if (isBin) 
    oneByteInitialize (FALSE) ;
  for (ii = 0 ; ii < arrayMax (targets) ; ii++)
    {
      sxwt = arrp (targets, ii, SXWQ) ;
      isDownTrack = sxwt->isDown ? 1 : -1 ;
      hits = sxwt->hits ;
      if (! sxwt->hits && !sxwt->upHits && !sxwt->downHits)
	continue ;

      if (sx->wiggleDir)
	{
	  int iHH ;
	  char *wNam = 0, *wNam2 = 0, *sNam = 0 ;
	  ACEOUT ao = 0 ;

	  for (iHH = 0 ; iHH < 3 ;iHH++)
	    {
	      switch (iHH)
		{
		case 0: 
		  hits = sxwt->hits ;
		  wNam = hprintf (h, "%s/%s/%s.slx"
				  , sx->wiggleDir
				  , sx->trackName
				  , dictName (dict, sxwt->t)
				  ) ;
		  wNam2 = hprintf (h, "%s/%s.slx"
				   , sx->trackName
				   , dictName (dict, sxwt->t)
				   ) ;		  
		  sNam = "" ;
		  break ;
		case 1: 
		  hits = sxwt->downHits ; 
		  wNam = hprintf (h, "%s/%s/%s.f.slx"
				  , sx->wiggleDir
				  , sx->trackName
				  , dictName (dict, sxwt->t)
				  ) ;
		  wNam2 = hprintf (h, "%s/%s.f.slx"
				   , sx->trackName
				   , dictName (dict, sxwt->t)
				   ) ;
		  sNam = ".f" ;
		  break ;
		case 2: 
		  hits = sxwt->upHits ; 
		  wNam = hprintf (h, "%s/%s/%s.r.slx"
				  , sx->wiggleDir
				  , sx->trackName
				  , dictName (dict, sxwt->t)
				  ) ;
		  wNam2 = hprintf (h, "%s/%s.r.slx"
				   , sx->trackName
				   , dictName (dict, sxwt->t)
				   ) ;
		  sNam = ".r" ;
		  break ;
		}
	      if (!hits || !keySetMax (hits))
		continue ;
	      ao = aceOutCreateToFile (wNam, "w", h) ;
	      if (!ao)
		continue ;

	      if (sx->period == 0) /* export a value each time the coverage changes */
		{
		  for (i = oldn = 0 ; hits && i < keySetMax (hits) ; i++)
		    {
		      n = keySet (hits, i) ;
		      if (1 || n != oldn)
			{
			  if (1 || n || sx->mrnaWiggle)
			    freeOutf ("%d\t%d\n"
				     , i, n /* , dictName (dict, sxwt->map), sxwt->a1 + i - 1 */
				     ) ;
			}
		      oldn = n ;
		    }
		  if (oldn > 0 && sx->mrnaWiggle)
		    freeOutf ("%d\t%d\n"
			     , i, 0 /* , dictName (dict, sxwt->map), sxwt->a1 + i - 1 */
			     ) ;
		}
	      else /* export a value at regular intervals */
		{
		  for (oldn = 0, i = 1 - (sxwt->a1 % sx->period) ; hits && i < (int)keySetMax (hits) ; i += sx->period)
		    {
		      if(i<0) continue ;
		      n = keySet (hits, i) ;
		      if (!n && !oldn) continue ;
		      if (!oldn && i > 1)
			freeOutf ("%d\t%d\n"
				 , i - sx->period, oldn /* , dictName (dict, sxwt->map), sxwt->a1 + i - sx->period */
				 ) ;
		      freeOutf ("%d\t%d\n"
			       , i, n /* , dictName (dict, sxwt->map), sxwt->a1 + i  */
			       ) ;
		      oldn = n ;
		    }
		}
	  
	      ac_free (ao) ;
	      freeOutf ("Sequence %s\nSolexa %s%s %s\n\n"
			, ac_protect (dictName (dict, sxwt->t), h)
			, sx->trackName
			, sNam
			, wNam2
			) ;
	    }
	}
      else
	{
	  freeOutf ("%s %s\n"
		    , sx->mrnaWiggle ? "mRNA" : "Sequence"
		    , ac_protect (dictName (dict, sxwt->t), h)
		    ) ;
	  if (isBin) 
	    {
	      bName = ac_protect (hprintf(h,  "%s/%s.sxb", dictName (dict, sxwt->t), sx->trackName), h) ;
	      freeOutf ("Solexa %s %s\n\n", sx->trackName, bName);
	      freeOutf ("Binary %s\n", bName) ;
	      for (i = 0 ; i < keySetMax (hits) ; i++)
		{
		  cc = oneByteEncode (keySet(hits, i)) ;
		  j = 40 + cc ; if (j < 126) j = 126 ;
		}
	    }
	  else
	    {	 
	      int iHH, j ;

	      for (iHH = 0 ; iHH < 3 ;iHH++)
		{
		  switch (iHH)
		    {
		    case 0: hits = sxwt->hits ; isDownTrack = 0 ; break ;
		    case 1: hits = sxwt->downHits ; isDownTrack = 1 ; break ;
		    case 2: hits = sxwt->upHits ; isDownTrack = -1 ; break ;
		    }
		  if (!hits || !keySetMax (hits))
		    continue ;
		  for (i = oldn = 0 ; i < keySetMax (hits) ; i++)
		    {
		      n = keySet (hits, i) ;
		      if (n != oldn)
			{
			  if (!sx->wiggleFeature)
			    freeOutf ("Solexa  %s %d %d\n"
				      , sx->trackName, i, n
				      ) ;
			  else if (n)
			    { 
			      for (j = i + 2 ; j < keySetMax (hits) && keySet (hits, j) == n ; j++) {} ; j-- ;
			      freeOutf ("Feature %s %d %d %d %s\n"
					, sx->wiggleFeature
					, isDownTrack == -1 ? j - 1 : i
					, isDownTrack == -1 ? i - 1 : j 
					, n, sx->trackName
					) ;
			    }
			}
		      oldn = n ;
		    }
		  if (oldn > 0 && sx->mrnaWiggle)
		    freeOutf ("Solexa %s %d %d\n", sx->trackName, i, 0) ;
		}
	      freeOutf ("\n") ;
	    }
	}
    }
  ac_free (h) ;
  return ;
} /*sxTargetWiggleExport */

/*************************************************************************************/

static void sxTargetWiggle (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  DICT *dict = dictHandleCreate (100000, h) ;
  Array targets =  arrayHandleCreate (1000, SXWQ, h) ;

  if (!sx->trackName) 
    messcrash ("missign argument -trackName") ;
  if (! sx->mrnaWiggle)
    sxTargetWiggleGetTargets(sx, targets, dict) ;
  sxTargetWiggleGetHits (sx, targets, dict, sx->stranded, h) ;
  sxTargetWiggleExport (sx, targets, dict) ;

  ac_free (h) ;
  return ;
} /* sxTargetWiggle */


/*************************************************************************************/
/* In exome study, compute the relative dosage of each gene */
static void sxGeneDosage (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  int MXRUN = 130 ;
  Array aa = arrayHandleCreate (20000, float, h) ; 
  Array kbs = arrayHandleCreate (1000000, float, h) ;
  Array lns = arrayHandleCreate (20000, float, h) ;
  Array covers = arrayHandleCreate (20000, float, h) ;
  Array gMedian = arrayHandleCreate (20000, float, h) ;
  Array rMedian = arrayHandleCreate (200, float, h) ;
  KEYSET badGenes =  keySetHandleCreate (h) ; 
  KEYSET g2c = keySetHandleCreate (h) ;  /* gene chromosome */
  KEYSET g2x = keySetHandleCreate (h) ;  /* gene coordinate on the chromosome */
  DICT *geneDict = dictHandleCreate (20000, h) ;
  DICT *runDict = dictHandleCreate (200, h) ;
  DICT *chromDict = dictHandleCreate (30, h) ;
  ACEIN ai = aceInCreate (0, 0, h) ;
  ACEOUT ao = aceOutCreate (sx->outFileName, ".txt", 0, h) ;
  const char *ccp ;
  int i, n, x, pass ;
  int maxRunDebug = 100000 ;  /* may be used to limit the size of the exported table */
  int maxGeneDebug = 10000000 ; 
  int chrom, gene, run, nRuns = 0, nGoodGenes = 0 ;
  BOOL inside = FALSE ;
  float kb, globalScale = 1 ;
  double T ;

  if (sx->exomeDosage || sx->variantDosage)
    aceInSpecial (ai, "\n") ;

  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#') continue ;

      if (sx->variantDosage)
	{
	  char buf[1000], cc, *cp, *cq ;
	  float snp, cover ;

	  if (! strcmp(ccp, "Run"))
	    {
	      aceInStep (ai, '\t') ;
	      ccp = aceInWord (ai) ;	
	      dictAdd (runDict, ccp, &run) ;
	      continue ;
	    }
	  if (! strncmp(ccp, "Ghs", 3))
	    {
	      dictAdd (runDict, ccp, &run) ;
	      continue ;
	    }
	  
	  dictAdd (geneDict, ccp, &gene) ;
	  strncpy (buf, ccp, 999) ;
	  for (cq = cp = buf ; *cq && *cq != ':' ; cq++) ;
	  cc = *cq ; *cq = 0 ;
	  dictAdd (chromDict, cp, &chrom) ; 
	  keySet (g2c, gene) = chrom ;
	  *cq++ = cc ;
	  if (cc)
	    {
	      for (cp = cq ; *cq >= '0' && *cq <= '9' ; cq++) ;
	      cc = *cq ; *cq = 0 ;
	      x = atoi(cp) ;  
	      keySet (g2x, gene) = x ;
	    } 
	  array (lns, gene, float) = 1 ;

	  aceInStep (ai, '\t') ;
	  ccp = aceInWord (ai) ;	
	  dictAdd (runDict, ccp, &run) ;
	  
	  if (dictMax(runDict) > MXRUN - 5)
	    messcrash ("Too many runs, please recompile") ;

	  aceInStep (ai, '\t') ;
	  aceInFloat (ai, &snp) ;
	  aceInStep (ai, '\t') ;
	  aceInFloat (ai, &cover) ;
	  kb =  cover ? 100 * snp/cover : 0 ;
	  if (kb < 10) kb = 0 ;
	  else if (kb > 90) kb = 100 ;
	  else kb = 50 ;
	  if (cover < 10) kb = -1 ;
	  array (kbs, MXRUN * gene + run, float) = kb ;

          continue ;
	}
      if (sx->exomeDosage)
	{
	  if (! strcmp(ccp, "Run"))
	    {
	      aceInStep (ai, '\t') ;
	      ccp = aceInWord (ai) ;	
	      dictAdd (runDict, ccp, &run) ;
	      continue ;
	    }
	  if (! strncmp(ccp, "Ghs", 3))
	    {
	      dictAdd (runDict, ccp, &run) ;
	      continue ;
	    }

	  if (strncmp(ccp, "level_",6)) continue ;

	  aceInStep (ai, '\t') ;
	  ccp = aceInWord (ai) ;	
	  dictAdd (runDict, ccp, &run) ;	  
	  if (dictMax(runDict) > MXRUN - 5)
	    messcrash ("Too many runs, please recompile") ;

	  aceInStep (ai, '\t') ;	  
	  ccp = aceInWord (ai) ;
	  gene = 0 ;
	  
	  if (1) /* per gene */
	    {
	      dictAdd (geneDict, ccp, &gene) ;
	      aceInStep (ai, '\t') ;	  
	      ccp = aceInWord (ai) ;
	    }
	  else  /* per pseudo exon */
	    {
	      aceInStep (ai, '\t') ;	  
	      ccp = aceInWord (ai) ;
	      dictAdd (geneDict, ccp, &gene) ;
	    }
	  aceInStep (ai, '\t') ;	  
	  ccp = aceInWord (ai) ;
	  dictAdd (chromDict, ccp, &chrom) ;
	  keySet (g2c, gene) = chrom ;
	  aceInStep (ai, '\t') ;	  
	  aceInInt (ai, &i) ; /* start */
	  aceInStep (ai, '\t') ;	  
	  aceInInt (ai, &x) ; /* end  */

	  keySet (g2x, gene) = (i+x)/2 ;

	  aceInStep (ai, '\t') ;	  
	  aceInInt (ai, &x) ; /* length */
          array (lns, gene, float) += x ;
	  aceInStep (ai, '\t') ;	  
	  aceInInt (ai, &x) ;  /* threshold, drop it */
	  aceInStep (ai, '\t') ;	  
	  aceInInt (ai, &x) ; /* length covered, drop it */
	  aceInStep (ai, '\t') ;	  
	  aceInFloat (ai, &kb) ; /* coverage */

	  array (kbs, MXRUN * gene + run, float) += kb/1000 ;
	  continue ;
	}
      if (! strcmp (ccp, "Gene"))
	{
	  ccp = aceInWord (ai) ;
	  gene = 0 ;
	  dictFind (geneDict, ccp, &gene) ; /* recognize the genes */
	  inside = TRUE ;
	}
      else if (! inside)
	{
	  /* register the genes in the mapping order */
	  dictAdd (chromDict, ccp, &chrom) ;
	  aceInStep (ai, '\t') ;
	  aceInInt (ai, &x) ;
	  aceInStep (ai, '\t') ;
	  ccp = aceInWord (ai) ;
	  dictAdd (geneDict, ccp, &gene) ;
	  keySet (g2c, gene) = chrom ;
	  keySet (g2x, gene) = x ;
	}
      else if (! strcmp (ccp, "Run_U"))
	{       
	  /* register the kb per gene, per run */
	  ccp = aceInWord (ai) ;
	  dictAdd (runDict, ccp, &run) ;
	  if (dictMax(runDict) > MXRUN - 5)
	    messcrash ("Too many runs, please recompile") ;
	  for (i = 3 ; i <= 7 ;  i++) /* drop columns 3 to 7 */
	    ccp = aceInWord (ai) ; 
	  aceInFloat (ai, &kb) ;
 
	  if (gene > 0 && gene < maxGeneDebug && run < maxRunDebug)
	    array (kbs, MXRUN * gene + run, float) = kb ;
         }
    }

  /* we now have a matrix of coverage of each gene in each run
     we want to compute the theoretical value, a bit like in the chi-square test
     and measure the relative coverage as
         100 * observed * global-total / (sum of row * sum of column)
     but to renove wild points, we replace the sum of each row and each column
     by the median * the number of elements
  */

  /* compute the median for each gene */
  for (gene = 1 ; gene <= dictMax (geneDict) ; gene++)
    {
      aa = arrayReCreate (aa, 1000, float) ; 
      for (run = 1 ; run < maxRunDebug && run <= dictMax (runDict) ; run++)
	array (aa, run - 1, float) = array (kbs, MXRUN * gene + run, float) ;
      arraySort (aa, floatOrder) ;
      n = arrayMax (aa) ;
      kb = array (aa, n/2, float) ;
      array (gMedian, gene, float) = n * kb ;
      array (covers, gene, float) = n * kb/(1+array(lns,gene,float)) ;
      nRuns = n ;
    }
  
  /* compute the coverage of each gene and find the limits */
  {
    float c1, c2, cover ;
    
    arraySort (covers, floatOrder) ;
    n = arrayMax(covers) ;
    c1 = array (covers, 30 * n/100.0, float) ;
    c2 = array (covers, 98 * n/100.0, float) ;
    
    T = 0 ;
    for (gene = 1 ; gene <= dictMax (geneDict) ; gene++)
      {
	cover = array (gMedian, gene, float)/(1+array(lns,gene,float)) ;
	if ((sx->variantDosage && kb >= 0) || (cover >= c1 && cover <= c2))
	  T +=  nRuns * kb ;
	else
	  keySet (badGenes, gene) =  1 ;
      }
  }
  /* compute the median for each run excluding exome segments on X and Y */
  for (run = 1 ; run <= dictMax (runDict) ; run++)
    {
      aa = arrayReCreate (aa, 1000, float) ; 
      for (i = 0, gene = 1 ; gene <= dictMax (geneDict) ; gene++)
	{
	  ccp = dictName (chromDict, keySet (g2c, gene)) ;
	  if (! keySet (badGenes, gene) && ! strstr (ccp, "X") && ! strstr (ccp, "Y"))
	    array (aa, i++, float) = array (kbs, MXRUN * gene + run, float) ;
	}
      arraySort (aa, floatOrder) ;
      n = arrayMax (aa) ;
      kb = array (aa, n/2, float) ;
      array (rMedian, run, float) = n * kb ;
      T +=  n * kb ;
      nGoodGenes = n ;
    }
  T /= 2 ;  /* half of (sum lines + sum columns ) */

  for (pass = 0 ; pass < (sx->variantDosage ? 1 : 2) ; pass++)
    {
      globalScale = 1 ;
      /* regle de trois */  
      if (pass == 1)
	{
	  aa = arrayReCreate (aa, 1000, float) ; i = 0 ;
	  for (gene = 1 ; gene <= dictMax (geneDict) ; gene++)
	    if (! keySet (badGenes, gene))
	      for (run = 1 ; run <= dictMax (runDict) ; run++)
		{
		  kb = 100 * T * array (kbs, MXRUN * gene + run, float) / 
		    (1 + array (gMedian, gene, float) * array (rMedian, run, float) ) ;
		  array (kbs, MXRUN * gene + run, float) = kb ; 
		  array (aa, i++, float) = kb ;
		}
	  arraySort (aa, floatOrder) ;
          globalScale = array (aa, i/2, float) ;
	  if (globalScale < 1) globalScale = 1 ;
	  globalScale = 100/globalScale ;    /* this will recenter the global median to 100 */
	}

      if (pass)
	ao = aceOutCreate (sx->outFileName, ".index.txt", 0, h) ;
      /* caption line: the name of the runs */
      aceOutf (ao, "# %s\n", timeShowNow ()) ;
      aceOutf (ao, pass ? "# Number of kb aligned in each gene" : "# Relative number of kb aligned in each gene") ;
      aceOutf (ao, "\nChromosome\tPosition\tGene\tNumber of runs with index < 5\tNumber of runs with index < 40\tNumber of runs with index > 200\tNumber of runs  with index > 500\tMedian") ;
      for (run = 1 ; run <= dictMax (runDict) ; run++)
	aceOutf (ao, "\t%s", dictName (runDict, run)) ;
      
      if (1)
	{
	  /* global stat lines */

	  aceOutf (ao, "\n\t\t\t\tMedian kb aligned per genes\t\t\t", x) ;
	      for (run = 1 ; run <= dictMax (runDict) ; run++)
		aceOutf (ao, "\t%.1f", array (rMedian, run, float)/nGoodGenes) ;
	      
	  for (x = 50 ; x >= 0 ; x -= 10) 
	    {
	      if (x == 0) x = 5 ;
	      aceOutf (ao, "\n\t\t\t\tNumber of genes with index < %d\t\t\t", x) ;
	      for (run = 1 ; run <= dictMax (runDict) ; run++)
		{
		  for (n = 0, gene = 1 ; gene <= dictMax (geneDict) ; gene++)
		    if (! keySet (badGenes, gene))
		      {
			kb = globalScale * array (kbs, MXRUN * gene + run, float) ;
			if (kb > -1 && kb < x)
			  n++ ;
		      }
		  aceOutf (ao, "\t%d", n) ;
		}
	    }
	  for (x = 200 ; x<= 1000 ; x += 100)
	    {
	      aceOutf (ao, "\n\t\t\t\tNumber of genes with index > %d\t\t\t", x) ;
	      for (run = 1 ; run <= dictMax (runDict) ; run++)
		{
		  for (n = 0, gene = 1 ; gene <= dictMax (geneDict) ; gene++)
		    if (! keySet (badGenes, gene))
		      {
			kb = globalScale * array (kbs, MXRUN * gene + run, float) ;
			if (kb > -1 && kb > x)
			  n++ ;
		      }
		  aceOutf (ao, "\t%d", n) ;
		}
	    }
	}
      aceOutf (ao, "\n") ;
      /* exort one line per gene */
      for (gene = 1 ; gene <= dictMax (geneDict) ; gene++)
	{
	  if (keySet (badGenes, gene))
	    continue ;

	  aceOutf (ao, "%s\t%d\t%s"
		   , dictName (chromDict, keySet (g2c, gene))
		   , keySet (g2x, gene)
		   , dictName (geneDict, gene)
		   ) ;
	
	  for (x = 5 ; x <= 40 ; x += 35) 
	    {
	      for (n = 0, run = 1 ; run <= dictMax (runDict) ; run++)
		{
		  kb = globalScale * array (kbs, MXRUN * gene + run, float) ;
		  if (kb > -1 && kb < x)
		  n++ ;
		}
	      aceOutf (ao, "\t%d", n) ;
	    }
	  for (x = 200 ; x<= 500 ; x += 300)
	    {
	      for (n = 0, run = 1 ; run <= dictMax (runDict) ; run++)
		{
		  kb = globalScale * array (kbs, MXRUN * gene + run, float) ;
		  if (kb > -1 && kb > x)
		  n++ ;
		}
	      aceOutf (ao, "\t%d", n) ;
	    }

	  aceOutf (ao, "\t%.1f", array (gMedian, gene, float)/nRuns) ;
	  /* exporte the dosage = regle de trois */  
	  for (run = 1 ; run < maxRunDebug && run <= dictMax (runDict) ; run++)
	    {
	      kb = globalScale * array (kbs, MXRUN * gene + run, float) ;
	      aceOutf (ao, "\t%.1f", kb) ;
	    }
	  aceOut (ao, "\n") ;
	}
      aceOut (ao, "\n") ;	
    }
  ac_free (h) ;
}  /* sxGeneDosage */
 
/*************************************************************************************/
/* In stranded case, use the E LR FR wiggles to locate transcription start and end sites */

static void sxFlagTranscriptEnds (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  Array w_LF = 0, w_LR = 0, w_RF = 0, w_RR = 0 ;
  char *cp ;
  int nn = 0 ;
  int min = sx->minimalSupport ;
  int minRatio = 5 ;
  int isDown ;
  Array ww = arrayHandleCreate (100000, SXX, h) ;

  if (!  sx->wiggleFileName)
    messcrash ("Option -flagTranscriptEnds requires option -wiggle fileName") ;
  
  cp = hprintf (h, "%s.ELF.BF.gz", sx->wiggleFileName) ;
  if (! (w_LF = sxNewExonsGetWiggle (sx, &SOLEXA_STEP, cp, h)))
    messcrash ("-flagTranscriptEnds: missing ELF wiggle file %s\n", cp) ;
      
  cp = hprintf (h, "%s.ELR.BF.gz", sx->wiggleFileName) ;
  if (! (w_LR = sxNewExonsGetWiggle (sx, &SOLEXA_STEP, cp, h)))
    messcrash ("-flagTranscriptEnds: missing ELR wiggle file %s\n", cp) ;
      
  cp = hprintf (h, "%s.ERF.BF.gz", sx->wiggleFileName) ;
  if (! (w_RF = sxNewExonsGetWiggle (sx, &SOLEXA_STEP, cp, h)))
    messcrash ("-flagTranscriptEnds: missing ERF wiggle file %s\n", cp) ;
      
  cp = hprintf (h, "%s.ERR.BF.gz", sx->wiggleFileName) ;
  if (! (w_RR = sxNewExonsGetWiggle (sx, &SOLEXA_STEP, cp, h)))
    messcrash ("-flagTranscriptEnds: missing ERR wiggle file %s\n", cp) ;
      
  for (isDown = 1 ; isDown >= 0 ; isDown--)
    {
      int k, ii, iMax ;
      Array w_f = isDown ? w_LF : w_LR ;
      Array w_r = isDown ? w_RF : w_RR ;
      SXW *sf, *sr ;
      int uu, vv, u[3], v[3], oldx = 0 ;
      int typeL = isDown ? 0 : 2 ;
      int typeR = isDown ? 1 : 3 ;

      memset (u, 0, sizeof(u)) ;
      memset (v, 0, sizeof(v)) ;

      iMax = arrayMax (w_f) ;
      ii = arrayMax (w_r) ;
      if (iMax > ii) {   array (w_r, iMax - 1, SXW).m = 0 ; }
      else { iMax = ii ; array (w_f, iMax - 1, SXW).m = 0 ; }

      for (ii = k = 0,uu = vv = 0,  sf = arrp (w_f, ii, SXW), sr =   arrp (w_r, ii, SXW) ; ii < iMax ; ii++, k++, sf++, sr++)
	{
	  k = k % 3 ;
	  if (sf->n && sr->n && sf->a1 != sr->a1)
	    messcrash ("ii = %d,  sf->a1 = %d != sr->a1 = %d", ii, sf->a1, sr->a1) ;
 
	  uu -=  u[k] ; u[k] = sf->n ; uu += u[k] ;
	  vv -=  v[k] ; v[k] = sr->n ; vv += v[k] ;
          if (uu > min && uu > vv * minRatio)
	    {
	      SXX *sxx = arrayp (ww, nn++, SXX) ;
	      sxx->type = typeL + 1 ; sxx->m = uu ; sxx->n = vv ; sxx->a1 = oldx ; sxx->a2 = oldx ;
	    }
          if (vv > min && vv > uu * minRatio)
	    {
	      SXX *sxx = arrayp (ww, nn++, SXX) ;
	      sxx->type = typeR + 1 ; sxx->m = vv ; sxx->n = uu ; sxx->a1 = oldx ; sxx->a2 = oldx ;
	    }
	  oldx = sf->a1 ;
	}
    }

  arraySort (ww, sxxOrder) ; /* sort by a1 */
  /* export the results */
  if (1)
    {
      const char *types[] = { "ELF", "ERF", "ELR", "ERR" } ;
      ACEOUT ao = aceOutCreate (sx->outFileName, ".endFlags.ace", sx->gzo, h) ;
      int type, ii, iMax = arrayMax (ww) ;
      SXX *sxx, *old = 0 ;
      aceOutf (ao, "Sequence %s\n", sx->sxxChromosome ? sx->sxxChromosome : "unknownChrom") ;
      for (type = 1 ; type < 5 ; type ++) 
	for (ii = 0, old = 0, sxx = arrp (ww, ii, SXX) ; ii < iMax ; ii++, sxx++)
	  {
	    if (sxx->type != type) continue ;
	    if (old && sxx->a1 == old->a1 + 10)
	      { sxx->m += old->m ; sxx->n += old->n ; sxx->a2 = sxx->a1 ; sxx->a1 = old->a1 ; old->type = 0 ; } 
	    old = sxx ;
	  }

      for (ii = 0, sxx = arrp (ww, ii, SXX) ; ii < iMax ; ii++, sxx++)
	if (sxx->type > 2)
	  { int a = sxx->a1 ; sxx->a1 = sxx->a2 ; sxx->a2 = a ; }

     for (ii = 0, sxx = arrp (ww, ii, SXX) ; ii < iMax ; ii++, sxx++)
	if (sxx->type)
	  aceOutf (ao, "Feature %s\t%d\t%d\t%d\t%s_%.3f_%d\n", types[sxx->type - 1], sxx->a1, sxx->a2, sxx->m, types[sxx->type - 1], sxx->n/(1.0*sxx->m), sxx->n) ;

      aceOutf (ao, "\n") ;
    }

  fprintf (stderr, "sxFlagTranscriptEnds exported %d flags\n", nn) ;

  ac_free (h) ;
} /*  sxFlagTranscriptEnds */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (void)
{
  fprintf (stderr, 
	   "// Usage: geneelements \n"
	   "// Example:  geneelements -cmp -tissue brain\n"  
	   "//   -title text : a title, reexported after the date in the output\n"
	   "//   -i fileName : optional: .hits input file, defaults to stdin\n"
	   "//   -o fileName : optional: output file, defaults to stdout\n"
	   "//   -gzi : gunzip the input file on the fly\n"
	   "//   -cmp : compare different kind of tag lists, hard coded\n"
	   "//   -wordProfile : distribution of (position word) i.e. error profiles\n"
	   "// -ignore fileName : name of a file of ignored tags, one name per line\n"
	   "//           any tags listed in this file will be ignored when reading the tags file\n"
	   "//           This is very useful to perform a recurvive search for primers in the profile\n"
	   
	   "//   -prefix count_err_file : count_err_file is tab delimited probe/err, stdin is fasta, \n"
	   "//           output gives for each 6bp prefix the total number of sequences aligned with n errors\n"
	   "//   -prefixq count_err_file : count homopolymer at either end, fastq file\n"
	   "//   -manip : ILM,LifeTech...\n"
	   "//   -tissue : Brain, UHR...\n"
	   "//   -lane : 1,2..., manip.tissue.lane is used to locate hits files and fasta files\n"
	   "//   -exp [ -plot | -count | -venn | -intron]: expression analysis, std input is an ace file\n"
	   "//   -saturation nn, count how many targets appear in the first k/total element os a list of a 1 column list of targets\n"
	   "//   -support -plot | -count | -venn | -intron]: expression analysis, std input is an ace file\n"
	   "//   -nmWiggle suffix : create a wiggle for all NMs ending like the suffix\n"
	   "// -cosmidWiggle target [-stranded] [-delta int] [-hitFile hits] -trackName name (-wiggleDir name |  -feature name)\n"
	   "//       Read stdin or -hitFile, export a wiggle per target,\n"
	   "//       target has 4 columns: target map a1 a2, target extends on map from  a1 to a2,\n"
	   "//       hits has 3 or 4 columns: map a1 a2 [fastc tag_name], matching the appelations in target\n"
	   "//          the optional 4th column is a multiplier\n"
	   "//       period: defaults to 0\n"
	   "//          if period==0: a value is exported each time there is a change in coverage\n"
	   "//          if period>=1: a value is exported at regular intervals delta\n"
           "//    -wiggleDir dirName:  the wiggle is exported in wiggle format in the externalDirectory,\n"
	   "//          the file is named trackName.slx, \n"
	   "//    -stranded, the data are exported in 2 separate files trackName.[fr].slx according to the sign of a1-a2,\n"
	   "//                           the corresponding aceFile is on stdout\n" 
	   "//    -feature name: the wiggle is exported as class:Sequence tagFeature 'name' coord coord number-of-tags trackName\n"
	   "//\n"
	   "// -geneDosage : dedicated code\n"
	   "// -exomeDosage : dedicated code\n"
	   "// -mapSDF sdfFilename  -spongeFile fileName : map opbjects relative to a set of spongeFile\n"
	   "//       the sdf and sponge file use the same format\n"
	   "//       maps the sdf objects hierarchically in the sponge hierarchical geometry\n"
	   "//       For example label introns donor/acceptors as coding/UTR/intronic/intergenic\n"
	   "// -sponge threshold -spongeFile fileName  [-sxxChromosome chrom] -wiggle f.bf[.gz]\n"
	   "//       threshold is an integer, say 12, below which positions in the wiggle file are skipped\n"
	   "//       sponge_file is a 6 column tab delimited file\n"
	   "//            object  segment  chromosome a1 a2 gene\n"
	   "//       object identifies a collection of segments, typically a transcript or  gene identifier, gene is the corresponding gene\n"
	   "//       segment is an sub-identifier, typically the ordinal number of each exon\n"
	   "//       chromosome, where the segment maps, should coincide with one of the targets of the wiggle file\n"
	   "//       a1, a2 1-based begin and end of the segment on the chromosome, a1>a2 for the reverse strand\n"
	   "//     Optionally one may specify a coma separated collection of sponge file names as for example\n"
	   "//        -spongeFile \"CDS:cdsFile,exon:exon_file,gene:gene_file,genome:genomeFile\"\n"
	   "//        in this case, the assigantion of the wiggle is recursive\n"
	   "//     Optionally restrict the sponge file to chromosome -sxxChromosome\n"
	   "//       f.bf a genomic wiggle exported by the wiggle program, in one of the possible formats: bf, bv, tabix  \n"
	   "// -intronSponge  threshold -mrnaRemapFile fileName [-bestHitFile fileName (defaults to stdin)]\n"
	   "//       Measure the support for all introns in the best hit file\n"
	   "//         we count the support of the donor and acceptor base, and traversing reads\n"
	   "//       mrnaRemapFile has 6 columns: mrna x1 x2 chromo a1 a2 (defining all exons)\n"
	   "//       bestHitFile is the file .bestg.gz describing score and alignemnt of all tags\n"
  
	   "//   -spongeWindowHisto -wiggle f.bf[.gz] [-width <int nn>] \n"
	   "//       read a wiggle file, gives the histo of tags per window of width nn bp [default 200]\n"
	   "//   -spongeWindowExclude fileName :  exclude the tags in these regions 'map\ta1\ta2'\n"
	   
	   "//\n"
	   "// -mrnaWiggle target [-hitFile file] [-delta int] -trackName name\n"
	   "//       Read stdin or -hitFile, export a wiggle per target,\n"
	   "//       Target has 1 column of mrna names, corresponding to col 1 of the hits file\n"
	   "//       hits has 3 columns: map a1 a2, matching the appelations in target\n"
	   "//       each tag is trimmed by delta bp (default 0) at both ends\n"
	   "//       A value is exported each time there is a change in coverage\n"
	   "//       the wiggle is exported in .ace format as class:mRNA tag:Solexa \n"
	   "//            i.e.: Wiggle 'trackName'   coord-on-mRNA     number-of-tags\n"
	   "//       hits (default stdin) is the output of \n"
	   "//           probeAlign -splice -seedOffet 1 -seedShift 50 ... | grep HIT (to find the approximate exons)\n"
	   "// -ventilate -run run_name -maxChromNameLength <int> \n"
	   "//       This step positions the files for the subsequent pipeline call to -newIntrons\n"
	   "//       The directory names are hard coded, sorry\n"
	   "//       Ventilate the tmp/PHITS_genome/$run_name/*overhangs.gz files into tmp/[OR,pA,SL]/$chrom/$run_name.clean.gz\n"
	   "//        chromosome names longer than maxChromNameLength (if non zero) are skipped\n"
	   "// -newIntrons -hitFile hits [-stranded] -seedLength n1 -overhangLength n2 -minimalSupport ns\n"
	   "//             [-non_gt_ag nother]  [-transsplicing n3]\n"
	   "//             [-intronMinLength lMin]  [-intronMaxLength lMax]  [-minX x1 -maxX x2]\n"
	   "//       hits (default stdin) is the output of -besthits option of this program \n"
	   "//           probeAlign -splice -seedOffet 1 -seedShift 50 ... | grep OR  (to find the donors acceptors)\n"
	   "//       n1 {default=12] is the required number of bases matching the genome\n"
	   "//       n2 [default= 8] is the required number of bases overhanging\n"
	   "//       ns [default=3]  is the minimal number of tag support at each the boundarybases overhanging\n"
	   "//       lMin [default=30]  minimum intron length to be reported \n"
	   "//       lMax [default=10,000]  maximum intron length to be reported \n"
	   "//       transsplicing [default FALSE], if set report transplicing events (please set high n1,n2,n3 thresholds)\n"
	   "//       nother [default=0], if set and support >= nother, report introns with non gt_ag or gc_ag boundaries\n"
	   "//       x1,x2 : if provided, hits outside [x1,x2] are rejected \n"
	   "// -newDoubleIntrons   -hitFile hits [-stranded] -t chrom.fasta \n"
	   "//                     -sxxNewIntronsFileName f.introns  -sxxDeUnoIntronsFileName f.deUno\n"
	   "//                     -seedLength n1 -overhangLength n2 -minimalSupport ns \n"
	   "//                     [-intronMinLength lMin]  [-intronMaxLength lMax] [-minX x1 -maxX x2]\n"
	   "//       Search for intron-exon-intron segments (or SL or pA)\n"
	   "//       hits (default stdin) is the output of -besthits option of this program \n"
	   "//       chrom.fasta is the target fasta file (chromosome or transcriptome) previously provided to probealign\n"
	   "//       f.introns a set of introns exported previoulsy by -newIntrons\n" 
	   "//       f.deUno is the .introns.gz file exported by clipalign\n"
	   "//       chr a chromosome sequence name as listed in f.slx and f.introns\n"
	   "//       n1 {default=12] minimal number of bases of the exon overhang matching the distal intron feet\n"
	   "//       n2 [default= 20] maximal number of bases reported on the overhanging segments\n"
	   "//       ns [default=1]  is the minimal number of tag support at each the boundarybases overhanging\n"
	   "//       lMin [default=30]  minimum intron length to be considered \n"
	   "//       lMax [default=10,000]  maximum intron length to be considered \n"
	   "//       x1,x2 : if provided, hits outside [x1,x2] are rejected \n"
	   "//       stranded : (use for stranded hits) only consider hits and introns from the same strand\n"
	   "// -newExons -wiggle_f f.BF[.gz] -wiggle_r r.BF[.gz] -wiggle_ns ns.bf[.gz]\n"
	   "//           -sxxNewIntronsFileName f.introns -sxxChromosome chr -minimalSupport ns \n"
	   "//       Search, on chromoosome chr, for exons bordering a candidate intron\n"
	   "//              supported at least ns times in the non stranded wiggle\n"
	   "//       f.bf a genomic wiggle exported the wiggle program, in one of the possible formats: bf, bv, tabix  \n"
	   "//       f.introns a set of introns exported previoulsy by -newIntrons\n" 
	   "//       chrom a chromosome sequence name as listed in f.bf and f.introns\n"
	   "// -strand -antistrand \n"
	   "//       search for new exons or intron only on the strand or the anti strand\n"
	   "// -flagTranscriptEnds -wiggle name -run run_name  -minimalSupport ns -sxxChromosome chrom\n"
	   "//       Open the wiggles files called name.ELF/ELR/ERF/ERR.BF.gz\n"
	   "//       Export all position where for the F and R starnd separately, L > 4R or R > 4L and L > ns\n"
	   "//       summing over a 30 bp moving window\n"
	   "//       L/R means left/right terminator and F/R forward reverse starnd of the genome\n"
	   );
  exit (1) ;
} /* usage */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  FILE *f = 0 ;
  int outlevel = 0 ;
  SX sx ;
  AC_HANDLE h = 0 ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&sx, 0, sizeof (SX)) ;
  sx.h = h ;

#ifdef MEM_TEST
  if (0) /* 1 byte wiggle compressor
	  * store any number up to 102477 in 1 byte with 1% precision
	  */
    {
      oneByteTest () ;
      return 0 ;
    }

  {
    mysize_t gb = 1<<30 ;
    int ii, jj ;
    void * v ;
    int imax = argc >= 2 ? atoi(argv[1]) : 0 ;

    fprintf (stderr, "imax=%d\n", imax) ;

    for (ii= 0, jj = 1; ii < imax ; jj *=1)
      {
	v = messalloc (jj*gb) ;
	ii += jj ;
	if (v)
	  fprintf (stderr, "Allocated %d Gb\n", ii) ;
	else
	  {
	    fprintf (stderr, "crashed while allocating another %d Gb\n", jj) ;	  
	    return 1 ;
	  }
      }
  }
  return 0 ;
#endif

  if (0)
    {
      showDBLI (0) ; /* for compiler happiness */
      showSponge (0) ;
      sxPowerTable (0) ;
    }

  /* optional temple argument */
  if (getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage () ;

  getCmdLineOption (&argc, argv, "-title", &(sx.title)) ;
  sx.cmp = getCmdLineOption (&argc, argv, "-cmp", 0) ;

  getCmdLineInt (&argc, argv, "-prefix", &sx.prefix) ;
  getCmdLineInt (&argc, argv, "-prefixq", &sx.prefixq) ;
  getCmdLineOption (&argc, argv, "-manip", &sx.manip) ;
  getCmdLineOption (&argc, argv, "-tissue", &sx.tissue) ;
  getCmdLineOption (&argc, argv, "-lane", &sx.lane) ;
  getCmdLineOption (&argc, argv, "-run", &sx.run) ;
  getCmdLineOption (&argc, argv, "-group", &sx.group) ;
  getCmdLineOption (&argc, argv, "-t", &sx.target) ;
  getCmdLineOption (&argc, argv, "-support", &sx.support) ;
  getCmdLineInt (&argc, argv, "-maxChromNameLength", &sx.maxChromNameLength) ;
  sx.exp = getCmdLineOption (&argc, argv, "-exp", 0) ;
  sx.ventilate = getCmdLineOption (&argc, argv, "-ventilate", 0) ;
  getCmdLineInt (&argc, argv, "-saturation", &sx.saturation) ;
  sx.count = getCmdLineOption (&argc, argv, "-count", 0) ;
  sx.nPlatform = getCmdLineOption (&argc, argv, "-nPlatform", 0) ;
  sx.hammer = getCmdLineOption (&argc, argv, "-hammer", 0) ;
  sx.hits2composite = getCmdLineOption (&argc, argv, "-hits2composite", 0) ;
  sx.plot = getCmdLineOption (&argc, argv, "-plot", 0) ;
  sx.venn = getCmdLineOption (&argc, argv, "-venn", 0) ;
  sx.wordProfile = getCmdLineOption (&argc, argv, "-wordProfile", 0) ;
  sx.intron = getCmdLineOption (&argc, argv, "-intron", 0) ;
  sx.geneDosage = getCmdLineOption (&argc, argv, "-geneDosage", 0) ;
  sx.exomeDosage = getCmdLineOption (&argc, argv, "-exomeDosage", 0) ;
  sx.variantDosage = getCmdLineOption (&argc, argv, "-variantDosage", 0) ;
  getCmdLineOption (&argc, argv, "-tissue", &(sx.tissue));
  getCmdLineOption (&argc, argv, "-ignore", &(sx.ignoredProbeFileName)) ;

  getCmdLineOption (&argc, argv, "-hitFile", &(sx.hitFile));

  getCmdLineOption (&argc, argv, "-nmWiggle", &sx.nmWiggle) ;
  sx.mrnaWiggle = getCmdLineOption (&argc, argv, "-mrnaWiggle", 0);
  sx.stranded = getCmdLineOption (&argc, argv, "-stranded", 0);
  sx.strand = getCmdLineOption (&argc, argv, "-strand", 0);
  sx.antistrand = getCmdLineOption (&argc, argv, "-antistrand", 0);
  sx.unique = getCmdLineOption (&argc, argv, "-u", 0);
  getCmdLineOption (&argc, argv, "-cosmidWiggle", &(sx.cosmidWiggle));
  if (getCmdLineOption (&argc, argv, "-strategy", &(sx.strategy)) &&
      ! strcasecmp (sx.strategy, "RNA_seq"))
    sx.RNA_seq = TRUE ;
  
  getCmdLineOption (&argc, argv, "-wiggleDir", &sx.wiggleDir) ;
  getCmdLineOption (&argc, argv, "-feature", &sx.wiggleFeature) ;
  getCmdLineOption (&argc, argv, "-trackName", &(sx.trackName));
  getCmdLineInt (&argc, argv, "-delta", &(sx.delta));
  getCmdLineInt (&argc, argv, "-period", &(sx.period));
  sx.newIntrons =  getCmdLineOption (&argc, argv, "-newIntrons", 0) ;
  getCmdLineInt (&argc, argv, "-non_gt_ag", &(sx.non_gt_ag)) ;
  getCmdLineInt (&argc, argv, "-transsplicing", &(sx.transsplicing)) ;
  sx.newDoubleIntrons =  getCmdLineOption (&argc, argv, "-newDoubleIntrons", 0) ;
  sx.newExons =  getCmdLineOption (&argc, argv, "-newExons", 0) ;
  getCmdLineOption (&argc, argv, "-wiggle", &sx.wiggleFileName) ;
  getCmdLineOption (&argc, argv, "-wiggle_ns", &sx.wiggleFileName_ns) ;
  getCmdLineOption (&argc, argv, "-wiggle_f", &sx.wiggleFileName_f) ;
  getCmdLineOption (&argc, argv, "-wiggle_r", &sx.wiggleFileName_r) ;
  getCmdLineOption (&argc, argv, "-sxxNewIntronsFileName", &sx.sxxNewIntronsFileName) ;
  getCmdLineOption (&argc, argv, "-sxxDeUnoIntronsFileName", &sx.sxxDeUnoIntronsFileName) ;
  /* 2016_01_16, i disable this command because the parser for sxxNewIntronsFileName in sxNewExonsGetIntron
     has changed be cause the format of the file EHITS.Test/14/introns.de_duo is completely modified
     also knownIntron is no longer used in scripts
     if we want to reuse it the format should be modified

  getCmdLineOption (&argc, argv, "-sxxKnownIntronsFileName", &sx.sxxKnownIntronsFileName) ;
  */
  getCmdLineOption (&argc, argv, "-sxxChromosome", &sx.sxxChromosome) ;

  getCmdLineInt (&argc, argv, "-intronMinLength", &sx.intronMinLength) ;
  getCmdLineInt (&argc, argv, "-intronMaxLength", &sx.intronMaxLength) ;
  getCmdLineInt (&argc, argv, "-seedLength", &sx.seedLength) ;
  getCmdLineInt (&argc, argv, "-overhangLength", &sx.overhangLength) ;
  getCmdLineInt (&argc, argv, "-minimalSupport", &sx.minimalSupport) ;
  sx.flagTranscriptEnds = getCmdLineOption(&argc, argv, "-flagTranscriptEnds", 0) ; 
  sx.exonLevel = sx.minimalSupport ? sx.minimalSupport : 10 ;
  sx.minimalIntronSupport = 1 ;
  if (sx.minimalSupport > 10) sx.minimalIntronSupport = sx.minimalSupport/10 ;
  getCmdLineInt(&argc, argv, "-minimalIntronSupport", &sx.minimalIntronSupport) ;
  if (sx.non_gt_ag && sx.non_gt_ag < 2 * sx.minimalIntronSupport)
    sx.non_gt_ag =  2 * sx.minimalIntronSupport ;
  if (sx.non_gt_ag && sx.non_gt_ag < 8) 
    sx.non_gt_ag =  8 ;
  getCmdLineInt(&argc, argv, "-minX", &sx.minX) ;
  getCmdLineInt(&argc, argv, "-maxX", &sx.maxX) ;
  getCmdLineInt(&argc, argv, "-saturation", &sx.saturation) ;
  getCmdLineInt(&argc, argv, "-sponge", &sx.sponge) ;
  getCmdLineInt(&argc, argv, "-width", &sx.width) ;
  getCmdLineOption(&argc, argv, "-spongeFile", &(sx.spongeFileName)) ;
  getCmdLineOption(&argc, argv, "-mapSDF", &(sx.sdfFileName)) ;
  getCmdLineInt(&argc, argv, "-intronSponge", &sx.intronSponge) ;
  sx.intronSpongeMerge = getCmdLineOption(&argc, argv, "-intronSpongeMerge", 0) ;
  getCmdLineOption (&argc, argv, "-mrnaRemapFile", &sx.mrnaRemapFileName) ;
  getCmdLineOption (&argc, argv, "-bestHitFile", &sx.bestHitFileName) ;
  getCmdLineOption (&argc, argv, "-genomeMask", &sx.genomeMaskFileName) ;
  sx.spongeWindow = getCmdLineOption (&argc, argv, "-spongeWindow", 0) ;
  getCmdLineOption (&argc, argv, "-spongeWindowExclude", &sx.spongeWindowExcludeFileName) ;
  getCmdLineOption (&argc, argv, "-split_mRNAs", &sx.splitMrnasFileName) ;

  fprintf (stderr, "// start: %s\n", timeShowNow()) ;

  if (getCmdLineOption (&argc, argv, "-o", &(sx.outFileName)))
    {
      if (sx.outFileName && ! sx.exomeDosage && ! sx.variantDosage && ! sx.geneDosage)  /* geneDosage is using aceOutCreate () */
	{
	  FILE *f = filopen (sx.outFileName, 0, "w") ;
	  if (f)
	    outlevel = freeOutSetFile (f) ;	
	}
    }
  if (!outlevel)
    outlevel = freeOutSetFile (stdout) ;	
  
  if (argc > 1)
    messcrash ("Unkwong argument %s, try -help ", argv[1]) ;

  if (sx.title) 
   {
     freeOutf ("// %s\n", timeShowNow()) ;
     freeOutf ("// %s\n", sx.title) ;
   }
  if (sx.ignoredProbeFileName)
    {
      int level = 0 ;
      char *cp ;
      FILE *f = filopen (sx.ignoredProbeFileName, 0, "r") ;
      if (f)
	{
	  sx.ignoredProbeDict = dictHandleCreate (10000, h) ;
	  level = freesetfile (f, 0) ;
	  while (freecard (level))
	    if ((cp = freeword ()))
	      {
		if (!strcasecmp (cp, "probe"))
		  cp = freeword () ;
		if (cp)
		  dictAdd (sx.ignoredProbeDict, cp, 0) ;
	      }
	  freeclose (level) ; /* will close f */
	}
    }
  if (sx.cmp)
    sxCompareCoverage (&sx) ; /* hard coded for SEQ-C 2008_12_14 */      
  else if (sx.wordProfile)
    sxWordProfile (&sx) ;   /* create a profile of positioned words */
  else if (sx.prefix || sx.prefixq) 
    sxPrefixDistrib (&sx) ;   /*  */
  else if (sx.nmWiggle)
    sxNmWiggle (&sx) ;   /* create a wiggle for all the NMs */
  else if (sx.mrnaWiggle || sx.cosmidWiggle)
    sxTargetWiggle (&sx) ;   /* create a wiggle for all the NMs */
  else if (sx.hits2composite) 
    sxHits2composite (&sx) ;   /* create a wiggle for all the NMs */
  else if (sx.exp)
    sxExp (&sx) ;   /* compare genes etc */
  else if (sx.saturation)
    sxSaturation (&sx) ;   /* compare genes etc */
  else if (sx.ventilate)
    sxVentilate (&sx) ;   /*  position the files for the subsequent pipeline call to -newIntrons */
  else if (sx.newIntrons)
    sxNewIntrons (&sx) ;   /* search for new introns in the output of probealign -slide */
  else if (sx.newDoubleIntrons)
    sxNewDoubleIntrons (&sx) ;   /* given hits and new introns, search for intron-exon-intron segments */
  else if (sx.newExons)
    sxNewExons (&sx) ;   /* search for XE XI XA using wiggles and a candidate file */
  else if (sx.sdfFileName)
    sxMapSdf (&sx) ;
  else if (sx.sponge)
    sxSponge (&sx) ;   /* count the support for exons introns */
  else if (sx.intronSponge)  /* 2013_09_02: possibly obsolete relative to MAGIC ii2a/ii2b */
    sxIntronSponge (&sx) ;   /* count the support for exons introns */
  else if (sx.intronSpongeMerge)
    sxIntronSpongeMerge (&sx) ;   /* merge the support for exons introns */
  else if (sx.spongeWindow && sx.wiggleFileName)
    sxSpongeWindowHisto (&sx) ;
  else if (sx.exomeDosage || sx.variantDosage || sx.geneDosage)
    sxGeneDosage (&sx) ;
  else if (sx.flagTranscriptEnds)
    sxFlagTranscriptEnds (&sx) ;
  else 
    usage () ;

  if (outlevel)
    freeOutClose (outlevel) ;
  if (f) filclose (f) ;
  ac_free (sx.h) ;

  {
    int mx ;
    messAllocMaxStatus (&mx) ; 
    fprintf (stderr, "// %s done: max memory %d Mb\n", timeShowNow(), mx) ;
  }

  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/*

bin/geneelements -sponge 1 -spongeFile ,Genome:tmp/METADATA/genome.ns.sponge -sxxChromosome 14 -wiggle tmp/WIGGLERUN/Rrn386/14/R.chrom.frns.u.BF.gz | grep Total


bin/geneelements -sponge 1 -spongeFile ,Discovery_genes:tmp/METADATA/magic.ns.gene.sponge,Genome:tmp/METADATA/genome.ns.sponge -sxxChromosome 14 -wiggle tmp/WIGGLERUN/Rrn386/14/R.chrom.frns.u.BF.gz | grep Total


bin/geneelements -sponge 1 -spongeFile ,RefSeq_genes:tmp/METADATA/RefSeq.ns.gene.sponge,AceView_genes:tmp/METADATA/av.ns.gene.sponge,Discovery_genes:tmp/METADATA/magic.ns.gene.sponge,Genome:tmp/METADATA/genome.ns.sponge -sxxChromosome 14 -wiggle tmp/WIGGLERUN/Rrn386/14/R.chrom.frns.u.BF.gz | grep Total

bin/geneelements -sponge 1 -spongeFile ,RefSeq_genes:tmp/METADATA/RefSeq.ns.gene.sponge,AceView_genes:tmp/METADATA/av.ns.gene.sponge,Genome:tmp/METADATA/genome.ns.sponge -sxxChromosome 14 -wiggle tmp/WIGGLERUN/Rrn386/14/R.chrom.frns.u.BF.gz | grep Total



bin/geneelements -sponge 1 -spongeFile ,Discovery_CDS:tmp/METADATA/magic_cds.ns.sponge,RefSeq_transcripts:tmp/METADATA/RefSeq.ns.sponge,AceView_transcripts:tmp/METADATA/av.ns.sponge,Discovery_transcripts:tmp/METADATA/magic.ns.sponge,RefSeq_genes:tmp/METADATA/RefSeq.ns.gene.sponge,AceView_genes:tmp/METADATA/av.ns.gene.sponge,Discovery_genes:tmp/METADATA/magic.ns.gene.sponge,Genome:tmp/METADATA/genome.ns.sponge -sxxChromosome 14 -wiggle tmp/WIGGLERUN/Rrn386/14/R.chrom.frns.u.BF.gz | grep Total


bin/geneelements -sponge 1 -spongeFile ,Discovery_CDS:tmp/METADATA/magic_cds.ns.sponge,RefSeq_transcripts:tmp/METADATA/RefSeq.ns.sponge,AceView_transcripts:tmp/METADATA/av.ns.sponge,Discovery_transcripts:tmp/METADATA/magic.ns.sponge,RefSeq_genes:tmp/METADATA/RefSeq.ns.gene.sponge,AceView_genes:tmp/METADATA/av.ns.gene.sponge,Discovery_genes:tmp/METADATA/magic.ns.gene.sponge,Genome:tmp/METADATA/genome.ns.sponge -sxxChromosome 14 -wiggle tmp/WIGGLERUN/Rrn386/14/R.chrom.frns.u.BF.gz | grep Total

*/
