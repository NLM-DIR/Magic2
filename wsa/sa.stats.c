/*
 * sa.stats.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module registers and exports detailled run statistics in tsf format
 */

#include "sa.h"

/**************************************************************/

void saCpuStatRegister (const char *nam, int agent, Array cpuStats, clock_t t1, clock_t t2, long int n)
{
  int s = arrayMax (cpuStats) ;
  CpuSTAT *bs = arrayp (cpuStats, s, CpuSTAT) ;
  bs->nB++ ;
  strncpy (bs->nam, nam, 30) ;
  bs->agent = agent ;
  bs->n += n ;
  bs->tA += t2 - t1 ;  /* active time */
  return ;
} /* saCpuStatRegister */

/**************************************************************/

void saCpuStatCumulate (Array aa, Array a)
{
  int i = arrayMax (aa) ;
  int j = arrayMax (a) ;
  for (int k = 0 ; k < j ; k++)
    array (aa, i + k, CpuSTAT) = array (a, k, CpuSTAT) ;
  return ;
} /* saCpuStatCumulate */

/**************************************************************/

static int cpuStatsOrder (const void *va, const void *vb)
{
  const CpuSTAT *up = va, *vp = vb ;
  int n = 0 ;

  n = strcmp (up->nam, vp->nam) ; if (n) return n ;
  n = up->agent - vp->agent ; if (n) return n ; 
  return 0 ;
} /* cpuStatsOrder */

/**************************************************************/

void saCpuStatExport (const PP *pp, Array stats)
{
  int i, j, iMax = arrayMax (stats) ;
  CpuSTAT *bt, *bs = arrayp (stats, 0, CpuSTAT) ;

  arraySort (stats, cpuStatsOrder) ;
  printf ("\n# Action\tAgent\tnB\tn\ttime (s)") ;
  if (1)
    {  /* details per agent, only useful for debugging and optimizing */
      for (i = 0 ; i < iMax ; i++, bs++)
	{
	  if (! bs->nB)
	    continue ;
	  for (j = i + 1, bt = bs + 1 ; j < iMax ; j++, bt++)
	    {
	      if (strcmp (bs->nam, bt->nam) || bs->agent != bt->agent)
		break ;
	      bs->nB += bt->nB ; bs->n += bt->n, bs->tA += bt->tA ; bt->nB = 0 ;
	    }
	  printf ("\n%s\t%d\t%d\t%ld\t%.0f", bs->nam, bs->agent, bs->nB, bs->n, bs->tA * 1.0 / CLOCKS_PER_SEC) ;
	}
      printf ("\n") ;
    }
  bs = arrayp (stats, 0, CpuSTAT) ;
  for (i = 0 ; i < iMax ; i++, bs++)
    {
      if (! bs->nB)
	continue ;
      for (j = i + 1, bt = bs + 1 ; j < iMax ; j++, bt++)
	{
	  if (strcmp (bs->nam, bt->nam))
	    break ;
	  if (bt->nB)
	    { bs->nB += bt->nB ; bs->n += bt->n, bs->tA += bt->tA ; bt->nB = 0 ; }
	}
      printf ("\n%s\t%d\t%d\t%ld\t%.0f", bs->nam, bs->agent, bs->nB, bs->n, bs->tA  * 1.0 / CLOCKS_PER_SEC) ;
    }
  printf ("\n") ;
  return ;
} /* saCpuStatExport */

/*************************************************************************************/

static int confirmedIntronsCountSites (const PP *pp, int run)
{
  int nn = 0, ii, iMax = arrayMax (pp->confirmedIntrons) ;
  INTRON *zp = iMax ? arrp (pp->confirmedIntrons, 0, INTRON) : 0 ;
  
  if (run)
    {
      for (ii = 0 ; ii < iMax ; ii++, zp++)
	if (zp->run == run)
	  {
	    int min = 3 ;
	    if (! zp->feet[0])
	      continue ;
	    else if (!strcmp (zp->feet, "gt_ag"))
	      min = 1 ;
	    else if (!strcmp (zp->feet, "gc_ag"))
	      min = 2 ;
	    if (zp->n + zp->nR >= min)
	      nn++ ;
	  }
    }
  else
    {
      int n = 0, chrom = 0, a1 = 0, a2 = 0, min = 3 ;
      for (ii = 0 ; ii < iMax ; ii++, zp++)
	{
	  if (! zp->feet[0])
	    continue ;
	  if (zp->chrom != chrom || zp->a1 != a1 || zp->a2 != a2)
	    {
	      if (n >= min)
		nn++ ;
	      if (!strcmp (zp->feet, "gt_ag"))
		min = 1 ;
	      else if (!strcmp (zp->feet, "gc_ag"))
		min = 2 ;
	      chrom = zp->chrom ; a1 = zp->a1 ; a2 = zp->a2 ;
	      n = 0 ;
	    }
	  n += zp->n + zp->nR ;
	}
      if (n >= min)
	nn++ ;
    }

  return nn ;
} /* confirmedIntronsCountSites */

/*************************************************************************************/

static int confirmedPolyAsCountSites (const PP *pp, int run)
{
  int nn = 0, ii, iMax = arrayMax (pp->confirmedPolyAs) ;
  POLYA *zp = iMax ? arrp (pp->confirmedPolyAs, 0, POLYA) : 0 ;

  if (run)
    {
      for (ii = 0 ; ii < iMax ; ii++, zp++)
	if (zp->run == run)
	  if (zp->n >= 3)
	    nn++ ;
    }
  else
    {
      int n = 0, chrom = 0, a1 = 0 ;
      for (ii = 0 ; ii < iMax ; ii++, zp++)
	{
	  if (zp->chrom != chrom || zp->a1 != a1)
	    {
	      if (n >= 3)
		nn++ ;
	      chrom = zp->chrom ; a1 = zp->a1 ; n = 0 ;
	    }
	  n += zp->n ;
	}
      if (n >= 3)
	nn++ ;
    }
  return nn ;
} /* confirmedPolyAsCountSites */

/*************************************************************************************/
/*************************************************************************************/

static void s2gSamStatsExports (const PP *pp, Array runStats)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *METHOD = pp->method ? pp->method : "01_SortAlign" ;
  Array aa = runStats ;
  RunSTAT *s0 = arrayp (aa, 0, RunSTAT) ;
  ACEOUT ao = aceOutCreate (pp->outFileName, ".s2g.samStats", 0, h) ;	
  const char *run = pp->runName ? pp->runName : (dictMax (pp->runDict) == 1 ? dictName (pp->runDict, 1) : "AllRuns") ;
  
  long int nRawReads = pp->nRawReads ? pp->nRawReads : s0->p.nReads ;
  long int nRawBases = pp->nRawBases ? pp->nRawBases : s0->p.nBase1 + s0->p.nBase2 ;
  
  aceOutf (ao, "%s\t%s\tPairs\t%ld\n", run, METHOD, s0->p.nPairs) ;
  aceOutf (ao, "%s\t%s\tAligned_pairs\t%ld\t%.2f%%\n", run, METHOD, s0->nPairsAligned, (100.0 * s0->nPairsAligned)/(nRawReads/2 + .000001)) ;
  aceOutf (ao, "%s\t%s\tnCompatible_pairs\t%ld\t%.2f%%\n", run, METHOD, s0->nCompatiblePairs, (100.0 * s0->nCompatiblePairs)/(s0->nPairsAligned + .000001)) ;
  aceOutf (ao, "%s\t%s\tnCircle_pairs\t%ld\t%.2f%%\n", run, METHOD, s0->nCirclePairs, (100.0 * s0->nCirclePairs)/(s0->nPairsAligned + .000001)) ;
  aceOutf (ao, "%s\t%s\tNon_compatible_pairs\t%ld\t%.2f%%\n", run, METHOD, s0->nIncompatiblePairs, (100.0 * s0->nIncompatiblePairs)/(s0->nPairsAligned + .000001)) ;
  

  aceOutf (ao, "%s\t%s\tReads\t%ld\n", run, METHOD, s0->p.nReads) ;  
  aceOutf (ao, "\n%s\t%s\tRawReads\t%ld\n", run, METHOD, nRawReads) ;
  if (s0->nClippedPolyA) aceOutf (ao, "\n%s\t%s\tClipped_PolyA\t%ld\n", run, METHOD, s0->nClippedPolyA) ;
  if (s0->nClippedPolyT) aceOutf (ao, "\n%s\t%s\tClipped_PolyT\t%ld\n", run, METHOD, s0->nClippedPolyT) ;
  if (s0->nClippedSL) aceOutf (ao, "\n%s\t%s\tClipped_SL\t%ld\n", run, METHOD, s0->nClippedSL) ;
  int nClipped = s0->nClippedAdaptor1L + s0->nClippedAdaptor2L + s0->nClippedAdaptor1R + s0->nClippedAdaptor2R ;
  if (nClipped) aceOutf (ao, "\n%s\t%s\tClipped_Adaptors\t%ld\n", run, METHOD, nClipped) ;
  aceOutf (ao, "%s\t%s\tPerfect_reads\t%ld\t%.2f%%\n", run, METHOD, s0->nPerfectReads, (100.0 * s0->nPerfectReads)/(nRawReads + .000001)) ;
  aceOutf (ao, "%s\t%s\tAlignedReads\t%ld\t%.2f%%\n", run, METHOD, s0->nMultiAligned[0], (100.0 * s0->nMultiAligned[0])/(nRawReads + .000001)) ;
  
  aceOutf (ao, "\n%s\t%s\tnAlignments\t%ld\t%.2f per aligned read\n", run, METHOD, s0->nAlignments, (1.0 * s0->nAlignments)/(s0->nMultiAligned[0] + .000001)) ;
  aceOutf (ao, "%s\t%s\tnMultiAligned\t%ld\t%.2f%%\n", run, METHOD, s0->nMultiAligned[0] - s0->nMultiAligned[1], (100.0 * (s0->nMultiAligned[0] - s0->nMultiAligned[1]))/(s0->p.nReads + .000001)) ;

  aceOutf (ao, "%s\t%s\tnBases\t%ld\n", run, METHOD, s0->p.nBase1 + s0->p.nBase2) ;
  aceOutf (ao, "\n%s\t%s\tRawBases\t%ld\n", run, METHOD, nRawBases) ;
  aceOutf (ao, "%s\t%s\tnAlignedBases\t%ld\t%.2f%%\n", run, METHOD, s0->nBaseAligned1 + s0->nBaseAligned2, 100.0 * (s0->nBaseAligned1 + s0->nBaseAligned2) / (nRawBases + .000001)) ;
  aceOutf (ao, "%s\t%s\tnErrorEvents\t%ld\t%.6f%%\n", run, METHOD, s0->nErr, (100.0 * s0->nErr)/(s0->nBaseAligned1 + s0->nBaseAligned2 + 0.00000001)) ;
  aceOutf (ao, "%s\t%s\tnMismatches_and_InDels\t%ld\t%.6f%%\n", run, METHOD, s0->nMID, (100.0 * s0->nMID)/(s0->nBaseAligned1 + s0->nBaseAligned2 + 0.00000001)) ;

  aceOutf (ao, "\n%s\t%s\tnPolyA_sites\t%ld\n", run, METHOD, confirmedPolyAsCountSites (pp, 0)) ;
  aceOutf (ao, "\n%s\t%s\tSupported_introns\t%ld\n", run, METHOD, confirmedIntronsCountSites (pp, 0)) ;
  aceOutf (ao, "\n%s\t%s\tnDoubleIntrons\t%ld\n", run, METHOD, arrayMax (pp->doubleIntrons)) ;
  aceOutf (ao, "%s\t%s\tIntron_supports\t%ld\t%ld\t%ld\t%.3f\n", run, METHOD
	   , s0->gt_ag_Support + s0->ct_ac_Support 
	   , s0->gt_ag_Support
	   , s0->ct_ac_Support
	   , s0->intronStranding
	   ) ;
  
  long int nUnaligned = nRawReads - s0->nMultiAligned[0] ;
  aceOutf (ao, "\n%s\t%s\tnMultiAligned %d times\t%ld\t%.2f%%\n", run, METHOD, 0
	   , nUnaligned
	   , 100.0 * nUnaligned / (nRawReads + .000001)
	   ) ;
  for (int j = 1 ; j < 11 ; j++)
    if (s0->nMultiAligned[j])
      {
	aceOutf (ao, "%s\t%s\tnMultiAligned %d times\t%ld\t%.2f%%\n", run, METHOD
		 , j, s0->nMultiAligned[j]
		 , 100.0 * s0->nMultiAligned[j] / (nRawReads + .000001)
		 ) ;
      }
  long int verif = 0 ;
  for (int j = 1 ; j < 11 ; j++)
    verif += s0->nMultiAligned[j] ;
  aceOutf (ao, "nReads = %ld , sum of multiAli = %ld, verif = %ld\n", nRawReads, verif, nRawReads - verif) ; 

  ac_free (h) ;
  return ;
} /* s2gSamStatsExport */

/**************************************************************/
/**************************************************************/
/* cumulate int global runStats the content of bb->runStats */
void saRunStatsCumulate (int run, PP *pp, BB *bb)
{
  Array aa = pp->runStats ;
  RunSTAT *up = arrayp (aa, run, RunSTAT) ;
  RunSTAT *vp = &(bb->runStat) ;
			    
  up->p.run = run ;
  up->p.nPairs += vp->p.nPairs ;
  up->nCompatiblePairs += vp->nCompatiblePairs ;
  up->nIncompatiblePairs += vp->nIncompatiblePairs ;
  up->nCirclePairs += vp->nCirclePairs ;
  up->nOrphans += vp->nOrphans ;
  up->n2ChromsPairs += vp->n2ChromsPairs ;
  up->nPairsAligned += vp->nPairsAligned ;
  up->p.nReads += vp->p.nReads ;
  up->p.nBase1 += vp->p.nBase1 ;
  up->p.nBase2 += vp->p.nBase2 ;
  up->lowEntropy += vp->lowEntropy ;
  up->tooShort += vp->tooShort ;
  up->lowEntropyBases += vp->lowEntropyBases ;
  up->tooShortBases += vp->tooShortBases ;
  up->nSupportedIntrons += vp->nSupportedIntrons ;
  up->nIntronSupportPlus += vp->nIntronSupportPlus ;
  up->nIntronSupportMinus += vp->nIntronSupportMinus ;
  up->gt_ag_Support += vp->gt_ag_Support ;
  up->ct_ac_Support += vp->ct_ac_Support ;
  
  for (int i = 0 ; i < 5 ; i++)
    up->p.ATGCN[i] += vp->p.ATGCN[i] ;
  for (int i = 0 ; i < 5 * LETTERMAX; i++)
    up->p.letterProfile1[i] += vp->p.letterProfile1[i] ;
  for (int i = 0 ; i < 5 * LETTERMAX; i++)
    up->p.letterProfile2[i] += vp->p.letterProfile2[i] ;
  for (int i = 0 ; i < 11 ; i++)
    up->nMultiAligned[i] += vp->nMultiAligned[i] ;
  for (int i = 0 ; i < 256 ; i++)
    up->nReadsAlignedPerTargetClass[i] += vp->nReadsAlignedPerTargetClass[i] ;
  for (int i = 0 ; i < 256 ; i++)
    up->nBasesAlignedPerTargetClass[i] += vp->nBasesAlignedPerTargetClass[i] ;

  for (int i = 0 ; i < OVERHANGMAX ; i++)
    up->overhangL1[i] += vp->overhangL1[i] ;
  for (int i = 0 ; i < OVERHANGMAX ; i++)
    up->overhangL2[i] += vp->overhangL2[i] ;

  for (int i = 0 ; i < OVERHANGMAX ; i++)
    up->overhangR1[i] += vp->overhangR1[i] ;
  for (int i = 0 ; i < OVERHANGMAX ; i++)
    up->overhangR2[i] += vp->overhangR2[i] ;

  up->nPerfectReads += vp->nPerfectReads ;
  up->nAlignments += vp->nAlignments ;
  up->nBaseAligned1 += vp->nBaseAligned1 ;
  up->nBaseAligned2 += vp->nBaseAligned2 ;
  up->nSupportedIntrons = saSupportedIntrons (pp, run) ;
  if (up->p.minReadLength == 0 || up->p.minReadLength > vp->p.minReadLength)
    up->p.minReadLength = vp->p.minReadLength ;
  if (up->p.maxReadLength < vp->p.maxReadLength)
    up->p.maxReadLength = vp->p.maxReadLength ;
  up->nErr += vp->nErr ;
  up->nMID += vp->nMID ;
  
  up->nClippedPolyA += vp->nClippedPolyA ;
  up->nClippedPolyT += vp->nClippedPolyT ;
  for (int i = 0 ; i < SLMAX ; i++)
    up->nClippedSls[i] += vp->nClippedSls[i] ;

  up->nClippedAdaptor1L += vp->nClippedAdaptor1L ;
  up->nClippedAdaptor1R += vp->nClippedAdaptor1R ;
  up->nClippedAdaptor2L += vp->nClippedAdaptor2L ;
  up->nClippedAdaptor2R += vp->nClippedAdaptor2R ;
  
  if (vp->errors)
    for (int i = 0 ; i < arrayMax (vp->errors) ; i++)
      array (up->errors, i, int) += array (vp->errors, i, int) ;
  for (int i = 0 ; i < 256 ; i++)
    {
      up->GF[i] += vp->GF[i] ;
      up->GR[i] += vp->GR[i] ;
    }
  if (!up->p.lengthDistribution)
    {
      up->p.lengthDistribution = arrayHandleCopy (vp->p.lengthDistribution, pp->h) ;
    }
  else
    {
      for (int i = 0 ; i < arrayMax (vp->p.lengthDistribution) ; i++)
	array (up->p.lengthDistribution, i , long int) += array (vp->p.lengthDistribution, i , long int) ;
    }
  if (!up->insertLengthDistribution)
    {
      up->insertLengthDistribution = arrayHandleCopy (vp->insertLengthDistribution, pp->h) ;
    }
  else
    {
      for (int i = 0 ; i < arrayMax (vp->insertLengthDistribution) ; i++)
	array (up->insertLengthDistribution, i , long int) += array (vp->insertLengthDistribution, i , long int) ;
    }
  return ;
} /* saRunStatsCumulate */
 
/**************************************************************/

static void saLetterProfileExport (const PP *pp, int run, RunSTAT *up)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *runName = dictName (pp->runDict, run) ;
  char *title = hprintf (h, "%s.letterProfile.tsf", runName) ;
  ACEOUT ao = aceOutCreate (pp->outFileName, title, 0, h) ;
  aceOutDate (ao, "##", hprintf (h, "Run %s Raw Letter Profile (before alignment) limited to %d letters", runName, LETTERMAX)) ;
  long int nn ;
  long int *aa = up->p.letterProfile1 ;
  int i, j ;

  for (j = 0, nn = 0 ; j < 5 ; j++)
    nn += aa[j] ;
  if (nn)
    {
      aceOutf (ao, "# Run.f\tPosition\tiiiiiifffff\tAny\tA\tT\tG\tC\tN\t%%A\t%%T\t%%G\t%%C\t%%N\n") ;
      for (i = 0 ; i < LETTERMAX ; i++)
	{
	  int i5 = 5 * i ;
	  for (j = 0, nn = 0 ; j < 5 ; j++)
	    nn += aa[i5 + j] ;
	  if (! nn)
	    break ;
	  aceOutf (ao, "%s.f\t%d\tiiiiiifffff\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n"
		   , runName
		   , i + 1, nn, aa[i5 + 0], aa[i5 + 1], aa[i5 + 2], aa[i5 + 3], aa[i5 + 4]
		   , 100.0 * aa[i5 + 0]/nn, 100.0 * aa[i5 + 1]/nn, 100.0 * aa[i5 + 2]/nn, 100.0 * aa[i5 + 3]/nn, 100.0 * aa[i5 + 4]/nn
		   ) ;
	}
    }

  aa = up->p.letterProfile2 ;
  for (j = 0, nn = 0 ; j < 5 ; j++)
    nn += aa[j] ;
  if (nn)
    {
      aceOutf (ao, "# Run.r\tPosition\tiiiiiifffff\tAny\tA\tT\tG\tC\tN\t%%A\t%%T\t%%G\t%%C\t%%N\n") ;
      for (i = 0 ; i < LETTERMAX ; i++)
	{
	  int i5 = 5 * i ;
	  for (j = 0, nn = 0 ; j < 5 ; j++)
	    nn += aa[i5 + j] ;
	  if (! nn)
	    break ;
	  aceOutf (ao, "%s.r\t%d\tiiiiiifffff\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n"
		   , runName
		   , i + 1, nn, aa[i5 + 0], aa[i5 + 1], aa[i5 + 2], aa[i5 + 3], aa[i5 + 4]
		   , 100.0 * aa[i5 + 0]/nn, 100.0 * aa[i5 + 1]/nn, 100.0 * aa[i5 + 2]/nn, 100.0 * aa[i5 + 3]/nn, 100.0 * aa[i5 + 4]/nn
		   ) ;
	}
    }

  aceOutf (ao, "\n\n") ;
  ac_free (h) ;
  return ;
} /* saLetterProfileExport */
 
/**************************************************************/

static BOOL readOneAdaptor (char adaptor[32], RunSTAT *up, BOOL isLeft, int pass, int *bestp, BOOL coded)
{
  int best = 0, best2, i, j, di, nn, nL = 0, good = 0 ;
  
  long int *aa = isLeft ? (pass ? up->overhangL2 : up->overhangL1) : (pass ? up->overhangR2 : up->overhangR1) ;
  /* find the best first letter */
  for (best = i = j = nn = 0 ; i < 4 ; i++)
    {
      int n = aa[150 * i + i] ;
      if (n > nn) { best = i ; nn = n ; }
    }
  /* cumulate in best the shifted histogram of the second letter of the best case */
  best2 = best ;
  di = 0 ;
  while (best2 == best && ++di < 6)
    {
      for (best2 = j = nn = 0 ; j < 5 ; j++)
	{
	  int n = aa[150 * best + 5 *di + j] ;
	  if (n > nn) { best2 = j ; nn = n ; }
	}
    }
  if (1) 
    { /* clean deconvoluted image */
      if (best2 != best)
	{ 
	  for (i = di ; i < 30 ; i++)
	    for (j = 0 ; j < 5 ; j++)
	      aa[150 * best + 5 * i + j] += aa[150 * best2 + 5 * (i-di) + j] ;
	  for (i = 0 ; i < di ; i++)
	    aa[150 * best + best + 5 *i] += aa[150 * best2 + best2] ;
	}
    }
  else
    {   /* standard blurred image, to present in the problem */
      di = 0 ;
      for (best2 =0 ; best2 < 4 ; best2++)
	if (best2 != best)
	  { 
	    for (i = di ; i < 30 ; i++)
	      for (j = 0 ; j < 5 ; j++)
		aa[150 * best + 5 * i + j] += aa[150 * best2 + 5 * (i-di) + j] ;
	    for (i = 0 ; i < di ; i++)
	      aa[150 * best + best + 5 *i] += aa[150 * best2 + best2] ;
	  }
    }
  if (1)
    {
      char *atgc = "ATGC" ;
      char ATGC[4] = {A_, T_, G_, C_} ;
      memset (adaptor, coded ? N_ : 'n' , 30) ; adaptor[30] = 0 ;
      
      for (i = 0 ; i < 30 ; i++)
	{
	  int i5 = best * 150 + i * 5 ;
	  for (j = 0, nn = 0 ; j < 5 ; j++)
	    nn += aa[i5 + j] ;
	  for (j = 0 ; j < 4 ; j++)
	    if (100 * aa[i5 + j] > 52 * nn && nn > 1000)
	      {
		if (j > 0) good++ ;
		nL += nn ;
		adaptor[i] = coded ? ATGC[j] : atgc[j] ;
	      }
	}
    }
  adaptor[31] = 0 ;
  if (bestp) *bestp = best ;
  if (nL < 20000)
    { memset (adaptor, 0, 30) ; good = 0 ; }
  return good >= 8 ? TRUE : FALSE ;
}

/**************************************************************/
/* read the adaptors for the frequency histogram of the overhangs */
BOOL saReadAdaptors (ADAPTORS *adaptors, RunSTAT *up, BOOL coded)
{
  BOOL ok = FALSE ;
  
  memset (adaptors, 0, sizeof(ADAPTORS)) ;
  ok |= readOneAdaptor (adaptors->a1L, up, TRUE, 0, 0, coded) ;
  ok |= readOneAdaptor (adaptors->a1R, up, FALSE, 0, 0, coded) ;
  if (up->p.nPairs)
    {
      ok |= readOneAdaptor (adaptors->a2L, up, TRUE, 1, 0, coded) ;
      ok |= readOneAdaptor (adaptors->a2R, up, FALSE, 1, 0, coded) ;
    }
  return ok ;
} /* saReadAdaptors */

/**************************************************************/

/* read1 then read2 */
static void saOverhangExport (const PP *pp, int run, RunSTAT *up, BOOL isLeft)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *runName = dictName (pp->runDict, run) ;
  char *title = hprintf (h, "%s.overhang.%s.tsf", runName, isLeft ? "5prime" : "3prime") ;
  long int nn = 0 ;
  long int *aa = 0 ;
  int i, j ;
  ACEOUT ao = aceOutCreate (pp->outFileName, title, 0, h) ;
  int best = 0 ;
  int pass, passMax = up->p.nPairs ? 2 : 1 ; 
  char *suffix ;
  char adaptor[32] = {0} ;
  
  for (pass = 0 ; pass < passMax ; pass++)
    {
      suffix = passMax == 1 ? "" : (pass == 0 ? "_R1" : "_R2") ;
      if (pass == 1) aceOutf (ao, "\n\n\n\n\n\n") ;
      aa = isLeft ? (pass ? up->overhangL2 : up->overhangL1) : (pass ? up->overhangR2 : up->overhangR1) ;
      aceOutDate (ao, "##", hprintf (h, "Run %s reverse complement of the %d prime unaligned overhang  profile limited to 32 letters", runName, isLeft ? 5 : 3)) ;

      best = 0 ;
      readOneAdaptor (adaptor, up, isLeft, pass, &best, FALSE) ;
      aceOutf (ao, "## Adaptor %s\n", adaptor) ;
      
      aceOutf (ao, "# Run.f\tPosition\tiiiiiifffff\tAny\tA\tT\tG\tC\tN\t%%A\t%%T\t%%G\t%%C\t%%N\n") ;
      for (i = 0 ; i < 30 ; i++)
	{
	  int i5 = best * 150 + i * 5 ;
	  for (j = 0, nn = 0 ; j < 5 ; j++)
	    nn += aa[i5 + j] ;
	  if (nn >= 100)
	    aceOutf (ao, "%s%s\t%d\tiiiiiifffff\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n"
		     , runName, suffix
		     , i + 1, nn, aa[i5 + 0], aa[i5 + 1], aa[i5 + 2], aa[i5 + 3], aa[i5 + 4]
		     , 100.0 * aa[i5 + 0]/nn, 100.0 * aa[i5 + 1]/nn, 100.0 * aa[i5 + 2]/nn, 100.0 * aa[i5 + 3]/nn, 100.0 * aa[i5 + 4]/nn
			 ) ;
	  else
	    aceOutf (ao, "%s%s\t%d\tiiiiii\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n"
		     , runName, suffix
		     , i + 1, nn, aa[i5 + 0], aa[i5 + 1], aa[i5 + 2], aa[i5 + 3], aa[i5 + 4]
		     ) ;
	  
	    }
    }
  
  aceOutf (ao, "\n\n") ;
  ac_free (h) ;
  return ;
} /* saOverhangExport */
  
/**************************************************************/
  
void saRunStatExport (const PP *pp, Array runStats)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (pp->outFileName, ".runStats.tsf", 0, h) ;
  int run, runMax = arrayMax (runStats) ;
  long int LD[7] ;

  memset (LD, 0, sizeof (LD)) ;
  
  aceOutDate (ao, "##", hprintf (h, "Run statistics on index %s, seed length %d, max repeats %d, error cost %d", pp->indexName, pp->seedLength, pp->maxTargetRepeats, pp->errCost)) ;
  for (run = 1 ; run < runMax ; run++)
    {
      RunSTAT *up = arrp (runStats, run, RunSTAT) ;
      saLetterProfileExport (pp, run, up) ;
      saOverhangExport (pp, run, up, TRUE) ;
      saOverhangExport (pp, run, up, FALSE) ;      
      if (up->p.nReads)
	{
	  const char *runNam = dictName (pp->runDict, run) ;
	  
	  aceOutf (ao, "%s\tPairs\ti\t%ld\n", runNam, up->p.nPairs) ;
	  aceOutf (ao, "%s\tAligned_pairs\tif\t%ld\t%.3f\n"
		   , runNam
		   , up->nPairsAligned
		   , 100.0 * up->nPairsAligned / (.000001 + up->p.nPairs)
		   ) ;
	  aceOutf (ao, "%s\tCompatible_pairs\tif\t%ld\t%.3f\n"
		   , runNam
		   , up->nCompatiblePairs
		   , 100.0 * up->nCompatiblePairs / (.000001 + up->nPairsAligned)
		   ) ;
	  aceOutf (ao, "%s\tCircle_pairs\tif\t%ld\t%.3f\n"
		   , runNam
		   , up->nCirclePairs
		   , 100.0 * up->nCirclePairs / (.000001 + up->nPairsAligned)
		   ) ;
	  aceOutf (ao, "%s\tNon_compatible_pairs\tif\t%ld\t%.3f\n"
		   , runNam
		   , up->nIncompatiblePairs
		   , 100.0 * up->nIncompatiblePairs / (.000001 + up->nPairsAligned)
		   ) ;
	  aceOutf (ao, "%s\tReads\ti\t%ld\n", runNam, up->p.nReads) ;
	  if (up->tooShort) aceOutf (ao, "%s\tToo_short_reads\ti\t%ld\n", runNam, up->tooShort) ;
	  if (up->lowEntropy)  aceOutf (ao, "%s\tLow_entropy_reads\ti\t%ld\n", runNam, up->lowEntropy) ;
	  if (up->tooShortBases) aceOutf (ao, "%s\tToo_short_bases\ti\t%ld\n", runNam, up->tooShortBases) ;
	  if (up->lowEntropyBases)  aceOutf (ao, "%s\tLow_entropy_bases\ti\t%ld\n", runNam, up->lowEntropyBases) ;
	  aceOutf (ao, "%s\tAligned_reads\tif\t%ld\t%.3f\n"
		   , runNam
		   , up->nMultiAligned[0]
		   , 100.0 * up->nMultiAligned[0] / (.000001 + up->p.nReads)
		   ) ;
	  aceOutf (ao, "%s\tPerfect_reads\tif\t%ld\t%.3f\n"
		   , runNam
		   , up->nPerfectReads
		   , 100.0 * up->nPerfectReads / (.000001 + up->p.nReads)
		   ) ;

	  aceOutf (ao, "%s\tBases\tiii\t%ld\t%ld\t%ld\n"
		   , runNam
		   , up->p.nBase1 + up->p.nBase2
		   , up->p.nBase1
		   , up->p.nBase2
		   ) ;
	  aceOutf (ao, "%s\tAligned_bases\tififif\t%ld\t%.3f\t%ld\t%.3f\t%ld\t%.3f\n"
		   , runNam
		   , up->nBaseAligned1 + up->nBaseAligned2
		   , 100.0*(up->nBaseAligned1 + up->nBaseAligned2)/(.0001 + up->p.nBase1 + up->p.nBase2)
		   , up->nBaseAligned1
		   , 100.0 * up->nBaseAligned1/(.000001 + up->p.nBase1)
		   , up->nBaseAligned2
		   , 100.0 * up->nBaseAligned2/(.000001 + up->p.nBase2)
		   ) ;

	  up->intergenic = up->wiggleCumul - up->cds - up->utr - up->intronic ;
	  aceOutf (ao, "%s\tGenome_length\ti\t%ld\n", runNam, pp->genomeLength) ;
	  aceOutf (ao, "%s\tCDS_utr_intronic_intergenic_Bases\tifififif\t%ld\t%.3f\t%ld\t%.3f\t%ld\t%.3f\t%ld\t%.3f\n"
		   , runNam
		   , up->cds
		   , 100.0*(up->cds)/(.000001 + up->cds + up->utr + up->intronic + up->intergenic)
		   , up->utr
		   , 100.0*(up->utr)/(.000001 + up->cds + up->utr + up->intronic + up->intergenic)
		   , up->intronic
		   , 100.0*(up->intronic)/(.000001 + up->cds + up->utr + up->intronic + up->intergenic)
		   , up->intergenic
		   , 100.0*(up->intergenic)/(.000001 + up->cds + up->utr + up->intronic + up->intergenic)
		   ) ;

	  aceOutf (ao, "%s\tATGCN\tiiiii\t%ld\t%ld\t%ld\t%ld\t%ld\n"
		   , runNam
		   , up->p.ATGCN[0]
		   , up->p.ATGCN[1]
		   , up->p.ATGCN[2]
		   , up->p.ATGCN[3]
		   , up->p.ATGCN[4]
		   ) ;


	  if (1)
	    {
	      ADAPTORS adaptors = {{0}} ;
	      saSetGetAdaptors (0, 0, &adaptors, run) ;
	      saReadAdaptors (&adaptors, up, FALSE) ;
	      if (1)
		{
		  if (up->nClippedAdaptor1L)
		    {
		      char *cp = adaptors.a1L ;
		      aceOutf (ao, "%s\tClipped_5prime_Adaptor_read1\tit\t%ld\t%s\n"
			       , runNam
			       , up->nClippedAdaptor1L
			       , cp 
			       ) ;
		    }
		  
		  if (up->nClippedAdaptor1R)
		    {
		      char *cp = adaptors.a1R  ;
		      aceOutf (ao, "%s\tClipped_3prime_Adaptor_read1\tit\t%ld\t%s\n"
			       , runNam
			       , up->nClippedAdaptor1R
			       , cp
			       ) ;
		    }
		  
		  if (up->nClippedAdaptor2L)
		    {
		      char *cp = adaptors.a2L ;
		      aceOutf (ao, "%s\tClipped_5prime_Adaptor_read1\tit\t%ld\t%s\n"
			       , runNam
			       , up->nClippedAdaptor2L
			       , cp
			       ) ;
		    }
		  
		  if (up->nClippedAdaptor2R)
		    {
		      char *cp = adaptors.a2R ;
		      aceOutf (ao, "%s\tClipped_3prime_Adaptor_read2\tit\t%ld\t%s\n"
			       , runNam
			       , up->nClippedAdaptor2R
			       , cp
			       ) ;
		    }
		}
	    }

	  for (int i = 0 ; i < SLMAX ; i++)
	    if (up->nClippedSls[i])
	      aceOutf (ao, "%s\tClipped_SL%d\tit\t%ld\t%s\n"
		       , runNam
		       , i, up->nClippedSls[i], "SL_sequence"
		       ) ;
	  if (up->nClippedPolyA)
	    aceOutf (ao, "%s\tClipped_polyA\ti\t%ld\n"
		     , runNam
		     , up->nClippedPolyA
		     ) ;
	  if (up->nClippedPolyT)
	    aceOutf (ao, "%s\tClipped_polyT\ti\t%ld\n"
		     , runNam
		     , up->nClippedPolyT
		     ) ;
	  
	  aceOutf (ao, "%s\tpolyA_sites\ti\t%ld\n"
		   , runNam
		   , confirmedPolyAsCountSites (pp, run)
		   ) ;
	  
	  aceOutf (ao, "%s\twiggleCumul\ti\t%ld\n", runNam, up->wiggleCumul) ;
	  aceOutf (ao, "%s\tErrors\tif\t%ld\t%.3f\n"
		   , runNam
		   , up->nErr
		   , 100.0 * up->nErr /(.0001 + up->nBaseAligned1 + up->nBaseAligned2)
		   ) ;
	  aceOutf (ao, "%s\tMismatches_and_InDels\tif\t%ld\t%.3f\n"
		   , runNam
		   , up->nMID
		   , 100.0 * up->nMID /(.0001 + up->nBaseAligned1 + up->nBaseAligned2)
		   ) ;
	  {{
	      int isRna = 0 ;
	      if (saSetGetAdaptors (0, &isRna, 0, run))
		aceOutf (ao, "%s\tStrategy\ti\t%d\n"
		   , runNam
		   , isRna
		   ) ;
	    }}		
	  aceOutf (ao, "%s\tIntron_supports\tiiif\t%ld\t%ld\t%ld\t%.2f%%\n"
		   , runNam
		   , up->gt_ag_Support + up->ct_ac_Support 
		   , up->gt_ag_Support
		   , up->ct_ac_Support
		   , up->intronStranding
		   ) ;
	  aceOutf (ao, "%s\tSupported_introns\ti\t%ld\n"
		   , runNam
		   , confirmedIntronsCountSites (pp, run)
		   ) ;
	  aceOutf (ao, "%s\tRaw_gt_ag_Support\ti\t%ld\n"
		   , runNam
		   , up->gt_ag_Support
		   ) ;
	  aceOutf (ao, "%s\tRaw_ct_ac_Support\ti\t%ld\n"
		   , runNam
		   , up->ct_ac_Support
		   ) ;
	  aceOutf (ao, "%s\tMin_read_length\ti\t%d\n", runNam, up->p.minReadLength) ;
	  aceOutf (ao, "%s\tMax_read_length\ti\t%d\n", runNam, up->p.maxReadLength) ;
	  if (up->p.lengthDistribution)
	    {
	      long int j = 0 ;
	      int mode = 0 ;
	      int i, iMax = arrayMax (up->p.lengthDistribution) ;
	      long int *xp = arrp (up->p.lengthDistribution, 0, long int) ;
	      for (i = 0 ; i < iMax ; i++, xp++)
		{
		  j += *xp ;
		  LD[6] += i * (*xp) ;
		  if (*xp > LD[5])
		    { LD[5] = *xp ; mode = i ; }
		}
	      LD[5] = mode ;
	      LD[6] = j ? LD[6]/j : 0 ; /* average */
	      /* construct the cumulated distrib */
	      xp = arrp (up->p.lengthDistribution, 0, long int) ;
	      for (i = 1 ; i < iMax ; i++, xp++)
		{
		  xp[0] += xp[-1] ;
		  if (! LD[0] && *xp >= j/100) LD[0] = i ;
		  if (! LD[1] && *xp >= 5*j/100) LD[1] = i ;
		  if (! LD[2] && *xp >= 50*j/100) LD[2] = i ;
		  if (! LD[3] && *xp >= 95*j/100) LD[3] = i ;
		  if (! LD[4] && *xp >= 99*j/100) LD[4] = i ;
		}
	      aceOutf (ao, "%s\tLength_distribution_1_5_50_95_99_mode_av\tiiiiiii", runNam) ;
	      for (i=0;i<7;i++)
		aceOutf (ao, "\t%ld", LD[i]) ;
	      aceOutf (ao, "\n") ;
	    }

	  if (up->insertLengthDistribution && arrayMax (up->insertLengthDistribution))
	    {
	      long int j = 0 ;
	      int mode = 0 ;
	      int i, iMax = arrayMax (up->insertLengthDistribution) ;
	      long int *xp = arrp (up->insertLengthDistribution, 0, long int) ;

	      memset (LD, 0, sizeof(LD)) ;
	      for (i = 0 ; i < iMax ; i++, xp++)
		{
		  j += *xp ;
		  LD[6] += i * (*xp) ;
		  if (*xp > LD[5])
		    { LD[5] = *xp ; mode = i ; }
		}
	      LD[5] = mode ;
	      LD[6] = j ? LD[6]/j : 0 ; /* average */
	      /* construct the cumulated distrib */
	      xp = arrp (up->p.lengthDistribution, 0, long int) ;
	      for (i = 1 ; i < iMax ; i++, xp++)
		{
		  xp[0] += xp[-1] ;
		  if (! LD[0] && *xp >= j/100) LD[0] = i ;
		  if (! LD[1] && *xp >= 5*j/100) LD[1] = i ;
		  if (! LD[2] && *xp >= 50*j/100) LD[2] = i ;
		  if (! LD[3] && *xp >= 95*j/100) LD[3] = i ;
		  if (! LD[4] && *xp >= 99*j/100) LD[4] = i ;
		}
	      aceOutf (ao, "%s\tInsertLength_distribution_1_5_50_95_99_mode_av\tiiiiiii", runNam) ;
	      for (i=0;i<7;i++)
		aceOutf (ao, "\t%ld", LD[i]) ;
	      aceOutf (ao, "\n") ;
	    }
	  
	  for (int ii = 1 ; ii < 2 ; ii++)
	    aceOutf (ao, "%s\tReads_Aligned_once\tif\t%ld\t%.3f\n"
		     , runNam
		     , up->nMultiAligned[ii]
		     , 100.0 * up->nMultiAligned[ii]/(.000001 + up->nMultiAligned[0]) 
		     ) ;
	  for (int ii = 2 ; ii < 11 ; ii++)
	    if (up->nMultiAligned[ii]) 
	      aceOutf (ao, "%s\tReads_multi_aligned__%d\tif\t%ld\t%.3f\n"
		       , runNam, ii
		       , up->nMultiAligned[ii]
		       , 100.0 * up->nMultiAligned[ii]/(.000001 + up->nMultiAligned[0])
		       ) ;
	  
	  for (int ii = 1 ; ii < 256 ; ii++)
	    {
	      if (up->nReadsAlignedPerTargetClass[ii])
		aceOutf (ao, "%s\tReads_aligned_in_class_%c\tif\t%ld\t%.3f\n"
			 , runNam, ii
			 , up->nReadsAlignedPerTargetClass[ii]
			 , 100.0 * up->nReadsAlignedPerTargetClass[ii]/(.000001 + up->nMultiAligned[0]) 
			 ) ;
	    }
	  for (int ii = 1 ; ii < 256 ; ii++)
	    {
	      if (up->nBasesAlignedPerTargetClass[ii])
		aceOutf (ao, "%s\tBases_aligned_in_class_%c\tif\t%ld\t%.3f\n"
			 , runNam, ii
			 , up->nBasesAlignedPerTargetClass[ii]
			 , 100.0 * up->nBasesAlignedPerTargetClass[ii]/(.000001 + up->nBasesAlignedPerTargetClass[0]) 
			 ) ;
	    }
	  for (int ii = 1 ; ii < 256 ; ii++)
	    {
	      if (up->nReadsAlignedPerTargetClass[ii])
		{
		  int f = up->GF[ii] ;
		  int r = up->GR[ii] ;
		  int t = f + r ;
		  if (t)		  
		    aceOutf (ao, "%s\tStranding_in_class_%c\tfii\t%.3f\t%ld\t%ld\n", runNam, ii
			     , 100.0 * f/t, f, r
			     ) ;
		}
	    }
	  if (up->gt_ag_Support + up->ct_ac_Support  > 10000)
	    {
	      aceOutf (ao, "%s\tStranding_in_class_I\tfii\t", runNam) ;
	      aceOutPercent (ao, (100.0 * up->gt_ag_Support) / (up->gt_ag_Support + up->ct_ac_Support)) ;
	      aceOutf (ao, "\t%ld\t%ld\n", up->gt_ag_Support, up->ct_ac_Support ) ;
	    }	  
	}
    }
  ac_free (h) ;

  s2gSamStatsExports (pp, runStats) ;

  return ;    
} /* saRunStatExport */

/**************************************************************/
