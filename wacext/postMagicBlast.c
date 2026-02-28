/*
 * authors: Danielle and Jean Thierry-Mieg, NCBI, 
 *  Copyright (C) D and J Thierry-Mieg 2020
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
 * This file is part of the AceView/Magic project developped by
 *	 D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *
 *  This should be linked with the acedb libraries
 *  available from http://www.acedb.org
 *
 *  Please direct all question about this code to
 *  mieg@ncbi.nlm.nih.gov

 * 15feb 2020
 * postmagicBlast.c
 *   Input: Output from MagicBlast (sam or tabular)
 *   Output: 
 *      best ali, ali.ace statistics, introns, SNPs 
 *      i.e. every interesting post-treatment we can think of
 *      using the additive .tsf format
 */
#define MALLOC_CHECK  
#define ARRAY_CHECK  

#include "../wac/ac.h"

typedef struct pmbStruct { 
  AC_HANDLE h ;
  
  const char *run, *targetClassName ; 
  
  const char *inFileList ; 
  const char *outFileName ; 

  const char *geneInfoFileName ;
  const char *splitMrnaFileName ;
  BOOL gzi, gzo ;
  BOOL sam, tabular, clipali ; /* inpt in SAM format, default:tabular */
  BOOL getIntrons ;
  BOOL getSnp ;
  BOOL getExpression ;
  BOOL hasPair ;
  BOOL showStats ;

  DICT *cloneDict, *targetDict, *target_classDict, *cigarDict, *intronDict, *snpDict, *runDict ;
  DICT *dnaDict, *geneDict ;
  BigArray alis, geneCounts, mrnaCounts, snps ;
  Array introns, geneInfo, mrnaInfo ;

  int bonus ;
  int strand ;
  int target_class ;
  int maxIntronLength ;
  int incompatibleTopology ;

  int alignedFragments ;
  int compatiblePairs ;
  int compatiblePairsInsideGene ;
  int compatiblePairsInGenome ;

  int r1Up, r1Down, r2Up, r2Down ;
} PMB ;

typedef struct aliStruct {
  int clone, mult, flag
    , target, target_class, gene
    , a1, a2, x1, x2, uu
    , score, pairScore
    , ln, ali, pairAli
    , nN, nErr
    , nintron, cigar, targetCigar
    , prefix, suffix
    , dna ; 
  BOOL isRead1 ;
} ALI ;
typedef struct intronStruct { int target, a1, a2, nr1, nr2 ; char feet[8] ;  } ITR ;
typedef struct countStruct {
  int tc, gene, mrna, seqs, tags, ali ;
} CNT ;
typedef struct geneInfoStruct {
  int target_class, mrna, gene, ln, gc, geneId ;
} GINFO ;
typedef struct snpStruct {
  int name, run, a1, a2, da, c, m, w, mp, wp, mm, wm, method, type, dx, suffix ;
} SNP ;


static int pmbTabularSmoke (PMB *pmb, ALI *up, vTXT txt) ;

/*************************************************************************************/

static int bestAliOrder (const void *va, const void *vb)
{
  const ALI *up = (const ALI *) va ;
  const ALI *vp = (const ALI *) vb ;
  int n ;

  n = up->clone - vp->clone ; if (n) return n ;
  n = (int)up->isRead1 - ((int)vp->isRead1) ; if (n) return -n ;
  n = up->pairScore - vp->pairScore ; if (n) return -n ;
  n = up->score - vp->score ; if (n) return -n ;
  n = up->target_class - vp->target_class ; if (n) return n ;
  n = up->target - vp->target ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;

  return 0 ;
} /* bestAliOrder */

/*************************************************************************************/
/* here all scores are alreay accepted, sort by gene/mrna/coords */
static int countOrder (const void *va, const void *vb)
{
  const ALI *up = (const ALI *) va ;
  const ALI *vp = (const ALI *) vb ;
  int n ;

  n = up->clone - vp->clone ; if (n) return n ;
  n = up->target_class - vp->target_class ; if (n) return n ;
  n = up->gene - vp->gene ; if (n) return n ;
  n = up->target - vp->target ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;

  n = (int)up->isRead1 - ((int)vp->isRead1) ; if (n) return -n ;
  n = up->pairScore - vp->pairScore ; if (n) return -n ;
  n = up->score - vp->score ; if (n) return -n ;


  return 0 ;
} /* countOrder */

/*************************************************************************************/
#ifdef JUNK

static int intronOrder (const void *va, const void *vb)
{
  const ITR *up = (const ITR *) va ;
  const ITR *vp = (const ITR *) vb ;
  int n ;

  n = up->target - vp->target ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;

  return 0 ;
} /* intronOrder */

#endif
/*************************************************************************************/

static void pmbPairScore (PMB *pmb)
{
  int i, j, iMax = bigArrayMax (pmb->alis) ;
  ALI *up, *vp ;
  int maxIntron = pmb->maxIntronLength ;

  bigArraySort (pmb->alis, bestAliOrder) ;

  for (i = 0, up = bigArrp (pmb->alis, i, ALI) ; i < iMax ; i++, up++) 
    {
      int clone = up->clone ;
      int isRead1 = up->isRead1 ;

      if (! up->pairScore)
	for (j = i + 1, vp = up + 1 ; j < iMax ; j++, vp++)
	  {
	    if (vp->clone != clone)
	      break ;
	    if (! vp->pairScore 
		&& vp->isRead1 != isRead1 
		&& vp->target == up->target
		&& vp->target_class == up->target_class
		&& ( 
		  (up->a1 < up->a2 && vp->a1 > vp->a2 && vp->a1 > up->a1 && vp->a1 - up->a1 < maxIntron)
		  || (up->a1 > up->a2 && vp->a1 < vp->a2 && vp->a1 < up->a1 && up->a1 - vp->a1 < maxIntron)
		     )
		)
	      {
		up->pairScore = vp->pairScore = up->score + vp->score ;
		up->pairAli = up->ali + vp->ali ;
	      }
	  }
    }

  return ;
} /* pmbPairScore */

/*************************************************************************************/

static void pmbSelect (PMB *pmb)
{
  long int ii, jj, kk, iiMax = bigArrayMax (pmb->alis) ;
  ALI *up, *vp, *wp ;
  int clone  ;
  int genome ;
  DICT *target_classDict = pmb->target_classDict ;
  int MX = 1000000 ;
  Array mrnaInfo = pmb->mrnaInfo ;
  int mrnaMax = mrnaInfo ? arrayMax (mrnaInfo) : 0 ;


  fprintf (stderr ,"pmbSelect before sort %ld %ld\n", iiMax, iiMax * pmb->alis->size) ;
  bigArraySort (pmb->alis, bestAliOrder) ;
  fprintf (stderr ,"pmbSelect after sort %ld\n", iiMax) ;
  dictAdd (target_classDict, "Z_genome", &genome) ;

  fprintf (stderr ,"pmbSelect old iMax = %ld\n", iiMax) ;
  if (iiMax)
    {
      for (ii = jj = 0, up = vp = bigArrp (pmb->alis, ii, ALI) ; ii < iiMax ; ii++, up++) 
	{
	  if (!up->score || (up->clone == vp->clone && up->isRead1 == vp->isRead1 && up->pairScore < vp->pairScore))
	    continue ;
	  if (up > vp) *vp = *up ;
	  vp++ ; jj++ ;
	}
      iiMax = bigArrayMax (pmb->alis) = jj ;
    }


  fprintf (stderr ,"pmbSelect new iMax = %ld\n", iiMax) ;
  if (iiMax)
    for (clone = 0, ii =  0, up = bigArrp (pmb->alis, ii, ALI) ; ii < iiMax ; ii++, up++) 
      {
	BOOL cpa = FALSE ;
	BOOL cpaGenome = FALSE ;
	BOOL cpaGene = FALSE ;
	int mrna = up->target ;
	
	if (! pmb->clipali && mrnaInfo && mrna && mrna < mrnaMax)
	  {
	    GINFO *mp = arrp (mrnaInfo, mrna, GINFO) ;	
	    up->gene = mp->gene ;
	  }
	if (up->clone == clone)
	  continue ;
	clone = up->clone ;
	
	for (vp = up, jj = ii ; jj < iiMax && vp->clone == clone && vp->isRead1 ; jj++, vp++)
	  for (wp = vp + 1, kk = jj + 1 ; kk < iiMax &&  wp->clone == clone ; kk++, wp++)
	    if (! wp->isRead1)
	      {
		if (vp->target == wp->target 
		    && 
		    (
		     (vp->a1 < vp->a2 && vp->a1 < wp->a1 && vp->a1 + MX > wp->a1 && wp->a1 > wp->a2) ||
		     (vp->a1 > vp->a2 && vp->a1 > wp->a1 && vp->a1 - MX < wp->a1 && wp->a1 < wp->a2) 
		     )
		    )
		  {
		    cpa = TRUE ;
		    if (wp->target_class == genome) cpaGenome = TRUE ;
		    if (wp->target_class && wp->target_class <= dictMax (target_classDict) && dictName (target_classDict, wp->target_class)[1] == 'T') cpaGene = TRUE ;
		  }
	      }
	
	pmb->alignedFragments += up->mult ;
	if (cpa) pmb->compatiblePairs += up->mult ;
	if (cpaGene) pmb->compatiblePairsInsideGene += up->mult ;
	if (cpaGenome) pmb->compatiblePairsInGenome += up->mult ;
      }
  return ;
} /* pmbSelect */

/*************************************************************************************/

static void pmbCount (PMB *pmb)
{
  BigArray alis = pmb->alis ;
  BigArray geneCounts = pmb->geneCounts ;
  BigArray mrnaCounts = pmb->mrnaCounts ;
  /*   Array geneInfo = pmb->geneInfo ; */
  Array mrnaInfo = pmb->mrnaInfo ;
  long int ii, jj, iMax = bigArrayMax (alis) ;
  ALI *ap0, *ap ;
  GINFO *mp ;
  vTXT namTxt = vtxtHandleCreate (pmb->h) ;
  int nr1 = 0, nr2 = 0 ;

  if ((pmb->getExpression || pmb->introns || pmb->getSnp) && iMax)
    {
        bigArraySort (alis, countOrder) ;
	fprintf (stderr, "// parseCount sort0 done: %s\n", timeShowNow()) ;
	
	for (ii = 0, ap0 = bigArrayp (alis, 0, ALI) ; ii < iMax ; ii++, ap0++)
	  {
	    int clone = ap0->clone, geneOld = 0, mrnaOld = 0 ;
	    for (jj = ii, ap = bigArrayp (alis, jj, ALI) ; jj < iMax && ap->clone == clone ; jj++, ap++)
	      {
		int mrna = ap->target ;
		
		if (ap->isRead1) 
		  nr1 += ap->mult ;
		else
		  nr2 += ap->mult ;


		if (! pmb->clipali && ap->cigar)
		  {
		    if (pmb->tabular)
		      ap->nErr = pmbTabularSmoke (pmb, ap, namTxt) ; /* isRead1 */
		  }
		if (ap->score == ap->pairScore || ap->isRead1)
		  if (pmb->getExpression && ap->score > 0 && mrna && mrna != mrnaOld && (pmb->clipali || mrna <arrayMax (mrnaInfo)))
		    {
		      int gene ;
		      CNT *mcp = bigArrayp (mrnaCounts, mrna, CNT) ;
		      
		      mrnaOld = mrna ;
		      mcp->seqs++ ;
		      mcp->tags += ap->mult ; 
		      mcp->ali += ap->ali * ap->mult ; 
		      if (! pmb->clipali)
			{
			  mp = arrayp (mrnaInfo, mrna, GINFO) ;	
			  gene = ap->gene = mp->gene ;
			}
		      else
			gene = ap->gene ;
		      mcp->gene = gene ;
		      mcp->mrna = mrna ;
		      if (gene && gene != geneOld)
			{
			  CNT *gcp = bigArrayp (geneCounts, gene, CNT) ;
			  gcp->gene = geneOld = gene ;
			  gcp->seqs++ ;
			  gcp->tags += ap->mult ; 
			  gcp->ali += ap->ali * ap->mult ; 
			}
		    }
	      }
	  }
	fprintf (stderr, "// parseCount sort1 done: %s\n", timeShowNow()) ;
	bigArraySort (alis, bestAliOrder) ;
	fprintf (stderr, "// parseCount sort2 done: %s\n", timeShowNow()) ;
	fprintf (stderr, "// nr1 = %d nr2 = %d \n", nr1, nr2) ;
    }
} /* pmbCount */

/*************************************************************************************/

static char *flipFeet (const char *old, char *new)
{
  new[0] = ace_lower(complementLetter(old[4])) ;
  new[1] = ace_lower(complementLetter(old[3])) ;
  new[2] = '_' ;
  new[3] = ace_lower(complementLetter(old[1])) ;
  new[4] = ace_lower(complementLetter(old[0])) ;
  new[5] = 0 ;
  return new ;
}

/*************************************************************************************/

static void pmbIntronsExport (PMB *pmb)
{
  int i, iMax = dictMax (pmb->intronDict) ;
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (pmb->outFileName, ".introns.tsf", pmb->gzo, h) ;
  ITR *ip ;
  DICT *targetDict = pmb->targetDict ;
  int isDown ;
  char newFeet[6] ;

  int nF = pmb->r1Down + pmb->r2Up ;
  int nR = pmb->r2Down + pmb->r1Up ;
  
  if (100 * nF > 55 * (nF + nR + 100)) isDown = 1 ;
  if (100 * nR > 55 * (nF + nR + 100)) isDown = -1 ;

  aceOutDate (ao, "###", "Introns reported by Magic-Blast") ;
  if (iMax)
    for (i = 1, ip = arrp (pmb->introns, 1, ITR) ; i <= iMax ; i++, ip++)
      {
	int a1 = ip->a1 ;
	int a2 = ip->a2 ;
	int n = ip->nr1 + ip->nr2 ;
	if (n > 0)
	  {
	    if (isDown == 1) /* flip the read2 */
	      {
		aceOutf (ao, "%s__%d_%d\t%s\tiiit\t%d\t%d\t%d\t%s\n", dictName (targetDict, ip->target), a1, a2, pmb->run, ip->nr1, ip->nr1, 0, ip->feet) ;
		aceOutf (ao, "%s__%d_%d\t%s\tiiit\t%d\t%d\t%d\t%s\n", dictName (targetDict, ip->target), a2, a1, pmb->run, ip->nr2, 0, ip->nr2, flipFeet(ip->feet, newFeet)) ;
	      }
	    else if (isDown == -1) /* flip the read1 */
	      {
		aceOutf (ao, "%s__%d_%d\t%s\tiiit\t%d\t%d\t%d\t%s\n", dictName (targetDict, ip->target), a1, a2, pmb->run, ip->nr2,  0, ip->nr2, ip->feet) ;
		aceOutf (ao, "%s__%d_%d\t%s\tiiit\t%d\t%d\t%d\t%s\n", dictName (targetDict, ip->target), a2, a1, pmb->run, ip->nr1,  ip->nr1, 0, flipFeet(ip->feet, newFeet)) ;
	      }
	    else
	      aceOutf (ao, "%s__%d_%d\t%s\tiiit\t%d\t%d\t%d\t%s\n", dictName (targetDict, ip->target), a1, a2, pmb->run,  n, ip->nr1, ip->nr2, ip->feet) ;
	  }
      }
  ac_free (h) ;
  fprintf (stderr, "Exported %d introns\n", iMax) ;
} /* pmbIntronsExport */

/*************************************************************************************/

static void pmbSnpsExpport (PMB *pmb)
{
  long int i, iMax = bigArrayMax (pmb->snps) ;
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (pmb->outFileName, ".snps.tsf", pmb->gzo, h) ;
  SNP *snp ;
  DICT *cloneDict = pmb->cloneDict ;
  DICT *snpDict = pmb->snpDict ;
  DICT *runDict = pmb->runDict ;

  aceOutDate (ao, "###", "Snps reported by Magic-Blast") ;
  aceOut (ao, " snp\rrun\tformat\ta1\ta2\tda\tc\tm\tw\tmp\twp\tmm\twm\tmethod\ttype\tdx\tsuffix\n") ;

  if (iMax)
    for (i = 1, snp = bigArrp (pmb->snps, 1, SNP) ; i <= iMax ; i++, snp++)
      {
	if (snp->a1)
	  {
	    aceOutf (ao, "%s\t%s\t", dictName (snpDict, snp->name),  dictName (runDict, snp->run)) ;
	    aceOut (ao, "tttiiiiiiiitttt") ;
	    aceOutf (ao, "\t%d\t%d\t%d", snp->a1, snp->a2, snp->da) ;   /* non additve, t format */
	    aceOutf (ao, "\t%d", snp->c) ; /* 7 additive counts */
	    aceOutf (ao, "\t%d\t%d", snp->m, snp->w) ; 
	    aceOutf (ao, "\t%d\t%d", snp->mp, snp->wp) ;
	    aceOutf (ao, "\t%d\t%d", snp->mm, snp->wm) ;
	    aceOutf (ao, "\t%s\t%s", dictName (cloneDict, snp->method), dictName (cloneDict, snp->type)) ;
	    aceOutf (ao, "\t%d\t%s\n", snp->dx, snp->suffix ? dictName (cloneDict, snp->suffix) : "-") ;
	  }
      }
  ac_free (h) ;
  fprintf (stderr, "Exported %ld snps\n", iMax) ;
} /* pmbSnpsExport */

/*************************************************************************************/

static void pmbGeneExpressionExport (PMB *pmb)
{
  BigArray aa = pmb->geneCounts ;
  int i, iMax = bigArrayMax (aa), nn = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (pmb->outFileName, ".genes.tsf", pmb->gzo, h) ;
  CNT *up ;
  DICT *geneDict = pmb->geneDict ;

  aceOutDate (ao, "###", "Gene expression reported by Magic-Blast") ;
  if (iMax)
    for (i = 1, up = bigArrp (aa, 1, CNT) ; i < iMax ; i++, up++)
      {
	if (up->tags && up->gene)
	  {
	    nn++ ;
	    aceOutf (ao, "%s\t%s\tiii\t%d\t%d\t%d\n", dictName (geneDict, up->gene), pmb->run,  up->seqs, up->tags, up->ali) ;
	  }
      }
  ac_free (h) ;
  fprintf (stderr, "Exported %d genes\n", nn) ;
} /* pmbGeneExpressionExport */

/*************************************************************************************/

static void pmbMrnaExpressionExport (PMB *pmb)
{
  BigArray aa = pmb->mrnaCounts ;
  int i, iMax = bigArrayMax (aa) ;
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (pmb->outFileName, ".mrnas.tsf", pmb->gzo, h) ;
  CNT *up ;
  DICT *mrnaDict = pmb->targetDict ;

  aceOutDate (ao, "###", "mRNA expression reported by Magic-Blast") ;
  if (iMax)
    for (i = 1, up = bigArrp (aa, 1, CNT) ; i < iMax ; i++, up++)
      {
	if (up->tags && up->mrna)
	  aceOutf (ao, "%s\t%s\tiii\t%d\t%d\t%d\n", dictName (mrnaDict, up->mrna), pmb->run,  up->seqs, up->tags, up->ali) ;
      }
  ac_free (h) ;
} /* pmbMrnaExpressionExport */

/*************************************************************************************/

static void pmbClassCount (PMB *pmb)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, ii, j, tc ;
  long int iMax = bigArrayMax (pmb->alis) ;
  ALI *up, *vp ;
  int clone = 0 , isRead1 = 0 ;
  int cMax = 1 + dictMax (pmb->target_classDict) ;
  int ns[cMax], nr[cMax], nali[cMax] ;
  int hns[cMax], hnr[cMax], hnali[cMax] ;
  int nrp[cMax], nrm[cMax], nra[cMax] ;
  int uuu[cMax][12] ;
  int perfect1 = 0 ;
  int perfect2 = 0 ;
  int nintronSupport = 0 ;
  Array firstBase1 = arrayHandleCreate (256, int, h) ;
  Array firstBase2 = arrayHandleCreate (256, int, h) ;
  Array lastBase1 = arrayHandleCreate (256, int, h) ;
  Array lastBase2 = arrayHandleCreate (256, int, h) ;
  Array nerrs1 = arrayHandleCreate (256, int, h) ;
  Array nerrs2 = arrayHandleCreate (256, int, h) ;
  Array aliLn = arrayHandleCreate (256, int, h) ;

  memset (ns, 0, sizeof(ns)) ;
  memset (nr, 0, sizeof(nr)) ;
  memset (nali, 0, sizeof(nali)) ;
  
  memset (hns, 0, sizeof(hns)) ;
  memset (hnr, 0, sizeof(hnr)) ;
  memset (hnali, 0, sizeof(hnali)) ;
  
  memset (nrp, 0, sizeof(nrp)) ;
  memset (nrm, 0, sizeof(nrm)) ;
  memset (nra, 0, sizeof(nra)) ;

  memset (uuu, 0, sizeof(uuu)) ;

  if (iMax)
    for (ii = 0, up = bigArrp (pmb->alis, ii, ALI) ; ii < iMax ; ii++, up++) 
      {
	BOOL ok = FALSE ;
	if (up->clone == clone && up->isRead1 == isRead1)
	  continue ;
	clone = up->clone ;
	isRead1 = up->isRead1 ;
	if (isRead1)
	  array (firstBase1, up->x1, int) += up->mult ;
	else
	  array (firstBase2, up->x1, int) += up->mult ;
	if (isRead1)
	  array (lastBase1, up->x2, int) += up->mult ;
	else
	array (lastBase2, up->x2, int) += up->mult ;
	{
	  int k = up->ali ;
	  int k1 = k/20 ;
	  if (k1>0) {k = k/k1 ; k = k * k1 ; }
	  array (aliLn, k, int) += up->mult ;
	}

	if (up->x1 == 1 && up->x2 == up->ln && up->nErr == 0)
	  {
	    if (isRead1)
	      perfect1 += up->mult ;
	    else
	      perfect2 += up->mult ;
	  }
	nintronSupport += up->mult * up->nintron ;
	
	if (isRead1)
	  array (nerrs1, up->nErr, int) += up->mult ;
	else
	  array (nerrs2, up->nErr, int) += up->mult ;
	
	for (tc = 0 ; tc < cMax ; tc++)
	  {
	    ok = FALSE ;
	    for (j = ii, vp = up ; j < iMax && vp->clone == clone && vp->isRead1 == isRead1 ; j++, vp++)
	      {
		if (vp->target_class == tc)
		  {
		    int k, s = 0 ;
		    ALI *wp ;
		    int u = vp->uu ;
		    
		    if (u > 10) u = 10 ;
		    ns[tc]++ ; nr[tc] += vp->mult ; nali[tc] += vp->ali ;
		    uuu[tc][u] += vp->mult ;
		    for (k = j, wp = vp ; k < iMax && wp->clone == clone && wp->isRead1 == isRead1 ; k++, wp++)
		      if (wp->target_class == tc)
			{
			  if (wp->a1 < wp->a2) s |= 0x1 ;
			  else s |= 0x2 ;
			}
		    if (s > 0 && s < 3 && pmb->strand == -1)
		      s = 3 - s ;   /* flip values 1 and 2 */
		    if (vp->isRead1)
		      switch (s)
			{
			case 1: nrp[tc] += vp->mult ; break ;
			case 2: nrm[tc] += vp->mult ; break ;
			case 3: nra[tc] += vp->mult ; break ;
			}
		    else
		      switch (s)
			{
			case 2: nrp[tc] += vp->mult ; break ;
			case 1: nrm[tc] += vp->mult ; break ;
			case 3: nra[tc] += vp->mult ; break ;
			}
		    
		    if (!ok)
		      {
			hns[tc]++ ; hnr[tc] += vp->mult ; hnali[tc] += vp->ali ;
		      }
		    break ;
		  }
		ok = TRUE ;
	      }
	  }
      }
  
  if (cMax && pmb->showStats)
    {
      ACEOUT ao = aceOutCreate (pmb->outFileName, ".stats.tsf", pmb->gzo, h) ;
      
      if (! ao)
	messcrash ("Cannot open outputfile %s.stats.tsf", pmb->outFileName ? pmb->outFileName : "./") ;
      aceOutDate (ao, "###", "Lane stats exported by postMagicBlast.c tabular format") ;

      for (tc = 0 ; tc < cMax ; tc++)
	{
	  if (hns[tc] > 0)
	aceOutf (ao, "h_ALI__%s\t%s\tiii\t%d\t%d\t%d\n", tc ? dictName (pmb->target_classDict, tc) : "zero", pmb->run, hns[tc], hnr[tc], hnali[tc]/1000) ;
	  if (ns[tc] > 0)
	    aceOutf (ao, "nh_ALI__%s\t%s\tiii\t%d\t%d\t%d\n", tc ? dictName (pmb->target_classDict, tc) : "zero", pmb->run, ns[tc], nr[tc], nali[tc]/1000) ;
	  if (nr[tc] > 10000)
	    aceOutf (ao, "Stranding__%s\t%s\tiiii\t%d\t%d\t%d\t%d\t//%.2f\t%.2f\t%.2f\n", tc ? dictName (pmb->target_classDict, tc) : "zero", pmb->run
		     , nr[tc], nrp[tc], nrm[tc], nra[tc], 100.0*nrp[tc]/nr[tc], 100.0*nrm[tc]/nr[tc], 100.0*nra[tc]/nr[tc]) ;
	  if (nr[tc] > 0)
	    {
	      int j ;
	      aceOutf (ao, "Unicity__%s\t%s\t10i", tc ? dictName (pmb->target_classDict, tc) : "zero", pmb->run) ;
	      for (j = 1 ; j <= 10 ; j++)
		aceOutf (ao, "\t%d", uuu[tc][j]) ;
	      aceOutf (ao, "\n") ;
	    }
	}
      
      
      array (firstBase2, arrayMax (firstBase1), int) += 0 ;
      array (firstBase1, arrayMax (firstBase2), int) += 0 ;
      for (i = 1 ; i < arrayMax (firstBase1) ; i++)
	{
	  int n1 = array (firstBase1, i, int) ;
	  int n2 = array (firstBase2, i, int) ;
	  if (n1 + n2)
	    {
	      if (pmb->hasPair)
		aceOutf (ao, "First_base_aligned__%d\t%s\tii\t%d\t%d\n", i, pmb->run, n1, n2) ;
	      else
		aceOutf (ao, "Last_base_aligned__%d\t%s\ti\t%d\n", i, pmb->run, n1) ;
	    }
	}
      
      array (lastBase2, arrayMax (lastBase1), int) += 0 ;
      array (lastBase1, arrayMax (lastBase2), int) += 0 ;
      for (i = 1 ; i < arrayMax (lastBase1) ; i++)
	{
	  int n1 = array (lastBase1, i, int) ;
	  int n2 = array (lastBase2, i, int) ;
	  if (n1 + n2)
	    {
	      if (pmb->hasPair)
		aceOutf (ao, "Last_base_aligned__%d\t%s\tii\t%d\t%d\n", i, pmb->run, n1, n2) ;
	      else
		aceOutf (ao, "Last_base_aligned__%d\t%s\ti\t%d\n", i, pmb->run, n1) ;
	    }
	}
      
      if (pmb->hasPair)  /* export the same number twwice: float then int in the ace database */
	aceOutf (ao, "Perfect_reads\t%s\tii\t%d\t%d\n", pmb->run, perfect1 + perfect2, perfect1 + perfect2) ;
      else
	aceOutf (ao, "Perfect_reads\t%s\tii\t%d\n", pmb->run, perfect1, perfect1) ;
      
      if (pmb->hasPair)
	{
	  aceOutf (ao, "Aligned_fragments\t%s\ti\t%d\n", pmb->run, pmb->alignedFragments) ;
	  aceOutf (ao, "Compatible_pairs\t%s\ti\t%d\n", pmb->run, pmb->compatiblePairs) ;
	  aceOutf (ao, "Compatible_pairs_inside_gene\t%s\ti\t%d\n", pmb->run, pmb->compatiblePairsInsideGene) ;
	  aceOutf (ao, "Compatible_pairs_in_genome\t%s\ti\t%d\n", pmb->run, pmb->compatiblePairsInGenome) ;
	  aceOutf (ao, "Incompatible_topology\t%s\ti\t%d\n", pmb->run, pmb->incompatibleTopology) ;
	}
      
      array (nerrs2, arrayMax (nerrs1), int) += 0 ;
      array (nerrs1, arrayMax (nerrs2), int) += 0 ;
      for (i = 0 ; i < arrayMax (nerrs1) ; i++)
	{
	  int n1 = array (nerrs1, i, int) ;
	  int n2 = array (nerrs2, i, int) ;
	  if (n1 + n2)
	    {
	      if (pmb->hasPair)
		aceOutf (ao, "Count_mismatch__%d\t%s\ti\t%d\n", i, pmb->run, n1 + n2) ;
	      else
		aceOutf (ao, "Count_mismatch__%d\t%s\ti\t%d\n", i, pmb->run, n1) ;
	    }
	}
      
      for (i = 0 ; i < arrayMax (aliLn) ; i++)
	{
	  int n1 = array (aliLn, i, int) ;
	  if (n1)
	    aceOutf (ao, "Aligned_length__%d\t%s\ti\t%d\n", i, pmb->run, n1) ;
	}
      
      if (0)  /* bad format */
	aceOutf (ao, "N_intron_support\t%s\ti\t%d\n",  pmb->run, nintronSupport) ;
    }

  ac_free (h) ;
  return ;
} /* pmbClassCount */

/*************************************************************************************/

static void pmbParseGeneInfo (PMB *pmb)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (pmb->geneInfoFileName, FALSE, h) ;
  DICT *mrnaDict = pmb->targetDict ;
  DICT *geneDict = pmb->geneDict ;
  Array aaG = pmb->geneInfo ;
  Array aaM = pmb->mrnaInfo ;

  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      const char *ccp ;
      int tc, ln, gc, gene, geneId, mrna ;

      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#')
	continue ;

      mrna = gene = geneId = tc = gc = ln = 0 ;
      if (0)
	{
	  if (ccp && ccp[0] != '-')
	    dictAdd (mrnaDict, ccp, &tc) ;
	  aceInStep (ai, '\t') ;
	  ccp = aceInWord (ai) ;
	}
      if (! ccp || *ccp == '#' || *ccp == '-')
	continue ;
      dictAdd (mrnaDict, ccp, &mrna) ;
      aceInStep (ai, '\t') ;
      aceInInt (ai, &ln) ;
      aceInStep (ai, '\t') ;
      aceInInt (ai, &gc) ;
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (ccp && strcmp (ccp, "-"))
	dictAdd (geneDict, ccp, &gene) ;
      aceInStep (ai, '\t') ;
      aceInInt (ai, &geneId) ;

      if (mrna && gene)
	{
	  GINFO *up = arrayp (aaG, gene, GINFO) ;
	  up->gene = gene ;
	  up->geneId = geneId ;
	  if (ln > up->ln) up->ln = ln ;
	  up->gc = gc ;	  
	  up->target_class = tc ;	  

	  up = arrayp (aaM, mrna, GINFO) ;
	  up->gene = gene ;
	  up->ln = ln ;
	  up->gc = gc ;	  
	  up->target_class = tc ;	  
	}
    }

  ac_free (h) ;
  return ;
} /* pmbParseGeneInfo */

/*************************************************************************************/

static void pmbParseSplitMrna (PMB *pmb)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (pmb->splitMrnaFileName, FALSE, h) ;

  ac_free (ai) ;
  ac_free (h) ;
  return ;
} /* pmbParseSplitMrna */

/*************************************************************************************/

static void pmbParseSam (ACEIN ai, PMB *pmb)
{
  DICT *cloneDict = pmb->cloneDict ;
  DICT *targetDict = pmb->targetDict ;
  DICT *cigarDict = pmb->cigarDict ;
  ALI *up ;
  int clone, flag, target, a1, cigar, ln, mult ;
  int target_class = pmb->target_class ;
  BigArray alis = pmb->alis ;
  long int nn = bigArrayMax (alis) ; ;
  int hasPair
    , read1Reversed
    , read2Reversed
    , isRead2
    ;

  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      const char *ccp, *ccq ;

      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#' || *ccp == '@')
	continue ;
      dictAdd (cloneDict, ccp, &clone) ;

      mult = 0 ;
      ccq = ccp + strlen (ccp) - 1 ;
      while (ccq > ccp && *ccq != '#')
	ccq-- ;
      if (*ccq == '#')
	{
	  ccq++ ;
	  while (*ccq >= '0' && *ccq <= '9')
	    { mult = 10 * mult + (*ccq - '0') ; ccq ++ ; }
	}
      if (mult == 0)
	mult = 1 ;

      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &flag) || !flag)
	continue ;

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#' || *ccp == '@')
	continue ;
      dictAdd (targetDict, ccp, &target) ;

      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &a1) || a1 < 1)
	continue ;

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* mapQual ignore */
      
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#' || *ccp == '@')
	continue ;
      dictAdd (cigarDict, ccp, &cigar) ;

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* mateTarget */

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* matePos */

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* templateLength */

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* fasta */
      ln = strlen (ccp) ;

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* phredQual */


      up = bigArrayp (alis, nn++, ALI) ;
      up->clone = clone ;
      up->mult = mult ;
      up->flag = flag ;
      up->target = target ;
      up->target_class = target_class ;
      up->a1 = a1 ;
      up->a2 = up->a1 + ln ;
      up->x1 = 1 ;
      up->x2 = up->x1 + ln ;
      up->ali = up->x2 - up->x1 + 1 ;

      up->cigar = cigar ;

      while (aceInStep (ai, '\t') && (ccp = aceInWord (ai)))
	{
	  if (ccp && ! strncmp (ccp, "NH:i:", 5))   /* unicity */
	    up->uu = atoi (ccp+5) ;
	  
	  if (ccp && ! strncmp (ccp, "AS:i:", 5))
	    up->score = atoi (ccp+5) ;
	  
	  if (ccp && ! strncmp (ccp, "NM:i:", 5))
	    up->nErr = atoi (ccp+5) ;

	  if (ccp && ! strncasecmp (ccp, "TC:a:", 5))
	    dictAdd (pmb->target_classDict, ccp+5, &(up->target_class)) ;
	}

      pmb->hasPair = hasPair = flag & 0x1 ;
      if (hasPair)
	{
	  /*
	    read2Unali = flag & 0x8 ;
	  */
	  read2Reversed = flag & 0x20 ;

	  up->isRead1 = (flag & 0x40 ? TRUE : FALSE) ;
	  isRead2 = flag & 0x80 ;
	}

      read1Reversed = flag & 0x10 ;

      /*
	read1Unali = flag & 0x4 ;
	secondaryAli = flag & 0x100 ;
	badQuality = flag & 0x200 ;
	isDuplicate = flag & 0x400 ;
	otherAli = flag & 0x800 ;
      */

      if ((up->isRead1 && read1Reversed) || (isRead2 && read2Reversed))
	{ int a0 = up->a1 ; up->a1 = up->a2 ; up->a2 = a0 ; }
    }
  

  pmbPairScore (pmb) ;
} /* pmbParseSam */

/*************************************************************************************/

static int cigarMult (char *cp, int *np)
{ 
  int n = 0, k = 0 ;
  while (*cp >= '0' && *cp <= '9')
    { n = 10*n + *cp - '0' ; k++ ; cp++ ; }
  *np = n ;
  return k ;
} /* cigarMult */

/*************************************************************************************/

static int pmbTabularSmoke (PMB *pmb, ALI *up, vTXT txt)
{
  AC_HANDLE h = ac_new_handle () ;
  DICT *cloneDict = pmb->cloneDict ;
  DICT *snpDict = pmb->snpDict ;
  DICT *runDict = pmb->runDict ;
  DICT *targetDict = pmb->targetDict ;
  char *cp ;
  char *cigar = strnew (dictName (pmb->cigarDict, up->cigar), h) ;
  int n, nerr = 0 ;
  int a1, a2 ;
  int isDown = (up->a1 < up->a2 ? 1 : -1) ;
  BOOL isReadF = TRUE ;

  if (isDown != up->isRead1)
    isReadF = FALSE ;

  up->nintron = 0 ;
  a1 = up->a1 ;
  for (cp = cigar ; *cp ; cp++)
    {
      cp += cigarMult (cp, &n) ;
      a1 += isDown * n ;
      if (*cp == 0)
	break ;
      else if (*cp == '%') /* deletion */
	{
	  int nd = 0 ;
	  cp++ ;
	  cp += cigarMult (cp, &nd) ;
	  if (*cp != '%') /* deletion end */
	    messcrash ("Bad deletion sign ^ in %s\n", cigar) ;
	  if (pmb->getSnp)
	    {
	      SNP s, *snp ;
	      char type[8], suffix[8] ;
	      const char *nam1 ;
	      char nam2[8], nam3[8] ;
	      
	      memset (&s, 0, sizeof (SNP)) ;

	      type[0] = 'D';
	      type[1] = 'e' ;
	      type[2] = 'l' ;
	      type[3] = cp[0] ;
	      type[4] = 0 ;
	      s.da = -nd ; s.dx = nd ;
	      suffix[0] = cp[0] ;
	      suffix[1] = 0 ;
	      
	      nam1 = "Del" ;
	      nam2[0] = 'x' ; nam2[1] = cp[0] ; nam2[2] = 0 ;
	      nam3[0] = 'x' ; nam3[1]= 0 ;
	    

	      if (! isDown && suffix[0])
		suffix[0] = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)suffix[0]])] ; 
	      if (isDown)
		{ s.a1 = a1 - 1 ; s.a2 = a1 + nd ;}
	      else
		{ s.a1 = a1 + 1 ; s.a2 = a1 - nd ;}
	      
	      dictAdd (pmb->cloneDict, type, &s.type) ;
	      if (s.da) dictAdd (cloneDict, suffix, &s.suffix) ;

	      vtxtClear (txt) ;
	      vtxtPrintf (txt, "%s:%d:%s:%s:%s", dictName (targetDict, up->target), a1, nam1, nam2, nam3) ;

	      dictAdd (snpDict, vtxtPtr (txt), &s.name) ;
	      dictAdd (cloneDict, "MB", &s.method) ;
	      dictAdd (runDict, pmb->run, &s.run) ;
	      
	      if (isReadF)
		s.m = s.mp = up->mult ;
	      else
		s.m = s.mm = up->mult ;
	      s.c = s.m + s.w ;
	      
	      snp = bigArrayp (pmb->snps, s.name, SNP) ;
	      
	      s.c += snp->c ;
	      s.m += snp->m ;
	      s.w += snp->w ;
	      s.mp += snp->mp ;
	      s.wp += snp->wp ;
	      s.mm += snp->mm ;
	      s.wm += snp->wm ;


	      *snp = s ;


	    }

	}
      else if (*cp == '_') /* insertion */
	{
	  int ni = 0 ;
	  cp++ ;
	  cp += cigarMult (cp, &ni) ;
	  if (*cp != '_') /* deletion end */
	    messcrash ("Bad insertion sign _ in %s\n", cigar) ;
	  nerr += ni ;
	  a1 += isDown * ni ;
	}
      else if (*cp == '^') /* intron */
	{
	  int ni = 0 ;
	  char feet[6] ;
	  cp++ ;
	  feet[0] = 0 ;

	  if (*cp >= 'a'&& *cp <='z') /* feet are reported */
	    {
	      feet[0] = cp[0] ;
	      feet[1] = cp[1] ;
	      cp+= 2 ;
	      cp += cigarMult (cp, &ni) ;
	      feet[2] = '_' ;
	      feet[3] = cp[0] ;
	      feet[4] = cp[1] ;
	      feet[5] = 0 ;
	      cp+= 2 ; ni += 4 ;
	    }
	  else
	    cp += cigarMult (cp, &ni) ;


	  if (*cp != '^') /* intron end */
	    messcrash ("Bad intron sign % in %s\n", cigar) ;
	  up->nintron++ ;
	  a2 = a1 + isDown * ni ;
	  if (pmb->getIntrons)
	    {
	      char intronName [1024] ;
	      ITR *ip ;
	      int k ;
	      int b1 = a1, b2 = a2 ;
	      if (isDown == 1)
		{ b1 = a1 ; b2 = a2 - 1 ; }
	      else 
		{ b1 = a1 ; b2 = a2 + 1 ; }

	      if (! isReadF)   { int b0 = b1 ; b1 = b2 ; b2 = b0 ; }

	      intronName[1023] = 0 ;
	      sprintf (intronName, "%d__%d_%d", up->target, b1, b2) ;
	      if (intronName[1023])
		messcrash ("intronName length overflow because in file %s\n// the target name is too long: %s\n"
			   , pmb->inFileList ? pmb->inFileList : "stdin"
			   , dictName (pmb->targetDict, up->target)
			   ) ;
	      dictAdd (pmb->intronDict, intronName, &k) ;
	      ip = arrayp (pmb->introns, k, ITR) ;
	      ip->a1 = b1 ; ip->a2 = b2 ; ip->target = up->target ;
	      if (up->isRead1)
		{
		  ip->nr1 += up->mult ;
		  if (! strcmp(feet, "gt_ag"))
		    pmb->r1Down += up->mult ;
		  else if (! strcmp(feet, "ct_ac"))
		    pmb->r1Up += up->mult ;
		}
	      else
		{
		  ip->nr2 += up->mult ;
		  if (! strcmp(feet, "gt_ag"))
		    pmb->r2Down += up->mult ;
		  else if (! strcmp(feet, "ct_ac"))
		    pmb->r2Up += up->mult ;
		}

	      strcpy (ip->feet, feet) ;
	    }
	  a1 += isDown * ni ;
	}
      else if (strchr ("ATGCN-", cp[0]) && strchr ("ATGCN-", cp[1]))
	{ /* substitution or 1 base indel */
	  if (pmb->getSnp)
	    {
	      SNP s, *snp ;
	      char type[8], suffix[8] ;
	      const char *nam1 ;
	      char nam2[8], nam3[8] ;
	      
	      memset (&s, 0, sizeof (SNP)) ;
	      if (cp[0] != '-' && cp[1] != '-')
		{ 
		  type[0] = cp[0] ;
		  type[1] = '2' ;
		  type[2] = cp[1] ;
		  type[3] = 0 ;
		  s.da = s.dx = 0 ;
		  suffix[0] = 0 ;

		  nam1 = "Sub" ;
		  nam2[0] = cp[0] ; nam2[1] = 0 ;
		  nam3[0] = cp[1] ; nam3[1]= 0 ;
		}
	      else if (cp[0] == '-' && cp[1] != '-')
		{ 
		  type[0] = 'I' ;
		  type[1] = 'n' ;
		  type[2] = 's' ;
		  type[3] = cp[1] ;
		  type[4] = 0 ;
		  s.da = s.dx = 1 ;
		  suffix[0] = cp[1] ;
		  suffix[1] = 0 ;

		  nam1 = "Ins" ;
		  nam2[0] = 'x' ; nam2[1] = 0 ;
		  nam3[0] = 'x' ; nam3[1] = cp[1] ; nam3[2] = 0 ;
		}
	      else if (cp[0] != '-' && cp[1] == '-')
		{ 
		  type[0] = 'D';
		  type[1] = 'e' ;
		  type[2] = 'l' ;
		  type[3] = cp[0] ;
		  type[4] = 0 ;
		  s.da = -1 ; s.dx = 1 ;
		  suffix[0] = cp[0] ;
		  suffix[1] = 0 ;

		  nam1 = "Del" ;
		  nam2[0] = 'x' ; nam2[1] = cp[0] ; nam2[2] = 0 ;
		  nam3[0] = 'x' ; nam3[1]= 0 ;
		}

	      if (! isDown && suffix[0])
		suffix[0] = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)suffix[0]])] ; 
	      if (isDown)
		{ s.a1 = a1 - 1 ; s.a2 = a1 + 1 ; if (s.da == 1) s.a2-- ; }
	      else
		{ s.a1 = a1 + 1 ; s.a2 = a1 - 1 ; if (s.da == 1) s.a2++ ; }
	      
	      dictAdd (pmb->cloneDict, type, &s.type) ;
	      if (s.da) dictAdd (cloneDict, suffix, &s.suffix) ;

	      vtxtClear (txt) ;
	      vtxtPrintf (txt, "%s:%d:%s:%s:%s", dictName (targetDict, up->target), a1, nam1, nam2, nam3) ;

	      dictAdd (snpDict, vtxtPtr (txt), &s.name) ;
	      dictAdd (cloneDict, "MB", &s.method) ;
	      dictAdd (runDict, pmb->run, &s.run) ;
	      
	      if (isReadF)
		s.m = s.mp = up->mult ;
	      else
		s.m = s.mm = up->mult ;
	      s.c = s.m + s.w ;
	      
	      snp = bigArrayp (pmb->snps, s.name, SNP) ;
	      
	      s.c += snp->c ;
	      s.m += snp->m ;
	      s.w += snp->w ;
	      s.mp += snp->mp ;
	      s.wp += snp->wp ;
	      s.mm += snp->mm ;
	      s.wm += snp->wm ;


	      *snp = s ;
	    }
	  nerr++ ; cp++ ;
	  a1 += isDown ;
	}
      else
	messcrash ("Bad sub in %s\n", cigar) ;
    }

  up->nErr = nerr ;
  ac_free (h) ;
  return nerr ;
} /* pmbTabularSmoke */

/*************************************************************************************/

static void pmbParseTabular (ACEIN ai, PMB *pmb)
{
  DICT *cloneDict = pmb->cloneDict ;
  DICT *targetDict = pmb->targetDict ;
  DICT *target_classDict = pmb->target_classDict ;
  DICT *cigarDict = pmb->cigarDict ;
  ALI *up ;
  int  mult ;
  BigArray alis = pmb->alis ;
  long int nn = bigArrayMax (alis) ;

  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      char *ccp, *ccq ;
      int bonus = 0 ;

      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#' || *ccp == '@')
	continue ;

      up = bigArrayp (alis, nn++, ALI) ;
      up->isRead1 = TRUE ; 
      ccq = ccp + strlen (ccp) - 2 ;
      if (pmb->hasPair && (ccq[0] == '.' && (ccq[1] == '1' || ccq[1] == '2')))
	{ 
	  if (ccq[1] == '2')
	    { up->isRead1 = FALSE ; pmb->hasPair = TRUE ; }
	  ccq[0] = 0 ;
	}
      

      dictAdd (cloneDict, ccp, &(up->clone)) ;

      mult = 0 ;
      ccq = ccp + strlen (ccp) - 1 ;
      while (ccq > ccp && *ccq != '#')
	ccq-- ;
      if (*ccq == '#')
	{
	  ccq++ ;
	  while (*ccq >= '0' && *ccq <= '9')
	    { mult = 10 * mult + (*ccq - '0') ; ccq ++ ; }
	}
      if (mult == 0)
	mult = 1 ;
      up->mult = mult ;

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#' || *ccp == '@')
	continue ;
      dictAdd (targetDict, ccp, &(up->target)) ;

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* percent identity ignore */

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* notused */
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* notused */
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* notused */

      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &(up->x1)) || up->x1 < 1)
	continue ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &(up->x2)) || up->x2 < 1)
	continue ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &(up->a1)) || up->a1 < 1)
	continue ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &(up->a2)) || up->a2 < 1)
	continue ;

      if (up->x1 > up->x2)
	{ 
	  int a0 = up->a1 ; up->a1 = up->a2 ; up->a2 = a0 ; 
	  int x0 = up->x1 ; up->x1 = up->x2 ; up->x2 = x0 ; 
	}

      up->ali = up->x2 - up->x1 + 1 ;

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* notused */
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* notused */

      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &(up->score)) || up->score < 1)
	continue ;

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* query strand plus/minus */
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* ref strand plus/minus */

      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &(up->ln)))  /* read-length */
	continue ;

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* BTOP = tabular-cigar */
      if (!ccp)
	continue ;
      dictAdd (cigarDict, ccp, &(up->cigar)) ;

      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &(up->uu)))  /* unicity */
	continue ;

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* notused */
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* compartment */

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* left overhang */
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* right overhang */

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* mate ref */
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* mate start */

      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &(up->pairScore)))  /* pairScore */
	continue ;
      if (! up->pairScore)
	up->pairScore = up->score ;

      aceInStep (ai, '\t') ;
      up->target_class = pmb->target_class ;
      ccp = aceInWord (ai) ; /* target_class */
      if (ccp)
	dictAdd (target_classDict, ccp, &(up->target_class)) ;

      aceInStep (ai, '\t') ;
      aceInInt (ai, &(bonus)) ; 

      up->score += bonus + pmb->bonus ; 
      up->pairScore += bonus + pmb->bonus ; 

      nn++ ;
    }
  
} /* pmbParseTabular */

/*************************************************************************************/
/*************************************************************************************/

static void pmbParseClipAli (ACEIN ai, PMB *pmb)
{
  DICT *cloneDict = pmb->cloneDict ;
  DICT *targetDict = pmb->targetDict ;
  DICT *target_classDict = pmb->target_classDict ;
  DICT *cigarDict = pmb->cigarDict ;
  DICT *geneDict = pmb->geneDict ;
  DICT *dnaDict = pmb->dnaDict ;
  ALI *up ;
  BigArray alis = pmb->alis ;
  long int nn = bigArrayMax (alis) ;

  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      char *ccp, *ccq ;

      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#' || *ccp == '@')
	continue ;

      up = bigArrayp (alis, nn++, ALI) ;
      up->isRead1 = TRUE ; 
      ccq = ccp + strlen (ccp) - 1 ;
      if (ccq[0] == '>'  || ccq[0] == '<')
	{ 
	  pmb->hasPair = TRUE ; 
	  if (ccq[0] == '<')
	    { up->isRead1 = FALSE ;}
	  ccq[0] = 0 ;
	}

      dictAdd (cloneDict, ccp, &(up->clone)) ;

      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &(up->score)) || up->score < 1)
	continue ;
      up->score += pmb->bonus ;  /* line bonus is already counted in the aligner */
      up->pairScore = up->score ;
      
      aceInStep (ai, '\t') ;
      aceInInt (ai, &(up->mult)) ;
      
      aceInStep (ai, '\t') ;
      aceInInt (ai, &(up->ln)) ; /* length to align after clipping */
      
      aceInStep (ai, '\t') ;
      aceInInt (ai, &(up->ali)) ;
      
      aceInStep (ai, '\t') ;
      aceInInt (ai, &(up->x1)) ;
      
      aceInStep (ai, '\t') ;
      aceInInt (ai, &(up->x2)) ;
      
      aceInStep (ai, '\t') ;
      up->target_class = pmb->target_class ;
      ccp = aceInWord (ai) ; /* target_class */
      if (ccp && strcmp (ccp,"-"))
	dictAdd (target_classDict, ccp, &(up->target_class)) ;
      
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (ccp && strcmp (ccp, "-"))
	dictAdd (geneDict, ccp, &(up->gene)) ;
      
      aceInStep (ai, '\t') ;
      aceInInt (ai, &(up->uu)) ;

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (1 && ccp && !strncmp (ccp,"MRNA:",5))
	ccp += 5 ;
      if (! ccp || ! *ccp)
	{ up->score = 0 ; continue ; }
      dictAdd (targetDict, ccp, &(up->target)) ;
      
      aceInStep (ai, '\t') ;
      aceInInt (ai, &(up->a1)) ;
      
      aceInStep (ai, '\t') ;
      aceInInt (ai, &(up->a2)) ;
      
      aceInStep (ai, '\t') ;
      aceInInt (ai, &(up->nN)) ;
      
      aceInStep (ai, '\t') ;
      aceInInt (ai, &(up->nErr)) ;
      
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (ccp && strcmp (ccp, "-"))
	dictAdd (cigarDict, ccp, &(up->cigar)) ;
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (ccp && strcmp (ccp, "-"))
	dictAdd (cigarDict, ccp, &(up->targetCigar)) ;
      
      ccp = aceInWord (ai) ;
      if (ccp && strcmp (ccp, "-"))
	dictAdd (dnaDict, ccp, &(up->prefix)) ;
      ccp = aceInWord (ai) ;
      if (ccp && strcmp (ccp, "-"))
	dictAdd (dnaDict, ccp, &(up->suffix)) ;

#ifdef JUNK
      ccp = aceInWord (ai) ;
      if (ccp && strcmp (ccp, "-"))
	dictAdd (dnaDict, ccp, &(up->targetPrefix)) ;
      ccp = aceInWord (ai) ;
      if (ccp && strcmp (ccp, "-"))
	dictAdd (dnaDict, ccp, &(up->targetSuffix)) ;
      
      
      /* if pair == false, do not compute dPair, but if known reexport it */
      deltaPair = 0 ;
      aceInStep (ai, '-') ; aceInInt (ai, & deltaPair)  ; aceInStep (ai, '\t') ;
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
      if (ccp && ! strcasecmp (ccp, "chain") &&
	  aceInInt (ai, &chain))
	{
	  up->chain = chain ;  
	  aceInStep (ai, '\t') ; aceInInt (ai, &(up3->c1)) ;
	  aceInStep (ai, '\t') ; aceInInt (ai, &(up3->c2)) ;
	}
#endif
      

      nn++ ;
    }
} /* pmbParseClipAli */

/*************************************************************************************/

static void pmbInit (PMB *pmb)
{
  AC_HANDLE h = pmb->h ;
  pmb->cloneDict = dictHandleCreate (100000, h) ;
  pmb->targetDict = dictHandleCreate (1000, h) ;
  pmb->target_classDict = dictHandleCreate (16, h) ;
  pmb->cigarDict = dictHandleCreate (1000, h) ;
  pmb->runDict = dictHandleCreate (100, h) ;
  pmb->dnaDict = dictHandleCreate (100000, h) ;
  pmb->geneDict = dictHandleCreate (100000, h) ;
  pmb->intronDict = dictHandleCreate (1000, h) ;

  if (pmb->getExpression || pmb->geneInfoFileName)
    {
      pmb->geneInfo = arrayHandleCreate (100000, GINFO, h) ;
      pmb->mrnaInfo = arrayHandleCreate (100000, GINFO, h) ;
      pmb->geneCounts = bigArrayHandleCreate (100000, CNT, h) ; 
      pmb->mrnaCounts = bigArrayHandleCreate (100000, CNT, h) ; 
    }
  if (pmb->getSnp)
    {
      pmb->snps = bigArrayHandleCreate (100000, SNP, h) ; 
      pmb->snpDict = dictHandleCreate (100000, h) ;
    }

  pmb->alis = bigArrayHandleCreate (100000, ALI, h) ; 
  if (pmb->getIntrons)
    pmb->introns = arrayHandleCreate (100000, ITR, h) ; 

  dictAdd (pmb->target_classDict, "0_SpikeIn", 0) ;
  dictAdd (pmb->target_classDict, "1_DNASpikeIn", 0) ;
  dictAdd (pmb->target_classDict, "A_mito", 0) ;
  dictAdd (pmb->target_classDict, "B_rrna", 0) ;
  dictAdd (pmb->target_classDict, "D_transposon", 0) ;

  dictAdd (pmb->target_classDict, "ET_av", 0) ;
  dictAdd (pmb->target_classDict, "KT_RefSeq", 0) ;
  dictAdd (pmb->target_classDict, "QT_smallRNA", 0) ;

  dictAdd (pmb->target_classDict, "Z_genome", 0) ;

  dictAdd (pmb->target_classDict, "b_bacteria", 0) ;
  dictAdd (pmb->target_classDict, "v_virus", 0) ;
  dictAdd (pmb->target_classDict, "z_gdecoy", 0) ;

  if (pmb->targetClassName)
    dictAdd (pmb->target_classDict, pmb->targetClassName, &(pmb->target_class)) ;
} /* pmbInit */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// postMagicBlast: post treatment of the MagicBLAST output\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, Tom Madden and Greg Boratyn, NCBI, October 2020, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "// Given a fasta, fastc or fastq file, aligned using MagicBlast\n"
	    "//    import the results\n"
	    "//    export all useful post treatment we may want in .tsf format\n"
	    "//\n"
	    "// Syntax:\n"
	    "// postMagicBlast -run run [options]\n"
	    "//   -run MANDATORY name of the run being analyzed\n"
	    "// All other parameters are optional and may be specified in any order\n"
	    "//\n"
	    "// INPUT:\n"
	    "//   -i alignment_file   [default: stdin]\n"
	    "//     Single file name, or coma separated list of files\n"
	    "//   -gzi : input files are gzipped\n"
	    "//     Usefull if parsing stdin\n"
	    "//     Input files named *.gz are automatically gunzipped\n"
	    "//   -tabular [default] :  input file format\n"
	    "//     Input files were generated by magicblast -tabular\n"
	    "//   -clipali :\n"
	    "//     Input files were generated by magic/clipalign\n"
	    "//   -sam :\n"
	    "//     Input files are in SAM format\n"
	    "//     CAVEAT: The output is more complete when using tabular input\n"
	    "//   -forward : strand plus\n"
	    "//   -reverse : strand minus\n"
	    "//     single-end sequencing: Orientation of the reads\n"
	    "//     paired-end sequencing: Orientation of the first reads of the pairs\n"
	    "//     Illumina short read sequencing is usually -reverse, nanopore is usually -forward\n"
	    "//     If neither strand is specified, the program makes a guess based on introns and transcripts alignments\n"
	    "//   -target_class tc : target class KT_RefSeq Z_Genome ...\n"
	    "//     The default target class can be specified as a parameter on the command line\n"
	    "//     For each alignment in tabular format it is expected as column 26\n"
	    "//                        in SAM format, it is expected as TC:a:target_class_name\n"
	    "//     If the second letter is T, the target is recognized as a transcript and used to guess the strand of the reads\n"
	    "//     The target_class specific counts are reported in the statistics\n"
	    "//   -bonus <int>:  bonus value. i.e. 1, -10\n"
	    "//     The global bonus can be specified as a parameter on the command line\n"
	    "//     Foreach alignment in tabular format the bonus can be given as column 27\n"
	    "//                       in SAM format, the bonus can be given as BN:i:integer\n"
	    "//     The eventual global and line specific bonuses are added to the line score\n"
	    "//     Adding a bonus to the score of the alignment is useful i.e. to favor ribosomal over genome alignments\n"
	    "// OUTPUT:\n"
	    "//   -o output_file_prefix\n"
	    "//      all exported files will be named output_file_prefix.action\n" 
            "//   -gzo : gzip all output files\n"
	    "// ACTIONS:\n"
	    "//   -stats [default]: report alignment statistics in tsf format\n" 
	    "//   -no_stats : do not report alignment statistics\n"
	    "//   -introns : report introns seen in genome alignment in tsf format\n" 
	    "//   -expression : report transcripts and gene counts in tsf format\n"
	    "//      -geneInfo fileName : mandatory if -expression is requested\n"
	    "//         Tab delimited, 5 columns: target_class, transcript_name, length, gc_percent, gene, gene_id\n"
	    "//         Column 2 and 5 are mandatory, other columns can be void, columns 3,4,6 contain integers\n" 
	    "//         Lines staring with # are ignored\n"
	    "//         The transcript names must match the targets (subjects) given in the alignment file\n"
	    "//      -split split_MRNA_file_name : optional\n"
	    "//         A way to attribute parts of long transcripts to different genes, in test\n"
	    "// Caveat:\n"
	    "//   Caption lines starting with # are ignored\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  qcsummary -help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  PMB pmb ;
  AC_HANDLE h = 0 ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&pmb, 0, sizeof (PMB)) ;
  pmb.h = h ;
  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  pmb.gzi = getCmdLineBool (&argc, argv, "-gzi") ;
  pmb.gzo = getCmdLineBool (&argc, argv, "-gzo") ;

  if (getCmdLineBool (&argc, argv, "-forward"))
    pmb.strand = 1 ;
  if (getCmdLineBool (&argc, argv, "-reverse"))
    {
      if (pmb.strand == 1)
	{
	  fprintf (stderr, "Arguments -forward and -reverse are incompatible, please try postMagicBlast -help\n") ;
	  exit (1) ; 
	}
      pmb.strand = -1 ;
    }

  pmb.sam = getCmdLineBool (&argc, argv, "-sam") ;
  pmb.hasPair = getCmdLineBool (&argc, argv, "-pair") ;
  pmb.clipali = getCmdLineBool (&argc, argv, "-clipali") ;
  pmb.tabular = getCmdLineBool (&argc, argv, "-tabular") ;
  if (pmb.sam + pmb.tabular + pmb.clipali > 2)
    {
      fprintf (stderr, "Arguments  -clipali -sam and -tabular are incompatible, please try postMagicBlast -help\n") ;
      exit (1) ; 
    }
  pmb.getIntrons = getCmdLineBool (&argc, argv, "-introns") ;
  pmb.getSnp = getCmdLineBool (&argc, argv, "-snp") ;
  pmb.getExpression = getCmdLineBool (&argc, argv, "-expression") ;


  pmb.showStats = TRUE ; /* show stats by default */
  pmb.showStats = getCmdLineBool (&argc, argv, "-stats") ; /* accept the default valu */
  pmb.showStats = ! getCmdLineBool (&argc, argv, "-no_stats") ; /* do not show stats */
    
  getCmdLineOption (&argc, argv, "-run", &(pmb.run)) ;
  getCmdLineOption (&argc, argv, "-i", &(pmb.inFileList)) ;
  getCmdLineOption (&argc, argv, "-o", &(pmb.outFileName)) ;
  getCmdLineOption (&argc, argv, "-info", &(pmb.geneInfoFileName)) ;

  getCmdLineOption (&argc, argv, "-target_class", &(pmb.targetClassName)) ;
  getCmdLineInt (&argc, argv, "-bonus", &(pmb.bonus)) ;

  pmb.maxIntronLength = 1000000 ;
  getCmdLineInt (&argc, argv, "-maxIntronLn", &(pmb.maxIntronLength)) ;
  if (pmb.maxIntronLength < 0) 
    {
      fprintf (stderr, "Sorry, parameter -maxIntron must be positive, please try postMagicBlast -help\n") ;
      exit (1) ;
    }

  if (! pmb.run)
    {
      fprintf (stderr, "Sorry, missing parameter -run, please try postMagicBlast -help\n") ;
      exit (1) ;
    }
  if (argc > 1) usage (messprintf ("Unknown parameters %s", argv[1])) ;

  if (pmb.getExpression && ! pmb.clipali && ! pmb.geneInfoFileName)
    { 
      fprintf (stderr, "Missing -geneInfo parameter, please try postMagicBlast -help\n") ;
      exit (1) ;
    }

  fprintf (stderr, "// start: %s\n", timeShowNow()) ;
  pmbInit (&pmb) ;
  if (pmb.geneInfoFileName)
    pmbParseGeneInfo (&pmb) ;
  if (pmb.splitMrnaFileName)
    pmbParseSplitMrna (&pmb) ;


  if (1)
    {
      ACEIN ai = 0 ;
      int nn = 0 ;
      const char *error = 0 ;

      /* check that all files exist */
      if (pmb.inFileList)
	{
	  if (! aceInCheckList (pmb.inFileList, &error, h))
	    messcrash ("Bad parameter -f, %s"
		       ,  error
		       ) ;
	}
      else 
	pmb.inFileList = "-" ; /*  '-'  means stdin in aceInCreate */
 
      /* parse a validated list of files */
      while ((ai = aceInCreateFromList (ai, &nn, pmb.inFileList, pmb.gzi, h)))
	{
	  if (pmb.clipali)
	    pmbParseClipAli (ai, &pmb) ;
	  else if (pmb.sam)
	    pmbParseSam (ai, &pmb) ;
	  else  if (pmb.tabular)
	    pmbParseTabular (ai, &pmb) ;
	  else
	    { 
	      fprintf (stderr, "Please specify the input file format (-tabular, -sam or -clipali), please try postMagicBlast -help\n") ;
	      exit (1) ;
	    }
	}
      fprintf (stderr, "// parse done: %s, found %ld alis in %d files\n", timeShowNow(), bigArrayMax (pmb.alis), nn) ;
    }
    
    
  fprintf (stderr, "// select start: %s\n", timeShowNow()) ;
  pmbSelect (&pmb) ;
  fprintf (stderr, "// select done: %s\n", timeShowNow()) ;
  pmbCount (&pmb) ;
  fprintf (stderr, "// pmbCount done: %s\n", timeShowNow()) ;

  pmbClassCount (&pmb) ;
  
  if (pmb.getIntrons)
    pmbIntronsExport (&pmb) ;

  if (pmb.getExpression)
    {
      pmbGeneExpressionExport (&pmb) ;
      pmbMrnaExpressionExport (&pmb) ;
    }
  if (pmb.getSnp)
     pmbSnpsExpport (&pmb) ;

  {
    int mx ;
    messAllocMaxStatus (&mx) ;   
    fprintf (stderr, "// done: %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
  }

  ac_free (h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
} /* main */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

