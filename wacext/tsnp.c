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

 * SNP and InDel DETECTION

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

/* #define MALLOC_CHECK    */
/* #define ARRAY_CHECK  */
#include "ac.h"
#include "channel.h"
#include "query.h"
#include "bql.h"

/* snp chaining distance */

typedef struct tsnpCallerTable {
  AC_HANDLE h ;
  const char *outFileName ;
  BOOL gzi, gzo ;
  ACEIN ai ; 
  ACEOUT ao ;

  const char *target_class ;
  const char *target ;
  int  t1, t2 ;

  BOOL intron_only ;
  const char *referenceGenome ;
  const char *inFileOfFileList ;
  const char *inFileList ;
  const char *project ;
  const char *zone ;
  const char *filter ;
  const char *select ;
  const char *remap2genome ;
  const char *remap2genes ;
  BOOL countTypes ;
  BOOL force ;
 
  int minSnpCover, minSnpCount ;
  float minSnpFrequency ;
  int minIntron ;
  int maxIntron ;

  int max_threads ;
  int maxLanes ;

  int nAna ; /* default 4, number of analyzers launched in parallel wego */
  DICT *snpDict ;
  DICT *snpTypeDict ;
  DICT *geneDict ;
  DICT *runDict ;
  DICT *runStatDict ;
  DICT *targetDict ;
  DICT *chromDict ;
  DICT *target_classDict ;
  DICT *selectDict ;
  DICT *varTypeDict ;

  Array runs ;
  Array snps ;
  int snpDetected, snpExported ;
  KEYSET wdelta ;
  BOOL doubleReport ;

  BOOL mergeCounts ;
  BOOL dbReport ;
  BOOL dbSnpProfile ;
  BOOL dbTranslate ;
  BOOL dbGroupCount ;
  BOOL dbGGG ;
  BOOL dropMonomodal ;
  const char *dbName ;
  const char *wiggleDir ;
  AC_DB db ; 
  AC_TABLE runMetaData ;

  Array target2geneAtlas ;
  Array target2exonAtlas ;
  Array runProfiles ;
  BOOL makeWords ;

  int pure, high, mid, low, ref, wellCovered ; /* prevalence */
  CHAN *getLaneInChan, *getLaneOutChan, *analyzeChan ;
  CHAN *doneChan ;
} TSNP ;


typedef struct snpStruct {
  int snp ; /* self, needed after we sort the table */

  int type ; /* offset in typeDict <<= 3, then use or not split counts in bits 0x3 and slide  */
  int varType ; /* offset in varTypeDict */
  int target, a1, a2 ; /* target and position as read, can be edited by calling the database for the remapping */
  BOOL select ;

  int coverp, coverm, acoverp, acoverm ; /* coverage or mutant or wild type on plus/minus trand of template, using pairs */
  int mp, mm, wp, wm, awp, awm, amp, amm ; /* mutant, wild-type or donor, acceptor. p/m: read in forward or reverse direction */
  int topRun ;
  float minFrequency ; /* initialized at 1000 */
  float maxFrequency ; /* initialized at -1000 */
  float alleleFrequency ;

  int tested, rejected, notMeasurable, wellCovered ; /* prevalence */
  int pure, high, mid, low, ref ;
  int xpure, xhigh, xmid, xlow, xref ;
  int ppure, phigh, pmid, plow, pref ;
  float pureC, highC, midC, lowC, refC ;
  float xpureC, xhighC, xmidC, xlowC, xrefC ;
  float ppureC, phighC, pmidC, plowC, prefC ;
  int wiggleProblemHigh ;
  int wiggleProblemLow ;
  Array counts ;
} SNP ;

typedef struct runStruct {
  int run ; /* self, needed after we sort the table */
  int pure, high, mid, low, ref, tested, rejected, notMeasurable, wellCovered ; /* prevalence */
  int pureC, highC, midC, lowC, refC, testedC, rejectedC, notMeasurableC, wellCoveredC ; /* prevalence */
  int xpure, xhigh, xmid, xlow, xref, xtested, xrejected, xnotMeasurable, xwellCovered ; /* prevalence */
  int ppure, phigh, pmid, plow, pref, ptested, prejected, pnotMeasurable, pwellCovered ; /* prevalence */
  int ppureC, phighC, pmidC, plowC, prefC, ptestedC, prejectedC, pnotMeasurableC, pwellCoveredC ; /* prevalence */
  int *types ;
  int *typesR ;
  int *typesC ;
  KEYSET r2g ; /* list of groups g of which r is a member */
  KEYSET g2r ; /* list of runs in this group */
} RC ;

static BOOL tsnpGetMonomodal (TSNP *tsnp, AC_OBJ Snp);

/*************************************************************************************/
/*************************************************************************************/
 
static int tsnpSnpA1Order  (const void *va, const void *vb)
{
  const SNP *up = (const SNP *)va, *vp = (const SNP*)vb ;
  int n ;
  
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;

  return 0 ;
} /* tsnpSnpA1Order */

/*************************************************************************************/
#define MAXTYPE 6
/* TYPE up->counts(run *MAXTYPE + i)  ::  0/1: var-donor/accp, 2/3 :ref d/a, 4/5: wiggle d/a */
static int tsnpSnpParseOne (TSNP *tsnp, ACEIN ai)
{
  int nn = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  DICT *selectDict = tsnp->selectDict ;
  DICT *runDict = tsnp->runDict ;
  DICT *typeDict = dictHandleCreate (8, h) ;

  dictAdd (typeDict, "Var", 0) ;
  dictAdd (typeDict, "VarDonor", 0) ;
  dictAdd (typeDict, "VarAccep", 0) ;

  dictAdd (typeDict, "Ref", 0) ;
  dictAdd (typeDict, "RefDonor", 0) ;
  dictAdd (typeDict, "RefAccep", 0) ;

  while (aceInCard (ai))
    {
      char cutter, *cq, *cp = aceInWordCut (ai, "\t", &cutter) ;
      int snp, type, target, run, delta, subDelIns = 0 ;
      int a0, a1, a2, a3, da = 0 ;
      SNP *up ;

      if (! cp || *cp == '#' || *cp == '/')
	continue ;
    
      if (0)
	{ /* the tsf file has no zone */
	  if (tsnp->zone && strcasecmp (cp, tsnp->zone))
	    continue ;
	  
	  aceInStep (ai,'\t') ;
	  cp = aceInWordCut (ai, "\t", &cutter) ;
	  if (! cp || *cp == '#' || *cp == '/')
	    continue ;
	} 
      if (1)
	{
	  cq = strstr (cp, "__") ;
	  if (!cq)
	    continue ;
	  *cq = 0 ; cq += 2 ;
	  if (! dictFind (typeDict, cp, &type))
	    continue ;
	  cp = cq ;
	}
      if (selectDict && ! dictFind (selectDict, cp, 0))
	continue ;
      if (! strcmp (cp, "NC_045512:9693:DelIns_118_119::"))  /* HACK 2020_06_10, bad variant */
	continue ;

      if (dictAdd (tsnp->snpDict, cp, &snp))
	nn++ ; /* else it is already detected */
      da = subDelIns = 0 ;
      up = arrayp (tsnp->snps, snp, SNP) ; /* make room */

      subDelIns = 0 ;
      a1 = a2 = target = 0 ;
      if (1)
	{
	  char *cr ;
	  cr = strchr (cp, ':') ;
	  if (cr)
	    {
	      *cr++ = 0 ;
	      if (*cp)
		dictAdd (tsnp->targetDict, cp, &target) ;
	      cp = cr ;
	      cr = strchr (cp, ':') ;
	      da = 1 ;
	      if (cr)
		{
		  cr++ ;
		  if (cp)  /* cp point on the coordinate  target:type:coord: */
		    {
		      char cc = 0 ;
		      sscanf (cp, "%d%c", &a1, &cc) ;
		      if (cc != ':')
			continue ;
		      cp = cr ;
		      cr = strchr (cp, ':') ;
		      if (cr)
			{  /* expect :a:aT for a substitution of Del_12 */
			  cr++ ;
			  if (cp)     /* cp point on the coordinate  target:type:coord: */
			    {
			      if (! strncmp (cp, "Ins", 3))
				{
				  da = 0 ;
				  subDelIns = 6 ;
				  if (cp[3] == '_')
				    subDelIns = 2 ;
				}
			      else if (! strncmp (cp, "DelIns", 6))
				{
				  int kD = 0, kI = 0 ;
				  
				  cc = 0 ;
				  subDelIns = 7 ; da = 1 ; /* wild guess */
				  cp += 7 ;
				  if (sscanf (cq, "%d_%d%c", &kD, &kI, &cc) == 3 && cc == ':')
				    da = kD ;
				}
			      else 
				{
				  if (! strncmp (cp, "Sub", 3))
				    {
				      subDelIns = 1 ;
				      if (cp[3] == '_')
					subDelIns = 2 ;
				    }
				  if (! strncmp (cp, "Del", 3))
				    {
				      subDelIns = 4 ;
				      if (cp[3] == '_')
					subDelIns = 7 ;
				    }
				  cr = cp + 3 ;
				  if (*cr == '_')
				    {
				      cr++ ;
				      cp = cr ;
				      cc = 0 ;
				      sscanf (cp, "%d%c", &da, &cc) ;
				      if (cc != ':')
					continue ;
				    }
				}
			    }
			}
		      a2 = a1 + da + 1 ;
		    }
		}
	    }
	}

      aceInStep (ai, '\t') ;
      cp = aceInWordCut (ai, "\t", &cutter) ;
      
      if (! cp || *cp == '#' || *cp == '/')
	continue ;
      if (tsnp->mergeCounts)
	dictAdd (runDict, cp, &(run)) ;
      if (! dictFind (runDict, cp, &(run)))
	continue ;
      
      aceInStep (ai, '\t') ;  /* the format, expect iiii or ii */
      cp = aceInWordCut (ai, "\t", &cutter) ; /* the format of this tsf file */
      if (! cp || strncmp(cp,"ii",2))
	continue ;
      
      aceInStep (ai,'\t') ;
      aceInInt (ai, &a0) ;   /* forward counts */
      aceInStep (ai,'\t') ;
      if (!aceInInt (ai, &a1))
	continue ;  /* reverse counts (direction of the read, not of the template as in pair sequencing */
      /* in the next 2 coulumns, we optionally have access to the stand coverage */
      aceInStep (ai,'\t') ;
      aceInInt (ai, &a2) ;   /* forward counts */
      aceInStep (ai,'\t') ;
      aceInInt (ai, &a3) ;  /* reverse counts (direction of the read, not of the template as in pair sequencing */

      up = arrayp (tsnp->snps, snp, SNP) ;
      up->snp = snp ; /* self */
      up->a1 = a1 ;
      up->a2 = a2 ;
      up->target = target ;
      if (! up->counts)
	{
	  up->counts = arrayCreate (MAXTYPE * (8+dictMax (runDict)), KEY) ;

	  up->minFrequency = 1000 ;
	  up->maxFrequency = -1000 ;
	}
      

      up->type &= 3 ;
      switch (type)
	{
	case 1: /* Variant */
	  up->type |= 0x1 ;   /* split var counts */
	  keySet (up->counts, MAXTYPE * run + 1) = a0 + a1 ; /* double register */
	  delta = 0 ;
	  break ;
	case 2: /* VarDonor*/
	  delta = 0 ;
	  break ;
	case 3: /* VarAccep */
	  delta = 1 ;
	  break ;

	case 4: /* Ref */
	  up->type |= 0x2 ;   /* split var counts */
	  keySet (up->counts, MAXTYPE * run + 3) = a0 + a1 ; /* double register */
	  delta = 2 ;
	  break ;
	case 5: /* RefDonor */
	  delta = 2 ;
	  break ;
	case 6: /* RefAccep */
	  delta = 3 ;
	  break ;
	}
      up->type |= ((up->type & 0x3) +  8 * subDelIns) ;
      
      keySet (up->counts, MAXTYPE * run + delta) = a0 + a1 ; /* merge the strands */

      switch (type)
	{                      /* coverage = cumul of all runs */
	case 1: /* Variant */
	  up->amp += a0 ;
	  up->amm += a1 ;
	  up->acoverp += a2 ;
	  up->acoverm += a3 ;
	  /* fallthru */
	case 2: /* VarDonor*/
	  up->mp += a0 ;
	  up->mm += a1 ;
	  up->coverp += a2 ;
	  up->coverm += a3 ;
	  break ;
	case 3: /* VarAccep */
	  up->amp += a0 ;
	  up->amm += a1 ;
	  up->acoverp += a2 ;
	  up->acoverm += a3 ;
	  break ;
	case 4: /* Ref */
	  up->awp += a0 ;
	  up->awm += a1 ;
	  up->acoverp += a2 ;
	  up->acoverm += a3 ;
	  /* fallthru */
	case 5: /* RefDonor */
	  up->wp += a0 ;
	  up->wm += a1 ;
	  up->coverp += a2 ;
	  up->coverm += a3 ;
	  break ;
	case 6: /* RefAccep */
	  up->awp += a0 ;
	  up->awm += a1 ;
	  up->acoverp += a2 ;
	  up->acoverm += a3 ;
	  break ;
	}
    }

  ac_free (h) ;
  return nn ;
} /* tsnpSnpParseOne */

/*************************************************************************************/
/*************************************************************************************/

/* import the runs metadata from the database */

static int tsnpGetSelectedSnps (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tbl = 0 ;
  const char *errors = 0 ;
  const char *ccp ;
  char *cp, *cq, *sep ;
  int ir, nn = 0, target ;
  DICT *targetDict = tsnp->targetDict ;
  DICT *snpDict = tsnp->snpDict ;
  vTXT txt = vtxtHandleCreate (h) ;

  if (tsnp->select)
    {
      cp = strnew (tsnp->select, h) ;
      sep = " where (" ;
      while (cp)
	{
	  cq = strchr (cp, ',') ;
	  if (cq) *cq = 0 ;
	  vtxtPrintf (txt, "%s s#%s ", sep, cp) ;
	  sep = " or " ;
	  if (cq) 
	    cp = cq + 1 ;
	  else
	    cp = 0 ;
	}
      vtxtPrintf (txt, " ) ") ;
    }

  tbl = ac_bql_table (tsnp->db, hprintf (h, "select s, c, a1, a2 from s in ?variant %s , c in s->intmap, a1 in c[1], a2 in c[2] where a2", vtxtPtr (txt)), 0, 0, &errors, tsnp->h) ;

  if (tbl)
    for (ir = 1 ; ir < tbl->rows ; ir++)
      {
	SNP *up ;
	int snp ;
	int a1 = ac_table_int (tbl, ir, 2, 0) ;
	int a2 = ac_table_int (tbl, ir, 3, 0) ;
	ccp = ac_table_printable (tbl, ir, 0, 0) ;
	dictAdd (snpDict, ccp, &snp) ;
	ccp = ac_table_printable (tbl, ir, 1, 0) ;
	dictAdd (targetDict, ccp, &target) ;
	up = arrayp (tsnp->snps, snp, SNP) ;
	up->snp = snp ;
	up->target = target ;
	up->a1 = a1 ;
	up->a2 = a2 ;
	up->select = TRUE ;
	nn++ ;
      }

  ac_free (h) ;
  return nn ;
} /* tsnpGetSelectedSnps */

/*************************************************************************************/

static int hisMutIsMyMut (SNP *up, SNP *vp)
{
  if (up == vp)
    return 0 ;
  if (1) /* his sides cover my mut A, so he is +B or ++, and his +B and ++ counts reinforce my + count */
    {
      int a1 = up->a1 ;
      int a2 = up->a2 ;
      int b1 = vp->a1 ;
      int b2 = vp->a2 ;
      int aT = up->type / 8 ;
      int bT = vp->type / 8 ;

      if (aT == 1)
	{
	  switch (bT)
	    {
	    case 2: /* my single sub is inide a muti sub, they are most likely identical multi sub */
	      if (a1 >= b1 && a1 <= b2)
		return 0 ;
	    case 1:
	      if (a1 <= b1 + 14 && a2 >= b2 - 14)
		return 1 ; /* his mut is my ref */
	      break ;
	    case 3:
	    case 4:
	    case 5:
	    case 6:
	    case 7: /* base vp->a1 and vp->a2 are correct */
	      if (bT < 7  || b2 - b1 < 14)
		{
		  if (a2 <= b1 && a2 >= b1 - 14)
		    return 1 ; /* his mut is my ref */
		  if (a1 <= b2 + 14 && a1 >= b2)
		    return 1 ; /* his mut is my ref */
		}
	      break ;
	    }
	}
      if (aT == 2)
	{
	  switch (bT)
	    {
	    case 2: /* my single sub is inide a muti sub, they are most likely identical multi sub */
	    case 1:
	    case 3:
	    case 4:
	    case 5:
	    case 6:
	      if (b2 - b1 < 14 && 
		  (
		   (a1 > b2 && a2 < b2 + 14 ) ||
		   (a2 < b1 && a1 > b1 - 14 ) 
		   )
		  )
		return 1 ; /* his mut is my ref */
	      break ;
	    case 7:
	      break ;
	    }
	}
      if (aT >=3 && aT <= 7 && a2 - a1 < 14)
	{
	  switch (bT)
	    {
	    case 1:
	    case 2:
	    case 5:
	    case 6:
	    case 7: /* base vp->a1 and vp->a2 are correct */
	      if (a1 >= b2 - 14 && a2 <= b1 + 14)
		return 1 ; /* his mut is my ref */
	      break ;
	    case 3:
	    case 4:
	      if (a2 >= b1 - 14 && a2 <= b1)
		return 1 ; /* his mut is my ref */
	      if (a1 >= b2 && a1 <= b2 + 14)
		return 1 ; /* his mut is my ref */
	      break ;
	    }
	}
    }

  return 0 ; /* remainder, for example a complex multisub */
}

/*************************************************************************************/

static void tsnpSnpParse (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;

  tsnp->targetDict = dictHandleCreate (64, tsnp->h) ;
  tsnp->snpDict = dictHandleCreate (1000, tsnp->h) ;
  if (tsnp->inFileOfFileList)
    {
      vTXT txt = 0 ;
      char *sep = "" ;
      ACEIN ai = aceInCreate (tsnp->inFileOfFileList, 0, h) ;
      if (!ai)
	messcrash ("Cannot open -file_of_file_list %s\n", tsnp->inFileOfFileList) ;

      txt = vtxtHandleCreate (h) ;
      while (aceInCard (ai))
	{
	  char cutter, *cp ;
	  while ((cp = aceInWordCut (ai, " ,\t", &cutter)))
	    {
	      if (cp && *cp && *cp != '#')
		{
		  vtxtPrintf (txt, "%s%s",sep, cp) ;
		  sep = "," ;
		}
	    }
	}
      if (*sep)
	tsnp->inFileList = vtxtPtr (txt) ;
    }

  tsnp->runs = arrayHandleCreate (300, RC, tsnp->h) ;
  tsnp->snps = arrayHandleCreate (30000, SNP, tsnp->h) ;
  if (tsnp->inFileList)
    {
      ACEIN ai = 0 ;
      int nn = 0 ;
      const char *error = 0 ;

      /* check that all files exist */
      if (! aceInCheckList (tsnp->inFileList, &error, h))
	messcrash ("Bad parameter -f, %s"
		       ,  error
		       ) ;	
 
      /* parse a validated list of files */
      while ((ai = aceInCreateFromList (ai, &nn, tsnp->inFileList, tsnp->gzi, h)))
	{
	  tsnp->snpDetected += tsnpSnpParseOne (tsnp, ai) ;
	}
     }
  if (tsnp->select) tsnpGetSelectedSnps (tsnp) ;
  arraySort (tsnp->snps, tsnpSnpA1Order) ;

  /* coallesce the neighbours */
  if (arrayMax (tsnp->snps) > 1)
    {
      int snp, snpMax = arrayMax (tsnp->snps) ;
      SNP *up, *vp ;
      int ir, irMax = dictMax (tsnp->runDict) ;
      
      for (snp = 1, up = arrp (tsnp->snps, snp, SNP) ; snp < snpMax ; up++, snp++)
	{
	  if (up->counts) /* count problems */
	    for (ir = 1 ; ir <= irMax ;ir++)
	      {
		int s ;
		int a1 = up->a1 ;
		int a2 = up->a2 ;
		KEY *kp = arrayp (up->counts, ir * MAXTYPE, KEY) ;
		int m5 = kp[0] ; /* mutant */
		int m3 = kp[1] ; /* mutant */
		int r5 = kp[2] ; /* ref */
		int r3 = kp[3] ; /* ref */
		
		for (s = snp - 30, vp = up - 30 ; s < snp + 30 && s < snpMax ; s++, vp++)
		  {
		    if (s <= 0)
		      continue ;
		    if (vp->counts && vp->a1 >= a1 - 14 && vp->a1 <= a1 +14)
		      {
			KEY *kp2 = arrayp (vp->counts, ir * MAXTYPE, KEY) ;
			int m52 = kp2[0] ; /* variant, including self */
			int r52 = kp2[2] ; /* ref, including self */
			switch (hisMutIsMyMut (up, vp))
			  {
			  case 2:
			    m5 += m52 ;
			    break ;
			  case 1:
			    if (r5 < m52 + r52) r5 = m52 + r52 ; /* all the counts of my neighbous are my wild type, unless they are my mutant */
			    break ;
			  }
		      }
		    if (vp->counts && vp->a1 >= a2 - 14 && vp->a2 <= a1 +14)
		      {
			KEY *kp2 = arrayp (vp->counts, ir * MAXTYPE, KEY) ;
			int m32 = kp2[1] ; /* variant, including self */
			int r32 = kp2[3] ; /* variant, including self */
			switch (hisMutIsMyMut (up, vp))
			  {
			  case 2:
			    m3 += m32 ;
			    break ; 
			  case 1:
			    if (r3 < m32 + r32) r3 = m32 + r32 ; /* all the counts of my neighbous are my wild type, unless they are my mutant */
			    break ;
			  }
		      }
		  }
		kp[4] = (m5 > r5 ? 2*m5 : 2*r5+1) ;  /* union with future wiggle value */
		kp[5] = (m3 > r3 ? 2*m3 : 2*r3+1) ;  /* if the neighbour is wild, i am wild */
	      }
	}
    }
  /* accept the coallesced values */
  if (arrayMax (tsnp->snps) > 1)
    {
      int snp, snpMax = arrayMax (tsnp->snps) ;
      SNP *up ;
      int ir, irMax = dictMax (tsnp->runDict) ;
      
      for (snp = 1, up = arrp (tsnp->snps, snp, SNP) ; snp < snpMax ; up++, snp++)
	{
	  if (up->counts) /* count problems */
	    for (ir = 1 ; ir <= irMax ;ir++)
	      {
		KEY *kp = arrayp (up->counts, ir * MAXTYPE, KEY) ;
		int m5 = kp[0] ; /* mutant */
		int m3 = kp[1] ; /* mutant */
		int r5 = kp[2] ; /* ref */
		int r3 = kp[3] ; /* ref */
		int w5 = kp[4] ; /* ref */ if (w5 & 0x1) w5 = - (w5-1)/2 ; else w5 = w5/2 ;
		int w3 = kp[5] ; /* ref */ if (w3 & 0x1) w3 = - (w3-1)/2 ; else w3 = w5/2 ;
		if (w5 > m5) { kp[0] = w5 ;}
		if (-w5 > r5) { kp[2] = -w5 ;}
		if (w3 > m3) { kp[1] = w3 ;}
		if (-w3 > r3) { kp[3] = -w3 ;}
		kp[4] = m3 + r3 ; /* corrected estimation of the donor site coverage */
		kp[5] = m5 + r5 ; /* acceptor = mutant-variant + wild-reference + others */
	      }
	}
    }
		
  ac_free (h) ;

  return  ;
} /* tsnpSnpParse */

/*************************************************************************************/
/*************************************************************************************/

static void tsnpWiggleParseOne (TSNP *tsnp, int run, int target)
{
  AC_HANDLE h = ac_new_handle () ;
  char *fNam = hprintf (h, "%s/%s/%s/R.chrom.frns.u.BF.gz", tsnp->wiggleDir, dictName (tsnp->runDict, run), dictName (tsnp->targetDict, target)) ;
  ACEIN ai = aceInCreate (fNam, 0, h) ;
  int line = 0, start = 0, step = 0, x = 0 ;
  float w ;
  Array ww = 0 ;
  char *cp ;

  if (ai)
    { /* parse the whole file, it is simpler and unavoidable anyway. We get a tmp sorted wiggle, we can then access it from the unsorted snp array */
      line++ ;
      ww = arrayHandleCreate (1000000, float, h) ;
      while (aceInCard (ai))
	{
	  if (!step)
	    {
	      if (line > 12) break ; /* badly formatted file */
	      cp = aceInWord (ai) ;
	      if (!cp || strcmp (cp, "fixedStep"))
		continue ;
	      cp = aceInWord (ai) ;
	      if (!cp)
		continue ;
	      cp = strstr (cp,"chrom=") ;
	      if (!cp)
		continue ;
	      cp += 6 ;
	      if (!cp || strcmp (cp, dictName (tsnp->targetDict, target)))
		continue ;
	      
	      cp = aceInWord (ai) ;
	      if (!cp)
		continue ;
	      cp = strstr (cp,"start=") ;
	      if (!cp)
		continue ;
	      cp += 6 ;
	      if (sscanf (cp, "%d", &start) < 1 || start<= 0)
		continue ;
	      
	      cp = aceInWord (ai) ;
	      if (!cp)
		continue ;
	      cp = strstr (cp,"step=") ;
	      if (!cp)
		continue ;
	      cp += 5 ;
	      if (sscanf (cp, "%d", &step) < 1 || step<= 0)
		continue ;
	      continue ;
	    }
	  x++ ; w = 0 ;
	  if (aceInFloat (ai, &w))
	    array (ww, x, float) = w ;	  
	}
    }

  if (ww && arrayMax (ww)) /* we have a wiggle, with known start and step, we scan the snps and update the wiggle values at a1 and a2 */
    {
      int kk, pass, snp, snpMax = arrayMax (tsnp->snps) ;
      int zz[2] ;
      int ww0[2] ;
      SNP *up ;
      int stop = start + step * arrayMax (ww) ;
      float wwww[5], w0, dw, dw0, z ;
      if (!step) step = 1 ;
      for (snp = 1, up = arrp (tsnp->snps, 0, SNP) ; snp < snpMax ; up++, snp++)
	{
	  KEY *kp = 0 ;
	  if (! up->select && tsnp->minSnpFrequency && (up->maxFrequency < 0 || up->maxFrequency < tsnp->minSnpFrequency))
	    continue ;
	  if (up->a1 < start+10 || up->a1 > stop - 10)
	    continue ;
	  if (! up->counts)
	    {
	      up->counts = arrayCreate (MAXTYPE * (1+dictMax (tsnp->runDict)), KEY) ;
	      
	      up->minFrequency = 1000 ;
	      up->maxFrequency = -1000 ;
	    }
	  kp = arrayp (up->counts, run * MAXTYPE, KEY) ; 
	  memset (wwww, 0, sizeof(wwww)) ;
	  ww0[0] = kp[4] ; /* was set as sum of negighbours */
	  ww0[1] = kp[5] ; /* was set as sum of neighbours */
	  for (kk = 0 ; kk < 2 ; kk++) /* donor acceptor */
	    {
	      int i, a1, a01 = kk == 0 ? up->a1 : up->a2 ;
	      if (! kk || up->a1 >= up->a2 + 3 || up->a1 <= up->a2 - 3)  /* reuse the same wwww */
		{
		  memset (wwww, 0, sizeof(wwww)) ;
		  for (pass = 0 ; pass < 5 ; pass++)
		    {
		      switch (pass)
			{
			case 0: a1 = a01 ; break ;
			case 1: a1 = a01 - 1 ; break ;
			case 2: a1 = a01 + 1 ; break ;
			case 3: a1 = a01 - 10 ; break ;
			case 4: a1 = a01 + 10 ; break ;
			}
		      if (a1)
			{
			  int b1 = (a1 - start)/step ;
			  float db = a1 - start - step * b1 ; db /= 1.0 * step ; /* how much we need of next position */
			  if (db >= 1) 
			    messcrash ("a1=%d b1=%d step=%d start=%d db=%f should be < 1, file %s", up->a1, step*b1, step, start, db, fNam) ;
			  if (b1 < arrayMax (ww) - 1)
			    {
			      z = (1 - db) * arr (ww, b1, float) + db * arr (ww, b1+1, float) ;
			      wwww[pass] = z ;
			    }
			}
		    }
		}
	      /* select best estimate, closest to the word counts */
	      if (! up->counts)
		up->counts = arrayCreate (MAXTYPE * (1+dictMax (tsnp->runDict)), KEY) ;
	      w0 = ww0[kk] ;
	      dw0 = 1000000000 ; /* horrible value */
	      z = 0 ;
	      for (i = 0 ; i < 5 ; i++) /* try direct positions */
		{
		  dw = wwww[i] - w0 ; if (dw < 0) dw = -dw ; if (dw0 > dw) { dw0 = dw ; z = wwww[i] ; }
		}
	      for (i = 1 ; i < 4 ; i += 2) /* try averages */
		{
		  dw = (wwww[i] + wwww[i+1])/2 - w0 ; if (dw < 0) dw = -dw ; if (dw0 > dw) { dw0 = dw ; z = (wwww[i] + wwww[i+1])/2 ; }
		}
	      zz[kk] = z ;
	    }
	  z = zz[0] + zz[1] ;
	  w0 = ww0[0] + ww0[1] ;
	  if (z > 20 && z > 2 * w0) 
	    { zz[0] = - zz[0] ; up->wiggleProblemHigh++ ; }
	  if (w0 > 20 && 2 * z < w0)
	    { zz[0] = w0/2 ; zz[1] = -z/2 ;  up->wiggleProblemLow++ ; }
	  kp[4] = zz[0] ; /* register best estimate and replace count gathering */
	  kp[5] = zz[1] ; /* register best estimate and replace count gathering */
	}
    }
  
  ac_free (h) ;
  return  ;
} /* tsnpWiggleParseOne */

/*************************************************************************************/

static void tsnpWiggleParse (TSNP *tsnp)
{
  int run, runMax = dictMax (tsnp->runDict) ;
  int target, targetMax = dictMax (tsnp->targetDict) ;

  if (! tsnp->wiggleDir)
    return ; 

  tsnp->wdelta = keySetHandleCreate (tsnp->h) ;
  for (run = 1 ; run <= runMax ; run++)
    for (target = 1 ; target <= targetMax ; target++)
      tsnpWiggleParseOne (tsnp, run, target) ;

  if (keySetMax (tsnp->wdelta))
    {
      AC_HANDLE h = ac_new_handle () ;
      ACEOUT ao = aceOutCreate (tsnp->outFileName, ".snp_wiggle_delta.txt", 0, h) ;
      int i, iMax = keySetMax (tsnp->wdelta) ;

      aceOutDate (ao, "###", "Histo of delta frequency relative to the wiggles") ;
      aceOutf (ao, " ## delta = 100 * w / (v + r), where w is the wiggle coverage and v and r and the reported numbers of 31-mers representing at the same position the reference (r) and the variant (v)\n") ;
      aceOut (ao, "# delta\tNumber of occurences\n") ;
      for (i = 0 ; i <= iMax ; i++)
	aceOutf (ao, "%d\t%d\n", i, keySet (tsnp->wdelta, i)) ;
      ac_free (h) ;
    }

  return  ;
} /* tsnpWiggleParse */


/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

typedef struct remapStruct { int target_class, gene, target, x1, x2, chrom, a1, a2 ; } RMP ;

/*****************/

static int atlasOrder1 (const void *va, const void *vb)
{
  const RMP *a = (const RMP *)va, *b = (const RMP *)vb ;
  int n ;
  n = a->target_class - b->target_class ; if (n) return n ;
  n = a->target - b->target ; if (n) return n ;
  n = a->x1 - b->x1 ; if (n) return n ;
  n = a->x2 - b->x2 ; if (n) return n ;
  n = a->chrom - b->chrom ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;

  return 0 ;
} /* atlasOrder1 */

/*****************/

static int atlasOrder2 (const void *va, const void *vb)
{
  const RMP *a = (const RMP *)va, *b = (const RMP *)vb ;
  int n ;
  n = a->chrom - b->chrom ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  n = a->target_class - b->target_class ; if (n) return n ;
  n = a->target - b->target ; if (n) return n ;
  n = a->x1 - b->x1 ; if (n) return n ;
  n = a->x2 - b->x2 ; if (n) return n ;

  return 0 ;
} /* atlasOrder1 */

/*************************************************************************************/
/* a tampon in bp, to allocate the snp to the promotor region of the gene */
#define PROMOTOR_DELTA 500

static int tsnpCreateAtlas (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai ;
  Array atlas = 0, geneAtlas = 0, map, geneMap ;
  const char *ccp ;
  RMP *up ;
  int nn = 0, x1, x2, a1, a2, target, chrom, gene, target_class ;

  x1 = x2 = a1 = a2 = 0 ;
  if (! tsnp->target_classDict)
    {
      tsnp->targetDict = dictHandleCreate (4, tsnp->h) ;
      tsnp->target_classDict = dictHandleCreate (4, tsnp->h) ;
      tsnp->chromDict = dictHandleCreate (64, tsnp->h) ;
      tsnp->target2geneAtlas = arrayHandleCreate (4, Array, tsnp->h) ;
      tsnp->target2exonAtlas = arrayHandleCreate (4, Array, tsnp->h) ;
    }

  if (tsnp->remap2genome)
    ai = aceInCreate (tsnp->remap2genome, FALSE, h) ;
  else
    ai = aceInCreate (tsnp->remap2genes, FALSE, h) ;
  aceInSpecial (ai, "\t\n") ;

  while (ai && aceInCard (ai))
    {
      ccp = aceInWord (ai) ; /* target_class */
      if (! ccp || *ccp == '#')
	continue ;
      if (tsnp->target_class && strcmp (ccp, (tsnp->target_class)))
	continue ;
      ccp = "any" ; /* ignore target class for the moment */
      dictAdd (tsnp->target_classDict, ccp, &target_class) ;
      atlas = array (tsnp->target2exonAtlas, target_class, Array) ;
      geneAtlas = array (tsnp->target2geneAtlas, target_class, Array) ;
      if (! atlas)
	{
	  atlas = arrayHandleCreate (100000, Array, tsnp->h) ;
	  array (tsnp->target2exonAtlas, target_class, Array) = atlas ;
	}
      if (! geneAtlas)
	{
	  geneAtlas = arrayHandleCreate (10000, Array, tsnp->h) ;
	  array (tsnp->target2geneAtlas, target_class, Array) = geneAtlas ;
	}
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* target */
      if (! ccp || *ccp == '#')
	continue ; 
      dictAdd (tsnp->targetDict, ccp, &target) ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &x1))  /* mRNA exon coordinates */
	continue ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &x2))  /* mRNA exon coordinates */
	continue ;
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* chrom */
      if (! ccp || *ccp == '#')
	continue ; 
      dictAdd (tsnp->chromDict, ccp, &chrom) ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &a1))  /* mRNA exon coordinates */
	continue ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &a2))  /* mRNA exon coordinates */
	continue ;
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* gene */
      if (! ccp || *ccp == '#')
	continue ;
      dictAdd (tsnp->targetDict, ccp, &gene) ;
      if (1)
	{
	  map = array (atlas,  tsnp->remap2genome ? target : chrom, Array) ;
	  if (! map)
	    map = array (atlas, tsnp->remap2genome ? target : chrom, Array) = arrayHandleCreate (8, RMP, tsnp->h) ;
	  geneMap = array (geneAtlas,  tsnp->remap2genome ? target : chrom, Array) ;
	  if (! geneMap)
	    geneMap = array (geneAtlas, tsnp->remap2genome ? gene : chrom, Array) = arrayHandleCreate (8, RMP, tsnp->h) ;
	  nn++ ;
	  up = arrayp (map, arrayMax (map), RMP) ;
	  up->target = target ;
	  up->chrom = chrom ;
	  if (tsnp->remap2genome || a1 < a2)
	    {
	      up->x1 = x1 ; up->x2 = x2 ;
	      up->a1 = a1 ; up->a2 = a2 ;
	    }
	  else
	    {
	      up->x1 = x2 ; up->x2 = x1 ;
	      up->a1 = a2 ; up->a2 = a1 ;
	    }

	  up =  arrayMax (geneMap) ? arrayp (geneMap, arrayMax (geneMap) - 1, RMP) : 0 ;
	  if (! up || up->target != gene)
	    up = arrayp (geneMap, arrayMax (geneMap), RMP) ;
	  up->target = gene ;
	  up->chrom = chrom ;
	  if (a1 < a2)
	    {
	      a1 -= PROMOTOR_DELTA ;
	      if (! up->a1 || a1 < up->a1) up->a1 = a1 ; 
	      if (a2 > up->a2) up->a2 = a2 ;
	      up->x1 = 1 ; up->x2 = up->a2 - up->a1 + 1 ;
	    }
	  else
	    {
	      a1 += PROMOTOR_DELTA ;
	      if (! up->a1 || a2 < up->a1) up->a1 = a2 ; 
	      if (a1 > up->a2) up->a2 = a1 ;
	      up->x2 = 1 ; up->x1 = up->a2 - up->a1 + 1 ;
	    }
	}
    }

  for (a1 = 0 ; atlas && a1 < arrayMax (atlas) ; a1++)
    {
      map = array (atlas, a1, Array) ;
      if (map)
	arraySort (map, tsnp->remap2genome ? atlasOrder1 : atlasOrder2) ;
    }
  for (a1 = 0 ; geneAtlas && a1 < arrayMax (geneAtlas) ; a1++)
    {
      map = array (geneAtlas, a1, Array) ;
      if (map)
	{
	  arraySort (map, tsnp->remap2genome ? atlasOrder1 : atlasOrder2) ;
	  arrayCompress (map) ;
	}
    }
   
  fprintf (stderr, "tsnpCreateAtlas found %d exons in file %s\n"
	   , nn
	   , tsnp->remap2genome ? tsnp->remap2genome :  tsnp->remap2genes
	   ) ;

  ac_free (h) ;
  return nn ;
} /* tsnpCreateAtlas */

/*************************************************************************************/
/* 
 * We import all the relevant data from the ZZ database
 * Remap the transcript variants into genome coordinates
 */

static BOOL tsnpRemap1Do (TSNP *tsnp, int mrna, int x1, int x2, int *chromp, int *a1p, int *a2p, int *strandp)
{
  int ii, target_class = 0 ;
  RMP *up ;
  Array atlas, map ;
  
  dictAdd (tsnp->target_classDict, "any", &target_class) ;
  if (tsnp->target2exonAtlas && 
      target_class < arrayMax (tsnp->target2exonAtlas) && 
      (atlas =  arr (tsnp->target2exonAtlas, target_class, Array)) &&
      mrna < arrayMax(atlas) &&
      (map = array (atlas, mrna, Array))
      )
    {
      for (ii = 0, up = arrp (map, 0, RMP) ; ii < arrayMax (map) ; ii++, up++)
	{
	  if (up->x1 <= x1 && up->x2 >= x1)
	    {
	      if (up->a1 < up->a2)
		{ *a1p = up->a1 + x1 - up->x1 ; *a2p = up->a1 + x2 - up->x1 ; *strandp = 1 ; *chromp = up->chrom ; return TRUE ; }
	      else
		{ *a1p = up->a1 - x1 + up->x1 ; *a2p = up->a1 - x2 + up->x1 ; *strandp = - 1 ; *chromp = up->chrom ; return TRUE ; }
	    }
	  else if (x1 < up->x1)
	    break ;
	}
    }
  return FALSE ;
} /* tsnpRemap1Do */

/*************************************************************************************/
/* brs_snp.tsf only contains the name of the SNP, the Counts and the tag Found_in_mRNA
 * here we add details deduced from the name
 */
static void tsnpRemap0 (TSNP *tsnp)
{
  AC_HANDLE  h1 = 0, h = ac_new_handle () ;
  AC_ITER iter ;
  AC_OBJ variant = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  int a1, a2, nn1 = 0, nn2 = 0 ;
  char *cp, *snpNam, *seqNam, *posBuf, *typeBuf, *Abuf, *Bbuf ;
  const char *errors = 0 ;
  
  if (tsnp->db)
    {
      if (tsnp->force)
	{
	  if (tsnp->remap2genome)
	    iter = ac_query_iter (tsnp->db, TRUE, "find variant Found_in_mRNA", 0, h) ;
	  if (tsnp->remap2genes)
	    iter = ac_query_iter (tsnp->db, TRUE, "find variant Found_in_genome", 0, h) ;
	}
      else
	{
	  if (tsnp->remap2genome)
	    iter = ac_query_iter (tsnp->db, TRUE, "find variant Found_in_mRNA && ! IntMap ", 0, h) ;
	  if (tsnp->remap2genes)
	    iter = ac_query_iter (tsnp->db, TRUE, "find variant Found_in_genome && ! mRNA ", 0, h) ;
	}
      while (ac_free (variant), ac_free (h1), variant = ac_iter_obj (iter))
	{
	  h1 = ac_new_handle () ;

	  nn1++ ;
	  snpNam = strnew (ac_name (variant), h1) ;
	  cp = seqNam = snpNam ;
	  cp = strchr (cp, ':') ; if (! cp) continue ;
          *cp++ = 0 ;
	  posBuf = cp ;
	  cp = strchr (cp, ':') ; if (! cp) continue ;
	  *cp++ = 0 ;
	  typeBuf = cp ;
	  a1 = 0 ;
	  if (sscanf (posBuf, "%d", &a1) != 1) continue ;
	  if (a1 < 1)  continue  ;
	  cp = strchr (cp, ':') ; if (! cp) continue ;
	  *cp++ = 0 ;
	  Abuf = cp ;
	  cp = strchr (cp, ':') ; if (! cp) continue ;
	  *cp++ = 0 ;
	  Bbuf = cp ;
	  cp = strchr (cp, ':') ; if (cp) *cp = 0 ; /* in case there are additional fields beyond target:234:Sub:A:G */
	  if (! seqNam[0]) continue ;
	  if (! Bbuf[0])   continue ;
	  if (! Abuf[0])   continue ;
	  if (Bbuf[0] == '~')   continue ;
	  if (Abuf[0] == '~')   continue ;

	  vtxtPrintf (txt, "\nVariant \"%s\"\n", ac_name (variant)) ;
	  vtxtPrintf (txt, "-D Type\n-D Seq\n-D VCF_hg38\n-D Site\n") ;
	  if (tsnp->remap2genes)
	    {
	      vtxtPrintf (txt, "Found_in_genome\n") ;
	      vtxtPrintf (txt, "Parent_sequence \"%s\"\n", seqNam) ;
	    }
	  if (tsnp->remap2genome)
	    {
	      vtxtPrintf (txt, "Found_in_mRNA\n") ;
	      vtxtPrintf (txt, "-D mRNA\nmRNA \"%s\"\n", seqNam) ;
	    }



	  if (!strcasecmp (typeBuf, "Sub")) 
	    {
	      a2 = a1 + 1 ; a1 = a1 - 1 ;
	      vtxtPrintf (txt, "Typ \"%c>%c\"\n", Abuf[0], Bbuf[0]) ;
	      vtxtPrintf (txt, "%c2%c\n", Abuf[0], Bbuf[0]) ;
	    }
	  else if (!strcasecmp (typeBuf, "Ins")) 
	    {
	      int k = strlen (Bbuf) ;
	      a2 = a1 + 1 ;
	      if (k > 0)
		vtxtPrintf (txt, "Typ Ins%s\n", Bbuf + 1) ;
	      if (k > 2)
		vtxtPrintf (txt, "Multi_insertion %d %s\n", k - 1, Bbuf+1) ;
	      else if (k > 1)
		vtxtPrintf (txt, "Ins%c\n", Bbuf[1]) ;
	    }
	  else if (!strcasecmp (typeBuf, "Del")) 
	    {
	      int k = strlen (Abuf) ;
	      a2 = a1 + k ;
	      if (k > 0)
		vtxtPrintf (txt, "Typ Del%s\n", Abuf + 1) ;
	      if (k > 2)
		vtxtPrintf (txt, "Multi_deletion %d %s\n", k - 1, Abuf+1) ;
	      else if (k > 1)
		vtxtPrintf (txt, "Del%c\n", Abuf[1]) ;
	    }
	  else
	    continue ;
	  nn2++ ;
	  if (tsnp->remap2genes)
	    {
	      vtxtPrintf (txt, "Site\nFound_in_genome\n") ;
	      vtxtPrintf (txt, "Parent_sequence \"%s\"\n", seqNam) ;
	      vtxtPrintf (txt, "IntMap \"%s\" %d %d\n", seqNam, a1, a2) ;
	    }
	  if (tsnp->remap2genome)
	    {
	      vtxtPrintf (txt, "Site\nFound_in_mRNA\n") ;
	      vtxtPrintf (txt, "mRNA \"%s\" %d %d\n", seqNam, a1, a2) ;
	    }

	  vtxtPrintf (txt, "\n") ;
	}
      ac_parse (tsnp->db, vtxtPtr (txt), &errors, 0, h) ; 
    }
  fprintf(stderr, "tsnpRemap0 found %d variants remapped %d\n", nn1, nn2) ;

  if (errors && *errors) fprintf(stderr, "tsnpRemap0 parsing error %s\n", errors) ;
  ac_free (h1) ;
  ac_free (h) ;
  return ;
} /* tsnpRemap0 */

/*************************************************************************************/
/* scan the VariantDB acedb database
 * add the remap info 
 * Remap the transcript variants into genome coordinates
 */
static int tsnpRemap1 (TSNP *tsnp)
{
  AC_HANDLE  h1 = 0, h = ac_new_handle () ;
  AC_ITER iter ;
  AC_OBJ variant = 0 ;
  AC_TABLE mrnaTable = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  int x1, x2, a1, a2, strand, chrom, mrna, nn1 = 0, nn2 = 0 ;
  const char *mm ;
  const char *errors = 0 ;
  
  if (tsnp->db)
    {
      if (tsnp->force)
	iter = ac_query_iter (tsnp->db, TRUE, "find variant Found_in_mRNA", 0, h) ;
      else
	iter = ac_query_iter (tsnp->db, TRUE, "find variant Found_in_mRNA && ! IntMap ", 0, h) ;
      while (ac_free (variant), ac_free (h1), variant = ac_iter_obj (iter))
	{
	  h1 = ac_new_handle () ;
	  nn1++ ;
	  mrnaTable = ac_tag_table (variant, "mRNA", h1) ;
	  if (mrnaTable)
	    {
	      int da = 0 ;
	      a1 = a2 = 0 ;
	      x1 = ac_table_int (mrnaTable, 0, 1, 0) ;
	      x2 = ac_table_int (mrnaTable, 0, 2, 0) ;

	      if (ac_has_tag (variant, "Substitution"))
		{
		}
	      else if (ac_has_tag (variant, "Deletion"))
		{
		  da = ac_tag_int (variant, "Multi_deletion", 1) ;
		  x2 = x1 + da + 1 ;
		}
	      else if (ac_has_tag (variant, "Insertion"))
		{
		  x2 = x1 + 1 ;
		}
	      mm = ac_table_printable (mrnaTable, 0, 0, 0) ;
	      if (mm && dictFind (tsnp->targetDict, mm, &mrna) && tsnpRemap1Do (tsnp, mrna, x1, x2, &chrom, &a1, &a2, &strand))
		{
		  nn2++ ;
		  vtxtPrintf (txt, "Variant %s\n", ac_protect (ac_name (variant), h1)) ;
		  /* fix inconsistent coordinates */
		  vtxtPrintf (txt, "mRNA %s %d %d\n", ac_protect (dictName (tsnp->targetDict, mrna), h1), x1, x2) ;
		  vtxtPrintf (txt, "IntMap %s %d %d\n\n", dictName (tsnp->chromDict, chrom), a1, a2) ; 
		}
	    }
	}
      ac_parse (tsnp->db, vtxtPtr (txt), &errors, 0, h) ; 
    }
  fprintf(stderr, "tsnpRemap1 found %d variants remapped %d\n", nn1, nn2) ;

  if (errors && *errors) fprintf(stderr, "tsnpRemap1 parsing error %s\n", errors) ;
  ac_free (h1) ;
  ac_free (h) ;
  return nn2 ;
} /* tsnpRemap1 */

/*************************************************************************************/
/* 
 * We import all the relevant data from the ZZ database
 *   Remap the genome variants into transcript coordinates
 *   Remap to the genebox without giving precise coordinates
 */

static BOOL tsnpRemap2geneBoxDo (TSNP *tsnp, int chrom, int pos, int *JJp, int *geneBoxp, int *x1p, int *strandp)
{
  int ii, target_class = 0 ;
  RMP *up ;
  Array atlas, map ;
  
  dictAdd (tsnp->target_classDict, "any", &target_class) ;
  if (tsnp->target2geneAtlas && 
      target_class < arrayMax (tsnp->target2geneAtlas) && 
      (atlas =  arr (tsnp->target2geneAtlas, target_class, Array)) &&
      chrom < arrayMax(atlas) &&
      (map = array (atlas, chrom, Array))
      )
    {
      for (ii = *JJp, up = arrp (map, ii, RMP) ; ii < arrayMax (map) ; ii++, up++)
	{
	  if (up->a1 <= pos && up->a2 >= pos)
	    {
	      if (up->x1 < up->x2)
		{ *x1p = up->x1 + pos - up->a1 - PROMOTOR_DELTA ; *strandp = 1 ; *geneBoxp = up->target ; *JJp = ii + 1 ; return TRUE ; }
	      else
		{ *x1p = up->x1 - pos + up->a1  - PROMOTOR_DELTA ; *strandp = - 1 ; *geneBoxp = up->target ; *JJp = ii + 1 ;  return TRUE ; }
	    }
	  else if (pos < up->a1)
	    break ;
	}
    }
  return FALSE ;
} /* tsnpRemap2Do */

static BOOL tsnpRemap2Do (TSNP *tsnp, int chrom, int pos, int *JJp, int *mrnap, int *x1p, int *strandp)
{
  int ii, target_class = 0 ;
  RMP *up ;
  Array atlas, map ;
  
  dictAdd (tsnp->target_classDict, "any", &target_class) ;
  if (tsnp->target2exonAtlas && 
      target_class < arrayMax (tsnp->target2exonAtlas) && 
      (atlas =  arr (tsnp->target2exonAtlas, target_class, Array)) &&
      chrom < arrayMax(atlas) &&
      (map = array (atlas, chrom, Array))
      )
    {
      for (ii = *JJp, up = arrp (map, ii, RMP) ; ii < arrayMax (map) ; ii++, up++)
	{
	  if (up->a1 <= pos && up->a2 >= pos)
	    {
	      if (up->x1 < up->x2)
		{ *x1p = up->x1 + pos - up->a1 ; *strandp = 1 ; *mrnap = up->target ; *JJp = ii + 1 ; return TRUE ; }
	      else
		{ *x1p = up->x1 - pos + up->a1 ; *strandp = - 1 ; *mrnap = up->target ; *JJp = ii + 1 ;  return TRUE ; }
	    }
	  else if (pos < up->a1)
	    break ;
	}
    }
  return FALSE ;
} /* tsnpRemap2Do */

/*************************************************************************************/
/* scan the VariantDB acedb database
 * add the remap info 
 *   Remap the genome variants into transcript coordinates\n"
 */
static int tsnpRemap2 (TSNP *tsnp)
{
  AC_HANDLE  h1 = 0, h = ac_new_handle () ;
  AC_ITER iter ;
  AC_OBJ variant = 0 ;
  AC_TABLE intMapTable = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  int pos, pos2, x1 = 0, strand, chrom, mrna = 0, geneBox = 0, nn1 = 0, nn2 = 0, nn3 = 0, JJ=0 ;
  const char *chromNam ;
  const char *errors = 0 ;
  
  if (tsnp->force)
    iter = ac_query_iter (tsnp->db, TRUE, "find variant found_in_genome ", 0, h) ;
  else
    iter = ac_query_iter (tsnp->db, TRUE, "find variant found_in_genome && ! geneBox", 0, h) ;
  while (ac_free (variant), ac_free (h1), variant = ac_iter_obj (iter))
    {
      h1 = ac_new_handle () ;
      nn1++ ;
      intMapTable = ac_tag_table (variant, "IntMap", h1) ;
      if (intMapTable)
	{
	  chromNam = ac_table_printable (intMapTable, 0, 0, 0) ;
	  pos = ac_table_int (intMapTable, 0, 1, 0) ;
	  pos2 = ac_table_int (intMapTable, 0, 2, 0) ;
	  if (pos2 == 1 || pos2 > pos) pos2 = 1 ;
	  else if (pos2 == -1 || pos2 < pos) pos2 = 1 ;
	  
	  if (pos && dictFind (tsnp->chromDict, chromNam, &chrom))
	    {
	      JJ = 0 ;
	      if (! ac_has_tag (variant, "mRNA") && ! ac_has_tag (variant, "Found_in_mRNA"))
		while (tsnpRemap2Do (tsnp, chrom, pos, &JJ, &mrna, &x1, &strand))
		  {
		    nn2++ ;
		    vtxtPrintf (txt, "Variant %s\n", ac_protect (ac_name (variant), h1)) ;
		    vtxtPrintf (txt, "mRNA %s %d %d\n\n", dictName (tsnp->targetDict, mrna), x1, strand * pos2) ; 
		  }
	      JJ = 0 ;
	      while (tsnpRemap2geneBoxDo (tsnp, chrom, pos, &JJ, &geneBox, &x1, &strand))
		{
		  nn2++ ;
		  vtxtPrintf (txt, "Variant %s\n", ac_protect (ac_name (variant), h1)) ;
		  vtxtPrintf (txt, "GeneBox %s %d %d\n\n", dictName (tsnp->targetDict, geneBox), x1, strand * pos2) ; 
		}
	    }
	}
    }
  fprintf(stderr, "tsnpRemap2 found %d variants remapped %d, located %d genes\n", nn1, nn2, nn3) ;
  if (vtxtPtr (txt))
    {
      ACEOUT ao = aceOutCreate ("toto", 0, 0, h) ;
      aceOut (ao,  vtxtPtr (txt)) ;
      ac_free (ao) ;
    }
  ac_parse (tsnp->db, vtxtPtr (txt), &errors, 0, h) ; 
  if (*errors) fprintf(stderr, "tsnpRemap parsing error %s\n", errors) ;
  ac_free (h1) ;
  ac_free (h) ;
  return nn2 ;
} /* tsnpRemap2 */

/*************************************************************************************/
/*************************************************************************************/
/* import the runs metadata from the database */

static int tsnpGetRunList (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *errors = 0 ;
  AC_TABLE tbl = tsnp->runMetaData = ac_bql_table (tsnp->db, hprintf (h, "select r, t, st from p in class project where p == \"%s\", r in p->run where r ISA runs, t in r->title, st in r->Sorting_title", tsnp->project), 0, "st", &errors, tsnp->h) ;
  int ir ;
  const char *ccp ;

  if (! tbl)
    messcrash ("No run belong to project:%s in database \n", tsnp->project ? tsnp->project : "unspecified") ;
  tsnp->runDict = dictHandleCreate (tbl->rows + 10, tsnp->h) ;

  /* ATTENTION: dans la table BQL il ne faut avoir que des data uniques associees au run pour avoir toujours exactement 1 ligne par run, sinon on a une erreur d'attribution des nombres qui eux utilisent runDict pour la memorisation */
  for (ir = 0 ; ir < tbl->rows ; ir++)
    {
      ccp = ac_table_printable (tbl, ir, 0, 0) ;
      if (! ccp || ! dictAdd (tsnp->runDict, ccp, 0))
	messcrash ("Asynchrony it is necessary to ensure that the dict and the table agree on the numbering of the runs ") ;
      if (strcmp (ccp, dictName(tsnp->runDict, ir+1)))
	messcrash ("Asynchrony it is necessary to ensure that the dict and the table agree on the numbering of the runs ") ;
    }

  ac_free (h) ;
  return dictMax (tsnp->runDict) ;
} /* tsnpGetRunList */

/*************************************************************************************/

static void tsnpReportFrequency (ACEOUT ao, ACEOUT ao2, TSNP *tsnp, int line, SNP *up, AC_KEYSET ks, BOOL intertwinFrequencyCounts)
{
  DICT *runDict = tsnp->runDict ;
  int ir, irMax = runDict ? dictMax (runDict) : 0 ;
  AC_TABLE tbl = tsnp->runMetaData ;
  BOOL hasW = tsnp->wiggleDir ? TRUE : FALSE ;
  
  if (ao)
    aceOutf (ao, "\t") ;
  if (line) /* titles */
    {
      if (ao)
	{
	  aceOutf (ao, "\tVariant") ;
	  for (ir = 1 ; ir <= irMax ;ir++)
	    {
	      aceOutf (ao, "\t%s", ac_table_printable (tbl, ir - 1, line-1, "")) ;
	      if (intertwinFrequencyCounts)
		aceOutf (ao, "\t%s counting v:variant, r:reference, w:wiggle at 5' and 3' sites if different", ac_table_printable (tbl, ir - 1, line-1, "")) ;
	    }
	}
    }
  else
    {
      if (ao)
	aceOutf (ao, "\t%s", dictName (tsnp->snpDict, up->snp)) ;
      for (ir = 1 ; ir <= irMax ;ir++)
	{
	  int m5 = 0, m3 = 0, r5 = 0, r3 = 3, w5 = 0, w3 = 0, z3, z5 ;
	  char wiggleProblem = ' ' ;
	  float f5, f3, f ;
	  f3 = f5 = f = -20 ;
	  if (up->counts)
	    {
	      KEY *kp = arrayp (up->counts, MAXTYPE*ir, KEY) ;
	      m5 = kp[0] ;
	      m3 = kp[1] ;
	      r5 = kp[2] ;
	      r3 = kp[3] ;
	      w5 = kp[4] ; if (w5 < 0) {  wiggleProblem = '#' ; w5 = -w5 ; }
	      w3 = kp[5] ; 
	      if (w3 < 0)  /* wiggle too low */
		{ 
		  wiggleProblem = '@' ; 
		  z3 = z5 = w5  ; /* true neighbour coverage */
		  w3 = w5 ;
		} 
	      else
		{
		  z3 = m3 + r3 ; if (w3 > 2 * z3) z3 = w3 ;
		  z5 = m5 + r5 ; if (w5 > 2 * z5) z5 = w5 ;
		}
	      f5 = (z5 >= 20 ?  100.0 * m5 /z5 : -20) ;
	      f3 = (z3 >= 20 ?  100.0 * m3 /z3 : -20) ;
	      f  = (z3 + z5 >= 40 ?  100.0 * (m3 + m5) /(z3 + z5) : -20) ;
	      if (ao2)
		{ /* tsf export of same data */
		  aceOutf (ao2, "%s\t%s" , dictName (tsnp->snpDict, up->snp), dictName(tsnp->runDict, ir)) ;
		  aceOutf (ao2, "\t10i\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n", m5,m3,r5,r3,w5,w3, z3,z5,f3, f5) ;
		}
	    }
	  if (ao)
	    {
	      aceOut (ao, "\t") ;
	      if (0 && (f3  > f5+10 || f3 < f5 - 10))
		{ aceOutPercent (ao, f5) ; aceOut (ao, " ") ; aceOutPercent (ao, f3) ; }
	      else
		aceOutPercent (ao, f) ;
	      if (intertwinFrequencyCounts)
		{ 
		  { /* ATTENTION duplicated in ReportCounts */
		    if (tsnp->doubleReport)
		      {
			if (up->type & 0x1)
			  aceOutf (ao, "\t%cv:%d", wiggleProblem, m5) ;
			else
			  aceOutf (ao, "\t%cv:%d_%d", wiggleProblem, m5,m3) ;
			if (up->type & 0x2)
			  aceOutf (ao, " r:%d", r5) ;
			else
			  aceOutf (ao, " r:%d_%d", r5,r3) ;
			if (hasW)
			  {
			    if ((up->type & 0x3)  && w3 != w5)
			      aceOutf (ao, " w:%d_%d", w5, w3) ;
			    else
			      aceOutf (ao, " w:%d", (w3+w5)/2) ;
			  }
		      }
		    else
		      {
			aceOutf (ao, "\t%cv:%d r:%d", wiggleProblem, (m5 + m3), (r5 + r3)/2) ;
			if (hasW)
			  aceOutf (ao, " w:%d", (w5 + w3)/2) ;
		      }
		    
		  }
		}
	    }
	}
    }
} /* tsnpReportFrequency */

/*************************************************************************************/

static void tsnpReportCounts (ACEOUT ao, TSNP *tsnp, int line, SNP *up, AC_KEYSET ks, BOOL best)
{
  DICT *runDict = tsnp->runDict ;
  int ir, irMax = dictMax (runDict) ;
  AC_TABLE tbl = tsnp->runMetaData ;
  BOOL hasW = tsnp->wiggleDir ? TRUE : FALSE ;

  if (1 || ! best) aceOutf (ao, "\t") ;
  if (line) /* titles */
    {
      if (best)
	{
	  aceOut (ao, "\t% variant reads 100*v/(v+ref)");
	  aceOut (ao, "\tLowest frequency\tHighest frequency") ;
	  aceOutf (ao, "\tTop run\t%s", "Top run counts v:variant, r:reference, w:wiggle at 5' and 3' sites if different") ;
	}
      else
	for (ir = 1 ; ir <= irMax ;ir++)
	  aceOutf (ao, "\t%s counting v:variant, r:reference, w:wiggle at 5' and 3' sites if different", ac_table_printable (tbl, ir - 1, line-1, "")) ;
    }
  else
    {
      int topRun = up->topRun ;
      if (best)
	{
	  int w, m ;
	  float f = -20 ;
	  
	  if (0)
	    {
	      if (up->awm > 5 * up->wm)  up->wm = up->awm ;
	      if (up->wm > 5 * up->awm)  up->awm = up->wm ;
	      if (up->awp > 5 * up->wp)  up->wp = up->awp ;
	      if (up->wp > 5 * up->awp)  up->awp = up->wp ;
	    }
	  w = up->wp + up->wm ; m = up->mp + up->mm ; f = (w + m >= 20 ? 100.0 * m / ((float)(m + w)) : -20) ;
	  aceOut (ao, "\t") ;      aceOutPercent (ao, f) ;

	  f = up->minFrequency < 200 ? up->minFrequency : -tsnp->minSnpCover ; aceOut (ao, "\t") ;      aceOutPercent (ao, f) ;
	  f = up->maxFrequency > -tsnp->minSnpCover  ? up->maxFrequency : -tsnp->minSnpCover ; aceOut (ao, "\t") ;      aceOutPercent (ao, f) ;
	}

      if (best && !topRun)
	aceOutf (ao, "\t\t") ;
      else
	for (ir = 1 ; ir <= irMax ;ir++)
	  { 
	    int m5 = 0, m3 = 0, r5 = 0, r3 = 0, w5 = 0, w3 = 0, z3, z5 ;
	    char wiggleProblem = ' ' ;

	    if (best)
	      {
		if (ir != topRun)
		  continue ;
		aceOutf (ao, "\t%s", dictName(tsnp->runDict, topRun)) ;			
	      }
	    if (up->counts)
	      {
		KEY *kp = arrayp (up->counts, MAXTYPE*ir, KEY) ;
		m5 = kp[0] ;
		m3 = kp[1] ;
		r5 = kp[2] ;
		r3 = kp[3] ;
		w5 = kp[4] ; if (w5 < 0) {  wiggleProblem = '#' ; w5 = -w5 ; }
		w3 = kp[5] ; 
		if (w3 < 0)  /* wiggle too low */
		  { 
		    wiggleProblem = '@' ; 
		    z3 = z5 = w5  ; /* true neighbour coverage */
		    w3 = w5 ;
		  } 
		else
		  {
		    z3 = m3 + r3 ; if (w3 > 2 * z3) z3 = w3 ;
		    z5 = m5 + r5 ; if (w5 > 2 * z5) z5 = w5 ;
		  }
	      }
	    
	    /* ATTENTION duplicated in ReportFrequency */
	    if (tsnp->doubleReport)
		{
		  if (up->type & 0x1)
		    aceOutf (ao, "\t%cv:%d", wiggleProblem, m5) ;
		  else
		    aceOutf (ao, "\t%cv:%d_%d", wiggleProblem, m5,m3) ;
		  if (up->type & 0x2)
		    aceOutf (ao, " r:%d", r5) ;
		  else
		    aceOutf (ao, " r:%d_%d", r5,r3) ;
		  if (hasW)
		  {
		    if ((up->type & 0x3) && w3 != w5)
		      aceOutf (ao, " w:%d_%d", w5, w3) ;
		    else
		      aceOutf (ao, " w:%d", (w3+w5)/2) ;
		  }
		}
	      else
		{
		  aceOutf (ao, "\t%cv:%d r:%d", wiggleProblem, (m5 + m3), (r5 + r3)/2) ;
		  if (hasW)
		    aceOutf (ao, " w:%d", (w5 + w3)/2) ;
		}
	  }
    }
} /* tsnpReportCounts */

/*************************************************************************************/

static void tsnpGetGlobalCounts (TSNP *tsnp)
{
  DICT *runDict = tsnp->runDict ;
  int ir, irMax = dictMax (runDict) ;
  int snp, snpMax = arrayMax (tsnp->snps) ;
  SNP *up ;
  int minSnpCount = tsnp->minSnpCount ;
  int minSnpCover = tsnp->minSnpCover ;

  if (snpMax > 1)
    for (snp = 1, up = arrp (tsnp->snps, snp, SNP) ; snp < snpMax ; snp++, up++)
      {
	float alleleFrequency = 0 ;
	int nAlleleFrequency = 0 ;

	up->wiggleProblemHigh = 0 ;
	up->wiggleProblemLow = 0 ;
	for (ir = 1 ; ir <= irMax ;ir++)
	  {
	    float f = -minSnpCover ;
	    if (up->counts)
	      {
		float m5, m3, r5, r3, w5, w3, f3, f5, z3, z5 ;
		KEY *kp = arrayp (up->counts, MAXTYPE*ir, KEY) ;
		
		m5 = kp[0] ;
		m3 = kp[1] ;
		r5 = kp[2] ;
		r3 = kp[3] ;
		w5 = kp[4] ; if (w5 < 0) { w5 = -w5 ; }
		w3 = kp[5] ; 
		if (w3 < 0)  /* wiggle too low */
		  { 
		    z3 = z5 = w5  ; /* true neighbour coverage */
		    w3 = w5 ;
		  } 
		else
		  {
		    z3 = m3 + r3 ; if (w3 > 2 * z3) z3 = w3 ;
		    z5 = m5 + r5 ; if (w5 > 2 * z5) z5 = w5 ;
		  }
		if (0 && w3 + w5 == 0)
		  {
		    if (r3 > 5 * r5) r5 = kp[2] = r3 ;
		    if (r5 > 5 * r3) r3 = kp[3] = r5 ;
		  }

		if (m3 + m5 + r3 + r5 >= 2 * minSnpCover)
		  {
		    up->wellCovered++ ; 
		    f = f3 = f5 = -1000 ;
		    z3 = m3 + r3 ; if (w3 > 2 * z3) z3 = w3 ;
		    z5 = m5 + r5 ; if (w5 > 2 * z5) z5 = w5 ;
		    if ((m5 >= minSnpCount || m5 == 0) && m5 + r5 >= 20) f5 = 100.0 * m5 /z5 ;
		    if ((m3 >= minSnpCount || m3 == 0) && m3 + r3 >= 20) f3 = 100.0 * m3 /z3 ;
		    
		    /* take the max */
		    f = (m3 + m5 == 0 || m3 + m5 >= minSnpCount ? 100.0 * (m5 + m3) /(z3 + z5) : f) ;
		    if (f < f5) f = f5 ;
		    if (f < f3) f = f3 ;
		    if (f >= 95)        up->pure++ ; 
		    else if (f >= 80)   up->high++ ; 
		    else if (f >= 20)   up->mid++  ; 
		    else if (f >= 5)    up->low++  ; 
		    else                up->ref++  ; 
		    
		    if (f >= 0) { alleleFrequency += 2 * f ; nAlleleFrequency++ ; }
		    
		    if (up->maxFrequency < f)
		      { up->maxFrequency = f ; up->topRun = ir ; }
		    if (f >= 0 && up->minFrequency > f)
		      up->minFrequency = f ; 
		    if (minSnpCount  && up->maxFrequency < -minSnpCount && m3 + m5 > 0 && m3 + m5 < minSnpCount) 
		      { up->maxFrequency = -minSnpCount ;  up->topRun = ir ; }
		  }
	      }
	  }
	if (nAlleleFrequency)
	  up->alleleFrequency = alleleFrequency / nAlleleFrequency ;
      }
} /* tsnpGetGlobalCounts */

/*************************************************************************************/

static void tsnpReportGlobalCounts (ACEOUT ao, TSNP *tsnp, int line, SNP *up, AC_KEYSET ks)
{
  aceOutf (ao, "\t") ;
  if (line) /* titles */
    {
      aceOut (ao, "\tReads supporting the variant on strand plus\tReads supporting the variant on strand minus") ;
      aceOut (ao, "\tReads supporting the reference  on strand plus\tReads supporting the reference on strand minus") ;


      aceOut (ao, "\tCoverage on plus strand of template\tCoverage on minus strand of template") ;
    }
  else
    {
      if (up->type & 0x1) aceOutf (ao, "\t%d", up->mp) ; else aceOutf (ao, "\t%d_%d", up->mp,up->amp) ;
      if (up->type & 0x1) aceOutf (ao, "\t%d", up->mm,up->amm) ;  else aceOutf (ao, "\t%d_%d", up->mm,up->amm) ;
      if (up->type & 0x2) aceOutf (ao, "\t%d", up->wp) ; else aceOutf (ao, "\t%d_%d", up->wp,up->awp) ;
      if (up->type & 0x2) aceOutf (ao, "\t%d", up->wm) ; else aceOutf (ao, "\t%d_%d", up->wm,up->awm) ;

      if ((up->type & 0x3) == 0x3) aceOutf (ao, "\t%d", up->coverp) ; else      aceOutf (ao, "\t%d_%d", up->coverp, up->acoverp) ;
      if ((up->type & 0x3)== 0x3) aceOutf (ao, "\t%d", up->coverm) ; else      aceOutf (ao, "\t%d_%d", up->coverm, up->acoverm) ;
    }
} /* tsnpReportGlobalCounts */

/*************************************************************************************/

static void tsnpReportPrevalence (ACEOUT ao, TSNP *tsnp, int line, SNP *up, AC_KEYSET ks)
{
  aceOutf (ao, "\t") ;
  if (line == 1) /* titles */
    {
      aceOutf (ao, "\tRuns not measurable (not covered at least 20 times)\tMeasurable (covered at least 20 times)\tReference (<5%%)\tLow (5%%-20%%)\tMid (20%%-80%%)\tHigh (80%%-95%%)\tPure (over 95%%)") ;
      aceOut (ao, "\tAllele frequency in cohort (assuming diploidy)");
    }
  else if (line)
    aceOut (ao, "\t\t\t\t\t\t\t\t") ;
  else
    {
      int runMax = dictMax (tsnp->runDict) ;
      
      aceOutf (ao, "\t%d", runMax - up->wellCovered) ;
      aceOutf (ao, "\t%d", up->wellCovered) ;
      aceOutf (ao, "\t%d", up->ref) ;
      aceOutf (ao, "\t%d", up->low) ;
      aceOutf (ao, "\t%d", up->mid) ;
      aceOutf (ao, "\t%d", up->high) ;
      aceOutf (ao, "\t%d", up->pure) ;
      aceOutf (ao, "\t%.2f", up->alleleFrequency) ;
    }
} /* tsnpReportPrevalence */

/*************************************************************************************/

static void tsnpReportWiggleProblem (ACEOUT ao, TSNP *tsnp, int line, SNP *up, AC_KEYSET ks)
{
  BOOL hasW = tsnp->wiggleDir ? TRUE : FALSE ;

  if (line == 1) /* titles */
    {
      if (hasW)
	aceOutf (ao, "\tWiggle too high\tWiggle too low") ;
      else
	aceOutf (ao, "\t\t") ;
    }
  else if (line)
    aceOut (ao, "\t\t") ;
  else
    {
      if (hasW) 
	aceOutf (ao, "\t%d\t%d", up->wiggleProblemHigh, up->wiggleProblemLow) ;
      else
	aceOut (ao, "\t\t") ;
    }
} /* tsnpReportWiggleProblem */

/*************************************************************************************/

static void tsnpReportCoding (ACEOUT ao, TSNP *tsnp, int line, SNP *up, AC_KEYSET ks)
{
  aceOutf (ao, "\t") ;
  if (line == 1) /* titles */
    {
      aceOutf (ao, "\tHGVS_name but shifted 5' to be compatible with VCF") ;
    }
  else if (line)
    aceOut (ao, "\t") ;
  else 
    {
      AC_HANDLE h = ac_new_handle () ;
      const char *errors = 0 ;
      AC_TABLE tbl = ac_bql_table (tsnp->db, "select s, map, gNam, rNam, pNam, pNam2, tit from s in @, map in s->IntMap, gNam in s->gName, rNam in s->rName, pNam in s->pName, pNam2 in s->Observed__protein_sequence[3], tit in s->title", ks, 0, &errors, h) ;      
      if (tbl)
	{
	  const char *map = ac_table_printable (tbl, 0,1, "") ;
	  const char *gNam = ac_table_printable (tbl, 0,2, "") ;
	  const char *rNam = ac_table_printable (tbl, 0,3, 0) ;
	  const char *pNam = ac_table_printable (tbl, 0,4, "") ;
	  const char *pNam2 = ac_table_printable (tbl, 0,5, "") ;
	  const char *title = ac_table_printable (tbl, 0,6, 0) ;
	  
	  if (title)
	    aceOutf (ao, "\t%s", title) ;
	  else
	    {
	      aceOutf (ao, "\t%s:", map) ;
	      aceOutf (ao, "%s", gNam) ;
	      if (rNam) aceOutf (ao, ", %s", rNam) ;
	      if (pNam[0] || pNam2[0]) aceOutf (ao, " (%s)[%s]", pNam,pNam2) ;
	    }
	}
      else
	aceOut (ao, "\t") ;

      ac_free (h) ;
    }
} /* tsnpReportCoding */

/*************************************************************************************/

static void tsnpReportSnippet (ACEOUT ao, TSNP *tsnp, int line, SNP *up, AC_KEYSET ks)
{
  aceOutf (ao, "\t") ;
  if (line == 1) /* titles */
    {
      aceOutf (ao, "\tGenomic 51mer reference > variant") ;
      aceOutf (ao, "\tProtein 17aa reference > variant") ;
    }
  else if (line)
    aceOut (ao, "\t\t") ;
  else 
    {
      AC_HANDLE h = ac_new_handle () ;
      const char *errors = 0 ;
       AC_TABLE tbl = ac_bql_table (tsnp->db, "select s, pref, pvar, gref, gvar from s in @, pref in s->Reference_protein_sequence[2], pvar in s->Observed__protein_sequence[2], gref in s->Reference_genomic_sequence, gvar in s->Observed__genomic_sequence", ks, 0, &errors, h) ;      
      if (tbl)
	{
	  const char *pR, *pV, *gR, *gV ;
	  
	  pR = ac_table_printable (tbl, 0,1,"") ;
	  aceOutf (ao, "\t%s", pR) ;
	  pV = ac_table_printable (tbl, 0,2,0) ;
	  if (pR && pV && !strcmp (pR,pV))
	    aceOutf (ao, " =") ;
	  else if (pV)
	    aceOutf (ao, " > %s", pV) ;
	  gR = ac_table_printable (tbl, 0,3,"") ;
	  aceOutf (ao, "\t%s", gR) ;
	  gV = ac_table_printable (tbl, 0,4,0) ;
	  if (gV) aceOutf (ao, " > %s", gV) ;
	}
      else
	aceOut (ao, "\t\t") ;

      ac_free (h) ;
    }
} /* tsnpReportSnippet */

/*************************************************************************************/

static void tsnpReportMethod (ACEOUT ao, TSNP *tsnp, int line, SNP *up, AC_KEYSET ks)
{
  if (line == 1) /* titles */
    {
      aceOutf (ao, "\tMethod") ;
    }
  else if (line)
    aceOut (ao, "\t") ;
  else 
    {
      AC_HANDLE h = ac_new_handle () ;
      const char *errors = 0 ;
      AC_TABLE tbl = ac_bql_table (tsnp->db, "select s, m from s in @, m in s=>method", ks, 0, &errors, h) ;      
      if (tbl)
	{
	  aceOutf (ao, "\t%s", ac_table_printable (tbl, 0,1,"")) ;
	}
      else
	aceOut (ao, "\t") ;

      ac_free (h) ;
    }
} /* tsnpReportMethod */

/*************************************************************************************/

static void tsnpReportVCF (ACEOUT ao, TSNP *tsnp, int line, SNP *up, AC_KEYSET ks)
{
  aceOutf (ao, "\t") ;
  if (line == 1) /* titles */
    {
      aceOutf (ao, "\tto\tTarget\tVCF coordinate: substituted base or base before indel\tReference\tVariant\tModified bases") ;
    }
  else if (line)
    aceOut (ao, "\t\t\t\t\t\t") ;
  else if (up)
    {
      AC_HANDLE h = ac_new_handle () ;
      const char *errors = 0 ;

      AC_TABLE tbl = ac_bql_table (tsnp->db, "select s, map, a1, a2, w1, w2, typ, ms, md, mi, from s in @, map in s->IntMap, v in s#VCF, a1 in map[1], a2 in map[2], w1 in v[2], w2 in v[3], typ in s->typ, ms in s->Multi_substitution, md in s->Multi_deletion, mi in s->Multi_insertion", ks, 0, &errors, h) ;      
      if (tbl)
	{
	  const char *ccp, *typ ;
	  int a1 = ac_table_int (tbl, 0, 2, 0) ;
	  int ms = ac_table_int (tbl, 0, 7, 0) ;
	  int a2 = ac_table_int (tbl, 0, 3, 0) ;
	  /*
	  int md = ac_table_int (tbl, 0, 8, 0) ;
	  int mi = ac_table_int (tbl, 0, 9, 0) ;
	  */
	  aceOutf (ao, "\t%d\t%s", a2, ac_table_printable (tbl, 0,1,"")) ;
	  aceOutf (ao, "\t%d", a1);

	  ccp = ac_table_printable (tbl, 0,4,"") ;  aceOutf (ao, "\t%s", ccp) ;
	  ccp = ac_table_printable (tbl, 0,5,"") ;  aceOutf (ao, "\t%s", ccp) ;
	  typ =  ac_table_printable (tbl, 0, 6, "")   ; /* type */
	  if (typ)
	    {
	      const char *ccp = strchr (typ, '>') ;
	      int nn = ccp ? ccp - typ : 0 ;
	      if (! nn && ms && typ[ms] == '2') nn = ms ;	      
	      if (! nn)
		aceOutf (ao, "\t%s", typ) ;
	      else
		{
		  int i ;
		  aceOutf (ao, "\t%c%d%c", ace_upper(typ[0]), a1, ace_upper(typ[nn+1])) ;
		  for (i = 1 ; i < nn ; i++)
		    if (ace_upper(typ[i]) !=  ace_upper(typ[nn+1+i]))
		      aceOutf (ao, ",%c%d%c", ace_upper(typ[i]), a1+i, ace_upper(typ[nn+1+i])) ;
		}
	    }
	}
      else
	aceOut (ao, "\t\t\t\t\t\t") ;

      ac_free (h) ;
    }
} /* tsnpReportVCF */

/*************************************************************************************/

static void tsnpReportLine (ACEOUT ao, ACEOUT ao2, TSNP *tsnp, int line, SNP *up)
{
  AC_HANDLE h = 0 ;
  AC_KEYSET ks = 0 ;
  BOOL intertwinFrequencyCounts = TRUE ;

  if (up && tsnp->db)
    {
      h = ac_new_handle () ;
      {
	const char *ccp = dictName (tsnp->snpDict, up->snp) ;
	if (tsnp->db)
	  ks = ac_dbquery_keyset (tsnp->db, hprintf (h, "find Variant IS \"%s\"", ccp), h) ;
      }
    }
 
  if (ao)
    {
      if (1) tsnpReportMethod (ao, tsnp, line, up, ks) ;
      if (1) tsnpReportWiggleProblem (ao, tsnp, line, up, ks) ;
      if (1) tsnpReportVCF (ao, tsnp, line, up, ks) ;
      if (1) tsnpReportCoding (ao, tsnp, line, up, ks) ;
      if (1) tsnpReportSnippet (ao, tsnp, line, up, ks) ;
      if (1) tsnpReportPrevalence (ao, tsnp, line, up, ks) ;
      if (1) tsnpReportCounts (ao, tsnp, line, up, ks, TRUE) ;
      if (1) tsnpReportGlobalCounts (ao, tsnp, line, up, ks) ;
    }
  if (1) tsnpReportFrequency (ao, ao2, tsnp, line, up, ks, intertwinFrequencyCounts);
  if (ao)
    {
      if (! intertwinFrequencyCounts) tsnpReportCounts (ao, tsnp, line, up, ks, FALSE) ;
      aceOut (ao, "\tZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ\n") ;
    }
  ac_free (h) ;
} /* tsnpReportLine */

/*************************************************************************************/
/* the rejected sites count as rejected in the run only if they are measurable
 * and never count as non measurable
 */
static void tsnpCountRejectedSnps (TSNP *tsnp, SNP *up)
{
  KEYSET counts = up->counts ;
  int ir, irMax = dictMax (tsnp->runDict) ;
  int minSnpCover = tsnp->minSnpCover ;
  if (! counts)
    return ;

  for (ir = 1 ; ir <= irMax ;ir++)
    {
      RC *rc = arrayp (tsnp->runs, ir, RC) ;
      float m5, m3, r5, r3 ;
      KEY *kp = arrayp (up->counts, MAXTYPE*ir, KEY) ;
      
      m5 = kp[0] ;
      m3 = kp[1] ;
      r5 = kp[2] ;
      r3 = kp[3] ;
     
      if (m3 + m5 + r3 + r5 >= 2 * minSnpCover)
	{
	  int k = up->varType ; 
	  if (k < 0 || k > dictMax(tsnp->varTypeDict))
	    k = 0 ;
	  rc->rejected++ ;
	  if (! rc->typesR)
	    rc->typesR = halloc (dictMax(tsnp->varTypeDict)+4 * sizeof(int), tsnp->h) ;
	  rc->typesR[k]++ ;
	}
    }
} /* tsnpCountRejected */

/*************************************************************************************/
/* the rejected sites are not rejected !
 * ditinguish well covered, ... up to pure
 */
static void tsnpCountAcceptedSnps (TSNP *tsnp, SNP *up)
{
  KEYSET counts = up->counts ;
  int ir, irMax = dictMax (tsnp->runDict) ;
  int minSnpCover = tsnp->minSnpCover ;
  int minSnpCount = tsnp->minSnpCount ;
  if (! counts)
    return ;
  for (ir = 1 ; ir <= irMax ;ir++)
    {
      RC *rc = arrayp (tsnp->runs, ir, RC) ;
      float m5, m3, r5, r3, f, f3, f5, w3, w5, z3, z5 ;
      KEY *kp = arrayp (up->counts, MAXTYPE*ir, KEY) ;
      int k ;

      if (! rc->types)
	rc->types = halloc (dictMax(tsnp->varTypeDict)+4 * sizeof(int), tsnp->h) ;
      if (! rc->typesC)
	rc->typesC = halloc (dictMax(tsnp->varTypeDict)+4 * sizeof(int), tsnp->h) ;
      if (! rc->typesR)
	rc->typesR = halloc (dictMax(tsnp->varTypeDict)+4 * sizeof(int), tsnp->h) ;

      m5 = kp[0] ;
      m3 = kp[1] ;
      r5 = kp[2] ;
      r3 = kp[3] ;
      w5 = kp[4] ;
      w3 = kp[5] ;
     
      rc->tested++ ;
      if (m3 + m5 + r3 + r3 < 2 * minSnpCover)
	{ rc->notMeasurable++ ; continue ; }
      rc->wellCovered++ ;

      f = f3 = f5 = -1000 ;
      z3 = m3 + r3 ; if (w3 > 2 * z3) z3 = w3 ;
      z5 = m5 + r5 ; if (w5 > 2 * z5) z5 = w5 ;
      if ((m5 >= minSnpCount || m5 == 0) && m5 + r5 >= 20) f5 = 100.0 * m5 /z5 ;
      if ((m3 >= minSnpCount || m3 == 0) && m3 + r3 >= 20) f3 = 100.0 * m3 /z3 ;
      
      /* take the max */
        f = (m3 + m5 == 0 || m3 + m5 >= 2 * minSnpCount ? 100.0 * (m5 + m3) /(z3 + z5) : f) ;
      if (f < f5) f = f5 ;
      if (f < f3) f = f3 ;
      if (f >= 95)        rc->pure++ ; 
      else if (f >= 80)   rc->high++ ;
      else if (f >= 20)   rc->mid++ ; 
      else if (f >= 5)    rc->low++ ;
      else                rc->ref++ ;
		    
      k = up->varType ;
      if (k < 0 || k > dictMax(tsnp->varTypeDict))
	k = 0 ;
      rc->types[k] ++ ;
      rc->typesC[k] += (m3 + m5)/2 ; 
    }
} /* tsnpCountAccepted */

/*************************************************************************************/

static int tsnpReport (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii, nn = 0, snpMax  = arrayMax (tsnp->snps) ;
  DICT *snpDict = tsnp->snpDict ;

  ACEOUT ao = aceOutCreate (tsnp->outFileName, ".snp_frequency_table.txt", tsnp->gzo, h) ;
  ACEOUT ao2 = tsnp->outFileName ? aceOutCreate (tsnp->outFileName, ".snp_counts.tsf", tsnp->gzo, h) : 0 ;
  aceOutDate (ao, "####", "SNPs and variations table") ;
  if (ao2)
    {
      aceOutDate (ao2, "####", "SNPs and variations count variants, reference and wiggle cover at donor and acceptor sites :\n") ;
      aceOutf (ao2, "# SNP\tRun\t10i\tVar5\tVar3\tRef5\tRef3\tCover5\tCover3\tWiggle5\tWiggel3\tFreq5\tFreq3\n") ;
    }
  tsnp->doubleReport = TRUE ;

  aceOutf (ao, "# Title") ; tsnpReportLine (ao, ao2, tsnp, 2, 0) ;
  aceOutf (ao, "# Sorting_title") ; tsnpReportLine (ao, ao2, tsnp, 3, 0) ;
  aceOutf (ao, "# Run") ; tsnpReportLine (ao, ao2, tsnp, 1, 0) ;
    
  for (ii = 1 ; ii < snpMax  ; ii++)
    {
      SNP *up = arrayp (tsnp->snps, ii, SNP) ;

      if (tsnp->t1 && (up->a1 > tsnp->t2 || up->a2 < tsnp->t1 || up->a1 < tsnp->t1 - 1000 || up->a2 > tsnp->t2 + 1000)) continue ;

      if (! up->select && tsnp->minSnpFrequency > 0 && (up->maxFrequency < 0 || up->maxFrequency < tsnp->minSnpFrequency))
	continue ;
      if (up->coverp + up->coverm + up->acoverp + up->acoverm < 2 * tsnp->minSnpCover)
	continue ;
      if (up->mp + up->mm + up->amp + up->amm < 2 * tsnp->minSnpCount)
	continue ;
      if (! up->select && tsnp->dropMonomodal && up->pure + up->high + up->mid == 0)
	{
	  tsnpCountRejectedSnps (tsnp, up) ;
	  continue ;
	}
      tsnpCountAcceptedSnps (tsnp, up) ;
      aceOutf (ao, "%s", dictName (snpDict, up->snp)) ;
      tsnpReportLine (ao, ao2, tsnp, 0, up) ;
      nn++ ;
    }
  ac_free (h) ;
  return nn ;
} /* tsnpReport */

/*************************************************************************************/

static int tsnpMergeExport (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii, nn = 0, snpMax  = arrayMax (tsnp->snps) ;

  ACEOUT ao = aceOutCreate (tsnp->outFileName, ".snp_frequency.tsf", tsnp->gzo, h) ;
  aceOutDate (ao, "####", "SNPs counts and frequency table") ;
  aceOutf (ao, "# SNP\tRun\t10i\tVar5\tVar3\tRef5\tRef3\tCover5\tCover3\tWiggle5\tWiggel3\tFreq5\tFreq3\n") ;

  for (ii = 1 ; ii < snpMax  ; ii++)
    {
      SNP *up = arrayp (tsnp->snps, ii, SNP) ;

      if (tsnp->t1 && (up->a1 > tsnp->t2 || up->a2 < tsnp->t1 || up->a1 < tsnp->t1 - 1000 || up->a2 > tsnp->t2 + 1000)) continue ;

      if (! up->select && tsnp->minSnpFrequency > 0 && (up->maxFrequency < 0 || up->maxFrequency < tsnp->minSnpFrequency))
	continue ;
      if (up->coverp + up->coverm + up->acoverp + up->acoverm < 2 * tsnp->minSnpCover)
	continue ;
      if (up->mp + up->mm + up->amp + up->amm < 2 * tsnp->minSnpCount)
	continue ;
      if (! up->select && tsnp->dropMonomodal && up->pure + up->high + up->mid == 0)
	{
	  continue ;
	}
      tsnpReportLine (0, ao, tsnp, 0, up) ;
      nn++ ;
    }
  ac_free (h) ;
  return nn ;
} /* tsnpMergeExport */

/*************************************************************************************/
/*************************************************************************************/
/* export snp stats per run that will go to the Ali object then to qc_ssummary */
static void tsnpReportRuns (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (tsnp->outFileName, "snp_runs.ace", tsnp->gzo, h) ;
  int i, ir, irMax = arrayMax (tsnp->runs) ;
  DICT *varTypeDict = tsnp->varTypeDict ;
  int iMax = dictMax (varTypeDict) ;
  RC *rc ;

  if (irMax)
    for (ir = 1, rc = arrp (tsnp->runs, ir, RC) ; ir < irMax ; rc++, ir++)
      {
	aceOutf (ao, "Ali \"%s\"\n", dictName (tsnp->runDict, ir)) ;
	aceOut (ao, "-D SNP\n") ;
	
	if (rc->tested)
	  {
	    aceOutf (ao, "Tested_sites %d any\n", rc->tested) ;
	    aceOutf (ao, "Rejected_sites %d any\n", rc->rejected) ;
	    aceOutf (ao, "Not_measurable_sites %d any\n", rc->notMeasurable) ;
	  }
	if (rc->wellCovered)
	  {
	    aceOutf (ao, "Measured_sites Genomic %d any %d reference %d low %d mid %d high %d pure\n"
		     , rc->wellCovered
		     , rc->ref
		     , rc->low
		     , rc->mid
		     , rc->high
		     , rc->pure
		     ) ;	     
	    aceOutf (ao, "Support_by_sites Genomic %%.0f any %%.0f reference %%.0f low %%.0f mid %%.0f high %%.0f pure\n"
		     , rc->wellCoveredC
		     , rc->refC
		     , rc->lowC
		     , rc->midC
		     , rc->highC
		     , rc->pureC
		     ) ;	     
	  }
	
	if (rc->xwellCovered)
	  aceOutf (ao, "Measured_sites Exonic %d any %d reference %d low %d mid %d high %d pure\n"
		   , rc->xwellCovered
		   , rc->xref
		   , rc->xlow
		   , rc->xmid
		   , rc->xhigh
		   , rc->xpure
		   ) ;
	if (rc->pwellCovered)
	  {
	    aceOutf (ao, "Measured_sites Protein_changing %d any %d reference %d low %d mid %d high %d pure\n"
		     , rc->pwellCovered
		     , rc->pref
		     , rc->plow
		     , rc->pmid
		     , rc->phigh
		     , rc->ppure
		     ) ;	    
	    aceOutf (ao, "Support_by_sites Protein_changing %%.0f any %%.0f reference %%.0f low %%.0f mid %d high %%.0f pure\n"
		     , rc->pwellCoveredC
		     , rc->prefC
		     , rc->plowC
		     , rc->pmidC
		     , rc->phighC
		     , rc->ppureC
		     ) ;	    
	  }
	for (i = 1 ; i <= iMax ; i++)
	  {
	    int n = rc->types[i] ;
	    int nR = rc->typesR[i] ;
	    float nS = rc->typesC[i] ;
	    if (n)
	      aceOutf (ao, "SNP_profile any %s %d %d %.0f\n", dictName (varTypeDict, i), n, nR, nS) ;
	  }
	for (i = 0 ; i < 1 ; i++)
	  {
	    int n = rc->types[i] ;
	    int nR = rc->typesR[i] ;
	    float nS = rc->typesC[i] ;
	    if (n)
	      aceOutf (ao, "SNP_profile any %s %d %d %.0f\n", "Other", n, nR, nS) ;
	  }

	aceOut (ao, "\n") ;
      }
  ac_free (h) ;
  return ;
} /* tsnpReportRuns */

/*************************************************************************************/
/*************************************************************************************/
static DICT *tsnpMakeVarTypeDict (AC_HANDLE h) ;
static void tsnpDbReport (TSNP *tsnp)
{
  tsnp->varTypeDict = tsnpMakeVarTypeDict (tsnp->h) ;

  tsnpSnpParse (tsnp) ;
  if (tsnp->snps && arrayMax (tsnp->snps))
    {
      tsnpGetGlobalCounts (tsnp) ;
      tsnpWiggleParse (tsnp) ;
      tsnp->snpExported = tsnpReport (tsnp) ;
      tsnpReportRuns (tsnp) ;
    }
} /* tsnpDbReport */
  
/*************************************************************************************/
/*************************************************************************************/
static DICT *tsnpMakeVarTypeDict (AC_HANDLE h) ;
static void tsnpMergeCounts (TSNP *tsnp)
{
  tsnp->runDict = dictHandleCreate (256, tsnp->h) ;
  tsnp->varTypeDict = tsnpMakeVarTypeDict (tsnp->h) ;
  tsnp->minSnpCount = tsnp->minSnpFrequency = 0 ;
  tsnpSnpParse (tsnp) ;
  if (tsnp->snps && arrayMax (tsnp->snps))
    {
      tsnpGetGlobalCounts (tsnp) ;
      tsnpMergeExport (tsnp) ;
    }
} /* tsnpMergeCounts */
  
/*************************************************************************************/
/*************************************************************************************/
typedef struct snpRunCounts { KEY run ; int m, c ;} RMC ;

static int rmcOrder (const void *va, const void *vb)
{
  const RMC *up = va ;
  const RMC *vp = vb ;
  int n = up->run - vp->run ;
  return n ;
} /* rmcOrder */

/*************************************************************************************/

static int tsnpDbRun2groups (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;
  DICT *runDict = dictHandleCreate (200, tsnp->h) ;
  AC_DB db = tsnp->db ;
  Array runs = tsnp->runs ;
  const char *errors = 0 ;
  char *qq = hprintf (h, "Find Project \"%s\" ; >run", tsnp->project) ;
  AC_TABLE tbl = ac_bql_table (db, qq, 0, 0, &errors, h) ;
  int ir, irMax = tbl ? tbl->rows : 0 ;

  tsnp->runStatDict = runDict ;
  for (ir = 0 ; ir < irMax ; ir++) 
    dictAdd (runDict, ac_table_printable (tbl, ir, 0, "NULL"), 0) ;

  if (! tsnp->runs)
    runs = tsnp->runs = arrayHandleCreate (256, RC, tsnp->h) ;
  for (ir = 0 ; ir < irMax ; ir++) 
    {
      AC_HANDLE h1 = ac_new_handle () ;
      AC_OBJ Run = ac_table_obj (tbl, ir, 0, h1) ;
      int run ;
      if (dictFind (runDict, ac_name(Run), &run))
	{
	  AC_TABLE tbl2 = ac_obj_bql_table (Run, ">>Group", 0, &errors, h1) ;
	  int jrMax = tbl2 ? tbl2->rows : 0 ;
      
	  if (jrMax)
	    {
	      RC *rc = arrayp (runs, run, RC) ;
	      KEYSET r2g = rc->r2g ;
	  
	      if (! r2g)
		r2g = rc->r2g = keySetHandleCreate (tsnp->h) ;
	      for (int g = 0, k =0, jr = 0 ; jr < jrMax ; jr++)
		if (dictFind (runDict, ac_table_printable (tbl2, jr, 0, 0), &g))
		  keySet (r2g, k++) = g ;
	    }
	}
      ac_free (h1) ;
    }
  
  ac_free (h) ;
  return dictMax (runDict) ;
} /* tsnpDbRun2groups */

/*************************************************************************************/
/* export a tsf file per run and group counting supported snps */
static void tsnpDbSnpProfile (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_HANDLE h1 = 0 ;
  AC_DB db = tsnp->db ;
  AC_ITER iter ;
  AC_OBJ Snp = 0 ;
  Array rmcs = arrayHandleCreate (300, RMC, h) ;
  Array stats = arrayHandleCreate (1000000, int, h) ;
  DICT *runDict = tsnp->runStatDict ;
  DICT *dict = dictHandleCreate (10000, h) ;
  DICT *gNameDict = dictHandleCreate (10000, h) ;
  int minSnpCover = tsnp->minSnpCover ;
  int minSnpCount = tsnp->minSnpCount ;
  int minSnpFrequency = tsnp->minSnpFrequency ;
  char *allSubs = "A>G_T>C_G>A_C>T_A>T_T>A_G>C_C>G_A>C_T>G_G>T_C>A" ;
  int nSites = 0, nIter = 0 ;
  tsnp->varTypeDict = tsnpMakeVarTypeDict (h) ;
  int nnn = 0 ;
  char typ2[256] ;

  memset (typ2, 0, sizeof (typ2)) ;
  if (tsnp->select)
    iter = ac_query_iter (db, TRUE, hprintf (h, "find variant IS  \"%s\" ", tsnp->select), 0, h) ;
  else
    iter = ac_query_iter (db, TRUE, "Find Variant typ && BRS_counts", 0, h) ;
  
  while (ac_free (Snp), ac_free (h1), Snp = ac_iter_obj (iter))
    {
      AC_TABLE tbl = 0 ;
      int irmc = 0, ir, irMax = dictMax (runDict) ;
      RMC *rmc ;  
      const char *typ = ac_tag_printable (Snp, "Typ", 0) ;
      int i, k, p = 0 ;
      const char *gName = ac_tag_printable (Snp, "gName", 0) ;

      nIter++ ;
      if (! gName || ! dictAdd (gNameDict, gName, 0))
	continue ;
      if (! typ)
	continue ;
      
      nSites++ ;
      h1 = ac_new_handle () ;
      tbl = ac_tag_table (Snp, "BRS_counts", h1) ;

      if (typ[1]== '>')
	{ 
	  char *cp = strstr (allSubs, typ) ;
	  p = 29 ; memcpy (typ2+3, typ, 3) ; typ2[6] = 0 ; 
	  if (cp) { p = (cp - allSubs) ; p = 10 + p/4 ; }
	}
      else
	{
	  k = strlen (typ) ; 
	  if (k > 6) k = 6 ; 
	  k -= 3 ;
	  if (k <= 0)
	    continue ;
	  
	  if (! strncmp (typ, "Del", 3))
	    { 
	      p = 20 + 10*k + 3 ;
	      for (i = 0 ; i < k ; i++)
		typ2[3+i] = '-' ;
	      for (i = 0 ; i < k ; i++)
		typ2[3+k+i] = typ[i+3] ;
	      typ2[3+2*k] = 0 ;
	    }
	  else if (! strncmp (typ, "Ins", 3))
	    { 
	      p = 20 + 10*k + 1 ;
	      for (i = 0 ; i < k ; i++)
		typ2[3+i] = '+' ;
	      for (i = 0 ; i < k ; i++)
		typ2[3+k+i] = typ[i+3] ;
	      typ2[3+2*k] = 0 ;
	    }
	  else if (! strncmp (typ, "Dim", 3))
	    { 
	      p = 20 + 10*k + 4 ;
	      typ2[3] = '*' ;
	      for (i = 0 ; i < k ; i++)
		typ2[3+1+i] = '-' ;
	      for (i = 0 ; i < k ; i++)
		typ2[3+1+k+i] = typ[i+3] ;
	      typ2[3+2*k+1] = 0 ;
	    }
	  else if (! strncmp (typ, "Dup", 3))
	    { 
	      p = 20 + 10*k + 2 ;
	      typ2[3] = '*' ;
	      for (i = 0 ; i < k ; i++)
		typ2[3+1+i] = '+' ;
	      for (i = 0 ; i < k ; i++)
		typ2[3+1+k+i] = typ[i+3] ;
	      typ2[3+2*k+1] = 0 ;
	    }
	}
      typ2[0] = '0' + p/10 ; 
      typ2[1] = '0' + p%10 ; 
      typ2[2] = '_' ;
      for (char *cp = typ2 + 3 ; *cp ; cp++)
	*cp = ace_lower (*cp) ;
      
      for (ir = 0 ; ir < irMax ; ir++)
	{
	  const char *runNam = ac_table_printable (tbl, ir, 0, 0) ;
	  int run = 0 ;
	  if (dictFind (runDict, runNam, &run))
	    {
	      int c = ac_table_int (tbl, ir, 1, 0) ;
	      int m = ac_table_int (tbl, ir, 2, 0) ;
	      RC *rc = arrayp (tsnp->runs, run, RC) ;
	      KEYSET r2g = rc->r2g ;
	      int igMax = r2g ? keySetMax (r2g) : 0 ;
	      rmc = arrayp (rmcs, irmc++, RMC) ;
	      rmc->run = run ;
	      rmc->c = c ;
	      rmc->m = m ;
	      for (int ig = 0 ; ig < igMax ; ig++)
		{
		  rmc = arrayp (rmcs, irmc++, RMC) ;
		  rmc->run = keySet (r2g, ig) ;
		  rmc->c = c ;
		  rmc->m = m ;
		}
	    }
	}	    
      arraySort (rmcs, rmcOrder) ;
      for (int i = 0 ; i < irmc ; i++)
	{
	  int c = 0, m = 0 ;
	  RMC *rmc = arrp (rmcs, i, RMC) ;
	  KEY run = rmc->run ;

	  for (int j = i ; j < irmc && rmc->run == run ; j++, rmc++)
	    { c += rmc->c ; m += rmc->m ; i = j ; }
	  if (c >= minSnpCover && m >= minSnpCount && 100*m >= minSnpFrequency * c)
	    {
	      int *kp, k = 0, t = -1 ;
	      char buf[276] ;

	      snprintf (buf, 275, "%05d_%s", run, "Measurable") ;
	      dictAdd (dict, buf, &k) ;
	      kp = arrayp (stats, k, int) ;
	      (*kp)++ ;

	      if (ac_has_tag (Snp, "monomodal"))
		{
		  snprintf (buf, 275, "%05d_%s", run, "Rejected_sites") ;
		  dictAdd (dict, buf, &k) ;
		  kp = arrayp (stats, k, int) ;
		  (*kp)++ ;
		  continue ;
		}

	      snprintf (buf, 275, "%05d_TYPE__%s", run, typ2) ;
	      dictAdd (dict, buf, &k) ;
	      kp = arrayp (stats, k, int) ;
	      (*kp)++ ;
	      t = -1 ;
	      if (1)
		{
		  /* char *pTypes = {"Pure","High","Mid","Low","Ref"} ; */
		  if (100 * m > 95*c) t = 4 ;
		  else if (100 * m > 80*c) t = 3 ;
		  else if (100 * m > 20*c) t = 2 ;
		  else if (100 * m > 5*c) t = 1 ;
		  else if (100 * m > 95*c) t = 0 ;
		}
	      if (t >= 0)
		{
		  snprintf (buf, 275, "%05d_Genomic__%d", run, t) ;
		  dictAdd (dict, buf, &k) ;
		  kp = arrayp (stats, k, int) ;
		  (*kp)++ ;	     
		} 
	      t = -1 ;
	      if (ac_has_tag (Snp, "Coding") && !ac_has_tag (Snp, "Synonymous"))
		{
		  /* char *pTypes = {"Pure","High","Mid","Low","Ref"} ; */
		  if (100 * m > 95*c) t = 4 ;
		  else if (100 * m > 80*c) t = 3 ;
		  else if (100 * m > 20*c) t = 2 ;
		  else if (100 * m > 5*c) t = 1 ;
		  else if (100 * m > 95*c) t = 0 ;
		}
	      if (t >= 0)
		{
		  snprintf (buf, 275, "%05d_Protein_changing__%d", run, t) ;
		  dictAdd (dict, buf, &k) ;
		  kp = arrayp (stats, k, int) ;
		  (*kp)++ ;	     
		} 
	    }
	}
    }

  /* export in tsf format */
  if (1)
    {
      ACEOUT ao = aceOutCreate (tsnp->outFileName, ".snp_profile.tsf", tsnp->gzo, h) ;
      int run = 0, kMax = arrayMax (stats) ;

      for (run = 1 ; run <= dictMax (runDict) ; run++)
	aceOutf (ao, "%s\t%s\ti\t%d\n", dictName (runDict, run), "Tested_sites", nSites) ; 
      
      for (int k = 1 ; k < kMax ; k++)
	{
	  int nn = array (stats, k, int) ;
	  if (nn)
	    {
	      const char *ccp = dictName (dict, k) ;
	      sscanf (ccp, "%d_", &run) ;
	      while (*ccp != '_') ccp++ ; 
	      ccp++ ;
	      if (!strcmp (ccp, "Measurable"))
		{ 
		  aceOutf (ao, "%s\t%s\ti\t%d\n", dictName (runDict, run), ccp, nn) ; 
		  nn = nSites - nn ; ccp = "Not_measurable_sites" ; 
		}
	      aceOutf (ao, "%s\t%s\ti\t%d\n", dictName (runDict, run), ccp, nn) ;
	      nnn++ ; 
	    }
	}
 
      fprintf (stderr, "tsnpDbSnpProfile exported %d/%d snps in file %s\n"
	       , nnn, nIter, aceOutFileName (ao)
	       ) ;
    }
  
  ac_free (h) ;
  return ;
} /* tsnpDbSnpProfile */

/*************************************************************************************/
/*************************************************************************************/
/* slide deletions and insertions, speclial treat for position 76 */

static BOOL tsnpSlide (const char *dna, int dnaLn, int a10, int da, int *dxp, int *slidep)
{
  int i, dx, dy, a1 = a10, xx ; 
  BOOL isDim = FALSE ;
  
  /* slide left */
  dx = 0 ; i = a1 - 1 ;
  while (i > 0 && dna[i] == dna[i + da])
    { i-- ; dx++ ; }
  a1 -= dx ;
  /* slide right */
  dy = 0 ; i = a1 ;
  while (i + da < dnaLn - 1 && dna[i] == dna[i + da])
    { i++ ; dy++ ; }
  if (dy >= da)
    isDim = TRUE ;
  xx = 75 ;
  if (a1 < xx && a1 + dy >= xx && da > 1000)
    { dy = xx - a1 ; a1 += dy ; }

  *slidep = dy ;
  *dxp = a1 - a10 ;
  return isDim ;
} /* tsnpSlide */

/*************************************************************************************/
/* slide deletions and insertions, speclial treat for position 76 */
/* buf is the most 5' insert in the mRNA 
 * bufG is its reverse complement (most 3' insert in the genome assuming negative strand
 * bufS is the most 3' insert in the mRNA (thus equal to buf rolled to th right on the mRNA) 
 * bufGS is its reverse complement (thus equal to bufG rolled to the right on the genome)
 */

static BOOL tsnpSlideDup (const char *dna, int dnaLn, int a10, int da, int *dxp, int *slidep, const char *insert0, char *buf, char *bufS, char *bufG, char *bufGS)
{
  int i, j, k, dx, dy, a1 = a10 ;
  BOOL isDup = FALSE ;
  char insert[da+1] ;

  if (strlen (insert0) != da)
    for (k = 0 ; k < da ; k++)
      insert[k] = 'n' ;
  else
    for (k = 0 ; k < da ; k++)
      insert[k] = insert0[k] ;
  insert[da] = 0 ;

  /* slide left */
  dx = 0 ; i = a1 - 1 ; j = da - 1 ;
  while (i > 0 && dna[i] == ace_lower(insert [j]))
    { i-- ; dx++ ; j-- ; if (j < 0) j += da ; }
  a1 -= dx ;
  /* copy the rotated buffer */
  for (k = 0 ; k < da ; k++)
    buf[k] = insert[(j+1+k+da) % da] ;
  buf[k] = 0 ;
  /* slide right */
  dy = 0 ; i = a1 ; j = 0 ;
  while (i < dnaLn - 1 && dna[i] == ace_lower(buf[j]))
    { i++ ; dy++ ; j = (j+1) % da ; }
  if (dy >= da)
    isDup = TRUE ;

  for (k = 0 ; k < da ; k++)
    buf[k] = ace_upper (buf[k]) ;
  for (k = 0 ; k < da ; k++)
    bufS[k] = buf[(k + dy) % da] ;
  for (k = 0 ; k < da ; k++)
    {
      bufG[k] = ace_upper (complementLetter (buf[da-k-1])) ;
      bufGS[k] = ace_upper (complementLetter (bufS[da-k-1])) ;
    }
  buf[da] = bufS[da] = bufG[da] = bufGS[da] = 0 ;
  *slidep = dy ;
  *dxp = a1 - a10 ;
  return isDup ;
} /* tsnpSlideDup */

/*************************************************************************************/

static const char *Rdna = 0 ;
static Array Gdna = 0 ;
static int RdnaLn = 0 ;
static int myMrna = 0 ;
static int g2m = 0 ;
static int dnaLn = 0 ;
BOOL fromMrna = FALSE ;

static char myRdna (int a) 
{
  int m = a + (fromMrna ? 0 : g2m) ;
  int cc = 'n' ;
  
  if (m >=0 && m < RdnaLn)
    {
      cc = Rdna[m] ;
      if (cc == 't') cc = 'u' ;
    }
  return cc ;
} /* myRdna */

/*************************************************************************************/

static int myTranslate (char *buf, char tBuf1[], char tBuf3[], int m1, int p1, AC_HANDLE h) 
{
  Array aa = arrayHandleCreate (strlen (buf) + 8, char, h) ;
  char * translationTable = pepGetTranslationTable(myMrna, 0) ; 
  char cc, *cp1, *cp2, *cB1, *cB3 ;
  int frame = 0, j, i = 50 + strlen (buf) ;
  char *tbuf = halloc (i, h) ;
  
  array(aa, i, char) = 0 ;
  arrayMax(aa) = i ;
  cp1 = arrp (aa, 0, char) ; cp2 = buf ;
  for (i = 0 ; *cp2 ; cp2++)
    { cc = dnaEncodeChar [(int)*cp2] ; if (cc) { *cp1++ = cc ; i++ ; }}
  *cp1 = 0 ; 
  if (i >= arrayMax(aa)) messcrash ("n importe quoi dans myTranslate") ;
  arrayMax (aa) = i ;
  memset (tbuf, 0, i) ;
  tbuf[0] = 0 ;
  i = m1 - p1 ;
  i = i % 3 ;
  i = 15 + 2 - i ;
  i = i % 3 ;
  cB1 = tBuf1 ;   cB3 = tBuf3 ; 
  for (j = 0 ; j < i ; j++)
    *cB3++ = '-' ;
  if (arrayMax (aa) > 2)
    for (cp1 = arrp (aa, i, char)  ; i < arrayMax (aa) - 2 ; i+= 3 , cp1 += 3 )
      {
	int j ;
	char cc = e_codon (cp1, translationTable) ;
	const char *pep = pepShortName[(int)cc] ;
	if (i>=19 && i <=21) cc = ace_lower (cc) ;
	*cB1++ = cc ;
	if (pep[0] == 'X') pep = "---" ;
	if (cc == '*') pep = "Ter" ;
	for (j = 0 ; j < 3 ; j++)
	  *cB3++ = *pep++ ;
      }
  *cB1 = *cB3 = 0 ;
  return frame ;
} /* myTranslate */

/*************************************************************************************/
/* dna is a char string, starting at pos, locate the next stop 
 * pos is in bio convention, first base is called 1
 * met is distances relative to pos
 * if met == 0, pos is a met codon
 * if no met is detected, the code returns met == -1
 */
static BOOL locateMet (const char *dna,  const char *translationTable, int pos, int *met)
{
  int i, j, n = strlen (dna) ;
  const char *cp ;
  Array dna2 = arrayHandleCreate (n + 3, char, 0) ;

  array (dna2, n - 1, char) = 0 ;
  memcpy (arrp (dna2, 0, char), dna, n) ;
  dnaEncodeArray (dna2) ;

  /* the new met may me way outside the snippet, we must translate the whole mrna */
  *met = -1 ;  /* stop not found */
  if (pos > 0)
    for (i = pos - 1, j = 0, cp = arrp (dna2, i, char) ; 
	 i < arrayMax (dna2) - 2 ; i+= 3, cp += 3, j++)
      if (e_codon (cp, translationTable) == 'M')
	{ *met = j ; break ; }

  ac_free (dna2) ;
  return TRUE ;
} /* locateMet */

/*************************************************************************************/
/* dna is a char string, starting at pos1 or pos2, locate the next stop 
 * pos1, pos2 are in bio convention, first base is called 1
 * stop1/2 are distances relative to pos1/2
 * if stop1 == 0, pos1 is a stop codon
 * if no stop is detected, the code returns stop1 == -1
 */
static BOOL locateStop (const char *dna,  const char *translationTable, int pos1, int pos2, int *stop1, int *stop2)
{
  int i, j, n = strlen (dna) ;
  const char *cp ;
  Array dna2 = arrayHandleCreate (n + 3, char, 0) ;

  array (dna2, n - 1, char) = 0 ;
  memcpy (arrp (dna2, 0, char), dna, n) ;
  dnaEncodeArray (dna2) ;

  /* the new stop may me way outside the snippet, we must translate the whole mrna */
  *stop1 = *stop2 = -1 ;  /* stop not found */
  if (pos1 > 0)
    for (i = pos1 - 1, j = 0, cp = arrp (dna2, i, char) ; 
	 i < arrayMax (dna2) - 2 ; i+= 3, cp += 3, j++)
      if (e_codon (cp, translationTable) == '*')
	{ *stop1 = j ; break ; }
  if (pos2 > 0)
    for (i = pos2 - 1, j = 0, cp = arrp (dna2, i, char) ; 
	 i < arrayMax (dna2) - 2 ; i+= 3, cp += 3, j++)
      if (e_codon (cp, translationTable) == '*')
	{ *stop2 = j ; break ; }
  
  ac_free (dna2) ;
  return TRUE ;
} /* locateStop */

/*************************************************************************************/
/* if i kill Met1, then use Leu2_Met124del  (i.e. we now use Met14 as our new Met)
 * stops are caller Ter
 * extension: the stop has a substitution :  Ter110GlyextTer31
 * (Arg123LysfsTer34)  paranthesis are not needed if we have experimental RNA or peptides 
 */
static void tsnpSetPName (vTXT txt, KEY product, char *pR1, char *pV1, char *pRef, char *pVar, int m1, int p1, int frame, int fs) 
{
  char *cp, *cq ;
  char * translationTable = pepGetTranslationTable(myMrna, 0) ; 
  char *sep = "" ;
  char buf[6];
  char buf0[6];
  char buf2[6];
  char buf12[6];
  char buf22[6];
  char idem[4] ;
  int i=0, j, k = 0, m ;
  m1 = (m1 - p1 + (fs ? 1 : 0)) / 3 ; m1++ ;
  
  if (1) /* sub del ins */
    {
      if (1) vtxtPrintf (txt, "\npName p.%s:", name(product)) ;
      cp = pRef ; cq = pVar ;
      while (*cp == '-') { cp++; cq++ ; i++; }
      while (*cp)
	{
	  m = m1 - 7+(i+2)/3 ;
	  if (m == m1)
	    {  
	      for (j = 0 ; j < 3 ; j++)
		idem[j] = cp[j] ;
	      idem[3] = 0 ;
	    }
	  if (!strncmp (cp, "Ter",3) || strncmp (cp, cq,3))
	    {
	      k++ ;
	      for (j = 0 ; j < 3 ; j++)
		buf0[j] = cp >= pRef+3 ? cp[j-3] : '-' ;
	      buf0[3] = 0 ;
	      for (j = 0 ; j < 3 ; j++)
		buf[j] = cp[j] ;
	      buf[3] = 0 ;
	      for (j = 0 ; j < 3 ; j++)
		buf2[j] = cq[j] ;
	      buf2[3] = 0 ;
	      for (j = 0 ; j < 3 ; j++)
		buf12[j] = cp[j+3] ;
	      buf12[3] = 0 ;
	      for (j = 0 ; j < 3 ; j++)
		buf22[j] = cq[j+3] ;
	      buf22[3] = 0 ;
	      if (m == 1 && ! strncmp (buf, "Met", 3) && strncmp (buf2, "Met", 3))
		{ 
		  int x = 3*(m-1) + p1 ; /* first base of modified triplet */
		  int met = 0 ;
		  int dx = 1 + (fs < 0 ? (-fs-1)/3 : 0) ;
		  
		  locateMet(Rdna, translationTable, x + 3*dx, &met) ;
		  if (met + dx == 0)
		    {
		      vtxtPrintf (txt, "%sMet2del", sep) ;
		      vtxtPrintf (txt, "\nMet1_lost Met2del") ;
		    }
		  else if (met + dx > 0)
		    {
		      vtxtPrintf (txt, "%s%s2_Met%ddel", sep, buf12, met+dx+1) ;
		      vtxtPrintf (txt, "\nMet1_lost \"%s2_Met%ddel\"", buf12, met+dx+1) ;
		    }
		  else
		    {
		      int stop1 = 0, stop2 = 0 ;
		      locateStop(Rdna, translationTable, x, 0, &stop1, &stop2) ;
		      vtxtPrintf (txt, "%sMet1_Ter%ddel", sep, stop2 + 1) ;
		      vtxtPrintf (txt, "\nMet1_lost \"Met1%s_NoOtherMet\"", buf2) ;
		    }
		  break ;
		}
	      else if (strncmp (buf, "Ter", 3) && ! strncmp (buf2, "Ter", 3))
		{ 
		  vtxtPrintf (txt, "%s%s%dTer", sep, buf,m) ;
		  vtxtPrintf (txt, "\nAA_to_Stop \"%s%dTer\"", buf,m) ;
		  break ;
		}
	      else if (! strncmp (buf, "Ter", 3) && strncmp (buf2, "Ter", 3))
		{ 
		  int x = 3*(m-1) + p1 ; /* first base of modified triplet */
		  int stop1 = 0, stop2 = 0 ;
		  int dx = fs < 0 ? (-fs + 2)/3 : 0 ; dx++ ;
		  locateStop (Rdna, translationTable, 0*(x + 3*dx), x-fs+3*dx, &stop1, &stop2) ;
		  if (stop2>0)
		    {
		      vtxtPrintf (txt, "%s%s%d%sextTer%d", sep, buf,m,buf2,stop2+1) ;
		      vtxtPrintf (txt, "\nStop_to_AA \"%s%d%sextTer%d\"",buf,m,buf2,stop2+1) ;
		      break ;
		    }
		  else
		    {
		      vtxtPrintf (txt, "%s%s%d%sextTer?", sep, buf,m,buf2) ;
		      vtxtPrintf (txt, "\nStop_to_AA \"%s%d%sextTer?\"", buf,m,buf2) ;
		      break ;
		    }
		}
	      else if (! strncmp (buf, "Ter", 3) && ! strncmp (buf2, "Ter", 3) && m > 1)
		{
		  vtxtPrintf (txt, "%s%s%d=", sep, buf,m) ;
		  vtxtPrintf (txt, "\nSynonymous \"%s%d=\"", buf,m) ;
		  break ;
		}
	    }
	  if (strncmp (cp, cq,3))
	    {
	      if (fs == 0)
		{
		  vtxtPrintf (txt, "%s%s%d%s", sep, buf, m, buf2) ;
		  vtxtPrintf (txt, "\nAA_substitution \"%s%d%s\"", buf, m, buf2) ;
		  break ;
		}
	      else
		{
		  int mk = vtxtMark (txt) ;
		  if (fs == -3 && !strcmp (buf2,buf12)) /* del in frame */
		    vtxtPrintf (txt, "%ddel%s", m, buf) ;
		  else if (fs == -3) /* del in frame */
		    vtxtPrintf (txt, "%s%d_%s%ddelins%s", buf,m,buf12,m+1, buf2) ;
		  else if (fs == 3 && !strcmp (buf,buf22)) /* ins in frame */
		    vtxtPrintf (txt, "%s%d_%s%dins%s", buf0,m-1,buf,m, buf2) ;
		  else if (fs == 3)
		    vtxtPrintf (txt, "%s%ddelins%s%s",buf, m, buf2, buf22) ;
		  else if (fs <0 && (-fs % 3 == 0))
		    vtxtPrintf (txt, "%d_%ddel", m, m - fs/3 -1) ;
		  else if (fs >0 && (fs % 3 == 0))
		    vtxtPrintf (txt, "%d_%dins(%d)", m, m + 1, fs/3) ;
		  else if (fs %3)
		    {
		      char *ter = strstr (cq, "Ter") ;
		      int x = 3*(m-1) + p1 ; /* first base of modified triplet */
		      int stop1 = 0, stop2 = 0 ;
		      int dx = fs < 0 ? (-fs + 2)/3 : 0 ;

		      if (ter)
			vtxtPrintf (txt, "%s%d%sfsTer%d", buf, m, buf2, 1+(ter-cq)/3) ;
		      else if (Rdna && RdnaLn > x+3*dx && RdnaLn > x+3*dx-fs)
			{ /* the snippet is too short, look in the mrna */
			  /* compute now the location of the stops 
			   * in the complete mRNA and modified mRNA 
			   */
			  locateStop (Rdna, translationTable, x+3*dx, x - fs+3*dx, &stop1, &stop2) ;
			  if (stop2 > 0)
			    vtxtPrintf (txt, "%s%d%sfsTer%d", buf, m, buf2, stop2 + 1 + dx) ;
			  else
			    vtxtPrintf (txt, "%s%d%sfsTer?", buf, m, buf2) ;
			  if (stop2 > stop1)
			    vtxtPrintf (txt, "\nExtension %d AA", stop2 - stop1) ;
			}
		      else
			vtxtPrintf (txt, "%s%d%sfs", buf,m,buf2) ;
		    }
		  if (fs % 3 == 0)
		    {
		      char *cr = strnew (vtxtAt (txt, mk), 0) ;
		      vtxtPrintf (txt, "\nFrame_preserving_indel \"%s\"", cr) ;
		      ac_free (cr) ;
		    }
		  else
		    {
		      char *cr = strnew (vtxtAt (txt, mk), 0) ;
		      vtxtPrintf (txt, "\nFrameshift \"%s\"", cr) ;
		      ac_free (cr) ;
		    }
		  break ; 
		}
	    }
	  cp += 3 ; cq += 3 ; i += 3 ;
	}
      if (k == 0) 
	{
	  vtxtPrintf (txt, "%s%d=", idem,m1) ;
	  vtxtPrintf (txt, "\nSynonymous %s%d=", idem, m1) ;
	}
    }
  /* clean up the pR1 pV1 */
  if (m1 > 0)
    {
      cp = strchr (pR1 + 6, '*') ;
      if (cp) *++cp = 0 ;
      cp = strchr (pV1 + 6, '*') ;
      if (cp) *++cp = 0 ;
      cp = strstr (pRef + 18, "Ter") ;
      if (cp) cp[3] = 0 ;
      cp = strstr (pVar + 18, "Ter") ;
      if (cp) cp[3] = 0 ;
    }
  return ; 
}  /* tsnpSetPName */

/*************************************************************************************/
/* biologits, like Socrates, never heard of zero */
static int socrate (int x)
{
  return x > 0 ? x : x - 1 ; 
} /* socrate */

/*************************************************************************************/

static int tsnpSetGName (vTXT txt, TSNP *tsnp, AC_OBJ Snp, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl = ac_tag_table (Snp, "IntMap", h) ;
  AC_TABLE Rtbl = ac_tag_table (Snp, "mRNA", h) ;
  KEY seq = ac_table_key (tbl, 0, 0, 0) ;
  KEY mrna = ac_tag_key (Snp, "mRNA", 0) ;
  static KEY oldSeq = 0 ;
  static KEY oldMrna = 0 ;
  static const char *dna = 0 ;
  BOOL ok = FALSE ;
  char bufR[52], bufV[52], bufVG[52] ;
  char RbufR[52], RbufV[52] ;
  char bufSeqTitle[128] ;  
  fromMrna = ac_has_tag (Snp, "Found_in_mRNA") ;
  /*
    KEY product = 0 ; 
    int mrnaA1 = 0 ;   // a coord of the first base of the mRNA in the map 
  int productX1 = 0 ; // x coord of the first base of the coding region in the mRNA 
  */
  memset (bufR, 0, sizeof (bufR)) ;
  memset (bufV, 0, sizeof (bufV)) ;
  memset (RbufR, 0, sizeof (bufR)) ;
  memset (RbufV, 0, sizeof (bufV)) ;

  socrate (0) ; /* for compiler happiness */

  if (Rtbl)
    {
      AC_OBJ Mrna = ac_table_obj (Rtbl, 0, 0, h) ;
      if (Mrna)
	{
	  const char *ccp = ac_tag_printable (Mrna, "Gene", 0) ;
	  if (ccp)
	    vtxtPrintf (txt, "Gene \"%s\"\n", ccp) ;
	  ac_free (Mrna) ;
	}
    }
  
  bufSeqTitle[0]= 0 ; 
  if (seq) 
    {
      AC_OBJ Seq = ac_table_obj (tbl, 0, 0, h) ;
      const char *ccp = Seq ?  ac_tag_printable (Seq, "Title", name(seq)) : 0 ;
      if (ccp && *ccp)
	strncpy (bufSeqTitle, ccp, 126) ;
    }

  if (! seq)
    {
      if (dna != Rdna) ac_free (dna) ;
      dna = 0 ;
      dnaLn = 0 ;
      oldSeq = 0 ;      
      ac_free (Gdna) ;
    }
  else if (seq != oldSeq && ! fromMrna)
    {
      AC_OBJ Seq = ac_tag_obj (Snp, "Parent_sequence", h) ;

      ac_free (dna) ;
      oldSeq = seq ;
	
      dna = Seq ? ac_obj_dna (Seq, h0) : 0 ;
      dnaLn = dna ? strlen (dna) : 0 ;
      ac_free (Gdna) ;
      Gdna = dnaGet (seq) ;
    }

  if (! mrna)
    {
      if (RdnaLn) ac_free (Rdna) ;
      Rdna = dna ;
      RdnaLn = 0 ;
      oldMrna = 0 ;
      if (!seq) 
	return 0 ;
    }
  else if (mrna != oldMrna)
    {
      AC_OBJ Mrna = ac_tag_obj (Snp, "mRNA", h) ;

      if (RdnaLn) ac_free (Rdna) ;
      oldMrna = mrna ;
      Rdna = Mrna ? ac_obj_dna (Mrna, h0) : 0 ;
      RdnaLn = Rdna ? strlen (Rdna) : 0 ;
      ac_free (Mrna) ;
    }
  myMrna = mrna  ;
  if (fromMrna)
    {
      dna = Rdna ;
      dnaLn = RdnaLn ;
    }

  if (! dna)
    return 0 ;
  if (tbl && tbl->cols >= 3)
    {

      const char *ccp = ac_table_printable (tbl, 0, 0, "xxx") ;
      int i, j ;
      int fs = 0 ; /* frameshift */
      int dda = 0 ; /* use to get the translation frame */
      int m1 = ac_table_int (Rtbl, 0, 1, 0) ;
      int m2 = ac_table_int (Rtbl, 0, 2, 0) ;
      int da = 1 ;
      int a1 = ac_table_int (tbl, 0, 1, 0) ;
      int a2 = ac_table_int (tbl, 0, 2, 0) ;
      char **sub, *subs[] = {"A2T","A2G","A2C","T2A","T2G","T2C","G2A","G2T","G2C","C2A","C2T","C2G",0} ;
      char **del, *dels[] = {"DelA","DelT","DelG","DelC", 0} ;
      char **ins, *inss[] = {"InsA","InsT","InsG","InsC", 0} ;
      g2m = RdnaLn ? m1 - a1 : 0 ;
      
      vtxtPrintf (txt, "%s\n", a1 < a2 ? "Forward" : "Reverse") ;
      if (a2 == 1 || a2 == -1)
	{
	  char *cp = strnew (ac_name(Snp), 0), *cq, *cr ;
	  cq = strchr (cp, ':') ;
	  if (cq)
	    {
	      cq++ ; 
	      cr = strchr (cq, '_') ;
	      if (cr)
		{
		  da = 0 ;
		  cr++ ;
		  while (*cr >= '0' && *cr <= '9')
		    { da = 10 * da + (*cr - '0') ; cr++ ; }
		}
	    }	  
	  if (a2 == 1)
	    a2 = a1 + da + 1 ;
	  else
	    a2 = a1 - (da + 1) ;
	  if (m2 == 1)
	    m2 = m1 + da + 1 ;
	  else
	    m2 = m1 - (da + 1) ;
	  ac_free (cp) ;
	}      
      if (strstr (ac_name(Snp), ":Sub"))
	{
	  if (! ac_has_tag (Snp, "Multi_substitution"))
	    {
	      for (sub = subs ; *sub ; sub++)
		{
		  if (ac_has_tag (Snp, *sub))
		    {
		      int am1 = fromMrna ? m1 : a1 ;
		      int a3 = a1 < a2 ? a1 + 1 : a1 - 1 ;
		      char buf[4] ;
		      if  (a1 < a2)
			{
			  buf[0] = ace_upper (sub[0][0]) ;
			  buf[1] = '>' ;
			  buf[2] = ace_upper(sub[0][2]) ;
			  buf[3] = 0 ;
			}
		      else
			{
			  buf[0] = ace_upper (complementLetter(sub[0][0])) ;
			  buf[1] = '>' ; 
			  buf[2] = ace_upper (complementLetter(sub[0][2])) ;
			  buf[3] = 0 ;
			}

		      ok = TRUE ;
		      vtxtPrint (txt, "-D Substitution\n") ; /* cleanup */
		      vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Base %d %c becomes %c\"\n", name (seq), a1, a2, a3, buf[0], buf[2]) ;
		      vtxtPrintf (txt,"Typ %c>%c\n%s\n", sub[0][0], sub[0][2], *sub) ; /* reinstate */ 
		      if (a1 < a2)
			{
			  vtxtPrintf (txt, "gName \"%s:g.%d%s\"\n", bufSeqTitle, a1+1, buf) ;
			  vtxtPrintf (txt, "VCF %s %d %c %c\n", name(seq), a1+1, ace_upper (sub[0][0]), ace_upper (sub[0][2])) ;
			}
		      else
			{
			  vtxtPrintf (txt, "gName \"%s:g.%d%s\"\n", bufSeqTitle, a1-1, buf) ;
			  vtxtPrintf (txt, "VCF %s %d %c %c\n", name(seq), a2 + 1
				      , ace_upper (complementLetter(sub[0][0]))
				      , ace_upper (complementLetter(sub[0][2]))
				      ) ;
			}
		      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d%c>%c\"\n", name(mrna), am1 + 1, sub[0][0], sub[0][2]) ;
		      
		      dda = 0 ; /* -1 for Socrates, but +1 because the sub is between a1 and a2, at a1+1 */

		      for (i = am1 - 21 + dda, j = 0 ; i < am1 + 30 + dda && i < dnaLn ; i++)
			if (i > 0 && j < 51) 
			  {
			    if (i == am1)
			      { 
				bufR[j] = ace_upper (dna[i]) ;
				RbufR[j] = ace_upper (myRdna(i)) ;
				bufV[j] = sub[0][2] ;
				RbufV[j] = sub[0][2] ;
			      }
			    else
			      {
				bufV[j] = bufR[j] = dna[i] ;
				RbufV[j] = RbufR[j] = myRdna(i) ;
			      }
			    j++ ;
			  }
		      RbufR[j] = bufR[j] = 0 ;
		      RbufV[j] = bufV[j] = 0 ;
		    }
		}
	    }
	  else
	    {
	      tbl = ac_tag_table (Snp, "Multi_substitution", h) ;
	      if (tbl->cols >= 3)
		{
		  int k,  da = a2 - a1 - 1 ;
		  const char *ccR = ac_table_printable (tbl, 0, 1, "") ; 
		  const char *ccV = ac_table_printable (tbl, 0, 2, "") ; 
		  int am1 = fromMrna ? m1 : a1 ;
		  int am2 = fromMrna ? m2 : a2 ;

		  vtxtPrint (txt, "-D Substitution\n") ; /* cleanup */
		  if (a1 < a2)
		    vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Bases %d to %d (%d bases) are modified\"\n", name (seq), a1, a2, a1+1, a2-1, da) ;
		  else
		    vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Bases %d to %d (%d bases) are modified\"\n", name (seq), a2, a1, a2+1, a1-1, da) ;
		  vtxtPrintf (txt, "Multi_substitution %d %s %s\n", da, ccR, ccV) ; /* reinstate */
		  if (ccV && strlen (ccV) == da)
		    {
		      ok = TRUE ;
		      if (a1 < a2)
			{
			  vtxtPrintf (txt, "VCF %s %d %s %s\n", name(seq), a1+1, ccR, ccV) ;
			  vtxtPrintf (txt, "gName \"%s:g.%d_%ddelins%s\"\n", bufSeqTitle, a1+1,a2 -1, ccV) ;
			}
		      else
			{
			  int i, n = strlen (ccV) ;
			  char buf[n+1] ;

			  for (i = 0 ; i < n ; i++)
			    buf[i] = complementLetter (ccV[n-1-i]) ;
			  buf[i] = 0 ;
			  vtxtPrintf (txt, "VCF %s %d %s %s\n", a2 + 1, name(seq), ccR, buf) ;
			  vtxtPrintf (txt, "gName \"%s:g.%d_%ddelins%s\"\n", bufSeqTitle, a2+1,a1 -1, buf) ;
			}

		      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddelins%s\"\n", name(mrna),am1+1,am2-1, ccV) ;
		      vtxtPrintf (txt, "Typ \"%s2%s\"\n", ccR, ccV) ;
		      dda = 0 ; /* -1 for Socrates, but +1 because the sub is between a1 and a2, starting at a1+1 */
		      
		      for (i = am1 - 21 + dda, j = k = 0 ; i < am1 + 30 + dda && i < dnaLn ; i++)
			if (i > 0 && j < 51)
			  {
			    if (i >= am1 -1 && i <= am2 - 2)
			      { 
				bufR[j] = ace_upper (dna[i]) ;
				if (a1 < a2)
				  bufV[j] = ace_upper (ccV[i - a1 + 1]) ;
				else
				  bufV[j] = ace_upper (complementLetter(ccV[da - (i - a1) - 1 + 1])) ;
				RbufR[j] = ace_upper (myRdna(i)) ;
				RbufV[j] = ace_upper (ccV[i - a1 + 1]) ;
				if (bufR[j] == bufV[j])
				  bufR[j] = bufV[j] = dna[i] ; /* lower the common bases */
				if (RbufR[j] == RbufV[j])
				  RbufR[j] = RbufV[j] = myRdna(i) ; /* lower the common bases */
			      }
			    else
			      {
				bufV[j] = bufR[j] = dna[i] ;
				RbufV[j] = RbufR[j] = myRdna(i) ;
			      }
			    j++ ;
			  }
		      RbufR[j] = bufR[j] = 0 ;
		      RbufV[j] = bufV[j] = 0 ;
		    }
		}
	    }
	}	    
      
      else if (ac_has_tag (Snp, "Deletion"))   /* happens because dimA is under DelA and should be under Diminution */
	{
	  if (! ac_has_tag (Snp, "Multi_deletion"))   /* happens because dimA is under DelA and should be under Diminution */
	    {
	      for (del = dels ; *del ; del++)
		if (ac_has_tag (Snp, *del))
		  {
		    int da = 1 ;
		    int dx = 0 ;
		    int slide = 0 ;
		    int am2, am1 = fromMrna ? m1 : a1 ;
		    BOOL isDim = tsnpSlide (dna, dnaLn, am1, da, &dx, &slide) ;
		    if (a1 < a2) { a1 += dx ; a2 += dx ;}
		    else { a1 -= dx ; a2 -= dx ;}
		    if (fromMrna) { m1 += dx ; m2 += dx ; }
		    am1 = fromMrna ? m1 : a1 ;
		    am2 = am1 + 2 ;
		    fs = -1 ; 

		    if (mrna) vtxtPrintf (txt, "mRNA %s %d %d \"Base %c %d is deleted\"\n", name (mrna), m1, m2, (*del)[3], m1+1) ;
		    if (a1 < a2)
		      vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Base %c %d is deleted\"\n", name (seq), a1, a2, (*del)[3], a1+1) ;
		    else
		      vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Base %c %d is deleted\"\n", name (seq), a1, a2, ace_upper(complementLetter(dna[am1])), a1  - 1) ;
		    vtxtPrintf (txt, "-D Sliding\n") ;
		    vtxtPrint (txt, "-D Deletion\n") ; /* cleanup */
		    vtxtPrintf (txt,"%s\n", *del) ; /* reinstate */
		    
		    
		    ok = TRUE ;
		    if (isDim)
		      {
			if (a1 < a2)
			  vtxtPrintf (txt, "gName \"%s:g.%ddim%c\"\n", bufSeqTitle, a1+1+slide, ace_upper(dna[am1])) ;
			else
			  vtxtPrintf (txt, "gName \"%s:g.%ddim%c\"\n", bufSeqTitle, a2+1, ace_upper(complementLetter(dna[am1]))) ;
			if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%ddim%c\"\n", name(mrna), am1+1+slide, ace_upper(dna[am1])) ;
			vtxtPrintf (txt, "Typ \"Dim%c\"\n", ace_upper(dna[am1])) ;
			vtxtPrintf (txt, "Diminution\nIn_repeat\n") ;
		      }
		    else
		      {
			if (a1 < a2)
			  vtxtPrintf (txt, "gName \"%s:g.%ddel%c\"\n", bufSeqTitle, a1+1, ace_upper(dna[am1])) ;
			else
			  vtxtPrintf (txt, "gName \"%s:g.%ddel%c\"\n", bufSeqTitle, a2+1, ace_upper(complementLetter(dna[am2+slide]))) ;
			if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%ddel%c\"\n", name(mrna), am1+1+slide, ace_upper(dna[am1])) ;
			vtxtPrintf (txt, "Typ \"Del%c\"\n", ace_upper(dna[am1])) ;
		      }
		    if (slide) 
		      vtxtPrintf (txt, "Sliding %d\n", slide) ;
		    if (a1 < a2)
		      vtxtPrintf (txt, "VCF %s %d %c%c %c\n", name(seq),  a1,  ace_lower(dna[am1-1]), ace_upper(dna[am1]), ace_lower(dna[am1-1])) ;
		    else
		      vtxtPrintf (txt, "VCF %s %d %c%c %c\n", name(seq),  a2 - slide,  complementLetter(dna[am2+slide-1]), ace_upper(complementLetter(dna[am2+slide-2])), complementLetter(dna[am2+slide-1])) ;

		    am1 = fromMrna ? m1 : a1 ;
		    for (i = am1 - 21, j = 0 ; i < am1 + 30 + da && i + da < dnaLn ; i++)
		      if (i > 0 && j < 51)
			{
			  if (i < am1)
			    {
			      bufV[j] = bufR[j] = dna[i] ;
			      RbufV[j] = RbufR[j] = myRdna(i) ;
			    }
			  else if (i == am1)
			    {
			      bufR[j] = ace_upper (dna[i]) ;
			      bufV[j] = '_' ; 
			      RbufR[j] = ace_upper (myRdna(i)) ;
			      RbufV[j] = '_' ; 
			    }
			  else if (i > am1)
			    {
			      bufR[j] =  bufV[j] = dna[i] ;
			      RbufV[j] = RbufR[j] = myRdna(i) ;
			    }
			  j++ ;
			}
		    RbufR[j] = bufR[j] = 0 ;
		    RbufV[j] = bufV[j] = 0 ;
		  }
	    }
	  else
	    {
	      tbl = ac_tag_table (Snp, "Multi_deletion", h) ;
	      if (tbl)
		{
		  int da = (a1 < a2 ? a2 - a1 - 1 : a1 - a2 - 1) ;
		  int dx = 0, k, kk ;
		  int slide = 0 ;
		  int am1 = fromMrna ? m1 : a1 ;
		  int am2 = fromMrna ? m2 : a2 ;
		  BOOL isDim = tsnpSlide (dna, dnaLn, am1, da, &dx, &slide) ;
		  char bufN[25] ;
		  char bufVS[30] ;
		  char bufVGS[30] ;
		  /* 		  char bufVC[51] ; */
		  if (a1 < a2) { a1 += dx ; a2 += dx ;}
		  else { a1 -= dx ; a2 -= dx ;}
		  if (fromMrna) { m1 += dx ; m2 += dx ; }
		  am1 = fromMrna ? m1 : a1 ;
		  am2 = fromMrna ? m2 : a2 ;
		  fs = -da ;
		  vtxtPrintf (txt, "-D Sliding\n") ;
		  vtxtPrint (txt, "-D Deletion\n") ; /* cleanup */
		  
		  if (a1 < a2)
		    vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Bases %d to %d (%d bases) are deleted\"\n", name (seq), a1, a2, a1+1, a2-1, da) ;
		  else
		    vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Bases %d to %d (%d bases) are deleted\"\n", name (seq), a1, a2, a2+1, a1-1, da) ;
		  ok = TRUE ;
		  /* temporarily store the deleted bases in bufV */
		  for (i = am1, j = 0 ; j < 51 && i < am2 - 1 && i < dnaLn ; i++)
		    bufV[j++] = ace_upper(dna[i]) ;
		  bufV[j] = 0 ;
		  for (i = am1, j = 0 ; j < 51 && i < am2 - 1 && i < dnaLn ; i++)
		    bufVS[j++] = ace_upper(dna[i+slide]) ;
		  bufVS[j] = 0 ;
		  if (a1 > a2)
		    {
		      for (i = 0 ; i <  da ; i++)
			bufVG[i] = ace_upper(dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)bufV[da - 1 - i]])]) ;
		      bufVG[i] = 0 ;
		      for (i = 0 ; i <  da ; i++)
			bufVGS[i] = ace_upper(dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)bufVS[da - 1 - i]])]) ;
		      bufVGS[i] = 0 ;
		    }
		  vtxtPrintf (txt,"Multi_deletion %d %s\n", da, bufV) ; /* reinstate */ 
		  if (da <= 20)
		    {
		      if (isDim)
			{
			  if (a1 < a2)
			    {
			      vtxtPrintf (txt, "VCF %s %d %c%s %c\n", name(seq), a1, ace_lower(dna[am1-1]), bufV, ace_lower(dna[am1-1])) ;
			      vtxtPrintf (txt, "gName \"%s:g.%d_%ddim%s\"\n", bufSeqTitle, a1+1+slide,a2-1+slide, bufVS) ;
			    }
			  else
			    {
			      char cc = ace_lower(complementLetter(dna[am2+slide-1])) ;
			      vtxtPrintf (txt, "VCF %s %d %c%s %c\n", name(seq), a2-slide, cc, bufVGS, cc) ; 
			      vtxtPrintf (txt, "gName \"%s:g.%d_%ddim%s\"\n", bufSeqTitle, a2+1,a1-1, bufVG) ;
			    }
			  if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddim%s\"\n", name(mrna), am1 + 1 + slide, am2 - 1 + slide, bufVS) ;
			  vtxtPrintf (txt, "Typ \"Dim%s\"\n", bufV) ;
			  vtxtPrintf (txt, "Diminution\n") ;
			}
		      else 
			{
			  if (a1 < a2)
			    {
			      vtxtPrintf (txt, "VCF %s %d %c%s %c\n", name(seq), a1, ace_lower(dna[am1-1]), bufV, ace_lower(dna[am1-1])) ;
			      vtxtPrintf (txt, "gName \"%s:g.%d_%ddel%s\"\n", bufSeqTitle, a1+1+slide,a2-1+slide, bufVS) ;
			      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddel%s\"\n", name(mrna), am1 + 1 + slide, am2 - 1 + slide, bufVS) ;
			    }
			  else
			    {
			      char cc = ace_lower(complementLetter(dna[am2+slide-1])) ;
			      vtxtPrintf (txt, "VCF %s %d %c%s %c\n", name(seq), a2-slide, cc, bufVGS, cc) ; 
			      vtxtPrintf (txt, "gName \"%s:g.%d_%ddel%s\"\n", bufSeqTitle, a2+1,a1-1, bufVG) ;
			      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddel%s\"\n", name(mrna), am2 + 1 - slide, am1 - 1 - slide, bufVS) ;
			    }

			  vtxtPrintf (txt, "Typ \"Del%s\"\n", bufV) ;
			}
		      vtxtPrintf (txt, "Multi_deletion %d %s\n", da, bufV) ;
		    }
		  else
		    {
		      vtxtPrintf (txt, "Typ \"Del_%d\"\n", da) ;
		      if (a1 < a2)
			{
			  vtxtPrintf (txt, "VCF %s %d %c%d %c\n", name(seq), a1, ace_lower(dna[am1-1]), da, ace_lower(dna[am1-1])) ;
			  vtxtPrintf (txt, "gName \"%s:g.%d_%ddel(%d)\"\n", bufSeqTitle, a1+1,a2-1, da) ;
			}
		      else
			{
			  char cc = ace_lower(complementLetter(dna[am2+slide-1])) ;
			  vtxtPrintf (txt, "VCF %s %d %c%d %c\n", name(seq), a1, cc, da, cc) ;
			  vtxtPrintf (txt, "gName \"%s:g.%d_%ddel(%d)\"\n", bufSeqTitle, a1+1,a2-1, da) ;
			}

		      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddel(%d)\"\n", name(mrna), am1 + 1 + slide, am2 - 1 + slide, da) ;
		      vtxtPrintf (txt, "Multi_deletion %d\n", da) ;
		    }
		  if (slide) 
		    vtxtPrintf (txt, "Sliding %d\n", slide) ;
		  bufN[0] = 0 ;
		  if (da ==  1) sprintf (bufN, "_") ; 
		  else if (da ==  2) sprintf (bufN, "__") ; 
		  else if (da ==  3) sprintf (bufN, "___") ; 
		  else if (da ==  4) sprintf (bufN, "____") ; 
		  else if (da ==  5) sprintf (bufN, "_____") ; 
		  else if (da ==  6) sprintf (bufN, "______") ; 
		  else if (da ==  7) sprintf (bufN, "_______") ; 
		  else if (da ==  8) sprintf (bufN, "________") ; 
		  else if (da ==  9) sprintf (bufN, "_________") ; 
		  else if (da == 10) sprintf (bufN, "__________") ; 
		  else if (da == 11) sprintf (bufN, "___________") ; 
		  else if (da == 12) sprintf (bufN, "____________") ; 
		  else if (da > 12) sprintf (bufN, "_%d_", da) ; 
		  
		  am1 = fromMrna ? m1 : a1 ;
		  am2 = fromMrna ? m2 : a2 ;
		  for (i = am1 - 21, k = kk = j = 0 ; i < am1 + 30 && i < dnaLn ; i++)
		    if (i > 0 && k < 51)
		      {
			if (i < am1)
			  {
			    bufV[k] = bufR[j] = dna[i] ;
			    RbufV[k] = RbufR[j] = myRdna(i) ;
			    k++ ; j++ ;
			  }
			else if (i >= am1 && i < am2 - 1)
			  {
			    bufR[j] = ace_upper (dna[i]) ;
			    RbufR[j] = ace_upper (myRdna(i)) ;
			    j++ ;
			    if (bufN[kk]) 
			      { RbufV[k] = bufV[k] = bufN[kk++] ; k++ ; } 
			    else
			      {
				bufV[k]  = dna[i+da-kk] ;		    
				RbufV[k] = myRdna(i+da-kk) ;
				k++ ;
			      }
			  }
			else if (i >= am2 -1)
			  {
			    if (j < 51) { bufR[j] = dna[i] ; RbufR[j] = myRdna(i) ;j++ ; }
			    bufV[k]  = dna[i+da-kk] ;		    
			    RbufV[k] = myRdna(i+da-kk) ;
			    k++ ;
			  }
		      }
		  RbufR[j] = bufR[j] = 0 ;
		  RbufV[k] = bufV[k] = 0 ;
		}
	    }
	}
      else if (ac_has_tag (Snp, "Insertion"))
	{  
	  if (! ac_has_tag (Snp, "Multi_insertion"))   /* happens because DupA is under InsA and should be under Duplication */
	    {
	      for (ins = inss ; *ins ; ins++)
		if (ac_has_tag (Snp, *ins))
		  {
		    int k ;
		    int da = 1 ;
		    int dx = 0 ;
		    BOOL isDup = FALSE ;
		    char buf[da+1] ;
		    char bufS[da+1] ;
		    char bufG[da+1] ;
		    char bufGS[da+1] ;
		    int slide = 0 ;
		    int am1 = fromMrna ? m1 : a1 ;
		    int am2 = fromMrna ? m2 : a2 ;
		    isDup = tsnpSlideDup (dna, dnaLn, am1, da, &dx, &slide, (*ins) + 3, buf, bufS, bufG, bufGS) ;
		    if (a1 < a2) { a1 += dx ; a2 += dx ;}
		    else { a1 -= dx ; a2 -= dx ;}
		    if (fromMrna) { m1 += dx ; m2 += dx ; }
		    am1 = fromMrna ? m1 : a1 ;
		    am2 = fromMrna ? m2 : a2 ;
		    ok = TRUE ;
		    ok = TRUE ;
		    fs = 1 ;

		    if (a1 < a2)
		      vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Base %c is inserted on plus strand between base %d and %d\"\n", name (seq), a1, a2, (*ins)[3], a1, a2) ;
		    else
		      vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Base %c is inserted on minus strand between base %d and %d\"\n", name (seq), a1, a2, ace_upper(complementLetter((*ins)[3])), a1, a2) ;
		    if (a1 < a2)
		      vtxtPrintf (txt, "VCF %s %d %c %c%c\n", name(seq), a1, ace_lower(dna[am1-1]), ace_lower(dna[am1-1]), (*ins)[3]) ;
		    else
		      vtxtPrintf (txt, "VCF %s %d %c %c%c\n", name(seq), a2-slide, ace_lower(complementLetter(dna[am1+slide])), ace_lower(complementLetter(dna[am1+slide])), ace_upper(complementLetter( (*ins)[3]))) ;
		    vtxtPrintf (txt, "-D Sliding\n") ;
		    vtxtPrint (txt, "-D Insertion\n") ; /* cleanup */
		    vtxtPrintf (txt,"%s\n", *ins) ; /* reinstate */
		    if (isDup)
		      {
			if (a1 < a2)
			  vtxtPrintf (txt, "gName \"%s:g.%ddup%c\"\n", bufSeqTitle, a1+slide, (*ins)[3]) ;
			else
			  vtxtPrintf (txt, "gName \"%s:g.%ddup%c\"\n", bufSeqTitle, a2, ace_upper(complementLetter ((*ins)[3]))) ;
			vtxtPrintf (txt, "Typ \"Dup%c\"\n", (*ins)[3]) ;
			if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%ddup%c\"\n", name(mrna), am1+slide, (*ins)[3]) ;
			vtxtPrintf (txt, "Duplication\nIn_repeat\n") ;
		      }
		    else
		      {
			if (a1 < a2)
			  vtxtPrintf (txt, "gName \"%s:g.%d_%dins%c\"\n", bufSeqTitle, a1+slide, a2+slide, (*ins)[3]) ;
			else
			  vtxtPrintf (txt, "gName \"%s:g.%d_%dins%c\"\n", bufSeqTitle, a2, a1, ace_upper(complementLetter ((*ins)[3]))) ;
			vtxtPrintf (txt, "Typ \"Ins%c\"\n", (*ins)[3]) ;
			if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%dins%c\"\n", name(mrna), am1+slide, am2+slide, (*ins)[3]) ;
		      }
		      
		    if (slide) 
		      vtxtPrintf (txt, "Sliding %d\n", slide) ;
		    
		    for (i = am1 - 21, j = k = 0 ; i < am1 + 30 && i < dnaLn ; i++)
		      if (i > 0 && j < 51) 
			{
			  bufR[j] = dna[i] ; RbufR[j] = myRdna(i) ;  j++ ;
			  if (k < 51) { bufV[k] = dna[i] ; RbufV[k] = myRdna(i) ; k++ ; }
			  if (i == am1-1)
			    { 
			      bufR[j] = '^' ; RbufR[j] = '^' ; j++ ;
			      ccp = buf ; while (*ccp) { bufV[k] = ace_upper(*ccp) ; RbufV[k] = ace_upper(*ccp) ; k++ ; ccp++ ; }
			    }
			}
		    RbufR[j] = bufR[j] = 0 ;
		    RbufV[j] = bufV[j] = 0 ;
		  }
	    }
	  else
	    {
	      tbl = ac_tag_table (Snp, "Multi_insertion", h) ;
	      if (tbl)
		{
		  int am1, am2, da = ac_table_int (tbl, 0, 0, 0) ;
		  const char *ccp = ac_table_printable (tbl, 0, 1, 0) ; 
		  if (0 && ! ccp) {
		    ccp = ac_name (Snp) ;
		    if (ccp) ccp = strchr (ccp, ':') ;
		    if (ccp) ccp = strchr (ccp+1, ':') ;
		    if (ccp) ccp = strchr (ccp+1, ':') ;
		    if (ccp) ccp = strchr (ccp+1, ':') ;
		    if (ccp) ccp = ccp+2 ;
		    da = strlen (ccp) ;
		  }
		  if (a2 - a1 == 1 || a2 - a1 == -1)
		    {
		      int i, k, dx = 0 ;
		      BOOL isDup = FALSE ;
		      char buf[da+1], bufN[15] ;
		      char bufS[da+1] ;
		      char bufG[da+1] ;
		      char bufGS[da+1] ;
		      int slide = 0 ;
		      am1 = fromMrna ? m1 : a1 ;
		      am2 = fromMrna ? m2 : a2 ;

		      if (ccp && da == strlen (ccp))
			isDup = tsnpSlideDup (dna, dnaLn, am1, da, &dx, &slide, ccp, buf, bufS, bufG, bufGS) ;
		      else
			{
			  int i ;
			  for (i = 0 ; i < da ; i++)
			    buf[i] = bufG[i] = 'n' ;
			  buf[i] = bufG[i] = 0 ;
			  ok = FALSE ;
			}

		      if (a1 < a2) { a1 += dx ; a2 += dx ;}
		      else { a1 -= dx ; a2 -= dx ;}
		      if (fromMrna) { m1 += dx ; m2 += dx ; }
		      fs = da ;
		      am1 = fromMrna ? m1 : a1 ;
		      am2 = fromMrna ? m2 : a2 ;

		      if (a1 < a2)
			vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"%d bases %s are inserted on plus strand between base %d and %d\"\n", name (seq), a1, a2, da, buf, a1, a2) ;
		      else
			vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"%d bases %s are inserted on minus strand between base %d and %d\"\n", name (seq), a1, a2, da, bufG, a1, a2) ;
		      
		      if (a1 < a2)
			vtxtPrintf (txt, "VCF %s %d %c %c%s\n", name(seq), a1, ace_lower(dna[am1-1]), ace_lower(dna[am1-1]), buf) ;
		      else
			{
			  char cc = ace_lower(complementLetter(dna[am1+slide+2])) ;
			  vtxtPrintf (txt, "VCF %s %d %c %c%s\n", name(seq), a2-slide, cc, cc, bufGS) ;
			}

		      vtxtPrintf (txt, "-D Sliding\n") ;
		      vtxtPrint (txt, "-D Insertion\n") ; /* cleanup */
		      vtxtPrintf (txt,"Multi_insertion %d\n", da) ; /* reinstate */ 
		      
		      if (isDup)
			{
			  vtxtPrintf (txt, "Typ \"Dup%s\"\n", buf) ;
			  if (a1 < a2)
			    {
			      vtxtPrintf (txt, "gName \"%s:g.%d_%ddup%s\"\n", bufSeqTitle, a1 + slide -da + 1, a1 + slide, bufS) ;
			      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddup%s\"\n", name(mrna), am1+slide-da+1, am1+slide, bufS) ;
			    }
			  else
			    {
			      vtxtPrintf (txt, "gName \"%s:g.%d_%ddup%s\"\n", bufSeqTitle, a2 -da + 1, a2, bufG) ;
			      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddup%s\"\n", name(mrna), am1+slide-da+1, am1+slide, bufS) ;
			    }
			  vtxtPrintf (txt, "Duplication\n") ;
			}
		      else
			{
			  vtxtPrintf (txt, "Typ \"Ins%s\"\n", buf) ;
			  if (a1 < a2)
			    {
			      vtxtPrintf (txt, "gName \"%s:g.%d_%dins%s\"\n", bufSeqTitle,a1+slide,a2+slide, bufS) ;
			      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%dins%s\"\n", name(mrna), am1 + slide, am2 + slide, bufS) ;
			    }
			  else
			    {
			      vtxtPrintf (txt, "gName \"%s:g.%d_%dins%s\"\n", bufSeqTitle,a2,a1, bufG) ;
			      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%dins%s\"\n", name(mrna), am1 + slide, am2 + slide, bufS) ;
			    }
			}
		      vtxtPrintf (txt,"Multi_insertion %d %s\n", da, buf) ; /* reinstate */ 
		      if (slide) 
			vtxtPrintf (txt, "Sliding %d\n", slide) ;
		      if (da ==  1) sprintf (bufN, "^") ; 
		      else if (da ==  2) sprintf (bufN, "^^") ; 
		      else if (da ==  3) sprintf (bufN, "^^^") ; 
		      else if (da ==  4) sprintf (bufN, "^^^^") ; 
		      else if (da ==  5) sprintf (bufN, "^^^^^") ; 
		      else if (da ==  6) sprintf (bufN, "^^^^^^") ; 
		      else if (da ==  7) sprintf (bufN, "^^^^^^^") ; 
		      else if (da ==  8) sprintf (bufN, "^^^^^^^^") ; 
		      else if (da ==  9) sprintf (bufN, "^^^^^^^^^") ; 
		      else if (da == 10) sprintf (bufN, "^^^^^^^^^^") ; 
		      else if (da == 11) sprintf (bufN, "^^^^^^^^^^^") ; 
		      else if (da == 12) sprintf (bufN, "^^^^^^^^^^^^") ; 
		      else if (da > 12)  sprintf (bufN, "-%d-", da) ; 
		      
		      ok = TRUE ;
		      
		      for (i = am1 - 21, j = k = 0 ; i < am1 + 30 && i < dnaLn ; i++)
			if (i > 0 && j < 51)
			  {
			    if (i < slide)
			      { bufV[k] = bufR[j] = dna[i] ; RbufV[k] = RbufR[j] = dna[i] ; j++ ; k++ ; }
			    if (i == slide)
			      {
				if (da)
				  {
				    int kk ;
				    for (kk = 0 ; kk < da && kk < 30 ; kk++)
				      { 
					bufR[j] = RbufR[j] = '^'; j++ ;
					bufV[k] = RbufV[k] = ace_upper(buf[kk]) ; k++ ; 
				      }
				  }
			      }
			    if (i >= am1)
			      {
				if (j < 51) { bufR[j] = dna[i] ; RbufR[j] = myRdna(i) ; j++ ; }
				if (k < 51) { bufV[k] = dna[i] ; RbufV[k] = myRdna(i) ; k++ ; }
			      }
			  }
		      RbufR[j] = bufR[j] = 0 ;
		      RbufV[k] = bufV[k] = 0 ;
		    }
		}
	    }
	}
      else if (ac_has_tag (Snp, "DelIns"))
	{
	  int am1 = 0, am2 = 0 ;
	  tbl = ac_tag_table (Snp, "DelIns", h) ;
	  if (tbl)
	    {
	      int dD = ac_table_int (tbl, 0, 0, 0) ;
	      int dI = ac_table_int (tbl, 0, 1, 0) ;
	      const char *ccD = ac_table_printable (tbl, 0, 2, 0) ; 
	      const char *ccI = ac_table_printable (tbl, 0, 3, 0) ; 

	      fs = dI - dD ;

	      vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Bases %d to %d (%d bases) are replaced by %d bases\"\n", name (seq), a1, a2, a1, a2, dD, dI) ;
	      if (tbl->cols >= 4 && ccI)
		{
		  int i, j, k ;
		  ok = TRUE ;
		  if (dD < 20 && dI < 20)
		    {
		      am1 = fromMrna ? m1 : a1 ;
		      vtxtPrintf (txt, "Typ \"DelIns%s>%s\"\n", ccD, ccI) ;
		      vtxtPrintf (txt, "VCF %s %d %c%s %c%s\n", name(seq), a1, ace_lower(dna[am1-1]), ccD, ace_lower(dna[am1-1]), ccI) ;
		      vtxtPrintf (txt, "gName \"%s:g.%d_%ddelins%s\"\n", bufSeqTitle, a1,a2, ccI) ;
		      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddelins%s\"\n", name(mrna), am1+1, am2-1, ccI) ;
		    }
		  else
		    {
		      vtxtPrintf (txt, "Typ \"DelIns_%d_%d\"\n", dD, dI) ;
		      vtxtPrintf (txt, "VCF %d %c%s %c%s\n", name(seq), a1, ace_lower(dna[am1-1]), ccD, ace_lower(dna[am1-1]), ccI) ;
		      vtxtPrintf (txt, "gName \"%s:g.%d_%ddelins(%d,%d)\"\n", bufSeqTitle, a1,a2, dD, dI) ;
		      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddelins(%d,%d)\"\n", name(mrna), am1+1,am2-1, dD, dI) ;
		    }
		  dda = -1 ;
		  am1 = fromMrna ? m1 : a1 ;
		  am2 = fromMrna ? m2 : a2 ;
		  for (i = am1 - 21 + dda, j = k = 0 ; i < am1 + 30 + dda && i < dnaLn ; i++)
		    if (i > 0 && j < 51) 
		      {
			if (i < am1 - 1 || i >= am2-1) 
			  {
			    if (j < 51)  { bufR[j] = dna[i] ; RbufR[j] = myRdna(i) ; j++ ; }
			    if (k < 51)  { bufV[k] = dna[i] ; RbufV[k] = myRdna(i) ; k++ ; }
			  }
			else if (i == am1 - 1) /* copy the 2 buffers */
			  {
			    ccp = ccD ;
			    while (*ccp && j < 51)
			      { bufR[j] = RbufR[j] = ace_upper(*ccp++) ; j++ ; }
			    ccp = ccI ;
			    while (*ccp && k < 51)
			      { bufV[k] = RbufV[k] = ace_upper(*ccp++) ; k++ ; }
			  }
		      }
		  RbufR[j] = bufR[j] = 0 ;
		  RbufV[k] = bufV[k] = 0 ;
		}
	    }
	}
      
      if (! fromMrna)
	{
	  if (bufR[0]) vtxtPrintf (txt, "Reference_genomic_sequence %s\n", bufR) ;
	  if (0 && bufV[0]) vtxtPrintf (txt, "-D Observed_genomic_sequence %s\n", bufV) ; /* old spelling */
	  if (bufV[0]) vtxtPrintf (txt, "Observed__genomic_sequence %s\n", bufV) ;
	}

      vtxtPrintf (txt, "-D Exonic\n") ;
      if (RdnaLn)
	{
	  AC_OBJ Mrna = ac_tag_obj (Snp, "mRNA", h) ;
	  int p1 = 0, p2 = 0 ;
	  KEY product = 0 ;
	  if (Mrna)
	    {
	      AC_TABLE pp = ac_tag_table (Mrna, "Product", h) ;
	      int ir ;
	      BOOL is5 = FALSE , is3 = FALSE ;
	      for (ir = 0 ; pp && ir < pp->rows ; ir++)
		{
		  p1 = ac_table_int (pp, ir, 1 ,0) ;
		  p2 = ac_table_int (pp, ir, 2 ,0) ;
		  if (m1 > p1 && m1 < p2)
		    { product = ac_table_key (pp, ir, 0, 0) ; break ; }
		  if (m1 < p1)
		    is5 = TRUE ;
		  if (m1 > p2)
		    is3 = TRUE ;
		}
	      if (! product && is5)
		vtxtPrintf (txt, "UTR_5Prime\n") ;
	      else if (! product && is3)
		vtxtPrintf (txt, "UTR_3Prime\n") ;
	      if (! pp || ! pp->rows)
		vtxtPrintf (txt, "Non_coding_transcript\n") ;
	    }
	  if (RbufR[0] && RbufV[0])
	    {
	      /* transmit the coord of the conserved hook base written at offset 21 of the snippet */
	      char pR1[60], pV1[60] ;
	      char pR3[60], pV3[60] ;
	      char *cp, *cq ;
	      int frame, k ;
	      vTXT pp = vtxtHandleCreate (h) ;

	      k = strlen (RbufR) ;
	      for (cp = RbufR + k - 1 ; *cp == 'n' && k > 0 ; cp--, k--)
		*cp = 0 ;
	      k = strlen (RbufV) ;
	      for (cp = RbufV + k - 1 ; *cp == 'n' && k > 0 ; cp--, k--)
		*cp = 0 ;

	      k = strlen (RbufR) ;
	      if (k < 52)
		{
		  for (cp = RbufR + k, cq = bufR + k ; *cq && k < 51 ; cp++, cq++, k++)
		    *cp = *cq ;
		  *cp = 0 ;
		}
	      k = strlen (RbufV) ;
	      if (k < 52)
		{
		  for (cp = RbufV + k, cq = bufV + k ; *cq && k < 51 ; cp++, cq++, k++)
		    *cp = *cq ;
		  *cp = 0 ;
		}
	      
	      /* in RNA, use U and possibly lower the common letters */
	      for (cp = RbufR ; *cp ; cp++)
		if (*cp == 'T') *cp = 'U' ;
	      for (cp = RbufV ; *cp ; cp++)
		if (*cp == 'T') *cp = 'U' ;
	      if (0)
		for (cp = RbufR, cq = RbufV ; *cp ; cp++, cq++) 
		  if (ace_lower(*cp) == ace_lower (*cq)) *cp = *cq = ace_lower (*cp) ;
	      
	      vtxtPrintf (txt, "Reference_RNAexon_sequence %s\n", RbufR) ;
	      vtxtPrintf (txt, "Observed__RNAexon_sequence %s\n", RbufV) ;

	      if (product)
		{
		  frame = myTranslate (RbufR, pR1, pR3, m1, p1, h) ;
		  myTranslate (RbufV, pV1, pV3, m1, p1, h) ;
		  
		  for (k = 0, cp = pR1, cq = pV1 ; *cp && *cq ; cp++, cq++)
		    if (ace_lower(*cp) != ace_lower (*cq)) 
		      { *cp = ace_lower (*cp) ; *cq = ace_lower (*cq) ; }
		    else
		      { *cp = ace_upper (*cp) ; *cq = ace_upper (*cq) ; }
		  
		  tsnpSetPName (pp, product, pR1, pV1,pR3, pV3, m1, p1, frame, fs) ;
		  cp = vtxtPtr (pp) ;
		  if (pR1[0]) { char *cq = pR1; while (*cq == 'X') cq++ ; vtxtPrintf (txt,   "Reference_protein_sequence %s %s\n", pR3, cq) ; }
		  if (pV1[0]) { char *cq = pV1; while (*cq == 'X') cq++ ; vtxtPrintf (txt,   "Observed__protein_sequence %s %s %s\n", pV3, cq, cp ? cp : "") ; }
		}
	    }
	}

      if (tsnp->project)
	{
	  vtxtPrintf (txt, "-D Monomodal %s\n", tsnp->project) ;
	  if (tsnpGetMonomodal (tsnp, Snp))
	    vtxtPrintf (txt, "Monomodal %s\n", tsnp->project) ;
	}
      if (0)
	{
	  char PbufR [51] ;
	  char PbufV [51] ;
	  char geneBuf[1000] ; 
	  KEY product = 0, geneBox = 0 ;
	  int r1 = 0 ;

	  memset (PbufR, 0, sizeof (PbufR)) ;
	  memset (PbufV, 0, sizeof (PbufV)) ;
	  memset (geneBuf, 0, sizeof (geneBuf)) ;
	  
	  if (product)
	    {
	      int px1 = 0, pa1 = 0, dx = pa1 % 3 ;
	      sprintf (geneBuf, "%s %s %d", name (geneBox), name (product), px1) ;
	      memcpy (PbufV, RbufV + dx, r1 - dx) ;
	      memcpy (PbufR, RbufR + dx, r1 - dx) ;

	      /* translate PbufV, PbufR */
	    }
	  if (bufR[0]) vtxtPrintf (txt, "Reference_sequence %s %s %s\n", RbufR, PbufR, geneBuf)  ;
	  if (0 && bufV[0]) vtxtPrintf (txt, "-D Observed__sequence \n") ; /* old spelling */
	  if (bufV[0]) vtxtPrintf (txt, "Observed__sequence %s %s %s\n", RbufV, PbufV, geneBuf) ;
	}
    }

  return ok ? 1 : 0 ;
} /* tsnpSetGName */

/*************************************************************************************/
/* enter all types of modifs in a logical order */
static DICT *tsnpMakeVarTypeDict (AC_HANDLE h) 
 {
   DICT *dict = dictHandleCreate (256, h) ;
   char *cp, *cq ;
   /*
   char *Types = "Any,Substitution,Transition,Transversion,Insertion,Deletion,Double insertion,Double deletion,Triple insertion,Triple deletion,Other,A>G,T>C,G>A,C>T,A>T,T>A,G>C,C>G,A>C,T>G,G>T,C>A,Ins A,Ins T,Ins G,Ins C,Del A,Del T,Del G,Del C,Ins AA,Ins TT,Ins GG,Ins CC,Ins AG,Ins CT,Ins AC,Ins GT,Ins TG,Ins CA,Ins TC,Ins GA,Ins AT,Ins TA,Ins GC,Ins CG,Del AA,Del TT,Del GG,Del CC,Del AG,Del CT,Del AC,Del GT,Del TG,Del CA,Del TC,Del GA,Del AT,Del TA,Del GC,Del CG,Ins AAA,Ins TTT,Ins GGG,Ins CCC,Ins AAT,Ins ATT,Ins AAG,Ins CTT,Ins AAC,Ins GTT,Ins TTA,Ins TAA,Ins TTG,Ins CAA,Ins TTC,Ins GAA,Ins GGA,Ins TCC,Ins GGT,Ins ACC,Ins GGC,Ins GCC,Ins CCA,Ins TGG,Ins CCT,Ins AGG,Ins CCG,Ins CGG,Ins ATA,Ins TAT,Ins ATG,Ins CAT,Ins ATC,Ins GAT,Ins AGA,Ins TSNP,Ins AGT,Ins ACT,Ins AGC,Ins GCT,Ins ACA,Ins TGT,Ins ACG,Ins CGT,Ins TAG,Ins CTA,Ins TAC,Ins GTA,Ins TGA,Ins TCA,Ins TGC,Ins GCA,Ins TCG,Ins CGA,Ins GAG,Ins CTC,Ins GAC,Ins GTC,Ins GTG,Ins CAC,Ins GCG,Ins CGC,Ins CAG,Ins CTG,Del AAA,Del TTT,Del GGG,Del CCC,Del AAT,Del ATT,Del AAG,Del CTT,Del AAC,Del GTT,Del TTA,Del TAA,Del TTG,Del CAA,Del TTC,Del GAA,Del GGA,Del TCC,Del GGT,Del ACC,Del GGC,Del GCC,Del CCA,Del TGG,Del CCT,Del AGG,Del CCG,Del CGG,Del ATA,Del TAT,Del ATG,Del CAT,Del ATC,Del GAT,Del AGA,Del TSNP,Del AGT,Del ACT,Del AGC,Del GCT,Del ACA,Del TGT,Del ACG,Del CGT,Del TAG,Del CTA,Del TAC,Del GTA,Del TGA,Del TCA,Del TGC,Del GCA,Del TCG,Del CGA,Del GAG,Del CTC,Del GAC,Del GTC,Del GTG,Del CAC,Del GCG,Del CGC,Del CAG,Del CTG,Ambiguous" ; 
   */
   char *Types = strnew ("Substitution,Transition,Transversion,Deletion,Single_deletion,Double_deletion,Triple_deletion,Deletion_4_30,Long_deletion,Insertion,Single_insertion,Double_insertion,Triple_insertion,Insertion_4_30,Long_insertion,a>g,t>c,g>a,c>t,a>t,t>a,g>c,c>g,a>c,t>g,g>t,c>a,+a,*+a,+t,*+t,+g,*+g,+c,*+c,-a,*-a,-t,*-t,-g,*-g,-c,*-c,++aa,++tt,++gg,++cc,++ag,++ct,++ac,++gt,++tg,++ca,++tc,++ga,++at,++ta,++gc,++cg,*++aa,*++tt,*++gg,*++cc,*++ag,*++ct,*++ac,*++gt,*++tg,*++ca,*++tc,*++ga,*++at,*++ta,*++gc,*++cg,--aa,--tt,--gg,--cc,--ag,--ct,--ac,--gt,--tg,--ca,--tc,--ga,--at,--ta,--gc,--cg,*--aa,*--tt,*--gg,*--cc,*--ag,*--ct,*--ac,*--gt,*--tg,*--ca,*--tc,*--ga,*--at,*--ta,*--gc,*--cg,+++aaa,+++ttt,+++ggg,+++ccc,+++aat,+++att,+++aag,+++ctt,+++aac,+++gtt,+++tta,+++taa,+++ttg,+++caa,+++ttc,+++gaa,+++gga,+++tcc,+++ggt,+++acc,+++ggc,+++gcc,+++cca,+++tgg,+++cct,+++agg,+++ccg,+++cgg,+++ata,+++tat,+++atg,+++cat,+++atc,+++gat,+++aga,+++tsnp,+++agt,+++act,+++agc,+++gct,+++aca,+++tgt,+++acg,+++cgt,+++tag,+++cta,+++tac,+++gta,+++tga,+++tca,+++tgc,+++gca,+++tcg,+++cga,+++gag,+++ctc,+++gac,+++gtc,+++gtg,+++cac,+++gcg,+++cgc,+++cag,+++ctg,*+++aaa,*+++ttt,*+++ggg,*+++ccc,*+++aat,*+++att,*+++aag,*+++ctt,*+++aac,*+++gtt,*+++tta,*+++taa,*+++ttg,*+++caa,*+++ttc,*+++gaa,*+++gga,*+++tcc,*+++ggt,*+++acc,*+++ggc,*+++gcc,*+++cca,*+++tgg,*+++cct,*+++agg,*+++ccg,*+++cgg,*+++ata,*+++tat,*+++atg,*+++cat,*+++atc,*+++gat,*+++aga,*+++tsnp,*+++agt,*+++act,*+++agc,*+++gct,*+++aca,*+++tgt,*+++acg,*+++cgt,*+++tag,*+++cta,*+++tac,*+++gta,*+++tga,*+++tca,*+++tgc,*+++gca,*+++tcg,*+++cga,*+++gag,*+++ctc,*+++gac,*+++gtc,*+++gtg,*+++cac,*+++gcg,*+++cgc,*+++cag,*+++ctg,---aaa,---ttt,---ggg,---ccc,---aat,---att,---aag,---ctt,---aac,---gtt,---tta,---taa,---ttg,---caa,---ttc,---gaa,---gga,---tcc,---ggt,---acc,---ggc,---gcc,---cca,---tgg,---cct,---agg,---ccg,---cgg,---ata,---tat,---atg,---cat,---atc,---gat,---aga,---tsnp,---agt,---act,---agc,---gct,---aca,---tgt,---acg,---cgt,---tag,---cta,---tac,---gta,---tga,---tca,---tgc,---gca,---tcg,---cga,---gag,---ctc,---gac,---gtc,---gtg,---cac,---gcg,---cgc,---cag,---ctg,*---aaa,*---ttt,*---ggg,*---ccc,*---aat,*---att,*---aag,*---ctt,*---aac,*---gtt,*---tta,*---taa,*---ttg,*---caa,*---ttc,*---gaa,*---gga,*---tcc,*---ggt,*---acc,*---ggc,*---gcc,*---cca,*---tgg,*---cct,*---agg,*---ccg,*---cgg,*---ata,*---tat,*---atg,*---cat,*---atc,*---gat,*---aga,*---tsnp,*---agt,*---act,*---agc,*---gct,*---aca,*---tgt,*---acg,*---cgt,*---tag,*---cta,*---tac,*---gta,*---tga,*---tca,*---tgc,*---gca,*---tcg,*---cga,*---gag,*---ctc,*---gac,*---gtc,*---gtg,*---cac,*---gcg,*---cgc,*---cag,*---ctg", 0) ; 

   for (cp = Types, cq = strchr (Types, ',') ; cp ; cp = cq ? cq + 1 : 0, cq = cp ? strchr (cp, ',') : 0) 
     {
       char cc = 0 ;
       if (cq)
	 { 
	   cc = *cq ; 
	   *cq = 0 ; 
	 }
       dictAdd (dict, cp, 0) ;
       if (cq) *cq = cc ;
     }
   ac_free (Types) ;
   return dict ;
 }

/*************************************************************************************/
/* establish VCF name, gName, rNam, pNam */
static void tsnpDbTranslate (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = tsnp->outFileName ? aceOutCreate (tsnp->outFileName, ".translate.ace", tsnp->gzo, h) : 0 ;
  AC_DB db = tsnp->db ;
  AC_ITER iter ;
  int nn = 0, nm = 0, nt = 0 ;
  AC_OBJ Snp = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  
  tsnp->varTypeDict = tsnpMakeVarTypeDict (h) ;

  if (tsnp->select)
    iter = ac_query_iter (db, TRUE, hprintf (h, "find variant IS  \"%s\" ", tsnp->select), 0, h) ;
  else
    iter = ac_query_iter (db, TRUE, "Find Variant", 0, h) ;

  /* NOT DONE: split the 729 multi sub and force create the individual pieces */
  while (ac_free (Snp), Snp = ac_iter_obj (iter))
    {
      nn++ ;
      vtxtClear (txt) ;
      vtxtPrintf (txt, "Variant \"%s\"\n", ac_name (Snp)) ;
      tsnpSetGName (txt, tsnp, Snp, h) ;

      vtxtPrint (txt,"\n") ;
      if (ao)
	aceOutf (ao, "%s\n", vtxtPtr (txt)) ;
      else
	{
	  const char *errors = 0 ;
	  ac_parse (db, vtxtPtr (txt), &errors, 0, h) ; 
	  if (errors && *errors)
	    fprintf(stderr, "tsnpDbTranslate parsing error %s\n", errors) ;
	  else
	    nt++ ;
	}
    }

  fprintf (stderr, "Found %d SNPs, remapped %d, translated %d\n", nn, nm, nt) ;

  ac_free (h) ;
} /* tsnpDbTranslate */
  
/*************************************************************************************/
/*************************************************************************************/
/* count per group */
static void tsnpDbGroupCount (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_HANDLE h1 = 0 ;
  AC_DB db = tsnp->db ;
  AC_ITER iter ;
  ACEOUT ao = tsnp->outFileName ? aceOutCreate (tsnp->outFileName, ".group_count.ace", tsnp->gzo, h) : 0 ;
  int nn = 0, nnn = 0 ;
  AC_OBJ Snp = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  KEYSET gCounts = 0 ;
  const char *errors = 0 ;
  int runClass = 0 ;
  Array r2gs = arrayHandleCreate (200, KEYSET, h) ;
  int minSnpCover = tsnp->minSnpCover ;

  
  char *qq = "select r,g from r in ?runs, g in r >> group where g#SNP_add_counts" ;
  AC_TABLE tbl = ac_bql_table (db, qq, 0, 0, &errors, h) ;
  if (tbl)
    for (int ir = 0 ; ir < tbl->rows ; ir++)
      {
	KEY r = ac_table_key (tbl, ir, 0, 0) ;
	if (r)
	  runClass = class (r) ;
	r = KEYKEY (r) ;
	KEYSET r2g = array (r2gs, r, KEYSET) ;
	if (! r2g)
	  r2g = array (r2gs, r, KEYSET) = keySetHandleCreate (h) ;
	keySet (r2g, keySetMax (r2g)) = KEYKEY (ac_table_key (tbl, ir, 1, 0)) ;
      }
  if (! tbl)
    { ac_free (h) ; return ; }

  if (tsnp->select)
    iter = ac_query_iter (db, TRUE, hprintf (h, "find variant BRS_counts && IS  \"%s\" ", tsnp->select), 0, h) ;
  else
    iter = ac_query_iter (db, TRUE, "Find Variant BRS_counts", 0, h) ;

  /* NOT DONE: split the 729 multi sub and force create the individual pieces */
  while (ac_free (Snp), ac_free (h1),  Snp = ac_iter_obj (iter))
    {
      AC_TABLE tbl = 0 ;
      h1 = ac_new_handle () ;
      nn++ ;
      vtxtPrintf (txt, "\nVariant \"%s\"\n", ac_name (Snp)) ;

      gCounts = keySetHandleCreate (h1) ;
      tbl = ac_tag_table (Snp, "BRS_counts", h1) ;
      if (tbl)
	{
	  for (int ir = 0 ; ir < tbl->rows ; ir++)
	    {
	      KEY run = ac_table_key (tbl, ir, 0, 0) ;
	      KEYSET r2g = array (r2gs, KEYKEY (run), KEYSET) ;
	      if (r2g)
		for (int ig = 0 ; ig < keySetMax (r2g) ; ig++)
		  {
		    KEY g = keySet (r2g, ig) ;
		    for (int ic = 1 ; ic <= 7 ; ic++)
		      {
			int x = ac_table_int (tbl, ir, ic, 0) ;
			keySet (gCounts, 8*g + ic) += x ;
		      }
		  }
	    }
	}
      for (int ig = 0 ; ig < keySetMax (gCounts) ; ig += 8)
	{
	  int m ; /* cover , mutant */
	  float f = 0 ;
	  int c = keySet (gCounts, ig + 1) ;
	  KEY gr = KEYMAKE (runClass, ig/8) ;
	  
	  if (c < minSnpCover)
	    continue ;

	  vtxtPrintf (txt, "BRS_counts %s", name (gr)) ;
	  for (int ic = 1 ; ic <= 7 ; ic++)
	    vtxtPrintf (txt, " %d ", keySet (gCounts, ig + ic)) ;
	  m = keySet (gCounts, ig + 2) ;
	  f = c ? 100.0 * m/c : 0 ;
	  vtxtPrintf (txt, " Frequency %.2f\n", f) ;
	  nnn++ ;
	}
      ac_free (h1) ;
    }
  vtxtPrint (txt,"\n") ;
  if (ao)
    aceOutf (ao, "%s\n", vtxtPtr (txt)) ;
  else
    {
      ac_parse (db, vtxtPtr (txt), &errors, 0, h) ; 
      if (errors && *errors)
	fprintf(stderr, "tsnpDbTranslate parsing error %s\n", errors) ;
    }

  fprintf (stderr, "Found %d SNPs, exported %d counts\n", nn, nnn) ;

  ac_free (h) ;
} /* tsnpDbGroupCount */
  
/*************************************************************************************/
/*************************************************************************************/
/* look for snp in geneBox but not in transcript */
static int tsnpPotential_splice_disruption (TSNP *tsnp, ACEOUT ao)
{
  AC_HANDLE  h1 = 0, h2 = 0, h = ac_new_handle () ;
  AC_DB db = tsnp->db ;
  AC_ITER iter ;
  AC_OBJ variant = 0 ;
  AC_TABLE spl, iMap, viMap = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  int nn = 0, vPos ;

  iter = ac_query_iter (db, TRUE, "find variant geneBox && !mRNA", 0, h) ;
		
  while (ac_free (variant), ac_free (h1), vtxtClear (txt), variant = ac_iter_obj (iter))
    {
      AC_OBJ mrna = 0 ;
      BOOL done = FALSE ;
      AC_ITER iter2 = ac_objquery_iter (variant, ">geneBox ; >mrna COUNT Splicing > 2", h1) ;
      
      h1 = ac_new_handle () ; 
      viMap = ac_tag_table (variant, "IntMap", h1) ;
      vPos =  viMap ? ac_table_int (viMap, 0,1,0) : 0 ;
      
      
      while (ac_free (mrna), ac_free (h2), (! done) && (mrna = ac_iter_obj (iter2)))
	{
	  int jr, g1, g2, gPos ;
	  h2 = ac_new_handle () ; 
	  
	  spl = ac_tag_table (mrna, "Splicing", h2) ;
	  iMap = ac_tag_table (mrna, "IntMap", h2) ;
	  g1 = iMap ? ac_table_int (iMap, 0, 1, 0) : 0 ;
	  g2 = iMap ? ac_table_int (iMap, 0, 2, 0) : 0 ;
	  
	  if (g1 < g2) gPos = vPos - g1 + 1 ;
	  else  gPos = g1 - vPos + 1 ;
	  
	  for (jr = 1 ; g2 && jr < spl->rows - 1 ; jr++)
	    {
	      int a1 = ac_table_int (spl, jr, 0, 0) ;
	      int a2 = ac_table_int (spl, jr, 1, 0) ;
	      int y1 = ac_table_int (spl, jr, 2, 0) ;
	      int y2 = ac_table_int (spl, jr, 3, 0) ;
	      int da ;
	      
	      if (! strcasestr (ac_table_printable (spl, jr, 4, "toto"), "intron"))
		continue ;
	      da = gPos - a1 + 1 ; /* a1 is first base of intron */
	      if (da <= 0) da-- ; /* Plato */
	      switch (da)
		{
		case -16:
		case -1:
		case 1:
		case 2:
		case 3:
		case 5:
		  done = TRUE ;
		  nn++ ;
		  vtxtPrintf(txt, "Variant %s\nNear_donor %d\n\n"
			     , ac_protect (ac_name (variant), h2)
			     , da 
			     , ac_protect (ac_name (mrna), h2)
			     , a1 - 1, y1 /* last base of donor exon */
			     ) ;
		  break ;
		}
	      
	      da = gPos - a2  ; /* a2 is last base of intron */
	      if (da <= 0) da-- ; /* Plato */
	      switch (da)
		{
		case -2:
		case -1:
		  done = TRUE ;
		  nn++ ;
		  vtxtPrintf(txt, "Variant %s\nNear_acceptor %d %s %d %d\n\n"
			     , ac_protect (ac_name (variant), h2)
			     , da
			     , ac_protect (ac_name (mrna), h2)
			     , a2 + 1, y2 /* first base of acceptor exon */
			     ) ;
		  break ;
		}		      
	    }
	}
      /* edit this variant */
      vtxtPrintf (txt, "\n") ;
      if (ao)
	aceOutf (ao, "%s\n", vtxtPtr (txt)) ;
      else
	{
	  const char *errors = 0 ;
	  ac_parse (db, vtxtPtr (txt), &errors, 0, h1) ; 
	  if (errors && *errors)
	    messerror (errors) ;
	}
    }

  
  ac_free (h) ;
  return nn ;
} /* tsnpPotential_splice_disruption */

/***************/
/* probably obsolete, to be transferred into setPname */
static int tsnpCodingModif  (TSNP *tsnp)
{
  AC_HANDLE  h1 = 0, h = ac_new_handle () ;
  ACEOUT ao = tsnp->outFileName ? aceOutCreate (tsnp->outFileName, ".coding.ace", tsnp->gzo, h) : 0 ;
  int nn = 0, nnnn = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  vTXT gTxt = vtxtHandleCreate (h) ;
  vTXT cTxt = vtxtHandleCreate (h) ;
  vTXT pTxt = vtxtHandleCreate (h) ;
  vTXT qq = vtxtHandleCreate (h) ;
  vTXT geneboxes = vtxtHandleCreate (h) ;
  vTXT avgeneboxes = vtxtHandleCreate (h) ;
  /*   vTXT location = vtxtHandleCreate (h) ; */
  vTXT gsnippet = vtxtHandleCreate (h) ;
  vTXT snippet = vtxtHandleCreate (h) ;
  vTXT pSnippet = vtxtHandleCreate (h) ; 
  vTXT pType = vtxtHandleCreate (h) ;
  AC_DB db = tsnp->db ;
  AC_ITER iter ;
  AC_TABLE mrnas ;
  AC_OBJ variant = 0 ;
  const char *errors = 0 ;
  int ir ;
  BOOL debug = FALSE ;
  TSNP* tsnp0 = tsnp ;

  /*
    char *best = "where exists_tag p->best_product or exists_tag p->very_good_product" ;
    ir = ac_keyset_count (ac_dbquery_keyset (db, "query find product best_product", h)) ;
    if (ir < 2)
    best = "" ;

    AND IS \"MEX3B.aAug10:1726:InsGTG\"  
    \"TFAP2B.cAug10:99:A2G\" \"SYNGR2.fAug10-unspliced:98:C2G\"
  */

  if (1)
    iter = ac_query_iter (db, TRUE, "find variant mRNA && IS *  ", 0, h) ;
  else
    iter = ac_query_iter (db, TRUE, "find variant IS  \"NACA.aAug10:6048:Sub:T:C\" ", 0, h) ;

  while (ac_free (variant), ac_free (h1), vtxtClear (txt), variant = ac_iter_obj (iter))
    {
      BOOL is_Potential_splice_disruption = FALSE ;

      if (tsnp != tsnp0)
	invokeDebugger () ;
      nnnn++ ;
      if (debug && nnnn > 20)
	break ;
      h1 = ac_new_handle () ; 
      /* try to remap the variant in the mrna */
      vtxtPrintf (txt, "Variant %s\n-D Potential_splice_disruption\n", ac_protect (ac_name (variant), h1)) ;
      vtxtPrintf (txt, "Non_coding_transcript\n") ; /* default, will possible be overwritten since the schema is UNIQUE Coding/non_coding */
 
      /* export the nature of the Variant in more polite form + translation of proteins */
   
      vtxtClear (qq) ;
      vtxtPrintf (qq, "Colonne 1\n"
		  " Subtitle Column #1\n"
		  " Width 12\n"
		  " Optional\n"
		  " Visible\n"
		  " Class\n"
		  " Class Variant\n"
		  " Condition IS \"%s\"\n"
		  " \n"
		  " Colonne 2\n"
		  " Subtitle Column #2 mRNA\n"
		  " Width 12\n"
		  " Optional\n"
		  " Visible\n"
		  " Class\n"
		  " Class mRNA\n"
		  " From 1\n"
		  " Tag mRNA\n"
		  " \n"
		  " Colonne 3\n"
		  " Subtitle Column #3 x1\n"
		  " Width 12\n"
		  " Optional\n"
		  " Visible\n"
		  " Integer\n"
		  " Right_of 2\n"
		  " Tag  HERE \n"
		  " \n"
		  " Colonne 4\n"
		  " Subtitle Column #4 x2\n"
		  " Width 12\n"
		  " Mandatory\n"
		  " Visible\n"
		  " Integer\n"
		  " Right_of 3\n"
		  " Tag  HERE \n"
		  " \n"
		  " Colonne 5\n"
		  " Subtitle Column #5 Product\n"
		  " Width 12\n"
		  " Mandatory\n"
		  " Visible\n"
		  " Class\n"
		  " Class Product\n"
		  " From 2\n"
		  " Tag Product\n"
		  " Condition best_product OR very_good_product\n"
		  " \n"
		  " Colonne 6\n"
		  " Subtitle Column #6 p1\n"
		  " Width 12\n"
		  " Mandatory\n"
		  " Visible\n"
		  " Integer\n"
		  " Right_of 5\n"
		  " Tag  HERE \n"
		  " \n"
		  " Colonne 7\n"
		  " Subtitle Column #7 p2\n"
		  " Width 12\n"
		  " Mandatory\n"
		  " Visible\n"
		  " Integer\n"
		  " Right_of 6\n"
		  " Tag  HERE \n"
		  " \n"
		  " Colonne 8\n"
		  " Subtitle Column #8 mRNA sequence\n"
		  " Width 12\n"
		  " Mandatory\n"
		  " Visible\n"
		  " Class\n"
		  " Class DNA\n"
		  " From 2\n"
		  " Tag DNA\n"
		  " \n"
		  " Colonne 9\n"
		  " Subtitle Column #9 Variant typ\n"
		  " Width 12\n"
		  " Optional\n"
		  " Visible\n"
		  " Text\n"
		  " From 1\n"
		  " Tag Typ\n"
		  " \n"
		  " Colonne 10\n"
		  " Subtitle Column #10 chrom\n"
		  " Width 12\n"
		  " Optional\n"
		  " Visible\n"
		  " Class\n"
		  " Class Map\n"
		  " From 1\n"
		  " Tag IntMap\n"
		  " \n"
		  " Colonne 11\n"
		  " Subtitle Column #11 a1\n"
		  " Width 12\n"
		  " Optional\n"
		  " Visible\n"
		  " Integer\n"
		  " Right_of 10\n"
		  " Tag  HERE \n"
		  " \n"
		  " Colonne 12\n"
		  " Subtitle Column #12 a2\n"
		  " Width 12\n"
		  " Optional\n"
		  " Visible\n"
		  " Integer\n"
		  " Right_of 11\n"
		  " Tag  HERE \n"
		  " \n"
		  ,  ac_name(variant)
		  ) ; 
      mrnas = ac_tablemaker_table (db, vtxtPtr (qq), 0, ac_tablemaker_text , 0 , 0, &errors, h) ;


      if (! mrnas || ! mrnas->rows)
	{
	  /*
	    mrnas = ac_bql_table (db, messprintf("select v,m,pos,s, dna,typ from v in class \"Variant\" where v == \"%s\", m in v->mrna,pos in m[1],s in pos[1], dna in m->DNA where dna, typ in v->Typ", ac_name(variant),best), 0, 0, 0, h1) ;
	  */
     vtxtClear (qq) ;
      vtxtPrintf (qq, "Colonne 1\n"
		  " Subtitle Column #1\n"
		  " Width 12\n"
		  " Optional\n"
		  " Visible\n"
		  " Class\n"
		  " Class Variant\n"
		  " Condition IS \"%s\"\n"
		  " \n"
		  " Colonne 2\n"
		  " Subtitle Column #2 mRNA\n"
		  " Width 12\n"
		  " Optional\n"
		  " Visible\n"
		  " Class\n"
		  " Class mRNA\n"
		  " From 1\n"
		  " Tag mRNA\n"
		  " \n"
		  " Colonne 3\n"
		  " Subtitle Column #3 x1\n"
		  " Width 12\n"
		  " Optional\n"
		  " Visible\n"
		  " Integer\n"
		  " Right_of 2\n"
		  " Tag  HERE \n"
		  " \n"
		  " Colonne 4\n"
		  " Subtitle Column #4 x2\n"
		  " Width 12\n"
		  " Mandatory\n"
		  " Visible\n"
		  " Integer\n"
		  " Right_of 3\n"
		  " Tag  HERE \n"
		  " \n"
		  " Colonne 5\n"
		  " Subtitle Column #5 mRNA sequence\n"
		  " Width 12\n"
		  " Mandatory\n"
		  " Visible\n"
		  " Class\n"
		  " Class DNA\n"
		  " From 2\n"
		  " Tag DNA\n"
		  " \n"
		  " Colonne 6\n"
		  " Subtitle Column #6 Variant Typ\n"
		  " Width 12\n"
		  " Optional\n"
		  " Visible\n"
		  " Text\n"
		  " From 1\n"
		  " Tag Typ\n"
		  " \n"
		  " Colonne 7\n"
		  " Subtitle Column #7 Chrom\n"
		  " Width 12\n"
		  " Optional\n"
		  " Visible\n"
		  " Class\n"
		  " Class Map\n"
		  " From 1\n"
		  " Tag IntMap\n"
		  "\n"
		  " Colonne 8\n"
		  " Subtitle Column #8 a1\n"
		  " Width 12\n"
		  " Optional\n"
		  " Visible\n"
		  " Integer\n"
		  " Right_of 7\n"
		  " Tag  HERE \n"
		  "\n"
		  " Colonne 9\n"
		  " Subtitle Column #9 a2\n"
		  " Width 12\n"
		  " Optional\n"
		  " Visible\n"
		  " Integer\n"
		  " Right_of 8\n"
		  " Tag  HERE \n"
		  "\n"
		  ,  ac_name(variant)
		  ) ; 

	  mrnas = ac_tablemaker_table (db, vtxtPtr (qq), 0, ac_tablemaker_text , 0 , 0, &errors, h) ;
	}

      if (mrnas && mrnas->rows)
	{
	  int isFrameshift = 0 ;
	  int hasCDS = 0, iMrnaWithCDS = 0 ;
	  int delta = 20, strand, x1, x2 ;
	  AC_OBJ Mrna = 0 ;
	  const char *typ ;

	  /* try to locate a CDS */
	  for (ir = 0 ; !hasCDS && ir < mrnas->rows ; ir++)
	    if (ac_table_int (mrnas, ir, 2, 0) && ac_table_int (mrnas, ir, 6, 0)) 
	      { 
		int x1 = ac_table_int (mrnas, ir, 2, 0) ;
		int p1 = ac_table_int (mrnas, ir, 5, 0) ;
		int p2 = ac_table_int (mrnas, ir, 6, 0) ;

		if (x1 + 1 >= p1 + 1 && x1 + 1  <= p2)  /* the first edited base touches the CDS + stop */
		  { hasCDS = 5 ; iMrnaWithCDS = ir ; }
		else if (p2 > 0 && x2 - 1 >= p1 + 1 && x2 - 1  <= p2)  /* the last edited base touches the CDS + stop */
		  { hasCDS = 5 ; iMrnaWithCDS = ir ; }
		else if (p2 > 0 && x1 + 1 > p2)  /* the last edited base touches the CDS */
		  { hasCDS = 4 ; iMrnaWithCDS = ir ; }
		else if (hasCDS < 3 && p1 > 0 && x1 + 1 < p2)  /* the last edited base touches the CDS */
		  { hasCDS = 3 ; iMrnaWithCDS = ir ; }
		if (hasCDS == 4)
		  break ;
	      }
	  if (hasCDS == 3) 
	    vtxtPrintf(txt, "UTR_5prime\n") ;
	  else 	if (hasCDS == 4) 
	    vtxtPrintf(txt, "UTR_3prime\n") ;
	  if (hasCDS > 2)
	    { 
	      ir = iMrnaWithCDS ;
	      hasCDS = 2 ;
	    }
	  else
	    {
	      /* try to locate any dna */
	      for (ir = 0 ; !hasCDS && ir < mrnas->rows ; ir++)
		{
		  if (ac_table_key (mrnas, ir, 1, 0))
		    { hasCDS = 1 ;  break ; }
		} 
	    }
	  /* name the gene  and locate a Potential_splice_disruption */
	  if (mrnas && ir < mrnas->rows)
	    {
	      AC_OBJ mrna = ac_table_obj (mrnas, ir, 1, 0) ;
	      const char *ccp = mrna ? ac_tag_printable (mrna, "Gene", 0) : 0 ;
	      if (ccp)
		vtxtPrintf (txt, "Gene \"%s\"\n", ccp) ;
	      if (! is_Potential_splice_disruption)
		{
		  int jr ;
		  AC_TABLE spl = ac_tag_table (mrna, "Splicing", h1) ;
		  AC_TABLE iMap = ac_tag_table (mrna, "IntMap", h1) ;
		  int g1 = iMap ? ac_table_int (iMap, 0, 1, 0) : 0 ;
		  int g2 = iMap ? ac_table_int (iMap, 0, 2, 0) : 0 ;
		  int gPos =  ac_table_int (mrnas, ir, mrnas->cols > 10 ? 12 : 9, 0) ;
		  const char *typ = ac_table_printable (mrnas, ir, mrnas->cols > 10 ? 8 : 5 , "xxx") ;
		  

		  if (spl && gPos && typ)
		    {
		      BOOL isSub = ! strcmp (typ, "Substitution") ;

		      if (g1 < g2) gPos = gPos - g1 + 1 ; /* last exact base before error in pre-mRNA coords */
		      else  gPos = g1 - gPos + 1 ;
		      
		      for (jr = 1 ; g2 && jr < spl->rows - 1 ; jr++)
			{
			  int a1 = ac_table_int (spl, jr, 0, 0) ;
			  int a2 = ac_table_int (spl, jr, 1, 0) ;
			  int y1 = ac_table_int (spl, jr, 2, 0) ;
			  int y2 = ac_table_int (spl, jr, 3, 0) ;
			  int da ;
			  
			  if (! strcasestr (ac_table_printable (spl, jr, 4, "toto"), "intron"))
			    continue ;
			  da = gPos + 1 - a1 + 1 ; /* a1 is first base of intron */
			  if (da <= 0) da-- ; /* Plato */
			  if (isSub)
			    {
			      switch (da)
				{
				case -16:
				case -1:
				case 1:
				case 2:
				case 3:
				case 5:
				  vtxtPrintf(txt, "Near_donor %d %s %d %d\n"
					     , da
					     , ac_protect (ac_name (mrna), h1)
					     , a1 - 1, y1 /* last base of donor exon */
					     ) ;
			      break ;
				}
			    }
			  da = gPos + 1 - a2  ; /* a2 is last base of intron */
			  if (da <= 0) da-- ; /* Plato */
			  if (isSub) 
			    {
			      switch (da)
				{
				case -2:
				case -1:
				  vtxtPrintf(txt, "Near_acceptor %d\n"
					     , da
					     , ac_protect (ac_name (mrna), h1)
					     , a2 + 1, y2 /* first base of acceptor exon */
					     ) ;
				  break ;
				}		 
			    }     
			}
		      
		      is_Potential_splice_disruption = TRUE ;
		    }
		}
	      ac_free (mrna) ;	      
	    }
	  if (hasCDS) /* export a dna/peptide motif of a few letters */
	    {
	      int pos, u1, u2, delta1, delta2, i, j, ln ;
	      const char *dna, *cd1, *cd2 = 0, *cd22 = 0 ;
	      char cc, *cd3 ;
	      char buf1[1024], buf2[1024] ;
	      char tbuf1[1024], tbuf2[1024] ;
	      
	      if (debug) fprintf (stderr, "%s -> %s ->hasCDS=%d\n", ac_name(variant), ac_table_printable(mrnas, ir, 1, "toto"), hasCDS);
	      pos = ac_table_int (mrnas, ir, 2, 0) ;
	      strand = ac_table_int (mrnas, ir, 3, 1) ;
	      Mrna = ac_table_obj (mrnas, ir, 1, h1) ;
	      dna = ac_obj_dna (Mrna, h1) ;
	      typ = ac_table_printable (mrnas, ir, mrnas->cols > 10 ? 8 : 5 , "xxx") ;
	      x1 = ac_table_int (mrnas, ir, 5, 0) ;
	      x2 = ac_table_int (mrnas, ir, 6, 0) ;

	      u1 = u2 = pos ;
	      if (typ && (! strncasecmp(typ,"Del",3) || ! strncasecmp(typ,"Ins",3)))
		{
		  int k, k0 = strlen(typ) - 3 ; /* number of base inserted or deleted */

		  /* shift left on homopolymer */
		  for (u1 = pos - 1, cd1 = dna + u1 - 1, k = k0 - 1 ; u1 > 1 && ace_lower(typ[3+k]) == *cd1 ; u1--, cd1--, k = (k0 + k - 1) % k0) ;  
		  if (u1 >= 1 && u1 < pos) { u1++ ; cd1++ ;}
		  /* shift right on homopolymer
		  for (u2 = pos+1, cd2 = dna + u2 - 1, k = 0 ; ace_lower(typ[3+k]) == *cd2 ; u2++, cd2++, k = (k0 + k - 1) % k0) ; 
		  u2-- ; cd2-- ;
		  */
		}
	      /* extend by delta if possible */
	      ln = strlen (dna) ;
	      delta1 = u1 - 1 ; if (delta1 > delta) delta1 = delta ;
	      delta2 = ln - u2 ; if (delta2 > delta+3) delta2 = delta+3 ;
	      /* fall in frame */
	      if (hasCDS == 2 && pos >= 0*x1 && pos <= x2)
		{ 
		  int k ;
		
		  i = u1 - delta1 ;
		  k = (30000000 + i - x1) % 3 ;
		  delta1 += k ;
		  if (u1 - delta1 < 1) delta1 -= 3 ;
		}
	      /* copy in the buffer the delta1 extension */
	      for (i = u1 - delta1, j = 0, cd1 = dna + i - 1, cd3 = buf1 ; j < delta1 ; j++, i++, cd1++, cd3++)
		*cd3 = ace_lower(*cd1) ;
	      /* copy in the buffer the homopolymer */
	      for ( ; i < u2 ; j++, i++, cd1++, cd3++)
		*cd3 = ace_upper(*cd1) ;
	      if ( strncasecmp (typ, "Ins", 3))
		for ( ; i <= u2 ; j++, i++, cd1++, cd3++)
		  *cd3 = ace_upper(*cd1) ;
	      /* in case of multi deletion, write in capital letters */
	      if (! strncasecmp (typ, "del", 3))
		{
		  int k ;

		  k = pos - u1 + 1 ;
		  cd1 -= k ; cd3 -= k ; i -= k ; j -= k ;
		  for (k = pos - u1 + strlen(typ) - 3 ; k > 0 ; k--, j++, i++, cd1++, cd3++)
		    *cd3 = ace_upper(*cd1) ;
		}
	      /* copy in the buffer the delta2 extension */
	      for ( ; i <= u2 + delta2 && j < 1023 ; j++, i++, cd1++, cd3++)
		*cd3 = ace_lower(*cd1) ;
	      while (((cd3 - buf1)%3) && *cd1) /* complete to triplets */
		*cd3++ = ace_lower(*cd1++) ;
	      *cd3 = 0 ; /* zero terminate */

	      /* this time we export the modifed buffer */
	      /* copy in the buffer the delta1 extension */
	      for (i = u1 - delta1, j = 0, cd1 = dna + i - 1, cd3 = buf2 ; j < delta1 ; j++, i++, cd1++, cd3++)
		*cd3 = ace_lower(*cd1) ;
	      /* copy in the buffer the homopolymer */
	      for ( ; i < u2 ; j++, i++, cd1++, cd3++)
		*cd3 = ace_upper(*cd1) ;
	      if ( strncasecmp (typ, "Ins", 3))
		for ( ; i <= u2 ; j++, i++, cd1++, cd3++)
		  *cd3 = ace_upper(*cd1) ;
	      /* copy the mutation */
	      if (typ[1] == '2') /* substitution */
		{
		  cc =  typ[2] ;
		  if (strand == -1)
		    cc = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)cc])] ;
		  *(cd3 - 1) = ace_upper (cc) ; /* we already copied the original letter */
		}
	      else if (! strncasecmp (typ, "Ins", 3)) /* insert the insertion */
		{ 
		  if (1) *cd3++ = ace_upper(*cd1++) ; /* because the variant name reports a la VCF the base before the insertion */
		  if (strand > 0)
		    {
		      for (cd2 = typ+3; *cd2 ; cd2++, cd3++, j++) 
			{ 
			  char cc = *cd2 ;
			  *cd3 = ace_upper(cc) ; isFrameshift++ ; 
			  
			  /* specialize Ins? to Dup? if *cp3 is a repeated letter */
			  if (isFrameshift == 1 && cd2[1] == 0 &&
			      *cd3 == ace_upper(*(cd3 - 1)))
			    vtxtPrintf (txt, "Dup%c\n", *cd3) ;
			}
		    }
		  else
		    {
		      /* we must insert after the base reporting the insert */
		      *cd3++ = *cd1++ ; i++ ;
		      for (cd22 = typ+3; *cd22 ; cd22++) ;
		      for (--cd22 ; cd22 >= typ+3 ; cd22--, cd3++, j++)
			{ 
			  cc = *cd22 ;
			  if (strand == -1)
			    cc = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)cc])] ;
			  *cd3 = ace_upper(cc) ; isFrameshift++ ; 

			  /* specialize Ins? to Dup? if *cp3 is a repeated letter */
			  if (isFrameshift == 1 && typ[4] == 0 &&
			      *cd3 == ace_upper(*(cd3 - 1)))
			    vtxtPrintf (txt, "Dup%c\n", ace_upper(*cd3)) ;
			}
		    }
		}
	      else if (! strncasecmp (typ, "del", 3)) /* delete the deletion */
		{ 
		  if (1) *cd3++ = *cd1++ ; /* because the variant name reports a la VCF the base before the deletion */
		  cd3-- ; cd1-- ; 
		  for (cd2 = typ+3; *cd2 ; cd2++) 
		    { 
		      /* specialize Del? to Dim? if *cp3 is a repeated letter */
		      if (isFrameshift == 0 && cd1[0] == cd1[-1] && typ[4] == 0)
			vtxtPrintf (txt, "Dim%c\n", ace_upper(typ[3])) ;
		      cd1++ ; isFrameshift-- ; 
		    }
		}
	      else if (strchr (typ, '2')) /* multi_substitution */
		{
		  char *buf = strnew (typ, h1) ;
		  char *cp = strstr (buf, "2") ;
		  
		  *cp++ = 0 ;
		  if (strand == 1)
		    {
		      cd3-- ;
		      while (*cp)
			*cd3++ = ace_upper (*cp++) ;
		    }
		  else
		    {
		      char *cq = cp + strlen (cp) - 1 ;
		      cd3-- ;
		      while (cq >= cp)
			*cd3++ = ace_upper (dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)(*cq--)])])   ;
		    }
		}
	      
	      /* copy in the buffer the delta2 extension */
	      for ( ; i <= u2 + delta2 && j < 1023 ; j++, i++, cd1++, cd3++)
		*cd3 = ace_lower(*cd1) ;
	      while (((cd3 - buf2)%3) && *cd1) /* complete to triplets */
		*cd3++ = ace_lower(*cd1++) ;
	      *cd3 = 0 ; /* zero terminate */
	      if (0 && isFrameshift < 0)
		{ i = - isFrameshift / 3 ; *(cd3 - 3*i) = 0 ; }

#ifdef JUNK
	      /* adjust upper/lower */
	      buf = strnew (typ, h1) ;
	      tsnpIsSliding (buf, buf1, buf2, TRUE, TRUE, 0, 0) ;
	      /* export the buffer */
	      vtxtPrintf(txt, "Reference_sequence  %s\n", buf1) ;
	      vtxtPrintf(txt, "Observed__sequence   %s\n", buf2) ;
#endif
	      /* translate */
	      /* see www.hgvs.org den Dunnen and Antonarakis nov 15, 2014  HGV human genone variations */
	      if (hasCDS == 2 && pos >= 0 * x1 && pos <= x2)
		{
		  Array aa = arrayHandleCreate (64, char , h1) ;
		  char * translationTable = pepGetTranslationTable(ac_table_key (mrnas, ir, 0, 0), 0) ;
		  char *cp1, *cp2 ;
		  int dMet = 0, dStop = 0, stop1 = -1, stop2 = -1, ppos ;
		  int cc10, cc11, cc12, cc21, cc22, cc99 = 0 ;
		  int delta ;

		  i = strlen(buf1) ;
		  array(aa, i, char ) = 0 ;
		  arrayMax(aa) = i ;
		  cp1 = arrp (aa, 0, char) ; cp2 = buf1 ;
		  memcpy (cp1, cp2, i+1) ;
		  dnaEncodeArray (aa) ;
		  memset (tbuf1, 0, sizeof(tbuf1)) ;
		  memset (tbuf2, 0, sizeof(tbuf2)) ;
		  tbuf1[0] = 0 ; cp2 = tbuf1 ; stop1 = 0 ;
		  if (arrayMax (aa) > 2)
		    for (cp1 = arrp (aa, 0, char), i = 0  ; i < arrayMax (aa) - 2 ; i+= 3 , cp1 += 3 , cp2++)
		      {
			*cp2 = e_codon (cp1, translationTable) ;
			if (*cp2 == '*')
			  { *cp2 = 'X' ; stop1 = cp2 - tbuf1 ;  if (3*stop1 > delta1 + 3) { cp2++ ; break ;} } 
		      }
		  *cp2 = 0 ;

		  /* now the modified translation */
		  i = strlen(buf2) ;
		  array(aa, i, char) = 0 ;
		  arrayMax(aa) = i ;
		  cp1 = arrp (aa, 0, char) ; cp2 = buf2 ;
		  memcpy (cp1, cp2, i+1) ;
		  dnaEncodeArray (aa) ;
		  tbuf2[0] = 0 ; cp2 = tbuf2 ; stop2 = 0 ;
		  if (arrayMax (aa) > 2)
		    for (cp1 = arrp (aa, 0, char), i = 0
			   ; i < arrayMax (aa) - 2  ; i+= 3 , cp1 += 3 , cp2++)
		      {
			*cp2 = e_codon (cp1, translationTable) ;
			if (*cp2 == '*')
			  { *cp2 = 'X' ; stop2 = cp2 - tbuf2 ;  if (3*stop2> delta1 + 3) { cp2++; break ;}}
		      }
		  *cp2 = 0 ;


		  /* compute now the location of the stops in the complete mRNA and modified mRNA */
		  delta = (pos - x1) % 3 ;
		  locateStop (dna, translationTable
			      , pos - delta
			      , pos + 6 - isFrameshift - delta
			      , &stop1, &stop2
			      ) ;
		  if (stop2 >= 0) 
		    stop2 +=  (isFrameshift > 0 ? 2 : 2) ; /* we started on the second or third codon */

		  /* the 2 first codons of the modified sequence must be read off tbuf2, not off dna */
		  cc10 = tbuf1[-1+(pos - u1 + delta1)/3] ;
		  cc11 = tbuf1[(pos - u1 + delta1)/3] ;
		  cc12 = tbuf1[1+(pos - u1 + delta1)/3] ;
		  /*		  cc13 = tbuf1[2+(pos - u1 + delta1)/3] ; */
		  cc21 = tbuf2[(pos - u1 + delta1)/3] ;
		  cc22 = tbuf2[1+(pos - u1 + delta1)/3] ;

		  if (cc21 == 'X') stop2 = 0 ;
		  else if (cc22 == 'X') stop2 = 1 ;		      
		  if (stop1 >= 0 && stop2 >= 0 && stop1 != stop2)
		    dStop =  stop2 - stop1 ;

		  vtxtPrint (txt, "-D Coding\n") ;
		  /* interpret the protein independently of what happens in the dna */
		  /* the 2 first codon are interpreted explicitelly, downstream we only care about stops */

		  /* analyse variations of the first codon */
		  if (hasCDS == 2 && pos < x1 && pos <= x2 && cc21 == 'M' && isFrameshift %3 == 0)
		    {
		      int i, delta, n = strlen (dna) ;
		      const char *ccp ;
		      Array dna2 = 0 ;
		      
		      dMet = 0, cc99 = 0, delta = (x1 - pos) % 3 ; 
		      dna2 = arrayHandleCreate (n + 3, char, 0) ;
		      array (dna2, n - 1, char) = 0 ;
		      memcpy (arrp (dna2, 0, char), dna, n) ;
		      dnaEncodeArray (dna2) ;

		      i = pos + delta - 1 ; ccp = arrp (dna2, pos + delta - 1, char) ;
		      for ( ; i <= x1 ; i+= 3, ccp += 3, dMet++)
			if (e_codon (ccp, translationTable) == '*')
			  break ;
		      if (i > x1)
			{
			  const char *ccp = pepShortName[(int)e_codon(arrp (dna2, x1 + 3 - 1, char), translationTable)] ;
			  dMet = (x1 + 2 - pos) / 3 ; 
			  vtxtPrintf(txt, "Met1_gained Met1_%s2ins", ccp) ;
			  if (dMet > 1 && cc22) vtxtPrintf(txt, "%s", pepShortName[(int) cc22]) ;
			  for (i = pos + (delta ? delta - 3 : 0) - 1 + 6, ccp = arrp (dna2, i, char) ; i < x1 ; i += 3, ccp += 3)
			    {
			      const char *ccq = pepShortName[(int)e_codon (ccp, translationTable)] ;
			      vtxtPrintf(txt, "%s", ccq  && *ccq != 'X' ? ccq : "---") ;
			    }
			  vtxtPrint (txt, "\n") ;
			}
		      else 
			{ vtxtPrintf(txt, "UTR_5prime\n") ; dStop = dMet = 0 ; }
		    }
		  else if (hasCDS == 2 && pos < x1)
		    {
		      vtxtPrintf(txt, "UTR_5prime\n") ; dStop = dMet = 0 ;
		    }
		  else if (cc11 == 'M' && cc21 != 'M' && pos == x1)
		    {
		      /* we lost the Met, we need to find the next one */
		      dMet = dStop = 0 ;
		      if (cc22 == 'M')
			dMet = - 1 ;
		      else
			{
			  /*  locateMet (dna, translationTable, x1, &cc99, &dMet) ; */
			  messcrash ("locateMet modified") ;
			  if (dMet == 2)  cc99 = cc22 ;
			  dMet = - dMet ;
			}
		    }

		  if (pos < x2 && dStop + dMet)
		    vtxtPrintf (txt, "Extension %d AA\n", dStop + dMet) ;
		  ppos = 1 + (pos - x1) / 3 ; /* biological numbering of the codon containing the snp */

		  if (pos >= x1)
		    {
		      if (dMet < 0)
			{
			  if (dMet == 1)
			    vtxtPrintf (txt, "Met1del") ;
			  else if (cc99)
			    vtxtPrintf (txt, "Met1_lost %s2_M%ddel", pepShortName[(int)cc12], dMet) ;
			  else
			    vtxtPrintf (txt, "Met1_lost Met>%s", pepShortName[(int)cc21]) ;
			}
		      else if (
			       (cc11 == 'X' && cc21 == 'X') ||
			       (cc11 == cc21 && cc12 == 'X' && cc22 == 'X') ||
			       (isFrameshift == 0 && cc11 == cc21)
			       )
			vtxtPrintf (txt, "Synonymous %s%d=", cc11 == 'X' ? "Ter" : (cc11 ? pepShortName[(int)cc11] : ""), ppos) ;
		      else if (cc11 != cc21 && cc21 != 'X' && cc11 !='X' &&
			       (
				(cc12 == 'X' && cc22 == 'X') ||
				(isFrameshift == 0 && cc12 == cc22)
				)
			       )
			vtxtPrintf (txt, "AA_substitution %s%d%s"
				    , pepShortName[(int)cc11]
				    , ppos
				    , cc21 ? pepShortName[(int)cc21] : ""
				    ) ;
		      else if (cc11 == cc21 && cc12 != cc22 && stop1 == 2 && stop2 == 2)
			vtxtPrintf (txt, "AA_substitution %s%d%s"
				    , pepShortName[(int)cc12]
				    , 1 + ppos
				    , cc22 ? pepShortName[(int)cc22] : ""
				    ) ;
		      else if (cc11 == 'X')
			{
			  if (isFrameshift == 0)
			    vtxtPrint (txt, "Stop_to_AA") ;
			  else if (isFrameshift % 3 == 0)
			    vtxtPrint (txt, "Frame_preserving_indel ") ; 
			  else 
			    {
			      vtxtPrint (txt, "Frameshift") ;
			      if (pepShortName[(int)cc21])
				{
				  vtxtPrintf (txt, " Ter%dext%s"
					      , ppos
					      , pepShortName[(int)cc21]
					      ) ;
				  if (stop2 >= 0)
				    vtxtPrintf (txt, "Ter%d", stop2) ;
				  else
				    vtxtPrint (txt, "Ter") ; 
				}
			    }
			}
		      else if (cc11 == cc21 && cc12 == 'X')
			{
			  if (isFrameshift == 0)
			    vtxtPrint (txt, "Stop_to_AA") ;
			  else if (isFrameshift % 3 == 0)
			    vtxtPrint (txt, "Frame_preserving_indel ") ; 
			  else 
			    {
			      vtxtPrint (txt, "Frameshift") ;
			      if(cc22)
				{
				  vtxtPrintf (txt, " Ter%dext%s"
					      , 1 + ppos
					      , pepShortName[(int)cc22]
					      ) ;
				  if (stop2 >= 0)
				    vtxtPrintf (txt, "Ter%d", stop2) ;
				  else
				    vtxtPrint (txt, "Ter") ; 
				}
			    }
			}
		      else if (cc21 == 'X')
			{
			  if (isFrameshift == 0)
			    vtxtPrint (txt, "AA_to_stop") ;
			  else if (isFrameshift % 3 == 0)
			    vtxtPrint (txt, "Frame_preserving_indel ") ; 
			  else 
			    {
			      vtxtPrint (txt, "Frameshift") ;
			      vtxtPrintf (txt, " %s%dTer"
					  , pepShortName[(int)cc11]
					  , ppos
					  ) ;
			    }
			}
		      else if (cc11 == cc21 && cc22 == 'X')
			{
			  if (isFrameshift == 0)
			    vtxtPrint (txt, "AA_to_stop") ;
			  else if (isFrameshift % 3 == 0)
			    vtxtPrint (txt, "Frame_preserving_indel ") ; 
			  else 
			    {
			      vtxtPrint (txt, "Frameshift") ;
			      vtxtPrintf (txt, " %s%dTer"
					  , pepShortName[(int)cc12]
					  , 1 + ppos
					  ) ;
			    }
			}
		      else if (isFrameshift == 0) /* but we already know that 2 bases are affected */
			vtxtPrintf (txt, "AA_substitution %s%d_%s%ddelins%s"
				    , pepShortName[(int)cc11], ppos 
				    , pepShortName[(int)cc12], ppos + 1
				    , cc21 ?  pepShortName[(int)cc21] : ""
				    , cc22 ? pepShortName[(int)cc22] : ""
				    ) ;
		      else if (isFrameshift == -3) /*  remove triplet, dim or del */
			{
			  vtxtPrint (txt, "Frame_preserving_indel ") ; 
			  if (cc12 == cc21)
			    {
			      if (cc11 == cc12)  /* dim, but we must shift rigth as much as possible */
				{ 
				  int j = 0 ;
				  while (tbuf2[j+(pos - u1 + delta1)/3] == cc11) j++ ;
				  vtxtPrintf (txt, " %s%ddim", pepShortName[(int)cc12], j + ppos) ;
				}
			      else if (cc10 == cc11)
				vtxtPrintf (txt, " %s%ddim", pepShortName[(int)cc11], ppos) ;
			      else
				vtxtPrintf (txt, " %s%ddel", pepShortName[(int)cc11], ppos) ;
			    }
			  else if (cc11 == cc21)
			    vtxtPrintf (txt, " %s%d%ddel"
					, pepShortName[(int)cc12], ppos + 1
					) ;
			  else
			    vtxtPrintf (txt, " %s%d_%s%ddelins%s"
					, pepShortName[(int)cc11], ppos 
					, pepShortName[(int)cc12], ppos + 1
					, cc21 ? pepShortName[(int)cc21] : ""
					) ;
			}
		      else if (isFrameshift == 3) /* insert triplet, dup or ins */
			{
			  vtxtPrint (txt, "Frame_preserving_indel ") ; 
			  
			  if (cc11 == cc22)
			    {
			      if (cc21 == cc22)  /* dim, but we must shift rigth as much as possible */
				{ 
				  int j = 0 ;
				  while (tbuf1[j+1+(pos - u1 + delta1)/3] == cc11) j++ ;
				  vtxtPrintf (txt, " %s%ddup", pepShortName[(int)cc21], j + ppos) ;
				}
			      else if (cc10 == cc21)
				vtxtPrintf (txt, " %s%ddup", pepShortName[(int)cc10], -1 + ppos) ;
			      else
				vtxtPrintf (txt, " %s%d_%s%dins%s"
					    , pepShortName[(int)cc10], -1 + ppos
					    , pepShortName[(int)cc11], ppos
					    , cc21 ? pepShortName[(int)cc21] : "" 
					    ) ;
			    }
			  else if (cc11 == cc21)
			    vtxtPrintf (txt, "\nAA_substitution  %s%dins%s"
					, pepShortName[(int)cc11], ppos
					, cc22 ? pepShortName[(int)cc22] : "" 
					) ;
			  else
			    vtxtPrintf (txt, "\nAA_substitution  %s%ddelins%s%s"
					, cc11 ? pepShortName[(int)cc11] : "?"
					, ppos
					, cc21 ? pepShortName[(int)cc21] : "?"
					, cc22 ? pepShortName[(int)cc22] : "?" 
					) ;
			}	
		      else if (isFrameshift % 3) 
			{	 
			  if (cc11 == cc21)
			    {
			      vtxtPrintf (txt, "Frameshift  %s%d%sfs"
					  , cc12 ? pepShortName[(int)cc12] : ""
					  , ppos
					  , cc22 ? pepShortName[(int)cc22] : ""
					  ) ;
			      stop2-- ;
			    }
			  else
			    vtxtPrintf (txt, "Frameshift  %s%d%sfs"
					, cc11 ? pepShortName[(int)cc11] : "?"
					, ppos
					, cc21 ? pepShortName[(int)cc21]: ""
					) ;
			  if (stop2 >= 0)
			    vtxtPrintf (txt, "Ter%d", stop2) ;
			  else
			    vtxtPrintf (txt, "Ter") ;
			}
		      vtxtPrint (txt, "\n") ;
		    }

		  if (!  ac_table_printable (mrnas, iMrnaWithCDS, 1, " ") )
		    messcrash ("bizare") ;

		  /* set non-common letters of the protein in lower case */
		  /* set all etters to upper */
		  cp1 = tbuf1 ; cp2 = tbuf2 ;
		  for (cp1 = tbuf1 ; *cp1 ; cp1++) *cp1 = ace_upper (*cp1) ;
		  for (cp2 = tbuf2 ; *cp2 ; cp2++) *cp2 = ace_upper (*cp2) ;
		  /* keep common fist letters upper */
		  cp1 = tbuf1 ; cp2 = tbuf2 ;
		  while (*cp1 && *cp1 == *cp2) { cp1++ ; cp2++ ; } 
		   /* lower the rest */
		  if (*cp1) { while (*cp1) { *cp1 = ace_lower (*cp1) ; cp1++ ; }} 
		  if (*cp2) { while (*cp2) { *cp2 = ace_lower (*cp2) ; cp2++ ; }}
		  /* re-upper common terminal letters */
		  while (--cp1 >= tbuf1 && --cp2 >= tbuf2 && *cp1 == *cp2) 
		    { *cp1 = ace_upper (*cp1) ;  *cp2 = ace_upper (*cp2) ; }

		  if (0)
		    {
		      vtxtPrintf(txt, "Reference_sequence  %s %s %s %s %d\n"
				 , buf1, tbuf1
				 , ac_table_printable (mrnas, iMrnaWithCDS, 1, " ") 
				 , ac_table_printable (mrnas, iMrnaWithCDS, 4, " ") 
				 , 1 + (pos - x1) / 3
				 ) ;
		      vtxtPrintf(txt, "Observed__sequence   %s %s%s %s %s %d\n"
				 , buf2, tbuf2
				 , 0 && ((999999 + isFrameshift) % 3) ? "fs" : "" /* we no longer add fs, since we now put the peptide in lower case */
				 , ac_table_printable (mrnas, iMrnaWithCDS, 1, " ") 
				 , ac_table_printable (mrnas, iMrnaWithCDS, 4, " ") 
				 , 1 + (pos - x1) / 3
				 ) ;
		    }
		}
	    }
	}
      
      /* edit this variant */
      vtxtPrintf (txt, "\n") ;
      if (ao)
	aceOutf (ao, "%s\n", vtxtPtr (txt)) ;
      else
	{
	  ac_parse (db, vtxtPtr (txt), &errors, 0, h1) ; 
	  if (errors && *errors)
	    messerror (errors) ;
	}
      vtxtClear (txt) ;
      if (1)
	{
	  vtxtClear (gTxt) ;
	  vtxtClear (cTxt) ;
	  vtxtClear (pTxt) ;
	  vtxtClear (geneboxes) ;
	  vtxtClear (avgeneboxes) ;
	  vtxtClear (gsnippet) ;
	  vtxtClear (snippet) ;
	  vtxtClear (pSnippet) ;
	  vtxtClear (pType) ;
	  
	  /*
	    const char *ccp ;
	    ccp = ac_tag_printable (variant, "Typ", 0) ;
	    if (ccp && strlen (ccp) < 64)
	    {
	    int gMap = 0, gPos = 0 ;
	    char gType[64], cType[64] ; 
	    gType[0] = 0 ;
	    cType[0] = 0 ; 
	    
	    strcpy (gType, ccp) ;
	      snpPrettyNames (snp, variant, gTxt, cTxt, pTxt, &gMap, &gPos, gType, cType, pType, location, geneboxes, avgeneboxes, gsnippet, snippet, pSnippet) ;
	    }
	  */
	}
    }
  if (tsnp != tsnp0)
    invokeDebugger () ;

  tsnpPotential_splice_disruption (tsnp, ao) ;

  ac_free (h1) ;
  ac_free (h) ;

  return nn ;
} /* tsnpCodingModif */

/*************************************************************************************/
/*************************************************************************************/
#define SNPTYPEMAX 1024
static void tsnpTypeInit (TSNP *tsnp)
{
  char *snpTypes = strnew ("Genomic,Exonic,Protein_changing,A>G,T>C,G>A,C>T,A>T,T>A,G>C,C>G,A>C,T>G,G>T,C>A,InsA,InsT,InsG,InsC,DupA,DupT,DupG,DupC,DelA,DelT,DelG,DelC,DimA,DimT,DimG,DimC,InsAA,InsTT,InsGG,InsCC,InsAG,InsCT,InsAC,InsGT,InsTG,InsCA,InsTC,InsGA,InsAT,InsTA,InsGC,InsCG,DupAA,DupTT,DupGG,DupCC,DupAG,DupCT,DupAC,DupGT,DupTG,DupCA,DupTC,DupGA,DupAT,DupTA,DupGC,DupCG,DelAA,DelTT,DelGG,DelCC,DelAG,DelCT,DelAC,DelGT,DelTG,DelCA,DelTC,DelGA,DelAT,DelTA,DelGC,DelCG,DimAA,DimTT,DimGG,DimCC,DimAG,DimCT,DimAC,DimGT,DimTG,DimCA,DimTC,DimGA,DimAT,DimTA,DimGC,DimCG,InsAAA,InsTTT,InsGGG,InsCCC,InsAAT,InsATT,InsAAG,InsCTT,InsAAC,InsGTT,InsTTA,InsTAA,InsTTG,InsCAA,InsTTC,InsGAA,InsGGA,InsTCC,InsGGT,InsACC,InsGGC,InsGCC,InsCCA,InsTGG,InsCCT,InsAGG,InsCCG,InsCGG,InsATA,InsTAT,InsATG,InsCAT,InsATC,InsGAT,InsAGA,InsTCT,InsAGT,InsACT,InsAGC,InsGCT,InsACA,InsTGT,InsACG,InsCGT,InsTAG,InsCTA,InsTAC,InsGTA,InsTGA,InsTCA,InsTGC,InsGCA,InsTCG,InsCGA,InsGAG,InsCTC,InsGAC,InsGTC,InsGTG,InsCAC,InsGCG,InsCGC,InsCAG,InsCTG,DupAAA,DupTTT,DupGGG,DupCCC,DupAAT,DupATT,DupAAG,DupCTT,DupAAC,DupGTT,DupTTA,DupTAA,DupTTG,DupCAA,DupTTC,DupGAA,DupGGA,DupTCC,DupGGT,DupACC,DupGGC,DupGCC,DupCCA,DupTGG,DupCCT,DupAGG,DupCCG,DupCGG,DupATA,DupTAT,DupATG,DupCAT,DupATC,DupGAT,DupAGA,DupTCT,DupAGT,DupACT,DupAGC,DupGCT,DupACA,DupTGT,DupACG,DupCGT,DupTAG,DupCTA,DupTAC,DupGTA,DupTGA,DupTCA,DupTGC,DupGCA,DupTCG,DupCGA,DupGAG,DupCTC,DupGAC,DupGTC,DupGTG,DupCAC,DupGCG,DupCGC,DupCAG,DupCTG,DelAAA,DelTTT,DelGGG,DelCCC,DelAAT,DelATT,DelAAG,DelCTT,DelAAC,DelGTT,DelTTA,DelTAA,DelTTG,DelCAA,DelTTC,DelGAA,DelGGA,DelTCC,DelGGT,DelACC,DelGGC,DelGCC,DelCCA,DelTGG,DelCCT,DelAGG,DelCCG,DelCGG,DelATA,DelTAT,DelATG,DelCAT,DelATC,DelGAT,DelAGA,DelTCT,DelAGT,DelACT,DelAGC,DelGCT,DelACA,DelTGT,DelACG,DelCGT,DelTAG,DelCTA,DelTAC,DelGTA,DelTGA,DelTCA,DelTGC,DelGCA,DelTCG,DelCGA,DelGAG,DelCTC,DelGAC,DelGTC,DelGTG,DelCAC,DelGCG,DelCGC,DelCAG,DelCTG,DimAAA,DimTTT,DimGGG,DimCCC,DimAAT,DimATT,DimAAG,DimCTT,DimAAC,DimGTT,DimTTA,DimTAA,DimTTG,DimCAA,DimTTC,DimGAA,DimGGA,DimTCC,DimGGT,DimACC,DimGGC,DimGCC,DimCCA,DimTGG,DimCCT,DimAGG,DimCCG,DimCGG,DimATA,DimTAT,DimATG,DimCAT,DimATC,DimGAT,DimAGA,DimTCT,DimAGT,DimACT,DimAGC,DimGCT,DimACA,DimTGT,DimACG,DimCGT,DimTAG,DimCTA,DimTAC,DimGTA,DimTGA,DimTCA,DimTGC,DimGCA,DimTCG,DimCGA,DimGAG,DimCTC,DimGAC,DimGTC,DimGTG,DimCAC,DimGCG,DimCGC,DimCAG,DimCTG", 0) ;  

  if (! tsnp->snpTypeDict)
    {
      char *cp, *cq ;
      DICT *dict ;
      
      dict = tsnp->snpTypeDict = dictHandleCreate (SNPTYPEMAX, tsnp->h) ;

      cp = snpTypes ;
      while (cp)
	{
	  cq = strchr (cp, ',') ;
	  if (cq) *cq = 0 ;
	  dictAdd (dict, cp, 0) ;
	  cp = cq ? cq + 1 : 0 ;
	}
    }
  ac_free (snpTypes) ;

  return ;
} /* tsnpTypesInit */

/*************************************************************************************/

typedef struct snpProfileStruct {
  BOOL monomodal ;
  int ref, low, mid, high, pure, any ;
  int gRef, gLow, gMid, gHigh, gPure, gAny ;
  int xRef, xLow, xMid, xHigh, xPure, xAny ;
  int pcRef, pcLow, pcMid, pcHigh, pcPure, pcAny ;
  int typeN[SNPTYPEMAX], typeD[SNPTYPEMAX] ;
} SP ;

/*************************************************************************************/

static DICT *tsnpGetRuns (TSNP *tsnp)
{
  DICT *runDict = tsnp->runDict ; 

  if (! runDict)
    {
      const char *errors = 0 ;
      AC_HANDLE h = ac_new_handle () ;
      char *qq = hprintf (h, "Find project %s ; >run ; CLASS runs", tsnp->project) ;
      AC_TABLE tbl = ac_bql_table (tsnp->db, qq, 0, 0, &errors, h) ;
      int ir ;

      runDict = tsnp->runDict = dictHandleCreate (156, tsnp->h) ;
      if (tbl)
	for (ir = 0 ; ir < tbl->rows ; ir++)
	  dictAdd (runDict, ac_table_printable (tbl, ir, 0, "toto"), 0) ;
      ac_free (h)  ;
    }

  return runDict ;
} /* tsnpGetRuns */

/*************************************************************************************/

static BOOL tsnpGetMonomodal (TSNP *tsnp, AC_OBJ Snp)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ;

  int run ;
  const char *typ ;
  int t ;
  SP *sp2, sp0 ;
  DICT *runDict = tsnpGetRuns (tsnp) ;
  Array aa = tsnp->runProfiles ;

  if (! aa)
    aa = tsnp->runProfiles = arrayHandleCreate (128, SP, tsnp->h) ; 
  tsnpTypeInit (tsnp) ;

  sp2 = &sp0 ;
  memset (sp2, 0, sizeof (SP)) ;
  h1 = ac_new_handle () ;
  typ = ac_tag_printable (Snp, "Typ", 0) ;
  if (typ && dictFind (tsnp->snpTypeDict, typ, &t))
    {
      int ir ;
      AC_TABLE tbl = ac_tag_table (Snp, "BRS_Counts", h1) ;
      
      for (ir = 0 ; ir < tbl->rows ; ir++)
	{
	  int c = ac_table_int (tbl, ir, 1, 0) ;
	  int m = ac_table_int (tbl, ir, 2, 0) ;
	  
	  if (! dictFind (runDict, ac_table_printable (tbl, ir, 0, "toto"), &run))
	    continue ;

	  SP *sp = arrayp (aa, run, SP) ;

	  if (c >= tsnp->minSnpCover)
	    {
	      float f = 100.0 * m / c ;
	      

	      sp->typeN[t]++ ;
	      if (1)
		{
		  sp->any++ ;   /* this run any Snp */
		  if (f <= 5)
		    sp->ref++ ;
		  else if (f <= 20)
		    sp->low++ ;
		  else if (f <= 80)
		    sp->mid++ ;
		  else if (f <= 95)
		    sp->high++ ;
		  else 
		    sp->pure++ ;

		  sp2->any++ ;  /* this Snp in any run */
		  if (f <= 5)
		    sp2->ref++ ;
		  else if (f <= 20)
		    sp2->low++ ;
		  else if (f <= 80)
		    sp2->mid++ ;
		  else if (f <= 95)
		    sp2->high++ ;
		  else 
		    sp2->pure++ ;
		}
	      if (ac_has_tag (Snp, "Coding"))
		{
		  sp->pcAny++ ;
		  if (f <= 5)
		    sp->pcRef++ ;
		  else if (f <= 20)
		    sp->pcLow++ ;
		  else if (f <= 80)
		    sp->pcMid++ ;
		  else if (f <= 95)
		    sp->pcHigh++ ;
		  else 
		    sp->pcPure++ ;
		}
	      else if (ac_has_tag (Snp, "Exonic"))
		{
		  sp->xAny++ ;
		  if (f <= 5)
		    sp->xRef++ ;
		  else if (f <= 20)
		    sp->xLow++ ;
		  else if (f <= 80)
		    sp->xMid++ ;
		  else if (f <= 95)
		    sp->xHigh++ ;
		  else 
		    sp->xPure++ ;
		}
	      else
		{
		  sp->gAny++ ;
		  if (f <= 5)
		    sp->gRef++ ;
		  else if (f <= 20)
		    sp->gLow++ ;
		  else if (f <= 80)
		    sp->gMid++ ;
		  else if (f <= 95)
		    sp->gHigh++ ;
		  else 
		    sp->gPure++ ;
		}
	    }
	}

      
      if (10 * sp2->low > sp2->any && sp2->low > 5 &&
	  ! (sp2->ref > sp2->low + 10 && sp2->mid + sp2->pure > sp2->low + 10) &&
	  ! (sp2->pure > sp2->high)
	  )
	sp2->monomodal =  TRUE ;
      else
	sp2->monomodal = FALSE ;
      
      if (sp2->monomodal) 
	{  /* count them in the run */
	  for (ir = 0 ; ir < tbl->rows ; ir++)
	    {
	      int c = ac_table_int (tbl, ir, 1, 0) ;
	      int m = ac_table_int (tbl, ir, 2, 0) ;
	      
	      if (! dictFind (runDict, ac_table_printable (tbl, ir, 0, "toto"), &run))
		continue ;
	      
	      SP *sp = arrayp (aa, run, SP) ;
	      
	      if (c >= tsnp->minSnpCover && m >= tsnp->minSnpCount)
		sp->monomodal++ ;
	    }
	}	
    }

  ac_free (h) ;
  return sp2->monomodal ;
} /* tsnpGetMonomodal */

/*************************************************************************************/
/* export tsf file for this section */
static void tsnpExportProfile (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (tsnp->outFileName, ".snp_profile.tsf", tsnp->gzo, h) ;
  int run ;
  int typeMax = 0 ;
  int runMax ;
  DICT *runDict = tsnp->runDict ; 
  Array aa = tsnp->runProfiles ;

  tsnpTypeInit (tsnp) ;
  if (! runDict)
    runDict = tsnp->runDict = dictHandleCreate (156, tsnp->h) ;

  typeMax = dictMax (tsnp->snpTypeDict) ;
  runMax = aa ? arrayMax (aa) : 0 ;
  for (run = 1 ; run <= runMax ; run++)
    {
      SP *sp = arrayp (aa, run, SP) ;
      if (sp->any)
	{
	  int t ;
	  const char* runName = dictName (tsnp->runDict, run) ;
	  /* export the categories */
	  if (sp->any)
	    aceOutf (ao, "%s\tAny\tititititititt\t%d\tAny\t%d\treference\t%d\tlow\t%d\tmid\t%d\thigh\t%d\tpure\tmonomodal\n"
		     , runName, sp->gAny, sp->gRef,sp->gLow,sp->gMid,sp->gHigh,sp->gPure
		     , sp->monomodal ? "Monomodal" : ""
		     ) ; 
	  if (sp->gAny)
	    aceOutf (ao, "%s\tGenomic\titititititit\t%d\tAny\t%d\treference\t%d\tlow\t%d\tmid\t%d\thigh\t%d\tpure\n"
		     , runName, sp->gAny, sp->gRef,sp->gLow,sp->gMid,sp->gHigh,sp->gPure) ; 
	  if (sp->xAny)
	    aceOutf (ao, "%s\tExonic\titititititit\t\t%d\tAny\t%d\treference\t%d\tlow\t%d\tmid\t%d\thigh\t%d\tpure\n"
		     , runName, sp->xAny, sp->xRef,sp->xLow,sp->xMid,sp->xHigh,sp->xPure) ; 
	  if (sp->pcAny)
	    aceOutf (ao, "%s\tProtein_changing\titititititit\t\t%d\tAny\t%d\treference\t%d\tlow\t%d\tmid\t%d\thigh\t%d\tpure\n"
		     , runName, sp->pcAny, sp->pcRef,sp->pcLow,sp->pcMid,sp->pcHigh,sp->pcPure) ; 
	  
	  /* export the profile */
	  for (t = 1 ; t < typeMax ; t++)
	    {
	      if (sp->typeN[t])
		aceOutf (ao, "%s\t%d__%s\tii\t%d\t%d\n"
			 , runName, t, dictName(tsnp->snpTypeDict,t), sp->typeN[t], sp->typeD[t]) ;
	    }
	}
    }
  ac_free (h) ;
  return;

} /* tsnpExportProfile */

/*************************************************************************************/
/* analyse the rejected and non rejected snippets */
typedef struct gggStruct { int snp ; BOOL mono ; int run, ncounts, cover, mp, mm, wp, wm ; } GGG ;

static void tsnpGGG (TSNP *tsnp)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  ACEOUT ao = 0 ; 
  AC_ITER iter = 0 ;
  AC_OBJ variant = 0 ;
  int run ;
  DICT *runDict = tsnp->runDict ; 
  Array aa = arrayHandleCreate (64, SP, h) ; 
  DICT *gggDict1 = dictHandleCreate (1024, h) ;
  Array ggg1 = arrayHandleCreate (100000, GGG, h) ;  
  DICT *gggDict2 = dictHandleCreate (1024, h) ;
  Array ggg2 = arrayHandleCreate (100000, GGG, h) ;  
  DICT *gggDict3 = dictHandleCreate (1024, h) ;
  Array ggg3 = arrayHandleCreate (100000, GGG, h) ;  
  KEYSET ks1 = keySetHandleCreate (h) ;
  KEYSET ks2 = keySetHandleCreate (h) ;
  KEYSET ks3 = keySetHandleCreate (h) ;

  tsnpTypeInit (tsnp) ;
  if (! runDict)
    runDict = tsnp->runDict = dictHandleCreate (156, tsnp->h) ;


  if (tsnp->db)
    {
      iter = ac_query_iter (tsnp->db, TRUE, "find variant Typ && BRS_counts", 0, h) ;
      while (ac_free (variant), ac_free (h1), variant = ac_iter_obj (iter))
	{
	  const char *typ ;
	  int t ;
	  const char *dna = 0 ;
	  SP *sp = 0 ;
	  BOOL mono = ac_has_tag (variant, "Monomodal") ; 

	  h1 = ac_new_handle () ;
	  typ = ac_tag_printable (variant, "Typ", 0) ;
	  if (typ && dictFind (tsnp->snpTypeDict, typ, &t))
	    {
	      int ir ;
	      AC_TABLE tbl = ac_tag_table (variant, "BRS_Counts", h1) ;

	      for (ir = 0 ; ir < tbl->rows ; ir++)
		{
		  int c = ac_table_int (tbl, ir, 1, 0) ;
		  int m = ac_table_int (tbl, ir, 2, 0) ;
		  if (c > tsnp->minSnpCover)
		    {
		      float f = 100.0 * m / c ;

		      dictAdd (runDict, ac_table_printable (tbl, ir, 0, "toto"), &run) ;
		      sp = arrayp (aa, run, SP) ;
		      sp->typeN[t]++ ;
		      if (1)
			{
			  sp->any++ ;
			  if (f <= 5)
			    sp->ref++ ;
			  else if (f <= 20)
			    sp->low++ ;
			  else if (f <= 80)
			    sp->mid++ ;
			  else if (f <= 95)
			    sp->high++ ;
			  else 
			    sp->pure++ ;
			}
		      if (ac_has_tag (variant, "Coding"))
			{
			  sp->pcAny++ ;
			  if (f <= 5)
			    sp->pcRef++ ;
			  else if (f <= 20)
			    sp->pcLow++ ;
			  else if (f <= 80)
			    sp->pcMid++ ;
			  else if (f <= 95)
			    sp->pcHigh++ ;
			  else 
			    sp->pcPure++ ;
			}
		      else if (ac_has_tag (variant, "Exonic"))
			{
			  sp->xAny++ ;
			  if (f <= 5)
			    sp->xRef++ ;
			  else if (f <= 20)
			    sp->xLow++ ;
			  else if (f <= 80)
			    sp->xMid++ ;
			  else if (f <= 95)
			    sp->xHigh++ ;
			  else 
			    sp->xPure++ ;
			}
		      else
			{
			  sp->gAny++ ;
			  if (f <= 5)
			    sp->gRef++ ;
			  else if (f <= 20)
			    sp->gLow++ ;
			  else if (f <= 80)
			    sp->gMid++ ;
			  else if (f <= 95)
			    sp->gHigh++ ;
			  else 
			    sp->gPure++ ;
			}
		    }
		}
	    
	      /* knowing if monomodal is TRUE or FALSE, analyze the snippets */
	      /* get the dna */
	      if (ac_has_tag (variant, "Found_in_mRNA"))
		dna = ac_tag_printable (variant, "Reference_RNAexon_sequence", 0) ;
	      else if (ac_has_tag (variant, "Found_in_genone"))
		dna = ac_tag_printable (variant, "Reference_genome_sequence", 0) ;
	      if (dna && typ && strlen(typ) < 200)
		{
		  /* locate the Upper case letter */
		  const char *ccp = dna ;
		  int pos = 0 ;
		  char buf[256] ;
		  while (*ccp)
		    {
		      if (*ccp == ace_upper(*ccp))
			break ;
		      pos++ ; ccp++ ;
		    }
		  if (*ccp) /* we located the SNP */
		    {
		      /* grab the 3 letters in front */
		      if (pos > 3)
			{
			  int i, ir, snp ;
			  for (i = 0 ; i < 3 ; i++)
			    {
			      buf[i] = ace_upper(dna[pos-3+i]) ;
			      if (buf[i] == 'U') buf[i] = 'T' ;
			    }
			  buf[i] = ace_lower(dna[pos]) ;
			  if (buf[i] == 'u') buf[i] = 't' ;
			  buf[4] = ':' ;
			  memcpy (buf+5, typ, strlen(typ)+1) ;
			  dictAdd (gggDict1, buf, &snp) ; 
			  keySet (ks1, snp)++ ;
			  for (ir = 0 ; ir < tbl->rows ; ir++)
			    {
			      GGG *xp = arrayp (ggg1, arrayMax (ggg1), GGG) ;
			      xp->snp = snp ;
			      xp->ncounts++ ;
			      xp->mono |= mono ;
			      xp->run = ac_table_key (tbl, ir, 0, 0) ;
			      xp->cover = ac_table_int (tbl, ir, 1, 0) ;
			      xp->mp = ac_table_int (tbl, ir, 4, 0) ;
			      xp->wp = ac_table_int (tbl, ir, 5, 0) ;
			      xp->mm = ac_table_int (tbl, ir, 6, 0) ;
			      xp->wm = ac_table_int (tbl, ir, 7, 0) ;
			    }
			}
		      /* grab the 3 letters behind */
		      if (pos < strlen(dna) - 3)
			{
			  int i = 0, ir, snp ;
			  buf[0] = ace_lower(dna[pos]) ;
			  if (buf[0] == 'u') buf[0] = 't' ;
			  for (i = 1 ; i < 4 ; i++)
			    {
			      buf[i] = ace_upper(dna[pos+i]) ;
			      if (buf[i] == 'U') buf[i] = 'T' ;
			    }
			  buf[4] = ':' ;
			  memcpy (buf+4, typ, strlen(typ)+1) ;
			  dictAdd (gggDict2, buf, &snp) ; 
			  keySet (ks2, snp)++ ;
			  for (ir = 0 ; ir < tbl->rows ; ir++)
			    {
			      GGG *xp = arrayp (ggg2, arrayMax (ggg2), GGG) ;
			      xp->snp = snp ;
			      xp->mono |= mono ;
			      xp->ncounts++ ;
			      xp->run = ac_table_key (tbl, ir, 0, 0) ;
			      xp->cover = ac_table_int (tbl, ir, 1, 0) ;
			      xp->mp = ac_table_int (tbl, ir, 4, 0) ;
			      xp->wp = ac_table_int (tbl, ir, 5, 0) ;
			      xp->mm = ac_table_int (tbl, ir, 6, 0) ;
			      xp->wm = ac_table_int (tbl, ir, 7, 0) ;
			    }
			}
		      /* grab the 5 letters centered */
		      if (pos > 2 && pos < strlen(dna) - 3)
			{
			  int i, ir, snp ;
			  for (i = 0 ; i < 5 ; i++)
			    {
			      buf[i] = ace_upper(dna[pos+i-2]) ;
			      if (buf[i] == 'U') buf[i] = 'T' ;
			    }
			  buf[2] = ace_lower(dna[pos]) ;
			  if (buf[2] == 'u') buf[2] = 't' ;
			  buf[5] = ':' ;
			  memcpy (buf+6, typ, strlen(typ)+1) ;
			  dictAdd (gggDict3, buf, &snp) ; 
			  keySet (ks3, snp)++ ;
			  for (ir = 0 ; ir < tbl->rows ; ir++)
			    {
			      GGG *xp = arrayp (ggg3, arrayMax (ggg3), GGG) ;
			      xp->snp = snp ;
			      xp->ncounts++ ;
			      xp->mono |= mono ;
			      xp->run = ac_table_key (tbl, ir, 0, 0) ;
			      xp->cover = ac_table_int (tbl, ir, 1, 0) ;
			      xp->mp = ac_table_int (tbl, ir, 4, 0) ;
			      xp->wp = ac_table_int (tbl, ir, 5, 0) ;
			      xp->mm = ac_table_int (tbl, ir, 6, 0) ;
			      xp->wm = ac_table_int (tbl, ir, 7, 0) ;
			    }
			}
		    }
		}
	    }
	}
    }

  if (1)
    {
      ao = aceOutCreate (tsnp->outFileName, ".GGG.before.tsf", tsnp->gzo, h) ;
      aceOutDate (ao, "###", "Triplet before") ;
      aceOutf (ao, "# Type\tRun\tciiiit\tNsnps\tNcounts\tCover\tmp\twp\tmm\twm\tmonomodal\n") ;
      for (int i = 0 ; i < arrayMax (ggg1) ; i++)
	{
	  GGG *xp = arrayp (ggg1, i, GGG) ;
	  if (xp->snp)
	    aceOutf (ao, "%s\t%s\tciiiiiit\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"
		     , dictName (gggDict1, xp->snp)
		     , name (xp->run)
		     , keySet (ks1, xp->snp)
		     , xp->ncounts
		     , xp->cover, xp->mp, xp->wp, xp->mm, xp->wm
		     , xp->mono ? "monomodal" : "-"
		     ) ;
	}
    }
  if (1)
    {
      ao = aceOutCreate (tsnp->outFileName, ".GGG.after.tsf", tsnp->gzo, h) ;
      aceOutDate (ao, "###", "Triplet after") ;
      aceOutf (ao, "# Type\tRun\tciiiit\tNsnps\tNcounts\tCover\tmp\twp\tmm\twm\tmonomodal\n") ;
      for (int i = 0 ; i < arrayMax (ggg2) ; i++)
	{
	  GGG *xp = arrayp (ggg2, i, GGG) ;
	  if (xp->snp)
	    aceOutf (ao, "%s\t%s\tciiiiiit\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"
		     , dictName (gggDict2, xp->snp)
		     , name (xp->run)
		     , keySet (ks2, xp->snp)
		     , xp->ncounts
		     , xp->cover, xp->mp, xp->wp, xp->mm, xp->wm
		     , xp->mono ? "monomodal" : "-"
		     ) ;
	}
    }
  if (1)
    {
      ao = aceOutCreate (tsnp->outFileName, ".GGG.centered.tsf", tsnp->gzo, h) ;
      aceOutDate (ao, "###", "% letters centered") ;
      aceOutf (ao, "# Type\tRun\tciiiit\tNsnps\tNcounts\tCover\tmp\twp\tmm\twm\tmonomodal\n") ;
      for (int i = 0 ; i < arrayMax (ggg3) ; i++)
	{
	  GGG *xp = arrayp (ggg3, i, GGG) ;
	  if (xp->snp)
	    aceOutf (ao, "%s\t%s\tciiiiiit\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"
		     , dictName (gggDict3, xp->snp)
		     , name (xp->run)
		     , keySet (ks3, xp->snp)
		     , xp->ncounts
		     , xp->cover, xp->mp, xp->wp, xp->mm, xp->wm
		     , xp->mono ? "monomodal" : "-"
		     ) ;
	}
    }
			       
  ac_free (h) ;
  return;

} /* tsnpGGG */

/*************************************************************************************/
/*************************************************************************************/

typedef struct subsStruct { KEY snp, target ; int a1, a2 ; } SUBS ;

/*************************************************************************************/
 
static int tsnpSubsOrder  (const void *va, const void *vb)
{
  const SUBS *up = (const SUBS *)va, *vp = (const SUBS *)vb ;
  int n ;
  
  n = up->target - vp->target ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;

  return 0 ;
} /* tsnpSubsOrder */

/*************************************************************************************/

static int tsnpMakeWords_Any (TSNP *tsnp, ACEOUT ao, BOOL isMrna)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl = 0 ;
  char *qq ;
  KEYSET ksIn = 0 ;
  const char *zone = tsnp->zone ? tsnp->zone : "" ;
  const char *errors = 0 ;
  char *bqlQuery = 0 ;
  BQL_ITER *bqlIter = 0 ;
  int nn = 0, nn1 = 0, nSub = 0 ;
  Array subs = arrayHandleCreate (500, SUBS, h) ;
  SUBS *up, *vp ;
  if (isMrna)
    {
      qq = hprintf (h, "%s %s %s %s " 
		    , "Find variant ! Found_in_genome && mRNA" 
		    , tsnp->filter ? " && (" : ""
		    , tsnp->filter ? tsnp->filter : ""
		    , tsnp->filter ? ") " : ""
		    ) ;
      bqlQuery = hprintf (h, "%s %s"
			  , " select v,ref,var,del, mdel, ins, mins   "
			  , "from v in @, ref in v->Reference_RNAexon_sequence , var in v->Observed__RNAexon_sequence where ref and var, del in v->deletion, mdel in v->multi_deletion, ins in v->insertion, mins in v->multi_insertion"
			  ) ;
    }
  else
    {
      qq = hprintf (h, "%s %s %s %s"
		    , "Find variant Found_in_genome" 
		    , tsnp->filter ? " && (" : ""
		    , tsnp->filter ? tsnp->filter : ""
		    , tsnp->filter ? ") " : ""
		    ) ;

      bqlQuery = hprintf (h, "%s %s "
			  , " select v,ref,var,del, mdel, ins, mins, sub, target, a1,a2   "
			  , "from v in @, ref in v->Reference_genomic_sequence , var in v->Observed__genomic_sequence where ref and var, del in v->deletion, mdel in v->multi_deletion, ins in v->insertion, mins in v->multi_insertion, sub in v#substitution, target in v->IntMap, a1 in target[1], a2 in target[2]"
			  ) ;
    }
  ksIn = query (0, qq) ;
  if (keySetMax (ksIn))   bqlIter = bqlIterCreate (bqlQuery, ksIn, &errors, h) ;
  if (bqlIter) while ((tbl = bqlIterTable (bqlIter)))
    {
      KEY v = ac_table_key (tbl, 0, 0, 0) ;
      const char *ref = ac_table_printable (tbl, 0, 1, "xxx") ;
      const char *var = ac_table_printable (tbl, 0, 2, "xxx") ;
      const char *del = ac_table_printable (tbl, 0, 3, 0) ;
      const char *ins = ac_table_printable (tbl, 0, 5, 0) ;
      int mdel = ac_table_int (tbl, 0, 4, 0) ;
      int mins = ac_table_int (tbl, 0, 6, 0) ;
      KEY sub = ac_table_key (tbl, 0, 7, 0) ;
      KEY target = ac_table_key (tbl, 0, 8, 0) ;
      int a1 = ac_table_int (tbl, 0, 9, 0) ;
      int a2 = ac_table_int (tbl, 0, 10, 0) ;
      char buf[128] ;

      nn1++ ;
      if (ref && strlen (ref) > 30 && var && strlen (var) > 30)
	{
	  if (ins && mins == 0)
	    { 
	      memcpy (buf, ref + 6, 15) ;
	      buf[15]='@' ;
	      memcpy (buf+15, ref + 22, 16) ;
	    }
	  else if (ins && mins)
	    {
	      const char *cr = ref + strlen(ref) - 1 ;
	      memcpy (buf, ref + 6, 15) ;
	      while (cr > ref + 20 && *cr != '^') 
		cr-- ;
	      memcpy (buf+15, cr+1, 16) ;
	    }
	  else
	    memcpy (buf, ref + 5 + (del ? 1 : 0), 31) ;
	  buf[32] = 0 ;
	  aceOutf (ao, "%s\tRef__%s\t%s\n", buf, name(v), zone) ;
	  if (sub)
	    {
	      up = arrayp (subs, nSub++, SUBS) ;
	      up->snp = v ;
	      up->a1 = a1 ;
	      up->a2 = a2 ;
	      up->target = target ;
	    }
	  if (del && mdel == 0)
	    { 
	      memcpy (buf, var + 6, 15) ;
	      memcpy (buf+15, var + 22, 16) ;
	    }
	  else if (del && mdel)
	    {
	      const char *cr = var + strlen(var) - 1 ;
	      memcpy (buf, var + 6, 15) ;
	      while (cr > var + 20 && *cr != '_') 
		cr-- ;
	      memcpy (buf+15, cr+1, 16) ;
	    }
	  else
	    memcpy (buf, var + (ins ? 6 : 5), 31) ;
	  buf[32] = 0 ;
	  aceOutf (ao, "%s\tVar__%s\t%s\n", buf, name(v), zone) ;
	  nn++ ;
	}
      if (nn < 0) break ;
    }

  if (nSub > 1)
    {
      int ii, jj ;
      arraySort (subs, tsnpSubsOrder) ;
      
      for (ii = 0, up = arrp (subs, 0, SUBS) ; ii < nSub ; up++, ii++)
	{
	  KEY v = up->snp ;
	  AC_OBJ  Snp1 = ac_get_obj (tsnp->db, "Variant", name (v), 0) ; 
	  const char *obs1 = Snp1 ? ac_tag_text (Snp1, "Observed__genomic_sequence", 0) : 0 ;
	  const char *ref1 = Snp1 ? ac_tag_text (Snp1, "Observed__genomic_sequence", 0) : 0 ;
	  char buf[32] ;

	  if (! obs1 || !ref1)
	    continue ;
	  memset (buf, 0, sizeof(buf)) ;
	  for (vp = up + 1, jj = ii + 1 ; jj < nSub && vp->target == up->target && vp->a1 < up->a1 + 15 ; jj++, vp++)
	    {
	      AC_OBJ  Snp2 = ac_get_obj (tsnp->db, "Variant", name (vp->snp), 0) ; 
	      const char *obs2 = Snp2 ? ac_tag_text (Snp2, "Observed__genomic_sequence", 0) : 0 ;

	      if (0 && obs2)
		{
		  aceOutf (ao, "%s\tVar_%d__%s\t%s:%d\n", buf, vp->a1, name(v), zone) ;
		  aceOutf (ao, "%s\tRef_%d__%s\t%s:%d\n", buf, vp->a1, name(v), zone) ;
		}
	      ac_free (Snp2) ;
	    }
	  for (vp = up - 1, jj = ii - 1 ; jj >= 0 && vp->target == up->target && vp->a1 > up->a1 - 15 ; jj--, vp--)
	    {
	      AC_OBJ  Snp2 = ac_get_obj (tsnp->db, "Variant", name (vp->snp), 0) ; 
	      const char *obs2 = Snp2 ? ac_tag_text (Snp2, "Observed__genomic_sequence", 0) : 0 ;

	      if (0 && obs2)
		{
		  aceOutf (ao, "%s\tVar_%d__%s\t%s:%d\n", buf, vp->a1, name(v), zone) ;
		  aceOutf (ao, "%s\tRef_%d__%s\t%s:%d\n", buf, vp->a1, name(v), zone) ;
		}
	      ac_free (Snp2) ;
	    }
	  ac_free (Snp1) ;
	}
    }

  ac_free (ksIn) ;
  ac_free (h) ;
  fprintf (stderr, "# tsnpMakeWords_Sub isMrna %s exported %d/%d snps\n", isMrna ? "TRUE" : "FALSE", nn, nn1) ;
  return nn ;
} /* tsnpMakeWords_Any */

/*************************************************************************************/

static int tsnpMakeWords_Sub (TSNP *tsnp, ACEOUT ao, BOOL isMrna)
{
  AC_HANDLE h1 =0, h = ac_new_handle () ;
  AC_TABLE tbl = 0 ;
  char *qq ;
  KEYSET ksIn = 0 ;
  const char *zone = tsnp->zone ? tsnp->zone : "" ;
  const char *errors = 0 ;
  char *bqlQuery = 0 ;
  BQL_ITER *bqlIter = 0 ;
  int nn = 0, nn1 = 0 ;

  if (isMrna)
    {
      qq = hprintf (h, "%s %s %s %s " 
		    , "Find variant substitution && ! multi_substitution && ! Found_in_genome && mRNA" 
		    , tsnp->filter ? " && (" : ""
		    , tsnp->filter ? tsnp->filter : ""
		    , tsnp->filter ? ") " : ""
		    ) ;
      bqlQuery = hprintf (h, "%s %s"
			  , " select v,typ,dna,a1    "
			  , "from v in @, typ in v->typ, m in v->mRNA, a1 in m[1], dna in DNA(m, a1-15, a1+15) where dna"
			  ) ;
    }
  else
    {
      qq = hprintf (h, "%s %s %s %s"
		    , "Find variant substitution && ! multi_substitution && Found_in_genome" 
		    , tsnp->filter ? " && (" : ""
		    , tsnp->filter ? tsnp->filter : ""
		    , tsnp->filter ? ") " : ""
		    ) ;

      bqlQuery = hprintf (h, "%s %s "
			  , " select v,typ,dna    "
			  , "from v in @, typ in v->typ, chr in v->intMap, a1 in chr[1], seq in v->Parent_sequence, dna in DNA(seq, a1-15, a1+15) where dna"
			  ) ;
    }
  ksIn = query (0, qq) ;
  if (keySetMax (ksIn))   bqlIter = bqlIterCreate (bqlQuery, ksIn, &errors, h) ;
  if (bqlIter) while (ac_free (h1), (tbl = bqlIterTable (bqlIter)))
    {
      KEY v = ac_table_key (tbl, 0, 0, 0) ;
      const char *typ = ac_table_printable (tbl, 0, 1, "xxx") ;
      char cc = typ[2] ;
      const char *dna = ac_table_printable (tbl, 0, 2, 0) ;
      char buf[32] ;

      nn1++ ;
      if (dna && strlen (dna) == 31)
	{
	  aceOutf (ao, "%s\tRef__%s\t%s\n", dna, name(v), zone) ;
	  memcpy (buf, dna, 32) ;
	  buf[15] = cc ; 
	  aceOutf (ao, "%s\tVariant__%s\t%s\n", buf, name(v), zone) ;
	  nn++ ;
	}
      if (nn < 0) break ;
    }

  ac_free (ksIn) ;
  ac_free (h) ;
  fprintf (stderr, "# tsnpMakeWords_Sub isMrna %s exported %d/%d snps\n", isMrna ? "TRUE" : "FALSE", nn, nn1) ;
  return nn ;
} /* tsnpMakeWords_Sub */

/*************************************************************************************/

static int tsnpMakeWords_Del (TSNP *tsnp, ACEOUT ao, BOOL isMrna)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  AC_TABLE tbl = 0 ;
  char *qq ;
  KEYSET ksIn = 0 ;
  const char *zone = tsnp->zone ? tsnp->zone : "" ;
  const char *errors = 0 ;
  char *bqlQuery = 0 ;
  BQL_ITER *bqlIter = 0 ;
  int nn = 0, nn1 = 0 ;

  if (isMrna)
    {
      qq = "Find variant deletion && ! Found_in_genome && mRNA" ;
      bqlQuery = hprintf (h, "%s %s %s %s %s"
			  , " select v,dnaRD,dnaRA,dnaVD,dnaVA    "
			  , "from v in @, m in v->mRNA, a1 in m[1], a2 in m[2], dnaRD in DNA(m, a1-15, a1+15), dnaRA in DNA(m, a2-15, a2+15) where dnaRD, dnaVD in DNA(m, a1-14, a1), dnaVA in DNA(m, a2, a2+15)  where dnaVA"
			  , tsnp->filter ? " and (" : ""
			  , tsnp->filter ? tsnp->filter : ""
			  , tsnp->filter ? ") " : ""
			  ) ;
    }
  else
    {
      qq = "Find variant deletion && Found_in_genome" ;
      bqlQuery = hprintf (h, "%s %s %s %s %s"
			  , " select v,dnaRD,dnaRA,dnaVD,dnaVA    "
			  , "from v in @, chr in v->intMap, a1 in chr[1], a2 in chr[2], seq in v->Parent_sequence, dnaRD in DNA(seq, a1-15, a1+15) where dnaRD, dnaRA in DNA(seq, a2-15, a2+15) , dnaVD in DNA(seq, a1-14, a1), dnaVA in DNA(seq, a2, a2+15)  where dnaVA"
			  , tsnp->filter ? " and (" : ""
			  , tsnp->filter ? tsnp->filter : ""
			  , tsnp->filter ? ") " : ""
			  ) ;
    }
  ksIn = query (0, qq) ;
  if (keySetMax (ksIn)) bqlIter = bqlIterCreate (bqlQuery, ksIn, &errors, h) ;
  if (bqlIter) while (ac_free (h1), (tbl = bqlIterTable (bqlIter)))
    {
      KEY v = ac_table_key (tbl, 0, 0, 0) ;
      const char *dnaRD = ac_table_printable (tbl, 0, 1, "") ;
      const char *dnaRA = ac_table_printable (tbl, 0, 2, "") ;
      const char *dnaVD = ac_table_printable (tbl, 0, 3, "") ;
      const char *dnaVA = ac_table_printable (tbl, 0, 4, "") ;

      nn1++ ;
      if (dnaRD && strlen (dnaRD) == 31)
	{
	  aceOutf (ao, "%s\tRefDonor__%s\t%s\n", dnaRD, name(v), zone) ;
	  aceOutf (ao, "%s\tRefAccep__%s\t%s\n", dnaRA, name(v), zone) ;
	  aceOutf (ao, "%s%s\tVariant__%s\t%s\n", dnaVD,dnaVA,name(v), zone) ;
	  nn++ ;
	}
      if (nn < 0) break ;
    }

  ac_free (ksIn) ;
  ac_free (h) ;
  fprintf (stderr, "# tsnpMakeWords_Del isMrna %s exported %d/%d snps\n", isMrna ? "TRUE" : "FALSE", nn, nn1) ;
  return nn ;
} /* tsnpMakeWords_Del */

/*************************************************************************************/

static int tsnpMakeWords_Ins (TSNP *tsnp, ACEOUT ao, BOOL isMrna)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  AC_TABLE tbl = 0 ;
  char *qq ;
  KEYSET ksIn = 0 ;
  const char *zone = tsnp->zone ? tsnp->zone : "" ;
  const char *errors = 0 ;
  char *bqlQuery = 0 ;
  BQL_ITER *bqlIter = 0 ;
  int nn = 0, nn1 = 0 ;

  if (isMrna)
    {
      qq = "Find variant Insertion && ! Found_in_genome && mRNA" ;
      bqlQuery = hprintf (h, "%s %s %s %s %s"
			  , " select v,dnaR,dnaVD,dnaVA,v2    "
			  , "from v in @, m in v->mRNA, a1 in m[1], a2 in m[2], dnaR in DNA(m, a1-15, a1+15), dnaRD in DNA(m, a1-15, a1) where dnaRD, dnaRA in DNA(m, a1+1, a1+15), v1 in v->VCF[2], v2 in v1[1] where v2"
			  , tsnp->filter ? " and (" : ""
			  , tsnp->filter ? tsnp->filter : ""
			  , tsnp->filter ? ") " : ""
			  ) ;
    }
  else
    {
      qq = "Find variant Insertion && Found_in_genome" ;
      bqlQuery = hprintf (h, "%s %s %s %s %s"
			  , " select v,dnaR,dnaVD,dnaVA,v2    "
			  , "from v in @, chr in v->intMap, a1 in chr[1], a2 in chr[2], seq in v->Parent_sequence, dnaR in DNA(seq, a1-15, a1+15), dnaVD in DNA(seq, a1-15, a1) where dnaVD, dnaVA in DNA(seq, a1+1, a1+15), v1 in v->VCF[2], v2 in v1[1] where v2"
			  , tsnp->filter ? " and (" : ""
			  , tsnp->filter ? tsnp->filter : ""
			  , tsnp->filter ? ") " : ""
			  ) ;
    }
  ksIn = query (0, qq) ;
  if (keySetMax (ksIn))   bqlIter = bqlIterCreate (bqlQuery, ksIn, &errors, h) ;
  if (bqlIter) 
    while (ac_free (h1), (tbl = bqlIterTable (bqlIter)))
      {
	KEY v = ac_table_key (tbl, 0, 0, 0) ;
	const char *dnaR = ac_table_printable (tbl, 0, 1, "") ;
	const char *dnaVD = ac_table_printable (tbl, 0, 2, "") ;
	const char *dnaVA = ac_table_printable (tbl, 0, 3, "") ;
	const char *v2 = ac_table_printable (tbl, 0, 4, "") ;
	char buf[64], *cp ; const char *cq ;
	int i = 0 ;
	
	nn1++ ;
	i = 0 ; sscanf (v2, "%d", &i) ; if (i>0) continue ; /* we scanned a number */ 
	if (v2 && v2[0]) v2++ ; /* remove the hook base */
	if (dnaR && strlen (dnaR) == 31 && v2 && v2[0])
	  {
	    i = 0 ;
	    cp = buf ; cq = dnaVD ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; } 
	    cq = v2 ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; } 
	    cq = dnaVA ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; }
	    *cp = 0 ; buf[32] = 0 ;
	    if (strlen(buf)==31)
	      {
		aceOutf (ao, "%s\tRef__%s\t%s\n", dnaR, name(v), zone) ;
		aceOutf (ao, "%s\tVariant__%s\t%s\n", buf,name(v), zone) ;
		nn++ ;
	      }
	  }
	if (nn < 0) break ;
      }
  
  ac_free (ksIn) ;
  ac_free (h) ;
  fprintf (stderr, "# tsnpMakeWords_Ins isMrna %s exported %d/%d snps\n", isMrna ? "TRUE" : "FALSE", nn, nn1) ;
  return nn ;
} /* tsnpMakeWords_Ins */

/*************************************************************************************/

static int tsnpMakeWords_MultiSub (TSNP *tsnp, ACEOUT ao, BOOL isMrna)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  AC_TABLE tbl = 0 ;
  char *qq ;
  KEYSET ksIn = 0 ;
  const char *zone = tsnp->zone ? tsnp->zone : "" ;
  const char *errors = 0 ;
  char *bqlQuery = 0 ;
  BQL_ITER *bqlIter = 0 ;
  int nn = 0, nn1 = 0 ;

  if (isMrna)
    {
      qq = "Find variant (Multi_substitution || DelIns) && ! Found_in_genome && mRNA" ;
      bqlQuery = hprintf (h, "%s %s %s %s %s"
			  , " select v,x1,x2,,y,dnaRD,dnaRA,v1,v2    "
			  , "from v in @, x1 in v->DelIns, x2 in x1[1], y in v->Multi_substitution, m in v->mRNA, a1 in m[1], a2 in m[2], dnaR in DNA(m, a1-15, a1+15), dnaRD in DNA(m, a1-15, a1-1) where dnaRD, dnaRA in DNA(m, a2, a2+15), v1 in v->VCF[2], v2 in v1[1] where v2 "
			  , tsnp->filter ? " and (" : ""
			  , tsnp->filter ? tsnp->filter : ""
			  , tsnp->filter ? ") " : ""
			  ) ;
    }
  else
    {
      qq = "Find variant (Multi_substitution || DelIns)  && Found_in_genome" ;
      bqlQuery = hprintf (h, "%s %s %s %s %s"
			  , " select v,x1,x2,y,dnaRD,dnaRA,v1,v2    "
			  , "from v in @, x1 in v->DelIns, x2 in x1[1], y in v->Multi_substitution, chr in v->intMap, a1 in chr[1], a2 in chr[2], seq in v->Parent_sequence, dnaRD in DNA(seq, a1-15, a1-1), dnaRA in DNA(seq, a2, a2+15) where dnaRA, v1 in v->VCF[2], v2 in v1[1] where v2 "
			  , tsnp->filter ? " and (" : ""
			  , tsnp->filter ? tsnp->filter : ""
			  , tsnp->filter ? ") " : ""
			  ) ;
    }
  ksIn = query (0, qq) ;
  if (keySetMax (ksIn)) bqlIter = bqlIterCreate (bqlQuery, ksIn, &errors, h) ;
  if (bqlIter) 
    while (ac_free (h1), (tbl = bqlIterTable (bqlIter)))
      {
	KEY v = ac_table_key (tbl, 0, 0, 0) ;
	int x1 = ac_table_int (tbl, 0, 1, 0) ;
	int x2 = ac_table_int (tbl, 0, 2, 0) ;
	int y = ac_table_int (tbl, 0, 3, 0) ;
	
	const char *dnaRD = ac_table_printable (tbl, 0, 4, "") ;
	const char *dnaRA = ac_table_printable (tbl, 0, 5, "") ;
	const char *v1, *v2 ;
	char bufR[64], bufV[64], *cp ; const char *cq ;
	int i ;
	nn1++ ;
	
	v1 = ac_table_printable (tbl, 0, 6, "") ;
	i = 0 ; sscanf (v1, "%d", &i) ; if (i>0) continue ; /* we scanned a number */ 
	if (x1 && strlen (v1) != x1) continue ;
	if (y && strlen (v1) != y) continue ;

	v2 = ac_table_printable (tbl, 0, 7, "") ; /* v1, v2 are not stable, this is a bug of the library that must be fixed */
	i = 0 ; sscanf (v2, "%d", &i) ; if (i>0) continue ; /* we scanned a number */ 
	if (x2 && strlen (v2) != x2) continue ;
	if (y && strlen (v2) != y) continue ;

	if (1) 
	  {
	    i = 0 ;
	    cp = bufR ; cq = dnaRD ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; } 
	    v1 = ac_table_printable (tbl, 0, 6, "") ;
	    cq = v1 ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; } 
	    cq = dnaRA ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; }
	    *cp = 0 ; bufR[32] = 0 ;
	    
	    i = 0 ;
	    cp = bufV ; cq = dnaRD ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; } 
	    v2 = ac_table_printable (tbl, 0, 7, "") ;
	    cq = v2 ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; } 
	    cq = dnaRA ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; }
	    *cp = 0 ; bufV[32] = 0 ;
	    if (strlen(bufR)==31 && strlen(bufV)==31)
	      {
		aceOutf (ao, "%s\tRef__%s\t%s\n", bufR, name(v), zone) ;
		aceOutf (ao, "%s\tVariant__%s\t%s\n", bufV,name(v), zone) ;
		nn++ ;
	      }
	  }
	if (nn < 0) break ;
      }
  
  ac_free (ksIn) ;
  ac_free (h) ;
  fprintf (stderr, "# tsnpMakeWords_Multi isMrna %s exported %d/%d snps\n", isMrna ? "TRUE" : "FALSE", nn, nn1) ;
  return nn ;
} /* tsnpMakeWords_MultiSub */

/*************************************************************************************/

static void tsnpMakeWords (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (tsnp->outFileName, ".w31", tsnp->gzo, h) ;

  tsnpMakeWords_Any (tsnp, ao, TRUE) ;
  tsnpMakeWords_Any (tsnp, ao, FALSE) ;
  if (0)
    {
      tsnpMakeWords_Sub (tsnp, ao, TRUE) ;
      tsnpMakeWords_Sub (tsnp, ao, FALSE) ;
      tsnpMakeWords_Del (tsnp, ao, TRUE) ;
      tsnpMakeWords_Del (tsnp, ao, FALSE) ;
      tsnpMakeWords_Ins (tsnp, ao, TRUE) ; 
      tsnpMakeWords_Ins (tsnp, ao, FALSE) ;
      tsnpMakeWords_MultiSub (tsnp, ao, TRUE) ;  
      tsnpMakeWords_MultiSub (tsnp, ao, FALSE) ;  
    }
  ac_free (h) ;
} /* tsnpMakeWords */

/*************************************************************************************/
/*************************************************************************************/
/* compute for jafar sum of Nth pwers of the integers in formal variables */
static void getSumNthPowers (int N)
{
  int ii, P = 8 , s[P+1], a[N +1], b[N+1], d = 1 ;
  int M = 5 ;
  /* compute the actual partial sums up to order P */
  for (ii = 0 ; ii <= 8 ; ii++)
    {
      s[ii] = 0 ;
      if (ii > 0)
	{
	  int k = ii ;
	  for (int j=2 ; j <= N ; j++)
	    k *= ii ;
	  s[ii] = s[ii-1] + k ;
	}
    }
  switch (N)
    {
    default: return ;
    case 1:  /* n (n+1)/2 */
      for (a[0] = 1 ; a[0] <= M ; a[0]++)
	for (a[1] = 1 ; a[1] <= M ; a[1]++)
	  for (b[0] = -M ; b[0] <= M ; b[0]++)
	    for (b[1] = -M ; b[1] <= M ; b[1]++)
	      {
		int jj ;

		d = 2 * a[0] * a[1] ;
		if (d == 0)
		  continue ;
		for (jj = 0 ; jj < 8 ; jj++)
		  {
		    int x = (a[0] * jj + b[0]) * (a[1] * jj + b[1]) ;
		    if (x % d)
		      break ;
		    if (x != d*s[jj])
		      break ;
		  }
		if (jj == 8) /* success */
		  {
		    printf("+++++  formula is sum(i^N) = (%d n + %d) (%d n + %d)/%d\n", a[0], b[0], a[1], b[1], d) ;
		    exit (0) ;
		  }
	      }
      break ;
    case 2: /* n (n+1) (2n+1)/6 */
      for (a[0] = 1 ; a[0] <= M ; a[0]++)
	for (a[1] = 1 ; a[1] <= M ; a[1]++)
	  for (a[2] = 1 ; a[2] <= M ; a[2]++)
	    for (b[0] = -M ; b[0] <= M ; b[0]++)
	      for (b[1] = -M ; b[1] <= M ; b[1]++)
		for (b[2] = -M ; b[2] <= M ; b[2]++)
	      {
		int jj ;

		d = 3 * a[0] * a[1] * a[2] ;
		if (d == 0)
		  continue ;
		for (jj = 0 ; jj < 8 ; jj++)
		  {
		    int x = (a[0] * jj + b[0]) * (a[1] * jj + b[1]) * (a[2] * jj + b[2]) ;
		    if (x % d)
		      break ;
		    if (x != d*s[jj])
		      break ;
		  }
		if (jj == 8) /* success */
		  {
		    printf("+++++  formula is sum(i^N) = (%d n + %d) (%d n + %d) (%d n + %d)/%d\n", a[0], b[0], a[1], b[1], a[2], b[2], d) ;
		    exit (0) ;
		  }
	      }
      break ;
    case 3: /* n^2 (n+1)^2/4 */
      for (a[0] = 1 ; a[0] <= M ; a[0]++)
	for (a[1] = 1 ; a[1] <= M ; a[1]++)
	  for (a[2] = 1 ; a[2] <= M ; a[2]++)
	    for (a[3] = 1 ; a[3] <= M ; a[3]++)
	    for (b[0] = -M ; b[0] <= M ; b[0]++)
	      for (b[1] = -M ; b[1] <= M ; b[1]++)
		for (b[2] = -M ; b[2] <= M ; b[2]++)
		  for (b[3] = -M ; b[3] <= M ; b[3]++)
	      {
		int jj ;

		d = N + 1 ;
		for (int i = 0 ; i <= N ; i++)
		  d = d * a[i] ;
		if (d == 0)
		  continue ;
		for (jj = 0 ; jj < 8 ; jj++)
		  {
		    int x = 1 ;
		    for (int i = 0 ; i <= N ; i++)
		      x *= a[i] * jj + b[i] ;
		    if (x % d)
		      break ;
		    if (x != d * s[jj])
		      break ;
		  }
		if (jj == 8) /* success */
		  {
		    printf("+++++  formula is sum(i^N) = ") ;
		    for (int i = 0 ; i <= N ; i++)
		      printf (" (%d n + %d) ", a[i], b[i]) ;
		    printf (" / %d\n", d) ;
		    exit (0) ;
		  }
	      }
      break ;
    case 4:
      /* faux la vrai formule est n (n+1)(2n+1)(3n^2 + 3n -1)/30 */
      for (a[0] = 1 ; a[0] <= M ; a[0]++)
	for (a[1] = 1 ; a[1] <= M ; a[1]++)
	  for (a[2] = 1 ; a[2] <= M ; a[2]++)
	    for (a[3] = 1 ; a[3] <= M ; a[3]++)
	      for (a[4] = 1 ; a[4] <= M ; a[4]++)
		for (b[0] = -M ; b[0] <= M ; b[0]++)
		  for (b[1] = -M ; b[1] <= M ; b[1]++)
		    for (b[2] = -M ; b[2] <= M ; b[2]++)
		      for (b[3] = -M ; b[3] <= M ; b[3]++)
			for (b[4] = -M ; b[4] <= M ; b[4]++)
		      {
		int jj ;

		d = N + 1 ;
		for (int i = 0 ; i <= N ; i++)
		  d = d * a[i] ;
		if (d == 0)
		  continue ;
		for (jj = 0 ; jj < 8 ; jj++)
		  {
		    int x = 1 ;
		    for (int i = 0 ; i <= N ; i++)
		      x *= a[i] * jj + b[i] ;
		    if (x % d)
		      break ;
		    if (x != d * s[jj])
		      break ;
		  }
		if (jj == 8) /* success */
		  {
		    printf("+++++  formula is sum(i^N) = ") ;
		    for (int i = 0 ; i <= N ; i++)
		      printf (" (%d n + %d) ", a[i], b[i]) ;
		    printf (" / %d\n", d) ;
		    exit (0) ;
		  }
	      }
      break ;
    }
  printf ("Solution NOT found N = %d\n", N) ;
		
}

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/
/* run -run NA12878MOD -laneList tmp/TSNP/NA12878MOD/LaneList -t 20 -target_fasta TARGET/CHROMS/hs.chrom_20.fasta.gz -t1 25000001 -t2 30010000 -minSnpFrequency 18 -minSnpCover 10 -minSnpCount 4 -target_class Z_genome -o tata -maxLanes 4
 */

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  if (argc == 0)
  fprintf (stderr,
	   "// Usage: tsnp -db $DB -project $MAGIC -db_report -i[-f] f.val.txt \n"
	   "//      Post processing of .val.txt SNP counts exported by tricoteur -val\n"
	   "//      try: -h --help \n"
	   "// Example\n"
	   "//     tsnp -i f1.val.txt,f2.val.txt -db tmp/TSNP/$zone -p MMM -db_report\n"
	   "// Several sets of parameters must or may be specified\n"
	   "// SNP COUNTS TO BE ANALYZED: \n"
	   "//   -i f[,g[,h]] : coma separated list of files\n"
	   "//   -f file : file contains a list of files, coma or space or tab or line separated\n"
	   "//      ATTENTION all file names must be local or fully qualified starting at /\n"
	   "//      All the input files are of type .val.txt as exported by tricoteur -val\n"
	   "//      The format is very specific and cannot easily be replicated\n"
	   "//      So the present program can only postprocess tricoteur\n"
	   "// TARGET ZONE:\n"
	   "//   --target_class : target class (A_mito ... Z_genome) [default Z_genome]\n"
	   "//   -t targetName :  [optional] often a chromosome name\n"
	   "//   --t1 int --t2 int :  [optional] analyze only section [t1,t2] of the target\n"
	   "//    example:\n"
	   "//        --target_class  Z_genome -t chr7 --t1 2000000 --t2 3000000\n"
	   "//       Typically t2 - t1 = 10 Mbases, or if t1, t2 are not specified\n"
	   "//       the whole chromosome is analyzed .\n"
	   "// SNP FILTERS, DISCOVERY PHASE\n"
	   "//   --method arbitrary_name: the method is echoed in the output file\n" 
  	   "//   --minSnpCover integer : min coverage [default 10] \n"
	   "//   --minSnpCount integer: [default 4]\n"
	   "//   --minSnpFrequency float: [default 18] minimal MAF precentage\n"
	   "//   --intron : [default off] diifferentiate introns from deletions\n"
	   "//   --intron_only : just detect and report introns\n"
	   "//   --min_intron <int> : [default 30] min reported intron length\n"
	   "//   --max_intron <int> : [default 0]  max reported intron length\n"
	   "// Phase 3 actions: Analyse the .snp files\n"
	   "//   The analyses rely on a the existence of an acedb TSNP_DB/$zone database, aware of genes, transcripts and coding structures\n"
	   "//     and containing a copy of the metadata of the runs, originally hand constructed in MetaDB\n"
	   "//     The parameter --db points to this database, we recommend one database per zone, allowing parallelization\n"
	   "//   --db_remap2genome tmp/METADATA/mrnaRemap.gz  --db ACEDB [--force]\n"
	   "//      Remap the transcript variants into genome coordinates\n"
	   "//      --force : remap all variants, default: remmap on those lacking the IntMap tag\n"
	   "//   --db_remap2genes tmp/METADATA/mrnaRemap.gz --db ACEDB [--force]\n"
	   "//      Remap the genome variants into transcript coordinates\n"
	   "//      --force : remap all variants, default: remmap on those lacking the GeneBox tag\n"
	   "//   --db_translate --db ACEDB [--select snp_name]\n"
	   "//      Translate all [matching] mRNA variants (or genome variants remapped to mRNAs) if they map to a protein coding exon\n"
	   "//   --db_count_types  --db ACEDB  : count per run the number of SNP passing the cover/count/frequency thresholds of each type sub/indels/transition ...\n"
	   "//   --db_group_count  --db ACEDB  : add up the counts per group\n"
	   "// GENE FUSION\n"
	   "//   --target_class : target class (KT_RefSeq ET_av...) [default ET_av]\n"
	   "//   --min_GF integer : [default 5]  filter geneFusionFile on min support \n" 
	   "//   --minOverhang integer : [default 15] minimal number of bases\n"
	   "//   --geneFusion fileName: file of genefusions to be analysed\n"
	   "//      mrna1 a1 a2 mrna2 b1 b2 n (n supports for a jump from mrna1[position a2] to m2[b1]\n"
	   "//      Scan the hit file(s) report for each donor/acceptor read count that support\n"
	   "//         the proposed donor and goes to the acceptor\n"
	   "//         OR align locally OR jump locally OR jump elsewhere\n"
	   "// OUTPUT\n"
	   "//     The program exports a vcf file and several histograms\n"
	   "//   -o fileNamePrefix : output file name, equivalent to redirecting stdout\n"
  	   "//   --gzo : the output files should be gzipped\n"
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
  TSNP tsnp ;
  AC_HANDLE h = 0 ;
  char commandBuf [4000] ;

  freeinit () ; 
  messErrorInit (argv[0]) ;

  h = ac_new_handle () ;
  memset (&tsnp, 0, sizeof (TSNP)) ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
  tsnp.h = h ;

  if (argc < 2)
    usage (commandBuf, argc, argv) ;
  if (getCmdLineBool (&argc, argv, "-help") ||
      getCmdLineBool (&argc, argv, "--help")||
      getCmdLineBool (&argc, argv, "-h")
      )
    usage (commandBuf, 0, argv) ;
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

  if (0)
    {
      int N = 0 ;
      if (getCmdLineInt (&argc, argv, "--jaf", &N))
	getSumNthPowers (N) ;
      exit(0) ;
    }
  /* consume optional args */
  tsnp.gzi = getCmdLineBool (&argc, argv, "--gzi") ;
  tsnp.gzo = getCmdLineBool (&argc, argv, "--gzo") ;

  getCmdLineOption (&argc, argv, "-o", &(tsnp.outFileName)) ;
 
  /* TARGET ZONE */
  tsnp.target_class = 0 ; /* "Z_genome" ; */
  getCmdLineOption (&argc, argv, "--referenceGenome", &(tsnp.referenceGenome)) ;
  getCmdLineOption (&argc, argv, "--target_class", &(tsnp.target_class)) ;
  getCmdLineOption (&argc, argv, "-t", &(tsnp.target)) ;
  getCmdLineInt (&argc, argv, "--t1", &(tsnp.t1)) ;
  getCmdLineInt (&argc, argv, "--t2", &(tsnp.t2)) ;

  tsnp.intron_only = getCmdLineBool (&argc, argv, "--intron_only") ;

  tsnp.minIntron = 30 ; /* default */
  tsnp.maxIntron = 0 ; /* default : do not search introns */
  if (getCmdLineBool (&argc, argv, "--intron"))
    tsnp.maxIntron = 100000 ; /* default if we search introns */
 
  if (tsnp.t2 < tsnp.t1)
    {
      fprintf (stderr, "FATAL ERROR: zone limits expeted t1 < t2, got %d > %d, please try\n\tvariant_caller -help\n", tsnp.t1, tsnp.t2) ;
      exit (1) ;
    }
  if (tsnp.t1 < 0)
    {
      fprintf (stderr, "FATAL ERROR: zone limits expeted t1 > 0 got %d  please try\n\tvariant_caller -help\n", tsnp.t1) ;
      exit (1) ;
    }

  /* SNP FILTERS */
  tsnp.minSnpCover = 20 ;
  tsnp.minSnpCount = 4;
  tsnp.minSnpFrequency = 5 ;

  getCmdLineInt (&argc, argv, "--minSnpCover", &(tsnp.minSnpCover)) ;
  getCmdLineInt (&argc, argv, "--minSnpCount", &(tsnp.minSnpCount)) ;
  getCmdLineFloat (&argc, argv, "--minSnpFrequency", &(tsnp.minSnpFrequency)) ;


  getCmdLineOption (&argc, argv, "-i", &(tsnp.inFileList)) ;
  getCmdLineOption (&argc, argv, "--inFiles", &(tsnp.inFileList)) ;
  getCmdLineOption (&argc, argv, "-f", &(tsnp.inFileOfFileList)) ;
  getCmdLineOption (&argc, argv, "--fileList", &(tsnp.inFileOfFileList)) ;
  getCmdLineOption (&argc, argv, "--project", &(tsnp.project)) ;
  getCmdLineOption (&argc, argv, "--select", &(tsnp.select)) ;
  getCmdLineOption (&argc, argv, "-p", &(tsnp.project)) ;
  getCmdLineOption (&argc, argv, "--zone", &(tsnp.zone)) ;
  getCmdLineOption (&argc, argv, "--filter", &(tsnp.filter)) ;
  getCmdLineOption (&argc, argv, "--db_remap2genome", &tsnp.remap2genome) ;
  getCmdLineOption (&argc, argv, "--db_remap2genes", &tsnp.remap2genes) ;
  tsnp.countTypes = getCmdLineBool (&argc, argv, "--db_count_types") ;
  tsnp.force = getCmdLineBool (&argc, argv, "--force") ;

  tsnp.makeWords = getCmdLineBool (&argc, argv, "--makeWords") ;

  tsnp.nAna = 4 ;
  getCmdLineInt (&argc, argv, "--nAna", &(tsnp.nAna)) ;
  
  
  /********* REPORT **********/

  getCmdLineOption (&argc, argv, "-db", &(tsnp.dbName)) ;
  getCmdLineOption (&argc, argv, "--db", &(tsnp.dbName)) ;
  getCmdLineOption (&argc, argv, "--wiggleDir", &(tsnp.wiggleDir)) ;
  tsnp.max_threads = 4 ;
  getCmdLineInt (&argc, argv, "--max_threads", &tsnp.max_threads) ;
  
  if (tsnp.minSnpCount > tsnp.minSnpCover)
    tsnp.minSnpCount =  tsnp.minSnpCover ;

  tsnp.mergeCounts = getCmdLineBool (&argc, argv, "--merge") ;
  tsnp.dbTranslate = getCmdLineBool (&argc, argv, "--db_translate") ;
  tsnp.dbGroupCount = getCmdLineBool (&argc, argv, "--db_group_count") ;
  tsnp.dbGGG = getCmdLineBool (&argc, argv, "--db_GGG") ;
  tsnp.dropMonomodal = getCmdLineBool (&argc, argv, "--dropMonomodal") ;

  tsnp.dbReport = getCmdLineBool (&argc, argv, "--db_report") ;
  tsnp.dbSnpProfile = getCmdLineBool (&argc, argv, "--db_snp_profile") ;

  if (tsnp.dbName)
    {
      const char *errors ;
      tsnp.db = ac_open_db (tsnp.dbName, &errors);
      if (! tsnp.db)
	messcrash ("Failed to open db %s, error %s", tsnp.dbName, errors) ;
    }
  if (argc != 1)
    {
      fprintf (stderr, "unknown argument, sorry\n") ;
      usage (commandBuf, argc, argv) ;
    }
  filAddDir ("./") ;

  if (tsnp.dbReport || tsnp.makeWords)
    {
      if (! tsnp.db)
	{
	  fprintf (stderr, "-db_report requires -db SnpMetaDataDB, sorry, try -help\n") ;
	  exit (1) ;
	}
      if (! tsnp.project)
	{
	  fprintf (stderr, "-db_report requires -project $MAGIC, sorry, try -help\n") ;
	  exit (1) ;
	}
      if (! tsnpGetRunList (&tsnp))
	messcrash ("No run in -db %s belong to project %s", tsnp.dbName, tsnp.project) ;
    }
  if (tsnp.dbTranslate)
    {
      if (0 && ! tsnp.project)
	{
	  fprintf (stderr, "-db_translate requires -project $MAGIC, sorry, try -help\n") ;
	  exit (1) ;
	}
      if (! tsnp.db)
	{
	  fprintf (stderr, "-db_translate requires -db SnpMetaDataDB, sorry, try -help\n") ;
	  exit (1) ;
	}
    }
  if (tsnp.dbGroupCount)
    {
      if (! tsnp.project)
	{
	  fprintf (stderr, "-db_group_count requires -project $MAGIC, sorry, try -help\n") ;
	  exit (1) ;
	}
      if (! tsnp.db)
	{
	  fprintf (stderr, "-db_group_count requires -db SnpMetaDataDB, sorry, try -help\n") ;
	  exit (1) ;
	}
    }
  if (tsnp.dbGGG)
    {
      if (! tsnp.db)
	{
	  fprintf (stderr, "-db_GGG -db SnpMetaDataDB, sorry, try -help\n") ;
	  exit (1) ;
	}
    }

  /* parallelization */
   if (tsnp.max_threads < 4)
     tsnp.max_threads = 4 ;
   if (tsnp.max_threads < tsnp.nAna)
     tsnp.max_threads = tsnp.nAna ;
 
  /* check the absolute args */

  wego_max_threads (tsnp.max_threads) ;
  tsnp.doneChan = channelCreate (12, int, tsnp.h) ;

  if (tsnp.makeWords)
    {
      tsnpMakeWords (&tsnp) ;
    }
  if (tsnp.dbReport)
    {
      tsnpDbReport (&tsnp) ;
    }
  if (tsnp.mergeCounts)
    {
      tsnpMergeCounts (&tsnp) ;
    }
  if (tsnp.dbGGG)
    {
      if (1) tsnpExportProfile (&tsnp) ; /* export tsf file for this section */
      if (1) tsnpGGG (&tsnp) ; /* export tsf file for this section */
    }
  if (tsnp.remap2genome)
    {
      if (1) tsnpRemap0 (&tsnp) ; 
      if (1) tsnpCreateAtlas (&tsnp) ;
      if (1) tsnpRemap1 (&tsnp) ; 
    }
  if (tsnp.remap2genes)
    {
      if (1) tsnpRemap0 (&tsnp) ;
      if (1) tsnpCreateAtlas (&tsnp) ;
      if (1) tsnpRemap2 (&tsnp) ;    
    }
  if (tsnp.dbTranslate)
    {
      if (1) tsnpDbTranslate (&tsnp) ;
      if (0) tsnpCodingModif (&tsnp) ; /* probably obsolete */
      if (0) tsnpExportProfile (&tsnp) ; /* export tsf file for this section */
    }
  if (tsnp.dbGroupCount)
    {
      if (1) tsnpDbGroupCount (&tsnp) ;
    }
  if (tsnp.dbSnpProfile)
    {
      if (! tsnp.db)
	{
	  fprintf (stderr, "-db_snp_profile requires -db SnpMetaDataDB, sorry, try -help\n") ;
	  exit (1) ;
	}
      if (! tsnp.project)
	{
	  fprintf (stderr, "-db_snp_profile requires -project $MAGIC, sorry, try -help\n") ;
	  exit (1) ;
	}
      if (tsnpDbRun2groups (&tsnp))
	tsnpDbSnpProfile (&tsnp) ;
    }
  wego_flush () ;

  ac_free (tsnp.ao) ;
  if (1)
    { 
      int mx ;
      messAllocMaxStatus (&mx) ; 
      fprintf (stderr, "// minSnpCount %d minSnpCover %d minSnpFrequency %.1f%%\n"
	       , tsnp.minSnpCount
	       , tsnp.minSnpCover
	       , tsnp.minSnpFrequency
	       ) ;
      if (tsnp.snps)
	fprintf (stderr, "// SNP detected %d SEG reported %d\n"
		 , tsnp.snpDetected
		 , tsnp.snpExported
		 ) ;
     }
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderrn and to ensure all pipes are closed*/
  ac_free (tsnp.h) ;

  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
