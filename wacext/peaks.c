/*  File: peaks.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2026
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

/*
 * Ananlysis of the peaks visible in the wiggle coverage plots
 * from the genome data
 */

#define ARRAY_CHECK 

#include "ac.h"
#include "bitset.h"
#include "peaks.h"
#include "wiggle.h"
#include "wig.h"

typedef struct pkStruct {
  AC_HANDLE h ;
  AC_DB db ;
  Array runs ;
  DICT *chromDict ;

  int step ;
  
  const char *dbName ;
  const char *project ;
  const char *group ;
  const char *outFileName ;

  /* accepted on the command line so that existing scripts keep working, but not used (yet) */
  const char *orderBy ;
  const char *inFileName ;
  int min_support ;

  int min_peak_cover ;

  BOOL intersects ;
  BOOL makeGroupWiggles ;
  BOOL getStrandShifts ; 
  BOOL makePeaks ;
  BOOL makeShiftedNonStrandedPeaks ;
  BOOL gzo ;  
} PK ;

typedef struct runStruct {
  AC_HANDLE h ;
  AC_OBJ run ;
  int type ; /* 0: run, 1:group. 2:reference */
  int input ; /* 1 : is input (control) run */
  const char *name ;    /* run name, "Reference" for the reference run, run is 0 then */
  Array peaks ;
  int minCover ;
  int medianNoZero ;
  float alignedBases ;
} RPK ;

typedef struct peaksStruct {
  int chrom, fr, a1, a2, width, contact ;
  float height, avg, area ;
} PEAK ;

/**************************************************************/

static int peakOrder (const void *a, const void *b)
{
  const PEAK *up = a, *vp = b ;
  int z ;

  z = up->chrom - vp->chrom ; if (z) return z ;
  z = up->a1 - vp->a1 ; if (z) return z ;
  z = up->a2 - vp->a2 ; if (z) return z ;

  return 0 ;
} /* peakOrder*/

/***************************************************************/

/* merge the overlapping or adjacent peaks of an array sorted by peakOrder
 * returns the number of peaks left, the array is shortened accordingly */
static int peakMerge (Array peaks)
{
  int ii, nn, iMax =  arrayMax (peaks) ;
  PEAK *up, *vp ;
  
  if (iMax < 2)
    return iMax ;

  /* vp is the last merged peak, up scans the input starting at the second peak */
  for (ii = 1, nn = 1, vp = arrp (peaks, 0, PEAK), up = vp + 1 ; ii < iMax ; ii++, up++)
    {
      if (up->chrom > vp->chrom || up->a1 > vp->a2 + 1)
	{
	  vp++ ; nn++ ; *vp = *up ;
	}
      else  if (vp->a2 < up->a2)
	vp->a2 = up->a2 ;
    }
  
  arrayMax (peaks) = nn ;
  return nn ;
} /* peakMerge */

/***************************************************************/

static int pkCountOneIntersect (Array peaks1, Array peaks2, int type)
{
  int ii, jj, iMax =  arrayMax (peaks1), jMax =  arrayMax (peaks2) ;
  PEAK *up, *vp ;
  int nn = 0 ;
  
  up = arrp (peaks1, 0, PEAK) ; ii = 0 ;
  vp = arrp (peaks2, 0, PEAK) ; jj = 0 ;
  while (ii < iMax && jj < jMax)
    {
      if (up->chrom < vp->chrom) { up++ ; ii++ ; continue ;}
      if (up->chrom > vp->chrom) { vp++ ; jj++ ; continue ;}
      if (up->a2 < vp->a1) { up++ ; ii++ ; continue ;}
      if (up->a1 > vp->a2) { vp++ ; jj++ ; continue ;}
      nn++ ;
      /* these are the true number of intersects 
      if (up->a2 > vp->a2) { vp++ ; jj++ ; continue ;}
      if (up->a2 <= vp->a2) { up++ ; ii++ ; continue ;}
      */
      /* but we want to be assymetric */
      if (! type) up->contact++ ;
      up++ ; ii++ ;
    }
  return nn ;
} /* pkCountOneIntersect */

/****************************************************************/

static int pkGetRuns (PK *pk)
{
  AC_OBJ run = 0 ;
  int nr = 0 ;
  
  pk->runs = arrayHandleCreate (256, RPK, pk->h) ;
  
  if (pk->project)
    {
      const char *errors = 0 ;
      int ir ;
      AC_TABLE tt = 0 ;
      AC_HANDLE h = ac_new_handle () ; 
      char *qq = hprintf (h," select r,u,b from p in class project where p == \"%s\", r in p->run where ! r#Sublibrary_of, u in r#union_of, ali in r->ali, tg in ali->nh_Ali where tg==\"any\", b in tg[5]", pk->project) ;
      
      tt = ac_bql_table (pk->db, qq, 0, 0, &errors, h) ;
      for (ir = 0 ; tt && ir < (tt->rows) ; ir++)
	{
	  run = ac_table_obj (tt, ir, 0, pk->h) ;
	  if (run)
	    {
	      RPK *rc = arrayp (pk->runs, nr++, RPK) ;
	      rc->run = run ;
	      rc->name = ac_name (run) ;   /* lives as long as run, both on pk->h */
	      rc->type = ac_has_tag (run, "Union_of") ? 1 : 0 ;
	      rc->input = ac_has_tag (run, "ChIP_Input") ? 1 : 0 ;
	      rc->alignedBases = ac_table_int (tt, ir, 2, 0) ;
	      rc->h = pk->h ;
	      fprintf (stderr, "%s\t%s\n", ac_table_printable (tt,ir,0,"xx"), ac_table_printable (tt,ir,1,"xx")) ;
	    }
	}
      ac_free (h) ;
    }

  fprintf (stderr, " pkGetRuns got %d runs in project %s\n", nr, pk->project) ;

  return nr ;
} /* pkGetRuns */

/****************************************************************/

static void pkMakeGroupWiggles (PK *pk)
{
  const char *errors = 0 ;
  int ir ;
  AC_TABLE tt = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  const char *qq2 = pk->group ? pk->group : "*" ;

  for (int pass = 0 ; pass < 2 ; pass++)
    {
      AC_HANDLE h1 = 0, h2 = 0 ;
      char *qq ;
      if (pass == 0)
	qq = hprintf (h, "select g,r,sShift from g in ?Run where g->project == \"%s\" and g ~ \"%s\", r in g>>Union_of where r->project == \"%s\", ali in r->ali where ali, sShift in ali->Strand_shift", pk->project, qq2, pk->project) ;
      else
	qq = hprintf (h, "select g,r,sShift from g in ?Run where g->project == \"%s\" and g ~ \"%s\", r in g>>Sublibraries, ali in r->ali where ali, sShift in ali->Strand_shift", pk->project, qq2) ;
      
      
      tt = ac_bql_table (pk->db, qq, 0, 0, &errors, h1) ;
      fprintf (stderr, "pass=%d qq=%s\nMax=%d\n", pass, qq, tt ? tt->rows : -1) ;
      if (!tt || ! tt->rows)
	{ ac_free (h1) ; continue ; }
      for (int iChrom = 1 ; iChrom <= dictMax (pk->chromDict) ; iChrom++)
	{
	  const char *chrom = dictName (pk->chromDict, iChrom) ;
	  for (int ifr = 0 ; ifr < 2 ; ifr++)
	    {
	      const char *fr = ifr == 0 ? "f" : "r" ;
	      for (ir = 0 ; ir < tt->rows ; ir++)
		{
		  KEY gk = ac_table_key (tt, ir, 0, 0) ;
		  BOOL ok = TRUE ;
		  const char *gNam = ac_table_printable (tt, ir, 0, 0) ;
		  char *fNamG = hprintf (h1, "tmp/WIGGLEGROUP/%s/%s.%s.u.%s.az", gNam, gNam, chrom, fr) ;
		  if (filCheckName (fNamG, 0, "r"))
		    {
		      for ( ; ir < tt->rows ; ir++)
			{
			  KEY gk2 = ac_table_key (tt, ir, 0, 0) ;
			  if (gk2 != gk) break ;
			}
		      ir-- ;
		      continue ;
		    }
		  ac_free (h1) ;
		  h1 = ac_new_handle () ;
		  WIG wig = wigCreate (gNam, chrom, pk->step, h1) ;
		  
		  pk->step = wig->step ;
		  for ( ; ir < tt->rows ; ir++)
		    {
		      KEY gk2 = ac_table_key (tt, ir, 0, 0) ;
		      if (gk2 != gk) break ;
		      
		      const char *rNam = ac_table_printable (tt, ir, 1, 0) ;
		      int strand_shift = ac_table_int (tt, ir, 2, 0) ;
		      fprintf (stderr, "%s\t%s\t%s\t%s\tshift=%d\n", gNam, rNam, chrom, fr, strand_shift) ;
		      ac_free (h2) ;
		      h2 = ac_new_handle () ;
		      
		      char *fNam = hprintf (h2, "tmp/SA/%s/wiggles/%s.%s.u.%s.az", rNam, rNam, chrom, fr) ;
		      if (! filCheckName (fNam, 0, "r"))
			{
			  fNam = hprintf (h2, "tmp/WIGGLEGROUP/%s/%s.%s.u.%s.az", rNam, rNam, chrom, fr) ;
			  if (! filCheckName (fNam, 0, "r"))
			    {
			      fprintf (stderr, "Cannot find file %s\n", fNam) ;
			      ok = FALSE ;
			      continue ;
			    }
			}
		      fprintf (stderr, "Found file %s strand_shift=%d\n", fNam, (ifr % 2) == 0 ? strand_shift : -strand_shift) ;
		      wigParseAZ (wig, fNam, (ifr % 2) == 0 ? strand_shift : -strand_shift) ;
		    }
		  ir-- ;
		  
		  ac_free (h2) ;
		  
		  if (ok) /* export the group wiggle */
		    {
		      const char *fNam3 = "tmp/WIGGLEGROUP" ;
		      filCreateDir (fNam3) ;
		      fNam3 = hprintf (h1, "tmp/WIGGLEGROUP/%s", gNam) ;
		      filCreateDir (fNam3) ;
		      fNam3 = hprintf (h1, "tmp/WIGGLEGROUP/%s/%s.%s.u.%s", gNam, gNam, chrom, fr) ;
		      
		      wigExportAZ (wig, fNam3) ;
		      
		      const char *fNam4 = hprintf (h1, "tmp/WIGGLEGROUP/%s/%s.%s.u.%s.peaks", gNam, gNam, chrom, fr) ;
		      wigExportPeaks (wig, fNam4, pk->gzo) ;
		      
		      int minCover = pk->min_peak_cover ; /* let the system choose */
		      
		      wigFindPeaks (wig, minCover) ;
		      wigExportPeaks (wig, fNam4, FALSE) ;
		  }
		  
		} /* ifr */
	    }  /* chrom */
	  ac_free (h1) ;
	}
    }

  ac_free (h) ;
  
  return ;
} /* pkMakeGroupWiggles */

/****************************************************************/

static void pkMakeShiftedNonStrandedPeaks (PK *pk)
{
  const char *errors = 0 ;
  int ir ;
  AC_TABLE tt = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  const char *qq2 = pk->group ? pk->group : "*" ;
  char *qq = hprintf (h, "select r, dx from r in ?Run where r->project == \"%s\" and r ~ \"%s\", ali in r->ali, dx in ali->Strand_shift", pk->project, qq2) ;

  int minCover = pk->min_peak_cover ;

  tt = ac_bql_table (pk->db, qq, 0, 0, &errors, h) ;
  fprintf (stderr, "qq=%s\nMax=%d\n", qq, tt ? tt->rows : -1) ;
  if (!tt)
    { ac_free (h) ; return ; }
  
  for (ir = 0 ; ir < tt->rows ; ir++)
    {
      const char *rNam = ac_table_printable (tt,ir,0,"x") ;
      int dx = ac_table_int (tt, ir, 1, 0) ;
      
      BOOL ok = TRUE ;
      for (int iChrom = 1 ; iChrom <= dictMax (pk->chromDict) ; iChrom++)
	{
	  const char *chrom = dictName (pk->chromDict, iChrom) ;

	  AC_HANDLE h2 = ac_new_handle () ;
	  char *cpF = 0, *cpR = 0 ;

	  cpF = hprintf (h2, "tmp/SA/%s/wiggles/%s.%s.u.f.az", rNam, rNam, chrom) ;
	  cpR = hprintf (h2, "tmp/SA/%s/wiggles/%s.%s.u.r.az", rNam, rNam, chrom) ;
	  ok = filCheckName (cpF, 0, "r") && filCheckName (cpR, 0, "r") ;
	  if (! ok)
	    {
	      cpF = hprintf (h2, "tmp/WIGGLEGROUP/%s/%s.%s.u.f.az", rNam, rNam, chrom) ;
	      cpR = hprintf (h2, "tmp/WIGGLEGROUP/%s/%s.%s.u.r.az", rNam, rNam, chrom) ;
	      ok = filCheckName (cpF, 0, "r") ;
	    }

	  char *fNamPk = hprintf (h2, "tmp/Peaks.%d/%s/%s.%s.u.nsShifted.peaks", minCover, rNam, rNam, chrom) ;

	  if (ok)
	    {
	      WIG wigFR = wigCreate (rNam, chrom, 10, h2) ;
	      if (wigFR)
		{ /* accumulate the shifted data */
		  wigParseAZ (wigFR, cpF, dx/2) ;
		  if (cpR) wigParseAZ (wigFR, cpR, -dx/2) ;
		  wigFindPeaks (wigFR, minCover) ;
		  wigExportPeaks (wigFR, fNamPk, pk->gzo) ;
		}
	    }
	  ac_free (h2) ;
	}
    }
} /* pkMakeShiftedNonStrandedPeaks */

/****************************************************************/

static void pkMakePeaks (PK *pk)
{
  const char *errors = 0 ;
  int ir ;
  AC_TABLE tt = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  const char *qq2 = pk->group ? pk->group : "*" ;
  char *qq = hprintf (h, "select r, dx from r in ?Run where r->project == \"%s\" and r ~ \"%s\", ali in r->ali, dx in ali->Strand_shift", pk->project, qq2) ;
  tt = ac_bql_table (pk->db, qq, 0, 0, &errors, h) ;
  fprintf (stderr, "qq=%s\nMax=%d\n", qq, tt ? tt->rows : -1) ;
  if (!tt)
    { ac_free (h) ; return ; }
  
  for (ir = 0 ; ir < tt->rows ; ir++)
    {
      const char *rNam = ac_table_printable (tt,ir,0,"x") ;
      
      BOOL ok = TRUE ;
      for (int ifr = 0 ; ifr < 2 ; ifr++)
	{
	  char *fr = ifr ? "r" : "f" ;
	  for (int iChrom = 1 ; ok && iChrom <= dictMax (pk->chromDict) ; iChrom++)
	    {
	      const char *chrom = dictName (pk->chromDict, iChrom) ;
	      
	      AC_HANDLE h2 = ac_new_handle () ;
	      char *cpF = hprintf (h2, "tmp/SA/%s/wiggles/%s.%s.u.%s.az", rNam, rNam, chrom, fr) ;

	      char *fNamPk = hprintf (h2, "%s/%s.%s.u.%s.peaks", pk->outFileName, rNam, chrom, fr) ;
	      ok = filCheckName (cpF, 0, "r") ;
	      if (ok)
		{
		  WIG wig = wigCreate (rNam, chrom, 10, h2) ;
		  if (wig)
		    { /* accumulate the shifted data */
		      wigParseAZ (wig, cpF, 0) ;
		      wigFindPeaks (wig, 0) ;
		      wigExportPeaks (wig, fNamPk, pk->gzo) ;
		    }
		}
	      ac_free (h2) ;
	    }
	}
    }
} /* pkMakePeaks */

/****************************************************************/

static void pkGetStrandShifts (PK *pk)
{
  const char *errors = 0 ;
  int ir ;
  AC_TABLE tt = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  const char *qq2 = pk->group ? pk->group : "*" ;
  char *qq = hprintf (h, "select r, isG, isL  from r in ?Run where r->project == \"%s\" and r ~ \"%s\", isG in r#Union_of, isL in r#Sublibraries, ali in r->ali ", pk->project, qq2) ;
  tt = ac_bql_table (pk->db, qq, 0, 0, &errors, h) ;
  fprintf (stderr, "qq=%s\nMax=%d\n", qq, tt ? tt->rows : -1) ;
  if (!tt)
    { ac_free (h) ; return ; }
  vTXT txt = vtxtHandleCreate (h) ;
  
  for (ir = 0 ; ir < tt->rows ; ir++)
    {
      const char *rNam = ac_table_printable (tt,ir,0,"x") ;
      AC_HANDLE h1 = ac_new_handle () ;
      CST cst = cisTransCreate (5010, h1) ;

      const char *isGroup = ac_table_printable (tt,ir,1,0) ;
      const char *isLib = ac_table_printable (tt,ir,2,0) ;
      for (int iChrom = 1 ; iChrom <= dictMax (pk->chromDict) ; iChrom++)
	{
	  const char *chrom = dictName (pk->chromDict, iChrom) ;
	  AC_HANDLE h2 = ac_new_handle () ;
	  char *cpF, *cpR ;
	  if (isGroup || isLib)
	    {
	      cpF = hprintf (h2, "tmp/WIGGLEGROUP/%s/%s.%s.u.f.az", rNam, rNam, chrom) ;
	      cpR = hprintf (h2, "tmp/WIGGLEGROUP/%s/%s.%s.u.r.az", rNam, rNam, chrom) ;
	    }
	  else
	    {
	      cpF = hprintf (h2, "tmp/SA/%s/wiggles/%s.%s.u.f.az", rNam, rNam, chrom) ;
	      cpR = hprintf (h2, "tmp/SA/%s/wiggles/%s.%s.u.r.az", rNam, rNam, chrom) ;
	    }
	  BOOL ok = filCheckName (cpF, 0, "r") && filCheckName (cpR, 0, "r") ;
	  if (ok)
	    {
	      WIG wigF = wigCreate (0, 0, 10, h2) ;
	      WIG wigR = wigCreate (0, 0, 10, h2) ;
	      if (wigF && wigR)
		{
		  wigParseAZ (wigF, cpF, 0) ;
		  wigParseAZ (wigR, cpR, 0) ;
		  cisTransFeed (cst, wigF->step, wigF, wigR) ;
		}
	    }
	  ac_free (h2) ;
	}
      if (cst->nu)
	{
	  cisTransNormalize (cst) ;
	  char *fNam = isGroup ? hprintf (h1, "tmp/WIGGLEGROUP/%s/strand_shift", rNam) : hprintf (h1, "tmp/SA/%s/strand_shift", rNam) ;
	  ACEOUT ao = aceOutCreate (fNam, 0, FALSE, h1) ;
	  cisTransExport (cst, ao) ;
	  vtxtPrintf (txt, "Ali %s\nStrand_shift %d %.3f\n\n"
		      , rNam
		      , cst->bestShift, cst->bestUv
		      ) ;
	  fprintf (stderr, "%s", vtxtPtr (txt)) ;		   
	  
	  const char *fNam2 = hprintf (h1, "tmp/WIGGLEGROUP/%s/strand_shift.ace", rNam) ;
	  ACEOUT ao2 = aceOutCreate (fNam2, 0, pk->gzo, h1) ;
	  aceOut (ao2, vtxtPtr (txt)) ;
	  
	  vtxtClear (txt) ;
	}
      ac_free (h1) ;
    }
  ac_free (h) ;
  
  return ;
} /* pkGetStrandShifts */

/*************************************************************************************/

/* chromosome names come from $chromSetAll, separated by any blanks, tabs or newlines */
static void pkGetChroms (PK *pk)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *chroms = getenv ("chromSetAll") ;

  pk->chromDict = dictHandleCreate (32, pk->h) ;
  if (chroms)
    {
      char *buf = strnew (chroms, h) ;
      char *save = 0 ;
      char *cp ;

      for (cp = strtok_r (buf, " \t\n", &save) ; cp ; cp = strtok_r (0, " \t\n", &save))
	dictAdd (pk->chromDict, cp, 0) ;
    }
  
  ac_free (h) ;
  return ;
} /* pkGetChroms */

/*************************************************************************************/
static BOOL debug1 = FALSE ;

static void pkGetPeaks (PK *pk)
{
  AC_HANDLE h = ac_new_handle () ;
  int runMax = arrayMax (pk->runs) ;
  int chromMax = dictMax (pk->chromDict) ;
  int minCover = pk->min_peak_cover ;
  KEYSET histo = keySetHandleCreate (h) ;
  
  for (int ir = 0 ; ir < runMax ; ir++)
    {
      RPK *rc = arrayp (pk->runs, ir, RPK) ;
      const char *runName = rc->name ;
      int ii = 0 ;
      int minC = 0, nMedian = 0 ;
      PEAK *up ;
      rc->peaks = arrayHandleCreate (256, PEAK, rc->h) ;
      histo = keySetReCreate (histo) ;
      
      for (int chrom = 1 ; chrom <= chromMax ; chrom++)
	for (int fr = 0 ; fr < 1 ; fr++)
	  {
	    if (debug1 && chrom != 3) continue ;
	    AC_HANDLE h1 = ac_new_handle () ;
	    const char *chromName = dictName (pk->chromDict, chrom) ;
	    char *fNam = hprintf (h1, "tmp/SA/%s/wiggles/%s.%s.u.%s.peaks", runName, runName, chromName, fr==0 ? "nsShifted" : "r") ;
	    if (! filCheckName (fNam, 0, "r"))
	      fNam = hprintf (h1, "tmp/Peaks.%d/%s/%s.%s.u.%s.peaks", minCover, runName, runName, chromName, fr==0 ? "nsShifted" : "r") ;
	    if (! filCheckName (fNam, 0, "r"))
	      continue ;
	    ACEIN ai = aceInCreate (fNam, 0, h1) ;
	    if (! ai) continue ;
	    
	    while (aceInCard (ai))
	      {
		// # target	x1	x2	width	height	avg	area kb	
		//    chr5	4540	4666	127	6.0	5	0.6	
		
		char *cp = aceInWord (ai) ;
		if (!cp || !*cp)
		  continue ;
		
		if (!strcmp (cp, "##"))
		  {
		    cp = aceInPos (ai) ;
		    int k = 0, md = 0,  mc = 0 ;
		    k = sscanf (cp, "Median\t%d\tMinCover\t%d", &md, &mc) ;
		    if (k == 2 && md > 0)
		      { keySet (histo, nMedian++) = md ; minC += mc ; }
		    continue ;
		  }
		
		if (*cp == '#')
		  continue ;
		
		up = arrayp (rc->peaks, ii++, PEAK) ;
		up->chrom = chrom ;
		up->fr = fr ;
		aceInStep (ai, '\t') ; aceInNext (ai) ;
		aceInInt (ai, &up->a1) ;
		aceInStep (ai, '\t') ; aceInNext (ai) ;
		aceInInt (ai, &up->a2) ;
		aceInStep (ai, '\t') ; aceInNext (ai) ;
		aceInInt (ai, &up->width) ;
		aceInStep (ai, '\t') ; aceInNext (ai) ;
		aceInFloat (ai, &up->height) ;
		aceInStep (ai, '\t') ; aceInNext (ai) ;
		aceInFloat(ai, &up->avg) ;
		aceInStep (ai, '\t') ; aceInNext (ai) ;
		aceInFloat (ai, &up->area) ;
		if (debug1 && up->a1 > 100000)break ;
	      }
	    ac_free (h1) ;
	  }
      if (rc->peaks && arrayMax (rc->peaks))
	{
	  arraySort (rc->peaks, peakOrder) ;
	  peakMerge (rc->peaks) ;
	  if (nMedian)
	    {
	      keySetSort (histo) ;
	      rc->medianNoZero = keySet (histo, nMedian/2) ;
	      rc->minCover = minC/nMedian ;
	    }
	}
    }
  ac_free (h) ;
} /* pkGetPeaks */

/* Build a piecewise-constant total-coverage profile from a list of
 * (possibly overlapping) peaks.
 *
 * Input : pks, an Array of PEAK {int chrom, a1, a2, cover;}
 *         a1 <= a2, a2 inclusive, cover may be any (nonzero) int.
 * Output: a new Array of PEAK, one struct per maximal run of positions
 *         sharing the same NONZERO summed coverage, sorted by (chrom,a1),
 *         non-overlapping and adjacent runs of equal cover already merged.
 *
 * Method: difference-array sweep over coordinate events.
 *   each peak contributes  +cover at a1  and  -cover at a2+1
 *   sort events by (chrom, pos); sweep; the running sum is the total
 *   coverage on [pos, next_pos).  Emit segments where the sum is != 0.
 *
 * Complexity: O(n log n) time (the sort), O(n) extra memory, n = pMax.
 */

/*************************************************************************************/
/*************************************************************************************/

typedef struct { int chrom, pos, delta ; } CEVENT ;

/* sort events by chromosome, then position. delta order within a tie
 * is irrelevant: all deltas at the same (chrom,pos) are summed before
 * any segment starting at pos is emitted. */
static int ceventOrder (const void *va, const void *vb)
{
  const CEVENT *a = (const CEVENT *)va ;
  const CEVENT *b = (const CEVENT *)vb ;
  if (a->chrom != b->chrom) return a->chrom - b->chrom ;
  return a->pos - b->pos ;
}

/*************************************************************************************/
/* Returns a new Array(PEAK) allocated on handle h; caller owns it. */

static Array pkPeaksCoverage (Array pks, BOOL oneOrZero, int minCover, int maxCover, AC_HANDLE h)
{
  int ii, pMax = arrayMax (pks) ;
  Array events, out ;
  CEVENT *ev ;
  PEAK *up, *vp ;
  int nOut = 0 ;
  out = arrayHandleCreate (pMax ? pMax : 1, PEAK, h) ;
  if (!pMax)
    return out ;

  /* 1. build 2 events per peak. Skip zero-cover or empty peaks. */
  events = arrayCreate (2 * pMax, CEVENT) ;
  for (ii = 0 ; ii < pMax ; ii++)
    {
      up = arrp (pks, ii, PEAK) ;
      if (up->height == 0 || up->a2 < up->a1)
	continue ;
      ev = arrayp (events, arrayMax (events), CEVENT) ;
      ev->chrom = up->chrom ; ev->pos = up->a1 ;     ev->delta =  oneOrZero ? 1 : up->height ;
      ev = arrayp (events, arrayMax (events), CEVENT) ;
      ev->chrom = up->chrom ; ev->pos = up->a2 + 1 ; ev->delta = -(oneOrZero ? 1 : up->height) ;
    }
  
  /* 2. sort events by (chrom, pos). */
  arraySort (events, ceventOrder) ;
  
  /* 3. sweep. Between consecutive distinct positions the running sum is
   *    constant; emit it as a segment when nonzero. Merge is automatic
   *    because we only open a new output struct when cover changes. */
  {
    int eMax = arrayMax (events) ;
    int jj = 0 ;
    int state = 0 ;
    while (jj < eMax)
      {
	int chrom = arrp (events, jj, CEVENT)->chrom ;
	int sum = 0 ;         /* running coverage on this chromosome   */
	// if (chrom != 3) {jj++ ; continue ;}
	/* process every event on this chromosome */

	while (jj < eMax && arrp (events, jj, CEVENT)->chrom == chrom)
	  {
	    int pos = arrp (events, jj, CEVENT)->pos ;
	    // if (pos>100000) {jj++ ; continue ;}
	    /* find the coverage at this position */
	    while (jj < eMax &&
		   arrp (events, jj, CEVENT)->chrom == chrom &&
		   arrp (events, jj, CEVENT)->pos   == pos)
	      {
		sum += arrp (events, jj, CEVENT)->delta ;
		jj++ ;
	      }
	    if (debug1)
	      fprintf (stderr, "jj = %d pos = %d sum=%d\n", jj, pos, sum) ;
	    if (sum >= minCover)
	      {
		if (state == 0 || chrom != vp->chrom)
		  {
		    vp = arrayp (out, nOut++, PEAK) ;
		    vp->chrom = chrom ;
		    vp->a1    = pos ;
		    vp->height = 0 ;
		    state = 1 ;
		  }
		vp->a2 = pos ;   /* inclusive */
		vp->width = vp->a2 - vp->a1 + 1 ;
		vp->height = sum > vp->height ? sum : vp->height ;
		state = 1 ;
	      }
	    else
	      {
		if (state == 1)
		  {
		    state = 0 ;
		    vp->a2    = pos - 1 ;   /* inclusive */
		    vp->width = vp->a2 - vp->a1 + 1 ;
		    vp->height = sum > vp->height ? sum : vp->height ;
		  }
	      }
	  }
	/* by construction sum == 0 here (every + has its matching -),
	 * so there is no dangling run to close. */
      }
  }
  
  arrayDestroy (events) ;
  return out ;
} /* pkPeaksCoverage */

/*************************************************************************************/

static Array pkPeaksClip (Array peaks, int minCover, int maxCover, AC_HANDLE h)
{
  Array new = arrayHandleCopy (peaks, h) ;
  int ii, jj, iMax = peaks ? arrayMax (peaks) : 0 ;
  PEAK *vp, *wp ;
  
  vp = wp = arrayp (new, 0, PEAK) ;
  for (ii = jj = 0 ; ii < iMax ; ii++, vp++)
    {
      if (vp->height >= minCover && vp->height < maxCover)
	{
	  if (jj < ii) *wp = *vp ;
	  jj++ ; wp++ ;
	}
    }
  arrayMax (new) = jj ;
  return new ;
} /* pkPeaksClip */

/*************************************************************************************/
/* Create a wiggle of reference positions
 * each run contributes qualitativaly as one or zero over its own peaks
 * we then look for locations over 60% and call tem reference peaks
 */
static void pkCreateRefPeaks (PK *pk)
{
  AC_HANDLE h = ac_new_handle () ;
  int runMax = arrayMax (pk->runs) ;
  Array newPeaks = arrayHandleCreate (10000, PEAK, h) ;   /* scratch, freed with h */
  Array newInputPeaks = arrayHandleCreate (10000, PEAK, h) ;   /* scratch, freed with h */
  int nn = 0 ;
  int nnInput = 0 ;
  PEAK *up, *wp ;
  /* copy all peaks */
  int minCover = 0 ;
  int minInputCover = 0 ;
  for (int ir = 0 ; ir < runMax ; ir++)
    {
      RPK *rc1 = arrayp (pk->runs, ir, RPK) ;
      Array p1 = rc1->peaks ;
      if (!rc1->type && !rc1->input && p1 && arrayMax (p1))
	{
	  int pMax = arrayMax (p1) ;
	  up = pMax ? arrp (p1, 0, PEAK) : 0 ;
	  array (newPeaks, nn + pMax - 1, PEAK).width = 0 ;
	  wp = arrayp (newPeaks, nn, PEAK) ;
	  memcpy (wp, up, pMax * sizeof (PEAK)) ;
	  nn += pMax ;
	  minCover++ ;
	}
      if (!rc1->type && rc1->input && p1 && arrayMax (p1))
	{
	  int pMax = arrayMax (p1) ;
	  up = pMax ? arrp (p1, 0, PEAK) : 0 ;
	  array (newInputPeaks, nnInput + pMax - 1, PEAK).width = 0 ;
	  wp = arrayp (newInputPeaks, nnInput, PEAK) ;
	  memcpy (wp, up, pMax * sizeof (PEAK)) ;
	  nnInput += pMax ;
	  minInputCover++ ;
	}
    }
  arraySort (newPeaks, peakOrder) ;
  arraySort (newInputPeaks, peakOrder) ;
  Array newQPeaks4 = pkPeaksCoverage (newPeaks, TRUE, 3, minCover+1, pk->h) ;   /* kept by the reference run */

  Array newQPeaks1 = pkPeaksClip (newQPeaks4, (3.0*minCover)/4, minCover + 1, pk->h) ;
  Array newQPeaks2 = pkPeaksClip (newQPeaks4, minCover/2, (3.0*minCover)/4, pk->h) ;
  Array newQPeaks3 = pkPeaksClip (newQPeaks4, 3, minCover/2, pk->h) ;

  Array newInputPeaks4 = pkPeaksCoverage (newInputPeaks, TRUE, 3, minInputCover+1, pk->h) ;   /* kept by the reference run */

  Array newInputPeaks1 = pkPeaksClip (newInputPeaks4, (3.0*minInputCover)/4, minInputCover + 1, pk->h) ;
  Array newInputPeaks2 = pkPeaksClip (newInputPeaks4, minInputCover/2, (3.0*minInputCover)/4, pk->h) ;
  Array newInputPeaks3 = pkPeaksClip (newInputPeaks4, 3, minInputCover/2, pk->h) ;
      
  
  /* the reference is a new run appended after the real runs */
  RPK *rc = arrayp (pk->runs, runMax, RPK) ;
  rc->name = hprintf (pk->h, "Reference_frequent.%d_to_%d", 3*minCover/4, minCover) ;
  rc->type = 2 ;
  rc->peaks = newQPeaks1 ;
  rc->run = 0 ;
  rc->h = pk->h ;

  rc = arrayp (pk->runs, runMax + 1, RPK) ;
  rc->name = hprintf (pk->h, "Reference_common.%d_to_%d", minCover/2, 3*minCover/4 - 1) ;
  rc->type = 2 ;
  rc->peaks = newQPeaks2 ;
  rc->run = 0 ;
  rc->h = pk->h ;

  rc = arrayp (pk->runs, runMax + 2, RPK) ;
  rc->name = hprintf (pk->h, "Reference_specific.%d_to_%d", 3, minCover/2 - 1) ;
  rc->type = 2 ;
  rc->peaks = newQPeaks3 ;
  rc->run = 0 ;
  rc->h = pk->h ;

  rc = arrayp (pk->runs, runMax + 3, RPK) ;
  rc->name = hprintf (pk->h, "Reference_all.%d_to_%d", 3, minCover) ;
  rc->type = 2 ;
  rc->peaks = newQPeaks4 ;
  rc->run = 0 ;
  rc->h = pk->h ;

  rc = arrayp (pk->runs, runMax + 4, RPK) ;
  rc->name = hprintf (pk->h, "Input_frequent.%d_to_%d", 3*minInputCover/4, minInputCover) ;
  rc->type = 2 ;
  rc->peaks = newInputPeaks1 ;
  rc->run = 0 ;
  rc->h = pk->h ;

  rc = arrayp (pk->runs, runMax + 5, RPK) ;
  rc->name = hprintf (pk->h, "Input_common.%d_to_%d", minInputCover/2, 3*minInputCover/4 - 1) ;
  rc->type = 2 ;
  rc->peaks = newInputPeaks2 ;
  rc->run = 0 ;
  rc->h = pk->h ;

  rc = arrayp (pk->runs, runMax + 6, RPK) ;
  rc->name = hprintf (pk->h, "Input_specific.%d_to_%d", 3, minInputCover/2 - 1) ;
  rc->type = 2 ;
  rc->peaks = newInputPeaks3 ;
  rc->run = 0 ;
  rc->h = pk->h ;

  rc = arrayp (pk->runs, runMax + 7, RPK) ;
  rc->name = hprintf (pk->h, "Input_all.%d_to_%d", 3, minInputCover) ;
  rc->type = 2 ;
  rc->peaks = newInputPeaks4 ;
  rc->run = 0 ;
  rc->h = pk->h ;

  ac_free (h) ;
} /* pkCreateRefPeaks */

/*************************************************************************************/
/* Count how many peaks in each run intersect the peaks of any other run
 */
static void pkCountIntersects (PK *pk)
{
  AC_HANDLE h = ac_new_handle () ;
  int runMax = arrayMax (pk->runs) ;
  Array areaCumul = arrayHandleCreate (runMax, long, h) ;
  Array widthCumul = arrayHandleCreate (runMax, long, h) ;
  Array r2nPeaks = arrayHandleCreate (runMax, int, h) ;
  Array r2w = arrayHandleCreate (runMax, int, h) ;
  Array r2wN = arrayHandleCreate (runMax, int, h) ;
  Array r2wHistos = arrayHandleCreate (runMax, KEYSET, h) ;
  Array r2a = arrayHandleCreate (runMax, int, h) ;
  Array r2aN = arrayHandleCreate (runMax, int, h) ;
  Array r2aHistos = arrayHandleCreate (runMax, KEYSET, h) ;

  Array xx = arrayHandleCreate (runMax * runMax + 1, int, h) ;
  for (int ir = 0 ; ir < runMax ; ir++)
    {
      RPK *rc1 = arrayp (pk->runs, ir, RPK) ;
      Array p1 = rc1->peaks ;
      if (p1 && arrayMax (p1))
	{
	  int ii, iMax =  arrayMax (p1) ;
	  PEAK *up = arrp (p1, 0, PEAK) ;
	  for (ii = 0 ; ii < iMax ; ii++)
	    {
	      array (areaCumul, ir, long) += up->area ;
	      array (widthCumul, ir, long) += up->width ;
	    }

	  for (int jr = 0 ; jr < runMax ; jr++)
	    {
	      RPK *rc2 = arrayp (pk->runs, jr, RPK) ;
	      Array p2 = rc2->peaks ;
	      if (p2 && arrayMax (p2))
		array (xx, runMax * ir + jr, int) = pkCountOneIntersect (p1, p2, rc2->type) ;
	    }
	}
    }
  
  ACEOUT ao = aceOutCreate (pk->outFileName, ".peakIntersects.txt", 0, h) ;
  aceOutDate (ao, "###", "PeakIntersects") ;
  aceOutf (ao, "## Number of peaks at threshold %d seen in line j touching a peak of column i\n", pk->min_peak_cover) ;
  aceOutf (ao, "## The table is not exactly symmetric because a single peak in line i may touch 2 peaks in column i or vice versa\n") ;
   aceOutf (ao, "\tnPeaks") ;
  for (int jr = 0 ; jr < runMax ; jr++)
    {
      RPK *rc2 = arrayp (pk->runs, jr, RPK) ;
      aceOutf (ao, "\t%s", rc2->name) ;
    }
  aceOut (ao, "\nNumber of peaks\t") ;
  for (int jr = 0 ; jr < runMax ; jr++)
    {
      RPK *rc2 = arrayp (pk->runs, jr, RPK) ;
      aceOutf (ao, "\t%d", arrayMax (rc2->peaks)) ;
    }
  for (int ir = 0 ; ir < runMax ; ir++)
    {
      RPK *rc1 = arrayp (pk->runs, ir, RPK) ;
      aceOutf (ao, "\n%s\t%d", rc1->name, arrayMax (rc1->peaks)) ;
      for (int jr = 0 ; jr < runMax ; jr++)
	aceOutf (ao, "\t%d", array (xx, runMax*ir+jr, int)) ;
    }
  aceOut (ao, "\n\n") ;
  if (debug1)
    for (int ir = 0 ; ir < runMax ; ir++)
      {
	RPK *rc1 = arrayp (pk->runs, ir, RPK) ;
	Array p1 = rc1->peaks ;
	for (int ii = 0 ; ii < arrayMax (p1) ; ii++)
	  {
	    PEAK *up = arrp (p1, ii, PEAK) ;
	    aceOutf(ao, "%s\t%d\t%d\t%d\tw=%d\th=%.0f\n"
		    , rc1->name, ii
		    , up->a1, up->a2, up->width, up->height
		    ) ;
	  }
      }
  
  ao = aceOutCreate (pk->outFileName, ".contacts.txt", 0, h) ;
  aceOutDate (ao, "###", "Contacts") ;
  aceOutf (ao, "## Number of peaks in a run contacting the peaks of at least n other runs at threshold %d\n", pk->min_peak_cover) ;
  aceOutf (ao, "## The first column, n=0, is the total number of peaks in that run\nRun", pk->min_peak_cover) ;
  aceOutf (ao, "\tNumber of peaks\tExtent on genome in kb\tNumber of Megabases in aligned in these peaks\tAligned Megabases\t%% of bases aligned in the peaks") ;
  for (int jr = 1 ; jr < runMax ; jr++)
    aceOutf (ao, "\t%d", jr) ;
  for (int ir = 0 ; ir < runMax ; ir++)
    {
      RPK *rc1 = arrayp (pk->runs, ir, RPK) ;
      aceOutf (ao, "\n%s", rc1->name) ;
      KEYSET ks = keySetCreate () ;
      Array p1 = rc1->peaks ;
      int ipMax = p1 ? arrayMax (p1) : 0 ;
      PEAK *up = ipMax ? arrp (p1, 0, PEAK) : 0 ;
      /* direct counts */
      for (int i = 0 ; i < ipMax ; i++, up++)
	keySet (ks, up->contact)++ ;
      /* cumuls */
      for (int i = runMax - 1, cumul = 0 ; i >= 0 ; i--)
	{
	  cumul += keySet (ks, i) ;
	  keySet (ks, i) = cumul ;
	}
      aceOutf (ao, "\t%d", keySet (ks, 0)) ;
      aceOutf (ao, "\t%.0f\t%.3f\t%.3f\t%.3f"
	       , array (widthCumul, ir, long)/1000.0
	       , array (areaCumul , ir, long)/1000.0
	       , rc1->alignedBases/1000.0
	       , 100*array (areaCumul , ir, long)/rc1->alignedBases
	       ) ;
      for (int jr = 1 ; jr < runMax ; jr++)
	aceOutf (ao, "\t%d", keySet (ks, jr)) ;

      ac_free (ks) ;
    }
  aceOut (ao, "\n\n") ;




  int wMAX = 4000 ;
  int wStep = 50 ;
  int aMAX = 2000 ;
  int aStep = 10 ;
  for (int ir = 0 ; ir < runMax ; ir++)
    {
      KEYSET ksw = keySetHandleCreate (h) ;
      KEYSET ksa = keySetHandleCreate (h) ;
      RPK *rc1 = arrayp (pk->runs, ir, RPK) ;
      Array p1 = rc1->peaks ;
      int n, wMax = 0, wNMax = 0, aMax = 0, aNMax = 0, pMax = arrayMax (p1) ;
      PEAK *up = pMax ? arrp (p1, 0, PEAK) : 0 ;

      for (int j = 0 ; j < pMax ; j++, up++)
	{
	  int w = up->width ;
	  w = w > wMAX ? wMAX : w ;
	  w = (w + wStep - 1)/wStep ;
	  keySet (ksw, w)++ ;
	  n = keySet (ksw, w) ;
	  if (w < wMAX && n > wMax) { wMax = n ; wNMax = w ;}

	  int a = up->area ;
	  a = a > aMAX ? aMAX : a ;
	  a = (a + aStep - 1)/aStep ;
	  keySet (ksa, a)++ ;
	  n = keySet (ksa, a) ;
	  if (a < aMAX && n > aMax) { aMax = n ; aNMax = a ;}
	}
      
      array (r2w, ir, int) = wNMax ;
      array (r2wN, ir, int) = wMax ;
      array (r2a, ir, int) = aNMax ;
      array (r2aN, ir, int) = aMax ;
      array (r2nPeaks, ir, int) = pMax ;
      array (r2wHistos, ir, KEYSET) = ksw ;
      array (r2aHistos, ir, KEYSET) = ksa ;
    }
  
  ao = aceOutCreate (pk->outFileName, ".peak_width.histo.txt", 0, h) ;
  aceOutDate (ao, "###", "Peak_width.histo") ;
  aceOutf (ao, "## Histogram of peaks widths at threshold %d\n", pk->min_peak_cover) ;
  for (int jr = 0 ; jr < runMax ; jr++)
    {
      RPK *rc2 = arrayp (pk->runs, jr, RPK) ;
      aceOutf (ao, "\t%s", rc2->name) ;
    }
  for (int w = 0 ; w * wStep <= wMAX ; w++)
    {
      aceOutf (ao, "\n%d", w * wStep) ;
      for (int ir = 0 ; ir < runMax ; ir++)
	{
	  KEYSET ksw = array (r2wHistos, ir, KEYSET) ;
	  aceOutf (ao, "\t%d", keySet (ksw, w)) ;
	}
    }
  aceOut (ao, "\n") ;
  
  ao = aceOutCreate (pk->outFileName, ".peak_width.max.txt", 0, h) ;
  aceOutDate (ao, "###", "Peak_width.histo") ;
  aceOutf (ao, "## Per run, what is the most frequent width of the peaks\n") ;
  aceOutf (ao, "## Width is the most frequently observed peak-width\n") ;
  aceOutf (ao, "## N is the number of peaks with this width\n")  ;
  aceOutf (ao, "## Total is the total number of peaks seen in this run\n") ;
  aceOutf (ao, "#Run\tMost frequent width\tNumber of peaks thr_%d with this width\tTotal number of peaks thr_%d\tmedian\n", pk->min_peak_cover, pk->min_peak_cover) ;
  for (int ir = 0 ; ir < runMax ; ir++)
    {
      RPK *rc1 = arrayp (pk->runs, ir, RPK) ;
      aceOutf (ao, "\n%s", rc1->name) ;
      aceOutf (ao, "\t%d\t%d\t%d\t%d"
	       , wStep * array (r2w, ir, int)
	       , array (r2wN, ir, int)
	       , array (r2nPeaks, ir, int)
	       , rc1->medianNoZero
	       /*	       , rc1->minCover * rc1->medianNoZero */
	       ) ;
    }
  aceOut (ao, "\n") ;
  
  ao = aceOutCreate (pk->outFileName, ".peak_area.histo.txt", 0, h) ;
  aceOutDate (ao, "###", "Peak_area.max") ;
  aceOutf (ao, "## Histogram of peaks area (kb aligned in the peak) at threshold %d\n", pk->min_peak_cover) ;
  for (int jr = 0 ; jr < runMax ; jr++)
    {
      RPK *rc2 = arrayp (pk->runs, jr, RPK) ;
      aceOutf (ao, "\t%s", rc2->name) ;
    }
  for (int a = 0 ; a * aStep <= aMAX ; a++)
    {
      aceOutf (ao, "\n%d", a * aStep) ;
      for (int ir = 0 ; ir < runMax ; ir++)
	{
	  KEYSET ksa = array (r2aHistos, ir, KEYSET) ;
	  aceOutf (ao, "\t%d", keySet (ksa, a)) ;
	}
    }
  aceOut (ao, "\n") ;
  
  ao = aceOutCreate (pk->outFileName, ".peak_area.max.txt", 0, h) ;
  aceOutDate (ao, "###", "Peak_area.max") ;
  aceOutf (ao, "## Per run, what is the most frequent width of the peaks\n") ;
  aceOutf (ao, "## Area in kb is the most frequently observed peak-area\n") ;
  aceOutf (ao, "## N is the number of peaks with this area\n")  ;
  aceOutf (ao, "## Total is the total number of peaks seen in this run\n") ;
  aceOutf (ao, "#Run\tMost frequent peak area thr_%d in kb\tNumber of peaks at this maximal area\tTotal number of peaks thr_%d\n", pk->min_peak_cover, pk->min_peak_cover) ;
  for (int ir = 0 ; ir < runMax ; ir++)
    {
      RPK *rc1 = arrayp (pk->runs, ir, RPK) ;
      aceOutf (ao, "\n%s", rc1->name) ;
      aceOutf (ao, "\t%d\t%d\t%d"
	       , aStep * array (r2a, ir, int)
	       , array (r2aN, ir, int)
	       , array (r2nPeaks, ir, int)
	       ) ;
    }
  aceOut (ao, "\n") ;

  ac_free (h) ;
} /* pkCountIntersects */

/*************************************************************************************/


/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: peaks --db MetaDB -p project [options] [actions]\n"
	   "//      try: -h --help\n"
	   "// Objective: \n"
	   "//     Analyze the peaks of the coverage plots of a ChIP-seq experiment\n"
	   "// Required:\n"
	   "//    --db <dir> : acedb MetaDB directory\n"
	   "//    -p, --project <project> : project name\n"
	   "//    environment variable chromSetAll : blank separated chromosome names, i.e. source LIMITS\n"
	   "// General Options:\n"
	   "//    --group <pattern> : restrict to the runs or groups matching this pattern\n"
	   "//    --min_peak_cover <n> : minimal coverage of a peak\n"
	   "//    -o <outFilePrefix> : output file prefix, all exported files will start with that prefix\n"
	   "//    --gzo : the output files are gziped\n"
	   "//    -i, --min_support, -orderBy : accepted, see the code\n"
	   "// Actions, executed in this order:\n"
	   "//    --getStrandShifts\n"
	   "//    --makeGroupWiggles\n"
	   "//    --makePeaks\n"
	   "//    --makeShiftedNonStrandedPeaks\n"
	   "//    --intersects : build the reference peaks and count the intersects between the runs\n"
	   "//\n"
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
  char *cp ;
  PK pk ;
  AC_HANDLE h = 0 ;
  char commandBuf [4000] ;

  messErrorInit (argv[0]) ;

  h = ac_new_handle () ;
  memset (&pk, 0, sizeof (PK)) ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
  pk.h = h ;

  if (argc < 2)
    usage (commandBuf, argc, argv) ;
  if (getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)||
      getCmdLineOption (&argc, argv, "-h", 0)
      )
    usage (commandBuf, 1, argv) ;

  { 
    int ix ;
    for (ix = 0, cp = commandBuf ;  ix < argc && cp + strlen (argv[ix]) + 1 < commandBuf + 3900 ; cp += strlen (cp), ix++)
      sprintf(cp, "%s ", argv[ix]) ;
  }

  /* consume optional args */
  getCmdLineOption (&argc, argv, "-i", &(pk.inFileName)) ;
  getCmdLineOption (&argc, argv, "-o", &(pk.outFileName)) ;
  pk.gzo = getCmdLineBool (&argc, argv, "--gzo") ;

  getCmdLineInt (&argc, argv, "--min_support", &(pk.min_support)) ;
  getCmdLineInt (&argc, argv, "--min_peak_cover", &(pk.min_peak_cover)) ;
  getCmdLineOption (&argc, argv, "--project", &pk.project) ;
  getCmdLineOption (&argc, argv, "--db", &pk.dbName) ;
  getCmdLineOption (&argc, argv, "-p", &pk.project) ;
  getCmdLineOption (&argc, argv, "--group", &pk.group) ;
  pk.intersects = getCmdLineBool (&argc, argv, "--intersects") ;
  pk.makeGroupWiggles = getCmdLineBool (&argc, argv, "--makeGroupWiggles") ;
  pk.getStrandShifts = getCmdLineBool (&argc, argv, "--getStrandShifts") ;
  pk.makePeaks = getCmdLineBool (&argc, argv, "--makePeaks") ;
  pk.makeShiftedNonStrandedPeaks = getCmdLineBool (&argc, argv, "--makeShiftedNonStrandedPeaks") ;
  
  pk.orderBy = "Sorting_title" ;
  getCmdLineOption (&argc, argv, "-orderBy", &pk.orderBy) ;

  if (argc > 1) 
    usage ("Unknown argument", argc, argv) ;
  fprintf (stderr, "// %s start : \n", timeShowNow()) ;

  if (! pk.dbName)
    {
      fprintf (stderr, "Sorry, missing parameter --db, please try peaks -help\n") ;
      exit (1) ;
    }
  {
    const char *errors ;
    
    pk.db = ac_open_db (pk.dbName, &errors);
    if (! pk.db)
      {
	fprintf (stderr, "Failed to open db %s, error %s", pk.dbName, errors) ;
	exit (1) ;
      }
  }
  getCmdLineOption (&argc, argv, "-project", &pk.project) ;

  if (! pk.project)
    {
      fprintf (stderr, "Missing argument -p <project>\n") ;
      goto done ;
    }
    
  pkGetChroms (&pk) ;
  if (! dictMax (pk.chromDict))
    {
      fprintf (stderr, "Sorry the env variable $chromSetAll is not set, please try source LIMITS\n") ;
      goto done ;
    }

  if (pk.getStrandShifts)
    {
      pkGetStrandShifts (&pk) ;
    }

  if (pk.makeGroupWiggles) 
    {
      pkMakeGroupWiggles (&pk) ;
    }
  
  if (pk.makePeaks)
    {
      pkMakePeaks (&pk) ;
    }

  if (pk.makeShiftedNonStrandedPeaks)
    {
      pkMakeShiftedNonStrandedPeaks (&pk) ;
    }

  if (pk.intersects)
    {
      if (! pkGetRuns (&pk)) 
	{
	  fprintf (stderr, "Sorry, I found no run in project %s\n", pk.project) ; 
	  goto done ;
	}
      pkGetPeaks (&pk) ;
      pkCreateRefPeaks (&pk) ;
      pkCountIntersects (&pk) ;
    }
    
 done:
  ac_db_close (pk.db) ;
  
  {
    int mx ;
    messAllocMaxStatus (&mx) ; 
    fprintf (stderr, "// %s done: max memory %d Mb\n", timeShowNow(), mx) ;
  }

  return 0 ;
} /* main */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
