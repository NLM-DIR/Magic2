/*  File: sv.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2015
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
 * Ananlysis of the SV: strutural variations
 * using sequence overhangs and wiggles
 * from the genome data
 */

/* #define ARRAY_CHECK */

#include "ac.h"
#include "bitset.h"
#include "channel.h"
#include "remote_channel.h"

typedef struct svStruct {
  AC_HANDLE h ;
  Array runs ;
  DICT *classeDict ;
  DICT *runDict ;
  const char *runsFileName ;
  const char *inFileName ;
  const char *outFileName ;
  Stack info ;
  int circle, deletion, minideletion, microdeletion, duplication, miniduplication, microduplication, palindrome, inversion, translocation ;

  BOOL de_uno, de_duo, de_paribus, de_av ;
  int nRuns ;
  int min_support ;
  int seedLength ;
  int overhangLength ;

  BOOL de_wiggle ;

  CHAN *todo, *done ;
  int max_threads ;
  int max_tasks ;
  BOOL gzo ;
  
  BOOL serverSide, submit ;
  int clientId ;
  TASK *serverTask ;
  TASK *clientTask ;
} SV ;

typedef struct runStruct {
  int type ;
  int run ;
  int sex ; /* 0 -> autosome, 1 : Y, 2 : X */
  float cumul, cumulX, cumulY ;
  Array cover, geneCover ;
  int machine, sample, tissue, title, runId, sortingTitle, sortingTitle2, otherTitle ;
  int histo[101] ;
} RC ;

/* Y branch point */
#define BBWIDTH 15
#define BBSTEP 10
typedef struct yStruct {
  int classe ;
  int  chromA, chromB, a1, a2, b1, b2, type, sPlus, sMinus, width, downA, downB ;
  int v ; /* 1 when commmon part is left and V is right, -1 if V is on the left and common part on the 3' side */
  int bb0, bb1, bb2, bb3 ; /* bb0: common part, bb1 bb2 the 2 branches:  offsets in BB stem-array */
} YY ;

/* one of the 3 stems of the Y */
typedef struct stemStruct {
  int type ;  /* 0,1,2 for the stem is the b0, b1 or b2 of the Y */
  int run ;
  int chrom, a1 , a2 ; /* always oriented like the genome */
  float ok[BBWIDTH] ;  /* ok wiggle, oriented a1 towards a2 */
  float nu[BBWIDTH] ;  /* ok wiggle, oriented a1 towards a2 */
  float pp[BBWIDTH] ;  /* partial wiggle */
  int y[5] ;    /* list of y */
} BB ;

typedef BB *BBP ;

typedef struct svuStruct { int status, clientId ; char runName[256], chromName[256] ; } SDU ;
typedef struct swStruct {
  AC_HANDLE h ; 
  Array yys, bbs, bbps ;
  DICT *dict, *typeDict, *chromDict, *genePairDict ;
  KEYSET genePairKs ;
  int circleFeet ; 
} SW ;

/*************************************************************************************/
/*************************************************************************************/
/* aorder an arry of pointers into a BB array */
static int bbpOrder (const void *a, const void *b)
{
  const BB *up = *(const BB **)a, *vp = *(const BB **)b ;
  int z ;

  z = up->chrom - vp->chrom ; if (z) return z ;
  z = up->a1 - vp->a1 ; if (z) return z ;

  return 0 ;
} /* bbpOrder*/

/*************/

#ifdef JUNK

typedef struct breakPointStruct {
  int bk ;
  int chrom, a1, a2, ln ;
  BOOL isDown ;
  int type ;
  float cumul[2], median[2], u2u, u2nu ;
} BK ;

/*************/

static int bkOrder (const void *a, const void *b)
{
  const BK *up = (const BK *)a, *vp = (const BK *)b ;
  int z ;

  z = up->chrom - vp->chrom ; if (z) return z ;
  z = up->a1 - vp->a1 ; if (z) return z ;
  z = up->type - vp->type ; if (z) return z ;

  return 0 ;
} /* bkOrder*/

/*************************************************************************************/

  const char *project ;
  const char *title ;
  const char *inFileName ;
  const char *outFileName ; 
  const char *gtitleFileName ;
  const char *dbName = 0 ;
  Array bks ;
  AC_DB db ;


	   "//    -i filname : data file to be analyzed, default stdin"
	   "//    --gzi : please gunzip the input file, default action is the file is called *.gz\n"
	   "//      --project : MAGIC Project name\n"
	   "//      --db MetaDB_directory : Directory of the acedb meta-database of the project\n"  
	   "//      --gtitle gTitleFileName \n"

 
  getCmdLineOption (&argc, argv, "-i", &(sv.inFileName)) ;
  getCmdLineOption (&argc, argv, "-o", &(sv.outFileName)) ;

  sv.gzi = getCmdLineOption (&argc, argv, "--gzi", 0) ;

  getCmdLineOption (&argc, argv, "--db", &dbName) ;
  getCmdLineOption (&argc, argv, "--gtitle", &(sv.gtitleFileName)) ;
  getCmdLineOption (&argc, argv, "--title", &(sv.title)) ;
  getCmdLineOption (&argc, argv, "--project", &(sv.project)) ;

  sv.bks = arrayHandleCreate (100000, BK, h) ;

  if (sv.gtitleFileName) 
    svParseTitles (&sv) ;

  if (dbName)
    {
      char *errors = 0 ;
      
      sv.db = ac_open_db (dbName, &errors);
      if (! sv.db)
	messcrash ("Failed to open db %s, error %s", dbName, errors) ;
    }
  if (sv.db) 
    ac_db_close (sv.db) ;

/* only runs in these lists will be considered */
static void svParseTitles (SV *sv)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (sv->gtitleFileName, FALSE, h) ;
  int run ;
  DICT *dict = sv->runDict ;
  RC *rc ;
  Stack s = sv->info ;

  if (ai)
    {
      aceInSpecial (ai, "\n") ;
      while (aceInCard (ai))
	{
	  char cutter ;
	  const char *ccp = aceInWordCut (ai, "\t", &cutter) ;
	  if (! ccp) continue ;
	  ccp = ac_unprotect (ccp, h) ;
	  if ( ! dictFind (dict, ccp, &run))
	    continue ;
	  rc = arrayp (sv->runs, run, RC) ;

	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->machine = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->sample = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  /* rc->system = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ; */
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->tissue = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->title = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  /* rc->system2 = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ; */
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->runId = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  /* rc->machine2 = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ; */
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->sortingTitle = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->sortingTitle = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->sortingTitle2 = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->otherTitle = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	}
    }
  ac_free (h) ;
  return ;
} /* svParseTitles */
#endif

 /*************************************************************************************/
/* only runs in these lists will be considered */
static int svParseRunList (SV *sv)
{
  int nn = 0 ;
  const char *ff = 0 ;

  ff = sv->runsFileName ; 
  if (ff)
    {
      AC_HANDLE h1 = ac_new_handle () ;
      ACEIN ai= aceInCreate (ff, FALSE, h1) ;
      int run = 0 ;
      RC *rc ;
      const char *ccp ;
      DICT *dict = sv->runDict ;

      if (ai)
	while (aceInCard (ai))
	  if ((ccp = aceInWord (ai)))
	    {
	      dictAdd (dict, ccp, &run) ;
	      rc = arrayp (sv->runs, run, RC) ;
	      rc->run = run ;
	      nn++ ;
	    }
      ac_free (h1) ;
    } 
  fprintf (stderr, "// parsed %d runs\n", nn) ;
  sv->nRuns = nn ;
  return nn ;
} /* svParseRunList */

/*************************************************************************************/
/*************************************************************************************/

static int svYyOrder (const void *a, const void *b) 
{
  const YY *up = (const YY *) a ; 
  const YY *vp = (const YY *) b ; 
  int n ;
 
  n = up->classe - vp->classe ; if (n) return n ;
  n = up->chromA - vp->chromA ; if (n) return n ;
  n = up->downA - vp->downA ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;

  n = up->chromB - vp->chromB ; if (n) return n ;
  n = up->downB - vp-> downB ; if (n) return n ;
  n = up->b1 - vp->b1 ; if (n) return n ;

  n = up->type - vp->type ; if (n) return n ;

  return 0 ;
} /* svYyOrder */

/*************************************************************************************/

static int svYyCompress (Array yys, int min)
{ 
  int i, j, n = arrayMax (yys) ;
  YY *up, *vp ;
  
  if (n > 1) 
    {
      arraySort (yys, svYyOrder) ;
      
      for (i = 1, j = 0, up = arrp (yys, 1, YY), vp = up - 1 ; i < n ; i++, up++)
	{
	  if (svYyOrder (up, vp))
	    {
	      vp = arrp (yys, ++j, YY) ;
	      if (i > j) *vp = *up ;
	    } 
	  else /* same yy, increment the supports */
	    {
	      vp->sPlus  += up->sPlus ;
	      vp->sMinus += up->sMinus ;
	    }
	}
      n = arrayMax (yys) = j + 1 ;
      for (i = 0, j = 0, up = vp = arrp (yys, 0, YY) ; i < n ; i++, up++)
	{
	  if (up->sPlus + up->sMinus >= min)
	    {
	      if (i > j) *vp = *up ;
	      j++ ; vp++ ;
	    }
	}
      n = arrayMax (yys) = j ;
      
    }
  return n ;
} /* svYyCompress */

/*************************************************************************************/
#ifdef JUNK
/* exchangethe 2 segements, but do not change strand */
static void svSwitch (DICT *dict, YY *dd)
{
  int x ;

  x = dd->sMinus ; dd->sMinus = dd->sPlus ; dd->sPlus = x ;
  x = dd->chromA ; dd->chromA = dd->chromB ; dd->chromB = x ;
  x = dd->a1 ; dd->a1 = dd->b2 ; dd->b2 = x ;
  x = dd->a2 ; dd->a2 = dd->b1 ; dd->b1 = x ;
  x = dd->downA ; dd->downA = dd->downB ; dd->downB = x ; 
  if (0 && dd->type)
    {
      char tt[6], uu[6] ;
      
      strncpy (uu, dictName(dict, dd->type), 6) ;
      tt[0] = uu[3] ; tt[1] = uu[4] ; tt[2] = uu[2] ; tt[3] = uu[0] ; tt[4] = uu [1] ; tt[5] = 0 ;
      dictAdd (dict, tt, &(dd->type)) ;
    }
} /* svSwitch  */
#endif
/*************************************************************************************/

/* complement en place the name of an intron foot ct_ac -> gt->ag */
static void svSwitchComplement (DICT *dict, YY *dd)
{   
  int type = dd->type ;
  int x ;

  x = dd->sMinus ; dd->sMinus = dd->sPlus ; dd->sPlus = x ;
  x = dd->chromA ; dd->chromA = dd->chromB ; dd->chromB = x ;
  x = dd->a1 ; dd->a1 = dd->b2 ; dd->b2 = x ;
  x = dd->a2 ; dd->a2 = dd->b1 ; dd->b1 = x ;
  x = dd->downA ; dd->downA = dd->downB ; dd->downB = x ; 

  if (1 && type)
    {
      int i ;
      char tt[6], uu[6] ;
      
      strncpy (uu, dictName(dict, type), 6) ;
      tt[5]  = 0 ;
      for (i = 0 ; i < 5 ; i++)
	tt[i] = (uu[4-i] == '_' ? '_' : dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)uu[4-i]])]) ;
       dictAdd (dict, tt, &(dd->type)) ;
    }
  dd->downA *= -1 ;
  dd->downB *= -1 ;

  return ;
} /* svSwitchComplement */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

typedef struct ixpStruct { BOOL reverse, isDown, isDownA ; int type, chrom, chromA, dan, ddn, bestdx, dx1, exon1, exon2, u1, u2, du ; char foot[8] ; } IXP ;

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

/*************************************************************************************/

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
/* bp equivalent of the oligo complexity
 *
 * input is dna, x1/x2 in bio-coords
 */
int hookEntropy (const char *hook, Array ee, int *oldNnp)
{
  int ss = 0 ;
  int na = 0, nt = 0, ng = 0, nc = 0, nn ;
  const char *ccp ;
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
  if (nn > 1 && nn != *oldNnp)
    {
      double s = 0, log4 = log(4.0) ; /* divide by log4 to be in base equivalent */ 
      int j ;

      ee = arrayReCreate (ee, nn, int) ;
      for (j = 1 ; j <= nn ; j++)
        { s = j ; s /= nn ; array (ee, j, int) = - (int) (1000.0 * j * log(s)/log4 ) ; }
      *oldNnp = nn ;
      array (ee, 0, int) = 0 ;
    }

  ss = arr (ee, na, int) + arr (ee, nt, int) + arr (ee, ng, int) + arr (ee, nc, int) ;
  ss = (ss + 499)/1000 ;
  return ss ;
} /* hookEntropy */

/*************************************************************************************/
/* collate the number of tags supporting the same jump */
static int svNewIntronsCompress (Array donors, int limit)
{
  int ii, jj = 0 ;
  DA *up, *vp ;
  
  if (! donors || ! arrayMax (donors))
    return 0 ;
  
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
}  /* svNewIntronsCompress */

/*************************************************************************************/

static int svYy2Bb (const SV *sv, SDU *sdu, SW* sw)
{
  Array bbs = 0 ;
  int i, j, nn = arrayMax (sw->yys) ;
  YY *yy ;
  BB *bb0, *bb1, *bb2, *bb3 ;

  sw->bbs = bbs = arrayHandleCreate (3*nn, BB, sw->h) ;

  for (i = j = 0, yy = arrp (sw->yys, i, YY) ; i < nn ; i++, yy++)
    {
      /* start with bb3 because of possible bbs reallocation */
      yy->bb3 = j + 3 ; bb3 = arrayp (bbs, j+3, BB) ;
      yy->bb0 = j ; bb0 = arrayp (bbs, j, BB) ;
      yy->bb1 = j + 1 ; bb1 = arrayp (bbs, j+1, BB) ;
      yy->bb2 = j + 2 ; bb2 = arrayp (bbs, j+2, BB) ;
      j += 4 ;

      if (yy->a2 < yy->b1)
	{
	  bb0->chrom = yy->chromA ; 
	  bb1->chrom = yy->chromA ;
	  bb0->a2 = yy->a2 - BBSTEP ; bb0->a1 = yy->a2 - BBWIDTH * BBSTEP ; 
	  bb1->a1 = yy->a2 + BBSTEP ; bb1->a2 = yy->a2 + BBWIDTH * BBSTEP ; 
	  
	  bb2->chrom = yy->chromB ; 
	  bb3->chrom = yy->chromB ;
	  bb2->a2 = yy->b1 - BBSTEP ; bb2->a1 = yy->b1 - BBWIDTH * BBSTEP ; 
	  bb3->a1 = yy->b1 + BBSTEP ; bb3->a2 = yy->b1 + BBWIDTH * BBSTEP ; 
	}
      else
	{
	  bb2->chrom = yy->chromA ; 
	  bb3->chrom = yy->chromA ;
	  bb2->a2 = yy->a2 - BBSTEP ; bb2->a1 = yy->a2 - BBWIDTH * BBSTEP ; 
	  bb3->a1 = yy->a2 + BBSTEP ; bb3->a2 = yy->a2 + BBWIDTH * BBSTEP ; 
	  
	  bb0->chrom = yy->chromB ; 
	  bb1->chrom = yy->chromB ;
	  bb0->a2 = yy->b1 - BBSTEP ; bb0->a1 = yy->b1 - BBWIDTH * BBSTEP ; 
	  bb1->a1 = yy->b1 + BBSTEP ; bb1->a2 = yy->b1 + BBWIDTH * BBSTEP ; 
	}

      if (bb0->a1 < 0) bb0->a1 = 0 ;
      if (bb1->a1 < 0) bb1->a1 = 0 ;
      if (bb2->a1 < 0) bb2->a1 = 0 ;
      if (bb3->a1 < 0) bb3->a1 = 0 ;

      if (bb0->a2 < 0) bb0->a2 = 0 ;
      if (bb1->a2 < 0) bb1->a2 = 0 ;
      if (bb2->a2 < 0) bb2->a2 = 0 ;
      if (bb3->a2 < 0) bb3->a2 = 0 ;
    }
  nn = arrayMax (bbs) ;
  sw->bbps = arrayHandleCreate (nn, BBP, sw->h) ;
  if (nn) array (sw->bbps, nn - 1, BBP) = 0 ;
  for (i = 0 ; i < nn ; i++)
    arr (sw->bbps, i, BBP) = arrp (sw->bbs, i, BB) ;
  arraySort (sw->bbps, bbpOrder) ;
  return nn ;
} /* svYy2Bb  */

/*************************************************************************************/
 /* step found, alter the coordinates of all BB on that chrom */
static int svWiggleRegularizeCoodinates (SW* sw, int chrom, int start, int step, float  *scalep)
{
 BB **bbp = 0 ; 
 int ib = 0 ;
 int ibMax = arrayMax (sw->bbs) ;
 int ddx, dx = start % step ;
 /* int dd = BBSTEP % step ; */

 for (ib = 0,  bbp = arrp (sw->bbps, ib, BBP) ; ib < ibMax ; ib++, bbp++)
   {
     BB *bb = *bbp ;
     if (bb->chrom > chrom)
       break ;
     if (bb->chrom < chrom)
       continue ;
     ddx = (bb->a1 - dx) % step ;
     bb->a1 -= ddx ; bb->a2 -= ddx ; /* now we are in phase with the wiggle */
     if (bb->a1 < 0) bb->a1 = 0 ;  
     if (bb->a2 < 0) bb->a2 = 0 ;
   }
 return 1 ;
} /* svWiggleRegularizeCoodinates */

/*************************************************************************************/

static int svWiggleDoOne (const SV *sv, SDU *sdu, SW* sw, int chrom, int ok)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = 0 ;
  int nn = 0, state = 0, line = 0, x, ib = 0, nk = 0 ;
  float y, scale = 1 ;
  int start = 0, step = 10 ;
  char *command  = 0 ;
  BB **bbp = 0 ;
  int ibMax = arrayMax (sw->bbs) ;
  
  bbp = arrp (sw->bbps, ib, BBP) ;
  command = hprintf (h, "gunzip -c %s/tmp/WIGGLERUN/%s/%s/R.chrom.frns.%s.BF.gz" 
		     , sv->inFileName ? sv->inFileName : "."
		     , sdu->runName
		     , dictName (sw->dict, chrom) 
		     , ok == 2 ? "u" : (ok == 1 ? "nu" : "pp")
		     ) ; 
  if (0) command = "cat tatou" ;
  ai = aceInCreateFromPipe (command, "r", 0, h) ;
  
  aceInSpecial (ai, "\n") ;
  while (ai && aceInCard (ai))
    {
      char *cp ;
      line++ ;

      switch (state)
	{
	case 0:	
	  cp = aceInWord (ai) ;    /* probe id, count the multiplicity */
	  if (! cp)  continue ;
	  if (line > 20)
	    goto done ;
	  if (strcmp (cp, "fixedStep"))
	    continue ;
	  cp = strstr (aceInPos (ai), "start=") ;
	  if (! cp)  continue ;
	  if (sscanf (cp, "start=%d step=%d", &start, &step) != 2)
	    continue ;
	  x = start - step ;
	  state = 1 ; scale = 1 ;
	  if (BBSTEP % step)
	    svWiggleRegularizeCoodinates (sw, chrom, start, step, &scale) ;
	  /* was scale = step / (BBSTEP * BBSTEP) ; */
	  continue ;
	  break ;
	case 1:
	  x += step ; y = 0 ;
	  aceInFloat (ai, &y) ;
	  if (y <= 0)
	    continue ;
	  break ;	  
	}
      nk++ ;
      if (x >= 257870)
	x = x ;
      if (ib == ibMax && ib > 0)
	{ bbp-- ; ib-- ; }
      for ( ; ib > 0 && (*bbp)->a2 > x - 10 ;  bbp--, ib--)
	{} 
      for (; ib < ibMax && (*bbp)->a2 < x - BBSTEP ; ib++, bbp++)
	{}
      if (ib == ibMax) break ;
      for (; ib < ibMax && (*bbp)->a1 < x + BBSTEP ; ib++, bbp++)
	{
	  int a, j, k ;
	  BB *bb = *bbp ;

	  for (k = 0, a = bb->a1 ; a < bb->a2 && k < BBWIDTH ; k++, a += BBSTEP)
	    {
	      j = x - a ; if (j < 0) j = -j ;
	      if (j == 0)
		{
		  nn += y ;
		  if (ok == 2)
		    bb->ok[k] += y * scale ;
		  else if (ok == 1)
		    bb->nu[k] += y * scale ;
		  else
		    bb->pp[k] += y * scale ;
		}
#ifdef JUNK
	      if (j < BBSTEP), but now we regularized the coordinates */
		{
		  nn += y ;
		  if (ok == 2)
		    bb->ok[k] += (BBSTEP - j) * y * scale ;
		  else if (ok == 1)
		    bb->nu[k] += (BBSTEP - j) * y * scale ;
		  else
		    bb->pp[k] += (BBSTEP - j) * y * scale ;
		}
#endif
	    }
	}
    }
  fprintf (stderr, "x = %d nk=%d\n", x, nk) ;
 done :
  ac_free (h) ;
  return nn ;
} /* svWiggleDoOne */

/*************************************************************************************/

static int svWiggleDo (const SV *sv, SDU *sdu, SW* sw) 
{
  int i, oldChrom = 0, nn = 0 ;
  int iMax = arrayMax (sw->bbps) ;
  BB **bbp ;

  for (i = 0, bbp = arrp (sw->bbps, i, BB*) ; i < iMax ; i++, bbp++)
    {
      BB *bb = *bbp ;
      if (bb->chrom != oldChrom)
	{
	  nn += svWiggleDoOne (sv, sdu, sw, bb->chrom,  2) ;
	  nn += svWiggleDoOne (sv, sdu, sw, bb->chrom, 1) ;
	  nn += svWiggleDoOne (sv, sdu, sw, bb->chrom,  0) ;
	  oldChrom = bb->chrom ;
	}
    }
  wego_log (hprintf (0, "svWiggleDo found %d support in %d BB in run %s", nn, iMax, sdu->runName)) ;
  return nn ;
} /* svWiggleDo */

/*************************************************************************************/

static int  svDeDuoParseOverhangs (const SV *sv, SDU *sdu, SW *sw)
{
  AC_HANDLE h = ac_new_handle () ;
  IXP *up, *vp ;
  Array exportIntrons = 0 ; 
  int nOk0 = 0, nOk1 = 0, nOk2 = 0, nOkRejected = 0, i, j, ia = 0, id = 0, ipA = 0, iSL = 0
    , ja, k1, k2, bestk, bestk1, bestk2, bestdx, bestdx1, bestdx2, a1, x1, x2,  dx1, dx2, score, chrom, exon, hook, maxA, maxD, maxSL ;
  const char *exon1, *exon2, *hook1, *hook2 ;
  char *cp, buf[1001] ;
  const char *ccp, *ccq ;
  int overhangLength = sv->overhangLength ;
  int seedLength = sv->seedLength ; 
  int min_support = sv->min_support ;
  int entropy, mult, multT ;
  int microLength = 36, miniLength = 1000 ;
  char polyABuf[overhangLength + 1] ;
  int polyA, polyT, SL ;
  BOOL isDonor, isAcceptor, isHit, isPa, isPt, isDown, reverse ;
  BitSet dBitSet = 0, aBitSet = 0, truedBitSet = 0, trueaBitSet = 0 ;
  ACEIN ai = 0 ;
  ACEOUT ao = 0 ;
  DA *dd = 0, *da = 0 ;
  Array polyAs = arrayHandleCreate (100000, DA, h) ;
  Array SLs = arrayHandleCreate (100000, DA, h) ;
  Array donors = arrayHandleCreate (100000, DA, h) ;
  Array acceptors = arrayHandleCreate (100000, DA, h) ;
  Array ee = arrayHandleCreate (32, int, h) ;
  DICT *dict = sw->dict ;
  char *command  = 0 ;
  char aaaBuf[101] ;
  char f1[6], f2[6], foot[8] ; 
  int oldNn = -1 ;

  memset (aaaBuf, 'A', 100) ;
  memset (polyABuf,'a', overhangLength) ;
  polyABuf [overhangLength] = 0 ;
  dictAdd (dict, polyABuf, &polyA) ;
  memset (polyABuf,'t', overhangLength) ;
  polyABuf [overhangLength] = 0 ;
  dictAdd (dict, polyABuf, &polyT) ;
  
  command = hprintf (h, "gunzip -c %s/tmp/PHITS_virus/%s/*.overhangs.gz"
		     , sv->inFileName ? sv->inFileName : "."
		     , sdu->runName) ; 
  if (0) command = "cat tatou" ;
  ai = aceInCreateFromPipe (command, "r", 0, h) ;
  command = hprintf (h, ".%s.de_duo.txt", sdu->runName) ;
  if (1)
    {
      ao = aceOutCreate (sv->outFileName, command, sv->gzo, h) ;
      aceOutDate (ao, "###", "de-duo breakpoint analysis") ;
      aceOutf (ao, "## Requested minimal_support %d\n",  sv->min_support) ;
    }
  if (sv->clientId) 
    wego_log (hprintf (h, "## client %d, processing %s\n", sv->clientId, aceInFileName (ai))) ;
  aceInSpecial (ai, "\n") ;
  
  /* Parse the output of clipalign -splice overhangs file
   * Input should look like:
   * CHROMOSOME_I    2918821      probe  21      taggcttaggcttaggcttaggc caagcctaagcccaa DONOR   Forward
   */

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
      entropy = hookEntropy (cp, ee, &oldNn) ;
  
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

  wego_log (hprintf (0, "// svNewIntrons parsed %u donors %d acceptors %d SL %d polyA,  rejected %d overhangs with low entropy"
		       , arrayMax (donors), arrayMax (acceptors), arrayMax (SLs), arrayMax (polyAs), nOkRejected)) ;

  /* sort */
  if (0) showDA ("Donors raw", dict, donors) ;
  arraySort (donors, DAorderHook) ;
  if (0) showDA ("Donors sorted", dict, donors) ;
   wego_log (hprintf (0, "// svNewIntrons sorted the donors %u donors %u acceptors"
			, arrayMax (donors), arrayMax (acceptors))) ;

  if (0) showDA ("Acceptors raw", dict, acceptors) ;
  arraySort (acceptors, DAorderHook) ;
  if (0) showDA ("Acceptors sorted", dict, acceptors) ;

  wego_log (hprintf (0, "// svNewIntrons sorted the acceptors %u donors %u acceptors"
		       , arrayMax (donors), arrayMax (acceptors))) ;
  
  arraySort (SLs, DAorderHook) ;
  arraySort (polyAs, DAorderHook) ;

  wego_log (hprintf (0, "// svNewIntrons sorted %d Sls, %d polyA"
		       , arrayMax (SLs), arrayMax (polyAs))) ;
  
  /* merge identical boundaries, and keep those with enough support */
  maxD = svNewIntronsCompress (donors, min_support) ;
  maxA = svNewIntronsCompress (acceptors, min_support) ;
  maxSL = svNewIntronsCompress (SLs, min_support) ;
  /* int maxpA = svNewIntronsCompress (polyAs, min_support) ; */
  if (!maxA || ! maxD)
    goto done ;
  arraySort (donors, DAorderN) ;
  arraySort (acceptors, DAorderN) ;
  arraySort (polyAs, DAorderN) ;
  arraySort (SLs, DAorderN) ;

  if (0) showDA ("Donors compressed", dict, donors) ;
  if (0) showDA ("Acceptors compressed", dict, acceptors) ;
  wego_log (hprintf (0, "// svNewIntrons consolidated %u donors %u acceptors %u SL %u pA"
		       , arrayMax (donors), arrayMax (acceptors)
		       , arrayMax (SLs), arrayMax (polyAs)	   
		       )) ;
  
  {
    /* match the pairs and report the candidate introns */
    /* we do 3 pass 
     * pass 0: match donors to acceptors 
     * pass 1: match donors to donors to detect inversion and transloc
     * pass 1: match acceptors to acceptors to detect inversion and transloc
     */
    
    Array trueDonors = donors ;
    Array trueAcceptors = acceptors ;
    int pass, type, bestType, bestDu, iiStrands ; 
 
    arrayDestroy (exportIntrons) ;
    exportIntrons = arrayHandleCreate (10000, IXP, h) ;
    truedBitSet = bitSetCreate (maxD, h) ;
    trueaBitSet = bitSetCreate (maxA, h) ;
    for (pass = 0 ; pass < 4 ; pass++)
      {
	int ddChrom, oldDdChrom ;
	BOOL ddDown, oldDdDown ;
	int ddA1 ;
	if (0 && pass > 0) continue ;
	switch (pass)
	  {
	  case 0: 
	    donors = trueDonors ;       dBitSet = truedBitSet ;
	    acceptors = trueAcceptors ; aBitSet = trueaBitSet ;
	    break ;
	  case 1: 
	    donors = trueDonors ;    dBitSet = truedBitSet ;
	    acceptors = trueDonors ; aBitSet = truedBitSet ;
	    break ;
	  case 2: 
	    donors = trueAcceptors ;    dBitSet = trueaBitSet ;
	    acceptors = trueAcceptors ; aBitSet = trueaBitSet ;
	    break ;
	  }
	maxD = arrayMax (donors) ;
	maxA = arrayMax (acceptors) ;

	ddChrom = oldDdChrom = -1 ; ja = 0 ;
	ddDown = oldDdDown = FALSE ;
	for (id = ja = 0 ; id < maxD ; id++)
	  {
	    if (0 && bit (aBitSet, id)) continue ;
	    dd = arrp (donors, id, DA) ; 
	    if (dd->n < min_support) continue ;

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
	    for (ia = ja, da = ia < maxA ? arrp (acceptors, ia, DA) : 0 ; ia < maxA ; da++, ia++)
	      {  
		if (0 && bit (aBitSet, ia)) continue ;
		if (da->n < min_support) continue ;
		if (pass < 3) /* same strand, same chrom */
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
		if (0 && dd->a1 > 4468 && dd->a1 < 4644338 && da->a1 > -1 && da->a1 < 4)
		  invokeDebugger () ;
		if (da->chrom != ddChrom)
		  type = sv->translocation ; /* TRANSLOCATION */
		else if (pass > 0)  /* same chrom same strand */
		  {
		      {
			if (da->a1 > ddA1 - 12 && da->a1 < ddA1 + 12)
			  type = sv->palindrome ;   /* PALINDROME */
			else
			  type = sv->inversion ; /* INVERSION */
		      }
		  }
		else  if (da->a1 > ddA1)  /* [pass == 0 and same chrom same strand */
		  {
		    if (da->a1 < ddA1 + microLength)
		       type = sv->microdeletion ; /* MICRODELETION */
		    else if (da->a1 < ddA1 + miniLength)
		       type = sv->minideletion ; /* MINIDELETION */
		    else 
		      type = sv->deletion ; /* DELETION */
		  }
		else if (da->a1 < ddA1)
		  {
		    if (da->a1 > ddA1 - microLength)
		       type = sv->microduplication ; /* MICRODUPLICATION */
		    else if (da->a1 > ddA1 - miniLength)
		       type = sv->miniduplication ; /* MINIDUPLICATION */
		    else
		      type = sv->duplication ; /* DUPLICATION */
		  }

		if (!type || type > bestType)
		  continue ;
		switch (type)
		  {
		  case 0: 
		    continue ;
		  case 9:  /* translocation */
		    break ;
		  default:
		    {
		      int du =  da->a1 - dd->a1 ;
		      if (du < 0) du = -du ;

		      if (0 && du > bestDu)
		      continue ;
		    }
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

		/* blunt end contact would be dx1 + dx2 = 4, one base gap would be 6 */
		if (0 && dx1 + dx2 < 4)
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
		      { dx1 = dx2 = -1 ; break ; }
		  }
		if (dx1 + dx2 < 0)
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
                if (type == sv->deletion)
		  {
		    for (i = 0 ; i <= dx1 - 2 ; i++)
		      {
			k1 = k2 = -1 ; 
			ccp = exon1 + i  ; ccq = ccp+1 ;
			f2[3] = *ccp ; f2[4] = *ccq ;
			f1[1] = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)*ccp])] ;
			f1[0] = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)*ccq])] ;
			if (
			    (*ccp == 'a' && *ccq == 'c')   /*    d1 == "gt"  */
			    )
			  k1 += 5 ;
			if (
			    (*ccp == 'g' && *ccq == 'c')     /* || d1 == "gc"  */
			    )
			  k1 += 4 ;
			else 
			  {
			    if (*ccp == 'a' && *ccq == 'g')   /*    d1 == "ct"  */
			      k2 += 4 ;  /*    d1 == "ag"  */
			  }
			
			ccp = exon2 + dx1 - i - 2 ; ccq = ccp+1 ;
			f1[3] = *ccp ; f1[4] = *ccq ;
			f2[1] = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)*ccp])] ;
			f2[0] = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)*ccq])] ;
			if (i == 0) /* always initialize */
			  { bestk = k1 ; bestdx = i ; reverse = FALSE ; strcpy (foot, f1) ; foot[5]='1'; foot[6]=0 ;}
			if (*ccp == 'a' && *ccq == 'g')   /*    d1 == "ag"  */
			  k1 += 4 ;
			else 
			  {
			    if (
				(*ccp == 'a' && *ccq == 'c')     /*    d1 == "ac"  */
				)
			      k2 += 5 ;
			    if (
				(*ccp == 'g' && *ccq == 'c')     /* || d1 == "gc"  */
				)
			      k2 += 4 ;
			  }
			if ( k1 > bestk)
			  { bestk = k1 ; bestdx = i ; reverse = FALSE ; strcpy (foot, f1) ; foot[5]='1'; foot[6]=0 ;}
			else if ( k2 > bestk)
			  { bestk = k2 ; bestdx = i ; reverse = TRUE ; strcpy (foot, f2) ; foot[5]='2'; foot[6]=0 ;} 
			else if (i == dx1 - 2 && ! dd->isDown)
			  { bestdx = dx1 - 2 ;  reverse = FALSE ; strcpy (foot, f1) ; }
		      }
		    if (bestk == 4) type = sv->deletion ;
		    if (bestk1 < bestk) bestk1 = 0 ;
		    if (bestk2 < bestk) bestk2 = 0 ;
		    if (0 && bestk1 && bestk2 && bestk < 4) bestk2 = 0 ; /* no ambiguous sliding non gt_ag */
		  }
		/* acceptable intron boundaries */
		for (iiStrands = 0 ; iiStrands < 1 ; iiStrands++)
		  {
		    DA *mydd = dd, *myda = da ;
		    
		    /*
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
		    */
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
			char c ;
			up->exon1 = myda->exon ;
			up->exon2 = mydd->exon ;
			c = up->foot[0] ;   up->foot[0] = up->foot[3] ; up->foot[3] = c ;
			c = up->foot[1] ;   up->foot[1] = up->foot[4] ; up->foot[4] = c ;
			up->foot[6]='R' ; up->foot[7] = 0 ;
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
		    bitSet (dBitSet, id) ;
		    bitSet (dBitSet, ia) ;
		  }
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
	  if (up->type == sv->circle) invokeDebugger () ;
	  switch (up->type)
	    {
	    case 1: /* MICRODELETION */
	    case 2: /* MICRODELETION */
	    case 3: /* DELETION */
	      break ;
	      if (up->u1 < up->u2)
		{ up->u1++ ; up->u2-- ; up->du = up->u2 - up->u1 + 1 ; }
	      else
		{ up->u1-- ; up->u2++ ; up->du = up->u1 - up->u2 + 1 ; }
	      break ;
	    case 4: /* CIRCLE */
	    case 5: /* MICRODUPLICATION */
	    case 6: /* MINIDUPLICATION */
	    case 7: /* DUPLICATION */
	      if (up->u1 < up->u2)
		{ up->u1++ ; up->u2-- ; up->du = up->u2 - up->u1 + 1 ; 
		
		  if (up->u1 == 1) 
		    up->type = sv->circle ; /* CIRCLE */
		}
	      else
		{ up->u1-- ; up->u2++ ; up->du = up->u1 - up->u2 + 1 ;
		  if (up->u2 == 1)
		    up->type = sv->circle ; /* CIRCLE */
		}
	      break ;
	    case 8: /* PALINDROME */
	    case 9: /* INVERSION */
	      if (up->u1 < up->u2)
		{ 
		  up->u1-- ; up->u2++ ; up->du = up->u2 - up->u1 + 1 ; 
		}
	      else
		{ 
		  up->u1++ ; up->u2-- ; up->du = up->u1 - up->u2 + 1 ; 
		}
	      break ;
	    case 10: /* TRANSLOCATION */
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
	
	if (0) aceOutf (ao, "exon1=%d exon2=%d\n", up->exon1,up->exon2) ;
	
	if (pass == 0 && up->type == 0)
	  {
	    if (! up->reverse)
	      aceOutf (ao, "%s\t%s\t%s\t%09d\t%09d\tDONOR:%s\t%d\tACCEPTOR:%s\t%d\t%s\tdx1=%d  bestdx=%d  ln=%d\n"
		       , dictName (sv->classeDict, up->type)
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
	      aceOutf (ao, "%s\t%s\t%s\t%09d\t%09d\tDONOR:%s\t%d\tACCEPTOR:%s\t%d\t%s\tdx1=%d  bestdx=%d  ln=%d \n"
		       , dictName (sv->classeDict, up->type)
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
	    if (1)
	      {
		YY *yy = arrayp (sw->yys, arrayMax (sw->yys), YY) ;
		yy->classe = up->type ;
		yy->chromA = up->chrom ;
		yy->chromB = up->chromA ;
		yy->a1 = up->u1 - 40 ;
		yy->a2 = up->u1 ;
		yy->b1 = up->u2 ;
		yy->b2 = up->u2 + 40 ;
		yy->width = up->du ;
		yy->downA = up->isDown ;
		yy->downB = up->isDownA ;
		yy->sPlus = up->ddn ;
		yy->sMinus = up->dan ;
		dictAdd (sw->typeDict, up->foot, &yy->type) ;
	      }
	    else
	      {
		aceOutf (ao, "%s\t%s\t%s\t%09d\t%s\t%s\t%09d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n"
			 , dictName (sv->classeDict, up->type)
			 , dictName (dict, up->chrom)	
			 , up->isDown ? "Forward" : "Reverse"
			 , up->u1
			 , dictName (dict, up->chromA)	
			 , up->isDownA ? "Forward" : "Reverse"
			 , up->u2
			 , up->du
			 , up->foot
			 , sdu->runName
			 , up->ddn
			 , up->dan
			 , exon1 + up->bestdx + 2
			 , exon2 + up->dx1 - up->bestdx  
			 ) ;
	      }
	  }
	vp = up ;
      }
  }

   wego_log (hprintf (0, "// svNewIntrons Verified %d introns, tested %d, attempted %d, min_support %d start"
			, nOk2, nOk1, nOk0, min_support)) ;


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
      aceOutf (ao, "SL%d\t%s\t%s\t%09d\t%09d\tDONOR:%s\t%d\tACCEPTOR:%s\t%d\tdx1=%d\n"
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
	    aceOutf (ao, "pA\t%s\t%s\t%09d\t%09d\tpA:%s\t%d\n"
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
  
 done:
  ac_free (h) ;
  return nOk2 ;
} /* svDeDuoParseOverhangs */

/*************************************************************************************/

static int svDeUnoParseIntrons (const SV *sv, SDU *sdu, SW *sw)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = 0  ;
  Array yys = sw->yys ;
  DICT *dict = sw->dict ;
  int nn = 0 ;
  int microLength = 36, miniLength = 1000 ;
  char *command ;
  YY dd, *yy ;  
 
  command = hprintf (h, "gunzip -c %s/tmp/PHITS_virus/%s/*.introns.gz  %s/tmp/PHITS_mito/%s/*.introns.gz"
		     , sv->inFileName ? sv->inFileName : "."
		     , sdu->runName
		     , sv->inFileName ? sv->inFileName : "."
		     , sdu->runName) ; 
  ai = aceInCreateFromPipe (command, "r", 0, h) ;
    /*  ai = aceInCreate ("toto18", 0, h) ; */
  aceInSpecial (ai, "\n") ;
  wego_log (hprintf(h, "## client %d, processing %s\n", sv->clientId, aceInFileName (ai))) ;

  while (aceInCard (ai))
    {
      char *cp = aceInWord (ai) ;
      if (! cp || *cp == '#')
	continue ;

       /* jump column 1 */
      aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 2 */

      aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* col 3/4/5 first position */
      dictAdd (dict, cp, &dd.chromA) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &dd.a1) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &dd.a2) ; 

      aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* col 6/7/8 second position */
      dictAdd (dict, cp, &dd.chromB) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &dd.b1) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &dd.b2) ;
      dd.downA = (dd.a1 < dd.a2) ? 1 : -1 ; 
      dd.downB = (dd.b1 < dd.b2) ? 1 : -1 ; 

      aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* column 9: type */
      dictAdd (sw->typeDict, cp, &dd.type) ; 
      aceInStep (ai, '\t') ; aceInInt (ai, &dd.width) ;       /* col 10 */
      aceInStep (ai, '\t') ; aceInInt (ai, &dd.sPlus) ;  /* col 11 */
      dd.sMinus = dd.sPlus ;
      if (dd.sPlus < sv->min_support) 
	continue ;
      dd.width = 0 ;
  
      if (dd.downA == -1 && dd.downB == -1)
	svSwitchComplement (sw->typeDict, &dd) ;
      
      if (dd.chromA != dd.chromB)
	dd.classe = sv->translocation ;
      
      else if (dd.downA  != dd.downB)
	{
	  dd.classe = sv->inversion ; 
	  if (dd.a2 > dd.b1)
	    svSwitchComplement (sw->typeDict, &dd) ;
	  dd.width = dd.b1 - dd.a2 + 1 ;
	  if (
	      (dd.b2 > dd.a2 - 12 && dd.b2 < dd.a2 + 12) ||
	      (dd.b1 > dd.a1 - 12 && dd.b1 < dd.a1 + 12)
	      )
	    dd.classe = sv->palindrome ;   /* PALINDROME */
	}
      else  /* same chrom both down == 1 */
	{
	  if (dd.b1 > dd.a2)
	    {
	      dd.a2++ ; dd.b1-- ; dd.width = dd.b1 - dd.a2 + 1 ;  
	      if (dd.width < microLength) dd.classe =  sv->microdeletion ;
	      else if  (dd.width < miniLength) dd.classe = sv->minideletion ;
	      else dd.classe = sv->deletion ;
	    }	      
	  else
	    {
	      dd.width = dd.a2 - dd.b1 + 1 ;  
	      if (dd.type == sw->circleFeet)
		dd.classe = sv->circle ; 
	      else if (dd.width < microLength) dd.classe =  sv->microduplication ;
	      else if  (dd.width < miniLength) dd.classe = sv->miniduplication ;
	      else dd.classe = sv->duplication ;
	    }
	}
  
      yy = arrayp (yys, nn++, YY) ;
      *yy = dd ;
    }
  
  nn = svYyCompress (yys, sv->min_support) ; /* merges multiple supports of same breakpoint, returns new arrayMax (yys) */
  wego_log (hprintf (0, "## client %d, found %d breakpoints in run %s de_uno", sv->clientId, nn, sdu->runName)) ;
  ac_free (h) ;
  return nn ;
} /* svDeUnoParseIntrons */

/*************************************************************************************/

static void  svWiggleReport (const SV *sv, SDU *sdu, SW *sw, int unoDuo)
{
  int nn = 0, ii, i ;
  int min = sv->min_support ;
  YY *yy ; 
  DICT *dict = sw->dict ;
  Array yys = sw->yys ;
  char *command = hprintf (sw->h, ".%s.de_%s.txt", sdu->runName, unoDuo == 1 ? "uno" : "duo") ;
  ACEOUT ao = aceOutCreate (sv->outFileName, command, sv->gzo, sw->h) ;

  aceOutDate (ao, "###", hprintf (sw->h, "de-%s breakpoint analysis",  unoDuo == 1 ? "uno" : "duo")) ;

  aceOutf (ao, "## Requested minimal_support %d\n",  min) ;
  aceOutf (ao, "# Breakpoint type\tLeft element sequence\tStrand\tCoordinate\tRight element sequence\tStrand\tCoordinate\tReaaranged fragment length\tSequence at boundary\tRun\tNumber of sequences supporting the left element\tNumber of sequences supporting the right element") ;

  aceOutf (ao, "\tLeft break point\tRight breakpoint") ;

  for (ii = 1 ; ii <= 2 ; ii++)
    {
      aceOutf (ao, "\t") ;
      for (i = 0 ; i < BBWIDTH ; i++)
	aceOutf (ao, "\tB%d -%d", ii,  (i - BBWIDTH) * BBSTEP) ;
      aceOutf (ao, "\t") ;
      for (i = 0 ; i < BBWIDTH ; i++)
	aceOutf (ao, "\tB%d -%d", ii,  (i + 1) * BBSTEP) ;
    }
  
  nn = arrayMax (yys) ; /* superfluous, but makes the code more sturdy */
  /* verify that we have some anomalous reads */
   if (sv->de_wiggle)
     for (ii = 0, yy = arrp (yys, 0, YY) ; ii < nn ; ii++, yy++)
       {
	 BB *bb0, *bb1, *bb2, *bb3 ;
	 int i, np = 0, nok = 0 ;

	 bb0 = arrp (sw->bbs, yy->bb0, BB) ;	
	 bb1 = arrp (sw->bbs, yy->bb1, BB) ;	
	 bb2 = arrp (sw->bbs, yy->bb2, BB) ;	
	 bb3 = arrp (sw->bbs, yy->bb3, BB) ;	
	
	 for (i = 0 ; i < BBWIDTH ; i++)	
	   {
	     nok += bb0->ok[i] + bb1->ok[i] + bb2->ok[i] + bb3->ok[i] ; 
	     np += bb0->pp[i] + bb1->pp[i] + bb2->pp[i] + bb3->pp[i] ; 
	   }
	 if (20 * np < nok)
	   {
	     yy->sPlus *= -1 ;
	     yy->sMinus *= -1 ;
	   }
       }

  for (ii = 0, yy = arrp (yys, 0, YY) ; ii < nn ; ii++, yy++)
    {
      int ok ;
      char *subClasse[] = {"/multi mapped","/anomalous or partial","/uniquely mapped" } ;
      for (ok = (sv->de_wiggle ? 2 : 0)  ; ok >= 0 ; ok--)
	if (yy->sPlus + yy->sMinus >= min)
	  {
	    aceOutf (ao, "\n%s%s\t%s\t%s\t%09d\t%s\t%s\t%09d\t%d\t%s\t%s\t%d\t%d"
		     , dictName (sv->classeDict, yy->classe)
		     , sv->de_wiggle ? subClasse[ok] : ""
		     , dictName (dict, yy->chromA), (yy->downA == 1 ? "Forward" : "Reverse"), yy->a2
		     , dictName (dict, yy->chromA), (yy->downB == 1 ? "Forward" : "Reverse"), yy->b1
		     , yy->width
		     , yy->type ? dictName (sw->typeDict, yy->type) : "-"
		     , sdu->runName
		     , yy->sPlus 
		     , yy->sMinus 
		     ) ;
	    if (sv->de_wiggle)
	      {
		BB *bb ;
		
		bb = arrp (sw->bbs, yy->bb0, BB) ;	
		aceOutf (ao, "\t%d:%d", bb->a1, bb->a2) ;
		bb = arrp (sw->bbs, yy->bb1, BB) ;	
		aceOutf (ao, "::%d:%d", bb->a1,bb->a2) ;
		bb = arrp (sw->bbs, yy->bb2, BB) ;	
		aceOutf (ao, "\t%d:%d", bb->a1, bb->a2) ;
		bb = arrp (sw->bbs, yy->bb3, BB) ;	
		aceOutf (ao, "::%d:%d", bb->a1,bb->a2) ;
		
		bb = arrp (sw->bbs, yy->bb0, BB) ;
		aceOutf (ao, "\t") ;
		for (i = 0 ; i < BBWIDTH ; i++)
		  aceOutf (ao, "\t%.1f", ok == 2 ? bb->ok[i] : (ok == 1 ? bb->pp[i] : bb->nu[i])) ;
		
		bb  = arrp (sw->bbs, yy->bb1, BB) ;
		aceOutf (ao, "\t") ;
		for (i = 0 ; i < BBWIDTH ; i++)
		      aceOutf (ao, "\t%.1f", ok == 2 ? bb->ok[i] : (ok == 1 ? bb->pp[i] : bb->nu[i])) ;
		
		bb  = arrp (sw->bbs, yy->bb2, BB) ;
		aceOutf (ao, "\t") ;
		for (i = 0 ; i < BBWIDTH ; i++)
		  aceOutf (ao, "\t%.1f", ok == 2 ? bb->ok[i] : (ok == 1 ? bb->pp[i] : bb->nu[i])) ;
		
		bb  = arrp (sw->bbs, yy->bb3, BB) ;
		aceOutf (ao, "\t") ;
		for (i = 0 ; i < BBWIDTH ; i++)
		  aceOutf (ao, "\t%.1f", ok == 2 ? bb->ok[i] : (ok == 1 ? bb->pp[i] : bb->nu[i])) ;
	      }
	  }
    }
  aceOutf (ao, "\n") ;
} /* svWiggleReport */

/*************************************************************************************/

static int svDeUnoDo (const SV *sv, SDU *sdu)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
  SW sw ;

  memset (&sw, 0, sizeof (SW)) ;
  sw.h = h ;
  sw.dict = dictHandleCreate (200, h) ;  
  sw.typeDict = dictHandleCreate (200, h) ;
  dictAdd (sw.typeDict, "--_--", &sw.circleFeet) ; 

  sw.yys = arrayHandleCreate (10000, YY, h) ;

  nn = svDeUnoParseIntrons (sv, sdu, &sw) ;
  if (sv->de_wiggle && nn) 
    {
      svYy2Bb (sv, sdu, &sw) ; /* create the bb, and a sorted set of pointers */
      svWiggleDo  (sv, sdu, &sw) ; /* grab bb coverage from the wiggle file */
      svWiggleDo  (sv, sdu, &sw) ; /* grab wild type coverage covering the break point */
      svWiggleReport (sv, sdu, &sw, 1) ;
    }

  ac_free (h) ;
  return nn ;
} /* svDeUnoDo  */

/*************************************************************************************/

static int svDeDuoDo (const SV *sv, SDU *sdu)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
  SW sw ;

  memset (&sw, 0, sizeof (SW)) ;
  sw.h = h ;
  sw.dict = dictHandleCreate (200, h) ;
  sw.typeDict = dictHandleCreate (200, h) ;
  sw.chromDict = dictHandleCreate (200, h) ;
  dictAdd (sw.typeDict, "--_--", &sw.circleFeet) ; 
  sw.yys = arrayHandleCreate (10000, YY, h) ;

  nn = svDeDuoParseOverhangs (sv, sdu, &sw) ;
  if (sv->de_wiggle && nn) 
    {
      svYy2Bb (sv, sdu, &sw) ; /* create the bb, and a sorted set of pointers */
      svWiggleDo  (sv, sdu, &sw) ; /* grab bb coverage from the wiggle file */
      svWiggleReport (sv, sdu, &sw, 2) ;
    }

  ac_free (h) ;
  return nn ;
} /* svDeDuoDo  */

/*************************************************************************************/
/*************************************************************************************/
/* analysing breakpoints from pairs is planed in several steps
 * 0) select a library length, ln, and declare incompatible pairs 2*da > 3 ln
 * 0) split genome in blocks of length da
 * 1) scan the hits files, cumulate compatible and uncompatible pairs in each block
 * 2) flag block with high count and high ratio of uncompatible pairs
 * 3) for the flag blocks, hash them and construct a weighted list of at most 6 friends, including NULL friends
 * 4) for spotted pairs reaasemble the left and right overhang concensus
 * 5) overlap the concesus, giving the seqeunce of the hole, realign locally
 */
typedef struct paribusStruct {
  AC_HANDLE h ;
  int nBins ;
  int dmaxOk ; /* maximal distance to declare a pair is ok, 2*dmaxOk implies incompatible pair */ 
  Array aaa ; /* contains pLeft, pRight, pLeftOk, pRightOk  number of oriented pairs in each genome bin */
  Array badPairs ;
  BitSet pFlags, nFlags ;
  DICT *chromDict ;
} PBS ;

typedef struct paribusHitStruct {
  int chromA, a1, a2, chomB, b1, b2 ;
} PH ;



static PBS* svDeParibusInit (AC_HANDLE h)
{
  PBS *pbs = halloc (sizeof(PBS), h) ;
  int i ;
  memset (pbs, 0, sizeof(PBS)) ;

  pbs->dmaxOk = 300 ;
  pbs->nBins = 100000 ;
  pbs->aaa = arrayHandleCreate (256, Array, h) ;
  for (i = 0 ; i < 256 ; i++)
    array (pbs->aaa, i, Array) =  arrayHandleCreate (pbs->nBins, int, h) ;

  pbs->badPairs = arrayHandleCreate (100000, PH, h) ;

  pbs->pFlags = bitSetCreate (pbs->nBins, h) ;
  pbs->nFlags = bitSetCreate (pbs->nBins, h) ;
  pbs->chromDict = dictHandleCreate (256, pbs->h) ;
  return pbs ;
}

/*************************************************************************************/
/* phase 1: scan the hits files, cumulate compatible and uncompatible pairs in each block */
static void svDeParibusCountPairs (const SV *sv, SDU *sdu, PBS *pbs)
{
  AC_HANDLE h = ac_new_handle (), h0 = pbs->h ;
  int  dmaxOk = pbs-> dmaxOk, nBadP = 0, nBadN = 0, nOk = 0 ;
  char *command ;
  Array aa, aaa = pbs->aaa ;
  ACEIN ai = 0  ;
  command = hprintf (h, "gunzip -c %s/tmp/COUNT/%s/f2.*.hits.gz"
		      , sv->inFileName ? sv->inFileName : "."
		     , sdu->runName) ;

  ai = aceInCreateFromPipe (command, "r", 0, h) ; 
  if (ai)
     {
       aceInSpecial (ai, "\n") ;
       while (aceInCard (ai))
	 {
	   int mult, tMult, a1, a2, da, chrom, ali, toAli, type ;
	   char *cp = aceInWord (ai) ;
	   if (! cp || *cp == '#')
	     continue ;
	   /* jump column 1 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 2 */
	   aceInStep (ai, '\t') ; aceInInt (ai, &mult) ;
	   aceInStep (ai, '\t') ; aceInInt (ai, &toAli) ;
	   aceInStep (ai, '\t') ; aceInInt (ai, &ali) ;
	   
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 6 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 7 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* traget_class */
	   if (strstr (cp, "pike") || strstr (cp, "ecoy") ||  strstr (cp, "rrna")) 
	     continue ;
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 9 */
	   
	   aceInStep (ai, '\t') ; aceInInt (ai, &tMult) ; if (tMult > 1) continue ;
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* traget */
	   dictAdd (pbs->chromDict, cp, &chrom) ;
	   aceInStep (ai, '\t') ; aceInInt (ai, &a1) ;
	   aceInStep (ai, '\t') ; aceInInt (ai, &a2) ;
	   
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 14 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 15 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 16 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 17 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 18 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 19 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 20 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 21 */
	   
	   aceInStep (ai, '\t') ; aceInInt (ai, &da) ;
	   switch (da)
	     {          /* synchronize to hack with bestali.c snp.c and wiggle.c and sv.c */
	     case -7: /* links to RNA */
	     case -11: /* links to spike in */
	       continue ;
	     case -2: /* gene to close by genome */
	     case -3: /* too distant inside a single gene  NA: the distance given is the true disstance even if big */
	     case -5: /* 2 transcripts of same gene */
	     case -10: /* close genes in cis */
	       type = (a1 < a2 ? 0 : 1) ; /* good pair */
	       break ;
	     case -1: /* orphan */
	     case -4: /* too far on genome, NA: the distance given is the true disstance even if big */
	     case -6: /* distant genes */
	     case -8: /* gene mito */
	     case -9: /* bad gene genome */
	     case -12:
	     case -13:
	     case -14: /* bad topology */
	       type = (a1 < a2 ? 2 : 3) ;  /* incompatible pair */
	       break ;
	       
	     default:
	       if (da > 0 || da < -15 ) /* good pair */
		 type = (a1 < a2 ? 0 : 1) ;
	       else
		 type = (a1 < a2 ? 2 : 3) ;  /* incompatible pair */
	       if (da < 0) da = -da ;
	       if (da > dmaxOk && da < 2 * dmaxOk) continue ;
	       break ;
	     }
	   a2 = a2/dmaxOk ;
	   aa = array (aaa, 4 * chrom + type, Array) ;
	   if (! aa)
	     aa = array (aaa, 4 * chrom + type, Array) = arrayHandleCreate (100000, int, h0) ;
	   array (aa, a1/dmaxOk, int) += mult ;
	   if (type == 3)
	     nBadP += mult ;
	   else if (type == 2)
	     nBadN += mult ;
	   else
	     nOk += mult ;
	 }
     }
  wego_log (hprintf (0, "found %d %d incompatible reads, %d ok in file %s\n", nBadP, nBadN, nOk, aceInFileName(ai) ? aceInFileName(ai): "-")) ;
  ac_free (h) ;
} /* svDeParibusCountPairs */

/*************************************************************************************/
/* phase 2: flag block with high count and high ratio of uncompatible pairs */
static void svDeParibusFlagBlocks (const SV *sv, SDU *sdu, PBS *pbs) 
{
  Array aaa = pbs->aaa, aa[4] ;
  BitSet pFlags = pbs->pFlags ;
  BitSet nFlags = pbs->nFlags ;
  int ii, i,  nn, x, y, nok = 0 ;
  int min_support = sv->min_support ;

  for (ii = 0 ; ii < arrayMax (aaa) ; ii += 4)
    {
      aa[0] = array (aaa, ii, Array) ;
      aa[1] = array (aaa, ii+1, Array) ;
      aa[2] = array (aaa, ii+2, Array) ;
      aa[3] = array (aaa, ii+3, Array) ;
      if (! aa[0] || ! aa[1] || ! aa[2] || ! aa[3])
	continue ;

      for (nn = i = 0 ; i < 4 ; i++)
	if (nn < arrayMax (aa[i])) 
	  nn = arrayMax (aa[i]) ;
      for (i = 0 ; i < nn ; i++)
	{
	  x = i < arrayMax(aa[0]) ? arr (aa[0], i, int) : 0 ;  /* 5' good read */
	  y = i < arrayMax(aa[2]) ? arr (aa[2], i, int) : 0 ;  /* 5' incompatible read */
	  if (y > min_support && y > 2*x)
	    {
	      bitSet (pFlags, i) ;  /* il en faut un par chrom */
	      fprintf (stderr, "pFlag %s %d :: %d %d\n", dictName(pbs->chromDict, ii/4), i* pbs->dmaxOk, x, y) ;
	    }
	  else
	    nok++ ;
	  x = i < arrayMax(aa[1]) ? arr (aa[1], i, int) : 0 ;  /* 3' good read */
	  y = i < arrayMax(aa[3]) ? arr (aa[3], i, int) : 0 ;  /* 3' incompatible read */
	  if (y > min_support && y > 2*x)
	    {
	      bitSet (nFlags, i) ;   
	      fprintf (stderr, "nFlag %s %d :: %d %d\n", dictName(pbs->chromDict, ii/4), i* pbs->dmaxOk, x, y) ;
	    }
	  else
	    nok++ ;
	}
    }
  wego_log (hprintf (0, "svDeParibusFlagBlocks flagged %d and %d zones, found %d ok zones in run %s", 
		     bitSetCount (pFlags) , bitSetCount (nFlags), nok, sdu->runName
		     )
	    ) ;
} /* svDeParibusFlagBlocks  */

#ifdef JUNK
/*************************************************************************************/
/* phase 4: for spotted pairs reaasemble the left and right overhang concensus */
static void svDeParibusLocateBlockPairs (SV *sv, SDU *sdu, PBS *pbs) 
{
} /* svDeParibusLocateBlockPairs */

/*************************************************************************************/
/* phase 5: overlap the concesus, giving the seqeunce of the hole, realign locally */
static void svDeParibusLocateBreakPoints (SV *sv, SDU *sdu, PBS *pbs) 
{
} /* svDeParibusLocateBreakPoints */

/*************************************************************************************/
#endif

static int svDeParibusDo (const SV *sv, SDU *sdu)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
  PBS *pbs =  svDeParibusInit (h) ;
  
  wego_log (hprintf (0, "svDeParibusDo id=%d %s start"
		      , sdu->clientId, sdu->runName
		  )
	 ) ;
    svDeParibusCountPairs (sv, sdu, pbs) ;
    svDeParibusFlagBlocks (sv, sdu, pbs) ;
  /*
     svDeParibusLocateBlockPairs (sv, sdu, pbs) ;
    svDeParibusLocateBreakPoints (sv, sdu, pbs) ;
  */

    /*   aceOutf (ao, "\n") ; */
  return nn ;
} /* svDeParibusDo  */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/* phase 5: overlap the concesus, giving the seqeunce of the hole, realign locally */
static void svDeAvCountGenePairsLane (const SV *sv, SDU *sdu, SW *sw, const char *fNam)
{
  AC_HANDLE h = ac_new_handle ()  ;
  KEYSET ks = sw->genePairKs ;
  int nnn = 0 ;
  char geneBuf1[100], geneBuf2[100], mrnaBuf[100], readBuf[100], oldReadBuf[100] ;
  char *command ;
  ACEIN ai ;

  if (1)
    {
      command = hprintf (h, "gunzip -c %s", fNam) ;
      if (0) command = "cat totoLSAMP" ;
      wego_log ( command) ;
      ai = aceInCreateFromPipe (command, "r", 0, h) ; 
    }
  else
    {
      command = hprintf (h, "gunzip -c %s > %s.txt", fNam, fNam) ;
      wego_log ( command) ;
      system (command) ;
      command = hprintf (h, "%s.txt", fNam) ;
      ai = aceInCreateFromFile (command, "r", 0, h) ; 

    }
  oldReadBuf[0] = 0 ;
  if (ai)
     {
       aceInSpecial (ai, "\n") ;
       while (nnn < 10000 && aceInCard (ai))
	 {
	   int mult, tMult, a1, a2, da, ali, toAli, k, genePair ;
	   char *cp = aceInWord (ai) ;
	   if (! cp || *cp == '#')
	     continue ;
	   if (! strstr (aceInPos (ai), "links_to:"))
	     continue ;
	   if (!strcmp (oldReadBuf, cp)) continue ;
	   strncpy (readBuf, cp, 99) ;
	   /* jump column 1 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 2 */
	   aceInStep (ai, '\t') ; aceInInt (ai, &mult) ;
	   aceInStep (ai, '\t') ; aceInInt (ai, &toAli) ;
	   aceInStep (ai, '\t') ; aceInInt (ai, &ali) ;
	   
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 6 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 7 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* target_class */
	   if (strcmp (cp, "ET_av")) continue ;                              /* select av */
	   strncpy (oldReadBuf, readBuf, 99) ;
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; 
	   strncpy (geneBuf1, cp, 99) ;
	   aceInStep (ai, '\t') ; aceInInt (ai, &tMult) ; if (tMult > 1) continue ;
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* target */
	   strncpy (mrnaBuf, cp, 99) ;

	   aceInStep (ai, '\t') ; aceInInt (ai, &a1) ;
	   aceInStep (ai, '\t') ; aceInInt (ai, &a2) ;
	   
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 14 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 15 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 16 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 17 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 18 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 19 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 20 */
	   aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (!cp) continue ; /* jump column 21 */
	   
	   aceInStep (ai, '\t') ; da = 0 ; aceInInt (ai, &da) ;
	   if (da != -6) /* distant genes */
	     continue ;
	   cp = aceInWord (ai) ; if (!cp) continue ;
	   cp = strstr (cp, "links_to:") ;  if (!cp) continue ;
	   cp += 9 ;
	   strncpy (geneBuf2, cp, 99) ;
	   k = strcmp (geneBuf1, geneBuf2) ;
	   if (! k)
	     continue ;
	   {
	     char buf[200] ;
	     sprintf (buf, "%s\t%s" 
		      , k < 0 ? geneBuf1 : geneBuf2
		      , k < 0 ? geneBuf2 : geneBuf1
		      ) ;
	     dictAdd (sw->genePairDict, buf, &genePair) ;
	   }
	   keySet (ks, genePair) += mult ; nnn += mult ;
	 }
     }
  wego_log (hprintf (0, "found %d supports for %d gene pairs in file %s in run %s de_av\n", keySetMax (ks), nnn, aceInFileName(ai), sdu->runName)) ;
  if (0)
    {
      command = hprintf (h, "%s.txt", fNam) ;
      unlink (command) ;
    }
  ac_free (h) ;
} /* svDeAvCountGenePairsLane */

/*************************************************************************************/

static void svDeAvCountGenePairs (const SV *sv, SDU *sdu, SW *sw)
{
  AC_HANDLE h = ac_new_handle ()  ;
  char *command, *cp ;
  ACEIN ai ;

  command = hprintf (h, "ls %s/tmp/COUNT/%s/f2.*.hits.gz"
		     , sv->inFileName ? sv->inFileName : "."
		     , sdu->runName) ;
  ai = aceInCreateFromPipe (command, "r", 0, h) ;
  while (aceInCard (ai))
    {
      if ((cp = aceInWord (ai)))
	svDeAvCountGenePairsLane (sv, sdu, sw, cp) ;
    }
 ac_free (h) ;
} /* svDeAvCountGenePairs */

/*************************************************************************************/

static void svDeAvReportGenePairs (const SV *sv, SDU *sdu, SW *sw)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao ;
  char *command ;
  DICT *dict = sw->genePairDict ;
  KEYSET ks = sw->genePairKs ;
  int n, i, iMax = keySetMax (ks) ;
  int min = sv->min_support ;

  command = hprintf (h, ".%s.gene_pairs.txt", sdu->runName) ;
  ao = aceOutCreate (sv->outFileName, command, sv->gzo, h) ;
  aceOutDate (ao, "###", "de-av breakpoint analysis") ;
  aceOutf (ao, "## %s Requested minimal_support %d\n",  sdu->runName, sv->min_support) ;

  for (i = 0 ; i < iMax ; i++)
    {
      n = keySet (ks, i) ;
      if (n >= min)
	aceOutf (ao, "%s\t%d\n", dictName (dict, i), n) ;
    }

  aceOutDate (ao, "###", "de-av breakpoint analysis done ") ;
  ac_free (h) ;
} /* svReportGenePairs */

/*************************************************************************************/

static int svDeAvDo (const SV *sv, SDU *sdu)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
   SW sw ;

  memset (&sw, 0, sizeof (SW)) ;
  sw.h = h ;
  sw.genePairDict = dictHandleCreate (200, h) ;
  sw.genePairKs = keySetHandleCreate (h) ;
  sw.yys = arrayHandleCreate (10000, YY, h) ;

  wego_log (hprintf (0, "svDeAvDo id=%d %s start"
		      , sdu->clientId, sdu->runName
		  )
	 ) ;
  if (1)
    {
      svDeAvCountGenePairs (sv, sdu, &sw) ;
      svDeAvReportGenePairs (sv, sdu, &sw) ;
    }
  sleep (1) ;
  return nn ;
} /* svDeAvDo  */

/*************************************************************************************/
/*************************************************************************************/

static void svClientLoop (const void *vp)
{
  const SV *sv = (SV *) vp ;
  int nn = 0 ;
  SDU sdu ;

  wego_log (" svClientLoop") ;
  memset (&sdu, 0, sizeof (SDU)) ;
  /* get the request */
  while (channelMultiGet (sv->todo, &sdu, 1, SDU))
    {
      sdu.clientId = sv->clientId ;
      if (sv->de_uno) 
	nn += svDeUnoDo (sv, &sdu) ; 
      if (sv->de_duo) 
	nn += svDeDuoDo (sv, &sdu) ; 
      if (sv->de_paribus) 
	nn += svDeParibusDo (sv, &sdu) ; 
       if (sv->de_av) 
	nn += svDeAvDo (sv, &sdu) ; 
      channelPut (sv->done, &sdu, SDU) ;
    }
  
  wego_log (hprintf (0, "svDeUnoClientSide %d done", sv->clientId)) ;
  return ;
} /*  svClientLoop */

/*************************************************************************************/

static int svDeAnyLoop (SV *sv)
{
  int nn = 0 ;
  int run ;
  SDU sdu ;

  /* send the requests */
  wego_log ("svDeAnyLoop") ;
  for (run = 0 ; run < sv->nRuns ; run++)
    {
      memset (&sdu, 0, sizeof (SDU)) ;      
      strncpy (sdu.runName, dictName (sv->runDict, run + 1), 256) ;
      channelPut (sv->todo, &sdu,  SDU) ;
    }
  channelClose (sv->todo) ;

  /* wait for completion */ 
  if (! sv->clientId)
    {
      wego_log (" svDeAnyLoop waiting for completion") ;
      for (run = 0 ; run < sv->nRuns ; run++)
	{
	  channelMultiGet (sv->done, &sdu, 1, SDU) ;
	  wego_log (hprintf (0, "number %d done", run)) ;
	}
    }
  wego_log ("svDeAnyLoop done") ;
  return nn ;
} /*  svDeAnyLoop */

/*************************************************************************************/

static void svDeAny (SV *sv)
{
  BOOL debug = FALSE ;

  if (sv->serverSide)
    { 
      AC_HANDLE h = ac_new_handle () ;
      int ii, nn = sv->nRuns ;  
      int  max_tasks = nn < sv->max_tasks ? nn : sv->max_tasks ;

      sv->todo = taskCreateWriterChannel (sv->serverTask, "de_any_request", nn, SDU) ;
      sv->done = taskCreateReaderChannel (sv->serverTask, "de_any_reply", nn, SDU) ;
      if (debug)
	{
	  channelDebug ( sv->todo, TRUE, "todo") ;     
	  channelDebug ( sv->done, TRUE, "done") ;     
	}
      if (! sv->inFileName) sv->inFileName = "." ;
      if (! sv->outFileName) sv->outFileName = "." ;
      for (ii = 0 ; ii < max_tasks ; ii++)
	taskDispatch (0, hprintf (h, "%s %s %s %s %s --clientId %d --min_support %d --overhangLength %d -i %s -o %s"
				, sv->de_uno ? "--de_uno" : ""
				, sv->de_duo ? "--de_duo" : ""
				, sv->de_paribus ? "--de_paribus" : ""
				, sv->de_av ? "--de_av" : ""
				, sv->de_wiggle ? "--de_wiggle" : ""
				, ii + 1
				, sv->min_support
				, sv->overhangLength
				, sv->inFileName, sv->outFileName
				)
		    , h) ;
      
      svDeAnyLoop (sv) ;
      
      /* all done */
      ac_free (h) ;
    }
  else if (sv->clientId)
    {
      int nn = 1 ;
      wego_log (hprintf (0, "client %d want to crate channel nn=%d", sv->clientId, nn)) ;
      sv->todo = taskCreateReaderChannel (sv->clientTask, "de_any_request", nn, SDU) ;
      sv->done = taskCreateWriterChannel (sv->clientTask, "de_any_reply", nn, SDU) ;
      if (sv->todo && sv->done)
	{
	  if (debug)
	    {
	      channelDebug ( sv->todo, TRUE, "todo") ;     
	      channelDebug ( sv->done, TRUE, "done") ;     
	    }
	  svClientLoop (sv) ;
	  channelClose (sv->done) ;
	}
    }
  else
    {
      int ii, nn = sv->nRuns ;  
      int  max_threads = nn < sv->max_threads ? nn : sv->max_threads ;

      sv->todo = channelCreate (nn, SDU, 0) ; /* do not allocate on sv->h, because the main code may finish before the go routine */
      sv->done = channelCreate (nn, SDU, 0) ;

      /* create wego threads */
      for (ii = 0 ; ii < max_threads ; ii++)
	wego_go (svClientLoop, sv, SV) ;
      
      svDeAnyLoop (sv) ;
    }
} /*  svDeAny */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: sv -de_uno -... \n"
	   "//      try: -h --help --hits_file_caption \n"
	   "// Objective: \n"
	   "//     Analysize the .introms and .overhangs files exportred by clipalign\n"
	   "//     Report de_uno and de_duo breakpoint candidates\n"
	   "//     Validate the candidates by local reextension\n"
	   "//     Measure the reulting support and branchig ratios of each break point\n "
	   "// General Options:\n"
	   "//    -o outFilePrefix : output file prefix, all exported files will start with that prefix\n"
  	   "//    --gzo : the output file is gziped\n"
	   "//    -i rootDirName\n"
	   "// Actions:\n"
	   "//      --runs <file name> : list of runs, one run name per line, i.e. MetaDB/RunList\n"
	   "//      --min_support <int> : [default 3] reject breakwpoint with insufficient support\n"
	   "//      --de_uno : analyze the $rootDirName/PHITS_virus/$run/*.introns.gz files\n"
	   "//                 export $outFilePrefix/d1.$run.de_uno.txt\n"
	   "//      --de_duo : analyze the $rootDirName/PHITS_virus/$run/*.overhangs.gz files\n"
	   "//                 export $outFilePrefix/d1.$run.de_uno.txt\n"
	   "//            --seedLength n1 --overhangLength n2\n"
	   "//               n1 {default=25] is the required number of bases matching the genome\n"
	   "//               n2 [default= 8] is the required number of overhanging bases\n"
	   "//      --de_paribus : analyze the $rootDirName/tmp/COUNT/$run/*.hits.gz files\n"
	   "//               locate breakpoints from clusters incompatible pairs\n"
	   "//  Multi-threading:\n"
	   "//      --max_threads <int> : default 4\n"
	   "//          number of concurent threads running on the server side\n"
	   "//          this number should be at least 3, because the server manages the\n"
	   "//           communication layers as parallel threads\n"
	   "//  Multi_tasking:\n"
	   "//      By default, the program runs as a single process, however\n"
	   "//      the present program may also spawns copies of itself as parallel tasks\n"
	   "//      synchronized via the acedb remote-channels library\n"
	   "//      --multi_tasking <int> : allows the program to spawn client subtasks\n" 
	   "//            specifing the max number of tasks submitted in parallel\n"
	   "//      [+rChan_Id <int> +rChan_h my_host +rChan_p 123452 +rChan_s 65432 _tChan_c 0] : do not use\n"
	   "//            reserved parameters used in multi-tasking interprocess communications\n"  
	   "//      --local: [default]  the subtasks run on the local computer, using a system() call\n"
	   "//      --submit : the subtasksrun on the compute farm, provided scripts/submit is properly configured\n"

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
  SV sv ;
  AC_HANDLE h = 0 ;
  char commandBuf [4000] ;

  messErrorInit (argv[0]) ;

  h = ac_new_handle () ;
  memset (&sv, 0, sizeof (SV)) ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
  sv.h = h ;

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
  getCmdLineOption (&argc, argv, "-i", &(sv.inFileName)) ;
  getCmdLineOption (&argc, argv, "-o", &(sv.outFileName)) ;
  sv.gzo = getCmdLineOption (&argc, argv, "--gzo", 0) ;

  sv.de_uno = getCmdLineOption (&argc, argv, "--de_uno", 0) ;
  getCmdLineOption (&argc, argv, "--runs", &(sv.runsFileName)) ;
  sv.min_support = 3 ; /* default */
  getCmdLineInt (&argc, argv, "--min_support", &(sv.min_support)) ;
 
  sv.de_av = getCmdLineOption (&argc, argv, "--de_av", 0) ;

  sv.de_duo = getCmdLineOption (&argc, argv, "--de_duo", 0) ;
  sv.seedLength = 25 ; sv.overhangLength = 10 ; /* 10, we want to be sure we correctly jumped a single base indel */
  getCmdLineInt (&argc, argv, "--seedLength", &sv.seedLength) ;
  getCmdLineInt (&argc, argv, "--overhangLength", &sv.overhangLength) ;

  sv.de_paribus = getCmdLineOption (&argc, argv, "--de_paribus", 0) ;
  sv.de_wiggle = getCmdLineOption (&argc, argv, "--de_wiggle", 0) ;

  sv.max_threads = 1 ;
  getCmdLineInt (&argc, argv, "--max_threads", &(sv.max_threads)) ;
  getCmdLineInt (&argc, argv, "--multi_tasking", &(sv.max_tasks)) ;
  getCmdLineOption (&argc, argv, "--local", 0) ; /* consume default -- local arg */
  sv.submit = getCmdLineBool (&argc, argv, "--submit") ;
  getCmdLineInt (&argc, argv, "--clientId", &(sv.clientId)) ;

  sv.runDict = dictHandleCreate (1000, h) ;
  sv.runs = arrayHandleCreate (256, RC, h) ;
  sv.info = stackHandleCreate (10000, h) ;

  sv.classeDict = dictHandleCreate (256, h) ;
  dictAdd (sv.classeDict, "MICRODELETION", &sv.microdeletion) ;  /* ATTENTION: order of these words must stay synchronized with switch type below */
  dictAdd (sv.classeDict, "MINIDELETION", &sv.minideletion) ;
  dictAdd (sv.classeDict, "DELETION", &sv.deletion) ;
  
  dictAdd (sv.classeDict, "CIRCLE", &sv.circle) ; 
  dictAdd (sv.classeDict, "MICRODUPLICATION", &sv.microduplication) ;
  dictAdd (sv.classeDict, "MINIDUPLICATION", &sv.miniduplication) ;
  dictAdd (sv.classeDict, "DUPLICATION", &sv.duplication) ;

  dictAdd (sv.classeDict, "PALINDROME", &sv.palindrome) ;
  dictAdd (sv.classeDict, "INVERSION", &sv.inversion) ;

  dictAdd (sv.classeDict, "TRANSLOCATION", &sv.translocation) ;
  if (sv.clientId)
    { /* consume the communication layer arguments */
      sv.clientTask = taskClientInit (&argc, argv, h) ;
    }
  
  if (argc > 1) 
    usage ("Unknown argument", argc, argv) ;
  fprintf (stderr, "// %s start : max_threads=%d\n", timeShowNow(),sv.max_threads) ;

  if (sv.runsFileName) 
    svParseRunList (&sv) ;
  if (! sv.clientId && ! sv.nRuns)
    {
      fprintf (stderr, "Sorry --runs RunList is empty\n") ; 
      goto done ;
    }
      
  if (1)
    {
      int nn = 4 ;  /* we need many threads to handle the communications */
      int  max_threads = nn < sv.max_threads ? nn : sv.max_threads ;
      arrayReport (-2) ;
      wego_max_threads (max_threads + 8) ;
      wego_log ("// start wego_logs") ;
    }
  
  if (sv.max_tasks > 0)
    { 
      TASK_CONFIG config = TASK_LOCAL ; /* TASK_SUBMIT */;

      sv.serverSide = TRUE ;
      if (sv.submit)  config = TASK_SUBMIT ;
      sv.serverTask = taskServerInit (argv[0], config, h) ;
      if (sv.max_tasks < 1)
	sv.max_tasks = 1 ;
    }
  
  if (sv.de_uno || sv.de_duo || sv.de_paribus || sv.de_av)
    svDeAny (&sv) ;
  
 done:
  wego_flush () ;
  sleep (1) ; /* to allow closure of wego_log */
  ac_free (sv.h) ;
  
  if (1)
    { 
      int mx ;
      messAllocMaxStatus (&mx) ; 
      fprintf (stderr, "// %s done: max memory %d Mb\n", timeShowNow(), mx) ;
     }

  ac_free (h) ;
  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
