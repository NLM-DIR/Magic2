/*  File: dnamotifs.c
 *  Author: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg. 2002
 *-------------------------------------------------------------------
 * This file is part of the ACEVIEW PACKAGE
 *	Jean Thierry-Mieg (NCBI and CNRS France) mieg@ncbi.nlm.nih.gov
 *
 * Description:
 *   Align short sequences (probes/solexa) onto a large target (the genome)
 *
 * As far as I know the algo is completly novel.
 * The idea is to construct all 4^w chains or repeated words of length w.
 *   If the chains of length w are known, we find the chains of length w+1
 * by considering the previous chain and comparing one more letter 
 * on average 4 tests per base.
 *   We initialise the w=0 chains by linking every base to the next one
 * and the last base of each probe to the first base of the target.
 *
 * Exported functions:
 * HISTORY:
 *   2008_02_19: The code is a spin-off of chromorepeats developped to
 *               cope with the Solexa data
 *   
 * Created: Nov 2002 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */
#define PLN 250

#define debug 0 

#include "../wac/ac.h"
#include "../wh/acedna.h"
#include "../wh/dict.h"
#include "../wac/ac.h"
#include "../wh/aceio.h"

#define _Z 0x10
#define NMAX 200000000
#define NPMAX 800000
#define DIM3 1

static const char *justType = 0 ;
static BOOL justStrand, justAntistrand ;
static int NN = 0 ;
typedef struct profileStructure { float a[4] ; } PF ;
typedef struct probeStruct { int nam, ln, start, stop, shift, count, ignore ; BOOL isDown ; 
  double profileScalar, antiProfileScalar ; } PP ;

typedef struct lookStruct { 
  AC_HANDLE h ;
  ACEOUT ao ; DICT *dict ; Array dna, probes ; 
  int force, nIter, np, nbp,nIgnored, ntags
    , probeLengthMin, probeLengthMax, fold
    , hexamers, correl, doubleGas, centrage, width 
    , Rhexamers, RhexamersC1, RhexamersC2, RhexamersA1, RhexamersA2, RhexamersTours 
    , pepSplit
    ;
  int bccL, bccD, bccN, bccBestN, bccNfill, *bccCodes, *bccLk ;
  unsigned char **bccDistances, **bccLinks ;
  BOOL triplet, tripletq, homopolymer, homopolymerq, barCodeCreate, bccSlow ;
  const char *cFileName, *oFileName, *title, *profileFileName
    , *motif, *g26, *f1, *f2, *barCode ;
  DICT *ignoredProbeDict ;
} LOOK ;

/***************************************************************/
/***************************************************************/

static void gnuPlot (const char *fNam,  const char *title)
{
  vTXT txt = vtxtCreate () ;
  
  if (!title)
    title = fNam ;
  vtxtPrint (txt, "gnuplot -bg white << EOF\n") ;
  vtxtPrintf (txt, "set style data linespoints\n"
	      "set title '%s'\n"
	      "set terminal postscript color\n"
	      "set terminal postscript solid\n"
	      "set terminal postscript landscape\n"
	      "set output '%s.ps'\n"
	      "plot '%s' \n"
	      "quit\n"
	      , title, fNam, fNam) ;
  vtxtPrint (txt, "EOF\n") ;
  vtxtPrintf (txt, "gv %s.ps &\n", fNam) ;
  system (vtxtPtr (txt)) ;
  vtxtDestroy (txt) ;
  return ;
} /* gnuPlot */

/***************************************************************/

static void gnuPlotATGC (const char *fNam, const char *title)
{
  vTXT txt = vtxtCreate () ;
  
  if (!title)
    title = fNam ;
  vtxtPrint (txt, "gnuplot -bg white << EOF\n") ;
  vtxtPrintf (txt, "set style data linespoints\n"
	      "set title '%s'\n"
	      "set terminal postscript color\n"
	      "set terminal postscript solid\n"
	      "set terminal postscript landscape\n"
	      "set output '%s.ps'\n" 
              "plot '%s' using 1:3 title 'T', '%s' using 1:2 title 'A' , '%s' using 1:5 title 'C', '%s' using 1:4 title 'G'\n"
	      "quit\n"
	      , title, fNam, fNam, fNam, fNam, fNam) ;
  vtxtPrint (txt, "EOF\n") ;
  vtxtPrintf (txt, "gv %s.ps &\n", fNam) ;
  system (vtxtPtr (txt)) ;
  vtxtDestroy (txt) ;
  return ;
} /* gnuPlotATGC */

/***************************************************************/
/***************************************************************/
/* Order by fit to profile */
static int profileOrder (const void *a, const void *b)
{
  const PP *pa = (const PP*)a,  *pb = (const PP*)b ;
  double za, zb ;
  za = pa->profileScalar > pa->antiProfileScalar ? pa->profileScalar : pa->antiProfileScalar ;
  zb = pb->profileScalar > pb->antiProfileScalar ? pb->profileScalar : pb->antiProfileScalar ;

  return za > zb ? -1 : 1 ; /* big first */
} /* profileOrder */

/*************************************************************************************/

static void showCentrage (LOOK *look, const char *title)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii = 0 ;
  int nvalid = 0, xmin = 0, xmax = 0 ;
  PP *pp ;
  Array histo = arrayHandleCreate (2000, int, h) ;
  ACEOUT ao =  aceOutCreateToFile (look->cFileName, "w", h) ;
  
  if (!title) title = look->title ;
  if (!title) title = "Centrage" ;
  aceOutf (ao, "# %s\n", title) ;
  aceOutf (ao, "#Position\tNumber") ;
  for (ii = 0, pp = arrp (look->probes, 0, PP) ; ii < arrayMax (look->probes) ; ii++, pp++)
    {
      if (pp->ignore)
	continue ;
      array (histo, 1000 + pp->shift, int)++ ;
      if (!nvalid || pp->shift < xmin) xmin = pp->shift ;
      if (!nvalid || pp->shift > xmax) xmax = pp->shift ;
    nvalid++ ;
    }
  for (ii = xmin ; ii <= xmax ; ii++)
    aceOutf (ao, "\n%d\t%d", ii, array (histo, 1000 + ii, int));
	
  aceOutf (ao, "\n") ;
  ac_free (h) ;
  gnuPlot (look->cFileName, title) ;
  return ;
} /* showCentrage */

/*************************************************************************************/

static void showProfile (LOOK *look, Array profile, const char *title)
{
  PF *pf ;
  int ii = 0 ;
  int i, dim = 4 ;
  ACEOUT ao ;

  ao = aceOutCreateToFile (look->profileFileName, "w", 0) ;
  if (!title) title = look->title ;
  if (!title) title = look->oFileName ;
  if (!title) title = "Profile" ;
  aceOutf (ao, "%s\n", title ? title : "Profile") ;
  aceOutf (ao, "#Position\tA\tT\tG\tC") ;
  for (ii = 0, pf = arrp (profile, 0, PF) ; ii < arrayMax (profile) ; ii++, pf++)
    {
      aceOutf (ao, "\n%d", ii - PLN/2) ;
      for (i = 0 ; i < dim ; i++)
	aceOutf (ao, "\t%g", pf->a[i]) ;
    }
  aceOutf (ao, "\n") ; 
  ac_free (ao) ;
  gnuPlotATGC (look->profileFileName, title) ;

  return ;
} /* showProfile */

/*************************************************************************************/

static Array doGetProfile (LOOK *look, BOOL isMotif, BOOL isSymmetric, int *nnp)
{
  PP *pp ;
  PF *pf ;
  int ntot= 0, ii = 0, nused = 0, nvalid = 0 ; ;
  int i, dim = 4, x, pos, xc ;
  float z, *a ;
  Array aa = arrayHandleCreate (256, PF, look->h) ;
  Array dna = look->dna ;
  char cc ;

  pf = arrayp (aa, PLN - 1, PF) ; /* make room */
  /* collect the sequences */
  for (ii = 0, pp = arrp (look->probes, 0, PP) ; ii < arrayMax (look->probes) ; ii++, pp++)
    {
      if (
	  (!isMotif && pp->ignore) ||
	  (isMotif && pp->ignore != 2)
	  )
	continue ;
      nvalid++ ;
    }
  if (nvalid > 1) arraySort (look->probes, profileOrder) ;
  for (ii = 0, pp = arrp (look->probes, 0, PP) ; ii < arrayMax (look->probes) ; ii++, pp++)
    {
      if (
	  (!isMotif && pp->ignore) ||
	  (isMotif && pp->ignore != 2)
	  )
	continue ;
      nused++ ;
      if (0*nused > nvalid) break ;
      for (pos = -PLN/2 ; pos <= PLN/2 ; pos++)
	{
	  pf = arrp (aa, pos + PLN/2, PF) ;
	  a = pf->a ;
	  xc = pp->start + pp->ln/2 - pp->shift ;
	  x = xc + pos * (pp->isDown ? 1 : -1) ;
	  if (x >= pp->start && x < pp->stop)
	    {
	      cc = arr (dna, x, char) ;
	      if (!pp->isDown) 
		cc = complementBase(cc) ;
	      switch (cc)
		{
		case A_: i = 0 ; break ;
		case T_: i = 1 ; break ;
		case G_: i = 2 ; break ;
		case C_: i = 3 ; break ;
		default: i = -1 ; break ;
		}
	      if (i >= 0) 
		{
		  a[i]++ ; ntot++ ;
		  if (isSymmetric)
		    {
		      pf = arrp (aa, -pos + PLN/2, PF) ;
		      a = pf->a ;
		      a[i ^ 0x1]++ ;
		    }		    
		}
	    }
	}
    }
  if (nnp) *nnp = ntot ;
  fprintf (stderr, "%d bp included in this profile, %d per position\n", ntot, ntot/(PLN+1)) ;
  /* normalize */
  for (ii = 0, pf = arrp (aa, 0, PF) ; ii < arrayMax (aa) ; ii++, pf++)
    {
      a = pf->a ;
      for (i = 0, z = 0 ; i < dim ; i++)
	if (a[i] != 0) z += a[i]*a[i] ;
      z = z > 0 ? sqrt (z) : 1 ; 
      for (i = 0 ; i < dim ; i++)
	if (a[i] != 0) a[i] /= z ;
    }
  return aa ;
} /* doGetProfile */

/*************/

static Array getSymmetricProfile (LOOK *look, int *nnp)
{
  return doGetProfile (look, FALSE, TRUE, nnp) ;
} /* getProfile */

/*************/

static Array getProfile (LOOK *look, int *nnp)
{
  return doGetProfile (look, FALSE, FALSE, nnp) ;
} /* getProfile */

/*************/

static Array getMotifProfile (LOOK *look, int *nnp)
{
  return doGetProfile (look, TRUE, FALSE, nnp) ;
} /* getProfile */

/*************************************************************************************/
#ifdef JUNK
static void squishProfile (Array aa)
{
  PF *pf, *pf1, *pf2 ;
  int pos, ii ;
  int i, dim = 4 ;
  float z, zp1, zp2, zm1, zm2, *a, *a1, *a2 ;

  for (pos = 1 ; pos <= PLN/2 ; pos++)
    {
      pf1 = arrp (aa, pos + PLN/2, PF) ;
      pf2 = arrp (aa, -pos + PLN/2, PF) ;
      a1 = pf1->a ;
      a2 = pf2->a ;

      zp1 = a1[0] + a2[1] ; /* A + t */
      zm1 = a1[0] + a2[0] ; /* A + a */
      zp2 = a1[1] + a2[0] ; /* T + a */
      zm2 = a1[1] + a2[1] ; /* T + t */
      a1[0]=zp1 ; a1[1] = zp2 ;
      a1[0]=zm1 ; a1[1] = zm2 ;
      a2[0]=zm1 ; a2[1] = zm2 ;

      zp1 = a1[2] + a2[3] ; /* A + t */
      zm1 = a1[2] + a2[2] ; /* A + a */
      zp2 = a1[3] + a2[2] ; /* T + a */
      zm2 = a1[3] + a2[3] ; /* T + t */
      a1[2]=zp1 ; a1[3] = zp2 ;
      a1[2]=zm1 ; a1[3] = zm2 ;
      a2[2]=zm1 ; a2[3] = zm2 ;

    }

  /* normalize */
  for (ii = 0, pf = arrp (aa, 0, PF) ; ii < arrayMax (aa) ; ii++, pf++)
    {
      a = pf->a ;
      for (i = 0, z = 0 ; i < dim ; i++)
	if (a[i] != 0) z += a[i]*a[i] ;
      z = z > 0 ? sqrt (z) : 1 ; 
      for (i = 0 ; i < dim ; i++)
	if (a[i] != 0) a[i] /= z ;
    }
  return ;
} /* squishProfile */
#endif
/*************************************************************************************/

static double probeProbeScalar (char *ddna, 
				PP *pp1, int shift1, BOOL isDown1, 
				PP *pp2, int shift2, BOOL isDown2)
{
  int i, x1, x2, xc1, xc2, nn = 0 ;
  double z = 0 ;
  char cc1, cc2 ;
  
  xc1 = pp1->start + pp1->ln/2 - shift1 ;
  xc2 = pp2->start + pp2->ln/2 - shift2 ;
  for (i = -PLN/2 ; i <= PLN/2 ; i++)
    {
      x1 = xc1 + (isDown1 ? i : -i) ;
      x2 = xc2 + (isDown1 ? i : -i) ;
      if (x1 >= pp1->start && x1 < pp1->stop &&
	  x2 >= pp2->start && x2 < pp2->stop)
	{
	  cc1 = ddna[x1] ;
	  cc2 = ddna[x2] ;
	  if (isDown1) cc1 = complementBase(cc1) ;
	  if (isDown2) cc2 = complementBase(cc2) ;
	  if (cc1 == cc2) z++ ;
	  nn++ ;
	}
    }
  return nn ? z/nn : 0 ;
} /* probeProbeScalar */

/*************************************************************************************/

static double probeProfileScalar (char *ddna, PP *pp, int shift, BOOL probeDown, Array profile, BOOL profileDown)
{
  char cc ;
  int i, j, x, xc ;
  double z = 0 ;
  PF *pf ;
  
  xc = pp->start + pp->ln/2 - shift ;
  for (i = -PLN/2 ; i <= PLN/2 ; i++)
    {
      if (probeDown == profileDown)
	{
	  x = xc  + i ;
	  if (x >= pp->start && x < pp->stop)
	    cc = ddna[x] ;
	  else
	    continue ;
	}
      else
	{
	  x = xc - i ;
	  if (x >= pp->start && x < pp->stop)
	    {
	      cc = ddna[x] ;
	      cc = complementBase(cc) ;
	    }
	  else
	    continue ;
	}
      pf = arrp (profile, i + PLN/2, PF) ;
      switch (cc)
	{
	case A_: j = 0 ; break ;
	case T_: j = 1 ; break ;
	case G_: j = 2 ; break ;
	case C_: j = 3 ; break ;
	default: j = -1 ; break ;
	}
      if (j >= 0) 
	z += pf->a[j] ;
    }
  return z ;
} /* probeProfileScalar */

/*************************************************************************************/

static void allProbesProfileScalar (LOOK *look, Array profile) 
{
  int ii ;
  PP *pp ;
  char *ddna = arrp (look->dna, 0, char) ; 

  for (ii = 0, pp = arrp (look->probes, 0, PP) ; ii < arrayMax (look->probes) ; ii++, pp++)
    {
      pp->profileScalar = probeProfileScalar (ddna, pp, pp->shift, TRUE, profile, TRUE) ;
      pp->antiProfileScalar = probeProfileScalar (ddna, pp,  pp->shift, TRUE, profile, FALSE) ;
    }
} /* allProbesProfileScalar */

/*************************************************************************************/

static void orientProbes (LOOK *look, Array profile)
{
  int ntop = 0, nbottom = 0 ;
  int ii ;
  PP *pp ;

  for (ii = 0, pp = arrp (look->probes, 0, PP) ; ii < arrayMax (look->probes) ; ii++, pp++)
    {
      if (pp->ignore) continue ;
      if (pp->profileScalar >= pp->antiProfileScalar)
	{
	  ntop++ ;
	  pp->isDown = TRUE ;
	}
      else
	{
	  nbottom ++ ;
	  pp->isDown = FALSE ; 
	}
    }
  fprintf (stderr, "Found %d top peaks, %d bottom peaks\n ", ntop, nbottom) ;
} /* orientProbes */

/*************************************************************************************/

static int shiftProbes (LOOK *look, Array profile)
{
  PP *pp ;
  int ii = 0, s2 = 0, shift, shiftBest = 0 ;
  float z, zbest ;
  char *ddna = arrp (look->dna, 0, char) ;

  for (ii = 0, pp = arrp (look->probes, 0, PP) ; ii < arrayMax (look->probes) ; ii++, pp++)
    {
      if (pp->ignore) continue ;
      /* compute the saclar product profile shifted sequence and optimize the shift */
      for (shift = -100, zbest = -99999 ; shift < 100 ; shift++)
	{
	  z = probeProfileScalar (ddna, pp, shift, pp->isDown, profile, TRUE) ;
	  if (z > zbest || 
	      (z == zbest && s2 > shift*shift) 
	      )
	    { zbest = z ; shiftBest = shift ; s2 = shift * shift ; }
	}
      if (zbest > 0)
	{
	  if (debug) fprintf (stderr, "shifting probe %d by %d\n", ii, pp->shift+shiftBest) ;
	  pp->shift = shiftBest ;
	  if (pp->isDown) 
	    pp->profileScalar  = zbest ;
	  else
	    pp->antiProfileScalar  = zbest ;
	}
    }
  return 0 ;
} /* shiftProbes */

/*************************************************************************************/

static void showRelativeOrientation (LOOK *look, Array profile, int nn)
{
  int ii1, ii2 ;
  int  zp, zm ;
  PP *pp1, *pp2 ;
  char *ddna = arrp (look->dna, 0, char) ;

  arraySort (look->probes, profileOrder) ;
  for (ii1 = 0, pp1 = arrp (look->probes, 0, PP) ; ii1 < nn && ii1 < arrayMax (look->probes) ; ii1++, pp1++)
    {
      for (ii2 = 0, pp2 = arrp (look->probes, 0, PP) ; ii2 < nn && ii2 < arrayMax (look->probes) ; ii2++, pp2++)
	{
	  zp =  (PLN + 1) * probeProbeScalar (ddna, pp1, pp1->shift, TRUE, pp2, pp2->shift, TRUE) ; 
	  printf ("\t+%d", zp) ;
	}
      printf ("\n") ;
    }
  printf ("\n") ;
  for (ii1 = 0, pp1 = arrp (look->probes, 0, PP) ; ii1 < nn && ii1 < arrayMax (look->probes) ; ii1++, pp1++)
    {
      for (ii2 = 0, pp2 = arrp (look->probes, 0, PP) ; ii2 < nn && ii2 < arrayMax (look->probes) ; ii2++, pp2++)
	{
	  zm =  (PLN + 1) * probeProbeScalar (ddna, pp1, pp1->shift, TRUE, pp2, pp2->shift, FALSE) ;
	  printf ("\t-%d", zm) ;
	} 
      printf ("\n") ;
    }
  printf ("\n") ;
  for (ii1 = 0, pp1 = arrp (look->probes, 0, PP) ; ii1 < nn && ii1 < arrayMax (look->probes) ; ii1++, pp1++)
    printf("%s\t%d\t%s\n"
	   , pp1->isDown ? "+" : "-"
	   , pp1->shift
	   , dictName (look->dict, pp1->nam)
	   ) ;
  printf ("\n") ;
  return ;
} /* showRelativeOrientation */

/*************************************************************************************/

static void getRelativeOrientation (LOOK *look, Array profile, int nn)
{
  int ii1, ii2 ;
  int  zp, zm ;
  PP *pp1, *pp2 ;
  char *ddna = arrp (look->dna, 0, char) ;
  
  arraySort (look->probes, profileOrder) ;
  /* orient  nn probes relative to the top probe */
  for (ii1 = 0, pp1 = arrp (look->probes, 0, PP) ; ii1 <1 && ii1 < arrayMax (look->probes) ; ii1++, pp1++)
    {
      for (ii2 = 0, pp2 = arrp (look->probes, 0, PP) ; ii2 < nn && ii2 < arrayMax (look->probes) ; ii2++, pp2++)
	{
	  zp =  (PLN + 1) * probeProbeScalar (ddna, pp1, pp1->shift, TRUE, pp2, pp2->shift, TRUE) ; 
	  zm =  (PLN + 1) * probeProbeScalar (ddna, pp1, pp1->shift, TRUE, pp2, pp2->shift, FALSE) ;
	  pp2->isDown = zp > zm ? TRUE : FALSE ;
	}
    }
  /* construct a profile for the nn probes */
  ii1 = arrayMax (look->probes) ;
  arrayMax (look->probes)  = nn < ii1 ? nn : ii1 ;
  profile = getProfile (look, 0) ;
  arrayMax (look->probes) = ii1 ;

  /* orient all probes relative to this concensus */
  allProbesProfileScalar (look, profile) ; 
  orientProbes (look, profile) ;

  return ;
} /* getRelativeOrientation */

/*************************************************************************************/
/* restrict to exact central motif */
static void  restrictToMotif (LOOK *look, int seuil, int *nnp)
{
  PP *pp ;
  int ii, nn ;
  Array  profile = getMotifProfile (look, nnp) ;

  allProbesProfileScalar (look, profile) ;
  for (ii = 0, nn = 0, pp = arrp (look->probes, 0, PP) ; ii < arrayMax (look->probes) ; ii++, pp++)
    {
      if (
	  pp->ignore ||
	  (pp->isDown && pp->profileScalar < seuil) || 
	  (!pp->isDown && pp->antiProfileScalar < seuil)
	  )
	pp->ignore = 1 ;
      else 
	nn++ ;
    }
  if (nnp) *nnp = nn ;
  arrayDestroy (profile) ;
  return ;
} /* restrictToMotif */

/*************************************************************************************/

static void getMotifOrientation (LOOK *look, Array profile)
{
  int ii ;
  PP *pp ;
  
  for (ii = 0, pp = arrp (look->probes, 0, PP) ; ii < arrayMax (look->probes) ; ii++, pp++)
    pp->isDown = TRUE ;
  shiftProbes (look, profile) ; /* will compute pp->profileScalar */
  for (ii = 0, pp = arrp (look->probes, 0, PP) ; ii < arrayMax (look->probes) ; ii++, pp++)
    pp->isDown = FALSE ;
  shiftProbes (look, profile) ;  /* will compute pp->antiProfileScalar */

  orientProbes (look, profile) ;
  return ;
}

/*************************************************************************************/

static int dnaMotif (LOOK *look)
{
  int nn, nIter = look->nIter ;
  /* global statistics */
  Array profile = getMotifProfile (look, &nn) ;
  getMotifOrientation (look, profile) ;

  while (nIter--)
    {
      shiftProbes (look, profile) ;
      profile = getProfile (look, &nn) ;
    }
  if (0)
    {
      int ii ;
      for (ii = 0 ; ii < 2 ; ii++)
	{
	  allProbesProfileScalar (look, profile) ;
	  getRelativeOrientation (look, profile, 10) ;
	  nIter = 3 ;
	  while (nIter--)
	    {
	      shiftProbes (look, profile) ;
	      profile = getProfile (look, &nn) ;
	    }
	}
      showRelativeOrientation (look, profile, 20) ; 
    }
  if (look->force)
    restrictToMotif (look,look->force, &nn) ;
  printf ("Found %d sequences with the motif %s with %d exact letters\n", nn, look->motif, look->force);
  if (look->profileFileName && strcasecmp (look->profileFileName,"null"))
    {
      if (0)
	profile = getSymmetricProfile (look, 0) ;
      else
	profile = getProfile (look, 0) ;
      showProfile (look, profile, 0) ;
      if (look->cFileName) 
	showCentrage (look, "Centrage") ;
    }
  return 0 ;
} /* dnaMotif */

/*************************************************************************************/
/*************************************************************************************/

static int dnaMotifCorrel (LOOK *look)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0, np = 0, ii, ln, x, dx, bestdx, max = 0, dy, nsite1= 0, nsite2 = 0, ok = 0 ;
  int center, width = look->width ? look->width : 100 ;
  unsigned char *cp, **m1, *mm[256] ;
  Array correl = arrayHandleCreate (1001, int, h) ; 
  PP *pp ;
  unsigned char *ddna = arrp (look->dna, 0, unsigned char) ; 
  unsigned char *motif = (unsigned char *) strnew (look->motif, h) ;
  
  mm[0] = 0 ;
  for (cp = motif, ii = 0 ; ii < 254 && *cp ; cp++ )
    switch (*cp)
      {
      case '/': *cp = 0 ; ii++ ; break ;
      case '_': *cp = 'n' ; /* fall thru */
      default: 
	if (!mm[ii])  /* register a new motif */
	  {mm[ii] = cp ; mm[ii+1] = 0 ; }
	*cp = dnaEncodeChar[(int)*cp] ;
	break ;
      }
  
  for (ii = 0, pp = arrp (look->probes, 0, PP) ; ii < arrayMax (look->probes) ; ii++, pp++)
    {
      if (pp->ignore) continue ;
      x = 0 ; ok = 0 ;
      ln = pp->ln - x ;
      while (ln > 0)
	{
	  bestdx = dy = 0 ;
	  cp = ddna + pp->start + x ;
	  for (m1 = mm ; *m1 ; m1++)
	    {
	      dx = dnaPickMatch (cp, ln, *m1, 0, 0) ;
	      if (dx && (bestdx == 0 || bestdx > dx))
		{ bestdx = dx ; dy = strlen((char *)*m1)/2 ; }
	    }
	  if (!bestdx)
	    break ;
	  
	  center = x + bestdx - pp->ln/2  + dy ;
	  if (center >= -100 && center <= 100)
	    { nn++ ; if (!(ok & 0x1)) nsite1++ ; ok |= 1 ;  }
	  if (look->correl && x > 0 && bestdx <= width && 
	      (
	       (center >= -100 && center <= 100) ||
	       (center - dx >= -100 && center -dx <= 100)
	       )
	      )
	    { 
	      np++ ; array (correl, bestdx, int)++ ; 
	      if (bestdx > max) max = bestdx ; 
	      if (!(ok & 0x2))nsite2++;
	      ok |= 2;
	    }	
	  else if (look->centrage && center >= -width && center <= +width)
	    { np++ ; array (correl, width + center, int)++ ; } ;
	  x += bestdx ;
	  ln -= bestdx ;
	}
    }
  if (look->correl)
    {
      printf ("# Distance histogram between nearest \"%s\" neighbours: found\t%d\tmotifs in central 200 bp\t%d\tsites,\t%d\tpairs of motifs where one is in central 200bp in \t%d\tsites\n", look->motif, nn, nsite1, np, nsite2) ;
      for (ii = 1 ; ii < max + 1 ; ii++)
	printf ("%d\t%d\n", ii, arr (correl, ii, int)) ;
    }
  else
    {
      printf ("# Centering of %s motif\tfound\t%d\tmotifs in central 200 bp of\t%d\tsites\n", look->motif, nn, nsite1) ;
      for (ii = -width ; ii <= width ; ii++)
	printf ("%d\t%d\n", ii, arr (correl, width + ii, int)) ;
    }  
  ac_free (h) ;
  return nn ;
}  /* dnaMotifCorrel */

/*************************************************************************************/
/*************************************************************************************/


typedef struct tpStruct { int cc, n ; } TP ;

static int tpOrder (const void *va, const void *vb)
{
  const TP *a = (const TP *)va, *b = (const TP *)vb ;
  return b->n - a->n ; /* big first */
} /* tpOrder */

/*************************************************************************************/

static  void dnaTripletDistrib (LOOK *sx)
{
  TP *tp ;
  int NN = 100, i, j, k, iMax = 0, i1, i2, i3, total ; /* max length studied */
  int n = 2 ;  /* fastqfile has the sequence on line 2 modulo 4 */
  unsigned int cc ;
  Array aa = 0 ;
  int bb[65*NN], n1[4], n2[16], n3[64], s1[4], s2[4] ;
  char *cp ;
  char *val[64], buf[4], let[4] = {'A','T','G','C'} ;
  ACEIN ai = aceInCreateFromStdin (0, 0, sx->h) ;
  BOOL ok = TRUE ;

  aa = arrayCreate (256, TP) ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	{
	  cc = (i << 4) | (j << 2) | k ;
	  buf[0] = let[i] ;
	  buf[1] = let[j] ;
	  buf[2] = let[k] ;
	  buf[3] = 0 ;
	  val[cc] = strnew(buf, 0) ;
	}
  memset (n1, 0, sizeof(n1)) ;
  memset (n2, 0, sizeof(n2)) ;
  memset (s1, 0, sizeof(s1)) ;
  memset (s2, 0, sizeof(s2)) ;
  memset (n3, 0, sizeof(n3)) ;
  memset (bb, 0, sizeof(bb)) ;

  aceInSpecial (ai, "\n") ;
  i = i1 = i2 = i3 = 0 ;
  while (n++, aceInCard (ai))
    {
      cp = aceInWord(ai) ;
      if (!cp || !*cp)
	continue ;
      if (sx->tripletq)
	{
	  if (n%4==3 && sx->ignoredProbeDict && dictFind (sx->ignoredProbeDict, cp+1, 0))
	    ok = FALSE ;
	  if (n%4) continue ;
	  else { cc = 0 ; i = 0 ; }
	}
      if (sx->triplet && *cp == '>') 
	{ 
	  if (sx->ignoredProbeDict && dictFind (sx->ignoredProbeDict, cp+1, 0))
	    ok = FALSE ;
	  i = 0 ; cc = 0 ; 
	  i1 = i2 = i3 = 0 ;
	  continue ; 
	}
      for (; ok && *cp ; i++, cp++)
	{
	  i3 = i2 ; i2 = i1 ;
	  switch (*cp)
	    {
	    case 'A': case 'a': 
	      cc <<= 2 ; cc |= 0x0 ; n1[0]++ ; i1 = 0 ;
	      break ;
	    case 'T': case 't': 
	      cc <<= 2 ; cc |= 0x1 ; n1[1]++ ; i1 = 1 ;
	      break ;
	    case 'G': case 'g': 
	      cc <<= 2 ; cc |= 0x2 ; n1[2]++ ; i1 = 2 ;
	      break ;
	    case 'C': case 'c': 
	      cc <<= 2 ; cc |= 0x3 ; n1[3]++ ; i1 = 3 ; 
	      break ;
	    default: i = 2*NN ; break ;
	    }
	  if (i>0) n2[4*i2+i1]++ ;	  
	  if (i>1) n3[16*i3+4*i2+i1]++ ;
	  cc &= 0x3f ;
	  if (i > iMax && i < NN) iMax = i ;
	  if (i < NN)
	    bb[64*i + cc]++ ;
	}
      ok = TRUE ;
    }

  /* report */ 
  aceOutf (sx->ao, "\nSingle letter Counts") ;
  for (total = i = 0 ; i < 4 ; i++)
    total += n1[i] ;
  for (i = 0 ; i < 4 ; i++)
    aceOutf (sx->ao, "\n%c\t%d\t%.2f", let[i], n1[i], 100.0 * n1[i]/total) ;
 aceOutf (sx->ao, "\nTotal\t%d\t100\n", total) ;

  aceOutf (sx->ao, "\nDoublets: Read the doublet AG in column A:line G\n") ;
  aceOutf (sx->ao, "\tA\tT\tG\tC\tAny\tA\tT\tG\tC\tAny") ;
  for (total = i = 0 ; i < 4 ; i++)
    total += n1[i] ;

  for (i1 = 0 ; i1 < 4 ; i1++)
    {
      aceOutf (sx->ao, "\n%c", let[i1]) ;
      for (i2 = 0 ; i2 < 4 ; i2++)
	{
	  n = n2[4*i2+i1] ;
	  aceOutf (sx->ao, "\t%d", n) ;
	  s1[i1] += n ;
	  s2[i2] += n ;
	} 
      aceOutf (sx->ao, "\t%d", s1[i1]) ;
      for (i2 = 0 ; i2 < 4 ; i2++)
	{
	  aceOutf (sx->ao, "\t%.2f", 100.0*n2[4*i2+i1]/s1[i1]) ;
	}
      aceOutf (sx->ao, "\t%.2f", 100.0*s1[i1]/total) ;
    }
  aceOutf (sx->ao, "\nTotal") ;
  for (n = i2 = 0 ; i2 < 4 ; i2++)
    {
      aceOutf (sx->ao, "\t%d", s2[i2]) ;
      n += s2[i2] ;
    }
  aceOutf (sx->ao, "\t%d", n) ;
  for (i2 = 0 ; i2 < 4 ; i2++)
    aceOutf (sx->ao, "\t%.2f",100.0* s2[i2]/n) ;
  aceOutf (sx->ao, "\t100\n", n, total) ;
  

  aceOutf (sx->ao, "\n\nTriplet\tCounts\t%%") ;
  for (total = i3 = 0 ; i3 < 4 ; i3++)
    for (i2 = 0 ; i2 < 4 ; i2++)
      for (i1 = 0 ; i1 < 4 ; i1++)
	{
	  n =  n3[16*i3+4*i2+i1] ;
	  total += n ;
	}
  for (i3 = 0 ; i3 < 4 ; i3++)
    for (i2 = 0 ; i2 < 4 ; i2++)
      for (i1 = 0 ; i1 < 4 ; i1++)
	{
	  aceOutf (sx->ao, "\n%c%c%c", let[i3], let[i2], let[i1]) ;
	  n =  n3[16*i3+4*i2+i1] ;
	  aceOutf (sx->ao, "\t%d\t%.2f", n, 100.0*n/total) ;
	}
  aceOutf (sx->ao, "\nAny\t%d\t100\n", total) ;

  aceOutf (sx->ao, "\nTriplet Counts: Read the triplet AGC in column A:line GC\n") ;
  aceOutf (sx->ao, "Triplet") ;
  for (i3 = 0 ; i3 < 4 ; i3++)
    aceOutf (sx->ao, "\t%c", let[i3]) ;
  for (i2 = 0 ; i2 < 4 ; i2++)
    for (i1 = 0 ; i1 < 4 ; i1++)
      {
	aceOutf (sx->ao, "\n%c%c", let[i2], let[i1]) ;
	for (i3 = 0 ; i3 < 4 ; i3++)
	  {
	    n =  n3[16*i3+4*i2+i1] ;
	    aceOutf (sx->ao, "\t%d", n3[16*i3+4*i2+i1]) ;
	  }
      }


  aceOutf (sx->ao, "\n\nPositioned triplets") ;
  for (i = 2 ; i <= iMax ; i++)
    {
      aceOutf (sx->ao, "\n%d", i-1) ;
      arrayReCreate (aa, 64, TP) ;
      for (cc = 0 ; cc < 64 ; cc++)
	{
	  tp = arrayp (aa, cc, TP) ;
	  tp->cc = cc ;
	  tp->n = bb[64*i + cc] ;
	}
      arraySort (aa, tpOrder) ;
      for (j = 0 ; j < 12 ; j++)
	{
	  tp = arrayp (aa, j, TP) ;
	  if (tp->n)
	    aceOutf (sx->ao, "\t%d:%s", tp->n, val[tp->cc]) ;
	}
    }
  aceOutf (sx->ao, "\nMost frequent snakes") ;
  /* find the most frequent first three letters */
  arrayReCreate (aa, 64, TP) ;
  i = 2 ;
  for (cc = 0 ; cc < 64 ; cc++)
    {
      tp = arrayp (aa, cc, TP) ;
      tp->cc = cc ;
      tp->n = bb[64*i + cc] ;
    }
  arraySort (aa, tpOrder) ;
  /* report these and follow the highest frequency snake */
  for (i1 = 0 ; i1 < 8 ; i1++)
    {
      tp = arrayp (aa, i1, TP) ;
      aceOutf (sx->ao, "\n%s", val[tp->cc]) ;
      cc = tp->cc ;
      for (i2 = 3 ; i2 <= iMax ; i2++)
	{
	  cc = (cc<<2) & 0x3c ; n = -1 ;
	  for (i3 = 0, cc-- ; i3 < 4 ; i3++)
	    {
	      cc++ ;
	      if (bb[64*i2 + cc] > n) { n = bb[64*i2 + cc] ; i = i3 ; }
	    }
	  aceOutf (sx->ao, "%c", let[i]) ;
	  cc &= 0x3c ; cc += i ;
	}
    }
  aceOutf (sx->ao, "\n\n") ;
  
} /* sxTripletDistrib */

/*************************************************************************************/
/*************************************************************************************/

static  void dnaHomopolymerDistrib (LOOK *sx)
{
  int NN = 1000, nTags = 0, i, j, k, jMax = 0, n1[5*NN], n2[5*NN] ; /* max length studied */
  int n = 2 ;  /* fastqfile has the sequence on line 2 modulo 4 */
  char *cp, *cq, cc, *ll = "ATGCN" ;
  ACEIN ai = aceInCreateFromStdin (0, 0, sx->h) ;
  BOOL ok = TRUE ;

  memset (n1, 0, sizeof(n1)) ;
  memset (n2, 0, sizeof(n2)) ;

  aceInSpecial (ai, "\n") ;
  i = 0 ;
  while (n++, aceInCard (ai))
    {
      cp = aceInWord(ai) ;
      if (!cp || !*cp)
	continue ;
      if (sx->homopolymerq)
	{
	  if (n%4==3 && sx->ignoredProbeDict && dictFind (sx->ignoredProbeDict, cp+1, 0))
	    ok = FALSE ;
	  if (n%4) continue ;
	  else i = 0 ;
	}
      else if (*cp == '>') 
	{
	  if (sx->ignoredProbeDict && dictFind (sx->ignoredProbeDict, cp+1, 0))
	    ok = FALSE ;
	  i = 0 ; 
	  continue ;
	}
      nTags ++ ;
      for (k = 0 ; k < 5 ; k++)
	{
	  for (i = j = 0, cq = cp ; *cq ; i++, cq++)
	    {
	      cc = ace_upper (*cq) ;
	      if (cc == ll[k]) j++ ; else break ;
	    }
	  if (j >= NN) j = NN - 1 ;
	  if (ok)
	    n1[5*j+k]++ ;
	  if (j > jMax) jMax = j ;
	}
      for (k = 0 ; k < 5 ; k++)
	{
	  for (i = strlen(cp), j = 0, cq = cp + strlen(cp) - 1 ; i > 0 && *cq ; i--, cq--)
	    {
	      cc = ace_upper (*cq) ;
	      if (cc == ll[k]) j++ ; else break ;
	    }
	  if (j >= NN) j = NN - 1 ;
	  if (ok)
	    n2[5*j+k]++ ;
	  if (j > jMax) jMax = j ;
	}
      ok = TRUE ;
    }

  /* report */
  aceOutf (sx->ao, "Analysed %d sequence tags\n", nTags) ;
  for (k = 0 ; k < 5 ; k++)
    aceOutf (sx->ao, "\t5' %c", ll[k]) ;
  for (k = 0 ; k < 5 ; k++)
    aceOutf (sx->ao, "\t3' %c", ll[k]) ;
  for (j = 0 ; j <= jMax ; j++)
    {
      aceOutf (sx->ao, "\n%d", j) ;
      for (k = 0 ; k < 5 ; k++)
	aceOutf (sx->ao, "\t%d", n1[5*j+k]) ;
      for (k = 0 ; k < 5 ; k++)
	aceOutf (sx->ao, "\t%d", n2[5*j+k]) ;
    }
  aceOut (sx->ao, "\n\n") ;
} /* sxHomopolymerDistrib */

/*************************************************************************************/
/*************************************************************************************/
/* count all pairs of approximate gas motifs with at least a hemi-gas each */
static int dnaDoubleGas (LOOK *look)
{
  AC_HANDLE h = ac_new_handle () ;
  PP *pp ;
  int n1, n2, n, ii ;
  int x1, x2, dx, ln , X ;
  Array aa = arrayHandleCreate (2000, int, h) ;
  char *cp, *ddna = arrp (look->dna, 0, char) ; 
   
  aa = arrayHandleCreate (2000, int, h) ;
  for (ii = 0, pp = arrp (look->probes, 0, PP) ; ii < arrayMax (look->probes) ; ii++, pp++)
    {
      if (pp->ignore) continue ;
      ln = pp->ln ;
      for (dx = 1 ; dx < 50 ; dx++)
	for (x1 = 0 ; x1 < ln - 9 - dx ; x1++)  
	  {
	    cp = ddna + pp->start + x1 ;
	    n1 = 0 ;
	    if (cp[0] == T_ && cp[1] == T_ && cp[2] == C_) n1++ ;
	    if (cp[6] == G_ && cp[7] == A_ && cp[8] == A_) n1++ ;
	    for (x2 = x1+dx ; n1 > 0 && x2 < ln - 9 ; x2+= 10000000)
	      {  
		cp = ddna + pp->start + x2 ;
		n2 = 0 ;
		if (cp[0] == T_ && cp[1] == T_ && cp[2] == C_) n2++ ;
		if (cp[6] == G_ && cp[7] == A_ && cp[8] == A_) n2++ ;
		n = n1 + n2 ;
		if (n1 > 0 && n2 > 0)
		  {
		    X = n + 10*dx ;
		    array (aa, X, int)++ ;
		  }
	      }
	  }
    }
  aceOutf(look->ao,"Number of hemi-gas") ;
  for (dx = 1 ; dx < 50 ; dx++)
    aceOutf(look->ao,"\t%d", dx) ;
  for (n = 2 ; n <= 4 ; n++)
    {
      aceOutf(look->ao,"\n%d", n) ;
      for (dx = 1 ; dx < 50 ; dx++)
	aceOutf(look->ao,"\t%d", array(aa, n + 10*dx, int)) ;
    }
  aceOutf(look->ao,"\n") ;
  ac_free (h) ;
  return 0 ;
}  /* dnaDoubleGas */

/*************************************************************************************/
/*************************************************************************************/

static int dnaFoldMotifs (LOOK *look)
{
  int x, ii, i, j, nn=0, NN = 64, N1 = 0 ;
  int *hh ;
  PP *pp ;
  char *cp, *cq, cc,  *ddna = arrp (look->dna, 0, char) ; 

  hh = (int *)malloc (NN*NN* sizeof (int)) ;
  memset (hh, 0, NN*NN* sizeof (int)) ;
  for (ii = 0, pp = arrp (look->probes, 0, PP) ; ii < arrayMax (look->probes) ; ii++, pp++)
    {
      x = pp->start ;
      nn++ ;
      if (pp->ln > N1) N1 = pp->ln ;
      for (i = 0 ; i < NN && i < pp->ln ; i++)
	{
	  cp = ddna + x + i ;
	  hh[i + NN * i]++ ;
	  for (j = i+1 ; j < NN && j < pp->ln ; j++)
	    {
	      cq = ddna + x + j ;
	      cc = complementBase((*cq)) ;
	      if (*cp == *cq)
		hh[i + NN * j]++ ;
	      if (*cp == cc)
		hh[j + NN * i]++ ;
	    }
	}
    }
  printf ("\n") ;
  if (N1 > NN) N1 = NN ;
  printf ("Rate of identity (above the diagonal) and complementary (below) among %d motifs\n   ", nn) ;
  if (!nn) nn = 1 ;
  for (i = 0 ; i < N1 ; i++)
    {
      printf ("%d", i+1) ;
      for (j = 0 ; j < i ; j++)
	{
	  nn = hh[i + NN * i] ;
	  printf ("\t%d", 133*hh[i + NN * j]/nn - 33) ;
	}
      for (; j < N1 ; j++)
	{
	  nn = hh[j + NN * j] ;
	  printf ("\t%d", 133*hh[i + NN * j]/nn - 33) ;
	}
      printf ("\n") ;
    }

  printf ("\n") ;
  return nn ;
} /* dnaFoldMotifs */

/*************************************************************************************/
/*************************************************************************************/
typedef struct hexaStruct { const char *motif ; int n ; } HX ;
/***************************************************************/
/* Order hexamers by decreasing count */
static int hxOrder (const void *a, const void *b)
{
  const HX *pa = (const HX*)a,  *pb = (const HX*)b ;

  return pb->n - pa->n ;
} /* hxOrder */

/***************************************************************/

static int dnaCountHexamers (LOOK *look)
{
  int x1, x2,  ii, i, j, NN = look->hexamers ;
  PP *pp ;
  HX *hx ;
  Array dnaBuf = arrayCreate (NN+1, char) ;
  Array hh = arrayCreate (4096, HX) ;
  KEY k = 0, mask ;
  char *cp, *ddna = arrp (look->dna, 0, char) ; 

  mask =  (1<<(2*NN)) - 1 ;
  array (hh, mask, HX).n = 0 ; /* make room */
  aceOutf (look->ao, "Count all hexamers, cumulating the 2 strands\n") ;
  /* count all hexamers */
  for (ii = 0, pp = arrp (look->probes, 0, PP) ; ii < arrayMax (look->probes) ; ii++, pp++)
    {
      x1 = pp->start + pp->ln/2 - look->width/2  ;
      x2 = pp->start + pp->ln/2 + look->width/2  ;
      for (i = x1, j = 0, k = 0, cp = ddna + i ; i < x2 ; cp++, i++)
	{
	  if (i <  pp->start || i >=  pp->start + pp->ln) continue ;
	  k <<= 2 ;
	  switch (*cp)
	    {
	    case A_: k |= 0x0 ; j++ ; break ;
	    case T_: k |= 0x3 ; j++ ; break ;
	    case G_: k |= 0x1 ; j++ ; break ;
	    case C_: k |= 0x2 ; j++ ; break ;
	    default: j = 0 ;
	    }
	  if (j < NN)
	    continue ;
	  k &= mask ;
	  hx = arrp (hh, k, HX) ;
	  (hx->n) += pp->count ; hx->motif = cp - NN + 1 ;
	}
    }
  arraySort (hh, hxOrder) ;
  array (dnaBuf, NN-1, char) = 0 ;
  for (ii = 0 ; ii < arrayMax (hh) ; ii++)
    {
      hx = arrp (hh, ii, HX) ;
      if (!hx->n) break ;
      if (hx->motif)
	{
	  memcpy (dnaBuf->base, hx->motif, NN) ;
	  dnaDecodeArray (dnaBuf) ;
	}
      else
	memcpy (dnaBuf->base, "TOTOTO", NN) ;
      aceOutf (look->ao, "%s\t%d\n", arrp (dnaBuf, 0, char), hx->n) ;
    }
  return 0 ;
} /* dnaCountHexamers */


/*************************************************************************************/
/*************************************************************************************/
/* the idea is to loop, searching for n-mers most over expressed in
 * the c1-c2 central region relative to the a1-a2 global region
 * on the forward strand
 * eliminating recursively the sequences on the first hit
 */
typedef struct rhxStruct { const char* motif ; int kk ; int nc ; int na ; int ratio, bonus ; } RHX ;
/* Order hexamers by decreasing count */
static int rhxOrder (const void *a, const void *b)
{
  const RHX *pa = (const RHX*)a,  *pb = (const RHX*)b ;

  return pb->ratio - pa->ratio + pb->bonus - pa->bonus ;
} /* hxOrder */

static Array dnaCountRHexamersCollect (LOOK *look, Array aaOld, AC_HANDLE h)
{
  int NN = look->Rhexamers ;
  int C1 = look->RhexamersC1 ;
  int C2 = look->RhexamersC2 ;
  int A1 = look->RhexamersA1 ;
  int A2 = look->RhexamersA2 ;
  PP *pp ;
  char *cp, *ddna = arrp (look->dna, 0, char) ; 
  int x1, x2, ii, i, j, kk ;
  int pass = 0 ;
  RHX *rhx, *rhxOld ;
  KEY k = 0, kOld, mask ;
  Array aa = arrayHandleCreate (1<<NN, RHX, h) ;
  KEYSET isUsed = 0 ;
  int ncMax = 0 ;

  mask =  (1<<(2*NN)) - 1 ;
  if (! aaOld)
    {
      aaOld = arrayHandleCreate (mask+1, RHX, h) ;
      for (i = 0 ; i < mask ; i++)
	{
	  rhx = arrayp (aaOld, i, RHX) ;
	  rhx->kk = i ;
	}
    }
  else
    isUsed = keySetCreate () ;
  array (aa, mask, RHX).nc = 0 ; /* make room */
  /* count all hexamers */
  for (kk = 0 ; kk < arrayMax (aaOld) ; kk++)
    {
      rhxOld = arrp (aaOld, kk, RHX) ;
      if (rhxOld->nc > ncMax)
	ncMax = rhxOld->nc  ;
    }
  for (kk = 0 ; kk < arrayMax (aaOld) ; kk++)
    {
      rhxOld = arrp (aaOld, kk, RHX) ;
      kOld = rhxOld->kk ;
      if (rhxOld->nc < ncMax/100) continue ;      
      for (ii = 0, pp = arrp (look->probes, 0, PP) ; ii < arrayMax (look->probes) ; ii++, pp++)
	{
	  if (isUsed && keySet (isUsed, ii)) continue ;
	  for (pass = 0 ; pass < 2 ; pass++)
	    {
	      if (pass == 0)
		{
		  x1 = pp->start + C1 - 1 ;
		  x2 = pp->start + C2 ;
		}
	      else
		{
		  x1 = pp->start + A1 - 1 ;
		  x2 = pp->start + A2 ;
		}
	      for (i = x1, j = 0, k = 0, cp = ddna + i ; i < x2 ; cp++, i++)
		{
		  if (i <  pp->start || i >=  pp->start + pp->ln) continue ;
		  k <<= 2 ;
		  switch (*cp)
		    {
		    case A_: k |= 0x0 ; j++ ; break ;
		    case T_: k |= 0x3 ; j++ ; break ;
		    case G_: k |= 0x1 ; j++ ; break ;
		    case C_: k |= 0x2 ; j++ ; break ;
		    default: j = 0 ;
		    }
		  if (j < NN)
		    continue ;
		  k &= mask ;
		  if (isUsed && k != kOld) continue ;
		  
		  if (isUsed) keySet (isUsed, ii) = 1 ;
		  rhx = arrayp (aa, k, RHX) ;
		  rhx->motif = cp - NN + 1 ;
		  rhx->kk = k ;
		  if (pass == 0) 
		    {
		      rhx->nc += pp->count ;
		      rhx->ratio = 10000 ;
		    }
		  else 
		    {
		      rhx->na += pp->count ;
		      rhx->ratio = 10000.0 * (rhx->nc)/((double)rhx->na) ;
		      rhx->bonus = (ncMax ? 300 * rhx->nc/ncMax : 0) ;
		    }
		}
	    }
	}
      if (! isUsed)
	break ;
    }
  keySetDestroy (isUsed) ;
  arrayDestroy (aaOld) ;
  arraySort (aa, rhxOrder) ;
  return aa ;
} /* dnaCountRHexamersCollect */

/*************************************************************************************/

static int dnaCountRHexamers (LOOK *look)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii, cumulC, cumulA ;
  int NN = look->Rhexamers ;
  Array aa = 0 ;
  Array dnaBuf = arrayCreate (NN+1, char) ;
  RHX *rhx ;
  
  array (dnaBuf, NN, char) = 0 ;
  arrayMax(dnaBuf) = NN ; /* cheat */
  for (ii = 0 ; ii < look->RhexamersTours ; ii++)
    aa = dnaCountRHexamersCollect (look, aa, h) ; /* collect the most frequent n-mers */

  /* report */
  for (ii = cumulC = cumulA = 0 ; ii < arrayMax (aa) ; ii++)
    {
      rhx = arrayp (aa, ii, RHX) ;
      if (! rhx->motif) continue ;
      memcpy (arrp (dnaBuf, 0, char), rhx->motif, NN) ;
      dnaDecodeArray (dnaBuf) ;
      cumulC += rhx->nc ;
      cumulA += rhx->na ;
      aceOutf (look->ao, "%s\t%d\t%d\t%.2f\n"
	       , arrp (dnaBuf, 0, char), rhx->nc, rhx->na, rhx->ratio/100.0
	       ) ;
    }
  aceOutf (look->ao, "No signal\t%d\t%d\n"
	   , arrayMax(look->probes) - cumulC, arrayMax(look->probes) - cumulA
	   ) ;
  aceOutf (look->ao, "Total\t%d\t%d\n"
	   , arrayMax(look->probes), arrayMax(look->probes)
	   ) ;

  
  ac_free (h) ;
  return 1 ; 
} /* dnaCountRHexamers */

/*************************************************************************************/
/*************************************************************************************/
/* export useful statistics about power laws */
static int dnaPowerTable (LOOK *look)
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
	  aceOutf (look->ao, "\t%d",(int)n[i]) ;
      }
  return 0 ;
} /* dnaPowerTable */

/*************************************************************************************/
/*************************************************************************************/

static int dnaMotifBarDistanceDistrib (LOOK *look, DICT *dict, int limit)
{
  AC_HANDLE h = ac_new_handle () ;
  int n, n1, n2, nnn = 0, nk ;
  const char *ccp, *ccq, *cp, *cq ;
  KEYSET ks = keySetHandleCreate (h) ;
  ACEOUT ao = look->ao ;

  for (n1 = 1 ; n1 <= dictMax (dict) ; n1++)
    { 
      ccp = dictName (dict, n1) ;
      nk = strlen (ccp) ;
      for (n2 = 1 ; n2 <= n1 ; n2++)
	{
	  ccq = dictName (dict, n2) ;
	  n = 0 ; 
	  if (n1 != n2)
	    for (cp = ccp, cq = ccq ; *cp && *cq ; cp++, cq++)
	      if (*cp != *cq)n++ ;
	  if (n > 0 && 3*n < 2*nk && (limit == 0 || n < limit))
	    aceOutf (ao, "\t##\t%d:%s\t%d:%s\t%d\n", n1, ccp, n2, ccq, n) ;
	  keySet (ks, n)++ ; nnn++ ;
	}
    }
  n1 = keySet (ks, 0) ;
  for (n = 0 ; n < keySetMax(ks) ; n++)
    printf("%d\t%d\n", n, keySet (ks, n)) ;
  printf ("Total\t%d =? %d * %d/2 = %d\n\n", nnn,  n1, n1+1, n1 * (n1+1)/2) ;

  ac_free (h) ;
  return 0 ;
} /* dnaMotifBarDistanceDistrib */

/*************************************************************************************/

/* analysis of most efficient bar coding */
typedef struct codeStruct { BOOL indel, w[64] ; } CODE ;
#ifdef JUNK
static int codeOrder (const void *a, const void *b)
{
  const CODE *up = (const CODE *)a, *vp = (const CODE *)b ;
  int m, n ;
  
  for (m = 0 ; m < 64 ; m++)
    { n = up->w[m] - vp->w[m] ; if (n) return n ; }
  return 0 ;
} /* codeOrder */
#endif
/*************************************************************************************/

static BOOL dnaMotifBarCodeIndelResistant (DICT *dict, int ii, int jj, int sameMax)
{
  char buf1[10] ;
  const char *cp, *ccp = dictName (dict, ii)  ;
  const char *cq, *ccq = dictName (dict, jj)  ;
  int i, k, p ;
  int same ;
  char cc, *atgc = "ATGC" ;
  /* insert a C to the left of ii 
   * we get 7 bases, delete any of them and compare to jj
   */
  if (1)
    {
      buf1[0] = 'C' ;
      strcpy (buf1+1,ccp) ;
      for (p = 0 ; p < 7 ; p++) 
	{
	  for (cp = buf1, cq = ccq, i = 0, same = sameMax ; same > 0 && i < 6 ; i++, cp++, cq++)
	    {
	      if (i == p) cp++ ;
	      if (*cp != *cq) same-- ;
	    }
	  if (same > 0)
	    return FALSE ;
	}
    }
  /* insert an A to the right of ii 
   * we get 7 bases, delete any of them and compare to jj
   */
  if (1)
    {
      strcpy (buf1,ccp) ;
      buf1[6] = 'A' ; buf1[7] = 0 ;
      
      for (p = 0 ; p < 7 ; p++) 
	{  /* i<7 because we must include the final 0 in the comparison */ 
	  for (cp = buf1, cq = ccq, i = 0, same = sameMax ; same > 0 && i < 7 ; i++, cp++, cq++)
	    {
	      if (i == p) cp++ ;
	      if (*cp != *cq) same-- ;
	    }
	  if (same > 0)
	    return FALSE ;
	}
    }
  
  /* insert a C to the left of jj
   * we get 7 bases, delete any of them and compare to ii
   */
  if (1)
    {
      buf1[0] = 'C' ;
      strcpy (buf1+1,ccq) ;
      for (p = 0 ; p < 7 ; p++) 
	{
	  for (cp = buf1, cq = ccp, i = 0, same = sameMax ; same > 0 && i < 6 ; i++, cp++, cq++)
	    {
	      if (i == p) cp++ ;
	      if (*cp != *cq) same-- ;
	    }
	  if (same > 0)
	    return FALSE ;
	}
    }
  /* insert an A to the right of jj 
   * we get 7 bases, delete any of them and compare to ii
   */
  if (1)
    {
      strcpy (buf1,ccq) ;
      buf1[6] = 'A' ; buf1[7] = 0 ;

      for (p = 0 ; p < 7 ; p++) 
	{  /* i<7 because we must include the final 0 in the comparison */ 
	  for (cp = buf1, cq = ccp, i = 0, same = sameMax ; same > 0 && i < 7 ; i++, cp++, cq++)
	    {
	      if (i == p) cp++ ;
	      if (*cp != *cq) same-- ;
	    }
	  if (same > 0)
	    return FALSE ;
	}
    }
  /* insert an A to the right of jj 
   * we get 7 bases, delete any of them and compare to ii
   */
  if (1)
    {
      strcpy (buf1,ccq) ;
      buf1[6] = 'A' ; buf1[7] = 0 ;
      
      for (p = 0 ; p < 7 ; p++) 
	{  /* i<7 because we must include the final 0 in the comparison */ 
	  for (cp = buf1, cq = ccp, i = 0, same = sameMax ; same > 0 && i < 7 ; i++, cp++, cq++)
	    {
	      if (i == p) cp++ ;
	      if (*cp != *cq) same-- ;
	    }
	  if (same > 0)
	    return FALSE ;
	}
    }
  
  /* insert an  atgc anywhere in ii and compare the first 6 letters */
  if (1)
    {
      for (p = 0 ; p < 6 ; p++) /* position of the insertion */
	for (k = 0 ; k < 4 ; k++) /* nature of the insertion */
	  {
	    for (cc = atgc[k], cp = ccp, cq = ccq, i = 0, same = sameMax ; same > 0 && i < 6 ; i++, cp++, cq++)
	      {
		if (i == p) { cp-- ; if (*cq != cc) same-- ; }
		else  if (*cp != *cq) same-- ;
	      }
	    if (same > 0)
	      return FALSE ;
	  }
    }
  /* insert an  atgc anywhere in jj and compare the first 6 letters */
  if (1)
    {
      for (p = 0 ; p < 6 ; p++) /* position of the insertion */
	for (k = 0 ; k < 4 ; k++) /* nature of the insertion */
	  {
	    for (cc = atgc[k], cp = ccp, cq = ccq, i = 0, same = sameMax ; same > 0 && i < 6 ; i++, cp++, cq++)
	      {
		if (i == p) { cq-- ; if (*cp != cc) same-- ; }
		else  if (*cp != *cq) same-- ;
	      }
	    if (same > 0)
	      return FALSE ;
	  }
    }
  return TRUE ;
} /* dnaMotifBarCodeIndelResistant */

/*************************************************************************************/
/* verify if the code is resistant to indels */
#ifdef JUNK
static int dnaMotifBarCodeIndelDo (ACEOUT ao, Array aa, DICT *dict, int *nCodes)
{
  BOOL ok ;
  BOOL debug1 = FALSE ;
  int nn[65], nnn = 0, ii, i, j ;
  int sameMax = 3 ;
  CODE *up ;

  memset (nn, 0, sizeof(nn)) ;
  for (ii = 0 ; ii < arrayMax (aa) ; ii++)
    {
      up = arrp (aa, ii, CODE) ;
      ok = TRUE ;
      if (debug1)
	{
	  fprintf (stderr, "Code %d ::", ii) ;
	  for (i = 0 ; ok && i < 64 ; i++) 
	    fprintf (stderr, " %s", dictName(dict, up->w[i])) ;
	  fprintf (stderr, "\n") ;
	}
      for (i = 0 ; ok && i < 64 ; i++)
	if (up->w[i])
	  for (j = i + 1 ; j < 64 ; j++)
	    {
	      if (! up->w[j])
		continue ;
	      ok = dnaMotifBarCodeIndelResistant (dict, up->w[i], up->w[j], sameMax) ;
	      if (!ok && debug1)
		fprintf (stderr, "Failed i=%d j=%d,  %s %s\n", i, j, dictName(dict, up->w[i]), dictName(dict, up->w[j])) ;
	      if (! ok) 
		up->w[j] = 0 ;
	    }
      for (i = j = 0 ; i < 64 ; i++)
	if (up->w[i]) j++ ;
      up->indel = j ;
      nn[j]++ ;
    }   
  aceOutf (ao, "\n\nNumber of codes where the edit distance is 4 substitutions or one indel and %d substitutions\n", sameMax - 1) ;
  aceOutf (ao, "Number ofacceptable barcodes\tNumber of codes\n") ;
  for (i = 1 ; i <= 64 ; i++)
    {
      aceOutf (ao, "%d\t%d\n", i, nn[i]) ;
      if (nn[i]) nnn = i ;
    }
  aceOutf (ao, "\n\nCodes where the edit distance is 4 substitutions or one indel and %d substitutions\nCode", sameMax - 1) ; 
  for (i = 1 ; i <= 64 ; i++)
    aceOutf (ao, "\tCodebar_%d", i) ; 
  for (ii = j = 0 ; ii < arrayMax (aa) ; ii++)
    {
      up = arrp (aa, ii, CODE) ;
      if (up->indel != nnn)
	continue ;
      aceOutf (ao, "\n%d", ++j) ;
      for (i = 0 ; i < 64 ; i++)
	if (up->w[i])
	  aceOutf (ao, "\t%s", dictName (dict, up->w[i])) ;
    } 
  aceOutf (ao, "\n\n") ;
  *nCodes = nn[nnn] ;
  return nnn ;
} /* dnaMotifBarCodeIndelDo */
#endif
/*************************************************************************************/

static int dnaMotifBarCodeIndel (ACEOUT ao, DICT *dict) 
{
  int ii, jj ;
  KEYSET ks = keySetCreate () ;
  BOOL ok ;

  for (ii = 1 ; ii <= dictMax (dict) ; ii++)
    keySet (ks, ii) = 1 ;
  for (ii = 1 ; ii <= dictMax (dict) ; ii++)
    {
      if (!keySet (ks, ii)) continue ;
      for (jj = 1 ; jj <= dictMax (dict) ; jj++)
	{
	  if (ii == jj || !keySet (ks, jj))
	    continue ;
	  ok = dnaMotifBarCodeIndelResistant (dict, ii, jj, 7) ;
	  if (! ok)
	    keySet (ks, jj) = 0 ;
	}
    }
  for (jj = 0, ii = 1 ; ii <= dictMax (dict) ; ii++)
    if (keySet (ks, ii))
      jj++ ;

  fprintf (stderr, "Constructed %d different indel resistant words\n", jj) ;

  keySetDestroy (ks) ;

  return jj ;
} /* dnaMotifBarCodeIndel */

/*************************************************************************************/
#ifdef JUNK
static int dnaMotifBarCode3 (LOOK *look)
{
  AC_HANDLE h = ac_new_handle () ;
  char *codes1[] = { "ATT", "AGG", "ACC", "TGC", "TCG" } ;
  char *cp ;
  char buf[64], *cr, *codes2[65], *atgc = "ATGC" ;
  int n, nCodes ;
  int NN = 16 ;   /* 4 for the words at distance 3,  16 for the 16 words at distance 2 */
  ACEOUT ao = aceOutCreate (look->oFileName, ".barcode", 0, look->h) ;
  char *permut[7]= { "ATGC", "ATCG", "AGTC", "AGCT", "ACGT", "ACTG", 0 } ;
  int i, j, k, w, m, i1, i2, z, z1, z2 ;
  CODE *up ;
  DICT *dict = dictHandleCreate (10000, h) ;
  Array aa = arrayHandleCreate (100000, CODE, h) ;
  Array bb = arrayHandleCreate (65, int, h) ;

  if (0) dictAdd (dict, "toto", 0) ; /* avoid zero */
  for (z = 0, nCodes = 0 ; z < 4 * 4 * 4 ; z++)
    {
      z1 = z ;
      for (i = 2 ; i>= 0 ; i--)
	{ 
	  buf[i] = atgc[z1 & 0x3] ; 
	  z1 >>= 2 ;
	}
      buf[3] = 0 ;
      dictAdd (dict, buf, 0) ;
    }

  /* construct the double symbols */
  n = 0 ;
  for (i = 0 ; i < 65 ; i++) codes2[i] = halloc (7,h) ;
  strcpy (codes2[n++], "AAA") ;
  if (0)
    {
      NN = 4 ;
      strcpy (codes2[n++], "TTT") ;
      strcpy (codes2[n++], "GGG") ;
      strcpy (codes2[n++], "CCC") ;
    }
  else
    {
      NN = 16 ;
      for (i = 0 ; i < 5 ; i++)
	{
	  cp = codes1[i] ;
	  cr = codes2[n++] ;
	  cr[0] = cp[0] ; cr[1] = cp[1] ; cr[2] = cp[2] ; cr[3] = 0 ;  
	  cr = codes2[n++] ;
	  cr[0] = cp[1] ; cr[1] = cp[0] ; cr[2] = cp[2] ; cr[3] = 0 ;  
	  cr = codes2[n++] ;
	  cr[0] = cp[1] ; cr[1] = cp[2] ; cr[2] = cp[0] ; cr[3] = 0 ;  
	}
    }
  /* construct all 24 possible permutations */
  n = 0 ;
  for (z = 0, nCodes = 0 ; z < 24 * 24 * 24 ; z++)
    {
      if (0 && z > 1) break ;
      up = arrayp (aa, n++, CODE) ; nCodes++ ;
      for (m = 0 ; m < NN ; m++) /* NN words in the code */
	{
	  cr = codes2[m] ; z1 = 24 * z ;
	  for (i1 = 0 ; i1 < 3 ; i1++) /* 3 letters */
	    {
	      z1 /= 24 ; z2 = z1 % 24 ; j = z2 % 6 ; k = z2 / 4 ;
	      cp = permut[j] ;
	      for (i2 = 0 ; i2 < 4 ; i2++) /* 4 letters in the alphabet */ 
		if (cr[i1] == atgc[i2])
		  {  buf[i1] = cp [(i2 + k) % 4] ; break ; }
	    }
	  buf[i1] = 0 ;	      
	  dictAdd (dict, buf, &w) ;
	  up->w[m] = w ; 
	}
      /* sort the code */
      for (m = 0 ; m < NN ; m++)
	array (bb, m, int) = up->w[m] ;
      arraySort (bb, intOrder) ;
      for (m = 0 ; m < NN ; m++)
	up->w[m] = array (bb, m, int) ;	 
      if (nCodes % 1000000 == 0)
	{
	  arraySort (aa, codeOrder) ; 
	  arrayCompress (aa) ;
	  n = arrayMax (aa) ;
	}
    }
 if (0)
   for (n = 0 ; n < arrayMax (aa)  ; n++)
    {
      up = arrayp (aa, n, CODE) ; nCodes++ ;
      if (up->w[0]==1)
	{
	  aceOutf (ao, "\n%d", n) ;
	  for (m = 0 ; m < NN ; m++)
	    aceOutf(ao, "\t%s", up->w[m] ? dictName(dict, up->w[m]) : "XXX") ;
	}
    }
  arraySort (aa, codeOrder) ;
  aceOutf (ao, "Constructed %d 16 words codes\n", nCodes) ;
  if (0) for (n = 0 ; n < arrayMax (aa)  ; n++)
    {
      up = arrayp (aa, n, CODE) ; nCodes++ ;
      if (up->w[0]==1)
	{
	  aceOutf (ao, "\n%d", n) ;
	  for (m = 0 ; m < NN ; m++)
	    aceOutf(ao, "\t%s", up->w[m] ? dictName(dict, up->w[m]) : "XXX") ;
	}
    }
  aceOutf(ao, "\n") ;
  arrayCompress (aa) ;
  aceOutf (ao, "Constructed %d different 16 words codes using %d words\n", arrayMax (aa), dictMax (dict)) ;
  for (n = 0 ; n < arrayMax (aa)  ; n++)
    {
      up = arrayp (aa, n, CODE) ; nCodes++ ;
      aceOutf (ao, "\n%d", n) ;
      for (m = 0 ; m < NN ; m++)
	aceOutf(ao, "\t%s", up->w[m] ? dictName(dict, up->w[m]) : "XXX") ;
    }
  aceOutf(ao, "\n") ;
  return 0 ;
}  /* dnaMotifBarCode3 */
#endif
/*************************************************************************************/
/*************************************************************************************/
/* construct the double symbols  by letter permuttation */
static int dnaMotifBarCodeConstructDoubles1 (DICT *dict1, DICT *dict2)
{
  DICT *myDict = dictCreate (1000) ; 
  int ii, i, n1, n, pass ;
  const char *ccp ;
  char *cp1, buf1[256] ;
  char *cp2, buf2[256] ;
  char *cp3, buf3[256] ;
  char *cp4, buf4[256] ;

  /*  char *codes1[] = { "ATT", "AGG", "ACC", "TGC", "TCG", "NKK", "NLL", "NMM", "KLM", "KML", "ZXX", "ZYY", "ZUU", "XYU", "XUY", "SVV", "SWW", "SRR", "VWR", "VRW"} ; */
  /* construct the double symbols  by letter permuttation */
for (pass = 0 ; pass < 2 ; pass++)
  for (ii = 1 ; ii <= dictMax (dict1) ; ii++)
    {
      cp1 = buf1 ; 
      cp2 = buf2 ; 
      cp3 = buf3 ; 
      cp4 = buf4 ; 
      for (ccp = dictName (dict1, ii) ; *ccp ; ccp++)
	switch ((int)*ccp)
	  {
	  case 'A':
	    *cp1++ = 'A' ; *cp1++ = 'A' ;
	    *cp2++ = 'C' ; *cp2++ = 'G' ;
	    *cp3++ = 'T' ; *cp3++ = 'C' ;
	    *cp4++ = 'G' ; *cp4++ = 'T' ;
	    break ;
	  case 'T': 
	    *cp1++ = 'T' ; *cp1++ = 'T' ;
	    *cp2++ = 'A' ; *cp2++ = 'T' ;
	    *cp3++ = 'A' ; *cp3++ = 'G' ;
	    *cp4++ = 'A' ; *cp4++ = 'C' ;
	    break ;
	  case 'G': 
	    *cp1++ = 'G' ; *cp1++ = 'G' ;
	    *cp2++ = 'T' ; *cp2++ = 'A' ;
	    *cp3++ = 'G' ; *cp3++ = 'A' ;
	    *cp4++ = 'C' ; *cp4++ = 'A' ;
	    break ;
	  case 'C': 
	    *cp1++ = 'C' ; *cp1++ = 'C' ;
	    *cp2++ = 'G' ; *cp2++ = 'C' ;
	    *cp3++ = 'C' ; *cp3++ = 'T' ;
	    *cp4++ = 'T' ; *cp4++ = 'G' ;
	    break ;
	  }
      *cp1 = *cp2 = *cp3 = *cp4 = 0 ;
      dictAdd (myDict, buf1, 0) ;
      dictAdd (dict2, buf1, 0) ;
      dictAdd (dict2, buf2, 0) ;
      dictAdd (dict2, buf3, 0) ;
      dictAdd (dict2, buf4, 0) ;
      if (pass == 1)
	{
	  printf ("%s", buf1) ;

	  for (n1 = 100000, i = 1 ; i <= dictMax (myDict) ; i++)
	    {
	      n = 0 ;
	      for (ccp = dictName (myDict, i), cp1 = buf2 ; *ccp ; cp1++, ccp++)
		if (*cp1 != *ccp) n++ ;
	      if (n < n1) n1 = n ;
	    }
	  printf ("\t%d%s%s", n1, n1 < 8 ? "*" : "", buf2) ;
	  for (n1 = 100000, i = 1 ; i <= dictMax (myDict) ; i++)
	    {
	      n = 0 ;
	      for (ccp = dictName (myDict, i), cp1 = buf3 ; *ccp ; cp1++, ccp++)
		if (*cp1 != *ccp) n++ ;  
	      if (n < n1) n1 = n ;
	    }
	  printf ("\t%d%s%s", n1, n1 < 8 ? "*" : "", buf3) ;

	  for (n1 = 100000, i = 1 ; i <= dictMax (myDict) ; i++)
	    {
	      n = 0 ;
	      for (ccp = dictName (myDict, i), cp1 = buf4 ; *ccp ; cp1++, ccp++)
		if (*cp1 != *ccp) n++ ;  
	      if (n < n1) n1 = n ;
	    }
	  printf ("\t%d%s%s", n1, n1 < 8 ? "*" : "", buf4) ; 
	  if (n < n1) n1 = n ;

	  printf ("\n") ;	  
	}
    }
  
  
 ac_free (myDict) ; 
 return dictMax (dict2) ;
} /* dnaMotifBarCodeConstructDoubles1 */

/*************************************************************************************/
/* construct the double symbols by shifting the last letter */
#ifdef JUNK
static int dnaMotifBarCodeConstructDoubles2 (DICT *dict1, DICT *dict2)
{
  DICT *myDict = dictCreate (1000) ; 
  int ii, i, n1, n, n3, pass, nnn = 0 ;
  const char *ccp ;
  char *cp1, buf1[256] ;
  char *cp2, buf2[256] ;
  char *cp3, buf3[256] ;
  char *cp4, buf4[256] ;
  char *cq1, duf1[256] ;
  char *cq2, duf2[256] ;
  char *cq3, duf3[256] ;
  char *cq4, duf4[256] ;

  memset (buf1, 0, sizeof (buf1)) ;
  memset (buf2, 0, sizeof (buf2)) ;
  memset (buf3, 0, sizeof (buf3)) ;
  memset (buf4, 0, sizeof (buf4)) ;

  memset (duf1, 0, sizeof (duf1)) ;
  memset (duf2, 0, sizeof (duf2)) ;
  memset (duf3, 0, sizeof (duf3)) ;
  memset (duf4, 0, sizeof (duf4)) ;

  /*  char *codes1[] = { "ATT", "AGG", "ACC", "TGC", "TCG", "NKK", "NLL", "NMM", "KLM", "KML", "ZXX", "ZYY", "ZUU", "XYU", "XUY", "SVV", "SWW", "SRR", "VWR", "VRW"} ; */
  /* we shift the last letter */
  for (pass = 0 ; pass < 2 ; pass++)
    {
      for (ii = 1 ; ii <= dictMax (dict1) ; ii++)
	{
	  ccp = dictName (dict1, ii) ; 
	  
	  strcpy (buf1, ccp) ;
	  strcpy (buf2, ccp) ;
	  strcpy (buf3, ccp) ;
	  strcpy (buf4, ccp) ;
	  
	  memset (duf1, 0, sizeof (duf1)) ;
	  memset (duf2, 0, sizeof (duf2)) ;
	  memset (duf3, 0, sizeof (duf3)) ;
	  memset (duf4, 0, sizeof (duf4)) ;
	  
	  n = strlen (buf1) ; n3 = 2*n/3 ;
	  for (i = n - n3 ; i < n ; i++)
	    {
	      cp1 = buf1 + i ;
	      cp2 = buf2 + i ;
	      cp3 = buf3 + i ;
	      cp4 = buf4 + i ;
	      
	      switch ((int)*cp1)
		{
		case 'A':
		  *cp2 = 'T' ;
		  *cp3 = 'G' ;
		  *cp4 = 'C' ;
		  break ;
		case 'T':
		  *cp2 = 'G' ;
		  *cp3 = 'C' ;
		  *cp4 = 'A' ;
		  break ;
		case 'G':
		  *cp2 = 'C' ;
		  *cp3 = 'A' ;
		  *cp4 = 'T' ;
		  break ;
		case 'C':
		  *cp2 = 'A' ;
		  *cp3 = 'T' ;
		  *cp4 = 'G' ;
		  break ;
		}
	    }
	  
	  /* now we permut the double alphabets in each buffer on the remaining letters */
	  cq1 = duf1 ; 
	  cq2 = duf2 ; 
	  cq3 = duf3 ; 
	  cq4 = duf4 ; 

	  for (i = 0 ; i < n3 ; i++)
	    {
	      ccp = dictName (dict1, ii) + i ; 

	      switch ((int)*ccp)
		{
		case 'A':
		  *cq1++ = 'A' ; *cq1++ = 'A' ;
		  *cq2++ = 'C' ; *cq2++ = 'G' ;
		  *cq3++ = 'T' ; *cq3++ = 'C' ;
		  *cq4++ = 'G' ; *cq4++ = 'T' ;
		  break ;
		case 'T': 
		  *cq1++ = 'T' ; *cq1++ = 'T' ;
		  *cq2++ = 'A' ; *cq2++ = 'T' ;
		  *cq3++ = 'A' ; *cq3++ = 'G' ;
		  *cq4++ = 'A' ; *cq4++ = 'C' ;
		  break ;
		case 'G': 
		  *cq1++ = 'G' ; *cq1++ = 'G' ;
		  *cq2++ = 'T' ; *cq2++ = 'A' ;
		  *cq3++ = 'G' ; *cq3++ = 'A' ;
		  *cq4++ = 'C' ; *cq4++ = 'A' ;
		  break ;
		case 'C': 
		  *cq1++ = 'C' ; *cq1++ = 'C' ;
		  *cq2++ = 'G' ; *cq2++ = 'C' ;
		  *cq3++ = 'C' ; *cq3++ = 'T' ;
		  *cq4++ = 'T' ; *cq4++ = 'G' ;
		  break ;
		}
	    }
	  /* now we double the last letters without changing alphabet */
	  cp1 = buf1 + n3 ; 
	  cp2 = buf2 + n3 ; 
	  cp3 = buf3 + n3 ; 
	  cp4 = buf4 + n3 ; 
	  for ( ; i < n ; cp1++, cp2++, cp3++, cp4++, i++)
	    {
	      *cq1++ = *cp1 ; *cq1++ = *cp1 ;
	      *cq2++ = *cp2 ; *cq2++ = *cp2 ;
	      *cq3++ = *cp3 ; *cq3++ = *cp3 ;
	      *cq4++ = *cp4 ; *cq4++ = *cp4 ;
	    }
	  
	  
	  *cq1 = *cq2 = *cq3 = *cq4 = 0 ;
	  nnn += 4 ;
	  dictAdd (myDict, duf1, 0) ;
	  dictAdd (dict2, duf1, 0) ;
	  dictAdd (dict2, duf2, 0) ;
	  dictAdd (dict2, duf3, 0) ;
	  dictAdd (dict2, duf4, 0) ;

	  if (pass == 1)
	    {
	      printf ("%s", duf1) ;
	      
	      for (n1 = 100000, i = 1 ; i <= dictMax (myDict) ; i++)
		{
		  n = 0 ;
		  for (ccp = dictName (myDict, i), cp1 = duf2 ; *ccp ; cp1++, ccp++)
		    if (*cp1 != *ccp) n++ ;
		  if (n < n1) n1 = n ;
		}
	      printf ("\t%d%s%s", n1, n1 < 8 ? "*" : "", duf2) ;
	      for (n1 = 100000, i = 1 ; i <= dictMax (myDict) ; i++)
		{
		  n = 0 ;
		  for (ccp = dictName (myDict, i), cp1 = duf3 ; *ccp ; cp1++, ccp++)
		    if (*cp1 != *ccp) n++ ;  
		  if (n < n1) n1 = n ;
		}
	      printf ("\t%d%s%s", n1, n1 < 8 ? "*" : "", duf3) ;
	      
	      for (n1 = 100000, i = 1 ; i <= dictMax (myDict) ; i++)
		{
		  n = 0 ;
		  for (ccp = dictName (myDict, i), cp1 = duf4 ; *ccp ; cp1++, ccp++)
		    if (*cp1 != *ccp) n++ ;  
		  if (n < n1) n1 = n ;
		}
	      printf ("\t%d%s%s", n1, n1 < 8 ? "*" : "", duf4) ; 
	      if (n < n1) n1 = n ;
	      
	      printf ("\n") ;	  
	    }
	}
    }

  printf ("Created %d words, %d distinct\n", nnn, dictMax(dict2)) ;
  return dictMax (dict2) ;
} /* dnaMotifBarCodeConstructDoubles2 */
#endif
/*************************************************************************************/

static int dnaMotifBarCodeConstruct (LOOK *look)
{
  AC_HANDLE h = ac_new_handle () ;
  char *codes1[] = { "ATT", "AGG", "ACC", "TGC", "TCG" } ;
  char *cp ;
  char buf[64] ;
  int i ;
  DICT *dict3 = dictHandleCreate (10000, h) ;
  DICT *dict6 = dictHandleCreate (10000, h) ;
  DICT *dict12 = dictHandleCreate (10000, h) ;

  /* construct the 3 letter codes */
  dictAdd (dict3, "AAA", 0) ;
  for (i = 0 ; i < 5 ; i++)
    {
      cp = codes1[i] ;
      buf[0] = cp[0] ; buf[1] = cp[1] ; buf[2] = cp[2] ; buf[3] = 0 ;  
      dictAdd (dict3, buf, 0) ;
      buf[0] = cp[1] ; buf[1] = cp[0] ; buf[2] = cp[2] ; buf[3] = 0 ;  
      dictAdd (dict3, buf, 0) ;
      buf[0] = cp[1] ; buf[1] = cp[2] ; buf[2] = cp[0] ; buf[3] = 0 ;  
      dictAdd (dict3, buf, 0) ;
    }
  dnaMotifBarDistanceDistrib (look, dict3, 2) ;
  /* double the code */
  dnaMotifBarCodeConstructDoubles1 (dict3, dict6) ; 
  if (0) return 16 ;
  dnaMotifBarDistanceDistrib (look, dict6, 4) ;


  /* double the code */
  dnaMotifBarCodeConstructDoubles1 (dict6, dict12) ; 
  return 16 ;
  dnaMotifBarDistanceDistrib (look, dict12, 6) ;
}  /* dnaMotifBarCodeConstruct */

/*************************************************************************************/
/* construction of the linear span of a set of F4 generators */
/* the 13 generators are Z2 Z2 valued with 00 01 10 11 noted 0,1,2,3 */
static void dnaMotifBarCodeBogdanova_12_5  (LOOK *look)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *gen[] = {    "111110000000"
			 , "110001110000"
			 , "101001001100"
			 , "010100101010"
			 , "211001000001"
			 , "022001101000"
			 , "021200001001"
			 , "130102100001"
			 , "110020201000"
			 , "011121021001"
			 , "101100102200"
			 , "000121002020"
			 , "131121002002"
			 , 0
  } ;
  unsigned int gen2[13] ;
  unsigned int NN, x, z, z1, zmax, ii, i ; 
  const char *ccp, **ccpp ;
  char buf[256] ;
  DICT *dict = dictHandleCreate (4000, h) ;
  
  /* expand as binary numbers of length 24 */
  for (ii = NN = 0, ccpp = gen ; *ccpp ; ii++, ccpp++)
    {
      for (x = i = 0, ccp = *ccpp ; *ccp ; i++, ccp++)
	{ 
	  x <<= 2 ;
	  x |= ((*ccp - '0') & 0x3) ;
	}
      gen2[ii] = x ; 
      NN++ ; /* number of generators */
      if (i > 256) 
	messcrash ("dnaMotifBarCodeBogdanova i=%d, Please increase buf[] declaration", i) ;
    }
  /* construct 2^13 generators by Z2 Z2 addition , i.e. addition without carry over, i.e. XOR */
  zmax = (1 << (NN+1)) - 1 ;
  for (z = 0 ; z < zmax ; z++)
    {
      for (z1 = z, x = 0, ccpp = gen, ii = 0 ; *ccpp ; ii++, z1 >>= 1, ccpp++)
	if ((z1 & 0x1)) x ^= gen2[ii] ;
      for (i = 0, z1 = x, ccp = gen[0] ; *ccp ; i++, z1 >>=2, ccp++)
	switch ((int)(z1 & 0x3))
	  {
	  case 0: buf[i] = 'A' ; break ;
	  case 1: buf[i] = 'G' ; break ;
	  case 2: buf[i] = 'C' ; break ;
	  case 3: buf[i] = 'T' ; break ;
	  }
      buf[i] = 0 ;
      dictAdd (dict, buf, 0) ;
    }
  dnaMotifBarDistanceDistrib (look, dict, 5) ;
  {
    ACEOUT ao = aceOutCreate (look->oFileName, ".barCodes", FALSE, h) ;
    for (i = 1 ; i <= dictMax (dict) ; i++)
      aceOutf (ao, "%s\n", dictName (dict, i)) ;
    ac_free (ao) ;
  } 

  dnaMotifBarCodeIndel (look->ao, dict) ;

  ac_free (h) ;
  return ;
}  /* dnaMotifBarCodeBogdanova_12_5 */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/* Bogdanova (8320,4)
 *
 * Author : NTM
 * Created : 16/11/2011
 *

 print to stdout bogdanova's 8,4,320 code:
 320 codewords, each word of length 8, min distance 4.
 Taken in Bogdanova et al, Codes Crypto, 2001.

 Alphabet is GF(4),
 represented as 0, 1, 2==alpha, 3==alpha^2 where alpha is a primitive element,
 and in the output we arbitrarily print O==A, 1==T, 2==G, 3==C.
 
*/
/*************************************************************************************/
/* f4Minus: substraction x-y in GF(4) */
static int f4Minus(int x, int y)
{
  if ((x<0)||(x>3)||(y<0)||(y>3))
    messcrash ("function f4Minus called with bad args x and y: %d, %d", x, y) ;

  else if ((x==0) || (y==0))
    return (x+y) ;
  else if (x==y)
    return 0 ;
  else if (x==1)
    return (5-y) ;
  else if (x==2)
    {
      if (y==1)
	return 3 ;
      else if (y==3)
	return 1 ;
    }
  else // x==3 and (y==1 or y==2)
    return (3-y) ;
  return 0 ;
}

/*************************************************************************************/
/* f4Mult: multiply x*y in GF(4) */
static int f4Mult(int x, int y)
{
  if ((x<0)||(x>3)||(y<0)||(y>3))
   messcrash ("function f4Mult called with bad args x and y: %d, %d", x, y) ;

  else if ((x==0)||(y==0))
    return 0 ;
  else if ((x==1)||(y==1))
    return (x*y) ;
  else if (x==y)
    /* 2*2===3 and 3*3==2 */
    return (5-x) ;
  else /* 2*3==3*2==1 */
    return 1 ;
  
  return 0 ;
} /* f4Mult */

/*************************************************************************************/
/* fillWord: take a pointer to a mem area storing 8 ints (must be malloc'd);
   first 3 ints must be filled with some values in 0..3;
   this function filles the remaining 5 ints with the values so that
   word contains a valid codeword using bogdanova's method for the 8,4,320 code.
   s_word is a pointer to a vector of 5 ints: it is the S word to use.
*/
static void f4FillWord(int *word, int *s_word)
{
  int x1 = word[0] ;
  int x2 = word[1] ;
  int x3 = word[2] ;

  // x4 == s1 - x1 - 3*x2 - x3
  word[3] = f4Minus(f4Minus(f4Minus(s_word[0],x1),f4Mult(3,x2)),x3) ;
  // x5 == s2 - x1 - 2*x2
  word[4] = f4Minus(f4Minus(s_word[1],x1),f4Mult(2,x2)) ;
  // x6 == s3 - x1 - x2
  word[5] = f4Minus(f4Minus(s_word[2],x1),x2) ;
  // x7 == s4 - x3
  word[6] = f4Minus(s_word[3],x3) ;
  // x8 == s5 - x3
  word[7] = f4Minus(s_word[4],x3) ;
}

/*************************************************************************************/
/* f4PrintWords: given an S word, print all codewords in ATGC */
static void f4PrintWords(int *s_word)
{
  int x1, x2, x3, i ;
  
  int *word = messalloc(8*sizeof(int)) ;
  if (word==NULL)
    {
      fprintf(stderr, "no memory for word\n") ;
      exit(1) ;
    }

  for (x1=0; x1 < 4; x1++)
    for (x2=0; x2 < 4; x2++)
      for (x3=0; x3 < 4; x3++)
	{
	  word[0] = x1 ;
	  word[1] = x2 ;
	  word[2] = x3 ;
	  
	  f4FillWord(word, s_word) ;

	  for (i=0; i<8; i++)
	    {
	      switch (word[i])
		{
		case 0 :
		  printf("A") ;
		  break ;
		case 1 :
		  printf("T") ;
		  break ;
		case 2 :
		  printf("G") ;
		  break ;
		case 3 :
		  printf("C") ;
		  break ;
		default :
		  fprintf(stderr, "in f4PrintWords, bad value for word[%d]: %d\n", i, word[i]) ;
		  exit(1) ;
		}
	    }
	  // done printing this word
	  printf("\n") ;
	}

 messfree(word) ;

}

/*************************************************************************************/

static void dnaMotifBarCodeBogdanova_8_4 (LOOK *look)
{
  int *s_word = messalloc(5*sizeof(int)) ;
  if (s_word == NULL)
    {
      fprintf(stderr, "no memory for s_word\n") ;
      exit(1) ;
    }

  // clumsy to initialize s_word but I don't remember how to do this cleanly in C */
  s_word[0] = 0 ; s_word[1] = 0 ; s_word[2] = 2 ; s_word[3] = 0 ; s_word[4] = 0 ; 
  f4PrintWords(s_word) ;

  s_word[0] = 0 ; s_word[1] = 3 ; s_word[2] = 3 ; s_word[3] = 2 ; s_word[4] = 3 ; 
  f4PrintWords(s_word) ;

  s_word[0] = 1 ; s_word[1] = 1 ; s_word[2] = 0 ; s_word[3] = 2 ; s_word[4] = 2 ; 
  f4PrintWords(s_word) ;

  s_word[0] = 2 ; s_word[1] = 0 ; s_word[2] = 1 ; s_word[3] = 3 ; s_word[4] = 1 ; 
  f4PrintWords(s_word) ;

  s_word[0] = 3 ; s_word[1] = 2 ; s_word[2] = 1 ; s_word[3] = 1 ; s_word[4] = 2 ; 
  f4PrintWords(s_word) ;


  // all done, clean up
  messfree(s_word) ;

  return ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
#ifdef JUNK
static int dnaMotifBarCode (LOOK *look)
{
  AC_HANDLE h = ac_new_handle () ;
  /* codes12 are as read, i.e. as the complement of the template, gatcggaagagcacacgtctgaactccagtcacTGACCAatctcgtatgccgtcttctgctt 
   * they are numbered as in the documentation of Illumina */
  /* char *codes12[] ={"ATCACG", "CGATGT", "TTAGGC", "TGACCA", "ACAGTG", "GCCAAT", "CAGATC", "ACTTGA", "GATCAG", "TAGCTT", "GGCTAC", "CTTGTA", 0} ; */
  /* these are the same codes numbered as in the libraries of Wang */
  char *libs[] = { "A_1", "A_2", "B_1", "B_2", "C_1", "C_2", "D_1", "D_2", "E_1", "E_2", "F_1", "F_2", 0 } ;
  char *codes[] = { "ATCACG", "CAGATC", "ACTTGA", "CGATGT", "GATCAG", "TTAGGC", "TAGCTT", "TGACCA", "ACAGTG", "GGCTAC", "CTTGTA", "GCCAAT"
 , 0 } ;
  char *codes1[] = { "ATT", "AGG", "ACC", "TGC", "TCG", "NKK", "NLL", "NMM", "KLM", "KML", "ZXX", "ZYY", "ZUU", "XYU", "XUY", "SVV", "SWW", "SRR", "VWR", "VRW"} ;
  char *cp, *cq, **cpp, **cqq, **lib ;
  char buf[64], *cr ;
  char *codes2[65], *atgc = "ATGC" ;
  int n, nnn, nCodes ;
  ACEOUT ao = aceOutCreate (look->oFileName, ".barcode", 0, look->h) ;
  char *permut[7]= { "ATGC", "ATCG", "AGTC", "AGCT", "ACGT", "ACTG", 0 } ;
  int a1, a2, b1, b2, nN, exact, i, j, k, w, m, ii, i1, i2, z, z1, z2 ;
  Array a = 0, b = 0, err = 0 ;
  A_ERR *ep ;
  CODE *up ;
  DICT *dict = dictHandleCreate (10000, h) ;
  Array aa = arrayHandleCreate (100000, CODE, h) ;
  Array bb = arrayHandleCreate (65, int, h) ;

  if (0) dictAdd (dict, "toto", 0) ; /* avoid zero */
  for (z = 0, nCodes = 0 ; z < 4 * 4 * 4 * 4 * 4 * 4 ; z++)
    {
      z1 = z ;
      for (i = 5 ; i>= 0 ; i--)
	{ 
	  buf[i] = atgc[z1 & 0x3] ; 
	  z1 >>= 2 ;
	}
      buf[6] = 0 ;
      dictAdd (dict, buf, 0) ;
    }
  /* construct the double symbols */
  n = 0 ;
  for (i = 0 ; i < 65 ; i++) codes2[i] = halloc (7,h) ;
  strcpy (codes2[n++], "AAA") ;
  strcpy (codes2[n++], "NNN") ;
  strcpy (codes2[n++], "ZZZ") ;
  strcpy (codes2[n++], "SSS") ;
  for (i = 0 ; i < 20 ; i++)
    {
      cp = codes1[i] ;
      cr = codes2[n++] ;
      cr[0] = cp[0] ; cr[1] = cp[1] ; cr[2] = cp[2] ; cr[3] = 0 ;  
      cr = codes2[n++] ;
      cr[0] = cp[1] ; cr[1] = cp[0] ; cr[2] = cp[2] ; cr[3] = 0 ;  
      cr = codes2[n++] ;
      cr[0] = cp[1] ; cr[1] = cp[2] ; cr[2] = cp[0] ; cr[3] = 0 ;  
    }

  /* devellop into single base codes */
  codes2[64] = 0 ;
  for (i = 0 ; i < 64 ; i++)
    {
      cr = codes2[i] ; cr[6] = 0 ;
      for (j = 2 ; j>= 0 ; j--)
	switch ((int)cr[j])
	  {
	  case 'A': cr[2*j] = 'A' ; cr[2*j+1] = 'A' ; break ;
	  case 'T': cr[2*j] = 'T' ; cr[2*j+1] = 'T' ; break ;
	  case 'G': cr[2*j] = 'G' ; cr[2*j+1] = 'G' ; break ;
	  case 'C': cr[2*j] = 'C' ; cr[2*j+1] = 'C' ; break ;

	  case 'K': cr[2*j] = 'A' ; cr[2*j+1] = 'T' ; break ;
	  case 'L': cr[2*j] = 'T' ; cr[2*j+1] = 'A' ; break ;
	  case 'M': cr[2*j] = 'G' ; cr[2*j+1] = 'C' ; break ;
	  case 'N': cr[2*j] = 'C' ; cr[2*j+1] = 'G' ; break ;

	  case 'X': cr[2*j] = 'A' ; cr[2*j+1] = 'G' ; break ;
	  case 'Y': cr[2*j] = 'G' ; cr[2*j+1] = 'A' ; break ;
	  case 'Z': cr[2*j] = 'T' ; cr[2*j+1] = 'C' ; break ;
	  case 'U': cr[2*j] = 'C' ; cr[2*j+1] = 'T' ; break ;

	  case 'V': cr[2*j] = 'A' ; cr[2*j+1] = 'C' ; break ;
	  case 'W': cr[2*j] = 'C' ; cr[2*j+1] = 'A' ; break ;
	  case 'R': cr[2*j] = 'T' ; cr[2*j+1] = 'G' ; break ;
	  case 'S': cr[2*j] = 'G' ; cr[2*j+1] = 'T' ; break ;
	  }
    }

  /* verify that no other word can be chosen at distance 4 */ 
  for (i = 0 ; i < 64 ; i++)
    fprintf (stderr, " %s", codes2[i]) ;
  for (i1= z = 0 ; z < 4096 ; z++)
    {
      for (i = 0, z1 = z ; i < 6 ; z1/=4, i++)
	switch (z1 & 0x3)
	  {
	  case 0: buf[i] = 'A' ; break ;
	  case 1: buf[i] = 'T' ; break ;
	  case 2: buf[i] = 'G' ; break ;
	  case 3: buf[i] = 'C' ; break ;
	  }
      buf[7] = 0 ;

      for (i = i2 = 0, k = 999 ; i < 64 ; i++)
	{
	  cr = codes2[i] ; cp = buf ;
	  for (j = 0, n = 0 ; j < 6 ; cp++, cr++, j++)
	    if (*cr != *cp) n++ ;
	  if (n == 0)
	    i2 = 1 ;
	  else if (n < k)
	    k = n ;
	}
      if (i2 == 0 && k >= 4)
	fprintf (stderr, "Found %d new word %s\n", ++i1, buf) ;
      if (i2 == 1 && k < 4)
	fprintf (stderr, "###Found %d bad word %s\n", ++i1, buf) ;
    }


  aceOutf (ao, "Number of substitutions between the original barcodes\t", nCodes) ;
  if (0)
    for (cpp = libs ; *cpp ; cpp++)
      aceOutf (ao, "\t%s", *cpp) ;  
  aceOutf (ao, "\n\t") ;
  for (cpp = codes2 ; *cpp ; cpp++)
    aceOutf (ao, "\t%s", *cpp) ;
  for (cpp = codes2, lib = libs ; *cpp ; lib++, cpp++)
    { 
      aceOutf (ao, "\n\t%s", *cpp) ;
      for (cqq = codes2 ; *cqq ; cqq++)
	{
	  for (n = 0, cp = *cpp, cq = *cqq ; *cp && *cq ; cp++, cq++)
	   if (*cp != *cq)n++ ;
	  aceOutf (ao, "\t%d", n) ;
	}
    }
  aceOutf (ao, "\n\n") ;

  /* construct all 24 possible permutations */
  n = 0 ;
  for (z = 0, nCodes = 0 ; z < 24 * 24 * 24 * 24 * 24 * 24 ; z++)
    {
      if (0 && z > 1) break ;
      up = arrayp (aa, n++, CODE) ; nCodes++ ;
      for (m = 0 ; m < 64 ; m++) /* 24 words in the code */
	{
	  cr = codes2[m] ; z1 = 24 * z ;
	  for (i1 = 0 ; i1 < 6 ; i1++) /* 6 letters */
	    {
	      z1 /= 24 ; z2 = z1 % 24 ; j = z2 % 6 ; k = z2 / 4 ;
	      cp = permut[j] ;
	      for (i2 = 0 ; i2 < 4 ; i2++) 
		if (cr[i1] == atgc[i2])
		  {  buf[i1] = cp [(i2 + k) % 4] ; break ; }
	    }
	  buf[i1] = 0 ;	      
	  dictAdd (dict, buf, &w) ;
	  up->w[m] = w ; 
	}
      /* sort the code */
      for (m = 0 ; m < 64 ; m++)
	array (bb, m, int) = up->w[m] ;
      arraySort (bb, intOrder) ;
      for (m = 0 ; m < 64 ; m++)
	up->w[m] = array (bb, m, int) ;	 
      if (nCodes % 1000000 == 0)
	{
	  arraySort (aa, codeOrder) ; 
	  arrayCompress (aa) ;
	  n = arrayMax (aa) ;
	}
    }
  arraySort (aa, codeOrder) ;
  arrayCompress (aa) ;
  aceOutf (ao, "Constructed %d 64 words codes\n", nCodes) ;
  aceOutf (ao, "Constructed %d different 64 words codes using %d words\n", arrayMax (aa), dictMax (dict)) ;
  nnn = dnaMotifBarCodeIndelDo (ao, aa, dict, &nCodes) ;
  aceOutf (ao, "Constructed %d different indel resistant %d words codes using %d words\n", nCodes, nnn, dictMax (dict)) ;
      
  aceOutf (ao, "Substitutions between different barcodes\n\t") ;
  for (cpp = libs ; *cpp ; cpp++)
    aceOutf (ao, "\t%s", *cpp) ;  
  aceOutf (ao, "\n\t") ;
  for (cpp = codes ; *cpp ; cpp++)
    aceOutf (ao, "\t%s", *cpp) ;
  for (cpp = codes, lib = libs ; *cpp ; lib++, cpp++)
    { 
      aceOutf (ao, "\n%s\t%s", *lib, *cpp) ;
      for (cqq = codes ; *cqq ; cqq++)
	{
	  for (n = 0, cp = *cpp, cq = *cqq ; *cp && *cq ; n++, cp++, cq++)
	    buf[n] = (*cp != *cq ? *cq : '.') ;
	  buf[n] = 0 ;
	  aceOutf (ao, "\t%s", buf) ;
	}
    }
  aceOutf (ao, "\n\n") ;

  for (ii = j = 0 ; ii < arrayMax (aa) ; ii++)
    {
      up = arrp (aa, ii, CODE) ;
      if (up->indel != nnn)
	continue ;
      for (i = 0 ; i < 64 ; i++)
	if (up->w[i])
	  codes2[j++] = (char *)dictName (dict, up->w[i]) ;
      codes2[j] = 0 ;
      break ;
    }
  aceOutf (ao, "Number of substitutions between different indel resistant barcodes\nOne code is taken as an example but the pattern of distances applies to all of the %d codes\t", nCodes) ;
  if (0)
    for (cpp = libs ; *cpp ; cpp++)
      aceOutf (ao, "\t%s", *cpp) ;  
  aceOutf (ao, "\n\t") ;
  for (cpp = codes2 ; *cpp ; cpp++)
    aceOutf (ao, "\t%s", *cpp) ;
  for (cpp = codes2, lib = libs ; *cpp ; lib++, cpp++)
    { 
      aceOutf (ao, "\n\t%s", *cpp) ;
      for (cqq = codes2 ; *cqq ; cqq++)
	{
	  for (n = 0, cp = *cpp, cq = *cqq ; *cp && *cq ; cp++, cq++)
	   if (*cp != *cq)n++ ;
	  aceOutf (ao, "\t%d", n) ;
	}
    }
  aceOutf (ao, "\n\n") ;


  aceOutf (ao, "Number of indels between different barcodes borders by C------A\n\t") ;
  for (cpp = libs ; *cpp ; cpp++)
    aceOutf (ao, "\t%s", *cpp) ;  
  aceOutf (ao, "\n\t") ;
  for (cpp = codes ; *cpp ; cpp++)
    aceOutf (ao, "\t%s", *cpp) ;

  for (cpp = codes, lib = libs ; *cpp ; lib++, cpp++)
    { 
      aceOutf (ao, "\n%s\t%s", *lib, *cpp) ;
      for (cqq = codes ; *cqq ; cqq++)
	{
	  a = arrayReCreate (a, 12, char) ;
	  b = arrayReCreate (b, 12, char) ;

	  array (a, 0, char) = C_ ;
	  array (b, 0, char) = C_ ;
	  for (n = 1, cp = *cpp, cq = *cqq ; *cp && *cq ; n++, cp++, cq++)
	    {
	      array (a, n, char) = dnaEncodeChar [(int)*cp] ;
	      array (b, n, char) = dnaEncodeChar [(int)*cq] ;
	    }
	  for (; n<30 ; n++)
	    {
	      array (a, n, char) = A_ ;
	      array (b, n, char) = A_ ;
	    }
	  for (; n<40 ; n++)
	    {
	      array (a, n, char) = G_ ;
	      array (b, n, char) = G_ ;
	    }
	  a1 = b1 = 0 ; a2 = b2 = n ;
	  nN = exact = 0 ;
	  err = arrayReCreate (err, 12, A_ERR) ;
	  aceDnaTrackErrors (a, a1, &a2
			     , b, b1, &b2
			     , &nN, err, 1, 3, FALSE, &exact, TRUE) ;
	  aceOutf (ao, "\t") ;
	  if (arrayMax(err) < 3) aceOutf (ao, "%d", arrayMax(err)) ;
	}
    }
  aceOutf (ao, "\n\n") ;

  aceOutf (ao, "Indels between different barcodes borders by C------A\n\t") ;
  for (cpp = libs ; *cpp ; cpp++)
    aceOutf (ao, "\t%s", *cpp) ;  
  aceOutf (ao, "\n\t") ;
  for (cpp = codes ; *cpp ; cpp++)
    aceOutf (ao, "\t%s", *cpp) ;
  for (cpp = codes, lib = libs ; *cpp ; lib++, cpp++)
    { 
      aceOutf (ao, "\n%s\t%s", *lib, *cpp) ;
      for (cqq = codes ; *cqq ; cqq++)
	{
	  a = arrayReCreate (a, 12, char) ;
	  b = arrayReCreate (b, 12, char) ;

	  array (a, 0, char) = C_ ;
	  array (b, 0, char) = C_ ;
	  for (n = 1, cp = *cpp, cq = *cqq ; *cp && *cq ; n++, cp++, cq++)
	    {
	      array (a, n, char) = dnaEncodeChar [(int)*cp] ;
	      array (b, n, char) = dnaEncodeChar [(int)*cq] ;
	    }
	  for (; n<30 ; n++)
	    {
	      array (a, n, char) = A_ ;
	      array (b, n, char) = A_ ;
	    }
	  for (; n<40 ; n++)
	    {
	      array (a, n, char) = G_ ;
	      array (b, n, char) = G_ ;
	    }

	  a1 = b1 = 0 ; a2 = b2 = n ;
	  nN = exact = 0 ;
	  err = arrayReCreate (err, 12, A_ERR) ;
	  aceDnaTrackErrors (a, a1, &a2
			     , b, b1, &b2
			     , &nN, err, 1, 3, FALSE, &exact, TRUE) ;
	  aceOutf (ao, "\t") ;
	  if (arrayMax(err) < 3)
	    {
	      for (ep = arrp (err, 0, A_ERR), n = 0 ; n < arrayMax (err) ; ep++, n++)
		{
		  switch (ep->type)
		    {
		    case TROU: aceOutf (ao, "-%d%c", ep->iShort, dnaDecodeChar [(int)arr (b, ep->iShort, char)]) ; break ;
		    case INSERTION: aceOutf (ao, "+%d%c", ep->iShort, dnaDecodeChar [(int)arr (a, ep->iLong - 1, char)]) ; break ;
		    case ERREUR: aceOutf (ao, "%c%d%c",  dnaDecodeChar [(int)arr (a, ep->iLong, char)], ep->iShort, dnaDecodeChar [(int)arr (b, ep->iShort, char)]) ; break ;
		    default: break ;
		    }
		}
	    }
	}
    }
  aceOutf (ao, "\n\n") ;

  ac_free (ao) ;
  return 0 ;
} /* dnaMotifBarCode */
#endif

/*************************************************************************************/
/*************************************************************************************/

static int dnaMotifBarCodeTestList (LOOK *look)
{
  AC_HANDLE h = ac_new_handle () ;
  int n1 ;
  const char *ccp ;
  DICT *dict = dictHandleCreate (256, h) ;
  ACEIN ai = aceInCreate (look->barCode, FALSE, h) ;

  while (ai && aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (ccp)
	{
	  dictAdd (dict, ccp, &n1) ;
	  if (0)  printf ("Code %d\t%s\n", n1, ccp) ;
	}
    }
  dnaMotifBarDistanceDistrib (look, dict, 4) ;
  dnaMotifBarCodeIndel (look->ao, dict) ;
  ac_free (h) ;
  return 0 ;
} /* dnaMotifBarCodeTestList */

/*************************************************************************************/

/* 96 SOLID BarCodes 2011_10_22
  char *codes2[] ={"tctGTGTAAGAGGctg"
, "tctAGGGAGTGGTctg"
, "tctATAGGTTATActg"
, "tctGGATGCGGTCctg"
, "tctGTGGTGTAAGctg"
, "tctGCGAGGGACActg"
, "tctGGGTTATGCCctg"
, "tctGAGCGAGGATctg"
, "tctAGGTTGCGACctg"
, "tctGCGGTAAGCTctg"
, "tctGTGCGACACGctg"
, "tctAAGAGGAAAActg"
, "tctGCGGTAAGGCctg"
, "tctGTGCGGCAGActg"
, "tctGAGTTGAATGctg"
, "tctGGGAGACGTTctg"
, "tctGGCTCACCGCctg"
, "tctAGGCGGATGActg"
, "tctATGGTAACTGctg"
, "tctGTCAAGCTTTctg"
, "tctGTGCGGTTCCctg"
, "tctGAGAAGATGActg"
, "tctGCGGTGCTTGctg"
, "tctGGGTCGGTATctg"
, "tctAACATGATGActg"
, "tctCGGGAGCCCGctg"
, "tctCAGCAAACTTctg"
, "tctAGCTTACTACctg"
, "tctGAATCTAGGGctg"
, "tctGTAGCGAAGActg"
, "tctGCTGGTGCGTctg"
, "tctGGTTGGGTGCctg"
, "tctCGTTGGATACctg"
, "tctTCGTTAAAGGctg"
, "tctAAGCGTAGGActg"
, "tctGTTCTCACATctg"
, "tctCTGTTATACCctg"
, "tctGTCGTCTTAGctg"
, "tctTATCGTGAGTctg"
, "tctAAAAGGGTTActg"
, "tctTGTGGGATTGctg"
, "tctGAATGTACTActg"
, "tctCGCTAGGGTTctg"
, "tctAAGGATGATCctg"
, "tctGTACTTGGCTctg"
, "tctGGTCGTCGAActg"
, "tctGAGGGATGGCctg"
, "tctGCCGTAAGTGctg"
, "tctATGTCATAAGctg"
, "tctGAAGGCTTGCctg"
, "tctAAGCAGGAGTctg"
, "tctGTAATTGTAActg"
, "tctGTCATCAAGTctg"
, "tctAAAAGGCGGActg"
, "tctAGCTTAAGCGctg"
, "tctGCATGTCACCctg"
, "tctCTAGTAAGAActg"
, "tctTAAAGTGGCGctg"
, "tctAAGTAATGTCctg"
, "tctGTGCCTCGGTctg"
, "tctAAGATTATCGctg"
, "tctAGGTGAGGGTctg"
, "tctGCGGGTTCGActg"
, "tctGTGCTACACCctg"
, "tctGGGATCAAGCctg"
, "tctGATGTAATGTctg"
, "tctGTCCTTAGGGctg"
, "tctGCATTGACGActg"
, "tctGATATGCTTTctg"
, "tctGCCCTACAGActg"
, "tctACAGGGAACGctg"
, "tctAAGTGAATACctg"
, "tctGCAATGACGTctg"
, "tctAGGACGCTGActg"
, "tctGTATCTGGGCctg"
, "tctAAGTTTTAGGctg"
, "tctATCTGGTCTTctg"
, "tctGGCAATCATCctg"
, "tctAGTAGAATTActg"
, "tctGTTTACGGTGctg"
, "tctGAACGTCATTctg"
, "tctGTGAAGGGAGctg"
, "tctGGATGGCGTActg"
, "tctGCGGATGAACctg"
, "tctGGAAAGCGTTctg"
, "tctAGTACCAGGActg"
, "tctATAGCAAAGCctg"
, "tctGTTGATCATGctg"
, "tctAGGCTGTCTActg"
, "tctGTGACCTACTctg"
, "tctGCGTATTGGGctg"
, "tctAAGGGATTACctg"
, "tctGTTACGATGCctg"
, "tctATGGGTGTTTctg"
, "tctGAGTCCGGCActg"
, "tctAATCGAAGAGctg"
, 0} ;
*/
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/* construct a bar code with bccN words of length bccL at distance bccD
   by iteration
   The method should find one such code if it exists, 
   there should be no false negative results
   However for large values of bccL we expect the run time to diverge
*/

const char *bccName (int code, int NL)
{
  char c ;
  static char buf[1024] ;
  int j ;

  for (j = 0 ; j < NL ; j++)
    {
      switch ((code >> ((NL - j - 1)<<1)) & 0x3)
	{
	case 0: c = 'A' ; break ;
	case 1: c = 'T' ; break ;
	case 2: c = 'G' ; break ;
	case 3: c = 'C' ; break ;
	}
      buf[j] = c ;
    }
  buf[j] = 0 ;
  return buf ;
} /* bccName */

/*************************************************************************************/

static void dnaMotifBccFill (LOOK *look, int code, int ii, int NL, int DD, int NN)
{
  unsigned char *dd = look->bccDistances[ii] ;
  unsigned char *mm = look->bccLinks[ii] ;
  int i, j, d, cc ;

  look->bccNfill++ ;
  if (!dd)
    {
      dd = look->bccDistances[ii] = (unsigned char *) messalloc (NN) ;
      mm = look->bccLinks[ii] = (unsigned char *) messalloc (NN) ;
    }
  if (! ii)
    memset (dd, NL, NN) ;
  else
    {
      memcpy (dd, look->bccDistances[ii - 1], NN) ;
      memcpy (mm, look->bccLinks[ii - 1], NN) ;
    }

  if (ii >= look->bccBestN)
    look->bccBestN = ii + 1 ;
  look->bccCodes[ii] = code ;
  look->bccLk[ii] = mm[code] ;

  dd[code] = 0 ;
  for (i = 0 ; i < NN ; i++)
    {
      cc = code ^ i ;
      if (dd[i] > DD - 1)
	{
	  for (j = d = 0 ; cc && j < NL ; j++, cc >>= 2)
	    if (cc & 0x3) d++ ;
	  if (dd[i] > d) dd[i] = d ;      
	  if (d == DD && mm[i] < 250) 
	    mm[i]++ ;
	}
    }
  return ;
} /* dnaMotifBccFill */

/*************************************************************************************/

static BOOL dnaMotifBccLoop (LOOK *look, int ii, int NL, int DD, int NN, int NC)
{
  int i, code, codeMax = NN - 1, maxD = 0, d, DD1, mmMax = 0 ;
  unsigned char *dd, *mm ;

  dd = look->bccDistances[ii - 1] ;
  mm = look->bccLinks[ii - 1] ;
  for (code = 0 ; code < NN ; code++)
    {
      d = dd[code] ;
      if (d  && d > maxD)
	maxD = d ;
      if (dd[code] && mm[code] >= mmMax) 
	{ codeMax = code ; mmMax = mm[code] ; }
    }
  if (debug)
    for (code = 0 ; code < NN ; code++)
      if (dd[code] && mm[code] == mmMax)
	{
	  for (i = 0 ; i < ii ; i++)  fprintf(stderr, "  ") ;
	  fprintf(stderr, "%d : %s=%d   mmMax=%d\n", ii,  bccName(code, NL), code, mmMax) ;
	  break ; 
	}
  for (DD1 = DD ; DD1 <= DD ; DD1++)
    for (code = 0 ; code <= codeMax ; code++)
      {
	d = dd[code] ;
	if (d < DD ||  mm[code] < mmMax)
	  continue ;
	if (dd[code] == DD1)
	  {
	    if (debug) 
	      {
		for (i = 0 ; i < ii ; i++)  fprintf(stderr, "  ") ;
		fprintf (stderr, "%d : %s  %d ==?  %d\n", ii, bccName(code, NL), mm[code], mmMax) ;
	      }
	    if (ii == NC - 1)
	      {
		look->bccCodes[ii] = code ;
		look->bccLk[ii] = mm[code] ;
		look->bccBestN = NC ;
		return TRUE ; /* success */
	      }
	    dnaMotifBccFill (look, code, ii, NL, DD, NN) ;
	    if (dnaMotifBccLoop (look, ii+1, NL, DD, NN, NC))
	      return TRUE ;
	  }
      }
  if (! look->bccSlow && maxD < DD)
    return TRUE ;
  return FALSE ; /* no more word at distance DD */
} /* dnaMotifBccLoop */

/*************************************************************************************/

static void dnaMotifBccConstruct  (LOOK *look)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (look->oFileName, ".code", 0, h) ;
  int i, code ;
  int NL = look->bccL ;
  int NN = (1 << (NL << 1)) ;
  int DD = look->bccD ;
  DICT *dict = dictHandleCreate (look->bccN, h) ;

  look->bccDistances = halloc (NN * sizeof (char *), h) ;
  look->bccLinks = halloc (NN * sizeof (char *), h) ;
  look->bccCodes = halloc (look->bccN * sizeof (int), h) ;
  look->bccLk = halloc (look->bccN * sizeof (int), h) ;

  for (i = code = 0 ; i < look->bccD ; i++)
    { code <<=2 ; code |= 0x1 ; }
  dnaMotifBccFill (look, 0, 0, NL, DD,NN) ;
  dnaMotifBccFill (look, code, 1, NL, DD, NN) ;
  dnaMotifBccLoop (look, 2, NL, DD, NN, look->bccN) ;
  fprintf (stderr, "Found %d codes of length %d at distance %d in %d iterations\n"
	   , look->bccBestN, look->bccL, look->bccD, look->bccNfill
	   ) ;
  for (i = 0 ; i < look->bccBestN ; i++)
    {
     code = look->bccCodes[i] ;
     aceOutf (ao, "%d\t%s\t%d\n", i + 1, bccName(code, NL), look->bccLk[i]) ;
     dictAdd (dict, bccName(code, NL), 0) ;
    }
	
  dnaMotifBarDistanceDistrib (look, dict, look->bccD) ;
  ac_free (ao) ;
  ac_free (h) ;
  return ;
}  /* dnaMotifBarCodeConstruct */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
 
static int parseFastaFile (ACEIN ai, LOOK *look)
{
  int np, nb = 0, newNam, count = 1, ntags = 0 ;
  PP *seq = 0 ;
  BOOL last = FALSE ;
  const char *ccp ;
  char *cp, *cr ;
  Array dna = look->dna ;
  Array seqs = look->probes ;

  np = arrayMax (seqs) ;
  look->probeLengthMin = ACEDB_MAXINT ;

  if (look->motif)
    {
      dictAdd (look->dict, "motif", &newNam) ;

      seq = arrayp (seqs, np++, PP) ;
      seq->nam = newNam ;
      seq->ignore = 2 ;
      seq->isDown = TRUE ;
      seq->start = NN ;
      for (ccp = look->motif ; *ccp; ccp++)
	array (dna, NN++, char) = dnaEncodeChar [(int)*ccp] ;
      seq->ln = NN ;
      seq->stop = NN ;  /* first _Z following the probe */
      array (dna, NN++, char) = 0 ; /* 00terminate each sequence */
      array (dna, NN++, char) = 0 ;
    }
  while (TRUE)
    {
      if (aceInCard (ai))
	{
	  cp = aceInWord (ai) ;
	  if (!cp || !*cp)
	    continue ;
	}
      else
	last = TRUE ;
      if (last || *cp == '>')
	{
	  /* register the probe and
	   * add a null base to propagate correctly bacwards
	   * the clause w[i] cannot be repeated if (w-1)[i+1] is not
	   */
	  if (!last)
	    {
	      dictAdd (look->dict, cp+1, &newNam) ;	  
	      count = fastcMultiplicity (cp, 0, 0) ;
	    }

	  if (seq)
	    {
	      array (dna, NN++, char) = 0 ; /* 00terminate each sequence */
	      array (dna, NN++, char) = 0 ;
	      seq->stop = NN - 1 ;
	      ntags += seq->count ;
	      nb += seq->count * seq->ln ;
	    }

	  if (seq && seq->ln)
	    {
	      if (seq->ln < look->probeLengthMin)
		look->probeLengthMin = seq->ln ;
	      if (seq->ln > look->probeLengthMax)
		look->probeLengthMax = seq->ln ;
	    }
	  
	  if (seq && look->hexamers && !justStrand)
	    { /* add the strand inverted probe */
	      char *cr ;
	      int j, n0 = NN ;
	      array (dna, NN + seq->ln + 2, char) = 0 ; /* make room */
	      for (cr = arrp (dna, seq->stop-2, char), j = seq->ln ; j-- ; cr--)
		array (dna, NN++, char) = complementBase(*cr) ;
	      array (dna, NN++, char) = 0 ;
	      array (dna, NN++, char) = 0 ;
	      seq = arrayp (seqs, np++, PP) ;
	      seq->nam = (seq-1)->nam ;
	      seq->isDown = FALSE ;
	      seq->ln = (seq-1)->ln ;
	      seq->start = n0 ;
	      seq->stop = NN - 1 ;  /* first _Z following the probe */
	    }
	  if (last)
	    break ;
	  seq = arrayp (seqs, np++, PP) ;
	  seq->nam = newNam ;
	  seq->count = count ;
	  if (justType && ! pickMatch (dictName(look->dict,newNam), justType))
	     { seq->ignore = 1 ; look->nIgnored++ ; }
	  seq->isDown = TRUE ;
	  seq->start = NN ;
	}
      else
	{
	  seq->ln += strlen(cp) ;
	  for (cr = cp ; *cr ; cr++)
	    array (dna, NN++, char) = dnaEncodeChar [(int)*cr] ;
	}
    }
  look->np = np ;
  look->ntags += ntags ;
  look->nbp = nb ;
  return np ;
} /* parseFastaFile  */

/*************************************************************************************/
/*************************************************************************************/

#include "topology.h"
static int dnaMotifg26 (LOOK *look)
{
  int nn = 0, nl = 0, nv = 0 ;
  const char *ccp ;
  int k1, k2 ;
  VERTEX *v ;
  LINK *w ;
  AC_HANDLE h = ac_new_handle () ;
  Array links = arrayHandleCreate (100000, LINK, h) ;
  Array vertices = arrayHandleCreate (100000, VERTEX, h) ;
  DICT *dict = dictHandleCreate (100000, h) ;
  ACEIN ai = aceInCreateFromFile (look->g26, "r", 0, h) ;

  while (ai && aceInCard (ai))
    {
      if ((ccp = aceInWord (ai)))
	{
	  dictAdd (dict, ccp, &k1) ;
	  if ((ccp = aceInWord (ai)))
	    {
	      dictAdd (dict, ccp, &k2) ;
	      v = arrayp (vertices, nv++, VERTEX) ;
	      w = arrayp (links, nl++, LINK) ;
	      w->a = k1 ; w->b = k2 ;
	      v->a = k1 ;
	    }
	}
    }
  arraySort (links, topoLinkOrder) ;
  arraySort (vertices, topoVertexOrder) ;
  arrayCompress (links) ;
  arrayCompress (vertices) ;
  nn = topoConnectedComponents (links, vertices) ;
  printf ("Found %u nodes and %u links forming %d clusters\n",
	  arrayMax (vertices), arrayMax(links), nn) ;
  ac_free (h) ;
  return nn ;
} /* dnaMotifg26 */

/*************************************************************************************/
/*************************************************************************************/
/* split the peptide on K or R except not on KP and not on RP
 * avoid incomplete fragments so just consider ]R,R] or ]R,stop] or [met,R]
 */
static int dnaMotifPepSplit (LOOK *look)
{
  int ii, n1 = 0, n2 = 0 ;
  BOOL ok ;
  char cc, cc1, *cp, *cq, *cp0 ;
  char oldSeq[1001], newSeq[1001] ;
  AC_HANDLE h = ac_new_handle () ;
  ACEIN aBig, aRef ;
  DICT *dict1 = dictHandleCreate (1000000, h) ;
  DICT *dict2 = dictHandleCreate (1000000, h) ;
  vTXT txt = vtxtHandleCreate (h) ;

  memset (oldSeq, 0, sizeof (oldSeq)) ;
  memset (newSeq, 0, sizeof (newSeq)) ;

  aBig = aceInCreateFromFile (look->f1, "r", 0, h) ;
  if (! aBig)
    messcrash ("Cannot open -f1 file %s", look->f1) ;
  aRef = aceInCreateFromFile (look->f2, "r", 0, h) ;
  if (! aRef)
    messcrash ("Cannot open -f2 file %s", look->f2) ;
  
  ok = 1 ;
  while (ok)
    {
      if (aceInCard (aRef))
	cp = aceInWord (aRef) ;
      else
	{ cp = 0 ; ok = 0 ; }
      if (!cp || ! *cp || *cp == '>')
	{
	  n1++ ;
	  /* add in the dict all n-mers of the previous sequence */ 
	  vtxtPrint (txt, "*") ; /* add a terminal stop for all RefSeq 
				    This is probably abusive
but the effect is to export less new petides, so at worst it will
generate false negatives
				 */
	  cp = cp0 = vtxtPtr (txt) ;
	  while (cp && *cp)
	    {
	      for (cq = cp ; *cq ; cq++)
		if (*cq == '*' || ((*cq == 'R' || *cq == 'K') && *(cq+1) != 'P'))
		  break ;
	      cc = 0 ;
	      if (*cq)
		{ cq++ ; cc = *cq ; *cq = 0 ; }
	      if (cp && *cp)
		dictAdd (dict1, cp, 0) ;
	      *cq = cc ;
	      cp = cq ;
	    }
	  vtxtClear (txt) ;
	  continue ;
	}
      vtxtPrint (txt, cp) ;
    }
  fprintf (stderr, "Found %d sequences, %d dictinct cleaved peptides in file %s\n"
	   , n1, dictMax (dict1), look->f2) ;

  ok = 1 ;
  while (ok)
    {
      if (aceInCard (aBig))
	cp = aceInWord (aBig) ;
      else
	{ cp = 0 ; ok = 0 ; }
      if (cp && *cp == '>')
	strncpy (newSeq, cp+1, 1000) ;
      if (!cp || ! *cp || *cp == '>')
	{
	  n2++ ;
	  /* add in the dict all n-mers of the previous sequence */ 

	  cp = cp0 = vtxtPtr (txt) ;
	  cq = cp ? strstr (cp+1, "*") : 0 ;
	  if (cq && cq < cp + strlen(cp) - 1) cp = 0 ; /* abandon peptides with an internal stop */
	  if (cp && cp[0] == '*' && cp[1] != 'M') cp = 0 ; /* abandon non met supposed start */
	  while (cp && *cp)
	    {
	      for (cc1 = 0, cq = cp ; *cq ; cq++)
		if (*cq == '*' || ((*cq == 'R' || *cq == 'K') && *(cq+1) != 'P'))
		  { cc1 = *cq ; break ; }
	      cc = 0 ;
	      if (*cq)
		{ cq++ ; cc = *cq ; *cq = 0 ; }
	      if (cp && *cp &&     /* forget very very short sequences */
		  !dictFind (dict1, cp, 0) &&        /* reject known sequences */
		  cp > cp0 &&                          /* reject is not strating on NH2-complete or or cleavage */
		  cc1 &&     /* end on stop or cleavage */
		  ! strstr(cp, "X") &&   ! strstr(cp, "-")
		  )
		{
		  int u = dictAdd (dict2, cp, 0) ;
		  if (u && strlen(cp) <= 3)
		    {
		      *cq = cc ; fprintf (stderr, "%s\t%s\n", cp, cp0) ;
		    }
		}
	      *cq = cc ;
	      cp = cq ;
	    }
	  vtxtClear (txt) ;
	  continue ;
	}
      vtxtPrint (txt, cp) ;
    }

  fprintf (stderr, "\n# REFSEQ\n") ;
  if (0)
    for (ii = 1 ; ii <= dictMax (dict1) ; ii++)
      aceOutf (look->ao, "##%s#\n", dictName (dict1, ii)) ;

  fprintf (stderr, "\n# NEW PEPTIDES\n") ;
  for (ii = 1 ; ii <= dictMax (dict2) ; ii++)
    aceOutf (look->ao, "%s\n", dictName (dict2, ii)) ;

  fprintf (stderr, "Found %d sequences, exported new %d cleaved peptides,  in file %s\n"
	   , n2, dictMax (dict2), look->f1) ;


  ac_free (h) ;
  return n2 ;
} /* dnaMotifPepSplit */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *err)
{
  fprintf (stderr, 
	   "// Usage: dnamotif -p fasta_file -... \n"
	   "// Example:  dnamotif -p seq.fasta -motif cgatcg/gcatgc\n"  
	   "// -p : the name of a fasta file to analyse \n"
	   "//      (tested with 1 Million sequences, upper limit depends on harware)\n"
	   "// -filter filter : keep only the sequence whose name match this filter, ? and * wild chars are recognized\n"
	   "// -o fileName : outfile name, same as redirecting stdout\n"
	   "// -title 'title' : This title is given inside the postscript output files\n"
	   
	   "//\n"  
	   "// Simple histograms:\n"
	   "//   -triplet : distribution of triplets of letters in a fasta file\n"
	   "//   -tripletq : distribution of triplets of letters from fastq file\n"
	   "//   -homopolymer : count homopolymer at either end, fasta file\n"
	   "//   -homopolymerq : count homopolymer at either end, fastq file\n"
	   "//\n"  
	   "//   -hexamers [nn: default 6] [-strand]: count all nn-mer motifs in '-width' central bases\n"
	   "//   -Rhexamers n -c1 c1 -c2 c2 -a1 a1 -a2 a2 -txt tours: recursivelly find the n-mers most expressed in segment c1-c2 relative to a1-a2\n"
	   "//\n"  
	   "// Motif analysis\n"  
	   "//   -motif atgcgct : motif used to seed the profile or in correl, _ means n, / means multi motif ex: atgcnnat/ttc___gcc \n"
	   "//   -force int : number of centered matches to the seed required\n"
	   "//   -profile fileName : a motif profile is generated, exported there, then gnuPlotted in fileName.ps\n"

	   "//   -fold : try to fold the motif(s)\n"
	   "//   -centrage : export the centrage histo of the motif(s)\n"
  
	   "//   -correl : export the correlation function of the motif(s) listed as -motif x/y...\n"
	   "//   -doubleGas : specialized distribution of distances between double gas motifs)\n"
	   "//   -width nn :  width in bp, used in centrage/correl\n"

	   "// -g26 : hard coded g26 analysis\n"
	   "// -barCode : hard coded bar code analysis\n"
	   "// -barCodeCreate -bccL <int> -bccD <int> -bccN <int> [-bccSlow]\n"
	   "//    construct by iteration an error correcting code on ATGC\n"
	   "//    with at most bccN words (default 1024) of length bccL at distance bccD\n" 
	   "//    In fast mode, we recover half of Bogdanova table\n"
	   "//    In -bccSlow mode, we recover several others but obscenelly slow\n"

	   "// -pepSplit int file1 file2\n"
	   "//     file1 is a petide fasta file, example all_pep_for_mass-spec.fasta\n"
	   "//     file2 is a petide fasta file, example refseq.peptide.fasta\n"
	   "//     the program exports a difference peptide fasta file\n"
	   "//     excluding from file1 the n-mers found in file2\n"
	   "//     The idea is to facilitate a mass-spec reanalysis of the orphan profiles\n"
	   ) ;
  
  if (err)
    fprintf (stderr, "\n// ERROR: %s\n", err) ;
  exit (1) ;
}

/*************************************************************************************/

static void myTest (void)
{
  DICT *dict = dictCaseSensitiveHandleCreate (10000, 0) ;
  char *cp, buf[1001], *cr ;
  const char *zz[100000] ;
  const char *cq ;
  int n = 0, nn = 0, i, k = 0, zn[100000] ;
  long int nc = 0 ;

  while ((cp = fgets (buf, 1000, stdin)))
    {
      nn++ ;
      while (*cp== ' ' || *cp == '\t') cp++ ;
      cr = cp + strlen(cp) - 1 ;
      while (cr > cp && (*cr == ' '  || *cr == '\n' || *cr == '\t'))
	*cr-- = 0 ;
      nc += strlen (cp) ;
      dictAdd (dict, cp, &n) ;
      cq = dictName (dict, n) ;
      if (cp && cq && strcmp (cq, cp))
	messcrash ("nn=%d n= %d \ncp = #%s#  \ncq = #%s#\n", nn, n, cp, cq) ;

      if ((nn % 10000000) == 0)
	{
	  if (0) printf ("nn=%d n= %d cp =%s  cq = %s", nn, n, cp, cq) ;
	  printf (" ok %d Millions words %ld M char\n", nn/1000000, nc/1000000) ;
	  dictStatus (dict) ;
	}
      if ((nn % 1000000) == 0)
	{ 
	  k = nn/1000000 ;
	  zz[k] = cq ; zn[k] = n ;
	  for (i = 1 ; i <= k ; i++)
	    if (zz[i] != dictName (dict, zn[i]))
	      messcrash ("pointer has moved i=%d n = %d zn[i]=%d\n", i, n, zn[i]) ; 
	  printf (".") ; fflush (stdout) ; 
	}
    }
  printf ("done %d\n", nn) ;
  dictStatus (dict) ;
  exit (0) ;
} /* myTest */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = 0 ;
  int mx = 0, nFound ;
  const char *pFileName = 0 ;
  const char *oFileName = 0 ;
  LOOK look ;

  freeinit () ;
  if (0) myTest () ;

  memset (&look, 0, sizeof(look)) ;
  look.h = h ;
  look.dict = dictHandleCreate (1000000, h) ;
  look.probes = arrayHandleCreate (100000, PP, h) ;
  look.dna = arrayHandleCreate (1000000, char, h) ;
  
  /* open and parse the probe and the target file */
  getCmdLineOption (&argc, argv, "-g26", &look.g26) ;
  getCmdLineOption (&argc, argv, "-filter", &justType) ;
  getCmdLineOption (&argc, argv, "-motif", &look.motif) ;
  getCmdLineOption (&argc, argv, "-cFileName", &look.cFileName) ;
  getCmdLineOption (&argc, argv, "-title", &(look.title)) ;
  getCmdLineOption (&argc, argv, "-profile", &(look.profileFileName)) ;
  getCmdLineInt (&argc, argv, "-width", &look.width) ;
  getCmdLineInt (&argc, argv, "-nIter", &look.nIter) ;
  if (getCmdLineInt (&argc, argv, "-pepSplit", &look.pepSplit))
    {
      if (getCmdLineOption (&argc, argv, "-f1", &look.f1)  &&
	  getCmdLineOption (&argc, argv, "-f2", &look.f2)
	  ) ;
      else
	usage ("Missing argumnets -f1 -f2, needed for -pepSplit") ;
    }

  getCmdLineOption (&argc, argv, "-barCode", &look.barCode) ;
  if ((look.barCodeCreate = getCmdLineOption (&argc, argv, "-barCodeCreate", 0)))
    {
      look.bccSlow = getCmdLineOption (&argc, argv, "-bccSlow", 0) ;
      if (! getCmdLineInt (&argc, argv, "-bccL", &look.bccL) ||
	  ! getCmdLineInt (&argc, argv, "-bccD", &look.bccD) ||
	  (! getCmdLineInt (&argc, argv, "-bccN", &look.bccN) && look.bccSlow) ||
	  look.bccL <= 0 || (look.bccSlow && look.bccN <= 0) || look.bccD <= 0
	  )
	messcrash ("-barCodeCreate requires positive values for -bccL -bccD -bccN, try -help") ;
      if (! look.bccN)
	look.bccN = 1024 ;
    }
  look.triplet = getCmdLineOption (&argc, argv, "-triplet", 0) ;
  look.tripletq = getCmdLineOption (&argc, argv, "-tripletq", 0) ;
  look.homopolymer = getCmdLineOption (&argc, argv, "-homopolymer", 0) ;
  look.homopolymerq = getCmdLineOption (&argc, argv, "-homopolymerq", 0) ;

  /* force the presence of that many letters */
  getCmdLineInt (&argc, argv, "-force", &look.force) ;
  getCmdLineInt (&argc, argv, "-mx", &mx) ;
  look.doubleGas = getCmdLineOption (&argc, argv, "-doubleGas", 0) ;
  look.correl = getCmdLineOption (&argc, argv, "-correl", 0) ;
  look.fold = getCmdLineOption (&argc, argv, "-fold", 0) ;
  if (!getCmdLineInt (&argc, argv, "-hexamers", &look.hexamers)) 
    {
      if (getCmdLineOption (&argc, argv, "-hexamers", 0))
	look.hexamers = 6 ;
    }

  look.centrage = getCmdLineOption (&argc, argv, "-centrage", 0) ;
  justStrand = getCmdLineOption (&argc, argv, "-strand", 0) ;
  justAntistrand = getCmdLineOption (&argc, argv, "-antistrand", 0) ;

  if (getCmdLineInt (&argc, argv, "-Rhexamers", &look.Rhexamers)) 
    {
      look.RhexamersTours = 8 ;
      getCmdLineInt (&argc, argv, "-rxt", &look.RhexamersTours) ;

      if (getCmdLineInt (&argc, argv, "-c1", &look.RhexamersC1) &&
	  getCmdLineInt (&argc, argv, "-c2", &look.RhexamersC2) &&
	  getCmdLineInt (&argc, argv, "-a1", &look.RhexamersA1) &&
	  getCmdLineInt (&argc, argv, "-a2", &look.RhexamersA2)
	  ) ;
      else
	messcrash ("-Rhexamers takes 5 numeric parameters: motif lentgth, coordinate of central segment, coordinate of external segment") ;
      justStrand = TRUE ; justAntistrand = FALSE ;
    }


  if (! look.barCode && ! look.barCodeCreate && !look.g26 && !look.pepSplit &&
      !getCmdLineOption (&argc, argv, "-p", &pFileName))
    usage ("Missing '-p ' option, the name of a fasta file is not specified, use -p - to analyse stdin") ;

  if (pFileName)
    {
      if (!strcmp (pFileName, "-"))
	ai = aceInCreateFromStdin (FALSE, 0, h) ;
      else
	ai = aceInCreateFromFile (pFileName, "r", 0, h) ;
      if (!ai)
	usage ("Sorry i cannot open the probe fasta file") ;
      else
	parseFastaFile (ai, &look) ;
      ac_free (ai) ;
    }


  getCmdLineOption (&argc, argv, "-o", &oFileName) ;
  look.ao =  aceOutCreate (oFileName, "", FALSE, look.h) ;
  look.oFileName = oFileName ;

  fprintf (stderr, "// dnamotif starts %s, %d ignored, %d sequences, %d tags, filter %s, %d bp, length %d to %d bp\n"
	   , timeShowNow()
	   , look.nIgnored
	   , look.np - look.nIgnored, look.ntags
	   , justType ? justType : "none"
	   , look.nbp
	   , look.probeLengthMin, look.probeLengthMax
	   ) ;
  if (mx)  look.probeLengthMax = mx ;
  if (look.barCode)
    {
      if (look.bccL == 12 && look.bccD == 5) 
	dnaMotifBarCodeBogdanova_12_5 (&look) ;
      else if (look.bccL == 8 && look.bccD == 4) 
	dnaMotifBarCodeBogdanova_8_4 (&look) ;
      else 
	dnaMotifBarCodeTestList (&look) ;

      if (0) dnaMotifBarCodeConstruct (&look) ;
    }
  else if (look.barCodeCreate)
    {
      dnaMotifBccConstruct (&look) ;
    }
  else if (look.pepSplit)
    nFound = dnaMotifPepSplit (&look) ;
  else if (look.g26)
    nFound = dnaMotifg26 (&look) ;
  else if (look.profileFileName)
    nFound = dnaMotif (&look) ;
  else if (look.doubleGas)
    nFound = dnaDoubleGas (&look) ;
  else if (look.correl || look.centrage)
    nFound = dnaMotifCorrel (&look) ;
  else if (look.fold)
    nFound = dnaFoldMotifs (&look) ;
  else if (look.hexamers)
    nFound = dnaCountHexamers (&look) ;
  else if (look.Rhexamers)
    nFound = dnaCountRHexamers (&look) ;
  else if (look.triplet || look.tripletq) 
    dnaTripletDistrib (&look) ;   /* distrib of triplets in a fastq file */
  else if (look.homopolymer || look.homopolymerq)
    dnaHomopolymerDistrib (&look) ;   /* distrib of triplets in a fastq file */
  else if (0)
    dnaPowerTable (&look) ;
  fprintf (stderr, "// dnamotif done %s, using seed %s. Found = %d motifs \n ", timeShowNow (), look.motif ? look.motif : "none", nFound) ;

  ac_free (look.h) ;
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/* motif danielle IRF4: 
        acaggaagtgagaa
           GGAAG 
   dnamotif -motif  acaggaagtgagaa -p peaks.IRF4mIL21.C.fasta
         gggggngnnacaggaagtgagaaggggggggggg  
*/
/*
 foreach i (3 4 5 6 7 8 9 10 11 12)
   foreach j (3 4 5 6 7 8 9 10 11 12)
     if($i >= $j && ($i <= 10 || $j >= 5) && ! -e CODES/code.L$i.D$j.code) then
       waligner/scripts/submit CODES/code.L$i.D$j "dnamotif -barCodeCreate -bccL $i -bccD $j -bccN 40000 -o CODES/code.L$i.D$j" 32G
     endif
   end
 end

 set toto="codes.txt"
 date > $toto
 foreach j (3 4 5 6 7 8 9 10 11 12)
   echo -n "\t$j" >> $toto
 end
 foreach i (3 4 5 6 7 8 9 10 11 12)
   echo -n "\n$i" >> $toto
   foreach j (3 4 5 6 7 8 9 10 11 12)
     if (-e code.L$i.D$j.code) then
       cat code.L$i.D$j.code | gawk '{n++}END{printf("\t%d",n);}' >> $toto
     else
       echo -n "\t0" >> $toto
     endif
   end
 end
 echo >> $toto
 cat $toto
 echo $toto
 

*/
