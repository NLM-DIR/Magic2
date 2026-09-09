/*  File: peaks.c
 *  Implementation of peaks.h, see that file for the method.
 *
 *  Pipeline
 *    1  sigma      robust noise level, MAD of successive differences
 *    2  background morphological opening (running min then running max)
 *                  followed by a moving average
 *    3  fs         lightly smoothed copy of f, used for DETECTION only
 *    4  components maximal runs of  fs > background + n*sigma,
 *                  then extended outwards down the flanks
 *    5  merge tree union-find over the sublevel sets of fs inside each
 *                  component -> every local max with its prominence
 *    6  selection  prominence > max (n*sigma, frac*height)
 *    7  domains    consecutive kept apices are split at the deepest
 *                  point between them
 *    8  window     half-mean fixed point iteration on g = f - localBase,
 *                  areas taken from a prefix sum of the RAW data
 *
 *  Everything internal is double; only the input is int.
 */

#include "peaks.h"
#include <math.h>
#include "ac.h"
/*======================================================================*/
/*=========================== small helpers ============================*/
/*======================================================================*/

static int pkDoubleOrder (const void *va, const void *vb)
{
  double a = *(const double *) va, b = *(const double *) vb ;
  return a < b ? -1 : (a > b ? 1 : 0) ;
}

static double pkMedian (Array a)   /* a : Array of double, destroyed order */
{
  int n = arrayMax (a) ;
  if (!n) return 0 ;
  arraySort (a, pkDoubleOrder) ;
  return (n & 1) ? arr (a, n/2, double)
                 : 0.5 * (arr (a, n/2 - 1, double) + arr (a, n/2, double)) ;
}

/*----------------------------------------------------------------------*/
/* Noise level measured AT THE SCALE OF THE PEAKS, not at the scale of
 * one bin.
 *
 * A peak of length L does not compete against bin-to-bin noise, it
 * competes against whatever the background does over a stretch of L
 * bins.  For white noise the two are equivalent: averaging L bins
 * shrinks the fluctuation as sigma/sqrt(L), which is the assumption
 * behind the A/sqrt(L) matched filter.  On a smoothed track, or on deep
 * data where the background itself undulates, the assumption fails:
 * adjacent bins are correlated, the windowed mean does not shrink as
 * fast, and a bin-level estimate understates the real variability by a
 * large factor.
 *
 * We therefore average the residual f - background over non overlapping
 * windows of L bins and take the MAD of that, rescaled by sqrt(L):
 *
 *    sigma = 1.4826 * median |mean_L (f - b)| * sqrt(L)
 *
 * On white noise the sqrt(L) cancels the averaging exactly and this
 * returns the bin level sigma, so nothing changes for well behaved data.
 * The windows containing peaks give large values, which the median
 * ignores as long as peaks occupy well under half the profile.
 */

static double pkNoiseSigma (Array fd, Array bg, int L)
{
  int i, nd = 0, n = arrayMax (fd) ;
  double sigma ;
  Array s = 0, d = 0 ;

  if (L < 2) L = 2 ;
  if (L > n / 8) L = n / 8 ;
  if (L < 2 || n < 32) return 1.0 ;

  s = arrayCreate (n + 1, double) ;
  array (s, 0, double) = 0 ;
  for (i = 0 ; i < n ; i++)
    array (s, i+1, double) = arr (s, i, double)
      + arr (fd, i, double) - arr (bg, i, double) ;

  d = arrayCreate (n / L + 2, double) ;
  for (i = 0 ; i + L <= n ; i += L)
    array (d, nd++, double) =
      fabs (arr (s, i+L, double) - arr (s, i, double)) / L ;

  if (nd < 16) { arrayDestroy (s) ; arrayDestroy (d) ; return 1.0 ; }
  sigma = 1.4826 * pkMedian (d) * sqrt ((double) L) ;
  arrayDestroy (s) ;
  arrayDestroy (d) ;

  /* f is integer: a noise level below one count is meaningless. */
  if (sigma < 1) sigma = 1 ;
  return sigma ;
}

/*----------------------------------------------------------------------*/
/* Sliding window minimum or maximum in O(n) with a monotonic deque.
 * Window is centred, of total width 2*(w/2)+1, clipped at the edges.
 */

static void pkRunExtremum (Array src, Array dst, int w, BOOL wantMin)
{
  int i, j, head = 0, tail = 0, n = arrayMax (src), half = w / 2 ;
  Array dq = arrayCreate (n + 1, int) ;

  if (half < 1) half = 1 ;
  for (i = 0 ; i < n + half ; i++)
    {
      if (i < n)
	{
	  double v = arr (src, i, double) ;
	  while (tail > head)
	    {
	      double u = arr (src, arr (dq, tail-1, int), double) ;
	      if (wantMin ? (u >= v) : (u <= v)) tail-- ;
	      else break ;
	    }
	  array (dq, tail++, int) = i ;
	}
      j = i - half ;               /* the output slot now complete */
      if (j >= 0)
	{
	  while (arr (dq, head, int) < j - half) head++ ;
	  array (dst, j, double) = arr (src, arr (dq, head, int), double) ;
	}
    }
  arrayDestroy (dq) ;
}

/*----------------------------------------------------------------------*/
/* Centred moving average, O(n) via a prefix sum. */

static void pkMovingAverage (Array src, Array dst, int w)
{
  int i, n = arrayMax (src), half = w / 2, i0, i1 ;
  Array s = arrayCreate (n + 1, double) ;

  if (half < 1) { for (i = 0 ; i < n ; i++) array (dst, i, double) = arr (src, i, double) ;
                  arrayDestroy (s) ; return ; }
  array (s, 0, double) = 0 ;
  for (i = 0 ; i < n ; i++)
    array (s, i+1, double) = arr (s, i, double) + arr (src, i, double) ;
  for (i = 0 ; i < n ; i++)
    {
      i0 = i - half ; if (i0 < 0) i0 = 0 ;
      i1 = i + half ; if (i1 > n - 1) i1 = n - 1 ;
      array (dst, i, double) =
	(arr (s, i1+1, double) - arr (s, i0, double)) / (i1 - i0 + 1) ;
    }
  arrayDestroy (s) ;
}

/*----------------------------------------------------------------------*/
/* Slow background: a morphological opening (erode then dilate) sits
 * under every structure narrower than the window, the moving average
 * then removes the staircase left by the two extremum filters.
 */

static void pkBackground (Array fd, Array bg, int w)
{
  int i, n = arrayMax (fd) ;
  double delta ;
  Array e = arrayCreate (n, double) ;
  Array d = arrayCreate (n, double) ;
  Array r = arrayCreate (n, double) ;

  pkRunExtremum (fd, e, w, TRUE) ;     /* erosion  */
  pkRunExtremum (e,  d, w, FALSE) ;    /* dilation */
  pkMovingAverage (d, bg, w) ;

  /* The opening is a LOWER envelope: it tracks the bottom of the noise,
   * not its centre, and would push half the samples above any threshold.
   * Lift it by the median residual, which the peaks cannot bias as long
   * as they occupy well under half the samples - the whole premise here.
   */
  for (i = 0 ; i < n ; i++)
    array (r, i, double) = arr (fd, i, double) - arr (bg, i, double) ;
  delta = pkMedian (r) ;
  for (i = 0 ; i < n ; i++)
    arr (bg, i, double) += delta ;

  arrayDestroy (e) ;
  arrayDestroy (d) ;
  arrayDestroy (r) ;
}

/*======================================================================*/
/*======================= merge tree / prominence ======================*/
/*======================================================================*/

typedef struct pkNodeStruct {
  int    apex ;        /* argmax of this branch                          */
  int    saddle ;      /* where it died, -1 for the root                 */
  double birth ;       /* fs at the apex                                 */
  double death ;       /* fs at the saddle                               */
  int    parent ;      /* union find                                     */
  BOOL   isRoot ;
  BOOL   kept ;
} PKNODE ;

typedef struct { int x ; double v ; } PKXV ;

static int pkXVOrder (const void *va, const void *vb)   /* value descending */
{
  const PKXV *a = (const PKXV *) va, *b = (const PKXV *) vb ;
  if (a->v > b->v) return -1 ;
  if (a->v < b->v) return  1 ;
  return a->x - b->x ;
}

static int pkIntOrder (const void *va, const void *vb)
{ return *(const int *)va - *(const int *)vb ; }

static int pkPeakOrder (const void *va, const void *vb)  /* by position */
{ return ((const PEAK *)va)->apex - ((const PEAK *)vb)->apex ; }

static int pkFind (Array nodes, int c)
{
  int r = c, nx ;
  while (arrp (nodes, r, PKNODE)->parent != r) r = arrp (nodes, r, PKNODE)->parent ;
  while (arrp (nodes, c, PKNODE)->parent != r)
    { nx = arrp (nodes, c, PKNODE)->parent ;
      arrp (nodes, c, PKNODE)->parent = r ;
      c = nx ;
    }
  return r ;
}

/*----------------------------------------------------------------------*/
/* Sublevel set merge tree of fs restricted to [i0,i1].
 *
 * Sweep the samples from the highest down.  Each new sample either starts
 * a branch, extends one, or joins two: in that last case the younger of
 * the two branches dies right there and its prominence is fixed once and
 * for all.  This is 0-dimensional persistent homology, and in 1-D it is
 * exactly topographic prominence.  O(n log n), the log coming from the
 * sort alone.
 *
 * New nodes are appended to nodes[], their indices returned in the
 * caller's array; label[] is scratch of length >= i1+1.
 */

static void pkMergeTree (Array fs, int i0, int i1, Array nodes, Array label)
{
  int i, j, n = i1 - i0 + 1, last = i0, c0 = arrayMax (nodes) ;
  Array ord = arrayCreate (n, PKXV) ;

  for (i = i0 ; i <= i1 ; i++)
    {
      PKXV *p = arrayp (ord, i - i0, PKXV) ;
      p->x = i ; p->v = arr (fs, i, double) ;
      array (label, i, int) = -1 ;
    }
  arraySort (ord, pkXVOrder) ;

  for (j = 0 ; j < n ; j++)
    {
      int l, r ;
      i = arrp (ord, j, PKXV)->x ;
      last = i ;
      l = (i > i0 && arr (label, i-1, int) >= 0) ? pkFind (nodes, arr (label, i-1, int)) : -1 ;
      r = (i < i1 && arr (label, i+1, int) >= 0) ? pkFind (nodes, arr (label, i+1, int)) : -1 ;

      if (l < 0 && r < 0)                     /* birth of a new maximum */
	{
	  int c = arrayMax (nodes) ;
	  PKNODE *nd = arrayp (nodes, c, PKNODE) ;
	  nd->apex = i ; nd->saddle = -1 ;
	  nd->birth = nd->death = arr (fs, i, double) ;
	  nd->parent = c ; nd->isRoot = FALSE ; nd->kept = FALSE ;
	  array (label, i, int) = c ;
	}
      else if (r < 0)      array (label, i, int) = l ;
      else if (l < 0)      array (label, i, int) = r ;
      else if (l == r)     array (label, i, int) = l ;
      else                                    /* two branches meet here */
	{
	  int surv, dying ;
	  PKNODE *nd ;
	  if (arrp (nodes, l, PKNODE)->birth >= arrp (nodes, r, PKNODE)->birth)
	    { surv = l ; dying = r ; }
	  else
	    { surv = r ; dying = l ; }
	  nd = arrp (nodes, dying, PKNODE) ;
	  nd->death = arr (fs, i, double) ;
	  nd->saddle = i ;
	  nd->parent = surv ;
	  array (label, i, int) = surv ;
	}
    }

  /* the one branch that never died owns the whole component */
  if ((int) arrayMax (nodes) > c0)
    {
      int root = pkFind (nodes, arr (label, last, int)) ;
      PKNODE *nd = arrp (nodes, root, PKNODE) ;
      nd->death = arr (fs, last, double) ;   /* lowest sample of the region */
      nd->saddle = last ;
      nd->isRoot = TRUE ;
    }
  arrayDestroy (ord) ;
}

/*======================================================================*/
/*========================= the half-mean window =======================*/
/*======================================================================*/

/* Local baseline of one peak: the straight line joining the two valleys
 * that bound its domain.  Isolated peaks get a nearly flat line close to
 * the global background; a shoulder sitting on the flank of a big
 * neighbour gets a steep line, so that its area is not inflated by the
 * neighbour.
 */

typedef struct {
  int    vL, vR ;      /* the two valleys                                */
  double bL, bR ;      /* baseline value there                           */
} PKBASE ;

static double pkBaseAt (const PKBASE *b, double x)
{
  if (b->vR == b->vL) return b->bL ;
  return b->bL + (b->bR - b->bL) * (x - b->vL) / (double) (b->vR - b->vL) ;
}

/* Integral of  f - baseline  over [x3,x4], inclusive.
 * The f part is O(1) out of the prefix sum, the baseline part is the
 * closed form of an arithmetic progression.
 */

static double pkArea (Array sum, const PKBASE *b, int x3, int x4)
{
  double L = x4 - x3 + 1 ;
  double sf = arr (sum, x4+1, double) - arr (sum, x3, double) ;
  double sb ;

  if (b->vR == b->vL)
    sb = L * b->bL ;
  else
    sb = L * b->bL + (b->bR - b->bL) / (double) (b->vR - b->vL)
                     * (L * 0.5 * (x3 + x4) - L * b->vL) ;
  return sf - sb ;
}

/*----------------------------------------------------------------------*/
/* Fixed point of   g(x3) = g(x4) = mean(g on [x3,x4]) / 2
 *
 * Start from the whole domain and let the two ends walk in.  The endpoint
 * test is done on the smoothed profile so that noise does not chop the
 * window short, the mean is taken from the raw integral so that the area
 * stays unbiased.  Two-cycles are possible on discrete data, they are
 * detected and stopped.
 */

static void pkHalfMeanWindow (Array fs, Array sum, const PKBASE *b,
			      int apex, int x1, int x2,
			      int *x3p, int *x4p, double *areap,
			      double *xLp, double *xRp)
{
  int it, x3 = x1, x4 = x2, p3 = -1, p4 = -1 ;
  double area = 0, target = 0, gA, gB ;

  for (it = 0 ; it < 40 ; it++)
    {
      int n3, n4 ;
      area = pkArea (sum, b, x3, x4) ;
      if (area <= 0) { x3 = x4 = apex ; area = pkArea (sum, b, x3, x4) ; break ; }
      target = 0.5 * area / (x4 - x3 + 1) ;

      n3 = apex ;
      while (n3 > x1 && arr (fs, n3-1, double) - pkBaseAt (b, n3-1) >= target) n3-- ;
      n4 = apex ;
      while (n4 < x2 && arr (fs, n4+1, double) - pkBaseAt (b, n4+1) >= target) n4++ ;

      if (n3 == x3 && n4 == x4) break ;        /* fixed point */
      if (n3 == p3 && n4 == p4) { x3 = n3 ; x4 = n4 ; break ; }  /* 2-cycle */
      p3 = x3 ; p4 = x4 ;
      x3 = n3 ; x4 = n4 ;
    }
  area = pkArea (sum, b, x3, x4) ;
  target = 0.5 * area / (x4 - x3 + 1) ;

  /* refine both ends to sub-sample precision by linear interpolation of
   * the crossing g = target between the last sample in and the first out
   */
  *xLp = x3 - 0.5 ;
  if (x3 > x1)
    {
      gA = arr (fs, x3,   double) - pkBaseAt (b, x3) ;
      gB = arr (fs, x3-1, double) - pkBaseAt (b, x3-1) ;
      if (gA > gB) *xLp = x3 - (gA - target) / (gA - gB) ;
    }
  *xRp = x4 + 0.5 ;
  if (x4 < x2)
    {
      gA = arr (fs, x4,   double) - pkBaseAt (b, x4) ;
      gB = arr (fs, x4+1, double) - pkBaseAt (b, x4+1) ;
      if (gA > gB) *xRp = x4 + (gA - target) / (gA - gB) ;
    }
  if (*xRp < *xLp) { *xLp = x3 - 0.5 ; *xRp = x4 + 0.5 ; }

  *x3p = x3 ; *x4p = x4 ; *areap = area ;
}

/*----------------------------------------------------------------------*/
/* Local noise level.  See PEAKPARAMS.poissonSigma. */

static double pkSigmaAt (double sigma, Array bg, double bMed, int i, BOOL poisson)
{
  double b ;
  if (! poisson) return sigma ;
  b = arr (bg, i, double) ;
  if (b < 1) b = 1 ;
  return sigma * sqrt (b / bMed) ;
}

/*======================================================================*/
/*=============================== driver ===============================*/
/*======================================================================*/

void peakParamsDefault (PEAKPARAMS *pp)
{
  if (!pp) return ;
  memset (pp, 0, sizeof (PEAKPARAMS)) ;
  pp->baseWindow       = 0 ;      /* auto */
  pp->smooth           = 1 ;
  pp->sigma            = 0 ;      /* auto */
  pp->nSigmaThreshold  = 3.0 ;
  pp->nSigmaProminence = 3.0 ;
  pp->fracProminence   = 0.20 ;
  pp->noiseScale       = 0 ;      /* auto */
  pp->poissonSigma     = FALSE ;
  pp->minWidth         = 1 ;
  pp->minSnr           = 0 ;
}

/*----------------------------------------------------------------------*/
/* Walk down the flank away from the thresholded component, so that the
 * window is free to open below the detection threshold, which was only
 * ever a device to find the peak, not to measure it.  We stop at the
 * first local minimum, at the background, or at the neighbour.
 */

static int pkExtendLeft (Array fs, Array bg, int i0, int limit)
{
  int j = i0 ;
  while (j > limit
	 && arr (fs, j-1, double) <= arr (fs, j, double)
	 && arr (fs, j, double) > arr (bg, j, double))
    j-- ;
  return j ;
}

static int pkExtendRight (Array fs, Array bg, int i1, int limit)
{
  int j = i1 ;
  while (j < limit
	 && arr (fs, j+1, double) <= arr (fs, j, double)
	 && arr (fs, j, double) > arr (bg, j, double))
    j++ ;
  return j ;
}

/*----------------------------------------------------------------------*/

Array findPeaksFull (Array aa, PEAKPARAMS *pp,
		     Array *backgroundp, double *sigmap, AC_HANDLE h)
{
  PEAKPARAMS par ;
  Array peaks = arrayHandleCreate (32, PEAK, h) ;
  Array fd = 0, fs = 0, bg = 0, sum = 0, nodes = 0, label = 0, keep = 0 ;
  int i, n, i0, i1, w, L ;
  double bMed = 1 ;
  double sigma ;

  if (!aa || !arrayMax (aa)) { if (sigmap) *sigmap = 0 ; return peaks ; }
  n = arrayMax (aa) ;

  if (pp) par = *pp ; else peakParamsDefault (&par) ;

  /*--- 0. raw profile as doubles, and its prefix sum -----------------*/
  fd  = arrayCreate (n, double) ;
  sum = arrayCreate (n + 1, double) ;
  array (sum, 0, double) = 0 ;
  for (i = 0 ; i < n ; i++)
    {
      array (fd, i, double) = arr (aa, i, unsigned int) ;
      array (sum, i+1, double) = arr (sum, i, double) + arr (fd, i, double) ;
    }

  /*--- 1. slow background --------------------------------------------*/
  w = par.baseWindow > 2 ? par.baseWindow : (n / 16 > 9 ? n / 16 : 9) ;
  bg = arrayHandleCreate (n, double, backgroundp ? h : 0) ;
  pkBackground (fd, bg, w) ;

  /*--- 2. noise, measured at the scale of the peaks -------------------*/
  L = par.noiseScale > 1 ? par.noiseScale : (w / 64 > 4 ? w / 64 : 4) ;
  if (L > 256) L = 256 ;
  sigma = par.sigma > 0 ? par.sigma : pkNoiseSigma (fd, bg, L) ;
  if (sigmap) *sigmap = sigma ;

  /* reference background level, for the local sigma; subsampled, since
   * only its median is needed and the profile may be very long
   */
  {
    Array sub = arrayCreate (n / 64 + 2, double) ;
    int j = 0 ;
    for (i = 0 ; i < n ; i += 64) array (sub, j++, double) = arr (bg, i, double) ;
    bMed = pkMedian (sub) ;
    if (bMed < 1) bMed = 1 ;
    arrayDestroy (sub) ;
  }

  /*--- 3. detection copy ---------------------------------------------*/
  fs = arrayCreate (n, double) ;
  if (par.smooth > 0) pkMovingAverage (fd, fs, 2 * par.smooth + 1) ;
  else for (i = 0 ; i < n ; i++) array (fs, i, double) = arr (fd, i, double) ;

  nodes = arrayCreate (64, PKNODE) ;
  label = arrayCreate (n, int) ;
  keep  = arrayCreate (32, int) ;

  /*--- 4..8 one connected component at a time ------------------------*/
  i = 0 ;
  while (i < n)
    {
      int c0, c, nKept, k, e0, e1, prev ;

      while (i < n &&
	     arr (fs, i, double) <= arr (bg, i, double) + par.nSigmaThreshold
	     * pkSigmaAt (sigma, bg, bMed, i, par.poissonSigma))
	i++ ;
      if (i >= n) break ;
      i0 = i ;
      while (i < n &&
	     arr (fs, i, double) >  arr (bg, i, double) + par.nSigmaThreshold
	     * pkSigmaAt (sigma, bg, bMed, i, par.poissonSigma))
	i++ ;
      i1 = i - 1 ;

      /* open the region downhill on both sides, below the threshold */
      e0 = pkExtendLeft  (fs, bg, i0, 0) ;
      e1 = pkExtendRight (fs, bg, i1, n - 1) ;

      /*--- 5. merge tree of the extended region ----------------------*/
      c0 = arrayMax (nodes) ;
      pkMergeTree (fs, e0, e1, nodes, label) ;

      /*--- 6. which local maxima are real peaks ----------------------*/
      keep = arrayReCreate (keep, 32, int) ;
      for (c = c0 ; c < (int) arrayMax (nodes) ; c++)
	{
	  PKNODE *nd = arrp (nodes, c, PKNODE) ;
	  double prom   = nd->birth - nd->death ;
	  double height = nd->birth - arr (bg, nd->apex, double) ;
	  double need   = par.nSigmaProminence
	    * pkSigmaAt (sigma, bg, bMed, nd->apex, par.poissonSigma) ;
	  if (par.fracProminence * height > need) need = par.fracProminence * height ;
	  /* The root is NOT exempt: a noise blip that happens to cross the
	   * threshold is a component too, and its root must earn its place
	   * exactly like any other branch.  A component may yield 0 peaks.
	   */
	  if (prom >= need)
	    { nd->kept = TRUE ; array (keep, arrayMax (keep), int) = c ; }
	}
      nKept = arrayMax (keep) ;
      if (!nKept) continue ;
      arraySort (keep, pkIntOrder) ;   /* by node index, not yet by position */
      {
	Array ap = arrayCreate (nKept, PKXV) ;
	for (k = 0 ; k < nKept ; k++)
	  { PKXV *p = arrayp (ap, k, PKXV) ;
	    p->x = arr (keep, k, int) ;
	    p->v = - arrp (nodes, p->x, PKNODE)->apex ;   /* ascending apex */
	  }
	arraySort (ap, pkXVOrder) ;
	for (k = 0 ; k < nKept ; k++) array (keep, k, int) = arrp (ap, k, PKXV)->x ;
	arrayDestroy (ap) ;
      }

      /*--- 7. cut the region at the deepest point between apices -----*/
      prev = e0 ;
      for (k = 0 ; k < nKept ; k++)
	{
	  PKNODE *nd = arrp (nodes, arr (keep, k, int), PKNODE) ;
	  int x1 = prev, x2, j, jmin ;
	  int apex = nd->apex ;

	  if (k + 1 < nKept)
	    {
	      int nextApex = arrp (nodes, arr (keep, k+1, int), PKNODE)->apex ;
	      jmin = apex ;
	      for (j = apex ; j <= nextApex ; j++)
		if (arr (fs, j, double) < arr (fs, jmin, double)) jmin = j ;
	      x2 = jmin ;
	    }
	  else
	    x2 = e1 ;

	  /*--- 8. local baseline, then the half-mean window ----------*/
	  {
	    PKBASE b ;
	    PEAK *pk ;
	    int x3, x4 ;
	    double area, xL, xR ;

	    b.vL = x1 ; b.bL = arr (fs, x1, double) ;
	    b.vR = x2 ; b.bR = arr (fs, x2, double) ;

	    pkHalfMeanWindow (fs, sum, &b, apex, x1, x2, &x3, &x4, &area, &xL, &xR) ;

	    if (x4 - x3 + 1 < par.minWidth) { prev = x2 ; continue ; }

	    pk = arrayp (peaks, arrayMax (peaks), PEAK) ;
	    memset (pk, 0, sizeof (PEAK)) ;
	    pk->apex  = apex ;
	    pk->fApex = arr (aa, apex, unsigned int) ;
	    pk->x1 = x1 ; pk->x2 = x2 ;
	    pk->x3 = x3 ; pk->x4 = x4 ;
	    pk->xLeft = xL ; pk->xRight = xR ;
	    pk->width = xR - xL ;
	    pk->base   = pkBaseAt (&b, apex) ;
	    pk->height = arr (fd, apex, double) - pk->base ;
	    pk->prominence = nd->birth - nd->death ;
	    pk->area = area ;
	    pk->totalArea = pkArea (sum, &b, x1, x2) ;
	    pk->snr = area
	      / (pkSigmaAt (sigma, bg, bMed, apex, par.poissonSigma)
		 * sqrt (pk->width > 0 ? pk->width : 1)) ;
	    pk->blendedLeft  = (x3 <= x1 && k > 0) ;
	    pk->blendedRight = (x4 >= x2 && k + 1 < nKept) ;

	    if (par.minSnr > 0 && pk->snr < par.minSnr)
	      arrayMax (peaks) -= 1 ;
	  }
	  prev = x2 ;
	}
    }

  arraySort (peaks, pkPeakOrder) ;

  if (backgroundp) *backgroundp = bg ; else arrayDestroy (bg) ;
  arrayDestroy (fd) ;
  arrayDestroy (fs) ;
  arrayDestroy (sum) ;
  arrayDestroy (nodes) ;
  arrayDestroy (label) ;
  arrayDestroy (keep) ;

  return peaks ;
}

/*----------------------------------------------------------------------*/

Array findPeaksWithParams (Array aa, PEAKPARAMS *pp, AC_HANDLE h)
{ return findPeaksFull (aa, pp, 0, 0, h) ; }

Array findPeaks (Array aa, AC_HANDLE h)
{ return findPeaksFull (aa, 0, 0, 0, h) ; }

/*======================================================================*/
/*=================== pooled cross-chromosome noise ====================*/
/*======================================================================*/
/*----------------------------------------------------------------------*/

static void peakShow (ACEOUT ao, Array peaks, const char *target, int step, int posMin, int minCover)
{
  int i ;
  aceOutf (ao, "# %d peak%s\n", peaks ? (int) arrayMax (peaks) : 0,
	   peaks && arrayMax (peaks) == 1 ? "" : "s") ;
  aceOutf (ao, "# target\tx1\tx2\twidth\tmax cover\tmax position\taverage cover\taligned bp\tbase\tpromin\tsnr\n") ;
  for (i = 0 ; peaks && i < arrayMax (peaks) ; i++)
    {
      PEAK *p = arrp (peaks, i, PEAK) ;
      float av = (float)p->area / ((p->x4 - p->x3 + 1) * step) ;
      if (av < 30)
	continue ;
      aceOutf (ao, "%s\t%d\t%d\t%d"
	       , target
	       , p->x3  * step + posMin - step/2, p->x4  * step + posMin + step/2
	       , (p->x4 - p->x3 + 1) * step
	       ) ;
      aceOutf (ao, "\t%.1f\t%.1f\t%.1f\t%.1f"
	       , (float)p->fApex / step 
	       , (float)p->apex * step + posMin
	       , av
	       , (float)p->area 
	       ) ;
      aceOutf (ao, "\t%.1f\t%.1f\t%.1f\n" 
	       , p->base/step, p->prominence/step, p->snr
	       ) ;
    }
} /* peakShow */

/**************************** End of File ******************************/

void peaksCreateExport (ACEOUT ao, const char *target, int posMin, int step, int minCover, Array cc)
{
  AC_HANDLE h = ac_new_handle () ;
  PEAKPARAMS pp ;
  Array peaks = 0 ;
  peakParamsDefault (&pp) ;
  pp.baseWindow = 30000 / step ;
  pp.smooth     = (step >= 10) ? 0 : 1 ;
  pp.minWidth   = 3 ;
  pp.minCover = minCover ;
  pp.minSnr     = 0 ;            /* filter afterwards, not here */
	
  peaks = findPeaksWithParams (cc, &pp, h) ;

  fprintf (stderr, "// %s : %d bins, sigma %.2f, baseWindow %d bins, minCover %d, %d peaks\n",
           "chrom", arrayMax (cc), pp.sigma, pp.baseWindow, pp.minCover, arrayMax (peaks)) ;

  peakShow (ao, peaks, target, step, posMin, minCover) ;

  ac_free (h) ;   /* frees peaks and bg, both allocated on h */
} /* sxWiggleExportMultiPeaks */

