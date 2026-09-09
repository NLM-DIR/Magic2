/*  File: peaks.h
 *  Multi-scale peak detection in a 1-D integer profile.
 *
 *  Model: f(x) >= 0, x in [0, arrayMax(aa)[, f(x) = array(aa,x,int)
 *         looks like a slowly varying background + noise + a set of peaks.
 *
 *  Two decisions are kept strictly separate:
 *
 *    1) HOW MANY peaks   -> topographic prominence (persistence of the
 *                           sublevel-set merge tree).  A local maximum is
 *                           a peak if the deepest valley separating it from
 *                           any higher maximum is deep enough.
 *
 *    2) HOW WIDE each is -> matched filter with a top-hat kernel against
 *                           white noise, i.e. maximise
 *
 *                                J = A / sqrt(L),   A = integral of g on the
 *                                                   window, L = its length
 *
 *                           whose stationarity condition is the implicit
 *                           "half-mean" rule
 *
 *                                g(x3) = g(x4) = (1/2) * mean(g on [x3,x4])
 *
 *                           g being f minus the peak's own local baseline.
 *
 *  Doing 1) before 2) is essential: J alone happily merges two Gaussians
 *  up to ~8.4 sigma apart into one fat window.
 *
 *  Calibration of the half-mean width (exact, noise free):
 *      rectangle of width W        ->  W            (exact, as expected)
 *      Gaussian of width sigma     ->  2.80 sigma   (FWHM = 2.355 sigma)
 *      triangle of base 2W         ->  1.333 W
 *      Lorentzian of HWHM gamma    ->  2.78 gamma
 *
 *  Author: (c) 2026, free to use.
 */

#ifndef DEF_PEAKS_H
#define DEF_PEAKS_H

#include "regular.h"   /* pulls in array.h : Array, BOOL, AC_HANDLE ... */
#include "aceio.h"

/*----------------------------------------------------------------------*/
/* One detected peak.  Positions are sample indices, inclusive.          */
/*    x1 <= x3 <= apex <= x4 <= x2                                       */
/*----------------------------------------------------------------------*/

typedef struct peakStruct {
  int    apex ;        /* position of the local maximum                  */
  int    fApex ;       /* raw f value there                              */

  int    x1, x2 ;      /* domain of the peak, valley to valley           */
  int    x3, x4 ;      /* optimal half-mean window  (the answer)         */

  double xLeft ;       /* x3 refined to sub-sample precision             */
  double xRight ;      /* x4 refined to sub-sample precision             */
  double width ;       /* xRight - xLeft                                 */

  double base ;        /* local baseline evaluated at the apex           */
  double height ;      /* fApex - base                                   */
  double prominence ;  /* depth of the deepest separating valley         */

  double area ;        /* integral of f - base over [x3,x4]              */
  double totalArea ;   /* integral of f - base over [x1,x2]              */
  double snr ;         /* area / (sigma * sqrt(width))                   */

  BOOL   blendedLeft ; /* window clipped by a valley shared with a peak  */
  BOOL   blendedRight ;/* idem on the right; width/area are then biased  */
} PEAK ;

/*----------------------------------------------------------------------*/
/* Tuning.  peakParamsDefault() fills in sane values; every field that   */
/* is left <= 0 is auto-computed from the data.                          */
/*----------------------------------------------------------------------*/

typedef struct peakParamsStruct {
  int minCover ;           /* default 10: kills counting fluctuations and low coverage */
  int    baseWindow ;      /* Morphological opening window used to build the
			    * slow background.  THE ONE PARAMETER THAT REALLY
			    * MATTERS: it must be clearly WIDER than the widest
			    * real peak, or the background climbs into that peak
			    * and splits it in two.  Wider is safe as long as it
			    * stays narrower than the background variations.
			    * <=0 : auto = max/16, a guess, override it.      */
  int    smooth ;          /* half width of the moving average used for
			    * DETECTION only; areas are always integrated
			    * on the raw data.  0 : no smoothing           */
  double sigma ;           /* noise level.  <=0 : robust MAD estimate     */
  double nSigmaThreshold ; /* detection threshold = background + n*sigma  */
  double nSigmaProminence ;/* a valley must be this deep, in sigma units  */
  double fracProminence ;  /* ... or this fraction of the peak height,
			    * whichever is larger.  0.20 is a good default */
  int    noiseScale ;      /* the noise is measured over windows of this
			    * many bins, which should be of the order of a
			    * typical peak width: a peak of length L
			    * competes against the background fluctuation
			    * at scale L, not at scale one bin.  This
			    * matters whenever the track is smoothed or
			    * deep enough that the background undulates.
			    * <=0 : auto = baseWindow/64, clamped to
			    * [4,256]                                     */
  BOOL   poissonSigma ;    /* FALSE : one sigma for the whole profile.
			    * TRUE  : sigma follows the local background as
			    *   sigma(x) = sigma * sqrt (b(x) / median b)
			    * which is what the physics of read pileup gives.
			    * If the noise atom is a read of length R arriving
			    * as a Poisson process, and f is the aligned bases
			    * per bin, then sigma = sqrt (R * f): Poisson, but
			    * inflated by the read length because one spurious
			    * event deposits R correlated bases, not one.
			    * Turn this on when the background varies by more
			    * than a factor of a few across the profile.     */
  int    minWidth ;        /* reject peaks whose window is narrower       */
  double minSnr ;          /* reject peaks below this signal to noise     */
} PEAKPARAMS ;

void peakParamsDefault (PEAKPARAMS *pp) ;

/*----------------------------------------------------------------------*/
/* Main entry points.                                                   */
/*                                                                      */
/* aa : Array of int, the profile f.                                    */
/* h  : handle the returned Array is allocated upon, may be 0.          */
/*                                                                      */
/* Returns a new Array of PEAK, sorted by increasing apex position,      */
/* possibly of length 0, never null.                                     */
/*----------------------------------------------------------------------*/

Array findPeaks (Array aa, AC_HANDLE h) ;
Array findPeaksWithParams (Array aa, PEAKPARAMS *pp, AC_HANDLE h) ;

/*----------------------------------------------------------------------*/
/* Pooled noise estimation across several profiles.
 *
 * The noise level is a property of the sample, not of the chromosome.
 * Estimating it separately per chromosome misfires badly on the nearly
 * empty ones (chrY in a female sample, small contigs), where the median
 * absolute difference is zero and sigma collapses to its floor.
 *
 * The absolute successive differences are small non-negative integers,
 * so the pooled median is obtained exactly from a histogram, in one
 * linear scan per chromosome and constant memory.  Give each thread its
 * own histogram and merge them:
 *
 *    Array hh = peakNoiseCreate (h) ;
 *    for each chromosome, in parallel :
 *        Array hi = peakNoiseCreate (h) ;
 *        peakNoiseAccumulate (hi, aa) ;
 *    ... then, serially :
 *        peakNoiseMerge (hh, hi) ;
 *    par.sigma = peakNoiseSigma (hh) ;
 *
 * Exclude the mitochondrion and any contig whose depth is wildly out of
 * line with the rest; run those with their own sigma.
 */
Array  peakNoiseCreate (AC_HANDLE h) ;
void   peakNoiseAccumulate (Array histo, Array aa) ;
void   peakNoiseMerge (Array dest, Array src) ;
double peakNoiseSigma (Array histo) ;             /* median, q = 0.50 */

/* The median estimator assumes independent adjacent samples.  It returns
 * the floor of 1.0 when the data are quantised (a float wiggle truncated
 * to int) or pre-smoothed, because most covered adjacent pairs are then
 * exactly equal.  peakNoiseShow() diagnoses that in one line; when it
 * happens, take a high quantile instead, or set PEAKPARAMS.sigma by hand.
 */
double peakNoiseSigmaQ (Array histo, double q) ;  /* q = 0.75 or 0.90 */
void   peakNoiseShow (Array histo, FILE *fo) ;

/* Same, but also hands back the background it subtracted (Array of
 * double, same length as aa) and the noise level it used.  Pass 0 for
 * anything you do not want.  *backgroundp is allocated on h.
 */
Array findPeaksFull (Array aa, PEAKPARAMS *pp,
		     Array *backgroundp, double *sigmap, AC_HANDLE h) ;

/* Human readable dump, one line per peak. */
//static void peakShow (ACEOUT ao, Array peaks, const char *target, int step, int posMin, int minCover) ;
void peaksCreateExport (ACEOUT ao, const char *target, int posMin, int step, int minCover, Array cc) ;
#endif /* DEF_PEAKS_H */

/**************************** End of File ******************************/
