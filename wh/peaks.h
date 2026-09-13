/*  File: peaks.h
 *  Nested, multi-level threshold peak caller for a 1-D integer profile.
 *
 *  Model: f(x) >= 0, x in [0, arrayMax(aa)[, f(x) = array(aa,x,uint).
 *
 *  Definition of a peak at level "cover":
 *      - (x1,x2) is a maximal run with f(x) >= cover for every x in it
 *        (a connected component of the thresholded set), AND
 *      - area(x1,x2) = sum_{x=x1}^{x2} f(x)  >  z * (x2-x1+1) * cover
 *        i.e. the component's own average exceeds z times the level
 *        that defined it, not just barely-clears-the-bar flat plateaux.
 *
 *  Levels are visited at  cover_k = minCover * z^k,  k = 0,1,2,...
 *  Components nest by construction: the thresholded set only shrinks
 *  as cover grows, so a level-(k+1) component is always contained in
 *  the level-k component it was found inside (whether or not that
 *  containing component itself qualified as a peak).  This means each
 *  level only needs to be searched *inside* the previous level's
 *  components, never over the whole profile again -- the recursion
 *  below exploits exactly that, so the total work stays close to
 *  linear in practice instead of linear-per-level.
 *
 *  Recursion stops on its own: once cover_k exceeds the local maximum
 *  of a component, no sub-component is found and that branch returns
 *  with no further cost.
 *
 *  A component is only ever recorded (and returned in the PEAK table)
 *  the first time its exact shape appears: as cover rises a component
 *  often persists unchanged for several levels before it finally
 *  shrinks or splits, and re-recording an unchanged interval at every
 *  intervening level would only inflate the table with duplicates of
 *  the same information, not add any.
 *
 *  No baseline, no smoothing, no noise model, no merge tree: this is
 *  the fast/cheap version, to be compared for speed and memory against
 *  the fuller prominence-based peaks.c before deciding on thresholds.
 *
 *  Author: (c) 2026, free to use.
 */

#ifndef DEF_PEAKS_H
#define DEF_PEAKS_H

#include "regular.h"   /* pulls in array.h : Array, BOOL, AC_HANDLE ... */
#include "aceio.h"

/*----------------------------------------------------------------------*/
/* One detected peak, in bin (array-index) coordinates, inclusive.       */
/*----------------------------------------------------------------------*/

typedef struct peakStruct {
  int    x1, x2 ;   /* bin coordinates of the component, inclusive       */
  int    level ;    /* 0 at minCover, 1 at minCover*z, 2 at minCover*z^2 */
  double cover ;    /* the threshold that produced this component       */
  double area ;     /* sum f(x) over [x1,x2]                            */
  double avg ;       /* area / (x2-x1+1)                                 */
} PEAK ;

#define PEAK_DEFAULT_Z 1.5   /* try 1.2 .. 3 */

/*----------------------------------------------------------------------*/
/* Core call.
 *   aa          : Array of unsigned int, the profile f
 *   minCover    : starting threshold (level 0)
 *   z           : level multiplier, cover_k = minCover * z^k ; z must be > 1
 *   levelCountp : if non-zero, receives a new Array of int (on h), one
 *                 entry per level, index = level.  Each entry counts
 *                 EVERY component that satisfied the shape rule at that
 *                 level, including ones that are an unchanged carry-over
 *                 from the level below and therefore do NOT get a fresh
 *                 row in the returned table.  This is the count you want
 *                 for choosing minCover / z; the returned table is what
 *                 you want for enumerating peaks without duplication.
 *   h           : handle the returned Arrays are allocated on, may be 0
 * Returns a new Array of PEAK, sorted by (x1 ascending, level ascending),
 * possibly of length 0, never null.
 */
Array peakCaller (Array aa, int minCover, double z, Array *levelCountp, AC_HANDLE h) ;

/*----------------------------------------------------------------------*/
/* If ao is non zero, export the full peak table,
 * if aoLevels is non-zero, the distribution of how many peaks were found at each level.
 * Coordinates are rescaled by "step" and offset by "posMin", matching
 * the original wiggle's genomic coordinates.
 */
void peaksExport (ACEOUT ao, ACEOUT aoLevels, const char *target,
		  int minCover, int posMin, int step, double z, Array peaks, Array levelCount) ;

/*----------------------------------------------------------------------*/
/* Convenience wrapper: calls peakCaller then peaksExport and frees
 * everything.  Pass aoLevels = 0 if you only want the peak table.
 */
void peaksCreateExport (ACEOUT ao, ACEOUT aoLevels, const char *target, int posMin, int step,
			int minCover, double z, Array aa) ;

#endif /* DEF_PEAKS_H */

/**************************** End of File ******************************/
