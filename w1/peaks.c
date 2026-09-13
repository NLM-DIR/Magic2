/*  File: peaks.c
 *  Implementation of peaks.h, see that file for the method.
 */

#include "peaks.h"
#include "ac.h"

/*----------------------------------------------------------------------*/
/* Prefix sum of the raw profile, for O(1) area queries.  One double
 * array of length n+1: this is the only O(n) buffer this module needs.
 */

static Array pkPrefixSum (Array aa)
{
  int i, n = arrayMax (aa) ;
  Array sum = arrayCreate (n + 1, double) ;

  array (sum, 0, double) = 0 ;
  for (i = 0 ; i < n ; i++)
    array (sum, i+1, double) = arr (sum, i, double) + arr (aa, i, unsigned int) ;
  return sum ;
}

/*----------------------------------------------------------------------*/
/* Scan [x1,x2] for maximal runs with f(x) >= cover.  A run's shape stays
 * exactly the same for every cover up to and including the minimum
 * value of f inside it: nothing breaks until cover exceeds that
 * minimum, at which point the bin(s) achieving it drop out and the run
 * must shrink or split.  That minimum is tracked for free while doing
 * the one scan that finds the run in the first place, so instead of
 * stepping cover by cover*z and re-scanning the same bins over and
 * over just to confirm nothing changed, we jump cover straight to the
 * first geometric level that exceeds the minimum.  The shape-rule test
 * (area > z*length*cover) is then just arithmetic on the fixed
 * (area,length) of the persisting run, so the intervening geometric
 * levels are filled into the level-count histogram at O(1) each, with
 * no re-scan of the bin array at all.  A real re-scan only happens
 * once the jump lands on a cover that can actually change the shape --
 * i.e. exactly when there is something genuinely new to find.
 *
 * A run is only ADDED TO THE PEAK TABLE the first time its shape
 * appears (see isNew below): re-confirming an unchanged run at higher
 * levels would just duplicate the same row.
 */

static void pkRecurse (Array aa, Array sum, int x1, int x2,
		       double cover, double z, int level,
		       Array peaks, Array levelCount)
{
  int i = x1 ;

  while (i <= x2)
    {
      int j, compEnd, length ;
      unsigned int minVal, v ;
      double area, curCover ;
      int curLevel ;
      BOOL isNew, stillPassing ;

      while (i <= x2 && arr (aa, i, unsigned int) < cover) i++ ;
      if (i > x2) break ;

      minVal = arr (aa, i, unsigned int) ;
      j = i ;
      while (j <= x2 && (v = arr (aa, j, unsigned int)) >= cover)
	{
	  if (v < minVal) minVal = v ;
	  j++ ;
	}
      compEnd = j - 1 ;
      length = compEnd - i + 1 ;
      area = arr (sum, compEnd + 1, double) - arr (sum, i, double) ;

      isNew = (level == 0) || i != x1 || compEnd != x2 ;

      /* this exact shape is valid for every cover in [cover, minVal];
       * walk the geometric levels in that range at O(1) each, with no
       * further scanning of the bins.  curCover/curLevel are always
       * advanced through to just past minVal before recursing below --
       * skipping that would call pkRecurse again with an unchanged
       * (cover,level) and recurse forever -- so only the levelCount
       * bookkeeping, not the advance itself, is skipped once the shape
       * rule stops passing (by monotonicity, fixed area/length can
       * never pass again once cover has grown past the point it fails).
       */
      curCover = cover ;
      curLevel = level ;
      stillPassing = TRUE ;
      while (curCover <= (double) minVal)
	{
	  if (stillPassing && area > z * length * curCover)
	    {
	      array (levelCount, curLevel, int) += 1 ;
	      if (isNew && curLevel == level)
		{
		  PEAK *pk = arrayp (peaks, arrayMax (peaks), PEAK) ;
		  pk->x1 = i ; pk->x2 = compEnd ;
		  pk->level = level ;
		  pk->cover = cover ;
		  pk->area  = area ;
		  pk->avg   = area / length ;
		}
	    }
	  else
	    stillPassing = FALSE ;
	  curCover *= z ;
	  curLevel++ ;
	}

      /* a real re-scan only happens once curCover exceeds minVal, i.e.
       * exactly when the shape can actually be different
       */
      pkRecurse (aa, sum, i, compEnd, curCover, z, curLevel, peaks, levelCount) ;

      i = j ;
    }
}

/*----------------------------------------------------------------------*/

Array peakCaller (Array aa, int minCover, double z, Array *levelCountp, AC_HANDLE h)
{
  Array peaks = arrayHandleCreate (256, PEAK, h) ;
  Array levelCount = arrayHandleCreate (32, int, h) ;
  Array sum ;
  int n ;

  if (!aa || ! (n = arrayMax (aa)))
    { if (levelCountp) *levelCountp = levelCount ; else arrayDestroy (levelCount) ;
      return peaks ; }
  if (minCover < 1 || z <= 1)
    messcrash ("peakCaller: need minCover >= 1 and z > 1") ;

  sum = pkPrefixSum (aa) ;
  pkRecurse (aa, sum, 0, n - 1, (double) minCover, z, 0, peaks, levelCount) ;
  arrayDestroy (sum) ;

  if (levelCountp) *levelCountp = levelCount ; else arrayDestroy (levelCount) ;
  return peaks ;
} /* peakCaller */

/*----------------------------------------------------------------------*/

void peaksExport (ACEOUT ao, ACEOUT aoLevels, const char *target,
		  int minCover, int posMin, int step, double z, Array peaks, Array levelCount)
{
  int i, n = peaks ? arrayMax (peaks) : 0 ;

  if (ao)
    {
      aceOutf (ao, "# %d distinct peak shape%s\n", n, n == 1 ? "" : "s") ;
      aceOutf (ao, "# target\tx1\tx2\twidth\tlevel\tcover\tarea\tavg\n") ;
      for (i = 0 ; i < n ; i++)
	{
	  PEAK *p = arrp (peaks, i, PEAK) ;

	  aceOutf (ao, "%s\t%d\t%d\t%d\t%d\t%.1f\t%.1f\t%.1f\n",
		   target,
		   p->x1 * step + posMin, p->x2 * step + posMin,
		   (p->x2 - p->x1 + 1) * step,
		   p->level, p->cover, p->area * step, p->avg) ;
	}
    }

  if (aoLevels && levelCount)
    {
      aceOutf (aoLevels, "# target\tminCover\tnPeaks\n") ;
      for (i = 0 ; i < arrayMax (levelCount) ; i++)
	{
	  aceOutf (aoLevels, "%s\t%d\t%d\n", target, minCover, arr (levelCount, i, int)) ;
	  minCover *= z ;
	}
    }
  /* levelCount and peaks are owned by whichever AC_HANDLE allocated them
   * in peakCaller (see peaksCreateExport below) -- not destroyed here,
   * or ac_free(h) would double-free them
   */
} /* peaksExport */

/*----------------------------------------------------------------------*/

void peaksCreateExport (ACEOUT ao, ACEOUT aoLevels, const char *target, int posMin, int step,
			int minCover, double z, Array aa)
{
  AC_HANDLE h = ac_new_handle () ;
  Array levelCount = 0 ;
  Array peaks = peakCaller (aa, minCover, z, &levelCount, h) ;

  fprintf (stderr, "// %s : %d bins, minCover %d, z %.2f, %d distinct peak shapes\n",
	   target, aa ? arrayMax (aa) : 0, minCover, z,
	   (int) arrayMax (peaks)) ;

  peaksExport (ao, aoLevels, target, minCover, posMin, step, z, peaks, levelCount) ;

  ac_free (h) ;   /* frees peaks */
} /* peaksCreateExport */

/**************************** End of File ******************************/
