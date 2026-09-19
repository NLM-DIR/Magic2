/*  File: peaks.c
 *  Implementation of peaks.h, see that file for the method.
 */

#include "ac.h"
#include "peaks.h"

typedef struct peakStruct {
  int x1, x2 ;   /* array coordinates */
  int yMax, area, level ;
} PEAK ;

/***********************************************************************/
/*
  bin/wiggle -f tmp/SA/SRR3740166/wiggles/SRR3740166.NC_050103.1.u.fr -I AZ -multiPeaks 2 -O COUNT -o tmp/Peaks/SRR3740166/SRR3740166.NC_050103.1.u.fr
*/
static void peakFind (Array aa, Array peaks, int minCover) 
{
  int iMax = arrayMax (aa), iPeak = 0 ;
  double z1 = 1.5 ;
  unsigned int *yp ;

  if (minCover < 5) minCover = 5 ;
  for (int i = 0 ; i < iMax ; i++)   /* scan whole wiggle */
    {
      int j, x1 = i, x2 = i, dx ;
      int yMax = 0 ;
      int area = 0 ;
      int areaBelow = 0 ;
      
      yp = arrp (aa, i, unsigned int) ;      /* passing over minCover */
      for (j = i ; j < iMax && *yp < minCover ; yp++, j++)
	areaBelow += *yp ;
      x1 = x2 = j ; area = 0 ;
      
      for ( ; j < iMax && *yp >= minCover ; yp++, j++)
	{    /* untill we fall below */
	  yMax = (*yp > yMax ? *yp : yMax) ;
	  area += *yp ;
	}
      x2 = j - 1 ;
      dx = x2 - x1 + 1 ;

      if (area > z1 * dx * minCover)
	{    /* register */
	  PEAK *pk = iPeak ? arrayp (peaks, iPeak - 1, PEAK) : 0 ;
	  if (pk && areaBelow > .8 * (x1 - pk->x2) * minCover &&
	      (area + pk->area + areaBelow) > z1 * (x2 - pk->x1) * minCover
	      )
	    { /* extend previous peak */
	      pk->x2 = x2 ;
	      pk->yMax = (yMax > pk->yMax ? yMax : pk->yMax) ;
	      pk->area += area + areaBelow ;
	    }
	  else
	    {
	      pk = arrayp (peaks, iPeak++, PEAK) ;
	      pk->x1 = x1 ;
	      pk->x2 = x2 ;
	      pk->yMax = yMax ;
	      pk->area  = area ;
	    }
	}
      i = j - 1 ;
    }
} /* peakFind */

/**********************************************************/

void peaksExport (ACEOUT ao, ACEOUT aoLevels,
		  Array peaks, int minCover, 
		  const char *target, int posMin, int step,
		  unsigned int median, unsigned int medianNoZero
		  )
{
  int i, n = peaks ? arrayMax (peaks) : 0 ;
  KEYSET histo = 0 ;
  
  if (aoLevels)
    {
      histo = keySetCreate () ;
      for (i = 0 ; i < n ; i++)
	{
	  PEAK *pk = arrp (peaks, i, PEAK) ;
	  int ln = pk->x2 - pk->x1 + 1 ;
	  float y = pk->area / ln ;
	  int k = utMainPart (y) ;
	  int k1 = 0 ;
	  int k2 = 10 ;
	  int k0 = k ;
	  while (k0 >= k2) {k1+=3; k0 /= 10 ;}
	  k1 += (k0 >= 2 ? 1 : 0) + (k0 >= 5 ? 1 : 0) ;
	  keySet (histo, k1) += 1 ;
	  pk->level = k ;
	}
    }
  if (ao)
    {
      aceOutf (ao, "# %d distinct peak shape%s threshold\t%.1f\n", n, n == 1 ? "" : "s", (float)minCover/step) ;
      aceOutf (ao, "# target\tx1\tx2\twidth\theight\tavg\tarea kb\tlevel\n") ;
      for (i = 0 ; i < n ; i++)
	{
	  PEAK *pk = arrp (peaks, i, PEAK) ;
	  int ln = step * (pk->x2 - pk->x1 + 1) ;
	  aceOutf (ao, "%s\t%d\t%d\t%d\t%.1f\t%.0f\t%.1f\t%d\n",
		   target,
		   pk->x1 * step + posMin - step/2,
		   pk->x2 * step + posMin + step/2,
		   ln, 
		   (float)pk->yMax / step,  /* yMax = nb of bases in one bin */
		   (float)pk->area / ln, (float)pk->area/1000.0, pk->level/step
		   ) ;
	}
    }

  if (aoLevels)
    {
      aceOutf (aoLevels, "# target\tminCover\tnPeaks\t median=%u medianNoZero=%u\n", median, medianNoZero) ;
      int k1 = 1, dk = 0 ;
      for (i = 0 ; i < keySetMax (histo) ; i++, dk = (dk + 1) % 3)
	{
	  aceOutf (aoLevels, "%s\t%d\t%d\n", target, k1, keySet (histo,i)) ;
	  switch (dk)
	    {
	    case 0:
	    case 2:
	      k1 *= 2 ;
	      break ;
	    case 1:
	      k1 /= 2 ; k1 *= 5 ;
	      break ;
	    }
	}
    }

  keySetDestroy (histo) ;
  /* levelCount and peaks are owned by whichever AC_HANDLE allocated them
   * in peakCaller (see peaksCreateExport below) -- not destroyed here,
   * or ac_free(h) would double-free them
   */
} /* peaksExport */

/***********************************************************************/

void peaksCreateExport (ACEOUT ao, ACEOUT aoLevels,
			Array aa ,    // array of unsigned int
			const char *target,
			int posMin, int step, // x = i * step + posMin
			int minCover   // default 3 median (aa) no zero
			)
{
  AC_HANDLE h = ac_new_handle () ;
  Array peaks = arrayHandleCreate (0x1 << 15, PEAK, h) ;  
  unsigned int median = arrayUintMedian (aa, FALSE) ;
  unsigned int medianNoZero = arrayUintMedian (aa, TRUE) ;
  minCover = 3 * medianNoZero ;

  peakFind (aa, peaks, minCover) ;
  
  fprintf (stderr, "// %s : found %d peaks at minCover %d\n",
	   target, arrayMax (peaks), minCover) ;

  peaksExport (ao, aoLevels, peaks, minCover, target, posMin, step, median, medianNoZero) ;

  ac_free (h) ;   /* frees peaks */
} /* peaksCreateExport */

/***********************************************************************/
/**************************** End of File ******************************/
/***********************************************************************/
