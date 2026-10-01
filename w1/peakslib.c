
/*  File: peaks.c
 *  Implementation of peaks.h, see that file for the method.
 */

#include "ac.h"
#include "peaks.h"
#include "wiggle.h"

typedef struct peakStruct {
  int x1, x2 ;   /* array coordinates */
  int yMax ;
  long int area, level ;
} PEAK ;

/**************************************************************/
/**************************************************************/
#define CSTMAGIC 5441927

static void cisTransCheck (CST cst, const char *caller)
{
  if (! cst)
    messcrash ("%s called on null cst", caller) ;
  if (cst->magic == 0)
    messcrash ("%s freed cst", caller) ;
  if (cst->magic != CSTMAGIC)
    messcrash ("%s called on bad SCR", caller) ;
} /* cisTransCheck */

/**************************************************************/

static void cisTransFinalize (void *vp)
{
  CST cst = (CST) vp ;
  cisTransCheck (cst, "cisTransFinalize") ;
  cst->magic = 0 ;
  free (cst->uu) ;
  free (cst->uv) ;
} /* scratchFinalize */

/**************************************************************/
/* Finally compute the best shift between the 2 strands */
void cisTransNormalize (CST cst)
{
  cisTransCheck (cst, "cisTransNormalize") ;
  int nu = cst->nu ;
  if (cst->nu > 0)
    {
      double u1 = cst->u1, v1 = cst->v1 ;
      cst->u2 -= u1 * u1 / nu ;
      cst->v2 -= v1 * v1 / nu ;
      double z = sqrt (cst->u2 * cst->v2) ;
      for (int dx = 0 ; dx < cst->dxMax ; dx++)
	{
	  cst->uv[dx] -= u1 * v1 / nu ; cst->uv[dx] /= z ;
	  cst->uu[dx] -= u1 * u1 / nu ; cst->uu[dx] /= cst->u2 ;
	  if (cst->bestUv < cst->uv[dx])
	    {
	      cst->bestUv = cst->uv[dx] ;
	      cst->bestShift = dx ;
	    }
	}
      cst->u1 /= nu ;
      cst->v1 /= nu ;
      cst->u2 /= nu ;
      cst->v2 /= nu ;
      cst->nu = -nu ;  /* to be sure we do not renormalize twice */
    }
  return ;
} /* cisTransNormalize */

/**************************************************************/
/* Accumulate the counts chromosome per chromosome  */
void cisTransCumulate (CST cst, int step, BigArray aaf, BigArray aar)
{
  long int ii, iMax = aaf ? bigArrayMax (aaf) : 0 ; 
  long int jj, jMax = aar ? bigArrayMax (aar) : 0 ;

  /* make room */
  WIGGLEPOINT *z1p = bigArrp(aaf, 0, WIGGLEPOINT) ;
  WIGGLEPOINT *z2p = bigArrp(aar, 0, WIGGLEPOINT) ;

  cisTransCheck (cst, "cisTransCumulate") ;
  ii = jj = 0 ;

  int nu = 0, dx0 = (z2p->x - z1p->x)/step ;
  double z, u1 = 0, v1 = 0, u2 = 0, v2 = 0 ;
  double *uu = cst->uu, *uv = cst->uv ;
  if (dx0 > 0) 
    { z1p += dx0 ; ii += dx0 ; }
  if (dx0 < 0) 
    { z2p -= dx0 ; jj -= dx0 ; }
  for (; ii < iMax  && jj < jMax ; z1p++, z2p++, ii++, jj++)
    { 
      z = z1p->y ; u2 += z*z ; u1 += z ;  nu++ ;
      if (z > 0.0001) 
	for (int dx = 0 ; dx < cst->dxMax ; dx++)
	  {
	    if (ii + dx < iMax) uu[dx] += z * ((z1p+dx)->y) ;
	    if (jj + dx < jMax) uv[dx] += z * ((z2p+dx)->y) ;
	  }
      z = z2p->y ; v2 += z*z ; v1 += z ;
    }
  cst->nu += nu ;
  cst->u1 += u1 ;
  cst->v1 += v1 ;
  cst->u2 += u2 ;
  cst->v2 += v2 ;

  return ;
} /* cisTransCumulate */

/**************************************************************/

CST cisTransCreate (int dxMax, AC_HANDLE h)
{
  int n1 = sizeof (struct cstStruct) ;
  int nn = dxMax * sizeof(double) ;
  CST cst = (CST) handleAlloc (cisTransFinalize, h, n1) ;
  memset (cst, 0, n1) ;
  
  cst->magic = CSTMAGIC ;
  cst->dxMax = dxMax ;
  cst->uu = malloc (nn) ;
  cst->uv = malloc (nn) ;
  memset (cst->uu, 0, nn) ;
  memset (cst->uv, 0, nn) ;

  return cst ;
} /* cisTransCreate */

/**************************************************************/
/**************************************************************/
/*
  bin/wiggle -f tmp/SA/SRR3740166/wiggles/SRR3740166.NC_050103.1.u.fr -I AZ -multiPeaks 2 -O COUNT -o tmp/Peaks/SRR3740166/SRR3740166.NC_050103.1.u.fr
*/
static void peakFind (Array aa, Array peaks, int minCover) 
{
  int iMax = arrayMax (aa), iPeak = 0 ;
  double z1 = 1.5 ;
  unsigned int *yp ;
  for (int i = 0 ; i < iMax ; i++)   /* scan whole wiggle */
    {
      int j, x1 = i, x2 = i, dx ;
      int yMax = 0 ;
      long int area = 0 ;
      long int areaBelow = 0 ;
      
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
		  unsigned int median, unsigned int medianNoZero,
		  int scale
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
	  float y = pk->area / (scale * ln) ;
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
	  int ln = (pk->x2 - pk->x1 + 1) ;
	  aceOutf (ao, "%s\t%d\t%d\t%d\t%.1f\t%.0f\t%.1f\t%d\n",
		   target,
		   step * pk->x1 + posMin - step/2,
		   step * pk->x2 + posMin + step/2,
		   step * ln, 
		   (float)pk->yMax / (step * scale),  /* yMax = nb of bases in one bin */
		   (float)pk->area / (ln * step * scale), (float) pk->area/(1000.0 * scale), pk->level/step
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
			int minCover,   // default 3 median (aa) no zero
			int scale      // divide aa values by scale
			)
{
  AC_HANDLE h = ac_new_handle () ;
  Array peaks = arrayHandleCreate (0x1 << 15, PEAK, h) ;  
  unsigned int median = arrayUintMedian (aa, FALSE) ;
  unsigned int medianNoZero = arrayUintMedian (aa, TRUE) ;
  if (minCover <= 0)
    {
      minCover = 5 * medianNoZero ; // was 3
      if (minCover < 5) minCover = 5 ;
    }
  if (scale < 1) scale = 1 ;
  peakFind (aa, peaks, scale * minCover) ;
  
  fprintf (stderr, "// %s : found %d peaks at minCover %d  medianNoZero %.2f\n",
	   target, arrayMax (peaks), minCover, (float)medianNoZero/scale) ;

  peaksExport (ao, aoLevels, peaks, minCover, target, posMin, step, median, medianNoZero, scale) ;

  ac_free (h) ;   /* frees peaks */
} /* peaksCreateExport */

/***********************************************************************/
/**************************** End of File ******************************/
/***********************************************************************/
