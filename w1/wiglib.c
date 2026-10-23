#include "../wsa/sa.h"
#include "wig.h"
#include "wiggle.h"
#include "fastItoA.h"
#include <zlib.h>
#include <fcntl.h>  // for O_WRONLY if using write()
#include <unistd.h> // for write()

typedef struct peakStruct {
  int x1, x2 ;   /* array coordinates */
  int yMax ;
  long int area, level ;
} WPEAK ;

/*************************************************************************************/

static inline int fast_itoa_nl(char *buf, int val)
{
    char *p = buf;

    if (val < 0) {
        *p++ = '-';
        val = -val;
    }

    if (val == 0) {
        *p++ = '0';
    } else {
        char tmp[10];
        int i = 0;
        while (val > 0) {
            tmp[i++] = '0' + val % 10;
            val /= 10;
        }
        while (i--) *p++ = tmp[i];
    }

    *p++ = '\n';
    return p - buf;           // bytes written
}

/**************************************************************/
/**************************************************************/
#define WIGMAGIC 2131927

static void wigCheck (WIG wig, const char *caller)
{
  if (! wig)
    messcrash ("%s called on null wig", caller) ;
  if (wig->magic == 0)
    messcrash ("%s freed wig", caller) ;
  if (wig->magic != WIGMAGIC)
    messcrash ("%s called on bad WIG", caller) ;
} /* wigCheck */

/**************************************************************/

static void wigFinalize (void *vp)
{
  WIG wig = (WIG ) vp ;
  wigCheck (wig, "wigFinalize") ;
  wig->magic = 0 ;
  ac_free (wig->aa) ;
  ac_free (wig->peaks) ;
} /* wigFinalize */

/**************************************************************/

WIG wigCreate (const char *runName, const char *chromName, int step, AC_HANDLE h)
{
  int n1 = sizeof (struct wigStruct) ;
  WIG wig = handleAlloc (wigFinalize, h, n1) ;
  memset (wig, 0, n1) ;
  
  wig->magic = WIGMAGIC ;
  wig->runName = runName ? runName : "xxx" ;
  wig->chromName = chromName ? chromName : "zzz" ;
  wig->step = step ? step : 1 ;
  wig->scale = 1 ;
  wig->aa = arrayCreate (64, unsigned int) ;

  return wig ;
} /* wigCreate */

/**************************************************************/
/* accumulate data in the wig */
void wigParseAZ (WIG wig, const char *fNam, int strand_shift)
{
  AC_HANDLE h = ac_new_handle () ;
  AZZ *az = wigAzOpen (fNam, h) ;
  int nn = 0 ;
  if (! az)
    messcrash ("Cannot open file: %s", fNam) ;
  if (az->posMax <= 0)
    { ac_free (h) ; return ; }
  if (wig->nBp == 0) wig->step = 0 ;
  if (wig->step && wig->step != az->step)
    messcrash ("makeGroupWiggles cannot mix wiggles with different step %d != %d\n", wig->step, az->step) ;
  wig->step = az->step ;
    
  wigAzZone (az, 0, az->posMax, strand_shift, wig, 0, &nn, &(wig->nBp), TRUE) ;
  fprintf (stderr, "// Parsed %d positions, %ld bases, in file %s\n", nn, wig->nBp, fNam) ;
  ac_free (h) ;
} /* wigParseAz */

/**************************************************************/
/**************************************************************/

int wigSetShift (WIG wig, int shift)  /* returns previous value */
{
  wigCheck (wig, "wigSetShift") ;
  int old = wig->shift ;
  wig->shift = shift ;
  return old ;
} /* wigSetShift */

/**************************************************************/
/**************************************************************/
/* the peaks are given in raw aa offsets */

int wigFindPeaks (WIG wig, int minCover)
{
  wigCheck (wig, "wigFindPeaks") ;
  Array aa = wig->aa ;
  int iMax = arrayMax (aa), iPeak = 0 ;
  int scale = wig->scale ;
  double z1 = 1.5 ;
  unsigned int *yp ;
  Array peaks = wig->peaks = arrayReCreate (wig->peaks, 1024, WPEAK) ;

  //  unsigned int median ; median = wig->median = arrayUintMedian (aa, FALSE) ;
  unsigned int medianNoZero ;  medianNoZero = wig->medianNoZero = arrayUintMedian (aa, TRUE) ;
  if (minCover <= 0)
    minCover = 10 ;
  wig->minCover = minCover ;
  int mm = medianNoZero > scale ? medianNoZero : scale ;
  minCover = minCover * mm ;


  
  if (scale < 1) scale = 1 ;

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
	  WPEAK *pk = iPeak ? arrayp (peaks, iPeak - 1, WPEAK) : 0 ;
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
	      pk = arrayp (peaks, iPeak++, WPEAK) ;
	      pk->x1 = x1 >= 0 ? x1 : 0 ;
	      pk->x2 = x2 ;
	      pk->yMax = yMax ;
	      pk->area  = area ;
	    }
	}
      i = j - 1 ;
    }
  return arrayMax (peaks) ;
} /* wigFindPeaks */

/**************************************************************/

static void wigExportPeaksDo (ACEOUT ao, WIG wig, AC_HANDLE h)
{
  Array peaks = wig->peaks ;
  int i, n = arrayMax (peaks) ;
  int step = wig->step ;
  int scale = wig->scale ? wig->scale : 1 ;
  int shift = wig->shift ;
  double totalCoverInPeaks = 0 ;
  double totalAreaInPeaks = 0 ;
  
  aceOutDate (ao, "###", hprintf (h, "Found %d Peaks in  run %s, target %s, at min height %.2f"
				  , n, wig->runName, wig->chromName, ((float)wig->minCover)/scale)
	      ) ;
  aceOutf (ao, "## Median\t%d\tMinCover\t%d\tPeaks\t%d\n"
	   , wig->medianNoZero/scale
	   , wig->minCover
	   , n
	   ) ;

  aceOutf (ao, "# target\tx1\tx2\twidth\theight\tavg\tarea kb\n") ;
  for (i = 0 ; i < n ; i++)
    {
      WPEAK *pk = arrp (peaks, i, WPEAK) ;
      int ln = (pk->x2 - pk->x1 + 1) ;
      int c1 = step * pk->x1 + shift - (step-1)/2 + wig->posMin ;
      int c2 = step * pk->x2 + shift + step/2 + wig->posMin ;
      if (c1 < 1) c1 = 1 ;
      if (c2 < c1) continue ;     /* peak entirely before position 1 */
      aceOutf (ao, "%s\t%d\t%d\t%d\t%.1f\t%.1f\t%.1f\n",
	       wig->chromName,
	       c1, c2,
	       step * ln, 
	       (float)pk->yMax / (scale),  /* yMax = nb of bases in one bin */
	       (float)pk->area / (ln * scale), (float) pk->area/(1000.0 * scale)
	       ) ;
      totalCoverInPeaks += step * ln ;
      totalAreaInPeaks += pk->area/scale ;
    }

  aceOutf (ao, "## totalCoverInPeaks\t%.3f\n"
	   , totalCoverInPeaks/scale
	   , ((float)wig->minCover)
	   , n
	   ) ;
  aceOutf (ao, "## totalAreaInPeaks\t%.3f\n"
	   , totalAreaInPeaks
	   ) ;
  aceOutf (ao, "## genomeLength\t%.3f\n"
	   , wig->chromLength
	   ) ;
  aceOutf (ao, "## totalAlignedBases\t%.3f\n"
	   , wig->bases
	   ) ;

} /* wigExportPeaksDo */

/**************************************************************/

void wigExportPeaks (WIG wig, const char *fName, BOOL gzo)
{
  AC_HANDLE h = ac_new_handle () ;

  ACEOUT ao = aceOutCreate (fName, ".peaks", gzo, h) ;
  wigCheck (wig, "wigExportPeaks") ;

  if (ao && wig && wig->peaks)
    wigExportPeaksDo (ao, wig,  h) ;
  ac_free (h) ;
  return ;
} /*  wigExportPeaks */

/**************************************************************/
/**************************************************************/
/* feed in a WIG an array of WSEG alignment segments STEP=1 */
void wigFeed (WIG wig, BigArray wsegs, int scale)
{
  wigCheck (wig, "wigFeed") ;
  Array aa = wig->aa ;
  wig->scale = scale ;
  wig->step = 1 ;
  
  int j1, jj, jMax = bigArrayMax (wsegs) ;
  WSEG *wseg ;
  int posMin = 2000000000 ; /* 2G */
  int posMax = 0 ;
  if (!jMax) return ;

  if (wig->nBp)
    messcrash ("wigFeed can only be called once, otherwise we should x-shift the data in a correct way") ;
  posMin = INT_MAX ; posMax = 0 ;
  for (j1 = 0, wseg = bigArrp (wsegs, 0, WSEG) ; j1 < jMax ; j1++, wseg++)
    {
      posMin = wseg->weight > 0 && wseg->pos < posMin ? wseg->pos : posMin ;
      posMax = wseg->weight > 0 && wseg->pos + wseg->ln > posMax ? wseg->pos + wseg->ln : posMax ;
    }
  posMin = posMax ? posMin :  0 ;
  posMax++ ; 
  wig->posMin = posMin ;
  wig->iMin = 0 ;
  wig->iMax = posMax - posMin ;
  wig->nBp = 0 ;
  
  array (aa, wig->iMax - 1, unsigned int) = 0 ; /* make room */
  
  for (jj = 0, wseg = bigArrp (wsegs, jj, WSEG)  ; jj < jMax ; jj++, wseg++)
    {
      int x1 = wseg->pos ;
      int x2 = wseg->pos + wseg->ln ;
      if (wseg->weight > 0)
	{
	  wig->nBp += (long)wseg->weight * wseg->ln ;
	  for (int i = x1 ; i < x2 ;  i++)
	    array (aa, i - posMin, unsigned int) += wseg->weight ;
	}
    }
} /* wigFeed */

/**************************************************************/

void wigCompress (WIG wig, int step)
{
  Array aa = wig->aa ;
  int ii, jj = 0 ;
  int scale = wig->scale ? wig->scale : 1 ;
  int iMax = arrayMax (aa) ;
  int posMin = wig->posMin ; /* coordinate of offset 0 at step 1 */
  int newPosMin = 0 ;
  int newMin = 0 ;
  int newMax = 0 ;
  int demiScale = scale / 2 ;

  if (step == wig->step)
    { /* just divide by scale and find the true ends */
      unsigned int *up, *vp ;
      wig->nBp = 0 ;
      for (ii = jj = 0 , up = vp = arrp (aa, 0, unsigned int) ; ii < iMax ; ii++, up++ )
	{
	  unsigned int cumul = *up ;
	  if (newMax || 2*cumul > scale) /* alreday initialised or positive first value */
	    {
	      *vp = (cumul + demiScale)/scale ;
	      wig->nBp += *vp ;
	      cumul = 0 ;
	      newPosMin = newMax ? newPosMin : ii * step + posMin ;
	      newMin = newMax ? newMin : jj ;
	      newMax = *vp ? jj + 1 : newMax ;
	      jj++ ; vp++ ;
	    }
	}
    }
  else
    {  
      if (wig->step > 1)
	messcrash ("wigCompress can only be called on step=1 wig") ;
      if (step < 1)
	messcrash ("wigCompress called with step < 1: step = %d", step) ;

      int dx = posMin % step ;
      if (wig->nBp > demiScale)   /* oldScale, if false the wiggle is effectively empty of significant data */
	{
	  int ii, k, iMax = arrayMax (aa) ;
	  unsigned int *up, *vp, cumul = 0 ;
	  wig->nBp = 0 ;
	  k = (dx - step/2 + step) % step ;
	  for (ii = jj = 0 , up = vp = arrp (aa, 0, unsigned int) ; ii < iMax ; ii++, up++, k = (k + 1) % step)
	    {
	      cumul += *up ;
	      if (k == 0 || ii == iMax - 1)
		{
		  if (newMax || 2*cumul > scale) /* alreday initialised or positive first value */
		    {
		      *vp = (cumul + demiScale)/scale ;
		      wig->nBp += *vp ;
		      cumul = 0 ;
		      newPosMin = newMax ? newPosMin : ii + posMin ;
		      newMin = newMax ? newMin : jj ;
		      newMax = *vp ? jj + 1 : newMax ;
		      jj++ ; vp++ ;
		    }
		}
	    }
	}
    }

  wig->scale = 1 ;   /* new scale */
  wig->step = step ;
  wig->iMin = newMin ;  /* expect zero, the array should be cor4ectly compressed */
  wig->iMax = newMax ;
  wig->posMin = newPosMin / step ; wig->posMin *= step ;
  arrayMax (wig->aa) = newMax ;
} /* wigCompress */

/**************************************************************/

void wigCountGenes (WIG wig, BigArray geneB, BigArray geneC, long int *cumuls, char wigStrand)
{
  Array aaa = wig->aa ;
  long int ib, ibMax = bigArrayMax (geneB) ;
  int jMax = arrayMax (aaa) ;
  GBX *gb = bigArrayp (geneB, 0, GBX) ;
  unsigned int *xp = arrayp (aaa, 0, unsigned int) ;
  int step = wig->step ;
  
  /* the candidate gene segments that may cover position x are
   * not earlier than the first segment igOld  covering the previous position
   * not later than the fist gene ig1 starting after x
   */
  for (ib = 0, gb = bigArrp (geneB, ib, GBX) ; ib < ibMax ; ib++, gb++)
    {
      BOOL sameStrand = (gb->strand == wigStrand ? TRUE : FALSE) ;
      BOOL master = gb->flag & 0x10 ? TRUE : FALSE ;
      BOOL intronic = gb->flag & 0x1 ? TRUE : FALSE ;
      BOOL ambiguous = gb->flag & 0x8 ? TRUE : FALSE ;
      int friends = ambiguous && gb->friends > 0 ? gb->friends : 1 ;

      if (sameStrand && geneC && gb->a2 >= wig->posMin)
	{
	  int pos0 = wig->posMin ;
	  int j1 = (gb->a1 - pos0) / step ;
	  int j2 = (gb->a2 - pos0) / step ;
	  for (int j = (j1 > 0 ? j1 : 0) ; j <= j2 && j < jMax ; j++)
	    {
	      int weight = xp[j] ; 
	      if (weight)
		{
		  if (gb->gene)
		    {
		      int gg = gb->gene ;
		      if (! intronic)
			bigArray (geneC, 2*gg, long int) += weight/friends ;
		      bigArray (geneC, 2*gg + 1, long int) += weight/friends ;
		    }
		  if (master) cumuls[gb->flag & 0x7] += weight ;
		}
	    }
	}
    }
} /* wigCountGenes */

/**************************************************************/

#define BUF_SIZE (1 * 1024 * 1024)  // 1 MB - tune based on your system/disk and << pthread_create stack default size 4M 
static void wigExportDo (WIG wig, const char *fNam,  BOOL bigWig, const char *indexName, BOOL gzo)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = 0 ;
  char *writeBuffer = halloc(BUF_SIZE, h) ;
  int bufPos = 0 ;
  int step = wig->step ;
  const char *chromNam = wig->chromName ;

  Stack s = stackHandleCreate (1024, h) ;
  gzFile gzf = 0 ;

  if (wig->scale > 1)
    messcrash ("Call wigCompress before calling wigExport to reset scale=1") ;
  if (bigWig)
    {
      /* -clip avoids a crash if we go out of chromSizes */
      char *chromSizeNam = hprintf (h, "%s", indexName) ;
      char *cmd = hprintf (h, "wigToBigWig -clip stdin %s %s.bw", chromSizeNam, fNam) ;
      ao = aceOutCreateToPipe (cmd, h) ; // the bw format allows direct access
    }
  if (! ao) /* export as BF */ 
    {
      if (! gzo)
	ao = aceOutCreate (fNam, ".BF", 0, h) ;
      else  // compressing the BF files is extremelly good 100X
	{
	  ao = aceOutCreateToStack (s,h) ; 
	  gzf = gzopen (hprintf (h, "%s.gz", fNam), "wb") ;
	  if (! gzf) messcrash ("wigExportDo failed to gzopen(%s)\n", fNam) ;
	}
    }
  if (ao)
    {
      int start = wig->posMin ; if (start < 1) start = 1 ;
      aceOutDate (ao, "##", "wiggle") ; 
      aceOutf (ao, "track type=wiggle_0\n") ;
      aceOutf (ao, "fixedStep chrom=%s start=%d step=%d span=%d\n", chromNam, start, step, step) ;
      if (gzf)
	{
	  char *cp = stackText (s, 0) ;
	  int k = strlen (cp) ;
	  gzwrite (gzf, cp, k) ;
	}
    }
  Array aaa = wig->aa ;
  unsigned int *xp = arrayp (aaa, 0, unsigned int) ;
  for (int j = 0, jMax = arrayMax(aaa) ; j < jMax ; j++)
    {
      unsigned int w = xp[j] ;
      if (ao)
	{
	  char buf[32] ;
	  int k = fast_itoa_nl(buf, (int)(0.49999 + (double) w / step)) ;
	  
	  // Flush if buffer full
	  if (bufPos + k + 1 >= BUF_SIZE)
	    {
	      {
		if (gzf)
		  gzwrite (gzf, writeBuffer, bufPos) ;
		else if (ao)
		  aceOutBinary(ao, writeBuffer, bufPos) ;
	      }
	      bufPos = 0 ;
	    }
	  // Copy to buffer
	  memcpy(writeBuffer + bufPos, buf, k) ;
	  bufPos += k ;
	}
    }
  
  // Final flush
  if (bufPos > 0)
    {
      if (gzf)
	gzwrite (gzf, writeBuffer, bufPos) ;
      else if (ao)
	aceOutBinary(ao, writeBuffer, bufPos) ;
    }
  if (gzf)  
    {
      if (gzclose(gzf) != Z_OK)
	messcrash("gzclose failed");   /* important: the trailer/flush happens here */
    }
  ac_free (ao) ; /* ensure flushing the pipe is need be */
  ac_free (h) ;
  return ;
} /*  wigExportDo */

/**************************************************************/

void wigExportAZ (WIG wig, const char *fNam)
{
  AC_HANDLE h = ac_new_handle () ;
  wigCheck (wig, "wigExportAZ") ;
  if (wig->iMax)
    {
      wigCompress (wig, wig->step) ;
      if (wig->iMax)
	wigAzWrite (fNam, wig, h) ;
    }
  
  ac_free (h) ;
  return ;
} /*  wigExportAZ */

/**************************************************************/

void wigExportBF (WIG wig, const char *fNam, BOOL gzo)
{
  wigCheck (wig, "wigExportBF") ;
  if (wig->iMax)
    wigCompress (wig, wig->step) ;
  if (wig->iMax)
    wigExportDo (wig, fNam, FALSE, 0, gzo) ;

  return ;
} /*  wigExportBF */

/**************************************************************/

void wigExportBigWig (WIG wig, const char *fName, const char *indexName)
{
  wigCheck (wig, "wigExportBF") ;
  if (wig->iMax)
    wigCompress (wig, wig->step) ;
  if (wig->iMax)
    {
      if (! indexName)
	messcrash ("wigExportBigWig : null indexName (file of chrom sizes)") ;
      if (! filCheckName (indexName, 0, "r"))
	messcrash ("wigExportBigWig : cannot read file of chrom sizes : %s\n", indexName) ;
      if (! isExecutableInPath ("wigToBigWig"))   /* see w1/utils.c */
	wigExportDo (wig, fName, FALSE, 0, TRUE) ;    
      else
	wigExportDo (wig, fName, TRUE, indexName, FALSE) ;
    }
  return ;
} /*  wigExportBF */

/**************************************************************/
/**************************************************************/
/**************************************************************/
/*************************************************************************************/
/*************************************************************************************/
#define AZMAGIC 55182836
/* 256k positions == 2560 kb is step = 10 */

/**************************************************************/

static void azFinalize (void *vp)
{
  AZZ *az = vp ;
  if (!az) return ;
  if (az->magic != AZMAGIC)
    messcrash ("azFinalize received corrupt array->magic");
  az->magic = 0 ;
  if (az->h) ac_free (az->h) ;
  az->h = 0 ;
} /* azFinalize */

/**************************************************************/
#ifdef JUNK
static Array wigAzRegularTiling (AZZ *az, Array wPoints, AC_HANDLE h)
{
  Array aa = arrayHandleCreate (az->xMax, unsigned int, h) ;
  int iw, iwMax = arrayMax (wPoints) ;
  int step = az->step ;
  WIGGLEPOINT *wp = iwMax ? arrp (wPoints, 0, WIGGLEPOINT) : 0 ;
  for (iw = 0 ; iw < iwMax ; iw++, wp++)
    if (wp->x + step -1 >= 0)
      array (aa, (wp->x + step - 1)/step, int) = wp->y > 0 ? wp->y : 0 ;
  return aa ;
} /* wigAzRegularTiling */
#endif
/**************************************************************/
#define WIGBMAX 18
AZZ *wigAzWrite (const char *fName, WIG wig, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  int jMax = 0x1 << WIGBMAX ;
  Array byte0 = arrayHandleCreate (jMax, unsigned char, h) ;
  Array byte1 = arrayHandleCreate (jMax, unsigned char, h) ;
  Array byte2 = arrayHandleCreate (jMax, unsigned char, h) ;
  Array byte3 = arrayHandleCreate (jMax, unsigned char, h) ;
  Array aa = wig->aa ;
  int step = wig->step ;
  int iMax = wig->iMax ;
  int posMin = wig->posMin ;
  //  int posMax = wig->posMin + iMax * step ;

  AZZ *az = handleAlloc (azFinalize, h0, sizeof(AZZ)) ;
  az->magic = AZMAGIC ;
  az->h = ac_new_handle () ;
  az->fName = strnew (fName, az->h) ;

  if (az->fName)
    { /* remove .az suffix */
      char *cp = az->fName ;
      while (*cp) cp++ ;
      if (cp >= az->fName + 3) cp -= 3 ;
      if (! strcasecmp (cp, ".az")) *cp = 0 ;
      if (! *az->fName) messcrash ("wigAzOpen received a bad file name just called .az") ;
    }
  filremove (hprintf (h, "%s.az", fName), 0) ;

  az->target = strnew (wig->chromName, h0) ;
  az->step = wig->step ;
 #ifdef JUNK
  Array wPoints = 0 ;
  if (! aa) aa = wigAzRegularTiling (az, wPoints, h) ; 

  if (wPoints)
    {
      int wMax = arrayMax(aa) ;
      WIGGLEPOINT *wp = wMax ? arrp (wPoints, 0, WIGGLEPOINT) : 0 ;
      az->xMin = posMin = wp ? wp->x : 0 ;
      wp = wMax ? arrp (wPoints, wMax - 1, WIGGLEPOINT) : 0 ;
      posMax = posMin + (wMax - 1) * step ;
    }
 #endif
  posMin /= step ;
  int xMax = iMax + posMin ; 
  az->NB = xMax ? ((xMax - 1) >> WIGBMAX) + 1 : 0 ;
  az->xMax = iMax ;
  int xMin = 0 ;
  
  if (! filCreateDir (az->fName))
    messcrash ("wigAzWrite cannot create the %s directory", az->fName) ;

  ACEOUT ao = aceOutCreate (az->fName, ".az", FALSE, h) ;
  if (! ao) messcrash ("wigAzWrite cannot create %s/header", az->fName) ;
  
  aceOutf (ao, "Step\t%d\n", az->step) ;
  aceOutf (ao, "Span\t%d\n", 1) ;
  aceOutf (ao, "PosMin\t%d\n", xMin * step) ;  // always zero
  aceOutf (ao, "PosMax\t%d\n", xMax * step) ;
  aceOutf (ao, "xMax\t%d\n", xMax) ;
  aceOutf (ao, "NB\t%d\n", az->NB) ;
  aceOutf (ao, "BMAX\t%d\n", WIGBMAX) ;
  if (az->target) aceOutf (ao, "Target\t%s\n", wig->chromName) ;
	  
  aceOutDestroy (ao) ;

  int ii = 0 ;  /* scan aa */
  for (int iB = 0 ; iB < az->NB ; iB++)
    {
      unsigned int zMin = 0, zMax = 0 ;
      int nBytes = 0 ;
      unsigned int dz ;
      unsigned int  bMin = iB << WIGBMAX ;    /* offsets in this block */
      unsigned int  bMax = (iB + 1) << WIGBMAX ;

      for (int i = ii  ; i < iMax && i + posMin < bMax ; i++)
	if (i + posMin >= bMin)
	  {
	    unsigned int z = array (aa, i, unsigned int) ;
	    /* if (z < zMin) zMin = z ;   zMin on purpose was initialized at zero, so it never happens */
	    if (z > zMax) zMax = z ;
	  }
      
      nBytes = 0 ;
      dz = zMax - zMin ;
      if (dz) { nBytes++ ; dz >>= 8 ; }
      if (dz) { nBytes++ ; dz >>= 8 ; }
      if (dz) { nBytes++ ; dz >>= 8 ; }
      if (dz) { nBytes++ ; dz >>= 8 ; }
      
      char *azFileName = hprintf (h, "%s/wig.%d.az.gz",  az->fName, iB) ;
      gzFile gzf = gzopen (azFileName, "wb"); 
      if (! gzf)
	messcrash("wigAzWrite cannot open file %s",  azFileName) ;
      gzprintf (gzf, "%u\t%u\t%d\t%d\t%d\t%.3f\n", zMin, zMax, nBytes, wig->chromLength, wig->bases) ;

      array (byte0, jMax -1, unsigned char) = 0 ;
      array (byte1, jMax -1, unsigned char) = 0 ;
      array (byte2, jMax -1, unsigned char) = 0 ;
      array (byte3, jMax -1, unsigned char) = 0 ;

      if (nBytes >= 1) memset (arrp (byte0, 0, unsigned char), 0, jMax * sizeof (unsigned char))  ;
      if (nBytes >= 2) memset (arrp (byte1, 0, unsigned char), 0, jMax * sizeof (unsigned char))  ;
      if (nBytes >= 3) memset (arrp (byte2, 0, unsigned char), 0, jMax * sizeof (unsigned char))  ;
      if (nBytes >= 4) memset (arrp (byte3, 0, unsigned char), 0, jMax * sizeof (unsigned char))  ;


      for (  ; ii < iMax && ii + posMin < bMax ; ii++)
	if (ii + posMin >= bMin)
	  {
	    unsigned int z = array (aa, ii, unsigned int) ;
	    unsigned int dz = z - zMin;
	    int jj = ii + posMin - bMin ;

	    if (nBytes >= 1) { array (byte0, jj, unsigned char) = dz & 0xff ; dz >>= 8 ; }
	    if (nBytes >= 2) { array (byte1, jj, unsigned char) = dz & 0xff ; dz >>= 8 ; }
	    if (nBytes >= 3) { array (byte2, jj, unsigned char) = dz & 0xff ; dz >>= 8 ; }
	    if (nBytes >= 4) { array (byte3, jj, unsigned char) = dz & 0xff ; }
	  }
      
      unsigned int k = 0 ;
      if (nBytes >= 1)  k += gzwrite(gzf, arrp (byte0, 0, unsigned char), jMax) ;
      if (nBytes >= 2)  k += gzwrite(gzf, arrp (byte1, 0, unsigned char), jMax) ;
      if (nBytes >= 3)  k += gzwrite(gzf, arrp (byte2, 0, unsigned char), jMax) ;
      if (nBytes >= 4)  k += gzwrite(gzf, arrp (byte3, 0, unsigned char), jMax) ;
      
      if (k != (unsigned)(nBytes * jMax))
	messcrash("gzwrite k=%u != nBytes*jMax=%u", k, (unsigned) (nBytes  * jMax)) ;  
      if (gzclose(gzf) != Z_OK)
	messcrash("gzclose failed");   /* important: gzclose flushes the data */
    }

  ac_free (h) ;
  return az ;
} /* wigAzWrite */


/**************************************************************/

AZZ *wigAzOpen (const char *fName, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  
  AZZ *az = handleAlloc (azFinalize, h0, sizeof(AZZ)) ;
  az->magic = AZMAGIC ;
  az->h = ac_new_handle () ;
  az->step = 0 ;
  az->fName = strnew (fName, az->h) ;
  az->xMax  = 0 ;
  az->NB = 0 ;
  az->cNB = -1 ;
  az->bMax = WIGBMAX ;
  
  if (az->fName)
    { /* remove .az suffix */
      char *cp = az->fName ;
      while (*cp) cp++ ;
      if (cp >= az->fName + 3) cp -= 3 ;
      if (! strcasecmp (cp, ".az")) *cp = 0 ;
      if (! *az->fName) messcrash ("wigAzOpen received a bad file name just called .az") ;
    }

  ACEIN ai = aceInCreate  (hprintf (az->h, "%s.az", az->fName), FALSE, h) ;
  if (ai)
    {
      aceInSpecial (ai, "\n") ;
      while (aceInCard (ai))
	{
	  char *cp = aceInWord (ai) ;
	  if (! cp || *cp == '#')
	    continue ;
	  if (! strcmp (cp, "Step"))
	    { aceInStep (ai, '\t') ; aceInInt (ai, &az->step) ; }
	  else if (! strcmp (cp, "PosMin"))
	    { aceInStep (ai, '\t') ; aceInInt (ai, &az->posMin) ; }
	  else if (! strcmp (cp, "PosMax"))
	    { aceInStep (ai, '\t') ; aceInInt (ai, &az->posMax) ; }
	  else if (! strcmp (cp, "xMax"))
	    { aceInStep (ai, '\t') ; aceInInt (ai, &az->xMax) ; }
	  else if (! strcmp (cp, "NB"))
	    { aceInStep (ai, '\t') ; aceInInt (ai, &az->NB) ; }
	  else if (! strcmp (cp, "BMAX"))
	    { aceInStep (ai, '\t') ; aceInInt (ai, &az->bMax) ; }
	  else if (! strcmp (cp, "Target"))
	    { aceInStep (ai, '\t') ; cp = aceInWord (ai) ; if (cp) az->target = strnew (cp, az->h) ; }
	}
      if (az->step <= 0)
	messcrash ("Bad step=%d <=0, in wigAzOpen(%s)\n", az->step, fName) ;
    }
  else
    { ac_free (az) ; az = 0 ; }
  ac_free (h) ;
  
  return az ;
} /* azOpen */

/**************************************************************/

static BOOL azParseOneBlock (AZZ *az, int nb)
{
  if (az->cNB == nb)
    return TRUE ;
  AC_HANDLE h = ac_new_handle () ;
  char *azFileName = hprintf (h, "%s/wig.%d.az.gz", az->fName, nb) ;
  gzFile gzf = gzopen (azFileName, "rb") ;
  if (gzf)
    {
      int jMax = 0x1 << az->bMax ;
      unsigned int zMin, zMax ;
      int nBytes, chromLength ;
      float bases ;
      Array byte0 = arrayHandleCreate (jMax, unsigned char, h) ;
      Array byte1 = arrayHandleCreate (jMax, unsigned char, h) ;
      Array byte2 = arrayHandleCreate (jMax, unsigned char, h) ;
      Array byte3 = arrayHandleCreate (jMax, unsigned char, h) ;
      if (! az->cache) az->cache = arrayHandleCreate (jMax, int, az->h) ;
      array (az->cache, jMax - 1, int) = 0 ;  /* make room */
      az->cNB = nb ;
      az->c1 = ((long)nb) << az->bMax ;
      az->c2 = (nb+1L) << az->bMax ;
      az->chromLength = chromLength ;
      az->bases = bases ;

      array (byte0, jMax -1, unsigned char) = 0 ;
      array (byte1, jMax -1, unsigned char) = 0 ;
      array (byte2, jMax -1, unsigned char) = 0 ;
      array (byte3, jMax -1, unsigned char) = 0 ;

      /* header line: must match writer's gzprintf("%u\t%u\t%d\t%d\n", ...) */
      {
        char headerLine[256] ;
        if (! gzgets (gzf, headerLine, sizeof(headerLine)))
          messcrash ("wigAzRead cannot read header of %s", azFileName) ;
        if (sscanf (headerLine, "%u\t%u\t%d\t%d\t%f", &zMin, &zMax, &nBytes, &chromLength, &bases) < 4)
          messcrash ("wigAzRead bad header in %s: %s", azFileName, headerLine) ;
      }
      
      unsigned int k = 0 ;
      if (nBytes >= 1)  k += gzread (gzf, arrp (byte0, 0, unsigned char), jMax) ;
      if (nBytes >= 2)  k += gzread (gzf, arrp (byte1, 0, unsigned char), jMax) ;
      if (nBytes >= 3)  k += gzread (gzf, arrp (byte2, 0, unsigned char), jMax) ;
      if (nBytes >= 4)  k += gzread (gzf, arrp (byte3, 0, unsigned char), jMax) ;
      if (k != (unsigned)(nBytes * jMax))
        messcrash ("gzread k=%u != nBytes*jMax=%u in %s",
                   k, (unsigned)(nBytes * jMax), azFileName) ;

      if (gzclose (gzf) != Z_OK)
        messcrash ("gzclose failed on %s", azFileName) ;
      
      for (int jj = 0 ; jj < jMax ; jj++)
        {
          unsigned int z = 0 ;
          unsigned int dz = 0 ;
          if (nBytes >= 4) { dz += array (byte3, jj, unsigned char) ; dz <<= 8 ; }
          if (nBytes >= 3) { dz += array (byte2, jj, unsigned char) ; dz <<= 8 ; }
          if (nBytes >= 2) { dz += array (byte1, jj, unsigned char) ; dz <<= 8 ; }
          if (nBytes >= 1) { dz += array (byte0, jj, unsigned char) ; }
          z  = dz + zMin ;
          array (az->cache, jj , unsigned int) = z ;
        }
    }
  ac_free (h) ;
  return gzf ? TRUE : FALSE ;
} /* azParseOneBlock */

/**************************************************************/
/* x1, x2 are given as true coordinates */
BOOL wigAzZone (AZZ *az, int x1, int x2, int shift, WIG wig, Array wPoints, int *nPosp, long int *nBpp, BOOL cumul)
{
  if (! az)
    messcrash ("wigAzZone called on null az") ;
  if (az->magic != AZMAGIC)
    messcrash ("wigAzZone called on corrupted az") ;
  if (! wig && ! arrayExists (wPoints))
    messcrash ("wigAzZone called on null wPointsd Array") ;
  if (! nBpp)
    messcrash ("wigAzZone called on null nBpp") ;
  if (x2 <= x1)
    messcrash ("wigAzZone called x1>=x2: x1=%d x2=%d", x1, x2) ;
  if (x1 < 0)
    messcrash ("wigAzZone called x1 < 0: x1 = %d", x1) ;
  int step = az->step ;
  int posMin = 0 ;
  Array aa = wig ? wig->aa : 0 ;
  if (wig) wPoints = 0 ;
  
  if (shift > 0) shift += step - 1 ;   /* we prefer to over-reach a little when we divide by step */
  if (shift < 0) shift -= (step - 1) ;
  shift = shift/step  ;  
  
  /* x1, x2 are coordinates in the WIG, a1,a2 are coordinates in the block */
  /* all numberrs are divisible by step */
  if (cumul)
    { posMin = 0 ; x1 /= step ; x2 /= step ; if (wig) wig->posMin = posMin ; }
  else
    { posMin = x1 + step * shift ; shift = 0 ; x2 = (x2 - x1)/ step ; x1 = 0 ; if (wig) wig->posMin = posMin ; }
  int absShift = shift > 0 ? shift : -shift ;
  posMin /= step ;   /* now every number is an offSet */

  if (wig)
    {
      array (aa, x2 + absShift, unsigned int) += 0 ;
      if (! cumul)
	{
	  memset (arrp (aa, 0, unsigned int), 0, (x2 + absShift + 1) * sizeof(unsigned int)) ;
	  arrayMax (aa) = x2 + absShift + 1 ;
	  wig->iMin = 0 ; wig->iMax = arrayMax (aa) ; 
	  if (nPosp) *nPosp = arrayMax (wig ? aa : wPoints) ;
	  *nBpp = arrayMax (aa) ; /* this is not really the number of bases, but a proxy */
	}
    }
  else
    {
      array (wPoints, x2 + absShift, WIGGLEPOINT).x += 0 ;
      if (! cumul)
	{
	  memset (arrp (wPoints, 0, int), 0, (x2 + absShift + 1) * sizeof(WIGGLEPOINT)) ;
	  arrayMax (wPoints) = x2 + absShift + 1 ;
	  *nBpp = arrayMax (wPoints) ; /* this is not really the number of bases, but a proxy */
	}
    }
  
  if (x1 < 0)
    x1 = 0 ;  
  if (x2 + posMin > az->xMax)
    x2 = az->xMax - posMin ;
  
  int yMin = az->posMin / step ;
  
  for (int x = x1 ; x < x2 ; )
    {
      int y = x + posMin - yMin ;
      int y2 = x2 + posMin - yMin ; 
      int nb = y  >> az->bMax ;
      if (azParseOneBlock (az, nb))
	{
	  int c1 = az->c1 ;
	  int c2 = az->c2 ;
	  int dc1 = y > c1 ? y - c1 : 0 ;
	  int dc2 = y2 < c2 ? y2 - c1 : c2 - c1 ;
	  /* copy [dc1,dc2[ of the block to  [x, x + dc2 - dc1[ of the WIG */
	  if (wig)
	    {
	      wig->chromLength = az->chromLength ;
	      if (! cumul) /* just copy the blocks */
		{
		  memcpy (arrp (aa, x, unsigned int), arrp (az->cache, dc1, unsigned int), (dc2 - dc1)*sizeof (unsigned int)) ;
		  wig->bases = az->bases ;
		}
	      else
		{
		  wig->bases += az->bases ;
		  for (int i = dc1 ; i < dc2 ; i++)
		    {
		      int xx = az->c1 + i - posMin + shift ;
		      if (xx >= 0 && xx < arrayMax (aa))
			{
			  unsigned int *up = arrp (aa, xx, unsigned int) ;
			  *up += arr (az->cache, i, unsigned int) ;
			  *nBpp += arr (az->cache, i, unsigned int) ; 
			}
		    }
		}
	      
	    }
	  else
	    {
	      for (int i = dc1 ; i < dc2 ; i++)
		{
		  int xx = az->c1 + i - posMin + shift ;
		  if (xx >= 0 && xx < arrayMax (wPoints))
		    {
		      WIGGLEPOINT *wp = arrayp (wPoints, xx, WIGGLEPOINT) ;
		      wp->x = (xx + posMin + shift)  * step ;
		      if (cumul)
			wp->y += arr (az->cache, i, unsigned int) ;
		      else
			wp->y = arr (az->cache, i, unsigned int) ;
		      *nBpp += arr (az->cache, i, unsigned int) ; 
		    }
		}
	    }
	  x += (dc2 - dc1) ;
	}
      else
	return FALSE ;
    }
  if (wig) { wig->iMin = 0 ; wig->iMax = arrayMax (aa) ;}
  if (nPosp) *nPosp = arrayMax (wig ? aa : wPoints) ;    
  
  return TRUE ;
} /* azZone */
     
/**************************************************************/
/**************************************************************/
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
  if (cst->uu) free (cst->uu) ;
  if (cst->uv) free (cst->uv) ;
} /* scratchFinalize */

/**************************************************************/
/* Finally compute the best shift between the 2 strands */
void cisTransNormalize (CST cst)
{
  cisTransCheck (cst, "cisTransNormalize") ;
  long int nu = cst->nu ;
  int dxMax = cst->dxMax / cst->step ;
  if (cst->u2 < 1) cst->nu = 0 ;  
  if (cst->nu > 0)
    {
      /* remove mean values, possibly not needed */
         double u1 = cst->u1, v1 = cst->v1 ; 
	 
	 cst->u2 -= u1 * u1 / nu ;
	 cst->v2 -= v1 * v1 / nu ;
      
      double z = (cst->u2 > 0 && cst->v2 > 0) ? sqrt (cst->u2 * cst->v2) : 0 ;
      cst->bestUv = 0 ;
      cst->bestShift = 0 ;
      for (int dx = 0 ; dx < dxMax ; dx++)
	{
	  /* remove mean values, possibly not needed */
	  cst->uv[dx] -= u1 * v1 / nu ;
	  cst->uu[dx] -= u1 * u1 / nu ;
	  
	  cst->uv[dx] = z > 0.0001 ? cst->uv[dx] / z : 0 ;
	  cst->uu[dx] /= cst->u2 ;
	  if (cst->bestUv <= cst->uv[dx]  && cst->uv[dx] > 0)
	    {
	      cst->bestUv = cst->uv[dx] ;
	      cst->bestShift = dx * cst->step ;
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
/* Accumulate the counts from several chromosomes */
void cisTransCumulate (CST genomeCst, CST chromCst)
{
  CST cg = genomeCst, cc = chromCst ; // shorter nicknames
  
  cisTransCheck (cg, "cisTransCumulate") ;
  cisTransCheck (cc, "cisTransCumulate") ;
  if (! cg->nu)
    {
      int dxMax = cg->dxMax / cc->step ;
      int nn = dxMax * sizeof(double) ;
      cg->step = cc->step ;
      cg->uu = malloc (nn) ;
      cg->uv = malloc (nn) ;
      if (! cg->uu || ! cg->uv)
	messcrash ("malloc failed in cisTransCreate dxMax=%d\n", dxMax) ;
      memset (cg->uu, 0, nn) ;
      memset (cg->uv, 0, nn) ;
    }
  if (cg->step != cc->step)
    messcrash ("cisTransCumulate step mismatch %d != %d\n", cg->step, cc->step) ;
  if (cg->dxMax != cc->dxMax)
    messcrash ("cisTransCumulate dxMax mismatch %d != %d\n", cg->dxMax, cc->dxMax) ;
  cg->nu += cc->nu ;
  cg->u1 += cc->u1 ;
  cg->u2 += cc->u2 ;
  cg->v1 += cc->v1 ;
  cg->v2 += cc->v2 ;

  int dxMax = cg->dxMax / cg->step ;
  for (int dx = 0 ; dx < dxMax ; dx++)
    {
      cg->uu[dx] += cc->uu[dx] ;
      cg->uv[dx] += cc->uv[dx] ;
    }
  return ;
} /* cisTransCumulate */

/**************************************************************/
/**************************************************************/
/* Accelerate the code by accessing 4 memories at once
 * uu[dx] += z0*a[dx] + z1*a[dx+1] + z2*a[dx+2] + z3*a[dx+3] ;
 * same for uv with b
 */
static inline void kernel4 (double *restrict uu, double *restrict uv,
                            const double *restrict a, const double *restrict b,
                            double z0, double z1, double z2, double z3, int n)
{
  for (int dx = 0 ; dx < n ; dx++)
    {
      uu[dx] += z0*a[dx] + z1*a[dx+1] + z2*a[dx+2] + z3*a[dx+3] ;
      uv[dx] += z0*b[dx] + z1*b[dx+1] + z2*b[dx+2] + z3*b[dx+3] ;
    }
} /* kernel4 */

/**************************************************************/

static inline void kernel1 (double *restrict uu, double *restrict uv,
                            const double *restrict a, const double *restrict b,
                            double z, int n)
{
  for (int dx = 0 ; dx < n ; dx++)
    {
      uu[dx] += z * a[dx] ;
      uv[dx] += z * b[dx] ;
    }
} /* kernel1 */

/**************************************************************/

/* Feed the counts chromosome per chromosome  */
void cisTransFeed (CST cst, int step, WIG waF, WIG waR) 
{
  Array aaf = waF->aa ;
  Array aar = waR->aa ;
  int ii, iMax = aaf ? arrayMax (aaf) : 0 ; 
  int jj, jMax = aar ? arrayMax (aar) : 0 ;

  /* make room */
  cisTransCheck (cst, "cisTransFeed") ;
  ii = jj = 0 ;
  if (step <= 0)
    messcrash ("cisTransFeed called on bad step= %d <= 0", step) ;
  if (cst->step && cst->step != step)
    messcrash ("cisTransFeed can only feed wiggles computed bat the same step: previous=%d != new=%d", cst->step, step) ;
  cst->step = step ;
  int dxMax = cst->dxMax / step ;
  if (! cst->uu)
    {
      int nn = dxMax * sizeof(double) ;
      cst->step = cst->step ;
      cst->uu = malloc (nn) ;
      cst->uv = malloc (nn) ;
      if (! cst->uu || ! cst->uv)
	messcrash ("malloc failed in cisTransCreate dxMax=%d\n", dxMax) ;
      memset (cst->uu, 0, nn) ;
      memset (cst->uv, 0, nn) ;
    }
  long int nu = 0 ;
  double u1 = 0, v1 = 0, u2 = 0, v2 = 0 ;
  double *uu = cst->uu, *uv = cst->uv ;
  int ddx = waR->posMin - waF->posMin ;
  if (ddx > 0) ii = ddx/step ;
  else jj = (-ddx)/step ;
  
  if (ii >= iMax) return ;
  if (jj >= jMax) return ;
  
  /* convert once to double (note: arrp pointers start at the chromosome origin) */
  unsigned int *z1u = arrp (aaf, 0, unsigned int) ;
  unsigned int *z2u = arrp (aar, 0, unsigned int) ;
  double *A = malloc (iMax * sizeof(double)) ;
  double *B = malloc (jMax * sizeof(double)) ;
  if (!A || !B) messcrash ("malloc failed in cisTransFeed") ;
  for (int k = 0 ; k < iMax ; k++) A[k] = (double)(int) z1u[k] ;  /* int cast: fast cvtdq2pd */
  for (int k = 0 ; k < jMax ; k++) B[k] = (double)(int) z2u[k] ;

  /* the scalar sums: cheap, keep them in their own simple loop (vectorizes well) */
  {
    int n = (iMax - ii < jMax - jj) ? iMax - ii : jMax - jj ;
    for (int k = 0 ; k < n ; k++)
      {
        double z = A[ii+k], y = B[jj+k] ;
        u1 += z ; u2 += z*z ; v1 += y ; v2 += y*y ;
      }
    nu += n ;
  }

  /* bulk: blocks of 4 rows, full dxMax width guaranteed in range */
  int i0 = ii, j0 = jj ;
  int nn ;
  for (nn = 0 ; i0 + nn + 4 + dxMax <= iMax && j0 + nn + 4 + dxMax <= jMax ; nn += 4)
    {
      const double *a = A + i0 + nn, *b = B + j0 + nn ;
      double z0 = a[0], z1 = a[1], z2 = a[2], z3 = a[3] ;
      if (z0 != 0 || z1 != 0 || z2 != 0 || z3 != 0)       /* keep the sparsity skip */
        kernel4 (uu, uv, a, b, z0, z1, z2, z3, dxMax) ;
    }

  /* tail: with the shrinking dx range */
  for (int k = nn ; i0 + k < iMax && j0 + k < jMax ; k++)
    {
      double z = A[i0+k] ;
      if (z != 0)
        {
          int n = dxMax ;
          if (iMax - (i0+k) < n) n = iMax - (i0+k) ;
          if (jMax - (j0+k) < n) n = jMax - (j0+k) ;
          kernel1 (uu, uv, A + i0 + k, B + j0 + k, z, n) ;
        }
    }

  free (A) ; free (B) ;
  cst->nu += nu ;
  cst->u1 += u1 ;
  cst->v1 += v1 ;
  cst->u2 += u2 ;
  cst->v2 += v2 ;

  return ;
} /* cisTransFeed */

/**************************************************************/
/**************************************************************/

void cisTransExport (CST cst, ACEOUT ao)
{
  /*
    int *limitp, limits[] = {1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, -1} ; //, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 10000000, -1
  */
  cisTransCheck (cst, "cisTransExport") ;
  int step=cst->step ;
  
  if (ao)
    {
      aceOutDate (ao, "##", "Autocorrelation of the wiggle on the top strand, used as control\n") ;
      aceOutf (ao, "## Trans correlation of the 2 strand showing the average lag of the minus strand wiggle\n") ;
      int bestShift = cst->bestShift ;
      double bestUv =  cst->bestUv ;
      if (bestUv < .3)  /* do not provoke a shift on a weak bestUv */
	{ bestShift = 0 ; bestUv = cst->uv[0] ; }
      aceOutf (ao, "## bestShift=%d\t%.2f\n", bestShift, cst->bestUv) ;
      aceOutf (ao, "# Distance\tCis autocorrelation\tTrans correlation\t\n") ;
      
      int dxMax = cst->dxMax / step ;
      for (int dx = 0 ; dx < dxMax ; dx++)
	aceOutf (ao, "%d\t%g\t%g\n", step*dx, cst->uu[dx], cst->uv[dx]) ;
      aceOutf (ao, "\n\n") ;
    }  
} /* cisTransExport */

/**************************************************************/

CST cisTransCreate (int dxMax, AC_HANDLE h)
{
  int n1 = sizeof (struct cstStruct) ;
  CST cst = (CST) handleAlloc (cisTransFinalize, h, n1) ;
  memset (cst, 0, n1) ;
  
  cst->magic = CSTMAGIC ;
  cst->dxMax = dxMax ;

  return cst ;
} /* cisTransCreate */

/**************************************************************/
/**************************************************************/
/**************************************************************/


