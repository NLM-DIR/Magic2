/*
 * sa.wiggle.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module implements all operations related
 * to the "wiggle" coverage plots
 *   construction of the different types of wiggle
 *   evaluation of gene expression
 *   annotating transcripts starts and ends
 *
 */
#define ARRAY_CHECK

#define WIGGLETYPEMAX 2 /* strand */
#include "sa.h"
#include "wiggle.h"
#include "fastItoA.h"
#include <fcntl.h>  // for O_WRONLY if using write()
#include <unistd.h> // for write()


typedef struct wigPosStruct {
  unsigned int pos ; unsigned short ln ; unsigned short weight ; unsigned int dummy ;
} __attribute__((aligned(16))) WP ;

/**************************************************************/

static int wpOrder (const void *va, const void *vb)
{
  const WP *up = va ;
  const WP *vp = vb ;
  int n ;
  n = (up->pos > vp->pos ? 1 : (up->pos < vp->pos ? -1 : 0)) ; if (n) return n ;

  return 0 ;
} /* wiggleOrder */

/*************************************************************************************/

static void wiggleCumulate (BigArray aaa, BigArray aa)
{
  WP *up, *vp ;
  long int iMax = bigArrayMax (aa), iiMax = bigArrayMax (aaa) ;
  if (iMax)
    {
      saSort (aa, 1) ;
      up = bigArrayp (aaa, iMax + iiMax - 1, WP) ; /* make room */
      up = bigArrp (aaa, iiMax, WP) ;
      vp = bigArrp (aa, 0, WP) ;
      memcpy (up, vp, iMax * sizeof (WP)) ;
      vp = bigArrp (aa, iMax -1, WP) ;
      if ((int)vp->pos < 1)
	vp->pos = 1 ;
      if ((int)vp->pos < 1)
	messcrash ("bad negative value in wiggleCumulate pos=%d ln=%d\n", vp->pos, vp->ln) ;
    }
  return ;
} /* wiggleCumulate */
  
/**************************************************************/
#define STEP1 1  
static int wiggleCreate (const PP *pp, BB *bb)
{
  ALIGN *ap ;
  BOOL wiggleEnds = pp->wiggleEnds ;
  long int ii, iMax = bigArrayMax (bb->aligns) ;
  int nsw = 0 ;
  int chrom = 0 ;
  BigArray wig = 0, wigL = 0, wigR = 0, wigP = 0, wigNU = 0 ;
  int chromMax = dictMax (pp->bbG.dict) + 1 ;
  Array wiggles = bb->wiggles = arrayHandleCreate (2 * chromMax, BigArray, bb->h) ;
  Array wigglesL = bb->wigglesL = arrayHandleCreate (2 * chromMax, BigArray, bb->h) ;
  Array wigglesR = bb->wigglesR = arrayHandleCreate (2 * chromMax, BigArray, bb->h) ;
  Array wigglesP = bb->wigglesP = arrayHandleCreate (2 * chromMax, BigArray, bb->h) ;
  Array wigglesNU = bb->wigglesNU = arrayHandleCreate (2 * chromMax, BigArray, bb->h) ;
  const int step = (STEP1 == 1 ? 1 : pp->wiggle_step) ;  /* examples s=10, 5, 1 */
  const int demiStep = (step - 1)/2 ;
  const int endLength = 30 / step ;

  for (ii = 0, ap = bigArrp (bb->aligns, 0, ALIGN) ; ii < iMax ; ap++, ii++)
    {
      int mult = ap->nTargetRepeats ? ap->nTargetRepeats : 1 ;
      int weight = 720/mult ;
      int targetClass = ap->targetClass ;

      switch (targetClass)
	{
	case 'G':
	case 'M':
	case 'C':
	  break ;
	default:
	  continue ;
	}
      if (weight)
	{
	  BOOL isRead2 = ap->read & 0x1 ? TRUE : FALSE ;
	  int a1 = ap->a1 ;
	  int a2 = ap->a2 ;
	  int x1 = (ap->x1 == ap->chainX1 ? ap->x1 : 0) ;
	  int x2 = (ap->x2 == ap->chainX2 ? ap->x2 : ap->readLength) ;
	  int w1 = (a1 + demiStep)/step ;
	  int w2 = (a2 + demiStep)/step ;
	  int aChrom = ap->chrom ^ (ap->read & 0x1) ; 
	  if (0) aChrom ^= 0x1 ; /* negative run */

	  /*  chr 0  read 0 faux,
	   *  chr 0 read 1 ok
	   *  chr 1 read 1 faux
	   *  chr 1 read 0 ok
	   */
	  if (0 && !(ap->chrom & 0x1)) continue ;
	  if (0 && (ap->read & 0x1)) continue ;
	  if (0 &&  ((ap->chrom ^ ap->read) & 0x1)) continue ; 
	  if (chrom != aChrom)
	    {
	      nsw++ ;
	      wig = 0 ;
	      chrom = aChrom ;
	      if (*dictName (pp->bbG.dict, chrom >> 1) == 'G' ||
		  *dictName (pp->bbG.dict, chrom >> 1) == 'M' ||
		  *dictName (pp->bbG.dict, chrom >> 1) == 'C'
		  )
		{
		  wig = array (wiggles, chrom, BigArray) ;
		  if (! wig)
		    wig = array (wiggles,  chrom, BigArray) = bigArrayHandleCreate (100000, WP, bb->h) ;
		  wigP = array (wigglesP, chrom, BigArray) ;
		  if (! wigP)
		    wigP = array (wigglesP,  chrom, BigArray) = bigArrayHandleCreate (100000, WP, bb->h) ;
		  wigNU = array (wigglesNU, chrom, BigArray) ;
		  if (! wigNU)
		    wigNU = array (wigglesNU,  chrom, BigArray) = bigArrayHandleCreate (100000, WP, bb->h) ;
		}
	      if (wiggleEnds)
		{
		  if (*dictName (pp->bbG.dict, chrom >> 1) == 'G' ||
		      *dictName (pp->bbG.dict, chrom >> 1) == 'M' 
		      )
		    {
		      wigL = array (wigglesL, chrom, BigArray) ;
		      if (! wigL)
			wigL = array (wigglesL,  chrom, BigArray) = bigArrayHandleCreate (100000, WP, bb->h) ;
		      wigR = array (wigglesR, chrom, BigArray) ;
		      if (! wigR)
			wigR = array (wigglesR,  chrom, BigArray) = bigArrayHandleCreate (100000, WP, bb->h) ;
		    }
		}
	    }

	  if (pp->wiggleEnds && x1 == ap->leftClip + 1)
	    {
	      if (wigL && wigR && w1 < w2)
		{
		  BigArray wig = isRead2 ? wigR : wigL ;
		  WP *wp = bigArrayp (wig, bigArrayMax (wig), WP) ;
		  wp->pos = w1 ; wp->ln = endLength ; wp->weight = weight ;
		}
	      if (wigL && wigR && w1 > w2 && w1 >= endLength)
		{
		  BigArray wig = isRead2 ? wigR : wigL ;
		  WP *wp = bigArrayp (wig, bigArrayMax (wig), WP) ;
		  wp->pos = w1 - endLength ; wp->ln = endLength ; wp->weight = weight ;
		}
	    }
	  
	  if (wigP && w1 > 10 && w1 > endLength &&  x1 > ap->leftClip + 25)
	    {
	      WP *wp = bigArrayp (wigP, bigArrayMax (wigP), WP) ;
	      if (w1 < w2) { wp->pos = w1 - endLength ; wp->ln = endLength ; }
	      else { wp->pos = w1 ; wp->ln = endLength ; }
	      wp->weight = weight ;
	    }
	  if (wigP && w2 > 10 && w2 > endLength &&  x2 < ap->rightClip - 25)
	    {
	      WP *wp = bigArrayp (wigP, bigArrayMax (wigP), WP) ;
	      if (w1 < w2) { wp->pos = w2 ; wp->ln = endLength ; }
	      else { wp->pos = w2 - endLength ; wp->ln = endLength ; }
	      wp->weight = weight ;
	    }
	  if (wig && ap->nTargetRepeats == 1)
	    {
	      WP *wp = bigArrayp (wig, bigArrayMax (wig), WP) ;
	      if (w1 > w2) { int w0 = w1 ; w1 = w2 ; w2 = w0 ; }
	      wp->pos = w1 ; wp->ln = w2 - w1 + 1 ; wp->weight = weight ;
	      if (w1 < 0) messcrash ("bad w1 in wiggleCreate") ;
	    }
	  if (wigNU && ap->nTargetRepeats > 1)
	    {
	      WP *wp = bigArrayp (wigNU, bigArrayMax (wigNU), WP) ;
	      if (w1 > w2) { int w0 = w1 ; w1 = w2 ; w2 = w0 ; }
	      wp->pos = w1 ; wp->ln = w2 - w1 + 1 ; wp->weight = weight ;
	    }
	}
    }

  return  arrayMax (wiggles) ;
} /* wiggleCreate */

/*************************************************************************************/

static long int wigCumul (BigArray wig)
{
  long int ii, nn = 0 ;
  WP *wp = wig ? bigArrp (wig, 0, WP) : 0 ;
  long int iMax = wig ? bigArrayMax (wig) : 0 ;
  
  for (ii = 0 ; ii < iMax ; ii++, wp++)
    nn += wp->ln * wp->weight ;
  return nn / 72 ;
}

void saWiggleCumulate (const PP *pp, BB *bb)
{
  Array ppWiggles = 0 ;
  Array bbWiggles = 0 ;
  int chromMax = dictMax (pp->bbG.dict) + 1 ;
  int iwMax = wiggleCreate (pp, bb) ; /* max number of bb->wiggles */
  BigArray *ap1, *ap0 ;

  if (iwMax > 2 * chromMax) messcrash ("iwMax too large ?") ;

  for (int type = 0 ; type < (pp->wiggleEnds ? 5 : 1) ; type++)
    {
      switch (type)
	{
	case 0:
	  ppWiggles = pp->wiggles ;
	  bbWiggles = bb->wiggles ;
	  break ;
	case 1:
	  ppWiggles = pp->wigglesL ;
	  bbWiggles = bb->wigglesL ;
	  break ;
	case 2:
	  ppWiggles = pp->wigglesR ;
	  bbWiggles = bb->wigglesR ;
	break ;
	case 3:
	  ppWiggles = pp->wigglesP ;
	  bbWiggles = bb->wigglesP ;
	break ;
	case 4:
	  ppWiggles = pp->wigglesNU ;
	  bbWiggles = bb->wigglesNU ;
	break ;
	}
      
      iwMax = bbWiggles ? arrayMax (bbWiggles) : 0 ;
      if (! iwMax) continue ;
      if (iwMax > 2 * chromMax)
	messcrash ("Too many bbWigggles type %d :: %d >= %d", type, iwMax, 2*chromMax) ;
      ap0 = arrayp (ppWiggles, 0, BigArray) ;
      for (int iw = 0 ; iw < iwMax ; iw++)
	{
	  BigArray aa = array (bbWiggles, iw, BigArray) ;
	  
	  if (aa)
	    {
	      BigArray aaa = array (ppWiggles, 2 * bb->run * chromMax + iw, BigArray) ;
	      if (! aaa)
		aaa = array (ppWiggles, 2 * bb->run * chromMax + iw, BigArray) = bigArrayHandleCreate (10000, WP, pp->h) ;
	      long int naaa1 = 0 ;
	      if (0) naaa1 = wigCumul (aaa) ;
	      wiggleCumulate (aaa, aa) ;
	      if (0)
		{
		  long int naa = wigCumul (aa) ;
		  long int naaa2 = wigCumul (aaa) ;
		  fprintf (stderr, "%ld + %ld = %ld verif %ld\n"
			   , naa, naaa1, naaa2, naaa2-naa-naaa1) ;
		}
	    }
	}
      ap1 = arrp (ppWiggles, 0, BigArray) ;
      if (ap1 != ap0)
	messcrash ("pp->wiggles was relocalized which is not allowed here because of multithreading") ;

    }
  
  return ;
} /* saWiggleCumulate */

/**************************************************************/

#define BUF_SIZE (1 * 1024 * 1024)  // 4MB - tune based on your system/disk and << pthread_create stack default size 4M 

#ifdef JUNK
static void exportToFile (const char *fNam)
{
  
  /* Loop over your millions of numbers */
  for (;;/* each unsigned int n */)
    {
      int n = 0 ;
      int len = fast_itoa_nl_table(buf + pos, n) ;
      pos += len ;
      break ;
      if (pos >= BUF_SIZE - 12) {  // Flush if near full (safety margin for max line)
	if (fwrite(buf, pos, 1, fp) != 1) { /* error */ }
	pos = 0;
      }
    }

  /* Flush remainder */
    if (pos > 0)
      {
        if (fwrite(buf, pos, 1, fp) != 1) { /* error */ }
      }

    // Optional: fflush(fp); or fsync(fileno(fp)); if you need immediate disk write
    fclose(fp);
}
#endif


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


static void wiggleExportOne (const PP *pp, int nw, int type)
{
  Array wiggles = 0 ;
  BigArray wig = 0 ;
  int chromMax = dictMax (pp->bbG.dict) + 1 ;
  int run = nw / (2 * chromMax) ;
  int chrom = (nw % (2 * chromMax)) ;
  long int np = array(pp->runStats, run, RunSTAT).gt_ag_Support ;
  long int nm = array(pp->runStats, run, RunSTAT).ct_ac_Support ;
  BOOL doFlip = (pp->antiStrand || (! pp->strand && nm > 3 && 100L * nm > 80L * (nm+np))) ;
  char flip = doFlip ? 0x1 : 0x0 ;  ;  
  char strand = ( nw & 0x1) ^ flip ? 'r' : 'f' ;
  long int ii, iMax = 0 ;
  long int cumul = 0, cumuls[8] = {0} ;
  unsigned int pos0 ;
  BigArray geneC = 0 ;
  BigArray geneB = 0 ;
  const int wiggle_step = pp->wiggle_step ;
  const int demiStep = (wiggle_step - 1)/2 ;
  const char *typeNam ;
  char wigStrand = (strand == 'f' ? 0x0 : 0x1) ;

  if (0 && chrom != 2) return ;
  switch (type)
    {
    case 0:

      typeNam = (strand == 'f' ? "u.f" : "u.r") ;
      wiggles = pp->wiggles ;

      if (1)
	{
	  geneB = pp->geneBoxes ? array (pp->geneBoxes, chrom >> 1, BigArray) : 0 ;
	  geneC = pp->geneCounts ? array (pp->geneCounts, nw, BigArray) : 0 ;
	}
      
      break ;
      
    case 1:
      typeNam = (strand == 'f' ? "u.ELF" : "u.ERR") ;
      wiggles = pp->wigglesL ;
      break ;
      
    case 2:
      typeNam = (strand == 'f' ? "u.ERF" : "u.ELR") ;
      wiggles = pp->wigglesR ;
      break ;
      
    case 3:
      typeNam = (strand == 'f' ? "pp.f" : "pp.r") ;		 
      wiggles = pp->wigglesP ;
      break ;
      
    case 4:
      typeNam = (strand == 'f' ? "nu.f" : "nu.r") ;		 
      wiggles = pp->wigglesNU ;
      break ;
      
    }
  wig = array (wiggles, nw, BigArray) ;
  iMax = wig ? bigArrayMax (wig) : 0 ;
  if (! iMax) return ;
  
  if (0)
    bigArraySort (wig, wpOrder) ;
  else
    saSort (wig, 1) ;
    
  if (wig && iMax)
    {
      AC_HANDLE h = ac_new_handle () ;
      WP *wp0, *wp = bigArrp (wig, iMax - 1, WP) ;
      unsigned int posMax = wp->pos ;
      Array a = arrayHandleCreate (posMax + 1000, int, h) ;
      unsigned int *xp = arrayp (a, posMax, unsigned int) ;
      Array aAZ= 0 ;

      if (!pp->bigWig && pp->wigAZ && arrayMax(a))
	aAZ = arrayHandleCreate (arrayMax (a)/wiggle_step + 1, unsigned int, h) ;
      
      wp = bigArrp (wig, 0, WP) ;
      wp0 = bigArrp (wig, 0, WP) ;
      pos0 = wp0->pos ; pos0 -= pos0 % wiggle_step ;
      for (ii = 0 ; ii < iMax ; ii++, wp++)
	{
	  if (wp->weight)
	    {
	      xp = arrayp (a, wp->pos + wp->ln - pos0, unsigned int) ;
	      xp -= wp->ln ;
	      posMax = wp->pos + wp->ln - 1 ;
	      for (int i = 0 ; i < wp->ln ; i++)
		xp[i] += wp->weight ;
	    }
	}

      if (arrayMax(a))
	{
#define BUFFER_SIZE (2 * 1024 * 1024)  // 2MB
	  char *writeBuffer = malloc(BUFFER_SIZE) ;
	  int bufPos = 0 ;
	  int iAAZ = 0 ;
	  const char *chromNam = dictName (pp->bbG.dict, chrom >> 1) + 2 ;
	  const char *runNam = dictMax (pp->runDict) < run || ! run ? "runX" : dictName (pp->runDict, run) ;
	  ACEOUT ao = 0 ;
	  Stack s = stackHandleCreate (1024, h) ;
	  gzFile gzf = 0 ;
	  if (! ao && pp->bigWig)
	    {
	      /* -clip avoids a crash if we go out of chromSizes */
	      char *chromSizeNam = hprintf (h, "%s/targetSizes.txt", pp->indexName) ;
	      char *fNam = hprintf (h, "%s/wiggles/%s.%s.%s.bw", pp->outFileName, runNam, chromNam, typeNam) ;
	      char *cmd = hprintf (h, "wigToBigWig -clip stdin %s %s", chromSizeNam, fNam) ;
	      ao = aceOutCreateToPipe (cmd, h) ; // the bw format allows direct access
	    }
	  if (! ao && pp->wigBF)
	    {
	      const char *chromNam = dictName (pp->bbG.dict, chrom >> 1) + 2 ;
	      const char *runNam = dictMax (pp->runDict) < run || ! run ? "runX" : dictName (pp->runDict, run) ;
	      char *fNam = hprintf (h, "%s/wiggles/%s.%s.%s.BF%s", pp->outFileName, runNam, chromNam, typeNam, pp->gzo ? ".gz" : "") ;
	      if (! pp->gzo)
		ao = aceOutCreate (fNam, 0, 0, h) ;
	      else  // compressing the BF files is extremelly good 100X
		{
		  ao = aceOutCreateToStack (s,h) ; 
		  gzf = gzopen (fNam, "wb") ;
		}
	    }
	  if (ao)
	    {
	      aceOutDate (ao, "##", "wiggle") ; 
	      aceOutf (ao, "track type=wiggle_0\n") ;
	      aceOutf (ao, "fixedStep chrom=%s start=%d step=%d\n", chromNam, pos0, wiggle_step) ;
	      if (gzf)
		{
		  char *cp = stackText (s, 0) ;
		  int k = strlen (cp) ;
		  gzwrite (gzf, cp, k) ;
		}
	    }

      	  xp = arrayp (a, 0, unsigned int) ;
	  for (int localCumul = 0, j = 0, jMax = arrayMax(a) ; j < jMax ; j++)
	    {
	      unsigned int w = xp[j] ;
	      cumul += w ;
	      localCumul += w ;
	      if ((j + demiStep) % wiggle_step == 0)
		{
		  if (ao)
		    {
		      char buf[32] ;
		      int k = fast_itoa_nl(buf, localCumul / 720) ;
		      
		      // Flush if buffer full
		      if (bufPos + k + 1 >= BUFFER_SIZE)
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
		  else if (aAZ)
		    array (aAZ, iAAZ++, unsigned int) = localCumul / 720 ;
		  
		  localCumul = 0 ;
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
	  if (aAZ)
	    {
	      const char *chromNam = dictName (pp->bbG.dict, chrom >> 1) ;
	      const char *runNam = dictName (pp->runDict, run) ;
	      
	      AC_HANDLE h1 = ac_new_handle () ;
	      char *fNam = hprintf (h1, "%s/wiggles/%s.%s.%s.AZ", pp->outFileName, runNam, chromNam + 2, typeNam) ;
	      wigAzWrite (fNam, chromNam + 2, aAZ, 0, wiggle_step, pos0, posMax, h1) ;
	      ac_free (h1) ;
	    }
	  free(writeBuffer) ;
	}
      if (!pp->bigWig && pp->wigAZ && arrayMax(a))
      if (arrayMax(a) && geneB)
	{
	  long int ib, ibMax = bigArrayMax (geneB) ;
	  int jMax = arrayMax (a) ;
	  GBX *gb = bigArrayp (geneB, 0, GBX) ;
      	  xp = arrayp (a, 0, unsigned int) ;
	  
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
	      
	      if (sameStrand)
		{
		  int j1 = gb->a1 - pos0 ;
		  int j2 = gb->a2 - pos0 ;
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
	}

      ac_free (h) ;
    }
  if (type == 0 || type == 4)
    {
      arr (pp->wiggleCumuls, nw, long int) += cumul ;
      arr (pp->cdss, nw, long int) += cumuls[4];
      arr (pp->utrs, nw, long int) += cumuls[2] ;
      arr (pp->intronics, nw, long int) += cumuls[1] ;
      arr (pp->intergenics, nw, long int) += cumul - cumuls[1] - cumuls[2] - cumuls[4] ;
    }
  return ;
} /* wiggleExportOne */

/**************************************************************/

void wiggleExportAgent (const void *vp)
{
  const PP *pp = vp ;
  int nw ;
  
  while (channelGet (pp->wwChan, &nw, int))
    {
      wiggleExportOne (pp, nw, 0) ;
      if (pp->wiggleEnds)
	{
	  wiggleExportOne (pp, nw, 1) ;
	  wiggleExportOne (pp, nw, 2) ;
	  wiggleExportOne (pp, nw, 4) ;  /* non unique */
	  wiggleExportOne (pp, nw, 3) ;  /* partiel */
	}
      channelPut (pp->wwDoneChan, &nw, int) ;
    }
  channelCloseSource (pp->wwDoneChan) ;
  return ;
} /* wiggleExportAgent */

/**************************************************************/
/**************************************************************/

typedef struct gcStruct {
  int gene, run ;
  long int boxCount, exonCount ; 
} __attribute__((aligned(16))) GC ;

/**************************************************************/

static int gcOrder (const void *va, const void *vb)
{
  const GC *up = va ;
  const GC *vp = vb ;
  long int n ;
  
  n = up->gene - vp->gene ; if (n) return n > 0 ? 1 : -1 ;
  n = up->run - vp->run ; if (n) return n > 0 ? 1 : -1 ;

  return 0 ;
} /* gcOrder */

/**************************************************************/

static float geneIndex (const PP *pp, GC *gc)
{
  int run = gc->run ;
  RunSTAT *rs = arrayp (pp->runStats, run, RunSTAT) ;
  double wall = 1 ;
  double logDeux = log((double)2.0) ;
  BOOL isLow = FALSE ;
  float index = 0 ;  /* depends on total counts */
  double damper = 1000 ; /* damper in bp */
  int geneLength = 3000 ;
  long int average_read_ln = (rs->p.nBase1 + rs->p.nBase2) /(1L + rs->p.nReads) ;
  int ln = geneLength, ln0 ;
  int accessibleLength = 5000 ;
  double z, abp, bp, genomicKb ;

  
  if (rs->accessibleLength8kb > 0)
    accessibleLength = rs->accessibleLength8kb ;
  else if (rs->accessibleLength8kb > 0)
    accessibleLength = rs->accessibleLength8kb ;

  if (ln < 20)
    ln = 1000 ;
  if (average_read_ln < 35)    /* we cannot usefully map shorter reads */
    average_read_ln = 35 ;
  if (ln > accessibleLength)
    ln = accessibleLength ;
  
  /*
    if (gx->isMicroRNA) { ln = 1000 ;  average_read_ln = 35 ; }
    if (gx->isMA)
    {
    ln = ln0 = average_read_ln = gx->isMA ;
    }
  */

  if (rs->nBaseAligned1 + rs->nBaseAligned2 < (0x1 << 23)) /* 8M */
    return -999 ;
  if (gc->boxCount < 100000)  /* 100 kbases */
    return -999 ;
  
  if (rs->p.nPairs)
    wall *= 2 ;



#ifdef JUNK
  z = 5.0/100.0 ;          /* dromadaire */
      genomicKb = rc->genomicKb - z * rc->targetKb ;
      if (genomicKb < 0) 
	genomicKb = 0 ;
      genomicKb = gx->seaLevel * genomicKb / 3 ;  
      if (rc->DGE || gx->DGE || rc->pA == 1 || (!rc->pA && gx->pA == 1)) genomicKb /= 10 ;      /* total, 90% go to introns and should not be attributed to genomic contamination */

      if (rc->intergenicKb)
	genomicKb = rc->intergenicKb * gx->seaLevel/3 ; /* was 0.4 seaLevel * rc->intergenicKb, sellevel default == 3 */
      if (rc->intergenicKb > rc->prunedTargetKb)  /* 2020_01_27 we are facing a DNA-seq experiment */
	genomicKb *= 2 ;
      noise_bp = ln0 * genomicKb/gx->genomeLengthInKb ;  /* here use the true length of the sequence, do not staturate at 3kb  */
     
      noise_tags = noise_bp/average_read_ln ;
      if (noise_tags < 10)
	noise_bp = 3 * noise_tags ;
      else if (noise_tags > 500)
	noise_bp = 1.2828 * noise_tags ;
      else
	noise_bp = noise_tags  + 6.3245 * sqrt (noise_tags) ;
      noise_bp *= average_read_ln ;
      noise_tags = noise_bp/average_read_ln ;
    }

  /* mieg: 2011_03_14 we try to substract the genomic noise
   * note that we should measure the effective length rather than rely on annotations
   * CF Double stranded or star-seq normalisation is supposed to yield a more uniform coverage
   * and may also influence the index is some way or other
   */
  index = -999 ;

  bp = 1000 * dc->kb - noise_bp ;  
  abp = 1000 * (rc->prunedTargetKb + rc->intergenicKb * gx->seaLevel/3.0) ; /* 2020_01_27 add the intergenic tags to the denominator */
  bonusTags = dc->tags - noise_tags - wall ;

  /****************** new method 2013_04_13
   * we separate damping the (quantized)  number of tags 
   * from the denominator, normalizing the size of the genes and the run
   * for large counts, the formula is not modified
   * for zero counts, the low value is around 7 for 60M reads (Rhs844)
   */

  if (! (LEMING && isIntron)) /* new method 2013_04_13*/
    {
      double tt ;
      int lf ;

      if (bonusTags < 0)
	{ 
	  isLow = TRUE ;
	  *isLowp = TRUE ;
	}
      else
	*isLowp = FALSE ;

      if (gx->isMA)
	lf = 0 ;
      else if (rc->fragmentLength)
	lf = rc->fragmentLength ;
      else
	lf = 300 ; /* assumed length of the fragment, we should have measured it */ 

      z = 1.0E12/ ((ln - lf > lf ? ln - lf : lf) * abp) ;
      
      if (1) /* excellent change: OS and Fav d not change, EF HREF HROS go up one point */
	{ /* 2013_07_31 revert to previous method */
	  bp = 1000 * dc->kb - noise_bp ;  
	  tt = bp/average_tag_ln ; z *= average_tag_ln ;
	  if (tt < 0) tt = 0 ; 
	  damper = 5 ;
	  if (1) damper = 3 ; /* SEQC2 test 2021_09_03 */
	  damper = 4 ; /* 2023_09_21, so if tt=0 log((tt+sqrt(damper +tt^2)/2)=log(1)=0 */
	}

      if(gx->isMA) damper = 9 ; /* because we count a local coverage, not a number of reads */
      index = (log(( tt + sqrt (damper  + tt * tt))/((double)2.0)) + log(z))/logDeux ;	 
      if (! rc->zeroIndex)
	{
	  /* mieg 2023_09_21: rc->zero cannot depend on the particular gene 
	   * i choose 2kb
	   */
	  z = 1.0E12/ (2000 * abp) ; z *= average_tag_ln ; 
	  rc->zeroIndex = (log(sqrt(damper)/2.0) + log(z))/logDeux ; 
	}
      index += rc->shift_expression_index ;
      return index ;
    }

  index += rc->shift_expression_index ;



#endif
  
  return index ;
} /* geneIndex */

/**************************************************************/

static GeneCounts wiggleExportGeneCounts (const PP *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  int nw, wMax = arrayMax (pp->wiggles) ;
  int chromMax = dictMax (pp->bbG.dict) + 1 ;
  BigArray allGeneC ;
  long int igc = 0, jgc, igcMax = 0, nnn = 0 ;
  GC *gc, *gc2 ; 
  char tBuf[25] ;
  GeneCounts gcs = {0} ;
  int ngenes = 0 ;
  
  fprintf (stderr, "%s: start geneCounts export\n", timeBufShowNow (tBuf)) ;
  
  allGeneC = bigArrayHandleCreate (1000000, GC, h) ;
  for (nw = 0 ; nw < wMax ; nw++)
    {
      BigArray geneC = array (pp->geneCounts, nw, BigArray) ;
      int run = nw / (2 * chromMax) ;
      long int gCMax = geneC ? bigArrayMax (geneC) : 0 ;
      int gene ;
      long int *xp ;
      if (gCMax)
	for (gene = 0, xp = bigArrp (geneC, 0, long int) ; gene < gCMax ; gene += 2, xp+=2)
	  if (*xp > 0)
	    {
	      gc = bigArrayp (allGeneC, igcMax++, GC) ;
	      gc->gene = gene >> 1 ;
	      gc->run = run ;
	      gc->exonCount += xp[0] ;
	      gc->boxCount += xp[1] ;
	    }
    }
  if (igcMax)
    {
      bigArraySort (allGeneC, gcOrder) ;
      for (igc = jgc = 0, gc = bigArrp (allGeneC, 0, GC), gc2 = gc ; igc < igcMax ; igc++, gc++)
	{
	  if (gc2->gene != gc->gene || gc2->run != gc->run)
	    {
	      gc2++ ; jgc++ ;
	      if (gc2 < gc) *gc2 = *gc ;
	    }
	  else if (gc2 < gc)
	    {
	      gc2->exonCount += gc->exonCount ;
	      gc2->boxCount += gc->boxCount ;
	    }
	}
      bigArrayMax (allGeneC) = jgc ;

      
      ACEOUT ao = aceOutCreate (pp->outFileName, ".geneCounts.tsf", pp->gzo, h) ;
      aceOutDate (ao, "##", "Gene counts") ;
      aceOutf (ao, "#Gene\tRun\tFormat\tGene_coverage\tExons_Xocoverage\n") ;

      Array geneCoords = pp->geneCoords ;
      GBX *gCo = 0 ;
      int gbxMax = geneCoords ? arrayMax (geneCoords) : 0 ;
      
      for (igc = jgc = 0, gc = bigArrp (allGeneC, 0, GC), gc2 = gc ; igc < igcMax ; igc++, gc++)
	if (gc->gene && gc->boxCount)
	  {
	    nnn+= gc->boxCount ; ngenes++ ;
	    gCo = gc->gene < gbxMax ? arrp (geneCoords, gc->gene, GBX) : 0 ;
	    aceOutf (ao, "%s\t%s\ttiifii\t%s\t%d\t%d\t%.2f\t%ld\t%ld\n"
		     , dictName (pp->geneDict, gc->gene)
		     , dictName (pp->runDict, gc->run)
		     , gCo && gCo->chrom ? dictName (pp->bbG.dict, gCo->chrom) : "-"
		     , gCo ? gCo->a1 : 0, gCo ? gCo->a2 : 0
		     , geneIndex (pp, gc)
		     , gc->boxCount/720, gc->exonCount/720
		     ) ;
	  }
    }
  
  fprintf (stderr, "%s: stop geneCounts export %d genes with total count %ld\n", timeBufShowNow (tBuf), ngenes, nnn/720) ;
  ac_free (h) ;
  return gcs ;
} /* wiggleExportGeneCounts */

/**************************************************************/

static void wiggleExportWiggleStats (PP *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  int wMax = arrayMax (pp->wiggles) ;
  int chromMax = dictMax (pp->bbG.dict) + 1 ;
  long int nnn = 0 ;
  ACEOUT ao = aceOutCreate (pp->outFileName, ".wiggleCumuls.tsf", pp->gzo, h) ;
  aceOutDate (ao, "##", "wiggleCumuls in Million Bases") ;
  aceOutDate (ao, "#", "Target\tRun\tFormat\tCumul\tCDS\tUTR\tIntronic\tIntergenic in Bases") ;
  static int n720 = 720 ;
  
  /* export the wiggle cumuls per run. chromosome, strand */
  for (int nw = 0 ; nw < wMax ; nw++)
    {
      int run = nw / (2 * chromMax) ;
      int chrom = (nw % (2 * chromMax)) ;
      char strand = ( nw & 0x1) ? 'r' : 'f' ;
      long int cumul = array (pp->wiggleCumuls, nw, long int) ;
      long int cds = array (pp->cdss, nw, long int) ;
      long int utr = array (pp->utrs, nw, long int) ;
      long int intronic = array (pp->intronics, nw, long int) ;
      long int intergenic = array (pp->intergenics, nw, long int) ;

      if (cumul)
	{
	  const char *chromNam = dictName (pp->bbG.dict, chrom >> 1) ;
	  const char *runNam = dictName (pp->runDict, run) ;

	  RunSTAT *rc = arrayp (pp->runStats, run, RunSTAT) ;
	  rc->wiggleCumul += cumul ;
	  rc->cds += cds ;
	  rc->utr += utr ;
	  rc->intronic += intronic ;
	  rc->intergenic += intergenic ;
	  nnn += cumul ;
	  aceOutf (ao, "%s.%c\t%s\tiiiii\t%ld\t%ld\t%ld\t%ld\t%ld\n"
		   , chromNam, strand
		   , runNam
		   , cumul / n720 
		   , cds/n720, utr/n720, intronic / n720 , intergenic / n720 
		   ) ;	       
	}
    }

  /* export the cumul per run and globally */
  pp->wiggleCumul = 0 ;
  pp->cds = 0 ;
  pp->utr = 0 ;
  pp->intronic = 0 ;
  pp->intergenic = 0 ;
  for (int run = 1 ; run <= dictMax (pp->runDict) ; run++)
    {
      RunSTAT *rc = arrayp (pp->runStats, run, RunSTAT) ;

      rc->wiggleCumul /= n720 ;
      rc->utr /= n720 ;
      rc->cds /= n720 ;
      rc->intronic /= n720 ;
      rc->intergenic /= n720 ;
      
      pp->wiggleCumul += rc->wiggleCumul ;
      pp->cds += rc->cds ;
      pp->utr += rc->utr ;
      pp->intronic += rc->intronic ;
      pp->intergenic += rc->intergenic ;

      const char *runNam = dictName (pp->runDict, run) ;
      
      aceOutf (ao, "%s\t%s\tiiiii\t%ld\t%ld\t%ld\t%ld\t%ld\n"
	       , "Any"
	       , runNam
	       , rc->wiggleCumul
	       , rc->cds
	       , rc->utr
	       , rc->intronic
	       , rc->intergenic
	       ) ;
    }
  
  ac_free (h) ;
  return ;
} /* wiggleExportWiggleStats */

/**************************************************************/
/**************************************************************/

GeneCounts saWiggleExport (PP *pp, int nAgents)
{
  AC_HANDLE h = ac_new_handle () ;
  int wMax = arrayMax (pp->wiggles) ;
  BOOL debug = FALSE ;
  char tBuf[25] ;
  GeneCounts gcs = {0} ;

  nAgents = nAgents > 32 ? 32 : nAgents ;
  
  pp->wiggleCumuls = arrayHandleCreate (wMax, long int, h) ;
  pp->cdss = arrayHandleCreate (wMax, long int, h) ;
  pp->utrs = arrayHandleCreate (wMax, long int, h) ;
  pp->intronics = arrayHandleCreate (wMax, long int, h) ;
  pp->intergenics = arrayHandleCreate (wMax, long int, h) ;

  /* initialize now to prevent any race condition on arrayMax */
  arrayMax (pp->wiggleCumuls) = wMax - 1 ;
  arrayMax (pp->cdss) = wMax - 1 ;
  arrayMax (pp->utrs) = wMax - 1 ;
  arrayMax (pp->intronics) = wMax - 1 ;
  arrayMax (pp->intergenics) = wMax - 1 ;
  
  if (pp->geneBoxes)
    pp->geneCounts = arrayHandleCreate (wMax, Array, h) ;

  int k = 0, n = 0 ;
  for (int nw = 0 ; nw < wMax ; nw++)
    {
      Array wig = arr (pp->wiggles, nw, Array) ;
      if (wig)
	n++ ;
    }
  fprintf (stderr, "%s: start exportation of  %d wiggles\n", timeBufShowNow (tBuf), n) ;
  /* parallelize: open a channel and start agents */
  pp->wwChan = channelCreate (wMax + 1, int, h) ;
  channelDebug (pp->wwChan, debug, "wwChan") ;
  pp->wwDoneChan = channelCreate (wMax+1, int, h) ;
  channelDebug (pp->wwDoneChan, debug, "wwChan") ;

  BOOL doChan = TRUE ;

  channelSources (pp->wwDoneChan, nAgents) ;
  if (doChan)
    for (int ii = 0 ; ii < nAgents ; ii++)
      {
	pp->agent = ii ;
	wego_go (wiggleExportAgent, pp, PP) ; 
      }
  
  /* load the channel to start execution */ 

  k = n = 0 ;
  for (int nw = 0 ; nw < wMax ; nw++)
    {
      Array wig = arr (pp->wiggles, nw, Array) ;
      if (wig)
	{
	  if (pp->geneCounts)
	    {
	      array (pp->geneCounts, nw, BigArray) = bigArrayHandleCreate (2048, long int, h) ;
	    }
	  if (doChan)
	    channelPut (pp->wwChan, &nw, int) ;
	  else
	    {
	      fprintf (stderr, "%s: wiggle export %d\n", timeBufShowNow (tBuf), nw) ;
	      wiggleExportOne (pp, nw, 0) ;
	      if (pp->wiggleEnds)
		{
		  wiggleExportOne (pp, nw, 1) ;
		  wiggleExportOne (pp, nw, 2) ;
		  wiggleExportOne (pp, nw, 4) ;  /* non unique */
		  wiggleExportOne (pp, nw, 3) ;  /* partiel */
		}
	    }
	}
    }
  channelClose (pp->wwChan) ;

  /* synchronize */
  if (doChan)
    while (channelGet (pp->wwDoneChan, &k, int))
      n++ ;

  fprintf (stderr, "%s: stop wiggle export\n", timeBufShowNow (tBuf)) ;

  if (pp->geneBoxes)
    gcs = wiggleExportGeneCounts (pp) ;

  wiggleExportWiggleStats (pp) ;
  
  ac_free (h) ;
  return gcs ;
}

/**************************************************************/
/**************************************************************/
/**************************************************************/
 
