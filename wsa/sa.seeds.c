/*
 * sa.seeds.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module implements all operations related
 * to the creation and sorting of the sequence seeds
*/

#include "sa.h"

/**************************************************************/

static void showCws (const PP *pp, BB *bb, BigArray cws)
{
  long int ii, iMax = bigArrayMax (cws), jj ;
  char buf[17], bufR[17] ;
  const long unsigned int mask26 = (1L << 26) - 1 ;
  
  for (ii = jj = 0 ; ii < iMax && jj < 50 ; ii++)
    {
      CW *cw = bigArrp (cws, ii, CW) ;
      if (1)
	{
	  if (!bb && ((cw->intron >> 31) & 0x1)==0)
	    continue ;  /* select the intron words */
	  else
	    {
	      if (cw->pos >= 0 && (cw->pos < 80 || cw->pos > 110)) continue ;
	      if (cw->pos <  0 && (-cw->pos < 80 || -cw->pos > 110)) continue ;
	    }
	}
      jj++ ;
      for (int i = 0 ; i < 16 ; i++)
	{
	  unsigned int z = (cw->seed >> (30 - 2 * i)) & 0x3 ;
	  switch (z)
	    {
	    case 0x0: buf[i] = 'A' ; break ;
	    case 0x1: buf[i] = 'C' ; break ;
	    case 0x2: buf[i] = 'G' ; break ;
	    case 0x3: buf[i] = 'T' ; break ;
	    }
	}
      buf[16] = 0 ;
      for (int i = 0 ; i < 16 ; i++)
	{
	  unsigned int z = (cw->seed >> (30 - 2 * i)) & 0x3 ;
	  switch (z)
	    {
	    case 0x0: bufR[15-i] = 'T' ; break ;
	    case 0x1: bufR[15-i] = 'G' ; break ;
	    case 0x2: bufR[15-i] = 'C' ; break ;
	    case 0x3: bufR[15-i] = 'A' ; break ;
	    }
	}
      bufR[16] = 0 ;
 
      if (1 || (cw->nam  == 44878))
	{
	  int da1 = 999, da = 0 ;
	  if (!bb && (ii % 256) < 2)
	    continue ;
	  if ((cw->intron >> 31) & 0x1)
	    {
	      da1 =  cw->intron & 0xf ;
	      da  =  ((cw->intron >> 4) & mask26) ;
	    }
	  printf (".. r=%d\t%d\t%u\t%s\t%s\t%s\tii=%ld\tda1=%d\tda=%d\n"
		  , cw->nam,cw->pos, cw->seed
		  , bb ? dictName(bb->dict, cw->nam >> 2) : dictName(pp->bbG.dict, cw->nam >> 1)
		  , buf, bufR, ii
		  , da1, da
		  ) ;
	}
    }
  printf ("......... max %ld\n", iMax) ;
} /* showCws */

/**************************************************************/

static int knownIntronOrder (const void *va, const void *vb)
{
  const GBX *up = va ;
  const GBX *vp = vb ;
  int n ;

  /* chrom order */
  n = up->chrom - vp->chrom ; if (n) return n ;
  /* strand order is implied by the parity of up->chrom */
  /* n = (up->a1 > up->a2) - (vp->a1 > vp->a2) ; if (n) return n ; */
  /* a1 order */
  n = up->a1 - vp->a1 ; if (n) return n ;
  return 0 ;
} /* knownIntronOrder */

/**************************************************************/

int saCodeIntronSeeds (PP *pp, BB *bbG)
{
  BOOL debug = FALSE ;
  CW *restrict cw  ;  
  const unsigned char *restrict cp ;
  Array dnas = bbG->dnas ;
  int k, p ;
  long int nn = 0 ;
  const int wLen = pp->seedLength ;
  const int nHidden = wLen > 16 ? wLen - 16 : 0 ;
  const int nHidden2 = nHidden << 1 ;
  const int nShift = 62 ;
  int NN = pp->nIndex ;
  const long unsigned int maskN = NN - 1 ;
  const long unsigned int mask32 = 0xffffffff ; /* 4 bytes integer */
  const long unsigned int maskSeedLn = (1L << 2*wLen) - 1 ;
  const unsigned int mask26 = (1L << 26) - 1 ;
  BigArray aa = pp->knownIntrons ;
  long int ii, iMax = bigArrayMax (aa) ;
  GBX *restrict upx = 0 ;
  int chrom = 0, a1, a2, da, v1, v2, dv, w1, w2, dw ;
  Array dna = 0 ;
  BOOL isIntronDown ;
  
  bigArraySort (aa, knownIntronOrder) ;
  
  upx = iMax ? bigArrp (aa, 0, GBX) : 0 ;
  for (ii = 0 ; ii < iMax ; ii++, upx++)
    {
      long unsigned int w, wr ;
      
      if (upx->chrom != chrom)
	{
	  chrom = upx->chrom >> 1 ;
	  dna = array (dnas, chrom, Array) ;
	}

      a1 = upx->a1 ;           /* 1-based Gt_ag position */
      a2 = upx->a2 ;        /* 1-based gt_aG position */
      isIntronDown = (chrom & 0x1) ? FALSE : TRUE ;
      
      da = a2 - a1 + 1 ;
      if (da >= (0x1 << 26))
	continue ; /* we use only 25 bits to code the intron length */

      /* check if left exon is shorter than 16 */
      v1 = v2 = 0 ; dv = 16 ;
      if (ii > 0)
	{
	  GBX *restrict vpx = upx - 1 ;
	  if (vpx->chrom == upx->chrom)
	    {
	      if (vpx->a1 < vpx->a2 && vpx->a2 > a1 - 17 && vpx->a2 < a1)
		{
		  v1 = vpx->a1 ; v2 = vpx->a2 ;
		  dv = a1 - v2 - 1 ; /* length of left exon */
		}
	    }
	}
      
      /* check if right exon is shorter than 16 */
      w1 = w2 = 0 ; dw = 16 ;
      if (ii < iMax - 1)
	{
	  GBX *restrict wpx = upx + 1 ;
	  if (wpx->chrom == upx->chrom)
	    {
	      if (wpx->a1 < wpx->a2 && wpx->a1 < a2 + 17 && wpx->a1 > a2)
		{
		  w1 = wpx->a1 ; w2 = wpx->a2 ;
		  dw = w1 - a2 - 1 ; /* length of right exon */
		}
	    }
	}
      
      /* Construct a 32 bases, 64 bits long unsigned integer representing wLen bases (A_, T_, G_, C_) using 2 bits per base 
       * first 15 bases belong to the left exon, upstream of Gt_ag
       * It does not matter if the annotated left or rigth exon is shorter than 15 bases
       */
      int dw0 = dw, dv0 = dv ;
      for (int pass = 0 ; pass < 4 ; pass++)
	{
	  /* In case of short exon  use either the dv and dw, or one of them or none */
	  switch (pass)
	    {
	    case 0: /* consider the up intron in isolation, forget the v and w introns */
	      dv = dw = 16 ;
	      break ;
	    case 1: /* consider the vp intron if it exists */
	      dv = dv0 ; dw = 16 ;
	      if (dv == 16) continue ;
	      break ;
	    case 2: /* consider the wp intron if it exists */
	      dv = 16 ; dw = dw0 ;
	      if (dw == 16) continue ;
	      break ;
	    case 3: /* consider both the vp and the wp intron if they both exit */
	      dv = dv0 ; dw = dw0 ;
	      if (dv == 16 || dw == 16) continue ;
	      break ;
	    } 
	  p = 0 ;
	  w = wr = 0 ;
	  
	  if (dv < 16) /* import some bases from 2 exons above */
	    {
	      cp = arrp (dna, v1 - (16 - dv) - 1, unsigned char) ;
	      for (int i = 0 ; i < 16 - dv ; i++, cp++)
		{
		  w <<= 2 ; wr >>= 2 ;
		  p++ ;
		  switch ((int)*cp)
		    {  /* alphabetic order and XOR complement */
		    case 0: p = 0 ; break ;
		    case A_:             wr |= 0x3L << nShift ; break ;
		    case C_: w |= 0x1L ; wr |= 0x2L << nShift ; break ;
		    case G_: w |= 0x2L ; wr |= 0x1L << nShift ; break ;
		    case T_: w |= 0x3L ;                   break ;
		      
		    default: p = 0 ;
		    }
		}
	    }
	  
	  /* import dv bases from the left exon */
	  cp = arrp (dna, a1 - dv - 1, unsigned char) ;
	  for (int i = 0 ; i < dv ; i++, cp++)
	    {
	      w <<= 2 ; wr >>= 2 ;
	      p++ ;
	      switch ((int)*cp)
		{  /* alphabetic order and XOR complement */
		case 0: p = 0 ; break ;
		case A_:             wr |= 0x3L << nShift ; break ;
		case C_: w |= 0x1L ; wr |= 0x2L << nShift ; break ;
		case G_: w |= 0x2L ; wr |= 0x1L << nShift ; break ;
		case T_: w |= 0x3L ;                   break ;
		  
		default: p = 0 ;
		}
	    }
	  
	  /* import dw bases from the right exon */
	  cp = arrp (dna, a2, unsigned char) ;
	  for (int i = 0 ; i < dw ; i++, cp++)
	    {
	      w <<= 2 ; wr >>= 2 ;
	      p++ ;
	      switch ((int)*cp)
		{  /* alphabetic order and XOR complement */
		case 0: p = 0 ; break ;
		case A_:             wr |= 0x3L << nShift ; break ;
		case C_: w |= 0x1L ; wr |= 0x2L << nShift ; break ;
		case G_: w |= 0x2L ; wr |= 0x1L << nShift ; break ;
		case T_: w |= 0x3L ;                   break ;
		  
		default: p = 0 ;
		}
	    }
	  /* import 16 - dw  bases from the following exon */
	  if (dw < 16)
	    {
	      cp = arrp (dna, w2, unsigned char) ; /* a2 is G in bio coords, so cp is the first base of the second exon */
	      for (int i = 0 ; i < 16 - dw  ; i++, cp++)
		{
		  w <<= 2 ; wr >>= 2 ;
		  p++ ;
		  switch ((int)*cp)
		    {  /* alphabetic order and XOR complement */
		    case 0: p = 0 ; break ;
		    case A_:             wr |= 0x3L << nShift ; break ;
		    case C_: w |= 0x1L ; wr |= 0x2L << nShift ; break ;
		    case G_: w |= 0x2L ; wr |= 0x1L << nShift ; break ;
		    case T_: w |= 0x3L ;                   break ;
		      
		    default: p = 0 ;
		    }
		}
	    }
	  
	  if (p == 32)   /* 4,5,6,7,8,9,10,11 = 8 possibilites, codes on 3 bits 21-24 */
	    for (int i1 = 4 - (wLen - 16)  ; i1 < 13 ; i1++) /* search every word with at least 4 bases on each side */
	      {
		if (i1 < 0) continue ;
		long unsigned int x =   (w >> (2 * (32 - wLen - i1))) & maskSeedLn ;
		long unsigned int xr =  (wr >> 2 * (i1)) & maskSeedLn ;
		int dx1 = 16 - i1 ;     /* number of letters in first exon */
		BOOL  minus  = (x > xr) ; /* word is antistrand to the exon-exon junction */
		long unsigned int word = ( minus ? xr : x ) ;
		k = word & maskN ;   /* right bits are used to select the relevant table */
		word >>= nHidden2 ;  /* remove excessive bits, by construction of NN, they are stored in k */
		word &= mask32 ;     /* masked down to 32 bits integer */
		cw = bigArrayp (bbG->cwsN[k], bigArrayMax(bbG->cwsN[k]), CW) ;
		cw->nam = (chrom << 1) | (minus ? 0x1 : 0) ;
		cw->seed = word ;
		cw->pos = a1 ; /* in introns, a1 cannot be zero */
		cw->intron =                  /* bit 30 alone, isExon */
		  (0x1 << 31)                 /* bit 31: isIntron */
		  | (isIntronDown ? (0x1 << 30) : 0) /* bit 30: intron strand */
		  | ((da & mask26)  << 4)            /* bits 5-30 : intron length */
		  | (dx1 & 0xf)          /* bits 0-4 : number of seed letters in the left exon, dx1 < 16 */
		  
		  ;
		nn++ ;
	      }
	}
    }
  fprintf (stderr, "+++++++ Coded %ld intron seeds\n", nn) ;
  if (debug && NN == 1)
    {
      BigArray uu = bigArrayCreate (nn, CW) ;
      CW *up, *vp = bigArrayp (bbG->cwsN[0], bigArrayMax (bbG->cwsN[0]) - nn, CW) ;
      for (int i = 0 ; i < nn ; i++, vp++)
	{
	  up = bigArrayp (uu, i, CW) ;
	  *up = *vp ;
	}
      showCws (pp, 0, uu) ;
    }
  return nn ;
} /* saCodeIntronSeeds */


/**************************************************************/
/* if we create N index, suppose N=16, we can mask log_4(NN), i.e. 2, letters
 * effectively indexing 18-mers in 32 bits words
 * using N index also accelerates sorting and searching
 * and consumes less memory for the sorting
 */
/*
actctatcacccaggctggagtgcagtggtgccatctcggctcactgcaacctccacctcccaggttcaatcgattctcctgcctcagcctcccgagtagctgggattataggcacccgccaccatgcccggctaatttttatattt
 actctatcacccaggctggagtgcagtggtgccatctcggctcactgcaacctccacctcccaggttcaatcgattctcctgcctcagcctcccgagtagctgggattataggcacccgccaccatgcccggctaatttttatattt
  actctatcacccaggctggagtgcagtggtgccatctcggctcactgcaacctccacctcccaggttcaatcgattctcctgcctcagcctcccgagtagctgggattataggcacccgccaccatgcccggctaatttttatattt
*/
void saCodeSequenceSeeds (const PP *pp, BB *bb, int step, BOOL isTarget) 
{
  int NN = pp->nIndex ;
  BigArray cwsN[NN] ;
  CW *restrict cw  ;  
  const unsigned char *restrict cp ;
  Array dnas = bb->dnas ;
  int k, ia, iaMax = arrayMax (dnas), nSeeds = 0 ;
  long int dMax = bigArrayMax (bb->globalDna) / (NN * step) ;
  const int wLen = pp->seedLength ;
  const int nHidden = wLen > 16 ? wLen - 16 : 0 ;
  const int nHidden2 = nHidden << 1 ;
  const int shiftB = 64 - 2 * wLen ;
  const long unsigned int maskN = NN - 1 ;
  const long unsigned int mask32 = 0xffffffff ; /* 4 bytes integer */
  const long unsigned int maskSeedLn = (1L << 2*wLen) - 1 ;
  int minEntropy = pp->minEntropy ;
  int minLength = pp->minLength ;
  
  memset (cwsN, 0, sizeof (cwsN)) ;
  if (minEntropy < 0) /* default value */
    {
      minEntropy = -minEntropy ;
      if (minEntropy > bb->runStat.p.maxReadLength / 2)
	minEntropy = bb->runStat.p.maxReadLength / 2 ;
    }
  if (minLength < 0) /* default value */
    {
      minLength = -minLength ;
      if (minLength > bb->runStat.p.maxReadLength / 2)
	minLength = bb->runStat.p.maxReadLength / 2 ;
    }
  if (step < 1)
    messcrash ("codeWordsDo received step = %d < 1", step) ;
  if (dMax >> 31)
    messcrash ("codeWordsDo received ln=%ls, step=%d, ln/step=%ld > 2G", bb->length, step, dMax) ;

  if (NN > 1)  /* Avoid reallocation of marginal size fluctuations */
    dMax *= 1.2 ;
  for (k = 0 ; k < NN ; k++)
    cwsN[k] = bigArrayHandleCreate (dMax, CW, bb->h) ;
  
  bb->step = step ;
  for (ia = 1 ; ia < iaMax ; ia++)
    {
      Array dna = arr (dnas, ia, Array) ;
      if (! dna)
	continue ;
      int ii, jj, p = 0, cStep, robin ;
      int iMax = arrayMax (dna) ;
      long unsigned int w, wr ;
      const int nShift = 62 ;
      if (1) step = iMax < 60 ? 1 : bb->step ;
      if (iMax < minLength)
	{
	  bb->runStat.tooShort++ ;
	  bb->runStat.tooShortBases += iMax ;
	  continue ;
	}
      if (minEntropy > 0 && oligoEntropy (arrp (dna, 0, unsigned char), iMax, minEntropy) < minEntropy)
	{
	  bb->runStat.lowEntropy++ ;
	  bb->runStat.lowEntropyBases += iMax ;
	  continue ;
	}
      w = wr = 0 ;
      /* each q vector is used round-robin style, so no copying is needed */
      BOOL qminus[step] ;
      long unsigned int qw[step], qwr[step], qx[step], qxr[step], qz[step], qhz[step] ;
      memset (qminus, 0, sizeof (qminus)) ;
      memset (qw, 0, sizeof (qw)) ;
      memset (qwr, 0, sizeof (qwr)) ;
      memset (qx, 0, sizeof (qx)) ;
      memset (qxr, 0, sizeof (qxr)) ;
      memset (qz, 0, sizeof (qz)) ;
      memset (qhz, 0, sizeof (qhz)) ;

      cp = arrp (dna, 0, unsigned char) ;

      for (ii = 0, jj = -wLen + 1, cStep = -wLen + 1, p = 0, robin = 0 ; ii < iMax ; ii++, jj++, cp++, cStep++, robin = (robin + 1) % step) 
	{
	  /* construct a 2*wLen bits long integer representing wLen bases (A_, T_, G_, C_) using 2 bits per base */
	  w <<= 2 ; wr >>= 2 ;
	  p++ ;
	  switch ((int)*cp)
	    {  /* alphabetic order and XOR complement */
	    case 0: p = 0 ; break ;
	    case A_:             wr |= 0x3L << nShift ; break ;
	    case C_: w |= 0x1L ; wr |= 0x2L << nShift ; break ;
	    case G_: w |= 0x2L ; wr |= 0x1L << nShift ; break ;
	    case T_: w |= 0x3L ;                   break ;
	      
	    default: p = 0 ;
	    }

	  qw[robin] = w ;
	  qwr[robin] = wr ;
	  qx[robin] = qw[robin] & maskSeedLn ;
	  qxr[robin] =  (qwr[robin] >> shiftB) & maskSeedLn ;
	  qminus[robin]  = (qx[robin] > qxr[robin]) ;
	  qz[robin] = qminus[robin] ? qxr[robin] : qx[robin] ;
	  qhz[robin] = (qz[robin] ^ (qz[robin] << 14)) ^ (qz[robin] << 23) ; 
	  if (
	      p < wLen ||
	      (((  w >> 2) << 2) == (w << 2) ) || /* avoid homopolymers */
	      (((  w >> 2) << 2) == (w << 2) )    /* avoid homodimers */
	      )
	    qz[robin] = 0 ;

	  
	  if (cStep == step) cStep = 0 ;
	  if (!cStep)
	    {
	      long int word = 0 ;
	      BOOL minus, qok = FALSE ;
	      int dx = 0 ;
	      
	      for (int i = step - 1 ; !qok && i >= 0 ; i--)
		if (qz[robin])
		  qok = TRUE ;
	      if (! qok)
		{ cStep = -1 ; continue ; } /* we hope to export the next word */
	      long int besthz = qhz[robin] ;
	      int bestk = robin ; dx = 0 ;
	      int myMax = (step > jj ? jj : step) ;
	      for (int i = 1 ; i < myMax ; i++)
		{
		  int k = (robin - i + step ) % step ;
		  if (qz[k] && qhz[k] < besthz)
		    { besthz = qhz[k]; bestk = k ; dx = i ; }
		}
	      word = qz[bestk] ;
	      minus = qminus[bestk] ;

	      k = word & maskN ;   /* right bits are used to select the relevant table */
	      word >>= nHidden2 ;  /* remove excessive bits, by construction of NN, they are stored in k */
	      word &= mask32 ;     /* masked down to 32 bits integer */
	      cw = bigArrayp (cwsN[k], bigArrayMax(cwsN[k]), CW) ;
	      /* store the sign in the name, semantically implying that we treat the 2 strands as independent chromosomes */
	      cw->nam = (ia << 1) | (minus ? 0x1 : 0) ;
	      cw->seed = word ;
	      cw->pos = jj - dx + 1 ;  /* bio coordinates of the first base of the seed */
	      cw->intron = 0 ;
	      nSeeds++ ;
	      cStep = dx ;
	      if (0)
		fprintf (stderr, "codeWordsDo p=%d k=%d pos=%d seed=%d\n", p, k, cw->pos,cw->seed) ;
	    }
	}
    }
  bb->cwsN = halloc (NN * sizeof(BigArray), bb->h) ;

  if (1)
    fprintf (stderr, "Agent %d lane %d allocated %d seeds\n", bb->readerAgent, bb->lane, nSeeds) ;
  for (k = 0 ; k < NN ; k++)
    bb->cwsN[k] = cwsN[k] ;
  if (0)
    {
      int mem = 0, mx = 0 ; /* megaBytes */
      messAllocStatus (&mem) ;
      messAllocMaxStatus (&mx) ;
      fprintf (stderr, "=== CodeWordDo finally allocated %d Mb, max %d Mb\n", mem, mx) ;
    }
  return ;
}  /* saCodeSequenceSeeds */

