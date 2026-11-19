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
#include "sa.seedEncoder.h"
/**************************************************************/

static void showCws (const PP *pp, BB *bb, BigArray cws)
{
  long int ii, iMax = bigArrayMax (cws), jj ;
  char buf[17], bufR[17] ;
  const unsigned long long mask26 = (1ULL << 26) - 1 ;
  
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
	  printf (".. r=%u\t%u\t%u\t%s\t%s\t%s\tii=%ld\tda1=%d\tda=%d\n"
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

static int dnaOrder (const void *va, const void *vb)
{
  const Array *up = va ;
  const Array *vp = vb ;
  int n ;

  if (! up)
    return vp ? -1 : 0 ;
  if (! vp)
    return 1 ;
  if (! *up)
    return *vp ? -1 : 0 ;
  if (! *vp)
    return 1 ;
  int iMax = arrayMax (*up) ;
  int jMax = arrayMax (*vp) ;
  n = iMax - jMax ; if (n) return -n ;  // short first

  unsigned const char *restrict cp = arrp (*up, 0, unsigned char) ;
  unsigned const char *restrict cq = arrp (*vp, 0, unsigned char) ;
  for (int i = 0 ; i < iMax ; i++, cp++, cq++)
    {
      n = *cp < *cq ? -1 : (*cp > *cq ? 1 : 0) ;
      if (n) return n ;
    }
  return 0 ;
} /* knownIntronOrder */

/**************************************************************/

static int knownIntronOrder (const void *va, const void *vb)
{
  const GBX *up = va ;
  const GBX *vp = vb ;
  int n ;

  /* chrom order */
  n = (up->chrom > vp->chrom) - (up->chrom < vp->chrom) ; if (n) return n ; 
  /* strand order is implied by the parity of up->chrom */
  /* a1 order */
  n = (up->a1 > vp->a1) - (up->a1 < vp->a1) ; if (n) return n ; 

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
  int NN = pp->nIndex ;
  const unsigned long long maskN = (unsigned long long)NN - 1 ;
  const unsigned long long mask32 = 0xffffffff ; /* 4 bytes integer */
  const unsigned long long maskSeedLn = (1ULL << 2*wLen) - 1 ;
  const unsigned int mask26 = (1ULL << 26) - 1 ;
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
      unsigned long long  w, wr ;
      
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
      
      /* Construct a 32 bases, 64 bits unsigned long long representing wLen bases (A_, T_, G_, C_) using 2 bits per base 
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
		  unsigned int base = (*cp & 0x0F) ;
		  
		  w  |= w_table[base] ;    /* construct the shifted seed */
		  wr |= wr_table[base] ;    /* construct the shifted reversed seed */
		  p   = p & valid_base[base] ;
		}
	    }
	  
	  /* import dv bases from the left exon */
	  cp = arrp (dna, a1 - dv - 1, unsigned char) ;
	  for (int i = 0 ; i < dv ; i++, cp++)
	    {
	      w <<= 2 ; wr >>= 2 ;
	      p++ ;
	      unsigned int base = (*cp & 0x0F) ;
	      
	      w  |= w_table[base] ;    /* construct the shifted seed */
	      wr |= wr_table[base] ;    /* construct the shifted reversed seed */
	      p   = p & valid_base[base] ;
	    }
	  
	  /* import dw bases from the right exon */
	  cp = arrp (dna, a2, unsigned char) ;
	  for (int i = 0 ; i < dw ; i++, cp++)
	    {
	      w <<= 2 ; wr >>= 2 ;
	      p++ ;
	      unsigned int base = (*cp & 0x0F) ;
	      
	      w  |= w_table[base] ;    /* construct the shifted seed */
	      wr |= wr_table[base] ;    /* construct the shifted reversed seed */
	      p   = p & valid_base[base] ;
	    }
	  /* import 16 - dw  bases from the following exon */
	  if (dw < 16)
	    {
	      cp = arrp (dna, w2, unsigned char) ; /* a2 is G in bio coords, so cp is the first base of the second exon */
	      for (int i = 0 ; i < 16 - dw  ; i++, cp++)
		{
		  w <<= 2 ; wr >>= 2 ;
		  p++ ;
		  unsigned int base = (*cp & 0x0F) ;
		  
		  w  |= w_table[base] ;    /* construct the shifted seed */
		  wr |= wr_table[base] ;    /* construct the shifted reversed seed */
		  p   = p & valid_base[base] ;
		}
	    }
	  
	  if (p == 32)   /* 4,5,6,7,8,9,10,11 = 8 possibilites, codes on 3 bits 21-24 */
	    for (int i1 = 4 - (wLen - 16)  ; i1 < 13 ; i1++) /* search every word with at least 4 bases on each side */
	      {
		if (i1 < 0) continue ;
		unsigned long long x =   (w >> (2 * (32 - wLen - i1))) & maskSeedLn ;
		unsigned long long xr =  (wr >> 2 * (i1)) & maskSeedLn ;
		int dx1 = 16 - i1 ;     /* number of letters in first exon */
		BOOL  minus  = (x > xr) ; /* word is antistrand to the exon-exon junction */
		unsigned long long word = ( minus ? xr : x ) ;
		k = word & maskN ;   /* right bits are used to select the relevant table */
		word >>= nHidden2 ;  /* remove excessive bits, by construction of NN, they are stored in k */
		word &= mask32 ;     /* masked down to 32 bits integer */
		cw = bigArrayp (bbG->cwsN[k], bigArrayMax(bbG->cwsN[k]), CW) ;
		cw->nam = (chrom << 1) | (minus ? 0x1 : 0) ;
		cw->seed = word ;
		cw->pos = a1 ; /* in introns, a1 cannot be zero */
		cw->intron =                  /* bit 30 alone, isExon */
		  (0x1u << 31)                 /* bit 31: isIntron */
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

static void saCodeSequenceSeedsStep1 (const PP *pp, BB *bb)
{
  int NN = pp->nIndex ;
  BigArray cwsN[NN] ;
  CW *restrict cw  ;  
  const unsigned char *restrict cp ;
  Array dnas = bb->dnas ;
  int k, ia, iaMax = arrayMax (dnas), nSeeds = 0 ;
  long int dMax = bigArrayMax (bb->globalDna) / NN ;
  const int wLen = pp->seedLength ;
  const int nHidden = wLen > 16 ? wLen - 16 : 0 ;
  const int nHidden2 = nHidden << 1 ;
  const int shiftB = 64 - 2 * wLen ;
  const unsigned long long maskN = (unsigned long long)NN - 1 ;
  const unsigned long long mask32 = 0xffffffff ;
  const unsigned long long maskSeedLn = (1ULL << 2*wLen) - 1 ;
  long int nReads1 = 0 ;
  int minEntropy = pp->minEntropy ;
  int minLength = pp->minLength ? pp->minLength : pp->seedLength ;
  int minAli = pp->minAli ;

  memset (cwsN, 0, sizeof (cwsN)) ;
  if (dMax >> 31)
    messcrash ("saCodeSequenceSeedsStep1 ln=%ld > 2G", dMax) ;

  dMax *= 1.2 ;  /* avoid reallocation of marginal size fluctuations */
  for (k = 0 ; k < NN ; k++)
    cwsN[k] = bigArrayHandleCreate (dMax, CW, bb->h) ;

  bb->step = 1 ;
  for (ia = 1 ; ia < iaMax ; ia++)
    {
      Array dna = arr (dnas, ia, Array) ;
      if (! dna)
        continue ;
      int ii, jj, p = 0 ;
      int iMax = arrayMax (dna) ;
      unsigned long long w = 0, wr = 0 ;

      if (iMax < minLength)
        {
          bb->runStat.p.tooShort++ ;
          bb->runStat.p.tooShortBases += iMax ;
          continue ;
        }
      if (minEntropy > 0)
	{
	  int entropy = oligoEntropy (arrp (dna, 0, unsigned char), iMax, minEntropy) ;
	  if (entropy < minEntropy && entropy < iMax/2)
	    {
	      bb->runStat.p.lowEntropy++ ;
	      bb->runStat.p.lowEntropyBases += iMax ;
	      continue ;
	    }
	}

      cp = arrp (dna, 0, unsigned char) ;

      for (ii = 0, jj = -wLen + 1 ; ii < iMax ; ii++, jj++, cp++)
        {
          w <<= 2 ; wr >>= 2 ;
          p++ ;

          __builtin_prefetch(cp + 64, 0, 1) ;
          unsigned int base = (*cp & 0x0F) ;

          w  |= w_table[base] ;
          wr |= wr_table[base] ;
          p   = p & valid_base[base] ;

          if (p < wLen)
            continue ;

          unsigned long long x  =  w & maskSeedLn ;
          unsigned long long xr = (wr >> shiftB) & maskSeedLn ;
          BOOL minus = (x > xr) ;
          unsigned long long z  = minus ? xr : x ;

          /* skip homopolymers and homodimers */
          if (
              (((w >> 2) << 2) == (w << 2)) ||
              (((w >> 4) << 4) == (w << 4))
             )
            continue ;

          k = z & maskN ;
          z >>= nHidden2 ;
          z &= mask32 ;

          cw = bigArrayp (cwsN[k], bigArrayMax(cwsN[k]), CW) ;
          cw->nam    = (ia << 1) | (minus ? 0x1 : 0) ;
          cw->seed   = z ;
          cw->pos    = jj + 1 ;   /* bio coordinates of first base of seed */
          cw->intron = 0 ;
          nSeeds++ ;
        }
    }

  bb->cwsN = halloc (NN * sizeof(BigArray), bb->h) ;

  long int avDna = bb->runStat.p.nBase1 /(nReads1 + 0.0001) ;
  if (minAli <= 0)
    {
      minAli = avDna / 2 ;
      if (minAli > 30) minAli = 30 ;
    }
  bb->runStat.p.minAli = minAli ;

  if (pp->debug)
    fprintf (stderr, "Agent %d lane %d allocated %d seeds\n", bb->readerAgent, bb->lane, nSeeds) ;
  for (k = 0 ; k < NN ; k++)
    bb->cwsN[k] = cwsN[k] ;

  return ;
}  /* saCodeSequenceSeedsStep1 */

/**************************************************************/

static BOOL saReadClipRight (BB *bb, int read, Array dna, unsigned char adaptor[32], int mult)
{
  int xv, xv2 ;
  unsigned char *dna3 = (unsigned char *)strnew (arrp (dna, 0, char), 0) ;
  BOOL clipped = FALSE ;
  
  int n = arrayMax (dna) ;
  int n0 = ustrlen (adaptor) ;
  int n1 = n0 ;
  xv = 0 ;
  if (n1 > 12)
    {
      while ((xv2 = dnaPickMatch (dna3 + xv, n - xv, adaptor, 2, 2)))
	xv += xv2 ;
      while (!xv && n1 > 12) /* try to shorten the adaptor and reextend the hit */
	{
	  int n2 = n1 / 2, n3 ; 
	  int e1 = 2, e2 = 2 ;
	  char cc = adaptor[n2] ;
	  adaptor[n2] = 0 ;
	  while ((xv2 = dnaPickMatch (dna3 + xv, n - xv, adaptor, 2, 2)))
	    xv += xv2 ;
	  adaptor[n2] = cc ;
	  if (xv) /* verify that the end of the read matches the vector */
	    {
	      n3 = n - xv + 1;
	      cc = adaptor[n3] ;
	      adaptor[n3] = 0 ;
	      if (n - xv < 12) e1 = e2 = 0 ;
	      xv2 = dnaPickMatch (dna3 + xv - 1, n - xv + 1, adaptor, e1, e2) ;
	      adaptor[n3] = cc ;
	      if (! xv2)
		{ xv = 0 ; n1 = 0 ; } /* reject */
	    }
	  n1 /= 2 ;
	}
    }
  else
    {
      while ((xv2 = dnaPickMatch (dna3 + xv, n - xv, adaptor, 0, 0)))
	xv += xv2 ;
    }
  if (xv && n - xv >= 8)
    {
      /*
	sx->exitAdaptorSeqFound++ ;
	sx->exitAdaptorTagFound += sx->count ;
	sx->exitAdaptorBpFound += sx->count * (n - xv + 1) ;
	sx->exitAdaptorBpClipped += (n - xv + 1) ;
      */
      clipped = TRUE ;
      unsigned char *cp = arrp (dna, xv-1, unsigned char) ; /* first unaligned base */
      long int *overH = (read & 0x1 ? bb->runStat.overhangR2 : bb->runStat.overhangR1) ;

      int first = -1 ;
      int cc = *cp ;
      switch (cc)
	{
	case A_:
	  first = 0 ;
	  break ;
	case T_:
	  first = 1 ;
	  break ;
	case G_:
	  first = 2 ;
	  break ;
	case C_:
	  first = 3 ;
	  break ;
	default:
	  first = -1 ;
	}

      if (first >= 0)
	for (int i = 0 ; i < 30 && *cp ; i++, cp++)
	  {
	    cc = *cp ;
	    int ii = 150 * first + 5 * i ;
	    switch (cc)
	      {
	      case A_:
		overH[ii + 0]+= mult ;
		break ;
		case T_:
		  overH[ii + 1]+= mult ;
		  break ;
	      case G_:
		overH[ii + 2]+= mult ;
		break ;
	      case C_:
		overH[ii + 3]+= mult ;
		break ;
		default:
		  overH[ii + 4]+= mult ;
		  break ;
	      }
	  }

      cp = arrp (dna, xv-1, unsigned char) ;
      *cp = 0 ; arrayMax (dna) = xv-1 ;
    }

  ac_free (dna3) ;
  return clipped ;
} /* saReadClipRight */

/**************************************************************/

  static BOOL saReadClip (BB *bb, int read, Array dna, unsigned char adaptor[32], BOOL rightClip, int mult)
{
  if (rightClip)
    return saReadClipRight (bb, read, dna, adaptor, mult) ;
  else
    ;
  return FALSE ;
} /* saReadClip */

/**************************************************************/

void saCodeSequenceSeeds (const PP *pp, BB *bb, int step)
{
  if (step == -1)  /* the new code is slower in Elapsed  and User time (/usr/bin/time -f 'E %E U %U' */
    return  saCodeSequenceSeedsStep1 (pp, bb) ;
  int NN = pp->nIndex ;
  BigArray cwsN[NN] ;
  CW *restrict cw  ;  
  const unsigned char *restrict cp ;
  Array dnas = bb->dnas ;
  int k, ia, iaMax = arrayMax (dnas), nSeeds = 0 ;
  long int dMax = (bb->globalDna ? bigArrayMax (bb->globalDna) : bb->length) / (NN * step) ;
  const int wLen = pp->seedLength ;
  const int nHidden = wLen > 16 ? wLen - 16 : 0 ;
  const int nHidden2 = nHidden << 1 ;
  const int shiftB = 64 - 2 * wLen ;
  const  unsigned long long maskN = (unsigned long long)NN - 1 ;
  const unsigned long long mask32 = 0xffffffffULL ; /* 4 bytes integer */
  const unsigned long long maskSeedLn = (1ULL << 2*wLen) - 1 ;
  int iMaxMax = 0 ;
  long int nReads1 = 0 ;
  long int nReads2 = 0 ;
  int minEntropy = pp->minEntropy ;
  int minLength = pp->minLength ;
  ADAPTORS adaptors = {{0}} ;
  BOOL newMultiplicity = FALSE ;
  
  if (pp->fastc)
    bb->multiplicity = arrayHandleCreate (iaMax, int, bb->h) ;
  memset (cwsN, 0, sizeof (cwsN)) ;
  if (step < 1)
    messcrash ("codeWordsDo received step = %d < 1", step) ;
  if (dMax >> 31)
    messcrash ("codeWordsDo received ln=%lu, step=%d, ln/step=%ld > 2G", bb->length, step, dMax) ;

  if (NN > 1)  /* Avoid reallocation of marginal size fluctuations */
    dMax *= 1.2 ;
  for (k = 0 ; k < NN ; k++)
    cwsN[k] = bigArrayHandleCreate (dMax, CW, bb->h) ;

  if (1)
    while (! saSetGetAdaptors (0, &(bb->isRna), &adaptors, bb->run, &(bb->errCost), &(bb->runStat.p.minAli)))
      {
	if (bb->lane <= 1)
	{
	  bb->runStat.p.minAli = pp->seedLength ;
	  break ; /* pilot block proceeds */
	}
	sleep (1) ;                /* all other blocks wait */
      }
  
  bb->step = step ;
  bb->runStat.p.nBase1 = 0 ;
  bb->runStat.p.nBase2 = 0 ;
  bb->runStat.p.nClippedBase1 = 0 ;
  bb->runStat.p.nClippedBase2 = 0 ;
  bb->runStat.p.nReads = 0 ;
  bb->runStat.p.nPairs = 0 ;
  
  for (ia = 1 ; ia < iaMax ; ia++)
    {
      Array dna = arr (dnas, ia, Array) ;
      if (! dna)
	continue ;

      int mult = 1 ;
      if (bb->multiplicity)
	{
	  const char *cm0, *cm ;
	  cm0 = cm = dictName (bb->dict, ia>>1) ;
	  while (*cm) cm++ ;
	  while (cm > cm0 && *cm != '#') cm-- ;
	  if (cm > cm0 && *cm == '#')
	    mult = atoi (cm + 1) ;
	  if (mult <= 0) mult = 1 ;
	  array (bb->multiplicity, ia, int) = mult ;
	}

      int iMax = arrayMax (dna) ;
      bb->runStat.p.nBase1 += (ia & 0x1) ? 0 : iMax * mult ;
      bb->runStat.p.nBase2 += (ia & 0x1) ? iMax * mult : 0 ;
      if (ia & 0x1)
	{
	  unsigned char *cp = arrp (dna, 0, unsigned char) ;
	  for (int i = 0 ; i < iMax && i < LETTERMAX ; i++, cp++)
	    bb->runStat.p.letterProfile2[5*i + natgc[(int)*cp]]+= mult ;
	}
      else
	{
	  unsigned char *cp = arrp (dna, 0, unsigned char) ;
	  for (int i = 0 ; i < iMax && i < LETTERMAX ; i++, cp++)
	    bb->runStat.p.letterProfile1[5*i + natgc[(int)*cp]]+= mult ;
	}

      if (iMax >= minLength)
	{
	  if (ia & 0x1)
	    {
	      if (adaptors.a2R[0])
		bb->runStat.p.nClippedAdaptor2R += saReadClip (bb, ia, dna, adaptors.a2R, TRUE, mult) ? mult : 0 ;
	      if (adaptors.a2L[0])
		bb->runStat.p.nClippedAdaptor2L += saReadClip (bb, ia, dna, adaptors.a2L, FALSE, mult) ? mult : 0 ;
	    }
	  else
	    {
	      if (adaptors.a1R[0])
		bb->runStat.p.nClippedAdaptor1R += saReadClip (bb, ia, dna, adaptors.a1R, TRUE, mult) ? mult : 0 ;
	      if (adaptors.a1L[0])
		bb->runStat.p.nClippedAdaptor1L += saReadClip (bb, ia, dna, adaptors.a1L, FALSE, mult) ? mult : 0 ;
	    }
	}

      iMax = arrayMax (dna) ;
      iMaxMax = iMaxMax > iMax ? iMaxMax : iMax ;
      nReads1 += (ia & 0x1) ? 0 : mult ;
      nReads2 += (ia & 0x1) ? mult : 0 ;
      if (1)
	{
	  unsigned char *cp = arrp (dna, 0, unsigned char) ;
	  for (int i = 0 ; i < iMax ; i++, cp++)
	    bb->runStat.p.NATGC[natgc[(int)*cp]]++ ;
	}
      
      if (iMax < minLength) 
	{
	  bb->runStat.p.tooShort += mult ;
	  bb->runStat.p.tooShortBases += iMax * mult ;
	  arrayMax (dna) = 0 ;
	  continue ;
	}
      if (minEntropy > 0)
	{
	  int entropy = oligoEntropy (arrp (dna, 0, unsigned char), iMax, minEntropy) ;
	  if (entropy < minEntropy && entropy < iMax/2)
	    {
	      bb->runStat.p.lowEntropy += mult ;
	      bb->runStat.p.lowEntropyBases += iMax * mult ;
	      arrayMax (dna) = 0 ;
	      continue ;
	    }
	}

    }

  if (0)
    {
      RunSTAT up = bb->runStat ;
      fprintf (stderr, "lane=%d clipped=%ld adap=%d\n", bb->lane, up.p.nClippedAdaptor1R, adaptors.a1R[0]) ;
    }

  if (1 && ! bb->multiplicity && nReads2 == 0 && iMaxMax < 130 && ! pp->keepReadNames)
    {
      DICT *newDict = dictHandleCreate (nReads1, bb->h) ;
      /* remove dnas of read2 and empty dnas */
      int ii, jj ;
      char buf[64] ;
      bb->multiplicity = arrayHandleCreate (nReads1, int, bb->h) ;
      newMultiplicity = TRUE ;
      for (ii = jj = 0 ; ii < iaMax ; ii += 2)
	{
	  Array d = array (dnas, ii, Array) ;
	  if (!ii || (d && arrayMax (d)))
	    array (dnas, jj++, Array) = d ;
	}
      iaMax = arrayMax (dnas) = jj ;
      arraySort (dnas, dnaOrder) ;
      for (ii = jj = 0 ; ii < iaMax ; ii++)
	{
	  Array di = array (dnas, ii, Array) ;
	  int k = 1 ;
	  for (int jk = ii+1 ; jk < iaMax ; jk++)
	    {
	      Array dj = array (dnas, jk, Array) ;
	      if (0 || dnaOrder (&di, &dj))
		break ;
	      k++ ;
	    }
	  array (dnas, jj, Array) = di ;
	  sprintf (buf, "s.%d#%d", jj, k) ;
	  dictAdd (newDict, buf, 0) ;
	  array (bb->multiplicity, 2*jj, int) = k ;
	  array (bb->multiplicity, 2*jj+1, int) = k ;
	  jj++ ;
	  ii += k - 1 ;
	}
      bb->dict = newDict ;
      iaMax = jj ;
      for (ii = iaMax - 1 ; ii >= 0 ; ii--)
	{
	  Array d = array (dnas, ii, Array) ;
	  array (dnas, 2*ii + 0, Array) = d ;
	  array (dnas, 2*ii + 1, Array) = 0 ;
	}
      iaMax = arrayMax (dnas) = 2 * iaMax ;
    }

  bb->runStat.p.nReads = nReads1 + nReads2 ;
  bb->runStat.p.nPairs += nReads2 ;
  

  iaMax = arrayMax (dnas) ;
  for (ia = 1 ; ia < iaMax ; ia++)
    {
      Array dna = arr (dnas, ia, Array) ;
      if (! dna)
	continue ;

      int mult = 1 ;
      if (newMultiplicity)
	mult = array (bb->multiplicity, ia, int) ;
      else if (bb->multiplicity)
	{
	  const char *cm0, *cm ;
	  cm0 = cm = dictName (bb->dict, ia>>1) ;
	  while (*cm) cm++ ;
	  while (cm > cm0 && *cm != '#') cm-- ;
	  if (cm > cm0 && *cm == '#')
	    mult = atoi (cm + 1) ;
	  if (mult <= 0) mult = 1 ;
	  array (bb->multiplicity, ia, int) = mult ;
	}

      int ii, jj, p = 0, cStep, robin ;
      int iMax = arrayMax (dna) ;
      unsigned long long  w, wr ;
      if (1) step = iMax < 60 ? 1 : bb->step ;

      bb->runStat.p.nClippedBase1 += (ia & 0x1) ? 0 : iMax * mult ;
      bb->runStat.p.nClippedBase2 += (ia & 0x1) ? iMax * mult : 0 ;
      if (ia & 0x1)
	{
	  unsigned char *cp = arrp (dna, 0, unsigned char) ;
	  for (int i = 0 ; i < iMax && i < LETTERMAX ; i++, cp++)
	    bb->runStat.p.clippedLetterProfile2[5*i + natgc[(int)*cp]]+= mult ;
	}
      else
	{
	  unsigned char *cp = arrp (dna, 0, unsigned char) ;
	  for (int i = 0 ; i < iMax && i < LETTERMAX ; i++, cp++)
	    bb->runStat.p.clippedLetterProfile1[5*i + natgc[(int)*cp]]+= mult ;
	}
      
      if (0)
	{
	  if (iMax < minLength) 
	    {
	      bb->runStat.p.tooShort += mult ;
	      bb->runStat.p.tooShortBases += iMax * mult ;
	      continue ;
	    }
	  if (minEntropy > 0)
	    {
	      int entropy = oligoEntropy (arrp (dna, 0, unsigned char), iMax, minEntropy) ;
	      if (entropy < minEntropy && entropy < iMax/2)
		{
		  bb->runStat.p.lowEntropy += mult ;
		  bb->runStat.p.lowEntropyBases += iMax * mult ;
		  continue ;
		}
	    }
	}
      w = wr = 0 ;
      /* each q vector is used round-robin style, so no copying is needed */
      BOOL qminus[step] ;
      unsigned long long qw[step], qwr[step], qx[step], qxr[step], qz[step], qhz[step] ;
      memset (qminus, 0, sizeof (qminus)) ;
      memset (qw, 0, sizeof (qw)) ;
      memset (qwr, 0, sizeof (qwr)) ;
      memset (qx, 0, sizeof (qx)) ;
      memset (qxr, 0, sizeof (qxr)) ;
      memset (qz, 0, sizeof (qz)) ;
      memset (qhz, 0, sizeof (qhz)) ;

      cp = arrp (dna, 0, unsigned char) ;

      for (ii = 0, jj = -wLen + 1, cStep = -wLen + 1, p = 0, robin = 0 ; ii < iMax ; ii++, jj++, cp++, cStep++, robin = (robin + 1 < step ? robin + 1 : 0))
	{
	  /* construct a 2*wLen bits long integer representing wLen bases (A_, T_, G_, C_) using 2 bits per base */
	  w <<= 2 ; wr >>= 2 ;
	  p++ ;

	  __builtin_prefetch(cp + 64, 0, 1) ;
	  unsigned int base = (*cp & 0x0F) ;
	  
	  w  |= w_table[base] ;    /* construct the shifted seed */
	  wr |= wr_table[base] ;    /* construct the shifted reversed seed */
	  p   = p & valid_base[base] ;

	  qw[robin] = w ;
	  qwr[robin] = wr ;
	  qx[robin] = qw[robin] & maskSeedLn ;
	  qxr[robin] =  (qwr[robin] >> shiftB) & maskSeedLn ;
	  qminus[robin]  = (qx[robin] > qxr[robin]) ;
	  qz[robin] = qminus[robin] ? qxr[robin] : qx[robin] ;
	  qhz[robin] = (qz[robin] ^ (qz[robin] << 14)) ^ (qz[robin] << 23) ;
	  /*  possibly better but incompatible: try
	    unsigned long long h = qz[robin] ;
	    h ^= h << 14 ; h ^= h >> 7 ; h ^= h << 23 ;
	    qhz[robin] = h ;
	  */
	  if (
	      p < wLen ||
	      (((  w >> 2) << 2) == (w << 2) ) || /* avoid homopolymers */
	      (((  w >> 4) << 4) == (w << 4) )    /* avoid homodimers */
	      )
	    qz[robin] = 0 ;

	   
	  if (cStep == step) cStep = 0 ;
	  if (!cStep)
	    {
	      unsigned long long word = 0 ;
	      int dx = 0 ;
	      BOOL minus, qok = FALSE ;
	      
	      for (int i = 0 ; !qok && i < step ; i++)
		qok = (qz[i] != 0) ;
	      if (! qok)
		{ cStep = -1 ; continue ; } /* we hope to export the next word */
	      unsigned long long besthz = qhz[robin] ;
	      int bestk = robin ; dx = 0 ;
	      int myMax = (step > jj ? jj : step) ;
	      for (int i = 1 ; i < myMax ; i++)
		{
		  int k = robin - i ; k = (k < 0 ? k + step : k) ; /* was:  k = (robin - i + step) % step ; */
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
		fprintf (stderr, "codeWordsDo p=%d k=%d pos=%u seed=%u wLen=%d\n", p, k, cw->pos,cw->seed, wLen) ;
	    }
	}
    }
  bb->cwsN = halloc (NN * sizeof(BigArray), bb->h) ;


#ifdef JUNK
  long int avDna = bb->runStat.p.nBase1 /(nReads1 + 0.0001) ;
  if (minAli <= 0)
    {
      minAli = avDna / 2 ;
      if (minAli > 30) minAli = 30 ;
    }
  bb->runStat.p.minAli = minAli ;
#endif
  // fprintf (stderr, "minEntropy=%d minAli=%d avDna=%ld\n", minEntropy, minAli, avDna) ;

  if (pp->debug)
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

