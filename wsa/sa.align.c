/*
 * sa.align.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module extends seed hits to full alignments
 */
/*
  #define ARRAY_CHECK
*/

#include "sa.h"
#define MAXJUMP 3
#define MAXJUMP2 8

/**************************************************************/
/**************************************************************/
/* a0 = a1 - x1 is the putative position of base 1 of the read 
 * It also works for the negative strand (a1 < 0, x1 > 0).
 */
static int countChromOrder (const void *va, const void *vb)
{
  const COUNTCHROM *up = va ;
  const COUNTCHROM *vp = vb ;
  int n ;
  n = up->weight - vp->weight ; if (n) return -n ;
  n= up->seeds - vp->seeds ; if (n) return -n ;
  n= up->chrom - vp->chrom ; if (n) return n ;
  n= up->a1 < up->a2 ? 1 : -1 ;
  
  return n ;
} /* countChromOrder */

/**************************************************************/

int saAlignOrder (const void *va, const void *vb)
{
  const ALIGN *up = va ;
  const ALIGN *vp = vb ;
  int n ;
  n = up->read - vp->read ; if (n) return n ;
  n = up->chainScore - vp->chainScore ; if (n) return -n ;
  n = up->chain - vp->chain ; if (n) return n ;
  n = up->chrom - vp->chrom ; if (n) return n ;
  n = up->x1 - vp->x1  ; if (n) return n ;
  n = up->x2 - vp->x2  ; if (n) return -n ;
  n = up->nErr - vp->nErr ; if (n) return n ;
  n = up->a1 - vp->a1  ; if (n) return n ;
  n = up->a2 - vp->a2  ; if (n) return n ;
  n = up->nTargetRepeats - vp->nTargetRepeats ; if (n) return 1 ;
  return 0 ;
} /* saAlignOrder */

/**************************************************************/
/* a0 = a1 - x1 is the putative position of base 1 of the read 
 * It also works for the negative strand (a1 < 0, x1 > 0).
 */
static int hitReadPosOrder (const void *va, const void *vb)
{
  const HIT *up = va ;
  const HIT *vp = vb ;
  int n, n1, n2 ;

  n = ((up->read > vp->read) - (up->read < vp->read)) ; if (n) return n ;
  n = ((up->chrom > vp->chrom) - (up->chrom < vp->chrom)) ; if (n) return n ; 
  n = ((up->x1 > vp->x1) - (up->x1 < vp->x1)) ; if (n) return n ;
  n1 = up->a1 + (up->x1 >> NSHIFTEDTARGETREPEATBITS) ;
  n2 = vp->a1 + (vp->x1 >> NSHIFTEDTARGETREPEATBITS) ;
  n = n1 - n2 ;
  return n ;
} /* hitReadPosOrder */

/**************************************************************/

static void showCountChroms (Array countChroms)
{
  int ii, iMax = countChroms ? arrayMax (countChroms) : 0 ;
  for (ii = 0 ; ii < iMax && ii < 3000 ; ii++)
    {
      COUNTCHROM *zp = arrp (countChroms, ii, COUNTCHROM) ;
      if (zp->chrom)
	printf (".. %d:  weight=%.2f\tseeds=%d\t1:%d/2:%d/4:%d/8:%d/16:%d/32:%d\tindex %d %d\tchrom=%d\tpos=%u %u\tda=%u\tx= %u %u\n"
		, ii, zp->weight, zp->seeds
		, zp->seed1
		, zp->seed2
		, zp->seed4
		, zp->seed8
		, zp->seed16
		, zp->seed32
		, zp->i1, zp->i2
		, zp->chrom
		, zp->a1, zp->a2
		, zp->a2 > zp->a1 ? zp->a2 - zp->a1 : zp->a1 - zp->a2
		, zp->x1, zp->x2
		) ;
    }
} /* showCountChroms */

/**************************************************************/
/* merge align->errors into bb->errors */
/* sorry, we use this extremelly obfuscated way to encode the errors
 * rather than a simple clear utilisation of vp->type + actual bases
 * just to maintain a small max value for the variable type.
 * The intention is to help the CPU cache by maintaining a very small aa array
 *
 * There is also a complex way to choose the strand
 * If we use paired end sequencing, we want to switch the strand of read2
 * In addition, if we sequence RNA (as opposed to DNA) and
 * we map antistrand to a known gene, we flip read1 and again read2
 */
static void mergeErrors (Array aa, Array bb, unsigned int flip)
{
  int i, iMax = arrayMax (bb) ;
  A_ERR *vp = 0 ;

  if (iMax)
    for (i = 0, vp =arrp (bb, 0, A_ERR) ; i < iMax ; i++, vp++)
      {
	int type = 0 ;
	switch (vp->type)
	  {
	  case ERREUR : /* use 4 right bits type < 15 */
	    switch (vp->baseLong)
	      {
	      case A_: break ;
	      case T_: type = 0xc ;break ;
	      case G_: type = 0x4 ; break ;
	      case C_: type = 0x8 ; break ;
	      }
	    switch (vp->baseShort)
	      {
	      case A_: break ;
	      case T_: type |= 0x3 ;break ;
	      case G_: type |= 0x1 ; break ;
	      case C_: type |= 0x2 ; break ;
	      }
	    type = 0xf & (type ^ flip) ; /* flip the last 4 bits : A<->T   G<->C */
	    break ;
	  case INSERTION:
	    switch (vp->baseShort)
	      {
	      case A_: type = 0x0 ; break ;
	      case T_: type = 0x3 ; break ;
	      case G_: type = 0x1 ; break ;
	      case C_: type = 0x2 ; break ;
	      }
	    type = 0x3 & (type ^ flip) ; /* flip the last 4 bits : A<->T   G<->C */
	    type |= 0x10 ;
	    break ;
	  case TROU:
	    switch (vp->baseLong)
	      {
	      case A_: type = 0x0 ; break ;
	      case T_: type = 0xc ;break ;
	      case G_: type = 0x4 ; break ;
	      case C_: type = 0x8 ; break ;
	      }
	    type = 0xc & (type ^ flip) ; /* flip the last 4 bits : A<->T   G<->C */
	    type |= 0x20 ;
	    break ;
	  case INSERTION_DOUBLE:
	    type = 0x32 ;
	    break ;
	  case INSERTION_TRIPLE:
	    type = 0x33 ;
	    break ;
	  case TROU_DOUBLE:
	    type = 0x38 ;
	    break ;
	  case TROU_TRIPLE:
	    type = 0x3c ;
	    break ;
	  default:
	    break ;
	  }
	if (type)
	  (array (aa, type, int)) += 1 ;
      }
  return ;
} /* mergeErrors */

/**************************************************************/

static BOOL alignExtendHit (Array dna, Array dnaG, Array dnaGR, Array err
			    , BOOL isDown, int chromLength
			    , int *a1p, int *a2p, int *x1p, int *x2p
			    , int errCost
			    , BOOL isIntron
			    , int errMax, int minAli, int maxJump
			    )
{
  int nN = 0, dx = 0, nerr = 0 ;
  int x1 = *x1p, x2 = *x2p ;
  int a1 = *a1p, a2 = *a2p ;
  arrayMax (err) = 0 ;

  errMax = isIntron ? 0 : errMax ;
  
  

  if (! isDown)
    {
      Array dummy = dnaG ; dnaG = dnaGR ; dnaGR = dummy ;
      a1 = chromLength - a1 + 1 ; a2 = chromLength - a2 + 1 ;
      *a1p = a1 ; *a2p = a2 ;
    }

  arrayMax (err) = 0 ;
  aceDnaDoubleTrackErrors (dna, &x1, &x2, TRUE /* bio coordinates, extend = TRUE */
			   , dnaG, dnaGR, &a1, &a2
			   , &nN, err, maxJump, errMax, TRUE, 0) ;
  if (x1 > *x1p || x2 < *x2p)
    return FALSE ;
  dx = x2 - x1 + 1 ;
  if (*a1p < *a2p && (a1 > a2 || a1 > *a1p || a2 < *a2p))
    return FALSE ;
  if (*a1p > *a2p && (a1 < a2 || a1 < *a1p || a2 > *a2p))
    return FALSE ;
  /* reclip */
  /* clip left errors */
  nerr = arrayMax (err) ;
  if (nerr)
    {
      int i, j, y1 = x1 - 1, y2 = x2 - 1 ;  /* natural coordinates */
      A_ERR *up = arrp (err, 0, A_ERR) ;
      A_ERR *vp = up ;
      BOOL wonder = TRUE ;
      
      for (i = j = 0 ; i < nerr ; i++, up++)
	{
	  if (wonder && errCost > (up->iShort - y1)) /* clip this error */
	    {
	      x1 = up->iShort + 1 ; /* +1 for bio coordinates */
	      a1 = up->iLong + 1 ;  /* +1 for bio coordinates */
	      switch (up->type)
		{
		case AMBIGUE: break ;
		case INSERTION_TRIPLE: x1 += 3 ; break ;
		case TROU_TRIPLE: a1 += 3 ;break ;
		case INSERTION_DOUBLE: x1 += 2 ; break ;
		case TROU_DOUBLE: a1 += 2 ; break ;
		case INSERTION: x1++ ; break ;
		case  TROU: a1++ ; break ;
		default : a1++ ; x1++ ; break ;
		}
	      y1 = x1 - 1 ; /* natural coordinates */
	      continue ;
	    }
	  else  /* register all remaining errors */
	    {	      
	      wonder = FALSE ;
	      if (j < i) *vp = *up ;
	      j++ ; vp++ ;
	    }
	}
      nerr =  arrayMax (err) = j ; /* some 5' errors are dropped */
	    
      wonder = TRUE ;
      if (nerr)
	for (i = nerr - 1, up = arrp (err, i, A_ERR) ; i >= 0 ; i--, up--)
	  {
	    if (wonder && errCost > y2 - up->iShort) /* clip this error */
	      {
		nerr-- ;
		x2 = up->iShort - 1 + 1 ; /* +1 for bio coordinates */
		a2 = up->iLong  - 1 + 1 ; /* +1 for bio coordinates */
		y2 = x2 - 1 ; /* natural coordinates */
		continue ;
	      }
	  }
      arrayMax(err) = nerr ; /* some 3' errors are dropped */
    }
  
  if (x1 > *x1p || x2 < *x2p)
    return FALSE ;
  if (*a1p < *a2p && (a1 > a2 || a1 > *a1p || a2 < *a2p))
    return FALSE ;
  if (*a1p > *a2p && (a1 < a2 || a1 < *a1p || a2 > *a2p))
    return FALSE ;

  if (! isDown)
    {
      a1 = chromLength - a1 + 1 ; a2 = chromLength - a2 + 1 ;
    }
  
  *x1p = x1 ; *x2p = x2 ;
  *a1p = a1 ; *a2p = a2 ;
  dx = x2 - x1 + 1 ;

  minAli = isIntron ? 8 : 22 ;
  if (dx < minAli)
    return FALSE ;
  if (errMax >= 0 &&  nerr > errMax)
    return FALSE ;
  if (dx < nerr * errCost)
    return FALSE ;

  return TRUE ;
} /* alignExtendHit */

/**************************************************************/
/**************************************************************/

static int alignFormatLeftOverhang (const PP *pp, BB *bb, ALIGN *up, Array dna, Array dnaG, Array dnaGR, ADAPTORS *adaptors, int i0, int iMax)
{
  int x1 = up->x1 ;
  int slBonus = 4 ;
  int leftClip = up->leftClip ;  /* already initialized to jump5 */
  BOOL doUpper = FALSE ;
  
  if (x1 > 1)
    {
      int nT = 0 ;
      int dx = x1 - 1 ;
      int Dx = dx ;
      char buf[32] ;
      memset (buf, 0, 32) ;
      if (dx > 30) { dx = 30 ; }
      char *cp = arrp (dna, x1 - 2, char) ; /* last unaligned base */
      long int *overH = (up->read & 0x1 ? bb->runStat.overhangL2 : bb->runStat.overhangL1) ;
      if (up->targetClass == 'G')
	{
	  int first = -1 ;
	  int cc = complementBase(*cp) ;
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
	      goto done ;
	    }

	  /* search for a polyT */
	  for (int i = 0 ; i < dx ; i++, cp--)
	    {
	      cc = complementBase(*cp) ;
	      buf[i] = cc ;
	      switch (cc)
		{
		case A_:
		  nT++ ;
		  break ;
		}
	    }
	  
	  /* may be we should look for the adaptor before Dx - dx */
	  if (bb->isRna >= 0 && nT >= 7 && 10 * nT >= 9 * dx && Dx < dx + 10)
	    {
	      POLYA *zp = arrayp (bb->confirmedPolyAs, arrayMax (bb->confirmedPolyAs), POLYA) ;
	      zp->chrom = up->chrom ^ 0x1 ;
	      zp->a1 = up->a1 - (up->chrom & 0x1 ? -1 : 1) ;
	      zp->run = bb->run ;
	      zp->n = 1 ;
	      doUpper = TRUE ;
	      bb->runStat.nClippedPolyT++ ;
	      leftClip = x1 - 1 ;
	      goto done ; /* do not polute the search for the adaptor with the presence of a polyT */
	    }
	  
	  cp = arrp (dna, x1 - 2, char) ; /* last unaligned base */
	  for (int i = 0 ; i < dx ; i++, cp--)
	    {
	      cc = complementBase(*cp) ;
	      int ii = 150 * first + 5 * i ;
	      switch (cc)
		{
		case A_:
		  overH[ii + 0]++ ;
		  break ;
		case T_:
		  overH[ii + 1]++ ;
		  break ;
		case G_:
		  overH[ii + 2]++ ;
		  break ;
		case C_:
		  overH[ii + 3]++ ;
		  break ;
		default:
		  overH[ii + 4]++ ;
		  break ;
		}
	    }

	  /* try to recognize the SLs */
	  if (bb->isRna >= 0 && pp->isWorm && pp->SLs)
	    {
	      int iSlMax = arrayMax (pp->SLs) ;
	      for (int iSl = 0 ; iSl < iSlMax ; iSl++)
		{	      
		  const char *sl = array (pp->SLs, iSl, char *) ;
		  int iMax = strlen (sl) ;
		  for (int di = 0 ; di < 6 ; di++)
		    {
		      int i, n = 0 ;
		      if (iMax > 30) iMax = 30 ;
		      
		      cp = arrp (dna, x1 - 2 + di, char) ; /* x1 - 2 =  first unaligned base */
		      for (i = 0 ; i < dx + di && i < 30 && i < iMax ; i++)
			n +=  (cp[-i] == sl[iMax - i - 1] ? 1 : 0) ;
		      if (n > 6 && 10 * n >= 9 * i)
			{
#ifdef JUNK			  
			  /* check that the prefix including an extra donor 'gt' 
			   * is not present in the genome
			   */
			  int kkk = ccp ? strlen ((const char*)ccp) : 0 ;
			  BOOL inGenome = FALSE ;
			  if (kkk < OVLN)
			    {
			      sprintf ((char *)buf2, "%s%c%c", ccp, G_, T_ | C_) ;
			      inGenome = dnaPickMatch (genome, genomeLn, buf2, 0, 0) ;
			    }
#endif
			  doUpper = TRUE ;
			  bb->runStat.nClippedSls[iSl]++ ;
			  leftClip = x1 - 1 + di ;
			  up->x1 += di ;
			  up->a1 += (up->a1 < up->a2 ? di :- di) ;
			  for (ALIGN *vp = up ; vp->read == up->read && vp->chain == up->chain ; vp++)
			    {
			      vp->score += slBonus - di ;
			      vp->chainScore += slBonus - di ;
			      vp->chainX1 += di ;
			      vp->chainA1 += (up->a1 < up->a2 ? di :- di) ;
			    }
			  SLS *zp = arrayp (bb->confirmedSLs, arrayMax (bb->confirmedSLs), SLS) ;
			  zp->chrom = up->chrom ;
			  zp->a1 = up->a1 ;
			  zp->run = (bb->run << 4) | (iSl & 0xf) ;
			  zp->n = 1 ;
			  goto done ;
			}
		    }	
		}
	    }
	  
	  /* try to recognize the known adaptor */
	  const char *adaptor = (up->read & 0x1 ? adaptors->a2L : adaptors->a1L) ;
	  if (adaptor[0])
	    {
	      int jj, iAmax = 1 ; /* arrayMax (aa) ; */
	      for (int iA = 0 ; iA < iAmax ; iA++)
		{	      
		  int iMax = strlen (adaptor) ;
		  for (int di = 0 ; di < 6 ; di++)
		    {
		      int i, n = 0 ;
		      if (iMax > 30) iMax = 30 ;
		      
		      cp = arrp (dna, x1 - 2 + di, char) ; /* x1 - 2 =  first unaligned base */
		      for (i = 0 ; i < dx + di && i < 30 && i < iMax ; i++)
			n +=  (cp[-i] == complementBase(adaptor[i]) ? 1 : 0) ;
		      if (n > 6 && 10 * n >= 9 * i)
			{
			  ALIGN *vp ;
			  
			  doUpper = TRUE ;
			  if (up->read & 0x1) 
			    bb->runStat.nClippedAdaptor2L++ ;
			  else
			    bb->runStat.nClippedAdaptor1L++ ;
			  leftClip = x1 - 1 + di ;
			  up->x1 += di ;
			  up->ali -= di ;
			  up->a1 += (up->a1 < up->a2 ? di :- di) ;
			  cp = arrp (dna, up->x1 - 1, char) ;
			  if (1 && di > 0) /* extend the inaligned buffer */
			    {
			      for (int i = 30 ; i >= di ; i--)
				buf[i] = buf[i- di] ;
			      for (int i = 0 ; i < di ; i++)
				buf[i] = complementBase(cp[-i - 1]) ; 
			      buf[31] = 0 ;
			    }
			    
			  for (vp = up, jj = i0 ; jj < iMax && vp->read == up->read && vp->chain == up->chain ; jj++, vp++)
			    {
			      vp->score -= di ;
			      vp->chainScore -= di ;
			      vp->chainAli -= di ;
			      vp->chainX1 += di ;
			      vp->chainA1 += (up->a1 < up->a2 ? di :- di) ;
			    }
			  goto done ;
			}
		    }	
		}
	    }
	}
    done:
      for (int i = 0 ; i < dx ; i++)
	buf[i] = dnaDecodeChar[(int)buf[i]] ;
      if (doUpper)
	bufferToUpper (buf) ;
      if (strlen (buf) > 0)
	  dictAdd (bb->dict, buf, &up->leftOverhang) ;
    }
  up->leftClip = leftClip ;
  return leftClip ;
} /* alignFormatLeftOverhang */

/**************************************************************/

static int alignFormatRightOverhang (const PP *pp, BB *bb, ALIGN *up, Array dna, Array dnaG, Array dnaGR, ADAPTORS *adaptors, int ii0, int iiMax)
{
  int x2 = up->x2 ;
  int ln = arrayMax (dna) ; /* length to align */
  int rightClip = ln ;
  BOOL doUpper = FALSE ;
  int slBonus = 4 ;
  
  if (x2 < ln)
    {
      int cc, nA = 0 ;
      int dx = ln - x2 ;
      int Dx = dx ;
      char buf[32] ;
      memset (buf, 0, 32) ;
      if (dx > 30) { dx = 30 ; }
      char *cp = arrp (dna, x2, char) ; /* first unaligned base */
      long int *overH = (up->read & 0x1 ? bb->runStat.overhangR2 : bb->runStat.overhangR1) ;
      
      for (int i = 0 ; i < dx ; i++, cp++)
	{
	  cc = *cp ;
	  buf[i] = cc ;
	  switch (cc)
	    {
	    case A_:
	      nA++ ;
	      break ;
	    }
	}
      
      if (ii0 == iiMax - 1 || up[1].read != up[0].read)
	/* isLast, proceed and look for polyA and adaptors */ ;
      else
	{ /* even if we change chain, we do not want i think to seach a polyA or adaptor */
	  if (dx > up[1].x1 - up[0].x2 - 1)
	    {
	      dx = up[1].x1 - up[0].x2 - 1 ;
	      if (dx >= 0 && dx < 32)
		buf[dx] = 0 ;
	    }	  
	  goto done ;  
	}

      if (1)
	{
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
	      goto done ;
	    }


	  /* search for a polyA */
	  if (up->targetClass == 'G')
	    {
	      /* may be we should look for the adaptor before Dx - dx */
	      if (bb->isRna >= 0 && nA >= 7 && 10 * nA >= 9 * dx && Dx < dx + 10)
		{
		  POLYA *zp = arrayp (bb->confirmedPolyAs, arrayMax (bb->confirmedPolyAs), POLYA) ;
		  zp->chrom = up->chrom ;
		  zp->a1 = up->a2 + (up->chrom & 0x1 ? -1 : 1) ;
		  zp->run = bb->run ;
		  zp->n = 1 ;
		  doUpper = TRUE ;
		  bb->runStat.nClippedPolyA++ ;
		  rightClip = x2 ;
		  goto done ; /* do not polute the search for the adaptor with the presence of a polyA */
		}
	    }
	  
	  cp = arrp (dna, x2, char) ; /* first unaligned base */
	  for (int i = 0 ; i < dx ; i++, cp++)
	    {
	      cc = *cp ;
	      int ii = 150 * first + 5 * i ;
	      switch (cc)
		{
		case A_:
		  overH[ii + 0]++ ;
		  break ;
		case T_:
		  overH[ii + 1]++ ;
		  break ;
		case G_:
		  overH[ii + 2]++ ;
		  break ;
		case C_:
		  overH[ii + 3]++ ;
		  break ;
		default:
		  overH[ii + 4]++ ;
		  break ;
		}
	    }
	  
	  /* try to recognize the SLs */
	  if (bb->isRna >= 0 && pp->isWorm && pp->SLs)
	    {
	      int iSlMax = arrayMax (pp->SLs) ;
	      for (int iSl = 0 ; iSl < iSlMax ; iSl++)
		{	      
		  const char *sl = array (pp->SLs, iSl, char *) ;
		  int iMax = strlen (sl) ;
		  for (int di = 0 ; di < 6 ; di++)
		    {
		      ALIGN *vp  ;
		      int i, n = 0 ;
		      if (iMax > 30) iMax = 30 ;
		      
		      cp = arrp (dna, x2 - di, char) ; /* x2 =  first unaligned base */
		      for (i = 0 ; i < dx + di && i < iMax ; i++)
			n +=  (cp[i] == sl[i] ? 1 : 0) ;
		      if (n > 6 && 10 * n >= 9 * i)
			{
#ifdef JUNK			  
			  /* check that the prefix including an extra donor 'gt' 
			   * is not present in the genome
			   */
			  int kkk = ccp ? strlen ((const char*)ccp) : 0 ;
			  BOOL inGenome = FALSE ;
			  if (kkk < OVLN)
			    {
			      sprintf ((char *)buf2, "%s%c%c", ccp, G_, T_ | C_) ;
			      inGenome = dnaPickMatch (genome, genomeLn, buf2, 0, 0) ;
			    }
#endif
			  doUpper = TRUE ;
			  bb->runStat.nClippedSls[iSl]++ ;
			  rightClip = x2 + 1 - di ;
			  up->x2 -= di ;
			  up->a2 -= (up->a1 < up->a2 ? di :- di) ;
			  for (vp = up ; vp->read == up->read && vp->chain == up->chain ; vp--)
			    {
			      vp->score += slBonus - di ;
			      vp->chainScore += slBonus - di ;
			      vp->chainX2 -= di ;
			      vp->chainA2 -= (up->a1 < up->a2 ? di :- di) ;
			    }
			  SLS *zp = arrayp (bb->confirmedSLs, arrayMax (bb->confirmedSLs), SLS) ;
			  zp->chrom = up->chrom ;
			  zp->a1 = up->a2 ;
			  zp->run = (bb->run << 4) | (iSl & 0xf) ;
			  zp->n = 1 ;
			  
			  goto done ;
			}
		    }	
		}
	    }
	}
	      
      /* try to recognize the known adaptor */
      const char *adaptor = (up->read & 0x1 ? adaptors->a2R : adaptors->a1R) ;
      if (adaptor[0])
	{
	  int iAmax = 1 ;
	  for (int iA = 0 ; iA < iAmax ; iA++)
	    {	      
	      int iMax = strlen (adaptor) ;
	      for (int di = 0 ; di < 6 ; di++)
		{
		  int jj, i, n = 0 ;
		  if (iMax > 30) iMax = 30 ;
		  
		  cp = arrp (dna, x2 - di, char) ; /* x2 =  first unaligned base */
		  for (i = 0 ; i < dx + di && i < iMax ; i++)
		    n +=  (cp[i] == adaptor[i] ? 1 : 0) ;
		  if (n > 6 && 10 * n >= 9 * i)
		    {
		      ALIGN *vp ;
		      doUpper = TRUE ;
		      if (up->read & 0x1) 
			bb->runStat.nClippedAdaptor2R++ ;
		      else
			bb->runStat.nClippedAdaptor1R++ ;
		      rightClip = x2 - di ;
		      up->x2 -= di ;
		      up->ali -= di ;
		      up->a2 -= (up->a1 < up->a2 ? di :- di) ;
		      cp = arrp (dna, up->x2, char) ;
		      if (1 && di > 0) /* extend the unaligned buffer */
			{
			  for (int i = 30 ; i >= di ; i--)
			    buf[i] = buf[i- di] ;
			  for (int i = 0 ; i < di ; i++)
			    buf[i] = cp[i] ;
			  buf[31] = 0 ;
			}
		      
		      for (vp = up, jj = ii0 ; jj >= 0 && vp->read == up->read && vp->chain == up->chain ; jj--, vp--)
			{
			  vp->score -= di ;
			  vp->chainScore -= di ;
			  vp->chainAli -= di ;
			  vp->chainX2 -= di ;
			  vp->chainA2 -= (up->a1 < up->a2 ? di :- di) ;
			}
		      goto done ;
		    }
		}	
	    }
	}
    
    done:
      for (int i = 0 ; i < dx ; i++)
	buf[i] = dnaDecodeChar[(int)buf[i]] ;
      if (doUpper)
	bufferToUpper (buf) ;
      if (strlen (buf) > 0)
	dictAdd (bb->dict, buf, &up->rightOverhang) ;
    }
  
  up->rightClip = rightClip ;
  return rightClip ;
} /* alignFormatRightOverhang */

/**************************************************************/

static void alignFormatErrors (const PP *pp, BB *bb, ALIGN *up, Array dna, Array dnaG, Array dnaGR, int read)
{
 A_ERR *ep = arrp (up->errors, 0, A_ERR) ;
  int ii, nerr = arrayMax (up->errors) ;
  vTXT txt1 = bb->txt1 ;
  vTXT txt2 = bb->txt2 ;
  char *sep ;
  int xShort, xLong  ;
  BOOL isUp = (up->a1 > up->a2) ;
  const BOOL debug = FALSE ;

  if (debug)
    aceDnaShowErr (up->errors) ;
  vtxtClear (txt1) ;
  vtxtClear (txt2) ;

  for (ii = 0, sep = "" ; ii < nerr ; ii++, ep++, sep = ",")
   {
      xShort = ep->iShort + 1 ;
      xLong = ep->iLong + 1 ;

      if (xShort <= 0 || xLong <= 0 || xShort >= arrayMax (dna))
	continue ;
      switch (ep->type)
	{
	  case TYPE80:
	    {
	      char cc1a, cc2a ;
	      
	      cc1a = arr (dnaG, (isUp ? xLong - 0 : xLong - 2), unsigned char) ;
	      cc2a = arr (dnaG, (isUp ? xLong - 1 : xLong - 1), unsigned char) ;
	      
	      vtxtPrintf (txt1,"%s%d:%c%c>oo"
			  , sep
			  , xShort - 1
			  , isUp ? dnaDecodeChar[(int)cc2a] : dnaDecodeChar[(int)cc1a]
			  , isUp ? dnaDecodeChar[(int)cc1a] : dnaDecodeChar[(int)cc2a]
			  ) ;
	      vtxtPrintf (txt2,"%s%d:%c%c>oo"
			  , sep
			  , isUp ? xLong - 1 : xLong - 1
			  , isUp ? dnaDecodeChar[(int)complementBase(cc1a)] : dnaDecodeChar[(int)cc1a]
			  , isUp ? dnaDecodeChar[(int)complementBase(cc2a)] : dnaDecodeChar[(int)cc2a]
			  ) ;
	      
	    }
	    break ;
	    
	case AMBIGUE:
	case ERREUR:
	  {
	    char ccS, ccL, ccSR, ccLR ;
	    /* int xLongR = arrayMax (dnaG) - xLong + 1 ; */
	    
	    ccS = ep->baseShort ;
	    ccL = ep->baseLong ;
	    
	    ccSR = isUp ? complementBase(ccS) : ccS ; 
	    ccLR = isUp ? complementBase(ccL) : ccL ; 

	    ccS = dnaDecodeChar[(int)ccS] ;
	    ccL = dnaDecodeChar[(int)ccL] ;
	    ccSR = dnaDecodeChar[(int)ccSR] ;
	    ccLR = dnaDecodeChar[(int)ccLR] ;
	    
	    vtxtPrintf(txt1, "%s%d:%c>%c"
		       , sep
		       , xShort 
		       , ccLR, ccS
		       ) ;

	    vtxtPrintf(txt2, "%s%d:%c>%c"
		       , sep
		       , xLong 
		       , ccL, ccSR
		       ) ;
	  }  
	  break ;
	case TROU: 
	  {
	    char *ss = "-", ccL, ccLR ;
	    
	    ccL = arr(dnaG, xLong - 1, unsigned char) ;
	    ccLR = isUp ? complementBase(ccL) : ccL ; 
	    
	    ccL = dnaDecodeChar[(int)ccL] ;
	    ccLR = dnaDecodeChar[(int)ccLR] ;

	    vtxtPrintf (txt1, "%s%d:%s%c"
			, sep
			, xShort
			, ss
			, ccLR
			) ;
	    vtxtPrintf (txt2,"%s%d:%s%c"
			, sep
			, xLong 
			, ss 
			, ccL 
			) ;
	    
	  }
	  break ;

	case TROU_DOUBLE:
	  {
	    char *ss = "--", cc1L, cc2L, cc1LR, cc2LR ;
	    
	    cc1L = arr(dnaG, xLong - 1, unsigned char) ;
	    cc2L = arr(dnaG, xLong - 1 + 1, unsigned char) ;
	    cc1LR = isUp ? complementBase(cc2L) : cc1L ; 
	    cc2LR = isUp ? complementBase(cc1L) : cc2L ; 
	    
	    cc1L = dnaDecodeChar[(int)cc1L] ;
	    cc2L = dnaDecodeChar[(int)cc2L] ;
	    cc1LR = dnaDecodeChar[(int)cc1LR] ;
	    cc2LR = dnaDecodeChar[(int)cc2LR] ;

	    vtxtPrintf (txt1, "%s%d:%s%c%c"
			, sep
			, xShort + (isUp ? 1 : 0)
			, ss
			, cc1LR, cc2LR
			) ;
	    vtxtPrintf (txt2, "%s%d:%s%c%c"
			, sep
			, xLong 
			, ss 
			, cc1L, cc2L
			) ;
	    
	  }
	  break ;
	case TROU_TRIPLE:
	  {
	    char *ss = "---", cc1L, cc2L, cc3L, cc1LR, cc2LR, cc3LR ;
	    
	    cc1L = arr(dnaG, xLong - 1, unsigned char) ;
	    cc2L = arr(dnaG, xLong - 1 + 1, unsigned char) ;
	    cc3L = arr(dnaG, xLong - 1 + 2, unsigned char) ;
	    cc1LR = isUp ? complementBase(cc3L) : cc1L ; 
	    cc2LR = isUp ? complementBase(cc2L) : cc2L ;
	    cc3LR = isUp ? complementBase(cc1L) : cc3L ; 
	    
	    cc1L = dnaDecodeChar[(int)cc1L] ;
	    cc2L = dnaDecodeChar[(int)cc2L] ;
	    cc3L = dnaDecodeChar[(int)cc3L] ;
	    cc1LR = dnaDecodeChar[(int)cc1LR] ;
	    cc2LR = dnaDecodeChar[(int)cc2LR] ;
	    cc3LR = dnaDecodeChar[(int)cc3LR] ;

	    vtxtPrintf (txt1, "%s%d:%s%c%c%c"
			, sep
			, xShort + (isUp ? 1 : 0)
			, ss
			, cc1LR, cc2LR, cc3LR
			) ;
	    vtxtPrintf (txt2, "%s%d:%s%c%c%c"
			, sep
			, xLong 
			, ss 
			, cc1L, cc2L, cc3L
			) ;
	  }
	  break ;
	case INSERTION: 
	  {
	    char *ss = "+", cc1S, cc1SR ;
	    
	    cc1S = arr(dna, xShort - 1, unsigned char) ;
	    cc1SR = isUp ? complementBase(cc1S) : cc1S ;
	    
	    cc1S = dnaDecodeChar[(int)cc1S] ;
	    cc1SR = dnaDecodeChar[(int)cc1SR] ;

	    vtxtPrintf (txt1, "%s%d:%s%c"
			, sep
			, xShort + (isUp ? 0 : 0)
			, ss
			, cc1S
			) ;
	    vtxtPrintf (txt2, "%s%d:%s%c"
			, sep
			, xLong
			, ss 
			, cc1SR
			) ;
	    
	  }
	  break ;

	case INSERTION_DOUBLE:
	  {
	    char *ss = "++", cc1S, cc2S, cc1SR, cc2SR ;

	    xShort = isUp ? xShort - 1 : xShort ; 
	    cc1S = arr(dna, xShort - 1, unsigned char) ;
	    cc2S = arr(dna, xShort - 1 + 1, unsigned char) ;
	    cc1SR = isUp ? complementBase(cc2S) : cc1S ; 
	    cc2SR = isUp ? complementBase(cc1S) : cc2S ;
	    
	    cc1S = dnaDecodeChar[(int)cc1S] ;
	    cc2S = dnaDecodeChar[(int)cc2S] ;
	    cc1SR = dnaDecodeChar[(int)cc1SR] ;
	    cc2SR = dnaDecodeChar[(int)cc2SR] ;

	    vtxtPrintf (txt1, "%s%d:%s%c%c"
			, sep
			, xShort + (isUp ? 0 : 0)
			, ss
			, isUp ? cc1S : cc1S
			, isUp ? cc2S : cc2S
			) ;
	    vtxtPrintf (txt2, "%s%d:%s%c%c"
			, sep
			, xLong
			, ss 
			, cc1SR, cc2SR
			) ;
	    
	  }
	  break ;
	  
	case INSERTION_TRIPLE:
	  {
	    char *ss = "+++", cc1S, cc2S, cc3S, cc1SR, cc2SR, cc3SR ;
	    xShort = isUp ? xShort - 2 : xShort ; 
	    cc1S = arr(dna, xShort - 1, unsigned char) ;
	    cc2S = arr(dna, xShort - 1 + 1, unsigned char) ;
	    cc3S = arr(dna, xShort - 1 + 2, unsigned char) ;
	    cc1SR = isUp ? complementBase(cc3S) : cc1S ; 
	    cc2SR = isUp ? complementBase(cc2S) : cc2S ;
	    cc3SR = isUp ? complementBase(cc1S) : cc3S ; 
	    
	    cc1S = dnaDecodeChar[(int)cc1S] ;
	    cc2S = dnaDecodeChar[(int)cc2S] ;
	    cc3S = dnaDecodeChar[(int)cc3S] ;
	    cc1SR = dnaDecodeChar[(int)cc1SR] ;
	    cc2SR = dnaDecodeChar[(int)cc2SR] ;
	    cc3SR = dnaDecodeChar[(int)cc3SR] ;

	    vtxtPrintf (txt1, "%s%d:%s%c%c%c"
			, sep
			, xShort + (isUp ? 0 : 0)
			, ss
			, isUp ? cc1S : cc1S
			, isUp ? cc2S : cc2S
			, isUp ? cc3S : cc3S
			) ;
	    vtxtPrintf (txt2, "%s%d:%s%c%c%c"
			, sep
			, xLong
			, ss 
			, cc1SR, cc2SR, cc3SR
			) ;
	    
	  }
	  break ;
      
	default:
	  vtxtPrintf(txt2, "%s%d:%c>%c"
		     , sep
		     , xLong 
		     , 'Z', 'Z'
		     ) ;
	  
	  vtxtPrintf(txt1, "%s%d:%c>%c"
		     , sep
		     , xShort 
		     , 'Z', 'Z'
		     ) ;
	}
    }

  #ifdef JUNK
  int nk1 = strlen (vtxtPtr (txt1)) ;
  int nk2 = strlen (vtxtPtr (txt2)) ;
  if (nk1 > 200)
    invokeDebugger () ;
  if (nk2 > 200)
    invokeDebugger () ;
  if (bb->lane == 169)
    {
      nnE++ ;
      fprintf (myErrFile, "nnE=%d\tread=%d::%s\n%s\n\n", nnE, read, vtxtPtr (txt1), vtxtPtr (txt2)) ;
      fflush (myErrFile) ;
    }
  #endif
  
  dictAdd (bb->errDict, vtxtPtr (txt1), &up->errShort) ;
  dictAdd (bb->errDict, vtxtPtr (txt2), &up->errLong) ;
  
  return ;
} /* alignFormatErrors */

/**************************************************************/

static void alignClipErrorLeft (ALIGN *vp, int errCost)
{
  A_ERR *ep ;
  int i, iMax = arrayMax (vp->errors) ;
  int bestMax = -1 ;
  int x1 = vp->x1 ;
  
  if (iMax)
    {
      for (i = 0, ep = arrp (vp->errors, i, A_ERR) ; i < iMax ; i++, ep++)
	{
	  int dx = ep->iShort - x1 ;
	  if (dx < (i+1) * errCost)
	    { bestMax = i ; x1 = ep->iShort ; }
	}
      if (bestMax > -1)
	{
	  ep = arrp (vp->errors, bestMax, A_ERR) ;
	  if (0 && ep->iShort > vp->x2 - 5)
	    return ;
	  vp->x1 = ep->iShort + 2 ; /* last exact base */
	  if (vp->a1 < vp->a2)
	    vp->a1 = ep->iLong + 2 ;
	  else
	    {
	      vp->a1 = ep->iLong - 0 ;
	    }
	  switch (ep->type)
	    {   
	    case INSERTION:
	      vp->x1++ ;
	      break ;
	    case INSERTION_DOUBLE:
	      vp->x1 += 2 ;
	      break ;
	    case INSERTION_TRIPLE:
	      vp->x1 += 3 ;
	      break ;
	    case TROU:
	      vp->a1 += ((vp->chrom & 0x1) ? -1: 1) ;
	      break ;
	    case TROU_DOUBLE:
	      vp->a1 += ((vp->chrom & 0x1) ? -2: 2) ;
	      break ;
	    case TROU_TRIPLE:
	      vp->a1 += ((vp->chrom & 0x1) ? -3 : 3) ;
	      break ;
	    default:
	      break ;
	    }
	  int k = bestMax + 1 ;
	  for (i = 0, ep = arrp (vp->errors, i, A_ERR) ; i < iMax - k ; i++, ep++)
	    *ep = *(ep + k) ;
	  arrayMax (vp->errors) -= k ;	    
	}
    }
  
  
  return ;
} /* alignClipErrorLeft */

/**************************************************************/

static void alignClipErrorRight (ALIGN *vp, int errCost)
{
  A_ERR *ep ;
  int i, iMax = arrayMax (vp->errors) ;
  int bestMax = iMax ;
  int x2 = vp->x2 ;
  
  if (iMax)
    {
      for (i = iMax - 1, ep = arrp (vp->errors, i, A_ERR) ; i >= 0 ; i--, ep--)
	{
	  int dx = x2 - ep->iShort ;
	  if (dx < (bestMax - i) * errCost)
	    { bestMax = i ; x2 = ep->iShort; }
	}
      if (bestMax < iMax)
	{
	  ep = arrp (vp->errors, bestMax, A_ERR) ;
	  vp->x2 = ep->iShort ; /* last exact base */
	  vp->a2 = ep->iLong ;
	  arrayMax (vp->errors) = bestMax ;
	}
    }

  return ;
} /* alignClipErrorRight */

/**************************************************************/
/* locate again the chains to eliminate the chain == -1 killed exons
 * ATTENTION there is a second version below  alignLocateChains1 to avoid killed chain zero 
 */
static int alignLocateChains (Array bestAp, Array aa, int myRead)
{
  int i1, i2, jj, iMax = arrayMax (aa) ;
  ALIGN *up, *vp ;
  
  bestAp = arrayReCreate (bestAp, keySetMax (bestAp), KEY) ;
  if (iMax)
    {
      for (i1 = i2 = jj = 0, up = vp = arrp (aa, 0, ALIGN) ; i1 < iMax ; i1++, up++)
	{
	  if (up->chain > 0 && up->a1 != up->a2)
	    {
	      if (vp < up) *vp = *up ;
	      if (vp->read == myRead && vp->chain > jj)
		{
		  jj = vp->chain ;
		array (bestAp, jj, int) = i2 ;
		}
	      vp++ ; i2++ ;
	    }
	}
      iMax = arrayMax (aa) = i2 ;
    }
      
  up = arrayp (aa, iMax, ALIGN) ;
  memset (up, 0, sizeof (ALIGN)) ; /* force a null record */
  arrayMax (aa) = iMax ;

  return iMax ;
} /* alignLocateChains */

/**************************************************************/
/* locate again the chains to eliminate the chain <= 0 killed exons
 * ATTENTION there is a second version above  alignLocateChains to only avoid killed chain == -1 exons
 */
static int alignLocateChains1 (Array bestAp, Array aa, int myRead)
{
  int i1, i2, jj, iMax = arrayMax (aa) ;
  ALIGN *up, *vp ;
  
  bestAp = arrayReCreate (bestAp, keySetMax (bestAp), KEY) ;
  if (iMax)
    {
      for (i1 = i2 = jj = 0, up = vp = arrp (aa, 0, ALIGN) ; i1 < iMax ; i1++, up++)
	{
	  if (up->chain > 0 && up->a1 != up->a2)
	    {
	      if (vp < up) *vp = *up ;
	      if (vp->read == myRead && vp->chain > jj)
		{
		  jj = vp->chain ;
		array (bestAp, jj, int) = i2 + 1 ;
		}
	      vp++ ; i2++ ;
	    }
	}
      iMax = arrayMax (aa) = i2 ;
    }
      
  up = arrayp (aa, iMax, ALIGN) ;
  memset (up, 0, sizeof (ALIGN)) ; /* force a null record */
  arrayMax (aa) = iMax ;

  return iMax ;
} /* alignLocateChains */

/**************************************************************/

static void alignAdjustIntrons (const PP *pp, BB *bb, Array bestAp, Array aa)
{
  ALIGN *up, *vp, *wp ;
  int ii, jj, chromA = 0 ;
  Array dnaG = 0 ;
  int chain = 0 ;
  
  for (int ic = 1 ; ic < arrayMax (bestAp) ; ic++)
    {
      /* adjust introns */
      int k = array (bestAp, ic, int) ;
      if (!k) continue ;
      up = vp = arrp (aa, k - 1, ALIGN) ; 
      int chain = up->chain ;
      if (up->chrom != chromA)
	{
	  chromA = up->chrom ;
	  dnaG = arr (pp->bbG.dnas, chromA >> 1, Array) ;
	}

      
      wp = (vp[1].chain == chain ? vp + 1 : 0) ;
      while (vp && wp)
	{
	  if (vp > up && vp->x1 > wp->x1)
	    {
	      memset (vp, 0, sizeof (ALIGN)) ; vp->chain = -1 ;
	      saIntronsOptimize (bb, up, wp, dnaG) ;
	    }
	  else
	    saIntronsOptimize (bb, vp, wp, dnaG) ;
	  if (wp->chain != -1)   /* happens if the exons were merged */
	    vp = wp ;
	  wp = (vp[0].chain && wp[1].chain == vp[0].chain ? wp + 1 : 0) ;
	}
    }
  
  for (ii = jj = 0, up = vp = arrp (aa, ii, ALIGN) ; ii < arrayMax (aa) ; up++, ii++)
    {
      if (up->chain > 0)
	{
	  if (vp < up) *vp = *up ;
	  jj++ ; vp++ ;
	}
    }
  arrayMax (aa) = jj ;
  
  for (int ii = 0 ; ii < arrayMax (aa) ; ii++)
    {
      up = arrayp (aa, ii, ALIGN) ; 
      if (up->chain != chain)
	{
	  chain = up->chain ;
	  vp = up ;
	  up->chainX1 = up->x1 ;
	  up->chainA1 = up->a1 ;
	}
      for (wp = vp ; wp <= up ; wp++)
	{
	  wp->chainX1 = vp->chainX1 ;
	  wp->chainX2 = up->x2 ;
	  wp->chainA2 = up->a2 ;
	  wp->chainAli = wp->chainX2 - wp->chainX1 + 1 ;
	}
    }
  
} /* alignAdjustIntrons */

/**************************************************************/

static void alignAdjustExons (const PP *pp, BB *bb, Array bestAp, Array aa, Array dna, int maxJump, int maxJump2)
{
  ALIGN *up, *vp ;
  int chromA = 0, chain = 0 ;
  AC_HANDLE h1 = ac_new_handle () ;
  Array dnaR = 0 ;
  Array dnaG = 0 ;
  

  for (int ic = 1 ; ic < arrayMax (bestAp) ; ic++)
    {
      /* adjust   exons */
      int k = array (bestAp, ic, int) ;
      if (!k) continue ;
      up = vp = arrp (aa, k - 1, ALIGN) ; 

      chain = up->chain ;
      if (up->chrom != chromA)
	{
	  chromA = up->chrom ;
	  dnaG = arr (pp->bbG.dnas, chromA >> 1, Array) ;
	}
      if (0 && up->a1 == 99016239)
	invokeDebugger () ;
      
      /* count all errors */
      int nErr = 0 ;
      int da = 0, nAli = 0 ;
      int lnShort = arrayMax (dna) ;
      for (vp = up ; vp->chain == chain || vp->chain == -1 ; vp++)
	{
	  da = vp->a2 - vp->a1 ;
	  vp->nErr = vp->errors ? arrayMax (vp->errors) : 0 ;
	  nErr += vp->nErr ;
	  nAli += da > 0 ? da + 1 : -da + 1 ;
	}

      if (nErr || nAli < lnShort)
	{
	  KEYSET ks = keySetHandleCreate (h1) ;
	  Array dnaI = 0, dnaShort = dna, errors = arrayHandleCreate (20, A_ERR, bb->h) ;
	  int ie, ia = 0, da ;
	  int ln = arrayMax (dnaG) ;
	  char *cp, *cq ;
	  A_ERR *ep ;
	  ALIGN zp ;
	  BOOL isDown ;
	  memset (&zp, 0, sizeof (A_ERR)) ;
		  
	  /* reconstruct the image of the transcript */

	  isDown = (up->chrom & 0x1) ? FALSE : TRUE ;
      	  nAli = nAli + 200 ;

	  dnaI = arrayHandleCreate (nAli + 1, char, h1) ;
	  array (dnaI, nAli, char) = 0 ; /* add a terminal zero */
	  arrayMax (dnaI) = nAli ;
	  if (! isDown)
	    {
	      int nvp = 0 ;
	      if (! dnaR)
		{
		  dnaR = dnaHandleCopy (dna, h1) ;
		  reverseComplement (dnaR) ;
		}
	      dnaShort = dnaR ;
	      for (vp = up ; vp->chain == chain ; vp++)
		{
		  int dummy = vp->a1 ; vp->a1 = vp->a2 ; vp->a2 = dummy ;
		  dummy = vp->x1 ; vp->x1 = lnShort - vp->x2 + 1 ; vp->x2 = lnShort - dummy + 1 ;
		  nvp++ ;
		}
	      for (int i = 0 ; i < nvp/2 ; i++)
		{
		  ALIGN wp = up[i] ;
		  up[i] = up[nvp - i - 1] ;
		  up[nvp -i - 1] = wp ;
		}
	    }
	  if (1)
	    {
	      int jj = up->a1 > 100 ? 100 : up->a1 - 1 ;
	      keySet (ks, ia++) = up->a1 - jj ;
	      cp = arrp (dnaI, 0, char) ;
	      if (jj)
		{
		  cp = arrayp (dnaI, jj, char) ;
		  cp = arrp (dnaI, 0, char) ;
		  cq = arrp (dnaG, up->a1 - 1 - jj, char) ;
		  memcpy (cp, cq, jj) ;
		  cp += jj ;
		}
	      memset (&zp, 0, sizeof(zp)) ;
	      zp.x1 = up->x1 ;
	      zp.a1 = jj + 1 ;
	      for (vp = up ; vp->chain == chain || vp->chain == -1 ; vp++)
		{
		  da = vp->a2 - vp->a1 + 1 ;
		  if (vp->chain == -1 || da < 1) continue ;
		  keySet (ks, ia++) = vp->a1 - jj - 1 ;
		  cp = arrayp (dnaI, jj + da, char) ;
		  cq = arrp (dnaG, vp->a1 - 1 , char) ;
		  cp = arrayp (dnaI, jj, char) ;
		  memcpy (cp, cq, da) ;
		  jj += da ;
		  cp += da ;
		  if (vp[1].chain == chain)
		    {
		      int du = vp[1].x1 - vp->x2 - 1 ;
		      int du0 = du ;
		      cp = arrayp (dnaI, jj + du + 6, char) ;
		      cp = arrp (dnaI, jj, char) ;
		      BOOL isDonor = TRUE ;
		      if (du > 0)
			{ /* in the overlap OR the bases form the donor and acceptor exon */
			  cq = arrp (dnaG, vp->a2 - 1, char) ;
			  char *cr = arrp (dnaG, vp[1].a1 - 1, char) ;
			  if (cr[-2] == A_ && cr[-1] == G_ && cq[du + 1] == G_ && cq[du+2] == T_)
			    { /* copy the extension of the donor exon */
			      memcpy (cp, cq+1, du) ;
			    }
			  else if (du >= 1 && cr[-2] == A_ && cr[-1] == G_ && cq[du + 1 -1] == G_ && cq[du+2-1] == T_)
			    { /* copy the extension of the donor exon with deletion */
			      du-- ;
			      memcpy (cp, cq+1, du) ;
			    }
			  else if (du >= 2 && cr[-2] == A_ && cr[-1] == G_ && cq[du + 1 -2] == G_ && cq[du+2-2] == T_)
			    { /* copy the extension of the donor exon with double deletion */
			      du-=2 ;
			      memcpy (cp, cq+1, du) ;
			    }
			  else if (du >=3 && cr[-2] == A_ && cr[-1] == G_ && cq[du + 1 -3] == G_ && cq[du+2-3] == T_)
			    { /* copy the extension of the donor exon with triple deletion */
			      du-=3 ;
			      memcpy (cp, cq+1, du) ;
			    }
			  else if (du >= 1 && cr[-2] == A_ && cr[-1] == G_ && cq[du + 2] == G_ && cq[du+3] == T_)
			    { /* copy the extension of the donor exon with insertion */
			      du++ ;
			      memcpy (cp, cq+1, du) ;
			    }
			  else if (du >= 1 && cr[-2] == A_ && cr[-1] == G_ && cq[du + 3] == G_ && cq[du+4] == T_)
			    { /* copy the extension of the donor exon with double insertion */
			      du+=2 ;
			      memcpy (cp, cq+1, du) ;
			    }
			  else if (du >=1 && cr[-2] == A_ && cr[-1] == G_ && cq[du + 4] == G_ && cq[du+5] == T_)
			    { /* copy the extension of the donor exon with triple insertion */
			      du+=3 ;
			      memcpy (cp, cq+1, du) ;
			    }
			  else if (cq[1] == G_ && cq[2] == T_ && cr[- du - 2] == A_ && cr[-du - 1] == G_)
			    { /* copy the extension of the acceptor exon */
			      isDonor = FALSE ;
			      /* do not copy it will be automatic since we shift vp[1] below */
			    }
			  else if (du >= 0 && cq[1] == G_ && cq[2] == T_ && cr[- du - 2] == A_ && cr[-du - 1] == G_)
			    { /* copy the extension of the acceptor exon with insertion */
			      isDonor = FALSE ;
			      /* do not copy it will be automatic since we shift vp[1] below */
			    }
			  else if (du >= 0 && cq[1] == G_ && cq[2] == T_ && cr[- du - 3] == A_ && cr[-du - 2] == G_)
			    { /* copy the extension of the acceptor exon with insertion */
			      isDonor = FALSE ;
			      du++ ;
			      /* do not copy it will be automatic since we shift vp[1] below */
			    }
			  else if (du >= 0 && cq[1] == G_ && cq[2] == T_ && cr[- du - 4] == A_ && cr[-du - 3] == G_)
			    { /* copy the extension of the acceptor exon with insertion */
			      isDonor = FALSE ;
			      du+=2 ;
			      /* do not copy it will be automatic since we shift vp[1] below */
			    }
			  else if (du >= 0 && cq[1] == G_ && cq[2] == T_ && cr[- du - 1] == A_ && cr[-du - 0] == G_)
			    { /* copy the extension of the acceptor exon with deletion */
			      isDonor = FALSE ;
			      du-- ;
			      /* do not copy it will be automatic since we shift vp[1] below */
			    }
			  else if (du >= 0 && cq[1] == G_ && cq[2] == T_ && cr[- du - 0] == A_ && cr[-du + 1] == G_)
			    { /* copy the extension of the acceptor exon with deletion */
			      isDonor = FALSE ;
			      du-=2 ;
			      /* do not copy it will be automatic since we shift vp[1] below */
			    }
			  else if (du >= 0 && cq[1] == G_ && cq[2] == T_ && cr[- du + 1] == A_ && cr[-du + 2] == G_)
			    { /* copy the extension of the acceptor exon with deletion */
			      isDonor = FALSE ;
			      du-=3 ;
			      /* do not copy it will be automatic since we shift vp[1] below */
			    }
			  else
			    memset (cp, N_, du) ;
			  if (1 || isDonor)
			    {
			      jj += du0 - du ;
			      cp += du0 - du ;
			      vp->a2 += du0 - du ; vp->x2 += du0 - du ;
			    }
			  else
			    du = du0 ;
			  vp[1].x1 -= du ; vp[1].a1 -= du ;
			}
		    }
		}
	      vp-- ; ia-- ;
	      zp.x2 = vp->x2 ;
	      zp.a2 = jj ;

	      if (vp->x2 < lnShort)
		{
		  cp = arrp (dnaI, jj, char) ;
		  da = ln > vp->a2 + 100 ? 100 : ln - vp->a2 ;
		  if (da)
		    {
		      cp = arrayp (dnaI, jj + da, char) ; /* may reallocate */
		      cp = arrp (dnaI, jj, char) ;
		      cq = arrp (dnaG, vp->a2, char) ;
		      memcpy (cp, cq, da ) ;
		      jj += da ;
		    }
		  arrayMax (dnaI) = jj ;
		}
	      
	      /* align the read on the genomic image of the transcript */
	      zp.errors = errors ;
	      arrayMax (zp.errors) = 0 ;
	      if (1)
	      {
		int x2 = zp.x2,  a2 = zp.a2 ;
		aceDnaDoubleTrackErrors (dnaShort, &(zp.x1), &(zp.x2), TRUE   /* isDown = TRUE */
					 , dnaI, 0, &(zp.a1), &(zp.a2)
					 , 0, zp.errors, maxJump, -3, TRUE, 0) ; /* bio coordinates, extend right */
		if (zp.x2 < x2 - 50)
		  {
		    zp.x2 = x2 ; zp.a2 = a2 ;
		    arrayMax (zp.errors) = 0 ;
		    aceDnaDoubleTrackErrors (dnaShort, &(zp.x1), &(zp.x2), TRUE   /* isDown = TRUE */ 
					     , dnaI, 0, &(zp.a1), &(zp.a2)
					     , 0, zp.errors, maxJump2, -1, FALSE, 0) ; /* bio coordinates, jump 8 but do not extend */
		  }
	      }
	      alignClipErrorLeft (&zp, pp->errCost) ;
	      alignClipErrorRight (&zp, pp->errCost) ;
	      /* remap */
	      int ja, dda = 0 ;
	      for (vp = up, ja = 1 ; vp->chain == chain || vp->chain == -1 ; vp++)
		{
		  int dz = keySet (ks, ja)  ;
		  int j = 0 ;
		  int del = 0, ins = 0, sub = 0 ;
		  da = vp->a2 - vp->a1 + 1 ;
		  if (vp->chain == -1 || da < 1) continue ;
		  if (vp->a1 > zp.a2 + dz)
		    { vp->chain = -1 ; continue ; }
		  if (vp->a1 < zp.a1 + dz) { vp->a1 = zp.a1 + dz ; vp->x1 = zp.x1 ; }
		  if (vp->a2 > zp.a2 + dz) { vp->a2 = zp.a2 + dz ; vp->x2 = zp.x2 ; }
		  
		  if (vp->errors)
		    arrayMax (vp->errors) = 0 ; 
		  for (ie = 0 ; ie < arrayMax (zp.errors) ; ie++)
		    {
		      ep = arrp (zp.errors, ie, A_ERR) ;
		      int es = ep->iShort + 1 ;
		      if (es >= vp->x1 && es <= vp->x2)
			{
			  if (! vp->errors)
			    vp->errors = arrayHandleCreate (8, A_ERR, bb->h) ;
			  A_ERR *eq = arrayp (vp->errors, j++, A_ERR) ;
			  *eq = *ep ;
			  eq->iLong += dz ; /* + dda ; */
			  switch (ep->type)
			    {
			    case INSERTION: dda++ ; ins++ ; break ;
			    case INSERTION_DOUBLE: dda+=2 ; ins += 2 ; break ;
			    case INSERTION_TRIPLE: dda+=3 ; ins += 3 ; break ;
			    case TROU: dda-- ; del++ ; break ;
			    case TROU_DOUBLE: dda-=2 ; del += 2 ; break ;
			    case TROU_TRIPLE: dda-=3 ; del += 3 ; break ;
			    case ERREUR: sub++ ; break ;
			    default : break ;
			    }
			}
		    }
		  int dd = (vp->x2 - vp->x1 + 1 + del) - ((vp->a1 < vp->a2 ? vp->a2 - vp->a1 : vp->a1 - vp->a2) + 1 + ins) ;

		  if (dd < 0) vp->x2 += dd ;
		  if (dd > 0) vp->a2 -= dd ;

		  ja++ ;
		  vp->nErr = vp->errors ? arrayMax (vp->errors) : 0 ;
		  vp->nMID = del + ins + sub ;
		}
	    }
	  /*
	  for (vp = up ; vp->chain == chain ; vp++)
	    {
	      if (vp->x1 < 1)
		{ int dx = 1 - vp->x1 ; vp->x1 += dx ; vp->a1 += (vp->a1 < vp->a2 ? dx : -dx) ;}
	      if (vp->x2 > arrayMax (dna))
		{ int dx = vp->x2 - arrayMax(dna) ; vp->x2 -= dx ; vp->a2 -= (vp->a1 < vp->a2 ? dx : -dx) ;}
	    }
	  */
	  if (! isDown)
	    {
	      int nvp = 0 ;
	      for (vp = up ; vp->chain == chain ; vp++)
		{
		  int dummy = vp->a1 ; vp->a1 = vp->a2 ; vp->a2 = dummy ;
		  dummy = vp->x1 ; vp->x1 = lnShort - vp->x2 + 1 ; vp->x2 = lnShort - dummy + 1 ;
		  nvp++ ;

		  /* flip the error positions */
		  if (vp->nErr)
		    {
		      int i ;
		      for (i = 0, ep = arrp (vp->errors, 0, A_ERR) ; i < vp->nErr ; i++, ep++)
			{
			  ep->iShort = lnShort - ep->iShort - 1 ;
			  if (0 && ep->baseShort > 17)
			    invokeDebugger () ;
			  ep->baseShort = complementBase(ep->baseShort & 0xf) ;
			  ep->sens = -1 ;
			}
		    }
		}
	      for (int i = 0 ; i < nvp/2 ; i++)
		{
		  ALIGN wp = up[i] ;
		  up[i] = up[nvp - i - 1] ;
		  up[nvp-i - 1] = wp ;
		}
	    }
	}
    }

  ac_free(h1) ;
  
} /* alignAdjustExons */

/**************************************************************/
/* Dynamic programming of path score */

static void alignSelectBestDynamicPath (const PP *pp, BB *bb, Array aaa, Array aa, Array dna, int chromA, Array dnaG, Array dnaGR, Array bestAp, int maxJump, int maxJump2)
{
  AC_HANDLE h = 0 ;
  int ii, jj, i1, i2, iMax ;
  ALIGN *up, *vp, *wp ;
  int i02 = 0 ;
  int chainAli = 0 ;
  int chainScore = 0 ;
  int bestChainScore = 0 ;
  int bestI1 = 0 ;
  int chain = 0 ;
  int maxIntron = pp->maxIntron ;
  int minAli = pp->minAli ;
  Array aaNew = 0 ;
  int errCost = pp->errCost ;
  int bigErrCost = 8 ; /* errCost ; */
  int myRead = arrp (aa, 0, ALIGN)->read ;
  iMax = arrayMax (aa) ;
  arraySort (aa, saAlignOrder) ;

  /* create scores */
  if (iMax)
    {
      for (ii = jj = 0, up = vp = arrp (aa, 0, ALIGN) ; ii < iMax ; ii++, up++)
	{
	  int ali = up->x2 - up->x1 + 1 ;
	  int score = ali - up->nErr * errCost ;
	  
	  if (score > 0)
	    {
	      up->chain = 0 ;
	      up->chainScore = 0 ;
	      up->ali = ali ;
	      up->id = jj + 1 ;
	      up->score = up->ali - up->nErr * bigErrCost ;
	      if (vp < up) *vp = *up ;
	      if (1)
		{
		  wp = vp - 1 ;
		  if (jj == 0 || vp->chrom != wp->chrom || vp->a1 != wp->a1 || vp->a2 != wp->a2 || vp->x1 != wp->x1 || vp->x2 != wp->x2)
		    {
		      vp++ ; jj++ ;
		    }
		  else if (jj)
		    {
		      if (vp->donor && ! wp->donor)
			wp->donor = vp->donor ;
		      if (vp->acceptor && ! wp->acceptor)
			wp->acceptor = vp->acceptor ;
		      }
		}
	    }
	}
    }

  iMax = arrayMax (aa) = jj ;
  if (! iMax) return ;
  
  /* for each exon moving left to right, construct the chains */
  i2 = 0 ; vp = arrp (aa, 0, ALIGN) ; /* preposition */
  
  if (0 && bb->lane == 169 && myRead == 41873)
    invokeDebugger() ;

  for (i1 = 0, up = arrp (aa, 0, ALIGN) ; i1 < iMax ; i1++, up++)
    {
      int chrom = up->chrom ;
      int x1 = up->x1 ;
      int x2 = up->x2 ;
      int a1 = up->a1 ;
      int a2 = up->a2 ;
      int bestPrevious = 0 ;
      int bestPreviousScore = 0 ;
      int bestPreviousAli = 0 ;
      BOOL isDown = ! (chrom & 0x1) ;
      BOOL foundI2 = FALSE ;
      
      i2 = i02 ; vp = arrp (aa, i2, ALIGN) ; /* preposition */
      for (foundI2 = FALSE ; pp->splice && i2 < iMax ; i2++, vp++)
	{
	  if (i1 == i2)
	    continue ;
	  if (vp->chrom > chrom || (vp->chrom == chrom && vp->x1 > x2 - maxJump2))
	    break ;
	  if (vp->chrom < chrom || vp->x2 < x1 - maxJump2)
	    {
	      if (!foundI2) i02 = i2 + 1 ;
	      continue ;
	    }  
	  foundI2 = TRUE ;
	  if (vp->chrom == chrom
	      && vp->x2 >= x1 - 8 && vp->x2 < x2 && vp->x1 < x2
	      &&
	      (
	       ( isDown && vp->a1 < a2 && vp->a2 + maxIntron > a1 && vp->a2 - vp->x2 < a1 - x1 + maxJump2) ||
	       ( ! isDown && vp->a1 > a2 && vp->a2 - maxIntron < a1 && vp->a2 + vp->x2 > a1 + x1  - maxJump2)
	       )
	      )
	    {
	      int dx = vp->x2 - x1 + 1 > 0 ? vp->x2 - x1 + 1 : 0 ;
	      if (vp->chainScore - dx > bestPreviousScore)
		{
		  bestPreviousScore = vp->chainScore - dx ;
		  bestPreviousAli = vp->chainAli - dx ;
		  bestPrevious = i2 + 1 ;
		}
	      }
	}
      
      up->chainScore = up->score + bestPreviousScore ;
      up->chainAli = up->ali + bestPreviousAli ;
      up->previous = bestPrevious ;
      if (up->chainScore > bestChainScore)
	{
	  bestChainScore = up->chainScore ;
	  bestI1 = i1 ;
	}
    }
  
  /* attribute a chain number recursivelly to the best paths */
  chain = 0 ;
  i1 = bestI1 ;
  h = ac_new_handle () ;
  aaNew = arrayHandleCreate (iMax, ALIGN, h) ;

  while (bestChainScore > 0)
    {
      /* find the top of the chain */
      up = arrp (aa, bestI1, ALIGN) ;
      chainScore = up->score ;
      chainAli = up->chainAli ;
      while (up->previous)
	{
	  vp = arrp (aa, up->previous - 1, ALIGN) ;
	  vp->next = up->id ;
	  int dx = vp->x2 - up->x1 + 1 ;
	  chainScore += vp->score + (dx> 0 ? -dx : 0) ;
	  chainAli += vp->ali + (dx> 0 ? -dx : 0) ;

	  i1 = up->previous ;
	  up = vp ;
	}
      /* register this chain in aaNew, kill it in aa */
      chain++ ;
      int j0 = 0, jj = arrayMax (aaNew) ;
      array (bestAp, chain, int) = jj ;  /* best chain */
      wp = 0 ;
      while (up)
	{
	  up->chainScore =  chainScore ;
	  up->chainAli =  chainAli ;
	  vp = arrayp (aaNew, jj++, ALIGN) ;
	  *vp = *up ;
	  
	  vp->id = jj ;
	  if (j0)
	    {
	      vp[-1].next = jj ;
	      vp->previous = j0 ;
	    }
	  else
	    {
	      wp = vp ;
	      wp->chainX1 = up->x1 ;
	    }
	  j0 = jj ;
	  up->chain = vp->chain = chain ;
	  wp->chainX2 = vp->x2 ;
	  up = (up->next ? arrayp (aa, up->next - 1, ALIGN) :0) ;
	}
      /* filter on chainX1/chainX2 */
      BOOL ok = TRUE ;
      if (wp->chainX2 - wp->chainX1 + 1 < minAli)
	ok = FALSE ;
      for (int ic = 1 ; ok && ic < chain ; ic++)
	{
	  int iw = arr (bestAp, ic, int) ;
	  up = arrp (aaNew, iw, ALIGN) ;
	  int z1 = (up->chainX1 > wp->chainX1 ? up->chainX1 : wp->chainX1) ;
	  int z2 = (up->chainX2 < wp->chainX2 ? up->chainX2 : wp->chainX2) ;
	  int dz = z2 - z1 ;
	  int du = up->chainX2 - up->chainX1 ;
	  int dw = wp->chainX2 - wp->chainX1 ;
	  int tc = *dictName(pp->bbG.dict,wp->chrom >> 1) ;
	  if (
	      (2 * dz > du || 2 * dz > dw) /* significant overlap */
	      && (up->chainScore > 2 * wp->chainScore   || up->chainScore > wp->chainAli + pp->bonus[tc]) /* 10: allow for class bonus */
	      )
	    ok = FALSE ;
	}
      if (! ok)
	{  /* destroy this chain in aaNew */
	  int iw = arr (bestAp, chain, int) ;
	  arrayMax (aaNew) = iw ;
	  arrayMax (bestAp) = chain ;
	  chain-- ;

	  /* flag the bad chain */
	  up = arrp (aa, bestI1, ALIGN) ;
	  up->chain = -1 ;
	  while (up->previous)
	    {
	      up = arrp (aa, up->previous - 1, ALIGN) ;
	      up->chain = -1 ;
	    }	  
	}
      
      /* edit the score of the other exons
       * and look for new best score
       * this applies even for destroyed chains 
       */
      bestI1 = 0 ;
      bestChainScore = 0 ;
      for (i1 = 0, up = arrp (aa, 0, ALIGN) ; i1 < iMax ; i1++, up++)
	{
	  if (up->chain) /* allready used */
	    continue ;
	  if (up->previous && (vp = arrp (aa, up->previous - 1, ALIGN)) && !vp->chain)
	    up->chainScore = up->score + vp->chainScore ;
	  else
	    { up->chainScore = up->score ; up->previous = 0 ; } /* disconnect */
	  if (up->chainScore > bestChainScore)
	    { bestChainScore = up->chainScore ; bestI1 = i1 ; }
	}
    }

  /* transfer the sorted chains back in aa */
  iMax = arrayMax (aaNew) ;
  if (iMax)
    for (i1 = i2 = 0, up = arrp (aaNew, 0, ALIGN) ; i1 < iMax ; i1++, up++)
      {
	if (up->chain > 0)
	  {
	    vp = arrp (aa, i2++, ALIGN) ;
	    *vp = *up ;
	  }
      }
  iMax =   arrayMax (aa) = i2 ;
  arrayMax (bestAp) = 0 ;
  if (iMax)
    for (i1 = i2 = jj = 0, up = vp = arrp (aa, 0, ALIGN) ; i1 < iMax ; i1++, up++)
      {
	if (1)
	  {
	    if (vp < up) *vp = *up ;
	    if (vp->chain > jj)
	      {
		jj = vp->chain ;
		array (bestAp, jj, int) = i2 ;
	      }
	    vp++ ; i2++ ;
	  }
      }
  iMax = arrayMax (aa) = i2 ;      /* readjust */
  
  /* adjust introns and scores
   * both operations need to know the genome 
   */
  
  /* locate the chains */
  iMax = alignLocateChains1 (bestAp, aa, myRead) ;

  if (iMax)
    {
      /* adjust introns */
      alignAdjustIntrons (pp, bb, bestAp, aa) ;
      iMax = alignLocateChains1 (bestAp, aa, myRead) ;
      /* adjust exons */
      alignAdjustExons (pp, bb, bestAp, aa, dna, maxJump, maxJump2) ;
      iMax = alignLocateChains (bestAp, aa, myRead) ;
      
      /* Compute the clean chain score */
      bestChainScore = 0 ;
      for (int ic = 1 ; ic < arrayMax (bestAp) ; ic++)
	{
	  ii =array (bestAp, ic, int) ;
	  up = arrp (aa, ii, ALIGN) ; 
	  
	  int tc = *dictName(pp->bbG.dict,up->chrom >> 1) ;
	  int chain = up->chain ;
	  int chainX1 = up->x1 ;
	  int chainX2 = up->chainX2 ;
	  
	  int chainAli = up->chainAli = 0 ;
	  int chainErr = up->chainErr = 0 ;
	  int chainMID = up->chainMID = 0 ;
	  int chainScore = pp->bonus[tc] ;
	  int chainA1 = up->a1 ;
	  int chainA2 = up->a2 ;
	  
	  for (vp = up, jj = ii ; jj < iMax && vp->chain == chain ; jj++, vp++)
	    {
	      vp->targetClass = tc ;
	      vp->ali = vp->x2 - vp->x1 + 1 ;
	      vp->score = vp->ali - errCost * vp->nErr ;
	      if (vp->score < -10)
		continue ;
	      chainX2 = vp->x2 ;
	      chainAli += vp->ali ;
	      chainErr += vp->nErr ;
	      chainMID += vp->nMID ;
	      chainScore += vp->score ;
	      
	      if (chainA1 > vp->a1) chainA1 = vp->a1 ;
	      if (chainA1 > vp->a2) chainA1 = vp->a2 ;
	      if (chainA2 < vp->a1) chainA2 = vp->a1 ;
	      if (chainA2 < vp->a2) chainA2 = vp->a2 ;
	      
	    }
	  if (chainAli > arrayMax (dna))
	    chainAli = arrayMax (dna) ;

	  /* filter */
	  if (chainScore < pp->minScore ||
	      chainAli < pp->minAli ||
	      100 * chainAli < pp->minAliPerCent * arrayMax (dna) ||
	      100 * chainErr > pp->errRateMax * chainAli
	      )
	    chainScore = chainAli = chainErr = chainMID = 0 ;
	  
	  if (bestChainScore < chainScore)
	    bestChainScore = chainScore ;
	  
	  /* set the chain values in all exons */
	  if (up->a1 > up->a2)
	    { int dummy = chainA1 ; chainA1 = chainA2 ; chainA2 = dummy ; }
	  for (vp = up, jj = ii ; jj < iMax && vp->chain == chain ; jj++, vp++)
	    {
	      vp->chainScore = chainScore ;
	      vp->chainAli = chainAli ;
	      vp->chainErr = chainErr ;
	      vp->chainMID = chainMID ;
	      vp->chainX1 = chainX1 ;
	      vp->chainX2 = chainX2 ;
	      vp->chainA1 = chainA1 ;
	      vp->chainA2 = chainA2 ;
	    }
	  up = vp - 1 ; ii = jj - 1 ;
	}
    }
  
  /* clean up the destroyed chains and adjust the chain numbers */
  iMax = arrayMax (aa) ;

  if (iMax)
    {
      arraySort (aa, saAlignOrder) ;
      
      int newChain = 0 ;
      for (i1 = i2 = chain = 0,  up = vp = arrp (aa, 0, ALIGN) ; i1 < iMax ; i1++, up++)
	{
	  if (up->chain && up->chainScore > 0 && up->score > -10)
	    {
	      int k = 0 ;
	      if (1)
		{
		  for (int i3 = 1 ; up->score && i3 <= newChain ; i3++)
		    {
		      wp = arrp (aa, arr (bestAp, i3, int), ALIGN) ;
		      if (wp->chainScore > up->chainScore)
			{
			  int z1 = (up->chainX1 > wp->chainX1 ? up->chainX1 : wp->chainX1) ;
			  int z2 = (up->chainX2 < wp->chainX2 ? up->chainX2 : wp->chainX2) ;
			  int dz = z2 - z1 ;
			  int du = up->chainX2 - up->chainX1 ;
			  int dw = wp->chainX2 - wp->chainX1 ;
			  
			  if  (2 * dz > du || 2 * dz > dw) /* significant overlap */
			    {
			      int c = up->chain ;
			      for (wp = up ; i1 < iMax && wp->chain == c ; wp++)
				{ k++ ; wp->chainScore = wp->chain = wp->score = 0 ; }
			      i1 += k - 1 ; up += k - 1 ;
			    }
			}
		      
		    }
		  if (up->score)
		    {
		      int c = up->chain ;
		      newChain++ ;
		      array (bestAp, newChain, int) = i2 ;
		      for (wp = up ; i1 < iMax && wp->chain == c ; wp++)
			if (wp->score)
			  {
			    wp->chain = newChain ;
			    if (vp < wp) *vp = *wp ;
			    
			    i2++ ; vp++ ;
			    i1++ ; up++ ;
			  }
		      i1-- ; up-- ;
		    }		
		}
	    }
	}
      arrayMax (bestAp) = newChain ? newChain + 1 : 0 ;
      arrayMax (aa) = i2 ;
    }
  
  /* register the alignments */
  int kMax = arrayMax (aaa) ;
  iMax = arrayMax (aa) ;
  if (iMax)
    {
      up = arrp (aa, 0, ALIGN) ;
      /*   bitSet (bb->isAligned, up->read) ; */
      for (ii = 0 ; ii < iMax ; ii++, up++)
	{
	  if (up->score > -10)
	    {
	      vp = arrayp (aaa, kMax++, ALIGN) ;
	      *vp = *up ;
	      vp->nChains = arrayMax (bestAp) ;
	    }
	}
    }
  iMax = alignLocateChains (bestAp, aaa, myRead) ;  
  ac_free (h) ;
  
  return ;
} /* alignSelectBestDynamicPath */

/**************************************************************/
/* Establish chain scores, select best */
static void  alignDoRegisterOnePair (const PP *pp, BB *bb, BigArray aaa, Array aa, int read, Array bestAp, ADAPTORS *adaptors)
  
{
  ALIGN *ap, *vp ;
  int ii ;
  int iMax = alignLocateChains1 (bestAp, aa, read) ;  
  int nChains = 0 ;
  char allTc[256] ;
  Array dna = 0, dna1 = 0, dna2 = 0 ;
  int chromA = 0 ;
  Array dnaG = 0, dnaGR = 0 ;
  int read1 = read & (~0x1) ;
  int read2 = read | 0x1 ;

  dna1 = arr (bb->dnas, read1, Array) ;
  dna2 = arr (bb->dnas, read2, Array) ;
  
  /* create chains */
  /* overhangs */
  iMax = arrayMax (aa) ;
  if (iMax)
    for (ii = 0, ap = arrp (aa, ii, ALIGN) ; ii < iMax ; ii++, ap++)
      if (read == ap->read)
	{
	  dna = ap->read & 0x1 ? dna2 : dna1 ;
	  ap->leftClip = ((ap->read & 0x1) ? bb->rc.jump5r2 : bb->rc.jump5r1) ;
	  if (ap->leftClip > arrayMax (dna))
	    ap->leftClip = 0 ;
	  ap->rightClip = arrayMax (dna) ;
	  if (ap->x1 > 1 && ap->x1 == ap->chainX1)
	    {
	      if (ap->chrom != chromA)
		{
		  chromA = ap->chrom ;
		  dnaG = arr (pp->bbG.dnas, chromA >> 1, Array) ;
		  dnaGR = arr (pp->bbG.dnasR, chromA >> 1, Array) ;
		}
	      ap->leftClip = alignFormatLeftOverhang (pp, bb, ap, dna, dnaG, dnaGR, adaptors, ii, iMax) ;
	      
	    }
	  if (ap->x2 == ap->chainX2 && ap->x2 < arrayMax (dna))	
	    {
	      if (ap->chrom != chromA)
		{
		  chromA = ap->chrom ;
		  dnaG = arr (pp->bbG.dnas, chromA >> 1, Array) ;
		  dnaGR = arr (pp->bbG.dnasR, chromA >> 1, Array) ;
		}
	      ap->rightClip = alignFormatRightOverhang (pp, bb, ap, dna, dnaG, dnaGR, adaptors, ii, iMax) ;
	    }
	}

  /* format the errors */
  if (0 && bb->lane == 169 && read == 41873)
    invokeDebugger() ;
  /*
    run -x Aligners/011_SortAlignG5R5/IDX.T2T.18.31 --maxTargetRepeats 31 -i titi.f.fasta+titi.r.fasta --align --method 011_SortAlignG5R5 --run FrontalCortex_CHP_Chimpanzee -o toto4 --step 5 --numactl --nB 1 --nA 1
  */
  
  iMax = arrayMax (aa) ;
  if (iMax)
    for (ii = 0, ap = arrp (aa, ii, ALIGN) ; ii < iMax ; ii++, ap++)
      if (read == ap->read)
	{
	  if (arrayExists (ap->errors))
	    {
	      unsigned int flip = 0 ;
	      
	      dna = ap->read & 0x1 ? dna2 : dna1 ;
	      if (ap->chrom != chromA)
		{
		  chromA = ap->chrom ;
		  dnaG = arr (pp->bbG.dnas, chromA >> 1, Array) ;
		  dnaGR = arr (pp->bbG.dnasR, chromA >> 1, Array) ;
		}
	      
	      if (bb->runStat.p.nPairs && (ap->read & 0x1))
		flip = 0x0f ; /* will flip last 4 bits */
	      if (arrayMax (ap->errors))
		mergeErrors (bb->errors, ap->errors, flip) ;
	      if (arrayMax (ap->errors))
		alignFormatErrors (pp, bb, ap, dna, dnaG, dnaGR, read) ;	  
	      if (! pp->sam && ! pp->bam)
		arrayDestroy (ap->errors) ;
	    }
	}
  
  /* global statistics */
  /* stranding : once per target class */
  memset (allTc, 0, sizeof (allTc)) ;
  
  for (int ic = 0 ; ic < arrayMax (bestAp) ; ic++)
    {
      int k = array (bestAp, ic, int) ;
      if (! k) continue ;
      ap = arrp (aa, k - 1, ALIGN) ;
      if (read == ap->read)
	{
	  int a1 = ap->a1 ;
	  int a2 = ap->a2 ;
	  int tc = ap->targetClass ;
	  BOOL s = a1 < a2 ;

	  if (ap->chrom & 0x1) s = !s ;
	  if (! tc || allTc[tc])
	    continue ;
	  allTc[tc] = 1 ;
	  if (ap->read & 0x1)
	    s = !s ;
	  if (ap->chrom & 0x1)
	    s = !s ;
	  if (s)
	    bb->runStat.GF[tc]++ ;
	  else
	    bb->runStat.GR[tc]++ ;
	  
	  bb->runStat.nReadsAlignedPerTargetClass[tc]++ ;
	  bb->runStat.nBasesAlignedPerTargetClass[tc]+=ap->chainX2 - ap->chainX1 + 1 ;
	}
    }
  /* increase the block stats */
  nChains = 0 ;
  if (arrayMax (bestAp))
    {
      int tc0 = 0 ;
      nChains = 0 ;
      memset (allTc, 0, sizeof (allTc)) ;
      
      ap = 0 ;
      for (int ic = 1 ; ic < arrayMax (bestAp) ; ic++)
	{
	  int k = array (bestAp, ic, int) ;
	  if (! k) continue ;
	  vp = arrp (aa, k - 1, ALIGN) ;
	  if (! vp->chainScore) continue ;
	  int tc = vp->targetClass ;
	  if (read == vp->read)
	    {
	      if (! ap)
		{
		  ap = vp ;
		  tc0 = ap->targetClass ;
		  bb->nAli++ ;
		  bb->runStat.nReadsAlignedPerTargetClass[0]++ ;
		  bb->runStat.nBasesAlignedPerTargetClass[0]+=ap->chainX2 - ap->chainX1 + 1 ;
		  bb->runStat.nMultiAligned[0]++ ;
		  nChains = 1 ;
		  bb->runStat.nAlignments++ ;
		  if (ap->chainErr == 0 && ap->leftClip + ap->chainAli >= ap->rightClip)
		    {
		      bb->runStat.nPerfectReads++ ;
		      if (0)
			fprintf (stderr, "PERFECT\t%s%s\n", dictName(bb->dict, (ap->read)>>1), (ap->read & 0x1 ? "<" : ">")) ;
		    }
		}
      
	      if (ap != vp && tc == tc0)
		{ /* count multiali only in main class */
		  int z1 = (ap->chainX1 > vp->chainX1 ? ap->chainX1 : vp->chainX1) ;
		  int z2 = (ap->chainX2 < vp->chainX2 ? ap->chainX2 : vp->chainX2) ;
		  int dz = z2 - z1 ;
		  int du = ap->chainX2 - ap->chainX1 ;
		  int dv = vp->chainX2 - vp->chainX1 ;
		  
		  if (2 * dz > du || 2 * dz > dv) /* significant overlap */
		    {
		      bb->runStat.nAlignments++ ;
		      nChains++ ;
		      continue ;
		    }
		}

	      if (tc == tc0)
		{
		  bb->runStat.nErr += vp->chainErr ;
		  bb->runStat.nMID += vp->chainMID ;
		  bb->aliDx += vp->chainAli ;
		  bb->aliDa += vp->chainAli ;
		  if (vp->read & 0x1)
		    bb->runStat.nBaseAligned2 += vp->chainAli ;
		  else
		    bb->runStat.nBaseAligned1 += vp->chainAli ;
		}
	    }
	}
      if (nChains)
	bb->runStat.nMultiAligned[nChains > 10 ? 10 : nChains]++ ;
    }

  /* measure the insert lengths */
#ifdef KUNK
  if (0 && read == read2 && kMax)
    {  /* Too complex to calculate because of the introns */
      int da = 0 ;
      ALIGN *wp ;
      wp = 0 ;
      ap = arrp (aa, 0, ALIGN) ;
      for (ii = kMax - 1, vp = bigArrayp (aaa, ii, ALIGN) ; ii >=0 && ap->read == read2 ; ii--, vp--)
	;
      vp++ ; ii++ ;
      if (vp->read == read2) wp = vp ;
      for (ii = ii - 1, vp = vp - 1 ; ii >=0 && ap->read == read2 ; ii--, vp--)
	;
      vp++ ; ii++ ;
      if (vp->read == read1) 
	{ /* vp, wp are psotionned on the best pair and we want to count only once */
	  if (vp->targetClass == wp->targetClass && vp->chrom == wp->chrom)
	    {
	      if (vp->chainA1 < vp->chainA2 &&
		  vp->chainA1 < wp->chainA1 &&
		  wp->chainA1 > wp->chainA2 &&
		  vp->chainA2 < wp->chainA1 &&
		  wp->chainA2 > vp->chainA1
		  )
		da = wp->chainA2 - vp->chainA2 + 1 ;
	      else if (vp->chainA1 > vp->chainA2 &&
		       vp->chainA1 > wp->chainA1 &&
		       wp->chainA1 < wp->chainA2
		       vp->chainA2 > wp->chainA1 &&
		       wp->chainA2 < vp->chainA1
		       )
		da = vp->chainA2 - wp->chainA2 + 1 ;
	    }
	  /* if (da < 0) the 2 aligned reads are overlapping */
	  if (da < 100000 &&
	      da > 0)
	    {  /* WRONG this includes the size of the jumped introns */
	      da /= 10 ;
	      if (da > 10000) da = 10000 ;
	      array (bb->insertLengthDistribution, da, long int)++ ;
	    }
	}
    }
#endif
  
  /* register the introns */
  int intronMaxOld = arrayMax (bb->confirmedIntrons) ;
  iMax = arrayMax (aa) ;
  if (iMax)
    {
      ALIGN *bp ;
      int jj ;
      int intronBonus = 0 ;
      
      ap = arrp (aa, 0, ALIGN) ;
      /*   bitSet (bb->isAligned, ap->read) ; */
      for (ii = 0 ; ii < iMax ; ii++, ap++)
	if (read == ap->read)
	  {
	    for (jj = ii + 1, bp = ap + 1 ; jj < iMax && ap->chrom == bp->chrom && bp->read == read && bp->chain == ap->chain && bp->x1 == ap->x2 + 1 ; jj++, bp++)
	      { /* found one intron */
		INTRON *zp = arrayp (bb->confirmedIntrons, arrayMax (bb->confirmedIntrons), INTRON) ;
		int chrom = ap->chrom ;
		int a1 = ap->a1 ;
		int a2 = ap->a2 ;
		int b1 = bp->a1 ;
		    
		BOOL isReadDown = a1 < a2 ? TRUE : FALSE ;

		zp->run = bb->run ;

		if (chrom != chromA)
		  {
		    chromA = chrom ;
		    dnaG = arr (pp->bbG.dnas, chromA >> 1, Array) ;
		    dnaGR = arr (pp->bbG.dnasR, chromA >> 1, Array) ;
		  }

		if (read & 0x1)
		  {
		    if (isReadDown)
		      { 
			zp->a1 = a2 + 1 ;
			zp->a2 = b1 - 1 ;
			bb->nIntronSupportPlus++ ;
			const char *cp = arrp (dnaG, zp->a1 - 1, char) ;
			zp->feet[4] = dnaDecodeChar[(int)complementBase(cp[0])] ;
			zp->feet[3] = dnaDecodeChar[(int)complementBase(cp[1])] ;
			zp->feet[2] = '_' ;
			cp = arrp (dnaG, zp->a2 - 2, char) ;
			zp->feet[1] = dnaDecodeChar[(int)complementBase(cp[0])] ;
			zp->feet[0] = dnaDecodeChar[(int)complementBase(cp[1])] ;
			zp->feet[5] = 0 ;
		      }
		    else
		      {
			zp->a1 = a2 - 1 ;
			zp->a2 = b1 + 1 ;
			bb->nIntronSupportMinus++ ;
			const char *cp = arrp (dnaG, zp->a1 - 1, char) ;
			zp->feet[4] = dnaDecodeChar[(int)cp[0]] ;
			zp->feet[3] = dnaDecodeChar[(int)cp[-1]]  ;
			zp->feet[2] = '_' ;
			cp = arrp (dnaG, zp->a2 - 1, char) ;
			zp->feet[1] = dnaDecodeChar[(int)cp[1]] ;
			zp->feet[0] = dnaDecodeChar[(int)cp[0]] ;
			zp->feet[5] = 0 ;
		      }
		    int a0 = zp->a1 ; zp->a1 = zp->a2 ; zp->a2 = a0 ;
		    zp->n = ap->chain ;
		    zp->chrom = chrom ^ 0x1;
		  }
		else
		  {
		    if (isReadDown)
		      { 
			zp->a1 = a2 + 1 ;
			zp->a2 = b1 - 1 ;
			bb->nIntronSupportPlus++ ;
			const char *cp = arrp (dnaG, zp->a1 - 1, char) ;
			zp->feet[0] = dnaDecodeChar[(int)cp[0]] ;
			zp->feet[1] = dnaDecodeChar[(int)cp[1]] ;
			zp->feet[2] = '_' ;
			cp = arrp (dnaG, zp->a2 - 2, char) ;
			zp->feet[3] = dnaDecodeChar[(int)cp[0]] ;
			zp->feet[4] = dnaDecodeChar[(int)cp[1]] ;
			zp->feet[5] = 0 ;
		      }
		    else
		      {
			zp->a1 = a2 - 1 ;
			zp->a2 = b1 + 1 ;
			bb->nIntronSupportMinus++ ;
			const char *cp = arrp (dnaG, zp->a1 - 1, char) ;
			zp->feet[0] = dnaDecodeChar[(int)complementBase(cp[0])] ;
			zp->feet[1] = dnaDecodeChar[(int)complementBase(cp[-1])] ;
			zp->feet[2] = '_' ;
			cp = arrp (dnaG, zp->a2 - 1, char) ;
			zp->feet[3] = dnaDecodeChar[(int)complementBase(cp[1])] ;
			zp->feet[4] = dnaDecodeChar[(int)complementBase(cp[0])] ;
			zp->feet[5] = 0 ;
		      }
		    zp->n = ap->chain ;
		    zp->chrom = chrom ;
		    
		    if (! strcmp (zp->feet, (chrom & 0x1 ? "gt_ag" : "gt_ag")))
		      {
			bb->runStat.gt_ag_Support++ ;
			if (bb->isRna >= 0 && bb->runStat.gt_ag_Support > bb->runStat.ct_ac_Support)
			  intronBonus++ ;
		      }
		    if (! strcmp (zp->feet, (chrom & 0x1 ? "ct_ac" : "ct_ac")))
		      {
			bb->runStat.ct_ac_Support++ ;
			if (bb->isRna >= 0 && bb->runStat.gt_ag_Support < bb->runStat.ct_ac_Support)
			  intronBonus++ ;
		      }
		    if (bb->isRna < 0)
		      intronBonus++ ; /* all discontinuities lower the global score */
		  }
		break ;
	      }
	  }

      intronBonus *= 4 * (bb->isRna >= 0 ? 1 : -1) ;
      
      ap = arrp (aa, 0, ALIGN) ;
      if (intronBonus)
	for (ii = 0 ; ii < iMax ; ii++, ap++)
	  if (read == ap->read)
	    {
	      ap->score += intronBonus ;
	      ap->chainScore += intronBonus ;
	    }
    }
  /* register the alignments */
  long int kMax = bigArrayMax (aaa) ;
  iMax = arrayMax (aa) ;
  if (iMax)
    {
      ap = arrp (aa, 0, ALIGN) ;
      /*   bitSet (bb->isAligned, ap->read) ; */
      for (ii = 0 ; ii < iMax ; ii++, ap++)
	if (read == ap->read)
	  {
	    vp = bigArrayp (aaa, kMax++, ALIGN) ;
	    *vp = *ap ;
	    vp->nChains = nChains ;
	    vp->nTargetRepeats = nChains ;
	  }
    }

     /* register the double introns */
  int intronMax = arrayMax (bb->confirmedIntrons) ;
  if (pp->introns)
    {
      
      for (ii = intronMaxOld ; ii < intronMax -1 ; ii++)
	{
	  INTRON *zp1 = arrayp (bb->confirmedIntrons, ii, INTRON) ;
	  INTRON *zp2 = arrayp (bb->confirmedIntrons, ii + 1, INTRON) ;
	  
	  if (zp1->chrom == zp2->chrom && zp1->n == zp2->n &&
	      (
	       (zp1->a1 < zp1->a2 && zp1->a2 < zp2->a1 && zp2->a1 < zp2->a2) ||
	       (zp1->a1 > zp1->a2 && zp1->a2 > zp2->a1 && zp2->a1 > zp2->a2) 
	       ) &&
	      ! strcmp (zp1->feet, "gt_ag") && ! strcmp (zp2->feet, "gt_ag") 
	      ) /* same chain chain */
	    {
	      DOUBLEINTRON *zzp = arrayp (bb->doubleIntrons, arrayMax (bb->doubleIntrons), DOUBLEINTRON) ;
	      zzp->chrom = zp1->chrom ;
	      zzp->run = bb->run ;
	      zzp->n = 1 ;
	      zzp->a1 = zp1->a1 ; zzp->a2 = zp1->a2 ;
	      zzp->b1 = zp2->a1 ; zzp->b2 = zp2->a2 ;
	      memcpy (zzp->feet1, zp1->feet, 6) ;
	      memcpy (zzp->feet2, zp2->feet, 6) ;
	    }
	}
      
  /* reset the count which was overlaoded with chain */
      for (ii = intronMaxOld ; ii < intronMax ; ii++)
	{
	  INTRON *zp1 = arrayp (bb->confirmedIntrons, ii, INTRON) ;
	  zp1->n = 1 ;
	}
    }
  return ;
  
} /* alignDoRegisterOnePair */

/**************************************************************/
static void alignDoOneRead (const PP *pp, BB *bb
			    , Array aaa, BigArray hits
			    , Array aa, Array err, Array bestAp
			    , int maxJump, int maxJump2)
{   
  BOOL debug = FALSE ;
  AC_HANDLE h = ac_new_handle () ;
  HIT * restrict hit ;
  ALIGN *ap ;
  long int ii, iMax = bigArrayMax (hits), kMax = 0 ;
  int a1, a2, x1, x2 ;
  int b1, b2, y1, y2, ha1, readOld = 0, chromOld = 0, readA = 0, chromA = 0, read1 = 0, iiGood = 0 ;
  BOOL isDownOld = TRUE ;
  Array dna = 0, dnaG = 0, dnaGR = 0 ;
  int errMax = pp->errMax ; /* 999999 ; */
  int chromLength = 0 ;
  /*
    HIT * restrict h1 ;
    int r1 = 0, nh1 = 0, chrom1 = 0 ;
  */
  int errCost = pp->errCost ;
  /*   unsigned int uu = 0 ; */
  int donor = 0, acceptor = 0 ;
  int intronBonus  = 1 ; /* this is a coodinate bonus, not a score bonus */
  
  int nTargetRepeats  = 1 ;
  int nTargetRepeatsOld = 0 ;
  const int nTRmask = (0x1 << NTARGETREPEATBITS) - 1 ;

  for (ii = 0, hit = bigArrp (hits, 0, HIT) ; ii < iMax ; ii++, hit++)
    {
      int read = hit->read ;
      int chrom = hit->chrom ;
      BOOL isDown = TRUE ;
      BOOL isIntron = ((hit->x1  >> NTARGETREPEATBITS )  & 0x7) ? TRUE : FALSE ;

      if (! read || ! chrom)
	continue ;
      if (ii < iMax  && ! memcmp (hit, hit + 1, sizeof (HIT)))
	continue ;
      if (read != read1)
	{
	  read1 = read ;
	  if (arrayMax (aa))
	    { /* create chain scores */
	      alignSelectBestDynamicPath (pp, bb, aaa, aa, dna, chromA, dnaG, dnaGR, bestAp, maxJump, maxJump2) ;
	    }
	  arrayMax (aa) = kMax = 0 ;
	}

      if (read != readA)
	{ readA = read ; dna = arr (bb->dnas, read, Array) ; }
      if (chrom != chromA)
	{
	  chromA = chrom ;
	  dnaG = arr (pp->bbG.dnas, chrom >> 1, Array) ;
	  dnaGR = arr (pp->bbG.dnasR, chrom >> 1, Array) ;
	  chromLength = arrayMax (dnaG) ;
	}

      x1 = hit->x1 ;
      nTargetRepeats = (x1 & nTRmask) ;
      if (arrayMax (dna) > 200 &&  nTargetRepeats > 10)
	continue ;
      x1 = x1 >> NTARGETREPEATBITS ;
      BOOL isIntronDown = (x1 >> 2) & 0x1 ;
      isDown = (chrom & 0x1)  ? FALSE : TRUE ;
      donor = x1 & 0x1 ;
      acceptor = x1 & 0x2 ;
      x1 = x1 >> 3 ; x2 = x1 + 1 ;    /* bio coordinates */   
      if (isDown)   /* plus strand of the genome */
	{
	  a1 =
	    hit->a1
	    + x1
	    + (isIntron ? intronBonus : 0) 
	    - 1 ;        /* compensate avoid zero */
	  a2 = a1 + 1 ;
	  if (donor) donor = a1 + 2 ; /* first base of intron */
	  if (acceptor) acceptor = a1 - 1 ; /* last base of intron */
	}
      else   /* minus strand of the genome */
	{
	  a1 =
	    hit->a1
	    - x1
	    + (isIntron ? intronBonus : 0) 
	    - 1 ;        /* compensate avoid zero */
	  a2 = a1 - 1 ;
	  if (donor) donor = a1 - 2 ; /* first base of intron */
	  if (acceptor) acceptor = a1 + 1 ; /* last base of intron */
	}

      if (! isIntronDown)
	{ donor = - donor ; acceptor = - acceptor ; }
      if (1 && read == readOld && chrom == chromOld && isDown == isDownOld &&
	  x1 >= y1 && x2 <= y2 && hit->a1 < ha1 + 3 &&   
	  (
	   (isDown && a1 >= b1 && a2 <= b2) ||
	   (! isDown && a1 <= b1 && a2 >= b2)
	   )
	  &&
	  nTargetRepeats >= nTargetRepeatsOld >> 1
	  )
	{
	  if (debug) fprintf (stderr, "Hit %ld\tr=%d\t%d\t%d\tc=%d\t%d\t%d\tDoublet of %d\t%s\t%d\n", ii, read, x1, x2, chrom, a1, a2, iiGood, dictName (pp->bbG.dict, chrom >> 1), hit->a1) ;
	  hit->read = 0 ;  /* remove doublets */
	  if (kMax)
	    {
	      if (donor)
		ap->donor = donor ;
	      if (acceptor)
		ap->acceptor = acceptor ;
	    }
	}
      else 
	{
	  chromOld = 0 ;
	  int a0 = a1, x0 = x1 ;
	  if (debug) fprintf (stderr, "Hit %ld\tr=%d\t%d\t%d\tc=%d\t%d\t%d\tbefore align\t%s\t%u\n"
			      , ii, read, x1, x2, chrom, a1, a2, dictName (pp->bbG.dict, chrom >> 1), hit->a1) ;
	  if (alignExtendHit (dna, dnaG, dnaGR, err, isDown, chromLength, &a1, &a2, &x1, &x2, errCost, isIntron, errMax, 22, maxJump))
	    {
	      if (debug) fprintf (stderr, "Hit %ld\tr=%d\t%d\t%d\tc=%d\t%d\t%d\tAccepted\t%s, u=%u, nErr=%d\n"
				  , ii, read, x1, x2, chrom, a1, a2
				  , dictName (pp->bbG.dict, chrom >> 1)
				  , hit->a1
				  , arrayMax (err)
				  ) ;
	      ap = arrayp (aa, kMax++, ALIGN) ;
	      memset (ap, 0, sizeof (ALIGN)) ;
	      ap->read = read ;
	      ap->chrom = chrom ;
	      ap->a0 = a0 ;
	      ap->a1 = a1 ;
	      ap->a2 = a2 ;
	      ap->x0 = x0 ;
	      ap->x1 = x1 ;
	      ap->x2 = x2 ;
	      ap->nTargetRepeats = nTargetRepeats ;
	      ap->donor = donor ;
	      ap->acceptor = acceptor ;
	      ap->readLength = arrayMax (dna) ;
	      ap->nErr = arrayMax (err) ;
	      if (ap->nErr)
		ap->errors = arrayHandleCopy (err, bb->h) ;
	      readOld = read ;
	      chromOld = chrom ;
	      isDownOld = isDown ;
	      nTargetRepeatsOld = nTargetRepeats ;
	      b1 = a1 ; b2 = a2 ; y1 = x1 ; y2 = x2 ; ha1 = hit->a1 ; iiGood = ii ;
	    }
	  else
	    {
	      if (debug) fprintf (stderr, "Hit %ld\tr=%d\t%d\t%d\tc=%d\t%d\t%d\tRejected\tu=%d\n", ii, read, x1, x2, chrom, a1, a2, hit->a1) ;
	      hit->read = 0 ; /* remove false positive */
	    }
	}
    }
  if (arrayMax (aa))
    { /* create chain scores */
      alignSelectBestDynamicPath (pp, bb, aaa, aa, dna, chromA, dnaG, dnaGR, bestAp, maxJump, maxJump2) ;
    }
  ac_free (h) ;
  return ;
} /* alignDoOneRead */

/**************************************************************/
/**************************************************************/
typedef struct pairStruct { ALIGN *up, *vp ; int score, chrom, a1, a2 ; } PAIR ;
static int pairOrder (const void *va, const void *vb)
{
  const PAIR *up = va ;
  const PAIR *vp = vb ;
  int n ;
  n = up->score - vp->score ; if (n) return -n ;
  n = up->chrom - vp->chrom ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;

  return 0 ;
} /* pairOrder */

/**************************************************************/
static void alignDoOnePair (const PP *pp, BB *bb
			    , BigArray aaaa, BigArray hits
			    , Array aa, Array err, ADAPTORS *adaptors)
{
  AC_HANDLE h = ac_new_handle () ;
  HIT *hit ;
  int ii, iMax = arrayMax (hits) ;
  int read1, read2 ;
  Array aaa = arrayHandleCreate (128, ALIGN, h) ;
  Array bestAp1 = arrayHandleCreate (8, int, h) ;
  Array bestAp2 = 0 ;
  int maxJump = MAXJUMP ;
  int maxJump2 = MAXJUMP2 ;

  /* maxJump = maxJump2 = 1 ; */ 
  read1 = bigArr (hits, 0, HIT).read ;
  read2 = bigArr (hits, iMax -1, HIT).read ;

  if (read1 != read2)
    {
      int iMax1 ;
      for (ii = 0, hit = bigArrp (hits, 0, HIT) ; ii < iMax ; hit++, ii++)
	if (hit->read == read2)
	  break ;
      arrayMax (hits)  = iMax1 = ii ;
      alignDoOneRead (pp, bb, aaa, hits, aa, err, bestAp1, maxJump, maxJump2) ;
      arrayMax (hits) = iMax ;
      if (iMax > iMax1)
	for (ii = iMax1, hit = bigArrp (hits, ii, HIT) ; ii < iMax ; hit++, ii++)
	  *(hit - iMax1) = *hit ;
      arrayMax (hits) = iMax - iMax1 ;
      arrayMax (aa) = arrayMax (err) = 0 ;
      bestAp2 = arrayHandleCreate (8, int, h) ;
      iMax1 = arrayMax (aaa) ;
      alignDoOneRead (pp, bb, aaa, hits, aa, err, bestAp2, maxJump, maxJump2) ;
    }
  else
    alignDoOneRead (pp, bb, aaa, hits, aa, err, bestAp1, maxJump, maxJump2) ;

  iMax = arrayMax (aaa) ;
  array (aaa, iMax, ALIGN).read = 0 ; /* impose a zero terminal record */
  arrayMax (aaa) = iMax ;
  
  if (bestAp2) /* we have a pair, good  example polyA_B_1 read 144/145*/
    {
      ALIGN *up, *vp ;
      int iMax1 = arrayMax (bestAp1) ;
      int iMax2 = arrayMax (bestAp2) ;

      BOOL hasPair = FALSE, hasGoodPair = FALSE, hasCirclePair = FALSE ;

      if (iMax1 && iMax2)
	{
	  Array pairs = arrayHandleCreate (iMax1 * iMax2 , PAIR, h) ;
	  PAIR *px ;
	  int i1, i2 ;
	  int jj = 0 ;

	  hasPair = TRUE ;
	  for (i1 = 1 ; i1 < iMax1 ; i1++)
	    for (i2 = 1 ; i2 < iMax2 ; i2++)
	      {
		int j1 = arr (bestAp1, i1, int) ;
		int j2 = arr (bestAp2, i2, int) ;

		up = arrp (aaa, j1, ALIGN) ;
		vp = arrp (aaa, j2, ALIGN) ;

		if ((up->chrom ^ vp->chrom) == 0x1)
		  {
		    if ((up->chrom & 0x1) == 1) { ALIGN *zp = up ; up = vp ; vp = zp ;}
		    int da = vp->chainA1 - up->chainA1 ;
		    int db = vp->chainA2 - up->chainA2 ;

		    if (da > 0 && db < 1000000) /* true pair */
		      {
			hasGoodPair = TRUE ;
			px = arrayp (pairs, jj++, PAIR) ;
			px->up = up ;
			px->vp = vp ;
			px->chrom = up->chrom ;
			px->a1 = up->a1 ;
			px->a2 = vp->a1 ;
			px->score = up->chainScore + vp->chainScore ;
		      }
		  }
	      }

	  int jMax = jj ;
	  if (jMax)
	    {
	      ALIGN *wp, *zp ;
	      int m ;
	      arraySort (pairs, pairOrder) ;
	      PAIR *px0 = arrayp (pairs, 0, PAIR), *qx = 0, *rx ;
	      int bestScore = px0->score, kk = 0 ;
	      Array aaa1 = arrayHandleCreate (arrayMax (aaa), ALIGN, h) ;
	      
	      for (jj = 0, px = px0 ; jj < jMax && px->score == bestScore ; jj++, px++)
		{
		  if (qx && px->up == qx->up && px->a2 > qx->a2)
		    continue ; /* eliminate vp2 in   ---> <--- <--- config */
		  BOOL ok = TRUE ;
		  for (m = jj + 1, rx = px+1 ; /* eliminate up1 in   ---> ---> <--- config */
		       ok && m < jMax && rx->a1 < px->a2 && rx->chrom == px->chrom && rx->score == bestScore ;
		       m++, rx++ )
		    ok = FALSE ;
		  if (!ok)
		    continue ;
		  qx = px ;

		  hasPair = TRUE ;
		  up = px->up ;
		  for (vp = up ; vp->chain == up->chain && vp->read == up->read ; vp++)
		    {
		      wp = px->vp ;
		      vp->pairScore = px->score ;
		      vp->mateChrom = wp->chrom ;
		      vp->mateA1 = wp->chainA1 ;
		      vp->mateA2 = wp->chainA2 ;

		      zp = arrayp (aaa1, kk++, ALIGN) ;
		      *zp = *vp ;
		    }	  


		  up = px->vp ;
		  for (vp = up ; vp->chain == up->chain && vp->read == up->read ; vp++)
		    {
		      wp = px->up ;
		      vp->pairScore = px->score ;
		      vp->mateChrom = wp->chrom ;
		      vp->mateA1 = wp->chainA1 ;
		      vp->mateA2 = wp->chainA2 ;

		      zp = arrayp (aaa1, kk++, ALIGN) ;
		      *zp = *vp ;
		    }	  
		}

	      arrayMax (aaa) = 0 ;
	      if (kk)
		{
		  zp = arrayp (aaa, kk, ALIGN) ;
		  memcpy (arrp (aaa, 0, ALIGN), arrp (aaa1, 0, ALIGN), kk * sizeof (ALIGN)) ;
		}
	    }
	}
      if (hasPair)
	{
	  bb->runStat.nPairsAligned++ ;
	  if (hasGoodPair)
	    bb->runStat.nCompatiblePairs++ ;
	  else if (hasCirclePair)
	    bb->runStat.nCirclePairs++ ;
	  else
	    bb->runStat.nIncompatiblePairs++  ;
	}
    }


  /* alignDoSelectBestPair */
  iMax = arrayMax (aaa) ;
  if (iMax)
    {
      arraySort (aaa, saAlignOrder) ;
      arrayCompress (aaa) ;
      alignDoRegisterOnePair (pp, bb, aaaa, aaa, read1, bestAp1, adaptors) ;
      if (bestAp2) alignDoRegisterOnePair (pp, bb, aaaa, aaa, read2, bestAp2, adaptors) ;
    }
  ac_free (h) ;
} /* alignDoOnePair */

/*
  run -x Aligners/011_SortAlignG5R5/IDX.hs2013.18.81 --maxTargetRepeats 81 -I Fasta/RNA_PolyA_B_1/RNA_PolyA_B_1.Config --align --method 011_SortAlignG5R5  --run RNA_PolyA_B_1  -o RESULTS/011_SortAlignG5R5/RNA_PolyA_B_1/RNA_PolyA_B_1 --step 5 --nRawReads 313433604 --nRawBases 47328474204 --numactl --nB 1
read1 = 144 line=5884
*/
  
/**************************************************************/

void saAlignDo (const PP *pp, BB *bb)
{
  AC_HANDLE h = 0 ;
  BOOL redo = FALSE ;
  int pass = 0 ;
  HIT * restrict hit ;
  HIT *h1, *h2 ;
  long int ii, jj, iMax ;
  BigArray hits = 0 ;
  BigArray hits2 = 0 ;
  Array aa = 0 ;
  Array err = 0 ;
  BigArray aaa = 0 ;
  Array countChroms = 0 ;
  int n = NTARGETREPEATBITS ;
  int mask = (1 << n) - 1 ;
  ADAPTORS adaptors = {{0}} ;

  /* in pilot block, bb->isRna defaults to 0=RNA */
  bb->isRna = 0 ;
  if (1) 
    while (! saSetGetAdaptors (0, &(bb->isRna), &adaptors, bb->run))
      {
	if (bb->lane <= 1) break ; /* pilot block proceeds */
	sleep (1) ;                /* all other blocks wait */
      }

 secondPass:
  iMax = bigArrayMax (bb->hits) ;
  h = ac_new_handle () ;
  hits = bigArrayHandleCreate (256, HIT, h) ;
  hits2 = bigArrayHandleCreate (256, HIT, h) ;  
  aa = arrayHandleCreate (128, ALIGN, h) ;
  err = arrayHandleCreate (256, A_ERR, h) ;
  aaa = bigArrayHandleCreate (iMax, ALIGN, h) ;
  countChroms = arrayHandleCreate (256, COUNTCHROM, h) ;
  bb->confirmedSLs = arrayHandleCreate (6400, SLS, bb->h) ;
  bb->confirmedPolyAs = arrayHandleCreate (64000, POLYA, bb->h) ;
  bb->confirmedIntrons = arrayHandleCreate (64000, INTRON, bb->h) ;
  bb->doubleIntrons = arrayHandleCreate (64000, DOUBLEINTRON, bb->h) ;
  /*
    bb->isAligned = bitSetHandleCreate (bb->nSeqs, bb->h) ;
  */
#ifdef JUNK
  if (1)
    {
      for (ii = 0, hit = iMax ? bigArrp (bb->hits, 0, HIT) : 0 ; ii < iMax ; ii++, hit++)
	fprintf (stderr, "BBHITS\t%ld\t%d\t%d\t%d\t%d\n"
		 , ii
		 , hit->read
		 , hit->chrom
		 , hit->a1
		 , hit->x1
		 ) ;
    }
#endif

  for (ii = 0, hit = iMax ? bigArrp (bb->hits, 0, HIT) : 0 ; ii < iMax ; ii++, hit++)
    {
      int nn = 1, read = hit->read, pair = read >> 1 ;
      for (jj = ii + 1, h1 = hit + 1 ; jj < iMax && (h1->read >> 1) == pair ; jj++, h1++)
	nn++ ;
      if (nn >= 1) /* this read has n+1 hit */
	{ /* create  a copy of the hits of that read */
	  h2 = bigArrayp (hits, nn - 1, HIT) ; /* make room */
	  bigArrayMax (hits) = nn ;
	  h2 = bigArrayp (hits, 0, HIT) ; /* make room */
	  memcpy (h2, hit, nn * sizeof(HIT)) ;

	  /* bb->hits, hence its slice hits, are in hitPairOrder : pair, chrom, position */

	  if (1)
	    {
	      int chrom = 0, mult ;
	      int k, kk = 0, a1 = 0 ;
	      HIT *up  ;
	      COUNTCHROM *zp, *zp0 ;
	      arrayMax (countChroms) = 0 ;

	      /* establish zones */
	      for (k = kk = 0, up = bigArrayp (hits, 0, HIT) ; k < nn ; up++, k++)
		{
		  BOOL isIntron = ((up->x1  >> NTARGETREPEATBITS )  & 0x7) ? TRUE : FALSE ;
		  if (bb->isRna < 0 && isIntron) continue ;
		  if (up->chrom != chrom || up->a1 > a1 + 1000000)
		    {
		      a1 = up->a1 ;

		      zp = arrayp (countChroms, kk++, COUNTCHROM) ;
		      zp->i1 = zp->i2 = k ;
		      zp->a1 = zp->a2 = up->a1 + (up->x1 >> NSHIFTEDTARGETREPEATBITS) ;
		      zp->chrom = chrom = up->chrom ;
		    }
		  else
		    {  /* hits to this chrom are between [zp->a1,zp->a2] */
		      zp->i2 = k ;
		      zp->a2 = up->a1 + (up->x1 >> NSHIFTEDTARGETREPEATBITS) ;
		      a1 = up->a1 ;
		    }
		}

	      /* compute scores as number of seeds in the read x coordinates */
	      if (kk > 0)
		{
		  for (k = 0, zp = arrayp (countChroms, 0, COUNTCHROM)  ; k < kk ; k++, zp++)
		    {
		      int ii, i1 = zp->i1, i2 = zp->i2 + 1 ;
		      int x0 = -999999999 ;

		      zp->seeds = 0 ;
		      zp->seed1 = 0 ;
		      zp->seed2 = 0 ;
		      zp->seed4 = 0 ;
		      zp->seed8 = 0 ;
		      zp->seed16 = 0 ;
		      zp->seed32 = 0 ;
		      zp->weight = 0 ;
		      
		      bigArraySortSlice (hits, i1, i2, hitReadPosOrder) ;
		      for (ii = i1, up = bigArrp (hits, ii, HIT) ; ii < i2 ; ii++, up++)
			{
			  int x1 = up->x1 >> NSHIFTEDTARGETREPEATBITS ;
			  if (x1 > x0)
			    {
			      x0 = x1 ;
			      mult = up->x1 & mask ;
			      /* zp->weight += mult < 4 ? 720/mult : 720/(4 * mult * mult) ; */
			      if (!zp->seeds) zp->x1 = zp->x2 = x1 ;
			      else zp->x2 = x1 ;
			      zp->seeds++ ;
			      if (mult == 1)       { zp->seed1++ ;  zp->weight += 4 ; }
			      else if (mult == 2)  { zp->seed2++ ;  zp->weight += 3 ; }
			      else if (mult <= 4)  { zp->seed4++ ;  zp->weight += 2 ; }
			      else if (mult <= 8)  { zp->seed8++ ;  zp->weight += 1 ; }
			      else if (mult <= 16) { zp->seed16++ ; zp->weight += 1 ; }
			      else if (mult <= 32) { zp->seed32++ ; zp->weight += 1 ; }
			      if (0) zp->weight = zp->seeds ;
			    }
			}
		    }
		  arraySort (countChroms, countChromOrder) ;
		  /* copy the relevant hits */
		  bigArrayMax (hits2) = 0 ;

		  int mm = 0, k8 = 1 ;
		  zp = zp0 = arrayp (countChroms, 0, COUNTCHROM) ;
		  if (zp->seed1 < zp->seed2) k8 = 2 ;
		  if (zp->seed1 < zp->seed4) k8 = 4 ;
		  if (zp->seed1 < zp->seed8) k8 = 8 ;
		  if (zp->seed1 < zp->seed16) k8 = 8 ;
		  if (zp->seed1 < zp->seed32) k8 = 8 ;
		  if (0 && zp->seed1 < 2 * zp->seed8) k8 = 8 ;
		  for (k = 0 ; k < kk ; k++, zp++)
		    {
		      /* keep at most 2 chromosomes */
		      if (
			  ( k > k8 && zp->weight < zp0->weight && zp->seeds < 3) ||
			  (k >= k8 && 3 * zp->weight < zp0->weight) ||
			  (k < k8 && 4 * zp->weight < zp0->weight)
			  )
			break ;
		      for (int i = zp->i1 ; i <= zp->i2 ; i++)
			{
			  HIT *up = bigArrp (hits, i, HIT) ;
			  BOOL isIntron = ((up->x1  >> NTARGETREPEATBITS )  & 0x7) ? TRUE : FALSE ;
			  if (bb->isRna < 0 && isIntron) continue ;
			  HIT *vp = bigArrayp (hits2, mm++, HIT) ;
			  *vp = *up ;
			  vp->chrom ^= (vp->read & 0x1) ;
			}
		    }
		  if (k + 4 < kk) arrayMax(countChroms) = k + 4 ;
		  if (0) arrayMax(countChroms) = 1 ;
		  arrayMax (hits2) = mm ;
		  arrayMax (aa) = arrayMax (err) = 0 ;
		  /* switch chroms and reorder */
		  if (0) hits2 = bb->hits ;
		  bb->gpu += saSort (hits2, 2) ; /* hitReadOrder */
		  if (0)  showCountChroms (countChroms) ;

		  alignDoOnePair (pp, bb, aaa, hits2, aa, err, &adaptors) ;
		}
	    }
	}
      hit += nn - 1 ; ii += nn - 1 ;
    }
  
  bb->aligns = bigArrayHandleCopy (aaa, bb->h) ; /* resize */
  if (bb->lane <= 5)
    {
      /* guess RNA/DNA on the concept that 4000 aligned bases should give us one intron support */
      int isRna = 1 ;
      long int baseAligned = bb->runStat.nBaseAligned1 + bb->runStat.nBaseAligned2 ;
      long int iSupport = bb->runStat.gt_ag_Support + bb->runStat.ct_ac_Support ;
      if (0) iSupport = bb->runStat.nIntronSupportPlus + bb->runStat.nIntronSupportMinus ;
      if (baseAligned > 2000000 && 4000 * iSupport < baseAligned)
	isRna = -1 ;
      if (saReadAdaptors (&adaptors, &(bb->runStat), TRUE))
	saSetGetAdaptors (1, &isRna, &adaptors, bb->run) ;
      else
	saSetGetAdaptors (1, &isRna, 0, bb->run) ;
      if (1)
	{
	  if (! redo && ((isRna < 0 && bb->isRna >= 0) || (isRna > 0 && bb->isRna < 0)))
	    redo = TRUE ;
	  else
	    redo = FALSE ;
	}
      if (0 && redo == 0) redo = 1 ;
      bb->isRna = isRna ;
      fprintf (stderr, "SETGET lane %d isRna=%d isupport= %ld %ld   gc_ag=%d %d  baseAli = %ld\n"
	       , bb->lane, isRna
	       , bb->runStat.nIntronSupportPlus, bb->runStat.nIntronSupportMinus
	       , bb->runStat.gt_ag_Support, bb->runStat.ct_ac_Support
	       ,  bb->runStat.nBaseAligned1 + bb->runStat.nBaseAligned2
	       ) ;
    }
  
  if (1 && pass==0 && redo)
    {
      PSD p = bb->runStat.p ;

      memset (&bb->runStat, 0, sizeof (RunSTAT)) ;
      bb->runStat.p = p ;

      bb->runStat.insertLengthDistribution = arrayHandleCreate (1024, long int, bb->h) ;

      bb->nAli = bb->aliDx = bb->aliDa = 0 ;
      bigArrayMax (aaa) = 0 ;

      pass++ ;
      goto secondPass ;
    }
  ac_free (h) ;
  return ;
} /* saAlignDo */

/**************************************************************/

#ifndef YANN
void saAlign (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  char tBuf[25] ;
  long int nnn = 0 ;
  clock_t  t1, t2 ;

  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->oaChan, &bb, BB))
    {
      if (pp->align && bb.hits)
	{
	  if (pp->debug) printf ("--- %s: Start align %lu seeds\n", timeBufShowNow (tBuf), bigArrayMax (bb.hits)) ;

	  t1 = clock () ;
	  saAlignDo (pp, &bb) ;
	  nnn += bb.nAli ;
	  t2 = clock () ;
	  
	  saCpuStatRegister ("7.Align_r", pp->agent, bb.cpuStats, t1, t2, bb.aligns ? bigArrayMax (bb.aligns) : 0) ;
	  if (pp->debug) printf ("--- %s: Stop align %lu ali, %lu mismatches\n", timeBufShowNow (tBuf), bb.nAli, bb.nerr) ;
	}
      channelPut (pp->aeChan, &bb, BB) ;
    }

  int n = channelCount (pp->plChan) ;
  if (pp->debug) printf ("..... close aeChan at %d,  found %ld ali\n", n, nnn) ;

  channelCloseSource (pp->aeChan) ;

  return ;
} /* saAlign */
#endif

/********************************************************************/
