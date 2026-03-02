/*
 * sa.introns.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module implements all operations related
 * to intron support and discovery
 */

#include "sa.h"

/**************************************************************/
/**************************************************************/

static void showIntrons (Array aa)
{
  if (aa && arrayMax (aa))
    {
      INTRON *zp ;
      int ii, iMax = arrayMax (aa) ;

      for (ii = 0, zp = arrp (aa, 0, INTRON) ; ii < iMax ; ii++, zp++)
	{
	  fprintf (stderr, "%d\trun %d\tchrom %d\ta1 %d\ta2 %d\tn %d\tnR %d\t%s\n"
		   , ii, zp->run, zp->chrom
		   , zp->a1, zp->a2
		   , zp->n, zp->nR, zp->feet
		   ) ;
	}
    }
} /* showIntrons */

/**************************************************************/

static int confirmedPolyAOrder (const void *va, const void *vb)
{
  const POLYA *up = va ;
  const POLYA *vp = vb ;
  int n ;

  /* chrom order */
  n = up->chrom - vp->chrom ; if (n) return n ;
  /* pos order */
  n = up->a1 - vp->a1 ; if (n) return n ;
  /* run order */
  n = up->run - vp->run ; if (n) return n ;
  return 0 ;
} /* confirmedPolyAOrder */

/**************************************************************/

static int confirmedSLsOrder (const void *va, const void *vb)
{
  const SLS *up = va ;
  const SLS *vp = vb ;
  int n ;

  /* chrom order */
  n = up->chrom - vp->chrom ; if (n) return n ;
  /* pos order */
  n = up->a1 - vp->a1 ; if (n) return n ;
  /* run order */
  n = up->run - vp->run ; if (n) return n ;
  return 0 ;
} /* confirmedSLsOrder */

/**************************************************************/

static int confirmedIntronsOrder (const void *va, const void *vb)
{
  const INTRON *up = va ;
  const INTRON *vp = vb ;
  int n ;

  /* chrom order */
  n = up->chrom - vp->chrom ; if (n) return n ;
  /* pos order */
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  /* run order */
  n = up->run - vp->run ; if (n) return n ;
  /* feet order, this clears the nR lines */
  n = up->feet[0] - vp->feet[0] ; if (n) return -n ;
  return 0 ;
} /* confirmedIntronOrder */

/**************************************************************/

static int doubleIntronsOrder (const void *va, const void *vb)
{
  const DOUBLEINTRON *up = va ;
  const DOUBLEINTRON *vp = vb ;
  int n ;

  /* chrom order */
  n = up->chrom - vp->chrom ; if (n) return n ;
  /* pos order */
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  n = up->b1 - vp->b1 ; if (n) return n ;
  n = up->b2 - vp->b2 ; if (n) return n ;
  /* run order */
  n = up->run - vp->run ; if (n) return n ;
  /* feet order, this clears the nR lines */
  n = up->feet1[0] - vp->feet1[0] ; if (n) return -n ;
  return 0 ;
} /* doubleIntronOrder */

/**************************************************************/

static int confirmedPolyAsCompress (Array aaa) 
{
  int ii, jj = 0, kk = 1, iMax = arrayMax (aaa) ;
  arraySort (aaa, confirmedPolyAOrder) ;
  POLYA *up, *vp, *wp ;

  if (iMax)
    for (ii = 0,  jj = 0, up = arrp (aaa, 0, POLYA), vp = up ; ii < iMax ; ii += kk, up += kk)
      {
	for (wp = up + 1, kk = 1 ; ii + kk < iMax ; kk++, wp++)
	  if (up->a1 == wp->a1 && up->chrom == wp->chrom && up->run == wp->run)
	    {
	      up->n += wp->n ; wp->n = 0 ;
	    }
	  else
	    break ;
	if (vp < up) *vp = *up ;
	vp++ ; jj++ ;
      }
  arrayMax (aaa) = jj ;
  return jj ;
} /* confirmedPolyAsCompress */

/**************************************************************/

static int confirmedSLsCompress (Array aaa) 
{
  int ii, jj = 0, kk = 1, iMax = arrayMax (aaa) ;
  arraySort (aaa, confirmedSLsOrder) ;
  SLS *up, *vp, *wp ;

  if (iMax)
    for (ii = 0,  jj = 0, up = arrp (aaa, 0, SLS), vp = up ; ii < iMax ; ii += kk, up += kk)
      {
	for (wp = up + 1, kk = 1 ; ii + kk < iMax ; kk++, wp++)
	  if (up->a1 == wp->a1 && up->chrom == wp->chrom && up->run == wp->run)
	    {
	      up->n += wp->n ; wp->n = 0 ;
	    }
	  else
	    break ;
	if (vp < up) *vp = *up ;
	vp++ ; jj++ ;
      }
  arrayMax (aaa) = jj ;
  return jj ;
} /* confirmedSLsCompress */

/**************************************************************/

static int confirmedIntronsCompress (Array aaa) 
{
  int ii, jj = 0, kk = 1, iMax = arrayMax (aaa) ;
  arraySort (aaa, confirmedIntronsOrder) ;
  INTRON *up, *vp, *wp ;

  if (iMax)
    for (ii = 0,  jj = 0, up = arrp (aaa, 0, INTRON), vp = up ; ii < iMax ; ii += kk, up += kk)
      {
	if (! up->feet[0])
	  continue ;
	for (wp = up + 1, kk = 1 ; ii + kk < iMax ; kk++, wp++)
	  if (up->a1 == wp->a1 && up->a2 == wp->a2 && up->chrom == wp->chrom && up->run == wp->run)
	    {
	      up->n += wp->n ; wp->n = 0 ;
	      up->nR += wp->nR ; wp->nR = 0 ;
	    }
	  else
	    break ;
	if (vp < up) *vp = *up ;
	vp++ ; jj++ ;
      }
  arrayMax (aaa) = jj ;
  return jj ;
} /* confirmedIntronsCompress */

/**************************************************************/

static int doubleIntronsCompress (Array aaa) 
{
  int ii, jj = 0, kk = 1, iMax = arrayMax (aaa) ;
  arraySort (aaa, doubleIntronsOrder) ;
  DOUBLEINTRON *up, *vp, *wp ;

  if (iMax)
    for (ii = 0,  jj = 0, up = arrp (aaa, 0, DOUBLEINTRON), vp = up ; ii < iMax ; ii += kk, up += kk)
      {
	if (! up->feet1[0])
	  continue ;
	for (wp = up + 1, kk = 1 ; ii + kk < iMax ; kk++, wp++)
	  if (up->a1 == wp->a1 && up->a2 == wp->a2 &&
	      up->b1 == wp->b1 && up->b2 == wp->b2 &&
	      up->chrom == wp->chrom && up->run == wp->run)
	    {
	      up->n += wp->n ; wp->n = 0 ;
	      up->nR += wp->nR ; wp->nR = 0 ;
	    }
	  else
	    break ;
	if (vp < up) *vp = *up ;
	vp++ ; jj++ ;
      }
  arrayMax (aaa) = jj ;
  return jj ;
} /* doubleIntronsCompress */

/**************************************************************/

void saPolyAsCumulate (PP *pp, BB *bb) 
{
  Array aaa = pp->confirmedPolyAs, aa = bb->confirmedPolyAs ;
  int iMax = aa ? arrayMax (aa) : 0 ;

  if (iMax)
    {
      if (! aaa)
	aaa = pp->confirmedPolyAs = arrayHandleCopy (bb->confirmedPolyAs, pp->h) ;
      else
	{
	  int jj = arrayMax (aaa)  ;
	  arrayp (aaa, jj + iMax - 1, POLYA)->n = 0 ; /*  make room */
	  memcpy (arrp (aaa, jj, POLYA), arrp (aa, 0, POLYA), iMax * sizeof (POLYA)) ;
	  confirmedPolyAsCompress (aaa) ;
	}
    }
  return ;
} /* saPolyAsCumulate */

/**************************************************************/

void saSLsCumulate (PP *pp, BB *bb) 
{
  Array aaa = pp->confirmedSLs, aa = bb->confirmedSLs ;
  int iMax = aa ? arrayMax (aa) : 0 ;

  if (iMax)
    {
      if (! aaa)
	aaa = pp->confirmedSLs = arrayHandleCopy (bb->confirmedSLs, pp->h) ;
      else
	{
	  int jj = arrayMax (aaa)  ;
	  arrayp (aaa, jj + iMax - 1, POLYA)->n = 0 ; /*  make room */
	  memcpy (arrp (aaa, jj, POLYA), arrp (aa, 0, POLYA), iMax * sizeof (POLYA)) ;
	  confirmedSLsCompress (aaa) ;
	}
    }
  return ;
} /* saSLsCumulate */

/**************************************************************/

void saIntronsCumulate (PP *pp, BB *bb) 
{
  Array aaa = pp->confirmedIntrons, aa = bb->confirmedIntrons ;
  int iMax = arrayMax (aa) ;
  
  if (iMax)
    {
      int jj = arrayMax (aaa) ;
      arrayp (aaa, jj + iMax - 1, INTRON)->n = 0 ; /*  make room */
      memcpy (arrp (aaa, jj, INTRON), arrp (aa, 0, INTRON), iMax * sizeof (INTRON)) ;
      confirmedIntronsCompress (aaa) ;
    }
  return ;
} /* saIntronsCumulate */

/**************************************************************/

void saDoubleIntronsCumulate (PP *pp, BB *bb) 
{
  Array aaa = pp->doubleIntrons, aa = bb->doubleIntrons ;
  int iMax = arrayMax (aa) ;
  
  if (iMax)
    {
      int jj = arrayMax (aaa) ;
      arrayp (aaa, jj + iMax - 1, DOUBLEINTRON)->n = 0 ; /*  make room */
      memcpy (arrp (aaa, jj, DOUBLEINTRON), arrp (aa, 0, DOUBLEINTRON), iMax * sizeof (DOUBLEINTRON)) ;
      doubleIntronsCompress (aaa) ;
    }
  return ;
} /* saDoubleIntronsCumulate */

/**************************************************************/

int saSupportedIntrons (const PP *pp, int run)
{
  int nn = 0 ;
  Array aa = pp->confirmedIntrons ;
  int ii, iMax = arrayMax (aa) ;
  INTRON *up ;

  if (iMax)
    for (ii = 0, up = arrp (aa, 0, INTRON) ; ii < iMax ; ii++, up++)
      if (! run || up->run == run) nn++ ;

  return nn ;
} /* saSupportedIntrons */

/**************************************************************/
/* sliding intron, clip errors in overlap r636 r452 */
void saIntronsOptimize (BB *bb, ALIGN *vp, ALIGN *wp, Array dnaG)  
{
  A_ERR *epX, *epY ;
  int x2 = vp->x2 ;
  int y1 = wp->x1 ;
  int nEx = vp->nErr ;
  int nEy = wp->nErr ;
  int i, j, nE ;
  int bestN, bestI = -1, bestJ = -1 ;
  int dy = vp->x2 - wp->x1 + 1 ;  /* dy > 0 recouvrement, dy < 0: trou dans le read */
  int da = vp->a2 - wp->a1 ;      
  BOOL isDown = (vp->a1 < vp->a2 ? TRUE : FALSE ) ;
  if (da < 5 && da > -5 && dy < 5 && dy > -5 && vp->chrom == wp->chrom)
    {
      /* merge the 2 alignments */
      vp->a2 = wp->a2 ;
      vp->x2 = wp->x2 ;
      vp->nErr += wp->nErr ;
      vp->errors = vp->errors ? vp->errors : wp->errors ;
      if (vp->errors && vp->nErr)
	epX = arrayp (vp->errors, vp->nErr - 1, A_ERR) ; /* the extra errors from wp are not copied, hence set to zero */
      if (! vp->nErr)
	{
	  vp->nErr = 1 ;
	  vp->errors = arrayHandleCreate (20, A_ERR, bb->h) ;
	  epX = arrayp (vp->errors, 0, A_ERR) ;
	  epX->iShort = wp->x1 ;
	  epX->iLong = wp->a1 ;
	  epX->type = ERREUR ;
	}	  
      memset (wp, 0, sizeof (ALIGN)) ; nEy = 0 ;
      wp->chain = -1 ;
      return ;
    }
  if (vp->x1 > wp->x1)
    {
      int dx = vp->x1 - wp->x1 ;
      wp->x1 = vp->x1 ;
      wp->a1 += (isDown ? dx : -dx ) ;  /* correct only if there are no indel in [wp->a1,wp->a1+dx] segment */
    }
  dy = vp->x2 - wp->x1 + 1 ;
  if (dy > 0 && nEx + nEy > 0)
    {
      int zX =x2 , zY ;
      int zEx = 0, zEy = nEy ;
      epX = epY = 0 ; bestI = 0 ; bestJ = -1 ; nE = 0 ; 
      for (i = 0 ; i < nEx ; i++)
	{  /* count the vp errors that cannot be clipped */
	  epX = arrp (vp->errors, i, A_ERR) ;
	  int zX = epX->iShort ; /* last good base */
	  if (zX  < y1 - 1) /* bio coords */
	    { bestI = i + 1 ; zEx++ ; }
	  else
	    break ;
	}
      nE = zEx + zEy ;
      bestN = nE ;
      int cI = bestI, cJ = bestJ ;
      int cY1 = y1, cX2 = zX ;
      
      /* if we clip both ali at bestY1 <= z <= bestX2, we keep nE errors
       * let us try to find a better position
       */
      while (TRUE)
	{
	  epX = (cI >=0 && cI < nEx ? arrp (vp->errors, cI, A_ERR) : 0) ;
	  epY = (cJ >=0 && cJ < nEy ? arrp (wp->errors, cJ, A_ERR) : 0) ;
	  if (epX && epY && epX->iLong == epY->iLong && epX->iShort == epY->iShort)
	    {
	      /* merge the 2 alignments */
	      vp->a2 = wp->a2 ;
	      vp->x2 = wp->x2 ;
	      cI++ ;
	      for (int i = cJ + 1 ; i < nEy ; i++, cI++)
		array (vp->errors, cI++, A_ERR) = array (wp->errors, i, A_ERR) ;
	      vp->nErr = arrayMax (vp->errors) = i ;
	      memset (wp, 0, sizeof (ALIGN)) ;
	      wp->chain = -1 ; 
	      return ;
	    }

	  /* run -x titi4/IDX -i titi4/s9324.fasta -o titi4/S --debug */
	  /* find the next zY (candidate bestY1) */
	  zY = cY1 ;
	  j = cJ + 1 ;
	  if (j < nEy)
	    {
	      epY = arrp (wp->errors, j, A_ERR) ;
	      zY = epY->iShort + 2  ; /* first good base starting from the right */
	      switch (epY->type)
 		{   /* we need to adjust the coordinates */
		case INSERTION:
		  break ;
		case INSERTION_DOUBLE:
		  zY+= 1 ;
		  break ;
		case INSERTION_TRIPLE:
		  zY += 2 ;
		  break ;
		case TROU:
		case TROU_DOUBLE:
		case TROU_TRIPLE:
		  zY -- ;
		  break ;
		default:
		  break ;
		}

	      if (zY <= cX2 + 1)
		{  /* eating this Y error is favorable */
		  zEy-- ;
		  nE = zEx + zEy ;
		  if (nE < bestN) { bestN = nE ; bestJ = j ; bestI = i ; }
		  cJ = j ; cY1 = zY ;
		  continue ;
		}
	    }	  
	  /* advance on zX */
	  if (cJ >= 0 && cJ < nEy)
	    {
	      zX = cX2 ;
	      i = cI + 1 ;
	      if (i < nEx)
		{
		  epX = arrp (vp->errors, i, A_ERR) ;
		  zX = epX->iShort ; /* last good base starting from the left */
		  /* not favorable but maybe we canadvance on Y */
		  zEx = i - 1 ; cI++ ; cX2 = zX ; continue ;
		}
	    }
	  /* we cannot advance on X and we already have advanced on Y, so we are done */
	  break ;
	}

      if (nEx && bestI <= nEx - 1)  /* we can clip the vp errors */
	{
	  if (bestI == -1) bestI = 0 ;
	  if (bestI >= 0)
	    {
	      epX = arrp (vp->errors, bestI, A_ERR) ;
	      vp->x2 = epX->iShort ;  /* last exact base bio coords */
	      /* vp->a2 = epX->iLong + (isDown ? 0 : 2) ; */
	      int zA = epX->iLong + (isDown ? 0 : -1) ;
	      vp->a2 = (vp->chrom & 0x1 ?  arrayMax(dnaG) - zA + 1 : zA) ;
	    }
	  nEx = bestI ;
	  vp->nErr = arrayMax (vp->errors) = nEx ;
	}

      if (nEy && bestJ >= 0) /* if bestJ == -1, wp->x1/a1 is well positioned */
	{ /* we must clip these errors */
	  epY = arrp (wp->errors, bestJ, A_ERR) ;
	  int zA = epY->iLong  ;
	  int zY = epY->iShort ; /* we need to find the first exact */
	  int dA = 0, dY = 0 ;
	  switch (epY->type)
	    {   
	    case INSERTION:
	      dA-- ;
	      break ;
	    case INSERTION_DOUBLE:
	      dY++ ;
	      dA-- ;
	      break ;
	    case INSERTION_TRIPLE:
	      dY += 2 ;
	      dA-- ;
	      break ;
	    case TROU:
	      dY -- ;
	      break ;
	    case TROU_DOUBLE:
	      dY-- ;
	      dA++ ;
	      break ;
	    case TROU_TRIPLE:
	      dY -- ;
	      dA += 2 ;
	      break ;
	    default:
	      break ;
	    }
	  wp->x1 = zY + dY + 1 ;  /* bio coords */
	  wp->a1 = (isDown ?  zA + dA + 1 : arrayMax(dnaG) - zA - dA) ;
	  wp->nErr -= bestJ + 1 ;
	  if (wp->nErr)
	    for (j = 0, epY = arrp (wp->errors, 0, A_ERR) ; j < wp->nErr ; epY++, j++)
	      *epY = *(epY + bestJ + 1) ;
	  arrayMax (wp->errors) = wp->nErr ;
	}
    }
  /*
  111M8590N187M696N213M1I21M1I 28M1S1071N 3M1I45M1I25M
  111M8590N187M696N213M1I21M1I 28M1I1071N     48M1I25M
  111M8590N187M696N210M1I21M1I 31M1071N   1I45M1I28M
  */
  BOOL foundDonor = FALSE ;
  BOOL foundAcceptor = FALSE ;
  int donor = 0, acceptor = 0 ;
  BOOL isReadDown = TRUE ;
  dy = vp->x2 - wp->x1 + 1 ;
  
  if (bb->isRna < 0)
    goto done ;

  /* trim if possible the vp->x2 end of the first exon on a known 'donor' */
  if (vp->a1 < vp->a2)
    {
      isReadDown = TRUE ;
      donor = vp->donor ;
      acceptor = wp->acceptor ;
      if (donor < 0)
	{ donor = - donor ; acceptor = -acceptor ; }
      else if (donor == 0 && acceptor < 0)
	{ acceptor = - acceptor ; }
      if (donor > acceptor)
	donor = acceptor = 0 ;
      
      if (dy > 0 && donor && vp->a2 >= donor && vp->a2 <= donor + dy - 1) 
	{  /* move back to the canonical donor site */
	  dy = vp->a2 - donor + 1 ;
	  vp->x2 -= dy ;
	  vp->a2 -= dy ;
	}
      if (donor && vp->a2 == donor - 1)
	foundDonor = TRUE ;
    }
  else
    {
      isReadDown = FALSE ;
      donor = vp->donor ;
      acceptor = wp->acceptor ;

      if (donor < 0)
	{ donor = - donor ; acceptor = -acceptor ; }
      if (donor == 0 && acceptor < 0)
	{ acceptor = - acceptor ; }
      if (donor && donor < acceptor)
	donor = acceptor = 0 ;
      
      if (dy > 0 && donor && vp->a2 <= donor && vp->a2 >= donor - dy + 1)
	
	{  /* move back to the canonical donor site */
	  dy = donor - vp->a2 + 1 ;
	  vp->x2 -= dy ;
	  vp->a2 += dy ;
	}      
      if (donor && vp->a2 == donor + 1)
	foundDonor = TRUE ;
    }

  /* alternativelly  trim if possible the wp->x2 start of the second  exon on a known 'acceptor' */
  dy = vp->x2 - wp->x1 + 1 ;
  if (isReadDown && acceptor)
    {
      if (dy > 0 && wp->a1 < acceptor + 1 && wp->a1 + dy >= acceptor + 1)
	{
	  dy = acceptor + 1 - wp->a1 ;
	  wp->x1 += dy ;
	  wp->a1 += dy ;
	}
      if (wp->a1 == acceptor + 1)
	{
	  foundAcceptor = TRUE ;
	  dy = vp->x2 - wp->x1 + 1 ;
	  if (dy > 0) /* trim vp->x2 */
	    {
	      vp->x2 -= dy ;
	      vp->a2 -= dy ;
	    }
	}
    }
  else if (!isReadDown && acceptor)
    {
      if (dy > 0 && wp->a1 > acceptor - 1 && wp->a1 - dy <= acceptor - 1)
	{
	  dy = wp->a1 - acceptor + 1 ;
	  wp->x1 += dy ;
	  wp->a1 -= dy ;
	}
      if (wp->a1 == acceptor - 1)
	{
	  foundAcceptor = TRUE ;
	  dy = vp->x2 - wp->x1 + 1 ;
	  if (dy > 0) /* trim vp->x2 */
	    {
	      vp->x2 -= dy ;
	      vp->a2 += dy ;
	    }
	}
    }

  /* no annotated junction, look for gt_ag */
  dy = vp->x2 - wp->x1 + 1 ;
  if (dy > 0 && isReadDown && ! foundDonor && !foundAcceptor  && wp->a1 > 3 && vp->a2 > dy)  
    {
      unsigned char *cp = arrp (dnaG, vp->a2 - dy, unsigned char) ; /* the base just after vp->a2 - dy */
      unsigned char *cq = arrp (dnaG, wp->a1 - 3, unsigned char) ; /* the base 2 bases before wp->a1 - dy */
      int bestI = -1 ;

      if (bb->runStat.gt_ag_Support >= bb->runStat.ct_ac_Support)
	{     /* favor gt_ag over ct_ac */
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == C_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == C_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == G_ && cq[i+1] == C_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok1 ; }
	}
      else
	{ /* favor ct_ac over gt_ag */
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == C_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == G_ && cq[i+1] == C_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == C_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cq[i] == A_ && cq[i+1] == C_)
	      { bestI = i ; goto ok1 ; }
	}
    ok1:
      if (bestI >= 0)
	{
	  dy -= bestI ;
	  vp->x2 -= dy ;
	  vp->a2 -= dy ;
	  donor = vp->a2 + 1 ;
	  
	  dy = bestI ;
	  wp->x1 += dy ;
	  wp->a1 += dy ;
	  acceptor = wp->a1 - 1 ;
	}	    
    }
  else if (dy > 0 && !isReadDown && ! foundDonor && !foundAcceptor && wp->a1 > dy && vp->a2 > 3) 
    {
      /* move backwards on the genome */
      unsigned char *cp = arrp (dnaG, wp->a1 - dy, unsigned char) ; /* the base just after vp->a2 - dy */
      unsigned char *cq = arrp (dnaG, vp->a2 - 3, unsigned char) ; /* the base 2 bases before wp->a1 - dy */
      int bestI = -1 ;
      if (bb->runStat.gt_ag_Support < bb->runStat.ct_ac_Support)
	{     /* favor gt_ag over ct_ac */
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == C_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == C_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == G_ && cq[i+1] == C_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok2 ; }
	}
      else
	{ /* favor ct_ac over gt_ag */
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == C_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == G_ && cq[i+1] == C_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == C_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cq[i] == A_ && cq[i+1] == C_)
	      { bestI = i ; goto ok2 ; }
	}
    ok2:
      if (bestI >= 0)
	{
	  dy -= bestI ;
	  wp->x1 += dy ;
	  wp->a1 -= dy ;
	  donor = vp->a2 - 1 ;
	  
	  dy = bestI ;
	  vp->x2 -= dy ;
	  vp->a2 += dy ;
	  acceptor = wp->a1 + 1 ;
	}	    
    }
  
 done:   /* In any case trim the second exon a bouts francs */
  dy = vp->x2 - wp->x1 + 1 ;
  if (dy > 0)
    { 
      wp->x1 += dy ;
      wp->a1 += (wp->a1 < wp->a2 ? dy : -dy) ;
    }

  dy = vp->x2 - wp->x1 + 1 ;
  da = vp->a2 - wp->a1 ;

  if (da < 5 && da > -5 && dy < 5 && dy > -5 && vp->chrom == wp->chrom)
    {
      /* merge the 2 alignments */
      vp->a2 = wp->a2 ;
      vp->x2 = wp->x2 ;
      vp->nErr += wp->nErr ;
      vp->errors = vp->errors ? vp->errors : wp->errors ;
      if (vp->errors && vp->nErr)
	epX = arrayp (vp->errors, vp->nErr - 1, A_ERR) ;
      memset (wp, 0, sizeof (ALIGN)) ;
      wp->chain = -1 ;
      return ;
    }
} /* saIntronsOptimize */

/**************************************************************/
/**************************************************************/

void saPolyAsExport (PP *pp, Array aaa)
{
  int iMax = aaa ? confirmedPolyAsCompress (aaa) : 0 ;
  
  if (iMax)
    {
      AC_HANDLE h = ac_new_handle () ;
      ACEOUT ao = aceOutCreate (pp->outFileName, ".polyA.tsf", 0, h) ;
      POLYA *up ;
      long int ii ;
      
      aceOutf (ao, "### PolyAd support in tsf format: chrom__a1,  run, i, strand, nb of reads supportingt the polyA\n") ;
      aceOutf (ao, "### Call bin/tsf -i %s -I tsf -O table -o my_table.txt to reformat this file into an excell compatible tab delimited table\n",
	       aceOutFileName (ao)
	       ) ;
      aceOutf (ao, "PolyA\tRun\tti\tStrand\tn\n") ;
      for (ii = 0, up = arrp (aaa, ii, POLYA) ; ii < iMax ; ii++, up++)
	{
	  int min = 1 ;
	  if (up->n >= min)
	    aceOutf (ao, "%s__%d\t%s\tti\t%s\t%d\n"
		     , dictName (pp->bbG.dict, up->chrom >> 1) 
		     , up->a1
		     , up->run ? dictName (pp->runDict, up->run) : "runX"
		     , up->chrom & 0x1 ? "-" : "+"
		     , up->n
		     ) ;
	}
      
      ac_free (h) ;
    }
} /* saPolyAsExport */

/**************************************************************/
/**************************************************************/

void saSLsExport (PP *pp, Array aaa)
{
  int iMax = aaa ? confirmedSLsCompress (aaa) : 0 ;
  
  if (iMax)
    {
      AC_HANDLE h = ac_new_handle () ;
      ACEOUT ao = aceOutCreate (pp->outFileName, ".SL.tsf", 0, h) ;
      SLS *up ;
      long int ii ;
      
      aceOutf (ao, "### SL support in tsf format: chrom__a1,  run, i, strand, nb of reads supportingt the polyA\n") ;
      aceOutf (ao, "### Call bin/tsf -i %s -I tsf -O table -o my_table.txt to reformat this file into an excell compatible tab delimited table\n",
	       aceOutFileName (ao)
	       ) ;
      aceOutf (ao, "SL\tRun\tti\tStrand\tn\n") ;
      for (ii = 0, up = arrp (aaa, ii, SLS) ; ii < iMax ; ii++, up++)
	{
	  int min = 1 ;
	  if (up->n >= min)
	    aceOutf (ao, "%s__%d___SL%d\t%s\tti\t%s\t%d\n"
		     , dictName (pp->bbG.dict, up->chrom >> 1) 
		     , up->a1, up->run  & 0xf
		     , up->run >> 4 ? dictName (pp->runDict, up->run >> 4) : "runX"
		     , up->chrom & 0x1 ? "-" : "+"
		     , up->n
		     ) ;
	}
      
      ac_free (h) ;
    }
} /* saSLsExport */

/**************************************************************/
/**************************************************************/

static char *flipFeet (char *feet)
{
  char buf[6] ;
  memcpy (buf, feet, 6) ;
  feet[0] = complementLetter(buf[4]) ;
  feet[1] = complementLetter(buf[3]) ;
  feet[3] = complementLetter(buf[1]) ;
  feet[4] = complementLetter(buf[0]) ;

  return feet ;
}
  
/**************************************************************/

void saIntronStranding (PP *pp, Array aa)
{
  INTRON *zp, *zpR ;
  int runMax = dictMax (pp->runDict) + 1 ;
  int nGt_ag[runMax] ;
  int nCt_ac [runMax] ;
  int cGood[runMax][12] ;
  int cOther [runMax][12] ;
  int k  ;
  float minS = 100 ;
  int run, ii, iMax = arrayMax (aa) ;
  float *s0 = pp->runStranding ;
  
  memset (nGt_ag, 0, sizeof (nGt_ag)) ;
  memset (nCt_ac, 0, sizeof (nCt_ac)) ;
  memset (cGood, 0, sizeof (cGood)) ;
  memset (cOther, 0, sizeof (cOther)) ;
  
  if (iMax)
    for (ii = 0, zp = arrp (aa, ii, INTRON) ; ii < iMax ; ii++, zp++)
      {
	for (char *cp = zp->feet ; *cp ; cp++)
	*cp = ace_lower ((int)*cp) ;
	if (! strcmp (zp->feet, "gt_ag"))
	  {
	    k = zp->n + zp->nR ;
	    nGt_ag[zp->run]++ ;
	    cGood[zp->run][k < 12 ? k : 11]++ ;
	  }
	else if (! strcmp (zp->feet, "ct_ac"))
	  {
	    k = zp->n + zp->nR ;
	    nCt_ac[zp->run]++ ;
	    cGood[zp->run][k < 12 ? k : 11]++ ;
	  }
	else
	  {
	    k = zp->n + zp->nR ;
	    cOther[zp->run][k < 12 ? k : 11]++ ;
	  }
      }
  
  /* check ratio of Good/Other to find the threshold */
  for (run = 0 ; run < runMax ; run++)
    {
      /* cumul by coverage */
      for (k = 1 ; k < 12 ; k++)
	{	
	  cGood[run][k] += cGood[run][k-1] ;
	  cOther[run][k] += cOther[run][k-1] ;
	}
    }
  for (run = 0 ; run < runMax ; run++)
    {
      s0[run] = array(pp->runStats, run, RunSTAT).intronStranding = 100.0 * (nGt_ag[run] + 1.0) / (nGt_ag[run] + nCt_ac[run] + .0001) ;
      
      if (pp->strand)
	s0[run] = 100 ;
      else if (pp->antiStrand)
	s0[run] = 0 ;
      if (s0[run] < minS && nCt_ac[run])
	minS = s0[run] ;
    }
  
  if (minS < 70) /* flip needed */
    {
      for (ii = 0, zp = arrp (aa, ii, INTRON) ; ii < iMax ; ii++, zp++)
	{
	  if (s0[zp->run] < 40)
	    {  /* flip the whole run */
	      int a0 = zp->a1 ; zp->a1 = zp->a2 ; zp->a2 = a0 ;
	      flipFeet (zp->feet) ;
	    }
	  else if (s0[zp->run] < 60)  
	    {  /* non stranded case, choose for every intron */
	      if (! strcmp (zp->feet, "ct_ac") || ! strcmp (zp->feet, "ct_gc"))
		{
		  int a0 = zp->a1 ; zp->a1 = zp->a2 ; zp->a2 = a0 ;
		  flipFeet (zp->feet) ;
		}
	    }
	}
    }
  
  /* compute the anti counts */
  if (iMax)
    {
      array (aa, 2*iMax -1, INTRON).n = 0 ;
      for (ii = 0, zp = arrp (aa, ii, INTRON), zpR = zp + iMax ; ii < iMax ; ii++, zp++, zpR++)
	{
	  *zpR = *zp ;
	  zpR->a1 = zp->a2 ; zpR->a2 = zp->a1 ; 
	  zpR->nR = zp->n ; zpR->n = zp->nR ; zpR->feet[0] = 0 ;
	}
      /* merged counts and antiCounts */
      for (run = 0 ; run < runMax ; run++)
	{
	  array(pp->runStats, run, RunSTAT).nIntronSupportPlus = 0 ;
	  array(pp->runStats, run, RunSTAT).nIntronSupportMinus = 0 ;
	}
      iMax = confirmedIntronsCompress (aa) ;
      for (ii = 0, zp = arrp (aa, ii, INTRON); ii < iMax ; ii++, zp++)
	{
	  array(pp->runStats, zp->run, RunSTAT).nIntronSupportPlus +=  (s0[zp->run] < 40 ? zp->nR : zp->n) ;
	  array(pp->runStats, zp->run, RunSTAT).nIntronSupportMinus += (s0[zp->run] < 40 ? zp->n : zp->nR) ;
	  array(pp->runStats, 0, RunSTAT).nIntronSupportPlus += (s0[zp->run] < 40 ? zp->nR : zp->n) ;
	  array(pp->runStats, 0, RunSTAT).nIntronSupportMinus += (s0[zp->run] < 40 ? zp->n : zp->nR) ;
	}
      for (run = 0 ; run < runMax ; run++)
	{
	  int np = array(pp->runStats, run, RunSTAT).nIntronSupportPlus ;
	  int nm = array(pp->runStats, run, RunSTAT).nIntronSupportMinus ;
	  array(pp->runStats, run, RunSTAT).intronStranding = 100.0 * np / (np + nm + 0.0000001) ;
	}
    }
  return ;
} /* saIntronStranding */

/**************************************************************/

void saIntronsExport (PP *pp, Array aaa)
{
  int iMax = aaa ? confirmedIntronsCompress (aaa) : 0 ;
  showIntrons (0) ; /* for compiler happiness */
  
  if (iMax)
    {
      AC_HANDLE h = ac_new_handle () ;
      ACEOUT ao = aceOutCreate (pp->outFileName, ".introns.tsf", 0, h) ;
      INTRON *up ;
      long int ii ;
      
      aceOutf (ao, "### Intron support in tsf format: chrom__a1_a2,  run, iit (format for 2 integers 1 text), nb of reads supportingt the intron, number antistrand, intron feet\n") ;
      aceOutf (ao, "### Call bin/tsf -i %s -I tsf -O table -o my_table.txt to reformat this file into an excell compatible tab delimited table\n",
	       aceOutFileName (ao)
	       ) ;
      aceOutf (ao, "#Intron\tRun\tformat\tn\tnR\tfeet\n") ;
      for (ii = 0, up = arrp (aaa, ii, INTRON) ; ii < iMax ; ii++, up++)
	{
	  int min = 3 ;
	  
	  if (! up->feet[0])
	    continue ;
	  if (!strcasecmp (up->feet, "gt_ag"))
	   min = 1 ;
	  else if (!strcasecmp (up->feet, "gc_ag"))
	    min = 2 ;

	  if (up->n + up->nR >= min)
	    aceOutf (ao, "%s__%d_%d\t%s\tiit\t%d\t%d\t%s\n"
		     , dictName (pp->bbG.dict, up->chrom >> 1) + 2, up->a1, up->a2
		     , dictName (pp->runDict, up->run)
		     , up->n, up->nR
		     , up->feet
		     ) ;
	}
      
      ac_free (h) ;
    }
} /* saIntronsExport */

/**************************************************************/

void saDoubleIntronsExport (PP *pp, Array aaa)
{
  int iMax = aaa ? doubleIntronsCompress (aaa) : 0 ;
  float *s0 = pp->runStranding ;
  
  if (iMax)
    {
      AC_HANDLE h = ac_new_handle () ;
      ACEOUT ao = aceOutCreate (pp->outFileName, ".double_introns.tsf", 0, h) ;
      DOUBLEINTRON *up ;
      long int ii ;
      

      aceOutf (ao, "### Double_intron support in tsf format: chrom__a1_a2__b1_b2,  run, iitt (format for 2 integers 2 texts), nb of reads supportingt the intron, number antistrand, intron feet\n") ;
      aceOutf (ao, "### Call bin/tsf -i %s -I tsf -O table -o my_table.txt to reformat this file into an excell compatible tab delimited table\n",
	       aceOutFileName (ao)
	       ) ;
      aceOutf (ao, "Double_intron\tRun\tiit\tn\tnR\tfeet1\tfeet2\n") ;
      for (ii = 0, up = arrp (aaa, ii, DOUBLEINTRON) ; ii < iMax ; ii++, up++)
	{
	  int min = 2 ;
	  if (! up->feet1[0])
	    continue ;
	  if (s0[up->run] < 40 || /* flip whole run */
	      (s0[up->run] < 40 && (! strcasecmp (up->feet1, "ct_ac") || ! strcasecmp (up->feet2, "ct_ac")))
	      ) /* flip needed */
	    {
	      if (!strcasecmp (up->feet1, "ct_ac") && !strcasecmp (up->feet2, "ct_ac"))
		min = 1 ;
	      else if (!strcasecmp (up->feet1, "ct_gc") && !strcasecmp (up->feet2, "ct_ac"))
		min = 1 ;
	      else if (!strcasecmp (up->feet1, "ct_ac") && !strcasecmp (up->feet2, "ct_gc"))
		min = 1 ;
	      if (up->n >= min)
		aceOutf (ao, "%s__%d_%d___%d_%d\t%s\tiitt\t%d\t%d\t%s\t%s\n"
			 , dictName (pp->bbG.dict, up->chrom >> 1) + 2, up->b2, up->b1, up->a2, up->a1
			 , dictName (pp->runDict, up->run)
			 , up->n, up->nR
			 , flipFeet (up->feet2)
			 , flipFeet (up->feet1)
			 ) ;
	    }
	  else 
	    {
	      if (!strcasecmp (up->feet1, "gt_ag") && !strcasecmp (up->feet2, "gt_ag"))
		min = 1 ;
	      else if (!strcasecmp (up->feet1, "gc_ag") && !strcasecmp (up->feet2, "gt_ag"))
		min = 1 ;
	      else if (!strcasecmp (up->feet1, "gt_ag") && !strcasecmp (up->feet2, "gc_ag"))
		min = 1 ;
	      if (up->n >= min)
		aceOutf (ao, "%s__%d_%d___%d_%d\t%s\tiitt\t%d\t%d\t%s\t%s\n"
			 , dictName (pp->bbG.dict, up->chrom >> 1) + 2, up->a1, up->a2, up->b1, up->b2
			 , dictName (pp->runDict, up->run)
			 , up->n, up->nR
			 , up->feet1
			 , up->feet2
			 ) ;
	    }
	}
      
      ac_free (h) ;
    }
} /* saDoubleIntronsExport */

/**************************************************************/
/**************************************************************/
/**************************************************************/

