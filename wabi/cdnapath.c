/*  File: cdnapath.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg. 1997
 *-------------------------------------------------------------------
 * This file is part of the ACEMBLY extension to 
        the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
   Align cDNA
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  9 15:18 1998 (fw)
 * Created: Thu Dec  9 00:01:50 1997 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */
#define CHRONO


#include "freeout.h"
#include "keyset.h"
#include "cdna.h"
#include "cdnapath.h"
#include "bitset.h"
#include "query.h"
static BOOL isFinalPath = FALSE ;

static void* myMagic = 0 ;

typedef struct cPathStruct *PP ;
typedef struct cPathStruct { Array hits, values ; BitSet bb, dead, no ; int searchRepeats ;
                             long (*oneValue) (PP pp, int iOld, int ii, int *offSetCostp) ;
                             int (*isCompatible) (PP pp, int ii, int jj) ;
                             char *title ; int a1, a2, nh, pMax ; void **magic ; } PPstruct ;
typedef struct cPathValueStruct {int p ; long value ; int offSetCost ; } VALUE ;

#define cPathDestroy(_pp) (cPathDoDestroy(_pp), pp=0)
static void cPathMakeOne (PP pp, int ii) ;
static void cPathShow (PP pp) ;

/***********************************************************************/

static PP cPathCreate (char *title, Array hits, 
		       long (*oneValue) (PP pp, int iOld, int ii, int *offSetCostp), 
		       int (*isCompatible) (PP pp, int ii, int jj))
{
  PP pp = (PP) messalloc (sizeof (struct cPathStruct)) ;
  
  pp->magic = &myMagic ;
  pp->title = title ;
  pp->hits = hits ;
  pp->values = 0 ;
  pp->nh = arrayMax(hits) ;
  pp->bb = bitSetCreate (2 * pp->nh, 0) ;
  pp->dead = bitSetCreate (2 * pp->nh, 0) ;
  pp->no = 0 ;

  pp->oneValue = oneValue ;
  pp->isCompatible = isCompatible ;
  pp->pMax = 0 ;

  return pp ;
} /* cPathCreate  */

/***********************************************************************/

static void cPathDoDestroy (PP pp)
{
  if (!pp)
    return ;
  if (pp->magic != &myMagic)
    messcrash ("Bad magic in cPathDestroy") ;
  pp->magic = 0 ;
  
  pp->hits = 0 ;  /* do not destroy pp->hits, not mine */
  arrayDestroy (pp->values) ;
  bitSetDestroy (pp->bb) ;
  bitSetDestroy (pp->dead) ;
  bitSetDestroy (pp->no) ;
  messfree (pp) ;
} /* cPathDoDestroy */

/***********************************************************************/

static void cPathMakeAll (PP pp)
{
  int ii = pp->nh ;

  while (ii--)
    cPathMakeOne (pp, ii) ;
} /* cPathMakeAll */

/***********************************************************************/
/***********************************************************************/
static int MAXOFFSETCOST = 360000 ; /* 12 error */
static int OFFSETCOST =  10 ;

static int cPathValueDecreassingOrder (const void *a, const void *b)
{
  const VALUE *va = (const VALUE *)a ;
  const VALUE *vb = (const VALUE *)b ;
  
  long dv = va->value - vb->value ;
  int doff = va->offSetCost - vb->offSetCost ;

  if (doff > MAXOFFSETCOST)
    doff = MAXOFFSETCOST ;
  if (doff < - MAXOFFSETCOST)
    doff = - MAXOFFSETCOST ;

  if (
      (va->offSetCost > 0 && vb->offSetCost < 0) ||
      (va->offSetCost < 0 && vb->offSetCost > 0)
      )
    doff = 0 ;
  if (va->offSetCost > 0 || vb->offSetCost > 0) /* distance in pair 3' 5' */
    return - (dv - doff) ;
  
  return - (dv - doff) ; /* favor small position on cosmid (marked as negative) 
			 * i changed the sign to dv - doff 2008_09_15, to please unc-116 
			 */
} /* cPathValueDecreassingOrder */

/***********************************************************************/
/***********************************************************************/

static long cPathEvaluateOne (PP pp, int p, int* offSetCostp)
{
  int ii, iiOld ;
  long value = 0 ;
  int nh = pp->nh ;
  BitSet b = pp->bb ;

  for (ii = 0, iiOld = -1 ; ii < nh ; ii++)
    {
      if (bit(b, p * nh + ii))
	{
	  value += pp->oneValue (pp, iiOld, ii, offSetCostp) ;
	  iiOld = ii ;
	}
    }  
  return value ;
} /* cPathEvaluateOne */

/***********************************************************************/

static void cPathEvaluateAll (PP pp)
{
  int p = pp->pMax ;
  VALUE *vp ;

  pp->values = arrayReCreate (pp->values, pp->pMax, VALUE) ;
  while (p--)
    {
      vp = arrayp (pp->values, p, VALUE) ;
      vp->p = p ;
      vp->value = cPathEvaluateOne (pp, p, &(vp->offSetCost)) ;
    }
  arraySort (pp->values, cPathValueDecreassingOrder) ;
} /* cPathEvaluateAll */

/***********************************************************************/
/************* algorithmic utilities ***********************************/

static int cPathCount (PP pp, int ii)
{
  int p = pp->pMax ;
  int nh = pp->nh ;
  BitSet b = pp->bb ;
  int n = 0 ;
  while (p--)
    if (bit(b, p * nh + ii))
      n++ ;

  return n ; 
} /* cPathCount */

/***********************************************************************/
static BOOL cPathIsCognate (PP pp, int ii, int jj)
{
  int p ;
  int nh = pp->nh ;
  BitSet b = pp->bb ;
  BitSet d = pp->dead ;

  p = pp->pMax ;
  while (p--) /* jj is allready in my descendants */
    if (bit(b, p * nh + ii) && bit(b, p * nh + jj))
	return TRUE ;
  p = pp->pMax ;
  while (p--)  /* i know i hate jj */
    if (bit(d, p * nh + ii) && bit(b, p * nh + jj))
	return TRUE ;

  return FALSE ;
} /* cPathIsCognate */

/***********************************************************************/

static BOOL cPathIsAntiCognate (PP pp, int jj)
{
  int p = pp->pMax ;
  int nh = pp->nh ;
  BitSet b = pp->bb ;
  BitSet no = pp->no ;

  while (p--)
    if (bit(no, p) && bit(b, p * nh + jj))
	return TRUE ;

  return FALSE ;
} /* cPathIsAntiCognate */

/***********************************************************************/

static void cPathAddNonCompatible (PP pp, int ii)
{
  int p = pp->pMax ;
  int nh = pp->nh ;
  BitSet b = pp->bb ;
  BitSet no = pp->no ;

   while (p--)
     if (bit(b, p * nh + ii))
       bitSet (no, p) ;
} /* cPathAddNonCompatible */

/***********************************************************************/
/* create a new path with ii as only member */
static void cPathNewPath (PP pp, int ii)
{
  int p = pp->pMax ;
  int nh = pp->nh ;
  BitSet b = pp->bb ;
 
  bitSet (b, p * nh + ii) ; 
  pp->pMax++ ;  
  
  bitUnSet (pp->bb, 2 * pp->nh * pp->pMax) ; /* make room */
  bitUnSet (pp->dead, 2 * pp->nh * pp->pMax) ; /* make room */
}

/***********************************************************************/
/* ii is compatible with jj and inherits all its paths */
static void cPathInheritPath (PP pp, int ii, int jj)
{
  int p ;
  int nh = pp->nh ;
  BitSet b = pp->bb ;
  BitSet d = pp->dead ;

  p = pp->pMax ;
  while (p--)  /* inherit live ones */
    if (bit(b, p * nh + jj))
      bitSet (b, p * nh + ii) ;
  
  p = pp->pMax ;
  while (p--)  /* inherit dead ones */
    if (bit(d, p * nh + jj))
      bitSet (d, p * nh + ii) ;
} /* cPathInheritPath */

/***********************************************************************/

static void cPathAddPathDownstream (PP pp, int startHit, int pOld, int pNew)
{
  int jj ;
  int nh = pp->nh ;
  BitSet b = pp->bb ;

  for (jj = startHit ; jj < pp->nh ; jj++)
    if (bit (b, pOld * nh + jj))
	bitSet (b, pNew * nh + jj) ;
} /* cPathAddPathDownstream */

/***********************************************************************/

static BOOL cPathAddNeeded (PP pp, int ii, int startHit, int pOld)
{
  int jj ;
  int nh = pp->nh ;
  BitSet b = pp->bb ;

  for (jj = startHit ; jj < pp->nh ; jj++)
    if (bit (b, pOld * nh + jj) &&
	! cPathIsCognate (pp, ii, jj))
      return TRUE ;
  return FALSE ;
} /* cPathAddNeeded */

/***********************************************************************/

static void cPathDuplicate (PP pp, int ii, int jj)
{
  int p, oldMax = pp->pMax, nn;
  int nh = pp->nh ;
  BitSet b = pp->bb ;

  for (p = 0, nn = 0 ; p < oldMax ; p++)
    {
      if (bit(b, p * nh + jj) &&
	  cPathAddNeeded (pp, ii, jj, p)) 
	{
	  nn++ ;
	  cPathAddPathDownstream (pp, jj, p, pp->pMax) ;
	  cPathNewPath (pp, ii) ; /* increments pp->pMax */
	}
    }
  /* because of bestSubPath selection addition should occur just once */
  if(nn > 100)
    fprintf (stderr, "anomaly in cPathDuplicate %s nn = %d\n", pp->title, nn)  ;
} /* cPathDuplicate */

/***********************************************************************/

static void cPathSelectBestSubPath (PP pp, int ii, int oldpMax)
{
  int p, pMax = pp->pMax ;
  int nh = pp->nh ;
  BitSet b = pp->bb ; 
  BitSet d = pp->dead ; 
  int pBest = -1 ;
  long value, bestValue = -1 ;
  int offSetCost = 0, bestOffSetCost = 0 ;

  for (p = 0 ; p < pMax ; p++)  
    {
      if (bit(b, p * nh + ii))
	{
	  offSetCost = 0 ;
	  value = cPathEvaluateOne (pp, p, &offSetCost) ;
	  if (pBest == -1)
	    { pBest = p ; bestValue = value ; bestOffSetCost = offSetCost ; }
	  else
	    {  
	      long dv = bestValue - value ;
	      int doff = bestOffSetCost - offSetCost ;
	      
	      if (doff > MAXOFFSETCOST)
		doff = MAXOFFSETCOST ;
	      if (doff < - MAXOFFSETCOST)
		doff = - MAXOFFSETCOST ;
	      if (
		  (bestOffSetCost > 0 && offSetCost < 0) ||
		  (bestOffSetCost < 0 && offSetCost > 0)
		  )
		doff = 0 ;
	      if (
		  (bestOffSetCost > 0 && dv - doff < 0) || /* distance in pair 3' 5' */
		  (bestOffSetCost < 0 && dv + doff < 0)  /* favor small offset (marked as negative) */
		  )
		{ pBest = p ; bestValue = value ; bestOffSetCost = offSetCost ; }
	    }
	}
    }

  for (p = 0 ; p < pMax ; p++)  
    {
      if (p != pBest && bit(b, p * nh + ii))
	{
	  int jj ;
	  for (jj = ii ; jj < nh ; jj++)
	    {
	      bitUnSet (b, p * nh + jj) ;
	      bitUnSet (d, p * nh + jj) ;
	      if (p < oldpMax)
		break ; /* mieg 2007_01_11 do not kill pre-exiting lower paths */
	    }
	}
    }	
} /* cPathSelectBestSubPath */
  
/***********************************************************************/
/***********************************************************************/
/* algorithm real work is here */
static void cPathMakeOne (PP pp, int ii)
{
  int jj, nNew = 0 ;
  int oldpMax = pp->pMax ;

  bitUnSet (pp->bb, 2 * pp->nh * (pp->pMax + 1)) ; /* make room */
  bitUnSet (pp->dead, 2 * pp->nh * (pp->pMax + 1)) ; /* make room */

  pp->no = bitSetReCreate (pp->no, 2 * (pp->pMax + 1)) ;
  for (jj = ii+1 ; jj < pp->nh ; jj++)
    {
      if (cPathIsCognate (pp, ii, jj))   /* pair ii jj was already evaluated */
	continue ;
      switch (pp->isCompatible (pp, ii, jj))
	{
	case 0:     /* start a new path */ 
	  cPathAddNonCompatible (pp, jj) ; /* signal that future path hates all paths with jj */
	  continue ;
	case 1:     /* huge intron, accept extension and at the same time split */
          if (!nNew++) cPathNewPath (pp, ii) ; /* new */
	  cPathDuplicate (pp, ii, jj) ;        /* odl and union */
	  break ;
	case 2:     /* extend previous path */
	  if (cPathIsAntiCognate (pp, jj))
	    cPathDuplicate (pp, ii, jj) ;   /* old and union */
	  else
	    cPathInheritPath (pp, ii, jj) ; /* just union: inherit in ii the path of jj */
	  break ;
	} 
    }
  switch (cPathCount (pp, ii))
    {
    case 0:  cPathNewPath (pp, ii) ; break ;  /* new singlet */
    case 1:  break ;
    default: cPathSelectBestSubPath (pp, ii, oldpMax) ; break ;
    }
} /* cPathMakeOne */

/***********************************************************************/
/******************  public interface  *********************************/
/***********************************************************************/
/* this function is for debugging */
static void cPathShow (PP pp)
{
  int ii = 0 ;
  int p, pMax = pp->pMax ;
  int nh = pp->nh ;
  BitSet b = pp->bb ;
  HIT *up ;
  VALUE *vv ;

  if (!pp || pp->magic != &myMagic)
    return ;
  printf ("\n a1=%d a2=%d searchRepeats=%d\n", pp->a1, pp->a2,  pp->searchRepeats) ;
  for (ii = 0 ; ii < pp->nh ; ii++)
    {
      up = arrp (pp->hits, ii, HIT) ;
      printf ("%2d: %s %6d %6d   %5d %5d  ", ii, name(up->est), up->a1, up->a2, up->x1, up->x2) ;
      for (p = 0 ; p < pMax ; p++)
	if (bit (b, p * nh + ii)) 
	  printf (" %d", p) ;
      printf ("\n") ;
    }

  if (pp->values)
    for (ii = 0 ; ii < pp->pMax; ii++)
      {
	vv = arrp (pp->values, ii, VALUE) ;
	printf ("path %2d :: v = %ld  offset = %d\n", vv->p, vv->value, vv->offSetCost) ;
      }
} /* cPathShow */

/*********************************************************************/

static int getPleaseAcceptNegativeIntron (void)
{
  static int style = -1 ;

  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone ; AcceptNegativeIntron") ;
      style = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    style = 0 ; /* default */

  return style ;
} /* getPleaseAcceptNegativeIntron */

/************* biological part *****************************************/
/* 0: non compatible
 * 1: huge intron, accept extension and at the same time split
 * 2: extend previous path
 */
#define BIGINTRON 1000000

static int cPathIsCompatible (PP pp, int ii, int jj)
{
  HIT *up = arrp (pp->hits, ii, HIT) ;
  HIT *vp = arrp (pp->hits, jj, HIT) ;
  BOOL isUp = up->x1 > up->x2 ? TRUE : FALSE ; /* up gene */
  float cover = 1.5 ; /* tried 1.5 then 1.0, 1.0 is bad in nm_020364.1, 2 is bad in 4f413 */
  int da, dx, dx1, dx2 ;
  int negativeGap = getPleaseAcceptNegativeIntron ()  ;
  da = vp->a1 - up->a2 ; /* a gap == intron */
  dx = isUp ? up->x2 - vp->x1 - 1 : vp->x1 - up->x2 - 1 ; /* positif : gap, negatif : intron overshooting */
  dx1 = up->x2 - up->x1 ;
  dx2 = vp->x2 - vp->x1 ;
  if (dx1 < 0) dx1 = - dx1 ;
  if (dx2 < 0) dx2 = - dx2 ;

#if 0
  if (da + dx + 40 < 0  ||  /* x overlap larger than a gap */ /* apr 23 changed 30 to 40 */
      (dx == 0 && da < 5) ||
      (dx != 0 && da - dx < 30) )    /* x gap  larger than a gap, except in exact case  */
    return 0 ;   /* non compatible exons */
#else
  /* changed nov 25 2003 i think signs were all wrong
     in perfect case dx = 0, da - delta = intron size
     if we have a gap of g bp, we get   dx = g > 0,  da = delta + g
     if we overshhoot the intron by n we have dx = -n < 0, da = delta - n
     so to evaluate incompatibility, we always wish to use da - dx
  */
  if (isFinalPath && negativeGap)
    {
      /* micro intron but accept 12 bp deletions in est  */
      if (0 && dx == 0 && da < 5 && da > -5 && dx1 < 20 && dx2 > 20)
	return 0 ;   /* non compatible exons */
      if (da  - dx  < -200) /* x overlap larger than a reasonable genome gap */
	return 0 ;   /* non compatible exons */
    }
  else
    {
      if (dx == 0 && da < 5)/* micro intron but accept 12 bp deletions in est  */
	return 0 ;   /* non compatible exons */
      if (dx && da  - dx  < -30) /* x overlap larger than a gap */
	return 0 ;   /* non compatible exons */
    }
#endif
  if ((!isUp && ( -cover * dx > up->x2 - up->x1 ||   /* was 1.5 was 2.0 */
		  -cover * dx > vp->x2 - vp->x1)) ||
      ( isUp && (-cover * dx > up->x1 - up->x2 || 
		 -cover * dx > vp->x1 - vp->x2)))
    return 0 ; /* intersect >= all was 2/3 (was half until aug 17) the length */

  if (da > BIGINTRON)
    return 1 ;
  return 2 ;
} /* cPathIsCompatible */

/***********************************************************************/

static long cPathOneValue (PP pp, int iiOld, int ii, int *offSetCostp)
{
  long vv = 0 ;
  int da = 0, di = 0, gap = 0, ne = 0 ;
         /* gain = nn *ecost - ng * gapcost - ni *icost */
  int exonCost = 10000, errorCost = 30000, gapCost = 3000, intronCost = 10 ;
  int offSetCost = 0 ;
  int crossOverCost = 0 ;
  BOOL isDown ; /* up or down gene */
  HIT *up, *vp ;

  MAXOFFSETCOST = 6 * errorCost ; /* 6 error */
  OFFSETCOST = 10 ; /* pp->searchRepeats ? 50 : 10 ; */
  if (pp->searchRepeats == 2)
    {
      gapCost *= 100 ;
      intronCost *= 100 ;
      OFFSETCOST  *= 10 ;
    }

  if (ii < 0 || ii >= arrayMax(pp->hits))
    messcrash ("Bad ii = %d in cPathOneValue", ii) ;
  vp = arrp (pp->hits, ii, HIT) ; 
  if (iiOld >= 0)
    {
      if (iiOld >= arrayMax(pp->hits))
	messcrash ("Bad iiOld = %d in cPathOneValue", iiOld) ;
      up = arrp (pp->hits, iiOld, HIT) ; 
    }
  else
    up = 0 ;
  
  isDown = vp->x1 < vp->x2 ? TRUE : FALSE ; 
  /* in case of repeated gene, try to squeeze both half of the gene to the same place */

  if (!up) /* vp is first exon, measure the offsetcost and the exon size */
    {
      /* exon length and errors */
      ne = vp->nerr ;
      if (isDown)
	da = vp->x2 - vp->x1 ;
      else
	da = vp->x1 - vp->x2 ;

      /* offset */
      if (vp->x1 > vp->x2 && pp->a1) /* we are in an up facing 5' read */
	{
	  /* the 200 is because the (unclipped) 3' read may start a bit below the top of the 5' read */
	  if (vp->a1 > pp->a1 - 200)
	    offSetCost = OFFSETCOST * (vp->a1 - pp->a1 + 200) ;  /* malus, squeeze close to pp->a1 */
	  else /* i am on the wrong side of my 3' read */
	    offSetCost = MAXOFFSETCOST ;
	}
      else if (vp->x1 < vp->x2 && pp->a2) /* we are in an down facing 5' read */
	{
	  /* the 200 is because the (unclipped) 3' read may start a bit below the top of the 5' read */
	  if (vp->a1 < pp->a2 + 200)
	    offSetCost = OFFSETCOST * (pp->a2 + 200 - vp->a1) ;  /* malus, squeeze close to pp->a2 */
	  else /* i am on the wrong side of my 3' read */
	    offSetCost = MAXOFFSETCOST ;
	}
      else /* not in a pair, add weigth towards the end of the cosmid */
	offSetCost = - vp->a1 ; /* negative: differ to relative cost between 2 paths */

      /* cross-over cost, abandonned */
      if ( 0 && /* pas bon, cela tue tout le path et pas juste le mauvais exon */
	  ((pp->a1 && vp->x1 > vp->x2 && vp->a1 < pp->a1 - 200) ||
	  (pp->a1 && vp->x1 < vp->x2 && vp->a1 > pp->a1 + 200)
	   ))
	crossOverCost = 150 * errorCost ;
    }
  else  /* up and vp, measure introns overlaps and and gaps */
    {  
      ne = vp->nerr ;
      di =  vp->a1 > up->a2 ? vp->a1 - up->a2 : 0 ; /* intron */
      if (isDown)
	{ 
	        /* da = exons: add x coverage but beware of overlaps */
	  da = vp->x2 - vp->x1 ;
	  if (vp->x1 < up->x2)   /* mrna overlap */
	    da -= up->x2 - vp->x1 ;
	  else                   /* mrna gap */
	    gap = vp->x1 - up->x2 ;
	}
      else
	{
	  da = vp->x1 - vp->x2 ;
	  if (vp->x1 > up->x2)   /* mrna overlap */
	    da -= vp->x1 - up->x2 ;
	  else                   /* mrna gap */
	    gap = up->x2 - vp->x1 ;
	}
    }     
  if (!isFinalPath) /* before fix introns */
    {
      ne -= 8 ; /* we expect about 8 errors due to intron leaking */
      if (vp->x1 < vp->x2)
	{
	  if (vp->x1 == vp->clipTop) ne += 4 ;
	  if (vp->x2 == vp->clipEnd) ne += 4 ;
	}
      if (vp->x1 > vp->x2)
	{
	  if (vp->x2 == vp->clipTop) ne += 4 ;
	  if (vp->x1 == vp->clipEnd) ne += 4 ;
	}
      if (ne < 0) ne = 0 ;
    }
  if (di > BIGINTRON) di = BIGINTRON ;
  vv = (long)da * exonCost - ne * errorCost - gap * gapCost - di * intronCost - crossOverCost ;
  if (offSetCost) *offSetCostp = offSetCost ;
  return vv ;
} /* cPathOneValue */

/***********************************************************************/
/* remove hits present 20 times */
static BOOL  cPathSimplifyPath (Array hits)
{
  HIT *up, *vp, *wp ; 
  int n, j1, j2, jMax = arrayMax(hits) ;
  int dx, dx0, dx00 ; /* useful for genomic clones */
  
  /* look for silly repeats */
  if (jMax > 1)
    for (j1 = jMax - 1, up = arrp(hits, j1, HIT) ; j1 >= 0 ; j1--, up--) 
      {
	n = 0 ; /* number of inclusions */
	if (! up->est)
	  continue ;
	wp = up ;
	dx0 = dx00 = up->x2 - up->x1 ;
	if (dx00 < 100) dx00 = 0 ;
	for (j2 = j1 - 1, vp = up - 1 ;  n < 10 && j2 > 0 ; j2--, vp--)
	  {
	    if (! vp->est)
	      continue ;
	    if (vp->x1 < up->x1 - 5)
	      break ;
	    if (vp->x2 > up->x2 + 5) 
	      continue ;
	    dx = vp->x2 - vp->x1 ;
	    if (dx00 > 20 * dx)
	      vp->est = 0 ;
	    if (vp->est) 
	      {
		n++ ;
		if (vp->nerr + dx0 > wp->nerr + dx) wp = vp ;
	      }
	  }
	for (j2 = j1 + 1, vp = up + 1 ; n < 10 && j2 < jMax ; j2++, vp++) 
	  {
	    if (! vp->est)
	      continue ;
	    if (vp->x1 > up->x2)
	      break ;
	    /* vp->x1 >= up->x1 is garanteed */
	    if (vp->x2 > up->x2 + 5)  
	      continue ; 
	    dx = vp->x2 - vp->x1 ;
	    if (dx00 > 20 * dx)
	      vp->est = 0 ;
	    if (vp->est) 
	      {
		n++ ;
		if (vp->nerr + dx0 > wp->nerr + dx) wp = vp ;
	      }
	  }
	if (n >= 6)  /* discard overwhelming repeats except first 6 */
	  wp->est = 0 ;
      }
  
  /* register happy few */
  j2 = 0 ;
  if (jMax)
    {
      for (j1 = j2 = 0, up = arrp(hits, j1, HIT) ; j1 < jMax ; j1++, up++) 
	{ 
	  if (up->est)
	    {
	      if (j1 != j2)
		arr (hits, j2, HIT) = arr(hits, j1, HIT) ;  
	      j2++ ;
	    }
	}
      arrayMax (hits) = j2 ;
    }
  return j2 < jMax ? TRUE : FALSE ; /* di i simplify ? */
} /* cPathSimplifyPath */

/***********************************************************************/
/*
 * a1>0 : indicates the top base of the down looking 3' read of my current 5' 
 * a2>0 : indicates the final base of the up looking 3' read of my current 5' 
 */
static void cPathSelectOneBestPath (Array hits, Array estHits, int a1, int a2, int searchRepeats)
{
  PP pp = 0 ;
  char *title = "toto" ;
  
  if (arrayMax (estHits) > 1)
    {
      cDNASwapX(estHits) ; 
      arraySort (estHits, cDNAOrderByX1) ; 
      
      while (cPathSimplifyPath (estHits)) ;
      
      cDNASwapA(estHits) ; 
      arraySort (estHits, cDNAOrderByA1) ; 
    }

  if (arrayMax (estHits))
    {
      if (hits && arrayMax(hits))
	title = name(arr (hits, 0, HIT).est) ; /* used for debugging */
      pp = cPathCreate (title, estHits, cPathOneValue, cPathIsCompatible) ;
      pp->a1 = a1 ; pp->a2 = a2 ; pp->searchRepeats = searchRepeats ;
      cPathMakeAll (pp) ;
      cPathEvaluateAll (pp) ; /* sorts */
      if (0) cPathShow (pp) ; /* for debugging */

      if (pp->values && arrayMax(pp->values))
	{
	  int p = arrp(pp->values, 0, VALUE)->p ;
	  int ii, jj ;
	  BitSet b = pp->bb ;
	  int nh = pp->nh ;
	  
	  for (ii = 0, jj = arrayMax(hits) ; ii < arrayMax(estHits) ; ii++)
	    {
	      if (bit (b, p * nh + ii)) /* register happy few */
		array (hits, jj++, HIT) = arr (estHits, ii, HIT) ;
	    }
	}
      cPathDestroy (pp) ;
    }
} /* cPathSelectOneBestPath */

/*********************************************************************/

static int cDNAOrderByA1Path (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va, *vp = (const HIT *)vb ;

  if (up->cDNA_clone != vp->cDNA_clone)
    return up->cDNA_clone - vp->cDNA_clone ;  /* the cDNA_clone */
  else if (up->x1 < up->x2 && vp->x1 > vp->x2)
    return -1 ;
  else if (up->x1 > up->x2 && vp->x1 < vp->x2)
    return 1 ;  /* down read first */
  else if (up->est != vp->est)
    return up->est - vp->est ;  /* the ESt */
  else if (up->a1 != vp->a1)
    return up->a1 - vp->a1 ;
  else if (up->a2 != vp->a2)
    return up->a2 - vp->a2 ;      
  else  
    return up->x1 - vp->x1 ;    
} /* cDNAOrderByA1Path */

/***********************************************************************/
/* select for each est the best path */
void cPathSelectBestPath (Array hits, BOOL final, int searchRepeats)
{
  Array oldHits = 0, estHits = 0 ;
  int a1 = 0, a2 = 0, ii, jj = 0, j1 ;
  KEY est = 0 ;
  HIT *up, *vp ;

  isFinalPath = final ;
  
  if (!hits || !arrayMax(hits)) 
    return ;
  cDNASwapA(hits) ; 
  arraySort (hits, cDNAOrderByA1Path) ; 

  oldHits = arrayCopy (hits) ;
  arrayMax(hits) = 0 ;
  for (ii = 0, jj = 0 ; ii < arrayMax(oldHits) ; ii++)
    {
      up = arrp (oldHits, ii, HIT) ;
      if (!est)
	{
	  est = up->est ; 
	  estHits = arrayReCreate (estHits, arrayMax(oldHits), HIT) ;
	  jj = 0 ;
	}
      if (est != up->est)
	{
	  est = up->est ;
	  /* run the graph theory algo on one est at a time */
	  /* a1 remembers the start of the down-looking EST */
	  /* so that the up looking EST will be selected down of a1 */
	  j1 = arrayMax (hits) ;
	  cPathSelectOneBestPath (hits, estHits, a1, a2, searchRepeats) ;
#ifndef JUNK
	  /* case i want down read first */
	  a1 = a2 = 0 ;
	  if (j1 < arrayMax (hits) && 
	      (vp = arrp (hits, j1, HIT) ) &&
	      up->cDNA_clone == vp->cDNA_clone &&
	      up->x1 > up->x2 && vp->x1 < vp->x2)
	    a1 = vp->a1 ;
#else
	  /* case i want 3' read first */
	  a1 = a2 = 0 ;
	  if (j1 < arrayMax (hits) &&  
	      (vp = arrp (hits, j1, HIT) ) &&
	      vp->reverse && !up->reverse &&
	      up->cDNA_clone == vp->cDNA_clone)
	    {
	      if (vp->x1 < vp->x2)  /* top exon of 3' looking down new alignment */
		  a1 = vp->a1 ;
	      else
		{
		  j1 = arrayMax (hits) - 1 ;
		  vp = arrp (hits, j1, HIT) ; /* last exon of 3' looking down new alignment */
		  a2 = vp->a2 ;
		}
	    }
#endif
	  arrayReCreate (estHits, arrayMax(oldHits), HIT) ;
	  jj = 0 ;
	} 
      /* commented out 2007_02_14 i do not understand these lines
      if ( getPleaseDoRepeats (0) && / oct 2004, stupid since i do not verify if 3' is better than 5' 
	  (
	   (a1 && up->x1 > up->x2 && up->a1 < a1 - 200) ||
	   (a1 && up->x1 < up->x2 && up->a1 > a1 + 200)
	   )
	   ) ;
	   else
      */
	array (estHits, jj++, HIT) = arr (oldHits, ii, HIT) ;
    } 
  if (jj) cPathSelectOneBestPath (hits, estHits, a1, a2, searchRepeats) ;
  arrayDestroy (estHits) ;
  arrayDestroy (oldHits) ;
} /* cPathSelectBestPath */

/***********************************************************************/
/***********************************************************************/
