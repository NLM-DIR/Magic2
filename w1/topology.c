/*  File: topology.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@frmop11.bitnet
 *
 * Description:
 * Exported functions:
    topoConnectedComponents()    
 * HISTORY:
 * Last edited: Mar  6 17:11 1996 (srk)
 * Created: Fri Oct 30 19:08:43 1992 (mieg)
 *-------------------------------------------------------------------
 */

/*  $Id: topology.c,v 1.6 2020/05/25 00:59:21 mieg Exp $ */
 
#include "ac.h"
#include "topology.h"

/********************************/

static int topoFirstKeyOrder (const void*a, const void*b)
{
  return ((const LINK*) a) -> a -  ((const LINK*) b)-> a ;
}

/**************************************/

static void topoDoubleUnorientedLinks(Array links)
{
  register int i = arrayMax(links) ; 
  register LINK *x , *y ;

  arrayExtend(links, 2*i) ;
  arrayMax(links) = 2*i ;
  x = arrp(links, 0, LINK) ;
  y = arrp(links, i, LINK) ;
  while(i--)
    { x -> group = y -> group = 0 ;
      *y = *x ;  /* anything else, i.e. type */
      y -> a = x -> b ; 
      y -> b = x -> a ;
      x++ , y++ ;
    }
  arraySort(links, topoFirstKeyOrder) ;
}

/**************************************/

static void topoEnumerateVertex(Array links, Array vx)
{
  register int i , n = arrayMax(links) , j = 0;
  register KEY k = 0 ;

  arrayMax(vx) = 0 ;
  for(i = 0 ; i < n ; i++)
    if ( k != arrp(links,i,LINK) -> a)
      {  k = arrp(links,i,LINK) -> a ;
	 array(vx, j, VERTEX).a = k ;
	 array(vx, j, VERTEX).group = 0 ;
	 array(vx, j, VERTEX).pos = i ;
	 j++ ;
       }
}

/**************************************/

static Array VX, LK;
static int GP ;      /* number of the connected component */

/******************/

static int topoFindVertex(KEY k)
{
  register int i = arrayMax(VX) ;
  while(i--)
    if(arrp(VX,i,VERTEX)->a == k)
      return i ;
  messcrash("Topofindvertex cannot find %d", k) ;

  return 0 ; /* for compiler happiness */
}

/******************/

static void topoNCC2(int n)
{ KEY a = array(VX,n,VERTEX).a ;
  int pos = array(VX,n,VERTEX).pos ;

  array(VX,n,VERTEX).group = GP ;
  for(;pos < arrayMax(LK)  && array(LK, pos,LINK).a == a;pos++)
    { n = topoFindVertex(array(LK, pos,LINK).b) ;
      if(!arr(VX,n,VERTEX).group)
	  topoNCC2(n) ;
      else if(arr(VX,n,VERTEX).group != GP )
	messcrash
	  ("topoNCC2 : Double numbering, sorry") ;
    }
}

/******************/

static void topoNumberLinks(Array links)
{ KEY a ;
  register int i = arrayMax(links), n ;
  while(i--)
    {  a = arr(links,i,LINK).a ;
       n = topoFindVertex(a) ;
       arr(links,i,LINK).group = arr(VX,n,VERTEX).group ;
     }
}

/*****************/

static int topoNumberConnectedComponents(Array links, Array vx)
{
  register int n = arrayMax(vx) ;

  VX = vx ; LK = links; GP = 0 ;
  while(n--)
    if(! arr(vx,n,VERTEX).group)
      { GP++ ;
	topoNCC2(n) ;
      }
  return GP ;
}

/**************************************************************/

int topoVertexOrder (const void *v1, const  void *v2)
{
  int dx ;
  const VERTEX *w1 = (const VERTEX *)v1, *w2 = (const VERTEX *)v2 ;
  if (w1->a < w2->a) return -1 ;
  if (w1->a > w2->a) return 1 ;
  dx = w1->group - w2->group ;
  if (dx) return dx ;
  dx = w1->pos - w2->pos ;
  if (dx) return dx ;
  return 0 ;
} /* topoVertexOrder */

/**************************************************************/

int topoLinkOrder (const void *v1, const  void *v2)
{
  int dx ;
  const LINK *w1 = (const LINK *)v1, *w2 = (const LINK *)v2 ;
  if (w1->a < w2->a) return -1 ;
  if (w1->a > w2->a) return 1 ;
  if (w1->b < w2->b) return -1 ;
  if (w1->b > w2->b) return 1 ;
  dx = w1->group - w2->group ;
  if (dx) return dx ;
  dx = w1->type - w2->type ;
  if (dx) return dx ;
  return 0 ;
} /* topoLinkOrder */

/**************************************************************/
/************* Public Routines ********************************/

int topoConnectedComponents(Array links, Array vx)
{
  Array cpLinks = arrayCopy(links) ;
  int i ;

  if(!arrayExists(links) || !arrayExists(vx))
    messcrash("topoConnectedComponents called with bad array") ;

  topoDoubleUnorientedLinks(cpLinks) ;
  topoEnumerateVertex(cpLinks, vx) ;
  i = topoNumberConnectedComponents(cpLinks, vx) ;
  arrayDestroy(cpLinks) ;
  topoNumberLinks(links) ;
  return i ;
}

/**************************************************************/
/**************************************************************/
/* chonoSort a 2D matrix */

typedef struct wiStruct { double w ; int i ; } WI ;

static int wiOrder (const void *a, const void *b) 
{
  const WI *up = (const WI *)a, *vp = (const WI *)b ;
  
  if (up->w < vp->w) return 1 ; /* heawy wins */
  if (up->w > vp->w) return -1 ;
  return up->i - vp->i ;
} /* wiOrder */

/**********************************************************************/

static BOOL topoChronoSortLines (Array aa, int xy, int iMax, int jMax, KEYSET lines, Array ww, int pass)
{
  BOOL moved = FALSE ;
  BOOL debug = FALSE ;
  WI *wi ;
  int i, j, k, k1, dx, dy ;
  double *zf = arrp (aa, 0, double) ;
  double z, z1, zMin, zMax , buf[iMax * jMax] ;
  int lnBuf [xy == 0 ? iMax : jMax] ;

  switch (xy)
    {
    case 0: 
      dx = jMax ;
      dy = 1 ;
      break ;
    case 1:
      dx = 1 ;
      dy = jMax ;
      i = iMax ; iMax = jMax ; jMax = i ;
      break ;
    }
      
  z = 0 ;
  arrayMax (ww) = iMax ;
  wi = arrp (ww, 0, WI) ;
  for (i = 0 ; i < iMax ; i++)
    {
      wi[i].i = -1 ;  
      wi[i].w = 0 ;
    }

  zMin = 30 ; zMax = 100 ;
  if (pass < 0)
    {zMin = 50 ; zMax = 90 ; }
  if (pass < 0)
    { zMin = 50 ; zMax = 70 ; }
  if (pass < 0)
    { zMin = 50 ; zMax = 70 ; }
  if (pass < 0)
    { zMin = 50 ; zMax = 80 ; }

  for (i = k = 0 ; i < iMax ; i++, k += dx)
    {
      int j0 = 0 ;
      z = 0 ;
      for (j = k1 = 0 ; j < jMax ; j++, k1 += dy)
	{
	  z1 = zf[k+k1] ;
	  if (z1 > zMax) z1 = zMax + (z1-zMax)/10 ;
	  if (z1 >= zMin)
	    z1 = 1 + (pass > 100 ? z1/30000 : 0) ;
	  else
	    z1 = 0 ; 
	  if (z1 && ! z)
	    { j0 = j + 10 ; z += 10*(exp (1.0*log(jMax - j))) * z1 ;  }
	  else
	    z += (j0 - j) * z1 ; 
	} 
      wi = arrp (ww, i, WI) ;
      wi->w = z ;
      wi->i = i ;
    }
  if (debug)
    {
      int kk = 0 ;
      fprintf (stderr, "#A pass %d xy=%d ", pass, xy) ; 
      for (i = k = 0 ; i < iMax ; i++, k += dx)
	{
	  wi = arrp (ww, i, WI) ;
	  fprintf (stderr, " %d", (int) wi->w) ;
	  if (wi->w > 1) kk++ ;
	}
      fprintf (stderr, " kk=%d\n" ,kk) ;
    }

  arraySort (ww, wiOrder) ;
  if (debug)
    {
      int kk = 0 ;
      fprintf (stderr, "#B pass %d xy=%d ", pass, xy) ; 
      for (i = k = 0 ; i < iMax ; i++, k += dx)
	{
	  wi = arrp (ww, i, WI) ;
	  fprintf (stderr, " %d", (int) wi->w) ;
	  if (wi->w > 1) kk++ ;
	}
      fprintf (stderr, " kk=%d\n" ,kk) ;

    }
  wi = arrp (ww, 0, WI) ;
  for (i = 0 ; i < iMax ; wi++, i++)
    {
      k = wi->i ;  
      if (k > iMax - 1)
	messcrash ("bad value") ;
    }

#ifdef JUNK
  /* allow this block to give the system a chance to
   * escape false minima
   * 2020_03_11
   * works but does not seem really useful on the example of crn who_is_who
   */ 
  if (0)
    {
      wi = arrp (ww, 0, WI) ;
      if (iMax > 1 && (pass % 3) == 0)
	{
	  for (i = 0 ; i < iMax - 1 ; wi++, i++)
	    {
	      int k1 = wi->i ;
	      
	      if (k1 >=0)
		{
		  extern double randdouble (void) ;
		  double z = randdouble () ;
		  if (z < .1)
		    {
		      int j = iMax * 5 * z ;
		      if (j < iMax && j != i)
			{
			  int k2 = k1 + 1 ;
			  if (k2 > 0)
			    {
			      moved = TRUE ;
			      if (k1 > iMax - 1)
				messcrash ("bad value") ;
			      if (k2 > iMax - 1)
				messcrash ("bad value") ;
			      wi->i= k2 ;
			      arrp (ww, j, WI)->i = k1 ;
			      i++ ;
			    }
			}
		    }
		}
	    }
	  
	  if (debug)
	    {
	      int kk = 0 ;
	      fprintf (stderr, "#C pass %d xy=%d ", pass, xy) ; 
	      for (i = k = 0 ; i < iMax ; i++, k += dx)
		{
		  wi = arrp (ww, i, WI) ;
		  fprintf (stderr, " %d", (int) wi->w) ;
		  if (wi->w > 1) kk++ ;
		}
	      fprintf (stderr, " kk=%d\n" ,kk) ;
	    }
	}
    }
#endif

  wi = arrp (ww, 0, WI) ;
  for (i = 0 ; i < iMax ; wi++, i++)
    {
      k = wi->i ;  
      if (k > iMax - 1)
	messcrash ("bad value") ;
    }

  memcpy (buf, zf, iMax * jMax * sizeof(double)) ;
  memcpy (lnBuf, arrp (lines, 0, KEY), iMax * sizeof(KEY)) ;

 if (debug)
    {
      for (i = k = 0 ; i < iMax ; i++,wi++, k += dx)
	{
	  for (j = k1 = 0 ; j < jMax ; j++, k1 += dy)
	    fprintf (stderr, "\t%d", (int)zf[k+k1]) ;
	  fprintf (stderr, "\n") ;
	}
    }	     

  wi = arrp (ww, 0, WI) ;
  for (i = k = 0 ; i < iMax ; i++,wi++, k += dx)
    {
      int i2 = wi->i ;
      if (i2 >= 0)
	{
	  keySet (lines, i) = lnBuf[i2] ; 
	  if (i2 != i)
	    {
	      int k2 = i2 * dx ;
	      moved = TRUE ;
	      for (j = k1 = 0 ; j < jMax ; j++, k1 += dy)
		zf[k+k1] = buf[k2+k1] ;
	    }
	}
    }
  
  if (debug)
    {
      for (i = k = 0 ; i < iMax ; i++,wi++, k += dx)
	{
	  for (j = k1 = 0 ; j < jMax ; j++, k1 += dy)
	    fprintf (stderr, "\t%d", (int)zf[k+k1]) ;
	  fprintf (stderr, "\n") ;
	}
    }	     

  return moved ;
} /* topoChronoSortLines */

/**********************************************************************/
/* sort a 2D matrix recursivelly by the total weights of its lines and columns
 * the order is retured as 'lines' and 'cols'
 */
BOOL topoChronoOrder (Array aa, KEYSET lines, KEYSET cols)
{
  AC_HANDLE h = ac_new_handle() ;
  int i, j, iMax, jMax, pass ;
  BOOL moved = TRUE ;
  Array ww = 0 ;
  BOOL isSym = FALSE ;

  if (! aa)
    messcrash ("topoChronoOrder received a null array") ;
  if (! lines)
    messcrash ("topoChronoOrder received a null lines keyset") ;
  if (! cols)
    messcrash ("topoChronoOrder received a null cols keyset") ;
  if  (aa->size != sizeof(double))
    messcrash ("topoChronoOrder did not receive an array of double: size=%d != sizeof(double) = %d", aa->size, sizeof(double)) ;
  if (arrayMax (aa) != keySetMax (lines) * keySetMax(cols))
    messcrash ("topoChronoOrder received non matching sizes (max(aa)=%d) != (max(lines=%d) * max(cols)=%d)", arrayMax(aa), keySetMax (cols), keySetMax(lines)) ;


  i = iMax = keySetMax (lines) ;
  j = jMax = keySetMax (cols) ;

  for (i = 0 ; i < iMax ; i++)
    {
      int run1 = keySet (lines, i) ;
      for (j = 0 ; j < jMax ; j++)
	{
	  int run2 = keySet (cols, j) ;
	  if (run1 != run2)
	    isSym = FALSE ;
	}
    }

  /* original order */
  while (i--)  keySet (lines, i) = i ;
  while (j--)  keySet (cols, j) = j ;
  /* after sorting lines[2]=7 means that lines 2 was called 7 in lines0
   * so the correct metadata of line i =  title[lines[i]]
   */
  i = iMax >  jMax ? iMax : jMax ;
  ww = arrayHandleCreate (i, WI, h) ;
  array (ww, i - 1, WI).w = 0 ; /* make room */
    
  for (pass = 0 ; moved && pass < 3000 ; pass++)
    {
      moved = topoChronoSortLines (aa, 0, iMax, jMax, lines, ww, pass) ;
      if (0 && isSym) /* FAUX, il faut remanier la matrices, pas juste la table des cosl */
	for (j = 0 ; j < jMax ; j++)
	  keySet (cols, j) = keySet (lines, j) ;
      else
	moved |= topoChronoSortLines (aa, 1, iMax, jMax, cols, ww, pass) ;
 
      if (! moved && i < 3) { moved = TRUE ; }
      if (! moved && i < 100) { i = 100 ; moved = TRUE ; }
    }
  fprintf (stderr, "topo stop at pass=%d\n", pass) ;
  ac_free (h) ;
  return TRUE ;
} /* topoChronoOrder */

/**********************************************************************/
/**********************************************************************/
