 /*  File: polynome.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg 2024
 *-------------------------------------------------------------------
 * This file is part of a toolbox for formal calculus in QFT
 *
 * Description:
 *    Handle Polynomes: 
 *                 show, add, mutiply, exponentiate
 *    Implements t'Hooft dimensional regularized integrals
 *                 by defining momenta and denominators 
 *    The elments of the polynomes are
 *          complex numbers, 
 *          rational fractions (to avoid rounding errors)
 *          symbols (x,y,...), 
 *          Grassman variables (i,j,k...)
 *          Metric with dummy indices and epsilon symbols
 *          Pauli matrices sigma and sigmaBar SIG. SIGB
 *          Momenta with a vector index, 
 *          Specialized propagator denominators 1/k^2(k+p)^2(k+p+q)^2
 *    Handle matrices of polynomes PMX
 *                 show, add, mutiply, exponentiate, determinant.
 * Exported functions:
 *              See Header file: polynome.h
 * HISTORY:
 * Last edited:
 * Created: Feb 2024 (mieg)
 *-------------------------------------------------------------------
 */

#include "polynome.h"
static int POLMAGIC = 153724322 ;
static int TTMAGIC = 0 ;                    /* zero means we did not overwrite, but no need to set it */
static int PMXMAGIC = 871537243 ;

/***********************************************************************************************************************************************/
/***************************************************************** security     ******************************************************************/
/***********************************************************************************************************************************************/

static void polCheck (POLYNOME pp)
{
  if (pp)
    {
      if  (pp->magic != POLMAGIC)
	messcrash ("Bad POLMAGIC") ;
      if (pp->tt.magic != TTMAGIC)
	messcrash ("Bad TTMAGIC") ;  
    }
} /* polCheck */

/***********************************************************************************************************************************************/
static void pmxCheck (PMX px)
{
  if (px)
    {
      int iMax = px->N * px->N ;

      if  (px->magic != PMXMAGIC)
	messcrash ("Bad PMXMAGIC") ;
      if (px->N <= 0)
	messcrash ("Bad pmx N=%d in %s", px->N, px->title) ;
      for (int i = 0 ; i < iMax ; i++)
	polCheck (px->pp[i]) ;
    }
} /* pmxCheck */

/***********************************************************************************************************************************************/
/***************************************************************** Utilities  ******************************************************************/
/***********************************************************************************************************************************************/

short firstDummyIndex = 'a' ;

short newDummyIndex (void)
{
  short cc = firstDummyIndex++ ;
  
  if (cc >= INDEXMAX)
    messcrash ("too many dummy indices") ;
  return (short)cc  ;
}


/* Warning : we have no provision to store a big array in the database */
 /***********************************************************************************************************************************************/
/* look for all indices, 
 * edit = FALSE: count occurences of a,b,c...
 * edit = TRUE: replace by a continuous list liberating the high values
 */
static void reduceIndicesTtDo (short *cp, KEYSET ks, BOOL edit)
{  
  int i ;
  
  for (i = 0 ; *cp && i<GMAX ; cp++, i++)
    {
      if (! edit)
	keySet (ks, *(short *)cp) ++ ;
      else
	*cp = keySet (ks, *(short *)cp) ;
    }
  return ;
} /* reduceIndices */

/***********************************************************************************************************************************************/

static void reduceIndicesTt (POLYNOME pp, KEYSET ks, BOOL edit)
{
  polCheck (pp) ;
  reduceIndicesTtDo (pp->tt.g, ks, edit) ;
  reduceIndicesTtDo (pp->tt.gg, ks, edit) ;
  reduceIndicesTtDo (pp->tt.sigma, ks, edit) ;
  reduceIndicesTtDo (pp->tt.sigB, ks, edit) ;
  reduceIndicesTtDo (pp->tt.eps, ks, edit) ;
  reduceIndicesTtDo (pp->tt.mm[0], ks, edit) ;
  reduceIndicesTtDo (pp->tt.mm[1], ks, edit) ;
  reduceIndicesTtDo (pp->tt.mm[2], ks, edit) ;
  reduceIndicesTtDo (pp->tt.mm[3], ks, edit) ;
  polCheck (pp) ;
} /* reduceIndicesTt */

/***********************************************************************************************************************************************/

static void reduceIndicesDo (POLYNOME pp, KEYSET ks, BOOL edit)
{
  if (! pp) 
    return ;
  polCheck (pp) ;
  reduceIndicesDo (pp->p1, ks, edit) ;
  reduceIndicesDo (pp->p2, ks, edit) ;
  if (pp->tt.type)
    reduceIndicesTt (pp, ks, edit) ;
} /* reduceIndices */

/***********************************************************************************************************************************************/

static int reduceIndicesIsProduct (POLYNOME pp)
{
  int i, j, n ;
  KEYSET ks = keySetCreate () ;
  KEYSET ks2 = keySetCreate () ;
  KEYSET ks3 = keySetCreate () ;

  polCheck (pp) ;
  /* count the occurence of all indices */
  reduceIndicesDo (pp, ks, FALSE) ;
  /* contract the list by dropping unused values */
  for (i = j = 0 ; 'a' + i < keySetMax (ks) ; i++)
    {
      n = keySet (ks, 'a' + i) ;
      if (n == 1)    /* block and do not rename the free indices present once */
	{
	  keySet (ks2, 'a' + i)  = 'a' + i ;
	  keySet (ks3, 'a' + i)  = 1 ;       /* block it */ 
	  if (firstDummyIndex <= 'a' + i)
	    firstDummyIndex = 'a' + i + 1 ;
	}
      else if (n > 2)
	{
	  showPol (pp) ;
	  messcrash ("no index (%c) should be repeated %d > 2 times", 'a'+i, n) ;
	}
    }
  for (i = 0 ; 'a' + i < keySetMax (ks) ; i++)
    {
      n = keySet (ks, 'a' + i) ;
      if (n == 2)
	{
	  for (j = 0 ; ; j++)
	    if (keySet (ks3, 'a' + j)  == 0)  /* first unblocked index  */
	      {
		keySet (ks2, 'a' + i)  = 'a' + j ;
		keySet (ks3, 'a' + j)  = 1 ;
		if (firstDummyIndex <= 'a' + j)
		  firstDummyIndex = 'a' + j + 1 ;
		break ;
	      }
	}
    }
  /* rename all indices */
  reduceIndicesDo (pp, ks2, TRUE) ;

  keySetDestroy (ks) ; 
  keySetDestroy (ks2) ; 
  keySetDestroy (ks3) ; 
  return 'a' + j ;
} /* reduceIndicesIsProduct */

/***********************************************************************************************************************************************/

POLYNOME reduceIndices (POLYNOME pp)
{
  static int level = 0 ;

  polCheck (pp) ;
  level++ ;
  if (0 && level == 1)
    {
      pp = squareMomentaCleanUp (pp) ;
    }
  if (pp && pp->isSum)
    {
      if (pp->p1) 
	pp->p1 = reduceIndices (pp->p1) ;
      if (pp->p2) 
	pp->p2 = reduceIndices (pp->p2) ;
    }
  if (pp && ! pp->isSum) reduceIndicesIsProduct (pp) ;
  level-- ;
  if (level == 0 && pp)
    {
      if (! pp->isSum)
	reduceIndicesIsProduct (pp) ;
    }
  return pp ;
} /* reduceIndices */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static int polOrder (const void *va, const void *vb)
{
  const POLYNOME p1 = *(POLYNOME *)va ;
  const POLYNOME p2 = *(POLYNOME *)vb ;
  double complex z1 = p1->tt.z ;
  double complex z2 = p2->tt.z ;
  int id1 = p1->tt.Id2 ;
  int id2 = p2->tt.Id2 ;
  int s ;

  polCheck (p1) ;
  polCheck (p2) ;
  p1->tt.z = 0 ;       
  p2->tt.z = 0 ;
  p1->tt.Id2 = 0 ;
  p2->tt.Id2 = 0 ;
  memset (p1->tt.freeIndex, 0, SMAX) ;
  memset (p2->tt.freeIndex, 0, SMAX) ;
  s = memcmp (&(p1->tt), &(p2->tt), sizeof(TT)) ;
  p1->tt.z = z1 ;       
  p2->tt.z = z2 ;
  p1->tt.Id2 = id1 ;
  p2->tt.Id2 = id2 ;       

  return s ;
}

/***********************************************************************************************************************************************/

static void sortPolGetSum (POLYNOME pp, Array aa)
{
  if (! pp)
    return ;
  if (pp->isSum)
    {
      sortPolGetSum (pp->p1, aa) ;
      sortPolGetSum (pp->p2, aa) ;
    }
  else if (pp)
    array (aa, arrayMax (aa), POLYNOME) = pp ;
  return ;
}

/***********************************************************************************************************************************************/

static POLYNOME sortReduceSum (Array aa, AC_HANDLE h)
{
  POLYNOME pp = 0 ;
  int ii, jj ;

  if (!aa)
    return 0 ;

  if (arrayMax (aa) > 1)
    arraySort (aa, polOrder) ;

  for (ii = 0 ; ii + 1 < arrayMax (aa) ; ii++)
    {	/* add identical terms */				
      int jj ;
      POLYNOME q1 = arr (aa, ii, POLYNOME) ; 
      if (! q1 || ! q1->tt.type)
	continue ;
      for (jj = ii + 1 ; jj < arrayMax (aa) ; jj++)
	{
	  POLYNOME q2 = arr (aa, jj, POLYNOME) ; 
	  if (q2)
	    {
	      int n = polOrder (&q1, &q2) ;
	      if (n == 0)
		{
		  q1->tt.z += q2->tt.z ;
		  arr (aa, jj, POLYNOME)  = 0 ;
		}
	      else if (n < 0)
		break ;
	    }
	}
    }
  for (ii = jj = 0 ; ii < arrayMax (aa) ; ii++)
    {  /* elimimate the killed terms */
      POLYNOME qq = arr (aa, ii, POLYNOME) ; 
      if (qq && (qq->isSum || qq->isProduct || qq->tt.z != 0))
	arr (aa, jj++, POLYNOME) = qq ;
    }
  arrayMax (aa) = jj ;
  
  if (arrayMax (aa) == 0)
    pp = 0 ;
  else if (arrayMax (aa) == 1)
    pp = arr (aa, 0, POLYNOME) ;
  else 
    {
      POLYNOME qq ;
      
      pp =  newPolynome (h) ;
      for (qq = pp, ii = 0 ; ii + 1 < arrayMax (aa) ; ii++)
	{
	  qq->isSum = TRUE ;
	  qq->p1 = arr (aa, ii, POLYNOME) ;
	  if (ii < arrayMax (aa) - 2)
	    { qq->p2 = newPolynome (h) ; qq = qq->p2 ; }
	}
      qq->p2 = arr (aa, ii, POLYNOME) ;
    }
  return pp ;
} /* sortReduceSum */
    
/***********************************************************************************************************************************************/

POLYNOME sortPol (POLYNOME pp)
{
  if (! pp)
    return 0 ;
  else if (pp->isSum)
    {
      int ii = 0, jj = 0 ;
      Array aa = arrayCreate (32, POLYNOME) ;
      
      sortPolGetSum (pp, aa) ;
      if (0)
	{
	  for (ii = jj = 0 ; ii < arrayMax (aa) ; ii++)
	    {
	      POLYNOME qq = arr (aa, ii, POLYNOME) ; 
	      qq = sortPol (qq) ;
	      if (qq)
		arr (aa, jj++, POLYNOME) = qq ;
	    }
	  arrayMax (aa) = jj ;
	}
      pp = sortReduceSum (aa, pp->h) ;
      if (0) arrayDestroy (aa) ;
    }
  else if (pp->isProduct)
    {
      pp->p1 = sortPol (pp->p1) ;
      pp->p2 = sortPol (pp->p2) ;
      if (pp->p1 && pp->p2)
	return pp ;
      return 0 ;
    }
  else if (! pp->tt.type ||  cabs (pp->tt.z) < minAbs)
    return 0 ;

  return pp ;
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static BOOL ppIsNumber (POLYNOME pp)
{
  TT tt = pp->tt ;
  BOOL ok = tt.type ;
  if (tt.x[0])  ok = FALSE ;
  if (tt.mm[0][0])  ok = FALSE ;
  if (tt.mm[1][0])  ok = FALSE ;
  if (tt.mm[2][0])  ok = FALSE ;
  if (tt.mm[3][0])  ok = FALSE ;
  if (tt.g[0])  ok = FALSE ;
  if (tt.gg[0])  ok = FALSE ;
  if (tt.sigma[0])  ok = FALSE ;
  if (tt.sigB[0])  ok = FALSE ;
  if (tt.eps[0])  ok = FALSE ;
  return ok ;
} /* ppIsNumber */

/*******************************************************************************************/

static void cleanTtSqrt (TT *ttp)
{
  if (ttp->sqrt1 == 0) ttp->sqrt1 = 1 ;
  if (ttp->sqrt2 == 0) ttp->sqrt2 = 1 ;
  
  ttp->sqrti = (ttp->sqrti + 40000)%4 ;
  switch (ttp->sqrti)
    {
    case 0:  break ;
    case 1:  break ;
    case 2: ttp->sqrti = 0 ; ttp->z *= 1.0I ; break ;
    case 3: ttp->sqrti = 1 ; ttp->z *= 1.0I ; break ;
    }

  if (*ttp->theta)
    {
      BOOL more = TRUE ;
      int sign = 1 ;
      while (more)
	{
	  more = FALSE ;
	  for (char *cp = ttp->theta ; ! more && *cp ; cp++)
	    for (char *cq = cp + 1 ; !more && *cq ; cq++)
	      {
		if (*cp == *cq)
		  { ttp->z = 0 ; return ; }
		if (*cp > *cq)
		  { char cc = *cp ; *cp = *cq ; *cq = cc ; more = TRUE ; sign = -sign ; break ;}
	      }
	}
      ttp->z *= sign ;
    }

  if (ttp->sqrt1 * ttp->sqrt2 > 1)
    {
      int i ;
      /* reduce to common denominator */
      for (i = 2 ; i <= ttp->sqrt1 && i <= ttp->sqrt2 ; i++)
	{
	  if ((ttp->sqrt1 % i == 0) && (ttp->sqrt2 % i == 0))
	    {
	      ttp->sqrt1 /= i ;
	      ttp->sqrt2 /= i ;
	      i = 1 ;
	      continue ;
	    }
	}
      /* remove squares */
      for (i = 2 ; i*i <= ttp->sqrt1 ; i++)
	{
	  if (ttp->sqrt1 % (i*i) == 0)
	    {
	      ttp->sqrt1 /= (i*i)  ;
	      ttp->z *= i ;
	      i = 1 ;
	      continue ;
	    }
	}
      for (i = 2 ; i*i <= ttp->sqrt2 ; i++)
	{
	  if (ttp->sqrt2 % (i*i) == 0)
	    {
	      ttp->sqrt2 /= (i*i)  ;
	      ttp->z /= i ;
	      i = 1 ;
	      continue ;
	    }
	}
    }
} /* cleanTtSqrt */


/***********************************************************************************************************************************************/

static int indexTtSort (short *cp, int dx, int sign)
{
  int i, j, k, ss = 1 ;
  BOOL modif = TRUE ;

  if (! *cp)
    return 1 ;
  
  while (modif)
    {   /* sort inside the goups of length dx (i.e. g_ba -> g_ab  eps_acbd -> - eps_abcd */
      modif = FALSE ;
      for (i = 0 ; i + dx - 1 < GMAX && cp[i + dx - 1] ; i += dx)
	for (j = 0 ; j < dx -1 ; j++)
	  if (cp[i+j+1]  && cp[i+j] > cp[i+j+1])
	    { k = cp[i+j] ; cp[i+j] = cp[i+j+1] ; cp[i+j+1] = k ; ss *= sign ; modif = TRUE ; }
    }
  
  modif = TRUE ;
  while (modif)
    {   /* sort the blocks (i.e. g_cdab -> g_abcd */
      modif = FALSE ;
      for (i = 0 ; i + 2 * dx - 1 < GMAX && cp[i] && cp[i+dx] ; i += dx)
	if (cp[i] > cp[i+dx])
	  for (j = 0 ; j < dx ; j++)
	    { k = cp[i+j] ; cp[i+j] = cp[i+j+dx] ; cp[i+j+dx] = k ; modif = TRUE ; }
    }
  return ss ;
} /* indextTtSort */

/***********************************************************************************************************************************************/
/* Einstein contraction rules */
static POLYNOME contractTtIndices (POLYNOME pp)
{
  AC_HANDLE h = pp->h ;
  double complex zz = 1 ;
  TT tt = pp->tt ;
  int kkkk = 0 ;
  /* sort and search repeated pair of indices inside the metric itself */
  if (tt.type)
    {
      int ii, i, j, k = 0 ;
      short *g, *s ;
      BOOL ok = FALSE ;

      while (! ok)
	{
	  ok = TRUE ;
	  kkkk++ ;
	  pp->tt = tt ;
	  if (0) 
	    {
	      printf ("inside kkk=%d\t", kkkk) ;
	      showPol (pp) ;
	    }
	  /* sort the indices */
	  for (i = 0 ; i < 4 ; i++)
	    tt.z *= indexTtSort (tt.mm[i], 1, 1) ;
	  tt.z *= indexTtSort (tt.g, 2, 1) ;
	  tt.z *= indexTtSort (tt.eps, 4, -1) ;
	  /* simplify repeated k indices */
	  for (i = 0 ; i < GMAX -1 ; i++)
	    if (tt.mm[0][i] && tt.mm[0][i] == tt.mm[0][i+1] && tt.denom[0]) 
	      {
		tt.denom[0]-- ;   /* divide by k^2 top and bottom */
		for (j = 0 ; i+j < GMAX -2 ; j++)
		  tt.mm[0][i+j] = tt.mm[0][i+j+2] ;
		tt.mm[0][GMAX-2] = 0 ; 
		tt.mm[0][GMAX-1] = 0 ;
		i-- ;  /* scan again the same position */
	      }

	  /* search repeated indices in a single metric */ 
	  for (i = 0, g = tt.g ; g[i] && i < GMAX ; i+=2)
	    { 
	      if (g[i] == g[i+1])
		{ zz *= 4 ; ok = FALSE ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; g[k++] = 0; g[k++] = 0 ; ok = FALSE ; }
	      if (g[i] > g[i+1]) /* switch : the Lorentz metric is Abelian */
		{ short cc = g[i] ; g[i] = g[i+1] ; g[i+1] = cc ; ok = FALSE ; }
	    }
	  if (! ok) continue ;
	  
	  /* search repeated indices between a metric and a sigma */
	  for (i = 0, g = tt.g ; ok && g[i] && i < GMAX - 2 ; i+=2)
	    for (j = 0, s = tt.sigma ; ok && s[j] && j < GMAX ; j++)
	      { 
		if (g[i] == s[j])
		  { s[j] = g[i+1] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
		else if (g[i+1] == s[j])
		  { s[j] = g[i] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
	      }
	  if (! ok) continue ;
	
	  /* search repeated indices between a metric and a sigma-bar */
	  for (i = 0, g = tt.g ; ok && g[i] && i < GMAX - 2 ; i+=2)
	    for (j = 0, s = tt.sigB ; ok && s[j] && j < GMAX ; j++)
	      { 
		if (g[i] == s[j])
		  { s[j] = g[i+1] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
		else if (g[i+1] == s[j])
		  { s[j] = g[i] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
	      }
	  if (! ok) continue ;
	    
	  /* search repeated indices between a metric and a momentum */
	  for (ii = 0 ; ii < 3 ; ii++)
	    for (i = 0, g = tt.g ; ok && g[i] && i < GMAX - 2 ; i+=2)
	      for (j = 0, s = tt.mm[ii] ; ok && s[j] && j < GMAX ; j++)
		{ 
		  if (g[i] == s[j])
		    { s[j] = g[i+1] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; g[k] = g[k+1] = 0 ; ok = FALSE ; }
		  else if (g[i+1] == s[j])
		    { s[j] = g[i] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; g[k] = g[k+1] = 0 ; ok = FALSE ; }
		}
	  if (! ok) continue ;
		  
	  /* search repeated indices between a metric and an epsilon */
	  for (i = 0, g = tt.g ; ok && g[i] && i < GMAX - 2 ; i+=2)
	    for (j = 0, s = tt.eps ; ok && s[j] && j < GMAX ; j++)
	      { 
		if (g[i] == s[j])
		  { s[j] = g[i+1] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; g[k] = g[k+1] = 0 ; ok = FALSE ; }
		else if (g[i+1] == s[j])
		  { s[j] = g[i] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; g[k] = g[k+1] = 0 ; ok = FALSE ; }
	      }
	  if (! ok) continue ; /* will sort the modified epsilon symbol */
		  
	  /* search repeated indices inside an epsilon */
	  for (i = 0, s = tt.eps ; ok && i < GMAX && s[i] ; i+= 4)
	    for (j = 0 ; j <= 2 ; j++)
	      for (k = j+1 ; k <= 3 ; k++)
		{
		  if (s[i+j] == s[i+k])
		    {
		      tt.z = 0 ;
		      pp->tt = tt ;
		      return pp ;
		    }
		}
	  
	  /* sort indices inside an epsilon : guaranteed by the above call to indexTtSort */

	  /* search repeated indices between pairs of epsilon */
	  for (int pass = 4 ; pass > 0 ; pass--)  /* only use n=0,2,4 in first pass */
	    for (i = 0, s = tt.eps ; 1 && ok && i < GMAX && s[i] ; i += 4)
	      for (j = i + 4 ; ok && j < GMAX && s[j] ; j += 4)
		{
		  int kkk[4], lll[4] ;
		  int n = 0, k, l, kk ;
		  short e = 0, f = 0 ;
		  short e1 = 0, e2 = 0, e3 = 0, f1 = 0, f2 = 0, f3 = 0 ;
		  for (k = 0 ; k < 4 ; k++)
		    { kkk[k] = lll[k] = 0 ; }
		  for (k = 0 ; k < 4 ; k++)
		    {
		      for (l = 0 ; l < 4 ; l++)
			if (s[i+k] == s[j+l])
			  { n++ ; kkk[k] = l+1 ; lll[l] = k + 1 ; }
		    }
		  if (n != pass)
		    continue ;
		  switch (n)
		    {
		    case 0:
		    case 10:
		      break ;
		    case 1:
		      kk = 0 ;
		      /* identify the non repeated indices */
		      for (k = 0 ; k < 4 ; k++)
			{
			  if (kkk[k] == 0)
			    {
			      if (! e1)
				e1 = s[i+k] ;
			      else if (! e2)
				e2 = s[i+k] ;
			      else
				e3 = s[i+k] ;
			    }
			  else
			    kk += k + kkk[k] - 1 ;
			}
		      for (k = 0 ; k < 4 ; k++)
			if (lll[k] == 0)
			  {
			    if (!f1)
			      f1 = s[j+k] ;
			    else if (!f2)
			      f2 = s[j+k] ;
			    else
			      f3 = s[j+k] ;
			  }
		      
		      g = tt.g ; l = 0 ;
		      while (*g) { g++ ; l++ ;}
		      *g++ = e1 ;
		      *g++ = f1 ;
		      *g++ = e2 ;
		      *g++ = f2 ;
		      *g++ = e3 ;
		      *g++ = f3 ;
		      *g++ = 0 ;
		      
		      tt.z *= 1 ; /* we contracted 1 index */
		      if (kk % 2) /* adjust the sign */
			tt.z *= -1 ;
		      
		      /* clean up the epsilons */
		      for (k = j ; k < GMAX - 4 ; k++)
			s[k] = s[k+4] ;
		      for (; k < GMAX - 4 ; k++)
			s[k] = 0 ;
		      for (k = i ; k < GMAX - 4 ; k++)
			s[k] = s[k+4] ;
		      for (; k < GMAX - 4 ; k++)
			s[k] = 0 ;
		      ok = FALSE ;
		      
		      /* sixplicate the polynome and anisymmetrize */
		      POLYNOME p1 = newScalar (1, h) ;
		      POLYNOME p2 = newScalar (1, h) ;
		      POLYNOME p3 = newScalar (1, h) ;
		      POLYNOME p4 = newScalar (1, h) ;
		      POLYNOME p5 = newScalar (1, h) ;
		      POLYNOME p6 = newScalar (1, h) ;
		      
		      p1->tt = tt ;
		      
		      p2->tt = tt ;
		      p2->tt.z *= -1 ;
		      p2->tt.g[l++] = e1 ;
		      p2->tt.g[l++] = f1 ;
		      p2->tt.g[l++] = e2 ;
		      p2->tt.g[l++] = f3 ;
		      p2->tt.g[l++] = e3 ;
		      p2->tt.g[l++] = f2 ;
		      l -= 6 ;
		      
		      p3->tt = tt ;
		      p3->tt.z *= -1 ;
		      p3->tt.g[l++] = e1 ;
		      p3->tt.g[l++] = f2 ;
		      p3->tt.g[l++] = e2 ;
		      p3->tt.g[l++] = f1 ;
		      p3->tt.g[l++] = e3 ;
		      p3->tt.g[l++] = f3 ;
		      l -= 6 ;
		      
		      p4->tt = tt ;
		      p4->tt.z *= 1 ;
		      p4->tt.g[l++] = e1 ;
		      p4->tt.g[l++] = f2 ;
		      p4->tt.g[l++] = e2 ;
		      p4->tt.g[l++] = f3 ;
		      p4->tt.g[l++] = e3 ;
		      p4->tt.g[l++] = f1 ;
		      l -= 6 ;
		      
		      p5->tt = tt ;
		      p5->tt.z *= 1 ;
		      p5->tt.g[l++] = e1 ;
		      p5->tt.g[l++] = f3 ;
		      p5->tt.g[l++] = e2 ;
		      p5->tt.g[l++] = f1 ;
		      p5->tt.g[l++] = e3 ;
		      p5->tt.g[l++] = f2 ;
		      l -= 6 ;
		      
		      p6->tt = tt ;
		      p6->tt.z *= -1 ;
		      p6->tt.g[l++] = e1 ;
		      p6->tt.g[l++] = f3 ;
		      p6->tt.g[l++] = e2 ;
		      p6->tt.g[l++] = f2 ;
		      p6->tt.g[l++] = e3 ;
		      p6->tt.g[l++] = f1 ;
		      l -= 6 ;
		      
		      POLYNOME ppp3[] = {p1, p2, p3, p4, p5, p6, 0} ;
		      POLYNOME p7 = polMultiSum (pp->h, ppp3) ;
		      ok = FALSE ;
		      *pp = *p7 ;
		      return pp ;
		      break ;
		    case 2:
		      kk = 0 ;
		    /* identify the non repeated indices */
		      for (k = 0 ; k < 4 ; k++)
			{
			  if (kkk[k] == 0)
			    {
			      if (! e1)
				e1 = s[i+k] ;
			      else
				e2 = s[i+k] ;
			    }
			  else
			    kk += k + kkk[k] - 1 ;
			}
		      for (k = 0 ; k < 4 ; k++)
			if (lll[k] == 0)
			  {
			    if (!f1)
			      f1 = s[j+k] ;
			    else
			      f2 = s[j+k] ;
			}
		      g = tt.g ; l = 0 ;
		      while (*g) { l++ ; g++ ;}
		      *g++ = e1 ;
		      *g++ = f1 ;
		      *g++ = e2 ;
		      *g++ = f2 ;
		      *g++ = 0 ;
		      
		      tt.z *= 2 ; /* we contracted 2 indices */
		      if (kk % 2) /* adjust the sign */
			tt.z *= -1 ;
		      /* clean up the epsilons */
		      for (k = j ; k < GMAX - 4 ; k++)
			s[k] = s[k+4] ;
		      for (; k < GMAX - 4 ; k++)
			s[k] = 0 ;
		      for (k = i ; k < GMAX - 4 ; k++)
			s[k] = s[k+4] ;
		      for (; k < GMAX - 4 ; k++)
			s[k] = 0 ;
		      ok = FALSE ;

		      /* duplicate the polynome and antisymmetrize */
		      p1 = newScalar (1,h) ;
		      p2 = newScalar (1,h) ;
		      
		      p1->tt = tt ;
		      p2->tt = tt ;
		      p2->tt.g[l++] = e1 ;
		      p2->tt.g[l++] = f2 ;
		      p2->tt.g[l++] = e2 ;
		      p2->tt.g[l++] = f1 ;
		      
		      p2->tt.z *= -1 ;
		      p3 = polSum (p1, p2,h) ;
		      ok = FALSE ;
		      *pp = *p3 ;
		      return pp ;
		      
		      break ;
		    case 3:
		      kk = 0 ;
		      /* identify the non repeated indices */
		      for (k = 0 ; k < 4 ; k++)
			{
			  if (kkk[k] == 0)
			    e = s[i+k] ;
			  else
			    kk += k + kkk[k] - 1 ;
			}
		      for (k = 0 ; k < 4 ; k++)
			if (lll[k] == 0)
			  f = s[j+k] ;
		      g = tt.g ;
		      while (*g) g++ ;
		      *g++ = e ;
		      *g++ = f ;
		      *g++ = 0 ;
		      tt.z *= 6 ; /* we contracted 3 indices */
		      if (kk % 2 == 1) /* adjust the sign */
			tt.z *= -1 ;
		      
		      /* clean up the epsilons */
		      for (k = j ; k < GMAX - 4 ; k++)
			s[k] = s[k+4] ;
		      for (; k < GMAX - 4 ; k++)
			s[k] = 0 ;
		      for (k = i ; k < GMAX - 4 ; k++)
			s[k] = s[k+4] ;
		      for (; k < GMAX - 4 ; k++)
			s[k] = 0 ;
		      ok = FALSE ;
		      break ;
		    case 4:
		      tt.z *= 24 ; /* we contracted 4 indices */
		      for (k = j ; k < GMAX - 4 ; k++)
			s[k] = s[k+4] ;
		      for (; k < GMAX - 4 ; k++)
			s[k] = 0 ;
		      for (k = i ; k < GMAX - 4 ; k++)
			s[k] = s[k+4] ;
		      for (; k < GMAX - 4 ; k++)
			s[k] = 0 ;
		      pp->tt = tt ;
		      ok = FALSE ;
		      break ;
		    }
		}
	  if (! ok) continue ;
	  

	  /* search repeated indices between pairs of pauli matrices and epsilon */
	  for (i = 0, g = tt.sigma ; ok && g[i] && i < GMAX ; i++)
	    for (j = 0, s = tt.eps ; ok && s[j] && j < GMAX ; j+= 4)
	      {
		int n = 0, k ;
		for (k = 0 ; k < 4 ; k++)
		  if (g[i] == s[j+k] || g[i+1] == s[j+k])
		      n++ ;
		if (n == 2) /* replace the pauli matrices and eliminate the epsilon */
		  {
		    POLYNOME p1, p2, p3 ;
		    int e = 0, f = 0, kk = 0 ;
		    tt.z *= 1 ;     /* I in Minkovski */
		    if (0 && g[i] > g[i+1])
		      tt.z *= -1 ;  /* 2025_04)08: wrong was breaking the tensor contribution to the fermion-vector ward identity */
		    for (k = 0 ; k < 4 ; k++)
		      if (g[i] != s[j+k] && g[i+1] != s[j+k])
			{
			  if (e == 0)
			    e = s[j+k] ;
			  else
			    f = s[j+k] ;
			  kk += k ;
			}
		    if ((i+kk) % 2 == 0)
		      tt.z *= -1 ;
		    g[i] = e ;
		    g[i+1] = f ;

		    for (k = j ; k < GMAX - 4 ; k++)
		      s[k] = s[k+4] ;
		    for (; k < GMAX - 4 ; k++)
		      s[k] = 0 ;

		    /* duplicate the polynome and antisymmetrize */
		    p1 = newScalar (1,h) ;
		    p2 = newScalar (1,h) ;

		    p1->tt = tt ;
		    p2->tt = tt ;
		    p2->tt.sigma[i] = f ;
		    p2->tt.sigma[i+1] = e;
		    p2->tt.z *= -1 ;
		    p3 = polSum (p1, p2,h) ;
		    ok = FALSE ;
		    *pp = *p3 ;
		    return pp ;
		    break ;
		  }
	      }
	  if (! ok) continue ;



	  /* search repeated indices between pairs of pauliBar matrices and epsilon */
	  for (i = 0, g = tt.sigB ; ok && g[i] && i < GMAX ; i++)
	    for (j = 0, s = tt.eps ; ok && s[j] && j < GMAX ; j+= 4)
	      {
		int n = 0, k ;
		for (k = 0 ; k < 4 ; k++)
		  if (g[i] == s[j+k] || g[i+1] == s[j+k])
		      n++ ;
		if (n == 2) /* replace the pauli matrices and eliminate the epsilon */
		  {
		    POLYNOME p1, p2, p3 ;
		    int e = 0, f = 0, kk = 0 ;
		    tt.z *= -1 ;    /* -I in Minkovski */
		    if (0 && g[i] > g[i+1])
		      tt.z *= -1 ; /* 2025_04)08: wrong was breaking the tensor contribution to the fermion-vector ward identity */
		    for (k = 0 ; k < 4 ; k++)
		      if (g[i] != s[j+k] && g[i+1] != s[j+k])
			{
			  if (e == 0)
			    e = s[j+k] ;
			  else
			    f = s[j+k] ;
			  kk += k ;
			}
		    if ((i+kk) % 2 == 0)
		      tt.z *= -1 ;
		    g[i] = e ;
		    g[i+1] = f ;

		    for (k = j ; k < GMAX - 4 ; k++)
		      s[k] = s[k+4] ;
		    for (; k < GMAX - 4 ; k++)
		      s[k] = 0 ;

		    p1 = newScalar (1,h) ;
		    p2 = newScalar (1,h) ;

		    p1->tt = tt ;
		    p2->tt = tt ;
		    p2->tt.sigB[i] = f ;
		    p2->tt.sigB[i+1] = e;
		    p2->tt.z *= -1 ;
		    p3 = polSum (p1, p2,h) ;
		    ok = FALSE ;
		    *pp = *p3 ;
		    return pp ;
		    break ;
		  }
	      }
	  if (! ok) continue ;


	  /* search pairs of repeated indices between epsilon and momenta */
	  for (ii = 0 ; ii < 3 ; ii++)
	    for (i = 0, s = tt.eps ; ok && s[i] && i < GMAX ; i+= 4)
	      for (j = 0 ; j <= 3 ; j++)
		{
		  int k, l, m ; short *t, *u ;
		  for (k = 0, t = tt.mm[ii] ; ok && t[k] && k < GMAX ; k++)		
		    if (s[i+j] == t[k])
		      for (l = 0, u = tt.mm[ii] ; ok && s[l] && l < GMAX ; l++)		
			for (m = 0 ; m <= 3 ; m++)
			  if (k != l && s[i+m] == u[l])
			    {
			      tt.z = 0 ;
			      pp->tt = tt ;
			      return pp ;
			    }
		}
	  /* search repeated indices in a pair of metrics, do one modif at a time */
	  for (i = 0, g = tt.g ; ok && g[i] && i < GMAX - 2 ; i+=2)
	    for (j = i + 2 ; ok && g[j] && j < GMAX ; j+= 2)
	      { 
		if (g[i] == g[j])
		  { g[i] = g[j+1] ; for (k = j ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
		else if (g[i] == g[j+1])
		  { g[i] = g[j] ; for (k = j ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
		else if (g[i+1] == g[j])
		  { g[i+1] = g[j+1] ; for (k = j ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
		else if (g[i+1] == g[j+1])
		  { g[i+1] = g[j] ; for (k = j ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
	      }
	  if (! ok) continue ;
	  
	  /* search repeated indices between a metric and epsilon */
	  for (i = 0, g = tt.g ; ok && g[i] && i < GMAX - 2 ; i+=2)
	    for (j = 0, s = tt.eps ; ok && s[j] && j < GMAX ; j++)
	      { 
		if (g[i] == s[j])
		  { s[j] = g[i+1] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; g[k] = g[k+1] = 0 ; ok = FALSE ; }
		if (g[i+1] == s[j])
		  { s[j] = g[i] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; g[k] = g[k+1] = 0 ; ok = FALSE ; }
	      }
	  if (! ok) continue ;
	
	      
	  /* sort alphabetically the momenta, they are Abelian */
	  for (ii = 0 ; ii < 3 ; ii++)
	    for (i = 0, s = tt.mm[ii] ; ok && s[i] && i < GMAX ; i++)
	      for (j = i + 1 ; s[j] && j < GMAX ; j++)
		if (s[j] < s[i])
		  { short cc = s[i] ; s[i] = s[j] ; s[j] = cc ; }
	  
	  /* search sigma or sigB repeated indices
	   * sigma_a sigB_a = 4, AND  s_a s_b s_a = -2 s_b AND  abEab = 4 E AND s_a sB_b s_c sB_a = 4 g_bc
	   */
	  for (i = 0, s = (tt.sigma[0] ? tt.sigma : tt.sigB) , g = tt.g ; ok && s[i] && i < GMAX - 3 ; i++)
	    { 
	      int u = 0, v = 0 ;
	      if (s[i] == s[i+1])
		{ zz *= 4 ; u = 0 ; v = 2 ; ok = FALSE ; }
	      else if (s[i] == s[i+2])
		{ zz *= -2 ; s[i] = s[i+1] ; u = 1 ; v = 2 ;  ok = FALSE ; }
	      else if (s[i] == s[i+3])
		{ zz *= 4 ;
		  k = strlen ((const char *)g) ; if (k > GMAX-3) messcrash ("k=%d too large", k) ;
		  g[k] = s[i+1]; g[k+1] = s[i+2] ; g[k+2] = 0 ;
		  u = 0 ; v = 4 ;  ok = FALSE ;
		}
	      else if (s[i] == s[i+4] && s[i+1] == s[i+5]  && s[i+2] == s[i+6])
		{ zz *= -32 ; s[i] = s[i+3] ; u = 1 ; v = 6 ; ok = FALSE ; }
	      if (! ok)
		{
		  for (k = i + u ; k < GMAX - v ; k++)
		    s[k] = s[k+v] ;
		  for (; k < GMAX ; k++)
		    s[k] = 0 ;
		}
	    }
	  if (! ok) continue ;
	}
    }
  tt.z *= zz ;
  pp->tt = tt ;
  return pp ;
}

POLYNOME contractIndices (POLYNOME pp)
{
  POLYNOME p1, p2 ;
  static int kkk ;
  
  if (!pp)
    return 0 ;
  p1 = pp->p1 ;
  p2 = pp->p2 ;
  p1 = contractIndices (p1) ;
  p2 = contractIndices (p2) ;

  kkk++ ;
  if (0)
    {
      printf ("outside kkk=%d\t", kkk) ;
      showPol (pp) ;
    }
  
  if (pp->tt.type)
    {
      pp = contractTtIndices (pp) ;
      if (0)
	{
	  printf ("..............outside kkk=%d\t", kkk) ;
	  if (pp) showPol (pp) ;
	}
      if (pp->tt.type)
	{
	  if (pp->tt.z == 0)
	    return 0 ;
	}
      else
	pp = contractIndices (pp) ;
    }
  return pp ;
}

static int slen (short *sp)
{
  int n = 0 ;
  while (*sp++)
    n++ ;
  return n ;
}

/* in a product of monomes, the list of symbols must merge */
static int contractTTProducts (POLYNOME pp, POLYNOME p1, POLYNOME p2)
{
  short *u, *v, *w ;
  short buf[GMAX] ;
  int ii, i, j ;
  TT tt = pp->tt ;
  TT t1 = p1->tt ;
  TT t2 = p2->tt ;

  /* merge numbers */

  cleanTtSqrt (&t1) ;
  cleanTtSqrt (&t2) ;
  tt.z = t1.z * t2.z ;
  tt.sqrt1 = t1.sqrt1 * t2.sqrt1 ;
  tt.sqrt2 = t1.sqrt2 * t2.sqrt2 ;
  tt.sqrti = t1.sqrti + t2.sqrti ;
  {
    char *cp = tt.theta ;
    char *cp1 = t1.theta ;
    char *cp2 = t2.theta ;
    while (*cp1)
      *cp++ = *cp1++ ;
    while (*cp2)
      *cp++ = *cp2++ ;
    *cp = 0 ;
  }
  cleanTtSqrt (&tt) ;

  
  if (cabs (tt.z) < minAbs)
    {
      if (0)
	{
	  ac_free (p1) ;
	  ac_free (p2) ;
	}
      pp->isFlat = FALSE ;
      return 0 ;
    }

  /* merge symbols */
  if (1)
    {
      int i = 0 ;
      char *u = tt.x ;
      char *v = t1.x ;
      char *w = t2.x ;
      while (*v)
	{ *u++ = *v++ ; i++ ; }
      while (*w)
	{ *u++ = *w++ ; i++ ; }
      for (;i < GMAX ; i++)
	*u++ = 0 ;
      tt.N = t1.N + t2.N ;
    }
  if (1)   /* symbols are commuting, needed for the SU21 determinant (2024_2), but wrong when exponentiating formal matrices (2023) */
    {
      BOOL ok = FALSE ;

      while (! ok)
	{
	  ok = TRUE ;
	  for (char *u = tt.x ; *u && ok ; u++)
	    for (char *v = u + 1 ; *v && ok ; v++)
	      if (*u > *v) 
		{
		  char cc = *u ; *u = *v ; *v = cc ; ok = FALSE ; 
		}
	}
    }
  
  /* count the chi to know if at the end we get one or zero x */
  if (1)
    {
      int nx = 0, s = 1, i = 0 ;
      char *v, *u = tt.x ;
      while (*u)
	if (*u++ == 'x') 
	  {
	    nx++ ;
	    v = u ;
	    while (*v)
	      {
		/* evaluate the sign search the sign */
		if (*v >= 'i' && *v < 'm')
		  s = -s ; 
		v++ ;
	      }
	  }
      tt.z *= s ;
      /* copy the non chi symbols */
      u = v = tt.x ; 
      while (*v)
	{
	  if (*v != 'x')
	    { *u++ = *v ; i++ ; }
	  v++ ;
	}
      if (nx %2)
	{ *u++ = 'x' ; i++ ; }
      for (;i < GMAX ; i++)
	*u++ = 0 ;

      /* kill the odd square */
      u = tt.x ; 
      while (*u)
	{
	  if (u[0] == u[1] && u[0] >= 'i' && u[0] <= 'm')
	    {
	      tt.z = 0 ;
	      tt.type = 1 ;
	    }
	  u++ ;
	}
    }

  /* merge denoms */
  for (i = 0 ; i < 4 ; i++)
    {
      if (tt.denom[i] != t1.denom[i] + t2.denom[i])
	pp->isFlat = FALSE ;
      tt.denom[i] = t1.denom[i] + t2.denom[i] ;
    }

  /* merge metrics */
  u = t1.g ; v = t2.g ; w = tt.g ;
  memcpy (buf, w, GMAX) ;
  i = slen (u) ; j = slen (v) ;
  if (i+j)
    {
      if (i+j >= GMAX) messcrash ("Merging too long metrics: %s %s", u, v) ;
      while ((*w++ = *u++)) ;  
      w-- ; 
      while ((*w++ = *v++)) ;
    }
  if (memcmp (buf, w, GMAX))
    pp->isFlat = FALSE ;

  u = t1.gg ; v = t2.gg ; w = tt.gg ;
  i = slen (u) ; j = slen (v) ;
  memcpy (buf, w, GMAX) ;
  if (i+j)
    {
      if (i+j >= GMAX) messcrash ("Merging too long metrics: %s %s", u, v) ;
      while ((*w++ = *u++)) ;
      w-- ; while ((*w++ = *v++)) ;
    }
  if (memcmp (buf, w, GMAX))
    pp->isFlat = FALSE ;

  /* merge epsilon */
  u = t1.eps ; v = t2.eps ; w = tt.eps ;
  memcpy (buf, w, GMAX) ;
  i = slen (u) ; j = slen (v) ;
  if (i+j)
    {
      if (i+j >= GMAX) messcrash ("Merging too long metrics: %s %s", u, v) ;
      while ((*w++ = *u++)) ;  
      w-- ; 
      while ((*w++ = *v++)) ;
    }
  if (memcmp (buf, w, GMAX))
    pp->isFlat = FALSE ;

  /* merge momenta */
  for (ii = 0 ; ii < 4 ; ii++)
    {
      memcpy (buf, w, GMAX) ;
      u = t1.mm[ii] ; v = t2.mm[ii] ; w = tt.mm[ii] ;
      i = slen (u) ; j = slen (v) ;
      if (i+j)
	{
	  memset (w, 0, GMAX) ;
	  if (i+j >= GMAX) messcrash ("Merging too long metrics: %s %s", u, v) ;
	  while ((*w++ = *u++)) ;
	  w-- ; while ((*w++ = *v++)) ;
	}
      if (memcmp (buf, w, GMAX))
	pp->isFlat = FALSE ;
    }

  /* merge matrices */
  if (t1.Id2 || t2.Id2)
    tt.Id2 = 1 ;
  if (t1.sigma[0] && t1.sigB[0])
    messcrash ("cannot have a sigma sigma product, should be sigma sigma->bar") ;
  else if (t2.sigma[0] && t2.sigB[0])
    messcrash ("cannot have a sigma sigma product, should be sigma sigma->bar") ;
  else if (! t1.sigma[0] && ! t1.sigB[0])
    { 
      w = tt.sigma ; v = t2.sigma ; 
      memcpy (buf, w, GMAX) ;
      while ((*w++ = *v++)) ; 
      if (memcmp (buf, w, GMAX))
	pp->isFlat = FALSE ;

      w = tt.sigB ;  v = t2.sigB  ; 
      memcpy (buf, w, GMAX) ;
      while ((*w++ = *v++)) ; 
      if (memcmp (buf, w, GMAX))
	pp->isFlat = FALSE ;
    }
  else if (! t2.sigma[0] && ! t2.sigB[0])
    {
      w = tt.sigma ; u = t1.sigma ;
      memcpy (buf, w, GMAX) ;
      while ((*w++ = *u++)) ; 
      if (memcmp (buf, w, GMAX))
	pp->isFlat = FALSE ;

      w = tt.sigB ;  u = t1.sigB  ;
      memcpy (buf, w, GMAX) ;
      while ((*w++ = *u++)) ; 
      if (memcmp (buf, w, GMAX))
	pp->isFlat = FALSE ;
    }
  else if (t1.sigma[0] && (slen(t1.sigma) % 2) && t2.sigma[0])
    messcrash ("cannot have a (odd-sigma) sigma product, should be (odd-sigma) sigmaBar") ;
  else if (t1.sigma[0] && (slen(t1.sigma) % 2 == 0) && t2.sigB[0])
    messcrash ("cannot have a (even-sigma) sigmaBar product, should be (even-sigma) sigma") ;
  else if (t1.sigB[0] && (slen(t1.sigB) % 2) && t2.sigB[0])
    messcrash ("cannot have a (odd-sigmaBar) sigmaBar product, should be (odd-sigmaBar) sigma") ;
  else if (t1.sigB[0] && (slen(t1.sigB) % 2 == 0) && t2.sigma[0])
    messcrash ("cannot have a (even-sigmaBar) sigma product, should be (even-sigmaBar) sigmaBar") ;
 else if (t1.sigma[0])
    { 
      u = t1.sigma ; i = slen (u) ; w = tt.sigma ;
      memcpy (buf, w, GMAX) ;
      if (i % 2) 
	{
	  v = t2.sigB ; j = slen (v) ;
	  if (i+j >= GMAX) messcrash ("Merging too long metrics: %s %s", u, v) ;
	  while ((*w++ = *u++)) ;  
	  w-- ; 
	  while ((*w++ = *v++)) ;
	}
      else
	{
	  v = t2.sigma; j = slen (v) ;
	  if (i+j >= GMAX) messcrash ("Merging too long metrics: %s %s", u, v) ;
	  while ((*w++ = *u++)) ;
	  w-- ; while ((*w++ = *v++)) ;
	}
      if (!tt.sigma[0])
	tt.Id2 = 1 ;
      if (memcmp (buf, w, GMAX))
	pp->isFlat = FALSE ;
    }
  else if (t1.sigB[0])
    {
      u = t1.sigB ; i = slen (u) ; w = tt.sigB ;
      memcpy (buf, w, GMAX) ;
      if (i % 2 == 0) 
	{
	  v = t2.sigB ; j = slen (v) ;
	  if (i+j >= GMAX) messcrash ("Merging too long metrics: %s %s", u, v) ;
	  while ((*w++ = *u++)) ;
	  w-- ; while ((*w++ = *v++)) ;
	}
      else 
	{
	  v = t2.sigma; j = slen (v) ;
	  if (i+j >= GMAX) messcrash ("Merging too long metrics: %s %s", u, v) ;
	  while ((*w++ = *u++)) ;
	  w-- ; while ((*w++ = *v++)) ;
	}
      if (!tt.sigB[0])
	tt.Id2 = 1 ;
      if (memcmp (buf, w, GMAX))
	pp->isFlat = FALSE ;
    }

  tt.type = 1 ;
  pp->tt = tt ;
  return tt.type ;
}
  
/*******************************************************************************************/

POLYNOME contractProducts (POLYNOME pp)
{
  POLYNOME p1, p2 ;
  BOOL debug = FALSE ;
  static int nn= 0 ;
  if (!pp)
    return 0 ;
  polCheck (pp) ;
  nn++ ;
  p1 = pp->p1 ;
  p2 = pp->p2 ;
  p2 = pp->p2 = contractProducts (p2) ;
  p1 = pp->p1 = contractProducts (p1) ;

  if (debug)  showPol(pp) ;
  if (0 && pp->isProduct && p1 && p2 && p1->tt.type && p2->isProduct && p2->p1 && p2->p1->tt.type)
      contractTTProducts (pp, p1, p2->p1) ;


  if (pp->isSum && p1 && !p2)
    {
      *pp = *p1 ;
      pp->isFlat = FALSE ;
    }

  if (pp->isSum && p2 && !p1)
    {
      *pp = *p2 ;
      pp->isFlat = FALSE ;
    }

  if (pp->isSum && p1 && p2 && ! p1->tt.type && p2->tt.type)
    { /* addition is Abelian */
      POLYNOME q = pp->p1 ;
      pp->p1 = pp->p2 ;
      pp->p2 = q ;
      pp->isFlat = FALSE ;
    }
  
  if (pp->isSum && p1 && p2 && p1->tt.type && p2->tt.type)
    {
      double complex z1 = p1->tt.z ;
      double complex z2 = p2->tt.z ;
      int s ;

      cleanTtSqrt (&(p1->tt)) ;
      cleanTtSqrt (&(p2->tt)) ;
      
      s = polOrder (&p1, &p2) ;

      if (s == 0)
	{
	  *pp = *p1 ;
	  pp->isFlat = FALSE ;
	  pp->tt.z = z1 + z2 ;
	}
      else if (s > 0)
	{ /* addition is Abelian */
	  POLYNOME q = pp->p1 ;
	  pp->p1 = pp->p2 ;
	  pp->p2 = q ;
	  pp->isFlat =  FALSE ;
	}
    }
  

  if (pp->isSum && p1 && p2 && p1->tt.type && p2->isSum && p2->p1 && p2->p1->tt.type)
    {
      cleanTtSqrt (&(p1->tt)) ;
      cleanTtSqrt (&(p2->p1->tt)) ;

      double complex z1 = p1->tt.z ;
      double complex z2 = p2->p1->tt.z ;
      int s ;

      p1->tt.z = 0 ;       
      p2->p1->tt.z = 0 ;       
      s = memcmp (&(p1->tt), &(p2->p1->tt), sizeof(TT)) ;
      p1->tt.z = z1 ;       
      p2->p1->tt.z = z2 ;       

      if (s == 0)
	{
	  *pp = *p2 ;
	  pp->isFlat = FALSE ;
	  pp->p1->tt.z = z1 + z2 ;
	}
      else if (s > 0)
	{ /* addition is Abelian */
	  POLYNOME q = pp->p1 ;
	  pp->p1 = pp->p2->p1 ;
	  pp->p2->p1 = q ;
	  pp->isFlat =  FALSE ;
	}
    }
  

  if (pp->isProduct && !p1 && !p2)
    pp = 0 ;

  if (pp->isProduct && p1 && p2 && p1->tt.type && p2->tt.type)
    {
      pp->isProduct = FALSE ;
      pp->isFlat = TRUE ;
      contractTTProducts (pp, p1, p2) ;
      if (0)
	{
	  ac_free (p1) ;
	  ac_free (p2) ;
	}
      pp->p1 = pp->p2 = 0 ;
      if (debug) showPol(pp) ;
      if (pp->tt.z == 0)
	pp = 0 ;
    }
  polCheck (pp) ;
  return pp ;
}
  
/*******************************************************************************************/
/*******************************************************************************************/

static KEYSET polynomeKs = 0 ;
static void checkPolynome (POLYNOME pp) 
{
  static int level = 0 ;
  static int nn ;
  int i ;

  if (! pp)
    return ;
  polCheck (pp) ;
  if (level == 0)
    {
      if (! polynomeKs)
	polynomeKs = keySetCreate () ;
      nn = 0 ;
    }
  level++ ;
  for (i = 0 ; i < nn ; i++)
    if (keySet (polynomeKs, i) == pp->id) 
      messcrash ("Duplicate node in polynome")  ;
  keySet (polynomeKs, nn++) = pp->id ; 
	
  if (pp->p1 == pp) messcrash ("pp == pp->p1 in checkPolynome") ;
  if (pp->p2 == pp) messcrash ("pp == pp->p2 in checkPolynome") ;
  if (pp->p1 && pp->p1 == pp->p2) messcrash ("pp->p1 == pp->p2 in checkPolynome") ;
  if (pp->p1 && pp->p2 && pp->p1->p2 == pp->p2) messcrash ("pp->p1->p2 == pp->p2 in checkPolynome") ;
  if (pp->p1) checkPolynome (pp->p1) ;
  if (pp->p2) checkPolynome (pp->p2) ;

  level-- ;
}

/*******************************************************************************************/
/* Flatten a polynome */
static POLYNOME expandDo (POLYNOME pp, int force)
{
  POLYNOME p1, p2 ;
  BOOL debug = FALSE ;
  AC_HANDLE h = pp ? pp->h : 0 ;
  
  if (!pp)
    return 0 ;
  polCheck (pp) ;

  if (force) pp->isFlat = FALSE ;
  p1 = pp->p1 ;
  p2 = pp->p2 ;
  if (p1 == pp) messcrash ("pp == pp->p1 in expandDo") ;
  if (p2 == pp) messcrash ("pp == pp->p2 in expandDo") ;
  if (p1 && p1 == p2) messcrash ("pp->p1 == pp->p2 in expandDo") ;

  if (pp->tt.type && cabs (pp->tt.z) < minAbs)
    return 0 ;
  if (p1 && p1->tt.type && cabs (p1->tt.z) < minAbs)
    p1 = pp->p1 = 0 ;
  if (p2 && p2->tt.type && cabs (p2->tt.z) < minAbs)
    p2 = pp->p2 = 0 ;
  
  if (pp->tt.type && cabs (pp->tt.z) < minAbs)
    { polFree (pp) ; return 0 ; }
  if (p1 && p1->tt.type && cabs (p1->tt.z) < minAbs)
    { polFree (p1) ; p1 = pp->p1 = 0 ; }
  if (p2 && p2->tt.type && cabs (p2->tt.z) < minAbs)
    { polFree (p2) ; p2 = pp->p2 = 0 ; }
  
  if (! pp)
    return 0 ;

  if (pp->isProduct)
    {
      if (!p1 || !p2)
	return 0 ;
    }
  if (pp->isSum)
    {
      if (!p1 && !p2)
	return 0 ;
      if (!p2)
	return expandDo (pp->p1, force) ;
      if (!p1)
	return expandDo (pp->p2, force) ;
    }

  if (pp->p2 && (force || ! pp->p2->isFlat))
    p2 = pp->p2 = expandDo (pp->p2, force) ;
  if (pp->p1 && (force || ! pp->p1->isFlat))
    p1 = pp->p1 = expandDo (pp->p1, force) ;
  if (debug) checkPolynome (pp) ;

  if (pp->isProduct) /* check again */
    {
      if (!p1 || !p2)
	return 0 ;
    }
  if (pp->isSum)
    {
      if (!p1 && !p2)
	return 0 ;
      if (!p2)
	return expandDo (pp->p1, force) ;
      if (!p1)
	return expandDo (pp->p2, force) ;
    }

  if (pp->isProduct && p1 && p1->isSum)
    {
      POLYNOME q1 = newPolynome (h) ;
      POLYNOME q2 = newPolynome (h) ;
      pp->isSum = TRUE ;
      pp->isProduct = FALSE ;
      q1->isProduct = TRUE ;
      q2->isProduct = TRUE ;
      q1->p2 = p2 ;
      q2->p2 = polCopy (p2, h) ;
      q1->p1 = p1->p1 ;
      q2->p1 = p1->p2 ;
      pp->p1 = q1 ;
      pp->p2 = q2 ;
      if (debug) checkPolynome (pp) ;
      if (pp->p2 && (force || ! pp->p2->isFlat))
	p2 = pp->p2 = expandDo (q2, force) ;
      if (pp->p1 && (force || ! pp->p1->isFlat))
	p1 = pp->p1 = expandDo (q1, force) ;
      if (debug) checkPolynome (pp) ;
    }
  if (pp->isProduct && p2 && p2->isSum)
    {
      POLYNOME q1 = newPolynome (h) ;
      POLYNOME q2 = newPolynome (h) ;
      pp->isSum = TRUE ;
      pp->isProduct = FALSE ;
      q1->isProduct = TRUE ;
      q2->isProduct = TRUE ;
      q1->p1 = p1 ;
      q2->p1 = polCopy (p1, h) ;
      q1->p2 = p2->p1 ;
      q2->p2 = p2->p2 ;
      if (q2 && q2->p1 && q2->p1 == q2->p2) messcrash ("q2->p1 == q2->p2 in expandDo") ;
      p1 = pp->p1 = q1 ;
      p2 = pp->p2 = q2 ;
      if (debug) checkPolynome (pp) ;
      if (pp->p2 && (force || ! pp->p2->isFlat))
	p2 = pp->p2 = expandDo (q2, force) ;
      if (pp->p1 && (force || ! pp->p1->isFlat))
	p1 = pp->p1 = expandDo (q1, force) ;
      if (debug) checkPolynome (pp) ;
    }
  if (pp->isSum && p1 && p2 && p1->isSum && p1->p2 && p2->isSum)
    {
      if (debug) checkPolynome (pp) ;
      p2 = pp->p2 = expandDo (pp->p2, force) ;
      if (debug) checkPolynome (p2) ;
      p1 = pp->p1 = expandDo (pp->p1, force) ;
      if (debug) checkPolynome (p1) ;
    }

  if (pp->isSum && p1 && p2 && p1->isSum && p1->p2 && p2->isSum)
    {
      pp->p1 = p1->p1 ;
      p1->p1 = p1->p2 ;
      p1->p2 = pp->p2 ;
      pp->p2 = p1 ;
      p1 = pp->p1 ; p2 = pp->p2 ;
      if (debug) checkPolynome (pp) ;
    }

  if (pp->isSum && p1 && p1->tt.type && p2 && p2->isSum && p2->p1 && p2->p1->tt.type)
    {
      POLYNOME p3 = pp->p2->p1 ;
      double complex z1 = p1->tt.z ;
      double complex z3 = p3->tt.z ;
      int s ;
      p1->tt.z = 0 ;       
      p3->tt.z = 0 ;       
      s = memcmp (&(p1->tt), &(p3->tt), sizeof(TT)) ;
      p1->tt.z = z1 ;       
      p3->tt.z = z3 ;       

      if (s == 0)
	{ /* add p1 inside p2, skip p1, return p2 */
	  p3->tt.z = z1 + z3 ;
	  p2->isFlat = FALSE ;
	  return expandDo (pp->p2, force) ;
	}
      else if (s > 0)
	{ /* addition is Abelian */
	  pp->p2->p1 = p1 ;
	  p1 = pp->p1 = p3 ;
	  pp->isFlat = FALSE ;
	  pp->p2->isFlat = FALSE ;
	}
      p1 = pp->p1 ; p2 = pp->p2 ;
    }
  

  if (! pp)
    return 0 ;

  if (pp->tt.type)
    cleanTtSqrt (&(pp->tt)) ;
    
  if (pp->tt.type && cabs (pp->tt.z) < minAbs)
    {
      return 0 ;
    }
  if (! pp->tt.type && !p1 && !p2)
    {
      return 0 ;
    }
  pp->isFlat = TRUE ;
  polCheck (pp) ;
  if (debug) checkPolynome (pp) ;
  return pp ;
} /* expandDo */

/*************************************************************************************************/

POLYNOME expand (POLYNOME pp)
{
  int force = 1 ;
  polCheck (pp) ;

  pp = sortPol (pp) ;
  if (pp)
    {
      int nn = 12 ;
      pp->isFlat = FALSE ;
      while (pp && ! pp->isFlat && nn-- > 0)
	{
	  pp = expandDo (pp, force) ;
	  pp = contractProducts (pp) ;
	  pp = contractIndices (pp) ;
	  pp = sortPol (pp) ;
	  force = FALSE ;
	}
    }
  polCheck (pp) ;
  return pp ;
}


/*******************************************************************************************/
/*******************************************************************************************/

static void gprintf (char *prefix, short *sp)
{
  printf (" %s_", prefix) ;
  sp-- ;
  while (*++sp)
    printf ("%c", (char)(*sp % 256)) ;
  printf (" ") ;
}

  static void showTT (POLYNOME pp)
{
  TT tt ;
  int i ; 
  
  if (!pp) 
    return ;
  polCheck (pp) ;
  tt = pp->tt ;  
  if (!tt.type)
    return ;
  
  cleanTtSqrt (&tt) ;

  if (1 && ! ppIsNumber (pp))
    {
      if (cabs (tt.z + 1) < minAbs)
	printf (" -") ;
      else if (cabs (tt.z - 1) > minAbs)
	{ tt.z = nicePrint ("", tt.z) ; }
    }
  else
    { tt.z = nicePrint ("", tt.z) ; }
  if (tt.sqrt1 > 1 || tt.sqrt2 > 1)
    {
      printf ("sqrt(%d",tt.sqrt1) ;
      if (tt.sqrt2 > 1) printf ("/%d",tt.sqrt2) ;
      printf (")") ;
    }
  switch (tt.sqrti)
    {
    case 0 :  break ;
    case 1:  printf(" eta ") ; break ;
    }


  if (*tt.x) printf (" %s ", tt.x) ;
  if (*tt.g) gprintf ("g", tt.g) ;
  if (*tt.gg) gprintf ("gg", tt.gg) ;
  if (*tt.eps) gprintf ("epsilon", tt.eps) ;

  if (*tt.sigma) gprintf ("s", tt.sigma) ;
  if (*tt.sigB) gprintf ("sB", tt.sigB) ;

  if (tt.mm[0][0]) gprintf ("k", tt.mm[0]) ;
  if (tt.mm[1][0]) gprintf ("p", tt.mm[1]) ;
  if (tt.mm[2][0]) gprintf ("q", tt.mm[2]) ;
  if (tt.mm[3][0]) gprintf ("r", tt.mm[3]) ;

  i = 0 ;
  if      (tt.denom[0]) printf (" %c k^%d ", i++ ? ' ' : '/', 2 * tt.denom[0]) ;
  if      (tt.denom[1]) printf (" %c (k+p)^%d ", i++ ? ' ' : '/', 2 * tt.denom[1]) ;
  if      (tt.denom[2]) printf (" %c (k+p+q)^%d ", i++ ? ' ' : '/', 2 * tt.denom[2]) ;
  if      (tt.denom[3]) printf (" %c (k+p+q+r)^%d ", i++ ? ' ' : '/', 2 * tt.denom[3]) ;

  if (*tt.theta) printf(" theta(%s) ", tt.theta) ;
  return ;
} /* showTT */

/*******************************************************************************************/

void showPol (POLYNOME pp)
{
  POLYNOME p1, p2 ;
  static int nn = 0, level = 0 ;

  polCheck (pp) ;
   if (!pp)
     {
       if (level == 0) printf ("(NULL) #### zero term ###\n") ;
       return ;
     }
   if (! level) nn = 0 ;
   level++ ;
   p1 = pp->p1 ;
   p2 = pp->p2 ;
   if (p1 == pp) messcrash ("pp == pp->p1 in showPol") ;
   if (p2 == pp) messcrash ("pp == pp->p2 in showPol") ;
   if (p1 && p1 == p2) messcrash ("pp->p1 == pp->p2 in showPol") ;

   if (pp->isProduct)
     { printf ("(");  showPol (p1) ; printf (")*(");  showPol(p2) ; printf(")") ; }
   else if (pp->isSum)
     { printf ("(");  showPol (p1) ; printf (")+") ; if (p2 && !p2->isSum) printf("(");  showPol(p2) ; if (p2 && !p2->isSum) printf(")") ; }
   else if (pp->tt.type)
     { showTT (pp) ; nn++ ; }
   level-- ;

   if (! level) printf ("  ### %d term%s ###\n", nn, nn > 1 ? "s" : "") ;
} /* showPol */

/***********************************************************************************************************************************************/
/******************************************* New Polynomes of all types ************************************************************************/
/***********************************************************************************************************************************************/

void polFree (POLYNOME pp)
{
  if (pp && pp->id)
    {
      pp->id = 0 ;
      pp->magic = 0 ;
      /*
	ac_free (pp->p1) ;
	ac_free (pp->p2) ;
      */
       /* do NOT       
	* ac_free (pp->h) ; 
	* pp may be a component of a larger pp with same h 
	*/
    }
  return ;
}

static void polFinalise (void *vp)
{
  POLYNOME pp = (POLYNOME) vp ;
  polFree (pp) ;
}

/***********************************************************************************************************************************************/

POLYNOME newPolynome (AC_HANDLE h)
{
  static int id = 0 ;
  int n = sizeof (struct polynomeStruct) ;
  POLYNOME pp = (POLYNOME) handleAlloc (polFinalise, h, n) ;
  memset (pp, 0, n) ;
  pp->id = ++id ;
  pp->h = h ;
  pp->magic = POLMAGIC ;
  return pp ;
}

/*************************************************************************************************/

POLYNOME polCopy (POLYNOME p1, AC_HANDLE h)
{
  POLYNOME pp = 0 ; 

  polCheck (p1) ;
  if (p1 && ! p1->isSum && ! p1->isProduct && ! p1->tt.z) p1 = 0 ;
  if (p1)
    {
      int id ;
      pp = newPolynome (h) ;
      id = pp->id ;
      memcpy (pp, p1, sizeof (struct polynomeStruct)) ;
      pp->id = id ;
      pp->h = h ;
      if (pp->p1)
	pp->p1 =  polCopy (pp->p1, h) ;
      if (pp->p2)
	pp->p2 =  polCopy (pp->p2, h) ;
    }
  return pp ;
}

/*************************************************************************************************/

POLYNOME newScalar (double complex z, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type = 1 ;
  p->tt.z = z ;
  return p ;
}

/*************************************************************************************************/

POLYNOME newG (short mu, short nu, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.g[0] = mu ;
  p->tt.g[1] = nu ;
  return p ;
}

/*************************************************************************************************/

POLYNOME newEpsilon (short a, short b, short c, short d, AC_HANDLE h)
{
  POLYNOME pp ;
  
  pp = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = 1.0 ;
  pp->tt.eps[0] = a ;
  pp->tt.eps[1] = b ;
  pp->tt.eps[2] = c ;
  pp->tt.eps[3] = d ;
  return pp ;
}

/*************************************************************************************************/
/* newAG (... ,  0) antisymmetric link 1/2(ac bd - ad bc). Optionally adding the i epsilon 
 * newAG (... , +1) is the self-dual projector P+ = 1/4 ( ac bd - ad bc) + i/2 epsilon(abcd)
 * newAG (... , -1) is the self-dual projector P- = 1/4 ( ac bd - ad bc) - i/2 epsilon(abcd)
 *     WE have (P+)^2 = (P+),   (P-)^2 = (P-),   (P+)(P-) = (P-)(P+) = 0 
 */
POLYNOME newAG (short a, short b, short c, short d, int parity, AC_HANDLE h)
{
  POLYNOME pp, ppp[4] ;

  pp = ppp[0] = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = 1.0/2.0 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = c ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = d ;

  pp = ppp[1] = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = -1.0/2.0 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = d ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = c ;

  ppp[2] = ppp[3] = 0 ;

  if (parity)
    {
      ppp[0]->tt.z /= 2 ;
      ppp[1]->tt.z /= 2 ;

      pp = ppp[2] = newEpsilon (a, b, d, c, h) ;
      pp->tt.z = (parity ) / 4.0 ;             /* I * parity,  in Minkovski */
      if (parity == 2 || parity == -2)
	return pp ;
    }
  
  pp = polMultiSum (h, ppp) ;
  return pp ; 
} /* newAG */

/*************************************************************************************************/

POLYNOME newAG6 (short a, short b, short c, short d, short e, short f, AC_HANDLE h)
{
  POLYNOME pp, ppp[7] ;

  pp = ppp[0] = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = 1 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = d ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = e ;
  pp->tt.g[4] = c ;
  pp->tt.g[5] = f ;

  pp = ppp[1] = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = -1 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = d ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = f ;
  pp->tt.g[4] = c ;
  pp->tt.g[5] = e ;

  pp = ppp[2] = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = -1 ;
  pp->tt.g[0] = b ;
  pp->tt.g[1] = e ;
  pp->tt.g[2] = a ;
  pp->tt.g[3] = f ;
  pp->tt.g[4] = c ;
  pp->tt.g[5] = d ;

  pp = ppp[3] = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = -1 ;
  pp->tt.g[0] = c ;
  pp->tt.g[1] = f ;
  pp->tt.g[2] = a ;
  pp->tt.g[3] = e ;
  pp->tt.g[4] = b ;
  pp->tt.g[5] = d ;

  pp = ppp[4] = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = 1 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = e ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = f ;
  pp->tt.g[4] = c ;
  pp->tt.g[5] = d ;

  pp = ppp[5] = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = 1 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = f ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = d ;
  pp->tt.g[4] = c ;
  pp->tt.g[5] = e ;

  ppp[6] = 0 ;
  
  pp = polMultiSum (h, ppp) ;
  return pp ; 
} /* newAG6 */

/*************************************************************************************************/

POLYNOME newK (short cc, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.mm[0][0] = cc ;
  return p ;
}

/*************************************************************************************************/

POLYNOME newP (short cc, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.mm[1][0] = cc ;
  return p ;
}

/*************************************************************************************************/

POLYNOME newQ (short cc, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.mm[2][0] = cc ;
  return p ;
}

/*************************************************************************************************/

POLYNOME newR (short cc, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.mm[3][0] = cc ;
  return p ;
}

/*************************************************************************************************/

POLYNOME newPQR (int pqr, short mu, AC_HANDLE h)
{
  POLYNOME p0 = pqr >=0 ? newK (mu, h) : 0 ;
  POLYNOME p1 = pqr >=1 ? newP (mu, h) : 0 ;
  POLYNOME p2 = pqr >=2 ? newQ (mu, h) : 0 ;
  POLYNOME p3 = pqr >=3 ? newR (mu, h) : 0 ;
  POLYNOME ppp[5] = {p0, p1, p2, p3, 0} ;

  return polMultiSum (h, ppp) ;
}

/*************************************************************************************************/

POLYNOME newSigma (short cc, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.Id2 = 1 ;
  p->tt.sigma[0] = cc ;
  return p ;
}

/*************************************************************************************************/

POLYNOME newSigB (short cc, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type = 1 ;
  p->tt.z = 1 ;
  p->tt.Id2 = 1 ;
  p->tt.sigB[0] = cc ;
  return p ;
} /* newSigB */

/*************************************************************************************************/

POLYNOME newTheta (char *cp, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  strncpy (p->tt.theta, cp, GMAX-2) ;
  return p ;
}

/*************************************************************************************************/

POLYNOME newSymbol (char *cp, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.N = strlen (cp) ;
  strncpy (p->tt.x, cp, GMAX-2) ;
  return p ;
}

/*************************************************************************************************/

static POLYNOME polDoSum (POLYNOME p1, POLYNOME p2, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  polCheck (p1) ;
  polCheck (p2) ;
  p->p1 = p1 ;
  p->p2 = p2 ;
  if (p1->tt.type && p2->tt.type && p1->tt.Id2 + p2->tt.Id2 == 1)
    messcrash ("Cannot add aPauli matrix to a number\n") ;;
  p->isSum = TRUE ;
  return p ;
} /* polSum */

/*************************************************************************************************/

POLYNOME polSum (POLYNOME p1, POLYNOME p2, AC_HANDLE h)
{
  if (p1 && ! p1->isSum && ! p1->isProduct && ! p1->tt.z) p1 = 0 ;
  if (p2 && ! p2->isSum && ! p2->isProduct && ! p2->tt.z) p2 = 0 ;
  polCheck (p1) ;
  polCheck (p2) ;
  if (p1)
    p1 = polCopy (p1, h) ;
  if (p2)
    p2 = polCopy (p2, h) ;
  if (p1 && p2)
    return polDoSum (p1, p2, h) ;
  else if (p1)
    return p1 ;
  else if (p2)
    return p2 ;
  return 0 ;
} /* polSum */

/*************************************************************************************************/

static POLYNOME polDoProduct (POLYNOME p1, POLYNOME p2, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->p1 = p1 ;
  p->p2 = p2 ;
  p->isProduct = TRUE ;
  polCheck (p1) ;
  polCheck (p2) ;
  return p ;
} /* polProduct */

/*************************************************************************************************/

POLYNOME polProduct (POLYNOME p1, POLYNOME p2, AC_HANDLE h)
{
  if (p1 && ! p1->isSum && ! p1->isProduct && ! p1->tt.z) p1 = 0 ;
  if (p2 && ! p2->isSum && ! p2->isProduct && ! p2->tt.z) p2 = 0 ;
  polCheck (p1) ;
  polCheck (p2) ;
  if (p1 && p2)
    {
      p1 = polCopy (p1, h) ;
      p2 = polCopy (p2, h) ;
      return polDoProduct (p1, p2, h) ;
    }
  return 0 ;
} /* polProduct */

/*************************************************************************************************/

BOOL polScale (POLYNOME pp, double complex z)
{
  if (! pp)
    return 0 ;

  polCheck (pp) ;

  if (pp->isSum)
    {
      polScale (pp->p1, z) ;
      polScale (pp->p2, z) ;
    }
  else if (pp->isProduct)
    {
      
      polScale (pp->p1, z) ;
    }
  else
    pp->tt.z *= z ;
  
  return TRUE ;
} /* polScale */

/*************************************************************************************************/

POLYNOME polMultiSum (AC_HANDLE h, POLYNOME ppp[])
{
  POLYNOME pp = 0, p1, p2 ;
  int i = -1 ;

  while (i++, ppp[i])
    ppp[i] = polCopy(ppp[i], h) ;
  if (i > 0)
    {
      pp = ppp[--i] ; 
      while (i > 0) 
	{
	  p2 = pp ;
	  p1 = ppp[--i] ;
	  pp = polDoSum (p1, p2, h) ;
	}
    }
  return pp ;
} /* polMultiSum */

/*************************************************************************************************/

POLYNOME polMultiProduct (AC_HANDLE h, POLYNOME ppp[])
{
  POLYNOME pp, p1, p2 ;
  int i = 0 ;

  while (ppp[i]) 
    {
      polCheck (ppp[i]) ;
      i++ ;
    }
  if (i <= 1) return polCopy (ppp[0], h) ;

  pp = ppp[--i] ;
  while (i > 0)
    {
      p2 = pp ;
      p1 = ppp[--i] ;
      pp = polProduct (p1, p2, h) ;
    }
  return pp ;
} /* polMultiProduct */

/*************************************************************************************************/
/******************************* t'hooft integral **************************************************/
static POLYNOME dimIntegralDo (POLYNOME pp, int pass, AC_HANDLE h) ;

/* incomplete, this only works for pairs of sigma, we need the cases4,6,8 ...
 * which create polynomes in gg, not monomes 
 */
static POLYNOME pauliTraceTT (POLYNOME pp)
{
  AC_HANDLE h = pp->h ;
  BOOL epsilon = TRUE ;
  TT tt = pp->tt ; 
  int iss ;
  short *s = tt.sigma ; 
  short *sb = tt.sigB ;
  int parity = 1 ; 

  polCheck (pp) ;
  pp->isFlat = FALSE ;

  if (s[0] && sb[0])
    messcrash ("FATAL ERROR: Computing the trace of a monome where sigma=%s and sigB=% are both present\n", s, sb) ; 
  if (sb[0])
    { s = sb ; parity = -1 ; }
  iss = slen (s) ;
  if (iss % 2)
    { tt.z = 0 ; return 0 ; }
  if (iss == 0)
    {
      if (tt.Id2)
	{
	  tt.Id2 = 0 ;
	  pp->tt.z *= 2 ; /* trace (identity) = 2 */
	}
    }
  else if (iss == 2)
    { 
      short *g = tt.g ;
      while (*g) g++ ;
      while ((*g++ = *s++)) ;
      memset (tt.sigma , 0, SMAX) ;
      memset (tt.sigB , 0, SMAX) ;
      pp->tt = tt ;

      tt.Id2 = 0 ;
      pp->tt.z *= 2 ; /* trace (identity) = 2 */
    }
  else if (iss == 4)
    { 
      int i, n, N = 4, NN = 3 ;
      short S[N] ;

      memcpy (S, s, N*sizeof (short)) ;
      pp->tt.z *= 2 ; /* trace (identity) = 2 */
      pp->tt.Id2 = 0 ;
      short *gg = tt.g ;
      int k = slen (gg) ;
      int ek = slen (tt.eps) ;
      POLYNOME ppp[NN+2] ;
      char *z[3] = { "abcd", "acbd", "adbc"} ;
      for (n = 0 ; n < NN ; n++)
	{               /* we need N products of type g_ab g_cd g_ef, then we zero terminate the list */
	  int i ;
	  ppp[n] = polCopy (pp, h) ;
	  memset (ppp[n]->tt.sigma , 0, SMAX) ;
	  memset (ppp[n]->tt.sigB , 0, SMAX) ;
	  if (n%2) ppp[n]->tt.z *= -1 ;   /* alternate signs */
	  for (i = 0 ; i < N ; i++)
	    {
	      ppp[n]->tt.g[k+i] = S[z[n][i] - 'a'] ;
	    }
	}
      if (epsilon)
	{
	  ppp[n] = polCopy (pp,h ) ;
	  memset (ppp[n]->tt.sigma , 0, SMAX) ;
	  memset (ppp[n]->tt.sigB , 0, SMAX) ;
	  for (i = 0 ; i < N ; i++)
	    {
	      ppp[n]->tt.eps[ek+i] = S[z[0][i] - 'a'] ;
	    }
	  ppp[n]->tt.eps[ek+i] = 0 ;
	  ppp[n]->tt.z *= parity ;
	  n++ ;
	}
      ppp[n++] = 0 ; /* zero terminate the list */	
      pp = polMultiSum (h, ppp) ;
    }
  else /* iss even and > 4 */
    { 
      int ii, i, j, k, N = 0 ;
      int NN = iss*iss*iss ;   /* max number of terms */
      POLYNOME ppp[NN] ;

      short *gg = tt.g ;
      int ng = slen (gg) ;
      int neps = slen (tt.eps) ;

      /* eliminate the zeroth and the ith pauli matrix and create a g term */
      for (N = 0, ii = 1 ; ii < iss ; ii++)
	{
	  int m ;
	  POLYNOME p1 = newScalar (tt.z, h) ;
	  short *s1 ;
	  p1->tt = tt ;
	  p1->tt.g[ng] = s[0] ;
	  p1->tt.g[ng+1] = s[ii] ;

	  if (ii % 2 == 0)
	    { p1->tt.z *= -1 ; }
	  p1->tt.sigma[0] = p1->tt.sigB[0] = 0 ;
	  s1 = (parity == 1 ? p1->tt.sigma : p1->tt.sigB) ;
	  for (m = 0 ; m < ii-1 ; m++)
	    s1[m] = s[m+1] ;
	  for (m = ii - 1 ; m < iss ; m++)
	    s1[m] = s[m+2] ;
	  p1 = pauliTraceTT (p1) ;
	  ppp[N++] = p1 ; 
	}
      /* eliminate the zeroth and three other pauli matrix and create an epsilon term */
      if (epsilon)
	for (i = 1 ; i < iss ; i++)
	  for (j = i + 1 ; j < iss ; j++)
	    for (k = j + 1 ; k < iss ; k++)
	      {
		int m ;
		POLYNOME p1 = newScalar (1, h) ;
		short *s1 ;
		
		p1->tt = tt ;
		p1->tt.z *= parity ;
		if ((i+k+j) % 2 == 1)
		  p1->tt.z *= -1 ;
		p1->tt.eps[neps+0] = s[0] ;
		p1->tt.eps[neps+1] = s[i] ;
		p1->tt.eps[neps+2] = s[j] ;
		p1->tt.eps[neps+3] = s[k] ;
		
		p1->tt.sigma[0] = p1->tt.sigB[0] = 0 ;
		s1 = (parity == 1 ? p1->tt.sigma : p1->tt.sigB) ;
		for (m = 0 ; m < i-1 ; m++)
		  s1[m] = s[m+1] ;
		for (m = i - 1 ; m < j - 2 ; m++)
		  s1[m] = s[m+2] ;
		for (m = j - 2 ; m < k - 3 ; m++)
		  s1[m] = s[m+3] ;
		for (m = k - 3 ; m < iss ; m++)
		  s1[m] = s[m+4] ;
		p1 = pauliTraceTT (p1) ;
		ppp[N++] = p1 ; 
	      }
      ppp[N] = 0 ;	
      /* add up all the contractions, since we always ue index zero, we are not overcounting */
      if (N >= NN)
	messcrash ("Too many terms iss=%d NN = %d N=%d", iss, N, NN) ;
      pp = polMultiSum (h, ppp) ;
    }
  return pp ;
} /* pauliTraceTT */

/*******************************************************************************************/

POLYNOME pauliTrace (POLYNOME pp, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  POLYNOME p1, p2 ;
  static int level = 0 ;
  if (!pp)
    return 0 ;
  polCheck (pp) ;

  if (level == 0)
    {
      pp = polCopy (pp, h) ;
      pp = expand (pp) ;
    }
  level++ ;

  if (pp->isSum)
    {
      p1 = pp->p1 ;
      p2 = pp->p2 ;
      if (p1) p1 = pp->p1 = pauliTrace (p1, h) ;
      if (p2) p2 = pp->p2 = pauliTrace (p2, h) ;
    }

  if (pp->tt.type)
    {
      TT tt = pp->tt ;
      short *s = tt.sigma ; 
      short *sb = tt.sigB ; 
      if (s[0] && sb[0]) messcrash ("Cannot have sigma=%s and sigmaBar=%s in the same monome", s, sb) ;
      if (tt.Id2)
	{
	  pp->isFlat = FALSE ;
	  pp = contractTtIndices (pp) ;
	  pp = pauliTraceTT (pp) ;
	}
      else if (tt.sigma[0] || tt.sigB[0])
	{
	  messcrash ("FATAL ERROR: Computing the trace of a monome with sigma, but zero Id2  sigma=%s and sigB=% are both present\n"
		     , pp->tt.sigma, pp->tt.sigB
		     ) ;
	}
      if (pp && pp->tt.type && cabs (pp->tt.z) < minAbs)
	{ pp = 0 ; }
    }

  level-- ;
  if (pp && level >= 0)
    pp = expand (pp) ;
  pp = polCopy (pp, h0) ;
  ac_free (h) ;
  return pp ;
} /* pauliTrace */
  
/*************************************************************************************************/

static POLYNOME killMomenta (POLYNOME pp)
{
  polCheck (pp) ;
  if (pp->p1) pp->p1 = killMomenta (pp->p1) ;
  if (pp->p2) pp->p2 = killMomenta (pp->p2) ;
  if (pp->isProduct)
    {
      if (! pp->p1 || ! pp->p2)
	return 0 ;
    }
  else if (pp->isSum)
    {
      if (! pp->p1)
	{ pp->p1 = pp->p2 ; pp->p2 = 0 ;}
      if (! pp->p2)
	pp = pp->p1 ;
    }
  else if (! pp->tt.type)
    return 0 ;
  else 
    {
      int i ;
      for (i = 1 ; i < 4 ; i++) /* do not kill k , just p,p,r */
	{
	  if (pp->tt.mm[i][0])
	    return 0 ;
	  pp->tt.denom[0] += pp->tt.denom[i] ; 
	  pp->tt.denom[i] = 0 ;
	}
      if (cabs (pp->tt.z) < minAbs)
	return 0 ;
    }
  if(pp) pp->isFlat = FALSE ;
  return pp ;
}

/*************************************************************************************************/
/* compute the dervative of (1/(k + p + q)^2*order) relative to (pqr)_mu */ 
static POLYNOME newDeriveDenom (POLYNOME p0, int pqr, int mu, AC_HANDLE h)
{
  TT tt = p0->tt ;
  POLYNOME ppp[6], pp = newPolynome (h) ;
  int j, nn = 0 ;
  BOOL debug = FALSE ;

  polCheck (p0) ;
  memset (ppp, 0, sizeof(ppp)) ;
  
  for (j = pqr ; j < 4 ; j++)
    {
      if (tt.denom[j]) /* ( k+p)^2 */
	{
	  int i, k, n ;
	  POLYNOME w, vv [5] ;
	  memset (vv, 0, sizeof(vv)) ;
	  if (j >=0) vv[0] = newK (mu, h) ;
	  if (j >=1) vv[1] = newP (mu, h) ;
	  if (j >=2) vv[2] = newQ (mu, h) ;
	  if (j >=3) vv[3] = newR (mu, h) ;
	  vv[4] = 0 ;
	  for (k = 0 ; k < 4 ; k++)
	    if (vv[k])
	      {
		for (i=0 ; i < 4 ; i++)
		  vv[k]->tt.denom[i] = tt.denom[i] ; /* passive factors in the denom */ 
		n = tt.denom[j] ;   /* the factor we derived */
		vv[k]->tt.denom[j] = n + 1 ;
		vv[k]->tt.z *= (-2 * n ) ;
	      }
	  w = polMultiSum (h,vv) ; /* 1/(k+p+q)^n  -> w = (k+p+q) = derivee du denominateur */
	  ppp[nn++] = w ;
	}
    }
  pp = 0 ;
  if (nn)
    {
      pp = polMultiSum (h, ppp) ;
      pp->isFlat = FALSE ;
      if (debug) checkPolynome (pp) ;  
    }
  return pp ;  
} /* newDeriveDenom */

/*************************************************************************************************/
/* partial derivative of a polynome with respect to p_mu or q_mu or r_mu */
static POLYNOME derivePdo (POLYNOME pp, int pqr, int mu, AC_HANDLE h)
{
  static     POLYNOME empty = 0 ;
  TT tt ;
  BOOL hasDenom = FALSE ;
  BOOL debug = FALSE ;

  if (! pp) return 0 ;
  polCheck (pp) ;
  if (debug) checkPolynome (pp) ;  
  if (pp->isSum)
    { /* linearity */
      pp->p1 = derivePdo (pp->p1, pqr, mu,h) ;
      pp->p2 = derivePdo (pp->p2, pqr, mu,h) ;
      pp->isFlat = FALSE ;
      return pp ;
    }
  if (pp->isProduct)
    {
      POLYNOME q1 = polCopy (pp->p1, h) ;
      POLYNOME q2 = polCopy (pp->p2, h) ;

      q1 = derivePdo (q1, pqr, mu,h) ;
      q2 = derivePdo (q2, pqr, mu,h) ;
      
      if (!q1 && ! q2)
	pp = 0 ; 
      else if (q1 && ! q2)
	pp->p1 = q1 ; 
      else if (! q1 && q2)
	pp->p2 = q2 ;
      else
	{ /* Leibnitz */
	  POLYNOME r1 = polProduct (q1, pp->p2, h) ;
	  POLYNOME r2 = polProduct (pp->p1, q2, h ) ;
	  pp->isSum = TRUE ;
	  pp->isProduct = FALSE ;
	  pp->p1 = r1 ; 
	  pp->p2 = r2 ;
	}
      if (pp) pp->isFlat = FALSE ;
      if (debug) checkPolynome (pp) ;  
      return pp ;
    }

  /* free object */
  if (pqr <1 || pqr > 3)
      messcrash ("You can only partial derive with respect to 1:p, 2:q, 3:r, not %d\n", pqr) ;
    
  /* derive the numerator */
  if (! empty)
    empty = newPolynome (h) ;

  
  tt = pp->tt ;

  if (! hasDenom)
    {
      int i ;
      for (i = pqr ; i < 4 ; i++)
	if (tt.denom[i])
	  hasDenom = TRUE ;
    }
  if (! tt.mm[pqr][0] && ! hasDenom)
    return 0 ;
  else if (tt.mm[pqr][0] && ! hasDenom)
    {
      short *u0 = tt.mm[pqr] ;
      int i, k, iMax = slen (u0)  ;
      
      if (iMax == 1) /* simplest case, just add a g_munu and suppress the p_mu */
	{
	  pp = polCopy (pp, h) ;
	  short *v = tt.g, *w = tt.mm[pqr] ;
	  v += slen (v) ;
	  v[0] = mu ; v[1] = tt.mm[pqr][0] ; v[2] = 0 ;
	  memset (w, 0, SMAX) ;	
	  if(v[0]==v[1]) tt.z /= 4.0 ;  /* because we want delta_{aa} not g_{aa} */
	  
	  pp->tt = tt ;	  pp->isFlat = FALSE ;
	  if (debug) checkPolynome (pp) ;  
	  return pp ;
	}
      else 
	{
	  POLYNOME qq[iMax+1] ;
	  for (k = 0 ; k < iMax ; k++)
	    {
	      short *v,*w = tt.mm[pqr] ;
	      /* copy the original polynome */
	      qq[k] = newPolynome (h) ;
	      qq[k]->tt.type = 1 ;
	      qq[k]->tt = tt ;
	      /* replace one dependence on p_alpha by g_mu_alpha */
	      v = qq[k]->tt.mm[pqr] ; 
	      for (i = k ; i < iMax ; i++)
		v[i] = v[i+ 1] ; 
	      v = qq[k]->tt.g ;
	      v += slen (v) ;
	      v[0] = mu ; v[1] = u0[k] ; v[2] = 0 ;
	      memset (w, 0, SMAX) ;
	      if(v[0]==v[1]) qq[k]->tt.z /= 4.0 ;  /* because we want delat_{aa} not g_{aa} */
	    }
	  qq[iMax] = 0 ;
pp = polMultiSum (h, qq) ;
	  if (debug) checkPolynome (pp) ;  
	}
    }
  else if (hasDenom)
    { /* construct the product (f'/g + fg'/g^2) */
      POLYNOME q3, q1prime, q1 = newScalar (1,h) ;
      POLYNOME q4, q2prime, q2 = newScalar(1,h) ;
      int i ;

      /* contruct f: copy and remove the denom */
      q1->tt = tt ;
      for (i = 0 ; i < 4 ; i++)
	q1->tt.denom[i] = 0 ;
      /* contruct g: copy the denom, but certainly not the matrices */
      for (i = 0 ; i < 4 ; i++)
	q2->tt.denom[i] = tt.denom[i] ;

      /* derive */
      q1prime = derivePdo (q1, pqr, mu, h) ;
      q2prime = newDeriveDenom(q2, pqr, mu, h) ;

      /* construct the 2 products (q1 and q1prime both commute) */
      q3 = polProduct (q1, q2prime,h) ;
      q4 = q1prime ? polProduct (q2, q1prime,h) : 0 ;
      pp = q4 ? polSum (q3,q4,h) : q3 ;
      if (debug) checkPolynome (pp) ;  
      /*
      ac_free (q1) ;      
      ac_free (q2) ;
      ac_free (q3) ;
      ac_free (q4) ;
      ac_free (q5) ;
      ac_free (q1prime) ;
      ac_free (q2prime) ;
      */
    }  
  if (debug) checkPolynome (pp) ;  
  if (pp) pp->isFlat = FALSE ;
  return pp ;
} /* derivePdo */

/*************************************************************************************************/

static POLYNOME deriveP (POLYNOME pp, int pqr, int mu, AC_HANDLE h)
{
  BOOL debug = FALSE ;

  polCheck (pp) ;
  pp = derivePdo (pp, pqr, mu,h) ;
  if (0 && pp) pp = killMomenta (pp) ;
  pp = expand (pp) ;
  if (debug) showPol (pp) ;
  pp = contractProducts (pp) ;
  if (debug) showPol (pp) ;
  pp = contractIndices (pp) ;
  if (debug) showPol (pp) ;
  if (pp) pp->isFlat = FALSE ;
  return pp ;
} /* deriveP */

/*************************************************************************************************/
/*************************************************************************************************/

static POLYNOME dimIntegrateByPart (POLYNOME pp) ;
/* dimensional integral, only works on a flat expression of total degree 4 */
static POLYNOME dimIntegralMonome (POLYNOME pp, int state, short *kk, int *np, int pass)
{
  int i, j, k ;
  BOOL debug = FALSE ;

  if (! pp) return 0 ;
  polCheck (pp) ;
  AC_HANDLE h = pp->h ;
		  
  if (debug) checkPolynome (pp) ;  
  if (debug)
    showPol(pp) ;
  switch (state)
    {
    case 0:
      if (1)
	{
	  int nd = 0 ;
	  short buf[GMAX] ;
	  
	  memset (buf, 0, sizeof(buf)) ;
	  dimIntegralMonome (pp, 1, buf, &nd, pass) ;
	  if (debug) checkPolynome (pp) ;  
	  k = 4 + slen (buf) - 2 * nd ;
	  if (k < 0) /* convergent integral */
	    return 0 ;
	  if (! pass)
	    return pp ;
	  else if (k == 1)
	    return dimIntegrateByPart (pp) ;
	  else if (k == 2)
	    {
	      pp = dimIntegrateByPart (pp) ;
	      if (! pp)
		return 0 ;
	      if (pp->isProduct)
		pp->p1->tt.z /= 2 ;
	      else if (pp->isSum)
		{
		  POLYNOME qq = newScalar (.5,h) ;
		  pp = polProduct (qq, pp,h) ;
		}
	      else
		pp->tt.z /= 2 ;
	      return pp ;
	    }
	  else if (k)
	    messcrash ("This integral should have been differentiated, it has k-order %d > 0", k) ;

	  dimIntegralMonome (pp, 2, 0, 0, pass) ; /* clean up all denoms and all k dependencies */
	  if (debug) checkPolynome (pp) ;  
	  /* transform the polynome into a sum of g_munu */
	  k = slen (buf) ;
	  if (k) 
	    {
	      /* eliminate repeated indices    k^2 / k^4 (k+p)^2 = 1/k^2 (k+p)^2 */
	      for (i = 0 ; i < k - 1 ; i++)
		for (j = i + 1 ; j < k ; j++)
		  if (buf[i] == buf[j])
		    buf[i] = buf[j] = 0 ;
	      for (i = j = 0 ; i < k ; i++)
		if (buf[i]) { buf[j] = buf[i] ; j++ ; }
	      buf[j] = 0 ;
	    }
	  k = slen (buf) ; /* simplified length */
	  if (k == 0)
	    {
	    }
	  else if (k == 2)
	    { 
	      int n, N = 2, NN = 1 ;
	      POLYNOME qqq[NN+1] ;
	      short *cp ;
	      pp->tt.z /= 4 ;
	      for (n = 0 ; n < NN ; n++)
		{                             /* we need N products of type g_ab g_cd g_ef, then we zero terminate the list */
		  int i ;
		  qqq[n] = polCopy (pp, h) ;
		  cp = qqq[n]->tt.g ;
		  while (*cp) cp++ ;
		  
		  for (i = 0 ; i < N ; i++)
		    *cp++ = buf[i] ;
		  *cp++ = 0 ;
		}
	      qqq[n] = 0 ; /* zero terminate the list */	
	      pp = polMultiSum (h, qqq) ;
	    }	  
	  else if (k == 4)
	    { 
	      int n, N = 4, NN = 3 ;
	      POLYNOME qqq[NN+1] ;
	      short *cp ;
	      char *z[3] = {
		"abcd","acbd","adbc"
	      } ;
	      pp->tt.z /= 24 ;
	      for (n = 0 ; n < NN ; n++)
		{                             /* we need N products of type g_ab g_cd g_ef, then we zero terminate the list */
		  int i ;
		  qqq[n] = polCopy (pp, h) ;
		  cp = qqq[n]->tt.g ;
		  while (*cp) cp++ ;

		  for (i = 0 ; i < N ; i++)
		    *cp++ = buf[z[n][i] - 'a'] ;
		  *cp++ = 0 ;
		}
	      qqq[n] = 0 ; /* zero terminate the list */	
	      pp = polMultiSum (h, qqq) ;
	    }	  
	  else if (k == 6)
	    { 
	      int n, N = 6, NN = 15 ;
	      POLYNOME qqq[NN+1] ;
	      short *cp ;
	      char *z[15] = {
		"abcdef","abcedf","abcfde",
		"acbdef","acbedf","acbfde",
		"adbcef","adbecf","adbfce",
		"aebcdf","aebdcf","aebfcd",
		"afbcde","afbdce","afbecd"
	      } ;
	      pp->tt.z /= 192 ; 
	      for (n = 0 ; n < NN ; n++)
		{                             /* we need N products of type g_ab g_cd g_ef, then we zero terminate the list */
		  int i ;
		  qqq[n] = polCopy (pp, h) ;
		  cp = qqq[n]->tt.g ;
		  while (*cp) cp++ ;

		  for (i = 0 ; i < N ; i++)
		    *cp++ = buf[z[n][i] - 'a'] ;
		  *cp++ = 0 ;
		}
	      qqq[n] = 0 ; /* zero terminate the list */	
	      pp = polMultiSum (h, qqq) ;
	    }	  
	  else if (k == 8)
	    { 
	      int n, N = 8, NN = 105 ;
	      POLYNOME qqq[NN+1] ;
	      short *cp ;
	      char *z[105] = {
		"abcdefgh","abcdegfh","abcdehfg",
		"abcedfgh","abcedgfh","abcedhfg",
		"abcfdegh","abcfdgeh","abcfdheg",
		"abcgdefh","abcgdfeh","abcgdhef",
		"abchdefg","abchdfeg","abchdgef",
		
		"acbdefgh","acbdegfh","acbdehfg",
		"acbedfgh","acbedgfh","acbedhfg",
		"acbfdegh","acbfdgeh","acbfdheg",
		"acbgdefh","acbgdfeh","acbgdhef",
		"acbhdefg","acbhdfeg","acbhdgef",
		
		"adcbefgh","adcbegfh","adcbehfg",
		"adcebfgh","adcebgfh","adcebhfg",
		"adcfbegh","adcfbgeh","adcfbheg",
		"adcgbefh","adcgbfeh","adcgbhef",
		"adchbefg","adchbfeg","adchbgef",
		
		"aecdbfgh","aecdbgfh","aecdbhfg",
		"aecbdfgh","aecbdgfh","aecbdhfg",
		"aecfdbgh","aecfdgbh","aecfdhbg",
		"aecgdbfh","aecgdfbh","aecgdhbf",
		"aechdbfg","aechdfbg","aechdgbf",
		
		"afcdebgh","afcdegbh","afcdehbg",
		"afcedbgh","afcedgbh","afcedhbg",
		"afcbdegh","afcbdgeh","afcbdheg",
		"afcgdebh","afcgdbeh","afcgdheb",
		"afchdebg","afchdbeg","afchdgeb",
		
		"agcdefbh","agcdebfh","agcdehfb",
		"agcedfbh","agcedbfh","agcedhfb",
		"agcfdebh","agcfdbeh","agcfdheb",
		"agcbdefh","agcbdfeh","agcbdhef",
		"agchdefb","agchdfeb","agchdbef",
		
		"ahcdefgb","ahcdegfb","ahcdebfg",
		"ahcedfgb","ahcedgfb","ahcedbfg",
		"ahcfdegb","ahcfdgeb","ahcfdbeg",
		"ahcgdefb","ahcgdfeb","ahcgdbef",
		"ahcbdefg","ahcbdfeg","ahcbdgef"
	      } ;
	      if (4) pp->tt.z /= 1920 ;
	      for (n = 0 ; n < NN ; n++)
		{                             /* we need N products of type g_ab g_cd g_ef, then we zero terminate the list */
		  int i ;
		  qqq[n] = polCopy (pp, h) ;
		  cp = qqq[n]->tt.g ;
		  while (*cp) cp++ ;

		  for (i = 0 ; i < N ; i++)
		    *cp++ = buf[z[n][i] - 'a'] ;
		  *cp++ = 0 ;
		}
	      qqq[n] = 0 ; /* zero terminate the list */	
	      pp = polMultiSum (h, qqq) ;
	    }	  
	  else 
	    messcrash ("Sorry, i cannot yet integrate k^10/k^14\n") ;
	}
      if (debug) checkPolynome (pp) ;  
      break ;

    case 1:
      if (pp->isSum)
	messcrash ("dimIntegralMonome found a sum hidden under a product, the polynome is not flat") ;
      if (pp->isProduct)
	{
	  dimIntegralMonome (pp->p1, 1, kk, np, pass) ;
	  dimIntegralMonome (pp->p2, 1, kk, np, pass) ;
	}
      if (pp->tt.type)
	{
	  int i ;
	  short *u, *v ;
	  
	  for (i = 0 ; i < 4 ; i++)
	    *np += pp->tt.denom[i] ;  /* power of k in 1/(k+...)^2 */
	  u = kk ; u += slen(u) ; 
	  v = pp->tt.mm[0] ;        /* numerator k_mu k_nu ... */
	  while ((*u++ = *v++)) ;
	  if (kk[GMAX-1])
	    messcrash ("Overflow while collecting the k[] vectors ") ;
	}
      break ;
    case 2:  /* clean up the k dependency */
      if (pp->isProduct)
	{
	  dimIntegralMonome (pp->p1, 2, 0, 0, pass) ;
	  dimIntegralMonome (pp->p2, 2, 0, 0, pass) ;
	}
      if (pp->tt.type)
	{
	  for (i = 0 ; i < 4 ; i++)
	    pp->tt.denom[i] = 0 ;   /* power of k in 1/(k+...)^2 */
	  memset (pp->tt.mm[0], 0, SMAX) ; /* numerator k_mu k_nu ... */
	}
      break ;
    }
  if (pp) pp->isFlat = FALSE ;
  if (debug) checkPolynome (pp) ;
  if (debug)
    showPol(pp) ;
  return pp ;
}  /* dimIntegralMonome */
  
/*************************************************************************************************/

static POLYNOME dimIntegrateByPart (POLYNOME pp) 
{
  AC_HANDLE h = pp->h ;
  POLYNOME ppp[4] ;
  int nn = 0 ;
  int i, j ;
  short mu = newDummyIndex () ;
  static int level = 0 ;
  BOOL debug = FALSE ;

  polCheck (pp) ;
  level++ ;
  for (j = 1 ; j <= 3 ; j++)
    {
      POLYNOME p1d = 0, p1 = polCopy (pp,h) ;
      POLYNOME p2d = 0, p2 = newScalar (1,h) ;
      POLYNOME puu = 0, pvv = 0 ;
 
     if (debug)
       {
	 int i ;
	 for (i = 0 ; i < level ; i++)
	   printf ("###") ;
	 printf ("# level %d before\n", level) ;
	 showPol (pp) ;
       }
      for (i = 0 ; i < 4 ; i++)
	{
	  p2->tt.denom[i] = p1->tt.denom[i] ;
	  p1->tt.denom[i] = 0 ; 
	}
      for (i = 0 ; i < GMAX ; i++)
	{
	  p2->tt.mm[0][i] = p1->tt.mm[0][i] ;
	  p1->tt.mm[0][i] = 0 ;
	}

      if (0) p1d = deriveP (p1, j, mu,h) ;
      p2d = deriveP (p2, j, mu,h) ;
      if (p1d)
	{
	  POLYNOME pm = newScalar (1.0, h) ;
	  pm->tt.mm[j][0] = mu ;
	  puu = dimIntegralDo (p1d, 1, h) ;
	  if (puu) puu = polProduct (pm, puu,h) ;
	}
      if (p2d)
	{
	  POLYNOME pm = newScalar (1.0, h) ;
	  pm->tt.mm[j][0] = mu ;
	  pvv = dimIntegralDo (p2d, 1, h) ;
	  if (pvv) pvv = polProduct (pm, pvv, h) ;
	}
      if (puu) puu = polProduct (p2, puu, h) ;
      if (pvv) pvv = polProduct (p1, pvv,h) ;
      if (puu && pvv) 
	ppp[nn++] = polSum (puu, pvv, h) ;
      else if (puu)
	ppp[nn++] = puu ;
      else if (pvv)
	ppp[nn++] = pvv ;
     if (debug)
       {
	 int i ;
	 for (i = 0 ; i < level ; i++)
	   printf ("###") ;
	 printf ("# level %d after\n", level) ;
	 showPol (pvv) ;
       }
    }
  ppp[nn] = 0 ;

  if (nn > 1)
    pp = polMultiSum (h, ppp) ;
  else if (nn == 1)
    pp = ppp[0] ;
  else
    pp = 0 ;

  level-- ;
  return pp ;
}

/*************************************************************************************************/
/* dimensional integral, only works on a flat expression of total degree 4 */
static POLYNOME dimIntegralDo (POLYNOME pp, int pass, AC_HANDLE h)
{  
  BOOL debug = FALSE ;
  static int level = 0 ;
  static int nnn = 0 ;
  if (! pp) return 0 ;
  level++ ;
  if (level == -1) firstDummyIndex += 4 ;
  pp = polCopy (pp, h) ;
  pp = expand (pp) ; 
  pp = contractProducts (pp) ;
  contractIndices (pp) ;
  if (! pp) { level-- ; return 0 ; }
 
  if (debug && level <=2) 
    {
      int i = level ;
      while (i--)
	printf ("##") ;
      showPol(pp) ;
    }
  if (pp->isSum)
    { /* linearity */
      int ii = 0, jj = 0 ;
      Array aa = arrayCreate (32, POLYNOME) ;
      
      sortPolGetSum (pp, aa) ;
      for (ii = jj = 0 ; ii < arrayMax (aa) ; ii++)
	{
      	  POLYNOME qq = arr (aa, ii, POLYNOME) ; 
	  nnn++ ;
	  qq = dimIntegralDo (qq, pass, h) ;
	  if (qq)
	    arr (aa, jj++, POLYNOME) = qq ;
	}
      arrayMax (aa) = jj ;
      pp = sortReduceSum (aa, h) ;
   }
  else
    {
      pp = dimIntegralMonome (pp, 0, 0, 0, pass) ;
    }
  if (pp) 
    {
      pp->isFlat = FALSE ;
      pp = expand (pp) ;
    }
  if (debug && level <=2) 
    {
      int i = level ;
      while (i--)
	printf ("XX") ;
      showPol(pp) ;
    }
  level-- ;
  return pp ;
}   /* dimIntegralDo */



/*************************************************************************************************/
/*************************************************************************************************/

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME squareMomentaCleanUpDo (POLYNOME pp, short alpha) 
{
  pp = expand (pp) ;

  if (! pp)
    return 0 ;

  if (pp->isSum && pp->p1)
    pp->p1 = squareMomentaCleanUpDo (pp->p1, alpha) ;
  if (pp->isSum && pp->p2)
    pp->p2 = squareMomentaCleanUpDo (pp->p2, alpha) ;
  if (pp->tt.type == 1)
    {
      if (pp->tt.mm[1][0] && pp->tt.mm[1][0] == pp->tt.mm[1][1])
	pp->tt.mm[1][0] = pp->tt.mm[1][1] = alpha ;
    }
  pp = expand (pp) ;

  return pp ;
}

POLYNOME squareMomentaCleanUp (POLYNOME pp)
{
  short alpha = newDummyIndex() ;
  return squareMomentaCleanUpDo (pp, alpha) ;
}

/***********************************************************************************************************************************************/
#ifdef JUNK
static POLYNOME squareMomentaKill (POLYNOME pp)
{
  if (! pp)
    return 0 ;

  if (pp->isSum && pp->p1)
    pp->p1 = squareMomentaKill (pp->p1) ;
  if (pp->isSum && pp->p2)
    pp->p2 = squareMomentaKill (pp->p2) ;
  if (pp->tt.type == 1)
    {
      if (pp->tt.mm[1][0] && pp->tt.mm[1][0] == pp->tt.mm[1][1])
	pp->tt.z = 0 ;
    }
  pp = expand (pp) ;

  return pp ;
}
#endif
/***********************************************************************************************************************************************/

POLYNOME momentaCleanUp (POLYNOME pp, short alpha) 
{
  AC_HANDLE h = pp->h ;
  int nn = 0 ;
  POLYNOME ppp[4] ;
  POLYNOME p1, p2, p3, p4, p5 ;
  polCheck (pp) ;
  pp = expand (pp) ;
  p1 = polCopy (pp, h) ;
  p2 = deriveP (p1, 1, alpha, h) ;
  p3 = expand (p2) ;
  p4 = contractIndices (p3) ;
  p5 = expand (p4) ;
  if (p5)
    ppp[nn++] = polProduct (newP(alpha,h), p5, h) ;

  p1 = polCopy (pp, h) ;
  p2 = deriveP (p1, 2, alpha, h) ;
  p3 = expand (p2) ;
  p4 = contractIndices (p3) ;
  p5 = expand (p4) ;
  if (p5)
    ppp[nn++] = polProduct (newQ(alpha,h), p5, h) ;
  p1 = polCopy (pp, h) ;
  p2 = deriveP (p1, 3, alpha, h) ;
  p3 = expand (p2) ;
  p4 = contractIndices (p3) ;
  p5 = expand (p4) ;
  if (p5)
    ppp[nn++] = polProduct (newR(alpha,h), p5, h) ;
  ppp[nn++] = 0 ;

  return polMultiSum (h, ppp) ;
}
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

POLYNOME dimIntegral (POLYNOME p0)
{
  AC_HANDLE h ;
  POLYNOME pp ;
  BOOL debug = FALSE ;
  if (! p0) return 0 ;
  h = p0->h ;
  polCheck (p0) ;
  pp = polCopy (p0, h) ;
  pp = expand (pp) ;
    if (debug) showPol(pp) ;
  freeIndex (pp) ;
    if (debug) showPol(pp) ;
  pp = contractProducts (pp) ;
  if (debug) showPol(pp) ;
  freeIndex (pp) ;
    if (debug) showPol(pp) ;
    contractIndices (pp) ;
    if (debug) showPol(pp) ;
    freeIndex (pp) ;
  if (debug) showPol(pp) ;
  pp = dimIntegralDo (pp, 0, h) ;
    if (debug) showPol(pp) ;
    freeIndex (pp) ;
    
  if (debug) showPol(pp) ;
  pp = dimIntegralDo (pp, 1, h) ;
    if (debug) showPol(pp) ;
    freeIndex (pp) ;

  if (debug) showPol(pp) ;
  
  return pp ;
} /* dimIntegral */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

/***********************************************************************************************************************************************/
/* given a list of tensor indices, fuse and eliminate repated index */
static BOOL freeIndexFuse (short *top, short **cpp)
{
  int ii, k, kMax = 0 ;
  short *cp, nn[INDEXMAX] ;
  memset (nn, 0, sizeof (nn)) ;

  for (ii = 0, cp = cpp[ii] ; cp ; cp = cpp[++ii])
    while (*cp)
      {
	int k = *cp++ ;
	nn[k]++ ;
	if (nn[k] > 2)
	  messcrash ("Index %c is repeated %d times in this tensor", k, nn[k]) ;
	else if (nn[k] == 1)
	  { if (k > kMax)
	      kMax = k ;
	  }
      }
  for (ii = 0, k = 1  ; k <= kMax ; k++)
    if (nn[k] == 1)
      top[ii++] = k ;
  top[ii] = 0 ;
  if (ii >= GMAX)
    messcrash ("too many indicies in this tensor %s", top) ;
  return TRUE ;

} /* freeIndexFuse */

/***********************************************************************************************************************************************/

static POLYNOME bbCleanUpDo (POLYNOME pp, short a, short b, short c, short d, BOOL *okp)
{
  if (!pp)
    return 0 ;
  
  polCheck (pp) ;

  if (pp->tt.type && cabs (pp->tt.z) < minAbs)
    { return 0 ; }
  
  if (pp->isSum)
    {
      pp->p1 = bbCleanUpDo (pp->p1, a, b, c, d, okp) ;
      pp->p2 = bbCleanUpDo (pp->p2, a, b, c, d, okp) ;
    }
  else if (pp->tt.type)
    {
      TT tt = pp->tt ;
      if (tt.eps[0])
	{
	  short *s = tt.eps ;
	  int n = slen (s) ;
	  int pass = 0 ;
	  int i, j, kk, k1, k2 ;

	  for (pass = 0 ; pass < 2 ; pass++)
	    {
	      short A = a, B = b ;
	      if (pass == 1)
		{ A = c ; B = d ; }

	      for (i = 0 ; i < n ; i+= 4)
		{
		  for (j = kk = k1 = k2 = 0 ; j < 4 ; j++)
		    {
		      if (s[i+j] == A)
			{ k1 = j + 1 ; kk++ ; }
		      if (s[i+j] == B)
			{ k2 = j + 1 ; kk++ ; }
		    }
		  if (kk == 2) /* i can replace this epsilon by a (gg - gg) */
		    {
		      short e = 0, f = 0 ;
		      
		      for (j = 0 ; j < 4 ; j++)
			{
			  if (j + 1 != k1 && j + 1 != k2)
			    {
			      if (e == 0)
				e = s[i+j] ;
			      else
				f = s[i+j] ;
			    }
			  if (e && f) /* replace eps(a,b,e,f) by 2g (aebf) */
			    {
			      short *cp ;
			      tt.z *= 2 ;
			      if (pass == 1) /* flip the epsilon sign */
				tt.z *= -1 ;
			      if ((k1 + k2) % 2 == 0)  /* non contiguous ab indices inside the epsilon */
				tt.z *= -1 ;
			      for (j = 0 ; s[i+j+4] ; j++)
				s[i+j] = s[i+j+4] ;
			      for (; i+j < GMAX ; j++)
				s[i+j] = 0 ;
			      cp = tt.g ;
			      while (*cp) cp++ ;
			      *cp++ = A ;
			      *cp++ = e ;
			      *cp++ = B ;	      
			      *cp++ = f ;
			      *cp++ = 0 ;
			      
			      *okp = FALSE ;
			      pp->tt = tt ;
			      return pp ;
			    }
			}
		    }
		}		
	    }
	}

      if (tt.g[0])
	{
	  short *s = tt.g ;
	  int n = slen (s) ;
	  int pass = 0 ;
	  int i1, j1, i2, j2 ;

	  if (n == 4)
	    {
	      short A = a, B = b ;
	      short C = c, D = d ;
	      int k = 0 ;
	      for (i1 = 0 ; i1 < n ; i1 ++)
		{
		  short x = s[i1] ;
		  if (x == A || x == B || x == C || x == D)
		    k++ ;
		}
	      if (k == 4)
		return 0 ;
	    }
	  for (pass = 0 ; pass < 2 ; pass++)
	    {
	      short A = a, B = b ;
	      short e = 0, f = 0 ;
	      if (pass == 1)
		{ A = c ; B = d ; }

	      for (i1 = 0 ; i1 < n ; i1 += 2)
		{
		  for (j1 = 0 ; j1 < 2 ; j1++)
		    if (s[i1+j1] == A)
		      {
			e = s [i1+1-j1] ;
			if (e == B)
			  return 0 ;
			for (i2 = 0 ; i2 < n ; i2 += 2)
			  {
			    for (j2 = 0 ; j2 < 2 ; j2++)
			      if (s[i2+j2] == B)
				{
				  f = s [i2+1-j2] ;
				  if (f == A)
				    return 0 ;
				  /* bingo if e>f we switch a and b */
				  if (e > f)
				    {
				      tt.z *= -1 ;
				      s [i1+j1] = B ;
				      s [i2+j2] = A ;
				      *okp = FALSE ;
				      pp->tt = tt ;
				      return pp ;
				    }
				}
			  }
		      }
		}
	    }
	}
    }
  
  polCheck (pp) ;
  return pp ;
} /* bbCleanUpDo */

/*******************************************/

POLYNOME bbCleanUp (POLYNOME pp, short a, short b, short c, short d)
{
  BOOL ok = FALSE ;
  
  if (!pp)
    return 0 ;
  polCheck (pp) ;
  while (pp && ! ok)
    {
      ok = TRUE ;
      pp = expand (pp) ;
      pp = bbCleanUpDo (pp, a, b, c, d, &ok) ;
    }
  
  return pp ;
} /* bbCleanUp */

/***********************************************************************************************************************************************/

static int scmp (short *cp, const short *cq)
{
  while (*cp || *cq) 
    {
      int n = *cp++ - *cq++ ;
      if (n) return n ;
    }
  
  return 0 ;
}

static void scpy (short *cp, const short *cq)
{
  while ((*cp++ = *cq++)) ;
  return ; 
}

BOOL freeIndex (POLYNOME pp)
{
  POLYNOME p1, p2 ;
  
  if ( !pp)
    return TRUE ; /* no problem */ ;
  polCheck (pp) ;
  p1 = pp->p1 ;
  p2 = pp->p2 ;

  memset (pp->tt.freeIndex, 0, SMAX) ;
  if (pp->isProduct)
    {
      freeIndex (p1) ;
      freeIndex (p2) ;
      if (! p1 && ! p2)
	return TRUE ;
      if (p1 && ! p2)
	{
	  scpy (pp->tt.freeIndex, p1->tt.freeIndex) ;
	  return TRUE ;
	}
      if (! p1 && p2)
	{
	  scpy (pp->tt.freeIndex, p2->tt.freeIndex) ;
	  return TRUE ;
	}
      if (p1 && p2)
	{
	  short *cp = pp->tt.freeIndex ;
	  short *cp1 = p2->tt.freeIndex ;
	  short *cp2 = p1->tt.freeIndex ;
	  short *cpp[] = {cp1, cp2, 0 } ;

	  freeIndexFuse (cp, cpp) ;
	}
      return TRUE ; /* a product may have any set of free indicies */
    }
  
  else if (pp->isSum)
    {
      freeIndex (p1) ;
      freeIndex (p2) ;
      if (! p1 && ! p2)
	return TRUE ;
      if (p1 && ! p2)
	{
	  scpy (pp->tt.freeIndex, p1->tt.freeIndex) ;
	  return TRUE ;
	}
      if (! p1 && p2)
	{
	  scpy (pp->tt.freeIndex, p2->tt.freeIndex) ;
	  return TRUE ;
	}
      if (p1 && p2)
	{
	  short *cp = pp->tt.freeIndex ;
	  short *cp1 = p1->tt.freeIndex ;
	  short *cp2 = p2->tt.freeIndex ;

	  if (scmp (cp1, cp2))
	    messcrash ("Index not repeated in a sum  #%s#  #%s#", p1->tt.freeIndex, p2->tt.freeIndex) ;
	  scpy (cp, cp1) ;
	}
      return TRUE ; /* a product may have any set of free indicies */
    }
  else if (pp->tt.type)
    {
      int jj, n = 0 ;
      short *cp = pp->tt.freeIndex ;
      short *cpp[12] ;
      static int pass = 0 ;

      pass++ ;
      if (*pp->tt.g) cpp[n++] = pp->tt.g ;
      if (*pp->tt.sigma) cpp[n++] = pp->tt.sigma ;
      if (*pp->tt.sigB) cpp[n++] = pp->tt.sigB ;
      if (*pp->tt.eps) cpp[n++] = pp->tt.eps ;
      for (jj = 0 ; jj < 4 ; jj++)
	if (*pp->tt.mm[jj]) cpp[n++] = pp->tt.mm[jj] ;
      cpp[n++] = 0 ;
      if (n > 1)
	freeIndexFuse (cp, cpp) ;
    }
  
  return TRUE ;
}

/***********************************************************************************************************************************************/
/*************************************************************************************/
/***************** Polynome Matrices *************************************************/

void pmxShow (PMX pmx)
{
  int N = pmx ? pmx->N : 0 ;
  pmxCheck (pmx) ;
  printf("#### PolynomeMatrix %s\n", pmx->title) ;
  for (int i = 0 ; i < N ; i++)
    {
      char *sep = " ##" ;
      for (int j = 0 ; j < N ; j++)
	{
	  POLYNOME p = pmx->pp[N*i+ j] ;
	  if (p)
	    {
	      printf("%s i=%d j=%d\t", sep, i, j)  ;
	      sep = "  #" ;
	      showPol (p) ;
	    }
	}
    }
  printf("\n") ;
} /* pmxShow */

/*************************************************************************************/

static void pmxFinalize (void *vp)
{
  PMX pmx = (PMX) vp ;
  if (pmx && pmx->id)
    {
      pmx->magic = 0 ;
      pmx->id = 0 ;
      ac_free (pmx->h) ;
    }
} /* pmxFinalize */

/*************************************************************************************/

PMX pmxCreate (int N, char *title, AC_HANDLE h)
{
  static int id = 0 ;
  int n = sizeof (struct pmxStruct) ;
  PMX pmx = (PMX) handleAlloc (pmxFinalize, h, n) ;
  memset (pmx, 0, n) ;
  pmx->h = ac_new_handle () ;
  pmx->N = N ;
  pmx->title = title ? strnew (title, pmx->h) : "No_title" ;
  pmx->id = ++id ;
  pmx->magic = PMXMAGIC ;
  pmx->pp = halloc ((N*N+1)*sizeof(POLYNOME), h) ;
  memset (pmx->pp, 0, (N*N+1)*sizeof(POLYNOME)) ;
  return pmx ;
} /* pmxCeate () */

/*************************************************************************************/

PMX pmxCopy (PMX pmx, char *title, AC_HANDLE h) 
{
  int N = pmx->N ;
  PMX px = pmxCreate (pmx->N, title, h) ;
  pmxCheck (pmx) ;
  for (int i = 0 ; i < N ; i++)
    for (int j = 0 ; j < N ; j++)
      {
	POLYNOME p = pmx->pp[N*i + j] ;
	px->pp[N*i + j] = polCopy (p, pmx->h)  ;
     }
  return px ;  
} /* pmxCopy */

/*************************************************************************************/

BOOL pmxExpand (PMX pmx)
{
  int N = pmx->N ;
  AC_HANDLE h1 = ac_new_handle () ;
  pmxCheck (pmx) ;


  for (int i = 0 ; i < N ; i++)
    for (int j = 0 ; j < N ; j++)
      {
	POLYNOME p = pmx->pp[N*i + j] ;
	if (p)
	  {
	    POLYNOME q = expand (p) ;  /* we can kill the previous p */
	    pmx->pp[N*i + j] = polCopy (q, h1) ;
	  }
     }
  ac_free (pmx->h) ;
  pmx->h = h1 ;
  return TRUE ;
} /* pmxExpand */

/*************************************************************************************/
/* Add a zero terminated set of PMX */
PMX pmxMultiSum (PMX *pmxs, char *title, AC_HANDLE h)
{
  int N = 0, nn = 0 ;
  PMX *pmx, pmx1 ;
  for (pmx = pmxs ; *pmx ; pmx++)
    {
      PMX px = *pmx ;
      pmxCheck (px) ;
      nn++ ;
      if (!N)
	N = px->N ;
      if (N != px->N)
	messcrash ("Inhomogenous call to pmxAdd N=%d != px->N = %d in %s", N, px->N, px->title) ;
    }
  pmx1 = pmxCreate (N, title, h) ;
  
  for (int i = 0 ; i < N ; i++)
    for (int j = 0 ; j < N ; j++)
      {
	int n = 0 ;
	POLYNOME qq[nn+1] ;
	for (pmx = pmxs ; *pmx ; pmx++)
	  {
	    PMX px = *pmx ;
	    POLYNOME p = px->pp[N*i + j] ;
	    if (p)
	      qq[n++] = p ;
	  }
	if (n)
	  {
	    qq[n] = 0 ;
	    pmx1->pp[N*i + j] = polMultiSum (pmx1->h, qq) ;
	  }
      }
  return pmx1 ;
} /* pmxMultiSum */

/*************************************************************************************/

PMX pmxSum (PMX pmx1, PMX pmx2, char *title, AC_HANDLE h)
{
  PMX pmxs[3] ;
  pmxs[0] = pmx1 ;
  pmxs[1] = pmx2 ;
  pmxs[2] = 0 ;
  return pmxMultiSum (pmxs, title, h) ;
} /* pmxSum */

/*************************************************************************************/
/*************************************************************************************/
/* create a polynome multiple of the identity matrix */
PMX pmxScalar (int N, char *title, POLYNOME p,  AC_HANDLE h)
{
  PMX pmx = pmxCreate (N, title, h) ;
  for (int i = 0 ; i < N ; i++)
    pmx->pp[N*i + i] = polCopy (p, pmx->h) ;
  return pmx ;
} /* pmxScalar */

/*************************************************************************************/
/* fill a pmx with a copy of a given polynome scaled according to the complex matrix zz 
 * the size of zz is checked by requesting a magic value at last position
 */
BOOL pmxSet (PMX pmx, POLYNOME p,  double complex *zz)
{
  pmxCheck (pmx) ;
  int N = pmx->N ;

  ac_free (pmx->h) ;
  pmx->h = ac_new_handle () ;
  if (zz[N*N] != -1)
    messcrash ("As a substitute to size checking, pmxSet requires zz[N*N] = -1") ;
  for (int i = 0 ; i < N ; i++)
    for (int j = 0 ; j < N ; j++)
      {
	double complex z = zz[N*i + j] ;
	if (z != 0) 
	  {
	    POLYNOME q = polCopy (p, pmx->h) ;
	    polScale (q, z) ;
	    pmx->pp[N*i + j] = q ;
	  }
	else
	  pmx->pp[N*i + j] = 0 ;
      } 
  
  return TRUE ;
} /* pmxScalar */

/*************************************************************************************/
/* multiply a pair of PMX */
PMX pmxProduct (PMX pmx1, PMX pmx2,  char *title, AC_HANDLE h)
{
  int N = 0 ;
  PMX pmx ;
  AC_HANDLE h1 = ac_new_handle () ;
  
  if (!pmx1 || ! pmx2)
    messcrash ("Bad call to pmxProduct, null pmx") ;
  pmxCheck (pmx1) ;
  pmxCheck (pmx2) ;

  N = pmx1->N ;
  if (N != pmx2->N)
    messcrash ("Bad call to pmxProduct, unequal dimensions N1=%d N2=%2 %s %s", N, pmx2->N, pmx1->title, pmx2->title) ;

  pmx = pmxCreate (N, title, h) ;
  for (int i = 0 ; i < N ; i++)
    for (int j = 0 ; j < N ; j++)
      {
	int kk = 0 ;
	POLYNOME ab[N+1] , pp = 0 ;
	for (int k = 0 ; k < N ; k++)
	  {
	    POLYNOME a = pmx1->pp[N*i + k] ;	    
	    POLYNOME b = pmx2->pp[N*k + j] ;
	    if (a && b)
		ab[kk++] = polProduct (a, b, h1) ;
	  }
	ab[kk] = 0 ;
	if (kk)
	  pp = polMultiSum (pmx->h, ab) ;
	pmx->pp [N*i + j] = pp ;
      }
  if (1) ac_free (h1) ;
  return pmx ;
} /* pmxProduct */

/*************************************************************************************/
/* multiply a zero terminated set of PMX */
PMX pmxMultiProduct (PMX *pmxs, char *title, AC_HANDLE h)
{
  int N = 0, nn = 0 ;
  PMX pmx = 0, *ppx ;
  for (ppx = pmxs ; *ppx ; ppx++)
    {
      PMX px = *ppx ;
      pmxCheck (px) ;
      nn++ ;
      if (!N)
	N = px->N ;
      if (N != px->N)
	messcrash ("Inhomogenous call to pmxAdd N=%d != px->N = %d in %s", N, px->N, px->title) ;
    }
  switch (nn)
    {
    case 0: break ;
    case 1: pmx = pmxCopy (pmxs[0],  title, h) ; break ;
    case 2: pmx =pmxProduct (pmxs[0], pmxs[1], title, h) ; break ;
    default:
      {
	AC_HANDLE h1 = ac_new_handle () ;
	pmx =pmxProduct (pmxs[0], pmxs[1], title, h1) ;
	for (int k = 2 ; k < nn - 1 ; k++)
	  pmx =pmxProduct (pmx, pmxs[k], title, h1) ;
	pmx =pmxProduct (pmx, pmxs[nn-1], title, h) ;
	ac_free (h1) ;
	break ;
      }
    }
  return pmx ;
} /* pmxMultiProduct */

/*************************************************************************************/
  /* Taylor expansion of exp(pmx)  up to order level */
PMX pmxExponential (PMX pmx, char *title, int level, AC_HANDLE h)
{
  int N = pmx->N ;
  AC_HANDLE h1 = ac_new_handle () ;
  POLYNOME p ; 
  PMX ppp, px0 ;
  PMX pmxs[level+2] ;
  double z = 1 ;

  pmxCheck (pmx) ;
  p = newScalar (1, h1) ;
  px0 = pmxScalar (N, title, p, h1) ;
  pmxs[0] =  px0 ;  /* unit matrix */

  for (int k = 1 ; k <= level ; k++)
    {
      PMX px  = pmxProduct (pmxs[k-1], pmx, title, h1) ;

      z /= k ;
      p = newScalar (z, h1) ;
      px0 = pmxScalar (N, title, p, h1) ; /* diag 1/k matrix */
      px  = pmxProduct (px0, px, title, h1) ;
      if (1) pmxExpand (px) ;
      pmxs[k] = px ;
    }
  pmxs[level+1] = 0 ;
       
  ppp = pmxMultiSum (pmxs, title, h) ;
  pmxExpand (ppp) ;
  if (1) ac_free (h1) ;
  return ppp ;
} /* pmxExponential */

/*************************************************************************************/
/*************************************************************************************/
/* create a pmx removing line ii and column 0 */
static PMX pmxMinor (PMX pmx, int ii, AC_HANDLE h)
{
  int N = pmx->N ;
  PMX px = pmxCreate (N-1, "minor", h) ;
  pmxCheck (pmx) ;
  for (int i = 0, i1 = 0 ; i < N ; i++)
    {
      if (i != ii)
	{
	  for (int j = 1 ; j < N ; j++)
	    {
	      POLYNOME p = pmx->pp[N*i + j] ;
	      if (p)
		p = polCopy (p, px->h) ;
	      px->pp[(N-1)*i1 + (j-1)] = p ; 
	    }
	  i1++ ;
	}
    }
  return px ;
} /* pmxMinor */

/*************************************************************************************/

static POLYNOME pmxDeterminantDo (PMX pmx, POLYNOME pAbove, AC_HANDLE h) 
{
  POLYNOME p1, p2, p3, det = 0 ;
  int N = pmx->N, n = 0 ; 
  POLYNOME pp[N+1] ;
  static int pass = 0 ;
  AC_HANDLE h1 = ac_new_handle () ;

  pass++ ;
  pmxCheck (pmx) ;
  switch (N)
  {
  case 1:
    p1 = pmx->pp[0] ;
    if (p1)
      {
	p2 = polProduct (pAbove, p1, h) ;
	det = expand (p2) ;
      }
    break ;
  default:
    for (int i = 0 ; i < N ; i++)
      {
	POLYNOME p1 = pmx->pp[N * i] ;
	p2 = polProduct (pAbove, p1, h1) ;
	p3 = expand (p2) ;
	if (p3)
	  {
	    PMX px = pmxMinor (pmx, i, h1) ;
	    POLYNOME p4 = pmxDeterminantDo (px, p3, h1) ;
	    if (p4)
	      {
		if (i % 2)
		  polScale (p4, -1) ; 
		pp[n++] = p4 ;
	      }
	  }
      }
    if (n)
      {
	pp[n] = 0 ;
	p1 = polMultiSum (h1, pp) ;
	p2 = expand (p1) ;
	det = polCopy (p2, h) ; 
      }
  }
  ac_free (h1) ;
  return det ;
} /* pmxDeterminantDo */

/*************************************************************************************/

POLYNOME pmxDeterminant (PMX pmx, AC_HANDLE h)
{
  POLYNOME det = 0 ;
  pmxCheck (pmx) ;
  int N = pmx->N ;
  switch (N)
    {
    case 1:
      det = polCopy (pmx->pp[0], h) ;
      break ;
    default:
      {  
	AC_HANDLE h1 = ac_new_handle () ;
	POLYNOME pAbove = newScalar (1, h1) ;
	POLYNOME p1 = pmxDeterminantDo (pmx, pAbove, h1) ;
	det = polCopy (p1, h) ;
	ac_free (h1) ;
      }
    }
  return  det ;
} /* pmxDeterminant */

/*************************************************************************************/
/*************************************************************************************/

/* hello */
