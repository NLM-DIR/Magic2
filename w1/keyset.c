/*  File: keyset.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 ** public :  keySetAND OR XOR MINUS Order
 **   In addition keyset.h defines as macros                 
 **   keySetInsert/Remove, Create/Destroy , Max, Sort, Find  
 **  consider also w7/ksetdisp.c                             

        A KEYSET is an orderd Array of KEYs                  

 * Exported functions:
 * HISTORY:
 * Last edited: Feb 15 21:19 1994 (rd)
 * Created: Fri Jun  5 16:18:58 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: keyset.c,v 1.19 2016/09/09 23:29:00 mieg Exp $ */

#include "ac.h"
#include "bigarray.h"
#define KEYSIZE sizeof(KEY)

/**************************************************************/
      /* Pure chronologic order */

int keySetOrder (const void *a, const void *b)
{
  return
    (*(const KEY*)a) > (*(const KEY*)b)  ?
      1 : (*(const KEY*)a) == (*(const KEY*)b) ? 0 : -1 ;
}

/*************************************************************************************/
/* keep only the keys repeated min times */
unsigned int keySetNCompress (KEYSET ks, unsigned int min)
{
  unsigned int i, j, nn ;
  KEY old = 0, new ; ;
  
  for (i = j = nn = 0 ; i < keySetMax (ks) ; i++)
    {
      new = keySet (ks, i) ;
      if (new == old) 
	nn++ ;
      else
	{
	  if (nn >= min)
	    keySet (ks, j++) = old ;
	  nn = 0 ;
	  old = new ; 
	}
    }
  if (nn >= min)
    keySet (ks, j++) = old ;
  keySetMax (ks) = j ;

  return j ;
} /* keySetNCompress */

/**************************************************************/
  /* This routine returns */
  /* the intersect (logical AND) of the first and second KEYSET */

KEYSET keySetAND(KEYSET x, KEYSET y)
{ 
  register unsigned int i = 0 , j = 0, k= 0 ;
  KEYSET z = keySetCreate();

  if(!x || !y )
    return z ;
  
  if ((x && x->size != KEYSIZE) || (y && y->size != KEYSIZE))
     messcrash ("keySetOR called on non KEY Array") ;

  while((i<keySetMax(x)) && (j<keySetMax(y)))
    { 
                      /*success, skip further in the 3 index */
      if(arr(x,i,KEY) == arr(y,j,KEY))
	{ array(z,k++,KEY) = arr(x,i++,KEY) ; 
	  j++ ;
	}
      else      /*recall that every index is in increasing order*/
	(keySetOrder(arrp(x,i,KEY),arrp(y,j,KEY)) < 0) ?  i++ : j++ ;
    }
 return z ;
}

/**************************************************************/
  /* This routine returns */
  /* the union (logical OR) of the first and second KEYSET */

KEYSET keySetOR(KEYSET x, KEYSET y)
{ 
  register unsigned int i, j, k ;
  KEYSET z = 0 ;

  if(!x && !y )
    return  keySetCreate() ;
  
  if ((x && x->size != KEYSIZE) || (y && y->size != KEYSIZE))
     messcrash ("keySetOR called on non KEY Array") ;

  if(!x)
    { 
      z = arrayCopy(y) ;
      return z ;
    }

  if(!y)
    { 
      z = arrayCopy(x) ;
      return z ;
    }

  z = keySetCreate() ;
  i = j = k = 0 ;
  while((i < keySetMax(x)) && (j < keySetMax(y)))
    { 
                      /*no doubles, skip further in the 3 index */
      if(arr(x,i,KEY) == arr(y,j,KEY))
	{ array(z,k++,KEY) = arr(x,i++,KEY) ; 
	  j++ ;
	}
      else      /*recall that every index is in increasing order*/
	if (keySetOrder(arrp(x,i,KEY),arrp(y,j,KEY)) < 0) 
	  array(z,k++,KEY) =  arr(x,i++,KEY) ;
	else
	  array(z,k++,KEY) =  arr(y,j++,KEY) ;
    }

  if(i == keySetMax(x))
    while (j < keySetMax(y))
      array (z, k++, KEY) =  arr (y, j++, KEY) ;
  else if (j == keySetMax(y))
    while(i < keySetMax(x))
      array (z, k++, KEY) =  arr (x, i++, KEY) ;
  
  return z ;
}

/**************************************************************/
  /* This routine returns */
  /* the exclusive union (logical XOR) of the first and second KEYSET */

KEYSET keySetXOR(KEYSET x, KEYSET y)
{ 
  register unsigned int i = 0 , j = 0, k= 0 ;
  KEYSET z = 0 ;


  if(!x && !y )
    return keySetCreate() ; 
  
  if ((x && x->size != KEYSIZE) || (y && y->size != KEYSIZE))
     messcrash ("keySetXOR called on non KEY Array") ;

  if(!x)
    { 
      z = arrayCopy(y) ;
      return z ;
    }

  if(!y)
    { 
      z = arrayCopy(x) ;
      return z ;
    }

  z = keySetCreate() ;
  while((i < keySetMax(x)) && (j < keySetMax(y)))
    { 
                /*discard the intersect */
      if(arr(x,i,KEY) == arr(y,j,KEY))
	{ i++; j++ ;}	
      else      /*recall that every index is in increasing order*/
	{ 
	  if (keySetOrder(arrp(x,i,KEY),arrp(y,j,KEY)) < 0)
	    array(z,k++,KEY) =  arr(x,i++,KEY) ;
	  else
	    array(z,k++,KEY) =  arr(y,j++,KEY) ;
	}
    }

 if(i == keySetMax(x))
   while(j<keySetMax(y))
     array(z,k++,KEY) =  arr(y,j++,KEY) ;
 else if(j == keySetMax(y))
   while(i<keySetMax(x))
     array(z,k++,KEY) =  arr(x,i++,KEY) ;

 return z ;
}

/**************************************************************/
  /* This routine returns */
  /* the difference (logical AND NOT) of the first and second KEYSET */

KEYSET keySetMINUS(KEYSET x, KEYSET y)
{ 
  register unsigned int i = 0 , j = 0, k= 0 ;
  KEYSET z = 0 ;

  if(!x)
    return keySetCreate();
  
  if(!y)
    { 
      z = arrayCopy(x) ;
      return z ;
    }

  if ((x && x->size != KEYSIZE) || (y && y->size != KEYSIZE))
     messcrash ("keySetMINUS called on non KEY Array") ;

  z = keySetCreate();
  while((i<keySetMax(x)) && (j<keySetMax(y)))
     { 
                      /*no doubles, skip further in x y */
      if(arr(x,i,KEY) == arr(y,j,KEY))
	{ i++ ; j++ ;
	}
      else      /*recall that every index is in increasing order*/
	if (keySetOrder(arrp(x,i,KEY),arrp(y,j,KEY)) < 0) 
	  array(z,k++,KEY) =  arr(x,i++,KEY) ;
	else
	  j++ ;
    }

  while(i<keySetMax(x))
     array(z,k++,KEY) =  arr(x,i++,KEY) ;
 
 return z ;
}

/**************************************************************/

/*  input: kss is a set of sets
 * output: kss is a set of exclusive sets
 *         all pairwise intersects are now empty
 */

unsigned int keySetPartition (Array kss)
{
  unsigned int ii, jj, ni, nj, nij, imax ;
  KEYSET ki, kj, kij, kinew, kjnew ;

  if (! arrayExists (kss))
    messcrash ("keySetPartition received a wrong array, should be a set of sets") ;
  for (ii = 0 ; ii < arrayMax (kss) ; ii++)
    {
      ki = arr (kss, ii, KEYSET) ;
      if (!ki ||
	  ! (ni = keySetMax (ki))
	  )
	continue ;
      if (! keySetExists (ki))
	messcrash ("keySetPartition received a wrong array, should be a set of sets") ;
      imax =  arrayMax (kss) ; /* no use to compare to the new pieces */
      for (jj = ii + 1 ; jj < imax ; jj++)
	{
	  kj = arr (kss, jj, KEYSET) ;
	  if (!kj ||
	      ! (nj = keySetMax (kj))
	      )
	    continue ;
	  kij = keySetAND (ki, kj) ;
	  nij = keySetMax (kij) ;
	  if (!nij) ;
	  else if (nij == ni && nij == nj)
	    nj = keySetMax (kj) = 0 ;
	  else
	    {
	      kinew = keySetMINUS (ki, kij) ;
	      kjnew = keySetMINUS (kj, kij) ;
	      arr (kss, ii, KEYSET) = kij ;
	      arr (kss, jj, KEYSET) = kinew ;
	      array (kss, arrayMax (kss), KEYSET) = kjnew ;
	      keySetDestroy (ki) ; keySetDestroy (kj) ;
	      ki = arr (kss, ii, KEYSET) ;
	      kj = arr (kss, jj, KEYSET) ;
	      kij = 0 ; /* avoid destruction */
	      ni = keySetMax (ki) ;
	      nj = keySetMax (kj) ;
	    }
	  keySetDestroy (kij) ;
	}
    }
  for (ii = jj = 0 ; ii < arrayMax (kss) ; ii++)
    {
      ki = arr (kss, ii, KEYSET) ;
      if (!ki ||
	  ! (ni = keySetMax (ki))
	  )
	keySetDestroy (ki) ;
      else
	{
	  if (jj < ii)
	    arr (kss, jj, KEYSET) = ki ;
	  jj++ ;
	}
    }
  arrayMax (kss) = jj ;

  return arrayMax (kss) ;    
} /* keySetPartition */

/******************************************************************/
/* this algorithm is from Knuth and gives smaller rounding errors */
/******************************************************************/

static int floatVarianceSorted (Array a, float *medianp, float *averagep, float *sigmap)
{
  float *zp1, *zp2 ;
  double delta = 0, shiftedMean = 0, M2 = 0, median = 0, sigma = 0 ;
  int n, n1, n2, N = arrayMax (a) ;

  n = N/2 ;
  if (N == 0)
    messcrash (" floatVarianceSorted received an empty array ") ;

  if (2 * n == arrayMax(a))
    {
      median = (arr (a, n - 1, float) + arr (a, n, float)) / 2 ;
    }
  else
    {
      median = arr (a, n, float) ;
    }
  n1 = n ; n2 = n - 1 ;
  if (N == 1)
    {
      shiftedMean = sigma = 0 ; /* shiftedMean is centralized on the median */
    }
  else
    {
      /* starting the additions from the center stabilizes the computation */
      for (zp1 = arrp (a, n1, float), zp2 = arrp (a, n2, float), n = 1 ; n1 < N ; n1++, zp1++, n2--, zp2--)
	{
	  if (n1 < N)
	    {
	      delta = *zp1 - shiftedMean - median ;  /* subtracting the median every where stabilizes the calculation */
	      shiftedMean = shiftedMean + delta/n ;
	      /*  M2 = M2 + delta * (*zp1 - shiftedMean - median) ; */
	      M2 = M2 + delta * delta * (n - 1.0) / n ;
	      n++ ;
	    }
	  if (n2 >= 0)
	    {
	      delta = *zp2 - shiftedMean - median ;  /* subtracting the median every where stabilizes the calculation */
	      shiftedMean = shiftedMean + delta/n ;
	      /*  M2 = M2 + delta * (*zp2 - shiftedMean - median) ; */
	      M2 = M2 + delta * delta * (n - 1.0) / n ;
	      n++ ;
	    }
	}
      /* sigma =  sqrt (M2/(N - 1)) ;*/
      sigma =  sqrt (M2/(N)) ;
    }
  if (medianp) *medianp = median ;
  if (sigmap) *sigmap = sigma ;
  if (averagep) *averagep = shiftedMean + median ;

  return N ;
} /* floatVarianceSorted */

/**************************************************************/

int floatVariance (Array a, float *medianp, float *averagep, float *sigmap)
{
  int N ;
  Array b = 0 ;

  if (medianp) *medianp = 0 ;
  if (averagep) *averagep = 0 ;
  if (sigmap) *sigmap = 0 ;
  
  if (! arrayExists (a) || !arrayMax (a))
    return 0 ;

  if (a->size != sizeof(float))
    messcrash ("floatVariance not called on an array of float") ;

  b = arrayCopy (a) ;
  arraySort (b, floatOrder) ;
  
  N = floatVarianceSorted (a, medianp, averagep, sigmap) ;
  arrayDestroy (b) ;
  
  return N ;
} /* floatVariance */

/**************************************************************/
/**************************************************************/

static int bigFloatVarianceSorted (BigArray a, float *medianp, float *averagep, float *sigmap)
{
  float *zp1, *zp2 ;
  double delta, shiftedMean = 0, M2 = 0, median, sigma ;
  long int n, n1, n2, N = arrayMax (a) ;

  n = N/2 ;
  if (N == 0)
    messcrash (" bigFloatVarianceSorted received an empty array ") ;

  if (2 * n == bigArrayMax(a))
    {
      median = (bigArr (a, n - 1, float) + bigArr (a, n, float)) / 2 ;
    }
  else
    {
      median = bigArr (a, n, float) ;
    }
  n1 = n ; n2 = n - 1 ;
  if (N == 1)
    {
      shiftedMean = sigma = 0 ; /* shiftedMean is centralized on the median */
    }
  else
    {
      /* starting the additions from the center stabilizes the computation */
      for (zp1 = bigArrp (a, n1, float), zp2 = bigArrp (a, n2, float) ; n1 < N ; n1++, zp1++, n2--, zp2--)
	{
	  if (n1 < N)
	    {
	      delta = *zp1 - shiftedMean - median ;  /* subtracting the median every where stabilizes the calculation */
	      shiftedMean = shiftedMean + delta/n ;
	      M2 = M2 + delta * (*zp1 - shiftedMean - median) ;
	    }
	  if (n2 >= 0)
	    {
	      delta = *zp2 - shiftedMean - median ;  /* subtracting the median every where stabilizes the calculation */
	      shiftedMean = shiftedMean + delta/n ;
	      M2 = M2 + delta * (*zp2 - shiftedMean - median) ;
	    }
	}
      /* sigma =  sqrt (M2/(N - 1)) ;*/
      sigma =  sqrt (M2/(N)) ;
    }
  if (sigmap) *sigmap = sigma ;
  if (averagep) *averagep = shiftedMean + median ;
 
  return N ;
} /* bigFloatVarianceSorted */

/**************************************************************/
/* for i < LL  bigArray(a, i, float) contains the number of times i has been seen
 * for i >= LL bigArray(a, i, float) enumerates the values
 * this is a very strong optimisiation when measuring huge SNP tables
 */

int bigFloatVariance (BigArray a,  int LL, float *medianp, float *averagep, float *sigmap)
{
  int N ;
  long int i, j, n = 0 ;
  BigArray b = 0 ;

  if (medianp) *medianp = 0 ;
  if (averagep) *averagep = 0 ;
  if (sigmap) *sigmap = 0 ;
  
  if (! bigArrayExists (a) || ! bigArrayMax (a))
    return 0 ;

  if (a->size != sizeof(float))
    messcrash ("floatVariance not called on an array of float") ;

  for (i = 0 ; i < LL ; i++)
    n +=  bigArray (a, i, float) ;
  if (n < 0)
    return 0 ;
  b = bigArrayCreate (n + bigArrayMax (a), float) ;
  bigArray (b, n + bigArrayMax (a), float) = 0 ; /* make room */

  /* register the small values stored by multiplicities in a */
  for (i = 0, j = 0 ; j < LL ; j++)
    {
      n  = bigArray (a, j, float) ; if (n < 0) n = 0 ;
      while (n--)
	bigArray (a, i++, float) = j ;
    }
  /* register the larger values stored as a list in a */
  n =  bigArrayMax (a) ;
  for (j = LL ; j < n ; j++)
    {
      bigArray (a, i++, float) = j ;
    }

  bigArraySort (b, floatOrder) ;
  
  N = bigFloatVarianceSorted (a, medianp, averagep, sigmap) ;
  bigArrayDestroy (b) ;
  
  return N ;
} /* bigFloatVariance */

/******************************************************************/
/* this algorithm is from Knuth and gives smaller rounding errors */
/******************************************************************/

static int doubleVarianceSorted (Array a, double *medianp, double *averagep, double *sigmap)
{
  double *zp1, *zp2 ;
  double delta = 0, shiftedMean = 0, M2 = 0, median = 0, sigma = 0 ;
  int n, n1, n2, N = arrayMax (a) ;

  n = N/2 ;
  if (N == 0)
    messcrash (" doubleVarianceSorted received an empty array ") ;

  if (2 * n == arrayMax(a))
    {
      median = (arr (a, n - 1, double) + arr (a, n, double)) / 2 ;
    }
  else
    {
      median = arr (a, n, double) ;
    }
  n1 = n ; n2 = n - 1 ;
  if (N == 1)
    {
      shiftedMean = sigma = 0 ; /* shiftedMean is centralized on the median */
    }
  else
    {
      /* starting the additions from the center stabilizes the computation */
      for (zp1 = arrp (a, n1, double), zp2 = arrp (a, n2, double), n = 1 ; n1 < N ; n1++, zp1++, n2--, zp2--)
	{
	  if (n1 < N)
	    {
	      delta = *zp1 - shiftedMean - median ;  /* subtracting the median every where stabilizes the calculation */
	      shiftedMean = shiftedMean + delta/n ;
	      /*  M2 = M2 + delta * (*zp1 - shiftedMean - median) ; */
	      M2 = M2 + delta * delta * (n - 1.0) / n ;
	      n++ ;
	    }
	  if (n2 >= 0)
	    {
	      delta = *zp2 - shiftedMean - median ;  /* subtracting the median every where stabilizes the calculation */
	      shiftedMean = shiftedMean + delta/n ;
	      /*  M2 = M2 + delta * (*zp2 - shiftedMean - median) ; */
	      M2 = M2 + delta * delta * (n - 1.0) / n ;
	      n++ ;
	    }
	}
      /* sigma =  sqrt (M2/(N - 1)) ;*/
      sigma =  sqrt (M2/(N)) ;
    }
  if (medianp) *medianp = median ;
  if (sigmap) *sigmap = sigma ;
  if (averagep) *averagep = shiftedMean + median ;

  return N ;
} /* doubleVarianceSorted */

/**************************************************************/

int doubleVariance (Array a, double *medianp, double *averagep, double *sigmap)
{
  int N ;
  Array b = 0 ;

  if (medianp) *medianp = 0 ;
  if (averagep) *averagep = 0 ;
  if (sigmap) *sigmap = 0 ;
  
  if (! arrayExists (a) || !arrayMax (a))
    return 0 ;

  if (a->size != sizeof(double))
    messcrash ("doubleVariance not called on an array of double") ;

  b = arrayCopy (a) ;
  arraySort (b, doubleOrder) ;
  
  N = doubleVarianceSorted (a, medianp, averagep, sigmap) ;
  arrayDestroy (b) ;
  
  return N ;
} /* doubleVariance */

/**************************************************************/
/* divariance: look for a dromadaire
 * find the point which splits the data in 2 parts of minimal combined variance
 *  retain this minimal variance
 *
 * SIDE-EFFECT  a is sorted
 */
BOOL diVariance (Array a, float *splitp, float *avp, float *sigmap, float *av1p, float *sigma1p, float *av2p, float *sigma2p, float *sigmaCombinedp)
{
  int i, iMax = arrayMax (a) ;
  int NXa, NYa, NXb, NYb, limit = 1, pass ;
  float median = 0, split = 0 ;
  unsigned int delta ;
  double z = 0, z1, z2, bestz = 0 ;
  double X1a, X2a, Y1a, Y2a, za = 0, z1a, z2a ;
  double X1b, X2b, Y1b, Y2b, zb = 0, z1b, z2b ;
  
  *splitp = *avp = *sigmap = *av1p = *av2p = *sigma1p =  *sigma2p = *sigmaCombinedp = 0 ;
  if (! arrayExists (a) || !arrayMax (a))
    return 0 ;

  if (a->size != sizeof(float))
    messcrash ("floatVariance not called on an array of float") ;

  arraySort (a, floatOrder) ;
  floatVarianceSorted (a, &median, avp, sigmap) ;

  *av1p = *av2p = *avp ;
  *sigma1p = *sigma2p = *sigmap ;
  /* start at median and by dicothomy find optimal variance */
  if (iMax < 20) return 0 ;

  /* start at center */
  limit = iMax / 2 ;
  /* compute the 2 variances */
  for (pass = 0, delta = limit/2 ; delta > 1 ; pass++, delta >>= 1)
    {
      NXa = NYa = 0 ; X1a = X2a = Y1a = Y2a = 0 ; /* a way to compute the derivative */
      NXb = NYb = 0 ; X1b = X2b = Y1b = Y2b = 0 ;
      for (i = 0 ; i < iMax ; i++)  /* number 1 to nRuns, easier for the binomial coef */
	{
	  z = array (a, i, float) ;
	  if (i == limit) split = z ;
	  if (i < limit) { NXa++ ; X1a += z ; X2a += z * z ; }
	  else { NYa++ ; Y1a += z ; Y2a += z * z ; }

	  if (i <= limit) { NXb++ ; X1b += z ; X2b += z * z ; }
	  else { NYb++ ; Y1b += z ; Y2b += z * z ; }
	}
      if (NXa < 1) { NXa = 1 ; }   if (NYa < 1) { NYa = 1 ; }  
      X1a /= NXa ; X2a /= NXa ; Y1a /= NYa ; Y2a /= NYa ;
      z1a = X2a - X1a*X1a ; z2a = Y2a - Y1a*Y1a ; za = z1a + z2a ;

      if (NXb < 1) { NXb = 1 ;}  if (NYb < 1) { NYb = 1 ;} 
      X1b /= NXb ; X2b /= NXb ; Y1b /= NYb ; Y2b /= NYb ;
      z1b = X2b - X1b*X1b ; z2b = Y2b - Y1b*Y1b ; zb = z1b + z2b ;

      z = (za + zb)/2 ;
      if (2*delta + 1 >= limit || z < bestz)
	{ bestz = z ; *splitp = split ; *av1p = X1a ; *av2p = Y1a ; *sigma1p = z1a ; *sigma2p = z2a ; }

      if (za < zb) limit -= delta ;
      else limit += delta ;
      z = z1a - z2a ; if (z<0)z = -z ; if (2000 * z < z1a + z2a) break ;
    }
  z1 = *sigma1p ; *sigma1p = sqrt (z1) ;
  z2 = *sigma2p ; *sigma2p = sqrt (z2) ;
  z = *sigmap ;
  *sigmaCombinedp = sqrt ((z1 + z2) * 1.3759691970109549) ; /* Pi / (2*Pi - 4) */

  /* Explanation of the constant Pi/ (2Pi - 4)

    The correctly normalized Gaussian is (1/sqrt(2*Pi)) exp (-x^2/2)
    The surface is 1, the average is 0 and the variance is 1 = <x^2> since <x> = 0

    If we take the halh Gaussian, 
    the (even power) surface is 1/2 and the even is by symmetry sum(x^2) = 1/2, so <x^2>=1
    But the average must be computed as  
        {1/surface == 2} integral {1/sqrt(2 Pi) (dx x exp(-x^2/2))} 
           = 2/sqrt(2Pi) { intregral(d(exp(-x^2/2)))==1}
           = 2/sqrt(2 Pi)
    So the variance of the half Gaussian is <x^2> - (<x>)^2 = 1 - 2/Pi
    The sum of the 2 halves is   2 - 4/Pi == inverse of Pi/(2 Pi - 4)      QED
  */

  return z > 1.3 * (*sigmap) ? TRUE : FALSE ;
} /* diVariance */

/**************************************************************/
/**************************************************************/

/* Solves, y = a x + b , in a and b   a defaults to *ap  if too few data given */

BOOL linearRegression (Array xy, double *ap, double *bp, double *rp, double *wp)
{
  int i, max ;
  double r, w, xm,ym,xym , x2m, y2m ,b, a;

  if(!xy || !(max = arrayMax(xy)))
    {
      messout("linearRegression received a null array") ;
      *ap = * bp = *rp = *wp  = 0 ;
      return FALSE ;
    }
  if(xy->size != sizeof(POINT2D))
    messcrash(
  "linearRegression received a wrong type of Array, should be  a pair of doubles ") ;

  i = max ;
  xm = ym = xym = 0 ;
  while (i--)
    {
      xm += arr(xy,i,POINT2D).x ; ym += arr(xy,i,POINT2D).y;
    }
  xm /= max ; ym /= max;
  x2m = y2m =xym =0;
  i = max ;
  while (i--)
    {
      xym += (arr(xy,i,POINT2D).x-xm) *( arr(xy,i,POINT2D).y-ym) ;
      x2m +=  (arr(xy,i,POINT2D).x-xm) *( arr(xy,i,POINT2D).x-xm) ;
      y2m +=  (arr(xy,i,POINT2D).y-ym) *( arr(xy,i,POINT2D).y-ym) ;
    }
  if(max >3 && x2m>0 && y2m >0)
    {
      r = xym / sqrt((double)x2m * y2m) ;
      if (r >= 1)    /* mhmp 28.05.99 */
	w = 0.999 ;
      else
   	w = r ;
      i = max - 3 ;
      w = .5 * log((1+w) / (1-w)) * sqrt((double) i) ;
    }
  else
    r = w = 0 ;

  i = max ;
  a = 0 ;
  while (i--)
    a += (arr(xy,i,POINT2D).x - xm) * (arr(xy,i,POINT2D).y - ym) ;

  if (x2m)
    { a /= x2m ; }
  else
    a = *ap ;  /* default */

  b = ym - a*xm ;
  
  *ap = a ;
  *bp = b ;
  *rp = r ;
  *wp = w ;
  
  return TRUE ;
}

/*************************************************************************/
/*************************************************************************/

