#include "ac.h"
#include <complex.h>
#include "matrix.h"
#include "polynome.h"
#ifdef JUNK8
/* Create july 2026
 *
 * Calculations supporting the BIM SU(2/1) coset space paper
 *
 * The 8 complex matrices are enumerated, and we then verify the commutators and the anomalies
 * we have explicitly the leptons, the quarks, and several indecomposable reps

 * A: Define the su(2/1) matrices K0 in the fundamental quartet rep Kac(a=0,y=0)
 *    construct all the structure constants and the Kac Killing
      K0 is an array of 12 matrices (0,1,2,3) even, (4,5,6,7) odd (8) = chi   (9) not used
 * B: Construct the Matrices Kac(a,y) for the grid a=0-5, y = 0-5
      verify the commutators
      compute the metric
      KKK is the array of all reps labelled as KKK[n], with n=100*a + y (no collisions possible)
 * C: Compute the Casimirs C2 and C3 as matrices
      Compute the casimir tensor C3_{YYY} C3_{Yhh}, verify that they are proportional to Tr(Y)
      Verify that they are diagonal, compute their eigen values,
      verify C3 = (y-1) C2; C2 = (y-a)(y+a+2) Tr(Y)
 * 
 */

typedef struct mmStruct { int s, d, d1, d2, d3, d4 ; MX mx ; } SMM ;
typedef struct fabcStruct {int a,b,c,sign; complex float z; const char *title ;} FABC ;
typedef struct kacStruct { int a, b ; int d, d1, d2,  d3, d4 ; MX chi ; MX SMM[12] ; AC_HANDLE h ; } KAC ;

static SMM smmCreate (AC_HANDLE h, char *title, int s, int d, int d1, int d2, int d3, int d4)
    {
      SMM mm ;
      mm.mx = mxCreate (h,  title, MX_INT, d, d, 0) ;
      mm.s = s ;
      mm.d = d ;
      mm.d1 = d1 ;
      mm.d2 = d2 ;
      mm.d3 = d3 ;
      mm.d4 = d4 ;

      return mm ;
    }

static void smmShow (SMM *mm)
{
  mxShow (mm->mx) ;
  return ;
}

// Easy to use supermatrix element, notice the transposition compensating a transposotion in the MX lib
#define MMM(_mm,_i,_j) _mm->mx[_mm->d * (_j) + (-i)]

typedef struct bimStruct {
  AC_HANDLE h ;
  int a, y ;   /* user defined Dynkin numbers */
  int aMax, yMax ; /* user defined Dynkin grid */
  KAC kac3 ; /* fundamental rep, filled by hand */
  Array FABC  ;      /* Array of structure constants, built by using MU3 */
  Array kacArray ; /* array of Kac modules, labelled by 100a+y */
} BIM ;

static KAC* bimConstructKacModule (int a, int b, BOOL show,  AC_HANDLE h) 
{
  KAC *kac = (KAC *) halloc (sizeof(KAC), h) ;
  
  int d1 = a+1, d2 = a+2, d3 = a, d4 = a+1 ;
  int d = d1 + d2 + d3 + d4 ;
  int xx[d*d] ;
  
  
  kac->h = h ;
  kac->a = a ; kac->b = b ;
  kac->d = d ;
  kac->d1 = d1 ;
  kac->d2 = d2 ;
  kac->d3 = d3 ;
  kac->d4 = d4 ;

  kac->MM = (MX *) halloc (12 * sizeof (MX), h) ;

  kac->chi = mxCreate (h,  "chi", MX_INT, d, d, 0) ;
  // Hypercharge
  kac->MM[0] = mxCreate (h,  "muY", 0, d, d1, d2, d3, d4) ;
  // sl(2) Chevalley basis
  kac->MM[1] = mxCreate (h,  "muE", 0, d, d1, d2, d3, d4) ;
  kac->MM[2] = mxCreate (h,  "muF", 0, d, d1, d2, d3, d4) ;
  kac->MM[3] = mxCreate (h,  "muH", 0, d, d1, d2, d3, d4) ;
  // the 4 odd matrices
  kac->MM[4] = mxCreate (h,  "muU", 1, d, d1, d2, d3, d4) ;
  kac->MM[5] = mxCreate (h,  "muV", 1, d, d1, d2, d3, d4) ;
  kac->MM[6] = mxCreate (h,  "muW", 1, d, d1, d2, d3, d4) ;
  kac->MM[7] = mxCreate (h,  "muX", 1, d, d1, d2, d3, d4) ;
  // K1={U,V} and K2={W,X}
  kac->MM[8] = mxCreate (h,  "muK1", 0, d, d1, d2, d3, d4) ;
  kac->MM[9] = mxCreate (h,  "muK2", 0, d, d1, d2, d3, d4) ;
  // U' and W' carrying the b coefficient of the odd raising operators
  kac->MM[10] = mxCreate (h,  "muUb", 1, d, d1, d2, d3, d4) ;
  kac->MM[11] = mxCreate (h,  "muWb", 1, d, d1, d2, d3, d4) ;

  /* Y hypercharge  Y = diag (a,a...a/a+1,...a+1) */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    xx[d * i + i] = 2*b - a ;
  for (i = 0 ; i < d2 + d3 ; i++)
    {
      j = d1 + i ;
      xx[d * j + j] = 2*b -a - 1 ;
    }
  for (i = 0 ; i < d4 ; i++)
    {
      j = d1 + d2 + d3 + i ;
      xx[d * j + j] = 2*b -a - 2 ;
    }
  mxSet (kac->MM[0], xx) ;
  if (show) mxShow(kac->MM[0]) ;
  
  /* even Cartan operator H = diag (a, a-2 .... -a in each SU(2) sector */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    xx[d*i +i] = d1 - 1 - 2*i ;
  for (i = 0 ; i < d2 ; i++)
    xx[d*(d1+i) + d1+i) = d2 - 1 - 2*i ;
  for (i = 0 ; i < d3 ; i++)
    xx[d *(d1 + d2 + i) + d1 + d2 + i] = d3 - 1 - 2*i ;
  for (i = 0 ; i < d4 ; i++)
    xx[d*(d1 +d2 + d3 +) + d1 + d2 + d3 + i] = d4 - 1 - 2*i ;
  mxSet (kac->MM[3].mx, xx) ;
  if (show) mmShow(kac->MM[3]) ;
 
  /* odd Cartan operator K1 = diag (0,-1,-2,...-a/-1,-2...-a) */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 1 ; i < d1 ; i++)
    {
      xx[d * i + i] = -i ;
      j = d1 + i - 1 ;
      xx[d * j + j] = -i ;
    }
  mxSet (kacMM[8], xx) ;
  if (show) mxShow(kacMM[8]) ;
  
  /* odd Cartan operator K2 = diag (-a,...-2,-1, 0/--a,...-2,-1) */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < a ; i++)
    {
      xx[d * i + i] = -a + i ;
      j = d1 + i ;
      xx[d * j + j] = -a + i ;
    }
  mxSet (kacMM[9], xx) ;
  if (show) mxShow(kacMM[9]) ;
  
  
   /* even raising operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 1 ; i < d1 ; i++)
    xx[d * i + i - 1] = i * (a  - i + 1) ;
  for (i = 1 ; i < d2 ; i++)
    xx[d * (d1 + i) + d1 + i - 1] =  i * (a - i + 2) ;
  for (i = 1 ; i < d3 ; i++)
    xx[d * (d1 + d2 + i) + d1 + d2 + i - 1] =  i * (a - i ) ;
  for (i = 1 ; i < d4 ; i++)
    xx[d * (d1 + d2 + d3 + i) + d1 + d2 + d3 + i - 1] =  i * (a - i + 1) ;
  mxSet (kacMM[1], xx) ;
  if (show) mxShow(kacMM[1]) ;
  
  /* even lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 - 1 ; i++)
    xx[d * i + i + 1] = 1 ;
  for (i = 0 ; i < d2 - 1 ; i++)
    xx[d * (d1 + i) + d1 + i + 1] =  1 ;
  for (i = 0 ; i < d3 - 1 ; i++)
    xx[d * (d1 + d2 + i) + d1 + d2 + i + 1] =  1 ;
  for (i = 0 ; i < d4 - 1 ; i++)
    xx[d * (d1 + d2 + d3 + i) + d1 + d2 + d3 + i + 1] =  1 ;
  mxSet (kacMM[2], xx) ;
  if (show) mxShow(kacMM[2]) ;

 /* odd raising operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (d1+i) + i + 1] = -1 ;
    }
  mxSet (kacMM[4], xx) ;
  if (show) mxShow(kacMM[4]) ;
  
  /* b part of odd raising operator */  
  int s = 1 ;

  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    {

      j = d1 ;
      xx[d * (j + i) + i] = s*1 ;
    }
  for (i = 0 ; i < d3 ; i++)
    {
      j = d1 + d2 ;
      xx[d * (j + i ) + i + 1 ] = s * ( -1) ;
      j = d1 + d2 + d3 ;
      xx[d * (j+i) + i+d1+d2] = 1 ;
    }
  for (i = 0 ; i < d4 ; i++)
    {
      j = d1 + d2 + d3 ;
      xx[d * (j + i ) + d1 + i + 1 ] = 1  ;
    }
  mxSet (kacMM[10], xx) ;
  if (show) mxShow(kacMM[10]) ;

 /* odd lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (i + 1)  + d1 + i] = i+1 ;
    }
  mxSet (kacMM[5], xx) ;
  if (show) mxShow(kacMM[5]) ;

 /* other lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (i)  + d1 + i] = 1 ;
    }
  mxSet (kacMM[7], xx) ;
  if (show) mxShow(kacMM[7]) ;



  /* other raising operator */  
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (d1+i) + i] = -a + i ;
    }
  mxSet (kacMM[6], xx) ;
  if (show) mxShow(kacMM[6]) ;

  /* b dependent other raising operator */  
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (d1+i) + i] = -a + i ;
    }
  mxSet (kacMM[11], xx) ;
  if (show) mxShow(kacMM[11]) ;

  return kac ;
}

#ifdef JUNK

static MX *KasimirConstructAtypicMatrices (KAS *kas)
{
  MX muY, muE, muF, muH, muU, muV, muW, muX ; /* the 8 generators of SU(2/1) in the Chevalley basis */
  MX muK1, muK2 ; /* the combinations Y + H and Y - H, with must match {U,V} and {W,X} */
  MX muUb, muWb ; /* the b dependent part of the odd raisng matrices U and W */
  MX *mu ;
  int d, d1 = 0, d2 = 0, d3 = 0, d4 = 0 ;
  int i, j ;
  int a = kas->a ;
  AC_HANDLE h = kas->h ;
  int s = 1 ;
  mu = kas->mu = (MX *) halloc (12 * sizeof (MX), kas->h) ;
  kas->chi = 1 ;
  
   /* atypic 1 */
  d1 = a + 1 ;
  d2 = a ;
  
  d = d1 + d2 ;
  kas->d = d ;
  kas->d1 = d1 ;
  kas->d2 = d2 ;
  kas->d3 = 0 ;
  kas->d4 = 0 ;
      
  
  muY = mxCreate (h,  "muY", MX_INT, d, d, 0) ;
  muH = mxCreate (h,  "muH", MX_INT, d, d, 0) ;
  muE = mxCreate (h,  "muE", MX_INT, d, d, 0) ;
  muF = mxCreate (h,  "muF", MX_INT, d, d, 0) ;
  muU = mxCreate (h,  "muU", MX_INT, d, d, 0) ;
  muV = mxCreate (h,  "muV", MX_INT, d, d, 0) ;
  muW = mxCreate (h,  "muW", MX_INT, d, d, 0) ;
  muX = mxCreate (h,  "muX", MX_INT, d, d, 0) ;
  muK1 = mxCreate (h,  "muK1", MX_INT, d, d, 0) ;
  muK2 = mxCreate (h,  "muK2", MX_INT, d, d, 0) ;
  muUb = mxCreate (h,  "muUb", MX_INT, d, d, 0) ;
  muWb = mxCreate (h,  "muWb", MX_INT, d, d, 0) ;
  
  int xx[d*d] ;

  /* even Cartan operator H = diag (a, a-2 .... -a in each SU(2) sector */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    xx[d * i + i] = d1 - 1 - 2*i ;
  for (i = 0 ; i < d2 ; i++)
    {
      j = d1 + i ;
      xx[d * j + j] = d2 - 1 - 2 * i ;
    }
  mxSet (muH, xx) ;
  if (kas->show) mxShow(muH) ;
 
  /* Y hypercharge  Y = diag (a,a...a/a+1,...a+1) */
 
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    xx[d * i + i] = -a ;
  for (i = 0 ; i < d2 + d3 ; i++)
    {
      j = d1 + i ;
      xx[d * j + j] = -a - 1 ;
    }
  for (i = 0 ; i < d4 ; i++)
    {
      j = d1 + d2 + d3 + i ;
      xx[d * j + j] = -a - 2 ;
    }
  mxSet (muY, xx) ;
  if (kas->show) mxShow(muY) ;
  
  /* odd Cartan operator K1 = diag (0,-1,-2,...-a/-1,-2...-a) */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 1 ; i < d1 ; i++)
    {
      xx[d * i + i] = -i ;
      j = d1 + i - 1 ;
      xx[d * j + j] = -i ;
    }
  mxSet (muK1, xx) ;
  if (kas->show) mxShow(muK1) ;
  
  /* odd Cartan operator K2 = diag (-a,...-2,-1, 0/--a,...-2,-1) */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < a ; i++)
    {
      xx[d * i + i] = -a + i ;
      j = d1 + i ;
      xx[d * j + j] = -a + i ;
    }
  mxSet (muK2, xx) ;
  if (kas->show) mxShow(muK2) ;
  
  
   /* even raising operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 1 ; i < d1 ; i++)
    xx[d * i + i - 1] = i * (a  - i + 1) ;
  for (i = 1 ; i < d2 ; i++)
    xx[d * (d1 + i) + d1 + i - 1] =  i * (a - i ) ;
  
  mxSet (muE, xx) ;
  if (kas->show) mxShow(muE) ;
  
  /* even lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d - 1 ; i++)
    {
      if (i != d1 - 1) xx[d * i + i + 1] = 1 ;
    }
  mxSet (muF, xx) ;
  if (kas->show) mxShow(muF) ;

 /* odd raising operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (d1+i) + i + 1] = -1 ;
    }
  mxSet (muU, xx) ;
  if (kas->show) mxShow(muU) ;
  
  /* b part of odd raising operator */  
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    {

      j = d1 ;
      xx[d * (j + i) + i] = s*1 ;
    }
  for (i = 0 ; i < d3 ; i++)
    {
      j = d1 + d2 ;
      xx[d * (j + i ) + i + 1 ] = s * ( -1) ;
      j = d1 + d2 + d3 ;
      xx[d * (j+i) + i+d1+d2] = 1 ;
    }
  for (i = 0 ; i < d4 ; i++)
    {
      j = d1 + d2 + d3 ;
      xx[d * (j + i ) + d1 + i + 1 ] = 1  ;
    }
  mxSet (muUb, xx) ;
  if (kas->show) mxShow(muUb) ;

 /* odd lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (i + 1)  + d1 + i] = i+1 ;
    }
  mxSet (muV, xx) ;
  if (kas->show) mxShow(muV) ;


  /* other raising operator */  
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (d1+i) + i] = -a + i ;
    }
  mxSet (muW, xx) ;
  if (kas->show) mxShow(muW) ;

  /* b dependent other raising operator */  
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (d1+i) + i] = -a + i ;
    }
  mxSet (muWb, xx) ;
  if (kas->show) mxShow(muWb) ;

 /* other lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (i)  + d1 + i] = 1 ;
    }
  
  mxSet (muX, xx) ;
  if (kas->show) mxShow(muX) ;
  
  mu[0] = muY ; mu[1] = muE ; mu[2] = muF ; mu[3] = muH ;
  mu[4] = muW ; mu[5] = muX ; mu[6] = muU ; mu[7] = muV ; 
  mu[8] = muK1 ; mu[9] = muK2 ;
  mu[10] = muUb ; mu[11] = muWb ;

  
  return mu ;
} /* KasimirConstructAtypicMatrices */

/***********************************************************************************************************************************************/
static MX *KasimirConstructAntiMatrices (KAS *kas)
{
  MX muY, muE, muF, muH, muU, muV, muW, muX ; /* the 8 generators of SU(2/1) in the Chevalley basis */
  MX muK1, muK2 ; /* the combinations Y + H and Y - H, with must match {U,V} and {W,X} */
  MX muUb, muWb ; /* the b dependent part of the odd raisng matrices U and W */

  MX *mu ;
  int d, d1 = 0, d2 = 0, d3 = 0, d4 = 0 ;
  int i, j, k ;
  int a = kas->a ;
  int b = kas->b ;
  AC_HANDLE h = kas->h ;
  mu = kas->mu = (MX *) halloc (12 * sizeof (MX), kas->h) ;

  kas->chi = 1 ;
  
  if (1)
    {
      d1 = a + 1 ;
      d2 = a + 2 ;
    }

  d = d1 + d2 ;
  kas->d = d ;
  kas->d1 = d1 ;
  kas->d2 = d2 ;
  kas->d3 = 0 ;
  kas->d4 = 0 ;
      
  
  muY = mxCreate (h,  "muY", MX_INT, d, d, 0) ;
  muH = mxCreate (h,  "muH", MX_INT, d, d, 0) ;
  muE = mxCreate (h,  "muE", MX_INT, d, d, 0) ;
  muF = mxCreate (h,  "muF", MX_INT, d, d, 0) ;
  muU = mxCreate (h,  "muU", MX_INT, d, d, 0) ;
  muV = mxCreate (h,  "muV", MX_INT, d, d, 0) ;
  muW = mxCreate (h,  "muW", MX_INT, d, d, 0) ;
  muX = mxCreate (h,  "muX", MX_INT, d, d, 0) ;
  muK1 = mxCreate (h,  "muK1", MX_INT, d, d, 0) ;
  muK2 = mxCreate (h,  "muK2", MX_INT, d, d, 0) ;
  muUb = mxCreate (h,  "muUb: Ub", MX_INT, d, d, 0) ;
  muWb = mxCreate (h,  "muWb: Wb", MX_INT, d, d, 0) ;
  
  int xx[d*d] ;
 
  /* even Cartan operator H = diag (a, a-2 .... -a in each SU(2) sector */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0, k = a ; i < d2 ; k -= 2, i++)
    {
      if ( i < d1)
	xx[d * i + i] = k ;
      j = d1 + i ;
      xx[d * j + j] = k + 1 ;
    }
  mxSet (muH, xx) ; 
  if (kas->show) mxShow(muH) ;


  /* Y hypercharge  Y = diag (a, a-2 .... -a in each SU(2) sector */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    xx[d * i + i] = 2*b -a ;
  for (i = 0 ; i < d2 + d3 ; i++)
    {
      j = d1 + i ;
      xx[d * j + j] = 2*b -a - 1 ;
    }
  for (i = 0 ; i < d4 ; i++)
    {
      j = d1 + d2 + d3 + i ;
      xx[d * j + j] = 2*b - a - 2 ;
    }
  mxSet (muY, xx) ;
  if (kas->show) mxShow(muY) ;

  /* odd Cartan operator K = diag (a,...2,1,/ a,...2,1,0) */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      j = d1 + i  ;
      xx[d * j + j] = d1 - i ;
      if (i < d1)
	xx[d * i + i] = d1 - i ;
    }
  mxSet (muK1, xx) ;
  if (kas->show) mxShow(muK1) ;

  /* odd Cartan operator K = diag (1,2,...a/0,1,2...a) */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      j = d1 + i  ;
      xx[d * j + j] = i ;
      if (i < d1)
	xx[d * i + i] = i + 1;
    }
  mxSet (muK2, xx) ;
  if (kas->show) mxShow(muK2) ;
  
  /* even raising operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 1 ; i < d1 ; i++)
    xx[d * i + i - 1] = i * (a + 1 - i ) ;
  for (i = 1 ; i < d2 ; i++)
    xx[d * (d1 + i) + d1 + i - 1] =  i * (a - i + 2 ) ;
  mxSet (muE, xx) ;
  if (kas->show) mxShow(muE) ;

  /* even lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d - 1 ; i++)
    {
      if (i != d1 - 1) xx[d * i + i + 1] = 1 ;
    }
  mxSet (muF, xx) ;
  if (kas->show) mxShow(muF) ;

  /* odd raising operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    {
      xx[d * (d1+i) + i] = 1 ;
    }
  mxSet (muU, xx) ;
  if (kas->show) mxShow(muU) ;
  
/* odd lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    {
      xx[d * (i )  + d1 + i] = (d1 - i) ;
    }
  mxSet (muV, xx) ;
  if (kas->show) mxShow(muV) ;

  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    {
      xx[d * (d1+i+1) + i] = -1 - i ;
    }
  mxSet (muW, xx) ;
  if (kas->show) mxShow(muW) ;

  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    {
      xx[d * (i)  + d1 + i + 1] = -1 ;
    }
  mxSet (muX, xx) ;
  if (kas->show) mxShow(muX) ;
  
  mu[0] = muY ; mu[1] = muE ; mu[2] = muF ; mu[3] = muH ;
  mu[4] = muW ; mu[5] = muX ;  mu[6] = muU ; mu[7] = muV ; 
  mu[8] = muK1 ; mu[9] = muK2 ;
  mu[10] = muUb ; mu[11] = muWb ;

  return mu ;
} /* KasimirConstructAntiMatrices */

/***********************************************************************************************************************************************/

static MX *KasimirConstructTypicMatrices (KAS *kas, BOOL show)
{
  MX muY, muE, muF, muH, muU, muV, muW, muX ; /* the 8 generators of SU(2/1) in the Chevalley basis */
  MX muK1, muK2 ; /* the combinations Y + H and Y - H, with must match {U,V} and {W,X} */
  MX muUb, muWb ; /* the b dependent part of the odd raisng matrices U and W */
  MX *mu ;
  int d, d1 = 0, d2 = 0, d3 = 0, d4 = 0 ;
  int i, j ;
  int a = kas->a, b = kas->b ;
  int s = 1 ;  /* scaling U V K1 K2 */
  AC_HANDLE h = kas->h ;
  mu = kas->mu = (MX *) halloc (12 * sizeof (MX), kas->h) ;

  kas->show = show ;
  kas->chi = 1 ;
  s = a + 1 ; 
  kas->scale = s * s ;
  
  if (1)
    {
      d1 = d4 = a + 1 ;
      d2 = a + 2 ;
      d3 = a  ;
    }

  d = d1 + d2 + d3 + d4 ;
  int xx[d*d] ;
  const int *xx1 = messalloc (d*d*sizeof(int)) ;
  const int *xx2 = messalloc (d*d*sizeof(int)) ;
      
  kas->d = d ;
  kas->d1 = d1 ;
  kas->d2 = d2 ;
  kas->d3 = d3 ;
  kas->d4 = d4 ;
  
  
  muY = mxCreate (h,  "muY: Y Hypercharge", MX_INT, d, d, 0) ;
  muH = mxCreate (h,  "muH: Even SU(2) Cartan operator", MX_INT, d, d, 0) ;
  muE = mxCreate (h,  "muE: Even raising operator", MX_INT, d, d, 0) ;
  muF = mxCreate (h,  "muF: Even lowering operator", MX_INT, d, d, 0) ;
  muU = mxCreate (h,  "muU: Odd raising operator", MX_INT, d, d, 0) ;
  muV = mxCreate (h,  "muV: Odd lowering operator", MX_INT, d, d, 0) ;
  muW = mxCreate (h,  "muW: Other odd raising operator", MX_INT, d, d, 0) ;
  muX = mxCreate (h,  "muX: Other odd lowering operator", MX_INT, d, d, 0) ;
  muK1 = mxCreate (h,  "muK1: K1 = {U,V}", MX_INT, d, d, 0) ;
  muK2 = mxCreate (h,  "muK2: K2 = {W,X}", MX_INT, d, d, 0) ;
  muUb = mxCreate (h,  "muUb: Ub", MX_INT, d, d, 0) ;
  muWb = mxCreate (h,  "muWb: Wb", MX_INT, d, d, 0) ;
  
  /* even Cartan operator H = diag (a, a-2 .... -a in each SU(2) sector */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    xx[d * i + i] = d1 - 1 - 2*i ;
  for (i = 0 ; i < d2 ; i++)
    {
      j = d1 + i ;
      xx[d * j + j] = d2 - 1 - 2 * i ;
    }
  for (i = 0 ; i < d3 + d3 ; i++)
    {
      j = d1 + d2 + i ;
      xx[d * j + j] = d3 - 1 - 2 * i ;
    }
  for (i = 0 ; i < d4 ; i++)
    {
      j = d1 + d2 + d3 + i ;
      xx[d * j + j] = d4 - 1 - 2 * i ;
    }
  mxSet (muH, xx) ;
  if (kas->show) mxShow(muH) ;

  /* Y hypercharge  Y = diag (a+1,a+1...a+1/a,a,....a) */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    xx[d * i + i] = 2*b - a ;
  for (i = 0 ; i < d2 + d3 ; i++)
    {
      j = d1 + i ;
      xx[d * j + j] = 2*b -a - 1 ;
    }
  for (i = 0 ; i < d4 ; i++)
    {
      j = d1 + d2 + d3 + i ;
      xx[d * j + j] = 2*b -a - 2 ;
    }
  mxSet (muY, xx) ;
  if (kas->show) mxShow(muY) ;

  
  /* odd Cartan operator K1 = diag (a,...2,1,/ a,...2,1,0) */
  /* odd Cartan operator K2 = diag (1,2,...a/0,1,2...a) */
  memset (xx, 0, sizeof (xx1)) ;
  memset (xx, 0, sizeof (xx2)) ;
  mxValues (muY, &xx1, 0, 0) ;
  mxValues (muH, &xx2, 0, 0) ;

    
  for (i = 0 ; i < d * d ; i++)
    xx[i] =  (xx1[i] + xx2[i])/2 ;
  mxSet (muK1, xx) ;
  if (kas->show) mxShow(muK1) ;

  for (i = 0 ; i < d * d ; i++)
    xx[i] =  (xx1[i] - xx2[i])/2 ;
  mxSet (muK2, xx) ;
  if (kas->show) mxShow(muK2) ;
  

  /* even raising operator E */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 1 ; i < d1 ; i++)
    {
      j = 0 ;
      xx[d * (j + i) + j + i - 1] = i * (d1 - i) ;
    }
  for (i = 1 ; i < d2 ; i++)
    {
      j = d1 ;
      xx[d * (j + i) + j + i - 1] = i * (d2 - i) ;
    }
  for (i = 1 ; i < d3  ; i++)
    {
      j = d1 + d2 ;
      xx[d * (j + i) + j + i - 1] = i * (d3 - i) ;
    }
  for (i = 1 ; i < d4  ; i++)
    {
      j = d1 + d2 + d3 ;
      xx[d * (j + i) + j + i - 1] = i * (d4 - i) ;
    }
  mxSet (muE, xx) ;
  if (kas->show) mxShow(muE) ;

  /* even lowering operator F */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 1 ; i < d1 ; i++)
    {
      j = 0 ;
      xx[d * (j + i - 1) + j + i] = 1 ;
    }
  for (i = 1 ; i < d2 ; i++)
    {
      j = d1 ;
      xx[d * (j + i - 1) + j + i] = 1 ;
    }
  for (i = 1 ; i < d3  ; i++)
    {
      j = d1 + d2 ;
      xx[d * (j + i - 1) + j + i] = 1 ;
    }
  for (i = 1 ; i < d4  ; i++)
    {
      j = d1 + d2 + d3 + 0 ;
      xx[d * (j + i - 1) + j + i] = 1 ;
    }
  mxSet (muF, xx) ;
  if (kas->show) mxShow(muF) ;
  
  /* odd raising operator */  
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    {

      j = d1 ;
      xx[d * (j + i) + i] = s*b ;
    }
  for (i = 0 ; i < d3 ; i++)
    {
      j = d1 + d2 ;
      xx[d * (j + i ) + i + 1 ] = s * ( a+1 -b) ;
      j = d1 + d2 + d3 ;
      xx[d * (j+i) + i+d1+d2] = b ;
    }
  for (i = 0 ; i < d4 ; i++)
    {
      j = d1 + d2 + d3 ;
      xx[d * (j + i ) + d1 + i + 1 ] = b - s ;
    }
  mxSet (muU, xx) ;
  if (kas->show) mxShow(muU) ;

  /* b part of odd raising operator */  
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    {

      j = d1 ;
      xx[d * (j + i) + i] = s*1 ;
    }
  for (i = 0 ; i < d3 ; i++)
    {
      j = d1 + d2 ;
      xx[d * (j + i ) + i + 1 ] = s * ( -1) ;
      j = d1 + d2 + d3 ;
      xx[d * (j+i) + i+d1+d2] = 1 ;
    }
  for (i = 0 ; i < d4 ; i++)
    {
      j = d1 + d2 + d3 ;
      xx[d * (j + i ) + d1 + i + 1 ] = 1  ;
    }
  mxSet (muUb, xx) ;
  if (kas->show) mxShow(muUb) ;

  /* odd lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    {
      j = 0 ;
      xx[d * (j+i) + i+d1] = s - i ;
    }
  for (i = 0 ; i < d4 ; i++)
    {
      j = d1 ;
      xx[d * (j+i+1) + i+d1+d2+d3] = s * (i+1) ;
    }
  for (i = 0 ; i < d3 ; i++)
    {
      j = 0 ;
      xx[d * (j+i+1) + i+d1+d2] = - (i+1) ;
      j = d1 + d2 ;
      xx[d * (j+i) + i+d1+d2+d3] = s * ( a - i) ;
    }
  mxSet (muV, xx) ;
  if (kas->show) mxShow(muV) ;
  if (0) exit (0) ;

  /* odd other raising operator */
  muW = KasCommut (muE, muU, -1, kas) ;
  muW->name = "muW" ;
  if (kas->show) mxShow(muW) ;
  muWb = KasCommut (muE, muUb, -1, kas) ;
  muWb->name = "muWb" ;
  if (kas->show) mxShow(muWb) ;
  
  /* odd other oweringing operator */
  muX = KasCommut (muV, muF, -1, kas) ;
  if (kas->show) mxShow(muX) ;


  mu[0] = muY ; mu[1] = muE ; mu[2] = muF ; mu[3] = muH ;
  mu[4] = muW ; mu[5] = muX ;  mu[6] = muU ; mu[7] = muV ; 
  mu[8] = muK1 ; mu[9] = muK2 ;
  mu[10] = muUb ; mu[11] = muWb ;

  return mu ;
} /* KasimirConstructTypicMatrices */

/***********************************************************************************************************************************************/

static void KasimirCheckSuperTrace (KAS *kas)
{
  int i, j, ii ;
  int d = kas->d ;
  int d1 = kas->d1 ;
  int d2 = kas->d2 ;
  int d3 = kas->d3 ;
  int d4 = kas->d4 ;

  for (ii = 0 ; ii < 10 ; ii++)
    {
      int n = 0 ;
      MX mu = kas->mu[ii] ;
      const int *xx ;

      if (! mu)
	continue ;
      mxValues (mu, &xx, 0, 0) ;
      for (i = 0 ; i < d1 ; i++)
	n +=  xx[d*i + i] ;
      for (i = 0 ; i < d2 ; i++)
	{
	  j = d1 + i ;
	  n -=  xx[d*j + j] ;
	}
      for (i = 0 ; i < d3 ; i++)
	{
	  j = d1 + d2 + i ;
	  n -=  xx[d*j + j] ;
	}
      for (i = 0 ; i < d4 ; i++)
	{
	  j = d1 + d2 + d3 + i ;
	  n +=  xx[d*j + j] ;
	}
      n *= kas->chi ;
      if (n)
	messcrash ("Str(%s) = %d\n", mu->name, n) ;
    }
  return ;
} /* KasimirCheckSuperTrace */

/***********************************************************************************************************************************************/

static MX KasCommut (MX a, MX b, int sign, KAS *kas)
{
  MX p = mxMatMult (a, b, kas->h) ;
  MX q = mxMatMult (b, a, kas->h) ;
  MX r = mxCreate (kas->h, "r", a->type, kas->d, kas->d, 0) ;

  r = sign == 1 ? mxAdd (r, p, q, kas->h) : mxSubstract (p, q, kas->h) ;
  
  return r ;
}

/***********************************************************************************************************************************************/

static MX KasCheck (LC *up, KAS *kas)
{
  int d = kas->d ;
  int dd = kas->d * kas->d ;
  MX a = kas->mu[up->a] ;
  MX b = kas->mu[up->b] ;
  MX c = kas->mu[up->c] ;
  MX r = mxCreate (kas->h, "r", MX_INT, d,d,0) ;
  MX s = mxCreate (kas->h, "s", MX_INT, d,d,0) ;
  MX t = mxCreate (kas->h, "t", MX_INT, d,d,0) ;
  const int *xx ;
  int yy [dd] ;
  int i, k ;

  MX ab = KasCommut (a, b, up->s, kas) ;
  int scale = (up->s == 1 ? (kas->scale ? kas->scale : 1) : 1) ;

  mxValues (c, &xx, 0, 0) ;
  for (i = 0 ; i < dd ; i++)
    yy[i] = up->n * scale * xx[i] ;
  mxSet (s, yy) ;
  t = mxSubstract (ab, s, kas->h) ;
  mxValues (t, &xx, 0, 0) ;
  for (i = k = 0 ; i < dd ; i++)
    k += xx[i] * xx[i] ;
  if (k > 0)
    {
      mxShow (ab) ;
      mxShow (c) ;
      messcrash ("\nKasChect Failed [%s,%s] = %d * %s\n", a->name, b->name, up->n, c->name) ;
    }
  else if (0)
    {
      if (up->s == -1 && up->n == 0)
	printf ("[%s,%s] = 0\n",  a->name, b->name) ;
      else if (up->s == 1 && up->n == 0)
	printf ("{%s,%s} = 0\n",  a->name, b->name) ;
      else if (up->s == -1 && up->n == 1)
	printf ("[%s,%s] = %s\n",  a->name, b->name, c->name) ;
      else if (up->s == -1 && up->n == 1)
	printf ("[%s,%s] = - %s\n",  a->name, b->name, c->name) ;
      else if (up->s == -1)
	printf ("[%s,%s] = %d %s\n",  a->name, b->name, up->n, c->name) ;

      else if (up->s == 1 && up->n == 1)
	printf ("{%s,%s} = %s\n",  a->name, b->name, c->name) ;
      else if (up->s == 1 && up->n == 1)
	printf ("{%s,%s} = - %s\n",  a->name, b->name, c->name) ;
      else if (up->s == 1)
	printf ("{%s,%s} = %d %s\n",  a->name, b->name, up->n, c->name) ;
    }
  return r ;
}

/***********************************************************************************************************************************************/

static MX KasCheckR16 (KAS *kas, MX a, MX b, MX c, int scale, int sign)
{
  int d = kas->d ;
  int dd = kas->d * kas->d ;
  MX r = mxCreate (kas->h, "r", MX_COMPLEX, d,d,0) ;
  MX s = mxCreate (kas->h, "s", MX_COMPLEX, d,d,0) ;
  MX t = mxCreate (kas->h, "t", MX_COMPLEX, d,d,0) ;
  const float complex  *xx ;
  float complex yy [dd] ;
  int i, k ;

  MX ab = KasCommut (a, b, sign, kas) ;

  mxValues (c, 0, 0, &xx) ;
  for (i = 0 ; i < dd ; i++)
    yy[i] = scale * xx[i] ;
  mxSet (s, yy) ;
  if (kas->xiPrime)
    s = mxMatMult (kas->chi16, s, kas->h) ;
  t = mxSubstract (ab, s, kas->h) ;
  mxValues (t, 0, 0, &xx) ;
  for (i = k = 0 ; i < dd ; i++)
    {
      float complex z = xx[i] * xx[i] ;
      float y = creal(z)*creal(z) + cimag(z)*cimag(z) ;
      if (y > 1/100.0)
	{
	  mxNiceShow (ab) ;
	  mxNiceShow (s) ;
	  messcrash ("\nKasCheckR12 Failed [%s,%s] = %d * %s\n", a->name, b->name, scale, c->name) ;
	}
    }
  return r ;
}

/***********************************************************************************************************************************************/

static void KasimirCheckCommutators (KAS *kas)
{
  LC *up, *XXX ;
  LC Su2XXX[] = {
    {1,2,3,1,-1},
    {3,1,1,2,-1},
    {3,2,2,-2,-1},
    {0,0,0,0,0}
  } ;

  LC Su21XXX[] = {
    {1,2,3,1,-1},
    {3,1,1,2,-1},
    {3,2,2,-2,-1},

    {0,1,1,0,-1},
    {0,2,2,0,-1},
    {0,3,3,0,-1},
    
    {0,4,4,1,-1},
    {0,5,5,-1,-1},
    {0,6,6,1,-1},
    {0,7,7,-1,-1},
    
    {3,4,4,1,-1},
    {3,5,5,-1,-1},
    {3,6,6,-1,-1},
    {3,7,7,1,-1},

    {1,4,4,0,-1},
    {1,5,7,-1,-1},
    {1,6,4,1,-1},
    {1,7,7,0,-1},

    {2,4,6,1,-1},
    {2,5,5,0,-1},
    {2,6,6,0,-1},
    {2,7,5,-1,-1},

    {4,4,9,0,1},
    {4,5,9,1,1},
    {4,6,1,0,1},
    {4,7,1,-1,1},

    {5,4,9,1,1},
    {5,5,9,0,1},
    {5,6,2,-1,1},
    {5,7,2,0,1},


    {6,4,1,0,1},
    {6,5,2,-1,1},
    {6,6,8,0,1},
    {6,7,8,1,1},

    {7,4,1,-1,1},
    {7,5,2,0,1},
    {7,6,8,1,1},
    {7,7,8,0,1},
    
    {0,0,0,0,0}
  } ;

  LC OSp21XXX[] = {
    {1,2,3,1,-1},
    {3,1,1,2,-1},
    {3,2,2,-2,-1},

    {3,4,4,1,-1},
    {3,5,5,-1,-1},
    {1,4,4,0,-1},
    {1,5,4,1,-1},
    {2,4,5,1,-1},
    {2,5,5,0,-1},

    {4,4,1,2,1},
    {5,5,2,-2,1},
    {4,5,3,-1,1},

    {0,0,0,0,0}
  } ;

  /* check that K1 = Y+H */
  if (! kas->isOSp)
    {
      int i, d = kas->d ;
      const int *xxY = messalloc (d*d*sizeof(int)) ;
      const int *xxH = messalloc (d*d*sizeof(int)) ;
      const int *xxK1 = messalloc (d*d*sizeof(int)) ;
      const int *xxK2 = messalloc (d*d*sizeof(int)) ;
      
      mxValues (kas->mu[0], &xxY, 0, 0) ;
      mxValues (kas->mu[3], &xxH, 0, 0) ;
      mxValues (kas->mu[8], &xxK1, 0, 0) ;
      mxValues (kas->mu[9], &xxK2, 0, 0) ;
	    
      for (i = 0 ; i < d * d ; i++)
	if (2*xxK1[i] !=  xxY[i] + xxH[i])
	  messcrash ("K1=(Y+H)/2 failed for i=%d\n",i) ;
      for (i = 0 ; i < d * d ; i++)
	if (2*xxK2[i] !=  xxY[i] - xxH[i])
	  messcrash ("K2=(Y-H)/2 failed for i=%d\n",i) ;
    }

  XXX = kas->isSU2 ? Su2XXX : (kas->isOSp ? OSp21XXX : Su21XXX) ;

  
  for (up = XXX ; up->s ; up++)
    KasCheck (up, kas) ;
  printf ("SUCCESS (a=%d, 0) all comutators have been verified\n", kas->a) ;
  return ;
} /* KasimirCheckCommutators */

/***********************************************************************************************************************************************/
/* compute the identities essociated with the exponentiation of the supergroup */
static int myGorelikTrace (KAS *kas, MX a)
{
  const int *xx ;
  int tr = 0 ;
  int d = kas->d ;

  mxValues (a, &xx, 0, 0) ;
  for (int k = 0 ; k < d ; k++)
    tr +=  xx[d*k + k] ;
  return tr ;
}
static int GorelikTrace (KAS *kas, int i, int j, int k, int l,  AC_HANDLE h)
{
  MX a = kas->mu[i] ;
  MX b = kas->mu[j] ;
  MX u = kas->mu[k] ;
  MX v = kas->mu[l] ;
  MX c1 = mxMatMult (a, b, h) ;
  MX c2 = mxMatMult (u, v, h) ;	    
  MX c4 = mxMatMult (c1, c2, h) ;

  int t1 = myGorelikTrace (kas, c1) ;
  int t2 = myGorelikTrace (kas, c2) ;
  int t4 = myGorelikTrace (kas, c4) ;
  int t6 = (8*t4 - t1*t2) ;
  
  printf ("...... t1=%d t2=%d    t4=%d  t6=%d\n",t1,t2,t4, t6) ;
  return t6 ;
}
static void  SuperGroup (KAS *kas)
{
  AC_HANDLE h = ac_new_handle () ;
  const int *xx ;
  int nY = 0 ;
  int d = kas->d ;
  int d1 = kas->d1 ;
  int d2 = kas->d2 ;
  int d3 = kas->d3 ;
  int s = kas->scale ;
  int NN = kas->NN ;
  MX Y = kas->mu[0] ;
  static int pass = 0 ;
  int ok = 0, ok4 = 0 ;

  if (pass++)
    {
      mxValues (Y, &xx, 0, 0) ;
      for (int k = 0 ; k < d ; k++)
	nY +=  xx[d*k+ k]  ;
      if (NN > 1)
	nY /= NN ;
      /* Permanent 2 */
      for (int i = 4 ; i < 8 ; i++)
	for (int j = i ; j < 8 ; j++)
	  {
	    
	    int tr= 0,  str =0 ;
	    MX a = kas->mu[i] ;
	    MX b = kas->mu[j] ;
	    int NN = kas->NN ;
	    MX c = mxMatMult (a, b, h) ;
	    
	    mxValues (c, &xx, 0, 0) ;
	    for (int k = 0 ; k < d ; k++)
	      tr +=  xx[d*k + k] ;
	    for (int k = 0 ; k < d ; k++)
	      {
		tr +=  xx[d*k + k] ;
		str += (k < d1 || k >= d1 + d2 + d3 ? xx[d*k + k] : - xx[d*k + k]) ;
	      }
	    str *= kas->chi ;
	    if (s > 1)
	      { tr /= s ; str /= s; }
	    if (NN > 1)
	      { tr /= NN ; str /= NN ; }
	    
	    if (! tr && ! str)
	      continue ;
	    if (2*tr != nY)
	      printf ("=== ERROR Tr(%d,%d)=%d   STr()=%d Tr(Y)=%d\n",i,j,tr,str,nY) ;
	    else
	      {
		ok++ ;
		printf ("=== SUCESS Tr(%d,%d)=%d   STr()=%d Tr(Y)=%d\n",i,j,tr,str,nY) ;
	      }
	  }
      if (1)
	{
	  printf ("=== Tr(ij+ji) == Tr(Y) in %d cases a=%d b=%d\n", ok, kas->a, kas->b) ;
	  /* Gorelik Trace 4 */
	  int i=4, j=5, k=6, l=7 ;
	  int t1 = GorelikTrace (kas, i,j,k,l, h) ;
	  int t2 = GorelikTrace (kas, i,j,l,k,h) ;
	  int t3 = GorelikTrace (kas, i,k,j,l,h) ;
	  int t4 = GorelikTrace (kas, i,k,l,j,h) ;
	  int t5 = GorelikTrace (kas, i,l,j,k,h) ;
	  int t6 = GorelikTrace (kas, i,j,k,j,h) ;
	  int tr = t1 + t2 + t3 + t4 + t5 + t6 -2 ;
	  int str = 0 ;
	  int NN = kas->NN ;
	  int s = kas->scale ;
	  
	  if (s > 1)
	    { tr /= s*s ; str /= s*s; }
	  if (NN > 1)
	    { tr /= NN ; str /=  NN ; }
	  
	  if (2*tr != nY)
	    printf ("=== ERROR Tr(%d %d %d %d)=%d   STr()=%d Tr(Y)=%d\n",i,j,k,l,tr,str,nY) ;
	  else
	    {
	      ok4++ ;
	      printf ("=== SUCESS Tr(%d %d %d %d)=%d   STr()=%d Tr(Y)=%d\n",i,j,k,l,tr,str,nY) ;
	    }
	
	  printf ("=== Tr(ij+ji) == Tr(Y) in %d cases a=%d b=%d\n", ok, kas->a, kas->b) ;
	  printf ("=== Tr4(ijkl) == Tr(Y) in %d cases a=%d b=%d\n", ok4, kas->a, kas->b) ;
	}
    }
  ac_free (h) ;
} /* SuperGroup */

/***********************************************************************************************************************************************/

static void PolyDeterminant (PMX uvexp)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn[4] ;
  int N = uvexp->N ;
  PMX zz = 0 ;
  POLYNOME dd = 0 ;
  char *title = 0 ;
  printf ("Poly determinant\n") ;
  for (int i = 0 ; i < 4 ; i++)
    {
      nn[0] = i ;
      for (int j = 0 ; j < 4 ; j++)
	{
	  if ((j-i) == 0)
	    continue ;
	  nn[1] = j ; 
	  for (int k = 0 ; k < 4 ; k++)
	    {
	      if ((k-i) * (k-j) == 0)
		continue ;
	      nn[2] = k ;
	      for (int l = 0 ; l < 4 ; l++)
		{
		  if ((l-i) * (l-j) * (l-k)  == 0)
		    continue ;
		  nn[3] = l ;
		  title = hprintf (h, "...Poly determinant [%d%d%d%d]",i,j,k,l) ;
		  zz = pmxCopy (uvexp, title, h) ;
		  for (int m = 0 ; m < 4 ; m++)
		    for (int n = 0 ; n < 4 ; n++)
		      {
			/* we can exchange the lines but NOT the columns because we compute 
			 * the determinant by expanding column by column
			 */
			zz->pp[N*m + n] = polCopy (uvexp->pp[N*nn[m] + n], h) ;
		      }
		  printf ("...... polyDeterminant[%s]\n", title) ;
		  dd = pmxDeterminant (zz, h) ;
		  showPol (dd) ;
		}
	    }
	}
    }
  
  ac_free (h) ;
} /* PolyDeterminant */

/***********************************************************************************************************************************************/

static void  SuperGroupExpMap (KAS *kas)
{
  AC_HANDLE h = ac_new_handle () ;
  int d = kas->d ;
  int b = kas->b ;
  
   /****** U ******/
  /* grab the matrix */
  MX muU1 = kas->mu[6] ;
  MX muUb = kas->mu[10] ;
  const int *xxU1 ;
  const int *xxUb ;
  mxValues (muU1, &xxU1, 0, 0) ;
  mxValues (muUb, &xxUb, 0, 0) ;

  mxShow (muU1) ;
  mxShow (muUb) ;
  /* create a Polynome Matrix */
  complex zU1[d*d+1] ;
  complex zUb[d*d+1] ;
  zU1[d*d] = -1 ;  zUb[d*d] = -1 ; /* size check */

  for (int i = 0 ; i < d*d ; i++)
    {
      zU1[i] = xxU1[i] - b * xxUb[i] ;
      zUb[i] = xxUb[i] ;
    }
  PMX pU1 = pmxCreate (d, "pU1", h) ;
  PMX pUb = pmxCreate (d, "pUb", h) ;
  POLYNOME pu1 = newTheta ("u", h) ;
  POLYNOME pub = newSymbol ("b", h) ;
  pu1->tt.theta[0] = 'u' ;
  pub->tt.theta[0] = 'u' ;
  pmxSet (pU1, pu1, zU1) ;
  pmxSet (pUb, pub, zUb) ;
  pmxShow (pU1) ;
  pmxShow (pUb) ;
  PMX U = pmxSum (pU1, pUb, "U", h) ;
  pmxShow (U) ;
  

   /****** W ******/
  /* grab the matrix */
  MX muW1 = kas->mu[4] ;
  MX muWb = kas->mu[11] ;
  const int *xxW1 ;
  const int *xxWb ;
  mxValues (muW1, &xxW1, 0, 0) ;
  mxValues (muWb, &xxWb, 0, 0) ;

  mxShow (muW1) ;
  mxShow (muWb) ;
  /* create a Polynome Matrix */
  complex zW1[d*d+1] ;
  complex zWb[d*d+1] ;
  zW1[d*d] = -1 ;  zWb[d*d] = -1 ; /* size check */

  for (int i = 0 ; i < d*d ; i++)
    {
      zW1[i] = xxW1[i] - b * xxWb[i] ;
      zWb[i] = xxWb[i] ;
    }
  PMX pW1 = pmxCreate (d, "pW1", h) ;
  PMX pWb = pmxCreate (d, "pWb", h) ;
  POLYNOME pw1 = newTheta ("w", h) ;
  POLYNOME pwb = newSymbol ("b", h) ;
  pw1->tt.theta[0] = 'w' ;
  pwb->tt.theta[0] = 'w' ;
  pmxSet (pW1, pw1, zW1) ;
  pmxSet (pWb, pwb, zWb) ;
  pmxShow (pW1) ;
  pmxShow (pWb) ;
  PMX W = pmxSum (pW1, pWb, "W", h) ;
  pmxShow (W) ;
  

   /****** V ******/
  /* grab the matrix */
  MX muV1 = kas->mu[7] ;
  const int *xxV1 ;
  mxValues (muV1, &xxV1, 0, 0) ;

  mxShow (muV1) ;
  /* create a Polynome Matrix */
  complex zV1[d*d+1] ;
  zV1[d*d] = -1 ; /* size check */

  for (int i = 0 ; i < d*d ; i++)
    {
      zV1[i] = xxV1[i] ;
    }
  PMX V = pmxCreate (d, "V", h) ;
  POLYNOME pv1 = newTheta ("v", h) ;
  pv1->tt.theta[0] = 'v' ;
  pmxSet (V, pv1, zV1) ;
  pmxShow (V) ;
  
   /****** X ******/
  /* grab the matrix */
  MX muX1 = kas->mu[5] ;
  const int *xxX1 ;
  
  mxValues (muX1, &xxX1, 0, 0) ;

  mxShow (muX1) ;
  /* create a Polynome Matrix */
  complex zX1[d*d+1] ;
  zX1[d*d] = -1 ; /* size check */

  for (int i = 0 ; i < d*d ; i++)
    {
      zX1[i] = xxX1[i] ;
    }
  PMX X = pmxCreate (d, "X", h) ;
  POLYNOME px1 = newTheta ("x", h) ;
  px1->tt.theta[0] = 'x' ;
  pmxSet (X, px1, zX1) ;
  pmxShow (X) ;
  
  /****** det exp UVWX ******/

  PMX uvwxSet[] = {U, V, W, X, 0} ;
  PMX uvwx = pmxMultiSum (uvwxSet, "u+v+w+x", h) ;
  pmxShow (uvwx) ;

  PMX uvexp = pmxExponential (uvwx, "exp(u+v+w+x)", 6, h) ;
  pmxShow (uvexp) ;


  printf ("Matrix polyodering determinant\n") ;
  PolyDeterminant (uvexp) ;
  printf ("Matrix determinant\n") ;
  POLYNOME dd = pmxDeterminant (uvexp, h) ;
  showPol (dd) ;
      

  exit (0) ;

  ac_free (h) ;
} /* SuperGroupExpMap */

/***********************************************************************************************************************************************/

static void  KasimirLowerMetric (KAS *kas)
{
  MX gg ;
  int i, j, k, k1 ;
  float  n, yy[100], zz, zscale ;
  static  float yyAdjoint[100] ;
  AC_HANDLE h = ac_new_handle () ;
  static BOOL firstPass = TRUE ;
  BOOL isAdjoint = (firstPass && kas->NN == 0 && kas->a == 1 && kas->b == 1) ? TRUE : FALSE ;
 
  firstPass = FALSE ;
  
  gg = kas->gg = mxCreate (kas->h,  "gg", MX_FLOAT, 10, 10, 0) ;

  printf ("Metric gg:: ") ;
  memset (yy, 0, sizeof (yy)) ;
  for (i = 0 ; i < 8 ; i++)
    for (j = 0 ; j < 8 ; j++)
      {
	int d = kas->d ;
	int d1 = kas->d1 ;
	int d2 = kas->d2 ;
	int d3 = kas->d3 ;
	int s = kas->scale ;
	MX a = kas->mu[i] ;
	MX b = kas->mu[j] ;
	int NN = kas->NN ;
	
	if (!a || !b) continue ;
	MX c = mxCreate (h, "c", MX_INT, d, d, 0) ;
	const int *xx ;
	
	c = mxMatMult (a, b, h) ;
	mxValues (c, &xx, 0, 0) ;
	n = 0 ;
	for (k = 0 ; k < d ; k++)
	  {
	    k1 = NN ? k % (d/NN) : k ;
	    n += (k1 < d1 || k1 >= d1 + d2 + d3 ? xx[d*k + k] : - xx[d*k + k]) ;
	  }
	n *= kas->chi ;
	if (s > 1 && i>=4 && j>=4)
	  n /= s ;
	if (NN)
	  n /= NN ;
	yy [10*i + j] = n/2.0 ;
	if (isAdjoint)
	  yyAdjoint [10*i + j] = n/2.0;
	if (n != 0)
	  printf (" %d:%d=%.2f ",i,j,n/2.0) ;
      }
  printf ("\n") ;
  if (! isAdjoint && ! kas->isCycle)
    {
      float z0 = yyAdjoint[0] ;
      float a = kas->a, b = kas->b ;
      float alpha_adjoint = -2 ;
      float alpha ;
      int N = kas->NN ;
      
      if (b == 0 && N == 0)
	alpha = a * (a+1)/2.0 ;
      else if (b == a+1 && N == 0)
	alpha =  - b * (b+1)/2.0 ;
      else
	alpha = - (a+1) ;
      
      if (z0 != -alpha_adjoint)
	{
	  printf ("ERROR in lower metric adjoint g_yy = %.2f, expected %.2f\n", z0, alpha_adjoint) ;
	  exit (1) ;
	}
      if (yy[0] != -alpha)
	{
	  printf ("ERROR in lower metric (a=%d,b=%d)  g_yy = %.2f, expected %.2f\n", kas->a,kas->b, yy[0], -alpha) ;
	  exit (1) ;
	}
      if (yy[33] != alpha)
	{
	  printf ("ERROR in lower metric (a=%d,b=%d)  g_33 = %.2f, expected %.2f\n", kas->a,kas->b, yy[0], alpha) ;
	  exit (1) ;
	}
      
      zscale = yy[0]/z0 ;
      
      
      for (i = 0 ; i < 8 ; i++)
	for (j = 0 ; j < 8 ; j++)
	  {
	    zz = yy[10*i + j] ;
	    z0 = yyAdjoint[10*i + j] ;
	    if (zz != zscale * z0)
	      {
		printf ("ERROR in lower metric non uniform scale at i=%d j=%d  zz=%g z0=%g zscale=%g\n", i,j,zz,z0,zscale) ;
		exit (1) ;
	      }
	  }
      printf ("SUCCESS all lower metric entries scale up relative to the adjoint by a factor %g\n", zscale) ;
    }

  mxSet (gg, yyAdjoint) ;
  ac_free (h) ;
  return  ;
} /* KasimirMetric */

/***********************************************************************************************************************************************/

static void KasimirUpperMetric (KAS *kas)
{
  MX gg, GG ;
  int i, j ;
  const float *xx ;
  float yy[100] ;

  gg = kas->gg ;
  GG = kas->GG = mxCreate (kas->h,  "gg", MX_FLOAT, 10, 10, 0) ;

  mxValues (gg,0, &xx, 0) ;
  memcpy (yy, xx, sizeof (yy)) ;
 
  printf ("Metric GG:: ") ;
  for (i = 0 ; i < 8 ; i++)
    for (j = 0 ; j < 8 ; j++)
      {
	float z = yy[10*i + j] ;
	if (i>=4 && i < 8)
	  z *= -1 ;
	yy[10*i + j] = z ? 1/z : 0 ;
	if (z)
	  printf (" %d:%d=%.2f",i,j,1/z) ;
      }
  mxSet (GG, yy) ;
  
  printf ("\n") ;
  return ;
} /* KasimirMetric */

/***********************************************************************************************************************************************/

static void KasimirOperatorK2 (KAS *kas)
{
  int i, j, k ;
  int d = kas->d ;
  int a = kas->a ;
  int b = kas->b ;
  AC_HANDLE h = ac_new_handle () ;
  const float *xx ;
  const int *yy ;
  float zz [d*d], dz ;
  int s = kas->scale ;
  
  memset (zz, 0, sizeof (zz)) ;
  mxValues (kas->GG, 0, &xx, 0) ;
  for (i = 0 ; i < 8 ; i++)
    for (j = 0 ; j < 8 ; j++)
      if (xx[10*i + j])
	{
	  MX a = kas->mu[i] ;
	  MX b = kas->mu[j] ;

	  if (!a || !b)
	    messcrash ("uninit generrator %d %d in KAS2",i,j) ;
	  MX c = mxMatMult (a, b, h) ;
	  float z = xx[10*i + j] ;

	  mxValues (c, &yy, 0, 0) ;
	  if (s > 1 && i>= 4 && j >= 4)
	    z /= s ;
	  for (k = 0 ; k < d*d ; k++)
	    zz[k] += z * yy[k]/2.0 ;
	}

  MX kas2 = kas->kas2 = mxCreate (kas->h,  "KAS2", MX_FLOAT, d, d, 0) ;
  /* compute the casimir using the fixed adjoint metric */
  dz = 2*b * (b - a - 1)/(a+1.0) ; /* natural metric STr(ab) in the same rep */
  dz = b * (b - a - 1) ; /* fixed  metric STr(ab) in the adjoint rep */
  mxSet (kas2, zz) ;
  if (kas->show && a<4) mxNiceShow (kas2) ;
  if (dz != 0)
    for (i = 0 ; i < d*d ; i++)
      zz[i] /= dz ;

  if (dz == 0 && (zz[0]*zz[0]) > 1.0/10000)
    messcrash ("ERROR, non zero atypic KAS2 %f\n", zz[0]) ;
  if (dz != 0 && ((zz[0]-1)*(zz[0]-1)) > 1/1000.0)
    messcrash ("ERROR, KAS2 != b(b-a-1) bad ratio=%f should be 1\n", zz[0]) ;


  if (dz != 0)
    printf ("SUCCESS Quadratic super-Casimir operator (a=%d,b=%d)  KAS2 = b(b-a-1) * %f\n", kas->a, kas->b, zz[0]) ;
  else
    printf ("SUCCESS Quadratic super-Casimir operator (a=%d,b=%d) ATYPIC  KAS2 = %f expected 0\n", kas->a, kas->b, zz[0]) ;


  ac_free (h) ;
  return ;
} /* KasimirOperatorK2 */

/***********************************************************************************************************************************************/

static void GhostKasimirOperatorMinus (KAS *kas)
{
  AC_HANDLE h = ac_new_handle () ;
  int d = kas->d ;
  MX u = kas->mu[4] ;
  MX v = kas->mu[5] ;
  MX w = kas->mu[6] ;
  MX x = kas->mu[7] ;

  MX uv = mxMatMult (u,v, h) ;
  MX vu = mxMatMult (v,u, h) ;
  MX wx = mxMatMult (w,x, h) ;
  MX xw = mxMatMult (x,w, h) ;

  MX p =  mxCreate (kas->h,  "p", MX_INT, d, d, 0) ;
  MX q =  mxCreate (kas->h,  "q", MX_INT, d, d, 0) ;
  mxAdd (p, uv, wx, h) ;
  mxAdd (q, vu, xw, h) ;
  MX r = mxMatMult (p,q, h) ;

  int zz [d*d], dz = kas->scale * kas->scale ;
  if (dz)
    {
      const int *xx1 = messalloc (d*d*sizeof(int)) ;
      mxValues (r, &xx1, 0, 0) ;
      memset (zz, 0, sizeof (zz)) ;
      int i ;
      for (i = 0 ; i < d*d ; i++)
	zz[i] = xx1[i]/dz ;
      mxSet (r, zz) ;
    }
  printf( "Ghost Casimir Minus\n") ;
  mxNiceShow(r) ;
  
  ac_free (h) ;
  return ;
} /* GhostKasimirOperatorMinus */

/***********************************************************************************************************************************************/

  static void GhostKasimirOperatorXtilde2 (KAS *kas)
{
  int i, j, k, l, m1 ;
  int d = kas->d ;

  return ;
  AC_HANDLE h = ac_new_handle () ;
  const float *xx ;
  const int *yy ;
  float zz [d*d], dz ;
  MX kas2 = kas->CHI = mxCreate (kas->h,  "Ghost-Casimir2", MX_FLOAT, d, d, 0) ;
  BOOL isAdjoint = (kas->NN == 0 && kas->a == 1 && kas->b == 0) ? TRUE : FALSE ;
  
  memset (zz, 0, sizeof (zz)) ;
  /* mxValues (kas->GG, 0, &xx, 0) ; */
  for (i = 4 ; i < 8 ; i++)
    for (j = 4 ; j < 8 ; j++)
      if (1 || xx[10*i + j])
	{
	  MX a = kas->mu[i] ;
	  MX b = kas->mu[j] ;
	  float z = 0 ;
	  BOOL ok = FALSE ;
	  
	  if (!a || !b)
	    continue ;
	  MX c = mxMatMult (a, b, h) ;
	  /* float z = xx[10*i + j] ; */
	  /* if (z>0)z=1;else z=-1; */
	  if (i == 4 && j == 5) z = 1 ;
	  else if (i == 5 && j == 4) z = -1 ;
	  else if (i == 6 && j == 7) z = 1 ;
	  else if (i == 7 && j == 6) z = -1 ;
	  else if (0) continue ;
	  if (kas->scale) z /= kas->scale ;

	  mxValues (c, &yy, 0, 0) ;
	  for (k = 0 ; k < d*d ; k++)
	    {
	      zz[k] -= z * yy[k] ;
	      if (yy[k] * yy[k] > 0)
		ok = TRUE ;
	    }
	  if (0 && ok)
	    {
	      printf ("************************** X2 i = %d j = %d sign=%.2f\n", i, j, -z) ;
	      mxNiceShow (c) ;
	    }
	}

  mxSet (kas2, zz) ;
  if (0) mxNiceShow (kas2) ;
  if (0) memset (zz, 0, sizeof (zz)) ;
  for (i = 4 ; i < 8 ; i++)
    for (j = 4 ; j < 8 ; j++)
      for (k = 4 ; k < 8 ; k++)
	for (l = 4 ; l < 8 ; l++)
	  {
	    int jj = (i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l);
	    BOOL ok = FALSE ;
	    
	    if (jj)
	      {
		MX a = kas->mu[i] ;
		MX b = kas->mu[j] ;
		MX c = kas->mu[k] ;
		MX d1 = kas->mu[l] ;
		float z = 1 ;
		
		if (!a || !b || !c || !d)
		  continue ;
		MX e = mxMatMult (a, b, h) ;
		MX f = mxMatMult (e,c, h) ;
		MX g = mxMatMult (f,d1, h) ;

		if (kas->scale) z /= (kas->scale * kas->scale) ;
		mxValues (g, &yy, 0, 0) ;
		for (m1 = 0 ; m1 < d*d ; m1++)
		  {
		    zz[m1] += z * (jj>0  ? yy[m1] : -yy[m1]) ;
		    if (yy[m1] * yy[m1] > 0)
		      ok = TRUE ;
		  }
		if (0 &&  ok)
		  {
		    printf ("*** X2 i = %d j = %d  k=%d  l=%d sign=%d\n", i, j,k,l,jj> 0 ? 1 : -1) ;
		    mxNiceShow (g) ;
		  }
	      }
	  }


  for (i = 0 ; i < d*d ; i++)
    zz[i] /= 6.0 ;
  
  int a = kas->a, b = kas->b ;
  dz = b * (b - a - 1) ;
  /* dz1 = -6*(2*b -a - 1)*(2*b - 1) ; */
  if (0) mxNiceShow (kas2) ;
  mxSet (kas2, zz) ;
  if (0) mxNiceShow (kas2) ;
  if (dz != 0)
    for (i = 0 ; i < d*d ; i++)
      zz[i] /= dz ;


  if (dz == 0 && (zz[0]*zz[0]) > 1.0/10000)
    messcrash ("ERROR, non zero ghost casimir %f\n", zz[0]) ;
  if (dz != 0 && zz[0] != 1.0)
    messcrash ("ERROR, ghost casimir != b(b-a-1) bad ratio=%f should be 1\n", zz[0]) ;

  if (dz == 0)
    printf ("\nSUCCESS Ghost Casimir operator Xtilde2 (a=%d,b=%d) ATYPIC %f  expect 0\n", kas->a, kas->b, zz[0]) ;
  else
    printf ("SUCCESS Ghost-Casimir operator Xtilde2 (a=%d,b=%d) expect = b * (b - a - 1) * %.3f\n", kas->a, kas->b, zz[0]) ;

  if (kas->show && kas->a<4) mxNiceShow (kas2) ;
  if (0 && ! isAdjoint) exit (0) ;
  ac_free (h) ;
  return ;
} /* GhostKasimirOperatorXtilde2 */

/***********************************************************************************************************************************************/

static void GhostKasimirOperatorXtilde2New (KAS *kas)
{
  int i, j, k, l, m1 ;
  int d = kas->d ;

  if (! kas->show)
    return ;
  AC_HANDLE h = ac_new_handle () ;
  const float *xx ;
  const int *yy ;
  float zz [d*d], dz ;
  float zC4 = kas->zC4 ;
  const float *GG ;
  MX XT2 = kas->CHI = mxCreate (kas->h,  "Ghost-Casimir2new", MX_FLOAT, d, d, 0) ;
  BOOL isAdjoint = (kas->NN == 0 && kas->a == 1 && kas->b == 0) ? TRUE : FALSE ;
  
  zC4 = 1 ; /* fixed scale (the calculation gives -1) */

  memset (zz, 0, sizeof (zz)) ;
  mxValues (kas->GG, 0, &GG, 0) ;
  for (i = 4 ; i < 8 ; i++)
    for (j = 4 ; j < 8 ; j++)
      if (1 || xx[10*i + j])
	{
	  MX a = kas->mu[i] ;
	  MX b = kas->mu[j] ;
	  float z = 0 ;
	  BOOL ok = FALSE ;
	  z = GG[10*j + i] ;  /* contract in direct order g^{ji} i j, that is UV = VU since g_UV = g^UV = 1 */	  
	  if (z==0) continue ;
	  z = -z ;
	  if (!a || !b)
	    continue ;
	  MX c = mxMatMult (a, b, h) ;

	  if (kas->scale) z /= kas->scale ;

	  mxValues (c, &yy, 0, 0) ;
	  for (k = 0 ; k < d*d ; k++)
	    {
	      zz[k] += z * yy[k] ;
	      if (yy[k] * yy[k] > 0)
		ok = TRUE ;
	    }
	  if (0 && ok)
	    {
	      printf ("************************** X2 i = %d j = %d sign=%.2f\n", i, j, -z) ;
	      mxNiceShow (c) ;
	    }
	}

  mxSet (XT2, zz) ;
  if (1 && kas->show && kas->a<4) mxNiceShow (XT2) ;
  if (0 && kas->show) memset (zz, 0, sizeof (zz)) ;
  
  if (0)   memset (zz, 0, sizeof (zz)) ;
  for (i = 4 ; i < 8 ; i++)
    for (j = 4 ; j < 8 ; j++)
      for (k = 4 ; k < 8 ; k++)
	for (l = 4 ; l < 8 ; l++)
	  {
	    int jj = (i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l);
	    BOOL ok = FALSE ;
	    
	    if (jj)
	      {
		MX a = kas->mu[i] ;
		MX b = kas->mu[j] ;
		MX c = kas->mu[k] ;
		MX d1 = kas->mu[l] ;
		float z = 1 ;
		z = z ;
		if (!a || !b || !c || !d)
		  continue ;
		MX e = mxMatMult (a, b, h) ;
		MX f = mxMatMult (e,c, h) ;
		MX g = mxMatMult (f,d1, h) ;

		if (kas->scale) z /= (kas->scale * kas->scale) ;
		mxValues (g, &yy, 0, 0) ;
		for (m1 = 0 ; m1 < d*d ; m1++)
		  {
		    zz[m1] += z * (jj>0  ? yy[m1] : -yy[m1]) ;
		    if (yy[m1] * yy[m1] > 0)
		      ok = TRUE ;
		  }
		if (0 && kas->show && ok)
		  {
		    printf ("*** X2 i = %d j = %d  k=%d  l=%d sign=%d\n", i, j,k,l,jj> 0 ? 1 : -1) ;
		    if (0) mxNiceShow (g) ;
		  }
	      }
	  }

  mxSet (XT2, zz) ;
  if (1 && kas->show && kas->a<4) mxNiceShow (XT2) ;
  
  /* we already added the 2 terms */
  for (i = 0 ; i < d*d ; i++)
    zz[i] *= 1/6.0 ;           /* divide by 6 */
  
  int a = kas->a, b = kas->b ;
  dz = b * (b - a - 1) ;
  /* dz1 = -6*(2*b -a - 1)*(2*b - 1) ; */
  mxSet (XT2, zz) ;
  if (dz != 0)
    for (i = 0 ; i < d*d ; i++)
      zz[i] /= dz ;


  if (kas->show && kas->a<4) mxNiceShow (XT2) ;

  if (1 && dz == 0 && (zz[0]*zz[0]) > 1.0/10000)
    messcrash ("ERROR, non zero ghost casimir %f\n", zz[0]) ;
  if (1 && dz != 0 && zz[0] != 1.0)
    messcrash ("ERROR, ghost casimir != b(b-a-1) bad ratio=%f should be 1\n", zz[0]) ;

  if (dz == 0)
    printf ("\nSUCCESS Ghost Casimir operator Xtilde2New (a=%d,b=%d) ATYPIC %f  expect 0\n", kas->a, kas->b, zz[0]) ;
  else
    printf ("SUCCESS Ghost-Casimir operator Xtilde2New (a=%d,b=%d) expect = b * (b - a - 1) * %.3f zC4=%.2f\n", kas->a, kas->b, zz[0],zC4) ;


  if (0 && ! isAdjoint) exit (0) ;
  ac_free (h) ;
  return ;
} /* GhostKasimirOperatorXtilde2New */

/***********************************************************************************************************************************************/

static void GhostKasimirOperatorXtilde3 (KAS *kas)
{
  int i, j, k, l, m, m1 ;
  int d = kas->d ;

  if (! kas->show)
    return ;

  AC_HANDLE h = ac_new_handle () ;
  const float *CCC ;
  const float *C5 ;
  const int *yy ;
  float zz [d*d], dz ;
  MX XT3 = mxCreate (kas->h,  "Ghost-Chisimir3", MX_FLOAT, d, d, 0) ;
  BOOL isAdjoint = (kas->NN == 0 && kas->a == 1 && kas->b == 0) ? TRUE : FALSE ;


  if (0)
    mxValues (kas->CCC, 0, &CCC, 0) ;
  else
    mxValues (kas->CCCGhost, 0, &CCC, 0) ;
  memset (zz, 0, sizeof (zz)) ;
  /* mxValues (kas->GG, 0, &xx, 0) ; */
  for (i = 0 ; i < 8 ; i++)
    for (j = 0 ; j < 8 ; j++)
      for (k = 0 ; k < 8 ; k++)
	if (1)
	{
	  int n ;
	  MX a = kas->mu[i] ;
	  MX b = kas->mu[j] ;
	  MX c = kas->mu[k] ;
	  float z = 0 ;
	  z = CCC[100*i + 10*j + k] ;
	  if (z == 0)
	    continue ;
	  BOOL ok = FALSE ;
	  z = 3 ;
	  if (!a || !b || !c)
	    continue ;

	  n = 0 ;
	  if (i < 4) n++ ;
	  if (j < 4) n++ ;
	  if (k < 4) n++ ;
	  if (n != 1)
	    continue ;
	  if (i < 4 && j > k) z = -z ;
	  if (j < 4 && i > k) z = -z ;
	  if (k < 4 && i > j) z = -z ;
	  
	  MX e = mxMatMult (a, b, h) ;
	  MX f = mxMatMult (e, c, h) ;

	  if (kas->scale) z /= kas->scale ;
	  
	  mxValues (f, &yy, 0, 0) ;
	  if (i*j*k == 0)
	  for (n = 0 ; n < d*d ; n++)
	    {
	      zz[n] += z * yy[n] ;
	      if (yy[n] * yy[n] > 0)
		ok = TRUE ;
	    }
	  if (1 && i*j*k != 0 && ok)
	    {
	      printf ("***#######*********************** X2 i = %d j = %d k=%d sign=%.2f\n", i, j, k, z) ;
	      mxNiceShow (f) ;
	    }
	}

  mxSet (XT3, zz) ;
  printf ("#$#$#$#$#$ GHOST Chisimir 3\n") ;
  if (1) mxNiceShow (XT3) ;
  mxNiceShow (kas->CCC) ;



  
  mxValues (kas->C5, 0, &C5, 0) ;
  if (0) memset (zz, 0, sizeof (zz)) ;
  /* mxValues (kas->GG, 0, &xx, 0) ; */
  for (i = 0 ; i < 8 ; i++)
    for (j = 0 ; j < 8 ; j++)
      for (k = 0 ; k < 8 ; k++)
	for (l = 0 ; l < 8 ; l++)
	  for (m = 0 ; m < 8 ; m++)
	    if (1)
	      {
		BOOL ok = FALSE ;
		int s, n, myA = 10, myI = 0 ;
		MX a = kas->mu[i] ;
		MX b = kas->mu[j] ;
		MX c = kas->mu[k] ;
		MX e = kas->mu[l] ;
		MX f = kas->mu[m] ;
		float z = 0 ;
				
		if (!a || !b || !c || !e || !f)
		  continue ;

		/* order of even operator does not count */
		if (i < 4) i -= 100 ; 
		if (j < 4) j -= 100 ; 
		if (k < 4) k -= 100 ; 
		if (l < 4) l -= 100 ; 
		if (m < 4) m -= 100 ; 
		
		/* cut the product in 2 pieces to avoid integer overflow */
		s = (i-j)*(i-k)*(i-l)*(i-m)*(j-k) ;
		
		if (s > 0) s = 1 ;
		else if (s < 0) s = -1 ;
		else s = 0 ;
		
		s = s*(j-l)*(j-m)*(k-l)*(k-m)*(l-m) ;
		
		/* reset the even indices before issuing a continue */
		if (i < 0) i += 100 ;
		if (j < 0) j += 100 ;
		if (k < 0) k += 100 ;
		if (l < 0) l += 100 ;
		if (m < 0) m += 100 ;
		
		if (s == 0) continue ;
		if (s > 0) s = 1 ;
		else s = -1 ;

		n = 0 ;
		if (i < 4) { myA = i ; myI = 1 ; n++ ; }
		if (j < 4) { myA = j ; myI = 2 ; n++ ; }
		if (k < 4) { myA = k ; myI = 3 ; n++ ; }
		if (l < 4) { myA = l ; myI = 4 ; n++ ; }
		if (m < 4) { myA = m ; myI = 5 ; n++ ; }
		
		if (n != 1 || myA == 10)
		  continue ;
		if (myI != 1 && myI != 3 && myI != 5)
		  continue ;
		if (1 && myA != 0)
		  continue ;
		if (0 && j+k != 9) continue ;
		
		z = C5[myA] ;
		z = -1/12.0 ;
		if (z == 0)
		  continue ;
		
		MX u = mxMatMult (a, b, h) ;
		MX v = mxMatMult (u, c, h) ;
		MX w = mxMatMult (v, e, h) ;
		MX x = mxMatMult (w, f, h) ;
		
		if (kas->scale) z /= (kas->scale * kas->scale) ;
		
		mxValues (x, &yy, 0, 0) ;
		for (n = 0 ; n < d*d ; n++)
		  {
		    zz[n] += 24  * s * z * yy[n] ;
		    if (yy[n] * yy[n] > 0)
		      ok = TRUE ;
		  }
		if (ok && yy[6] > 0)
		  {
		    printf ("***#######*********************** X2 i = %d j = %d k=%d l=%d m=%d sign=%.2f\n", i, j, k, l, m, z) ;
		    mxNiceShow (x) ;
		  }
	      }

  mxSet (XT3, zz) ;
  printf ("#$#$#$#$#$ GHOST Chisimir 3\n") ;
  if (1) mxNiceShow (XT3) ;
  mxNiceShow (kas->C5) ;

  return ;
  
  if (1)
    {

      MX Y = kas->mu[0] ;
      
      MX a = kas->mu[4] ;
      MX b = kas->mu[5] ;

      MX c = kas->mu[6] ;
      MX d1 = kas->mu[7] ;

      MX ab = mxMatMult (a,b, h) ;
      MX ba = mxMatMult (b,a, h) ;

      MX cd = mxMatMult (c,d1, h) ;
      MX dc = mxMatMult (d1,c, h) ;

      const int *yyY ;
      const int *yyab ;
      const int *yyba ;
      const int *yycd ;
      const int *yydc ;

      mxValues (Y, &yyY, 0, 0) ;
      mxValues (ab, &yyab, 0, 0) ;
      mxValues (ba, &yyba, 0, 0) ;
      mxValues (cd, &yycd, 0, 0) ;
      mxValues (dc, &yydc, 0, 0) ;
		
      for (m1 = 0 ; m1 < d*d ; m1++)
	{
	  int m2 = m1 % d ;
	  zz[m1] =  yyY[m2 + d * m2] *   (yyab[m1] - yyba[m1])*( yycd[m1] - yydc[m1]) ;
	  if (kas->scale) zz[m1] /= (kas->scale * kas->scale) ;
	}
    }
  for (i = 0 ; i < 0 ; i++)
    for (j = 4 ; j < 8 ; j++)
      for (k = 4 ; k < 8 ; k++)
	for (l = 4 ; l < 8 ; l++)
	  {
	    int jj = (i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l);
	    BOOL ok = FALSE ;
	    
	    if (jj)
	      {
		MX a = kas->mu[i] ;
		MX b = kas->mu[j] ;
		MX c = kas->mu[k] ;
		MX d1 = kas->mu[l] ;
		float z = 1 ;
		
		if (!a || !b || !c || !d)
		  continue ;
		MX e = mxMatMult (a, b, h) ;
		MX f = mxMatMult (e,c, h) ;
		MX g = mxMatMult (f,d1, h) ;

		if (kas->scale) z /= (kas->scale * kas->scale) ;
		mxValues (g, &yy, 0, 0) ;
		for (m1 = 0 ; m1 < d*d ; m1++)
		  {
		    zz[m1] += z * (jj>0  ? yy[m1] : -yy[m1]) ;
		    if (m1==0 && yy[m1] * yy[m1] > 0)
		      ok = TRUE ;
		  }
		if (1 && ok)
		  {
		    printf ("*** X2 i = %d j = %d  k=%d  l=%d sign=%d\n", i, j,k,l,jj> 0 ? 1 : -1) ;
		    mxNiceShow (g) ;
		  }
	      }
	  }


  for (i = 0 ; i < d*d ; i++)
    zz[i] /= 1.0 ;
  
  int a = kas->a, b = kas->b ;
  dz = b * (b - a - 1) * (2*b - a - 1) ;
  /* dz1 = -6*(2*b -a - 1)*(2*b - 1) ; */
  mxSet (XT3, zz) ;
  if (dz != 0)
    for (i = 0 ; i < d*d ; i++)
      zz[i] /= dz ;


  if (0 && dz == 0 && (zz[0]*zz[0]) > 1.0/10000)
    messcrash ("ERROR, non zero ghost casimir %f\n", zz[0]) ;
  if (0 && dz != 0 && zz[0] != 1.0)
    messcrash ("ERROR, ghost casimir != b(b-a-1) bad ratio=%f should be 1\n", zz[0]) ;

  if (0 && dz == 0)
    printf ("\nSUCCESS Ghost Casimir operator Xtilde3 (a=%d,b=%d) ATYPIC %f  expect 0\n", kas->a, kas->b, zz[0]) ;
  else
    printf ("######QUESTION  Ghost-Casimir operator Xtilde3 (a=%d,b=%d) expect = b * (b - a - 1) * (2b - a - 1) = %d\n", kas->a, kas->b, b * (b-a-1)*(2*b-a-1)) ;

  if (kas->show && kas->a<4) mxNiceShow (XT3) ;
  if (0 && ! isAdjoint) exit (0) ;
  ac_free (h) ;
  return ;
} /* GhostKasimirOperatorXtilde3 */

/***********************************************************************************************************************************************/
/* Casimir proposed by Peter, july 28 */
static void KasimirOperatorK4 (KAS *kas)
{
  int d = kas->d ;
  AC_HANDLE h = ac_new_handle () ;
  MX U, V, W, X ;
  int ii ;
  
  for (ii = 0 ; ii < 2 ; ii++)
    {
      if (ii==0)
	{
	  U = kas->mu[6] ;
	  V = kas->mu[7] ;
	  W = kas->mu[4] ;
	  X = kas->mu[5] ;
	}
      else
	{
	  V = kas->mu[6] ;
	  U = kas->mu[7] ;
	  X = kas->mu[4] ;
	  W = kas->mu[5] ;
	}
      
      MX Y = kas->mu[0] ;
      
      MX WX = mxMatMult (W,X,h) ;
      MX UWX = mxMatMult (U,WX,h) ;
      MX WXU = mxMatMult (WX,U,h) ;
      MX uwx =  mxCreate (h, "[U,WX]", MX_INT, d, d, 0) ;
      uwx = mxSubstract (UWX, WXU, h) ;
      MX Vuwx = mxMatMult (V, uwx,h) ;
      MX uwxV = mxMatMult (uwx,V,h) ;
      MX vuwx = mxCreate (h, "{V,[U,WX]}", MX_INT, d, d, 0) ;
      vuwx = mxAdd (vuwx, Vuwx,uwxV, h) ;
      
      MX Y2 = mxMatMult (Y,Y,h) ;
      MX Y2WX = mxMatMult (Y2,WX,h) ;
      MX UY2WX = mxMatMult (U,Y2WX,h) ;
      MX Y2WXU = mxMatMult (Y2WX,U,h) ;
      MX uy2wx =  mxCreate (h, "[U,Y2WX]", MX_INT, d, d, 0) ;
      uy2wx = mxSubstract (UY2WX, Y2WXU, h) ;
      MX Vuy2wx = mxMatMult (V, uy2wx,h) ;
      MX uy2wxV = mxMatMult (uy2wx,V,h) ;
      MX vuy2wx = mxCreate (h, "{V,[U,Y2WX]}", MX_INT, d, d, 0) ;
      vuy2wx = mxAdd (vuy2wx, Vuy2wx,uy2wxV, h) ;
        
      mxNiceShow (vuwx) ;
      mxNiceShow (vuy2wx) ;
    }
  
  ac_free (h) ;
  return ;
} /* KasimirOperatorK4 */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static void  KasimirLower3tensor (KAS *kas, BOOL isGhost)
{
  int i, j, k, i1, scale ;
  float yy[1000] ;
  static  float yyAdjoint[1000] ;
  static  float yyAdjointGhost[1000] ;
  float zz, zscale = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  MX ccc ;
  int mx0 = 0 ;
  int mx1 = 8 ;
  static BOOL firstPass = TRUE ;
  static BOOL firstPassGhost = TRUE ;
  BOOL isAdjoint = (kas->NN >= 0 && kas->a == 1 && kas->b == 0) ? TRUE : FALSE ;

  if (isGhost)
    {
      if (!kas->cccGhost)
	kas->cccGhost = mxCreate (kas->h,  "cccGhost", MX_FLOAT, 10, 10, 10, 0) ;
      ccc = kas->cccGhost ;
      if (0)
	if (! isAdjoint || ! firstPassGhost)
	  goto done ;
      if (isAdjoint)
	firstPassGhost = FALSE ;
    }
  else
    {
      if (! kas->ccc)
	kas->ccc = mxCreate (kas->h,  "ccc", MX_FLOAT, 10, 10, 10, 0) ;
      ccc = kas->ccc ;
      if (0)
	if (! isAdjoint || ! firstPass)
	  goto done ;
      if (isAdjoint)
	firstPass = FALSE ;
    }
  
  printf ("Lower ccc:: ") ;

  memset (yy, 0, sizeof (yy)) ;
  for (i = mx0 ; i < mx1 ; i++)
    for (j = mx0 ; j < mx1 ; j++)
      for (k = mx0 ; k < mx1 ; k++)
      {
	int d = kas->d ;
	int d1 = kas->d1 ;
	int d2 = kas->d2 ;
	int d3 = kas->d3 ;
	int s ;
	MX a = kas->mu[i] ;
	MX b = kas->mu[j] ;
	MX c = kas->mu[k] ;
	MX u,v,x,y ;
	MX z = mxCreate (h, "z", MX_INT, d, d, 0) ;
	const int *xx ;
	float zz1 ;	
	u = mxMatMult (a, b, h) ;
	v = mxMatMult (u, c, h) ;
	x = mxMatMult (a, c, h) ;
	y = mxMatMult (x, b, h) ;

	if (1 && isGhost  && i<4 && j<4 && k<4)
	  continue ;
	if (j >= 4 && j <= 7 && k >= 4 && k <= 7)
	  s = -1 ;
	else
	  s = 1 ;
	if (i > 40)
	  s = -s ;
	if (s == -1)
	  z = mxSubstract (v, y, h) ;
	else
	  z = mxAdd (z, v, y, h) ;
	mxValues (z, &xx, 0, 0) ;
	zz = 0  ; /* STr (a[i,j]) */
	zz1 = 0 ; /*  Tr (a[i,j]) */
	for (i1 = 0 ; i1 < d ; i1++)
	  {
	    int NN = kas->NN ;
	    int dd2 = NN ? d/NN : d ;
	    int i2 = i1 % dd2 ; 
	    zz += (i2 < d1 || i2 >= d1 + d2 + d3 ? xx[d*i1 + i1] : - xx[d*i1 + i1]) ;
	    zz1 += (i2 < d1 || i2 >= d1 + d2 + d3 ? xx[d*i1 + i1] : + xx[d*i1 + i1]) ;
	  }
	zz *= kas->chi/2.0 ;

	scale = (i>=4 || j >= 4 || k >= 4) ? kas->scale : 0 ;
	if (scale != 0)
	  { zz /= scale ; zz1 /= scale ; }
	
	yy [100*i + 10*j + k] = zz ;
	if (0)
	  {
	    if ((i+j+k==0)|| (zz != 0))
	      printf ("LOWERcccSTr (%d%d%d a=%d b=%d)=%g\n ",i,j,k,kas->a,kas->b,zz) ;
	    if ((i+j+k==0)|| (zz1 != 0))
	      printf ("LOWERcccTr (%d%d%d a=%d b=%d)=%g\n ",i,j,k,kas->a,kas->b,zz1) ;
	  }
      }
  if (! isAdjoint || ! firstPassGhost)
    goto done ;
  if (!firstPass && ! isAdjoint && ! isGhost)
    {
      float z0 = yyAdjoint[0] ;
      zscale = yy[0]/z0 ;
      
      for (i = mx0 ; i < mx1 ; i++)
	for (j = mx0 ; j < mx1 ; j++)
	  for (k = mx0 ; k < mx1 ; k++)
	    
	    {
	      zz = yy[100*i + 10*j + k] ;
	      z0 = yyAdjoint[100*i + 10*j + k] ;
	      if (zz != zscale * z0)
		{
		  printf ("ERROR in lower3tensor at i=%d j=%d k=%d zz=%g z0=%g zscale=%g\n", i,j,k,zz,z0,zscale) ;
		  exit (1) ;;
		}
	    }
      printf ("SUCCESS all lower 3 tensor scale up by a factor %g\n", zscale) ;
    }

  mxNiceShow (ccc) ;
  /* the lower 3 tensor scales (a,b) relative to the lepton (a=1,b=0) by a factor s=(a+1)(2b-a-1) = (a+1)(y-1) = 1/4  Tr(Y)
   * for the quarks b=2/3,a=0  s=1/3, really -1/3 because we start on a right state, hence BIM lepton + 3 quarks = 0
   * whereas as operrators C_3(lepton)==0 (atypic) c_3(quarks) non zero
   */
  if (isGhost)
    memcpy (yyAdjointGhost, yy, sizeof (yy)) ;
  else
    memcpy (yyAdjoint, yy, sizeof (yy)) ;
 done:
  if (isGhost)
    mxSet (ccc, yyAdjointGhost) ;
  else
    mxSet (ccc, yyAdjoint) ;
 
  ac_free (h) ;
  return  ;
} /* KasimirLower3tensor */

/***********************************************************************************************************************************************/

static void  KasimirLower4tensor (KAS *kas)
{
  int i, j, k, l, i1, scale ;
  float zz, zc4 = 0 ;
  static float zc4Adjoint = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  int mx0 = 4 ;
  int mx1 = 8 ;
  static BOOL firstPass = TRUE ;
  BOOL isAdjoint = (kas->NN >= 0 && kas->a == 1 && kas->b == 0) ? TRUE : FALSE ;

  if (!firstPass || !isAdjoint)
    {
      kas->zc4 = zc4Adjoint ;
      return ;
    }
  /* we do not ned to compute c4: it is always antisymmetrized, so it is dual to a scalar zc4
     ccc = kas->c4 = mxCreate (kas->h,  "ccc", MX_FLOAT, 10, 10, 10, 10, 0) ;
  */

  firstPass = FALSE ;
  printf ("Lower c4:: ") ;

  for (i = mx0 ; i < mx1 ; i++)
    for (j = mx0 ; j < mx1 ; j++)
      for (k = mx0 ; k < mx1 ; k++)
	for (l = mx0 ; l < mx1 ; l++)
	  {
	    int d = kas->d ;
	    int d1 = kas->d1 ;
	    int d2 = kas->d2 ;
	    int d3 = kas->d3 ;
	    int s ;
	    MX a = kas->mu[i] ;
	    MX b = kas->mu[j] ;
	    MX c = kas->mu[k] ;
	    MX e = kas->mu[l] ;
	    MX u,v,w ;
	    const int *xx ;
	    
	    s = (i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l) ;
	    if (s == 0) continue ;
	    if (s > 0) s = 1 ;
	    else s = -1 ;
	    
	    u = mxMatMult (a, b, h) ;
	    v = mxMatMult (u, c, h) ;
	    w = mxMatMult (v, e, h) ;
	    
	    mxValues (w, &xx, 0, 0) ;
	    
	    /* compute the supertrace */
	    zz = 0 ;
	    for (i1 = 0 ; i1 < d ; i1++)
	      {
		int NN = kas->NN ;
		int dd2 = NN ? d/NN : d ;
		int i2 = i1 % dd2 ; 
		zz += (i2 < d1 || i2 >= d1 + d2 + d3 ? xx[d*i1 + i1] : - xx[d*i1 + i1]) ;
	      }
	    zz *= s * kas->chi ;
	    
	    scale = kas->scale * kas->scale ; /* we use 4 odd operators */
	    if (scale != 0)
	      zz /= scale ;
	    
	    /*
	      yy [1000*i + 100*j + 10*k + l] = zz ;
	      if (isAdjoint)
	      yyAdjoint [1000*i + 100*j + 10*k + l] = zz ;
	    */
	    zc4 += zz ;
	    if (zz != 0)
	      printf ("C3(%d%d%d%d)=%g ",i,j,k,l,zz) ;
	  }

  zc4Adjoint = kas->zc4 = zc4/4 ;
  
  ac_free (h) ;
  return  ;
} /* KasimirLower4tensor */

/***********************************************************************************************************************************************/

static void  KasimirLower5tensor (KAS *kas)
{
  int i, j, k, l, m, p,  i1, scale ;
  float yy[10] ;
  static  float yyAdjoint[10] ;
  float zz, zscale = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  MX c5 ;
  int mx0 = 0 ;
  int mx1 = 8 ;
  static BOOL firstPass = TRUE ;
  BOOL isAdjoint = (kas->NN >= 0 && kas->a == 1 && kas->b == 0) ? TRUE : FALSE ;

  if (! kas->c5)
    kas->c5 = mxCreate (kas->h,  "c5", MX_FLOAT, 10, 0) ;
  c5 = kas->c5 ;
  if (!firstPass || !isAdjoint)
    goto done ;

  if (isAdjoint)
    firstPass = FALSE ;
  
  printf ("Lower c5:: ") ;

  memset (yy, 0, sizeof (yy)) ;
  for (i = mx0 ; i < mx1 ; i++)
    for (j = mx0 ; j < mx1 ; j++)
      for (k = mx0 ; k < mx1 ; k++)
	for (l = mx0 ; l < mx1 ; l++)
	  for (m = mx0 ; m < mx1 ; m++)
	    {
	      int myA = 0 ;
	      int d = kas->d ;
	      int d1 = kas->d1 ;
	      int d2 = kas->d2 ;
	      int d3 = kas->d3 ;
	      int s ;
	      MX a = kas->mu[i] ;
	      MX b = kas->mu[j] ;
	      MX c = kas->mu[k] ;
	      MX e = kas->mu[l] ;
	      MX f = kas->mu[m] ;
	      MX u,v,w,x ;
	     const int *xx ;
	     
	     /* count the even operator, keep one */
	      s = 0 ;
	      if (i < 4) s++ ;
	      if (j < 4) s++ ;
	      if (k < 4) s++ ;
	      if (l < 4) s++ ;
	      if (m < 4) s++ ;
	      
	      if (s != 1)
		continue ;
	      
	      /* order of even operator does not count */
	      if (i < 4) i -= 100 ; 
	      if (j < 4) j -= 100 ; 
	      if (k < 4) k -= 100 ; 
	      if (l < 4) l -= 100 ; 
	      if (m < 4) m -= 100 ; 
	      
	      /* cut the product in 2 pieces to avoid integer overflow */
	      s = (i-j)*(i-k)*(i-l)*(i-m)*(j-k) ;

	      if (s > 0) s = 1 ;
	      else if (s < 0) s = -1 ;
	      else s = 0 ;
	      
	      s = s*(j-l)*(j-m)*(k-l)*(k-m)*(l-m) ;

	      /* reset the even indices before issuing a continue */
	      if (i < 0) i += 100 ;
	      if (j < 0) j += 100 ;
	      if (k < 0) k += 100 ;
	      if (l < 0) l += 100 ;
	      if (m < 0) m += 100 ;

	      if (s == 0) continue ;
	      if (s > 0) s = 1 ;
	      else s = -1 ;

	      if (i < 4) myA = i ;
	      if (j < 4) myA = j ;
	      if (k < 4) myA = k ;
	      if (l < 4) myA = l ;
	      if (m < 4) myA = m ;


	      
	      u = mxMatMult (a, b, h) ;
	      v = mxMatMult (u, c, h) ;
	      w = mxMatMult (v, e, h) ;
	      x = mxMatMult (w, f, h) ;
	      
	      mxValues (x, &xx, 0, 0) ;
	      
	      /* compute the supertrace */
	      zz = 0 ;
	      for (i1 = 0 ; i1 < d ; i1++)
		{
		  int NN = kas->NN ;
		  int dd2 = NN ? d/NN : d ;
		  int i2 = i1 % dd2 ; 
		  zz += (i2 < d1 || i2 >= d1 + d2 + d3 ? xx[d*i1 + i1] : - xx[d*i1 + i1]) ;
		}
	      zz *= s * kas->chi ;
	      
	      scale = kas->scale * kas->scale ; /* we use 4 odd operators */
	      if (scale != 0)
		zz /= scale ;

	      zz /= 24 ;
	      yy [myA] += zz ;
	      yyAdjoint [myA] += zz ;
	    }
  if (!firstPass && ! isAdjoint)
    {
      float z0 = yyAdjoint[0] ;
      zscale = yy[0]/z0 ;
      
      for (p = 0 ; p < 4 ; p++)
	{
	  zz = yy[p] ;
	  z0 = yyAdjoint[p] ;
	  if (zz != zscale * z0)
	    {
	      printf ("ERROR in lower5tensor at i=%d zz=%g z0=%g zscale=%g\n", i,zz,z0,zscale) ;
	      exit (1) ;;
	    }
	}
      printf ("SUCCESS all lower 5 tensor scale up by a factor %g\n", zscale) ;
    }

  /* the lower 3 tensor scales (a,b) relative to the lepton (a=1,b=0) by a factor s=(a+1)(2b-a-1)
   * for the quarks b=2/3,a=0  s=1/3, really -1/3 because we start on a right state, hence BIM lepton + 3 quarks = 0
   * whereas as operrators C_3(lepton)==0 (atypic) c_3(quarks) non zero
   */
 done:
    mxSet (c5, yyAdjoint) ;
    mxNiceShow (c5) ; 
  ac_free (h) ;
  return  ;
} /* KasimirLower5tensor */

/***********************************************************************************************************************************************/

static void  KasimirUpper3tensor (KAS *kas)
{
  int i, j, k ;
  float yy[1000] ;
  float yyGhost[1000] ;
  AC_HANDLE h = ac_new_handle () ;
  MX CCC = kas->CCC = mxCreate (kas->h,  "CCC", MX_FLOAT, 10, 10, 10, 0) ;
  MX CCCGhost = kas->CCCGhost = mxCreate (kas->h,  "CCCGhost", MX_FLOAT, 10, 10, 10, 0) ;
  const float *GG ;
  const float *ccc ;
  const float *cccGhost ;

  printf ("Upper CCC:: ") ;
  mxValues (kas->GG, 0,  &GG, 0) ;
  mxValues (kas->ccc, 0, &ccc, 0) ;
  mxValues (kas->cccGhost, 0, &cccGhost, 0) ;
  memset (yy, 0, sizeof (yy)) ;
  memset (yyGhost, 0, sizeof (yyGhost)) ;
  for (i = 0 ; i < 8 ; i++)
    for (j = 0 ; j < 8 ; j++)
      for (k = 0 ; k < 8 ; k++)
	{
	  int a, b, c ; /* dummy indices */
	  float  z = 0 ;
	  float  zGhost = 0 ;
	  for (a = 0 ; a < 8 ; a++)
	    for (b = 0 ; b < 8 ; b++)
	      for (c = 0 ; c < 8 ; c++)
		{
		       z += GG[10*i + a] * GG[10*j + b] * GG[10*k + c] * ccc[100*a + 10 * b + c] ;
		  zGhost += GG[10*i + a] * GG[10*j + b] * GG[10*k + c] * cccGhost[100*a + 10 * b + c] ;
		}

	  if (i>4 || j>4 || k>4) { z = -z ; zGhost = - zGhost ; }
	  if (kas->show && z != 0)
	    printf (" %d:%d:%d=%.2f ::ghost %.2f",i,j,k,z,zGhost) ;
	  yy[100*i + 10*j + k] += z ;
	  yyGhost[100*i + 10*j + k] += zGhost ;
	}
  mxSet (CCC, yy) ;
  mxSet (CCCGhost, yyGhost) ;
  ac_free (h) ;
  return  ;
} /* KasimirUpper3tensor */

/***********************************************************************************************************************************************/

static void  KasimirUpper4tensor (KAS *kas)
{
  if (kas->zc4)
    kas->zC4 = 1/kas->zc4 ;
  return ;
} /* KasimirUpper4tensor */

/***********************************************************************************************************************************************/
/* this is really a single index tensor a */
static void  KasimirUpper5tensor (KAS *kas)
{
  int i ;

  if (! kas->c5) return ;
  float yy[10] ;
  AC_HANDLE h = ac_new_handle () ;
  MX CCC = kas->C5 = mxCreate (kas->h,  "C5", MX_FLOAT, 10, 0) ;
  const float *GG ;
  const float *ccc ;

  printf ("Upper C5:: ") ;
  mxValues (kas->GG, 0,  &GG, 0) ;
  mxValues (kas->c5, 0, &ccc, 0) ;
  memset (yy, 0, sizeof (yy)) ;
  for (i = 0 ; i < 8 ; i++)
    {
      int a ; /* dummy indices */
      float  z = 0 ;
      for (a = 0 ; a < 8 ; a++)
	{
	  z += GG[10*i + a] * ccc[a] ;
	}
      if (kas->show && z != 0)
	printf (" %d=%.2f",i,z) ;
      yy[i] += z * kas->zC4 ;
    }
  mxSet (CCC, yy) ;
  mxNiceShow (CCC) ;
  ac_free (h) ;
  return  ;
} /* KasimirUpper5tensor */

/***********************************************************************************************************************************************/

static void KasimirUpperTensor (KAS *kas)
{
  KasimirLower3tensor (kas, FALSE) ;
  KasimirLower3tensor (kas, TRUE) ;
  KasimirLower4tensor (kas) ;
  KasimirLower5tensor (kas) ;

  KasimirUpper3tensor (kas) ;
  KasimirUpper4tensor (kas) ;
  KasimirUpper5tensor (kas) ;
  return ;
} /* KasimirUppertensor */

/***********************************************************************************************************************************************/

static void KasimirOperatorK3 (KAS *kas)
{
  int i, j, k, m ;
  int d = kas->d ;
  AC_HANDLE h = ac_new_handle () ;
  const float *xx ;
  const int *yy ;
  float z, zz [d*d] ;
  int mx0 = 0 ;
  int mx1 = 8 ;
  int s = kas->scale ;
  float zexpected ;
  int a = kas->a ;
  int b = kas->b ;
  
  memset (zz, 0, sizeof (zz)) ;
  mxValues (kas->CCC, 0, &xx, 0) ;
  for (i = mx0 ; i < mx1 ; i++)
    for (j = mx0 ; j < mx1 ; j++)
      for (k = mx0 ; k < mx1 ; k++)
	if (xx[100*i + 10*j + k])
	  {
	    MX a = kas->mu[i] ;
	    MX b = kas->mu[j] ;
	    MX c = kas->mu[k] ;
	    MX u = mxMatMult (a, b, h) ;
	    MX v = mxMatMult (u, c, h) ;
	    float z = 1, n = xx[100*i + 10*j + k] ;
	    
	    mxValues (v, &yy, 0, 0) ;
	    if (s > 1 && (i>= 4 || j>=4 || k >= 4))
	      z = 1.0/s ;
	      
	    for (m = 0 ; m < d*d ; m++)
	      zz[m] += z * n * yy[m]/6 ;
	  }

  MX kas3 = kas->kas3 = mxCreate (kas->h,  "KAS3", MX_FLOAT, d, d, 0) ;
  mxSet (kas3, zz) ;

  z = b * (b - a - 1) ;
  zexpected = 4 * (b - a)  * (b - a - 1) * (2*b - a - 1) * (2*b + a - 1)  ;
  zexpected = b * (b - a - 1) * (2*b - a - 1)  ;  /* using a fixed (adjoint) C_{abc} lifted using G^{ab} also adjoint, i.e. a fixed operator for all reps */
  if (z == 0)
    printf ("\nSUCCESS Cubic super-Casimir operator KAS3 (a=%d,b=%d) ATYPIC %f  expect 0\n", kas->a, kas->b, zz[0]) ;
  else if (0 ||  (2*zz[0] - zexpected)*(2*zz[0] - zexpected) < .1)
    printf ("\nSUCCESS Cubic super-Casimir operator KAS3 (a=%d,b=%d) = %f, zexpected= b  * (b - a - 1) * (2*b - a - 1)/2 = %f  = z * %f\n", kas->a, kas->b, zz[0] , zexpected/2, 2*zz[0]/zexpected) ;
  else
    messerror ("\nCubic super-Casimir operator KAS3 (a=%d,b=%d) z = %f expect b(b-a-1)(2b - a -1)/2 =  %f\n", kas->a, kas->b, zz[0], zexpected/2.0) ;
  
  if (kas->show && kas->a<6) mxNiceShow (kas3) ;

  ac_free (h) ;
  return ;
} /* KasimirOperatorK3 */

/***********************************************************************************************************************************************/

static void QFTscalar (KAS *kas)
{
  int i, j, k, l ;
  if (! kas->show)
    return ;

  AC_HANDLE h = ac_new_handle () ;
  int d = kas->d ;
  float zz[d*d] ;
  float scale = kas->scale ;

  if (scale == 0)
    scale = 1 ;
  memset (zz, 0, sizeof (zz)) ;
  printf ("In the 4 scalar vertex we want to compute Tr(ijkl(1+chi)/2) symmetrized in ik and jl\n") ;
  for (i = 4 ; i < 8 ; i++)
    for (j = 4 ; j < 8 ; j++)
      for (k = 4 ; k < 8 ; k++)
	for (l = 4 ; l < 8 ; l++)
	  {
	    if (k<i) continue ;
	    if (l<j) continue ;
	    MX a = kas->mu[i] ;
	    MX b = kas->mu[j] ;
	    MX c = kas->mu[k] ;
	    MX d = kas->mu[l] ;
	    MX mm1[] = {a,b,c,d,0} ;
	    MX mm2[] = {a,d,c,b,0} ;
	    MX mm3[] = {c,b,a,d,0} ;
	    MX mm4[] = {c,d,a,b,0} ;
	    MX z1 = mxMatMultiProduct (h, mm1) ;
	    MX z2 = mxMatMultiProduct (h, mm2) ;
	    MX z3 = mxMatMultiProduct (h, mm3) ;
	    MX z4 = mxMatMultiProduct (h, mm4) ;
	    MX zz ;
	    MX zz1[] = {z1,z2,z3,z4,0} ;
	    MX zz2[] = {z1,z3,0} ;
	    MX zz3[] = {z1,z2,0} ;
	    if ((i-k)*(j-l) != 0) zz = mxMultiSum (h, zz1) ;
	    else if ((i-k) != 0) zz = mxMultiSum (h, zz2) ;
	    else if ((j-l) != 0) zz = mxMultiSum (h, zz3) ;
	    else  zz = z1 ;
	    float complex z = mxMatTrace (zz) ; 
	    float y = creal(z)*creal(z) + cimag(z)*cimag(z) ;
	    if (y > 1/1000) 
	      {
		float u = creal(z), v = cimag(z) ;
		if (u*u < 0.01) u = 0 ;
		if (v*v < 0.01) v = 0 ;
		printf("%d %d %d %d -> %.2f + i %.2f\n", i,j,k,l,u, v) ;
	      }
	  }


  printf(" In the scalar psi-psi diagram g^ij (i j) 4 + 4 g^{ji}{i j} should look like 4\n") ;
  
  MX w = mxCreate (h, "wave function", MX_FLOAT, kas->d, kas->d, 0) ;
  mxValues (kas->GG, 0, &GG, 0) ; 
    for (i = 0 ; i < 8 ; i++)
      for (j = 0 ; j < 8 ; j++)
      {
	float z = GG[10*i + j] ;
	if (z)
	  {
	    int m ;
	    MX a = kas->mu[i] ;
	    MX b = kas->mu[j] ;

	    MX c = mxMatMult (a, b, h) ;
	    MX e = mxMatMult (c, K, h) ;
	    MX f = mxMatMult (K, c, h) ;
	    
	    MX g ;
	    if ( i >= 40)
	      continue ;
	    if (i >= 14)
	      g = mxSubstract (e, f,h) ;
	    else
	      g = mxAdd (0,e, f,h) ;
	    mxValues (g, &yy, 0, 0) ;
	    for (m = 0 ; m < d*d ; m++)
	      zz[m] += z * yy[m] / scale ;
	    printf ("QFT i=%d j=%d\n", i, j) ;
	    if (0)
	      {
		mxNiceShow (a) ;
		mxNiceShow (b) ;
		mxNiceShow (c) ;
		mxNiceShow (e) ;
		mxNiceShow (f) ;
	      }
	    mxNiceShow (g) ;
	  }
      }
    
    mxSet (w, zz) ;
    mxNiceShow (K) ;
    mxNiceShow (w) ;
    
    ac_free (h) ;
} /* QFTscalar */

  /***********************************************************************************************************************************************/

static void Kasimirs (int a, int b, BOOL show)
{
  KAS kas ;
  memset (&kas, 0,sizeof(KAS)) ;
  kas.h = ac_new_handle () ;
  kas.a = a ;    /* Kac Dynkin weights of the heighest weight */
  kas.b = b ;
  kas.show = show ;
  kas.isOSp = FALSE ;
  AC_HANDLE h = kas.h ;
  
  if (a>=2000)
    KasimirConstructSU2Matrices (&kas) ;
  else if (a>=1000)
    KasimirConstructOSp1_2Matrices (&kas) ;
  else if (a>0 && b == 0)
    KasimirConstructAtypicMatrices (&kas) ;
  else if (a >= 0 && b == a + 1)
    KasimirConstructAntiMatrices (&kas) ;
  else
    KasimirConstructTypicMatrices (&kas, show) ;

  KasimirCheckSuperTrace (&kas) ;
  KasimirCheckCommutators (&kas) ;

  KasimirLowerMetric (&kas) ;
  if (0 && show) exit (0) ;
  
  
  KasimirUpperMetric (&kas) ;
  KasimirUpperTensor (&kas) ;
  
  KasimirOperatorK2 (&kas) ;
  GhostKasimirOperatorXtilde2 (&kas) ;
  GhostKasimirOperatorXtilde2New (&kas) ;
  GhostKasimirOperatorMinus (&kas) ;
  
  if (0) GhostKasimirOperatorXtilde3 (&kas) ;
  if (1 && kas.show) 
    {
      QFTscalar (&kas) ;
      exit (0) ;
    }
  if (0) KasimirOperatorK4 (&kas) ;
  if (0) return ;

  MX qmuH = kas.mu[3] ;
  MX qmuX = kas.mu[6] ;
  int d = kas.d ;

  if (kas.show)
    {
      printf ("Verify that the casimir commutes with H\n") ;
      MX CKX = mxMatMult (kas.kas2, qmuH, h) ;
      MX CXK = mxMatMult (qmuH, kas.kas2, h) ;
      MX Com =  mxCreate (h, "[casimir,H]", MX_COMPLEX,d,d, 0) ;
      Com = mxSubstract (CKX, CXK, h) ;
      if (kas.show)
	mxNiceShow (Com) ;
      
      printf ("Verify that the casimir commutes with X\n") ;
      MX CKX2 = mxMatMult (kas.kas2, qmuX, h) ;
      MX CXK2 = mxMatMult (qmuX, kas.kas2, h) ;
      MX Com2 =  mxCreate (h, "[casimir,X]", MX_COMPLEX,d,d, 0) ;
      Com = mxSubstract (CKX2, CXK2, h) ;
      if (kas.show)
	mxNiceShow (Com2) ;
      
      printf ("Verify that the S-casimir anticommutes with XU and YV\n") ;
      MX SCKX = mxMatMult (kas.CHI, qmuX, h) ;
      MX SCXK = mxMatMult (qmuX, kas.CHI, h) ;
      MX SCom =  mxCreate (h, "{S-casimir,X}", MX_COMPLEX, d,d, 0) ;
      SCom = mxAdd (SCom, SCKX, SCXK, h) ;
      if (kas.show) mxNiceShow (SCom) ;
      
      printf ("Compute the square of the S-casimir\n") ;
      MX SC2 = mxMatMult (kas.CHI,kas.CHI, h) ;
      if (0) SC2 = mxLinearCombine (SC2, 1, SC2, -1, kas.kas2, h) ;
      SC2->name = "S-Casimir square" ;
      if (kas.show) mxNiceShow (SC2) ;
      
      printf ("Compute the product of the Casimir by the S-casimir Q^3 = Q\n") ;
      MX SC3 = mxMatMult (kas.kas2, kas.CHI, h) ;
      if (0) SC3 = mxLinearCombine (SC3, 1, SC3, -1, kas.CHI, h) ;
      SC3->name = "S-Casimir cube" ;
      if (kas.show) mxNiceShow (SC3) ;
    }
  KasimirUpperTensor (&kas) ;
  
  if (show)
    KasimirOperatorK3 (&kas) ;
  SuperGroup (&kas) ;
  if (show)
    SuperGroupExpMap (&kas) ;

} /* Kasimirs */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static void GhostKasimirOperatorR16 (KAS *kas)
{
  int i, j, k, l, m1 ;
  int d = kas->d ;
  AC_HANDLE h = ac_new_handle () ;
  const complex float *xx ;
  const complex float *yy ;
  complex float zz [d*d] ;
  MX XT2 = kas->CHI = mxCreate (kas->h,  "Ghost-CasimirR16", MX_COMPLEX, d, d, 0) ;

  memset (zz, 0, sizeof (zz)) ;
  if (1)
    {
      MX a = kas->Rmu[4] ;
      MX b = kas->Rmu[5] ;
      MX c = mxMatMult (a, b, h) ;
      MX d = mxMatMult (b, a, h) ;
      MX e = mxSubstract (c, d, h) ;
      a = kas->Rmu[7] ;
      b = kas->Rmu[6] ;
      c = mxMatMult (a, b, h) ;
      d = mxMatMult (b, a, h) ;
      MX f = mxSubstract (c, d, h) ;
      XT2 = mxAdd (XT2, e, f, h) ;
    }

  if (1) mxNiceShow (XT2) ;
  
  mxValues (XT2, 0, 0, &xx) ;
  for (m1 = 0 ; m1 < d*d ; m1++)
    zz[m1] = -2 * xx[m1] ;

  for (i = 4 ; i < 8 ; i++)
    for (j = 4 ; j < 8 ; j++)
      for (k = 4 ; k < 8 ; k++)
	for (l = 4 ; l < 8 ; l++)
	  {
	    int jj = (i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l);
	    if (jj)
	      {
		MX a = kas->Rmu[i] ;
		MX b = kas->Rmu[j] ;
		MX c = kas->Rmu[k] ;
		MX d1 = kas->Rmu[l] ;
		if (!a || !b || !c || !d1)
		  continue ;
		MX e = mxMatMult (a, b, h) ;
		MX f = mxMatMult (e,c, h) ;
		MX g = mxMatMult (f,d1, h) ;

		mxValues (g, 0, 0, &yy) ;
		for (m1 = 0 ; m1 < d*d ; m1++)
		  zz[m1] += (jj>0  ? yy[m1] : -yy[m1]) ;
	      }
	  }

  mxSet (XT2, zz) ;
  printf ("\nSUCCESS Ghost Casimir operator R16 computed\n") ;
  if (1) mxNiceShow (XT2) ;


  ac_free (h) ;
  return ;
} /* GhostKasimirOperatorR16 */

/***********************************************************************************************************************************************/

static void BBB (void)
{
  complex double z, b[4][4], B[4][4], bb[4][4], BB[4][4] ;
  int i, j, k, l ;
  memset (b, 0, sizeof (b)) ;   /* self dual lower index 2-tensor */
  memset (B, 0, sizeof (b)) ;   /* raise both indices */
  memset (bb, 0, sizeof (bb)) ; /* anti-self-dual lower index 2-tensor */
  memset (BB, 0, sizeof (bb)) ; /* aise both indices */


  for (i = 1 ; i < 4 ; i++)
    {
      randint() ;
      b[0][i] = randint() % 100 ;    
      b[i][0] = - b [0][i] ;
      j = (i - 1 + 1) % 3 + 1 ;    
      k = (i - 1 + 2) % 3 + 1 ;
      b[j][k] = I * b[0][i];
      b[k][j] = - b[j][k] ;
    }
  for (i = 1 ; i < 4 ; i++)
    for (j = 1 ; j < 4 ; j++)
      B[i][j] = b[i][j] ;
  for (i = 1 ; i < 4 ; i++)
    {
      B[0][i] = - b[0][i] ;
      B[i][0] = - b[i][0] ;
    }

  for (i = 1 ; i < 4 ; i++)
    {
      bb[0][i] = randint() % 10 ;    
      bb[i][0] = - bb [0][i] ;
      j = (i - 1 + 1) % 3 + 1 ;    
      k = (i - 1 + 2) % 3 + 1 ;
      bb[j][k] = -I * bb[0][i];
      bb[k][j] = - bb[j][k] ;
    }
  for (i = 1 ; i < 4 ; i++)
    for (j = 1 ; j < 4 ; j++)
      BB[i][j] = bb[i][j] ;
  for (i = 1 ; i < 4 ; i++)
    {
      BB[0][i] = - bb[0][i] ;
      BB[i][0] = - bb[i][0] ;
    }

  /* compute the traces */
  z = 0 ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      z += b[i][j] * B[i][j] ;
  printf ("b:B=%f + I %f\n",creal(z), cimag(z)) ;

  /* compute the traces */
  z = 0 ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      z += b[i][j] * BB[i][j] ;
  printf ("b:BB=%f + I %f\n",creal(z), cimag(z)) ;

  /* compute the traces */
  z = 0 ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      z += bb[i][j] * BB[i][j] ;
  printf ("bb:BB=%f + I %f\n",creal(z), cimag(z)) ;

  /* compute the traces */
  z = 0 ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++) 
	for (l = 0 ; l < 4 ; l++) 
	  z += b[i][j] * B[j][k] * b[k][l] * B[l][i] ;
  printf ("b.B.b.B=%f + I %f\n",creal(z), cimag(z)) ;

  /* compute the traces */
  z = 0 ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++) 
	for (l = 0 ; l < 4 ; l++) 
	  z += 4*b[i][j] * BB[j][k] * b[k][l] * BB[l][i] ;
  printf ("4 b.BB.b.BB=%f + I %f\n",creal(z), cimag(z)) ;

  /* compute the traces */
  z = 0 ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++) 
	for (l = 0 ; l < 4 ; l++) 
	  z += 4*b[i][j] * B[j][k] * bb[k][l] * BB[l][i] ;
  printf ("4 b.B.bb.BB=%f + I %f\n",creal(z), cimag(z)) ;

   /* compute the traces */
  z = 0 ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++) 
	for (l = 0 ; l < 4 ; l++) 
	  z += b[i][j] * B[j][i] * bb[k][l] * BB[l][k] ;
  printf ("b:B*bb:BB=%f + I %f\n",creal(z), cimag(z)) ;
} /* BBB */

/***********************************************************************************************************************************************/

static void KasimirR16 (void)
{
  AC_HANDLE h = ac_new_handle () ;
  KAS kas0, *kas = &kas0 ;
  int ii, i ;
  int d = 16 ;
  int dd = d * d ;
  MX *mu, *Rmu, chi, chiP, chiM, xi, xiP, xiM ;
  MX *Tmu, Tp, Tm ;
  MX M ;
  int xx[dd], xxP[dd], xxM[dd] ;
  float complex zz[dd] ;
  float complex zzP[dd] ;
  float complex zzM[dd] ;
  const float complex *zz1, *zz2 ; 
  BOOL xiPrime = FALSE     ;
  
  memset (kas,0, sizeof (KAS)) ;
  kas->h = h ;
  kas->d = d ;
  kas->xiPrime = xiPrime ;
  
  mu = kas->Rmu = (MX *) halloc (10 * sizeof (MX), kas->h) ;
  for (ii = 0 ; ii < 10 ; ii++)
    mu[ii] = mxCreate (h,  messprintf ("mu[%d]", ii) , MX_COMPLEX, d, d, 0) ;
  Rmu = kas->Rmu = (MX *) halloc (10 * sizeof (MX), kas->h) ;
  for (ii = 0 ; ii < 10 ; ii++)
    Rmu[ii] = mxCreate (h,  messprintf ("Rmu[%d]", ii) , MX_COMPLEX, d, d, 0) ;
  Tmu = kas->Rmu = (MX *) halloc (10 * sizeof (MX), kas->h) ;
  for (ii = 0 ; ii < 10 ; ii++)
    Tmu[ii] = mxCreate (h,  messprintf ("Tmu[%d]", ii) , MX_COMPLEX, d, d, 0) ;
  
  kas->chi16 =  chi = mxCreate (h,  "chi", MX_INT, d, d, 0) ;
  chiP = mxCreate (h,  "chiP", MX_INT, d, d, 0) ;
  chiM = mxCreate (h,  "chiM", MX_INT, d, d, 0) ;
  xi = mxCreate (h,  "xi", MX_COMPLEX, d, d, 0) ;
  xiP = mxCreate (h,  "xiP", MX_COMPLEX, d, d, 0) ;
  xiM = mxCreate (h,  "xiM", MX_COMPLEX, d, d, 0) ;

  /* chi and xi matrices */
  memset (xx, 0, sizeof(xx)) ;
  memset (xxP, 0, sizeof(xx)) ;
  memset (xxM, 0, sizeof(xx)) ;
  memset (zz, 0, sizeof(zz)) ;
  memset (zzP, 0, sizeof(zz)) ;
  memset (zzM, 0, sizeof(zz)) ;

  for (ii = 0 ; ii < 16 ; ii++)
    {
      if (ii % 4 == 0)
	{
	  xx [d*ii + ii] = -1 ;
	  xxM [d*ii + ii] = 1 ;
	  zz [d*ii + ii] = I ;
	  zzM [d*ii + ii] = I ;
	}
      else if (ii % 4 == 3)
	{
	  xx [d*ii + ii] = -1 ;
	  xxM [d*ii + ii] = 1 ;
	  zz [d*ii + ii] = xiPrime ? -I : I ;
	  zzM [d*ii + ii] = xiPrime ? -I : I ;
	}
      else
	{
	  xx [ d*ii + ii] = 1 ;
	  xxP [d*ii + ii] = 1 ;
	  zz [ d*ii + ii] = xiPrime ? -1 : 1 ;
	  zzP [ d*ii + ii] = xiPrime ? -1 : 1 ;
	}
    }
  mxSet (chi, xx) ;
  mxSet (chiP, xxP) ;
  mxSet (chiM, xxM) ;
  mxSet (xi, zz) ;
  mxSet (xiP, zzP) ;
  mxSet (xiM, zzM) ;


  /* Y matrices L0 */
  memset (zz, 0, sizeof(zz)) ;
  Rmu[0] = mxCreate (h,  "Y", MX_COMPLEX, d, d, 0) ;
  for (ii = 0 ; ii < 4 ; ii++)
    {
      if (ii == 1)
	zz[d * ii + ii] = -I ;
      else if (ii == 2)
	zz[d * ii + ii] = -I ;
      else if (ii == 3)
	zz[d * ii + ii] = -2.0I ;
    }
  for (ii = 4 ; ii < 16 ; ii++)
    {
      if (ii % 4 == 0)
	zz[d * ii + ii] = 4.0I/3.0 ;
      else if (ii % 4 == 1)
	zz[d * ii + ii] = I/3.0 ;
      else if (ii % 4 == 2)
	zz[d * ii + ii] = I/3.0 ;
      else if (ii % 4 == 3)
	zz[d * ii + ii] = -2.0I/3.0 ;
    }
  mxSet (Rmu[0], zz) ;

  /* sl(2) matrices L3 */
  memset (zz, 0, sizeof(zz)) ;
  Rmu[3] = mxCreate (h,  "Rmu[3]", MX_COMPLEX, d, d, 0) ;
  for (ii = 0 ; ii < 16 ; ii++)
    {
      if (ii % 4 == 1)
	zz[d * ii + ii] = I ;
      else if (ii % 4 == 2)
	zz[d * ii + ii] = -I ;
    }
  mxSet (Rmu[3], zz) ;

  /* sl(2) matrices L1 */
  memset (zz, 0, sizeof(zz)) ;
  Rmu[1] = mxCreate (h,  "Rmu[1]", MX_COMPLEX, d, d, 0) ;
  for (ii = 0 ; ii < 16 ; ii++)
    {
      if (ii % 4 == 1)
	{
	  zz[d * ii + ii + 1] = I ;
	  zz[d * (ii + 1) + ii] = I ;
	}
    }
  mxSet (Rmu[1], zz) ;

  /* sl(2) matrices L2 */
  memset (zz, 0, sizeof(zz)) ;
  Rmu[2] = mxCreate (h,  "Rmu[2]", MX_COMPLEX, d, d, 0) ;
  for (ii = 0 ; ii < 16 ; ii++)
    {
      if (ii % 4 == 1)
	{
	  zz[d * ii + ii + 1] = 1 ;
	  zz[d * (ii + 1) + ii] = -1 ;
	}
    }
  mxSet (Rmu[2], zz) ;

  /* matrix L8 = (L0 + L3)/2 */
  mxValues (Rmu[0], 0, 0, &zz1) ;
  mxValues (Rmu[3], 0, 0, &zz2) ;

  memset (zz, 0, sizeof(zz)) ;
  Rmu[8] = mxCreate (h,  "Rmu[8]", MX_COMPLEX, d, d, 0) ;
  for (i = 0 ; i < dd ; i++)
    zz[i] = (zz1[i] - zz2[i]) ;
  mxSet (Rmu[8], zz) ;

  memset (zz, 0, sizeof(zz)) ;
  Rmu[9] = mxCreate (h,  "Rmu[9]", MX_COMPLEX, d, d, 0) ;
  for (i = 0 ; i < dd ; i++)
    zz[i] = (zz1[i] + zz2[i]) ;
  mxSet (Rmu[9], zz) ;


  
  mxNiceShow (chi) ;
  mxNiceShow (xi) ;
  mxNiceShow (Rmu[0]) ;
  for (ii = 1 ; ii < 4 ; ii++)
    mxNiceShow (Rmu[ii]) ;
  for (ii = 8 ; ii < 10 ; ii++)
    mxNiceShow (Rmu[ii]) ;
    
  M = KasCommut (Rmu[1], Rmu[2], -1, kas) ;
  M->name = "[1,2]" ;
  mxNiceShow (M) ;
  
  M = KasCommut (Rmu[3], Rmu[1], -1, kas) ;
  M->name = "[3,1]" ;
  mxNiceShow (M) ;
  
  M = KasCommut (Rmu[3], Rmu[2], -1, kas) ;
  M->name = "[3,2]" ;
  mxNiceShow (M) ;

  KasCheckR16 (kas, Rmu[1], Rmu[2], Rmu[3], 2, -1) ;
  KasCheckR16 (kas, Rmu[2], Rmu[3], Rmu[1], 2, -1) ;
  KasCheckR16 (kas, Rmu[3], Rmu[1], Rmu[2], 2, -1) ;
  
  KasCheckR16 (kas, Rmu[0], Rmu[1], Rmu[2], 0, -1) ;
  KasCheckR16 (kas, Rmu[0], Rmu[2], Rmu[2], 0, -1) ;
  KasCheckR16 (kas, Rmu[0], Rmu[3], Rmu[2], 0, -1) ;

  printf ("Success for all even-even commutators\n") ;
  
  /* sl(2/1,R) odd matrix L6 */
  memset (zz, 0, sizeof(zz)) ;
  if (0) Rmu[6] = mxCreate (h,  "Rmu[6]", MX_COMPLEX, d, d, 0) ;
  for (ii = 0 ; ii < 4 ; ii++)
    {
      if (ii % 4 == 0)
	{
	  zz[d * ii + ii + 1] = 0 ;
	  zz[d * (ii + 1) + ii] = 0 ;
	}
      else if (ii % 4 == 2)
	{
	  zz[d * ii + ii + 1] = 1 ;
	  zz[d * (ii + 1) + ii] = 1 ;
	}
    }
  for (ii = 4 ; ii < 16 ; ii++)
    {
      if (ii % 4 == 0)
	{
	  zz[d * ii + ii + 1] = -sqrt(2.0/3.0) ;
	  zz[d * (ii + 1) + ii] = sqrt(2.0/3.0) ; ;
	}
      else if (ii % 4 == 2)
	{
	  zz[d * ii + ii + 1] = sqrt(1.0/3.0) ;
	  zz[d * (ii + 1) + ii] = sqrt(1.0/3.0) ; ;
	}
    }
  mxSet (mu[6], zz) ;
  mxNiceShow (mu[6]) ;
  Rmu[6] = mxMatMult (xi, mu[6], kas->h) ;
  Rmu[6]->name = "Rmu[6]" ;
  mxNiceShow (Rmu[6]) ;
  mxNiceShow (xi) ;
  
  KasCheckR16 (kas, Rmu[6], Rmu[6], Rmu[9], -1, 1) ;

  if (0) exit (0) ;
  
  Rmu[7] = KasCommut (Rmu[3], Rmu[6], -1, kas) ;
  Rmu[7]->name = "Rmu[7]" ;
  Rmu[4] = KasCommut (Rmu[1], Rmu[6], -1, kas) ;
  Rmu[4]->name = "Rmu[4]" ;
  Rmu[5] = KasCommut (Rmu[3], Rmu[4], -1, kas) ;
  Rmu[5]->name = "Rmu[5]" ;
  
  mxNiceShow (Rmu[4]) ;
  mxNiceShow (Rmu[5]) ;
  mxNiceShow (Rmu[6]) ;
  mxNiceShow (Rmu[7]) ;
  
  KasCheckR16 (kas, Rmu[3], Rmu[6], Rmu[7], 1, -1) ;
  KasCheckR16 (kas, Rmu[3], Rmu[7], Rmu[6], -1, -1) ;
  KasCheckR16 (kas, Rmu[3], Rmu[4], Rmu[5], 1, -1) ;
  KasCheckR16 (kas, Rmu[3], Rmu[5], Rmu[4], -1, -1) ;

  KasCheckR16 (kas, Rmu[0], Rmu[6], Rmu[7], -1, -1) ;
  KasCheckR16 (kas, Rmu[0], Rmu[7], Rmu[6], 1, -1) ;
  KasCheckR16 (kas, Rmu[0], Rmu[4], Rmu[5], 1, -1) ;
  KasCheckR16 (kas, Rmu[0], Rmu[5], Rmu[4], -1, -1) ;

  KasCheckR16 (kas, Rmu[1], Rmu[6], Rmu[4], 1, -1) ;
  KasCheckR16 (kas, Rmu[1], Rmu[7], Rmu[5], -1, -1) ;
  KasCheckR16 (kas, Rmu[1], Rmu[4], Rmu[6], -1, -1) ;
  KasCheckR16 (kas, Rmu[1], Rmu[5], Rmu[7], 1, -1) ;

  KasCheckR16 (kas, Rmu[2], Rmu[4], Rmu[7], -1, -1) ;
  KasCheckR16 (kas, Rmu[2], Rmu[5], Rmu[6], -1, -1) ;
  KasCheckR16 (kas, Rmu[2], Rmu[6], Rmu[5], 1, -1) ;
  KasCheckR16 (kas, Rmu[2], Rmu[7], Rmu[4], 1, -1) ;

  printf ("Success for all even-odd commutators\n") ;

    
  KasCheckR16 (kas, Rmu[4], Rmu[4], Rmu[8], -1, 1) ;
  KasCheckR16 (kas, Rmu[5], Rmu[5], Rmu[8], -1, 1) ;
  KasCheckR16 (kas, Rmu[6], Rmu[6], Rmu[9], -1, 1) ;
  KasCheckR16 (kas, Rmu[7], Rmu[7], Rmu[9], -1, 1) ;

  KasCheckR16 (kas, Rmu[4], Rmu[5], Rmu[9], 0, 1) ;
  KasCheckR16 (kas, Rmu[4], Rmu[6], Rmu[2], 1, 1) ;
  KasCheckR16 (kas, Rmu[4], Rmu[7], Rmu[1], -1, 1) ;
  
  KasCheckR16 (kas, Rmu[5], Rmu[4], Rmu[9], 0, 1) ;
  KasCheckR16 (kas, Rmu[5], Rmu[6], Rmu[1], -1, 1) ;
  KasCheckR16 (kas, Rmu[5], Rmu[7], Rmu[2], -1, 1) ;
  
  KasCheckR16 (kas, Rmu[6], Rmu[4], Rmu[2], 1, 1) ;
  KasCheckR16 (kas, Rmu[6], Rmu[5], Rmu[1], -1, 1) ;
  KasCheckR16 (kas, Rmu[6], Rmu[7], Rmu[1], 0, 1) ;
  
  KasCheckR16 (kas, Rmu[7], Rmu[4], Rmu[1], -1, 1) ;
  KasCheckR16 (kas, Rmu[7], Rmu[5], Rmu[2], -1, 1) ;
  KasCheckR16 (kas, Rmu[7], Rmu[6], Rmu[1], 0, 1) ;
  
  printf ("Success for all odd-odd anti-commutators\n") ;

  if (0)
    {  /* i do not know how to twist */
      /* Construct the twisted matrices */
      Tp = mxMatMult (xiP, Rmu[4], kas->h) ;
      Tm = mxMatMult (xiM, Rmu[4], kas->h) ;
      Tmu[4] = mxAdd (Tmu[4], Tp, Tm, kas->h) ;
      
      Tp = mxMatMult (xiP, Rmu[5], kas->h) ;
      Tm = mxMatMult (xiM, Rmu[5], kas->h) ;
      Tmu[5] = mxSubstract (Tp, Tm, kas->h) ;
      
      Tp = mxMatMult (xiP, Rmu[6], kas->h) ;
      Tm = mxMatMult (xiM, Rmu[6], kas->h) ;
      Tmu[6] = mxSubstract (Tp, Tm, kas->h) ;
      
      Tp = mxMatMult (xiP, Rmu[7], kas->h) ;
      Tm = mxMatMult (xiM, Rmu[7], kas->h) ;
      Tmu[7] = mxAdd (Tmu[7], Tp, Tm, kas->h) ;
      
      mxNiceShow (Tmu[6]) ;
      
      KasCheckR16 (kas, Tmu[6], Tmu[6], Rmu[9], -1, 1) ;
      KasCheckR16 (kas, Rmu[3], Tmu[6], Tmu[7], 1, -1) ;
      
      printf ("Success for all Tmu commutators\n") ;
    }

  for (ii = 0 ; ii < 10 ; ii++)
    kas->Rmu[ii] = Rmu[ii] ;
  GhostKasimirOperatorR16 (kas) ;
  
  ac_free (h) ;
} /* KasimirR16 */

/***********************************************************************************************************************************************/
/*****  SU(2/1) representation theory. This is used by, but does not depend on the analysis above of the Feynman diagrams **********************/
/*****  Scalar anomaly paper , indecomposable representations submited to Arxiv and JHEP in My 20, 2020 ****************************************/
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/


#define NTYPES 12

MX *neq[NTYPES] ;
MX *marcu[NTYPES] ;
MX *Marcu[NTYPES] ;
MX nchiT[NTYPES] ;
MX nchiS[NTYPES] ;
MX nchiL[NTYPES] ;
MX nchiR[NTYPES] ;

int ss[] = {4,4,4, 8,8,8, 8,8,8, 8,8,8} ;
MX nn[10], ee[10], qq[10], N2[10], E2[10], Q2[10], N2a[10], E2a[10], Q2a[10], N2b[10], E2b[10], Q2b[10] ;
MX nnmarcu[10], eemarcu[10], eeMarcu[10], qqmarcu[10] ;
MX chiT, chiS, chiL, chiR ;
MX chiT2, chiS2, chiL2, chiR2 ;
MX SG[4], SB[4] ;
float gg[4][4] = {{-1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}} ;
complex float eps[4][4][4][4] ;
complex float PP[4][4][4][4] ;
complex float PM[4][4][4][4] ;

static int c3Mask = 0 ;
static int SU3 = 0 ;
static int myType = -1 ;
static MX muHermite (MX m, AC_HANDLE h)
{
  MX m1 = 0 ;
  m1 = mxMatTranspose (m1, m, h) ;
  mxConjugate (m1, m1, h) ;
  
  return m1 ;
} /* muHermite */

/*************************************************************************************************/

static void muStructure (void)
{
  AC_HANDLE h = ac_new_handle () ;
  int a,b,c, t ;
  MX mm1, mm2, mm3, mm4 ;
  FABC *f ;
  FABC ff[] = {
    {1,2,3,-1,2.0I, "[m_1,m_2] = 2i m_3"},
    {2,3,1,-1,2.0I, "[m_2,m_3] = 2i m_1"},
    {3,1,2,-1,2.0I, "[m_3,m_1] = 2i m_2"},

    {1,4,7,-1, I,"\n  [m_1,m_4] =  i m_7"},
    {1,5,6,-1,-I,    "[m_1,m_5] = -i m_6"},
    {1,6,5,-1, I,    "[m_1,m_6] =  i m_5"},
    {1,7,4,-1,-I,    "[m_1,m_7] = -i m_4"},

    {2,4,6,-1, I,"\n  [m_2,m_4] =  i m_6"},
    {2,5,7,-1, I,    "[m_2,m_5] =  i m_7"},
    {2,6,4,-1,-I,    "[m_2,m_6] = -i m_4"},
    {2,7,5,-1,-I,    "[m_2,m_7] = -i m_5"},

    {3,4,5,-1, I,"\n  [m_3,m_4] =  i m_5"},
    {3,5,4,-1,-I,    "[m_3,m_5] = -i m_4"},
    {3,6,7,-1,-I,    "[m_3,m_6] = -i m_7"},
    {3,7,6,-1, I,    "[m_3,m_7] =  i m_6"},

    {0,4,5,-1, I,"\n  [m_0,m_4] =  i m_5"},
    {0,5,4,-1,-I,    "[m_0,m_5] = -i m_4"},
    {0,6,7,-1, I,    "[m_0,m_6] =  i m_7"},
    {0,7,6,-1,-I,    "[m_0,m_7] = -i m_6"},

    {4,4,9, 1,1,"\n  {m_4,m_4} = -m_9"},
    {5,5,9, 1,1,    "{m_5,m_5} = -m_9"},

    {6,6,8, 1,1,   "{m_6,m_6} = -m_8"},
    {6,7,0, 1, 0,    "{m_6,m_7} = 0  "},
    {7,7,8, 1,1,   "{m_7,m_7} = -m_8"},
    {4,5,0, 1, 0,"\n  {m_4,m_5} = 0  "},
    {6,7,0, 1, 0,    "{m_6,m_7} = 0  "},

    {4,6,1, 1, -1,"\n  {m_4,m_6} = m_1"},
    {5,6,2, 1, -1,    "{m_5,m_6} = m_2"},
    {4,7,2, 1,1,"\n  {m_4,m_7} = -m_2"},
    {5,7,1, 1, -1,    "{m_5,m_7} = m_1"},


    /*

    */
    {-1,0,0,0,0}
  } ;

  for (f = ff ; f->a >=0 ; f++)
    {
      a = f->a ;
      b = f->b ;
      c = f->c ;

      printf ("# %s\t", f->title) ;
      for (t = 0 ; t < NTYPES ; t++)
	{
	  double z ;
	  if (0 && t != 11) continue ;
	  mm1  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
	  mm2  = mxCreate (h,  "mm2", MX_COMPLEX, ss[t], ss[t], 0) ;
	  mm3  = mxCreate (h,  "mm3", MX_COMPLEX, ss[t], ss[t], 0) ;
	  mm4  = mxCreate (h,  "mm4", MX_COMPLEX, ss[t], ss[t], 0) ;

	  mm1 = mxMatMult (neq[t][a], neq[t][b], h) ;
	  mm2 = mxMatMult (neq[t][b], neq[t][a], h) ;
	  mm3 = mxLinearCombine (mm3, 1, mm1, f->sign, mm2, h) ;
	  mm4 = mxLinearCombine (mm4, 1, mm3, f->z,  neq[t][c], h) ;
	  
	  if (0 && a == 2 && b == 7 && t == 7)
	    {
	      printf ("\n# mu(%d) type %d\n", a, t) ;
	      mxNiceShow  (neq[t][a]) ;
	      printf ("\n# mu(%d) type %d\n", b, t) ;
	      mxNiceShow  (neq[t][b]) ;
	      printf ("\n# [a,b] type %d\n", t) ;
	      mxNiceShow (mm3) ;
	      printf ("\n# [a,b] should be equal to neq[%d][%d]\n", t,c) ;
	      mxNiceShow (neq[t][c]) ;
	      printf ("\n# norm");
	    }
	  z = mxFNorm (mm4) ;
	  if (z < .0000001) z = 0 ;
	  printf ("\t%.2g", z) ;
	  }
	printf ("\n") ;
      }
 
  if (1)
    {
      mxNiceShow (neq[6][6]) ;
      mxNiceShow (neq[6][7]) ;
    }
  ac_free (h) ;
  return ;
} /* muStructure */

/*************************************************************************************************/

static void muSigma (AC_HANDLE h)
{
  int i, j, k, l, m, n ;
  float z ;

  memset (eps, 0, sizeof(eps)) ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	for (l = 0 ; l < 4 ; l++)
	  { /* checked, this is correct */
	    n =(i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l) ;
	    if (n > 0)
	      eps[i][j][k][l] = n = 1 ;
	    else if (n < 0)
	      eps[i][j][k][l] = n = -1 ;
	    if (0) if (n) printf ("espison(%d,%d,%d,%d) = %d\n", i,j,k,l,n) ;
	  }


  complex float sg0[] = {1,0,0,1} ;
  complex float sb0[] = {-1,0,0,-1} ;
  complex float sg1[] = {0,1,1,0} ;
  complex float sg2[] = {0,I,-I,0} ;
  complex float sg3[] = {1,0,0,-1} ;

  SG[0] = mxCreate (h, "Sigma_0", MX_COMPLEX, 2, 2, 0) ;
  SG[1] = mxCreate (h, "Sigma_1", MX_COMPLEX, 2, 2, 0) ;
  SG[2] = mxCreate (h, "Sigma_2", MX_COMPLEX, 2, 2, 0) ;
  SG[3] = mxCreate (h, "Sigma_3", MX_COMPLEX, 2, 2, 0) ;

  SB[0] = mxCreate (h, "SB_0", MX_COMPLEX, 2, 2, 0) ;
  SB[1] = mxCreate (h, "SB_1", MX_COMPLEX, 2, 2, 0) ;
  SB[2] = mxCreate (h, "SB_2", MX_COMPLEX, 2, 2, 0) ;
  SB[3] = mxCreate (h, "SB_3", MX_COMPLEX, 2, 2, 0) ;

  mxSet (SG[0], sg0) ;
  mxSet (SB[0], sb0) ;
  mxSet (SG[1], sg1) ;
  mxSet (SB[1], sg1) ;
  mxSet (SG[2], sg2) ;
  mxSet (SB[2], sg2) ;
  mxSet (SG[3], sg3) ;
  mxSet (SB[3], sg3) ;

  printf ("### Verify that the sigma sigma-bar obey the Clifford algebra sg_i sb_j + sg_j sb_i = 2 g_ij Identity[2] \n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      {
	AC_HANDLE h = ac_new_handle () ;
	MX mmm[3] ; 
	MX mm1 =  mxCreate (h, "ij", MX_COMPLEX, 2, 2, 0) ;
	MX mm2 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	MX mm3 =  mxCreate (h, "ij+ji", MX_COMPLEX, 2, 2, 0) ;
	MX mm4 =  mxCreate (h, "zero", MX_COMPLEX, 2, 2, 0) ;

	mm1 = mxMatMult (SG[i], SB[j], h) ;
	mm2 = mxMatMult (SG[j], SB[i], h) ;

	mmm[0] = SG[i] ;
	mmm[1] = SB[j] ;
	mmm[2] = 0 ;
	mm1 = mxMatListMult (h, mmm) ;
	if (0)
	  {
	    printf ("###### sigma sbar : i=%d j=%d \n", i, j) ;
	    mxNiceShow (mm1) ;
	  }

	mmm[0] = SG[j] ;
	mmm[1] = SB[i] ;
	mmm[2] = 0 ;
	mm2 = mxMatListMult (h, mmm) ;
	if (0) mxNiceShow (mm2) ;
	mm3 = mxLinearCombine (mm3, 1,mm1, 1, mm2, h) ; 
	mm4 = mxLinearCombine (mm4, 1,mm3, -2*gg[i][j], SG[0], h) ; 
	if (0)
	  {
	    mxNiceShow (mm1) ;
	    mxNiceShow (mm2) ;	
	    mxNiceShow (mm3) ;	
	    mxNiceShow (mm4) ;
	  }
	z = mxFNorm(mm4) ;
	if (z > .001)
	  printf ("###### sigma sbar : i=%d j=%d {i,j} = 2 g_ij Id :: verif %g\n",i,j, z) ;
	
	ac_free (h) ;	
      }
  if (0) exit (0) ;

  printf ("### Verify that the sigma sigma-bar obey the Clifford algebra sby_i sg_j + sb_j sg_i = 2 g_ij Identity[2] \n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      {
	AC_HANDLE h = ac_new_handle () ;

	MX mm1 =  mxCreate (h, "ij", MX_COMPLEX, 2, 2, 0) ;
	MX mm2 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	MX mm3 =  mxCreate (h, "ij+ji", MX_COMPLEX, 2, 2, 0) ;
	MX mm4 =  mxCreate (h, "zero", MX_COMPLEX, 2, 2, 0) ;

	mm1 = mxMatMult (SB[i], SG[j], h) ;
	mm2 = mxMatMult (SB[j], SG[i], h) ;
	mm3 = mxLinearCombine (mm3, 1,mm1, 1, mm2, h) ; 

	mm4 = mxLinearCombine (mm4, 1,mm3, -2*gg[i][j], SG[0], h) ; 
	if (0)
	  {
	    mxNiceShow (mm1) ;
	    mxNiceShow (mm2) ;	
	    mxNiceShow (mm3) ;	
	    mxNiceShow (mm4) ;
	  }
	z = mxFNorm(mm4) ;
	if (z > .001)
	  printf ("###### sbar sigma : i=%d j=%d {i,j} = 2 g_ij Id :: verif %g\n",i,j, z) ; 
	
	ac_free (h) ;	
      }

  /* check the projectors */
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	for (l = 0 ; l < 4 ; l++)
	  {
	    PP[i][j][k][l] = (gg[i][k]*gg[j][l] - gg[i][l]*gg[j][k] + I * eps[i][j][k][l])/4.0 ;
	    PM[i][j][k][l] = (gg[i][k]*gg[j][l] - gg[i][l]*gg[j][k] - I * eps[i][j][k][l])/4.0 ;
	  }

  printf("### Verify that PP is a projector PP^2 = PP\n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	for (l = 0 ; l < 4 ; l++)
	  {
	    float z ;
	    complex float z2 = 0, z1 = PP[i][j][k][l] ;
	    int a, b ;
	    for (a = 0 ; a < 4 ; a++)
	      for (b = 0 ; b < 4 ; b++)
		z2 += PP[i][j][a][b] *gg[a][a] * gg[b][b] * PP[a][b][k][l] ;
	    z = cabsf (z2 - z1) ;
	    if (z > minAbs)
	      printf("PP PP - PP not zero ijkl = %d %d %d %d  zz=%g\n", i,j,k,l,z) ;
	  }

  printf("### Verify that PM is a projector PM^2 = PM\n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	for (l = 0 ; l < 4 ; l++)
	  {
	    float z ;
	    complex float z2 = 0, z1 = PM[i][j][k][l] ;
	    int a, b ;
	    for (a = 0 ; a < 4 ; a++)
	      for (b = 0 ; b < 4 ; b++)
		z2 += PM[i][j][a][b] *gg[a][a] * gg[b][b] * PM[a][b][k][l] ;
	    z = cabsf (z2 - z1) ;
	    if (z > minAbs)
	      printf("PM PM - PM not zero ijkl = %d %d %d %d  zz=%g\n", i,j,k,l,z) ;
	  }

  printf("### Verify that PP is a projector PP PM = 0\n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	for (l = 0 ; l < 4 ; l++)
	  {
	    float z, z2 = 0, z1 = 0 ;
	    int a, b ;
	    for (a = 0 ; a < 4 ; a++)
	      for (b = 0 ; b < 4 ; b++)
		z2 += PP[i][j][a][b] *gg[a][a] * gg[b][b] * PM[a][b][k][l] ;
	    z = fabsf (z2 - z1) ;
	    if (z > minAbs)
	      printf("PP PM  not zero ijkl = %d %d %d %d  zz=%g\n", i,j,k,l,z) ;
	  }

  printf("### Verify that SG SB = PP SG SB\n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      {
	AC_HANDLE h = ac_new_handle () ;

	MX mm1 =  mxCreate (h, "ij", MX_COMPLEX, 2, 2, 0) ;
	MX mm2 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	MX mm3 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	int a, b ;

	mm2 = mxMatMult (SG[i], SB[j], h) ;
	mm3 = mxMatMult (SG[j], SB[i], h) ;
	mm1 = mxLinearCombine (mm1, 0.5,mm2, -0.5, mm3, h) ;
	if (0)
	  {
	    printf ("## i=%d j=%d ::\n", i, j) ;
	    mxNiceShow (mm1) ;
	  }
	mm2 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;

	for (a = 0 ; a < 4 ; a++)
	  for (b = 0 ; b < 4 ; b++)
	    {
	      MX mm ;
	      mm =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	      mm = mxMatMult (SG[a], SB[b], h) ;
	      mm3 =  mxCreate (h, "mm2", MX_COMPLEX, 2, 2, 0) ;
	      mm3 = mxLinearCombine (mm3, 1,mm2, PP[i][j][a][b], mm, h) ; 
	      mm2 = mm3 ;
	    }
	mm3 =  mxCreate (h, "zero", MX_COMPLEX, 2, 2, 0) ;
	if (0) mxNiceShow (mm2) ;
	mm3 = mxLinearCombine (mm3, 1,mm1, -1*gg[i][i]*gg[j][j], mm2, h) ;
	z = mxFNorm(mm3) ;
	if (z > minAbs)
	  printf ("###### S_i sb_j not equal PP s sb: sbar i=%d j=%d z = %f\n", i,j,z) ;
	
	ac_free (h) ;	
      }


  printf("### Verify that SB SG = PM SB SG\n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      {
	AC_HANDLE h = ac_new_handle () ;

	MX mm1 =  mxCreate (h, "ij", MX_COMPLEX, 2, 2, 0) ;
	MX mm2 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	MX mm3 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	int a, b ;

	mm2 = mxMatMult (SB[i], SG[j], h) ;
	mm3 = mxMatMult (SB[j], SG[i], h) ;
	mm1 = mxLinearCombine (mm1, 0.5,mm2, -0.5, mm3, h) ;
	if (0)
	  {
	    printf ("## i=%d j=%d ::\n", i, j) ;
	    mxNiceShow (mm1) ;
	  }
	mm2 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;

	for (a = 0 ; a < 4 ; a++)
	  for (b = 0 ; b < 4 ; b++)
	    {
	      MX mm ;
	      mm =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	      mm = mxMatMult (SB[a], SG[b], h) ;
	      mm3 =  mxCreate (h, "mm2", MX_COMPLEX, 2, 2, 0) ;
	      mm3 = mxLinearCombine (mm3, 1,mm2, PM[i][j][a][b], mm, h) ; 
	      mm2 = mm3 ;
	    }
	mm3 =  mxCreate (h, "zero", MX_COMPLEX, 2, 2, 0) ;
	if (0) mxNiceShow (mm2) ;
	mm3 = mxLinearCombine (mm3, 1,mm1, -1*gg[i][i]*gg[j][j], mm2, h) ;
	z = mxFNorm(mm3) ;
	if (z > minAbs)
	  printf ("###### S_i sb_j not equal PP s sb: sbar i=%d j=%d z = %f\n", i,j,z) ;
	
	ac_free (h) ;	
      }

  printf("### Verify that Tr(SG_i SB_j SG_k SB_l = 2 * (g_ijg_kl - g_ik_g_jl+gil_gjk + I eps_ijkl)\n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	for (l = 0 ; l < 4 ; l++)
	  {
	    float z ;
	    float complex z1, z2 ;
	    AC_HANDLE h = ac_new_handle () ;
	    MX mmm[5] ;

	    MX mm2 =  mxCreate (h, "ji", MX_COMPLEX, 2, 3, 0) ;
	    MX mm3 =  mxCreate (h, "ji", MX_COMPLEX, 3, 5, 0) ;
	    MX mm4 =  0 ;

	    mm4 = mxMatMult (mm2, mm3, h) ;



	    mmm[0] = mm2 ;
	    mmm[1] = mm3 ;
	    mmm[2] = 0 ;
	    mmm[4] = 0 ;


	    mm4 = mxMatListMult (h, mmm) ;               


	    mmm[0] = SG[i] ;
	    mmm[1] = SB[j] ;
	    mmm[2] = SG[k] ;
	    mmm[3] = SB[l] ;
	    mmm[4] = 0 ;


	    mm4 = mxMatListMult (h, mmm) ;               
	    z1 = mxMatTrace (mm4) ;
	    z2 = 2*(gg[i][j]*gg[k][l] - gg[i][k]*gg[j][l] + gg[i][l]*gg[j][k] + I*eps[i][j][k][l]) ;
	    z = cabsf (z2 - z1) ;
	    if (z > minAbs)
	      printf ("###### Trace (sigma ijkl) not equal gg - gg + gg + i epsilon: i=%d j=%d k=%d l=%d z = %f\n", i,j,k,l,z) ;
	    
	    ac_free (h) ;	
	  }

  printf("### Verify that Tr(SB_i SG_j SB_k SG_l = 2 * (g_ijg_kl - g_ik_g_jl+gil_gjk - I eps_ijkl)\n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	for (l = 0 ; l < 4 ; l++)
	  {
	    float z ;
	    float complex z1, z2 ;
	    AC_HANDLE h = ac_new_handle () ;
	    MX mmm[5] ;
	    MX mm4 =  0 ;
	    
	    mmm[0] = SB[i] ;
	    mmm[1] = SG[j] ;
	    mmm[2] = SB[k] ;
	    mmm[3] = SG[l] ;
	    mmm[4] = 0 ;
		
	    mm4 = mxMatListMult (h, mmm) ;               
	    z1 = mxMatTrace (mm4) ;
	    z2 = 2*(gg[i][j]*gg[k][l] - gg[i][k]*gg[j][l] + gg[i][l]*gg[j][k] - I*eps[i][j][k][l]) ;
	    z = cabsf (z2 - z1) ;
	    if (z > minAbs)
	      printf ("###### Trace (sb sg ijkl) not equal 2 *( gg - gg + gg - i epsilon): i=%d j=%d k=%d l=%d z = %f\n", i,j,k,l,z) ;
	    
	    ac_free (h) ;	
	  }


  printf("### Verify that Tr(SG_i^6 = 2 * (ggg 15 terms + i g epsilon 15 terms)\n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	for (l = 0 ; l < 4 ; l++)
	  for (m = 0 ; m < 4 ; m++)
	    for (n = 0 ; n < 4 ; n++)
	      {
		float z ;
		float complex z1, z2 ;
		AC_HANDLE h = ac_new_handle () ;
		MX mmm[7] ;
		MX mm4 =  0 ;	
		
		mmm[0] = SG[i] ;
		mmm[1] = SB[j] ;
		mmm[2] = SG[k] ;
		mmm[3] = SB[l] ;
		mmm[4] = SG[m] ;
		mmm[5] = SB[n] ;
		mmm[6] = 0 ;
		
		
		mm4 = mxMatListMult (h, mmm) ;               
		z1 = mxMatTrace (mm4) ;						     
		z2 = 2*(
			+ gg[i][j]*gg[k][l]*gg[m][n] - gg[i][j]*gg[k][m]*gg[l][n] + gg[i][j]*gg[k][n]*gg[l][m]
			- gg[i][k]*gg[j][l]*gg[m][n] + gg[i][k]*gg[j][m]*gg[l][n] - gg[i][k]*gg[j][n]*gg[l][m]
			+ gg[i][l]*gg[j][k]*gg[m][n] - gg[i][l]*gg[j][m]*gg[k][n] + gg[i][l]*gg[j][n]*gg[k][m]
			- gg[i][m]*gg[j][k]*gg[l][n] + gg[i][m]*gg[j][l]*gg[k][n] - gg[i][m]*gg[j][n]*gg[k][l]
			+ gg[i][n]*gg[j][k]*gg[l][m] - gg[i][n]*gg[j][l]*gg[k][m] + gg[i][n]*gg[j][m]*gg[k][l]
			) ;
		z2 += 2*(
			 + gg[i][j]*I*eps[k][l][m][n]
			 - gg[i][k]*I*eps[j][l][m][n]
			 + gg[i][l]*I*eps[j][k][m][n]
			 - gg[i][m]*I*eps[j][k][l][n]
			 + gg[i][n]*I*eps[j][k][l][m]

			 + gg[j][k]*I*eps[i][l][m][n]
			 - gg[j][l]*I*eps[i][k][m][n]
			 + gg[j][m]*I*eps[i][k][l][n]
			 - gg[j][n]*I*eps[i][k][l][m]

			 + gg[k][l]*I*eps[i][j][m][n]
			 - gg[k][m]*I*eps[i][j][l][n]
			 + gg[k][n]*I*eps[i][j][l][m]

			 + gg[l][m]*I*eps[i][j][k][n]
			 - gg[l][n]*I*eps[i][j][k][m]

			 + gg[m][n]*I*eps[i][j][k][l]			 
			 ) ;
		z = cabsf (z2 - z1) ;
		if (z > minAbs)
		  {
		    printf ("###### Trace (sigma ijkl) not equal gg - gg + gg + i epsilon: i=%d j=%d k=%d l=%d m=%d n=%d z = %f\n", i,j,k,l,m,n,z) ;
		    exit (1) ;
		  }
		ac_free (h) ;	
	      }

  printf("### Verify that Tr(SB_i^6 = 2 * (ggg 15 terms - i g epsilon 15 terms)\n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	for (l = 0 ; l < 4 ; l++)
	  for (m = 0 ; m < 4 ; m++)
	    for (n = 0 ; n < 4 ; n++)
	      {
		float z ;
		float complex z1, z2 ;
		AC_HANDLE h = ac_new_handle () ;
		MX mmm[7] ;
		MX mm4 =  0 ;	
		
		mmm[0] = SB[i] ;
		mmm[1] = SG[j] ;
		mmm[2] = SB[k] ;
		mmm[3] = SG[l] ;
		mmm[4] = SB[m] ;
		mmm[5] = SG[n] ;
		mmm[6] = 0 ;
		
		
		mm4 = mxMatListMult (h, mmm) ;               
		z1 = mxMatTrace (mm4) ;						     
		z2 = 2*(
			+ gg[i][j]*gg[k][l]*gg[m][n] - gg[i][j]*gg[k][m]*gg[l][n] + gg[i][j]*gg[k][n]*gg[l][m]
			- gg[i][k]*gg[j][l]*gg[m][n] + gg[i][k]*gg[j][m]*gg[l][n] - gg[i][k]*gg[j][n]*gg[l][m]
			+ gg[i][l]*gg[j][k]*gg[m][n] - gg[i][l]*gg[j][m]*gg[k][n] + gg[i][l]*gg[j][n]*gg[k][m]
			- gg[i][m]*gg[j][k]*gg[l][n] + gg[i][m]*gg[j][l]*gg[k][n] - gg[i][m]*gg[j][n]*gg[k][l]
			+ gg[i][n]*gg[j][k]*gg[l][m] - gg[i][n]*gg[j][l]*gg[k][m] + gg[i][n]*gg[j][m]*gg[k][l]

			) ;
		z2 += -2*(
			 + gg[i][j]*I*eps[k][l][m][n]
			 - gg[i][k]*I*eps[j][l][m][n]
			 + gg[i][l]*I*eps[j][k][m][n]
			 - gg[i][m]*I*eps[j][k][l][n]
			 + gg[i][n]*I*eps[j][k][l][m]

			 + gg[j][k]*I*eps[i][l][m][n]
			 - gg[j][l]*I*eps[i][k][m][n]
			 + gg[j][m]*I*eps[i][k][l][n]
			 - gg[j][n]*I*eps[i][k][l][m]

			 + gg[k][l]*I*eps[i][j][m][n]
			 - gg[k][m]*I*eps[i][j][l][n]
			 + gg[k][n]*I*eps[i][j][l][m]

			 + gg[l][m]*I*eps[i][j][k][n]
			 - gg[l][n]*I*eps[i][j][k][m]

			 + gg[m][n]*I*eps[i][j][k][l]			 
			 ) ;
		z = cabsf (z2 - z1) ;
		if (z > minAbs)
		  {
		    printf ("###### Trace (sigma ijkl) not equal gg - gg + gg + i epsilon: i=%d j=%d k=%d l=%d m=%d n=%d z = %f\n", i,j,k,l,m,n,z) ;
		    exit (1) ;
		  }
		ac_free (h) ;	
	      }


  printf("### Verify that Tr(SG_i SB_j SG_k SB_m SG_k SB_n = 2 * (ggg 15 terms + i g epsilon 15 terms)\n") ;
  int N = 2 ;
  for (i  = 0 ; i < N ; i++)
    for (j = 0 ; j < N ; j++)
      for (k = 0 ; k < N ; k++)
	for (l = 0 ; l < N ; l++)
	  for (m = 0 ; m < N ; m++)
	    for (n = 0 ; n < N ; n++)
	      for (o = 0 ; o < N ; o++)
		for (p = 0 ; p < N ; p++)
		  {
		    int x[9] ;
		    AC_HANDLE h = ac_new_handle () ;
		    x[0] = i ;
		    x[1] = i ;
		    x[2] = i ;
		    x[3] = i ;
		    x[4] = i ;
		    x[5] = i ;
		    x[6] = i ;
		    x[7] = i ;

		    int z = 0 ;
		    for (int i1 = 0 ; i1 < 8 ; i1++)
		      {
			z = 1 - 2 * (i%2) ;
			for (int i2 = i1+1 ; i2 < 8 ; i2++)
			  {
			    z *= gg[x[i1]][x[i2]] ;
			    for (int i3 = i1 + 1 ; z > 0 && i3 < 8 ; i3++)
			      {
				if (i3 == i2) continue ;
				z *=-1 ;
				for (int i4 = i3 + 1 ; i4 < 8 ; i4++)
				  {
				    if (i4 == i2) continue ;
				    z *=-1 ;
				    z *= gg[x[i3]][x[i4]] ;
				    z *=-1 ;
				    for (int i5 = i3 + 1 ; z > 0 && i5 < 8 ; i5++)
				      {
					if (i5 == i2 || i5 == i4) continue ;
					for (int i6 = i5 + 1 ; i6 < 8 ; i6++)
					  {
					    if (i6 == i2 || i6 == i4) continue ;
					    z *= gg[x[i5]][x[i6]] ;
					    for (int i7 = i5 + 1 ; z > 0 && i7 < 8 ; i7++)
					      {
						if (i7 == i6 || i7 == i4 || i7 == i2) continue ;
						for (int i8 = i7 + 1 ; i8 < 8 ; i8++)
						  {
						    if (i8 == i2 || i8 == i4 || i8 == i6) continue ;
						    z *= gg[x[i7]][x[i8]] ;
						  }
					      }
					  }
				      }
				  }
			      }
			  }
		      }
		    MX mmm[9] ;
		    MX mm4 =  0 ;	
		    
		    mmm[0] = SB[i] ;
		    mmm[1] = SG[j] ;
		    mmm[2] = SB[k] ;
		    mmm[3] = SG[l] ;
		    mmm[4] = SB[m] ;
		    mmm[5] = SG[n] ;
		    mmm[6] = SB[o] ;
		    mmm[7] = SG[p] ;
		    mmm[8] = 0 ;
		    
		    
		    mm4 = mxMatListMult (h, mmm) ;               
		    complex float z1 = mxMatTrace (mm4) ;
		    z1 -= 2*z ;
		    if (cabsf(z1) > .1)
		      messcrash ("error") ;
		    ac_free (h) ;
		  }


  return ;
} /* muSigma */

/*************************************************************************************************/

static void muInit (AC_HANDLE h)
{
  int t, i ;
  float s2 = sqrt (2) ;
  float s3 = sqrt (3) ;
  MX mm =  mxCreate (h, "mm", MX_COMPLEX, 4, 4, 0) ;
  MX mm2 = 0 ;

  complex float mu1[] = {0,0,0,0, 0,0,1,0, 0,1,0,0, 0,0,0,0} ;
  complex float mu2[] = {0,0,0,0, 0,0,-I,0, 0,I,0,0, 0,0,0,0} ;
  complex float mu3[] = {0,0,0,0, 0,1,0,0, 0,0,-1,0, 0,0,0,0} ;

  complex float mu0n[] = {1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,-1} ;
  complex float mu0e[] = {0,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0, -2} ;
  complex float mu0SU3[] = {0,0,0,0, 0,1/s3,0,0, 0,0,1/s3,0, 0,0,0, -2/s3} ;
  complex float mu0q[] = {4/3.0,0,0,0, 0,1/3.0,0,0, 0,0,1/3.0,0, 0,0,0,-2/3.0} ;

  complex float mu8n[] = {1,0,0,0, 0,1,0,0, 0,0,-1,0, 0,0,0,-1} ;
  complex float mu8e[] = {0,0,0,0, 0,0,0,0, 0,0,-2,0, 0,0,0,-2} ;
  complex float mu8q[] = {4/3.0,0,0,0, 0,4/3.0,0,0, 0,0,-2/3.0,0, 0,0,0,-2/3.0} ;

  complex float mu9n[] = {1,0,0,0, 0,-1,0,0, 0,0,1,0, 0,0,0,-1} ;
  complex float mu9e[] = {0,0,0,0, 0,-2,0,0, 0,0,0,0, 0,0,0,-2} ;
  complex float mu9q[] = {4/3.0,0,0,0, 0,-2/3.0,0,0, 0,0,4/3.0,0, 0,0,0,-2/3.0} ;

  complex float mu4n[] = {0,0,-1/s2,0, 0,0,0,1/s2, 1/s2,0,0,0, 0,1/s2,0,0} ;
  complex float mu5n[] = {0,0,I/s2,0, 0,0,0,-I/s2, I/s2,0,0,0, 0,I/s2,0,0} ;
  complex float mu6n[] = {0,1/s2,0,0, -1/s2,0,0,0, 0,0,0,1/s2, 0,0,1/s2,0} ;
  complex float mu7n[] = {0,-I/s2,0,0, -I/s2,0,0,0, 0,0,0,-I/s2, 0,0,I/s2,0} ;

  complex float mu4e[] = {0,0,0,0, 0,0,0,1, 0,0,0,0, 0,1,0,0} ;
  complex float mu5e[] = {0,0,0,0, 0,0,0,-I, 0,0,0,0, 0,I,0,0} ;
  complex float mu6e[] = {0,0,0,0, 0,0,0,0, 0,0,0,1, 0,0,1,0} ;
  complex float mu7e[] = {0,0,0,0, 0,0,0,0, 0,0,0,-I, 0,0,I,0} ;

  complex float mu4q[] = {0,0,-s2/s3,0, 0,0,0,1/s3, s2/s3,0,0,0, 0,1/s3,0,0} ;
  complex float mu5q[] = {0,0,I*s2/s3,0, 0,0,0,-I/s3, I*s2/s3,0,0,0, 0,I/s3,0,0} ;
  complex float mu6q[] = {0,s2/s3,0,0, -s2/s3,0,0,0, 0,0,0,1/s3, 0,0,1/s3,0} ;
  complex float mu7q[] = {0,-I*s2/s3,0,0, -I*s2/s3,0,0,0, 0,0,0,-I/s3, 0,0,I/s3,0} ;

  complex float xT[] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1} ;
  complex float xS[] = {-1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,-1} ;
  complex float xL[] = { 0,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,0} ;
  complex float xR[] = { 1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1} ;


  complex float marcu0n[] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1} ;
  complex float marcu4n[] = {0,0,1/s2,0, 0,0,0,1/s2, -1/s2,0,0,0, 0,1/s2,0,0} ;
  complex float marcu5n[] = {0,0,I/s2,0, 0,0,0,I/s2, I/s2,0,0,0, 0,-I/s2,0,0} ;
  complex float marcu6n[] = {0,-1/s2,0,0, 1/s2,0,0,0, 0,0,0,1/s2, 0,0,1/s2,0} ;
  complex float marcu7n[] = {0,-I/s2,0,0, -I/s2,0,0,0, 0,0,0,I/s2, 0,0,-I/s2,0} ;

  
  complex float marcu0e[] = {2*s2/3,0,0,0, 0,2*s2/3,0,0, 0,0,2*s2/3,0, 0,0,0,2*s2/3} ;
  complex float marcu4e[] = {0,0,0,0, 0,0,0,1, -1,0,0,0, 0,1,0,0} ;
  complex float marcu5e[] = {0,0,0,0, 0,0,0,-I, -I,0,0,0, 0,I,0,0} ;
  complex float marcu6e[] = {0,0,0,0, 1,0,0,0, 0,0,0,1, 0,0,1,0} ;
  complex float marcu7e[] = {0,0,0,0, I,0,0,0, 0,0,0,-I, 0,0,I,0} ;

  complex float Marcu4e[] = {0,0,-2,0, 0,0,0,1, 0,0,0,0, 0,1,0,0} ;
  complex float Marcu5e[] = {0,0,2*I,0, 0,0,0,-I, 0,0,0,0, 0,I,0,0} ;
  complex float Marcu6e[] = {0,2,0,0, 0,0,0,0, 0,0,0,1, 0,0,1,0} ;
  complex float Marcu7e[] = {0,-2*I,0,0,0,0,0,0, 0,0,0,-I, 0,0,I,0} ;

  /* ERROR in eq appendix H.2 of Scalar paper:
   * in the paper we should replace sqrt(2) by -sqrt(2) in equation H.2
   * there is probably a related error of sign in H.3
   * No conclusion is modified
   */
  
  complex float marcu0q[] = {2*s2/3,0,0,0, 0,2*s2/3,0,0, 0,0,2*s2/3,0, 0,0,0,2*s2/3} ;
  complex float marcu4q[] = {0,0,1/s3,0, 0,0,0,s2/s3, -1/s3,0,0,0, 0,s2/s3,0,0} ;
  complex float marcu5q[] = {0,0,I/s3,0, 0,0,0,I*s2/s3, I/s3,0,0,0, 0,-I*s2/s3,0,0} ;
  complex float marcu6q[] = {0,-1/s3,0,0, 1/s3,0,0,0, 0,0,0,s2/s3, 0,0,s2/s3,0} ;
  complex float marcu7q[] = {0,-I/s3,0,0, -I/s3,0,0,0, 0,0,0,I*s2/s3, 0,0,-I*s2/s3,0} ;

  
  chiT = mxCreate (h, "chiT", MX_COMPLEX, 4, 4, 0) ;
  chiS = mxCreate (h,  "chi", MX_COMPLEX, 4, 4, 0) ;
  chiL = mxCreate (h, "chiL", MX_COMPLEX, 4, 4, 0) ;
  chiR = mxCreate (h, "chiR", MX_COMPLEX, 4, 4, 0) ;

  mxSet (chiT, xT) ;
  mxSet (chiS, xS) ;
  mxSet (chiL, xL) ;
  mxSet (chiR, xR) ;

  for (t = 0 ; t < 3 ; t++)
    {
      nchiT[t] = chiT ;
      nchiS[t] = chiS ;
      nchiL[t] = chiL ;
      nchiR[t] = chiR ;
    }
  neq[0] = nn ;
  neq[1] = ee ;
  neq[2] = qq ;

  marcu[0] = nnmarcu ;
  marcu[1] = eemarcu ;
  Marcu[1] = eeMarcu ;
  marcu[2] = qqmarcu ;

  for (i = 0 ; i < 10 ; i++)
    {
      nn[i] = mxCreate (h, messprintf ("nn_%d", i), MX_COMPLEX, 4, 4, 0) ;
      ee[i] = mxCreate (h, messprintf ("ee_%d", i), MX_COMPLEX, 4, 4, 0) ;
      qq[i] = mxCreate (h, messprintf ("qq_%d", i), MX_COMPLEX, 4, 4, 0) ;
      marcu[0][i] = mxCreate (h, messprintf ("marcu_%d", i), MX_COMPLEX, 4, 4, 0) ;
      marcu[1][i] = mxCreate (h, messprintf ("marcu_%d", i), MX_COMPLEX, 4, 4, 0) ;
      Marcu[1][i] = mxCreate (h, messprintf ("Marcu_%d", i), MX_COMPLEX, 4, 4, 0) ;
      marcu[2][i] = mxCreate (h, messprintf ("marcu_%d", i), MX_COMPLEX, 4, 4, 0) ;
    }
  for (t = 0 ; t < 3 ; t++)
    {
      mxSet (neq[t][1], mu1) ;
      mxSet (neq[t][2], mu2) ;
      mxSet (neq[t][3], mu3) ;
    }
  mxSet (nn[0], mu0n) ;
  if(SU3 == 0)
    mxSet (ee[0], mu0e) ;
  else
    mxSet (ee[0], mu0SU3) ;
  mxSet (qq[0], mu0q) ;

  mxSet (nn[8], mu8n) ;
  mxSet (ee[8], mu8e) ;
  mxSet (qq[8], mu8q) ;

  mxSet (nn[9], mu9n) ;
  mxSet (ee[9], mu9e) ;
  mxSet (qq[9], mu9q) ;

  mxSet (nn[4], mu4n) ;
  mxSet (nn[5], mu5n) ;
  mxSet (nn[6], mu6n) ;
  mxSet (nn[7], mu7n) ;

  mxSet (ee[4], mu4e) ;
  mxSet (ee[5], mu5e) ;
  mxSet (ee[6], mu6e) ;
  mxSet (ee[7], mu7e) ;

  mxSet (qq[4], mu4q) ;
  mxSet (qq[5], mu5q) ;
  mxSet (qq[6], mu6q) ;
  mxSet (qq[7], mu7q) ;


  mxSet (marcu[0][0], marcu0n) ;
  mxSet (marcu[0][4], marcu4n) ;
  mxSet (marcu[0][5], marcu5n) ;
  mxSet (marcu[0][6], marcu6n) ;
  mxSet (marcu[0][7], marcu7n) ;
    
  mxSet (marcu[1][0], marcu0e) ;
  mxSet (marcu[1][4], marcu4e) ;
  mxSet (marcu[1][5], marcu5e) ;
  mxSet (marcu[1][6], marcu6e) ;
  mxSet (marcu[1][7], marcu7e) ;
  mxSet (marcu[1][0], marcu0e) ;
  mxSet (Marcu[1][4], Marcu4e) ;
  mxSet (Marcu[1][5], Marcu5e) ;
  mxSet (Marcu[1][6], Marcu6e) ;
  mxSet (Marcu[1][7], Marcu7e) ;
  

    
  mxSet (marcu[2][0], marcu0q) ;
  mxSet (marcu[2][4], marcu4q) ;
  mxSet (marcu[2][5], marcu5q) ;
  mxSet (marcu[2][6], marcu6q) ;
  mxSet (marcu[2][7], marcu7q) ;
    
  if (1) 
    {
      mxNiceShow (qq[1]) ;
      mxNiceShow (qq[2]) ;
      mxNiceShow (qq[3]) ;
      
      mxNiceShow (ee[6]) ;
      mxNiceShow (ee[7]) ;
      mxNiceShow (ee[0]) ;
      mxNiceShow (qq[6]) ;
      mxNiceShow (qq[7]) ;
      mxNiceShow (qq[0]) ;

      mxNiceShow (marcu[0][0]) ;
      mxNiceShow (marcu[0][4]) ;
      mxNiceShow (marcu[0][5]) ;
      mxNiceShow (marcu[0][6]) ;
      mxNiceShow (marcu[0][7]) ;

	    
      mm = mxMatMult (ee[4],ee[4],h) ;
      mxNiceShow (mm) ;
      
      mm2 = mxMatMult (ee[6],ee[6],h) ;
      mxNiceShow (mm2) ;
    }
} /* muInit */

/*************************************************************************************************/

static MX muComposeMatrix (MX mm, MX a00, MX a01, MX a10, MX a11, complex float x00, complex float x01, complex float x10, complex float x11)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, iMax = 4, iiMax = 8 ;
  int di, dj ;
  const complex float *zc ;
  complex float zz[64] ;

  memset (zz, 0, sizeof (zz)) ;
  if (a00)
    {
      mxValues (a00, 0, 0, &zc) ;
      di = 0 ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = x00 * zc[iMax * i + j] ;
    }
  if (a01)
    {
      mxValues (a01, 0, 0, &zc) ;
      di = 0 ; dj = 4 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = x01 * zc[iMax * i + j] ;
    }
  if (a10) 
    {
      mxValues (a10, 0, 0, &zc) ;
      di = 4 ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = x10 * zc[iMax * i + j] ;
    }
  if (a11)
    {
      mxValues (a11, 0, 0, &zc) ;
      di = 4 ; dj = 4 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = x11 * zc[iMax * i + j] ;
    }

  mxSet (mm, zz) ;
  ac_free (h) ;
  
  return mm ;
} /* muComposeMatrix */

/*************************************************************************************************/

static MX muComposeIntMatrix (int d, MX mm, MX a00, MX a01, MX a10, MX a11, int x00, int x01, int x10, int x11)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, iMax = d, iiMax = 2*d ;
  int di, dj ;
  const int *zi ;
  int zz[4*d*d] ;

  memset (zz, 0, sizeof (zz)) ;
  if (a00)
    {
      mxValues (a00, &zi, 0, 0) ;
      di = 0 ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = x00 * zi[iMax * i + j] ;
    }
  if (a01)
    {
      mxValues (a01, &zi, 0,0) ;
      di = 0 ; dj = d ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = x01 * zi[iMax * i + j] ;
    }
  if (a10) 
    {
      mxValues (a10, &zi, 0, 0) ;
      di = d ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = x10 * zi[iMax * i + j] ;
    }
  if (a11)
    {
      mxValues (a11, &zi, 0, 0) ;
      di = d ; dj = d ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = x11 * zi[iMax * i + j] ;
    }

  mxSet (mm, zz) ;
  ac_free (h) ;
  
  return mm ;
} /* muComposeIntMatrix */

/*************************************************************************************************/

static MX muTripleComposeIntMatrix (int d, MX mm, MX a00, MX a01, MX a02, MX a10, MX a11, MX a12, MX a20, MX a21, MX a22)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, iMax = d, iiMax = 3*d ;
  int di, dj ;
  const int *zi ;
  int zz[9*d*d] ;

  memset (zz, 0, sizeof (zz)) ;
  if (a00)
    {
      mxValues (a00, &zi, 0, 0) ;
      di = 0 ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = zi[iMax * i + j] ;
    }
  if (a01)
    {
      mxValues (a01, &zi, 0,0) ;
      di = 0 ; dj = d ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = zi[iMax * i + j] ;
    }
  if (a02)
    {
      mxValues (a02, &zi, 0,0) ;
      di = 0 ; dj = 2*d ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = zi[iMax * i + j] ;
    }
  if (a10) 
    {
      mxValues (a10, &zi, 0, 0) ;
      di = d ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = zi[iMax * i + j] ;
    }
  if (a11)
    {
      mxValues (a11, &zi, 0, 0) ;
      di = d ; dj = d ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = zi[iMax * i + j] ;
    }
  if (a12)
    {
      mxValues (a12, &zi, 0, 0) ;
      di = d ; dj = 2*d ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = zi[iMax * i + j] ;
    }



   if (a20) 
    {
      mxValues (a20, &zi, 0, 0) ;
      di = 2*d ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = zi[iMax * i + j] ;
    }
  if (a21)
    {
      mxValues (a21, &zi, 0, 0) ;
      di = 2*d ; dj = d ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = 2*zi[iMax * i + j] ;
    }
  if (a22)
    {
      mxValues (a22, &zi, 0, 0) ;
      di = 2*d ; dj = 2*d ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = zi[iMax * i + j] ;
    }

  mxSet (mm, zz) ;
  ac_free (h) ;
  
  return mm ;
} /* muTripleComposeIntMatrix */

/*************************************************************************************************/

static MX muMarcuComposeIntMatrix (int NN, int d, int ii, MX mm, MX mu, MX nu, BOOL new)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, iMax = d, iiMax = NN*d ;
  int marcu, di, dj ;
  const int *zi ;
  int zz[NN*NN*d*d] ;


  memset (zz, 0, sizeof (zz)) ;
  if (mu) /* install in the white diagonals the triangular block copies of mu */
    {
      int w = 1 ;
      int diag ;
      int diagMax = nu ? (new ? 1 : NN) : 1 ; /* if nu==0, just populate the block diagonal: good for SU(2) */

      mxValues (mu, &zi, 0, 0) ;
      for (diag = 0 ; diag < diagMax ; diag += 2)
	{
	  w = 1 ;
	  if (diag == 4) w = 24 ;
	  for (marcu = 0 ; marcu < NN - diag ; marcu++)
	    {
	      di = d * (diag + marcu) ; dj = d * (marcu) ;
	      if (diag == 0) w = 1 ;
	      if (diag == 2 && marcu) w = 4*w ; /* gamma */
	  w = w ;
	      for (i = 0 ; i < iMax ; i++)
		for (j = 0 ; j < iMax ; j++)
		  zz[iiMax * (i + di) + (j + dj)] = w * zi[iMax * i + j] ;
	    }
	}
    }
  if (nu && (!new || ii == 4 || ii == 6)) /* install in the black diagonals the triangular block copies of nu */
    {
      int w = 1 ;
      int diag ;
      int diagMax = nu ? (new  ? 2 : NN) : 1 ; /* if nu==0, just populate the block diagonal */
      mxValues (nu, &zi, 0, 0) ;
      
      for (diag = 1 ; diag < (new ? 2 : diagMax) ; diag += 2)
	{
	  w = w ;
	  if (diag == 3) w = 4 ;
	  for (marcu = 0 ; marcu < NN - diag ; marcu++)
	    {
	      di = d * (diag + marcu) ; dj = d * (marcu) ;
	      if (diag == 1 && marcu) w = 1 * w ;
	      if (diag == 3 && marcu) w = 8 * w ; /* to be determined */
	  w = w ;
	      for (i = 0 ; i < iMax ; i++)
		for (j = 0 ; j < iMax ; j++)
		  zz[iiMax * (i + di) + (j + dj)] = w * zi[iMax * i + j] ;
	    }
	}
    }
  mxSet (mm, zz) ;
  ac_free (h) ;
  
  return mm ;
} /* muMarcuComposeIntMatrix */

/*************************************************************************************************/

/* extract the Hermitian part of a matrix */
static MX muBiHK (MX a, int sign, AC_HANDLE h)
{
  MX mm = mxCreate (h, "muBiHK", MX_COMPLEX, 4, 4, 0) ;
  if (a)
    {
      int i, j, iMax = 4 ;
      MX at = mxMatTranspose (0, a, h)  ;
      const complex float *za ;
      const complex float *zat ;
      complex float zz[16] ;
      
      memset (zz, 0, sizeof (zz)) ;
      mxValues (a, 0, 0, &za) ;
      mxValues (at, 0, 0, &zat) ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iMax * (i) + (j)] = 0.5 * (za [iMax * i + j] + sign * conj(zat [iMax * i + j])) ; 
      mxSet (mm, zz) ;

    }
  return mm ;
} /* muBiHK */

/*************************************************************************************************/

/* extract the anti-Hermitian part of a matrix */
static MX muBiH (MX a, AC_HANDLE h)
{
  return muBiHK (a, 1, h) ;
}
static MX muBiK (MX a, AC_HANDLE h)
{
  return muBiHK (a, -1, h) ;
}

/*************************************************************************************************/
/* Mixing 2 famillies, using a pair of angles
 * alpha and beta
 * alpha is the Hermitian angle, it concerns the down quarks
 * beta is the anti-Hermitian angle, it conscerns the up quarks
 *
 * if alpha = beta, or in the electron case
 *   this is just a change of variables global to all the right states
 *   and the representation remains decomposable
 * theta = alpha - beta   could hopefully be the cabbibo angle
 *   It describes the misalignment of te up/c qarks relative to the down/s quarks
 *   We verify here that the mix matrices represent SU(2/1)
 *   We need to verify that the representatin is indecomposable
 *   It has 2 highest weights u_R and c_R
 * and seems to share d_R + s_R with a same phase ? 
 */
static MX muBiComposeMatrix (MX mm, MX a00, MX a01, MX a10, MX a11, int x00, int x01, int x10, int x11)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, iMax = 4, iiMax = 8 ;
  int di, dj ;
  float pi = 3.1415926535 ;
  float alpha = 1*pi/6 ;
  float beta = 1*pi/4 ;

  complex float x00K = x00 * cos (beta) ;
  complex float x01K = x01 * sin (beta) ;
  complex float x10K = x10 * sin (beta) ;
  complex float x11K = x11 * cos (beta) ;
  complex float x00H = x00 * cos (alpha) ;
  complex float x01H = x01 * sin (alpha) ;
  complex float x10H = x10 * sin (alpha) ;
  complex float x11H = x11 * cos (alpha) ;

  const complex float *zcH ;
  const complex float *zcK ;
  complex float zz[64] ;
  MX a00H = muBiH (a00, h) ;
  MX a01H = muBiH (a01, h) ;
  MX a10H = muBiH (a10, h) ;
  MX a11H = muBiH (a11, h) ;
  MX a00K = muBiK (a00, h) ;
  MX a01K = muBiK (a01, h) ;
  MX a10K = muBiK (a10, h) ;
  MX a11K = muBiK (a11, h) ;

  memset (zz, 0, sizeof (zz)) ;
  if (a00)
    {
      mxValues (a00H, 0, 0, &zcH) ;
      mxValues (a00K, 0, 0, &zcK) ;
      di = 0 ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  {
	    zz[iiMax * (i + di) + (j + dj)] = 
	      x00H * zcH[iMax * i + j] +
	      x00K * zcK[iMax * i + j] ;
	  }
    }
  if (a01)
    {
      mxValues (a01H, 0, 0, &zcH) ;
      mxValues (a01K, 0, 0, &zcK) ;
      di = 0 ; dj = 4 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  {
	    zz[iiMax * (i + di) + (j + dj)] = 
	      x01H * zcH[iMax * i + j] +
	      x01K * zcK[iMax * i + j] ;
	  }
    }
  if (a10) 
    {
      mxValues (a10H, 0, 0, &zcH) ;
      mxValues (a10K, 0, 0, &zcK) ;
      di = 4 ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  {
	    zz[iiMax * (i + di) + (j + dj)] = 
	      x10H * zcH[iMax * i + j] +
	      x10K * zcK[iMax * i + j] ;
	  }
    }
  if (a11)
    {
      mxValues (a11H, 0, 0, &zcH) ;
      mxValues (a11K, 0, 0, &zcK) ;
      di = 4 ; dj = 4 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  {
	    zz[iiMax * (i + di) + (j + dj)] = 
	      x11H * zcH[iMax * i + j] +
	      x11K * zcK[iMax * i + j] ;
	  }
    }

  mxSet (mm, zz) ;
  ac_free (h) ;
  
  return mm ;
} /* muBiComposeMatrix */


/*************************************************************************************************/
/* construct the rotated 8x8 mattrices */
static void muInit2 (AC_HANDLE h)
{
  int i, t ;

  chiT2 = mxCreate (h, "chiT", MX_COMPLEX, 8, 8, 0) ;
  chiS2 = mxCreate (h, "chi", MX_COMPLEX, 8, 8, 0) ;
  chiL2 = mxCreate (h, "chiL", MX_COMPLEX, 8, 8, 0) ;
  chiR2 = mxCreate (h, "chiR", MX_COMPLEX, 8, 8, 0) ;

  muComposeMatrix (chiT2, chiT, 0, 0, chiT, 1, 0, 0, 1) ;
  muComposeMatrix (chiS2, chiS, 0, 0, chiS, 1, 0, 0, 1) ;
  muComposeMatrix (chiL2, chiL, 0, 0, chiL, 1, 0, 0, 1) ;
  muComposeMatrix (chiR2, chiR, 0, 0, chiR, 1, 0, 0, 1) ;

  for (t = 3 ; t < NTYPES ; t++)
    {
      nchiT[t] = chiT2 ;
      nchiS[t] = chiS2 ;
      nchiL[t] = chiL2 ;
      nchiR[t] = chiR2 ;
    }

  for (i = 0 ; i < 10 ; i++)
    { 
      N2[i] = mxCreate (h, messprintf ("N2_%d", i), MX_COMPLEX, 8, 8, 0) ;
      E2[i] = mxCreate (h, messprintf ("E2_%d", i), MX_COMPLEX, 8, 8, 0) ;
      Q2[i] = mxCreate (h, messprintf ("Q2_%d", i), MX_COMPLEX, 8, 8, 0) ;

      N2a[i] = mxCreate (h, messprintf ("N2a_%d", i), MX_COMPLEX, 8, 8, 0) ;
      E2a[i] = mxCreate (h, messprintf ("E2a_%d", i), MX_COMPLEX, 8, 8, 0) ;
      Q2a[i] = mxCreate (h, messprintf ("Q2a_%d", i), MX_COMPLEX, 8, 8, 0) ;

      N2b[i] = mxCreate (h, messprintf ("N2b_%d", i), MX_COMPLEX, 8, 8, 0) ;
      E2b[i] = mxCreate (h, messprintf ("E2b_%d", i), MX_COMPLEX, 8, 8, 0) ;
      Q2b[i] = mxCreate (h, messprintf ("Q2b_%d", i), MX_COMPLEX, 8, 8, 0) ;
    }
  neq[3] = N2 ;
  neq[4] = E2 ;
  neq[5] = Q2 ;

  neq[6] = N2a ;
  neq[7] = E2a ;
  neq[8] = Q2a ;
  
  neq[9] = N2b ;
  neq[10] = E2b ;
  neq[11] = Q2b ;
  
  /* even matrices, same block diagonal */
  for (t = 0 ; t < 3 ; t++)
    for (i = 0 ; i < 10 ; i++)
      {
	if (i > 3 && i < 8) continue ;
	muComposeMatrix (neq[t+3][i], neq[t][i], 0, 0, neq[t][i], 1, 0, 0, 1) ;
	muComposeMatrix (neq[t+6][i], neq[t][i], 0, 0, neq[t][i], 1, 0, 0, 1) ;
	muComposeMatrix (neq[t+9][i], neq[t][i], 0, 0, neq[t][i], 1, 0, 0, 1) ;
      }

  /* odd matrices block diagonal */
  for (i = 4 ; i < 8 ; i++)
    for (t = 0 ; t < 3 ; t++)
      muComposeMatrix (neq[t+3][i], neq[t][i], 0, 0, neq[t][i], 1, 0, 0, 1) ;

  /* odd matrices block diagonal + bottom corner */
    for (t = 0 ; t < 3 ; t++)
      {
	muComposeMatrix (neq[t+6][4], neq[t][4], 0, neq[t][5], neq[t][4], 1, 0, 1, 1) ;
	muComposeMatrix (neq[t+6][5], neq[t][5], 0, neq[t][4], neq[t][5], 1, 0, -1, 1) ;
	muComposeMatrix (neq[t+6][6], neq[t][6], 0, neq[t][7], neq[t][6], 1, 0, 1, 1) ;
	muComposeMatrix (neq[t+6][7], neq[t][7], 0, neq[t][6], neq[t][7], 1, 0, -1, 1) ;
      }      
  /* odd matrices block diagonal + top corner */
    for (t = 0 ; t < 3 ; t++)
      {
	muBiComposeMatrix (neq[t+9][4], neq[t][4], neq[t][5], neq[t][5], neq[t][4],1,1,1,1) ;
	muBiComposeMatrix (neq[t+9][5], neq[t][5], neq[t][4], neq[t][4], neq[t][5],1,-1,-1,1) ;
	muBiComposeMatrix (neq[t+9][6], neq[t][6], neq[t][7], neq[t][7], neq[t][6],1,1,1,1) ;
	muBiComposeMatrix (neq[t+9][7], neq[t][7], neq[t][6], neq[t][6], neq[t][7],1,-1,-1,1) ;
      }      
  return ; 
} /* muInit2 */

/*************************************************************************************************/
  
  exit (0) ;
    printf("#### extract the OSp(2/1) sub-superalgebbra Lepton Cabibbo \n") ;
    if (1)
      { 
	/* we extract the generators F Y H X E of OSp(2/1) from the generators of SU(2/1)
	 * for SU(2/1) 
	 *         X = (4 + i5)/2, Y = (4-i5)/2, Z = (6+i7)/2, T = (6 - i7)/2
	 * Now we extract the OSp generators by projection of the odd generators of the eightfold way adjoint of SU(2/1) on the SU(2) axis
         *         OX = (X+T)/2    OY=(Z-Y)/2
	 * Finally
	 *         OE= OX OX,   OF = - OY OY, OH = - OX OY - OY OX
	 *  we now write that as a program
	*/
	MX mx =  mxCreate (h, "m1", MX_COMPLEX, 8, 8, 0) ;
	MX my =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
	MX mz =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
	MX mt =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;

	MX OX =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
	MX OY =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
	MX OH =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
			
	mx = mxLinearCombine (mx, 0.5, E2a[4], 0.5I, E2a[5], h) ;
	my = mxLinearCombine (my, 0.5, E2a[4], -0.5I, E2a[5], h) ;
	mz = mxLinearCombine (mz, 0.5, E2a[6], 0.5I, E2a[7], h) ;
	mt = mxLinearCombine (mt, 0.5, E2a[6], -0.5I, E2a[7], h) ;

	OX = mxLinearCombine (OX, 1, mx, 1, mt, h) ;
	OY = mxLinearCombine (OY, 1, mz, -1, my, h) ;
	MX OE = mxMatMult (OX, OX, h) ;
	MX OF = mxMatMult (OY, OY, h) ;
	OF = mxLinearCombine (OF, -1, OF, 0, OF, h) ;
	OH = mxLinearCombine (OH, -1, mxMatMult (OX, OY, h), -1, mxMatMult (OY, OX, h), h) ;

	OX->name = "OX" ;
	OY->name = "OY" ;
	OE->name = "OE" ;
	OF->name = "OF" ;
	OH->name = "OH" ;

	mxNiceShow (OX) ;
	mxNiceShow (OY) ;
	mxNiceShow (OH) ;
	mxNiceShow (OE) ;
	mxNiceShow (OF) ;

	/* Casimir  HH + 2 EF + 2 Fe + 2 XY - 2 YX */
	MX K1 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX K2 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX K3 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX K4 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX K5 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX KK =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	
	MX KHH = mxMatMult (OH, OH, h) ;
	MX KEF = mxMatMult (OE, OF, h) ;
	MX KFE = mxMatMult (OF, OE, h) ;
	MX KXY = mxMatMult (OX, OY, h) ;
	MX KYX = mxMatMult (OY, OX, h) ;
	K1 = mxAdd (K1, KEF, KFE, h) ;
	K2 = mxAdd (K2, KHH, K1, h) ;
	K3 = mxAdd (K3, K1, K2, h) ;
	K4 = mxSubstract (KYX, KXY, h) ;
	KK = mxAdd (K5, K3, K4, h) ;

	const complex float *zz4 ;
	complex float zz45[64] ;
	mxValues (K4, 0, 0, &zz4) ;
	memcpy (zz45, zz4, sizeof (zz45)) ;
	for (i = 0 ; i < 8 ; i++)
	  zz45[8*i + i] -= .5 ;
	mxSet (K4, zz45) ;

	KK->name = "Q_Casimir" ;
	KHH->name = "Q_Casimir HH" ;
	K1->name = "Q_Casimir EF-FE" ;
	K2->name = "Q_Casimmir HH + EF+ FE " ;
	K3->name = "Q_Casimmir HH + 2 EF+ FE " ;
	K4->name = "Q_Casimir XY - YX -1/2" ;
	if (0)
	  {
	    mxNiceShow (KHH) ;
	    mxNiceShow (K1) ;
	    mxNiceShow (K2) ;
	    mxNiceShow (K3) ;
	  }
	mxNiceShow (K4) ;
	mxNiceShow (KK) ;

	printf ("Verify that the casimir commutes with X and Y 1\n") ;
	MX CKX = mxMatMult (KK, OX, h) ;
	MX CXK = mxMatMult (OX, KK, h) ;
	MX Com =  mxCreate (h, "[casimir,X]", MX_COMPLEX, 8, 8, 0) ;
	Com = mxSubstract (CKX, CXK, h) ;
	mxNiceShow (Com) ;
	
	printf ("Verify that the S-simir anticommutes with X and Y 2\n") ;
	MX SCKX = mxMatMult (K4, OX, h) ;
	MX SCXK = mxMatMult (OX, K4, h) ;
	MX SCom =  mxCreate (h, "{S-casimir,X}", MX_COMPLEX, 8, 8, 0) ;
	SCom = mxAdd (SCom, SCKX, SCXK, h) ;
	mxNiceShow (SCom) ;
	
	printf ("####### OSp(2/1) Lepton-Cabibbo representation done\n") ;
      }
    printf("#### extract the OSp(2/1) sub-superalgebbra Quark Cabibbo \n") ;
    if (1)
      { 
	/* we extract the generators F Y H X E of OSp(2/1) from the generators of SU(2/1)
	 * for SU(2/1) 
	 *        X = (4 + i5)/2, Y = (4-i5)/2, Z = (6+i7)/2, T = (6 - i7)/2
	 * Now we extract the OSp generators by projection of the odd generators of the eightfold way adjoint of SU(2/1) on the SU(2) axis
         *         OX = (X+T)/2    OY=(Z-Y)/2
	 * Finally
	 *         OE= OX OX,   OF = - OY OY, OH = - OX OY - OY OX
	 * we now write that as a program
	*/
	MX mx =  mxCreate (h, "m1", MX_COMPLEX, 8, 8, 0) ;
	MX my =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
	MX mz =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
	MX mt =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;

	MX OX =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
	MX OY =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
	MX OH =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
			
	mx = mxLinearCombine (mx, 0.5, Q2a[4], 0.5I, Q2a[5], h) ;
	my = mxLinearCombine (my, 0.5, Q2a[4], -0.5I, Q2a[5], h) ;
	mz = mxLinearCombine (mz, 0.5, Q2a[6], 0.5I, Q2a[7], h) ;
	mt = mxLinearCombine (mt, 0.5, Q2a[6], -0.5I, Q2a[7], h) ;

	OX = mxLinearCombine (OX, 1, mx, 1, mt, h) ;
	OY = mxLinearCombine (OY, 1, mz, -1, my, h) ;
	MX OE = mxMatMult (OX, OX, h) ;
	MX OF = mxMatMult (OY, OY, h) ;
	OF = mxLinearCombine (OF, -1, OF, 0, OF, h) ;
	OH = mxLinearCombine (OH, -1, mxMatMult (OX, OY, h), -1, mxMatMult (OY, OX, h), h) ;

	OX->name = "QX" ;
	OY->name = "QY" ;
	OE->name = "QE" ;
	OF->name = "QF" ;
	OH->name = "QH" ;

	mxNiceShow (OX) ;
	mxNiceShow (OY) ;
	mxNiceShow (OH) ;
	mxNiceShow (OE) ;
	mxNiceShow (OF) ;
	
	/* Casimir  HH + 2 EF + 2 Fe + 2 XY - 2 YX */
	MX K1 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX K2 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX K3 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX K4 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX K5 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX KK =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	
	MX KHH = mxMatMult (OH, OH, h) ;
	MX KEF = mxMatMult (OE, OF, h) ;
	MX KFE = mxMatMult (OF, OE, h) ;
	MX KXY = mxMatMult (OX, OY, h) ;
	MX KYX = mxMatMult (OY, OX, h) ;
	K1 = mxAdd (K1, KEF, KFE, h) ;
	K2 = mxAdd (K2, KHH, K1, h) ;
	K3 = mxAdd (K3, K1, K2, h) ;
	K4 = mxSubstract (KYX, KXY, h) ;
	KK = mxAdd (K5, K3, K4, h) ;

	if (0)
	  {
	    const complex float *zz4 ;
	    complex float zz45[64] ;
	    mxValues (K4, 0, 0, &zz4) ;
	    memcpy (zz45, zz4, sizeof (zz45)) ;
	    for (i = 0 ; i < 8 ; i++)
	      zz45[8*i + i] -= 0.5 ;
	    mxSet (K4, zz45) ;
	  }
	else
	  K4 = mxLinearCombine (K4, 1, K4, -.25, KK, h) ;
	KK->name = "Casimir" ;
	KHH->name = "Q_Casimir HH" ;
	K1->name = "Q_Casimir EF-FE" ;
	K2->name = "Q_Casimmir HH + EF+ FE " ;
	K3->name = "Q_Casimmir HH + 2 EF+ FE " ;
	K4->name = "Q_Casimir XY - YX -1/2" ;
	if (0)
	  {
	    mxNiceShow (KHH) ;
	    mxNiceShow (K1) ;
	    mxNiceShow (K2) ;
	  }
	mxNiceShow (K3) ;

	mxNiceShow (K4) ;
	mxNiceShow (KK) ;

	printf ("Verify that the casimir commutes with HHH\n") ;
	MX CKXH = mxMatMult (KK, OH, h) ;
	MX CXKH = mxMatMult (OH, KK, h) ;
	MX ComH =  mxCreate (h, "[casimir,H]", MX_COMPLEX, 8, 8, 0) ;
	Com = mxSubstract (CKXH, CXKH, h) ;
	mxNiceShow (Com) ;
	
	printf ("Verify that the casimir commutes with XXX\n") ;
	MX CKX = mxMatMult (KK, OX, h) ;
	MX CXK = mxMatMult (OX, KK, h) ;
	MX Com =  mxCreate (h, "[casimir,X]", MX_COMPLEX, 8, 8, 0) ;
	Com = mxSubstract (CKX, CXK, h) ;
	mxNiceShow (Com) ;
	
	printf ("Verify that the S-casimir anticommutes with X and Y 3\n") ;
	MX SCKX = mxMatMult (K4, OX, h) ;
	MX SCXK = mxMatMult (OX, K4, h) ;
	MX SCom =  mxCreate (h, "{S-casimir,X}", MX_COMPLEX, 8, 8, 0) ;
	SCom = mxAdd (SCom, SCKX, SCXK, h) ;
	mxNiceShow (SCom) ;
	
	printf ("Compute the square of the S-casimir\n") ;
	MX SC2 = mxMatMult (K4, K4, h) ;
	SC2 = mxLinearCombine (SC2, 8/9.0, SC2, -1, KK, h) ;
	SC2->name = "S-Casimir square" ;
	mxNiceShow (SC2) ;
	
	printf ("Compute the product of the Casimir by the S-casimir Q^3 = Q\n") ;
	MX SC3 = mxMatMult (KK, K4, h) ;
	SC3 = mxLinearCombine (SC3, 1/2.00, SC3, -1, K4, h) ;
	SC3->name = "S-Casimir cube" ;
	mxNiceShow (SC3) ;
	
	printf ("####### OSp(2/1) Quark-Cabibbo representation done\n") ;
      }


/*************************************************************************************************/

/* construct the triple Marcu matrices where the cartan subalgebra is non diagonal with U,V non zero as before 2022_05_05 */
static void muInitNMarcuOld (int a, int b, int NN)
{
  KAS kas, kas2, kasQ ;
  AC_HANDLE h = ac_new_handle () ;
  int i, d ;
  MX *mu, *nu, *QQ ;
  MX qmuY, qmuH, qmuE, qmuF, qmuU, qmuV, qmuW, qmuX, qmuK1, qmuK2 ;

  memset (&kas, 0, sizeof(KAS)) ;
  memset (&kas2, 0, sizeof(KAS)) ;
  memset (&kasQ, 0, sizeof(KAS)) ;
  
  kas.a = a ;
  kas.b = b ;
  kas.h = h ;
  kas2.a = a ;
  kas2.b = b + 1 ; /* b+1 ; */
  kas2.h = h ;
  kasQ.h = h ;

  Kasimirs(1,1, FALSE) ;
  Kasimirs(1,0, FALSE) ;
  KasimirConstructTypicMatrices (&kas, FALSE) ;
  KasimirConstructTypicMatrices (&kas2, FALSE) ;
  kasQ.NN = NN ;
  kasQ.show = TRUE ;
  mu = kas.mu ;
  nu = kas2.mu ; 
  QQ = kasQ.mu = (MX *) halloc (10 * sizeof (MX), kas.h) ;

  kasQ.a = kas.a ;
  kasQ.b = kas.b ;
    
  kasQ.d = d = NN * kas.d ;
  kasQ.d1 = kas.d1 ;
  kasQ.d2 = kas.d2 ;
  kasQ.d3 = kas.d3 ;
  kasQ.d4 = kas.d4 ;
  kasQ.chi = kas.chi ;
	
  kasQ.scale = kas.scale ;
  QQ[0] = qmuY = mxCreate (h,  "qmuY", MX_INT, d, d, 0) ;
  QQ[3] = qmuH = mxCreate (h,  "qmuH", MX_INT, d, d, 0) ;
  QQ[1] = qmuE = mxCreate (h,  "qmuE: E", MX_INT, d, d, 0) ;
  QQ[2] = qmuF = mxCreate (h,  "qmuF", MX_INT, d, d, 0) ;
  QQ[6] = qmuU = mxCreate (h,  "qmuU", MX_INT, d, d, 0) ;
  QQ[7] = qmuV = mxCreate (h,  "qmuV", MX_INT, d, d, 0) ;
  QQ[4] = qmuW = mxCreate (h,  "qmuW", MX_INT, d, d, 0) ;
  QQ[5] = qmuX = mxCreate (h,  "qmuX", MX_INT, d, d, 0) ;
  QQ[8] = qmuK1 = mxCreate (h,  "qmuK1: K1 = {U,V}", MX_INT, d, d, 0) ;
  QQ[9] = qmuK2 = mxCreate (h,  "qmuK2: K2 = {W,X}", MX_INT, d, d, 0) ;
  
  
  /* even and odd matrices, same block diagonal, use nu in the bottom left */
  if (1) /* flip sign of nu[5] */
  {
    nu[5] = mxLinearCombine (nu[5], -1, nu[5], 0, nu[5], h) ;
    nu[7] = mxLinearCombine (nu[7], -1, nu[7], 0, nu[7], h) ;
  }
  for (i = 1 ; i < 4 ; i++) /* block diagonal SU(2) */
  {
    muMarcuComposeIntMatrix (NN, kas.d, i, QQ[i], mu[i],0, FALSE) ;
  }
  for (i = 4 ; i < 8 ; i++) /* triangular odd generators */
  {
    muMarcuComposeIntMatrix (NN, kas.d, i, QQ[i], mu[i],nu[i], FALSE) ;
  }

  QQ[8] = qmuK1 = KasCommut(qmuU,qmuV,1,&kasQ) ;
  QQ[9] = qmuK2 = KasCommut(qmuW,qmuX,1,&kasQ) ;
  QQ[0] = qmuY = mxLinearCombine (qmuY, 1, qmuK1, 1, qmuK2, h) ;
  
  if (1) /* rescale */
  {
    int dd = kasQ.d, d2 = dd*dd ;
     int i, yy[dd*dd], s2 = kasQ.scale ;
     const int *xx ;


     mxValues (qmuK1, &xx, 0, 0) ;
    for (i = 0 ; i < d2 ; i++)
      yy[i] = xx[i]/s2 ;
    mxSet (qmuK1, yy) ;
     mxValues (qmuK2, &xx, 0, 0) ;
    for (i = 0 ; i < d2 ; i++)
      yy[i] = xx[i]/s2 ;
    mxSet (qmuK2, yy) ;

    mxValues (qmuY, &xx, 0, 0) ;
    for (i = 0 ; i < d2 ; i++)
      yy[i] = xx[i]/s2 ;
    mxSet (qmuY, yy) ;
  }

  printf ("###### Marcu\n") ;
  for (i = 0 ; i < 1 ; i++)
    niceIntShow (QQ[i]) ;

  
  MX zUV = mxMatMult (qmuU, qmuV, h) ;
  MX zVU = mxMatMult (qmuV, qmuU, h) ;
  niceIntShow (qmuY) ;
  niceIntShow (qmuU) ;
  niceIntShow (qmuV) ;
  niceIntShow (qmuW) ;
  niceIntShow (qmuX) ;
  niceIntShow (zUV) ;
  niceIntShow (zVU) ;
  
    
  KasimirCheckCommutators (&kasQ) ;

  KasimirLowerMetric (&kasQ) ;
  KasimirUpperMetric (&kasQ) ;
  KasimirUpperTensor (&kasQ) ;
      
  KasimirOperatorK2 (&kasQ) ;
  GhostKasimirOperatorXtilde2 (&kasQ) ;
  GhostKasimirOperatorXtilde2New (&kasQ) ;
  if (0) GhostKasimirOperatorXtilde3 (&kasQ) ;
  GhostKasimirOperatorMinus (&kasQ) ;
  
  if (0) KasimirOperatorK4 (&kasQ) ;

	printf ("Verify that the casimir commutes with H  4\n") ;
	MX CKXH = mxMatMult (kasQ.kas2, qmuH, h) ;
	MX CXKH = mxMatMult (qmuH, kasQ.kas2, h) ;
	MX ComH =  mxCreate (h, "[casimir,X]", MX_COMPLEX, d, d, 0) ;
	ComH = mxSubstract (CKXH, CXKH, h) ;
	mxNiceShow (ComH) ;
	
	printf ("Verify that the casimir commutes with X  4\n") ;
	MX CKX = mxMatMult (kasQ.kas2, qmuX, h) ;
	MX CXK = mxMatMult (qmuX, kasQ.kas2, h) ;
	MX Com =  mxCreate (h, "[casimir,X]", MX_COMPLEX, d, d, 0) ;
	Com = mxSubstract (CKX, CXK, h) ;
	mxNiceShow (Com) ;
	
	printf ("Verify that the S-casimir anticommutes with X and Y 5\n") ;
	MX SCKX = mxMatMult (kasQ.CHI, qmuX, h) ;
	MX SCXK = mxMatMult (qmuX, kasQ.CHI, h) ;
	MX SCom =  mxCreate (h, "{S-casimir,X}", MX_COMPLEX, d, d, 0) ;
	SCom = mxAdd (SCom, SCKX, SCXK, h) ;
	mxNiceShow (SCom) ;
	
	printf ("Compute the square of the S-casimir\n") ;
	MX SC2 = mxMatMult (kasQ.CHI,kasQ.CHI, h) ;
	SC2->name = "S-Casimir square" ;
	mxNiceShow (SC2) ;
	
	printf ("Compute the product of the Casimir by the S-casimir Q^3 = Q\n") ;
	MX SC3 = mxMatMult (kasQ.kas2, kasQ.CHI, h) ;
	SC3->name = "S-Casimir cube" ;
	mxNiceShow (SC3) ;

	if(1)
	  {
	    KasimirUpperTensor (&kasQ) ;
	  }
	if (0 && kasQ.show)
	  KasimirOperatorK3 (&kasQ) ;

  exit(0) ;
  return ;
} /* muInitNMarcuOld */

/*************************************************************************************************/
/* construct the triple Marcu matrices where the cartan subalgebra is non diagonal with U,V non zero as before 2022_05_05 */
static void muInitNMarcu (int a, int b, int NN)
{
  KAS kas, kas2, kasQ ;
  AC_HANDLE h = ac_new_handle () ;
  int i, d ;
  MX *mu, *nu, *QQ ;
  MX qmuY, qmuH, qmuE, qmuF, qmuU, qmuV, qmuW, qmuX, qmuK1, qmuK2 ;

  memset (&kas, 0, sizeof(KAS)) ;
  memset (&kas2, 0, sizeof(KAS)) ;
  memset (&kasQ, 0, sizeof(KAS)) ;
  
  kas.a = a ;
  kas.b = b ;
  kas.h = h ;
  kas2.a = a ;
  kas2.b = b + 1 ; /* b+1 ; */
  kas2.h = h ;
  kasQ.h = h ;

  Kasimirs(1,1, FALSE) ;
  Kasimirs(1,0, FALSE) ;
  KasimirConstructTypicMatrices (&kas, FALSE) ;
  KasimirConstructTypicMatrices (&kas2, FALSE) ;
  kasQ.NN = NN ;
  kasQ.show = TRUE ;
  mu = kas.mu ;
  nu = kas2.mu ; 
  QQ = kasQ.mu = (MX *) halloc (10 * sizeof (MX), kas.h) ;

  kasQ.a = kas.a ;
  kasQ.b = kas.b ;
    
  kasQ.d = d = NN * kas.d ;
  kasQ.d1 = kas.d1 ;
  kasQ.d2 = kas.d2 ;
  kasQ.d3 = kas.d3 ;
  kasQ.d4 = kas.d4 ;
  kasQ.chi = kas.chi ;
	
  kasQ.scale = kas.scale ;
  QQ[0] = qmuY = mxCreate (h,  "qmuY", MX_INT, d, d, 0) ;
  QQ[3] = qmuH = mxCreate (h,  "qmuH", MX_INT, d, d, 0) ;
  QQ[1] = qmuE = mxCreate (h,  "qmuE: E", MX_INT, d, d, 0) ;
  QQ[2] = qmuF = mxCreate (h,  "qmuF", MX_INT, d, d, 0) ;
  QQ[6] = qmuU = mxCreate (h,  "qmuU", MX_INT, d, d, 0) ;
  QQ[7] = qmuV = mxCreate (h,  "qmuV", MX_INT, d, d, 0) ;
  QQ[4] = qmuW = mxCreate (h,  "qmuW", MX_INT, d, d, 0) ;
  QQ[5] = qmuX = mxCreate (h,  "qmuX", MX_INT, d, d, 0) ;
  QQ[8] = qmuK1 = mxCreate (h,  "qmuK1: K1 = {U,V}", MX_INT, d, d, 0) ;
  QQ[9] = qmuK2 = mxCreate (h,  "qmuK2: K2 = {W,X}", MX_INT, d, d, 0) ;
  
  
  /* even and odd matrices, same block diagonal, use nu in the bottom left */
  if (1) /* flip sign of nu[5] */
  {
    nu[5] = mxLinearCombine (nu[5], -1, nu[5], 0, nu[5], h) ;
    nu[7] = mxLinearCombine (nu[7], -1, nu[7], 0, nu[7], h) ;
  }
  for (i = 1 ; i < 4 ; i++) /* block diagonal SU(2) */
  {
    muMarcuComposeIntMatrix (NN, kas.d, i, QQ[i], mu[i],0, TRUE) ;
  }
  if (1) /* redefine the nu matrices */
    {
      int a = kas.a + 1 ;
      int d = kas.d ;
      int d1 = kas.d1 ;
      int d2 = kas.d2 ;
      int d3 = kas.d3 ;
      int i, j, yy[d*d] ;
      
      memset (yy, 0, sizeof(yy)) ;
      nu[6] = mxCreate (h,  "qnuU", MX_INT, d, d, 0) ;
      for (i = 0, j = d1 ; i < d1 ; i++, j++)
	yy[d*j + i] = a ; 
      for (i = 1, j = d1 + d2 ; i < d1 ; i++, j++)
	yy[d*j + i] = -a ; 
      for (i = d1+1, j = d1 + d2 + d3 ; i < d1+d2 ; i++, j++)
	yy[d*j + i] = 1 ; 
      for (i = d1+d2, j = d1 + d2 + d3 ; i < d1+d2+d3 ; i++, j++)
	yy[d*j + i] = 1 ; 
      mxSet (nu[6], yy) ;
      mxShow(kas.mu[1]) ;
      mxShow(nu[6]) ;
      nu[4] = KasCommut(kas.mu[1],nu[6],-1,&kas) ;
      mxShow(nu[4]) ;
      if (0)   exit (0) ;
    }
  for (i = 4 ; i < 8 ; i++) /* triangular odd generators */
  {
    muMarcuComposeIntMatrix (NN, kas.d, i, QQ[i], mu[i],nu[i], TRUE) ;
  }

  QQ[8] = qmuK1 = KasCommut(qmuU,qmuV,1,&kasQ) ;
  QQ[9] = qmuK2 = KasCommut(qmuW,qmuX,1,&kasQ) ;
  QQ[0] = qmuY = mxLinearCombine (qmuY, 1, qmuK1, 1, qmuK2, h) ;
  
  if (1) /* rescale */
  {
    int dd = kasQ.d, d2 = dd*dd ;
     int i, yy[dd*dd], s2 = kasQ.scale ;
     const int *xx ;


     mxValues (qmuK1, &xx, 0, 0) ;
    for (i = 0 ; i < d2 ; i++)
      yy[i] = xx[i]/s2 ;
    mxSet (qmuK1, yy) ;
     mxValues (qmuK2, &xx, 0, 0) ;
    for (i = 0 ; i < d2 ; i++)
      yy[i] = xx[i]/s2 ;
    mxSet (qmuK2, yy) ;

    mxValues (qmuY, &xx, 0, 0) ;
    for (i = 0 ; i < d2 ; i++)
      yy[i] = xx[i]/s2 ;
    mxSet (qmuY, yy) ;
  }

  printf ("###### Marcu\n") ;
  for (i = 0 ; i < 1 ; i++)
    mxNiceShow (QQ[i]) ;

  
  MX zUV = mxMatMult (qmuU, qmuV, h) ;
  MX zVU = mxMatMult (qmuV, qmuU, h) ;
  mxNiceShow (qmuY) ;
  mxNiceShow (qmuU) ;
  mxNiceShow (qmuV) ;
  mxNiceShow (qmuW) ;
  mxNiceShow (qmuX) ;
  mxNiceShow (zUV) ;
  mxNiceShow (zVU) ;
  
    
  KasimirCheckCommutators (&kasQ) ;

  KasimirLowerMetric (&kasQ) ;
  KasimirUpperMetric (&kasQ) ;
  KasimirUpperTensor (&kasQ) ;
      
  KasimirOperatorK2 (&kasQ) ;
  GhostKasimirOperatorXtilde2 (&kasQ) ;
  GhostKasimirOperatorXtilde2New (&kasQ) ;
  if (0) GhostKasimirOperatorXtilde3 (&kasQ) ;
  GhostKasimirOperatorMinus (&kasQ) ;

  
  if (0) KasimirOperatorK4 (&kasQ) ;

	printf ("Verify that the casimir commutes with H  4\n") ;
	MX CKXH = mxMatMult (kasQ.kas2, qmuH, h) ;
	MX CXKH = mxMatMult (qmuH, kasQ.kas2, h) ;
	MX ComH =  mxCreate (h, "[casimir,X]", MX_COMPLEX, d, d, 0) ;
	ComH = mxSubstract (CKXH, CXKH, h) ;
	mxNiceShow (ComH) ;
	
	printf ("Verify that the casimir commutes with X  4\n") ;
	MX CKX = mxMatMult (kasQ.kas2, qmuX, h) ;
	MX CXK = mxMatMult (qmuX, kasQ.kas2, h) ;
	MX Com =  mxCreate (h, "[casimir,X]", MX_COMPLEX, d, d, 0) ;
	Com = mxSubstract (CKX, CXK, h) ;
	mxNiceShow (Com) ;
	
	printf ("Verify that the S-casimir anticommutes with X and Y 5\n") ;
	MX SCKX = mxMatMult (kasQ.CHI, qmuX, h) ;
	MX SCXK = mxMatMult (qmuX, kasQ.CHI, h) ;
	MX SCom =  mxCreate (h, "{S-casimir,X}", MX_COMPLEX, d, d, 0) ;
	SCom = mxAdd (SCom, SCKX, SCXK, h) ;
	mxNiceShow (SCom) ;
	
	printf ("Compute the square of the S-casimir\n") ;
	MX SC2 = mxMatMult (kasQ.CHI,kasQ.CHI, h) ;
	SC2->name = "S-Casimir square" ;
	mxNiceShow (SC2) ;
	
	printf ("Compute the product of the Casimir by the S-casimir Q^3 = Q\n") ;
	MX SC3 = mxMatMult (kasQ.kas2, kasQ.CHI, h) ;
	SC3->name = "S-Casimir cube" ;
	mxNiceShow (SC3) ;

	if(1)
	  {
	    KasimirUpperTensor (&kasQ) ;
	  }
	if (1 && kasQ.show)
	  KasimirOperatorK3 (&kasQ) ;

  exit(0) ;
  return ;
} /* muInitNMarcu */

/*************************************************************************************************/
/* construct the 1 > 3 > 1 < <(1) cycle */
static KAS *cycle (int a, int b)
{
  KAS *kas ;
  AC_HANDLE h = ac_new_handle () ;
  int i, j,  d = 8 ;
  MX *mu ;
  MX muY, muH, muE, muF, muU, muV, muW, muX, muK1, muK2 ;
  int xx[d*d] ;
  const int *xx1 ;
  const int *xx2 ;

  kas = halloc (sizeof(KAS), h) ;
  memset (kas, 0, sizeof(KAS)) ;
  kas->a = 0 ;
  kas->b = 1 ;
  kas->h = h ;
  kas->d = d ;
  kas->d1 = 3 ;   /* states y=0 (universal donor),plus y=2, -2  scalar of the triplets */
  kas->d2 = 2 ;
  kas->d3 = 2 ;
  kas->d4 = 1 ;  /* state 0 universal sink */
  kas->isCycle = TRUE ;
  kas->show = TRUE ;
  mu = kas->mu = (MX *) halloc (10 * sizeof (MX), h) ;

  muY = mxCreate (h,  "muY: Y Hypercharge", MX_INT, d, d, 0) ;
  muH = mxCreate (h,  "muH: Even SU(2) Cartan operator", MX_INT, d, d, 0) ;
  muE = mxCreate (h,  "muE: Even raising operator", MX_INT, d, d, 0) ;
  muF = mxCreate (h,  "muF: Even lowering operator", MX_INT, d, d, 0) ;
  muU = mxCreate (h,  "muU: Odd raising operator", MX_INT, d, d, 0) ;
  muV = mxCreate (h,  "muV: Odd lowering operator", MX_INT, d, d, 0) ;
  muW = mxCreate (h,  "muW: Other odd raising operator", MX_INT, d, d, 0) ;
  muX = mxCreate (h,  "muX: Other odd lowering operator", MX_INT, d, d, 0) ;
  muK1 = mxCreate (h,  "muK1: K1 = {U,V}", MX_INT, d, d, 0) ;
  muK2 = mxCreate (h,  "muK2: K2 = {W,X}", MX_INT, d, d, 0) ;


  memset (xx, 0, sizeof(xx)) ;
  xx[3*d+3] = 1 ;
  xx[4*d+4] = -1 ;
  xx[5*d+5] = 1 ;
  xx[6*d+6] = -1 ;
  mxSet (muH, xx) ;

  memset (xx, 0, sizeof(xx)) ;
  xx[1*d+1] = 2 ;
  xx[2*d+2] = -2 ;
  xx[3*d+3] = 1 ;
  xx[4*d+4] = 1 ;
  xx[5*d+5] = -1 ;
  xx[6*d+6] = -1 ;
  mxSet (muY, xx) ;

  /* even operators */
  memset (xx, 0, sizeof(xx)) ;
  xx[3*d+4] = 1 ;
  xx[5*d+6] = 1 ;
  mxSet (muF, xx) ;

  memset (xx, 0, sizeof(xx)) ;
  xx[4*d+3] = 1 ;
  xx[6*d+5] = 1 ;
  mxSet (muE, xx) ;

  /* odd operators */
  memset (xx, 0, sizeof(xx)) ;
  xx[0*d+5] = -a ;
  xx[1*d+3] = 1 ;
  xx[6*d+2] = -1 ;
  xx[4*d+7] = b ;
  mxSet (muV, xx) ;

  memset (xx, 0, sizeof(xx)) ;
  xx[0*d+4] = -11-a ;
  xx[3*d+1] = 1 ;
  xx[2*d+6] = 1 ;
  xx[5*d+7] = -17-b ;
  mxSet (muU, xx) ;

  /* odd other raising operator */
  muW = KasCommut (muE, muU, -1, kas) ;
  muW->name = "muW" ;
  
  /* odd other oweringing operator */
  muX = KasCommut (muV, muF, -1, kas) ;
  muX->name = "muX" ;
  
  /* odd Cartan operator K1 = diag (a,...2,1,/ a,...2,1,0) */
  /* odd Cartan operator K2 = diag (1,2,...a/0,1,2...a) */
  muK1 = KasCommut (muU, muV, 1, kas) ;
  muK2 = KasCommut (muW, muX, 1, kas) ;
  muK1->name = "muK1" ;
  muK2->name = "muK2" ;
  mxValues (muK1, &xx1, 0, 0) ;
  mxValues (muK2, &xx2, 0, 0) ;
  memset (xx, 0, sizeof(xx)) ;
  for (i = 0 ; i < d ; i++)
    for (j = 0 ; j < d ; j++)
      xx[i*d + j] = (xx1[i*d + j] + xx2[i*d + j]) ;
  mxSet (muY, xx) ;

  mu[0] = muY ; mu[1] = muE ; mu[2] = muF ; mu[3] = muH ;
  mu[4] = muW ; mu[5] = muX ;  mu[6] = muU ; mu[7] = muV ; 
  mu[8] = muK1 ;
  mu[9] = muK2 ;

  for (i = 0 ; i < 10 ; i++)
    mxShow (mu[i]) ;
  
  KasimirCheckCommutators (kas) ;
  return kas ;
} /* cycle */

/*************************************************************************************************/

static void marcuCycle (int nn, int a, int b)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii, d = 8, d2 = 16 ;
  MX *mu ;
  MX muY, muH, muE, muF, muU, muV, muW, muX, muK1, muK2 ;
  int xx[d2*d2] ;
  const int *xx1 ;
  const int *xx2 ;
  const int *xx3 ;
  KAS *kas, *kas1, *kas2, *kas3 ;
  int nn0 = nn ;
  int i, j ;
  
  kas = kas1 = cycle (a, b) ;
  if (nn0 > 1)
    {
      kas2 = cycle (a+1, b+2) ;
      kas3 = cycle (7, 13) ;
      
      
      kas = halloc (sizeof(KAS), h) ;
      memset (kas, 0, sizeof(KAS)) ;
      kas->a = 0 ;
      kas->b = 1 ;
      kas->isCycle = TRUE ;
      kas->h = h ;
      kas->d = 2 * kas1->d ;
      kas->d1 = 2 * kas1->d1 ;
      kas->d2 = 2 * kas1->d2 ;
      kas->d3 = 2 * kas1->d3 ;
      kas->d4 = 2 * kas1->d4 ;
      
      kas->show = TRUE ;
      mu = kas->mu = (MX *) halloc (10 * sizeof (MX), h) ;
      
      muY = mxCreate (h,  "muY: Y Hypercharge", MX_INT, d2, d2, 0) ;
      muH = mxCreate (h,  "muH: Even SU(2) Cartan operator", MX_INT, d2, d2, 0) ;
      muE = mxCreate (h,  "muE: Even raising operator", MX_INT, d2, d2, 0) ;
      muF = mxCreate (h,  "muF: Even lowering operator", MX_INT, d2, d2, 0) ;
      muU = mxCreate (h,  "muU: Odd raising operator", MX_INT, d2, d2, 0) ;
      muV = mxCreate (h,  "muV: Odd lowering operator", MX_INT, d2, d2, 0) ;
      muW = mxCreate (h,  "muW: Other odd raising operator", MX_INT, d2, d2, 0) ;
      muX = mxCreate (h,  "muX: Other odd lowering operator", MX_INT, d2, d2, 0) ;
      muK1 = mxCreate (h,  "muK1: K1 = {U,V}", MX_INT, d2, d2, 0) ;
      muK2 = mxCreate (h,  "muK2: K2 = {W,X}", MX_INT, d2, d2, 0) ;
      
      mu[0] = muY ; mu[1] = muE ; mu[2] = muF ; mu[3] = muH ;
      mu[4] = muW ; mu[5] = muX ;  mu[6] = muU ; mu[7] = muV ; 
      mu[8] = muK1 ;
      mu[9] = muK2 ;
      
      for (ii = 0 ; ii < 10 ; ii++)
	{
	  memset (xx, 0, sizeof(xx)) ;
	  mxValues (kas1->mu[ii], &xx1, 0, 0) ;
	  mxValues (kas2->mu[ii], &xx2, 0, 0) ;
	  mxValues (kas3->mu[ii], &xx3, 0, 0) ;
	  for (i = 0 ; i < d ; i++)
	    for (j = 0 ; j < d ; j++)
	      {
		xx[i*d2 + j] = xx1[i*d + j]  ;
		xx[(i+d)*d2 + (j+d)] = xx2[i*d + j]  ;
		if (ii == 4 || ii == 6)
		  xx[i*d2 + (j+d)] = xx3[i*d + j] ;
		if (ii == 5 || ii == 7)
		  xx[i*d2 + (j+d)] = -xx3[i*d + j] ;
	      }
	  mxSet (mu[ii], xx) ;
	  if (0) mxShow (mu[ii]) ;
	}
      
      kas->mu[8] = muK1 = KasCommut (muU, muV, 1, kas) ;
      kas->mu[9] = muK2 = KasCommut (muW, muX, 1, kas) ;
      mxValues (muK1, &xx1, 0, 0) ;
      mxValues (muK2, &xx2, 0, 0) ;
      memset (xx, 0, sizeof(xx)) ;
      for (i = 0 ; i < d2 ; i++)
	for (j = 0 ; j < d2 ; j++)
	  xx[i*d2 + j] = (xx1[i*d2 + j] + xx2[i*d2 + j]) ;
      mxSet (muY, xx) ;
      
      mxShow (muK1) ;
      mxShow (muK2) ;
      mxShow (muY) ;
      mxShow (muU) ;
      mxShow (muV) ;
      
      if (0)
	for (ii = 0 ; ii < 10 ; ii++)
	  mxShow (mu[ii]) ;
    }
  KasimirCheckCommutators (kas) ;

  KasimirLowerMetric (kas) ;
  KasimirUpperMetric (kas) ;
  KasimirUpperTensor (kas) ;
  
  KasimirOperatorK2 (kas) ;
  GhostKasimirOperatorXtilde2 (kas) ;
  GhostKasimirOperatorXtilde2New (kas) ;
  GhostKasimirOperatorMinus (kas) ;
 if (0)  GhostKasimirOperatorXtilde3 (kas) ;
  exit (0) ;
  return ;
}

/*************************************************************************************************/
/*************************************************************************************************/

void muConjugate (AC_HANDLE h)
{
  int i, t ;
  for (t = 0 ; t < NTYPES ; t++)
    for (i = 0 ; i < 10 ; i++)
      neq[t][i] = muHermite (neq[t][i], h) ;
} /* muConjugate */

/*************************************************************************************************/

static void casimir2 (const char *title)
{
  AC_HANDLE h = ac_new_handle () ;
  MX mm1 = 0, mm2 = 0, chi = 0, casimir = 0 ;
  int i, j, t ;
  float complex z ;
  float a, b ;

  printf ("%s\n", title) ;
  for (t = 0 ; t < NTYPES ; t++)
    {
      casimir  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
      printf ("# Casimir2 Type %d\n", t) ;
      for (i = 0 ; i < 8 ; i++) 
	for (j = 0 ; j < 8 ; j++)
	{
	  if (i < 4 && j >= 4)
	    continue ;
	  if (i >= 4 && j < 4)
	    continue ;

	  /* compute the coefficient g_ab */
	  chi = nchiS[t] ;
	  mm1  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
	  mm2  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
	  mm1 = mxMatMult (neq[t][i], neq[t][j], h) ;
	  mm2 = mxMatMult (chi, mm1, h) ;
	  z = mxMatTrace (mm2) ;
	  z = z/2.0 ;
	  a = creal (z) ;
	  b = cimag (z) ;
	  
	  if (a*a + b*b < .01)
	    continue ;
	  casimir = mxLinearCombine (casimir, 1, casimir, z, mm1, h) ;
	  if (0 && i == 3 && j == 3)
	    mxNiceShow (casimir) ;
	}
      mxNiceShow (casimir) ;
    }
  printf ("\n") ;
  ac_free (h) ;
} /* casimir2 */

/*************************************************************************************************/

static void casimir3 (const char *title, BOOL isHyper)
{
  AC_HANDLE h = ac_new_handle () ;
  MX mm1 = 0, mm1a = 0, mm1b = 0, mm2 = 0, mm3 = 0, mm4 = 0, chi = 0, casimir = 0 ;
  int i, j, k, t ;
  float complex z ;
  float a, b ;

  printf ("%s\n", title) ;
  for (t = 0 ; t < NTYPES ; t++)
    {
      if (myType != -1 && t != myType)
	continue ;
      casimir  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
      printf ("# Casimir3 Type %d\n", t) ;
      for (i = 0 ; i < 8 ; i++) 
	for (j = 0 ; j < 8 ; j++)
	  for (k = 0 ; k < 8 ; k++)
	    {
	      switch (c3Mask)
		{
		case 333:   /* abc => expect zero */
		  if (i*j*k == 0)
		    continue ;
		  if (i > 3 || j > 3 || k > 3)
		    continue ;
		  break ;
		case 888:  /* 888 */
		  if (i+j+k > 0)
		    continue ;
		  break ;
		case 833:  /* 8aa */
		  if (i*j*k > 0)
		    continue ;
		  if (i+j+k == 0 || i > 3 || j > 3 || k > 3)
		    continue ;
		  break ;
		case 844:   /* 8ij */
		  if (i*j*k > 0)
		    continue ;
		  if (i+j+k ==0 || i*(i-4) < 0 || j*(j-4) < 0 || k*(k-4) < 0 )
		    continue ;
		  break ;
		case 344:   /* aij */
		  if (i*j*k == 0)
		    continue ;
		  if (i < 4 && j < 4 && k < 4)
		    continue ;
		  break ;
		case 444:   /* ijk => expect zero */
		  if (i < 4 || j < 4 || k < 4)
		    continue ;
		  break ;
		}

	      /* compute the coefficient t_abc */
	      chi = SU3 ? nchiT[t] : nchiS[t] ;
	      mm1  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
	      mm1a  = mxCreate (h,  "mm2", MX_COMPLEX, ss[t], ss[t], 0) ;
	      mm1b  = mxCreate (h,  "mm2", MX_COMPLEX, ss[t], ss[t], 0) ;
	      mm2  = mxCreate (h,  "mm2", MX_COMPLEX, ss[t], ss[t], 0) ;
	      mm3  = mxCreate (h,  "mm3", MX_COMPLEX, ss[t], ss[t], 0) ;
	      mm1a = mxMatMult (neq[t][j], neq[t][k], h) ;
	      mm1b = mxMatMult (neq[t][k], neq[t][j], h) ;
	      z = 1 ;
	      if (SU3 == 0 && ( j>=4 && k >= 4))
		z = -1 ;
	      mm1 =  mxLinearCombine (mm1, 1, mm1a, z, mm1b , h) ;
	      if (0)
		{
		  mxNiceShow (mm1) ;	
		  mxNiceShow (neq[t][i]) ;	
		}
	      mm2 = mxMatMult (neq[t][i], mm1, h) ;
	      mm3 = mxMatMult (chi, mm2, h) ;
	      z = isHyper ? mxMatTrace (mm2) : mxMatTrace (mm3) ;

	      z *= 9 ;
	      z /= 8 ;
	      if (SU3 == 0) z *= 1 ;
	      if (0) 
		mxNiceShow (mm3) ;	

	      if (t%3 == 2) z *= 3 ; /* 3 quark colors */
	      if (SU3 == 0 && i*j*k == 0 && ! isHyper) z = -z ;/* g^00 and (g^00)cube == -1 */
	      a = creal (z) ;
	      b = cimag (z) ;
	      
	      if (0 && a*a + b*b < .000001) 
		continue ;
	      if (0)
		{
		  printf ("# mm1 Type %d [%d,%d,%d] %.2f %.2f\n", t,i,j,k, a,b) ;
		  mxNiceShow (mm2) ;
		}
	      mm4 = casimir ;
	      casimir  = mxCreate (h,  "casimir3", MX_COMPLEX, ss[t], ss[t], 0) ;
	      casimir = mxLinearCombine (casimir, 1, mm4, z, mm2, h) ;

	      if (SU3 == 1 && j==7 && k == 7)
		{
		  printf ("# Casimir3 Type %d [%d,%d,%d]\n", t,i,j,k) ;
		  mxNiceShow (casimir) ;
		}
	    }
      mxNiceShow (casimir) ;
    }
  printf ("\n") ;
  ac_free (h) ;
} /* casimir3 */

/*************************************************************************************************/

static void mu2p (const char *title)
{
  AC_HANDLE h = ac_new_handle () ;
  MX mm1 = 0, mm2 = 0, chi ;
  int i, j, t, pass, ok ;
  BOOL debug = FALSE ;
  float complex zz[NTYPES] ;

  printf ("%s\n", title) ;
  printf ("# Index\t\tN  \te  \tq  \tf \tCheck\tN2   \tE2   \tQ2   \tf2   \tCheck2\tN2a   \tE2a   \tQ2a   \tf2a   \tCheck2a\tN2b   \tE2b   \tQ2b   \tf2b   \tCheck2b") ;
  for (i = 0 ; i < 8 ; i++) 
    for (j = 0 ; j < 8 ; j++)
      for (pass = 0 ; pass < 2 ; pass++)
	{
	  if (pass == 0)
	    ok = 0;
	  if (pass == 1 && ok == 0)
	    continue ;
	  
	  if (pass == 1)
	    printf ("\n(%d,%d)\t", i, j) ;
	  for (t = 0 ; t < NTYPES ; t++)
	    {
	      float complex z = 0 ;
	      float a, b ;

	      if (i < 4)
		chi = nchiS[t] ;
	      else
		chi = nchiL[t] ;

	      if (debug) mxNiceShow (neq[t][i]) ;
	      if (debug) mxNiceShow (neq[t][j]) ;
	      mm1  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
	      mm2  = mxCreate (h,  "mm2", MX_COMPLEX, ss[t], ss[t], 0) ;

	      mm1 = mxMatMult (chi, neq[t][i], h) ;
	      mm2 = mxMatMult (mm1, neq[t][j], h) ;
	      if (debug) mxNiceShow (mm1) ;
	      if (debug) mxNiceShow (mm2) ;
	      z = mxMatTrace (mm2) ;
	      zz[t] = z ;

	      a = creal (z) ;
	      b = cimag (z) ;
	      if (pass == 0)
		{
		  if (a != 0 || b != 0)
		    ok = 1;
		}
	      else
		{
		  nicePrint ("\t", z) ;
		    if (t%3 == 2)
		      {
			nicePrint ("\t", zz[t-1] + 3 * zz[t]) ;
			nicePrint ("\t", -4 * zz[t-2] + zz[t-1] + 3 * zz[t]) ;
		      }
		}
	    }
	}
  printf ("\n") ;
  ac_free (h) ;
} /* mu2p */

/*************************************************************************************************/

static void mu3p (const char *title, int type)
{
  AC_HANDLE h = ac_new_handle () ;
  MX chi1, chi2 ;
  MX mm1 = 0, mm2 = 0, mm3 = 0 ;
  int i, j, k, t, pass, ok, sign ;
  float complex zz[NTYPES] ;
  BOOL debug = FALSE ;

  printf ("\n%s\n", title) ;
  printf ("# Index\t\tN  \te  \tq  \tf \tCheck\tN2   \tE2   \tQ2   \tf2   \tCheck2\tN2   \tE2   \tQ2   \tf2   \tCheck2a\tN2b   \tE2b   \tQ2b   \tf2b   \tCheck2b") ;
  for (i = 0 ; i < 8 ; i++) 
    for (j = 0 ; j < 8 ; j++)
      for (k = j ; k < 8 ; k++)
	{
	  switch (type)
	    {
	    case 0: /* f-abc */
	    case 1: /* d-abc */
	      if (j < i || i > 3 || j > 3 || k > 3)
		continue ;
	      break ;
	    case 2: /* f-aij vector-scalar */
	    case 20: /* f-aij vector-scalar */
	    case 21: /* f-aij vector-scalar */
	    case 22: /* d-aij vector-scalar */
	    case 23: /* d-aij vector-scalar */
	      if (i > 3 || j <= 3 || k <= 3)
		continue ;
	      break ;
	    case 4: /* f-aij vector-scalar anomaly */
	      if (i > 3 || j <= 3 || k <= 3)
		continue ;
	      break ;
	    case 3: /* f-abi f-ijk should vanish */
	      if (i <= 3 && (j > 3 || k <= 3))
		continue ;
	      if (i > 3 && (j <= 3 || k <= 3))
		continue ;
	      break ;
	    }
	  for (pass = 0 ; pass < 2 ; pass++)
	    {	      
	      if (pass == 0)
		ok = 0;
	      if (pass == 1 && ok == 0)
		continue ;
	      
	      if (pass == 1)
		printf ("\n(%d,%d,%d)\t", i, j, k) ;
	      for (t = 0 ; t < NTYPES ; t++)
		{
		  float complex z = 0 ;
		  float a, b ;
		  if (debug) mxNiceShow (neq[t][i]) ;
		  if (debug) mxNiceShow (neq[t][j]) ;
		  switch (type)
		    {
		    case 0: /* f-abc symmetrize in mu-nu, skew in bc, use trace:  (L+R) (abc - acb) */
		      chi1 = nchiT[t] ;
		      chi2 = nchiT[t] ;
		      sign = -1 ;
		      break ;
		    case 1: /* d-abc skew-symmetrize in mu-nu, sym in bc, use super trace:  (L-R) (abc + acb) */
		      chi1 = nchiS[t] ;
		      chi2 = nchiS[t] ;
		      sign = 1 ;
		      break ;
		    case 2: /* f-aij use Laij - Raji */
		      chi1 = nchiL[t] ;
		      chi2 = nchiR[t] ;
		      sign = -1 ;
		      break ;
		    case 20: /* f-aij use Trace aij - aji, expect zero */
		      chi1 = nchiT[t] ;
		      chi2 = nchiT[t] ;
		      sign = -1 ;
		      break ;
		    case 21: /* f-aij use STrace aij - aji, expect zero */
		      chi1 = nchiS[t] ;
		      chi2 = nchiS[t] ;
		      sign = -1 ;
		      break ;
		    case 22: /* d-aij use Trace aij - aji, expect zero */
		      chi1 = nchiT[t] ;
		      chi2 = nchiT[t] ;
		      sign = +1 ;
		      break ;
		    case 23: /* d-aij use STrace aij - aji, expect zero */
		      chi1 = nchiS[t] ;
		      chi2 = nchiS[t] ;
		      sign = +1 ;
		      break ;
		    case 3: /* should be null */
		      chi1 = nchiT[t] ;
		      chi2 = nchiT[t] ;
		      sign = 0 ;
		      break ;
		    case 4: /* should be null */
		      chi1 = nchiT[t] ;
		      chi2 = nchiT[t] ;
		      sign = -1 ;
		      break ;
		    }

		  mm1  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
		  mm2  = mxCreate (h,  "mm2", MX_COMPLEX, ss[t], ss[t], 0) ;
		  mm3  = mxCreate (h,  "mm3", MX_COMPLEX, ss[t], ss[t], 0) ;

		  mm1 = mxMatMult (chi1, neq[t][i], h) ;
		  mm2 = mxMatMult (mm1, neq[t][j], h) ;
		  mm3 = mxMatMult (mm2, neq[t][k], h) ;
		  if (debug) mxNiceShow (mm1) ;
		  if (debug) mxNiceShow (mm2) ;
		  z = mxMatTrace (mm3) ;
		  
		  mm1 = mxMatMult (chi2, neq[t][i], h) ;
		  mm2 = mxMatMult (mm1, neq[t][k], h) ;
		  mm3 = mxMatMult (mm2, neq[t][j], h) ;
		  z += sign * mxMatTrace (mm3) ;
		  
		  zz[t] = z ; /* memorize, to be able to compute the Family e + 3*q */
		  a = creal (z) ;
		  b = cimag (z) ;
		  if (pass == 0)
		    {
		      if (a != 0 || b != 0)
			ok = 1;
		    }
		  else
		    {
		      nicePrint ("\t", z) ;
		      if (t%3 == 2)  /* compute the family vertex */
			{
			  nicePrint ("\t", zz[t-1] + 3 * zz[t]) ;
			  nicePrint ("\t", -4 * zz[t-2] + zz[t-1] + 3 * zz[t]) ;
			}
		    }
		}
	    }
	}
  printf ("\n\n") ;
  ac_free (h) ;
} /* mu3p */

/*************************************************************************************************/

static float complex tetraTrace (MX chi, int t, int i, int j, int k, int l)
{
  float complex z = 0 ;
  MX mm1, mm2, mm3, mm4 ;
  AC_HANDLE h = ac_new_handle () ;

  mm1  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
  mm2  = mxCreate (h,  "mm2", MX_COMPLEX, ss[t], ss[t], 0) ;
  mm3  = mxCreate (h,  "mm3", MX_COMPLEX, ss[t], ss[t], 0) ;
  mm4  = mxCreate (h,  "mm4", MX_COMPLEX, ss[t], ss[t], 0) ;

  mm1 = mxMatMult (chi, neq[t][i], h) ;
  mm2 = mxMatMult (mm1, neq[t][j], h) ;
  mm3 = mxMatMult (mm2, neq[t][k], h) ;
  mm4 = mxMatMult (mm3, neq[t][l], h) ;
  z = mxMatTrace (mm4) ;
  
  ac_free (h) ;
  return z ;
} /* tetraTrace */

/******************/

static void mu4p (const char *title, int type)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, k, l, t, a, b, c, d, pass, ok, mult ;
  float complex zz[NTYPES] ;
  
  printf ("\n%s\n", title) ;
  printf ("# Index  \t\tN  \te  \tq  \tf \tCheck\tN2   \tE2   \tQ2   \tf2   \tCheck2\tN2   \tE2   \tQ2   \tf2   \tCheck2\tN2a   \tE2a   \tQ2a   \tf2a   \tCheck2a\tN2b   \tE2b   \tQ2b   \tf2b   \tCheck2b") ;
  for (a = 0 ; a < 6 ; a+=1) 
    for (b = 0 ; b < 8 ; b+=1)
      for (c = 0 ; c < 8 ; c+=1)
	for (d = 0 ; d < 8 ; d+=1)
	  {
	    if (0 && (a-b)*(a-c)*(a-d)*(b-c)*(b-d)*(c-d) == 0)
	      continue ;
	    mult = 1 ;
	    switch (type)
	      {
	      case 0: /* K-abcd 4-vectors, [ab] [cd] usual */
	      case 1: /* K-abcd 4-vectors, [ab] {cd} should be zero */
		/* case {ab} {cd} vanishes because of the g-mu,nu symmetries */
	      case 2: /* anomaly, epsilon mu,nu,tho,sigma, use STr and fully anti-sym in [abcd] */
		if (a > 3 || b > 3 || c > 3 || d > 3)
		  continue ;
		if (a > b || c > d)
		  continue ;
		break ;
	      case 3: /* K=abij 2-vectors 2-scalars */
		i = c ; j = d ;
		if (a > 3 || b > 3 || i < 4 || j < 4)
		  continue ;
		if (a > b || i > j)
		  continue ;
		break ;
	      case 4: /* K=ijkl 4-scalars */
		i = a ; j = b ; k = c ; l = d ;
		if (i < 4 || j < 4 || k < 4 || l < 4)
		  continue ;
		if (i>j || k > l)
		  continue ;
		if (i<j)
		  mult *= 2 ;
		if (k<l)
		  mult *= 2 ;
		break ;
	      }
	    for (pass = 0 ; pass < 2 ; pass++)
	      {	      
		if (pass == 0)
		  ok = 0;
		if (pass == 1 && ok == 0)
		  continue ;
		
		if (pass == 1)
		  printf ("\n(%d,%d,%d,%d)\t", a, b, c, d) ;
		for (t = 0 ; t < 3 && t < NTYPES ; t++)
		  {
		    float complex z = 0 ;
		    switch (type)
		      {
			
		      case 0: /* K-abcd 4-vectors, [ab] [cd] usual use trace and skew symmetrize in (ab) and in (cd) */
			z = 0 ;
			z += tetraTrace (nchiT[t], t, a, b, c, d) ;
			z -= tetraTrace (nchiT[t], t, a, b, d, c) ;
			z -= tetraTrace (nchiT[t], t, b, a, c, d) ;
			z += tetraTrace (nchiT[t], t, b, a, d, c) ;
			break ;
		      case 1: /* K-abcd 4-vectors, [ab] {cd}  use trace expect zero */
			/* case {ab} {cd} vanishes because of the g-mu,nu symmetries */
			z = 0 ;
			z += tetraTrace (nchiT[t], t, a, b, c, d) ;
			z += tetraTrace (nchiT[t], t, a, b, d, c) ;
			z -= tetraTrace (nchiT[t], t, b, a, c, d) ;
			z -= tetraTrace (nchiT[t], t, b, a, d, c) ;
			break ;
		      case 2: /* anomaly, epsilon mu,nu,tho,sigma, use STr and fully anti-sym in [abcd] */
			z = 0 ;
			z += tetraTrace (nchiS[t], t, a, b, c, d) ;
			z -= tetraTrace (nchiS[t], t, a, b, d, c) ;
			z -= tetraTrace (nchiS[t], t, a, c, b, d) ;
			z += tetraTrace (nchiS[t], t, a, c, d, b) ;
			z += tetraTrace (nchiS[t], t, a, d, b, c) ;
			z -= tetraTrace (nchiS[t], t, a, d, c, b) ;
			break ;
		      case 3: /* K-abij 2-vectors, 2-scalars terme direct {ab}(Lij+Rji) */
			/* if we use STrace everywhere, we get zero on lepton + quarks */
			z = 0 ;
			z += tetraTrace (nchiL[t], t, a, b, i, j) ;
			z += tetraTrace (nchiL[t], t, b, a, i, j) ;
			z += tetraTrace (nchiR[t], t, a, b, j, i) ;
			z += tetraTrace (nchiR[t], t, b, a, j, i) ;
			
			/* K-abij 2-vectors, 2-scalars terme croise Laibj + Rajbi */
			z += -2 * tetraTrace (nchiL[t], t, a, i, b, j) ;
			z += -2 * tetraTrace (nchiR[t], t, a, j, b, i) ;
			break ;
		      case 4: /* K-ijkl, 4 scalars symmetrize in {kl} : L(ikjl + iljk) */
			      /* use Strace => zero, use nchiR == 1/2 Trace  => Higgs potential */
			z = 0 ;
			z += tetraTrace (nchiS[t], t, i, k, j, l) ;
			z += tetraTrace (nchiS[t], t, i, l, j, k) ;
			z += tetraTrace (nchiS[t], t, j, k, i, l) ;
			z += tetraTrace (nchiS[t], t, j, l, i, k) ;
		      }
		    z = mult * z ;
		    zz[t] = z ; /* memorize, to be able to compute the Family e + 3*q */
		    if (pass == 0)
		      {
			if (creal (z * conj(z)) > .1)
			  ok = 1;
		      }
		    else
		      {
			nicePrint ("\t", z) ;
			if (t == 2 || t == 5 || t == 8)  /* compute the family vertex */
			  {
			    nicePrint ("\t", (zz[1] + 3 * zz[2])*.3/.8) ;
			    nicePrint ("\t", -4 * zz[0] + zz[1] + 3 * zz[2]) ;
			  }
		      }
		  }
	      }
	  }
  
  printf ("\n\n") ;
  ac_free (h) ;
} /* mu4p */


/*************************************************************************************************/
/*************************************************************************************************/

/* check the non Abelian expansion exp(a)exp(b)exp(-b) = exp (b + [a,b] + [a,[a,b]]/2! + [a,[a[a,b]]]/3! ...) */
static POLYNOME expPol (POLYNOME pp, int NN, int sign, AC_HANDLE h)
{
  int i, fac = 1 ;
  POLYNOME ppp, p[NN+2] ;

  pp = expand (pp) ;
  if (1)
    {
      POLYNOME q2 ;
      q2 = polCopy (pp,h) ;
      if (0)
	{
	  printf (".Q2...... expPol") ;
	  showPol (q2) ;
	}
      q2 = limitN (q2, NN-1) ;
      if (0)
	{
	  printf (".Q2..... expPol") ;
	  showPol (q2) ;
	}
    }

  p[0] = newScalar (1,h) ;
  for (i = 1 ; i <= NN ; i++)
    {
      if (i==1)
	p[i] = polProduct (p[i-1], pp,h) ;
      else
	{
	  POLYNOME q1, q2 ;
	  q1 = polCopy (p[i-1],h) ;
	  q1 = limitN (q1, NN-1) ;
	  q2 = polCopy (pp,h) ;
	  if (0)
	    {
	      printf (".Q2..... expPol") ;
	      showPol (q2) ;
	    }
      q2 = limitN (q2, NN-i+1) ;
      if (0)
	{
	  printf (".QQ2..... expPol") ;
	  showPol (q2) ;
	}
      p[i] = polProduct (q1, q2,h) ;
	}
      if (0)
	{
	  printf (".A...... expPol[x^%d]", i) ;
	  showPol (p[i]) ;
	}
      p[i] = expand (p[i]) ;
      if (0)
	{
	  printf (".B...... expPol[x^%d]", i) ;
	  showPol (p[i]) ;
	}
      p[i] = limitN (p[i], NN) ;
      if (0)
	{
	  printf (".C...... expPol[x^%d]", i) ;
	  showPol (p[i]) ;
	}
    }
  p[i] = 0 ;
  for (i = 1 ; i <= NN ; i++)
    {
      fac *= sign * i ;
      polScale (p[i], 1.0/fac)  ;
    }
  ppp = polMultiSum (h,p) ;
  ppp = expand (ppp) ;
  ppp = expand (ppp) ;
  return ppp ;
}

static POLYNOME superCommutator (POLYNOME p1, POLYNOME p2, AC_HANDLE h)
{
  if (!p1 || !p2)
    return 0 ;

  if (p1 && p1->isSum)
    {
      POLYNOME r1 = superCommutator (p1->p1, p2,h) ;
      POLYNOME r2 = superCommutator (p1->p2, p2,h) ;
      
      return polSum (r1, r2,h) ;
    }
  if (p2 && p2->isSum)
    {
      POLYNOME r1 = superCommutator (p1, p2->p1,h) ;
      POLYNOME r2 = superCommutator (p1, p2->p2,h) ;

      return polSum (r1, r2,h) ;
    }

  POLYNOME r1 = polProduct (p1, p2,h) ;
  POLYNOME r2 = polProduct (p2, p1,h) ;
  POLYNOME r3 ;
  
  int sign = -1 ;
  char *u  = r1->tt.x ;
  while (*u)
    {
      char *v  = r2->tt.x ;
      while (*v)
	{
	  if (*u >= 'i' && *u < 'm' && *v == 'x')
	    sign = -sign ;
	  if (*v >= 'i' && *v < 'm' && *u == 'x')
	    sign = -sign ;
	  v++ ;
	}
    }
  
  r3 = r2 ;
  if (sign == -1)
    {
      polScale (r3, -1) ;
    }
  return expand (polSum (r1, r3,h)) ;
} /* superCommutator */

static POLYNOME repeatedSuperCommutator (POLYNOME p1, POLYNOME p2, int NN, AC_HANDLE h)
{

  POLYNOME p3 = p2 ;

  if (NN < 1)
    messcrash ("NN=%d < 1 in repeatedSuperCommutator", NN) ;
  if (NN > 1)
    p3 = repeatedSuperCommutator (p1, p2, NN - 1,h) ;
  return superCommutator (p1, p3,h) ;
} /* repeatedSuperCommutator */

static void superExponential (int NN, int type, int typeb, AC_HANDLE h)
{
  POLYNOME pp, ss, qa,  qb, qc, qa2, qb2,  rr, p[6], q[6], r[6], pa, pb, pc, pa2, pb2 ;

  char *a = "a" ;
  char *b = "b" ;
  char *c = "c" ;
  int n ;
  

  switch (type)
    {
    case 1: a = "i" ; b = "j" ; break ;
    case 2: a = "i" ; b = "ax" ; break ;
    case 3: a = "a" ; b = "i" ; break ;
    default: a = "a" ; b = "b" ; break ;
    }
  
  if (0)
    {
      qa = newScalar (2,h) ;
      qb = newSymbol ("ii",h) ;
      qa = polProduct (qa, qb,h) ;
      showPol (qa) ;
      qb = expand (qa) ;
      showPol (qb) ;
      exit (0) ;
    }


  qa = newSymbol (a,h) ;
  qb = newSymbol (b,h) ;


  p[0] = expPol (qa, NN, 1,h) ;
  printf (" exp(%s) = ", a) ;
  showPol (p[0]) ;

  p[1] = expPol (qb, NN, 1,h) ;
  printf (" exp(%s) = ", b) ;
  showPol (p[1]) ;
  p[2] = expPol (qa, NN, -1,h) ;
  printf (" exp(-%s) = ", a) ;
  showPol (p[2]) ;
  p[3]= 0 ;

  pp = polMultiProduct (h,p) ;
  pp = expand (pp) ;
  pp = expand (pp) ;
  pp = limitN (pp, NN) ;
  pp = expand (pp) ;
  printf (" exp(%s)exp(%s)exp(-%s) = ", a, b, a) ;
  showPol (pp) ;

  r[0] = qb ;
  r[1] = superCommutator (qa, qb,h) ;
  printf ("\n\n[%s,%s] =", a, b) ;
  showPol (r[1]) ;
  polScale (r[1], 1) ;

  int fac = 1 ;
  for (n = 2 ; n <  NN ; n++)
    {
      fac *= n ;
      r[n] = repeatedSuperCommutator (qa, qb, n,h) ;
      printf ("\n\nn=%d [%s,.. [%s,%s]..] =", n, a, a, b) ;
      showPol (r[n]) ;
      polScale (r[n], 1.0/fac) ;
    }

  r[NN] = 0 ;
  rr = polMultiSum (h,r) ;
  printf (" %s + [%s,%s] =", b, a, b) ;  
  showPol (rr) ;
  
  rr = expPol (rr, NN, 1,h) ;
  rr = expand (rr) ;
  rr = limitN (rr, NN) ;
  rr = expand (rr) ;
  printf ("                     exp( %s + [%s,%s]) =\n", b, a, b) ;  r[2]= 0 ;
  showPol (rr) ;
  showPol (pp) ;


  printf ("\n\nexp(%s)exp(%s)exp(%s) - exp( %s + [%s,%s]) =", a, b, a, b, a, b) ;


  polScale (rr, -1) ;
  ss = polSum (pp, rr,h) ;
  ss = expand (ss) ;
  ss = expand (ss) ;


  showPol (ss) ;

  switch (typeb)
    {
    case 1: a="a" ; b = "b" ; c = "k" ; break ;
    case 2: a="a" ; b = "ix" ; c = "c" ; break ;
    case 3: a="a" ; b = "ix" ; c = "k" ; break ;
    case 4: a="i" ; b = "b" ; c = "c" ; break ;
    case 5: a="i" ; b = "b" ; c = "k" ; break ;
    case 6: a="i" ; b = "jx" ; c = "c" ; break ;
    case 7: a="i" ; b = "jx" ; c = "k" ; break ;

    case 8: a="a" ; b = "bx" ; c = "k" ; break ;
    case 9: a="a" ; b = "i" ; c = "c" ; break ;
    case 10: a="a" ; b = "i" ; c = "k" ; break ;
    case 11: a="i" ; b = "bx" ; c = "c" ; break ;
    case 12: a="i" ; b = "bx" ; c = "k" ; break ;
    case 13: a="i" ; b = "j" ; c = "c" ; break ;
    case 14: a="i" ; b = "j" ; c = "k" ; break ;

    default: a="a" ; b = "b" ; c = "c" ; break ;
    }

  qa = newSymbol (a,h) ;
  qb = newSymbol (b,h) ;
  qc = newSymbol (c,h) ;
  
  pa = expPol (qa, NN, 1,h) ;
  pb = expPol (qb, NN, 1,h) ;
  pc = expPol (qc, NN, 1,h) ;
  pa2 = expPol (qa, NN, -1,h) ;
  pb2 = expPol (qb, NN, -1,h) ;

  pp = polProduct (pa2, pc,h) ;   pp = limitN (pp, NN) ;
  pp = polProduct (pp, pa,h) ;   pp = limitN (pp, NN) ;
  pp = polProduct (pb2, pp,h) ;   pp = limitN (pp, NN) ;
  pp = polProduct (pp, pb,h) ;   pp = limitN (pp, NN) ;

  pp = polProduct (pa, pp,h) ;   pp = limitN (pp, NN) ;
  pp = polProduct (pp, pa2,h) ;   pp = limitN (pp, NN) ;
  pp = polProduct (pb, pp,h) ;   pp = limitN (pp, NN) ;
  pp = polProduct (pp, pb2,h) ;   pp = limitN (pp, NN) ;

  printf (".............Holonomy\n") ;
  showPol(pp) ;

  p[0] = polProduct (qa, qb,h) ;
  p[1] = polProduct (qb, qa,h) ;
  polScale (p [1], -1) ;
  ss = polSum (p[0], p[1],h) ; /* commutator [a,b] */
  ss = expand(ss) ;
  printf (".............[%s,%s]\n",a,b) ;
  showPol (ss) ;


  fac = 1 ;
  r[0] = qc ;
  for (n = 1 ; n <  NN ; n++)
    {
      fac *= -n ;
      r[n] = repeatedSuperCommutator (ss, qc, n,h) ;
      printf ("\n\nn=%d [%s,.. [%s,%s]..] =", n, "[]","[[]]", c) ;
      r[n] = limitN (r[n], NN) ;
      polScale (r[n], 1.0/fac) ;
      showPol (r[n]) ;
    }

  r[NN] = 0 ;
  rr = polMultiSum (h,r) ;
  printf ("............... iterated commutator\n") ;
  showPol (rr) ;
  rr = expPol (rr, NN, 1,h) ;
  printf ("...............exp (minus iterated commutator)\n") ;
  showPol (rr) ;

  polScale (rr, -1) ;
  ss = polSum (pp, rr,h) ;
  ss = expand (ss) ;
  printf ("............... holonomy - exp (-[])\n") ;
  showPol (ss) ;
      
  
  exit (0) ;
  
  rr = expand (rr) ;
  showPol (rr) ;
  rr = limitN (rr, NN) ;
  showPol (rr) ;
  exit (0) ;
  
  polScale (rr, -1) ;
  ss = polSum (rr, pp,h) ;
  ss = expand (ss) ;
  showPol (ss) ;
  exit (0) ;
  
  printf (" exp(%s) = ", b) ;
  showPol (p[0]) ;
  p[1] = expPol (qa, NN, 1,h) ;
  printf (" exp(%s) = ", a) ;
  showPol (p[1]) ;
  p[2] = expPol (qc, NN, 1,h) ;
  printf (" exp(%s) = ", c) ;
  showPol (p[2]) ;
  p[3] = expPol (qa2, NN, 1,h) ;
  printf (" exp(-%s) = ", a) ;
  showPol (p[3]) ;
  p[4] = expPol (qb2, NN, 1,h) ;
  printf (" exp(-%s) = ", b) ;
  showPol (p[4]) ;

  p[5] = 0 ;
  pp = polMultiProduct (h,p) ;
  pp = expand (pp) ;
  pp = limitN (pp, NN) ;
  showPol(pp) ;
  polScale (pp, -1) ;
  q[1] = pp ;


  p[0] = polProduct (qa,qb,h) ;
  p[1] = polProduct (qb,qa,h) ;
  showPol (p[1]) ;
  polScale (p[1], -1) ;
  showPol (p[1]) ;

  p[2] = qc ;
  p[3] = 0 ;
  q[2] = polMultiSum (h,p) ;
  q[2] = expand (q[2]) ;
  showPol (q[2]) ;

  q[2] = expPol(q[2], NN, 1,h) ;
  q[2] = expand (q[2]) ;
  showPol (q[2]) ;


  polScale (q[2], -1) ;
  q[3] = 0 ;
  
  pp = polMultiSum (h,q) ;
  pp = expand (pp) ;
  showPol (pp) ;
  }

void pmxSwap (PMX pmx)
{
  int N = pmx ? pmx->N : 0 ;

  POLYNOME q[N*N] ;
  int m, n, sw[N] ;
  for (int i = 0 ; i < N ; i++)
    sw[i] = i ;
  if (N < 4) messcrash ("pmxSwap N=%d < 4", N) ;
  sw[0] = 0 ; sw[1] = 2 ;
  sw[2] = 1 ; sw[3] = 3 ;

  for (int i = 0 ; i < N ; i++)
    for (int j = 0 ; j < N ; j++)
      q[N*i + j] = pmx->pp[N*i + j] ;
  
  for (int i = 0 ; i < N ; i++)
    for (int j = 0 ; j < N ; j++)
      {
	m = sw[i] ;
	n = sw[j] ;
	pmx->pp[N*i + j] = q [N*m + n] ;
      }  
} /* pmxSwap */

/*************************************************************************************/

static void THETA (void)
{
  POLYNOME p1, p2, p3, p21a, p21b, p12a, p12b, p11, p22, p21, p12, det1, det2, det3 ;
  AC_HANDLE h = ac_new_handle () ;

  fprintf (stderr, "### Formal calculations with Grassman variables\n") ;
  fprintf (stderr, "### The hope is to show that det(supergroup SU(2/1)) = 1+alpha Tr(Y)\n") ;


  if (0)   /* check signs in products of grassman */
    {
      char buf[5] ;
      int n = 0, i, j, k, l ;
      for (i = 0 ; i < 4 ; i++)
	for (j = 0 ; j < 4 ; j++)
	  for (k = 0 ; k < 4 ; k++)
	    for (l = 0 ; l < 4 ; l++)
	      {
		buf[0] = 'a'+i ;
		buf[1] = 'a'+j ;
		buf[2] = 'a'+k ;
		buf[3] = 'a'+l ;
		buf[4] = 0 ;
		POLYNOME p = newTheta (buf, h) ;
		POLYNOME p1 = polCopy (p, h) ;
		POLYNOME p2 = expand (p1) ;
		if (p2)
		  {
		    n++ ; 
		    printf ("# %d\t%s\t", n, buf) ;
		    showPol(p2) ;
		  }
	      }
      exit (0) ;		    
    }
  if (0)   /* check signs in determinants */
    {
      PMX px = pmxCreate (2, "test", h) ;
      complex zz[] = {0,1, 1,0, -1} ;
      POLYNOME p = newScalar (1, h) ;
      pmxSet (px, p, zz) ;
      pmxShow (px) ;
      POLYNOME d = pmxDeterminant (px, h) ;
      showPol (d) ; 
      expand (d) ;
      showPol (d) ;
      exit (0) ;		    
    }
  if (0)   /* check signs in determinants */
    {
      PMX px = pmxCreate (4, "test", h) ;
      complex zz[] = {0,1,0,0,  1,0,0,0, 0,0,0,1, 0,0,1,0,   -1} ;
      POLYNOME p = newScalar (1, h) ;
      pmxSet (px, p, zz) ;
      pmxShow (px) ;
      POLYNOME d = pmxDeterminant (px, h) ;
      showPol (d) ;
      expand (d) ;
      showPol (d) ;
      exit (0) ;		    
    }
  p1 = newScalar (1,h) ;
  p2 = newSymbol ("B",h) ;
  strcpy (p2->tt.theta, "vu") ;
  p3 = newSymbol ("B",h) ;
  strcpy (p3->tt.theta, "uv") ;
  p11 = polSum (p1, p2,h) ;
  p22 = polSum (p1, p3,h) ;
  printf("\n### p11 \n")  ;
  showPol (p11) ;
  printf("\n### p22 \n")  ;
  showPol (p22) ;

  p21a = newSymbol ("b",h) ;
  p21a->tt.theta[0] = 'u' ;
  p21a->tt.sqrti = 1 ;
  p21a->tt.z = 1 ;

  p21b = newSymbol ("b",h) ;
  p21b->tt.theta[0] = 'v' ;
  p21b->tt.sqrti = 1 ;
  p21b->tt.z = I ;

  p21 = polSum (p21a, p21b,h) ;
  printf("\n### p21 \n")  ;
  showPol (p21) ;

  p12a = newSymbol ("b",h) ;
  p12a->tt.theta[0] = 'u' ;
  p12a->tt.sqrti = 1 ;
  p12a->tt.z = 1 ;

  p12b = newSymbol ("b",h) ;
  p12b->tt.theta[0] = 'v' ;
  p12b->tt.sqrti = 1 ;
  p12b->tt.z = -I ;

  p12 = polSum (p12a, p12b,h) ;
  printf("\n### p12 \n")  ;
  showPol (p12) ;

  det1 = polProduct (p11, p22,h) ;
  det2 = polProduct (p21, p12,h) ;
  det2->tt.z *= -1 ;
  printf("\n### det1 \n")  ;
  showPol (det1) ;
  printf("\n### det2 \n")  ;
  showPol (det2) ;

  det1 = expand (det1) ;
  det1 = expand (det1) ;
  det2 = expand (det2) ;
  printf("\n### det1 \n")  ;
  showPol (det1) ;
  printf("\n### det2 \n")  ;
  showPol (det2) ;
  
  det2 = expand(det2) ;
  showPol (det2) ;


  det3 = polSum (det1, det2,h) ;

  printf("\n### determinant det(p11,p12)(p21,p22) \n")  ;
  showPol (det3) ;

  printf("\n##########Test the matrix system\n") ;

  BOOL test = FALSE ;
   /****** U ******/
  complex zu1a[] = {0, 1, 0, 0,   1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,   -1} ;
  complex zu1b[] = {0, 1, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,   -1} ;
  PMX u1 = pmxCreate (4, "u1", h) ;
  POLYNOME pu1 = newTheta ("u", h) ;
  if (! test) pu1->tt.x[0] = 'b' ; 
  pu1->tt.sqrti = 1 ;
  if (1) pmxSet (u1, pu1, zu1a) ;
  pmxShow (u1) ;
  
  complex zu2a[] = {0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 1,   0, 0, 1, 0,    -1} ;
  complex zu2b[] = {0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 1,   0, 0, 0, 0,    -1} ;
  PMX u2 = pmxCreate (4, "u2", h) ;
  POLYNOME pu2 = newTheta ("u", h) ;
  if (! test) pu2->tt.x[0] = 'c' ;
  if (test) pu2->tt.z = 1.I ;
  pu2->tt.sqrti = 1 ;
  pmxSet (u2, pu2, zu2a) ;
  pmxShow (u2) ;

  PMX u = pmxSum (u1, u2, "u", h) ;
  pmxShow (u) ;
  
  /****** V ******/
  complex zv1a[] = {0, -1.I, 0, 0,   1.I, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,     -1} ;
  complex zv1b[] = {0, 0, 0, 0,   1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,     -1} ;
  PMX v1 = pmxCreate (4, "v1", h) ;
  POLYNOME pv1 = newTheta ("v", h) ;
  if (! test) pv1->tt.x[0] = 'b' ;
  pv1->tt.sqrti = 1 ;
  pmxSet (v1, pv1, zv1a) ;
  pmxShow (v1) ;

  complex zv2a[] = {0, 0, 0, 0,  0, 0, 0, 0,   0, 0, 0, -1.I,   0, 0, 1.I, 0,    -1} ;
  complex zv2b[] = {0, 0, 0, 0,  0, 0, 0, 0,   0, 0, 0, 0,      0, 0, 1, 0,      -1} ;
  PMX v2 = pmxCreate (4, "v2", h) ;
  POLYNOME pv2 = newTheta ("v", h) ;
  if (! test) pv2->tt.x[0] = 'c' ;
  if (test) pv2->tt.z = 1.I ;
  pv2->tt.sqrti = 1 ;
  pmxSet (v2, pv2, zv2a) ;
  pmxShow (v2) ;

  PMX v = pmxSum (v1, v2, "v", h) ;
  pmxShow (v) ;

  /****** W ******/
  complex zw1a[] = {0, 0, -1, 0,   0, 0, 0, 0,  -1, 0, 0, 0,   0, 0, 0, 0,      -1} ;
  complex zw1b[] = {0, 0, -1, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,      -1} ;
  PMX w1 = pmxCreate (4, "w1", h) ;
  POLYNOME pw1 = newTheta ("w", h) ;
  if (!test) pw1->tt.x[0] = 'b' ;
  pw1->tt.sqrti = 1 ;
  pmxSet (w1, pw1, zw1a) ;
  pmxShow (w1) ;
  
  complex zw2a[] = {0, 0, 0, 0,   0, 0, 0, 1,  0, 0, 0, 0,    0, 1, 0, 0,      -1} ;
  complex zw2b[] = {0, 0, 0, 0,   0, 0, 0, 1,  0, 0, 0, 0,    0, 0, 0, 0,      -1} ;
  PMX w2 = pmxCreate (4, "w2", h) ;
  POLYNOME pw2 = newTheta ("w", h) ;
  if (!test) pw2->tt.x[0] = 'c' ;
  if (test) pw2->tt.z = 1.I ;
  pw2->tt.sqrti = 1 ;
  pmxSet (w2, pw2, zw2a) ;
  pmxShow (w2) ;

  PMX w = pmxSum (w1, w2, "w", h) ;
  pmxShow (w) ;
  
  /****** X ******/
  complex zx1a[] = {0, 0, 1.I, 0,   0, 0, 0, 0,  -1.I, 0, 0, 0,   0, 0, 0, 0,      -1} ;
  complex zx1b[] = {0, 0, 0, 0,   0, 0, 0, 0,       -1, 0, 0, 0,   0, 0, 0, 0,      -1} ;
  PMX x1 = pmxCreate (4, "x1", h) ;
  POLYNOME px1 = newTheta ("z", h) ;
  if (!test) px1->tt.x[0] = 'b' ;
  px1->tt.sqrti = 1 ;
  pmxSet (x1, px1, zx1a) ;
  pmxShow (x1) ;

  complex zx2a[] = {0, 0, 0, 0,   0, 0, 0, -1.I,  0, 0, 0, 0,  0, 1.I, 0, 0,      -1} ;
  complex zx2b[] = {0, 0, 0, 0,   0, 0, 0, 0 ,    0, 0, 0, 0,  0, 1, 0, 0,      -1} ;
  PMX x2 = pmxCreate (4, "x2", h) ;
  POLYNOME px2 = newTheta ("z", h) ;
  if (!test) px2->tt.x[0] = 'c' ;
  if (test) px2->tt.z = 1.I ;
  px2->tt.sqrti = 1 ;
  pmxSet (x2, px2, zx2a) ;
  pmxShow (x2) ;

  PMX x = pmxSum (x1, x2, "x", h) ;
  pmxShow (x) ;
 
  /****** UVWX ******/

  PMX uvwxSet[] = {u, v, w, x, 0} ; /* {u,v,w,x,0} ; */
  PMX uvwxSet1[] = {u1, v1, w1, x1, 0} ; /* {u,v,w,x,0} ; */
  PMX uvwxSet2[] = {u2, v2, w2, x2, 0} ; /* {u,v,w,x,0} ; */ 

  PMX uvwx1 = pmxMultiSum (uvwxSet1, "u1+v1+w1+x1", h) ;
  PMX uvwx2 = pmxMultiSum (uvwxSet2, "u2+v2+w2+x2", h) ;
  PMX uvwx = pmxMultiSum (uvwxSet, "u+v+w+x", h) ;
  pmxShow (uvwx) ;
  pmxSwap (uvwx) ;
  PMX uvexp1 = pmxExponential (uvwx1, "exp(u1+v1+w1+x1)", 6, h) ;
  PMX uvexp2 = pmxExponential (uvwx2, "exp(u2+v2+w2+x2)", 6, h) ;
  PMX uvexp = pmxExponential (uvwx, "exp(u+v+w+x)", 6, h) ;
  pmxShow (uvexp1) ;
  pmxShow (uvexp2) ;
  pmxShow (uvexp) ;
  pmxShow (u1) ;
  pmxShow (v1) ;
  pmxShow (w1) ;
  pmxShow (x1) ;


 
  printf ("Matrix determinant\n") ;
  POLYNOME dd = pmxDeterminant (uvexp, h) ;
  showPol (dd) ;
  dd = expand (dd) ; 
  dd = expand (dd) ;
  dd = expand (dd) ;
  showPol (dd) ;
      
  ac_free (h) ;
  exit (0) ;
} /* THETA */

#endif

/*************************************************************************************/
/***************************** Public interface *************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// bimSU21: Construction of su(2/1) representations and the SU(2/1) coset space\n"
	    "// Authors: Jean Thierry-Mieg, NCBI, 2026\n"
	    "// Purpose\n"
	    "// Construct the matrices of irreducible and indecomposable representations\n"
	    "// Construct the Casimirs, super Casimirs, Gorelik ghost Casimir\n"
	    "// Verify all assertions of the BIM SU(2/1) coset paper on a grid of actual examples\n"
	    "// over-saturating the degre of all equations, hence proved\n"
	    "// \n"
	    "//\n"
	    "// Syntax:\n"
	    "// bimSU21 [options]\n"
	    "//   [] [-h] [-help] [--help] : this message\n"
	    "// A: Representations\n"
	    "//   su21 -a <int> -b <int> [-N <int>]\n"
	    "//     export the matrices, Casimirs and verifications for the module with \n"
	    "//     Dynkin lables (a,b), a positive integer, b signed integer\n"
	    "//     Number of generations N (N >= 2)\n"
	    "//       In theory, b can be any complex number,\n"
	    "//     for numerical convenience, we restrict here to positive signed integers\n"
	    "//     but the formulas like the Casimir eigen values are anlytic in b\n"
	    "//       When a or N are large, many outputs are suppressed, try first a<=3, N<=3\n"
	    "//\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  dna2dna --help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/
/*************************************************************************************/
#endif // JUNK8
int main (int argc, const char **argv)
{
  AC_HANDLE h = ac_new_handle () ;

#ifdef JUNK88
  BIM bim ;

  bim.h = h ;
  freeinit () ;

  
  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  int a = 0, b = 0 ;
  getCmdLineInt (&argc, argv, "-a", &a) ;   /* even Dynkin number */
  getCmdLineInt (&argc, argv, "-b", &b) ;   /* odd Dynkin number */


  if (a==-1) /* a test */
    {
      /* eigen values of the cubic super casimir Kas3, scaled by (a+1)^2 */
      /* they were computed by this program called with params " su21 -a a -b b" */
      int z,z1, a, b ;
      int xx[8][8] = {
		     { 0, 0, 0, 0, 0, 0, 0, 0} ,
		     { 0, 0, -8, -48, -144, -1920, -4200, -8064},
		     {72, 0, -8, 0, -24,-128, -48, -240},
		     {600, 192, 0, -48, -24, 0, -48, -240},
		     {2352,1152,400,0,-144,-128,-48,0},
		     {6480, 3840,1960,720,0,-320, -360, -240},
		     {14520,9600,5832,3072,1176,  0, -600, -768},
		     {28392,20160,13552,8400,4536,1792,0,-1008}
      } ;
      for (a = 0 ; a < 7 ; a++)
	for (b = 1 ; b < 8 ; b++)
	  {
	    /* the polynome z gives the eigen values and reported in the paper su21rep.tex with jarvis */
	    z = 4 * b * (b - a -1) *( 2*b - a - 1) * (2*b - a- 1) ;
	    z1 = z ? z : 1 ;
	    printf ("a=%d b=%d x=%d z=%d x/z=%.2f\n", a, b, xx[b][a], z,  xx[b][a]*1.0/z1) ;

	  }
      if (1) exit (0) ;
    }

  if (a < 0)
    usage ("SU(2) Dynkin weigth a should be a positiver integer") ; 
  
  bim.MU3 = bimConstructKacModule (0, 0,TRUE, h) ;
  if ((a || b))
    {
      /* 2021_03_18 
       * construct the 8 matrices for the generic irreps of su(2/1) with h.w. (a,b)
       * verify all commutations relations
       * compute the casimir tensors and operators 

       */
      bimConstructKacModule (a, b, TRUE, h) ;

      /*
	Kasimirs (1,1, FALSE) ;
	Kasimirs (1,0, FALSE) ;
	Kasimirs (a,b, TRUE) ;
      */
      exit (0) ;
    }
  /* always init, otherwise the gcc linker is unhappy */
  //  if (0) muInit (h) ;   /* init the 4x4 matrices */


  /* Verifications des traces sur la theorie des groupes pour l'article sur les anomalies scalaires */
  if (getCmdLineBool (&argc, argv, "-G"))
    {
      //      muConjugate (h) ;
      
      
      printf ("########## Compute the relevant traces of products of 2,3,4 SU(2/1) matrices\n") ;
      printf ("########## In each case, the trace is computed for the neutral representation (N), then for leptons (e), quarks (q) and family (e+3*q)\n") ;
      printf ("########## The observation is that leptons and quarks have anomalous traces, but they compensate each other\n") ;
      printf ("########## The family trace, one lepton +  quarks, is proportional to the neutral trace\n") ;
      printf ("########## In the last column, we check that S = e + 3*q - 4*n == 0\n") ;
      
      printf ("########## Verify the commutators,   all computed norms should vanish\n");
      // muStructure () ;
    
      
      // if (0) mu2p ("######### Metric\n# For the even generators (a,b=0123), compute the Super-Trace: STr(ab)\n# For the odd generators (i=4567), compute the Left trace: LTr(ij)\n We hope to find the SU(2/1) Super-Killing metric") ;
      
      //if (1) casimir2 ("######### Casimir 2\n# 1/2 g^AB mu_A mu_B,   we hope to find a diagonal matrix") ;
      
      // if (1) casimir3 ("######### Super Casimir 3\n# 1/6 d^ABC mu_A mu_B mu_C,   we hope to find a diagonal matrix", FALSE) ;
      //if (0) casimir3 ("######### Hyper Casimir 3\n# 1/6 d^ABC mu_A mu_B mu_C,   we hope to find a diagonal matrix", TRUE) ;
      
      exit (0) ;
    }

  /* superalgebra Jacobi indentities */
  
  if (0) // jacobi
    {
      /*
	muInit (0) ;
	
      if (0) mu3p ("######### Triple Vector Vertex\n# Lie algebra f-abc vertex,\n# compute the trace anti-symmetrized in bc: Tr(a[bc])\n# we hope to find the Lie algebra f-123 = 4i", 0) ;
      
      if (0) mu3p ("######### Adler-Bardeen Anomalous Triple Vector Vertex\n# d-abc anomalous vertex\n# compute the super-trace symmetrized in bc: STr(a{bc})\n# The anomaly should vanish", 1) ;
	  if (0) mu3p ("######### Vector Scalar Vertex\n# since  i and j are oriented, do not symmetrized in i,j but use LTr(aij)-RTr(aji)\n# We hope to find the super-algebra d-aij\n", 2) ;
      if (0) mu3p ("######### Vector Scalar Vertex\n# use Trace (aij - aji), expect zero in f=famille\n", 20) ;
      if (0) mu3p ("######### Vector Scalar Vertex STr measure\n# use SuperTrace (aij - aji), expect zero in f=famille\n", 21) ;
      if (0) mu3p ("######### Vector Scalar Vertex Tr measure\n# use Trace (aij - aji), expect irregularities\n", 22) ;
      if (0) mu3p ("######### Vector Scalar Vertex STr vertex\n# use STrace (aij + aji), expect universal d_aij\n", 23) ;
      if (0) mu3p ("######### The other types of triple vertices, i.e. f-abi and f-ijk should be zero because they do not conserve the even/odd grading\n", 3) ;
      if (0) mu3p ("######### Vector scalar anomaly, Tr (a [ij]) should vanish\n", 4) ;
     
      
      
      
      printf ("\n######### Four vector vertices\n# The 3 types of (abcd) symmetrisations are implied by the trace on the Pauli matrices of the Fermion loop\n") ;
      if (0) mu4p ("#########  K-abcd 4 vectors\n# [ab] [cd]: standard Lie Algebra vertex g_mn f^m_ab f^n_cd", 0) ;
      if (0) mu4p ("#########  K-abcd 4 vectors\n# [ab] {cd} should vanish", 1) ;
      if (0) mu4p ("#########  K-abcd 4 vectors anomaly\n# [abcd]", 2) ;
      
      printf ("\n######### Two vectors, 2 scalars vertices\n# The scalars are oriented, so we do not symmetrize on (ij)\n") ;
      if (0) mu4p ("#########  K-abij 2 vectors, 2 scalars\n# abij: Symmetize in {ab}, use Lij+Rji\n# Then add the K-aibj Symmetrize in {ab}, use (-2)(L.i.j+R.j.i)", 3) ;
      
      printf ("\n######### Four scalars\n# The scalars are oriented,{ij} incoming, {kl} outcoming\n") ;
      if (1) mu4p ("#########  K-ijkl Symmetrize in {ij} and {kl}, use Likjl + Liljk", 4) ;
      */
    }
#endif
    return 0 ;
}

  

