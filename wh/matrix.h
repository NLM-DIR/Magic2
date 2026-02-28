/*  File: matrix.h
 *  Author: J.Thierry-Mieg, 2018 (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) Public domain, no reservations, no guarantee  of any kind
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: header for matrix.h
 *              A general utility for handling vector, matrices and tensors of any rank
 *              The matrices are typed as BOOL, int, float complex
 *              with an accuracy of 8/16/32/64 bits
 *              In mixed mode operations, the lower accuracy matrix is promoted before any calculation
 *
 *              The rank and dimensions of a matrix  are fixed at creation, they are not elastic.
 *              The matrices form a ring, with Abelian addition and associative multiplication.
 *              Trace and multiplication are implemented as index contractions (with matching rankension)
 *              Generalized transposition and reshapings are implemented.
 *              When aadding a p-tensor to a q-tensor (p>q), broadcastiing is implied, generalizing the rule
 *              that is M is a matric and c a number, M+c has the shape of M with all cells increased by c
 *    
  * Exported functions:
 *              mx... Create, Destroy, Trace, Add, Multiply, Contract, Reshape
 * CAVEAT:
 *              Matrix INVERSION is NOT supported
 * HISTORY:
 * Created: Aug 24 2014 (mieg)
 *-------------------------------------------------------------------
 */

#ifndef DEF_MATRIX_H
#define DEF_MATRIX_H

#include "complex.h"

#ifdef RESTRICT 
#define Restrict restrict
#else
#define Restrict 
#endif

#ifndef MATRIX_IS_DEFINED
#define MATRIX_IS_DEFINED

/* SIYNCHRONIZE VERY CAREFULLY WITH matrix.h matrix.c matrix_slices.c */
/* Type of a matrix */
typedef enum {
  MX_NULL = 0, MX_BOOL, MX_INT, MX_FLOAT, MX_COMPLEX
} MX_TYPE ;

/**********************************************************************/
/* Type of the matrix. i.e. MX_INT, MX_FLOAT ....
 *   
 * Rank of the matrix. i.e. its number of indexes
 * 0: a scaler, 1: a vector, 2: a matrix, p>2 a tensor
 *
 * Shape of the matrix.
 *    A list of integers giving the size for each rank
 * If mxRank(mx)  == 3, iterate i0, i1, i2 according to
 *     0 <= i2 < mxShapes[2]
 * 0: a scalar, 1: a vector, 2: a matrix, p>2 a tensor
 */
/**********************************************************************/
#define MAXRANK 64
typedef struct mxStruct *MX ;
typedef struct mxStruct { 
  const char *name ;
  const MX_TYPE type ;
  const int rank ;
  const int shapes[MAXRANK] ;
  const int slice0, slice1 ;
} MM ;
#endif

/* After type, give a zero termintaed list of integers
 * For example: mxCreate (h, "test", MX_FLOAT, 7, 2, 0) ;
 * create a 7 by 2 matrix
 * 7 will be the fastest moving index
 * The product of 2 matrices contracts by default the first index
 * example (Einstein notation, the i repeated index is summed)
 *      a(j,k,l,m) = b(i,j,k, l) c(i,m)
 * therefore in a product the dim of the first index must match
 *
 * ALTERNATIVE CALL
 *      int shapes[MAXRANK] = {3,4,0,...} ;
 *   or int shapes = other_tensor->shapes ;
 *    mxCreate (h, "test", MX_FLOAT, -999,  shapes);
 *    creates a tensor with the given shapes
 */ 
MX mxCreate(AC_HANDLE h,const char *name, MX_TYPE type, ...) ;

void uMxDestroy (void *va) ;
#define mxDestroy(_mx) {uMxDestroy(_mx);(_mx) = 0;}


/* Rank of the matrix. i.e. its number of indexes
 * 0: a scaler, 1: a vector, 2: a matrix, p>2 a tensor
 */
int mxRank (MX mx) ;
MX_TYPE mxType (MX mx) ;

/* Shape of the matrix.
 * A list of integers giving the size for each rank
 * If mxRank(mx)  == 3, iterate i0, i1, i2 according to
 *     0 <= i2 < mxShape[2]
 * 0: a scalar, 1: a vector, 2: a matrix, p>2 a tensor
 */
const int *mxShape (MX mx) ;
MX  mxTranspose (MX at, MX a, int g1, int g2, AC_HANDLE h) ; /* if at==0, at is created, returns at = a.transposed */
MX  mxFlipTranspose (MX at, MX a, int g1, int g2, AC_HANDLE h) ; /* if at==0, at is created, returns at = a.transposed AND mirror image around central point */

BOOL mxCheck (MX a, const char *nam) ; /* check matrix existence */

float complex nicePrint (const char *prefix, float complex z) ; /* show on stdout */
void mxShow (MX a) ;        /* show on stderr */
void mxNiceShow (MX a) ;        /* show on stderr */
void mxSet (MX a, const void *x) ; /* x is memcpy in a->z, it must have the correct type and size, no check */


/* copy matrix b in matrix a, 
 * if a exists, a ann be must have same shape and type
 * copy content of b in a 
 */
MX mxCopy (MX a, MX b, AC_HANDLE h);

/* reshape a matrix and or change its type
 * if a == NULL, creates a as needed
 * else check a type is equal or more precise than b type
 * check that the product of the a dims is equalt to the product of the b dims
 * returns a ;
 */
MX mxReshape(MX a, MX b, MX_TYPE type, int rank, int i0, int i1, int i2, int i3, AC_HANDLE h) ;

/* retype a matrix to a more accurate type (i.e. promote int to float or to complex)
 * if a == NULL, cretates a as needed
 * else check a type is equal or more precise than b type
 * returns a ;
 */ 
MX mxRetype (MX a, MX b, MX_TYPE type, AC_HANDLE h) ;
/* complex conjugate every element, do not change the shape
 * if a == NULL, cretates a as needed
 * else check a type is equal or more precise than b type
 * returns a = b.conjugate ;
 */ 
MX mxConjugate (MX a, MX b, AC_HANDLE h) ;  

void mxRandomInit (MX a) ;
void mxClear (MX a) ; /* reset data to zero, keep shape and type */
/* Map a function f on each element of a matrix
 * returns af = f(a)
 *   if af == 0, af is created
 *   if af == a, the filtering is applied en place.
 *   function f is applied in parallel to each element of a
 */
MX mxMap (MX af, MX a, MX cache, void f(MX, MX, MX), AC_HANDLE h) ;

/***************************************************************************************/
/* Standard Rank 2 Matrix operation */
/* Matrix addition : 
 * returns a = b + c 
 *   if a == 0, a is created
 *   a, b, c must have same rank and compatible shapes
 */
MX mxAdd (MX a, MX b, MX c, AC_HANDLE h) ;       /* A = B + C */
MX mxSubstract (MX b, MX c, AC_HANDLE h) ;       /* returns  B - C */
/* general linear combination:   A = beta B + gamma C */
MX mxLinearCombine (MX a, complex float beta, MX b, complex float gamma, MX c, AC_HANDLE h) ;

MX mxMultiSum (AC_HANDLE h, MX *ms) ; /* returns A+B+...0 */
/* Standard Matrix product 
 * returns a =  matrix multiply (b , c)
 *   if a == 0, a is created
 *   the b->shape[1] must equal c->shape[0]
 *   a,b,c->rank must be equal to 2
 *   a->shape[0] must be equalt to b->shape[0]
 *   a->shape[1] must be equalt to c->shape[1]
 *
 * the matrix product is as usual the contraction of second b index with the first c index
 */
MX mxMatMult (MX b, MX c, AC_HANDLE h) ;
MX mxMatMultiProduct (AC_HANDLE h, MX *ms) ; /* returns ABC...0 */
MX mxMatTranspose (MX a, MX b, AC_HANDLE h) ; /* standard transposition of a rank2 matrix */
MX mxMatTransposeConjugate (MX a, MX b, AC_HANDLE h) ; /* standard hermitian conjugate of a rank2 matrix */
float complex  mxMatTrace (MX a) ;  /* standard trace of a rank2 square matrix */
BOOL mxMatDeterminant (MX a, int *ip, float *fp, complex float *cp) ;
BOOL mxMatCofactor (MX a, int *ip, float *fp, complex float *cp) ; /* cofactor of element i,j */
MX mxMatInverse (MX a, AC_HANDLE h) ; /* return NULL or the inverse matrix if available */
int mxIntDeterminant (int *aa, int rank) ;
int mxIntInverse (int *ai, int *aa, int rank) ;   /* inverse a matrix given as an int buffer */
void  mxIntMult (int *a, int *b, int *c, int rank) ; /* A=B*C, no-verif integer matrix product */

/***************************************************************************************/

/* single value scalar matrix */
MX mxScalarMatrix (const char *name, float x, AC_HANDLE h) ;
float complex mxTopCorner (MX a) ;
void mxValues (MX mx, const int **zip, const float **zfp, const complex float **zcp) ;

/* Matrix elementwise multiplication 
 * returns ((a)) = ((b * c)) 
 *   if a == 0, a is created
 *   a, b, c must have same rank and same shapes
 */
MX mxElementWiseMultiplication (MX a, MX b, MX c, AC_HANDLE h) ;

/* Matrix dot-multiplication : 
 * returns a = b . c 
 *   if a == 0, a is created
 *   the last index of b must have the same dimension as the first index of c
 * recall that r = mxCreate (type,1,3,0,0,0,h) retuns a 3 dimensional row vector 
 * recall that c = mxCreate (type,2,1,7,0,0,h) retuns a 7 dimensional column vector 
 * recall that m = mxCreate (type,2,3,7,0,0,h) retuns a 3 rows,7 columns matrix
 *  i.e.    r.m is a 7 dim row vector
 *          m.v is a 3 dim column vector
 *          r.m.v is a scalar, i.e. a matrix of mx(type,0,0,0,0)
 */

MX mxContractFirstIndex (MX a, MX b, MX c, AC_HANDLE h) ; 
MX mxContractLastIndex (MX a, MX b, MX c, AC_HANDLE h) ;
MX mxContractFirstNIndices (int n, MX a, MX b, MX c, AC_HANDLE h) ;
MX mxSumFirstIndex (MX a, MX b, AC_HANDLE h) ;  /* the resulting matrix loses its first index */
MX mxSumLastIndex (MX a, MX b, AC_HANDLE h) ;  /* the resulting matrix loses its last index */
MX mxDot (MX a, MX b, MX c, AC_HANDLE h) ; /* alias of  mxContractFirstIndex */
MX mxD2Dot (MX a, MX b, MX c, AC_HANDLE h) ; /* alias of  mxContractFirst2Indices */
MX mxD3Dot (MX a, MX b, MX c, AC_HANDLE h) ; /* alias of  mxContractFirst3Indices */
MX mxConvolution (MX a, MX b, MX W, AC_HANDLE h) ;
MX mxMaxPooling (MX a, MX b, MX W, AC_HANDLE h) ;
MX mxMaxPoolingWithCache (MX a, MX b, MX W, MX cache, AC_HANDLE h) ;
MX mxMaxPoolingBack (MX a, MX b, MX c, MX W, AC_HANDLE h) ;

/* mxMatListMult
 *   multiply a zero terminated list of rank-4 matrices
 * example:   m = mxMaltListMult (h, m1,m2,m3...m77, 0) 
 *    returns m = m1*m2*m3*...*m77, allocated on handle h
*/

MX mxMatListMult (AC_HANDLE h0, MX *ap) ;


/* Matrix product contraction
 * returns a =  contract (b , c) on dimemsion g1,g2 
 *   if a == 0, a is created
 *   the b->shape[g1] must equal c->shape[g2]
 *   a->rank must be equal to b->rank + c->rank - 2
 *   a->shape must be equalt to the union (b->shape, c->shape) removing g1 g2 from b c
 *
 * the dot product is implemented as the contraction of the last b index with the first c index
 */
MX mxContract (MX a, MX b, MX c, int g1, int g2, AC_HANDLE h) ;

/* returns sum_over_ijk... { a[i,j,k...] b[i,j,k...] }, where a and b have same shape */
/* if nGoodp, count the number of coincidences of hihest value in a and b when looping on last (slowest) index */
float complex mxFullContraction (MX a, MX b, int *nGoodp, int *nBadp) ;
/* Frobenius norm: sum aij aij conjugate over all matrix cells */
double mxFNorm (MX a) ; 

/* slice the matrix on its slowest index */
MX mxSliceLastIndex (MX a, MX b, int s0, int s1) ;

/* activation functions */
typedef MX (*MX_ACTIVATION_FUNCTION)(MX af, MX a, MX cache, AC_HANDLE h) ;

MX mxSoftMax (MX af, MX a, MX cache, AC_HANDLE h) ;
MX mxSigmoid (MX af, MX a, MX cache, AC_HANDLE h) ;
MX mxRelu (MX af, MX a, MX cache, AC_HANDLE h) ;
MX mxElu (MX af, MX a, MX cache, AC_HANDLE h) ;
MX mxL3elu (MX af, MX a, MX cache, AC_HANDLE h) ;
MX mxTanh (MX af, MX a, MX cache, AC_HANDLE h) ;
MX mxSoftRelu (MX af, MX a, MX cache, AC_HANDLE h) ;
MX mxLog (MX af, MX a, BOOL isMinus, AC_HANDLE h) ;
BOOL mxTest (void) ; 

#endif
