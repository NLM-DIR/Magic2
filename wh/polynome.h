#ifndef _POLYNOME_H
#define _POLYNOME_H

#include "ac.h"
#include "matrix.h"

#define minAbs 0.0000001
#define INDEXMAX 1000
#define GMAX 256
#define SMAX (256*sizeof(short))
typedef struct termStruct {
  int type ;
  double complex z ;   /* complex scalar multiplier. If zero, the whole TT is NULL */
  int sqrt1,sqrt2 ;        /* holds rational sqrt to avoid rounding errors,    sqrt1=2,  sqrt2=3   => z = sqrt(2/3) exactly */
  short sigma[GMAX] ; /* sigma     matrices : non-commutative list of index "ab" means sigma_a sigma-bar_b */
  short sigB[GMAX] ;  /* sigma-bar matrices : non-commutative list of index "ab" means sigma-bar_a sigma_b */
  char x[GMAX]  ; /* a,b,  i,j  x(meaning chi)  symbol to exponentiate */
  short sqrti ; /* exp(i pi/4) */
  char theta[GMAX] ; /* grassman indices, they anticommute and square to zero */
  short g[GMAX] ;     /* Lorentz metric */
  short gg[GMAX] ;    /* group metric */
  short eps[GMAX] ; /* espislon anti symmetric set of n times 4 indices */ 
  short mm[4][GMAX] ;      /*(k p q r)_mu momenta :  "1 ab" means the product p_a p_b */
                          /*  "1 a"  "2 b" means the product p_a q_b */
  int  denom[4] ;     /* number of terms of the form 1/k^, 1/(k+p)^2, 1/(k+p+q)^2, 1/(k+p+q+r)^2 */
  int Id2 ; /* Pauli identity matrix, needed its value is 2 when we trace */
  short freeIndex[GMAX] ;
  int N ; /* Taylor degree in x  symbol */
  int magic ;
} TT ;

typedef struct polynomeStruct *POLYNOME ;
struct polynomeStruct {
  AC_HANDLE h ;
  BOOL isFlat, isSum, isProduct ;
  POLYNOME p1, p2 ;   /* p1, p2 are always allocated on p->h, so extracting p->p1 is not a public operation, so p1, p2 are really private */
  TT tt ;
  int id ;
  int magic ;
} ;

typedef struct pmxStruct *PMX ;
struct pmxStruct {
  AC_HANDLE h ;
  int N ;
  char *title ;
  POLYNOME *pp ;   /* all pp polynomes are copied on pmx->h, extracting them i forbidden without a copy , the should be private */
  TT tt ;
  int id ;
  int magic ;
} ;


extern short firstDummyIndex ;
short newDummyIndex (void) ;


POLYNOME reduceIndices (POLYNOME pp) ;


POLYNOME expand (POLYNOME pp) ;
void showPol (POLYNOME pp) ;
POLYNOME squareMomentaCleanUp (POLYNOME pp) ;
BOOL freeIndex (POLYNOME pp) ;
void polFree (POLYNOME pp) ;

/**** New Polynomes of all types ***/
POLYNOME newPolynome (AC_HANDLE h) ;
POLYNOME polCopy (POLYNOME p1, AC_HANDLE h) ;
POLYNOME newScalar (complex double z, AC_HANDLE h) ;
POLYNOME newG (short mu, short nu, AC_HANDLE h) ;
POLYNOME newEpsilon (short a, short b, short c, short d, AC_HANDLE h) ;
/* newAG (... ,  0) antisymmetric link 1/2(ac bd - ad bc). Optionally adding the i epsilon 
 * newAG (... , +1) is the self-dual projector P+ = 1/4 ( ac bd - ad bc) + i/2 epsilon(abcd)
 * newAG (... , -1) is the self-dual projector P- = 1/4 ( ac bd - ad bc) - i/2 epsilon(abcd)
 *     WE have (P+)^2 = (P+),   (P-)^2 = (P-),   (P+)(P-) = (P-)(P+) = 0 
 */
POLYNOME newAG (short a, short b, short c, short d, int parity, AC_HANDLE h) ;
POLYNOME newAG6 (short a, short b, short c, short d, short e, short f, AC_HANDLE h) ;
POLYNOME newK (short cc, AC_HANDLE h) ;
POLYNOME newP (short cc, AC_HANDLE h) ;
POLYNOME newQ (short cc, AC_HANDLE h) ;
POLYNOME newR (short cc, AC_HANDLE h) ;
POLYNOME newPQR (int pqr, short mu, AC_HANDLE h) ;
POLYNOME newSigma (short cc, AC_HANDLE h) ;
POLYNOME newSigB (short cc, AC_HANDLE h) ;
POLYNOME newSymbol (char *cp, AC_HANDLE h) ;
POLYNOME newTheta (char *cp, AC_HANDLE h) ;

BOOL polScale (POLYNOME pp, double complex z) ;  /* en place */
POLYNOME polSum (POLYNOME p1, POLYNOME p2, AC_HANDLE h) ;
POLYNOME polProduct (POLYNOME p1, POLYNOME p2, AC_HANDLE h) ;
POLYNOME polMultiSum (AC_HANDLE h, POLYNOME ppp[]) ;
POLYNOME polMultiProduct (AC_HANDLE h, POLYNOME ppp[]) ;

POLYNOME dimIntegral (POLYNOME pp) ;
POLYNOME pauliTrace (POLYNOME pp, AC_HANDLE h) ;
POLYNOME momentaCleanUp (POLYNOME pp, short alpha)  ;
POLYNOME bbCleanUp (POLYNOME pp, short a, short b, short c, short d) ;
POLYNOME contractIndices (POLYNOME pp) ;
POLYNOME sortPol (POLYNOME pp) ;
POLYNOME contractProducts (POLYNOME pp) ;

void pmxShow (PMX pmx) ;
PMX pmxCreate (int N, char *title, AC_HANDLE h) ;
PMX pmxCopy (PMX pmx, char *title, AC_HANDLE h) ;
PMX pmxScalar (int N, char *title, POLYNOME p,  AC_HANDLE h) ;
PMX pmxSum (PMX pmx1, PMX pmx2, char *title, AC_HANDLE h) ;
PMX pmxMultiSum (PMX *pmxs, char *title, AC_HANDLE h) ;
PMX pmxProduct (PMX pmx1, PMX pmx2, char *title, AC_HANDLE h) ;
PMX pmxMultiProduct (PMX *pmxs, char *title, AC_HANDLE h) ;
PMX pmxExponential (PMX pmx, char *title, int level, AC_HANDLE h) ;
BOOL pmxExpand (PMX pmx) ;   /* en place */

POLYNOME pmxDeterminant (PMX pmx, AC_HANDLE h) ;
BOOL pmxSet (PMX pmx, POLYNOME p,  double complex *zz) ;
/* fill a pmx with a copy of a given polynome scaled according to the complex matrix zz 
 * to weakly check that the size of pmx and zz are equal,  zz[pmx->N^2]=-1 is required
 */

#endif
