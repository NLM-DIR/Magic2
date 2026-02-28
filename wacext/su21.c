#include "ac.h"
#include <complex.h>
#include "matrix.h"
#include "polynome.h"

/* Create june 2020
 * edited 27 july 2021
 *
 * My personnal way of computing the renormalisability of gauge theories
 *
 * A diagram is a hooking of subdiagrams, like propagators and vertex
 * It is an ordered product of functions of g_ij, epsilon_ijkl, sigma_i, sigma_bar_i
 * and has dependence of momenta (p,q,r), at the numerator and at the denominator
 * In addition, it contains group matrices and structure constants.
 *
 * The idea is to represent each term as a structure containing such operators
 * We only want to compute the pole part
 * 
 * To reduce the problem to a numeric expression, we estimate the degree in p,q,r
 * and compute the correct number of formal derivatives with repect to p_mu, q_nu, r_rho.
 * Then the loop integral reduces to a polynome in g_ab, p, q, (constr
 * Separately we know the overall degree in sigma, sigma_bar, multiply by sigma and trace
 * We extract the overall epsilon by multiplying by an additional epsilon.
 * 
 * expand: allows to flatten the products of polynomes
 * contractIndices: sort and reduces the Einstein repeated indices and the epsilons
 * PauliTrace:  transforms the ordered products of Pauli matrices into a polynome in g_mn and epsilon
 *
 * Z2: Construct the Feynman diagram of a propagator given the field types
 *        by hooking the relevant Feynam rules
 * Z3, Z4: Idem for 3-fields and 4-fields vertices
 * Ward: Construct the sum of all diagrams that are expected to give a Ward identity 
 *
 * Testing: We must check all subcomponents calculations
 *          and check that we recover the Yang-Mill Ward identites
 *          The idea is that the progam can rerun all checks each time it is modified
 *          So each component comes with its internal checks.
 *
 * Purpose: Prove that QAD is renormalisable.

 * A separate part of the program is meant to check the consistency of reps of SU(2/1)
 * The 8 complex matrices are enumerated, and we then verify the commutators and the anomalies
 * we have explicitly the leptons, the quarks, and several indecomposable reps
 */


/***********************************************************************************************************************************************/
/* find the polynome equal to the sum n=0 to N of the n^k */
static void powerSum (void)
{
  int i, k ;
  double ii, ss[12] ;
  int N = 200000 ;
  int a, b, c, d, e ;
  BOOL ok, done ;

  memset (ss, 0, sizeof (ss)) ;

  for (i=0 ; i < N ; i++)
    {
      ii = 1 ;
      if (20*i % N == 0)       printf ("%d", i) ;
      for (k = 1 ; k < 6 ; k++)
	{
	  ii *= i ;
	  ss[k] += ii ;
	  if (20*i % N == 0) printf ("\t%d\t%.3f", k, exp((k+1)*log(i) - log(ss[k]))) ;
	}
      if (20*i % N == 0) 
	printf ("\n") ;
    }

  printf ("\n\nTEST  sum (n)\n") ;
  N = 10 ;
  memset (ss, 0, sizeof (ss)) ;
  for (a = 0, done = FALSE ; !done && a <= 6 ; a++)
    {
      int s = 0 ;
      for (i=0, ok = TRUE ; ok && i < N ; i++)
	{
	  for (ii = 1, k = 1 ; k <= 1 ; k++)
	    ii *= i ;
	  if (1)
	    {
	      int s2 ;

	      s += ii ;
	      s2 = i * (i+a) ;
	      printf ("test k=%d a=%d i=%d, s=%d s2=%d\n", k, a,i, s,s2) ;
	      if ((a+1) * s != s2)
		ok = FALSE ;
	    }
	}
      done = ok ;
    }
  if (done)
    printf ("Success: k=%d a=%d N=%d\n", k, a, N) ;
  else
    printf ("Failed  k=%d a=%d N=%d\n", k, a, N) ;

  printf ("\n\nTEST  sum (n^2)\n") ;
  N = 10 ;
  memset (ss, 0, sizeof (ss)) ;
  for (a = 0, done = FALSE ; !done && a <= 6; a++)
    for (b = 1 ; !done && b <= 6 ; b++)
      for (c = 0; !done && c <= 6 ; c++)
    {
      int s = 0 ;
      for (i=0, ok = TRUE ; ok && i < N ; i++)
	{
	  ii = 1 ;
	  for (k = 1 ; k <= 2 ; k++)	
	    ii *= i ;
	  if (1)
	    {
	      int s2, sm ;

	      s += ii ;
	      s2 = i * (i+a) * (b * i + c) ;
	      sm = s2 / ((a+1)*(b+c)) ;   
	      printf ("test k=%d a=%d b=%d c=%d i=%d, s=%d sm=%d s2=%d\n", k, a, b, c, i, s, sm, s2) ;
	      if ( (a+1) * (c+b) * s != s2)
		ok = FALSE ;
	    }
	}
      if (ok)
	{
	  done = ok ;
	  goto done2 ;
	  break ;
	}
    }
 done2:
  if (done)
    printf ("Success: k=%d a=%d b=%d c=%d N=%d\n", k, a, b, c, N) ;
  else
    printf ("Failed  k=%d a=%d N=%d\n", k, a, N) ;


  printf ("\n\nTEST  sum (n^3)\n") ;
  N = 10 ;
  memset (ss, 0, sizeof (ss)) ;
  for (a = 0, done = FALSE ; !done && a <= 6; a++)
    for (b = 1 ; !done && b <= 6 ; b++)
      for (c = 0; !done && c <= 6 ; c++)
	for (d = 1 ; !done && d <= 6 ; d++)
	  for (e = 0; !done && e <= 6 ; e++)
	    {
	      int s = 0 ;
	      for (i=0, ok = TRUE ; ok && i < N ; i++)
		{
		  ii = 1 ;
		  for (k = 1 ; k <= 3 ; k++)	
		    ii *= i ;
		  if (1)
		    {
		      int s2, sm ;
	      
		      s += ii ;
		      s2 = i * (i+a) * (b * i + c) * (d * i + e) ;
		      sm = s2 / ((a+1)*(b+c))*(d + e) ;   
		      printf ("test k=%d a=%d b=%d c=%d d=%d e=%d  i=%d, s=%d sm=%d s2=%d\n", k, a, b, c, d, e,  i, s, sm, s2) ;
		      if ( (a+1) * (c+b) * (d + e) * s != s2)
			ok = FALSE ;
		    }
		}
      if (ok)
	{
	  done = ok ;
	  goto done3 ;
	  break ;
	}
    }
 done3:
  if (done)
    printf ("Success: k=%d a=%d b=%d c=%d d=%d e=%d N=%d\n", k, a, b, c,d, e,  N) ;
  else
    printf ("Failed  k=%d a=%d N=%d\n", k, a, N) ;

  return ; 
}  /* powerSum */

/***********************************************************************************************************************************************/
/**************************************** Quantum Field theory rules ***************************************************************************/
/***********************************************************************************************************************************************/
#ifdef JUNK

/* k/sig to g */
static int indexTrace (short *old, short *new, int sign, AC_HANDLE h)
{
  short *cp, buf[GMAX] ;
  int i, j, k ;
  static int level = 0 ;
  Array a = arrayNandleCreate () ;
  level++ ;
  for (i = 1 ; i < GMAX && old[i] ; i++) // which index shall i contract with index 0 
    {
      for (k = 0, j = 1 ; j < GMAX && old[j])
	if (i != j)
	  buf[k++] = old[j] ;
      new[0] = old[0] ;
      new[1] = old[i] ;
      indexTrace (buf, new, sign) ;
    }
  level-- ;

  return 0 ;
}
#endif

/*******************************************************************************************/
/* kill monomes with Taylor symbol degree > NN */ 
static POLYNOME limitN (POLYNOME pp, int NN)
{
  static int level = 0 ;
  if (!pp)
    return 0 ;

  if (! level)
    pp = expand (pp) ;
  if (!pp)
    return 0 ;


  level++ ;
  if (pp->isSum)
    {
      pp->p1 = limitN (pp->p1, NN) ;
      pp->p2 = limitN (pp->p2, NN) ;
    }
  if (! pp->isProduct)
    {
      if (pp->tt.N > NN)
	pp->tt.z = 0 ;
    }

  level-- ;

  if (! level)
    pp = expand (pp) ;
  return pp ;
} /* limitN */

/*******************************************************************************************/
/*******************************************************************************************/
/***********************************************************************************************************************************************/
/* New vertex A B HB and A H BB */
static POLYNOME vertex_A_B_HB (short mu, short a, short b, int mm[4], AC_HANDLE h) /* A_mu B_a_b, momentum of the incoming photon */
{
  POLYNOME pp, ppp[6] ; 
  int nn = 0 ; 
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  int z = 0 ; /* -1 */
  BOOL useProjector = TRUE ;
  if (1) useProjector = FALSE ;
	   ;
  nn = 0 ;
  if (! useProjector) { c = a ; d = b ; }
  if (mm[0]) { pp = newK (c,h) ; pp->tt.z = mm[0] ; ppp[nn++] = pp ; }
  if (mm[1]) { pp = newP (c,h) ; pp->tt.z = mm[1] ; ppp[nn++] = pp ; }
  if (mm[2]) { pp = newQ (c,h) ; pp->tt.z = mm[2] ; ppp[nn++] = pp ; }
  if (mm[3]) { pp = newR (c,h) ; pp->tt.z = mm[3] ; ppp[nn++] = pp ; }
  ppp[nn++] = 0 ;

  pp = polMultiSum (h, ppp) ;
  nn = 0 ;
  ppp[nn++] = newScalar (I, h) ;
  ppp[nn++] = pp ;
  ppp[nn++] = newG (mu, d, h) ;
  if (useProjector) ppp[nn++] = newAG (a,b,c,d,z,h) ;
  ppp[nn++] = 0 ;

  
  pp = polMultiProduct (h,ppp) ;

  return pp ; 
} /* vertex_A_B_HB */

/**************************************************/

static POLYNOME vertex_A_H_BB (short mu, short a, short b, int mm[4], AC_HANDLE h) /* momentum of the incoming mu-photon */
{
  POLYNOME pp, ppp[6] ; 
  int nn = 0 ; 
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  int z = 0 ;
  BOOL useProjector = TRUE ;
  if (1) useProjector = FALSE ;
  nn = 0 ;
  if (! useProjector) { c = a ; d = b ; }
  if (mm[0]) { pp = newK (c,h) ; pp->tt.z = mm[0] ; ppp[nn++] = pp ; }
  if (mm[1]) { pp = newP (c,h) ; pp->tt.z = mm[1] ; ppp[nn++] = pp ; }
  if (mm[2]) { pp = newQ (c,h) ; pp->tt.z = mm[2] ; ppp[nn++] = pp ; }
  if (mm[3]) { pp = newR (c,h) ; pp->tt.z = mm[3] ; ppp[nn++] = pp ; }
  ppp[nn++] = 0 ;

  pp = polMultiSum (h, ppp) ;
  nn = 0 ;
  ppp[nn++] = newScalar (I,h) ;
  ppp[nn++] = pp ;
  ppp[nn++] = newG (mu, d,h) ;
  if (useProjector) ppp[nn++] = newAG (a,b,c,d,z,h) ;
  ppp[nn++] = 0 ;
  
  pp = polMultiProduct (h,ppp) ;
  
  return pp ;
} /* vertex_A_H_BB  */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/*   A_mu diffusion of B_ab -> BB_cd */
static POLYNOME vertex_A_B_BB (short mu, short a, short b, short c, short d, int  mm1[4], int  mm2[4], AC_HANDLE h)   /* B(mm1)->BB(mm2) in the orientation of the BB->B line */
{
  POLYNOME pp, ppp[6] ; 
  int i, nn = 0 ; 

  for (i = 0 ; i < 4 ; i++)
    {
      if (mm1[i])
	{
	  pp = newScalar (1 * mm1[i],h) ;
	  pp->tt.mm[i][0] = a ;
	  pp->tt.g[0] = mu ;
	  pp->tt.g[1] = c ;
	  pp->tt.g[2] = b ;
	  pp->tt.g[3] = d ;
	  ppp[nn++] = pp ;	  
	}
      if (mm2[i])
	{
	  pp = newScalar (-1 * mm2[i],h) ;
	  pp->tt.mm[i][0] = c ;
	  pp->tt.g[0] = mu ;
	  pp->tt.g[1] = a ;
	  pp->tt.g[2] = b ;
	  pp->tt.g[3] = d ;
	  ppp[nn++] = pp ;	  
	}
    }
  ppp[nn++] = 0 ;	  
  pp = polMultiSum (h,ppp) ;
  pp->tt.z *= I ;

  return pp ;
} /* vertex_A_B_BB */

/***********************************************************************************************************************************************/
/*   A_mu(p) diffusion of A_nu(q)\, A_rho(r) */
static POLYNOME vertex_A_A_A (short a, short b, short c, int  p[4], int  q[4], int r[4], AC_HANDLE h)   /* all 3 are incoming momenta */
{
  POLYNOME p1, ppp[20] ;
  int nn = 0 ;
  
  for (int i = 0 ; i < 4 ; i++)
    {
      if (p[i])
	{
p1 = ppp[nn++] = newG (a,b,h) ; 
	  p1->tt.mm[i][0] = c ;
	  p1->tt.z = p[i] ;
	  p1->tt.z *= I ;
	  p1 = ppp[nn++] = newG (c,a,h) ; 
	  p1->tt.mm[i][0] = b ;
	  p1->tt.z = -p[i] ;
	  p1->tt.z *= I ;
	}
      if (q[i])
	{
	  p1 = ppp[nn++] = newG (b,c,h) ; 
	  p1->tt.mm[i][0] = a ;
	  p1->tt.z = q[i] ;
	  p1->tt.z *= I ;
	  p1 = ppp[nn++] = newG (a,b,h) ; 
	  p1->tt.mm[i][0] = c ;
	  p1->tt.z = -q[i] ;
	  p1->tt.z *= I ;
	}
      if (r[i])
	{
	  p1 = ppp[nn++] = newG (c,a,h) ; 
	  p1->tt.mm[i][0] = b ;
	  p1->tt.z = r[i] ;
	  p1->tt.z *= I ;
	  p1 = ppp[nn++] = newG (b,c,h) ; 
	  p1->tt.mm[i][0] = a ;
	  p1->tt.z = -r[i] ;
	  p1->tt.z *= I ;
	}

    }
  ppp[nn] = 0 ;
  POLYNOME pp = polMultiSum (h, ppp) ;
  pp = expand (pp) ;
  return pp ;
}

/**********************************************************************************************************************************************/
/*   A_mu diffusion of H -> HB */
static POLYNOME vertex_A_H_HB (short mu, int mm[4], AC_HANDLE h)  /* 2k+p = (2,1,0,0) : sum of the momentum of the H->HB in that direction */
{
  POLYNOME pp, ppp[6] ; 
  int i, nn = 0 ; 

  for (i = 0 ; i < 4 ; i++)
    if (mm[i])
      {
	pp = newScalar (mm[i],h) ;
	pp->tt.mm[i][0] = mu ;
	ppp[nn++] = pp ;	  
      }
  ppp[nn++] = 0 ;	  
  pp = polMultiSum (h,ppp) ;
  pp->tt.z *= I ;
  return pp ;
} /* vertex_A_H_HB */

/**********************************************************************************************************************************************/
/*   A_mu diffusion of c  -> cB ghost */
static POLYNOME vertex_A_c_cB (short mu, int mm[4], AC_HANDLE h)  /* k+p = (1,1,0,0) : incoming momentum of the cB */
{
  POLYNOME pp, ppp[6] ; 
  int i, nn = 0 ; 

  for (i = 0 ; i < 4 ; i++)
    if (mm[i])
      {
	pp = newScalar (mm[i],h) ;
	pp->tt.z *= -I ;        /* -1 because of the incoming momentum */
	pp->tt.mm[i][0] = mu ;
	ppp[nn++] = pp ;	  
      }
  ppp[nn++] = 0 ;	  
  pp = polMultiSum (h,ppp) ;
  return pp ;
} /* vertex_A_c_cB */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/* This vertex is by itelf self dual */
static POLYNOME vertex_B_PsiR_PsiLB (short a, short b, AC_HANDLE h)
{
  int X = 0 ; /* -1 B is anti-self-dual */
  short mu = newDummyIndex() ;
  short nu = newDummyIndex() ;
  POLYNOME p1 = newSigB (mu,h) ;
  POLYNOME projector = newAG(a,b,mu,nu,X,h) ;
  p1->tt.sigB[1] = nu ;
  p1->tt.z = .5 ;
  p1->tt.z *= I ;
  p1->tt.sqrt1 = 1 ;
  p1->tt.sqrt2 = 1 ;

  return polProduct (projector, p1,h) ; 
}

/***********************************************************************************************************************************************/
/* This vertex is by itelf anti self dual */
static POLYNOME vertex_BB_PsiL_PsiRB (short a, short b, AC_HANDLE h)
{
  int X = 0 ; /* 1 : Bbar is self-dual */ ;
  short mu = newDummyIndex() ;
  short nu = newDummyIndex() ;
  POLYNOME p1 = newSigma (mu,h) ;
  POLYNOME projector = newAG(a,b,mu,nu,X,h) ;
  p1->tt.sigma[1] = nu ;
  p1->tt.z = .5 ;
  p1->tt.z *= I ;
  p1->tt.sqrt1 = 1 ;
  p1->tt.sqrt2 = 1 ;

  return polProduct (projector, p1,h) ; 
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME vertex_A_PsiR_PsiRB (short mu, AC_HANDLE h)
{
  POLYNOME p = newSigma (mu,h) ;
  p->tt.z *= I ;
  return p ;
}

/***********************************************************************************************************************************************/

static POLYNOME vertex_A_PsiL_PsiLB (short mu, AC_HANDLE h)
{
POLYNOME p = newSigB (mu,h) ;
  p->tt.z *= I ;
  return p ;
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/


static POLYNOME vertex_H_PsiR_PsiLB (AC_HANDLE h)
{
  POLYNOME p = newScalar (1,h) ;
  p->tt.z *= I ;  /* 2*I/3 */
  p->tt.sqrt1 = 1 ;
  p->tt.sqrt2 = 1 ;
  return p ;
}

/***********************************************************************************************************************************************/

static POLYNOME vertex_HB_PsiL_PsiRB (AC_HANDLE h)
{
  POLYNOME p = newScalar (1,h) ;
  p->tt.z *= I ;
  p->tt.sqrt1 = 1 ;
  p->tt.sqrt2 = 1 ;
  return p ; 
}

/***********************************************************************************************************************************************/

static POLYNOME prop_BB_B (short mu, short nu, short rho, short sig, int pqr, AC_HANDLE h)
{
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  int z = -1 ; /* 0: no epsilon, 1:self dual, -1:anti self, 2 just epsilon */
  int u = -1 ;
  int pqrD = pqr ;
  int pqrN = pqr ;
  if (pqr == 99)   /* contruct the  lagrangian */
    { z = 0 ; pqrD = pqrN = 0 ; u = 1 ; }
  if (pqr == 20) 
    { pqrD = 2 ; pqrN = 0 ; }
  POLYNOME p1 = newAG (mu,nu,a,b, z,h) ; /* z */
  if (0) p1 = newEpsilon (mu, nu, a,b,h) ;
  POLYNOME p2 = newPQR (pqrN, a,h) ;
  POLYNOME p3 = newPQR (pqrN, c,h) ;
  POLYNOME p4 = newG  (b, d,h) ;
  POLYNOME p5 = newAG (c,d,rho,sig, -z,h) ; /* -z */
  if (0) p5 = newEpsilon (c,d,rho,sig,h) ;
  POLYNOME pp, ppp[] = {p1,p2,p3, p4, p5, 0} ; 
  /* POLYNOME pp, ppp[] = {p1,p31,p41,p5, 0} ;  */

  p4->tt.denom[pqrD] = 2 ;
  p4->tt.z *= 16 ;
  if (0) p4->tt.z *= u*I ;
  pp = polMultiProduct (h, ppp) ;

  return pp ;
}

/***********************************************************************************************************************************************/

static POLYNOME prop_AA (short mu, short nu, int pqr, AC_HANDLE h)
{
  POLYNOME pp = newG (mu,nu, h) ;
  pp->tt.denom[pqr] = 1 ;
  pp->tt.z *= I ;
  return pp ;
}

/***********************************************************************************************************************************************/
#ifdef JUNK
static POLYNOME prop_AR (short mu, short nu, int pqr, AC_HANDLE h)
{
  POLYNOME pp = newG (mu,nu, h) ;
  pp->tt.denom[pqr] = 1 ;
  pp->tt.z = -I ; /* negative norm */
  return pp ;
}
#endif
/***********************************************************************************************************************************************/

static POLYNOME prop_HB_H (int pqr, AC_HANDLE h) 
{
  POLYNOME pp = newScalar (1.0, h) ;
  pp->tt.denom[pqr] = 1 ;
  pp->tt.z *= I ;
  return pp ;
}

/***********************************************************************************************************************************************/

static POLYNOME prop_cB_c (int pqr, AC_HANDLE h) 
{
  POLYNOME pp = newScalar (1.0, h) ;
  pp->tt.denom[pqr] = 1 ;
  pp->tt.z *= I ;
  return pp ;
}

/***********************************************************************************************************************************************/

static POLYNOME prop_PsiRB_PsiR (int pqr, AC_HANDLE h)
{
  short cc = newDummyIndex () ;
  POLYNOME p1 = newSigB (cc, h) ;
  p1->tt.z *= I ;
  p1->tt.denom[pqr] = 1 ;
  POLYNOME p2 = newPQR (pqr,cc,h) ;
  POLYNOME pp = contractIndices(polProduct (p1, p2,h)) ;
  return contractProducts (pp) ;
}

/***********************************************************************************************************************************************/

static POLYNOME prop_PsiLB_PsiL (int pqr, AC_HANDLE h)
{
  short cc = newDummyIndex () ;
  POLYNOME p1 = newSigma (cc,h) ;
  p1->tt.z *= I ;
  p1->tt.denom[pqr] = 1 ;
  POLYNOME p2 = newPQR (pqr,cc,h) ;
  POLYNOME pp = contractIndices(polProduct (p1, p2,h)) ;
  return contractProducts (pp) ;
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME pauliCleanUp (POLYNOME pp, short w)
{
  if (! pp)
    return 0 ;
  if (pp->isSum)
    {
      pp->p1 = pauliCleanUp (pp->p1, w) ;
      pp->p2 = pauliCleanUp (pp->p2, w) ;
      pp = expand (pp) ;
    }
  else if (pp->tt.type)
    {
      AC_HANDLE h = pp->h ;
      POLYNOME p1 = pp->tt.sigma[0] || (pp->p1 && pp->p1->tt.sigma[0]) ? newSigB (w,h) : newSigma (w,h) ;
      POLYNOME p11 = pp->tt.sigma[0] || (pp->p1 && pp->p1->tt.sigma[0]) ? newSigma (w,h) : newSigB (w,h) ;
      POLYNOME p2 = polProduct (p1, pp,h) ;
      POLYNOME p4 = pauliTrace (p2, 0) ;
      POLYNOME p6 = contractIndices (p4) ;
      POLYNOME p7 = polProduct (p11, p6,h) ;
      POLYNOME p8 = expand (p7) ;
      POLYNOME p9 = contractIndices (p8) ;
      POLYNOME p10 = expand (p9) ;
      if (p10) p10->tt.z /= 2 ;
      p10 = reduceIndices (p10) ;
      return p10 ;
    }
  return pp ;
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME Z2_AA__loopH  (const char *title) 
{
  AC_HANDLE h = ac_new_handle () ;
  firstDummyIndex = 'a' ;
    
  short mu = newDummyIndex () ;
  short nu = newDummyIndex () ;
  int ppv[4] = {2,1,0,0} ; /* 2k + p : vertex */
  
  POLYNOME p1 = vertex_A_H_HB (mu, ppv,h) ; /*(2k + p)_mu */
  POLYNOME p2 = prop_HB_H (1,h) ;   /* (1/(p+^2 */
  POLYNOME p3 = vertex_A_H_HB (nu, ppv,h) ;
  POLYNOME p4 = prop_HB_H (0,h) ;   /* (1/(k)^2 */
  POLYNOME ppp[] = {p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(polMultiProduct (h, ppp)) ;

  printf ("%s\n", title) ;
  showPol (pp) ;

  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  squareMomentaCleanUp (pp) ;
  printf ("### Z2 AA loop H expect ::  1/3 (p_ab - g_ab p^2)\n") ;
  showPol (pp) ;

  printf ("DONE %s\n", title) ;

  return pp ;
} /* Z2_AA__loopH */

/***********************************************************************************************************************************************/

static POLYNOME Z2_AA__loopA  (const char *title) 
{
  AC_HANDLE h = ac_new_handle () ;
  firstDummyIndex = 'a' ;
    
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;
  
  int ppv1a[4] = {0,1,0,0} ;  /* p : vertex */
  int ppv1c[4] = {-1,-1,0,0} ;  /* k + p : vertex */
  int ppv1f[4] = {1,0,0,0} ;  /* k : vertex */
  int ppv2d[4] = {1,1,0,0} ;  /* p : vertex */
  int ppv2b[4] = {0,-1,0,0} ;  /* k + p : vertex */
  int ppv2e[4] = {1,0,0,0} ;  /* k : vertex */

	
  POLYNOME p10 = newScalar (.5,h) ; /* real loop, symmetric */
  POLYNOME p1 = vertex_A_A_A (a,c,f,ppv1a,ppv1c,ppv1f,h) ;
  POLYNOME p2 = prop_AA (c,d,1,h) ;   /* (1/(pk)^2 */
  POLYNOME p3 = vertex_A_A_A (d,b,e,ppv2d,ppv2b,ppv2e,h) ;
  POLYNOME p4 = prop_AA (e,f,0,h) ;   /* (1/(k) */
  POLYNOME ppp[] = {p10,p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(polMultiProduct (h, ppp)) ;

  printf ("%s\n", title) ;
  showPol (pp) ;

  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  squareMomentaCleanUp (pp) ;
  showPol (pp) ;
  printf ("### Z2 AA loop A expect ::  3/2 (p_ab - g_ab p^2) : =3/2(vector) -1/6(ghost) = -5/3 (t'Hooft)\n") ;
  showPol (pp) ;

  printf ("DONE %s\n", title) ;

  return pp ;
} /* Z2_AA__loopA */

/***********************************************************************************************************************************************/

static POLYNOME Z2_AA__loopGhost  (const char *title) 
{
  AC_HANDLE h = ac_new_handle () ;
  firstDummyIndex = 'a' ;
    
  short mu = newDummyIndex () ;
  short nu = newDummyIndex () ;
  int ppv1[4] = {-1,0,0,0} ; /* k  : vertex */
  int ppv2[4] = {-1,-1,0,0} ; /* k + p : vertex */

  POLYNOME p10 = newScalar (-1,h) ; /* Ghost loop */
  POLYNOME p1 = vertex_A_c_cB (mu, ppv2,h) ; /*(2k + p)_mu */
  POLYNOME p2 = prop_cB_c (1,h) ;   /* (1/(k+p)^2 */
  POLYNOME p3 = vertex_A_c_cB (nu, ppv1,h) ; 
  POLYNOME p4 = prop_cB_c (0,h) ;   /* (1/(k)^2 */
  POLYNOME ppp[] = {p10,p1,p2,p3,p4,0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  pp = contractIndices(pp) ;

  printf ("%s\n", title) ;
  showPol (pp) ;

  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  pp = squareMomentaCleanUp (pp) ;
  showPol (pp) ;
  printf ("### Z2 AA loop ghost expect ::  1/6 p_ab + 1/12g_ab p^2\n") ;
  showPol (pp) ;

  printf ("DONE %s\n", title) ;

  return pp ;
} /* Z2_AA__loopGhost */

/***********************************************************************************************************************************************/

static POLYNOME Z2_BB__loopPsi  (const char *title) 
{
  AC_HANDLE h = ac_new_handle () ;
  firstDummyIndex = 'a' ;
    
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;

  POLYNOME p10 = newScalar (-1,h) ; /* Ghost loop */
  POLYNOME p1 = vertex_BB_PsiL_PsiRB (a,b,h) ;
  POLYNOME p2 = prop_PsiLB_PsiL (1,h) ;   /* (1/(k)^2 */
  POLYNOME p3 = vertex_B_PsiR_PsiLB (c,d,h) ;
  POLYNOME p4 = prop_PsiRB_PsiR (0,h) ;   /* (1/(k)^2 */
  POLYNOME ppp[] = {p10,p1,p2,p3,p4,0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("%s\n", title) ;
  showPol (pp) ;

  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;
  
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  pp = pauliTrace (pp, 0) ;
  showPol (pp) ;
  pp = squareMomentaCleanUp (pp) ;
  pp = reduceIndices (pp) ;
  pp = expand (pp) ;
  printf ("### Z2 Tensor avec loop PsiB_L Psi_L expect ::  je_sais_pas , on trouve 4/3 \n") ;
  showPol (pp) ;
  /* bbCleanUp is very fragile, i do not understand the order,
   * but this seems to work, we know the final symmetry in abcd 
   */
  pp = bbCleanUp (pp, b,a, d, c) ;
  showPol (pp) ;

  printf ("\n### this is the current definition of the propagator squezed between projectors expect 16\n") ;
  POLYNOME p5 = prop_BB_B (a, b, c, d, 0,h) ;
  showPol (p5) ;
  p5 = expand (p5) ;
  showPol (p5) ;
  p5 = bbCleanUp (p5, a, b, c, d) ;
  showPol (p5) ;
  printf ("DONE %s\n", title) ;

  
  return pp ;
} /* Z2_BB__loopPsi */

/***********************************************************************************************************************************************/

static POLYNOME Z2_HH__loopPsi  (const char *title) 
{
  AC_HANDLE h = ac_new_handle () ;
  firstDummyIndex = 'a' ;

  POLYNOME p10 = newScalar (-1,h) ; /* Ghost loop */
  POLYNOME p1 = vertex_HB_PsiL_PsiRB (h) ;
  POLYNOME p2 = prop_PsiLB_PsiL (1,h) ;   /* (1/(k)^2 */
  POLYNOME p3 = vertex_H_PsiR_PsiLB (h) ;
  POLYNOME p4 = prop_PsiRB_PsiR (0,h) ;   /* (1/(k)^2 */
  POLYNOME p5 = prop_HB_H (0,h) ;
  POLYNOME ppp[] = {p10,p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(polMultiProduct (h,ppp)) ;

  printf ("%s\n", title) ;
  showPol (p5) ;

  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  pp = pauliTrace (pp, 0) ;
  showPol (pp) ;
  pp = squareMomentaCleanUp (pp) ;
  pp = reduceIndices (pp) ;
  pp = expand (pp) ;
  printf ("### Z2 scalar avec loop PsiB_L Psi_L expect ::  je_sais_pas \n") ;
  showPol (pp) ;
  printf ("DONE %s\n", title) ;

  return pp ;
} /* Z2_HH__loopPsi */

/***********************************************************************************************************************************************/

static POLYNOME Z2_AA__loopPsi  (const char *title) 
{
  AC_HANDLE h = ac_new_handle () ;
  firstDummyIndex = 'a' ;
    
  short mu = newDummyIndex () ;
  short nu = newDummyIndex () ;

  POLYNOME p10 = newScalar (-1,h) ; /* Ghost loop */
  POLYNOME p1 = vertex_A_PsiL_PsiLB (mu,h) ;
  POLYNOME p2 = prop_PsiLB_PsiL (0,h) ;   /* (1/(k)^2 */
  POLYNOME p3 = vertex_A_PsiL_PsiLB (nu,h) ;
  POLYNOME p4 = prop_PsiLB_PsiL (1,h) ;   /* (1/(k+p)^2 */
  POLYNOME p5 = prop_AA (mu, nu, 0,h) ;
  POLYNOME ppp[] = {p10,p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(polMultiProduct (h,ppp)) ;

  printf ("%s\n", title) ;
  p5 = expand (p5) ;
  p5 = expand (p5) ;
  showPol (p5) ;
  showPol (pp) ;

  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  pp = dimIntegral (pp) ;
  showPol (pp) ;
  pp = pauliTrace (pp, 0) ;
  showPol (pp) ;
  pp = squareMomentaCleanUp (pp) ;
  showPol (pp) ;
  pp = reduceIndices (pp) ;
  pp = expand (pp) ;
  printf ("### Z2 Photon avec loop PsiB_L Psi_L expect :: 2/3 (p_ab - g_ab p^2) \n") ;
  showPol (pp) ;
  printf ("DONE %s\n", title) ;

  return pp ;
} /* Z2_AA__loopPsi */

/***********************************************************************************************************************************************/

static POLYNOME Z2_AA__loopB (const char *title) 
{
  AC_HANDLE h = ac_new_handle () ;
  firstDummyIndex = 'a' ;

  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;
  short g = newDummyIndex () ;
  short h1 = newDummyIndex () ;
  short i = newDummyIndex () ;
  short j = newDummyIndex () ;

  int mm1[4] = {1,1,0,0} ;    /* k+p */
  int mm2[4] = {1,0,0,0} ;    /* k */

  POLYNOME p0 = newScalar (1,h) ;
  POLYNOME p1 = vertex_A_B_BB (a,c,d,i,j,mm1,mm2,h) ;
  POLYNOME p2 = prop_BB_B (c,d,e,f,1,h) ;   /* (1/(k+p)^2 */
  POLYNOME p3 = vertex_A_B_BB (b,g,h1,e,f,mm2,mm1,h) ;
  POLYNOME p4 = prop_BB_B (g,h1,i,j,0,h) ;   /* (1/(k)^2 */
  POLYNOME ppp[] = {p0, p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(polMultiProduct (h, ppp)) ;

  printf ("%s\n", title) ;
  showPol (pp) ;
  freeIndex (pp) ;
  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  if (0)
    {
      
      printf ("freeIndex\n") ;
      freeIndex (pp) ;
      showPol (pp) ;
    }
  
  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  pp = squareMomentaCleanUp (pp) ;
  showPol(pp) ;
  pp = expand (pp) ;
  pp = expand (pp) ;
  
  printf ("### Z2 Photon avec loop BB_B expect :: x (p_ab - g_ab p^2, we get 2p_ab - 4 g_ab, this has to be wrong \n") ;
  showPol(pp) ;
  printf ("DONE %s\n\n", title) ;

  return pp ;
} /* Z2_AA__loopB */

/***********************************************************************************************************************************************/

static POLYNOME Z2_AA__loopHB (const char *title) 
{    
  AC_HANDLE h = ac_new_handle () ;
  firstDummyIndex = 'a' ;

  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;

  int mm1[4] = {0,1,0,0} ; /*  p  */
  int mm2[4] = {0,-1,0,0} ; /* -p  */

  POLYNOME p1 = prop_HB_H (1,h) ;   /* (1/(k+p)^2 */
  POLYNOME p2 = vertex_A_B_HB (a, c, d, mm1,h) ; 
  POLYNOME p3 = prop_BB_B (c,d,e,f, 0,h) ;   /* (1/(k)^2 */
  POLYNOME p4 = vertex_A_H_BB (b, e, f, mm2,h) ; /*(2k + p)_mu */

  POLYNOME ppp[] = {p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(polMultiProduct (h,ppp)) ;

  printf ("%s\n", title) ;
  showPol (pp) ;

  if (1)
    {
      printf ("Expand\n") ;
      pp = expand (pp) ;
      showPol (pp) ;
    }
  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  pp = squareMomentaCleanUp (pp) ;
  showPol(pp) ;
  pp = expand (pp) ;
  printf ("### Z2 Photon avec loop H_B expect zero (convergent because the derivatives on A  do not affect k)\n") ;
  showPol(pp) ;
  printf ("DONE %s\n\n", title) ;

  return pp ;
} /* Z2_AA__loopHB */

/***********************************************************************************************************************************************/

static POLYNOME Z2_HH__Aunder (const char *title) 
{
  AC_HANDLE h = ac_new_handle () ;
  firstDummyIndex = 'a' ;
    
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  int ppv[4] = {1,2,0,0} ; /* 2p + k : vertex */
  
  POLYNOME p1 = vertex_A_H_HB (b, ppv,h) ; /*(2k + p)_mu */
  POLYNOME p2 = prop_HB_H (1,h) ;   /* (1/(p+k)^2 */
POLYNOME p3 = vertex_A_H_HB (a, ppv,h) ;
  POLYNOME p4 = prop_AA (a,b,0,h) ;   /* (1/(k)^2 */
  POLYNOME ppp[] = {p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(polMultiProduct (h,ppp)) ;

  printf ("%s\n", title) ;
  showPol (pp) ;

  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  pp = squareMomentaCleanUp (pp) ;
  printf ("### Z2 scalar avec vector under, expect x p^2 ?, ZERO IN SU(1/1)\n") ;
  showPol (pp) ;
  printf ("DONE %s\n", title) ;

  return pp ;
} /* Z2_HH__Aunder */

/***********************************************************************************************************************************************/

static POLYNOME Z2_ghost__Aunder (const char *title) 
{
  AC_HANDLE h = ac_new_handle () ;
  firstDummyIndex = 'a' ;
    
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  int ppva[4] = {-1,-1,0,0} ; /* -p - k : vertex */
  int ppvb[4] = {0,1,0,0} ;   /* p : vertex */
  
  POLYNOME p1 = vertex_A_c_cB (b, ppva,h) ; /*(2k + p)_mu */
  POLYNOME p2 = prop_cB_c (1,h) ;   /* (1/(p+k)^2 */
  POLYNOME p3 = vertex_A_c_cB (a, ppvb,h) ; 
  POLYNOME p4 = prop_AA (a,b,0,h) ;   /* (1/(k)^2 */
  POLYNOME ppp[] = {p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(polMultiProduct (h,ppp)) ;

  printf ("%s\n", title) ;
  showPol (pp) ;

  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  pp = squareMomentaCleanUp (pp) ;
  printf ("### Z2 ghost avec vector under, expect (-1/2) p^2 ?, ZERO IN SU(1/1)\n") ;
  showPol (pp) ;
  printf ("DONE %s\n", title) ;

  return pp ;
} /* Z2_ghost__Aunder */

/***********************************************************************************************************************************************/

static POLYNOME Z2_HH__loopAB (const char *title) 
{
  AC_HANDLE h = ac_new_handle () ;
  firstDummyIndex = 'a' ;
    
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;
  int m1[4] = {-1, -1,0,0} ; /* -p - k : incoming A momentum */
  int m2[4] = {1,1,0,0} ; /* p + k : vertex */
  
  POLYNOME p1 = vertex_A_H_BB (a,c,d,m1,h) ; /*(2k + p)_mu */
  POLYNOME p2 = prop_BB_B (c,d,e,f,0,h) ;   /* (1/(p+k)^2 */
  POLYNOME p3 = vertex_A_B_HB (b,e,f, m2,h) ; 
  POLYNOME p4 = prop_AA (b,a,1,h) ;   /* (1/(k)^2 */
  POLYNOME ppp[] = {p1,p2,p3,p4,0} ; 
  /* POLYNOME ppp[] = {p2,p3,0} ; */

  POLYNOME pp = contractIndices(polMultiProduct (h,ppp)) ;

  printf ("%s\n", title) ;
  showPol (pp) ;
  freeIndex (pp) ;
  
  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  pp = squareMomentaCleanUp (pp) ;
  printf ("### Z2 scalar avec new ABH loop, expect x p^2 ?\n") ;
  showPol (pp) ;
  printf ("DONE %s\n", title) ;

  return pp ;
} /* Z2_HH__loopAB */

/***********************************************************************************************************************************************/

static POLYNOME Z2_BB__loopAH (const char *title) 
{
  AC_HANDLE h = ac_new_handle () ;
  firstDummyIndex = 'a' ;
    
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;
  short g = newDummyIndex () ;
  short h1 = newDummyIndex () ;
  short i = newDummyIndex () ;
  short j = newDummyIndex () ;
  int m1[4] = {1,0,0,0} ; /*  k : vertex */
  int m2[4] = {-1,0,0,0} ; /* -k : vertex */


  
  /*  POLYNOME p0 = newAG(a,b,g,h1,1,h) ; */
  POLYNOME p1 = vertex_A_H_BB (e,a,b,m1,h) ; /* k incoming */
  POLYNOME p2 = prop_HB_H (1,h) ;   /* (1/(p+k)^2 */
  POLYNOME p3 = vertex_A_B_HB (f,c,d, m2,h) ; /* k outgoing */
  POLYNOME p4 = prop_AA (f,e,0,h) ;   /* (1/(k)^2 */
  /*   POLYNOME p5 = newAG(c,d,i,j,-1,h) ; */
  POLYNOME ppp[] = {p1,p2,p3,p4,0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;

  printf ("%s\n", title) ;
  showPol (pp) ;
  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  pp = squareMomentaCleanUp (pp) ;
  showPol (pp) ;

  pp = expand (pp) ;
  printf ("### Z2 tensor avec new scalar-vector loop, before double projector \n") ;
  showPol (pp) ;

  exit (0) ;
  
  POLYNOME p02 = newAG(a,b,g,h1,1,h) ;
  POLYNOME p52 = newAG(i,j,c,d,-1,h) ;
  POLYNOME ppp2[] = {p02,pp,p52,0} ;
  pp = polMultiProduct (h,ppp2) ;
  showPol (pp) ;
  pp = expand (pp) ;
  showPol (pp) ;
  pp = reduceIndices(pp) ;
  showPol (pp) ;
  pp = expand (pp) ;
  showPol (pp) ;
  pp = expand (pp) ;

  printf ("### Z2 tensor avec new scalar-vector loop, expect je sais pas \n") ;
  showPol (pp) ;
  exit (0) ;

  return pp ;
} /* Z2_BB__loopAH */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME Z2_PsiL__B_Psi (const char *title) 
{
  AC_HANDLE h = ac_new_handle () ;
  short mu = newDummyIndex () ;
  short nu = newDummyIndex () ;
  short rho = newDummyIndex () ;
  short sigma = newDummyIndex () ;

  short w = newDummyIndex () ;
    
  POLYNOME p4 = vertex_BB_PsiL_PsiRB (mu,nu,h) ;
  POLYNOME p2 = prop_BB_B (mu,nu,rho,sigma, 0,h) ;   /* (1/(k)^2 */
  POLYNOME p3 = prop_PsiRB_PsiR (1,h) ; /* (1/(k+p)^2 */
  POLYNOME p1 = vertex_B_PsiR_PsiLB (rho,sigma,h) ;
  POLYNOME ppp[] = {p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(polMultiProduct (h,ppp)) ;
  
  printf ("%s\n",title) ;
  showPol (pp) ;

  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  printf ("pauli cleanup\n") ;
  printf ("### Z2 Psi left avec B_mu_nu under, expect ::  je_sais_pas * p-slash\n") ;
  printf ("...... Pauli cleanUp \n") ;
  pp = pauliCleanUp (pp, w) ;
  showPol(pp) ;

  printf ("Z2_  done\n\n") ;

  return pp ;
}

/***********************************************************************************************************************************************/

static POLYNOME Z2_PsiL__A_Psi (const char *title)
{
  AC_HANDLE h = ac_new_handle () ;
  short mu = newDummyIndex () ;
  short nu = newDummyIndex () ;
  short w = newDummyIndex () ;

  POLYNOME p1 = vertex_A_PsiL_PsiLB (mu,h) ;
  POLYNOME p2 = prop_AA (mu,nu, 0,h) ;   /* (1/(k)^2 */
  POLYNOME p3 = prop_PsiLB_PsiL (1,h) ; /* (1/(k+p)^2 */
  POLYNOME p4 = vertex_A_PsiL_PsiLB (nu,h) ;
  POLYNOME ppp[] = {p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(polMultiProduct (h,ppp)) ;
  
  printf ("%s\n",title) ;
  showPol (pp) ;

  pp = expand (pp) ;
  if (1) showPol (pp) ;

  pp = dimIntegral (pp) ;
  printf ("### Z2 Psi left avec A_mu under expect  p-slash\n") ;
  showPol (pp) ;
  pp = pauliCleanUp (pp, w) ;
  pp = contractIndices (pp) ;
  showPol (pp) ;
  return pp ;
}

/***********************************************************************************************************************************************/

static POLYNOME Z2_PsiL__H_Psi (const char *title)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL debug = FALSE ;
  POLYNOME p1 = vertex_HB_PsiL_PsiRB (h) ;
  POLYNOME p2 = prop_HB_H (0,h) ;          /* (1/(k)^2 */ 
  POLYNOME p3 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p4 = vertex_H_PsiR_PsiLB (h) ;
  POLYNOME ppp[] = {p1,p2,p3,p4,0} ;
  short w = newDummyIndex () ;
  
  POLYNOME pp = contractIndices(polMultiProduct (h,ppp)) ;
  
  printf ("%s\n", title) ;
  showPol (pp) ;
  
  expand (pp) ;

  if (debug) showPol (pp) ;
  pp = dimIntegral (pp) ;
  printf ("### Z2 Psi left avec H under, expect  1/2 p-slash\n") ;  
  showPol (pp) ;
  pp = pauliCleanUp (pp, w) ;
  showPol (pp) ;

  return pp ;
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/* check some integrals */
static void scpy (short *cp, const short *cq)
{
  while ((*cp++ = *cq++)) ;
  return ; 
}

static void sncpy (short *cp, const short *cq, int n)
{
  while (n-- > 0)
    *cp++ = *cq++ ;
  return ; 
}

static void tcpy (short *cp, const char *cq)
{
  while ((*cp++ = *cq++)) ; 
  return ; 
}

static void Thooft (void) 
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;
  short g = newDummyIndex () ;
  short h1 = newDummyIndex () ;

  POLYNOME pp, qq, rr ;

  pp = newPolynome (h) ;
  pp->tt.z = 1 ;
  pp->tt.type = 1 ;
  pp->tt.denom[0] = 1 ;
  pp->tt.denom[1] = 1 ;
  printf ("### Int { 1/k^2 (k+p)^2, expect I1\n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;

  pp = newPolynome (h) ;
  pp->tt.z = 1 ;
  pp->tt.type = 1 ;
  qq = newPolynome (h) ;
  qq->tt.z = 1 ;
  qq->tt.type = 1 ;

  pp->tt.mm[0][0] = a ;
  pp->tt.mm[0][1] = b ;
  scpy (qq->tt.g, pp->tt.mm[0]) ;
  pp->tt.denom[0] = 1 ;
  pp->tt.denom[1] = 1 ;
  printf ("### Int { k^a k^b/k^2 (k+p)^2, expect p_ab/3 - g_ab p^2/12\n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  printf ("# trace it, expect zero\n") ;
  rr = polProduct (pp, qq, h) ;
  rr = expand (rr) ;
  rr = reduceIndices (rr) ;
  rr = expand (rr) ;
  showPol (rr) ;

  pp = newPolynome (h) ;
  pp->tt.z = 1 ;
  pp->tt.type = 1 ;
  qq = newPolynome (h) ;
  qq->tt.z = 1 ;
  qq->tt.type = 1 ;

  pp->tt.mm[0][0] = a ;
  pp->tt.mm[0][1] = b ;
  scpy (qq->tt.g, pp->tt.mm[0]) ;
  pp->tt.denom[0] = 1 ;
  pp->tt.denom[1] = 1 ;
  pp->tt.denom[2] = 1 ;
  printf ("### Int { k^a k^b/k^2 (k+p)^2(k+p+q)^2, expect g_ab/4\n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  rr = polProduct (pp, qq,h) ;
  rr = expand (rr) ;
  showPol (rr) ;


  pp = newPolynome (h) ;
  pp->tt.z = 1 ;
  pp->tt.type = 1 ;
  qq = newPolynome (h) ;
  qq->tt.z = 1 ;
  qq->tt.type = 1 ;

  pp->tt.mm[0][0] = a ;
  pp->tt.mm[0][1] = b ;
  scpy (qq->tt.g, pp->tt.mm[0]) ;
  pp->tt.denom[0] = 1 ;
  pp->tt.denom[1] = 2 ;
  printf ("### Int { k^a k^b/k^2 (k+p)^4, expect g_ab/4\n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  rr = polProduct (pp, qq,h) ;
  rr = expand (rr) ;
  showPol (rr) ;


  pp = newPolynome (h) ;
  pp->tt.z = 1 ;
  pp->tt.type = 1 ;
  qq = newPolynome (h) ;
  qq->tt.z = 1 ;
  qq->tt.type = 1 ;

  pp->tt.mm[0][0] = a ;
  pp->tt.mm[0][1] = b ;
  pp->tt.mm[0][2] = c ;
  pp->tt.mm[0][3] = d ;
  scpy (qq->tt.g, pp->tt.mm[0]) ;
  pp->tt.denom[0] = 2 ;
  pp->tt.denom[1] = 2 ;
  printf ("### Int { k^a k^b k^c k^d/k^4 (k+p)^4, expect g_ab g_cd/24\n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  rr = polProduct (pp, qq, h) ;
  rr = expand (rr) ;
  rr = expand (rr) ;
  showPol (rr) ;

  pp = newPolynome (h) ;
  pp->tt.z = 1 ;
  pp->tt.type = 1 ;
  qq = newPolynome (h) ;
  qq->tt.z = 1 ;
  qq->tt.type = 1 ;

  pp->tt.mm[0][0] = a ;
  pp->tt.mm[0][1] = b ;
  pp->tt.mm[0][2] = c ;
  pp->tt.mm[0][3] = d ;
  pp->tt.mm[0][4] = e ;
  pp->tt.mm[0][5] = f ;
  scpy (qq->tt.g, pp->tt.mm[0]) ;
  pp->tt.denom[0] = 2 ;
  pp->tt.denom[1] = 3 ;
  printf ("### Int { k^a k^b k^c k^d k^e k^f/k^4 (k+p)^6, expect g_ab g_cd/192\n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  rr = polProduct (pp, qq, h) ;
  showPol (rr) ;
  rr = expand (rr) ;
  rr = expand (rr) ;
  showPol (rr) ;

  pp = newPolynome (h) ;
  pp->tt.z = 1 ;
  pp->tt.type = 1 ;
  qq = newPolynome (h) ;
  qq->tt.z = 1 ;
  qq->tt.type = 1 ;

  pp->tt.mm[0][0] = a ;
  pp->tt.mm[0][1] = b ;
  pp->tt.mm[0][2] = c ;
  pp->tt.mm[0][3] = d ;
  pp->tt.mm[0][4] = e ;
  pp->tt.mm[0][5] = f ;
  pp->tt.mm[0][6] = g ;
  pp->tt.mm[0][7] = h1 ;
  scpy (qq->tt.g, pp->tt.mm[0]) ;
  pp->tt.denom[0] = 3 ;
  pp->tt.denom[1] = 3 ;
  pp->tt.z = 1920 ;
  printf ("### Int { k^a k^b k^c k^d k^e k^f k^g k^h/k^6 (k+p)^6, expect g_ab g_cd (i multiplied by 1920\n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  rr = polProduct (pp, qq, h) ;
  rr = expand (rr) ;
  rr = expand (rr) ;
  showPol (rr) ;
 

 pp = newPolynome (h) ;
  pp->tt.z = 1 ;
  pp->tt.type = 1 ;
  qq = newPolynome (h) ;
  qq->tt.z = 1 ;
  qq->tt.type = 1 ;

  pp->tt.mm[0][0] = a ;
  pp->tt.mm[0][1] = b ;
  pp->tt.mm[0][2] = c ;
  pp->tt.mm[0][3] = d ;
  pp->tt.mm[0][4] = e ;
  sncpy (qq->tt.g, pp->tt.mm[0],4) ;
  pp->tt.denom[0] = 1 ;
  pp->tt.denom[1] = 1 ;
  pp->tt.denom[2] = 1 ;
  pp->tt.denom[3] = 1 ;
  printf ("### Int { k^a k^b k^c k^d k^e /k^2 (k+p)^2 (k+p+q)^2 (k+p+q+r)^2, expect p_a g_bc g_de/ ?\n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  rr = polProduct (pp, qq, h) ;
  showPol (rr) ;
  rr = expand (rr) ;
  showPol (rr) ;
  rr = expand (rr) ;
  showPol (rr) ;

  ac_free (h) ;
  return ;
} /* tHooft */

/***********************************************************************************************************************************************/
/* check self duality projectors */
static POLYNOME Hodge (void) 
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;
  short g = newDummyIndex () ;
  short h1 = newDummyIndex () ;
  short i = newDummyIndex () ;
  short j = newDummyIndex () ;
  short mu = newDummyIndex () ;
  short nu = newDummyIndex () ;
  short rho = newDummyIndex () ;
  short sig = newDummyIndex () ;

  POLYNOME pp, ppp[6], pP1, pP2, pM2,pG1, pG2 ;
  /*
    int mm1[4] = {1,1,0,0} ; //
    int mm2[4] = {-1,-1,0,0} ; //  p + k 
  */

  if (1)
    {
      pP1 = newAG (mu,nu,a,b,1,h) ;
      /*  pM1 = newAG (mu,nu,a,b,-1, h) ; */
      pP2 = newAG (rho,sig,a,b,1, h) ;
      pM2 = newAG (rho,sig,a,b,-1, h) ;
      pG1 = newAG (mu,nu,a,b,0, h) ;
      pG2 = newAG (rho,sig,a,b,0, h) ;
    }
  if (0)
    {
      pp = newG (mu,nu,h) ;
      tcpy (pp->tt.eps, "abcdabcd") ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      
      pp = newG (mu,nu,h) ;
      tcpy (pp->tt.eps, "abcdabdc") ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      
      pp = newG (mu,nu,h) ;
      tcpy (pp->tt.eps, "abcdacbd") ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      
      pp = newG (mu,nu,h) ;
      tcpy (pp->tt.eps, "abcdabce") ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      
      pp = newG (mu,nu,h) ;
      tcpy (pp->tt.eps, "abcdaebc") ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      
      
      pp = newG (mu,nu,h) ;
      tcpy (pp->tt.eps, "abcdabec") ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
    }
  
  if (0)
    {
      ppp[0] = newG ('a','b',h) ;
      ppp[1] = 0 ;
      pp = polMultiSum (h,ppp) ;
      showPol (pp) ;
      pp = sortPol (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      ppp[0] = newG ('a','b',h) ;
      ppp[1] = newG ('c','d',h) ;
      ppp[2] = 0 ;
      pp = polMultiSum (h,ppp) ;
      showPol (pp) ;
      pp = sortPol (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      ppp[0] = newG ('a','b',h) ;
      ppp[1] = newG ('c','d',h) ;
      ppp[2] = newG ('a','b',h) ;
      ppp[3] = 0 ;
      pp = polMultiSum (h,ppp) ;
      showPol (pp) ;
      pp = sortPol (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

    }
  if (1)
    {
      ppp[0] = newG ('a','b',h) ;
      ppp[1] = newG ('c','d',h) ;
      ppp[2] = newG ('a','b',h) ;
      ppp[3] = newG ('c','d',h) ;
      ppp[4] = 0 ;
      pp = polMultiSum (h,ppp) ;
      showPol (pp) ;
      pp = sortPol (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      exit (0) ;

      ppp[0] = newG ('a','b',h) ;
      ppp[1] = newG ('c','d',h) ;
      ppp[2] = newG ('a','b',h) ;
      ppp[3] = newG ('c','d',h) ;
      ppp[4] = newG ('a','b',h) ;
      ppp[5] = 0 ;
      pp = polMultiSum (h,ppp) ;
      showPol (pp) ;
      pp = sortPol (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
    }

  if (1)
    {
      printf ("Is PP a projector: compute PP^2\n") ;
      ppp[0] = pP1 ;
      ppp[1] = pP2 ;
      ppp[2] = 0 ;
      pp = polMultiProduct (h,ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;

      showPol (pp) ;
  
      printf ("Is PG a projector: compute PG PG\n") ;
      ppp[0] = pG1 ;
      ppp[1] = pG2 ;
      ppp[2] = 0 ;
      pp = polMultiProduct (h,ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      pp = contractProducts (pp) ;
      contractIndices (pp) ;
      pp = expand (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      
      printf ("Is PP * PM zero: compute PP PM\n") ;
      ppp[0] = pP1 ;
      ppp[1] = pM2 ;
      ppp[2] = 0 ;
      pp = polMultiProduct (h,ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      pp = contractProducts (pp) ;
      contractIndices (pp) ;
      pp = expand (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      
      printf ("Test : compute PP PG\n") ;
      ppp[0] = pP1 ;
      ppp[1] = pG2 ;
      ppp[2] = 0 ;
      pp = polMultiProduct (h,ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      
      printf ("Test : compute PG eps\n") ;
      ppp[0] = pG1 ;
      ppp[1] = newEpsilon (rho,sig,a,b,h) ;
      ppp[2] = 0 ;
      pp = polMultiProduct (h,ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      printf ("Test : compute B propagator\n") ;
      ppp[0] = newAG (a,b,e,f,-1,h) ; 
ppp[1] = prop_BB_B (e,f,g,h1,0,h) ;
      ppp[2] = newAG (g,h1,c,d,1,h) ; 
      ppp[3] = 0 ;
      pp = polMultiProduct (h,ppp) ;
      showPol (pp) ;
      pp = reduceIndices (pp) ;
      pp = expand (pp) ;
      pp = reduceIndices (pp) ;
      pp = expand (pp) ;
      pp = reduceIndices (pp) ;
      pp = expand (pp) ;
      pp = sortPol (pp) ;
      pp = expand (pp) ;
      pp = expand (pp) ;
      
      showPol (pp) ;
    }

  if (1)
    {
      printf ("Test : compute eps eps\n") ;
      ppp[0] = newEpsilon (mu,nu,a,b,h) ;
      ppp[1] = newEpsilon (rho,b,a,sig,h) ;
      ppp[0]->tt.z = I/4 ;
      ppp[1]->tt.z = I/4 ;
      ppp[2] = 0 ;
      pp = polMultiProduct (h,ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;

      showPol (pp) ;
    }
  
  if (1)
    {
      printf ("Test : compute chiral projectors PP PP = PP\n") ;
      ppp[0] = newAG (a,b,c,d,1,h) ;
      ppp[1] = newAG (c,d,e,f,1,h) ;
      ppp[2] = 0 ;
      pp = polMultiProduct (h,ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      printf ("Test : compute chiral projectors PP PP PP = PP\n") ;
      ppp[0] = newAG (a,b,c,d,1,h) ;
      ppp[1] = newAG (c,d,e,f,1,h) ;
      ppp[2] = newAG (e,f,g,h1,1,h) ;
      ppp[3] = 0 ;
      pp = polMultiProduct (h,ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      printf ("Test : compute chiral projectors PP PP PP PP = PP\n") ;
      ppp[0] = newAG (a,b,c,d,1,h) ;
      ppp[1] = newAG (c,d,e,f,1,h) ;
      ppp[2] = newAG (e,f,g,h1,1,h) ;
      ppp[3] = newAG (g,h1,i,j,1,h) ;
      ppp[4] = 0 ;
      pp = polMultiProduct (h,ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      printf ("Test : compute chiral projectors PP PM = 0\n") ;
      ppp[0] = newAG (a,b,c,d,1,h) ;
      ppp[1] = newAG (c,d,e,f,-1,h) ;
      ppp[2] = 0 ;
      pp = polMultiProduct (h,ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      printf ("Test : compute chiral projectors PM PP = 0\n") ;
      ppp[0] = newAG (a,b,c,d,-1,h) ;
      ppp[1] = newAG (c,d,e,f,1,h) ;
      ppp[2] = 0 ;
      pp = polMultiProduct (h,ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      printf ("Test : compute chiral projectors PP PM PP = PP\n") ;
      ppp[0] = newAG (a,b,c,d,1,h) ;
      ppp[1] = newAG (c,d,e,f,-1,h) ;
      ppp[2] = newAG (e,f,g,h1,1,h) ;
      ppp[3] = 0 ;
      pp = polMultiProduct (h,ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      printf ("Test : double epsilon, no common index\n") ;
      pp = newG (a,b,h) ;
      tcpy (pp->tt.eps, "cdefghij") ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
    }
  ac_free (h) ;
  exit (0) ;
  return pp ;
} /* Hodge */

/***********************************************************************************************************************************************/

static POLYNOME Z2_BB__Aunder (const char *title) 
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;
  short g = newDummyIndex () ;
  short h1 = newDummyIndex () ;
  short i = newDummyIndex () ;
  short j = newDummyIndex () ;
  short k = newDummyIndex () ;
  short l = newDummyIndex () ;
  short mu = newDummyIndex () ;
  short nu = newDummyIndex () ;

  POLYNOME pp, ppp[12] ;
  int mm1[4] = {0,1,0,0} ;   /* p + k */
  int mm2[4] = {1,0,0,0} ; /* p + k */
  int mm3[4] = {0,1,0,0} ;   /* p */

  printf ("\n\n%s\n", title) ;

  if (1)
    { 
      ppp[0] = newAG (a,b,i,j, -1,h) ;
      ppp[0]->tt.z = 384 ;
      ppp[1] = prop_AA (mu, nu, 1,h) ; /* 1/(p+k^2 */
      ppp[2] = vertex_A_B_BB (mu, i, j, e, f, mm1, mm2,h) ;
      ppp[3] = prop_BB_B (e,f,g,h1,0,h) ; /* kk/(k)^4 */
      ppp[4] = vertex_A_B_BB (nu, g, h1, k, l, mm2, mm3,h) ;
      ppp[5] = newAG (k,l,c,d, 1,h) ;
      ppp[6] = newScalar (384,h) ;
      ppp[7] = 0 ;

      pp = polMultiProduct (h, ppp) ;
      printf ("########### 384 Z2 B with transient A TRES FAUX resultat en p^4 au lieu de p^2 \n") ;  
      showPol (pp) ;
      pp = expand (pp) ;
      pp = dimIntegral (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      pp = squareMomentaCleanUp (pp) ;
      showPol(pp) ;
      pp = reduceIndices (pp) ;
      pp = expand (pp) ;
      showPol(pp) ;
    }  
  return pp ;
} /* Z2_BB__Aunder */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
static POLYNOME Z3_AHH__loopPsiL (const char *title)
{
  AC_HANDLE h = ac_new_handle () ;
  short mu = newDummyIndex () ;
  /* calcul plus simple, on calcule successivement q=r=0, r=p=0, p=q=0,  et on additionne */

  /* set q == 0 */
  POLYNOME p10 = newScalar (-1,h) ; /* Fermion loop */
  POLYNOME p11 = prop_PsiLB_PsiL (0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p12 = vertex_A_PsiL_PsiLB (mu,h) ;
  POLYNOME p13 = prop_PsiLB_PsiL (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p14 = vertex_HB_PsiL_PsiRB (h) ;
  POLYNOME p15 = prop_PsiRB_PsiR (2,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p16 = vertex_H_PsiR_PsiLB (h) ;
  POLYNOME pppP[] = {p10,p11,p12,p13,p14,p15,p16,0} ;

  POLYNOME PP = contractIndices(polMultiProduct (h,pppP)) ;

  printf ("%s\n", title) ;
  printf ("############# Z3 A H HB with loop psiL (psiL touches A)\n") ;
  showPol (PP) ;

  PP = dimIntegral (PP) ;
  if (0) showPol (PP) ;
  printf ("# trace the pauli matrix :\n ") ;
  PP = pauliTrace (PP, 0) ;
  if (0) showPol (PP) ;
  PP = expand(PP) ;
  PP = reduceIndices (PP) ;
  PP = expand(PP) ;
  printf ("### Z3 A H HB with psi loop, expect (-p-2q\n") ;
  showPol (PP) ;

  return PP ;
} /* Z3_A_HHH__loopPsiL */

/***********************************************************************************************************************************************/

static POLYNOME Z3_ABH__loopPsiL (const char *title)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  /* calcul plus simple, on calcule successivement q=r=0, r=p=0, p=q=0,  et on additionne */

  /* set q == 0 */
  POLYNOME p10 = newScalar (-1,h) ; /* Fermion loop */
  POLYNOME p11 = prop_PsiLB_PsiL (0,h) ;   /* (1/(k)^2 */ 

  POLYNOME p12 = vertex_B_PsiR_PsiLB (b,c,h) ;
  POLYNOME p13 = prop_PsiRB_PsiR (2,h) ;   /* (1/(k+p+q)^2 */
  POLYNOME p14 = vertex_HB_PsiL_PsiRB (h) ;

  POLYNOME p15 = prop_PsiLB_PsiL (1,h) ;   /* (1/(k+p)^2 */
  POLYNOME p16 = vertex_A_PsiL_PsiLB (a,h) ;


  POLYNOME pppP[] = {p10,p11,p12,p13,p14,p15,p16,0} ;
  
POLYNOME PP = contractIndices(polMultiProduct (h,pppP)) ;

  printf ("\n\n%s\n", title) ;
  printf ("############# Z3 A B HB with psi loop R\n") ;
  showPol (PP) ;


  PP = dimIntegral (PP) ;
  if (0) showPol (PP) ;
  printf ("# trace the pauli matrix :\n ") ;
  PP = pauliTrace (PP, 0) ;
  if (0) showPol (PP) ;
  PP = expand(PP) ;
  PP = reduceIndices (PP) ;
  PP = expand(PP) ;
  printf ("### Z3 A B HB with psi loop L, expect (-p: we derive the incoming vector\n") ;
  showPol (PP) ;

  return PP ;
} /* Z3_ABH__loopPsiL */

/***********************************************************************************************************************************************/

static POLYNOME Z3_ABB__loopPsiL (const char *title)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  
  int mm1[4] = {0, 1, 0, 0} ;
  int mm2[4] = {0, 0, 1, 0} ;
  
  /* set q == 0 */
  POLYNOME p10 = newScalar (-1,h) ; /* Fermion loop */
  POLYNOME p11 = prop_PsiLB_PsiL (0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p12 = vertex_A_PsiL_PsiLB (a,h) ;
  POLYNOME p13 = prop_PsiLB_PsiL (2,h) ;   /* (1/(k+p+q)^2 */
  POLYNOME p14 = vertex_B_PsiR_PsiLB (b,c,h) ;
  POLYNOME p15 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p16 = vertex_BB_PsiL_PsiRB (d,e,h) ;
  POLYNOME p20 = vertex_A_B_BB (a, b,c,d,e, mm1, mm2,h) ;
  
  POLYNOME ppp[] = {p10, p11,p12,p13,p14,p15,p16,0} ;

  POLYNOME PP = polMultiProduct (h,ppp) ;

  p20 = expand (p20) ;
  
  printf ("\n\n%s\n", title) ;
  printf ("############# Z3 A B BB raw vertex\n") ;
  showPol (p20) ;
  printf ("############# Z3 A B BB with psi loop L\n") ;
  PP = expand(PP) ;
  PP = expand(PP) ;
  showPol (PP) ;

  PP = dimIntegral (PP) ;
  if (0) showPol (PP) ;
  PP = expand(PP) ;
  PP = reduceIndices (PP) ;
  PP = expand(PP) ;
  if (1) showPol (PP) ;
  PP = sortPol (PP) ;
  if (1) showPol (PP) ;
  printf ("# trace the pauli matrix :\n ") ;

  PP = pauliTrace (PP, 0) ;
  if (0) showPol (PP) ;
  PP = expand(PP) ;
  if (1) PP = reduceIndices (PP) ;
  if (1) showPol (PP) ;
  PP = expand(PP) ;
  PP = expand(PP) ;
  showPol (PP) ;

  exit (0) ;
  
  ppp[0] = newSigma (b,h) ;
  ppp[0]->tt.sigma[1] = c ;
  ppp[1] = PP ;
  ppp[2] = newSigma (d,h) ;
  ppp[2]->tt.sigB[1] = e ;
  ppp[3] = 0 ;
  
  printf ("### Z3 A B BB with psi loop L, expect 2p+q\n") ;
  PP = polMultiProduct (h,ppp) ;
  PP = expand(PP) ;
  showPol (PP) ;

  exit (0) ;
  return PP ;
} /* Z3_ABB__loopPsiL */

/***********************************************************************************************************************************************/

static POLYNOME Z3_ABB__loopPsiR (const char *title)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;

  POLYNOME p0 = newSigma (a,h) ;
  tcpy (p0->tt.sigma, "abcdefdh")  ;
  printf ("\n\n%s\n", title) ;
  showPol (p0) ;
  p0 = expand(p0) ;
  showPol (p0) ;


  /* set q == 0 */
  POLYNOME p10 = newScalar (-1,h) ; /* Fermion loop */
  POLYNOME p11 = prop_PsiRB_PsiR (0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p12 = vertex_BB_PsiL_PsiRB (d,e,h) ;
  POLYNOME p13 = prop_PsiLB_PsiL (2,h) ;   /* (1/(k+p+q)^2 */
  POLYNOME p14 = vertex_B_PsiR_PsiLB (b,c,h) ;
  POLYNOME p15 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k+p)^2 */
  POLYNOME p16 = vertex_A_PsiR_PsiRB (a,h) ;
  
  POLYNOME pppP[] = {p10, p11,p12,p13,p14,p15,p16,0} ;

  POLYNOME PP = contractIndices(polMultiProduct (h,pppP)) ;

  printf ("############# Z3 A B BB with right psi loop L\n") ;
  showPol (PP) ;

  PP = expand(PP) ;
  if (1) showPol (PP) ;
  PP = dimIntegral (PP) ;
  if (1) showPol (PP) ;
  PP = expand(PP) ;
  if (1) showPol (PP) ;
  PP = reduceIndices (PP) ;
  PP = expand(PP) ;
  if (1) showPol (PP) ;
  printf ("# trace the pauli matrix :\n ") ;
  PP = pauliTrace (PP, 0) ;
  if (0) showPol (PP) ;
  PP = expand(PP) ;
  PP = reduceIndices (PP) ;
  PP = expand(PP) ;
  printf ("### Z3 A B BB with psi loop R, expect 2p+q\n") ;
  showPol (PP) ;

  return PP ;
} /* Z3_ABB__loopPsiR */

/***********************************************************************************************************************************************/

static POLYNOME Z3_AAA__loopPsiL (const char *title)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short n = newDummyIndex () ;

  /* calcul plus simple, on calcule successivement q=r=0, r=p=0, p=q=0,  et on additionne */

  /* set q == 0 */
  POLYNOME p10 = newScalar (-1,h) ; /* Fermion loop */
  POLYNOME p11 = vertex_A_PsiL_PsiLB (a,h) ;
  POLYNOME p12 = prop_PsiLB_PsiL (0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p13 = vertex_A_PsiL_PsiLB (c,h) ;
  POLYNOME p14 = prop_PsiLB_PsiL (2,h) ;   /* (1/(k+p+q)^2 */
  POLYNOME p15 = vertex_A_PsiL_PsiLB (b,h) ;
  POLYNOME p16 = prop_PsiLB_PsiL (1,h) ;   /* (1/(k+p)^2 */ 



  POLYNOME pppP[] = {p10, p11,p12,p13,p14,p15,p16,0} ;

  POLYNOME PP = contractIndices(polMultiProduct (h,pppP)) ;

  printf ("\n\n%s\n", title) ;
  printf ("############# Z3 A A A with psi loop\n") ;
  showPol (PP) ;

  PP = dimIntegral (PP) ;
  showPol (PP) ;
  printf ("# trace the pauli matrix :\n ") ;
  PP = pauliTrace (PP, 0) ;
  showPol (PP) ;
  PP = momentaCleanUp (PP, n) ;
  PP = expand (PP) ;
  printf ("### Z3 A A A with psi loop, expect g_bc (p+2q)_rho + .. + ..\n") ;
  showPol (PP) ;


  return PP ;
} /* Z3_AAA__loopPsiL */

/***********************************************************************************************************************************************/

static POLYNOME Z3_AAA__loopH (const char *title)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short n = newDummyIndex () ;

  int pa[4] = {2,1,1,0} ;
  int pb[4] = {2,1,0,0} ;
  int pc[4] = {2,2,1,0} ;
  
  /* calcul plus simple, on calcule successivement q=r=0, r=p=0, p=q=0,  et on additionne */

  /* set q == 0 */
  POLYNOME p10 = newScalar (1,h) ; /* scalar loop */
  POLYNOME p11 = vertex_A_H_HB (a, pa,h) ;
  POLYNOME p12 = prop_HB_H (0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p13 = vertex_A_H_HB (c, pc,h) ;
  POLYNOME p14 = prop_HB_H (2,h) ;   /* (1/(k+p+q)^2 */
  POLYNOME p15 = vertex_A_H_HB (b, pb,h) ;
  POLYNOME p16 = prop_HB_H (1,h) ;   /* (1/(k+p)^2 */ 



  POLYNOME pppP[] = {p10, p11,p12,p13,p14,p15,p16,0} ;

  POLYNOME PP = contractIndices(polMultiProduct (h,pppP)) ;

  printf ("\n\n%s\n", title) ;
  printf ("############# Z3 A A A with complex scalar loop\n") ;
  showPol (PP) ;

  PP = dimIntegral (PP) ;
  showPol (PP) ;
  printf ("# trace the pauli matrix :\n ") ;
  PP = pauliTrace (PP, 0) ;
  showPol (PP) ;
  PP = momentaCleanUp (PP, n) ;
  PP = expand (PP) ;
  printf ("### Z3 A A A with complex scalar loop, expect g_bc (p-q)_rho + .. + ..\n") ;
  showPol (PP) ;


  return PP ;
} /* Z3_AAA__loopH */

/***********************************************************************************************************************************************/

static POLYNOME Z3_AAA__loopGhost (const char *title)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short n = newDummyIndex () ;

  /* calcul plus simple, on calcule successivement q=r=0, r=p=0, p=q=0,  et on additionne */

  /* set q == 0 */
  POLYNOME p10 = newScalar (-1,h) ; /* Ghost loop */
  POLYNOME p11 = vertex_A_PsiL_PsiLB (a,h) ;
  POLYNOME p12 = prop_PsiLB_PsiL (0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p13 = vertex_A_PsiL_PsiLB (c,h) ;
  POLYNOME p14 = prop_PsiLB_PsiL (2,h) ;   /* (1/(k+p+q)^2 */
  POLYNOME p15 = vertex_A_PsiL_PsiLB (b,h) ;
  POLYNOME p16 = prop_PsiLB_PsiL (1,h) ;   /* (1/(k+p)^2 */ 



  POLYNOME pppP[] = {p10, p11,p12,p13,p14,p15,p16,0} ;

  POLYNOME PP = contractIndices(polMultiProduct (h,pppP)) ;

  printf ("\n\n%s\n", title) ;
  printf ("############# Z3 A A A with ghost loop\n") ;
  showPol (PP) ;

  PP = dimIntegral (PP) ;
  showPol (PP) ;
  printf ("# trace the pauli matrix :\n ") ;
  PP = pauliTrace (PP, 0) ;
  showPol (PP) ;
  PP = momentaCleanUp (PP, n) ;
  PP = expand (PP) ;
  printf ("### Z3 A A A with ghost loop, expect g_bc (p+2q)_rho + .. + ..\n") ;
  showPol (PP) ;


  return PP ;
} /* Z3_AAA__loopGhost */

/***********************************************************************************************************************************************/

static POLYNOME Z3_AAA__loopA (const char *title)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;
  short g = newDummyIndex () ;
  short h1 = newDummyIndex () ;
  short i = newDummyIndex () ;


  int pa[4] = {0,-1,-1,0} ;
  int pb[4] = {0,1,0,0} ;
  int pc[4] = {0,0,1,0} ;
  int pd[4] = {-1,0,0,0} ;
  int pe[4] = {1,0,0,0} ;
  int pf[4] = {-1,-1,0,0} ;
  int pg[4] = {1,1,0,0} ;
  int ph[4] = {-1,-1,-1,0} ;
  int pi[4] = {1,1,1,0} ;

  /* calcul plus simple, on calcule successivement q=r=0, r=p=0, p=q=0,  et on additionne */

  /* set q == 0 */
  POLYNOME p10 = newScalar (1,h) ; /* Vector loop */
  POLYNOME p11 = vertex_A_A_A (a,d,i,pa,pd,pi,h) ;
  POLYNOME p12 = prop_AA (d,e,0,h) ;
  POLYNOME p13 = vertex_A_A_A (b,f,e,pb,pf,pe,h) ;
  POLYNOME p14 = prop_AA (f,g,1,h) ;
  POLYNOME p15 = vertex_A_A_A (c,h1,g,pc,ph,pg,h) ;
  POLYNOME p16 = prop_AA (h1,i,2,h) ;



  POLYNOME pppP[] = {p10, p11,p12,p13,p14,p15,p16,0} ;

  POLYNOME PP = contractIndices(polMultiProduct (h,pppP)) ;

  printf ("\n\n%s\n", title) ;
  printf ("############# Z3 A A A with vector loop\n") ;
  showPol (PP) ;

  PP = dimIntegral (PP) ;
  showPol (PP) ;
  printf ("# trace the pauli matrix :\n ") ;
  PP = pauliTrace (PP, 0) ;
  showPol (PP) ;
  PP = momentaCleanUp (PP, 'n') ;
  PP = expand (PP) ;
  printf ("### Z3 A A A with psi loop, expect g_bc (p+2q)_rho + .. + ..\n") ;
  showPol (PP) ;


  return PP ;
} /* Z3_AAA__loopA */

/***********************************************************************************************************************************************/

static POLYNOME Z3_AAA__loopB (const char *title)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short n = newDummyIndex () ;

  /* calcul plus simple, on calcule successivement q=r=0, r=p=0, p=q=0,  et on additionne */

  /* set q == 0 */
  POLYNOME p10 = newScalar (-1,h) ; /* Ghost loop */
  POLYNOME p11 = vertex_A_PsiL_PsiLB (a,h) ;
  POLYNOME p12 = prop_PsiLB_PsiL (0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p13 = vertex_A_PsiL_PsiLB (c,h) ;
  POLYNOME p14 = prop_PsiLB_PsiL (2,h) ;   /* (1/(k+p+q)^2 */
  POLYNOME p15 = vertex_A_PsiL_PsiLB (b,h) ;
  POLYNOME p16 = prop_PsiLB_PsiL (1,h) ;   /* (1/(k+p)^2 */ 



  POLYNOME pppP[] = {p10, p11,p12,p13,p14,p15,p16,0} ;

  POLYNOME PP = contractIndices(polMultiProduct (h,pppP)) ;

  printf ("\n\n%s\n", title) ;
  printf ("############# Z3 A A A with tensor loop\n") ;
  showPol (PP) ;

  PP = dimIntegral (PP) ;
  showPol (PP) ;
  printf ("# trace the pauli matrix :\n ") ;
  PP = pauliTrace (PP, 0) ;
  showPol (PP) ;
  PP = momentaCleanUp (PP, n) ;
  PP = expand (PP) ;
  printf ("### Z3 A A A with tensor loop, expect g_bc (p+2q)_rho + .. + ..\n") ;
  showPol (PP) ;


  return PP ;
} /* Z3_AAA__loopB */

/***********************************************************************************************************************************************/

static POLYNOME Z3_A_H_BB__loopABH (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;
  short g = newDummyIndex () ;
  short h1 = newDummyIndex () ;
  short i = newDummyIndex () ;

  /* set p=q == 0 in the vertex and propagators since the integral is k^4/k^8 */
#ifdef JUNK
  /* true momenmta  */
  int mmA1[4] = {0,1,0,0} ; /* p-incoming on photon lne at bottom left of triangle*/
  int mmA2[4] = {-1,-1,-1,0} ; /* p-incoming on photon lne at bottom right of triangle */
  int mmA3[4] = {1,1,1,0} ; /* p-incoming on photon lne at top of triangle */
#else
  /* simplified momenmta  */
  int mmA1[4] = {0,1,0,0} ; /* p-incoming on photon lne at bottom left of triangle*/
  int mmA2[4] = {-1, 0, 0, 0} ; /* p-incoming on photon lne at bottom right of triangle */
  int mmA3[4] = {1,0,0,0} ; /* p-incoming on photon lne at top of triangle */
#endif
    

  POLYNOME p1 = vertex_A_B_HB (a,h1,i,mmA1,h) ;  
  POLYNOME p2 = vertex_A_H_BB (d,f,g,mmA2,h) ;
  POLYNOME p3 = vertex_A_H_BB (e,b,c,mmA3,h) ;
  POLYNOME p12 = prop_BB_B (f,g,h1,i,0,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p23 = prop_AA (d,e,0,h) ;         /* (1/(k+p+q)^2 */
  POLYNOME p31 = prop_HB_H (0,h) ;           /* (1/(k)^2 */ 



  POLYNOME pppP[] = {p1,p12,p2,p23,p3,p31,0} ;

  POLYNOME PP = contractIndices(polMultiProduct (h,pppP)) ;

  showPol (PP) ;

  PP = dimIntegral (PP) ;
  showPol (PP) ;
  
  short n = newDummyIndex () ;
  PP = momentaCleanUp (PP, n) ;
  PP = expand (PP) ;

  printf ("### Z3  A_H_BB with psi boson loop, expect g_ac (p)_b + .. + ..\n") ;
  showPol (PP) ;


  return PP ;
} /* Z3_H_B_BB__loopABH */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME ZF_H_PsiR_PsiLB (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;

  int mm[4] = {1,0,0,0} ;

  POLYNOME p1 = vertex_A_PsiR_PsiRB (b,h) ;
  POLYNOME p2 = prop_PsiRB_PsiR (0,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p3 = vertex_B_PsiR_PsiLB (e,f,h) ;

  POLYNOME p4 = prop_BB_B (c,d,e,f, 0,h) ; /* 1/k^2 */
  POLYNOME p5 = prop_AA (a,b,0,h) ; /* 1/k^2 */
  POLYNOME p6 = vertex_A_B_HB (a,c,d,mm,h) ;
  POLYNOME ppp[7] = {p1,p2,p3,p4,p5,p6,0} ;
  POLYNOME pp = polMultiProduct (h,ppp) ;


  printf ("Z3 New vertex A/H/B contrib to vertex phi-psi-psi : \n") ;
  showPol(p6) ;
  showPol(pp) ;
  pp = expand (pp) ;
  showPol(pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;


  return pp ;
} /* ZF_H_PsiR_PsiLB */

/***********************************************************************************************************************************************/

static POLYNOME Z3_A_PsiL_PsiLB__Hunder (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;

  POLYNOME p1 = vertex_HB_PsiL_PsiRB (h) ;
  POLYNOME p2 = prop_PsiRB_PsiR (0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p3 = vertex_A_PsiR_PsiRB (a,h) ;
  POLYNOME p4 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p5 = vertex_H_PsiR_PsiLB (h) ;
  POLYNOME p6 = prop_HB_H (2,h) ; /* 1/(k+p+q)^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex A psi with H under: expect 1/2 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;


  return pp ;
} /* Z3_A_PsiL_PsiLB__Hunder */

/***********************************************************************************************************************************************/

static POLYNOME Z3_A_PsiL_PsiLB__Aunder (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;

  POLYNOME p1 = vertex_A_PsiL_PsiLB (b,h) ;
  POLYNOME p2 = prop_PsiLB_PsiL (0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p3 = vertex_A_PsiL_PsiLB (a,h) ;
  POLYNOME p4 = prop_PsiLB_PsiL (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p5 = vertex_A_PsiL_PsiLB (c,h) ;
  POLYNOME p6 = prop_AA (b,c,2,h) ; /* 1/(k+p+q)^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex A psi with A under: expect 1 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;


  return pp ;
} /* Z3_A_PsiL_PsiLB__Aunder */

/***********************************************************************************************************************************************/

static POLYNOME Z3_A_PsiL_PsiLB__Bunder (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;

  POLYNOME p5 = vertex_BB_PsiL_PsiRB (b,c,h) ;
  POLYNOME p2 = prop_PsiRB_PsiR (0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p3 = vertex_A_PsiR_PsiRB (a,h) ;
  POLYNOME p4 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p1 = vertex_B_PsiR_PsiLB (d,e,h) ;
  POLYNOME p6 = prop_BB_B (b,c,d,e,2,h) ; /* 1/(k+p+q)^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex A psi with B under: expect je sais pas \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;


  return pp ;
} /* Z3_A_PsiL_PsiLB__Bunder */

/***********************************************************************************************************************************************/

static POLYNOME Z3_A_PsiL_PsiLB__Hover (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  int ppv[4] = {2,1,1,0} ;
  POLYNOME p1 = vertex_H_PsiR_PsiLB (h) ;
  POLYNOME p2 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p3 = vertex_HB_PsiL_PsiRB (h) ;
  POLYNOME p4 = prop_HB_H (0,h) ; /* 1/(k)^2 */
  POLYNOME p5 = vertex_A_H_HB (a, ppv,h) ;
  POLYNOME p6 = prop_HB_H (2,h) ; /* 1/(k+p+q)^2 */

  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex A psi with H over: expect 1/2 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;


  return pp ;
} /* Z3_A_PsiL_PsiLB__Hover */

/***********************************************************************************************************************************************/

static POLYNOME Z3_A_PsiL_PsiLB__Ghostover (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  int ppv[4] = {2,1,1,0} ;
  POLYNOME p1 = vertex_H_PsiR_PsiLB (h) ;
  POLYNOME p2 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p3 = vertex_HB_PsiL_PsiRB (h) ;
  POLYNOME p4 = prop_HB_H (0,h) ; /* 1/(k)^2 */
  POLYNOME p5 = vertex_A_H_HB (a, ppv,h) ;
  POLYNOME p6 = prop_HB_H (2,h) ; /* 1/(k+p+q)^2 */

  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex A psi with Ghost  over: expect 1/2 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;


  return pp ;
} /* Z3_A_PsiL_PsiLB__Ghostover */

/***********************************************************************************************************************************************/

static POLYNOME Z3_A_PsiL_PsiLB__Aover (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  
  /*
    int ppva[4] = {0, -1, -1, 0} ;
  int ppvb[4] = {-1, 0, 0, 0} ;
  int ppve[4] = {1, 1, 1, 0} ;
  */
  int ppva[4] = {0, 0, 0, 0} ;
  int ppvb[4] = {-1, 0, 0, 0} ;
  int ppve[4] = {1, 0, 0, 0} ;


  POLYNOME p1 = vertex_A_PsiL_PsiLB (d,h) ;
  POLYNOME p2 = prop_PsiLB_PsiL (1,h) ;   /* (1/(k+p)^2 */
  POLYNOME p3 = vertex_A_PsiL_PsiLB (c,h) ;
  POLYNOME p4 = prop_AA (b,c,0,h) ; /* 1/(k)^2 */
  POLYNOME p5 = vertex_A_A_A (a,b,e,ppva,ppvb,ppve,h) ;
  POLYNOME p6 = prop_AA (d,e,2,h) ; /* 1/(k+p+q)^2 */

  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex A psi with A over: expect 1 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;


  return pp ;
} /* Z3_A_PsiL_PsiLB__Aover */

/***********************************************************************************************************************************************/

static POLYNOME Z3_A_PsiL_PsiLB__Bover (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;
  short g = newDummyIndex () ;
  short h1 = newDummyIndex () ;
  short i = newDummyIndex () ;
  
  /*
    int ppvcd[4] = {-1,0,0,0} ;
    int ppvhi[4] = {1,1,1,0} ;  // the diagram is log divergent, we do not need the external p,q 
  */
  int ppvcd[4] = {-1,0,0,0} ;
  int ppvhi[4] = {1,0,0,0} ;
      
  POLYNOME p1 = vertex_B_PsiR_PsiLB (f,g,h) ;
  POLYNOME p2 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k+p)^2 */
  POLYNOME p3 = vertex_BB_PsiL_PsiRB (d,e,h) ;
  POLYNOME p4 = prop_BB_B (d,e,b,c,0,h) ; /* 1/(k)^2 */
  POLYNOME p5 = vertex_A_B_BB (a,b,c,h1,i,ppvcd,ppvhi,h) ;
  POLYNOME p6 = prop_BB_B (h1,i,f,g,20,h) ; /* 1/(k+p+q)^2 */

  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;
  /*   POLYNOME ppp2[7] = {p1, p2, p5, p6, 0} ; */

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex A psi with tensor over: expect 1 \n") ;
  showPol (pp) ;
  pp = expand (pp) ;
  if (0) showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;


  return pp ;
} /* Z3_A_PsiL_PsiLB__Bover */

/***********************************************************************************************************************************************/

static POLYNOME Z3_B_PsiR_PsiLB__BBA (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short mu = newDummyIndex () ;
  short nu = newDummyIndex () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  
  short f = newDummyIndex () ;

  int mm1[4] = {0,0,0,0} ;
  int mm2[4] = {1,0,0,0} ;


  POLYNOME p1 = prop_PsiRB_PsiR (0,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p2 = vertex_A_PsiR_PsiRB (a,h) ;
  POLYNOME p3 = prop_AA (a,b,0,h) ; /* 1/k^2 */
  POLYNOME p4 = vertex_A_B_BB (b,mu,nu,c,d,mm1,mm2,h) ;
  POLYNOME p5 = prop_BB_B (c,d,e,f, 0,h) ; /* 1/k^2 */
  POLYNOME p6 = vertex_B_PsiR_PsiLB (e,f,h) ;
  POLYNOME ppp[7] = {p1,p2,p3,p4,p5,p6,0} ;
  POLYNOME pp = polMultiProduct (h,ppp) ;


  printf ("Z3 Classic vertex A/B contrib to vertex B_ab -psi-psi : \n") ;
  showPol(p6) ;
  showPol(pp) ;
  pp = expand (pp) ;
  showPol(pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;
  return pp ;
} /* Z3_B_PsiR_PsiLB__BBA */

/***********************************************************************************************************************************************/

static POLYNOME Z3_H_PsiR_PsiLB__Aunder (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;

  POLYNOME p1 = vertex_A_PsiL_PsiLB (b,h) ;
  POLYNOME p2 = prop_PsiLB_PsiL (2,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p3 = vertex_H_PsiR_PsiLB (h) ;
  POLYNOME p4 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k)^2 */  
  POLYNOME p5 = vertex_A_PsiR_PsiRB (a,h) ;
  POLYNOME p6 = prop_AA (a,b,0,h) ;   /* (1/(k)^2 */ 
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex H psi Aunder: expect 4 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3_H_PsiR_PsiLB__Aunder */

/***********************************************************************************************************************************************/

static POLYNOME Z3_H_PsiR_PsiLB__HAH (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;
  int kk[4] = {-1,0,0,0} ;

  POLYNOME p1 = vertex_H_PsiR_PsiLB (h) ;
  POLYNOME p2 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p3 = vertex_A_PsiR_PsiRB (f,h) ;
  POLYNOME p4 = prop_HB_H (2,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p5 = prop_AA (e,f,0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p6 = vertex_A_H_HB (e,kk,h) ;
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 New vertex H psi HAB: expect -3 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3_H_PsiR_PsiLB__HAH */

/***********************************************************************************************************************************************/

static POLYNOME Z3_H_PsiR_PsiLB__HHA (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;
  int kk[4] = {-1,0,0,0} ;

  
  POLYNOME p1 = vertex_A_PsiR_PsiRB (f,h) ;
  POLYNOME p2 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p3 = vertex_H_PsiR_PsiLB (h) ;
  POLYNOME p4 = prop_HB_H (0,h) ;    /* (1/(k)^2 */ 
  POLYNOME p5 = prop_AA (e,f,2,h) ;  /* (1/(k+p+q)^2 */ 
  POLYNOME p6 = vertex_A_H_HB (e,kk,h) ;
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 New vertex H psi HAB: expect -3 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3_H_PsiR_PsiLB__HAA */

/***********************************************************************************************************************************************/

static POLYNOME Z3_H_PsiR_PsiLB__HAB (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;
  int kk[4] = {-1,0,0,0} ;

  POLYNOME p1 = vertex_B_PsiR_PsiLB (c,d,h) ;
  POLYNOME p2 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p3 = vertex_A_PsiR_PsiRB (f,h) ;
  POLYNOME p4 = prop_BB_B (a,b,c,d,2,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p5 = prop_AA (e,f,0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p6 = vertex_A_H_BB (e,a,b,kk,h) ; /* 1/k^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 New vertex H psi HAB: expect -3 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3_H_PsiR_PsiLB__HAB */

/***********************************************************************************************************************************************/

static POLYNOME Z3_H_PsiR_PsiLB__HBA (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;
  int kk[4] = {1,1,1,0} ;

  POLYNOME p1 = vertex_A_PsiL_PsiLB (f,h) ;
  POLYNOME p2 = prop_PsiLB_PsiL (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p3 = vertex_B_PsiR_PsiLB (c,d,h) ;

  POLYNOME p4 = prop_BB_B (a,b,c,d,0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p5 = prop_AA (e,f,2,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p6 = vertex_A_B_HB (e,a,b,kk,h) ; /* 1/k^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 New vertex H psi HBA: expect -3 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3_H_PsiR_PsiLB__HBA */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME Z3_A_PsiR_PsiLB__Hunder (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;

  POLYNOME p1 = vertex_A_PsiL_PsiLB (b,h) ;
  POLYNOME p2 = prop_PsiLB_PsiL (2,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p3 = vertex_H_PsiR_PsiLB (h) ;
  POLYNOME p4 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k)^2 */  
  POLYNOME p5 = vertex_A_PsiR_PsiRB (a,h) ;
  POLYNOME p6 = prop_HB_H (0,h) ;   /* (1/(k)^2 */ 
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex A psi Hunder: expect 4 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3_H_PsiR_PsiLB__Hunder */

/***********************************************************************************************************************************************/

static POLYNOME Z3_A_PsiR_PsiLB__Aunder (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;

  POLYNOME p1 = vertex_A_PsiL_PsiLB (b,h) ;
  POLYNOME p2 = prop_PsiLB_PsiL (2,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p3 = vertex_H_PsiR_PsiLB (h) ;
  POLYNOME p4 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k)^2 */  
  POLYNOME p5 = vertex_A_PsiR_PsiRB (a,h) ;
  POLYNOME p6 = prop_AA (a,b,0,h) ;   /* (1/(k)^2 */ 
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex H psi Aunder: expect 4 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3_H_PsiR_PsiLB__Aunder */

/***********************************************************************************************************************************************/

static POLYNOME Z3_A_PsiR_PsiLB__Bunder (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;

  POLYNOME p1 = vertex_A_PsiL_PsiLB (b,h) ;
  POLYNOME p2 = prop_PsiLB_PsiL (2,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p3 = vertex_H_PsiR_PsiLB (h) ;
  POLYNOME p4 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k)^2 */  
  POLYNOME p5 = vertex_A_PsiR_PsiRB (a,h) ;
  POLYNOME p6 = prop_AA (a,b,0,h) ;   /* (1/(k)^2 */ 
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex H psi Aunder: expect 4 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3_H_PsiR_PsiLB__Bunder */

/***********************************************************************************************************************************************/

static POLYNOME Z3_A_PsiR_PsiLB__HAH (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;
  int kk[4] = {-1,0,0,0} ;

  POLYNOME p1 = vertex_H_PsiR_PsiLB (h) ;
  POLYNOME p2 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p3 = vertex_A_PsiR_PsiRB (f,h) ;
  POLYNOME p4 = prop_HB_H (2,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p5 = prop_AA (e,f,0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p6 = vertex_A_H_HB (e,kk,h) ;
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 New vertex H psi HAB: expect -3 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3_H_PsiR_PsiLB__HAH */

/***********************************************************************************************************************************************/

static POLYNOME Z3_H_PsiR_PsiLB__HHdoublon (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;
  short f = newDummyIndex () ;
  int kk[4] = {-1,0,0,0} ;

  
  POLYNOME p1 = vertex_A_PsiR_PsiRB (f,h) ;
  POLYNOME p2 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p3 = vertex_H_PsiR_PsiLB (h) ;
  POLYNOME p4 = prop_HB_H (0,h) ;    /* (1/(k)^2 */ 
  POLYNOME p5 = prop_AA (e,f,2,h) ;  /* (1/(k+p+q)^2 */ 
  POLYNOME p6 = vertex_A_H_HB (e,kk,h) ;
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 New vertex H psi HAB: expect -3 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3_H_PsiR_PsiLB__HAA */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME Z3_A_c_cB__Aunder (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  int ppva[4] = {-1, -1, -1, 0} ;
  int ppvb[4] = {-1, 0, 0, 0} ;
  int ppvc[4] = {0, 0, -1, 0} ;

  POLYNOME p1 = vertex_A_c_cB (c, ppvc,h) ;
  POLYNOME p2 = prop_AA (b, c, 1,h) ;
  POLYNOME p3 = vertex_A_c_cB (b, ppvb,h) ;
  POLYNOME p4 = prop_cB_c (0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p5 = vertex_A_c_cB (a, ppva,h) ;
  POLYNOME p6 = prop_cB_c (2,h) ; /* 1/(k+p+q)^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex A c cB ghost with A under: expect 1/2 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;


  return pp ;
} /* Z3_A_c_cB__Aunder */

/***********************************************************************************************************************************************/

static POLYNOME Z3_A_c_cB__Aover (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;

  int ppva[4] = {0, -1, -1, 0} ;
  int ppvb[4] = {-1, 0, 0, 0} ;
  int ppvc[4] = {1, 1, 0, 0} ;
  int ppvd[4] = {0, 0, -1, 0} ;
  int ppve[4] = {1, 1, 1, 0} ;


  
  POLYNOME p1 = vertex_A_c_cB (d, ppvd,h) ;
  POLYNOME p2 = prop_cB_c (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p3 = vertex_A_c_cB (c, ppvc,h) ;
  POLYNOME p4 = prop_AA (b, c, 0,h) ;
  POLYNOME p5 = vertex_A_A_A (a,b,e,ppva,ppvb,ppve,h) ;
  POLYNOME p6 = prop_AA (d, e, 2,h) ;
  
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex A c cB ghost with A over: expect 1/2 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;


  return pp ;
} /* Z3_A_c_cB__Aover */

/***********************************************************************************************************************************************/

static POLYNOME Z3_A_H_HB__Aunder (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  int ppva[4] = {-2, -1, -1, 0} ;
  int ppvb[4] = {-1, 1, 0, 0} ;
  int ppvc[4] = {-1, -1, -2, 0} ;

  POLYNOME p1 = vertex_A_H_HB (c, ppvc,h) ;
  POLYNOME p2 = prop_AA (b, c, 1,h) ;
  POLYNOME p3 = vertex_A_H_HB (b, ppvb,h) ;
  POLYNOME p4 = prop_HB_H (0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p5 = vertex_A_H_HB (a, ppva,h) ;
  POLYNOME p6 = prop_HB_H (2,h) ; /* 1/(k+p+q)^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex A H HB scalar with A under: expect 1/2 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;


  return pp ;
} /* Z3_A_H_HB__Aunder */

/***********************************************************************************************************************************************/

static POLYNOME Z3_A_H_HB__Aover (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;

  int ppva[4] = {0, -1, -1, 0} ;
  int ppvb[4] = {-1, 0, 0, 0} ;
  int ppvc[4] = {1, 2, 0, 0} ;
  int ppvd[4] = {1, 1, -1, 0} ;
  int ppve[4] = {1, 1, 1, 0} ;

  POLYNOME p1 = vertex_A_H_HB (d, ppvd,h) ;
  POLYNOME p2 = prop_HB_H (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p3 = vertex_A_H_HB (c, ppvc,h) ;
  POLYNOME p4 = prop_AA (b, c, 0,h) ;
  POLYNOME p5 = vertex_A_A_A (a,b,e,ppva,ppvb,ppve,h) ;
  POLYNOME p6 = prop_AA (d, e, 2,h) ;
  
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex A H HB scalar with A over: expect 1/2 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;


  return pp ;
} /* Z3_A_H_HB__Aover */

/***********************************************************************************************************************************************/

static POLYNOME Z3_B_PsiR_PsiLB__BAH (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  int kk[4] = {-1,0,0,0} ;

  POLYNOME p1 = vertex_H_PsiR_PsiLB (h) ;
  POLYNOME p2 = prop_PsiRB_PsiR (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p3 = vertex_A_PsiR_PsiRB (d,h) ;
  POLYNOME p4 = prop_HB_H (2,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p5 = prop_AA (c,d,0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p6 = vertex_A_B_HB (c,a,b,kk,h) ; /* 1/k^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 New vertex B psi BAH: expect sB_ab/2 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3_B_PsiR_PsiLB__BAH */

/***********************************************************************************************************************************************/

static POLYNOME Z3_B_PsiR_PsiLB__BHA (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  int kk[4] = {1,1,1,0} ;

  POLYNOME p1 = vertex_A_PsiL_PsiLB (d,h) ;
  POLYNOME p2 = prop_PsiLB_PsiL (1,h) ;   /* (1/(k+p)^2 */ 
  POLYNOME p3 = vertex_H_PsiR_PsiLB (h) ;

  POLYNOME p4 = prop_HB_H (0,h) ;   /* (1/(k)^2 */ 
  POLYNOME p5 = prop_AA (c,d,2,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p6 = vertex_A_B_HB (c,a,b,kk,h) ; /* 1/{k+p+q)^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 New vertex B psi BAH: expect sB_ab/2 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3_B_PsiR_PsiLB__BHA */

/***********************************************************************************************************************************************/

static POLYNOME Z3_B_PsiR_PsiLB__Aunder (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  POLYNOME  pp = 0 ;
  
  if (1)
    {
      POLYNOME p1 = vertex_A_PsiL_PsiLB (d,h) ;
      POLYNOME p2 = prop_PsiLB_PsiL (0,h) ;   /* (1/(k+p+q)^2 */ 
      POLYNOME p3 = vertex_B_PsiR_PsiLB (a,b,h) ;
      POLYNOME p4 = prop_PsiRB_PsiR (0,h) ;   /* (1/(k)^2 */ 
      POLYNOME p5 = vertex_A_PsiR_PsiRB (c,h) ;
      POLYNOME p6 = prop_AA (c,d,0,h) ; /* 1/(k+p)^2 */
      POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

      pp = polMultiProduct (h,ppp) ;
    }
  else
    {
pp = newScalar (1,h) ;
      tcpy (pp->tt.sigma, "efgh") ;
      tcpy (pp->tt.g, "afbg") ;
    }
  printf ("Z3 Classic vertex B psi with A under: expect 0 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3_B_PsiR_PsiLB__A_under */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME Z3H_A_PsiR_PsiLB (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;

  POLYNOME p1 = vertex_H_PsiR_PsiLB (h) ;
  POLYNOME p2 = prop_PsiLB_PsiL (0,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p3 = vertex_A_PsiL_PsiLB (a,h) ;
  POLYNOME p4 = prop_PsiLB_PsiL (0,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p5 = vertex_HB_PsiL_PsiRB (h) ;
  POLYNOME p6 = prop_HB_H (0,h) ; /* 1/k^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex A psi with H under: expect 0 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3H_A_PsiR_PsiLB */

/***********************************************************************************************************************************************/

static POLYNOME Z3A_A_PsiR_PsiLB (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;

  POLYNOME p1 = vertex_A_PsiR_PsiRB (b,h) ;
  POLYNOME p2 = prop_PsiRB_PsiR (0,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p3 = vertex_A_PsiR_PsiRB (a,h) ;
  POLYNOME p4 = prop_PsiRB_PsiR (0,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p5 = vertex_A_PsiR_PsiRB (c,h) ;
POLYNOME p6 = prop_AA (b,c,0,h) ; /* 1/k^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex A psi with A under: expect 0 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3A_A_PsiR_PsiLB */

/***********************************************************************************************************************************************/

static POLYNOME Z3B_A_PsiR_PsiLB (void)
{
  AC_HANDLE h = ac_new_handle () ;
  short a = newDummyIndex () ;
  short b = newDummyIndex () ;
  short c = newDummyIndex () ;
  short d = newDummyIndex () ;
  short e = newDummyIndex () ;

  POLYNOME p1 = vertex_B_PsiR_PsiLB (b,c,h) ;
  POLYNOME p2 = prop_PsiLB_PsiL (0,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p3 = vertex_A_PsiL_PsiLB (a,h) ;
  POLYNOME p4 = prop_PsiLB_PsiL (0,h) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p5 = vertex_BB_PsiL_PsiRB (d,e,h) ;
  POLYNOME p6 = prop_BB_B (d,e,b,c,0,h) ; /* 1/k^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;
  short w = newDummyIndex () ;
  POLYNOME pp = polMultiProduct (h,ppp) ;
  printf ("Z3 Classic vertex A psi with B under: expect 0 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;
  pp = expand (pp) ;
  pp = contractIndices (pp) ;
  showPol(pp) ;

  showPol (pp) ;

  printf ("Z3 TEST \n") ;
  p1 = newSigma ('a',h) ;
  tcpy (p1->tt.sigma, "bdecbde") ;
  showPol(p1) ;
  p1 = contractIndices (p1) ;
  showPol(p1) ;
  p1 = pauliCleanUp (p1, w) ;
  showPol (p1) ;
 
  return pp ;
} /* Z3B_A_PsiR_PsiLB */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static BOOL polynomeTest (void)
{
  AC_HANDLE h = ac_new_handle () ;
  POLYNOME p1, p2, p3, p4, pp, gac, sga, sgb, sbc , sbd ;
  BOOL debug = FALSE ;
  short alpha ;
  sga = newSigma ('a',h) ;
  sgb = newSigma ('b',h) ;
  sbc = newSigB ('c',h) ;
  sbd = newSigB ('d',h) ;
  gac = newG ('a','c',h) ;

p1 = polSum (sga, sgb,h) ;
 p2 = polSum (sbc, sbd,h) ;
 p4 = polProduct (p1,p2,h) ;
 p3 = polProduct (gac, p4,h) ;

  if (debug) showPol (p3) ;
  p4 = expand (p3) ;
  if (debug) showPol (p4) ;
  contractIndices (p4) ;
  if (debug) showPol (p4) ;
  pp = contractProducts (p4) ;
  if (debug) showPol (p4) ;
  contractIndices (p4) ;
  if (debug) showPol (p4) ;
  pauliTrace (p4, 0) ;
  if (debug) showPol (p4) ;
  contractIndices (p4) ;

  printf ("Test  g_ac (s_ + s_b) * (s_c + s_d)\n") ;
  if (debug) showPol (p4) ;



  printf ("Test  g_ab p_a p_c\n") ;
  p1 = newG ('a','b',h) ;
  p2 = newSigB ('a',h) ;
  p3 = newSigB ('c',h) ;
  p4 = polProduct (p1,p2,h) ;
  pp = polProduct (p4,p3,h) ;
  if (debug) showPol (pp) ;
  pp = contractProducts (pp) ;
  contractIndices (pp) ;
  if (debug) showPol (pp) ;

  p1 = newK ('a',h) ;
  p1->tt.mm[0][1] = 'b' ;
  p1->tt.mm[0][2] = 'c' ;
  p1->tt.mm[0][3] = 'd' ;
  p1->tt.denom[0] = 4 ;
  showPol(p1) ;
p4 = polCopy(p1,h) ;
  p2 = dimIntegral (p1) ;
  showPol(p2) ;
  p3 = polProduct (newAG('a','b','e','f',0,h),p2,h) ;
  showPol(p3) ;
  p3 = expand(p3) ;
  showPol(p3) ;
  p3 = polProduct (newG('e','f',h),p3,h) ;
  showPol(p3) ;
  p3 = expand(p3) ;
  showPol(p3) ;

  showPol(p4) ;
  p3 = polProduct (newG('c','d',h),p4,h) ;
  showPol(p3) ;
  p2 = dimIntegral (p3) ;
  showPol(p2) ;

  printf ("##### Test des integrales  /(k+p)^2 expect 1\n") ;
  p1 = newScalar (1,h) ;
  p1->tt.mm[0][0] = newDummyIndex () ;
  p1->tt.denom[0] = 1 ;
  p1->tt.denom[1] = 0 ;
  p1->tt.denom[2] = 1 ;
  showPol(p1) ;
  p2 = dimIntegral (p1) ;
  showPol(p2) ;
  alpha = newDummyIndex () ;
  p3 = momentaCleanUp (p2, alpha) ;
  printf ("##### expect -(p+q)/2 \n") ;
  showPol(p3) ;

  printf ("##### test k3/k^6 cases\n") ;
  p1 = newScalar (1,h) ;
  p1->tt.mm[0][0] = newDummyIndex () ;
  p1->tt.mm[0][1] = newDummyIndex () ;
  p1->tt.mm[0][2] = newDummyIndex () ;
  p1->tt.denom[0] = 1 ;
  p1->tt.denom[1] = 1 ;
  p1->tt.denom[2] = 1 ;


  showPol(p1) ;

  p2 = dimIntegral (p1) ;
  showPol(p2) ;
  alpha = newDummyIndex () ;
  p3 = momentaCleanUp (p2, alpha) ;
  printf ("##### check please \n") ;
  showPol(p3) ;

  printf ("##### pauli traces\n") ;
  p1 = newSigma ('a',h) ;
  p1->tt.sigma[0] = 'a' ;
  p1->tt.sigma[1] = 'b' ;

  showPol (p1) ;
  p2 = pauliTrace (p1, 0) ;
  showPol (p2) ;

  p1 = newSigma ('a',h) ;
  p1->tt.sigma[0] = 'a' ;
  p1->tt.sigma[1] = 'b' ;
  p1->tt.sigma[2] = 'c' ;
  p1->tt.sigma[3] = 'd' ;

  showPol (p1) ;
  p2 = pauliTrace (p1, 0) ;
  showPol (p2) ;

  p1 = newSigma ('a',h) ;
  p1->tt.sigma[0] = 'a' ;
  p1->tt.sigma[1] = 'b' ;
  p1->tt.sigma[2] = 'c' ;
  p1->tt.sigma[3] = 'd' ;
  p1->tt.sigma[4] = 'e' ;
  p1->tt.sigma[5] = 'f' ;

  showPol (p1) ;
  p2 = pauliTrace (p1, 0) ;
  showPol (p2) ;
  

  printf ("#########################################\n") ;
  firstDummyIndex = 'a' ;
  if (1)
    {
      short a, b, c, d ;
      pp = newScalar (1/2.0,h) ; /* we derive twice in p_ab */
      pp->tt.mm[0][0] = a = newDummyIndex () ;
      pp->tt.mm[0][1] = b = newDummyIndex () ;
      pp->tt.denom[0] = 1 ;
      pp->tt.denom[1] = 1 ;
      if (0)
	{
	  pp->tt.g[0] = c = newDummyIndex () ;
	  pp->tt.g[1] = d = newDummyIndex () ;
	}
      printf ("# test des integrales  (g_ac k_bd /k^2 (k+p)^2)\n") ;
      showPol(pp) ;
      pp = dimIntegral (pp) ;
      printf ("# Expect 1/3 p_ab - 1/12 g_ab p_cc\n") ;
      showPol(pp) ;
      printf ("# trace it and expect zero\n") ;
      p1 = newScalar (1,h) ;
      p1->tt.g[0] = a ;
      p1->tt.g[1] = b ;
      p2 = polProduct (p1, pp,h) ;
      showPol(p1) ;
      p2 = expand(p2) ;
      p2 = contractIndices(p2) ;
      p2 = expand(p2) ;
      showPol(p2) ;
      p2 = squareMomentaCleanUp (p2) ;
      showPol(p2) ;
    }
  
  
  firstDummyIndex = 'a' ;
  if (0)
    {
      short a,b,c,d ;
      pp = newScalar (1,h) ;
      pp->tt.mm[0][0] = a = newDummyIndex () ;
      pp->tt.mm[0][1] = b = newDummyIndex () ;
      pp->tt.mm[0][2] = c = newDummyIndex () ;
      pp->tt.mm[0][3] = d = newDummyIndex () ;
      pp->tt.denom[0] = 1 ;
      pp->tt.denom[1] = 1 ;
      printf ("# test des integrales k_abcd / k^2 (k+p)2\n") ;
      showPol(pp) ;
      pp = dimIntegral (pp) ;
      showPol(pp) ;
      printf ("#  FAUX on voit des indicies libres nouveaux\n") ;
    }
  firstDummyIndex = 'a' ;
  if (1)
    {
      short a,b,c,d ;
      pp = newScalar (1,h) ;
      pp->tt.mm[0][0] = a = newDummyIndex () ;
      pp->tt.mm[0][1] = b = newDummyIndex () ;
      pp->tt.mm[0][2] = c = newDummyIndex () ;
      pp->tt.mm[0][3] = d = newDummyIndex () ;
      pp->tt.denom[0] = 2 ;
      pp->tt.denom[1] = 2 ;
      printf ("# test des integrales k_abcd / k^4 (k+p)4\n") ;
      showPol(pp) ;
      pp = dimIntegral (pp) ;
      showPol(pp) ;
      printf ("# trace it, expect 1: CORRECT\n") ;
      p1 = newScalar (1,h) ;
      p1->tt.g[0] = a ;
      p1->tt.g[1] = b ;
      p1->tt.g[2] = c ;
      p1->tt.g[3] = d ;
      p2 = polProduct (p1, pp,h) ;
      showPol(p1) ;
      p2 = expand(p2) ;
      p2 = contractIndices(p2) ;
      p2 = expand(p2) ;
      showPol(p2) ;
    }

  firstDummyIndex = 'a' ;
  if (1)
    {
      short a,b,c,d ;
      pp = newScalar (1,h) ;
      pp->tt.mm[0][0] = a = newDummyIndex () ;
      pp->tt.mm[0][1] = b = newDummyIndex () ;
      pp->tt.mm[0][2] = c = newDummyIndex () ;
      pp->tt.mm[0][3] = d = newDummyIndex () ;
      pp->tt.denom[0] = 2 ;
      pp->tt.denom[1] = 1 ;
      printf ("# test des integrales k_abcd / k^4 (k+p)2\n") ;
      showPol(pp) ;
      pp = dimIntegral (pp) ;
      showPol(pp) ;
      printf ("# FAUX je veux du p_ab - g_ab p^2 ou un trc du genre , pas du pur p^2\n") ;
      printf ("# trace it, expect 1: CORRECT\n") ;
      p1 = newScalar (1,h) ;
      p1->tt.g[0] = a ;
      p1->tt.g[1] = b ;
      p1->tt.g[2] = c ;
      p1->tt.g[3] = d ;
      p2 = polProduct (p1, pp,h) ;
      showPol(p1) ;
      p2 = expand(p2) ;
      p2 = contractIndices(p2) ;
      p2 = expand(p2) ;
      showPol(p2) ;
    }

  firstDummyIndex = 'a' ;
  if (1)
    {
      short a,b,c,d,e,f ;
pp = newScalar (1,h) ;
      pp->tt.mm[0][0] = a = newDummyIndex () ;
      pp->tt.mm[0][1] = b = newDummyIndex () ;
      pp->tt.mm[0][2] = c = newDummyIndex () ;
      pp->tt.mm[0][3] = d = newDummyIndex () ;
      pp->tt.mm[0][4] = e = newDummyIndex () ;
      pp->tt.mm[0][5] = f = newDummyIndex () ;
      pp->tt.denom[0] = 2 ;
      pp->tt.denom[1] = 2 ;
      pp->tt.denom[2] = 1 ;
      printf ("# test des integrales k_abcdefgh / k^4 (k+p)4 (k+p+q)^4\n") ;
      showPol(pp) ;
      pp = contractIndices (pp) ;
      showPol(pp) ;
      pp = dimIntegral (pp) ;
      showPol(pp) ;
      printf ("# trace it, expect 1: CORRECT\n") ;
      p1 = newScalar (1,h) ;
      p1->tt.g[0] = a ;
      p1->tt.g[1] = b ;
      p1->tt.g[2] = c ;
      p1->tt.g[3] = d ;
      p1->tt.g[4] = e ;
      p1->tt.g[5] = f ;
      showPol(p1) ;
      p1 = contractIndices (p1) ;
      showPol(p1) ;

p2 = polProduct (p1, pp,h) ;

      showPol(p2) ;
      p2 = expand(p2) ;
      showPol(p2) ;
      p2 = contractIndices(p2) ;
      p2 = expand(p2) ;
      showPol(p2) ;
    }

  firstDummyIndex = 'a' ;
  if (1)
    {
      short a,b,c,d ;
      pp = newScalar (1,h) ;
      pp->tt.eps[0] = a = newDummyIndex () ;
      pp->tt.eps[1] = b = newDummyIndex () ;
      pp->tt.eps[2] = c = newDummyIndex () ;
      pp->tt.eps[3] = d = newDummyIndex () ;
      printf ("# test des espilon\n") ;
      showPol(pp) ;
      p1 = newScalar (1,h) ;
      p1->tt.g[0] = a ;
      p1->tt.g[1] = b ;
      showPol(p1) ;
      p2 = polProduct (p1, pp, h) ;
      p2 = expand(p2) ;
      p2 = contractIndices(p2) ;
      p2 = expand(p2) ;
      printf ("# g_ab eps_abcd  expect zero \n") ;
      showPol(p2) ;

      p3 = newScalar (1,h) ;
      p3->tt.mm[1][0] = c ;
      p3->tt.mm[1][1] = a ;
      showPol(p3) ;
      p2 = polProduct (p3, pp, h) ;
      p2 = expand(p2) ;
      p2 = contractIndices(p2) ;
      p2 = expand(p2) ;
      printf ("# (p)_ab eps_abcd  expect zero \n") ;
      showPol(p2) ;

      p3->tt.mm[1][0] = b ;
      p2 = polProduct (p3, p1, h) ;
      p2 = expand(p2) ;
      p2 = contractIndices(p2) ;
      p2 = expand(p2) ;
      printf ("# (p)_ab g_ab  expect p_aa \n") ;
      showPol(p2) ;
    }
  return 0 ;
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/*****  SU(2/1) representation theory. This is used by, but does not depend on the analysis above of the Feynman diagrams **********************/
/*****  Casimir studies with Pater Jarvis, mars 2021 **********************************************/
/***** Scalar anomaly paper is below **********************************************/
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

typedef struct kasStruct { MX *mu, *Rmu, *QQ, gg, GG, kas2, CHI, ccc, CCC, cccGhost, CCCGhost, c4, c5, C4, C5, kas3, chi16 ; int a, b, d, d1, d2, d3, d4, chi, scale, NN ; BOOL isOSp, isSU2, isCycle, show, xiPrime ; float zc4, zC4 ; AC_HANDLE h ; } KAS ;
typedef struct comtpStruct { int a, b, c, n, s ; } LC ;
static MX KasCommut (MX a, MX b, int sign, KAS *kas) ;
  
static MX *KasimirConstructSU2Matrices (KAS *kas)
{
  MX  muE, muF, muH  ; /* the 5 generators of OSp(2/1) in the Chevalley basis */
  MX *mu ;
  int d, d1 = 0, d2 = 0 ;
  int i ;
  int a ;
  AC_HANDLE h = kas->h ;
  mu = kas->mu = (MX *) halloc (10 * sizeof (MX), kas->h) ;

  kas->chi = 1 ;
  kas->isSU2 = TRUE ;
  
  a = kas->a = kas->a - 2000 ;
  d1 = a + 1 ;
  d2 = 0 ;
  d = d1 + d2 ;

  kas->d = d = d1 + d2 ;
  kas->d1 = d1 ;
  kas->d2 = d2 ;
  kas->d3 = 0 ;
  kas->d4 = 0 ;
      
  
  muH = mxCreate (h,  "muH", MX_INT, d, d, 0) ;
  muE = mxCreate (h,  "muE", MX_INT, d, d, 0) ;
  muF = mxCreate (h,  "muF", MX_INT, d, d, 0) ;
  
  int xx[d*d] ;

  /* even Cartan operator H = diag (a, a-2 .... -a in each SU(2) sector */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    xx[d * i + i] = d1 - 1 - 2*i ;
  mxSet (muH, xx) ;
  mxShow (muH) ;
 
   /* even raising operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 1 ; i < d1 ; i++)
    xx[d * i + i - 1] = i * (a  - i + 1) ;
  
  mxSet (muE, xx) ;
  mxShow (muE) ;
  
  /* even lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d - 1 ; i++)
    {
      if (i != d1 - 1) xx[d * i + i + 1] = 1 ;
    }
  mxSet (muF, xx) ;
  mxShow (muF) ;

  mu[0] =0 ; mu[1] = muE ; mu[2] = muF ; mu[3] = muH ;
  mu[4] = 0 ; mu[5] = 0 ; mu[6] = 0 ; mu[7] = 0 ;
  mu[8] = 0 ;
  mu[9] = 0 ;

  return mu ;
} /* KasimirConstructSU2Matrices */

static MX *KasimirConstructOSp1_2Matrices (KAS *kas)
{
  MX  muE, muF, muH, muS, muT  ; /* the 5 generators of OSp(2/1) in the Chevalley basis */
  MX *mu ;
  int d, d1 = 0, d2 = 0 ;
  int i, j ;
  int a ;
  AC_HANDLE h = kas->h ;
  mu = kas->mu = (MX *) halloc (10 * sizeof (MX), kas->h) ;
  kas->chi = 1 ;
  kas->isOSp = TRUE ;
  
  a = kas->a = kas->a - 1000 ;
  d1 = a + 1 ;
  d2 = a ;
  d = d1 + d2 ;

  kas->d = d = d1 + d2 ;
  kas->d1 = d1 ;
  kas->d2 = d2 ;
  kas->d3 = 0 ;
  kas->d4 = 0 ;
      
  
  muH = mxCreate (h,  "muH", MX_INT, d, d, 0) ;
  muE = mxCreate (h,  "muE", MX_INT, d, d, 0) ;
  muF = mxCreate (h,  "muF", MX_INT, d, d, 0) ;
  muS = mxCreate (h,  "muS", MX_INT, d, d, 0) ;
  muT = mxCreate (h,  "muT", MX_INT, d, d, 0) ;
  
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
  mxShow (muH) ;
 
   /* even raising operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 1 ; i < d1 ; i++)
    xx[d * i + i - 1] = i * (a  - i + 1) ;
  for (i = 1 ; i < d2 ; i++)
    xx[d * (d1 + i) + d1 + i - 1] =  i * (a - i ) ;
  
  mxSet (muE, xx) ;
  mxShow (muE) ;
  
  /* even lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d - 1 ; i++)
    {
      if (i != d1 - 1) xx[d * i + i + 1] = 1 ;
    }
  mxSet (muF, xx) ;
  mxShow (muF) ;

 /* odd raising operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (i+1)  + d1 + i] = (i + 1) ;
      xx[d * (d1 + i)  + i] = d2 - i ;
    }
  mxSet (muS, xx) ;
  mxShow (muS) ;
  
 /* odd lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (i)  + d1 + i] = -1 ;
      xx[d * (d1 + i)  + i + 1] = 1 ;
    }
  mxSet (muT, xx) ;
  mxShow (muT) ;
  
  mu[0] =0 ; mu[1] = muE ; mu[2] = muF ; mu[3] = muH ;
  mu[4] = muS ; mu[5] = muT ; mu[6] = 0 ; mu[7] = 0 ;
  mu[8] = 0 ;
  mu[9] = 0 ;

  return mu ;
} /* KasimirConstructOSp1_2Matrices */

/***********************************************************************************************************************************************/

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

#ifdef JUNK
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
#endif
    
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

typedef struct fabcStruct {int a,b,c,sign; complex float z; const char *title ;} FABC ;
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

  #ifdef JUNK

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

#endif
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
#ifdef JUNK
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
#endif
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
#ifdef JUNK
  

  
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
#endif

/*************************************************************************************************/
#ifdef JUNK
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
#endif
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

/*************************************************************************************/
/***************************** Public interface *************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// su21: Construction of su(2/1) representations and Feynman diagrams\n"
	    "// Authors: Jean Thierry-Mieg, NCBI, 2020-, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "// Construct the matrices of irreducible and indecomposable representations\n"
	    "// Construct the Casimirs, super Casimirs, Gorelik ghost Casimir\n"
	    "// \n"
	    "// Also compute the anomalies and the feynman diagrams supporting my JHEP su(2/1) papers\n"
	    "//\n"
	    "// Syntax:\n"
	    "// su21 [options]\n"
	    "//   [] [-h] [-help] [--help] : this message\n"
	    "// A: Representations\n"
	    "//   su21 -a <int> -b <int> [-N <int>]\n"
	    "//     export the matrices, Casimirs and verifications for the module with \n"
	    "//     Dynkin lables (a,b), a positive integer, b signed integer\n"
	    "//     Number of generations N (N >= 2)\n"
	    "//       In theory, b can be any complex number,\n"
	    "//     for numerical convenience, we restrict here to signed integers\n"
	    "//     but the formulas like the Casimir eigen values are anlytic in b\n"
	    "//       When a or N are large, many outputs are suppressed, try first a<=3, N<=3\n"
	    "//\n"
	    "// B: Feynman diagrams\n"
	    "//   Not documented, sorry: check the source code !\n"
	    "//\n"
	    "// THETA: -theta  SU(2/1) supergroup det = 1 ?\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  dna2dna --help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  AC_HANDLE h = ac_new_handle () ;

  freeinit () ;

  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  getCmdLineInt (&argc, argv, "-c3Mask", &c3Mask) ;
  getCmdLineInt (&argc, argv, "-t", &myType) ;
  /*   BOOL SU3 = getCmdLineBool (&argc, argv, "-su3") ; */
  int NN = 0 ;
  int CYCLE = 0 ;
  
  getCmdLineInt (&argc, argv, "-N", &NN) ; /* Number of generations >= 2 */
  getCmdLineInt (&argc, argv, "-NN", &NN) ; /* synonim */
  BOOL jacobi = getCmdLineBool (&argc, argv, "-jacobi") ;
    
  if (getCmdLineBool (&argc, argv, "-R16"))  /* lpto-quark real rep */
    {
      KasimirR16 () ;
      return 0 ;
    }

  if (getCmdLineBool (&argc, argv, "-bbb"))  /* self dual tensors potential */
    {
      BBB ();
      return 0 ;
    }

  if (getCmdLineBool (&argc, argv, "-theta"))  /* su(2/1) group determinant */
    {
      THETA ();
      return 0 ;
    }
      
  getCmdLineInt (&argc, argv, "-cycle", &CYCLE) ;

  int a = 0, b = 0 ;
  int king = 0 ;

  getCmdLineInt (&argc, argv, "-a", &a) ;
  getCmdLineInt (&argc, argv, "-b", &b) ;
  getCmdLineInt (&argc, argv, "-king", &king) ;

  if (0)
    {      /* a test to time the and evaluate the ranfloat function */
      int i = a ;
      int n = 0, nn = 0 ;
      float z ;
      
      while (i--)
	{
	  z = randfloat () ;
	  nn++ ;
	  if (z > .5) n++ ;
	}
      printf("#RAND nn=%d n=%d f=%g\n", nn , n , 2.0*n/nn - 1.0) ;
      exit (0) ;

    }

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
  if (NN != 0 && NN < 2)
    usage ("The number of generations N should be an integer >= 2") ;
  

  if (0)
    { /* sum of fibonnaci numbers, Euler set of problems as D programming language example,  2021_09_03 */
      int n[2] = {1,2} ;
      int j, f=0 ;
      int N = 4000000 ;
      long int s = 3 ;

      for (j=0;f<N;j=1-j)
	{
	  f=n[j]+=n[1-j] ;
	  s+= f;
	}
      printf ("f=%d s=%ld\n",n[j],s-f) ;
      exit(0) ;
    }


  if (CYCLE)
    {
      marcuCycle (NN, a, b) ;
      exit (0) ;
    }

  if (king > 0) /* a test */
    { /* check the non Abelian expansion exp(a)exp(b)exp(-b) = exp (b + [a,b] + [a,[a,b]]/2! + [a,[a[a,b]]]/3! ...) */
      superExponential (king, a, b, h) ;
      exit (0) ;
    }

  if (getCmdLineBool (&argc, argv, "-powerSum"))
    { /* find the value of sum n^k */
      powerSum() ;
      exit (0) ;
    }

  if (!NN && (a || b))
    {
      /* 2021_03_18 
       * construct the 8 matrices for the generic irreps of su(2/1) with h.w. (a,b)
       * verify all commutations relations
       * compute the casimir tensors and operators 

       */
      Kasimirs (1,1, FALSE) ;
      Kasimirs (1,0, FALSE) ;
      Kasimirs (a,b, TRUE) ;
      exit (0) ;
    }
  /* always init, otherwise the gcc linker is unhappy */
  if (0) muInit (h) ;   /* init the 4x4 matrices */
  if (0) muInit2 (h) ;  /* init the 2-families 8x8 rotated matrices */
  if (NN >= 2) muInitNMarcu (a,b, NN) ;  /* init the 2-families 8x8 marcu indecomposable matrices */
  /* verification numerique directe de traces de matrices de pauli */ 
  if (0) { muSigma (h) ; exit (0);}

  /* Verifications des traces sur la theorie des groupes pour l'article sur les anomalies scalaires */
  if (getCmdLineBool (&argc, argv, "-G"))
    {
      muConjugate (h) ;
      
      
      printf ("########## Compute the relevant traces of products of 2,3,4 SU(2/1) matrices\n") ;
      printf ("########## In each case, the trace is computed for the neutral representation (N), then for leptons (e), quarks (q) and family (e+3*q)\n") ;
      printf ("########## The observation is that leptons and quarks have anomalous traces, but they compensate each other\n") ;
      printf ("########## The family trace, one lepton +  quarks, is proportional to the neutral trace\n") ;
      printf ("########## In the last column, we check that S = e + 3*q - 4*n == 0\n") ;
      
      printf ("########## Verify the commutators,   all computed norms should vanish\n");
      muStructure () ;
    
      
      if (0) mu2p ("######### Metric\n# For the even generators (a,b=0123), compute the Super-Trace: STr(ab)\n# For the odd generators (i=4567), compute the Left trace: LTr(ij)\n We hope to find the SU(2/1) Super-Killing metric") ;
      
      if (1) casimir2 ("######### Casimir 2\n# 1/2 g^AB mu_A mu_B,   we hope to find a diagonal matrix") ;
      
      if (1) casimir3 ("######### Super Casimir 3\n# 1/6 d^ABC mu_A mu_B mu_C,   we hope to find a diagonal matrix", FALSE) ;
      if (0) casimir3 ("######### Hyper Casimir 3\n# 1/6 d^ABC mu_A mu_B mu_C,   we hope to find a diagonal matrix", TRUE) ;
      
      exit (0) ;
    }

  if (getCmdLineBool (&argc, argv, "-F"))
    {
      if (0)
	{
	  POLYNOME pp = newAG6 (0,0,0,0,0,0,h) ;
	  showPol (pp) ;
	  firstDummyIndex = 'a' ;
	  printf ("============t'Hooft, check some integrals\n") ;
	  if (0) Thooft () ;
	  
	  printf ("============Hodge, check some projectors\n") ;
	  if (1) Hodge () ;

	  printf ("============polynomeTest\n") ;
	  if (0) polynomeTest () ;
	  exit (0) ;
	}

      /* dimIntegral test, expect 1 each time */
      if (0)
	{
	  printf ("============dimIntegralTest: expect 1 each time\n") ;
	  firstDummyIndex = 'a' ;
	  char a = newDummyIndex () ;
	  char b = newDummyIndex () ;
	  char c = newDummyIndex () ;
	  char d = newDummyIndex () ;
	  char e = newDummyIndex () ;
	  char f = newDummyIndex () ;
	  char g = newDummyIndex () ;
	  char h1 = newDummyIndex () ;

	  POLYNOME p1 = newG (a, b, h) ;
	  POLYNOME p2 = newScalar (1.0, h) ;
	  p2->tt.denom[0] = 1 ;
	  p2->tt.denom[1] = 1 ;
	  p2->tt.denom[2] = 1 ;
	  p2->tt.mm[0][0] = a ;
	  p2->tt.mm[0][1] = b ;
	  printf ("============ k^ab/k^2(k+p)^2(k+p+q)^2\n") ;
	  showPol (p2) ;
	  p2 = dimIntegral (p2) ;
	  showPol (p2) ;
	  POLYNOME p3 = polProduct (p1, p2, h) ;
	  showPol (p3) ;
	  p3 = expand (p3) ;
	  showPol (p3) ;
 
	  p1 = newG (a, b, h) ;
	  p1->tt.g[2] = c ;
	  p1->tt.g[3] = d ;
	  p2 = newScalar (1.0, h) ;
	  p2->tt.denom[0] = 1 ;
	  p2->tt.denom[1] = 1 ;
	  p2->tt.denom[2] = 2 ;
	  p2->tt.mm[0][0] = a ;
	  p2->tt.mm[0][1] = b ;
	  p2->tt.mm[0][2] = c ;
	  p2->tt.mm[0][3] = d ;
	  printf ("============ k^abcd/k^2(k+p)^2(k+p+q)^4\n") ;
	  showPol (p2) ;
	  p2 = dimIntegral (p2) ;
	  showPol (p2) ;
	  p3 = polProduct (p1, p2, h) ;
	  showPol (p3) ;
	  p3 = expand (p3) ;
	  showPol (p3) ;
 
	  p1 = newG (a, b, h) ;
	  p1->tt.g[2] = c ;
	  p1->tt.g[3] = d ;
	  p1->tt.g[4] = e ;
	  p1->tt.g[5] = f ;
	  p2 = newScalar (1.0, h) ;
	  p2->tt.denom[0] = 1 ;
	  p2->tt.denom[1] = 2 ;
	  p2->tt.denom[2] = 2 ;
	  p2->tt.mm[0][0] = a ;
	  p2->tt.mm[0][1] = b ;
	  p2->tt.mm[0][2] = c ;
	  p2->tt.mm[0][3] = d ;
	  p2->tt.mm[0][4] = e ;
	  p2->tt.mm[0][5] = f ;
	  printf ("============ k^abcdef/k^2(k+p)^4(k+p+q)^4\n") ;
	  showPol (p2) ;
	  p2 = dimIntegral (p2) ;
	  showPol (p2) ;
	  p3 = polProduct (p1, p2, h) ;
	  showPol (p3) ;
	  p3 = expand (p3) ;
	  showPol (p3) ;
 
	  p1 = newG (a, b, h) ;
	  p1->tt.g[2] = c ;
	  p1->tt.g[3] = d ;
	  if (1)
	    {
	      p1->tt.g[4] = e ;
	      p1->tt.g[5] = f ;
	      p1->tt.g[6] = g ;
	      p1->tt.g[7] = h1 ;
	    }
	  else
	    {
	      p1->tt.eps[0] = e ;
	      p1->tt.eps[1] = f ;
	      p1->tt.eps[2] = g ;
	      p1->tt.eps[3] = h1 ;
	    }
	  p2 = newScalar (1.0, h) ;
	  p2->tt.denom[0] = 2 ;
	  p2->tt.denom[1] = 2 ;
	  p2->tt.denom[2] = 2 ;
	  p2->tt.mm[0][0] = a ;
	  p2->tt.mm[0][1] = b ;
	  p2->tt.mm[0][2] = c ;
	  p2->tt.mm[0][3] = d ;
	  p2->tt.mm[0][4] = e ;
	  p2->tt.mm[0][5] = f ;
	  p2->tt.mm[0][6] = g ;
	  p2->tt.mm[0][7] = h1 ;
	  printf ("============ k^abcdefgh/k^4(k+p)^4(k+p+q)^4\n") ;
	  showPol (p2) ;
	  p2 = dimIntegral (p2) ;
	  showPol (p1) ;
	  showPol (p2) ;
	  p3 = polProduct (p1, p2, h) ;
	  showPol (p3) ;
	  p3 = expand (p3) ;
	  showPol (p3) ;
 
	  exit (0) ;
	}
	  
      /* projector tests */
      if (0)
	{
	  firstDummyIndex = 'a' ;
	  char a = newDummyIndex () ;
	  char b = newDummyIndex () ;
	  char c = newDummyIndex () ;
	  char d = newDummyIndex () ;
	  char e = newDummyIndex () ;
	  char f = newDummyIndex () ;
	  char g = newDummyIndex () ;
	  char h1 = newDummyIndex () ;
	  
	  char i = newDummyIndex () ;
	  char j = newDummyIndex () ;
	  int z ;

	  if (0)
	    for (z = -1 ; z < 3 ; z++)
	      {
		printf ("\n\n\n@@@@@@@@@ Projector tests z = %d\n", z) ;
		POLYNOME p1 = newAG (a, b, c, d, z,h) ;
		POLYNOME p2 = newAG (c, d, e, f, z,h) ;
		POLYNOME p3 = newAG (e, f, g, h1, z,h) ;
		POLYNOME p4 = newAG (g, h1, i, j, z,h) ;
		POLYNOME ppp[] = {p1, p2, p3, p4, 0} ;
		POLYNOME pp = polMultiProduct (h,ppp) ;
		showPol (pp) ;
		pp = expand (pp) ;
		showPol (pp) ;
		pp = bbCleanUp (pp, a, b, i, j) ;
		showPol (pp) ;
	      }
	  
	  if (1)  /* is the propagator the inverse f the lagrangian */
	      {
		printf ("\n\n\n@@@@@@@@@ Projector tests z = %d\n", z) ;
		POLYNOME p1 = prop_BB_B (a, b, c, d, 99,h) ;
		POLYNOME p2 = prop_BB_B (c, d, e, f, 0,h) ;
		POLYNOME p3 = prop_BB_B (e, f, g, h1, 99,h) ;
		POLYNOME p4 = newScalar (1,h) ; /* compensate the values in the current propagators */
		POLYNOME ppp[] = {p1, p2, p3, p4, 0} ;
		showPol (p1) ;
		showPol (p2) ;
		showPol (p3) ;
		POLYNOME pp = polMultiProduct (h,ppp) ;
		showPol (pp) ;
		pp = expand (pp) ;
		showPol (pp) ;
		pp = bbCleanUp (pp, a, b, g, h1) ;
		showPol (pp) ;
	      }
	  exit (0) ;
	  
	  for (z = -1 ; z < 3 ; z++)
	    {
	      printf ("\n\n\n@@@@@@@@@ Projector tests z = %d\n", z) ;
	      POLYNOME p1 = newAG (a, b, c, d, z,h) ;
	      POLYNOME p2 = newScalar (1,h) ;
	      tcpy (p2->tt.sigma, "ecd") ;
	      POLYNOME pp = polProduct (p1, p2,h) ;
	      showPol (pp) ;
	      pp = expand (pp) ;
	      showPol (pp) ;
	    }
	  for (z = -1 ; z < 3 ; z++)
	    {
	      printf ("\n\n\n@@@@@@@@@ Projector tests z = %d\n", z) ;
	      POLYNOME p1 = newAG (a, b, c, d, z,h) ;
	      POLYNOME p2 = newScalar (1,h) ;
	      tcpy (p2->tt.sigB, "cd") ;
	      POLYNOME pp = polProduct (p1, p2,h) ;
	      showPol (pp) ;
	      pp = expand (pp) ;
	      showPol (pp) ;
	    }

	  
	  for (z = -1 ; z < 3 ; z++)
	    {
	      printf ("\n\n\n@@@@@@@@@ Projector tests A_psiB_psi B under z = %d\n", z) ;
	      POLYNOME p1 = newAG (a, b, c, d, z,h) ;
	      POLYNOME p2 = newScalar (1,h) ;
	      tcpy (p2->tt.sigma, "abfefcd") ;
	      POLYNOME pp = polProduct (p1, p2,h) ;
	      showPol (pp) ;
	      pp = expand (pp) ;
	      showPol (pp) ;
	    }

	  for (z = -1 ; z < 3 ; z++)
	    {
	      printf ("\n\n\n@@@@@@@@@ Projector tests A_psiB_psi B under z = %d\n", z) ;
	      POLYNOME p1 = newAG (a, b, e, h1, z,h) ;
	      POLYNOME p2 = newAG (c, d, g, h1, -z,h) ;
	      POLYNOME p3 = newScalar (1,h) ;
	      tcpy (p3->tt.sigma, "abgfecd") ;
	      POLYNOME ppp[] = {p1, p2, p3, 0} ;
	      POLYNOME pp = polMultiProduct (h,ppp) ;
	      showPol (pp) ;
	      pp = expand (pp) ;
	      showPol (pp) ;
	    }

	  
	  printf ("\n\n\n@@@@@@@@@ Projector tests DONE\n") ;

	  if (1) /* verify the eps eps contractions */
	    {
	      POLYNOME pp, ppp[12] ;
	      
	      ppp[0] = newEpsilon(a,c,f,b,h) ;
	      ppp[1] = newEpsilon(a,b,f,e,h) ;
	      ppp[2] = 0 ;
	      pp = polMultiProduct (h,ppp) ;
	      showPol (pp) ;
	      pp = expand (pp)  ;
	      showPol (pp) ;
	      pp = expand (pp)  ;
	      showPol (pp) ;
	      exit (0) ;
	    }
	  
	  if (0) /* verify the eps eps contractions */
	    {
	      POLYNOME pp = newEpsilon(a,b,i,f,h) ;
	      tcpy(pp->tt.eps,"abefcdei") ;
	      showPol (pp) ;
	      pp = expand (pp)  ;
	      showPol (pp) ;
	      exit (0) ;
	      
	      pp = newEpsilon(a,b,i,f,h) ;
	      tcpy(pp->tt.eps,"abhfcdhi") ;
	      showPol (pp) ;
	      pp = expand (pp)  ;
	      showPol (pp) ;
	      
	      pp = newEpsilon(a,b,i,f,h) ;
	      tcpy(pp->tt.g,"eh") ;
	      tcpy(pp->tt.eps,"abefcdhi") ;
	      showPol (pp) ;
	      pp = expand (pp)  ;
	      showPol (pp) ;
	      
	      pp = newEpsilon(a,b,i,f,h) ;
	      tcpy(pp->tt.eps,"abefcdhe") ;
	      showPol (pp) ;
	      pp = expand (pp)  ;
	      showPol (pp) ;
	      
	      
	      pp = newEpsilon(a,b,i,f,h) ;
	      tcpy(pp->tt.g,"ei") ;
	      tcpy(pp->tt.eps,"abefcdhi") ;
	      showPol (pp) ;
	      pp = expand (pp)  ;
	      showPol (pp) ;
	      
	      pp = newEpsilon(a,b,i,f,h) ;
	      tcpy(pp->tt.g,"eh") ;
	      tcpy(pp->tt.mm[1],"if") ;
	      tcpy(pp->tt.eps,"abefcdhi") ;
	      showPol (pp) ;
	      pp = expand (pp)  ;
	      showPol (pp) ;
	      
	      pp = newEpsilon(a,b,i,f,h) ;
	      tcpy(pp->tt.g,"ei") ;
	      tcpy(pp->tt.mm[1],"hf") ;
	      tcpy(pp->tt.eps,"abefcdhi") ;
	      showPol (pp) ;
	      pp = expand (pp)  ;
	      showPol (pp) ;
	      
	      pp = newEpsilon(a,b,i,f,h) ;
	      tcpy(pp->tt.mm[1],"if") ;
	      tcpy(pp->tt.eps,"abefcdei") ;
	      showPol (pp) ;
	      pp = expand (pp)  ;
	      showPol (pp) ;
	      
	      pp = newEpsilon(a,b,i,f,h) ;
	      tcpy(pp->tt.mm[1],"hf") ;
	      tcpy(pp->tt.eps,"abefcdhe") ;
	      showPol (pp) ;
	      pp = expand (pp)  ;
	      showPol (pp) ;
	      
	    }
	  exit (0) ;
	}

      if (0) /* verify some Pauli contractions */
	{
	  firstDummyIndex = 'a' ;
	  if (0) /* verify the Pauli trace */
	    {
	      POLYNOME pp = newSigma('a',h) ;
	      tcpy (pp->tt.sigma,"abcd") ;
	      showPol (pp) ;
	      pp = pauliTrace(pp, 0) ;
	      showPol (pp) ;
	      pp = newSigma('a',h) ;
	      tcpy (pp->tt.sigma,"abcdef") ;
	      showPol (pp) ;
	      pp = pauliTrace(pp, 0) ;
	      showPol (pp) ;
	      pp = newSigma('a',h) ;
	      tcpy (pp->tt.sigma,"abcdefgh") ;
	      showPol (pp) ;
	      pp = pauliTrace(pp, 0) ;
	      showPol (pp) ;
	      exit (0) ;
	    }
	  
	  if (0)
	    {
	      printf ("g_ab s_abcdef)\n") ;
	      POLYNOME p1 = newSigma ('e',h) ;
	      tcpy (p1->tt.sigma,"abcdef") ;
	      POLYNOME p2 = newG ('a','b',h) ;
	      POLYNOME ppp[] = {p1,p2,0} ;
	      POLYNOME pp = polMultiProduct (h,ppp) ;
	      pp = expand (pp) ;
	      showPol (pp) ;
	      printf ("Tr (g_ab s_abcdef)\n") ;
	      pp = pauliTrace(pp, 0) ;
	      pp = expand (pp) ;
	      showPol (pp) ;

	      printf ("Tr (s_abcdef)\n") ;
	      p1 = pauliTrace (p1, 0) ;
	      showPol (p1) ;
	      POLYNOME ppp2[] = {p1,p2,0} ;
	      pp = polMultiProduct (h,ppp2) ;
	      printf ("g_cd Tr(s_abcdef)\n") ;
	      pp = expand (pp) ;
	      showPol (pp) ;

	      exit (0) ;
	    }

	  if (0)
	    {
	      printf ("(g_cd s_abcdefgh)\n") ;
	      POLYNOME p1 = newSigma ('a',h) ;
	      tcpy (p1->tt.sigma,"abcdefgh") ;
	      POLYNOME p2 = newG ('a','b',h) ;
	      POLYNOME ppp[] = {p1,p2,0} ;
	      POLYNOME pp = polMultiProduct (h,ppp) ;
	      pp = expand (pp) ;
	      showPol (pp) ;
	      printf ("Tr (g_cd s_abcdefgh)\n") ;
	      pp = pauliTrace (pp, 0) ;
	      pp = expand (pp) ;
	      showPol (pp) ;

	      printf ("\n\nTr (s_abcdefgh)\n") ;
	      p1 = pauliTrace (p1, 0) ;
	      showPol (p1) ;
	      p2 = newG ('a','b',h) ;
	      p2->tt.z *= -1 ;
	      POLYNOME pp2, pp3, ppp2[] = {p1,p2,0} ;
	      pp2 = polMultiProduct (h,ppp2) ;
	      printf ("g_cd Tr(s_abcdefgh)\n") ;
	      pp2 = expand (pp2) ;
	      showPol (pp2) ;
	      pp2 = pauliTrace (pp2, 0) ;
	      printf ("pp2\n") ;
	      showPol (pp2) ;
	      printf ("pp\n") ;
	      showPol (pp) ;
	      pp3 = polSum (pp, pp2,h) ;
	      pp3 = expand (pp3) ;
	      printf ("pp\n") ;
	      showPol (pp) ;
	      printf ("pp2\n") ;
	      showPol (pp2) ;
	      printf ("pp3\n") ;
	      showPol (pp3) ;
	    }

	  if (0)
	    {
	      printf ("\n\n-g_ac g_ac\n") ;
	      POLYNOME p1 = newScalar (1,h) ;
	      POLYNOME p2 = newScalar (1,h) ;
	      tcpy (p1->tt.g,"ac") ;
	      tcpy (p2->tt.g,"acbdefgh") ;
	      POLYNOME pp = polProduct (p1, p2,h) ;
	      showPol (pp) ;
	      pp = expand (pp) ;
	      showPol (pp) ;
	    }
	  
	  if (0)
	    {
	      printf ("\n\ns_abca\n") ;
	      POLYNOME p1 = newScalar (1,h) ;
	      tcpy (p1->tt.sigma,"abcd") ;
	      POLYNOME p2 = newScalar (-1,h) ;
	      tcpy (p2->tt.sigma,"acbd") ;
	      POLYNOME pp = polSum (p1, p2,h) ;
	      showPol (pp) ;
	      pp = expand (pp) ;
	      showPol (pp) ;
	      p1 = newScalar (1,h) ;
	      tcpy (p1->tt.g,"ad") ;
	      pp = polProduct (p1, pp,h) ;
		  showPol (pp) ;
		  pp = expand (pp) ;
		  showPol (pp) ;
	    }
      
	  exit (0) ;
	}

      /* verify the B propagator */
      if (0)
	{
	    char a = newDummyIndex () ;
	    char b = newDummyIndex () ;
	    char c = newDummyIndex () ;
	    char d = newDummyIndex () ;
	    char e = newDummyIndex () ;
	    char f = newDummyIndex () ;
	    char g = newDummyIndex () ;
	    char h1 = newDummyIndex () ;
	    char i = newDummyIndex () ;
	    char j = newDummyIndex () ;
	    char k = newDummyIndex () ;

	    printf ("Simplification of the propagator\n") ;
	    POLYNOME p1 = newAG (a,b,e,f,1,h) ;
	    POLYNOME p2 = newEpsilon (e,f,g,h1,h) ;
	    p2->tt.mm[0][0] = g ;
	    POLYNOME p3 = newEpsilon (i,h1,j,k,h) ;
	    p3->tt.mm[0][0] = i ;
	    POLYNOME p4 = newAG (j,k,c,d, -1,h) ;
	    POLYNOME ppp [] =  {p1,p2,p3,p4,0} ;
	    POLYNOME pp = polMultiProduct (h,ppp) ;
	    showPol (pp) ;
	    pp = expand (pp) ;
	    showPol (pp) ;
	    pp = reduceIndices (pp) ;
	    showPol (pp) ;
	    pp = expand (pp) ;
	    showPol (pp) ;
	    
	    exit (0) ;
	}
      
      /* Fermion propagator */
      
      if (0) 
	{
	  printf ("\n\n\n@@@@@@@@@ Fermion propagator scalar under\n") ;
	  Z2_PsiL__H_Psi ("######### Fermion propagator, Scalar under\n") ;
		  
	  firstDummyIndex = 'a' ;
	  Z2_PsiL__A_Psi ("######### Fermion propagator, Vector under\n") ;
	  	  
	  firstDummyIndex = 'a' ;
	  Z2_PsiL__B_Psi ("######### Fermion propagator, Tensor under\n") ;
	  	  
	  printf ("\n\n\n@@@@@@@@@ Fermion propagator DONE\n\n") ;
	  exit (0) ;
	}

      /* pure gauge theory, coupling of the Vector to the Fermion in the presence of scalar/vector/tensor under */
      
      if (0) 
	{
	  printf ("\n\n\n@@@@@@@@@ Classic Ward identity : A_PsiB_Psi scalar under\n") ;
	  firstDummyIndex = 'a' ;
	  if (1) Z3_A_PsiL_PsiLB__Hunder () ; 
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ Classic Ward identity : A_PsiB_Psi scalar over\n") ;
	  if (1) Z3_A_PsiL_PsiLB__Hover () ; 

	  printf ("\n\n\n@@@@@@@@@ Classic Ward identity : A_PsiB_Psi Vector under\n") ;
	  firstDummyIndex = 'a' ;
	  if (1) Z3_A_PsiL_PsiLB__Aunder () ;
	  printf ("\n\n\n@@@@@@@@@ Classic Ward identity : A_PsiB_Psi Vector over\n") ;
	  if (1) Z3_A_PsiL_PsiLB__Aover () ;
	     
	  printf ("\n\n\n@@@@@@@@@ Classic Ward identity : A_PsiB_Psi tensor under\n") ;
	  firstDummyIndex = 'a' ;
	  if (1) Z3_A_PsiL_PsiLB__Bunder () ;
	  printf ("\n\n\n@@@@@@@@@ Classic Ward identity : A_PsiB_Psi tensor over\n") ;
	  firstDummyIndex = 'a' ;
	  if (1) Z3_A_PsiL_PsiLB__Bover () ;   /* BUG: expand loops forever */
	  
	  printf ("\n\n\n@@@@@@@@@ Classic Ward identity : A_PsiB_Psi DONE\n") ; 


	  exit (0) ;
	}

      if (0)
	{
	  printf ("\n\n\n@@@@@@@@@ Classic Ward identity counted on ghosts : A_cB_c\n") ;
	  firstDummyIndex = 'a' ;
	  if (1) Z2_AA__loopA ("######### Vector propagator, Vector loop, null in su(1/1) \n") ;
	  firstDummyIndex = 'a' ;
	  if (1) Z2_AA__loopGhost ("######### Vector propagator, Ghost loop, null in su(1/1) \n") ;
	  firstDummyIndex = 'a' ;
	  if (1) Z2_ghost__Aunder ("######### Ghost propagator, Vector under\n") ;
	  firstDummyIndex = 'a' ;
	  printf ("Z2 done\n") ;
	  if (1) Z3_A_c_cB__Aunder () ;
	  firstDummyIndex = 'a' ;
	  if (1) Z3_A_c_cB__Aover () ;
	  exit (0) ;
	}

      if (0)
	{
	  printf ("\n\n\n@@@@@@@@@ Classic Ward identity counted on scalars : A_cB_c\n") ;
	  firstDummyIndex = 'a' ;
	  if (1) Z2_HH__Aunder ("######### Scalar propagator, Vector under\n") ;
	  firstDummyIndex = 'a' ;
	  printf ("Z2 done\n") ;
	  if (1) Z3_A_H_HB__Aunder () ;
	  firstDummyIndex = 'a' ;
	  if (1) Z3_A_H_HB__Aover () ;
	  exit (0) ;
	}
      
      /* Boson propagators Fermion loops*/
      if (1)
	{
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ Boson propagators, Fermion loops */\n") ;

	  if (1) Z2_HH__loopPsi ("######### Scalar propagator, Fermion loop\n") ;
	  if (1) Z2_AA__loopPsi ("######### Vector propagator, Fermion loop\n") ;
	  if (1) Z2_BB__loopPsi ("######### Tensor propagator, Fermion loop\n") ; 

	  printf ("\n\n\n@@@@@@@@@ Boson propagators Fermion loops DONE\n") ;
	  exit (0) ;
	}
      
      /* Boson propagators Boson loops*/
      if (0)
	{
	  POLYNOME p1 = newScalar (1,h) ;
	  POLYNOME p2 = newScalar (1,h) ;
	  POLYNOME pp = newScalar (1,h) ;

	  /* Adkp Befp cdkn efmn   Ap-cn  Bp-mn   */
	  tcpy (p1->tt.eps, "adkpbefpcdknefmn") ;
	  tcpy (pp->tt.eps, "adkpcdknbefpefmn") ;
	  if (1) tcpy (p1->tt.eps, "adkpcdkn") ;
	  if (1) tcpy (p2->tt.eps, "befpefmn") ;
	  if (0)
	    {
	      showPol(p1) ;
	      p1 = expand (p1) ;
	      showPol (p1) ;
	      showPol(p2) ;
	      p2 = expand (p2) ;
	      showPol (p2) ;
	      POLYNOME p3 = polProduct (p1,p2,h) ;
	      showPol(p3) ;
	      p3 = expand (p3) ;
	      showPol (p3) ;
	    }
	  if (0) tcpy (pp->tt.eps, "akmnbckn") ;
	  showPol(pp) ;
	  pp = expand (pp) ;
	  showPol (pp) ;
	  
	  exit (0) ;
	}
      
      if (0)
	{
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ Boson propagators, Boson loops */\n") ;

	  if (1)
	    {
	      if (1) Z2_AA__loopH ("######### Vector propagator, Complex Scalar loop, -1/3 null in su(1/1) \n") ;
	      if (1) Z2_AA__loopGhost ("######### Vector propagator, Ghost loop, -1/6 null in su(1/1) \n") ;
	      if (1) Z2_AA__loopA ("######### Vector propagator, Vector loop, -3/2null in su(1/1) \n") ;
	      if (1) Z2_AA__loopB ("######### Vector propagator, Tensor loop, ?, null in su(1/1) \n") ;
	      if (0) Z2_AA__loopHB ("######### Vector propagator, Tensor-Scalar 0, loop, expect 0\n") ;
	      exit (0) ;
	    }
	  
	  if (1)
	    {
	      if (0) Z2_HH__Aunder ("######### Scalar propagator, Vector-under, 0 in SU(1/1)\n") ;
	      if (1) Z2_HH__loopAB ("######### Scalar propagator, Vector-Tensor loop\n") ;
	    }
	  
	  if (1)
	    {

	      if (0) Z2_BB__Aunder ("######### Tensor propagator, Vector-under, 0 in SU(1/1)\n") ;
	      if (1) Z2_BB__loopAH ("######### Tensor propagator, NEW Vector-Scalar loop\n") ;
	    }
	  printf ("\n\n\n@@@@@@@@@ Boson propagators Boson loops DONE\n") ;
	  exit (0) ;
	}
      
      /* coupling of the scalar to the Fermions, influenced by the scalar/vector/tensor */
      if (0)
	{
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ New Ward identity  H_PsiB_Psi Aunder vertex\n") ;
	  if (1) Z3_H_PsiR_PsiLB__Aunder () ;  
	  
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ New Ward identity  H_PsiB_Psi HAB vertex\n") ;
	  if (1) Z3_H_PsiR_PsiLB__HAB () ;  
	  firstDummyIndex = 'a' ; 
	  printf ("\n\n\n@@@@@@@@@ New Ward identity  H_PsiB_Psi HBA vertex\n") ; 
	  if (1) Z3_H_PsiR_PsiLB__HBA () ; 
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ New Ward identity  H_PsiB_Psi HHA vertex\n") ; 
	  if (1) Z3_H_PsiR_PsiLB__HHA () ; 
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ New Ward identity  H_PsiB_Psi HAH vertex\n") ; 
	  if (1) Z3_H_PsiR_PsiLB__HAH () ; 

	  printf ("\n\n\n@@@@@@@@@ New H-PsiB-Psi Ward identity  DONE\n") ; 


	  exit (0) ;
	}

      /* coupling of the vector to the Fermions, influenced by the scalar/vector/tensor */
      if (0)
	{
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ Classic Ward identity  A_PsiB_Psi Hunder vertex\n") ;
	  if (1) Z3_A_PsiR_PsiLB__Hunder () ;  
	  
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ Classic Ward identity  A_PsiB_Psi Aunder vertex\n") ;
	  if (1) Z3_A_PsiR_PsiLB__Aunder () ;  
	  
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ Classic Ward identity  A_PsiB_Psi Bunder vertex\n") ;
	  if (1) Z3_A_PsiR_PsiLB__Bunder () ;  
	  
#ifdef JUNK
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ New Ward identity  A_PsiB_Psi HAH vertex\n") ;
	  if (1) Z3_A_PsiR_PsiLB__HAH () ;  
	  firstDummyIndex = 'a' ; 
	  printf ("\n\n\n@@@@@@@@@ New Ward identity  A_PsiB_Psi AAA vertex\n") ; 
	  if (1) Z3_H_PsiR_PsiLB__AAA () ; 
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ New Ward identity  H_PsiB_Psi HHA vertex\n") ; 
	  if (1) Z3_H_PsiR_PsiLB__HHA () ; 
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ New Ward identity  H_PsiB_Psi HAH vertex\n") ; 
	  if (1) Z3_H_PsiR_PsiLB__HAH () ; 

	  printf ("\n\n\n@@@@@@@@@ New H-PsiB-Psi Ward identity  DONE\n") ; 

#endif
	  exit (0) ;
	}

      /* coupling of the tensor to the Fermions, influenced by the scalar/vector/tensor */
      if (0)
	{
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ New Ward identity  B_PsiB_Psi BAH vertex\n") ;
	  if (1) Z3_B_PsiR_PsiLB__BAH () ;  
	  
	  firstDummyIndex = 'a' ; 
	  printf ("\n\n\n@@@@@@@@@ New Ward identity  B_PsiB_Psi BHA vertex\n") ; 
	  if (1) Z3_B_PsiR_PsiLB__BHA () ; 

	  firstDummyIndex = 'a' ; 
	  printf ("\n\n\n@@@@@@@@@ New Ward identity  B_PsiB_Psi Aunder vertex\n") ;
	  if (1) Z3_B_PsiR_PsiLB__Aunder () ;  

	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ New Ward identity  B_PsiB_Psi BBA vertex\n") ; 
	  if (1) Z3_B_PsiR_PsiLB__BBA () ; 

	  printf ("\n\n\n@@@@@@@@@ New B-PsiB-Psi Ward identity  DONE\n") ; 


	  exit (0) ;
	}

      /* scalar/vector/tensor vertex, Boson loop */
      if (0)
	{
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ New Ward identity  A_H_BB boson loop\n") ;
	  if (1) Z3_A_H_BB__loopABH () ;  
	  printf ("\n\n\n@@@@@@@@@ New A_H_BB Boson loop  DONE\n") ; 


	  exit (0) ;
	}
      /* vector interactions with the scalar-vector-tensor in the presence of a Fermion loop */
      if (0)
	{
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ Vector-Boson vertex, Fermion loops */\n") ;

	  firstDummyIndex = 'a' ;
	  if (1) Z3_AHH__loopPsiL ("######### Vector-Scalar-Scalar, Fermion loop\n") ;
	  firstDummyIndex = 'a' ;
	  if (1) Z3_AAA__loopPsiL ("######### Vector-Vector-Vector, Fermion loop\n") ;
	  firstDummyIndex = 'a' ;
	  if (0) Z3_ABB__loopPsiL ("######### Vector-Tensor-Tensor, Left Fermion loop\n") ;
	  firstDummyIndex = 'a' ;
	  if (0) Z3_ABB__loopPsiR ("######### Vector-Tensor-Tensor, Right Fermion loop\n") ;
	  firstDummyIndex = 'a' ;
	  if (1) Z3_ABH__loopPsiL ("######### Scalar_Vector-Tensor, Fermion loop\n") ;

	  printf ("\n\n\n@@@@@@@@@ Boson propagators Fermion loops DONE\n") ;
	  exit (0) ;
	}
      /* triple vector interactions with boson loops */
      if (0)
	{
	  firstDummyIndex = 'a' ;
	  printf ("\n\n\n@@@@@@@@@ Triple Vector vertex, Boson loops */\n") ;

	  firstDummyIndex = 'a' ;
	  if (1) Z3_AAA__loopPsiL ("######### Triple Vector, Fermion loop\n") ;
	  firstDummyIndex = 'a' ;
	  if (1) Z3_AAA__loopH ("######### Triple Vector, Scalar loop\n") ;
	  firstDummyIndex = 'a' ;
	  if (1) Z3_AAA__loopGhost ("######### Triple Vector, Ghost loop\n") ;
	  firstDummyIndex = 'a' ;
	  if (1) Z3_AAA__loopA ("######### Triple Vector, vector loop\n") ;
	  firstDummyIndex = 'a' ;
	  if (1) Z3_AAA__loopB ("######### Triple Vector, tensor loop\n") ;

	  exit (0) ;
	}
      


      if (0)
	{
	  printf ("\n\n\n@@@@@@@@@ epsilon tests \n") ;
	  
	  POLYNOME p1 = newScalar (1,h) ;

	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdaefg") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;

	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdabed") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;

	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdabcd") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;
	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdabdc") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;
	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdacbd") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;
	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "dbcadacb") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;
	  
	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdabce") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;
	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdabed") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;
	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdaecd") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;
	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdbecd") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;
	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdbedc") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;
	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdbced") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;
	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdbcde") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;
	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdacde") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;


	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdabef") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;
	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdaefb") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;
	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdaecf") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;

	  p1 = newScalar (1,h) ; tcpy (p1->tt.eps, "abcdaefg") ;
	  showPol (p1) ; p1 = expand (p1) ; showPol (p1) ;

	  exit (0) ;	  
	}

      if (0)
        {
	  printf("\n### check that skew pairs of pauli are self dual \n") ;

	  char a = newDummyIndex () ;
	  char b = newDummyIndex () ;
	  char c = newDummyIndex () ;
	  char d = newDummyIndex () ;
	  char e = newDummyIndex () ;
	  char f = newDummyIndex () ;

	  POLYNOME r1 = newAG (a,b,c,d, 1, h) ;
	  POLYNOME r2 = newAG (a,b,c,d, -1, h) ;
	  POLYNOME r3 = newAG (c,d,e,f, 1, h) ;
	  POLYNOME r4 = newAG (c,d,e,f, -1, h) ;
	  POLYNOME s1 = newSigma (c, h) ;
	  s1->tt.sigma[1] = d ;
	  POLYNOME s2 = newSigB (c, h) ;
	  POLYNOME p5 = 0 ;

	  s2->tt.sigB[1] = d ;
	  
	  printf ("\nP^2\n") ;
	  p5 = polProduct (r1,r3,h) ;
	  showPol (p5) ;
	  p5 = expand (p5) ;
	  showPol (p5) ;
	  
	  printf ("\nPM\n") ;
	  p5 = polProduct (r1,r4,h) ;
	  showPol (p5) ;
	  p5 = expand (p5) ;
	  showPol (p5) ;
	  
	  printf ("\nM^2\n") ;
	  p5 = polProduct (r2,r4,h) ;
	  showPol (p5) ;
	  p5 = expand (p5) ;
	  showPol (p5) ;
	  
	  printf ("\nMP\n") ;
	  p5 = polProduct (r2,r3,h) ;
	  showPol (p5) ;
	  p5 = expand (p5) ;
	  showPol (p5) ;
	  
	  printf ("\nP.Sig.SigB\n") ;
	  p5 = polProduct (r1,s1,h) ;
	  showPol (p5) ;
	  p5 = expand (p5) ;
	  showPol (p5) ;
	  
	  printf ("\nM.Sig.SigB\n") ;
	  p5 = polProduct (r2,s1,h) ;
	  showPol (p5) ;
	  p5 = expand (p5) ;
	  showPol (p5) ;
	  
	  printf ("\nP.SigB.Sig\n") ;
	  p5 = polProduct (r1,s2,h) ;
	  showPol (p5) ;
	  p5 = expand (p5) ;
	  showPol (p5) ;
	  
	  printf ("\nM.SigB.Sig\n") ;
	  p5 = polProduct (r2,s2,h) ;
	  showPol (p5) ;
	  p5 = expand (p5) ;
	  showPol (p5) ;
	  exit (0) ;
	}

      if (0)
	{
	  printf ("\n\n\n@@@@@@@@@ Boson propagators, contraction tests */\n") ;
	  if (1)
	    {
	      firstDummyIndex = 'a' ;
	      char a = newDummyIndex () ;
	      char b = newDummyIndex () ;
	      char c = newDummyIndex () ;
	      char d = newDummyIndex () ;
	      char e = newDummyIndex () ;
	      char f = newDummyIndex () ;
	      
	      POLYNOME p1 = newAG (a,b,c,b,1,h) ;
	      POLYNOME p2 = newAG (c,d,a,d,-1,h) ;
	      POLYNOME ppp[] = {p1, p2, 0} ;
	      POLYNOME p3 = polMultiProduct (h,ppp) ;
	      showPol (p3) ;
	      p3 = expand (p3) ;
	      showPol (p3) ;
	      
	      p1 = newAG (a,b,c,d,1,h) ;
	      p2 = newAG (d,a,b,c,-1,h) ;
	      POLYNOME ppp2[] = {p1, p2, 0} ;
	      p3 = polMultiProduct (h,ppp2) ;
	      showPol (p3) ;
	      p3 = expand (p3) ;
	      showPol (p3) ;
	      
	      p3 = newAG (a,b,a,b,1,h) ;
	      showPol (p3) ;
	      p3 = expand (p3) ;
	      showPol (p3) ;
	      
	      p3 = newAG (a,b,c,b,0,h) ;
	      showPol (p3) ;
	      p3 = expand (p3) ;
	      showPol (p3) ;
	      
	      p1 = newPQR (0, a,h) ;
	      p2 = newPQR (0, b,h) ;
	      p3 = newPQR (0, c,h) ;
	      POLYNOME p4 = newPQR (0, d,h) ;
	      p1->tt.denom[0] = 2 ;
	      p1->tt.denom[1] = 1 ;
	      POLYNOME ppp3[] = {p1, p2, p3, p4, 0} ;
	      p3 = polMultiProduct (h,ppp3) ;
	      showPol (p3) ;
	      p3 = dimIntegral (p3) ;
	      p3 = expand (p3) ;
	      showPol (p3) ;
	      
	      p1 = newAG (a,c,e,f,1,h) ;
	      p2 = newAG (b,d,e,f,-1,h) ;
	      POLYNOME ppp4[] = {p1, p3, p2, 0} ;
	      p3 = polMultiProduct (h,ppp4) ;
	      showPol (p3) ;
	      p3 = expand (p3) ;
	      p3 = squareMomentaCleanUp (p3) ;
	      showPol (p3) ;
	       
	    }
	  exit (0) ;
	}

      if (jacobi)
	{
	  printf ("\n\n\n@@@@@@@@@ 4-tensor vertex, hoping to find zero to match the *F F lagrangian 2024_03\n") ;
	  if (1)
	    {
	      firstDummyIndex = 'a' ;
	      /* external tensor indices (circular trace or B^2 C^2 */
	      char a = newDummyIndex () ;
	      char b = newDummyIndex () ;
	      char c = newDummyIndex () ;
	      char d = newDummyIndex () ;
	      /* pauli tensor indices */
	      char m = newDummyIndex () ;
	      char n = newDummyIndex () ;
	      char o = newDummyIndex () ;
	      char p = newDummyIndex () ;
	      char q = newDummyIndex () ;
	      char r = newDummyIndex () ;
	      char s = newDummyIndex () ;
	      char t = newDummyIndex () ;
	      POLYNOME B1 = newAG (a,b,m,n,-1,h) ;
	      POLYNOME C1 = newAG (b,c,o,p,1,h) ;
	      POLYNOME B2 = newAG (c,d,q,r,-1,h) ;
	      POLYNOME C2 = newAG (d,a,s,t,1,h) ;

	      POLYNOME BB1 = newAG (a,b,m,n,-1,h) ;
	      POLYNOME CC1 = newAG (c,d,o,p,1,h) ;
	      POLYNOME BB2 = newAG (a,b,q,r,-1,h) ;
	      POLYNOME CC2 = newAG (c,d,s,t,1,h) ;

	      /* internal lines */
	      char i = newDummyIndex () ;
	      char j = newDummyIndex () ;

	      /* Loop with crossing internal lines */
	      printf ("\n.... BCBC circular trace\n");
	      POLYNOME loop = newSigma (m, h) ;
	      int ii, *ip, idx[] = {m,n,i,o,p,j,q,r,i,s,t,j,0} ;
	      for (ii = 0, ip = idx ; *ip ; ii++, ip++)
		loop->tt.sigma[ii] = *ip ;
	      showPol (loop) ;
	      
	      POLYNOME ppp[] = {B1,C1,B2,C2,loop, 0} ;
	      POLYNOME p1 = polMultiProduct (h,ppp) ;
	      showPol (p1) ;
	      p1 = expand (p1) ;
	      showPol (p1) ;
	      p1 = pauliTrace (p1, 0) ;
	      showPol (p1) ;
	      p1 = expand (p1) ;
	      showPol (p1) ;

	      printf ("\n.... BCBC : BB CC trace\n");
	      POLYNOME loop2 = newSigma (m, h) ;
	      
	      int idx2[] = {m,n,i,o,p,j,q,r,i,s,t,j,0} ;
	      for (ii = 0, ip = idx2 ; *ip ; ii++, ip++)
		loop2->tt.sigma[ii] = *ip ;
	      showPol (loop2) ;
	      
	      POLYNOME ppp2[] = {BB1,CC1,BB2,CC2,loop2, 0} ;
	      p1 = polMultiProduct (h,ppp2) ;
	      showPol (p1) ;
	      p1 = expand (p1) ;
	      showPol (p1) ;
	      p1 = pauliTrace (p1, 0) ;
	      showPol (p1) ;
	      p1 = expand (p1) ;
	      showPol (p1) ;
	      
	    }
	  exit (0) ;
	}
    
    }
#ifdef JUNK
  
  /* triangles which only exist in the non-Abelian (non su(1/1) case */
  if (1) Z3_AHH_loopAAH ("######### Vector-scalar-scalar vertex, Scalar_below-vector-vector loop\n") ;
  if (1) Z3_AHH_loopAAB ("######### Vector-scalar-scalar vertex, Tensor_below-vector-vector loop\n") ;
  if (1) Z3_AHH_loopHHA ("######### Vector-scalar-scalar vertex, Vector_below-scalar-scalar loop\n") ;
  if (1) Z3_AHH_loopBBA ("######### Vector-scalar-scalar vertex, Vector_below-tensor-tensor loop\n") ;
  if (1) Z3_AHH_loopBHA ("######### Vector-scalar-scalar vertex, Vector_below-tensor-scalar loop\n") ;
  if (1) Z3_AHH_loopHBA ("######### Vector-scalar-scalar vertex, Vector_below-scalar-tensor loop\n") ;
  
  
  
#endif
  
  /* superalgebra Jacobi indentities */
  
  if (jacobi)
    {
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

    }
  return 0 ;
}

  
