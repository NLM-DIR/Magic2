#include    <stdio.h>

int firstDivisor (int nn)
{
#define NN 10000              /* max size of prime table */
  static int primeList [NN] ; /* lazy table of known primes */
  static int maxKnown = 0 ;
  int i, p, q ;               /* running integers */
  int testing = 1 ;           /* 1: show prime list when new number computed */

  if (!maxKnown) /* initialise */
    {
      memset (primeList, 0, NN * sizeof (int)) ;
      maxKnown = 1 ;
      primeList [0] = 2 ;
    }
  while (1) /* lazy fill of primeList */
    /* NB if p is our largest known prime, there
       is always another prime q < p*p
    */
    {
      p = primeList [maxKnown - 1] ;
      if (p*p > nn) 
	break ;
      for (q = p + 1 ; firstDivisor (q) ; q++) ;
      if (testing)
	printf ("%6d is prime\n", q) ;
      primeList [maxKnown++] = q ;
      if (maxKnown == NN)
	{ 
	  printf ("Overflow in firstDivisor nn=%d, maxKnown=%d q = %d\n", 
		  nn, maxKnown, q) ;
	  exit (1) ;
	}
    }
  /* testing */
  for (i = 0 ; i < maxKnown ; i++)
    {
      p = primeList [i] ;
      if (p*p > nn) 
	return 0 ; /* nn is prime, I return 0 for ease of use */
      if ( nn == p * (nn/p))
	return p ; /* not a prime */
    }
  return 0 ; /* nn is prime, I return 0 for ease of use */
}

int main (int argc, char *argv[])
{
  int k, p ;

  k = 8 ; p = firstDivisor (k) ;
  printf ( "first divisor of %d is %d\n\n", k, p) ;

  k = 301 ; p = firstDivisor (k) ;
  printf ( "first divisor of %d is %d\n\n", k, p) ;
  k = 1027 ; p = firstDivisor (k) ;
  printf ( "first divisor of %d is %d\n\n", k, p) ;
  exit (0) ;

}
/*hello from kaa */
