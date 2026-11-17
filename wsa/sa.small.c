/*
 * analysis if short RNA

 * 1: scan and count all 12-mers in 100k first reads
 * 2: extend the most frequent words, the winner is the adaptor
 * 3: dna2dna -rightClipOn
 * 4: create a table:   dna_word:atgc  BOOL:ref_Id/inRun  char*:ID/run_name
 * 5: sort the table, compress:  dna_word ref_Id(extended via new_1234) run/count//run/count//run/count...
 * 6: sort by total count, map the top of the table via a hash table by scanning the genome

 */

#include "sa.h"

/* Find the adaptor. in short reads we expect to read it in most sequences */
BOOL saCountKmers (Array adaptor, Array dnas, int run, int KK)
{
  AC_HANDLE h = ac_new_handle () ;
  if (KK > 12) messcrash ("saCountKmers can only count up to KK=12") ;

  Array kMers = arrayHandleCreate (0x1 << (2*KK), int, h) ;
  BOOL ok = FALSE ;
  int ii, iMax = arrayMax (dnas) ;
  unsigned int w = 0, wMax = 0, n = 0,  nMax = 0 ;
  char d2a[4], d2[256] = {0} ;
  d2[A_] = 0 ;
  d2[T_] = 3 ;
  d2[G_] = 1 ;
  d2[C_] = 2 ;

  d2a[0] = 'a' ;
  d2a[3] = 't' ;
  d2a[1] = 'g' ;
  d2a[2] = 'c' ;
  
  /* scan all sequences */
  for (ii = 0 ; ii < iMax ; ii++)
    {
      Array dna = array (dnas, ii, Array) ;
      int j, jMax = dna ? arrayMax(dna) - KK : 0 ;

      if (jMax < 0)
	continue ;
      /* count the KK-mers */
      unsigned char *cp = arrp (dna, 0, unsigned char) ;
      for (j = 0 ; j < KK - 1 ; j++)
	w = (w << 2) | d2[(int)*cp] ;
      for ( ; j < jMax ; j++)
	{
	  w = (w << 2) | d2[(int)*cp] ;
	  n = array (kMers, w, int) + 1 ;
	  array (kMers, w, int) = n ;
	  if (n > nMax) { nMax = n ; wMax = w ; }
	}
    }

  arrayMax (adaptor) = 0 ;
  if (nMax > arrayMax (dnas)/2)
    {
      ok = TRUE ;			  
      /* construct the best word */
      char cc ;
      int nBest ;      

      w = wMax ;
      n = nMax = nBest = array (kMers, w, int) ;
      
      unsigned int w1 = w ;
      char buf[KK+1] ;
      buf[KK] = 0 ;
      for (int i = 0 ; i < KK ; i++)
	{ buf [KK - i - 1] = d2a[w1 & 0x3] ; w1 >>= 2 ; }
      char *vv = strnew (buf, h) ;
      unsigned int aw, tw, gw, cw ;
      int naw, ntw, ngw, ncw ;
      
      /* extend left */
      while (nBest > nMax/2)
	{
	  aw = tw = gw = cw = (w >> 2) ;
	  aw |= 0x0 << (2 * KK) ; naw = array (kMers, aw, int) ; 
	  tw |= 0x3 << (2 * KK) ; ntw = array (kMers, aw, int) ; 
	  gw |= 0x1 << (2 * KK) ; ngw = array (kMers, aw, int) ; 
	  cw |= 0x2 << (2 * KK) ; ncw = array (kMers, aw, int) ; 
	  int nBest = naw ; cc = 'a' ;
	  if (ntw > nBest) { nBest = ntw ; cc = 't' ; } 
	  if (ngw > nBest) { nBest = ngw ; cc = 'g' ; } 
	  if (ncw > nBest) { nBest = ncw ; cc = 'c' ; } 
	  if (nBest < nMax / 2) break ;
	  vv = hprintf (h, "%c%s", cc, vv) ;
	}
      
      /* extend right */
      w = wMax ;
      n = nMax = nBest = array (kMers, w, int) ;
      while (nBest > nMax/2)
	{
	  unsigned int aw, tw, gw, cw ;
	  int n1 = array (kMers, w, int) ;
	  
	  aw = tw = gw = cw = (w << 2) ;
	  aw |= 0x0 ; naw = array (kMers, aw, int) ; 
	  tw |= 0x3 ; ntw = array (kMers, aw, int) ; 
	  gw |= 0x1 ; ngw = array (kMers, aw, int) ; 
	  cw |= 0x2 ; ncw = array (kMers, aw, int) ; 
	  int nBest = naw ; cc = 'a' ;
	  if (ntw > nBest) { nBest = ntw ; cc = 't' ; } 
	  if (ngw > nBest) { nBest = ngw ; cc = 'g' ; } 
	  if (ncw > nBest) { nBest = ncw ; cc = 'c' ; } 
	  if (nBest < n1 / 2) break ;
	  vv = hprintf (h, "%s%c", vv, cc ) ;
	}
      
      fprintf (stderr, "Adaptor : %s %d\n", vv, nMax) ;
      int i ; char *cp ;
      for (i = 0, cp = vv ; *cp ; i++, cp++)
	array (adaptor, i, char) = *cp ;
      dnaEncodeArray (adaptor) ;
    }
  
  ac_free (h) ;
  return ok ;
} /* saCountKMers */

	  

	  
