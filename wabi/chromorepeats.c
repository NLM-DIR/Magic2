/*  File: chromorepeats.c
 *  Author: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg. 2002
 *-------------------------------------------------------------------
 * This file is part of the ACEMBLY extension to 
        the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk
 *	Jean Thierry-Mieg (NCBI and CNRS France) mieg@ncbi.nlm.nih.gov
 *
 * Description:
   Find all repeates of length 1 to w in a whole genome
 * Exported functions:
 * HISTORY:
 *   several algo to count the repeats in the complete genome
 *
 *   chromoCleanRepeats is a compact form
 *   chromoRepeats is the full fledged form
 *     these are based on jumping to next position repeating current
 *     working recursively from repeat length zero to max
 *   chromoRepeatSort is a totally different algo
 *   
 * Created: Nov 2002 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */
#define MAXWORD 100

#ifndef ACEMBLY

#include <stdio.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include "../wh/acedna.h"
#else

#include "acedb.h"
#include "bitset.h"
#include "../wabi/chromorepeats.h"
#include "../whooks/sysclass.h"
#include "../whooks/systags.h"
#include "dna.h"
#include "a.h"
#include "cdna.h"
#include "bs.h"
#include "lex.h"
#include "../whooks/classes.h"
#include "../whooks/tags.h"

/***************************************************************/
/***************************************************************/

static Array  chromoRepeatsTabulate (int NN, int w, unsigned int *zz)
{ 
  unsigned int j, n, *zp, i, i0 ;
  BitSet bb = bitSetCreate (NN, 0) ; /* resets bb to zero */
  Array rr = arrayCreate (NN, unsigned char) ;
  /* store the number of repeats of len wMax starting at each position */
  
  array (rr, NN - 1, int) = 0 ;
  /* for each start point
     if not repeated -> *rr = 0
     if repeated, circulate and count, then circulate, set *rr and flag
  */
  for (i = i0 = 0, zp = zz ; i < NN - w ; i++, zp = zz + i)
    {
      if (bit(bb,i))
	continue ;
      if (!*zp)
	{ array (rr, i, int) = 0 ;  continue ; }
      /* count the repeats */ 
      n = 0 ; 
      while (*zp) { n++ ; j = *zp ; zp = zz + j ; } 
      /* flag the repeats */
      zp = zz + i ;  array (rr, i, int) = n ; bitSet (bb,i) ;
      while (*zp) { j = *zp ; array (rr, j, unsigned char) = n > 255 ? 255 : n ; bitSet (bb,j) ; zp = zz + j ; }
    }
  bitSetDestroy (bb) ; 
  return rr ;
} /* chromoRepeatsTabulate */

#endif /* ACEMBLY */

/***************************************************************/
/* count the dna repeats of length up to wMax */
/* theoretical expectation if dna was fully random */

static long chromoRepeatsTheory (long nn, int w)
{
  double p, q, z1, z2, z3, r ;
  
  z1 = exp (w * log((double)4.0)) ;
  p = 1.0/z1 ; q = 1.0 - p ;
  z2 = exp (nn * log (q)) ;
  z3 = nn * p * exp ((nn - 1) * log (q)) ;
  r = z1 * (1 - z2 - z3) ;
   
  return r >= 0 ? (long) (r + .49) : 0 ;
}

/* may only works on a 64 bits machine */

/***********/
/* most compact algorithm, actually 50% faster than the next one */
static void* chromoCleanRepeats (unsigned int *probes, int nProbe, char *dna, long NN, int wMax, int doTabulate)
{
  int w ;
  long int i, isNew, nRepeats, nNewRepeats, nComparisons, nJumps ;
  unsigned int *zz, *zp, *tp, z1, mask = 0x80000000 ;
  char cc ;
  void *rr = 0 ;
  time_t t0 = time (0)  ;

  zz = (unsigned int *) malloc (NN * sizeof(int)) ;
  /* NN is the length of the dna */
  /* zz[i] is the offset of the first repeat of [dna[i]..dna[i+w-1] */
  /* jump is the next position which may have a repeat */
  printf("\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n\n"
	 , "word length","# types","# repeats","# bounces", "# jumps", "Theory", "Time") ;
  nRepeats = nNewRepeats = 0 ;
  for (zp = zz, i = 1 ; i < NN ; i++)
    { 
      switch(dna[i])
	{ 
	case 'a': case 't': case 'g': case 'c': 
#ifdef ACEMBLY
	case A_: case T_: case G_: case C_:
#endif
  /* initialisation, all 'atgc' words of length zero are repeated at next 'atgc' position (zz) */
  /* all 'atgc' position must be tested for repeats (jump) */
	  *zp = i ; zp = zz + i ;
	  nRepeats++ ;
	  break ;
	default: /* immediately skip ambiguous letters like n or - */
	  zz [i] = 0 ;
	  break ;
	}
    }
  *zp = 0 ;  /* except the last one */
  if (nProbe > 1)
    {
      z1 = probes[nProbe-1] + 1 ;
      for (i = 0 ; i < nProbe - 1 ; i++)
	{
	  zz[probes[i]] = z1 ;
	  zz[probes[i] + 1] = 0 ;
	}
    }
  if (0) for (i = 0, zp = zz ; i < NN  ; i++, zp++)
    fprintf (stderr, "%ld\t%d\n", i, *zp) ;
  printf ("\t%d\t%d\t%ld\t%ld\t%ld\t%d\t%ds\n", 0, 4, nRepeats, NN, NN, 4, (int)(time(0) - t0)) ;
  /* recursion */
  /* nRepeats is the number of repeats of length w */
  /* nNewRepeats is the number of different repeats of length w */
  /* nComparisons is the number of comparison needed during this recursion */
  for (w = 1 ; nRepeats > 0 && w <= wMax ; w++)
    {
      nRepeats = nNewRepeats = nComparisons = nJumps = 0 ;
      isNew = 1 ;
      if (0) for (i = 0, zp = zz ; i < NN  ; i++, zp++)
	fprintf (stderr, "%ld\t%d\n", i, *zp) ;
      for (i = 0, zp = zz ; i < NN - w ; i++, zp++)
	{
	  /* search if  word w[i] of length w starting at i is repeated downwards */
	  z1 = *zp ;
	  if (!z1)      /* w[i] cannot be repeated if (w-1)[i] is not */
	    continue ;
	  if (w > 1 && !*(zp+1)) 
	    { *zp = 0 ; continue ; } /* w[i] cannot be repeated if (w-1)[i+1] is not */

	  cc = dna[i + w - 1 ] ;  /* just compare a single new letter */
	  if (z1 & mask)
	    { 
	      *zp = z1 = z1 & ~mask ; /* clean up the sign flag */
	      isNew = 0 ;   /* there exists a repeat upstream */
	    }
	  else 
	    isNew = 1 ;
	  tp = zz + z1 ;

	  while (1)
	    {
	      nComparisons++ ;
	      if (cc == dna[w - 1 +  z1])    /* found the next repeat */
		{ 
		  nRepeats++ ; /* count all repeats */ 
		  nNewRepeats += isNew ;
		  /* if (0) printf("    %d:%d\n", i,*tp) ;  if you want, show details */
		  *zp = z1 ;    /* i.e. zz[i] now points to *tp */
		  z1 = *tp ;
		  if (z1 > 0) *tp |= mask ;  /* flag the downstream repeat as not new */

		  break ;        /* done: found repeat */
		} 
	      z1 = *tp & ~mask ;
	      if (z1) 
		tp = zz + z1 ;  /* loop: try next candidate */
	      else
		{ *zp = 0 ; break ; }       /* done: no repeat */
	    }
	}
      printf ("\t%d\t%ld\t%ld\t%ld\t%ld\t%ld\t%ds\n", w, nNewRepeats, nRepeats, nComparisons, nJumps
	      , chromoRepeatsTheory(NN, w)
	      , (int)(time(0) - t0)) ;
    }

#ifdef ACEMBLY
  if (doTabulate && wMax > 0 && w == wMax+1)
    rr = (void *) chromoRepeatsTabulate (NN, w, zz) ;
#endif

  free (zz) ;
  return rr ;

} /* chromoCleanRepeats */

/***********/
/* we try now to work by half base, to have less bouncing */

static void* chromoHalfBaseRepeats (char *dna, int NN, int wMax, int doTabulate)
{
  int i, isNew, w, nRepeats, nNewRepeats, nComparisons, nJumps ;
  unsigned int *zz, *zp, *tp, z1, mask = 0x80000000 ;
  char cc ;
  void *rr = 0 ;
  time_t t0 = time (0)  ;

  zz = (unsigned int *) malloc (NN * sizeof(int)) ;
  /* NN is the length of the dna */
  /* zz[i] is the offset of the first repeat of [dna[i]..dna[i+w-1] */
  /* jump is the next position which may have a repeat */
  printf("\n%18s%18s%12s%12s%12s%8s%8s\n\n"
	 , "word length","# types","# repeats","# bounces", "# jumps", "Theory", "Time") ;
  nRepeats = nNewRepeats = 0 ;
  for (zp = zz, i = 1 ; i < NN ; i++)
    { 
      switch(dna[i])
	{ 
	case 'a': case 't': case 'g': case 'c': 
	case A_: case T_: case G_: case C_:
  /* initialisation, all 'atgc' words of length zero are repeated at next 'atgc' position (zz) */
  /* all 'atgc' position must be tested for repeats (jump) */
	  *zp = i ; zp = zz + i ;
	  nRepeats++ ;
	  break ;
	default: /* immediately skip ambiguous letters like n or - */
	  zz [i] = 0 ;
	  break ;
	}
    }
  *zp = 0 ;  /* except the last one */
  printf ("%18d%18d%12d%12d%12d%8d%8ds\n", 0, 4, nRepeats, NN, NN, 4, (int)(time(0) - t0)) ;
  /* recursion */
  /* nRepeats is the number of repeats of length w */
  /* nNewRepeats is the number of different repeats of length w */
  /* nComparisons is the number of comparison needed during this recursion */
  for (w = 1 ; nRepeats > 0 && w <= 2 * wMax ; w++)
    {
      nRepeats = nNewRepeats = nComparisons = nJumps = 0 ;
      isNew = 1 ;
      for (i = 0, zp = zz ; i < NN - (w+1)/2 ; i++, zp++)
	{
	  /* search if  word w[i] of length w starting at i is repeated downwards */
	  z1 = *zp ;
	  if (!z1)      /* w[i] cannot be repeated if (w-1)[i] is not */
	    continue ;
	  if (w > 2 && !*(zp+1)) 
	    { *zp = 0 ; continue ; } /* w[i] cannot be repeated if (w-1)[i+1] is not */

	  cc = dna[i + (w+1)/2 - 1 ] ; /* next letter */
	  if (w & 0x1)
	    switch (cc)
	      {
	      case A_: case G_: cc = A_ | G_ ; break ;
	      case T_: case C_: cc = T_ | C_ ; break ;
	      default: cc = N_ ; break ;
	      }

	  if (z1 & ~mask)
	    { 
	      *zp = z1 = z1 & ~mask ; /* clean up the sign flag */
	      isNew = 0 ;   /* there exists a repeat upstream */
	    }
	  else 
	    isNew = 1 ;
	  tp = zz + z1 ;

	  while (1)
	    {
	      nComparisons++ ;
	      if (cc & dna[(w+1)/2 - 1 +  z1])    /* found the next half repeat */
		{ 
		  nRepeats++ ; /* count all repeats */ 
		  nNewRepeats += isNew ;
		  /* if (0) printf("    %d:%d\n", i,*tp) ;  if you want, show details */
		  *zp = z1 ;    /* i.e. zz[i] now points to *tp */
		  z1 = *tp & ~mask ;

		  break ;        /* done: found repeat */
		} 
	      z1 = *tp & ~mask ; ;
	      if (z1) 
		tp = zz + z1 ;  /* loop: try next candidate */
	      else
		{ *zp = 0 ; break ; }       /* done: no repeat */
	    }
	}
      printf ("%18d%18d%12d%12d%12d%8ld%8ds\n", w, nNewRepeats, nRepeats, nComparisons, nJumps
	      , chromoRepeatsTheory(NN, (w+1)/2)
	      , (int)(time(0) - t0)) ;
    }

#ifdef ACEMBLY
  if (doTabulate && wMax > 0 && w == wMax+1)
    rr = (void *) chromoRepeatsTabulate (NN, w, zz) ;
#endif

  free (zz) ;
  return rr ;

} /* chromoHalfBaseRepeats */

#ifdef ACEMBLY

/***************************************************************/
 /***************************************************************/
/***********/
/* in addition to the algorithm above, we include a jump table
 * which accelerates the system by jumping useless areas
 * the jump could also be used to jump from the end of each est into the genome
 * which of course means we should reinitialise a search at the
 * top of each est
 *
 * actually by correct initialisation of zz, we may obtain the same effect
 * without the need for jump[]
 */

static Array chromoRepeats (char *dna, int NN, int wMax, BOOL doTabulate)
{
  unsigned int i, i0, w, nRepeats, nNewRepeats, nComparisons, nJumps, *zz, *jump, *zp, *tp, *jp ;
  char cc ;
  Array rr = 0 ;
  BitSet bb = bitSetCreate (NN, 0) ;
  time_t t0 = time (0)  ;

  zz = messalloc (NN * sizeof(unsigned int)) ;
  jump = messalloc (NN * sizeof(unsigned int)) ;
  /* NN is the length of the dna */
  /* zz[i] is the offset of the first repeat of [dna[i]..dna[i+w-1] */
  /* jump is the next position which may have a repeat */
  printf("\n%18s%18s%12s%12s%12s%8s%8s\n\n"
	 , "word length","# types","# repeats","# bounces", "# jumps", "Theory", "Time") ;
  nRepeats = nNewRepeats = 0 ;
  for (jp = jump, zp = zz, i = 1 ; i < NN ; i++)
    { 
      switch(dna[i])
	{ 
	case 'a': case 't': case 'g': case 'c': 
	case A_: case T_: case G_: case C_:
  /* initialisation, all 'atgc' words of length zero are repeated at next 'atgc' position (zz) */
  /* all 'atgc' position must be tested for repeats (jump) */
	  *jp = i ; jp = jump + i ;
	  *zp = i ; zp = zz + i ;
	  nRepeats++ ;
	  if (nNewRepeats < 4 && !bit(bb, dna[i]))
	    { nNewRepeats++ ; bitSet (bb, dna[i]) ; }
	  break ;
	default: /* immediatly skip ambiguous letters like n or - */
	  jump[i] = 0 ; zz [i] = 0 ;
	  break ;
	}
    }
  *zp = 0 ; *jp = NN ; /* except the last one */
  printf ("%18d%18d%12d%12d%12d%8d%8ds\n", 0, nNewRepeats, nRepeats, NN, NN, 4, (int)(time(0) - t0)) ;

  /* recursion */
  /* nRepeats is the number of repeats of length w */
  /* nNewRepeats is the number of different repeats of length w */
  /* nComparisons is the number of comparison needed during this recursion */
  for (w = 1 ; nRepeats > 0 && w <= wMax ; w++)
    {
      nRepeats = nNewRepeats = nComparisons = nJumps = 0 ;
      bb = bitSetReCreate (bb, NN) ; /* resets bb to zero */
      for (i = i0 = 0, zp = zz ; i < NN - w ; i = jump[i], zp = zz + i)
	{
	  nJumps++ ;
	  /* search if  word w[i] of length w starting at i is repeated downwards */
	  if (w > 1 && !*(zp+1)) 
	    {
	      *zp = 0 ; /* w[i] cannot be repeated if (w-1)[i+1] is not */
	    }
	  if (*zp)    /* w[i] cannot be repeated if (w-1)[i] is not */
	    { tp = zp ; cc = dna[i + w - 1 ] ; } /* just compare a single new letter */
	  while (*zp)
	    {
	      nComparisons++ ;
	      if (cc == dna[*tp + w - 1])    /* found the next repeat */
		{ 
		  nRepeats++ ; /* count all repeats */ 
		  if (0) printf("    %d:%d\n", i,*tp) ; /* if you want, show details */
		  *zp = *tp ;    /* i.e. zz[i] now points to *tp */
		  if (!bit (bb, i)) 
		    nNewRepeats++ ;     /* no upstream repeat, since not flagged previously */
		  bitSet (bb, *tp) ;  /* flag the downstream repeat */
		  break ;        /* done: found repeat */
		} 
	      else if (*tp) 
		tp = zz + *tp ;  /* loop: try next candidate */
	      if (!*tp)
		*zp = 0 ;        /* done: no repeat */
	    }
	  if (!*zp)             /* we will not need to test position i at stage w+1 */
	    jump[i0] = jump[i] ;  /* but will jump direct from i0 to successor of i */ 
	  else 
	    i0 = i ;    /* establish i as the successor of i0 at next round */
	}
      printf ("%18d%18d%12d%12d%12d%8ld%8ds\n", w, nNewRepeats, nRepeats, nComparisons, nJumps
	      , chromoRepeatsTheory(NN, w)
	      , (int)(time(0) - t0)) ;
    }

  messfree (jump) ;
  bitSetDestroy (bb) ; 
  
  if (doTabulate && wMax > 0 && w == wMax+1)
    rr = chromoRepeatsTabulate (NN, w, zz) ;

  messfree (zz) ; 
  return rr ;
} /* chromoRepeats */

/***************************************************************/
/***********/
/* in addition to the algorithm above, 
 * we include a prefix of length p
 * We create a table giving the successive positions of all
 * the occurences of that prefix in the genome, so the length
 * of this table is NP =~ length(genome)/4^p, which can be small
 * except for a few saturated prefixes
 * we then do everything as above but working of this shorter table
 * so the memory requirement is now 
 *      genome * 2 bits + (genome/4^p)*(3 integers, zz, jump, offset)
 */
typedef struct prStuct { unsigned int z, jump, x ;} PRS ;
static Array chromoPrefixRepeats (char *dna, int NN, int wMax, BOOL doTabulate, char *prefix)
{
  int i, i0, x, NP, w, z1 = 0, nRepeats, nNewRepeats, nComparisons, nJumps ;
  int prefixLength = strlen (prefix) ;
  char cc = 0 ;
  Array rr = 0 ;
  BitSet bb = 0 ;
  PRS *pp0, *pp = 0, *pp1, *prs ;
  Array aa = 0 ;
  time_t t0 = time(0) ;

  /* NN is the length of the dna */
  /* NP is the number of repeats of the prefix in the dna */
  /* zz[i] is the offset of the first repeat of [dna[p[i].x]..dna[p[i].x+w-1] */
  /* jump is the next position which may have a repeat */

  /* step 0: initialise the table for the given prefix */
  aa = arrayCreate (NN/( 1 << (2*prefixLength - 1)), PRS) ; /* elastic self reallocated table */

  printf("Prefix %s\n\n%18s%18s%12s%12s%12s%8s%8s\n\n"
	 , prefix, "word length","# types","# repeats","# bounces","# jumps", "Theory","Time") ;
  nRepeats = nNewRepeats = 0 ;
  for (i = 0, x = 0 ; x < NN - prefixLength + 1 ; x++)
    { 
      if (strncmp (dna + i, prefix, prefixLength))
	continue ;
      pp = arrayp (aa, i++, PRS) ;
      pp->x = x ; /* offset of the jth occurence of the prefix in the genome */
      pp->z = i ; /* all 'prefix' words are repeated at next 'prefix' position */
      pp->jump = i ; /* there is no hole */
    }
  pp->z = 0 ; pp->jump = NN ; /* the last occurence is the last one ! */
  nRepeats = NP = i ;
  prs = arrayp (aa, 0, PRS) ; /* now that aa has reached its full size, its address is final */

  printf ("%18d%18d%12d%12d%12d%8d%8ds\n", 0, 1, nRepeats, NN, 1, 1, (int)(time(0) - t0)) ;

  /* step 2: compute recursivelly the repeats of length 1, 2,   ... w */
  /* recursion */
  /* nRepeats is the number of repeats of length w */
  /* nNewRepeats is the number of different repeats of length w */
  /* nComparisons is the number of comparison needed during this recursion */

  for (w = prefixLength  + 1 ; nRepeats > 0 && w <= wMax ; w++)
    {
      nRepeats = nNewRepeats = nComparisons = nJumps = 0 ;
      bb = bitSetReCreate (bb, NP) ; /* resets bb to zero */
      for (i = 0, pp = pp0 = prs ; i < NP ; i = pp->jump, pp = prs + i)
	{
	  /* search if  word w[i] of length w starting at i is repeated downwards */
	  /* w[i] cannot be repeated if (w-1)[i+1] is not, no loger applies */
	  pp1 = pp ;
	  nJumps++ ;
	  if (pp->z)    /* w[i] cannot be repeated if (w-1)[i] is not */
	    { z1 = pp->z ; pp1 = prs + z1 ; cc = dna[pp->x + w - 1 ] ; } /* just compare a single new letter */
	  while (pp1->z)
	    {
	      nComparisons++ ;
	      if (cc == dna[pp1->x + w - 1])    /* found the next repeat */
		{ 
		  nRepeats++ ; /* count all repeats */ 
		  if (0) printf("    %d:%d\n", i, pp1->x) ; /* if you want, show details */
		  pp->z = z1 ;    /* i.e. pp now points to pp1 */
		  if (!bit (bb, i)) 
		    nNewRepeats++ ;     /* no upstream repeat, since not flagged previously */
		  bitSet (bb, z1) ;  /* flag the downstream repeat */
		  break ;        /* done: found repeat */
		} 
	      else if (pp1->z) 
		{ z1 = pp1->z ; pp1 = prs + z1 ; }  /* loop: try next candidate */
	      if (!pp1->z)
		pp->z = 0 ;        /* done: no repeat */
	    }
	  if (!pp->z)             /* we will not need to test position i at stage w+1 */
	    pp0->jump = pp->jump ;  /* but will jump direct from pp0 to successor of pp */ 
	  else 
	    pp0 = pp ;    /* establish pp as the successor of pp0 at next round */
	}
      printf ("%18d%18d%12d%12d%12d%8ld%8ds\n", w, nNewRepeats, nRepeats, nComparisons, nJumps
	      , chromoRepeatsTheory(NN, w - prefixLength)
	      , (int)(time(0) - t0)) ;
    }
  if (doTabulate && wMax > 0 && w == wMax+1)
    { unsigned int j, n ;
    
    /* store the number of repeats of len wMax starting at each position */
      rr = arrayCreate (NP, int) ;
      array (rr, NP - 1, int) = 0 ;
      bb = bitSetReCreate (bb, NP) ; /* resets bb to zero */
      /* for each start point
	 if not repeated -> *rr = 0
	 if repeated, circulate and count, then circulate, set *rr and flag
      */
      for (i = i0 = 0, pp = prs ; i < NP ; i++, pp++)
	{
	  if (bit(bb,i))
	    continue ;
	  if (!pp->z)
	    { array (rr, i, int) = 0 ;  continue ; }
	  /* count the repeats */ 
	  n = 0 ; 
	  while (pp->z) { n++ ; j = pp->z ; pp = prs + j ; }
	  /* flag the repeats */
	  pp = prs + i ;  array (rr, i, int) = n ; bitSet (bb,i) ;
	  while (pp->z) { j = pp->z ; array (rr, j, int) = n ; bitSet (bb,j) ; pp = prs + j ; }
	}
    }
  bitSetDestroy (bb) ;
  arrayDestroy (aa) ;

  return rr ;
} /* chromoPrefixRepeats */

/***************************************************************/
/***************************************************************/
/* implement the HM Martinez algorithm 
 * which has to do with sorting all words 
 */

static char *myDna = 0 ; 

static int chromoRepeatOrder (const void *va, const void *vb)
{
  const unsigned int ia = *(const unsigned int*)va, ib = *(const unsigned int*)vb ;
  char *cp, *cq ;
  int nn = 30 ;
  cp = myDna + ia ;
  cq = myDna + ib ;
  while (nn--)
    {
      if (*cp != *cq) return *cp - *cq ;
      cp++ ; cq++ ;
    }
  return 0 ;  
}

static void chromoRepeatSort (char *dna, long NN, int wMax)
{
  unsigned int i, w, *zp, *tp ;
  unsigned int ia, ib ;
  char *cp, *cq ;
  Array zz ;

  myDna = dna ; /* used in sort */

  zz = arrayCreate (NN, unsigned int) ;
  array (zz, NN-1, unsigned int) = 0 ;
  tp = (unsigned int*) messalloc (wMax * sizeof (unsigned int)) ;
  printf("\n%18s%18s%12s%12s%12s\n\n", "word length","# types","# repeats","# bounces","Theory") ;
  
  for (zp = arrp(zz, 0, unsigned int), i = 0 ; i < NN - wMax ; zp++, i++)
    *zp = i ;
  arraySort (zz, chromoRepeatOrder) ;
  for (i = 0 ; i + 1 < arrayMax(zz) ; i++)
    {
      ia = arr(zz, i, unsigned int) ;
      ib = arr(zz, i+1, unsigned int) ;
      cp = myDna + ia ;
      cq = myDna + ib ;
      for (w = 0 ; *cp == *cq && w < wMax ; cp++, cq++, w++)
	*(tp+w) += 1 ;
    }
  for (w = 0 ; w < wMax ; w++)
    printf("%18d%18d%12d%12d%12ld\n", w+1, 0, *(tp+w), 0,  chromoRepeatsTheory(NN, w+1)) ;
  messfree (tp) ;
  arrayDestroy (zz) ;
}

/***************************************************************/
/***************************************************************/

typedef struct offStruct { KEY key ; unsigned int len, offset ; } OFS ;

int  chromoRepeatsKeyset (KEYSET ks)
{
  mysize_t ii ; KEY key ;
  unsigned long int total = 0 ;
  Stack s = 0 ;
  Array dna, ofs = 0 ;
  char *cp ;
  OFS *up ;
  
  if (!ks || !keySetMax (ks))
    return 0 ;
  ofs = arrayCreate (keySetMax(ks), OFS) ;
  for (ii = total = 0 ; ii < keySetMax(ks) ; ii++)
    {
      key = keySet (ks, ii) ;
      if ((dna = dnaGet (key)))
	{
	  if (!s)
	    s = stackCreate (arrayMax(dna)) ;
	  /* dnaDecodeArray (dna) ; */
	  up = arrayp (ofs, ii, OFS) ;
	  up->key = key ; 
	  up->len = arrayMax(dna) ;
	  total += up->len ;
	  if (ii > 0) up->offset = (up - 1)->offset + (up - 1)->len ;
	  catBinary (s, arrp(dna, 0, char), up->len) ; 
	  cp = stackText(s,0) ;
	  printf("Adding %s (%d) total length %lu\n", name(key),up->len, total) ;
	  arrayDestroy (dna) ;
	}
    }
  cp = stackText(s,0) ;
  printf("total length %ld\n", total) ;
  if (1) 
    {
      OBJ Chrom, Cosmid = 0 ;
      int j, x1, x2 ;
      BSunit *uu ;
      Array units = arrayCreate (1000, BSunit) ;
      KEY rkey, gKey;
      Array aa = 0, aa2 = 0 ;
      int rLength = 18 ;
      Array rr = chromoCleanRepeats (0, 0, cp, strlen(cp), rLength, TRUE) ;

      if (rr)
	{
	  for (ii = 0 ; ii < arrayMax (ofs) ; ii++)
	    {
	      up = arrp (ofs, ii, OFS) ;
	      aa = arrayReCreate (aa, up->len, unsigned char) ;
	      array (aa, up->len - 1, unsigned char) = 0 ;
	      memcpy (aa->base, arrp(rr, up->offset, int), sizeof(unsigned char) * up->len) ;
	      if (keyFindTag (up->key, str2tag("Genomic")) || 
		  class (up->key) == _VmRNA)
		{ /* store directly in this object */
		  lexaddkey (messprintf ("%s.i100", name(up->key)), &rkey, _VOligoRepeat) ;
		  arrayStore (rkey, aa, "i") ;
		}
	      else if (keyFindTag (up->key, _Subsequence))
		{ /* store in its genomic descendant */
		  Chrom = bsCreate (up->key) ;
		  
		  if (Chrom && 
		      bsGetArray (Chrom, _Subsequence, units, 3))
		    for (j = 0 ; j < arrayMax(units) ; j+= 3)
		      {
			uu = arrp (units, j, BSunit) ;
			x1 = uu[1].i ; x2 = uu[2].i ;
			gKey = uu[0].k ;
			if (x1 < x2 &&
			    keyFindTag (gKey, str2tag("Genomic")))
			  {
			    lexaddkey (messprintf ("%s.i100", name(gKey)), &rkey, _VOligoRepeat) ;
			     aa2 = arrayReCreate (aa2, x2 - x1 + 1, unsigned char) ;
			     array (aa2, (x2 - x1 + 1) - 1, unsigned char) = 0 ;
			     memcpy (aa2->base, arrp(aa, x1 - 1, unsigned char), sizeof(unsigned char) * (x2 - x1)) ; 
			     arrayStore (rkey, aa2, "c") ;
			     if ((Cosmid = bsUpdate (gKey)))
			       {
				 bsAddKey (Cosmid, str2tag("OligoRepeat"), rkey) ;
				 bsAddData (Cosmid, _bsRight, _Int, &rLength) ;
				 bsSave (Cosmid) ;
			       }
			  }
		      }
		  
		  bsDestroy (Chrom) ;
		}		
	    }
	  messfree (rr) ;
	}
      arrayDestroy (units) ;
      arrayDestroy (aa) ;
      arrayDestroy (aa2) ;
    }
  else
    {
      char prefix[] = { 'a', 't', 'g', 'c', 'g', 'c', 0} ; /* some prefix */
      
      printf ("Clean repeats start %s\n", timeShowNow()) ;
      chromoCleanRepeats (0, 0, cp, total, MAXWORD, FALSE) ; 
      printf (" Clean repeats done %s\n", timeShowNow()) ;
      if (0)
	{printf ("Half repeats start %s\n", timeShowNow()) ;
	chromoHalfBaseRepeats (cp, total, MAXWORD, FALSE) ; 
	printf (" Half repeats done %s\n", timeShowNow()) ;
	printf ("Fast repeats start %s\n", timeShowNow()) ;
	chromoRepeats (cp, total, MAXWORD, FALSE) ; 
	printf (" Fast repeats done %s\n", timeShowNow()) ;
	}
      if (0) chromoPrefixRepeats (cp, total, 100, FALSE, prefix) ;
    }
  if (0) chromoRepeatSort (cp, total, 30) ;
  stackDestroy (s) ;
  arrayDestroy (ofs) ;
  return 1 ;
}

#else

/* compile and run as follows
 *
gcc -O4 -o crep chromorepeats.c -lm
crep < chrom3.fasta # a fasta file with a single sequence

* this is 15% faster than when run inside acedb
*/

#define NMAX 300000000
#define NPMAX 2000
int main (int argc, char **argv)
{
  char *cp, *buff ;
  unsigned int i, nn = NMAX, n1, n2 = 0, nProbe=0 ;
  int total = 0 ;
  unsigned int probes[NPMAX] ;
  unsigned int probeLength[NPMAX] ;

  cp = buff = malloc (NMAX) ;
  while (fgets (cp, nn, stdin))
    {
      if (*cp == '>')
	{
	  /* regiter the probe and
	     add a null base to propagate correctly bacwards
	   * the clause w[i] cannot be repeated if (w-1)[i+1] is not
	   */
	  if (cp > buff)
	    {
	      probes[nProbe++] = n2 - 1 ;
	      *cp++ = 0 ; nn-- ; n2++ ;
	    }
	  continue ;
	}
      n1 = strlen (cp) - 1  ; /* remove the \n */
      cp += n1 ;
      nn -= n1 ;
      n2 += n1 ;
      if (n2 > NMAX)
	{
	  fprintf (stderr,
		   "// There were %d bp > %d,\n"
		   "//  please edit the source code or parse less dna"
		   , n2, NMAX) ;
	  exit (1) ;
	}
      if (nProbe > NPMAX)
	{
	  fprintf (stderr,
		   "// There were %d probes > %d,\n"
		   "//  please edit the source code or parse less dna"
		   , nProbe, NPMAX) ;
	  exit (1) ;
	}
    }
  for (i = 0, cp = buff ; i < n2 ; i++, cp++)
    *cp = tolower(*cp) ;
    
  printf ("Clean repeats start\n") ;
  chromoCleanRepeats (probes, nProbe, buff, n2, probeLengthMax ? probeLengthMax : wMax, 0) ;
  printf (" Clean repeats done\n") ;

  free (buff) ;
}

#endif

