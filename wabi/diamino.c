/*  File: diamino.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg. 2001
 *-------------------------------------------------------------------
 * This file is part of the ACEMBLY extension to 
        the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
     Given an abstract set of letters (i.e. a set of vertex of an unoriented graph)
     with compatibility rules we construct all maximal words such are
     all the letters of a given word are mutually compatible.
     Maximal means that no word in the resulting list is contained in a longer
     word of that same list.

 * Exported functions: 
 ****   Array diaminoCreate (int lMax, void *vp, DCF isCompatible, DCF isConnected, int diamMax)
  Input:
     lMax        :: The number of letters,
     vp          :: a void* passed back to the compatibility function.
     isCompatible:: a function to test compatibility of a pair of letters,
     isConnected :: a function to test connectivity of a pair of letters,
     diamMax     :: if > 0, an absolute limit on number of words to create
  Result:
                an Array          of DIAMWORD
                (list of words)   (each word)  

 ****  void diaminoTest(void) : a  built in test
 ****  void diaminoDestroy (Array aa) : destroy the result structure

 * HISTORY:
 * Last edited: Dec  9 15:18 1998 (fw)
 * Created: Thu Dec  9 00:01:50 1997 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */

#define CHRONO

#include "cdna.h"
#include "cdnatr.h"
#include "bitset.h"
#include "diamino.h"

/***********************************************************************/
/**********************************************************************************/
/* Diamino algorithm (danielle) */
/* we try to construct all words of maximal length such that any pair of
   letter within the word is pair compatible

   init: we initiliaze and list all antilinks
   run: we factor the graph in connected components
        for each component we separate the most highly connected vertex
        and iterate
*/

typedef struct diamStruct 
{ int lMax ;           /* user supplied number of letters */
  void *vp ;           /* user supplied call back handle */
  DCF isCompatible;    /* user supplied compatibility function */
  DCF isConnected;     /* user supplied connectivity function */
  Array links ;        /* array of bitset to flag non-compatibilities */
  Array used ;         /* allocated here to avoid create/destroy in recursion */
  /* array contig ;       numbers the connected component of each letter */
  /* int ctgMax ;         number of connected components */
  /* Array diams ;        final results */     
  AC_HANDLE h ;     /* everything is allocated on that handle */
  int DIAMMAX ;
  void **magic ;
} DIAM ;

#define diamDestroy(_diam) (diamDoDestroy(_diam), diam=0)
static void diamShow (int lMax, Array diams, int target) ;
static void* diamMagic = 0 ;
static Array diamRun (DIAM *diam, BitSet mask, int level) ;

/**********************************************************************************/
/**********************************************************************************/

static DIAM *diamCreate (void *vp, int lMax, DCF isCompatible, DCF isConnected, int diamMax) 
{
  AC_HANDLE h = handleCreate () ;
  DIAM *diam = (DIAM*) halloc (sizeof (struct diamStruct), h) ;

  diam->h = h ;
  diam->lMax = lMax ;
  diam->vp = vp ;
  diam->isCompatible = isCompatible ;
  diam->isConnected = isConnected ;
  diam->used = arrayHandleCreate (lMax, KEY, h) ;
  array (diam->used, lMax - 1, KEY) = 0 ; /* make room */
  diam->DIAMMAX = diamMax > 0 ? diamMax : ACEDB_MAXINT  ;
  diam->magic = &diamMagic ;
  return diam ;
}

/**********************************************************************************/

static void diamDoDestroy (DIAM *diam)
{
  AC_HANDLE h = 0 ;
  
  if (diam && diam->h)
    { h = diam->h ; messfree (h) ; }
}

/**********************************************************************************/
/* construct antilinks */
static void diamLinks (DIAM* diam)
{
  int lMax = diam->lMax ;
  BitSet b1, b2 ;
  Array links ;
  int i1, i2 ;
  AC_HANDLE h = diam->h ;

  diam->links = links = arrayHandleCreate (lMax, BitSet, h) ;

  i1 = lMax ;
  while (i1--)
    {
      b1 = array (links, i1, BitSet) = bitSetCreate (lMax, h) ;
      bitUnSet (b1, lMax - 1) ; /* make room */
    }

  for (i1 = 0 ; i1 < lMax - 1 ; i1++)
    {
      b1 = arr (links, i1, BitSet) ;
      for (i2 = i1 + 1 ; i2 < lMax ; i2++)
	{ 
	  b2 = arr (links, i2, BitSet) ;
	  if (!(diam->isCompatible)(diam->vp, i1, i2))
	    {
	      bitSet (b1, i2) ;
	      bitSet (b2, i1) ;
	    }
	}
    }
}

/**********************************************************************************/
/* number the connected components */
static int diamContigs (int lMax, DIAM* diam, Array contig, BitSet mask, Array cardinal)
{
  BitSet b1 ;
  Array links = diam->links ;
  KEYSET used ;
  int i, i1, i2, c1, c2, cc, ctgMax ;

  used = diam->used = arrayReCreate (diam->used, lMax, KEY) ;
  for (i = 0 ; i < lMax ; i++)
    {  /* initialy, each letter is isolated in its own contig */
      arr (contig, i, int) = bit (mask, i) ? -1 : i  ;
    }

  for (i1 = 0 ; i1 < lMax - 1 ; i1++)
    if (!bit (mask, i1))
      {
	b1 = arr (links, i1, BitSet) ;
	for (i2 = i1 + 1 ; i2 < lMax ; i2++) 
	  if (!bit (mask, i2))	
	    if (bit (b1, i2))
	      {
		c1 = arr (contig, i1, int) ;
		c2 = arr (contig, i2, int) ;
		if (c1 > c2) { cc = c1 ; c1 = c2 ; c2 = cc ; }
		
		for (i = 0 ; i < lMax ; i++)
		  if (!bit (mask, i))
		    if (arr (contig, i, int) == c2)
		      arr (contig, i, int) = c1 ;
	      }
      }
  /* find occupied contigs */

  for (i = 0 ; i < lMax ; i++)
    if (!bit (mask, i))
      keySet (used, arr (contig, i, int)) = 1 ;

  /* attribute minimal number set */
  for (i1 = i2 = 0 ; i1 < lMax ; i1++)
    if (keySet(used, i1)) keySet (used, i1) = i2++ ;
  ctgMax = i2 ;

  /* renumber contigs */
   for (i = 0 ; i < lMax ; i++) 
     if (!bit (mask, i))
       {
	 i1 = arr (contig, i, int) ; /* in 2 lines because C does not guarantee order */
	 i2 = keySet (used, i1) ;
	 arr (contig, i, int) = i2 ;
	 arr (cardinal, i2, int)++ ;
       } 
  return ctgMax ;
}

/**********************************************************************************/
/**********************************************************************************/
/* component of tensor product: (ab)(cd) -> abcd)  */
static void diamDoProduct (int lMax, DIAMWORD *ma, DIAMWORD *mb, DIAMWORD *mc)
{
  BitSet ba = ma->bb, bb = mb->bb, bc = mc->bb ;
  int i  = lMax, p = mc->p ;

  while (i--)
    if (bit(ba,i) || bit(bb,i))
      { bitSet (bc, i) ; p++ ;}
  mc->p = p ;
  mc->discard = FALSE ;
}

/**********************************************************************************/
/* Tensor product  (a + bc)(d + e) -> ad + ae + bcd + bce */
static Array diamProduct (DIAM *diam, Array aa, Array bb)
{
  int ia, ib, jj, iaMax, ibMax ;
  Array product ;
  DIAMWORD *ma, *mb, *mc ;
  
  iaMax = aa ? arrayMax (aa) : 0 ;
  ibMax = bb ? arrayMax (bb) : 0 ;
  
  product = arrayCreate (iaMax * ibMax, DIAMWORD) ; 
  for (ia = 0, jj = 0 ; jj < diam->DIAMMAX && ia < iaMax ; ia++)
    { 
      ma = arrayp (aa, ia, DIAMWORD) ; 
      if (!ma->discard)
	for (ib = 0 ; jj < diam->DIAMMAX && ib < ibMax ; ib++)
	  {
	    mb = arrayp (bb, ib, DIAMWORD) ; 
	    if (!mb->discard)
	      {
		mc = arrayp (product, jj++, DIAMWORD) ;
		mc->p = 0 ;
		mc->bb = bitSetCreate (diam->lMax, 0) ;
		diamDoProduct (diam->lMax, ma, mb, mc) ;
	      }
	  }
    }
  return product ;
}

/**********************************************************************************/
/* Addition: (a + b) + (c + d) -> (a + b + c + d) */
static void diamSum (DIAM *diam, Array aa, Array bb)
{
  int ia, ib, ibMax ;
  DIAMWORD mb ;

  ia = arrayMax (aa) ;
  ibMax = arrayMax (bb) ;
  
  for (ib = 0 ; ib < ibMax ; ib++)
    { 
      mb = array (bb, ib, DIAMWORD) ; 
      if (ia < diam->DIAMMAX && !mb.discard)
	array (aa, ia++, DIAMWORD) = mb ;
      else
	bitSetDestroy (mb.bb) ;
    }
  arrayMax(bb) = 0 ;
}

/**********************************************************************************/

static void diamAddLetter (DIAM *diam, Array diams, int let)
{
  int nn = arrayMax(diams) ;
  DIAMWORD *mm ;

  if (nn == 0) /* create a void word */
    {
      mm = arrayp (diams, 0, DIAMWORD) ; 
      mm->p = 1 ;
      bitSetDestroy (mm->bb) ;
      mm->bb = bitSetCreate (diam->lMax, 0) ;
      bitUnSet (mm->bb, diam->lMax - 1) ;
      bitSet (mm->bb, let) ;
    }
  else
    while (nn--)
      {
	mm = arrayp (diams, nn, DIAMWORD) ;
	if (bit (mm->bb, let))
	  messcrash ("bit already set in diamAddLetter") ;
	bitSet (mm->bb, let) ; 
	(mm->p)++ ;
      }
}

/**********************************************************************************/
/* treat a connected graph with several vertex */
static void diamFilterDiams (DIAM *diam, Array diams, int ii)
{
  BitSet bb, ll = arr (diam->links, ii, BitSet) ;
  DIAMWORD *mm ;
  int i, nn = arrayMax(diams), lMax = diam->lMax ;
  BOOL ok ;

   while (nn--)
      {
	mm = arrayp (diams, nn, DIAMWORD) ;
	bb = mm->bb ;
	ok = FALSE  ;
	i = lMax ;
	while (i--)
	  if (bit (ll,i) && bit (bb,i))
	    ok = TRUE ;
	if (!ok)
	  mm->discard = TRUE ;
      }
}

/**********************************************************************************/
/* treat a connected graph with several vertex */
static Array diamRecurse (DIAM *diam, int ctg, Array contig, BitSet oldMask, int level)
{
  int lMax = diam->lMax ;
  BitSet mask = bitSetCreate (lMax, 0) ;
  int i, j, n ;
  int besti = 0 , bestn = 0 ;
  BitSet bb, ll ;
  Array diamsWith = 0 , diamsWithout = 0 ;
  BOOL debug = FALSE ;

  /* mask all but the relevant contig  */
  if (debug) printf ("diamRecurse(level %d):: using ", level) ;
  for (i = 0 ; i < lMax ; i++)
    if (bit (oldMask,i) || arr (contig, i, int) != ctg)
      bitSet (mask, i) ;
    else
      {
	bitUnSet (mask, i) ;
	if (debug) printf ("%d ",i) ;
      }
  if (debug) printf ("\n") ;
  /* count the number of links of each vertex, select best one */
  for (i = 0 ; i < lMax ; i++)
    if (!bit (mask, i)) 
      {
	bb = arr (diam->links, i, BitSet) ;
	for (n = j = 0 ; j < lMax ; j++)  
	  if (!bit (mask, j))
	    if (bit(bb, j)) n++ ;
	if (n > bestn) 
	  { bestn = n ; besti = i ; }
      }
  if (debug) printf ("highest node %d: %d links\n", besti, bestn) ;
  
  if (!bestn)
    messcrash ("diamRecurse: a connected component should always have at least one link") ;
  
  bitSet (mask, besti) ;
  /* first,  avoid besti, and use all other letters */
  diamsWithout = diamRun (diam, mask, level + 1) ; 
  if (debug) printf ("diamsWithout : %d words\n", arrayMax(diamsWithout)) ;
  diamFilterDiams (diam, diamsWithout, besti) ; /* each member must have a link to besti */

  /* then, use besti, and mask out all letters connected to it */
  ll = arr (diam->links, besti, BitSet) ;
  for (j = 0 ; j < lMax ; j++)  
    if (bit(ll, j)) 
      bitSet (mask, j) ;
  if (arrayMax(diamsWithout) < diam->DIAMMAX)
    {
      diamsWith = diamRun (diam, mask, level + 1) ;
      if (debug) printf ("diamsWith : %d words\n", arrayMax(diamsWith)) ;
      diamAddLetter (diam, diamsWith, besti) ;
      
      diamSum (diam, diamsWith, diamsWithout) ; /* adds into diamsWith */
      diaminoDestroy (diamsWithout) ;
    }
  else
    diamsWith = diamsWithout ;
  bitSetDestroy (mask) ;
  return diamsWith ;
}

/**********************************************************************************/
/* break the graph into connected commponents, and return product of individual results */
static Array diamEvaluate (DIAM *diam, int ctgMax, Array contig, BitSet mask, Array cardinal, int level) 
{
  Array oldDiams, diams  = arrayCreate (12, DIAMWORD), subDiams ;
  int ctg, j, lMax = diam->lMax ;
  BOOL debug = FALSE ;

  if (debug)  printf ("diamEval(%d)\n", level) ;
  for (ctg = 0 ; arrayMax(diams) < diam->DIAMMAX && ctg < ctgMax ; ctg++)
    switch (arr (cardinal, ctg, int))
      {
      case 0: 
	messcrash ("empty contig in diamEvaluate") ; break ;
      case 1: 
	for (j = 0 ; j < lMax && arr (contig, j, int) != ctg; j++) ;
	if (j >= lMax)
	  messcrash ("no member to ctg %d in diamEvaluate", ctg) ;
	diamAddLetter (diam, diams, j) ;
	if (debug)
	  {
	     printf ("adding letter %d contig %d\n", j, ctg) ;
	     diamShow (lMax, diams, 0) ;
	  }
	break ;
      default:
	subDiams = diamRecurse (diam, ctg, contig, mask, level) ;
	if (debug)
	  {
	     printf ("subContig %d\n", ctg) ;
	     diamShow (lMax, diams, 0) ;
	  }
	if (arrayMax(diams))
	  {
	    oldDiams = diams ;
	    diams = diamProduct (diam, oldDiams, subDiams) ;
	    if (debug)
	      {
		printf ("After product\n") ;
		diamShow (lMax, diams, 0) ;
	      }
	    diaminoDestroy (oldDiams) ;
	    diaminoDestroy (subDiams) ;
	  }
	else
	  {
	    diaminoDestroy (diams) ;
	    diams = subDiams ;
	  }
	break ;
      }
  if (debug)
    {
      printf ("Final result of level %d\n", level) ;
      diamShow (lMax, diams, 0) ;
    }
  if (debug)  printf ("diamEval(%d) ends\n", level) ;
  return diams ;
}

/**********************************************************************************/
/**********************************************************************************/

static int diamSplitOrder (const void *a, const void *b)
{
  BitSet bba = ((DIAMWORD *)a)->bb, bbb = ((DIAMWORD *)b)->bb ;
  unsigned long int i, nn = bigArrayMax(bba) ;
  unsigned int *ia, *ib ;

  if (nn)
    {
      ia = bigArrp (bba, 0, unsigned int) ;
      ib = bigArrp (bbb, 0, unsigned int) ;
      for (i = 0 ; i < nn ; ia++, ib++, i++)
	{
	  if (*ia < *ib) return -1 ;
	  if (*ia > *ib) return 1 ;
	}
    }
  return 0 ;
}

/**********************************************************************************/
/* because of breaking i now have diams containing one another, clean them out */
static Array diamCompress (Array oldDiams, DIAM *diam)
{
  Array diams = arrayCreate (arrayMax(oldDiams) , DIAMWORD) ;
  BitSet bba, bbb ;
  int ii, jj, kk, isIn,isOk ;
  DIAMWORD *dmw, *dmw2 ;
  unsigned int *ia, *ib ;
 
  for (ii = kk = 0 ; ii < arrayMax(oldDiams) ; ii++)
    {
      unsigned long int i, nn ;
      dmw = arrp (oldDiams, ii, DIAMWORD) ;
      bba = dmw->bb ;
      nn = bigArrayMax(bba) ;
      for (isIn = 0, jj = ii + 1 ; !isIn && jj < arrayMax(oldDiams) ; jj++)
	{
	  dmw2 = arrp (oldDiams, jj, DIAMWORD) ;
	  bbb= dmw2->bb ;
	  ia = bigArrp (bba, 0, unsigned int) ;
	  ib = bigArrp (bbb, 0, unsigned int) ;
	  for (i = 0, isOk = 1 ; isOk &&  i < nn ; ia++, ib++, i++)
	    {
	      if (*ia & (~(*ib))) /* one bit of a is not set in b */
		isOk = 0 ;  /* a not included in b */
	    }
	  if (isOk)
	    isIn = 1 ;
	}
      if (!isIn)
	{
	  dmw2 = arrayp (diams, kk++, DIAMWORD) ;
	  *dmw2 = *dmw ; dmw->bb = 0 ; /* do not destroy */
	}
      else
	bitSetDestroy (dmw->bb) ;
    } 
  diaminoDestroy (oldDiams) ;
  return diams ;
}

 /**********************************************************************************/
/* break the graph into connected commponents, and return product of individual results */
static Array diamSplitDisconnectedParts (Array oldDiams, DIAM *diam)
{
  int lMax = diam->lMax ;
  KEYSET pp = arrayCreate (lMax, KEY) ;
  Array diams = arrayCreate (arrayMax(oldDiams) , DIAMWORD) ;
  BitSet bb, bb2 ;
  int ii, jj = 0, i, j, k, pi, pj, pmax ;
  DIAMWORD *dmw, *dmw2 ;

  for (ii = jj = 0 ; ii < arrayMax(oldDiams) ; ii++)
    {
      dmw = arrp (oldDiams, ii, DIAMWORD) ;
      bb = dmw->bb ;
      /* number the connected components of this word */
      pp = keySetReCreate (pp) ;
      for (i = pmax = 0 ; i < lMax ; i++)
	if (bit (bb, i))
	  {
	    pi = keySet (pp, i) ;
	    if (!pi)
	      pi = keySet (pp, i) = ++pmax ;
	    for (j = i + 1 ; j < lMax ; j++)
	      if (bit (bb, j))
		{ 
		  pj = keySet (pp, j) ;
		  if (pi == pj ||
		      !(diam->isConnected)(diam->vp, i, j))
		    continue ;
		  if (!pj)
		     keySet (pp, j) = pi ;
		  else
		    for (k = 0 ; k < lMax ; k++)
		      if (bit (bb, k) &&  keySet (pp, k) == pj)
			keySet (pp, k) = pi ;
		}
	  }
      /* report in different new diamwords the different connected components */
      for (pi = 1 ; pi <= pmax ; pi++)
	{
	  dmw2 = arrayp (diams, jj++, DIAMWORD) ;
	  dmw2->discard = FALSE ;
	  dmw2->p = 0 ;
	  bb2 = dmw2->bb = bitSetCreate (lMax, 0) ;
	  
	  for (i = 0 ; i < lMax ; i++)
	    if (bit (bb, i) &&  keySet (pp, i) == pi)
	      { bitSet (bb2, i) ; dmw2->p++ ; }
	}
    }
  keySetDestroy (pp) ;
  diaminoDestroy (oldDiams) ;
  arraySort (diams, diamSplitOrder) ;
  diams = diamCompress (diams, diam) ;

  return diams ;
}

/**********************************************************************************/

static Array diamRun (DIAM *diam, BitSet mask, int level)
{
  int lMax = diam->lMax, ctgMax  ;
  Array contig = arrayCreate (lMax, int) ;
  Array cardinal = arrayCreate (lMax, int) ;
  Array diams = 0 ;

  array (contig, lMax - 1, int) = 0 ; /* make room */
  array (cardinal, lMax - 1, int) = 0 ; /* make room */

  ctgMax = diamContigs (lMax, diam, contig, mask, cardinal) ;
  diams = diamEvaluate (diam, ctgMax, contig, mask, cardinal, level) ;
  arrayDestroy (cardinal) ;
  arrayDestroy (contig) ;

  return diams ;
}

/**********************************************************************************/
/* invoke via the debugger */
static void diamShow (int lMax, Array diams, int target)
{ 
  DIAMWORD *mm ;
  int ii, j, dMax = diams ? arrayMax(diams) : 0 ;

  switch (dMax)
    {
    case 0:
      printf ("Empty diamino\n") ;
      break ;
    case 1:
      printf ("Single word in diamino\n") ;
      break ;
    default:
      printf ("%d words in diamino\n", dMax) ;
      break ;
    }

  for (ii = 0 ; ii < dMax ; ii++)
    {
      mm = arrp (diams, ii, DIAMWORD) ; 
      if (target && ! bit(mm->bb,target))
	continue ;
      printf ("%2d:: %d ", ii + 1, mm->p) ;
      if (mm->discard) printf ("--- ") ;
      else printf ("+++ ") ;
      for (j = 0 ; j < lMax ; j++)
	if (bit(mm->bb,j)) printf("%d ",j);
      printf ("\n") ;
    }
  if (FALSE) diamShow (lMax, diams, 0) ; /* to please the compiler */
}

/**********************************************************************************/
/**********************************************************************************/

static BOOL diamTestCompatible (void *vp, int i1, int i2)
{
  int n = 100*i1 + i2 ;

  switch (n)
    {
    case 104: case 204: case 401: case 402:
      return FALSE ;
    }
  return TRUE ;
}

/**********************************************************************************/

static BOOL diamTestConnected (void *vp, int i1, int i2)
{
  return TRUE ;
}

/**********************************************************************************/
/**********************************************************************************/
/*  public function */

void diaminoTest (void)
{
  Array aa = 0 ;
  aa = diaminoCreate (NULL, 5, diamTestCompatible, diamTestConnected, 0) ;
  diaminoDestroy (aa) ;
}

/**********************************************************************************/

Array diaminoCreate (void *vp, int lMax, DCF isCompatible, DCF isConnnected, int diamMax)
{
  Array diams = 0 ;
  BitSet mask = bitSetCreate (lMax, 0);
  DIAM *diam = diamCreate (vp, lMax, isCompatible, isConnnected, diamMax) ;
		   
  diamLinks (diam) ; 
  diams = diamRun (diam, mask, 0) ;

  if (0) diamShow (lMax, diams, 0) ;
  if (!diamMax || arrayMax(diams) < diamMax)
    diams = diamSplitDisconnectedParts (diams, diam) ;
  if (0) diamShow (lMax, diams, 0) ;

  diamDestroy (diam) ;
  bitSetDestroy (mask) ;
  return diams ;
}

/**********************************************************************************/

void diaminoDestroy (Array aa)
{
  int ii = aa ? arrayMax(aa) : 0 ; /* to be sure we do not leak empty mm->bb */
  DIAMWORD *mm ;

  while (ii--)
    {
      mm = arrayp (aa,ii, DIAMWORD) ;
      bitSetDestroy (mm->bb) ;
    }
  arrayDestroy (aa) ;
}

/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/

