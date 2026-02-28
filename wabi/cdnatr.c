/*  File: cdnatr.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg. 2001
 *-------------------------------------------------------------------
 * This file is part of the ACEMBLY extension to 
        the ACEDB genome database package, written by
 *         Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *        Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
   Align cDNA
 * Exported functions:
 * HISTORY:
 * Created: Mars 2001 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */
#define CHRONO


#include "acedb.h"
#include "bitset.h"
#include "bs.h"
#include "cdna.h"
#include "cdnatr.h"
#include "diamino.h"
#include "makemrna.h"
#include "freeout.h"
#include "chrono.h"
#include "query.h"
#include "lex.h"

BOOL isComposite = FALSE ;

/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/

/* A shadow contig, 
   is a set of shadows with at least one exon in common

   An individual shadows will be the letters of the 
   danielle algo

   A shadow (SH) corresponds to one or a set of cdna_clones
   together with yes/no bitSets of genome pieces used as exons or introns

   The coordinates of the exons/introns flagged by the bitSets are held in chain
*/

typedef struct chainStruct { int x ; BOOL isUp ; BOOL isGap ; } CHAIN ;

/**********************************************************************************/

static BOOL s2mConnectedShadows (SH *sh1, SH *sh2, int bMax)
{
  BitSet yes1, yes2, ghost1, ghost2, no1, no2 ;
  int i = bMax, j1, j2 ;
  BOOL isDown = TRUE ;

  if (sh1->a2 < sh2->a1 || sh2->a2 < sh1->a1)
    return FALSE ;
  yes1 = sh1->yes ;
  yes2 = sh2->yes ;
  no1 = sh1->no ;
  no2 = sh2->no ;
  ghost1 = sh1->ghost ;
  ghost2 = sh2->ghost ;

  /* simple case: exon is shared */
  i = bMax ;
  while (i--)
    if ((bit(yes1,i) && bit(yes2,i)))
      return TRUE ;

  /* ghost touches exon */
  i = bMax ;
  while (i--)
    if (bit(ghost1,i) && bit(yes2,i)) 
      {
        i = bMax ;
        while (i--)
          if ((bit(yes1,i) && bit(no2,i)) ||
              (bit(yes2,i) && bit(no1,i)))
            return FALSE ;  /* incompatibility */
        /* no incompatibility, check that we touch seriously */
        /* check orientation of ghost */
        i = bMax ;
        for (i = 0 ; i < bMax ; i++)
          {
            if (bit(yes1,i)) { isDown = TRUE ; break ; }
            if (bit(ghost1,i)) { isDown = FALSE ; break ; }
          }
        if (isDown) /* go down, wait for yes1, start counting */
          {
            for (i = 0 ; i < bMax ; i++) 
              if (bit(yes1,i))
                break ;
            for (j1 = j2 = 0 ; i < bMax ; i++)
              {
                if (bit(yes2,i)) j1++ ;
                if (bit(ghost1,i) && bit(yes2,i))  j2++ ;
              }
            if (j2 > 0 && (j2 >= 2 || (j1 == 1 && j2 == 1)))
              return TRUE ;
            return FALSE ;
          }
        else /* go up, wait for yes1, start counting */
          {
            for (i = bMax - 1 ; i >= 0 ; i--)
              if (bit(yes1,i))
                break ;
            for (j1 = j2 = 0 ; i >= 0 ; i--)
              {
                if (bit(yes2,i)) j1++ ;
                if (bit(ghost1,i) && bit(yes2,i))  j2++ ;
              }
            if (j2 > 0 && (j2 >= 2 || (j1 == 1 && j2 == 1)))
              return TRUE ;
            return FALSE ;
          }
        return FALSE ;
      }
  /* ghost touches exon */
  i = bMax ;
  while (i--)
    if (bit(ghost2,i) && bit(yes1,i)) 
      {
        i = bMax ;
        while (i--)
          if ((bit(yes1,i) && bit(no2,i)) ||
              (bit(yes2,i) && bit(no1,i)))
            return FALSE ;  /* incompatibility */
        /* no incompatibility, check that we touch seriously */
        /* check orientation of ghost */
        i = bMax ;
        for (i = 0 ; i < bMax ; i++)
          {
            if (bit(yes2,i)) { isDown = TRUE ; break ; }
            if (bit(ghost2,i)) { isDown = FALSE ; break ; }
          }
        if (isDown) /* go down, wait for yes1, start counting */
          {
            for (i = 0 ; i < bMax ; i++)
              if (bit(yes2,i))
                break ;
            for ( j1 = j2 = 0 ; i < bMax ; i++)
              {
                if (bit(yes1,i)) j1++ ;
                if (j1 > 0 && bit(ghost2,i) && bit(yes1,i))  j2++ ;
              }
            if (j2 > 0 && (j2 >= 2 || (j1 == 1 && j2 == 1)))
              return TRUE ;
            return FALSE ;
          }
        else /* go up, wait for yes1, start counting */
          { 
            for (i = bMax - 1 ; i >= 0 ; i--)
              if (bit(yes2,i))
                break ;
            for (j1 = j2 = 0 ; i >= 0 ; i--)
              {
                if (bit(yes1,i)) j1++ ;
                if (j1 > 0 && bit(ghost2,i) && bit(yes1,i))  j2++ ;
              }
            if (j2 > 0 && (j2 >= 2 || (j1 == 1 && j2 == 1)))
              return TRUE ;
            return FALSE ;
          }
        return FALSE ;
      }
  return FALSE ;
}

/**********************************************************************************/

static BOOL s2mCompatibleShadows (SH *sh1, SH *sh2, int bMax)
{
  BitSet yes1, yes2, no1, no2, suspect1, suspect2 ;
  int i = bMax, j, end3p1, end3p2, begin3p1, begin3p2 ;
  
  if (sh1->a2 < sh2->a1 || sh2->a2 < sh1->a1)
    return FALSE ;

  yes1 = sh1->yes ; no1 = sh1->no ;
  yes2 = sh2->yes ; no2 = sh2->no ;
  suspect1 = sh1->suspect ; suspect2 = sh2->suspect ;
  end3p1 = sh1->end3p ;
  end3p2 = sh2->end3p ;
  begin3p1 = sh1->begin3p ;
  begin3p2 = sh2->begin3p ;

  i = bMax ;
  while (i--)
    if (
        (!bit(suspect1,i) && !bit(suspect2,i)) &&
        (
         (
          (bit(yes1,i) && bit(no2,i)) ||
          (bit(yes2,i) && bit(no1,i)) ||
          (bit(suspect1,i) && bit(no2,i)) ||
          (bit(suspect2,i) && bit(no1,i))
          )
         )
        )
      return FALSE ;

  i = bMax ;
  while (i--)
    {
      if (bit(suspect1,i)) /* check i hit an exon on both side */
        {
          for (j = i-1 ; j >= 0 ; j--)
            if (bit(yes1,j) && ! bit(yes2,j))
              return FALSE ;
          for (j = i+1 ; j < bMax ; j++)
            if (bit(yes1,j) && ! bit(yes2,j))
              return FALSE ;
        }
      if (bit(suspect2,i)) /* check i hit an exon on both side */
        {
          for (j = i-1 ; j >= 0 ; j--)
            if (bit(yes2,j) && ! bit(yes1,j))
              return FALSE ;
          for (j = i+1 ; j < bMax ; j++)
            if (bit(yes2,j) && ! bit(yes1,j))
              return FALSE ;
        }
    }
  /* example XQ756 */
  if (end3p1)   /* after a 3prime end i want no intron or gap followed by exon */
    for (i = end3p1 + 1 ; i < bMax ; i++)
      if (bit(no2,i))  /*  an intron */
        {
          for (; i < bMax ; i++)
            if (bit(yes2,i)) 
              return FALSE ;
          break ;
        }
  if (end3p2)   /* after a 3prime end i want no intron followed by exon */
    for (i = end3p2 + 1 ; i < bMax ; i++)
      if (bit(no1,i))  /* after an intron */
        {
          for (; i < bMax ; i++)
            if (bit(yes1,i)) 
              return FALSE ;
          break ;
        }

  if (begin3p1)   /* after a 3prime begin i want no intron followed by exon */
    for (i = begin3p1 - 1 ; i >= 0 ; i--)
      if (bit(no2,i))  /* after an intron */
        {
          for (; i >= 0 ; i--)
            if (bit(yes2,i)) 
              return FALSE ;
          break ;
        }
  if (begin3p2)   /* after a 3prime begin i want no intron followed by exon */
    for (i = begin3p2 - 1 ; i >= 0 ; i--)
      if (bit(no1,i))  /* after an intron */
        {
          for (; i >= 0 ; i--)
            if (bit(yes1,i)) 
              return FALSE ;
          break ;
        }

  return TRUE ;
}

/**********************************************************************************/

static BOOL s2mIdenticalShadows (SH *sh1, SH *sh2, int bMax)
{
  BitSet yes1, yes2, no1, no2 ;
  int i = bMax ;
  
  yes1 = sh1->yes ; no1 = sh1->no ;
  yes2 = sh2->yes ; no2 = sh2->no ;

  while (i--)
    if (
        (bit(yes1,i) != bit(yes2,i)) ||
        (bit(no1,i) !=  bit(no2,i))
        )
      return FALSE ;
  return TRUE ;
}

/**********************************************************************************/

static void s2mKillShadow (SH *sh)
{
  if (keySetExists (sh->clones))
    {
      keySetDestroy (sh->clones) ;
      bitSetDestroy (sh->yes) ;
      bitSetDestroy (sh->no) ;
      bitSetDestroy (sh->ghost) ;
      bitSetDestroy (sh->suspect) ;
    }
}

/**********************************************************************************/

static void s2mKillContig (SC *sc)
{ 
  if (keySetExists (sc->sh.clones))
    {
      keySetDestroy (sc->ks) ;
      keySetDestroy (sc->sh.clones) ;
      bitSetDestroy (sc->sh.yes) ;
      bitSetDestroy (sc->sh.no) ;
      bitSetDestroy (sc->sh.ghost) ;
      bitSetDestroy (sc->sh.suspect) ;
    }
}

/**********************************************************************************/

static void s2mDoMergeShadows (SH *sh1, SH *sh2, int bMax)
{
  BitSet yes1, yes2, no1, no2, ghost1, ghost2, suspect1, suspect2 ;
  int i = bMax ;
  
  yes1 = sh1->yes ; no1 = sh1->no ; ghost1 = sh1->ghost ; suspect1 = sh1->suspect ;
  yes2 = sh2->yes ; no2 = sh2->no ; ghost2 = sh2->ghost ; suspect2 = sh2->suspect ;

  while (i--)
    {
      if (bit(suspect2,i) && !bit(yes1,i) && !bit(no1,i)) bitSet (suspect1,i) ;
      if (bit(suspect1,i) && !bit(suspect2,i) && (bit(yes2,i) || bit(no2,i)))
        bitUnSet (suspect1,i) ;
      if (bit(yes2,i)) bitSet (yes1,i) ;
      if (bit(no2,i)) bitSet (no1,i) ;
      if (bit(ghost2,i)) bitSet (ghost1,i) ;
    }

  for (i = 0 ; i < keySetMax (sh2->clones) ; i++)
    keySetInsert (sh1->clones, keySet (sh2->clones, i)) ;
  
  if (sh2->begin3p &&
      (!sh1->begin3p || sh1->begin3p > sh2->begin3p))
    sh1->begin3p  = sh2->begin3p ;

  if (sh2->end3p &&
      (!sh1->end3p || sh1->end3p < sh2->end3p)) 
    sh1->end3p  = sh2->end3p ;
  if (sh1->a1 > sh2->a1) sh1->a1 = sh2->a1 ;
  if (sh1->a2 < sh2->a2) sh1->a2 = sh2->a2 ;
}

/**********************************************************************************/

static void s2mDoMergeContigs (SC *sc1, SC *sc2, int bMax)
{
  int i = bMax ;

  s2mDoMergeShadows (&(sc1->sh), &(sc2->sh), bMax) ;

  for (i = 0 ; i < keySetMax (sc2->ks) ; i++)
    keySetInsert (sc1->ks, keySet (sc2->ks, i)) ;
}

/**********************************************************************************/

static void s2mCompressShadows (Array shadows)
{
  int ii, jj ;
  SH* sh1, *sh2 ;

  for (ii = 0, jj = 0 ; ii < arrayMax(shadows) ; ii++)
    {
      sh1 = arrp (shadows, ii, SH) ;
      if (sh1->clones)
        {
          if (jj < ii)
            {
              sh2 = arrp (shadows, jj, SH) ;
              *sh2 = *sh1 ;
              sh1->clones = 0 ; sh1->yes = 0 ; sh1->no = 0 ; sh1->ghost = 0 ; sh1->suspect = 0 ;
            }
          jj++ ;
        }
    }
  arrayMax(shadows) = jj ;
}

/**********************************************************************************/

static void s2mCompressContigs (Array contigs)
{
  int ii, jj ;
  SC* sc1, *sc2 ;

  for (ii = 0, jj = 0 ; ii < arrayMax(contigs) ; ii++)
    {
      sc1 = arrp (contigs, ii, SC) ;
      if (sc1->ks)
        {
          if (jj < ii)
            {
              sc2 = arrp (contigs, jj, SC) ;
              *sc2 = *sc1 ;
              sc1->sh.clones = 0 ; sc1->ks = 0 ; sc1->sh.yes = 0 ; sc1->sh.no = 0 ;
            }
          jj++ ;
        }
    }
  arrayMax(contigs) = jj ;
}

/**********************************************************************************/

static void s2mShowOneShadow (int ii, SH *sh, KEY targetClone)
{
  KEYSET ks = 0 ;
  int i ;
 
  if (sh->clones)
    {
      if (targetClone)
        {
          for (i = 0 ; i < keySetMax (sh->clones) ; i++)
            if (keySet (sh->clones, i) == targetClone)
              break ;
          if (i >= keySetMax (sh->clones))
            return ; 
        }
      ks = keySetCopy (sh->clones) ;
      arraySort (ks, keySetAlphaOrder) ;
      printf ("%d:: ln=%d  a1=%d a2=%d %d sh = ", ii, sh->a2 - sh->a1, sh->a1, sh->a2, keySetMax (ks)) ;
      for (i = 0 ; i < keySetMax (ks) ; i++)
        printf ("%s ", name(keySet (ks, i))) ;
      printf ("\n") ;
      keySetDestroy (ks) ;
      
      printf ("  yes :: ") ;
      for (i = 0 ; i < sh->bMax ; i++)
        if (bit (sh->yes, i)) printf ("%d ", i) ;
      printf ("\n") ;
      printf ("   no :: ") ;
      for (i = 0 ; i < sh->bMax ; i++)
        if (bit (sh->no, i)) printf ("%d ", i) ;
      printf ("\n") ;
      printf ("ghost :: ") ;
      for (i = 0 ; i < sh->bMax ; i++)
        if (bit (sh->ghost, i)) printf ("%d ", i) ;
      printf ("\n") ;
      printf ("suspect :: ") ;
      for (i = 0 ; i < sh->bMax ; i++)
        if (bit (sh->suspect, i)) printf ("%d ", i) ;
      printf ("\n") ;
    }
}

static void s2mShowShadows (Array shadows, char *title)
{
  int ii ;
  KEY clo = 0 ;

  if (! title || ! lexword2key (title, &clo, _VcDNA_clone))
    clo = 0 ;
  if (shadows)
    {
      printf ("Currently %d %s\n", arrayMax(shadows), title) ;
      for (ii = 0 ; ii < arrayMax(shadows) ; ii++) 
        s2mShowOneShadow (ii, arrp (shadows, ii, SH), clo) ;
    }
  else
    printf ("sorry shadows == 0\n") ;
}

static void s2mShowContigs (Array contigs, char *title, int qual)
{
  int ii ;
  KEY clo = 0 ;

  if (! title || ! lexword2key (title, &clo, _VcDNA_clone))
    clo = 0 ;

  printf ("s2mShowContigs qual=%d Currently %d %s\n", qual, arrayMax(contigs), title) ;
  for (ii = 0 ; ii < arrayMax(contigs) ; ii++) 
    {
      printf("contig %d ma(ks)=%d ", ii, arrayMax(arr(contigs, ii, SC).ks)) ;
      s2mShowOneShadow (ii, &(arr(contigs, ii, SC).sh), clo) ;
    }
}

/**********************************************************************************/
/**********************************************************************************/

static void* s2mMagic = 0 ;
#define s2mDestroy(_s2m) (s2mDoDestroy(_s2m), s2m=0)
static S2M* s2mCreate (int type, KEY cosmid, KEY cosmid1, Array dnaD, Array dnaR, Array plainHits, Array geneHits, Array linkPos)
{
  AC_HANDLE h = handleCreate () ;
  S2M *s2m = (S2M*) halloc (sizeof (struct s2mStruct), h) ;

  s2m->h = h ;
  s2m->type = type ;
  s2m->cosmid = cosmid ;
  s2m->cosmid1 = cosmid1 ;
  s2m->dnaD = dnaD ;
  s2m->dnaR = dnaR ;
  s2m->plainHits = plainHits ;
  s2m->geneHits = geneHits ;
  s2m->gmrnas = arrayHandleCreate (12, SMRNA, s2m->h) ; 
  s2m->linkPos = linkPos ;
  s2m->magic = &s2mMagic ;
  return s2m ;
}

/**********************************************************************************/

static void s2mDoDestroy (S2M *s2m)
{
  if (s2m->h)
    {
      AC_HANDLE h = s2m->h ;
      messfree (h) ; /* in 2 lines because mesfree resets h to zero ! */
    }
}

/**********************************************************************************/
/**********************************************************************************/
/* a chain is a continuous set of sections of type exon, intron, stop, gap */
static int chainOrder (const void *a, const void *b)
{
  const CHAIN *x = (const CHAIN*)a , *y = (const CHAIN*)b ;
  if (x->isUp != y->isUp)
    return x->isUp ? 1 : -1 ;
  return x->x - y->x ;
}

/**********************************************************************************/

static void s2mMakeChain (S2M *s2m, BOOL debug)
{
  Array chain, geneHits ;
  HIT *hh, *hh1 ;
  int ii, jj ;
  BOOL isUp ;
  CHAIN *ch = 0 ;

  geneHits = s2m->geneHits ;
  chain = s2m->chain = arrayHandleCreate (200, CHAIN, s2m->h) ;

  /* beware of fuzzy introns */
  for (ii = 0 ; ii < arrayMax(geneHits) ; ii++)
    {
      int ok ;
      hh = arrp (geneHits, ii, HIT) ;
      if (!class(hh->est)) /* ghost */
        continue ;
      if (! (hh->type & gJ))
        continue ;
      ok = 0 ;
      for (jj = ii - 1, hh1 = hh - 1 ; jj >= 0 ; hh1--, jj--)
        {
          if (!(hh1->type & gX))
            continue ;
          if (hh1->a2 > hh->a1 - 4 && hh1->a2 < hh->a1 + 4)
            { hh->a1 = hh1->a2 + 1 ; ok++ ; break ; }
        }
      for (jj = ii + 1, hh1 = hh + 1 ; jj  < arrayMax(geneHits) ; hh1++, jj++)
        {
          if (!(hh1->type & gX))
            continue ;
          if (hh1->a1 > hh->a2 - 4 && hh1->a1 < hh->a2 + 4)
            { hh->a2 = hh1->a1 - 1 ; ok++ ; break ; }
        }
      if (ok < 2) /* do not leave floating fuzzy introns around they destroy the est */
        hh->est = 0 ;
    }

  for (ii = 0, jj = 0 ; ii < arrayMax(geneHits) ; ii++)
    {
      hh = arrp (geneHits, ii, HIT) ;
      if (!class(hh->est)) /* ghost */
        continue ;
      isUp = hh->reverse ;
      ch = arrayp (chain, jj++, CHAIN) ;
      ch->isUp = isUp ; ch->x = hh->a1 - 1 ;
      ch = arrayp (chain, jj++, CHAIN) ;
      ch->isUp = isUp ; ch->x = hh->a2 ;
    }

  if (debug)
    {
      printf ("Chain:: ") ;
      for (ii = 0 ; ii < arrayMax(chain) ; ii++)
        { 
          ch = arrayp (chain, ii, CHAIN) ;
          printf ("%c%d ", ch->isUp ? '-' : '+', ch->x) ;
        }
      printf ("\n") ;
    }
  arraySort (chain, chainOrder) ;
  arrayCompress (chain) ;
  s2m->bMax = arrayMax(chain) ;

  if (debug)
    {
      printf ("Chain:: ") ;
      for (ii = 0 ; ii < arrayMax(chain) ; ii++)
        { 
          ch = arrayp (chain, ii, CHAIN) ;
          printf ("%c%d ", ch->isUp ? '-' : '+', ch->x) ;
        }
      printf ("\n") ;
    }
  return ;
} /* s2mMakeChain */

/**********************************************************************************/

/* a shadow is the projection as a bit set of a single clone */
static void s2mMakeOneShadow (S2M *s2m, BOOL isUp, int j0, int j1, KEYSET coverage)
{
  HIT *vp ;
  SH *sh ;
  int spy, ii, i, j, a1, a2, c1, c2, x1, x2, bMax = s2m->bMax, j10, j11, limit ;
  Array chain = s2m->chain ;
  CHAIN *ch ;
  static BitSet nono = 0 ;
  BitSet yes, no, ghost, suspect ;
  BOOL iFirst, ignoreSL0 = TRUE, ignorePolyA = FALSE, debug = FALSE ;
  int iLast, blankFirst, blankLast ;
  KEY oldEst ;
  KEYSET thisCoverage = 0 ;

  vp = arrp (s2m->plainHits, j0, HIT) ;
  if (vp->type & gDroppedGene)
    return ;

  thisCoverage = arrayCreate (bMax, KEY) ;
  j10 = j1 ; j11 = j0 ;
  for (ii = j0 ; ii < j1 ; ii++)
    {   
      vp = arrp (s2m->plainHits, ii, HIT) ; 
      if (vp->type & (gX | gFuseToGhost)) /* start on exon */
        { j10 = ii ; break ; }
    }
  for (ii = j10 ; ii < j1 ; ii++)
    {   
      vp = arrp (s2m->plainHits, ii, HIT) ; 
      if (vp->type & (gX | gFuseToGhost)) /* end on exon */
        j11 = ii ;
    }
  if (j11 < j10)  /* only ghosts */
    return ;
 
  sh = arrayp (s2m->shadows, arrayMax(s2m->shadows), SH) ;
  sh->bMax = bMax ;
  sh->isUp = isUp ;
  sh->clones = arrayHandleCreate (8, KEY, s2m->h) ;
  keySet (sh->clones, 0) = vp->cDNA_clone ;

  nono =  bitSetReCreate (nono, bMax) ;
  sh->yes = yes = bitSetCreate (bMax, s2m->h) ;
  sh->no = no = bitSetCreate (bMax, s2m->h) ;
  sh->ghost = ghost = bitSetCreate (bMax, s2m->h) ;
  sh->suspect = suspect = bitSetCreate (bMax, s2m->h) ;

  oldEst = 0 ; ii = j0 ; /* run separately on each est to handle correctly the SL and polyA */
  blankFirst = 0 ; blankLast = bMax ; 
 newEst:
  oldEst = 0 ;
  iFirst = TRUE ; iLast = bMax ;
  for (; ii < j1 ; ii++)
    {   
      vp = arrp (s2m->plainHits, ii, HIT) ; 
      if (oldEst && vp->est != oldEst)
        { vp-- ; break ; }

      oldEst = vp->est ;
      /*
        up = ii > j0 ? vp - 1 : 0 ;      
        wp = ii < j1 - 1 ? vp + 1 : 0 ;  
      */
      a1 = vp->a1 ; a2 = vp->a2 ;
      for (i = 0, c1 = 1, ch = arrp(chain, 0, CHAIN) ; i < bMax ; i++, ch++)
        {
          c2 = ch->x ;
          
          x1 = a1 > c1 ? a1 : c1 ;
          x2 = a2 < c2 ? a2 : c2 ;
          limit = (vp->nerrAll ? (vp->nerrAll < 5 ? 10 : 5 + 2*vp->nerrAll) : 5) ;
          spy=0 ;
          if (x1 <= x2 - limit ||                  /* big unambiguous overlap */
              (a1 < c1 - 5 && a2 > c2 + 5) ||      /* the c segment is clearly within the a segment */
              ( 
               x2 - x1 == c2 - c1 &&          /* small exact */
               ! ((vp->type & gA) &&      /* but not a polyA wobble */
                  (x2 - x1 < 5 || vp->nerr) && /* which is shorther longer than 5 or noisy */
                  (                              
                   (!vp->reverse &&  x2 > a2 - limit) ||
                   (vp->reverse &&  x1 < a1 + limit) 
                   )
                  ) &&
               (spy+=1) &&
               ! ( /* a small 5p extension of less than 3 bp */
                   x2 - x1 == c2 - c1 &&
                   x2 - x1 <= 2 &&
                   ((ii == j10 && !vp->reverse && vp->a2 > x2) || (ii == j11 && vp->reverse && vp->a1 < x1)) &&
                   !(vp->type & gS) &&
                   (spy+=2)/* except if a sl itself */
                  ) &&
               (spy+=10)
               ) ||
              (
               x2 - x1 < limit && /* very short probably team jumped */ 
               !vp->reverse &&         
               (
                (vp->x1 < vp->x2 && a2 == c2) ||
                (vp->x1 > vp->x2 && a1 == c1) 
                ) &&
               (spy+=100)
               ) ||
              (
               (   /* very short zone completely inside the clone exon, so is exon */
                (vp->type & gX) && a1 < x1  && a2 > x2  /* 4L677, there is an SL0 at every base  */
               ) && 
               (spy+=200)
               )
              )
            {
              if (debug)
                printf ("est %s spy=%d x1=%d x2=%d a1=%d a2=%d\n    ",
                        name(vp->est), spy, vp->x1, vp->x2, vp->a1, vp->a2) ; 

              if (TRUE &&
                  vp->type & gSuspect) 
                {
                  /* 2005_04_19:
                   * wish make suspect non compatible with big introns  over 30 bp 
                   * so i simply do not flag short suspects
                   * 2005_06_12 -> FALSE exemple 1D214
                   * je supprime ce flagging
                   * en fait gSuspect signifie suspected_internal_deletion
                   * et donc compatible a exon et introns alike
                   */
                  if (x2 - x1 > 30)
                    bitSet (suspect, i) ;
                }
              else if (vp->type & gMicro) 
                {
                  /* 2006_12_07 see BG189081 on S_Y_3 human
                   */
                  bitSet (suspect, i) ;
                }
              else if (! class(vp->est) && (vp->type & gFuseToGhost) )           /* ghost */
                {
                  if (x1 <= x2 - 20  && x1 <= a2 - 100) /* CF rps-10 */
                    { bitSet (ghost, i) ; bitSet (yes, i) ; }
                }
              else
                {
                  if (vp->type & (gX | gFuseToGhost)) 
                    { 
                      if (debug)
			printf("yes %d ",i) ;
                      bitSet (yes, i) ; keySet(thisCoverage,i) += x2 - x1 + 1 ;
                     }/* exon */
                  if (vp->type & gI) 
                    {
                      if (debug) 
			printf("no %d ",i)  ;
                      bitSet (no, i) ;  /* intron */
                      if (x2 - x1 > 50)
                        bitSet (nono, i) ;  /* absolute intron */
                    }
                  if (!isUp && iFirst && (vp->type & gS))                   /* transpliced leader sl */
                     blankFirst = i ;
                  if (((vp->type & gReal5p) || !ignoreSL0) && 
                      !isUp && iFirst &&  (vp->type & gS0))   /* transpliced leader sl0 */
                     blankFirst = i ;
                  if (isUp && iFirst && !ignorePolyA && (vp->type & gA))                   /* polyA */
                    { sh->begin3p = i ;  if ((vp->type & gReal3p)) blankFirst = i ; }
                  iFirst = FALSE ; iLast = i ;
                  if (debug) printf("\n") ;
		  if (!sh->a2 || vp->a1 < sh->a1) sh->a1 = vp->a1 ;
		  if (!sh->a2 || vp->a2 > sh->a2) sh->a2 = vp->a2 ;
                }
            }
          c1 = c2 + 1 ;
        }
    }

  if (oldEst)
    {
      if (!isUp && !ignorePolyA && (vp->type & gA))                   /* polyA */
        { sh->end3p = iLast ; if ((vp->type & gReal3p)) blankLast = iLast ; }
      if (isUp && (vp->type & gS))                   /* SL1 */
        blankLast = iLast ;
      if (((vp->type & gReal5p) || !ignoreSL0)  && isUp && (vp->type & gS0))                   /* SL0 */
        blankLast = iLast ;
      goto newEst ;
    } 
  
  /* in case the other part of the clone contradicts a blanking, move the blanking */
  for (j = blankFirst - 1 ; j >= 0 ;  j--)  
    if (bit (yes, j)) 
      blankFirst = j ;
  for (j = 0 ; j < blankFirst ; j++)  
    bitSet (no, j) ;
  for (j = blankLast + 1 ; j < bMax ; j++)  
    if (bit (yes, j)) 
      blankLast = j ;        
  for (j = blankLast + 1 ; j < bMax ; j++)  
    bitSet (no, j) ;
  for (j = sh->end3p + 1 ; j < bMax ; j++)
    if (bit (yes, j))
      sh->end3p   = 0 ; /* ignore an end contradicted by the other EST */
  for (j = sh->begin3p - 1 ; j >= 0 ; j--)
    if (bit (yes, j))
      sh->begin3p   = 0 ; /* ignore a begin contradicted by the other EST */
  /* because of allowed wobble, some clones may be undecidable */
  for (i = 0 ; i < bMax ; i++)
    if (bit (yes, i) && bit (nono, i))
      { bitUnSet (yes, i) ; }
  for (i = 0 ; i < bMax ; i++)
    if (bit (yes, i) && bit (no, i))
      { bitUnSet (yes, i) ; bitUnSet (no, i) ; }
    else
      {
        if (keySet (coverage, i) < keySet (thisCoverage, i))
          keySet (coverage, i) = keySet (thisCoverage, i) ;
      }

  keySetDestroy (thisCoverage) ;
} /* s2mMakeOneShadow */

/**********************************************************************************/

static void s2mMakeShadows (S2M *s2m, BOOL debug, BOOL *hasReadPairs, BOOL *hasFusedClones)
{
  HIT *hh ;
  Array hits = s2m->plainHits ;
  KEY old = 0, oldEst = 0 ;
  BOOL oldUp = FALSE ;
  int ii, jj, bMax = s2m->bMax, c1, cc, cover ;
  KEYSET coverage = 0 ;
  CHAIN *ch ; SH*sh ;

  coverage = arrayCreate (bMax, KEY) ;

  s2m->shadows = arrayHandleCreate (32, SH, s2m->h) ;
  for (ii = jj = 0 ; ii < arrayMax(hits) ; ii++)
    {
      hh = arrp (s2m->plainHits, ii, HIT) ;
      if (old && (old != hh->cDNA_clone || hh->est != oldEst || hh->reverse != oldUp))
        { 
          if (old == hh->cDNA_clone)
            *hasReadPairs = TRUE ;
          s2mMakeOneShadow (s2m, oldUp, jj, ii, coverage) ;
          jj = ii ;
        }
      if (old != hh->cDNA_clone &&
          keyFindTag ( hh->cDNA_clone, str2tag ("Fuse_to_clone")))
        *hasFusedClones = TRUE ;
      old = hh->cDNA_clone ;
      oldUp = hh->reverse ;
      oldEst = hh->est ;
    }

  if (old)
    s2mMakeOneShadow (s2m, oldUp, jj, ii, coverage) ;

  if (arrayMax(s2m->chain))
    for (ii = 0, c1 = 0, ch = arrp(s2m->chain, ii, CHAIN) ; ii < bMax ; ch++, ii++)
      {
        cover = keySet(coverage , ii), cc = ch->x - c1 ;
        
        if (cc > 3 * cover && cover < 100)
          for (jj = 0, sh = arrp (s2m->shadows, 0, SH) ; jj < arrayMax(s2m->shadows)  ; sh++, jj++)
            {
              if (debug && bit(sh->yes, ii))
                { printf("destroying a yes jj=%d x1=%d x2=%d", jj, c1, ch->x) ;}               
              bitUnSet (sh->yes, ii) ; 
            } 
        c1 = ch->x  ;
      }
  keySetDestroy (coverage) ;
  
  return ;
}  /* s2mMakeShadows */

/**********************************************************************************/

static void s2mMergeCloneShadows (S2M *s2m)
{
  int i, ii, jj, imin,imax,  bMax = s2m->bMax, ok ;
  SH *sh1, *sh2 ;
  Array shadows = s2m->shadows ;
  KEY clone = 0 ;
  BitSet suspect = 0 ;

  /* merge reads from same clone */
  for (ii = 0 ; ii + 1 < arrayMax(shadows) ; ii++)
    {
      sh1 = arrp (shadows, ii, SH) ;
      if (!sh1->clones)
        continue ; 
      clone = keySet (sh1->clones, 0) ; ok = 0 ;
      for (jj = ii + 1 ; jj < arrayMax(shadows) ; jj++)
        {
          sh2 = arrp (shadows, jj, SH) ;
          if (sh1->isUp == sh2->isUp && 
              sh2->clones &&
              keySet (sh2->clones, 0) == clone)
            {
              if (! ok)
                {
                  ok = 1 ;
                  imin = 1 ; imax = 0 ;
                  if (!suspect)
                    suspect = bitSetCreate (bMax, 0) ;
                  for (i = 0 ; i < bMax ; i++)
                    bitSet (suspect, i) ; /* flag all */
                  /* unflag the sh1 region */
                  for (i = 0 ; i < bMax ; i++)
                    if (bit (sh1->yes, i)) 
                      { imin = i ; break ; }  /* begin */
                  for (i = bMax - 1 ; i >= 0 ; i--)
                    if (bit (sh1->yes, i)) 
                      { imax = i ; break ; } /* end */
                  for (i = imin ; i <= imax ; i++)
                    if (! bit (sh1->suspect, i)) 
                      bitUnSet (suspect, i) ; /* exclude EST gaps */
                }
              /* unflag the sh2 region */
              imin = 1 ; imax = 0 ;
              for (i = 0 ; i < bMax ; i++)
                if (bit (sh2->yes, i)) 
                  { imin = i ; break ; }  /* begin */
              for (i = bMax - 1 ; i >= 0 ; i--)
                if (bit (sh2->yes, i)) 
                  { imax = i ; break ; } /* end */
              for (i = imin ; i <= imax ; i++)
                if (! bit (sh2->suspect, i)) 
                  bitUnSet (suspect, i) ; /* exclude EST gaps */
              s2mDoMergeShadows (sh1, sh2, s2m->bMax) ;
              s2mKillShadow (sh2) ;
            }
        }
      if (ok &&
          FALSE) /* june 4 2005  
                  * this code was introduced april 19 2005 to prevent
                  * the hooking of a floating est inside a 5' 3' gap
                  * unless the floating EST touches both ends of the gap
                  * but then in the human genome we decided to remove clone pairs
                  * so now we are restricted to the worm where we do not
                  * want to create alternatives like bar-1 if i remove the
                  * prediction and complete mrna
                  * so i jump this code
		  * if i open this code, bar-1 beomes totally crazy
                  */
        {
          /* remove the suspect flag outside the exon region */  
          for (i = 0 ; i < bMax ; i++)
            if (bit (sh1->yes, i) && ! bit (sh1->ghost, i)) break ;
            else bitUnSet (suspect, i) ;
          for (i = bMax - 1 ; i >= 0 ; i--)
            if (bit (sh1->yes, i) && ! bit (sh1->ghost, i) ) break ;
            else bitUnSet (suspect, i) ;
          /* we leave as suspect the declared gaps inside an est and the region between est */
          for (i = 0 ; i < bMax ; i++)
            if (bit (suspect, i)) bitSet (sh1->suspect, i) ;
        }
    }

  /* register happy few */
  s2mCompressShadows (shadows) ;
  bitSetDestroy (suspect) ;

  return ; 
}  /* s2mMergeCloneShadows */

/**********************************************************************************/

static void s2mMergeFusedCloneShadows (S2M *s2m)
{
  int ii, jj ;
  SH *sh1, *sh2 ;
  Array shadows = s2m->shadows ;
  KEYSET ks1, ks2 ;

  /* merge reads from fused clones */
  for (ii = 0 ; ii + 1 < arrayMax(shadows) ; ii++)
    {
      sh1 = arrp (shadows, ii, SH) ;
      if (!sh1->clones)
        continue ; 
      ks1 = query (sh1->clones, ">fuse_to_clone") ;
      if (keySetMax (ks1))
        for (jj = ii + 1 ; jj < arrayMax(shadows) ; jj++)
          {
            sh2 = arrp (shadows, jj, SH) ;
            if (sh1->isUp == sh2->isUp && 
                sh2->clones)
              {
                ks2 = keySetAND (ks1, sh2->clones) ;
                if (keySetMax (ks2))
                  {
                    s2mDoMergeShadows (sh1, sh2, s2m->bMax) ;
                    s2mKillShadow (sh2) ;
                  }
                keySetDestroy (ks2) ;
              }
        }
      keySetDestroy (ks1) ;
    }

  /* register happy few */
  s2mCompressShadows (shadows) ;

  return ; 
}  /* s2mMergeFusedCloneShadows */

/**********************************************************************************/

static void s2mMergeIdenticalShadows (S2M *s2m)
{
  int ii, jj ;
  SH *sh1, *sh2 ;
  Array shadows = s2m->shadows ;

  /* merge identical */
  for (ii = 0 ; ii + 1 < arrayMax(shadows) ; ii++)
    {
      sh1 = arrp (shadows, ii, SH) ;
      if (!sh1->clones)
        continue ;
      for (jj = ii + 1 ; jj < arrayMax(shadows) ; jj++)
      {
        sh2 = arrp (shadows, jj, SH) ;
        if (!sh2->clones)
          continue ;

        if (s2mIdenticalShadows (sh1, sh2, s2m->bMax))
          {
            s2mDoMergeShadows (sh1, sh2, s2m->bMax) ;
            s2mKillShadow (sh2) ;
          }
      }
    }

  /* register happy few */
  s2mCompressShadows (shadows) ;

  return ; 
}  /* s2mMergeidenticalShadows */

/**********************************************************************************/

static void s2mCompleteShadows (S2M *s2m)
{
  int ii, jj, bMax = s2m->bMax, c1, cc ;
  SH *sh ;
  CHAIN *ch ;
  Array shadows = s2m->shadows ;

  /* short gap between 2 exon counts as exon */
  for (jj = 0, c1 = 0, ch = arrp(s2m->chain, jj, CHAIN) ; jj < bMax - 1 ; ch++, jj++)
    {
      cc = ch->x - c1 ;
      if (jj > 0 && cc < 9)
        {
          for (ii = 0 ; ii < arrayMax(shadows) ; ii++)
            {
              sh = arrp (shadows, ii, SH) ;
              if (!sh->clones)
                continue ;
              if (bit(sh->yes, jj - 1) && bit(sh->yes, jj + 1) && !bit(sh->yes, jj) && !bit(sh->no, jj))
                bitSet (sh->yes, jj) ;
            }
        }
        c1 = ch->x  ;
    }
  return ; 
}  /* s2mCompleteShadows */

/**********************************************************************************/

static int s2mIncludedShadows (SH *sh1, SH *sh2, int bMax, BOOL mergePolyA)
{
  BitSet yes1, yes2, no1, no2, suspect1, suspect2  ;
  int i = bMax, a1, a2, da = 0 ;
  int bad = 0 ;
  BOOL ok ;

  yes1 = sh1->yes ; no1 = sh1->no ;
  yes2 = sh2->yes ; no2 = sh2->no ;
  suspect1 = sh1->suspect ; suspect2 = sh2->suspect ;
 
  if (!mergePolyA && 
      (
       sh1->end3p != sh2->end3p  ||
       sh1->begin3p != sh2->begin3p
       )
      )
    return 0 ;
      
  if (!s2mCompatibleShadows (sh1, sh2, bMax))
    return 0 ;
   /* test if 1 is included in 2 */
  i = bMax ; while (i--)
    if ( 
        (!bit(suspect1,i) && !bit(suspect2,i)) &&
        (
         (bit(yes1,i) && !bit(yes2,i)) ||
         (bit(no1,i)  && !bit(no2,i))
         )
        )
      { bad++ ; break ; } /* 1 is not included in 2 */

  /* test if 2 is included in 1 */
  if (bad)
    {
      i = bMax ; while (i--)
        if ( 
            (!bit(suspect1,i) && !bit(suspect2,i)) &&
            (
             (bit(yes2,i) && !bit(yes1,i)) ||
             (bit(no2,i)  && !bit(no1,i))
             )
            )
          return 0 ;
    }
  /* test if they touch by at least one exon */
  da = 0 ;
  ok = FALSE ; 
  i = bMax ; while (i--)
    if ( 
        (!bit(suspect1,i) && !bit(suspect2,i)) &&
        bit(yes2,i) && bit(yes1,i)
        )
      ok = TRUE ;
  if (ok)
    {
      a1 = sh1->a1 > sh2->a1 ? sh1->a1 : sh2->a1 ;
      a2 = sh1->a2 < sh2->a2 ? sh1->a2 : sh2->a2 ;
      if (a1 < a2) da = a2 - a1 ;
      else da = 0 ;
    }
  return da ;
} /* s2mIncludedShadows */

/**********************************************************************************/

static int shadowLengthOrder (const void *va, const void *vb)
{
  SH *sha = (SH *)va ;
  SH *shb = (SH *)vb ;
  return (sha->a2 - sha->a1) - (shb->a2 - shb->a1) ;
}  /* shadowLengthOrder */

/**********************************************************************************/

static void s2mMergeIncludedShadows (S2M *s2m, BOOL mergePolyA)
{
  int ii, jj, da, bestDa, bestJj ;
  SH *sh1, *sh2 ;
  Array shadows = s2m->shadows ;

  arraySort (shadows, shadowLengthOrder) ;
  /* merge included */
  for (ii = 0 ; ii + 1 < arrayMax(shadows) ; ii++)
    {
      sh1 = arrp (shadows, ii, SH) ;
      if (!sh1->clones)
        continue ;
      bestDa = 0 ; bestJj = ii ;
      for (jj = ii + 1 ; jj < arrayMax(shadows) ; jj++)
	{
	  sh2 = arrp (shadows, jj, SH) ;
	  if (!sh2->clones)
	    continue ;
	  
	  da = s2mIncludedShadows (sh1, sh2, s2m->bMax, mergePolyA) ;
	  if (da && da >= bestDa)
	    { bestDa = da ; bestJj = jj ; }
	}
      if (bestJj > ii)
	{
	  sh2 = arrp (shadows, bestJj, SH) ;
	  s2mDoMergeShadows (sh2, sh1, s2m->bMax) ;
	  s2mKillShadow (sh1) ;
	}
    }

  /* register happy few */
  s2mCompressShadows (shadows) ;

  return ; 
}  /* s2mMergeIncludedShadows */

/**********************************************************************************/

static void s2mVerifyContig (S2M *s2m, SC *sc)
{
  SH *sh ;
  int i, ii ;
  KEY clone ;
  BOOL debug = FALSE ;
  KEYSET allClones = 0 ;

 lao:
  allClones = keySetReCreate (allClones) ;
  if (debug)
    printf (" s2mVerifyContig\n") ;
  for (ii = 0 ; ii < arrayMax(sc->ks) ; ii++)
    {
      sh = arrp (s2m->shadows, keySet (sc->ks,ii), SH) ; 
      if (keySetExists (sh->clones))
        {
          if (debug) printf("v:%d: ", ii) ;
          for (i = 0 ; i < keySetMax (sh->clones)  ; i++)
            {
              clone = keySet (sh->clones, i) ;
              keySetInsert (allClones, clone) ;
              if (debug) printf(" %s ", name(clone)) ;
            }
          if (debug) printf("\n") ;
        }
    }
  if (debug) printf ("\n max(sc->sh.clones)=%d cumulatedMax=%d\n", keySetMax(sc->sh.clones), keySetMax(allClones)) ;
  if (keySetMax(sc->sh.clones) != keySetMax(allClones))
    {
      if (debug) messcrash ("wrong count in s2mVerifyContig") ;
      else { debug = TRUE ; goto lao ; }
    }

  keySetDestroy (allClones) ;

}

/**********************************************************************************/
  /* create contigs of shadows by moving left to right along the chain
   * and merging shadows with a common exon
   */

static void s2mMakeContigs (S2M *s2m)
{
  SH *sh ; SC *sc, *sc0 ;
  Array contigs = 0 ;
  Array shadows = s2m->shadows ;
  int ii, jj, jMax, sMax = arrayMax(shadows), bMax = s2m->bMax  ;
  BOOL debug = FALSE ;

  contigs = s2m->contigs = arrayHandleCreate (12, SC, s2m->h) ;

  if (debug) printf ("s2mmakecontigs max(s2m->shadows)=%d", sMax) ;
  for (ii = jj = 0 ; ii < sMax ; ii++)  /* happens because of dropEsts */
    {
      sh = arrp (shadows, ii, SH) ;
      if (!keySetExists (sh->clones) || 
          !keySetMax (sh->clones))
        continue ;
      if (jj < ii)
        arr (shadows, jj, SH) = arr (shadows, ii, SH) ;
      jj++ ;
    }
  sMax = arrayMax(shadows) = jj ; 
  if (debug) printf("  dropping to max(s2m->shadows)=%d", sMax) ;
  if (debug) printf ("make contigs: %d shadows ", jj) ;
  for (ii = 0, jMax = 0 ; ii < sMax ; ii++)  /* for each shadow */
    {
      sh = arrp (shadows, ii, SH) ;
      sc0 = 0 ;
      for (jj = 0 ; jj < jMax ; jj++)  /* try to merge in every scontig */
        {
          sc = arrp (contigs, jj, SC) ;
          if (!sc->sh.clones)
            continue ;
          if (s2mConnectedShadows (&(sc->sh), sh, bMax))
            {
              if (!sc0) /* merge shadow in this first contig */
                {
                  s2mDoMergeShadows (&(sc->sh), sh, bMax) ; /* no kill */
                  keySet (sc->ks, keySetMax(sc->ks)) = ii ;
                  sc0 = sc ;
                }
              else  /* merge the second contig to the first good one */
                {
                  s2mDoMergeContigs (sc0, sc, bMax) ;
                  s2mKillContig (sc) ;
                }                  
            }
        }
      if (!sc0) /* create a new contig */
        {
          sc = arrayp (contigs, jMax++, SC) ;
          sc->s2m = s2m ;
          sc->isUp = sh->isUp ;
          sc->sh.clones = arrayHandleCreate (12, KEY,  s2m->h) ;
          sc->ks = arrayHandleCreate (12, KEY, s2m->h) ;
          keySet (sc->ks, 0) = ii ;
          sc->sh.bMax = s2m->bMax ;
          sc->sh.yes = bitSetCreate (bMax, s2m->h) ;
          sc->sh.no = bitSetCreate (bMax, s2m->h) ;
          sc->sh.ghost = bitSetCreate (bMax, s2m->h) ;
          sc->sh.suspect = bitSetCreate (bMax, s2m->h) ;
          sc->sh.a1 = sh->a1 ;
          sc->sh.a1 = sh->a2 ;
          s2mDoMergeShadows (&(sc->sh), sh, bMax) ;
          s2mVerifyContig (s2m, sc) ;
        }
      else 
        s2mVerifyContig (s2m, sc0) ;
    }

  /* register happy few */
  s2mCompressContigs (contigs) ;
  if (debug) printf( "  %d contigs\n", arrayMax(contigs)) ;
}  /* s2mMakeContigs */

/**********************************************************************************/
/**********************************************************************************/

static BOOL scCompatible (void *vp, int i1, int i2)
{
  SC *sc = (SC*) vp ;
  SH *sh1, *sh2 ;
  BOOL ok, debug = FALSE ;

  sh1 = arrp (sc->s2m->shadows, keySet(sc->ks,i1), SH) ;
  sh2 = arrp (sc->s2m->shadows, keySet(sc->ks,i2), SH) ;
  
  if (sh1->clones && sh2->clones)
    {
      ok = s2mCompatibleShadows (sh1, sh2, sc->s2m->bMax) ;
      if (debug)
        printf ("scCompatible %s %s %s\n", ok ? "yes" : "no ",
                name(keySet(sc->sh.clones,i1)), name(keySet(sc->sh.clones,i2))) ;
      return ok ;
    }
  else
    return FALSE ;
}

/**********************************************************************************/

static BOOL scConnected (void *vp, int i1, int i2)
{
  SC *sc = (SC*) vp ;
  SH *sh1, *sh2 ;
  BOOL ok, debug = FALSE ;

  sh1 = arrp (sc->s2m->shadows, keySet(sc->ks,i1), SH) ;
  sh2 = arrp (sc->s2m->shadows, keySet(sc->ks,i2), SH) ;
  
  if (sh1->clones && sh2->clones)
    {
      ok = s2mConnectedShadows (sh1, sh2, sc->s2m->bMax) ;
      if (debug)
        printf ("scConnected %s %s %s\n", ok ? "yes" : "no ",
                name(keySet(sc->sh.clones,i1)), name(keySet(sc->sh.clones,i2))) ;
      return ok ;
    }
  else
    return FALSE ;
}

/**********************************************************************************/

static Array s2mDiamino (S2M *s2m, SC *sc, int *diamMax)
{
  int i, ii, scMax = arrayMax(sc->ks) ;
  int bMax = s2m->bMax ;
  Array diams = diaminoCreate (sc, arrayMax(sc->ks), scCompatible, scConnected, *diamMax) ;
  DIAMWORD *mm ;
  Array mrnas = 0 ;
  SC *mrna ; SH *sh2 ;
  BOOL debug = FALSE ;

  mrnas = arrayHandleCreate (arrayMax(diams), SC, s2m->h) ;
  if ( *diamMax && arrayMax(diams) >= *diamMax)
    {
      arrayDestroy (mrnas) ;
      messout ("error: s2mDiamino created %d > %d diams with %d clones in cosmid %s",
                 arrayMax(diams), *diamMax, keySetMax(sc->ks), name(s2m->cosmid)) ;
     }
  else
    for (ii = 0 ; ii < arrayMax(diams) ; ii++)
      {
        mrna = arrayp (mrnas, ii, SC) ;
        mrna->s2m = s2m ;
        mrna->ks = arrayHandleCreate (20, KEY, s2m->h) ;
        mrna->sh.clones = arrayHandleCreate (20, KEY, s2m->h) ;
        mrna->sh.yes = bitSetCreate (bMax, s2m->h) ;
        mrna->sh.no = bitSetCreate (bMax, s2m->h) ;
        mrna->sh.ghost = bitSetCreate (bMax, s2m->h) ;
        mrna->sh.suspect = bitSetCreate (bMax, s2m->h) ;
        mrna->sh.bMax = s2m->bMax ;
        
        mm = arrp (diams, ii, DIAMWORD) ;
        for (i = 0 ; i < scMax ; i++)
          {
            if (bit (mm->bb, i))
              {
                sh2 = arrp (s2m->shadows, keySet (sc->ks, i), SH) ;
                if (sh2->clones)
                  s2mDoMergeShadows (&(mrna->sh) , sh2, bMax) ;
              }
          }
#ifdef DEBUG
        /* verify */
        for (i = 0 ; i < bMax ; i++)
          if (bit (mrna->sh.yes, i) && bit (mrna->sh.no, i))
            messerror ("error in s2mdiamino") ;
        
        for (jj = ii+1 ; jj < arrayMax(diams) ; jj++)
          {
            m2 = arrp (diams, jj, DIAMWORD) ;
            for (i = 0 ; i < scMax ; i++)
              if (bit (mm->bb, i) != bit (m2->bb, i))
                break ;
            if (i >= scMax)
              messerror ("diam %d %d are identical", ii, jj) ;
          }
#endif        
      }

  *diamMax = arrayMax(diams) ;
  diaminoDestroy (diams) ;
      
  if (debug && mrnas)
    s2mShowContigs (mrnas, "mrnas obtained via diamino", -1) ;
  return mrnas ;
}

/**********************************************************************************/
/* each transcript is described by a HIT array */
static void mrnaFindAlternatives (SC *sc, Array mrnaShadows, Array allYes, Array allNo)
{
  int *np, ii, jj, type ;
  int bMax = sc->sh.bMax, mMax = arrayMax (mrnaShadows) ;
  BitSet yes, no ;

  /* count how many mrna use bit ii as exon intron */
  for (jj = 0 ; jj < arrayMax(mrnaShadows) ; jj++)
    {
      yes = arr (mrnaShadows, jj, SC).sh.yes ;
      no = arr (mrnaShadows, jj, SC).sh.no ;
      for (ii = 0 ; ii < bMax ; ii++)
        {
          if (bit(yes,ii)) arr (allYes, ii, int)++ ;
          if (bit(no,ii)) arr (allNo, ii, int)++ ;
        }
    }
  /* any contiguous set of exon intron not used maximally is alternative */
  for (ii = 0 ; ii < bMax ; ii++)
    if (arr (allYes, ii, int)) /* entering an exon */
      {
        type = 0 ;
        for (jj = ii, np = arrp (allYes, jj, int) ; jj < bMax && *np ; np++, jj++)
          if (*np < mMax) { type = 1 ; break ; }
        if (type)
          for (jj = ii, np = arrp (allYes, jj, int) ; jj < bMax && *np ; np++, jj++)
            *np = -1 ;
        ii = jj ; /* position ii on a !bit */
      }
  for (ii = 0 ; ii < bMax ; ii++)
    if (arr (allNo, ii, int)) /* entering an intron */
      {
        type = 0 ;
        for (jj = ii, np = arrp (allNo, jj, int) ; jj < bMax && *np ; np++, jj++)
          if (*np < mMax) { type = 1 ; break ; }
        if (type)
          for (jj = ii, np = arrp (allNo, jj, int) ; jj < bMax && *np ; np++, jj++)
            *np = -1 ;
        ii = jj ; /* position ii on a !bit */
      }
}

/**********************************************************************************/

static void s2mMakeOneTranscript (S2M *s2m, SC *sc, SC *mm, SMRNA *smrna, Array allYes, Array allNo)
{
  HIT *up = 0 ;
  int ii, iMrna , b1, b2, type, oldType,  isAlternative ;
  int bMax = s2m->bMax ;
  Array chain = s2m->chain ;
  BitSet yes, no ;
  CHAIN *c1, *c2 ;  /* mieg 2023_11_23 */

  smrna->clones = arrayHandleCreate (keySetMax(mm->sh.clones), KEY, s2m->h) ;
  for (ii = 0 ; ii < keySetMax(mm->sh.clones) ; ii++)
    keySet (smrna->clones, ii) = keySet(mm->sh.clones, ii) ;

  /* find first and last used bit */
  yes = mm->sh.yes, no = mm->sh.no ;
  b1 = -1 ; b2 = -2 ;
  for (ii = 0 ; ii < bMax ; ii++)
    if (bit (yes,ii))
      {
        if (b1 == -1) b1 = ii ;
        b2 = ii ;
      }
   
  /* report all used exon intron gap */
  if (0) s2mShowOneShadow (100, &(mm->sh), 0) ;
  for (iMrna = oldType = 0, c1 = arrp (chain, 0, CHAIN),  ii = 1 ; ii < bMax ; ii++)
    {
      c2 = arrp (chain, ii, CHAIN) ;
      type = isAlternative = 0 ;
      if (bit (yes,ii))
        {
          type = gX ;
          if (arr (allYes, ii, int) < 0) type |= gB ;            
        }
      else if (ii >= b1 && ii <= b2 && bit (no, ii))
        {
          type = gI ;
          if (arr (allNo, ii, int) < 0) type |= gB ;
        }              
      else if (ii > b1 && ii < b2)
        type = gGap ;
      if (type)
        {
          if (type & oldType & (gX | gI | gGap))
            {
              up->a2 = c2->x ; /* extend */
              up->type |= type ;
            }
          else
            {
              up = arrayp (smrna->hits, iMrna++, HIT) ;
              up->a1 = c1->x + 1 ; up->a2 = c2->x ;
              up->type = type ;
            }
        }
      oldType = type ;
      c1 = c2 ;
    }  
  
}  /* mrnaMakeOneTranscript */

/**********************************************************************************/
/* import detailled flags */
/* mieg feb12 2003, was s2m->geneHits just before
 * problem is that geneHits does not have all the clone at any give intron 
 * but plainHits is missing some flagging 
*/

static BOOL s2mFlagOneTranscript (Array clones, Array sHits, Array plainHits, Array gHits)
{
  int ii, jj, jj2, lastjj2 ;
  HIT *up, *vp, *vp2 ;
  Array gHitsSorted = 0 ;
  Array plainHitsSorted = 0 ;
  BOOL found = FALSE ;

  if (!arrayMax(sHits) || !arrayMax(gHits)) return FALSE ;

  gHitsSorted = arrayCopy(gHits) ;
  cDNASwapA (gHitsSorted) ;
  arraySort (gHitsSorted, cDNAOrderGloballyByA1) ;

  plainHitsSorted = arrayCopy(plainHits) ;
  cDNASwapA (plainHitsSorted) ;
  arraySort (plainHitsSorted, cDNAOrderGloballyByA1) ;

  up = arrp (sHits, 0, HIT) ;
  for (ii = 0 ; ii < arrayMax (sHits) ; ii++, up++)
    {
      if (up->type & gX)
        {
	  up->type &= ~gReal5p ; /* tilda NOT */
	  up->type &= ~gReal3p ; /* tilda NOT */
          up->type &= ~gA ; /* tilda NOT */
          up->type &= ~gS ; /* tilda NOT */
          up->type &= ~gS0 ; /* tilda NOT */
          up->type &= ~gCompleteCDS ; /* tilda NOT */

          for (jj2 = 0, vp2 = arrp (plainHitsSorted, jj2, HIT); 
               jj2 < arrayMax(plainHitsSorted) ; vp2++, jj2++)
	    {
	      if (keySetFind (clones, vp2->cDNA_clone, 0))
		{
		  if ((vp2->type & gReal3p) &&
		      vp2->a1 < up->a2 && vp2->a2 > up->a1)
		    up->type |= gReal3p ; 
		  if ((vp2->type & gA) &&
		      vp2->a1 < up->a2 && vp2->a2 > up->a1)
		    up->type |= gA ;
		  if ((vp2->type & gReal5p) &&
		      vp2->a1 < up->a2 && vp2->a2 > up->a1)
		    up->type |= gReal5p ;
		  if ((vp2->type & gCompleteCDS) &&
		      vp2->a1 < up->a1+60 && vp2->a2 > up->a1)
		    up->type |= gCompleteCDS ; 
		  if (vp2->reverse)
		    {
		      if ((vp2->type & gS0) &&
			  vp2->a2 > up->a2 - 60 && vp2->a1 < up->a1)
			up->type |= gS0 ; 
		      if ((vp2->type & gS) &&
			  vp2->a2 == up->a2)
			up->type |= gS ; 
		    }
		  else
		    {
		      if ((vp2->type & gS0) &&
			  vp2->a1 < up->a1+60 && vp2->a2 > up->a1)
			up->type |= gS0 ; 
		      if ((vp2->type & gS) &&
			  vp2->a1 == up->a1)
			up->type |= gS ; 
		    }
		}
            }
        }
      if (up->type & gI) 
        {
          up->type |= gJ ; /* tilda NOT */
          for (jj2 = 0, vp2 = arrp (plainHitsSorted, jj2, HIT); 
               jj2 < arrayMax(plainHitsSorted) ; vp2++, jj2++)
            {
              if ((vp2->type & gI) && !(vp2->type & gJ) &&
                  vp2->a1 == up->a1 && vp2->a2 == up->a2)
                { up->type &= ~gJ ; break ; }
            }
        }
    }

  if (0) /* obsolete 2009_06_17 */
    {
      up = arrp (sHits, 0, HIT) ;
      vp = arrp (gHitsSorted, 0, HIT) ;
      lastjj2 = 0 ;
      for (ii = jj = 0 ; ii < arrayMax (sHits) && jj < arrayMax(gHitsSorted) ; )
	{
	  if (!class (vp->est))  { vp++ ; jj++ ; continue ;} /* ghost */
	  if (!keySetFind (clones, vp->cDNA_clone, 0))
	    { /* try to find the same zone in plainhitssorted */
	      /* we must restart each time on vp->a1 */
	      found = FALSE ;
	      if (lastjj2 < arrayMax(plainHitsSorted))
		{
		  for (jj2 = lastjj2,  vp2 = arrp (plainHitsSorted, jj2, HIT) ;
		       vp2->a1 < vp->a1 && jj2 < arrayMax(plainHitsSorted) ;
		       vp2++, jj2++) ;
		  lastjj2 = jj2 ;
		  for ( ; jj2 < arrayMax(plainHitsSorted) &&
			  vp2->a1 == vp->a1 && !found ; vp2++, jj2++)
		    if (keySetFind (clones, vp2->cDNA_clone, 0))
		      found = TRUE ;
		}
	      if (!found)
		{ vp++ ; jj++ ; continue ; }
	    }          
	  if (up->a1 < vp->a1) 
	    {
	      if (up->type & gI) up->type |= gJ ; /* since i did not find its justification */
	      up++ ; ii++ ; continue ;
	    }
	  if (up->a1 > vp->a1 || up->a2 > vp->a2) { vp++ ; jj++ ; continue ;}
	  if (up->type & gI) 
	    {
	      if (up->a2 == vp->a2 && (vp->type & gI)) 
		{ if (vp->type & gJ) up->type |= gJ ; }
	      else
		up->type |= gJ ; /* since i did not find its justification */
	      up++ ; ii++ ; vp++ ; jj++ ;
	      continue ;
	    }
	  /* up++ ; ii++ ; */  vp++ ; ; jj++ ;    
	} 
    }

  arrayDestroy (gHitsSorted) ;
  arrayDestroy (plainHitsSorted) ;

  up = arrp (sHits, 0, HIT) ; 
  for (ii = 0 ; ii < arrayMax (sHits) ; up++, ii++)
    if (up->type & gX)
      return TRUE ;
  return FALSE ;
}

/**********************************************************************************/
/* make an array of transcripts */
static Array s2mMakeTranscripts (S2M *s2m, SC *sc, Array mrnaShadows) 
{
  Array smrnas = arrayHandleCreate (12, SMRNA, s2m->h) ;
  SMRNA *smrna ;
  int ii, jj ;
  Array allYes = 0, allNo = 0 ;
  SC *mm ;

  allYes = arrayCreate (s2m->bMax, int) ; array (allYes, s2m->bMax -1, int) = 0 ; 
  allNo = arrayCreate (s2m->bMax, int) ; array (allNo, s2m->bMax -1, int) = 0 ; 
  mrnaFindAlternatives (sc, mrnaShadows, allYes, allNo) ;

  for (ii = jj = 0 ; ii < arrayMax(mrnaShadows) ; ii++)
    { 
      mm = arrp (mrnaShadows, ii, SC) ;
      if (!keySetMax(mm->sh.clones))
        continue ; /* may happen after dropEst */
      smrna = arrayp (smrnas, jj++, SMRNA) ;
      smrna->hits = arrayHandleCreate (s2m->bMax, HIT, s2m->h) ;

      s2mMakeOneTranscript (s2m, sc, mm, smrna, allYes, allNo) ;
      if (!s2mFlagOneTranscript (smrna->clones, smrna->hits, 
                                 s2m->plainHits, s2m->geneHits) ||
          !arrayMax(smrna->hits))
        { arrayDestroy (smrna->hits) ; jj-- ; arrayMax(smrnas)-- ; continue ; }

      if (0)
        {
          printf ("mrnaMakeTranscripts: ii = %d\n", ii) ;
          showHits(smrna->hits) ;
        }
    }
  arrayDestroy (allYes) ;
  arrayDestroy (allNo) ;
  return smrnas ;
}

/*********************************************************************/

static BOOL mrnaFilterGene (SC *sc)
{
  int g1, g2, x = 0 ;
  HIT *hh ;
  Array linkPos = sc->s2m->linkPos ;

  g1 = sc->a1 ; g2 = sc->a2 ;

  /* keep only genes which are more left than the begin of cosmid2 */
  if (linkPos && arrayMax(linkPos))
    {
      hh = arrp(linkPos, arrayMax(linkPos) - 1, HIT) ;
      x = hh->x1  < hh->x2 ? hh->x1 : hh->x2 ;
      if (g1 >= x && g2 >= x)
        return FALSE ;  /* drop this gene */
    }
  return TRUE ;
}

/**********************************************************************************/

static SMRNA *mrnaMakeGene (S2M *s2m, SC *sc, Array mrnaShadows, Array smrnas)
{
  int ii, jj, ig,  g1, g2 ;
  HIT *up, *vp ;
  SMRNA *smrna, *gmrna ;

  gmrna = arrayp (s2m->gmrnas, arrayMax(s2m->gmrnas), SMRNA) ;
  gmrna->hits = arrayHandleCreate (20, HIT, s2m->h) ;

  /* find the global origin of the gene */
  g1 = g2 = -1 ;
  for (ii = 0; ii < arrayMax(smrnas) ; ii++)
    {
      smrna = arrp (smrnas, ii, SMRNA) ;
      for (jj = 0 ; jj < arrayMax(smrna->hits) ; jj++)
        {
          up = arrp (smrna->hits, jj, HIT) ;
          if (!ii && !jj)
            { 
              g1 = up->a1 < up->a2 ? up->a1 : up->a2 ;
              g2 = up->a1 > up->a2 ? up->a1 : up->a2 ;
            }
          if (g1 > up->a1) g1 = up->a1 ;
          if (g1 > up->a2) g1 = up->a2 ;
          if (g2 < up->a1) g2 = up->a1 ;
          if (g2 < up->a2) g2 = up->a2 ;
        }
      /*
        printf ("mrnaMakeGene original: ii = %d\n", ii) ;
        showHits(mrna) ;
      */
    }

  /* reset the origins and accumulate in gene */
  if (sc->isUp)
    { int g0 = g1 ; g1 = g2 ; g2 = g0 ; }
  sc->a1 = g1 ; sc->a2 = g2 ; 

  for (ii = ig = 0; ii < arrayMax(smrnas) ; ii++)
    {
      smrna = arrp (smrnas, ii, SMRNA) ;
      
      for (jj = 0 ; jj < arrayMax(smrna->hits) ; jj++)
        {
          up = arrp (smrna->hits, jj, HIT) ;
          if (!sc->isUp)
            {
              up->reverse = FALSE ;
              up->a1 = up->a1 - g1 + 1 ;
              up->a2 = up->a2 - g1 + 1 ;
            }
          else
            { 
              up->reverse = TRUE ;
              up->a1 = g1 - up->a1 + 1 ;
              up->a2 = g1 - up->a2 + 1 ;
            }
          
          vp = arrayp (gmrna->hits, ig++, HIT) ;
          *vp = *up ;
        }
      cDNASwapA (smrna->hits) ; 
      arraySort (smrna->hits, cDNAOrderByA1) ;
      smrna->a1 = arrp (smrna->hits, 0, HIT)->a1 ;
      smrna->a2 = arrp (smrna->hits, arrayMax(smrna->hits) -1, HIT)->a2 ;
      for (jj = 0 ; jj < arrayMax(smrna->hits) ; jj++)
        {
          up = arrp (smrna->hits, jj, HIT) ;
          up->a1 = up->a1 - smrna->a1 + 1 ;
          up->a2 = up->a2 - smrna->a1 + 1 ;
        }
    }
  cDNASwapA (gmrna->hits) ; 
  arraySort (gmrna->hits, cDNAOrderByA1) ;
  arrayCompress (gmrna->hits) ; 
  
  for (ii = jj = 0 ; ii < arrayMax(gmrna->hits) - 1 ; ii++)
    {
      up = arrp (gmrna->hits, ii, HIT) ;
      vp = up + 1 ;
      if (up->type && vp->type && up->a1 == vp->a1 && up->a2 == vp->a2)
        {
          if ((up->type & gGap) && !(vp->type & gGap))
            up->type = 0 ;
          if ((up->type & gJ) && !(vp->type & gJ))
            up->type = 0 ;
          if (!(up->type & gGap) && (vp->type & gGap))
            vp->type = 0 ;
          if (!(up->type & gJ) && (vp->type & gJ))
            vp->type = 0 ;
        }
    }
  /* keep happy few */
  for (ii = jj = 0 ; ii < arrayMax(gmrna->hits) ; ii++)
    {
      up = arrp (gmrna->hits, ii, HIT) ;
      if (!up->type) continue ;
      if (jj < ii)
        {
          vp = arrp (gmrna->hits, jj, HIT) ;
          *vp = *up ;
        }
      jj++ ;
    }
  arrayMax(gmrna->hits) = ii ;
 
  return gmrna ;
}

/**********************************************************************************/

static void mrnaMakeEstHits (S2M *s2m, SC *sc, SMRNA *gmrna) 
{
  HIT *up, *vp ;
  int ii, jj, g1 = sc->a1, amin, amax ; 
  Array allHits, hits ;
  KEYSET clones = sc->sh.clones ;

  allHits = s2m->plainHits ;
  gmrna->estHits = hits = arrayHandleCreate (12 * keySetMax (clones) , HIT, s2m->h) ;
  cDNASwapA (allHits) ;
  
  /* normal clones */
  amin = ACEDB_MAXINT ;
  amax = -1 ;
  for (ii = jj = 0, up = arrp (allHits, ii, HIT) ; ii < arrayMax(allHits) ; up++, ii++)
    {
      if (keySetFind (clones, up->cDNA_clone, 0))
        {
          vp = arrayp (hits, jj++, HIT) ;
          *vp = *up ;
          if (up->a1 < amin) amin = up->a1 ;
          if (up->a2 > amax) amax = up->a2 ;
          if (sc->isUp)
            {
              vp->a1 = g1 - vp->a1 + 1 ;
              vp->a2 = g1 - vp->a2 + 1 ;
            }
          else
            {
              vp->a1 = vp->a1 - g1 + 1 ;
              vp->a2 = vp->a2 - g1 + 1 ;
            }
        }
    } 
  /* now recover the tracking clones */
  for (ii = 0, up = arrp (allHits, ii, HIT) ; ii < arrayMax(allHits) ; up++, ii++)
    {
      if ((up->type & gDroppedGene) && 
          up->a1 < amax - 30 && up->a2 > amin + 30)
        {
          keySetInsert (clones, up->cDNA_clone) ;
          vp = arrayp (hits, jj++, HIT) ;
          *vp = *up ;
          if (sc->isUp)
            {
              vp->a1 = g1 - vp->a1 + 1 ;
              vp->a2 = g1 - vp->a2 + 1 ;
            }
          else
            {
              vp->a1 = vp->a1 - g1 + 1 ;
              vp->a2 = vp->a2 - g1 + 1 ;
            }
        }
    } 
  cDNASwapA (hits) ;
  arraySort (hits, cDNAOrderByA1) ;
}

/**********************************************************************************/
#ifdef JUNK
static void s2mReviseDebutExon (S2M *s2m, SC *sc, int i0)
{
  int ii = 0, a1 ;
  HIT *hh, *h0 ;
  KEY est ;

  h0 = arrp (s2m->plainHits, i0, HIT) ;
  a1 = h0->a1 ; est = h0->est ;
  for (ii = 0 ; ii < arrayMax(s2m->plainHits) ; ii++)
    {
      hh = arrp (s2m->plainHits, ii, HIT) ;
      if (hh->type & gDroppedGene)
        continue ;

      if (! (hh->type & gX) || !(hh->type & gD))
        continue ; 
      if (est == hh->est)
        continue ;

      if (hh->a1 < a1 - 10 || hh->a1 > a1 + 10)
        continue ;
      
    }
}

/**********************************************************************************/

static void s2mReviseFinExon (S2M *s2m, SC *sc, int i0)
{
  int ii = 0, a2;
  HIT *hh, *h0 ;
  KEY est ;
 
  h0 = arrp (s2m->plainHits, i0, HIT) ;
  a2 = h0->a2 ; est = h0->est ;

  for (ii = 0 ; ii < arrayMax(s2m->plainHits) ; ii++)
    {
      hh = arrp (s2m->plainHits, ii, HIT) ;
      if (hh->type & gDroppedGene)
        continue ;

      if (! (hh->type & gX) || !(hh->type & gF))
        continue ;
      if (est == hh->est)
        continue ;

      if (hh->a2 < a2 - 10 || hh->a2 > a2 + 10)
        continue ;
    }
}

/**********************************************************************************/

static void s2mReviseJumps (S2M *s2m, SC *sc)
{
  int ii = 0 ;
  HIT *hh ;
  BOOL debug = FALSE ;

  if (debug)
    showHits (s2m->plainHits) ;

  for (ii = 0 ; ii < arrayMax(s2m->plainHits) ; ii++)
    {
      hh = arrp (s2m->plainHits, ii, HIT) ;
      if (hh->type & gDroppedGene)
        continue ;
      if (! (hh->type & gX))
        continue ;
      if ( ! (hh->type & gD))
         s2mReviseDebutExon (s2m, sc, ii) ;
      if ( ! (hh->type & gF))
         s2mReviseFinExon (s2m, sc, ii) ;
    }
  return ;
}
#endif
/**********************************************************************************/

static void s2mDropEsts (S2M *s2m, SC *sc, int contg, int bestQual, KEYSET clonesWithBadIntron)
{
  int ii, jj, i, j ;
  KEY clone, est ;
  SH *sh ;
  BOOL debug = FALSE, noMrna = TRUE, isBad ;
  Array shadows = s2m->shadows ;
  KEY _Is_est = str2tag ("Is_est") ;
  KEYSET ks = 0, allClones = keySetCreate (), goodClones = keySetCreate () ;

  if (0 && bestQual <= 2) debug=TRUE ;
  if (debug)
    printf ("\n// s2mDropEsts contig %d qual %d max(ks) = %d ",  contg, bestQual, arrayMax(sc->ks)) ;
  
  for (ii = 0 ; ii < arrayMax(sc->ks) ; ii++)
    {
      sh = arrp (shadows, keySet (sc->ks,ii), SH) ; 
      if (keySetExists (sh->clones))
        {
          if (debug) printf("\n v:%d: ", ii) ;
          for (i = 0 ; i < keySetMax (sh->clones)  ; i++)
            {
              clone = keySet (sh->clones, i) ;
              keySetInsert (allClones, clone) ;
              if (debug) printf(" %s ", name(clone)) ;
            }
          if (debug) printf("\n") ;
        }
    }
  if (debug) printf ("\n max(sc->ks)=%d cumulatedMax=%d\n", keySetMax(sc->sh.clones), keySetMax(allClones)) ;
  
  
  for (ii = jj = 0 ; ii < arrayMax(allClones) ; ii++)
    {
      clone = keySet (allClones, ii) ;
      if (keySetFind (clonesWithBadIntron, clone, 0))
        isBad = TRUE ;
      else
        isBad = FALSE ;
      if (debug) printf(" ??:%s ", name(clone)) ;
      if ((est = keyGetKey (clone, _Read)))
        {
          /* always keep the NM */
          if (keyFindTag (est, _Ref_Seq) && !keyFindTag (est, _Is_est))
            { 
              noMrna = FALSE ; keySet (goodClones, jj++) = clone ;
              if (debug) printf(" ok1:%s ", name(clone)) ;
            }
          /* keep mrna >= 3 if bestQual > 0 */
          else if (isBad)
            { 
              ks = queryKey (est
                             , messprintf("(Ref_mRNA && !IS_est && From_cosmid:9 <= %d) || "
                                          "(!Ref_mRNA && !IS_est && From_cosmid:9 <= %d)"
                                          , bestQual > 3 ? bestQual - 2 : 2
                                          , bestQual > 4 ? bestQual - 3 : 1 
                                          )) ;
              if (keySetMax (ks))
                { 
                  keySet (goodClones, jj++) = clone ;
                  if (debug) printf(" ok2:%s ", name(clone)) ;
                }
              else
                {
                  if (debug) printf(" bad2 ") ;
                }
              keySetDestroy (ks) ;
            }
          /* keep mrna >= 3 if bestQual > 0 */
          else if (bestQual > 0 && keyFindTag (est, _Ref_mRNA) && !keyFindTag (est, _Is_est))
            { 
              ks = queryKey (est
                             , messprintf("(! other && From_cosmid:9 <= %d) ||"
                                          " (other && From_cosmid:9 <= %d)"
                                          , bestQual > 3 ? bestQual : 3
                                          , bestQual > 3 ? bestQual - 1 : 2 
                                          )) ;
              if (keySetMax (ks))
                { 
                  noMrna = FALSE ;
                  keySet (goodClones, jj++) = clone ;
                  if (debug) printf(" ok3:%s ", name(clone)) ;
                }
              else
                {
                  if (debug) printf(" bad3 ") ;
                }
              keySetDestroy (ks) ;
            }
	  else if (keyFindTag (est, _Composite))
	    {
	      int n = 1 << (2 * (11 - bestQual)) ;
	      if (n < 2) n = 1000000 ;
              ks = queryKey (est
			     , messprintf ("Tags > %d || Composite > %d || COUNT intron > 2"
					   , n, n ) 
			     ) ;
              if (keySetMax (ks))
                { 
                  noMrna = FALSE ;
                  keySet (goodClones, jj++) = clone ;
                  if (debug) printf(" ok3:%s ", name(clone)) ;
                }
              else
                {
                  if (debug) printf(" bad30 ") ;
                }
              keySetDestroy (ks) ;


	    }
          /* keep est if bestQual > 0 */
          else if (bestQual > 0)
            { 
              ks = queryKey (est
                             , messprintf("(! other && gt_ag && From_cosmid:9 <= %d) ||"
                                          "(! other && From_cosmid:9 <= %d) ||"
                                          " (other && From_cosmid:9 <= %d)"
                                          , bestQual
                                          , bestQual - 2
                                          , bestQual - 3
                                          )) ;
              if (keySetMax (ks))
                { 
                  keySet (goodClones, jj++) = clone ;
                  if (debug) printf(" ok4:%s ", name(clone)) ;
                }
              else
                {
                  if (debug) printf(" bad4 ") ;
		  /*
                  if (debug && !keySetMax (ks) && !strcasecmp(name(est),"CA310331"))
                    {
                      OBJ Est1 = bsCreate(est) ;
                      BOOL b1, b2, b3 ;
                      
                      b1 = bsFindTag (Est1, _Other) ;
                      b2 = bsFindTag (Est1, _gt_ag) ;
                      b3 = bsFindTag (Est1, str2tag("From_cosmid")) ;
                      invokeDebugger() ;
                      
                      bsDestroy (Est1) ;
                    }
		  */
                }
              keySetDestroy (ks) ;
            }
          /* bestQual <= 0, just keep the mrna, mrna est are already quality 3, est 1 */
          else if (keyFindTag (est, _Ref_mRNA) && !keyFindTag (est, _Is_est))
            { 
              ks = queryKey (est
                             , messprintf("(! other && From_cosmid:9 <= %d) ||"
                                          " (other && From_cosmid:9 <= %d)"
                                          , bestQual + 2
                                          , bestQual + 1
                                          )) ;
              if (keySetMax (ks))
                { 
                  noMrna = FALSE ;
                  keySet (goodClones, jj++) = clone ;
                  if (debug) printf(" ok4:%s ", name(clone)) ;
                }
              keySetDestroy (ks) ;
            }
        }
    }
  if (debug) printf ("\ns2mDropEsts %d clones %d goodClones nomRNA=%d\n",   
                     arrayMax(sc->sh.clones), keySetMax(goodClones), noMrna) ;
  for (ii = jj = 0 ; ii < arrayMax(sc->ks) ; ii++)
    {
      j = 0 ;
      sh = arrp (shadows, keySet (sc->ks,ii), SH) ; 
      if (keySetExists (sh->clones))
        {
          BOOL ok = FALSE ;
          for (i = j = 0 ; i < keySetMax (sh->clones)  ; i++)
            {
              clone = keySet (sh->clones, i) ;
              if (keySetFind (goodClones, clone, 0))
                ok = TRUE ;
            }
                        
          for (i = j = 0 ; i < keySetMax (sh->clones)  ; i++)
            {
              clone = keySet (sh->clones, i) ;
              if (0 && !strcmp(name(clone),"AF052793"))
                invokeDebugger () ;
              
              if (debug) printf("testing ii=%d qual=%d i=%d %s ",  ii, bestQual, i, name(clone)) ;
              if (keySetFind (goodClones, clone, 0))
                { 
                  if (debug)
                    printf (" keeping %s \n",  name(clone)) ;
                  keySet (sh->clones, j++) = clone ;
                }
              else if (ok)
                { 
                  if (debug)
                    printf (" keeping also %s \n",  name(clone)) ;
                  keySet (sh->clones, j++) = clone ;
                }
              else
                {
                  OBJ Clone = 0 ;
                  if (debug)
                    printf (" dropping %s \n",  name(clone)) ;
                  if ((Clone = bsUpdate (clone)))
                    {
                      bsAddTag (Clone, str2tag ("Ignored_in_gene")) ;
                      bsSave (Clone) ;
                    }
                }
            }
          arrayMax (sh->clones) = j ;
        }
    }
  
  keySetDestroy (allClones) ;
  keySetDestroy (goodClones) ;
  return ;
} /* s2mDropEsts */

/**********************************************************************************/
/**********************************************************************************/
#define MAX_SHADOW 160
KEYSET mrnaCreate (int type, KEY cosmid, KEY cosmid1,
                   Array dnaD, Array dnaR, Array plainHits, Array geneHits, 
                   Array linkPos, KEYSET clipTops, KEYSET clipEnds,
                   KEYSET clonesWithBadIntron)
{
  S2M *s2m = s2mCreate (type, cosmid, cosmid1, dnaD, dnaR, plainHits, geneHits, linkPos) ;
  int ii, ig = 0, diamMax = 0, DIAMMAX = 1800, bestQual ;
  Array smrnas = 0 ;
  SMRNA *gmrna = 0 ;
  KEYSET genes = keySetCreate () ;
  BOOL debug = FALSE, debug1 = FALSE ;
  static int nCall22 = 0 ; /* for debugging */
  extern int mrnaPleaseNoAbandon(void) ;
  BOOL isDropped = FALSE ;
  BOOL hasReadPairs = FALSE, hasFusedClones = FALSE ;
  if (debug1)
    {
      freeOutf ("Entering mrnaCreate %s nCall = %d\n", name(cosmid), nCall22) ;
      freeOutf ("MC:%s ", timeShowNow()) ;
    }
  chrono ("mrnaCreate") ;
  s2mMakeChain (s2m, debug) ;
    /* make shadow sets with a single clone per shadow */ 
  s2mMakeShadows (s2m, debug, &hasReadPairs, &hasFusedClones ) ;  
  if (debug) s2mShowShadows (s2m->shadows, "make-shadows") ;
  if (debug1)  freeOutf ("mrnaCreate: %d reads = initial shadows\n", arrayMax (s2m->shadows)) ;
  if (hasReadPairs)
    {
      s2mMergeCloneShadows(s2m) ;     /* merge shadows from single clone */
      if (debug) s2mShowShadows (s2m->shadows, "clone-shadows") ;
      if (debug1)  freeOutf ("mrnaCreate: %d clone-shadows\n", arrayMax (s2m->shadows)) ;
    }
  if (hasFusedClones) /* already taken care of by the ghost system, this would be an alternative */
    {
      s2mMergeFusedCloneShadows(s2m) ;     /* merge shadows from fused clones */
      if (debug) s2mShowShadows (s2m->shadows, "fused clone shadows") ;
    }
  s2mMergeIdenticalShadows(s2m) ;     /* merge identical shadows */
  if (debug) s2mShowShadows (s2m->shadows, "merge-identical shadows") ;
  if (debug1)  freeOutf ("mrnaCreate: %d merge-identical shadows\n", arrayMax (s2m->shadows)) ;
  s2mCompleteShadows (s2m) ;  /* a short gap between exon is counted as exon */
  if (debug) s2mShowShadows (s2m->shadows, "gap-completed shadows") ;
  if (debug1)  freeOutf ("mrnaCreate: %d reads = initial shadows\n", arrayMax (s2m->shadows)) ;
  if (debug1)  freeOutf ("mrnaCreate: %d merge-included-nonPolyA shadows\n", arrayMax (s2m->shadows)) ;
  s2mMergeIncludedShadows(s2m, TRUE) ;     /* merge included shadows, merge polyA */
  if (debug) s2mShowShadows (s2m->shadows, "merge-included-polyA shadows") ;
  if (debug1)  freeOutf ("mrnaCreate: %d merge-included-PolyA shadows\n", arrayMax (s2m->shadows)) ; 
  if (debug1)  freeOutf ("mrnaCreate: %d final shadows\n", arrayMax (s2m->shadows)) ;
  /* check we do not make too many mergedshadows per contig */
  bestQual = 9 ;
  isDropped = TRUE ;
  s2mMakeContigs (s2m) ;      /* group shadows sharing an exon */
  if (debug) s2mShowContigs (s2m->contigs, "contigs", bestQual) ;
  while (isDropped)
    {
      bestQual-- ; isDropped = FALSE ;
      if (bestQual < 1) /* we are down to mrna and est quality 1 */
        {
          messerror ("error: mrnaCreate: %s problem1 with contig  in %s, after bestQual == 1",
                     timeShowNow(), name(cosmid)) ;
          break ;
        }
      /* if some contigs are too large drop some shadows in those contigs and iterate */
      for (ii = 0 ; ii < arrayMax (s2m->contigs) ; ii++)
        { 
          SC *sc = arrp (s2m->contigs, ii, SC) ;

          if (0 && debug) messout ("mrnaCreate1: qual=%d contig %d max(sc->ks) = %d", bestQual, ii,  arrayMax(sc->ks)) ;
          if (0 && debug1) freeOutf ("mrnaCreate1: qual=%d contig %d max(sc->ks) = %d", bestQual, ii,  arrayMax(sc->ks)) ;
          if (arrayMax(sc->ks) >  MAX_SHADOW && /* was 160 before aug 17 2005 */
              ! mrnaPleaseNoAbandon())         /* worse case in worm is 112 117 120 141 */
            {
              messout ("warning: mrnaCreate: %s s2mDropEsts1 q=%d in %s, max(sc->ks) = %d > %d OR diamMax=%d > %d",
                       timeShowNow(),  bestQual, name(cosmid), arrayMax(sc->ks), MAX_SHADOW, diamMax, DIAMMAX-2) ;
              if (debug) s2mShowContigs (s2m->contigs, "phase0a contigs ",bestQual) ;
              s2mDropEsts (s2m, sc, ii, bestQual, clonesWithBadIntron) ; isDropped = TRUE ;
              if (debug) s2mShowContigs (s2m->contigs, "phase0b contigs ",bestQual) ;
              isDropped = TRUE ;
            }
        }
      if (isDropped)
        {
          if (debug) s2mShowContigs (s2m->contigs, "phase1a contigs ",bestQual) ;
          s2mMakeContigs (s2m) ;      /* group shadows sharing an exon */
          if (debug) s2mShowContigs (s2m->contigs, "phase1b contigs",bestQual) ;
        }
    }
  
  /* check we do not make too complex diaminos per contig */
  bestQual++ ; isDropped = TRUE ;
  while (isDropped)
    { 
      bestQual-- ; isDropped = FALSE ;
      if (bestQual < -1)
        {
          messerror ("error: mrnaCreate: %s problem2 with contig  in %s, after bestQual == -1",
                     timeShowNow(), name(cosmid)) ;
          break ;
        }
      for (ii = 0 ; ii < arrayMax (s2m->contigs) ; ii++)
        { 
          SC *sc = arrp (s2m->contigs, ii, SC) ;

          if (debug)
            messout ("mrnaCreate2: qual=%d contig %d max(sc->ks) = %d", bestQual, ii,  arrayMax(sc->ks)) ;
          if (debug1)
            freeOutf ("mrnaCreate2: qual=%d contig %d max(sc->ks) = %d", bestQual, ii,  arrayMax(sc->ks)) ;
          if (arrayMax(sc->ks) >  MAX_SHADOW && 
              ! mrnaPleaseNoAbandon())         /* dropEst abovbe should have reduced max(sc->ks) */
            continue ;
          diamMax = DIAMMAX ;
          chrono ("s2mDiamino") ; 
          sc->s2mContigs = s2mDiamino (s2m, sc, &diamMax) ; 
          chronoReturn () ;
          
          if (!sc->s2mContigs || diamMax > DIAMMAX - 2)
            {
              messout ("warning: mrnaCreate: %s s2mDropEsts2 q=%d in %s, contig %d max(sc->ks) = %d  diamMax=%d > %d",
                       timeShowNow(),  bestQual, name(cosmid), ii, arrayMax(sc->ks), diamMax, DIAMMAX-2) ;
              s2mDropEsts (s2m, sc, ii, bestQual, clonesWithBadIntron) ; isDropped = TRUE ;
              sc->s2mContigs = 0 ;
              isDropped = TRUE ;
            }
        }
      if (isDropped)
        { 
          s2mMakeContigs (s2m) ;      /* group shadows sharing an exon */
          if (debug) s2mShowContigs (s2m->contigs, "phase 2b contigs",bestQual) ;
        }
    }

  for (ii = 0 ; ii < arrayMax (s2m->contigs) ; ii++)
    { 
      SC *sc = arrp (s2m->contigs, ii, SC) ;
      
      if (s2m->bMax && sc->s2mContigs)
        {
	  KEY gene ;
          smrnas = s2mMakeTranscripts (s2m, sc, sc->s2mContigs) ;
          if (!arrayMax (smrnas))
            continue ; /* feb 2004, could be messcrash, happen in rare low quality cases */
          gmrna = mrnaMakeGene (s2m, sc, sc->s2mContigs, smrnas) ;
          mrnaMakeEstHits (s2m, sc, gmrna) ;
          if (s2m->type == 2 && !mrnaFilterGene (sc))
            continue ;
	  gene = makeMrnaGene (s2m, sc, gmrna, smrnas, clipTops, clipEnds, linkPos) ;
	  if (gmrna->genes)
	    {/* composite method */
	      for (int i = 0 ; i < keySetMax (gmrna->genes) ; i++)
		keySet (genes, ig++) = keySet (gmrna->genes, i) ;
	    }
	  else
	    keySet (genes, ig++) = gene ;
        }
    }
  s2mDestroy (s2m) ;

  chronoReturn () ;
  if (debug1)
    {
      freeOutf ("Leaving mrnaCreate %s\n", name(cosmid)) ;
    }
  return genes ;
} /* mrnaCreate */

/**********************************************************************************/
/**********************************************************************************/
