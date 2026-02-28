/*  File: cdnaflipgene.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg. 1997
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

#define ARRAY_CHECK
#define MALLOC_CHECK


#include "acedb.h"
/* for webquery */
#include "lex.h"
#include "a.h"
#include "session.h"
/* for dna */
#include "dna.h"
#include "cdna.h"
#include "query.h"
#include "freeout.h"
#include "longtext.h"
#include "../whooks/systags.h"
#include "../whooks/tags.h"
#include "../whooks/classes.h"
#include "parse.h"
#include "pick.h"
#include "dict.h"
#include "vtxt.h"
#include "session.h"
#include "makemrna.h"

/***************************************************************/

void cDNAEliminateDeadMrnas (void)
{
  /* eliminate dead bodies */
  KEYSET ks = query (0, "Find mrna ! from_gene && ! from_prediction ; {IS *} $| {>dna} $| {>product} $| {>product ; >peptide } $| {>RefSeqMaker} $| {>RefSeqMaker; >dna}") ;
  keySetKill (ks) ; 
} /* cDNAEliminateDeadMrnas */

/***************************************************************/

static BOOL isGoodGene (KEY tg)
{
  BOOL ok = FALSE ;

  if (keyFindTag (tg, _gt_ag) ||
      keyFindTag (tg, _gc_ag))
    ok = TRUE ;
  return ok ;
}

/***************************************************************/

static int cDnaDoFlipRead (KEY read, KEYSET kFlipped)
{ 
  OBJ Read = 0 ;
  int nn = 0, x=0, nct_ac = 0, ngt_ag = 0, ngc_ag = 0 ;

  if (!keySetFind (kFlipped, read, 0))
    {
      if ((Read = bsUpdate (read)))
        {
          if (
              (bsFindTag (Read, _mForward) &&
               bsFindTag (Read, _Forward)
               ) ||
              (bsFindTag (Read, _mReverse) &&
               bsFindTag (Read, _Reverse)
               )        
              ) ;
          else
            {
              bsGetData (Read, _gt_ag, _Int, &ngt_ag) ;
              bsGetData (Read, _gc_ag, _Int, &ngc_ag) ;
              bsGetData (Read, _ct_ac, _Int, &nct_ac) ;
              if (ngt_ag + ngc_ag == 0 || ngt_ag + ngc_ag < nct_ac)
                {
		  bsAddTag (Read, _Flipped) ;
                  if (bsFindTag (Read, _Forward))
                    bsAddTag (Read, _Reverse) ;
                  else
                    { bsAddTag (Read, _Forward) ; x = 1 ;}
                  if (bsFindTag (Read, _Intron_boundaries))
                    { bsRemove (Read) ; x = 0 ; }
                  bsSave (Read) ;
                  freeOutf("// Flipping Read %s now becomes %s\n", 
                           name(read), x ? "Forward" : "Reverse") ;
                  nn = 1;
                }
              keySetInsert (kFlipped, read) ;
            }
        }
      bsSave (Read) ;
    }
  return nn ;
}

/***************************************************************/

static int  cDnaDoFlipGene (KEY tg, KEYSET kFlipped, KEY heros)
{
  int ii, nn = 0 ;
  KEYSET reads = queryKey (tg, ">Read") ;

  freeOutf("// Flipping Gene %s => %d reads to fit into %s\n", name(tg), keySetMax(reads), 
           heros ? name(heros) : "narcisse") ;
  for (ii = 0 ; ii < keySetMax (reads) ; ii++)
    nn += cDnaDoFlipRead (keySet (reads, ii), kFlipped) ;
    
  keySetDestroy (reads) ;


  return nn ;
}

/***************************************************************/

static int cDnaFlipTest (KEY tg1, KEY tg2, KEYSET kFlipped, BOOL countClones)
{
  Array units = 0, hits1 = 0, hits2 = 0 ;
  AC_HANDLE handle = 0 ;
  OBJ Tg1 = 0, Tg2 = 0 ;
  KEY map1 = 0, map2 = 0 ;
  int ii, jj, a1, a2, b1, b2, nn = 0 ;
  BSunit *uu ;
  HIT *h1, *h2 ;
  BOOL ok = FALSE, ok2 ;

  if  (isGoodGene  (tg2))
    return FALSE ;

  if (0) freeOutf("//comparing %s %s\n", name(tg1), name(tg2)) ;

  /* initialise */
  handle = handleCreate () ;
  units= arrayHandleCreate (64, BSunit, handle) ;
  hits1 = arrayHandleCreate (64, HIT, handle) ;
  hits2 = arrayHandleCreate (64, HIT, handle) ;

  Tg1 = bsCreate (tg1) ;
  Tg2 = bsCreate (tg2) ;
  
  if (!isGoodGene (tg1)) /* gene without gt_ag intron */
    {   /* turn around those genes with less than 1/3 of the clones */
      ii = jj = 0 ;
      if (bsGetArray (Tg1, _cDNA_clone, units, 1))
        ii = arrayMax (units) ;
      if (bsGetArray (Tg2, _cDNA_clone, units, 1))
        jj = arrayMax (units) ;
      if (ii > 3 * jj) 
        nn = cDnaDoFlipGene (tg2, kFlipped, tg1) ;
      goto abort ; 
    }

  /* get relative exons coordinates */
  if (bsGetArray (Tg1, _Splicing, units, 3))
    for (ii = jj = 0 ; ii < arrayMax (units) ; ii += 3)
      {
        uu = arrp (units, ii, BSunit) ;
        h1 = arrayp (hits1, jj++, HIT) ;
        h1->a1 = uu[0].i ;
        h1->a2 = uu[1].i ;
        h1->gene = tg1 ;
        h1->est = uu[2].k ;
      }

  if (bsGetArray (Tg2, _Splicing, units, 3))
    for (ii = jj = 0 ; ii < arrayMax (units) ; ii += 3)
      {
        uu = arrp (units, ii, BSunit) ;
        h2 = arrayp (hits2, jj++, HIT) ;
        h2->a1 = uu[0].i ;
        h2->a2 = uu[1].i ;
        h2->gene = tg1 ;
        h2->est = uu[2].k ;
      }

  if (!arrayMax (hits1) || ! arrayMax (hits2))
    goto abort ;

  /* get absolute exons coordinates */
  if (bsGetKey (Tg1, _IntMap, &map1) &&
      bsGetData (Tg1, _bsRight, _Int, &a1) &&
      bsGetData (Tg1, _bsRight, _Int, &b1) &&
      bsGetKey (Tg2, _IntMap, &map2) &&
      bsGetData (Tg2, _bsRight, _Int, &a2) &&
      bsGetData (Tg2, _bsRight, _Int, &b2))
    {
      for (ii = 0, h1 = arrayp (hits1, 0, HIT) ; ii < arrayMax (hits1) ; ii++, h1++)
        { 
          if (a1 < b1)
            { h1->x1 = a1 + h1->a1 - 1 ; h1->x2 = a1 + h1->a2 - 1 ; h1->reverse = FALSE ; }
          else
            { h1->x1 = a1 - h1->a1 + 1 ; h1->x2 = a1 - h1->a2 + 1 ; h1->reverse = TRUE ; }
        }
      for (ii = 0, h2 = arrayp (hits2, 0, HIT) ; ii < arrayMax (hits2) ; ii++, h2++)
        { 
          if (a2 < b2)
            { h2->x1 = a2 + h2->a1 - 1 ; h2->x2 = a2 + h2->a2 - 1 ; h2->reverse = FALSE ; }
          else
            { h2->x1 = a2 - h2->a1 + 1 ; h2->x2 = a2 - h2->a2 + 1 ; h2->reverse = TRUE ; }
        }
    }

  /* look for full inclusion of tg2 in tg1 */
  cDNASwapX (hits1) ;
  arraySort (hits1, cDNAOrderByX1) ;
  cDNASwapX (hits2) ;
  arraySort (hits2, cDNAOrderByX1) ;

  ii = jj = 0 ;   
  h1 = arrayp (hits1, 0, HIT) ; h2 = arrayp (hits2, 0, HIT) ; 
  ok = TRUE ;
  for (jj= 0 ; ok && jj < arrayMax (hits2) ; h2++, jj++)
    {
      if (!strstr (name(h2->est), "xon"))
        continue ;
      if (h2->x2 - h2->x1 < 12)
        continue ;
      ok2 = FALSE ;
      for (; ii < arrayMax (hits1) ; h1++, ii++)
        { 
          if (!strstr (name(h1->est), "xon"))
            continue ;
          if (h1->x1 < h2->x1 + 200 && h1->x2 > h2->x2 - 200) /* was 50, may 2 2004 */
            { ok2 = TRUE ; break ; }
        }
      if (!ok2) 
        ok = FALSE ;
    }
    
  if (ok)
    nn = cDnaDoFlipGene (tg2, kFlipped, tg1) ;

 abort:
  handleDestroy (handle) ;
  bsDestroy (Tg1) ;
  bsDestroy (Tg2) ;
  return nn ;
}

/***************************************************************/

static int cDnaFlipGene (KEY tg, KEYSET kFlipped)
{
  BOOL done = FALSE ;
  int ii, n, nn = 0 ;
  KEYSET ks = 0 ;

  if (keyFindTag (tg, _ct_ac) &&
      !keyFindTag (tg, _gt_ag))
    { done = TRUE ; cDnaDoFlipGene (tg, kFlipped, 0) ; }

  if (!done)
    {  
      ks = queryKey (tg, ">Antisens_to") ;
      
      for (ii = 0 ; ii < keySetMax (ks) ; ii++)
        {
          n = cDnaFlipTest (tg, keySet (ks, ii), kFlipped, FALSE) ;
          if (n) done = TRUE ; 
          nn += n ;
        }
      keySetDestroy (ks) ;
    }
  if (!done && keyFindTag (tg, _ct_ac))
    {
      OBJ Tg = 0 ;
      int nct = 0, ngt = 0 ;

      if ((Tg = bsCreate (tg)))
        {
          bsGetData (Tg, _gt_ag, _Int, &ngt) ;
          bsGetData (Tg, _ct_ac, _Int, &nct) ;
          bsDestroy (Tg) ;
        }      
      if (nct > 2 * ngt)
        { done = TRUE ; cDnaDoFlipGene (tg, kFlipped, 0) ; }        
    }
  
  return nn ;
}

/***************************************************************/

int cDnaFlipGeneKeySet (KEYSET ks) 
{
  int i, nn = 0 ;
  KEYSET ks1 = 0, ks2 = 0, kEst = 0, kFlipped = keySetCreate ()  ;
  KEY *kp ;

  freeOutf ("cDnaFlipGeneKeySet phase 1\n") ;
  if (ks)
    ks1 = ks ;
  else
    {
      ks1 = query (0, "Find tg") ;
      kEst = query (0, "Find EST ct_ac && ! gt_ag  && ! manualStrand") ;
    }
  ks2 = query (ks1, "CLASS Transcribed_gene && Antisens_to") ;
  if (arrayMax(ks2))
    for (i = 0, kp = arrp (ks2, 0, KEY) ; i < keySetMax(ks2) ; kp++, i++)
      {
        nn += cDnaFlipGene (*kp, kFlipped) ;
        /* printf (" cDnaFlipGeneKeySet i=%d nn = %d\n", i, nn) ; */
      }
        
  if (ks1 != ks)
    keySetDestroy (ks1) ;
  keySetDestroy (ks2) ;
  
  freeOutf ("// cDnaFlipGeneKeySet phase 2 max(kest)=%d\n", kEst ? arrayMax(kEst) : 0) ;
  if (kEst && arrayMax(kEst))
    for (i = 0, kp = arrp (kEst, 0, KEY) ; i < keySetMax(kEst) ; kp++, i++)
      nn += cDnaDoFlipRead (*kp, kFlipped) ;
  
  keySetDestroy (kEst) ;

  kEst = query (0, "Find Read (Forward && mReverse) || (reverse && mForward)") ;
  if (kEst && arrayMax(kEst))
    for (i = 0, kp = arrp (kEst, 0, KEY) ; i < keySetMax(kEst) ; kp++, i++)
      cDnaDoFlipRead (*kp, kFlipped) ;  /* do not count these, they are final */
  

  keySetDestroy (kFlipped) ;
  
  freeOutf("// Flipped %d  Read\n", nn) ;
  return nn ;
}

/***************************************************************/
/***************************************************************/
/************** Strategy, imported from clone MainClone **************/

static int getPleaseMinIntronSize (void)
{
  static int minIntronSize = -1 ;

  if (minIntronSize == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && MinIntronSize") ;
      minIntronSize = keySetMax (ks) ? 1 : -1 ;
      
      if (minIntronSize == 1)
        {
          OBJ Clone = bsCreate (keySet (ks,0)) ;
          
          if (Clone)
            {
              bsGetData (Clone, str2tag("MinIntronSize"), _Int, &minIntronSize) ;
              bsDestroy (Clone) ;
            }
        }
      keySetDestroy (ks) ;
    }
  if (minIntronSize == -1)
    minIntronSize = 0 ; /* default */

  return minIntronSize ;
} /* getPleaseMinIntronSize */

/***************************************************************/
/***************************************************************/

static int cDnaFlagOneSuspectSalvage (OBJ TG, KEYSET clones, KEY clone)
{
  static Array introns = 0 ;
  BSunit *uu ;
  int ii, oldX1 = 0, oldX2 = 0, x1, x2, nn = 0, n, ngt = 0, nct = 0, ncomp = 0;
  KEY type, oldType = 0, clo ;
  BOOL foundClo = FALSE ;
  KEYSET reads = queryKey (clone, ">Read") ;
  OBJ Est = 0 ;

  for (ii = 0 ; ii < keySetMax (reads) ; ii++)
    if ((Est = bsCreate (keySet(reads, ii))))
      {  
        if (bsGetData (Est, _gt_ag, _Int, &n)) ngt += n ;
        if (bsGetData (Est, _gc_ag, _Int, &n)) ngt += n ;
        if (bsGetData (Est, _Composite, _Int, &n)) ncomp += n ;
	if (bsGetData (Est, _ct_ac, _Int, &n)) nct += n ;
	if (bsFindTag (Est, _Ref_Seq)) nct += 100000 ; /* salvage */
        bsDestroy (Est) ;
      }
  keySetDestroy (reads) ;

  if (ncomp > 2)
    {
      reads = queryKey (bsKey (TG), ">Read Composite && gt_ag") ;
      ii = arrayMax (reads) ;
      if (ii)
	{
	  int jj = 0 ;
	  Array aa = arrayCreate (ii, int) ;
	  
	  while (ii--)
	    {
	      if ((Est = bsCreate (keySet(reads, ii))))
		{  
		  if (bsGetData (Est, _Composite, _Int, &n)) array (aa, jj++, int) = n ;
		  bsDestroy (Est) ;
		}
	    }

	  arraySort (aa, intOrder) ;
	  if (! jj || ncomp >= array (aa, jj/2, int))
	    return 2 ;
	}
      keySetDestroy (reads) ;
    }

  if (nct > 9999 && nct >= ngt) /* do not flag as deletion a reversed clone */
    return 2 ;
  if (nct > 0 && nct >= ngt) /* do not flag as deletion a reversed clone */
    return 1 ;

  /* given that we do not like this intron
   * we wonder if the suspect clone is the only support of another gt_ag intron
   * note that if we successively kill several clones
   * we mot realize that we are killing another good intron
   */
  introns = arrayReCreate (introns, 64, BSunit) ;
  if (bsGetArray (TG,  str2tag ("Intron_boundaries"), introns, 5))
    for (nn = ii = 0 ; ii < arrayMax(introns) ; ii += 5)
      {
        uu = arrp (introns, ii, BSunit) ;
        type = uu[0].k ;
	/*         nb = uu[1].i ; */
        x1 = uu[2].i ;
        x2 = uu[3].i ;
        clo = uu[4].k ;
        
        if (type != _gt_ag && type != _gc_ag)
          continue ;
        if (type != oldType || x1 != oldX1 || x2 != oldX2) 
          { 
            if (foundClo && !nn) /* clo is needed to support this intron */
              return 1 ;
            nn = 0 ; foundClo = FALSE ;
          }
        oldType = type ; oldX1 = x1 ; oldX2 = x2 ;
        if (!nn)
          {
            if (clo == clone) 
              foundClo = TRUE ;
            else if (!keySetFind (clones, clo, 0)) /* somebody else */
              nn++ ;
          }
      }
  if (foundClo && !nn) /* clo is needed to support this intron */
    return 1 ;
  return 0 ;
} /* cDnaFlagOneSuspectSalvage */

/***************************************************************/

int cDnaFlagOneSuspectIntrons (KEY tg, BOOL doIgnore, KEYSET result)
{
  int nn = 0, ii, jj, i1, i2 ;
  int limit = 3 ; /* default */
  OBJ TG = 0 ;
  HIT *up ;
  Array suspects = 0 ; 
  Array others = arrayCreate (64, BSunit) ;
  Array small = arrayCreate (64, BSunit) ;
  Array fuzzies = arrayCreate (64, BSunit) ;
  Array at_acs = arrayCreate (64, BSunit) ;
  Array ct_acs = arrayCreate (64, BSunit) ;
  Array afroms = arrayCreate (64, BSunit) ;
  Array spls = arrayCreate (64, BSunit) ;
  BSunit *uu, *vv ;
  KEY clone ;
  KEYSET clones = 0, badClones = 0, horribleClones = 0, badClones2 = keySetCreate () ;
  KEY _Ignore_this_clone = str2tag ("Ignore_this_clone") ;
  KEY _Ignore_this_clone_automatic = str2tag ("Ignore_this_clone_automatic") ;
  KEY _Duplicate_clone = str2tag ("Duplicate_clone") ;
  int minIntronSize = getPleaseMinIntronSize () ;

  clones = queryKey (tg, ">cdna_clone") ;
  ii = keySetMax (clones) ;
  keySetDestroy (clones) ;
  ii /= 30 ;
  if (ii > limit) limit = ii ;

  badClones = queryKey (tg, ">cdna_clone ; (double_fuzzy || Suspected_internal_deletion) && ! Ignore_this_clone && !Ignore_this_clone_automatic") ;
  if ((TG = bsCreate (tg)))
    {
      bsGetArray (TG, _Other, others, 5) ;
      bsGetArray (TG, str2tag("at_ac"), at_acs, 5) ;
      bsGetArray (TG, _ct_ac, ct_acs, 5) ;
      bsGetArray (TG, str2tag ("Small_deletion"),  small, 4) ;
      bsGetArray (TG, _Fuzzy, fuzzies, 5) ;
      bsGetArray (TG, _Splicing, spls, 5) ;
      bsGetArray (TG, _Assembled_from, afroms, 5) ;
      
      /* we cat the at-ac on the others */
      for (ii = 0, jj = arrayMax(others) ; ii < arrayMax (at_acs) ; ii += 5)
        {
          uu = arrp (at_acs, ii, BSunit) ;
          clone = uu[3].k ; 
          if (keyFindTag (clone, _Ignore_this_clone) ||
              keyFindTag (clone, _Ignore_this_clone_automatic))
            continue ;
          array (others, jj + 4, BSunit).i = 0 ; /* make room */
          vv = arrayp (others, jj, BSunit) ; jj += 5 ;
          vv[0].i = 0 ; 
          vv[1] = uu[0] ;
          vv[2] = uu[1] ;
          vv[3] = uu[2] ;
          vv[4].k = clone ;
        }

      /* we cat the ct-ac on the others */
      for (ii = 0, jj = arrayMax(others) ; ii < arrayMax (ct_acs) ; ii += 5)
        {
          uu = arrp (ct_acs, ii, BSunit) ;
          clone = uu[3].k ; 
          if (keyFindTag (clone, _Ignore_this_clone) ||
              keyFindTag (clone, _Ignore_this_clone_automatic))
            continue ;
          array (others, jj + 4, BSunit).i = 0 ; /* make room */
          vv = arrayp (others, jj, BSunit) ; jj += 5 ;
          vv[0].i = 0 ; 
          vv[1] = uu[0] ;
          vv[2] = uu[1] ;
          vv[3] = uu[2] ;
          vv[4].k = clone ;
        }

      /* we cat the micro_gc on the others */
      for (ii = 0, jj = arrayMax(others) ; ii < arrayMax (small) ; ii += 4)
        {
          uu = arrp (small, ii, BSunit) ;
          if (uu[2].i - uu[1].i + 1 >=  minIntronSize)
            continue ;
          clone = uu[3].k ; 
          if (keyFindTag (clone, _Ignore_this_clone) ||
              keyFindTag (clone, _Ignore_this_clone_automatic))
            continue ;
          array (others, jj + 4, BSunit).i = 0 ; /* make room */
          vv = arrayp (others, jj, BSunit) ; jj += 5 ;
          vv[0].i = 0 ; 
          vv[1].i = 0 ;
          vv[2] = uu[1] ;
          vv[3] = uu[2] ;
          vv[4].k = clone ;
        }

      /* we only include the non gt_ag gc_ag fuzzy cases */
      for (ii = 0, jj = arrayMax(others) ; ii < arrayMax (fuzzies) ; ii += 5)
        {
          uu = arrp (fuzzies, ii, BSunit) ;
          clone = uu[3].k ; 
          if (keyFindTag (clone, _Ignore_this_clone) ||
              keyFindTag (clone, _Ignore_this_clone_automatic))
            continue ;
          if (!uu[4].s ||
              (uu[4].s && !strstr(uu[4].s, "gt_ag") && !strstr(uu[4].s, "gc_ag"))
               )
            {
              array (others, jj + 4, BSunit).i = 0 ; /* make room */
              vv = arrayp (others, jj, BSunit) ; jj += 5 ;
              vv[0].i = 0; 
              vv[1] = uu[0] ;
              vv[2] = uu[1] ;
              vv[3] = uu[2] ;
              vv[4].k = clone ;
            }
        }

      /* we now check all the different types of suspect clones at once */
      for (ii = 0 ; ii < arrayMax (others) ; ii += 5)
        {
          uu = arrp (others, ii, BSunit) ;
          uu[0].i = 1 ;  /* prepare to block reiterations */
        }

      for (ii = 0 ; ii < arrayMax (others) ; ii += 5)
        {
          uu = arrp (others, ii, BSunit) ;
          if (!uu[0].i) /* block reiterations */
            continue ;
          clones = keySetReCreate (clones) ;
          for (i1 = i2 = 0 ; ii + 5*i1 + 3 < arrayMax (others) ; i1++)
            if (uu[2].i == uu[2 + 5*i1].i && uu[3].i == uu[3 + 5*i1].i)
              { 
                uu[5*i1].i = 0 ; /* case done */ /* block reiterations */
                clone  = uu[4 + 5*i1].k ;
                keySetInsert (clones, clone) ;
                if (!keyFindTag (uu[4 + 5*i1].k,_Duplicate_clone))
                  i2++ ;
              }
          if (1) /* no limit if clone is covered, limit applies only to non covered cases */
            {
              for (i1 = 0 ; i1 < keySetMax (clones) ; i1++)
                {
                  clone = keySet (clones, i1) ;        
                  if (!keyFindTag (clone, _Manual_no_internal_deletion))
                    {
                      if (
                          (minIntronSize &&  uu[3].i - uu[2].i + 1 < minIntronSize && i2 < 3*limit) ||
                          i2 <= limit)  /* good quasi unique clone sustaining the intron */
                        keySetInsert (badClones, clone) ;  /* store the bad clone */
                      if (!suspects)                     /* store the message */
                        suspects = arrayCreate (12, HIT) ;
                      up = arrayp (suspects, arrayMax(suspects), HIT) ;
                      up->cDNA_clone = clone ;
                      for (jj = 0 ; jj < arrayMax (afroms) ; jj += 5)
                        {
                          vv = arrp (afroms, jj, BSunit) ;
                          /* do not break, this table is not ordered on vv[1] == hit->a2 */
                          if (vv[1].i >= uu[2].i - 5 && vv[1].i <= uu[2].i + 4 && /* fuzzy case */
                              (clone == keyGetKey (vv[2].k, _cDNA_clone))
                              )
                            {
                              up->est = vv[2].k ;
                              up->x1 = vv[4].i ;
                              /* do not break i may have to mark both 3p and 5p reads */
                            }                    
                        } 
                      
                      /* check that we will not create a gap because another read will help us */
                      {
                        int ok = 0 ;
                        BSunit *vv2 ;
                        int jj2 ;
                        
                        for (jj = 0 ; ok <= i2 && jj < arrayMax (afroms) ; jj += 5)
                          {
                            vv = arrp (afroms, jj, BSunit) ;
                            if (vv[0].i <= uu[2].i  && vv[1].i >= uu[2].i - 5) 
                              for (jj2 = jj ; ok <= i2 && jj2 < arrayMax (afroms) ; jj2 += 5) 
                                {
                                  vv2 = arrp (afroms, jj2, BSunit) ;
                                  if (vv2[2].k == vv[2].k && /* same EST */
                                      vv2[0].i <= uu[3].i + 5  && vv2[1].i > uu[3].i) 
                                    {
                                      KEY clo = keyGetKey (vv[2].k, _cDNA_clone) ;
                                      if (!clo || !keyFindTag (clo,_Duplicate_clone))
                                        ok++ ;
                                    }
                                }
                          }
                        if (ok > i2 &&   /* so we are covered */
                            (i2 < 3 || i2 < ok - i2)
                            )
                          keySetInsert (badClones2, clone) ;  /* this clone is covered */
                      }
                    }
                }
            }
        }
      /* add the clones with a gap */
      {
        KEYSET gapClones = 0, ks = 0 ;
        gapClones = queryKey (tg, "total_gap_length ; >mrna gap_length ; >cdna_clone ; >read gap_length ; >cdna_clone") ;
        ks = badClones ;
        badClones = keySetOR (ks, gapClones) ;
        keySetDestroy (ks) ; 
        ks = badClones2 ;
        badClones2 = keySetOR (ks, gapClones) ;
        keySetDestroy (ks) ; 
        keySetDestroy (gapClones) ;
      }

      if (doIgnore)
	{
	  KEY est ;
	  KEYSET ests = queryKey (tg, ">read ; Fuzzy OR Other") ;
	  OBJ Est = 0 ;

	  for (i1 = 0 ; i1 < keySetMax (ests) ; i1++)
	    {
	      est = keySet (ests, i1) ;
	      if ((Est = bsCreate (est)))
		{
		  int good = 0, bad = 0 ;
		  clone = 0 ;
		  bsGetKey (Est, _cDNA_clone, &clone) ;
		  bsGetArray (Est, _Intron_boundaries, others, 3) ;
		  for (ii = jj = 0 ; jj < arrayMax(others) ; jj += 3)
		    {
		      vv = arrayp (others, jj, BSunit) ;
		      if (vv[0].k == _gt_ag)
			good += vv[1].i ;
		      else if (vv[0].k == _Other)
			bad += vv[2].i ;
		      else if (vv[0].k == _Fuzzy)
			bad += vv[1].i ;	    
		    }

		  if (! bsFindTag (Est, _Ref_Seq) && 
		      ((bad > 1 && 3*bad > good) || bad > good)
		      )
		    {
		      if (!horribleClones)
			horribleClones = keySetCreate () ;	
		      keySetInsert (horribleClones, clone) ;
		      keySetInsert (badClones, clone) ;
		    }
		  bsDestroy (Est) ;
		}		  
	    }
	  keySetDestroy (ests) ;
	  ests = queryKey (tg, ">mrna ; fuzzy || other || gap ; COUNT cdna_clone = 1 ; >cdna_clone") ;
	  for (i1 = 0 ; i1 < keySetMax (ests) ; i1++)
	    {
	      clone = keySet (ests, i1) ;
	      keySetInsert (badClones, clone) ;
	    }	      
	}
	      
      /* actualy flag those clones that do not bring a good gt-ag */
      if (keySetMax (badClones))
        { 
          OBJ Clone = 0 ;
          int salvage ;
          BOOL ok ;

          for (i1 = 0 ; i1 < keySetMax (badClones) ; i1++)
            {
              ok = FALSE ;
              clone = keySet (badClones, i1) ;
	      if (horribleClones && keySetFind (horribleClones, clone, 0))
		salvage = 0 ;
	      else
		salvage = cDnaFlagOneSuspectSalvage (TG, badClones, clone) ;
              if (salvage < 2 &&
                  (Clone = bsUpdate (clone)))
                {
                  if (! bsFindTag (Clone, _Ignore_this_clone) &&
                      ! bsFindTag (Clone, _Ignore_this_clone_automatic) &&
                      ! bsFindTag (Clone, _Manual_no_internal_deletion))
                    {
                      if (doIgnore && !salvage)
                        { 
                          ok = TRUE ;
                          bsAddTag (Clone, _Ignore_this_clone_automatic) ;
                          bsAddData (Clone, _bsRight, _Text, "Automatically added by FlagIgnoreIntrons2") ;
                          bsAddKey (Clone, _Colour, _PALERED) ;
                        }
                      if (keySetFind (badClones2, clone, 0))
                        {
                          ok = TRUE ;
                          bsAddTag (Clone, _Suspected_internal_deletion) ;
                          for (ii = 0 ; suspects && ii < arrayMax (suspects) ; ii++)
                            {
                              up = arrp (suspects, ii, HIT) ;
                              if (up->cDNA_clone == clone && up->est)
                                {
                                  bsAddKey (Clone, _Suspected_internal_deletion, up->est) ;
                                  bsAddData (Clone, _bsRight, _Int, &(up->x1)) ;
                                }                          
                            }
                        }
                      if (ok)
                        {
                          freeOutf ("// Flagged clone %s as %s in tg %s\n"
                                    , name(clone)
                                    , doIgnore ? "ignore_this_clone_automatic" : "suspect" 
                                    , name(tg)) ;
                          nn++ ;
                        }
                      if (result)
                        keySet (result, keySetMax (result)) = clone ;
                    }
                  bsSave (Clone) ;
                  /* do not break i may have to mark both 3p and 5p reads */
                } 
            }
        }
      bsDestroy (TG) ;
    }
  
  arrayDestroy (suspects) ;
  arrayDestroy (others) ;
  arrayDestroy (small) ;
  arrayDestroy (fuzzies) ;
  arrayDestroy (at_acs) ;
  arrayDestroy (ct_acs) ;
  arrayDestroy (afroms) ;
  arrayDestroy (spls) ;
  keySetDestroy (clones) ;
  keySetDestroy (horribleClones) ;
  keySetDestroy (badClones) ;
  keySetDestroy (badClones2) ;

  return nn ;
} /* cDnaFlagOneSuspectIntrons */

/***************************************************************/

int cDnaFlagSuspectIntrons (KEYSET genes, BOOL doIgnore) 
{
  int ii, jj, nn = 0, ng = 0, n1 ;
  KEYSET ks1 ;

  ks1 = query (genes, "CLASS Transcribed_gene ; Other || at_ac || Fuzzy || Small_deletion || Total_gap_length ") ;

  for (ii = jj = 0 ; ii < keySetMax(ks1) ; ii++)
    {
      n1 = cDnaFlagOneSuspectIntrons (keySet (ks1,ii), doIgnore, 0) ;
      if (n1 > 0) { nn += n1 ; ng++ ; keySet (genes, jj++) = keySet (ks1,ii) ; }
    }
  keySetMax (genes) = jj ;
  keySetDestroy (ks1) ;
  if (!nn)  freeOutf("\n") ;
  freeOutf("// Flagged %d suspect introns %sin %d genes\n", nn
           , doIgnore ? "as ignore_this_clone " : "", ng) ;
  return jj ;
} /* cDnaFlagSuspectIntrons */

/***************************************************************/
/***************************************************************/
   
int cDnaExportGeneboxPrimers (KEYSET ks0) 
{
  int i, ii, nn = 0, a1, a2, b1, b2 ;
  OBJ Cosmid = 0 ;
  KEY gene, cosmid ;
  Array dna = 0, p1 = 0, p2 = 0 ;
  KEYSET ks = 0 ;
  int SIZE = 30 ;
  BOOL isUp = FALSE ;

  if (!ks0 || ! arrayMax (ks0)) ks = query (0, "Find Gene") ;
  else ks = query (ks0, "CLASS Gene") ;

  p1 = arrayCreate (SIZE + 2, char) ;
  p2 = arrayCreate (SIZE + 2, char) ;

  array (p1, SIZE, char) = 0 ; /* make zero terminated room */
  arrayMax (p1) = SIZE ;
  array (p2, SIZE, char) = 0 ; /* make zero terminated room */
  arrayMax (p2) = SIZE ;

  for (ii = 0 ; ii < keySetMax(ks) ; ii++)
    {
      gene = keySet (ks, ii) ;
      cosmid = keyGetKey (gene, str2tag ("Genomic_sequence")) ;
      if (!cosmid ||
          !(Cosmid = bsCreate (cosmid)))
        continue ;

       if (0 && bsFindKey (Cosmid, _Genes, gene) &&
           bsGetData (Cosmid, _bsRight, _Int, &b1) &&
           bsGetData (Cosmid, _bsRight, _Int, &b2) &&
           (dna = dnaGet (gene)))
        {
          if (dna && arrayMax(dna) > SIZE + 6)
            {
              a1 = 1 ; a2 = arrayMax(dna) - 1 ;
              for (i = 0 ; i < SIZE ; i++)
                arr (p1, i, char) = dnaDecodeChar [(int)arr (dna, a1 + i - 1, char)] ;
              for (i = 0 ; i < SIZE ; i++)
                arr (p2, i, char) = dnaDecodeChar [(int)complementBase(arr (dna, a2 - i - 1, char))] ;
              
              freeOutf ("Primer gbx_%s_f\nForward\nMotif %s\nGene %s\n\n", name(gene), arrp(p1, 0, char), name(gene)) ;
              freeOutf ("Primer gbx_%s_r\nReverse\nMotif %s\nGene %s\n\n", name(gene), arrp(p2, 0, char), name(gene)) ;
              freeOutf ("Gene %s\nPrevious_cosmid %s %d %d\n\n", name (gene), name(cosmid), b1, b2) ;
            }
          else
            freeOutf ("Gene %s // ERROR dnaMax = %d\n\n",  name (gene), arrayMax(dna)) ;
          arrayDestroy (dna) ;
        }          

      if (1 && bsFindKey (Cosmid, _Genes, gene) &&
          bsGetData (Cosmid, _bsRight, _Int, &a1) &&
          bsGetData (Cosmid, _bsRight, _Int, &a2) &&
          (dna = dnaGet (cosmid)))
        { 
          if (0) printf ("method2 %d ==? %d  a1=%d a2=%d\n",
                  dna ? arrayMax(dna) : 0, a2 - a1, a1, a2) ;
          b1 = a1 ; b2 = a2 ;
          if (a1 > a2) { isUp = TRUE ; i = a1 ; a1 = a2 ; a2 = i ; }
          else isUp = FALSE ;

          if (a2 < a1 + SIZE + 1) a2 = a1 + SIZE + 1 ;
          if (a1 < a2 && a1 >= 0 && a2 < arrayMax(dna) && a2 - a1 > SIZE)
            {
              for (i = 0 ; i < SIZE ; i++)
                arr (p1, i, char) = dnaDecodeChar [(int)arr (dna, a1 + i - 1, char)] ;
              for (i = 0 ; i < SIZE ; i++)
                arr (p2, i, char) = dnaDecodeChar [(int)complementBase(arr (dna, a2 - i - 1, char))] ;

              if (isUp) { Array p3 = p1 ; p1 = p2 ; p2 = p3 ; }
              freeOutf ("Primer gbx_%s_f\nForward\nMotif %s\nGene %s\n\n", name(gene), arrp(p1, 0, char), name(gene)) ;
              freeOutf ("Primer gbx_%s_r\nReverse\nMotif %s\nGene %s\n\n", name(gene), arrp(p2, 0, char), name(gene)) ;
              freeOutf ("Gene %s\nPrevious_cosmid %s %d %d\n\n", name (gene), name(cosmid), b1, b2) ;
            }
          arrayDestroy (dna) ;
        }
      bsDestroy (Cosmid) ;
    }
  
  keySetDestroy (ks) ;
  arrayDestroy (dna) ;
  arrayDestroy (p1) ;
  arrayDestroy (p2) ;
  return nn ;
}

/***************************************************************/

static int cdnaTagFloating3pGenes (KEYSET tgs0)
{
  int ii, jj, a1, a2, b1, b2, nn = 0 ;
  OBJ Tg0, Tg ;
  KEYSET tgs = 0, reads0 = 0, reads1 = 0, readsBoth = 0 ;
  KEY chrom1, chrom2, tg, tg0;
  BOOL ok ;

  for (ii = 0 ; ii < keySetMax (tgs0) ; ii++)
    {
      tg0 = keySet (tgs0, ii) ;
      ok = FALSE ;
      if ((Tg0 = bsCreate (tg0)))
        {
          if (bsGetKey (Tg0, _IntMap, &chrom1) &&
              bsGetData (Tg0, _bsRight, _Int, &a1) &&
              bsGetData (Tg0, _bsRight, _Int, &a2))
            ok = TRUE ;
          bsDestroy (Tg0) ;
        }
      if (!ok) continue ;
      keySetDestroy (tgs) ;
      keySetDestroy (reads0) ;
      reads0 = queryKey (tg0, ">read") ;
      tgs = queryKey (tg0, "{>cdna_clone ; {>parent_clone} SETOR {>child_clone} ; > from_gene} SETMINUS {>to_be_fused_with}") ;
      for (jj = 0 ; jj < keySetMax (tgs) ; jj++)
        {
          ok = FALSE ;
          tg = keySet (tgs, jj) ;
          if (tg == tg0)
            continue ;
          if ((Tg = bsCreate (tg)))
            {
              if (bsGetKey (Tg, _IntMap, &chrom2) &&
                  bsGetData (Tg, _bsRight, _Int, &b1) &&
                  bsGetData (Tg, _bsRight, _Int, &b2))
                ok = TRUE ;
              bsDestroy (Tg) ;
            }
          if (!ok) continue ;
          if (chrom1 != chrom2 ||
              ((b2 - b1) > 0 &&  (a2 - a1) < 0) ||
              ((b2 - b1) < 0 &&  (a2 - a1) > 0)
              )
            continue ;
          if (b1 < b2 && (a2 < b1 - 3000 || a1 > b2 + 3000))
            continue ;
          if (b2 < b1 && (a1 < b2 - 3000 || a2 > b1 + 3000))
            continue ;
          keySetDestroy (reads1) ;
          reads1 = queryKey (tg, ">read") ;
          keySetDestroy (readsBoth) ;
          readsBoth = keySetAND (reads0, reads1) ;
          if (keySetMax (readsBoth))
            continue ;
          if ((Tg0 = bsUpdate (tg0)))
            {
              bsAddKey (Tg0, str2tag ("To_be_fused_with"), tg) ;
              bsSave (Tg0) ;
              nn++ ;
            }
        }
    }
  keySetDestroy (tgs) ;
  keySetDestroy (reads0) ;
  keySetDestroy (reads1) ;
  keySetDestroy (readsBoth) ;
  return nn ;
} /* cdnaTagFloating3pGenes */

/***************************************************************/

static void cdnaTagSheddedGenesOneCluster (KEYSET tgs)
{
  KEYSET genes = 0, tgs1 = 0, principalTgs = 0, clones, tgPairs = keySetCreate () ;
  KEY tg, tg1, principalTg = 0 ;
  int i, j, ii, n, nMax = 0 ;
  vTXT txt = vtxtCreate () ;
  Array aa = arrayCreate (12, HIT) ;
  Array bb = arrayCreate (12, HIT) ;
  HIT *up, *vp ;
  
  /* select the principal genes */
  for (ii = 0 ; ii < keySetMax (tgs) ; ii++)
    {
      tg = keySet (tgs, ii) ;
      clones = queryKey (tg, ">cDNA_clone") ;
      n = keySetMax (clones) ;
      if (n > nMax)
        { principalTg = tg ; nMax = n ; }
      keySetDestroy (clones) ;
    }
  if (!principalTg)
    return ;
  principalTgs = query (tgs, messprintf ("cdna_clone = nm* ||  gt_ag || gc_ag || COUNT cdna_clone > %d", nMax/4, nMax/4)) ;
  if (keySetMax (principalTgs) > 1)
    {
      tgs1 = principalTgs ;
      principalTgs = query (tgs1, "cdna_clone = nm* || gt_ag || gc_ag") ;
      
      keySetDestroy (tgs1) ;
    }
  
  if (0 && keySetMax (principalTgs) > 1) 
    {
      /* no longer necessary since we re-use the clone pairs in the construction */
      /* avoid re splitting clone pairs */
      KEYSET ks ;
      BOOL found = TRUE, found2 = FALSE ;
      while (found)
        {
          for (ii = 0, found = FALSE ; ii < keySetMax ( principalTgs) ; ii++)
            {
              tg = keySet (principalTgs, ii) ;
              if (!tg) continue ;
              ks = queryKey (tg, ">cdna_clone ; {IS *} SETOR {>child_clone} SETOR {>parent_clone} ;>from_gene") ;
              for (i = 0 ; i < keySetMax (ks) ; i++)
                {
                  tg1 = keySet (ks, i) ;
                  /* remove gene friends from the principal gene list */
                  for (j = ii+1 ; j < keySetMax ( principalTgs) ; j++)
                    if (keySet (principalTgs, j) == tg1)
                      { keySet (principalTgs, j)  = 0 ; found = TRUE ; found2 = TRUE ; break ; }
                }
              keySetDestroy (ks) ;
            }
        }
      if (found2) /* keep happy few */
        {
          for (ii = i = 0 ; i < keySetMax (principalTgs) ; i++)
            {
              if (keySet (principalTgs, i))
                {
                  if (ii < i) keySet (principalTgs, ii) = keySet (principalTgs, i) ;
                  ii++ ;
                }
            }
          keySetMax (principalTgs) = ii ;
        }          
    }
  if (keySetMax (principalTgs) > 1)
    {
      /* regroup overlapping gene pairs */
      /* since we are inside the to be fused with system, we are not unduly grouping the cloud genes */
      OBJ Tg ;
      int a1, a2, jj, jj2, da, db, dc, c1, c2 ;
      KEY map ;
      
      /* collect coords */
      for (ii = jj = jj2 = 0 ; ii < keySetMax (tgs) ; ii++)
        {
          tg = keySet (tgs, ii) ;
          if ((Tg = bsCreate (tg)))
            {
              if (bsGetKey (Tg, _IntMap, &map) &&
                  bsGetData (Tg, _bsRight, _Int, &a1) &&
                  bsGetData (Tg, _bsRight, _Int, &a2)
                  )
                {
                  up = arrayp (bb, jj++, HIT) ;
                  up->gene = tg ;
                  if (a1 < a2)
                    { up->a1 = a1 ; up->a2 = a2 ; }
                  else
                    { up->a1 = a2 ; up->a2 = a1 ; }
                  up->x1 = up->a2 - up->a1 -1 ; /* true length */
                  if (keySetFind (principalTgs, tg, 0))
                    {
                      vp = arrayp (aa, jj2++, HIT) ;
                      *vp = *up ;
                    }                    
                }
              bsDestroy (Tg) ;
            }
        }
      arraySort (aa, cDNAOrderGloballyByA1) ; /* coords of principals */
      arraySort (bb, cDNAOrderGloballyByA1) ; /* coords of all */
      
      /* compare */
      for (ii = 0, up = arrp (aa, ii, HIT) ; ii < keySetMax (aa) ; up++, ii++)
        {
          da = up->a2 - up->a1 ; /* current length */ 
          if (up->x2) continue ;
          for (jj = ii+1, vp = up+1 ; jj < keySetMax (aa) ; vp++, jj++)
            {
              if (vp->x2) continue ;
              db = vp->a2 - vp->a1 ;
              c1 = up->a1 > vp->a1 ? up->a1 : vp->a1 ;
              c2 = up->a2 < vp->a2 ? up->a2 : vp->a2 ;
              dc = c2 - c1 ;
              if (dc <= 0) break ;
              /* if overlap, create a virtual large box */
              if (3*dc > 2*da || 3*dc > 2*db)
                { 
                  if (up->x1 < vp->x1)
                    { up->gene = vp->gene ; up->x1 = vp->x1 ; }
                  if (up->a2 < vp->a2)
                    { up->a2 = vp->a2 ; da = up->a2 - up->a1 ; }
                  vp->x2 = 1 ;
                }
            }          
        }
      keySetMax (principalTgs) = 0 ;
      for (ii = i = 0, up = arrp (aa, ii, HIT) ; ii < keySetMax (aa) ; up++, ii++)
        if (!up->x2)
          keySet (principalTgs, i++) = up->gene ;
    }
  
  keySetSort (principalTgs) ;
  keySetCompress (principalTgs) ;
  if (! keySetMax (principalTgs))
    keySet (principalTgs, 0) = principalTg ;
  principalTg = keySet (principalTgs, 0) ;
  /* tag the others */
  /* principalGene = keyGetKey (principalTg, _Gene) ; */
  tgs1 = keySetMINUS (tgs, principalTgs) ;
  
  for (ii = 0 ; ii < keySetMax (tgs1) ; ii++)
    {
      HIT *wp = aa && arrayMax (aa) ? arrp (aa, 0, HIT) : 0 ;
      tg = keySet (tgs1, ii) ;
      tg1 = principalTg ;
      /* locate coords and corresponding best */
      if (wp && keySetMax (principalTgs) > 1)
        {
          for (i = 0, up = arrp (bb, i, HIT) ; i < keySetMax (bb) ; i++, up++)
            if (up->gene == tg) 
              {
                int a1 = up->a1, a2 = up->a2, b1, b2, c1, c2, bestDb, bestDc ;
                b1 = wp->a1 ; b2 = wp->a2 ;
                c1 = a1 > b1 ? a1 : b1 ; c2 = a2 < b2 ? a2 : b2 ;
                bestDc = c2 - c1 ; bestDb = b2 - b1 ; 
                for (j = 1, vp = wp+1 ; j < keySetMax (aa) ; j++, vp++)
                  {
                    if (vp->x2) continue ; /* in aa, but eliminated from principalTgs */
                    b1 = vp->a1 ; b2 = vp->a2 ;
                    c1 = a1 > b1 ? a1 : b1 ; c2 = a2 < b2 ? a2 : b2 ;
                    if (c2 - c1 > bestDc || (c2 - c1 == bestDc && b2 - b1 < bestDb))
                      { wp = vp ; bestDc = c2 - c1 ; bestDb = b2 - b1 ; }
                  }
                tg1 = wp->gene ;
                break ;
              }
        }
      vtxtPrintf (txt, "Transcribed_gene \"%s\"\nShedded_from \"%s\"\n\n"
                  , name (tg)
                  , name (tg1)) ;
    }
  for (ii = 0 ; ii < keySetMax ( principalTgs) ; ii++)
    {
      tg = keySet ( principalTgs, ii) ;
      vtxtPrintf (txt, "Transcribed_gene \"%s\"\n-D Shedded_from\n\n"
                  , name (tg)) ;
      
    }
  parseBuffer (vtxtPtr (txt), 0) ;

  vtxtDestroy (txt) ;
  keySetDestroy (tgs1) ;
  keySetDestroy (genes) ;
  keySetDestroy (principalTgs) ;
  keySetDestroy (tgPairs) ;
  return ;
} /* cdnaTagSheddedGenesOneCluster */

/***************************************************************/

int cdnaTagSheddedTranscribedGenes (KEYSET tgs0)
{
  int nn = 0 ;
  KEYSET tgs, ks0, ks1 ;
  KEY tg ;

  nn = cdnaTagFloating3pGenes (tgs0) ;

  tgs = query (tgs0, "to_be_fused_with") ;
  while (keySetMax (tgs))
    {
      tg = keySet (tgs, 0) ; /* get first candidate */
      ks1 = queryKeyTransitiveClosure (tg, str2tag ("to_be_fused_with")) ;
      cdnaTagSheddedGenesOneCluster (ks1) ;
      ks0 = keySetMINUS (tgs, ks1) ;
      keySetDestroy (tgs) ;
      keySetDestroy (ks1) ;
      tgs = ks0 ; ks0 = 0 ;
    }

  keySetDestroy (tgs) ;
  
  ks1 = query (tgs0, "Shedded_from") ;
  nn = keySetMax (ks1) ;
  keySetDestroy (ks1) ;

  return nn ; 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/


double fitResults(int chrom,double X1plusX2half)
{
        int i;
        double res,xDegree;
        double fitCoef[6][10]={
                {
                -17.73438,
                -7.91388E-6,
                1.21376E-11,
                -4.59871E-18,
                8.31368E-25,
                -7.32405E-32,
                1.77929E-39,
                1.90459E-46,
                -1.50777E-53,
                3.20403E-61
                },
                {
                -19.21157,
                8.56921E-6,
                -2.96269E-12,
                3.26213E-20,
                4.12595E-25,
                -1.42216E-31,
                2.20934E-38,
                -1.81749E-45,
                7.70021E-53,
                -1.32396E-60
                },
                {
                -27.37335,
                3.29293E-6,
                -7.22335E-12,
                8.74387E-18,
                -3.68703E-24,
                7.89053E-31,
                -9.59259E-38,
                6.71265E-45,
                -2.51818E-52,
                3.92031E-60
                },
                {
                -27.40893,
                -1.92325E-6,
                9.50652E-12,
                -3.50665E-18,
                5.75935E-25,
                -4.48769E-32,
                8.74399E-40,
                1.02247E-46,
                -6.9016E-54,
                1.28634E-61
                },
                {
                -20.04515,
                -3.94894E-6,
                4.36283E-12,
                -9.44693E-19,
                9.46066E-26,
                -5.08353E-33,
                1.80178E-40,
                -6.86314E-48,
                2.49302E-55,
                -4.02849E-63
                },
                {
                -19.503,
                -4.21145E-7,
                3.09507E-13,
                1.02339E-18,
                -4.87583E-25,
                1.00242E-31,
                -1.12285E-38,
                7.09936E-46,
                -2.3738E-53,
                3.25756E-61
                }
        };
        if(chrom>=(int)'1' && chrom<=(int)'5')chrom-=(int)'0';
        else if(chrom=='x' || chrom=='X')chrom=6;

        
        xDegree=1;
        for(i=0;i<10;i++){
                res+=fitCoef[chrom-1][i]*xDegree;
                xDegree*=X1plusX2half;
        }
        return xDegree;
}

/*
 [4/22/2002 16:27 "/Graph1" (2452386)]
Polynomial Regression for chrom1_gpos:
Y = A + B1*X + B2*X^2 + B3*X^3 + B4*X^4 + B5*X^5 
    + B6*X^6 + B7*X^7 + B8*X^8 + B9*X^9

Parameter   Value   Error
------------------------------------------------------------
A   -17.73438   0.47364
B1  -7.91388E-6 2.26494E-6
B2  1.21376E-11 2.96338E-12
B3  -4.59871E-18    1.69778E-18
B4  8.31368E-25 1E-20
B5  -7.32405E-32    1E-20
B6  1.77929E-39 1E-20
B7  1.90459E-46 1E-20
B8  -1.50777E-53    1E-20
B9  3.20403E-61 1E-20
------------------------------------------------------------

R-Square(COD)   SD  N   P
------------------------------------------------------------
0.9975  0.45087 96  <0.0001
------------------------------------------------------------


[4/22/2002 16:27 "/Graph2" (2452386)]
Polynomial Regression for chrom2_gpos:
Y = A + B1*X + B2*X^2 + B3*X^3 + B4*X^4 + B5*X^5 
    + B6*X^6 + B7*X^7 + B8*X^8 + B9*X^9

Parameter   Value   Error
------------------------------------------------------------
A   -19.21157   0.50893
B1  8.56921E-6  2.71877E-6
B2  -2.96269E-12    3.87751E-12
B3  3.26213E-20 2.28352E-18
B4  4.12595E-25 1E-20
B5  -1.42216E-31    1E-20
B6  2.20934E-38 1E-20
B7  -1.81749E-45    1E-20
B8  7.70021E-53 1E-20
B9  -1.32396E-60    1E-20
------------------------------------------------------------

R-Square(COD)   SD  N   P
------------------------------------------------------------
0.9975  0.39463 73  <0.0001
------------------------------------------------------------


[4/22/2002 16:27 "/Graph3" (2452386)]
Polynomial Regression for chrom3_gpos:
Y = A + B1*X + B2*X^2 + B3*X^3 + B4*X^4 + B5*X^5 
    + B6*X^6 + B7*X^7 + B8*X^8 + B9*X^9

Parameter   Value   Error
------------------------------------------------------------
A   -27.37335   0.47308
B1  3.29293E-6  1.92888E-6
B2  -7.22335E-12    2.74522E-12
B3  8.74387E-18 1.72417E-18
B4  -3.68703E-24    1E-20
B5  7.89053E-31 1E-20
B6  -9.59259E-38    1E-20
B7  6.71265E-45 1E-20
B8  -2.51818E-52    1E-20
B9  3.92031E-60 1E-20
------------------------------------------------------------

R-Square(COD)   SD  N   P
------------------------------------------------------------
0.99868 0.35974 95  <0.0001
------------------------------------------------------------


[4/22/2002 16:27 "/Graph4" (2452386)]
Polynomial Regression for chrom4_gpos:
Y = A + B1*X + B2*X^2 + B3*X^3 + B4*X^4 + B5*X^5 
    + B6*X^6 + B7*X^7 + B8*X^8 + B9*X^9

Parameter   Value   Error
------------------------------------------------------------
A   -27.40893   0.52985
B1  -1.92325E-6 2.73456E-6
B2  9.50652E-12 3.10587E-12
B3  -3.50665E-18    1.55023E-18
B4  5.75935E-25 1E-20
B5  -4.48769E-32    1E-20
B6  8.74399E-40 1E-20
B7  1.02247E-46 1E-20
B8  -6.9016E-54 1E-20
B9  1.28634E-61 1E-20
------------------------------------------------------------

R-Square(COD)   SD  N   P
------------------------------------------------------------
0.99695 0.48532 82  <0.0001
------------------------------------------------------------


[4/22/2002 16:27 "/Graph5" (2452386)]
Polynomial Regression for chrom5_gpos:
Y = A + B1*X + B2*X^2 + B3*X^3 + B4*X^4 + B5*X^5 
    + B6*X^6 + B7*X^7 + B8*X^8 + B9*X^9

Parameter   Value   Error
------------------------------------------------------------
A   -20.04515   0.19751
B1  -3.94894E-6 7.75596E-7
B2  4.36283E-12 9.32725E-13
B3  -9.44693E-19    4.54891E-19
B4  9.46066E-26 1E-20
B5  -5.08353E-33    1E-20
B6  1.80178E-40 1E-20
B7  -6.86314E-48    1E-20
B8  2.49302E-55 1E-20
B9  -4.02849E-63    1E-20
------------------------------------------------------------

R-Square(COD)   SD  N   P
------------------------------------------------------------
0.9994  0.18618 66  <0.0001
------------------------------------------------------------


[4/22/2002 16:28 "/Graph6" (2452386)]
Polynomial Regression for chromx_gpos:
Y = A + B1*X + B2*X^2 + B3*X^3 + B4*X^4 + B5*X^5 
    + B6*X^6 + B7*X^7 + B8*X^8 + B9*X^9

Parameter   Value   Error
------------------------------------------------------------
A   -19.503 2.39823
B1  -4.21145E-7 5.87152E-6
B2  3.09507E-13 5.07405E-12
B3  1.02339E-18 2.17928E-18
B4  -4.87583E-25    1E-20
B5  1.00242E-31 1E-20
B6  -1.12285E-38    1E-20
B7  7.09936E-46 1E-20
B8  -2.3738E-53 1E-20
B9  3.25756E-61 1E-20
------------------------------------------------------------

R-Square(COD)   SD  N   P
------------------------------------------------------------
0.9964  0.82115 95  <0.0001
------------------------------------------------------------


*/

