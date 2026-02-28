/*  File: dnaalign.c
 *  Author: Ulrich Sauvage (ulrich@kaa.crbm.cnrs-mop.fr)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  9 15:16 1998 (fw)
 * Created: Wed Jun  1 20:40:35 1994 (ulrich)
 *-------------------------------------------------------------------
 */


/* @(#)dnaalign.c	1.18 11/12/96 */

#include <ctype.h>
#include "acedb.h"
#include "topology.h"
#include "dnaalign.h"
#include "query.h"
#include "dna.h"
#include "a.h"
#include "pick.h"
#include "plot.h"
#include "lex.h"
#include "bs.h"
#include "session.h"
#include "tags.h"
#include "chrono.h"
#include "systags.h"
#include "classes.h"
#include "display.h"
#include "aligntools.h"
#include "call.h"
#include "restriction.h"
#include "acembly.h"
#include "mytime.h"

typedef struct { KEY k1, k2 ; } KEY_KEY ;

static KEYSET ksstep = 0 ;
static char u[9] ;
static KEY dnaAlignPaires (DEFCPT look, KEY key1, KEY key2, int *taux) ;

extern void fMapReDrawWindow (void) ;

/***************************************************************/
/*************************** Divers ****************************/
/***************************************************************/

void dnaAlignInit (void)
{ u[A_] = 0 ;
  u[T_] = 3 ;
  u[G_] = 1 ;
  u[C_] = 2 ;
  alignToolsInit () ;
}

/***************************************************************/

void dnaAlignDestroy (DEFCPT look)
{ alignToolsDestroy (look) ;
  dnaDispGraph (0, 0) ; /* should destroy ksstep */
  if (keySetExists(ksstep))
    keySetDestroy(ksstep) ;
  if (look) look->step = 0 ;
}

/***************************************************************/

void dnaAlignForget(DEFCPT look, KEY key)
{
  monDnaForget(look, key) ;
} 

/***************************************************************/

char* dnaAlignDecodeOligo(KEY key)
{ unsigned int rac ;
  char a[13] ;
  char *cp ;

  rac = (unsigned int)key ;
  rac |= (1 << 24) ;
  cp = a + 12 ;
  while(cp--, rac > 1)
    { switch(rac & 3)
	{
	case 0:
	  *cp = (char)A_ ;
	  break ;
	case 1:
	  *cp = (char)G_ ;
	  break ;
	case 2:
	  *cp = (char)C_ ;
	  break ;
	case 3:
	  *cp = (char)T_ ;
	  break ;
	}
      rac >>= 2 ;
    }
  return dnaDecodeString(a) ;
}

/***************************************************************/
/*
int dnaAlignGiveNbBad(void)
{
  return keySetExists(badSeq) ? keySetMax(badSeq) : 0 ;
}
*/
/***************************************************************/

void dnaAlignAddSeqIn (DEFCPT look, KEY key)
{ KEY key1 ;
  KEYSET ks ;

  if (!look || !(ks = look->def))
    return ;
  if (class(key) == _VDNA && monDnaGet (look, 0, key))
    keySet(ks, keySetMax(ks)) = key ;
  else if (class(key) == _VSequence && lexReClass(key, &key1, _VDNA) &&
	   monDnaGet (look, 0, key1))
    keySet(ks, keySetMax(ks)) = key1 ;
  else
    messcrash("Bad class(key) or DNA too short in dnaAlignAddSeqIn") ;
}

/***************************************************************/

static int dnaAlignMaillonOrder(const void *a, const void *b)
{ return
    ((const KEY_KEY*)a)->k2 > ((const KEY_KEY*)b)->k2 ? -1 :
      ((const KEY_KEY*)a)->k2 == ((const KEY_KEY*)b)->k2 ? 0 : 1 ;
}

/***************************************************************/

static void dnaAlignSortMaillon (DEFCPT look)
{ int i, j = 0, imax ;
  Array temp, tempbis, maillon = look->maillon ;
  TREE_DEF *ksmp ;
  KEY_KEY *keyp ;

  temp = arrayCreate(128, KEY_KEY) ;
  imax = arrayMax(maillon) ;
  ksmp = arrp(maillon, 0, TREE_DEF) - 1 ;
  while(ksmp++, imax--)
    { i = keySetMax(ksmp->ks) ;
      while(i--)
	{ arrayp(temp, j, KEY_KEY)->k1 = keySet(ksmp->ks, i) ;
	  tempbis = monDnaGet (look, 0, keySet(ksmp->ks, i)) ;
	  arrp(temp, j++, KEY_KEY)->k2 = arrayMax(tempbis) ;
	}
      keySetDestroy(ksmp->ks) ;
      keySetDestroy(ksmp->x) ;
    }
  arraySort(temp, dnaAlignMaillonOrder) ;
  arrayReCreate(maillon, 1, TREE_DEF) ;
  ksmp = arrayp(maillon, 0, TREE_DEF) ;
  ksmp->ks = keySetCreate() ;
  ksmp->x = keySetCreate() ;
  imax = arrayMax(temp) ;
  keyp = arrp(temp, 0, KEY_KEY) ;
  for (i = 0 ; i < imax ; keyp++, i++)
    { keySet(ksmp->ks, i) = keyp->k1 ;
      keySet(ksmp->x, i) = keyp->k2 ;
    }
  arrayDestroy(temp) ;
}

/***************************************************************/

void dnaAlignTryAll(DEFCPT look)
{ int i ;
  KEYSET ksp ;

  dnaAlignSortMaillon (look) ;
  ksp = arrayp(look->maillon, 0, TREE_DEF)->x ;
  i = keySetMax(ksp) ;
  while(i--)
    keySet(ksp, i) = 0 ;
  dnaAlignCptSegments(look) ;
}

/***************************************************************/

void dnaAlignNewTryPaires(DEFCPT look)
{ int imax, i, j, taux, tauxmax = look->taux ;
  KEY key, key1, key2, *keyp, *keyq ;
  KEYSET def = look->def ;
  Array maillon = look->maillon ;

  dnaAlignSortMaillon (look) ;
  imax = keySetMax(arrp(maillon, 0, TREE_DEF)->ks) ;
  for (i = 0 ; i < imax - 1 ; i++) /* tour sur maillon stop quand les sequences trop petites */
    { if (keySet(arrp(maillon, 0, TREE_DEF)->x, i) < 1000)
	break ;
      if (!(key1 = keySet(arrp(maillon, 0, TREE_DEF)->ks, i)))
	continue ;
      for (j = i + 1 ; j < imax ; j++)
	{ if (!(key2 = keySet(arrp(maillon, 0, TREE_DEF)->ks, j)))
	    continue ;
	  key = dnaAlignPaires(look, key1, key2, &taux) ;
	  if (key && taux <= tauxmax)
	    { key1 = key ;
	      keySet(arrp(maillon, 0, TREE_DEF)->ks, j) = 0 ;
	      j = i ;
	    }
	}
      keySet(def, keySetMax(def)) = key1 ;
    }
  j = keySetMax(def) ;
  while(i < imax) /* add the last one */
    if ((key = keySet(arrp(maillon, 0, TREE_DEF)->ks, i++)))
      keySet(def, j++) = key ;
  ksstep = keySetReCreate(ksstep) ; /* display */
  keyp = arrp(def, 0, KEY) ;
  array(ksstep, j - 1, KEY) = 0 ;
  keyq = arrp(ksstep, 0, KEY) ;
  while(j--)
    lexReClass(*keyp++, keyq++, _VSequence) ;
  if (look->display)
    dnaDispGraph(ksstep, look->tour) ;
}

/***************************************************************/
/************************** Get  Data **************************/
/***************************************************************/

static void dnaMakeMotif (DEFCPT look, Associator adna, int nbolig)
{ int i, j, max, frac, nd, nf, index = 0 ;
  unsigned int rac, gtr, n ;
  int debut, fin ;
  char *vp, *vp0 = 0 ;
  KEY key ;
  KEYSET dna = look->def, mar = look->mar ;
  Array mondna ;
  char *v0, *vq ;

  i = arrayMax(dna) ;
  frac = ((nbolig + 1) / 2) * 2 ;
  while (i--)
    {
      mondna = monDnaGet (look, 0, keySet(dna, i)) ;
      if (!mondna) continue ;
      debut = - 1 ;
      fin = arrayMax(mondna) ;
      max = fin - 1 ;
      n = nbolig ;
      nd = nf = 0 ;
      while (n--)
	{
	  if (n % 2)
	    { j = max * (frac - nf++) / frac ;
	      if (j >= fin)
		j = fin - 1 ;
	      if (j < DODECA || 
		  !alignToolsMakeOligo(mondna, j, 0, DODECA, &rac, &fin))
		continue ;
	    }
	  else
	    {
	      j = max * nd++ / frac ;
	      if (j <= debut)
		j = debut + 1 ;
	      if (j > max - DODECA ||
		  !alignToolsMakeOligo(mondna, j, max, DODECA, &rac, &debut))
		continue ;
	    }
	  v0 = vp0 + rac ;
	  if (!(assFind(adna, v0, &vp)))
	    {
	      key = (KEY)rac ;
	      vq = vp0 + key ;
	      keySet(mar, index++) = key ;
	      assInsert(adna, v0, vq) ;
	      gtr = dnaComplement(rac, NB_OCTET) ;
	      v0 = vp0 + gtr ;
	      if (gtr != rac) /* same association for one oligo and its complement */
		assInsert(adna, v0, vq) ;
	    }
	}
    }
  keySetSort(mar) ;
}

/***************************************************************/

void dnaAlignGetData (DEFCPT look)
/* KEYSET dna, KEYSET mar, Array marin, Array defin, int nbolig */
{ 
  int i, j, n, max, nbolig, newnumber ;
  unsigned int mask = 1 << 24 ; /* pour avoir des dodecameres */
  KEY dummy ; 
  KEYSET ksp, ksq ; 
  KEYSET ks, *kstmp, dna, mar ;
  char *cdep, buf[128] ;
  unsigned int test ;
  BOOL depart ;
  Array histdna = 0, mondna, marin, defin ;
  char *vp, *vp0 = 0 ;
  Associator adna = 0 ;
#ifndef NON_GRAPHIC 
  int histtot ;
  void *vdummy ;
  KEYSET aa = 0 ; 
  KEY key, key1 ;
#endif

   /* protection de mondnaget */
  if (look->def)
    { ks = keySetCopy (look->def) ;
      j = 0 ;
      for (i = 0 , j = 0 ;  i < keySetMax (ks) ; i++)
	{ mondna = monDnaGet (look, 0, keySet(ks, i)) ;
	  if (mondna)
	    keySet (look->def, j++) = keySet(ks, i) ;
	}
      keySetMax (look->def) = j ;
    }

  dna = look->def ;
  mar = look->mar ;
  marin = look->marInDef ;
  defin = look->defInMar ;
  nbolig = (int)look->whatDis ;


#ifndef NON_GRAPHIC /* il faut passer par load file qui cree le keyset */
  if (!keySetMax(dna))
    {
      if (!look->rejected)
	{
	  if (!keySetActive(&aa, &vdummy))
	    { messout("First select a keyset containing sequences") ;
	      return ;
	    }
	  j = 0 ;
	  for (i = 0 ; i < keySetMax(aa) ; i++)
	    { key = keySet(aa,i) ;
	      if (class(key) == _VDNA && monDnaGet (look, 0, key))
		keySet(dna,j++) = key ;
	      else if (class(key) == _VSequence &&
		       lexReClass(key, &key1, _VDNA) && monDnaGet (look, 0, key1))
		keySet(dna,j++) = key1 ;
	    }
	  if (!keySetMax(dna))
	    { messout("First select a keyset containing dna") ;
	      return ;
	    }
	  look->tour = 0 ;
	}
      else /* a corriger pour que cela soit impossible de faire Iterate dans ce cas */
	{ messout("Please ReInsert Bad Sequences") ;
	  return ;
	}
    }
#endif /* attention au cas que des bads ! attention au contexte */
  keySetSort(dna) ;
  keySetCompress (dna) ;
  i = keySetMax(dna) ;
  adna = assBigCreate(2*i) ;
  if (nbolig > 1 && i * nbolig > 15000)
    { newnumber = 15000 / i ;
      if (!newnumber)
	newnumber = 1 ;
      if (look->manEntry)
	{ sprintf(buf, "You are asking for too many oligos (%d), please confirm", nbolig) ;
	ici:
	  if (messPrompt(buf, messprintf("%d", newnumber), "i"))
	    { freeint(&j) ;
	      if (j > nbolig)
		goto ici ;
	      newnumber = j ;
	    }
	}
      nbolig = newnumber ;
      look->whatDis = (unsigned char)nbolig ;
      sprintf(look->nboligo, "%d", nbolig) ;
    }
  chrono("dnaMakeMotif") ;
  dnaMakeMotif(look, adna, nbolig) ;
  if (!(j = keySetMax(mar)))
    { messout("I can't find an oligo") ;
      chronoReturn() ;
      return ;
    }
  chronoReturn () ;
  
  chrono("Matrice Transposee") ;
  kstmp = arrayp(marin, i - 1, KEYSET) + 1 ;
  while(kstmp--, i--)
    *kstmp = keySetCreate() ;
  kstmp = arrayp(defin, j - 1, KEYSET) + 1 ;
  while(kstmp--, j--)
    *kstmp = keySetCreate() ;
  i = keySetMax(dna) ;
#ifndef NON_GRAPHIC
  histdna = arrayCreate(1000, int) ;
  histtot = 0 ;
#endif
  while (i--)/* Creation de la matrice transposee */
    { mondna = monDnaGet (look, 0, keySet(dna,i)) ;
      max = arrayMax(mondna) ;
#ifndef NON_GRAPHIC      
      array(histdna, max, int)++ ;
      histtot += max ;
#endif
      test = 1 ;
      depart = FALSE ;
      cdep = arrp(mondna, 0, char) - 1 ;
      ksp = arr(marin, i, KEYSET) ;
      n = 0 ;
      while (cdep++, max--)
	{
	  switch(*cdep)
	    {
	    case A_: case T_: case G_: case C_:
	      test = ((test << 2) | u[(int)*cdep]) ;
	      break ;
	    default:
	      depart = FALSE ;
	      test = 1 ;
	      continue ;
	    }
	  if (!depart)
	    { if (test & (mask >> 2))
		depart = TRUE ;
	      continue ;
	    }
	  test &= (mask - 1) ;
	  if (!assFind(adna, vp0 + test, &vp))
	    continue ;
	  dummy = vp - vp0 ;
	  if (!keySetFind(mar, dummy, &j))
	    messcrash("Problem keySetFind in dnaAlignGetData") ;
	  keySet(ksp, n++) = j ; /* def i contient marqueur j */
	  ksq = arr(defin, j, KEYSET) ;
	  keySet(ksq, keySetMax(ksq)) = i ; /* marqueur j dans def i */
	}
    }
  assDestroy(adna) ;
  i = arrayMax(marin) ;
  kstmp = arrp(marin, i - 1, KEYSET) + 1 ;
  while(kstmp--, i--)
    { keySetSort(*kstmp) ;
      keySetCompress(*kstmp) ;
    }
  i = arrayMax(defin) ;
  kstmp = arrp(defin, i - 1, KEYSET) + 1 ;
  while(kstmp--, i--)
    { keySetSort(*kstmp) ;
      keySetCompress(*kstmp) ;
    }
  chronoReturn() ;
#ifdef NON_GRAPHIC
  arrayDestroy (histdna) ;
#else
  if (look->display)
    plotHisto(messprintf("Number of sequences per length (total length %d)", 
			 histtot), histdna) ;
  else 
    arrayDestroy (histdna) ;
#endif
}

/***************************************************************/
/***************************************************************/

static KEYSET dnaAlignFilter (DEFCPT look, KEYSET ks)
{ KEYSET ks1 = keySetCreate () ;
  KEY key ;
  int i, j = 0 ;

  for (i = 0 ; i < keySetMax(ks) ; i++)
    { key = keySet(ks, i) ; 
      if (class (key) == _VSequence)
	lexReClass (key, &key, _VDNA) ;
      if (class(key) == _VDNA && monDnaGet (look, 0, key))
	keySet(ks1, j++) = key ;
    }
  keySetSort (ks1) ;
  keySetCompress (ks1) ;
  return ks1 ;
}

/***************************************************************/

int dnaAlignLoad (DEFCPT look, char *cp)
{ KEYSET ks = 0, ks1, ks2, ks3 ;
  int i = 0 ;
  BOOL detruit = TRUE ;
  KEYSET baseCallNewScf (void) ;

  if (!cp)
    ks = baseCallNewScf () ;
  else if (!strncmp (cp, "-active", 7))
    { 
      if (keySetExists (look->actif))
	{ ks = look->actif ;
	  detruit = FALSE ;
	}
      else 
	{ messout ("Please select a keySet") ;
	  return 0 ;
	}
    }
  else
    ks = query (0, cp) ;
  ks2 = query (ks, ">Subsequence") ;
  ks3 = keySetOR (ks, ks2) ;
  keySetDestroy (ks2) ;
  if (look->def)
    { keySetSort (look->def) ;
      i = keySetMax (look->def) ;
      ks1 = keySetOR (look->def, ks3) ;
      keySetDestroy (look->def) ;
      look->def = dnaAlignFilter (look, ks1) ;
      keySetDestroy (ks1) ;
    }
  else
    { if (detruit)   /* kludge pour le fuse -> on ne veut pas avoir le meme link */
	{ ks2 = query (ks, "CLASS Assembly") ;
	  if (keySetMax (ks2) == 1)
	    defCptChangeLook (look, 0, keySet (ks2, 0)) ;
	  keySetDestroy (ks2) ;
	}
      look->def = dnaAlignFilter (look, ks3) ;
    }
  if (detruit)
    keySetDestroy (ks) ;
  keySetDestroy (ks3) ;
  return keySetMax (look->def) - i ;
}

/***************************************************************/
/************************ Display Tools ************************/
/***************************************************************/

#ifdef NON_GRAPHIC
void dnaDispGraph(Array ks, int tour)
{ if (ks)
    { keySetSort(ks) ;
      keySetCompress(ks) ;
    }
}

/***************************************************************/

#else
static void *dispHandle = 0 ;

void dnaDispGraph (Array ks, int tour)
{ static Graph dispGraph = 0 ;

  if (!ks)
    { if (graphExists (dispGraph))
	{ graphActivate(dispGraph) ;
	  graphDestroy() ;
	}
      return ;
    }
  keySetSort (ks) ;
  keySetCompress (ks) ;
  if (graphExists (dispGraph))
    {
      keySetShow (ks, dispHandle) ;
      graphRetitle(messprintf("Step %d", tour)) ;
    }
  else
    { 
      dispGraph = displayCreate(DtKeySet) ;
      graphRetitle(messprintf("Step %d", tour)) ;
      dispHandle = keySetShow (ks, 0) ;
    }
}
#endif /* non graphic */

/***************************************************************/

    /* these are sets in class DNA */
static void dnaAlignRemoveBad(DEFCPT look)
{
  int i = look->tour - 3 ;
  KEY *ksp, *ksq ;
  KEYSET ksInter = 0, ksInter2 = 0, ksInter3 = 0 ;

  if (look->step < 3 || i < 0)
    return ;
  if (!(look->rejected))
    look->rejected = keySetCreate() ;
  keySetSort(look->def) ;
  keySetCompress(look->def) ;
  ksInter = query(look->def, messprintf
		  ("IS = _segment_%d.* AND IS > _segment_%d.%d", look->id, look->id, i)) ;
  ksInter2 = keySetMINUS(look->def, ksInter) ;
  ksInter3 = keySetOR(ksInter2, look->rejected) ;
  keySetDestroy(ksInter2) ;
  keySetDestroy(look->rejected) ;
  look->rejected = ksInter3 ;
  keySetReCreate(look->def) ;
  i = keySetMax(ksInter) ;
  if (i)
    { keySet(look->def, i - 1) = 0 ;
      ksp = arrp(look->def, 0, KEY) ;
      ksq = arrp(ksInter, 0, KEY) ;
      while(i--)
	*ksp++ = *ksq++ ;
    }
  keySetDestroy(ksInter) ;
}

/***************************************************************/

/* utilisation interne on garde badSeq */
static void dnaAlignDoInsertBad(KEYSET badSeq, KEYSET ksDna)
{
  int i ;
  KEYSET ksInter = 0 ;
  KEY *ksp, *ksq ;

  if (!badSeq)
    return ;
  keySetSort(ksDna) ;
  keySetCompress(ksDna) ;
  ksInter = keySetOR(ksDna, badSeq) ;
  ksDna = keySetReCreate(ksDna) ;
  i = keySetMax(ksInter) ;
  keySet(ksDna, i - 1) = 0 ;
  ksp = arrp(ksDna, 0, KEY) ;
  ksq = arrp(ksInter, 0 , KEY) ;
  while(i--)
    *ksp++ = *ksq++ ;
  keySetDestroy(ksInter) ;
}

/***************************************************************/

/* utilisation interne pour ksstep */
static void dnaAlignDoInsertBadSeq(KEYSET badSeq, KEYSET ksSeq)
{
  int i ;
  KEYSET ksInter = 0 ;
  KEY *ksp, *ksq ;

  ksInter = query (ksSeq, ">DNA") ;
  dnaAlignDoInsertBad(badSeq, ksInter) ;
  ksSeq = keySetReCreate(ksSeq) ;
  i = keySetMax(ksInter) ;
  keySet(ksSeq, i - 1) = 0 ;
  ksp = arrp(ksSeq, 0, KEY) ;
  ksq = arrp(ksInter, 0 , KEY) ;
  while(i--)
    lexReClass (*ksq++, ksp++, _VSequence) ;
  keySetDestroy(ksInter) ;
}

/***************************************************************/

/* utilisation public on recupere et on detruit badSeq */
/* sets in class Sequence */
void dnaAlignInsertBad (DEFCPT look)
{ int i ;
  KEY *ksp, *ksq ;

  dnaAlignDoInsertBad (look->rejected, look->def) ;
  keySetDestroy (look->rejected) ;
  ksstep = keySetReCreate (ksstep) ;
  i = keySetMax (look->def) ;
  keySet(ksstep, i - 1) = 0 ;
  ksp = arrp (ksstep, 0, KEY) ;
  ksq = arrp (look->def, 0 , KEY) ;
  while(i--)
    lexReClass (*ksq++, ksp++, _VSequence) ;
  dnaDispGraph(ksstep, look->tour) ;
}

/***************************************************************/
/************************ Align Segment ************************/
/***************************************************************/

static void dnaNewSeg(TREE_DEF ksmi, TREE_DEF *ksm, int *j, BOOL neverCut)
{
  int ksmax, xmax, segmax, i ;
  static int last, nb ;

  if (!keySetExists(ksm->ks))
    { ksm->ks = keySetCreate() ;
      ksm->x = keySetCreate() ;
      last = 0 ;
      nb = 0 ;
    }
  segmax = keySetMax(ksm->ks) ;
  if (*j == last + 1)
    nb++ ;
  else 
    nb = 1 ;
  if (nb < 4 || neverCut)
    { keySet(ksm->ks, segmax) = keySet(ksmi.ks, *j) ;
      keySet(ksm->x, segmax) = keySet(ksmi.x, *j) ;
      last = *j ;
      (*j)++ ;
    }
  else
    { ksmax = keySetMax(ksmi.ks) ;
      xmax = keySet(ksmi.x, *j) ;
      for (i = *j ; i < ksmax && keySet(ksmi.x, i) >= xmax ; i++)
	{ keySet(ksm->ks, segmax) = keySet(ksmi.ks, i) ;
	  keySet(ksm->x, segmax++) = keySet(ksmi.x, i) ;
	}
      *j = i ; /* indice du dernier rejete + 1 */
    }
}

/***************************************************************/

static void dnaShift(Array dna, int *debut)
{ int i, nb = 1000 ;
  char *cp, *cq ;

  i = arrayMax(dna) ;
  if ((*debut + nb) < 0)
    nb = -(*debut) + 500 ;
  *debut += nb ;
  array(dna, i + nb - 1, char) = 0 ;
  cp = arrp(dna, i, char) ;
  cq = cp + nb ;
  while (cp--, cq--, i--)
    *cq = *cp ;
  cq++ ;
  while(cq--, nb--)
    *cq = 0 ;
}

/***************************************************************/

static void dnaAddPiece(Array mydna, Array dnasup, Array erreur, int sens,
			int *x1, int *x2, int *shift)
{ int i, errpos, errtype, nerrmax = 0, nerr ;
  int debut, fin, taille ;
  static char bu[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } ;
  char base[16] ;
  char *cp, *cp0, *cq ;

  chrono("dnaAssBilan") ;

  if (!bu[2])
    { i = 16 ;
      while(i--)
	bu[i] = 0 ;
      bu[8] = 8 ;
      bu[4] = 4 ;
      bu[2] = 2 ;
      bu[1] = 1 ;
      if (A_ != 1 || T_ != 2 || G_ != 4 || C_ != 8)
	messcrash("Incompatible values of A_, T_, G_, C_ in dnaalign.c") ;
    }
  else if (!bu[3])
    { i = 15 ;
      while(i--)
	bu[i] = i ;
    }
  debut = *x1 ;
  if (debut < 0)
    {
      *shift = - debut ;
      dnaShift(mydna, &debut) ;
      *shift += debut ;
    }
  else *shift = 0 ;
  if (erreur)
    nerrmax = arrayMax(erreur) ;
  nerr = 0 ; errtype = 0 ;
  if (nerrmax)
    { if ((fin = *shift))
	{ i = nerrmax ;
	  while(i--)
	    arr(erreur, i, A_ERR).iLong += fin ;
	}
      errpos = arr(erreur, nerr, A_ERR).iLong ;
      errtype = arr(erreur, nerr++, A_ERR).type ;
    }
  else
    errpos = -1 ;
  taille = arrayMax(dnasup) ;
  i = 16 ;
  if (sens == 1) /* dnasup dans le sens direct */
    { cq = arrp(dnasup, 0, char) - 1 ;
      while(i--)
	base[i] = bu[i] ;
    }
  else /* dnasup dans le sens reverse */
    { cq = arrp(dnasup, taille - 1, char) + 1 ;
      while(i--)
	base[i] = complementBase(bu[i]) ;
    }
  fin = debut + taille + 2 * nerrmax ;
  if (fin >= arrayMax(mydna))
    array(mydna, fin, char) = 0 ;
  cp0 = arrp(mydna, 0, char) ;
  cp = cp0 + debut - 1 ;
  while(cp++, cq += sens, taille--)
    {
      if (cp - cp0 == errpos)
	{ switch(errtype)
	    {
	    case ERREUR: case AMBIGUE:
	      break ;
	    case TROU: /* trou dans la sequence ajoutee */
	      cq -= sens ;
	      taille++ ;
	      break ;
	    case TROU_DOUBLE:
	      cq -= 2 * sens ;
	      taille += 2 ;
	      break ;
	    case TROU_TRIPLE:
	      cq -= 3 * sens ;
	      taille += 3 ;
	      break ;
	    case INSERTION: /* trou dans le bilan */ /* ancien apres */
	      cq += sens ;
	      if (taille)
		taille-- ;
	      break ;
	    case INSERTION_DOUBLE:
	      cq += 2 * sens ;
	      taille -= 2 ;
	      if (taille < 0)
		taille = 0 ;
	      break ;
	    case INSERTION_TRIPLE:
	      cq += 3 * sens ;
	      taille -= 3 ;
	      if (taille < 0)
		taille = 0 ;
	      break ;
	    }
	  if (nerr < nerrmax)
	    { errpos = arr(erreur, nerr, A_ERR).iLong ;
	      errtype = arr(erreur, nerr++, A_ERR).type ;
	    }
	  else errpos = - 1 ;
	}
      else if (!*cp || *cp == N_)
	*cp = base[(int)*cq] ;
    }
  *x1 = debut ;
  *x2 = cp - 1 - cp0 ;
}

/***************************************************************/

static void doCptSegment(DEFCPT look, TREE_DEF seg, TREE_DEF *ksnew, Stack bilan, 
			 Array mydna, int taux, int *nbseg, int *indice, 
			 int *x1, int *x2, BOOL neverCut)
{ KEY key ;
  int xTree ;
  Array mondna ;
  Array erreur = 0 ;
  int y1, y2, sens, monmax, shift, nbessai = taux > 15 ? 17 : taux + 2 ;
  static int shifttot = 0 ;
  static int *min, *max ;
  BOOL control = FALSE ;

  key = keySet(seg.ks, *indice) ;
  xTree = keySet(seg.x, *indice) ;
  if (!*indice) /* pour prendre tout l'arbre en cas de coupure */
    xTree = 0 ; /* a modifier en faisant des vrais sous arbres */
  mondna = monDnaGet (look, 0, key) ;
  monmax = arrayMax(mondna) - 1 ;
  if (!*indice)
    { y1 = 2000 ;
      sens = 1 ;
      min = x1 ;
      max = x2 ;
      *min = y1 ;
      *max = y1 + monmax ;
    }
  else
    {
      if (!(erreur = alignToolsMatch(mydna, *min, *max, *x1, *x2, mondna, 0, monmax,
				     0, monmax, &y1, &y2, &sens, 0, taux, 
				     nbessai, control)))
	{ dnaNewSeg(seg, ksnew, indice, neverCut) ;
	  arrayDestroy(erreur) ;
	  return ;
	}
    }
  dnaAddPiece(mydna, mondna, erreur, sens, &y1, &y2, &shift) ;
  arrayDestroy(erreur) ;
  (*nbseg)++ ;
  if (shift)
    { *min += shift ;
      *max += shift ;
      shifttot += shift ;
      push(bilan, shift, int) ;
    }
  shift = shifttot ;
  if (y1 < *min)
    *min = y1 ;
  if (y2 > *max)
    *max = y2 ;
  (*indice)++ ;
  while(*indice < keySetMax(seg.ks) && xTree <= keySet(seg.x, *indice))
    doCptSegment(look, seg, ksnew, bilan, mydna, taux, nbseg, indice, &y1, &y2, neverCut) ;
  y1 = y1 + shifttot - shift + 1 ; /* + 1 pour etre en coordonnee bio 1 -> max */
  y2 = y2 + shifttot - shift + 1 ; /* et non plus info 0 -> max - 1 */
  if (sens == 1)
    { push(bilan, y2, int) ;
      push(bilan, y1, int) ;
    }
  else
    { push(bilan, y1, int) ;
      push(bilan, y2, int) ;
    }
  push(bilan, key, KEY) ;
  push(bilan, 1 << 30, int) ;
}

/***************************************************************/

static TREE_DEF cptThisSegments (DEFCPT look, TREE_DEF seg, Stack bilan, Array mydna, 
				 int taux, int *nbseq, BOOL neverCut)
{ int x1, x2, ind = 0, i ;
  char *cp, *cq ;
  TREE_DEF ksnew ;

  ksnew.ks = ksnew.x = 0 ;
  *nbseq = 0 ;
  doCptSegment(look, seg, &ksnew, bilan, mydna, taux, nbseq, &ind, &x1, &x2, neverCut) ;
  if (*nbseq == 1) /* pas de changement donc suite inutile */
    return ksnew ;
  cp = arrp(mydna, 0, char) ;
  cq = cp + x1 ;
  i = x2 - x1 + 1 ;
  arrayMax(mydna) = i ;
  push(bilan, - x1, int) ;
  while(i--)
    *cp++ = *cq++ ;
  return ksnew ;
}

/***************************************************************/

void dnaAlignCptSegments(DEFCPT look)
{ /* attention taux d'erreur a choisir dans defcpt */
  int i, imax, max, istep, segstep, shift, shift2, shift3, taux, clip ;
  TREE_DEF ksi, ksnew ;
  KEY mykey, mykey2, mykey3, dnakey3 ;
  KEYSET dna ;
  char buff[20] ;
  Array mydna = 0, segments, units = 0 ;
  int nochange, maxu ;
  Stack dnaBilan = 0 ;
  OBJ obj, obj2 = 0 ;

  messStatus("Assembling Segments") ;
  if (!look)
    return ;
  if (!look->tour) /* efface les _segments */
    alignToolsDestroy_Segs(look->id) ;
  look->assembling = TRUE ;
  dna = look->def ;
  segments = look->maillon ;
  taux = look->taux ;
  look->step++ ;
  look->tour++ ;
  if (keySetExists(ksstep))
    keySetDestroy(ksstep) ;
  ksstep = keySetCreate() ;
  istep = segstep = 0 ;
  for (i = 0 ; i < arrayMax(segments) ; i++)
    { imax = arrayMax(segments) ;
      ksi = arr(segments, i, TREE_DEF) ;
      nochange = 0 ;
      if (keySetMax(ksi.ks) > 1)
	{ dnaBilan = stackCreate(6 * keySetMax(ksi.ks)) ;
	  mydna = arrayCreate(5000, char) ;
	  ksnew = cptThisSegments(look, ksi, dnaBilan, mydna, taux, &nochange, FALSE) ;
	  if (ksnew.ks)
	    array(segments, imax, TREE_DEF) = ksnew ;
	}
      if (!nochange || nochange == 1) /* pas de changement dans le segment */
	{ mykey = keySet(ksi.ks, 0) ;
	  keySet(dna, i) = mykey ;
	  lexReClass(mykey, &mykey2, _VSequence) ;
	  keySet(ksstep, istep++) = mykey2 ;
	  if (nochange) /* destroy if create and not the previous mydna */
	    arrayDestroy(mydna) ;
	}
      else /* nouveau segment */
	{ sprintf(buff, "_segment_%d.%d.%d", look->id, look->tour, segstep++) ;
	  lexaddkey(buff, &mykey, _VDNA) ;
	  keySet(dna, i) = mykey ;
	  max = arrayMax(mydna) ;
	  monDnaInsert (look, mykey, mydna) ;
	  lexaddkey(buff, &mykey2, _VSequence) ;
	  keySet(ksstep, istep++) = mykey2 ;
	  if ((obj = bsUpdate(mykey2)))
	    { bsAddKey(obj, _DNA, mykey) ;
	      bsAddData(obj, _bsRight, _Int, &max) ;
	      if (bsFindTag(obj, _Assembled_into))
		bsRemove(obj) ;
	      shift = 0 ;
	      maxu = 0 ;
	      units = arrayCreate (500, BSunit) ;
	      while(nochange)
		{ 
		  shift2 = pop(dnaBilan, int) ;
		  if (shift2 == (1 << 30))
		    { dnakey3 = pop (dnaBilan, KEY) ;
		      lexReClass(dnakey3, &mykey3, _VSequence) ;
		      array (units, maxu++, BSunit).k = mykey3 ;
		      shift3 = pop(dnaBilan, int) + shift ;
		      array (units, maxu++, BSunit).i = shift3 ;
		      shift3 = pop(dnaBilan, int) + shift ;
		      array (units, maxu++, BSunit).i = shift3 ;
		      if (!(obj2 = bsCreate (mykey3)) || 
			  !bsGetData (obj2, _Clipping, _Int, &clip))
			clip = 1 ;
		      bsDestroy (obj2) ;
		      array (units, maxu++, BSunit).i = clip ;
		      clip += arrayMax (monDnaGet (look, mykey2, dnakey3)) - 1 ;
		      array (units, maxu++, BSunit).i = clip ;
		      nochange-- ;
		    }
		  else
		    shift += shift2 ;
		}
	      if (!bsAddArray (obj, _Assembled_from, units, 5))
		messerror ("Pb in constructing segment") ;
	      bsSave(obj) ;
	      arrayDestroy (units) ;
	    }
	}
      stackDestroy(dnaBilan) ;
      keySetDestroy (ksi.ks) ;
      keySetDestroy (ksi.x) ;
    }
  if (look)
    { dnaAlignRemoveBad (look) ;
      dnaAlignDoInsertBadSeq (look->rejected, ksstep) ;
      arrayDestroy (look->maillon) ;
    }
  if (look->display)
    dnaDispGraph (ksstep, look->tour) ;
}

/***************************************************************/

static void dnaAlignUpdateSubSeqPos (OBJ obj, KEY contig, int x, KEY tag)
{ Array seq ;
  OBJ obj2 = bsCreate (contig) ;
  int i, y ;
  BSunit *u ;

  seq = arrayCreate (50, BSunit) ;
  if (bsFindTag (obj2, tag) && bsFlatten (obj2, 5, seq))
    for (i = 0 ; i < arrayMax (seq) ; i += 5)
      {	u = arrp (seq, i, BSunit) ;
	if (bsFindKey (obj, tag, u[0].k))
	  { y = u[1].i + x ;
	    bsAddData (obj, _bsRight, _Int, &y) ;
	    y = u[2].i + x ;
	    bsAddData (obj, _bsRight, _Int, &y) ;
	    y = u[3].i ;
	    if (y)
	      { bsAddData (obj, _bsRight, _Int, &y) ;
		y = u[4].i ;
		bsAddData (obj, _bsRight, _Int, &y) ;
	      }
	  }
      }
  bsDestroy (obj2) ;
  arrayDestroy (seq) ;
}

/***************************************************************/

static void dnaAlignRemoveTag (KEY key, KEY tag)
{ OBJ obj = bsUpdate (key) ;

  if (obj)
    { if (bsFindTag (obj, tag))
	bsRemove (obj) ;
      bsSave (obj) ;
    }
}

/***************************************************************/

void dnaAlignAddKeySetIn (KEY link, KEY contig, KEYSET reads, int taux)
{ int i, j, nochange = 0, maxdna, shift, shift2, shift3, x, clip ;
  KEY dnacontig, key, dummy, *keyq, dnakey ;
  TREE_DEF ksi, ksnew ;
  Stack dnaBilan = 0 ;
  OBJ obj, obj2 = 0 ;
  Array mydna, newseq = 0 ;
  DEFCPT look = defCptGetLook (link) ;

  i = keySetMax (reads) ;
  if (!i)
    return ;
  messStatus ("Assembling Reads") ;
  dnaAlignRemoveTag (contig, _Previous_contig) ;
  ksi.ks = arrayCreate (i + 1, KEY) ;
  ksi.x =  arrayCreate (i + 1, KEY) ;
  if (!lexReClass (contig, &dnacontig, _VDNA))
    messcrash ("Bad Type of Key contig") ;
  keySet (ksi.ks, 0) = dnacontig ;
  keyq = arrp (reads, 0, KEY) - 1 ;
  j = 1 ; /* already one seq in the keySet */
  while (keyq++, i--)
    if (lexReClass (*keyq, &dummy, _VDNA) && monDnaGet (look, 0, dummy))
      keySet (ksi.ks, j++) = dummy ;
  array (ksi.x, j - 1, KEY) = 0 ;
  dnaBilan = stackCreate (6 * j) ;
  mydna = arrayCreate (20000, char) ;
  ksnew = cptThisSegments (look, ksi, dnaBilan, mydna, taux, &nochange, TRUE) ;
  if (nochange == 1)
    { arrayDestroy (mydna) ;
      goto abort ;
    }
  monDnaForget (look, dnacontig) ;
  maxdna = arrayMax (mydna) ;
  dnaStoreDestroy (dnacontig, mydna) ;
  newseq = arrayCreate (5 * nochange, BSunit) ;
  i = 0 ;
  shift = 0 ;
  while (nochange)
    { shift2 = pop (dnaBilan, int) ;
      if (shift2 == 1 << 30)
	{ dnakey = pop (dnaBilan, KEY) ;
	  lexReClass (dnakey, &key, _VSequence) ;
	  array (newseq, i++, BSunit).k = key ;
	  shift3 = pop (dnaBilan, int) + shift ;
	  array (newseq, i++, BSunit).i = shift3 ;
	  shift3 = pop (dnaBilan, int) + shift ;
	  array (newseq, i++, BSunit).i = shift3 ;
	  if (!(obj2 = bsCreate (key)) || 
	      !bsGetData (obj2, _Clipping, _Int, &clip))
	    clip = 1 ;
	  bsDestroy (obj2) ;
	  array (newseq, i++, BSunit).i = clip ;
	  clip += arrayMax (monDnaGet (look, contig, dnakey)) - 1 ;
	  array (newseq, i++, BSunit).i = clip ;
	  nochange-- ;
	}
      else
	shift += shift2 ;
    }
  key = arr (newseq, 0, BSunit).k ;
  x = arr (newseq, 1, BSunit).i - 1 ; /* then newpos = oldpos + x */
  if (key != contig)
    messcrash ("It should be the contig") ;
  if (!(obj = bsUpdate (contig)))
    goto abort ;
  if (bsFindKey (obj, _DNA, dnacontig))
    bsAddData (obj, _bsRight, _Int, &maxdna) ;
  if (x)
    { dnaAlignUpdateSubSeqPos (obj, contig, x, _Assembled_from) ;
      dnaAlignUpdateSubSeqPos (obj, contig, x, _Aligned) ;
    } /* first one is the contig */
  bsAddArray (obj, _Assembled_from, newseq, 5) ;
  bsAddTruncatedArray (obj, _Previous_contig, newseq, 5, 3) ;
  bsSave (obj) ;
  reads = keySetReCreate (reads) ;
  if ((x = keySetMax (ksnew.ks)))
    for (i = 0 ; i < x ; i++)
      keySet (reads, i) = keySet (ksnew.ks, i) ;
 abort:
  stackDestroy (dnaBilan) ;
  keySetDestroy (ksnew.ks) ;
  keySetDestroy (ksnew.x) ;
  keySetDestroy (ksi.ks) ;
  keySetDestroy (ksi.x) ;
  arrayDestroy (newseq) ;
}

/***************************************************************/
/*********************** Align Consensus ***********************/
/***************************************************************/

static KEYSET dnaAlignFilterContig (DEFCPT look, KEYSET ks, char *cp)
{ int i, i1, i2, i3 ;
  KEYSET ks1, ks2, ks3, ks4 ;
  KEY key, *keyp ;
  OBJ obj = 0 ;

  if (!ks || !keySetMax (ks))
    { ks1 = keySetCreate () ;
      return ks1 ;
    }
  ks1 = keySetCreate () ; i1 = 0 ;
  ks2 = keySetCreate () ; i2 = 0 ;
  ks3 = keySetCreate () ; i3 = 0 ;
  i = keySetMax (ks) ;
  keyp = arrp (ks, 0, KEY) - 1 ;
  while (keyp++, i--)
    { key = *keyp ;
      if ((class (key) != _VSequence) && !lexReClass (key, &key, _VSequence))
	continue ;
      if ((obj = bsCreate (key)))
	{
	  if (bsFindTag (obj, _Assembled_from)) /* i e contig */
	    keySet (ks1, i1++) = key ;
	  else if (bsFindTag (obj, _Subsequence))    /* i e Assembly */
	    keySet (ks2, i2++) = key ;
	  else
	    keySet (ks3, i3++) = key ;
	  bsDestroy (obj) ;
	}
    }
  keySetSort (ks1) ;
  keySetCompress (ks1) ;
  keySetSort (ks2) ;
  keySetCompress (ks2) ;
  keySetSort (ks3) ;
  keySetCompress (ks3) ;
  ks4 = query (ks2, ">Subsequence") ;
  keySetDestroy (ks2) ;
  ks2 = keySetOR (ks1, ks4) ;	/* contigs */
  keySetDestroy (ks1) ;
  keySetDestroy (ks4) ;
  ks1 = query (ks2, ">Assembled_from") ; /* reads in contigs */
  ks4 = keySetMINUS (ks3, ks1) ; /* new reads */
  keySetDestroy (ks3) ;
  keySetDestroy (ks1) ;
  ks3 = dnaAlignMakeSubSequence (look->link, ks4, cp) ;
  keySetDestroy (ks4) ;
  ks1 = keySetOR (ks2, ks3) ;
  keySetDestroy (ks2) ;
  keySetDestroy (ks3) ;
  return ks1 ;
}

/***************************************************************/

KEY dnaAlignDoMakeSuperLink (KEYSET ks, char *cp)
{ int i ;
  KEY link, dummy ;
  KEYSET ks1 = 0 ;
  char *cq ;
  DEFCPT look = 0 ;

  i = keySetMax (ks) ;
  if (!i)
    return 0 ;
  look = defCptGetLook (1) ;
  if (!cp || !*cp || lexword2key (cp, &dummy, _VSequence))
    cq = dnaAlignNewName (cp) ;
  else
    cq = cp ;
  lexaddkey (cq, &link, _VSequence) ;
  ks1 = dnaAlignFilterContig (look, ks, cq) ;
  defCptDestroyLook (1) ;
  if (!keySetMax (ks1))
    { keySetDestroy (ks1) ;
      return 0 ;
    }
  contigLengthSort (ks1) ;
  alignToolsAdjustLink (link, ks1, 0) ;
  keySetDestroy (ks1) ;
  if (cq != cp) messfree (cq) ;
  return link ;
}

/***************************************************************/

static void dnaAlignMakeSuperLink (DEFCPT look, KEYSET ks)
{ int i ;
  KEY link ;

  i = keySetMax(ks) ;
  if (!i)
    return ;
  lexaddkey (messprintf("_Assembly_%d.%d", look->id, look->tour), &link, _VSequence) ;
  defCptChangeLook (look, look->link, link) ;
  contigLengthSort (ks) ;
  alignToolsAdjustLink (look->link, ks, 0) ;
  keySetSort (ks) ;
#ifndef NON_GRAPHIC
  if (look->display)
    display(link, 0, 0) ;
#endif
}

/***************************************************************/

Array blyDnaGet (KEY link, KEY contig, KEY mykey)
{ static DEFCPT look = 0 ;
  
  if (!look ||  look->magic != 727830152 /* see defcpt.c */
      || look->link != link) 
    look = defCptGetLook (link) ;
  return monDnaGet (look, contig, mykey) ;
}

/***************************************************************/

void dnaAlignFixSegConsensus(DEFCPT look, KEYSET segment)
{ KEY segDna, segSeq, targetDna, targetSeq ;
  KEYSET dna = look->def ;
  int i, j, jstep = 0 ;
  char buff[64] ;
  OBJ obj ;

  look->tour++ ;
  look->assembling = TRUE ;
  i = keySetMax(segment) ;
  if (!keySetExists(ksstep))
    ksstep = 0 ;
  ksstep = keySetReCreate(ksstep) ;
  for (j = 0 ; j < i ; j++)
    { segDna = keySet(segment, j) ;
      lexReClass (segDna, &segSeq, _VSequence) ;
      if ((obj = bsCreate(segSeq)))
	{ if (bsFindTag(obj, _Assembled_from))
	    { sprintf(buff, "_segment_%d.%d.%d", look->id, look->tour, jstep++) ;
	      lexaddkey(buff, &targetDna, _VDNA) ;
	      lexaddkey(buff, &targetSeq, _VSequence) ;
	      doForceAssembleSeg (look, segSeq, targetSeq, 0, -2, 2, 0) ;
	      dnaAlignCleanLeftArrows (0, segSeq, targetSeq, 1, 1) ;
	      segDna = targetDna ;
	      segSeq = targetSeq ;
 	    }
	  bsDestroy(obj) ;
	}
      keySet(dna, j) = segDna ;
      keySet(ksstep, j) = segSeq ;
    }
  if (look->display)
    dnaDispGraph (ksstep, look->tour) ;
  dnaAlignMakeSuperLink (look, ksstep) ;
}

/***************************************************************/

void dnaAlignFixContig (KEY link, KEY key)
{
  OBJ obj = 0 ;
  int i ; char *cp ;
  KEY dnaKey ;
  Array dna = 0 ;
  DEFCPT look = 0 ;

  if (!link)
    look = defCptGetLook (1) ;
  else
    look = defCptGetLook (link) ;
  if (class(key) != _VSequence)
    lexReClass(key, &key, _VSequence) ;
  if ((obj = bsCreate(key)))
    { 
      if (!bsFindTag(obj, _Assembled_from))
	{ bsDestroy(obj) ;
	return ;
	}
      if (!bsGetKey(obj, _DNA, &dnaKey) || !(dna = monDnaGet (look, 0, dnaKey)))
	{ bsDestroy(obj) ;
	  return ;
	}
      bsDestroy(obj) ;
      cp = arrp(dna, 0, char) - 1 ;
      i = arrayMax(dna) ;
      while(cp++, i--) /* kludge pour corriger les n apres edition */
	if (!*cp || *cp == N_)
	  *cp = (char)G_ ;	/* probalement meilleur que A_ */
      doForceAssembleSeg(look, key, key, 0, -2, 2, 1) ;	/* was 0, -2, 0, 1 */
    }
  if (!link)
    defCptDestroyLook (look->link) ; /* look->link should be 1 */
}

/***************************************************************/
#ifdef NON_GRAPHIC
void dnaAlignFixActKeyset(void)
{
}
#else
void dnaAlignFixActKeyset(void) /* fonction a corriger */
{ int i ;
  KEYSET aa ;
  KEY key ;
  void *dummy ;
  OBJ obj = 0 ;
  DEFCPT look = defCptGetLook (0) ;

  if (!keySetActive(&aa, &dummy))
    { messout("First select a keyset containing sequences") ;
      return ;
    }
  messStatus("Fixing KeySet") ;
  for (i = 0 ; i < keySetMax(aa) ; i++)
    { key = keySet(aa, i) ;
      if (class(key) == _VDNA)
	lexReClass(key, &key, _VSequence) ;
      if ((obj = bsCreate(key)))
	{ if (!bsFindTag(obj, _Assembled_from))
	    { bsDestroy(obj) ;
	      continue ;
	    }
	  bsDestroy(obj) ;
	  doForceAssembleSeg (look, key, key, 0, -2, 0, 1) ;
	}
    }
  keySetShow(aa, dummy) ;
}
#endif
/***************************************************************/
/*
static void dnaAlignRemoveSeq(KEY cons, KEY seq, KEYSET new)
{ Array nbt ;
  int i, j ;
  OBJ obj ;


}
*/
/***************************************************************/
/*********************** Other Functions ***********************/
/***************************************************************/

typedef struct  { KEY key ; int ll ;} KEY_LEN ;

static int contigLengthOrder (const void *a, const void *b)
{ int 
    x1 = ((const KEY_LEN*)a)->ll, 
    x2 = ((const KEY_LEN*)b)->ll ; 

   return x2 - x1 ;
}

void contigLengthSort (KEYSET ks)
{ KEY_LEN * kl ;
  Array aa = arrayCreate (50, KEY_LEN) ;
  OBJ obj ;
  int x, i, j = 0 ; KEY key, dummy ;

  i = keySetMax (ks) ;
  while (i--)
    { key = keySet (ks, i) ;
      obj = bsCreate (key) ;
      if (obj)
	{ x = 300 ; /* default */
	  if (bsGetKey (obj, _DNA, &dummy))
	    bsGetData (obj, _bsRight, _Int, &x) ;
	  bsDestroy (obj) ;
	}
      kl = arrayp (aa, j++, KEY_LEN) ;
      kl->key = key ; kl->ll = x ;
    }
  arraySort (aa, contigLengthOrder) ;
  i = keySetMax (ks) ;
  while (i--)
    keySet(ks, i) = arrp (aa, i, KEY_LEN)->key ;
  arrayDestroy (aa) ;
}

/***************************************************************/

KEYSET dnaAlignMakeSubSequence (KEY link, KEYSET reads, char *cp)
{ Array dna, dna2 ;
  char buf[128], buf2[256] ;
  KEYSET ks ;
  OBJ obj = 0 ;
  KEY *kp, dummy, dnaKey, key, key2, dnaKey2 ;
  int i, j = 1, k = 0, max, nn = 1, clip = 1 ;
  DEFCPT look = defCptGetLook (link) ;

  strncpy (buf, cp, 100) ;
  i = keySetMax (reads) ;
  ks = keySetCreate () ;
  if (!i) return ks ;
  kp = arrp (reads, 0, KEY) - 1 ;
  while (kp++, i--)
    { key = *kp ; 
      if (!lexReClass (key, &dnaKey, _VDNA) || !(dna = monDnaGet (look, 0, dnaKey)))
	continue ;
      if ((obj = bsCreate (key)))
	{ if (!bsGetData (obj, _Clipping, _Int, &clip))
	    clip = 1 ;
	  bsDestroy (obj) ;
	}
      else clip = 1 ; /* should not happen */
      dna2 = dnaCopy (dna) ;
      max = arrayMax (dna2) ;
      while (lexword2key(messprintf("%s.%d", buf, j), &dummy, _VSequence)) 
	j++ ;
      sprintf (buf2, "%s.%d", buf, j) ;
      lexaddkey (buf2, &key2, _VSequence) ;
      lexaddkey (buf2, &dnaKey2, _VDNA) ;
      j++ ;
      dnaStoreDestroy (dnaKey2, dna2) ;
      obj = bsUpdate (key2) ;
      bsAddKey (obj, _DNA, dnaKey2) ;
      bsAddData (obj, _bsRight, _Int, &max) ;
      bsAddKey (obj, _Assembled_from, key) ;
      bsAddData (obj, _bsRight, _Int, &nn) ;
      bsAddData (obj, _bsRight, _Int, &max) ;
      bsAddData (obj, _bsRight, _Int, &clip) ;
      clip += max - 1 ;
      bsAddData (obj, _bsRight, _Int, &clip) ;
      bsSave (obj) ;
      keySet (ks, k++) = key2 ;
    }
  return ks ;
}

/***************************************************************/

static void dnaAlignDoSaveAs (DEFCPT look, char *cp)
{ KEY key, dummy, dnaKey, *kp, link ;
  KEYSET ks = 0, ks1 = 0, ks2 = 0 ;
  int i, j;
  char *cq0, buf[256] ;

  if (!look->link)
    { messout ("no consensus to save") ;
      return ;
    }
  strncpy(buf, cp, 128) ;  /* be sure we have room left in name for numbers */
  if (class (look->link) != _VSequence) /* save without fix */
    { 
      if (look->def && keySetMax (look->def))
	{ ks = keySetCreate () ; j = 0 ;
	  for (i = 0 ; i < keySetMax (look->def) ; i++)
	    { key = keySet (look->def, i) ;
	      if (lexReClass (key, &key, _VSequence))
		keySet (ks, j++) = key ;
	    }
	  keySetSort (ks) ;
	  ks1 = query (ks, "NOT Vector") ;
	  keySetDestroy (ks) ;
	  lexaddkey (buf, &link, _VSequence) ;
	  defCptChangeLook (look, look->link, link) ;
	  contigLengthSort (ks1) ;
	  alignToolsAdjustLink (look->link, ks1, 0) ;
	  keySetSort (ks1) ;
	}
      else
	{ messout ("no consensus to save") ;
	  return ;
	}
    }
  else
    { ks1 = queryKey (look->link, ">Subsequence ; NOT Vector") ;
      lexAlias (&(look->link), messprintf("%s", buf), FALSE, FALSE) ;
    }
  ks = query (ks1, "Assembled_from AND IS _seg*") ; /* a tester absolument */
  ks2 = keySetMINUS (ks1, ks) ;
  keySetDestroy (ks1) ;

  cq0 = buf + strlen(buf) ; *cq0++ = '.' ; *cq0 = 0 ;
  j = 1 ;
  contigLengthSort (ks) ;
  i = keySetMax (ks) ;
  kp = arrp (ks, 0, KEY) - 1 ;
  while (kp++, i--)
    { key = *kp ; dnaKey = 0 ; lexReClass (key, &dnaKey, _VDNA) ;
      while (lexword2key(messprintf("%s%d", buf, j), &dummy, _VSequence)) j++ ; 
      lexAlias (&key, messprintf("%s%d", buf, j), FALSE, FALSE) ;
      if (dnaKey) lexAlias (&dnaKey, messprintf("%s%d", buf, j), FALSE, FALSE) ;
    }
  if (keySetMax (ks2))
    { ks1 = dnaAlignMakeSubSequence (look->link, ks2, cp) ;
      keySetDestroy (ks2) ;
      keySetSort (ks) ;
      keySetSort (ks1) ;
      ks2 = keySetOR (ks, ks1) ;
      alignToolsAdjustLink (look->link, ks2, 0) ;
    }
  look->assembling = FALSE ;
  keySetDestroy (ks) ;
  keySetDestroy (ks1) ;
  keySetDestroy (ks2) ;
  if (look->display)
    { dnaDispGraph(ksstep, look->tour) ;
#ifndef NON_GRAPHIC
      if (display (look->link, 0, 0))
	fMapReDrawWindow () ;
#endif
    }
}

/***************************************************************/

BOOL dnaAlignCopyContig (KEY key, KEY *kp, char *cp, BOOL left)
{ 
  int j = 0, k, max, x1, x2, clip ;
  KEY seqKey, dnaKey, nKey, ndnaKey, key2, dummy ;
  BSunit *u ;
  Array seq = 0, dna ;
  OBJ obj = 0, obj2 = 0 ;
  char buff[512] ;

  if (strlen (cp) > 250)
    { messerror ("Char cp too long in dnaAlignCopyContig") ;
      return FALSE ;
    }
  if (class (key) == _VSequence)
    seqKey = key ;
  else if (class (key) == _VDNA)
    lexReClass (key, &seqKey, _VSequence) ;
  else return FALSE ;
  obj = bsCreate (seqKey) ;
  if (!obj || !bsGetKey (obj, _DNA, &dnaKey) ||
      !bsGetData (obj, _bsRight, _Int, &max))
    { bsDestroy (obj) ;
      return FALSE ;
    }
  seq = arrayCreate (200, BSunit) ;
  if (!bsFindTag (obj, _Assembled_from) || !bsFlatten (obj, 5, seq))
    { bsDestroy (obj) ;
      arrayDestroy (seq) ;
      return FALSE ;
    }
  dna = dnaGet (dnaKey) ;
  while (lexword2key (messprintf ("%s.%d", cp, j), &dummy, _VSequence))
    j++ ;
  sprintf (buff, "%s.%d", cp, j) ;
  lexaddkey (buff, &nKey, _VSequence) ;
  lexaddkey (buff, &ndnaKey, _VDNA) ;
  dnaStoreDestroy (ndnaKey, dna) ;
  obj2 = bsUpdate (nKey) ;
  bsAddKey (obj2, _DNA, ndnaKey) ;
  bsAddData (obj2, _bsRight, _Int, &max) ;
  for (k = 0 ; k < arrayMax (seq) ; k += 5)
    { u = arrp (seq, k, BSunit) ;
      key2 = u[0].k ;
      x1 = u[1].i ;
      x2 = u[2].i ;
      bsAddKey (obj2, _Assembled_from, key2) ;
      bsAddData (obj2, _bsRight, _Int, &x1) ;
      bsAddData (obj2, _bsRight, _Int, &x2) ;
      clip = u[3].i ;
      if (clip)
	{ bsAddData (obj2, _bsRight, _Int, &clip) ;
	  clip = u[4].i ;
	  bsAddData (obj2, _bsRight, _Int, &clip) ;
	}
    }
  if (left)
    {
      if (bsGetArray (obj, _Aligned, seq, 5))
	bsAddArray (obj2, _Aligned, seq, 5) ;
      if (bsGetArray (obj, _Previous_contig, seq, 3))
	bsAddArray (obj2, _Previous_contig, seq, 3) ;
    }
  bsSave (obj2) ;
  bsDestroy (obj) ;
  arrayDestroy (seq) ;
  *kp = nKey ;
  return TRUE ;
}

/***************************************************************/

static void dnaAlignCopy (DEFCPT look, char *cp, BOOL order)
{ OBJ Link = 0, obj1 = 0, obj2 = 0 ;
  Array seq = 0, dna = 0, unit = 0, bilan = 0 ;
  KEYSET ks ;
  mytime_t tim ;
  KEY dummy, key, key2, key3, dnakey1, dnakey2 ;
  int i, j, k, max, x1, x2, ii = 0, jj  = 0, clip ;
  char buf[255], buff[512] ;
  BSunit *u, *v ;

  strncpy (buf, cp, 128) ;
  unit = arrayCreate (100, BSunit) ;
  Link = bsCreate (look->link) ;
  if (!Link || ! bsFindTag (Link, _Subsequence) || !bsFlatten (Link, 3, unit))
    { bsDestroy (Link) ;
      goto abort ;
    }
  bsDestroy (Link) ;
  v = arrp (unit, 0, BSunit) ;
  bilan = arrayCreate (100, BSunit) ;
  ks = keySetCreate () ;
  j = 1 ;
  lexaddkey (buf, &key, _VSequence) ;
  tim = timeParse ("today") ;
  Link = bsUpdate (key) ;
  bsAddKey (Link, _Derived_from, look->link) ;
  bsAddData (Link, _bsRight, _DateType, &tim) ;
  bsSave (Link) ;
  defCptChangeLook (look, look->link, key) ; /* modify look in associator */
  for (i = 0 ; i < arrayMax (unit) ; v += 3, i += 3)
    { obj1 = bsCreate (v->k) ;
      if (!obj1 || !bsGetKey (obj1, _DNA, &dnakey1) ||
	  !bsGetData (obj1, _bsRight, _Int, &max))
	{ bsDestroy (obj1) ;
	  continue ;
	}
      seq = arrayReCreate (seq, 200, BSunit) ;
      if (!bsFindTag (obj1, _Assembled_from) || !bsFlatten (obj1, 5, seq))
	{ bsDestroy (obj1) ; /* a revoir ? */
	  continue ;
	}
      dna = dnaGet (dnakey1) ;
      while (lexword2key (messprintf ("%s.%d", buf, j), &dummy, _VSequence))
	j++ ;
      sprintf (buff, "%s.%d", buf, j) ;
      lexaddkey (buff, &key2, _VSequence) ;
      lexaddkey (buff, &dnakey2, _VDNA) ;
      j++ ;
      dnaStoreDestroy (dnakey2, dna) ;
      obj2 = bsUpdate (key2) ;
      bsAddKey (obj2, _DNA, dnakey2) ;
      bsAddData (obj2, _bsRight, _Int, &max) ;
      for (k = 0 ; k < arrayMax (seq) ; k += 5)
	{ u = arrp (seq, k, BSunit) ;
	  key3 = u[0].k ;
	  x1 = u[1].i ;
	  x2 = u[2].i ;
	  bsAddKey (obj2, _Assembled_from, key3) ;
	  bsAddData (obj2, _bsRight, _Int, &x1) ;
	  bsAddData (obj2, _bsRight, _Int, &x2) ;
	  clip = u[3].i ;
	  if (clip)
	    { bsAddData (obj2, _bsRight, _Int, &clip) ;
	      clip = u[4].i ;
	      bsAddData (obj2, _bsRight, _Int, &clip) ;
	    }
	}
      if (bsGetArray (obj1, _Aligned, seq, 5))
	bsAddArray (obj2, _Aligned, seq, 5) ;
      if (bsGetArray (obj1, _Previous_contig, seq, 3))
	bsAddArray (obj2, _Previous_contig, seq, 3) ;
      bsSave (obj2) ;
      bsDestroy (obj1) ;
      keySet (ks, ii++) = key2 ;
      array (bilan, jj++, BSunit).k = key2 ;
      array (bilan, jj++, BSunit).i = (v+2)->i - (v+1)->i ;
    }
  if (order)
    alignToolsAdjustLink (look->link, 0, bilan) ;
  else
    { contigLengthSort (ks) ;
      alignToolsAdjustLink (look->link, ks, 0) ;
    }
  arrayDestroy (bilan) ;
  keySetDestroy (ks) ;
 abort:
  arrayDestroy (unit) ;
  arrayDestroy (seq) ;
}  

/***************************************************************/

char *dnaAlignNewName (char *cp)
{ char *ctmp, *cnew, *cq ;
  KEYSET ks = 0 ;
  int i = 1 ;
  KEY dummy ;

  ctmp = strnew (cp ? cp : "a", 0) ;
  cq = ctmp + strlen(ctmp) - 1 ;
  while (cq > ctmp && (isdigit((int)*cq) || ispunct ((int)*cq))) *cq-- = 0 ;
  while (TRUE)
    { while (lexword2key (messprintf ("%s%d", ctmp, i), &dummy, _VSequence) ||
	     lexword2key (messprintf ("%s%d.", ctmp, i), &dummy, _VSequence))
	i++ ;
      ks = query (0, messprintf(">? Sequence IS %s%d OR IS %s%d.*", ctmp, i, ctmp, i)) ;
      if (!keySetMax (ks))
	break ;
      keySetDestroy (ks) ;
      i++ ;
    }
  keySetDestroy (ks) ;
  cnew = strnew (messprintf ("%s%d", ctmp, i), 0) ;
  messfree (ctmp) ;
  return cnew ;
}

/***************************************************************/

void dnaAlignSave (DEFCPT look, char *cp, BOOL order)
{ KEYSET ks = 0 ;
  char *cnew ;

  if (!sessionGainWriteAccess ())
    { 
      messout ("Sorry, you cannot gain write access.") ;
      return ;
    }
  ks = query (0, messprintf(">? Sequence IS %s OR IS %s.*", cp, cp)) ;
  if (keySetMax (ks))
    cnew = dnaAlignNewName (cp) ;
  else
    cnew = cp ;
  keySetDestroy (ks) ;
  if (look->assembling)
    dnaAlignDoSaveAs (look, cnew) ;
  else
    dnaAlignCopy (look, cnew, order) ;
  if (cnew != cp) messfree (cnew) ;
  sessionDoSave (TRUE) ;
}

/***************************************************************/

void dnaAlignSaveAs (DEFCPT look, char **cpp)
{ char *cp ;

  if (!messPrompt("Save Your Assembly as", "contig", "w"))
    return ;
  cp = freeword() ;
  *cpp = messalloc (strlen(cp) + 1) ;
  strcpy (*cpp, cp) ;
  dnaAlignSave (look, *cpp, FALSE) ;
}

/***************************************************************/

char *dnaAlignSaveDefault (DEFCPT look)
{ char *cp = dnaAlignNewName ("a") ;
  dnaAlignSave (look, cp, FALSE) ;
  return cp ;
}

/***************************************************************/

Array dnaAlignCptErreur(Array longDna, Array shortDna,
			int *debutp, int *finp, int *shdebutp, int *shfinp)
{
/*  if (*debutp < 0 || *shdebutp < 0 || *shfinp < 0)
    invokeDebugger () ; */
  return alignToolsCptErreur(longDna, shortDna, debutp, finp, shdebutp, shfinp) ;
}

/***************************************************************/
/* prend des cles dna ou sequence et retourne une cle de la meme classe
   que la cle 1. Attention les deux cles doivent exister */
/* assemble si elle trouve une zone de 20 bases identiques dans les deux
   sequences. Retourne le taux d'erreur trouve. */
/* Attention ne doit etre utilisee que pour assembler des contigs pas de
sequences isolees */
static KEY dnaAlignPaires (DEFCPT look, KEY key1, KEY key2, int *taux)
{ int i, j, max1, max2, x = 0, dx = 15, psta, pend, sens, x1, x2, shift ;
  int y1, y2, n, class ;
  char *cp, *cq, buff[12] ;
  static int num = 0 ;
  Array dna1, dna2, newdna, erreur = 0 ;
  KEY key, keydna, key1dna, key2dna ;
  OBJ obj ;

  class = class(key1) ;
  if (class(key1) == _VDNA)
    { key1dna = key1 ;
      if (!lexReClass(key1, &key1, _VSequence))
	return 0 ; /* messerror ("key sequence does not exist") ; */
    }
  else if (!lexReClass(key1, &key1dna, _VDNA))
    return 0 ; /* messerror ("Key dna does not exist") ; */
  if (class(key2) == _VDNA)
    { key2dna = key2 ;
      if (!lexReClass(key2, &key2, _VSequence))
	return 0 ; /* messerror ("key sequence does not exist") ; */
    }
  else if (!lexReClass(key2, &key2dna, _VDNA))
    return 0 ; /* messerror ("Key dna does not exist") ; */
  dna1 = monDnaGet (look, 0, key1dna) ;
  dna2 = monDnaGet (look, 0, key2dna) ;
  if (!dna1 || !dna2 || !(max1 = arrayMax(dna1)) || !(max2 = arrayMax(dna2)))
    return 0 ; /* messerror ("Wrong keys passed") ; */
  if (max2 > max1) /* I want the longest dna is the first */
    { key = key2 ; key2 = key1 ; key1 = key ;
      newdna = dna2 ; dna2 = dna1 ; dna1 = newdna ;
      i = max2 ; max2 = max1 ; max1 = i ;
    }
  for (x = 0 ; x < max2 - dx ; x += 10)
    if ((erreur = alignToolsMatch(dna1, 0, max1 - 1, 0, max1 - 1,
				 dna2, x, x + dx, x, x + dx,
				 &psta, &pend, &sens, 1, 0, 1, FALSE)))
      break ;
  if (!erreur)
    return 0 ;
  x1 = x2 = 0 ; y1 = max1 - 1 ; y2 = max2 - 1 ;
  if (sens == -1)
    x += dx ;
  localCptErreur(dna1, x1, y1, psta, dna2, x2, y2, x, sens, &n, &x2, &y2, &i, erreur) ;
  *taux = (arrayMax(erreur) - n) * 100 / i ;
  newdna = dnaCopy(dna1) ;
  dnaAddPiece(newdna, dna2, erreur, sens, &x2, &y2, &shift) ;
  arrayDestroy (erreur) ;
  x1 += shift ; y1 += shift ; /* x1, y1 position de dna1 dans bilan */
  i = x1 < x2 ? x1 : x2 ; /* debut du dna */
  j = y1 > y2 ? y1 : y2 ; /* fin du dna */
  x1 -= i - 1 ; /* elimination des zeros et passage coordonne bio */
  x2 -= i - 1 ;
  y1 -= i - 1 ;
  y2 -= i - 1 ;
  n = j - i + 1 ;
  arrayMax(newdna) = n ;
  if (i)
    { cp = arrp(newdna, 0, char) ;
      cq = arrp(newdna, i, char) ;
      while(n--)
	*cp++ = *cq++ ;
    }
  n = arrayMax(newdna) ;
  if (sens == -1)
    { x = x2 ;
      x2 = y2 ;
      y2 = x ;
    }
  sprintf(buff, "_paire%d", num++) ;
  lexaddkey(buff, &key, _VSequence) ;
  lexaddkey(buff, &keydna, _VDNA) ;
  if ((obj = bsUpdate(key)))
    { bsAddKey(obj, _DNA, keydna) ;
      bsAddData(obj, _bsRight, _Int, &n) ;
      if (bsFindTag(obj, _Assembled_from))
	bsRemove(obj) ;
      if (bsFindTag(obj, _Assembled_into))
	bsRemove(obj) ;
      bsAddKey(obj, _Assembled_from, key1) ;
      bsAddData(obj, _bsRight, _Int, &x1) ;
      bsAddData(obj, _bsRight, _Int, &y1) ;
/*    x = 1 ;  clips inutile pour les contigs
      bsAddData (obj, _bsRight, _Int, &x) ;
      x += arrayMax (dna1) - 1 ;
      bsAddData (obj, _bsRight, _Int, &x) ; */
      bsAddKey(obj, _Assembled_from, key2) ;
      bsAddData(obj, _bsRight, _Int, &x2) ;
      bsAddData(obj, _bsRight, _Int, &y2) ;
/*    x = 1 ;
      bsAddData (obj, _bsRight, _Int, &x) ;
      x += arrayMax (dna2) ;
      bsAddData (obj, _bsRight, _Int, &x) ; */
      bsSave(obj) ;
    }
  else
    messcrash("I can't Create Obj") ;
  dnaStoreDestroy(keydna, newdna) ;
  if (class == _VSequence)
    return key ;
  else
    return keydna ;
}

/***************************************************************/

Array dnaAlign(Array longDna, Array shortDna, int *debutp, int *finp, 
	       int *shdebutp, int *shfinp)
{ int x1, nsho, nlon, dn = 40, pols, pole, pos ;
  int sens, nerrmax, taux ;
  Array erreur = 0 ;

  nsho = arrayMax(shortDna) ;
  nlon = arrayMax(longDna) ;
  if (!nsho)
    return 0 ;  
  *shdebutp = 0 ;
  *debutp = 0 ;
  *shfinp = nsho - 1 ;
  *finp = nlon - 1 ;
  if (dn > nsho - 1) dn = nsho - 1 ;
  for (nerrmax = 5 ; nerrmax < 20 ; nerrmax += 3)
    { taux = (100 * nerrmax) / dn ;
      for (x1 = (nsho > 90 ? 40 : 0) ; x1 < nsho - dn - 1 ; x1 += 20)
	if ((erreur = alignToolsMatch(longDna, 0, nlon - 1, 0, nlon - 1,
				     shortDna, x1, x1 + dn, x1, x1 + dn, 
				     &pols, &pole, &sens, 1, taux, 1, 
				     FALSE)))
	  break ;
      if (erreur)
	{ arrayDestroy(erreur) ;
	  if (sens == 1)
	    pos = x1 ;
	  else
	    pos = x1 + dn ;
	  erreur = alignToolsMakeErreur(longDna, debutp, finp, pols, 
					shortDna, 0, nsho - 1, pos, sens) ;
	  break ;
	}
    }
  return erreur ;
}

/***************************************************************/

 /* renvoie dans topp, endp, la position dans long des bases x1 x2 de court */
BOOL dnaAlignForceMatch (Array shortDna, int x1, int x2, 
			 Array longDna, int y1, int y2, 
			 int *topp, int *endp, int *sensp)
{ BOOL result = FALSE ;
  Array err = 0 ;

  if ((err = alignToolsMatch(longDna, y1, y2, y1, y2, shortDna, x1, x2, x1, x2,
			    topp, endp, sensp, 0, 10, 8, FALSE)))
    { arrayDestroy(err) ;
      result = TRUE ;
    }
  return result ;
}

/***************************************************************/

static int tripletOrder (const void *a, const void *b)
{ if (((const PAIR*)a)->u1 != ((const PAIR*)b)->u1)
    return
      ((const PAIR*)a)->u1 - ((const PAIR*)b)->u1 ;
  else
    return
      ((const PAIR*)a)->u2 - ((const PAIR*)b)->u2 ;
}

/***************************************************************/

Array dnaAlignMatchTriplet (Array dna1, int x1, int y1, int *n1, Array dna2,
			    int x2, int y2, int *n2, BOOL direct)
{ Array result = arrayCreate(100, PAIR), inter1 = 0, inter2 = 0 ;
  unsigned int i, j, k, oli, ilo, win = 0x0000003f ; /* trimere */
  int n, m, *np, *mp, nn = 0;

  *n1 = *n2 = 0 ;
  i = 4 ;
  while(i--)
    { j = 4 ;
      while(j--)
	{ if (j == i)
	    continue ;
	  k = 4 ;
	  while(k--)
	    { if (k == i || k == j)
		continue ;
	      oli = ((i << 4) | (j << 2) | k) ;
	      if ((inter1 = alignToolsFindShortOligo (dna1, x1, y1, win, 3, 
						     oli)))
		*n1 += arrayMax(inter1) ;
	      if (direct)
		ilo = oli ;
	      else
		ilo = (((~k & 3) << 4) | ((~j & 3) << 2) | (~i & 3)) ;
	      if ((inter2 = alignToolsFindShortOligo (dna2, x2, y2, win, 3,
						     ilo)))
		*n2 += arrayMax(inter2) ;
	      if (!inter1 || !inter2)
		{ arrayDestroy(inter1) ;
		  arrayDestroy(inter2) ;
		  continue ;
		}
	      n = arrayMax(inter1) ;
	      np = arrp(inter1, 0, int) - 1 ;
	      m = arrayMax(inter2) ;
	      mp = arrp(inter2, 0, int) - 1 ;
	      while(mp++, m--)
		*mp += 1 ;
	      while(np++, n--)
		{ *np += 1 ;
		  m = arrayMax(inter2) ;
		  mp = arrp(inter2, 0, int) ;
		  while(m--)
		    { arrayp(result, nn, PAIR)->u1 = *np ;
		      arrp(result, nn++, PAIR)->u2 = *mp++ ;
		    }
		}
	      arrayDestroy(inter1) ;
	      arrayDestroy(inter2) ;
	    }
	}
    }
  if (arrayMax(result) > 1)
    arraySort(result, tripletOrder) ;
  return result ;
}

/***************************************************************/
/* Attention on recoit les coordonnees Jean (1ere incluse derniere excluse */
void dnaAlignRecale (Array longDna, int *xl, int *yl, Array shortDna,
		     int xs, int ys)
{ Array err = 0, dna ;
  int deb, fin, sens = 1,
      a1, a2, amax, taux,  max = arrayMax(longDna) ;
  BOOL isUp ;

  a1 = *xl ; a2 = *yl ;
  amax = arrayMax (longDna) ;
  isUp = a1 < a2 ? FALSE : TRUE ;
  if (isUp)
    { 
      dna = dnaCopy (longDna) ;
      reverseComplement (dna) ;
      a1 = amax - a1 - 1 ;
      a2 = amax - a2 - 1 ;
    }
  else
    dna = longDna ;

  deb = a1 > 50 ? a1 - 50 : 0 ;
  fin = a2 < max - 50 ? a2 + 49 : max - 1 ; /* passage last inclus */

  taux = 2 ;  
  while (taux < 65 && 
	 !(err = alignToolsMatch (dna, deb, fin, deb, fin, shortDna, xs, ys,
				  xs, ys, &deb, &fin, &sens, 1, taux, 5, 
				  TRUE)))
    taux *= 2 ;
  
  if (!err)
    goto done ;
  if (sens !=  1)
    { /* printf ("I have complemented the sequence\n") ; */
      goto done ;
    }
  
  a1 = deb ;
  a2 = fin + 1 ; /* retour coordonnees Jean */
  
  if (isUp)
    { a1 = amax - a1 - 1 ;
      a2 = amax - a2 - 1 ;
    }
  *xl = a1 ; *yl = a2 ;

 done:
  arrayDestroy(err) ;
  if (isUp)
    arrayDestroy (dna) ;
}
  
/***************************************************************/
typedef struct { int u1, u2, g ; } FIT ;
typedef struct { int x1, x2, x3, a1, a2, a3, ct, ce ; KEY  key, dnaKey ; } READ ;

static int fitOrder (const void *a, const void *b) 
{ int 
    x1 = ((const FIT *)a)->u1 , y1 = ((const FIT *)b)->u1 ,
    x2 = ((const FIT *)a)->u2 , y2 = ((const FIT *)b)->u2 ;
  return x1 < y1 ? -1 : (x1 == y1 ? x2 - y2 :  1 ) ;
}

/***************************************************************/
/*
static Array dnaAlignMatchTriplet (Array dna1, int x1, int x2,
			     Array dna2, int y1, int y2, BOOL direct)
{
  Array aa = arrayCreate(10, int) ;
  int i = 0, n, err, delta ;

  n = randint() % 20 ;
  for (i = 0 ; i < 2*n ; i++)
    { delta = randint() % 15 ;
      delta = randint() % 8 ;
      array(aa, i, int) = delta ;
      array(aa, i+ 1 , int) = delta + err ;
    }
  return aa ;
}
*/
/***************************************************************/

int fMapTraceFindBestGroup (Array fit, BOOL direct, int wobble, int taille)
{ int 
    ddx2 = wobble * wobble + 1, dx1, dx2, dx,
    bestg, found, i, j, v1, g, gmax, ng, nbestg ;
  FIT *f, *f1 ;
  
  /* now make groups */
  arraySort (fit, fitOrder) ;

  found = 1 ; g = 0 ;
  while (found)
    { found = 0 ;
      for (i = 0 ; i < arrayMax(fit) ; i++)
	{ f = arrp(fit, i, FIT) ;
	  if (f->g)
	    continue ;
	  found = 1 ;
	  g++ ;
	  f->g = g ;
	  v1 = -1 ;
	  for (j = i + 1 ; j < arrayMax(fit) ; j++)
	    { f1 = arrp(fit, j, FIT) ;
	      if (f1->g || f1->u1 == v1)
		continue ;
	      v1 = f1->u1 ;
	      dx1 = f1->u1 - f->u1 ; 
	      dx2 = f1->u2 - f->u2 ;
	      if (direct)
		{ dx = dx1 - dx2 ;
		  if (dx1 >= taille && dx2 >= taille  && dx*dx < ddx2) 
		    { f1->g = g ;
		      f = f1 ;
		    }
		}
	      else
		messcrash("Jean doit une bierre a Ulrich") ;
	    }
	}
    }

  /* find best group */
  gmax = g + 1 ;
  bestg = 0 ; nbestg = 0 ;
  for (g = 0 ; g < gmax ; g++)
    { ng = 0 ;
      for (i = 0, f = arrp (fit, 0, FIT)  ; i < arrayMax(fit) ; i++, f++)
	if (f->g == g) ng++ ;
      if (ng > nbestg)
	{ nbestg = ng ; bestg = g ; }
    }
/*
  if (nbestg > 5)
    for (i = 0 , f = arrp (fit, 0, FIT)  ; i < arrayMax(fit) ; i++, f++)
      { printf( "\n group %d, pair %d  %d ", f->g, f->u1, f->u2) ;
	if (f->g == bestg)
	  { n++ ; }
      }
*/
  	
  return bestg ;
}

/***************************************************************/

int contigFindMatch (KEY contig1, KEY contig2, Array dnaContig1, Array dnaContig2,
		     Array r1, Array r2, BOOL direct, Array *fitp)
{ int bestg, ntr1, ntr2, n, i, u1, u2, i1, i2 ;
  Array dna1, dna2, aa = 0, fit = arrayCreate (100, FIT) ;
  READ *fp1, *fp2 ; PAIR *pp ;
  FIT *f ;
  
  i1 = arrayMax(r1) ;
  fp1 = arrp(r1, 0, READ) -1 ;
  while (fp1++, i1--)
    { i2 = arrayMax(r2) ;
      fp2 = arrp(r2, 0, READ) -1 ;
      dna1 = dnaGet (fp1->dnaKey) ;
      fp1->x3 = arrayMax(dna1) ;
      
      while (fp2++, i2--)
	{
	  dna2 = dnaGet (fp2->dnaKey) ;
	  fp2->x3 = arrayMax(dna2) ;
	  aa = dnaAlignMatchTriplet (dna1, fp1->x1, fp1->x3, &ntr1, 
				     dna2, fp2->x1, fp2->x3, &ntr2, direct) ;
	  arrayDestroy(dna2) ;
	  i = arrayMax(aa) ; pp = arrp(aa, 0, PAIR) - 1 ;
	  while(pp++, i--)
	    { u1 = pp->u1 ;
	      u2 = pp->u2 ;
	      f = arrayp(fit, arrayMax(fit), FIT) ;
	      if (fp1->a1 < fp1->a2)
		f->u1 = u1 + fp1->a1 - fp1->x1 ;
	      else
		f->u1 = - u1 + fp1->a1 + fp1->x1 ;
	      if (fp2->a1 < fp2->a2)
		f->u2 = u2 + fp2->a1 - fp2->x1 ;
	      else
		f->u2 = - u2 + fp2->a1 + fp2->x1 ;
	      f->g = 0 ;
	    }
	}
      arrayDestroy(dna1) ;
    }
  arrayDestroy (aa) ;

  bestg = fMapTraceFindBestGroup (fit, direct, 2,3) ;
  i = arrayMax(fit) ; n = 0 ;
  f = arrp(fit, 0, FIT) - 1 ;
  while (f++, i--)
    if (f->g == bestg) 
      { n++ ;
	
      }
  return n ;
}

/***************************************************************/
   /* if zone > 0, assemble in ends, else anywhere */
BOOL dnaAlignAsmbPaireDna (Array dna1, Array dna2, int taille, int max, int zone, 
			   int nn0, int *v1, int *v2, int *sens, BOOL tryboth)
{ int i, n, nn, nnt = 0, u1, u2, gg, tour = 2 ;
  Array mm = 0, fit = 0, dna2r = dna2, dna2rr = 0 ;
  BOOL ok = FALSE ;
  FIT *f ;

  while (tour--)
    { mm = alignToolsMakeShortMatch (dna1, dna2r, 1, taille, max, zone) ;
/*    printf("MakeShort got %d %d-mers", mm ? arrayMax(mm)/ 2 : 0, taille) ; */
      if (mm && arrayMax(mm))
	{ gg = dnaAlignMakeGroups (mm, &fit, TRUE, taille) ; 
	  nn = 0 ;
	  i = arrayMax(fit) ;
	  f = arrp(fit, 0, FIT) - 1 ;
	  while (f++, i--)
	    { if (f->g == gg)
		nn++ ;
	    }
/*        printf("%d en direct\n", nn) ; */
	  if (nn > nn0 && nn > nnt)
	    { n = 0 ;
	      i = arrayMax(fit) ;
	      f = arrp(fit, 0, FIT) - 1 ;
	      while (f++, i--)
		{ if (f->g == gg)
		    n++ ;
		}
	      n /= 2 ;
	      i = arrayMax(fit) ;
	      f = arrp(fit, 0, FIT) - 1 ;
	      u1 = u2 = -1 ;
	      while (f++, i--)
		if (f->g == gg)
		  { u1 = f->u1 ;
		    u2 = f->u2 ;
		    if (!n--)
		      break ;
		  }
	      if ( zone <= 0 ||
		   (u1 <= zone && u2 + zone >= arrayMax (dna2)) ||
		    (u2 <= zone && u1 + zone >= arrayMax (dna1)) )		   
		{ nnt = nn ;
		  if (tour)
		    { *sens = 1 ; /* reverse */
		      *v1 = u1 ;
		      *v2 = u2 ;
		    }
		  else
		    { *sens = -1 ; /* reverse */
		      *v1 = u1 ;
		      *v2 = arrayMax(dna2) ; /* 2 lines, beware of unsigne d int */
		      *v2 -= u2 + 1 ;
		    }
		  ok = TRUE ;
		}
	    }
	  arrayDestroy (fit) ;
	  arrayDestroy (mm) ;
	  if (ok && !tryboth)
	    break ;
	  if (!dna2rr)
	    { 
	      dna2r = dna2rr = dnaCopy (dna2) ;
	      reverseComplement (dna2r) ;
	    }
	}
    }
  arrayDestroy (dna2rr) ;
  return ok ;
}

/***************************************************************/
/* obj ouvert en ecriture et non sauve ici */
static void dnaAlignAddPrevious (OBJ obj, KEY key, int v1, int v2)
{ int i, sens, x1, x2 ;
  KEY key1 ;
  OBJ petit = 0 ;
  Array prev = 0 ;

  petit = bsCreate (key) ;
  prev = arrayCreate (10, BSunit) ;
  if (!petit || !bsFindTag (petit, _Previous_contig) || !bsFlatten (petit, 3, prev))
    goto abort ;
  sens = v1 < v2 ? 1 : -1 ;
  for (i = 0 ; i < arrayMax (prev) ; i += 3)
    { key1 = arrp (prev, i, BSunit)->k ;
      x1 = arrp (prev, i+1, BSunit)->i ;
      x2 = arrp (prev, i+2, BSunit)->i ;
      if (bsAddKey (obj, _Previous_contig, key1))
	{
	  if (sens == 1)
	    { x1 += v1 - 1 ; x2 += v1 - 1 ; }
	  else
	    { x1 = v1 - x1 + 1 ; x2 = v1 - x2 + 1 ; }
	  bsAddData (obj, _bsRight, _Int, &x1) ;
	  bsAddData (obj, _bsRight, _Int, &x2) ;
	}
    }
 abort:
  arrayDestroy (prev) ;
  bsDestroy (petit) ;
}
/***************************************************************/
/* pre & 1 => previous de key1 ; pre & 2 => previous de key2 */
KEY  dnaAlignAsmbPaire21 (KEY key1, KEY key2, int id, int isTour, int taille,
			  int pre, char *nm)
{ int i, nn0, uu, u1, u2, sens = 1, x, v1, v2, w1, w2, p1, p2 ;
  KEY key, seqKey, dnaKey ;
  OBJ obj, Contig = 0 ;
  BSunit *u ;
  char *cp, *cq ;
  Array units = 0, dna1, dna2, dna ;
  char nom[48] ;
  BOOL isBout = FALSE ;

  if (taille == -1)
    isBout = TRUE ;
  if (taille < 3) taille = 3 ;
  switch (taille)
    { 
    case 7:
      nn0 = 5 ;
      break ;
    case 5:
    default:
      nn0 = 35/taille ;
      break ;
    }

  if ((class(key1) != _VDNA && class(key1) != _VSequence) ||
      (class(key2) != _VDNA && class(key2) != _VSequence))
    { messout("Please, choose a paire of sequences") ;
      return 0 ;
    }
  lexReClass(key1, &key1, _VDNA) ;
  lexReClass(key2, &key2, _VDNA) ;

  dna1 = dnaGet(key1) ; if (!dna1) return 0 ;
  dna2 = dnaGet(key2) ; if (!dna2) { arrayDestroy (dna1); return 0 ; }

  if (isTour)
    sprintf(nom, "_segment_%d.%d", id, isTour) ;
  else if (nm && *nm)
    strncpy (nom, nm, 48) ;
  else
    { strncpy (nom, name(key1), 48) ;
      cp = nom + strlen (nom) - 1 ;
      while (cp > nom && *cp == '.') *cp-- = 0 ;
    }
  if (dnaAlignAsmbPaireDna (dna1, dna2, taille, 2000, 500, nn0, &u1, &u2, &sens, FALSE))
    goto ok ;
  if (isBout &&  baseCallTrackBout (dna1, dna2,  &u1, &u2, &sens))
    goto ok ;

  arrayDestroy (dna1) ;
  arrayDestroy (dna2) ;
  return 0 ;
 ok:
  i = 0 ;
  while (++i)
    { if (!lexword2key(messprintf("%s.%d", nom, i),
		       &key, _VSequence))
	break ;
    }
  lexaddkey(messprintf("%s.%d", nom, i), &dnaKey, _VDNA) ;
  lexaddkey(messprintf("%s.%d", nom, i), &seqKey, _VSequence) ;

  if (sens == -1)
    { 
      reverseComplement (dna2) ;
      u2 = arrayMax (dna2) - u2 - 1 ;
    }
  i = arrayMax(dna1) + arrayMax(dna2) ;
  dna = arrayCreate (i, char) ;
  array(dna, i - 1, char) = 0 ;
/* copy the longest sequence left of the cut */
  cp = arrp (dna, 0, char)  ;
  if (u1 > u2)
    { cq = arrp(dna1, 0, char) ;
      i = uu = u1 ;
      while (i--) *cp++ = *cq++ ;
    }
  else
    { cq = arrp(dna2, 0, char) ;
      i = uu = u2 ;
      while (i--) *cp++ = *cq++ ;
    }
/* now the right longest one */
  if (arrayMax(dna1) - u1 > arrayMax(dna2) - u2)
    { cq = arrp(dna1, u1, char) ;
      i = arrayMax(dna1) - u1 ;
      arrayMax(dna) = uu + i ;
      while (i--) *cp++ = *cq++ ;
    }
  else
    { cq = arrp(dna2, u2, char) ;
      i = arrayMax(dna2) - u2 ;
      arrayMax(dna) = uu + i ;
      while (i--) *cp++ = *cq++ ;
    }
  v1 = uu - u1 + 1 ; w1 = uu - u2 + 1 ;/* no zero */
  v2 = v1 + arrayMax(dna1) - 1 ;
  w2 = w1 + arrayMax(dna2) - 1 ;
  if (sens == -1)
    { i = w1 ; w1 = w2 ; w2 = i ;
    }
  i = arrayMax(dna) ;
  dnaStoreDestroy(dnaKey, dna) ;

  if (i && (obj = bsUpdate(seqKey)))
    { bsKill (obj) ;
      obj = bsUpdate(seqKey) ;
      if (bsAddKey(obj, _DNA, dnaKey))
	bsAddData(obj, _bsRight, _Int, &i) ;
      lexReClass(key1, &key1, _VSequence) ;
      lexReClass(key2, &key2, _VSequence) ;
      if (pre & 1)
	dnaAlignAddPrevious (obj, key1, v1, v2) ;
      else
	{ bsAddKey (obj, _Previous_contig, key1) ;
	  bsAddData (obj, _bsRight, _Int, &v1) ;
	  bsAddData (obj, _bsRight, _Int, &v2) ;
	}
      if (pre & 2)
	dnaAlignAddPrevious (obj, key2, w1, w2) ;
      else
	{ bsAddKey (obj, _Previous_contig, key2) ;
	  bsAddData (obj, _bsRight, _Int, &w1) ;
	  bsAddData (obj, _bsRight, _Int, &w2) ;
	}

      units = arrayCreate (90, BSunit) ;
      if ((Contig = bsCreate(key1)))
	{ if (bsFindTag (Contig, _Assembled_from) &&
	      bsFlatten (Contig, 5, units))
	    { for (i = 0 ; i < arrayMax(units) ; i += 5)
		{ u = arrp(units,i,BSunit) ;
		  key = u[0].k ;
		  p1 = u[1].i ; p2 = u[2].i ;
		  p1 += v1 - 1 ; p2 += v1 - 1 ;
		  bsAddKey (obj, _Assembled_from, key) ;
		  bsAddData (obj, _bsRight, _Int, &p1) ;
		  bsAddData (obj, _bsRight, _Int, &p2) ;
		  x = u[3].i ;
		  if (x)
		    { bsAddData (obj, _bsRight, _Int, &x) ;
		      x = u[4].i ;
		      bsAddData (obj, _bsRight, _Int, &x) ;
		    }
		}
	    }
	  bsDestroy (Contig) ;
	}
      if ((Contig = bsCreate(key2)))
	{ if (bsFindTag (Contig, _Assembled_from) &&
	      bsFlatten (Contig, 5, units))
	    { for (i = 0 ; i < arrayMax(units) ; i += 5)
		{ u = arrp(units,i,BSunit) ;
		  key = u[0].k ;
		  p1 = u[1].i ; p2 = u[2].i ;
		  
		  if (sens == 1)
		    { p1 += w1 - 1 ; p2 += w1 - 1 ; }
		  else
		    { p1 = w1 - p1 + 1 ; p2 = w1 - p2 + 1 ;
		    }
		  bsAddKey (obj, _Assembled_from, key) ;
		  bsAddData (obj, _bsRight, _Int, &p1) ;
		  bsAddData (obj, _bsRight, _Int, &p2) ;
		  x = u[3].i ;
		  if (x)
		    { bsAddData (obj, _bsRight, _Int, &x) ;
		      x = u[4].i ;
		      bsAddData (obj, _bsRight, _Int, &x) ;
		    }
		}
	    }
	  bsDestroy (Contig) ;
	}
      bsSave(obj) ;
    }
  arrayDestroy (units) ;
  arrayDestroy (dna1) ;
  arrayDestroy (dna2) ;
  return seqKey ;
}

/***************************************************************/

Array dnaAlignCompareDna (Array dna1, Array dna2, int *x1, int *x2, int *sens, BOOL countN)
{ int mymax, max1, max2, nn, lg, pp, u1, u2, nbn, recou ;
  Array erreur = 0 ;

  nbn = countN ? 1 : 0 ; /* to get the ambiguities or not */
  max1 = arrayMax (dna1) ;
  max2 = arrayMax (dna2) ;
  mymax = max1 > max2 ? max1 : max2 ;

  if (mymax > 500)
    { lg = 640 ; pp = 7 ; nn = 5 ;
      while (lg < mymax)
	{ pp++ ; lg *= 4 ;
	}
      if (pp > 12) pp = 12 ;	/* saturation du makeShortOligo */
    }
  else
    { pp = 7 ; nn = 5 ; }
  if (dnaAlignAsmbPaireDna (dna1, dna2, pp, 5000, 0, nn, &u1, &u2, sens, TRUE))
    { erreur = arrayCreate (50, A_ERR) ;
      localCptErreur (dna1, 0, max1 - 1, u1, dna2, 0, max2 - 1, u2, *sens,
		      &nbn, x1, x2, &recou, erreur) ;
      (*x1)++ ; (*x2)++ ; /* bio Algebra */
    }
  return erreur ;
}

/***************************************************************/
/* compare les contigs */
BOOL dnaAlignCompare (KEY key1, KEY key2, Array dna1, Array dna2)
{ int x1, x2, xx, sens ;
  OBJ obj = 0 ;
  Array erreur = 0 ;

  if ((erreur = dnaAlignCompareDna (dna1, dna2, &x1, &x2, &sens, FALSE)))
    { arrayDestroy (erreur) ;
      if ((obj = bsUpdate (key1)))
	{ bsAddKey (obj, _Aligned, key2) ;
	  if (sens == - 1)
	    { xx = x1 ; x1 = x2 ; x2 = xx ;
	    }
	  if (bsFindKey (obj, _Aligned, key2))
	    { bsAddData (obj, _bsRight, _Int, &x1) ;
	      bsAddData (obj, _bsRight, _Int, &x2) ;
	      x1 = 1 ;
	      x2 = arrayMax (dna2) ;
	      bsAddData (obj, _bsRight, _Int, &x1) ;
	      bsAddData (obj, _bsRight, _Int, &x2) ;
	    }
	  bsSave (obj) ;
	}
      return TRUE ;
    }
  return FALSE ;
}

/***************************************************************/

int dnaAlignMakeGroups (Array mm, Array *fitp, BOOL direct, int taille)
{ Array ff = arrayCreate (arrayMax(mm), FIT) ;
  int i, j, gg ;
  FIT *f ;

  for (i = 0 , j = 0 ; i < arrayMax(mm) ; j++, i += 2)
    { f = arrayp (ff, j, FIT) ;
      f->u1 = arr(mm, i, int) ;
      f->u2 = arr(mm, i + 1, int) ;
      f->g = 0 ;
    }

  gg = fMapTraceFindBestGroup (ff, direct, 2, taille) ;
/*for (i = 0 ; i < arrayMax(ff) ; i++)
    f = arrayp (ff, i, FIT) ; */
  *fitp = ff ;
  return gg ;
}

/***************************************************************/
#ifndef NON_GRAPHIC
void dnaAlignAsmbPaire (DEFCPT look, KEY key1, KEY key2)
{ int i ;
  KEYSET aa ;
  KEY key ;

  if ((class(key1) != _VDNA && class(key1) != _VSequence) ||
      (class(key2) != _VDNA && class(key2) != _VSequence))
    messout("Please, choose a paire of sequences") ;
  key = dnaAlignPaires(look, key1, key2, &i) ;
  if (!key)
    { messout("Sorry, I have not found a join") ;
      return ;
    }
  messout("I found less than %d %% error", i + 1) ;
  aa = keySetCreate() ;
  keySet(aa, 0) = key1 ;
  keySet(aa, 1) = key2 ;
  keySet(aa, 2) = key ;
  keySetSort(aa) ;
  if (class(key) == _VDNA)
    lexReClass(key, &key, _VSequence) ;
  doForceAssembleSeg (look, key, key, 0, -2, 2, 0) ;
  display (key, 0, 0) ;
  displayCreate(DtKeySet) ;
  graphRetitle("Paire") ;
  keySetShow(aa, 0) ;
}

/***************************************************************/
/* attention mn2 doit etre initialise correctement */
static BOOL dnaAlignUnJoinContig (KEY contigKey, Array contig1, Array pvc1, Array contig2,
				  Array pvc2, int where, int *mx1, int *mn2)
{ int i, j, a1, a2, max1 = 0, maxpc = 0, min2, minpc, min20 ;
  Array units = 0 ;
  KEYSET ks1 = 0, ks2 = 0, ks3 = 0, ks4 = 0 ;
  KEY key ;
  BSunit *u ;
  READ *fp ;
  OBJ obj = bsCreate (contigKey) ;

  if (!obj)
    return FALSE ;
  min2 = min20 = minpc = *mn2 ;
  units = arrayCreate (128, BSunit) ;
  ks1 = keySetCreate () ;
  ks2 = keySetCreate () ;
  if (bsFindTag (obj, _Previous_contig) && bsFlatten (obj, 3, units))
    { for (i = 0 ; i < arrayMax(units) ; i += 3)
	{ u = arrp(units, i, BSunit) ;
	  key = u[0].k ; a1 = u[1].i ; a2 = u[2].i ;
	  if ((a1 < a2 && (where - a1 + 1) > (a2 - where)) ||
	      (a1 > a2 && (where - a2 + 1) > (a1 - where)))
	    { keySet (ks1, keySetMax (ks1)) = key ;
	      fp = arrayp (pvc1, arrayMax (pvc1), READ) ;
	      fp->key = key ;
	      fp->a1 = a1 ;
	      fp->a2 = a2 ;
	      if (a1 > maxpc)
		maxpc = a1 ;
	      if (a2 > maxpc)
		maxpc = a2 ;
	    }
	  else
	    { keySet (ks2, keySetMax (ks2)) = key ;
	      fp = arrayp (pvc2, arrayMax (pvc2), READ) ;
	      fp->key = key ;
	      fp->a1 = a1 ;
	      fp->a2 = a2 ;
	      if (a1 < minpc)
		minpc = a1 ;
	      if (a2 < minpc)
		minpc = a2 ;
	    }
	}
    }
  else goto abort ;
  keySetSort (ks1) ;
  keySetSort (ks2) ;
  ks3 = query (ks1, ">Assembled_from") ;
  ks4 = query (ks2, ">Assembled_from") ;
  minpc -= 100 ; /* in case of bizare realignment */
  maxpc += 100 ;
  if (bsFindTag (obj, _Assembled_from) && bsFlatten (obj, 5, units))
    { for (i = 0 ; i < arrayMax(units) ; i += 5)
	{ u = arrp(units, i, BSunit) ;
	  key = u[0].k ; a1 = u[1].i ; a2 = u[2].i ;
	  if (a1 < minpc && a2 < minpc)
	    { fp = arrayp (contig1, arrayMax (contig1), READ) ;
	      fp->key = key ;
	      fp->a1 = a1 ;
	      fp->a2 = a2 ;
	      fp->ct = u[3].i ;
	      fp->ce = u[4].i ;
	      if (a1 > max1)
		max1 = a1 ;
	      if (a2 > max1)
		max1 = a2 ;
	    }
	  else if (a1 > maxpc && a2 > maxpc)
	    { fp = arrayp (contig2, arrayMax (contig2), READ) ;
	      fp->key = key ;
	      fp->a1 = a1 ;
	      fp->a2 = a2 ;
	      fp->ct = u[3].i ;
	      fp->ce = u[4].i ;
	      if (a1 < min2)
		min2 = a1 ;
	      if (a2 < min2)
		min2 = a2 ;
	    }
	  else if (keySetFind (ks3, key, &j))
	    { fp = arrayp (contig1, arrayMax (contig1), READ) ;
	      fp->key = key ;
	      fp->a1 = a1 ;
	      fp->a2 = a2 ;
	      fp->ct = u[3].i ;
	      fp->ce = u[4].i ;
	      if (a1 > max1)
		max1 = a1 ;
	      if (a2 > max1)
		max1 = a2 ;
	    }
	  else if (keySetFind (ks4, key, &j))
	    { fp = arrayp (contig2, arrayMax (contig2), READ) ;
	      fp->key = key ;
	      fp->a1 = a1 ;
	      fp->a2 = a2 ;
	      fp->ct = u[3].i ;
	      fp->ce = u[4].i ;
	      if (a1 < min2)
		min2 = a1 ;
	      if (a2 < min2)
		min2 = a2 ;
	    }
	  else messerror ("Canot find where to put %s ; dnaAlignUnJoinContig", name (key)) ;
	}
    }
  else goto abort ;
  *mx1 = max1 ;
  *mn2 = min2 ;
 abort:
  keySetDestroy (ks1) ;
  keySetDestroy (ks2) ;
  keySetDestroy (ks3) ;
  keySetDestroy (ks4) ;
  arrayDestroy (units) ;
  bsDestroy (obj) ;
  if (max1 * (min2 - min20) == 0)
    return FALSE ;
  return TRUE ;
}
/*  max1 = 0 ; min20 = min2 ; min2 is length of contigDna */

/***************************************************************/
/* attention mn2 doit etre initialise correctement */
static BOOL dnaAlignDoCutContig (KEY contigKey, Array contig1, Array contig2,
				 int where, int *mx1, int *mn2)
{ int i, a1, a2, max1 = 0, min2, min20 ;
  KEY key ;
  Array units = 0 ;
  BSunit *u ;
  READ *fp ;
  OBJ obj = bsCreate (contigKey) ;

  if (!obj)
    return FALSE ;
  min2 = min20 = *mn2 ;
  units = arrayCreate (128, BSunit) ;
  if (bsFindTag (obj, _Assembled_from) && bsFlatten (obj, 5, units))
    { for (i = 0 ; i < arrayMax(units) ; i += 5)
	{ u = arrp(units, i, BSunit) ;
	  key = u[0].k ; a1 = u[1].i ; a2 = u[2].i ;
	  if ((a1 < a2 && (where - a1 + 1) > (a2 - where)) ||
	      (a1 > a2 && (where - a2 + 1) > (a1 - where)))
	    { fp = arrayp (contig1, arrayMax (contig1), READ) ;
	      fp->key = key ;
	      fp->a1 = a1 ;
	      fp->a2 = a2 ;
	      fp->ct = u[3].i ;
	      fp->ce = u[4].i ;
	      fp->dnaKey = 0 ;
	      if (a1 > max1)
		max1 = a1 ;
	      if (a2 > max1)
		max1 = a2 ;
	    }
	  else
	    { fp = arrayp (contig2, arrayMax (contig2), READ) ;
	      fp->key = key ;
	      fp->a1 = a1 ;
	      fp->a2 = a2 ;
	      fp->ct = u[3].i ;
	      fp->ce = u[4].i ;
	      fp->dnaKey = 0 ;
	      if (a1 < min2)
		min2 = a1 ;
	      if (a2 < min2)
		min2 = a2 ;
	    }
	}
    }
  bsDestroy(obj) ;
  arrayDestroy (units) ;
  if (max1 * (min2 - min20) == 0)
    return FALSE ;

  *mx1 = max1 ;
  *mn2 = min2 ;
  return TRUE ;
}

/***************************************************************/
/* prend une cle sequence et coupe (ou defait un join) a la position determinee */

BOOL dnaAlignCutContig (KEY contigKey, int where, 
			KEY *kp1, KEY *kp2, int *mx1, int *mx2, char action)
{ int i, a1, a2, max1 = 0, max2, min2 ;
  OBJ obj = 0 ;
  KEY seqKey1, seqKey2, key1d, key2d, dnaKey ;
  READ *fp = 0 ;
  BOOL result = FALSE ;
  Array contig1 = 0, contig2 = 0, pvc1 = 0, pvc2 = 0, dna1, dna2, contigDna ;
  char buff[40], *cp, *cq ;

  obj = bsCreate (contigKey) ;
  if (!obj ||
      !bsGetKey (obj, _DNA, &dnaKey) ||
      ! bsGetData (obj, _bsRight, _Int, &min2) ||
      where < 0 || where > min2 - 2 ||
      !(contigDna = dnaGet (dnaKey)))
    { bsDestroy(obj) ;
      messerror ("Bad value for cutting position") ;
      return FALSE ;
    }
  bsDestroy (obj) ;

  contig1 = arrayCreate (48, READ) ;
  contig2 = arrayCreate (48, READ) ;
  switch (action)
    {
    case 'c': /* attention min2 doit etre initialise correctement */
      if (!dnaAlignDoCutContig (contigKey, contig1, contig2, where, &max1, &min2))
	{ messout ("Failed to cut the contig") ;
	  goto abort ;
	}
      break ;
    case 'u': /* attention min2 doit etre initialise correctement */
      pvc1 = arrayCreate (24, READ) ;
      pvc2 = arrayCreate (24, READ) ;
      if (!dnaAlignUnJoinContig (contigKey, contig1, pvc1, contig2, pvc2, where, &max1, &min2))
	{ messout ("Failed to unjoin the contig") ;
	  goto abort ;
	}
      break ;
    default:
      goto abort ;
    }

/* search for a pair of subContig names */
  i = 0 ;
  while(++i)
    if (!lexword2key
	(messprintf("%s_%d", name(contigKey), i), &seqKey1, _VSequence))
      break ;
  lexaddkey 
    (messprintf("%s_%d", name(contigKey), i), &seqKey1, _VSequence) ;
  lexaddkey(messprintf("%s_%d", name(contigKey), i), &key1d, _VDNA) ;
  i++ ;
  sprintf (buff, "%s_%d", name(contigKey), i) ;
  lexaddkey(buff, &seqKey2, _VSequence) ;
  lexaddkey(buff, &key2d, _VDNA) ;

/* create and store dna1 */
  dna1 = arrayCreate(max1, char) ;
  array(dna1, max1 - 1, char) = 0 ;
  cp = arrp(dna1, 0, char) ;
  cq = arrp(contigDna, 0, char) ;
  i = max1 ;
  while(i--)
    *cp++ = *cq++ ;
  dnaStoreDestroy(key1d, dna1) ;

/* create and store dna2 */
  i = max2 = arrayMax(contigDna) - min2 + 1 ;
  dna2 = arrayCreate(i, char) ;
  array(dna2, i - 1, char) = 0 ;
  cp = arrp(dna2, 0, char) ;
  cq = arrp(contigDna, min2 - 1, char) ;
  while(i--)
    *cp++ = *cq++ ;
  dnaStoreDestroy (key2d, dna2) ;

/* reconstruct subcontig1 */
  obj = bsUpdate(seqKey1) ;
  bsAddKey (obj, _DNA, key1d) ;
  bsAddData (obj, _bsRight, _Int, &max1) ;
  i = arrayMax(contig1) ;
  if (i) fp = arrp(contig1, 0, READ) - 1 ;
  while(fp++, i--)
    { bsAddKey (obj, _Assembled_from, fp->key) ;
      bsAddData (obj, _bsRight, _Int, &fp->a1) ;
      bsAddData (obj, _bsRight, _Int, &fp->a2) ;
      if (fp->ct)
	{ bsAddData (obj, _bsRight, _Int, &fp->ct) ;
	  bsAddData (obj, _bsRight, _Int, &fp->ce) ;
	}
    }
  if (pvc1 && (i = arrayMax (pvc1)))
    { fp = arrp (pvc1, 0, READ) - 1 ;
      while (fp++, i--)
	{ bsAddKey (obj, _Previous_contig, fp->key) ;
	  bsAddData (obj, _bsRight, _Int, &fp->a1) ;
	  bsAddData (obj, _bsRight, _Int, &fp->a2) ;
	}
    }
  bsSave(obj) ;

/* reconstruct subcontig2 */
  obj = bsUpdate(seqKey2) ;
  bsAddKey (obj, _DNA, key2d) ;
  bsAddData (obj, _bsRight, _Int, &max2) ;
  i = arrayMax(contig2) ;
  if (i) fp = arrp(contig2, 0, READ) - 1 ;
  while(fp++, i--)
    { bsAddKey (obj, _Assembled_from, fp->key) ;
      a1 = fp->a1 - min2 + 1 ; a2 = fp->a2 - min2 + 1 ;
      bsAddData (obj, _bsRight, _Int, &a1) ;
      bsAddData (obj, _bsRight, _Int, &a2) ;      
      if (fp->ct)
	{ bsAddData (obj, _bsRight, _Int, &fp->ct) ;
	  bsAddData (obj, _bsRight, _Int, &fp->ce) ;
	}
    }
  if (pvc2 && (i = arrayMax (pvc2)))
    { fp = arrp (pvc2, 0, READ) - 1 ;
      while (fp++, i--)
	{ bsAddKey (obj, _Previous_contig, fp->key) ;
	  a1 = fp->a1 - min2 + 1 ; a2 = fp->a2 - min2 + 1 ;
	  bsAddData (obj, _bsRight, _Int, &a1) ;
	  bsAddData (obj, _bsRight, _Int, &a2) ;
	}
    }
  bsSave(obj) ;

/* done */
  *kp1 = seqKey1 ; *mx1 = max1 ;
  *kp2 = seqKey2 ; *mx2 = max2 ;
  defCptForget (0, contigKey) ;
  result = TRUE ;

 abort:
  arrayDestroy (contig1) ;
  arrayDestroy (pvc1) ;
  arrayDestroy (contig2) ;
  arrayDestroy (pvc2) ;
  arrayDestroy (contigDna) ;
  return result ;
}

/***************************************************************/

void dnaAlignFindRepeat (Array source,
			 Array dna, Array color, int x1,  int x2, int taille)
{ int max ;
  Associator oligo = 0 ;

  if (taille < 4)
    { 
      taille = 4 ;
    }
  if (taille > 15)
    { 
      taille = 15 ;
    }
  dnaRepaint (color) ;
  max = arrayMax (source) ;
  if (x1 < 0) x1 = 0 ;
  if (x2 >= max) x2 = max - 1 ;
  if (x2 < x1 + 20) return ;
  oligo = assBigCreate (x2 - x1 + 100) ;
  alignToolsMakeSuccesOligo (source, oligo, x1, x2, taille) ;
  alignToolsFindColorOligo (dna, oligo, color, taille) ;
  assDestroy (oligo) ;
}

/***************************************************************/

void dnaAlignCompareClip (KEY link)
{ int i, ted1, ted2, x, y ;
  KEYSET ks = 0 ;
  KEY *keyp ;
  OBJ obj = 0 ;
  Array ae, ag, af, av ;

  ks = queryKey (link, ">Subsequence ; >Assembled_from") ;
  ae = arrayCreate (1000, int) ;
  ag = arrayCreate (1000, int) ;
  af = arrayCreate (1000, int) ;
  av = arrayCreate (1000, int) ;
  i = keySetMax (ks) ;
  keyp = arrp (ks, 0, KEY) - 1 ;
  while (keyp++, i--)
    { obj = bsCreate (*keyp) ;
      if (!obj || !bsFindTag (obj, _Old_Clipping) ||
	  !bsGetData (obj, _bsRight, _Int, &ted1) ||
	  !bsGetData (obj, _bsRight, _Int, &ted2))
	{ bsDestroy (obj) ;
	  continue ;
	}
      if (bsFindTag (obj, _Excellent_upto) &&
	  bsGetData (obj, _bsRight, _Int, &x))
	{ y = x - ted2 + 500 ;
	  if (y >= 0)
	    array (ae, y, int)++ ;
	}
      if (bsFindTag (obj, _Good_upto) &&
	  bsGetData (obj, _bsRight, _Int, &x))
	{ y = x - ted2 + 500 ;
	  if (y >= 0)
	    array (ag, y, int)++ ;
	}
      if (bsFindTag (obj, _Fair_upto) &&
	  bsGetData (obj, _bsRight, _Int, &x))
	{ y = x - ted2 + 500 ;
	  if (y >= 0)
	    array (af, y, int)++ ;
	}
      if (bsFindTag (obj, _Vector_Clipping) &&
	  bsGetData (obj, _bsRight, _Int, &x))
	{ y = x - ted1 + 100 ;
	  if (y >= 0)
	    array (av, y, int)++ ;
	}
      bsDestroy (obj) ;
    }
  keySetDestroy (ks) ;
  plotHisto ("Excellent / Ted (centre 500)", ae) ;
  plotHisto ("Good / Ted (centre 500)", ag) ;
  plotHisto ("Fair / Ted (centre 500)", af) ;
  plotHisto ("Vector Top / Ted (centre 100)", av) ;
}
#endif /* NON_GRAPHIC */
/***************************************************************/
/* Verifie qu'il est possible de trouver au moins mini dodecameres dans la sequence clippee */
BOOL dnaAlignCheckSequence (Array dna, int cTop, int cEnd,int mini)
{ int i = mini, x = cTop - 1 ;
  unsigned int dummy ;

  while (x < cEnd && i-- && alignToolsMakeOligo (dna, x + 1, cEnd, 12, &dummy, &x)) ;
  if (i == -1) /* ok */
    return TRUE ;
  else return FALSE ;
}

/***************************************************************/

void dnaAlignAdjustLink (KEY link)
{ Array order = 0, units = 0 ;
  int i, j = 0 ;
  BSunit *u ;
  OBJ obj = bsCreate (link) ;

  if (!obj)
    return ;
  units = arrayCreate (100, BSunit) ;
  if (!bsGetArray (obj, _Subsequence, units, 3))
    goto abort ;
  order = arrayCreate (100, BSunit) ;
  for (i = 0 ; i < arrayMax (units) ; i += 3)
    { u = arrp (units, i, BSunit) ;
      array (order, j++, BSunit).k = u[0].k ;
      array (order, j++, BSunit).i = u[2].i - u[1].i ;
    }
  alignToolsAdjustLink (link, 0, order) ;
 abort:
  arrayDestroy (order) ;
  arrayDestroy (units) ;
  bsDestroy (obj) ;
}

/***************************************************************/

BOOL dnaAlignGetClip (KEY link, KEY contig, KEY key, int *ctp, int *cep)
{ int dummy, max ;
  OBJ Contig = 0, Read = 0 ;
  BOOL ok = FALSE ;
  KEY dnaKey ;
  KEYSET ks1, ks2 , ks3 ;

  if (!(Read = bsCreate (key)) || 
      !bsGetKey (Read, _DNA, &dnaKey) ||
      !bsGetData (Read, _bsRight, _Int, &max) ||
      max < 1)
    goto abort ;
	   
  if (!contig && link)
    { ks1 = queryKey (key, "> Assembled_into") ;
      ks2 = queryKey (link,"> Subsequence") ;
      ks3 = keySetAND (ks1, ks2) ;

      keySetDestroy (ks1) ;
      keySetDestroy (ks2) ;
      contig = keySetMax(ks3) > 0 ? keySet(ks3,0) : 0 ;
      keySetDestroy (ks3) ;
    }
  if (contig)
    {
      Contig = bsCreate(contig) ;
      if (Contig && bsFindKey (Contig, _Assembled_from, key) &&
	  bsGetData(Contig, _bsRight, _Int, &dummy) &&
	  bsGetData(Contig, _bsRight, _Int, &dummy) &&
	  bsGetData(Contig, _bsRight, _Int, ctp) &&
	  bsGetData(Contig, _bsRight, _Int, cep) &&
	  0 < *ctp && *ctp < *cep && *cep <=  max
	  )
	ok = TRUE ;
    }
  else
    { 
      if (bsGetData (Read, _Clipping, _Int, ctp) && 
	  bsGetData (Read, _bsRight, _Int, cep) &&
	  0 < *ctp && *ctp < *cep && *cep <=  max) ;
      else
	{ *ctp = 1 ; *cep = max ; }
      ok = TRUE ;
    }

 
 abort:
  bsDestroy (Contig) ;
  bsDestroy (Read) ;
  return ok ;
}

/***************************************************************/

static BOOL abiAlign (KEY target, KEY key, KEY tag)
{ int a1, a2, x1, x2, sens, atp ;
  Array bigDna = dnaGet (target), shDna = dnaGet (key) ;
  OBJ obj = 0 ;
  BOOL resul = FALSE ;
  
  if (!bigDna || !shDna || !dnaAlignGetClip (0, 0, key, &x1, &x2))
    goto abort ;
  a1 = 0 ;
  a2 = arrayMax (bigDna) - 1 ;
  if (!dnaAlignForceMatch (shDna, x1 - 1, x2 - 1, bigDna, a1, a2, 
			   &a1, &a2, &sens))
    goto abort ;
  a1++ ; a2++ ; /* bio algebra */
  if (sens == -1)
    { atp = a1 ; a1 = a2 ; a2 = atp ;
    }
  resul = TRUE ;
  if ((obj = bsUpdate (target)))
    { bsAddKey (obj, tag, key) ;
      bsAddData (obj, _bsRight, _Int, &a1) ;
      bsAddData (obj, _bsRight, _Int, &a2) ;
      bsAddData (obj, _bsRight, _Int, &x1) ;
      bsAddData (obj, _bsRight, _Int, &x2) ;
      bsSave (obj) ;
    }
 abort:
  arrayDestroy (bigDna) ;
  arrayDestroy (shDna) ;
  return resul ;
}

/***************************************************************/

static BOOL newAbiAlign (KEY target, KEY key, KEY tag)
{ int a1, a2, x1, x2, sens, atp ;
  Array bigDna, tpDna, shDna = 0 ;
  OBJ obj = 0 ;
  BOOL resul = FALSE ;

  bigDna = dnaGet (target) ; tpDna = dnaGet (key) ;
  if (!bigDna || !tpDna || !dnaAlignGetClip (0, 0, key, &x1, &x2))
    goto abort ;
  if (x2 > arrayMax (tpDna))
    x2 = arrayMax (tpDna) ; /* kludge : Pb in Clip & quality after deletions */
  a1 = 0 ;
  a2 = arrayMax (bigDna) - 1 ;
  x1-- ; /* Jean algebra (i e x2 exclue) */
  shDna = arrayTruncatedCopy (tpDna, x1, x2) ;
  arrayDestroy (tpDna) ;
  if (!dnaAlignCompareDna (bigDna, shDna, &a1, &a2, &sens, FALSE))
    goto abort ;
  x1++ ; /* bio */
/* a1++ ; a2++ ;  bio algebra ( utile seulement avec abiAlign pas avec newAbiAlign */
  if (sens == -1)
    { atp = a1 ; a1 = a2 ; a2 = atp ;
    }
  resul = TRUE ;
  if ((obj = bsUpdate (target)))
    { bsAddKey (obj, tag, key) ;
      bsAddData (obj, _bsRight, _Int, &a1) ;
      bsAddData (obj, _bsRight, _Int, &a2) ;
      bsAddData (obj, _bsRight, _Int, &x1) ;
      bsAddData (obj, _bsRight, _Int, &x2) ;
      bsSave (obj) ;
    }
 abort:
  arrayDestroy (bigDna) ;
  arrayDestroy (shDna) ;
  return resul ;
}

/***************************************************************/

BOOL dnaAlignContigForcePair (KEY seqKey, KEY target, KEYSET petitKs)
{ int i, j, iPetit, a1, a2, a3, a4 ;
  Array aa = 0, bb = 0 ;
  OBJ obj = 0 ;
  KEY petit ;
  BOOL gotIt = FALSE ;

  iPetit = keySetMax (petitKs) ;
  while (iPetit--)
    { petit = keySet (petitKs, iPetit) ;
      if (!petit)
	continue ;
      if (abiAlign (target, petit, _Aligned)) /* success */
	{ gotIt = TRUE ;
	  keySet (petitKs, iPetit) = 0 ; /* prevent recursions */
	  obj = bsUpdate (target) ;
	  if (bsFindKey (obj, _Aligned, petit) &&
	      bsGetData (obj, _bsRight, _Int, &a1) &&
	      bsGetData (obj, _bsRight, _Int, &a2) &&
	      bsGetData (obj, _bsRight, _Int, &a3) &&
	      bsGetData (obj, _bsRight, _Int, &a4))
	    { bsAddKey (obj, _Assembled_from, petit) ;
	      bsAddData (obj, _bsRight, _Int, &a1) ;
	      bsAddData (obj, _bsRight, _Int, &a2) ;
	      bsAddData (obj, _bsRight, _Int, &a3) ;
	      bsAddData (obj, _bsRight, _Int, &a4) ;
	    }
	  bsSave (obj) ;
	}
    }

  if (gotIt)
    { alignToolsAdjustContigSize (seqKey, target) ;
      aa = arrayCreate (100, BSunit) ;
      bb = arrayCreate (100, BSunit) ;
      if ((obj = bsCreate (seqKey)) && bsFindTag (obj, _Subsequence) &&
	  bsFlatten (obj, 3, bb))
	{ j = 0 ;
	  for (i = 0 ; i < arrayMax (bb) ; i += 3)
	    { array (aa, j++, BSunit).k = arr (bb, i, BSunit).k ;
	      array (aa, j++, BSunit).i = arr (bb, i+2, BSunit).i - arr (bb, i+1, BSunit).i ;
	    }
	  alignToolsAdjustLink (seqKey, 0, aa) ;
	}
      bsDestroy (obj) ;
      arrayDestroy (aa) ;
      arrayDestroy (bb) ;      

    }
  return gotIt ;
}

/***************************************************************/

static void dnaAlignDoCleanLeftArrows (KEY key)
{ OBJ obj = 0 ;

  if ((obj = bsUpdate (key)))
    { if (bsFindTag (obj, _Previous_contig))
	bsRemove (obj) ;
      if (bsFindTag (obj, _Aligned))
	bsRemove (obj) ;
      bsSave (obj) ;
    }
}

/***************************************************************/

void dnaAlignCleanLeftArrows (KEY link, KEY oldCont, KEY newCont, BOOL seton, int choix)
{ int i ;
  static KEYSET ks = 0 ;
  KEY key ;
  static BOOL onoff = FALSE ;

  if (link && seton)
    { onoff = TRUE ;
      if (keySetExists (ks))
	keySetDestroy (ks) ;
      ks = queryKey (link, ">Subsequence") ;
    }
  if (!seton)
    { keySetDestroy (ks) ;
      onoff = FALSE ;
    }
  if (!onoff)
    return ;

  switch (choix)
    {
    case 1:
      if (oldCont && keySetFind (ks, oldCont, &i))
	dnaAlignDoCleanLeftArrows (newCont) ;
      break ;
    case 2:
      if ((i = keySetMax (ks)))
	while (i--)
	  { key = keySet (ks, i) ;
	    dnaAlignDoCleanLeftArrows (key) ;
	  }
      break ;
    case 3:
      dnaAlignDoCleanLeftArrows (newCont) ;
      break ;
    }
}

/***************************************************************/

BOOL dnaAlignDoReInsertLoners (KEY seqKey, KEYSET read)
{ int i, j, k ;
  KEYSET ks, ks1, ks2 ;
  KEY *keyp ;
  OBJ obj = 0 ;
  Array sub = 0, bilan = 0 ;
  BOOL gotIt = FALSE ;

  ks1 = queryKey (seqKey, ">Subsequence") ;
  ks2 = query (ks1, ">Assembled_from") ;
  ks = keySetMINUS (read, ks2) ;
  keySetDestroy (ks2) ;
  i = keySetMax (ks1) ;
  keyp = arrp (ks1, 0, KEY) - 1 ;
  while (keyp++, i--)
    if (dnaAlignContigForcePair (seqKey, *keyp, ks))
      gotIt = TRUE ;
  keySetDestroy (ks) ;
  keySetDestroy (ks1) ;

/* Reinsertion des loners comme contigs */
  ks1 = queryKey (seqKey, ">Subsequence") ; /* not the same object */
  ks2 = query (ks1, ">Assembled_from") ;
  ks = keySetMINUS (read, ks2) ;
  keySetDestroy (ks1) ;
  keySetDestroy (ks2) ;
  ks1 = dnaAlignMakeSubSequence (seqKey, ks, name (seqKey)) ;
  keySetDestroy (ks) ;
  if (keySetMax (ks1))
    gotIt = TRUE ;
  sub = arrayCreate (50, BSunit) ;
  if ((obj = bsCreate (seqKey)) && bsFindTag (obj, _Subsequence) &&
      bsFlatten (obj, 3, sub))
    { bilan = arrayCreate (100, BSunit) ;
      j = 0 ;
      for (k = 0 ; k < arrayMax (sub) ; k += 3)
	{ array (bilan, j++, BSunit).k = arr (sub, k, BSunit).k ;
	  array (bilan, j++, BSunit).i = arr (sub, k+2, BSunit).i - arr (sub, k+1, BSunit).i ;
	}
      for (k = 0 ; k < keySetMax (ks1) ; k++)
	{ array (bilan, j++, BSunit).k = keySet (ks1, k) ;
	  array (bilan, j++, BSunit).i = 1 ;
	}
    }
  arrayDestroy (sub) ;
  bsDestroy (obj) ;
  keySetDestroy (ks1) ;
  if (arrayMax (bilan))
    alignToolsAdjustLink (seqKey, 0, bilan) ;
  arrayDestroy (bilan) ;
  return gotIt ;
}

/***************************************************************/

BOOL dnaAlignReInsertLoners (KEY seqKey, char tag)
{ KEYSET read = 0, ks = 0 ;
  BOOL resul ;

  switch (tag)
    { case 'r':
	if (!messPrompt 
	    ("Which loners do you want to reinsert in the existing contigs ?\n"
	     "The default \"Find Read *\" will reinsert all of them","Find Read *", "t"))
	  return FALSE ;
	read = query (0, freeword()) ;
	break ;
      case 's':
	if (!messPrompt 
	    ("Which Subclones do you want to reinsert in the existing contigs ?\n"
	     "The default \"Find Subclone *\" will reinsert all of them","Find Subclone *", "t"))
	  return FALSE ;
	ks = query (0, freeword()) ;
	read = query (ks, ">Read") ;
	keySetDestroy (ks) ;
	break ;
      }
  resul = dnaAlignDoReInsertLoners (seqKey, read) ;
  keySetDestroy (read) ;
  return resul ;
}

/***************************************************************/

BOOL dnaAlignAssemblyCompare (KEY key1, KEY key2)
{ int i, j ;
  KEY seqKey1, seqKey2 ;
  KEYSET ks1, ks2 ;
  Array dna1, dna2 ;
  BOOL result = FALSE ;

  ks1 = queryKey (key1, ">Subsequence") ;
  if (!keySetMax (ks1)) /* no subsequences => contig or read */
    keySet (ks1, 0) = key1 ;
  ks2 = queryKey (key2, ">Subsequence") ;
  if (!keySetMax (ks2))
    keySet (ks2, 0) = key2 ;
 
  i = keySetMax (ks1) ;
  while (i--)
    { seqKey1 = keySet (ks1, i) ;
      if (!seqKey1 || !(dna1 = dnaGet (seqKey1)))
	continue ;
      if (arrayMax (dna1) < 3000)
	{ arrayDestroy (dna1) ;
	  keySet (ks1, i) = 0 ;
	  continue ;
	}
      j = keySetMax (ks2) ;
      while (j--)
	{ seqKey2 = keySet (ks2, j) ;
	  if (!seqKey2 || !(dna2 = dnaGet (seqKey2)))
	    continue ;
	  if (arrayMax (dna2) < 3000)
	    { arrayDestroy (dna2) ;
	      keySet (ks2, j) = 0 ;
	      continue ;
	    }
	  if (dnaAlignCompare (seqKey1, seqKey2, dna1, dna2))
	    result = TRUE ;
	  arrayDestroy (dna2) ;
	}
      arrayDestroy (dna1) ;
    }
  keySetDestroy (ks1) ;
  keySetDestroy (ks2) ;
  return result ;
}

/***************************************************************/

void doAssembleAllTraces (KEY target, KEYSET ks, char fonc)
{
  KEYSET  ksBad = 0, ks3 = 0 ;
  int i ; /* , a1, a2, x1, x2 ; */
  KEY kk, link = 0, wild, key, tag = _Aligned, *keyp, *keyq ;
  OBJ obj = 0 ;

  /* New Link */
  ks3 = keySetCreate () ;
  keySet (ks3, 0) = target ;
  link = dnaAlignDoMakeSuperLink (ks3, name (target)) ;
  keySetDestroy (ks3) ;
  if (link)
    { ks3 = queryKey (link, ">Subsequence") ;
      target = keySet (ks3, 0) ;
      keySetDestroy (ks3) ;
      tag = _Assembled_from ;
    }
  /* target est un contig maintenant */

  i = keySetMax (ks) ;
  while (i--)
    { kk = keySet (ks, i) ;
      lexReClass (kk, &key, _VSequence) ;
      if (!(obj = bsCreate(key)))
	continue ;
      if (target)
	wild = target ;
      else if (!bsGetKey(obj, _Fragment_of, &wild))
	{ bsDestroy (obj) ;
	  continue ;
	}
      bsDestroy(obj) ;

      if (target)
	wild = target ;
      switch (fonc)
	{
	case 'o':
	  if (!abiAlign (wild, key, tag))
	    { if (!ksBad)
		ksBad = keySetCreate() ;
	      keySet(ksBad, keySetMax(ksBad)) = key ;
	    }
	  break ;
	case 'n':default:
	  if (!newAbiAlign (wild, key, tag))
	    { if (!ksBad)
		ksBad = keySetCreate() ;
	      keySet(ksBad, keySetMax(ksBad)) = key ;
	    }
	  break ;
	}
    }
  if (ksBad)
    { 
      if (!link)
	dnaDispGraph (ksBad, 0) ;
      else
	{ ks3 = dnaAlignMakeSubSequence (link, ksBad, name (link)) ;
	  i = arrayMax (ks3) ;
	  if (i)
	    { keyp = arrayp (ks3, i, KEY) ; /* make room */
	      keyq = keyp + 1 ;
	      while (keyp--, keyq--, i--)
		*keyq = *keyp ;
	      *keyq = target ;
	      alignToolsAdjustLink (link, ks3, 0) ;
	    }
	  keySetDestroy (ks3) ;
	  keySetDestroy (ksBad) ;
	}
    }

#ifndef NON_GRAPHIC 
  if (ksBad && !link)
    dnaDispGraph (ksBad, 0) ;
  if (link)
    wild = link ;
  display (wild, 0, 0) ;
#endif
  keySetDestroy(ks) ;
}

/***************/

#ifndef NON_GRAPHIC
void assembleAllTraces (char fonc)
{ KEYSET ks = 0, aa ;
  int i ; 
  KEY target = 0 ;
  char *cp ;
  void *dummy ;

  sessionGainWriteAccess () ;
  if (!keySetActive(&aa, &dummy))
    { messout("First select a keyset containing dna") ;
      return ;
    }
  ks = query (aa, "CLASS Sequence") ;
  i = keySetMax(ks) ;
  if (!i)
    { messout("First select a keyset containing Sequence") ;
      return ;
    }
    
  if (!messQuery (messprintf("I will try to realign %d traces, should I proceed ?", i)))
    return ;
  
  if (messPrompt("Choose the Target or cancel to use Fragment_of", "", "wz"))
    { cp = freeword() ;
      if (!lexword2key(cp, &target, _VSequence))
	{ messout("Unknown Target") ;
	  return ;
	}
    }
  messStatus ("Aligning fragments" ) ;
  doAssembleAllTraces(target, ks, fonc) ;
}
/***************************************************************/

void oldAssembleAllTraces (void)
{ assembleAllTraces ('o') ;
}

void newAssembleAllTraces (void)
{ assembleAllTraces ('n') ;
}

#endif 
/***************************************************************/
/***************************************************************/
/***************************************************************/

/************************ Pour les stats ***********************/
/* attention taille doit etre < 16 */
static unsigned int dnaAlignShortComplement(unsigned int rac, int taille)
{ unsigned int flag = 3, ret = 0 ;

  while(taille--)
    { ret |= (((~(rac & flag)) & flag) << (2 * taille)) ;
      rac >>= 2 ;
    }
  return ret ;
}

/***************************************************************/

void dnaAlignMakeSpecialMotif (DEFCPT look, Associator adna, int nbolig, int taille)
{ int i, j, max, frac, nd, nf, index = 0 ;
  unsigned int rac, gtr, n ;
  int debut, fin ;
  char *vp, *vp0 = 0 ;
  KEY key ;
  KEYSET dna = look->def, mar = look->mar ;
  Array mondna ;
  char *v0, *vq ;

  i = arrayMax(dna) ;
  frac = ((nbolig + 1) / 2) * 2 ;
  while (i--)
    {
      mondna = monDnaGet (look, 0, keySet(dna, i)) ;
      if (!mondna) continue ;
      debut = - 1 ;
      fin = arrayMax(mondna) ;
      max = fin - 1 ;
      n = nbolig ;
      nd = nf = 0 ;
      while (n--)
	{
	  if (n % 2)
	    { j = max * (frac - nf++) / frac ;
	      if (j >= fin)
		j = fin - 1 ;
	      if (j < taille || 
		  !alignToolsMakeOligo(mondna, j, 0, taille, &rac, &fin))
		continue ;
	    }
	  else
	    {
	      j = max * nd++ / frac ;
	      if (j <= debut)
		j = debut + 1 ;
	      if (j > max - taille ||
		  !alignToolsMakeOligo(mondna, j, max, taille, &rac, &debut))
		continue ;
	    }
	  v0 = vp0 + rac ;
	  if (!(assFind(adna, v0, &vp)))
	    {
	      key = (KEY)rac ;
	      vq = vp0 + key ;
	      keySet(mar, index++) = key ;
	      assInsert(adna, v0, vq) ;
	      gtr = dnaAlignShortComplement(rac, taille) ;
	      v0 = vp0 + gtr ;
	      if (gtr != rac) /* same association for one oligo and its complement */
		assInsert(adna, v0, vq) ;
	    }
	}
    }
  keySetSort(mar) ;
}

/***************************************************************/
