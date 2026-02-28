/*  File: aligntools.c
 *  Author: Ulrich Sauvage (ulrich@kaa.crbm.cnrs-mop.fr)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.ycam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: May 22 14:31 1996 (ulrich)
 * Created: Thu Jun  2 13:04:59 1994 (ulrich)
 *-------------------------------------------------------------------
 */

/* %W% %G% */


#include "acembly.h"
#include "topology.h"
#include "interval.h"
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

#define READ_MIN 40 /* taille minimale des dna alignes (cf monDnaGet) */
#define FIXTOOLMAG 950424121

typedef struct { int pos1, pos2, sens ;} MATCH ;

typedef struct FIXTOOLSTUFF { int magic ;
			      int min, max ;
			      BOOL isFixed ;
			      KEY target, seg, dnakey ;
			      Array seq, dnad, dnar ;
			      Array err, nbt, mod ;
			      Array dummy, alig ;
			    } *FIXTOOL ;

static char u[9] ;
static unsigned int revCompByte[256] ;
static int nbcpterror = 0 ;
static int localCptReport = 0 ;
static void monDnaRemove (DEFCPT look) ;

/***************************************************************/
/*************************** Divers ****************************/
/***************************************************************/

void alignToolsInit (void)
{ unsigned int i, j, k ;
  int n ;

  i = 256 ;
  while(i--)
    { j = 0;
      n = 4 ;
      while(n--)
	{
	  k = (i & (3 << (2*n))) >> (2*n) ;
	  j |= (k << (2*(3 - n))) ;
	}
      revCompByte [i] = (~j) & 255 ;
    }
  i = 9 ; while (i--) u[i] = 0 ;
  u[A_] = 0 ;
  u[T_] = 3 ;
  u[G_] = 1 ;
  u[C_] = 2 ;
}

/***************************************************************/

void alignToolsDestroy (DEFCPT look)
{ if (look)
    monDnaRemove (look) ;
}

/***************************************************************/

unsigned int dnaComplement(unsigned int lambda, int noct)
{ unsigned char c ;
  unsigned int mu = 0 ;
  int n = noct ;

  if (noct <= 0 || noct > sizeof(int))
    messcrash("invalid value noct (%d) in dnaComplement", noct) ;

  if (n>0) while (n--)
    {
      c = lambda & 255 ;
      lambda >>= 8 ;
      mu <<= 8 ;
      mu |= revCompByte[(int)c] ;
    }
  return mu ;
}

/***************************************************************/
static int reportDna = 0, reportBase = 0 ;

void monDnaReportOld (int *nDnap, int *nBasep)
{ *nDnap = reportDna ; *nBasep = reportBase ;
}

Array monDnaGet (DEFCPT look, KEY contig, KEY key)
{ char  *cp, *cq, *cp0 = 0 ;
  Array d ;
  OBJ Seq = 0, Contig = 0 ;
  int cliptop = 0, clipend = 0, i, cVecTop = 0, cVecEnd = 0 ;
  char mask = 15 ; /* Pour nettoyer les flags sur le dna */
  KEY seqKey = 0, dnaKey = 0 ;
  
  if (class (key) ==  _VDNA)
    dnaKey = key ;
  else if (dnaSubClass (key, &dnaKey))
    seqKey = key ;
  else
    return 0 ;

  cp = cp0 + dnaKey ;
  if (!look->assDnaGet)
    look->assDnaGet = assBigCreate (10000) ;

  if (look->assDnaGet && assFind(look->assDnaGet, cp, &cq))
    { d = (Array) cq ;
      if (arrayExists (d))
	return d ;
      invokeDebugger () ;    /* kludge: should not happen */
      assRemove (look->assDnaGet, cp) ;
    }
  
  d = dnaGet(dnaKey) ;
  if (!d)
    return 0 ;
  
  if (!dnaReClass (key, &seqKey)) 
    return d ;
  
  if (contig && class(contig) == _VDNA &&
      !dnaReClass (contig, &contig))
    messcrash ("Bad contig in monDnaGet") ;

  if (!contig && look->link)
    { KEYSET ks1, ks2, ks3 ;
      
      ks1 = queryKey (look->link, "> Subsequence") ;
      ks2 = queryKey (seqKey, "> Assembled_into") ;
      ks3 = keySetAND (ks1, ks2) ;
      
      if (keySetMax(ks3) == 1)
	contig = keySet (ks3, 0) ;
      keySetDestroy (ks1) ;
      keySetDestroy (ks2) ;
      keySetDestroy (ks3) ;
    }
      
  if (seqKey)
    { 
      if (contig && (Contig = bsCreate (contig)) && bsFindKey (Contig, _Assembled_from, seqKey) &&
	  bsGetData (Contig, _bsRight, _Int, &i) && bsGetData (Contig, _bsRight, _Int, &i) &&
	  bsGetData (Contig, _bsRight, _Int, &cliptop) &&
	  bsGetData (Contig, _bsRight, _Int, &clipend)) ; /* Eureka */
      else if ((Seq = bsCreate (seqKey)) && bsGetData (Seq, _Clipping, _Int, &cliptop) &&
	       bsGetData (Seq, _bsRight, _Int, &clipend)) ;
      else if ((Contig || (Contig = bsCreate (look->link))) && bsFindKey (Contig, _Assembled_from, seqKey)&&
	  bsGetData (Contig, _bsRight, _Int, &i) && bsGetData (Contig, _bsRight, _Int, &i) &&
	  bsGetData (Contig, _bsRight, _Int, &cliptop) &&
	  bsGetData (Contig, _bsRight, _Int, &clipend)) ; /* Eureka bis*/
      else
	{ cliptop = 0 ; /* i e no clip */
	  clipend = arrayMax (d) ;
	}
    }
  else 
    { cliptop = 0 ; /* should not happen */
      clipend = arrayMax (d) ;
    }
  if (Seq && bsGetData (Seq, _Vector_Clipping, _Int, &cVecTop))
    bsGetData (Seq, _bsRight, _Int, &cVecEnd) ;
  bsDestroy (Contig) ;
  bsDestroy (Seq) ;

  if (cliptop < 0) cliptop = 0 ;
  if (cliptop < cVecTop)
    cliptop = cVecTop ;
  if (clipend > arrayMax(d))
    clipend = arrayMax(d) ;
  if (cVecEnd && cVecEnd < clipend)
    clipend = cVecEnd ;
  if (clipend < cliptop)
    goto abort ; /* rejected vector */
  if (cliptop)
    cliptop-- ; /* decalage ieme base <-> array(i-1) */
  if (!dnaAlignCheckSequence (d, cliptop, clipend, 4)) /* check if we may find 4 dodecamers */
    goto abort ;
  arrayMax(d) = clipend - cliptop ; /* clip changes the Max */
  i = arrayMax(d) ;
  if (i < READ_MIN)
    goto abort ;
  cp = arrp (d, 0, char) ;
  cq = cp + cliptop ;
  while (i--)
    *cp++ = (mask & *cq++) ; /* to clean the flag */
  *cp = 0 ;
  reportDna++ ; reportBase += arrayMax (d) ;
  cp = cp0 + dnaKey ;
  if (look->assDnaGet)
    assInsert (look->assDnaGet, cp, (void*)d) ;

  return d ;
 abort:
  arrayDestroy (d) ;
  return 0 ;             /* reject */
}

/***************************************************************/

void monDnaInsert(DEFCPT look, KEY key, Array dna)
{ Array dnaInter, dummy = 0 ;
  char *kp = 0 ;

  dnaInter = dnaCopy(dna) ;
  dnaStoreDestroy(key, dnaInter) ;
/* may i need to forget ? */
  if (!look->assDnaGet)
    look->assDnaGet = assBigCreate (10000) ;
  if (look->assDnaGet && assFind (look->assDnaGet, kp + key, &dummy))
    { if (dummy != dna && arrayExists (dummy))
	{ messerror ("double insert in mondnaget") ; 
	  arrayDestroy (dummy) ;
	}
      assRemove (look->assDnaGet, kp + key) ;
    }
  if (look->assDnaGet)
    assInsert (look->assDnaGet, kp + key, (void*)dna) ;
  reportDna++ ; reportBase += arrayMax (dna) ;
}

/***************************************************************/

void monDnaForget(DEFCPT look, KEY key)
{ void *vp = 0 ;
  char *kp = 0 ;
  Array d ;
  KEY mykey ;

  if (class(key) != _VDNA)
    { if (!dnaSubClass(key, &mykey))
	return ;
    }
  else
    mykey = key ;
  kp += mykey ;
  if (!look || !look->assDnaGet || !assFind (look->assDnaGet, kp, &vp))
    return ;
  d = (Array)vp ;
  if (!arrayExists(d))
    messcrash("See Ulrich in dnaAlignForget") ;
  reportDna-- ; reportBase -= arrayMax (d) ;
  arrayDestroy(d) ;
  if (!assRemove(look->assDnaGet, kp))
    messcrash("in dnaAlignForget : Pb with assRemove") ;
}

/***************************************************************/

static void monDnaRemove (DEFCPT look)
{ void *vp, *wp ;
  Array d ;

  if (!look || !assExists(look->assDnaGet))
    return ;
  vp = 0 ;
  while(assNext(look->assDnaGet, &vp, &wp))
    { d = (Array)wp ;
      if (arrayExists(d))
	arrayDestroy(d) ;
    }
  assDestroy (look->assDnaGet) ;
}

/***************************************************************/

void alignToolsDestroy_Segs (int id)
{ KEY key ;
  OBJ obj ;
  int i ;
  KEYSET ks = 0, ksd = 0 ;

  if (id)
    { ks = query (0, messprintf
		  ("FIND Sequence IS _segment_%d.* | IS _Assembly_%d.*", id, id)) ;
      ksd = query (0, messprintf
		   ("FIND DNA IS _segment_%d.* | IS _Assembly_%d.*", id, id)) ;
    }
  else
    { ks = query (0, "FIND Sequence IS _segment* | IS _Assembly*") ;
      ksd = query (0, "FIND DNA IS _segment* | IS _Assembly*") ;
    }
  i = keySetMax (ks) ;
  while (i--)
    { key = keySet (ks, i) ;
      if ((obj = bsCreate (key)))
	{ bsDestroy (obj) ;
	  obj = bsUpdate (key) ;
	  bsKill (obj) ;
	}
    }
  keySetDestroy (ks) ;
  i = keySetMax (ksd) ;
  while (i--)
    { key = keySet (ksd, i) ;
      arrayKill(key) ;
    }
  keySetDestroy (ksd) ;
}

/***************************************************************/

/* attention il faut au plus des hexadecameres */
/* Pour le moment toutes autres fonctions sont sur dodecameres */
BOOL alignToolsMakeOligo(Array dna, int debut, int fin, int taille, 
			  unsigned int *olig, int *pos)
{ char *cp, *cp0 ;
  int i, n = taille, max, decal, sens = 1 ;
  unsigned int rac = 0, mask = (1 << (taille * 2))  - 1 ;
  
  cp0 = arrp(dna, 0, char) ;
  max = arrayMax(dna) ;
  if (max && debut >= max && debut < max + 10) debut = max - 1 ; /* very frequent */

    
  if (debut >= max) debut = max - 1 ;
  if (debut < 0) debut = 0 ;
  
  if (max && fin >= max && fin < max + 10) fin = max - 1 ; /* very frequent */
  if (fin >= max || fin < 0)
    { /* printf("// AlignToolsMakeOligo : bad value for fin = %d > max = %d\n", fin, max) ; */
      fin = max - 1 ;
    }


  if (fin < debut)
    sens = - 1 ;
  if (sens == 1) /* oligo en montant */
    { cp = arrp(dna, debut, char) - 1 ;
      i = fin - debut + 1 ;
      if (i>0) while(cp++, i-- && n)
	{ switch(*cp & 0x0f)
	    {
	    case A_: case T_: case G_: case C_:
	      rac = ((rac << 2) | u[(int)*cp]) ;
	      n-- ;
	      if (!n && i) /* pour eviter de rejeter une sequence a cause d'un poly A ou T */
		{ rac &= mask ;
		  if (!rac || rac == mask)
		    n++ ;
		}
	      break ;
	    default:
	      rac = 0 ;
	      n = taille ;
	      break ;
	    }
	}
      cp -= taille ;
    }
  else /* oligo en descendant */
    { cp = arrp(dna, debut, char) + 1 ;
      i = debut - fin + 1 ;
      decal = (taille - 1) * 2 ;
      if (i>0) while(cp--, i-- && n)
	{ switch(*cp & 0x0f)
	    {
	    case A_: case T_: case G_: case C_:
	      rac = ((rac >> 2) | (u[(int)*cp] << decal)) ;
	      n-- ;
	      if (!n && i) /* pour eviter de rejeter une sequence a cause d'un poly A ou T */
		{ rac &= mask ;
		  if (!rac || rac == mask)
		    n++ ;
		}
	      break ;
	    default:
	      rac = 0 ;
	      n = taille ;
	      break ;
	    }
	}
      cp++ ;
    }
  *olig = rac & mask ;
  *pos = cp - cp0 ; /* pos en coordonne array : 0 -> max - 1 */
  if (n || *olig == 0 || *olig == mask)
    return FALSE ;
  else
    return TRUE ;
}

/***************************************************************/

static void alignToolsFindOligo (Array dna, int xd, int xf, Array matching,
				 int pos2, unsigned int match)
{ int i, j, pos, n ;
  unsigned int mask = 1 << 24 ; /* pour avoir des dodecameres */
  char *cp ;
  unsigned int matchcomp, rac = 0 ;

  matchcomp = dnaComplement(match, 3) ; /* dodecamere */
  i = xf - xd + 1 ;
  if (i<0 || xd < 0 || xf > arrayMax(dna)) 
    messcrash ("i=%d in alignToolsFindOligo",i) ;
  cp = arrp(dna, xd, char) - 1 ;
  pos = xd ;
  n = 0 ;
  while (pos++, cp++, i--)
    {
       switch (*cp)
	{
	case A_: case T_: case G_: case C_:
	  rac = ((rac << 2) | u[(int)*cp]) ;
	  n++ ;
	  break ;
	default:
	  n = 0 ;
	  continue ;
	}
       if (n >= 12)
	 {
	   rac &= mask - 1 ;
	   if (rac == match)
	     { j = arrayMax(matching) ;
	       array(matching, j, MATCH).pos1 = pos - 12 ; /* attention dodecameres */
	       arr(matching, j, MATCH).pos2 = pos2 ;
	       arr(matching, j, MATCH).sens = 1 ;
	     } 
	   if (rac == matchcomp)
	     { j = arrayMax(matching) ;
	       array(matching, j, MATCH).pos1 = pos - 1 ;
	       arr(matching, j, MATCH).pos2 = pos2 ;
	       arr(matching, j, MATCH).sens = - 1 ;
	     }
	 }
     }
}

/***************************************************************/
/* find avec taille variable mais sens impose ! */
/* xd et xf positions incluses dans le dna pour la zone de calcul */
Array alignToolsFindShortOligo(Array dna, int xd, int xf, 
				      unsigned int win, int taille,
				      unsigned int olig)
{ int i, n, j = 0 ;
  char *cp, *cp0 ;
  unsigned int rac = 0 ;
  Array result = 0 ;

  i = xd < 0 ? -xd : 0 ;
  cp0 = arrp(dna, 0, char) ;
  cp = cp0 + xd + i ;
  i = xf - xd + 1 ;
  n = 0 ;
  result = arrayCreate(30, int) ;
  if (i>0) while(i--)
    {
      switch(*cp)
	{
	case A_: case T_: case G_: case C_:
	  rac = ((rac << 2) | u[(int)*cp++]) ;
	  n++ ;
	  break ;
	default:
	  cp++ ;
	  n = 0 ;
	  continue ;
	}
      if (n >= taille)
	{
	  rac &= win ;
	  if (rac == olig)
	    array(result, j++, int) = cp - cp0 - taille ;
	}
    }
/*  if (!j)
    arrayDestroy(result) ; */
 /* mieg probleme si destroy dnas calling routine */
/* Pourquoi ne pas detruire le tableau s'il est vide */
  return result ;
}

/***************************************************************/
/* attention taille doit etre < 16 */
static unsigned int shortComplement(unsigned int rac, int taille)
{ unsigned int flag = 3, ret = 0 ;

  while(taille--)
    { ret |= (((~(rac & flag)) & flag) << (2 * taille)) ;
      rac >>= 2 ;
    }
  return ret ;
}

/***************************************************************/

static BOOL invariantOligo (char *cp0, int taille)
{ char *cp, *cq ;
  int i, j ;
  char buf[64] ;

  if (taille > 16) messcrash("taille > 16 in invariantoligo") ;
  
  cp = buf ; cq = cp0 ; i = taille ;
  while (i--) *cp++ = *cq++ ;

  cq = cp0 ; i = taille ; /* second copy */
  while (i--) *cp++ = *cq++ ;
    
  for (i = 1 ; i < taille ; i++)
    { cp = &buf[i] ; cq = cp0 ;
      j = taille ;
      while (j--)
	if (*cp++ != *cq++)
	  goto ok1 ;
      return FALSE ;
    ok1:
      continue ;
    }
  return TRUE ;
}

/***************************************************************/
/* ordre normal long, court */
Array alignToolsMakeShortMatch(Array dna1, Array dna2, int sens, int taille,
			       int max, int zone)
{ Array paire = 0, dna ;
  BOOL reverse = FALSE ;
  char *cp, *cp0, *asp = 0, *vp ;
  Associator exist = 0 ;
  int i, j, n, nb = 0, pos, pos1, pos2, assmax, saut = 0, jump ;
  unsigned int acgt = 0 , tgca, win ;

  if (taille > 12)
    messcrash("Dodecamere Max") ;
  win = (1 << (2 * taille)) - 1 ;
  if (arrayMax(dna1) < arrayMax(dna2))
    { reverse = TRUE ;
      dna = dna1 ;
      dna1 = dna2 ;
      dna2 = dna ; 
      dna = 0 ;
    }
  jump = zone > arrayMax (dna2) ? arrayMax (dna2) : zone ;
  if (jump > 500)
    saut = jump/500 ;
  jump = saut ;
  exist = assBigCreate(2 * arrayMax(dna2)) ;
  cp0 = arrp(dna2, 0, char) ;
  cp = cp0 - 1 ;
  if (zone <= 0 || 2 * zone > arrayMax(dna2))
    { j = 1 ;
      i = arrayMax(dna2) ;
    }
  else
    { j = 2 ;
      i = zone ;
    }
  n = 0 ;
  while(j--) /* boucle pour ne prendre que les bouts */
    { if (i>0) while(cp++, i--)
	{
	  switch(*cp)
	    {
	    case A_: case C_: case G_: case T_:
	      acgt = ((acgt << 2) | u[(int)*cp]) ;
	      n++ ; jump++ ;
	      break ;
	    default:
	      n = 0 ;
	      break ;
	    }
	  if (n >= taille && jump > saut)
	    {
	      acgt &= win ;
	      if (!acgt || !invariantOligo(cp - taille + 1, taille))
		continue ;
	      pos = cp - cp0 - taille + 1 ;
	      tgca = acgt << 8 ;
	      assmax = 254 ;
	      while(assmax-- && assFind(exist, asp + tgca, &vp))
		tgca++ ;
	      assInsert(exist, asp + tgca, asp + pos) ;
	      jump = 0 ;
	    }
	}
      i = zone ;
      n = 0 ;
      cp = cp0 + arrayMax(dna2) - zone - 1 ;
    }
  paire = arrayCreate(1000, int) ;
  i = arrayMax(dna1) ;
  n = 0 ;
  cp0 = arrp(dna1, 0, char) ;
  cp = cp0 - 1 ;
  if (zone <= 0 || 2 * zone > arrayMax(dna1))
    { j = 1 ;
      i = arrayMax(dna1) ;
    }
  else
    { j = 2 ;
      i = zone ;
    }
  while(j--)
    { if (i>0) while(cp++, i--)
	{
	  switch(*cp)
	    {
	    case A_: case C_: case G_: case T_:
	      acgt = ((acgt << 2) | u[(int)*cp]) ;
	      n++ ;
	      break ;
	    default:
	      n = 0 ;
	      break ;
	    }
	  if (n >= taille)
	    {
	      acgt &= win ;
	      if (sens == 1)
		tgca = acgt ;
	      else
		tgca = shortComplement(acgt, taille) ;
	      tgca <<= 8 ;
	      if (!acgt || !invariantOligo(cp - taille + 1, taille) || 
		  !assFind(exist, asp + tgca, &vp))
		continue ;
	      if (sens == 1)
		pos1 = cp - cp0 - taille + 1 ;
	      else
		pos1 = cp - cp0 ;
	      assmax = 254 ;
	      while(assmax-- && assFind(exist, asp + tgca, &vp))
		{ tgca++ ;
		  pos2 = vp - asp ;
		  if (reverse)
		    { array(paire, nb++, int) = pos2 ;
		      array(paire, nb++, int) = pos1 ;
		    }
		  else
		    { array(paire, nb++, int) = pos1 ;
		      array(paire, nb++, int) = pos2 ;
		    }
		  if (nb > max)
		    break ;
		}
	    }
	  if (nb > max)
	    break ;
	}
      i = zone ;
      n = 0 ;
      cp = cp0 + arrayMax(dna1) - zone - 1 ;
    }
  assDestroy (exist) ;
  if (!arrayMax(paire))
    arrayDestroy(paire) ;
  return paire ;
}

/***************************************************************/

void alignToolsMakeSuccesOligo (Array dna, Associator oligo, int x1, 
				int x2, int taille)
{ char *cp, *cp0, *cq, *cr ;
  int i, n = taille ;
  unsigned int rac = 0, gtr ;
  void *vp ;

  i = x2 - x1 ;
  cp0 = arrp (dna, 0, char) + taille + 1 ;
  cp = arrp (dna, x1, char) - 1 ;
  if (i>0) while (cp++, i--)
    { switch (*cp & 0x0f)
	{ 
	case A_: case C_: case G_: case T_:
	  rac = ((rac << 2) | u[(int)*cp]) ;
	  n-- ;
	  break ;
	default:
	  rac = 0 ;
	  n = taille ;
	  break ;
	}
      if (!n)
	{ cq = (char *)0 + rac ;
	  if (!assFind (oligo, cq, &vp))
	    { cr = (char *)0 + (cp - cp0) ;
	      if (cq) assInsert (oligo, cq, cr) ;
	      gtr = shortComplement (rac, taille) ;
	      if (gtr != rac)
		{ cq = (char *)0 + gtr ;
		  if (cq) assInsert (oligo, cq, cr) ;
		}
	    }
	  n = taille ;
	  rac = 0 ;
	}
    }
}

/***************************************************************/


#define TINT_HIGHLIGHT1 0x01	/* highest priority, dna highlighting */
#define TINT_HIGHLIGHT2 0x02	/* highlight friends */
#define TINT_RED        0x04 
#define TINT_LIGHTGRAY  0x08 
#define TINT_MAGENTA    0x10 
#define TINT_CYAN       0x20 
#define TINT_LIGHTGREEN 0x40 
#define TINT_YELLOW     0x80 

void alignToolsFindColorOligo (Array dna, Associator oligo, Array color,
			       int taille)
{ char *cp, *cq, *cp0, *cq0, *vp1, *vp2, *cr, mycolor[5] ;
  int i, j, n = 0, mask = (1 << (2 * taille)) - 1, col ;
  unsigned int rac = 0 , gtr ;
  void *vp ;

  cr = mycolor ;
  *cr++ = TINT_RED ;
  *cr++ =  TINT_LIGHTGRAY ;
  *cr++ =  TINT_MAGENTA ;
  *cr++ =  TINT_CYAN ;
  *cr++ =  TINT_LIGHTGREEN ;
  cp = arrp (dna, 0, char) - 1 ;
  cp0 = arrp (dna, 0, char) ;
  cq0 = arrp (color, 0, char) ;
  i = arrayMax (dna) ;
  if (i != arrayMax (color))
    messcrash ("arrayMax (color) not good") ;
  while (cp++, i--)
    { switch (*cp & 0x0f)
	{
	case A_: case C_: case G_: case T_:
	  rac = ((rac << 2) | u[(int)*cp]) ;
	  n++ ;
	  break ;
	default:
	  n = 0 ;
	  break ;
	}
      if (n >= taille)
	{ rac &= mask ;
	  gtr = shortComplement (rac, taille) ;
	  vp1 = (char *)0 + rac ;
	  vp2 = (char *)0 + gtr ;
	  if (assFind (oligo, vp1, &vp) || assFind (oligo, vp2, &vp))
	    { cq = cq0 + (cp - cp0) ;
	      col = (char *)vp - (char *)0 ;
	      col /= 50 ; col %= 5 ;
	      j = taille ;
	      while (j--)
		*cq-- |=  mycolor[col] ;
	    }
	}
    }
}

/***************************************************************/
/************************** Erreurs ****************************/
/***************************************************************/
/* reCreates errArray */
/* reCreates errArray */
void localCptErreur(Array longDna, int xl1, int xl2, int pol,
		    Array shortDna, int xs1, int xs2, int pos, int sens,
		    int *nbN, int *start, int *stop, int *recouv, Array errArray)
{ 
  int 
    err = 0, i, j, i1=0, j1=0, j2=0, type, 
    xsd, xsf, all = *nbN, n = 0, pass,
    bt, bx, b1 ;  /* best type best dx */
  char *cp, *cq, *cp1=0, *cq1=0, *cp0, *cq0, base[16] ;
  A_ERR *newErr ;

  if (!arrayExists (errArray))
    messcrash ("localCptErreur receive a non Existing Array") ;
  if (xl2 >= arrayMax (longDna))
    { /* messerror ("xlong2 depasse le dna") ; */
      xl2 = arrayMax (longDna) - 1 ;
    }
  nbcpterror++ ;
  localCptReport++ ;
  i = 16 ;
  if (sens == 1)
    while(i--)
      base[i] = i ;
  else
    while(i--)
      base[i] = complementBase(i) ;
  cp0 = arrp(longDna, 0, char) ;
  cq0 = arrp(shortDna, 0, char) ;
  cp = cp0 + pol ;
  cq = cq0 + pos ;
  i = pol - xl1 + 1 ;
  j = sens == 1 ? pos - xs1 + 1 : xs2 - pos + 1 ;
  xsd = sens == 1 ? xs1 : xs2 ;
  xsf = sens == 1 ? xs2 : xs1 ;
  while(i && j) /* recherche de la position de depart */
    { if (i < 0) /* may happen in DNA _Aligned to a contig */
	break ;
      i-- ; j-- ; /* decrementation des deux ou aucun */
      if (*cp && *cq && !(*cp & base[(int)*cq]))
	{ type = 6 ; /* AMBIGUE non teste donc 0 non teste */
	  pass = 1 ; bt = 0 ; bx = 10000 ;
	  while(--type)
	    { switch (type)
		{
		case ERREUR:
		  cp1 = cp - 1 ;
		  cq1 = cq - sens ;
		  i1 = i ;
		  j1 = j ;
		  if (!i1 || !j1)
		    continue ;
		  j2 = GO_AFT_ERREUR ;
		  break ;
		case TROU:
		  cp1 = cp - 1 ;
		  cq1 = cq ;
		  i1 = i ;
		  if (!i1)
		    continue ;
		  j1 = j + 1 ;
		  if (*cp1 & base[(int)*cq1])
		    j2 = GO_AFT_SINGLE ;
		  else
		    continue ;
		  break ;
		case TROU_DOUBLE:
		  cp1 = cp - 1 ;
		  cq1 = cq ;
		  if (*cp1 & base[(int)*cq1])
		    continue ;  /* no double single was tried */
		  cp1-- ;
		  i1 = i - 1 ;
		  if (i1 <= 0)
		    continue ;
		  j1 = j + 1 ;
		  j2 = GO_AFT_DOUBLE ;
		  break ;
		case INSERTION:
		  cp1 = cp ;
		  cq1 = cq - sens ;
		  i1 = i + 1 ;
		  j1 = j ;
		  if (!j1)
		    continue ;
		  if (*cp1 & base[(int)*cq1])
		    j2 = GO_AFT_SINGLE ;
		  else
		    continue ;
		  break ;
		case INSERTION_DOUBLE:
		  cp1 = cp ;
		  cq1 = cq - sens ;
		  if (*cp1 & base[(int)*cq1]) /* no double single was tried */
		    continue ;
		  cq1 -= sens ;
		  i1 = i + 1 ;
		  j1 = j - 1 ;
		  if (j1 <= 0)
		    continue ;
		  j2 = GO_AFT_DOUBLE ;
		  break ;
		}
	      b1 = 0 ;
	      if (j2) while(b1++ < 500 && i1-- && j1-- && j2--)
		{ 
		  if (!*cp1 || !*cq1 || *cp1 == N_ || *cq1 == N_)
		    j2++ ;
		  else if (!(*cp1 & base[(int)*cq1]))
		    { if (type == INSERTION || type == INSERTION_DOUBLE)
			{ j2 += 2 ;
			  if (pass || j2 > 24) /* 24 errors is enough ! */
			    break ;
/*		insertion seulement si ok
   	  if (*cp1 & base[(int)*(cq1 - sens)])
			    { cq1 -= sens ; j1-- ; }
*/

			  if (*(cp1 - 1) & base[(int)*(cq1 - sens)]) ;  /* ponctuel */
			  else if (type == INSERTION_DOUBLE || /* avance coute que coute */
				   (*cp1 & base[(int)*(cq1 - sens)]) )
			    { j2++ ; cq1 -= sens ; j1-- ; }
			}
		      else if (type == TROU || type == TROU_DOUBLE)
			{ j2 += 2 ;
			  if (pass || j2 > 24) /* 24 errors is enough ! */
			    break ;
			  if (*(cp1 - 1) & base[(int)*(cq1 - sens)]) ; /* ponctuel */
			  else if (type == TROU_DOUBLE || /* avance coute que coute */
				   (*(cp1 - 1) & base[(int)*cq1]) )
			    { j2++ ; cp1-- ; i1-- ; }
			}
		      else
			break ;
		    } 
		  cp1-- ;
		  cq1 -= sens ;
		}
	      if (j2 < 2 && (j2 == -1 || j1 == -1 || i1 == -1)) /* so the error is jumped */
		break ;
	      if (!pass && b1 < bx) { bt = type ; bx = b1 ; }
	      if (pass && type == 3)
		{ pass = 0 ; type = 6 ; } /* was type = 5 */
	    }
	  if (!type)
	    {
	      if (!pass)
		type = bt ;
	      else
		type = ERREUR ; 
	    }
	  switch(type)
	    {
	    case TROU_DOUBLE:
	      cp -= 2 ;
	      i -= 2 ;
	      break ;
	    case TROU:
	      cp-- ;
	      i-- ;
	      break ;
	    case INSERTION_DOUBLE:
	      cq -= 2 * sens ;
	      j -= 2 ;
	      break ;
	    case INSERTION:
	      cq -= sens ;
	      j-- ;
	      break ;
/* en fait ambigue non testes mais correspond a : rien d'autre ne marche */
	    case ERREUR:case AMBIGUE:
	      break ;
	    }
	}
      cp-- ;
      cq -= sens ;
    }      /* en fin de boucle je suis 1 en dessous des derniers char testes */
  cp++ ;   /* positionnement sur le depart */
  cq += sens ;
  if (!i && xl1 != cp - cp0)
    messcrash("erreur dans les indices cp dans CptErreur") ;
  if (!j && xsd != cq - cq0)
    messcrash("erreur dans les indices cq dans CptErreur") ;
  *start = xl1 + i - j ;
  errArray = arrayReCreate(errArray, 50, A_ERR) ;
  cp = cp0 + *start ;
  cq = cq0 + xsd ;
  j -= i ; /* it may be < 0 in _Aligned */
  if (j > 0)
    {
      if (all && j) /* ambigue en dehors de longDna */
	{ cp-- ; cq -= sens ;
	  while (cp++, cq += sens, j--)
	    if (!*cq || *cq == N_)
	      { n++ ;
		newErr = arrayp(errArray, err++, A_ERR) ;
		newErr->type = AMBIGUE ;
		newErr->iLong = cp - cp0 ;
		newErr->iShort = cq - cq0 ;
		newErr->baseLong = *cp ;
		newErr->baseShort = *cq ;
		newErr->sens = sens ;
	      }
	}
      else if (j)
	{ cp += j ;
	  cq += j * sens ;
	}
    }
  *recouv = cp0 - cp ; /* ce n'est pas la valeur finale */
  i = xl2 - (cp - cp0) + 1 ;
  if (sens == 1)
    j = xsf - (cq - cq0) + 1 ;
  else
    j = cq - cq0 - xsf + 1 ;
  while(i && j)  /* calcul du tableau d'erreur */
    { i-- ; j-- ;
      if (all && (!*cq || *cq == N_)) /* enregistre ambigue */
	{ n++ ;
	  newErr = arrayp(errArray, err++, A_ERR) ;
	  newErr->type = AMBIGUE ;
	  newErr->iLong = cp - cp0 ;
	  newErr->iShort = cq - cq0 ;
	  newErr->baseLong = *cp ;
	  newErr->baseShort = *cq ;
	  newErr->sens = sens ;
	}
      if (*cp && *cq && !(*cp & base[(int)*cq]))
	{ type = 6 ; /* Ambigue non erreur donc 0 non teste */
	  pass = 1 ; bt = 0 ; bx = 10000 ;
	  while(--type)
	    { switch (type)
		{
		case ERREUR:
		  cp1 = cp + 1 ;
		  cq1 = cq + sens ;
		  i1 = i ;
		  j1 = j ;
		  if (!i1 || !j1)
		    continue ;
		  j2 = GO_AFT_ERREUR ;
		  break ;
		case TROU:
		  cp1 = cp + 1 ;
		  cq1 = cq ;
		  i1 = i ;
		  if (!i1)
		    continue ;
		  j1 = j + 1 ;
		  if (*cp1 & base[(int)*cq1])
		    j2 = GO_AFT_SINGLE ;
		  else
		    continue ;
		  break ;
		case TROU_DOUBLE:
		  cp1 = cp + 1 ;
		  cq1 = cq ;
		  if (*cp1 & base[(int)*cq1])
		    continue ;  /* no double single was tried */
		  cp1 ++ ;
		  i1 = i - 1 ;
		  if (i1 <= 0)
		    continue ;
		  j1 = j + 1 ;
		  j2 = GO_AFT_DOUBLE ;
		  break ;
		case INSERTION:
		  cp1 = cp ;
		  cq1 = cq + sens ;
		  i1 = i + 1 ;
		  j1 = j ;
		  if (!j1)
		    continue ;
		  if (*cp1 & base[(int)*cq1])
		    j2 = GO_AFT_SINGLE ;
		  else
		    continue ;
		  break ;
		case INSERTION_DOUBLE:
		  cp1 = cp ;
		  cq1 = cq + sens ;
		  if (*cp1 & base[(int)*cq1]) /* no double single was tried */
		    continue ;
		  cq1 += sens ;
		  i1 = i + 1 ;
		  j1 = j - 1 ;
		  if (j1 <= 0)
		    continue ;
		  j2 = GO_AFT_DOUBLE ;
		  break ;
		}
	      b1 = 0 ;
	      if (j2) while(b1++ < 500 && i1-- && j1-- && j2--)
		{ 
		  if (!*cp1 || !*cq1 || *cp1 == N_ || *cq1 == N_)
		    j2++ ;
		  else if (!(*cp1 & base[(int)*cq1]))
		    { 
		      if (type == INSERTION || type == INSERTION_DOUBLE)
			{ j2 += 2 ;
			  if (pass || j2 > 24) /* 24 errors is enough ! */
			    break ;
			  if (*(cp1 + 1) & base[(int)*(cq1 + sens)]) ;
			  else if (type == INSERTION_DOUBLE || /* avance coute que coute */
				   (*cp1 & base[(int)*(cq1 + sens)]) )
			    { j2++ ; cq1 += sens ; j1-- ; }
			}
 		      else if (type == TROU || type == TROU_DOUBLE)
			{ j2 += 2 ;
			  if (pass || j2 > 24) /* 24 errors is enough ! */
			    break ;
			  if (*(cp1 + 1) & base[(int)*(cq1 + sens)]) ; /* ponctuel */
			  else if (type == TROU_DOUBLE || /* avance coute que coute */
				   (*(cp1 + 1) & base[(int)*cq1]) )
			    { j2++ ; cp1++ ; i1-- ; }
			}
		      else /* erreur */
			{ j2 += 2 ;
			  if (pass || j2 > 12) /* 12 errors is enough ! */
			    break ;
			}
		    } 
		  cp1++ ;
		  cq1 += sens ;
		}
	      if (j2 < 2 && (j2 == -1 || j1 == -1 || i1 == -1)) /* so the error is jumped */
		break ;
	      if (!pass && b1 < bx) { bt = type ; bx = b1 ; }
	      if (pass && type == 3)
		{ type = 6 ; pass = 0 ;
		}
	    }
	  if (!type)
	    {
	      if (!pass)
		type = bt ;
	      else
		type = ERREUR ; 
	    }
/* erreur car rien n'a marche... mais ambigue testes avant */
	  newErr = arrayp(errArray, err++, A_ERR) ;
	  newErr->type = type ;
	  newErr->iLong = cp - cp0 ;
	  newErr->iShort = cq - cq0 ;
	  newErr->sens = sens ;
	  switch(type)
	    {
	    case TROU_DOUBLE:
	      newErr->baseLong = *cp ;
	      newErr->baseShort = '#' ;
	      cp += 2 ;
	      i -= 2 ;
	      break ;
	    case TROU:
	      newErr->baseLong = *cp ;
	      newErr->baseShort = '*' ;
	      cp++ ;
	      i-- ;
	      break ;
	    case INSERTION_DOUBLE:
	      newErr->baseLong = *cp ;
	      newErr->baseShort = *cq ;
	      cq += 2 * sens ;
	      j -= 2 ;
	      break ;
	    case INSERTION:
	      newErr->baseLong = *cp ;
	      newErr->baseShort = *cq ;
	      cq += sens ;
	      j-- ;
	      break ;
	    case ERREUR:
	      newErr->baseLong = *cp ;
	      newErr->baseShort = *cq ;
	      break ;
	    }
	}
      cp++ ;
      cq += sens ;
    }
  cp-- ; /* on sort un au dessus du dernier teste */
  cq -= sens ;
  *recouv += cp - cp0 ;
  *stop = cp - cp0 + j ; /* si short depasse */
  if (all && j) /* ambigue en dehors de longDna */
    { cp-- ; cq -= sens ;
      while (cp++, cq += sens, j--)
	if (!*cq || *cq == N_)
	  { n++ ;
	    newErr = arrayp(errArray, err++, A_ERR) ;
	    newErr->type = AMBIGUE ;
	    newErr->iLong = cp - cp0 ;
	    newErr->iShort = cq - cq0 ;
	    newErr->baseLong = *cp ;
	    newErr->baseShort = *cq ;
	    newErr->sens = sens ;
	  }
    }
  if (all) *nbN = n ; /*renvoie le nb d'ambiguite si on demande le tableau complet */
}

/***************************************************************/

Array alignToolsMatch(Array longDna, int xmin, int xmax, int x1, int x2,
		      Array shortDna, int ymin, int ymax, int y1, int y2,
		      int *debut, int *end, int *sens, int all, int taux, int nbessai,
		      BOOL control) /* control pour rester dans la zone */
{ Array matching = 0, erreur = 0, erreur2 = arrayCreate (50, A_ERR), atemp = 0 ;
  MATCH *test ;
  unsigned int oligo ;
  int deb, fin, i, start, stop, over = 1, nberr = 1, nberr2, over2, nbN ;
  int zonedeb = x1 - 10, zonefin = x2 + 10, mv ; /* 10 est un kludge pour pouvoir depasser */
  BOOL biger ;

  if (x2 - x1 >= y2 - y1)
    { biger = TRUE ;
      deb = y1 ;
      fin = y2 ;
    }
  else 
    { biger = FALSE ;
      deb = x1 ;
      fin = x2 ;
      control = FALSE ; /* kludge to realign biger a corriger dans les appels */
    }
  mv = (y2 - y1)/(2 * nbessai + 1) ;
  while(nbessai--)
    { if (deb >= fin - 11)
	break ;
      matching = arrayReCreate(matching, 10, MATCH) ; /* 12 car dodecameres (a changer ?) */
      if (biger)
	{ if (alignToolsMakeOligo(shortDna, deb, fin, 12, &oligo, &deb))
	    alignToolsFindOligo(longDna, x1, x2, matching, deb, oligo) ;
	  if (alignToolsMakeOligo(shortDna, fin, deb, 12, &oligo, &fin))
	    alignToolsFindOligo(longDna, x1, x2, matching, fin, oligo) ;
	  if (!arrayMax(matching))
	    { deb += 4 + mv ; fin += 7 - mv ;
/*	      deb += 7 ; fin += 4 ; */
/*	      deb++ ; fin += 10 ; */
	      continue ;
	    }
	}
      else
	{ if (alignToolsMakeOligo(longDna, deb, fin, 12, &oligo, &deb))
	    alignToolsFindOligo(shortDna, y1, y2, matching, deb, oligo) ;
	  if (alignToolsMakeOligo(longDna, fin, deb, 12, &oligo, &fin))
	    alignToolsFindOligo(shortDna, y1, y2, matching, fin, oligo) ;
	  if (!(i = arrayMax(matching)))
	    { deb += 4 + mv ; fin += 7 - mv ;
/*	      deb += 7 ; fin += 4 ; */
/*	      deb++ ; fin += 10 ; */
	      continue ;
	    }
	  test = arrp(matching, 0, MATCH) - 1 ;
	  while(test++, i--)
	    { over2 = test->pos1 ;
	      test->pos1 = test->pos2 ;
	      test->pos2 = over2 ;
	    }
	}
      i = arrayMax(matching) ;
      test = arrp(matching, 0, MATCH) - 1 ;
      while(test++, i--)
	{ nbN = all ;
	  localCptErreur(longDna, xmin, xmax, test->pos1, shortDna, ymin, ymax,
			 test->pos2, test->sens, &nbN, &start, &stop, &over2, erreur2) ;
	  nberr2 = arrayMax(erreur2) - nbN ;
	  if (control && (start < zonedeb || stop > zonefin))
	    continue ;
	  if ((100 * nberr2 <= taux * over2) && (over2 >= OVERLAP_MIN) &&
	      (nberr2 * over < nberr * over2))
	    { atemp = erreur ;
	      erreur = erreur2 ;
	      if (atemp)
		erreur2 = atemp ;
	      else
		erreur2 = arrayCreate (50, A_ERR) ;
	      over = over2 ;
	      nberr = nberr2 ;
	      *debut = start ;
	      *end = stop ;
	      *sens = test->sens ;
	    }
	}
      if (erreur)
	break ;
/*    deb++ ; fin += 10 ;  */
      deb += 4 + taux ; fin += 7 - taux ;
    }
  arrayDestroy (erreur2) ;
  arrayDestroy (matching) ;
  return erreur ;
}

/***************************************************************/
/************************ Fix Consensus ************************/
/***************************************************************/

static int allErrOrder(const void *a, const void *b)
{
  int i = ((A_ERR*) a)->iLong - ((A_ERR*) b)->iLong ;
  if (i)
    return i ;
  i = (int)(((A_ERR*) a)->type - ((A_ERR*) b)->type) ; 
  if (i)
    return i ;
  i = (int)(((A_ERR*) a)->baseShort - ((A_ERR*) b)->baseShort) ; 
  return i ;
}

/***************************************************************/

static void doFindAssemble (OBJ Whole, DEFCPT look, FIXTOOL tool, KEY seg,
			    int a, int b, int ctop, KEY fath, int let)
{ KEY part ;
  int x, y, nn, i, top = 0 ;
  OBJ obj  = bsCreate(seg) ;
  BSMARK bsk = 0 ;
  KEY dna ;
  Array temp ;

  if (!obj)
    return ;
  bsk = bsMark (obj, 0) ; /* allocates bsk */
  if (bsGetKey (obj, _Assembled_from, &part))
    { if (let)
	{ if (let == 1)
	    { bsMark(obj, bsk) ;
	      if (bsGetKey(obj, _DNA, &dna) && 
		  bsGetData(obj, _bsRight, _Int, &nn))
		{ i = arrayMax(tool->dummy) ;
		  array(tool->dummy, i++, BSunit).k = seg ;
		  if (a < b)
		    { x = a + 1 ; y = a + nn ; } /* bio algebra */
		  else
		    { x = a + 1 ; y = a + 2 - nn ; } /* bio algebra */
		  array(tool->dummy, i++, BSunit).i = x ;
		  array(tool->dummy, i, BSunit).i = y ;
		}
	      bsGoto(obj, bsk) ;
	    }
	  let-- ;
	}
      do
	{ bsMark (obj, bsk) ;
	  if (bsGetData (obj, _bsRight, _Int, &x) &&
	      bsGetData (obj, _bsRight, _Int, &y))
	    { x-- ; y-- ; /* info algebra */
	      if (!bsGetData (obj, _bsRight, _Int, &top))
		top = 0 ;
	      if (a < b)
		doFindAssemble (Whole, look, tool, part, a + x, a + y, top, seg, let) ;
	      else
		doFindAssemble (Whole, look, tool, part, a - x, a - y, top, seg, let) ;
	    }
	  bsGoto (obj, bsk) ;
	} while (bsGetKey (obj, _bsDown, &part)) ; /* do while construct */
      if (seg != tool->target && bsGetKey (obj, _DNA, &dna))
	monDnaForget (look, dna) ;
      if (fath == tool->target && bsFindKey (Whole, _Assembled_from, seg))
	bsRemove (Whole) ;
    }
  else if (seg != tool->target &&
	   bsGetKey (obj, _DNA, &dna) && 
	   (temp = monDnaGet (look, fath, dna)) && (nn = arrayMax(temp)))
    { 
      if (!ctop)
	{ if (!bsGetData (obj, _Clipping, _Int, &top))
	    top = 1 ;
	}
      else 
	top = ctop ;
      bsAddKey (Whole, _Assembled_from, seg) ;
      if (a < b)
	{ x = a + 1 ; y = a + nn ; } /* bio algebra */
      else
	{ x = a + 1 ; y = a + 2 - nn ; } /* bio algebra */
      bsAddData (Whole, _bsRight, _Int, &x) ;
      bsAddData (Whole, _bsRight, _Int, &y) ;
      bsAddData (Whole, _bsRight, _Int, &top) ;
      top += nn - 1 ;
      bsAddData (Whole, _bsRight, _Int, &top) ;
    }
  bsDestroy (obj) ;
  bsMarkFree (bsk) ;
}

/***************************************************************/

static void toolCleanErr (Array err, Array lgDna, Array shDna, BOOL isUp)
{ int i, n = arrayMax(err), start, stop ;
  A_ERR *eq, *ep, *epMax ;
  char *cs, cc ;
  int sens = isUp ? -1 : 1 ;

  if (!n) return ;
  if (isUp)
    { ep = arrp (err, 0, A_ERR) - 1 ;
      n-- ;
      if (n>0) while (ep++, n--)
	if ((ep->type == INSERTION) && ((ep+1)->type == ERREUR) &&
	    (ep->iShort == (ep+1)->iShort + 1))
	  { i = ep->type ; ep->type = (ep+1)->type ; (ep+1)->type = i ;
	    i = ep->iLong ; ep->iLong = (ep+1)->iLong ; (ep+1)->iLong = i + 1;
	  }
    }
  n = arrayMax (err) ;
  ep = arrp(err, 0, A_ERR) - 1 ;
  while (ep++, n--)
    { 
      switch (ep->type)
	{
	case TROU: case TROU_DOUBLE:
	  if (ep->iLong < arrayMax(lgDna))
	    { cs = arrp(lgDna, ep->iLong, char) ;
	      cc = *cs ; i = arrayMax(lgDna) - ep->iLong ;
	      while (i-- && *(++cs) == cc)
		{ ep->iShort += sens ;
		  ep->iLong++ ;
		}
	      if (ep->type == TROU_DOUBLE && !isUp)
		ep->iShort -- ;
	    }
	  break ;
	case INSERTION: case INSERTION_DOUBLE:
	  start = ep->iShort ;
	  epMax = arrp (err, arrayMax (err) - 1, A_ERR) + 1 ;
	  if (!isUp && ep->iShort < arrayMax (shDna))
	    { cs = arrp (shDna, ep->iShort, char) ;
	      cc = *cs ; i = arrayMax(shDna) - ep->iShort ;
	      while (i-- && *(++cs) == cc)
		{ ep->iShort++ ;
		  ep->iLong++ ;
		}
	      if (ep->type == INSERTION_DOUBLE)
		ep->iLong-- ;
	    }
	  else if (isUp && ep->iShort < arrayMax(shDna))
	    { cs = arrp (shDna, ep->iShort, char) ;
	      cc = *cs ; i = ep->iShort ;
	      while (i-- && *(--cs) == cc && ep->iShort > 0)
		{ ep->iShort-- ;
		  ep->iLong++ ;
		}
	    }
	  stop = ep->iShort ;
	  eq = ep ;
	  if (isUp)
	    while (++eq < epMax && eq->iShort < start && 
		   eq->iShort >= stop)
	      eq->iShort++ ;
	  else
	    while (++eq < epMax && eq->iShort > start && 
		   eq->iShort <= stop)
	      eq->iShort-- ;
	  break ;
	default: /* point error or N */
	  break ;
	}
    }
}

/***************************************************************/
#define KLUDGE
static Array adjustAssembleSeg (FIXTOOL tool, Array shDna, 
				int *u1, int *u2)
{ int i, v1, v2, v3, sens, end = arrayMax(shDna) - 1, taux, x1, x2 ;
  Array err = 0, longDna ;
  BOOL isUp = FALSE ;
  A_ERR *ep, *eq, er ;
  int FENETRE = 30 ;
  int iniu1 = *u1, iniu2 = *u2 ;

  v3 = arrayMax(tool->dnad) - 1 ;
  if (*u1 > *u2)
    { isUp = TRUE ;
      *u1 = v3 - *u1 ;
      *u2 = v3 - *u2 ;
      longDna = tool->dnar ;
    }
  else
    longDna = tool->dnad ;
  v1 = *u1 > FENETRE ? *u1 - FENETRE : 0 ;
  v2 = *u2 < v3 - FENETRE ? *u2 + FENETRE : v3 ;
  taux = 4 ;
  while (taux < 20 && !(err = alignToolsMatch (longDna, v1, v2, v1, v2, shDna, 0, end, 0, end,
					       u1, u2, &sens, 1, taux, 5, TRUE)) && sens != 1)
    { arrayDestroy (err) ;
      taux *= 4 ;
    }
#ifdef KLUDGE
  if (!err)
    err = dnaAlignCompareDna (longDna, shDna, &x1, &x2, &sens, FALSE) ;
  if (!err)
    { *u1 = iniu1 ; /* no realignment -> return with no changes */
      *u2 = iniu2 ;
      return 0 ;
    }
#else
  if (!err)
    if (!dnaAlignAsmbPaireDna (longDna, shDna, 7, 500, -1, 
			       5, &x1, &x2, &sens, FALSE) ||
	x2 < v1 || x2 > v2)
      { /* invokeDebugger() ;  */
	*u1 = iniu1 ; /* no realignment -> return with no changes */
	*u2 = iniu2 ;
	return 0 ;
      }
    else
      { err = arrayCreate (30, A_ERR) ;
	localCptErreur(longDna, v1, v2, x1,
		       shDna, 0, end, x2,
		       sens, &nbN, u1, u2, &recouv, err) ;
      }
#endif
  if (isUp)
    { *u1 = v3 - *u1 ;
      *u2 = v3 - *u2 ;
      i = arrayMax (err) ;
      if (i)
	{ ep = arrp (err, 0, A_ERR) - 1 ;
	  while (ep++, i--)
	    { switch (ep->type)
		{
		case INSERTION: case INSERTION_DOUBLE:
		  ep->iLong-- ;
		  break ;
		case TROU: case TROU_DOUBLE:
		  ep->iShort-- ;
		  break ;
		default: /* point error or N */
		  break ;
		}
	      ep->iLong = v3 - ep->iLong ;
	      ep->sens = - 1 ;
	    }
	  i = arrayMax (err) ;
	  ep = arrp (err, 0, A_ERR) ;
	  eq = arrp (err, i - 1, A_ERR) ;
	  while (ep < eq)
	    { er = *ep ; *ep++ = *eq ; *eq-- = er ; }
	}
    }
  toolCleanErr (err, tool->dnad, shDna, isUp) ;
  return err ;
}

/***************************************************************/

static void adjustAssembly (DEFCPT look, FIXTOOL tool)
{ int x1, x2, z1, z2, i, j, k, *ip, min, max ;
  BSunit *u ;
  KEY seg, dnaKey ;
  OBJ obj2 ;
  Array shDna, err = 0 ;
  A_ERR *ep, *eq ;

  tool->err = arrayCreate (256, A_ERR) ;
  min = tool->min ;
  max = tool->max ;
  for (k = 0 ; k < arrayMax (tool->seq) ; k += 4)
    { u = arrp (tool->seq, k, BSunit) ;
      seg = u[0].k ;
      x1 = u[1].i - 1 ; /* info */
      x2 = u[2].i - 1 ;
      if ((x1 < min && x2 < min) || (x1 > max && x2 > max))
	continue ;
      if ((obj2 = bsCreate (seg)) &&
	  bsGetKey (obj2, _DNA, &dnaKey) &&
	  (shDna = monDnaGet (look, tool->target, dnaKey)))
	{ err = adjustAssembleSeg (tool, shDna, &x1, &x2) ;
	  if (err && arrayMax(err))
	    { i = arrayMax(err) ;
	      ep = arrp(err, i - 1, A_ERR) ;
	      j = arrayMax(tool->err) + i ;
	      eq = arrayp(tool->err, j - 1, A_ERR) ;
	      while(i--)
		*eq-- = *ep-- ;
	    }
	  if (err)
	    { 
	      if (x1 < x2)
		{ z1 = x1 >= 0 ? x1 : 0 ;
		  z2 = x2 < arrayMax(tool->nbt) ? x2 : arrayMax (tool->nbt) - 1 ;
		}
	      else
		{ z1 = x2 >= 0 ? x2 : 0 ;
		  z2 = x1 < arrayMax(tool->nbt) ? x1 : arrayMax (tool->nbt) - 1 ;
		}
	      for (i = z1, ip = arrp (tool->nbt, i, int) ; i <= z2 ; i++, ip++)
		(*ip)++ ;
	      x1++ ; x2++ ; /* no zero */
	      u[1].i = x1 ;
	      u[2].i = x2 ;
	    }
	  arrayDestroy(err) ;
	}
      bsDestroy (obj2) ;
    }
}

/***************************************************************/
/* x1 incluse, x2 exclue y1 vrai max de new (derniere base + 1) */
static void toolShortCopy (Array old, Array new, Array mod, 
			   int x1, int x2, int *y1)
{ char *cp, *cq, *cq0 ;
  int *ip, i ;

  i = x2 - x1 ;
  if (i <= 0) return ;
  array (new, *y1 + i - 1, char) = 0 ;
  ip = arrp (mod, x1, int) ;
  cp = arrp (old, x1, char) ;
  cq0 = arrp (new, 0, char) ;
  cq = cq0 + *y1 ;
  *y1 += i ;
  if (i>0) while (i--)
    { *ip++ = cq - cq0 ;
      *cq++ = *cp++ ;
    }
}

/***************************************************************/

static void fixDnaConsensus (FIXTOOL tool, BOOL first)
{ int i, j, max, safemax, maxerr = arrayMax (tool->err), n, nt, ntbis ;
  int *ip, minc = arrayMax (tool->dnad), maxc = - 1, maxdna = 0 ;
  Array dna2 = dnaCopy (tool->dnad) ;
  char *safe, *cp, *cp0, *cq, *cq0, cr=0 ;
  BOOL fixtest = FALSE ;
  int errBase[9], errType[9] ;
  A_ERR *ep, *eq, *epMax ;
  int sens ;

  tool->isFixed = FALSE ;
  if (!maxerr)
    return ;
  if (tool->min > tool->max)
    return ;
  max = arrayMax (tool->dnad) ;
  tool->mod = arrayCreate (max, int) ;
  array (tool->mod, max - 1, int) = 0 ;
  safemax = max + 200 ;
  safe = arrayp (tool->dnad, safemax, char) ;  /* make room */
  if (tool->min)
    toolShortCopy (dna2, tool->dnad, tool->mod, 0, tool->min, &maxdna) ;
  cp = arrp (tool->dnad, maxdna, char) ;
  cq = arrp (dna2, tool->min, char) ;
  cp0 = arrp (tool->dnad, 0, char) ;
  cq0 = arrp (dna2, 0, char) ;
  if (!cq0) messcrash("") ;
  ip = arrp (tool->mod, tool->min, int) ;
  ep = arrp (tool->err, 0, A_ERR) ;
  epMax = arrp (tool->err, maxerr - 1, A_ERR) + 1 ;
  while (ep < epMax && ep->iLong < tool->min) ep++ ;
  n = ep < epMax ? ep->iLong : -1 ;
  for (i = tool->min ; i < tool->max ; i++)
    { if (i != n)
        { *ip++ = cp - cp0 ;
	  *cp++ = *cq++ ;
	  continue ;
	}
      if (cq - cq0 > i + 1000 || cq < cq0 - 100)
	invokeDebugger () ;
      j = 9 ;
      fixtest = TRUE ;
      while(j--)
	errBase[j] = errType[j] = 0 ;
      eq = ep ;
      while(eq < epMax && eq->iLong == n)
	{ errType[eq->type]++ ;
	  if (eq->type == ERREUR)
	    { if (eq->sens == 1)
		errBase[(int)eq->baseShort]++ ;
	      else
		errBase[(int)complementBase(eq->baseShort)]++ ;
	    }
	  eq++ ;
	}
      nt = arr(tool->nbt, i, int) ;
      ntbis = nt ; /* pour eviter de choisir les ambigues */
      if (first && nt > 1)
	nt-- ; /* kludge pour accepter la majorite relative au premier tour */
      j = errType[TROU] + errType[TROU_DOUBLE] ;
      if (errType[INSERTION] + errType[INSERTION_DOUBLE] > nt / 2)
	{ 
	  if (errType[INSERTION_DOUBLE] > nt / 2)
	    while(ep < epMax && ep->type != INSERTION_DOUBLE) ep++ ;
	  else
	    while(ep < epMax && ep->type != INSERTION) ep++ ;
	}
      else if (j > nt / 2)
	{ 
	  if (errType[TROU_DOUBLE] > nt / 2)
	    while(ep < epMax && ep->type != TROU_DOUBLE) ep++ ;
	  else
	    while(ep < epMax && ep->type != TROU) ep++ ;
	}
      else if (errType[ERREUR] > (nt - j - errType[AMBIGUE]) / 2)
	{ nt -= j + errType[AMBIGUE] ;
	  while(ep < epMax && ep->type != ERREUR) ep++ ;
	  if (*cq && *cq != N_)
	    errBase[(int)*cq] = nt - errBase[A_] - errBase[T_] - errBase[G_] - errBase[C_] ;
	  if ( *cq && *cq != N_ && errBase[(int)*cq] >= errBase[A_] && 
	      errBase[(int)*cq] >= errBase[T_] && errBase[(int)*cq] >= errBase[G_] && 
	      errBase[(int)*cq] >= errBase[C_])
	    fixtest = FALSE ;
	  else
	    { if (errBase[A_] >= errBase[T_] && 
		  errBase[A_] >= errBase[G_] && errBase[A_] >= errBase[C_])
		cr = (char) A_ ; /* cast a cause de la dec */
	      if (errBase[T_] >= errBase[A_] && 
		  errBase[T_] >= errBase[G_] && errBase[T_] >= errBase[C_])
		cr = (char) T_ ;
	      if (errBase[G_] >= errBase[A_] && errBase[G_] >= errBase[T_] && 
		  errBase[G_] >= errBase[C_])
		cr = (char) G_ ;
	      if (errBase[C_] >= errBase[A_] && errBase[C_] >= errBase[T_] && 
		  errBase[C_] >= errBase[G_])
		cr = (char) C_ ;
	    }
	}
      else if (errType[AMBIGUE] == ntbis - j) /* que des n kludge post edition */
	{ while(ep < epMax && ep->type != AMBIGUE) ep++ ;
	  if (!*cq || *cq == N_) /* correction deja faite */
	    fixtest = FALSE ;
	}
      else fixtest = FALSE ; /* pas de correction a faire */
      if (fixtest)
	{ sens = ep->sens ;
	  tool->isFixed = TRUE ;
	  if (minc > cp - cp0)
	    minc = cp - cp0 ;
	  if (maxc < cp - cp0)
	    maxc = cp - cp0 ;
	  *ip++ = cp - cp0 ;
	  switch (ep->type)
	    {
	    case ERREUR:
	      *cp++ = cr ;
	      cq++ ;
	      break ;
	    case TROU:
	      cq++ ;
	      break ;   /* insert nothing in consensus */
	    case TROU_DOUBLE:
	      *ip++ = cp - cp0 ;
	      cq += 2 ;
	      while (ep < epMax && ep->iLong == i) ep++ ;
	      i++ ;    /* sinon on decale i par rapport a cq - cq0 */
	      break ;  /* insert nothing in consensus */
	    case INSERTION:
	      if (sens == 1)
		*cp++ = ep->baseShort ;
	      else
		*cp++ = complementBase(ep->baseShort) ;
	      *cp++ = *cq++ ;
	      break ;
	    case INSERTION_DOUBLE:
	      if (sens == 1)
		*cp++ = ep->baseShort ;
	      else
		*cp++ = complementBase(ep->baseShort) ;
	      *cp++ = (char)A_ ; /* juste pour retablir le bon nombre de base */
	      *cp++ = *cq++ ;
	      break ;
	    case AMBIGUE: /* kludge post edition de sequence */
	      *cp++ = 0 ;
	      cq++ ;
	      break ;
	    default:
	      *cp++ = *cq++ ;
	      break ;
	    }
	}
      else
	{ *ip++ = cp - cp0 ;
	  *cp++ = *cq++ ; /* pas de correction */
	}
   /* because of the insertions, we may run out of dna */
      if (safe - cp < max - i + 2)
	{ int ii = cp - cp0 ;
	  safe = arrayp (tool->dnad, safemax += 100 , char) ;  /* make room */
	  cp0 = arrp (tool->dnad, 0, char) ;
	  cp = cp0 + ii ;
	}
      while(ep < epMax && ep->iLong == i)
	ep++ ;
      n = ep < epMax ? ep->iLong : -1 ;
    }
  maxdna = cp - cp0 ;
  if (tool->max != cq - cq0)
    { /* invokeDebugger () ;*/
      tool->max = cq - cq0 ;
      if (tool->max > max)
	tool->max = max ;
    }
  if (tool->max < max)
    toolShortCopy (dna2, tool->dnad, tool->mod, tool->max, max, &maxdna) ;
  arrayMax(tool->dnad) = maxdna ;
  tool->min = minc > tool->min + 25 ? minc - 25 : tool->min ; /* inclu */
  tool->max = maxc < maxdna - 25 ? maxc + 25 : maxdna ; /* non inclu */
  arrayDestroy(dna2) ;
}

/***************************************************************/
/* profondeur du tableau 3 ou 4 */
static void alignToolsAddPreviousSeq (DEFCPT look, FIXTOOL tool, Array preSeq, 
				      KEY tag, int prof)
{ int i, imax = arrayMax(preSeq), x, y, z1 = 0, z2 ;
  OBJ Whole = bsUpdate(tool->target) ;
  KEY seg, segDna ;
  BSunit *unit ;
  Array shDna, err ;

  if (!imax || !Whole)
    { bsDestroy (Whole) ;
      return ;
    }
  unit = arrp(preSeq, 0, BSunit) ;
  for (i = 0 ; i < imax ; i += prof)
    { seg = (unit++)->k ;
      x = (unit++)->i - 1 ; /* info algebra */
      y = (unit++)->i - 1 ; /* info algebra */
      if (prof == 4)
	z1 = (unit++)->i ;
      if (tool->target == seg)
	{ if (bsFindKey(Whole, _Assembled_from, seg))
	    bsRemove(Whole) ;
	  continue ;
	}
      dnaSubClass (seg, &segDna) ;
      shDna = monDnaGet (look, tool->target, segDna) ;
      if (shDna)
	{ err = adjustAssembleSeg(tool, shDna, &x, &y) ;
	  arrayDestroy(err) ;
	  x++ ; y++ ; /* passage en coordonnees bio */
	  bsAddKey (Whole, tag, seg) ;
	  bsAddData (Whole, _bsRight, _Int, &x) ;
	  bsAddData (Whole, _bsRight, _Int, &y) ;
	  if (z1)
	    { bsAddData (Whole, _bsRight, _Int, &z1) ;
	      z2 = z1 + arrayMax (shDna) - 1 ;
	      bsAddData (Whole, _bsRight, _Int, &z2) ;
	      z1 = 0 ;
	    }
	}
    }
  bsSave(Whole) ;
}

/***************************************************************/

static void fixToolDestroy (FIXTOOL tool)
{ if (arrayExists (tool->seq))
    arrayDestroy (tool->seq) ;
/*if (arrayExists (tool->dnad)) // is suppose to be in the associator for monDnaGet
    arrayDestroy (tool->dnad) ; */
  if (arrayExists (tool->dnar))
    arrayDestroy (tool->dnar) ;
  if (arrayExists (tool->err))
    arrayDestroy (tool->err) ;
  if (arrayExists (tool->nbt))
    arrayDestroy (tool->nbt) ;
  if (arrayExists (tool->mod))
    arrayDestroy (tool->mod) ;
  if (arrayExists (tool->dummy))
    arrayDestroy (tool->dummy) ;
  if (arrayExists (tool->alig))
    arrayDestroy (tool->alig) ;
  messfree (tool) ;
}

/***************************************************************/

static void alignToolsAdjustPos (FIXTOOL tool)
{ int i, j, max ;
  BSunit *u ;

  max = arrayMax (tool->mod) ;
  for (i = 0 ; i < arrayMax (tool->seq) ; i += 4)
    { u = arrp (tool->seq, i, BSunit) ;
      j = u[1].i ; /* bio algebra */
      if (j > max) j = max ; /* kludge pour Jean */
      if (j < 1) j = 1 ;
      u[1].i = arr (tool->mod, j - 1, int) + 1 ; /* read info + -> bio */
      j = u[2].i ;
      if (j > max) j = max ;
      if (j < 1) j = 1 ;
      u[2].i = arr (tool->mod, j - 1, int) + 1 ;
    }
}

/***************************************************************/

void doForceAssembleSeg (DEFCPT look, KEY seg, KEY target, int min, int max,
			int let, int cut)
{ OBJ obj = 0, Whole = 0 ;
  KEY dnaKey1, dnaKey2 ;
  int i, j ;
  Array dna, temp = 0 ;
  BOOL first = TRUE ;
  FIXTOOL tool ;
  BSunit *u ;

  messStatus ("Fixing") ;
  nbcpterror = 0 ;
  tool = (FIXTOOL)messalloc(sizeof(struct FIXTOOLSTUFF)) ;
  tool->magic = FIXTOOLMAG ;
  tool->seg = seg ;
  tool->target = target ;
  tool->min = min ;
  tool->max = max ;
  if ((obj = bsCreate (seg)))
    { if (bsGetKey (obj, _DNA, &dnaKey1))
	{ 
	  if (seg != target)
	    { if ((dna = dnaGet(dnaKey1)))
		{ lexaddkey (name(target), &dnaKey2, _VDNA) ;
		  dnaStoreDestroy (dnaKey2, dna) ;
		  tool->dnakey = dnaKey2 ;
		}
	    }
	  else tool->dnakey = dnaKey1 ;
	}
      bsDestroy (obj) ;
    }
  if (let)
    tool->dummy = arrayCreate (10, BSunit) ;
  if ((Whole = bsUpdate (target)))
    { doFindAssemble (Whole, look, tool, seg, 0, 1, 0, seg, let) ;
      bsSave (Whole) ;
      alignToolsAdjustContigSize (look->link, target) ;
      tool->seq = arrayCreate (1000, BSunit) ;
      if (!(Whole = bsCreate (target)) ||
	  !bsFindTag (Whole, _Assembled_from) || 
	  !bsFlatten (Whole, 4, tool->seq))
	goto abort ;
      bsDestroy (Whole) ;

      tool->dnad = monDnaGet (look, tool->target, tool->dnakey) ;
      if (!tool->dnad) goto abort ;
      tool->dnar = dnaCopy (tool->dnad) ;
      reverseComplement (tool->dnar) ;

      alignToolsAddPreviousSeq (look, tool, tool->seq, _Assembled_from, 4) ;
      alignToolsAdjustContigSize (look->link, target) ;
      arrayDestroy (tool->seq) ;
      arrayDestroy (tool->dnar) ;
    }
  else goto abort ;
  if ((obj = bsCreate (seg)))
    { if (bsFindTag (obj, _Previous_contig))
	{ temp = arrayCreate (10, BSunit) ;
	  bsFlatten (obj, 3, temp) ;
	  if (!tool->dummy)
	    tool->dummy = temp ;
	  else
	    { i = arrayMax (tool->dummy) ;
	      for (j = 0 ; j < arrayMax (temp) ; j += 3)
		{ u = arrp (temp, j, BSunit) ;
		  array (tool->dummy, i++, BSunit).k = u[0].k ;
		  array (tool->dummy, i++, BSunit).i = u[1].i ;
		  array (tool->dummy, i++, BSunit).i = u[2].i ;
		}
	      arrayDestroy (temp); 
	    }
	}
      if (bsFindTag (obj, _Aligned))
	{ tool->alig = arrayCreate (10, BSunit) ;
	  bsFlatten (obj, 4, tool->alig) ;
	}
      bsDestroy (obj) ;
    }
  if (!(Whole = bsCreate (target))) goto abort ;
  tool->seq = arrayCreate (1000, BSunit) ;
  if (!bsFindTag (Whole, _Assembled_from) || 
      !bsFlatten (Whole, 4, tool->seq))
    goto abort ;
  bsDestroy (Whole) ;

  tool->dnad = monDnaGet (look, tool->target, tool->dnakey) ;
  arrayDestroy (tool->dnar) ;
  if (!tool->dnad)		/* may happen if we canot find oligos in the contig */
    goto abort ;

#ifdef JUNK
  trackFixContig(look->link, target, dnaKey2, tool->seq, tool->dnad) ;
  dnaStoreDestroy (tool->dnakey, dnaCopy (tool->dnad)) ;
#else
/*************  fix ulrich  ************/

  if (tool->max == -2) /* mean we want fix all */
    tool->max = arrayMax (tool->dnad) ;

  i = 6 ; /* NB MAX D'ITERATIONS */
  tool->isFixed = TRUE ;
  while (i-- && tool->isFixed)
    { /* messalloccheck () ; */
      j = arrayMax (tool->dnad) ;
      tool->nbt = arrayCreate (j, int) ;
      array (tool->nbt, j - 1, int) = 0 ;
      arrayDestroy (tool->dnar) ;
      tool->dnar = dnaCopy (tool->dnad) ;
      reverseComplement (tool->dnar) ;
      adjustAssembly (look, tool) ;
     /*  messalloccheck () ; */
      if (tool->err && arrayMax(tool->err))
	{ arraySort (tool->err, allErrOrder) ;
	    messStatus (messprintf ("Fixing %s", name(tool->target))) ;
	  fixDnaConsensus (tool, first) ;
	  first = FALSE ;
	  arrayDestroy (tool->err) ;
	  alignToolsAdjustPos (tool) ;
	  if (tool->max <= tool->min)
	    tool->isFixed = FALSE ;
	}
      else tool->isFixed = FALSE ;
      arrayDestroy (tool->nbt) ;
    }
  dnaStoreDestroy (tool->dnakey, dnaCopy (tool->dnad)) ;


/*********  fin du fix ******************/
#endif

/*
  printf ("nombre d'iteration %d (%s)\n", 5 - i, name (tool->target)) ;
  printf ("Appel localCptError %d (nb seq : %d) ; total %d\n", nbcpterror,
	  arrayMax (tool->seq)/4, localCptReport) ;
*/
  arrayDestroy (tool->dnar) ;
  tool->dnar = dnaCopy (tool->dnad) ;
  reverseComplement (tool->dnar) ;

  if (tool->seq)
    alignToolsAddPreviousSeq(look, tool, tool->seq, _Assembled_from, 4) ;
  if (tool->dummy)
    alignToolsAddPreviousSeq(look, tool, tool->dummy, _Previous_contig, 3) ;
  if (tool->alig)
    alignToolsAddPreviousSeq(look, tool, tool->alig, _Aligned, 4) ;
  alignToolsAdjustContigSize (look->link, target) ;
 abort:
  fixToolDestroy (tool) ;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/

Array alignToolsCptErreur(Array longDna, Array shortDna,
			  int *debutp, int *finp, int *shdebutp, int *shfinp)
{ int debut = *debutp, fin = *finp, shdeb, shfin, sens, recou, nn = 1 ;
  Array erreur = arrayCreate (50, A_ERR) ;

  if (debut > fin)
    messerror ("bad sign for longDna") ;
  if (*shdebutp < *shfinp)
    {
      sens = 1 ;
      shdeb = *shdebutp ;
      shfin = *shfinp ;
    }
  else
    {
      sens = -1 ;
      shdeb = *shfinp ;
      shfin = *shdebutp ;
    }
  localCptErreur (longDna, debut, fin, debut, shortDna, shdeb, shfin, 
		  *shdebutp, sens, &nn, debutp, finp, &recou, erreur) ;
  return erreur ;
}

/***************************************************************/

Array alignToolsMakeErreur(Array longDna, int *debutp, int *finp, int pol,
			   Array shortDna, int shdebut, int shfin, int pos,
			   int sens)
{ Array erreur = arrayCreate (50, A_ERR) ;
  int nn = 1, recou ;
  
  localCptErreur(longDna, *debutp, *finp, pol, shortDna, shdebut,
		 shfin, pos, sens, &nn, debutp, finp, &recou, erreur) ;
  if (sens == - 1) /* inversion des coordonnees pour indiquer reverse */
    { nn = *debutp ;
      *debutp = *finp ;
      *finp = nn ;
    }
  return erreur ;
}

/***************************************************************/
/* fine adjust : return FALSE if no change done in the object */
BOOL alignToolsAdjustContigSize (KEY link, KEY key)
{ KEY dnaKey = 0, keymin = 0, keymax = 0 ;
  Array seq = 0, dna = 0, newdna = 0, subdna = 0 ;
  int i, j, min, max, maxdna, newmax, x1, x2, x3 ;
  BSunit *u ;
  OBJ obj = 0, newobj = 0 ;

  char *cp=0, *cq=0 ;
  DEFCPT look = defCptGetLook (link) ;

  if (class (key) == _VDNA)
    { dnaKey = key ;
      dnaReClass (key, &key) ;
    }
  else if (!dnaReClass (key, &key))
    return FALSE ;
  if (!(obj = bsCreate (key)))
    return FALSE ;
  seq = arrayCreate (2000, BSunit) ;
  if (!bsGetKey (obj, _DNA, &dnaKey) || !(dna = monDnaGet (look, 0, dnaKey)) ||
      !bsFindTag (obj, _Assembled_from) || !bsFlatten (obj, 5, seq))
    goto abort ;
  maxdna = arrayMax (dna) ;
  max = -10000000 ; min = 10000000 ;
  u = arrp (seq, 0, BSunit) ;
  for (i = 0 ; i < arrayMax (seq) ; i += 5, u += 5)
    { j = u[1].i ;
      if (j > max) 
	{ max = j ;
	  keymax = u[0].k ;
	}
      if (j < min) 
	{ min = j ;
	  keymin = u[0].k ;
	}
      j = u[2].i ;
      if (j > max) 
	{ max = j ;
	  keymax = u[0].k ;
	}
      if (j < min) 
	{ min = j ;
	  keymin = u[0].k ;
	}
    }
  if (max < min) goto abort ;
  min-- ; /* coordonnes info */
  if (!min && max == maxdna)
    goto abort ;

  newmax = max - min ;
  newdna = arrayCreate (newmax, char) ;
  array (newdna, newmax - 1, char) = 0 ;
  if (min >= maxdna || newmax + min <= 0) 
    j = 0 ;
  else if (min > 0)
    { cq = arrp (dna, min, char) ;
      cp = arrp (newdna, 0, char) ;
      j = maxdna - min > newmax ? newmax : maxdna - min ;
    }
  else
    { cq = arrp (dna, 0, char) ;
      cp = arrp (newdna, - min, char) ;
      j = newmax + min > maxdna ? maxdna : newmax + min ;
    }
  if (j>0) while (j--)
    *cp++ = *cq++ ;
  if (max > maxdna)
    { j = max - maxdna ;
      if (bsFindKey (obj, _Assembled_from, keymax) && bsGetData (obj, _bsRight, _Int, &x1)
	  && bsGetData (obj, _bsRight, _Int, &x2))
	{ dnaSubClass (keymax, &keymax) ;
	  subdna = monDnaGet (look, key, keymax) ;
	  if (!subdna)
	    goto laba ;
	  cp = arrp (newdna, newmax - 1, char) ;
	  if (x1 < x2)
	    { x3 = arrayMax (subdna) - 1 ;
	      cq = arrp (subdna, x3, char) ;
	      if (j>0) while (j--)
		*cp-- = *cq-- ;
	    }
	  else
	    { cq = arrp (subdna, 0, char) ;
	      if (j>0) while (j--)
			 *cp-- = complementBase(*cq++) ;
	    }
	}
    }
  if (min < 0)
    { j = - min ;
      if (keymin &&
	  bsFindKey (obj, _Assembled_from, keymin) && 
	  bsGetData (obj, _bsRight, _Int, &x1) &&
	  bsGetData (obj, _bsRight, _Int, &x2))
	{ dnaSubClass (keymin, &keymin) ;
	  subdna = monDnaGet (look, key, keymin) ;
	  if (!subdna)
	    goto laba ;
	  cp = arrp (newdna, 0, char) ;
	  if (x1 < x2)
	    { cq = arrp (subdna, 0, char) ;
	      if (j>0) while (j--)
		*cp++ = *cq++ ;
	    }
	  else
	    { x3 = arrayMax (subdna) - 1 ;
	      cq = arrp (subdna, x3, char) ;
	      if (j>0) while (j--)
			 *cp++ = complementBase(*cq--) ;
	    }
	}
    }
 laba:
  monDnaForget (look, dnaKey) ;
  dnaStoreDestroy (dnaKey, newdna) ;
  if (!(newobj = bsUpdate (key)))
    goto abort ;
  if (bsFindTag (newobj, _Assembled_from))
    bsRemove (newobj) ;
  u = arrp (seq, 0, BSunit) ;
  for (i = 0 ; i < arrayMax (seq) ; i += 5, u += 5)
    { bsAddKey (newobj, _Assembled_from, u[0].k) ;
      j = u[1].i - min ;
      bsAddData (newobj, _bsRight, _Int, &j) ;
      j = u[2].i - min ;
      bsAddData (newobj, _bsRight, _Int, &j) ;
      j = u[3].i ;
      if (j)
	{ bsAddData (newobj, _bsRight, _Int, &j) ;
	  j = u[4].i ;
	  bsAddData (newobj, _bsRight, _Int, &j) ;
	}
    }
  if (bsFindTag (obj, _Previous_contig) && bsFlatten (obj, 3, seq))
    { if (bsFindTag (newobj, _Previous_contig))
	bsRemove (newobj) ;
      u = arrp (seq, 0, BSunit) ;
      for (i = 0 ; i < arrayMax (seq) ; i += 3, u += 3)
	{ bsAddKey (newobj, _Previous_contig, u[0].k) ;
	  j = u[1].i - min ;
	  bsAddData (newobj, _bsRight, _Int, &j) ;
	  j = u[2].i - min ;
	  bsAddData (newobj, _bsRight, _Int, &j) ;
	}
    }
  if (bsFindTag (obj, _Aligned) && bsFlatten (obj, 3, seq))
    { if (bsFindTag (newobj, _Aligned))
	bsRemove (newobj) ;
      u = arrp (seq, 0, BSunit) ;
      for (i = 0 ; i < arrayMax (seq) ; i += 3, u += 3)
	{ bsAddKey (newobj, _Aligned, u[0].k) ;
	  j = u[1].i - min ;
	  bsAddData (newobj, _bsRight, _Int, &j) ;
	  j = u[2].i - min ;
	  bsAddData (newobj, _bsRight, _Int, &j) ;
	}
    }
  bsSave (newobj) ;
  arrayDestroy (seq) ;
  bsDestroy (obj) ;
  return TRUE ;
 abort:
  arrayDestroy (newdna) ;
  bsDestroy (obj) ;
  arrayDestroy (seq) ;
  return FALSE ;
}

/***************************************************************/

static void alignToolsDoPurifyContig (KEY key, KEY seqKey, KEY dnaKey, Array seq, 
				      Array dna, int a1, int a2)
{ Array dna2 ;
  int i, x1, x2, clip ;
  BSunit *u ;
  char *cp, *cq ;
  OBJ obj, fat = bsCreate (key) ;

  if (!(obj = bsUpdate (seqKey)))
    { bsDestroy (fat) ;
      return ;
    }
  i = a2 - a1 ;
  dna2 = arrayCreate (i, char) ;
  array (dna2, i - 1, char) = 0 ;
  cp = arrp (dna2, 0, char) ;
  cq = arrp (dna, a1, char) ;
  if (i>0) while (i--)
    *cp++ = *cq++ ;
  i = arrayMax(dna2) ;
  if (bsFindTag (obj, _Assembled_from))
    bsRemove (obj) ;
  for (i = 0 ; i < arrayMax (seq) ; i += 5)
    { u = arrp (seq, i, BSunit) ;
      if (bsFindKey (fat, _Previous_contig, u[0].k))
	continue ;
      x1 = u[1].i ; x2 = u[2].i ;
      if (x1 > a1 && x1 <= a2)
	{ if (!(x2 > a1 && x2 <= a2))
	    messcrash ("impossible sequence entre deux contigs") ;
	  x1 -= a1 ; x2 -= a1 ;
	  bsAddKey (obj, _Assembled_from, u[0].k) ;
	  bsAddData (obj, _bsRight, _Int, &x1) ;
	  bsAddData (obj, _bsRight, _Int, &x2) ;
	  clip = u[3].i ;
	  if (clip)
	    { bsAddData (obj, _bsRight, _Int, &clip) ;
	      clip = u[4].i ;
	      bsAddData (obj, _bsRight, _Int, &clip) ;
	    }
	}
    }
  bsSave (obj) ;
  dnaStoreDestroy (dnaKey, dna2) ;
  bsDestroy (fat) ;
}

/***************************************************************/
/* pas de realignement l'objet est considere comme bon */
/* la cle de depart est toujours conservee ; renvoie le keyset des cles
   supplementaires necessaires */
typedef struct { int a1, a2 ; } COUPLE ;
KEYSET alignToolsPurifyContig (KEY link, KEY key)
{ Array nbt, seq, dna, dna2=0, contig ;
  BSunit *u ;
  COUPLE *cou ;
  int i, *ip, *ip0, a1, a2, z1, z2, j ;
  KEY dnaKey, seqKey ;
  KEYSET ks = 0 ;
  OBJ obj ;
  BOOL debut ;
  DEFCPT look = defCptGetLook (link) ;

  if (!(obj = bsCreate(key)))
    return 0 ;
  seq = arrayCreate (200, BSunit) ;
  if (!bsGetKey (obj, _DNA, &dnaKey) || 
      !bsFindTag (obj, _Assembled_from) || !bsFlatten (obj, 5, seq))
    goto abort ;
  if (!(dna = monDnaGet (look, 0, dnaKey)))
    goto abort ;
  nbt = arrayCreate (arrayMax(dna), int) ;
  for (i = 0 ; i < arrayMax(seq) ; i += 5)
    { u = arrp (seq, i, BSunit) ;
      if (bsFindKey (obj, _Previous_contig, u[0].k))
	continue ;
      a1 = u[1].i ; a2 = u[2].i ;
      z1 = a1 < a2 ? a1 - 1 : a2 - 1 ;
      z2 = a1 < a2 ? a2 : a1 ;
      j = array (nbt, z2 - 1, int) ; /* side effect to eventually reallocate nbt */
      j = z2 - z1 ;
      ip = arrp (nbt, z1, int) - 1 ;
      if (j>0) while (ip++, j--)
	(*ip)++ ;
    }
  bsDestroy (obj) ;
  i = arrayMax (nbt) ;
  ip0 = arrp (nbt, 0, int) ;
  ip = ip0 - 1 ;
  debut = FALSE ;
  j = 0 ;
  contig = arrayCreate(10, COUPLE) ;
  if (i>0) while (ip++, i--)
    { if (*ip && !debut)
	{ debut = TRUE ;
	  arrayp (contig, j, COUPLE)->a1 = ip - ip0 ;
	}
      if (!(*ip) && debut)
	{ debut = FALSE ;
	  arrp (contig, j++, COUPLE)->a2 = ip - ip0 ;
	}
    }
  if (!(i = arrayMax (contig)))
    messcrash ("Pb in contig") ;
  if (!arrp (contig, i - 1, COUPLE)->a2) /* toujours a faire normalement vu const de nbt */
    arrp (contig, i - 1, COUPLE)->a2 = arrayMax (nbt) ;
  if (i == 1 && arrp(contig, 0, COUPLE)->a1 == 0 && 
      arrp(contig, 0, COUPLE)->a2 == arrayMax(dna))
    goto quit ;
  seqKey = key ;
  cou = arrp(contig, 0, COUPLE) - 1 ;
  if (i > 1)
    ks = keySetCreate() ;
  dna2 = dnaCopy (dna) ;
  monDnaForget (look, key) ;
  j = 0 ;
  if (i>0) while (cou++, i--)
    { alignToolsDoPurifyContig (key, seqKey, dnaKey, seq, dna2, cou->a1, cou->a2) ;
      if (i)
	{ while (++j)
	    if (!lexword2key(messprintf("%s_%d", name(key), j), &seqKey, _VSequence))
	      break ;
	  lexaddkey(messprintf("%s_%d", name(key), j), &seqKey, _VSequence) ;
	  lexaddkey(messprintf("%s_%d", name(key), j), &dnaKey, _VDNA) ;
	  keySet(ks, keySetMax(ks)) = seqKey ;
	}
    }
 quit:
  arrayDestroy (nbt) ;
  arrayDestroy (seq) ;
  arrayDestroy(contig) ;
  arrayDestroy(dna2) ;
  return ks ;
 abort:
  arrayDestroy(seq) ;
  bsDestroy (obj) ;
  return 0 ;
}

/***************************************************************/

void alignToolsAdjustLink (KEY link, KEYSET contigs, Array order)
{ int i, n = 1, max, x, y ;
  OBJ Link, obj = 0 ;
  KEY *keyp, dummy, key ;
  BSunit *unit ;

  if (!contigs && !order)
    return ;
  if (!(Link = bsUpdate (link)))
    return ;
  if (bsFindTag (Link, _Subsequence))
    bsRemove (Link) ;
  if (contigs)
    { i = arrayMax (contigs) ;
      if (!i) goto abort ;
      keyp = arrp (contigs, 0, KEY) - 1 ;
      while (keyp++, i--)
	{ if (!*keyp || !(obj = bsCreate (*keyp)))
	    continue ;
	  if (bsGetKey (obj, _DNA, &dummy) &&
	      bsGetData (obj, _bsRight, _Int, &max))
	    { bsAddKey (Link, _Subsequence, *keyp) ;
	      bsAddData (Link, _bsRight, _Int, &n) ;
	      n += max - 1 ;
	      bsAddData (Link, _bsRight, _Int, &n) ;
	      n += max > 3000 ? 1000 : 300 ;
	    }
	  bsDestroy (obj) ;
	}
    }
  else if (order)
    { for (i = 0 ; i < arrayMax (order) ; i += 2)
	{ unit = arrp (order, i, BSunit) ;
	  key = unit->k ;
	  unit++ ;
	  if (!key || !(obj = bsCreate (key)))
	    continue ;
	  if (bsGetKey (obj, _DNA, &dummy) &&
	      bsGetData (obj, _bsRight, _Int, &max))
	    { bsAddKey (Link, _Subsequence, key) ;
	      if (unit->i < 0)
		{ x = n + max - 1 ;
		  y = n ;
		}
	      else
		{ x = n ;
		  y = n + max - 1 ;
		}
	      bsAddData (Link, _bsRight, _Int, &x) ;
	      bsAddData (Link, _bsRight, _Int, &y) ;
	      n += max - 1 ;
	      n += max > 3000 ? 1000 : 300 ;
	    }
	  bsDestroy (obj) ;
	}
    }
 abort:
  bsSave (Link) ;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
