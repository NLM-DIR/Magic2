/*  File: abifix.c
 *  Author: Jean Thierry-Mieg (mieg@crick.wustl.edu)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jan 21 14:40 1999 (fw)
 * Created: Mon Dec  5 10:25:50 1994 (mieg)
 *-------------------------------------------------------------------
 */

/* @(#)abifix.c	1.11 11/12/96 */

#include "acedb.h"
#include "call.h"
#include "bs.h"
#include "dna.h"
#include "systags.h"
#include "sysclass.h"
#include "classes.h"
#include "tags.h"
#include "dnaalign.h" 
#include "query.h"
#include "pick.h"

#include "display.h"
#include "bindex.h"
#include "dict.h"

#include "acembly.h"
#include "basecall.h"
#include "freeout.h"
#include "makemrna.h"

/********************************************************************/

static int abiFixDoubleReads (KEY contig, int *ip)
{ OBJ Contig, obj ;
  Array dna = 0, consensus, consensusReverse, consensusDirect = 0, aa = 0, errArray = 0 ;
  int amax, i, j, nn, x, x1, x11, x2, x3, a1, a11, a2, a3 , oldx2 ;
  int vClipTop, vClipEnd, nmv = 0 ;
  KEY key, dnaKey, consensusKey ;
  A_ERR *ep ;
  BOOL isUp ;

  Contig = bsUpdate (contig) ;
  if (!Contig) 
    return 0 ;
  if (!bsFindTag (Contig, _Assembled_from) ||
      !bsGetKey (Contig, _DNA, &consensusKey) ||
      !(consensus = dnaGet(consensusKey)))
    { bsDestroy (Contig) ;
      return 0 ;
    }
  consensusDirect = consensus ;
  consensusReverse = dnaCopy (consensus) ;
  reverseComplement (consensusReverse) ;
  amax = arrayMax (consensus) ;

  aa = arrayCreate(3000, BSunit) ;
  bsFindTag (Contig, _Assembled_from) ;
  bsFlatten (Contig, 5, aa) ;

  for (i = 0 ; i < arrayMax(aa) ; i += 5)
    { key = arr (aa, i, BSunit).k ;
      a1 = arr (aa, i + 1, BSunit).i ;
      a2 = arr (aa, i + 2, BSunit).i ;
      x1 = arr (aa, i + 3, BSunit).i ;
      x2 = arr (aa, i + 4, BSunit).i ;

      isUp = a1 < a2 ? FALSE : TRUE ;
      if (isUp)
	{ a1 = amax - a1 - 1 ;
	  a2 = amax - a2 - 1 ;
	  consensus = consensusReverse ;
	}
      else
	consensus = consensusDirect ;
      obj = bsUpdate(key) ;
      if (!obj) continue ;
      
      if (!bsGetKey (obj, _DNA, &dnaKey) ||
	  bsFindTag (obj, _Assembled_from) ||
	  !(dna = dnaGet(dnaKey)))
	goto abort ;


      if (!bsFindTag (obj, _Old_Clipping))
	{ bsAddData (obj, _Old_Clipping, _Int, &x1) ;
	  bsAddData (obj, _bsRight, _Int, &x2) ;
	}
	
      oldx2 = x2 ;
      vClipTop = -1 ;
      if (bsGetData (obj, _Vector_Clipping, _Int, &vClipTop))
	{ vClipTop-- ;
	  bsGetData (obj, _bsRight, _Int, &vClipEnd) ;
	}

      a1-- ; x1-- ; /* no ZERO */
      errArray = 0 ;
/***************** adjust clipTop ********************/
      if (x1 < 35)
	{ x11 = x1 + 50 ;
	  a11 = a1 + 50 ;
	  if (x11 > x2) x11 = x2 ;
	  if (a11 > a2) a11 = a2 ;
	  errArray = dnaAlignCptErreur (consensus, dna, &a1, &a11, &x1, &x11) ;

	  j = arrayMax(errArray) ; 
	  if (j)
	    ep = arrp(errArray, 0, A_ERR) ;

	  if (!j || ep->iShort > x1 + 3)
	    goto topOk ;
	  if (j > 1) /* if j == 1 ep -> newClipEnd */
	    { nn = 0 ;
	      x = x1 ; j-- ;
	      while (nn++, j-- &&
		     (ep->iShort - x1 < 3 * nn) &&
		     (ep + 1)->iShort < ep->iShort + 4) ep++ ;
	    }
	  
	  switch (ep->type)
	    {
	    case TROU: 
	      x1 = ep->iShort ;
	      a1 = ep->iLong + 1 ;
	      break ;
	    case TROU_DOUBLE: 
	      x1 = ep->iShort ;
	      a1 = ep->iLong + 2 ;
	      break ;
	    case INSERTION : /*_AVANT: */
	      x1 = ep->iShort + 1 ;
	      a1 = ep->iLong ;
	      break ;
	    case INSERTION_DOUBLE: /* _AVANT: */
	      x1 = ep->iShort + 2 ;
	      a1 = ep->iLong ;
	      break ;
	    default:
	      x1 = ep->iShort + 1 ;
	      a1 = ep->iLong + 1 ;
	      break ;
	    }
	}
    topOk:
      arrayDestroy (errArray) ;
      /***************** adjust clipEnd ********************/

      x3 = arrayMax(dna) ;
      if (vClipTop != -1 &&
	  x3 > vClipEnd)
	x3 = vClipEnd ;

      a3 = a2 + x3 - x2 ;
      if (a3 > amax) 
	{
	  x3 = x3 - (a3 - amax) ;
	  a3 = amax ;
	}
      a3-- ; x3-- ;
      errArray = dnaAlignCptErreur (consensus, dna, &a1, &a3, &x1, &x3) ;
      a3++ ; x3++ ;
      
      ep = arrp(errArray, 0, A_ERR)  - 1 ;
      j = arrayMax(errArray) ; nn = 0 ;
      
      while (ep++, j--)
	{ x = ep->iShort ;
	  
	  if (x < x2 - 50) continue ;
	  nn++ ;
	  if (x > x2 && 
	      nn * 20 > x - x2 - 50)   /* % d'erreur = 1/20 */
	    { x3 = x ; a3 = ep->iLong ; break ; }
	}
      /* now remove all errors in last 10 bases */
      ep = arrp(errArray, 0, A_ERR)  - 1 ;
      j = arrayMax(errArray) ;
      while (ep++, j--)
	{ x = ep->iShort ;
	  if (x > x3 - 10) 
	    { x3 = x ; 
	      a3 = ep->iLong ;
	      break ; 
	    }
	}
      if (vClipTop != -1 &&
	  x3 > vClipEnd)
	{ a3 -= x3 - vClipEnd ;
	  x3 = vClipEnd ;
	}

	  /* because of the insertions I may be out of the consensus */
      if (a3 > amax)
	{ x3 = x3 - (a3 - amax) ;
	  a3 = amax ;
	}

      if (isUp)
	{ a1 = amax - a1 - 1 ;
	  a3 = amax - a3 - 1 ;
	}

      if (TRUE)
	{ x1++ ; 
	  if (x3 > 1000)
	    invokeDebugger() ;
	  bsAddData (obj, _Clipping, _Int, &x1) ;
	  bsAddData (obj, _bsRight, _Int, &x3) ;
	  nmv++ ; *ip += x3 - oldx2 ;
	  bsFindKey (Contig, _Assembled_from, key) ;
	  a1++ ;
	  bsAddData (Contig, _bsRight, _Int, &a1) ;
	  bsAddData (Contig, _bsRight, _Int, &a3) ;
	  bsAddData (Contig, _bsRight, _Int, &x1) ;
	  bsAddData (Contig, _bsRight, _Int, &x3) ;
	  defCptForget (0, key) ;
	}
      bsSave (obj) ;

    abort:
      arrayDestroy (dna) ;
      arrayDestroy (errArray) ;
      bsDestroy (obj) ;
    }

  bsSave (Contig) ;
  arrayDestroy (aa) ;
  arrayDestroy (consensusDirect) ;
  arrayDestroy (consensusReverse) ;
  return nmv ;
}

/********************************************************************/

#ifdef JUNK
typedef struct { int x1, x2, x3, oldX2, a1, a2, a3, nbErr ; 
		 KEY key ; Array dna, err ; } FIX ;

static void abiFixExtendTerminalReads (KEY contig, BOOL force)
{ OBJ Contig, obj ;
  Array dna = 0, consensus = 0, aa = 0, errArray = 0, fixes ;
  int 
    i, ii, i1, i2, j, nn, x, x1, x2, oldX2, 
    x3, a1 = 0, a2 = 0, a3, amax, endOther = 0 ;
  int vClipTop, vClipEnd ;
  KEY key, dnaKey, consensusKey ;
  A_ERR *ep ;
  FIX *fp, *bestfp ;
  char *cp, *cq ;

  Contig = bsUpdate (contig) ;
  if (!Contig) 
    return ;
  if (!bsFindTag (Contig, _Assembled_from) ||
      !bsGetKey (Contig, _DNA, &consensusKey) ||
      !(consensus = dnaGet(consensusKey)))
    { bsDestroy (Contig) ;
      return ;
    }

  aa = arrayCreate(3000, BSunit) ;
  bsFindTag (Contig, _Assembled_from) ;
  bsFlatten (Contig, 5, aa) ;

       /* look for terminal reads */
  fixes = arrayCreate (12, FIX) ;
  amax = arrayMax (consensus) ;
  for (ii = 0 , j = 0 ; ii < arrayMax(aa) ; ii += 5)
    { key = arr (aa, ii, BSunit).k ;
      a1 = arr (aa, ii + 1, BSunit).i ;
      a2 = arr (aa, ii + 2, BSunit).i ;
      x1 = arr (aa, ii + 3, BSunit).i ;
      x2 = arr (aa, ii + 4, BSunit).i ;
      
      obj = bsCreate(key) ;
      if (!obj || bsFindTag(obj, _Assembled_from))
	{ bsDestroy (obj) ;
	  continue ;
	}
      bsDestroy (obj) ;

      if (a1 > a2 &&
	  a1 >= endOther)
	endOther = a1 + 1 ;
	
      if (a1 > a2 || a2 < amax - 200)
	{ if (a2 > endOther)
	    endOther = a2 ;
	  continue ;
	}

      obj = bsUpdate(key) ;
      if (!obj) continue ;
      
      if (!bsGetKey (obj, _DNA, &dnaKey) ||
	  !(dna = dnaGet(dnaKey)))
	{ bsDestroy(obj) ;
	  continue ;
	}

      if (!bsFindTag (obj, _Old_Clipping))
	{ bsAddData (obj, _Old_Clipping, _Int, &x1) ;
	  bsAddData (obj, _bsRight, _Int, &x2) ;
	}
	
      vClipTop = -1 ;
      if (bsGetData (obj, _Vector_Clipping, _Int, &vClipTop))
	bsGetData (obj, _bsRight, _Int, &vClipEnd) ;

      oldX2 = 0 ;
      if (bsGetData (obj, _Old_Clipping, _Int, &x3))
	bsGetData (obj, _bsRight, _Int, &oldX2) ;

      a1-- ; x1-- ; /* no ZERO */


      if (x2 > 3000)
	invokeDebugger() ;
      bsDestroy (obj) ;
      fp = arrayp (fixes, j++, FIX) ;
      fp->key = key ;
      fp -> a1 = a1 ; fp->a2 = a2 ; fp->a3 = a2 ; /* at least */
      if (x2 > arrayMax(dna))
	x2 = arrayMax(dna) ;
      fp->x1 = x1 ; fp->x2 = x2 ; fp->oldX2 = oldX2 ;
      x3 = x2 + 100 ;
      if (x3 > arrayMax(dna))
	x3 = arrayMax(dna) ;
      if (vClipTop != -1 && x3 > vClipEnd)
	x3 = vClipEnd ;
      fp->x3 = x3 ;
      i = fp->a2 + fp->x3 - fp->x2 ;
      if (i > amax + 100)  /* do not over extend into rubbish */
	fp->x3 -= i - amax - 100 ;
      fp->dna = dna ;
    }

     /* count the errors in last 100 bases */

  arrayDestroy (aa) ;
  i = arrayMax(fixes) ;
  if (i < 1) goto abort ;

  fp = arrp (fixes, 0, FIX) - 1 ;
  while (fp++, i--)
    { key = fp->key ;
      a1 = fp->a1 ; a2 = fp->a2 - 1;
      dna = fp->dna ;
      x1 = fp->x1 ; x2 = fp->x2 - 1 ;
      x = amax - a2 ;
      x2 += x ; a2 += x ;
      if (x2 <= arrayMax(fp->dna))
	{ 
	  errArray = dnaAlignCptErreur (consensus, dna, &a1, &a2, &x1, &x2) ;
	  
	  if (errArray)
	    {
	      ep = arrp(errArray, 0, A_ERR)  - 1 ;
	      j = arrayMax(errArray) ; nn = 0 ;
	      
	      while (ep++, j--)
		{ x = ep->iLong ;
		  
		  if (x > amax - 100) 
		    nn ++ ;
		}
	      fp->nbErr = nn ;
	      arrayDestroy (errArray) ;
	    }
	}
      else
	fp->nbErr = 2 * amax  + 1000 ;
    }
  
      /* look for best sequence */
  i = arrayMax(fixes) ;
  fp = arrp (fixes, 0, FIX) - 1 ;
  bestfp = 0 ;
  j = 2 * amax ; /* plus infini */
  while (fp++, i--)
    { nn = fp->nbErr ;
      if (fp->a2 > amax - 50 &&
	  nn < j) 
	{ j = nn ; bestfp = fp ; }
    }
  
  /* extend consensus */
  if (bestfp)
    { if (bestfp->x3 > arrayMax (bestfp->dna))
	bestfp->x3 = arrayMax(bestfp->dna) ;
      i = bestfp->a2 + bestfp->x3 - bestfp->x2 ;
      if (i > arrayMax(consensus))
	array(consensus, i - 1, char) = 0 ;
      else
	arrayMax(consensus) = i ; /* This may shorten the consensus */
      bestfp->a3 = i ;
      cp = arrayp(consensus, bestfp->a2, char) ;
      cq = arrayp(bestfp->dna, bestfp->x2, char) ;
      i = bestfp->x3 - bestfp->x2 ;
      if (i>0)
	while (i--) *cp++ = *cq++ ;
      cp = arrayp(consensus, bestfp->a2, char) ;
      cq = arrayp(bestfp->dna, bestfp->x2, char) ;
      i = bestfp->x2  < bestfp->a2 ? bestfp->x2 : bestfp->a2 ;
      while (cp--, cq--, --i)
	if (*cp == 0) *cp = *cq ;
    }
  else 
    goto abort ;

  /* restudy the errors */
  amax = arrayMax (consensus) ;
  i = arrayMax(fixes) ;
  fp = arrp (fixes, 0, FIX) - 1 ;
  while (fp++, i--)
    { a1 = fp->a1 ; a2 = fp->a2 - 1 ;
      x1 = fp->x1 ; x2 = fp->x2 - 1 ;
      x3 = fp->x3 - 1 ;
      a3 = a2 + x3 - x2 ;
      if (a3 >= amax)
	{ x3 = x3 - (a3 - amax + 1) ;
	  a3 = amax - 1 ;
	}
      errArray = dnaAlignCptErreur (consensus, dna, &a1, &a3, &x1, &x3) ;      
      ep = arrp(errArray, 0, A_ERR)  - 1 ;
      j = arrayMax(errArray) ; nn = 0 ;
      
      while (ep++, j--)
	{ x = ep->iShort ;
	  
	  if (x < x2 - 100) continue ;
	  
	  if (x > x2 - 100) nn ++ ;
	  if (x > x2 &&  (nn * 20 > (x - x2 + 100)))  /*  4% d'erreur = 1/25 */
	    { x3 = x - 1 ; a3 = ep->iLong - 1 ; break ; }
	}
       /* now remove all errors in last 10 bases */
      ep = arrp(errArray, 0, A_ERR)  - 1 ;
      j = arrayMax(errArray) ; nn = 0 ;
       while (ep++, j--)
	{ x = ep->iShort ;
	  if (x < x3 - 10) continue ;
	  x3 = x - 1 ; a3 = ep->iLong - 1 ; break ; 
	}
      x3++ ; a3++ ; /* coordonne bio */
      if (a3 > amax)
	{ x3 = x3 - (a3 - amax) ;
	  a3 = amax ;
	}
     if (!force || fp != bestfp)
	{ fp->x3 = x3 ; fp->a3 = a3 ; fp->a1 = a1 ; fp->x1 = x1 ; }
      arrayDestroy (errArray) ;
    }
  
  /* check for a3 present twice */
  i = arrayMax(fixes) ;
  fp = arrp (fixes, 0, FIX) - 1 ;
  if (!force)
    {
      i1 = i2 = -1 ;      
      while (fp++, i--)
	{ if (fp->a3 >= i1)
	    { i2 = i1 ; i1 = fp->a3 ; }
	  else if (fp->a3 >= i2)
	    i2 = fp->a3 ;
	}
      if (i2 < endOther)
	i2 = endOther ;	
      if (arrayMax(fixes) == 1)
	{ fp = arrp (fixes, 0, FIX) ;
	  i2 = i1 ;
	  if (i2 > endOther)
	    { fp->x3 -= (i2 - endOther) ;
	      fp->a3 -= (i2 - endOther) ;
	      i2 = endOther ;
	    }
	  if (fp->x3 < fp->oldX2)
	    { i2 += fp->oldX2 - fp->x3 ;
	      fp->a3 += fp->oldX2 - fp->x3 ;
	      fp->x3 = fp->oldX2 ;
	    }
	}
    
      arrayMax (consensus) = i2 ;  /* second largest a3 */
      i = arrayMax(fixes) ;
      fp = arrp (fixes, 0, FIX) - 1 ;
      while (fp++, i--)
	{ if (fp->a3 > i2)
	    { fp->x3 = fp->x3 - (fp->a3 - i2) ;
	      fp->a3 = i2 ;
	    }
	}

      i = arrayMax(fixes) ;
      fp = arrp (fixes, 0, FIX) - 1 ;
      while (fp++, i--)
	{ if (fp->a3 == i2)
	    goto ok ;
	  if (fp->a1 <0)
	    messcrash("negative fp->a1") ;
	}
      if (i2 != endOther)
	messcrash ("no seg reaches a3") ;
    }
 ok:
  

  /* register */

  amax = arrayMax(consensus) ;
  dnaStoreDestroy(consensusKey, consensus) ;
  defCptForget (0, consensusKey) ;

  i = arrayMax(fixes) ;
  fp = arrp (fixes, 0, FIX) - 1 ;
  while (fp++, i--)
    { 
      if (fp->nbErr > amax)
	continue ;
      obj = bsUpdate(fp->key) ;
      x1 = fp->x1 ; x2 = fp->x2 ; x3 = fp->x3 ; 
      a1 = fp->a1 ; a2 = fp->a2 ; a3 = fp->a3 ;
/*
      printf("Terminals %s x1 = %d x2 = %d x3 = %d amax = %d a1 = %d a2 = %d a3 = %d \n",
	     name(fp->key), x1, x2, x3, amax, a1, a2, a3) ;
*/      
      if (a3 > amax)
	{ x3 = x3 - (a3 - amax) ;
	  a3 = amax ;
	}
      if (x3 > arrayMax(fp->dna))
	{ a3 -= x3 - arrayMax(fp->dna) ; /* on devrait restudy les erreurs */
	  x3 = arrayMax(fp->dna) ;
	}
      x1++ ; /* x2 ; x3 ; ulrich */
      bsAddData (obj, _Clipping, _Int, &x1) ;
      bsAddData (obj, _bsRight, _Int, &x3) ;
      bsFindKey (Contig, _Assembled_from, fp->key) ;
      a1++ ; /* a3 ; ulrich */
      bsAddData (Contig, _bsRight, _Int, &a1) ;
      bsAddData (Contig, _bsRight, _Int, &a3) ;
      bsAddData (Contig, _bsRight, _Int, &x1) ;
      bsAddData (Contig, _bsRight, _Int, &x3) ;
      defCptForget (0, fp->key) ;

      bsSave (obj) ;
    }
   /* clean up */
  if (bsFindKey(Contig, _DNA, consensusKey))
    bsAddData (Contig, _bsRight, _Int, &amax) ; 

 abort:
  bsSave (Contig) ;
  arrayDestroy (consensus) ;
  i = arrayMax(fixes) ;
  fp = arrp (fixes, 0, FIX) - 1 ;
  while (fp++, i--)
    { arrayDestroy (fp->dna) ;
      arrayDestroy (fp->err) ;
    }
  arrayDestroy (fixes) ;
}

/********************************************************************/

static void abiFixFlipAssembly (KEY contig)
{ OBJ Contig ;
  Array aa, consensus = 0 ;
  int 
    nn, i, a1, a2, clipt, clipe ;
  KEY key,  consensusKey ;

  Contig = bsUpdate (contig) ;
  if (!Contig) 
    return ;
  if (!bsFindTag (Contig, _Assembled_from) ||
      !bsGetKey (Contig, _DNA, &consensusKey) ||
      !(consensus = dnaGet(consensusKey)))
    { bsDestroy (Contig) ;
      return ;
    }

  aa = arrayCreate(3000, BSunit) ;
  bsFindTag (Contig, _Assembled_from) ;
  bsFlatten (Contig, 5, aa) ;
  bsFindTag (Contig, _Assembled_from) ;
  bsRemove (Contig) ;
  nn = arrayMax(consensus) + 1 ;
  for (i = 0 ; i < arrayMax(aa) ; i += 5)
    { key = arr (aa, i, BSunit).k ;
      a1 = arr (aa, i + 1, BSunit).i ;
      a2 = arr (aa, i + 2, BSunit).i ;
      clipt = arr (aa, i + 3, BSunit).i ;
      clipe = arr (aa, i + 4, BSunit).i ;
      a1 = nn - a1 ; a2 = nn - a2 ;
      bsAddKey (Contig, _Assembled_from, key) ;
      bsAddData (Contig, _bsRight, _Int, &a1) ;
      bsAddData (Contig, _bsRight, _Int, &a2) ;
      if (clipt)
	{ bsAddData (Contig, _bsRight, _Int, &clipt) ;
	  bsAddData (Contig, _bsRight, _Int, &clipe) ;
	}
    }
  arrayDestroy (aa) ;
  reverseComplement (consensus) ;
  dnaStoreDestroy (consensusKey, consensus) ;
  defCptForget (0, consensusKey) ;
  bsSave (Contig) ;
}
#endif
/********************************************************************/
/******************* public routines  *******************************/
/********************************************************************/

int abiFixDoubleContig (KEY link, KEY contig, int *ip)
{
  messStatus ("Double Stranding") ;
  return abiFixDoubleReads (contig, ip) ;
  /*
  abiFixExtendTerminalReads (contig, FALSE) ;

  abiFixFlipAssembly (contig) ;
  abiFixExtendTerminalReads (contig, FALSE) ;
  abiFixFlipAssembly (contig) ;
  if (link) 
    alignToolsAdjustContigSize (link, contig) ;
    */
}

/*************/

int abiFixDouble (KEY link, KEYSET ks, int *ip)
{ int n = 0, i = ks ? keySetMax(ks) : 0 ;
  KEY key, seq ;
  
  while (i--)
    { key = keySet(ks, i) ;
      if (class (key) == _VDNA)
	dnaReClass (key, &seq) ;
      else
	seq = key ;
      if (seq)
	n += abiFixDoubleContig (link, seq, ip) ;
    }
  return n ;
}

/********************************************************************/

void abiFixExtendContig (KEY link, KEY contig)
{ OBJ Contig, obj ;
  Array aa ;
  int 
    i, a1, a2, dx, amax, x1, x2, end, fair2, v1, v2, c1, c2 ;
  KEY key, dnaKey ;

  messStatus ("Extend") ;
  Contig = bsUpdate (contig) ;
  if (!Contig) 
    return ;
  if (!bsGetKey (Contig, _DNA, &dnaKey) ||
      !bsGetData (Contig, _bsRight, _Int, &amax) ||
      !bsFindTag (Contig, _Assembled_from))
    { bsDestroy (Contig) ;
      return ;
    }

  aa = arrayCreate(3000, BSunit) ;
  bsFindTag (Contig, _Assembled_from) ;
  bsFlatten (Contig, 5, aa) ;
  bsFindTag (Contig, _Assembled_from) ;
  bsRemove (Contig) ;

  for (i = 0 ; i < arrayMax(aa) ; i += 5)
    { key = arr (aa, i, BSunit).k ;
      a1 = arr (aa, i + 1, BSunit).i ;
      a2 = arr (aa, i + 2, BSunit).i ;
      c1 = arr (aa, i + 3, BSunit).i ;
      c2 = arr (aa, i + 4, BSunit).i ;

       if ((obj = bsUpdate (key)))
	 { if (bsGetKey (obj, _DNA, &dnaKey) &&
	       bsGetData (obj, _bsRight, _Int, &end) &&
	       bsGetData (obj, _Clipping, _Int, &x1) &&
	       bsGetData (obj, _bsRight, _Int, &x2))	       
	     { if (!bsFindTag (obj, _Old_Clipping))
		 { if (bsAddData (obj, _Old_Clipping, _Int, &x1))
		     bsAddData (obj, _bsRight, _Int, &x2) ;
		 }
	       if (c1 && c2) /* so dna was clipped in contig */
		 { x1 = c1 ;
		   x2 = c2 ;
		   if (bsGetData (obj, _Fair_upto, _Int, &fair2))
		     { if (fair2 < end)
			 end = fair2 ;
		     }
		   else
		     { if (end > x2 + 50)
			 end = x2 + 50 ;
		     }
		   if (bsGetData (obj, _Vector_Clipping, _Int, &v1) &&
		       bsGetData (obj, _bsRight, _Int, &v2) &&
		       v2 < end)
		     end = v2 ;
		   
		   dx = end - x2 ;
		   
		   if (a1 < a2)
		     { if (a2 + dx > amax)
			 { a2 += dx ; 
			   x2 += dx ;
			 }
		     }
		   else
		     { if (a2 - dx < 0)
			 { a2 -= dx ; 
			   x2 += dx ;
			 }
		     }
		 }
	     }
	   bsSave (obj) ;
	   defCptForget (0, key) ;
	 }
      
      bsAddKey (Contig, _Assembled_from, key) ;
      bsAddData (Contig, _bsRight, _Int, &a1) ;
      bsAddData (Contig, _bsRight, _Int, &a2) ;
      if (c1 && c2)
	{ bsAddData (Contig, _bsRight, _Int, &x1) ;
	  bsAddData (Contig, _bsRight, _Int, &x2) ;
	}
    }
  arrayDestroy (aa) ;
  bsSave (Contig) ;
  alignToolsAdjustContigSize (link, contig) ;
  dnaAlignFixContig (0, contig) ;
}

/********************************************************************/

void abiFixExtend (KEY link, KEYSET ks)
{ int i = ks ? keySetMax(ks) : 0 ;
  KEY key, seq ;

  while (i--)
    { key = keySet(ks, i) ;
      if (class(key) == _VDNA)
	dnaReClass (key, &seq) ;
      else if (class(key) == _VSequence)
	seq = key ;
      else
	continue ;
      
      abiFixExtendContig (link, seq) ;
    }
}

/********************************************************************/
 
#ifndef NON_GRAPHIC

static void dumpFinish (FILE *f, KEY contig)
{ OBJ obj = bsCreate (contig) ;
  Array aa ;
  static int nn = 1 ; /* unique seq id for finish */
  int i, a1, a2, pp, ll ;
  KEY key ;
  char buf[256], *cp, *cq ;

  if (!obj)
    return ;

  aa = arrayCreate (200, BSunit) ;
  if (bsFindTag (obj, _Assembled_from))
    { bsFlatten (obj, 3, aa) ;
      if (arrayMax (aa))
	  fprintf (f, "\nLEFT\nLEFT\n") ;	  
      for (i = 0 ; i < arrayMax(aa) ; i += 3)
	{ key = arr(aa, i, BSunit).k ;
	  a1 = arr(aa, i + 1, BSunit).i ;
	  a2 = arr(aa, i + 2, BSunit).i ;
	  if (a1 < a2)
	    { pp = a1 ; ll = a2 - a1 + 1 ;}
	  else
	    { pp = a2 ; ll = a2 - a1 - 1 ;}	

	    /* name is CX12.bp84a09.s1abi */
	  strncpy (buf, name(key), 254) ;
            /* clean up the cosmid name */
	  cp = buf ;
	  while (*cp && *cp != '.') cp++ ;
	  if (*cp == '.') *cp++ = 0 ;
	   /* clean up the ending */
	  cq = cp + strlen(cp) ;
	  while (cq > cp && !isdigit (*cq)) cq-- ;
	  if (cq > cp) *(cq + 1) = 0 ;
	  if (pp > 0)
	    fprintf (f, "%s %d %d %d 0 0 \n",
		     cp, nn++ , pp, ll) ;
	}
    }
  bsDestroy (obj) ;
  fprintf (f, "\n\n") ;
}

static FREEOPT finishTypeMenu[] =
{ {2, "Types"},
  {'x', "XL"},
   {'r', "RV"}
} ;

void abiParseFinish (FILE *f, Array aa, Array bb)
{ Array cc1, cc2, cc3 ;
  KEY kk,  type ;
  KEY new, contig ;
  int i, level = freesetfile(f,0) ;
  char *cp , buf[300] ;
  int a1, a2 ;
  OBJ obj ;

  while (freecard(level))
    {
      if (!freekey (&type, finishTypeMenu))
	continue ;
  /* get the new read */
      cp = freeword () ;
      lexaddkey (cp, &new, _VSequence) ;
  /* find the read we started from */
      if ((cp = freeword ()))
	strncpy (buf, cp,298) ;
      if (!freeint (&a1) || !freeint(&a2))
	continue ;
      cc1 = query (aa, messprintf("IS *%s*", cp)) ;
      cc2 = cc3 = 0 ;
      contig = 0 ;
  /* find the correct contig */
      for (i= 0 ; i < keySetMax(cc1) ; i++)
	{ kk = keySet(cc1, i) ;
	  keySetDestroy(cc2) ;
	  keySetDestroy(cc3) ;
	  cc2 = queryKey (kk, "> Assembled_into") ;
	  cc3 = keySetAND (cc2, bb) ;
	  if (keySetMax(cc3))
	    { contig = keySet(cc3, 0) ;
	      break ;
	    }
	}
      keySetDestroy(cc1) ;
      keySetDestroy(cc2) ;
      keySetDestroy(cc3) ;
      if (!contig)
	continue ;

      if ((obj = bsUpdate(contig)))
	{ bsAddKey (obj, _Assembled_from, new) ;
	  bsAddData (obj, _bsRight, _Int, &a1) ;
	  bsAddData (obj, _bsRight, _Int, &a2) ;
	  bsSave (obj) ;
	}
      if ((obj = bsUpdate(new)))
	{ bsAddTag (obj, _Proposed) ;
	  bsSave (obj) ;
	}
    }
}

void abiFixDoFinish (KEY link)
{
  KEYSET aa = 0, bb = 0, old ;
  int i, j ;
  FILE *f ;
  char *tmpName  = 0 , buf[100] , *cp ;
  static BOOL firstPass = TRUE ;

    /* get whole link and its hierarchy */
  i = 0 ;
  aa = bb = keySetCreate () ;
  aa = keySetCreate () ;
  keySet(aa, 0) = keySet (bb, i++) = link ;
  while (TRUE)
    { old = aa ;
      aa = query (old, "> Subsequence") ;
      keySetDestroy (old) ;
      if (!keySetMax(aa))
	break ;
      j = keySetMax (aa) ;
      while (j--)
	keySet(bb, i++) = keySet (aa, j) ;
    }
  
  keySetDestroy (aa) ;
  
  keySetSort (bb) ;
  keySetCompress (bb) ;
  

  if (!keySetMax(bb))
    { messout("No reads to finish") ;
      goto abort ;
    }

  aa = query (bb ,"> Assembled_from") ;
  strncpy (buf, name(keySet(aa,0)), 99) ;

  cp = buf ;
  while (*cp && *cp != '.') cp++ ;
  *cp = 0 ;

  f = filtmpopen(&tmpName,"w") ;
  if (!f)
    goto abort ;
  i = keySetMax(bb) ;
  while (i--)
    dumpFinish (f, keySet(bb, i)) ;
  filclose (f) ;

  if (firstPass)
    { firstPass = FALSE ;
      messout ("wscripts/finish \nContributed by Gabor Marth (wash U)") ;
    }
  messStatus ("Finishing") ;
  if (!callScript("ace2finish", messprintf ( "%s %s", buf, tmpName)))
    { if ((f = filopen(tmpName,"a","r")))
	{ externalFileDisplay ("Finish", f, 0) ;
	  filremove(tmpName,"a") ;
	}
      if ((f = filopen(tmpName,"c","r")))
	{ abiParseFinish (f, aa, bb) ;
	  /*      filremove(tmpName,"c") ;
	   */
	  callScript ("emacs", messprintf("%s.c &", tmpName)) ;
	}	
    }

  filtmpremove(tmpName) ;
  display (link, 0, 0) ;
 abort:
  keySetDestroy(aa) ;
  keySetDestroy(bb) ;
}

#endif /* !NON_GRAPHIC */

/*********************************************************/

int laneGlobalOrder  (const void *a, const void *b)
{ const LANE* la = (const LANE*)a,  *lb = (const LANE*) b ;

  BOOL
    ra = la->upSequence ,
    rb = lb->upSequence ;
  int 
    xa = ra ? la->x2 : la->x1 ,
    xb = rb ? lb->x2 : lb->x1 ;
/*
  if (ra && !rb)
    return -1 ;
  if (rb && !ra)
    return 1 ;
*/
  return
    xa < xb ? -1 : (xa == xb ? 0 : 1) ;
}

/********************************************************************/

static void traceSort (LOOK look)
{ arraySort (look->lanes, laneGlobalOrder) ;
}

static void reRegister (Array lanes, OBJ obj, int sens)
{ int ss ;
  int x1, x2, c1, c2,  iLane ;
  LANE *lane ;

  if (!obj) return ;
 
  for (iLane = 0 ; iLane < arrayMax(lanes) ; iLane++)
    { 
      lane = arrp (lanes, iLane, LANE) ;
      ss = lane->upSequence ? -sens : sens ;
      if (ss > 0)
	{
	  x1 = lane->x1 + 1 ;
	  x2 = lane->x2 + 2 ;
	}
      else
	{
	  x1 = lane->x1 + 1 ;
	  x2 = lane->x2 ;
	}
      c1 = lane->clipTop+ 1 ;
      c2 = lane->clipEnd ;
      
      if (lane->isAligned)
	{
	  if (!bsFindKey (obj, _Aligned, lane->key))
	    bsAddKey (obj, _Aligned, lane->key) ;
	}
      else
	{ 
	  if (!bsFindKey (obj, _Assembled_from, lane->key))
	    bsAddKey (obj, _Assembled_from, lane->key) ;
	}

      bsAddData (obj, _bsRight, _Int, &x1) ;
      bsAddData (obj, _bsRight, _Int, &x2) ;
      bsAddData (obj, _bsRight, _Int, &c1) ;
      bsAddData (obj, _bsRight, _Int, &c2) ;
    }
}


typedef struct {LANE* lane ; BOOL up, in ; int pos ;
		char *cp ; JUMP* jp ; } READ ;

void trackContig (LOOK look, int z1, int z2, BOOL whole)
{ int 
    i, dx0, dq, dx, dxw, pos, dxModif, maxModif, nrtot,
    n, nr, nrr, nRes, iLane, iLanesMax, iRead, nRead ;
  char *cp, *cq, *cq0, *cp1, *cq1, *cqMax, cc1, cc2, ccq ;
  Array 
    reads, result,
    consensus = look->dna , 
    lanes = look->lanes ; 
  int done[16], nq, sens ;
  float nn[16] ;
  JUMP *jp ;
  extern JUMP jumper[], reAligner[] ;
  LANE *lane ; READ *rr ;
  OBJ obj ;


  if (whole || z2 >= arrayMax (consensus))
    z2 = arrayMax (consensus) ;

  traceSort (look) ;
  i = 16 ; while (i--) done[i] = 0 ;
  i = 16 ; while (i--) nn[i] = 0 ;

  iLane = 0 ; iLanesMax = lanes ? arrayMax (lanes) : 0 ;
  if (whole && iLanesMax)
    { 
      lane = arrp (lanes, 0, LANE) ;
      z1 = lane->upSequence ? lane->x2 : lane->x1 ;
    }

  if (z2 < z1)
    return ;
 
  obj = bsUpdate (look->key) ;
  result = arrayCreate (2*arrayMax(look->dna) + 2000, char) ;
  reads = arrayCreate (50, READ) ; nRead = 0 ;
  pos = z1 ;
  if (z1 > 0 && !whole) /* copier le debut, mais pas si whole */
    for (nRes = 0 ; nRes < z1 ; nRes++)
      array (result, nRes, char) = array (consensus, nRes, char) ;
  if (z1 < 0)
    { /* shift inside obj */
      float dq = -z1 , zero = 0;
      if (obj)
	bsCoordShift (obj, _Structure, zero, dq, TRUE) ;
      /*  shilft lanes */
      iLane = iLanesMax ;
      while (iLane--)
	{ lane = arrp (lanes, iLane, LANE) ;
	  lane->x1 += dq ; lane->x2 += dq ;
	}
      /* shift the consensus */
      arrayDestroy (look->dnaR) ;
      i = arrayMax(consensus) ;
      array(consensus, i - 1 - z1, char) = N_ ; /* make room */
      while (i-- > 0)
	arr (consensus, i- z1, char) = arr (consensus, i, char) ;  
      i = -z1 ;
      while (i-- > 0)
	arr (consensus, i- z1, char) = N_ ;
      look->dnaR = dnaCopy (look->dna) ;
      reverseComplement (look->dnaR) ;
      z2 += -z1 ; z1 = 0 ;
    }
  
  iLane = 0 ;
  cq0 = arrp (consensus, 0, char) ;
  cqMax = cq0 + z2 ;
  nRes = z1 > 0 ? z1 : 0 ; 
  arrayMax (result) = 0 ; array (result, nRes, char) = 0 ; 
  pos = z1 ; /* may be negative */
  cq = cq0 + pos ;

  dq = 0 ;
  while (cq < cqMax || (nRead && whole)) 
    { 
      dx0 = 0 ; pos = cq - cq0 ;
      while (iLane < iLanesMax &&
	     (lane = arrp (lanes, iLane, LANE)))
	{ 
	  if (
	      (!lane->scf && !traceGetLane (look, lane)) || /* non usable */
	      lane->isAligned ||  /* uninteresting */
	      ((lane->upSequence && lane->x1 < pos) ||
	       (!lane->upSequence && lane->x2 < pos)) ||   /* too far down in consensus */
	      (!lane->upSequence &&  /* too far in the read */
	       pos - lane->x1  >= lane->clipEnd - lane->clipTop - 1) ||
	      (lane->upSequence &&    /* too far in the read */
	       pos - lane->x2  >= lane->clipEnd - lane->clipTop)
	      )
	   { iLane++ ; continue ; }  /* forget that one */

	   if ((lane->upSequence && lane->x2 > pos) ||
	       (!lane->upSequence && lane->x1 > pos)  )
	   break ; /* too early */
	   iLane++ ; /* avoid recursion */
	  /* insert */
	   nRead++ ;
	  for (i = 0 ; i < arrayMax(reads) ; i++)
	    { rr = arrp (reads, i, READ) ;
	      if (!rr->in) /* search empty spot */
		break ;
	    }

	  rr = arrayp (reads, i, READ) ;
	  rr->in = TRUE ;
	  rr->lane = lane ;
	  if (lane->upSequence)
	    { rr->up = TRUE ;
	      rr->pos = lane->clipEnd - pos + lane->x2 ;
	      rr->cp = arrp (lane->dna, rr->pos, char) ;
	    }
	  else
	    { rr->up = FALSE ;
	      rr->pos = lane->clipTop + pos - lane->x1 ;
	      rr->cp = arrp (lane->dna, rr->pos, char) ;
	    }
	  /* realign */
	  jp = reAligner ; jp-- ;
	  sens = rr->up ? -1 : 1 ;
          cp = rr->cp ;
	  while (jp++) /* 0, 0, 0, 0 = last always accepted */
	    { if (cq < cq0 || cq >= cqMax)
	       { 
		 if (jp->ok)
		   continue ;
		 else
		   break ;
	       }
	      n = 0 ; i = jp->ok ;
	      cp1 = cp + sens * ( jp->dcp - 1) ; cq1 = cq + jp->dcq ;
	      while (cp1 += sens, i--)
		{ cc1 = (*cp1 & 0x0F) ;
		  if (sens == -1) cc1 = complementBase((cc1)) ;
		  if (cc1 != (*cq1++ & 0x0F))
		    goto nextjp1 ;
		}
	      n = 0 ; i = jp->lng ;
	      cp1 = cp + sens * (jp->dcp - 1) ; cq1 = cq + jp->dcq ;
	      while (cp1 += sens, i--)
		{ cc1 = (*cp1 & 0x0F) ;
		  if (sens == -1) cc1 = complementBase((cc1)) ;
		  if (cc1 != (*cq1++ & 0x0F))
		    if (++n >= jp->ok)
		      goto nextjp1 ;
		}
	      goto foundjp1 ;
	    nextjp1:
	      continue ;	  
	    }
	foundjp1:
	  if (jp)
	    { dx = jp->dcp - jp->dcq ;
	      if (rr->up) dx = - dx ; 
	      rr->pos += dx ;
	      rr->cp += dx ;
	    }
	  if (!rr->up)
	    rr->lane->x1 = nRes - (rr->pos - lane->clipTop) ;
	  else
	    rr->lane->x2 = nRes + (rr->pos - lane->clipEnd) ;   
	}
      

      maxModif = dx0 =  0 ;
      i = 16 ; while (i--) done[i] = 0 ; cc1 = 0 ;
    lao:
      dx = 0 ; nr = nrr = nq = nrtot = 0 ; dxModif = 0 ;
      i = 16 ; while (i--) nn[i] = 0 ;
      iRead = arrayMax(reads) ; ccq = cq >= cq0  && cq < cqMax ? (*cq  & 0x0f) : cc1 ;
      if (ccq == N_) ccq = 0 ;
      while (iRead--)
	{ rr = arrp (reads, iRead, READ) ;
	  if (!rr->in)
	    continue ;
	  dxw = 10000 ; /* poids de cette base en nombre de base */
	  /* dxw := in the future: 55 70 90 100 110 130 145 */
          dxModif += dxw - 10000 ; nrr += 10000 ;
	  /* if (rr->excellent) */
	    nr = 10000  + (/* rr->excellent */  - rr->pos) ;  /* N_ counts like 1/3 */
	  rr->jp = 0 ;
	  cp = rr->cp ;
	  cc1 = *cp & 0x0f ; if (cc1 == N_) { cc1 = 0 ; nr -= 7000 ; }  /* N_ counts like 1/3 */
	  if (rr->up) cc1 = complementBase(cc1) ;
	  nn[(int)cc1] += nr ; /* this favors good quality i hope */
	  nrtot += nr ;
	  if (cc1 & ccq)  /* see the N_ */
	    { nq++ ; continue ; }
	  jp = jumper ; jp-- ;
	  sens = rr->up ? -1 : 1 ;
	  while (jp++) /* 1, 1, 0, 0 = last always accepted */
	    { 
	      if (cq < cq0 || cq >= cqMax)
	       { 
		 if (jp->ok)
		   continue ;
		 else
		   break ;
	       }
	      n = 0 ; i = jp->ok ;
	      cp1 = cp + sens * ( jp->dcp - 1) ; cq1 = cq + jp->dcq ;
	      while (cp1 += sens, i--)
		{ cc2 = (*cp1 & 0x0F) ;
		  if (sens == -1) cc2 = complementBase(cc2) ;
		    if (cc2 != (*cq1++ & 0x0F))
		      goto nextjp ;
		}
	      n = 0 ; i = jp->lng ;
	      cp1 = cp + sens * (jp->dcp - 1) ; cq1 = cq + jp->dcq ;
	      while (cp1 += sens, i--)
		{ cc2 = (*cp1 & 0x0F) ;
		  if (sens == -1) cc2 = complementBase(cc2) ;
		  if (cc2 != (*cq1++ & 0x0F))
		    if (++n > jp->ok)
		      goto nextjp ;
		}
	      goto foundjp ;
	    nextjp:
	      continue ;	  
	    }
	foundjp:
	  rr->jp = jp ;
	  if (rr->jp)
	    { 
	      if (jp->dcp > jp->dcq) dx += 10000 ; 
	      else if (jp->dcp < jp->dcq) dx -= 10000 ;
	      dxModif += 10000 * (jp->dcp - jp->dcq) ; 
	    }
	}
      if (!dx0) maxModif = dxModif ;
      if (2*dx > nrr && dx0 >= 0 && nrr * dx0 < maxModif 
	  && cq > cq0 && cq < cqMax)
	{ dx0++ ; dq += 1 ; *--cq  = cc1 ; goto lao ; }
      if (2*dx < - nrr && dx0 <= 0 && nrr * dx0 > maxModif && cq > cq0)
	{ dx0-- ; dq -= 1 ; cq++ ; goto lao ; }

      /* en moyenne on a maintenant le bon nombre de bases */
      if (2 * nn[(int)(ccq)] < nrtot)
	{ int f = 0 ;
	  cc1 = N_ ;
	  if (nn[A_] > f) { f = nn[A_] ; cc1 = A_ ; }
	  if (nn[T_] > f) { f = nn[T_] ; cc1 = T_ ; }
	  if (nn[G_] > f) { f = nn[G_] ; cc1 = G_ ; }
	  if (nn[C_] > f) { f = nn[C_] ; cc1 = C_ ; }
	  
	  if (ccq != cc1 && !done[(int)(cc1)])
	    { done[(int)(cc1)] = TRUE ; if (cq >= cq0 && cq < cqMax) *cq = cc1 ; goto lao ; }
	}
      /* eliminate finished reads */
      i = arrayMax (reads) ;
      while (i--)
	{ rr = arrp (reads, i, READ) ;
	  if (!rr->in) continue ;
	  if (!rr->up && rr->pos >= rr->lane->clipEnd - 1)
	    { nRead-- ; 
	      rr->in = FALSE ;
	      rr->lane->x2 = nRes + 1 ;
	      arrayDestroy (rr->lane->errArray) ; /* now obsolete */
	    }
	  else if (rr->up && rr->pos <= rr->lane->clipTop)
	    { nRead-- ;
	      rr->in = FALSE ;
	      rr->lane->x1 = nRes ;
	      arrayDestroy (rr->lane->errArray) ; /* now obsolete */ 
	    }
	}
      /* register new solution */
      array(result, nRes++, char) = ccq ;
      arrayMax (result) = nRes ;

      /* move forward in active reads */
      i = arrayMax(reads) ; /* cc = ccq ; */
      while (i--)
	{ rr = arrp (reads, i, READ) ;
	  if (!rr->in)
	    continue ;
	  if (rr->jp && rr->jp->dcp > rr->jp->dcq)
	    dx = 1 + rr->jp->dcp ; /* jump */
	  else if (rr->jp && rr->jp->dcp < rr->jp->dcq && nRead > 1)
	    dx = 0 ; /* wait */
	  else
	    dx = 1 ;
	  rr->cp += rr->up ? -dx : dx ;
	  rr->pos += rr->up ? -dx : dx ;
	}
      cq++ ;   pos = cq - cq0 ;
    }

  iLane = 0 ; z2 = pos ;
  if (obj)
    bsCoordShift (obj, _Structure, (float) pos, (float) dq, TRUE) ;
      /* eliminate finished reads */
  i = arrayMax (reads) ;
  while (i--)
    { rr = arrp (reads, i, READ) ;
      if (!rr->in) continue ;
      if (!rr->up)
      { nRead-- ; 
        rr->in = FALSE ;
        rr->lane->x2 += dq ;
        arrayDestroy (rr->lane->errArray) ; /* now obsolete */
      }
    else if (rr->up && rr->pos <= rr->lane->clipTop)
      { nRead-- ;
        rr->in = FALSE ;
        rr->lane->x1 += dq ;
        arrayDestroy (rr->lane->errArray) ; /* now obsolete */ 
      }
    }
  /* shift all other reads */
  while (iLane < iLanesMax &&
	 (lane = arrp (lanes, iLane++, LANE)))
    { if (lane->x1 > pos) lane->x1 +=  dq ;
      if (lane->x2 > pos) lane->x2 +=  dq ;
      arrayDestroy (lane->errArray) ; /* now obsolete, could be shifted en pos, dq */	    
      if (lane->x1 > z2) z2 = lane->x1 ;
      if (lane->x2 > z2) z2 = lane->x2 ;
    }
      
  if (obj)
    { reRegister (lanes, obj, look->sens) ;
      bsSave (obj) ;
    }
  if (z2 < 0) z2 = 0 ; /* possibly trim off the end */
  cq-- ; cqMax = cq0 + arrayMax (consensus) ;
  while (++cq < cq0 + z2) /* but copy over the explored region to the end of last read */
    array (result, nRes++, char) = cq >= cq0 && cq < cqMax ? *cq : N_ ;
  arrayDestroy (reads) ;
  arrayDestroy (look->dna) ;
  arrayDestroy (look->dnaR) ;
  look->dna = dnaCopy(result) ;
  look->dnaR = dnaCopy (look->dna) ;
  reverseComplement (look->dnaR) ;
  if (look->sens < 0)
    reverseComplement (result) ;
  if (look->dnaKey)
    dnaStoreDestroy (look->dnaKey, result) ;
  else 
    arrayDestroy (result) ;
}

void trackFixContig (KEY link, KEY contig, KEY contigDnaKey, Array af, Array dna)
{ int i, j ;
  LOOK look =(LOOK)  messalloc (sizeof (struct LookStruct)) ;
  BSunit *uu ;
  LANE *lane = 0 ;
  char *cp, *cq ;
       
  look->magic = (void*)1 ;
			
  look->key = contig ;
  look->dnaKey = contigDnaKey ;
  look->dna = dnaCopy (dna) ;
  look->dnaR = dnaCopy (dna) ;
  reverseComplement (look->dnaR) ;
  look->lanes = arrayCreate (arrayMax(af)/4, LANE) ;

  for (i = 0, j = 0 ; i < arrayMax(af) ; i += 4)
    { lane = arrayp (look->lanes, j++, LANE) ;
      uu = arrp (af, i, BSunit) ;
      lane->key = uu[0].k ; 
      if (uu[1].i < uu[2].i)
	{ lane->upSequence = FALSE ; 
	  lane->x1 = uu[1].i - 1 ;
	  lane->x2 = uu[2].i ;
	  lane->x3 = uu[2].i ;
	}
      else
	{ lane->upSequence = TRUE ; 
	  lane->x1 = uu[1].i - 1 ;
	  lane->x2 = uu[2].i - 2 ;
	  lane->x3 = uu[2].i - 2 ;
	}
      lane->dna = blyDnaGet (link, contig, lane->key) ; /* clipped dna */
      lane->clipTop = 0 ;
      lane->minBase = uu[3].i ; /* surcharge */
      lane->clipEnd = arrayMax (lane->dna) ;
      lane->scf = 1 ;
    }

  trackContig (look, 0, arrayMax(look->dna), TRUE) ;

  /* reregister for ulrich */

  cp = arrp (look->dna, 0, char) ;
  i = arrayMax(look->dna) ;
  array (dna, i,char) = 0 ;
  arrayMax(dna) = i ;
  cq = arrp (dna, 0, char) ;
  while (i--) *cq++ = *cp++ ;

  for (i = 0, j = 0 ; i < arrayMax(af) ; i += 4)
    { lane = arrayp (look->lanes, j++, LANE) ;
      uu = arrayp (af, i, BSunit) ;
      uu[0].k = lane->key ;
      if (lane->upSequence)
	{ uu[1].i = lane->x1 + 1 ;
	  uu[2].i = lane->x2 + 2;
	}
      else
	{ uu[1].i = lane->x1 + 1 ;
	  uu[2].i = lane->x2 ;
	}
      uu[3].i = lane->minBase ;   /* overloading */
    }

  arrayDestroy (look->dna) ;
  arrayDestroy (look->dnaR) ;
  arrayDestroy (look->lanes) ;
  look->magic = 0 ;
  messfree (look) ;
}

/********************* end of file **********************************/
/********************************************************************/
 
static BOOL doGetVector (KEY *vecp, KEY seq,  Array v1, Array v2, int *posp, int *posMaxp, BOOL *isCodagep)
{ KEYSET ks = 0 ;  KEY vec, newVec ;
  static KEY  _sv = 0, _lm, _lr, _Maximal_position ;
  int i ;
  OBJ Seq = 0, MainClone = 0, Vec = 0 ;
  char *cp ;
  static KEY mainClone = 0 ;
  char *codage = 0 ;

  if (!_sv)
    {
      lexaddkey ("Sequencing_Vector", &_sv, _VSystem) ;
      lexaddkey ("Left_Motif", &_lm, _VSystem) ;
      lexaddkey ("Right_Motif", &_lr, _VSystem) ;
      lexaddkey ("Maximal_position", &_Maximal_position, _VSystem) ;
    }
  if (!mainClone)
    { if ((ks = query (0, "Find Clone Main_Clone")) &&
	  keySetMax (ks) > 0)
	mainClone = keySet (ks, 0) ;
      else
	mainClone = 1 ;
      keySetDestroy (ks) ;
    }
      
  /* search the vector in the read then in main clone */
  vec = *vecp ;
lao:
  if (!(Seq = bsCreate (seq)))
    return FALSE ;
  codage = 0 ; newVec = 0 ;
  if (!vec)
    { 
      bsGetKey (Seq, _sv, &newVec) ;
      if (!newVec && mainClone && mainClone > 1)
	{
	  if ((MainClone = bsCreate (mainClone))) 
	    {
	      if (bsGetKey (MainClone, _sv, &newVec))
		bsGetData (MainClone, _bsRight, _Text, &codage) ;
	      bsDestroy (MainClone) ;
	    }
	}
    }
  else
    { 
      if (bsFindKey (Seq, _sv, vec))
	bsGetKey (Seq, _bsDown, &newVec) ;
      else if (mainClone && mainClone > 1)
	{
	  if ((MainClone = bsCreate (mainClone)))
	    {
	      if (bsFindKey (MainClone, _sv, vec))
		{
		  if (bsGetKey (MainClone, _bsDown, &newVec))
		    bsGetData (MainClone, _bsRight, _Text, &codage) ;
		}
	      bsDestroy (MainClone) ;
	    }
	}
    }
  bsDestroy (Seq) ;

  vec = newVec ;
  if (!vec)
    return FALSE ;
 
  Vec = bsCreate (vec) ;
  if (!Vec)
    goto lao ;

  if (!codage)
    bsGetData (Vec, str2tag("Codage"), _Text, &codage) ;

  if (codage && *codage && !pickMatch(name(seq), codage))
    {
      KEYSET vKs = queryKey (seq, codage) ;
      int nnvKs = keySetMax (vKs) ;

      keySetDestroy (vKs) ;
      if (!nnvKs)
	{
	  bsDestroy (Vec) ;
	  goto lao ;
	}
    }

  *isCodagep = codage ? TRUE : FALSE ;
  *posp = *posMaxp = 0 ;
  bsGetData (Vec, _Position, _Int, posp) ;
  bsGetData (Vec, _Maximal_position, _Int, posMaxp) ;
  v1 = arrayReCreate (v1, 20, char) ;
  v2 = arrayReCreate (v2, 20, char) ;
  if (bsGetData (Vec, _lm, _Text, &cp))
    { i = 0 ;
      while ((array (v1, i++, char) = dnaEncodeChar[(int)(*cp++)])) ;
      arrayMax(v1) -- ;
    }
  if (bsGetData (Vec, _lr, _Text, &cp))
    { i = 0 ;
      while ((array (v2, i++, char) = dnaEncodeChar[(int)(*cp++)])) ;
      arrayMax(v2) -- ;
    }
  bsDestroy (Vec) ;
  *vecp = vec ;
  return TRUE ;
}

/********************************************************************/
 
static BOOL doGetSl (KEY *vecp, Array v1)
{ 
  static KEYSET ksA = 0 ;
  KEY vec ;
  int i ;
  OBJ Vec = 0 ;
  char *cp ;

  vec = *vecp ;
  *vecp = 0 ;

  if (!ksA)
    {
      KEYSET ks = query (0, "Find Motif IS SL* && Match_sequence") ;
  
      if (keySetMax (ks))
	ksA = keySetAlphaHeap (ks, keySetMax (ks)) ;
      keySetDestroy (ks) ;
    }

  if (!ksA)
    return FALSE ;
  i = -1 ;
  if (vec) /* find next one */
    {
      for (i = 0 ; i < keySetMax (ksA) ; i++)
	if (vec == keySet (ksA, i))
	  break ;
      
      if (i >= keySetMax(ksA) - 1)
	return FALSE ;
    }
   
  vec = keySet (ksA, i+1) ;
 
  if ((Vec = bsCreate (vec)))
    {
      i = 0 ;
      v1 = arrayReCreate (v1, 20, char) ;

      if (bsGetData (Vec, str2tag ("Match_sequence"), _Text, &cp))
	{
	  while ((array (v1, i++, char) = dnaEncodeChar[(int)(*cp++)])) ;
	  arrayMax(v1) -- ;
	}
      if (i > 1)
	*vecp = vec ;
      bsDestroy (Vec) ;
    }

  return *vecp ? TRUE : FALSE ;  
}

/***********************/
/* force
   0x1 : external flag to overide the tag
   0x2: vector not found by exact code
   0x4: sl not found by exact code
*/
static int doTrackVector (KEY seq, Array v1, Array v2, int force) 
{ 
  Array dna = 0, dnaR = 0, vector, units ;
  BSunit *uu ;
  BOOL isCodage = FALSE, isLinker = FALSE ;
  OBJ obj, Tg, Vec ;
  char *cp, *cq, linker [256] ;
  KEY dnaKey, vec, tg ;
  int 
    ii, minLen = 8, nok, nerr, best, myvec, myok1, myerr1, oldVecPos,
    old, c1, c2, ok = 0, pos,
    pp = 0 , pp1, pp2, min, lg, lg1, nn, 
    tg1, /* beginning of current alignment to the genome */
    /*  tgerr=0, nb of errors in first exon */
    pos1, pos2, /* begin end of insert */ 
    posMax, posMax2,  /* never pi > posMax, careful after posMax2 = (pp1+posMax) / 2 */
    lg2 ;
  Array err = 0 ;
  int seuil = 3 ;
  
  if (class(seq) == _VDNA && !dnaReClass(seq, &seq))
    return 0 ;
  
  obj = bsCreate (seq) ;
  if (!obj)
    return 0 ;
  
  if (!bsGetKey (obj, _DNA, &dnaKey) ||
      (bsFindTag (obj, _Vector_Clipping) && !(force & 0x1)))
    goto abort ;
  
  linker[0] = 0 ;
  tg1 = -1 ;
  if (bsGetKey (obj, _From_gene, &tg) &&
      bsFindTag (obj, _Forward) &&
      (Tg = bsCreate (tg)))
    {
      units = arrayCreate (256, BSunit) ;
      if (bsGetArray (Tg, _Assembled_from, units, 6))
	for (ii = 0 ; ii < arrayMax (units) ; ii += 6)
	  {
	    uu = arrp (units, ii, BSunit) ;
	    if (uu[2].k == seq)
	      { tg1 = uu[3].i ; /* tgerr = uu[5].i ; */ break ; }
	  }
      bsDestroy (Tg) ;
      arrayDestroy (units) ;
    }

  dna = 0 ;
  dnaR = 0 ;

  if (!dna) 
    dna = dnaGet (dnaKey) ;
  if (!dna)
    goto abort ;
  pp1 = 1 ; pp2 = arrayMax (dna) ;

  oldVecPos = 0 ;
  bsGetData (obj, _Vector_clipping, _Int, &oldVecPos) ;
  vec = 0 ;
  myerr1 = myok1 = myvec = ok = 0 ; best = 0 ;
  while ((force & 0x4) && tg1 > 0 && !ok && doGetSl (&vec, v1))
    { 

      /* search pre_insert of sl at beginning of alignment */
      vector = v1 ;
      lg = old = arrayMax (vector) ;
      if (lg < 20)
	continue ; 
      for (pos = tg1 + 3 ; !ok && pos > tg1 - 6 && pos > 6 ; pos--)
	{
	  nn = lg ; nerr = nok = 0 ;
	  if (nn > pos) nn = pos ;
	  cp = arrp (vector, lg - 1, char) ;
	  cq = arrp (dna, pos, char) ;
	  while (nn-- >= 0 && nerr < 3)
	    {
	      if (pos - nok < oldVecPos)
		break ;
	      if (*cp-- == *cq--) nok++ ;
	      else { nok++ ; nerr++ ; }

	      if (nok >= 6 && nerr < 3 && nok - 6 * nerr > best) /* the vector matches the sequence */
		{ 
		  ok |= 1 ;
		  min = nerr ;
		  pp1 = pos + 1 ; /* first base of the sl */
		  myvec = vec ; 
		  myok1 = nok ;
		  myerr1 = nerr ;
		  best =  nok - 6 * nerr ;
		}
	    }
	}
    }
  if (ok)
    {
      freeOutf ("Sequence %s\n"
		, freeprotect (name(seq))
		) ;
      freeOutf ("Transpliced_to2 %s %d \"Match_at_gene_%d/%d\"\n\n"
		, freeprotect (name(myvec)), pp1  
		, myok1 - myerr1, myok1
		) ;
    }

  pp1 = 1 ; pp2 = arrayMax (dna) ;
  vec = 0 ; ok = 0 ;
  while (!ok && (force & 0x2) && doGetVector (&vec, seq, v1, v2, &pos1, &posMax, &isCodage))
    { 
      /* search pre_insert of vector among 50 bp at beginning of sequence */
      if (!dna) 
	dna = dnaGet (dnaKey) ;
      if (!dna)
	goto abort ;

      min = 100 ;
      vector = v1 ;
      lg = old = arrayMax (vector) ;
      if (isCodage && pos1 && pos1 < 50)
	for (pos = 0 ; min && pos + 8 < arrayMax (vector) ; pos++)
	  {
	    nn = pos ;
	    lg1 = arrayMax (vector) ; lg1 -= pos  ;
	    cp = arrp (vector, pos, char) ;
	    cq = arrp (dna, 0, char) ;
	    while (nn < arrayMax (vector) && (*cp++ == *cq++)) nn++ ;
	    if (nn ==  arrayMax (vector) && /* the vector exactly corresponds to the beginning of the sequence */
		(lg1 > tg1- 5) &&  /* otherwise , it is preferable to search at tg1, as below */
		(lg1 > 12 || (tg1 && tg1 > lg1 - 6))) /* be sure it is really the vector, not a match to the genome */
	      { 
		min = 0 ;
		ok |= 1 ;
		pp = lg1 + 1 ;
	      }
	  }
      /* search for the vector sequence just around the beginning of the alignment */
      reverseComplement (vector) ;
      if (!dnaR)
	{
	  dnaR = dnaCopy (dna) ;
	  reverseComplement (dnaR) ;
	}
      if (!ok && lg > 15 && tg1 >= 0)
	{
	  best = 0 ;
	  pp1 = (long)arrayMax (dnaR) - tg1 + 1 - 15 ;
	  if (pp1 < 0) pp1 = 0 ;
	  for (ii = 0, pos = pp1 ; min > 0 && ii < 50 && pos < arrayMax(dnaR) - 20; pos++, ii++) 
	    {
	      lg2 = lg ; pos2 = pos + lg ;
	      if (pos2 + 1 >= arrayMax(dnaR)) pos2 = arrayMax (dnaR) -1 ;
	      err = trackErrors (vector, 0, &lg2, dnaR, pos, &pos2, &nn, err, 0) ;
	      nerr = arrayMax (err) ;
	      if (
		  ((lg2 == lg && min > nerr) || (lg2 > minLen && lg2 - 5 * nerr > best) ) &&
		  ( nerr == 0 
		    || arr(err, 0, A_ERR).iShort > 1
		    || arr(err, 0, A_ERR).type ==	AMBIGUE						   
		    || arr(err, 0, A_ERR).type ==	ERREUR
		    )  /* avoid a jump on base 1, by scanning on pos we can find a non-jump */
		  )
	      { min = nerr ; pp = arrayMax(dna) - pos + 1 ; best = lg2 - 5 * nerr ; }
	      if (
		  (
		   (lg2 > minLen + 2 && min > 1 && (nerr == 1 || (nerr > 1 && arr(err, 1, A_ERR).iShort > minLen + 2))) ) &&
		  ( nerr == 0 
		    || arr(err, 0, A_ERR).iShort > 1
		    || arr(err, 0, A_ERR).type ==	AMBIGUE						   
		    || arr(err, 0, A_ERR).type ==	ERREUR
		    )  /* avoid a jump on base 1, by scanning on pos we can find a non-jump */
		  )
	      { min = 1 ; pp = arrayMax(dna) - pos + 1 ; best = minLen ; }
	    }
	  if (min < seuil)
	    ok |= 2 ;
	  else
	    pp = 0 ;
	}

      posMax2 = (pos1 + posMax)/2 ;
      if (isCodage && (ok || (pp < posMax2 && min < 2) || (pp < posMax && min < 1)))
	{ 
	  ok++ ;
	  pp1 = pp ;
	}
      else 
	pp1 = pp = 0 ;


      /* search end of post_insert_region of vector at end of sequence */
      min = 100 ; nn = 0 ;
      vector = v2 ;
      lg = old = arrayMax (vector) ; /* to thrhshold at 10, should include 6 AAAAAA before specific sequence */
      if (0 && isCodage &&  lg > 15 /* tg2 not yet defined tg2 > 0 */ )
	{
	  nn = 0 ;
	  cp = arrp (vector, 0, char) ;
	  cq = arrp (dna, arrayMax(dna) -  pos, char) ;
	  while (nn < pos && ( *cp++ == *cq++)) nn++ ;
	  if (nn ==  pos) /* the vector exactly corresponds to the beginning of the sequence */
	    { 
	      min = 0 ;
	      cp = arrp (vector, 0, char) ;
	      while (*cp++ == A_) pos-- ;
	      if (pos >= 4)
		pp = pp2 = arrayMax(dna) -  pos ;
	    }
	}
      lg = arrayMax (vector) ; 
      if (0 && isCodage && min > 0 && lg > 20)
	for (pos = pp1 > 1 ? pp1 : pos1 ; pos + lg < arrayMax(dna) && min > 0 && nn < 10 ; pos++)
	  { 
	    lg2 = lg ; pos2 = pos + lg ;
	    err = trackErrors (vector, 0, &lg2, dna, pos, &pos2, &nn, err, 0) ;
	    nerr = arrayMax (err) ;
	    if (lg2 == lg && min > nerr) 
	      { min = nerr ; pp = pos ; }
	    if (nerr > 2*seuil && min < seuil)
	      break ;
	  }
      if (min < seuil) /*  && arrp (err,0, A_ERR)->iShort != 0) */
	{ 
	  if (pp2 > pp) pp2 = pp ;
	  ok++ ;
	}
    }

   

  if (!(force & 0x2) &&
      bsGetData (obj, _Vector_clipping, _Int, &pp1) &&
      bsGetData (obj, _bsRight, _Int, &nok) &&
      bsGetKey (obj, _bsRight, &vec))
   
 /* eat up linker stretches just downstream of the existing vector */   
  if (vec && pp1 > 1 && 
      (Vec = bsCreate (vec)))
    {
      char *mylinker = 0 ;
      
      if (bsGetData (Vec,str2tag("Followed_by_linker"), _Text, &mylinker) &&
	  strlen (mylinker) >= 3)
	{
	  pp = pp1 ;
	  nok = 0 ;
	  cp = arrp (dna, pp-1, char) ;
	  cq = mylinker ;
	  while (*cq && *cp == dnaEncodeChar [(int)(*cq)])
	    { cp++ ; cq++ ; }
	  if (!*cq) /* all linker was eaten */
	    { 
	      ok++ ;
	      isLinker = TRUE ;
	      strncpy (linker, mylinker, 255) ;
	      pp1 = pp + strlen (linker) ; /* plato */
	    } 
	}      
      bsDestroy (Vec) ;
    }
  
  /* eat up TTT stretches just downstream of the existing vector */
  if (pp1 && keyFindTag (vec,str2tag("Followed_by_polyT")))
    {
      pp = pp1 ;
      nok = 0 ;
      cp = arrp (dna, pp, char) ;
      while (pp--, pp >= 0 && *cp-- == T_) ;
      pp++ ; cp++ ; /* reposition just before the T */
      while (pp++, *++cp == T_) nok++ ;
      if (nok >= 6)
	{ 
	  ok++ ;
	  pp1 = pp + 1 ; /* plato */
	} 
    }      

  /* eat up AAA stretches just downstream of the vector */
  if (pp1 && keyFindTag (vec,str2tag("Followed_by_polyT")))
    {
      pp = pp1 ;
      nok = 0 ;
      cp = arrp (dna, pp, char) ;
      while (pp--, pp >= 0 && *cp-- == A_) ;
      pp++ ; cp++ ; /* reposition just before the A */
      while (pp++, *++cp == A_) nok++ ;
      if (nok >= 8)
	{ 
	  ok++ ;
	  pp1 = pp + 1 ; /* plato */
	} 
    }

  arrayDestroy (dna) ;
  arrayDestroy (dnaR) ;
  
  if (ok)
    { 
      c1 = -1 ; c2 = 1000000 ;
      bsGetData (obj, _Clipping, _Int, &c1) ;
      bsGetData (obj, _bsRight, _Int, &c2) ;
      if (c2 > pp2) c2 = pp2 ; 
      if (c1 < 0 || pp1 > c1 ) c1 = pp1 ;
 
      freeOutf ("Sequence %s\n"
		, freeprotect (name(seq))
		) ;
      freeOutf ("Vector_clipping %d %d %s %s\nClipping %d %d\n\n"
		, pp1, pp2,  freeprotect (name(vec))
		, isLinker ? linker : "" 
		, c1, c2
		) ;   
    }
  else if (0)
    { /* Prevent recursive search */
      freeOutf ("Sequence %s\nVector_Clipping\n\n"
		, freeprotect (name(seq))
		) ;
    }
  defCptForget (0, dnaKey) ;
  
 abort:
  arrayDestroy (err) ;
  bsDestroy (obj) ;
  return ok ;
}
 
/***********************/

static int doTrackVectorExact (KEY seq, Array v1, Array v2, int force, BOOL sl) 
{ 
  Array dna = 0, dnaR = 0, vector, units ;
  BSunit *uu ;
  OBJ obj, Tg ;
  char *cp, *cq ;
  BOOL isCodage ;
  KEY dnaKey = 0, vec, myvec = 0, tg ;
  int 
    ii, myerr1, myerr2, nok, myok1, myok2, nerr, best = 0,minLen = 30,
    old, c1, c2, ok = 0,ok2 = 0,  n, pos,
    pp1, pp2, min, lg, pos1, posMax
    ;
  Array err = 0 ;
  
  if (class(seq) == _VDNA && !dnaReClass(seq, &seq))
    return 0 ;
  
  obj = bsCreate (seq) ;
  if (!obj)
    return 0 ;
  
  if (!bsGetKey (obj, _DNA, &dnaKey))
    goto abort ;
  
  if (bsFindTag (obj, str2tag("Vector_Clipping")) && !force)
    goto searchSl ;
  
  if (0)
    { 
      /* int  tg1 = -1  :  beginning of current alignment to the genome */
      /* int  tgerr : nb of errors in first exon */

      if (bsGetKey (obj, _From_gene, &tg) &&
	  bsFindTag (obj, _Forward) &&
	  (Tg = bsCreate (tg)))
	{
	  units = arrayCreate (256, BSunit) ;
	  if (bsGetArray (Tg, _Assembled_from, units, 6))
	    for (ii = 0 ; ii < arrayMax (units) ; ii += 6)
	      {
		uu = arrp (units, ii, BSunit) ;
		if (uu[2].k == seq)
		  { /* tg1 = uu[3].i ;  tgerr = uu[5].i ; */ break ; }
	      }
	  bsDestroy (Tg) ;
	  arrayDestroy (units) ;
	}
    }
  
  dna = dnaGet (dnaKey) ;
  if (!dna)
    goto abort ;
  
  pp1 = 1 ; pp2 = arrayMax (dna) ;
  vec = 0 ;
  ok = 0 ;
  min = 100 ;
  myok1 = myok2 = myerr1 = myerr2 = 0 ;
  while (min && doGetVector (&vec, seq, v1, v2, &pos1, &posMax, &isCodage))
    { 
      /* search exact motif anywhere in the sequence */
      vector = v1 ;
      lg = old = arrayMax (vector) ; best = 0 ;
      minLen = isCodage ? 20 : 30 ;
      if (lg < minLen)
	continue ;
      for (pos = minLen ; min && pos < arrayMax (dna) ; pos++)
	{
	  n = lg ; nerr = nok = 0 ;
	  cp = arrp (vector, lg - 1, char) ;
	  cq = arrp (dna, pos, char) ;
	  while (n-- && nerr < 5)
	    {
	      if (*cp-- == *cq--) nok++ ;
	      else { nok++ ; nerr++ ; }
	      if (nok >= minLen && nerr * 5 < nok)
		{
		  if (nok - 5 * nerr > best)
		   {
		     ok |= 1 ;
		     min = 0 ;
		     best = nok - 5 * nerr ;
		     pp1 = pos + 2 ;
		     myvec = vec ;
		     myok1 = nok ;
		     myerr1 = nerr ; 
		   }
		}
	    }
	  if (nok >= minLen && nerr * 5 < nok && nok - 5 * nerr > best) /* the vector matches the sequence */
	    { 
	      ok |= 1 ;
	      min = 0 ;
	      pp1 = pos + 2 ;
	      myvec = vec ;
	      myok1 = nok ;
	      myerr1 = nerr ;
	    }
	}
    }

  vec = 0 ;
  while (min && doGetVector (&vec, seq, v1, v2, &pos1, &posMax, &isCodage))
    { 
      /* search around the expected position using trackErrors */
      int dx, ddx, lg2, pos2, lg1, nn, pp = 0;
      
      vector = v1 ;
      lg = old = arrayMax (vector) ; best = 0 ;
      minLen = isCodage ? 20 : 30 ;
      if (lg < minLen || !pos1)
	continue ;
      if (posMax <= 0) { posMax = pos1 + 35 ; }
      reverseComplement (vector) ;
      if (!dnaR)
	{
	  dnaR = dnaCopy (dna) ;
	  reverseComplement (dnaR) ;
	}
      pos1 = (long)arrayMax (dna) - pos1  + 1 ;
      posMax = (long)arrayMax (dna) - posMax  + 1 ;
      for (dx = 0 ; min > 0 && dx < 35 ; dx++)
	for (ddx = -1 ; min > 0 && ddx < 2 ; ddx += 2)
	  { 
	    pos = pos1 + ddx*dx ;
	    if (ddx == 1 && dx == 0) continue ;
	    if (pos <= lg ||  pos > posMax || pos + lg >= arrayMax(dna))
	      continue ;
	    lg2 = lg ; pos2 = pos + lg ; lg1 = 0 ;
	    err = trackErrors (vector, lg1, &lg2, dnaR, pos, &pos2, &nn, err, 0) ;
	    nerr = arrayMax (err) ;
	    if (
		((lg2 == lg && min > nerr) || (lg2 > minLen && lg2 - 5 * nerr > best)) &&	
		( nerr == 0  
		  || arr(err, 0, A_ERR).iShort > 1
		  || arr(err, 0, A_ERR).type ==	AMBIGUE						   
		  || arr(err, 0, A_ERR).type ==	ERREUR
		  ) /* avoid a jump on base 1, by scanning on pos we can find a non-jump */
		)
	      { min = nerr ; pp = arrayMax(dna) - pos + 1 ; best = lg2 - 5 * nerr ; }
	  }
      if (best > 0)
	{
	  ok |= 1 ;
	  pp1 = pp ;
	  myvec = vec ;
	  myok1 = lg ;
	  myerr1 = min ;
	  min = 0 ; /* no need to look for another vector */
	}
    }

  vec = 0 ;
  min = 100 ;
  while (min && doGetVector (&vec, seq, v1, v2, &pos1, &posMax, &isCodage) && v2)
    { 
      /* search exact motif de sortie  anywhere in the sequence */
      vector = v2 ;
      lg = old = arrayMax (vector) ;
      minLen = isCodage ? 20 : 30 ;
      if (lg < minLen)
	continue ;
      for (pos = 0 ; min && pos + lg < arrayMax (dna) ; pos++)
	{
	  n = lg ; nerr = nok = 0 ; best = 0 ;
	  cp = arrp (vector, 0, char) ;
	  cq = arrp (dna, pos, char) ;
	  while (n-- && nerr < 5)
	    {
	      if (*cp++ == *cq++) nok++ ;
	      else { nok++ ; nerr++ ; }
	      if (nok >= minLen && nerr * 5 < nok)
		{
		  if (nok - 5 * nerr > best)
		   {
		     ok |= 2 ;
		     min = 0 ;
		     best = nok - 5 * nerr ;
		     pp2 = pos ;
		     myvec = vec ;
		     myok2 = nok ;
		     myerr2 = nerr ; 
		   }

		}
	    }
	  if (nok >= minLen && nerr * 5 < nok && nok - 5 * nerr > best) /* the vector matches the sequence */
	    { 
	      ok |= 2 ;
	      min = 0 ;
	      pp2 = pos ;
	      myok2 = nok ;
	      myvec = vec ;
	      myerr2 = nerr ;
	    }
	}
    }
  if (ok)
    {
      char *ok_text[] = { "NULL", "Exact_start", "Exact_stop", "Exact"} ;

      c1 = -1 ; c2 = 1000000 ;
      bsGetData (obj, _Clipping, _Int, &c1) ;
      bsGetData (obj, _bsRight, _Int, &c2) ;
      if (c2 > pp2) c2 = pp2 ; 
      if (c1 < 0 || pp1 > c1 ) c1 = pp1 ;

      freeOutf ("Sequence %s\n"
		, freeprotect (name(seq))
		) ;
      if (myerr1 + myerr2 == 0)
	freeOutf ("Vector_clipping %d %d %s %s\nClipping %d %d\n\n"
		  , pp1, pp2,  freeprotect (name(myvec)), ok_text[ok]
		  , c1, c2
		  ) ;
      else
	freeOutf ("Vector_clipping %d %d %s \"Match_%d/%d__%d/%d\"\nClipping %d %d\n\n"
		  , pp1, pp2,  freeprotect (name(myvec))
		  , myok1 - myerr1, myok1, myok2 - myerr2, myok2
		  , c1, c2
		  ) ;      
    }
  else if (0)
    {
      freeOutf ("Sequence %s\n"
		, freeprotect (name(seq))
		) ;
      freeOutf ("Vector_Clipping\n\n") ; /* Prevent recursive search */
    }
  
 searchSl:
  if (!dna)
    dna = dnaGet (dnaKey) ;
  if (!dna)
    goto abort ;
 
  pp1 = 1 ; pp2 = arrayMax (dna) ;
  vec = 0 ;
  ok2 = 0 ; myerr1 = myok1 = 0 ;
  min = 100 ;
  while (min && doGetSl (&vec, v1))
    { 
      /* search exact motif anywhere in the sequence */
      vector = v1 ;
      lg = old = arrayMax (vector) ;
      if (lg < 20)
	continue ;
      for (pos = 20 ; min && pos < arrayMax (dna) ; pos++)
	{
	  n = lg ; nerr = nok = 0 ;
	  cp = arrp (vector, lg - 1, char) ;
	  cq = arrp (dna, pos, char) ;
	  while (n-- && nerr < 3)
	    {
	      if (*cp-- == *cq--) nok++ ;
	      else { nok++ ; nerr++ ; }
	    }
	  if (nok >= 20 && nerr < 3 && nerr < min) /* the vector matches the sequence */
	    { 
	      ok2 |= 1 ;
	      min = nerr ;
	      pp1 = pos + 1 ; /* first base of the sl */
	      myvec = vec ; 
	      myok1 = nok ;
	      myerr1 = nerr ;
	    }
	}
    }

  if (ok2)
    {
      char *ok_text[] = { "NULL", "Exact" } ;

      freeOutf ("Sequence %s\n"
		, freeprotect (name(seq))
		) ;
      if (myerr1 == 0)
	freeOutf ("Transpliced_to2 %s %d %s\n\n"
		  , freeprotect (name(myvec)), pp1,  ok_text[ok2]
		  ) ;
      else
	freeOutf ("Transpliced_to2 %s %d \"Match_%d/%d\"\n\n"
		  , freeprotect (name(myvec)), pp1  
		  , myok1 - myerr1, myok1
		  ) ;

    }

  defCptForget (0, dnaKey) ;

 abort:
  arrayDestroy (err) ;
  arrayDestroy (dnaR) ;
  bsDestroy (obj) ;
  return ok + 2 * ok2 ;
} /* doTrackVectorNew */
 
/***********************/

int trackVector (KEYSET ks, KEY key, BOOL force)
{ 
  int ok, force1, myforce, n = 0, i = ks ? keySetMax (ks) : 0 ; 
  KEY *kp = i ? arrp(ks, 0, KEY) - 1 : 0 ;
  Array 
    v1 = arrayCreate (20, char),
    v2 = arrayCreate (20, char) ;

  if (!ks && key) { i =  1; kp = &key ; kp-- ; }
  force1 = force ? 0x7 : 0x6 ;
  while (kp++, i--)
    {
      ok = 0 ; myforce = force1 ;
      switch (doTrackVectorExact (*kp, v1, v2, myforce, TRUE))
	{
	case 1: ok = 1 ; myforce  = force ? 0x5 : 0x4 ; break ; /* vector found */
	case 2: ok = 2 ; myforce  = force ? 0x3 : 0x2 ; break ; /* sl found found */
	case 3: ok = 3 ; myforce = 0 ; break ; /* both found */
	}
      if (myforce)
	switch (doTrackVector (*kp, v1, v2, myforce))
	  {
	  case 1: ok |= 1 ; myforce  = force ? 0x5 : 0x4 ; break ; /* vector found */
	  case 2: ok |= 2 ; myforce  = force ? 0x3 : 0x2 ; break ; /* sl found found */
	  case 3: ok |= 3 ; myforce = 0 ; break ; /* both found */
	  }
      switch (ok)
	{
	  case 1: n++ ; break ;
	  case 2: n++ ; break ;
	  case 3: n++ ; break ;
	}
    }

  arrayDestroy (v1) ;
  arrayDestroy (v2) ;
  return n ;
}
 
/***********************/
/* locate the largest zone without 30% N over next 20 N */
/* we use the mountain search algo */
static int doGetBestRange (Array dna , int *x1p, int *x2p)
{
  int 
    bestRange = 0,
    slope = 0,
    delta = 30, /* i.e. 10 non compensated N */
    x1,  /* end previous zone */
    y1, y2, /* tmp end of current zone */
    z1, z2 ; /* current pos */
  char *cp ;

  *x1p = *x2p = 0 ;
  if (!dna || ! arrayMax(dna))
    return 0 ;
  x1 = - 1 ; y1 = y2 = z1 = z2 = 0 ; slope = -1 ;
  for (z1 = 0, cp = arrp (dna, 0, char) ; z1 < arrayMax(dna) ; z1++, cp++)
    {
      switch (*cp)
	{
	case A_: case T_: case G_: case C_: z2 -= 1 ; break ;
	default: z2 += 3 ; break ; /* we wish 30% N or ambiguous letter */
	}
      if (slope == -1 && z2 < y2) /* z becomes our tmp mini */
	{ y1 = z1 ; y2 = z2 ; }
      else if (slope == 1 && z2 > y2) /* z becomes our tmp maxi */
	{ y1 = z1 ; y2 = z2 ; }
      else if (slope == -1 && z2 > y2 + delta) /* y is a true minimum */
	{ 
	  if (y1 - x1 > bestRange) /* best range with few N */
	    { bestRange = y1 - x1 ; *x1p = x1 + 2 ; *x2p = y1 + 1 ; }
	  x1 = y1 ; /* x2 = y2 ; */ slope = 1 ;
	  y1 = z1 ; y2 = z2 ;
	}
      else if (slope == 1 && z2 < y2 - delta) /* y is a true maximum */
	{ 
	  x1 = y1 ; /* x2 = y2 ; */ slope = -1 ;
	  y1 = z1 ; y2 = z2 ;
	}
    }
  if (slope == -1 && y1 - x1 > bestRange) /* best range with few N */
    { bestRange = y1 - x1 ; *x1p = x1 + 2 ; *x2p = y1 + 1 ; }

  return bestRange ;
} /* doGetBestRange */

/***********************/

static int doTrackBadQuality (KEY seq, int nMin, int delta)
{
  BOOL isBad = FALSE ;
  Array dna = dnaGet (seq) ;
  int i, j, max, na, nc, nt, ng, x1, x2 ;
  char *cp ;
  int clipTop, clipEnd ;
  int nn[100] ; /* rolling position of last NMAX=6 + 1 n */

  max = dna ? arrayMax(dna) : 0 ;
  clipTop = 1 ; clipEnd = max ;
  doGetBestRange (dna, &x1, &x2) ;
  if (x1 > 1 || x2 < max)
    {
      OBJ Seq = bsUpdate (seq) ;
      KEY vector = 0 ;
      int ii = 0 ;

      if (Seq)
	{
	  if (bsGetData (Seq, _Vector_Clipping, _Int, &clipTop) &&
	      bsGetData (Seq, _bsRight, _Int, &clipEnd))
	    bsGetKey (Seq, _bsRight, &vector) ;
	  
	  if (clipTop < x1)
	    { ii = 1 ; clipTop = x1 ; }
	  if (clipEnd > x2)
	    { ii = 1 ; clipEnd = x2 ; }
	  if (ii)
	    {
	      bsAddData (Seq, _Vector_Clipping, _Int, &clipTop) ;
	      bsAddData (Seq, _bsRight, _Int, &clipEnd) ;
	      if (vector)
		bsAddKey (Seq, _bsRight, vector) ;
	      printf ("Clipping x1=%d x2=%d len=%d est=%s\n", clipTop, clipEnd, max, name(seq)) ;
	    }	  
	  bsSave (Seq) ;
	}
    }
  if (max > clipEnd)
    max = clipEnd ;
  /* with cp on 7th n, there is only 6 n up tocp ! */
  j = nMin + 1 ; while (j--) nn[j] = clipTop - 1 ; 

  if (max > 600) max = 600 ;
  if (dna)
    {
      isBad = TRUE ;
      na = nt = nc = ng = 0 ;

      if (max > clipTop)
	for (i = clipTop, cp = arrp (dna, clipTop, char) ; 
	     isBad && i <= max ; cp++, i++)
	  {
	    if (i < max)  /* we enter the loop once after the last letter */
	      switch (*cp)
		{  /* do not trust repeats longer than 4 of a kind */
		case A_: na++ ; nt = nc = ng = 0 ; if (na < 5) continue ; else break ;
		case T_: nt++ ; na = nc = ng = 0 ; if (nt < 5) continue ; else break ;
		case C_: nc++ ; na = nt = ng = 0 ; if (nc < 5) continue ; else break ;
		case G_: ng++ ; na = nt = nc = 0 ; if (ng < 5) continue ; else break ;
		}
	    if (i - nn[nMin] > delta)  /* rock  */
	      isBad = FALSE ;
	    else                         /* roll */
	      { j =  nMin ;  while (j--) nn[j + 1] = nn[j] ;  nn[0] = i ; }
	  }
    }
  
  arrayDestroy (dna) ;
  if (isBad)
    {
      OBJ Seq = bsUpdate (seq) ;
      KEY tag = str2tag ("Bad_quality") ;
      
      cp = messprintf("At least %dn inbest %d bp", nMin, delta) ;
      if (Seq)
	{
	  bsAddTag (Seq, tag) ;
	  if (bsFindTag (Seq, tag))
	    bsAddData (Seq, _bsRight, _Text, cp) ;
	  bsSave (Seq) ;
	}
    }

  return isBad ;
}

/***********************/

int trackBadQuality (KEYSET ks)
{ 
  int n = 0, i, nn = 6, delta = 100 ;
  char *cp ;

  if ((cp = freeword()) && !strcmp (cp, "-n") && freeint (&nn) && nn > 0) ;
  else nn = 6 ;
  if ((cp = freeword()) && !strcmp (cp, "-d") && freeint (&delta) && delta > 50) ;
  else delta = 100 ;

  if (ks) 
    for (i = 0 ; i < keySetMax(ks) ; i++)
      if (doTrackBadQuality (keySet (ks, i), nn, delta))
	n++ ;
  return n ;
}
 
/***********************/
  /* search for at least 8 contiguous matches */
BOOL baseCallTrackBout (Array dna1, Array dna2, int *x1p, int *x2p, int *sensp)
{ int i, j, x, p, mn, y1, y2 ;
  int max1 = arrayMax (dna1), max2 = arrayMax (dna2) ;
  char *cp, *cq, c ;
  
  if (max1 < 152 || max2 < 152) return FALSE ;
/* esaie 1 2 dans l'ordre */
  mn = y1 = y2 = 0 ; 
  for (x = max1 - 150 ; x < max1 - 8 ; x++)
    { p = 0 ; i = x ; j = 0 ;
      cp = arrp (dna1, i, char) ;
      cq = arrp (dna2, j, char) ;
      
      while (i++ < max1 && j++ < max2)
	{ c = *cp++ & *cq++ & 0x0f ;
	  if (c == N_) continue ;
	  else if (c) p++ ;
	  else 
	    { if (p > mn) 
		{ mn = p ; y1 = x ; y2 = 0 ; *sensp = 1 ; }
	      p = 0 ;
	    }
	}
      if (p > mn) 
	{ mn = p ; y1 = x ; y2 = 0 ; *sensp = 1 ; }
    }
	      
/* esaie 2 1 dans l'ordre */

  for (x = max2 - 150 ; x < max2 - 8 ; x++)
    { p = 0 ; i = 0 ; j = x ;
      cp = arrp (dna1, i, char) ;
      cq = arrp (dna2, j, char) ;
      
      while (i++ < max1 && j++ < max2)
	{ c = *cp++ & *cq++ & 0x0f ;
	  if (c == N_) continue ;
	  else if (c) p++ ;
	  else 
	    { if (p > mn) 
		{ mn = p ; y1 = 0 ; y2 = x ; *sensp = 1 ; }
	      p = 0 ;
	    }
	}
      if (p > mn) 
	{ mn = p ; y1 = 0 ; y2 = x ; *sensp = 1 ; }
    }
	      
/* esaie 1 anti-2 dans l'ordre */

  for (x = max1 - 150 ; x < max1 - 8 ; x++)
    { p = 0 ; i = x ; j = max2 - 1 ;
      cp = arrp (dna1, i, char) ;
      cq = arrp (dna2, j, char) ;
      
      while (i++ < max1 && j--)
	{ c = *cp++ & complementBase(*cq-- & 0x0f) ;
	  if (c == N_) continue ;
	  else if (c) p++ ;
	  else 
	    { if (p > mn) 
		{ mn = p ; y1 = x ; y2 = max2 - 1 ; *sensp = -1 ; }
	      p = 0 ;
	    }
	}
      if (p > mn) 
	{ mn = p ; y1 = x ; y2 = max2 - 1 ; *sensp = -1 ; }
    }
	      
/* esaie anti-2 1 dans l'ordre */

  for (x = 150 ; x > 8 ; x--)
    { p = 0 ; i = 0 ; j = x ;
      cp = arrp (dna1, i, char) ;
      cq = arrp (dna2, j, char) ;
      
      while (i++ < max1 && j--)
	{ c = *cp++ & complementBase (*cq-- & 0x0f) ;
	  if (c == N_) continue ;
	  else if (c) p++ ;
	  else 
	    { if (p > mn) 
		{ mn = p ; y1 = 0 ; y2 = x ; *sensp = -1 ; }
	      p = 0 ;
	    }
	}
      if (p > mn) 
	{ mn = p ; y1 = 0 ; y2 = x ; *sensp = -1 ; }
    }

  if (mn >= 12)
    { *x1p = y1 ; *x2p = y2 ; return TRUE ; }
  return FALSE ;
}

/***********************/
/***********************/
/*********************************************************************/
/*********************************************************************/

static void fixGenePolyA (KEY est, OBJ Gene, int a2g, int x2g, int z1g, int z2g, Array dnaCosmid, Array dnaEst)
{
  KEY cosmid, _na = 0 ;
  OBJ Est = 0 ;
  int ia, ja, nta = 0, pa, ppa, a1, a2, x1, x2, z1, z2 ;
  int iapure, ntapure, papure, dummy, clipStart, clipEnd, manual_pa ;
  Array dna = 0 ;
  BOOL isDown = TRUE ;
  char *cp, *cq ;

  lexword2key ("Number_of_terminal_A", &_na, _VSystem) ;
  if (!bsIsTagInClass (_VSequence, _na))
    _na = 0 ;

  if (est)
    {
      clipEnd = arrayMax (dnaEst) ;
      if (keyFindTag (est, _Vector_clipping) &&
	  (Est = bsCreate (est)))
	{
	  if (bsGetData (Est, _Vector_clipping, _Int, &dummy))
	    bsGetData (Est, _bsRight, _Int, &clipEnd) ;
	  bsDestroy (Est) ;
	}
      if (clipEnd > arrayMax (dnaEst))
	clipEnd = arrayMax (dnaEst) ;
      clipStart = 1 ;
      /* search for tail of nnnnnn */
      if (1)
	{
	  int i1 = clipEnd, iN = 10 ;

	  cp = arrp (dnaEst, clipEnd - 1, char) ;
	  while (iN > 0 && i1-- >0)
	    switch (*cp--)
	      {
	      case N_: if (iN < 10) iN++ ; break ;
	      default: iN-- ; break ;
	      }
	  if (clipEnd - i1 > 18)
	    clipEnd = i1 + 10 ;
	  if ((Est = bsUpdate (est)))
	    {
	      bsGetData (Est, _Vector_clipping, _Int, &clipStart) ;
	      bsAddData (Est, _Vector_clipping, _Int, &clipStart) ;
	      bsAddData (Est, _bsRight, _Int, &clipEnd) ;
	      bsSave (Est) ;
	    }
	}
      ntapure = -1 ;
      if (_na && keyFindTag (est, _Forward))
	{
	  pa = clipEnd ;
	  cp = arrp (dnaEst, pa - 1, char) + 1 ;
	  ia = iapure = papure = nta = 0 ; ntapure = -1 ; ja = 6 ; ppa = 0 ;
	  while (cp--, pa-- > 0)
	    {
	      if (*cp == A_)
		{ 
		  ia++ ; if (ntapure == -1) { iapure++ ; papure = pa ; }
		  if (ja < 6) 
		    ja++ ; 
		  else
		    { 
		      nta = ia ;
		      ppa = pa ; /* last trustable AA */
		    }
		}
	      else if (*cp & A_)
		{ if (ntapure == -1) ntapure = iapure ; }
	      else 
		{ ja-=2 ; if (ntapure == -1) ntapure = iapure ; }
	      if (ja <= 0)
		break ;
	    }
	  if (nta > 10)
	    pa = ppa ;
	  else if (ntapure > 5)
	    pa = papure ;
	  else
	    nta = pa = 0 ;
	  if (nta > 5 && pa > 200 && 
	      (Est = bsUpdate (est)))
	    {
	      bsAddData (Est, _na, _Int, &nta) ;
	      if (bsGetKey (Est, str2tag("Best_cosmid_alignment"), &cosmid) &&
		  bsGetData (Est, _bsRight, _Int, &a1) &&
		  bsGetData (Est, _bsRight, _Int, &a2) &&
		  bsGetData (Est, _bsRight, _Int, &x1) &&
		  bsGetData (Est, _bsRight, _Int, &x2) &&
		  pa + ia < x2 - 30) ;
	      else
		bsAddData (Est, _PolyA_after_base, _Int, &pa) ;
	      bsSave (Est) ;
	    }	
	}

      if (_na && keyFindTag (est, _Forward) 
	  && ! keyFindTag (est, _mForward)  
	  && ! keyFindTag (est, _gt_ag) 
	  && ! keyFindTag (est, _gc_ag) 
	  && ! keyFindTag (est, _PolyA_after_base)
	  && ! keyFindTag (est, _aForward)) /* is there an initila polyT ? */
	{
	  pa = clipStart ;
	  cp = arrp (dnaEst, pa - 1, char) ;
	  ia = nta = 0 ; ja = 3 ; ppa = 0 ;
	  while (cp++, pa++ > 0)
	    {
	      if (*cp == T_)
		{ 
		  ia++ ;
		  if (ja < 3) 
		    ja++ ; 
		  else
		    { 
		      nta = ia ;
		      ppa = pa ; /* last trustable AA */
		    }
		}
	      else if (*cp & T_) ;
	      else 
		ja-- ;
	      if (ja <= 0)
		break ;
	    }
	  if (nta > 12)
	    pa = ppa ;
	  else
	    nta = pa = 0 ;
	  if (nta > 12 && pa < clipStart + nta + 20 && 
	      (Est = bsUpdate (est)))
	    {
	      bsAddData (Est, _na, _Int, &nta) ;
	      if (bsGetKey (Est, str2tag("Best_cosmid_alignment"), &cosmid) &&
		  bsGetData (Est, _bsRight, _Int, &a1) &&
		  bsGetData (Est, _bsRight, _Int, &a2) &&
		  bsGetData (Est, _bsRight, _Int, &x1) &&
		  bsGetData (Est, _bsRight, _Int, &x2) &&
		  pa - ia > x1 + 5) ; 
	      else
		{
		  bsAddTag (Est, _Reverse) ;
		  bsAddData (Est, _PolyA_after_base, _Int, &pa) ;
		}
	      bsSave (Est) ;
	    }	
	}

      if (_na && keyFindTag (est, _Reverse))
	{
	  pa = 0 ;
	  cp = arrp (dnaEst, 0, char) - 1 ;
	  ia = nta = 0 ; ja = 3 ; ppa = 0 ;
	  while (cp++, pa++ < clipEnd)
	    {
	      if (*cp == T_)
		{ 
		  ia++ ;
		  if (ja < 3) 
		    ja++ ; 
		  else
		    { 
		      nta = ia ;
		      ppa = pa ; /* last trustable AA */
		    }
		}
	      else if (*cp & T_) ;
	      else 
		ja-- ;
	      if (ja <= 0)
		break ;
	    }
	  if (nta > 8)
	    pa = ppa ;
	  else
	    nta = pa = 0 ;

	  if (nta >= 8 &&
	      (Est = bsUpdate (est)))
	    {
	      bsAddData (Est, _na, _Int, &nta) ;
	      if (bsGetKey (Est, str2tag("Best_cosmid_alignment"), &cosmid) &&
		  bsGetData (Est, _bsRight, _Int, &a1) &&
		  bsGetData (Est, _bsRight, _Int, &a2) &&
		  bsGetData (Est, _bsRight, _Int, &x1) &&
		  bsGetData (Est, _bsRight, _Int, &x2) &&
		  pa - ia > x1 + 30) ;
	      else
		{
		  pa++ ;
		  bsAddData (Est, _PolyA_after_base, _Int, &pa) ;
		}
	      bsSave (Est) ;
	    }	
	}
      pa = manual_pa = 0 ;
      if ((Est = bsUpdate (est)))
	{
	  isDown = bsFindTag (Est, _Forward) ;
	  if (bsGetData (Est, _Manual_polyA, _Int, &manual_pa))
	    bsAddData (Est, _PolyA_after_base, _Int, &manual_pa) ;
	  bsGetData (Est, _PolyA_after_base, _Int, &pa) ;
	  if ((bsFindTag (Est, _Forward) && pa < 200) ||
	      (bsFindTag (Est, _Reverse) && pa > 100))
	    pa = 0 ;
	  bsSave (Est) ;
	}

      if (!Gene || !dnaEst || manual_pa || ntapure >= 12) ;      
      else 
	{
	  a1 = a2 = x1 = x2 = 0 ;
	  a2 = a2g ; x2 = x2g ;
	  if (pa && ((isDown && pa > x2 + 5) ||  (!isDown && pa < x2 - 5)))
	    ;  /* not correctly aligned, do not touch the polya */
	  else if (x2 && dnaCosmid)
	    {
	      z1 = z1g ; z2 = z2g ;
	      dna = dnaCosmid ;
	      if (dna)
		{
		  nta = 0 ;
		  a2 += z1 - 1 ;
		  if (a2 >= 0 && a2 < arrayMax(dna) &&
		      x2 > 1 && x2 + 1 < arrayMax(dna)
		      ) /* compute nta */
		    {
		      if (isDown)
			{
			  cp = arrp (dna, a2, char) ; /* on ., first base post alignment */
			  cq = arrp (dnaEst, x2, char) ;
			  cp-- ; cq-- ; /* last aligned base */
			  while (*cq == A_) { x2-- ; a2-- ; cp-- ; cq-- ; }
			  cq++ ; x1 = x2 ; /* first A or . of xxxAA. stretch */
			  ia = 0 ; ja = 2 ; nta = 0 ;
			  while (*cq && ja > 0 &&
				 x2 < arrayMax(dnaEst)) 
			    { 
			      x2++ ; 
			      if (*cq == A_) { nta++ ; if (ja < 1) ja++ ; }
			      else ja-- ;
			      cq++ ;
			    }
			}
		      else
			{
			  cp = arrp (dna, a2, char) ; /* on ., first base post alignment */
			  cq = arrp (dnaEst, x2-2, char) ;
			  cp-- ; cq++ ;  /* last aligned base */
			  while (*cq == T_) { x2++ ; a2-- ; cp-- ; cq++ ; }
			  cq-- ; x1 = x2 ;  /* first base post A */
			  ia = 0 ; ja = 2 ; nta = 0 ;
			  while (*cq && ja > 0 &&
				 x2 > 1)
			    { 
			      x2-- ;
			      if (*cq == T_) { nta++ ; if (ja < 1) ja++ ; }
			      else ja-- ;
			      cq-- ;
			    }
			}
		    }
	
		  if (1)
		    {
		      if (nta > 9 && x1 != pa)
			printf ("Moving Gene %s Est %s, %d to %d z2=%d\n", 
				bsName(Gene), name (est), pa, x1, z2) ;
		      if ((Est = bsUpdate (est)))
			{
			  if (nta > 9)
			    {
			      bsAddData (Est, _na, _Int, &nta) ;
			      bsAddData (Est, _PolyA_after_base, _Int, &x1) ;
			    }
			  else
			    { 
			      if (bsFindTag (Est, _PolyA_after_base))
				bsRemove (Est) ;
			      if (bsFindTag (Est, _na))
				bsRemove (Est) ;
			    }
			  bsSave (Est) ;
			}
		    }
		}
	    }
	}
    } 

} /* fixGenePolyA */

/*********************************************************************/

void fixPolyA (KEYSET ks0)
{
  KEYSET ests = query (ks0, "Is_read") ;
  KEYSET estsNoGene = query (ests, "! from_gene") ;
  KEYSET ks1, genes = query (ests, ">From_gene") ;
  int ii, jj, ia, a2, x1, x2, z1, z2 ;
  Array aa = 0, dna = 0, dnaEst = 0 ;
  KEY est, gene, cosmid ;
  OBJ Gene = 0 ;
  BSunit *uu ;

  ks1 = keySetCreate () ;
  aa = arrayCreate (300, BSunit) ; 
  for (ii = 0 ; ii < keySetMax (estsNoGene) ; ii++)
    {
      est = keySet (estsNoGene, ii) ;
      dnaEst = dnaGet (est) ;
      if (dnaEst && arrayMax(dnaEst) > 10)
	fixGenePolyA (est, 0, 0, 0, 0, 0, 0, dnaEst) ;
      arrayDestroy (dnaEst) ;
    }
  for (jj = 0 ; jj < keySetMax (genes) ; jj++)
    {
      gene = keySet (genes, jj) ;
      if ((Gene = bsCreate (gene)))
	{
	  z1 = z2 = 0 ; dna = 0 ; cosmid = 0 ;
	  if (bsGetKey (Gene, _Genomic_sequence, &cosmid) &&
	      bsGetArray (Gene, _Covers, aa, 4) &&
	      arrayMax(aa) > 3)
	    {
	      uu = arrp (aa, 0, BSunit) ;
	      z1 = uu[2].i ; z2 = uu[3].i ;
	      if (cosmid && (dna = dnaGet (cosmid)))
		{
		  if (z1 > z2)
		    {
		      int i = arrayMax(dna) ;
		      
		      reverseComplement (dna) ;
		      z1 = i - z1 + 1 ; z2 = i - z2 + 1 ;
		    }
		}
	    }
	  if (!dna) { cosmid = 0 ; z1 = z2 = 0 ; }
	  bsGetArray (Gene, _Assembled_from, aa, 5) ;
	  a2 = x1 = x2 = 0 ;
	  /* find last base of alignment of each est */
	  for (ia = (long)arrayMax (aa) - 5 ; ia > 0 ; ia -= 5)
	    {
	      int ii ;

	      uu = arrp (aa, ia, BSunit) ;
	      est = uu[2].k ;  a2 = uu[1].i ; x2 = uu[4].i ;
	      if (keySetFind (ests, est, &ii))
		{
		  if (!keySet (ks1, ii))
		    {
		      keySet (ks1, ii) = 1 ;
		      dnaEst = dnaGet (est) ;
		      if (dnaEst && arrayMax(dnaEst) > 10)
			fixGenePolyA (est, Gene, a2, x2, z1, z2, dna, dnaEst) ;
		      arrayDestroy (dnaEst) ;
		    }
		}
	    }
	  bsDestroy (Gene) ;
	  arrayDestroy (dna) ;
	}
    }
  arrayDestroy (aa) ; 
  keySetDestroy (ks1) ;
  keySetDestroy (genes) ;
  keySetDestroy (ests) ;
  keySetDestroy (estsNoGene) ;
} /* fixPolyA */

/*********************************************************************/
/*********************************************************************/

static void fixGeneVector (KEY est, OBJ Gene, Array aag, int z1g, int z2g, Array dnaCosmid, Array dnaEst)
{
  KEY vector ;
  OBJ Est = 0 ;
  int ia, a1, a2, x1, x2, z1, clipTop, clipEnd ;
  BSunit *uu ;
  Array dna = 0 ;
  BOOL isDown = TRUE ;
  char *cp, *cq ;

  if (est && Gene && keyFindTag (est, _Forward) 
      && (Est = bsCreate (est))) /* ensure isDown in calling function */
    {
      clipTop = 1 ; clipEnd = -1 ; vector = 0 ;
      if (keyFindTag (est, _Vector_clipping))
	{
	  if (bsGetData (Est, _Vector_clipping, _Int, &clipTop) &&
	      bsGetData (Est, _bsRight, _Int, &clipEnd))
	    bsGetKey (Est, _bsRight, &vector) ;
	}
      if (clipTop < 1)
	clipTop = 1 ;

      bsDestroy (Est) ;

      if (clipTop > 1 && Gene && aag)
	{
	  a1 = a2 = x1 = x2 = 0 ;
	  for (ia = 0 ; ia + 4 < arrayMax (aag) ; ia += 5)
	    {
	      uu = arrp (aag, ia, BSunit) ;
	      if (est == uu[2].k && /* find the alignment on the genome of the first base */
		  clipTop >= uu[3].i && clipTop <= uu[4].i)
		{
		  a1 = uu[0].i ; a2 = uu[1].i ; 
		  x1 = uu[3].i ; x2 = uu[4].i ; 
		  break ;
		}
	    }
	  if (x1 && dnaCosmid)
	    {
	      z1 = z1g ; dna = dnaCosmid ;
	      if (dna)
		{
		  cq = 0 ; a1 += z1 - 1 ;
		  a1 += clipTop - x1 - 2 ; x1 = clipTop - 2 ; /* corresponding first bp inside the vector */
		  if (a1 >= 0 && a1 < arrayMax(dna))
		    {
		      if (clipEnd == -1 || clipEnd > arrayMax (dnaEst))
			clipEnd = arrayMax (dnaEst) ;

		      cp = arrp (dnaEst, clipTop - 2, char) ;
		      cq = arrp (dna, a1, char) ; /* first aligned base left of clipTop */
		      if (isDown)
			{
			  ia = 0 ;
			  while (*cp == *cq  &&
				 x1 > 0 &&
				 a1 > 0)
			    {ia++ ; a1-- ; x1-- ; cp-- ; cq-- ; }
			 
			  if (ia)
			    {
			      printf ("Moving Gene %s Est %s, %d to %d\n", 
				      bsName(Gene), name (est), clipTop, clipTop - ia) ;
			      if ((Est = bsUpdate (est)))
				{
				  if (bsFindTag (Est,  _Vector_clipping))
				    bsRemove (Est) ;
				  clipTop -= ia ;
				  if (clipTop < 4)
				    clipTop = 1 ;
				  if (clipTop == 1 || ia > 6)
				    vector = 0 ;
				  bsAddData (Est, _Vector_clipping, _Int, &clipTop) ;
				  bsAddData (Est, _bsRight, _Int, &clipEnd) ;
				  if (vector)
				    bsAddKey (Est,  _bsRight, vector) ;
				  bsSave (Est) ;
				}
			    }
			}
		    }
		}
	    }
	}
      /* search for initial polyA */
      if (keyFindTag (est, str2tag ("gt_ag")) || 
	  keyFindTag (est, str2tag ("gc_ag"))
	  )
	{ /* but only if we are absolutely sure of the orientation */
	  cp = arrp (dnaEst, 0, char) ;
	  ia = 0 ;
	  while (*cp++ == A_) { ia++ ; }
	  if (ia > 6 &&
	      (Est || (Est = bsUpdate (est))))
	    {
	      x1 = 1 ;
	      bsAddData (Est, str2tag ("Initial_polyA"), _Int, &x1) ;
	      bsAddData (Est, _bsRight, _Int, &ia) ;
	    }
	  if (clipTop > 1 &&
	      (Est || (Est = bsUpdate (est))))
	    {
	      cp = arrp (dnaEst, 0, char) ;
	      ia = 0 ;
	      while (*cp++ == A_) { ia++ ; }
	      if (ia > 6)
		{
		  x1 = clipTop ;
		  bsAddData (Est, str2tag ("Initial_polyA"), _Int, &x1) ;
		  bsAddData (Est, _bsRight, _Int, &ia) ;
		}
	    }
	  cp = arrp (dnaEst, 0, char) ;
	  ia = 0 ;
	  while (*cp++ == T_) { ia++ ; }
	  if (ia > 6 &&
	      (Est || (Est = bsUpdate (est))))
	    {
	      x1 = 1 ;
	      bsAddData (Est, str2tag ("Initial_polyT"), _Int, &x1) ;
	      bsAddData (Est, _bsRight, _Int, &ia) ;
	    }
	  if (clipTop > 1 &&
	      (Est || (Est = bsUpdate (est))))
	    {
	      cp = arrp (dnaEst, 0, char) ;
	      ia = 0 ;
	      while (*cp++ == T_) { ia++ ; }
	      if (ia > 6)
		{
		  x1 = clipTop ;
		  bsAddData (Est, str2tag ("Initial_polyT"), _Int, &x1) ;
		  bsAddData (Est, _bsRight, _Int, &ia) ;
		}
	    }
	  bsSave (Est) ;
	}
    }
} /* fixVector */

/*********************************************************************/

void fixVector (KEYSET ks0)
{
  KEYSET ests = query (ks0, "Is_read && Forward") ;
  KEYSET estsNoGene = query (ests, "! from_gene") ;
  KEYSET ks1, genes = query (ests, ">From_gene") ;
  int ii, jj, ia, z1, z2 ;
  Array aa = 0, dna = 0, dnaEst = 0 ;
  KEY est, gene, cosmid ;
  OBJ Gene = 0 ;
  BSunit *uu ;

  ks1 = keySetCreate () ;
  aa = arrayCreate (300, BSunit) ; 
  for (ii = 0 ; ii < keySetMax (estsNoGene) ; ii++)
    {
      est = keySet (estsNoGene, ii) ;
      dnaEst = dnaGet (est) ;
      if (dnaEst && arrayMax(dnaEst) > 10)
	fixGeneVector (est, 0, 0, 0, 0, 0, dnaEst) ;
      arrayDestroy (dnaEst) ;
    }
  for (jj = 0 ; jj < keySetMax (genes) ; jj++)
    {
      gene = keySet (genes, jj) ;
      if ((Gene = bsCreate (gene)))
	{
	  z1 = z2 = 0 ; dna = 0 ; cosmid = 0 ;
	  if (bsGetKey (Gene, _Genomic_sequence, &cosmid) &&
	      bsGetArray (Gene, _Covers, aa, 4) &&
	      arrayMax(aa) > 3)
	    {
	      uu = arrp (aa, 0, BSunit) ;
	      z1 = uu[2].i ; z2 = uu[3].i ;
	      if (cosmid && (dna = dnaGet(cosmid)))
		{
		  if (z1 > z2)
		    {
		      int i = arrayMax(dna) ;
		      
		      reverseComplement (dna) ;
		      z1 = i - z1 + 1 ; z2 = i - z2 + 1 ;
		    }
		}
	    }
	  if (!dna) { cosmid = 0 ; z1 = z2 = 0 ; }
	  bsGetArray (Gene, _Assembled_from, aa, 5) ;
	
	  /* find last base of alignment of each est */
	  for (ia = (long)arrayMax (aa) - 5 ; ia > 0 ; ia -= 5)
	    {
	      int  ii ;

	      uu = arrp (aa, ia, BSunit) ;
	      est = uu[2].k ; /*   x2 = uu[4].i ; */
	      if (keySetFind (ests, est, &ii))
		{
		  if (!keySet (ks1, ii))
		    {
		      keySet (ks1, ii) = 1 ;
		      dnaEst = dnaGet (est) ;
		      if (dnaEst && arrayMax(dnaEst) > 10)
			fixGeneVector (est, Gene, aa, z1, z2, dna, dnaEst) ;
		      arrayDestroy (dnaEst) ;
		    }
		}
	    }
	  bsDestroy (Gene) ;
	  arrayDestroy (dna) ;
	}
    }
  arrayDestroy (aa) ; 
  keySetDestroy (ks1) ;
  keySetDestroy (genes) ;
  keySetDestroy (ests) ;
  keySetDestroy (estsNoGene) ;
} /* fixVector */

/*********************************************************************/
/*********************************************************************/
/* returns the number of A to unzip, or to zip
 * if a2 has 8 A above and 13 below, the system returns nA = zip==FALSE ? 8 : -13 
 * so a2 - nA is the position wanted
 */
static int abiFixUnzip (Array dna, int a2, BOOL isDown, BOOL zip)
{
  const char *cp ;
  int nA = 0, amax = arrayMax(dna) ;

  if (!isDown)
    a2 = -a2 ;
  
  /* start at a2 and count the A */
  if (dna && a2 < amax && a2 >= 1)   
    {
      if (zip)
	{
	  cp = arrp (dna, a2-1, char) ; /* last base of mRNA */
	  if(isDown)
	    {
	      while (*cp == A_ && a2 >= 1 && a2 < amax)
		{ 
		  nA-- ;
		  a2++ ; cp++ ;
		}
	    }
	  else
	    {
	      while (*cp == T_ && a2 >= 1 && a2 < amax)
		{ 
		  nA-- ;
		  a2-- ; cp-- ; 
		}
	    }
	}
      else
	{
	  cp = arrp (dna, a2-1, char) ; /* last base of mRNA */
	  if(isDown)
	    {
	      while (*cp == A_ && a2 >= 1 && a2 < amax)
		{ 
		  nA++ ;
		  a2-- ; cp-- ;
		}
	    }
	  else
	    {
	      while (*cp == T_ && a2 >= 1 && a2 < amax)
		{ 
		  nA++ ;
		  a2++ ; cp++ ; 
		}
	    }
	}

    }
  return nA ;
}  /* abiFixUnzip */

/*********************************************************************/
/*********************************************************************/

static BOOL abiFixIsGenomeArich (Array dna, BOOL isDown, int a2, int *nAp, int *nAllp, BOOL isKim)
{
  BOOL isRich = FALSE ;
  char *cp, cc, AA = A_, TT = T_ ;
  int i, j, nn, nPure, nPureBest = 0 ; 
  int minAA[] = { 8, 9, 10, 10, 11, 12, 12, 13, 14, 14, 15, 16, 16, 17, 18, 99, 99, 0} ;
  /*            { 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 0}  */
  /* 8/8 and 9/9 will be removed if more than 2 clones 10/11 11/12 if more than 5 clones */
  int minAK[] = { 8, 9, 10, 10, 11, 12, 12, 13, 14, 14, 15, 16, 16, 17, 18, 99, 99, 0} ;
  /*            { 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 0}  */
  /* 8/8 and 9/9 will be removed if more than 2 clones 10/11 11/12 if more than 5 clones */
  int *minA = isKim ? minAK : minAA ;

  if (isKim) { AA = A_ | G_ ; TT = T_ | C_ ; }
  if (!isDown)
    a2 = -a2 ;
  
  if (!dna ||
      a2 + 25 > arrayMax (dna) ||
      a2 < 10
      )
    goto done ;

  /* start around a2 and look the rate of A */
  if (isDown)
    {
      cp = arrp (dna, a2-1, char) ; /* last base of mRNA */
      while (*cp & AA) cp-- ; /* eat back the A */
      cp++ ;  /* first base after mRNA */
      for (i = 1, nn = nPure = 0, isRich = FALSE ; !isRich && i <= 24 ; i++)
	{
	  cc = *cp++ ;
	  if (isKim && cc == G_) cc = A_ ;
	  switch ((int)cc)
	    {
	    case A_:
	      nn++ ; nPure++ ;
	      if (nPure > nPureBest) nPureBest = nPure ;
	      if (nPureBest >= 5 && nn >=8 && i >= 8 && nn >= minA[i-8])
		{ 
		  while (*cp++ & AA) { i++ ; nn++ ;}
		  while (1) 
		    { 
		      j=0 ; while(*cp++ & AA) j++ ;
		      if (j>=4)
			{
			  i += j + 1 ; nn += j ;
			}
		      else
			break ;
		    }
		  isRich = TRUE ; if(nAp) *nAp = nn ; if (nAllp) *nAllp = i ;
		}
	      break ;
	    default:
	      if (nPure > nPureBest) nPureBest = nPure ;
	      nPure = 0 ;
	      break ;
	    }
	  if (i == 6 && nn < 4) break ;
	}
    }
  else /* same on up strand, avoids complementing the whole cosmid */
    {
      cp = arrp (dna, a2 - 1, char) ; /* last base of mRNA */
      while (*cp & TT) cp++ ; /* eat back the A */
      cp-- ;  /* first base after mRNA */
      for (i = 1, nn = nPure = 0, isRich = FALSE ; !isRich && i <= 24 ; i++)
	{
	  cc = *cp-- ;
	  if (isKim && cc == C_) cc = T_ ;
	  switch ((int)cc)
	    {
	    case T_:
	      nn++ ; nPure++ ;
	      if (nPure > nPureBest) nPureBest = nPure ;
	      if (nPureBest >= 5 && nn >=8 && i >= 8 && nn >= minA[i-8])
		{ 
		  while (*cp-- & TT) { i++ ; nn++ ;}
		  while (1) 
		    {
		      j=0; while(*cp-- & TT) j++ ;
		      if (j>=4)
			{ i += j+1 ; nn += j ;}
		      else
			break ;
		    }
		  isRich = TRUE ; if(nAp) *nAp = nn ; if (nAllp) *nAllp = i ;
		}
	      break ;
	    default:
	      if (nPure > nPureBest) nPureBest = nPure ;
	      nPure = 0 ;
	      break ;
	    }
	  if (i == 6 && nn < 4) break ;
	}
    }

 done:
  return isRich ;
} /* abiFixIsGenomeArich */

/*********************************************************************/
typedef struct aaaStruct { int a1, a2, a22, m1, m2, ma1, ma2, x1, x2, a20, e1, e2, de2, nClones3, nClones5, isSl, nAClones, nAClones0, mClones3, mClones03, nA, nAll, tail, group ; KEY mrna, est, est0, clone ; int method ; int signal, aSignal, aRich ; BOOL top, iPriming, polyAPriming  ; unsigned char buffer[7] ; } AAA ;

static int abiFixHasPolyASignal (AAA *ap, Array gDna)
{
  BOOL ok = FALSE ;
  int i, ii, nn = 0, dx = 0 ;
  int MXLN = 32 ;
  /* favor a few motifs recursivelly */
  unsigned char *cp, buf[60], *buffer = ap->buffer ;

  /* liste modifiee le 2010_02_17 
   *
   AATAAA tataaa cataaa gataaa aatgaa tatgaa agtaaa attaaa catgaa aataca aatata aatgta actaaa aataag gatgaa ggtaaa aaaaaa aattaa aataat tgtaaa aagaaa tttaaa agtgaa gaaaaa aacaaa tctaaa aaataa caaaaa aaacaa aatcaa
   */

  unsigned char polyASignal[] = { 

    A_,A_,T_,A_,A_,A_, 0,
    A_,A_,T_,G_,A_,A_, 0,
    G_,A_,T_,A_,A_,A_, 0,
    C_,A_,T_,A_,A_,A_, 0,
    T_,A_,T_,A_,A_,A_, 0,
    T_,A_,T_,G_,A_,A_, 0,
    C_,A_,T_,G_,A_,A_, 0,
    A_,G_,T_,A_,A_,A_, 0,
    A_,A_,T_,A_,C_,A_, 0,
    G_,A_,T_,G_,A_,A_, 0,
    A_,A_,T_,A_,T_,A_, 0,
    A_,T_,T_,A_,A_,A_, 0,
    A_,A_,A_,T_,A_,A_, 0,
    A_,A_,T_,A_,A_,G_, 0,
    A_,C_,T_,A_,A_,A_, 0,
    G_,G_,T_,A_,A_,A_, 0,
    A_,A_,T_,G_,T_,A_, 0,
    A_,A_,T_,A_,A_,T_, 0,
    A_,A_,T_,T_,A_,A_, 0,
    T_,C_,T_,A_,A_,A_, 0,
    T_,G_,T_,A_,A_,A_, 0,
    A_,G_,T_,G_,A_,A_, 0,
    A_,A_,A_,A_,A_,A_, 0,
    A_,A_,G_,A_,A_,A_, 0,
    A_,A_,C_,A_,A_,A_, 0,
    T_,T_,T_,A_,A_,A_, 0,
    G_,A_,A_,A_,A_,A_, 0,
    C_,A_,A_,A_,A_,A_, 0,
    A_,A_,A_,C_,A_,A_, 0,
  
    0
  } ;
  
  ap->signal = 0 ;
  if (ap->a1 < ap->a2)
    {
      dx = 1 ;
      if (ap->a2 < MXLN)
	return 0 ;
      for (i = 0, cp = arrp (gDna, ap->a2 - MXLN, unsigned char) ; i < MXLN ; cp++, i++)
	buf [i] = *cp ;
    }
  else
    {
      dx = 3 ;
      if (-ap->a2 + MXLN > arrayMax (gDna))
	return 0 ;
      for (i = 0, cp = arrp (gDna, -ap->a2 + MXLN, unsigned char) ; i < MXLN ; cp--, i++)
	buf [i] = complementBase(*cp) ;
    }
  buf[MXLN] = 0 ;

  for (ii = 0 ; !ok && polyASignal[7*ii] > 0 ; ii++)
    {
      nn = dnaPickMatch (buf, MXLN, polyASignal + 7*ii, 0, 0) ;
      if (nn)
	{
	  if (ii == 0)
	    {
	      ok = TRUE ; ap->signal = 1 ; ap->aSignal = MXLN - nn + dx ;
	      strcpy ((char *)buffer, "AATAAA") ;
	    }
	  else
	    {
	      ok = TRUE ; ap->signal = 2 ; ap->aSignal = MXLN - nn + dx ;
	      for (i = 0 ; i < 6 ; i++)
		buffer[i] = dnaDecodeChar[(int)buf[nn+i-1]] ;
	      buffer[i] = 0 ;
	    }
	}
    }
  if (0 && ! ok) /* any one letter variant */
    {
      nn = dnaPickMatch (buf, MXLN, polyASignal, 1, 0) ;
      if (nn)
	{
	  ok = TRUE ; ap->signal = 2 ; ap->aSignal = MXLN - nn + dx ;
	  for (i = 0 ; i < 6 ; i++)
	   buffer[i] = dnaDecodeChar[(int)buf[nn+i-1]] ;
	  buffer[i] = 0 ;
	}
    }
  if (0 && ! ok)
    {
      nn = dnaPickMatch (buf, MXLN, polyASignal, 2, 0) ;
      if (nn)
	{
	  ap->signal = 2 ; ap->aSignal = MXLN - nn + dx ;
	  for (i = 0 ; i < 6 ; i++)
	   buffer[i] = dnaDecodeChar[(int)buf[nn+i-1]] ;
	  buffer[i] = 0 ;
	}
    }
  
  return nn > 0 ? MXLN + 1 - nn : 0 ; /* distance from end OR zero */
} /* abiFixhasPolyASignal */

/*********************************************************************/

static void showAAA (Array a, BOOL all)
{
  int ii ;
  AAA *ap ;

  if (arrayExists (a))
    {
      printf ("\n") ;
      for (ii = 0 ; ii < arrayMax (a) ; ii++)
	{
	  ap = arrp (a, ii, AAA) ;
	  if(ap->est == 0 && ! all) continue ;
	  printf ("%3d: %03d, a=%d %d   x=%d %d  %se=%d %d nA0=%d nA=%d nc5=%d nc3=%d %s%d %s%s %s %s %s %s tag:%d:%d\n"
		  , ii
		  , ap->group
		  , ap->a1, ap->a2
		  , ap->x1, ap->x2
		  , ap->est ? "*" : " "
		  , ap->e1, ap->e2
		  , ap->nAClones0
		  , ap->nAClones
		  , ap->nClones5
		  , ap->nClones3
		  , ap->aRich  ? "R" : ""
		  , ap->aRich 
		  , ap->signal ? (ap->signal == 1 ? "S" : "s" ) : ""
		  , ap->top ? "Top" : ""
		  , ap->tail ? "Tail" : ""
		  , ap->mrna ? name(ap->mrna) : "0"
		  , ap->est > 1 ? name(ap->est) : "0"
		  , ap->est0 ? name(ap->est0) : "0"
		  , ap->method 
		  , ap->mClones3
		  ) ;
	}
    }
  if (0) showAAA (a,0) ; /* for compiler happiness */
}

/*********************************************************************/

static int aaaMrnaOrder (const void *va, const void *vb)
{
  const AAA *a = (const AAA*)va ;
  const AAA *b = (const AAA*)vb ;

  if (a->mrna < b->mrna) return -1 ;
  else if (a->mrna > b->mrna) return 1 ;
  
  if (a->x2 !=  b->x2) return a->x2 - b->x2 ;
  if (a->de2 != b->de2) return a->de2 - b->de2 ;
  return 0 ;
}

/*********************************************************************/

static int aaaTgOrder (const void *va, const void *vb)
{
  const AAA *a = (const AAA*)va ;
  const AAA *b = (const AAA*)vb ;

  if (a->group !=  b->group) return a->group - b->group ;
  if (a->a2 !=  b->a2) return a->a2 - b->a2 ;
  if (a->de2 != b->de2) return a->de2 - b->de2 ;
  return a->x2 - b->x2 ;
}

/*********************************************************************/

static int aaaTgOrder20 (const void *va, const void *vb)
{
  const AAA *a = (const AAA*)va ;
  const AAA *b = (const AAA*)vb ;

  if (a->group !=  b->group) return a->group - b->group ;
  if (a->a20 !=  b->a20) return a->a20 - b->a20 ;
  if (a->de2 != b->de2) return a->de2 - b->de2 ;
  return a->x2 - b->x2 ;
}

/*********************************************************************/

static BOOL abiFixLabelGatherPolyA (KEY mrna, Array gDna, Array aa5,  Array aa3, KEYSET clones5, KEYSET clones3, DICT *dict)
{
  int ii, jj, a01, a1, a2, ma1=1, m1, m2, p2, nSl = 0 ;
  Array mDna = 0 ;
  AAA *ap, *ap1 ;
  KEY est, product ;
  Array units = 0 ;
  BSunit *uu ;
  BOOL hasProduct = FALSE, isDown = TRUE ;
  KEY _Valid5p = str2tag ("Valid5p") ;
  KEY _Ref_mRNA = str2tag ("Ref_mRNA") ;
  KEY _Ref_seq = str2tag ("Ref_seq") ;
  KEY _Is_partial = str2tag ("Is_partial") ;
  KEY _Complete_CDS = str2tag ("Complete_CDS") ;
  KEY _Internal_priming = str2tag ("Internal_priming") ;
  KEY _PolyA_after_base = str2tag ("PolyA_after_base") ;
  OBJ Mrna = bsUpdate (mrna) ;
  BOOL CHEAT = FALSE   ;  /* if TRUE: special trick to export for kris the global clustering of her polya flags */

  if (CHEAT)
    {
      mDna = dnaCopy (gDna) ;
    }
  else
    mDna = dnaGet (mrna) ;
  p2 = 1 ; /* by default end of the mRNA */
  units = arrayCreate (3000, BSunit) ;
  
  /* locate the mRNA on the cosmid */
  a01 = a1 = a2 = m1 = m2 = 0 ;
  if (bsGetArray (Mrna, _DNA, units, 2))
    {
      ii = 0 ;
      uu = arrp (units, ii, BSunit) ; 
    }
   if (bsGetArray (Mrna, _Covers, units, 5))
    {
      ii = arrayMax (units) - 5 ;
      uu = arrp (units, ii, BSunit) ; 
      a01 = a1 = uu[2].i ; /* on cosmid */
      a2 = uu[4].i ; 
      if (a1 > a2) isDown = FALSE ;
    }

  /* locate the last exon of the mRNA */
  if (bsGetArray (Mrna, _Splicing, units, 5))
    {
      ii = arrayMax (units) - 5 ;
      uu = arrp (units, ii, BSunit) ; 
      if (isDown)
	{ 
	  a2 = a1 + uu[1].i - 1 ; 
	  a1 = a1 + uu[0].i - 1 ; /* on cosmid */
	}
      else
	{
	  a2 = a1 - uu[1].i + 1 ; 
	  a1 = a1 - uu[0].i + 1 ; 
	}
      ma1 = uu[0].i ; /* on unspliced mrna */
      /* ma2 = uu[1].i ;  */
      m1 = uu[2].i ; /* on mrna */
      m2 = uu[3].i ; 
    }
  if (CHEAT) m2 = a2 = arrayMax(mDna) ;
  /* coords of reliable product on mRNA */
  if (bsGetArray (Mrna, _Product, units, 3))
    for (ii = 0 ; ii < arrayMax (units) ; ii+= 3)
      {
	uu = arrp (units, ii, BSunit) ;
	product = uu[0].k ;
	if (keyFindTag (product, _Best_product))
	  {
	    if (keyFindTag (product, _Good_product))
	      {
		if (keyFindTag (product, _COOH_complete))
		  { p2 = uu[2].i ; hasProduct = TRUE ; }
		if (keyFindTag (product, _NH2_complete))
		  { /*has5pProduct = TRUE ;*/ }
	      }
	    break ;
	  }
      }    
  if (hasProduct)
    {
      p2 += 3 ; /* eat the stop */
      if (p2 > m1 && p2 < m2)
	{
	  int da = p2 - m1 ;
	  m1 += da ;
	  ma1 += da ;
	  if (isDown) a1 += da ;
	  else a1 -= da ;
	}
    }
  if (CHEAT) hasProduct = FALSE ;
  /* gather the contributing 5p ESTs */
  jj = arrayMax (aa5) ;
  arrayMax (units) = 0 ;
  if (bsGetArray (Mrna, _Constructed_from, units, 5))
    for (ii = 0 ; ii < arrayMax (units) ; ii+= 5)
      {
	int e1, e2, x1, x2, v1 = 1, x, isSl = 0, n ;
	BOOL top = FALSE ;
	OBJ Est ;
	KEY clone, sl ;

	uu = arrp (units, ii, BSunit) ;
	x1 = uu[0].i ;  
	x2 = uu[1].i ; 
	e1 = uu[3].i ;
	e2 = uu[4].i ; /* est coords */
	est = uu[2].k ; 
	if (keyFindTag (est, _Is_AM)) 
	  continue ;
	Est = bsCreate (est) ;
	bsGetKey (Est, _cDNA_clone, &clone) ;
	if ((bsFindTag (Est, _Forward) &&
	     (
	      e1 < 5 || 
	      (bsGetData (Est, _Vector_clipping, _Int, &v1) && e1 < v1 + 5)
	      )
	     )
	    )
	  top = TRUE ;
	if (bsFindTag (Est, _Forward) &&
	    bsFindTag (Est, _Composite) &&
	    bsGetData (Est, _Composite, _Int, &n) &&
	    bsFindTag (Est, _Real_5prime)
	    )
	  { top = TRUE ; isSl = 1 ; nSl+= n; }
	if (bsFindTag (Est, _Forward) &&
	    bsGetKey (Est, _Transpliced_to, &sl) &&
	    bsGetData (Est, _bsRight, _Int, &x) &&
	    e1 < x + 5)
	  { top = TRUE ; isSl = 1 ; nSl++ ; }
	if (bsFindTag (Est, _Reverse) &&
	    bsGetKey (Est, _Transpliced_to, &sl) &&
	    bsGetData (Est, _bsRight, _Int, &x) &&
	    e1 < x + 5 && e1 > x - 5)
	  { top = TRUE ; isSl = 1 ; nSl++ ; }
	bsDestroy (Est) ;
	
	if (!top)
	  continue ;
	
	ap = arrayp (aa5, jj++, AAA) ;
	ap->mrna = mrna ;
	ap->est = est ;
	ap->est0 = est ;
	ap->clone = clone ;
	if (isDown) ap->a1 = a1 + x1 - m1 ;
	else  ap->a1 = -(a1 -x1 + m1) ;   /* negate coords, to simplify sorting */
	if (isDown) ap->ma1 = a1 ; /* last base of last exon on cosmid */
	else  ap->ma1 = -a1 ;   /* negate coords, to simplify sorting */
	ap->m1 = ma1 + x1 - m1 ; /* spliced mrna coord */
	ap->x1 = x1 ; /* spliced mrna coord */
	ap->x2 = x2 ;
	ap->e1 = e1 ;
	ap->e2 = e2 ;
	ap->top = top ;
	ap->isSl = isSl ;
	if (keySetInsert (clones5, clone)) /* only if new, since may be in a variant */
	  ap->nClones5 = 1 ;
      }
  /* gather the contributing 3p ESTs */
  jj = arrayMax (aa3) ;
  if (1)
    for (ii = 0 ; ii < arrayMax (units) ; ii+= 5)
      {
	int e1, e2, de2 = 0, x2, pA = 0, v1 = 1, tail = 0, n = 1 ;
	BOOL isMrna, iPriming, polyAPriming, isComposite = FALSE ;
	OBJ Est ;
	KEY clone ;

	uu = arrp (units, ii, BSunit) ;
	/* 	x1 = uu[0].i ;  */
	x2 = uu[1].i ; 
	e1 = uu[3].i ;
	e2 = uu[4].i ; /* est coords */
	if (x2 < m1) /* we are only interested in the terminal exon, downstream of the best/good CDS */
	  continue ;
	est = uu[2].k ; 
	isComposite = keyFindTag (est, _Composite) ;
	if (keyFindTag (est, _Is_AM)) 
	  continue ;
	Est = bsCreate (est) ;
	bsGetKey (Est, _cDNA_clone, &clone) ;
	bsGetData (Est, _Composite, _Int, &n) ;

	if (!strncmp(name(clone),"yk",2) && !strncmp(name(est),"GB",2))
	  { bsDestroy (Est); continue ; }
	/*
	*/
	iPriming = keyFindTag (clone, _Internal_priming) ;
	polyAPriming = keyFindTag (clone, _Primed_on_polyA) ;
	if (bsFindTag (Est, _PolyA_after_base))
	  {
	    if (bsFindTag (Est, _Reverse) &&
		bsFindTag (Est, _Composite)
		)
	      bsGetData (Est, _Composite, _Int, &n) ;
	    if (bsGetData (Est, _PolyA_after_base, _Int, &pA) &&
		(
		 (e1 < e2 && e2 > pA - 5 && e2 < pA + 5) || /* fully aligned */
		 (e1 > e2 && e2 > pA - 5 && e2 < pA + 5)
		 )
		)
	      {
		tail = 2 ;
		if (e1 < e2)
		  {
		    de2 = pA - e2 ;
		    if (de2 > 0)
		      { x2 += de2 ; e2 += de2 ; a2 += de2 ; de2 = 0; } /* on deroule */
		  }
		else
		  {
		    de2 = e2 - pA ;
		    if (de2 > 0)
		      { x2 += de2 ; e2 -= de2 ; a2 += de2 ; de2 = 0; } /* on deroule */
		  }
	      }
	  }

	if (! isComposite)
	  {
	    isMrna = bsFindTag (Est, _Ref_mRNA) &&
	      e2 > m1 + 10 &&
	      !bsFindTag (Est, _Ref_seq) &&
	      !bsFindTag (Est, _Is_partial) ;
	    if (0 && isMrna && bsFindTag (Est, _Complete_CDS))
	      tail = 2 ;
	    de2 = 0 ;
	    if (e1 > e2)
	      {
		v1 = 1 ;
		if (pA > 0) v1 = pA ;
		else
		  bsGetData (Est, _Vector_clipping, _Int, &v1) ;
		if (e2 > v1)  
		  de2 = e2 - v1 ; /* number of unaligned bases */
	      }
	    else if (isMrna) 
	      {
		int dummy, v2 ;
		v2 = 9999999 ;
		if (pA > 0) v2 = pA ;
		else
		  if (bsGetData (Est, _Vector_clipping, _Int, &dummy))
		    bsGetData (Est, _bsRight, _Int, &v2) ;
		if (e2 < v2)  
		  de2 = v2 - e2 ; /* number of unaligned bases */
	      }
	  }
	bsDestroy (Est) ;

	if (iPriming)
	  continue ;
	if (iPriming || de2 > (polyAPriming ? 12 : 6))
	  continue ;
	if (iPriming || de2 > (polyAPriming ? 12 : 6) || (e1 < e2))
	  continue ;
	if (iPriming || de2 > (polyAPriming ? 12 : 6) || (e1 < e2 && !tail && !isMrna))
	  continue ;
	if (iPriming || de2 > (polyAPriming ? 12 : 6) || (e1 < e2 && !tail && !isMrna) 
	    || (e1 > e2 && !tail && !polyAPriming))
	  continue ;
	if (iPriming || de2 > (polyAPriming ? 12 : 6) || (e1 < e2 && !tail && !isMrna) 
	    || (e1 > e2 && !tail && !polyAPriming))
	  continue ;
	if (0 &&               /* 0: on deroule en cooperatif */
	    e1 > e2 && polyAPriming) /* 1: on deroule a priori */
	  { x2 += de2 ; a2 += de2 ; de2 = 0; }
	ap = arrayp (aa3, jj++, AAA) ;
	ap->mrna = mrna ;
	ap->est = est ;
	ap->est0 = est ;
	ap->clone = clone ;
	if (isDown) ap->a2 = a1 + x2 - m1 ;
	else  ap->a2 = -(a1 -x2 + m1) ;   /* negate coords, to simplify sorting */
	ap->a20 = ap->a2 ;
	if (isDown) ap->ma2 = a2 ; /* last base of last exon on cosmid */
	else  ap->ma2 = -a2 ;   /* negate coords, to simplify sorting */
	ap->m2 = ma1 + x2 - m1 ; /* spliced mrna coord */
	ap->x2 = x2 ; /* current spliced mrna coord */
	ap->e1 = e1 ;
	ap->e2 = e2 ;
	ap->de2 = de2 ;
	ap->tail = tail ;
	ap->nAClones = tail ? 1 : 0 ;
	ap->nAClones0 = ap->nAClones ;
	ap->iPriming = iPriming ;
	ap->polyAPriming = polyAPriming ;
	if (keySetInsert (clones3, clone)) /* only if new, since may be in a variant */
	  {
	    ap->nClones3 = n ;
	    ap->mClones03 = n ;
	  }
      }

  /* gather valid3p from the Feature AAA stored on the cosmid */
  {
    KEYSET cosmids = queryKey (mrna, ">Genomic_sequence ; {Feature} SETOR {>parts} SETOR {>In_junction} SETOR {>In_junction ; > Parts } SETOR {>source} ; Feature") ;
    int iCosmid ;
    int g1, g01 = 0 ;
    KEY cosmid = keyGetKey (mrna, _Genomic_sequence) ;
    OBJ Cosmid = bsCreate (cosmid) ;    

    if (Cosmid && bsGetArray (Cosmid, _IntMap, units, 3))
      {
 	uu = arrp (units, 0, BSunit) ;
	g01 = uu[1].i ; 
      }
    bsDestroy (Cosmid) ;

    for (iCosmid = 0 ; g01 && cosmids && iCosmid < keySetMax (cosmids) ; iCosmid++)
      {
	cosmid = keySet (cosmids, iCosmid) ;
	Cosmid = bsCreate (cosmid) ;    
	
	if (bsGetArray (Cosmid, _IntMap, units, 3))
	  {
	    uu = arrp (units, 0, BSunit) ;
	    g1 = uu[1].i ; 
	  }
	if (Cosmid && bsGetArray (Cosmid, str2tag("Feature"), units, 5))
	  for (ii = 0 ; ii < arrayMax (units) ; ii+= 5)
	    {
	      int w1, w2, w3 ;
	      KEY method ;
	      const char *ccq ;
	    
	      uu = arrp (units, ii, BSunit) ;
	      method = uu[0].k ;   /* should be AAA        */
	      if (strncmp (name(method), "AAA_", 4))
		continue ;
	      if (!strncmp (name(method), "AAA_cDNA", 8))
		continue ;
	      w1 = uu[1].i + g1 - g01 ;       /* coord of the AAAA on the original cosmid */
	      w2 = uu[2].i + g1 - g01 ;       /* coord of the AAAA on the original cosmid */
	      w3 = uu[3].f ;       /* number of supports              */
	      ccq = uu[4].s ;
	      
	      if (!w3 || (isDown && w1 > w2) || (!isDown && w1 < w2) || !ccq) continue ;	      
	      if (! (isDown && w1 > a1 && w1 < a2 + 300) &&
		  ! (!isDown && w1 < a1 && w1 > a2 - 300)
		  )
		continue ;
	      
	      if (!isDown) w1 = -w1 ; /* negate coords, to simplify sorting */

	      
	      ap = arrayp (aa3, jj++, AAA) ;
	      ap->mrna = mrna ;
	      ap->est = 1 ;
	      ap->est0 = 0 ;
	      ap->mClones3 = w3 ;
	      ap->mClones03 = w3 ;
	      dictAdd (dict, ccq, &(ap->method)) ;
	      if (1)
		{
		  if ( ( 1 || strstr (name(method), "AAA_454_Kim")) &&
		       ! ap->aRich && 
		       abiFixIsGenomeArich (gDna, isDown, w1, 0, 0, TRUE))
		    ap->aRich = 1 ;
		}
	      
	      ap->a2 = w1 ;
	      ap->a20 = ap->a2 ;
	      
	      ap->ma2 = w1 ; /* on unspliced mrna */
	      if (isDown) ap->m2 = w1 - a01 + 1 ; /* spliced mrna coord */
	      else ap->m2 = a01 + w1 + 1 ;
	      
	      if (isDown) ap->x2 = w1 - a1 + m1 ; /* spliced mrna coord */
	      else ap->m2 = a1 + w1 + m1 ;
	    }
	bsDestroy (Cosmid) ;
      }
    keySetDestroy (cosmids) ;
  }

  arraySort (aa3, aaaMrnaOrder) ;
  if (0) showAAA (aa3,0) ;

  {
    int x = 1 ;
    
    if (bsFindTag (Mrna, _Valid5p))
      bsRemove (Mrna) ;
    if (bsFindTag (Mrna, _Transpliced_to) || bsFindTag (Mrna, str2tag("Aggregated_5p_clones")))
      {
	x = 1 ;
	bsAddData (Mrna, _Valid5p, _Int, &x) ;
	bsAddData (Mrna, _bsRight, _Int, &x) ;
	x = nSl ;
	bsAddData (Mrna, _bsRight, _Int, &x) ;
      }
  } /* if products) {} */
  
  
  
  /* unzip to first nonA base */
  if (arrayMax(aa3))
    {
      BOOL sortNeeded = FALSE ;
      BOOL doZip = FALSE ; /* 2024_09_24 was TRUE */ ;

      for (ii = arrayMax (aa3) - 1, ap = arrp (aa3,ii, AAA) ; ii >= 0 ; ii--, ap--)
	{
	  int dx = abiFixUnzip (gDna, ap->a2, isDown, doZip) ; /*returns number of A to unzip */
	  if (dx)
	    {
	      sortNeeded = TRUE ;
	      if (isDown)
		{
		  ap->x2 -= dx ; 
		  ap->m2 -= dx ;
		  ap->ma2 -= dx ;
		  ap->a2 -= dx ;
		  ap->a20 -= dx ;
		}
	      else
		{
		  ap->x2 -= dx ; 
		  ap->m2 -= dx ;
		  ap->ma2 -= dx ; /* since i negated the coords it is again minus */
		  ap->a2 -= dx ;
		  ap->a20 -= dx ;
		}
	    }
	}
      if (sortNeeded)
	{
	  arraySort (aa3, aaaMrnaOrder) ;   
	  for (ii = arrayMax (aa3) - 1, ap = arrp (aa3,ii, AAA) ; ii >= 0 ; ii--, ap--)
	    {
	      if (ap->est0)
		continue ;
	      for (jj = ii - 1, ap1 = ap - 1; jj >= 0 ; jj--, ap1--)
		{
		  if (ap->a2 != ap1->a2)
		    break ;
		  if (ap1->est0 || ap->method != ap1->method)
		    continue ;
		  ap->mClones3 += ap1->mClones3 ;
		  ap->mClones03 += ap1->mClones03 ;
		  ap1->est = 0 ;
		  ap1->method = 0 ;
		  ap1->mrna = 999999999 ;
		  ap1->mClones3 = 0 ;
		  ap1->mClones03 = 0 ;	   
		}
	    }
	}
    }

  
  arrayDestroy (units) ;
  if (mrna != 1)
    arrayDestroy (mDna) ;
  bsSave (Mrna) ;
  
  return isDown ;
} /* abiFixLabelGatherPolyA */

/*********************************************************************/

static BOOL getPleaseIsWorm (void)
{
  static int isWorm = -1 ;

  if (isWorm == -1)
    {
      KEYSET ks = query (0, "Find clone Strategy && species == worm") ;
      if (keySetMax (ks) > 0)
	isWorm = 1 ;
      else
	isWorm = 0 ;
      keySetDestroy (ks) ;
    }

  return isWorm ;
} /* getPleaseIsWorm */

/*********************************************************************/

static void abiFixLabelClusterPolyAPair (Array aa3, AAA *ap, AAA *ap1)
{
  AAA *ap2 ;
  int ii2, group = ap1->group ;
	
  /* regroup all the member of the lameduck ap1 group */
  for (ii2 = 0, ap2 = arrp (aa3, 0, AAA) ; ii2 < arrayMax (aa3) ; ii2++, ap2++)
    if (ap2->group == group)
      {
	ap2->group = ap->group ;
	if(ap->est)
	  ap2->est = 0 ; 
	ap2->a22 = ap->a22 ;
      }
  
  /* absorb  the ap1 counts in ap */
  if (ap->tail < ap1->tail)
    ap->tail = ap1->tail ; 
  ap->polyAPriming |= ap1->polyAPriming ;
  ap->nAClones += ap1->nAClones ;
  ap->nClones3 += ap1->nClones3 ;
  ap->mClones3 += ap1->mClones3 ;
  if (ap->a20 == ap1->a20)
    {
      ap->nAClones0 += ap1->nAClones0 ;
      ap->mClones03 += ap1->mClones03 ;
      ap1->nAClones0 = 0 ; ap1->mClones03 = 0 ;
    }
  ap1->est = 0 ; ap1->a22 = ap->a22 ;
} /* abiFixLabelClusterPolyAPair */

/*********************************************************************/

static void abiFixLabelClusterPolyA (Array aa5, Array aa3, Array gDna, BOOL isDown, int *nClonep)
{
  int ii, jj, a22, nRound, cluster ;
  AAA *ap, *ap1 ;
  BOOL isWorm = getPleaseIsWorm () ;
  BOOL sortNeeded ;
  arraySort (aa5, aaaTgOrder) ;
  arraySort (aa3, aaaTgOrder) ;
  /* check for genome A rich */
  if (!arrayMax (aa3))
    return ;
  
  /* extend to next cluster the unaligned sequences */
  if (arrayMax (aa3))
    {
      int de2, a2a, a2n, da ;
      
      for (ii = 0, ap = arrp (aa3, 0, AAA) ; ii < arrayMax (aa3) ; ii++, ap++)
	{
	  de2 = ap->de2 ;
	  if (de2 <= 0)
	    continue ;
	  a2n = ap->a2 ;
	  a2a = ap->a2 + (4*de2)/3 + 10 ;
	  for (jj = ii+1, ap1 = ap + 1 ; jj < arrayMax (aa3) && ap1->a2 < a2a  ; jj++, ap1++)
	    if (ap1->de2 < de2)
	      a2n = ap1->a2 ;
	  if (a2n > ap->a2 + ap->de2)
	    a2n = ap->a2 + ap->de2 ;
	  da = a2n - ap->a2 ;
	  ap->a2 += da ;
	  ap->m2 += da ; ap->x2 += da ; 
	  ap->de2 -= da ;
	}
    }
  arraySort (aa3, aaaTgOrder) ; /* sort again since we shifted */
  
  if (arrayMax (aa3))
    {
      int nA = 0, nAll = 0, b2 =  - arrayMax (gDna) - 1000 ;
      BOOL isRich = FALSE ;
      
      for (ii = 0, ap = arrp (aa3, 0, AAA) ; ii < arrayMax (aa3) ; ii++, ap++)
	{
	  if (ap->a2 != b2 )
	    {
	      if (abiFixIsGenomeArich (gDna, isDown, ap->a2, &nA, &nAll, FALSE))
		isRich = TRUE ;
	      else
		{
		  isRich = FALSE ;
		  if (isWorm &&
		      (strstr (name(ap->clone), "U454") || strstr (name(ap->clone), "Umarco") || (strstr (name(ap->clone), "PA") == name(ap->clone)))
		      )
		    isRich |=  abiFixIsGenomeArich (gDna, isDown, ap->a2, &nA, &nAll, TRUE) ; 
		}
	      b2 = ap->a2 ;
	    }
	  if (isRich)
	    { 
	      ap->aRich |= 1 ;
	      ap->nA = nA ; ap->nAll = nAll ; 
	    }
	}
    }
  
  /* clusterize the aRich and be demanding on number of A for large clusters */
  if (arrayMax (aa3))
    {
      int a2, nRich, nA, nAll, jj0 ;
      
      for (ii = jj0 = arrayMax (aa3) - 1, ap = arrp (aa3,ii, AAA) ; ii >= 0 ; ii--, ap--)
	if (ap->aRich == 1)
	  {
	    nA = nAll = nRich = 0 ;
	    a2 = ap->a2 ;
	    nA = ap->nA ; nAll = ap->nAll ; 
	    for (jj = jj0, ap1 = arrp (aa3,jj0, AAA) ; jj >= 0 ; jj--, ap1--)
	      {
		if (ap1->a2 > a2 + 10)
		  { if (jj0 > 0) jj0-- ; continue ; }
		if (ap1->a2 > a2 - 10)
		  { nRich++ ; if (ap1->nA > nA) { nA = ap1->nA ; nAll = ap1->nAll ; }}
		else
		  break ;
	      }
	    if (
		(nRich > 2 && nA < 10) || 
		(nRich > 5 && nA < 12)
		)
	      ap->aRich = 3 ;
	    if (!ap->nA && ap->de2 < 5) { ap->nA = nA ; ap->nAll = nAll ; }
	  }
      for (ii = jj0 = arrayMax (aa3) - 1, ap = arrp (aa3,ii, AAA) ; ii >= 0 ; ii--, ap--)
	if (ap->aRich == 3)
	  ap->aRich = 0 ;
    }

  /* inherit the A-rich tag around by -20 +10 bp for upgoing reads */
  if (arrayMax (aa3))
    {
      int nA = 0, nAll = 0, a2 = arrayMax (gDna) + 1000 ;
      for (ii = arrayMax (aa3) - 1, ap = arrp (aa3,ii, AAA) ; ii >= 0 ; ii--, ap--)
	if (ap->aRich == 1)
	  { a2 = ap->a2 ; nA = ap->nA ; nAll = ap->nAll ; }
	else if (ap->tail && ap->a2 > a2 - 10)
	  { ap->aRich = 2 ; ap->nA = nA ; ap->nAll = nAll ; }
	else if (!ap->tail && ap->a2 > a2 - 20)
	  { ap->aRich = 2 ; ap->nA = nA ; ap->nAll = nAll ; }
      a2 = - arrayMax (gDna) - 1000 ;
      for (ii = 0, ap = arrp (aa3, 0, AAA) ; ii < arrayMax (aa3) ; ii++, ap++)
	if (ap->aRich == 1)
	  { a2 = ap->a2 ; nA = ap->nA ; nAll = ap->nAll ; }
	else if (ap->tail && ap->a2 < a2 + 10)
	  { ap->aRich = 2 ; ap->nA = nA ; ap->nAll = nAll ; }
	else if (!ap->tail && ap->a2 < a2 + 20)
	  { ap->aRich = 2 ; ap->nA = nA ; ap->nAll = nAll ; }
    }
  
  
  /* export the a-rich flag */
  /* then ignore them unless we are in the worm */
  if (arrayMax (aa3))
    {
      for (ii = 0, ap = arrp (aa3, 0, AAA) ; ii < arrayMax (aa3) ; ii++, ap++)
	{
	  if ((ap->tail || ap->e1 > ap->e2) && ap->nA && 
	      ((isWorm && ap->aRich == 1) || (!isWorm && ap->aRich > 0)))
	    {
	      OBJ Clone = 0 ;
	      char buf[256];
	      
	      sprintf (buf, "%d/%d", ap->nA, ap->nAll) ;
	      if ((Clone = bsUpdate (ap->clone)))
		{
		  bsAddTag (Clone, _Internal_priming_on_A_rich) ;
		  bsAddData (Clone, _bsRight, _Text, buf) ;
		  bsSave (Clone) ;
		  if (nClonep) (*nClonep)++ ;
		}			
	    }
	  if (ap->aRich)
	    {
	      if (isWorm) ap->aRich = - ap->aRich ;
	      else ap->tail = 0 ;
	    }
	  else 
	    {
	      if (keyFindTag (ap->clone, _Internal_priming_on_A_rich))
		{
		  OBJ Clone = 0 ;
		  if ((Clone = bsUpdate (ap->clone)))
		    {
		      if (bsFindTag (Clone, _Internal_priming_on_A_rich))
			bsPrune (Clone) ;
		      bsSave (Clone) ;
		    }
		}
	    }
	}
    }

  /* clusterize */
  for (nRound = 0 ; nRound < 3 ; nRound++)
    {
      arraySort (aa3, aaaTgOrder) ; /* sort again since we shifted */
      if (nRound == 0) /* initialize the groups by position of ap->a2 */
	for (ii = 0, ap = arrp (aa3,ii, AAA) ; ii < arrayMax (aa3) ; ap++, ii++)
	  { ap->group = ii ; ap->a22 = ap->a2 ; }
      cluster = FALSE ;
      sortNeeded = FALSE ;

      if (arrayMax (aa3) > 1)   /* cluster polyA at same position */
	{
	  sortNeeded = FALSE ;
	  for (ii = arrayMax (aa3) - 1, ap = arrp (aa3,ii, AAA) ; ii >= 0 ; ii--, ap--)
	    if (ap->est && ap->aRich <= 0)
	      {
		a22 = ap->a22 ; 
		for (jj = ii - 1, ap1 = ap - 1 ; jj >= 0 ; jj--, ap1--) 
		  {
		    if (ap1->est && ap1->aRich <= 0 && ap1->a22 == a22)
		      {
			cluster = TRUE ;
			sortNeeded = TRUE ;
			abiFixLabelClusterPolyAPair (aa3, ap, ap1) ;
		      }
		    if (ap1->a22 > a22)
		      break ;
		  }
	      }
	  if (sortNeeded)
	    arraySort (aa3, aaaTgOrder) ; /* sort again since we shifted */
	}
      
      /*
       * Type 1 == polyA seen,   ap->nAclones > 0
       * Type 2 == deep seq,     ap->mClones3 > 0
       * Type 3 == yuji 3p clones, ap->nClones3 > 0 && ap->nAclones == 0
       */

      if (arrayMax (aa3) > 1)  /* cluster on type 1 (polyA seen) positions the clusters witout polyA if they are within 25 bp */
	{
	  sortNeeded = FALSE ;
	  for (ii = arrayMax (aa3) - 1, ap = arrp (aa3,ii, AAA) ; ii >= 0 ; ii--, ap--)
	    if (ap->est && (ap->nAClones ||  (nRound > 1 && ap->mClones3)) && ap->aRich <= 0)
	      {
		a22 = ap->a22 ; 
		for (jj = ii - 1 ; jj >= 0 ; jj--)
		  {
		    ap1 = arrp (aa3,jj, AAA) ;
		    if (! ap1->est) continue ;
		    if (ap1->nAClones || ap1->a22 < a22 - 25) break ;
		    if (ap1->aRich <= 0)
		      {
			if (ap1->a22 <= a22 && ap1->a22 >= a22 - (ap1->mClones3 ? 10 : 25))
			  {
			    cluster = TRUE ;
			    sortNeeded = TRUE ;
			    abiFixLabelClusterPolyAPair (aa3, ap, ap1) ;
			  }
		      }
		  }
		for (jj = ii + 1 ; jj < arrayMax (aa3) ; jj++)
		  {
		    ap1 = arrp (aa3,jj, AAA) ;
		    if (! ap1->est) continue ;
		    if (ap1->nAClones || ap1->a22 > a22 + 12) break ;
		    if (ap1->aRich <= 0)
		      {
			if (ap1->a22 >= a22 && ap1->a22 <= a22 + (ap1->mClones3 ? 5 : 12))
			  {
			    cluster = TRUE ;
			    sortNeeded = TRUE ;
			    abiFixLabelClusterPolyAPair (aa3, ap, ap1) ;
			  }
		      }
		  }
	      }
	  if (sortNeeded)
	    arraySort (aa3, aaaTgOrder) ; /* sort again since we shifted */
	}
      
      if (arrayMax (aa3) > 1)  /* cluster several type 1 or 2 if they are within 12 bp on the most 3p represented position */ 
	{
	  sortNeeded = FALSE ;
	  for (ii = arrayMax (aa3) - 1, ap = arrp (aa3,ii, AAA) ; ii >= 0 ; ii--, ap--)
	    if (ap->est && ap->nAClones + ap->mClones3 > 0 && ap->aRich <= 0)
	      {
		a22 = ap->a22 ;
		for (jj = ii-1, ap1 = ap - 1 ; jj >= 0 ; jj--, ap1--) 
		  if (ap1->a22 > a22 - 11)
		    {
		      if (ap1->est && ap1->aRich <= 0)
			{
			  sortNeeded = TRUE ;
			  cluster = TRUE ;
			  abiFixLabelClusterPolyAPair (aa3, ap, ap1) ;
			}
		    }
		  else
		    break ;
	      }
	  if (sortNeeded)
	    arraySort (aa3, aaaTgOrder) ; /* sort again since we shifted */
	}

      if (arrayMax (aa3) > 1 && nRound > 1) /* cluster type 3 on second iteration if not well supported */ 
	{
	  sortNeeded = FALSE ;
	  for (ii = arrayMax (aa3) - 1, ap = arrp (aa3,ii, AAA) ; ii >= 1 ; ii--, ap--)
	    if (ap->est && ap->aRich <= 0)
	      {
		a22 = ap->a22 ;
		for (jj = ii-1, ap1 = ap - 1 ; jj >= 0 ; jj--, ap1--)
		  if (ap1->est && ap1->aRich <= 0 && 
		      ! ap1->nAClones && ! ap1->mClones3 && ap1->nClones3 < 6 
		      )
		    {
		      if (ap1->a22 > a22 - 25)
			{
			  sortNeeded = TRUE ;
			  cluster = TRUE ;
			  /* group = ap1->group ; */
			  abiFixLabelClusterPolyAPair (aa3, ap, ap1) ;
			}
		      else
			break ;
		    }
	      }
	  if (sortNeeded)
	    arraySort (aa3, aaaTgOrder) ; /* sort again since we shifted */
	}
      
      if (arrayMax (aa3) > 1 && cluster)
	{
	  /* relocalize clusters on their best representative support */    
	  int a20, n, bestA20, bestN, oldGroup = -1 ;
	  
	  /* sort by cluster, then by original position */
	  arraySort (aa3, aaaTgOrder20) ;
	  
	  for (ii = arrayMax (aa3) - 1, ap = arrp (aa3,ii, AAA) ; ii >= 0 ; ii--, ap--)
	    if (ap->group != oldGroup)
	      {
		oldGroup = ap->group ;
		a20 = bestA20 = ap->a20 ; 
		n = 1 + 5 * ap->nAClones0 + ap->mClones03 ; 
		if (ap->tail) n+=2 ;
		bestN = n ; 
		for (jj = ii-1, ap1 = ap-1 ; jj >= 0 && ap1->group == ap->group ; jj--, ap1--)
		  {
		    if (ap1->a20 != a20) 
		      n = 0 ;
		    /* 2010_02_21, change the bonus of clones from 5*ap1->nAClones0  down to 2*ap1->nAClones0 */
		    n += 1 + 2*ap1->nAClones0 + ap1->mClones03 ; a20 = ap1->a20 ; 
		    if (ap1->tail) 
		      n+=2 ; 
		    if (n > bestN) 
		      { bestN = n ; bestA20 = a20 ; }
		  }
		for (jj = ii, ap1 = ap ; jj >= 0 && ap1->group == ap->group ; jj--, ap1--)
		  ap1->a22 = ap1->a2 = bestA20 ;  /* most supported coord inside the cluster */
	      }
	}
      arraySort (aa3, aaaTgOrder) ;
    } /* end while nRound loop */


  /* check for a signal close to the end */
  if (arrayMax(aa3))
    for (ii = 0, ap = arrp (aa3,ii, AAA) ;  ii < arrayMax (aa3) ; ii++, ap++)
      abiFixHasPolyASignal (ap, gDna) ;

  return ;
}  /* abiFixLabelClusterPolyA */

/*********************************************************************/

static int abiFixLabelReportPolyA (KEY mrna, Array aa5, Array aa3, DICT *dict)
{
  int ii, jj, n1 = 0, a0 = 0, a1, a2, m1, m2, nClones3, mrnaLength = 0 ;
  AAA *ap, *ap1 ;
  KEYSET ks = 0 ;
  BOOL hasOpenProduct5 = FALSE ;
  BOOL hasOpenProduct3 = FALSE ;
  KEY _Valid5p = str2tag ("Valid5p") ;
  KEY _Multi5p = str2tag ("Multi5p") ;
  KEY _Valid3p = str2tag ("Valid3p") ;
  KEY _Multi3p = str2tag ("Multi3p") ;
  OBJ Mrna = 0 ;
  Array units = arrayCreate (256, BSunit) ;
  BSunit *uu ;
  BOOL isDown = TRUE ;
  BOOL isWorm = getPleaseIsWorm () ;
  BOOL UTR_showAll = TRUE ; /* 2009_06_01 : should be FALSE except to export to KRIS a maximal number of polya flags */
  BOOL CHEAT = FALSE   ;  /* if TRUE: special trick to export for kris the global clustering of her polya flags */

  ks = queryKey (mrna, ">product Best_product && Good_product && !NH2_complete") ;
  if (keySetMax (ks))
    hasOpenProduct5 = TRUE ;
  keySetDestroy (ks) ;
  ks = queryKey (mrna, ">product Best_product && Good_product && !COOH_complete") ;
  if (keySetMax (ks))
    hasOpenProduct3 = TRUE ;
  keySetDestroy (ks) ;
  Mrna = bsUpdate (mrna) ;

  /* locate the mRNA on the cosmid */
  a1 = a2 = m1 = m2 = 0 ;
  if (bsGetArray (Mrna, _Covers, units, 5))
    {
      ii = arrayMax (units) - 5 ;
      uu = arrp (units, ii, BSunit) ; 
      mrnaLength = uu[0].i ; 
      a0 = a1 = uu[2].i ; /* on cosmid */
      a2 = uu[4].i ; 
      if (a1 > a2) isDown = FALSE ;
    }
 
  /* locate the last exon of the mRNA */
  if (bsGetArray (Mrna, _Splicing, units, 5))
    {
      ii = arrayMax (units) - 5 ;
      uu = arrp (units, ii, BSunit) ; 
      if (isDown)
	{ 
	  a2 = a1 + uu[1].i - 1 ; 
	  a1 = a1 + uu[0].i - 1 ; /* on cosmid */
	}
      else
	{
	  a2 = a1 - uu[1].i + 1 ; 
	  a1 = a1 - uu[0].i + 1 ; 
	  a1 = - a1 ; a2 = - a2 ; /* negate coordinates */
	  a0 = - a0 ;
	}
      m1 = uu[2].i ;
      m2 = uu[3].i ;
    }
  if (CHEAT) m2 = a2 = 30000000  ;

  /* report */
  {
    int x ;

    if (bsFindTag (Mrna, _Multi5p))
      bsRemove (Mrna) ;
    if (bsFindTag (Mrna, _Multi3p))
      bsRemove (Mrna) ;
    if (! hasOpenProduct5 && arrayMax(aa5))
      {
	/* locate first position */
	for (ii = arrayMax (aa5) - 1, ap = arrp (aa5,ii, AAA) ; ii >= 0 ; ii--, ap--)
	  if (ap->mrna == mrna) { break ; }
	
	for (ii = 0, ap = arrp (aa5, 0, AAA) ; ii < arrayMax (aa5) ; ii++, ap++)
	  if (ap->est && /* a cluster point */
	      ap->top &&
	      ap->isSl >= 3
	      )
	    { /* report if represented in this mrna */
	      for (jj = ii, ap1 = ap ; jj >= 0 && ap1->a1 == ap->a1 ; jj--, ap1--)
		if (ap1->mrna == mrna &&
		    (ap1->a1 == ap1->ma1 || ap1->a1 >= ap1->ma1 + 50)
		    )
		  {
		    {
		      bsAddData (Mrna, _Valid5p, _Int, &(ap1->m1)) ;
		      bsAddData (Mrna, _bsRight, _Int, &(ap1->x1)) ;
		      bsAddData (Mrna, _bsRight, _Int, &(ap->isSl)) ;
		      n1++ ;
		    }
		    break ;
		  }
	    }
      }
    if (! hasOpenProduct3 && arrayMax(aa3))
      {
	int pass, xTags, nTags ;
	/* locate last position */
	nTags = 0 ;
	for (ii = arrayMax (aa3) - 1, ap = arrp (aa3,ii, AAA) ; ii >= 0 ; ii--, ap--)
	  if (ap->mrna == mrna) { break ; }
	
	/* pass 0: count the total number of tags, 
	 * pass 1: discard position with less than 5% of the tags and no clones
	 */
	for (pass = 0 ; pass < 2 ; pass++)
	  for (ii = 0, ap = arrp (aa3, 0, AAA) ; ii < arrayMax (aa3) ; ii++, ap++)
	    if (ap->est && /* a cluster point */
		ap->aRich <= 0 &&
		(ap->signal || isWorm || ap->mClones03 > 10) &&
		(ap->nClones3 >= (ap->polyAPriming ? 1 : 3) || ap->tail || ap->mClones3 > (UTR_showAll ? 1 : 3))
		)
	      { /* report if represented in this mrna, or if matching the end of the mrna */
		/* count the supporting reads */
		BOOL atEnd = FALSE ;
		
		xTags = ap->nClones3 ;  
		if (pass == 0) 
		  {
		    nTags +=  xTags ;
		    continue ;
		  }
		if (xTags < nTags/20)
		  continue ;
		if (xTags == 1 && nTags > 10)
		  continue ;
		nClones3 = xTags ;
		if (!nClones3)
		  continue ;
		if  (!(ap->nClones3 >= (ap->polyAPriming ? 1 : 3)) && ! ap->tail && ap->nClones3 < 2 && !atEnd)
		  continue ; 
		x = ap->a22 - a0 + 1 ;
		/* reject flags stolen far downstream if we alrerady have one */
		if (x - mrnaLength > 10 && n1) 
		  continue ;
		/* reject flags stolen far downstream even if we do not yet have a flag */
		if (x - mrnaLength > 50) 
		  continue ;
		if (x < 0 || m2 + ap->a22 - a2 < 0) 
		  continue ;
		if (1)
		  { /* check we are inside the last exon */
		    BOOL ok = FALSE ;
		    if (bsGetArray (Mrna, _Splicing, units, 5))
		      {
			for (int i = arrayMax (units) - 5 ; !ok && i >= 0 && i < arrayMax (units) ; i += 5)
			  {
			    uu = arrp (units, i, BSunit) ; 
			    if (uu[4].i == _Exon &&
				uu[0].i <= x &&
				uu[1].i >= x)
			      ok = TRUE ;
			  }
		      }
		    if (ok)
		      {
			bsAddData (Mrna, _Valid3p, _Int, &x) ;
			x = m2 + ap->a22 - a2 ;
			bsAddData (Mrna, _bsRight, _Int, &x) ;
			
			bsAddData (Mrna, _bsRight, _Int, &nClones3) ;
			
			if (ap->signal && *ap->buffer)
			  {
			    if (1)
			      {
				bsAddData (Mrna, _bsRight, _Text, ap->buffer) ;
			      }
			    else
			      {
				if (ap->signal == 1)
				  bsAddData (Mrna, _bsRight, _Text, "AATAAA") ;
				else
				  bsAddData (Mrna, _bsRight, _Text, "Variant") ;
			      }
			    
			    x = ap->aSignal  ;
			    bsAddData (Mrna, _bsRight, _Int, &x) ;
			  }
			else
			  {
			    bsAddData (Mrna, _bsRight, _Text, "No_signal") ;
			    x = 0 ;
			    bsAddData (Mrna, _bsRight, _Int, &x) ;
			  }
		      }
		  }
		n1++ ;
	      }
      }
  } /* if products) {} */  
  
  bsSave (Mrna) ;
  arrayDestroy (units) ;
  
  return n1 ;
} /* abiFixLabelOnePolyA */

/*********************************************************************/

static void abiFixLabelRestoreManualTags (KEY tg)
{
  KEYSET ks ;
  int ii ;
  OBJ Clo ;
  
  ks = queryKey (tg, ">cdna_clone Internal_priming_manual && !Internal_priming") ;
  for (ii = 0 ; ii < keySetMax (ks) ; ii++)
    {
      if ((Clo = bsUpdate (keySet (ks, ii))))
	{
	  bsAddTag (Clo, _Internal_priming) ; 
	  bsSave (Clo) ; 
	}
    }
  keySetDestroy (ks) ;
} /* abiFixLabelRestoreManualTags */

/*********************************************************************/
/*
  $tacembly . <<EOF | tee -a d4.reverse.prelog
  table -o d4.polyAsuspect1.txt  $ici/tables/10.polyAsuspect1.def
  table -o d4.polyAsuspect2.txt  $ici/tables/10.polyAsuspect2.def
  table -o d4.polyAsuspect3.txt  $ici/tables/10.polyAsuspect3.def
  table -o d4.polyAsuspect4.txt  $ici/tables/10.polyAsuspect4.def
  quit
  EOF
  # suspect 1 finds reverse read starting inside the ORF
  # suspect 2 finds forward polyA read endding inside the ORF
  # suspect 3 marks as internal priming 3' clones reverse with no stop in product
  # suspect 4 removes polyA of forward read if a reverse read assembles further down
  
  $myawk  -F '\t' '/\"/ {printf ("Sequence %s\nGene_wall\n-D polyA_after_base\n\n",$9);}'  d4.polyAsuspect1.txt >!  d4.polyAsuspect1.ace
  $myawk  -F '\t' '/\"/ {printf ("Sequence %s\nGene_wall\n-D polyA_after_base\n\ncDNA_clone %s\nInternal_priming\n\n",$9, $10);}'  d4.polyAsuspect2.txt >!  d4.polyAsuspect2.ace
  $myawk  -F '\t' '/\"/ {printf ("cDNA_clone %s\nInternal_priming\n\nSequence %s\nGene_wall\n\n",$5,$4);}'  d4.polyAsuspect3.txt >!  d4.polyAsuspect3.ace
  $myawk  -F '\t' '/\"/ {printf ("Sequence %s\nGene_wall\n-D polyA_after_base\n\n",$4);}'  d4.polyAsuspect4.txt >!  d4.polyAsuspect4.ace
  
  #endif
*/  
/*********************************************************************/

int abiFixLabelPolyATg (KEY tg, int *nClonep, int *n0p, int *nMrnap)
{
  KEYSET mrnas = queryKey (tg, "{CLASS mRNA} SETOR {>Mrna}") ;
  int ii, n1 = 0 ;
  Array gDna = 0,  aa5 = arrayCreate (300, AAA), aa3 = arrayCreate (300, AAA) ;
  KEY cosmid = 0 ;
  KEYSET clones5 = keySetCreate () ;
  KEYSET clones3 = keySetCreate () ;
  DICT *dict = dictCreate (100) ;
  BOOL isDown = TRUE ;
  
  abiFixLabelRestoreManualTags (tg) ;
  if ((cosmid = keyGetKey (tg, _Genomic_sequence)))
    gDna = dnaGet (cosmid) ;
  if (!gDna)
    goto done ;

  for (ii = 0 ; ii < keySetMax (mrnas) ; ii++)
    isDown = abiFixLabelGatherPolyA (keySet (mrnas, ii), gDna, aa5, aa3, clones5, clones3, dict) ;
  abiFixLabelClusterPolyA (aa5, aa3, gDna, isDown, nClonep) ;
  for (ii = n1 = 0 ; ii < keySetMax (mrnas) ; ii++)
    n1 += abiFixLabelReportPolyA (keySet (mrnas, ii), aa5, aa3, dict) ;

 done:
  if (n0p) *n0p += keySetMax (mrnas) ;
  if (nMrnap) *nMrnap += n1 ;

  arrayDestroy (aa5) ;
  arrayDestroy (aa3) ;
  keySetDestroy (mrnas) ;
  keySetDestroy (clones5) ;
  keySetDestroy (clones3) ;
  dictDestroy (dict) ;
  arrayDestroy (gDna) ;

  return n1 ;
} /* abiFixLabelPolyATg */

/*********************************************************************/

int abiFixLabelPolyA (KEYSET ks0, int *nClonep, int *n0p, int *nMrnap)
{
  KEYSET tgs = query (ks0, "{CLASS transcribed_gene} SETOR {CLASS mRNA}") ;
  int ii ;

  *nClonep = *n0p = *nMrnap = 0 ;
  for (ii = 0 ; ii < keySetMax (tgs) ; ii++)
    abiFixLabelPolyATg (keySet (tgs, ii), nClonep, n0p, nMrnap) ;

  keySetDestroy (tgs) ;

  return *nMrnap ;
} /* abiFixLabelPolyA */

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
