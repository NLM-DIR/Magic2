/*  File: basecall.c
 *  Author: Jean Thierry-Mieg (mieg@kaa.crbm.cnrs-mop.fr)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 19 13:31 1998 (fw)
 * Created: Fri Dec 23 16:43:01 1994 (mieg)
 *-------------------------------------------------------------------
 */

/* @(#)basecall.c	1.15 5/23/97 */
#include "acedb.h"
#include "acembly.h"
#include "basecall.h"
#include "a.h"
#include "lex.h"
#include "tags.h"
#include "classes.h"
#include "session.h"
#include "freeout.h"
#include "plot.h"
#include "Read.h"   /* Staden's */
#include "ctf.h" 
#include "myNetwork.h"
#include "saucisse.h"
#include "seqIOCTF.h"
#include "acembly.h"  


int seqLeftCutoff () ;
int seqRightCutoff () ;

static void seq2ace(LANE *lane) ;
void findXclipping (LANE *lane) ;
static BOOL baseCallFindClips (LANE *lane, BOOL doPLot) ;

#ifdef NON_GRAPHIC
#include "graph.h"
static void fMapTraceForget(KEY key) { return ; }
#else
extern  void fMapTraceForget(KEY key) ;
#endif

/***************************************************************************/
/***************************************************************************/

BOOL baseCallGetCafEdits (LANE *lane) 
{ 
  Array 
    aa = 0 , dna = 0 , 
    base = 0, newBase = 0, basePos = 0, newBasePos = 0 ;
  OBJ obj = 0 ;
  int x, dx, ii, i, j, i0, j0, i1, i2, j1, j2 ;
  char c1, c2 ;
  BSunit *u ;

  if (!lane->key || !lane->dna ||
      !(obj = bsCreate(lane->key)))
    return FALSE ;
  
  aa = arrayCreate (400, BSunit) ;
  if (!bsGetArray (obj, _Align_to_SCF, aa,4))
    { bsDestroy (obj) ;
      return FALSE ;
    }

  bsDestroy (obj) ;

  dna = lane->dna ;
  base = lane->base ;
  i = arrayMax(lane->base) + 20 ;
  newBase = arrayCreate (i, char) ;
  basePos = lane->basePos ;
  newBasePos = arrayCreate (i, short) ;
  i0 = j0 = 0 ;
  for (ii = 0 ; ii < arrayMax (aa) ; ii += 4)
    {  u = arrp (aa, ii, BSunit) ;
       i1 = u[0].i - 1 ;
       i2 = u[1].i - 1 ;
       j1 = u[2].i - 1 ;
       j2 = u[3].i - 1 ;
       if (j2 - j1 != i2 - i1 ||
	   j2  < j1)
	 { messerror ("Align_to_scf %s error %d %d %d %d",
		      name(lane->key), i1, i2, j1, j2) ;
	   goto abort ;
	 }
       x = array (basePos, j0, short) ;
       dx = array (basePos, j1, short) - x ;
       if (i1 > i0) dx /= (i1 - i0) ;
       for (i = i0 ; i < i1 ; i++)
	 { array (newBasePos, i, short) = x ;
	   x += dx ;
	   c1 = array (dna, i, char) ;
	   array (newBase, i, char) = c1 | BC_ACE ;
	 }
       for (i = i1 , j = j1 ; i <= i2 ; i++, j++)
	 { array (newBasePos, i, short) = array (basePos, j, short) ;
	   c1 = array (dna, i, char) ;
	   c2 = array (base, j, char) ;
	   if ((c1 & 0xf) != (c2 &0xf))
	     c2 = c1 | BC_ACE ;
	   array (newBase, i, char) = c2 ;
	 }
       i0 = i2 + 1 ; j0 = j2 ;
     }
  arrayDestroy (aa) ;
  lane->basePos = newBasePos ;
  arrayDestroy (basePos) ;
  lane->base = newBase ;
  arrayDestroy (base) ;

  return TRUE ;
 abort:
  arrayDestroy (aa) ;
  arrayDestroy (newBase) ;
  arrayDestroy (newBasePos) ; 
  bsDestroy (obj) ;
  return FALSE ;
}

/***************************************************************************/
BOOL baseCallDumpDif (KEY key) 
{ 
/*
  Array bcA1 = 0, bcA2 = 0 ;
  BOOL done = FALSE ;
  OBJ obj = 0 ;
  KEY bcKey1, bcKey2 ;
  int n1, n2 ;
  int n, x ; 
  char b, ba, *bc, buf[256] ;

  if (!key || 
      !(obj = bsCreate(key)))
    return FALSE ;

  if (!bsGetKey (obj, _BaseCall, &bcKey1) ||
      !bsGetKey (obj, _OriginalBaseCall, &bcKey2))
    goto abort ;
   
  bcA1 = arrayGet (bsKey1, char, "c") ;
  bcA2 = arrayGet (bcKey2, char, "c") ;
  if (!bc1 || !bc2)
    goto abort ;

  n1 = arrayMax (bc1) ; n1 /= 2 ;
  n2 = arrayMax (bc2) ; n2 /= 2 ;
  x1 = x2 = 0 ;
      bc = arrp (bc1, 0, char) ;
      while (n--)
	{ x += *bc++ ;
	  b = *bc++ ;
	  ba = ace_upper(dnaDecodeChar [b &0xf]) ;
	  if (b & BC_LOW) ba = ace_lower (ba) ;
	  sprintf (buf, "%6d %c", x, ba) ;
	  if (b & BC_HAND) 
	    strcat (buf, " h") ;
	  if (b & BC_ACE) 
	    strcat (buf, " a") ;
	  strcat (buf, "\n") ;
	  freeOut (buf) ;
	}
  return TRUE ;
 abort:
  arrayDestroy (bcA1) ;
  arrayDestroy (bcA2) ;
  bsDestroy (obj) ;
  return done ;
*/
 return TRUE ;
}

/***************************************************************************/

BOOL baseCallDump (FILE* f, Stack s, KEY key) 
{ Array compressedBasecall = 0, quality = 0 ;
  int n, x, q ; 
  char  ba, buf[256] ;
  unsigned char b, *bc, *bq = 0 ;
  KEY qualityKey = 0 ;

  if (key)
    compressedBasecall = arrayGet (key, char, "c") ;
  if (!compressedBasecall || 
      !arrayMax(compressedBasecall))
    { arrayDestroy (compressedBasecall) ;
      return FALSE ;
    }

  if (lexReClass (key, &qualityKey, _VBaseQuality))
    quality = arrayGet (qualityKey, char, "c") ;
  if (quality && 
      arrayMax(quality) == arrayMax(compressedBasecall))
    bq = arrp (quality, 0, unsigned char) ;
  else
    bq = 0 ;
  if (compressedBasecall)
    { n = arrayMax (compressedBasecall) ;
      n /= 2 ;
      x = 0 ;
      bc = arrp (compressedBasecall, 0, unsigned char) ;
      while (n--)
	{ x += *bc++ ;
	  b = *bc++ ;
	  ba = ace_upper(dnaDecodeChar [b &0xf]) ;
	  if (b & BC_LOW) ba = ace_lower (ba) ;
	  
	  if (bq)
	    { q = *bq++ ;  
	      sprintf (buf, "%6d %c %5d ", x, ba, q) ;
	    }
	  else
	    sprintf (buf, "%6d %c", x, ba) ;

	  if (b & BC_HAND) 
	    strcat (buf, " h") ;
	  if (b & BC_ACE) 
	    strcat (buf, " a") ;
	  strcat (buf, "\n") ;
	  if (f)
	    fprintf (f, "%s", buf) ;
	  else if (s)
	    catText (s, buf) ; 
	  else 
	    freeOut (buf) ; 
	}
    }
  arrayDestroy (compressedBasecall) ;
  arrayDestroy (quality) ;
  return TRUE ;
}

/***************************************************************************/

BOOL baseCallParse (int level, KEY key) 
{ Array aa = arrayCreate (1000, char) ;
  Array quality = arrayCreate (1000, unsigned char) ;
  int n = 0, dx, x, x0 = 0 , nq = 0, nbase = 0 , q ; 
  char c, c1, *cp, flag ;
  KEY qualityKey = 0, seqKey = 0 , baseCallKey = key ;
  OBJ obj ;

  while (freecard(level))
    { flag = 0 ;
      if (!freeint(&x))
	{ if (!freeword())
	    break ;
	  messerror ("  basecallparse error at line %7d in %.25s : no position\n", 
		     freestreamline(level), name(key)) ;
	  goto abort ;
	}
      freenext() ;
      if (!(cp = freeword())) 
	{ messerror ("  basecallparse error at line %7d in %.25s, posotion %d  : no basecall\n", 
		     freestreamline(level), name(key), x) ;
	  goto abort ;
	}
      c = *cp ;
      if (c == ace_lower (c))
        flag |= BC_LOW ;
      c1 = dnaEncodeChar[((int)c) & 0x7f]  ;      /* accepted base codes */
      if (!c1)
	messerror (
		   "  baseCallparse error at line %7d in %.25s : bad char %d:%c\n", 
		   freestreamline(level), name(key), (int)c, c) ;
      freenext() ;
      nbase++ ;
      if (freeint (&q) && q >=0 && q < 255)
	array(quality, nq++, unsigned char) = q ;
      if ((cp = freeword()))
	switch (*cp)
	  {
	  case 'h': flag |= BC_HAND ; break ;
	  case 'a': flag |= BC_ACE ; break ;
	  }
      dx = x - x0 ;
    /* these will smoothe away after a few steps */
      if (dx <= 0) dx = 1 ;
      if (dx > 254) dx = 254 ;
      x0 += dx ;
      array(aa,n++,char) = dx ;
      array(aa,n++,char) = c1 | flag ;
    }
  arrayStore (key, aa, "c") ;
  arrayDestroy (aa) ;

  if (nq == nbase)
    { lexaddkey (name(key), &qualityKey, _VBaseQuality) ;
      arrayStore (qualityKey, quality, "c") ;
    }

  lexaddkey (name(key), &seqKey, _VSequence) ;
  if ((obj = bsUpdate (seqKey)))
    { bsAddKey (obj, _BaseCall, baseCallKey) ;
   /*    bsAddData (obj, _bsRight, _Int, &nbase) ; */
      if (nq == nbase)
	{ bsAddKey (obj, _Quality, qualityKey) ;
	  bsAddData (obj, _bsRight, _Int, &nq) ;
	}
      bsSave (obj) ;
    }
  arrayDestroy (quality) ;
  return TRUE ;

 abort:
  arrayDestroy (aa) ;
  arrayDestroy (quality) ;
  return FALSE ;
}

/***************************************************************************/
/***************************************************/

static void baseCallDoStore (LANE *lane, BOOL reserve)
{ KEY 
    key = lane->key , 
    baseCallKey = lane->baseCallKey ,
    basePositionKey = lane->basePositionKey ,
    baseQualityKey = lane->baseQualityKey ;
  OBJ obj ;
  Array base = lane->base, basePos = lane->basePos, baseCall = 0,
  baseQuality = lane->baseQuality, 
  pp1 = 0 ;
  int n ; short x ,dx, *bp ;
  char *b ; unsigned char *bc ;

  if (!arrayExists(base) || 
      !arrayExists(basePos) ||
      arrayMax(base) != arrayMax(basePos))
    messcrash("Problem in basecallStore") ;

  if (!baseCallKey)
    { obj = bsUpdate(key) ;
      if (!obj)
	return ; 
      lexaddkey (name(lane->key), &baseCallKey, _VBaseCall) ;
      if (reserve)
	{ KEY res ;
	  lexaddkey("OriginalBaseCall", &res, 0) ;
	  bsAddKey (obj, res, baseCallKey) ;
	}
      else
	{ lane->baseCallKey = baseCallKey ;
	  bsAddKey (obj, _BaseCall, baseCallKey) ;
	}
      bsSave (obj) ;
    }

  if (!basePositionKey)
    { obj = bsUpdate(key) ;
      if (!obj)
	return ; 
      lexaddkey (name(lane->key), &basePositionKey, _VBasePosition) ;
      lane->basePositionKey = basePositionKey ;
      bsAddKey (obj, _SCF_Position, basePositionKey) ;
      bsSave (obj) ;
    }

  if (!baseQualityKey)
    { obj = bsUpdate(key) ;
      if (!obj)
	return ; 
      lexaddkey (name(lane->key), &baseQualityKey, _VBaseQuality) ;
      lane->baseQualityKey = baseQualityKey ;
      bsAddKey (obj, _Quality, baseQualityKey) ;
      bsSave (obj) ;
    }

  n = arrayMax (base) ;
  baseCall = arrayCreate (2*n, char) ;
  array (baseCall, 2*n - 1, char) = 0 ;
  x = 0 ;
  b = arrp (base, 0, char) ;
  bp = arrp (basePos, 0, short) ;
  bc = arrp (baseCall, 0, unsigned char) ;
  while (n--)
    { dx = *bp++ - x ;
    /* these will smoothe away after a few steps */
      if (dx <= 0) dx = 1 ;
      if (dx > 254) dx = 254 ;
      *bc++ = dx ;
      x += dx ; 
      *bc++ = *b++ ;
    }

  arrayStore (baseCallKey, baseCall,"c") ;
  arrayDestroy (baseCall) ;

  n = arrayMax (basePos) ;
  pp1 = arrayCreate (n, char) ;
  array (pp1, n - 1, char) = 0 ;
  x = 0 ;
  bp = arrp (basePos, 0, short) ;
  bc = arrp (pp1, 0, unsigned char) ;
  while (n--)
    { dx = *bp++ - x ;
    /* these will smoothe away after a few steps */
      if (dx <= 0) dx = 1 ;
      if (dx > 127) dx = 127 ;
      *bc++ = dx ;
      x += dx ; 
    }
  arrayStore (basePositionKey, pp1,"c") ;
  arrayDestroy (pp1) ;

  if (baseQuality)
    arrayStore (baseQualityKey, baseQuality,"c") ;
}

/***************************************************************************/
/***************************************************************************/
/*********************************************************/

void laneDestroy (LANE *lane)
{ arrayDestroy (lane->errArray) ;
  arrayDestroy (lane->dna) ;
  arrayDestroy (lane->base) ;
  arrayDestroy (lane->basePos) ;
  arrayDestroy (lane->baseCall) ;
  arrayDestroy (lane->tags) ;
  freeSeq (lane->seq) ;
  memset(lane, 0, sizeof(struct LaneStruct)) ;
}

/*********************************************************/

/* problem: 377 abi machines byte saturate, resulting in silly drops
   of 255, i hack any drop of 150 followed soon by a rise of 150
   as a byte problem 
   */
static void hackSaturation (Read* seq)
{ int ii, i, j, nn = 0 ;
  TRACE *bp[4], *u, *v ;

  bp[0] = seq->traceA ;
  bp[1] = seq->traceC ;
  bp[2] = seq->traceG ;
  bp[3] = seq->traceT ;

  ii = 4 ; 
  while (ii--)
    { u = bp[ii] ;
      i = seqMax(seq)/2 ;
      nn = 0 ;
      while (u++, i--)
	{ 
	  if (!nn && *u > *(u - 1) + 30 && *(u + 1) + 120 < *u)
	    { j = 30 ; v = u ;
	      while (v++, j--)
		if (*v < *(v-1) - 30 && *(v+1) > *v + 120)
		  { nn = 1 ; *u -= 256 ; break ; }
	    }
	  else if (nn && *u < *(u - 1) - 30 && *(u + 1) > *u + 120)
	    { *u += 256 ; nn = 0 ; }
	  if (nn) *u += 256 ;
	}
    }
}

/*********************************************************/

static void renormalize (LANE *lane, BOOL doit)
{ int ii, i ;
  TRACE *bp[4], *sp, hMin = 0, hMax = 0 ;
  float f ;
  Read* seq = lane->seq ;
  int max = seqMax (seq) ;

  bp[0] = seq->traceA ;
  bp[1] = seq->traceC ;
  bp[2] = seq->traceG ;
  bp[3] = seq->traceT ;

  if (! doit) /*  un renormalize */
    {
      f = lane->traceRenormalize ;
      hMin = lane->traceShift ;
      if (f != 0) /* otherwise nothing to do */
	{ 
	  for (ii=0 ; ii<4 ; ii++)
	    { 
	      sp = bp[ii] ;
	      for (i = 0; i < max ; i++, sp++)
		*sp = (*sp)/f + hMin ;
	    }
	   lane->traceRenormalize = 0 ; /* avoid recursion */
	}
    }
  else  /* do renormalize */
    {
      hMin = hMax = 0 ;
      for (ii=0 ; ii<4 ; ii++)
	{ sp = bp[ii] ;
	for (i = 0; i < max ; i++, sp++)
	  { if (hMax < *sp) hMax = *sp ;
	  if (hMin > *sp) hMin = *sp ;
	  }
	}
      
      f = hMax - hMin ; 
      if (f > 16) /* otherwise it is too silly */
	{ 
	  f = 254.0/f ;
	  for (ii=0 ; ii<4 ; ii++)
	    { 
	      sp = bp[ii] ;
	      for (i = 0; i < max ; i++, sp++)
		*sp = (*sp - hMin) * f ;
	    }
	  lane->traceRenormalize = f ;
	  lane->traceShift = hMin ;
	}
      else
	lane->traceRenormalize = 0 ;
    }
}

/*********************************************************/

static BOOL doGetTrace (LANE *lane, char *myname, char style)
{ 
  FILE *ff = 0 ;

  if (lane->seq)
    freeSeq (lane->seq) ;
  
  if (myname)
    ff = fopen (myname, "rb") ;
  else
    { ff = stdin ; myname = "stdin" ; }
  if (!ff)
    return FALSE ;
  
  if (style == 's' || style == 'c')
    lane->seq = fread_reading (ff, myname, 0) ;
  
  if (myname && strcmp(myname,"stdin") && ff) fclose (ff) ;
  
  if (!lane->seq)
    return FALSE ;
  
  if (seqMax(lane->seq) <= 0)
    { messout ("Sorry, trace file %s contains a sequence of length 0",
	       myname) ;
    freeSeq (lane->seq) ;
    return FALSE ;
    }
  hackSaturation (lane->seq) ;
  renormalize (lane, TRUE) ;
  
  return TRUE ;
}

/*********************************************************/
  
static BOOL baseCallGetScfSeq (LANE *lane)
{
  char *cp, *cq, *myname = 0 ;

  if (lane->key == 1)   /* scf2ctf */
    return doGetTrace (lane, 0, 's') ;

  cq = baseCallNameGuess (lane->key, 1) ;
  if ((cp = getenv ("SCF_DATA")))
    { myname = filName (messprintf ("%s/%s", cp, cq), 0, "r") ;
      if (!myname)
	myname = filName (messprintf ("%s/%s", cp, cq), "scf", "r") ;
    }
  if (!myname)
    myname = filName (messprintf ("SCF/%s", cq), 0, "r") ;
  if (!myname)
    myname = filName (messprintf ("SCF_DATA/%s", cq), "scf", "r") ;
  
  if (!myname || !doGetTrace (lane, myname, 's'))
    return FALSE ;

  return TRUE ;
}

/*********************************************************/
  
static BOOL baseCallGetCtfSeq (LANE *lane)
{
  char *cp, *cq, *myname = 0 ;

  if (lane->key == 2)   /* ctf2scf */
    return doGetTrace (lane, 0, 'c') ;

  cq = baseCallNameGuess (lane->key, 2) ;
  if ((cp = getenv ("CTF_DATA")))
    { myname = filName (messprintf ("%s/%s", cp, cq), 0, "r") ;
      if (!myname)
	myname = filName (messprintf ("%s/%s", cp, cq), "ctf", "r") ;
    }
  if (!myname)
    myname = filName (messprintf ("CTF/%s", cq), "ctf", "r") ;
  if (!myname)
    myname = filName (messprintf ("CTF/%s", cq), 0, "r") ;
  
 
 if (!myname || !doGetTrace (lane, myname, 'c'))
    return FALSE ;
  return TRUE ;
}

/*********************************************************/

BOOL baseCallGetSeq (LANE *lane)
{ 
  if (lane->scf > 2)
    return lane->scf > 3 ? TRUE : FALSE ;

  if (lane->key != 1 && baseCallGetCtfSeq (lane))
    { lane->scf = 5 ; return TRUE ; }
  if (lane->key != 2 && baseCallGetScfSeq (lane))
    { lane->scf = 4 ; return TRUE ; }
  else
    { lane->scf = 3 ; return FALSE ; }
}

/*********************************************************/
/* construct base and basePos using dna, useful for sequences with no available trace */
BOOL baseCallFakePos (LANE *lane)
{
  Array 
    base = 0, basePos = 0 ;
  int n ; unsigned short x , *bp ;
  char *b ; unsigned char *b1 ;

  if (lane->base) return TRUE ;
  if (lane->scf < 3 || !lane->dna) return FALSE ;

  n = arrayMax (lane->dna) ;
  
  base = arrayCreate (n, char) ;
  array (base, n - 1, char) = 0 ;
  basePos = arrayCreate (n, short) ;
  array (basePos, n - 1, short) = 0 ;
  b = arrp (base, 0, char) ;
  bp = arrp (basePos, 0,  unsigned short) ;
  b1 = arrp (lane->dna, 0, unsigned char) ;
  x = -9 ;
  while (n--)
    { 
      x += 10 ;
      *bp++ = x ;
      *b++ = *b1++ ;
    }
  lane->basePos = basePos ;
  lane->base = base ;
  lane->maxPos = x ;
  
  return TRUE ;
}

/**********************************************************/
  /* Extract DNA from SCF files */
BOOL ctfExportKey (KEY key)
{ 
  LANE *lane = (LANE*) messalloc (sizeof (struct LaneStruct)) ;
  BOOL done = FALSE ;
  FILE *ff = 0 ;

  cDNAAlignInit () ;

  lane->key = key ;   /* seq == 1: scf2ctf from stdin */
                      /* seq == 2: ctf2scf from stdin */ 

  if (key != 1)
    {
      ff = filopen (messprintf("CTF/%s", name(key)), "ctf", "wb") ;
      if (!ff)
	return FALSE ;
    }
  else
    ff = stdout ;

  if (baseCallGet (lane) &&
      baseCallGetSeq(lane) &&  
      lane->scf >= 4)
    {
      done =  fwrite_reading (ff, lane->seq, TT_CTF) == 0 ? TRUE : FALSE ; /* staden convention */
    }
  laneDestroy (lane) ;
  messfree (lane) ;
  if (key != 1)
    filclose (ff) ;

  return done ;
}

/**********************************************************/
  /* Extract DNA from SCF files */
BOOL scfExportKey (KEY key)
{ 
  LANE *lane = (LANE*) messalloc (sizeof (struct LaneStruct)) ;
  BOOL done = FALSE ;
  FILE *ff = 0 ;

  cDNAAlignInit () ;

  lane->key = key ;   /* seq == 1: scf2scf from stdin */
                      /* seq == 2: scf2scf from stdin */ 

  if (key != 1)
    {
      ff = filopen (messprintf("SCF/%s", name(key)), "scf", "wb") ;
      if (!ff)
	return FALSE ;
    }
  else
    ff = stdout ;

  if (baseCallGet (lane) &&
      baseCallGetSeq(lane) &&  
      lane->scf >= 4)
    {
      done = fwrite_reading (ff, lane->seq, TT_SCF) == 0 ? TRUE : FALSE ; /* staden convention */
    }
  laneDestroy (lane) ;
  messfree (lane) ;
  if (key != 1)
    filclose (ff) ;

  return done ;
}

/**********************************************************/

void scfTraceExportKeySet (KEYSET ks)
{
  int i ; 
  if (!keySetExists(ks) || !keySetMax(ks))
    return ;
  
  for (i = 0 ; i < keySetMax(ks) ; i++)
    scfExportKey (keySet(ks,i)) ;
}

/**********************************************************/

  /* Extract DNA from SCF files */
static BOOL ctfFillSaucisse (Saucisse saucisse, KEY seq, int predictionMode)
{ 
  LANE *lane = (LANE*) messalloc (sizeof (struct LaneStruct)) ;
  BOOL done = FALSE ;
  Array a = 0, b = 0 ;
  int n ;
  short s, *sp ; unsigned char *cp ;
  extern Array ctfDecorrelate (Read *read, int predictionMode) ; /* do not put in .h, this is a  test code */

  cDNAAlignInit () ;

  lane->key = seq ;
  if (baseCallGet (lane) &&
      baseCallGetSeq(lane) &&  
      lane->scf >= 4)
    {
      a = ctfDecorrelate (lane->seq,predictionMode) ;
      n = arrayMax (a) ;
      b = arrayCreate (n, unsigned char) ;
      array (b, n - 1, unsigned char) = 0 ; 
      cp = arrp (b, 0, unsigned char) ;
      sp = arrp (a, 0, short) ;
      while (n--)
	{ 
	  s = *sp++ + 128 ;
	  *cp++ = s > 0 && s < 255 ? s : 255 ; 
	}
      saucisseFill (saucisse, b) ;
      arrayDestroy (a) ;
      arrayDestroy (b) ;
    }
  laneDestroy (lane) ;
  messfree (lane) ;
  return done ;
}

/**********************************************************/

void ctfTraceExportKeySet (KEYSET ks)
{
  int i ; 
  Saucisse saucisse = 0 ;
  int mode = 0 ; /* 0: export,  1: profile */
  int predictionMode = -1 ;  /* auto choice in decorelattor */

  if (!keySetExists(ks) || !keySetMax(ks))
    return ;
  
  switch (mode)
    {
    case 0: /* normal operations */
      for (i = 0 ; i < keySetMax(ks) ; i++)
	ctfExportKey (keySet(ks,i)) ;
      break ;

    case 1: /* saussage mode, to profile the compressor */
      saucisse =  saucisseCreate (2000, 8, FALSE) ;
      for (i = 0 ; i < keySetMax(ks) ; i++)
	ctfFillSaucisse (saucisse, keySet(ks,i), predictionMode) ;
      saucisseShow (saucisse, 1) ;
      saucisseDestroy (saucisse) ;
      break ;
    }  
}

/***************************************************************************/

static BOOL doNotStore = FALSE ;

BOOL baseCallGet (LANE *lane)
{ KEY 
    dnaKey = 0,
    baseCallKey = lane->baseCallKey ,
    positionKey = lane->basePositionKey ,
    qualityKey = lane->baseQualityKey ;
  Array 
    base = 0, basePos = 0, compressedBasecall = 0 ,
    pp1 = 0, qq1 = 0, dna1 ;
  int n ;  short x , y, *bp ; 
  unsigned short *bp2 ;
  char *b, *b1, c, *bp1 ;
  unsigned char *bc ;

  if (lane->base)
    return TRUE ;
  if (baseCallKey)
    compressedBasecall = arrayGet (baseCallKey, char, "c") ;

  if (compressedBasecall)
    { n = arrayMax (compressedBasecall) ;
      n /= 2 ;
      base = arrayCreate (n, char) ;
      array (base, n - 1, char) = 0 ;
      basePos = arrayCreate (n, short) ;
      array (basePos, n - 1, short) = 0 ;
      x = 0 ;
      b = arrp (base, 0, char) ;
      bp = arrp (basePos, 0, short) ;
      bc = arrp (compressedBasecall, 0, unsigned char) ;
      while (n--)
	{ x += *bc++ ;
	  *bp++ = x ;
	  *b++ = *bc++ ;
	}
    }
  else if (qualityKey && positionKey && dnaKey &&
	   (pp1 = arrayGet (positionKey, char, "c")) &&
	   (qq1 = arrayGet (qualityKey, unsigned char, "c")) &&
	   (dna1 = arrayGet (dnaKey, char, "c")))
    { n = arrayMax (dna1) ;
      if (!n)
	{ arrayDestroy (dna1) ;
	  arrayDestroy (pp1) ;
	  arrayDestroy (qq1) ;
	  return FALSE ;
	}
      basePos = arrayCreate (n, short) ;
      array (basePos, n - 1, short) = 0 ;
      b = arrp (base, 0, char) ;
      bp = arrp (basePos, 0, short) ;
      bp1 = arrayp (pp1, 0, char) ;
      x = - 1000 ; y = 0 ;
      while (n--)
	{ y += *bp1 ;
	  if (x < y)
	    x = y ;
	  else
	    x++ ;
	  *bp++ = x ;
	}
	 
      base = dna1 ; dna1 = 0 ;
      b = arrp (base, 0, char) ;
      bp1 = arrayp (qq1, 0, char) ;
      n = arrayMax (qq1) ;
      while (n--)
	{ if (*bp1++ < 10)
	    *b |= BC_LOW ;
	  b++ ;
	}
      arrayDestroy (pp1) ;
      arrayDestroy (qq1) ;
    }
  else
    { if (!baseCallGetSeq (lane))
	return FALSE ;
      n = seqMaxBase(lane->seq) ;
      if (!n)
	{ freeSeq (lane->seq) ;
	  return FALSE ;
	}
      if (TRUE) /* exists if read from CTF format */
	{
	  BOOL needConversion = FALSE ;

	  switch (lane->scf)
	    {
	    case 4: 
	    case 5: needConversion = TRUE ; break ; ;
	    }


	  base = arrayCreate (n, char) ;
	  array (base, n - 1, char) = 0 ;
	  basePos = arrayCreate (n, short) ;
	  array (basePos, n - 1, short) = 0 ;
	  b = arrp (base, 0, char) ;
	  bp = arrp (basePos, 0, short) ;
	  b1 = lane->seq->base ;
	  bp2 = lane->seq->basePos ;
	  

	  x = - 10002 ;
	  while (n--)
	    { if (x < *bp2)
	      x = *bp2++ ;
	    else
	      { x++ ; bp2++ ; }
	    *bp++ = x ;
	    
	    c = *b1++ ;
	    if (needConversion)
	      switch (c)
	      {
	      case 'A': c = A_ ; break ;
	      case 'T': c = T_ ; break ;
	      case 'G': c = G_ ; break ;
	      case 'C': c = C_ ; break ;
	      case 'a': c = BC_LOW | A_ ; break ;
	      case 't': c = BC_LOW | T_ ; break ;
	      case 'g': c = BC_LOW | G_ ; break ;
	      case 'c': c = BC_LOW | C_ ; break ;
	      default : c = BC_LOW | N_ ; break ;
	      }
	    *b++ = c ;
	    }
	  x = seqMax(lane->seq) ;
      
	}
    }
  lane->basePos = basePos ;
  lane->base = base ;
  lane->maxPos = x ;
  
  if (lane->base && lane->dna)
    findXclipping (lane) ;
  if (!doNotStore && !compressedBasecall)
    { baseCallGetCafEdits (lane) ;
      baseCallDoStore (lane, TRUE) ;
    }
  arrayDestroy (compressedBasecall) ;
  return TRUE ;
}

/*****************************************************/

void baseCallStore (LANE *lane)
{ baseCallDoStore (lane, FALSE) ;
}

/*********************************************************/

char* baseCallNameGuess (KEY key, int type)
{ 
  static KEY mainClone = 0 ;
  OBJ obj = bsCreate (key) ;
  static char* buffer = 0 ;
  char *cp, *cq ;

  messfree (buffer) ;

  if (!mainClone)
    { 
      KEYSET ks = query (0,"Find Clone Main_Clone") ;
      if (keySetMax(ks) == 1)
	mainClone = keySet (ks, 0) ;
      else 
	mainClone = 1 ;
      keySetDestroy (ks) ;
    }

  if (!_CTF_File)
    cDNAAlignInit() ;
  if (!obj || 
      (type == 1 && !bsGetData(obj, _SCF_File, _Text, &cp)) ||
      (type == 2 && !bsGetData(obj, _CTF_File, _Text, &cp)))
    {
      cp = name (key) ;
      if (mainClone != 1)
	{
	  int i ;
	  cq = name (mainClone) ;
	  i = strlen (cq) ;
	  if (i && strlen (cp) > i && !strncmp (cp, cq, i) && 
	      *(cp + i) == '.')
	    cp += i + 1 ;
	}
    }
  bsDestroy (obj) ;

  buffer = strnew (cp, 0) ;
  return buffer ;
}
    
/*****************************************************/
/*****************************************************/

Array baseCallBase (KEY key)
{ Array aa ;
  static LANE *lane = 0 ;

  if (!lane)
    lane = (LANE*) messalloc (sizeof (struct LaneStruct)) ;

  lane->key = key ;
  baseCallGet(lane) ;
  aa = lane->base ;
  lane->base = 0 ;
  laneDestroy (lane) ;

  return aa ;
}

/***************************************************/
/***************************************************/

#ifdef JUNK
 /* realign the read onto its scf datafile */
static void importExtension (LOOK look, LANE *lane)
{ int n, top, end, sens, i, j, a1, a2, x1, x2 ;
  char *cp ;
  Array 
    dna = lane->dna, abi = 0, base2 = 0 , basePos2 = 0 ,
    errArray = 0  ;
  A_ERR *ep ;
  
  return ;

  top = lane->clipTop ;
  end = lane->clipEnd ;
  
  if (top)  /* clipping already known */
    return ;
  
        /* make a copy and realign it */
  abi = arrayCopy (lane->base) ;
  n = arrayMax(abi) ;
  cp = arrp (abi, 0, char) - 1 ;
  while (cp++, n--) *cp &= 0x0f ; /* unflag */
  
  if (!dnaAlignForceMatch (abi, 0, arrayMax(lane->base) -1,
			   lane->dna, 0, arrayMax (lane->dna) - 1, 
			   &top, &end, &sens) || sens == -1)
    { invokeDebugger () ;
      arrayDestroy (abi) ;
      return ;
    }

      /* import the clipped regions from the base call */       
  arrayDestroy (abi) ;
  if (top > 0 || end < arrayMax (dna)) /* then i should import */
    { abi = arrayCreate (arrayMax(lane->base) + 50, char) ;
      for (j = 0, n = 0 ; n < top ; n++)
	arr (abi, j++, char) = arr (lane->base, n, char) ;
      for (i = 0 ; i < arrayMax(lane->dna) ; i++)
	array (abi, j++, char) = arr (lane->dna, n, char) ;
      for (n = end + 1 ; n < arrayMax(lane->base) ; n++)
	array (abi, j++, char) = arr (lane->base, n, char) ;
      
      arrayDestroy(lane->dna) ;
      lane->dna = abi ; abi = 0 ;/* reset to avoid destroy(abi) */
    }
  
  lane->clipTop = top ;
  lane->clipEnd = end ;
  lane->clipExtend = end + 20 ;
  if (lane->clipExtend > arrayMax(lane->dna))
    lane->clipExtend = arrayMax(lane->dna) ;
  /* include the edits inside the base call itself */

  a1 = x1 = 0 ; 
  a2 = arrayMax(lane->base) - 1 ;
  x2 = arrayMax(lane->dna) - 1 ;
  errArray = 
    dnaAlignCptErreur (lane->base, lane->dna, 
		       &a1, &a2, &x1, &x2) ;

  ep = arrp (errArray, 0, A_ERR) ;
  n = arrayMax(errArray) ;
  base2 = arrayCreate (arrayMax(lane->base) + 50, char) ;
  basePos2 = arrayCreate (arrayMax(lane->base) + 50, short) ;

  arrayDestroy (errArray) ;
  arrayDestroy (abi) ;
  arrayDestroy (base2) ;
  arrayDestroy (basePos2) ;
  
   /* forget in related displays */
  while (i--)
    array (dna, j++, char) = *cp++ ;
  fMapTraceForget(lane->key) ;
  dnaAlignForget (lane->key) ;
}
#endif
/************************************************************/

BOOL traceGetLane (LOOK look, LANE *lane)
{ char *cp ;
  int z1, z2, x1 = lane->x1, x2 = lane->x2 ;
  Array laneDna = 0 ; 
  KEY dnaKey, key = lane->key,
    baseCallKey = 0, baseQualityKey = 0, basePositionKey = 0 ;
  OBJ obj = 0 ;
#ifndef NON_GRAPHIC
  KEY myCol ;
  extern int fMapQueryColor (KEY key) ;
#endif

  if (lane->scf) return TRUE ;
  lane->scf = 1 ; /* prevent recursion */
  lane->key = key ;
  obj = bsCreate (key) ;
  if(!obj)
    return FALSE ;
    
  cp = name(key) + strlen (name(key)) ;
  while (cp > name(key) && *cp != '.') cp-- ;
  if (bsGetKey (obj, _DNA, &dnaKey))
    laneDna = dnaGet (dnaKey) ;
  lane->color = 0 ;
#ifndef NON_GRAPHIC
  if (bsGetKeyTags (obj, _Colour, &myCol))
    lane->color = myCol - _WHITE + WHITE ;
  if ((myCol = fMapQueryColor (key)))
    lane->color = myCol ;
#endif
  if (bsGetKey (obj, _BaseCall, &baseCallKey))
    lane->baseCallKey = baseCallKey ;
  if (bsGetKey (obj, _SCF_Position, &basePositionKey))
    lane->basePositionKey = basePositionKey ;
  if (bsGetKey (obj, _Quality, &baseQualityKey))
    lane->baseQualityKey = baseQualityKey ;
  if (bsGetData (obj, _Hand_Clipping, _Int, &z1) &&
      bsGetData (obj, _bsRight, _Int, &z2))
    { lane->handClipTop = z1 - 1 ; lane->handClipEnd = z2 ; }
  else
    lane->handClipTop = lane->handClipEnd = 0 ;
  if (!lane->clipEnd)
    { 
      if (bsGetData (obj, _Clipping, _Int, &z1) &&
	  z1 >= 0 &&
	  bsGetData (obj, _bsRight, _Int, &z2))
	{ 
	  lane->clipTop = z1 - 1 ; /* top included */
	  lane->clipEnd = z2 ;   /* end excluded */
	  if (laneDna && lane->clipEnd > arrayMax (laneDna))
	    z2 = lane->clipEnd = arrayMax (laneDna) ;
	  lane->clipExtend = z2 + 50 ;
	  if (laneDna && lane->clipExtend > arrayMax (laneDna))
	    lane->clipExtend = arrayMax (laneDna) ;
	}
      else
	{
	  lane->clipTop = 0 ;
	  lane->clipEnd = lane->clipExtend = 
	    laneDna ? arrayMax (laneDna) : 0 ;
	}
    }
  if (bsGetData (obj, _Vector_Clipping, _Int, &z1) &&
      z1 >= 0 &&
      bsGetData (obj, _bsRight, _Int, &z2))
    { 
      lane->vectorTop = z1 - 1 ; /* top included */
      lane->vectorEnd = z2 ;   /* end excluded */
      if (laneDna && lane->vectorEnd > arrayMax (laneDna))
	z2 = lane->vectorEnd = arrayMax (laneDna) ;
    }
  else
    lane->vectorTop = lane->vectorEnd = -1 ;
  
  if (!bsFindTag (obj, _Assembly_tags))
    lane->hasTag = 1 ;  /* evaluated to absent */
  bsDestroy(obj) ;
  
  lane->dnaKey = dnaKey ;
  lane->hide = look ? look->hide : 0 ;
  
  if (x1 < x2)
    { lane->upSequence = FALSE ;
      lane->x3 = x2 + (lane->clipExtend - lane->clipEnd) ;
    }
  else
    { lane->upSequence = TRUE ;
      lane->x3 = x2  - (lane->clipExtend - lane->clipEnd) ;
    }

  lane->scf = 2 ;  /* success */
  if (lane->dna) messcrash ("anomaly in traceGetlane") ;
  lane->dna = laneDna ;
  lane->baseCall = 0 ;
  
  if (look)
    laneMakeErrArray(look, lane) ;
  lane->laneMax = 0 ; /* to force a call to laneHeight */

  return TRUE ;
}

/***************************************************/
/***************************************************/
static KEY _New_Read = 0 ;
  /* Extract DNA from SCF files */
static BOOL newScf2dna (KEY seq)
{ OBJ obj ;
  KEY dnaKey = 0 ;
  int dx = 0, clipTop = 1, clipEnd = 1;
  LANE *lane = (LANE*) messalloc (sizeof (struct LaneStruct)) ;
  
  lane->key = seq ;
  if (!baseCallGet (lane) ||
      lane->scf < 4)
    { laneDestroy (lane) ;
      messfree (lane) ;
      return FALSE ;
    }

  seq2ace(lane) ;                            /* cree lane->dna */
  baseCallStore (lane) ;  /*  register basecall */
  lexaddkey (name(seq), &dnaKey, _VDNA) ;     /* register it */
  dx = arrayMax (lane->dna) ;

  if ((obj = bsUpdate (lane->key)))
    { bsAddKey (obj, _DNA, dnaKey) ;
      bsAddData (obj, _bsRight, _Int, &dx) ;
      bsSave (obj) ;
    }

  lane->clipTop = 0 ;
  if (lane->dna) lane->clipEnd = arrayMax(lane->dna) ;
  else lane->clipEnd = 1 ;
  
  if (baseCallUnclipLane (lane, 'G', &dx) &&   /* evaluate quality */
      (obj = bsCreate (lane->key)) && bsGetData (obj, _Clipping, _Int, &clipTop))
    bsGetData (obj, _bsRight, _Int, &clipEnd) ;
  else
    { clipTop = 1 ; clipEnd = dx ;}
  bsDestroy (obj) ;

  clipTop-- ;
  if (!dnaAlignCheckSequence (lane->dna, clipTop, clipEnd,1)) /* rejected */
    if ((obj = bsUpdate (lane->key)))
      { clipTop++ ; clipEnd = clipTop + 1 ;
	bsAddData (obj, _Clipping, _Int, &clipTop) ;
	bsAddData (obj, _bsRight, _Int, &clipEnd) ;
	bsAddData (obj, _Vector_Clipping, _Int, &clipTop) ;
	bsAddData (obj, _bsRight, _Int, &clipEnd) ;
	bsSave (obj) ;
      }
  dnaStoreDestroy (dnaKey, lane->dna) ;  lane->dna = 0 ;
  laneDestroy (lane) ;
  messfree (lane) ;
  return TRUE ;
}

/***************************************************/
  /* eliminate  New_read */

static void clearNewRead (void) 
{ OBJ obj = 0 ;
  KEYSET ks = query (0, ">? New_Read") ; 
  int i = keySetMax (ks) ;

  while (i--)
    if ((obj = bsUpdate(keySet(ks, i))))
      { if (bsFindTag (obj, _New_Read))
	  bsRemove (obj) ;
	bsSave (obj) ;
      }
  keySetDestroy (ks) ;
}

KEYSET  baseCallNewScf(void)
{ KEYSET 
    ks = query (0, ">? New_Scf"), 
    ks1 = keySetCreate () ;
  int i = keySetMax (ks), j = 0 ;
  OBJ obj ; KEY *kp ;

  if (!_New_Read)
    lexaddkey ("New_Read", &_New_Read, _VSystem) ;
      
  sessionGainWriteAccess ()  ;      /*  Get if possible */

  for (i = 0 ; i < keySetMax(ks) ; i++)
    {
      if (!isWriteAccess())
	break ;
      printf("%s\n", name(keySet(ks, i))) ;
      if (newScf2dna (keySet(ks, i)))
	keySet (ks1, j++) = keySet (ks, i) ;
      if ((! (j%1000)))
	{ sessionClose(1) ;  sessionGainWriteAccess () ;}
    }
  keySetDestroy (ks) ;

  if (keySetMax(ks1))
    { clearNewRead () ;
      i = keySetMax (ks1) ;
      kp = arrp (ks1, 0, KEY) - 1 ;
      while (kp++, i--)
	if ((obj = bsUpdate (*kp)))
	  { bsAddTag (obj, _New_Read) ;
	    bsSave (obj) ;
	  }
    }
  return ks1 ;
}

/***************************************************/
/***************************************************/

BOOL baseCorrel(Read* seq1, int x1, BOOL direct,
	       Read* seq2, int x2, int ll, int nstep, int step, int *bestdxp)
{ 
  int
    x11 = x1 - ll/2, x21 = x2 - ll/2 ,
    x12 = x1 + ll/2, x22 = x2 + ll/2 ,
    dx, dxmax = nstep * step,
    i,
    max1 = seqMax(seq1) ,
    max2 = seqMax(seq2) , n, di = 1 ;
  double  s11, s22, s12 ;
  float z, bestz ;
  TRACE *bp1[4], *bp2[4], *ip1[4] , *ip2[4] ;
  int t1, t2 ;
  
  if (x11 < dxmax || x21 < dxmax || 
      x12 + dxmax >= max1 || x22 + dxmax >= max2)
    return FALSE ;

  bp1[0] = seq1->traceA ;
  bp1[1] = seq1->traceT ;
  bp1[2] = seq1->traceG ;
  bp1[3] = seq1->traceC ;
  
  if (direct)
    { bp2[0] = seq2->traceA ;
      bp2[1] = seq2->traceT ;
      bp2[2] = seq2->traceG ;
      bp2[3] = seq2->traceC ;
    }
  else
    { bp2[1] = seq2->traceA ;
      bp2[0] = seq2->traceT ;
      bp2[3] = seq2->traceG ;
      bp2[2] = seq2->traceC ;
      x21 = x22 ; di = -1 ;
    }

      
  bestz = - 1 ;
  for (dx = -dxmax ; dx < dxmax ; dx += step)
    {
      i = 4 ; while (i--) ip1[i] = bp1[i] + x11 ;
      i = 4 ; while (i--) ip2[i] = bp2[i] + x21 + dx ;
  
      n = ll ;
      s12 = s11 = s22 = 0 ;
      while (n--)
	{ i = 4 ; 
	  while(i--) 
	    { t1 = *ip1[i]++ ;
	      t2 = *ip2[i] ; ip2[i] += di ;
              s11 += t1 * t1 ; s22 += t2 * t2 ; s12 += t1 * t2 ;
	    }
	}
      if (s12 > 0)
	{ z = s12*s12 ; /* /(s11*s22) ; */
	  if (z > bestz)
	    { *bestdxp = dx ; bestz = z ;
	    }
	}
    }
  return bestz > 0 ? TRUE : FALSE ;
}

/***************** energy ***********************************/
/*
static Array seqFindEnergy(Read* seq, int min, int max)
{
  TRACE *a, *t, *g, *c ; int i ; int e ; 
  Array energy ;

  if (min < 0 || max > seqMax(seq))
       return 0 ;
 
  energy = arrayCreate(seqMax(seq), int) ;
  
  a = seq->traceA + min ; 
  t = seq->traceT + min ;  
  g = seq->traceG + min ;  
  c = seq->traceC + min ;  
  
  array(energy,max-min, int) = 0 ; // make room 
  for (i= min; i< max ; i++, a++, t++, g++, c++)
    {
      e = (int)*a * *a  + (int)*t * *t + (int)*g * *g +(int)*c * *c ;
      arr(energy, i - min, int) = e ;

    }
  return energy ;
}
*/
/*************/
 /* carre de la derivee de l'enveloppe */
 /* corrigee pour le bruit de fond */
int seqEnergyOfDerivee (Read* seq, int min, int max)
{
  TRACE *a, *t, *g, *c, a0, t0, g0, c0, y, y2 ; 
  int x, i ;
  double e = 0, s1 = 0, s2 = 0, dd ;

  if (!seq || min < 0 || max > seqMax(seq))
       return 0 ;
 
  a = seq->traceA + min ; 
  t = seq->traceT + min ;  
  g = seq->traceG + min ;  
  c = seq->traceC + min ;  
  
  a0 = *a ; t0 = *t ; g0 = *g ; c0 = *c ;

  for (i= min; i< max ; i++, a++, t++, g++, c++)
    { x = 0 ; y = 0 ; y2 = 0 ;
      s1 += *a * *a + *t * *t + *g * *g  + *c * *c ;
      if (*a > y) { y = *a ; x = *a - a0 ; }
      if (*t > y) { y = *t ; x = *t - t0 ; }
      if (*g > y) { y = *g ; x = *g - g0 ; }
      if (*c > y) { y = *c ; x = *c - c0 ; }
      if (*a != y && *a > y2) { y2 = *a ; }
      if (*t != y && *t > y2) { y2 = *t ; }
      if (*g != y && *g > y2) { y2 = *g ; }
      if (*c != y && *c > y2) { y2 = *c ; }
      e += x * x ; if (y2) s2 += (y - y2) * (y - y2) ;
      a0 = *a ; t0 = *t ; g0 = *g ; c0 = *c ;
    }
  if (!s1) s1 = 1 ;
  dd = e * e / s1 ; 
  dd = dd * s2/s1 ;
  return (int) dd ;
}

/***************************************************/
/***************************************************/
/*********** extrema ************************************/

static void findMaxima(TRACE *bp, int max, int delta, Array xx)
{ short d = delta ; TRACE lasty = *bp, *bp1 ;
  int i, i1, j = 0, lastx = 0 ;
  BOOL up = FALSE ;
  
  for (i = 0 ; i < max ; i++, bp++) 
    {
      if (up)
	{ if (*bp < lasty - d)
	    { up = FALSE ;
	      bp1 = bp - i + lastx ; i1 = 0 ;
	      while (*++bp1 == lasty) i1++ ;
	      lastx += i1/2 ;  /* middle of a saturation */
	      array(xx, j, Extremum).x = lastx ;
	      arr(xx, j, Extremum).y = lasty ;
	      j++ ;
	      lastx = i ; lasty = *bp ;
	    }
	  if (*bp > lasty)
	    { 
	      lastx = i ; lasty = *bp ;
	    }
	}
      else
	{ if (*bp > lasty + d)
	    { up = TRUE ;
/*	don t store the minima
              array(xx, j, Extremum).x = lastx ;
	      arr(xx, j, Extremum).y = lasty ;
	      j++ ;
*/
	      lastx = i ; lasty = *bp ;
           }
	  if (*bp < lasty)
	    { 
	      lastx = i ; lasty = *bp ;
	    }
	}
    }
}

/************/

static BOOL find4Extrema(LANE *lane, Array aExtrema[])
{ TRACE *bp[4], yy, y[4] ;
  int x, i, i1, j, k, max, delta ;  
  Extremum *e1, *e2 ;

  delta = 3 ;
  if(!lane->seq)
    return FALSE ;

  
  max = seqMax(lane->seq) ;

  bp[0] = lane->seq->traceA ; 
  bp[1] = lane->seq->traceG ; 
  bp[2] = lane->seq->traceC ; 
  bp[3] = lane->seq->traceT ; 

  for (i=0 ; i<4 ; i++)
    { aExtrema[i] = arrayReCreate(aExtrema[i], 
				  100, Extremum) ;
      findMaxima(bp[i], max, delta, aExtrema[i]) ;
    }
  for (i=0 ; i<4 ; i++)
    { 
      for (j = 0 , k = 0, e1 = e2 = arrp(aExtrema[i], 0, Extremum) ;
	   j < arrayMax(aExtrema[i]) ; j++, e1++)
	{ x = e1->x ;
	  for (i1 = 0 ; i1 < 4 ; i1++)
	    y[i1] = *(bp[i1] + x)  ;
	  yy = y[i] ;
	  for (i1 = 0 ; i1 < 4 ; i1++)
	    if (yy < y[i1])
	      goto forget ;
	  if (e1 != e2)
	    *e2 = *e1 ; 
	  e2++ ; k++ ;
	forget:
	  continue ;
	}
      arrayMax(aExtrema[i]) = k ;
    }  
  return TRUE ;
}


static void baseCallAdd (LANE *lane, Array bc, int dx, int *debut, int fin)
{ int reject, 
    n, ni, i, ibc, i1, j1, x, x0 = 0, x1, x2, y, yy, dx1, dx2, ddx ;
  BASECALL *bb = 0, *bb1 = 0;
  TRACE *bp[4], *tp ;
  static int lastdebut, lastx2 ;

  bp[0] = lane->seq->traceA ; 
  bp[1] = lane->seq->traceG ; 
  bp[2] = lane->seq->traceC ; 
  bp[3] = lane->seq->traceT ; 

  if (*debut > lastdebut)
    { ibc = *debut ;
      while (ibc > 0 && arrp(bc, ibc, BASECALL)->x > lastx2)
	ibc-- ;
      *debut = ibc > 0 ? ibc : 0 ;
    }
  lastdebut = *debut ;
  for (ibc = *debut ; 
       ibc < fin && ibc < arrayMax(bc) - 1 ; ibc++)  /* will stop sur avant dernier */
    { bb = arrp (bc, ibc, BASECALL) ;
      x1 = bb->x ; x2 = (bb + 1)->x ; 
      ni = (x2 - x1 + .3  * dx) / dx ; /* nb of steps in interval */
      if (ni < 2)
	continue ;

      yy = 0 ; j1 = 0 ;

      reject = 7 ;
    tryotherpeak:
      x = x1 + (x2 - x1) / ni ;
      for (i = x - .35 * dx ; i <= x + .35 * dx ; i++)
	{ i1 = 4 ;
	  while (i1--) /* find highest trace */
	    { if (i1 == reject) continue ;
	      y = *(bp[i1] + i) ;
	      if (yy < y) { yy = y ; j1 = i1 ; x0 = i ; }
	    }
	}
      i1 = ibc - 16 ; if (i1 < 0) i1 = 0 ;
      for (y = yy, n = 1, i = i1 ; i < arrayMax(bc) && n < 5 ; i++)
	{ if (arrp(bc, i, BASECALL)->t == j1)
	    { y += arrp(bc, i, BASECALL)->y ; n++ ; }
	}
      y /= n ; /* hauteur moyenne des pics  */
      if (yy * 4 <  y) /* not high enough to keep it */
	{ if (reject == 7)
	    { reject = j1 ;
	      goto tryotherpeak ; /* may be second peak is good */
	    }
	  if (x2 - x1 < 2.11 * dx)
	    continue ; /* no real peak here */
	  j1 = reject ; /* well, there is nothing better and i want one */
	}
      if (!yy || yy * 8 <  y) /* really not high enough to keep it */
	continue ;
      /* check that curvature is negative */
      
      i = x - .33 * dx ; /* centre on middle not on highest point */
      tp = bp[j1] + i ; x = i - 1 ;
      dx1 = *tp - *(tp - 1) ;
      i = .66  * dx ;
      ddx = 0 ;
      while (tp++, x++, i--)
	{ dx2 = dx1 ;
	  dx1 = *tp - *(tp - 1) ;
	  if (dx1 -  dx2 < ddx)
		{ ddx = dx1 - dx2 ; x0 = x ; }
	}
      
      if (ddx < 0 || ni > 2  || (x2 -  x1 > 2.1 * dx))
	{
	  /* insertion */
	  i = arrayMax(bc) ;
	  if (ibc > i - 2) messcrash("ibc too big in addBase") ;
	  bb1 = arrayp (bc, i, BASECALL) ; /* increases arrayMax */
	      while (--i > ibc) { *bb1 = *(bb1 - 1) ; bb1-- ;}
	  bb = arrp (bc, ibc, BASECALL) ;
	  bb1 = bb + 1 ;
	  bb1->x = x0 ; bb1->y = yy ; 
	  bb1->t = j1 ; bb1->flag = BC_LOW ;
	  if (bb->x >= bb1->x) invokeDebugger() ;
	}
    }
  lastx2 = bb->x - 10 ;
}

static void baseCallSuppress (Array bc, int dx, int debut, int fin)
{ int n, i, ibc, i1, imax, t, j1, x1, x2, y, yy, ymin ;
  BASECALL *bb, *bb0, *bb1, *bb2, *bbmin, *bbmax ;

  for (ibc = debut + 1 ; ibc < fin && ibc < arrayMax(bc) - 1 ; ibc++)  /* will stop sur avant dernier */
    { bb = bb0 = arrp (bc, ibc, BASECALL) ;
      bb1 = bb2 = bb ;
      t = bb->t ; imax = arrayMax(bc) ;
      i = ibc ; while (i-- > 0 && bb1->t == t) bb1-- ; 
      i = ibc ; while (++i < imax && bb2->t == t) bb2++ ;
      
      x1 = bb1 ->x ; x2 = bb2->x ; n = bb2 - bb1 - 1 ;
      if (x2 - x1 < (n + .6) * dx)
	{ yy = bb0->y ; j1 = bb0->t ;
	  i1 = ibc - 16 ; if (i1 < 0) i1 = 0 ;
	  for (y = yy, n = 1, i = i1 ; i < arrayMax(bc) && n < 5 ; i++)
	    { if (arrp(bc, i, BASECALL)->t == j1)
		{ y += arrp(bc, i, BASECALL)->y ; n++ ; }
	    }
	  y /= n ; /* hauteur moyenne des pics  */
	  bbmin = bb0 ; ymin = bbmin->y ;
	  for (bb = bb1 + 1 ; bb < bb2 ; bb++)
	    if (bb->y < ymin)
	      { bbmin = bb ;  ymin = bbmin->y ; }

	  if (ymin * 2 < y) /* low enough to remove */
	    { /* deletion */ 
              i = arrayMax(bc) - 1 ;
	      bb1 = bbmin ; bbmax = arrp (bc, i, BASECALL) ;
	      if (bb1 <  arrp (bc, 0, BASECALL))
		messcrash("nnanna nerr") ;
	      while (bb1 < bbmax) { *bb1 = *(bb1 + 1) ; bb1++ ; }
	      arrayMax(bc)-- ;
	    }
	  else
	    bb0->flag |= BC_LOW ;
	}
    }
}

static int getDx (KEY key, Array bc, int dx0, int deb, int fin)
{ static Array hh = 0 ;
#ifndef NON_GRAPHIC
  static int tour = 0 ;
#endif

  int gg[] = 
    { 2, 6, 13, 27, 48, 72, 92, 100, 92, 72, 48, 27, 13, 6, 2 } ;
  int i, i1, icc, j, z, *ggp, *ih, dx, a, x1, x2, x3, y ;
  BASECALL *bb, *bb1 ;

  hh = arrayReCreate (hh, 200, int) ;
  if (fin > arrayMax(bc))
    fin = arrayMax(bc) ;
  if (fin < deb + 80) return dx0 ;

  for ( i = deb, bb = arrp(bc, i, BASECALL) ;
       i < fin - 4 ; bb++, i++)
    for (bb1 = bb + 1, j = 1 ; j < 4 ; bb1++, j++)
      { z = bb1->x - bb->x ;
	icc = 15 ; i1 = z + 7 ; if (i1 < 0) i1 = 0 ;
	ggp = &gg[0] ; ih = arrayp(hh, i1, int) ; 
	while (icc-- && i1--)
	    *ih-- += *ggp++ ; 
      }

  a = dx0/3 ; x1 = x2 = x3 = 0 ; 
  array(hh, 2*dx0, int) = 0 ; /* zeiterbiftleiter */
  y = 0 ;
  for (i = dx0 - a, ih = arrayp(hh, i, int) ; i < dx0 + a ; ih++, i++)
    if (*ih > y) { y = *ih ; x1 = i ; }
  y = 0 ;
  for (i = 2*x1 - a, ih = arrayp(hh, i, int) ; i < 2*x1 + a ; ih++, i++)
    if (*ih > y) { y = *ih ; x2 = i ; }
  y = 0 ;
  for (i = x1 + x2 - a, ih = arrayp(hh, i, int) ; i < x1 + x2 + a ; ih++, i++)
    if (*ih > y) { y = *ih ; x3 = i ; }
  
  dx = (9*x1 + 3*x2 + x3) / 18 ;

#ifndef NON_GRAPHIC
  if (tour++ < 0)
    plotHisto(messprintf("%s: zone %d %d",name(key), deb, fin), 
	      arrayCopy(hh)) ;
#endif

  return dx ;
}

BOOL findBaseCall (LANE *lane)
{ int x, mx[4], i, jbc, i1, dx, x1, x2 ;
  Extremum *ep[4] ;
  Array bc, aExtrema[4] ; 
  BASECALL *bb ;
  TRACE *bp[4] ; 

  if (arrayExists(lane->baseCall))
    return TRUE ;
  if (lane->scf < 3 && !baseCallGetSeq (lane))
    return FALSE ;
  bc = lane->baseCall = arrayCreate (1000, BASECALL) ;
  jbc = 0 ;

  bp[0] = lane->seq->traceA ; 
  bp[1] = lane->seq->traceG ; 
  bp[2] = lane->seq->traceC ; 
  bp[3] = lane->seq->traceT ; 
  i = 4 ;
  while (i--)
    aExtrema[i] = arrayCreate (1000, Extremum) ;

  if (!find4Extrema(lane, aExtrema))
    return FALSE ;

  i = 4 ;
  while (i--)
    { mx[i] = arrayMax (aExtrema[i]) ;
      ep[i] = arrp (aExtrema[i], 0, Extremum) ;
    }
 encore:
  x = 1000000 ; i1 = - 1 ;
  i = 4 ;
  while (i--)
    { if (x > ep[i]->x && mx[i] > 0)
	{ x = ep[i]->x ; i1 = i ; }
    }
  bc = lane->baseCall ; /* this seems to be needed, 5/2000 */
  if (i1 >= 0)
    { bb = arrayp(bc, jbc++, BASECALL) ;
      bb->t = i1 ; bb->x = ep[i1]->x ; bb->y = ep[i1]->y ; 
      if (bb->y < 16)
	bb->flag |= BC_LOW ;
      else 
	{ 
	  i = 4 ;
	  while (i--)
	    if (i != i1 && 
		( (*(bp[i] + bb->x) * 4) >  2 * bb->y))
	      { bb->flag |= BC_LOW ; ; break ; }
	}
      
      ep[i1]++ ; mx[i1]-- ;
      if (jbc > 1 && bb->x < (bb - 1)->x)
      invokeDebugger() ;
      goto encore ;
    }
  i = 4 ;
  while (i--) 
    arrayDestroy (aExtrema[i]) ;

  i = arrayMax(bc) > 300 ? 300 : arrayMax(bc) - 1 ;
  if (!i)
    return FALSE ;
  dx = arr(bc, i, BASECALL).x / i ;

  for (x1 = 0 ; x1 < arrayMax(bc) ; x1 += 95) 
    { dx = getDx (lane->key, bc, dx, x1, x1 + 100) ;
      x2 = x1 + 100 ;
      if (dx < 8) break ;
      baseCallAdd (lane, bc, dx, &x1, x2) ;
      baseCallSuppress (bc, dx, x1, x2) ;
    }
  return TRUE ;
}

/***************************************************/
/***************************************************/

static void cleanErr (Array err, Array dnaLong, Array dnaShort, BOOL isUp)
{ int i, ii, j, n = arrayMax(err), nn, start, stop ;
  A_ERR *eq, *ep, *epMax ;
  char *cl, *cs, cc ;
  int sens = isUp ? -1 : 1 ;
  ALIGN_ERROR_TYPE type ;

  if (n < 2) return ;
  if (isUp)
    { ep = arrp (err, 0, A_ERR) - 1 ;
      n-- ;
      while (ep++, n--)
	if ((ep->type == INSERTION) && ((ep+1)->type == ERREUR) &&
	    (ep->iShort == (ep+1)->iShort + 1))
	  { type = ep->type ; ep->type = (ep+1)->type ; (ep+1)->type = type ;
	    i = ep->iLong ; ep->iLong = (ep+1)->iLong ; (ep+1)->iLong = i + 1;
	  }
    }

  n = arrayMax (err) ;
  ep = arrp(err, 0, A_ERR) - 1 ; ii = -1 ;
  while (ep++, ii++, n--)
    { 
      switch (ep->type)
	{
	case AMBIGUE:
	  nn = arrayMax (err) - 1 ;
	  if (isUp)
	    { 
	      if (ii > 0 &&
		  (ep - 1)->type == INSERTION &&
		  (ep - 1)->iShort < arrayMax(dnaShort) &&
		  (ep - 1)->iLong >= 0 &&
		  (ep - 1)->iLong < arrayMax(dnaLong))
		{ cl = arrp(dnaLong, (ep - 1)->iLong, char) ; 
		  cs = arrp(dnaShort, (ep - 1)->iShort, char) ;
		  i = j = 0 ;
		  cc = *cl ; while (*cl++ == cc) i++ ;
		  if (sens == -1)
		    { cc = complementBase(cc) ;
		      while (*cs-- == cc) j++ ;
		    }
		  else
		    while (*cs++ == cc) j++ ; 
		  if (i >= j && (ep - 1)->iShort == ep->iShort + j)
		    { 
		      eq = ep - 2 ;
		      i = n + 1 ;
		      while (eq++, i-- > 0) *eq = *(eq + 1) ;
		      arrayMax(err) -= 1 ;
		      ep-- ; ii-- ;
		      ep->type = INSERTION ;
		      ep->iLong++ ;
		    }
		}
	      if (ii < nn  &&
		  (ep + 1)->type == INSERTION &&
		  (ep + 1)->iShort < arrayMax(dnaShort) &&
		  (ep + 1)->iLong >= 0 &&
		  (ep + 1)->iLong < arrayMax(dnaLong))
		{ cl = arrp(dnaLong, (ep + 1)->iLong, char) ; 
		  cs = arrp(dnaShort, (ep + 1)->iShort, char) ;
		  i  = j = 0 ;
		  cc = *cl ; while (*cl++ == cc) i++ ;
		  if (sens == -1)
		    { cc = complementBase(cc) ;
		      while (*cs-- == cc) j++ ;
		    }
		  else
		    while (*cs++ == cc) j++ ; 
		  if (i >= j && (ep + 1)->iShort == ep->iShort + j)
		    { 
		      eq = ep - 1 ;
		      i = n  ;
		      while (eq++, i-- > 0) *eq = *(eq + 1) ;
		      arrayMax(err) -= 1 ;
		      ep-- ; ii-- ;
		      ep->type = INSERTION ;
		      ep->iLong++ ;
		    }
		}
	    }
	  else
	    { if ( ii < nn  &&
		  (ep + 1)->type == INSERTION &&
		  ep->iShort < arrayMax(dnaShort) &&
		  ep->iLong >= 0 &&
		  ep->iLong < arrayMax(dnaLong))
/*		  (ep + 1)->type == INSERTION_DOUBLE) plus difficile */
		{ cl = arrp(dnaLong, ep->iLong, char) ; 
		  cs = arrp(dnaShort, ep->iShort, char) + sens ;
		  i = j = 0 ;
		  cc = *cl ; while (*cl++ == cc) i++ ;
		  if (sens == -1)
		    { cc = complementBase(cc) ;
		      while (*cs-- == cc) j++ ;
		    }
		  else
		    while (*cs++ == cc) j++ ; 
		  if (i == j && (ep + 1)->iShort == ep->iShort + i * sens)
		    { ep->type = INSERTION ;
		      eq = ep ;
		      i = n - 1 ;
		      while (eq++, i-- > 0) *eq = *(eq + 1) ;
		      arrayMax(err) -= 1 ;
		      n-- ;
		    }
		}
	    }
	  break ;
	case INSERTION: case INSERTION_DOUBLE:
	  start = ep->iShort ;
	  epMax = arrp (err, arrayMax (err) - 1, A_ERR) + 1 ;
	  if (n >= 0 && !isUp && ep->iShort >= 0 &&
	      ep->iShort < arrayMax(dnaShort))
	    { cs = arrp(dnaShort, ep->iShort, char) ;
	      cc = *cs ; i = arrayMax(dnaShort) - ep->iShort ;
	      while (i-- && *(++cs) == cc)
		{ ep->iShort++ ;
		  ep->iLong++ ;
		}
	      if (ep->type == INSERTION_DOUBLE)
		ep->iLong-- ;
	    }
	  else if (n >= 0 && isUp && ep->iShort >= 0 &&
		   ep->iShort < arrayMax(dnaShort))
	    { cs = arrp(dnaShort, ep->iShort, char) ;
	      cc = *cs ; i = ep->iShort ;
	      while (i--  > 0 && *(--cs) == cc && ep->iShort > 0)
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
	  /* condense insertion-error, nov 1 , 2003 */
	  nn = arrayMax (err) - 1 ;
	  if (ii < nn &&
	      ep->type == INSERTION && n &&
	      (ep+1)->type == ERREUR &&
	      ep->iLong == (ep+1)->iLong &&
	      ep->iShort - 1 == (ep+1)->iShort)
	    { 
	      nn = n - 1 ;
	      eq = ep ;
	      while (++eq, nn-- > 0)
		*eq = *(eq + 1) ;
	      n-- ;
	      arrayMax (err) -= 1 ;
	    }
	  if (ep->type == INSERTION_DOUBLE && ii < nn)
	    { nn = n + 1 ;
	      eq = arrayp (err, arrayMax (err), A_ERR) + 1 ;
	      ep = arrp (err, ii,A_ERR) ; /* may have been reallocated */
	      while (eq--, nn-- > 0)
		*eq = *(eq - 1) ;
	      ep->type = INSERTION ;
	      ep->baseShort = W_ ;
	      if (!isUp)
		ep->iShort-- ;
	      ep++ ; ii++ ;
	      ep->type = INSERTION ;
	      if (isUp)
		ep->iShort++ ;
	      ep->baseShort = W_ ;
	      ep++ ; ii++ ; if (n) n-- ;
	    }
	  break ;
#ifdef JUNK
	case INSERTION: 
	  if (ii < nn &&
	      (ep+1)->type == TROU)
	    {
	      cc = arr (dnaShort, ep->iShort +1, char) ;
	      for (i = ep->iShort + 1 ; i < (ep+1)->iShort ; i++)
		if (cc != arr(dnaShort, i, char))
		  goto idouble ;
	       for (i = ep->iLong ; i <= (ep+1)->iLong ; i++)
		 if (cc != arr(dnaLong, i, char))
		   goto idouble ;
	       ep->type = ERREUR ;
	       ep->baseShort = cc ;
	       (ep+1)->type = 0 ;  (ep+1)->baseShort = cc ;
	       break ;
	    }
	idouble:
	case INSERTION_DOUBLE:
	  start = ep->iShort ;
	  epMax = arrp (err, arrayMax (err) - 1, A_ERR) + 1 ;
	  if (n >= 0 && !isUp && ep->iShort >= 0 &&
	      ep->iShort < arrayMax(dnaShort))
	    { cs = arrp(dnaShort, ep->iShort, char) ;
	      cc = *cs ; i = arrayMax(dnaShort) - ep->iShort ;
	      while (i--  > 0 && *(++cs) == cc)
		{ ep->iShort++ ;
		  ep->iLong++ ;
		}
	      if (ep->type == INSERTION_DOUBLE)
		ep->iLong-- ;
	    }
	  else if (n >= 0 && isUp && ep->iShort >= 0 &&
		   ep->iShort < arrayMax(dnaShort))
	    { cs = arrp(dnaShort, ep->iShort, char) ;
	      cc = *cs ; i = ep->iShort ;
	      while (i--  > 0 && *(--cs) == cc && ep->iShort > 0)
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
	  nn = arrayMax (err) - 1 ;
	  if (ep->type == INSERTION_DOUBLE && ii < nn)
	    { nn = n + 1 ;
	      eq = arrayp (err, arrayMax (err), A_ERR) + 1 ;
	      ep = arrp (err, ii,A_ERR) ; /* may have been reallocated */
	      while (eq--, nn-- > 0)
		*eq = *(eq - 1) ;
	      ep->type = INSERTION ;
	      ep->baseShort = W_ ;
	      if (!isUp)
		ep->iShort-- ;
	      ep++ ; ii++ ;
	      ep->type = INSERTION ;
	      if (isUp)
		ep->iShort++ ;
	      ep->baseShort = W_ ;
	      ep++ ; ii++ ; if (n) n-- ;
	    }
	  break ;
#endif
	case TROU: case TROU_DOUBLE:
	  if (n >= 0 && ep->iLong >= 0 && ep->iLong < arrayMax(dnaLong))
	    { cs = arrp(dnaLong, ep->iLong, char) ;
	      cc = *cs ; i = arrayMax(dnaLong) - ep->iLong ;
	      while (i--  > 0 && *(++cs) == cc)
		{ ep->iShort += sens ;
		  ep->iLong++ ;
		}
	      if (ep->type == TROU_DOUBLE && !isUp)
		ep->iShort -- ;
	    }
	  break ;
	case ERREUR:
	  break ;
	default:
	  break ;
	}
    }

}

/*********************************************************/
static void fuseErr (Array err1, Array err2)
{ int i, n1 = arrayMax(err1), n2 = arrayMax(err2) ;
  for (i = 0 ; i < n2 ; i++)
    array(err1, n1 + i, A_ERR) = arr(err2, i, A_ERR) ;
}

static Array makeErr2 (Array longDna, Array shortDna, 
		       int *x1p, int *x2p, int *c1p, int *c2p, int new)
{ Array err = 0, err2 = 0 ;
  int nn = 0, c1, c2, d2, e2, x1, x2, y2, z2 ;

  c1 = *c1p ; x1 = *x1p ;
  switch (new)
    {
    case 5:  /* favor insert */
      err = trackErrors (shortDna, c1, c2p, longDna, x1, x2p, &nn, err, 2) ; 
      break ;
    case 4:  /* favor delete */
      err = trackErrors (shortDna, c1, c2p, longDna, x1, x2p, &nn, err, 1) ; 
      break ;
    case 3:  /* normal */
      err = trackErrors (shortDna, c1, c2p, longDna, x1, x2p, &nn, err, 0) ; 
      break ;
    case 2:  /* try all */
      c2 = d2 = e2 = *c2p ; x2 = y2 = z2 = *x2p ;  /* do notpass err to trackErrors, we need 2 simultaneous err */
      err = trackErrors (shortDna, c1, &c2, longDna, x1, &x2, &nn, 0, 0) ; 
      if (4 * arrayMax(err) <  (*x2p - x1))
	{ *c2p = c2 ; *x2p = x2 ;  break ; }
      err2 = trackErrors (shortDna, c1, &d2, longDna, x1, &y2, &nn, 0, 1) ; 
      if (arrayMax(err2) < arrayMax(err))
	{  *c2p = d2 ; *x2p = y2 ; arrayDestroy (err) ; err = err2 ; err2 = 0 ; *c2p = c2 ; *x2p = x2 ; break ; }
      err2 = trackErrors (shortDna, c1, &e2, longDna, x1, &z2, &nn, 0, 2) ; 
      if (arrayMax(err2) < arrayMax(err))
	{  *c2p = e2 ; *x2p = z2 ; arrayDestroy (err) ; err = err2 ; err2 = 0 ; *c2p = c2 ; *x2p = x2 ; break ; } 
      *c2p = c2 ; *x2p = x2 ;
      break ;
    default:
      err = dnaAlignCptErreur (longDna, shortDna, x1p, x2p, c1p, c2p) ;
      break ;
    }

  return err ;
}

/*********************************************************/

Array baseCallCptErreur (Array dnaD, Array dnaR, Array dnaCourt, BOOL isUp,
			 int x1, int x2, int x3, 
			 int *clipTopp, int *cEndp, int *cExp, int mode)
{ int amax, ctop, c2, c3, i ;
  Array err1, err2, dnaLong ;
  A_ERR *ep, *eq, ee ;

  amax = arrayMax (dnaD) ;

  if (isUp)
    { x1 = amax - x1 - 1 ;
      x2 = amax - x2 - 1 ;
      x3 = amax - x3 - 1 ;
      dnaLong = dnaR ;
    }
  else
    dnaLong = dnaD ;

  ctop = *clipTopp ;
  x2-- ;  c2 = *cEndp -1 ;
  err1 = makeErr2 (dnaLong, dnaCourt, &x1, &x2, &ctop, &c2, mode) ;
  *cEndp = c2 = c2 + 1 ; x2++ ;  /* may have moved a bit */      
  *clipTopp = ctop ;
  c3 = *cExp - 1 ;  
  if (x2 < x3 && c2 < c3)
    { err2 = makeErr2 (dnaLong, dnaCourt, &x2, &x3, &c2, &c3, mode) ;
      x3++ ;
      *cExp = c3 + 1 ;
      fuseErr (err1, err2) ;  /* accumulates in err1 */
      arrayDestroy (err2) ;
    }
  if (isUp)
    { x1 = amax - x1 - 1 ;
      x2 = amax - x2 - 1 ;
      x3 = amax - x3 - 1 ;
      i = arrayMax(err1) ;
      if (i)
	{ ep = arrp (err1, 0, A_ERR) - 1 ;
	  while (ep++, i--)
	    { 
	      switch (ep->type)
		{
		case INSERTION:
		case INSERTION_DOUBLE:
		  ep->iLong-- ;
		  break ;
		case TROU:
		case TROU_DOUBLE:
		  ep->iShort-- ;
		  break ;
		case AMBIGUE:
		case ERREUR:
		  break ;
		default:
		  break ;
		}
	      ep->iLong = amax - ep->iLong - 1 ;
	    }
	  i = arrayMax(err1) ;
	  ep = arrp (err1, 0, A_ERR) ;
	  eq = arrp (err1, i - 1, A_ERR) ;
	  while (ep < eq)
	    { ee = *ep ; *ep++ = *eq ; *eq-- = ee ; }
	}
    }
  cleanErr(err1, dnaD, dnaCourt, isUp) ;
  return err1 ;
}

/*********************************************************/

void  laneMakeErrArray (LOOK look, LANE *lane)
{
  BOOL isUp = lane->x1 < lane->x2 ? FALSE : TRUE ;

  arrayDestroy (lane->errArray) ;
  if (!lane->dna)
    goto abort ;

  lane->errArray = baseCallCptErreur
    (look->dna, look->dnaR, lane->dna, isUp,
     lane->x1, lane->x2, lane->x3,
     &lane->clipTop, &lane->clipEnd, &lane->clipExtend, 
     2 + lane->favorDelete) ;

abort:
    if (!lane->errArray)
      lane->errArray = arrayCreate (12, A_ERR) ;
  return ;
}

/*******************************************************************************/

static BOOL patchError (Array dna, LANE *lane, A_ERR *ep, Array err)
{ int 
    t, i, nb, ns, jmax,
    c0, cmax,  /* in dna */
    s0, s1, s2, /* in lane base */
    j0, /* in basecall */
    x0, x1, x2 ;
  Array 
    bc = lane->baseCall, 
    base = lane->base, 
    basePos = lane->basePos ;
  char cc , *cp ;
  short *posp ;
  BASECALL *bb, *bb1, *bb2 ;
  A_ERR *ep1 ;

  switch (ep->type)
    {
    case TROU_DOUBLE:  /* attention, used later as a flag to do nothing */
    case INSERTION_DOUBLE: 
      return FALSE ;
    default:
      break ;
    }

  c0 = ep->iLong ; 
  s0 = s1 = s2 = ep->iShort ;
  s1-- ; s2++ ;
  if (s1 < 5 || s2 > arrayMax(basePos) - 10 ||
      c0 < 0 || c0 >= arrayMax (dna))
    return FALSE ;
  x1 = arr(basePos, s1, short) ;
  x0 = arr(basePos, s0, short) ;
  x2 = arr(basePos, s2, short) ;
  
    /* search for the last base left of x0 */
  j0 = 0 ; jmax = arrayMax (bc) ;
  bb = arrp (bc, j0, BASECALL) ;
  while (j0 < jmax && bb->x < x0)
    { j0 += 20 ; bb += 20 ; }
  if (j0 > 0)
    { bb -= 20 ; j0 -= 20 ; }
  while (j0 < jmax && bb->x < x0)
    { j0 ++ ; bb++ ;}
  if (j0 > 0 && bb->x > x0 )
    { bb -- ; j0-- ;}

        /* identify the base in bb notation */
  cmax = arrayMax(dna) ;
  if (ep->type == INSERTION)
    cc = ep->baseShort ;
  else
    cc =  arr(dna, c0, char) & 0x0f ;
  switch (cc)
    { 
    case A_: t = 0 ; break ;
    case G_: t = 1 ; break ;
    case C_: t = 2 ; break ;
    case T_: t = 3 ; break ;
    default: return FALSE ;
    }

        /* locate preciselly the correct base */
 /* find closest base */
  if (bb->x - x0 > x0 - (bb - 1)->x) bb-- ;
  else if (x0 - bb->x  > (bb + 1)->x - x0) bb++ ;
  switch (ep->type)
    {
    case AMBIGUE: /* find closest base is sufficient */
      break ;
    case ERREUR:
    case INSERTION:
      if (bb->t != t && (bb+1)->t == t && bb->x < x0) bb++ ;
      else if (bb->t != t && (bb-1)->t == t && bb->x > x0) bb-- ;
      break ;
    case TROU: /* trou are always referred downstream */
      if (bb->t != t && (bb - 1)->t == t)  bb-- ;
	break ;
    default: /* doubles , eliminated above */
      break ;
    }
  bb1 = bb ;

        /* count identical bases in area */
  nb = 0 ; ns = 0 ;
  i = c0 ;
  while (i > 0 && cc == ((arr(dna, i - 1, char)) & 0x0f) ) i-- ;
  while (i < cmax && cc == ((arr(dna, i, char)) & 0x0f)) { i++; nb++ ; }
      
  if (bb1->t == t) while ((bb1 - 1)->t == t) bb1-- ;
  bb2 = bb1 ; while (bb2->t == t) { ns++ ; bb2++ ; } 

  if (nb == ns)   /* success, bb is last base of correct type */
    switch (ep->type)
      {
      case AMBIGUE: 
      case ERREUR: 
	{ 
	  arr (base, s0, char) = BC_ACE | cc ;
	  /* reposition  after fixing */
	  s1 = s2 = s0 ; /* position de l'insert */
	  while (((arr (base, s1 - 1, char)) & 0x0f) == cc) s1-- ;
	  while (((arr (base, s2 + 1, char)) & 0x0f) == cc) s2++ ;
	  if (arr (basePos, s1 - 2, short) < (bb1 - 1)->x &&
	      arr (basePos, s1 + nb + 2, short) > (bb1 + nb + 1)->x )
	    { 
	      for (i = -1 ; i < nb + 2 ; i++)
		arr (basePos, s1 + i, short) = (bb1 + i)->x ;
	    }
 	  ep->type = TROU_DOUBLE ;
	  return TRUE ;
	}
      case INSERTION_DOUBLE: /*_AVANT: */
	return FALSE ;
      case INSERTION: /*_AVANT: */
	/* reposition, I have one less so i will delete base bb */
	arr (basePos, s1, short) = x1 + (x0 - x1)/3 ;
	arr (basePos, s2, short) = x2 + (x0 - x2)/3 ;
	/* delete */
	cp = arrp (base, s0, char) ;
	posp = arrp (basePos, s0, short) ;
	i = arrayMax(base) - s0 - 1 ;
	while (i--)
	  { *cp = *(cp + 1) ; *posp = *(posp + 1) ; posp++, cp++ ; }
	arrayMax(base)-- ;
	arrayMax(basePos)-- ;
	if (s0 < lane->clipTop)
	  lane->clipTop-- ;
	if (s0 < lane->clipEnd)
	  lane->clipEnd-- ;
	if (s0 < lane->clipExtend)
	  lane->clipExtend-- ;
	/* update the error table */
	ep1 = arrp(err, 0, A_ERR) - 1 ;
	i = arrayMax(err) ;
	while (ep1++, i--)
	  if (ep1->iShort > s0)
	    ep1->iShort-- ;
	ep->type = TROU_DOUBLE ;
	return TRUE ;

    case TROU:
                /* insert */
      cp = arrayp (base, arrayMax(base), char) ;
      posp = arrayp (basePos, arrayMax(basePos), short) ;
      i = arrayMax(base) - 1 ; /* note that arrayMax has now shifted */
      if ( (((arr (base, s0, char)) & 0x0f) != cc) &&
	  (((arr (base, s0 - 1, char)) & 0x0f) == cc))
	{ s0-- ; }
      s1 = s2 = s0 ; s1-- ; s2++ ;
      while (i > s2)
	{ *cp = *(cp - 1) ; *posp = *(posp - 1) ; posp--, cp-- ; i-- ;}
      array (base, s2, char) = cc | BC_ACE ;
      /* reposition,  i have one more so i will insert right of s0, at s2 */
      x1 = arr(basePos, s1, short) ;
      x0 = arr(basePos, s0, short) ;
      x2 = arr(basePos, s2, short) ;     
      arr (basePos, s0, short) = x1 + (x2 - x1)/3 ;
      arr (basePos, s2, short) = x2 - (x2 - x1)/3 ;
      
      /* reposition  after inserting */
      s1 = s2 = s0 + 1 ; /* position de l'insert */
      while (((arr (base, s1 - 1, char)) & 0x0f) == cc) s1-- ;
      while (((arr (base, s2 + 1, char)) & 0x0f) == cc) s2++ ;
      if (arr (basePos, s1 - 2, short) < (bb1 - 1)->x &&
	  arr (basePos, s1 + nb + 2, short) > (bb1 + nb + 1)->x )
	{ 
	  for (i = -1 ; i < nb + 2 ; i++)
	    arr (basePos, s1 + i, short) = (bb1 + i)->x ;
	}
               /* update the error table */
      ep1 = arrp(err, 0, A_ERR) - 1 ;
      i = arrayMax(err) ; s2 = s0 + 1 ; /* position de l'insert */
      while (ep1++, i--)
	if (ep1->iShort > s2)
	  ep1->iShort++ ;	
      if (s0 < lane->clipTop)
	lane->clipTop++ ;
      if (s0 <= lane->clipEnd)
	lane->clipEnd++ ;
      if (s0 <= lane->clipExtend)
	lane->clipExtend++ ;
      ep->type = TROU_DOUBLE ;
      return TRUE ;

      default: /* doubles , eliminated above */
	break ;
	
      }

  /* try now a neural net system for the remaining ? */
  if (ep->type == AMBIGUE)
    { char c0 = nnBaseCall (lane, s0) ;
      if ((c0 & 0x0f) == cc)
	{ array (base, s0, char) = (c0 | BC_ACE | BC_LOW) ;
	  /* reposition  after fixing */
	  s1 = s2 = s0 ; /* position de l'insert */
	  while (((arr (base, s1 - 1, char)) & 0x0f) == cc) s1-- ;
	  while (((arr (base, s2 + 1, char)) & 0x0f) == cc) s2++ ;
	  if (arr (basePos, s1 - 2, short) < (bb1 - 1)->x &&
	      arr (basePos, s1 + nb + 2, short) > (bb1 + nb + 1)->x )
	    { 
	      for (i = -1 ; i < nb + 2 ; i++)
		arr (basePos, s1 + i, short) = (bb1 + i)->x ;
	    }
 	  ep->type = TROU_DOUBLE ;
	  return TRUE ;
	}
    }
  return FALSE ;
}

/***************************************************/ 

static void seq2ace(LANE *lane) 
{ char *cp ;
  int n = lane->base ? arrayMax (lane->base) : 0 ;
  if (lane->scf < 3 || n <= 0)
    return ;
  arrayDestroy (lane->dna) ;
  defCptForget (0, lane->key) ;

  lane->dna = dnaCopy (lane->base) ;

  cp = arrp(lane->dna, 0, char) - 1 ;
  while (cp++, n--)
    *cp &= 0x0F ;
}

/***************************************************/

void laneEditSave (LOOK look, LANE *lane, int pos, int nb)
{ OBJ obj ;
  int x, old1, old2, i ;
  Array a = 0, units = 0 ;
  /*  static int countEdits = 0 ; */
  static int max = 0 ;
  KEY dnaKey ;
  BSunit *u ;
  LANE *lane1 ;

  if (lane->scf < 3)
    return ;
  seq2ace (lane) ;

  fMapTraceForget(lane->key) ;

  dnaSubClass(lane->key, &dnaKey) ;
  a = dnaCopy (lane->dna) ;
  dnaStoreDestroy (dnaKey, a) ;
	  
  if (look && look->mode != CDNA && (obj = bsUpdate(lane->key)))
    { if (!bsFindTag (obj, _Old_Clipping) &&
	  bsGetData (obj, _Clipping, _Int, &old1) &&
	  bsGetData (obj, _bsRight, _Int, &old2))      
	{ bsAddData (obj, _Old_Clipping, _Int, &old1) ;
	  bsAddData (obj, _bsRight, _Int, &old2) ;
	}
      units = arrayCreate (20, BSunit) ;
      if (pos && nb && bsFindTag (obj, _Clips) && bsFlatten (obj, 3, units))
	{ u = arrp (units, 0, BSunit) ;
	  for (i = 0 ; i < arrayMax (units) ; u += 3, i += 3)
	    { 
	      if ((x = u[1].i))
		{ if (pos < x)
		    x += nb ;
		  bsAddData (obj, u[0].k, _Int, &x) ;
		  if ((x = u[2].i))
		    { if (pos < x)
		      x += nb ;
		    bsAddData (obj, _bsRight, _Int, &x) ;
		    }
		}
	    }
	} /* le vecteur clip est ajoute ailleurs ? */
      else /* ? possible ? */
	{ x = lane->clipTop + 1 ;
	  bsAddData (obj, _Clipping, _Int, &x) ;
	  x = lane->clipEnd ;
	  bsAddData (obj, _bsRight, _Int, &x) ;
	  if (lane->vectorTop > 0 && lane->vectorEnd > 0)
	    { x = lane->vectorTop + 1 ;
	      bsAddData (obj, _Vector_Clipping, _Int, &x) ;
	      x = lane->vectorEnd ;
	      bsAddData (obj, _bsRight, _Int, &x) ;
	    }
	}
      bsSave (obj) ;
      arrayDestroy (units) ;
    }
  if (look && look->mode == CDNA && (obj = bsUpdate(lane->key)))
    { units = arrayCreate (20, BSunit) ;
      if (pos && nb && bsFindTag (obj, _Clips) && bsFlatten (obj, 3, units))
	{ u = arrp (units, 0, BSunit) ;
	  for (i = 0 ; i < arrayMax (units) ; u += 3, i += 3)
	    { 
	      if ((x = u[1].i))
		{ if (pos < x)
		    x += nb ;
		  bsAddData (obj, u[0].k, _Int, &x) ;
		  if ((x = u[2].i))
		    { if (pos < x)
		      x += nb ;
		    bsAddData (obj, _bsRight, _Int, &x) ;
		    }
		}
	    }
	}   
      if (pos && nb && bsGetData (obj, _PolyA_after_base, _Int, &x) && x > pos)
	{ x += nb ; bsAddData (obj, _PolyA_after_base, _Int, &x) ; }
      if (lane->vectorTop > 0 || lane->vectorEnd > 0)
	{ x = lane->vectorTop + 1 ;
	bsAddData (obj, _Vector_Clipping, _Int, &x) ;
	x = lane->vectorEnd ;
	if (x <= 0) x = arrayMax (lane->dna) ;
	bsAddData (obj, _bsRight, _Int, &x) ;
	}
      if (lane->qualityEnd > 0)
	{
	  if (bsGetData (obj, _Clipping, _Int, &x))
	    {
	      x = lane->qualityEnd ;
	      bsAddData (obj, _bsRight, _Int, &x) ;
	    }
	}

      bsSave (obj) ;
      arrayDestroy (units) ;
    }
  findXclipping (lane) ; /* to readjust the clipping */

  if (look && look && look->mode != CDNA)
    { max = arrayMax (look->dna) ;
      if (look->mode != CDNA && (obj = bsUpdate(look->key)))
	{ if (bsFindKey (obj, _Assembled_from, lane->key))
	    { 
	      if (look->sens < 0)
		{ x = max - lane->x1 ; 
		  bsAddData (obj, _bsRight, _Int, &x) ;
		  x = max - lane->x2 + 1 ;
		  bsAddData (obj, _bsRight, _Int, &x) ;
		}
	      else
		{ x = lane->x1 + 1 ;
		  bsAddData (obj, _bsRight, _Int, &x) ;
		  x = lane->x2 + 2 ;
		  bsAddData (obj, _bsRight, _Int, &x) ;
		}
	      x = lane->clipTop + 1 ;
	      bsAddData (obj, _bsRight, _Int, &x) ;
	      x = lane->clipEnd ;
	      bsAddData (obj, _bsRight, _Int, &x) ;
	    }
	  if (lane->isAligned && bsFindKey (obj, _Aligned, lane->key))
	    {
	      if (look->sens < 0)
		{ x = max - lane->x1 ;
		  bsAddData (obj, _bsRight, _Int, &x) ;
		  x = max - lane->x2 + 1 ;
		  bsAddData (obj, _bsRight, _Int, &x) ;
		}
	      else
		{ x = lane->x1 + 1 ;
		  bsAddData (obj, _bsRight, _Int, &x) ;
		  x = lane->x2 + 2 ;
		  bsAddData (obj, _bsRight, _Int, &x) ;
		}
	    }
	  bsSave (obj) ;
	}    
    }
  if (look) laneMakeErrArray(look, lane) ;
  defCptForget (0, lane->key) ;
  if (lane->scf >= 4)
    baseCallStore (lane) ;

  if (look && look->mode == CDNA  && pos) /* destroy everything to recover point edits */
    {
      i = arrayMax(look->lanes) ;
      while (i--)
	{
	  lane1 = arrp (look->lanes, i, LANE) ;
	  if (lane == lane1 || lane->key != lane1->key) continue ;
	  if (lane1->clipTop > pos) lane1->clipTop += nb ;
	  if (lane1->clipEnd > pos) lane1->clipEnd += nb ;
	  /*
	  if (lane1->x1 > pos) lane1->x1 += nb ;
	  if (lane1->x2 > pos) lane1->x2 += nb ;
	  if (lane1->x3 > pos) lane1->x3 += nb ;
	  */
	  if (lane1->vectorTop > pos) lane1->vectorTop += nb ;
	  if (lane1->vectorEnd > pos) lane1->vectorEnd += nb ;
	  /*
	  if (lane1->dna)
	    {
	      arrayDestroy (lane1->dna) ;
	      lane1->dna = dnaCopy (lane->dna) ;
	    }
	  laneMakeErrArray(look, lane1) ;
	  */
	  arrayDestroy (lane1->errArray) ;
	  arrayDestroy (lane1->dna) ;
	  arrayDestroy (lane1->base) ;
	  arrayDestroy (lane1->basePos) ;
	  arrayDestroy (lane1->baseQuality) ;
	  arrayDestroy (lane1->baseCall) ;
	  lane1->scf = 0 ;
	  traceGetLane (look, lane1) ;
	}
      if (nb && (obj = bsUpdate (look->gene)))
	{
	  units = arrayCreate (20, BSunit) ;
	  if (bsGetArray (obj, _Assembled_from, units, 5))
	    { 
	      u = arrp (units, 0, BSunit) ;
	      for (i = 0 ; i < arrayMax (units) ; u += 5, i += 5)
		{ 
		  if (lane->key != u[2].k) continue ;
		  if ((x = u[3].i))
		    { if (pos < x)
		      u[3].i += nb ;
		    }
		  if ((x = u[4].i))
		    { if (pos < x)
		      u[4].i += nb ;
		    }
		}
	    }
	  bsAddArray (obj,_Assembled_from, units, 5) ;
	  bsSave (obj) ;
	  arrayDestroy (units) ;
	}
    }
}

/***************************************************/

int baseCallPatchLane (Array dnaDirect, Array dnaReverse,
			LANE *lane)
{ BOOL f1, f2 = FALSE, isUp ;
  Array dna, err = 0 ;
  A_ERR *ep ;
  int nn = 0, i, tour, x1, x3, ctop, c3, amax = arrayMax(dnaDirect) ;

  if (!baseCallGet (lane))      
    return FALSE ;
  if (! arrayExists(lane->baseCall) &&
      ! findBaseCall (lane))
    return FALSE ;
 
  if (! arrayMax(lane->baseCall))
    return FALSE ;

  isUp = lane->x1 > lane->x2 ? TRUE : FALSE ;

  if (!isUp)
    dna = dnaDirect ;
  else
    { dna = dnaReverse ;
      if (!dna || lane->x1 > amax)
	return FALSE ;
      lane->x1 = amax - lane->x1 - 1 ; 
      lane->x2 = amax - lane->x2 - 1 ; 
      lane->x3 = amax - lane->x3 - 1 ; 
    }

  f1 = TRUE ; tour = 2 ; 
  while (tour-- && f1)
    { f1 = FALSE ;
      x1 = lane->x1 ; ctop = lane->clipTop ;
      x3 = lane->x3 - 1 ; c3 = lane->clipExtend - 1 ;  
      if (x3 >= arrayMax(dnaDirect))
	x3 = arrayMax(dnaDirect) - 1;
      err = dnaAlignCptErreur (dna, lane->dna, 
			       &x1, &x3,
			       &ctop, &c3) ;
      if (!err || ! arrayMax(err))
	break ;
  
      lane->x1 = x1 ;
      lane->x3 = x3 + 1 ;
      lane->clipTop = ctop ;
      lane->clipExtend = c3 + 1 ;
  
      for (i = 0, ep = arrp(err, 0, A_ERR) ; i < arrayMax(err) ; ep++, i++)
	{
	  if (ep->iLong >= arrayMax(dna))
	    break ;
           /* try only if err is reliable */
	  if (i < 6 || 
	      i > arrayMax(err) - 6 ||
	      (ep + 4)->iLong - (ep - 4)->iLong > 18)
	    if (patchError (dna, lane, ep, err))
	      { nn++ ; f1 = TRUE ; }
	}
      if (f1) seq2ace (lane) ;
      arrayDestroy (err) ;
      f2 |= f1 ;
    }

  if (isUp)
    { lane->x1 = amax - lane->x1 - 1 ; 
      lane->x2 = amax - lane->x2 - 1 ; 
      lane->x3 = amax - lane->x3 - 1 ; 
    }
  arrayDestroy (err) ;
  if (f2) monDnaForget (0, lane->key) ;
  return nn ;
}

/*****************************************************/

void baseCallRedoLaneBaseCall (LANE *lane, int *np)
{ 
  int i, dummy ;
  Array dna = 0 , bp = 0 ;
  BASECALL *bb ; char *d ; short *bpp ;
  OBJ obj ;
 
  bb = arrp (lane->baseCall, 0, BASECALL) - 1 ;
  i = arrayMax (lane->baseCall) ;
  dna = arrayCreate (i, char) ;
  bp = arrayCreate (i, short) ;
  array (dna, i - 1, char)= 0 ; /* make room */
  array (bp, i - 1, short)= 0 ; /* make room */
  d = arrp (dna, 0, char) - 1 ;
  bpp = arrp (bp, 0, short) - 1 ;
  while (d++, bpp++, bb++, i--)
    { 
      switch (bb->t)
	{
	case 0 : *d = A_ ; break ;
	case 1 : *d = G_ ; break ;
	case 2 : *d = C_ ; break ;
	case 3 : *d = T_ ; break ;
	default : *d = N_ ; break ;
	}
      *bpp = bb->x ;
    }
  arrayDestroy (lane->base) ;
  lane->base = dna ;  
  arrayDestroy (lane->basePos) ;
  lane->basePos = bp ;
  if (np) (*np)++ ; 
  seq2ace(lane) ;
  laneEditSave (0, lane, 0, 0) ;
  /* rm the clips */  
  if ((obj = bsUpdate(lane->key)))
    { if (bsFindTag (obj, _Clips))
       bsRemove (obj) ;
      bsSave (obj) ;
    } 
  baseCallUnclipLane(lane, (KEY)'G', &dummy) ;
  lane->vectorTop = lane->vectorEnd = 0 ;
  trackVector (0, lane->key, TRUE) ;
}

/***************************************************/ 

BOOL baseCallRedoBaseCall (KEY key, int *np)
{ static LANE *lane = 0 ;
  
  if (!lane)
    lane = (LANE*) messalloc (sizeof (struct LaneStruct)) ;
  lane->key = key ;

  if (!baseCallGet (lane) ||
      !findBaseCall (lane) ||
      !arrayMax(lane->baseCall))
    { 
      laneDestroy (lane) ;
      return 0 ;
    }

  baseCallRedoLaneBaseCall (lane, np) ;
  laneDestroy (lane) ;
  return TRUE ;
}

/***************************************************/ 

void findXclipping (LANE *lane)
{ int n, nBase, top, end, extend ;
  short *basePos ;

  if (!lane->base && !baseCallGet(lane))
    return ;   
  if (!lane->base)
    return ;
  if (lane->clipEnd > arrayMax (lane->dna))
    lane->clipEnd = arrayMax (lane->dna) ;
  top = lane->clipTop ;
  end = lane->clipEnd ;
  extend =  end + 20 ;
  if (extend >= arrayMax (lane->basePos))
    extend = arrayMax (lane->basePos) ;
  lane->clipExtend = extend ;
  lane->x3 = lane->x2 + (lane->x1 < lane->x2 ? extend - end : end - extend) ;

  n = 1 ; nBase = arrayMax (lane->basePos) - 3 ;
  if (arrayMax (lane->basePos) < 3)
    return ;
  basePos = arrp (lane->basePos, 2, short) ;

  if (top)
    { while(n < top && nBase)
	{ basePos++  ; nBase-- ; n++ ;
	}
      lane->xClipTop = (*(basePos - 1) + *(basePos - 2))/ 2 ;
    }
  else
    lane->xClipTop = 0 ;
  while(n < end && nBase)
    { basePos++  ; nBase-- ; n++ ;
    }
  if (lane->seq)
    {
      lane->xClipEnd = 
	nBase ?
	(*(basePos - 2) + *(basePos - 1))/ 2 :  seqMax (lane->seq) ;
      while(n < extend)
	{ basePos++  ; nBase-- ; n++ ;
	}
      lane->xClipExtend =  nBase ?
	(*(basePos - 2) + *(basePos - 1))/ 2 :  seqMax (lane->seq) ;
    }
}

/***************************************************/ 

#ifndef NON_GRAPHIC

static void laneDoPlotPeriodicite (LANE *lane, Array h2, Array h3)
{ Array 
    bc = lane->baseCall ;
  BASECALL *bb, *bb1 ;
  int 
    j, j1, deb, dx, dxm, dxold[10], fin ;

  if (lane->scf < 4)
    return ;

  j = 10 ;
  while (j--)
    dxold[j] = 0 ;

  dxm = 0 ;

  dxm = 0 ; fin = arrayMax(bc) ;
  for  (deb  = 1, bb1 = 0,  
	bb = arrp(bc, deb , BASECALL) ;
	deb < fin - 2  ; bb++, deb++)
    { if (bb->flag & BC_LOW) ; /* keep same bb1 */
      else /* bb is HIGH */
	{ bb1 = bb + 1 ; j1 = deb + 2 ;
	  while ((bb1->flag & BC_LOW) && j1 < fin) 
	    { j1++ ; bb1++ ;
	    }
	}
      if (!bb1)
	continue ;
      dx = bb1->x - bb->x ;
      j = 10 ;
      while (--j)
	dxold[j] = dxold[j-1] ;
      dxold[0] = dx ;
      dxm = dxm + dx - dxold[9] ;
      if (h2) array (h2, deb, int ) = dx ;
      if (h3) array (h3, deb, int ) = dxm ;
    }
}

/***************************************************/ 

void lanePlotDerivee (LOOK look, LANE *lane)
{ int i, x1, x2, mm, *ip ;
  Array 
    hh, 
    basePos = lane->basePos ;
  Read* seq = lane->seq ;

  if (lane->scf < 4)
    return ;
  
  hh = arrayCreate (1000, int) ;

  for (i = 11 ; i < arrayMax (basePos) - 11 ; i++)
     { x1 = arr (basePos, i - 10, short)  ;
       x2 = arr (basePos, i + 10, short)  ;
       array (hh, i, int) = seqEnergyOfDerivee (seq, x1, x2) ;
     }
  i = arrayMax (hh) ; mm = 0 ;
  ip = arrp (hh, 0, int) - 1 ;
  while (ip++, i--)
    if (*ip > mm) mm = *ip ;
     /* mark the clippings */
  array (hh, lane->clipTop, int) = mm ;
  array (hh, lane->clipEnd - 1, int) = mm ;


  plotHisto(messprintf("Energie de la derivee de %s",name(lane->key)), hh) ;
}

/***************************************************/ 

void lanePlotPeriodicite (LOOK look, LANE *lane)
{ Array 
    h2, h3 ;

  if (lane->scf < 4)
    return ;

  h2 = arrayCreate (1000, int) ;
  h3 = arrayCreate (1000, int) ;

  laneDoPlotPeriodicite (lane, h2, h3) ;
     /* mark the clippings */
  array (h2, lane->clipEnd - 1, int) = 0 ;
  array (h2, lane->clipTop, int) = 0 ;
  array (h3, lane->clipEnd - 1, int) = 0 ;
  array (h3, lane->clipTop, int) = 0 ;

  plotHisto(messprintf("Periodicite de  %s",
		       name(lane->key)), h2) ;
  plotHisto(messprintf("Periodicite moyennee de  %s",
		       name(lane->key)), h3) ;
  baseCallFindClips(lane, TRUE) ;
  lanePlotDerivee (look, lane) ;
}

#endif

/*****************************************************/
#ifdef JUNK
BOOL baseCall (KEY key)
{ OBJ obj ;
  LANE *lane ;
  
  if (!isWriteAccess ())	/* may occur is somebody else grabed it */
    { messout("Sorry, you do not have write access") ;
      return FALSE ;
    }

  if (!key ||
      class(key) != _VSequence ||
      !(obj = bsUpdate(key)))
    return FALSE ;
  bsSave (obj) ;

  lane = (LANE*) messalloc (sizeof (struct LaneStruct)) ;
  lane->key = key ;
  lane->x1 = 0 ; lane->x2 = 100 ;

  if (traceGetLane (0, lane))
    { findXclipping (lane) ; /* to adjust the clipping */
      if (aceBaseCallLane (lane))
	laneEditSave (0, lane, 0, 0) ;
      laneDestroy (lane) ;
    }
  messfree (lane) ;

  return FALSE ;
}
#endif
/***********************************************************/
/***********************************************************/
#ifdef JUNK
static BOOL baseCallFindClipsOld (LANE *lane)
{ Array bc ;

  int x, i, *ip, xe, xg, xf, xe1, xg1, xf1, maxDna, seuil ;
  OBJ obj = bsCreate (lane->key) ;
  Array h3 = 0 ;
  short *bpos ; 
  KEY dummy ;

  if (!obj) /* why bother */
    return FALSE ;
  if (bsFindTag (obj, _Excellent_upto))  /* Already done */
    { bsDestroy (obj) ; return TRUE ; }
  
  if (!traceGetLane (0, lane) || 
      !baseCallGetSeq(lane) ||
      !baseCallGet(lane) ||
      ! findBaseCall (lane) ||
      lane->scf < 4 ||
      arrayMax(lane->basePos) < 150 ||
      !bsGetKey (obj, _DNA, &dummy) ||
      !bsGetData (obj, _bsRight, _Int, &maxDna))
    { bsDestroy (obj) ; return FALSE ; }
 /* divide by 10, because i average 200 bases on 10 neighbours */
  seuil = (arr(lane->basePos, 150, short) - arr(lane->basePos, 50, short)) / 10 ;
  h3 = arrayCreate (1000, int) ;
  laneDoPlotPeriodicite (lane, 0, h3) ;
  xe = xg = xf = 0 ;
  for (i = 80, ip = arrp(h3, 80, int) ; i < arrayMax(h3) ; ip++, i++)
    { 
      if (!xe && *ip > 2 * seuil) xe = i ;
      if (!xg && *ip > 3 * seuil) xg = i ;
      if (!xf && *ip > 4 * seuil) xf = i ;
      if (xf && i < xf + 60 && *ip < 3 * seuil)
	xf = i ;
    }
  if (!xe) xe = arrayMax(lane->baseCall) - 2 ;
  if (!xg) xg = arrayMax(lane->baseCall) - 2 ;
  if (!xf) xf = arrayMax(lane->baseCall) - 2 ;
  arrayDestroy (h3) ;
  
  i = 0 ; bc = lane->baseCall ;
  x = arrp(bc, xe, BASECALL)->x ; xe1 = 0 ;
  for (bpos = arrp(lane->basePos, i, short); i < arrayMax(lane->basePos) ; bpos++, i++)
    if (x < *bpos) { xe1 = i - 1 ; break ; } 
  x = arrp(bc, xg, BASECALL)->x ; xg1 = 0 ;
  for (bpos = arrp(lane->basePos, i, short); i < arrayMax(lane->basePos) ; bpos++, i++)
    if (x < *bpos) { xg1 = i - 1 ; break ; } 
  x = arrp(bc, xf, BASECALL)->x ; xf1 = 0 ;
  for (bpos = arrp(lane->basePos, i, short); i < arrayMax(lane->basePos) ; bpos++, i++)
    if (x < *bpos) { xf1 = i - 1 ; break ; } 

  if (!xe1 || xe1 > maxDna) xe1 = maxDna ;
  if (!xg1 || xg1 > maxDna) xg1 = maxDna ;
  if (!xf1 || xf1 > maxDna) xf1 = maxDna ;
  
  /* now edit the lane */
  bsDestroy (obj) ;  /* destroy/update: some of above code needs update mode */
  obj = bsUpdate (lane->key) ;
  bsAddData (obj, _Excellent_upto, _Int, &xe1) ;
  bsAddData (obj, _Good_upto, _Int, &xg1) ;
  bsAddData (obj, _Fair_upto, _Int, &xf1) ;

  bsSave (obj) ;
  return TRUE ;
}
#endif

static BOOL baseCallFindClips (LANE *lane, BOOL doPlot)
{ BASECALL *bb ;
  int 
    i1, i2, nn, xe, xg, xf,  maxDna, maxbc, maxbase,
    iabi, ibc, ibc0, ngood, nabi, nbest, nbest200, xbest, xbegin ;
  OBJ obj = bsCreate (lane->key) ;
  Array h3 = 0 ;
  short *bp ; 
  KEY dummy ;

  if (!obj) /* why bother */
    return FALSE ;
  if (!doPlot && bsFindTag (obj, _Excellent_upto))  /* Already done  */
    { bsDestroy (obj) ;
      return TRUE ; 
    }

  if (!traceGetLane (0, lane) || 
      !baseCallGetSeq(lane) ||
      !baseCallGet(lane) ||
      ! findBaseCall (lane) ||
      lane->scf < 4 ||
      arrayMax(lane->basePos) < 150 ||
      !bsGetKey (obj, _DNA, &dummy) ||
      !bsGetData (obj, _bsRight, _Int, &maxDna))
    { bsDestroy (obj) ; return FALSE ; }

  h3 = arrayCreate (200, int) ;
  maxbase = arrayMax(lane->base) ;
  maxbc = arrayMax(lane->baseCall) ;

  /* new system: count nuber of good in best 200 bp stretch */
  ibc0 = 0 ; nabi = 200 ;
  nbest200 = 0 ;
  for (iabi = 40 ; iabi < maxbase - nabi && iabi < 400 - nabi ; iabi++)
    { bp = arrp(lane->basePos, iabi, short) ;
      i1 = *(bp) ;
      i2 = *(bp + nabi) ;
      ngood = 0 ; 
      bb = arrp(lane->baseCall, ibc0, BASECALL) ;
      ibc = ibc0 ;
      while (ibc < maxbc && bb->x < i1)
	{ ibc++ ; bb++ ; ibc0++ ; }
      while (ibc < maxbc && bb->x < i2)
	{ if (!(bb->flag & BC_LOW)) ngood ++ ;
	  ibc++ ; bb++ ;
	}
      if (ngood > nbest200)
	{ 
	  nbest200 = ngood ;
	}
      if (!ngood)
	break ;
    }

  /* first, i look for the best area */
  ibc0 = 0 ; nabi = 40 ;
  xbest = maxbase/4 ; nbest = 2*nabi ;
  for (iabi = 40 ; iabi < maxbase - nabi ; iabi++)
    { bp = arrp(lane->basePos, iabi, short) ;
      i1 = *(bp) ;
      i2 = *(bp + nabi) ;
      ngood = 0 ; 
      bb = arrp(lane->baseCall, ibc0, BASECALL) ;
      ibc = ibc0 ;
      while (ibc < maxbc && bb->x < i1)
	{ ibc++ ; bb++ ; ibc0++ ; }
      while (ibc < maxbc && bb->x < i2)
	{ if (!(bb->flag & BC_LOW)) ngood ++ ;
	  ibc++ ; bb++ ;
	}
      nn = nabi - ngood ;
      if (nn < nbest)
	{ xbest = iabi ;
	  nbest = nn ;
	}
      if (!nn)
	break ;
    }

  /* now i look for the end of the good zone */
  xe = xg = xf = 0 ; ibc0 = 0 ; nabi = 20 ;
  for (iabi = xbest ; iabi < maxbase - nabi ; iabi++)
    { bp = arrp(lane->basePos, iabi, short) ;
      i1 = *(bp) ;
      i2 = *(bp + nabi) ;
      ngood = 0 ; 
      bb = arrp(lane->baseCall, ibc0, BASECALL) ;
      ibc = ibc0 ;
      while (ibc < maxbc && bb->x < i1)
	{ ibc++ ; bb++ ; ibc0++ ; }
      while (ibc < maxbc && bb->x < i2)
	{ if (!(bb->flag & BC_LOW)) ngood ++ ;
	  ibc++ ; bb++ ;
	}
      nn = nabi - ngood ;
      if (nn < 0) nn = -nn ;
      if (!xe && 100*nn > 20*nabi) xe = iabi  ;
      if (!xg && (100*nn > 40*nabi || ( 100*nn > 30*nabi && iabi > xe + 30)))
	xg = iabi ;
      if (!xf && (100*nn > 50*nabi || ( 100*nn > 40*nabi && iabi > xg + 30)))
	xf = iabi ;
      if (nabi) array(h3, iabi, int) = 100 * nn/nabi ;
      if (!doPlot && xf) break ;
    }

  /* now i look for the beginning of the good zone */
  ibc0 = maxbc - 5 ; xbegin = 20 ; nabi = 20 ;
  /* if (xbest > 100) xbest = 100 ; */
  for (iabi = xbest ; iabi > 1 ; iabi--)
    { bp = arrp(lane->basePos, iabi, short) ;
      i1 = *(bp) ;
      i2 = *(bp + nabi) ;
      ngood = 0 ; 
      bb = arrp(lane->baseCall, ibc0, BASECALL) ;
      ibc = ibc0 ;
      while (ibc > 0 && bb->x > i2 + 2)
	{ ibc-- ; bb-- ; ibc0-- ; }
      while (ibc > 0 && bb->x > i1 -3 )
	{ if (!(bb->flag & BC_LOW)) ngood ++ ;
	  ibc-- ; bb-- ;
	}
      nn = nabi - ngood ;
      if (nn < 0) nn = -nn ;
      if (100*nn > 40*nabi)
	{ xbegin = iabi + 2 ;
	  break ;
	}
    }

  if (!xe) xe = arrayMax(lane->base) - 5 ;
  if (!xg) xg = arrayMax(lane->base) - 3 ;
  if (!xf) xf = arrayMax(lane->base) - 1 ;
  
  /* now edit the lane */
  bsDestroy (obj) ;  /* destroy/update: some of above code needs update mode */
  if ((obj = bsUpdate (lane->key)))
    { int ctop, cend ;

      bsAddData (obj, _Excellent_upto, _Int, &xe) ;
      bsAddData (obj, _Good_upto, _Int, &xg) ;
      bsAddData (obj, _Fair_upto, _Int, &xf) ;
      bsAddData (obj, str2tag("Best200"), _Int, &nbest200) ;
      ctop = 0 ; cend = 0 ;
      if (bsGetData (obj, _Clipping, _Int, &ctop))
	bsGetData (obj, _bsRight, _Int, &cend) ;
      xbegin++ ;  /* np zero */
      if (ctop < xbegin)
	{ bsAddData (obj, _Clipping, _Int, &xbegin) ;
	  if (cend)
	    bsAddData (obj, _bsRight, _Int, &cend) ;
	}
      
      bsSave (obj) ;
    }  

  /* insure that the dna is identical to that in the trace */
  if (lane->baseCall && arrayMax(lane->baseCall) > 0)
    {
      Array dnaLane = 0 ;
      Array dna = dnaGet (lane->key) ;
      char *cp, *cq ;
      int i, j, imax = arrayMax(lane->base) ;
      BOOL ok = FALSE ;
      
      dnaLane = arrayCreate (imax + 1, char) ;
      array (dnaLane, imax, char) = 0 ; arrayMax(dnaLane) = imax ; /* make zero terminated dna */

      cq = arrp (dnaLane, 0, char) ;
      cp = arrp(lane->base, 0, char) ;
      for (i = 0 ; i < imax ; cq++, cp++, i++)
	*cq = (*cp) & 0xf ;

      i = dna ? arrayMax(dna) : 0 ;
      j = arrayMax(dnaLane) ;
      
      if (i == j && i > 0)
	{
	  cp = arrp(dna, 0, char) ;
	  cq = arrp(dnaLane, 0, char) ;
	  while (i-- && *cp++ == *cq++) ;
	  if (i == -1) ok = TRUE ;
	}
      if (1)
	{
	  if (ok) printf("BaseCall %s ok\n", name(lane->key)) ;
	  else printf("Est %s restoring lost basecall editions\n", name(lane->key)) ;

	}
      if (!ok) 
	{
	  defCptForget (0, lane->key) ;
	  dnaStoreDestroy (lane->key, dnaLane) ;
	}
      else
	arrayDestroy (dnaLane) ;
      arrayDestroy (dna) ;
    }

#ifndef NON_GRAPHIC
  if (doPlot) 
    plotHisto ("pourcent de bon", h3) ;
  else
    arrayDestroy (h3) ;
#else
  arrayDestroy (h3) ;
#endif
  return TRUE ;
}

/***********************************************************/

BOOL baseCallUnclipLane (LANE *lane, KEY type, int *dxp)
{ OBJ obj ;
  int x, dx, dx1, clipTop = lane->clipTop, newClipTop = lane->clipTop, clipEnd = lane->clipEnd ;
  int vClipTop, vClipEnd ;

  if (!baseCallFindClips (lane, FALSE))
    return FALSE ;
  
  obj = bsUpdate (lane->key) ;  /* read after FindClips to see the results */
  if (!obj)
    return FALSE ;
  x = clipEnd ;
  if (!newClipTop)
    {
      bsGetData (obj, _Clipping, _Int, &newClipTop) ;
      newClipTop-- ;
    }
  switch (ace_upper(type & 0xff))
    {
    case 'E':  bsGetData (obj, _Excellent_upto, _Int, &x); break ;
    case 'G':  bsGetData (obj, _Good_upto, _Int, &x) ; break ;
    case 'F':  bsGetData (obj, _Fair_upto, _Int, &x) ; break ;
    case 'H':  
      if (bsGetData (obj, _Hand_Clipping, _Int, &newClipTop) &&
	  bsGetData (obj, _bsRight, _Int, &x))
	{ newClipTop-- ; }
      else 
	{ bsDestroy (obj) ;
          return FALSE ;
	}
      break ;
    }
  vClipTop = vClipEnd = -1 ;
  if (bsGetData (obj, _Vector_Clipping, _Int, &vClipTop) &&
      bsGetData (obj, _bsRight, _Int, &vClipEnd))
    if (x > vClipEnd)
      x = vClipEnd ;
  vClipTop-- ; /* info */
  dx = x - clipEnd ;
  if (lane->upSequence)
    lane->x2 -= dx ;
  else
    lane->x2 += dx ;
  lane->x3 = lane->x2 ;

  if (vClipTop > 5)
    newClipTop = vClipTop ;
  dx1 = newClipTop - clipTop ;
  if (lane->upSequence)
    lane->x1 -= dx1 ;
  else
    lane->x1 += dx1 ;
  lane->clipTop = newClipTop ;
  lane->clipEnd = lane->clipExtend = x ;
  
  x = lane->clipTop + 1 ; /* No zero */
  if (x<1) x = 1 ;
  bsAddData (obj, _Clipping, _Int, &x) ;
  x = lane->clipEnd ;
  if (x < 1) x = 1 ;
  bsAddData (obj, _bsRight, _Int, &x) ;

  bsSave (obj) ;
  *dxp += dx ;
  if (dx || dx1)
    defCptForget (0, lane->key) ;
  return dx != 0 ? TRUE : FALSE ;
}

/***********************************************************/

int baseCallUnclipContig (KEY contig, KEY type, int *dxp)
{ OBJ obj = 0 ;
  Array aa = 0 ;
  int i, nn = 0 ;
  static LANE *lane = 0 ;
  KEY dnaKey ;
  BSunit *uu ;

  if (class(contig) == _VDNA)
    { dnaKey = contig ;
      dnaReClass (dnaKey, &contig) ;
    }

  if (!lane)
    lane = (LANE*) messalloc (sizeof (struct LaneStruct)) ;
  aa = arrayCreate (100, BSunit) ;

  obj = bsCreate (contig) ;
  if (!obj || !bsGetArray (obj, _Assembled_from, aa, 5))
    goto abort ;
  bsDestroy (obj) ;

  messStatus ("Unclip") ;

  for (i = 0 ; i < arrayMax(aa) ; i+= 5)
    { uu = arrp (aa, i, BSunit) ;
      lane->key = uu[0].k ;
      lane->x1 = uu[1].i - 1 ;
      lane->x2 = uu[2].i ;
      if (lane->x1 <= lane->x2) 
	{ lane->upSequence = FALSE ; }
      else 
	{ lane->upSequence = TRUE ; lane->x2 -= 2 ; }
      lane->x3 = lane->x2 ;
      lane->clipTop = uu[3].i - 1 ;
      lane->clipEnd = uu[4].i ;
      if (baseCallUnclipLane (lane, type, dxp))
	{ nn++ ;
	  if (lane->upSequence)
	    { uu[1].i = lane->x1 + 1 ;
	      uu[2].i = lane->x2 + 2 ;
	    }
	  else
	    { uu[1].i = lane->x1 + 1 ;
	      uu[2].i = lane->x2 ;
	    }
	  uu[3].i = lane->clipTop + 1 ;
	  uu[4].i = lane->clipEnd ;
          
	  laneEditSave (0, lane, 0, 0) ;
	}
      laneDestroy (lane) ;
    }
  if ((obj = bsUpdate (contig)))
    { bsAddArray (obj, _Assembled_from, aa, 5) ;
      bsSave (obj) ;
    }
 abort:
  arrayDestroy (aa) ;
  bsDestroy (obj) ;
  return nn ;
}

/***********************************************************/

static int baseCallUnclipRead (KEY seq, KEY type, int *dxp)
{ static LANE *lane = 0 ;
  int done = 0 ;

  messStatus ("Unclip") ;
  if (!lane)
    lane = (LANE*) messalloc (sizeof (struct LaneStruct)) ;

  lane->key = seq ;
  dnaAlignGetClip (0, 0, seq, &lane->clipTop, &lane->clipEnd) ;
      
  lane->clipTop-- ;  /* bio */
  if (type)
    { if (baseCallUnclipLane (lane, type, dxp))
	done = 1 ;
    }
  else
    baseCallFindClips (lane, FALSE) ;
  laneDestroy (lane) ;
  return done ;
}

/***********************************************************/

int baseCallUnclipKeySet (KEY link, KEYSET ks, KEY type, int *dxp)
{ KEYSET ks1 = keySetCopy (ks), ks2, ks3, ks4, ks5 ; 
  KEY *kp ; int i, nn = 0, n ;

  i = keySetMax (ks1) ; kp = arrp(ks1, 0, KEY) - 1 ;
  while (kp++, i--)
    if (class(*kp) == _VDNA)
      dnaReClass (*kp, kp) ;

  keySetSort (ks1) ;
  ks2 = query (ks1, "CLASS Read") ;
  i = keySetMax (ks2) ;
  if (i)
    { kp = arrp(ks2, 0, KEY) - 1 ;
      while (kp++, i--)
	nn += baseCallUnclipRead (*kp, type, dxp) ;
    }
  ks3 = keySetMINUS (ks1, ks2) ;

  ks4 = query (ks3, ">Subsequence") ;
  ks5 = keySetOR (ks3, ks4) ;
  i = keySetMax (ks5) ; 
  if (i)
    { kp = arrp(ks5, 0, KEY) - 1 ;
      while (kp++, i--)
	{ n = baseCallUnclipContig (*kp, type, dxp) ;
	  if (n)
	    { nn += n ;
	      dnaAlignFixContig (link, *kp) ;
	    }
	}
      dnaAlignAdjustLink (link) ;
    }
  keySetDestroy (ks1) ;
  keySetDestroy (ks2) ;
  keySetDestroy (ks3) ;
  keySetDestroy (ks4) ;
  keySetDestroy (ks5) ;

  return nn ;
}

/***********************************************************/

int baseCallClipContig2Max (KEY link, int max, int *dxp)
{ int n = 0, i, j, dx, mv ;
  KEY *kp ;
  KEYSET ks = 0 ;
  Array aa ;
  BSunit *u ;
  OBJ obj = 0 ;

  if (!link)
    return n ;
  ks = queryKey (link, ">Subsequence") ;
  if (!(i = keySetMax (ks)))
    goto abort ;
  kp = arrp (ks, 0, KEY) - 1 ;
  aa = arrayCreate (100, BSunit) ;
  while (kp++, i--)
    { if ((obj = bsCreate (*kp)) && bsGetArray (obj, _Assembled_from, aa, 5))
	{ bsDestroy (obj) ;
	  u = arrp (aa, 0, BSunit) ;
	  mv = 0 ;
	  for (j = 0 ; j < arrayMax (aa) ; j += 5, u += 5)
	    if (u[4].i > max)
	      { dx = max - u[4].i ;
		u[4].i = max ;
		n++ ;
		mv = 1 ;
		*dxp += dx ;
		if (u[1].i < u[2].i)
		  u[2].i += dx ;
		else
		  u[2].i -= dx ;
		defCptForget (link, u[0].k) ;
	      }
	  if (mv && (obj = bsUpdate (*kp)))
	    { bsAddArray (obj, _Assembled_from, aa, 5) ;
	      bsSave (obj) ;
	      alignToolsAdjustContigSize (link, *kp) ;
	    }
	}
      bsDestroy (obj) ;
    }
  arrayDestroy (aa) ;
  if (n)
    dnaAlignAdjustLink (link) ;
 abort:
  keySetDestroy (ks) ;
  return n ;
}

/***********************************************************/
/************************ Tiling ***************************/
/***********************************************************/
typedef struct {KEY key ; int a, b, t, e, flag ;} TILE ;

static int tileOrder (const void *a, const void *b)
{ const TILE *t1 = a, *t2 = b ;

  return t1->a - t2->a ;
}

/***********************************************************/

static BOOL getNextCandidate (Array tiles, TILE *t1, TILE *t2, TILE **tr)
{ TILE *t, *tc ;
  int i, x1 = 0 ;

  i = arrayMax (tiles) ;
  if (i < 2) return FALSE ;
  t = arrp (tiles, 0, TILE) ;
  if (t1 < t || t2 <= t1 || (t2 - t) >= i)
    { messerror ("Appel incorect de getNextCandidate : t1 | t2 ne sont plus dans tiles") ;
      return FALSE ;
    }
  tc = *tr = t1 ;
  if (!t1->flag) x1 = t1->b ;
  for (t = tc + 1 ; t < t2 ; t++)
    if (!t->flag && t->b > x1)
      { x1 = t->b ; tc = t ; }
  if (tc == t1 && t1->flag) /* mean nothing match */
    return FALSE ;
  else 
    { *tr = tc ;
      return TRUE ;
    }
}

/***********************************************************/

static int doTile (Array tiles, KEY type, int *dxp)
{ int i, dx, f, x1, x2, vtop, vend ;
  TILE *t, *t1, *t2 ;
  OBJ obj = 0 ;
  int found = 0 ;

  if (arrayMax(tiles) < 2) 
    return 0 ;
  t1 = arrp (tiles, 0, TILE) ;
  t2 = arrp (tiles, 0, TILE) ;
  
  i = arrayMax(tiles) - 1 ;
  while (t2++, i--)
    { /* cherche la plus longue de t1 a t2 */
      while (getNextCandidate (tiles, t1, t2, &t) && (t2->a > t->b))
	{ vtop = vend = f = 0 ;
	  if ((obj = bsCreate (t->key)))
	    { if (bsGetData (obj, _Vector_Clipping, _Int, &vtop))
		bsGetData (obj, _bsRight, _Int, &vend) ;
	      if (bsGetData (obj, _Fair_upto, _Int, &f))
		if (vend && f > vend)
		  f = vend ;
	      bsDestroy (obj) ;
	      if (f > t->e)
		{ x1 = t2->a - t->b + 20 ; /* i would want that */
		  x2 = f - t->e ; /* i have that */
		  if (x1 - x2 > 30)
                    { t->flag = 1 ;
                      continue ;
                    }
		  dx = x1 < x2 ? x1 : x2 ;
		  t->e += dx ;
		  t->b += dx ;
		  *dxp += dx ;
		  found++ ;
		  if (t->b > t2->a)
		    break ;
		}
	      t->flag = 1 ;
	    }
	  bsDestroy (obj) ;
	}
    }
  return found ;
}

/***********************************************************/

int baseCallTileContig (KEY contig, KEY type, int *dxp)
{ int i, j, mx = 0, nmvd = 0 ;
  Array aa = arrayCreate(900, BSunit), tiles = arrayCreate (300, TILE) ;
  OBJ Contig ;
  KEY dummy ;
  BSunit *u ;
  TILE *t ;

  Contig = bsUpdate (contig) ;
  if (!Contig || 
      !bsGetKey (Contig, _DNA, &dummy) ||
      !bsGetData (Contig, _bsRight, _Int, &mx) ||
      mx < 20 ||
      !bsFindTag (Contig, _Assembled_from) ||
      !bsFlatten (Contig, 5, aa))
    goto abort ;

  if (bsFindTag (Contig, _Assembled_from)) 
    bsRemove (Contig) ;
  /* get down reads */
  for (i = 0, j = 0 ; i < arrayMax(aa) ; i += 5)
    { u = arrp(aa, i, BSunit) ;
      if (u[1].i < u[2].i)
	{ t = arrayp(tiles, j++, TILE) ;
	  t->key = u[0].k ; 
	  t->a = u[1].i - 1 ; t->b = u[2].i ;
	  t->t = u[3].i - 1 ; t->e = u[4].i ;
	  t->flag = 0 ;
	}
    }
  if (arrayMax(tiles))
    { arraySort (tiles, tileOrder) ;
      nmvd += doTile(tiles, type, dxp) ;                /* actual work */
      t = arrp(tiles, 0, TILE) - 1;
      i = arrayMax(tiles) ;
      while (t++, i--)
	{ defCptForget (0, t->key) ;
	  bsAddKey (Contig, _Assembled_from, t->key) ;
	  t->a++ ;
	  bsAddData (Contig, _bsRight, _Int, &(t->a)) ;
	  bsAddData (Contig, _bsRight, _Int, &(t->b)) ;
	  t->t++ ;
	  bsAddData (Contig, _bsRight, _Int, &(t->t)) ;
	  bsAddData (Contig, _bsRight, _Int, &(t->e)) ;
	}
    }
  tiles = arrayReCreate (tiles, 300, TILE) ;
  /* get up reads */
  for (i = 0, j = 0 ; i < arrayMax(aa) ; i += 5)
    { u = arrp(aa, i, BSunit) ;
      if (u[1].i > u[2].i)
	{ t = arrayp(tiles, j++, TILE) ;
	  t->key = u[0].k ; 
	  t->a = mx - u[1].i ; t->b = mx - u[2].i + 1 ;
	  t->t = u[3].i - 1 ; t->e = u[4].i ;
	  t->flag = 0 ;
	}
    }
  if (arrayMax(tiles))
    { arraySort (tiles, tileOrder) ;
      nmvd += doTile(tiles, type, dxp) ;                /* actual work */
      t = arrp(tiles, 0, TILE) -  1 ;
      i = arrayMax(tiles) ;
      while (t++, i--)
	{ defCptForget (0, t->key) ;
	  bsAddKey (Contig, _Assembled_from, t->key) ;
	  t->a = mx - t->a ; t->b = mx - t->b + 1 ;
	  bsAddData (Contig, _bsRight, _Int, &(t->a)) ;
	  bsAddData (Contig, _bsRight, _Int, &(t->b)) ;
	  t->t++ ;
	  bsAddData (Contig, _bsRight, _Int, &(t->t)) ;
	  bsAddData (Contig, _bsRight, _Int, &(t->e)) ;
	}
    }
    
 abort:
  bsSave (Contig) ;
  arrayDestroy (aa) ;
  arrayDestroy (tiles) ;

  return nmvd ;
}

/***********************************************************/

int baseCallTileContigs (KEYSET ks, KEY type, int *dxp)
{ KEYSET ks1 = keySetCopy (ks), ks2 ;
  KEY *kp ; int i, j = 0 ;

  i = keySetMax (ks1) ; kp = arrp(ks1, 0, KEY) - 1 ;
  while (kp++, i--)
    if (class(*kp) == _VDNA)
      dnaReClass (*kp, kp) ;

  ks2 = query (ks1, "{Assembled_from} $| {>Subsequence}") ;

  i = keySetMax (ks2) ; 
  if (i)
    { kp = arrp(ks2, 0, KEY) - 1 ;
      while (kp++, i--)
        j += baseCallTileContig (*kp, type, dxp) ;
    }
  i = keySetMax (ks2) ;
  keySetDestroy (ks1) ;
  keySetDestroy (ks2) ;

/*  return i ; */
  return j ;
}

/***********************************************************/
/***********************************************************/

BOOL baseCallPatchContig (KEY contig, int *nnp)
{ OBJ obj = bsUpdate(contig) ;
  Array aa, dna = 0 ;
  int n1, nn = 0, i, amax, clipEnd, x2, dx = 0 ;
  BOOL noAbort = TRUE ;
  KEY dnaKey ;
  LANE *lane ;
  BSunit *uu ;

  sessionGainWriteAccess() ;
  
  aa = arrayCreate (100, BSunit) ;
  if (!obj ||
      !bsGetArray (obj, _Assembled_from, aa, 5) ||
      !bsGetKey (obj, _DNA, &dnaKey) ||
      !(dna = dnaGet (dnaKey)))
    { bsDestroy (obj) ;
      arrayDestroy (aa) ;
      return TRUE ;
    }
  bsSave (obj) ;
 
  messStatus ("Autoedit") ;
  lane = (LANE*) messalloc (sizeof (struct LaneStruct)) ;
  
  for (i = 0 ; i < arrayMax(aa) ; i += 5)
    { if (messIsInterruptCalled())
	{ noAbort = FALSE ;
	  goto abort ;
	}
      uu = arrp (aa, i, BSunit) ;
      lane->key = uu[0].k ;
      lane->x1 = uu[1].i - 1 ;
      lane->x2 = uu[2].i ;
      if (lane->x1 <= lane->x2) 
	{ lane->upSequence = FALSE ; }
      else 
	continue ; /* voir plus bas */
      lane->x3 = lane->x2 ;
      lane->clipTop = uu[3].i - 1 ;
      lane->clipEnd = uu[4].i ;

      if (traceGetLane (0, lane))
	{ x2 = lane->x2 ;
	  clipEnd = lane->clipEnd ;
	  baseCallFindClips (lane,FALSE) ;
	  baseCallUnclipLane (lane, 'F', &dx) ; /* so i patch a longer stretch */
	  lane->x2-- ;
	  if (lane->clipEnd > arrayMax (lane->dna))
	    lane->clipEnd = arrayMax (lane->dna) ;
	  dnaAlignRecale(dna, &(lane->x1), &(lane->x2),
		     lane->dna, lane->clipTop, lane->clipEnd - 1 ) ;
	  lane->x2++ ;
	  n1 = baseCallPatchLane (dna, 0, lane) ;
	  lane->x2 = x2 ;
	  lane->clipEnd = clipEnd ;
	  if (n1)
	    { nn += n1 ;
	      findXclipping (lane) ; /* to adjust the clipping */
	      laneEditSave (0, lane, 0, 0) ;
	    }
	}
      laneDestroy (lane) ;
    }
  
  reverseComplement (dna) ;
  amax = arrayMax (dna) ;
  for (i = 0 ; i < arrayMax(aa) ; i += 5)
    { if (messIsInterruptCalled())
	{ noAbort = FALSE ;
	  goto abort ;
	}

      uu = arrp(aa, i, BSunit) ;
      lane->key = uu[0].k ;
      lane->x1 = uu[1].i - 1 ;
      lane->x2 = uu[2].i ;
      if (lane->x1 > lane->x2) 
	{ lane->upSequence = TRUE ; lane->x2 -= 2 ; }
      else 
	continue ; /* voir plus haut */
      lane->x3 = lane->x2 ;
      lane->clipTop = uu[3].i - 1 ;
      lane->clipEnd = uu[4].i ;
      lane->x1 = amax - lane->x1 - 1 ; 
      lane->x2 = amax - lane->x2 - 1 ; 
      lane->x3 = amax - lane->x3 - 1 ; 

      if (traceGetLane (0, lane))
	{ x2 = lane->x2 ;
	  clipEnd = lane->clipEnd ;
	  baseCallUnclipLane (lane, 'F', &dx) ;
	  lane->x2-- ;
	  dnaAlignRecale(dna, &(lane->x1), &(lane->x2),
		     lane->dna, lane->clipTop, lane->clipEnd - 1) ;
	  lane->x2++ ;
	  n1 = baseCallPatchLane (dna, 0, lane) ;
	  lane->x2 = x2 ;
	  lane->clipEnd = clipEnd ; /* faux du nombre d'insertion */
	  lane->x1 = amax - lane->x1 - 1 ; 
	  lane->x2 = amax - lane->x2 - 1 ; 
	  lane->x3 = amax - lane->x3 - 1 ; 
	  if (n1)
	    { nn += n1 ;
	      findXclipping (lane) ; /* to adjust the clipping */
	      laneEditSave (0, lane, 0, 0) ;
	    }
	}
      laneDestroy (lane) ;
    }
 
 abort: 
  arrayDestroy (aa) ;
  arrayDestroy (dna) ;

  messfree (lane) ;
/*dnaAlignFixContig (contig) ; */
  *nnp += nn ;
  return noAbort ; /* FALSE if F4 interupt, to break higher loops */
}

/*********************************************************/
/*********************************************************/
/*********************************************************/
  /* recognise multiplets of reads if their name only differ by endding */
 /* compare ks to the whole set of other reads */
void baseCallMakeSubclones (KEYSET ks)
{ KEY *kp, *kp1, seq, subclone ;
  OBJ obj = 0 ;
  int n, ii, i1 ;
  char *cp, buf[1000], buf1[1000] ;
  KEYSET kA = 0, ks1 = 0 ;

/* reanalyse the whole set of reads */

  ks1 = /* ks ? keySetCopy (ks) : */ query(0,"Find Read") ;
  kA = keySetAlphaHeap(ks1, keySetMax(ks1)) ;
  keySetDestroy (ks1) ;

  for (ii = 0 , kp = arrp (kA, 0, KEY) ; ii < arrayMax (kA) ;  kp++, ii++)
    { strcpy (buf, name(*kp)) ;
      cp = buf + strlen(buf) ;
      while (cp-- > buf && *cp != '.') ;
      if (cp > buf && *cp == '.' &&
	  strlen(cp + 1) &&
	  strlen(cp) < 12
	  )
	{ *cp = 0 ;
	  for (n = 1, kp1 = kp + 1 , i1 = ii + 1 ; i1 < keySetMax (kA) ; kp1++ , i1++)
	    { strcpy (buf1, name(*kp1)) ;
	      cp = buf1 + strlen(buf1) ;
	      while (cp-- > buf1 && *cp != '.') ;
	      if (cp > buf1 && *cp == '.' &&
		  strlen(cp + 1) &&
		  strlen(cp) < 12 &&
		  (*cp = 0, !strcmp (buf, buf1))
		  ) n++ ;
	      else
		break ;
	    }
	  if (n > 1) /* more than one in this group */
	    { lexaddkey (buf, &subclone, _VClone) ;
	      kp-- ; ii-- ; /* keep ii and kp synchrone */
	      while (ii++, kp++, n--)
		if ( dnaReClass(*kp, &seq) &&
		    (obj = bsUpdate (seq)))
		  { bsAddKey (obj, _Subclone, subclone) ;
		    bsSave (obj) ;
		  }
	      kp-- ; ii-- ; /* keep ii and kp synchrone */
	    }
	}
    }
  keySetDestroy (kA) ;
}

/*********************************************************/
/*********************************************************/

void baseCallTest (void)
{
  KEYSET ks = query (0, "Find EST yk*") ;
  OBJ obj ;
  char *cp, ctf [132], scf[132] ;
  int  i = keySetMax (ks) ;
  KEY key, *kp ; Array aa = 0 ;

  cDNAAlignInit () ;
  doNotStore = TRUE ;

  printf ("Found %d est\n", i) ;
  if (!i) return ;
  kp = arrp (ks, 0, KEY) - 1 ;
  while (kp++, i--)
    {
      key = *kp ;
      
      *ctf = 0 ; *scf = 0 ;

      if ((obj = bsCreate(key)))
	{
	  if (bsGetData (obj, _CTF_File, _Text, &cp))
	    strncpy (ctf, cp, 100) ;
	  if (bsGetData (obj, _SCF_File, _Text, &cp))
	    strncpy (scf, cp, 100) ;
	  bsDestroy (obj) ;
	}
      if (!*ctf && !*scf)
	{ printf ("NODATAFILE %s\n", name(key)) ; continue ; }
      if (!*ctf)
	{ printf ("NOCTFFILE %s %s\n", name(key), scf) ; continue ; }

      aa = baseCallBase (key) ;
      if (aa && arrayMax(aa) > 50)
	printf ("OK %s %d %s %s\n", name(key), arrayMax(aa), ctf, scf) ;
      else if (aa && arrayMax(aa) <= 50)
	printf ("NOLENGTH %s %d %s %s\n", name(key), arrayMax(aa), ctf, scf) ;
      else if (!aa)
	printf ("NOARRAY %s %s %s\n", name(key), ctf, scf) ;
      arrayDestroy (aa) ;
    }
}

/*********************************************************/
/*********************************************************/
/* search for the last base left of x0 */
static BOOL baseCallLocalQuality (LANE *lane, int baseShort)
{
  int i, j0, jmax, x0 = arr (lane->basePos, baseShort, short) ; 
  BASECALL *bb ;
  Array bc = lane->baseCall ;
  int isGood = 0 ;
  
  j0 = 0 ; jmax = arrayMax (bc) ;
  bb = arrp (bc, j0, BASECALL) ;
  while (j0 < jmax && bb->x < x0)
    { j0 += 20 ; bb += 20 ; }
  if (j0 > 20)
    { bb -= 20 ; j0 -= 20 ; }
  while (j0 < jmax && bb->x < x0)
    { j0 ++ ; bb++ ;}
  if (j0 > 0 && bb->x > x0 )
    { bb-- ; j0-- ;}
  for (i = 0 ; i < 5 && i + j0 >= 0 && i + j0 < jmax ; i++)
    {
      if ((bb+i)->flag & BC_LOW || (bb-i)->flag & BC_LOW) 
	break ;
      else
	isGood++ ;
    }
  return isGood > 4 ? TRUE : FALSE ;
} /* baseCallLocalQuality */

/*********************************************************/

static int baseCallFlagOneRnaEditing (KEY tg, LANE *lane, Array dna, Array dnaR, int a1, int a2)
{
  int i, ii, jj, nag = 0, b1, b2, c1, c2, x1, x2, eMax ;
  OBJ Tg = bsCreate (tg) ;
  Array dnaLong, dnaEst = 0, units = arrayCreate (500, BSunit), errors = arrayCreate (64, A_ERR) ;
  BSunit *uu ;
  BOOL isDown ;
  KEY est = lane->key ;
  A_ERR *ep ;

  if (!Tg)
    goto abort ;
  
  if (Tg && bsGetArray (Tg, _Assembled_from, units, 5))
    for (ii = jj = 0 ; ii < arrayMax(units) ; ii += 5)
      {
	uu = arrp(units, ii, BSunit) ;
	if (uu[2].k != est)
	  continue ;

	b1 = uu[0].i ; b2 = uu[1].i ; 
	x1 = uu[3].i ; x2 = uu[4].i ;

	dnaEst = dnaGet (est) ;
	if (a1 < a2) 
	  { 
	    c1 = a1 + b1 - 1 ; c2 = a1 + b2 - 1 ; 
	    isDown = x1 < x2 ? TRUE : FALSE ;
	    dnaLong = dna ;
	    errors = (Array) aceDnaDoubleTrackErrors (dnaEst, &x1, &x2, isDown,
						  dna, dnaR, &c1, &c2, 0, errors, 8, 0, FALSE, 0) ;
	  }
	else 
	  { 
	    c1 = a1 - b1 + 1 ; c2 = a1 - b2 + 1 ;
	    c1 = arrayMax(dna) - c1 + 1 ;
	    c2 = arrayMax(dna) - c2 + 1 ;
	    isDown = x1 < x2 ? TRUE : FALSE ;
	    dnaLong = dnaR ;
	    errors = (Array) aceDnaDoubleTrackErrors (dnaEst, &x1, &x2, isDown,
						      dnaR, dna, &c1, &c2, 0, errors, 8, 0, FALSE, 0) ;
	  }

	if (errors)
	  for (i = 0, ep = arrp (errors, i, A_ERR) ; i < arrayMax(errors) ; ep++, i++)
	    {
	      eMax = arrayMax(errors) ;
	      if (ep->type == ERREUR &&
		  arr (dnaLong, ep->iLong - 1, char) == A_ &&
		  ep->baseShort == G_ &&
		  ep->iShort > 30 && ep->iShort < 500 &&
		  (i == 0 || ep->iLong > (ep - 1)->iLong + 12) &&
		  (i == eMax - 1 || ep->iLong < (ep + 1)->iLong - 12) &&
		  baseCallLocalQuality (lane, ep->iShort) 
		  )
		{
		  OBJ Est = 0 ;
		  int x ;

		  Est = bsUpdate (est) ;
		  arr (dnaEst, ep->iShort - 1, char) = isDown ? R_ : Y_ ;
		  arr (lane->base, ep->iShort - 1, char) = BC_TAG | (isDown ? R_ : Y_) ;
		  x = ep->iLong ;
		  if (dnaLong == dnaR) x = arrayMax(dna) - x - 1 ;
		  if (a1 < a2) x = x - a1 + 1 ;
		  else { x = a1 - x + 1 ; } ;
		  if (Est)
		    {
		      if (bsIsTagInObj (Est, est, str2tag("RNA_editing")))
			{
			  bsAddData (Est, str2tag("RNA_editing"), _Int, &(ep->iShort)) ;
			  bsAddKey (Est, _bsRight, tg) ;
			  bsAddData (Est, _bsRight, _Int, &x) ;
			  bsSave (Est) ;
			}
		      nag++ ;
		    }
		}
	      else
		{ /* for debugging 
		  cc1 = A_ ;
		  cc2 = G_ ;
		  cc1 = arr (dnaLong, ep->iLong - 1, char) ;
		  cc2 = ep->baseShort ;
		  isGood = baseCallLocalQuality (lane, ep->iShort)  ;
		  ********/
		}
	    }
	if (nag)
	  dnaStoreDestroy (est, dnaEst) ;
	else
	  arrayDestroy (dnaEst) ;	  
      }
 abort:
  bsDestroy (Tg) ;
  arrayDestroy (errors) ;
  arrayDestroy (units) ;
  return nag ;
}

/*********************************************************/

BOOL baseCallFlagRnaEditing (KEY tg, int *nreadp, int *nagp, int *nagrp)
{
  KEYSET reads = queryKey (tg, ">Read") ;
  int ii, x, i, nag = 0, nagr = 0, nread = 0, a1 = -1, a2 = -1 ;
  LANE *lane = (LANE*) messalloc (sizeof (struct LaneStruct)) ; 
  KEY read, cosmid ;
  Array dna = 0, dnaR = 0 ;
  Array units = arrayCreate (500, BSunit) ;
  BSunit *uu ;
  OBJ Tg = bsCreate (tg) ;

  cosmid = keyGetKey (tg, str2tag("Genomic_sequence")) ;
  if (Tg && bsGetArray (Tg, str2tag("Covers"), units, 5))
    for (ii = 0 ; ii < arrayMax(units) ; ii += 5)
      {
	uu = arrp(units, ii, BSunit) ;
	a1 = uu[2].i ; a2 = uu[3].i ;
      }
  bsDestroy (Tg) ;

  dna = dnaGet (cosmid) ;
  if (!dna || a1 < 0 || a2 < 0 || a1 >= arrayMax(dna) || a2 >= arrayMax(dna))
    goto abort ;
  dnaR = dnaCopy (dna) ;
  reverseComplement (dnaR) ;

  for (i = 0; i < keySetMax(reads) ; i++)
    {
      read = keySet (reads, i) ;
      lane->key = read ;
      x = 0 ;
      
      if (baseCallGet (lane) &&
	  findBaseCall (lane) &&
	  arrayMax(lane->baseCall))
	{
	  nread++ ;
	  x = baseCallFlagOneRnaEditing (tg, lane, dna, dnaR, a1, a2) ;
	  if (x) 
	    {
	      nagr++ ;
	      if (1) baseCallStore (lane) ;
	    }
	  nag += x ;
	} 
      laneDestroy (lane) ;
    }

  if (nreadp) *nreadp += nread ;
  if (nagp) *nagp += nag ;
  if (nagrp) *nagrp += nagr ;

 abort:
  messfree (lane) ;
  keySetDestroy (reads) ;
  arrayDestroy (dna) ;
  arrayDestroy (dnaR) ;
  arrayDestroy (units) ;
  
  return nag > 0 ;
}

/*********************************************************/
/*********************************************************/

