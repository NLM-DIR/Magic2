/*  File: dnasubs.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 **  Packs an upacks dna arrays.                             
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 27 22:14 1998 (rd)
 * * Apr 23 17:14 1995 (rd): dnaGet() now gets from a Sequence,
 		recursively finding the DNA - complex code from fmap.
 * * Oct 23 20:16 1991 (mieg): Change + to n in decodeChar
 * Created: Wed Oct 23 18:10:21 1991 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: dnasubs.c,v 1.31 2017/02/15 20:36:25 mieg Exp $ */

#define CHRONO

#include "acedb.h"
#include "a.h"
#include "bs.h"
#include "dna.h"
#include "cdna.h"
#include "lex.h"
#include "../whooks/sysclass.h"
#include "../whooks/classes.h"
#include "../whooks/systags.h"
#include "../whooks/tags.h"
#include "pick.h"
#include "freeout.h"
#include "java.h"
#include "bindex.h"
#include "dump.h"
#include "chrono.h"
#include "acedna.h"

static Array dnaUnpackArray(Array pack, AC_HANDLE h) ;
static int dnaDoDumpFastAKeySet (KEYSET kSet, FILE *fil, Stack s, BOOL noMissmatch) ;
static int dnaDoAdd (Array a, KEY seq, int start, int stop, BOOL noMismatch, BOOL reportErrors) ;
static BOOL dnaStoreDestroy2 (KEY key, Array dna, BOOL isPure) ;

#ifdef KS_ADDED
  KEYSET ksAdded = 0 ;
#endif 


/* _VSequence is declared in classes.h/class.c */
static BOOL dnaInitialise (void)
{ 
  KEY key ;
  if (!lexword2key("?Sequence", &key, _VModel)  &&
      iskey(key) != 2)
    { 
      messerror("No model for class Sequence, please edit wspec/models.wrm") ;
      return FALSE ;
    }
  return TRUE ;
}

/************************************/

BOOL dnaDump (FILE* f, Stack buf, KEY key) 
{ 
  Array unpack, pack = arrayGet(key, char, "c") ;
  int level = 0 ;

  if (!pack)
    return FALSE ;
  unpack = dnaUnpackArray(pack, 0) ; 

  if (f)
    level = freeOutSetFile (f) ;
  else if (buf)
    level = freeOutSetStack (buf) ;
  dnaDoDump (unpack, 0, arrayMax(unpack), 9) ; /* can we remove the 9 spaces? */
  if (level)
    freeOutClose(level) ;

  if (pack != unpack)
    arrayDestroy(pack) ;
  arrayDestroy(unpack) ;
  
  return TRUE ;
}

/************************************/

BOOL javaDumpDNA(KEY key) 
{ Array dna, pack = arrayGet(key, char, "c") ;

  if (!pack)
    return FALSE ;

  freeOutf("?DNA?%s?\t?dna?",freejavaprotect(name(key)));

  dna = dnaUnpackArray(pack, 0) ;  /* dna may be packed on disk */
  dnaDecodeArray (dna) ;        /* decode from A_, T_ ... into ascii */
  array(dna,arrayMax(dna),char) = 0 ; /* ensure 0 terminated */
  --arrayMax(dna) ;
  freeOut (arrp(dna,0,char)) ;
  freeOut("?\n\n");

  if (pack !=dna) 
    arrayDestroy(pack);
  arrayDestroy(dna) ;
  return TRUE;
}

/**********************************************************/
   /* called also from dnacptfastaDump */

static BOOL dnaDumpCStyle (KEY key, Array dna, int from, int to, char style)
{
  int ii = to - from + 1 ;
  char buf[2] ;
  char *cp = name(key) ;

  if (!dna || from < 0 || from >= to || to >= arrayMax(dna))
    return FALSE ;
  if (style == 'C')
    dnaDecodeArray (dna) ;

  buf[0] = 'N' ; buf[1] = class(key) ;
  freeOutBinary (buf, 2) ;
  freeOutBinary (cp, strlen(cp) + 1) ; 

  buf[0] = '\n' ;
  buf[1] = 'c' ;
  freeOutBinary (buf,2) ; 
  
  freeOutBinary ((char*) &ii, 4) ; 

  buf[0] = '\n' ;
  buf[1] = 'D' ;
  freeOutBinary (buf,2) ; 

  freeOutBinary (arrp (dna, from, char), ii) ;

  buf[0] = '\n' ; 
  buf[1] = '#' ; 
  freeOutBinary (buf,2) ; 
  
  return TRUE ;
}

/************************************/

BOOL dnaDumpKeyCstyle (KEY key) 
{ 
  Array dna, pack = arrayGet(key, char, "c") ;
  int ii ;

  if (!pack)
    return FALSE ;

  dna = dnaUnpackArray(pack, 0) ;  /* dna may be packed on disk */
  array(dna,arrayMax(dna),char) = 0 ; /* ensure 0 terminated */
  --arrayMax(dna) ;

  ii = arrayMax(dna) ;
  dnaDumpCStyle (key, dna, 0, ii - 1, 'C') ;

  if (pack !=dna) 
    arrayDestroy(pack);
  arrayDestroy(dna) ;

  return TRUE;
}

/************************************/

BOOL dnaParse (int level, KEY key)
{
  char c = 0 ;
  unsigned char c1 = 0 ;
  Array dna = 0 ;
  KEY seqKey = 0, dnaKey = 0 ;
  BOOL ok = FALSE, isPure = FALSE ; /* safe bet */
  char *seqName = 0, *comment = 0 ;
  AC_HANDLE h = ac_new_handle () ;

  dna = dnaParseLevel (level, &c1, &seqName, &comment, h) ;

  /* check name at end of parsing, this allows to see more errors at once */
  if (!dnaReClass (key, &seqKey) || !dnaSubClass (key, &dnaKey))
    {
      messerror ("  dnaParse error at line %d in %s:%s : "
	     "No tag DNA in the model of this class",
	     freestreamline(level), className(key), name(key), c, c) ;
    }

  else if (c1 != 0xff) /* succesful parsing */
    {
      if (dna && arrayMax(dna))
	{
	  Array dna2 = dna ; /* do not copy, just rename */
	  dna = 0 ;   /* to avoid double destruction at end of routine */
	  if (dnaStoreDestroy2 (seqKey, dna2, isPure)) /* updates the sequence obj */
	    ok = TRUE ;
	  else
	    messerror (" failed to store dna %s", name(key)) ;	
	}
      else				/* an empty object is OK */
	ok = TRUE ;
    }

  /* failure */
  if (c1 == 0xff)
    messerror (
	    "  DNAparse error at line %7d in %.25s : bad char %d:%c\n", 
	    freestreamline(level), name(key), (int)c,  c ? c :'-') ; /* avoid printing char zero */

  arrayDestroy(dna) ;
  ac_free (h) ;
  return ok ;
}

/************************************/

BOOL dnaSubClass (KEY seq, KEY *dnaKeyp)   /* gets the parent class of dna, by default a Sequence */
{ 
  KEY dna = 0 ;
  if (class(seq) == _VDNA)
    dna = seq ; 
  else if (class(seq) == _VSequence)
    lexaddkey (name(seq), &dna, _VDNA) ;
  else if (bsIsTagInClass (class (seq), _DNA))
    lexaddkey(messprintf("%s:%s", className(seq), name(seq)), &dna, _VDNA) ;
  if (dna)
    { *dnaKeyp = dna ; return TRUE ; }
  return FALSE ;
}

BOOL dnaReClass (KEY dna, KEY *seqp)   /* gets the parent class of dna, by default a Sequence */
{
  KEY seq = 0 ;
  
  if (class(dna) == _VDNA && strlen(name(dna)))
    {
      char *cp, *buf = strnew (name(dna), 0) ;
      cp = buf ;
      for (cp=buf; *cp && *cp!= ':'; cp++) ;
      if (*cp==':')
	{
	  int cl = 0 ;
	  *cp = 0 ;
	  cl = pickWord2Class (buf) ;
	  if (cl && bsIsTagInClass (cl, _DNA) && strlen (cp+1))
	    lexaddkey (cp+1, &seq, cl) ;
	}
      if (!dnaInitialise())
	return FALSE ;
      if (!seq)
	lexaddkey (name(dna), &seq, _VSequence) ;
      messfree (buf) ;
    }
  else if (bsIsTagInClass (class(dna), _DNA))
    seq = dna ;
  if (seq)
    { *seqp = seq ; return TRUE ; }
  return FALSE ;
}

BOOL dnaInClass (int classe)   /* does this class supports dna ? */
{
  if (classe > 0 &&
      classe < 256 &&
      (classe == _VDNA || bsIsTagInClass (classe, _DNA)))
    return TRUE ;
  return FALSE ;
}

static BOOL seqDnaReClass (KEY seq, KEY *dnap)   /* get a dna with the same name if it exists */
{
  int nn = strlen(name(seq)) ;
  char buf [nn+100] ;

  *dnap = 0 ;
  if (! dnaInClass (class(seq)))
    return FALSE ;
  else if (class(seq) == _VSequence && strlen(name(seq)) &&
      lexword2key (name(seq), dnap, _VDNA))
    return TRUE ;
  else
    {
      sprintf (buf, "%s:%s", className (seq), name(seq)) ;
      if (lexword2key (buf, dnap, _VDNA))
	return TRUE ;
    }
  return FALSE ;
}

/*************************************/
/*************************************/
#define MAGIC_PACK 7 
#define MAGIC_PACK_ODD 8 
#define MAGIC_DOUBLE_PACK 6
#define MAGIC_NEW_DOUBLE_PACK 4

static int dnaPackingType(Array dna)
{ 
  int n = arrayMax(dna) ;
  char *cp = arrp(dna,0,char) ;

  if(n<4)   /* no packing */
    return 0 ;
  
  chrono ("dnaPackingType") ;

  while(n--)
    switch(*cp++)
      { 
      case A_: case T_: case G_: case C_:
	break ;
      default:   /* at least one ambiguous base */
	return 1 ;
      }

  chronoReturn () ;
  return 2 ;  /* double packing, only a t g c present */
}

/********************/

static void modifyLength (KEY key, KEY seq0, int n)
{ KEY seq = seq0;
  OBJ obj ;
  int n1 ;

  /* this addition is written for the sake of the server client
     it amounts to prevent XREFing while parsing from server */
  if (keyFindTag(seq,_DNA) && (obj = bsCreate (seq)))
    { if (bsFindKey (obj, _DNA, key) &&
	  bsGetData (obj, _bsRight, _Int, &n1) &&
	  n == n1)
      { bsDestroy (obj) ; return ; }
      bsDestroy (obj) ;
    }
  if ((obj = bsUpdate (seq)))
    { bsAddKey (obj, _DNA, key) ;
      bsAddData (obj, _bsRight, _Int, &n) ;
      bsSave (obj) ;
    }
}

/********************/
/* isPure is an optimisation if we know we just use atgc */
static BOOL dnaStoreDestroy2 (KEY key, Array dna, BOOL isPure)
{ 
  char *cp , *cq , c1, *base ;
  int m = dna ? arrayMax(dna) : 0 ;
  char dbp[16] ;
  KEY seqKey = 0, dnaKey = 0 ;

  if (! dnaReClass (key, &seqKey) ||
      ! dnaSubClass (key, &dnaKey))
    {   
      messerror ("dnaStoreDestroy called on a non-dna key: %s:%s", className(key), name(key)) ; 
      arrayDestroy (dna) ; 
      return FALSE ;
    }
 
  chrono ("dnaStoreDestroy2") ;
  if (dna) modifyLength (dnaKey, seqKey, arrayMax(dna)) ;
  dbp[A_] = 0 ;  dbp[G_] = 1 ;  dbp[C_] = 2 ;  dbp[T_] = 3 ;

  if (arrayMax (dna))
    switch(isPure ? 2 : dnaPackingType(dna) )
      {
      case 0:   /* no packing , no coding */
	dnaDecodeArray(dna) ;
	break ;
	
      case 1:  /* 2 bases per byte */
	c1  = (array(dna,0,char)  << 4) | array(dna,1,char) ;
	array(dna,0,char) = m%2 ? MAGIC_PACK_ODD : MAGIC_PACK ;
	array(dna,1,char) = c1 ;  /* first 2 bases */
	
	/* all the rest but possibly one */
	base =  arrp(dna,0,char) ;
	cp = cq =  arrp(dna,2,char) ;
	m -= 2 ;
	while(m>1)
	  {
	    *cp++ = ((*cq) << 4 ) | *(cq + 1) ;
	    cq += 2 ;
	    m -= 2 ;
	  }
	if(m)              /* last base in odd case */
	  {
	    *cp++ = *cq << 4 ;
	  }
	arrayMax(dna) = cp - base ; 
	break ;
	
      case 2:  /* 4 bases per byte */
	cq = arrp(dna,0,char) ;
	c1 =
	  (dbp[((int) *cq) & 0xff] << 6 ) |
	  (dbp[((int) *(cq+1)) & 0xff] << 4 ) |
	  (dbp[((int) *(cq+2)) & 0xff] << 2 ) |
	  dbp[(int) *(cq+3)] ;
	array(dna,0,char) =  MAGIC_NEW_DOUBLE_PACK ;
	array(dna,1,char) = m%4 ;
	array(dna,2,char) = c1 ;  /* first 4 bases */
	
	/* all the rest but possibly 3 */
	base =  arrp(dna,0,char) ;
	cp =  arrp(dna,3,char) ;
	cq =  arrp(dna,4,char) ;
	m -= 4 ;
	while(m>3)
	  {
	    *cp++ =
	      (dbp[((int) *cq) & 0xff] << 6 ) |
	      (dbp[((int) *(cq+1)) & 0xff] << 4 ) |
	      (dbp[((int) *(cq+2)) & 0xff] << 2 ) |
	      dbp[((int) *(cq+3)) & 0xff] ;
	    cq += 4 ;
	    m -= 4 ;
	  }
	
	if(m--)              /* last 3 bases */
	  { base-- ; /* to fix arrayMax, without using cp++ */
	  *cp = (dbp[((int) *cq++) & 0xff] << 6 ) ;
	  if(m--)
	    { *cp |= (dbp[((int) *cq++) & 0xff] << 4 ) ;
	    if(m--)
	      *cp |= (dbp[((int) *cq++) & 0xff] << 2) ;
	    }
	  }
	arrayMax(dna) = cp - base ; 
	break ;
      }
  
  arrayStore(dnaKey, dna, "c") ;
   arrayDestroy(dna) ;  
   /* because the array has been mangled by the compressor
      in a way it would be more clean to compress a copy of the dna array
      but it is most often ok to destroy the array so we win a lot
      of efficiency by letting the user copy himself if he needs to.
   */
   chronoReturn () ;
   return TRUE ;
}

BOOL dnaStoreDestroy (KEY key, Array dna)
{ 
  return dnaStoreDestroy2 (key, dna, FALSE) ; 
}

/********************************************************************/

static Array dnaUnpackArray(Array pack, AC_HANDLE h)
{ Array unpack ;
  char *cp, *cq ;
  static char undoublepack[] = {A_, T_, G_, C_} ; /* RD must be static to initialise on SGI */
  static char newundoublepack[] = {A_, G_, C_, T_} ; /* RD must be static to initialise on SGI */
  int m, n ;

  if(!pack)
    return 0 ;
  cp = arrp(pack,0,char) ;
  m = 0 ;
  switch(*cp)
    {
    case MAGIC_PACK_ODD: /* 2 bases per byte, odd total */
      m = -1 ; /* fall through */
    case MAGIC_PACK: /* MAGIC packed form */
        n = arrayMax(pack) ;
      if(n<=1)
	  return 0 ;

      m += 2*(n-1) ;  /* skip magic, then every char is 2 base except */
	              /* last char may be a single base in ODD case */
      unpack = arrayHandleCreate(m+1,char,h) ;  /* space for the zero */
      array(unpack,m,char) = 0 ; /* zero terminates */
      arrayMax(unpack) = m ; /* implies arrayMax = m */
      cp = arrp(pack,0,char) ;  /* so as to start decoding on byte 1 */
      cq = arrp(unpack,0,char) ;
      while(cp++, m--)
	{ *cq++ = (*cp >> 4 ) & (char)(15) ;  /* first half byte */
	       /* &0xf to ensure a left zero after right shift */
	  if(m--)
	    *cq++ = *cp & (char)(15) ;      /* second half byte */
	  else
	    break ;
	}
      return unpack ;

    case MAGIC_DOUBLE_PACK:  /* 4 bases per byte */
        n = arrayMax(pack) ;
      if(n<=2)           /* first byte is MAGIC,  */
	  return 0 ;

      m = array(pack,1,char) ;  /* second byte is max%4 */
      if (!m) m=4 ;  /* favorable case last byte contains 4 bases not zero */
      m = 4*(n-2) - (4-m);
               /* skip magic, residue, then every char is 4 base except */
	              /* last char which holds residue */
      unpack = arrayHandleCreate(m+1,char, h) ;  /* ensures zero terminated string */
      array(unpack,m-1,char) = 0 ; /* implies arrayMax = m */
      cp = arrp(pack,1,char) ; /* so as to start decoding on byte 2 */
      cq = arrp(unpack,0,char) ;
      while(cp++, m--)
	{ *cq++ = undoublepack[(*cp >> 6 ) & (char)(3)] ;  /* first quarter */
	  if(m--)
	    *cq++ =  undoublepack[(*cp >> 4 ) & (char)(3)] ; 
	  else
	    break ;
	  if(m--)
	    *cq++ =  undoublepack[(*cp >> 2 ) & (char)(3)] ; 
	  else
	    break ;
	  if(m--)
	    *cq++ =  undoublepack[(*cp ) & (char)(3)] ; 
	  else
	    break ;
	}
      return unpack ;

    case MAGIC_NEW_DOUBLE_PACK:  /* 4 bases per byte, new coding, such that ~ is complement */
        n = arrayMax(pack) ;
      if(n<=2)           /* first byte is MAGIC,  */
	  return 0 ;

      m = array(pack,1,char) ;  /* second byte is max%4 */
      if (!m) m=4 ;  /* favorable case last byte contains 4 bases not zero */
      m = 4*(n-2) - (4-m);
               /* skip magic, residue, then every char is 4 base except */
	              /* last char which holds residue */
      unpack = arrayHandleCreate(m+1,char,h) ;  /* ensures zero terminated string */
      array(unpack,m-1,char) = 0 ; /* implies arrayMax = m */
      cp = arrp(pack,1,char) ; /* so as to start decoding on byte 2 */
      cq = arrp(unpack,0,char) ;
      while(cp++, m--)
	{ *cq++ = newundoublepack[(*cp >> 6 ) & (char)(3)] ;  /* first quarter */
	  if(m--)
	    *cq++ =  newundoublepack[(*cp >> 4 ) & (char)(3)] ; 
	  else
	    break ;
	  if(m--)
	    *cq++ =  newundoublepack[(*cp >> 2 ) & (char)(3)] ; 
	  else
	    break ;
	  if(m--)
	    *cq++ =  newundoublepack[(*cp ) & (char)(3)] ; 
	  else
	    break ;
	}
      return unpack ;

    default:    /* uncoded char form, rare I hope */
      dnaEncodeArray(pack) ;
      return pack ;
    }
}

/******************************************************************/
/************** functions to get DNA from a KEY *******************/

static Array localDnaGet (KEY key, AC_HANDLE h)
{
  Array pack, unpack ;

  pack = arrayGet(key, char,"c") ;
  if(!pack)
    return 0 ;
  unpack = dnaUnpackArray(pack, h) ;
  if (pack != unpack)
    arrayDestroy(pack) ;

  return unpack ;
}

/*************************/

typedef struct { int x, y ; } ExonStruct ;
typedef struct { int x, y ; KEY k ; } Exon3Struct ;

static int dnaExonsOrder (const void *a, const void *b)
{
  const ExonStruct *ea = (const ExonStruct *)a ;
  const ExonStruct *eb = (const ExonStruct *)b ;
  return ea->x - eb->x ;
}
#ifdef FOR_DEBUGGING
static void  dnaExonsShow(Array a) 
{
  int i ;
  Exon3Struct *e ;
  for (i=0; i < arrayMax(a) ; i++)
    { 
      e = arrp (a, i, Exon3Struct) ;
      printf("%d %d %s\n",e->x, e->y, name(e->k)) ;
    }
}
#endif
void dnaExonsSort3 (Array a, int nn)
{		/* can't just arraySort, because a is in BSunits, not
		   pairs of BSunits */
  int i, n ;
  Array b = 0 ;

  n = arrayMax(a) / nn ;
  b = arrayCreate (n, Exon3Struct) ;
  for (i = 0 ; i < n ; ++i)
    { arrayp(b,i,Exon3Struct)->x = arr(a,nn*i,BSunit).i ;
      arrp(b,i,Exon3Struct)->y = arr(a,nn*i+1,BSunit).i ;
      if (nn == 3) arrp(b,i,Exon3Struct)->k = arr(a,nn*i+2,BSunit).k ;
    }
  arraySort (b, dnaExonsOrder) ;
  for (i = 0 ; i < n ; ++i)
    { arr(a,nn*i,BSunit).i = arrp(b,i,Exon3Struct)->x ;
      arr(a,nn*i+1,BSunit).i = arrp(b,i,Exon3Struct)->y ;
      if (nn == 3) arr(a,nn*i+2,BSunit).k = arrp(b,i,Exon3Struct)->k ;
    }
  arrayDestroy (b) ;
}

void dnaExonsSort (Array a)
{
  dnaExonsSort3 (a, 2) ;
}
/*************************/

static Array dnaDoGet (KEY key, BOOL noMismatch, AC_HANDLE h)
{ 
  Array a = 0 ;			/* the result */
  OBJ obj = 0 ;
  KEY dna ;
  int len ;

  if (class(key) == _VDNA)
    a = localDnaGet (key, h) ;

  else if (keyFindTag (key, _DNA) &&
	   seqDnaReClass (key, &dna) &&
	   (a = localDnaGet (dna, h))
	   )
    { obj = 0 ; /* to please the compiler */
#ifdef KS_ADDED
      if (!ksAdded) ksAdded = keySetCreate() ;
      keySet(ksAdded,keySetMax(ksAdded)) = key ;
#endif
    }


  else if (( !keyFindTag (key, _DNA) &&
	     !keyFindTag (key, _Source) &&
	     !keyFindTag (key, _Subsequence) &&
	     !keyFindTag (key, str2tag("SMAP"))) ||
	   !(obj = bsCreate(key)))
    return 0 ;
  
  else if (bsGetKey (obj, _DNA, &dna))
    { a = localDnaGet (dna, h) ;
#ifdef KS_ADDED
      if (!ksAdded) ksAdded = keySetCreate() ;
      keySet(ksAdded,keySetMax(ksAdded)) = key ;
#endif
    }

  else				/* nontrivial - must use dnaAdd() */
    { OBJ Seq = bsCreate (key) ;
      KEY seq = key ;
      KEY parent, dummy = 0 ;
      int start = 1, stop = 0, x1, x2; 
      Array s_children = 0 ;

	/* recurse up - done already if from fMapDisplay, but harmless */
      while (Seq)
	{
	  if (bsFindTag (Seq, str2tag("SMAP")) && 
	      !bsFindTag (Seq, str2tag("Splicing")) &&     /* stop the up recursion */
	      bsFindTag (Seq, str2tag("S_Parent")) &&
	      bsGetKeyTags (Seq, _bsRight, &dummy) &&
	      bsGetKey (Seq, _bsRight, &parent)
	      ) /* new tag2 recursive system */
	    { 
	      int i ;
	      BSunit *uu ;
	      /* find parent
	      if (bsFindTag (Seq, str2tag("Splicing")) ||  
		  !bsFindTag (Seq, str2tag("S_Parent")) ||
		  !bsGetKeyTags (Seq, _bsRight, &dummy) ||
		  !bsGetKey (Seq, _bsRight, &parent))
		break ;  stop the recursion */
	      bsDestroy (Seq) ;
	      /* find the child in the parent */
	      s_children = arrayCreate (128, BSunit) ;
	      if ((Seq = bsCreate(parent)) &&
		  bsGetArray (Seq, str2tag("S_Child"), s_children, 4))
		{
		  for (uu = 0, i = 0 ; i < arrayMax(s_children) ; i += 4)
		    {
		      uu = arrp(s_children, i, BSunit) ;
		      if (uu[1].k == seq) break ;
		      uu = 0 ;
		    }
		  if (uu)
		    {
		      x1 = uu[2].i ; x2 = uu[3].i ; 
		      /* finalize */
		      if (!stop)
			{ start = x1 ; stop = x2 ; }
		      else
			if (x1 < x2)
			  { start += x1-1 ; stop += x1-1 ; }
			else
			  { start = x1+1 - start ; stop = x1+1 - stop ; }
		      seq = parent ;
		      continue ;
		    }
		}
	    }

	  if (bsGetKey (Seq, _Source, &parent))
	    { 
	      bsDestroy (Seq) ; 
	      /* find the child in the parent */
	      if (!(Seq = bsCreate(parent)) ||
		  !bsFindKey (Seq, _Subsequence, seq) ||
		  !bsGetData (Seq, _bsRight, _Int, &x1) ||
		  !bsGetData (Seq, _bsRight, _Int, &x2))
		break ;/* stop the recursion */
	      /* finalize */
	      if (!stop)
		{ start = x1 ; stop = x2 ; }
	      else
		if (x1 < x2)
		  { start += x1-1 ; stop += x1-1 ; }
		else
		  { start = x1+1 - start ; stop = x1+1 - stop ; }
	      seq = parent ;
	      continue ;
	    }
	  else
	    break ;
	}

      bsDestroy (Seq) ;		/* OK if 0 */
      arrayDestroy (s_children) ;

      if (!stop)   /* no parent - find length from children */
	{ 
	  Array u = arrayCreate (12, BSunit) ;
	  int i ;
	  
	  stop = 0 ; start = 1 ;
	  if (bsGetArray (obj,str2tag("S_Child"), u, 4))
	    for (i = 0 ; i < arrayMax(u) ; i += 4)
	      { if (arr(u,i+2,BSunit).i > stop)
		  stop = arr(u,i+2,BSunit).i ;
		if (arr(u,i+3,BSunit).i > stop)
		  stop = arr(u,i+3,BSunit).i ;
                if (arr(u,i+2,BSunit).i < start)
		  start = arr(u,i+2,BSunit).i ;
		if (arr(u,i+3,BSunit).i < start)
		  start = arr(u,i+3,BSunit).i ;
	      }
	  if (bsGetArray (obj,_Subsequence, u, 3))
	    for (i = 0 ; i < arrayMax(u) ; i += 3)
	      { if (arr(u,i+1,BSunit).i > stop)
		  stop = arr(u,i+1,BSunit).i ;
		if (arr(u,i+2,BSunit).i > stop)
		  stop = arr(u,i+2,BSunit).i ;
                if (arr(u,i+1,BSunit).i < start)
		  start = arr(u,i+1,BSunit).i ;
		if (arr(u,i+2,BSunit).i < start)
		  start = arr(u,i+2,BSunit).i ;
	      }
	  arrayDestroy (u) ;
	  if (start < 1 || stop < start)  /* mieg */
	    { bsDestroy (obj) ;
	      return 0 ;
	    }
	  len = stop - start + 1 ; /* mieg */
	}
      else if (start < stop)
	len = stop - start + 1 ;
      else
	len = start - stop + 1 ;

      a = arrayCreate (len+1, char) ;    /* space for terminal 0 */
      array (a, len, char) = 0 ;   /* add terminal zero and make room */
      arrayMax(a) = len ;  /* restore max */

      /*
	if (dnaAdd2 (a, seq, 0, start, stop, noMismatch)  < 2)
	arrayDestroy(a) ;
      */
      if (dnaDoAdd (a, seq, start - 1, stop - 1, noMismatch, FALSE)  < 2)
	arrayDestroy(a) ;
      else
	{ 
	  if (bsFindTag (obj, _Source_Exons))
	    { Array b = arrayCreate (len+1, char) ;
	      Array u = arrayCreate (12, BSunit) ;
	      int i, j , mx = 0 , x, lastx = 0 ;
	      char *cp, *cq ;

	      bsFlatten (obj, 2, u) ;
	      dnaExonsSort (u) ;
	      for (i = 0 ; i < arrayMax(u) ; i += 2)
		{ 
		  x = arr(u, i, BSunit).i - 1 ;
		  j = arr(u, i+1, BSunit).i - arr(u, i, BSunit).i + 1 ;
		  if (x < 0 || 
		      j <= 0 || x + j  > arrayMax(a) ||
		      lastx > x
		      ) 
		    { arrayDestroy (b) ; break ; } /* mieg: a protection, really we should translate */
		  cq = arrayp(a, x, char) ; /* mieg, may be x >= len */
		  lastx = x + j ;  /* mieg: 2000_09_28 last base copied, never backtrack */
                  x = mx + j - 1 ;
		  array(b, x, char) = 0 ; /* make room, may reallocate */

		  cp = arrp(b, mx, char) ;
		  mx += j ;
		  while (j--)
		    *cp++ = *cq++ ;
		}
	      arrayDestroy (u) ;
	      arrayDestroy (a) ; 	/* replace a by b */
	      a = b ;
	    }
	}
    }

  bsDestroy (obj) ;
  if (!a)
    return 0 ;

	/* add a terminal zero, useful for unix calls, but not 
	   really correct because there can be internal 0's */
  array(a, arrayMax(a), char) = 0 ;
  --arrayMax(a) ;

  return a ;
}

/* mieg, june 2000, 
   Subsequence ?Sequence UNIQUE a1 UNIQUE a2 b1 UNIQUE b2
        // a1 a2 : coordinates of subseq in parent
        // b1 b2 : dna in parent default : a1 a2, must be inside a1 a2
	//  use b1 == -1 to use NO dna from this subsequence
  Exemple
Sequence Chrom
Subsequence c1 11 30  16 25
Subsequence c2 21 40  28 38
Subsequence c3 31 50  -1

>c1
AAAAATTTTTTTTTTCCCCC
>c2
GGGGGGGAAAAAAAAAAGGG
>c3
ATGATGATGGCTACTACTAA
*/

Array dnaGet (KEY key)
{
  return dnaDoGet (key, TRUE, 0) ;
}

Array dnaHandleGet (KEY key, AC_HANDLE h)
{
  return dnaDoGet (key, TRUE, h) ;
}

Array dnaGetWithErrors (KEY key)
{
  return dnaDoGet (key, FALSE, 0) ;
}
/*
  static int dnaAdd2 (Array a, KEY seq, int offset, int sstart, int sstop, BOOL noMismatch)
*/

static int dnaDoAdd (Array a, KEY seq, int start, int stop, BOOL noMismatch, BOOL reportErrors)
{  
  OBJ obj = 0 ;
  Array dna = 0 ;
  Array aa = 0 ;
  KEY subSeq, dnaKey ;
  int pos1, pos2, i ;
  char *cp, *cq, *cq0 ;
  int r1, result = 1 ; /* 0: mismatch, 1: absent,  2: success */
  BOOL isForwards = (stop >= start) ;
 
  if (!seq)
    return 1 ;

				/* seq itself has a DNA sequence */
  if (bIndexFind(seq,_DNA) &&
      (obj = bsCreate (seq)) &&
      bsGetKey (obj, _DNA, &dnaKey) && 
      (dna = localDnaGet (dnaKey, 0)))
    { cq = cq0 = arrayp(a, 0, char) ;

      if (start < 0)		/* first check bounds */
	{ if (!isForwards)
	    { bsDestroy (obj) ;
	      return 1 ;
	    }
	  cq -= start ;	/* remember start is -ve, so this adds */
	  start = 0 ;
	}
      if (start >= arrayMax(dna))
	{ if (isForwards)
	    { bsDestroy (obj) ;
	      return 1 ;
	    }
	  cq += start - arrayMax(dna) + 1 ;
	  start = arrayMax(dna) - 1 ;
	}
      if (stop < 0)
	{ if (isForwards)
	    { bsDestroy (obj) ;
	      return 1 ;
	    }
	  stop = 0 ;
	}
      if (stop >= arrayMax(dna))
	{ if (!isForwards)
	    { bsDestroy (obj) ;
	      return 1 ;
	    }
	  stop = arrayMax(dna) - 1 ;
	}

#ifdef KS_ADDED
      if (!ksAdded) ksAdded = keySetCreate() ;
      keySet(ksAdded,keySetMax(ksAdded)) = seq ;
#endif

      result = 2 ;
      if (isForwards)
	{ cp = arrp(dna, start, char) ;
	  for (i = start ; i <= stop ; ++i)
	    { 
	      if (*cp && *cq && !(*cp & *cq))  /* N_ is compatible with A_ and so on */
		{
		  /* Complain at most once for each bit of dna and only      */
		  /* if caller wants errors reported.                        */
		  if (reportErrors)
		    {
		      if (result)
			messerror ("sequence mismatch at %d in %s", i, name(seq)) ;
		      result = 0 ;
		    }
		}
	      *cq++ |= *cp++ ;
	    }
	}
      else
	{ cp = arrp(dna, start, char) ;
	  for (i = start ; i >= stop ; --i) /* RD added complementBase! */
	    { 
	      if (*cp && *cq && !(complementBase(*cp) & *cq))
		{
		  /* Complain at most once for each bit of dna and only      */
		  /* if caller wants errors reported.                        */
		  if (reportErrors)
		    {
		      if (result)
			messerror ("sequence mismatch at %d in %s", i, name(seq)) ;
		      result = 0 ;
		    }
		}
	      *cq++ |= complementBase(*cp--) ;  /* mix the 2 values */
	    }
	}

      arrayDestroy(dna) ;
      bsDestroy (obj) ;
      if (cq - cq0 > arrayMax (a) || cq < cq0)
	messcrash ("overflow in dnaGet: max = %d, length = %d",
		   arrayMax(a), cq - cq0) ;
      return result ;
    }

				/* reconstruct from subsequences */
  if (!obj &&
      (
       (!bIndexFind(seq,_Subsequence) && !bIndexFind(seq,str2tag("SMAP"))) ||
       !(obj = bsCreate (seq)))
      )

    goto end ;
  aa = arrayCreate (12, BSunit) ;
  if (bsGetArray (obj,str2tag("S_Child"), aa, 4))
    { 
      for (i=0 ; i<0 && i < arrayMax(aa) ; i += 4)
	{ 
	  pos1 = arr(aa,i+2,BSunit).i ;
	  pos2 = arr(aa,i+3,BSunit).i ;
	  subSeq = arr(aa,i+1,BSunit).k ;
	  if (keyFindTag (subSeq, str2tag("Splicing")))
	    continue ;
	  if (pos1 && pos2)
	    { --pos1 ; --pos2 ; r1 = 1 ;
	    if (pos1 < pos2)
	      { if ((start >= pos1 || stop >= pos1) && 
		    (start <= pos2 || stop <= pos2))
		r1 = dnaDoAdd(a, subSeq, start-pos1, stop-pos1, noMismatch, reportErrors);
	      }
	    else
	      { if ((start >= pos2 || stop >= pos2) && 
		    (start <= pos1 || stop <= pos1))
		r1 = dnaDoAdd(a, subSeq, pos1-start, pos1-stop, noMismatch, reportErrors);
	      }
	    switch (r1)
	      {
	      case 0: if (noMismatch) { result = 0 ; goto end ;}
		/* else fall through */
	      case 1: break ; /* sequence unknown */
	      case 2: result = 2 ; break ;  /* success */
	      }
	    }
	  else
	    messerror ("Subsequence coords missing for %s in %s", 
		       name (subSeq), name(seq)) ;
	}
    }
  
  if (bsGetArray (obj,_Subsequence, aa, 3))
    {
      for (i=0 ; i < arrayMax(aa) ; i += 3)
	{
	  pos1 = arr(aa,i+1,BSunit).i ;
	  pos2 = arr(aa,i+2,BSunit).i ;
	  subSeq = arr(aa,i,BSunit).k ;
	  if (pos1 && pos2)
	    { --pos1 ; --pos2 ; r1 = 1 ;
	    if (pos1 < pos2)
	      { if ((start >= pos1 || stop >= pos1) && 
		    (start <= pos2 || stop <= pos2))
		r1 = dnaDoAdd(a, subSeq, start-pos1, stop-pos1, noMismatch, reportErrors);
	      }
	    else
	      { if ((start >= pos2 || stop >= pos2) && 
		    (start <= pos1 || stop <= pos1))
		r1 = dnaDoAdd(a, subSeq, pos1-start, pos1-stop, noMismatch, reportErrors);
	      }
	    switch (r1)
	      {
	      case 0: if (noMismatch) { result = 0 ; goto end ;}
		/* else fall through */
	      case 1: break ; /* sequence unknown */
	      case 2: result = 2 ; break ;  /* success */
	      }
	    }
	  else
	    messerror ("Subsequence coords missing for %s in %s", 
		       name (subSeq), name(seq)) ;
	}
    }
      
 end:
  arrayDestroy (aa) ;
  bsDestroy (obj) ;

  return result ;
}

/**********************************************************/

/* dnaadd is used by the graphic code and knows base zero */
int dnaAdd (Array a, KEY seq, int start, int stop, BOOL noMismatch)
{ 
  BOOL reverse = FALSE ;
  int result ;

  if (start > stop)
    {
      int tmp = stop ; stop = start ; start = tmp ;
      reverse = TRUE ;
    } 
  
                /* dna2 est ecrit en coord bio */
  array (a, stop - start, char) = 0 ;
  /* result = dnaAdd2 (a, seq, offset, start + 1, stop + 1, noMismatch) ;  */
  result = dnaDoAdd (a, seq, start, stop, noMismatch, FALSE) ; 

  if (reverse)
    reverseComplement (a) ;

  return result ;
}

/**********************************************************/
/**********************************************************/

FILE *dnaFileOpen (void)
{
  static char fileName[FIL_BUFFER_SIZE],dirName[DIR_BUFFER_SIZE] ;

  return filqueryopen (dirName, fileName, "dna", "w",
		       "Choose a file for export in fasta format") ;
}

/**********************************************************/

static BOOL dnaDoDumpFastAKey (KEY key, FILE *fil, Stack s, BOOL noMissmatch, char style, int x1, int x2, BOOL noClassNam)
{ 
  Array a = dnaDoGet (key, noMissmatch,0) ;

  if (a)
    { 
      BOOL result = FALSE ;

      if (x1 > arrayMax(a)) 
	x1 = arrayMax (a) ;
      if (x2 > arrayMax(a)) 
	x2 = arrayMax (a) ;
      if (x1 == 0)
	x1 = 1 ;
      if (x2 == 0)
	x2 = arrayMax(a) ;
      if (x1 == -1)
	x1 = arrayMax(a) ;
      if (x2 == -1 || x2 == 0)
	x2 = arrayMax(a) ;
      if (x1 > x2)
	{ reverseComplement (a) ; x1 = arrayMax(a) - x1 + 1 ; x2 = arrayMax(a) - x2 + 1 ; }

      if (x1 != x2)  /* in this case i do not know the orientation, or both were outside 1/max */
	{
	  x1-- ; x2-- ;
	  if (style == 'C' || style == 'B')
	    {
	      dnaDumpCStyle (key, a, x1, x2, style) ;
	    }
	  else
	    {
	      if (class (key) == _VDNA || class (key) == _VSequence)
		{
		  KEY model = keyGetKey (key, str2tag("Model_of")) ;
		  KEY geneid = keyGetKey (key, str2tag("GeneId_pg")) ;
		  
		  if (model && geneid)
		    result = dnaDumpFastA (a, x1, x2,
					   messprintf ("%s|Gene|%s|GeneId|%s", name(key), name(model), name(geneid))
					   , fil, s) ;
		  else if (geneid)
		    result = dnaDumpFastA (a, x1, x2,
					   messprintf ("%s|GeneId|%s", name(key), name(geneid))
					   , fil, s) ;
		  else if (model)
		    result = dnaDumpFastA (a, x1, x2,
					   messprintf ("%s|Gene|%s", name(key),name(model))
					   , fil, s) ;
		  else
		    result = dnaDumpFastA (a, x1, x2,
					   messprintf ("%s",name(key))
					   , fil, s) ;
		}
	      else if (class (key) == _VmRNA)
		{
		  KEY gene = 0, tg = keyGetKey (key, _From_gene) ;
		  char geneid[1000] ;
		  if (! tg) tg = keyGetKey (key, _From_prediction) ;
		  gene = keyGetKey (key, _Gene) ; /* prefer the gene name which links directly to the mRNA */
		  if (tg && ! gene) gene = keyGetKey (tg, _Gene) ;
		  
		  geneid[0] = 0 ;
		  if (gene && keyFindTag (gene, _GeneId))
		    {
		      int i ;
		      char *cp = geneid, *cq, *buf1 = geneid + 900 ;
		      KEYSET ks = queryKey (gene, ">GeneId") ;
		      
		      for (i = 0 ; i < keySetMax(ks) && cp < buf1 ; i++)
			{
			  if (i) { cp-- ; *cp++ = ';' ; }
			  cq = name (keySet(ks,i));
			  while ((*cp++ = *cq++)) ;
 			}
		      keySetDestroy (ks) ;
		    }
		  
		  if (noClassNam)
		    {
		      if (geneid[0])
			result = dnaDumpFastA (a, x1, x2,
					       messprintf ("%s|Gene|%s|GeneId|%s", name(key),name(gene),geneid)
					       , fil, s) ;
		      else if (gene)
			result = dnaDumpFastA (a, x1, x2,
					       messprintf ("%s|Gene|%s", name(key),name(gene))
					       , fil, s) ;
		      else
			result = dnaDumpFastA (a, x1, x2,
					       messprintf ("%s", name(key))
					       , fil, s) ;
		    }
		  else
		    {
		      if (geneid[0])
			result = dnaDumpFastA (a, x1, x2,
					       messprintf ("%s:%s|Gene|%s|GeneId|%s", className(key), name(key),name(gene),geneid)
					       , fil, s) ;
		      else if (gene)
			result = dnaDumpFastA (a, x1, x2,
					       messprintf ("%s:%s|Gene|%s", className(key), name(key),name(gene))
					       , fil, s) ;
		      else
			result = dnaDumpFastA (a, x1, x2,
					       messprintf ("%s:%s", className(key), name(key))
					       , fil, s) ;
		    }
		}
	      else
		{
		  if (noClassNam)
		    {
		      result = dnaDumpFastA (a, x1, x2
					     , name(key)
					     , fil, s) ;
		    }
		  else
		    {
		      result = dnaDumpFastA (a, x1, x2
					     , messprintf ("%s:%s", className(key), name(key))
					     , fil, s) ;
		    }
		}
	    }
	}
      arrayDestroy (a) ;
      return result ;
    }
  else
    return FALSE ;
}

BOOL dnaDumpFastAKey (KEY key, FILE *fil, Stack s, char style)
{
  return dnaDoDumpFastAKey (key, fil, s, TRUE, style, 0, 0, 0) ;
}

BOOL dnaDumpFastAKeyWithErrors (KEY key, FILE *fil, Stack s, char style)
{
  return dnaDoDumpFastAKey (key, fil, s, FALSE, style, 0, 0, 0) ;
}

BOOL dnaZoneDumpFastAKey (KEY key, FILE *fil, Stack s, char style, int x1, int x2, BOOL withMissmatch, BOOL noClassNam)
{
  return dnaDoDumpFastAKey (key, fil, s, ! withMissmatch, style, x1, x2, noClassNam) ;
}

/**********************************************************/

static int dnaDoDumpFastAKeySet (KEYSET kSet, FILE *fil, Stack s, BOOL noMissmatch)
{
  KEYSET alpha ;
  int i, n = 0 ;

  if (!keySetExists(kSet) || ! keySetMax(kSet))
    return 0 ;

  alpha = keySetAlphaHeap (kSet, keySetMax(kSet)) ;
  for (i = 0 ; i < keySetMax(alpha) ; ++i)
    if (dnaDoDumpFastAKey (keySet(alpha, i), fil, s, noMissmatch, 0, 0, 0, 0))
      ++n ;
  keySetDestroy (alpha) ;

  messout ("// I wrote %d sequences", n) ;
  return n ;
}

int dnaDumpFastAKeySet (KEYSET kSet, FILE *fil, Stack s)
{
  return dnaDoDumpFastAKeySet (kSet, fil, s, TRUE) ;
}


/************************************************************/
/************************************************************/

/*
   save BaseQuality and BasePosition information for assemblies in two new
   classes, each an array class of unsigned char (for position this is
   a delta-packing).
	BaseQuality is a measure of base quality: unsigned char
   	BasePosition is the change in position in the trace
   Both autogenerate an entry in the Sequence object of the same name, 
   with the array length following the key like DNA. 
*/

BOOL baseQualityDump (FILE *f, Stack s, KEY k)
{ int level = 0 ;
  Array a ;
  int n, x, i ; 
  unsigned char *bc ;

  a = arrayGet (k, unsigned char, "c") ;
  if (!a || !arrayMax(a))
    { arrayDestroy (a) ;
      return FALSE ;
    }

  if (f)
    level = freeOutSetFile (f) ;
  else if (s)
    level = freeOutSetStack (s) ;

  n = arrayMax (a) ;
  bc = arrp (a, 0, unsigned char) ;
  while (n)
    { for (i = 0 ; i < 50 && n ; i++, --n)
	{ x = *bc++ ;
	  freeOutf ("%d ", x) ;
	}
      freeOut ("\n") ;
    }
  freeOut ("\n") ;

  arrayDestroy (a) ;
  if (level)
    freeOutClose (level) ;
  return TRUE ;
}



BOOL basePositionDump (FILE* f, Stack s, KEY k)
{ int level = 0 ;
  Array a ;
  int n, x, i ; 
  signed char *bc ;

  a = arrayGet (k, signed char, "c") ;
  if (!a || !arrayMax(a))
    { arrayDestroy (a) ;
      return FALSE ;
    }

  if (f)
    level = freeOutSetFile (f) ;
  else if (s)
    level = freeOutSetStack (s) ;

  n = arrayMax (a) ;
  bc = arrp (a, 0, signed char) ;
  x = 0 ;
  while (n)
    { for (i = 0 ; i < 50 && n ; i++, --n)
	{ x += *bc++ ;
	  freeOutf ("%d ", x) ;
	}
      freeOut ("\n") ;
    }
  freeOut ("\n") ;

  arrayDestroy (a) ;
  
  if (level)
    freeOutClose (level) ;
  return TRUE ;
}

BOOL baseQualityParse (int level, KEY key)
{ 
  Array a = arrayCreate (1000, unsigned char) ;
  int x, n = 0 ;

  if (class(key) != _VBaseQuality)
    messcrash ("baseQualityParse called on a non-BaseQuality key") ;

  while (freecard (level) && freeint (&x))
    { do 
	{ if (x < 0 || x > 255)
	    goto abort ;
	  array(a, n++, unsigned char) = x ; 
	} while (freeint (&x)) ;
      if (freeword())
	goto abort ;
    }
  if (freeword())
    goto abort ;

  if (arrayMax(a))
    { KEY seq ;
      OBJ obj ;

      lexaddkey (name(key), &seq, _VSequence) ;
      if ((obj = bsUpdate (seq)))
	{ KEY _Quality ;

	  lexaddkey ("Quality", &_Quality, 0) ;
	  bsAddKey (obj, _Quality, key) ;
	  bsAddData (obj, _bsRight, _Int, &n) ;
	  bsSave (obj) ;
	}

      arrayStore (key, a, "c") ;
    }

  arrayDestroy (a) ;
  return TRUE ;

 abort:
  messerror ("Error parsing BaseQuality %s at line %d (not an int 0-255)", 
	     name(key), freestreamline(level)) ;
  arrayDestroy (a) ;
  return FALSE ;
}

BOOL basePositionParse (int level, KEY key)
{ 
  Array a = arrayCreate (1000, signed char) ;
  int dx, x, old = 0, n = 0 ;

  if (class(key) != _VBasePosition)
    messcrash ("basePositionParse called on a non-BasePosition key") ;

  while (freecard (level) && freeint (&x)) 
    { do 
	{ dx = x - old ;
	  if (dx > 127) dx = 127 ;
	  if (dx < -128) dx = -128 ;
	  array (a, n++, signed char) = dx ;
	  old += dx ;
	} while (freeint (&x)) ;
      if (freeword ())
	goto abort ;
    }
  if (freeword ())
    goto abort ;
  
  if (arrayMax(a))
    { KEY seq ;
      OBJ obj ;

      lexaddkey (name(key), &seq, _VSequence) ;
      if ((obj = bsUpdate (seq)))
	{ KEY _SCF_Position ;

	  lexaddkey ("SCF_Position", &_SCF_Position, 0) ;
	  bsAddKey (obj, _SCF_Position, key) ;
	  bsAddData (obj, _bsRight, _Int, &n) ;
	  bsSave (obj) ;
	}

      arrayStore (key, a, "c") ;
    }

  arrayDestroy (a) ;
  return TRUE ;

 abort:
  messerror ("Non-integer parsing BasePosition %s at line %d", 
	     name(key), freestreamline(level)) ;
  arrayDestroy (a) ;
  return FALSE ;
}

/**********************************************************/
/**********************************************************/
  /* Given a sequence keyset
   * look for dnasequences whose descriptor is in keyset
   * in each, look recursively for the coding sequences,
   *  accumulate the codon usage and display.
   * by starting from top most Source (hopefully chromosome)
   * we are garanteed not to double count.
   */

Array dnaGetWordUsage (Array kSet, int step, unsigned long *nbp, unsigned long *nwp, 
		       BOOL create, const char *fName, int comb)
{
  Array dna = 0 ; KEY seq = 0 ;
  Array words = 0 ;
  void *vp ;
  int   i, j, ii, s, iiShow ;
  unsigned int n, nw, nN, nA, nT, nC, nG, nr, nseen, nsat ;
  unsigned int w = 0, pos, rest ;
  unsigned int w1 = 0, w2 = 0, wcombed = 0 ;
  unsigned int mask,  n1 ;
  unsigned char *cp, *cq ;
  FILE *f = 0 ;
  int max = 15 ;
  unsigned long tot, jl, jltot, nb = 0 ;
  unsigned long total[17] ;  /* max + 2, but must be a constant to compile */
  BOOL debug = FALSE ;

  if (fName && !create && (!filCheckName (fName,0,"rb")))
    return 0 ;
  if (step < 1) step = 1 ;
  if (comb < 0 || comb > 2)
    {
      messcrash ("// comb value %d, out of range", comb) ;
      return 0 ;
    }
  for (i=0; i <= max + 1; i++) total[i] = 0 ;
  nw = 1 << ((2*15) -1) ;  /* 500 Mbytes */

  if (!create)
    {
  /* it is actually extremelly costly to read a 500 Mb file
   * 5 minutes to write it, 3 minutes to read it
   * and often better to just create the data on the fly
   */
      if (!fName || !*fName || !(f = filopen (fName,0,"rb")))
	return 0 ;
   
      if (debug) messout ("//reading %s", fName) ;
      words = arrayCreate(nw,  unsigned char) ; 
      array(words, nw - 1, unsigned char) = 0 ;
      vp = arrp (words, 0, unsigned char) ;
      
      nr = fread(vp, 1, nw, f) ;
      filclose (f) ;
      if (debug) messout ("//reading %s done", fName) ;
      if (nr != nw)
	{
	  filremove (fName,0) ;
	  arrayDestroy (words) ;
	}
      w = 0 ; 

      for (tot = 0, ii = 0, cp = arrp (words, 0, unsigned char) ; ii < nw ; ii++, cp++)
	{
	  nA = (*cp & 0xf ) ;
	  nT = (*cp & 0xf0) ; nT >>= 4 ;
	  tot += nA + nT ;
	}
      freeOutf ("// words file contains %lu entries\n", tot) ;
      if (nwp) *nwp = tot ;
      return words ;
    }

  /* else create the array */
  
  if (fName && (filCheckName (fName,0,"rb")))
    {
      messout ("// fil %s allready exists, kill it first", fName) ;
      return 0 ;
    }
 
  if (fName && !(f = filopen (fName,0,"wb")))
    return 0 ;
 
  nw = 1 << 29 ;
  mask = (((unsigned int)1) << (2*15)) - 1 ;
  
  words = arrayCreate(nw,  unsigned char) ; 
  array(words, nw - 1, unsigned char) = 0 ;

  total[0] = ((unsigned long) nw) << 1 ; /* because nw is 2 half bytes */
  w = 0 ; nN = nA = nT = nC = nG = nb = nseen = nsat = 0 ;

  iiShow = keySetMax(kSet)/30 ;
  if (iiShow < 1) iiShow = 1 ;
  for (ii = 0 ; ii < keySetMax(kSet) ; ii++)
    if ((dna = dnaGetWithErrors(seq=keySet(kSet,ii))))
      { 
	cq = arrp(dna, 0, unsigned char) - 1 ;
	n = 0 ; i = arrayMax(dna) ; s = step ; j = 0 ;
	if ((ii % iiShow) == 0)
	  fprintf(stderr,"// Sequence %s ::\t %9d bp\n", name(keySet(kSet,ii)), i ) ;
	while (j++, n++, cq++, s--, i--)
	  {
	    nb++ ;
	    w1 <<= 2 ; w1 |= ((w2 >> 30) & 0x3) ; /* transfer 2 bits */
	    switch (*cq)
	      {  /* ATTENTION copied from B2[] in cdnainit.h, you must keep the 2 identical */
	      case A_: nA++ ; w2 <<= 2 ; w2 |= 0x0 ; break ;
	      case G_: nG++ ; w2 <<= 2 ; w2 |= 0x1 ; break ;
	      case C_: nC++ ; w2 <<= 2 ; w2 |= 0x2 ; break ;
	      case T_: nT++ ; w2 <<= 2 ; w2 |= 0x3 ; break ;
	      default: n = 0 ;  break ;
	      }
	    switch (comb)
	      {
	      case 0:
		w = w2 & mask ;
		break ;
	      case 1:
		/* this comb :: 110110110110110110110110110110110 forgets the wobbling base */
		wcombed = 
		  (((w2 >> 2) & 0xf) << 0) |
		  (((w2 >> 8) & 0xf) << 4) |
		  (((w2 >> 14) & 0xf) << 8) |
		  (((w2 >> 20) & 0xf) << 12) |
		  (((w2 >> 26) & 0xf) << 16) |
		  (((w1 >> 0) & 0xf) << 20) |
		  (((w1 >> 6) & 0xf) << 24) |
		  (((w1 >> 12) & 0xf) << 28) ;
		w = wcombed & mask ;
		break ;
	      case 2:
		/* this comb :: 110110110110110110110110110110110 forgets the wobbling base */
		wcombed = 
		  (((w2 >> 0) & 0xf) << 0) |
		  (((w2 >> 4) & 0xf) << 4) |
		  (((w2 >> 10) & 0xf) << 8) |
		  (((w2 >> 18) & 0xf) << 12) |
		  (((w2 >> 26) & 0xf) << 16) |
		  (((w1 >> 0) & 0xf) << 20) |
		  (((w1 >> 12) & 0xf) << 24) |
		  (((w1 >> 28) & 0xf) << 28) ;
		w = wcombed & mask ;
		break ;
	      }
	    if (s == 0)
	      {
		s = step ;
		if (j >= 15 && n < 15) nN++ ; /* this word is missed because of an N */
		if (n >= 15) 
		  { 
		    nseen++ ;               /* total indexed words */
		    pos = w >> 1 ; /* divide by 2 to be in half char */
		    rest = w & 0x1 ;
		    cp =  arrp(words, pos, unsigned char) ;
		    n1 = *cp ;
		    if (rest) n1 >>= 4 ;
		    n1 &= 0x0f ;
		    if (n1 == max) nsat++ ; /* total saturated words */
		    else if (n1 < max) 
		      {
			if (total[n1] > 0) total[n1]-- ;
			n1++ ; 
			total[n1]++ ;
			if (rest) 
			  { 
			    n1 <<= 4 ; n1 &= 0xf0 ; 
			    *cp &= 0x0f ; /* zero the left half byte */
			    *cp |= n1 ;
			  }
			else
			  {
			    n1 &= 0x0f ; 
			    *cp &= 0xf0 ; /* zero the right half byte */
			    *cp |= n1 ;
			  }
		      }
		  }
	      }
	  }
	arrayDestroy(dna) ;
      }
  if (nbp) *nbp = nb ; 

  jltot = 0 ;
  freeOutf ("// step = %d\t analysed %lu bp on the positive strand \n", step, nb) ;
  freeOutf ("// Nb of different words seen ii times\n") ;
  freeOutf ("// %3s\t%12s\t%12s\n", "ii", "nn", "ii*nn") ;
    for (n = 0, i = 0 ; i <= max  ; i++)
    {
      jl = i * total[i] ;
      n += total[i] ;
      jltot += jl ;
      freeOutf ("// %3d\t%12u\t%12u\n",i,total[i], jl) ;
    }  
  if (nwp) *nwp = jltot + nsat ;

  freeOutf ("// sat\t%12s\t%12lu\n", "-", nsat) ;
  freeOutf ("//   N\t%12s\t%12lu\n", "-", nN) ;
  freeOutf ("\n// sum\t%12lu\t%12lu\n", n, jltot+nsat+nN) ;
  freeOutf ("//verif\t%12u\t%12u\n", 2*nw, nseen+ nN) ;
  tot = nN ;
  tot = tot + (unsigned long) (nA + nC) + (unsigned long) (nG + nT) ; /* casting so i 2 lines */
  freeOutf ("\n// nN = %12u\n// nA = %12u\n// nT = %12u\n// nC = %12u\n// nG = %12u\n// tot= %12lu\n\n"
	    ,nN, nA, nT, nC, nG, tot) ;

  if (f &&  
      words &&
      arrayMax(words) == nw)
    {
      vp = arrp (words, 0, unsigned char) ;
      
      messout ("// Saving the comb=%d frequency table in %s", comb, fName) ;
      nr = fwrite (vp,1, nw, f) ;
      filclose (f) ;
      
      if (nr != nw)
	filremove (fName,0 );
    }
  
  return words ;
}

/*************************************************************/
/*************************************************************/
