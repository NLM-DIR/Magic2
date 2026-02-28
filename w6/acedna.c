/*  File: acedna.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 **  All the dna utilities not dependant on the acedb database
 **    genetic code, alignments, ...
 **    but not database storage
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 27 22:14 1998 (rd)
 * * Apr 23 17:14 1995 (rd): dnaGet() now gets from a Sequence,
 		recursively finding the DNA - complex code from fmap.
 * * Oct 23 20:16 1991 (mieg): Change + to n in decodeChar
 * Created: Wed Oct 23 18:10:21 1991 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: acedna.c,v 1.71 2020/10/06 11:56:53 mieg Exp $ */

#include "../wac/ac.h"
#include "dna.h"
#include "freeout.h"
/*********************************************************************/

void showDna (Array dna, int n)
{
  int i = n + 60 ;
  unsigned char *cp ;

  if (n < 0 || !arrayExists(dna) || dna->size != 1) return ;
  if (i > arrayMax(dna)) i = arrayMax(dna) ;
  cp = arrp(dna, n, unsigned char) ;
  for (; n < i ; n++, cp++) fprintf(stderr, "%c", dnaDecodeChar [*cp]) ;
  fprintf(stderr, "\n") ;
}

/*************************************************************************/

void reverseComplement (Array dna) 
{ char  c,  *cp = arrp(dna,0,char), *cq = cp + arrayMax(dna) - 1 ;

  while (cp < cq)
    { c = complementBase (*cp) ; 
      *cp++ = complementBase(*cq) ; 
      *cq-- = c  ;
    }
  if (cp == cq)
    *cp = complementBase(*cp) ;

  /* add a terminal zero, 
     I do it again here because we might have done
     dnaR = arrayCopy(dna) ; (which looses the terminal zero)
     please use rather dnaCopy which does not lose the terminal zero
    note that once the array is allocated at size max+1
     re-imposing the zero costs nothing
     but if it is allocated at size max, it implies a new array reallocation
     reverseComplement(dnaR) ;
  */
  if (! dna->lock)
    {
      array(dna, arrayMax(dna), char) = 0 ;
      --arrayMax(dna) ;
    }
}

/**************************************************************/
/*************************************/
    /* Decode string works in a static buffer, next 3 en place */
char * dnaDecodeString(char *cp)
{ static char buf[256] ;
  int i = 255 ;
  char *cq = buf ;

  while(i-- && *cp)
    *cq++ = dnaDecodeChar[((int)*cp++) & 0xf] ;
  *cq = 0 ;
  if(*cp)
    messout
("Warning, dnaDecodeString receives toolong a buffer.") ;
  return buf ;
}

/*************************************/

void dnaDecodeArray(Array a)
{ int n = arrayMax(a) ;
  char *cp = arrp(a,0,char) ;

  cp-- ;
  while(++cp, n--)
    *cp = dnaDecodeChar[((int)*cp) & 0xf] ; /*  [] is only 16 entires */
}

/*************************************/

void dnaDecodeExtendedArray(Array a)
{ int n = arrayMax(a) ;
  char *cp = arrp(a,0,char) ;

  cp-- ;
  while(++cp, n--)
    *cp = dnaDecodeExtendedChar[((int)*cp) & 0xff] ;/* rd added a . special to alignDumpKey */
}

/*************************************/

void rnaDecodeArray(Array a)
{ int n = arrayMax(a) ;
  char *cp = arrp(a,0,char) ;

  cp-- ;
  while(++cp, n--)
    *cp = rnaDecodeChar[((int)*cp) & 0x10] ;
}

/*************************************/

void dnaEncodeString(char *cp)
{ --cp ;
  while(*++cp)
    *cp = dnaEncodeChar[((int)*cp) & 0x7f] ;
}

/*************************************/

void dnaEncodeArray(Array a)
{ int n = arrayMax(a) ;
  char *cp = arrp(a,0,char) ;

  cp-- ;
  while(++cp, n--)
    *cp = dnaEncodeChar[((int)*cp) & 0x7f] ;
}

/*************************************/
/* transform an acedb dna array into a SOLID transition dna 
 * First base is kept as is, followed by transitions,
 * so the coordinates are unchanged
 * non ATGC bases are treated as repeats of the previous base
 */
void dnaSolidEncodeArray (Array a, BOOL isDown)
{ 
  int n = a ? arrayMax(a) - 1 : 0 ;
  char cc1, cc2, *cp ;

  if (!a || n <= 0) return ;
  cp = arrp(a,0,char) ;
  cc1 = *cp ;
  switch (cc1)
    {
    case A_ :
    case T_ :
    case G_ :
    case C_ : break ;
    default : cc1 = T_ ; break ;
    }
  while(++cp, n--)
    {
      cc2 = *cp ;
      switch (cc1)
	{
	case A_ :
	  switch (cc2)
	    {
	    case A_ : *cp = (isDown ? A_ : T_) ; break ;	
	    case T_ : *cp = (isDown ? T_ : A_) ; break ;
	    case G_ : *cp = (isDown ? G_ : C_) ; break ;
	    case C_ : *cp = (isDown ? C_ : G_) ; break ;
	    default: cc2 = cc1 ; *cp = '.' ; break ;
	    }
	  break ;
	case T_ :
	  switch (cc2)
	    {
	    case A_ : *cp = (isDown ? T_ : A_) ; break ;
	    case T_ : *cp = (isDown ? A_ : T_) ; break ;
	    case G_ : *cp = (isDown ? C_ : G_) ; break ;
	    case C_ : *cp = (isDown ? G_ : C_) ; break ;
	    default: cc2 = cc1 ; *cp = '.' ; break ;
	    }
	  break ;
	case G_ :
	  switch (cc2)
	    {
	    case A_ : *cp = (isDown ? G_ : C_) ; break ;
	    case T_ : *cp = (isDown ? C_ : G_) ; break ;
	    case G_ : *cp = (isDown ? A_ : T_) ; break ;
	    case C_ : *cp = (isDown ? T_ : A_) ; break ;
	    default:  cc2 = cc1 ; *cp = '.' ; break ;
	    }
	  break ;
	case C_ :
	  switch (cc2)
	    {
	    case A_ : *cp = (isDown ? C_ : G_) ; break ;
	    case T_ : *cp = (isDown ? G_ : C_) ; break ;
	    case G_ : *cp = (isDown ? T_ : A_) ; break ;	
	    case C_ : *cp = (isDown ? A_ : T_) ; break ;
	    default: cc2 = cc1 ; *cp = '.' ; break ;
	    }
	  break ;
	}
      if (cc2 != '.') cc1 = cc2 ;
    }
  return ;
} /* dnaSolidEncodeArray */

/*************************************/
/* transform back a SOLID transition dna into an acedb dna array
 * the double transformation will treat 
 * non ATGC bases as repeats of the previous base
 */
void dnaSolidDecodeArray (Array a, BOOL isDown)
{ 
  int n = a ? arrayMax(a) - 1 : 0 ;
  char *cp, *cq ;

  if (!a || n <= 0) return ;
  cp = cq = arrp(a,0,char) ;

  messcrash ("Sorry  dnaSolidDecodeArray is not yet programmed ") ;
  cp-- ; cq-- ;
  while(++cp, ++cq, n--)
    {
      switch (*cp)
	{
	case A_ :
	  switch (*cq)
	    {
	    case A_ : *cp = (isDown ? A_ : T_ ) ; break ;
	    case T_ : *cp = (isDown ? T_ : A_ ) ; break ;
	    case G_ : *cp = (isDown ? G_ : C_ ) ; break ;
	    case C_ : *cp = (isDown ? C_ : G_ ) ; break ;
	    default: *cq = *cp ; *cp = A_ ; break ;
	    }
	  break ;
	case T_ :
	  switch (*cq)
	    {
	    case A_ : *cp = (isDown ? T_ : A_ ) ; break ;
	    case T_ : *cp = (isDown ? A_ : T_ ) ; break ;
	    case G_ : *cp = (isDown ? C_ : G_ ) ; break ;
	    case C_ : *cp = (isDown ? G_ : C_ ) ; break ;
	    default: *cq = *cp ; *cp = A_ ; break ;
	    }
	  break ;
	case G_ :
	  switch (*cq)
	    {
	    case A_ : *cp = (isDown ? G_ : C_ ) ; break ;
	    case T_ : *cp = (isDown ? C_ : G_ ) ; break ;
	    case G_ : *cp = (isDown ? A_ : T_ ) ; break ;
	    case C_ : *cp = (isDown ? T_ : A_ ) ; break ;
	    default:  *cq = *cp ; *cp = A_ ; break ;
	    }
	  break ;
	case C_ :
	  switch (*cq)
	    {
	    case A_ : *cp = (isDown ? C_ : G_ ) ; break ;
	    case T_ : *cp = (isDown ? G_ : C_ ) ; break ;
	    case G_ : *cp = (isDown ? T_ : A_ ) ; break ;
	    case C_ : *cp = (isDown ? A_ : T_ ) ; break ;
	    default: *cq = *cp ; *cp = A_ ; break ;
	    }
	  break ;
	default :
	  *cp = A_ ; 
	  break ;
	}
    }
  return ;
} /* dnaSolidDecodeArray */

/************************************/

BOOL dnaDoDump (Array dna, int debut, int fin, int offset)
{ register int i, k, k1 ;
  register char *cp, *cq ;
  char buffer [4010] ;

  if (!dna || debut < 0 || 
      debut >= fin || fin > arrayMax (dna))
    return FALSE ;
  i = fin - debut ;
  cp = arrp(dna,debut,char) ;
  cq = buffer ;

  while(i > 0)
    { cq = buffer ;
      for (k=0 ; k < 4000/(50 + offset + 3) ; k++)
        if (i > 0)
          { for (k1 = offset ; k1-- ;)
	      *cq++ = ' ' ;
            k1 = 50 ;
            while (k1--  && i--)
              *cq++ = dnaDecodeChar[*cp++ & 0xff] ;
            *cq++ = '\n' ; *cq = 0 ;
          }
      freeOut (buffer) ;
    }
  return TRUE ;
}

/**********************************************************/
   /* called also from dnacptfastaDump */

BOOL dnaDumpFastA (Array dna, int from, int to, 
		   char *text, FILE* fil, Stack s)
{ 
  int level = 0 ;

  if (!arrayExists(dna))
    return FALSE ;
  
  if (from < 0)
    from = 0 ;
  ++to ;
  if (to > arrayMax(dna))
    to = arrayMax(dna) ;

  if (fil)
    level = freeOutSetFile (fil) ;
  else if (s)
    level = freeOutSetStack (s) ;
  /* else freeOut */
  if (text) 
    freeOutf (">%s\n", text) ;

  dnaDoDump (dna, from, to, 0) ;
  if (level)
    freeOutClose(level) ;
  return TRUE ;
}

/************************************/

Array dnaParseLevel (int level, unsigned char *resultp, char **seqNamep, char **commentp, AC_HANDLE h)
{
  char *cp, c = 0, c1 = 0 ;
  Array dna = 0 ;
  register int i = 0 ;
  BOOL isFasta = FALSE ;
  char cutter ;

  freecardback () ; /* thus i can parse the obj line again */
  freecard (level) ;
  cp = freeword () ;
  if (*cp++ == '>')
    {
      isFasta = TRUE ;
      if (!*cp) { cp = freeword() ; } /* > name is allowed */
      if (!cp || !*cp) { *resultp = 1 ; return 0 ; }
      *seqNamep = strnew (cp, h) ;
      cp = freeword() ;
      if (cp)
	{ cp = freepos() ; /* additional comments */
	  if (cp)  /* we could parse them as a comment in the sequence onject */
	    * commentp = strnew (cp, h) ;
	}
    }
  dna = arrayHandleCreate(5000,char, h) ;

  while (c=0, freecard(level))
    if ((cp = freewordcut("\n", &cutter)))
      {
	if (isFasta)
	  { 
	    if (*cp == '>')  /* fasta separator */
	      { freecardback () ;
		break ;
	      }
	    if (*cp == ';') /* fasta comment */
	      continue ;
	  }
	while ((c = *cp++))
	  { if ((c1 = dnaEncodeChar[((int)c) & 0x7f]))        /* accepted base codes */
	      array(dna,i++,char) = c1 ;
	    else
	      switch(c)
		{			/* accepted tabulations */
		case '-': case 'x':  case 'X': /* x is vector masking */
		  array(dna,i++,char) = 0 ;
		  break ;
		case '*':  /* phil green padding character, ignore */
		case '\n': case '\t': case ' ':
		case '0': case '1': case '2': case '3': case '4': 
		case '5': case '6': case '7': case '8': case '9': 
		  break;
		default:		/* error in the input file */
		  c1 = 0xff ; goto abort ;
		}
	  }
      }
    else
      break ;
 abort:
  *resultp = c1 ;
  return dna ;
}

/************************************************************/
/************************************************************/
/************************************************************/
/* copied from pickMatch and modified to use bit masks */
/* match to template
   
   returns 0 if not found
           1 + pos of first sigificant match if found
*/

int dnaPickMatch (const unsigned char *cp, int n, const unsigned char *tp, int maxError, int maxN)
{
  register const unsigned char *c=cp, *t=tp, *cs=cp, *cEnd = cp + n, *tEnd = tp + strlen((const char*)tp) ;
  register int  i = maxError, j = maxN ;
 
  if (!*cp || !*tp)
    return 0 ;

  if (n > 0 && tEnd > tp + n)
    tEnd = tp + n ;
  while (c < cEnd)
  {
    if (t >= tEnd)
      return  cs - cp + 1 ;
    if (*c == N_ && --j < 0)
      { t = tp ; c = ++cs ; i = maxError ; j = maxN ; }
    else if (!(*t++ & *c++) && (--i < 0))
      { t = tp ; c = ++cs ; i = maxError ; j = maxN ; }
  } 
  return t < tEnd ? 0 : cs - cp + 1 ;
}

/*************************************************************/
/*************************************************************/
/*************************************************************/
/**********************************************************************************/
/********************************************************************/
/*********************************************************************/

void aceDnaShowErr (Array err)
{
  int j ;
  A_ERR *ep ;
  
  if (arrayExists(err))
    {
      for (j = 0  ; j < arrayMax(err) ; j++)
	{
	  ep = arrp(err, j, A_ERR) ;
	  printf("%3d type ", j) ;
	  switch (ep->type)
	    {
	    case ERREUR: printf ("ERREUR ") ; break ;
	    case AMBIGUE: printf ("AMBIGUE ") ; break ;
	    case INSERTION: printf ("INSERTION ") ; break ;
	    case TROU: printf ("TROU ") ; break ;
	    case TROU_DOUBLE: printf ("TROU_DOUBLE") ; break ;
	    case TROU_TRIPLE: printf ("TROU_TRIPLE") ; break ;
	    case INSERTION_DOUBLE: printf ("INSERTION_DOUBLE") ; break ;
	    case INSERTION_TRIPLE: printf ("INSERTION_TRIPLE") ; break ;
	    default: printf ("autre ") ; break ;
	    }
	  printf ("iLong=%d iCourt = %d baseShort = %c baseLong = %c sens = %d\n",
		  ep->iLong, ep->iShort, dnaDecodeChar[(int)ep->baseShort],  dnaDecodeChar[(int)ep->baseLong], ep->sens ) ;
	}
      aceDnaShowErr (0) ; /* self reference trick to avoid compiler warning */
    }
}

/********************************************************************/

/* in a padded system, there should be no insert delete ? */
#ifdef JUNK
static JUMP paddedJumper [] = {
{ 1, 1, 0, 0 }   /* default is punctual */
} ;
#endif

static int useJumper = 0 ;

void aceDnaSetIlmJumper (BOOL ok)
{ useJumper = ok ? 1 : 0 ; }
void aceDnaSetSolidJumper (BOOL ok)
{ useJumper = ok ? 2 : 0 ; }
void aceDnaSetRocheJumper (BOOL ok)
{ useJumper = ok ? 3 : 0 ; }
void aceDnaSetPacBioJumper (BOOL ok)
{ useJumper = ok ? 4 : 0 ; }
void aceDnaSetNanoporeJumper (BOOL ok)
{ useJumper = ok ? 5 : 0 ; }
void aceDnaSetEditGenomeJumper (BOOL ok)
{ useJumper = ok ? 6 : 0 ; }

static JUMP jumpN [] = {{1000,0,0,0}} ;  /* a kludge, to jump ambiguities */
  /* dna1 = short, dna2 = long */

JUMP jumper [] = {  /* called from abifix.c */
 {1, 0, 10, 0},    /* insert in 1 */
 {0, 1, 10, 0},    /* trou in 1 */
 {1, 1, 12, 0},    /* ponctuel */

 {1, 0, 8, 0},    /* insert in 1 */
 {0, 1, 8, 0},    /* trou in 1 */
 {2, 0, 10, 0}, 
 {0, 2, 10, 0}, 
 {3, 0, 10, 0}, 
 {0, 3, 10, 0}, 
 {1, 0, 10, 2},    /* insert in 1 */
 {0, 1, 10, 2},    /* trou in 1 */
 {2, 0, 10, 2}, 
 {0, 2, 10, 2}, 
 {3, 0, 13, 2},
 {0, 3, 13, 2},
 {4, 0, 15, 3},
 {0, 4, 15, 3},
 {5, 0, 15, 3},
 {0, 5, 15, 3},
 {6, 0, 15, 3},
 {0, 6, 15, 3},
 {7, 0, 15, 3},
 {0, 7, 15, 3},
 {8, 0, 15, 3},
 {0, 8, 15, 3},
 {9, 0, 15, 3},
 {0, 9, 15, 3},
 {10, 0, 15, 3},
 {0, 10, 15, 3},
 {1, 1, 0, 0}    /* default is punctual */
} ;

JUMP editGenomeJumper [] = {  /* called from dna2dna -editGenome */
 {1, 0, 80, 1},    /* insert in 1 */
 {0, 1, 80, 1},    /* trou in 1 */
 {1, 1, 80, 1},    /* ponctuel */

 {2, 0, 80, 1}, 
 {0, 2, 80, 1}, 
 {3, 0, 80, 1},
 {0, 3, 80, 1},
 {1, 1, 0, 0}    /* default is punctual */
} ;

JUMP ilmJumper [] = {  /* called from abifix.c */
 {1, 0, 10, 0},    /* insert in 1 */
 {0, 1, 10, 0},    /* trou in 1 */
 {1, 1, 12, 0},    /* ponctuel */

 {1, 0, 8, 0},    /* insert in 1 */
 {0, 1, 8, 0},    /* trou in 1 */
 {2, 0, 10, 0}, 
 {0, 2, 10, 0}, 
 {3, 0, 10, 0}, 
 {0, 3, 10, 0}, 
 {1, 0, 10, 2},    /* insert in 1 */
 {0, 1, 10, 2},    /* trou in 1 */
 {2, 0, 10, 2}, 
 {0, 2, 10, 2}, 
 {3, 0, 13, 2},
 {0, 3, 13, 2},
 {1, 1, 0, 0}    /* default is punctual */
} ;


static JUMP solidJumper [] = {  /* called from abifix.c */
 {1, 0, 10, 0},    /* insert in 1 */
 {0, 1, 10, 0},    /* trou in 1 */
 {1, 1, 10, 0},    /* ponctuel */

 {2, 2, 10, 0},   /* ponctuel double as in SOLID */
 {2, 1, 10, 0},    /* insert in 1 as in SOLID */
 {1, 2, 10, 0},    /* trou in 1 as in SOLID */
 {3, 1, 12, 0},    /* insert double in 1 as in SOLID */
 {1, 3, 12, 0},    /* trou double in 1 as in SOLID */
 {4, 1, 14, 0},    /* insert triple in 1 as in SOLID */
 {1, 4, 14, 0},    /* trou triple in 1 as in SOLID */

 {1, 1, 8, 0},    /* ponctuel */
 {1, 0, 8, 0},    /* insert in 1 */
 {0, 1, 8, 0},    /* trou in 1 */
 {1, 0, 10, 2},    /* insert in 1 */
 {0, 1, 10, 2},    /* trou in 1 */
 {2, 0, 10, 2}, 
 {0, 2, 10, 2}, 
 {3, 0, 13, 2},
 {0, 3, 13, 2},
 {4, 0, 15, 3},
 {0, 4, 15, 3},
 {5, 0, 15, 3},
 {0, 5, 15, 3},
 {6, 0, 15, 3},
 {0, 6, 15, 3},
 {7, 0, 15, 3},
 {0, 7, 15, 3},
 {8, 0, 15, 3},
 {0, 8, 15, 3},
 {9, 0, 15, 3},
 {0, 9, 15, 3},
 {10, 0, 15, 3},
 {0, 10, 15, 3},
 {1, 1, 0, 0}    /* default is punctual */
} ;


static JUMP insertJumper [] = {
 {1, 1, 7, 0},    /* ponctuel */
 {1, 0, 5, 1},    /* insert in 1 */
 {0, 1, 2, 0},    /* trou in 1 */
 {1, 1, 3, 0},    /* ponctuel */
 {1, 0, 9, 2},    /* insert in 1 */
 {0, 1, 4, 2},    /* trou in 1 */
 {2, 0,11, 2}, 
 {0, 2, 6, 2}, 
 {3, 0, 13, 2},
 {0, 3, 13, 2},
 {4, 0, 15, 3},
 {0, 4, 15, 3},
 {5, 0, 15, 3},
 {0, 5, 15, 3},
 {6, 0, 15, 3},
 {0, 6, 15, 3},
 {7, 0, 15, 3},
 {0, 7, 15, 3},
 {8, 0, 15, 3},
 {0, 8, 15, 3},
 {9, 0, 15, 3},
 {0, 9, 15, 3},
 {10, 0, 15, 3},
 {0, 10, 15, 3},
 {1, 1, 0, 0}    /* default is punctual */
} ;

static JUMP deleteJumper [] = {
 {1, 1, 7, 0},    /* ponctuel */
 {1, 0, 2, 0},    /* insert in 1 */
 {0, 1, 5, 1},    /* trou in 1 */
 {1, 1, 3, 0},    /* ponctuel */
 {1, 0, 4, 2},    /* insert in 1 */
 {0, 1, 9, 2},    /* trou in 1 */
 {2, 0, 6, 2}, 
 {0, 2, 11, 2}, 
 {3, 0, 13, 2},
 {0, 3, 13, 2},
 {4, 0, 15, 3},
 {0, 4, 15, 3},
 {5, 0, 15, 3},
 {0, 5, 15, 3},
 {6, 0, 15, 3},
 {0, 6, 15, 3},
 {7, 0, 15, 3},
 {0, 7, 15, 3},
 {8, 0, 15, 3},
 {0, 8, 15, 3},
 {9, 0, 15, 3},
 {0, 9, 15, 3},
 {10, 0, 15, 3},
 {0, 10, 15, 3},
 {1, 1, 0, 0}    /* default is punctual */
} ;

/* 2020_10_04
 * j'ai deja aligne seqc2 mais je change le
 * jumper pacbio (voir pacbioJumperOld below)
 * pour le rendre identique
 * au jumper nanopore utilise en depuis juin 2020
 */
static JUMP pacbioJumper [] = {
 {1, 0, 10, 0},    /* insert in 1 */
 {0, 1, 10, 0},    /* trou in 1 */
 {1, 1, 12, 0},    /* ponctuel */

 {1, 0, 8, 0},    /* insert in 1 */
 {0, 1, 8, 0},    /* trou in 1 */
 {1, 1, 7, 0},    /* trou in 1 */
 {2, 0, 10, 0}, 
 {0, 2, 10, 0}, 
 {3, 0, 10, 0}, 
 {0, 3, 10, 0}, 
 {1, 0, 10, 2},    /* insert in 1 */
 {0, 1, 10, 2},    /* trou in 1 */
 {2, 0, 10, 2}, 
 {0, 2, 10, 2}, 
 {3, 0, 13, 2},
 {0, 3, 13, 2},
 {4, 0, 15, 3},
 {0, 4, 15, 3},
 {5, 0, 15, 3},
 {0, 5, 15, 3},
 {6, 0, 15, 3},
 {0, 6, 15, 3},
 {7, 0, 15, 3},
 {0, 7, 15, 3},
 {8, 0, 15, 3},
 {0, 8, 15, 3},
 {9, 0, 15, 3},
 {0, 9, 15, 3},
 {10, 0, 15, 3},
 {0, 10, 15, 3},
 {1, 1, 6, 0},     /* sub */
 {1, 0, 6, 0},    /* insert in 1 */
 {0, 1, 6, 0},    /* trou in 1 */
 {1, 1, 0, 0}    /* default is punctual */
} ;

#ifdef JUNK
static JUMP pacbioJumperOld
 [] = {  /* favor insert */
 {1, 1, 7, 0},    /* ponctuel */
 {0, 1, 2, 0},    /* trou in 1 */
 {1, 1, 3, 0},    /* ponctuel */
 {0, 1, 4, 2},    /* trou in 1 */
 {1, 0, 5, 1},    /* insert in 1 */
 {1, 0, 9, 2},    /* insert in 1 */
 {0, 2, 6, 2}, 
 {2, 0,11, 2}, 
 {0, 3, 8, 2},
 {3, 0, 13, 2},
 {0, 4, 12, 3},
 {4, 0, 15, 3},
 {5, 0, 15, 3},
 {0, 5, 15, 3},
 {6, 0, 15, 3},
 {0, 6, 15, 3},
 {7, 0, 15, 3},
 {0, 7, 15, 3},
 {8, 0, 15, 3},
 {0, 8, 15, 3},
 {9, 0, 15, 3},
 {0, 9, 15, 3},
 {10, 0, 15, 3},
 {0, 10, 15, 3},
 {1, 1, 0, 0}    /* default is punctual */
} ;
#endif

static JUMP nanoporeJumper [] = {
 {1, 0, 10, 0},    /* insert in 1 */
 {0, 1, 10, 0},    /* trou in 1 */
 {1, 1, 12, 0},    /* ponctuel */

 {1, 0, 8, 0},    /* insert in 1 */
 {0, 1, 8, 0},    /* trou in 1 */
 {1, 1, 7, 0},    /* trou in 1 */
 {2, 0, 10, 0}, 
 {0, 2, 10, 0}, 
 {3, 0, 10, 0}, 
 {0, 3, 10, 0}, 
 {1, 0, 10, 2},    /* insert in 1 */
 {0, 1, 10, 2},    /* trou in 1 */
 {2, 0, 10, 2}, 
 {0, 2, 10, 2}, 
 {3, 0, 13, 2},
 {0, 3, 13, 2},
 {4, 0, 15, 3},
 {0, 4, 15, 3},
 {5, 0, 15, 3},
 {0, 5, 15, 3},
 {6, 0, 15, 3},
 {0, 6, 15, 3},
 {7, 0, 15, 3},
 {0, 7, 15, 3},
 {8, 0, 15, 3},
 {0, 8, 15, 3},
 {9, 0, 15, 3},
 {0, 9, 15, 3},
 {10, 0, 15, 3},
 {0, 10, 15, 3},
 {1, 1, 6, 0},     /* sub */
 {1, 0, 6, 0},    /* insert in 1 */
 {0, 1, 6, 0},    /* trou in 1 */
 {1, 1, 0, 0}    /* default is punctual */
} ;

/* Laxist jumper to fill the holes in nanopore, 2019_02 
atatcttcaa gtcttcatcctgaatccttcaat    gaaatatgtggaa  gttcacc  caat  gagaactatgttgaa  gttacattatctttag gTatctcttaata taata Caatt atcttaaatggcatttagcatatacacttatttt
atatcttcaaagtcttcatcctgaatccttcaatatgagaaatatgtggaagagttc ccttcaatatgagaactatgttgaagagttacattatctttagtgcatctcttaatactaataaaaatttatcttaaatggcatttagcatatacacttatttt

*/

static JUMP laxJumper [] = {
 {1, 1, 3, 0},    /* ponctuel */
 {0, 1, 3, 0},    /* trou in 1 */
 {0, 2, 3, 0},    /* trou in 1 */
 {0, 3, 3, 0},    /* trou in 1 */
 {0, 4, 3, 0},    /* trou in 1 */
 {1, 1, 0, 0}    /* default is punctual */
} ;

JUMP reAligner [] = {  /* called from abifix.c */
 {0, 0, 5, 1},    /* accept current location */
 {1, 0, 10, 0},    /* insert in 1 */
 {0, 1, 10, 0},    /* trou in 1 */
 {1, 1, 5, 0},    /* ponctuel */
 {1, 0, 5, 1},    /* insert in 1 */
 {0, 1, 5, 1},    /* trou in 1 */
 {1, 0, 7, 2},    /* insert in 1 */
 {0, 1, 7, 2},    /* trou in 1 */
 {2, 0, 9, 2}, 
 {0, 2, 9, 2}, 
 {3, 0, 13, 2},
 {0, 3, 13, 2},
 {4, 0, 15, 3},
 {0, 4, 15, 3},
 {5, 5, 15, 3},
 {5, 0, 15, 3},
 {0, 5, 15, 3},
 {6, 6, 15, 3},
 {6, 0, 15, 3},
 {0, 6, 15, 3},
 {7, 7, 15, 3},
 {7, 0, 15, 3},
 {0, 7, 15, 3},
 {8, 8, 15, 3},
 {8, 0, 15, 3},
 {0, 8, 15, 3},
 {9, 9, 15, 3},
 {9, 0, 15, 3},
 {0, 9, 15, 3},
 {10, 10, 15, 3},
 {10, 0, 15, 3},
 {0, 10, 15, 3},
 {0, 0, 0, 0}    /* default is accept current location */
} ;

/* maxError: 0 stop at 5 grouped errors
 *          -1 extend for ever
 *           n stop at n+1 error
 *     n <= -2 stop on first error
*/
 
/**************************************************************/

#ifdef __SSE2__   /* use 128 bit intructions  */
#define EP128
#include <emmintrin.h> // SSE2
/*  aaaa aaaa aaaa aaaa  / aaaa aaaa aaaa aaga   -> 14 (number of exact matches)
static int first_non_equal_byte(unsigned char *cp, unsigned char *cq) {
    __m128i v1 = _mm_loadu_si128((__m128i*)cp); // aaaaaaaaaaaaaaaa
    __m128i v2 = _mm_loadu_si128((__m128i*)cq); // aaaaaaaaaaaaaaga
    __m128i cmp = _mm_cmpeq_epi8(v1, v2);       // 0xff for equal, 0x00 for non-equal
    int mask = _mm_movemask_epi8(cmp);          // 0xffff for equal, bit 0 for non-equal
    return mask == 0xffff ? 16 : __builtin_ctz(~mask); // First 0 bit
}
*/
#endif


Array aceDnaTrackErrors (Array  dna1, int pos1, int *pp1, 
			 Array dna2, int pos2, int *pp2, int *NNp, Array errArray, 
			 int maxJump, int maxError, BOOL doExtend, int *maxExactp, BOOL isErrClean)
{ 
  register const char *cp, *cq, *cp1, *cq1, *cp0, *cq0, *cpmax, *cqmax ;
  register char b ;  
  /* register int  *cpmax4, *cqmax4 ; */
  register int i, n ;
  int sens = 1, nerr = 0, nAmbigue = 0, nExact, maxExact = 0 ;
  register JUMP *jp ;
  JUMP *activeJumper = jumper ;
  A_ERR *ep ;
  /* int mask4 = 0xffffffff ; */
  /* static int nnn=0 ; */
  /* on one worm chromo we have 240 M calls to this function, so chrone may be too expansive */
  /* chrono("aceDnaTrackErrors") ; */

#ifdef JUNK
  /* pedantic runtime check if we cross compile */
#ifdef EP128       /* use 128 bit intructions  */
  int hasSSE2 = __builtin_cpu_supports("sse2") ;
#else
  int hasSSE2 = 0 ;
#endif
#endif
  
  switch (useJumper)
    {
    case 1: /* illumina */
      activeJumper = ilmJumper ;
      break ;
    case 2: /* SOLiD */
      activeJumper = solidJumper ;
      break ;
    case 3: /* Roche */
      activeJumper = nanoporeJumper ;
      break ;
    case 4: /* PacBio */
      activeJumper = pacbioJumper ;
      break ;
    case 5: /* Oxford Nanopore Technology ONT: minion, gridion promethion */
      activeJumper = nanoporeJumper ;
      if (0) activeJumper = laxJumper ;
      break ;
    case 6: /* editGenome */
      activeJumper = editGenomeJumper ;
      break ;
    default:
      activeJumper = jumper ;
      break ;
    }

  if (! errArray || ! isErrClean) 
    errArray = arrayReCreate (errArray, 5, A_ERR) ;
  if (!dna1 || ! dna2)
    messcrash ("bad call to  aceDnaTrackErrors dna1 OR dna2 == 0") ;
  if (! arrayMax(dna1) || ! arrayMax(dna2))
    return errArray ;
  cp0 = arrp (dna1, 0, char) ;
  cp = cp0 + pos1 - 1 ; 
  if (*pp1 > arrayMax(dna1))
    *pp1 = arrayMax(dna1) ;
  cpmax = cp0 + (doExtend ?  arrayMax(dna1) : *pp1) ;

  cq0 = arrp (dna2, 0, char) ;
  cq = cq0 + pos2 - 1 ; 
  cqmax =  cq0 +  arrayMax(dna2) ;

  /* cpmax4 = cpmax - 4 ; cqmax4 = cqmax - 4 ; */
  nExact = 0 ;
  while (++cq, ++cp < cpmax && cq < cqmax) /* always increment both */
    {
#ifdef EP128       /* use 128 bit intructions  */
      /* if (hasSSE2) pedantic costly run time check */
	while (cp < cpmax - 16 && cq < cqmax - 16)
	  {
	    __m128i v1 = _mm_loadu_si128((__m128i*)cp);
	    __m128i v2 = _mm_loadu_si128((__m128i*)cq);
	    __m128i cmp = _mm_cmpeq_epi8(v1, v2);
	    int mask = _mm_movemask_epi8(cmp);
	    if (mask != 0xffff)
	      {
		int pos = __builtin_ctz(~mask);
		cp +=  pos ; cq += pos ;
		break; /* Mismatch at cp/cq */
	      }
	    else
	      {
		cp += 16; cq += 16; 
		nExact += 16 ;
	      }
	  }
#endif      
      if (*cp == *cq)
	{
	  nExact++ ;
	  continue ;
	}
      if  (!(*cp & 0x0F) || (*cp & 0x0F) == N_) 
	{ jp = jumpN ; goto foundjp ; }
      if  (!(*cq & 0x0F) || (*cq & 0x0F) == N_) 
	{ jp = jumpN ; goto foundjp ; }
      if (*cp & *cq) /* other kind of ambiguous letter */
	{
	  nExact++ ;
	  continue ;
	}
      if (maxError == 0)
	break ;
      jp = activeJumper ;
      jp-- ;
      while (jp++) /* 1, 1, 0, 0 = last always accepted */
	{ 
	  if (!jp->lng)
	     goto foundjp ;	  
	  if (jp->dcp - jp->dcq > maxJump ||
	      jp->dcq - jp->dcp > maxJump)
	    continue ;
	  n = 0 ; i = jp->ok ;
	  cp1 = cp + jp->dcp ; cq1 = cq + jp->dcq ; 
	  if (cp1 >= cpmax || cq1 >= cqmax)
	    continue ;
	  b = 0 ; if (jp->dcp >  jp->dcq) { b = *cp1 ; } 
 	  while (i--)
	    {
	      if (cp1 >= cpmax || cq1 >= cqmax || !((*cp1++ & 0x0F) & (*cq1++ & 0x0F)))
		goto nextjp ; 
	      if (*cp1 != b) b = 0 ;
	    } 
	  n = 0 ; i = jp->lng ;
	  if (b)
	    while (1) /* keep searching while we see the same base */
	      {
		if (cp1 >= cpmax || cq1 >= cqmax || !((*cp1 & 0x0F) & (*cq1 & 0x0F)))
		  goto nextjp ; 
		if (*cp1 != b)
		  break ;
		i-- ; cp1++; cq1++ ;
	      }
	  if (i > 0)
	    {
	      cp1 = cp + jp->dcp ; cq1 = cq + jp->dcq ;
	      if (i + cp1 >= cpmax || i + cq1 >= cqmax)
		goto nextjp ;
	      while (i--)
		if (! ((*cp1++ & 0x0F) & (*cq1++ & 0x0F)))
		  if (++n > jp->ok)
		    goto nextjp ;
	    }
	  goto foundjp ;
	nextjp:
	  continue ;	  
	}
    foundjp:
      ep = arrayp(errArray, nerr++, A_ERR) ;
      if (maxExact < nExact) maxExact = nExact ;
      nExact = 0 ;
      switch (jp->dcp - jp->dcq)
	{ 
	case 1000:
	  ep->type = AMBIGUE ;
	  ep->iShort = cp - cp0 ;
	  ep->iLong = cq - cq0 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = *cp ;
	  nAmbigue++ ;
	  break ;
	case 0:
	  ep->type = ERREUR ;
	  ep->iShort = cp - cp0 ;
	  ep->iLong = cq - cq0 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = *cp ;
	  break ;
	case 1:
	  ep->type = INSERTION ;
	  ep->iShort = cp - cp0 ;
	  ep->iLong = cq - cq0 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = *cp ;
	  cp += sens * (1 - jp->dcq) ; cp -= sens ;
	  cq -= jp->dcq ; cq-- ;
	  break ;
	case 2:
	  ep->type = INSERTION_DOUBLE ;
	  ep->iShort = cp - cp0 ;
	  ep->iLong = cq - cq0 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = *cp ;
	  cp += sens * (2 - jp->dcq) ; cp -= sens ;
	  cq -= jp->dcq ; cq-- ;
	  break ;
	case -1:
	  ep->type = TROU ;
	  ep->iShort = cp - cp0 ;
	  ep->iLong = cq - cq0 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = '*' ;
	  cp -= sens * jp->dcp ; cp -= sens ;
	  cq += (1 - jp->dcp) ; cq-- ;
	  break ;
	case -2:
	  ep->type = TROU_DOUBLE ;
	  ep->iShort = cp - cp0 ;
	  ep->iLong = cq - cq0 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = '*' ;
	  cp -= sens * jp->dcp ; cp -= sens ;
	  cq += (2 - jp->dcp) ; cq-- ;
	  break ;
	default:
	  if ( jp->dcp > jp->dcq)
	    {
	      ep->type = INSERTION_TRIPLE ;
	      ep->iShort = cp - cp0 ;
	      ep->iLong = cq - cq0 ;
	      ep->sens = sens ;
	      ep->baseLong = *cq ;
	      ep->baseShort = *cp ;
	      cp += sens * (3 - jp->dcq) ; cp -= sens ;
	      cq -= jp->dcq ; cq-- ;
	    }
	  else
	    {
	      ep->type = TROU_TRIPLE ;
	      ep->iShort = cp - cp0 ;
	      ep->iLong = cq - cq0 ;
	      ep->sens = sens ;
	      ep->baseShort = '*' ; 
	      cp -= sens * jp->dcp ; cp -= sens ;
	      cq += (3 - jp->dcp) ; cq-- ;
	    }
	}
      if (/*doStop && */
	  arrayMax(errArray) > 5 &&
	  maxError != -1 &&
	  ep->iShort < (ep - 5)->iShort + 10)
	break ;
      /* we need a cast because we allow calls with maxError<0 to stop on first error */
      if (maxError >= 0 &&  maxError != -1 && (int)arrayMax(errArray) > maxError)
	break ;
    }
  *pp1 = cp - cp0 ;
  *pp2 = cq - cq0 ;
  if (NNp) *NNp = nAmbigue ;
  if (maxExact < nExact) maxExact = nExact ;
  /*  chronoReturn () ; */
  /*
    printf("\n nnn=%d maxJump=%d\n", ++nnn, maxJump) ;
    aceDnaShowErr (errArray) ;
  */
  if (maxExactp && *maxExactp < maxExact) *maxExactp = maxExact ;
  return errArray ;
} /* aceDnaTrackErrors */

/**************************************************************************/
/* dna1 is the EST or tag to be analysed, dna2 is the reference
 * if NNp != 0, report the N in the EST, but not in the reference
 */
Array aceDnaTrackErrorsBackwards (Array  dna1, int pos1, int *pp1, 
			 Array dna2, int pos2, int *pp2, int *NNp, Array err, 
			 int maxJump, int maxError, BOOL doExtend)
{ char *cp, *cq, *cp1, *cq1, *cp0, *cq0, *cpmin, *cqmin ;
  int i, n, sens = 1, nerr = 0, nAmbigue = 0 ;
  JUMP *jp ;
  A_ERR *ep ;
  Array errArray ;
  /* static int nnn=0 ; */
  /* on one worm chromo we have 240 M calls to this function, so chrone may be too expansive */
  /* chrono("aceDnaTrackErrors") ; */

  errArray = arrayReCreate (err, 5, A_ERR) ;
  cp0 = arrp (dna1, 0, char) ;
  cp = cp0 + pos1 + 1 ; 
  if (*pp1 < 0)
    *pp1 = 0 ;
  cpmin = doExtend ? cp0 : cp0 + *pp1 ;

  cq0 = arrp (dna2, 0, char) ;
  cq = cq0 + pos2 + 1 ; 
  cqmin = cq0 ;

  while (--cq, --cp >= cpmin && cq >= cqmin) /* always increment both */
    {
      if (*cp == *cq)
	continue ;
      if  (!(*cp & 0x0F) || (*cp & 0x0F) == N_) 
	{ jp = jumpN ; goto foundjp ; }
      if  (!(*cq & 0x0F) || (*cq & 0x0F) == N_) 
	{ jp = jumpN ; goto foundjp ; }
      if (*cp & *cq) /* other kind of ambiguous letter */
	continue ;
      if (maxError == 0)
	break ;
      jp = jumper ; jp-- ;
      while (jp++) /* 1, 1, 0, 0 = last always accepted */
	{ 
	  if (!jp->lng)
	     goto foundjp ;	  
	  if (jp->dcp - jp->dcq > maxJump ||
	      jp->dcq - jp->dcp > maxJump)
	    continue ;
	  n = 0 ; i = jp->ok ;
	  cp1 = cp - jp->dcp ; cq1 = cq - jp->dcq ;
	  while (i--)
	    if (cp1 < cpmin || cq1 < cqmin || !((*cp1-- & 0x0F) & (*cq1-- & 0x0F)))
	      goto nextjp ;
	  n = 0 ; i = jp->lng ;
	  cp1 = cp - jp->dcp ; cq1 = cq - jp->dcq ;
	  if (cp1 - i < cpmin || cq1 - i < cqmin)
	    goto nextjp ;
	  while (i--)
	    if (! ((*cp1-- & 0x0F) & (*cq1-- & 0x0F)))
	      if (++n > jp->ok)
		goto nextjp ;
	  goto foundjp ;
	nextjp:
	  continue ;	  
	}
    foundjp:
      ep = arrayp(errArray, nerr++, A_ERR) ;
      
      switch (jp->dcp - jp->dcq)
	{ 
	case 1000:
	  ep->type = AMBIGUE ;
	  ep->iShort = cp - cp0 ;
	  ep->iLong = cq - cq0 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = *cp ;
	  nAmbigue++ ;
	  break ;
	case 0:
	  ep->type = ERREUR ;
	  ep->iShort = cp - cp0 ;
	  ep->iLong = cq - cq0 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = *cp ;
	  break ;
	case 1:
	  ep->type = INSERTION ;
	  ep->iShort = cp - cp0 ;
	  ep->iLong = cq - cq0 + 1 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = *cp ;
	  cp -= sens * (1 - jp->dcq) ;
	  cq += jp->dcq ;
	  break ;
	case 2:
	  ep->type = INSERTION_DOUBLE ;
	  ep->iShort = cp - cp0 - 1 ;
	  ep->iLong = cq - cq0 + 1 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = *cp ;
	  cp -= sens * (2 - jp->dcq) ;
	  cq += jp->dcq ;
	  break ;
	case -1:
	  ep->type = TROU ;
	  ep->iShort = cp - cp0 + 1 ; /* +1 because we go backwards */
	  ep->iLong = cq - cq0 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = '*' ;
	  cp += sens * jp->dcp ;
	  cq -= (1 - jp->dcp) ; 
	  break ;
	case -2:
	  ep->type = TROU_DOUBLE ;
	  ep->iShort = cp - cp0 + 1 ; /* +1 because we go backwards */
	  ep->iLong = cq - cq0 - 1 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = '*' ;
	  cp += sens * jp->dcp ;
	  cq -= (2 - jp->dcp) ;
	  break ;
	default:
	  if ( jp->dcp > jp->dcq)
	    {
	      ep->type = INSERTION_TRIPLE ;
	      ep->iShort = cp - cp0 - 2 ;
	      ep->iLong = cq - cq0 + 1 ;
	      ep->sens = sens ;
	      ep->baseLong = *cq ;
	      ep->baseShort = *cp ;
	      cp -= 3*sens ;
	    }
	  else
	    {
	      ep->type = TROU_TRIPLE ;
	      ep->iShort = cp - cp0 + 1 ;
	      ep->iLong = cq - cq0 - 2 ;
	      ep->sens = sens ;
	      ep->baseLong = *cq ;
	      ep->baseShort = '*' ;
	      cq -= 3 ;
	    }
	}
      if (/*doStop && */
	  arrayMax(errArray) > 5 &&
	  maxError != -1 &&
	  ep->iShort > (ep - 5)->iShort - 10)
	break ;
      if (maxError &&  maxError != -1 && arrayMax(errArray) > maxError)
	break ;
    }
  *pp1 = cp - cp0 + 1 ;
  *pp2 = cq - cq0 + 1 ;
  if (NNp) *NNp = nAmbigue ;

  /*  chronoReturn () ; */
  /*
    printf("\n nnn=%d maxJump=%d\n", ++nnn, maxJump) ;
    showErr (errArray) ;
  */
  return errArray ;
} /* aceDnaTrackErrorsBackwards */

/***********************************************************************/
/* to be used when the est is goind up in the mrna */

static void fuseErr (Array err1, Array err2)
{ int i, j, n1 = arrayMax(err1), n2 = arrayMax(err2) ;
  for (i = j = 0 ; i < n2 ; i++)
    {
      if (n1 == 0 || memcmp (arrp (err1, n1 + j - 1, A_ERR), arrp(err2, i, A_ERR), sizeof (A_ERR)))
	{ array(err1, n1 + j, A_ERR) = arr(err2, i, A_ERR) ; j++ ; }
    }
}

/*********************************************************************/

static Array spliceTrackAllErrors (Array  dna1, int pos1, int *pp1, 
                                   Array dna2, int pos2, int *pp2, 
                                   int *NNp, Array err, int maxJump, int maxError, BOOL doExtend, int *maxExactp)
{
  int x1 = *pp1, x2 = *pp2, y1, y2, jj = 0 ;
  Array err2 = 0 ;

  err = aceDnaTrackErrors (dna1, pos1, &x1, dna2, pos2, &x2, NNp, err, maxJump, maxError, doExtend, maxExactp, TRUE) ;
  while (jj++ < 3 && x1 < *pp1 && x2 < *pp2)
    {
      y1 = *pp1 ; y2 = *pp2 ;
      err2 = aceDnaTrackErrors (dna1, x1, &y1, dna2, x2, &y2, NNp, err2, maxJump, maxError, doExtend, maxExactp, FALSE) ;
      x1 = y1 ; x2 = y2 ;
      fuseErr (err, err2) ;  /* accumulates in err1 */
    }
  *pp1 = x1 ; *pp2 = x2 ;

  arrayDestroy (err2) ;
  return err ;
}

/***************************************************************************************/

static void spliceTrackSlideErrors (Array errors, Array dnaLong, Array dnaShort)
{
  A_ERR *ep = arrp (errors, 0, A_ERR) ;
  int ii, nerr = arrayMax (errors) ;

  for (ii = 0 ; ii < nerr ; ii++, ep++)
    {
      int x1 = ep->iShort ;
      int a1 = ep->iLong ;
      int dx = 0 ;
      int da = 0 ;
      switch (ep->type)
	{
	case TYPE80:
	case AMBIGUE:
	case ERREUR:
	  continue ;
	  break ;
	case INSERTION: dx = 1 ; break ;
	case INSERTION_DOUBLE: dx = 2 ; break ;
	case INSERTION_TRIPLE: dx = 3 ; break ;
	case TROU: da = 1 ; break ;
	case TROU_DOUBLE: da = 2 ; break ;
	case TROU_TRIPLE: da = 3 ; break ;
	}
      /* push left on read and right on target target */
      if (dx > 0)
	{
	  x1 = ep->iShort = x1 + dx - 1 ;
	  a1 = ep->iLong = a1 + 1 ;
	}
      if (da > 0)
	{
	  x1 = ep->iShort = x1 - 1 ;
	  a1 = ep->iLong = a1 + 1 - da ;
	}
      char *cS = arrp (dnaShort, x1, char) ;  /* inserted base */
      char *cL = arrp (dnaLong, a1, char)  ; /* base after the insertion */
      while (*cS == complementBase(*cL))
	{ ep->iShort-- ; ep->iLong++ ; cS-- ; cL++ ; }
      /* ep->iShort += dx ; ep->iLong += da  ; cS++ ; cL-- ;  */
    }

  return ;
} /* spliceTrackSlideErrors  */

#ifdef JUNK
static void spliceTrackSlideErrorsOld (Array err, Array dnaLong, Array dnaShort)
{
  int i, n ;
  A_ERR *ep ;
  char *cs, cc ;

  n = arrayMax (err) ;
  ep = arrp(err, 0, A_ERR) - 1 ;
  while (ep++, n--)
    { 
      switch (ep->type)
        {
        case AMBIGUE:
        case ERREUR:
	case TYPE80:
          break ;
        case TROU:
          cs = arrp(dnaLong, ep->iLong - 1, char) ;
          cc = *cs ; i = arrayMax(dnaLong) - ep->iLong ;
          while (i-- && *(++cs) == cc)
            { 
              ep->iShort-- ;
              ep->iLong++ ;
            }
          ep->iShort-- ; /* bonus */
          break ;
        case TROU_TRIPLE:
	  ep->iShort-=2 ;
	  break ;
        case TROU_DOUBLE:
          cs = arrp(dnaLong, ep->iLong - 1, char) ;
          cc = *cs ; i = arrayMax(dnaLong) - ep->iLong ;
          while (i-- && *(++cs) == cc)
            { 
              ep->iShort-- ;
              ep->iLong++ ;
            }
          ep->iLong-- ; ep->iShort-- ; /* bonus */
          break ;
        case INSERTION: 
          cs = arrp(dnaLong, ep->iLong - 1, char) ;
          cc = ep->baseShort ; i = arrayMax(dnaLong) - ep->iLong ;
          while (i-- && *(++cs) == cc)
            { 
              ep->iShort-- ;
              ep->iLong++ ;
            }
          ep->iLong++ ; /* bonus long sliding */
          break ;
        case INSERTION_DOUBLE:        
          cs = arrp(dnaLong, ep->iLong - 1, char) ;
          cc = ep->baseShort ; i = arrayMax(dnaLong) - ep->iLong ;
	  if (0)           while (i-- && *(++cs) == cc)
            { 
              ep->iShort-- ;
              ep->iLong++ ;
            }
          ep->iLong++ ; ep->iShort++ ;  /* bonus double sliding */
          break ;
        case INSERTION_TRIPLE:
          cs = arrp(dnaLong, ep->iLong - 1, char) ;
          cc = ep->baseShort ; i = arrayMax(dnaLong) - ep->iLong ;
          if (0) while (i-- && *(++cs) == cc)
            { 
              ep->iShort-- ;
              ep->iLong++ ;
            }
	  ep->iLong++ ; ep->iShort+=2 ;  /* bonus double sliding */
	  break ;
        }
    }
}
#endif

/*************************************************************/
/*********************************************************************/
/* accepts biologists coordinates (1-based, extremities included) 
 * always give *a1p < *a2p
 * if (isDown) we extend in both directions
 * else the code to extend backwards is not written  
 */




Array aceDnaDoubleTrackErrors (Array  dna1, int *x1p, int *x2p, BOOL isDown,
			       Array dna2, Array dna2R, int *a1p, int *a2p, 
			       int *NNp, Array err, int maxJump, int maxError, 
			       BOOL doExtend, int *maxExactp)
{
  int i, j, nn, u1, u2, y1 = *x1p, b1 = *a1p ;
  A_ERR *ep, *ep2, *eptmp ;
  Array err1 = 0, err2 = 0 ;
  int maxErrorOld = maxError ;
  BOOL doExtendOld = doExtend ;
  
  if (isDown)      /* Plato */
    {
      if (*x1p < 1) *x1p = 1 ;
      if (*x2p > arrayMax(dna1)) *x2p = arrayMax(dna1) ;

      if (*a1p < 1) *a1p = 1 ;
      if (*a2p > arrayMax(dna2)) *a2p = arrayMax(dna2) ;

      if (maxErrorOld == -2)
	{ maxError = -2 ; doExtend = FALSE ; }
      if (maxErrorOld == -3)
	{ maxError = -2 ; doExtend = TRUE ; }
      if (maxError == -9)
	{
	  spliceTrackAllErrors (dna1, *x1p - 1, x2p, dna2, *a1p - 1, a2p, NNp, err, maxJump, maxError
				, FALSE, maxExactp) ;
	  return err ;
	}
      err1 = arrayCreate (5, A_ERR) ; /* smaller value 5 limits the number of memset 0 */
      err1 = spliceTrackAllErrors (dna1, *x1p - 1, x2p, dna2, *a1p - 1, a2p, NNp, err1, maxJump, maxError
				   , doExtend, maxExactp) ;
      
	
      if (maxErrorOld == -2)
	{ maxError = -2 ; doExtend = doExtendOld ; }
      if (maxErrorOld == -3)
	{ maxError = -3 ; doExtend = doExtendOld ; }

      u1 = *x1p - 1 ; u2 = *a1p - 1 ; /* C type coord of first base */
      if (err1 && arrayMax (err1))
	{
	  ep = arrp (err1, 0, A_ERR) ;
	  *a1p = ep->iLong  ;
	  *x1p = ep->iShort ;
	  /* *a1p *x2p is the bio coord base after the first error */
	  if (*x1p > 1 && *a1p > 1)
	    {
	      u1 = *x1p ; u2 = *a1p ; /* C type coord of first error */
	    }
	}
      if (maxErrorOld == -2)
	{ maxError = -2 ; doExtend = TRUE ; }
      if (maxErrorOld == -3)
	{ maxError = -2 ; doExtend = FALSE ; }
      if ((u1 > y1 && u2 >= b1) || doExtend)
	  {
	    int uz1 = u1 - 1, uz2 = u2 - 1 ;
	    u1 = y1 - 1; u2 = b1 - 1 ;
            err2 = arrayCreate (5, A_ERR) ;
	    err2 = aceDnaTrackErrorsBackwards (dna1, uz1, &u1, dna2, uz2, &u2, NNp, err2, maxJump, maxError, doExtend) ;
	    *x1p = u1 + 1 ; *a1p = u2 + 1 ;
	  }
      if (maxErrorOld == -2)
	{ maxError = -2 ; doExtend = doExtendOld ; }
      if (maxErrorOld == -3)
	{ maxError = -3 ; doExtend = doExtendOld ; }
      j = 0 ; 
      if (err2 && arrayMax (err2))
	for (i = arrayMax (err2) - 1 ; i >= 0 ; i--)
	  { 
	    ep = arrp (err2, i, A_ERR) ;
	    if (ep->iShort < *x1p - 1)
	      continue ;
	    if (i > 1 && ep->type == ERREUR)
	      {
		if (*x1p < ep->iShort + 1 && ep->iShort > (ep[-1].iShort-2))
		  { *x1p = ep->iShort + 1 ; *a1p = ep->iLong + 1 ; }
	      }
	  }
      if (err2 && arrayMax (err2))
	for (i = arrayMax (err2) - 1 ; i >= 0 ; i--)
	  { 
	    ep = arrp (err2, i, A_ERR) ;
	    if (ep->iShort < *x1p + 1)
	      continue ;
	    if ( (NNp && ep->baseShort == N_) || ep->type)
	      {
		ep2 = arrayp (err, j++, A_ERR) ;
		*ep2 = *ep ;
	      }
	  }
      for (i = 0 ; err1 && i < arrayMax (err1) ; i++)
	{ 
	  ep = arrp (err1, i, A_ERR) ;
	  if (ep->iShort > *x2p - 1)
	    continue ;
	  if ((NNp && ep->baseShort == N_) || ep->type)
	    {
	      ep2 = arrayp (err, j++, A_ERR) ;
	      *ep2 = *ep ;
	    }
	}
      arrayMax (err) = j ;
      arrayDestroy (err1) ;
      arrayDestroy (err2) ;
      if (err && arrayMax (err))
	{
	  int dx ;

	  for (i = arrayMax (err) - 1, dx = 0 ; i >= 0 ; i--)
	    { 
	      ep = arrp (err, i, A_ERR) ;
	      if (i == arrayMax (err) - 1)
		{
		  if ((*x2p > ep->iShort && *x2p < ep->iShort + 5) || (*a2p > ep->iLong && *a2p < ep->iLong +5) )
		    { *x2p = ep->iShort ; *a2p = ep->iLong ; }
		}
	      if (ep->iShort == *x2p - 1)
		{ 
		  (*x2p)-- ; (*a2p)-- ; dx++ ; 
		  switch (ep->type)
		    {
		    case TROU:  (*a2p)-- ; break ;
		    case TROU_DOUBLE:  (*a2p) -= 2 ; break ;
		    case TROU_TRIPLE:  (*a2p) -= 3 ;  break ;
		    case INSERTION: (*x2p)-- ; break ;
		    case INSERTION_DOUBLE: (*x2p) -= 2 ;  break ;
		    case INSERTION_TRIPLE: (*x2p) -= 3 ; break ;
		    default:  break ;
		    }
		}
	      else
		break ;
	    }
	  arrayMax (err) -= dx ; 
	}
    }
  
  else
    {

      if (doExtend) messcrash ("Sorry aceDnaDoubleTrackErrors (isDown=FALSE doExtend=TRUE) is not written") ;
      nn = arrayMax (dna2) ; 
      if (*x2p < 1) *x2p = 1 ;
      if (*x1p > arrayMax(dna1)) *x1p = arrayMax(dna1) ;

      if (*a2p < 1) *a2p = 1 ;
      if (*a1p > arrayMax(dna2)) *a1p = arrayMax(dna2) ;

      u1 = nn - *a1p + 1 ; u2 = nn - *a2p + 1 ; /* reversed biologist coordinates */
      
      err = spliceTrackAllErrors (dna1, *x1p - 1, x2p, dna2R, u1 - 1, &u2, NNp, err, maxJump, maxError, doExtend, maxExactp) ;
      
      *a2p = nn - u2 + 1 ; /* restore */
      
      /* remove the ambiguous errors */
      if (arrayMax(err) && ! NNp)
        {
          for (i = j = 0, ep = ep2 = arrp (err, 0, A_ERR) ; i < arrayMax(err) ; ep++, i++)
            {
	      if (ep->type) /* jump AMBIGUE */
		{
		  if (i>j) *ep2 = *ep ;
		  ep2++ ; j++ ;
		}
	    }
	  arrayMax (err) = j ;	  
	}
      /* complement the result */
      if (arrayMax(err))
        {
          for (i = 0, ep = arrp (err, i, A_ERR) ; i < arrayMax(err) ; ep++, i++)
            {
              /* ep is base-0 non-bio coords */
              ep->iLong = nn - ep->iLong - 1 ;
	      ep->baseShort = complementBase(ep->baseShort) ;
            }
        }
      /* reverse the table of errors */
      if (arrayMax(err) > 1)
        {
          nn = arrayMax (err) ;
          eptmp = arrayp (err, nn, A_ERR) ; /* a good tmp place */
          for (i = 0, ep = arrp (err, i, A_ERR) ; i < arrayMax(err)/2 ; ep++, i++)
            {
              ep2 = arrp (err, nn - 1 - i, A_ERR) ;
              *eptmp = *ep2 ; *ep2 = *ep ; *ep = *eptmp ;
            }
          
          arrayMax(err) = nn ; /* restore */
        }
    }
  if (!isDown && arrayMax(err))
    spliceTrackSlideErrors (err, dna2, dna1) ;
  return err ;
}

static int trackFirstJump (Array  dna1, int pos1, Array dna2, int pos2)
{ char *cp, *cq, *cp1, *cq1 ;
  int i, n ;
  JUMP *jp ;

  if (pos1 < 0 || pos1 > arrayMax(dna1))
    return 0 ;
  cp = arrp (dna1, pos1, char) ;
  if (pos2 < 0 || pos2 > arrayMax(dna2))
    return 0 ;
  cq = arrp (dna2, pos2, char) ;

  if (*cp == *cq)
    return 0 ;

  jp = reAligner ; jp-- ;
  while (jp++) /* 1, 1, 0, 0 = last always accepted */
    { n = 0 ; i = jp->ok ;
      cp1 = cp + jp->dcp ; cq1 = cq + jp->dcq ;
      while (i--)
	if ((*cp1++ & 0x0F) != (*cq1++ & 0x0F))
	  goto nextjp ;
      n = 0 ; i = jp->lng ;
      cp1 = cp + jp->dcp ; cq1 = cq + jp->dcq ;
      while (i--)
	if ((*cp1++ & 0x0F) != (*cq1++ & 0x0F))
	  if (++n >= jp->ok)
	    goto nextjp ;
      goto foundjp ;
    nextjp:
      continue ;	  
    }
 foundjp:

  return jp->dcq - jp->dcp ;
}

/********************************************************************/
/********************************************************************/

Array trackErrors (Array  dna1, int pos1, int *pp1, 
		   Array dna2, int pos2, int *pp2, int *NNp, Array err, int mode)
{ char *cp, *cq, *cp1, *cq1, *cp0, *cq0, *cpmax, *cqmax ;
  int i, n, sens = 1, nerr = 0, nAmbigue = 0 ;
  JUMP *jp ;
  A_ERR *ep ;
  Array errArray = arrayReCreate (err, 5, A_ERR) ;
  JUMP *myJumper ;

  switch (mode)
    {
    case 1: /* favor delete */
      myJumper = deleteJumper ;
      break ;
    case 2: /* favor insert */
      myJumper = insertJumper ;
      break ;
    default:
      myJumper = jumper ;
      break ;
    }

  cp0 = arrp (dna1, 0, char) ;
  cp = cp0 + pos1 ; cpmax = cp0 + 
    ( *pp1 < arrayMax(dna1)  ? *pp1 : arrayMax(dna1)) ;
  cq0 = arrp (dna2, 0, char) ;
  cq = cq0 + pos2 ; cqmax = cq0 + arrayMax(dna2) ;
  cp-- ; cq-- ;
  while (++cp < cpmax && ++cq < cqmax)
    {
      if (*cp == *cq)
	continue ;

#ifdef JUNK
      /* old system, never allow N_ in jumps */
      if  (!(*cp & 0x0F) || (*cp & 0x0F) == N_) 
	{ jp = jumpN ; goto foundjp ; }

      jp = jumper ; jp-- ;
      while (jp++) /* 1, 1, 0, 0 = last always accepted */
	{ n = 0 ; i = jp->ok ;
	  cp1 = cp + jp->dcp ; cq1 = cq + jp->dcq ;
	  while (i--)
	    if ((*cp1++ & 0x0F) != (*cq1++ & 0x0F))
	      goto nextjp ;
	  n = 0 ; i = jp->lng ;
	  cp1 = cp + jp->dcp ; cq1 = cq + jp->dcq ;
	  while (i--)
	    if ((*cp1++ & 0x0F) != (*cq1++ & 0x0F))
	      if (++n > jp->ok)
		goto nextjp ;
	  goto foundjp ;
	nextjp:
	  continue ;	  
	}
#endif /* JUNK */
      /* allow N_ */
      if  (!(*cp & 0x0F))
	{ jp = jumpN ; goto foundjp ; }

      jp = myJumper ; jp-- ;
      while (jp++) /* 1, 1, 0, 0 = last always accepted */
	{ n = 0 ; i = jp->ok ;
	  cp1 = cp + jp->dcp ; cq1 = cq + jp->dcq ;
	  while (i--)
	    {
	      if (cp1 >= cpmax || cq1 >= cqmax)
		goto nextjp ;
	      if (! ((*cp1 & 0x0F) & (*cq1++ & 0x0F)))
		goto nextjp ;
	      if ((*cp1++ & 0x0F) == N_)
		{
		  i++ ; 
		  if (++n > jp->ok)
		    goto nextjp ;
		}
	    }
	  n = 0 ; i = jp->lng ;
	  cp1 = cp + jp->dcp ; cq1 = cq + jp->dcq ;
	  while (i--)
	    {
	      if (cp1 >= cpmax || cq1 >= cqmax)
		goto nextjp ;
	      if (! ((*cp1 & 0x0F) & (*cq1++ & 0x0F)))
		if (++n > jp->ok)
		  goto nextjp ;
	      if ((*cp1++ & 0x0F) == N_)
		{
		  i++ ; 
		  if (++n > jp->ok)
		    goto nextjp ;
		}
	    }
	  goto foundjp ;
	nextjp:
	  continue ;	  
	}
    foundjp:
      ep = arrayp(errArray, nerr++, A_ERR) ;
      
      switch (jp->dcp - jp->dcq)
	{ 
	case 1000:
	  ep->type = AMBIGUE ;
	  ep->iShort = cp - cp0 ;
	  ep->iLong = cq - cq0 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = *cp ;
	  nAmbigue++ ;
	  break ;
	case 0:
	  ep->type = *cp == N_ ? AMBIGUE : ERREUR ;
	  ep->iShort = cp - cp0 ;
	  ep->iLong = cq - cq0 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = *cp ;
	  break ;
	case 1:
	  ep->type = INSERTION ;
	  ep->iShort = cp - cp0 ;
	  ep->iLong = cq - cq0 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = *cp ;
	  cp += sens ;
	  if ((*cp & 0x0f) == N_) { cp -= sens ; cq-- ; }
	  break ;
	case 2:
	  ep->type = INSERTION_DOUBLE ;
	  ep->iShort = cp - cp0 ;
	  ep->iLong = cq - cq0 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = *cp ;
	  cp += 2*sens ;
	  if ((*cp & 0x0f) == N_) { cp -= sens ; cq-- ; }
	  break ;
	case -1:
	  ep->type = TROU ;
	  ep->iShort = cp - cp0 ;
	  ep->iLong = cq - cq0 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = '*' ;
	  cq++ ;
	  if ((*cp & 0x0f) == N_) { cp -= sens ; cq-- ; }
	  break ;
	case -2:
	  ep->type = TROU_DOUBLE ;
	  ep->iShort = cp - cp0 ;
	  ep->iLong = cq - cq0 ;
	  ep->sens = sens ;
	  ep->baseLong = *cq ;
	  ep->baseShort = '*' ;
	  cq+=2 ;
	  if ((*cp & 0x0f) == N_) { cp -= sens ; cq-- ; }
	  break ;
	default:
	  if ( jp->dcp > jp->dcq)
	    {
	      ep->iShort = cp - cp0 ;
	      ep->iLong = cq - cq0 ;
	      ep->sens = sens ;
	      ep->baseShort = *cp ;
	      ep->baseLong = *cq ;
	      /* insertion accepted, try single insertion */
	      cp1 = cp + sens ; cq1 = cq ;
	      if ((*cp1 & 0x0F) & (*cq1 & 0x0F))
		{ ep->type = INSERTION ; cp += sens ; }
	      else
		{ ep->type = INSERTION_TRIPLE; cp += 3*sens ; }
	    }
	  else
	    {
	      ep->iShort = cp - cp0 ;
	      ep->iLong = cq - cq0 ;
	      ep->sens = sens ;
	      ep->baseShort = '*' ;
	      ep->baseLong = *cq ;
	      /* trou accepted, try single trou */
	      cp1 = cp + 1 ; cq1 = cq ;
	      if ((*cp1 & 0x0F) & (*cq1 & 0x0F))
		{ ep->type = TROU ; cq++ ; }
	      else
		{ ep->type = TROU_TRIPLE ; cq += 3 ; }
	    }
	  if ((*cp & 0x0f) == N_) { cp -= sens ; cq-- ; }
	}
    }
  *pp1 = cp - cp0 ;
  *pp2 = cq - cq0 ;
  *NNp = nAmbigue ;
  return errArray ;
}


/********************************************************************/

void newLocalCptErreur(Array longDna, int xl1, int xl2, int pol,
		    Array shortDna, int xs1, int xs2, int pos, int sens,
		    int *NNp, int *startp, int *stopp, int *recouvp, Array errArray)
{ 
  if (sens == 1)
    { xl1 += trackFirstJump (shortDna, xs1, longDna, xl1) ;
      trackErrors (shortDna, xs1, &xs2, 
		   longDna, xl1, &xl2, NNp, errArray, 0) ;
      *stopp = xl2 ; *startp = xl1 ;
      *recouvp = xl2 - xl1 ;
    }
  else
/* 
    oldLocalCptErreur (longDna, xl1, xl2, pol, shortDna, xs1, xs2,
		       pos, sens, NNp, startp, stopp, recouvp, errArray) ; 
*/ messcrash ("newLocalCptErreur pas pgm") ;
  return ;
}

/********************************************************************/
/********************************************************************/
/* melting temperature of a single stranded oligo acoording to Maniatis
 * LaDeana Hillier in OSp analysis.c uses the same formula
 *
 * input is dna, x1/x2 in bio-coords
 */
float oligoTm (Array dna, int x1, int x2, float* GC_rate) 
{
  float tm = 0 ;
  int nGC = 0, N = 0 ;
  unsigned char *ccp ;

  /*  calculate the  melting temperatures for a given length, and percent gc
   * using the formula from Maniatis
   *
   *  Tm = 62.3 + 0.41*(%G+C) - 500/N where N is the length of the sequence
   */
  if (x1 > 0 && dna && x2 <= arrayMax (dna))
    {
      ccp = arrp (dna , x1 - 1, unsigned char) - 1 ;
      while (ccp++, N++, x1++ <= x2) 
	if (*ccp == G_ || *ccp == C_) 
	    nGC++ ;
    }
  if (GC_rate) *GC_rate = N > 0 ? (100.0 * nGC)/(1.0 * N) : 0 ;
  tm = N > 0 ? 62.3 + (41.0 * nGC - 500)/N : 0 ;
  
  return tm ;
} /* oligoTm */

/********************************************************************/
/* bp equivalent of the oligo complexity
 *
 */
#ifdef JUNK
we now preder entropy 16, based on frequency of dimers
static int oligoEntropy4 (unsigned const char *dna, int ln, int minEntropy) 
{
  int x1 = 0 ;
  int ss = 0 ;
  int na = 0, nt = 0, ng = 0, nc = 0, nn ;
  unsigned const char *ccp ;
  static BOOL oldNn = -1 ;
  static Array ee = 0 ;
  /*  count all letters
   */
  if (ln > 0)
    {
      ccp = dna - 1 ;
      while (ccp++, x1++ < ln) 
	switch ((int)(*ccp))
	  {
	    case A_: na++ ; break ;
	    case T_: nt++ ; break ;
	    case G_: ng++ ; break ;
	    case C_: nc++ ; break ;
	    default:  break ;
	  }
    }

  nn = na + nt + ng + nc ;

  /* lazy calculation of all the logs */
  if (! ee || (nn > 1 && nn != oldNn))
    {
      double s = 0, log4 = log(4.0) ; /* divide by log4 to be in base equivalent */ 
      int j ;

      ee = arrayReCreate (ee, nn, int) ;
      for (j = 1 ; j <= nn ; j++)
        { s = j ; s /= nn ; array (ee, j, int) = - (int) (1000.0 * j * log(s)/log4 ) ; }
      array (ee, 0, int) = 0 ;
    }

  ss = arr (ee, na, int) + arr (ee, nt, int) + arr (ee, ng, int) + arr (ee, nc, int) ;
  ss = (ss + 499)/1000 ;
  return ss ;
} /* oligoEntropy4 */
#endif

/*********************************************************************/
/* bp equivalent of the oligo complexity seen in dimers */

static int oligoEntropy16 (unsigned const char *dna, int ln, int minEntropy) 
{
  int ss = 0 ;
  int n, nn = 0, dx = ln ;
  unsigned const char *ccp ;
  static Array ee = 0 ;
  int cc, bb[16] ;
  BOOL complement = FALSE ;

  memset (bb, 0, sizeof(bb)) ;

  if (dna)
    {
      /* count all non N dimers */
      
      /* if (dx > 1000) dx = 1000 ; */
      if (dx < 0) dx = 0 ;
      for (n = nn = 0, cc = 0, ccp = dna ; dx ; dx--, ccp++)
	{
	  cc = (cc << 2) & 0xf ; n++ ;
	  switch ((int)(*ccp))
	    {
	    case A_: case 'a': if ( complement) cc |= 0x1 ; break ;
	    case T_: case 't': if (!complement) cc |= 0x1 ; break ;
	    case G_: case 'g': if ( complement) cc |= 0x3 ; else cc |= 0x2 ; break ;
	    case C_: case 'c': if ( complement) cc |= 0x2 ; else cc |= 0x3 ; break ;
	    case FS_: if (ccp[1] == RS_) complement = TRUE ; n = 0 ; break ;
	    default: n = 0 ; break ;
	    }
	  if (n >= 2) { nn++ ; bb[cc]++ ;}
	}
      if (minEntropy && nn < minEntropy)
	return FALSE ;
    }
  /* lazy calculation of all the logs */
  if (nn < 128)
    {
      double s1 = 0 ;
      if (! ee || nn * nn + nn  >= arrayMax(ee))
	{
	  double s = 0, log16 = log(16.0) ; /* divide by log16 to be in base equivalent */ 
	  int n, j, n2 ;
	  nn = 128 ;	  
	  ee = arrayReCreate (ee, nn*nn + nn, int) ;
	  for (n2 = 1 ; n2 < nn ; n2 <<= 1) ;
	  for (n = 1 ; n <= n2 ; n++)
	    {
	      for (j = 1 ; j <= n ; j++)
		{ s = j ; s /= n ; array (ee, n*n + j, int) = - (int) (1000.0 * j * log(s)/log16 ) ; }
	      array (ee, n*n, int) = 0 ;
	    }
	}
      
      for (s1 = 0, n = 0 ; dna && n < 16 ; n++)
	{ 
	  s1 += arr (ee, nn * nn + bb[n], int) ; 
	  /* fprintf (stderr,  "   %d\t%d\t%.2f\n", n, bb[n], s1) ;  */
	}
      s1 = (s1 + 499)/1000 ;
      ss = s1 ;
    }
  else
    { 
      double s = 0, s1 = 0, log16 = log(16.0) ; /* divide by log16 to be in base equivalent */ 
      ss = 0 ;
      for (ss = 0, n = 0 ; dna && n < 16 ; n++)
	if (bb[n])
	  { 
	    s = bb[n] ; s1 -=  (1000.0 * s * log(s/nn)/log16 ) ; 
	    /* fprintf (stderr, "   %d\t%d\t%.2f\n", n, bb[n], s1) ; */
	  }
      s1 = (s1 + 499)/1000 ;
      ss = s1 ;
    }

  if (minEntropy && ss < minEntropy)
    return 0 ;
  return ss ;
} /* oligoEntropy16 */

/********************************************************************/


int oligoEntropy (unsigned const char *dna, int ln, int minEntropy) 
{
  return oligoEntropy16 (dna, ln, minEntropy) ;
}

/********************************************************************/
/********************************************************************/
/********************************************************************/
