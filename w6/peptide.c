/*  File: peptide.c
 * $Header: /home/mieg/aa/CVS_ACEVIEW/ace/w6/peptide.c,v 1.30 2016/11/10 01:22:46 mieg Exp $
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Created: Tue May 10 23:59:13 1994 (rd)
 *-------------------------------------------------------------------
Thanks Jean,

I got your file.  
I also update several other tables in the peptide.c file:
the molecular weight of Sec(U) is 168.05, so the molecularWeight[] is 
changed
to 168 for U, I also update Z, B, and X molecular weight.
PKa of NH and COOH of U is update according to that of C
I did not look for the exact number but U and C shuld be close.

The molecular weight of peptide can be simply calculated by
sum(Mw of amino acids) - (length-1)*Mw of Water
if you want to be more accurate, then molecularWeight[] should be
declared as float or double.
The peptide.c file is calculating the number of atoms (kind of
 unnecessary and slower).
Kemin Zhu

 */

/* $Id: peptide.c,v 1.30 2016/11/10 01:22:46 mieg Exp $ */

#include <ctype.h>
#include <wh/acedb.h>
#include <wh/a.h>
#include <wh/bs.h>
#include <wh/peptide.h>
#include <wh/lex.h>
#include <whooks/sysclass.h>
#include <whooks/classes.h>
#include <whooks/systags.h>
#include <whooks/tags.h>
#include <wh/dna.h>
#include <wh/java.h>
#include <wh/chrono.h>

#include "freeout.h"
#include <wh/pick.h>   /* pickWord2Class */
#include <wh/dump.h>
#define  SIZEOFINT (8 * sizeof (int))

/********************************************************/
/* _VProtein is declared in classes.h/class.c */

static BOOL pepInitialise (void)
{ 
  KEY key ;
  if (!lexword2key("?Protein", &key, _VModel)  &&
      iskey(key) != 2)
    { 
      messerror("No model for class Protein, please edit wspec/models.wrm") ;
      return FALSE ;
    }
  return TRUE ;
}

/****************/

char e_reverseCodon (const char* cp, const char *translationTable)
{
  char temp[3] ;

  temp[0] = complementBase(cp[2]) ;
  temp[1] = complementBase(cp[1]) ;
  temp[2] = complementBase(cp[0]) ;
  return e_codon (temp, translationTable) ;
}
/****************/

char e_antiCodon (const char* cp, const char *translationTable)
{
  char temp[3] ;

  temp[0] = complementBase(cp[0]) ;
  temp[1] = complementBase(cp[-1]) ;
  temp[2] = complementBase(cp[-2]) ;
  return e_codon (temp, translationTable) ;
}

/****************************************/

void pepDecodeArray (Array pep)
{
  int i = arrayMax(pep) ;
  char *cp = arrp(pep,0,char) ;

  while (i--)
    { *cp = pepDecodeChar[(int)*cp] ;
      ++cp ;
    }
}

void pepEncodeArray (Array pep)
{
  int i = arrayMax(pep) ;
  signed char *cp = arrp(pep,0, signed char) ;

  while (i--)
    { *cp = pepEncodeChar[*cp] ;
      ++cp ;
    }
}

/*****************************/

static void peptideDoDump (Array pep)
{
  register int i, k, k1, fin = arrayMax(pep) ;
  register char *cp, *cq ;
  char buffer [4100] ;

  i = fin ;
  cp = arrp(pep,0,char) ;
  cq = buffer ;

  while(i > 0)
    { cq = buffer ;
      memset(buffer,0,4100) ;
      for (k=0 ; k < 4000/50  ; k++)
        if (i > 0)
          { 
            k1 = 50 ;
            while (k1--  && i--)
              *cq++ = pepDecodeChar[*cp++ & 0xff] ;
            *cq++ = '\n' ; *cq = 0 ;
          }
      freeOut (buffer) ;
    }

  return;
} /* peptideDoDump */

/*****************************/

BOOL peptideDump (FILE* f, Stack buf, KEY k) 
{
  Array pep = 0 ;
  int level = 0 ;

  if (!(pep = arrayGet(k, char, "c"))) 
    return FALSE ;

 
  if (f)
    level = freeOutSetFile (f) ;
  else if (buf)
    level = freeOutSetStack (buf) ;

  peptideDoDump(pep) ;

  if (level)
    freeOutClose(level) ;

  arrayDestroy(pep);
  return TRUE ;
} /* peptideDump */

/************************************/

BOOL javaDumpPeptide(KEY key) 
{ Array pep = arrayGet(key, char, "c") ;
 
  if (!pep)
    return FALSE ;

  freeOutf("?Peptide?%s?\t?peptide?",freejavaprotect(name(key)));

  pepDecodeArray(pep) ;
  array(pep,arrayMax(pep),char) = 0 ; /* ensure 0 terminated */
  --arrayMax(pep) ;
  freeOut (arrp(pep,0,char)) ;
  freeOut("?\n\n");

  arrayDestroy(pep) ;

  return TRUE;
} /* javaDumpPeptide */

/************************************/

BOOL peptideParse (int level, KEY key)
{
  char *cp, c = 0 ;
  signed char c1 = 0 ;
  static Array pep = 0 ;
  KEY proteinKey = 0, pepKey = 0 ;
  int i = 0 ; 
  BOOL ok = FALSE ; /* safe bet */

  if (!pepInitialise())
    return FALSE ;

  pep = arrayReCreate (pep, 1000, char) ;

  while (freecard(level) && (cp = freeword())) do
    while ((c = *cp++))
      { c1 = pepEncodeChar[((int)c) & 0x7f] ;
	if (c1 >= 0)
	  array(pep,i++,signed char) = c1 ;
	else if (c1 == -1)
	  goto abort ;
      }
    while ((cp = freeword())) ;

  /* The CDS should include the stop codon so remove it from the translation.*/
  /* N.B. we could issue a warning if it doesn't but perhaps this is         */
  /* overkill.                                                               */
  if (c == '*')			/* removing trailing '*' */
    { --arrayMax(pep) ; --i ; }

 abort:
  /* check name at end of parsing, this allows to see more errors at once */
  if (! pepReClass (key, &proteinKey) ||
      ! pepSubClass (key, &pepKey))
    {
      messerror ("  peptideParse error at line %d in %s:%s : "
	     "No tag Peptide in the model of this class",
	     freestreamline(level), className(key), name(key), c, c) ;
    }

  else if (c1 >= 0) /* succesful parsing */
    {
      if (!i || peptideStore (proteinKey, pep)) /* check i because an empty object is OK */
	ok = TRUE ;
      else
	messerror (" failed to store peptide %s", name(key)) ;	
    }

  /* failure */
  if (c1 == -1)
    messerror ("  peptideParse error at line %7d in %.25s :  bad char %d:%c\n", 
	     freestreamline(level), name(key), (int)c, c ? c : '-') ; /* avoid printing char zero */
  arrayDestroy(pep) ;
  return ok ;
} /* peptideParse */

/************************************/
/************************************/

BOOL pepSubClass (KEY protein, KEY *pepKeyp)   /* gets the parent class of pep, by default a Proteinuence */
{ 
  KEY pep = 0 ;
  if (class(protein) == _VPeptide)
    pep = protein ; 
  else if (class(protein) == _VProtein)
    lexaddkey (name(protein), &pep, _VPeptide) ;
  else if (bsIsTagInClass (class (protein), _Peptide))
    {
      char buf[1000] ;
      strcpy(buf, className(protein)) ;
      strcat (buf,":") ;
      strcat (buf, name(protein)) ;
      lexaddkey(messprintf("%s:%s", className(protein), name(protein)), &pep, _VPeptide) ;
    }
  if (pep)
    { *pepKeyp = pep ; return TRUE ; }
  return FALSE ;
}

/************************************/

BOOL pepReClass (KEY pep, KEY *proteinp)   /* gets the parent class of pep, by default a Protein */
{
  KEY protein = 0 ;
  
  if (class(pep) == _VPeptide && strlen(name(pep)))
    {
      char *cp, *buf = strnew (name(pep), 0) ;
      cp = buf ;
      for (cp=buf; *cp && *cp!= ':'; cp++) ;
      if (*cp==':')
	{
	  int cl = 0 ;
	  *cp = 0 ;
	  cl = pickWord2Class (buf) ;
	  if (cl && bsIsTagInClass (cl, _Peptide) && strlen (cp+1))
	    lexaddkey (cp+1, &protein, cl) ;
	}
      if (!pepInitialise())
	return FALSE ;
      if (!protein)
	lexaddkey (name(pep), &protein, _VProtein) ;
      messfree (buf) ;
    }
  else if (bsIsTagInClass (class(pep), _Peptide))
    protein = pep ;
  if (protein)
    { *proteinp = protein ; return TRUE ; }
  return FALSE ;
}

/*******************************/

BOOL pepInClass (int classe)   /* does this class supports peptide ? */
{
  if (classe > 0 &&
      classe < 256 &&
      (classe == _VPeptide || bsIsTagInClass (classe, _Peptide)))
    return TRUE ;
  return FALSE ;
}

/*******************************/
/*
#define MID_VAL 10

int pepchecksum(Array a)
{
  int num,tot=0;
  register i,n=1,j=0;
  int max=-99,min=99;
  
  num = (int)(arrayMax(a)/4);

  for(i=0,j=0;j<num;i+=4,j++){
    tot += (int)arr(a,i,signed char)-MID_VAL;
    tot += ((int)arr(a,i+1,signed char)-MID_VAL)*10;
    tot += ((int)arr(a,i+2,signed char)-MID_VAL)*100;
    tot += ((int)arr(a,i+3,signed char)-MID_VAL)*1000;
  }
  
  do remainder 
  for(i=num*4;i<arrayMax(a);i++){
    tot += ((int)arr(a,i,signed char)-MID_VAL)*n;
    n *=10;
  }
  if(tot)
    return tot;
  else
    return 1;
}
*/
/*******************************/
int hashArray(Array a)
{
#define  SIZEOFINT (8 * sizeof (int))
  int n = 32;
  unsigned int i, h = 0 , j = 0,
    rotate = 13 ,
    leftover = SIZEOFINT - rotate ;

  chrono("hashText") ;  
  for(i=0;i<arrayMax(a);i++)
    { h = ace_upper(pepDecodeChar[(int)arr(a,i,signed char)]) ^  /* XOR*/
        ( ( h >> leftover ) | (h << rotate)) ; 
    }
  /* compress down to n bits */
  for(j = h, i = n ; i< SIZEOFINT; i+=n)
    j ^= (h>>i) ;
  chronoReturn() ;

  return j &  0xffffffff ;
}
/*******************************/

BOOL peptideStore (KEY key, Array pep)
{ 
  OBJ obj ;
  KEY proteinKey = 0, pepKey = 0 ;
  int len, len1 ;
  int checksum = 123, ck1 = 123;

  if (!pepInitialise())
    return FALSE ;

  if (! pepReClass (key, &proteinKey) ||
      ! pepSubClass (key, &pepKey))
    {
      messerror ("peptideStore called on a non-peptide key: %s:%s", className(key), name(key)) ; 
      return FALSE ;
    }

  if (!pep || !(len = arrayMax(pep)))
    { 
      arrayKill (key) ;
      if ((obj = bsCreate (proteinKey)))
	{
	  if (bsFindKey (obj, _Peptide, pepKey))
	    { 
	      bsDestroy (obj) ;
	      if ((obj = bsUpdate(proteinKey)))
		{
		  if (obj && bsFindKey (obj, _Peptide, pepKey))
		    bsRemove (obj) ;
		  bsSave (obj) ;
		}
	      else
		return FALSE ;
	    }
	  else
	    bsDestroy (obj) ;
	}
    }
  else
    { 
      arrayStore (pepKey, pep, "c") ;
      checksum = hashArray(pep);
      if ((obj = bsCreate (proteinKey)))
	{  /* check first for sake of client server speed */
	  if (bsFindKey (obj, _Peptide, pepKey)  &&
	      bsGetData (obj, _bsRight, _Int, &len1) &&
	      len1 == len &&
	      ( bsType(obj, _bsRight) != _Int ||
		( bsGetData (obj,_bsRight, _Int, &ck1) && ck1 == checksum)
		))
	    {  
	      bsDestroy (obj) ; return TRUE ;
	    } 
	  bsDestroy (obj) ;
	}
      if ((obj = bsUpdate (proteinKey)))
	{
	  bsAddKey (obj, _Peptide, pepKey) ;
	  bsAddData (obj, _bsRight, _Int, &len) ;
	  if(bsType(obj, _bsRight) == _Int)
	    bsAddData (obj,_bsRight, _Int, &checksum);
	  bsSave(obj) ;
	}
      else
	return FALSE ;
    }
  return TRUE ;
}

/**************************************/

/* This is the old peptideGet() which only translates the object if it is a  */
/* CDS. Its just a cover function for the newer peptideTranslate() which     */
/* will translate either just the CDS or the whole object.                   */
/*                                                                           */
Array peptideGet(KEY key)
{
  Array result = NULL ;

  result = peptideTranslate(key, TRUE) ; 

  return result ;
}

/*
* given a "sequence" record, translate the DNA into the appropriate
* protein.  It uses the correct genetic_code as defined in the
* "sequence" record passed in, or the human nuclear code by default.
* 
* The returned array is encoded using the numeric name for the protein.
* If the DNA is not a multiple of 3 bases, the last amino acid will be
* X.  If there are 0 bases (i.e. not any of tcag), the matching protein
* will also be X, even if the 0 base is in a "don't care" position.
*
*/

/* Translates the dna for a particular object to its peptide representation. */
/*                                                                           */
/* If CDS_only is TRUE, then the object will be checked for a CDS tag, and   */
/* and start/end positions for the CDS will be used to define the start and  */
/* end of peptide translation within the exons.                              */
/* If CDS_only is FALSE, then all of the exons for the object will be        */
/* translated.                                                               */
/* In addition, if the Start_not_found tag is found, this will be used to    */
/* frame shift the start of translation within the object.                   */
/*             WE MAY NEED TO ALLOW AN EXTRA ARGUMENT HERE TO ADJUST THE     */
/*             THE START OF TRANSLATION VIA FMAP INTERACTIVELY....           */
/*                                                                           */
/* There is much overlap in checking of CDS & Start_not_found here with that */
/* in fmapconvert() but its difficult to produce a set of routines that are  */
/* convenient to use in both functions.                                      */
/*                                                                           */
Array peptideTranslate(KEY key, BOOL CDS_only)
{
  OBJ obj = 0 ;
  int min, max ;
  Array dna = 0, pep = 0 ;
  char *tt = 0;
  int bases,pepmax,x,cc;

  if (class(key) == _VPeptide)
    {
      pep = arrayGet(key, char, "c") ;
      goto done ;
    }

  if (!(obj = bsCreate (key)))
    return 0 ;

  if (bsGetKey (obj, _Peptide, &key) && class(key) == _VPeptide)
    { bsDestroy (obj) ;
      pep = arrayGet(key, char, "c") ;
      if (!pep) return 0 ;
      goto done ;
    }

  if (bsGetKey (obj, _Corresponding_DNA, &key))	/* try that */
    { bsDestroy (obj) ;
      if (!(obj = bsCreate (key)))
	return 0 ;
    }

  if (!(dna = dnaGet (key)))
    { bsDestroy (obj) ; return 0 ; }
  min = 1 ; max = arrayMax(dna) ;
  if (CDS_only)
    {
      if (!bsFindTag (obj, _CDS))
	{ bsDestroy (obj) ; arrayDestroy (dna) ; return 0 ; }
      bsGetData (obj, _bsRight, _Int, &min) ;
      bsGetData (obj, _bsRight, _Int, &max) ;
    }

  bases = max - min + 1 ;
  x = pepmax = bases / 3;
  if (bases % 3)
    x++;
  pep = arrayCreate (x + 1, signed char) ;
  --min ;
  x = 0 ; cc = 0 ;

  /*
  * Look for relevant  genetic_code and translation table.
  */
  tt = pepGetTranslationTable (key, 0) ;

  while (x < pepmax)
    {
      cc = e_codon (arrp(dna, min, char),tt) ;
      array (pep, x, signed char) = pepEncodeChar[cc] ;
      x++;
      min += 3;
    }
  if (bases % 3)
    array(pep, x++, signed char) = pepEncodeChar['X'];
  else
    {
  /* The CDS should include the stop codon so remove it from the translation.*/
  /* N.B. we could issue a warning if it doesn't but perhaps this is         */
  /* overkill.                                                               */

      if (0 &&           /* 2023_11-15 : export the stop */
	  cc == '*') 
	--arrayMax(pep) ;
    }

 done:
  if (pep)
    {
      array(pep, arrayMax(pep), signed char) = 0 ; /* zero terminate */
      --arrayMax(pep) ;
    }

  arrayDestroy (dna) ;
  bsDestroy (obj) ;

  return pep ;
} /* peptideGet */

/*************************************************************/
/*************************************************************/

FILE *pepFileOpen (void)
{
  static char fileName[FIL_BUFFER_SIZE],dirName[DIR_BUFFER_SIZE] ;

  return filqueryopen (dirName, fileName, "pep", "w",
		       "Choose a file for export in fasta format") ;
} /* pepFileOpen */

/**********************************************************/

   /* called also from dnacptfastaDump */

static void pepDoDump (Array a, int from, int to)
{ 
  char *cp = 0, buf [55] ;
  register int i, j, k ;
  
   for (i = from ; i < to && ! cp ;)
    { 
      k = 0 ;
      for (j = 50 ; i < to && j-- ;)
	buf [k++] = pepDecodeChar[((int)arr(a, i++, char)) & 0xff] ;
      buf [k++] = '\n' ; buf[k++] = 0 ;
      cp = strstr (buf, "*") ;
      if (cp)
	{ *(++cp) = '\n' ; *(++cp) = 0 ; }
      freeOut(buf) ;
    }

   return;
} /* pepDoDump */

/**********************************************************/

BOOL pepDumpFastA (Array a, int from, int to, const char *text, FILE* fil, Stack s)
{ 
  int level = 0 ;

  if (!a)
    return FALSE ;

  /* ATTENTION: this function used to invoke the file chooser if fil=s=0
   * but this is uncompatible with the freeout strategy needed for aceserver
   * i changed this oct 1, 1997,  mieg 
   */

  if (from < 0)
    from = 0 ;

  ++to ;
  if (to > arrayMax(a))
    to = arrayMax(a) ;


  if (fil)
    level = freeOutSetFile (fil) ;
  else if (s)
    level = freeOutSetStack (s) ;
  /* else freeOut */

  freeOutf (">%s\n", text) ;
  pepDoDump(a,from,to) ;
  if (level)
    freeOutClose(level) ;

  return TRUE ;
} /* pepDumpFastA */

/**********************************************************/

static BOOL pepDumpCStyle (KEY key, Array pep, int from, int to)
{
  int ii = to - from + 1 ;
  char buf[2] ;
  char *cp = name(key) ;

  if (!pep || from < 0 || from >= to || to >= arrayMax(pep))
    return FALSE ;

  pepDecodeArray (pep) ;
  buf[0] = 'N' ; buf[1] = class(key) ;
  freeOutBinary (buf, 2) ;
  freeOutBinary (cp, strlen(cp) + 1) ; 

  buf[0] = '\n' ;
  buf[1] = 'c' ;
  freeOutBinary (buf,2) ; 
  
  freeOutBinary ((char*) &ii, 4) ; 

  buf[0] = '\n' ;
  buf[1] = 'p' ;
  freeOutBinary (buf,2) ; 

  freeOutBinary (arrp (pep, from, char), ii) ;

  buf[0] = '\n' ; 
  buf[1] = '#' ; 
  freeOutBinary (buf,2) ; 
  
  return TRUE ;
}
 
/**********************************************************/

BOOL peptideDumpKeyCstyle (KEY key)
{
  Array a = peptideGet (key) ;
    BOOL result = FALSE ;

  if (a)
    { 
      result = pepDumpCStyle (key, a, 0, arrayMax(a) - 1) ;
      arrayDestroy (a) ;
    }

  return result ;
}

/**********************************************************/

BOOL pepDumpFastAKey (KEY key, FILE *fil, Stack s, char style)
{ 
  Array a = peptideGet (key) ;
  BOOL result = FALSE ;

  if (a)
    { 
      if (style == 'C')
	result = pepDumpCStyle (key, a, 0, arrayMax(a) - 1) ;
      else
	result = pepDumpFastA (a, 0, arrayMax(a)-1, name(key), fil, s) ;
      arrayDestroy (a) ;
      
    }
  
  return result ;
}

/**********************************************************/

int pepDumpFastAKeySet (KEYSET kSet, FILE *fil, Stack s)
{
  KEYSET alpha ;
  int i, n = 0 ;

  if (!kSet || (!fil && !stackExists (s) &&!(fil = pepFileOpen())))
    return 0 ;

  alpha = keySetAlphaHeap (kSet, keySetMax(kSet)) ;
  for (i = 0 ; i < keySetMax(alpha) ; ++i)
    if (pepDumpFastAKey (keySet(alpha, i), fil, s, 0))
      ++n ;
  keySetDestroy (alpha) ;

  messout ("I wrote %d sequences", n) ;
  return n ;
}

/**********************************************************/
/**********************************************************/
/*
* Extended Codon translator
*
* This translates 3 dna bases into the equivalent protein.  Because the
* mapping is different in different places, you really need to use a map
* that is sensitive to where the DNA came from.
*
* The Sequence object in the database contains a field named Genetic_Code.
* That field refers to an object of type genetic_code that defines the
* mapping.  You need two steps to make the translation:
* 1) use pepGetTranslationTable() to translate the genetic_code key into
*    a char * that is used for the translation.  (returns null on error)
* 2) pass the char * to e_codon() along with the 3 bases to translate.
*
* You never free the value returned by pepGetTranslationTable(); it is
* cached internally.
*
*/

/*
* e_codon() - do the actual translation
* 
* s points at the 3 bytes of DNA code to translate, encoded as in dna.h
*
* map is a string containing the 64 byte codon map to use.  Get it
* from pepGetTranslationTable()
*
* The return value is the single character protein identifier or X if
* it cannot be known.
*
* ( It would be easy enough to return a list of all the proteins it
* could be, if that were of interest. )
*
* This function always examines all 64 possibilities.  You could
* optimize it by recognizing specific bit patterns and duplicating
* the loop bodies, but it probably isn't worth bothering.
*/
char e_codon (const char *s, const char *map)
{
  int x,y,z;
  char it=0;

  for (x=0 ; x<4 ; x++)
    if (s[0] & (1<<x))
      {
        for (y=0; y<4; y++)
	  if (s[1] & (1<<y))
	    {
	      for (z=0; z<4; z++)
		if (s[2] & (1<<z))
		  {
		    if (!it)
		      it = map[(x<<4)|(y<<2)|z] ;
		    if (map[(x<<4)|(y<<2)|z] != it)
		      return 'X';
		  }
	    }
      }
  if (!it) return 'X';	/*  I think this never happens */
  return it ;
} /* e_codon */

/*
* helper function for pepGetTranslationTable()
*/

static unsigned char dnaStrictEncodeChar[] =
{  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,

   0,   1,   0,   4,   0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   1,   0,   4,   0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
} ;

static int t_codon_code (char *s)
{
  unsigned char a, b, c;
    
  a = dnaStrictEncodeChar [(int)s[0]] ;
  b = dnaStrictEncodeChar [(int)s[1]] ;
  c = dnaStrictEncodeChar [(int)s[2]] ;
  if (!a-- || !b-- || !c--)
    return -1 ;
  return (a<<4) | (b<<2) | c ;
} /* t_codon_code */

/********************************/

static KEY pepGetGeneticCode (KEY key) 
{
  KEY gKey = 0 ;

  if (class (key) == _VGenetic_code)
    gKey = key ;
  else
    while (key && !gKey)
      {
	gKey = keyGetKey (key, _Genetic_code) ;
	if (!gKey)
	  {
	    key = keyGetKey (key, _Source) ; /* standard system */
	    if (!key)
	      key = keyGetTag2Key (key, _S_Parent) ; /* smap parent */
	  }
      }

  return gKey ;
} /* pepGetGeneticCode */

/********************************/

char *pepGetTranslationTable (KEY key, KEY *geneticCodep)
{
  OBJ obj = 0 ;
  char *translationTable = 0 ;
  KEY gKey = 0 ; int index = 0 ;
  static Array maps = 0 ;

  /*
   * Create maps.
   */
  if (!maps)
    maps = arrayCreate (64, char*) ;

  /*
   * find the correct genetic_code
   */
  if (!(gKey = pepGetGeneticCode (key)))
    goto done ;
  
  /*
   * see if we have already built a map for it
   */
  index = KEYKEY (gKey) ;
  translationTable = array (maps, index, char *) ;
  if (translationTable)
    goto done ;

   /*
   * we have not.  fetch the data components from the object
   * and create a map.
   */
  
  if ((obj = bsCreate(gKey)))
    {
      char *base1, *base2, *base3, *translation ;
      char result[4*4*4+1] ;
      int x,n ;

      /*
       * these 4 fields encode the translation table.
       */
      if (bsGetData (obj, str2tag("Translation"), _Text, &translation) &&
	  bsGetData (obj, str2tag("Base1"), _Text, &base1) &&
	  bsGetData (obj, str2tag("Base2"), _Text, &base2) &&
	  bsGetData (obj, str2tag("Base3"), _Text, &base3))
	{
	  /*
	   * if anything doesn't get filled abort, the answer is X
	   */
	  memset(result,0,4*4*4) ;
	  result[4*4*4] = 0 ;

	  /*
	   * If there are not 64 bytes of translation data, something
	   * is wrong.
	   */
          if (strlen(translation) != 64) goto bad;
          if (strlen(base1      ) != 64) goto bad;
          if (strlen(base2      ) != 64) goto bad;
          if (strlen(base3      ) != 64) goto bad;

	  /*
	   * walk through the combinations and do the conversion.  I don't
	   * assume that the bases are in a particular order, even though it
	   * looks like it from the database.
	   */
	  for (x=0 ; x<64 ; x++)
	    {
	      char xl[3] ;
	      xl[0] = base1[x] ;
	      xl[1] = base2[x] ;
	      xl[2] = base3[x] ;
	      n = t_codon_code (xl) ;
	      if (n >= 0)
		result[n] = translation[x] ;
	    }

 	  /*
	   * if we have not set all 64 translations, the genetic code is 
	   * corrupted.
	   */
	  if (strlen(result) != 64) goto bad;
	  
	  /*
	   * success, store the result
	   */
	  translationTable = array (maps, index, char*) = strnew(result, 0) ;
	}
      else
        {
	/* if any tags are missing from the object, the genetic code is bad.  */
	goto bad;
	}
    }
  
 done:
  bsDestroy(obj) ;
  /* default to the standard human nuclear genetic translationTable */
  if (!translationTable)
    translationTable = "KNKNIIMIRSRSTTTT*Y*YLFLF*CWCSSSSEDEDVVVVGGGGAAAAQHQHLLLLRRRRPPPP" ;
  
  if (gKey && geneticCodep)
    *geneticCodep = gKey ;
  return translationTable ;

bad:
  messerror("Genetic code %s is corrupted - assuming human nuclear code", name(gKey));
  translationTable =  array (maps, index, char*) =
    "KNKNIIMIRSRSTTTT*Y*YLFLF*CWCSSSSEDEDVVVVGGGGAAAAQHQHLLLLRRRRPPPP" ;
  goto done;

}  /* pepGetTranslationTable */ 

/***************** end of file *******************/
