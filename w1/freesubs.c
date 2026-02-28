/*  File: freesubs.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 * -------------------------------------------------------------------
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: free format input - record based
 * Exported functions: lots - see regular.h
 * HISTORY:
 * Last edited: May  9 12:18 2002 (edgrif)
 * * Apr 26 10:06 1999 (edgrif): Add code to trap NULL value for freedouble.
 * * Jan 14 13:31 1999 (fw): increased version umber to 1.1.2
 *              due to changes in the directory handling interface
 *              (Array)filDirectoryCreate -> (FilDir)filDirHandleCreate
 * * Dec  3 14:46 1998 (edgrif): Insert version macros for libfree.
 *    freecard ignores "\r" and "\n" under WIN32
 * * Sep 30 14:19 1998 (edgrif)
 * * Nov 27 12:30 1995 (mieg): freecard no longer stips \, freeword does
 *    also i added freeunprotect
 * Created: Sun Oct 27 18:16:01 1991 (rd)
 * CVS info:   $Id: freesubs.c,v 1.48 2020/05/30 16:50:30 mieg Exp $
 *-------------------------------------------------------------------
 */

#undef CHRONO

#include "regular.h"
#include "version.h"
#include <ctype.h>
#include "mytime.h" 
#include "call.h" 
#include "ac.h" 
#include "chrono.h" 
/* free package version and copyright string.    */
/*                                               */
#define FREE_TITLE   "Free library"
#define FREE_DESC    "Sanger Centre Informatics utilities library."

#define FREE_VERSION 1
#define FREE_RELEASE 1
#define FREE_UPDATE  2
#define FREE_VERSION_NUMBER  UT_MAKE_VERSION_NUMBER(FREE_VERSION, FREE_RELEASE, FREE_UPDATE)

/* UT_COPYRIGHT_STRING(FREE_TITLE, FREE_VERSION, FREE_RELEASE, FREE_UPDATE, FREE_DESC) */

 
BOOL isInteractive = TRUE ;      /* can set FALSE, i.e. in tace */

 
#define MAXSTREAM 80
#define MAXNPAR 80
typedef struct
  { FILE *fil ;
    char filName[1024] ;
    const char *text ;
    char special[24] ;
    int npar ;
    int parMark[MAXNPAR] ;
    int line ;
    BOOL isPipe ;
    int bufPos ;
    char *bufp ;
    char buff[4096] ;
  } STREAM ;
static STREAM   stream[MAXSTREAM] ;
static int	streamlevel ;
static FILE	*currfil ;	/* either currfil or currtext is 0 */
static const char	*currtext ;	/* the other is the current source */
static Stack    parStack = 0 ;
static int	maxcard = 4096 ;
static unsigned char *card = 0, *word = 0, *cardEnd, *pos ;
static Associator filAss = 0 ;
static char *bufp, *buff ;  /* line buffer in current stream */
static int bufPos = 0 ;

#define _losewhite    while (*pos == ' '|| (special[(int) '\t'] && *pos == '\t')) ++pos
#define _stepover(x)  (*pos == x && ++pos)

#if 0 /* use getc == old acedb method of 1989 */
#define _FREECHAR     (currfil ? getc (currfil) : *currtext++)
#define buffUngetc(_ch,_fil) (_fil ? ungetc(_ch,_fil) : currtext--)

#else /* use gets == new hopefully fatser method mieg 2003 */

/************************************/
/* mieg: july 2003: optimise parsing by using gets() */

static char bufRead (void)
{
  static int isEscape = 0 ;

  if (! fgets (buff, 4096, currfil)
      || feof(currfil)
      )
    {
      bufPos = 0 ;
      isEscape = 0 ;
      return 0 ;
    }
  else
    {      
      bufPos = strlen(buff) ;
      bufp = buff ;
      if (bufPos >=2 && buff[bufPos-1]=='\n' && buff[bufPos-2]=='\\')
	{ isEscape |= 0x2 ; bufPos -= 2 ; } /* effectivelly espaces the end of line */
      if (isEscape & 0x1) /* then jump next time the leading tab */
	{
	  if (*bufp == '\t')
	    { bufp++ ; bufPos-- ; }
	}
      if (isEscape & 0x2)
	isEscape = 0x1 ;
      else
	isEscape = 0 ;
      bufPos--; /* one char transmitted now */
      return *bufp++ ;
    }
}

#if 1 /* use gets, use macro : i.e. inline the next 3 functions */
/* i tried a macro to inline bufRead, but surprisingly it was actually a tiny bit slower on linus and on alpha */
#define _FREECHAR    (currfil ? (bufPos ? (bufPos--, *bufp++) : bufRead()) : *currtext++)
#define buffUngetc(_ch,_fil) (_fil ? (bufPos++,bufp--) : currtext--)
#else /* use gets, no-macro */


#define _FREECHAR FREECHAR()
static unsigned char  FREECHAR (void)
{
  unsigned char ch = 0 ;
  
  if (currfil)
    {
      if (bufPos)
	{
	  bufPos-- ;
	  ch = *bufp++ ;
	}
      else
	ch = bufRead() ;
    }
  else
    ch = *currtext++ ;
  return ch ;
}

static void buffUngetc (char ch, FILE *fil)
{
  if (currfil)
    {
      bufPos++ ;
      bufp-- ;
    }
  else
    --currtext ;
}

#endif /* use gets, no-macro */
#endif /* use getc/use gets */

/************************************/
 
void freeshutdown (void)
{
  messfree (card) ;
  messfree (word) ;
  stackDestroy (parStack) ;
  assDestroy (filAss) ;
}

/************************************/
extern void freeOut (const char *txt) ;

void freeinit (void)
{ static BOOL isInitialised = FALSE ;
  
  if (!isInitialised)
    { streamlevel = 0 ;
      currtext = 0 ;
      stream[streamlevel].fil = currfil = stdin ;
      stream[streamlevel].text = 0 ;
      freespecial ("\n\t\"\\/@%") ;
      card = (unsigned char *) messalloc (maxcard) ;
      cardEnd = &card[maxcard-1] ;
      pos = card ;
      word = (unsigned char *) messalloc (maxcard) ;
      filAss = assCreate () ;
      parStack = stackCreate (128) ;
      buff = stream[streamlevel].buff ;
      bufPos = 0 ; /* will force bufRead */
      isInitialised = TRUE ;
      freeOut (0) ; /* needed be the linker 2017_02_26 */
      /*  filAddDir (".") ; */
    }
}

/*******************/

/* Sometimes you may need to know if a function below you succeeded in       */
/* changing a stream level.                                                  */
int freeCurrLevel(void)
  {
  return streamlevel ;
  }

/*******************/

static unsigned char* freeExtend (unsigned char *in)
{	/* only happens when getting card */
  unsigned char *oldCard = card, *out ;

  maxcard *= 2 ;
  card = (unsigned char *) messalloc (maxcard) ;
  if (oldCard)     /* jtm june 22, 1992 */
    memcpy (card, oldCard, maxcard/2) ;
  cardEnd = &card[maxcard-1] ;
  out = in + (card - oldCard) ;
  messfree (oldCard) ;
  messfree (word) ;
  word = (unsigned char *) messalloc (maxcard) ;

  return out ;
}
 
/********************/

static char special[256] ;

void freespecial (char* text)
{
  if (!text)
    messcrash ("freespecial received 0 text") ;
  if (strlen(text) > 23)
    messcrash ("freespecial received a string longer than 23") ;
  if (text != stream[streamlevel].special)
    strcpy (stream[streamlevel].special, text) ;
  memset (special, 0, (mysize_t) 256) ;
  while (*text)
    special [((int) *text++) & 0xFF] = TRUE ;
  special[0] = TRUE ;
  special[(unsigned char)EOF] = TRUE ;		/* better have these to ensure streams terminate! */
}
 
/********************/
 
void freeforcecard (const char *string)
{ int level = freesettext (string, "") ;
  freespecial ("\"") ;
  freecard (level) ;
}
 
/********************/

char* freecard (int level)	/* returns 0 when streamlevel drops below level */
{ 
  register unsigned char *in,ch,*cp ;
  int kpar ;
  int isecho = FALSE ;		/* could reset sometime? */
  FILE *fil ;
  BOOL acceptShell, acceptCommand ;

restart :
  if (level > streamlevel)
    return 0 ;

  if (isecho)
    printf (!currfil ? "From text >" : "From file >") ;
  in = card - 1 ;

  acceptCommand = special[(int)'@'] ? TRUE : FALSE ;
  acceptShell = special[(int)'$']  ? TRUE : FALSE ;

  while (TRUE)
    {
      if (++in >= cardEnd)
	{ 
	  in = freeExtend (in) ;
	}

      *in = _FREECHAR ;
    lao:
      if (special[((int) *in) & 0xFF] && *in != '$' && *in != '@' )
	switch (*in)
	  {
#if defined(WIN32)
	  case '\r':
		  continue ; /* ignore carriage returns */ 
#endif
	  case '\"':
	    continue ;  /* will be treated by _losewhite */
	  case '\n':		/* == '\x0a' */
	  case ';':		/* card break for multiple commands on one line */
	    goto got_line ;
	  case (unsigned char) EOF:
	  case '\0':
	    if (currfil && currfil != stdin && currfil != stdout &&
		! stream[streamlevel].isPipe)
	      {
		long nn1, nn2 ;
		nn1 = ftell (currfil) ;
		fseek (currfil, 0, SEEK_END) ;
		nn2 = ftell (currfil) ;
		if (nn1 != nn2)
		  messerror ("File %s line %d freecard hit an EOF at %ld before true EOF %ld"
			     , stream[streamlevel].filName
			     , stream[streamlevel].line
			     , nn1, nn2) ;
	      }
	    freeclose(streamlevel) ;
	    goto got_line;
	  case '\t':     /* tabs should get rounded to 8 spaces */
	    if (isecho)	/* write it out */
	      putchar (*in) ;
            *in++ = ' ' ;
            while ((in - card) % 8)
              { if (in >= cardEnd)
                  in = freeExtend (in) ;
                *in++ = ' ' ;
              }
	    --in ;
	    continue ;  
	  case '/':		/* // means start of comment */
	    if ((ch = _FREECHAR) == '/')
	      { 
		while ((ch = _FREECHAR) != '\n' && ch && ch != (unsigned char)EOF) ;
		*in = ch ;
		goto lao ; /* mieg, nov 2001, to treat the EOF, this was breaking the server if sent a comment */
	      }
	    else
	      { 
		if (isecho) putchar (*in) ;
		buffUngetc (ch, currfil) ; /* push back ch */
	      }
	    break ;
	  case '%':		/* possible parameter */
	    --in ; kpar = 0 ;
	    while (isdigit (ch = _FREECHAR))
	      kpar = kpar*10 + (ch - '0') ;
	    if (kpar > 0 && kpar <= stream[streamlevel].npar)
	      for (cp = (unsigned char *) stackText (parStack, 
			     stream[streamlevel].parMark[kpar-1]) ; *cp ; ++cp)
		{ if (++in >= cardEnd)
                  in = freeExtend (in) ;
		  *in = *cp ;
		  if (isecho)
		    putchar (*in) ;
		}
	    else
	      messout ("Parameter %%%d can not be substituted", kpar) ;
	    if (++in >= cardEnd)
	      in = freeExtend (in) ;
	    *in = ch ; 
	    goto lao ; /* mieg */
	  case '\\':		/* escapes next character - interprets \n */
	    *in = _FREECHAR ;
	    if (*in == '\n')    /* fold continuation lines */
	      { if (isInteractive && !streamlevel)
		  printf ("  Continuation >") ;
		while ((ch = _FREECHAR) == ' ' || ch == '\t') ;
			/* remove whitespace at start of next line */
		buffUngetc (ch, currfil) ; /* push back ch */	
		stream[streamlevel].line++ ;
		--in ;
	      }
#if !defined(WIN32)
	    else if (*in == 'n') /* reinterpret \n as a format */
	      { *in = '\n' ; 
	      }
#endif
	    else  /* keep the \ till freeword is called */
	      { *(in+1) = *in ;
		*in = '\\' ;
		if (++in >= cardEnd)
                  in = freeExtend (in) ;
	      }
	    break ;
	  default:
	    messerror ("freesubs got unrecognised special character 0x%x = %c\n",
		     *in, *in) ;
	  }
      else
	{ if (iscntrl((int)*in) && *in != '\t' && *in != '\n') /* mieg 2015_10_12, replaced isprint by ! isctrl to allow parsing Göndör, Millán-Ariño, Németi and other accentuated characters */
	    --in ;
	  else if (isecho)	/* write it out */
	    putchar (*in) ;
	}
    }				/* while TRUE loop */
  
 got_line:
  
  chrono ("freecard") ;
  
  stream[streamlevel].line++ ;
  *in = 0 ;
  if (isecho)
    putchar ('\n') ;
  pos = card ;
  _losewhite ;
  if (acceptCommand && _stepover ('@'))        /* command file */
    { char *name ;
      if ((name = freeword ()) && 
	  (fil = filopen (name, 0, "r")))
	freesetfile (fil, (char*) pos) ;
      chronoReturn () ; 
      bufPos = 0 ; /* this card is over */
      goto restart ;
    }
  if (acceptShell && _stepover ('$'))        /* shell command */
    {
#if !defined(MACINTOSH)
      callSystem ((char*)pos) ;
#endif
  
      *card = 0 ; pos = card ; 
      bufPos = 0 ; /* this card is over */ 
      goto restart ;
    }
  chronoReturn () ;
  bufPos = 0 ; /* this card is over */
  return (char*) card ;
}
 
/************************************************/

void freecardback (void)    /* goes back one card */
{ stream[streamlevel].line-- ;
  freesettext ((char*) card, "") ;
  stream[streamlevel].line = stream[streamlevel - 1].line ;
}

/************************************************/

/* reads card from fil, has some embedded logic for dealing with
 * escaped newlines which doesn't use the freespecial() stuff.
 * I don't know the reason for this and probably it could be added
 * but this is now used exclusively by acediff.c
 */
BOOL freeread (FILE *fil)
{
  unsigned char ch, *in = card ;
  int  *line, chint ;
  
  /* freespecial ("\n\t\"\\@%") ; we use to have this but acediff call freeread and needs \t not to be in 'freespecial' */
  
  if (!assFind (filAss, fil, &line))
    {
      line = (int*) messalloc (sizeof (int)) ;
      assInsert (filAss, fil, line) ;
    }
  
  --in ;
  while (TRUE)
    {
      ++in ;
      if (in >= cardEnd)
	in = freeExtend (in) ;
      chint = getc(fil) ;
      if (ferror(fil))
	messerror ("chint was bad");
      *in = chint ;
      switch (*in)
        {
	case '\n' :
	  ++*line ;
	case (unsigned char) EOF :
	  goto got_line ;
	case '/' :		/* // means start of comment */
	  if ((ch = getc (fil)) == '/')
	    { 
	      while (getc(fil) != '\n' && !feof(fil)) 
		{} ;
	      ++*line ;
	      if (in > card)	/* // at start of line ignores line */
		goto got_line ;
	      else
		--in ; /* in = 0   unprintable, so backstepped */
	    }
	  else
	    ungetc (ch,fil) ;
	  break ;
	  
	case '\\' :		/* escape next character, NOTE fall through. */
	  *in = getc(fil) ;
	  if (*in == '\n')	/* continuation */
	    { ++*line ;
	      while (isspace (*in = getc(fil))) ;    /* remove whitespace */
	    }
	  else if (*in == '"' || *in == '\\' || *in == '/') /* escape for freeword */
	    { *(in+1) = *in ;
	      *in = '\\' ;
	      ++in ;
	    }
	  /* NB fall through - in case next char is nonprinting */
	  
	default:
	  {
	    /* mieg 2015_10_12, replaced isprint by ! isctrl to allow parsing Göndör, Millán-Ariño, Németi and other accentuated characters */
	    if (iscntrl (*in) && *in != '\t')	/* ignore control chars, e.g. \x0d */
	      --in ;
	  }
	}
    }
 got_line :
  *in = 0 ;
  pos = card ;
  _losewhite ;
  if (feof(fil))
    { assRemove (filAss, fil) ;
      messfree (line) ;
    }
  return (*pos != 0) || !feof(fil) ? TRUE : FALSE ;
}

int freeline (FILE *fil)
{ int *line ;

  if (assFind (filAss, fil, &line))
    return *line ;
  else
    return 0 ;
}
 
int freestreamline (int level)
{ 
  return stream[level].line  ;
}
 
/********************************************/


static void freenewstream (const char *parms)
{
  int nn, kpar ;

  stream[streamlevel].fil = currfil ;
  stream[streamlevel].text = currtext ;
  stream[streamlevel].bufPos = bufPos ;
  stream[streamlevel].bufp = bufp ;

  if (++streamlevel == MAXSTREAM)
    messcrash ("MAXSTREAM overflow in freenewstream") ;
  strcpy (stream[streamlevel].special, stream[streamlevel-1].special) ;

  stream[streamlevel].npar = 0 ;
  stream[streamlevel].line = 1 ;
  (stream[streamlevel].filName)[0] = 0 ;
  buff = stream[streamlevel].buff ;
  bufPos = 0 ; /* will force bufRead */

 if (!parms || !*parms)
    return ;                           /* can t abuse NULL ! */
  pos = (unsigned char *) parms ;			/* abuse freeword() to get parms */
  nn = strlen(parms) ;
  while (nn > maxcard - 1000)  /* allow very long parameters (dna in table-maker conditions) */
    freeExtend (card) ; /* doubles maxcard */
  for (kpar = 0 ; kpar < MAXNPAR && freeword () ; kpar++) /* read parameters */
    { stream[streamlevel].parMark[kpar] = stackMark (parStack) ;
      pushText (parStack, (char*) word) ;
    }

  stream[streamlevel].npar = kpar ;
  stream[streamlevel].isPipe = FALSE ;
  pos = card ;			/* restore pos to start of blank card */
  *card = 0 ;
}
 
int freesettext (const char *string, const char *parms)
{
  int nn = string ? strlen(string) : maxcard ;

  freenewstream (parms) ;

  if (nn > maxcard - 1000)
    {
      while (nn > maxcard - 1000) 
	freeExtend (card) ; /* doubles maxcard */
    }
  currfil = 0 ;
  currtext = string ;

  return streamlevel ;
}

int freesetfile (FILE *fil, const char *parms)
{
  freenewstream (parms) ;

  currfil = fil ;
  currtext = 0 ;

  return streamlevel ;
}

int freesetnamedfile (char *filName, FILE *fil, const char *parms)
{
  freesetfile (fil, parms) ;
  strncpy (stream[streamlevel].filName, filName ? filName : "", 1023) ;
  return streamlevel ;
}

int freesetpipe (FILE *fil, const char *parms)
{
  freenewstream (parms) ;

  currfil = fil ;
  currtext = 0 ;

  stream[streamlevel].isPipe = TRUE ;
  return streamlevel ;
}

void freeclose(int level)
{ 
  int kpar ;

  while (streamlevel >= level)
    { if (currfil && currfil != stdin && currfil != stdout)
	{
	  if (stream[streamlevel].isPipe)
	    pclose (currfil) ;
	  else
	    filclose (currfil) ;
	}
      for (kpar = stream[streamlevel].npar ; kpar-- ;)
	popText (parStack) ;
      --streamlevel ;
      currfil = stream[streamlevel].fil ;
      currtext = stream[streamlevel].text ;
      freespecial (stream[streamlevel].special) ;
      buff = stream[streamlevel].buff ;
      bufp = stream[streamlevel].bufp ;
      bufPos = stream[streamlevel].bufPos ;
    }
}

/************************************************/
/* freeword(), freewordcut() and freestep() are the only calls that
     directly move pos forward -- all others act via freeword().
   freeback() moves pos back one word.
*/
 
char *freeword (void)
{
  unsigned char *cw ;
 
  _losewhite ;             /* needed in case of intervening freestep() */

  if (special[(int) '"'] && _stepover ('"'))
    { for (cw = word ; !_stepover('"') && *pos ; *cw++ = *pos++)
	if (special[(int) '\\'] && _stepover('\\'))	/* accept next char unless end of line */
	  if (!*pos)
	    break ;
      _losewhite ;
      *cw = 0 ;
      return (char*) word ;	/* always return a word, even if empty */
    }

		/* default: break on space and \t, not on comma */
  for (cw = word ; isgraph (*pos) && *pos != '\t' ; *cw++ = *pos++)
    if (special[(int) '\\'] && _stepover('\\'))	/* accept next char unless end of line */
      if (!*pos)
	break ;
  _losewhite ;
  *cw = 0 ;
  return *word ? (char*) word : 0 ;
}
 
/************************************************/

#if defined(WIN32)

char *freepath (void)
{
  unsigned char *cw ;
 
  _losewhite ;             /* needed in case of intervening freestep() */

  if (_stepover ('"'))
    { for (cw = word ; !_stepover('"') && *pos ; *cw++ = *pos++)
	if (_stepover('\\'))	/* accept next char unless end of line */
	  if (!*pos)
	    break ;
      _losewhite ;
      *cw = 0 ;
      return (char*) word ;	/* always return a word, even if empty */
    }

  /* default: break on space, \t or end of line, not on comma
	 also, does not skip over backslashes which are assumed to be
	 MS DOS/Windows path delimiters */
  for (cw = word ; ( *pos == '\\' || isgraph (*pos) ) && *pos != '\t' ; *cw++ = *pos++) ;
  _losewhite ;
  *cw = 0 ;
  return *word ? (char*) word : 0 ;
}

#endif
 
/************************************************/
 
char *freewordcut (char *cutset, char *cutter)
        /* Moves along card, looking for a character from cut, which is a
           0-terminated char list of separators.
           Returns everything up to but not including the first match.
           pos is moved one char beyond the character.
           *cutter contains the char found, or if end of card is reached, 0.
        */
{ unsigned char *cc,*cw ;
 
  for (cw = word ; *pos ; *cw++ = *pos++)
    for (cc = (unsigned char *) cutset ; *cc ; ++cc)
      if (*cc == *pos)
        goto wcut ;
wcut:
  *cutter = *pos ;
  if (*pos)
    ++pos ;
  _losewhite  ;
  *cw-- = 0 ;
  /* lose trailing spaces in the word */
  while (cw >= word && *cw==' ')*cw--=0;
  return *word ? (char*) word : 0 ;
}
 
/************************************************/
 
void freeback (void)    /* goes back one word - inefficient but reliable */
 
 {unsigned char *now = pos ;
  unsigned char *old = pos ;
 
  pos = card ; _losewhite ;
  while  (pos < now)
   {old = pos ;
    freeword () ;
   }
  pos = old ;
 }
 
/************************************************/
/* Read a word representing an int from wherever free is pointing to.        */
/*                                                                           */
/* If there is no word OR the word cannot be converted to an int, reset      */
/*   the free pos, don't set the int param. and return FALSE                 */
/* If the word is "NULL" set int param to the POSIX "too small" int value    */
/*   and return TRUE                                                         */
/* Otherwise set the int param to the converted int and return TRUE          */
/*                                                                           */
/* Note that valid range of ints is    INT_MIN < valid < INT_MAX             */
/* otherwise UT_NON_INT doesn't work.                                        */
/*                                                                           */
BOOL freeint (int *p)
 
 {unsigned char *keep = pos ;
  unsigned char *cp ;
  int value = 0 ;
  BOOL isMinus = FALSE ;
 
  if (freeword ())
    { /*printf ("freeint got '%s'\n", word) ;*/
      cp = word ;
      if (!strcmp ((char*)cp, "NULL") || !strcasecmp ((char*)cp, "nan"))
	{ *p = UT_NON_INT ;
	  return TRUE ;
	}
      if (*cp == '-')
        { isMinus = TRUE ;
          ++cp ;
        }
      while (*cp)
      	{ if (*cp >= '0' && *cp <= '9')
      	    value = value*10 + (*cp++ - '0') ;
      	  else
      	    { pos = keep ;
      	      return FALSE ;
      	    }
      	 }
   	  *p = isMinus ? -value : value ;
      return (TRUE) ;
    }
  else
    { pos = keep ;
      return (FALSE) ;
    }
 }
 
/*****************************/
/* Read a word representing a float from wherever free is pointing to.       */
/*                                                                           */
/* If there is no word OR the word cannot be converted to a float, reset     */
/*   the free pos, don't set the float param. and return FALSE               */
/* If the word is "NULL" set float param to the POSIX "too small" float value*/
/*   and return TRUE                                                         */
/* Otherwise set the float param to the converted float and return TRUE      */
/*                                                                           */
/* Note that valid range of floats is:                                       */
/*                             -ve  -FLT_MAX < valid <  -FLT_MIN             */
/*                             +ve   FLT_MIN < valid <   FLT_MAX             */
/* otherwise UT_NON_FLOAT doesn't work as a range check for applications.    */
/*                                                                           */

BOOL freefloat (float *p)
{
  unsigned char *keep = pos ;
  float old = *p ;
  char dummy ; 
  double xx ;
 
  if (freeword ())
    { 
      if (!strcmp ((char*)word, "NULL") || !strcasecmp ((char*)word, "nan"))
	{ *p = UT_NON_FLOAT ;
	return TRUE ;
	}
      if (sscanf ((char*) word,"%lf%c",&xx,&dummy) == 1 &&
	  xx > -FLT_MAX  && xx < FLT_MAX 
	  )
	{ *p = xx ; return (TRUE) ; }
    }

  pos = keep ;
  *p = old ;
  return (FALSE) ;
} /* freefloat */
 
/* always move to the next cutter, return TRUE if a float was scanned */
BOOL freefloatcut (float *p, char *cutset, char *cutter)
{
  char *cp ;
  float old = *p ;
  char dummy ; 
 
  if ((cp = freewordcut (cutset, cutter)))
    { if (!strcmp ((char*)cp, "NULL") || !strcasecmp ((char*)cp, "nan"))
	{ *p = UT_NON_FLOAT ;
	  return TRUE ;
	}
      if (sscanf ((char*) cp,"%f%c",p,&dummy) == 1)
	return (TRUE) ;
    }
  else if (*cutter) /* nothing between cutters means NULL */
    {
      *p = UT_NON_FLOAT ;
      return TRUE ;
    }
 *p = old ;
  return (FALSE) ;
} /* freefloatcut */

/* always move to the next cutter, return TRUE if a float was scanned */
BOOL freedoublecut (double *p, char *cutset, char *cutter)
{
  char *cp ;
  double old = *p ;
  char dummy ; 
 
  if ((cp = freewordcut (cutset, cutter)))
    { if (!strcmp ((char*)cp, "NULL") || !strcasecmp ((char*)cp, "nan"))
	{ *p = UT_NON_DOUBLE ;
	  return TRUE ;
	}
      if (sscanf ((char*) cp,"%lf%c",p,&dummy) == 1)
	return (TRUE) ;
    }
  else if (*cutter) /* nothing between cutters means NULL */
    {
      *p = UT_NON_DOUBLE ;
      return TRUE ;
    }
 *p = old ;
  return (FALSE) ;
} /* freefloatcut */

/**************************************************/
/* Read a word representing a double from wherever free is pointing to.      */
/*                                                                           */
/* If there is no word OR the word cannot be converted to a double, reset    */
/*   the free pos, don't set the double param. and return FALSE              */
/* Otherwise set the double param to the converted double and return TRUE    */
/*                                                                           */
/* Note that valid range of doubles is:                                      */
/*                             -ve  -DBL_MAX < valid <  -DBL_MIN             */
/*                             +ve   DBL_MIN < valid <   DBL_MAX             */
/* otherwise UT_NON_DOUBLE doesn't work as a range check for applications.   */
/*                                                                           */
/* Note that because double is the largest number we can only check to see   */
/* if the converted number is equal to this, not >= as in the case of floats.*/
/*                                                                           */
BOOL freedouble (double *p)
{ 
  unsigned char *keep = pos ;
  double old = *p ;
  char dummy ;
 
  if (freeword () && (sscanf ((char*) word,"%lf%c",p,&dummy) == 1))
    return (TRUE) ;
  else
    { pos = keep ;
      *p = old ;
      return (FALSE) ;
    }
}
 
/*************************************************/
 
static int ambiguouskey;
 
BOOL freekey (KEY *kpt, FREEOPT *options)
{
  unsigned char  *keep = pos ;

  if (!freeword())
    return FALSE ;

  if (freekeymatch ((char*) word, kpt, options))
    return TRUE;
 
  if (ambiguouskey)
    messout ("Keyword %s is ambiguous",word) ;
  else if (word[0] != '?')
    messout ("Keyword %s does not match",word) ;
 
  pos = keep ;
  return FALSE ;
}
 
/*****************/
 
BOOL freekeymatch (char *cp, KEY *kpt, FREEOPT *options)
{
  char  *io,*iw ;
  FREEOPT *opt ;
  int  len, nopt ;
  KEY key = 0 ;
 
  ambiguouskey = FALSE;

  /* check for bad arguments */
  nopt = (int)options->key ;
  if (!nopt || !cp)
    return FALSE ;
 
  /* check for full words */
  opt = options ; len = strlen (cp) ;
  while (opt++, nopt--)
    { iw = cp ;
      io = opt->text ;

      if (!strncasecmp (iw, io, len) && (! io[len] || io[len] == ' '))
	{ 
	  *kpt = opt->key ; 
	  return TRUE ; 
	}
    }

  /* not a full word match */
  opt = options ; nopt = (int)options->key ;
  while (opt++, nopt--)
    { iw = cp ;
      io = opt->text ;

      while (ace_upper (*iw++) == ace_upper(*io++))
	if (!*iw)
	  { 
	    key = opt->key ; 
	    goto done ;
	  }
    }
  return FALSE ;

 done:
  while (opt++, nopt--)	/* check that later options are different */
    { 
      io = opt->text ;
      iw = (char*) word ;
      while (ace_upper (*iw++) == ace_upper (*io++))
	if (!*iw)
	  { ambiguouskey = TRUE;
	    return FALSE ;
	  }
    }
 
  *kpt = key ;
  return TRUE ;
}
 
/***************************************************/
  /* Return the text corresponding to the key */
char *freekey2text (KEY k, FREEOPT *o)  
{ int i = o->key ; char *title = o->text ;
  if (i<0)
    messcrash("Negative number of options in freekey2text") ;
  while (o++, i--) 
    if (o->key == k)
      return (o->text) ;
  return title ;
}

/***************************************************/
 
BOOL freeselect (KEY *kpt, FREEOPT *options)     /* like the old freemenu */
{
  if (isInteractive)
    printf ("%s > ",options[0].text) ;
  freecard (0) ;                       /* just get a card */
  if (isInteractive)
    while (freestep ('?'))            /* write out options list */
      { int i ;
	for (i = 1 ; i <= options[0].key ; i++)
	  printf ("  %s\n",options[i].text) ;
	printf ("%s > ",options[0].text) ;
	freecard (0) ;
      }
  return freekey (kpt,options) ;
}
  /* same but returns TRUE, -1, if stremlevel drops below level */
BOOL freelevelselect (int level, KEY *kpt, FREEOPT *options)     /* like the old freemenu */
{
  if (isInteractive)
    printf ("%s > ",options[0].text) ;
  if (!freecard (level)) /* try to get another card */
    { *kpt = (KEY)(-1) ;  
      return TRUE ;                       
    }

  if (isInteractive)
    while (freestep ('?'))            /* write out options list */
      { int i ;
	for (i = 1 ; i <= options[0].key ; i++)
	  printf ("  %s\n",options[i].text) ;
	printf ("%s > ",options[0].text) ;
	if (!freecard (level)) /* try to get another card */
	 { *kpt = (KEY)(-1) ;  
	   return TRUE ;                       
	 }
       }
  return freekey (kpt,options) ;
} 

/**************************************/
 
BOOL freequery (const char *query)
{
  if (isInteractive)
    { int retval, answer = 0 ;
      printf ("%s (y or n) ",query) ;
      answer = getchar () ;
      retval = (answer == 'y' || answer == 'Y') ? TRUE : FALSE ;
      while (answer != (unsigned char) EOF &&
	     answer != -1 && /* mieg: used not to break on EOF in pipes */
	     answer != '\n')
        answer = getchar () ;
      return retval ? TRUE : FALSE ;
    }
  else
    return TRUE ;
}
 
/**********/
 
BOOL freeprompt (const char *prompt, char *dfault, char *fmt)
{ 
  if (isInteractive)
    printf("%s ? > ",prompt);
  freecard (0) ;                       /* just get a card */
  if (freecheck (fmt))
    return TRUE ;
  else
    { messout ("input mismatch : format '%s' expected, card was\n%s",
	       fmt, card) ;
      return FALSE ;
   }
}
 
/*************************************/
 
int freefmtlength (char *fmt)
 
 {char *cp ;
  int length = 0 ;
 
  if (isdigit((int)*fmt))
   {sscanf (fmt,"%d",&length) ;
    return length ;
   }
 
  for (cp = fmt ; *cp ; ++cp)
    switch (*cp)
     {
case 'i' : case 'f' : case 'd' : length += 8 ; break ;
case 'w' : length += 32 ; break ;
case 't' : length += 80 ; break ;
case 'o' :
      if (*++cp)
        messcrash ("'o' can not end free format %s",fmt) ;
      length += 2 ; break ;
     }
 
  if (!length)
    length = 40 ;
  return length ;
 }
 
/****************/
 
BOOL freecheck (char *fmt)
        /* checks that whatever is in card fits specified format
           note that 't' format option changes card by inserting a '"' */
 {unsigned char *keep = pos ;
  union {int i ; float r ; double d ;}
          target ;
  char *fp ;
  unsigned char *start ;
  int nquote = 1 ;
 
  for (fp = fmt ; *fp ; ++fp)
    switch (*fp)
     {
case 'w' : if (freeword ()) break ; else goto retFALSE ;
case 'i' : if (freeint (&target.i)) break ; else goto retFALSE ;
case 'f' : if (freefloat (&target.r)) break ; else goto retFALSE ;
case 'd' : if (freedouble (&target.d)) break ; else goto retFALSE ;
case 't' :      /* must insert '"' and escape any remaining '"'s or '\'s */
      for (start = pos ; *pos ; ++pos)
        if (*pos == '"' || *pos == '\\')
          ++nquote ;
      *(pos+nquote+1) = '"' ;		/* end of line */
      for ( ; pos >= start ; --pos)
	{ *(pos + nquote) = *pos ;
	  if (*pos == '"' || *pos == '\\')
	    *(pos + --nquote) = '\\' ;
        }
      *start = '"' ;
      goto retTRUE ;
case 'z' : if (freeword ()) goto retFALSE ; else goto retTRUE ;
case 'o' :
      if (!*++fp) messcrash ("'o' can not end free format %s",fmt) ;
      freestep (*fp) ; break ;
case 'b' : break; /* special for graphToggleEditor no check needed  il */
default :
      if (!isdigit((int)*fp) && !isspace((int)*fp))
        messerror ("unrecognised char %d = %c in free format %s",
		   *fp, *fp, fmt) ;
     }
 
retTRUE :
  pos = keep ; return TRUE ;
retFALSE :
  pos = keep ; return FALSE ;
 }
 
/************************ little routines ************************/
 
BOOL freestep (char x)
 {return (*pos && ace_upper (*pos) == x && pos++) ? TRUE : FALSE  ;
 }
 
void freenext (void)
 {_losewhite ;
 }
 
char* freepos (void)		/* cheat to give pos onwards */
{ return (char*) pos ;
}

/*************************************************************************/
/*************************************************************************/

char* freeunprotect (const char *text)
{
  static char *buf = 0 ;
  char *cp, *cp0, *cq ;
  messfree (buf) ;
  buf = text ? strnew(text, 0) : strnew (" ", 0) ;

  /* remove external space and tabs and first quotes */
  cp = buf ;
  while (*cp == ' ' || *cp == '\t') cp++ ;
  if (*cp == '"') cp++ ;
  while (*cp == ' ' || *cp == '\t') cp++ ;

  cq = cp + strlen(cp) - 1 ;

  while (cq > cp && (*cq == ' ' || *cq == '\t')) *cq-- = 0 ;

  if (*cq == '"') /* remove one unprotected quote */
    {
      int i = 0 ; char *cr = cq - 1 ;
      while (cr > cp && *cr == '\\')
	{ i++ ; cr-- ; }
      if ( i%2 == 0)
	*cq-- = 0 ;  /* discard */
    }
  while (cq > cp && (*cq == ' ' || *cq == '\t')) *cq-- = 0 ;

  /* gobble the \ */
  cp0 = cq = cp-- ;
  while (*++cp)
    switch (*cp)
      {
      case '\\': 
	if (*(cp+1) == '\\') { cp++ ; *cq++ = '\\' ; break ;}
	if (*(cp+1) == '\n') { cp ++ ; break ; } /* skip backslash-newline */
	if (*(cp+1) == 'n') { cp ++ ; *cq++ = '\n' ; break ; }
	break ;
      default: *cq++ = *cp ;
      }
  *cq = 0 ;   /* terminate the string */
  return cp0 ;
}

/*************/

char* freeprotect_old (const char* text)	/* freeword will read result back as text */
{
  static Array a = 0 ;
  const char *ccp ;
  char *cq ;
  int base ;

		/* code to make this efficiently reentrant */

  if (a && arrayMax(a) && text >= arrp(a,0,char) && text < arrp(a,arrayMax(a)-1,char))
    { base = text - arrp(a,0,char) ;
      array (a, base+3*(1+strlen(text)), char) = 0 ; /* ensure long enough */
      text = arrp(a,0,char) + base ;            /* may have relocated */
      base += 1 + strlen(text) ;
    }
  else
    { a = arrayReCreate (a, 2*(1+strlen(text))+1, char) ;
      base = 0 ;
      array (a, 2*(1+strlen(text)), char) = 0 ; /* ensure long enough */
    }

  cq = arrp (a, base, char) ;
  *cq++ = '"' ;
  for (ccp = text ; *ccp ; )
    { 
      if (*ccp == '\\'  || *ccp == '"' || 		       /* protect these */
	  *ccp == '/' || *ccp == '%' || *ccp == ';' ||
	  *ccp == '\t' 
	  /* || *ccp == '\n' NO this is done on next line mieg:2017-11-21 */
	  )
	*cq++ = '\\' ;
      if (*ccp == '\n') 
	{*cq++ = '\\' ; *cq++ = 'n' ; ccp++ ; } /* do not acedump a \n, bad for other scripts */
      else
	*cq++ = *ccp++ ;
    }
  *cq++ = '"' ;
  *cq = 0 ;
  return arrp (a, base, char) ;
}

/**********/

char* freeprotect (const char* text)	/* freeword will read result back as text */
{
  static char *result = 0 ;

  if (result)
    messfree (result) ;

  result = ac_protect (text, 0) ;
  return result ;
}

/*************/

char* freejavaprotect (char* text)	/* freeword will read result back as text */
{
  static Array a = 0 ;
  char *cp, *cq ;
  int base ;

		/* code to make this efficiently reentrant */

  if (a && text >= arrp(a,0,char) && text < arrp(a,arrayMax(a),char))
    { base = text - arrp(a,0,char) ;
      array (a, base+3*(1+strlen(text)), char) = 0 ; /* ensure long enough */
      text = arrp(a,0,char) + base ;            /* may have relocated */
      base += 1 + strlen(text) ;
    }
  else
    { a = arrayReCreate (a, 128, char) ;
      base = 0 ;
      array (a, 2*(1+strlen(text)), char) = 0 ; /* ensure long enough */
    }

  cq = arrp (a, base, char) ;
  cp = text;
  while (*cp) 
    switch (*cp)
      {
      case '\n':
	*cq++ = '\\';
	*cq++ = 'n';
	cp++;
	break;
      case '\\': case '?':
	*cq++ = '\\' ;
	/* fall thru */
      default:
	*cq++ = *cp++;
      }
  *cq = 0 ;
  return arrp (a, base, char) ;
}

static QueryRoutine	  queryRoutine = 0 ;
static PromptRoutine	  promptRoutine = 0 ;

QueryRoutine messQueryRegister (QueryRoutine func)
{ QueryRoutine old = queryRoutine ; queryRoutine = func ; return old ; }

PromptRoutine messPromptRegister (PromptRoutine func)
{ PromptRoutine old = promptRoutine ; promptRoutine = func ; return old ; }

BOOL messPrompt (char *prompt, char *dfault, char *fmt)
{ 
  BOOL answer ;
  
  if (promptRoutine)
    answer = (*promptRoutine)(prompt, dfault, fmt) ;
  else
    answer = freeprompt (prompt, dfault, fmt) ;

  return answer ;
}
/*****************************/

BOOL messQuery (char *format,...)
{ 
  BOOL answer = FALSE ;
  Stack s = stackCreate (500) ;
  int len ;
  char *cp ;
  va_list ap;
  
  va_start(ap, format);
  len = utPrintfSizeOfArgList (format, ap) ;
  va_end (ap) ;
  
  stackExtend (s, len + 1) ;
  cp = stackText (s, 0) ;
  
  va_start(ap, format);
  len = vsprintf(cp, format, ap) ;
  va_end (ap) ;
  
  if (queryRoutine)
    answer = (*queryRoutine)(cp) ;
  else
    answer = freequery (cp) ;
  
  stackDestroy (s) ;
  return answer ;
}

/*****************************************************************/

/*********** end of file *****************/
 
 
