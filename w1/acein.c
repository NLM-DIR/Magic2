/*  File: acein.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
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
 *
 * Description: threadsafe version of freesubs
 * HISTORY:
 * Last edited: Apr  6 15:40 2001 (edgrif)
 * * May  8 13:12 2000 (edgrif): Added optimisation of aceInGetChar, see
 *              comments on aceInGetCharMacro() for further info.
 * * Jan  6 11:51 2000 (edgrif): Reinsert crash for exceeding MAXSTREAM.
 * * Nov  1 15:25 1999 (srk): added readline prompt
 * * Oct 18 13:35 1999 (fw): removed aceInSelect (unused)
 *              the proper function aceInLevelSelect
 * * Oct 18 10:48 1999 (fw): removed aceInRead/fi->filAss (unused)
 * * Oct 14 17:23 1999 (fw): rewrite constructor/destructor
 *              aceInSetFile/Pipe etc. now private
 * * Oct  6 13:34 1999 (fw): tidy up memory handling of AceInStruct
 *              add assertions
 * Created: Sun Aug 30 1999 (mieg)
 * CVS info:   $Id: acein.c,v 1.54 2020/06/12 20:20:11 mieg Exp $
 *-------------------------------------------------------------------
 */

#include <wh/regular.h>
#include <wh/aceio.h>		/* our public API */
#include <wh/version.h>
#include <wh/call.h>
#include <ctype.h>
#include <errno.h>
#include <wh/ac.h>
#include <stdlib.h>
#define NO_READLINE
#ifndef NO_READLINE
#include <readline/readline.h>
#endif /*NO_READLINE */

/* #define DEBUG_STREAMLEVEL */


/************************************************************/

enum {MAXSTREAM = 80,		/* highest streamlevel possible
				 * Number must be at least 1 */
      MAXNPAR = 80,
      MAXCARD = 1024,
      MAXSPECIAL = 24
} ;

typedef struct 
{
  FILE *fil ;
  const char *text ;
  const char *text_start ;					    /* Needed for streamPos on text stream */
  char *RLtext ;
  char *prompt ;
  char special[MAXSPECIAL] ;
  int npar ;
  int parMark[MAXNPAR] ;
  int line ;
  BOOL isPipe ;
  int streamlevel ;
  char *fname;			/* remember the file name, e.g. to debug */
} STREAM ;

static magic_t ACEIN_MAGIC = "ACEIN";

struct AceInStruct		/* complete the opaque type locally */
{
  magic_t *magic;		/* == &ACEIN_MAGIC */
  AC_HANDLE handle ;
  int streamlevel ;
  STREAM stream[MAXSTREAM + 1] ;
  int maxcard ;
  unsigned char special[256] ;
  const char *label ;
  char *card;
  char *pos;
  char *cardEnd;
  char *word ;
  Stack parStack ;
  Array buf ;
  int id;			/* count the number of created ACEIN's */
  BOOL debug;			/* to trace open/close */
  BOOL useReadLine;             /* private: using readline to get chars from STDIN */
  BOOL isInteractive;		/* public: clearing this inhibits use of readline
				   and stops prompts from being created
				   (set by the -noprompt option in tace 
				   and friends */
  FILE *curr_fil ;					    /* if curr level acein is a file, this */
							    /* is the FILE, NULL otherwise. */

};
 
/***********************/

#define _losewhite(fi)    while (*(fi->pos) == ' '|| (fi->special[(int) '\t'] && *(fi->pos) == '\t')) ++(fi->pos)
#define _stepover(fi,x)  (*(fi->pos) == x && ++(fi->pos))


/* 4_7 got input chars using a macro (_FREECHAR), this was replaced in 4_8   */
/* with a call to aceInGetChar(), this is too slow because it is called for  */
/* EVERY single char read. This has been corrected with two changes:         */
/*                                                                           */
/* - the below macro 'inlines' the aceInGetChar() code for files             */
/*                                                                           */
/* - curr_fil now gives direct access to the input file instead of having to */
/*   use massive indirection:  fi->stream[fi->streamlevel].fil)              */
/*                                                                           */
/* These changes brought 4_8 performance to within 8 or 9% of 4_7.           */
/*                                                                           */
#define aceInGetCharMacro(FI, RESULT)                   \
{                                                       \
  if (fi->curr_fil)                                     \
    {                                                   \
      unsigned char c = getc(fi->curr_fil);                       \
      if (c == (unsigned char)EOF)                                     \
	RESULT = '\0';                                  \
      else                                              \
	RESULT = (char)c;                               \
    }                                                   \
  else                                                  \
    RESULT =  aceInGetChar(FI) ;                        \
}

/***********************/

static ACEIN aceInDoCreate (AC_HANDLE handle);
static void aceInSetText (ACEIN fi, const char *string, const char *parms);
#ifndef NO_READLINE
static void aceInSetReadLine (ACEIN fi, char *parms);
#endif /*NO_READLINE */
static void aceInSetFile (ACEIN fi, FILE *fil, const char *fname, const char *parms);
static void aceInSetPipe (ACEIN fi, FILE *fil, const char *parms);
static void aceInClose(ACEIN fi, int level);
static void aceInExtend (ACEIN fi, char **pin) ;
static void aceInNewStream (ACEIN fi, const char *parms);
static void aceInFinalise (void *block);
static BOOL aceInExists (ACEIN fi);
static char aceInGetChar (ACEIN fi); /* get next char */
static void aceInUnGetChar (ACEIN fi, char ch); /* push back a char */
static int aceInKeyMatch (ACEIN fi, const char *cp, KEY *keyp, FREEOPT *options);
static ACEIN aceInCreateFromCompressedFile (const char *filename, const char *params, 
				     AC_HANDLE handle);

static char *doUnprotect(ACEIN fi, const char *text, BOOL just_quotes) ;

/******************************************************************
 ************************* public functions ***********************
 ******************************************************************/

int aceInWaitAndRetry (int retry)
{
  static int state = 0 ;
  switch (retry)
    {
    case 1: /* set retry */
      state = 1 ;
      break ;
    case 0: /* unset retry, very rare */
      state = 0 ;
      break ;
    }
  return state ;
}  /* aceInWaitAndRetry */

ACEIN aceInCreateFromFile (const char *filename, const char *spec, const char *params,
			   AC_HANDLE handle)
{
  FILE *fil = 0 ;
  ACEIN fi = NULL;

  if (!filename) messcrash("aceInCreateFromFile() - NULL filename");
  if (!spec) messcrash("aceInCreateFromFile() - NULL spec");
  if (spec[0] != 'r') messcrash("aceInCreateFromFile() - non-'r' spec");

  if (1)
    {
      int type = 0, n = strlen(filename);
      
      if (n > 4 && strcmp(filename+n-4, ".zip") == 0) type = 'z' ;
      if (n > 3 && strcmp(filename+n-3, ".gz") == 0) type = 'g' ;
      else if (n > 2 && strcmp(filename+n-2, ".Z") == 0) type = 'Z' ;
      
      if  (type)
       return aceInCreateFromCompressedFile (filename, params, handle);
    }
  if (1)
    /* uncompressed file */
    {
      int pass = aceInWaitAndRetry (-1) ?  6 : 2 ; /* up to 1 minute wait */
      int delay = 1 ;
      int failed = 0 ;

      while (! fil && pass--)
	{  /* loop in case of NFS delays */
	  fil = fopen (filename, (char*)spec);
	  if (!fil) fil = filopen (filename, "", (char*)spec);
	  if (fil)
	    break ;
	  sleep (delay) ;
	  delay *= 2 ; 
	  failed++ ;
	}

      if (fil)
	{
	  if (failed)
	    fprintf (stderr, "WARNING it took %d tries to open %s\n", failed + 1, filename) ;
	  fi = aceInDoCreate (handle);
	  
	  aceInSetFile (fi, fil, filename, params);	/* cannot fail */
	}
    }

  return fi;
} /* aceInCreateFromFile */


ACEIN aceInCreateFromURL (char *url, const char *spec, const char *params,
			  AC_HANDLE handle)
{
  ACEIN fi = NULL;
  
  if (!url) messcrash("aceInCreateFromURL() - NULL url");
  if (!spec) messcrash("aceInCreateFromURL() - NULL spec");
  if (spec[0] != 'r') messcrash("aceInCreateFromURL() - non-'r' spec");

  if (strcmp (url, "stdin://") == 0)
    {
      fi = aceInCreateFromStdin (FALSE, params, handle);
    }
  else if (strcmp (url, "readline://") == 0)
    {
      fi = aceInCreateFromStdin (TRUE, params, handle);
    }
  else if (strncmp (url, "text://",7) == 0)
    {
      messerror ("Cannot re-open an input that was "
		 "reading from a text-buffer");
    }
  else if (strncmp (url, "pipe://",7) == 0)
    {
      char *command;

      command = url + 7;
      fi = aceInCreateFromPipe (command, spec, params, handle);
    }
  else if (strncmp (url, "script://",9) == 0)
    {
      char *url_copy = strnew (url, 0);
      char *script;
      char *args;

      script = url_copy + 9;
      args = strstr (url_copy, "#");
      if (!args) messcrash("aceInCreateFromURL() - "
			   "script:// URL has no # character");
      *args = '\0';
      args++;

      fi = aceInCreateFromScriptPipe (script, args, params, handle);

      messfree (url_copy);
    }
  
  return fi;
} /* aceInCreateFromURL */


ACEIN aceInCreateFromChooser (const char *prompt, 
			      char *directory, char *filename,
			      char *extension,
			      const char *spec, AC_HANDLE handle)
     /* Note that extension is const char - it doesn't get changed
      * but the user can edit it and it is incoporated into the
      * filename-buffer */
{
  FILE *fil;
  ACEIN fi = NULL;
  char path[MAXPATHLEN];

  if (!directory) messcrash("aceInCreateFromChooser() - NULL directory");
  if (!filename) messcrash("aceInCreateFromChooser() - NULL filename");
  if (!spec) messcrash("aceInCreateFromChooser() - NULL spec");
  if (spec[0] != 'r') messcrash("aceInCreateFromChooser() - non-'r' spec");

  fil = filqueryopen (directory, filename, 
		      (char*)extension, (char*)spec, prompt);
  if (!fil)
    return NULL;		/* user clicked "Cancel" etc. */

  /* Assemble pathname :
   * Note, the assembled filename will be <directory>/<filename>
   * the filename will incorporate the extension, as it was
   * editable by the user. */
  sprintf (path, "%s%s%s", directory, SUBDIR_DELIMITER_STR, filename);
      
  fi = aceInCreateFromCompressedFile (path, "", handle);
  if (fi)
    {
      filclose (fil);		/* return value from 
				 * filechooser not needed */
    }
  else
    /* uncompressed file */
    {
      fi = aceInDoCreate (handle);
      
      aceInSetFile (fi, fil, path, ""); /* cannot fail */
    }

  return fi;
} /* aceInCreateFromChooser */


ACEIN aceInCreateFromPipe (char *command, const char *spec, const char *params,
			   AC_HANDLE handle)
{
  FILE *fil = 0 ;
  ACEIN fi = NULL;
      int pass = aceInWaitAndRetry (-1) ?  6 : 2 ; /* up to 1 minute wait */
  int delay = 1 ;
  int failed = 0 ;

  if (!command) messcrash("aceInCreateFromPipe() - NULL command");
  if (!spec) messcrash("aceInCreateFromPipe() - NULL spec");

     while (! fil && pass--)
    {
      fil = popen((const char*)command, spec); /* will be closed in 
						* aceInFinalise() */
      if (fil)
	break ;  
      sleep (delay) ;
      delay *= 2 ; 
      failed++ ;
    }

  if (fil)
    {
      if (failed)
	fprintf (stderr, "WARNING it took %d tries to open %s\n", failed + 1, command) ;
      fi = aceInDoCreate (handle);
      aceInSetPipe (fi, fil, params); /* cannot fail */

      {
	/* assemble url */
	Stack s = stackCreate(100); 
	pushText(s, "pipe://");
	catText(s, command);
	
	fi->stream[fi->streamlevel].fname = 
	  strnew(stackText(s, 0), fi->handle);
	fi->stream[fi->streamlevel].isPipe = TRUE ;

	stackDestroy (s);
      }
      if (fi->debug)
	printf ("DEBUG open stream %d(%d) : %s\n", 
		fi->id, fi->streamlevel, fi->stream[fi->streamlevel].fname);
    }

  return fi;
} /* aceInCreateFromPipe */


ACEIN aceInCreateFromScriptPipe (char *script, const char *args, const char *params,
				 AC_HANDLE handle)
{
  FILE *fil;
  ACEIN fi = NULL;

  if (!script) messcrash("aceInCreateFromScriptPipe() - NULL script");
  if (!args) messcrash("aceInCreateFromScriptPipe() - NULL args");

  fil = callScriptPipe(script, args); /* will be closed in 
				       * aceInFinalise() */

  if (fil)
    {
      fi = aceInDoCreate (handle);

      aceInSetPipe (fi, fil, params); /* cannot fail */

      {
	/* assemble url */
	Stack s = stackCreate(100);
	pushText(s, "script://");
	catText(s, script);
	catText(s, "#");
	catText(s, args);

	fi->stream[fi->streamlevel].fname = 
	  strnew (stackText(s, 0), fi->handle);
	stackDestroy (s);
      }

      if (fi->debug)
	printf ("DEBUG open stream %d(%d) : %s\n", 
		fi->id, fi->streamlevel, fi->stream[fi->streamlevel].fname);
    }

  return fi;
} /* aceInCreateFromScriptPipe */


ACEIN aceInCreateFromStdin (BOOL isInteractive, const char *params,
			    AC_HANDLE handle)
{
  /* Clearing isInteractive does two things.
     1) It stops prompts from being created.
     2) It inhibits the use of readline under all circumstances.
     
     If IsInteractive is true, readline is used iff stdin is a tty */
  
  ACEIN fi;
  
  fi = aceInDoCreate (handle);
  
  fi->isInteractive = isInteractive;
#ifndef NO_READLINE
  if (isInteractive && (isatty(fileno(stdin)) == 1))
    aceInSetReadLine (fi, params); /* cannot fail */
  else
#endif /*NO_READLINE */
    aceInSetFile (fi, stdin, "stdin://", params); /* cannot fail */

  return fi;
} /* aceInCreateFromStdin */


ACEIN aceInCreateFromText (const char *text, const char *params,
			   AC_HANDLE handle)
     /* NOTE: caller remains responsible for the text buffer */
{
  ACEIN fi = NULL;

  if (!text) messcrash("aceInCreateFromFile() - NULL text");

  fi = aceInDoCreate (handle);
  
  aceInSetText (fi, text, params); /* cannot fail */

  return fi;
} /* aceInCreateFromText */

/************************************************/


/* Reinitialise an acein struct with some text, use this if you want to reuse
 * an acein rather than constantly allocate/deallocate one. */
void aceInSetNewText(ACEIN fi, const char *text, const char *params)
{
  if (text != NULL || !aceInExists(fi))
    {
      
      aceInClose (fi, 1) ;					    /* close all levels */
    }

  aceInSetText(fi, text, params) ;

  return ;
}



/************************************************/

char *aceInGetURL (ACEIN fi)
{
  char *url = NULL;

  /* Return the filename of the first streamlevel on this stream.
   * If it has multiple levels, then we need to remember to 
   * file at the bottom of the stack to reporduce the same result.
   */
  if (!aceInEOF(fi))
    url = fi->stream[1].fname;

  return url;
} /* aceInGetURL */

/************************************************/

void uAceInDestroy (ACEIN fi)	/* use only aceInDestroy macro! */
{
  if (!aceInExists(fi))
    messcrash("aceInDestroy() - received invalid fi pointer");

  /* just a stub - you could just call messfree
   * cleanup is taken care of by finalisation */
  messfree (fi);
} /* uAceInDestroy */

/*******************/

BOOL aceInPrompt (ACEIN fi, const char *prompt)
{
  if (fi->stream[fi->streamlevel].prompt)
    messfree(fi->stream[fi->streamlevel].prompt);

  fi->stream[fi->streamlevel].prompt = strnew (prompt, fi->handle);
  
  
  /* NB. if isInteractive is clear, we don't want a prompt.
     if we are using readline readline  generates the prompt itself so
     we don't need to do it here.
  */
  if (fi->isInteractive && !fi->useReadLine && 
      fi->stream[fi->streamlevel].prompt &&
      *fi->stream[fi->streamlevel].prompt
      )
    fputs(fi->stream[fi->streamlevel].prompt, stdout);  /* fputs so no \n */
  aceInCard (fi);

  if (aceInEOF(fi))
    return FALSE;

  return TRUE;
} /* aceInPrompt */

BOOL aceInOptPrompt (ACEIN fi, int *keyp, FREEOPT *options)
     /* Using a prompt we get a selection from the user
      * fi - (only used if text-I/O)
      *       input used to get the choice from the user (e.g. STDIN)
      *
      * kpt - pointer in which the choice is returned, 
      *         refers to the option's KEY in the FREEOPT-options
      * options - FREEOPT array of possible choices
      *           The text of the first entry is the prompt-text
      *           Any option that would have the same effect as
      *           an EOF on the input should have the option-key (-1)
      *
      * RETURNS:
      *   TRUE - if a choice was made. The choice is returned
      *          in *keyp.
      *   FALSE - choice was ambiguous or didn't match anything in 
      *          options. *keyp is left untouched.
      */
{
  BOOL answer;

  messfree(fi->stream[fi->streamlevel].prompt);

  fi->stream[fi->streamlevel].prompt = (char*) halloc(strlen(options[0].text)+3,
					      fi->handle);
  strcpy(fi->stream[fi->streamlevel].prompt, options[0].text);
  strcat(fi->stream[fi->streamlevel].prompt, "> ");
  
  /* NB. if isInteractive is clear, we don't want a prompt.
     if we are using readline readline  generates the prompt itself so
     we don't need to do it here.
  */
  if (fi->isInteractive && !fi->useReadLine &&
      fi->stream[fi->streamlevel].prompt &&
      *fi->stream[fi->streamlevel].prompt
      )
    fputs(fi->stream[fi->streamlevel].prompt,stdout ); /* fputs so no \n */

  aceInCard (fi);
  
  if (aceInEOF(fi))
    /* input stream is finished (e.g. EOF character on STDIN) */
    { 
      *keyp = -1;  
      /* we return true to recognise the EOF as an explicit
       * command. In this case we set the option to -1
       * so the calling code can recognise this as the EOF
       * command, just like an explicit quit-command. */
      answer = TRUE;                       
    }
  else
    {
      KEY key;
      /* Get the next word in the card and match it against the options
       * Note: if there's no word in the card it'll return FALSE */
      answer = aceInKey (fi, &key, options) ;
      if (answer) *keyp = (int)key;
    }

  return answer;
} /* aceInOptPrompt */

/************/


/* Returns TRUE if acein is a tty and caller did _not_ set interactive FALSE,*/
/* returns FALSE otherwise.                                                  */
/*                                                                           */
/* Hence returns TRUE when:                                                  */
/*   1) readline is being used, readline can only be used with a tty).       */
/*   2) readline is not used but acein is a tty and the caller did not set   */
/*      interactive false.                                                   */
/* The corollary of this is that pipes will _not_ be interactive, even if    */
/* the acein was created with interactive TRUE, this is IMPORTANT for calling*/
/* tace from scripts, see tacemain.c where it decides whether to save the    */
/* database.                                                                 */
/*                                                                           */
BOOL aceInIsInteractive (ACEIN fi)
{
  BOOL interactive ;

  if (!aceInExists(fi))
    messcrash("aceInIsInteractive() - received invalid fi pointer");

  if (fi->useReadLine)
    interactive = TRUE ;
  else if (fi->isInteractive && fi->curr_fil && (isatty(fileno(fi->curr_fil)) == 1))
    interactive = TRUE ;
  else
    interactive = FALSE ;

  return interactive ;
} /* aceInIsInteractive */

/*******************/

BOOL aceInEOF (ACEIN fi)
     /* TRUE if we're at the end of the input */
     /* Notice that stream level maybe 0, but there is still one card
      * left, this happens when fi is text-based for instance.  */
{
  BOOL isEOF = FALSE;

  if (!aceInExists(fi))
    messcrash("aceInEOF() - received invalid fi pointer");

  if (fi->card == NULL)
    /* set to NULL by aceInCard, if streamlevel reaches 0 */
    isEOF = TRUE;
  else if (fi->useReadLine && fi->streamlevel == 0)
    /* The ZERO streamlevel was generated by receiving
     * EOF on an interactive STDIN. We have to query useReadLine
     * because that's decided upon whether it is a real tty-device
     * rather than what the user decides */
    /* This clause is important in deciding when the interactive 
     * stream ends, because fi->card wouldn't be NULL until 
     * aceInCard is called the next time round */
    isEOF = TRUE;

  return isEOF;
} /* aceInEOF */

/*******************/

void aceInSpecial (ACEIN fi, const char* text)
     /* The text contains all characters that have special meaning in aceInCard().
      *
      * Possible choices are "\n\t\";/%\\@$"
      *    The following actions will be performed upon inclusion
      *    of the characters in the special-text for a particular stream.
      * \n - NEWLINE will terminate card (when parsing line-based text-files).
      * \t - TAB will be expanded to 8 spaces
      * \" - protect text until next \"
      * ;  - NEWLINE can be used a line-break, if multiline input is given on one command-line
      * /  - COMMENT text after // until the enf-of-line will be ignored
      * %  - PARAMETER, e.g. %1 will be replaced by the first parameter
      * \\ - if last char on line - used to break up single line input over multiple lines
      *      and when "\n" (2 chars) should be treated as '\n' (1 newline char)
      */
{
  if (!aceInExists(fi))
    messcrash("aceInSpecial() - received invalid fi pointer");
  if (!text)
    messcrash ("aceInSpecial() - received NULL text") ;
  if (strlen(text) > 23)
    messcrash ("aceInSpecial() - received a string longer than 23") ;

  /* if (text != fi->stream[fi->streamlevel].special)  mieg 2011_09_26 it makes no sense to compare the address */
    strcpy (fi->stream[fi->streamlevel].special, text) ;
  memset (fi->special, 0, (mysize_t) 256) ;
  while (*text)
    fi->special [((int) *text++) & 0xFF] = TRUE ;
  fi->special[0] = TRUE ;     /* ensures EOF recognition */
  
  return;
} /* aceInSpecial */
 
/********************/
 
void aceInForceCard (ACEIN fi, const char *string)
{ 
  if (!aceInExists(fi))
    messcrash("aceInForceCard() - received invalid fi pointer");

  aceInSetText (fi, string, "") ;
  aceInCard (fi) ;

  return;
} /* aceInForceCard */
 
/********************/
 
char* aceInCard (ACEIN fi)
     /* returns NULL when EOF is reached */
{ 
  char *in;
  char ch;
  char *cp ;
  int kpar ;
  FILE *fil ;
  BOOL acceptShell, acceptCommand ;

  if (!aceInExists(fi))
    messcrash("aceInCard() - received invalid fi pointer");

restart :
  if (fi->streamlevel == 0)
    { 
      fi->card = 0;	/* important to recognise aceInEOF() */
      return 0 ;
    }

  in = fi->card ; --in ;

  acceptCommand = fi->special[(int)'@'] ? TRUE : FALSE  ;
  acceptShell = fi->special[(int)'$'] ? TRUE : FALSE  ;
 
  while (TRUE)
    { if (++in >= fi->cardEnd)
	aceInExtend (fi, &in) ;

      aceInGetCharMacro(fi, *in)
      /* returns 0 if end-of-text or end-of-file */

    lao:
      if (fi->special[((int) *in) & 0xFF] && *in != '$' && *in != '@' )
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

	  case '\0':
	    if (fi->stream[fi->streamlevel].fil &&
		fi->stream[fi->streamlevel].fil != stdin &&
		! fi->stream[fi->streamlevel].isPipe)
	      {
		long nn1, nn2 ;
		nn1 = ftell (fi->stream[fi->streamlevel].fil) ;
		fseek (fi->stream[fi->streamlevel].fil, 0, SEEK_END) ;
		nn2 = ftell (fi->stream[fi->streamlevel].fil) ;
		if (nn1 != nn2)
		  messerror ("File %s line %d aceincard hit an EOF at %ld before true EOF %ld"
			     , fi->stream[fi->streamlevel].fname ? fi->stream[fi->streamlevel].fname : "" 
			     , fi->stream[fi->streamlevel].line
			     , nn1, nn2) ;
	      }
	    
	    aceInClose(fi, fi->streamlevel) ;
	    /*
	     * the code returns an empty card for the last empty line of a file endding as usual on  \n
	     * it could seem preferable to return NULL, but that would change a 20 year old behaviour
	    if (! *fi->card) // mieg, 2010_01_06 
	      { 
		fi->card = 0 ;
		return 0 ;
	      }
	    */
	    goto got_line;

	  case '\t':     /* tabs should get rounded to 8 spaces */
            *in++ = ' ' ;
            while ((in - fi->card) % 8)
              { if (in >= fi->cardEnd)
                  aceInExtend (fi, &in) ;
                *in++ = ' ' ;
              }
	    --in ;
	    continue ;  

	  case '/':		/* // means start of comment */
            aceInGetCharMacro(fi, ch) ;
	    if (ch == '/')
	      { 
		do {		/* skip the rest of this line */
		  aceInGetCharMacro(fi, ch) ;
		} while (ch != '\n' && ch != '\0');

		goto got_line ;	/* it'll return "" - an empty string */
	      }
	    else
	      { 
		aceInUnGetChar (fi, ch);
		/*
		if (fi->stream[fi->streamlevel].fil)
		  ungetc (ch, fi->stream[fi->streamlevel].fil) ;
		else
		  --(fi->stream[fi->streamlevel].text) ;
		*/
	      }
	    break ;

	  case '%':		/* possible parameter */
	    --in ; kpar = 0 ;
	    while ( isdigit ((int)(ch = aceInGetChar(fi))) )
	      kpar = kpar*10 + (ch - '0') ;
	    if (kpar > 0 && kpar <= fi->stream[fi->streamlevel].npar)
	      for (cp = stackText (fi->parStack, 
				   fi->stream[fi->streamlevel].parMark[kpar-1]) ; *cp ; ++cp)
		{ if (++in >= fi->cardEnd)
		    aceInExtend (fi, &in) ;
		  *in = *cp ;
		}
	    else
	      messout ("Parameter %%%d cannot be substituted", kpar) ;
	    if (++in >= fi->cardEnd)
	      aceInExtend (fi, &in) ;
	    *in = ch ; 
	    goto lao ; /* mieg */

	  case '\\':		/* escapes next character - interprets \n */
	    *in = aceInGetChar(fi) ;
	    if (*in == '\n')    /* fold continuation lines */
	      { 
		while ((ch = aceInGetChar(fi)) == ' ' || ch == '\t') ;
			/* remove whitespace at start of next line */

		aceInUnGetChar (fi, ch);
		/*
		if (fi->stream[fi->streamlevel].fil)
		  ungetc (ch, fi->stream[fi->streamlevel].fil) ;
		else
		  --(fi->stream[fi->streamlevel].text) ;
		*/

		fi->stream[fi->streamlevel].line++ ;
		--in ;
	      }
#ifndef WIN32
	    else if (*in == 'n') /* reinterpret \n as a format */
	      {
		*in = '\n' ; 
	      }
#endif /* !WIN32 */
	    else  /* keep the \ till aceinword is called */
	      { 
		*(in+1) = *in ;
		*in = '\\' ;
		if (++in >= fi->cardEnd)
		  aceInExtend (fi, &in) ;
	      }
	    break ;

	  default:
	    messerror ("aceinsubs got unrecognised special character 0x%x = %c\n",
		       *in, *in) ;
	  }
      else
	{ if (iscntrl((int)*in) && *in != '\t' && *in != '\n') /* mieg 2015_10_12, replaced isprint by ! isctrl to allow parsing Göndör, Millán-Ariño, Németi and other accentuated characters */
	    --in ;
	}
    } /* while TRUE loop */
 
 got_line:
  fi->stream[fi->streamlevel].line++ ;
  *in = 0 ;
  fi->pos = fi->card ;
  _losewhite (fi) ;

  if (acceptCommand && _stepover (fi, '@'))        /* command file */
    { 
      char *name ;

      if ((name = aceInWord (fi))
	  && (fil = filopen (name, "", "r")))
	aceInSetFile (fi, fil, name, fi->pos) ;

      goto restart;
    }

  /* mieg - ?? can we accept shells in thread safe way ?? */
  if (acceptShell && _stepover (fi, '$'))        /* shell command */
    {
      callSystem ((char*)fi->pos) ;
      goto restart ;
    }
  if (fi->debug)
    fprintf (stderr, "aceInCard returns #%s#\n",  fi->card ?  fi->card : "zero") ; 

  return fi->card ;
} /* aceInCard */
 
/************************************************/

void aceInCardBack (ACEIN fi)    /* goes back one card */
{ 
  if (!aceInExists(fi))
    messcrash("aceInCardBack() - received invalid fi pointer");

  fi->stream[fi->streamlevel].line-- ;
  aceInSetText (fi, fi->card, "") ;
  fi->stream[fi->streamlevel].line = fi->stream[fi->streamlevel - 1].line ;
}

/************************************************/
 
int aceInStreamLine (ACEIN fi)
{ 
  if (!aceInExists(fi))
    messcrash("aceInStreamLine() - received invalid fi pointer");

  return fi->stream[fi->streamlevel].line  ;
}

/************************************************/
 
const char * aceInFileName (ACEIN fi)
{ 
  if (!aceInExists(fi))
    messcrash("aceInStreamFileName() - received invalid fi pointer");

  return fi->stream[fi->streamlevel].fname ;
}

/************************************************/
 
const char * aceInLabel (ACEIN fi)
{ 
  if (!aceInExists(fi))
    messcrash("aceInLabel() - received invalid fi pointer");

  return fi->label ;
}

/************************************************/
 
BOOL aceInStreamLength (ACEIN fi, long int *length)
{
  long int oldPos, endPos;

  if (!aceInExists(fi))
    messcrash("aceInStreamLength() - received invalid fi pointer");
  if (!length)
    messcrash("aceInStreamLength() - received NULL length pointer");

  if (fi->stream[fi->streamlevel].isPipe)
    /* can't determine length of streams from pipes */
    return FALSE;

  if (fi->stream[fi->streamlevel].text)
    /* stream comes from text-buffer, report its size */
    {
      *length = (long int)strlen(fi->stream[fi->streamlevel].text);
      return TRUE;
    }

  if (fi->stream[fi->streamlevel].fil == NULL)
    return FALSE;

  if (fi->stream[fi->streamlevel].fil == stdin)
    return FALSE;

  oldPos = ftell (fi->stream[fi->streamlevel].fil);
  if (oldPos == -1)
    return FALSE;		/* we should report ERRNO */

  if (fseek (fi->stream[fi->streamlevel].fil, 0L, SEEK_END) == -1)
    return FALSE;		/* we should report ERRNO */

  endPos = ftell (fi->stream[fi->streamlevel].fil);
  if (endPos == -1)
    return FALSE;		/* we should report ERRNO */

  if (fseek (fi->stream[fi->streamlevel].fil, oldPos, SEEK_SET) == -1)
    return FALSE;		/* we should report ERRNO */

  *length = endPos;
  return TRUE;
} /* aceInStreamLength */


BOOL aceInStreamPos (ACEIN fi, long int *pos)
{
  long int filPos;

  if (!aceInExists(fi))
    messcrash("aceInStreamLength() - received invalid fi pointer");
  if (!pos)
    messcrash("aceInStreamLength() - received NULL pos pointer");

  if (fi->stream[fi->streamlevel].isPipe)
    /* can't determine file-pos of streams from pipes */
    return FALSE;

  if (fi->stream[fi->streamlevel].text)
    {
      *pos = (long int)(fi->stream[fi->streamlevel].text - fi->stream[fi->streamlevel].text_start);
      return TRUE;
    }

  if (fi->stream[fi->streamlevel].fil == NULL)
    return FALSE;

  if (fi->stream[fi->streamlevel].fil == stdin)
    return FALSE;

  filPos = ftell (fi->stream[fi->streamlevel].fil);
  if (filPos == -1)
    return FALSE;		/* we should report ERRNO */

  *pos = filPos;

  return TRUE;
} /* aceInStreamPos */


/********************************************/
/* aceinword(), aceinwordcut() and aceinstep() are the only calls that
     directly move pos forward -- all others act via aceinword().
   aceinback() moves pos back one word.
*/
 
char *aceInWord (ACEIN fi)
{
  char *cw ;
 
  if (!aceInExists(fi))
    messcrash("aceInWord() - received invalid fi pointer");

  _losewhite (fi) ;         /* needed in case of intervening aceinstep() */

  if (fi->special[(int) '"'] && _stepover (fi, '"'))
    { 
      for (cw = fi->word ; !_stepover(fi, '"') && *(fi->pos) ; *cw++ = *(fi->pos)++)
	if (fi->special[(int) '\\'] && _stepover(fi, '\\'))	/* accept next char unless end of line */
	  if (*fi->pos == '\0')
	    break ;
      _losewhite (fi) ;
      *cw = 0 ;
      return (char*) fi->word ;	/* always return a word, even if empty */
    }
  
  /* default: break on space and \t, not on comma */
  for (cw = fi->word ; isgraph ((int)*(fi->pos)) && *(fi->pos) != '\t' ; 
       *cw++ = *(fi->pos)++)
    if (fi->special[(int) '\\'] && _stepover(fi, '\\'))	/* accept next char unless end of line */
      if (!*(fi->pos))
	break ;
  _losewhite (fi) ;
  *cw = 0 ;

  if (*(fi->word))
    return (char *)fi->word;

  return (char*)NULL;
}  /* aceInWord */
 
/************************************************/
/* cp word of in buffer, cp len -1 char and add a zero */
BOOL aceInWordScan (ACEIN fi, char *buffer, int len)
{
  char *cp = aceInWord (fi) ;
  int n = cp ? strlen (cp) : 0 ;

  if (!buffer || len < 0)
    messcrash ("bad call to aceInWordScan len=%d"
	       , len, buffer ? "" : ", buffer == 0") ;
  if (n > len - 1) n = len - 1;
  if (n >= 0) { strncpy (buffer, cp, n) ;   buffer[n] = 0 ; }
  return cp ? TRUE : FALSE ;
}

/************************************************/

char *aceInPath (ACEIN fi)
{
#ifdef __CYGWIN__
  char *cw ;
#endif /* __CYGWIN__ */

  if (!aceInExists(fi))
    messcrash("aceInPath() - received invalid fi pointer");

#ifdef __CYGWIN__
  _losewhite (fi) ;             /* needed in case of intervening aceinstep() */


  if (_stepover (fi, '"'))
    { for (cw = fi->word ; !_stepover(fi, '"') && *(fi->pos) ; *cw++ = *(fi->pos)++)
	if (_stepover(fi, '\\')) /* accept next char unless end of line */
	  if (!*(fi->pos))
	    break ;
      _losewhite(fi) ;
      *cw = 0 ;
      return (char*) fi->word ;	/* always return a word, even if empty */
    }

  /* default: break on space, \t or end of line, not on comma
	 also, does not skip over backslashes which are assumed to be
	 MS DOS/Windows path delimiters */
  for (cw = fi->word ; ( *(fi->pos) == '\\' || isgraph (*(fi->pos)) ) && *(fi->pos) != '\t' ; *cw++ = *(fi->pos)++) ;
  _losewhite(fi) ;
  *cw = 0 ;
#else  /* !__CYGWIN__*/
  aceInWord(fi) ;
#endif /* !__CYGWIN__ */


  if (*fi->word)
    {
      char *cp = filGetFullPath (fi->word, 0);
      if (cp)
	{
	  strncpy (fi->word, cp, fi->maxcard);
	  fi->word[fi->maxcard-1] = '\0';
	  messfree(cp);
	}

      return (char*)fi->word;
    }

  return (char*)NULL;
} /* aceInPath */

 
/************************************************/
 
char *aceInWordCut (ACEIN fi, const char *cutset, char *cutter)
     /* Moves along card, looking for a character from cut, which is a
      * 0-terminated char list of separators.
      * Returns everything up to but not including the first match.
      * pos is moved one char beyond the character.
      * *cutter contains the char found, or if end of card is reached, 0.
      */
{ 
  const char *cc ; char *cw ;
  
  if (!aceInExists(fi))
    messcrash("aceInWordCut() - received invalid fi pointer");

  for (cw = fi->word ; *fi->pos ; *cw++ = *fi->pos++)
    for (cc = cutset ; *cc ; ++cc)
      if (*cc == *fi->pos)
        goto wcut ;
wcut:
  if (cutter) *cutter = *fi->pos ;
  if (*fi->pos)
    ++(fi->pos) ;
  _losewhite (fi) ; 
  *cw-- = 0 ;
  /* lose trailing spaces in the word */
  while (cw >= fi->word && *cw==' ')*cw--=0;

  if (*fi->word)
    return (char*)fi->word;

  return (char*)NULL;
} /* aceInWordCut  */
 
/************************************************/
 
void aceInBack (ACEIN fi)
     /* goes back one word - inefficient but reliable */
{
  char *now = fi->pos ;
  char *old = fi->pos ;
  
  if (!aceInExists(fi))
    messcrash("aceInBack() - received invalid fi pointer");

  fi->pos = fi->card ; _losewhite (fi) ;
  while  (fi->pos < now)
    {
      old = fi->pos ;
     aceInWord (fi) ;
    }
  fi->pos = old ;

  return;
}  /* aceInBack */ 

/************************************************/
/* Read a word representing an int from wherever acein is pointing to.        */
/*                                                                           */
/* If there is no word OR the word cannot be converted to an int, reset      */
/*   the acein pos, don't set the int param. and return FALSE                 */
/* If the word is "NULL" set int param to the POSIX "too small" int value    */
/*   and return TRUE                                                         */
/* Otherwise set the int param to the converted int and return TRUE          */
/*                                                                           */
/* Note that valid range of ints is    INT_MIN < valid < INT_MAX             */
/* otherwise UT_NON_INT doesn't work.                                        */
/*                                                                           */
BOOL aceInInt (ACEIN fi, int *p)
{
  BOOL result = FALSE ;
  char *keep ;
  enum {DECIMAL_BASE = 10} ;
  char *endptr ;
  long int bigint ;
  
  if (!aceInExists(fi))
    messcrash("aceInInt() - received invalid fi pointer");

  keep = fi->pos ;
  if (aceInWord (fi))
    { /*printf ("aceInInt got '%s'\n", word) ;*/
      if (strcmp((const char*)fi->word, "NULL") == 0)
	{ 
	  *p = UT_NON_INT ;
	  result = TRUE ;
	}
      else
	{
	  errno = 0 ;
	  bigint = strtol((char *)fi->word, &endptr, DECIMAL_BASE) ;
	
	  if ((bigint == 0 && endptr == fi->word)	    /* first character wrong */
	      || (endptr != fi->word + strlen((const char*)fi->word)) /* some other character wrong */
	      || (errno == ERANGE)			    /* number too small/big for long int */
	      || (bigint <= INT_MIN || bigint >= INT_MAX))    /* number too small/big for int */
	    {
	      fi->pos = keep ;
	      return FALSE ;
	    }
	  else
	    {
	      *p = (int)bigint ;
	      result = TRUE ;
	    }
	}
    }
  else
    {
      fi->pos = keep ;
      result = FALSE ;
    }
  
  return result ;
} /* aceinInt */
 
/*****************************/
/* Read a word representing a float from wherever acein is pointing to.       */
/*                                                                           */
/* If there is no word OR the word cannot be converted to a float, reset     */
/*   the acein pos, don't set the float param. and return FALSE               */
/* If the word is "NULL" set float param to the POSIX "too small" float value*/
/*   and return TRUE                                                         */
/* Otherwise set the float param to the converted float and return TRUE      */
/*                                                                           */
/* Note that valid range of floats is:                                       */
/*                             -ve  -FLT_MAX < valid <  -FLT_MIN             */
/*                             +ve   FLT_MIN < valid <   FLT_MAX             */
/* otherwise UT_NON_FLOAT doesn't work as a range check for applications.    */
/*                                                                           */
BOOL aceInFloat (ACEIN fi, float *p)
{
  BOOL result = FALSE ;  
  char *keep ;
  double bigfloat = 0 ;
  char *endptr ;

  if (!aceInExists(fi))
    messcrash("aceInFloat() - received invalid fi pointer");

  keep = fi->pos ;
  if (aceInWord (fi))
    {
      if (strcmp ((const char*)fi->word, "NULL") == 0)
	{
	  *p = UT_NON_FLOAT ;				    /* UT_NON_FLOAT = -FLT_MAX */
	  result = TRUE ;
	}
      else
	{
	  errno = 0 ;
	  bigfloat = (double)strtod((const char*)fi->word, &endptr) ;

	  if ((bigfloat == +0.0 && endptr == fi->word)  /* first character wrong */
	      || (endptr != fi->word + strlen((const char*)fi->word)) /* some other character wrong */
	      || (errno == ERANGE)			    /* number too small/big for double */
	      || (bigfloat < 0 && (bigfloat >= -FLT_MIN || bigfloat <= -FLT_MAX))
	      || (bigfloat > 0 && (bigfloat <= FLT_MIN || bigfloat >= FLT_MAX)))
	    {						    /* number too small/big for float */
	      fi->pos = keep ;
	      result = FALSE ;
	    }
	  else
	    {
	      *p = (float)bigfloat ;
	      result = TRUE ;
	    }
	}
    }
  else
    {
      fi->pos = keep ;
      result = FALSE ;
    }

  return result ;
}  /* aceinFloat  */
 
/**************************************************/
/* Read a word representing a double from wherever acein is pointing to.      */
/*                                                                           */
/* If there is no word OR the word cannot be converted to a double, reset    */
/*   the acein pos, don't set the double param. and return FALSE              */
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
BOOL aceInDouble (ACEIN fi, double *p)
{ 
  BOOL result = FALSE ;  
  char *keep ;
  double bigfloat ;
  char *endptr ;

  if (!aceInExists(fi))
    messcrash("aceInDouble() - received invalid fi pointer");

  keep = fi->pos ;
  if (aceInWord(fi))
    {
      if (strcmp ((const char*)fi->word, "NULL") == 0)
	{
	  *p = UT_NON_DOUBLE ;		 /* UT_NON_DOUBLE = -DBL_MAX */
	  result = TRUE ;
	}
      else
	{
	  errno = 0 ;
	  bigfloat = strtod((const char*)fi->word, &endptr) ;

	  if ((bigfloat == +0.0 && endptr == fi->word)  /* first character wrong */
	      || endptr != fi->word + strlen((const char*)fi->word) /* some other character wrong */
	      || (errno == ERANGE)			    /* number too small/big for double */
	      || (bigfloat < 0 && (bigfloat == -DBL_MIN || bigfloat == -DBL_MAX))
	      || (bigfloat > 0 && (bigfloat == DBL_MIN || bigfloat == DBL_MAX)))
	    {						    /* number too small/big. */
	      fi->pos = keep ;
	      result = FALSE ;
	    }
	  else
	    {
	      *p = bigfloat ;
	      result = TRUE ;
	    }
	}
    }
  else
    {
      fi->pos = keep ;
      result = FALSE ;
    }

  return result ;
} /* aceInDouble */

/*************************************************/
 
BOOL aceInKey (ACEIN fi, KEY *keyp, FREEOPT *options)
{
  char *keep;
  int n = 0 ;

  if (!aceInExists(fi))
    messcrash("aceInKey() - received invalid fi pointer");

  keep = fi->pos ;

  if (!aceInWord(fi))
    return FALSE ;

  n = aceInKeyMatch (fi, fi->word, keyp, options) ;
  switch (n)
    {
    case 2:
      return TRUE;

    case 1:
      messout ("Keyword %s is ambiguous\n", fi->word) ;
      break ;

    case 0:
      messout ("Keyword %s does not match\n", fi->word) ;
      break ;
    }
  fi->pos = keep ;

  return FALSE ;
} /* aceInKey */
 
/*************************************************/
 
BOOL aceInBinary (ACEIN fi, void *buffer, long int size)
{
  long int num_bytes_read = 0 ;

  if (!aceInExists(fi) || ! fi->curr_fil)
    messcrash("aceInKey() - received invalid fi pointer");

  if (fi->curr_fil)
    {
      num_bytes_read = fread (buffer,1, size,fi->curr_fil);
      
      if (ferror (fi->curr_fil) != 0)
	return FALSE ; /* we should report errno */
      if (num_bytes_read != size)
	return FALSE ;
    }

  return TRUE ;
} /* aceInKey */
 
/*****************/

  /* Return the text corresponding to the key */
char *aceInKey2Text (KEY k, FREEOPT *o)  
{ 
  int i = o->key ; char *title = o->text ;

  if (i<0)
    messcrash("aceInKey2Text() - Negative number of options") ;
  while (o++, i--) 
    if (o->key == k)
      return (o->text) ;
  return title ;
} /* aceInKey2Text */

/***************************************************/
 
BOOL aceInQuery (ACEIN fi, ACEOUT fo, const char *query)
{
  if (!aceInExists(fi))
    messcrash("aceInQuery() - received invalid fi pointer");

  if (fi->debug)
    fprintf (stderr,"aceInQuery %s\n", query) ;
  if (aceInIsInteractive(fi))
    { 
      BOOL retval;
      char answer;

      aceOutf (fo, "%s (y or n) ",query) ;
      answer = aceInGetChar(fi) ;
      retval = (answer == 'y' || answer == 'Y') ? TRUE : FALSE ;
      while (answer != '\0' && answer != '\n')
        answer =  aceInGetChar(fi) ;

      return retval ;
    }

  return TRUE ;
}  /* aceInQuery */
 
/*************************************/
 
int aceInFmtLength (char *fmt)
{
  char *cp ;
  int length = 0 ;
  
  if (isdigit((int)*fmt))
    {
      sscanf (fmt,"%d",&length) ;
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
	  messcrash ("'o' can not end acein format %s",fmt) ;
	length += 2 ; break ;
      }
  
  if (!length)
    length = 40 ;
  return length ;
}  /* aceInFmtlength */

/****************/
 
BOOL aceInCheck (ACEIN fi, const char *fmt)
        /* checks that whatever is in card fits specified format
           note that 't' format option changes card by inserting a '"' */
{
  char *keep;
  union {int i ; float r ; double d ;}
  target ;
  const char *fp ;
  char *start ;
  int nquote = 1 ;
  
  if (!aceInExists(fi))
    messcrash("aceInCheck() - received invalid fi pointer");

  keep = fi->pos;

  for (fp = fmt ; *fp ; ++fp)
    switch (*fp)
      {
      case 'w' : if (aceInWord (fi)) break ; else goto retFALSE ;
      case 'i' : if (aceInInt (fi, &target.i)) break ; else goto retFALSE ;
      case 'f' : if (aceInFloat (fi, &target.r)) break ; else goto retFALSE ;
      case 'd' : if (aceInDouble (fi, &target.d)) break ; else goto retFALSE ;
      case 't' :      /* must insert '"' and escape any remaining '"'s or '\'s */
	for (start = fi->pos ; *fi->pos ; ++(fi->pos))
	  if (*(fi->pos) == '"' || *(fi->pos) == '\\')
	    ++nquote ;
	*(fi->pos + nquote + 1) = '"' ;		/* end of line */
	for ( ; fi->pos >= start ; --(fi->pos))
	  { *(fi->pos + nquote) = *(fi->pos) ;
	  if (*(fi->pos) == '"' || *(fi->pos) == '\\')
	    *(fi->pos + --nquote) = '\\' ;
	  }
	*start = '"' ;
	goto retTRUE ;
      case 'z' : if (aceInWord (fi)) goto retFALSE ; else goto retTRUE ;
      case 'o' :
	if (!*++fp) messcrash ("'o' can not end acein format %s",fmt) ;
	aceInStep (fi, *fp) ; break ;
      case 'b' : break; /* special for graphToggleEditor no check needed  il */
      default :
	if (!isdigit((int)*fp) && !isspace((int)*fp))
	  messerror ("unrecognised char %d = %c in acein format %s",
		     *fp, *fp, fmt) ;
      }
  
retTRUE:
  fi->pos = keep ; return TRUE ;
retFALSE:
  fi->pos = keep ; return FALSE ;
} /* aceInCheck */

/************************ little routines ************************/
 
BOOL aceInStep (ACEIN fi, char x)
{
  if (!aceInExists(fi))
    messcrash("aceInStep() - received invalid fi pointer");

  return (*(fi->pos) && ace_upper (*(fi->pos)) == x && (fi->pos)++) ? TRUE : FALSE  ;
}
 
void aceInNext (ACEIN fi)
{ _losewhite (fi) ;
}

char* aceInPos (ACEIN fi)
     /* cheat to give pos onwards - i.e. returns rest-of-line */
{
  if (!aceInExists(fi))
    messcrash("aceInPos() - received invalid fi pointer");

  return (char*) fi->pos ;
} /* aceInPos */

/*************************************************************************/
/*************************************************************************/

/* It's sometimes necessary to quote parts of the text so it doesn't get     */
/* interpreted too literally (i.e. you may want to keep \n as "\n"). The     */
/* aceInProtect() routine will do this. But then you will probably need to   */
/* unquote the text at some time.                                            */
/*                                                                           */
/* Quoting protects text by putting \" at the start and end of the text and  */
/* \\ in front of any special chars.                                         */
/*    - sometimes you want to completely reverse this - use aceInUnprotect() */
/*       to do this.                                                         */
/*    - sometimes you just want to remove the leading and trailing \", use   */
/*       aceInUnprotectQuote() to do this.                                   */
/*                                                                           */
char* aceInUnprotect(ACEIN fi, const char *text)
{
  char *result = NULL ;

  result = doUnprotect(fi, text, FALSE) ;

  return result ;
}

char *aceInUnprotectQuote(ACEIN fi, const char *text)
{
  char *result = NULL ;

  result = doUnprotect(fi, text, TRUE) ;

  return result ;
}

static char *doUnprotect(ACEIN fi, const char *text, BOOL just_quotes)
{
  char *ac_unprotectquote (const char *text, AC_HANDLE h) ;

  if (!aceInExists(fi))
    messcrash("aceInUnprotect() - received invalid fi pointer");

  
  if (!aceInExists(fi))
    messcrash("aceInProtect() - received invalid fi pointer");
  if (just_quotes)
    return ac_unprotectquote (text, fi->handle) ;
  else
    return ac_unprotect (text, fi->handle) ;
 } /* aceInUnprotect */

/****************************************/

char* aceInProtect (ACEIN fi, const char* text)	/* aceinword will read result back as text */
{
  if (!aceInExists(fi))
    messcrash("aceInProtect() - received invalid fi pointer");
  return ac_protect (text, fi->handle) ;
}  /* aceInProtect */

/****************************************/

char* aceInJavaProtect (ACEIN fi, char* text)
{
  Array a;
  char *cp, *cq ;
  int n = 2*(1+strlen(text)) ;
 
  if (!aceInExists(fi))
    messcrash("aceInJavaProtect() - received invalid fi pointer");

  a = fi->buf;
  if (a)
    a = arrayReCreate (a, n + 1, char) ;
  else
    a = fi->buf = arrayHandleCreate (n + 1, char, fi->handle) ;
  array (a, n, char) = 0 ; /* ensure long enough */
  
  cq = arrp (a, 0, char) ;
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

  return arrp (a, 0, char) ;
} /* aceInJavaProtect */

/******************************************************************
 ************************ private functions ***********************
 ******************************************************************/

static ACEIN aceInDoCreate (AC_HANDLE handle)
     /* generic constructor isn't public - we can only make an
      * ACEIN object from a file/pipe or the like */
{ 
  static int total_num_created = 0; /* for debug */
  ACEIN fi;

  fi = (ACEIN)halloc (sizeof (struct AceInStruct), handle);
  blockSetFinalise (fi, aceInFinalise);

  fi->magic = &ACEIN_MAGIC;
  fi->handle = handleCreate();	/* killed in aceInFinalise */
  fi->id = ++total_num_created;
  fi->streamlevel = 0 ;
  fi->stream[fi->streamlevel].fil = NULL ;
  fi->stream[fi->streamlevel].text = fi->stream[fi->streamlevel].text_start = NULL ;
  fi->stream[fi->streamlevel].RLtext = NULL ;
  fi->stream[fi->streamlevel].prompt = NULL ;
  fi->maxcard = MAXCARD ;
  aceInSpecial (fi, "\n\t\"\\/") ;
  fi->card = (char *)halloc (fi->maxcard, fi->handle) ;
  fi->cardEnd = fi->card + fi->maxcard - 1 ;
  fi->pos = fi->card ;
  fi->word = (char *)halloc (fi->maxcard, fi->handle) ;
  fi->parStack = stackHandleCreate (128, fi->handle) ;
  fi->useReadLine = FALSE;
  fi->isInteractive = FALSE;
  fi->curr_fil = NULL ;

#ifdef DEBUG_STREAMLEVEL
  fi->debug = TRUE;
#else  /* !DEBUG_STREAMLEVEL */
  fi->debug = FALSE;
#endif  /* !DEBUG_STREAMLEVEL */
  
  return fi ;
} /* aceInDoCreate */


static void aceInSetText (ACEIN fi, const char *string, const char *parms)
     /* switch fi to parse from a text-buffer */
{
  if (!aceInExists(fi))
    messcrash("aceInSetText() - received invalid fi pointer");

  aceInNewStream (fi, parms) ;

  fi->stream[fi->streamlevel].text = fi->stream[fi->streamlevel].text_start = string ;
  fi->stream[fi->streamlevel].fname = "text://" ;

  if (fi->debug)
    printf ("DEBUG open stream %d(%d) : %s\n", 
	    fi->id, fi->streamlevel, fi->stream[fi->streamlevel].fname);

  fi->curr_fil = NULL ;

  return ;
} /* aceInSetText */

#ifndef NO_READLINE
static void aceInSetReadLine (ACEIN fi, const char *parms)
{
  if (!aceInExists(fi))
    messcrash("aceInSetReadLine() - received invalid fi pointer");

  aceInNewStream (fi, parms) ;

  fi->stream[fi->streamlevel].text = "";
  fi->stream[fi->streamlevel].fname = strnew("readline://", fi->handle);
  fi->useReadLine = TRUE;

  if (fi->debug)
    printf ("DEBUG open stream %d(%d) : %s\n", 
	    fi->id, fi->streamlevel, fi->stream[fi->streamlevel].fname);

  fi->curr_fil = NULL ;

  return ;
} /* aceInSetReadLine */
#endif /*NO_READLINE */

static void aceInSetFile (ACEIN fi, FILE *fil, const char *fname, const char *parms)
     /* switch fi to read from the contents of a file */
{
  if (!aceInExists(fi))
    messcrash("aceInSetFile() - received invalid fi pointer");
  if (!fil) messcrash("aceInSetFile() - NULL fil");
  if (!fname) messcrash("aceInSetFile() - NULL fname");

  aceInNewStream (fi, parms) ;

  fi->stream[fi->streamlevel].fil = fil ;
  fi->stream[fi->streamlevel].fname = strnew(fname, fi->handle);

  if (fi->debug)
    printf ("DEBUG open stream %d(%d) : %s\n", 
	    fi->id, fi->streamlevel, fi->stream[fi->streamlevel].fname);

  fi->curr_fil = fil ;

  return ;
} /* aceInSetFile */

static void aceInSetPipe (ACEIN fi, FILE *fil, const char *parms)
     /* switch fi to read contents of a pipe - parse program results etc. */
{
  if (!aceInExists(fi))
    messcrash("aceInSetPipe() - received invalid fi pointer");
  if (!fil) messcrash("aceInSetPipe() - NULL fil");

  aceInNewStream (fi, parms) ;

  fi->stream[fi->streamlevel].fil = fil ;
  fi->stream[fi->streamlevel].isPipe = TRUE ;

  fi->curr_fil = fil ;

  return ;
} /* aceInSetPipe */

/*******************/

static void aceInFinalise (void *block)
     /* called when an ACEIN struct is cleaned up, either it is being
      * messfree'd explicitly or its parent handle is being destroyed */
{
  ACEIN fi = (ACEIN)block;

  if (!aceInExists(fi))
    messcrash("aceInFinalise() - received invalid block pointer");

  aceInClose (fi, 1);		/* close all levels */

  /* the handle contains all its associated stacks, buffers etc. */
  messfree(fi->handle);
  fi->magic = 0;		/* taint this memory area as freed
				 * future access to this memory block
				 * will trap aceInExists-assertions */
  return ;
} /* aceInFinalise */


static BOOL aceInExists (ACEIN fi)
     /* verify the validity of an ACEIN pointer */
{
  if (fi && fi->magic == &ACEIN_MAGIC)
    return TRUE;

  return FALSE;
} /* aceInExists */

static void aceInNewStream (ACEIN fi, const char *parms)
{
  int kpar ;

  if (fi->streamlevel == MAXSTREAM)
    messcrash("Maximum number of streams (%d) exceeded in ACEIN",  MAXSTREAM) ;
  else
    {
      fi->streamlevel++;

      memset (&(fi->stream[fi->streamlevel]), 0, sizeof(STREAM)) ;

      strcpy (fi->stream[fi->streamlevel].special, fi->stream[fi->streamlevel-1].special) ;

      fi->stream[fi->streamlevel].npar = 0 ;
      fi->stream[fi->streamlevel].line = 1 ;

      if (parms && strlen(parms) > 0)
	{
	  fi->pos = strnew (parms, 0) ; /* abuse aceinword() to get parms */
     
	  for (kpar = 0 ; kpar < MAXNPAR && aceInWord (fi) ; kpar++) /* read parameters */
	    {
	      fi->stream[fi->streamlevel].parMark[kpar] = stackMark (fi->parStack);
	      pushText (fi->parStack, fi->word) ;
	    }
     
	  fi->stream[fi->streamlevel].npar = kpar ;
	  fi->stream[fi->streamlevel].isPipe = FALSE ;
	  messfree (fi->pos) ;
	  fi->pos = fi->card ;			/* restore pos to start of blank card */
	  *(fi->card) = 0 ;
	}
    }

  return ;
} /* aceInNewStream */
 
/******************************/

static void aceInClose (ACEIN fi, int level)
{ 
  int kpar ;

  if (!aceInExists(fi))
    messcrash("aceInClose() - received invalid fi pointer");

  while (fi->streamlevel >= level)
    {
      if (fi->debug)
	printf ("DEBUG close stream %d(%d) : %s%s\n",
		fi->id, fi->streamlevel,
		fi->stream[fi->streamlevel].isPipe ? "PIPE " : "",
		fi->stream[fi->streamlevel].fname);


      if (fi->stream[fi->streamlevel].fil)
	{
	  if (fi->stream[fi->streamlevel].isPipe)
	    pclose (fi->stream[fi->streamlevel].fil) ;
	  else
	    filclose (fi->stream[fi->streamlevel].fil) ;
	}
      
      if (fi->stream[fi->streamlevel].RLtext)
	free(fi->stream[fi->streamlevel].RLtext);

      for (kpar = fi->stream[fi->streamlevel].npar ; kpar-- ;)
	popText (fi->parStack) ;

      --(fi->streamlevel) ;
      
      if (fi->streamlevel == 0 && !fi->stream[0].line) 
	fi->stream[0].line = fi->stream[1].line ;
    }

  /* Must reset for current level of acein                                   */
  if (fi->stream[fi->streamlevel].fil)
    fi->curr_fil = fi->stream[fi->streamlevel].fil ;
  else
    fi->curr_fil = NULL ;

  return;
} /* aceInClose */

/******************************/

static void aceInExtend (ACEIN fi, char **pin)
{	/* only happens when getting card */
  char *oldCard = fi->card ;
  long int offset = *pin - fi->card ;
  
  fi->maxcard *= 2 ;
  fi->card = (char *)halloc (fi->maxcard, fi->handle);
  if (oldCard)     /* jtm june 22, 1992 */
    memcpy (fi->card, oldCard, fi->maxcard/2) ;
  fi->cardEnd = fi->card + fi->maxcard - 1 ;
  *pin += (int)(fi->card - oldCard) ;
  *pin = fi->card + offset ;
  messfree (oldCard) ;
  messfree (fi->word) ;
  fi->word = (char *)halloc (fi->maxcard, fi->handle);

  return;
} /* aceInExtend */
 
static char aceInGetChar (ACEIN fi)
     /* read another character from 'fi'.
      * We have to make sure only to return valid 'char'-values.
      * therefore EOF has to be cast to NULL-character, because
      * EOF is typically -1, which is not a valid 'char'-value.
      */
{
  char ch;

  if (fi->curr_fil)
    {
      int c = getc(fi->curr_fil);
      if (c == EOF)
	ch = '\0';
      else
	ch = (char)c;
    }
  else if (*(fi->stream[fi->streamlevel].text))
    {
      ch = *(fi->stream[fi->streamlevel].text++);
    }
#ifndef NO_READLINE
  else if (fi->useReadLine)
    {
      if (fi->stream[fi->streamlevel].RLtext)
	{
	  /* reached on on readline: return /n and set up for new call */
	  free(fi->stream[fi->streamlevel].RLtext);
	  fi->stream[fi->streamlevel].RLtext = 0;
	  fi->stream[fi->streamlevel].text = "";
	  ch = '\n';
	}
      else
	{
	  /* new call to readline */
	  fi->stream[fi->streamlevel].text = 
	    fi->stream[fi->streamlevel].RLtext = 
	    readline(fi->stream[fi->streamlevel].prompt ? 
		     fi->stream[fi->streamlevel].prompt : "");
	  if (!fi->stream[fi->streamlevel].RLtext) /* EOF typed */
	    {
	      fi->stream[fi->streamlevel].text = "";
	      ch = '\0';
	    }
	  else if (!(*(fi->stream[fi->streamlevel].RLtext))) /* empty line typed */
	    { 
	      /* make sure we get called straight away again */
	      fi->stream[fi->streamlevel].text = ""; 
	      free(fi->stream[fi->streamlevel].RLtext);
	      fi->stream[fi->streamlevel].RLtext = 0;
	      ch = '\n';
	    }
	  else
	    { 
	      add_history(fi->stream[fi->streamlevel].text);
	      ch = *(fi->stream[fi->streamlevel].text++);
	    }
	}
    }
#endif /*NO_READLINE */
  else
    ch = '\0';
  
  return ch;
} /* aceInGetChar */
 
static void aceInUnGetChar (ACEIN fi, char ch)
     /* push back a character */
{
  if (fi->stream[fi->streamlevel].fil)
    /* NOTE: only 4 bytes of pushback are guaranteed */
    ungetc ((int)ch, fi->stream[fi->streamlevel].fil) ;
  else
    --(fi->stream[fi->streamlevel].text) ;

  return;
} /* aceInUnGetChar */


/* This routine compares *cp to the first word (i.e. the keyword) of each
 * option->text string. The routine can return the follwing values:
 * 
 *    0 : no match at all.
 *    1 : ambiguous, *cp matches more then one options keyword.
 *    2 : exact between *cp and an options keyword or just one inexact match,
 *        in this case *keyp is filled with the options->key for that match.
 */

#define KEY_UNDEFINED 0

static int aceInKeyMatch (ACEIN fi, const char *cp, KEY *keyp, FREEOPT *options)
{
  const char  *io, *iw ;
  int   nopt = (int)options->key ;
  KEY   key ;
  BOOL finished ;
  int match ;
 
  if (!aceInExists(fi) || !cp  || !keyp || !nopt)
    messcrash("aceInKeyMatch() - received bad parameters: %s%s%s%s",
	      !aceInExists(fi) ? "invalid fi pointer" : "",
	      !cp ? " NULL *cp pointer" : "",
	      !keyp ? " NULL *keyp pointer" : "",
	      !nopt ? " empty FREEOPT array" : "") ;

  key = KEY_UNDEFINED ;
  match = 0 ;						    /* default is no match. */
  finished = FALSE ;
  while (!finished)
    {
      iw = cp ;
      io = (++options)->text ;

      while (ace_upper (*iw++) == ace_upper(*io++))
	{
	  if (!*iw)
	    {
	      if (*io != ' ')
		{
		  if (key == KEY_UNDEFINED)
		    {
		      match = 2 ;			    /* possible single inexact match. */
		      key = options->key ;
		    }
		  else
		    {
		      match = 1 ;			    /* More than one inexact match. */
		    }
		}
	      else
		{
		  match = 2 ;				    /* exact match. */
		  finished = TRUE ;
		  key = options->key ;
		}
	      break ;
	    }
	}

      if (!--nopt)
	finished = TRUE ;
    }

  /* Only return key for exact match or SINGLE inexact match. */
  if (match == 2)
    *keyp = key ;

  return match ;
}

 
/************************************************************/

static ACEIN aceInCreateFromCompressedFile (const char *filename, const char *params, 
				     AC_HANDLE handle)
{
  int n;
  char type = '\0';
  ACEIN fi = NULL;
  char args[MAXPATHLEN+10];

  if (!filename) messcrash("aceInCreateFromCompressedFile() - NULL filename");

  /* try to open close the file, to vanquish tcp delays */
  if (1)
    {
      FILE *fil = 0 ;
      int pass = aceInWaitAndRetry (-1) ?  6 : 2 ; /* up to 1 minute wait */ 
      int delay = 1 ;
      int failed = 0 ;

      
      while (! fil && pass--)
	{  /* loop in case of NFS delays */
	  fil = fopen (filename, "rb") ;
	  if (!fil) fil = filopen (filename, "", "rb") ;
	  if (!fil)
	    break ;  
	  sleep (delay) ;
	  delay *= 2 ; 
	  failed++ ;
	}

      if (fil)
	{
	  if (failed)
	    fprintf (stderr, "WARNING it took %d tries to open %s\n", failed + 1, filename) ;
 	  fclose (fil) ;
	}
      else
	return NULL ;
    }


  n = strlen(filename);

  if (n > 4 && strcmp(filename+n-4, ".zip") == 0) type = 'z' ;
  if (n > 3 && strcmp(filename+n-3, ".gz") == 0) type = 'g' ;
  else if (n > 2 && strcmp(filename+n-2, ".Z") == 0) type = 'Z' ;

  switch (type)
    {
    case 'z':			/* .zip - winzip */
      sprintf (args, " -p %s", filename);
      fi = aceInCreateFromScriptPipe ("unzip", args, params, handle);
      if (!fi)
	fi = aceInCreateFromScriptPipe ("pkunzip", args, params, handle);
      break;

    case 'g':			/* .gz - gzip */
      sprintf (args, " -fdc %s", filename);
      fi = aceInCreateFromScriptPipe ("gzip", args, params, handle);
      if (!fi)
	{
	  sprintf (args, " -fc %s", filename);
	  fi = aceInCreateFromScriptPipe ("gunzip", args, params, handle);
	}
      break;

    case 'Z':			/* .Z - compress */
      sprintf (args, " -fdc %s", filename);
      fi = aceInCreateFromScriptPipe ("gzip", args, params, handle);
      if (!fi)
	{
	  sprintf (args, " -c %s", filename);
	  fi = aceInCreateFromScriptPipe ("gunzip", args, params, handle);
	}
      if (!fi)
	{
	  sprintf (args, " %s", filename);
	  fi = aceInCreateFromScriptPipe ("zcat", args, params, handle);
	}
      break;

    default:;			/* not a compressed file */
    }

  /* even if we recognised a compression extension, this may still
   * be NULL if the scriptPipe failed */
  return fi;
} /* aceInCreateFromCompressedFile */

/*************************************************************************************/
/*****************************    Utilities     **************************************/
/*************************************************************************************/

ACEIN aceInCreate (const char *inFileName, BOOL gzi, AC_HANDLE h)
{
  ACEIN ai = 0 ;
  char *label = 0, *cr = 0 ;

  cr = (inFileName ?  strchr (inFileName, ':') : 0) ;

  if (cr)
    {
      label = strnew (inFileName, h) ;
      cr = strchr (label, ':') ;
      *cr++ = 0 ;
      inFileName = cr ;
    }
  if (inFileName && inFileName[0] == '<')
    {
      ai = aceInCreateFromURL (hprintf(h,"pipe://%s",inFileName+1), "r", 0, h) ;
      if (ai)
	ai->label = label ;
      return ai ; 
    }
  if (inFileName && strcmp (inFileName, "-"))
    {
      char *rtype = "r" ;
      const char *ccp ;
      BOOL isBam = FALSE ;

      if (!gzi &&
	  (ccp = strstr (inFileName, ".gz")) &&
	  ! strcmp (ccp, ".gz")
	  )
	gzi = TRUE ;

      if (!gzi &&
	  (ccp = strstr (inFileName, ".bam")) &&
	  ! strcmp (ccp, ".bam")
	  )
	isBam = TRUE ;

      if (*inFileName == '/' || *inFileName == '.' || filCheckName (inFileName, 0, "r"))
	{
	  if (isBam)
	    ai = aceInCreateFromPipe (hprintf (h, "samtools view %s", inFileName), rtype, 0, h) ;
	  else if (gzi)
	    ai = aceInCreateFromPipe (hprintf (h, "gzip -dc %s", inFileName), rtype, 0, h) ;
	  else
	    ai = aceInCreateFromFile (inFileName, rtype, 0, h) ;
	}
      if (! ai)
	messerror ("file %s not found", inFileName) ;
    }
  else
    {
      char *rtype = "r" ;
      if (gzi)
	ai = aceInCreateFromPipe ("gunzip -c -", rtype, 0, h) ;
      else
	ai = aceInCreateFromStdin (0,0,h) ;
      inFileName = "stdin" ;
    }
  if (ai)
    {
      ai->label = label ;
      aceInSpecial (ai, "\n\\\"") ;
    }

  return ai ;
} /* aceInCreate */

/*************************************************************************************/
/* expects a valid comma delimited fileList : f1,f2,...
 * returns TRUE if all files can be opened
 * FALSE if not, and with a message in *error
 * Please ac_free(h) when done
 */
BOOL aceInCheckList (const char *fileList, const char **error, AC_HANDLE handle)
{
  BOOL ok = TRUE ;
  AC_HANDLE h = ac_new_handle () ;   
  char *cp, *cq, *cr, *list = strnew (fileList, h) ;
  int nn = 0 ;

  cp = cq = list ;
  while (ok && cp)
    {
      nn++ ;
      cq = strchr (cp, ',') ;
      if (cq) *cq++ = 0 ;

      if (0)
	{
	  cr = strchr (cp, ':') ;
	  if (cr)
	    {
	      *cr++ = 0 ;
	      cp = cr ;
	    }
	}
      if (! filCheckName (cp, 0, "r"))
	{
	  *error = hprintf (handle, "ERROR cannot open file %d : %s in list\n %s\n", nn, cp, fileList) ;
	  ok = FALSE ;
	  break ;
	}
      cp = cq ;
    }

  ac_free (h) ;
  return ok ;
} /* aceInCheckList */

/*************************************************************************************/
/* Loop on a file list
 * ACEIN ai = 0 ;
 * int nn = 0 ;
 * while (aceInCreateFromList (ai, &nn, "f1,f2,...", gzi, h) )
 *  { parse ai ; }
 * ac_free (h) ;
 *
 * expects a valid comma delimited fileList : f1,f2,...
 * ai MUST BE initialised to NULL
 * it is freed inside aceInCreateFromList
 * or a file cannot be opened.
 * Please ac_free(h) when done
 * but DO NOT edit nn, DO NOT ac_free the returned ai. 
 * returns ai
 */
ACEIN aceInCreateFromList (ACEIN ai, int *np, const char *fileList, BOOL gzi, AC_HANDLE handle)
{ 
  AC_HANDLE h = ac_new_handle () ;   
  char *cp, *cq, *cr, *label, *list = strnew (fileList, h) ;
  int nn = 0 ;
  
  if (! np)
    messcrash ("aceInCreateFromList received invalid &nn, fileList=%s, &nn must exist", fileList) ;
  if (ai)
    {
      if (ai->magic != &ACEIN_MAGIC)
	messcrash ("aceInCreateFromList received invalid ai (nn == %d), fileList=%s, please initialise ai=0 before first call and do not edit nn or ai", *np, fileList) ;
      ac_free (ai) ;
    }
  else
    *np = 0 ;
  
  cp = cq = list ;
  while (cp)
    {
      cq = strchr (cp, ',') ;
      if (cq) *cq++ = 0 ;

      if (0)
	{
	  cr = strchr (cp, ':') ;
	  if (cr)
	    {
	      *cr++ = 0 ;
	      label = strnew (cp, handle) ;
	      cp = cr ;
	    }
	}

      if (*np == nn)
	{
	  *np = nn + 1 ; /* for next call */
	  if (! (ai = aceInCreate (cp, gzi, handle)))
	    messcrash ("ERROR cannot open file %d : %s in list %s\n", nn, cp, fileList) ;
	  if (0) ai->label = label ;
	  ac_free (h) ;
	  return ai ;
	}
      nn++ ;
      cp = cq ;
    }
  ac_free (h) ;
  return 0 ;
} /* aceInCreateFromList */

/*************************************************************************************/
/*************************************************************************************/

