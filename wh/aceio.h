/*
 *  File aceio.h  : header file for ACEDB input/output module
 *              (formerly known as freesubs/freeout)
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
 * HISTORY:
 * Last edited: Dec  2 13:33 2002 (edgrif)
 * * Oct 18 13:51 1999 (fw): aceInKeyMatch now private to acein.c
 * * Oct 18 13:45 1999 (fw): aceInSelect removed
 * * Oct  6 11:38 1999 (fw): changed aceInCreate interface to include handle
 * Created: Sept 1999 (mieg)
 * CVS info:   $Id: aceio.h,v 1.22 2019/05/24 21:46:06 mieg Exp $
 *-------------------------------------------------------------------
 */


#ifndef ACEDB_ACEIO_DEF
#define ACEDB_ACEIO_DEF

/* constructors */
ACEIN aceInCreateFromFile (const char *filename, const char *spec, const char *params,
			   AC_HANDLE handle);
ACEIN aceInCreateFromChooser (const char *prompt, 
			      char *directory, char *filename,
			      char *extension,
			      const char *spec, AC_HANDLE handle);
ACEIN aceInCreateFromPipe (char *command, const char *spec, const char *params,
			   AC_HANDLE handle);
ACEIN aceInCreateFromScriptPipe (char *script, const char *args, const char *params,
				 AC_HANDLE handle);
ACEIN aceInCreateFromStdin (BOOL isInteractive, const char *params,
			    AC_HANDLE handle);
          /* use stdinIsInteractive() to determine if STDIN is a tty */

ACEIN aceInCreateFromText (const char *text, const char *params,
			   AC_HANDLE handle);

/* simple reads from stdin if inFileNme == NULL 
 * Please ac_free(h) when done
 */
ACEIN aceInCreate (const char *inFileName, BOOL gzi, AC_HANDLE h) ;

/* expects a valid comma delimited fileList : f1,f2,... or u:f1,v:f2...
 * returns TRUE if all files can be opened
 * FALSE if not, and with a message in *error
 * Please ac_free(h) when done
 *
 * The idea is to check a full file list before running a program
 * to avoid processing 6 large file and crash on a typo in the 7th file name
 */
BOOL aceInCheckList (const char *fileList, const char **error, AC_HANDLE handle);

/* Loop on a file list
 * ACEIN ai = 0 ;
 * int nn = 0 ;
 * while (aceInCreateFromList (ai, &nn, "u:f1,v:f2,...", gzi, h) )
 *  { parse ai ; }
 * ac_free (h) ;
 *
 * expects a valid comma delimited fileList: f1,f2,...
 * or a list of labelled file names:  u:f1,v:f2
 * in which case u, v are stripped and stored in aceInLabel(ai)
 *
 * ai MUST BE initialised to NULL
 * it is freed inside aceInCreateFromList
 * or a file cannot be opened.
 * Please ac_free(h) when done
 * but DO NOT edit nn, DO NOT edit or free the returned ai. 
 * returns ai
 */
ACEIN aceInCreateFromList (ACEIN ai, int *np, const char *fileList, BOOL gzi, AC_HANDLE handle) ;

/* destructor */
void uAceInDestroy (ACEIN fi);	/* just a stub for messfree */
#define aceInDestroy(fi) (uAceInDestroy(fi),fi=0)

/* call with 1 to set, 0 to unset, any other number to interogate
 * if set, aceInCreate will retry to open the file for 1 minute
 * if unset (defaut) aceInCreate will try 3 times within 4 seconds 
 * then if the file is not found aceIn returns NULL
 * All this is useful on disks with heavy delays
 * so we recommand setting after initialisation
 * and before launching a big parallele code
 * in particular wego_init sets the waiting period to 1 minute
 */ 
int aceInWaitAndRetry (int retry) ;

/******/

char *aceInGetURL (ACEIN fi);

BOOL aceInPrompt (ACEIN fi, const char *prompt); /* lie aceInCard, but sets
					    * prompt for readline stdin */
BOOL aceInOptPrompt (ACEIN fi, int *keyp, FREEOPT *options);

BOOL aceInIsInteractive (ACEIN fi); /* is the user interactivly in 
				     * control of this stream ? */

/* set which characters get special treatment in aceIn-package
 * possible choices are "\n\t;/%\\@$"
 * e.g. if you leave out @ it won't be able to open command-files 
 * and therefore the @ character will be returned literally */
void aceInSpecial (ACEIN fi, const char* text) ;

/* Use these calls to "quote" and unquote text, this allows some control     */
/* over how text will be interpreted when its read in.                       */
char *aceInProtect (ACEIN fi, const char* text);    /* aceinword will read result back as text */
char *aceInUnprotect (ACEIN fi, const char *text) ;		    /* Remove all quoting:  \"\\  removed */
char *aceInUnprotectQuote(ACEIN fi, const char *text) ;	    /* Remove only quotes:  \"  removed  */

/* Read text/data from the acein.                                            */
char *aceInCard (ACEIN fi) ;	/* get a whole line */
void aceInForceCard (ACEIN fi, const char *string) ; /* fake a line */
void aceInCardBack (ACEIN fi) ;
char *aceInWord (ACEIN fi) ;	/* get next word as text */
BOOL aceInWordScan (ACEIN fi, char *buffer, int len) ;   
/* copy word of in buffer, 
 * copies len -1 char and garantees a terminal zero 
*/

char *aceInPath (ACEIN fi) ;	/* get word as path/filename */
BOOL aceInInt (ACEIN fi, int *p) ; /* get word as int */
BOOL aceInFloat (ACEIN fi, float *p) ; /* get word as float */
BOOL aceInDouble (ACEIN fi, double *p) ; /* get word as double */

char *aceInWordCut (ACEIN fi, const char *cutset, char *cutter) ;
void aceInBack (ACEIN fi) ;


BOOL aceInKey (ACEIN fi, KEY *keyp, FREEOPT *options) ;
BOOL aceInQuery (ACEIN fi, ACEOUT fo, const char *query) ;
BOOL aceInCheck (ACEIN fi, const char *fmt) ;	/* checks if current card 
					 * matches format */
BOOL aceInStep (ACEIN fi, char x) ;
void aceInNext (ACEIN fi) ;
char* aceInPos (ACEIN fi) ;

BOOL aceInBinary (ACEIN fi, void *buffer, long int size) ;

BOOL aceInEOF (ACEIN fi);	/* TRUE if fi->streamlevel == 0 */

const char *aceInFileName (ACEIN fi) ; /* name of the file */
const char *aceInLabel (ACEIN fi) ;    /* label of a member of a file list */
int aceInStreamLine (ACEIN fi) ; /* current line in file/text */
BOOL aceInStreamLength (ACEIN fi, long int *length); /* total file/text-size */
BOOL aceInStreamPos (ACEIN fi, long int *pos); /* current file/text-pos */

/* constructors */
ACEOUT aceOutCreateToFile (const char *filename, const char *spec, 
			   AC_HANDLE handle);
ACEOUT aceOutCreateToChooser (const char *prompt, 
			      char *directory, char *filename,
			      char *extension,
			      const char *spec, AC_HANDLE handle);
ACEOUT aceOutCreateToMail (const char *address, AC_HANDLE handle);
ACEOUT aceOutCreateToStdout (AC_HANDLE handle);
ACEOUT aceOutCreateToStderr (AC_HANDLE handle);
ACEOUT aceOutCreateToStack (Stack s, AC_HANDLE handle);
ACEOUT aceOutCreateToGzippedFile (const char *filename, AC_HANDLE handle) ;
ACEOUT aceOutCreateToPipe (const char *command, AC_HANDLE handle) ;

ACEOUT aceOutCopy (ACEOUT source_fo, AC_HANDLE handle); /* only useful for stdout right now. */
void aceOutSetNewStack(ACEOUT fo, Stack s) ;

/* simple outputs to stdout if outFileNme == NULL */
ACEOUT aceOutCreate (const char *outFileName, const char *suffix, BOOL gzo, AC_HANDLE h) ; 
void aceOutDate (ACEOUT ao, const char *commentPrefix, const char *title)  ; /* export "prefix title\tdate\tfile=filename",  commentPefix would usually be  "#" or "//" */

				      
/* destructor */
void uAceOutDestroy (ACEOUT fo) ;   /* will close file or send mail */
#define aceOutDestroy(fo) (uAceOutDestroy(fo),fo=0)

/******/
const char *aceOutFileName(ACEOUT fo) ;   /* Retrieve filename, defaults to "stdout" */

char *aceOutGetURL (ACEOUT fo);

BOOL aceOutFlush (ACEOUT fo);/* FALSE if fo isn't file-based or 
				 * the fflush on fo->fil failed */

BOOL aceOutRewind (ACEOUT fo);	/* FALSE if fo isn't file-based or 
				 * the fseek on fo->fil failed */

BOOL aceOutStreamPos (ACEOUT fo, long int *posp);
/* determines the file-position (ftell) or stack-pos if possible */

/* these functions return errno, which is ESUCCESS if all is OK
 * calling code should check for ENOSPC and others */
int aceOutBinary(ACEOUT fo, const void *data, long int size) ;
int aceOut(ACEOUT fo, const char *simple_string) ;	    /* Use for simple text string. */
int aceOutf(ACEOUT fo, char *format,...) ;		    /* Use for formatted printing. */
int aceOutxy(ACEOUT fo, const char *text, int x, int y) ; 
void aceOutPercent (ACEOUT ao, float z) ; /* export a precent, increasing the number of decimal near 0%, 50% 100% */
int aceOutInt (ACEOUT ao, int n) ; /* fast equivalent of aceOutf(ao,"%d", n) */
int aceOutLongInt (ACEOUT ao, int n) ; /* fast equivalent of aceOutf(ao,"%ld", n) */
int aceOutLongFloat (ACEOUT ao, float z) ; /* fast equivalent of aceOutf(ao,"%f", z) */


int aceOutLine (ACEOUT fo) ;	  /* how many lines written */
long int aceOutByte (ACEOUT fo) ; /* how many bytes have been written */
int aceOutPos (ACEOUT fo) ;

/****** ACETMP - writing to temporary files using ACEOUT ************/

/* Create a tmpfile in the system tmp directory.                             */
ACETMP aceTmpCreate (const char *spec, AC_HANDLE handle) ;

/* Create a tmpfile in the given directory if it is an absolute path, or in  */
/* a subdir of the system tmp directory if its a relative path.              */
/* Currently aceTmpCreateSubdir() will not create any intermediate subdirs   */
/* specified in the subdir parameter, they must already exist. If the actual */
/* subdir does not exist, it is created.                                     */
ACETMP aceTmpCreateDir(const char *dir, const char *spec, AC_HANDLE handle) ;

void uAceTmpDestroy (ACETMP atmp);
#define aceTmpDestroy(atmp) (uAceTmpDestroy(atmp),atmp=0)

/* By default aceTmp will delete its temporary file when its destroyed, you  */
/* can use this call to turn this on and off. If you turn it off you must    */
/* first retrieve the filename so you can later delete the file.             */
void aceTmpNoRemove(ACETMP atmp, BOOL no_remove) ;

/* The user has access to the internal output stream */
ACEOUT aceTmpGetOutput (ACETMP atmp);

/* close the tmp-file prematurely */
void aceTmpClose (ACETMP atmp);

/* The filename can still be enquired even when the 
 * output stream has already finished */
char *aceTmpGetFileName (ACETMP atmp);


#endif  /* !ACEDB_ACEIO_DEF */


