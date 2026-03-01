/*
 *  File regular.h  : header file for libfree.a utility functions
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
 * Last edited: Sep 12 09:36 2001 (edgrif)
 * * May 23 12:06 2000 (edgrif): Fix SANgc07771
 * * Apr 26 10:07 1999 (edgrif): Added UTIL_NON_DOUBLE for freedouble check.
 * * Mar 22 14:47 1999 (edgrif): Added messSetMsgInfo, getSystemName calls,
 *              removed messErrorInit.
 * * Mar 18 10:53 1999 (edgrif): Moved some includes to mystdlib.h & put
 *              POSIX constants for invalid ints/floats.
 * * Jan 25 15:56 1999 (edgrif): Added getLogin from session.h
 * * Jan 21 15:27 1999 (edgrif): Added UtArraySize macro to determine array
 *              size at compile time.
 * * Jan 13 11:50 1999 (edgrif): Interface changed for messGetErrorProgram().
 * * Sep  9 16:54 1998 (edgrif): Add messErrorInit decl.
 * * Sep  9 14:31 1998 (edgrif): Add filGetFileName decl.
 * * Aug 20 11:50 1998 (rbrusk): AUL_FUNC_DCL
 * * Sep  3 11:50 1998 (edgrif): Add macro version of messcrash to give
 *              file/line info for debugging.
 * Created: 1991 (rd)
 *-------------------------------------------------------------------
 */

#ifndef DEF_REGULAR_H
#define DEF_REGULAR_H



/* #define MEM_DEBUG  */
/* #define MALLOC_CHECK  */
/* #define ARRAY_CHECK */

				/* library EXPORT/IMPORT symbols */
#if defined (WIN32)
#include "win32libspec.h"  /* must come before mystdlib.h...*/
#else
#define UTIL_VAR_DCL	extern
#define UTIL_VAR_DEF
#endif

#include <limits.h>
#include <strings.h>
#include <stdio.h>
#include <float.h> 
#include <stdlib.h>
#include "mystdlib.h"    /* contains full prototypes of system calls */
#include "ctype.h"

#if defined(WIN32)
#if defined(_DEBUG)
#define MEM_DEBUG /* must be defined here, before acelibspec.h */
#include <crtdbg.h>
#endif
UTIL_VAR_DCL char* linkDate ;
UTIL_VAR_DCL int isInteractive ;      /* can set FALSE, i.e. in tace */
#endif

#ifdef FALSE
  typedef int BOOL ;
#else
  typedef enum {FALSE=0,TRUE=1} BOOL ;
#endif

typedef struct _AC_HANDLE_STRUCT *AC_HANDLE ; /* opaque outside memsubs.c */
#include "aceiotypes.h"

typedef unsigned char UCHAR ; /* for convenience */

typedef unsigned int KEY ;

typedef void (*VoidRoutine)(void) ;
typedef void (*Arg1Routine)(void *arg1) ;

/* magic_t : the type that all magic symbols are declared of.
   They become magic (i.e. unique) by using the pointer
   to that unique symbol, which has been placed somewhere
   in the address space by the compiler */
/* type-magics and associator codes are defined at
   magic_t MYTYPE_MAGIC = "MYTYPE";
   The address of the string is then used as the unique 
   identifier (as type->magic or graphAssXxx-code), and the
   string can be used during debugging */
typedef char* magic_t;

/* Use this macro to determine the size of a declared array.                 */
#define UtArraySize(ARRAY) ((unsigned int) (sizeof(ARRAY) / sizeof(ARRAY[0])))

/* Values that mean "not a valid number", returned by freeint/float/double.  */
/* We use POSIX values to set these, meaning that for floats and ints we get */
/* a slightly reduced range of valid ints/floats.                            */
#define UT_NON_INT      INT_MIN
#define UT_NON_FLOAT   -FLT_MAX
#define UT_NON_DOUBLE  -DBL_MAX
#define UT_NON_BOOL (2)		/* always test b & TRUE , b & NON_BOOL */
#define UT_NON_FLAG (~0)		/* all bits - same as in flag.h */


/* Values used for rounding in a machine independant way the acedb floats */
/* the rounding is maintained using  sprintf (buffer, "%1.7g", obj->curr->n.f) */
#define ACE_FLT_MIN  1e-40 
#define ACE_FLT_RESOLUTION .25e-12

typedef struct freestruct
  { KEY  key ;
    char *text ;
  } FREEOPT ;

/*---------------------------------------------------------------------*/
/* The free package for reading from files/stdout, see freesubs.c      */
/*                                                                     */

void freeinit (void) ;
int freeCurrLevel(void) ;			    /* Returns current level. */
char* freecard (int level) ;	/* 0 if below level (returned by freeset*) */
void freecardback (void) ;  /* goes back one card */
void freeforcecard (const char *string);
int  freesettext (const char *string, const char *parms) ; /* returns level to be used in freecard () */
int  freesetfile (FILE *fil, const char *parms) ;
int freesetnamedfile (char *filName, FILE *fil, const char *parms) ;
int  freesetpipe (FILE *fil, const char *parms) ;  /* will call pclose */
void freeclose(int level) ; /* closes the above */
void freespecial (char *set) ;	/* set of chars to be recognized from "\n;/%\\@$" */
BOOL freeread (FILE *fil) ;	/* returns FALSE if EOF */
int  freeline (FILE *fil) ;	/* line number in file */
int  freestreamline (int level) ;/* line number in stream(level)*/
char *freeword (void) ;
void freeshutdown(void) ;

#if defined(WIN32)  /* A variation to correctly parse MS DOS/Windows pathnames */
  char *freepath (void) ;
#else	/* NOT defined(WIN32) */
#define freepath freeword  /* freeword() works fine if not in WIN32 */
#endif	/* defined(WIN32) */

/* always move to the next cutter, return TRUE if an int or float was scanned */
BOOL freeintcut (float *p, char *cutset, char *cutter) ;
BOOL freefloatcut (float *p, char *cutset, char *cutter) ;
BOOL freedoublecut (double *p, char *cutset, char *cutter) ;
char *freewordcut (char *cutset, char *cutter) ;

void freeback (void) ;		/* goes back one word */
BOOL freeint (int *p) ;
BOOL freefloat (float *p) ;
BOOL freedouble (double *p) ;

BOOL freekey (KEY *kpt, FREEOPT *options) ;
BOOL freekeymatch (char *text, KEY *kpt, FREEOPT *options) ;
void freemenu (void (*proc)(KEY), FREEOPT *options) ;
char *freekey2text (KEY k, FREEOPT *o)  ;  /* Return text corresponding to key */
BOOL freeselect (KEY *kpt, FREEOPT *options) ;
BOOL freelevelselect (int level,
				    KEY *kpt, FREEOPT *options);
void freedump (FREEOPT *options) ;
BOOL freestep (char x) ;
void freenext (void) ;
BOOL freeprompt (const char *prompt, char *dfault, char *fmt) ;/* gets a card */
BOOL freecheck (char *fmt) ;	/* checks remaining card fits fmt */
int  freefmtlength (char *fmt) ;
BOOL freequery (const char *query) ;
char *freepos (void) ;		/* pointer to present position in card */
char *freeprotect (const char* text) ; /* protect so freeword() reads correctly */
char* freeunprotect (const char *text) ; /* reverse of protect, removes \ etc */

extern const char UT_UPPER[128] ;
#define ace_upper(x)	(UT_UPPER[(x) & 0x7f])  /* table is only 128 long */

extern const char UT_LOWER[128] ;
#define ace_lower(x)	(UT_LOWER[(x) & 0x7f])

void utUnlimitResources(BOOL allow_user_abort) ;
int utPrintfSizeOfArgList (const char * formatDescription, va_list marker) ;
BOOL getCmdLineOption (int *argcp, const char **argv,
		       const char *arg_name, const char **arg_val) ; /* in utils.c */
BOOL getCmdLineText (AC_HANDLE h, int *argcp, const char **argv,
		       const char *arg_name, const char **arg_val) ; /* in utils.c */
BOOL getCmdLineBool (int *argcp, const char **argv, const char *arg_name) ;
BOOL getCmdLineInt (int *argcp, const char **argv, const char *arg_name, int *val) ;
BOOL getCmdLineLong (int *argcp, const char **argv, const char *arg_name, long int *val) ; 
BOOL getCmdLineFloat (int *argcp, const char **argv, const char *arg_name, float *val) ;

/**********************************************************************/
/******************** message routines - messubs.c ********************/
/**********************************************************************/

/* 'Internal' functions, do not call directly.                               */
void uMessSetErrorOrigin(const char *filename, int line_num) ;
void uMessCrash(char *format, ...) ;

/* External Interface.                                                       */
/* Note that messcrash is a macro and that it makes use of the ',' operator  */
/* in C. This means that the messcrash macro will only produce a single C    */
/* statement and hence can be used within brackets etc. and will not break   */
/* existing code, e.g.                                                       */
/*                     funcblah(messcrash("hello")) ;                        */
/* will become:                                                              */
/* funcblah(uMessSetErrorOrigin(__FILE__, __LINE__), uMessCrash("hello")) ;  */
/*                                                                           */

void messErrorInit (const char *progname) ; /* Record the
					 applications name for use
					 in error messages, etc */
const char *messGetErrorProgram (void) ; /* Returns the
						    application name */

char *messprintf (char *format, ...) ;	  
char *hprintf (AC_HANDLE h, char *format, ...) ;
				/* sprintf into (static!) string */
				/* !!!! beware the code is non reentrant */

void messbeep (void) ; /* make a beep */

void messout (char *format, ...) ;  /* simple message */
void messdump (char *format, ...) ; /* write to log file */
void messerror (const char *format, ...) ; /* error message and write to log file */
void messExit(char *format, ...) ;  /* error message, write to log file & exit */
#define messcrash   uMessSetErrorOrigin(__FILE__, __LINE__), uMessCrash
						  /* abort - but see below */
BOOL messQuery (char *text,...) ;	  /* ask yes/no question */
BOOL messPrompt (char *prompt, char *dfault, char *fmt) ;
	/* ask for data satisfying format get results via freecard() */

const char* messSysErrorText (void) ; 
	/* wrapped system error message for use in messerror/crash() */

int messErrorCount (void);
	/* return numbers of error so far */

BOOL messIsInterruptCalled (void);
	/* return TRUE if an interrupt key has been pressed */

/**** registration of callbacks for messubs ****/

typedef void (*OutRoutine)(const char*) ;
typedef BOOL (*QueryRoutine)(const char*) ;
typedef BOOL (*PromptRoutine)(const char*, char*, char*) ;
typedef BOOL (*IsInterruptRoutine)(void) ;

VoidRoutine	messBeepRegister (VoidRoutine func) ;
OutRoutine	messOutRegister (OutRoutine func) ;
OutRoutine	messDumpRegister (OutRoutine func) ;
OutRoutine	messErrorRegister (OutRoutine func) ;
OutRoutine	messExitRegister (OutRoutine func) ;
OutRoutine	messCrashRegister (OutRoutine func) ;
QueryRoutine	messQueryRegister (QueryRoutine func) ;
PromptRoutine	messPromptRegister (PromptRoutine func) ;
IsInterruptRoutine messIsInterruptRegister (IsInterruptRoutine func) ;

/**** routines to catch crashes if necessary, e.g. when acedb dumping ****/

#include <setjmp.h>

jmp_buf*	messCatchCrash (jmp_buf* ) ;
jmp_buf*	messCatchError (jmp_buf* ) ;
char*	messCaughtMessage (void) ;

  /* if a setjmp() stack context is set using messCatch*() then rather than
     exiting or giving an error message, messCrash() and messError() will
     longjmp() back to the context.
     messCatch*() return the previous value. Use argument = 0 to reset.
     messCaughtMessage() can be called from the jumped-to routine to get
     the error message that would have been printed.
  */

/********************************************************************/
/************** memory management - memsubs.c ***********************/
/********************************************************************/

#define NULL_HANDLE ((AC_HANDLE)0)

AC_HANDLE handleHandleCreate (AC_HANDLE handle) ;
#define handleCreate() handleHandleCreate(0)
#define handleDestroy(handle) messfree(handle)

#if defined(WIN32) && defined(_DEBUG)
#define MEM_DEBUG
#include <crtdbg.h>
#endif

#if !defined(MEM_DEBUG)

void *handleAlloc (void (*final)(void *), AC_HANDLE handle, mysize_t size) ;
    /* handleAlloc is deprecated, use halloc, and blockSetFinalize instead */
char *timeHandleShowNow (AC_HANDLE h) ;
void *halloc(mysize_t size, AC_HANDLE handle) ;
char *strnew(const char *old, AC_HANDLE handle) ;

#else		/* MEM_DEBUG from rbrusk */

void *halloc_dbg(mysize_t size, AC_HANDLE handle, const char *hfname, int hlineno) ;
void *handleAlloc_dbg(void (*final)(void *), AC_HANDLE handle, mysize_t size,
					  const char *hfname, int hlineno) ;
char *strnew_dbg (const char *old, AC_HANDLE handle, const char *hfname, int hlineno) ;

#define halloc(_s, _h) halloc_dbg(_s, _h, __FILE__, __LINE__)
#define handleAlloc(_f, _h, _s) handleAlloc_dbg(_f, _h, _s, __FILE__, __LINE__)
#define strnew(_o, _h) strnew_dbg(_o, _h, __FILE__, __LINE__)
#define messalloc_dbg(_size,_fname,_lineno) halloc_dbg(_size, 0, _fname, _lineno)

#endif

void blockSetFinalise(void *block, void (*final)(void *)) ;
void handleSetFinalise(AC_HANDLE handle, void (*final)(void *), void *arg) ;
void handleInfo (AC_HANDLE handle, int *number, mysize_t *size) ;
#define messalloc(size) halloc(size, 0)
void umessfree (void *cp) ;
#define messfree(cp)  ((cp) ? umessfree((void*)(cp)),(cp)=0,TRUE : FALSE)
void messalloccheck (void) ;	/* can be used anywhere - does nothing
				   unless MALLOC_CHECK set in messubs.c */
int messAllocStatus (int *np) ; /* returns the current number of  allocated buffers
				    *np is current allocated memory in Mb */
int messAllocMaxStatus (int *np) ; /* returns the current number of  allocated buffers
				    *np is maximal allocated memory in Mb */
void messAllocShutdown(void) ;

/********************************************************************/
/******** growable arrays and flexible stacks - arraysub.c **********/
/********************************************************************/

/* to be included after the declarations of AC_HANDLE etc. */
#include "array.h"

/********************************************************************/
/************** file opening/closing from filsubs.c *****************/
/********************************************************************/

void filAddPath (const char *path) ;	/* Adds a set of pathnames to the pathname stack */
void filAddDir (const char *dir) ;	/* Adds a single pathname to the pathname stack */

/* returns an absolute path string for dir in relation to user's CWD */
/* returns pointer to internal static */
char *filGetFullPath (char *dir, AC_HANDLE handle);

/* returns filename part of a pathname. */
/* returns pointer to internal static */
char *filGetFileName (const char *path, AC_HANDLE h) ;

/* returns the file-extension part of a path or file-name */
/* returns pointer to internal static */
char *filGetExtension(const char *path);

char *filName (const char *name, const char *ending, const char *spec) ;
BOOL filCheckName (const char *name, const char *ending, const char *spec) ;
char *filStrictName (const char *name, const char *ending, const char *spec) ;

/* determines time since last modification, FALSE if no file */
BOOL  filAge (const char *name, const char *ending,
			    int *diffYears, int *diffMonths, int *diffDays,
			    int *diffHours, int *diffMins, int *diffSecs);

FILE *filopen (const char *name, const char *ending, const char *spec) ;
FILE *filmail (const char *address) ;
void filclose (FILE* fil) ;

BOOL filremove (const char *name, const char *ending) ;

FILE *filtmpopen (char **nameptr, const char *spec) ;
BOOL filtmpremove (const char *name) ;
void filtmpcleanup (void) ;
BOOL filCopyFile(char *curr_name, char *new_name) ;
BOOL filCreateDir (const char *dirName) ; /* return TRUE if success */

/* file chooser */
typedef FILE* (*QueryOpenRoutine)(char*, char*, char*, const char*, const char*) ;
QueryOpenRoutine filQueryOpenRegister (QueryOpenRoutine newR);
		/* allow graphic file choosers to be registered */

FILE *filqueryopen (char *dirname, char *filname,
				  char *ending, const char *spec, const char *title);

	/* if dirname is given it should be DIR_BUFFER_SIZE long
	     and filname FILE_BUFFER_SIZE long
	   if not given, then default (static) buffers will be used */


/* directory access */
Array filDirectoryCreate (char *dirName,
					char *ending, 
					char *spec);

void filDirectoryDestroy (Array filDirArray);

/*******************************************************************/
/******************** multi threading signals **********************/

void memSetIsMultiThreaded (void) ; /* set, called by wego_init */
BOOL memIsMultiThreaded (void) ;    /* report */

/*******************************************************************/
/************* randsubs.c random number generator ******************/

double randfloat (void) ;
double randgauss (void) ;
int randint (void) ;
void randsave (int *a) ;
void randrestore (int *a) ;


/* Unix debugging.                                                           */
/* put "break invokeDebugger" in your favourite debugger init file */
/* this function is empty, it is defined in messubs.c used in
   messerror, messcrash and when ever you need it.
*/
void invokeDebugger(void) ;




/*******************************************************************/
/************* some WIN32 debugging utilities **********************/

#if defined (WIN32)
#if defined(_DEBUG)
/* See win32util.cpp for these functions */
const char *dbgPos( const char *caller, int lineno, const char *called ) ;
void WinTrace(char *prompt, unsigned long code) ;
void AceASSERT(int condition) ;
void NoMemoryTracking() ;
#else   /* !defined(_DEBUG) */
#define dbgPos(c,l,fil)   (const char *)(fil)
#endif	/* !defined(_DEBUG) */
#else	/* defined(WIN32) */
#define dbgPos(c,l,fil)   (const char *)(fil)
#endif	/* defined(WIN32) */

/*******************************************************************/

/*
* pseudoword features (w1/pseudoword.c)
*/

/*
*
* generate_pseudoword( buffer, number, phoneme_set )
*
* creates a pronouncable string that looks like it might be a real word
*
* buffer is a char pointer to return the result; 2147483647 is 
*	"riyemisukona" in pseudo-japanese and "lyskoplerkar" in
*	pseudo-english, so 15 bytes is enough space on 32 bit machines
*
* number is >= 0, and selects which word number to use.  The returned
*	words are actually just numbers written in a funny number base,
*	but because you don't know the number base, they look like
*	words and are easier to remember.
*
*	This was intended for use mostly with numbers less than 1 000 000,
*	where it produces 3 syllable words in pseudo-English and 4 
*	syllable words in pseudo-Japanese.
*
* phoneme_set is pseudoword_english to make a word that looks like it
*	might be English, or pseudoword_japanese to make a word that
*	looks like it might be Japanese.
*
* NOTE: There are many words that are common to both the English and
*	Japanese word sequences.  If you use words from both, you
*	must handle duplicates yourself.
*
* Bugs:  I created the Japanese pseudo-word sequence, but I do not know
* anything about the sounds of the Japanese language except what I read 
* on a few web sites yesterday.
*
*/

enum pseudoword_phoneme_set { pseudoword_english = 0, pseudoword_japanese = 1 };
int decode_pseudoword (char *cp, enum pseudoword_phoneme_set phoneme_set ) ;
char * generate_pseudoword (char *cp, int n, enum pseudoword_phoneme_set phoneme_set );

/*
* If you don't like single syllable words, do not use any number
* less than this:
*/
extern int 
	pseudoword_english_number_base,
	pseudoword_japanese_number_base;

#endif /* defined(DEF_REGULAR_H) */

/******************************* End of File **********************************/
