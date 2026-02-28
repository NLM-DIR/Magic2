/*  File: messubs.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
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
 *
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: low level messaging/exit routines, encapsulate vararg messages,
 *              *printf, crash handler, calling of application registered exit
 *              routines for cleaning up etc.
 *
 * Exported functions: see regular.h
 *
 * HISTORY:
 * Last edited: Jul  8 15:11 2003 (edgrif)
 * * May 23 10:18 2000 (edgrif): Fix SANgc07771
 * * Apr 13 16:10 2000 (edgrif): Attempt to improve messages output when
 *              messcrash recurses, totally rubbish before.
 * * Mar 22 14:41 1999 (edgrif): Replaced messErrorInit with messSetMsgInfo.
 * * Jan 27 09:19 1999 (edgrif): Put small orphan functions in utils.c
 * * Jan 13 11:20 1999 (edgrif): Small correction to Jeans changes to messGetErrorProgram
 *              and messGetErrorFile.
 * * Nov 19 13:26 1998 (edgrif): Removed the test for errorCount and messQuery
 *              in messerror, really the wrong place.
 * * Oct 22 15:26 1998 (edgrif): Replaced strdup's with strnew.
 * * Oct 21 15:07 1998 (edgrif): Removed messErrorCount stuff from graphcon.c
 *              and added to messerror (still not perfect), this was a new.
 *              bug in the message system.
 * * Sep 24 16:47 1998 (edgrif): Remove references to ACEDB in messages,
 *              change messExit prefix to "EXIT: "
 * * Sep 22 14:35 1998 (edgrif): Correct errors in buffer usage by message
 *              outputting routines and message formatting routines.
 * * Sep 11 09:22 1998 (edgrif): Add messExit routine.
 * * Sep  9 16:52 1998 (edgrif): Add a messErrorInit function to allow an
 *              application to register its name for use in crash messages.
 * * Sep  3 11:32 1998 (edgrif): Rationalise strings used as prefixes for
 *              messages. Add support for new messcrash macro to replace
 *              messcrash routine, this includes file/line info. for
 *              debugging (see regular.h for macro def.) and a new
 *              uMessCrash routine.
 * * Aug 25 14:51 1998 (edgrif): Made BUFSIZE enum (shows up in debugger).
 *              Rationalise the use of va_xx calls into a single macro/
 *              function and improve error checking on vsprintf.
 *              messdump was writing into messbuf half way up, I've stopped
 *              this and made two buffers of half the original size, one for
 *              messages and one for messdump.
 * * Aug 21 13:43 1998 (rd): major changes to make clean from NON_GRAPHICS
 *              and ACEDB.  Callbacks can be registered for essentially
 *              all functions.  mess*() versions continue to centralise
 *              handling of ... via stdarg.
 * * Aug 20 17:10 1998 (rd): moved memory handling to memsubs.c
 * * Jul  9 11:54 1998 (edgrif): 
 *              Fixed problem with SunOS not having strerror function, system
 *              is too old to have standard C libraries, have reverted to
 *              referencing sys_errlist for SunOS only.
 *              Also fixed problem with getpwuid in getLogin function, code
 *              did not check return value from getpwuid function.
 * * Jul  7 10:36 1998 (edgrif):
 *      -       Replaced reference to sys_errlist with strerror function.
 * * DON'T KNOW WHO MADE THESE CHANGES...NO RECORD IN HEADER....(edgrif)
 *      -       newformat added for the log file on mess dump.
 *      -       Time, host and pid are now always the first things written.
 *      -       This is for easier checking og the log.wrm with scripts etc.
 *      -       Messquery added for > 50 minor errors to ask if user wants to crash.
 *      -       Made user,pid and host static in messdump.
 * * Dec  3 15:52 1997 (rd)
 * 	-	messout(): defined(_WINDOW) =>!defined(NON_GRAPHIC)
 * * Dec 16 17:26 1996 (srk)
 * * Aug 15 13:29 1996 (srk)
 *	-	WIN32 and MACINTOSH: seteuid() etc. are stub functions
 * * Jun 6 10:50 1996 (rbrusk): compile error fixes
 * * Jun  4 23:31 1996 (rd)
 * Created: Mon Jun 29 14:15:56 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: messubs.c,v 1.22 2017/02/27 14:57:39 mieg Exp $ */

#include <assert.h>
#include <errno.h>
#include "regular.h"
#include "mytime.h"


/* This is horrible...a hack for sunos which is not standard C compliant.    */
/* to allow accessing system library error messages, will disappear....      */
#ifdef SUN
extern const char *sys_errlist[] ;
#endif


/* Mac has its own routine for crashing, see messcrash for usage.            */
#if !defined(MACINTOSH) 
extern void crashOut (char* text) ;
#endif



/* This buffer is used only by the routines that OUTPUT a message. Routines  */
/* that format messages into buffers (e.g. messprintf, messSysErrorText)     */
/* have their own buffers. Note that there is a problem here in that this    */
/* buffer can be overflowed, unfortunately because we use vsprintf to do     */
/* our formatting, this can only be detected after the event.                */
/*                                                                           */
/* Constraints on message buffer size - applicable to ALL routines that      */
/* format externally supplied strings.                                       */
/*                                                                           */
/* BUFSIZE:  size of message buffers (messbuf, a global buffer for general   */
/*           message stuff and a private ones in messdump & messprintf).     */
/* PREFIX:   length of message prefix (used to report details such as the    */
/*           file/line info. for where the error occurred.                   */
/* MAINTEXT: space left in buffer is the rest after the prefix and string    */
/*           terminator (NULL) are subtracted.                               */
/* Is there an argument for putting this buffer size in regular.h ??         */
/*                                                                           */
enum {BUFSIZE = 32768, PREFIXSIZE = 1024, MAINTEXTSIZE = BUFSIZE - PREFIXSIZE - 1} ;

/* Macro to format strings using va_xx calls 
 *  since it contains declarations, you must enclose it 
 * at the begining of a {} block and use the resulting
 * char *message
 * in the same block
 * each of the function using this macro define their own
 * static elastic buffer
 *   - there is no hard limit to the message size
 *   - BUT the buffer is static so the function is non reentrant
 *     and NOT thread safe
 */

#define ACFORMAT(_prefix)  static Stack s = 0 ; \
  int len, prefixLength ; \
  char *message ; \
  va_list ap; \
  s = stackReCreate (s, 500) ; \
\
  va_start(ap, format);\
  len = utPrintfSizeOfArgList (format, ap) ;\
  va_end (ap) ;\
\
  prefixLength = strlen(_prefix) ; \
  stackExtend (s, len + prefixLength + 1) ;\
  catText (s, _prefix) ;\
  message = stackText (s, 0) ;\
\
  va_start(ap, format);\
  len = vsprintf(message + prefixLength, format, ap) ;\
  va_end (ap) ;

/* Some standard defines for titles/text for messages:                       */
/*                                                                           */
#define ERROR_PREFIX "ERROR: "
#define EXIT_PREFIX "EXIT: "
#define CRASH_PREFIX_FORMAT "FATAL ERROR %s reported by %s at line %d: "
#define FULL_CRASH_PREFIX_FORMAT "FATAL ERROR %s reported by program %s, in file %s, at line %d: "
#if defined(MACINTOSH)
#define      SYSERR_FORMAT             "system error %d"
#else
#define      SYSERR_FORMAT             "system error %d - %s"
#endif
#define PROGNAME "The program"

/* messcrash now reports the file/line no. where the messcrash was issued    */
/* as an aid to debugging. We do this using a static structure which holds   */
/* the information and a macro version of messcrash (see regular.h), the     */
/* structure elements are retrieved using access functions.                  */
typedef struct _MessErrorInfo
  {
  char *progname ;				  /* Name of executable reporting error. */
  char *filename ;				  /* Filename where error reported */
  int line_num ;				  /* Line number of file where error
						     reported. */
  } MessErrorInfo ;

static MessErrorInfo messageG = {NULL, NULL, 0} ;

static int messGetErrorLine() ;
static const char *messGetErrorFile() ;


/* Keeps a running total of errors so far (incremented whenever messerror is */
/* called).                                                                  */
static int errorCount_G = 0 ;


/* Function pointers for application supplied routines that are called when  */
/* ever messerror or messcrash are called, enables application to take       */
/* action on all such errors.                                                */
static jmp_buf *errorJmpBuf = 0 ;
static jmp_buf *crashJmpBuf = 0 ;



/***************************************************************/
/********* call backs and functions to register them ***********/

static VoidRoutine	  beepRoutine = 0 ;
static OutRoutine	  outRoutine = 0 ;
static OutRoutine	  dumpRoutine = 0 ;
static OutRoutine	  errorRoutine = 0 ;
static OutRoutine	  exitRoutine = 0 ;
static OutRoutine	  crashRoutine = 0 ;
static IsInterruptRoutine isInterruptRoutine = 0 ;

VoidRoutine messBeepRegister (VoidRoutine func)
{ VoidRoutine old = beepRoutine ; beepRoutine = func ; return old ; }

OutRoutine messOutRegister (OutRoutine func)
{ OutRoutine old = outRoutine ; outRoutine = func ; return old ; }

OutRoutine messDumpRegister (OutRoutine func)
{ OutRoutine old = dumpRoutine ; dumpRoutine = func ; return old ; }

OutRoutine messErrorRegister (OutRoutine func)
{ OutRoutine old = errorRoutine ; errorRoutine = func ; return old ; }

OutRoutine messExitRegister (OutRoutine func)
{ OutRoutine old = exitRoutine ; exitRoutine = func ; return old ; }

OutRoutine messCrashRegister (OutRoutine func)
{ OutRoutine old = crashRoutine ; crashRoutine = func ; return old ; }


IsInterruptRoutine messIsInterruptRegister (IsInterruptRoutine func)
{ IsInterruptRoutine old = isInterruptRoutine ; isInterruptRoutine = func ; return old ; }



/***************************************************/
BOOL messIsInterruptCalled (void)
{
  if (isInterruptRoutine)
    return (*isInterruptRoutine)() ;

  /* unless a routine is registered, we assume no interrupt
     (e.g. F4 keypress in graph-window) has been called */
  return FALSE;
}


/* The message output routines.                                              */
/*                                                                           */
/*                                                                           */


/***************************************************/
void messbeep (void)
{
  if (beepRoutine)
    (*beepRoutine)() ;
  else
    { printf ("%c",0x07) ;  /* bell character, I hope */
      fflush (stdout) ;	/* added by fw 02.Feb 1994 */
    }
}


/*******************************/

void messout (char *format,...)
{
  if (outRoutine)
   {
     ACFORMAT("") ; /* some compilers do not like a 0 as parameter here */
     (*outRoutine)(message) ;
   }
  else
    {
      va_list ap; 
      va_start(ap, format);
      fprintf (stdout, format, ap) ;
      va_end (ap) ;
    }
}

/*****************************/


void  messPrompt2 (void)
{ 
  if (1) filclose (0) ;

 return ;
}

/*****************************/

/*****************************************************************/

void messdump (char *format,...)
{
  if (dumpRoutine)
    {
     ACFORMAT("") ; /* some compilers do not like a 0 as parameter here */
     (*dumpRoutine)(message) ;
    }
}


/*****************************************/


/* Access function for returning running error total.                        */
int messErrorCount (void) { return errorCount_G ; }


/* Output a non-fatal error message, for all messages a call to messdump is  */
/* made which may result in the message being logged. The single error count */
/* is also incremented so that functions can use this to check how many      */
/* errors have been recorded so far.                                         */
void messerror (const char *format, ...)
{
  char *prefix = ERROR_PREFIX ;
  ACFORMAT(prefix) ;

  /* always increment the error count.                                       */
  ++errorCount_G ;

  /* If application registered an error handler routine, call it.            */
  if (errorJmpBuf)
    longjmp (*errorJmpBuf, 1) ;

  /* Log the message.                                                        */
  messdump("%s", message) ;

  /* Now report the error to the user.                                       */
  if (errorRoutine)
    {
      if (errorCount_G < 30)
	(*errorRoutine)(message) ;
    }
  else
    fprintf (stderr, "%s\n", message) ;

  invokeDebugger () ;
}

/*******************************/

/* Use this function for errors that while being unrecoverable are not a     */
/* problem with the acedb code, e.g. if the user starts xace without         */
/* specifying a database.                                                    */
/* Note that there errors are logged but that this routine will exit without */
/* any chance to interrupt it (e.g. the crash routine in uMessCrash), this   */
/* could be changed to allow the application to register an exit handler.    */
/*                                                                           */
void messExit(char *format, ...)
{
  ACFORMAT (EXIT_PREFIX) ;

  if (exitRoutine)
    (*exitRoutine)(message) ;
  else
    fprintf (stderr, "%s\n", message) ;

  messdump ("%s", message) ;
  
  exit(EXIT_FAILURE) ;
  
  return ;					  /* Should never get here. */
}


/*******************************/

/* This is the routine called by the messcrash macro (see regular.h) which   */
/* actually does the message/handling and exit.                              */
/* This routine may encounter errors itself, in which case it will attempt   */
/* to call itself to report the error. To avoid infinite recursion we limit  */
/* this to just one reporting of an internal error and then we abort.        */
/*                                                                           */
void uMessCrash (char *format, ...)
{
  enum {MAXERRORS = 1} ;
  static int internalErrors = 0 ;
  static char prefix[1024] ;
  int rc ;
  extern void abort(void) ;

  /* Check for recursive calls and abort if necessary.                       */
  if (internalErrors > MAXERRORS) 
    {
      fprintf (stderr, "%s : fatal internal error, abort", 
	       messageG.progname);
      abort() ;
    }
  
  /* Construct the message prefix, adding the program name if possible.      */
  if (messGetErrorProgram() == NULL)
    rc = sprintf(prefix, CRASH_PREFIX_FORMAT, timeShowNow (), messGetErrorFile(), messGetErrorLine()) ;
  else
    rc = sprintf(prefix, FULL_CRASH_PREFIX_FORMAT,
		 messGetErrorProgram(), timeShowNow (), messGetErrorFile(), messGetErrorLine()) ;
  if (rc < 0) 
    {
      fprintf (stderr, "sprintf failed inside w1/messubs.c:uMessCrash") ;
      abort() ;
    }
    
  /* Format the message string.                                              */
  {
    ACFORMAT (prefix) ;
    
    
    if (crashJmpBuf)		/* throw back up to the function that registered it */
      longjmp(*crashJmpBuf, 1) ;
    internalErrors++ ;          /* must come after the longjmp
				 * if user registered a routine he may accepts loads of errors
				 * for example in dump.c
				 */
    messdump("%s", message) ;
    
    if (crashRoutine)
      (*crashRoutine)(message) ;
    else
      fprintf(stderr, "%s\n", message) ;
    
    invokeDebugger() ;
    
    exit(EXIT_FAILURE) ;
  }
  
  return ;					  /* Should never get here. */
} /* uMessCrash */


/******* interface to crash/error trapping *******/

jmp_buf* messCatchError (jmp_buf* neuf)
{
  jmp_buf* old = errorJmpBuf ;
  errorJmpBuf = neuf ;
  return old ;
}

jmp_buf* messCatchCrash (jmp_buf* neuf)
{
  jmp_buf* old = crashJmpBuf ;
  crashJmpBuf = neuf ;
  return old ;
}

static char messbuf[BUFSIZE] ;
char* messCaughtMessage (void) { return messbuf ; }

/* Message formatting routines.                                              */
/*                                                                           */
/*                                                                           */

/* This function writes into its own buffer,note that subsequent calls will  */
/* overwrite this buffer.                                                              */
/*                                                                           */


char *messprintf (char *format, ...)
{
  ACFORMAT ("") ;
  return message ;
}

/*****************************************************************************/

/* Used internally for formatting into a specified buffer.                   */
/* (currently only used as a cover function to enable us to use ACEFORMAT-   */
/* STRING from messSysErrorText)                                             */
static char *printToBuf(char *buffer, unsigned int buflen, char *format, ...)
{
  ACFORMAT ("") ;
  
  if (message && strlen (message) > buflen)
    {
      fprintf (stderr, "printToBuf() : "
	       "the message is longer than the provided buffer length (%d) : %s"
	       , buflen, message) ;
      
      invokeDebugger();
      exit (EXIT_FAILURE);
    }
  strncpy (buffer, message, buflen) ;
  
  return message ;
}

/* Return the string for a given errno from the standard C library.          */
/*                                                                           */
const char* messSysErrorText (void)
  {
  enum {ERRBUFSIZE = 2000} ;				    /* Should be enough. */
  static char errmess[ERRBUFSIZE] ;
  char *mess ;

#ifdef SUN
  /* horrible hack for Sunos/Macs(?) which are not standard C compliant */
  mess = printToBuf(&errmess[0], ERRBUFSIZE, SYSERR_FORMAT, errno, sys_errlist[errno]) ;
#elif defined(MACINTOSH)
  mess = printToBuf(&errmess[0], ERRBUFSIZE, SYSERR_FORMAT, errno) ;
#else
  mess = printToBuf(&errmess[0], ERRBUFSIZE, SYSERR_FORMAT, errno, strerror(errno)) ;
#endif

  return mess ? mess : "(null)" ; /* mieg may be 0 */
  }


char* messErrnoText (int error)
  {
  enum {ERRBUFSIZE = 2000} ;				    /* Should be enough. */
  static char errmess[ERRBUFSIZE] ;
  char *mess ;

#ifdef SUN
  /* horrible hack for Sunos/Macs(?) which are not standard C compliant */
  mess = printToBuf(&errmess[0], ERRBUFSIZE, SYSERR_FORMAT, error, sys_errlist[errno]) ;
#elif defined(MACINTOSH)
  mess = printToBuf(&errmess[0], ERRBUFSIZE, SYSERR_FORMAT, error) ;
#else
  mess = printToBuf(&errmess[0], ERRBUFSIZE, SYSERR_FORMAT, error, strerror(errno)) ;
#endif

  return mess ;
  }
void messSetMsgInfo(char *progname, char *progversion,
		    char *userid, char *machineid)
{
  return ;
}


/********************** crash file/line info routines ************************/
/* When the acedb needs to crash because there has been an unrecoverable     */
/* error we want to output the file and line number of the code that         */
/* detected the error. Here are the functions to do it.                      */
/*                                                                           */

/* Applications can optionally initialise the error handling section of the  */
/* message package, currently the program name can be set (argv[0] in the    */
/* main routine) as there is no easy way to get at this at run time except   */
/* from the main.                                                            */
/*                                                                           */
void messErrorInit(const char *progname)
  {

  if (progname != NULL) messageG.progname = strnew(progname, 0) ;

  return ;
  }

/* This function is called by the messcrash macro which inserts the file and */
/* line information using the __FILE__ & __LINE__ macros.                    */
/*                                                                           */
void uMessSetErrorOrigin(const char *filename, int line_num)
{
  
  assert(filename != NULL && line_num != 0) ;
  
  /* We take the basename here because __FILE__ can be a path rather than    */
  /* just a filename, depending on how a module was compiled.                */
  messageG.filename = strnew(filename, 0) ;
  
  messageG.line_num = line_num ;
}

/* mieg: protected these func against bad return, was crashing solaris server */
/* Access functions for message error data.                                  */
const char *messGetErrorProgram(void)
{
  return messageG.progname ?  messageG.progname : "programme_name_unknown"  ;
}  

static const char *messGetErrorFile(void)
{
  return messageG.filename ? messageG.filename  : "file_name_unknown" ; 
}  

static int messGetErrorLine(void)
{
  return messageG.line_num ;
}  

/*************************************************************************/
/**** end of file ****/

