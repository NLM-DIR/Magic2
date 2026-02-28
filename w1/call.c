/*  File: call.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: provides hooks to optional code, basically a call by
 			name dispatcher
	        plus wscripts/ interface
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 19 14:40 1998 (fw)
 * * Mar  3 15:51 1996 (rd)
 * Created: Mon Oct  3 14:05:37 1994 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: call.c,v 1.13 2017/06/19 20:54:32 mieg Exp $ */

/* #include "acedb.h" */
#include "call.h"
#include "mytime.h"
#include <ctype.h>		/* for isprint */
#include <sys/wait.h>

/************************* call by name package *********************/

typedef struct
{ const char *name ;
  CallFunc func ;
} CALL ;

static Array calls ;	/* array of CALL to store registered routines */

/************************************/

static int callOrder (const void *a, const void *b) 
{ return strcmp (((CALL*)a)->name, ((CALL*)b)->name) ; }

void callRegister (const char *name, CallFunc func)
{
  CALL c ;

  if (!calls)
    calls = arrayCreate (16, CALL) ;
  c.name = name ; c.func = func ;
  if (!arrayInsert (calls, &c, callOrder))
    messcrash ("Duplicate callRegister with name %s", name) ;
}

BOOL callExists (const char *name)
{
  CALL c ;
  int i ;

  c.name = name ;
  return (calls && arrayFind (calls, &c, &i, callOrder)) ? TRUE : FALSE  ;
}


#include <stdarg.h>   /* for va_start */

#if ! defined(GPLUSPLUS)
/*  i cannot compile that for C++ */
BOOL call (char *name, ...)
{
  va_list args ;
  CALL c ;
  int i ;

  c.name = name ;
  if (calls && arrayFind (calls, &c, &i, callOrder))
    { va_start(args, name) ;
      (*(arr(calls,i,CALL).func))(args) ;
      va_end(args) ;
      return TRUE ;
    }

  return FALSE ;
}
#endif

/***************** routines to run external programs *******************/

/* ALL calls to system() and popen() should be through these routines
** First, this makes it easier for the Macintosh to handle them.
** Second, by using wscripts as an intermediate one can remove system
**   dependency in the name, and even output style, of commands.
** Third, if not running in ACEDB it does not look for wscripts...
*/

static char *buildCommand (const char *dir, const char *script, const char *args)
{
  static Stack command = 0 ;
#ifdef ACEDB  /* until we resolve this bit, we have to include wscripts bit even ifndef ACEDB */
  char *cp ;
  static Stack s = 0 ;		/* don't use messprintf() - often used to make args */
  s = stackReCreate (s, 32) ; 
  if (!dir)
    {
      catText (s, "wscripts/") ; 
      catText (s, script) ;
      if ((cp = filName (stackText (s, 0), 0, "x")))
	script = cp ;  /* mieg else fall back on direct unix call */
    }
#endif

  command = stackReCreate (command, 128) ;
  if (dir)
    { catText (command, "cd ") ;
      catText (command, dir) ;
      catText (command, "; ") ;
    }
  catText (command, script) ;
  if (args)
    { catText (command, " ") ;
      catText (command, args) ;
    }
  return stackText (command, 0) ;
}

int callCdScript (const char *dir, const char *script, const char *args)
{
#if !defined(MACINTOSH)
  return callSystem (buildCommand (dir, script, args)) ;
#else
  return -1 ;
#endif
}

int callScript (const char *script, const char *args)
{
  return callCdScript (0, script, args) ;
}

FILE* callCdScriptPipe (const char *dir, const char *script, const char *args)
{
  char *command = buildCommand (dir, script, args) ;
  FILE *pipe ;
  int peek ;

#if !(defined(MACINTOSH) || defined(WIN32))
  pipe = popen (command, "r" ) ;
#elif defined(WIN32)
  pipe =  _popen (command, "rt")  ;
#else	/* defined(MACINTOSH) */
  return 0 ;
#endif	/* defined(MACINTOSH) */

  peek = fgetc (pipe) ;		/* first char from popen on DEC
				   seems often to be -1 == EOF!!!
				*/
  if (isprint(peek)) ungetc (peek, pipe) ;
#ifdef DEBUG
  printf ("First char on callCdScriptPipe is %c (0x%x)\n", peek, peek) ;
#endif
  return pipe ;
}

FILE* callScriptPipe (const char *script, const char *args)
{
  return callCdScriptPipe (0, script, args) ;
}

/* mieg 2007_01_17
 *
 * Man page on solaris says:
 * The following function can be used in a  multithreaded  pro-
 * cess  in  place  of  the  most  common  usage  of the Unsafe
 *system(3C) function:
 */ 
   
int callSystem (const char *cmd)
{
  FILE *p;
  
  if ((p = popen(cmd, "w")) == NULL)
    return (-1);
  return (pclose(p));
}

/********************************************************************/
/********************************************************************/
/* from Vahan Simonyan   2007_01_19 
 *
 * Problem: in SUSE linux pipes are unidirectional
 * although they are bi-directional in Solaris and other linuxes
 * the following function creates a bi-directional pipe
 *
 * Usage: 
 

int main(int argc, char * argv[], char * envp[])
{

        const char * inputBuf = ""
                "This text is to be directed \n"
                "to grep through pipe. \n"
                "Output of the pipe is supposed to \n"
                "contain lines which have \n"
                "the word Danielle inside. \n"
                "Jean won't be found here\n";
// ifndef ACEDB
        char outBuf[1024]; sprintf(outBuf,"not executed");
        popen2("grep Danielle", inputBuf, strlen(inputBuf)
                              , outBuf,sizeof(outBuf));
        printf("###############o\n%s\n################\n",outBuf);
// else
        ACEOUT fo = ... ; // rename to callPipe and use and acedb elastic outFlow 
	callPipe ("grep Danielle", inputBuf, strlen(inputBuf), fo) ;
#endif
}

*/
/********************************************************************/
#include <sys/types.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "aceio.h"

#define READ 0
#define WRITE 1

#ifdef OUTSIDE_ACEDB
int popen2 (const char *command, const char *inBuf, unsigned int inSize
	      char *outBuf, unsigned int outSize)
#else
int callPipe (const char *command, const char *inBuf, unsigned int inSize
		, ACEOUT fo)
#endif
{
  int parent2child[2] ;
  int child2parent[2] ;
  pid_t pid ;
  int status, inFl, outFl ;
  unsigned int nn = 1, nByteRead = 0, nByteWritten = 0 ;
  
  /* create pipes */
  if (pipe (parent2child) != 0 || pipe (child2parent) != 0)
      return -1000 ; /* error */
  
  /* fork your process */
  pid = fork () ;
  if (pid < 0) 
    return(int) pid ; /* error */
  
  /* successfull forking */
  else if (pid == 0) /* I am the child, fix the streams and execute */
    {      
      close (parent2child[WRITE]) ;
      dup2 (parent2child[READ], READ);    /* make fildes [READ] my effective stdin */
      close (child2parent[READ]) ;
      dup2 (child2parent[WRITE], WRITE);  /* make fildes[WRITE] my effective stdout */
      
      execl ("/bin/sh", "sh", "-c", command, NULL); /* execute the command */
      perror ("execl"); /* if there are any errors print them out */

      exit (0); /*  done */
    }
  
  /* here we are in the parent process */
  close (parent2child[READ]) ;
  close (child2parent[WRITE]) ;
  outFl = parent2child[WRITE] ; /* child reads where I write */
  inFl = child2parent[READ]   ; /* child writes where I read */
  
  /* transmit data to the child */
  if (inSize)
    nByteWritten = write (outFl, inBuf, inSize) ;
  close (outFl) ; /* we close here so the child doesn't wait for more input */
 
  /* read the results */

  if (inSize ==  nByteWritten)
    {
#ifdef OUTSIDE_ACEDB
      *outBuf = 0 ;
      nn = read (inFl, outBuf, outSize) ;
      outBuf[nn] = 0 ;
      nByteRead += nn ;
#else
      {
	char outBuf[100000] ;
	int outSize = 100000 ;
	while (nn)
	  {
	    nn = read (inFl, outBuf, outSize) ;
	    aceOutBinary (fo, outBuf, nn) ; 
	    nByteRead += nn ;
	  }
      }
#endif
    }
  /* final clean up */
  close (parent2child[WRITE]) ;
  close (child2parent[READ]) ;
  close (inFl) ;
  status = 0 ;
  waitpid (pid, &status, 0) ;
  return nByteRead ;
}

static const char *callPipeMagic = "callPipeMagic" ;

typedef struct callPipeStruct { const char *magic, *command ;; pid_t child ; int inFlow, outFlow ;} CALL_PIPE ;
BOOL callPipeDestroy (CALL_PIPE *pp)
{
  if (pp->magic == callPipeMagic)
    {
      pp->magic = 0 ;
      close (pp->inFlow) ;
      close (pp->outFlow) ;
      messfree (pp->command) ;
      /* kill child */
      return TRUE ;
    } 
  return FALSE ;
}

CALL_PIPE *callPipeCreate (const char *command)
{
  int parent2child[2] ;
  int child2parent[2] ;
  pid_t pid ;
  CALL_PIPE *pp = (CALL_PIPE *) messalloc (sizeof(CALL_PIPE)) ;

  pp->command = strnew (command, 0) ;
  pp->magic = callPipeMagic ;

  /* create pipes */
  if (pipe (parent2child) != 0 || pipe (child2parent) != 0)
      return 0 ;  /* error */
  
  /* fork your process */
  pid = fork () ;
  if (pid < 0) 
    return 0 ; /* error */
  
  /* successfull forking */
  else if (pid == 0) /* I am the child, fix the streams and execute */
    {      
      close (parent2child[WRITE]) ;
      dup2 (parent2child[READ], READ);    /* make fildes [READ] my effective stdin */
      close (child2parent[READ]) ;
      dup2 (child2parent[WRITE], WRITE);  /* make fildes[WRITE] my effective stdout */
      
      execl ("/bin/sh", "sh", "-c", command, NULL); /* execute the command */
      perror ("execl"); /* if there are any errors print them out */
      
      /* DO NOT exit */
    }

  /* here we are in the parent process */
  close (parent2child[READ]) ;
  close (child2parent[WRITE]) ;
  pp->outFlow = parent2child[WRITE] ; /* child reads where I write */
  pp->inFlow = child2parent[READ]   ; /* child writes where I read */
  return pp ;
}

int callPipeRun (CALL_PIPE *pp, const char *inBuf, unsigned int inSize, ACEOUT fo)
{
  CALL_PIPE *pp1 ;
  char outBuf[100000] ;
  int pass = 0, nIn = 1, nOut = 0, nByteRead = 0, outSize = 100000 ;

  if (!pp || pp->magic == callPipeMagic)
    messcrash ("Invalid handle passed to  callPipeRun") ;

  /* transmit data to the child */
 
  while (pass++ == 0 && (nOut = write (pp->outFlow, inBuf, inSize) < 0))
    {
      /* write failed, try once to reconstruct the pipe */
      callPipeDestroy (pp) ;
      close (pp->inFlow) ;
      close (pp->outFlow) ;
      pp1 = callPipeCreate (pp->command) ;
      messfree (pp->command) ;
      *pp = *pp1 ;
    }
  
  if (nOut >= 0)
      /* read the results */
    while (nIn)
      {
	nIn = read (pp->inFlow, outBuf, outSize) ;
	aceOutBinary (fo, outBuf, nIn) ; 
	nByteRead += nIn ;
      }
  return nByteRead ;
}

/********************************************************************/
/********************************************************************/
/********************************************************************/
