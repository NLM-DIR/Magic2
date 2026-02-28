/*  File: aceserver.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *-------------------------------------------------------------------
 * 
 * copied from aceserver.c and modified for the acetcp server
 */

#include "acedb.h"
#include "sysclass.h"
#include "session.h"
#include "pick.h"
#include "bs.h"
#include "systags.h"
#include "tags.h"
#include "call.h"
#include "mytime.h"
#include "aceclient.h"
#include "parse.h"
#include "freeout.h"
#include "version.h"
#include "command.h"
#include "keyset.h"
#include "query.h"
#include "a.h"
#include "dump.h"
#include "lex.h"
#include <signal.h>
#include <netinet/in.h>

#include <sys/types.h>
#include <sys/stat.h>		/* for umask */




/************************************************************/

#ifdef GIFACESERVER
#include "acedbgraph.h"
#endif

#define debug FALSE

static char usage[] = "Usage: aceserver database_dir [-nowrite] [-noswap] [port_number [params]]\n"
		      "  params are clientTimeout:serverTimeout:maxKbytes:autoSaveInterval\n"
		      "  defaults are port             = 12345\n"
		      "               clientTimeout    = 600 seconds\n"
		      "               serverTimeout    = 600 seconds\n"
		      "               maxKbytes        = 0 [no limit]\n"
		      "               autoSaveInterval = 600 seconds\n"
		      "  autoSave checks are only certain at serverTimeout intervals\n"
		      "  -nowrite :  prevents write access\n"
		      "  -noswap  :  calls mlockall() which on some platforms will keep the server in RAM\n"
                      "  -http    :  server will talk a mini http protocol, it will reply to GET but not to acedb clients\n"
		      "  example:  aceserver /local1/worm  20000100  1200:1200:100\n"

		      "            aceserver /local1/worm d20000100  1200:1200:100\n"
		;

/* reason for lack of -option notation is that inetd only allows 5 params in some 
   implementations 
*/

extern int isInteractive ;           /* in freesubs */
extern BOOL parseKeepGoing ;

/* Time outs etc. */
static int clientTimeOut = 600 ;
static int serverTimeOut = 600 ;
static int maxKbytes = 0;
static int autoSaveInterval = 600 ;

static FILE *out = 0 ;

/**************************************************************/
/* this option greatly increases the performances of the server */
/* prevent this process to swap, now and forever */
/* works on solaris and linux, may be not on all platforms */
#include <sys/mman.h>

static void aceCommandNoSwap (void)
{
  int n =  mlockall(MCL_FUTURE) ;
  if (n) messdump (messprintf ("-noswap: mloackall error %d\n", n)) ;
  else messdump (messprintf ("-noswap: mloackall accepted")) ;
}

/******************* INITIALISATION ***************************/
/**************************************************************/

static u_long port = 12345 ;


#if 0
  if (*encore == 3)   /* ace_in */
    { KEYSET ks ;
      if (!look->mayWrite)
	pushText (s,"// Sorry, you do not have write access.\n") ;
      else
	 { 
	   sessionGainWriteAccess() ; /* try to grab it */ 
	   if (!isWriteAccess ())	/* may occur is somebody else grabed it */ 
	     pushText (s,"// Sorry, you do not have write access.\n") ;
	   else
	     { int level = freesettext (question,"") ;
	       freespecial ("\n/\\\"\t") ;
	       ks = keySetCreate () ;
	       parseKeepGoing = TRUE ;
	       parseLevel (level, ks) ;
	       parseKeepGoing = FALSE ;
	       pushText (s, messprintf ("// Read %d objects \n", keySetMax (ks))) ;
	       keySetDestroy (ks) ;
	     }
	 }
      look->lastAccess = timeNow() ;
      goto fin ;
    }
#endif

extern FILE *tcp_log_file  ;

static void openOutFile (void)
{ char buf [12*1024] ;
  strcpy (buf, sessionFilName ("server", "log", 0)) ;
        
  out = filopen (buf, 0, "a") ;
  if (!out)
    out = filopen (buf, 0, "w") ;
  if (!out)
    out = stderr ;
  setbuf (out, NULL) ;
  tcp_log_file = out ;
  /* setbuf (stdout, NULL) ; this flat crashes under daemon on new linux */
}


/* Defines a routine to return the compile date of this file. */
UT_MAKE_GETCOMPILEDATEROUTINE()

/************************************************************/

int main (int argc, const char **argv)
{ 
  u_long p ;
  int n, t, ix ;
  char *cp ;
  const char *ccp ;
/* Declarations for inetd detection LDS 7/2/98 */
  /* struct sockaddr_in saddr; */
  /* int asize = sizeof (saddr); */
  extern void  acetcp_wait_for_client (u_long port) ;

#ifdef GIFACESERVER
  extern void (*gifEntry)(KEYSET,int,BOOL) ;	/* entry point in command.c */
  extern void gifControl (KEYSET ks, int level, BOOL isI)   ;

  gifEntry = gifControl ; 

  freeinit () ;			/* must come before graphInit */

  /* Initialise the graphics package and its interface for acedb.            */
  acedbAppGraphInit(&argc, (char **)argv) ;

  messErrorInit("gifaceserver");   /* Record program name */
#else
  messErrorInit("aceserver");   /* Record program name */
#endif


  out = stderr ; /* this static cannot be initialized at file scope in WIN32? */

  if (getCmdLineOption (&argc, argv, "-nowrite", 0))
    {
       sessionForbidWriteAccess () ;
    }
  if (getCmdLineOption (&argc, argv, "-http", 0))
    {
      extern BOOL TCPWEBSERVER ;
      TCPWEBSERVER = TRUE ;
    }
  if (getCmdLineOption (&argc, argv, "-noswap", 0))
    {
      aceCommandNoSwap () ;
    }
  if (getCmdLineOption (&argc, argv, "-port", &ccp))
    {
      if (sscanf (ccp, "%d", &ix) == 1 && ix >=0)
	port = ix ;
      else 
	{
	  fprintf (out, "%s", usage) ;
	  fprintf (stderr, "%s", usage) ;
	  exit (EXIT_FAILURE) ;
	}
    }
  n = 2 ;
  if (argc > 2 && (sscanf(argv[n],"%lu",&p) == 1))
    { port = p;
      n++;
    }
  if (argc > n)
    {
      ccp = argv[n] ;
      if (sscanf(ccp,"%d",&t) == 1)
	{ 
	  clientTimeOut = t ;
	  if ((cp = strstr(argv[n],":")) &&
	      sscanf(++cp,"%d",&t) == 1)
	    { serverTimeOut = t;
	    if ((cp = strstr(cp,":")) &&
		sscanf(++cp,"%d",&t) == 1)
	      { maxKbytes = t;
	      if ((cp = strstr(cp,":")) &&
		  sscanf(++cp,"%d",&t) == 1)
		autoSaveInterval = t;
	      }
	    }
	  n++;
	}
    }
  if (argc < 2 || argc > n)
    { fprintf (out, "%s", usage) ;
      fprintf (stderr, "%s", usage) ;
      exit (EXIT_FAILURE) ;
    }

  aceInit (argc>1 ? argv[1]: 0) ;  /* Chooses the database directory */

  openOutFile() ;
  cp = getenv ("HOST") ; if (!cp) cp = "(unknown)" ;
  fprintf(out,"\n\n#### TCP-Server starts %s\n", timeShowNow()) ;
  fprintf(out,"#### host=%s  port=%lu  ACEDB=%s\n",  
	  cp, port, sessionFilName(0,0,0)) ;
  fprintf(out,"#### clientTimeout=%d serverTimeout=%d maxKbytes=%d autoSaveInterval=%d\n",
	  clientTimeOut, serverTimeOut, maxKbytes, autoSaveInterval); 

  setsid () ;

  acetcp_wait_for_client(port) ;
  messcrash("Acedb network server error: main loop returned") ;

  return(EXIT_FAILURE) ;			  /* We should not reach here... */
}

/**************************************************/
/**************************************************/
 
 
 
 
