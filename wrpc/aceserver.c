/*  File: aceserver.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:      Network Query Server   
 *     Check authorizations and times-out

 * added back , because we need to know to correctly connect (mieg)
 * [rd 970722] remove "rpc.*" argv[0] test for isDaem, because that was
 * hardwired in anyway.

    checkLife keeps track of the number of active clients.  
   If a client is not active for clientTimeOut seconds, it will be disconnected.
   If there are no more active clients, the server itself will disconnect
   after serverTimeOut seconds.
 
   If a time_out is set to 0 (zero) expiration is deactivated

   Note that, when saving, you save all clients at once, there
   is NO WAY to save independantly a single client session.
 *
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  3 18:09 1998 (fw)
 * * Dec  2 16:32 1998 (edgrif): Correct decl. of main, add code to
 *              record build time of this module.
 * * Oct 15 10:57 1998 (edgrif): Add new graph/acedb initialisation call.
 * * Feb 25 23:00 1996 (mieg)
 * * Mar 97 (mieg) SIGPIPE handler
 * Created: Wed Nov 25 12:26:06 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: aceserver.c,v 1.19 2020/05/30 16:50:36 mieg Exp $ */

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
static void copy_wspec (BOOL encode) ;
/************************************************************/

#ifdef GIFACESERVER
#include "acedbgraph.h"
#endif

#define debug FALSE

static char usage[] = "Usage: aceserver database_dir [port_number [params]]\n"
		      "  params are clientTimeout:serverTimeout:maxKbytes:autoSaveInterval\n"
		      "  defaults are port             = 20000101\n"
		      "               clientTimeout    = 600 seconds\n"
		      "               serverTimeout    = 600 seconds\n"
		      "               maxKbytes        = 0 [no limit]\n"
		      "               autoSaveInterval = 600 seconds\n"
		      "  autoSave checks are only certain at serverTimeout intervals\n"
		      "  example:  aceserver /local1/worm  20000100  1200:1200:100\n"

		      "            aceserver /local1/worm d20000100  1200:1200:100\n"
		      "  note that 'port' is really an RPC program number\n" ;

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
static int nActiveClients = 0; 


static BOOL isDaemon = FALSE ;
static BOOL isSGI = FALSE ;
static Array wspecTar = 0 , wspecTarEncoded = 0 ;

#if defined(WIN32)
#define DCESERVER
extern void my_Signal( void (*callback)(int arg) ) ;
extern void my_Alarm(unsigned int seconds) ;
#define signal(S,C) my_Signal(C) /* ignore the SIGALRM type */
#define alarm(T) my_Alarm(T)

typedef unsigned long u_long ;
extern void stopDCEServer(void) ;
#endif


static FILE *out = 0 ;
static BOOL refuseNewClient = FALSE ;
static Associator clientAss = 0 ;

/************************************************************/

typedef struct LOOKSTUFF { int writeMagic, readMagic ; 
			   int clientId ;
			   char *magicReadFileName ;
			   char *magicWriteFileName ;
			   mytime_t lastAccess;
			   void *aceCommand ;
			   KEY lastCommand ;
			   int public ;
                           BOOL mayWrite ;
			 } *LOOK ;

/******************* INITIALISATION ***************************/
/**************************************************************/

static void clientDestroy (LOOK look)
{ char *cp, *cp0 = 0 ;

  fprintf (out, "\n%s Closing Client %d, %d active clients\n", 
	       timeShowNow(), look->clientId, --nActiveClients) ;
  cp = cp0 + look->clientId ;
  assRemove (clientAss, cp) ;
  aceCommandDestroy (look->aceCommand) ;
  messfree (look->magicReadFileName) ;
  messfree (look->magicWriteFileName) ;
  messfree (look) ;
}

static void wormClose(void)
{ extern void closePortMap(void);
  
#if defined(DCESERVER)
  stopDCEServer() ;
#else
  if (!isDaemon)  /* in daemon mode, do not close port mapper 
		  *  this is needed on alpha, may be machine dependant ?
		  */
    {
      fprintf(out,"\n\n%s #### Closing port map\n",timeShowNow()) ;
      closePortMap();	
    }
#endif

  if (isWriteAccess())
    aceQuit (TRUE);
  else
    aceQuit (FALSE);

  fprintf(out,"\n\n%s #### Server normal exit\n  A bientot \n\n",
	  timeShowNow()) ;
  if (out && out!= stderr)
    { filclose(out) ; out = 0 ; }

  exit(0) ;
} /* wormClose */

extern void wait_for_client (u_long port, BOOL isDaemon) ;
static u_long port = DEFAULT_PORT ;


/**************************************************************/
#if !defined(DCESERVER)

  /* SIGPIPE, generated if say the client types CTRL-C 
   * i hope resuming operations in this way is ok, not sure
   * mieg: march 97
   */
static void sigPipeError (int dummy)  /* dummy:  to please signal() prototypes */
{ 
  fprintf(out, "// SIGPIPE ERROR, hopefully on client side, ignore // %s\n", timeShowNow()) ;
  wait_for_client(port, isSGI && isDaemon ) ;
}
#endif
/**************************************************************/
  /* checkLife is called by alarm, 
   * but also after every query including from unauthorised clients
   */
static void checkLife (int dummy)  /* dummy:  to please signal() prototypes */
{ static BOOL dying = FALSE ;
  LOOK look;
  int dt ;
  char *cp = 0 ;
  mytime_t t = timeNow () ;
  static   mytime_t last_t = 0 ;

  if (refuseNewClient && !nActiveClients) 
    { fprintf(out,"// %s server shut down on request\n\n", timeShowNow()) ;
      wormClose () ;
    }

  if (dying && !nActiveClients)
    { if (timeDiffSecs (last_t, t, &dt) &&
	10*dt >= 9*serverTimeOut)  /* there seem to exist rounding errors, hence the necessity of 10/9 */
	{ fprintf(out,"// %s Reached timeout, server will shut down\n\n", timeShowNow()) ;
	  wormClose () ;
	}
      return ; /* do not update last_t in this case */
    }

  while (assNext (clientAss, &cp, &look)) /* need to loop over all clients */
    if (timeDiffSecs (look->lastAccess, t, &dt) &&
	dt > clientTimeOut) 
      { clientDestroy (look) ;
      fprintf(out, "// %s Destroying an idle client\n", timeShowNow()) ;
      cp = 0 ;    /* restart the search - needed to keep assNext() happy after assRemove */
      } 

      /* else printf("Not Destroying an idle client\n") ; */
  
  if (autoSaveInterval > 0)
    sessionTimeSave (autoSaveInterval, 0) ; /* saves unless recently done */
  signal (SIGALRM, checkLife) ;
  if (nActiveClients) 
    {
      if (clientTimeOut)
	alarm ((unsigned int) clientTimeOut) ;
      dying = FALSE ;
    } 
  else 
    { fprintf(out,"\nNo more active clients, starting to count-down: %d s // %s\n", 
	      serverTimeOut, timeShowNow()) ;
      if (serverTimeOut)
	alarm ((unsigned int)serverTimeOut) ;
      dying = TRUE ;
    }
  last_t = timeNow () ; /* saving may have taken several minutes */
}

/****************** Password file **********************/

static void setPassFile (LOOK look, 
			 int magicRead, int magicWrite,
			 int clientId)
{ 
  static char *magicWriteDir = 0 ;
  static char *magicReadDir = 0 ;
  static BOOL  firstPass = TRUE ;
  static int   public = 0 ;
  FILE *f = 0 ;
  char *cp ;
  int   level ;

  if (firstPass)
    { firstPass = FALSE ;
      cp = sessionFilName ("wspec/server", "wrm", "r") ;
      if (!cp ||
	  !(f = filopen (cp, 0, "r")))
	{ 
	  if (! cp) cp="no file";
	  fprintf (out, "Server cannot open server.wrm : %s", cp) ;
	  messcrash ("Server cannot open server.wrm : %s", cp) ;
	}
      level = freesetfile(f,0) ;
      while (freecard(level))
	{ cp = freeword() ;
	  if (cp)
	    { 
	      if (!strcmp (cp,"WRITE_ACCESS_DIRECTORY") &&
		  (cp = freeword()) )
	        magicWriteDir = strnew (cp, 0) ;
	      else if (!strcmp (cp,"READ_ACCESS_DIRECTORY") &&
		(cp = freeword()))
		{ 
		  if (!strcmp(cp, "PUBLIC"))
		    public = 2 ;
		  else if (!strcmp(cp, "RESTRICTED"))
		    public = 0 ;
		  else
		    { public = 1 ; magicReadDir = strnew (cp, 0) ; }
		}
	    }
	}
      freeclose (level) ;
    }
  
  if (!public && !magicWriteDir)
    { fprintf (out, "Nobody can read or write, please set in the file wspec/server.wrm\n"
	       "the variables READ_ACCESS_DIRECTORY and WRITE_ACCESS_DIRECTORY") ;
      messcrash ("Nobody can read or write, please set in the file wspec/server.wrm\n"
		 "the variables READ_ACCESS_DIRECTORY and WRITE_ACCESS_DIRECTORY") ;
    }

  look->magicReadFileName = 0 ;
  if (magicReadDir)
    { cp = filName (messprintf("%s/client.R.%d",magicReadDir,clientId), 0, "w") ;
      if (!cp ||
	  !(f = filopen(cp, 0, "w")))
	{ fprintf (out, "Server cannot create a magic file in the READ_ACCESS_DIRECTORY %s", magicReadDir) ;
	  messerror ("Server cannot create a magic file in the READ_ACCESS_DIRECTORY %s", magicReadDir) ; 
	  look->magicReadFileName = strnew ("totojunk", 0) ;
	}
      else
	{
	  fprintf(f, "%d\n", magicRead)  ;
	  filclose (f) ;
	  /*      system (messprintf ("chmod 640 %s", cp)) ; */
	  look->magicReadFileName = strnew (cp, 0) ;
	}
    }

  look->magicWriteFileName = 0 ;
  if (magicWriteDir)
    { cp = filName (messprintf("%s/client.W.%d",magicWriteDir,clientId), 0, "w") ;
      if (!cp ||
	  !(f = filopen(cp, 0, "w")))
	{ fprintf (out, "Server cannot create a magic file in the WRITE_ACCESS_DIRECTORY %s", magicWriteDir) ;
	  messerror ("Server cannot create a magic file in the WRITE_ACCESS_DIRECTORY %s", magicWriteDir) ;
	  look->magicWriteFileName = strnew ("totojunk", 0) ;
	}
      else
	{
	  fprintf(f, "%d\n", magicWrite)  ;
	  filclose (f) ;
	  look->magicWriteFileName = strnew (cp, 0) ;
	}
    }

  look->public = public ;
}

static void removePassFile (LOOK look)
{
  if ( look->magicReadFileName)
    filremove ( look->magicReadFileName, 0) ;
  if ( look->magicWriteFileName)
    filremove ( look->magicWriteFileName, 0) ;
}

/****************** MAIN LOOP **********************/

 /* By convention, each reply starts with a message ending on /n
    which will be deleted if aceclient runs silently rather than verbose 
    */

/* DWB - int maxChar, int *encore added to parameter list for processQueries */

Stack processQueries(int *clientIdp, int *clientMagic, char *question, int maxChar, int *encore)
{ KEY option ;
  char *cp ;
  static Stack s = 0 ;
  static int nClients = 0 , nTransactions = 0 ;
  int level, magic1, magic2, magic3, magicRead, magicWrite ;
  LOOK look ;

  if (debug) messerror("ProcessQueries: %s\n", question ? question : "null") ;

  /* DWB - don't let client ask for more than server limit */
  if ((maxKbytes > 0) && (maxChar > (maxKbytes*1024)))
    maxChar = maxKbytes*1024;
/* security system:
   the client is given magic1
   but should later return magic3
   which he contruct by reading the groupreadable file 
   containing magic2
*/
  /* JC1 switch of timeOut alarm , 
     we switch it back on at the end of this routine,
     in the call to checkLife */
  alarm((unsigned int)0);

  s = stackReCreate(s, 20000) ;
  nTransactions++ ;

  if (!*clientIdp)
    { 
      if (refuseNewClient)
	{ pushText (s, "Too late, I am dyyyiiing") ;
	  checkLife(0) ;
	  return s ;
	}

      look = (LOOK) messalloc(sizeof(struct LOOKSTUFF)) ;
      *clientIdp = look->clientId = ++nClients ;
      magic1 = assInt(look) ; 
      magic2 =  randint() ^ ((char*)(&cp) - (char*)0) ;
      magic3 = randint() ^  (int)(timeNow()) ;
      if (magic1 < 0) magic1 *= -1 ;
      if (magic2 < 0) magic2 *= -1 ;
      if (magic3 < 0) magic3 *= -1 ;
      if (!magic1) magic1 = 1 ;
      if (!magic2) magic2 = 1 ; 
      if (!magic3) magic3 = 1 ;
      *clientMagic = magic1 ;
      magicRead  = magic1 * magic2 % 73256171 ;
      magicWrite = magic1 * magic3 % 43532334 ;
      setPassFile (look, magic2, magic3, *clientIdp) ;
if (debug) messerror("magic1 = %d, magic2 = %d, magic3 = %d, magicRead = %d, magicWrite = %d\n",
	magic1, magic2, magic3, magicRead, magicWrite) ;
      look->readMagic = look->writeMagic = 0 ;
      if (look->magicReadFileName) look->readMagic = magicRead ;
      if (look->magicWriteFileName) look->writeMagic = magicWrite ;

      cp = (char *)0 + nClients ;
      assInsert (clientAss, cp, look) ;
      fprintf (out, "\n%sNew Client %d, %d active clients\n", timeShowNow(), 
	      *clientIdp, ++nActiveClients) ;

      if (look->magicWriteFileName) 
	{ look->mayWrite = TRUE ;
	  pushText (s, look->magicWriteFileName) ;
	}
      else
	{ look->mayWrite = FALSE ;
	  pushText (s, "NON_WRITABLE") ;
	}
      catText (s, " ") ;
      if (look->magicReadFileName) catText (s, look->magicReadFileName) ;
      else catText(s, look->public ? "PUBLIC" : "RESTRICTED") ;
      /* store time of last access. Used for the killing of 
         dead clients */
      look->lastAccess = timeNow() ;
#ifndef GIFACESERVER 
      look->aceCommand = (void*) aceCommandDoCreate (9,0,s) ;
#else
#ifdef ACEMBLY
      look->aceCommand = (void*) aceCommandDoCreate (29,0,s) ;
#else
      look->aceCommand = (void*) aceCommandDoCreate (25,0,s) ;
#endif
#endif
      checkLife(0);
if (debug) messerror(stackText(s,0)) ;
      return s ;
    }     
  else if ((cp = (char *)0 + *clientIdp), 
	   !assFind (clientAss, cp, &look))
    { messout("Aceserver received a bad client id = %d", *clientIdp) ;
      pushText(s, "//! Unauthorised access_1 : closing connection\n") ;
      checkLife(0);
      return 0 ;
    }

  if (debug) messerror ("ProcessQueries passed auth: %s\n", question ? question : "null") ;

  look->lastAccess = timeNow (); /* last access, to kill  dead clients */

  removePassFile (look) ; /* single guess */

  switch(look->public)
    {
    case 0: /* no read_access_dir, you must have write access */
     if ( look->mayWrite && look->writeMagic &&
	  *clientMagic == look->writeMagic)   /* ok, you can have read/write access */
	break ;
     goto nasty ;

    case 1:  /* there is a read_access_dir, you must read it */
      if ( look->mayWrite && look->writeMagic &&
	   *clientMagic == look->writeMagic)   /* ok, you can have read/write access */
	break ;
      if (*clientMagic == look->readMagic)   /* ok, you can have read access */
	 {
	   look->mayWrite = FALSE ;              /* mouse trap */
	   aceCommandNoWrite(look->aceCommand) ;
	   break ;
	 }
      goto nasty ;
      
    case 2: /* read access is public */
      if ( look->mayWrite && look->writeMagic &&
	   *clientMagic == look->writeMagic)   /* ok, you can write access */
	break ;
      look->mayWrite = FALSE ; look->public = 2 ;  /* mouse trap */
      aceCommandNoWrite(look->aceCommand) ;
      break ;
    }
  
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
	 
  fprintf (out, "%s Client %d (%d active). Query: %s\n", timeShowNow() ,  *clientIdp, nActiveClients, *encore ? "encore" : question) ;
  if (strncmp (question, "wspec",5))
    {
      level = freesettext (question, "") ;
      freespecial ("\"\t\\/@%") ;  /* forbid sub shells */
      option = aceCommandDoExecute (look->aceCommand, level, 
				    -1, *encore ? look->lastCommand : 0, maxChar) ;
      freeclose (level) ;
      *encore = 0 ;
      look->lastAccess = timeNow() ;
      
      switch (option)
	{
	case 'q':
	  clientDestroy(look) ;
	  *clientIdp = 0 ;
	  pushText(s, "// A bientot") ;
	  break ;
	  
	case 'm':
	case 'D':
	case 'B':
	case 'M':
	  look->lastCommand = option ;
	  *encore = 2 ;  /* automatic looping */
	  break ;
	  
	case 'W':  /* who */
	  if (!look->mayWrite)
	    { pushText (s,"// Sorry, you do not have write access.\n") ;
	    break ;
	    }
	  pushText(s, messprintf("// %d active Clients, %d transactions so far",
				 nActiveClients, nTransactions)) ;
	  break ;
	  
	case 'U':  /* shUtdown [now] */
	  if (!look->mayWrite)
	    { pushText (s,"// Sorry, you do not have write access.\n") ;
	    break ;
	    }
	  refuseNewClient = TRUE ;
	  cp = freeword() ;
	  if (cp && !strcasecmp(cp, "now")) 
	    wormClose () ;
	  else
	    {
	      pushText (s,"// The server will now refuse any new connection \n") ;
	      catText (s, "// and shutdown when the last client quits\n") ;
	      catText(s,messprintf("// Now %d active Client(s)\n",
				   nActiveClients)) ;
	    }
	  break ;
	}
    }
  else if (!strcmp (question, "wspec1"))
    {  /* export wspec 
	  does not work, i do not know how to untar at the client side
      catBinary (s, arrp (wspecTar, 0, char), arrayMax(wspecTar)) ;
      */
    }
  else if (!strcmp (question, "wspec"))
    {  /* export wspec */
      if (!wspecTarEncoded)
	copy_wspec (TRUE) ;
      /*
	i do not know how to recosntruct the non encoded file at the other hand, sorry
	copy_wspec (FALSE) ;
      */
      catBinary (s, arrp (wspecTarEncoded, 0, char), arrayMax(wspecTarEncoded)) ;
    }

/*  if (oldAnswer != stackText(s, 0))
    messout("Relocation: length(answer) = %d",  stackMark(s)) ;
    oldAnswer = stackText(s, 0) ;

*/
  
fin:
  checkLife(0);
  return s ;
  
nasty:
  fprintf (out,
	   "%s Unauthorised access by client %d, magic = %d != look->writeMagic = %d, look->readMagic = %d\n"
	   , timeShowNow(),
	   *clientIdp, *clientMagic, look->writeMagic, look->readMagic) ;
  clientDestroy (look) ;
  *clientIdp = 0 ; 
  pushText(s, "//! Unauthorised access_2: closing connection\n") ;
  goto fin ;
}

/**************************************************/
/********* ACEDB non graphic Query Server  ********/
/**************************************************/

static void copy_wspec (BOOL encode)
{
  FILE *f = 0 ;
  char c, *cp ;
  int n = 0 ;
  Array aa = arrayCreate (20000, char) ; 

  cp = filName("wspec",0,"r") ;
  if (!cp) 
    messcrash ("Server cannot open wspec, sorry") ;
  
  freeOut("\n") ;
  if (encode)
    f = popen(messprintf("cd %s ; tar chf - %s %s | uuencode server.wspec.tar",
			 cp, "cachesize.wrm cachesize.wrm database.wrm displays.wrm",
			 "layout.wrm models.wrm options.wrm psfonts.wrm subclasses.wrm xfonts.wrm"),"r") ;
  else
    f = popen(messprintf("cd %s ; tar chf - %s %s",
			 cp, "cachesize.wrm cachesize.wrm database.wrm displays.wrm",
			 "layout.wrm models.wrm options.wrm psfonts.wrm subclasses.wrm xfonts.wrm"),"r") ;
  while (c = getc(f), c  != (char)EOF)
    array (aa, n++, char) = c ;
  pclose (f) ;
  if (encode)
    wspecTarEncoded = aa ;
  else
    wspecTar = aa ;
}



static void daemonize(void)
{
#ifdef JUNK 
  /* on alpha htis creates an infinity of server process, bad ! */
 int pid ; /* detach from terminal */

  pid = fork() ;
 if (pid < 0)
   exit (0) ; /* messcrash("Cannot fork") ; */
 if (pid > 0)
   exit (0) ; /* parent suicides */
#endif

 setsid () ;   /* reset neutral environment */
 /* chdir ("/") ;  
	How can the server possibly work with chdir "/" here?  
	filopen() can't find my files, so how does it find the
	data files?
 */
 umask (0) ;
}

static void openOutFile (void)
{ char buf [12*1024] ;
  strcpy (buf, sessionFilName ("server", "log", 0)) ;
        
  out = filopen (buf, 0, "a") ;
  if (!out)
    out = filopen (buf, 0, "w") ;
  if (!out)
    out = stderr ;
  setbuf (out, NULL) ;
  /* setbuf (stdout, NULL) ; this flat crashes under daemon on new linux */
}


/* Defines a routine to return the compile date of this file. */
UT_MAKE_GETCOMPILEDATEROUTINE()

/************************************************************/

FILE *get_log_file_fp()
{
return out;
}

int main (int argc, char **argv)
{ 
  u_long p ;
  int n, t ;
  char *cp ;
  /* Declarations for inetd detection LDS 7/2/98 */
  struct sockaddr_in saddr;
  socklen_t asize = sizeof (saddr);


#ifdef GIFACESERVER
  extern void (*gifEntry)(KEYSET,int,BOOL) ;	/* entry point in command.c */
  extern void gifControl (KEYSET ks, int level, BOOL isI)   ;

  gifEntry = gifControl ; 

  freeinit () ;			/* must come before graphInit */

  /* Initialise the graphics package and its interface for acedb.            */
  acedbAppGraphInit(&argc, argv) ;

  messErrorInit("gifaceserver");   /* Record program name */
#else
  messErrorInit("aceserver");   /* Record program name */
#endif


  out = stderr ; /* this static cannot be initialized at file scope in WIN32? */

  /* This code determines whether we are running from inetd.
     It relies on the fact that if we are running under inetd,
     file descriptor 0 will be a socket.
     LDS 7/2/98 */
  if (getsockname(0, (struct sockaddr *)&saddr, &asize) == 0)
    isDaemon = TRUE ;
  else
    isDaemon = FALSE ;

#ifdef SGI
  /* in sgi case, i call rpc_init etc with diferent param (zero)
   * in the daemon case, 
   */
  isSGI = TRUE ;
#endif

  n = 2;
  if (argc > n && (!strcmp(argv[n],"-port")))
    { /* ignore, optional -port */
      n++;
    }
  if (argc > n && (sscanf(argv[n],"%lu",&p) == 1))
    { port = p; /* port = 2000700 ; */
      n++;
    }
  if (argc > n)
    {
      cp = argv[n] ;
      if (sscanf(cp,"%d",&t) == 1)
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

  if (debug) messerror ("AceServer: m1") ;
  isInteractive = TRUE ;
  clientAss  = assCreate() ;
  signal (SIGALRM, checkLife) ;
  if (debug) messerror ("AceServer: m2") ;
  /* start timer countDown */
  if (serverTimeOut)
    alarm((unsigned int)serverTimeOut);

  openOutFile() ;
  if (debug) messerror ("AceServer: m3") ;
  cp = getenv ("HOST") ; if (!cp) cp = "(unknown)" ;
  fprintf(out,"\n\n#### Server starts %s\n", timeShowNow()) ;
  fprintf(out,"#### host=%s  port=%lu  ACEDB=%s\n",  
	  cp, port, sessionFilName(0,0,0)) ;
  fprintf(out,"#### clientTimeout=%d serverTimeout=%d maxKbytes=%d autoSaveInterval=%d\n",
	  clientTimeOut, serverTimeOut, maxKbytes, autoSaveInterval); 

#if !defined(DCESERVER)
  signal (SIGPIPE, sigPipeError) ;
#endif

  if (debug) messerror ("AceServer: m4") ;
  daemonize() ;
  if (debug) messerror ("AceServer: m5") ;
  wait_for_client(port, isSGI && isDaemon) ;
  messcrash("Acedb network server error: main loop returned") ;

  return(EXIT_FAILURE) ;			  /* We should not reach here... */
}

/**************************************************/
/**************************************************/
 
 
 
 
