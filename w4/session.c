/*  File: session.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
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
 *	Jean Thierry-Mieg (CRBM du CNR, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 *   Session control, passwords, write access, locks etc.
 * Exported functions:
 * HISTORY:
 * Last edited: Apr 19 12:10 1999 (fw)
 * * Dec  3 14:43 1998 (edgrif): Change calls to new interface to aceversion.
 * * Nov 26 11:21 1998 (fw): removed waction stuff
 * * Oct 22 11:54 1998 (edgrif): Add action.h for init/close_action()
 *              THERE IS A PROBLEM WITH CLASHES WITH disk.h for BP..sigh.
 * * Oct 21 14:01 1998 (edgrif): Remove use of isGraphics to make graphOut
 *              do a printf, not needed now.
 * * Sep 11 09:56 1998 (rbsk) need strnew for filGetFullPath assigned to ACEDIR
 * * Sep 11 09:56 1998 (edgrif): Replace calls to messcrash with messExit
 *              where the problem is a user error.
 * * Sep 11 09:43 1998 (edgrif): Replace messcrash calls with messExit
 *              where this is appropriate. Remove from messages any mention
 *              of "setenv ACEDB" as this is now deprecated.
 * * Sep  3 11:14 1998 (edgrif): Remove "ACEDB" from fprintf string in
 *              acedbCrash, this is provided by uMessCrash now.
 * * Aug 21 14:49 1998 (rd): moved ruid/euid back here, along with getLogin
 *      and dump/crash handlers
 * * missing date/name
 *      - Moved the write access messages to the locking subroutines to cover
 *      - all uses. And removed time etc which is now automatically done.
 * * Sep  5 10:44 1997 (rd)
 *	-	WIN32 runs sessionFilName() like UNIX when defined(NON_GRAPHIC)
 *		WinTace user should specify ACEDB directory on the commandline?
 * * Jun  9 14:50 1996 (rd)
 *	-	extern uid_t euid, ruid moved to session.h
 * * Jun  6 17:41 1996 (rd)
 * * Jun  3 19:00 1996 (rd)
 * * Feb  2 11:46 1992 (mieg): interactive session control
 * * Jan 24 11:20 1992 (mieg): session fusion.
 * Created: Fri Jan 24 11:19:27 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: session.c,v 1.34 2020/05/30 16:50:32 mieg Exp $ */
/* #define CHRONO  */

#define DEF_BP /* avoid typedefing it as void* in disk.h */
typedef struct block *BLOCKP ;
typedef BLOCKP BP ;

#include "acedb.h"

#include "mytime.h"
#include <signal.h>           /* for autosave */
#include "freeout.h"
#include "lex.h"
#include "systags.h"
#include "sysclass.h"
#include "../wh/byteswap.h"		/* for calls in swapSuperBlock() */
#include "block.h"
#include "session.h"
#include "a.h"
#include "bs.h"
#include "disk.h"
#include "adisk.h"
#include "disk_.h"              /* defines BLOCKHEADER */
#include "session_.h"		/* uses BLOCKHEADER */
#include "lex_sess_.h"		/* lexSessionStart etc.. */

#include "pick.h"
#include "pref.h"		/* prefInit */
#include "model.h"		/* for readModels() */
#include "command.h"		/* for tStatus */
#include "help.h"		/* for helpSetDir during initilisation */
#include "cdna.h"               /* for cDNAAlignInit */
#include "chrono.h"
extern void CDTShutDown (void) ;
extern void scShutDown (void) ;
extern void condShutDown (void) ;

/***************** system includes **********************/

#include <errno.h>

#if  defined(MACINTOSH) || defined(MAC_X_BAD) || defined(WIN32) || defined(CYGWIN)

#define F_ULOCK 0
#define F_TLOCK 0

static BOOL suppressLock = TRUE ; /* don't use the lockf(3) mechanism */
BOOL lockf ( int lockFD, int mode1, int mode2 ) { return FALSE ; }

/* Stub functions for user ID's in non-multiuser systems
   Note: A WIN32 implementation for Windows NT may be possible */

#if !defined (MAC_X_BAD) && !defined(CYGWIN)
typedef int pid_t;
void	seteuid(uid_t uid) { return ; }


uid_t	getuid() { return 1 ; }
uid_t	geteuid() { return 1 ; }
int	gethostname(char *b, int l) { *b = 0; }

pid_t	getpid (void) { return 0; }
#endif
#else  /* UNIX */

#include <unistd.h>     /* extern lockf(); */
static BOOL suppressLock = FALSE ; /* use the lockf(3) mechanism */
/* on HP, for some reason, i did not get these which are innsys/unistd.h */
#ifndef F_TLOCK
#define F_TLOCK       2       /* Test and lock a region for exclusive use */
#endif

#ifndef F_ULOCK
#define F_ULOCK       0        /* Unlock a previously locked region */
#endif

#endif /* not MAC or WIN32 */

/************************************************************/

#ifdef IBM

struct passwd {
  char    *pw_name;
  char    *pw_passwd;
  long   pw_uid;
  long   pw_gid;
  char    *pw_gecos;
  char    *pw_dir;
  char    *pw_shell;
} ;
extern struct passwd *getpwuid (long uid) ;
extern struct passwd *getpwnam (const char *name) ;

#else  /* IBM */

#include <pwd.h>

#endif /* !IBM */

/************************************************************/

       /* globals declared in session.h */

/* uid_t	ruid = -1, euid = -1; now in utils.c */
SESSION thisSession = { 0 };
BOOL	swapData = FALSE ;

/************************************************************/
mytime_t sessionUserSession2Time(KEY userSession) { return 0 ; }

static BOOL dropWriteAccess = FALSE ;	/* once TRUE, write access
				   is forbidden and can't be regained
				   this is an odd behaviour */

static int saveMagic = 0 ; /* read/write verification */

 /* number of sessions to keep living before a session expires
    (if non-permanent and not readlocked) Default is 2, must be >= 1 */
static int sessionKeepAlive;	/* init'd in sessionInit */

static int writeAccess = 0 ;
static BLOCK theSuperBlock ;

static BOOL debug = FALSE;	    /* more log-file output */
static BOOL debugReadLock = FALSE ; /* to see the mesages failed to rename readlock file */

/************************************************************/


/* if set, externalServer will be called once on every new object, and
   always in query(0,..) and in grep 
   Registration of this pointer done via xCientInit
*/

KEYSET (*externalServer) (KEY key, char *query, char *grep, BOOL getNeighbours) ;
void   (*externalSaver) (KEY key, int action) ;
void   (*externalAlias) (KEY key, char *oldname, BOOL keepOldName) ;

void   (*externalCommit) (void) ;
BOOL   (*externalServerState) (int nn) ;
char*  (*externalServerName) (void) ;

/************************************************************/

static void sessionDoInit (void) ; /* session.c's own initialisation */
static char *getConfiguredDatabaseName (void) ;
static int  getConfiguredKeepAlive (void) ;

static int  sessionFindSon (Array sessionTree, KEY father, KEY *sonp) ;
static BOOL sessionDatabaseWritable (BOOL doTalk);
static BOOL checkSessionNumber(int type);
static BOOL sessionCheckUserName (BOOL doTalk);

static BOOL readlockCreateDir (void);
static BOOL readlockCreateFile (void);
static BOOL readlockUpdateFile (int new_session_num);
static Array readlockGetSessionNumbers (void);

/************************************************************/

/* this mechanism is used to register a function that redraws the 
   window layout or changes menus etc. when the write access changes */

static VoidRoutine writeAccessChangeRoutine = 0;

VoidRoutine writeAccessChangeRegister (VoidRoutine func)
{ 
  VoidRoutine old = writeAccessChangeRoutine ; 
  writeAccessChangeRoutine = func ; 
  return old ;
} /* writeAccessChangeRegister */

/*************************************************************/

/* this mechanism is used to reflect changes in the session structure
   in a graphical representation of session objects */

static VoidRoutine sessionChangeRoutine = 0;

VoidRoutine sessionChangeRegister (VoidRoutine func)
{ 
  VoidRoutine old = sessionChangeRoutine ; 
  sessionChangeRoutine = func ; 
  return old ;
} /* sessionChangeRegister */

/*************************************************************/

#ifndef ACEDB5
#ifdef ACEDB4
void swapSuperBlock(BLOCK *bp)	/* used also by disknew.c */
{
  bp->h.disk = swapDISK(bp->h.disk);
  bp->h.nextdisk = swapDISK(bp->h.nextdisk);
  bp->h.session = swapInt(bp->h.session);
  bp->h.key = swapKEY(bp->h.key);
  bp->gAddress = swapDISK(bp->gAddress);
  bp->mainRelease = swapInt(bp->mainRelease);
  bp->subDataRelease = swapInt(bp->subDataRelease);
  bp->subCodeRelease = swapInt(bp->subCodeRelease);
} /* swapSuperBlock */
#endif
#endif

/*************************************************************/

static void sessionSetPath (const char *ace_path)
{ 
#ifndef MACINTOSH
  if (ace_path)
    { filAddPath (ace_path) ;
    }
  if (getenv ("ACEDB"))			/* $ACEDB */
    { filAddPath (getenv ("ACEDB")) ;
    }
  if (getenv ("ACEDB_COMMON"))
    filAddDir (getenv("ACEDB_COMMON")) ;
#ifndef WIN32
  else
    filAddDir ("/usr/local/acedb") ;
#endif
  filAddDir (".") ;
#endif /* MACINTOSH */
} /* sessionSetPath */

/*************************************************************/

char *sessionFilName (char *name, const char *ending, const char *spec)
     /* general public function */
{
#ifdef MACINTOSH
  return filName (name, ending, spec) ;
#else  /* !MACTINTOSH */
  static Stack s = 0 ;
  static char *ACEDIR ;
  char *cp = 0 ;

  if (!s) 
    { 
      s = stackCreate (1024) ;
#ifdef WIN32
    retry:
#endif /* WIN32 */
      cp = 0 ;
      while (!cp)
	{
	  cp = filStrictName ("database", 0, "rd") ;
	  if (!cp)
	    {
#if (defined(WIN32) && !defined(NON_GRAPHIC)) /* graphical WIN32 only */
	      messout("Home database directory not found?\n"
		      "You must give a parent directory with subdirectory \"database\"\n"
		      "as a command line \"ACEDB=<pathname>\" argument\n"
		      "or by setting up an ACEDB database profile.") ;
	      getSessionParameters() ; /* see win32lib.cpp - can change path */
	      sessionSetPath(0) ;	/* see above */
	      goto retry ;
#else  /* all other UNIX or non-graphic WIN32 build */
	      messExit ("Home database directory not found - you must give a "
			"parent directory with subdirectories database/ and wspec/ "
			"as a command line argument.") ;
#endif /* UNIX or nongraph WIN32 */
	    }
	}
      cp[strlen(cp)-9] = 0 ;	/* remove "/database" */
      ACEDIR = filGetFullPath (cp,0) ;
      if (!ACEDIR)
	messExit ("Can not find absolute path name for ACEDB home directory") ;
      if (!filName (messprintf ("%s/wspec", ACEDIR), 0, "rd"))
	{
#if (defined(WIN32) && !defined(NON_GRAPHIC)) /* graphical WIN32 only */
	  messout("Home database directory not found?\n"
		  "You must give a parent directory with subdirectory \"database\"\n"
		  "as a command line \"ACEDB=<pathname>\" argument\n"
		  "or by setting up an ACEDB database profile.") ;
	  getSessionParameters() ; /* see win32lib.cpp - can change path */
	  sessionSetPath(0) ;	/* see above */
	  goto retry ;
#else  /* all other UNIX or non-graphic WIN32 build */
	  messExit ("Home database directory not found - you must give a "
		    "parent directory with subdirectories database/ and wspec/ "
		    "as a command line argument.") ;
#endif /* UNIX or nongraph WIN32 */
	}

      /* now set location of helpdirectory */
      helpSetDir (messprintf ("%s/whelp", ACEDIR));

      if (!getenv ("ACEDB_NO_BANNER"))
	printf ("// Database directory: %s\n", ACEDIR) ;
          /* // so that it does not pollute an acefile output when using acelib */
    }

  if (!name)
    return ACEDIR ; 
  
  stackClear(s) ;
  if (ACEDIR)			/* might still be NULL */
    pushText (s, ACEDIR) ;
  catText (s, "/") ;
  catText (s, name) ;
  if (ending)
    { catText (s, ".") ;
      catText (s, ending) ;
    }
  return filStrictName (stackText (s,0), 0, spec) ;
#endif /* !MACINTOSH */
} /* sessionFilName */

/*************************************************************/
/****************** getLogin() *******************************/

/**** getLogin() for real/effective user names 
***** original by 03.05.1992 kocab DKFZ  
*****/
#ifdef JUNK
now in utils.c
char *getLogin (BOOL isReal)
     /* general public function */
{
  /* isReal for real or effective can't fail */

#if defined(MACINTOSH)
  return "acedb" ;
#elif defined (WIN32)
  char *name = NULL;
  if ((name = getenv("USERNAME")) != NULL )
    return name ;   /* Windows NT has usernames */
  else
    return "acedb" ;
#else  /* all UNIX */
  static char *rname = 0 ;
  static char *ename = 0 ;

/* RD 980417: changed so getlogin() was last resort.  It can return
   "root" inappropriately where getpwuid(ruid)->pw_name returns the
   right answer (ddts problem SANgc02099)
*/

  if (!ename)
    {
      if (ruid == -1)
	{ ruid = getuid() ;
	euid = geteuid () ;
	}
      
      if (!rname)
	if (ruid)
	  if (getpwuid(ruid))
	    rname = strnew(getpwuid(ruid)->pw_name, 0) ;
      if (!rname)
	if (getlogin())
	  rname = strnew(getlogin(), 0) ;
      if (!rname)
	rname = "getLoginFailed" ;
      
      if (!ename)
	if (getpwuid(euid))
	  ename = strnew(getpwuid(euid)->pw_name, 0) ;
      if (!ename)
	ename = rname ;
      /* mieg: jan5 99, new problem
	 on my dec alpha, if i run as a daemon, getLogin fails totally
	 printf ("rname=%s ename =%s",rname?rname:"null",ename ? ename:"null") ; 
	 */
    }
  return isReal ? rname : ename ;
#endif /* UNIX */
} /* getLogin */
#endif /* JUNK */
/*************************************************************/
/*************************************************************/

static char *messcrashbuf = 0 ;	/* a memory buffer that is
				   allocated at initialization time
				   and free'd when acedb crashes, to
				   gain some extra memory for the crash
				   recovery procedures */
/******************************************************************/

void simpleCrashRoutine (void)
     /* The least amount of cleanup work we need to do,
	this function is the default for the messcrashroutine,
	more elaborate functions will typically save the current
	session or close windows etc. */
{ 
  sessionReleaseWriteAccess() ;	/* must remove lock file */
  readlockDeleteFile ();	/* release readlock */

  filtmpcleanup () ;		/* deletes any temporary files made */

  return;
} /* simpleCrashRoutine */

/******************************************************************/

VoidRoutine messcrashroutine = 0; /* as soon as necessary it will become
				     simpleCrashRoutine, the minimal
				     acedb clean-up routine. */
/* moved here from w1/messubs.c */
/* The meaning of this function has changed now.
   it is now only used for acedb code, where acedbCrash()
   is messCrashRegister'd with the message package.
   acedbCrash will call this messcrashroutine (init'd by aceInit)
   to allow the to perform specialized clean-up tasks, like
   release locks etc. */
/* keep using this name for now because used frequently in code */

static FILE *dumpFile = 0 ;

/******************************************************************/

static void acedbCrash (const char *text)
/* this function is called by messcrash to do an acedb-specific 
   crash-exit. It logs the error to the log file and calls the user
   defined cleanup routine messcrashroutine, which should most
   definitively call aceQuit to clean upo read/write locks */
{
  if (messcrashbuf)
    messfree (messcrashbuf) ;	/* frees up space that might be useful */

	/* if memory failure report usage info */

  if (!strncmp (text, "FATAL ERROR: Memory allocation failure", 38))
    { 
      int level ;

      if (dumpFile)		/* messcrash() should have opened this */
	{ level = freeOutSetFile (dumpFile) ;
	  tStatus () ;
	  freeOutClose (level) ; 
	  fflush (dumpFile) ;
	}
      else		   /* e.g. if can't write database/ directory */
	tStatus () ;
    }

  fprintf (stderr, "%s\n", text) ;		  /* for good measure */
  messout (strnew (text, 0)) ;			  /* strnew to prevent
						     messbuf tangle */


  fprintf (stderr, "\n// Sorry for the crash, please report it, "
	                 "if it is reproducible\n") ;

  /* the messcrashroutine should call aceQuit, to clean up,
     it may give the user a chance to save their work so far */
  if (messcrashroutine)
    {
      VoidRoutine the_routine = messcrashroutine ;
      messcrashroutine = 0;	/* avoids blind recursion */
      (*the_routine)();		/* call user cleanup routine */
    }

  invokeDebugger ();		/* after read/write locks are released */

  exit (EXIT_FAILURE);
} /* acedbCrash */

/************************************************************/

static void acedbExit (const char *text)
/* this is similar to acedbCrash, apart from the status report
   and the log-file reports, which make an error look really bad.
   This function is supposed to be a more graceful exit from acedb */
{
  fprintf (stderr, "%s\n", text) ;		  /* for good measure */
  messout (strnew (text, 0)) ;			  /* strnew to prevent
						     messbuf tangle */

  fprintf (stderr, "\n// Exit\n") ;

  /* the messcrashroutine should call aceQuit, to clean up,
     it may give the user a chance to save their work so far */
  if (messcrashroutine)
    {
      VoidRoutine the_routine = messcrashroutine ;
      messcrashroutine = 0;	/* avoids blind recursion */
      (*the_routine)();		/* call user cleanup routine */
    }

  invokeDebugger ();		/* after read/write locks are released */

  exit (EXIT_FAILURE);		/* although this is a graceful exit
				   something still went wrong */
} /* acedbExit */

/************************************************************/

static void acedbLog (const char *text)
     /* registered to be the messDumpRoutine for messdump */
{
  static FILE *dumpFile = 0 ;
  static char host[100];
  static int pid;
  static BOOL cantOpenFile = FALSE ;
  char timeBuf[25] ;

  if (cantOpenFile)
    return ;

  if (euid != ruid)		/* log file is in database/ directory */
    seteuid(euid);

  if (!dumpFile)
    {
      char *user = getLogin(TRUE) ;
      gethostname (host, 100) ;
      pid = getpid () ; 

      dumpFile = fopen (sessionFilName ("database/log", "wrm", 0), "a") ;
      if (!dumpFile)
	{ cantOpenFile = TRUE ;
	  return ;
	}

      fprintf (dumpFile, "\nNew start User:%s,  %s,  %s :: %s\n"
	       , user, aceGetVersionString(), aceGetLinkDateString()
	       , timeShow(timeNow(), timeBuf, 25)) ;
    }

				/* prepend time, host, pid to format */
  fprintf (dumpFile, "%s %s %d\t", timeShow(timeNow(), timeBuf, 25), host, pid) ;
  fprintf (dumpFile, "%s\n", text) ;
  fflush (dumpFile) ;

  if (euid != ruid)
    seteuid(ruid);
} /* acedbLog */

/************************************************************/

static void sessionInit2 (const char* ace_path)
     /* general public function */
{
  extern void qvInit(void) ;	/* in quovadis.c */

/* only ever called from aceInit()
   or in programs that don't use the acelib API yet
   it is only called in the beginning from the main function */
  KEY lastSession = 0;

#if !(defined(MACINTOSH) || defined(WIN32))
  euid = geteuid();
  ruid = getuid();
  if (euid != ruid)
  seteuid(ruid);

  if (euid != ruid && !getenv("ACEDB_NO_BANNER"))
    { printf("Acedb is running for user %s.\n", getLogin (TRUE));
      printf("Acedb is inheriting database write access from user %s.\n", 
	     getLogin (FALSE));
    }

  if (getgid() != getegid())
     messerror("WARNING. Your acedb executable has a permission "
	       "flag that inherits the group-id on execution (setgid). "
	       "This is pointless, and will stop acedb from spawning "
	       "external programs on systems with dynamic linking. "
               "Use executables with a set-user-id flag instead.");

#endif

  /* Macro KEYKEY relies on KEY beeing 4 bytes long  */
  /* Furthermore, it may be assumed implicitly elsewhere */

  if (sizeof(KEY) != 4)
    messcrash ("KEY size is %d - it must be 4 - rethink and recompile!",
	        sizeof(KEY)) ;

  if (__Global_Bat != KEYMAKE(1,1)) /* something is really wrong */
    messcrash 
      ("Hand edit systags.h because KEYMAKE(1,1) = %d != __Global_Bat",
       KEYMAKE(1,1)) ; 


  /* init these freeXXX first, they belong to the utilities-library,
     do all the acedb specific inits after that */
  sessionSetPath (ace_path) ;
  freeinit () ;
  freeOutInit () ;  


  
  /* prompt the system to check whether the file
     exists. If the database directory ace_path (just set above)
     is not recognised, the program will exit */
  sessionFilName("wspec/models.wrm", 0, "rd");
  
  /* this is panic-memory, we claim some memory here which 
     we'll never use. In case the program messcrash'es we free this
     chunk of storage, so hopefully we can exit gracefully and
     save our work before we have to exit */
  messcrashbuf = messalloc (12*1024) ;

  /* now register all the acedb specific crash/recovery routines */
  messCrashRegister (acedbCrash) ;
  messDumpRegister (acedbLog) ;
  messExitRegister (acedbExit) ;

  /* routine to simplyu clean up read/write locks */
  messcrashroutine = simpleCrashRoutine;
#ifdef CHRONO
  chronoStart () ;
#endif  
     /* The order of the initialisations is crucial */
  chrono ("preInit") ;
  blockInit() ;			/* Initialisation du cache */  
  lexInit() ;                   /* Initialisation du lexique */
  pickPreInit() ;		/* Enough SysClass definitions 
				   to read the session voc */
  sessionDoInit() ;             /* Checks ACEDB.wrm 
				   (and boostraps if neces.)
				   Sets the number thisSession.session
				   sets readlock for thisSession 
				   Reads in the lexiques etc, */

  lastSession = lexLastKey(_VSession);
  if (lastSession == 0)
    messcrash("The session class is empty, Reconstruct, sorry") ;
  sessionStart(lastSession) ;

  chronoReturn () ;
  chrono ("realInit") ;

  pickInit() ;                  /* Class definitions */
  prefInit();                   /* Preferences */
  pickCheckConsistency() ;      /* Check tag values and class hierarchy */
  lexAlphaClear () ;            /* Prevents pointless savings */
  sessionInitialize () ;	/* must come after pickInit */
  qvInit() ;			/* in quovadis.c */
  bIndexInit(BINDEX_AUTO) ;           /* will do only from disk */
  cDNAAlignInit () ;

  sessionKeepAlive = getConfiguredKeepAlive();

  {
    /* purge the readlocks directory of timed-out readlocks */
    Array tmp = readlockGetSessionNumbers();
    arrayDestroy (tmp);
  }

  chronoReturn () ;
  return;
} /* sessionInit2 */

/*************************************************************/
/******* Global way to start stop the acedb kernel ***********/
/*************************************************************/

static void defaultAceCrashRoutine (void)
 /* called by acedbCrash() in case of a messcrash and
    acedbExit() in case of a messExit to clean-up */
 /* this is usually overridden by the user to provide their
    own clean-up routine that may allow the user to save their work */
{ 
  aceQuit (FALSE);		/* do not save, just relase locks */

  fprintf (stderr,"// NB: changes made since the last explicit save "
	              "have NOT been saved\n") ;
} /* defaultAceCrashRoutine */

/*************************************************************/

BOOL aceInit (const char *ace_path)
{
  extern VoidRoutine messcrashroutine;

  sessionInit2 (ace_path) ;

  freesetfile(stdin,"") ; /* detach */

  /* called by acedbCrash (which is called by messcrash)
     this function will perform the aceQuit procedures
     in case of a crash to close the ace-context & database down */
  messcrashroutine = defaultAceCrashRoutine;

  return TRUE ;
}

/*********************************************************************/

BOOL aceQuit(BOOL commit)
     /* complementary call to aceInit() */
     /* in case of a messcrash, the defaultAceCrashRoutine()
	or (whatever messcrashroutine, there is at the moment)
	is called instead by acedbCrash(), to cleanup
	tmp-files, read/write locks etc. */
{
  sessionClose(commit) ;	/* close session,
				   save session if commit==TRUE,
				   release write access if necc.*/
  filtmpcleanup () ;		/* deletes any temporary files made */
  readlockDeleteFile ();	/* release readlock */

  sessionShutDown () ;
  
  return TRUE ;
}

/*********************************************************************/
/*********************************************************************/
/* definitelly release all memories before exit */

void sessionShutDown (void)
{
  int n0, n1 ;
  BOOL debug = FALSE ;

#ifdef CHRONO
  chronoStop () ;
  chronoReport () ;
#endif  
  messAllocStatus (&n1) ;
  
  if (debug) freeOutf ("%18s%8s%8s\n","// Layer","Freed","Used") ;
  if (debug) freeOutf ("%18s%8d%8d kb\n","// Shutdown start", 0, n1 >> 10) ; 
  n0 = n1 ;

  lexAlphaShutDown () ; 
  messAllocStatus (&n1) ; 
  if (debug) freeOutf ("%18s%8d%8d kb\n","// LexAlpha", (n0 - n1) >> 10, n1 >> 10) ;
  n0 = n1 ;

  lexShutDown () ; 
  messAllocStatus (&n1) ; 
  if (debug) freeOutf ("%18s%8d%8d kb\n","// Lex", (n0 - n1) >> 10, n1 >> 10) ;
  n0 = n1 ;

  bIndexShutDown () ;
  messAllocStatus (&n1) ; 
  if (debug) freeOutf ("%18s%8d%8d kb\n","// Bindex", (n0 - n1) >> 10, n1 >> 10) ;
  n0 = n1 ;

  bsMakePathShutDown () ; 
  messAllocStatus (&n1) ; 
  if (debug) freeOutf ("%18s%8d%8d kb\n","// MakePaths", (n0 - n1) >> 10, n1 >> 10) ;
  n0 = n1 ;

  cacheShutDown () ;
  messAllocStatus (&n1) ; 
  if (debug) freeOutf ("%18s%8d%8d kb\n","// Cache2", (n0 - n1) >> 10, n1 >> 10) ;
  n0 = n1 ;

#ifndef ACEDB5
  blockShutDown  () ;
#endif
  messAllocStatus (&n1) ; 
  if (debug) freeOutf ("%18s%8d%8d kb\n","// Blocks", (n0 - n1) >> 10, n1 >> 10) ;
  n0 = n1 ;

  BSshutdown () ;
  messAllocStatus (&n1) ; 
  if (debug) freeOutf ("%18s%8d%8d kb\n","// BS", (n0 - n1) >> 10, n1 >> 10) ;
  n0 = n1 ;

  BTshutdown () ;
  messAllocStatus (&n1) ; 
  if (debug) freeOutf ("%18s%8d%8d kb\n","// BT", (n0 - n1) >> 10, n1 >> 10) ;
  n0 = n1 ;

  OBJShutDown () ;
  messAllocStatus (&n1) ; 
  if (debug) freeOutf ("%18s%8d%8d kb\n","// OBJ", (n0 - n1) >> 10, n1 >> 10) ;
  n0 = n1 ;
  
  diskShutDown () ;
  messAllocStatus (&n1) ; 
  if (debug) freeOutf ("%18s%8d%8d kb\n","// OBJ", (n0 - n1) >> 10, n1 >> 10) ;
  n0 = n1 ;
  
  CDTShutDown () ;
  scShutDown () ;
  condShutDown () ;
  /*
    prefShutDown () ;
  */ 

  messfree (messcrashbuf) ;
  if (debug)  freeOutf ("\n") ;
  freeOutShutDown () ; /* cannot use freeout after freeshutdown, which anyway uses very little */  
  freeshutdown () ;
  messAllocShutdown () ;
  n0 = n1 ;
}

/*************************************************************/

BOOL I_own_the_process(void)
     /* private to session module */
{
#if !(defined(MACINTOSH) || defined(WIN32))
  return getuid() == geteuid() ? TRUE : FALSE  ;
#else
  return TRUE ;
#endif
}

/*************************************************************/
/************** Disk Lock System *****************************/
/*************************************************************/

static BOOL showWriteLock (void)
{
  static char* logName = "Don't-know-who" ;
  int n = 0, level ;
  FILE* lockFile = 0;
    
  lockFile = filopen (sessionFilName ("database/lock", "wrm", 0), 0, "r") ;
  if (!lockFile)
    return FALSE ;
    
  level = freesetfile (lockFile,"") ;
  if (freecard(level))
    { freenext() ;
      if (freeint(&n))
	{ freenext() ;
	  logName = freeword() ;
	}
    }
  freeclose (level) ;
  messout ("%s (session %d) already has write access", logName, n) ;
  return TRUE ;
} /* showWriteLock */

/**************************/

static int lockFd = -1 ;	/* lock file descriptor */

static void releaseWriteLock (void)
{

  if (euid != ruid)
    seteuid(euid);
  if (lockFd == -1)
    messcrash ("Unbalanced call to releaseWriteLock") ;

  if (!suppressLock)
    if (lockf (lockFd, F_ULOCK, 0)) /* non-blocking */
      messcrash ("releaseWriteLock cannot unlock the lock file");

  close (lockFd) ;
  if (!filremove(sessionFilName ("database/lock", "wrm", 0), 0))
    messerror ("Cannot unlock database "
	       "by removing database/lock.wrm (%s)", 
	       messSysErrorText()) ;
  lockFd = -1 ;
  if (euid != ruid)
    seteuid(ruid);
  messdump ("Write access dropped") ; 
} /* releaseWriteLock */

/**************************/

static BOOL setWriteLock (void)
{  
  unsigned int size ;
  char *cp ;

  if (lockFd != -1)
    messcrash ("Unbalanced call to setWriteLock") ;

  if (suppressLock)
    {

      /* check for existence of unlocked lock.wrm file
	 if it exists, then the database is considered locked
      */       
      lockFd = open(sessionFilName ("database/lock", "wrm", 0), O_RDONLY);
      if(lockFd != -1)
	{
	  /* the lock-file exists */
	  showWriteLock();
	  close(lockFd);
	  lockFd = -1;
	  return FALSE;
	}
    }
  if (euid != ruid)
    seteuid(euid);
#if  defined(WIN32) 
  lockFd = open (sessionFilName ("database/lock", "wrm", 0), 
		 O_RDWR | O_CREAT, S_IREAD | S_IWRITE ) ; 
#elif defined(MACINTOSH) || defined (MAC_X)
  /* 2008_08 this is wrong on mac 64, should be 0x0080 #define O_SYNC  0x04  hopefully stable ! including sys/vnode.h creates further problems */
  lockFd = open (sessionFilName ("database/lock", "wrm", 0), 
		 O_RDWR | O_CREAT | O_SYNC, 0644 ) ;
#else
  lockFd = open (sessionFilName ("database/lock", "wrm", 0),
		 O_RDWR | O_CREAT | O_SYNC, 0644) ;
#endif
  if (lockFd == -1)
    { 
      /* is there a permissions problem ? */
      messerror ("I can't open or create the lock file %s. (%s)", 
		 sessionFilName ("database/lock", "wrm", 0), 
		 messSysErrorText()) ;
      if (euid != ruid)
	seteuid(ruid);
      return FALSE ;
    }

  if (!suppressLock)
    if (lockf (lockFd, F_TLOCK, 0)) /* non-blocking */
      { 
	  if (errno == EACCES)
	    {
	      /* "Permission denied" - standard case
	       *  can't lock file - show who has locked locked it */
	      showWriteLock();
	    }
	  else
	    {
	      /* other error */
	      /* the errno ENOLCK - "No locks available"
		 hints to a possible misconfiguration 
		 of the lockdaemon, see man-page lockf(3) */
	      messout ("Sorry, I am unable to establish a lock "
		       "on file database/lock.wrm\n"
		       "(%s)\n"
		       "There may be a problem with the lock daemon "
		       "on this machine, esp. if running NFS.\n"
		       "To avoid the use of the lock daemon, "
		       "insert the word NOLOCK\n"
		       "in the file %s",
		       messSysErrorText(),
		       sessionFilName ("wspec/passwd", "wrm", 0));
	    }

	close(lockFd);
        lockFd = -1 ;
	if (euid != ruid)
	  seteuid(ruid);
	return FALSE ;
      }
	
  cp = messprintf ("%d %s\n\n", thisSession.session, getLogin(TRUE)) ;
  size = strlen(cp) + 1 ;
  write (lockFd, cp, size) ;
  if (euid != ruid)
    seteuid(ruid);
  messdump ("Write access") ; 
  return TRUE;
} /* setWriteLock */

/****************************************************************/
/****************************************************************/

static mytime_t  thisSessionStartTime = 0 ;

void sessionInitialize (void)
     /* private function in session module */
{
  thisSessionStartTime = timeNow() ;
  thisSession.userKey = 0 ; /* so sessionUserKey will reset */
  lexSessionStart () ;
} /* sessionInitialize */

/************************************************************/

void sessionStart (KEY fromSession)
     /* private function in session module */
{
  OBJ Session = 0 ;
  int  a3 = 0, b3 = 0, c3 = 0, h3 = 0,
    parent_index_version = 0 ;
  
  Session = bsCreate(fromSession) ;
  if(!Session )
    { messout("Sorry, I cannot find %s", name(fromSession)) ;
      return ;
    }
  if(!bsGetData(Session,_VocLex,_Int,&a3) ||
     !bsGetData(Session,_bsRight,_Int,&b3)  ||
     !bsGetData(Session,_bsRight,_Int,&c3)  ||
         /* Hasher only introduced in release 1.4  !bsGetData(Session,_bsRight,_Int,&h0)  || */
     !bsGetData(Session,_CodeRelease,
	    _Int,&thisSession.mainCodeRelease) ||
     !bsGetData(Session,_bsRight,
		_Int,&thisSession.subCodeRelease) ||
     !bsGetData(Session,_DataRelease,
	    _Int,&thisSession.mainDataRelease) ||
     !bsGetData(Session,_bsRight,
	    _Int,&thisSession.subDataRelease) ) 
    { messout("No address in this session") ;
      bsDestroy(Session) ;
      return ;
    }
  bsGetData(Session,_IndexVersion, 
	    _Int,&parent_index_version) ;
  bsDestroy(Session) ;

  thisSession.from = fromSession ;
  thisSession.upLink = fromSession ;
  thisSession.index_version = bIndexVersion(parent_index_version) ;

  lexOverLoad(__lexi3,(DISK)a3) ;
  lexOverLoad(__voc3,(DISK)c3) ;
  lexOverLoad(__lexh3,(DISK)h3) ;

  lexRead() ;

  return;
} /* sessionStart */

/*************************************************************/

static BOOL sessionOedipe (void) /* Kill the father of the present session */
{ KEY father = thisSession.upLink, grandFather, fPlus, fMinus ;
  OBJ Father ; 
  Array fPlusArray, fMinusArray ;

   if(!father || !(Father = bsUpdate(father)))
    return FALSE ;
  if(!bsGetKey(Father,_BatPlus,&fPlus) ||
     !bsGetKey(Father,_BatMinus,&fMinus) ) 
    { messout("No address in this father") ;
      bsDestroy(Father) ;
      return FALSE ;
    }

  if(!bsFindTag(Father,_Permanent_session))
    { messout("The previous session being delared permanent, I keep it") ;
      bsDestroy(Father) ;
      return FALSE ;
    }

  if (bsGetKey(Father, _Up_linked_to, &grandFather))
    { thisSession.upLink = grandFather ;
      bsAddData(Father, _Destroyed_by_session, _Int, &thisSession.session) ;
      bsSave(Father) ;
    }
  else
   { thisSession.upLink = 0 ;
     bsDestroy(Father) ;
   }  
  
  fPlusArray = arrayGet(fPlus,unsigned char,"c") ;
  fMinusArray = arrayGet(fMinus,unsigned char,"c") ;

  diskFuseOedipeBats(fPlusArray, fMinusArray) ; 

  arrayDestroy(fPlusArray) ;
  arrayDestroy(fMinusArray) ;
              /* Kill the father bats */
  arrayKill(fPlus) ;
  arrayKill(fMinus) ;

  saveAll() ;
  return TRUE ;
} /* sessionOedipe */

/******************/

static BOOL sessionFuse (KEY father, KEY son) 
{ 
  KEY fPlus, fMinus, sPlus, sMinus, grandFather ;
  OBJ Father = 0 , Son = 0 ; 
  int nPlus = 0, nMinus = 0 ; 

  Father = bsUpdate(father) ; 
  Son = bsUpdate(son) ;
  if(!Father || !Son || 
     !bsGetKey(Father,_BatPlus,&fPlus) ||
     !bsGetKey(Father,_BatMinus,&fMinus) ||
     bsFindTag(Father, _Destroyed_by_session) ||
     !bsGetKey(Son,_BatPlus,&sPlus) ||
     !bsGetKey(Son,_BatMinus,&sMinus) 
     )
    goto abort ;

  if (!diskFuseBats(fPlus, fMinus, sPlus, sMinus, &nPlus, &nMinus)) 
    goto abort ;

  /* Store the sons bat sizes */
  bsAddKey(Son, _BatPlus, sPlus) ;   /* to reposition */
  bsAddData(Son, _bsRight, _Int, &nPlus) ;
  bsAddKey(Son, _BatMinus, sMinus) ;
  bsAddData(Son, _bsRight, _Int, &nMinus) ;
  
  bsAddData(Father, _Destroyed_by_session, _Int, &thisSession.session) ;
  if (bsGetKey(Father, _Up_linked_to,&grandFather))
    bsAddKey(Son , _Up_linked_to, grandFather) ;

  bsSave(Father) ;
  bsSave(Son) ;
  saveAll() ;  
  return TRUE ;

 abort:
  bsSave(Son) ;
  bsSave(Father) ;
  return FALSE ;
} /* sessionFuse */

/************************************************************/

Array sessionTreeCreate(BOOL showDeads)
     /* private function in session module */
{ 
  KEY key = 0 , kk ; 
  int treeArrayMax ; 
  char *cp ;
  OBJ obj ;
  ST *st , *st1 ;
  int gen, i, j;
  BOOL found ;
  Array sessionTree;

  if (0) /* if we loose control of the Xterm, we cannot save */ 
    messStatus ("Reading sessions");

  sessionTree = arrayCreate(50, ST) ;
  
  /* Add all Session objects to the session tree */
  treeArrayMax = 0;
  while (lexNext(_VSession, &key))
    if ((obj = bsCreate(key)))
      {
	st = arrayp(sessionTree, treeArrayMax++, ST) ;
	st->key = key ;
	bsGetData(obj,_Session, _Int, &st->number);

	if (bsGetKey(obj,_Created_from, &kk))
	  st->father = st->ancester = kk ;
	if (bsGetKey(obj,_Up_linked_to, &kk))
	  st->ancester = kk ;
	else if (st->key != KEYMAKE(_VSession, 2))
	  st->isDestroyed = TRUE ; /* i.e. to old a release */
	if (bsFindTag(obj,_Destroyed_by_session))
	  st->isDestroyed = TRUE ;
	if (bsFindTag(obj,_Permanent_session))
	  st->isPermanent = TRUE ;

	if (bsGetData(obj,_User, _Text, &cp))
	  st->user = strnew(cp, 0);

	if (bsGetData(obj,_Date, _Text, &cp))
	  st->date = strnew(cp, 0);

	if (bsGetData(obj,_Session_Title, _Text, &cp))
	  st->title = strnew(cp, 0);
	
	bsDestroy(obj) ;
      }

  /* Add thisSession to the session tree 
     (the object for the current session is not a finalized object
     yet and therfore won't show up in the lexer up) */
  st = arrayp(sessionTree, treeArrayMax++, ST) ;
  st->key = _This_session ;
  st->father = st->ancester = thisSession.from ;

  /* links the st */
  for (i = 0, st = arrp(sessionTree,0,ST) ; i < arrayMax(sessionTree) ; st++, i++)
    for (j = 0, st1 = arrp(sessionTree,0,ST) ; j < arrayMax(sessionTree) ; st1++, j++)
      if(st->key == st1->ancester)
	st1->sta = st ;

  /* Count generations */
  if (showDeads)
    { 
      for (i = 0, st = arrp(sessionTree,0,ST); i < arrayMax(sessionTree) ; st++, i++)
	if (!st->ancester)
	  st->generation = 1 ;
    }
  else
    {
      for (i = 0, st = arrp(sessionTree,0,ST) ; i < arrayMax(sessionTree) ; st++, i++)
	if ((!st->isDestroyed) && ( !st->ancester || (st->sta && st->sta->isDestroyed)) )
	  st->generation = 1 ;
    }
  gen = 1 ;
  found = TRUE ;
  while (found)
    { found = FALSE ;
      if (showDeads)
	{
	  for (i = 0, st = arrp(sessionTree,0,ST) ; i < arrayMax(sessionTree) ; st++, i++)
	    if (st->sta && st->sta->generation == gen)
	      { st->generation = gen + 1 ;
		found = TRUE ;
	      }
	}
      else
	{
	  for (i = 0, st = arrp(sessionTree,0,ST) ; i < arrayMax(sessionTree) ; st++, i++)
	    if (!st->isDestroyed && st->sta && st->sta->generation == gen)
	      { st->generation = gen + 1 ;
		found = TRUE ;
	      }
	}
      gen++ ;
    }

  /* mark readlocked sessions, they won't get destroyed, because
     other processes still have that data cached */
  {
    Array readlockedSessions;
    int n, locked_num;

    if ((readlockedSessions = readlockGetSessionNumbers()))
      {
	for (n = 0; n < arrayMax(readlockedSessions); n++)
	  {
	    locked_num = arr (readlockedSessions, n, int);

	    for (i = 0, st = arrp(sessionTree,0,ST) ;
		 i < arrayMax(sessionTree) ; st++, i++)
	      if (st->number == locked_num - 1)
		st->isReadlocked = TRUE;
	  }
	  arrayDestroy (readlockedSessions);
      }
  }

  return sessionTree;
} /* sessionTreeConstruct */
    
/************************************************************/  

void sessionTreeDestroy(Array sessionTree)
     /* private function in session module */
{
  ST *st ;
  int i;
 
  for (i = 0, st = arrp(sessionTree, 0, ST);
       i < arrayMax(sessionTree); st++, i++)
    {
      if (st->user) messfree (st->user);
      if (st->date) messfree (st->date);
      if (st->title) messfree (st->title);
    }
  
  arrayDestroy (sessionTree);

  return;
} /* sessionTreeDestroy */

/************************************************************/  

static int sessionFindSon(Array sessionTree, KEY father, KEY *sonp)
{ 
  int i, n = 0;
  ST *st, *sta ;

  *sonp = 0 ;             /* Suicide */
  if (father == _This_session)
    return -2 ;
  if (!father)           /* Trivial */
    return -3 ;
  if (father <= KEYMAKE(_VSession, 2))
    return -6 ;

  for (i = 0, sta = arrp(sessionTree,0,ST) ; i < arrayMax(sessionTree) ; sta++, i++)
    if (sta->key == father)
      break ;

  if (sta->key != father)       /* not found error */
    return -5 ;
  
  if (sta->isDestroyed)          /* dead */
      return -3 ;
 
  for (i = 0, st = arrp(sessionTree,0,ST) ; i < arrayMax(sessionTree) ; st++, i++)
    if(!st->isDestroyed && st->sta == sta)
      { *sonp = st->key ;
	n++ ;
      }
  
  if (n == 1)            /* is st1 in my direct ancestors */
    { for (i = 0, st = arrp(sessionTree,0,ST) ; i < arrayMax(sessionTree) ; st++, i++)
	if (st->key == thisSession.from)
	  {
	    while (st)
	      { if (st == sta)
		  return 1 ;
		st = st->sta ;  
	      }
	    return -4 ;
	  }
      return -4 ;        /* not ancestor */
    }

  return n ;             /* number of sons */
} /* sessionFindSon */

/************************************************************/

static BOOL sessionKillAncestors (Array sessionTree, int n)
/* Kill all father over level n */
{
  KEY father = thisSession.upLink, son ;
  ST *st ;
  int i ;

  for (i = 0, st = arrp(sessionTree,0,ST);
       i < arrayMax(sessionTree); st++, i++)
    if (st->key == father)
      break ;

  if (father && st->key != father) 
    { 
      messerror ("sessionKillAncestors() - "
		 "could not find the father he wanted to murder.") ;
      return FALSE ;
    }
 
  while (st && --n)            /* Move up the ladder n times */
    st = st->sta ;
  if (!st || n)                /* Not many ancesters alive */
    return FALSE ;
  
  while (st && (st->isPermanent || 
		st->isReadlocked || 
		sessionFindSon(sessionTree, st->key, &son) != 1))
    st = st->sta ;            /* Try higher up */
    
  if (st && st->key && son)
    return sessionFuse(st->key, son) ;

  return FALSE ;
} /* sessionKillAncestors */

/*************************************************************/

static char* getConfiguredDatabaseName (void)
{ 
  static char dbname[32] = "";
  FILE *fil ;
  int level ;
  char *word ;

  if (sessionFilName("wspec/database", "wrm", "r") &&
      (fil = filopen (sessionFilName("wspec/database", "wrm", 0),0,"r")))
    {
      level = freesetfile (fil, 0) ;
      while (freecard (level))
	if ((word = freeword()) &&
	    !strcmp (word, "NAME") && (word = freeword()))
	  {
	    if (strlen (word) > 31)
	      messExit ("Error in wspec/database.wrm : "
			"the NAME is too long (31 characters max).") ;
	    else
	      strcpy (dbname, word) ;
	  }
      /* NOTE, that we will loop over the entire file, so the
	 free-package closes the file when it reaches the end */
    }

  return &dbname[0];
} /* getConfiguredDatabaseName */

/*************************************************************/

static int getConfiguredKeepAlive (void)
{ 
  FILE *fil ;
  int level ;
  char *word ;
  int keepAlive = 0;

  if (sessionFilName("wspec/database", "wrm", "r") &&
      (fil = filopen (sessionFilName("wspec/database", "wrm", 0),0,"r")))
    {
      level = freesetfile (fil, 0) ;
      while (freecard (level))
	if ((word = freeword()) &&
	    !strcmp (word, "SESSION_KEEP_ALIVE") && (freeint(&keepAlive)))
	  {
	    if (keepAlive < 2)
	      messExit ("Error in wspec/database.wrm : "
			"smallest legal value for SESSION_KEEP_ALIVE is 2") ;
	  }
      /* NOTE, that we will loop over the entire file, so the
	 free-package closes the file when it reaches the end */
    }
  if (keepAlive == 0)
    keepAlive = 2;	   /* default, if not defined in database.wrm */

  return keepAlive;
} /* getConfiguredKeepAlive */

/*************************************************************/

static void sessionDoInit(void)
{
  FILE *f ;

  if (sessionFilName ("database/ACEDB","wrm","r")) 
    { 
      /* this doesn't swap the superblock (write does) */
#ifndef ACEDB5
      dataBaseAccessInit();
      diskblockread (&theSuperBlock,1) ; 
#ifdef ACEDB4
      if (theSuperBlock.byteSexIndicator != 1)
	{ swapData = TRUE;
	  swapSuperBlock(&theSuperBlock);
	  if (!getenv("ACEDB_NO_BANNER"))
	    printf("// Database has non-native byte-order: swapping\n");
	}
      else
	{ 
	  swapData = FALSE;
	}
#endif /* ACEDB4 */
#else
      aDiskAccessInit(&theSuperBlock);
#endif /* ACEDB5 */
      if (theSuperBlock.mainRelease != aceGetVersion())
	{
	  if (getenv ("ACEDB_CHANGE_MAIN_RELEASE")) /* override */
	    theSuperBlock.mainRelease = aceGetVersion() ;
	  else
	    messExit
	      ("You are running %s of the code.\n"
	       "The database was created under release %d.\n"
	       "Either use the correct binary or upgrade the database.\n",
	       aceGetVersionString(), theSuperBlock.mainRelease) ;
	}
#ifdef ACEDB5
      thisSession.session = ++theSuperBlock.session ;
#else
      thisSession.session = ++theSuperBlock.h.session ;
#endif
      thisSession.gAddress = theSuperBlock.gAddress ;
      lexOverLoad (__lexi1,theSuperBlock.gAddress) ;
      lexReadGlobalTables () ;

      if (!readlockCreateFile())
	fprintf (stderr,
		 "Unable to set readlock for session %d.\n",
		 thisSession.session);

      return ;
    }
  /* else reinitialise the system */
  
  /* Only the acedb administrator is entitled to proceed
   * I check this by (a) getuid == geteuid
   *                 (b) trying to open passwd.wrm for writing
   * In the install documentation, it is stated that
   * passwd.wrm should be set chmod 644
   */
  
  if (!I_own_the_process () || 
      /* check whether we could write to the password file, but
	 only if it exists */
      !(sessionFilName ("wspec/passwd", "wrm", "r") &&
	sessionFilName ("wspec/passwd", "wrm", "a")))
    messExit
      ("Sorry, the database seems empty.  Only someone with "
       "write access to the password file can reinitialise it.");
  
  if (!messQuery ("The file %s does not exist, \n"
		  "indicating that the database is empty.\n"
		  "  Should I re-initialise the system?",
		  sessionFilName ("database/ACEDB", "wrm", 0)))
    messExit ("Database will remain uninitialized. Bye!") ;

#if !defined(MACINTOSH)
  if (!sessionCheckUserName (TRUE))  /* Parses the eventual NOLOCK,
					Wolf May 1994 */
    messExit ("Not authorized to initialize the database. Check wspec/passwd.wrm");
#endif

  if (!setWriteLock())
    /* this can happen if another process is also in the middle
       of the initialization process at the same time */
    messExit ("Cannot set write lock! Initialization failed!") ;

  /****************************************************************/
  /****** all permissions clear, boostrap an empty database *******/
  /****************************************************************/

  thisSession.session = 0;	/* special session number to be used
				   for this first readlock 
				   will be updated by sessionDoClose()
				   called at the end of this function */
  if (!readlockCreateFile())
    messExit ("Unable to set readlock in database/readlocks/");
  
  /****************************************************************/

  writeAccess = TRUE  ; /* implied */
  swapData = FALSE ; /* The database file must be in our order */
  thisSession.session = 1 ;
  strcpy (thisSession.name, "Empty_db") ;
  thisSession.gAddress = 0 ; /* so far */
  thisSession.mainCodeRelease = aceGetVersion() ; 
  thisSession.subCodeRelease = aceGetRelease() ; 
  thisSession.mainDataRelease = aceGetVersion() ; 
  thisSession.subDataRelease = 0 ; 
  thisSession.from = 0 ; 
  thisSession.upLink = 0 ; 

#ifndef ACEDB5
  theSuperBlock.h.disk = 1 ;
  theSuperBlock.h.nextdisk = 0 ;
  theSuperBlock.h.type = 'S' ;
  theSuperBlock.h.key  = 0 ;
  theSuperBlock.byteSexIndicator = 1;
#endif

  strcpy (theSuperBlock.dbName, getConfiguredDatabaseName());

  /* DataBase creation */
  if (!dataBaseCreate())
    messExit ("Database Creation failed, check write permissions etc.") ;


  messdump("\n**********************************************************\n");
  messdump("Making a new database file and starting at session 1 again");
  
  readModels ();

  /* application specific initialisations, in quovadis.wrm */
  diskPrepare ();

  sessionInitialize ();		/* do this last before saving */

  if (!(f = filopen (sessionFilName("database/ACEDB", "wrm", 0), 0, "w")))
    messExit ("Sorry, the operating system won't let me create %s",
	       sessionFilName("database/ACEDB", "wrm", 0)) ;
  fprintf (f, "# DESTROY this file to reinitialize the database\n") ;
  fprintf (f, "# BUT DO NOT EDIT THIS FILE\n") ;
  fprintf (f, "# It was created automatically by the acedb kernel on %s\n", timeShowNow ()) ;
  fprintf (f, "ACEDB_DISK_BLOCK_SIZE=%d\n",BLOC_SIZE/1024) ;
 
  filclose (f) ;
  
  f = filopen (sessionFilName("database/lock","wrm",0), 0, "w") ;
  fprintf (f, "The lock file\n") ;
  filclose (f) ;

#ifndef NEW_MODELS
  readModels () ; /* second read needed to register the _VModel obj */
#endif
  bIndexInit(BINDEX_AUTO) ;
  sessionDoClose () ;	/* real save without sessionInitialize
			   will also update readlock-file */
} /* sessionDoInit */

/************************************************************/

char *sessionDbName (void)
     /* general public function */
{
  /* returns the name of the database as gathered from the superblock,
     used for the window title and messages in data updates */
 if (theSuperBlock.mainRelease < 4 || theSuperBlock.subCodeRelease < 3)
    return 0 ;			/* not implemented before 4.3 */
  else if (*theSuperBlock.dbName)
    return theSuperBlock.dbName ;
  else
    return 0 ;
} /* sessionDbName */

/************************************************************/

static void sessionSaveSessionObject(void)
{
  /* store the Session Object */
  /* update the BAT block access table for thisSession */
  DISK d ;
  OBJ Session ;
  KEY kPlus, kMinus ;
  unsigned long int free2, plus, minus, used, ndread, ndwrite ;
  char timeBuf[25] ;


  Session = bsUpdate (thisSession.key) ;
  if(!Session)
    messcrash
      ("Session register fails to create the session object") ;

  if (thisSession.session == 1)
    bsAddTag(Session,_Permanent_session) ;

  bsAddData(Session,_User,_Text,thisSession.name) ;

  if (!bsFindTag (Session, _Date)) /* don't add date twice - it changes => save loop */
    {
      bsAddTag (Session, _Date) ;

      if (bsType (Session, _bsRight) == _Text)
	/* dates used to be _Text, we have to account for that */
	bsAddData (Session, _Date, _Text, timeShow(timeNow(), timeBuf, 25)) ;
      else
	/* with recent versions, the date is _DateType */
	{ 
	  mytime_t session_date = timeNow() ;
	  bsAddData (Session, _Date, _DateType, &session_date) ;
	}
    }

  if (*thisSession.title)
    bsAddData(Session,_Session_Title,_Text,thisSession.title) ;
  bsAddData(Session,_Session,_Int,&thisSession.session) ;
  if(thisSession.from)
    bsAddKey(Session,_Created_from,thisSession.from) ;
  if(thisSession.upLink)
      bsAddKey(Session,_Up_linked_to,thisSession.upLink) ;

  thisSession.mainCodeRelease = aceGetVersion() ; 
  thisSession.subCodeRelease = aceGetRelease() ; 
  bsAddData(Session,_CodeRelease,
	    _Int,&thisSession.mainCodeRelease) ;
  bsAddData(Session,_bsRight,
	    _Int,&thisSession.subCodeRelease) ;
  bsAddData(Session,_DataRelease,
	    _Int,&thisSession.mainDataRelease) ;
  bsAddData(Session,_bsRight,
	    _Int,&thisSession.subDataRelease) ;
  { int v = bIndexVersion(-1) ;
    if (v)
      bsAddData(Session,_IndexVersion,
		_Int,&v) ;
  }

  diskavail(&free2,&used,&plus,&minus,&ndread,&ndwrite) ;

  d = lexDisk(__Global_Bat) ;
  bsAddData(Session,_GlobalBat,_Int,&d) ; 
  bsAddData(Session,_bsRight,_Int,&used) ;

  d = lexDisk(__lexi1) ;
  bsAddData(Session,_GlobalLex,_Int,&d) ;

  d = lexDisk(__lexi2) ;
  bsAddData(Session,_SessionLex,_Int,&d) ;
  d = lexDisk(__lexa2) ;
  bsAddData(Session,_bsRight,_Int,&d) ;
  d = lexDisk(__voc2) ;
  bsAddData(Session,_bsRight,_Int,&d) ;
  d = lexDisk(__lexh2) ;
  bsAddData(Session,_bsRight,_Int,&d) ;

  d = lexDisk(__lexi3) ;
   bsAddData(Session,_VocLex,_Int,&d) ;
  d = lexDisk(__lexa3) ;
  bsAddData(Session,_bsRight,_Int,&d) ;
  d = lexDisk(__voc3) ;
  bsAddData(Session,_bsRight,_Int,&d) ;
  d = lexDisk(__lexh3) ;
  bsAddData(Session,_bsRight,_Int,&d) ;

  lexaddkey(messprintf("p-%d",thisSession.session),&kPlus,_VBat) ;
  lexaddkey(messprintf("m-%d",thisSession.session),&kMinus,_VBat) ;
 
  {
    int plus2 = (int) plus ;
    int minus2 = (int) plus ;

    bsAddKey(Session,_BatPlus,kPlus) ;
    bsAddData(Session,_bsRight,_Int,&plus2) ;
    
    bsAddKey(Session,_BatMinus,kMinus) ;
    bsAddData(Session,_bsRight,_Int,&minus2) ;
  }
  bsSave(Session) ;
} /* sessionSaveSessionObject */

/************************************************************/

static void sessionSaveUserSessionObject (void)
{ 
  /* complete the UserSession object */
  /* this might be called several times over 
     in the process of saveAll() */
  OBJ obj ; 
  KEY key, tag ;
  mytime_t timEnd = timeNow() ;
  
  lexSessionEnd () ;	/* writes new and touched keysets */
  
  sessionUserKey () ;
  obj = bsUpdate (thisSession.userKey) ;
  
  bsAddKey (obj, _Session, thisSession.key) ;
  bsAddData (obj, _User, _Text, thisSession.name) ;
  
  /* don't change dates again, once set, to avoid causing
     saveAll() to think that modifications have happened
     in the multiple saving loop */
  if (!bsFindTag (obj, _Start))
    {
      if (!thisSessionStartTime)
	messcrash ("Can't save UserSession object "
		   "with unititialized StartTime");
      bsAddData (obj, _Start, _DateType, &thisSessionStartTime) ;
    }
  if (!bsFindTag (obj, _Finish))
    bsAddData (obj, _Finish, _DateType, &timEnd) ;
  
  if (lexword2key (messprintf ("new-%s", name(thisSession.userKey)),
		   &key, _VKeySet))
    { 
      lexaddkey ("New", &tag, 0) ;
      bsAddKey (obj, tag, key) ;
    }

  if (lexword2key (messprintf ("touched-%s", name(thisSession.userKey)),
		   &key, _VKeySet))
    { 
      lexaddkey ("Touched", &tag, 0) ;
      bsAddKey (obj, tag, key) ;
    }
  
  bsSave (obj) ;
} /* sessionSaveUserSessionObject */

/************************************************************/

void sessionRegister(void)
     /* general public function */
{
  /* this is called several times over during saveAll() */
  OBJ ModelObj ;

  if (!isWriteAccess())
    return ;

#ifndef NEW_MODELS
  {
    OBJ obj1 = bsCreate (KEYMAKE (_VUserSession,0)) ;
    if (!obj1)
      {
	messout ("You are missing the model for UserSession.  "
	       "To bring your database up to date please Read Models.  "
	       "You will not be able to save until you do this.") ;
        return ;
      }
    bsDestroy (obj1) ;
  }
#endif


  strncpy (thisSession.name, getLogin(TRUE), 78) ;
  lexaddkey (messprintf ("s-%d-%s", 
			 thisSession.session, thisSession.name),
	     &thisSession.key, _VSession) ;

#ifdef NEW_MODELS
  { KEY ModelKey ;
    lexaddkey ("#Session", &ModelKey, _VModel) ;
    ModelObj = bsCreate(ModelKey) ;
  }
#else
  ModelObj = bsCreate(KEYMAKE(_VSession,0)) ;
#endif

  if(!ModelObj || !bsFindTag(ModelObj,_Up_linked_to))
    /* Always check on latest tag added to ?Session */
    { messerror("%s\n%s\n%s\n%s\n%s\n",
		"FATAL, an eror occured in sessionRegister",
		"Probably, your model for the session class is not correct",
		"and corresponds to an earlier release.",
		"Look for the latest version of the wspec/models.wrm file,",
		"read the model file and try again.") ;
      exit(1) ;
    }
  bsDestroy(ModelObj) ;


  sessionSaveSessionObject() ;

  sessionSaveUserSessionObject() ;

  /* will be reset for next new session (in sessionInitialize) */
  thisSessionStartTime = 0 ;
} /* sessionRegister */

/*************************************************************/

BOOL sessionDoClose (void)
     /* private function in session module */
{
  BOOL savedSomething = FALSE ;
  KEY son ;
  char timeBuf[25] ;

  if (!isWriteAccess())
    return FALSE ;

  /* Try to save all caches - return TRUE, if there
     had been modifications that were saved to disk*/
  if (saveAll())		
    {
      /* Yes we saved modifications - adjust the session structure*/

      /* I keep only a certain number of session alive 
	 (the sessionKeepAlive static int) */
      if (thisSession.session != 1)
	{
	  Array sessionTree;

	  sessionTree = sessionTreeCreate (FALSE); /* exclude dead sessions */

	  if (sessionKeepAlive != 1)
	    while (sessionKillAncestors(sessionTree, sessionKeepAlive))
	      {
		sessionTreeDestroy(sessionTree);
		sessionTree = sessionTreeCreate (FALSE);
	      }
	  else if (sessionFindSon(sessionTree, thisSession.upLink, &son) == 1)
	    sessionOedipe () ;

	  sessionTreeDestroy (sessionTree);
	}

                  /****************************/
                  /* increment session number */
                  /****************************/

      /* if you edit this code, edit also readsuperblock 
	 and sessionDestroy */
      theSuperBlock.gAddress = thisSession.gAddress = lexDisk(__lexi1) ;
      theSuperBlock.mainRelease = aceGetVersion() ;
      theSuperBlock.subCodeRelease = aceGetRelease() ;
      theSuperBlock.magic = saveMagic ;
      thisSession.from =  thisSession.key ;
      thisSession.upLink =  thisSession.key ;
#ifdef ACEDB5   
      theSuperBlock.session =  thisSession.session ;
#else
      theSuperBlock.h.session =  thisSession.session ;
#endif
      diskWriteSuperBlock(&theSuperBlock); /* writes and zeros the BATS*/

      messdump ("Save session  Session:%d Time:%s", 
		thisSession.session, timeShow(timeNow(), timeBuf, 25)) ;

      readlockUpdateFile (thisSession.session+1); /* update readlock 
						     for next session */

      thisSession.session++;	/* increase current session number */

      savedSomething = TRUE ;

      if (sessionChangeRoutine)
	(*sessionChangeRoutine)();
    }
  /* else savedSomething is still FALSE */

  sessionReleaseWriteAccess()  ;

  if (externalCommit) externalCommit() ;


  return savedSomething ;
} /* sessionDoClose */

/*************************************************************/

BOOL isWriteAccess(void)
     /* general public function */
{ return writeAccess ? TRUE : FALSE ;
} /* isWriteAccess */

/*************************************************************/

BOOL sessionGainWriteAccess (void)
     /* general public function */     
{
  /* gains write access for the user, if possible */

  if (writeAccess)		/* already has write access */
    return TRUE;

  /* we don't have write access, but can we grap it */

  if (dropWriteAccess ||        /* Cannot re-gain write access, once
				   permanently dropped.
				   Cannot gain write access in 
				   an older session */
      !sessionDatabaseWritable(TRUE) || /* will warn if impossible */
      !sessionCheckUserName (TRUE) || /* parses wspec/passwd.wrm
					 gives warning if it fails */
      !setWriteLock() ||	/* will warn if it fails */
      !checkSessionNumber (1))	/* NB after setWriteLock() */
    {
      if (debug)  /* for debugging */
	{
	  messdump ("gainWriteAccess Failed") ;
	  if (dropWriteAccess) messdump ("dropWriteAcces=TRUE ") ;
	  if (!sessionDatabaseWritable(TRUE)) messdump ("sessionDatabaseWritable=FALSE") ; 
	  if (!sessionCheckUserName(TRUE)) messdump ("sessionCheckUserName=FALSE") ; 
	  if (!setWriteLock()) messdump ("setWriteLock=FALSE") ;
	}
      return FALSE;
    }

  /* all checks passed */
  writeAccess = TRUE ;		/* was false before */

  /* therefor calls the routine that gets notified of the change,
     this is usually to re-draw buttons or change the menus, etc */
  if (writeAccessChangeRoutine)
    (*writeAccessChangeRoutine)();

  return TRUE;
} /* sessionGainWriteAccess */

/*************************************************************/

void sessionReleaseWriteAccess (void)
     /* general public function */     
{
  if (isWriteAccess())
    {
      releaseWriteLock () ;
      writeAccess = FALSE ;

      if (writeAccessChangeRoutine)
	(*writeAccessChangeRoutine)();
    }

  return;
} /* sessionReleaseWriteAccess */

/*************************************************************/

static BOOL sessionCheckUserName (BOOL doTalk)
{
#if defined(MACINTOSH)

    extern BOOL passPrompt ( void ) ;

    if(!passPrompt()) {
    	messout ("Sorry, your name was not found in the password file") ;
    	return FALSE ;
    }
    else
      return TRUE ;

#else  /* !MACINTOSH */
  int level ;
  FILE *fil ; 
  register char *cp, *name ;
  Stack uStack ;
  static int result = 0 ;

  /* check only once, since the userid of the process cannot change */
  switch (result)
    {
    case 0: result = 1 ; break ;
    case 1: return FALSE ;
    case 2: return TRUE ;
    }

/* Only users listed in wspec/passwd.wrm can get write access */

/* I cannot use cuserid, because I want the real userid and, 
   depending on systems, cuserid may return either the real or 
   the effective userid which is probably reset by setuid().
   
   On the other hand, Unix getlogin() is replaced by local 
   getLogin() because getlogin() returns NULL if the process 
   is not attached to a terminal */

  name = getLogin(TRUE) ;
  if (debug) messdump ("getLogin=%s",name?name:"name=0") ;
    
  /* use fopen to avoid warning from filopen() */
  fil = fopen (sessionFilName("wspec/passwd","wrm",0), "r") ;
  if (debug) messdump ("wspec=%s",sessionFilName("wspec/passwd","wrm",0)) ;
  if (!fil)
    {				
      messout ("The database has no file wspec/passwd.wrm. To enable "
	       "write access to the database, the usernames "
	       "of trusted administrators have to be listed in the "
	       "file : %s",
	       sessionFilName("wspec/passwd","wrm",0));
      return FALSE ;
    }
  fclose (fil);

  /* now use filopen, because the free-package will use filclose
     to close the file, once it read up until the end */
  fil = filopen (sessionFilName("wspec/passwd","wrm",0), 0, "r");

  uStack = stackCreate(32);

  level = freesetfile (fil, "") ;
  while (freecard (level))
    { freenext () ;
      if (!freestep ('#'))	/* skip # lines */
	if ((cp = freeword()))	/* read user name */
	  pushText(uStack,cp) ;
    }

  stackCursor(uStack, 0) ;
  while (!stackAtEnd (uStack))
    if (!strcmp ("NOLOCK", stackNextText(uStack)))
      suppressLock = TRUE ;

  stackCursor(uStack, 0) ;
  while (!stackAtEnd (uStack))
    if (!strcmp (name, stackNextText(uStack)))
      { stackDestroy (uStack) ;
        result = 2 ;
	return TRUE ;
      }
  stackDestroy(uStack) ;

  /* failed to find the name */
  if (doTalk)
    {
      if (strcmp (name, "getLoginFailed") == 0)
	messout ("On your particular machine, the operating system call "
		 "\"getlogin\" failed. To circumvent this problem, you "
		 "may uncomment the user name \"getLoginFailed\" in "
		 "the file wspec/passwd.wrm.\n"
		 "But from then on, everybody will have write access "
		 "unless you set up other security measures using "
		 "write permissions.") ;
      else
	messout ("%s is not a registered user.  To register a user "
		 "for write access add their login name "
		 "to wspec/passwd.wrm", name) ;
    }

  return FALSE ;
#endif /* !MACINTOSH */
}

/*************************************************************/

static BOOL sessionDatabaseWritable (BOOL doTalk)
{
  FILE *f;
  char *name = sessionFilName("database/test", 0, 0);

  if (euid != ruid)
    seteuid(euid);
  unlink(name);
  f = fopen (name, "a");
  if (f) fclose(f);
  unlink(name);
  if (euid != ruid)
    seteuid(ruid);
  if (!f)
    { 
#if defined(MACINTOSH) 
      if (doTalk)
	messout ("Please create a folder called 'database' and start over");
#else   
      if (doTalk)
	messout ("Sorry, the program does not have permission to write "
		 "files in %s", sessionFilName ("database", 0, 0)) ;
#endif
      return FALSE ;
    }
  return TRUE ;
} /* sessionDatabaseWritable */

/*************************************************************/

  /* If another process has written since I started
     either he has grabed the lock, and I can't write,
     or he may have released it but then he may have updated
     the disk.
     In this case, the session number on block 1 will be
     modified and I will know that my BAT is out of date
     */

static BOOL checkSessionNumber(int type)
{
  int mm, nn ;
#ifndef ACEDB5
  BLOCK bb ;
#endif

/* The call to diskFlushLocalCache forces the latest version of the  */
/* superblock to be read from the NFS server (for remote disks). This closes */
/* a concurrency  hole which can see this function succeed because it is */
/* checking a stale, client side cached copy of the superblock. */

  diskFlushLocalCache();

#ifdef ACEDB5
  nn = aDiskGetGlobalAdress() ; 
  return TRUE ; /* cheat !*/
#else
  diskblockread (&bb,1) ; 
#ifdef ACEDB4
  if (swapData)
    swapSuperBlock(&bb);
#endif
  nn = bb.gAddress ;
  mm = bb.magic ;
#endif
  if (type == 1 && thisSession.gAddress != nn)
    { messout ("Sorry, while you where reading the database, another process "
	       "updated it.\n//  You can go on reading, but to write, you must "
	       "quit and rerun.") ;
      releaseWriteLock() ;
      return FALSE ;
    }
  if (type == 2 && (saveMagic != mm || thisSession.gAddress != nn))
    { messcrash ("Sorry, I cannot save, while you where writing, another process "
	       "updated the database.\n"
	       "// Probably, on your machine locking is faulty or delayed\n"
	       "// If you are working on a single machine, please report this bug"
	       ) ;
    releaseWriteLock() ;
    return FALSE ;
    }
  return TRUE ;
} /* checkSessionNumber */

/*************************************************************/

void sessionForbidWriteAccess (void)
     /* general public function */
{
  /* definitive, can't re-gain write access ever again,
     Should we have a AllowWriteAccess function to match it ?? */

  readlockDeleteFile();		/* read locks no longer needed */
  dropWriteAccess = TRUE ;	/* make it permanent */

  /* the permanent loss of write access will change the
     behaviour of writeAccessPosibble(), and therefor may 
     change buttons or menus */
  if (writeAccessChangeRoutine)
    (*writeAccessChangeRoutine)();
 
  return;
} /* sessionForbidWriteAccess */

/*************************************************************/

BOOL writeAccessPossible(void)
     /* general public function */
{ 
  /* return TRUE if it would be possible to gain write access */
  static int status = 0;
  
  if (dropWriteAccess)
    return FALSE ;
  if (status)
    return status == 1 ? TRUE : FALSE ;
  
  if (sessionDatabaseWritable(FALSE) && sessionCheckUserName(FALSE))
    { status = 1;
      return TRUE;
    }
  else
    { status = 2;
      return FALSE;
    }
} /* writeAccessPossible */

/**************************************************************/

void sessionClose (BOOL doSave)
     /* general public function */
{ 
  sessionAutoSave (-3, 0) ;    /* synchronizes autosave, without activating it */
  if (doSave == FALSE)
    { 
      sessionReleaseWriteAccess() ;
    }
  else
    {
      if (sessionDoClose ()) /* returns TRUE if something was saved */
	{
	  sessionInitialize () ;	/* only init if we had saved */
	}
    }
  sessionAutoSave (-4, 0) ;     /* synchronizes autosave, without activating it */
} /* sessionClose */

/*****************************************************************/
/*****************************************************************/
/******* AutoSaving package, mieg: march 2000 ********************/
/*****************************************************************/
/* Saves if more than timeLapse seconds elapsed since last save 
 * or since timeZero if set
 * returns time since last save
 */

unsigned int sessionTimeSave (unsigned int timeLapse, int mode)
{ 
  /* use time(0) when possible (much faster) and timeNow() when saving */
  int dt = 0 ;
  time_t now = time (0) ;  /* in seconds */
  BOOL debug = FALSE ;
  static time_t timeZero = 0 ;

  switch (mode)
    {
    case 2:
      timeZero = now ; 
      break ;
    
    case 1:
      timeZero = 0 ;
      break ;
    
    case 0:
      if (timeZero)
	dt = now - timeZero ;
      else
	timeDiffSecs (thisSessionStartTime, timeNow(), &dt) ;
      if ((dt + 10 > timeLapse) &&
	  isWriteAccess()) /* avoid rounding problems, and never substract on unsigned int  */
	{
	  VoidRoutine old = writeAccessChangeRegister (0) ;  /* avoid redraws */

	  messStatus (messprintf("Auto-Saving every %d minutes",timeLapse/60));
	  messdump("Autosaving if %d s elapsed done, last save was %d s ago", timeLapse, dt) ;
	  timeZero = 0 ;
	  sessionClose (TRUE) ;
	  sessionGainWriteAccess () ; /* grab again */
	  
	  now = time (0) ;
	  messdump("Autosaving done") ;
	  timeDiffSecs (thisSessionStartTime, timeNow(), &dt) ;
	  writeAccessChangeRegister (old) ;
	}
      else if (debug)
	messdump("Autosaving  if %d s elapsed not needed, last %d s ago  ", 
		 timeLapse,  dt) ;
      break ;

    default:
      messcrash ("Wrong call to sessionTimeSave, mode = %d", mode) ;
      break ;
    }
  return dt ;
}

/*************************************************************/
/* alarm routine */
static void  sessionAutoSaveAlarm (int dummy)
{
  unsigned int 
    dt = 0,                                /* time since last save */ 
    dtAlarm = 0,                          /* delay till next alarm */
    s1 = sessionAutoSave (-1, 0) ,                /* in seconds */
    s2 = 0 ;                          /* ask value after saving */
  BOOL debug = FALSE ;

  alarm ((unsigned int) 0) ;             /* disables autosaving */
  if (s1 > 0)
    {
      if (debug) messdump("AutoSaveAlarm calling sessionTimeSave (%d,0)", s1, 0) ;
      dt = sessionTimeSave (s1, 0) ;          /* save if necessary */
      s2 = sessionAutoSave (-2, 0) ;
      if (s2 > 0)
	{
	  dtAlarm = s2 > dt ? s2 - dt : s2 ;
	  if (dtAlarm < 10) dtAlarm = 10 ;
	  signal (SIGALRM, sessionAutoSaveAlarm) ; /* self register */
	  alarm (dtAlarm) ;                          /* raise alarm */
	  if (debug)
	    messdump ("Autosaving alarm expected in %d s", dtAlarm) ;
	}
      else if (debug)
	messdump ("Autosaving alarm cancelled", dtAlarm) ;
    }
  else if (debug)
    messdump ("Autosaving alarm cancelled", dtAlarm) ;
  return ;
}

/*************************************************************/
/*
 * sessionAutoSave (int saveInterval, int inactivityInterval)
 *
 * Example
 *  From anywhere in the code call
 *    sessionAutoSave (1800, 300)
 * As a result the system will autosave every half hour
 * or at the end of the first 5 minutes of innactivity
 * or of course on explicit user request
 *
 ************** Public calls 
 * saveInterval >  0 : The maximal delay before a save
 * saveInterval =  0 : cancels
 * inactivityInterval : The maximal inactivity delay before a save
 * Activity is raised by the acedb kernel calls to bsCreate
 ************** Kernel calls 
 * negative values are private to sessionAutoSave/acedb-kernel
 *** saveInterval = -1 : returns maximal saveTime interval
 *** saveInterval = -2 : returns delay till next alarm
 *** saveInterval = -3 : ace kernel starts saving
 *** saveInterval = -4 : ace kernel finished saving
 *** saveInterval = -7 : returace kernel finished saving
 *** saveInterval = -8 : signals modif, called by bsUpdate
 *** saveInterval = -9 : signals activity, called by bsCreate
 ************** State
 *** 0: unset
 *** 1: saving while unset
 *** 2: saving while set
 *** 3: dead
 *** 4: inactif
 *** 5: actif
 ************* Algorithm
 * This routine is just a finite state machine, actual job done by 
 * sessionAutoSaveAlarm  which self registers in alarm() and calls
 * sessionTimeSave which saves if last save too old
 */

unsigned int sessionAutoSave (int saveInterval, int inactivityInterval)
{
  static int 
    state = 0,  
    longDelay = 0 ,
    shortDelay = 0 ;
  int noDelay = 1 ;
  BOOL debug = TRUE ;

  if (TRUE || /* may 8 2003, disable autosave for a few days to see if it generates X crashes */ 
      prefValue ("ACEDB_NO_AUTOSAVE")) return 0 ;
  if (saveInterval == -9)        /* record activity, called by all bsCreate */
    switch (state)
      {
      case 0:                  /* autoSave not set, ignore */
      case 2:                  /* saving while set , ignore */
      case 1:                  /* saving while unset, ignore */
	return 0 ;
      case 3:                  /* dead, ignore */
	return 3 ; 
      case 4:                  /* inactif */
	state = 5 ;            /* new activity occured */
      case 5:                  /* actif */
	return 5 ;             /* return actif */
      default:
	messcrash ("sessionAutoSave in wrong state %d", state) ;
	return state ;
      }
  
  if (saveInterval == -8)        /* record modif, called by all bsUpdate */
    switch (state)
      {
      case 0:                  /* autoSave not set, ignore */
      case 2:                  /* saving while set , ignore */
      case 1:                  /* saving while unset, ignore */
	return 0 ;
      case 3:                  /* raise from dead */
      if (debug) messdump("sessionAutoSave(-8) calling sessionTimeSave (0,2)") ;
	sessionTimeSave (0, 2) ;    /* synchronize */
	sessionAutoSaveAlarm (0) ;
      case 4:                  /* inactif */
	state = 5 ;            /* new activity occured */
      case 5:                  /* actif */
	return 5 ;             /* return actif */
      default:
	messcrash ("sessionAutoSave in wrong state %d", state) ;
	return state ;
      }
  
  if (saveInterval == -1)      /* gives maximal delay since previous save */
    switch (state)
      {
      case 0:                  /* autoSave not set */
      case 1:                  /* saving while unset, ignore */
      case 2:                  /* saving while set, ignore */
      case 3:                  /* dead, goto sleep */
	return 0 ;
      case 4:                  /* inactif, save immediately */
	return noDelay ;
      case 5:                  /* actif */
	state = 4 ;            /* sets to inactif, crucial switch of whole system */	
	return longDelay ;
      default:
	messcrash ("sessionAutoSave in wrong state %d", state) ;
	return state ;
      }
  
  if (saveInterval == -2)      /* gives delay till next alarm */
    switch (state)
      {
      case 0:                  /* autoSave not set */
      case 1:                  /* saving while unset, ignore */
      case 2:                  /* saving while set, ignore */
      case 3:                  /* dead, go to sleep */
	return 0 ;        
      case 4:                  /* inactif */
      case 5:                  /* actif */
	return shortDelay ;
      default:
	messcrash ("sessionAutoSave in wrong state %d", state) ;
	return state ;
      }
  
  if (saveInterval == -3)        /* saving starts */
    switch (state)
      {
      case 0:                  /* autoSave not set */
	state = 1 ;            /* set saving while unset */
      case 1:                  /* saving while unset, ignore */
	return 1 ;
      case 2:                  /* saving while set, ignore */
      case 3:                  /* dead */
      case 4:                  /* inactif */
      case 5:                  /* actif */
	state = 2 ;
	return  2 ;
      default:
	messcrash ("sessionAutoSave in wrong state %d", state) ;
	return state ;
      }
  
  if (saveInterval == -4)        /* saving ends */
    switch (state)
      {
      case 0:                  /* autoSave not set */
	return 0 ;
      case 1:                  /* saving while unset, revert to unset */
	state = 0 ;
      if (debug) messdump("sessionAutoSave(-4) calling sessionTimeSave (0,2)") ;
	sessionTimeSave (0, 2) ;   /* synchronize */
	return 0 ;
      case 2:                  /* saving while set, revert to dead */
      case 3:                  /* dead */
      case 4:                  /* inactif */	
	state = 4 ;   /* mhmp 05.06.02  state = 3 */
	return  3 ;
      case 5:                  /* actif */
	state = 3 ; 
	if (debug) messdump("sessionAutoSave(-5) calling sessionTimeSave (0,2)") ;
	sessionTimeSave (0, 2) ;   /* synchronize */
	sessionAutoSaveAlarm (0) ;  
	return  3 ;
      default:
	messcrash ("sessionAutoSave in wrong state %d", state) ;
	return state ;
      }

  if (saveInterval <= 0)         /*  zero or silly values, disables autosave */
    switch (state)
      {
      case 0:                  /* autoSave not set, ignore */
      case 1:                  /* saving while unset ignore */
	return 0 ;
      case 2:                  /* saving while set, unset with delay */
	state = 1 ;
	longDelay = shortDelay = 0 ;
	return 1 ;
      case 3:                  /* autosave set, unset */
      case 4:                 
      case 5:
	state = 0 ;	
	longDelay = shortDelay = 0 ;
	sessionAutoSaveAlarm (0) ; 
	return 1 ;
      default:
	messcrash ("sessionAutoSave in wrong state %d", state) ;
	return state ;
      }
  
  /* verify parameters are reasonable */
  if (saveInterval <= 3 * 60) 
    saveInterval = 180 ; /* 300 ; */
  if (inactivityInterval < 60)
    inactivityInterval = 60 ;
  
  if (shortDelay != inactivityInterval ||
      longDelay != saveInterval)
    state = 0 ;  /* reregistering needed */
  
  if (state != 5)
    {
      state = 5 ;        /* short wait needed */
      shortDelay = inactivityInterval ; 
      longDelay = saveInterval ;
      sessionAutoSaveAlarm (0) ;   /* register */
    }

  return state ;
}

/*****************************************************************/
/********** AutoSaving package end *******************************/
/*****************************************************************/

void sessionDoSave (BOOL reGainWriteAccess)  
     /* general public function */
{ 
  /* saves even if i had forgotten to take write access before,
     there can be no conflict, because if another user has written
     i won't be able to grab the write access */

  messStatus ("Saving") ;

  if (!isWriteAccess())  
    sessionGainWriteAccess () ; /* try to grab it */
  if (isWriteAccess())
    { sessionClose (TRUE) ;
      if (reGainWriteAccess) 
	sessionGainWriteAccess () ; /* grab again */
    }
}

/*************************************************************/

KEY sessionUserKey (void)
{ 
  char timeBuf[25] ;

  if (!thisSession.userKey)
   { 	/* next line ensures that the lex/voc for _VUserSession is
	   loaded, which appears to be necessary during booting
	   cf the double lexaddkey during bootstrap in lexsubs.c */
     lexword2key ("junk", &thisSession.userKey, _VUserSession) ;
     
     if (!thisSessionStartTime)
       /* this happens when this function is called very early on 
	  during bootstrap (during readModels at that stage) */ 
       thisSessionStartTime = timeNow();

     if (!lexaddkey (messprintf ("%s_%s",
				 timeShow(thisSessionStartTime, timeBuf, 25), 
				 getLogin(TRUE)), 
		     &thisSession.userKey, _VUserSession))
       { /* we've had that key already, now make a unique one */
	 int i = 0 ;
	 while (!lexaddkey (messprintf ("%s.%d_%s", 
					timeShow (thisSessionStartTime, timeBuf, 25),
					++i,
					getLogin(TRUE)), 
			    &thisSession.userKey, _VUserSession)) ;
       }
   }
  return thisSession.userKey ;	
/* ALPHA does not seem to keep global access reliably up to date */
}

/*************************************************************/
/****************************************************************/


/*************************************************************/
/********************* Readlock manager **********************/
/*************************************************************/

const int GRACE_PERIOD_HOURS = 12;	/* remove readlock files
					   older than that */

static char readlock_filename[MAXPATHLEN] = "";

/*************************************************************/
/* creates the directory database/readlocks/ and gives
   everybody write access to it. This directory will contain
   lock files, that specify which sessions are being looked at,
   so they aren't destroyed under an old process' feet */
/*************************************************************/

static BOOL readlockCreateDir (void)
{
  static int readlock_status = 0; /* 0 : not tried before
				     1 : everything OK
				     2 : had problems, don't try again */
  char *readlock_dir;

  if (readlock_status != 0)
    return (readlock_status == 1 ? TRUE : FALSE);

  if (sessionFilName("database/readlocks", 0, "wd"))
    {
      /* a writable directory already exists */
      readlock_status = 1;	/* no need to check again next time */
      return TRUE;			
    }

  readlock_dir = 
    strnew(messprintf ("%sreadlocks/",
		       sessionFilName("database/", 0, "rd")), 0);
  
  /* make the readlocks directory with rwxrwxrwx permissions,
     so that any other user can create and remove files from there */
  if (mkdir (readlock_dir, 0000777) == -1)
    {
      /* unsuccessful completion of command */
      /* This will happen, if we us this code on a database, which
	 the administrator (the one with all the permissions) hasn't
	 initialised yet - that would have created the directory
	 correctly, even for users to use, that normally only have
	 read-only permissions. */
      
      /* no big error-message, we'll have to silently accept it */
      messfree (readlock_dir);
      readlock_status = 2;	/* don't try (and fail) again */
      return FALSE;
    }

  /* I noticed that the mkdir(2) won't assign the right permissions
     at all times. That may have to do with different umask(2)
     settings. So we explicitly change the permissions again. */    
  if (chmod (readlock_dir, 0777) == -1)
    {
      messfree (readlock_dir);
      readlock_status = 2;	/* don't try (and fail) again */
      return FALSE;
    }

  
  messfree (readlock_dir);
  readlock_status = 1;		/* everything went OK */

  return TRUE;		/* managed to create directory alright */
} /* readLockCreateDir */



/*************************************************************/
/* for the given session, create the readlock file
   the file-name has the format <session>.<host>.<pid>
   and the file is in the directory database/readlocks/... */
/*************************************************************/
static BOOL readlockCreateFile (void)
{
  char host_name[100];
  int pid;
  int fd;
  char *buffer;
  char timeBuf[25] ;

  if (dropWriteAccess) return TRUE;

  if (!(readlockCreateDir()))
    return FALSE;

  gethostname (host_name, 100) ;
  pid = getpid () ; 
  
  sprintf (readlock_filename, "%s%d.%s.%d",
	   sessionFilName("database/readlocks/", 0, "wd"),
	   thisSession.session, host_name, pid);
	   
  /* create the readlockfile with rw-rw-rw- permission. */
  /* That makes it possible for other processes to
     remove expired readlock files that this proces may leave behind */

  fd = open(readlock_filename, O_RDWR | O_CREAT | O_SYNC, 0666);

  if (fd == -1)
    {
      strcpy (readlock_filename, "");
      return FALSE;
    }

  buffer = messprintf
    ("Readlock file to prevent destruction of "
     "sessions still in use by other processes\n"
     "Created: %s\n"
     "User: %s\n"
     "Program: %s\n"
     "Version: %s\n"
     "Linkdate: %s\n",
     timeShow(timeNow(), timeBuf, 25), 
     getLogin(TRUE), 
     messGetErrorProgram(), 
     aceGetVersionString(), 
     aceGetLinkDateString());
	   
  write (fd, buffer, strlen(buffer));
  close (fd);

  /* I noticed that the open(2) won't assign the right permissions
     at all times. That may have to do with different umask(2)
     settings. So we explicitly change the permissions again. */    
  if (chmod (readlock_filename, 0666) == -1)
    {
      strcpy (readlock_filename, "");      
      return FALSE;
    }

  return TRUE;
} /* readLockCreateFile */



/*************************************************************/
/* when you save a session, we get a new session number,
   so we rename the current readlock-file to contain
   the new session-number in the file name */
/*************************************************************/
static BOOL readlockUpdateFile (int new_session_num)
{
  char host_name[100];
  int pid;
  char new_filename[MAXPATHLEN] = "";

  if (dropWriteAccess) return TRUE;

  if (new_session_num == thisSession.session)
    messcrash("readlockUpdateFile() - new_session_num == thisSession.session");

  if (strlen(readlock_filename) == 0)
    return FALSE;

  gethostname (host_name, 100) ;
  pid = getpid () ; 

  sprintf (new_filename, "%s/%d.%s.%d",
	   sessionFilName("database/readlocks/", 0, "r"),
	   new_session_num, host_name, pid);
	   
  if (rename (readlock_filename, new_filename) == -1)
    {
      if (debugReadLock)
	messerror ("failed to rename readlock file %s to %s (%s)",
		   readlock_filename,
		   new_filename,
		   messSysErrorText());

      /* don't know how we lost it, but we must have a new readlock */
      return readlockCreateFile();
    }

  strcpy (readlock_filename, new_filename);

  return TRUE;
} /* readlockUpdateFile */

/*************************************************************/
/* called when the process exists, readlock-files that
   don't get destroyed this way because of an uncontrolled
   crash will be dealt by readlockIsSessionDestructible()
   if needed. */
/*************************************************************/
BOOL readlockDeleteFile (void)
{
  if (dropWriteAccess) return TRUE;

  if (!(readlockCreateDir()) || strlen(readlock_filename) == 0)
    return FALSE;

  if (unlink(readlock_filename) == -1)
    {
      if (debugReadLock)
	messerror ("failed to remove readlock file %s (%s)",
		 readlock_filename,
		 messSysErrorText());
      return FALSE;
    }

  return TRUE;
} /* readlockDeleteFile */


/*************************************************************/
/* We will create an array of session numbers that still have
   a valid readlock. Sessions with those numbers will then be
   treated with caution, because other processes are still
   looking at that data of those older sessions. */
/*************************************************************/
static Array readlockGetSessionNumbers (void)
{
  Array readlockedSessions;
  Array dirArray;
  char lockfile_name[MAXPATHLEN];
  char *cp1, *cp2, *cp3;
  char *dirEntry;
  int  file_session_num;
  char file_host_name[100];
  int  file_pid;
  int n, hourAge;

  if (!(readlockCreateDir()))
    return 0;			/* no readlocks being used */

  dirArray = filDirectoryCreate
    (sessionFilName("database/readlocks", 0, "rd"), "", "r");

  if (!dirArray)
    {
      /* this shouldn't really go wrong, because we've managed to
	 open or create the directory already */
      if (debugReadLock)
	messerror ("failed to open directory database/readlocks/ (%s)",
		   messSysErrorText());
      
      return 0;
    }

  if (arrayMax(dirArray) == 0)	/* nothing is readlocked */
    return 0;

  readlockedSessions = arrayCreate (arrayMax(dirArray), int);

  for (n = 0; n < arrayMax(dirArray); n++)
    /* start at 1 to skip . entry */
    {
      dirEntry = arr(dirArray, n, char*);

      strcpy (lockfile_name, 
	      sessionFilName ("database/readlocks/", 0, "rd"));
      strcat (lockfile_name, dirEntry);

      /* don't consider the readlock file of this process */
      if (strcmp (lockfile_name, readlock_filename) == 0)
	continue;
	
      /* extract the <sessionnum>.<hostname>.<pid> parts from filename */

      cp1 = dirEntry;
      /* find the first dot in the name */
      if (!(cp2 = strstr(cp1, "."))) continue;

      *cp2 = 0; 		/* cp1 now ends at the dot */
      ++cp2;			/* cp2 starts just after it */

      /* find the first dot from the end */
      cp3 = cp2 + strlen(cp2);
      while (cp3 > cp2 && *cp3 != '.') --cp3;

      if (cp3 == cp2) continue;

      *cp3 = 0;			/* cp2 now ends at the 2nd dot */
      ++cp3;			/* cp3 starts just after */

      if (!*cp1 || !*cp2 || !*cp3) continue;

      if (sscanf (cp1, "%d", &file_session_num) != 1) continue;
      if (sscanf (cp2, "%s", file_host_name) != 1) continue;
      if (sscanf (cp3, "%d", &file_pid) != 1) continue;

      hourAge = 0;

      if (!(filAge (lockfile_name, "", 0, 0, 0, &hourAge, 0, 0)))
	continue;		/* couldn't get age, don't consider
				   it a valid readlock then */

      if (hourAge < GRACE_PERIOD_HOURS) /* still locked then */
	{
	  int max = arrayMax(readlockedSessions);
	  /* get session-number from file-name */
	  /* terminate the filename at the position of the dot
	     then read a number from that first part of the string */
	  array(readlockedSessions, max, int) = file_session_num;
	}
      else			/* lock has expired */
	{
	  if (unlink (lockfile_name) == 0)
	    {
	      /* successfully removed */
	      messdump
		("Readlock manager : "
		 "session %d readlock of process %d on %s "
		 "is older than %d hours. readlock removed.",
		 file_session_num, file_pid, file_host_name,
		 GRACE_PERIOD_HOURS);
	    }
	}
    }

  filDirectoryDestroy (dirArray);

  return readlockedSessions;
} /* readlockGetSessionNumbers */

/******************************************************************/
/*********************** eof **************************************/
