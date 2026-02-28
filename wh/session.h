/*  File: session.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
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
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  1 13:14 1998 (fw)
 * * Sep 25 12:24 1998 (edgrif): Added prototype for getLogin()
 * * Richard Bruskiewich
 *		-	declare function prototypes for setuid() et al. for WIN32 & MAC's
 *		-	also can put extern declarations for ruid and euid here
 * * Jun  3 19:08 1996 (rd)
 * Created: Fri Sep 25 12:24:38 1998 (edgrif)
 *-------------------------------------------------------------------
 */

/* $Id: session.h,v 1.5 2008/11/18 02:53:39 mieg Exp $ */

                     /* session.h                             */
                     /* public functions of sessionsubs.c     */
                     /* handling the session.                 */

#ifndef DEFINE_SESSION_H
#define DEFINE_SESSION_H


/* On all but the macintosh the new acedb_version code is used to get the release numbers. */
/* When I've found a macintosh to test on I will replace the code on that to.              */
#if defined(MACINTOSH) 

#if defined(ACEDB5)
#define MAINCODERELEASE 5
#define SUBCODERELEASE 0
#define DATE_OF_CODE_RELEASE "Version 5 is UNRELIABLE,  DO NOT USE IT"
#else
#define MAINCODERELEASE 4
#define SUBCODERELEASE 6
#define DATE_OF_CODE_RELEASE "July 4, 1997"
#endif

#endif /* MACINTOSH */



#ifndef DEF_DISK
#include "disk.h"
#endif

    /* sizeof(SESSION) must remain smaller than BLOCKSIZE */
typedef struct {
                 int session ; /* unique identifier */
                 DISK gAddress ;
		 int mainCodeRelease, subCodeRelease,
		     mainDataRelease, subDataRelease ;
                 char name[80], title[255] ;
		 KEY key, from, upLink, userKey ;
                 int index_version ;
	       } SESSION ;

/************************************************************/

extern SESSION thisSession ;

BOOL aceInit (const char *ace_path) ;/* launches the acedb kernel */
BOOL  aceQuit(BOOL commit) ;         /* terminates and closes the disks */


void  sessionShutDown (void) ;
void  sessionRegister (void) ;
void  sessionClose (BOOL doSave) ;
unsigned int  sessionAutoSave (int s, int tInactivity) ; 
              /* Sets autosaving every s seconds
	       * or tInactivity  seconds after last call to sessionAutoSave
	       * if s == 0, disables autosaving 
	       * iff s < 0, just returns registered interval
	       * ATTENTION: grabs alarm() */
unsigned int sessionTimeSave (unsigned int s, int mode) ;
                        /* mode == 0 :saves and regains write access 
			 * if s seconds have  elapsed since last save 
			 * returns time since last save 
			 * mode == 1,2 private
			 */
char *sessionFilName (char *name, const char *ending, const char *spec) ;

BOOL  isWriteAccess (void) ;
void  sessionDoSave (BOOL reGainWriteAccess);   /* tries to save if 
						   write access
						   is allowed, will 
						   re-gain write
						   access if required. */
char *sessionDbName(void);	/* name of database,
				 used to title main-window etc. */

BOOL readlockDeleteFile (void);
BOOL sessionSetUserSessionUser(char *user);		    /* Username must be alphanumerics. */
mytime_t sessionUserSession2Time(KEY userSession);
mytime_t sessionStartTime(void);

/****************** write access functions *************************/

BOOL  sessionGainWriteAccess (void); /* Get if possible */

void  sessionReleaseWriteAccess(void); /* drop write lock and access */

void  sessionForbidWriteAccess (void); /* can never re-gain
					 write access again */

BOOL  writeAccessPossible(void);

VoidRoutine writeAccessChangeRegister (VoidRoutine func); /* register
							     a function
							     to be 
							     called,
							     if write
							     access
							     changes */

/*******************************************************************/

KEY   sessionUserKey (void) ;

char *getLogin (BOOL isReal) ;	/* also used in acedbgraph.c */

extern uid_t ruid, euid;

#if defined(MACINTOSH) || defined(WIN32)
void seteuid (uid_t uid) ;
uid_t getuid () ;
uid_t geteuid () ;
#endif        /* defined(MACINTOSH) || defined(WIN32) */

/************************************/
extern BOOL swapData;

/************************************/

/********************** ace server callbacks   *******************************/
/*                                                                           */
/* The server needs to do additional processing at various stages of session */
/* processing (e.g. init/termination) it does this via callbacks that it     */
/* registers at server start up, prior to making the call to sessionInit().  */
/*                                                                           */
/* NOTE the philosophy is that the server deserves to be a special case of   */
/* initialising acedb and so it is acceptable that it has its own callback   */
/* interface.                                                                */
/*                                                                           */
/* sessionRegisterServerCB() is called before sessionInit().                 */
/*                                                                           */
typedef void (*ServerCB)(void *server_pointer) ;
typedef void (*ServerMesgCB)(void *server_pointer, char *mesg) ;

typedef struct _ServerContextStruct
{
  union u_funcs
  {
    ServerCB func ;
    ServerMesgCB mesgfunc ;
  } funcs ;

  void *data ;
} ServerContextStruct ;


typedef struct _ServerSessionCBstruct
{
  ServerContextStruct init ;
  ServerContextStruct checkFiles ;
  /* struct messContextStruct mess_out ;			    Directly register with messout */
  /* struct messContextStruct mess_err ;			    Directly register with messerror */
  ServerContextStruct exit ;
  ServerContextStruct crash ;
  ServerContextStruct sighup ;
} ServerSessionCBstruct, *ServerSessionCB ;

void sessionRegisterServerCB(ServerSessionCB callbacks) ;
ServerSessionCB sessionGetServerCB(void) ;
void monDnaReport (int *nDnap, int *nBasep) ;

#endif	/* DEFINE_SESSION_H */
 
