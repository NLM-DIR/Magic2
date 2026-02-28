/*  File: xclientmain.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Main() of the complete program
        Provides Initialisation and Termination routines.
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 11 12:12 1998 (fw)
 * * Dec  2 15:54 1998 (edgrif): Add code to return build time of this module.
 * * Oct 15 10:50 1998 (edgrif): Add new call for graphics initialisation.
 * Created: Tue Oct 22 12:51:04 1991 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: xclientmain.c,v 1.4 2011/10/09 03:06:40 mieg Exp $ */

#include "acedb.h"

#include "display.h"
#include "acedbgraph.h"
#include "main.h"
#include "session.h"
#include "systags.h"
#include "banner.h"
#include "model.h"
#include "version.h"

/************************************************************/

static void xaceclientCrashRoutine (void) /* called by acedbCrash(),
				       which is called by messcrash()
				       special version for xace instead
				       of defaultAceCrashRoutine() */
{
  graphCleanUp () ;  /* kills all open windows,
			forces displayed objects to cache */

  writeAccessChangeRegister (0); /* avoid unnecessary re-draws
				    just before exiting */

  /* release read/write locks and clean up temp-files 
     is all done by aceQuit() */
  if (isWriteAccess() && 
      messQuery ("You did not save your work, should I ?"))
    aceQuit (TRUE);
  else
    aceQuit (FALSE);
   
  xClientClose () ;

  /* graphFinish () ;
     probably not needed anyway, it causes problems for now */

  if (!getenv("ACEDB_NO_BANNER"))
    printf ("\n\n A bientot\n") ;	/* use printf just as bannerWrite does */
} /* xaceclientCrashRoutine */

/*************************************************/

static BOOL xaceQuery (char* text)
     /* callback for messQuery */
     /* special version for graphical apps, that use messStatus */
{
  BOOL queryResult;

  messStatus("Please Answer") ;
  queryResult = graphQuery(text) ; 
  messStatus("Please wait") ;
  return queryResult;
}

/******************  ACEDB  ***********************/
/**************************************************/

/* Defines a routine to return the compile date of this file.                */
UT_MAKE_GETCOMPILEDATEROUTINE()


int main (int argc, char **argv)
     /* main function for xaceclient */
{
  extern VoidRoutine messcrashroutine;
  messErrorInit(argv[0]);   /* Record program name for crash messages */

  bannerWrite (bannerMainStrings ("xaceclient", TRUE, FALSE)) ;

  freeinit () ;			/* must come before graphInit */
  freespecial ("\n\t\"/%\\@") ;	/* No subshells ($ removed) */


  /* Initialise the graphics package and its interface for acedb */
  acedbAppGraphInit(&argc, argv) ;

  pickGetArgs (&argc, argv) ;   /* allow pick to extract starting class and template - srk */

  aceInit (argc>1 ? argv[1]: 0) ;

  xClientInit () ;

  messcrashroutine = xaceclientCrashRoutine; /* instead of
						defaultAceCrashRoutine,
						as set in aceInit() */

  bannerWrite (bannerDataInfoStrings ()) ; /* must come after aceInit */

  messQueryRegister(xaceQuery); /* graphQuery wrapped in messStatus */

  writeAccessChangeRegister (pickDraw);	/* re-draw if write access 
					   changes to adapt buttons 
					   and menus */

  modelChangeRegister (pickReReadClasses); /* change class menu
					      when models change */

  quitFunctionRegister (xaceclientCrashRoutine); /* same close actions
						    for quit/exit button
						    as for crash */

  pickCreate();			/* start main window */

  graphLoop (FALSE) ;		/**** MAIN LOOP:  does all the work ****/

  /* shouldn't get here */
  xaceclientCrashRoutine();
  return(EXIT_FAILURE) ;
}

/**************************************************/
/**************************************************/
 
 
