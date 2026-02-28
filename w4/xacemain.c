/*  File: xacemain.c
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
 * Last edited: Dec 21 11:42 1998 (edgrif)
 * * Dec  3 12:16 1998 (edgrif): Remove externs from within function.
 * * Dec  2 15:52 1998 (edgrif): Add code to record compile time of this module.
 * * Oct 15 10:49 1998 (edgrif): Add new call for graphics initialisation.
 * * Sep  9 16:51 1998 (edgrif): Add call to messErrorInit to record
 *              program name (xace probably) for use in messcrash calls.
 * Created: Tue Oct 22 12:51:04 1991 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: xacemain.c,v 1.7 2015/10/12 22:16:16 mieg Exp $ */

#include "acedb.h"

#include "display.h"
#include "acedbgraph.h"
#include "main.h"
#include "session.h"
#include "systags.h"
#include "banner.h"
#include "model.h"
#include "version.h"

extern VoidRoutine messcrashroutine;

/************************************************************/

static void xaceCrashRoutine (void) /* called by acedbCrash(),
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
   
  /* graphFinish () ;
     probably not needed anyway, it causes problems for now */

  if (!getenv("ACEDB_NO_BANNER"))
    printf ("\n\n A bientot\n") ;	/* use printf just as bannerWrite does */
} /* xaceCrashRoutine */

/*************************************************/

static BOOL xaceQuery (const char* text)
     /* callback for messQuery */
     /* special version for graphical apps, that use messStatus */
{
  BOOL queryResult;

  if ( graphActive ())  /* may fail if the windowing system is dead */
    {
      messStatus("Please Answer") ;
      queryResult = graphQuery(text) ; 
      messStatus("Please wait") ;
    }
  else  /* fall back on basic text interface */
    queryResult = freequery (text) ;
  return queryResult;
}


/* Defines a routine to return the compile date of this file.                */
UT_MAKE_GETCOMPILEDATEROUTINE()



/******************  ACEDB  ***********************/
/**************************************************/

int main (int argc, char **argv)
     /* main function used for xace and xacembly */
{
  messErrorInit(argv[0]);   /* Record program name for crash messages */

#ifdef ACEMBLY
  bannerWrite (bannerMainStrings ("xacembly", TRUE, TRUE)) ;
#else
  bannerWrite (bannerMainStrings ("xace", TRUE, FALSE)) ;
#endif



  /* Initialise the graphics package and its interface for acedb */


  freeinit () ;			/* must come before graphInit */
  freespecial ("\n\t\"/%\\@") ;	/* No subshells ($ removed) */
  acedbAppGraphInit(&argc, argv) ;
  aceInit (argc>1 ? argv[1]: 0) ;
  pickGetArgs (&argc, argv) ;   /* allow pick to extract starting class and template - srk */

  if (0)
    {
      char *cp1 = "Göndör" ;
      messout (cp1) ;
      exit (0) ;
    }


  messcrashroutine = xaceCrashRoutine;	/* instead of
					   defaultAceCrashRoutine,
					   as set in aceInit() */

  bannerWrite (bannerDataInfoStrings ()) ; /* must come after aceInit */

#ifdef ACEMBLY
  acemblyInit () ;
#endif

  messQueryRegister(xaceQuery); /* graphQuery wrapped in messStatus */

  writeAccessChangeRegister (pickDraw);	/* re-draw if write access 
					   changes to adapt buttons 
					   and menus */

  modelChangeRegister (pickReReadClasses); /* change class menu
					      when models change */

  quitFunctionRegister (xaceCrashRoutine); /* same close actions
					      for quit/exit button
					      as for crash */

  pickCreate();			/* start main window */

  graphLoop (FALSE) ;		/**** MAIN LOOP:  does all the work ****/

  /* shouldn't get here */
  xaceCrashRoutine();
  return(EXIT_FAILURE) ;
}

/**************************************************/
/**************************************************/
