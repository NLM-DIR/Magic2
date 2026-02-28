/*  File: graphxt.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Xt/Xaw level of X version of graph package
 * Exported functions: many
 * HISTORY:
 * Last edited: May 12 13:58 1999 (edgrif)
 * * May 12 13:56 1999 (edgrif): Added mapping event handler for busy
 *              cursor control, added code to popdown popup menus.
 * * Apr  7 10:12 1999 (edgrif): Fix SANgc03918, menu position problem.
 * * Mar 24 14:35 1999 (edgrif): Fix SANgc03641, see NOMENU flag in menuCall.
 * * Dec 21 10:43 1998 (edgrif): Make sure busy cursor is NOT shown for
 *              pointer box dragging type stuff.
 * * Dec 15 23:41 1998 (edgrif): Moved busy cursor stuff from event loop,
 *              to the nnnnCall routines.
 * * Nov 26 09:56 1998 (edgrif): Fixed small buglet in cut/paste.
 * * Nov 19 10:38 1998 (edgrif): Cleaned up X error handling: added handlers
 *              for X IO errors & Xt errors..both fatal. Added X11r6
 *              signal safe code for CNTL-C, made pre-X11r5 code
 *              as safe as possible.
 * * Nov 17 11:08 1998 (edgrif): Fred added stuff for SunOS not having
 *              XPointer, tidied up to put all SunOS stuff at top of file.
 * * Nov 12 14:25 1998 (edgrif): Added code to get the real size of the
 *              screen instead of guessing, commented it out because all
 *              the window sizes get messed up.
 * * Nov 12 14:09 1998 (edgrif): Replaced old cut/paste code to BUFFER0
 *              with standard cut/paste to XA_PRIMARY selection buffer.
 * * Oct 13 14:42 1998 (edgrif): Replace ACEDB defs with function calls.
 * * Sep 24 16:51 1998 (edgrif): Altered user initiated abort in intProc
 *              to be a messExit instead of messcrash.
 * * Feb 19 12:28 1996 (daz)
 * * Jan 10 16:50 1996 (rd)
 * * Sep 20 14:19 1992 (rd): menu correctly positioned when scrollbar
 * * Jan  2 02:29 1992 (rd): exposeCall now restores gActive
 * Created: Thu Jan  2 02:27:24 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: graphxt.c,v 1.13 2017/01/17 21:15:06 mieg Exp $ */

#include "regular.h"
#include <stdnoreturn.h>

#include "graphxt_.h"
#include "key.h"		/* for keyboard stuff */
/*#include "mytime.h"*/		/*mhmp 24.11.97*/

/* #include <X11/Xaw/Cardinals.h> */
#include <X11/Xaw/Viewport.h>
#include <X11/Xaw/Scrollbar.h>
#include <X11/Xaw/SimpleMenu.h>
#include <X11/Xaw/SmeBSB.h>
#include <X11/Xaw/AsciiText.h>
#include <X11/Xaw/Dialog.h>
#include <X11/Xaw/Command.h>

#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
#define  XK_MISCELLANY					    /* UUggghhh, what is this for... */
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <X11/Intrinsic.h>				    /* for XtTranslateCoords - srk */
#include <X11/Xatom.h>					    /* Next 3 all for selections support. */
#include <X11/Xmu/Atoms.h>
#include <X11/Xmu/StdSel.h>
#include <X11/cursorfont.h>
#include <signal.h>





/* SunOS is very back level and needs lots of help to compile stuff...       */
#ifdef SUN

/* RD fix for X11R4 on Sun unaccessed library calls missing in link */
void get_wmShellWidgetClass (void) { }
void get_applicationShellWidgetClass (void) { }

/* Stupid Sun does not even know about this basic X Window pointer type...   */
typedef char *XPointer;

#endif



extern XFontStruct * getDefaultFont(void);		    /* from graphxlib.c */
extern BOOL cmapDebug;					    /* from graphxlib.c */
extern int preferredVisual;				    /* from graphxlib.c */


/* globals within the X graph package */

float  gPixAspect = 1.0 ;

XtAppContext app_con;	/* application context */
Widget root_widget;	/* parent of all */
static Associator lastLeft = 0;

int menuBox ;	/* a global from graphxt_.h used in xtsubs.c */



/* There is a nasty problem in which the leaveCall callback for a widget
   gets called after it has been detroyed, thus putting a dangling pointer
   into lastLeft, having failed to find a fix for this in X, I've implemented
   a hack in which such  pointers get put into destroyed, and then ignored.
   It is NOT pretty, but it works - srk */
static Associator destroyed = 0 ;

/*******/

static Arg args[100];		/* convenience variable */

/* typeName used for widget instance names - 
   corresponds to enum GraphType for resource setting
*/
static char* typeName[] = {"plain",
			   "text_scroll",
			   "text_fit",
			   "map_scroll",
			   "pixel_scroll",
			   "text_full_scroll",
			   "pixel_fit",
			   "text_full_edit"
			  } ;

/* screenx, screeny are for positioning and sizing created graphs
   SUN and other common screens are 1152x900 (1.3x900 = 1170)
   so this is a kind of default. These sizes should be 
   obtained from the X server, see graphDevInit(), currently this is not
   not done because all the window sizes get messed up.
*/
static float screenx = 900 ;
static float screeny = 900 ;

static Associator assw2g = 0 ;

typedef struct LoopInfoStruct {
  Graph_ g;
  VoidCharRoutine reportActivity ;
  BOOL isBlocking ;
  BOOL isDone ;
  int retval ;
} LOOPINFO ;
static Array loopInfo = 0 ;

static void fakeResize (Graph_ g) ;
static BOOL activeSet (Widget w) ;

/* Should menu popup button do a grab ?                                      */
static BOOL noGrab = FALSE ;
static BOOL doGrab = FALSE ; /* mieg: i reverse the logic */

BOOL installOwnColorMap = FALSE;



/*                     Xt Selection support (new cut & paste)                */
/*                                                                           */

/* we use this structure to pass results from the paste callback function to */
/* graphPasteBuffer().                                                       */
typedef struct graphSelectResult_ {
  BOOL finished ;					    /* TRUE when selection finished. */
  char *text ;						    /* NULL if selection failed. */
} graphSelectResult ;

/* We use this to pass the Cut string between callback functions registered  */
/* by graphPostBuffer()                                                      */
static char *graphSelectCutString_G = NULL ;



/*                    Interrupt signal support.                              */
/*                                                                           */
/* NOTE: this is only fully supported in X Windows 11.6, previous releases   */
/*       could crash if the signal handler tried to issue X calls (which     */
/*       our signal handler does).                                           */
/*                                                                           */
/*       I have fixed this for X11r6 and attempted to put in a warning for   */
/*       earlier releases...there is no real fix for earlier releases.       */
/*                                                                           */
#if XtSpecificationRelease >= 6
static XtSignalId sigintID_G ;
#else
static Boolean inSignalRoutine_G = FALSE ;
#endif



/*****************************************************************************/
/* registered functions.                                                     */

/* General keypress etc.                                                     */
static void exposeCall(Widget w, XEvent* ev, String *params, Cardinal *num_params);
static void mouseCall(Widget w, XEvent* ev, String *params, Cardinal *num_params);
static void pickCall(Widget w, XEvent* ev, String *params, Cardinal *num_params);
static void middleCall(Widget w, XEvent* ev, String *params, Cardinal *num_params);
static void keyboardCall(Widget w, XEvent* ev, String *params, Cardinal *num_params);
static void resizeCall (Widget w, XEvent* ev, String *params, Cardinal *num_params) ;
static void propertyCall (Widget w, XEvent* ev, String *params, Cardinal *num_params) ;
static void subDevMapCall (Widget w, XtPointer client_data, XEvent *event, Boolean *propagate) ;
static void menuCall (Widget w, XEvent* ev, String *params, Cardinal *num_params) ;
static void activeCall (Widget w, XEvent* ev, String *params, Cardinal *num_params) { activeSet (w) ; }
static void destroyCallback (Widget w, XtPointer client_data, XtPointer call_data);
static void enterCall (Widget w, XEvent* ev, String *params, Cardinal *num_params) ;
static void leaveCall (Widget w, XEvent* ev, String *params, Cardinal *num_params) ;
static void returnCall (Widget w, XEvent* ev, String *params, Cardinal *num_params) { ; } /* a stub for now */
static void blockingRaise (void) ; /* mhmp 12.11.97 */

/* Cut/Paste support.                                                        */
static void getPasteData(Widget w, XtPointer our_ptr,
			 Atom *selection, Atom *type,
			 XtPointer value, unsigned long *length, int *format) ;
static Boolean sendCutData(Widget w,  Atom *selection, Atom *target,
			    Atom *type, XtPointer *value, unsigned long *length, int *format) ;
static void lostCutSelection(Widget w, Atom *selection) ;
static void finishedCutSend(Widget w, Atom *selection, Atom *target) ;


/* Interrupt signal support.                                                 */
static void intProc (int signum) ;
#if XtSpecificationRelease >= 6
void sigintCallback(XtPointer client_data, XtSignalId *id) ;
#endif



/*****************************************************************/




/* action translations for the composite widget (= graph itself) */

static String actionTrans =	/* RD added Shift<Btn3Down> */
    "<Expose>:           exposeCall()\n\
     <ConfigureNotify>:  resizeCall()\n\
     <PropertyNotify>:   propertyCall()\n\
     <EnterNotify>:      enterCall()\n\
     <LeaveNotify>:      leaveCall()\n\
     Shift<Btn1Down>:         pickCall()\n\
     <Btn1Down>:         pickCall()\n\
     <Btn1Motion>:       mouseCall(\"1\")\n\
     <Btn1Up>:           mouseCall(\"2\")\n\
     <Btn2Down>:         middleCall()\n\
     <Btn2Motion>:       mouseCall(\"4\")\n\
     <Btn2Up>:           mouseCall(\"5\")\n\
     Shift<Btn3Down>:    menuCall()\n\
     Ctrl<Btn3Down>:     menuCall()\n\
     <Btn3Down>:         menuCall()\n\
     <Key>:              keyboardCall()\n";


static String viewportTrans =
    "<EnterNotify>:      enterCall()\n\
     <LeaveNotify>:      leaveCall()\n";

static String textFullActionTrans =
    "#replace \n\
     Ctrl<Key>A:         beginning-of-line() \n\
     Ctrl<Key>B:         backward-character() \n\
     Ctrl<Key>D:         delete-next-character() \n\
     Ctrl<Key>E:         end-of-line() \n\
     Ctrl<Key>F:         forward-character() \n\
     Ctrl<Key>G:         multiply(Reset) \n\
     Ctrl<Key>H:         delete-previous-character() \n\
     Ctrl<Key>J:         newline-and-indent() \n\
     Ctrl<Key>K:         kill-to-end-of-line() \n\
     Ctrl<Key>L:         redraw-display() \n\
     Ctrl<Key>M:         newline() \n\
     Ctrl<Key>N:         next-line() \n\
     Ctrl<Key>O:         newline-and-backup() \n\
     Ctrl<Key>P:         previous-line() \n\
     Ctrl<Key>T:         transpose-characters() \n\
     Ctrl<Key>V:         next-page() \n\
     Ctrl<Key>W:         kill-selection() \n\
     Ctrl<Key>Y:         insert-selection(CUT_BUFFER1) \n\
     Ctrl<Key>Z:         scroll-one-line-up() \n\
     Meta<Key>B:         backward-word() \n\
     Meta<Key>F:         forward-word() \n\
     Meta<Key>K:         kill-to-end-of-paragraph() \n\
     Meta<Key>Q:         form-paragraph() \n\
     Meta<Key>V:         previous-page() \n\
     Meta<Key>Y:         insert-selection(PRIMARY, CUT_BUFFER0) \n\
     Meta<Key>Z:         scroll-one-line-down() \n\
     :Meta<Key>d:        delete-next-word() \n\
     :Meta<Key>D:        kill-word() \n\
     :Meta<Key>h:        delete-previous-word() \n\
     :Meta<Key>H:        backward-kill-word() \n\
     :Meta<Key>\\<:      beginning-of-file() \n\
     :Meta<Key>\\>:      end-of-file() \n\
     :Meta<Key>]:        forward-paragraph() \n\
     :Meta<Key>[:        backward-paragraph() \n\
     ~Shift Meta<Key>Delete:    delete-previous-word() \n\
      Shift Meta<Key>Delete:    backward-kill-word() \n\
     ~Shift Meta<Key>BackSpace: delete-previous-word() \n\
      Shift Meta<Key>BackSpace: backward-kill-word() \n\
     <Key>Right:         forward-character() \n\
     <Key>Left:          backward-character() \n\
     <Key>Down:          next-line() \n\
     <Key>Up:            previous-line() \n\
     <Key>Delete:        delete-previous-character() \n\
     <Key>BackSpace:     delete-previous-character() \n\
     <Key>Linefeed:      newline-and-indent() \n\
     <Key>Return:        newline() \n\
     <Key>:              insert-char() \n\
     <EnterNotify>:      enterCall()\n\
     <LeaveNotify>:      leaveCall()\n\
     <Btn1Down>:         select-start() \n\
     <Btn1Motion>:       extend-adjust() \n\
     <Btn1Up>:           extend-end(PRIMARY, CUT_BUFFER0) \n\
     <Btn2Down>:         insert-selection(PRIMARY, CUT_BUFFER0) \n\
     <Btn3Down>:         menuCall()\n\
     Shift<Btn3Down>:    menuCall()\n\
     Ctrl<Btn3Down>:     menuCall()\n";

/* Not yet implemented:
     Ctrl<Key>R:         search(backward) \n\
     Ctrl<Key>S:         search(forward) \n\
     Meta<Key>I:         insert-file() \n\
*/



#ifdef ED_G_NEVER_INCLUDE_THIS_CODE

/* CURRENTLY UNUSED SO DEFINED OUT....                                       */
static String textLineActionTrans =
    "#replace \n\
     Ctrl<Key>A:         beginning-of-line() \n\
     Ctrl<Key>B:         backward-character() \n\
     Ctrl<Key>D:         delete-next-character() \n\
     Ctrl<Key>E:         end-of-line() \n\
     Ctrl<Key>F:         forward-character() \n\
     Ctrl<Key>G:         multiply(Reset) \n\
     Ctrl<Key>H:         delete-previous-character() \n\
     Ctrl<Key>K:         kill-to-end-of-line() \n\
     Ctrl<Key>L:         redraw-display() \n\
     Ctrl<Key>M:         returnCall() \n\
     Ctrl<Key>T:         transpose-characters() \n\
     Ctrl<Key>W:         kill-selection() \n\
     Ctrl<Key>Y:         insert-selection(CUT_BUFFER1) \n\
     Meta<Key>B:         backward-word() \n\
     Meta<Key>F:         forward-word() \n\
     Meta<Key>Y:         insert-selection(PRIMARY, CUT_BUFFER0) \n\
     :Meta<Key>d:        delete-next-word() \n\
     :Meta<Key>D:        kill-word() \n\
     :Meta<Key>h:        delete-previous-word() \n\
     :Meta<Key>H:        backward-kill-word() \n\
     ~Shift Meta<Key>Delete:    delete-previous-word() \n\
      Shift Meta<Key>Delete:    backward-kill-word() \n\
     ~Shift Meta<Key>BackSpace: delete-previous-word() \n\
      Shift Meta<Key>BackSpace: backward-kill-word() \n\
     <Key>Right:         forward-character() \n\
     <Key>Left:          backward-character() \n\
     <Key>Delete:        delete-previous-character() \n\
     <Key>BackSpace:     delete-previous-character() \n\
     <Key>Return:        returnCall() \n\
     <Key>:              insert-char() \n\
     <Btn1Down>:         select-start() \n\
     <Btn1Motion>:       extend-adjust() \n\
     <Btn1Up>:           extend-end(PRIMARY, CUT_BUFFER0) \n\
     <Btn2Down>:         insert-selection(PRIMARY, CUT_BUFFER0) \n\
     <Btn3Down>:         menuCall()\n\
     Shift<Btn3Down>:    menuCall()\n\
     Ctrl<Btn3Down>:     menuCall()\n";

#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */


/* define here all the normal returns on things */

static XtActionsRec regular_actions[] = {
  {"exposeCall", exposeCall},
  {"resizeCall", resizeCall},
  {"mouseCall", mouseCall},
  {"pickCall", pickCall},
  {"middleCall", middleCall},
  {"keyboardCall", keyboardCall},
  {"activate", activeCall},
  {"promptReturn", promptReturn},
  {"propertyCall", propertyCall},
  {"menuCall", menuCall},
  {"enterCall", enterCall},
  {"leaveCall", leaveCall},
  {"returnCall", returnCall}
};

/**********************/

/* X Error handling support.                                                 */
/*                                                                           */
/* errorHandler   - reports X protocol and server errors.                    */
/*                                                                           */

extern void invokeDebugger (void) ;

int supressXError = 0;

/* Handle X Protocol/server errors.                                          */
/* NOTE: you must not issue any calls that will cause X Protocol requests    */
/* i.e. requests to the X server, in this routine.                           */
static int xErrorHandler (Display *display, XErrorEvent *err)
  {
  char msg_code[80] ;
  char msg_major[80] ;
  int result = 0 ;

  if (supressXError)
    {
    /* Application has asked to suppress X error reporting.                  */
    result = 0 ;
    }
#if XtSpecificationRelease < 6
  else if (inSignalRoutine_G == TRUE)
    {
    /* There is always a chance that the signal handler will interrupt an    */
    /* X call, then the server will report a protocol error....this only     */
    /* applies to X11r5 and earlier.                                         */
    fprintf(stderr,
	    "There has been an X Protocol error almost certainly caused"
	    "by you attempting to interrupt the program, you should ask for"
	    "your system to be upgraded to X Windows Version 11, release 6 to"
	    "avoid this problem in the future. You should now exit the application"
	    "to avoid potential further errors.");
    result = 0 ;
    }
#endif
  else
    {
    /* X Protocol/server error.                                              */
    XGetErrorText (display, err->error_code, msg_code, 80) ;
    XGetErrorDatabaseText (display, "GraphPackage", 
			   messprintf ("XRequestMajor.%d",err->request_code), 
			   "undefined", msg_major, 80) ;

    fprintf (stderr, "X error: code %d (%s), request (%d, %d) %s\n", 
	     err->error_code, msg_code, 
	     err->request_code, err->minor_code, msg_major) ;

    result = 0 ;
    }

  invokeDebugger () ;
  return (result) ;
  }

/* Handle errors from Xlib when it has detected an error in a system call,   */
/* often because the network has gone down or the server has died.           */
/* This routine must exit and there is almost certainly no point in trying   */
/* to Xlib calls.                                                            */
static int xIOErrorHandler(Display *display)
  {
  int rc ;

  rc = fprintf(stderr,
	       "Xlib has detected a system error, this may be a network problem,"
	       "the X server may have died, or there may a resource shortage on"
	       "your machine (e.g. no memory or disk space left). The application"
	       "cannot carry on and will now terminate.\n") ;

  invokeDebugger () ;					    /* Probably pointless... */

  exit(EXIT_FAILURE) ;

  return (rc) ;						    /* Never get here.. */
  }

/* Handle minor errors from Xt, OK to carry on working.                      */
static void xtWarningHandler (char *msg)
{
  if (msg && *msg)
    {
      fprintf (stderr, "Xt Warning: %s\n", msg) ;
      fflush (stderr) ;
    }
  
  invokeDebugger () ;
  
  return ;
}

/* Handle drastic errors in Xt (recursion in Form widgets for instance).     */
/* This routine must exit, Xt behaviour will be undefined if we try to carry */
/* on.                                                                       */
static noreturn void xtErrorHandler(String msg)
{
  fprintf(stderr, "Xt Error: %s\n", msg) ;
  
  invokeDebugger () ;
  
  exit(EXIT_FAILURE) ;
}



/* Interrupt signal support.                                                 */
/* Handle user interrupting the application (perhaps via CTRL-C)             */
/* The pre release 6 code is NOT safe, it can cause X protocol errors, we    */
/* attempt to trap this and report it.                                       */
/*                                                                           */
#if XtSpecificationRelease >= 6

static void intProc (int signum)
  {
  XtNoticeSignal(sigintID_G) ;

  return ;
  }

void sigintCallback(XtPointer client_data, XtSignalId *id)
  {

  if (graphQuery ("Do you want to abort?"))
    messExit("User initiated abort") ;

  return ;
  }

#else
/* Pre X11r6 code, set a global so we know if we were in this routine when   */
/* a protocol error occurred, most times we will get away with this, but     */
/* only 'most' times...                                                      */
static void intProc (int signum)
  {
  inSignalRoutine_G = TRUE ;

  if (graphQuery ("Do you want to abort?"))
    messExit("User initiated abort") ;

  inSignalRoutine_G = FALSE ;

  return ;
  }
#endif


static char  idleCall (void *vp)
{
  return gIdle (graphEventX, graphEventY) ? (char)1 : (char)0 ; 
}


void graphDevInit (int *argcptr, char **argv)
{
  uid_t savuid = geteuid(); /* preserve effective uid */
  int i = 0 ;
  static String fallback_resources[] = { "GraphPackage*Scrollbar.height: 14", 
					 "GraphPackage*Scrollbar.width: 14", NULL } ;

#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
  int screen_num ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */


  assw2g = assCreate();
  loopInfo = arrayCreate (4, LOOPINFO) ;

  /*
   * Zip through the arguments and see if the -install option is
   * present.  This indicates that we should install our own colormap
   */
  for(i=1;i<*argcptr;i++) {
        if (!strcmp(argv[i],"-install")) {
                installOwnColorMap = TRUE;
                for(i++;i<*argcptr;i++) {
                        argv[i-1] = argv[i];
                }
                argv[*argcptr-1] = 0;
                (*argcptr)--;
                break;
        }
  }
  i = 0;

  seteuid(getuid()); 
     /* set the effective uid to be the real uid here, so that 
        X access control can read the .Xauthority file of the
        user. Note that this executes before the user id munging
        in session.c, we restore the status quo ante so as not
        to upset that code. */


/*   XtSetArg(args[i], XtNmappedWhenManaged, False); i++; */
  root_widget = XtAppInitialize(&app_con, "GraphPackage",
                                NULL, ZERO, argcptr, argv, 
                                fallback_resources,
                                args, i) ;
/* creates an unrealized top widget to be parent for the rest */



#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
  /* Get the true size of the screen instead of guessing as before...        */
  /* Commented out because it messes up window sizes completely...           */
  screen_num = DefaultScreen(XtDisplay(root_widget)) ;
  screenx = DisplayWidth(XtDisplay(root_widget), screen_num) ;
  screeny= DisplayHeight(XtDisplay(root_widget), screen_num) ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */



  if (getenv("ACEDB_X_SYNCHRO"))
    XSynchronize (XtDisplay(root_widget), True) ; 


  /* Set up X error handlers, the default for all but xtWarningHandler is    */
  /* to issue a message and exit the application.                            */
  XSetErrorHandler(xErrorHandler) ;			    /* Handle X Protocol/server errors. */
  XSetIOErrorHandler(xIOErrorHandler) ;			    /* Handle system errors detected by Xlib. */
  XtAppSetWarningHandler(app_con, xtWarningHandler) ;	    /* Handle minor problems in Xt. */
  XtAppSetErrorHandler(app_con, xtErrorHandler) ;	    /* Handle major problems in Xt. */



  cmapDebug = (getenv("ACEDB_CMAP_DEBUG") != NULL);
  if (getenv("ACEDB_CMAP_VISUAL") != NULL)
    {
    /* They have suggested a visual that they want to use */
    sscanf(getenv("ACEDB_CMAP_VISUAL"),"%d",&preferredVisual);
    }
  XtVaSetValues (root_widget, XtNcolormap, getGlobalColorMap(), NULL) ;


  XtAppAddActions(app_con, regular_actions, XtNumber(regular_actions) );
/*  XawSimpleMenuAddGlobalActions(app_con) ; */


  /* Interrupt handling.                                                     */
  /* The X11r6 code is safe, the other is not.                               */
#if XtSpecificationRelease >= 6
  /* Register the signal callback and get the signal id.                     */
  sigintID_G = XtAppAddSignal(app_con, sigintCallback, NULL) ;
#endif
  
  /* Now register the signal handler itself.                                 */
  signal (SIGINT, intProc) ;


  /* Initialise the Xmu atom caching needed for cut/paste support, and set   */
  /* up atoms used in cut/paste that are not part of the standard X protocol.*/
  XmuInternAtom(XtDisplay(root_widget), XmuMakeAtom("NULL")) ;
  XA_TARGETS(XtDisplay(root_widget)) ;
  XA_COMPOUND_TEXT(XtDisplay(root_widget)) ;
  XA_TEXT(XtDisplay(root_widget)) ;


  /* mieg: i reverse the logic */
  noGrab = (getenv("ACEDB_X_NOGRAB") != NULL);
  doGrab = (getenv("ACEDB_X_GRAB") != NULL);

  if (!doGrab)  /* mieg: i reverse the logic */
    XtRegisterGrabAction (menuCall, True,
                          (unsigned)(ButtonPressMask | ButtonReleaseMask),
                          GrabModeAsync, GrabModeAsync) ;
  lastLeft = assCreate();
  destroyed = assCreate();
  seteuid(savuid); /* restore */

  return ;
}



void graphDevFinish(void)
  {
  if (gWriteColors) XStoreBuffer(XtDisplay(root_widget), 
				 (char*)&gWriteColors, 1, 5) ;

  /* ACEDB-GRAPH INTERFACE: acedb needs to stop communication with other     */
  /* programs using X specific methods.                                      */
  if (getGraphAcedbWcommFinish() != NULL) (getGraphAcedbWcommFinish())() ;

  XFreePixmap(XtDisplay(root_widget), GreyPixmap);
  XtDestroyWidget(root_widget) ;
  assDestroy(lastLeft);
  }

/********** creation ************/

/* comment for when I interface into the graph package:
   I could use the instance name as designated by the
   type--a type2string or something */

static BOOL sneakyTransient = FALSE ;

static XtGeometryResult GeometryManager(Widget w, 
                                        XtWidgetGeometry *request,
                                        XtWidgetGeometry *reply)
/* 
   This is a simple geometry manager to be used for the Composite
   widget.  It will allow a child widget to arbitrarily size and
   position itself inside its parent as long as it fits within
   parent's boundaries and allow any border width 
*/
{
  Widget self = w->core.parent;
  Position x = w->core.x;
  Position y = w->core.y;
  Dimension width = w->core.width;
  Dimension height = w->core.height;
  Dimension border_width = w->core.border_width;


  /* Viewport Widget is not allowed to resize so that it won't try to 
     accomodate the size of its child when one of the scrollbars is turned 
     off.  This is similar to the way Viewports are handled by Shell widget */

  if (XtClass(w)==viewportWidgetClass && 
      ((request->request_mode & CWWidth)||(request->request_mode & CWHeight)))
    return XtGeometryNo;

  if ((request->request_mode & CWWidth) && 
      x+request->width <= self->core.width &&
      request->width > 1 /* 2007_02_13 */
      )
    width = request->width;
  
  if ((request->request_mode & CWHeight) && 
      y+request->height <= self->core.height &&
      request->height > 1 /* 2007_02_13 */
    )
    height = request->height;
  
  if ((request->request_mode & CWX) && 
      request->x <= self->core.width - width &&
      request->x >= 0  /* 2007_02_13 */
      )
    x = request->x;
  
  if ((request->request_mode & CWY) && 
      request->y <= self->core.height - height &&
      request->y >= 0  /* 2007_02_13 */
      )
    y = request->y;

  if (request->request_mode & CWBorderWidth)
    border_width = request->border_width;

  XtConfigureWidget(w, x, y, width, height, border_width);

  return XtGeometryDone;
}

static SubDev YsubCreate (Graph_ graph, Widget parent,
                          int x, int y, int w, int h)
/* YsubCreate creates a subgraph of specified type */
{
  int i;
  SubDev subdev;

  subdev = (SubDev) messalloc (SUBDEV_SIZE);

  i=0;
  XtSetArg(args[i], XtNx, x); i++;
  XtSetArg(args[i], XtNy, y); i++;
  XtSetArg(args[i], XtNheight, h ); i++;
  XtSetArg(args[i], XtNwidth, w ); i++;
  XtSetArg(args[i], XtNforceBars, TRUE ); i++;
  XtSetArg(args[i], XtNborderWidth, 2); i++;
  XtSetArg(args[i], XtNborderPixmap, GreyPixmap); i++;
  XtSetArg(args[i], XtNcolormap, getGlobalColorMap()); i++;
  
  switch (graph->type)
    {
    case TEXT_FULL_EDIT:    /* scrollbars are managed by the widget */
    case PLAIN: case TEXT_FIT: case PIXEL_FIT:
      break ;
    case TEXT_SCROLL:
      XtSetArg(args[i], XtNallowVert, TRUE); i++ ;
      break ;
    case MAP_SCROLL:
      XtSetArg(args[i], XtNallowHoriz, TRUE); i++ ;
      break ;
    case PIXEL_SCROLL:
    case TEXT_FULL_SCROLL:
      XtSetArg(args[i], XtNallowVert, TRUE); i++ ;
      XtSetArg(args[i], XtNallowHoriz, TRUE); i++ ;
      break ;
    }

  subdev->viewport = XtCreateManagedWidget("viewport", 
                                           viewportWidgetClass,
                                           parent, args, i);

  XtAugmentTranslations(subdev->viewport,
			XtParseTranslationTable(viewportTrans));

  /* You have to use SubstructureNotifyMask to get a map call for ALL        */
  /* subgraphs...this involves some filtering of events in subDevMapCall     */
  XtAddEventHandler(subdev->viewport, SubstructureNotifyMask, False, subDevMapCall,
		    (XtPointer)graph) ;


  i=0;
  XtSetArg(args[i], XtNheight, h ); i++;
  XtSetArg(args[i], XtNwidth, w ); i++;
  XtSetArg(args[i], XtNcolormap, getGlobalColorMap()); i++;
  
  if (graph->type == TEXT_FULL_EDIT)
    {
      /* auto-fill is off by default, but can be turned on from menu */
      XtSetArg(args[i], XtNautoFill, True ); i++; 
      XtSetArg(args[i], XtNeditType, XawtextEdit ); i++;
      XtSetArg(args[i], XtNscrollVertical, XawtextScrollAlways ); i++;
      XtSetArg(args[i], XtNwrap, XawtextWrapLine ); i++;
      XtSetArg(args[i], XtNfont, getDefaultFont() ); i++;

      XtSetArg(args[i], XtNtranslations,
	       XtParseTranslationTable(textFullActionTrans) ); i++;
      subdev->base = XtCreateManagedWidget("full_text", asciiTextWidgetClass,
					   subdev->viewport, args, i);
    }
  else
    {
      XtSetArg(args[i], XtNtranslations,
	       XtParseTranslationTable(actionTrans) ); i++;
      subdev->base = XtCreateManagedWidget("composite", compositeWidgetClass,
					   subdev->viewport, args, i);
      /* Install the geometry manager */
      ((CompositeWidgetClass) (subdev->base->core.widget_class))
	->composite_class.geometry_manager = GeometryManager;
    }

  return subdev ;
}

static Dev YCreate (Graph_ graph, int x, int y, int w, int h)
/* YCreate creates a popup shell and puts a subgraph inside */
{
  int i;
  Dev dev;

  dev = (Dev) messalloc (DEV_SIZE)  ;

  i=0;
  XtSetArg(args[i], XtNmappedWhenManaged, True); i++;
  XtSetArg(args[i], XtNtitle, graph->name); i++;
  XtSetArg(args[i], XtNiconName, graph->name); i++;
  XtSetArg(args[i], XtNx, x ); i++;
  XtSetArg(args[i], XtNy, y ); i++;
  XtSetArg(args[i], XtNcolormap, getGlobalColorMap()); i++;

  dev->popup = XtCreatePopupShell(typeName[graph->type],
	               sneakyTransient ? overrideShellWidgetClass
				       : topLevelShellWidgetClass,
				  root_widget, args, i);

  XtAddCallback(dev->popup, 
		XtNdestroyCallback, destroyCallback, "hi folks, I'm dead");
  
  dev->subdev=YsubCreate(graph, dev->popup, 0, 0, w, h);

  return dev ;
}

void graphDevCreate (Graph_ graph, float x, float y, float w, float h)
{
  int ix = x*screenx ;
  int iy = y*screeny ;
  int iw = w*screenx ;
  int ih = h*screeny ;
  static BOOL isFirst = TRUE ;

  graph->dev = YCreate (graph, ix, iy, iw, ih) ;
  graph->subdev = graph->dev->subdev;

/*  ?no longer needed -- we used to have this to stop overdrawing when drawing a window
  graph->dev->isExposed = FALSE ;
*/

  assRemove(destroyed, graph->dev->popup);
  assInsert(assw2g, graph->dev->popup, graph);
  assInsert(assw2g, graph->dev->subdev->base, graph);
  assInsert(assw2g, graph->dev->subdev->viewport, graph);

  graph->w = iw;
  graph->h = ih;
				/* put window on screen */
  XtRealizeWidget (graph->dev->popup) ;
  XtPopup (graph->dev->popup, XtGrabNone) ;

  fakeResize (graph) ;		/* set up size state in graph structure */

  if (isFirst)
    {
    /* ACEDB-GRAPH INTERFACE: acedb needs to set up communication with other */
    /* programs using X specific methods.                                    */
    if (getGraphAcedbWcommInit() != NULL)
      (getGraphAcedbWcommInit())(XtDisplay(graph->dev->popup), /* not in graphDevInit cos */
				 XtWindow(graph->dev->popup)) ;	/* needs real window */

    XSetWindowColormap(XtDisplay(graph->dev->popup),
		       XtWindow(graph->dev->popup),getGlobalColorMap());

    isFirst = FALSE ;
    }

}  

void graphSubDevCreate (Graph_ graph, float x, float y, float w, float h)
{
  graph->subdev = YsubCreate (graph, graph->parent->subdev->base, 
			      x, y, w, h) ;
  
  assInsert(assw2g, graph->subdev->base, graph);
  assInsert(assw2g, graph->subdev->viewport, graph);
  assRemove(destroyed, graph->subdev->viewport); /* may re-use memory */
  assRemove(destroyed, graph->subdev->base); /* may re-use memory */

  graph->w = w;
  graph->h = h;
				/* put window on screen */
  XtRealizeWidget (graph->subdev->viewport) ;

  fakeResize (graph) ;		/* set up size state in graph structure */
}  

void graphDevDestroy (void)
{
  Dev dev = gDev ;

  while (graphLoopReturn (0)) ;

  if (assRemove (assw2g,dev->popup)) /* i.e. not looping via destroyCallback */
    { 
      assInsert(destroyed, dev->popup, 0);
      XtPopdown(dev->popup) ;
      XtDestroyWidget (dev->popup) ;
    }

  gActive->dev = gDev = 0;
  messfree(dev);
}

void graphSubDevDestroy (void)
{
  SubDev subdev = gSubDev;
  char*  data ;
  XImage *im ;
  
  while (graphLoopReturn (0)) ;
  
  if (assRemove (assw2g, subdev->viewport))
    { assInsert(destroyed, subdev->viewport, 0);
      XtDestroyWidget (subdev->viewport) ;
    }

  assRemove (assw2g, subdev->base) ;
  assRemove (assw2g, subdev->viewport) ;
  assRemove (lastLeft, subdev->base) ;
  assRemove (lastLeft, subdev->viewport) ;
  
  if (subdev->images)
    { data = 0 ; im = 0 ;
      while (assNext (subdev->images, &data, &im))
	XDestroyImage (im) ;
      assDestroy (subdev->images) ;
    }
  assDestroy (subdev->boxMenus) ;
  
  gActive->subdev = gSubDev = 0;
  messfree(subdev);
}

void graphMapBounds (float ux, float uy, float uw, float uh, float aspect)
{
  if (gActive->type != MAP_SCROLL)
    messcrash ("MapBounds called on invalid graph type %d",
	       gActive->type) ;

  gActive->ux = ux ;
  gActive->uy = uy ;
  gActive->uw = uw ;
  gActive->uh = uh ;
  gActive->aspect = aspect ;

  if (gDev && gSubDev)		/* reset width */
    { int i = 0 ;
      gActive->w = gActive->h * gActive->uw / 
		      (gActive->uh * gActive->aspect * gPixAspect) ;
      XtSetArg (args[i], XtNwidth, gActive->w ); i++;
      XtSetValues (gSubDev->base, args, i) ;
      graphFacMake () ;
    }
}

void graphScreenSize (float *sx, float *sy, float *fx, float *fy, int *px, int *py)
{
  int fw, fh ;  /* Font width and height */
  Screen *screen = XtScreen (root_widget) ;

  gFontInfo (0, &fw, &fh) ;

  if (sx) *sx = (float)screen->width/screenx ;
  if (sy) *sy = (float)screen->height/screeny ;

  if (fx) *fx = (float)screen->width/fw ;
  if (fy) *fy = (float)screen->height/fh ;

  if (px) *px = screen->width ;
  if (py) *py = screen->height ;
}

/* those values should be used to recreate the same window */
BOOL graphWindowSize (float *wx, float *wy, float *ww, float *wh)
{
  if (gSubDev && screenx && screeny) 
    {
      /* cast explicitelly since i do not know what a Dimension is */
      Dimension w, h ; Position x, y ;
      XtVaGetValues (gSubDev->viewport, XtNheight, &h, XtNwidth, &w, 
		     XtNx, &x, XtNy, &y, NULL) ;
      if (wx) *wx = x/screenx ; 
      if (ww) *ww = w/screenx ; 
      if (wy) *wy = y/screeny ; 
      if (wh) *wh = h/screeny ; 
      return TRUE ;
    }
  return FALSE ;
}


void graphTextBounds (int nx, int ny)
{
  if (gActive->type != TEXT_SCROLL && gActive->type != TEXT_FULL_SCROLL)
    messcrash ("textBounds called on invalid graph type %d",
	       gActive->type) ;

  gActive->uw = nx ;
  gActive->uh = ny ;

  if (gDev && gSubDev)		/* reset height */
    { int i = 0 ;
      gActive->h = ny * gActive->yFac ;
      XtSetArg (args[i], XtNheight, gActive->h ); i++;
      if (gActive->type == TEXT_FULL_SCROLL)
	{ gActive->w = nx * gActive->xFac ;
	  XtSetArg (args[i], XtNwidth, gActive->w ); i++;
	}
      XtSetValues (gSubDev->base, args, i) ;
    }
}

/* fake modif, used in forestdisp.c to prepare a longer printout
 * and restore old at the end
 */
float graphFakeBounds (float ny)
{
  float old = gActive->uh ;
  gActive->uh = ny ;

  if (gDev && gSubDev)	
    { int i = 0 ;
      gActive->h = ny * gActive->yFac ;
      XtSetArg (args[i], XtNheight, gActive->h ); i++;
      XtSetValues (gSubDev->base, args, i) ;
    }
  return old ;
}

void graphPixelBounds (int nx, int ny)
{
  if (gActive->type != PIXEL_SCROLL)
    messcrash ("pixelBounds called on invalid graph type %d",
	       gActive->type) ;

  gActive->uw = nx ;
  gActive->uh = ny ;

  if (gDev && gSubDev)
    { int i = 0 ;
      gActive->h = ny ;
      gActive->w = nx ;
      XtSetArg (args[i], XtNwidth, gActive->w ); i++;
      XtSetArg (args[i], XtNheight, gActive->h ); i++;
      XtSetValues (gSubDev->base, args, i) ;
    }
}

/*******************************/

static Graph_ blockGraph = 0 ;
static Widget wMess = 0 ;
static int loopLevel = 0 ;

static void destroyCallback (Widget w, XtPointer client_data, XtPointer call_data)
{
  if (activeSet(w))
    { assRemove (assw2g,gDev->popup) ;
      gDev->popup = 0 ;
      graphDestroy () ;
    }
}

/******* process all outstanding events in the queue *******/

void graphProcessEvents (void)
{
  XEvent ev ;

  while (XtAppPending (app_con))
    { XtAppNextEvent (app_con, &ev) ;
      XtDispatchEvent (&ev) ;
    }
}

/******* main loop call *************************************************/
/******* can now be recursive and either blocking or non-blocking *******/

	  
int graphLoop (BOOL isBlock)
{
  LOOPINFO *info ;
  XEvent ev ;
  int retval ;

  if (!gDev || !gSubDev)
    return 0 ;

  info = arrayp(loopInfo, ++loopLevel, LOOPINFO) ;
  info->g = gActive ;
  info->isBlocking = isBlock ;
  info->retval = 0 ;
  info->isDone = FALSE ;
  info->reportActivity = getGraphAcedbMainActivity()  ;

  if (isBlock)
    blockGraph = gActive ;

  while (!info->isDone)
    {
    XtAppNextEvent (app_con, &ev) ;
				/* OpenWindows exit - from Suzi */
    if (ev.type == ClientMessage &&
	ev.xclient.data.l[0] == myExitAtom) /* defined in graphxlib.c */
      {
      activeSet (XtWindowToWidget(XtDisplay(gSubDev->base), 
				  ev.xclient.window)) ;
      graphDestroy ();  /* WM wants a window closed... */
      }
    else
      { 
	XtDispatchEvent (&ev) ;
      }
    XtAppAddWorkProc(app_con, idleCall, NULL) ;

    if(!isBlock && !wMess && info->reportActivity)
      (info->reportActivity) (NULL) ;  
    }
  retval = info->retval ;

  info = arrp(loopInfo, --loopLevel, LOOPINFO) ;
  if (info->isBlocking)
    blockGraph = info->g ;
  else
    blockGraph = 0 ;

  return retval ;
}


#ifdef DEBUG
static void loopStatus (void)
{
  int i ;
  LOOPINFO *loop = arrp(loopInfo, 0, LOOPINFO) ;

  for (i = 1, ++loop ; i <= loopLevel ; ++i, ++loop)
    printf (" level %d, isDone %d, isBlocking %d, ret %d, graph %x\n",
	    i, loop->isDone, loop->isBlocking, loop->retval, loop->g) ;
}
#endif

BOOL graphLoopReturn (int retval)
{
  int i = loopLevel ;

#ifdef DEBUG
  printf ("LoopReturn %d on graph %x\n", retval, gActive) ;
  loopStatus() ;
#endif

  while (i)
    if (arrp(loopInfo, i, LOOPINFO)->g == gActive &&
	!arrp(loopInfo, i, LOOPINFO)->isDone)
      break ;
    else
      --i ;
  if (i)
    { arrp(loopInfo, i, LOOPINFO)->retval = retval ;
      while (i <= loopLevel)
	arrp(loopInfo, i++, LOOPINFO)->isDone = TRUE ;
#ifdef DEBUG
      printf ("Return TRUE\n") ; loopStatus() ;
#endif
      return TRUE ;
    }
#ifdef DEBUG
  printf ("Return FALSE\n") ; loopStatus() ;
#endif
  return FALSE ;
}

/*****************************************************************/

void graphGoto (float x, float y)  /* tries to center x,y in visible region */
{ 
  Widget sbar ;
  float shown ;

  if (!gDev || !gSubDev)
    return ;
  if ((sbar = XtNameToWidget (gSubDev->viewport,"horizontal")))
    {
      XtVaGetValues (sbar, XtNshown, &shown, NULL) ;

      x = (x - gActive->ux)/gActive->uw - 0.5*shown ;
      if (x < 0) x = 0 ;
      if (x > 1 - shown) x = 1 - shown ;

      shown = x ;
      XtCallCallbacks (sbar, XtNjumpProc, &shown) ;
    }
  if ((sbar = XtNameToWidget (gSubDev->viewport,"vertical")))
    {
      XtVaGetValues (sbar, XtNshown, &shown, NULL) ;

      y = (y - gActive->uy)/gActive->uh - 0.5*shown ;
      if (y < 0) y = 0 ;
      if (y > 1 - shown) y = 1 - shown ;

      shown = y ;
      XtCallCallbacks (sbar, XtNjumpProc, &shown) ;
    }
}

void graphWhere (float *x1, float *y1, float *x2, float *y2)
{
  Widget sbar ;
  Dimension dim ;
  float pos, shown ;

  if (!gDev) {
	/* DAZ - shouldn't quietly return with no error code */
	messcrash("call to graphWhere before gDev is available");
  }
  if ((sbar = XtNameToWidget (gSubDev->viewport,"horizontal")))
    { XtVaGetValues (sbar, XtNtopOfThumb, &pos, 
		     XtNshown, &shown, NULL) ;
      if (x1) *x1 = gActive->ux + pos * gActive->uw ;
      if (x2) *x2 = gActive->ux + (pos + shown) * gActive->uw ;
    }
  else
    { XtVaGetValues (gSubDev->viewport, XtNwidth, &dim, NULL) ;
      if (x1) *x1 = gActive->ux ;
      if (x2) *x2 = gActive->ux + dim / gActive->xFac ;
    }
  if ((sbar = XtNameToWidget (gSubDev->viewport,"vertical")))
    { XtVaGetValues (sbar, XtNtopOfThumb, &pos, 
		     XtNshown, &shown, NULL) ;
      if (y1) *y1 = gActive->uy + pos * gActive->uh ;
      if (y2) *y2 = gActive->uy + (pos + shown) * gActive->uh ;
    }
  else
    { XtVaGetValues (gSubDev->viewport, XtNheight, &dim, NULL) ;
      if (y1) *y1 = gActive->uy ;
      if (y2) *y2 = gActive->uy + dim / gActive->yFac ;
    }
}

void graphEvent (int action, float x, float y)
/* sends an event over the "wire" to active window
   currently can only handle printable ascii events
   and mouse events */
{
  XEvent *event = 0 ;
  Widget w ;
  static char buf[2] ;
   int tmp ;
  
  if (!gDev || !gSubDev)
    return ;
  w = gSubDev->base ;

  if(action < 0 || action > 127) 
    { messcrash ("graphEvent() can only handle actions between 0 and \
                127\n");
      return ;
    }
  
  if (action <= RIGHT_UP) /* means a mouse event */ {
    if (action % 3 == 1 ) 
      { fprintf (stderr, "Sorry, no drag events yet \n") ;
	return ;
      }

    /* Button up or down */
    event = (XEvent*) messalloc(sizeof (XButtonEvent) ) ;
    event->xbutton.type = ((tmp = action % 3) == 2 ? ButtonRelease :
			   ButtonPress);
    event->xbutton.button = ((tmp =  action / 3 )== 0 ? Button1 :
			     (tmp == 1 ? Button2 : Button3 )
			     ) ;
#ifdef DEBUG
    fprintf(stderr, "Button# = %d , action# = %d \n", tmp+1, action) ;
#endif
     event->xbutton.x = uToXabs(x) ;
     event->xbutton.y = uToYabs(y) ;
   }
 
  /* keyboard events */
  if(action > RIGHT_UP ) {
    
    event = (XEvent*) messalloc(sizeof (XKeyEvent)) ;
    event->xkey.type = KeyPress ;
    
    buf[0] = action ;
    buf[1] = 0 ;
    
    event->xkey.keycode = XKeysymToKeycode(XtDisplay(w),XStringToKeysym(buf)) ;
    /* KeyPressMask ; */
    event->xkey.x = uToXabs(x) ;
    event->xkey.y = uToYabs(y) ;
    
  }

  /* information common to all events */

  event->xany.send_event = TRUE ;
  event->xany.display = XtDisplay(w) ;
  event->xany.window = XtWindow(gSubDev->base) ;

  XSendEvent (XtDisplay(w), XtWindow(w), 
	     FALSE, KeyPressMask,		       
	     event) ;
#ifdef DEBUG
  fprintf(stderr, "ID %d, status %d, event =%d\n", XtDisplay(w),s, (int ) event);
#endif
  messfree(event) ;
}

void graphRetitle (const char *title)
{
  if (!gDev || !gSubDev) return ;
  
  XtSetArg(args[0], XtNtitle, title) ;
  XtSetArg(args[1], XtNiconName, title);
  XtSetValues(gDev->popup, args, 2) ;
}

/************/

void graphPop (void)
{
  if (gDev && gSubDev)
    { 
/*
      if (getenv("ACEDB_SCREEN_FIX"))
	XIconifyWindow (XtDisplay(gDev->popup),XtWindow(gDev->popup), 0) ;

      XtVaSetValues (gDev->popup, XtNiconic, FALSE, NULL) ;
*/
      XRaiseWindow (XtDisplay(gDev->popup),XtWindow(gDev->popup)) ;
      XMapWindow(XtDisplay(gDev->popup),XtWindow(gDev->popup)) ;
    }
}

/**********/

static void killMessage (Widget w, XtPointer client_data, XtPointer call_data)
{
  SubDev subdev = (SubDev) client_data ;

  XtPopdown(subdev->message) ;
  XtDestroyWidget (subdev->message) ;
  subdev->message = 0 ;
}

static void messageDestroy (Widget w, XtPointer client_data, 
			    XtPointer call_data)
{
  Graph_ g = (Graph_) client_data ;
  Graph old = graphActive () ;

  wMess = 0 ; /* at this point, gSubDev->message is already zero */

  if (graphExists (g->id) && 
      g->func[MESSAGE_DESTROY] && 
      graphActivate (g->id))
    { (*g->func[MESSAGE_DESTROY])() ;
      graphActivate (old) ;
    }
}

void graphMessage (char *text)
{
  Widget dialog, button ;
  Position x,y ;
  int n ;

  if (!gDev || !gSubDev || gSubDev->message)
    return ;

  XtTranslateCoords (gSubDev->viewport, 
		     (Position) 25, (Position) 20,
		     &x, &y) ;
  n = 0 ;
  XtSetArg (args[n], XtNx, x) ;	n++ ;
  XtSetArg (args[n], XtNy, y) ;	n++ ;
  XtSetArg(args[n], XtNcolormap, getGlobalColorMap()); n++;
  gSubDev->message = XtCreatePopupShell("message", 
					transientShellWidgetClass,
					gDev->popup, args, n) ;
  dialog = XtCreateManagedWidget ("textmessage", dialogWidgetClass,
				  gSubDev->message, NULL,0) ;
  XtVaSetValues (dialog, XtNlabel, uBrokenText (text, 40), NULL) ;
  XtVaSetValues (dialog, XtNcolormap, getGlobalColorMap(), NULL) ;
  button = XtCreateManagedWidget ("Remove", commandWidgetClass, 
				  dialog, args, 0) ;
  XtAddCallback (button, XtNcallback, killMessage, (XtPointer)gSubDev) ;
  XtAddCallback (gSubDev->message, XtNdestroyCallback, 
		 messageDestroy, (XtPointer)gActive) ;
  XtPopup (gSubDev->message,XtGrabNone) ;
  wMess = gSubDev->message ;
}

void graphUnMessage (void) 
{ 
  if (gDev && gSubDev && gSubDev->message)
    { 
      XtPopdown(gSubDev->message) ;
      XtDestroyWidget (gSubDev->message) ;
      gSubDev->message = 0 ;
    }
}

/*****************************************************************/

/* callbacks, associators and all that */

/* The following take care of the callback and action routine
 * interface to the graph package
 * 1. install/interface the graphs callbacks on the widgets actions
 */

float graphEventX, graphEventY ;
	/* for users to get at event position, in user coords */

static BOOL activeSet(Widget w)
{
  Graph_ g ;

  if (assFind(assw2g, w, &g))
    { graphActivate (g->id) ;
      return TRUE ;
    }
  else
    return FALSE ;
}

/*****************************************************************/

static Bool exposeCheck (Display* display, XEvent* ev0, char* arg)
{
  Window win = (Window) arg ;
  XAnyEvent* ev = (XAnyEvent*)ev0 ;

  return (ev->type == Expose) && (ev->window == win) ;
}

static void exposeCall (Widget w, XEvent* ev0, String *params, Cardinal *num_params)
{
  XEvent  ev1 ;
  XExposeEvent* exp0 = (XExposeEvent*)ev0 ;
  XExposeEvent* exp1 = (XExposeEvent*)&ev1 ;
  int xmin,xmax,ymin,ymax ;
  Graph gsave = gActive->id ;

  xmin = exp0->x ; xmax = xmin + exp0->width ;
  ymin = exp0->y ; ymax = ymin + exp0->height ;

  while (XCheckIfEvent (exp0->display, 
			&ev1, exposeCheck, 
			(char*)(exp0->window)))
    { if (xmin > exp1->x) 
	xmin = exp1->x ;
      if (xmax < exp1->x + exp1->width) 
	xmax = exp1->x + exp1->width ;
      if (ymin > exp1->y) 
	ymin = exp1->y ;
      if (ymax < exp1->y + exp1->height) 
	ymax = exp1->y + exp1->height ;
    }

  if (activeSet (w))
    { graphClipDraw (xmin,ymin,xmax,ymax) ;
      graphActivate (gsave) ;
    }
}





/* This routine calls the graph package every time a graph/subgraph          */
/* gets mapped, we do this by monitoring map events on the viewport widget   */
/* which is the 'parent' widget for all subgraphs of a graph.                */
static void subDevMapCall (Widget w, XtPointer client_data, XEvent* event, Boolean *propagate)
{
  /* We have to filter here for MapNotify events only, and then we need to   */
  /* filter for events that are for the viewport widgets and not their child */
  /* widgets.                                                                */ 
  if (event->type == MapNotify)
    {
      if (!activeSet (w))
	return ;

      /* When graphCreate is called the busy cursors are turned on, this is      */
      /* where they are turned off, i.e. once the window is shown on the screen. */
      if (gActive->window_being_constructed == TRUE)
	{
	  graphBusyCursorAll(FALSE) ;
	  gActive->window_being_constructed = FALSE ;
	}

    }
  
  return ;
}


static void resizeCall (Widget w, XEvent* ev, String *params, Cardinal *num_params)
{
  int width = ev->xconfigure.width ;
  int height = ev->xconfigure.height ;
  int i, isRedraw, isResize = FALSE ;
  Widget sbar ;

  if (!activeSet (w))
    return ;
  
  switch (gActive->type)
    {
    case PLAIN :
      isRedraw = (width <= gActive->w && height <= gActive->h &&
		  (width < gActive->w || height < gActive->h)) ;
      gActive->w = width ;
      gActive->h = height ;
      graphFacMake () ;
      if (isRedraw)
	graphRedraw () ;
      break ;
    case TEXT_FIT: case PIXEL_FIT:
      if (gActive->w != width || gActive->h != height)
	isResize = TRUE ;
      gActive->w = width ;
      gActive->h = height ;
      gActive->uw = width / gActive->xFac ;
      gActive->uh = height / gActive->yFac ;
      break ;
    case TEXT_SCROLL:
      if (gActive->w != width)
	isResize = TRUE ;
      gActive->w = width ;
      i=0;				
      XtSetArg (args[i], XtNheight, gActive->h ); i++;
      XtSetValues (w, args, i) ;
      break ;
    case PIXEL_SCROLL:
    case TEXT_FULL_EDIT:
    case TEXT_FULL_SCROLL:
      i=0;				
      XtSetArg (args[i], XtNwidth, gActive->w ); i++;
      XtSetArg (args[i], XtNheight, gActive->h ); i++;
      XtSetValues (w, args, i) ;
      break ;
    case MAP_SCROLL:
      if (height == gActive->h)	/* prevents recursion */
	break ;
      isResize = TRUE ;
      if ((sbar = XtNameToWidget (gSubDev->viewport,"horizontal")))
	{ float xpos,shown ;

	  XtVaGetValues (sbar, 
			 XtNtopOfThumb, &xpos, 
			 XtNshown, &shown,
			 NULL) ;
	  xpos *= height / (float) gActive->h ;
	  if (xpos > 1-shown)
	    xpos = 1-shown ;
	  XtCallCallbacks (sbar, XtNjumpProc, &xpos) ;
	}
      gActive->h = height ;
      gActive->w = gActive->h * gActive->uw / 
		      (gActive->uh * gActive->aspect * gPixAspect) ;
      i=0;				
      XtSetArg (args[i], XtNwidth, gActive->w ); i++;
      XtSetValues (w, args, i) ;
      graphFacMake () ;
      break ;
    }

  /* Application may have registered a callback routine for resizes, so set  */
  /* busy cursors and call it.                                               */
  if (isResize && gActive->func[RESIZE])
    {
    graphBusyCursorAll (TRUE);

    (*gActive->func[RESIZE])() ;

    graphBusyCursorAll (FALSE) ;
    }

  return ;
}

static void fakeResize (Graph_ g)
{
  XConfigureEvent ev ;
  ev.height = g->h ;
  ev.width = g->w ;
  resizeCall (g->subdev->base, (XEvent*)(&ev), 0, 0) ;
}

static void blockingRaise (void)
{  
  Widget w = blockGraph->dev->popup ;
  XRaiseWindow (XtDisplay(w), XtWindow(w)) ;
  XMapWindow(XtDisplay(w), XtWindow(w)) ;
}

static void blockingRaiseMess (void) /*mhmp 19.11.97 */
{  
  Widget w = wMess ;
  XRaiseWindow (XtDisplay(w), XtWindow(w)) ;
  XMapWindow(XtDisplay(w), XtWindow(w)) ;
}



/* Note that the Xevent will be one of MotionNotify or ButtonRelease, so it's */
/* important to check this (see translation tables at top of this file).      */
/*                                                                            */
static void mouseCall (Widget w, XEvent* ev, String *params, Cardinal *num_params)
{
  int myEv ;

  

  if (!activeSet (w))
    return ;

  if (gActive->block_mode == GRAPH_BLOCKABLE)
    if (blockGraph && !graphContainedIn(blockGraph->id))
      { blockingRaise () ; return ; }

  if (wMess) /*mhmp 19.11.97 */
    blockingRaiseMess () ;
  graphEventX = XtoUabs(ev->xbutton.x) ;
  graphEventY = YtoUabs(ev->xbutton.y) ;

  if (!sscanf (*params,"%d",&myEv) || myEv < 1 || myEv >= 9)
    messout ("Bad mouse event %s",*params) ;
  else if (gActive->func[myEv] )
    {
    /* Application has registered a callback for mouse clicks, so set busy   */
    /* cursor and call their routine.                                        */
    /* We are only doing this for non-motion events, otherwise it gets a bit */
    /* irritating.                                                           */
    if (ev->type != MotionNotify) graphBusyCursorAll (TRUE);

    (*gActive->func[myEv])(graphEventX, graphEventY); 

    if (ev->type != MotionNotify) graphBusyCursorAll(FALSE) ;
    }

}

static void pickCall (Widget w, XEvent* ev, String *params, Cardinal *num_params)
{
  if (!activeSet (w))
    return ;

  if (gActive->block_mode == GRAPH_BLOCKABLE)
    if (blockGraph && !graphContainedIn(blockGraph->id))
      { blockingRaise () ; return ; }

  if (wMess) /*mhmp 19.11.97 */
    blockingRaiseMess () ;

  graphEventX = XtoUabs(ev->xbutton.x) ;
  graphEventY = YtoUabs(ev->xbutton.y) ;

  gLeftDown (graphEventX, graphEventY) ; 
}

static void middleCall (Widget w, XEvent* ev, String *params, Cardinal *num_params)
{
  if (!activeSet (w))
    return ;

  if (gActive->block_mode == GRAPH_BLOCKABLE)
    if (blockGraph && !graphContainedIn(blockGraph->id))
     { blockingRaise () ; return ; }  

if (wMess) /*mhmp 19.11.97 */
    blockingRaiseMess () ;

  graphEventX = XtoUabs(ev->xbutton.x) ;
  graphEventY = YtoUabs(ev->xbutton.y) ;

  gMiddleDown (graphEventX, graphEventY) ; 
}

static void menuCall (Widget w, XEvent* ev, String *params, 
		      Cardinal *num_params)
{
  Box box = 0 ;
  Widget menu = 0 ;
  static Cursor RightArrow = 0 ;

  if (!RightArrow) 
    RightArrow =
      XCreateFontCursor(XtDisplay(root_widget), XC_draft_large);

#ifdef DEBUG
  printf ("menuCall(): widget %x (position %d, %d).\n",
	   w, ev->xbutton.x, ev->xbutton.y) ;
#endif
  if (!activeSet (w))
      return ;       
  /* only check for blocked graphs if they can be blocked */
  if (gActive->block_mode == GRAPH_BLOCKABLE)
    if (blockGraph && !graphContainedIn(blockGraph->id))
      { blockingRaise () ; return ; }

  if (wMess) /*mhmp 19.11.97 */
    blockingRaiseMess () ;

  graphEventX = XtoUabs(ev->xbutton.x) ;
  graphEventY = YtoUabs(ev->xbutton.y) ; 

  for (menuBox = gActive->nbox ; --menuBox ;)
    { box = gBoxGet (menuBox) ;

      if (box->flag & GRAPH_BOX_NOMENU_FLAG)		    /* Exclude these boxes ala nopick. */
	continue ;

      if (graphEventX >= box->x1 && graphEventX <= box->x2 && 
	  graphEventY >= box->y1 && graphEventY <= box->y2)
	break ;
    }

  if (menuBox == 0)
    menu = gSubDev->menu0 ;	/* background menu */
  else 
    {
      if (box->flag & GRAPH_BOX_MENU_FLAG)
	assFind (gSubDev->boxMenus,
		 assVoid(4*menuBox), &menu) ; /* use box menu */
      else
	menu = gSubDev->menu0 ; /* no menu specified,
				   use the background menu */
    }

  if (!menu)
    return ;

  /* Code to reposition menu.                                                */
  /* Note: we use XTranslateCoordinates because it avoids problems that      */
  /* XtTranslateCoordinates has with the Athena widget set.                  */
  {
    Dimension width, height, border_width ;
    int px, py, npx, npy ;
    Window child = None ;
    /* Find out where we are on the screen.                                  */
    px = ev->xbutton.x;
    py = ev->xbutton.y;
    XTranslateCoordinates(XtDisplay(w),
			  XtWindow(w), 
			  XRootWindowOfScreen(XtScreen(w)),
			  px, py, &npx, &npy,
			  &child) ;
    px = npx; py = npy;

    /* Set the menu position relative to position of pointer, we don't want  */
    /* the pointer sitting within the menu as this results in a menu item    */
    /* being highlighted and can lead to the user accidentally selecting a   */
    /* menu item.                                                            */
    /* Default position of pointer is top left of menu. if too close to      */
    /* bottom/right of screen we position menu further up/to left of pointer.*/
    /* Note the implicit assumption that the menu is NOT bigger than the     */
    /* screen, in 4_8 I will trap this problem when the menu is created.     */
    /*      Also, window positions do not include window border width, which */
    /* we need to take into account when positioning the menu.               */
    XtRealizeWidget (menu) ;	
    XtVaGetValues (menu,
		   XtNwidth, &width,
		   XtNheight, &height,
		   XtNborderWidth, &border_width,
		   NULL) ;

    /* Too far right ?  - put menu to left of pointer.                       */
    if (px + width > WidthOfScreen(XtScreen(menu)))
      {
        px = px - width - (border_width * 2) ;
      }

    /* Too far down ?   - put menu further up.                               */
    if (py + height > HeightOfScreen(XtScreen(menu)))
      {
	py = HeightOfScreen(XtScreen(menu)) - height ;
      }

    /* OK, position the menu.                                                */
    XtVaSetValues (menu, XtNx, px, XtNy, py, NULL) ;
  }


  XtPopupSpringLoaded (menu) ;

  if (doGrab)  /* mieg: i reverse the logic */
    /* Actively grab the pointer */
    XtGrabPointer(menu, True, 
		  (unsigned)
		  (ButtonPressMask | ButtonReleaseMask | PointerMotionMask), 
		  GrabModeAsync, GrabModeAsync,
		  None, RightArrow, ev->xbutton.time);
}


static void keyboardCall (Widget w, XEvent* ev, String *params, Cardinal *num_params)
{
  KeySym keysym ;
  char buffer[2] ;
  int kval ;

  if (!activeSet (w))
    return ;
  /* only check for blocked graphs if they can be blocked */
  if (gActive->block_mode == GRAPH_BLOCKABLE)
    if (blockGraph && !graphContainedIn(blockGraph->id))
      { blockingRaise () ; return ; }

  if (wMess) /*mhmp 19.11.97 */
    blockingRaiseMess () ;

  graphEventX = XtoUabs(ev->xkey.x) ;
  graphEventY = YtoUabs(ev->xkey.y) ;


  /* I do not believe the below is working as intended...we should be looking*/
  /* at the keysym first and then deciding what to do, not looking at the    */
  /* result of XLookupString first...I will have a go at sorting this out    */
  /* at sometime but not now....It is complicated by the fact that we don't  */
  /* do all keyboard stuff via here, some things like Cntl-Y are mapped      */
  /* separately.                                                             */
  if (XLookupString((XKeyEvent*) ev, buffer, 2, &keysym, NULL) )
    kval = buffer[0] ;
  else
    switch (keysym)
      {
      case XK_F1: 
      case XK_F10: 
      case 65386:
	help () ;
	return ;
      case XK_F2:
      case XK_F9:
	/*    case 65478: duplicate */
	graphPrint () ;  
	return ;

/* hard coding keys r1 r2 r3 r5 as arrow keys for the sun under openlook */
      case XK_Delete: kval = DELETE_KEY ; break ;
      case XK_BackSpace: kval = BACKSPACE_KEY ; break ;
      case 65491: 
      case XK_Up: kval = UP_KEY ; break ;
      case 65494: 
      case XK_Down: kval = DOWN_KEY ; break ;
      case 65490:
      case XK_Left: kval = LEFT_KEY ; break ;
      case 65492:
      case XK_Right: kval = RIGHT_KEY ; break ;
      case XK_Home: kval = HOME_KEY ; break ;
      case XK_End: kval = END_KEY ; break ;
      case 65379: kval = INSERT_KEY ; break ;
      case 65365: kval = PAGE_UP_KEY ; break ;
      case 65366: kval = PAGE_DOWN_KEY ; break ;
      case 0:			/* e.g. lower arrows on Sparc */
	switch (((XKeyEvent*)ev)->keycode)
	  {
	  case 27: kval = UP_KEY ; break ;
	  case 31: kval = LEFT_KEY ; break ;
	  case 35: kval = RIGHT_KEY ; break ;
	  case 34: kval = DOWN_KEY ; break ;
	  default: kval = 0 ; break ;
	  }
	break ;
      default: 
/*
 	if (keysym < 65500)
	  printf ("Got unmapped keysym %d (%d)\n",
		  keysym, XtGetActionKeysym (ev, NULL)) ;
*/
	kval = 0 ; 
        break ;
      }
/*
  printf ("buffer, keysym, kval are: %d %d %d\n", 
	  buffer[0],keysym,kval) ;
*/


  /* We don't really want a busy cursor to shown for every key press, e.g.   */
  /* just typing letters, but unfortunately it is not enough to just not     */
  /* do the busycursor call here because each key press causes an expose     */
  /* event which does do the busy cursor !!  To avoid polluting the expose   */
  /* handler, I've added a call to the busy cursor routines that will        */
  /* make the normal busy cursor call do nothing. This should be used with   */
  /* care and in paired calls to start/stop disabling (as below) otherwise   */
  /* the busy cursor stuff will get completely screwed up.                   */
  /* Note also that I have assumed that no one does something like sets up   */
  /* a translation for say "a" that launches some long running bit of code   */
  /* which locks up the application without a busy cursor...                 */
  /*                                                                         */
  if (kval && gActive->func[KEYBOARD])
    {
    if (((keysym >= XK_KP_Space) && (keysym <= XK_KP_9))    /* Keypad keys */
	|| ((keysym >= XK_space) && (keysym <= XK_asciitilde)) /* normal keyboard keys */
	|| keysym == XK_BackSpace || keysym == XK_Delete    /* delete keys */
	|| keysym == XK_Up || keysym == XK_Down	            /* up/down arrow keys */
	|| keysym == XK_Left || keysym == XK_Right)	    /* left/right arrow keys */
      graphSetBusyCursorInactive(TRUE) ;

    graphBusyCursorAll (TRUE);

    (*gActive->func[KEYBOARD])(kval, (GraphKeyModifier) ev->xkey.state) ;

    graphBusyCursorAll (FALSE) ;

    graphSetBusyCursorInactive(FALSE) ;
    }

  return ;
}


static void propertyCall (Widget w, XEvent* ev, String *params, Cardinal *num_params)
  {
  
  /* ACEDB-GRAPH INTERFACE: acedb needs to communicate with other programs   */
  /* using X specific methods.                                               */
  if (getGraphAcedbWprop() != NULL) (getGraphAcedbWprop())(ev) ;

  return ;
  }


static void enterCall (Widget w, XEvent* ev, String *params,
		       Cardinal *num_params)
{
  int i;
  Widget viewport;
  Widget parent;
  Widget last = 0, root, child;

  /* Only dispatch this call if:
     (1) We are actually the widget entered, not our child
     -- or --
     (2) We are the parent, but the child is a scrollbar
  */
  if (((ev->xcrossing.detail == NotifyVirtual) || 
       (ev->xcrossing.detail == NotifyNonlinearVirtual)) &&
      (child = XtWindowToWidget(ev->xcrossing.display, 
				ev->xcrossing.subwindow)) &&
      (XtClass(child)
	== scrollbarWidgetClass))
    return;

  /* unhighlight all pending */
  while (assNext(lastLeft, &last, &root))
    {
      if (root != root_widget)
	messcrash("Can't find the widget you came from.");

      if (last->core.being_destroyed )
	continue ;
      parent = XtParent(last);
      viewport = 
	(XtClass(last) == viewportWidgetClass) ? last :
	  (parent && parent != (void*)(-1) && XtClass(parent) == viewportWidgetClass) ? parent : NULL;
      
      if (viewport != NULL && !assFind(destroyed, viewport, 0))
	{
	  i = 0;
	  /* unhighlight lastLeft */
	  XtSetArg (args[i], XtNborderPixmap, GreyPixmap); i++;
	  XtSetValues (viewport, args, i) ;
	}
    }
  assClear(lastLeft);

  parent = XtParent(w);
  viewport = 
    (XtClass(w) == viewportWidgetClass) ? w :
      (parent && XtClass(parent) == viewportWidgetClass) ? parent : NULL;
  
  if (viewport != NULL && !assFind(destroyed, viewport, 0))
    {
      i = 0;
      /* highlight newly enetered graph */
      XtSetArg (args[i], XtNborderPixmap, XtUnspecifiedPixmap); i++;
      XtSetValues (viewport, args, i) ;
    }
}

static void leaveCall (Widget w, XEvent* ev, String *params,
			  Cardinal *num_params)
{ assInsert(lastLeft, w, root_widget); }

/**********************/

/* Look ahead for function key F4 */
/* Neil also tried #define KEYSYM_CTRL_C 65507 but
   1/ CTRL-C seems to cause a redraw and isn't as tidy as F4 (neil, 8 Jun 92)
   2/ the UNDO key should not be the same as CTRL-C anyway, because
   someone panicking to undo an operation may press the key repeatedly without
   wishing to exit the program
*/

BOOL graphInterruptCalled(void)
{ XEvent xev;
  KeySym keysym;
  char c[2];

  return 
    XCheckTypedEvent(XtDisplay(root_widget), KeyPress, &xev) &&
      !XLookupString((XKeyEvent *) &xev, c, 2, &keysym, NULL) &&
    (keysym == XK_F4 || keysym == XK_Escape)   ;
}




/*****************************************************************************/
/*                  Cut & Paste support.                                     */
/*                                                                           */
/* Now completely rewritten to use the Xt selection mechanism which provides */
/* better support for large data transfers and much else besides. This also  */
/* clears up a number of bugs with the original code which used CUT_BUFFER0  */
/* (e.g. memory leaks).                                                      */
/*                                                                           */
/* Data is now 'cut to/pasted from' the XA_PRIMARY selection buffer, this is */
/* the buffer used by most X clients these days.                             */
/*                                                                           */
/* Note that a better name for graphPostBuffer would have been graphCutBuffer*/
/* in line with the X usage of cut and paste.                                */
/*                                                                           */


/* Returns a pointer to the text pasted from the selection buffer into this  */
/* application.                                                              */
/*                                                                           */
/* Note that the basic method is that this routine is called when the user   */
/* presses keys that result in a paste into this application. This routine   */
/* then registers that it wants the text from the selection buffer and       */
/* registers a routine to receive that text. We then sit in this routine     */
/* dispatching events until our getPasteData routine has been called, in     */
/* other words, this is a blocking routine.                                  */
/*                                                                           */
char *graphPasteBuffer (void)
  {
  char *result = NULL ;
  static graphSelectResult selection = {FALSE, NULL} ;
  XEvent ev ;
  Widget paste_widget = gSubDev->base ;

  XtGetSelectionValue(paste_widget, XA_PRIMARY, XA_STRING, getPasteData,
		      (XtPointer)&selection, XtLastTimestampProcessed(XtDisplay(paste_widget))) ;

  /* Loop until getPasteData has been called and returned.                   */
  while (selection.finished == FALSE)
    {
    XtAppNextEvent (app_con, &ev) ;
    XtDispatchEvent (&ev) ;
    }
  selection.finished = FALSE ;				    /* Reset for next call. */

  if (selection.text != NULL) result = selection.text ;

  return (result) ;
  }

/* This callback routine actually gets the text from the selection buffer,   */
/* it then sets a flag to stop event processing in graphPasteBuffer and      */
/* returns the pasted text.                                                  */
/*                                                                           */
/* NOTE: I have turned off the error messages because they get in the way,   */
/*       it seems that the usual action by applications is to just ring the  */
/*       terminal bell or do nothing at all.                                 */
/*                                                                           */
static void getPasteData(Widget w, XtPointer our_ptr,
			 Atom *selection_buffer, Atom *type,
			 XtPointer value, unsigned long *length, int *format)
  {
  graphSelectResult *selection = (graphSelectResult *)our_ptr ;
  int volume = 100 ;

  if (*type == XA_STRING)
    {
    /* Requested data was successfully converted so we can return it.        */
    static char *buffer = NULL ;
    char *result ;

    /* Get rid of any old buffer and make a copy of the new selection.       */
    if (buffer != NULL) messfree(buffer) ;
    buffer = (char*)messalloc(*length + 1) ;
    result = strncpy(buffer, value, *length) ;
    if (result == NULL) messcrash("unexpected error in strncpy function: %s.", messSysErrorText()) ;

    selection->text = buffer ;				    /* return the string. */
    }
  else if (*type == XT_CONVERT_FAIL)
    {
    /* Bad news, timed out, there is something in the buffer but the app.    */
    /* providing the data has not responded.                                 */
#ifdef GRAPH_DEBUG
    messerror("paste of data failed because program that originally cut the data "
	      "has not responded quickly enough, please try the paste again.") ;
#endif

    XBell(XtDisplay(w), volume) ;
    selection->text = NULL ;
    }
  else if ((char *)value == NULL && *length == 0)
    {
    /* Bad news, no owner for selection any more, or owner cannot convert    */
    /* selected data into our requested type (XA_STRING), or perhaps even    */
    /* more likely, nobody has put anything in the selection buffer yet....  */
#ifdef GRAPH_DEBUG
    messerror("paste of data failed because program that originally cut the data "
	      "has terminated or is unable to convert the data into our requested type "
	      "of XA_STRING.") ;
#endif

    XBell(XtDisplay(w), volume) ;
    selection->text = NULL ;
    }
  else
    {
    /* Bad news, conversion failed, so notify user.                          */
#ifdef GRAPH_DEBUG
    messerror("paste of data failed because program has been sent data in a format "
	      "different to that which it requested (XA_STRING), please report this error.") ;
#endif

    XBell(XtDisplay(w), volume) ;
    selection->text = NULL ;
    }
    
  XtFree(value) ;					    /* Get rid of Xt's copy. */

  selection->finished = TRUE ;				    /* stops local event loop in
							       graphPasteBuffer() */

  return ;
  }





/* Sends the supplied text to the selection buffer.                          */
/*                                                                           */
/* Note that this routine is called when the user presses keys that result   */
/* in a request to send text from this application to the selection buffer.  */
/*         This routine then attempts to get ownership of the selection      */
/* buffer and if successful registers routines to handle sending the data to */
/* the buffer, clearing up when the send is finished and clearing up if      */
/* another program takes control of the selection buffer.                    */
/*                                                                           */
/* (I don't know if I have identified the cut widget correctly, but actually */
/* it doesn't matter too much for our purposes as long as the widget has     */
/* been active recently so that we get a reasonable value time as given by   */
/* XtLastTimestampProcessed.)                                                */
/*                                                                           */
void graphPostBuffer (char *text)
  {
  Boolean select_owner = FALSE ;
  Widget cut_widget = gSubDev->base ;			    /* The cut widget ?? */

  /* Try to get control of the selection buffer.                             */
  select_owner = XtOwnSelection(cut_widget, XA_PRIMARY,
				XtLastTimestampProcessed(XtDisplay(cut_widget)),
				sendCutData, lostCutSelection, finishedCutSend) ;
  
  /* It's unlikely but we can fail.                                          */
  if (select_owner != TRUE)
    {
    messerror("cut of data to selection clip board failed, another application "
	      "has taken control of the clip board, please try again.") ;
    graphSelectCutString_G = NULL ;
    }
  else
    {
    /* OK, we have control of the selection buffer, so make a 'global' copy  */
    /* of the string to be put into the selection buffer.                    */
    if (graphSelectCutString_G != NULL)
      {
      messfree(graphSelectCutString_G) ;
      graphSelectCutString_G = NULL ;
      }

    graphSelectCutString_G = strnew(text, 0) ;
    }

  return ;
  }


/* This routine is the one actually called to send the text to the selection */
/* buffer, it is only called when a request is made by an application to     */
/* get hold of the text from the buffer.                                     */
/*                                                                           */
static Boolean sendCutData(Widget w,  Atom *selection, Atom *target,
			   Atom *type, XtPointer *value, unsigned long *length, int *format)
  {
  Boolean result = FALSE ;
  Boolean debug = FALSE ;
  Display* d = XtDisplay(w);
  XSelectionRequestEvent* req = XtGetSelectionRequest(w, *selection, (XtRequestId)NULL);
							    /* Need req for event stuff... */

  if (*target == XA_TARGETS(d))
    {
    /* This is a standard selection request in which we return what data     */
    /* types we can deal with, we just say 'string'.                         */
    Atom* targetP;
    Atom* std_targets;
    unsigned long std_length;

    XmuConvertStandardSelection(w, req->time, selection, target, type,
				(XPointer*)&std_targets, &std_length, format);
    *value = XtMalloc(sizeof(Atom)*(std_length + 5));
    targetP = *(Atom**)value;
    *targetP++ = XA_STRING;
    *length = std_length + (targetP - (*(Atom **) value));
    memmove( (char*)targetP, (char*)std_targets, sizeof(Atom)*std_length);
    XtFree((char*)std_targets);
    *type = XA_ATOM;
    *format = 32;

    result = TRUE ;
    }
  else if (*target == XA_STRING || *target == XA_TEXT(d) || *target == XA_COMPOUND_TEXT(d))
    {
    /* We support the usual C style strings so deal with them.               */
    if (*target == XA_COMPOUND_TEXT(d)) *type = *target;
    else *type = XA_STRING;
    *length = strlen(graphSelectCutString_G) ;
    *value = graphSelectCutString_G ;
    *format = 8;

    result = TRUE ;
    }
  else
    {
    /* Attempt to deal with any selection types we don't understand but are  */
    /* standard for Xt and will be handled by Xmu/Xt.                        */
    if (XmuConvertStandardSelection(w, req->time, selection, target, type,
				    (XPointer*)value, length, format)) result = TRUE ;
    else if (debug)
      messerror("cannot convert cut data to the standard Xt XA_STRING data format, "
		"please make a note of the data you attempted to cut and report this "
		"failure to the application developers.") ;
    }

  return (result) ;
  }


/* If user cuts text from this application but then changes their mind and   */
/* cuts text from another application, then this routine is called to allow  */
/* us to clear up as our data is no longer wanted.                           */
/*                                                                           */
static void lostCutSelection(Widget w, Atom *selection)
  {

  /* Clear up our copy of the cut text.                                      */
  messfree(graphSelectCutString_G) ;
  graphSelectCutString_G = NULL ;

  return ;
  }


/* This routine is called after sendCutData and after the requesting program */
/* has received the data. In the end we do not need to do anything here, if  */
/* we only wanted to send the data once to the selection buffer we could     */
/* uncomment the XtDisownSelection() and messfree() calls. Currently we will */
/* continue to send the data as long as we own the selection buffer.         */
/*                                                                           */
static void finishedCutSend(Widget w, Atom *selection, Atom *target)
  {

#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
  /* Uncomment this and the transfer will only be done once.                 */
  XtDisownSelection(w, XA_PRIMARY, XtLastTimestampProcessed(XtDisplay(w))) ;
  messfree(graphSelectCutString_G) ;
  graphSelectCutString_G = NULL ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */

  return ;
  }

/********** end of file **********/

