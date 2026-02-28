/*  File: graphcon.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: device independent control routines for graph package
 * Exported functions: graphInit, Finish, Create, Destroy, Activate...
 * HISTORY:
 * Last edited: May 12 13:50 1999 (edgrif)
 * * May 12 13:49 1999 (edgrif): Added new graphSetBlockMode to help
 *              busy cursor setting, added new implementation graphSetBusyCursorAll.
 * * Jan  8 11:44 1999 (edgrif): Correct stupid error in busyCursorInactive_G
 * * Dec 15 23:26 1998 (edgrif): Added busy cursor to graphCreate.
 * * Dec  3 14:46 1998 (edgrif): Put in library version macros.
 * * Nov 19 15:03 1998 (edgrif): Fixed callback func. dec. for
 *              graphTextScrollEditor.
 * * Nov 18 13:51 1998 (edgrif): Added graph package version number.
 * * Oct 21 16:59 1998 (edgrif): Added a layer of message calls which
 *              then call device dependent code in xtsubs etc.. This
 *              allowed the addition of meaningful labels for buttons
 *              e.g. Continue, End etc. It has also reintroduced a
 *              layering that used to be in graph but had disappeared.
 * * Oct 21 11:27 1998 (edgrif): Removed outlandish isGraphics flag
 *              but also the internal isInitialised flag, documented this
 *              in wdoc/graph.html.
 * * Sep 11 09:30 1998 (edgrif): Add messExit function register.
 * * Aug 24 10:19 1998 (fw)
 *	-	graphDeleteContents(dying) code extracted from graphClear(), since the full
 *		graphClear() attempts (WIN32) illegal operations on the destroyed window device
 *	-	Changed "gActive" variables to "dying" midway through graphDestroy():
 *              seems more meaningful?
 * Jan 98  mieg: added in the call backs the value of e->box
       because of the garanteed ANSI-C way parameters are handled in C
       this is correct, the callback deposits on the stack box, val
       the old call-back func only dig for val, which is garanteed ok
 * * Jun 10 15:17 1996 (rd)
 * Created: Thu Jan  2 02:34:40 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: graphcon.c,v 1.7 2015/08/11 22:14:46 mieg Exp $ */
#include "regular.h"

#include "graph_.h"					    /* defines Graph, gXxxx externals, 
							       graphDevXxxx() */
#include "freeout.h"					    /* declares freeOutInit() */
#include "help.h"					    /* For helpOnRegister() */
#include "version.h"					    /* Library version utils */


typedef    BOOL (*EDIT_FUNC)() ;
typedef    BOOL (*FLOAT_EDIT_FUNC)(float,int) ;


/******* externals visible from elsewhere *******/


/* graph package version and copyright string.                     */
/*                                                                */
#define GRAPH_TITLE   "Graph library"
#define GRAPH_DESC    "Sanger Centre Informatics graph library for window control"

#define GRAPH_VERSION 1
#define GRAPH_RELEASE 1
#define GRAPH_UPDATE  1
#define GRAPH_VERSION_NUMBER  UT_MAKE_VERSION_NUMBER(GRAPH_VERSION, GRAPH_RELEASE, GRAPH_UPDATE)

/*
UT_COPYRIGHT_STRING(GRAPH_TITLE, GRAPH_VERSION, GRAPH_RELEASE, GRAPH_UPDATE, GRAPH_DESC)
*/



/* for user selection of legal foreground colours */
FREEOPT  graphColors[] =
{ 
  { 32,   "Colours" },
  { WHITE,        "White" },
  { BLACK,        "Black" },
  { GRAY,         "Gray" },
  { PALEGRAY,     "Pale Gray" },
  { LIGHTGRAY,    "Light Gray" },
  { DARKGRAY,     "Dark Gray" },
  { RED,          "Red" },
  { GREEN,        "Green" },
  { BLUE,         "Blue" },
  { YELLOW,       "Yellow" },
  { CYAN,         "Cyan" },
  { MAGENTA,      "Magenta" },
  { LIGHTRED,     "Light Red" },
  { LIGHTGREEN,   "Light Green" },
  { LIGHTBLUE,    "Light Blue" },
  { DARKRED,      "Dark Red" },
  { DARKGREEN,    "Dark Green" },
  { DARKBLUE,     "Dark Blue" },
  { PALERED,      "Pale Red" },
  { PALEGREEN,    "Pale Green" },
  { PALEBLUE,     "Pale Blue" }, 
  { PALEYELLOW,   "Pale Yellow" },
  { PALECYAN,     "Pale Cyan" },
  { PALEMAGENTA,  "Pale Magenta" },
  { BROWN,        "Brown" },
  { ORANGE,       "Orange" },
  { PALEORANGE,   "Pale Orange" },
  { PURPLE,       "Purple" },
  { VIOLET,       "Violet" },
  { PALEVIOLET,   "Pale Violet" },
  { CERISE,       "Cerise" },
  { MIDBLUE,      "Mid Blue" }
} ;

int graphContrastLookup[] = { 
  BLACK, /* WHITE */
  WHITE, /* BLACK */
  BLACK, /* LIGHTGRAY */
  WHITE, /* DARKGRAY */
  WHITE, /* RED */
  BLACK, /* GREEN */
  WHITE, /* BLUE */
  BLACK, /* YELLOW */
  BLACK, /* CYAN */
  BLACK, /* MAGENTA */
  BLACK, /* LIGHTRED */
  BLACK, /* LIGHTGREEN */
  BLACK, /* LIGHTBLUE */
  WHITE, /* DARKRED */
  WHITE, /* DARKGREEN */
  WHITE, /* DARKBLUE */
  BLACK, /* PALERED */
  BLACK, /* PALEGREEN */
  BLACK, /* PALEBLUE */
  BLACK, /* PALEYELLOW */
  BLACK, /* PALECYAN */
  BLACK, /* PALEMAGENTA */
  WHITE, /* BROWN */
  BLACK, /* ORANGE */
  BLACK, /* PALEORANGE */
  WHITE, /* PURPLE */
  BLACK, /* VIOLET */
  BLACK, /* PALEVIOLET */
  BLACK, /* GRAY */
  BLACK, /* PALEGRAY */
  BLACK, /* CERISE */
  WHITE, /* MIDBLUE */
  0 , 0  /* NUM_TRUECOLOURS, TRANSPARENT, for safety */
  };

/****** externals only for use within graph package *****/
/*                ----                                                       */
/*                 |                                                         */
/*          yeh...sure thing boss..                                          */
Graph_  gActive = 0 ;
Dev     gDev = 0 ;
SubDev  gSubDev = 0 ;
int     gIsCreating = FALSE ;   /* needed for careful activation switch
                                   during creation */
Box     gBox = 0 ;
Stack   gStk = 0 ;
char    *acefont = 0 ;

/****** statics for this file ******/

static char     *progname = NULL ;
static int      ngraph = 0 ;    /* number of graphs made by this process */
static Graph_   graphList = 0 ;
static Associator  graphIdAss = 0 ;

static BOOL busyCursorInactive_G = FALSE ;		    /* Can be used to disable the busy cursor */
							    /* call..a hack, see graphBusyCursor.... */


/***********************************************************/
/***************** initialisation and finish ****************/

/* Must be called before any other graph functions, if this routine is       */
/* called more than once the results are undefined.                          */
/*                                                                           */
void graphInit (int *argcptr, char **argv)
{
  int i ;

  /* If this is not true then it is a gross coding error by a developer.     */
  /* this sort of thing should be via assert calls.                          */
  if (NEVENTS != DESTROY+1)
    messcrash ("Code error: number of graph register events mismatch") ;

  progname = argv[0] ;

  graphIdAss = assCreate() ;

  /* ACEDB-GRAPH INTERFACE: initialise graph package overloadable            */
  /* functions.                                                              */
  initOverloadFuncs() ;

  graphDevInit (argcptr, argv) ;

  /* Register graph versions of message and file calls, note that graphics   */
  /* message routines overload the freeOut message routines, futher calls    */
  /* to freeOutInit() do nothing and therefore can't divert freeOutXxxx      */
  /* away from the graphics routines registered below                        */
  freeOutInit() ;
  messOutRegister (graphOut) ;
  messErrorRegister (graphError) ;
  messExitRegister (graphEnd) ;
  messCrashRegister (graphEnd) ;
  messQueryRegister (graphQuery) ;
  messPromptRegister (graphPrompt) ;

  /* long loops that check messIsInterruptCalled() 
     can now be terminated by pressing F4 */
  messIsInterruptRegister (graphInterruptCalled);

  /* register graphical filechooser */
  filQueryOpenRegister (graphQueryOpen) ;

  /* internal slim-line browser as default                                   */
  helpOnRegister(graphHtmlViewer);
			

  /* Custom font */
  for(i=1;i<*argcptr;i++)
    {
      if (!strcmp(argv[i],"-acefont"))
	{
	  if ((*argcptr - i) < 2)
	    {
	      /* no argument to -acefont option */
	      messExit ("No argument specified for -acefont option");
	    }
	  acefont = strnew(argv[i+1], 0);

	  /* clear both arguments from the list */
	  for(i+=2;i<*argcptr;i++)
	    argv[i-2] = argv[i];
	  argv[*argcptr-2] = 0;
	  (*argcptr) -= 2;
	  break;
        }
    }

  return ;
}

/*******/



/* Routines registered with the message package by graphIinit.               */
GRAPH_FUNC_DCL void graphOut(const char *text)
  {
  char *label = GRAPH_OUT_MSG ;

  graphWinOut(text, label) ;

  return ;
  }

GRAPH_FUNC_DCL void graphError(const char *text)
  {
  char *label = GRAPH_ERROR_MSG ;

  graphWinOut(text, label) ;

  return ;
  }

GRAPH_FUNC_DCL void graphEnd(const char *text)
  {
  char *label = GRAPH_END_MSG ;

  graphWinOut(text, label) ;

  return ;
  }


BOOL graphQuery(const char *text)
  {
  BOOL result ;
  
  result = graphWinQuery(text) ;

  return result ;
  }



BOOL graphPrompt(const char *prompt, char *dfault, char *fmt)
  {
  BOOL result ;
  
  result = graphWinPrompt (prompt, dfault, fmt) ;

  return result ;
  } 



/***********************/

void graphFinish ()
  {

  while (gActive)
    graphDestroy () ;

  graphDevFinish () ;

  assDestroy(graphIdAss) ;
  graphIdAss = 0 ;

  messOutRegister (0) ;
  messErrorRegister (0) ;
  messExitRegister (0) ;
  messCrashRegister (0) ;
  messQueryRegister (0) ;
  messPromptRegister (0) ;
  filQueryOpenRegister (0) ;
  helpOnRegister(0);

  return ;
  }

/***********************/

BOOL graphActivateChild (void)
{
  if (gActive->children)
    return  graphActivate (gActive->children->id) ;
  return FALSE ;
}

/***********************/

BOOL graphActivate (Graph gId)
{ Graph_ g ;
  char *cp0 = 0 ;

  if (gActive && gId == gActive->id)
    return TRUE ;
 
  if (!gId || !assFind(graphIdAss, cp0 + gId, &g))
    return FALSE ;

  if (!gIsCreating && gDev && gSubDev)
    graphDevActivate (FALSE) ;
  gActive = g ;
  gDev = gActive->dev ;
  gSubDev = gActive->subdev ;
  gStk = gActive->stack ;
  gBox = gActive->nbox ? gBoxGet (gActive->currbox) : 0 ;
  if (gDev && gSubDev)
    graphDevActivate (TRUE) ;

  return TRUE ;
}

Graph graphActive (void) { return gActive ? gActive->id : 0 ; }



/********************* create and destroy ************************/

static MENUOPT defaultMenu [] = { 
  { graphDestroy,	"Quit" },
  { help,		"Help" },
  { graphPrint,		"Print"},
  { 0, 0 }
} ;

void setTextAspect (Graph_ g)	/* used in WIN32 code */
{
  int dx,dy ;

  if (!gFontInfo (0, &dx, &dy))
    messcrash ("Can't get font info for default font") ;
			/* default font dimensions are device dependent */
  g->aspect = dy / (float) dx ;
  g->xFac = dx + 0.00001 ; /* bit to ensure correct rounding */
  g->yFac = dy + 0.00001 ;
}

static void typeInitialise (Graph_ g, int type)
	/* This is called before the window is generated, and
	   so can not refer to w,h.
	   After the window is generated the resize procedure
	   will be called once - so any initialisation for that
	   must be placed here.
	   After that xFac,yFac,xWin,yWin must be set, so that
	   drawing can take place.
	*/
{
  g->type = type ;

  switch (type)
    {
    case PLAIN: case MAP_SCROLL:
      g->ux = g->uy = 0.0 ;
      g->uw = g->uh = g->aspect = 1.0 ;
      break ;
    case TEXT_FULL_EDIT:/* default 80x25 text window */
    case TEXT_SCROLL:
    case TEXT_FULL_SCROLL:
      g->uw = 80 ;
      g->uh = 25 ;	/* note deliberate fall through here */
    case TEXT_FIT:
      g->ux = g->uy = 0 ;
      setTextAspect (g) ;
      g->xWin = g->yWin = 0 ;
      break ;
    case PIXEL_SCROLL: 
      g->uw = g->uh = 100 ;
    case PIXEL_FIT:
      g->ux = g->uy = 0 ;
      g->aspect = 1.0 ;
      g->xFac = g->yFac = 1.0 ;
      g->xWin = g->yWin = 0 ;
      break ;
    default: 
      messcrash ("Invalid graph type %d requested",type) ;
    }
}

static void positionCheck (float *px, float *py, float *pw, float *ph)
	/* tile if x,y are 0.0 - use old values if w,h are 0.0 */
{
  static float oldx = 0.0 ;
  static float oldw = 0.0 ;
  static float oldy = 0.0 ;
  float maxx, maxy ;

  graphScreenSize (&maxx, &maxy, 0, 0, 0, 0) ;

  if (!*py) *py = oldy ;
  if (!*px) 
    { *px = oldx + oldw ;
      if (*px + *pw/2 > maxx)
	{ *px = 0 ;
	  *py += 0.1 ;
        }
    }

  if (*pw < 0) *pw = 0 ;
  if (*pw > maxx) *pw = maxx ;

  if (*px < 0) *pw = 0 ;
  if (*px + *pw/2 > maxx) *px = maxx - *pw/2 ;

  if (*ph < 0) *ph = 0 ;
  if (*ph > maxy) *ph = maxy ;

  if (*py < 0) *py = 0 ;
  if (*py + *ph/2 > maxy) *py = maxy - *ph/2 ;

  oldx = *px ;
  oldy = *py ;
  oldw = *pw ;
}

Graph_ graphCreateStruct (int type)
{ char *cp0 = 0 ;
  AC_HANDLE handle = handleCreate () ;
  Graph_ graph = (Graph_) handleAlloc (0, handle, GRAPH_SIZE) ;
  int i ;

  graph->magic = GRAPH_MAGIC ;
  graph->id = ++ngraph ;
  graph->window_being_constructed = TRUE ;
  graph->block_mode = GRAPH_BLOCKABLE ;			    /* By default graphs are non-blocking */
                                                            /* but can themselves be blocked. */
  assInsert(graphIdAss, cp0 + graph->id, (void *)graph) ;
  graph->handle = handle ;
  graph->clearHandle = handleHandleCreate (handle) ;

  graph->boxes = arrayHandleCreate (64, BoxStruct, handle) ;
  graph->stack = stackHandleCreate (1024, handle) ;
  graph->boxstack = stackHandleCreate (32, handle) ;
  graph->bubbleInfo = arrayHandleCreate (64, BUBBLEINFO, handle) ;
  graph->bubbleDict = dictHandleCreate (1024, handle) ;
  dictAdd (graph->bubbleDict, "null", 0) ; /* so real entries are non-zero */
  graph->nbox = 0 ;
  graph->assoc = assHandleCreate(handle) ;
  typeInitialise (graph, type) ;

        /* initialise statics for the graph */
  graph->linewidth = 0.002 ;  /* mieg 2002: keep this value in synch with graphgd.c (same comment) */
  graph->pointsize = 0.005 ;
  graph->textheight = 0.0 ;
  graph->color = FORE_COLOR ;
  graph->textFormat = PLAIN_FORMAT ;
  for (i = 0 ; i < NEVENTS ; i++)
    graph->func[i] = 0 ;

  return graph ;
}


Graph graphCreate (int type, const char *nam, float x, float y, float w, float h)
{ 
  Graph_ graph;

#if !defined(DEBUG)
  if (type == TEXT_FULL_EDIT)
    { type = TEXT_FULL_SCROLL ;
      messout ("TEXT_FULL_EDIT window not yet available - using TEXT_FULL_SCROLL") ;
    }
#endif

  graph = graphCreateStruct (type) ;

  if (!nam)
    nam = messprintf ("%s%d", progname, ngraph) ;
  graph->name = (char *) handleAlloc (0, graph->handle, strlen(nam) + 1) ;
  strcpy (graph->name, nam) ;

  graph->nxt = graphList ;
  graphList = graph ;         /* add to list of graphs */
  graph->parent=0;
  graph->children=0;          /* initialize list of subgraphs */

  if (gDev && gSubDev)
    graphDevActivate (FALSE) ;  /* preparation for devCreate */
  gIsCreating = TRUE ;


  /* Actually create the new window                                          */
  positionCheck (&x, &y, &w, &h) ;
  graphDevCreate (graph,x,y,w,h) ;


  /* Set a busy cursor on all windows so that user will be warned not to use */
  /* the application until window is completely drawn, the cursors will be   */
  /* reset by graphxt.c which tests the window_being_constructed flag to see */
  /* if the cursor should be reset (n.b. mapping can occur for reasons other */
  /* than window creation, in which case we don't want to reset the cursors.)*/
  graphBusyCursorAll (TRUE);


  graphActivate (graph->id) ;
  gIsCreating = FALSE ;

  graphClear () ;                       /* initialises box, stack structures*/

#if !defined(MACINTOSH)
  graphMenu (defaultMenu) ;
#endif

  return graph->id ;
}


Graph graphSubCreate (int type, float ux, float uy, float uw, float uh)
{ 
  Graph_ graph;
  int x, y, w, h;


#if !defined(DEBUG)
  if (type == TEXT_FULL_EDIT)
    messcrash ("Graph type TEXT_FULL_EDIT is not yet implemented.\n");
#endif

  graph = graphCreateStruct (type) ;

  graph->name = 0;                         /* subgraphs are not named */
  
  graph->parent = gActive;

  graph->nxt = graph->parent->children;
  graph->parent->children = graph;

  if (gDev && gSubDev)
    graphDevActivate (FALSE) ;  /* preparation for devCreate */
  gIsCreating = TRUE ;
  graph->dev=graph->parent->dev;


  x = uToXabs(ux) ;
  if (x < 0) x = 0 ;
  if (x > graph->parent->w) x = graph->parent->w ;
  
  y = uToYabs(uy) ;
  if (y < 0) y = 0 ;
  if (y > graph->parent->h) y = graph->parent->h ;

  w = uToXrel(uw);
  h = uToYrel(uh);

  graphSubDevCreate (graph, x, y, w, h) ;


  graphActivate (graph->id) ;
  gIsCreating = FALSE ;
  

  /* Set a busy cursor on all windows so that user will be warned not to use */
  /* the application until window is completely drawn, the cursors will be   */
  /* reset by graphCBMap which tests the window_being_constructed flag to see*/
  /* if the cursor should be reset (n.b. mapping can occur for reasons other */
  /* than window creation, in which case we don't want to reset the cursors.)*/
  graphBusyCursorAll (TRUE);


  graphClear () ;                       /* initialises box, stack structures*/
  
#if !defined(MACINTOSH)
  graphMenu (defaultMenu) ;
#endif
  return graph->id ;
}

/*******************************/

void graphDestroy ()
{ char *cp0 = 0 ;
  Graph_ dying = gActive ;
  Graph_ g, *list, next ;
  static Associator ass = 0 ;            /* used to store set of dying graphs */

  if (!gActive)
    return ;
  if (!ass)
    ass = assCreate () ;
  if (!assInsert (ass,dying,&g))        /* dying already there */
    return ;

  g = dying->children;                  /* destroy all children */
  while (g)
    { graphActivate(g->id);
      next = g->nxt ;		/* save next pointer to survive
				   the graphDestroy of the current
				   linked list element */
      graphDestroy();
      g = next ;
    }

  graphActivate(dying->id);

  if (dying->subdev)                    /* device specific functions */
    graphSubDevDestroy();

  if (dying->dev && dying->parent==0)   /* outermost window */
    graphDevDestroy() ;

  devMenuDestroy () ;

  graphDeleteContents(dying) ;			/* kills anything attached to boxes */
						/* extracted from graphClear(); rbrusk 9/96 */
  /* Changed "gActive" variables to "dying" here; makes more sense? rbrusk 9/96 */
  if (dying->func[DESTROY])           	/* user registered functions */
    (*(dying->func[DESTROY]))() ;     	/* can not use graph, as gDev gone */

  if (dying->parent) 
    list = &dying->parent->children ; 
  else 
    list = &graphList ;

  if (dying == *list)               /* find in graph list and unlink */
    *list = dying->nxt ;
  else
    { 
      for (g = *list ; g && g->nxt != dying ; g = g->nxt)
	;
      if (!g)
        messcrash ("Dying graph not in graph list") ;
      g->nxt = dying->nxt ;
    }

  assRemove(graphIdAss, cp0 + dying->id) ;

  if (*list)
    graphActivate ((*list)->id) ;         /* must activate something else */
  else if (dying->parent)
    graphActivate (dying->parent->id);
  else
    { gActive = 0 ; gDev = 0 ; gSubDev = 0 ; gStk = 0 ; gBox = 0 ;
    }

  dying->magic = 0 ;
  assRemove (ass, dying) ;
  if (dying->handle)
    { AC_HANDLE h = dying->handle;
      dying->handle = 0; /* for possible recursion */
      handleDestroy (h) ;
    }
}

/*********************************/

BOOL graphExists (Graph gId)
{ Graph_ g ;
  char *cp0 = 0 ;
  
  return  
    gId && assFind(graphIdAss, cp0 + gId, &g) ;
}

/*********************************/

BOOL graphContainedIn (Graph parent)
/* TRUE if active graph is a subgraph (not necessarily immediate) of parent */
{ char *cp0 = 0 ;
  Graph_ dad;
  Graph_ ptr;

  if (!assFind(graphIdAss, cp0 + parent, &dad))
    return FALSE;

  if (gActive == dad)
    return TRUE;

  ptr = gActive->parent;
  while (ptr)
    {
      if (ptr == dad)
	return TRUE;
      ptr = ptr->parent;
    }
  return FALSE;
}

/*********************************/

Graph graphParent (void)
{ if (!gActive || !gActive->parent) 
    return 0 ;
  else
    return gActive->parent->id ;
}

/*********************************/

void graphCleanUp (void)  /* kill all windows except current active one */
{ 
  Graph_ gkill, g = graphList, gSave = gActive ;

  if (!gSave)
    return ;

  while (g)
    { gkill = g ;
      g = g->nxt ;
      
      if (gkill != gSave && graphActivate (gkill->id))
	graphDestroy() ;
    }
}



/* Graphs can be of three types for purposes of user interactions:           */
/*       GRAPH_BLOCKABLE - normal window (default), input blocked when a     */
/*                         modal dialog is posted.                           */
/*   GRAPH_NON_BLOCKABLE - normal window but input is not blocked when a     */
/*                         modal dialog is posted.                           */
/*        GRAPH_BLOCKING - modal dialog window, cannot itself be blocked,    */
/*                         there should only be one of these at a time.      */
/*                                                                           */
void graphSetBlockMode(GraphWindowBlockMode mode)
{

  if (mode != GRAPH_BLOCKABLE && mode != GRAPH_NON_BLOCKABLE && mode != GRAPH_BLOCKING)
    {
      messerror(GRAPH_INTERNAL_LOGIC_ERROR_MSG "graphSetBlockMode received invalid blocking mode. "
		GRAPH_RESTART_ERROR_MSG) ;
    }

 gActive->block_mode = mode ;
}



/**************** Code to implement busy cursors       ***********************/
/* See graphInternals.html for a discussion of how to use these calls.       */

/* graphBusyCursorAll will turn on/off the busy cursor on all currently      */
/* visible graph windows, if you add new busy cursor calls you MUST add them */
/* in 'on/off' paired calls otherwise the stack logic below will break. We   */
/* keep the stack because often cursor calls are nested and this can lead    */
/* to the busy cursor being turned off at inappropriate moments.             */
/*                                                                           */
void graphBusyCursorAll (BOOL on)
{ 
  static int cursor_stack = 0 ;

#ifdef GRAPH_DEBUG
  if (cursor_stack == 0)
    printf("\n") ;
  printf("Busy cursors, request is: %s,  stack is: %d\n", (on ? "ON" : "OFF"), cursor_stack) ;
#endif

  /* If we have a logic error, just reset and try to carry on.               */
  if (cursor_stack < 0)
    {
    messerror(GRAPH_INTERNAL_LOGIC_ERROR_MSG
	      "busy cursor stack out of sync, "
	      "the 'watch' cursor may be turned on/off at the wrong times. "
	       GRAPH_RESTART_ERROR_MSG) ;
    cursor_stack = 0 ;
    }
  else
    {
    /* Has graph disabled busy cursors (see graphSetBusyCursorInactive) ??   */
    if (busyCursorInactive_G == FALSE)
      {
	/* We always need to turn cursors on (a new window may have been       */
	/* created), we only turn them all off again when we return to the     */
	/* bottom of the stack.                                                */
	Graph_ g ;
	Graph_ oldActive = gActive ;
	void *v = 0 ;

	/* The logic here is that we only turn cursors off when we are at    */
	/* the bottom of the stack or for a blocking window after it has     */
	/* been constructed (the user should respond to this sort of window. */
	/* We always turn cursor on because any application request may      */
	/* create new windows that will need a busy cursor setting while     */
	/* they are created.                                                 */
	if (on == FALSE && cursor_stack > 1)
	  {
	    while (assNext (graphIdAss, &v, &g))
	      {
		gActive = g ;

		if (g->window_being_constructed == FALSE && g->block_mode == GRAPH_BLOCKING)
		  {
		    graphBusyCursor (on) ;
		  }

#ifdef GRAPH_DEBUG
		printf("Busy cursors, actually set: %s\n", (on ? "ON" : "OFF")) ;
#endif
	      }

	  }
	else
	  {
	    while (assNext (graphIdAss, &v, &g))
	      {
		gActive = g ;

		/* We do busy cursors for:                                         */
		/*     - all windows that are being constructed                    */
		/*     - all displayed windows that are blockable                  */
		if ((g->window_being_constructed == TRUE)
		    || (g->window_being_constructed == FALSE && g->block_mode == GRAPH_BLOCKABLE))

		  {
		    graphBusyCursor (on) ;
		  }

#ifdef GRAPH_DEBUG
		printf("Busy cursors, actually set: %s\n", (on ? "ON" : "OFF")) ;
#endif
	      }

	  }

	gActive = oldActive ;
      
      }

    /* Always adjust the stack.                                              */
    if (on)
      cursor_stack++ ;
    else
      cursor_stack-- ;
    }

  return ;
}


/* graphSetBusyCursorInactive controls whether the above routine will        */
/* actually do anything. This is necessary because for some application      */
/* callbacks we don't want to show a busy cursor, BUT the callbacks          */
/* themselves may call other graph routines that set the cursors (e.g.       */
/* redraw) and we need to stop all nested busy cursor calls.                 */
/* A good example would be dragging the cursor which would show a constantly */
/* flickering cursor unless we turned off cursor setting in this way.        */
void graphSetBusyCursorInactive(BOOL inactive)
  {
  busyCursorInactive_G = inactive ;

  return ;
  }




/**** routine to register the Help  ****/

char *graphHelp (char *item)
{ char *old = gActive->help ;
  if (item && *item)
    gActive->help = item ;

  return old ;
}

/**** routine to register functions for events (e.g. mouse) ****/

GraphFunc uGraphRegister (int i, GraphFunc proc) /* i is a GraphEvent */
{
  GraphFunc old = gActive->func[i] ;
  gActive->func[i] = proc ;
  return old ;
}

/**** button package ****/

int graphButton (const char* text, VoidRoutine func, float x, float y)
{ 
  int k = graphBoxStart () ;  

  graphText (text, x + XtoUrel(3), y + YtoUrel(2)) ;
  graphRectangle (x, y, gBox->x2 + XtoUrel(3), gBox->y2) ;
  gBox->flag |= GRAPH_BOX_BUTTON_FLAG ;	/* Added for mac CLH */
  graphBoxEnd () ;
  graphBoxDraw (k, BLACK, WHITE) ;
  if (!gActive->buttonAss)
    gActive->buttonAss = assHandleCreate (gActive->handle) ; 
  assInsert (gActive->buttonAss, assVoid(k*4), (void*) func) ;
  graphBoxInfo (k, 0, 
		strnew (messprintf ("BUTTON:\"%s\"", text), 
			gActive->handle)) ;
  return k ;
}

int graphColouredButton (const char *text,
			 ColouredButtonFunc func,
			 void *arg,
			 int fg,
			 int bg,
			 float x,
			 float y)
{ int k = graphBoxStart();
  
  graphText(text, x + XtoUrel(3), y + YtoUrel(2));
  graphRectangle(x, y, gBox->x2 + XtoUrel(3), gBox->y2);
  gBox->flag |= GRAPH_BOX_BUTTON_FLAG; /* Added for mac CLH */
  graphBoxEnd();
  graphBoxDraw (k, fg, bg);
  if (!gActive->buttonAss)
    gActive->buttonAss = assHandleCreate (gActive->handle);
  assInsert(gActive->buttonAss, assVoid(k*4), (void*) func) ;
  assInsert(gActive->buttonAss, assVoid(k*4 + 1), arg);
  assInsert(gActive->buttonAss, assVoid(k*4 + 2), assVoid(fg));
  assInsert(gActive->buttonAss, assVoid(k*4 + 3), assVoid(bg));
  graphBoxInfo (k, 0, 
		strnew (messprintf ("BUTTON:\"%s\"", text), 
			gActive->handle)) ;
  return k;
}

int graphButtons (MENUOPT *buttons, float x0, float y0, float xmax)
{
  float x = x0 ;
  float y = y0 ;
  int n=0, box = 0, len, ix, iy ;
  float ux, uy ;
  
  if (!buttons) return 0 ;

  graphTextInfo (&ix, &iy, &ux, &uy) ;
  uy += YtoUrel(8) ;

  while (buttons->text)
    { 
      len = strlen (buttons->text) ;
      
      if (!(len == 0 && buttons->f == menuSpacer))
	{ /* menu option is not a separator */
	  if (x + ux*len + XtoUrel(15) > xmax && x != x0)
	    { x = x0 ; y += uy ; }
	  box = graphButton (buttons->text, buttons->f, x, y) ;
	  x += ux*len + XtoUrel(15) ; 
	  ++n ;
	}
       ++buttons ;
    }

  return box + 1 - n ;
}

int graphColouredButtons (COLOUROPT *buttons, float x0, float *y0p, float xmax)
{
  float x = x0 ;
  float y = *y0p ;
  int n=0, box = 0, len, ix, iy ;
  float ux, uy ;
  COLOUROPT *old= buttons;

  if (!buttons) return 0 ;

  graphTextInfo (&ix, &iy, &ux, &uy) ;
  uy += YtoUrel(8) ;

  while (buttons->text) /* terminate on NULL text */
    { 
      len = strlen (buttons->text) ;
      
      if (x + ux*len + XtoUrel(15) > xmax && x != x0)
	{ x = x0 ; y += uy ; *y0p = y ; }
      box = graphColouredButton (buttons->text, 
				 buttons->f,
				 buttons->arg,
				 buttons->fg,
				 buttons->bg,
				 x, y) ;
      x += ux*len + XtoUrel(15) ; 
      ++n ;
      if (buttons->next == buttons || buttons->next == old) break;
      if (buttons->next) 
	buttons = buttons->next;
      else 
	buttons++;
    }  

  return box + 1 - n ;
}

/************** graph associator stuff ***************/

BOOL uGraphAssociate (void* in, void* out)
{
  if (gActive)
    return assInsert (gActive->assoc,in,out) ;
  else
    return FALSE ;
}

BOOL uGraphAssFind (void *in, void* *out)
{
  if (gActive)
    return assFind (gActive->assoc,in,out) ;
  else
    return FALSE ;
}

BOOL graphAssRemove (void* in)
{
  if (gActive)
    return assRemove (gActive->assoc,in) ;
  else
    return FALSE ;
}

/********** utility functions to let users get at handles **********/

AC_HANDLE graphHandle (void)  /* gets freed when graph is freed */
{ if (gActive)
    return gActive->handle ;
  else
    return 0 ;
}

AC_HANDLE graphClearHandle (void) /* freed when graph is cleared */
{ if (gActive)
    return gActive->clearHandle ;
  else
    return 0 ;
}

/*******************************************************************/
/****** graph editors - single item labelled textEntry boxes ******/

typedef struct {
  char *text ;
  char *format ;
  void *p ;
  int len ;
  int box,index; /* for radio buttons */
  union { BOOL (*i)(int,int) ; BOOL (*s)(char*,int) ;BOOL (*f)(float,int) ;} func ;
} EDITOR ;

static BOOL callNoError ;

static void callEditor (char *text, Graph_ g)
{ 
  EDITOR *e ;
  union {int i ; float f ; char *s ;} x ;
  int i ;

  if (!g->editors)
    messcrash ("callEditor() called without editors") ;

  for (i = 0 ; i < arrayMax(g->editors) ; ++i)
    if ((e = arrp(g->editors, i, EDITOR))->text == text)
      { if(e->format[0] != 'b') freeforcecard (text) ;
	if (!freecheck (e->format) && e->format[0] != 'b')
	  { messout ("Format does not check for entry \"%s\"", text) ;
	    callNoError = FALSE ;
	  }
	else
	  switch (*e->format)
	    { 
	    case 'b':
	      break;
	    case 'i':
	      freeint (&x.i) ;
	      if (e->func.i && !(*e->func.i)(x.i,e->box))
		{ messout ("Entry \"%s\" does not check", text) ;
		  callNoError = FALSE ;
		}
	      else
		*(int*)e->p = x.i ;
	      break ;
	    case 'f':
	      freefloat (&x.f) ;
	      if (e->func.f && !(*e->func.f)(x.f,e->box))
		{ messout ("Entry \"%s\" does not check", text) ;
		  callNoError = FALSE ;
		}
	      else
		*(float*)e->p = x.f ;
	      break ;
	    case 'w':
	    case 't':
	      x.s = freeword() ;
	      if (strlen (x.s) > e->len)
		{ messout ("Entry \"%s\" is too long (max %d)", x.s, e->len) ;
		  callNoError = FALSE ;
		}
	      else if (e->func.s && !(*e->func.s)(x.s,e->box))
		{ messout ("Entry \"%s\" does not check", x.s) ;
		  callNoError = FALSE ;
		}
	      else
		strncpy (e->text, x.s, e->len) ;
	      break ;
	    }
	return ;
      }

  messcrash ("callEditor called on bad text") ;
}

void callBackEditor (char *text) { callEditor (text, gActive) ; }

static int drawEditor (char *label, char *text, int len, int wlen,
		       float x, float y)
{
  graphText (label, x, y) ;
  x += strlen (label) + 0.5 ;
  return graphTextScrollEntry (text, len, wlen, x, y, callBackEditor) ;
}

static EDITOR *addEditor (char *text, void *p, char *format, int len)
{ 
  EDITOR *e ;

  if (!gActive->editors)
    gActive->editors = arrayHandleCreate (8, EDITOR, gActive->clearHandle) ;

  e = arrayp(gActive->editors, arrayMax(gActive->editors), EDITOR) ;
  e->text = text ;
  e->format = format ;
  e->p = p ;
  e->len = len ;

  return e ;
}

static void drawAltEditor(char *label, float x, float y, EDITOR *e)
{
  int box1,box2;

  graphText(label, x+1.5, y);
/*  x += strlen (label) + 0.5 ;*/
  
  box2 = graphBoxStart();
  graphFillArc(x+0.5, y+0.5, 0.4, 0, 360);
  graphBoxEnd();
 
  box1 = graphBoxStart();
  graphArc(x+0.5, y+0.5, 0.7, 0, 360);
  gBox->flag |= GRAPH_BOX_TOGGLE_FLAG;
  graphBoxEnd();
  
  graphAssociate(assVoid(box1),assVoid(arrayMax(gActive->editors)-1));

  graphBoxDraw(box1, BLACK, TRANSPARENT);
  if (*(BOOL*)e->p)
    graphBoxDraw(box2, BLACK, TRANSPARENT);
  else
    graphBoxDraw(box2, WHITE, TRANSPARENT);

  graphBoxSetPick(box2, FALSE);
}


/******* public routines ********/

void editorToggleChange(int box)
{
  EDITOR *e,*e1;
  int i;

  graphAssFind(assVoid(box), &i);
  e = arrp(gActive->editors,i,EDITOR);
  if(e->len < 1000){  /* if it's a toggle button */
    *(BOOL*)e->p = !*(BOOL*)e->p;
    
    graphBoxDraw(box, BLACK, TRANSPARENT);
    if (*(BOOL*)e->p)
      graphBoxDraw(box-1, BLACK, TRANSPARENT);
    else
      graphBoxDraw(box-1, WHITE, TRANSPARENT);
  }
  else{ /* if it's a radio button */
    for(i=0;i<arrayMax(gActive->editors);i++){
      e1 = arrp(gActive->editors,i,EDITOR);
      if(e1->len == e->len) /* i.e. same tag */
	graphBoxDraw(e1->box, WHITE, TRANSPARENT);
    }
    graphBoxDraw(e->box,BLACK,TRANSPARENT);
    graphAssFind(assVoid(e->len), &i);
    e1 = arrayp(gActive->editors,i,EDITOR);
    
    *(int*)e1->p = e->index; 
  }
}

void graphToggleEditor (char *label, BOOL *p, float x, float y)
{
  char *text = (char*) halloc (16, gActive->clearHandle) ;
  EDITOR *e = addEditor(text,(void*)p,"b", 0);

  strncpy (text,label, 15) ;

  drawAltEditor(label,x,y,e);
}

static void colourChange(KEY key,int box)
{
  int i;
  EDITOR *e;

  if(graphAssFind(assVoid(box), &i)){
    e = arrp(gActive->editors,i,EDITOR);
    *(int*)e->p = key;
    if(e->len != -2)
      graphBoxDraw(box,graphContrastLookup[*(int*)e->p], *(int*)e->p);
    else
      graphBoxDraw(e->box,BLACK,*(int*)e->p);
  }
}
void graphRedoColourBoxes()
{
  int i;
  EDITOR *e;

  for(i=0;i<arrayMax(gActive->editors);i++){
    e = arrp(gActive->editors,i,EDITOR);
    if(e->format[0] == 'b' && e->len == -1)
      graphBoxDraw(e->box,graphContrastLookup[*(int*)e->p], *(int*)e->p);
    else if(e->format[0] == 'b' && e->len == -2)
      graphBoxDraw(e->box,BLACK,*(int*)e->p);
  }
}

void graphColourEditor(char *label, char *text, int *p,float x, float y)
{
  EDITOR *e;
  char *text2 = (char*) halloc (strlen(label)+1, gActive->clearHandle) ;

  strcpy(text2,label);
  e = addEditor(text2,(void*)p,"b", 0);

  e->len = -1;
  e->box = graphBoxStart();
  graphRectangle(x,y,x+3+strlen(text),y+1);
  graphText(text,x+1.0,y);
  graphBoxEnd();
  graphBoxFreeMenu(e->box,colourChange,graphColors);
  if(text[0] != ' '){
    e->len = -1;
    graphBoxDraw(e->box,graphContrastLookup[*(int*)e->p],*(int*)e->p);
  }
  else{
    graphBoxDraw(e->box,BLACK,*(int*)e->p);    
    e->len = -2;
  }
  x += strlen(text) + 4;
  graphText (label,x,y);
  graphAssociate(assVoid(e->box),assVoid(arrayMax(gActive->editors)-1));

}			

/***********************************************************************************************************/
/*
*p will contain the index for the chosen radio button
 tag is the reference number for a radio button.
*/
int graphRadioCreate(char *text,int *p,float x, float y)
{
  int tag = 1000,i;
  EDITOR *e = addEditor(text,0,"b", 0);

  graphText(text,x,y);
  for(i=0;i<arrayMax(gActive->editors);i++){
    e = arrp(gActive->editors,i,EDITOR);
    if(e->format[0] == 'b' && e->len > 1000){
      tag = e->len +1;
    }
  }
  e = addEditor (text, (void*)p, "b", 0) ;
  e->len = tag;
  graphAssociate(assVoid(tag),assVoid(arrayMax(gActive->editors)-1));
  return tag;
}
/***********************************************************************************************************/
/*
tag ->   reference number for the radio button.
index -> this will be returned to e->p in the create routine if this radio is chosen.
*/
 
void graphAddRadioEditor(char *text, int tag, int index, float x, float y)
{
  EDITOR *e = addEditor(text,0,"b", 0);
  int box1;
  
  e->len = tag;

  graphText(text, x, y);
  x += strlen (text) + 0.5 ;
  
  e->text = 0;
  e->box  = graphBoxStart();
  graphFillArc(x+0.5, y+0.5, 0.4, 0, 360);
  graphBoxEnd();
 
  box1 = graphBoxStart();
  graphArc(x+0.5, y+0.5, 0.7, 0, 360);
  gBox->flag |= GRAPH_BOX_TOGGLE_FLAG;
  graphBoxEnd();
  
  graphAssociate(assVoid(box1),assVoid(arrayMax(gActive->editors)-1));

  e->index = index;
  graphBoxDraw(box1, BLACK, TRANSPARENT);
  if (e->index == 0)
    graphBoxDraw(e->box, BLACK, TRANSPARENT);
  else
    graphBoxDraw(e->box, WHITE, TRANSPARENT);

  graphBoxSetPick(e->box, FALSE);
}

void graphSetRadioEditor(int tag,int index)
{
  EDITOR *e;
  int i;

  for(i=0;i<arrayMax(gActive->editors);i++){
    e = arrp(gActive->editors,i,EDITOR);
    if(e->format[0] == 'b' && e->len == tag && e->index == index){
      editorToggleChange(e->box+1);
    }
  }
}

int graphIntEditor (char *label, int *p, float x, float y, 
		    BOOL (*checkFunc)(int))
{
  char *text = (char*) halloc (16, gActive->clearHandle) ;
  EDITOR *e = addEditor (text, (void*)p, "iz", 0) ;
  e->func.i = (EDIT_FUNC) checkFunc ;
  sprintf (text, "%d", *p) ;
  return  e->box = drawEditor (label, text, 15, 8, x, y) ;
}

int graphFloatEditor (char *label, float *p, float x, float y,
		    BOOL (*checkFunc)(float))
{
  char *text = (char*) halloc (16, gActive->clearHandle) ;
  EDITOR *e = addEditor (text, (void*)p, "fz", 0) ;
  e->func.f = (FLOAT_EDIT_FUNC) checkFunc ;
  sprintf (text, "%.4g", *p) ;
  return e->box = drawEditor (label, text, 15, 8, x, y) ;
}

int graphWordEditor (char *label, char *text, int len, float x, float y,
		    BOOL (*checkFunc)(char*))
{
  EDITOR *e = addEditor (text, 0, "wz", len) ;
  e->func.s = (EDIT_FUNC) checkFunc ;
  return  e->box = drawEditor (label, text, len, len < 16 ? len : 16, x, y) ;
}

int graphTextEditor (char *label, char *text, int len, float x, float y,
		    BOOL (*checkFunc)(char*))
{
  EDITOR *e = addEditor (text, 0, "tz", len) ;
  e->func.s = (EDIT_FUNC) checkFunc ;
  return  e->box = drawEditor (label, text, len, len < 16 ? len : 16, x, y) ;
}

int graphTextScrollEditor (char *label, char *text, int len, int wlen, float x, float y,
		    BOOL (*checkFunc)(char *text, int box))
{
  EDITOR *e = addEditor (text, 0, "tz", len) ;
  e->func.s = (EDIT_FUNC) checkFunc ;
  return  e->box = drawEditor (label, text, len, wlen, x, y) ;
}

#ifdef JUNK
int graphMultiLineEditor (char *label, char *text, int len, int width, int height, 
					 float x, float y,
					 BOOL (*checkFunc)(char* text, int box)) 
{
  GMLEDITOR *e ;
  
  if (!gActive->gmlEditors)
    gActive->gmlEditors = arrayHandleCreate (8, GMLEDITOR, gActive->clearHandle) ;
  
  e = arrayp(gActive->mlEditors, arrayMax(gActive->gmlEditors), GMLEDITOR) ;
  e->text = label ;
  e->text = text ;
  e->len = len ;
  e->len = width ;
  e->len = height ;
  e->func.s = (EDIT_FUNC) checkFunc ;
  return  e->box = gmleDraw (e) ;
}
#endif
BOOL graphCheckEditors (Graph graph, BOOL callCheckFunctions)
{ char *cp0 = 0 ;
  int i ;
  Graph old = 0 ;
  Graph_ g ;

  if (!graph || !assFind(graphIdAss, cp0 + graph, &g)) 
    return FALSE ;
  if (!g->editors) 
    return TRUE ;
  if (g != gActive) 
    old = gActive->id ;

  callNoError = TRUE ;
  for (i = 0 ; i < arrayMax(g->editors) ; ++i)
    callEditor (arrp(g->editors, i, EDITOR)->text, g) ;
  
  if (old) graphActivate (old) ;

  return callNoError ;
}

BOOL graphUpdateEditors (Graph graph)
{ char *cp0 = 0 ;
  int i ;
  Graph old = 0 ;
  Graph_ g ;
  EDITOR *e ;

  if (!graph || !assFind(graphIdAss, cp0 + graph, &g)) 
    return FALSE ;
  if (!g->editors) 
    return TRUE ;
  if (g != gActive) 
    { old = gActive->id ;
      graphActivate (g->id) ;
    }

  for (i = 0 ; i < arrayMax(g->editors) ; ++i)
    { e = arrp(g->editors, i, EDITOR) ;
      switch (*e->format)
	{ 
	case 'i': sprintf (e->text, "%d", *(int*)e->p) ; break ;
	case 'f': sprintf (e->text, "%.4g", *(float*)e->p) ; break ;
	}
      graphTextScrollEntry (e->text, 0, 0, 0, 0, 0) ;
    }

  if (old) graphActivate (old) ;

  return callNoError ;
}

/********* utility to draw a triangle for menu selection *******/

int graphMenuTriangle (BOOL filled, float x, float y)
{
  int i = 3, box;

  box = graphBoxStart();

  graphLine(x, y+0.25, x+1, y+0.25);
  graphLine(x+1, y+0.25, x+0.5, y+0.75);
  graphLine(x+0.5, y+0.75, x, y+0.25);

  if (filled)
    {
      Array temp = arrayCreate(2*i, float) ;

      array(temp, 0, float) = x;
      array(temp, 1, float) = y+0.25;
      array(temp, 2, float) = x+1;
      array(temp, 3, float) = y+0.25;
      array(temp, 4, float) = x+0.5;
      array(temp, 5, float) = y+0.75;
      graphPolygon(temp);
      arrayDestroy(temp);
    }

  graphBoxEnd();

  return box;
}

/************ end of file *****************/
