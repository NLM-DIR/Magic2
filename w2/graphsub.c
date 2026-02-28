/*  File: graphsub.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: graph package - device independent routines acting on
 *	the current active graph.
 * Exported functions: (many)
 * HISTORY:
 * Last edited: Mar 24 14:33 1999 (edgrif)
 * * Mar 24 14:31 1999 (edgrif): Fix SANgc03641, added graphBoxSetMenu()
 * * Mar  5 13:05 1999 (edgrif): Fix SANgc03582, missing call to busy cursor.
 * * Dec 15 23:30 1998 (edgrif): Added busy cursor to gLeftDown, gMiddleDown.
 * * Nov 12 13:36 1998 (edgrif): Added code to trap NULL return from
 *              graphPasteBuffer(), otherwise we bomb out.
 * * Sep 30 09:59 1998 (edgrif): Replace #ifdef ACEDB with a call to
 *              the acedb/graph interface. Correct invalid use of
 *              '? :' construct for output to freeOutf or fprintf.
 * * Apr 30 16:21 1997 (rd)
 *		-   In graphBoxStart(), removed newGraphBox()
 *		-	graphDeleteContents() extracted from graphClear() and invoked by
 *			both graphClear(), and in graphDestroy() in lieu of graphClear()
 *
 * * Jun 6 16:34 1996 (rbrusk):
 *	- graphTextEntry2() etc. now all "void (*fn)(char*)" callbacks
 * * Jun  5 01:37 1996 (rd)
 * * Jan  5 13:50 1994 (mieg): restore gActive after button action
 * * Aug 19 12:02 1992 (rd): polygon squares.
 * * Aug 19 12:01 1992 (mieg): textPtrPtr
 * * Feb 11 00:11 1992 (mieg): graphClassEntry -> graphCompletionEntry
 * * Jan 28 1992 (neil): save existing drag routines over a drag operation, 
                         restore them on button up;
                         (changed graphBoxDrag, upBox)
 * * Dec  1 19:21 1991 (mieg): textEntry recognition of arrow keys
 * * Dec  1 19:20 1991 (mieg): graphClassEntry: autocompletes on tab
 * * Oct 10 14:51 1991 (rd): cleaned graphBoxClear - now a toggle
 * * Oct  8 23:00 1991 (mieg): added graphBoxClear
 * * Oct  7 23:00 1991 (rd): added graphBoxShift
 *-------------------------------------------------------------------
 */

/* $Id: graphsub.c,v 1.11 2019/03/04 22:44:01 mieg Exp $ */

#include "regular.h"
#include "freeout.h"					    /* for freeOutF */
#include "graph_.h"     /* defines externals within graph package */
#include "key.h"

#if defined(MACINTOSH)
  #define isascii(x) ((x) < 128)
#endif
#include <ctype.h>

/*******************************/

void graphFacMake (void)	/* for PLAIN and MAP graphs */
{
  float xfac,yfac ;
  float asp ;

  xfac = gActive->w / gActive->uw ;
  yfac = gActive->h / gActive->uh ;

  if (gActive->aspect)
    { asp = gPixAspect * gActive->aspect ;
      if (xfac*asp > yfac)
	xfac = yfac/asp ;
      yfac = xfac * asp ;
    }
  else
    gActive->aspect = yfac/xfac ;

  gActive->xFac = xfac ;
  gActive->yFac = yfac ;
  gActive->xWin = uToXrel (gActive->ux) ;
  gActive->yWin = uToYrel (gActive->uy) ;
}

void gUpdateBox0 (void)
{
  Box box ;
  
  if (!gActive->nbox)
    return ;

  box = gBoxGet(0) ;
  box->x1 = gActive->ux ; box->x2 = gActive->ux + gActive->uw ;
  box->y1 = gActive->uy ; box->y2 = gActive->uy + gActive->uh ;
}

void graphPlainBounds (float ux, float uy, float uw, float uh, float aspect)
{                                 /* makes ->fac, ->xWin, ->yWin */
  gActive->ux = ux ;
  gActive->uy = uy ;
  gActive->uw = uw ;
  gActive->uh = uh ;
  gActive->aspect = aspect ;
  graphFacMake() ;
  gUpdateBox0 () ;
}

void graphTextInfo (int *dx, int *dy, float *cw, float *ch)
{
  if (!gFontInfo (uToYrel(gActive->textheight), dx, dy)) /* number of pixels of a character */
    messcrash ("Can't get font info for current font") ;
     /* default font dimensions are device dependent */
  if (cw && dx)
    *cw = XtoUrel(*dx) ; /* char width in user's coord */
  if (ch && dy)
    *ch = YtoUrel(*dy) ; /* char height in user's coord */
}

void graphFitBounds (int *nx, int *ny)
{
/* commented out fw  15.03 1994*/
/*  if (gActive->type != TEXT_FIT && gActive->type != TEXT_SCROLL)
    messcrash 
      ("FitBounds called on invalid graph type %d",gActive->type) ;
*/
  if (nx)
    *nx = gActive->w / gActive->xFac ;
  if (ny)
    *ny = gActive->h / gActive->yFac ;
}

/* other graphXxxBounds() are in graphdev.c since they actively
   change the size of the scrolled window, which is a device
   dependent operation.
*/

/********************* box control ********************************/

Box gBoxGet (int k)
{ 
  if (k < 0 || ! gActive->boxes || k >= arrayMax (gActive->boxes)) /* mieg 2105_04_05 was k >= gActive->nbox but it crashed occasionanly on cursor boxes */
    messcrash ("Cannot get box %d - index out of range",k) ;
  return arrp (gActive->boxes, k, BoxStruct) ;
}

int graphLastBox(void)
{
  return gActive->nbox - 1;
}

int graphBoxStart (void)
{
  Box box ;
  int k = gActive->nbox++ ;
  int col = k ? gBox->bcol : BACK_COLOR ;
#ifdef BETTER_BOX_SHIFT
  int parentid = k ? gBox->id : 0 ;
#endif

  box = arrayp (gActive->boxes, k, BoxStruct) ;	/* invalidates gBox */

  box->linewidth = gActive->linewidth ;
  box->pointsize = gActive->pointsize ;
  box->textheight = gActive->textheight ;
  box->fcol = gActive->color ;
  box->format = gActive->textFormat ;
  box->bcol = col;
  box->x1 = 1000000.0 ; box->x2 = -1000000.0 ;
  box->y1 = 1000000.0 ; box->y2 = -1000000.0 ;
  box->flag = 0 ;
#ifdef BETTER_BOX_SHIFT
  box->id = k ;
  box->parentid = parentid ;
#endif

  push (gStk, BOX_START, int) ; push (gStk, k, int) ;
  box->mark = stackMark (gStk) ;

  gBox = box ;
  if (k)
    push (gActive->boxstack, gActive->currbox, int) ;
  gActive->currbox = k ;

  return k ;
}

void graphBoxEnd (void)
{ Box box = gBox ;

  if (gActive->currbox == 0)
    messcrash ("Tried to end box 0") ;

  push (gStk, BOX_END, int) ;
  gActive->currbox = pop (gActive->boxstack, int) ;
  gBox = gBoxGet (gActive->currbox) ;

  if (box->x1 < gBox->x1) gBox->x1 = box->x1 ;
  if (box->x2 > gBox->x2) gBox->x2 = box->x2 ;
  if (box->y1 < gBox->y1) gBox->y1 = box->y1 ;
  if (box->y2 > gBox->y2) gBox->y2 = box->y2 ;
}

BOOL graphBoxClear (int k)
/* toggles appearance status of box and 
   redraws/clears box area as necessary */
{
  Box box = gBoxGet(k) ;

  if (!k ||                /* do not clear background box */
      k >= gActive->nbox)  /* avois stupid crash */
    return TRUE ;  

  box->x2 = -box->x2 ;		/* will not redraw if x2 < x1 */
  box->x1 = -box->x1 ;

  if (box->x1 > box->x2)	/* Clear the screen */
   {
     graphClipDraw (uToXabs(-box->x1), uToYabs(box->y1),
		    uToXabs(-box->x2), uToYabs(box->y2)) ;
     return FALSE;
   }
  else
    {
      graphBoxDraw (k, -1, -2) ;	/* expose box */
      return TRUE;
    }
}

BOOL graphBoxMarkAsClear (int k)
{
  Box box = gBoxGet(k) ;

  if (!k)  /* do not clear background box */
    return TRUE ;  

  box->x2 = -box->x2 ;		/* will not redraw if x2 < x1 */
  box->x1 = -box->x1 ;

  return (box->x1 > box->x2) ? FALSE : TRUE ;
}

void graphBoxDim (int k, float *x1, float *y1, float *x2, float *y2)
{ Box box = gBoxGet (k) ;
  if (x1) *x1 = box->x1 ; 
  if (x2) *x2 = box->x2 ;
  if (y1) *y1 = box->y1 ; 
  if (y2) *y2 = box->y2 ;
}

int graphBoxAt(float x, float y, float *rx, float *ry)
{ int i;
  Box box = 0;

  for(i = gActive->nbox; i--; )
    { box = gBoxGet(i);
      if (x >= box->x1 && x <= box->x2 &&
	  y >= box->y1 && y <= box->y2)
	break;
    }

  if (i != -1)
    {
      if (rx) *rx = x - box->x1;
      if (ry) *ry = y - box->y1;
    }

  return i;
}

void graphBoxSetPick (int k, BOOL pick)
{
  Box box = gBoxGet (k) ;

  if (pick)
    box->flag &= ~GRAPH_BOX_NOPICK_FLAG ;
  else
    box->flag |= GRAPH_BOX_NOPICK_FLAG ;
}


/* A new and maybe temp. function to set a flag to exclude a box from menu   */
/* click logic...if this flag is on, then the box is excluded altogether     */
/* from being selected by a right (menu) click.                              */
void graphBoxSetMenu (int k, BOOL menu)
{
  Box box = gBoxGet (k) ;

  if (menu)
    box->flag &= ~GRAPH_BOX_NOMENU_FLAG ;
  else
    box->flag |= GRAPH_BOX_NOMENU_FLAG ;
}

/******************************************/
/*********** bubble info stuff ************/

void graphBubbleDisplay (int box)
{
  BUBBLEINFO *inf ;
  int i ;

  if (!gActive || box <= 0 )
    return ;
  if (!gActive->bubbleDict || !gActive->bubbleInfo)
    return ;

  graphUnMessage()  ;
  i = arrayMax(gActive->bubbleInfo) ;
  inf = arrayp(gActive->bubbleInfo, 0, BUBBLEINFO) - 1 ;
  while (inf++, i--)
    {
      if (inf->box == box)
	{
	  graphMessage (messprintf("%s\n%s",
				   (inf->bName ? dictName (gActive->bubbleDict, inf->bName) : ""),
				   (inf->key ? name(inf->key) : ""))) ;
	  break ;
	}
    } 
}

/******************************************/

void graphBubbleInfo (int box, KEY key, char *ficheName, char *bubbleName) 
{ 
  BUBBLEINFO *bubble ;

  if (!gActive || box <= 0)
    return ;
  if (!gActive->bubbleDict || !gActive->bubbleInfo)
    return ;
  if (0 && strstr (bubbleName, "primed in"))
    invokeDebugger() ;
  bubble = arrayp(gActive->bubbleInfo, box, BUBBLEINFO) ;
  bubble->box = box ;
  bubble->key = key ;
  if (ficheName && *ficheName)
    dictAdd (gActive->bubbleDict, ficheName, &(bubble->fName)) ;
  if (bubbleName && *bubbleName)
    {
      char *cp = strnew (bubbleName, 0) ;
      char *cq, *cr ;

      for (cq = cp ; cq && *cq ; cq = cr)
	{
	  cr = cq ;
	  while ((cr = strchr (cr + 1, ',')))
	    if (cr > cq + 20 && cr[1] != '\n')
	      {
		*cr++ = '\n' ;
		break ;
	      }
	}
      dictAdd (gActive->bubbleDict,cp, &(bubble->bName)) ;
      messfree (cp)  ;
    }
}

/******************************************/

static void graphDoDumpBubble (void) /* dumps on active out */
{ 
  int i ;
  BUBBLEINFO *bubble ;
  Box box ;
  DICT *dict ;
  
  if (!gActive || !gActive->bubbleInfo || !(dict = gActive->bubbleDict) )
    return ;
  
  for (i = 0, bubble = arrayp(gActive->bubbleInfo, 0, BUBBLEINFO) ;
       i < arrayMax(gActive->bubbleInfo) ; bubble++, i++)
    {      
      if (!bubble->box)
	continue ;      
      
      box = gBoxGet(bubble->box) ;
      if (!box || box->x2 < box->x1 || box->y2 < box->y1) /* ignore empty or hidden boxes */
	continue ;
      
      freeOutf ("%d\t%d\t%d\t%d\t",               /* box coordinates in pixcell units */
		uToXabs(box->x1), uToYabs(box->y1), 
		uToXabs(box->x2), uToYabs(box->y2)) ;
      freeOutf ("%s\t", bubble->key ? name(bubble->key) : "0") ;
      freeOutf ("%s\t%", bubble->fName ? dictName (dict, bubble->fName) : "0") ;
      freeOutf ("%s\t%", bubble->bName ? dictName (dict, bubble->bName) : "0") ;
      freeOut ("\n");
    }
    
  return ;
}

/******************************************/

BOOL graphDumpBubble (char *fileName)
{
  FILE *fil = 0 ;
  BOOL level = 0 ;
  if (!fileName || !*fileName)
    return FALSE ;

  if (strcmp (fileName, "-"))
    {
      fil = filopen (fileName,"bubble", "w") ;
      if (fil)
	level = freeOutSetFile (fil) ;
      else
	return FALSE ;
    }
  graphDoDumpBubble () ;
  
  if (level)
    freeOutClose(level) ;
  if (fil)
    filclose (fil) ;

  return TRUE ;
}


/******************************************/
/************** box info stuff ************/

void graphBoxInfo (int box, KEY key, char *text)
{ 
  BOXINFO *inf ;

  if (!gActive)
    return ;
  if (!gActive->boxInfo)
    gActive->boxInfo = arrayHandleCreate (32, BOXINFO, gActive->clearHandle) ;

  inf = arrayp(gActive->boxInfo, box, BOXINFO) ;
  inf->key = key ;
  inf->text = text ;
}


void graphBoxInfoFile (FILE *fil)
  { 
  int i ;
  BOXINFO *inf ;
  Box box ;


  if (!gActive || !gActive->boxInfo) return ;

  for (i = 0 ; i < arrayMax(gActive->boxes) ; ++i)
    {
    inf = arrayp(gActive->boxInfo, i, BOXINFO) ;

    if (!(inf->key || inf->text)) continue ;


    box = gBoxGet(i) ;
    if (box->x2 < box->x1 || box->y2 < box->y1) /* ignore empty or hidden boxes */
	continue ;


    /* 2 Mar 1998 LS, so that boxes can go to gifaceserver */
    if (fil != NULL) fprintf (fil, "%d   %d %d %d %d",
			      i, uToXabs(box->x1), uToYabs(box->y1), 
			      uToXabs(box->x2), uToYabs(box->y2)) ;
    else freeOutf ("%d   %d %d %d %d",
		  i, uToXabs(box->x1), uToYabs(box->y1), 
		  uToXabs(box->x2), uToYabs(box->y2)) ;



    if (inf->key)
      { 

      /* ACEDB-GRAPH INTERFACE: if acedb has registered a function to print  */
      /* the acedb class then call it.                                       */
      if (getGraphAcedbClassPrint() != NULL) (getGraphAcedbClassPrint())(fil, inf->key) ;
      else
	{
	/* I have no idea why this is done again, it has already been done     */
	/* above and box is not used below...remove ??                         */
	box = gBoxGet(i) ;

	if (fil != NULL) fprintf (fil, "  %d", inf->key) ;
	else  freeOutf ("  %d", inf->key); 
	}

      }

    if (inf->text)
      {
      if (fil != NULL) fprintf (fil, "  %s", inf->text) ;
      else freeOutf ("  %s", inf->text);
      }

    if (fil != NULL) fprintf (fil, "\n") ;
    else freeOut ("\n");
    }


  return ;
  }




void graphBoxInfoWrite (void)
{ 
  FILE *fil ;

  if (!gActive || !gActive->boxInfo ||
      !(fil = filqueryopen (0, 0, "box", "w", "File for box information")))
    return ;

  graphBoxInfoFile (fil) ;
  filclose (fil) ;
}

/********************* text entry box *******************************/
/** structure definition must be here because used in gLeftDown() ***/

typedef struct EntryBoxStruct
  { int		magic ;
    char*	text ;
    char*	cp ;
    int		len ;
    int		box ;
    int		cursorBox ;
    int         insertMode ;
    void	(*action)(char*) ;
    GraphCompletionFunc  completionFunction ;
    struct EntryBoxStruct *nxt ;
    char*	winText ;
    int		winPos ;
    int		winLen ;
  } *EntryBox ;

static int EBOX_MAGIC = 187648 ;
static int TBOX_HANDLE ;

static int graphTextEntry2 (GraphCompletionFunc f, char* text, 
			    int len, int winlen,
			    float x, float y, void (*fn)(char*)) ;
static void graphPositionCursorFromEvent (EntryBox ebox) ;
static void entryBoxEntry (int key) ;

/********* box picking and dragging **********/

static void (*dragRoutine)(float*,float*,BOOL) ;
static GraphFunc leftDragR, leftUpR, middleDragR, middleUpR;
  /*saved over dragging operations, by graphBoxDrag, restored by upBox*/
static float oldx,oldy,xbase,ybase,xuser,yuser ;
static int draggedBox ;

void gifLeftDown (int x, int y) {
  gLeftDown(x/gActive->xFac, y/gActive->yFac);
}

void gLeftDown (float x, float y)
{
  int i = -1 ;
  Box box = 0 ;
  VoidRoutine buttonFunc ;
  typedef void (*MouseFunc)(double, double) ;
  typedef void (*PickFunc)(int, double, double) ;
  MouseFunc mfn ;
  PickFunc pfn ;
  Graph old ;
  void *arg;
  int bg, fg;

  for (i = gActive->nbox ; i-- ;)
    { box = gBoxGet (i) ;
      if (box->flag & GRAPH_BOX_NOPICK_FLAG)
	continue ;
      if (x >= box->x1 && x <= box->x2 && 
	  y >= box->y1 && y <= box->y2)
	break ;
    }

  if (i == -1) 	/* outside even the whole drawing area! */
    return ;

  if (box->flag & GRAPH_BOX_ENTRY_FLAG)	/* make this the active entry */
    { EntryBox ebox ;
      if (graphAssFind (&TBOX_HANDLE, &ebox))
	for ( ; ebox ; ebox = ebox->nxt)
	  if (i == ebox->box || i == ebox->cursorBox)
	    { graphTextEntry2 (0, ebox->text, 0, 0, 0, 0, 0) ;
	      graphPositionCursorFromEvent (ebox) ;
	    }
		/* NB - continue because user may need to do something */
    }
  else if (box->flag & GRAPH_BOX_TOGGLE_FLAG)
    editorToggleChange(i);
  else if (i && gActive->buttonAss && 
	   assFind (gActive->buttonAss, assVoid(i*4), &buttonFunc))
    { graphBoxDraw (i, WHITE, BLACK) ;
      old = graphActive() ;
				/* do the action */
      if (assFind(gActive->buttonAss, assVoid(i*4 + 1), &arg))
	{ ColouredButtonFunc cbf = (ColouredButtonFunc) buttonFunc ;
	graphBusyCursorAll (TRUE);
	cbf (arg) ; /* mieg, to please the IBM compiler */
	graphBusyCursorAll (FALSE) ;
	}
      else
	{
	graphBusyCursorAll (TRUE) ;
	(*buttonFunc)() ;
	graphBusyCursorAll (FALSE) ;
	}

	/* redraw button */
      if (graphActivate(old) &&
	  i < gActive->nbox && (box = gBoxGet(i)) &&
	  box->fcol == WHITE && box->bcol == BLACK &&
	  gActive->buttonAss &&
	  assFind (gActive->buttonAss, assVoid(i*4), &buttonFunc))
	{ fg = (assFind(gActive->buttonAss, assVoid(i*4 + 2), &arg)) 
	    ? assInt(arg) : BLACK ;
	  bg = (assFind(gActive->buttonAss, assVoid(i*4 + 3), &arg))
	    ? assInt(arg) : WHITE ;
	  graphBoxDraw (i, fg, bg) ;
	}
      return ;
    }

  /* If the user clicked the left mouse button on something, this is where   */
  /* the applications callback function gets called, we use busy cursors to  */
  /* warn the user we are busy.                                              */
  /* This applies to PICK & LEFT_DOWN.                                        */
  if (gActive->func[PICK])
    {
    pfn = (PickFunc)(gActive->func[PICK]) ;

    graphBusyCursorAll (TRUE);

    (*pfn)(i,(double)(x - box->x1),(double)(y - box->y1)) ;

    graphBusyCursorAll (FALSE) ;

    return ;
    }

  if (gActive->func[LEFT_DOWN])
    {
    mfn = (MouseFunc)(gActive->func[LEFT_DOWN]);

    graphBusyCursorAll (TRUE);
    
    (*mfn)((double)x,(double)y) ;

    graphBusyCursorAll (FALSE) ;

    }

  return ;
}


void gMiddleDown (float x, float y)
{
  typedef void (*MouseFunc)(double, double) ;
  MouseFunc mfn ;
  EntryBox ebox ;
  int i ;
  Box box = 0 ;

  for (i = gActive->nbox ; i-- ;)
    { box = gBoxGet (i) ;
      if (box->flag & GRAPH_BOX_NOPICK_FLAG)
	continue ;
      if (x >= box->x1 && x <= box->x2 && 
	  y >= box->y1 && y <= box->y2)
	break ;
    }
  if (box->flag & GRAPH_BOX_ENTRY_FLAG &&
      graphAssFind (&TBOX_HANDLE, &ebox))	/* make active and paste */
    { for ( ; ebox ; ebox = ebox->nxt)
	if (i == ebox->box || i == ebox->cursorBox)
	  { graphTextEntry2 (0, ebox->text, 0, 0, 0, 0, 0) ;
	    graphPositionCursorFromEvent (ebox) ;
	    entryBoxEntry (25) ; /* CTRL-Y is paste action */
	  }
      /* NB - continue because user may need to do something */
    }

  /* If the user clicked the middle mouse button on something, this is where */
  /* the applications callback function gets called, we use busy cursors to  */
  /* warn the user we are busy.                                              */
  if (gActive->func[MIDDLE_DOWN])
    {
    mfn = (MouseFunc)(gActive->func[MIDDLE_DOWN]);

    graphBusyCursorAll (TRUE);
    (*mfn)((double)x,(double)y) ;
    graphBusyCursorAll (FALSE) ;
    }


  return ;
}


BOOL gIdle (float x, float y)
{
  IdleFunc   ifn ;
  int i ;
  Box box = 0 ;
  BOOL done = TRUE ;

  for (i = gActive->nbox ; i-- ;)
    { box = gBoxGet (i) ;
      if (box->flag & GRAPH_BOX_NOPICK_FLAG)
	continue ;
      if (x >= box->x1 && x <= box->x2 && 
	  y >= box->y1 && y <= box->y2)
	{ x -= box->x1 ; y -= box->y1 ; }   /* go to box coords */
	break ;
    }

  /* If the user clicked the middle mouse button on something, this is where */
  /* the applications callback function gets called */
  if (gActive->func[IDLE])
    {
    ifn = (IdleFunc)(gActive->func[IDLE]);

    done = (*ifn)(i, (double)x,(double)y) ;
    }

  return done ;
}


static void dragBox (double x, double y)
{
  float oldxuser = xuser, oldyuser = yuser ;
 
  xbase += x-oldx ;
  ybase += y-oldy ;
  xuser = xbase;
  yuser = ybase;
  (*dragRoutine)(&xuser,&yuser,FALSE) ;
  if (xuser != oldxuser || yuser != oldyuser)
    { graphXorBox (draggedBox,oldxuser,oldyuser) ;
      graphXorBox (draggedBox,xuser,yuser) ;
    }
  oldx = x ; oldy = y ;
}

static void upBox (double x, double y)
{ 
  graphXorBox (draggedBox,xuser,yuser) ;
  xbase += x-oldx ;
  ybase += y-oldy ;
  (*dragRoutine)(&xbase,&ybase, TRUE) ;

  /*.....restore the saved methods*/
  graphRegister(LEFT_DRAG, leftDragR);
  graphRegister(LEFT_UP, leftUpR);
  graphRegister(MIDDLE_DRAG, middleDragR);
  graphRegister(MIDDLE_UP, middleUpR);

  return;
}

void graphBoxDrag (int k, void (*clientRoutine)(float*,float*,BOOL))
{
  Box box = gBoxGet(k) ;

  dragRoutine = clientRoutine ;
  draggedBox = k ;

  leftDragR=graphRegister(LEFT_DRAG, dragBox);
  leftUpR=graphRegister(LEFT_UP, upBox);
  middleDragR=graphRegister(MIDDLE_DRAG, dragBox);
  middleUpR=graphRegister(MIDDLE_UP, upBox);
    /*saving the existing methods for restoration by upBox*/

  xbase = xuser = box->x1 ; ybase = yuser = box->y1 ;
  oldx = graphEventX ; oldy = graphEventY ; 
  (*dragRoutine) (&xuser,&yuser,FALSE) ;

  graphXorBox (draggedBox,xuser,yuser) ;

  return;
}

/*****************************************************************/
/******* text entry box code - structure before gLeftDown() ******/

static void graphUnSetCursor(EntryBox ebox)
{ 
  graphBoxDraw (ebox->cursorBox, PALEGREEN, TRANSPARENT) ;
}

/******************************/

static void graphSetCursor(EntryBox ebox)
{ 
  float x1,x2,y1,y2 ;
  int offset = ebox->cp - ebox->text ;

  if (ebox->winText != ebox->text)
    { if (offset < 0)
	offset = 0 ;
      if (offset < ebox->winPos)
	ebox->winPos = offset ;
      else if (offset >= ebox->winPos + ebox->winLen)
	ebox->winPos = offset - ebox->winLen + 1 ;
      strncpy (ebox->winText, ebox->text + ebox->winPos, ebox->winLen) ;
    }
      
  graphBoxDim (ebox->box, &x1, &y1, &x2, &y2) ;
  graphBoxShift (ebox->cursorBox, 
		 x1 + UtextX*(ebox->cp - ebox->text - ebox->winPos),
		 y1) ;
  graphBoxDraw (ebox->cursorBox, RED, TRANSPARENT) ;
}

static void graphPositionCursorFromEvent (EntryBox ebox)
{ 
  float x1,x2,y1,y2 ;

  graphBoxDim (ebox->box, &x1, &y1, &x2, &y2) ;
  if (graphEventY >= y1 && graphEventY <= y2)
    { int n = (graphEventX - x1) + ebox->winPos ;
      if (n <= 0 || n > strlen(ebox->text))
	n = strlen(ebox->text) ;
      ebox->cp = ebox->text + n ;
    }

  graphSetCursor (ebox) ;
}

/******************************/

static void entryBoxEntry (int key)
{
  EntryBox ebox ;
  int n, oldbox ;
  char *cp, *cq, *cr ;
  Graph currGraph = graphActive () ;

  if (!graphAssFind (&TBOX_HANDLE, &ebox) ||
        ebox->magic != EBOX_MAGIC )
    messcrash ("entryBoxEntry() can't find entryBox") ;

  cp = ebox->cp ;
  if (strlen(ebox->text) > ebox->len)
    messcrash ("Over flow in entryBoxEntry") ;
  if (cp < ebox->text || cp >= ebox->text + ebox->len)
    cp = ebox->text + strlen(ebox->text) ;

  switch (key)
    {
    case '\t': /* auto complete */
      if (ebox->completionFunction)
	n = (*ebox->completionFunction)(ebox->text, ebox->len) ;
      else
	n = 0 ;
      cp = ebox->text + strlen(ebox->text) ;
      if (n != 1)
	break ;			/* note fall-through if one match only */
    case RETURN_KEY:
      ebox->cp = cp ;
      graphSetCursor(ebox) ;
      graphBoxDraw (oldbox = ebox->box, BLACK, RED) ;
      if (ebox->action)
	(*(ebox->action))(ebox->text) ;
      graphActivate (currGraph) ;  /* May work or not */
      if (graphAssFind (&TBOX_HANDLE, &ebox) &&  /* reget after action */
	  ebox->magic == EBOX_MAGIC &&
	  ebox->box == oldbox &&
	  ebox->box < gActive->nbox)
	{ Box testBox = gBoxGet (ebox->box) ;
	  /* test if ebox->action hasn't changed the entry boxes
	     colour to make it look inactive */
	  if (testBox->fcol == BLACK && testBox->bcol == RED)
	    /* action hasn't change it, so redraw it YELLOW */
	    graphBoxDraw (ebox->box, BLACK, YELLOW) ;
	}
      return ;   /* because action may have destroyed or graphRedrawn */
    case 12:			/* C-l center box on cursor */  
      if (cp > ebox->text + ebox->winLen/2)
	ebox->winPos = cp - ebox->text - ebox->winLen/2 ;
      else
	ebox->winPos = 0 ;
      break;
    case INSERT_KEY:
      ebox->insertMode = 1 - ebox->insertMode ;
      break ;
    case DELETE_KEY:
    case BACKSPACE_KEY:		/* delete one char to my left */
      if (cp == ebox->text)
	break ;
      cq = --cp ;
      while (*cq)
	{ *cq = *(cq + 1) ;
	  cq++ ;
	}
      break ;
    case 4:			/* C-d delete at cp */
      cq = cp ;
      while(*cq)
	{ *cq = *(cq + 1) ;
	  cq ++ ;
	}
      break ;
    case 23:			/* C-w delete a word leftwards at cp */
      if (cp == ebox->text)
	break ;
      cr = cq = cp ;
      while (*cq == ' ' && cq > ebox->text) cq-- ;
      while (*cq != ' ' && cq > ebox->text) cq-- ;
      cp = cq ;
      while ((*cq++ = *cr++)) ;	/* copy */
      while (*cq) *cq++ = 0 ;	/* clean up */
      break ;
    case LEFT_KEY:
    case 2:			/* C-b */
      if (cp > ebox->text)
	cp-- ;
      break ;
    case RIGHT_KEY:
    case 6:			/* C-f */
      if (*cp)
	cp++ ;
      if (cp > ebox->text + ebox->len && ebox->len > 0 )
	cp-- ;
      break ;
    case HOME_KEY:
    case 1:			/* C-a */
      cp = ebox->text ; 
      break ;
    case 11:			/* C-k kill from cp to the right */
      cq = ebox->text + ebox->len ;
      graphPostBuffer (cp) ;
      while (cq >= cp) 
	*cq-- = 0 ;
      break ;
    case 21:			/* C-u kill whole entry */
      cq = ebox->text + ebox->len ;
      cp = ebox->text ;
      graphPostBuffer (cp) ;
      while (cq >= cp) 
	*cq-- = 0 ;
      break ;
    case 25:			/* C-y yank buffer */
      cr = graphPasteBuffer() ;
      if (cr != NULL)					    /* paste can fail. */
	{
	if (ebox->insertMode)
	  while (*cr && strlen(ebox->text) < ebox->len)
	    { cq = ebox->text + strlen(ebox->text) ;
	    while (--cq >= cp) 
	      *(cq + 1) = *cq ;
	    *cp++ = *cr++ ;
	    }
	else
	  while (*cr && cp < ebox->text + ebox->len)
	    *cp++ = *cr++ ;
	}
      break ;
    case END_KEY:
    case 5:			/* C-e go to last char */
      cp = ebox->text + strlen(ebox->text) ;
      break ;
    default:
      if (!isascii(key) || !isprint(key))
	{ ebox->cp = cp ; 
	  return ; 
	}
      if (ebox->insertMode)
	{ if (strlen(ebox->text) < ebox->len)
	    { cq = ebox->text + strlen(ebox->text) ;
	      while (--cq >= cp) 
		*(cq + 1) = *cq ;
	      *cp++ = key ;
	    }
	}
      else
	if (cp < ebox->text + ebox->len)
	  *cp++ = key ;
    }

  ebox->cp = cp ;
  graphSetCursor (ebox) ;
  graphBoxDraw (ebox->box,-1,-1) ;
}  

static int graphTextEntry2 (GraphCompletionFunc f, char* text, 
			    int len, int winlen,
			    float x, float y, void (*fn)(char*))
{
  EntryBox ebox = 0, previous = 0, old = 0 ;
  int n = (int)x ;
  char *cp ;

  if (len < 0)
    messcrash ("Negative length in graphTextEntry") ;
  if (!text)
    messcrash ("Null text in graphTextEntry") ;
  cp = text - 1 ;  /* Clean the buffer */
  while (*++cp) ;
  while (cp < text + len)
    *cp++ = 0 ;

  graphRegister (KEYBOARD, entryBoxEntry) ;

  if (graphAssFind (&TBOX_HANDLE, &previous))
    { if (previous->text == text)
	{ n = strlen(text) ;
	  if (n > previous->len)
	    messcrash ("graphTextEntry2 text > limit %d :\n %s",
		       previous->len, text) ;
	  for ( ; n < previous->len ; n++)
	    text[n] = 0 ;  /* clean up text */
	  if (previous->cp > text + strlen(text))
	    previous->cp = text + strlen(text) ;

	  graphSetCursor (previous) ;
	  graphBoxDraw (previous->box, BLACK, YELLOW) ;
	  return previous->box ;
	}
      else
	{ graphUnSetCursor (previous) ;
	  graphBoxDraw (previous->box, BLACK, PALEGREEN) ;
	  graphAssRemove (&TBOX_HANDLE) ;
	  old = previous ;
	  for (ebox = old->nxt ; ebox && ebox->text != text ; 
	       ebox = ebox->nxt)
	    old = ebox ; 
	}
    }

  if (old && ebox)			/* merely unlink */
    old->nxt = ebox->nxt ;
  else				/* new box */
    { ebox = (EntryBox) handleAlloc (0, gActive->clearHandle, 
				     sizeof (struct EntryBoxStruct)) ;
      ebox->magic = EBOX_MAGIC ;
      ebox->text = text ;
      ebox->cp = text + strlen(text) ;
      ebox->len = len ;
      ebox->winLen = winlen ;
      if (winlen < len)
	ebox->winText = (char*) handleAlloc (0, gActive->clearHandle,
					     winlen+1) ;
      else
	ebox->winText = text ;
      ebox->winPos = 0 ;
      ebox->action = fn ;
      ebox->completionFunction = f ;
      ebox->insertMode = 1 ;
      
      ebox->box = graphBoxStart () ;
      graphTextPtr (ebox->winText,x,y,winlen) ;
      ebox->cursorBox = graphBoxStart () ;
      graphLine (x, y, x, y+UtextY) ;
      graphBoxEnd () ;		/* cursor */
      graphBoxEnd () ;		/* ebox */

      { Box box = gBoxGet(ebox->box) ;	/* set flag for gLeftDown */
	box->flag |= GRAPH_BOX_ENTRY_FLAG ;
	box = gBoxGet(ebox->cursorBox) ;
	box->flag |= GRAPH_BOX_ENTRY_FLAG ;
      }

      text[len-1] = 0 ;		/* for subsequent safety !!!!!!!!!!!!!*/
    }

  ebox->nxt = previous ;	/* 0 if previous not found */
  graphSetCursor (ebox) ;
  graphBoxDraw (ebox->box, BLACK, YELLOW) ;
  graphAssociate (&TBOX_HANDLE, ebox) ;
 
  return ebox->box ;
}

int graphTextEntry (char* text, int len, float x, float y, void (*fn)(char*))
{ return
    graphTextEntry2 (0, text, len, len, x, y, fn) ;
}

int graphCompletionEntry (GraphCompletionFunc f, char* text, int len, float x, float y, void (*fn)(char*))
{ return
    graphTextEntry2 (f, text, len, len, x, y, fn) ;
}

/* gha added next two procedures
   Richard rewwrote this because of cursor problems
 */
int graphTextScrollEntry (char* text, int len, int wlen, float x, float y, void (*fn)(char*))
{ 
  return
    graphTextEntry2 (0, text, len, wlen, x, y, fn) ;
}

int graphCompScrollEntry (GraphCompletionFunc f, char* text, int len, int wlen, float x, float y, void (*fn)(char*))
{ 
  return
    graphTextEntry2 (f, text, len, wlen, x, y, fn) ;
}

void graphEntryDisable (void)
{
  EntryBox ebox ;

  if (graphAssFind (&TBOX_HANDLE, &ebox))
    { graphUnSetCursor(ebox) ;
      graphBoxDraw (ebox->box, BLACK, PALEGREEN) ;
      graphRegister (KEYBOARD, 0) ;
    }
}

/********* clear box list, stacks etc - also used during initialisation ***/

void graphDeleteContents (Graph_ theGraph)  /* Used in graphDestroy() instead of graphClear() */
{
  EntryBox ebox ;

  if (!theGraph)
    return ;

  stackClear (theGraph->boxstack) ;
  stackClear (theGraph->stack) ;
  assClear (theGraph->buttonAss) ;

  /* clean up things hanging off clearHandle */
  handleDestroy (theGraph->clearHandle) ;
  theGraph->clearHandle = handleHandleCreate (theGraph->handle) ;
  theGraph->editors = 0 ;	/* zero dangling pointers */
  theGraph->boxInfo = 0 ;

  if (graphAssFind (&TBOX_HANDLE, &ebox))
    { graphAssRemove (&TBOX_HANDLE) ;
      graphRegister (KEYBOARD, 0) ;
    }
}

void graphClear ()
{
  if (!gActive)
    return ;

  graphDeleteContents (gActive) ;
  graphWhiteOut () ;

  gActive->nbox = 0 ;
  graphBoxStart () ;
  graphTextHeight (0.0) ;
  graphColor (BLACK) ;
  graphTextFormat (PLAIN_FORMAT) ;
  gActive->isClear = TRUE ;

  graphGoto (0.0,0.0) ;		/* RMD 28/5/92 */
}

/********************************* Low Level Graphics ************/

                        /* checks to increase box size */
#define CHECK4\
  if (x0 < x1) { if (x0 < gBox->x1) gBox->x1 = x0 ;\
		   if (x1 > gBox->x2) gBox->x2 = x1 ; }\
  else         { if (x1 < gBox->x1) gBox->x1 = x1 ;\
		   if (x0 > gBox->x2) gBox->x2 = x0 ; }\
  if (y0 < y1) { if (y0 < gBox->y1) gBox->y1 = y0 ;\
		   if (y1 > gBox->y2) gBox->y2 = y1 ; }\
  else         { if (y1 < gBox->y1) gBox->y1 = y1 ;\
		   if (y0 > gBox->y2) gBox->y2 = y0 ; }
#define CHECK3\
  if (x-r < gBox->x1) gBox->x1 = x-r ;\
  if (x+r > gBox->x2) gBox->x2 = x+r ;\
  if (y-r < gBox->y1) gBox->y1 = y-r ;\
  if (y+r > gBox->y2) gBox->y2 = y+r ;

#define CHECK_TEXT\
  if (x < gBox->x1)         gBox->x1 = x ;\
  if (x + len > gBox->x2)   gBox->x2 = x + len ;\
  if (y < gBox->y1)         gBox->y1 = y ;\
  if (y + UtextY > gBox->y2)     gBox->y2 = y + UtextY ;

void graphLine (float x0, float y0, float x1, float y1)

{ push (gStk, LINE, int) ;
  push (gStk, x0, float) ; push (gStk, y0, float) ;
  push (gStk, x1, float) ; push (gStk, y1, float) ;
  CHECK4
}

void graphRectangle (float x0, float y0, float x1, float y1)
{ push (gStk, RECTANGLE, int) ;
  push (gStk, x0, float) ; push (gStk, y0, float) ;
  push (gStk, x1, float) ; push (gStk, y1, float) ;
  CHECK4
}

void graphFillRectangle (float  x0, float y0, float x1, float y1)
{ push (gStk, FILL_RECTANGLE, int) ;
  push (gStk, x0, float) ; push (gStk, y0, float) ;
  push (gStk, x1, float) ; push (gStk, y1, float) ;
  CHECK4
}

void graphCircle (float x, float y, float r)
{ push (gStk, CIRCLE, int) ;
  push (gStk, x, float) ; push (gStk, y, float) ; push (gStk, r, float) ;
  CHECK3
}

void graphPoint (float x, float y)
{ float r = 0.5 * gActive->pointsize ;  /* needed for check3 */
  push (gStk, POINT, int) ;
  push (gStk, x, float) ; push (gStk, y, float) ;
  CHECK3
}

void graphText (const char *text, float x, float y)
{
  float len ;

  if (!gActive || !text) 
    return ;

  push (gStk, TEXT, int) ;
  push (gStk, x, float) ;
  push (gStk, y, float) ;
  pushText (gStk, text) ;

  len = strlen(text) * UtextX ;
  CHECK_TEXT
}

void graphTextUp (const char *text, float x, float y)
{
  float len ;

  if (!gActive || !text) 
    return ;

  push (gStk, TEXT_UP, int) ;
  push (gStk, x, float) ;
  push (gStk, y, float) ;
  pushText (gStk, text) ;

  len = strlen(text) ;
  if (x < gBox->x1) gBox->x1 = x ;
  if (x + 1 > gBox->x2) gBox->x2 = x + 1 ;
  if (y - len*0.6 < gBox->y1) gBox->y1 = y - len*0.6 ;
  if (y > gBox->y2) gBox->y2 = y ;
}

void graphTextPtr (const char *text, float x, float y, int st_len)
{
  float len ;

  if (!gActive || !text)
    return ;

  push (gStk, TEXT_PTR, int) ;
  push (gStk, x, float) ;
  push (gStk, y, float) ;
  push (gStk, text, const char*) ;

  len = st_len * UtextX ;
  CHECK_TEXT ;
} 

void graphTextPtrPtr (const char **text, float x, float y, int st_len)
{
  float len ;

  if (!text)
    return ;
  push (gStk, TEXT_PTR_PTR, int) ;
  push (gStk, x, float) ;
  push (gStk, y, float) ;
  push (gStk, text, const char**) ;

  len = st_len * UtextX ;
  CHECK_TEXT
} 

void graphTextExternal (const char **text, float x, float y, int * st_len)
{
  float len ;

  if (!text)
    return ;
  push (gStk, TEXT_EXTERNAL, int) ;
  push (gStk, x, float) ;
  push (gStk, y, float) ;
  push (gStk, text, const char**) ;
  push (gStk, st_len, int*) ;

  len = (*st_len) * UtextX ;
  CHECK_TEXT
} 

void graphPixels (char *pixels, int w, int h, int len,
		     float x0, float y0, float x1, float y1)
{
  push (gStk, PIXELS, int) ;
  push (gStk, x0, float) ; push (gStk, y0, float) ;
  push (gStk, pixels, char*) ;
  push (gStk, w, int) ; 
  push (gStk, h, int) ; 
  push (gStk, len, int) ;
  push (gStk, x1, float) ; push (gStk, y1, float) ;

  CHECK4
}

void graphPixelsRaw (char *pixels, int w, int h, int len,
		     float x0, float y0)
{
  float x1, y1 ;

  push (gStk, PIXELS_RAW, int) ;
  push (gStk, x0, float) ; push (gStk, y0, float) ;
  push (gStk, pixels, char*) ;
  push (gStk, w, int) ; 
  push (gStk, h, int) ; 
  push (gStk, len, int) ;

  x1 = x0 + XtoUrel(w) ;
  y1 = y0 + YtoUrel(h) ;
  CHECK4
}

void graphColorSquares (char *colors, float x, float y, int len, int skip, int *tints)
{
  if (len < 0)
    messcrash ("len < 0 in graphColorSquares") ;
  if (!len)
    return ;

  push (gStk, COLOR_SQUARES, int) ;
  push (gStk, x, float) ;
  push (gStk, y, float) ;
  push (gStk, colors, char*) ;
  push (gStk, len, int) ;
  push (gStk, skip, int) ;
  push (gStk, tints, int*) ;

  CHECK_TEXT
}

void graphPolygon (Array pts)
{
  float x, y ;
  int i ;

  push (gStk, POLYGON, int) ;
  push (gStk, arrayMax(pts)/2, int);

  for (i = 0 ; i + 1 < arrayMax(pts) ; i += 2)
    {
      x = arr (pts, i, float) ;
      y = arr (pts, i+1, float) ;
      push (gStk, x, float) ;
      push (gStk, y, float) ;
      if (x < gBox->x1) gBox->x1 = x ;
      if (y < gBox->y1) gBox->y1 = y ;
      if (x > gBox->x2) gBox->x2 = x ;
      if (y > gBox->y2) gBox->y2 = y ;
    }
}

void graphLineSegs (Array pts)
{
  float x, y ;
  int i ;

  push (gStk, LINE_SEGS, int) ;
  push (gStk, arrayMax(pts)/2, int);

  for (i = 0 ; i + 1 < arrayMax(pts) ; i += 2)
    { x = arr (pts, i, float) ;
      y = arr (pts, i+1, float) ;
      push (gStk, x, float) ;
      push (gStk, y, float) ;
      if (x < gBox->x1) gBox->x1 = x ;
      if (y < gBox->y1) gBox->y1 = y ;
      if (x > gBox->x2) gBox->x2 = x ;
      if (y > gBox->y2) gBox->y2 = y ;
    }
}

void graphArc (float x, float y, float r, float ang, float angDiff)
{ push (gStk, ARC, int) ;
  push (gStk, x, float) ; push (gStk, y, float) ; 
  push (gStk, r, float) ; 
  CHECK3
  push (gStk, ang, float) ; 
  push (gStk, angDiff, float) ; 
}

void graphFillArc (float x, float y, float r, float ang, float angDiff)
{ push (gStk, FILL_ARC, int) ;
  push (gStk, x, float) ; push (gStk, y, float) ; 
  push (gStk, r, float) ; 
  CHECK3
  push (gStk, ang, float) ; 
  push (gStk, angDiff, float) ; 
}

/********************** routines to set static properties ************/

float graphLinewidth (float x)
{ float old = gActive->linewidth ;
  if (x >= 0)
    { push (gStk, LINE_WIDTH, int) ;
      push (gStk, x, float) ;
      gActive->linewidth = x ;
    }
  return old ;
}

float graphPointsize (float x)
{float old = gActive->pointsize ;
  if (x >= 0)
    { push (gStk, POINT_SIZE, int) ;
      push (gStk, x, float) ;
      gActive->pointsize = x ;
    }
  return old ;
}

float graphTextHeight (float x)
{ float old = gActive->textheight ;
  float dx, dy ;
  if (x >= 0)
    { push (gStk, TEXT_HEIGHT, int) ;
      push (gStk, x, float) ;
      gActive->textheight = x ;
      graphTextInfo (&gActive->textX, &gActive->textY, &dx, &dy) ;
    }
  return old ;
}

int graphColor (int color)
{ int old = gActive->color ;

  if (color >= 0)
    { push (gStk, COLOR, int) ;
      push (gStk, color, int) ;
      gActive->color = color ;
    }
  return old ;
}

void graphBoxColor (int k, int fcol, int bcol)
{ 
  Box box = gBoxGet (k) ;
  if (fcol >= 0)
    box->fcol = fcol ;
  if (bcol >= 0)
    box->bcol = bcol ;
}

int graphTextFormat (int textFormat)
{ int old = gActive->textFormat ;

  if (textFormat >= 0)
    { push (gStk, TEXT_FORMAT, int) ;
      push (gStk, textFormat, int) ;
      gActive->textFormat = textFormat ;
    }
  return old ;
}

/******** routine to shift a box to a new origin ********/

#ifdef BETTER_BOX_SHIFT
static void adjustParentBounds(Box box)
{
  if (box->id) /* for all boxes other than box0 */
    { Box pBox = gBoxGet(box->parentid) ;

      if (box->x1 < pBox->x1) pBox->x1 = box->x1 ;
      if (box->x2 > pBox->x2) pBox->x2 = box->x2 ;
      if (box->y1 < pBox->y1) pBox->y1 = box->y1 ;
      if (box->y2 > pBox->y2) pBox->y2 = box->y2 ;
      adjustParentBounds (pBox) ;
    }
}
#endif

void graphBoxShift (int kbox, float xbase, float ybase)
{
  Box box = gBoxGet (kbox), subBox, box0 ;
  float dx = xbase - box->x1 ;
  float dy = ybase - box->y1 ;
  int nDeep = 1 ;
  int action, r, *lenp = 0;
  float t = 0 ;
  char* text = 0 ;

  stackCursor (gStk,box->mark) ; /* sets position to mark */
  while (!stackAtEnd (gStk))
    switch (action = stackNext (gStk,int))
      {
      case BOX_END:
	if (!--nDeep)
	  goto redraw ;                        /* exit point */
	break ;
      case BOX_START:
	r = stackNext (gStk, int) ; /* shift its coords */
	subBox = gBoxGet (r) ;
	subBox->x1 += dx ; subBox->x2 += dx ;
	subBox->y1 += dy ; subBox->y2 += dy ;
	++nDeep ;
	break ;
      case COLOR: case TEXT_FORMAT:
	r = stackNext (gStk,int) ; 
	break ;
      case LINE_WIDTH: case TEXT_HEIGHT: case POINT_SIZE:
	t = stackNext (gStk,float) ;
	break ;
      case LINE: case RECTANGLE: case FILL_RECTANGLE:
	stackNext (gStk,float) += dx ;
	stackNext (gStk,float) += dy ;
	stackNext (gStk,float) += dx ;
	stackNext (gStk,float) += dy ;
	break ;
      case PIXELS: case PIXELS_RAW:
	stackNext (gStk,float) += dx ;
	stackNext (gStk,float) += dy ;
	text = stackNext (gStk,char*) ;
	r = stackNext (gStk, int) ;
	r = stackNext (gStk, int) ;
	r = stackNext (gStk, int) ;
	if (action == PIXELS)
	  { stackNext (gStk,float) += dx ;
	    stackNext (gStk,float) += dy ;
	  }
	break ;
      case POLYGON : case LINE_SEGS:
	r = stackNext (gStk, int) ;
	if (r > 2)
	  while (r--)
	    { stackNext (gStk,float) += dx ;
	      stackNext (gStk,float) += dy ;
	    }
	break ;
      case CIRCLE: case POINT: 
      case TEXT: case TEXT_UP: case TEXT_PTR: case TEXT_PTR_PTR:case TEXT_EXTERNAL:
      case COLOR_SQUARES: case FILL_ARC: case ARC:
	stackNext (gStk,float) += dx ;
	stackNext (gStk,float) += dy ;
	switch (action)
	  {
	  case CIRCLE:
	    t = stackNext (gStk,float) ; 
	    break ;
	  case FILL_ARC: case ARC:
	    t = stackNext (gStk,float) ; 
	    t = stackNext (gStk,float) ; 
	    t = stackNext (gStk,float) ; 
	    break ;
	  case POINT: 
	    break ;
	  case TEXT: case TEXT_UP:
	    text = stackNextText (gStk) ;
	    break ;
	  case TEXT_PTR:
	    text = stackNext (gStk,char*) ;
	    break ;
	  case TEXT_PTR_PTR:
	    text = *stackNext (gStk,char**) ;
	    break ;
	  case TEXT_EXTERNAL:
	    text = *stackNext (gStk,char**) ;
	    lenp = stackNext (gStk,int *) ;
	    break ;
	  case COLOR_SQUARES:
	    text = stackNext (gStk,char*) ; /* colors */
	    r = stackNext (gStk, int) ;	/* len */
	    r = stackNext (gStk, int) ;	/* skip */
	    text = (char*) stackNext (gStk,int*) ; /* tints */
	    break ;
	  }
	break ;
      case IMAGE: 
	{ char *gim = stackNext (gStk, char*) ; 
	  gim = gim + 0 ;		/* to suppress compiler warning */
	}
	break ;
     default:
       messout ("Invalid draw action %d when shifting, %s %f %d\n",action, text ? text : "", t, *lenp) ;
      }

redraw:
  box->x1 = xbase ; box->x2 += dx ; box->y1 = ybase ; box->y2 += dy ;

#ifdef BETTER_BOX_SHIFT  /* Change parent box bounds, recursively */
  adjustParentBounds (box) ;
#else			 /* incorrect hack */
  box0 = gBoxGet(0) ;
  if (box->x1 < box0->x1) box0->x1 = box->x1 ;
  if (box->x2 > box0->x2) box0->x2 = box->x2 ;
  if (box->y1 < box0->y1) box0->y1 = box->y1 ;
  if (box->y2 > box0->y2) box0->y2 = box->y2 ;
#endif

  graphClipDraw (uToXabs(box->x1-dx)-1, uToYabs(box->y1-dy)-1, 
		 uToXabs(box->x2-dx)+1, uToYabs(box->y2-dy)+1) ;
  graphBoxDraw (kbox, -1, -1) ;
}

/******** end of file ********/
 
 
 
 
 
