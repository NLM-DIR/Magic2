/*  File: annot.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
    A general tool to annotate
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 28 18:56 1996 (ulrich)
 * Created: Thu Dec  9 00:01:50 1993 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */
#include "acedb.h"
#include "bs.h"
#include "graph.h"
#include "systags.h"
#include "session.h"
#include "bindex.h"
#include "annot.h"
#include "mytime.h"


typedef struct annotStruct { Graph subGraph ; KEY key ; void *magic ; Graph graph ; char remark[1001] ; } *ANNOT ;

static ANNOT annot = 0 ;
static BOOL annotSave () ;

static ANNOTATION *notes = 0 ;


static MENUOPT subAMenu[] =
{
  {graphDestroy , "Hello Quit"},
  /*   {graphHelp, "Help"}, */
  {graphPrint , "Hello Print"},
  {0, 0}
} ;
static MENUOPT annotMenu[] =
{
  {graphDestroy , "Quit"},
  /*   {graphHelp, "Help"}, */
  {graphPrint , "Print"},
  {0, 0}
} ;

static void annotDestroy (void)
{
  ANNOTATION* note = &notes[0] ;
  if (annot)
    {
      graphCheckEditors (annot->graph, 0) ;
      annotSave () ;
      annot->magic = 0 ;
      annot->graph = 0 ;
      messfree (annot) ;
    
      note-- ; while ((++note)->type)
	switch (note->type)
	  {
	  case _Text:
	    messfree (note->u.s) ;
	    break ;
	  case 0:
	    break ;
	  }
    }
}

BOOL annotSubPick (int box, double x, double y)
{
  messout (messprintf("annotSubPick (%d)", box)) ; return TRUE ;
}
BOOL annotSubKBD (int key)
{
  messout (messprintf("annotKBD (%c)", key)) ; return TRUE ;
}

BOOL annotRD (double x, double y)
{
  messout ("annotRD") ;
  return TRUE ;
}


static void annotDisplay (void)
{
  annotDestroy () ;
  annot = (ANNOT) messalloc (sizeof (struct annotStruct)) ;
  annot->magic = &annot ;    
  annot->graph = graphCreate (TEXT_FIT, "Annotator",.2, .2, .4, .6) ;

  annot->subGraph = graphSubCreate (TEXT_FULL_EDIT, 5.0, 5.0, 16.0, 16.0) ;

  graphRegister (PICK, annotSubPick) ;
  graphRegister (KEYBOARD, annotSubKBD) ;
  graphRegister (RIGHT_DOWN, annotRD) ;

  graphMenu (subAMenu) ;
  graphActivate (annot->graph) ;
  graphRegister (DESTROY , annotDestroy) ;
  graphMenu (annotMenu) ;
}

static void annotInit (KEY key) 
{
  char *cp ;
  OBJ Obj = bsCreate (key) ;
  ANNOTATION* note = &notes[0] ;
  annot->key = key ;
  
  note = &notes[0] ; note-- ; 
  while ((++note)->type)
    if (note->tagName && strlen(note->tagName))
      note->tag = str2tag (note->tagName) ;
  
  note = &notes[0] ; note-- ; 
  while ((++note)->type)
    if (note->tag) switch (note->type)
      {
	case _Text:
	  if (!note->u.s)
	    note->u.s = messalloc (note->length + 1) ;
	  memset (note->u.s, 0, note->length) ;
	  if (Obj && bsGetData (Obj, note->tag, _Text, &cp))
	    strncpy (note->u.s, cp, note->length) ;
	  break ;
	case _Int:
	  note->u.i = 0 ;
	  if (Obj) bsGetData (Obj, note->tag, _Int, &(note->u.i)) ;
	  break ;
	case _Unsigned:
	  note->u.i = 0 ;
	  if (Obj && bsFindTag (Obj, note->tag))
	    note->u.i = 1 ;
	  break ;
	case _Float:
	  note->u.f = 0 ;
	  if (Obj) bsGetData (Obj, note->tag, _Float, &(note->u.f)) ;
	  break ;
	case _DateType:
	  note->u.f = 0 ;
	  if (Obj) bsGetData (Obj, note->tag, _DateType, &(note->u.time)) ;
	  break ;
      default:
	note->u.i = bsFindTag (Obj, note->tag) ;
	break ;
      }    
  bsDestroy (Obj) ;
}

static BOOL annotSave ()
{
  ANNOTATION* note = &notes[0] ;
  OBJ Obj = annot ? bsUpdate (annot->key) : 0 ;
  
  if (!Obj)
    return FALSE ;
  sessionGainWriteAccess() ;
  graphCheckEditors (annot->graph, TRUE) ; 
  note-- ; while ((++note)->type)
    if (note->tag) switch (note->type)
      {
      case _Text:
	if (note->u.s && *(note->u.s))
	  bsAddData (Obj, note->tag, _Text, note->u.s) ;
	else if (bsFindTag (Obj, note->tag))
	  bsRemove (Obj) ;
	break ;
      case _Int:
	if (note->u.i)
	  bsAddData (Obj, note->tag, _Int, &(note->u.i)) ;
	else if (bsFindTag (Obj, note->tag))
	  bsRemove (Obj) ;
	break ;
      case _Unsigned:
	if (note->u.i)
	  bsAddTag (Obj, note->tag) ;
	else if (bsFindTag (Obj, note->tag))
	  bsRemove (Obj) ;
	break ;
      case _Float:
	if (note->u.f != 0.0)
	  bsAddData (Obj, note->tag, _Float, &(note->u.f)) ;
	else if (bsFindTag (Obj, note->tag))
	  bsRemove (Obj) ;
	break ;
      case _DateType:
	if (note->u.time)
	  bsAddData (Obj, note->tag, _DateType, &(note->u.time)) ;
	else if (bsFindTag (Obj, note->tag))
	  bsRemove (Obj) ;
	break ;
      default:
	if (note->u.i)
	  bsAddTag (Obj, note->tag) ;
	else if (bsFindTag (Obj, note->tag))
	  bsRemove (Obj) ;
	break ;
      }
  bsSave (Obj) ;
  return TRUE ;
}

static void annotDraw ()
{
  int line = 1 ;
  ANNOTATION* note = &notes[0] ;
  char timeBuf[25] ;

  graphClear () ;
  
  graphButton ("OK", graphDestroy, 30, line) ;
  graphText ("Gene:", 1, line) ;
  graphText (name(annot->key), 8, line++) ; line++ ;
  note-- ; while ((++note)->type)
    if (note->tag) switch (note->type)
      {
      case _Text:
	graphTextEditor (note->title, note->u.s, note->length, note->x, line++, 0) ;
	break ;
      case _Int:
	graphIntEditor (note->title, &(note->u.i), note->x, line++, 0) ;
	break ;
      case _Unsigned:
	graphIntEditor (note->title, &(note->u.i), note->x, line++, 0) ;
	break ;
      case _Float:
	graphFloatEditor (note->title, &(note->u.f), note->x, line++, 0) ;
	break ;
      case _DateType:
	graphText (note->title, note->x, line) ;
	graphText (timeShow (note->u.time, timeBuf, 25), 5 + strlen(note->title), line++) ;
	break ;
      default:
	{ 
	  BOOL bb = note->u.i ;
	  graphToggleEditor (note->title,  &bb, note->x, line++) ;
	  note->u.i = bb ;
	}
	break ;
      }
    else
      graphText (note->title, note->x, line++) ;
  graphRedraw () ;
}

static void subAnnotDraw (void)
{ 
  int nn = 1 ;
  int box = graphBoxStart () ;
  graphColor (RED) ;
  graphText ("subgraph", 2, 2) ;
  
  graphBoxEnd() ;
  graphBoxMenu (box, subAMenu) ;

  graphBoxDraw (box, BLACK, YELLOW) ;

  graphText (messprintf("Toto %d", nn),  8, 8) ;


  graphRedraw () ;
}

static void annotDraw2 ()
{
  int line = 1 ;
  static int nn = 1 ;

  Graph g = graphActive () ;

  graphClear () ;
  graphButton ("OK", graphDestroy, 30, line) ;
  graphText ("Gene:", 1, line) ;

  graphActivate (annot->subGraph) ;
  subAnnotDraw () ;
  graphActivate (g) ;
  graphText (messprintf ("Gene2: %d", nn), 2, line+20) ;
  graphRedraw () ;
}
 
void annotate (KEY key, ANNOTATION *an)
{
  notes = an ;
  if (!annot || !graphActivate (annot->graph)) 
    annotDisplay () ;
  if (annot->key != key)
    annotInit (key) ;
  graphPop () ;
  if (0) annotDraw () ; /* for compiler happiness */ 
  annotDraw2 () ;
}
