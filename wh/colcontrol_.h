/*  File: colcontrol_.h
 *  Author: Simon Kelley (srk@sanger.ac.uk)
 * Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: private header file for colcontrol.c and viewedit.c
 * HISTORY:
 * Last edited: Dec 17 13:57 1998 (fw)
 * * Oct 13 14:25 1998 (edgrif): Removed ACEDB #if, added a few things
 *              needed for acedb/graph interface.
 * Created: Wed 1 Mar 1995 (srk)
 *-------------------------------------------------------------------
 */

/*  $Id: colcontrol_.h,v 1.2 2007/03/21 21:03:14 mieg Exp $  */

#ifndef DEF_COLCONTROL__H
#define DEF_COLCONTROL__H

#include "colcontrol.h"		/* Our public header. */

BOOL instanceExists(MAPCONTROL map, COLPROTO proto);
void viewWindowCreate(void *dummy);
COLINSTANCE controlInstanceCreate(COLPROTO proto,
				  MAPCONTROL map,
				  BOOL displayed,
				  OBJ init,
				  char *name);

/* These were static in viewedit.c but need to be exposed for the acedb-     */
/* graphinterface.                                                           */
VIEWWINDOW currentView (char *callerFuncName) ;
void viewWindowDraw(VIEWWINDOW view) ;


struct ConversionRecord {   /* record of a conversion routine to be called */
  void *(*convertRoutine)(MAPCONTROL map, void *params);
  void *params;
  void *results;
};

#define EDITLEN 60

struct ViewWindowStruct {
  magic_t *magic;		/* == &VIEWWINDOW_MAGIC */
  COLCONTROL control;
  MAPCONTROL currentMap;
  Graph graph, helpGraph;
  Associator boxToArrow;
  Associator boxToInstance;
  Associator boxToProto;
  Associator boxToButton;
  FREEOPT *mapMenu;
  COLPROTO currentProto;
  COLINSTANCE currentInstance;
  int currentInstanceIndex;
  BOOL editInstance;
  AC_HANDLE handle;
  BOOL dirty;
  int suppressBox, submenusBox, cambridgeBox, hideBox;
  MENU menu;
  char buffer[EDITLEN+1];
} ;


/* ACEDB-GRAPH INTERFACE: Moved here from colcontrol.c itself, because we    */
/* need it for graph-acedb interface.                                        */
typedef struct SpacerPriv {
  int colour;
  float width;
} *SPACERPRIV;


#endif /* DEF_COLCONTROL__H */
