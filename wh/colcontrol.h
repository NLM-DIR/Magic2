/*  File: colcontrol.h
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: colcontrol.h - interface for colcontrol.c
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 21 13:54 1998 (fw)
 * * Oct  7 10:31 1998 (edgrif): Removed ACEDB defines, now some acedb
 *              fields are permanently in the structures.
 *-------------------------------------------------------------------
 */

/* $Id: colcontrol.h,v 1.2 2007/03/21 21:03:14 mieg Exp $  */

#ifndef _COLCONTROL_DEF
#define _COLCONTROL_DEF

#include "acedb.h"

/* ACEDB-GRAPH INTERFACE: This typedef is a hangover from the acedb stuff,   */
/* it is actually defined identically in bs.h. Since void* are guaranteed    */
/* not to lose data we are safe locally defining OBJ.                        */
#ifndef DEFINE_OBJ
#define DEFINE_OBJ
  typedef void* OBJ ;
#endif


/* forward definitions                                                       */
typedef struct MapControlStruct *MAPCONTROL ;
typedef struct InstanceStruct *COLINSTANCE;
typedef struct ProtoStruct *COLPROTO;
typedef struct ViewWindowStruct *VIEWWINDOW;


extern magic_t GRAPH2COLCONTROL_ASSOC;
extern magic_t GRAPH2COLCONTROL_LOCALS_ASSOC;

/* Several of the structures below contain 'acedb only' fields, for the      */
/* moment these have been left in here, if we wanted to we could put them    */
/* into extension records so that they did not appear in this header at all. */
/* This would prevent people messing about with them.                        */
/* I have moved all of the acedb fields to the end of their structures.      */


typedef struct ColControlStruct
{
  magic_t *magic;		/* == &COLCONTROL_MAGIC */
  Array instances; /* Column instances in this display */
  Array maps;      /* Maps in this display */
  MAPCONTROL currentMap;  
                  /* Which one draws the header and gets it menu displayed */
  Graph graph;     /* window we control */
  float graphWidth, graphHeight; /* stats about above */
  float topMargin, halfGraphHeight;
  float realTopMargin; /* for header hiding */
  AC_HANDLE handle; /* storage here */
  Array boxIndex;  /* maps boxes -> column instances */
  Array boxIndex2; /* maps boxes -> objects like SEGs */
  Array conversionRoutines; /* column prototype conversion routines. */
  BOOL  needRuler; /* Display a horizontal ruler on left button in box 0 */
  int lastHeaderBox; /* box no of first box for header */
  int activeBox;
  COLINSTANCE activeInstance; /* column instance it belongs to */
  VIEWWINDOW viewWindow;
  int fromBox;     /* the box it got drawn in */
  BOOL hideHeader ;
  Array mapTransitions; /* keep track of which columns are where */
  Associator bottomBoxes; /* To configure columns from the bottom. */
  MAPCONTROL thumbMap; /* used during drag in locator column */

  /* ACEDB-GRAPH INTERFACE: These fields are for use by acedb apps. ONLY.    */
  /*                                                                         */
  KEY from;						    /* Highlight this key on first drawing */
  KEY activeKey;					    /* inter-map communication */

} *COLCONTROL;



struct MapControlStruct
{
  magic_t *magic ;		/* == &MAPCONTROL_MAGIC */
  char *displayName ;		/* type of map fmap, gmap, etc */
  char *name;
  float min, max ;
  float mag, centre ;
  float fastScroll ;		/* x coord of locator */
  float scaleUnit ; 
  float magConf;    /* Configured magnification */
  BOOL hasProjectionLines ; /* Set to allow lines between scale and Locator. */
  BOOL hasCursor ;          /* Set to display cursor */
  struct 
    { int box ;
      BOOL drawn;
      float fac, offset, halfwidth, x ;
    } thumb ;
  struct
    { int box, pickBox ;
      BOOL drawn;
      int val, resolution ;
      float unit, x ;
      char text[32] ;
    } cursor ;
  COLPROTO *prototypes; /* Column prototypes for this map */
  Array protoArray; /* above collected together and copied */
  COLCONTROL control;  /* control structure for the window this map is in */
  MENUOPT *menu; /* top level menu for this map */
  void *look; /* polymorphic */
  void (*beforeDraw)(MAPCONTROL map);     /* called before cols drawn */
  float (*drawHeader)(MAPCONTROL map);     /* draw it */
  void (*headerPick)(MAPCONTROL map, int box, double x, double y); 
                                          /* pick on header */
  void (*convertTrans)(COLCONTROL control, MAPCONTROL map);
  void (*keyboard)(int k);
  
  BOOL needConvert;		/* Does this map need converting */
  Array permConversionRoutines;         /* All routines which need to be called
				        to convert this map */
  Array transConversionRoutines;
  int topMargin;
  COLOUROPT *buttons;
  COLOUROPT **buttonsAddHere;
  Array menusOnButtons;
  int colour;
  FREEOPT *viewMenu; /* temporary, pending new menus */
  AC_HANDLE handle ;
  BOOL suppressed; /* don't display this map  just now */
  BOOL cambridgeOptions; /* for where we agree to disagree */
  BOOL submenus; /* Only hang menus on boxes if this is TRUE */
  BOOL noButtons; /* set by colcontrol, used by mapdisplay, don't have buttons */

  /* ACEDB-GRAPH INTERFACE: These fields are for use by acedb apps. ONLY.    */
  /*                                                                         */
  KEY key;						    /* key of object displayed by this map */
  KEY viewKey;						    /* View used to make the columns */
  KEY parentKey;					    /* used for ManyMap to stop form losing map on recalculation */
  KEY viewTag;						    /* must be set in view for it to be OK */

} ;



struct InstanceStruct
{
  char *name;
  MAPCONTROL map;
  BOOL displayed;
  COLPROTO proto;
  AC_HANDLE handle;
  void (*drawInit)(COLINSTANCE);
  void (*draw)(COLINSTANCE, float *offset);
  BOOL (*configure)(COLINSTANCE);
  BOOL (*setSelectBox)(COLINSTANCE instance, int box, double x, double y);
  BOOL (*unSelectBox) (COLINSTANCE instance, int box);
  void (*doColour)(COLINSTANCE instance, int box);
  void (*followBox)(COLINSTANCE instance, int box, double x, double y);
  void (*pick)(COLINSTANCE instance, int box, double x, double y);
  void (*conversionRegister)(COLCONTROL control, COLINSTANCE instance);
  void (*configFinal)(COLINSTANCE instance, void *locals,BOOL okay);
  Graph configGraph;
  void *private;

  /* ACEDB-GRAPH INTERFACE: These fields are for use by acedb apps. ONLY.    */
  /*                                                                         */
  void (*save)(COLINSTANCE, OBJ);

} ;


struct ProtoStruct {
  void (*init)(COLPROTO proto, MAPCONTROL map);
  BOOL (*create)(COLINSTANCE instance, OBJ init);
  void (*destroy)(void *p);
  char *name;
  void *specialisation;
  BOOL unique;						    /* If set, only one instance of this 
							       column allowed */ 
  void *(*conRoutine)();				    /* how to convert the map */
  void *convertParams;					    /* params for above */
  char *helpText;					    /* What the column does , what the 
							       configparams do. */

  /* Fields after this are not part of the definition, they're maintenance */
  /* used by the code */
  void **convertResults ;				    /* pointer to pointer to convert results */ 
  KEY key;						    /* result of lexaddkey on name */
} ;


typedef struct { COLPROTO proto;
		 BOOL displayed;
		 char *name;
	       } COLDEFAULT[];


#define MAP2GRAPH(map,x) \
  ((map->mag * ((x) - map->centre)) + map->control->halfGraphHeight + \
                                      map->control->topMargin)
#define GRAPH2MAP(map,x) \
  ((((x) - map->control->halfGraphHeight - map->control->topMargin) / \
    map->mag) + map->centre)

#define MAP2WHOLE(map,x) \
  (1 + map->control->topMargin + ((x) - map->thumb.offset)*map->thumb.fac)

#define WHOLE2MAP(map,x) \
  (map->thumb.offset + ((x) - 1 - map->control->topMargin)/map->thumb.fac)


/******************/

COLCONTROL colControlCreate(BOOL isoldgraph, char *name, char* displayName);
COLCONTROL colControlBasicCreate(BOOL isoldgraph, int type, char *name, 
				 float x, float y, float w, float h);
MAPCONTROL mapControlCreate(COLCONTROL control, void (*finalise)(void *p)) ;
COLCONTROL currentColControl(char *callerFuncName);
MAPCONTROL currentMapControl(void);
Graph controlCreateConfig(COLINSTANCE instance, void *locals, char *text,
			  float width, float height) ; /* il */
void controlAddMap(COLCONTROL control, MAPCONTROL map) ;
void controlDeleteMap(COLCONTROL control, MAPCONTROL map);
void controlDeleteInstance(COLCONTROL control, COLINSTANCE instance);

void controlReadConfig(MAPCONTROL map, COLDEFAULT init);
void controlDestroy(COLCONTROL control);
BOOL controlSetColByName(MAPCONTROL map, char *name, int mode);
void controlDraw(void);
void controlOverlay(void);
void controlDrawControl(void *control);
void controlConvertRegister(COLCONTROL control);
void **controlConvert(Array conversionRoutines, MAPCONTROL map, void *(*conRoutine)(), void *params);
void controlConvertPerm(MAPCONTROL map);
void controlConvertTrans(MAPCONTROL map);

void controlAddButton(MAPCONTROL map, COLOUROPT *button, MENUOPT *menu);
float controlDrawButtons(MAPCONTROL map, float x, float y, float max);
							    /* returns the bottom of the bounding box */
void controlUnselectAll(void);
void controlSelectBox(int box);
void controlPrint(void);

void controlRegBox(COLINSTANCE instance, int box, void *private);
void *controlBoxRegd(COLINSTANCE instance, int box);
BOOL controlCheckY(MAPCONTROL map, float coord, float slop, float *ret);
void controlTruncatedLine(COLCONTROL control, float x1, float y1,
			  float x2, float y2);

void mapControlCursorSet (MAPCONTROL map, float x) ;
void mapControlCursorShift(MAPCONTROL map) ;
#define mapControlCursorPos(map) ((map)->cursor.val * (map)->cursor.unit)

extern struct ProtoStruct mapLocatorColumn;
extern struct ProtoStruct mapScrollerColumn;
extern struct ProtoStruct mapScaleColumn;
extern struct ProtoStruct spacerColumn;


#endif  /* _COLCONTROL_DEF */
