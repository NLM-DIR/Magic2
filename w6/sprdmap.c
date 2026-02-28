/*  File: sprdmap.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 **  The aim of this tool is to display on a unique graph as many
 **  system of maps as control by sprdctrl.c
 * Exported functions:
 * HISTORY:
 * Last edited: May 12 10:32 1999 (edgrif)
 * * May 12 10:32 1999 (edgrif): Remove isBlocked, redundant.
 * * Oct 15 16:44 1998 (edgrif): Change names of MAP2GRAPH & GRAPH2MAP
 *              as they clash with the graph headers versions.
 * Created: Wed May 27 11:40:21 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: sprdmap.c,v 1.3 2008/10/20 22:26:50 mieg Exp $ */

#include "acedb.h"

#include "display.h"
#include "spreaddisp_.h"

#include "lex.h"
#include "bs.h"
#include "tags.h"
#include "sysclass.h"
#include "classes.h"
#include "systags.h"

#include "query.h"
#if defined(applec)
#include "Math.h"
#endif

/************************************************************/

typedef struct segstruct SEG ;
struct segstruct
  { KEY key ;
    int line ;
    COL *c ;
    int x, xs ; float y, ys ;    
    int flag ;
    SEG *friends ;  /* Loop of related segs */
    int ibox ;
    KEY parent, grandParent ;
  } ;

typedef struct
  { SEG* seg ;		/* pointer OK: segs fixed when drawing */
    int  color ;
    float xs, ys ;
  } BOX ;


/* p1, p2 are used for contigs: spreadMap <-> pmap
                   for in_situ data for in_situ clones
*/

 
static void spreadMapPick (int box, double x , double y) ;
static void spreadMapDragCursor (float *x, float *y, BOOL isDone) ;
static void spreadMapMiddleDown (double x, double y) ;
static void spreadMapDrag (double x, double y) ;
static void spreadMapUp (double x, double y) ;
static void spreadMapConvert (SPREAD spread) ;
static void spreadMapDraw (SPREAD spread, KEY key) ;
static BOOL spreadMapFindLimits(SPREAD spread) ;
static void spreadMapResize (void) ;
static void spreadMapSelect (SPREAD spread, int box) ;
static void spreadMapFollow (SPREAD spread, double x, double y) ;
static void drawChromosomeBox (SPREAD spread) ;

static void spreadMapHighlight (void) ;
static void spreadMapWhole (void) ;

static void spreadMapChangeSymbolSize (void) ;

 /* void mapPrint(void) ; wrong sinc ei don t use the graph package directly */
static void myMapPrint(void) ;

static MENUOPT spreadMapMenu[] =
              { 
		{graphDestroy, "Quit"},
		{help,"Help"},
		{graphPrint, "Print Screen"},
		{myMapPrint, "Print Whole Map"},
		{displayPreserve,"Preserve"},
		{spreadMapHighlight,"Highlight Selected Objects"},
		{spreadMapChangeSymbolSize, "Symbol size"},
		 {0, 0}
	      } ;

#define FLAG_HIGHLIGHT 1
#define FLAG_INSITU 2
#define FLAG_NAME 4


/*****************************************/

/*****************************************/
/*********** action routines *************/
/**********************************/

static void spreadMapResize (void)
{
  SPREAD spread = currentSpread("spreadMapResize") ;

  spreadMapDraw (spread, 0) ;
}

/***********************************/

static void spreadMapPick (int box, double x, double y) 
{
  int i ;
  COL *c ;
  SPREAD spread = currentSpread("spreadMapPick") ;

  if (!box)
    return ;

  if (box == spread->cursorBox)
    graphBoxDrag (spread->cursorBox, spreadMapDragCursor) ;

  i = arrayMax(spread->colonnes) ;
  while(i--)
    { c = arrp(spread->colonnes,i,COL) ;
      if (box == c->subTitleBox ||
	  box == c->mapColBox)
	{
	  if (c != spread->activeMap)
	    { if (spread->activeMap->subTitleBox)
		graphBoxDraw(spread->activeMap->subTitleBox, BLACK, WHITE) ;
	      spread->activeMap = c ;
	      if (spread->activeMap->subTitleBox)
		graphBoxDraw(spread->activeMap->subTitleBox, BLACK, LIGHTGREEN) ;
	      drawChromosomeBox(spread) ;
	    }
	  return ;
	}
    }
				
  
  if (box == spread->activeMapBox)
    spreadMapFollow (spread, x, y) ;
  else
    spreadMapSelect (spread, box) ;
} /* spreadMapPick */

/*********************/

static void spreadMapClear (void)
{
  KEY curr ;
  SPREAD spread = currentSpread("spreadMapClear") ;
  
  if (spread->activeMapBox && arrp(spread->mapBoxes,spread->activeMapBox,BOX)->seg)
    curr = arrp(spread->mapBoxes,spread->activeMapBox,BOX)->seg->key ;
  else 
    curr = 0 ;
  
  spreadMapDraw (spread, curr) ;
  
  return;
} /* spreadMapClear */


static void spreadMapSelect (SPREAD spread, int box) 
{
  COL *c ;
  SEG *seg, *friend ;
  
  if (spread->activeMapBox)
    { graphBoxDraw (spread->activeMapBox, BLACK, LIGHTBLUE) ;
      seg = friend = arrp(spread->mapBoxes, spread->activeMapBox, BOX)->seg ;
      while (friend = friend->friends, friend != seg)
	if (friend->ibox)
	  graphBoxDraw (friend->ibox, BLACK, LIGHTBLUE) ;
    }

  if ((seg = arrp(spread->mapBoxes, box, BOX)->seg))
    { spread->activeMapBox = box ;
      graphBoxDraw (spread->activeMapBox, BLACK, RED) ;
      seg = friend = arrp(spread->mapBoxes, spread->activeMapBox, BOX)->seg ;
      while (friend = friend->friends, friend != seg)
	if (friend->ibox)
	  graphBoxDraw (friend->ibox, BLACK, LIGHTRED) ;

      c = seg->c ;
      if (c != spread->activeMap)
	{ if (spread->activeMap->subTitleBox)
	    graphBoxDraw(spread->activeMap->subTitleBox, BLACK, WHITE) ;
	  spread->activeMap = c ;
	  if (spread->activeMap->subTitleBox)
	    graphBoxDraw(spread->activeMap->subTitleBox, BLACK, LIGHTGREEN) ;
	  drawChromosomeBox(spread) ;
	}
    
      graphActivate(spread->graph) ;
      spreadSelectFromMap(spread, seg->line, c->colonne) ;
      graphActivate(spread->mapGraph) ;
    }
  else
    spread->activeMapBox = 0 ;

  return;
} /* spreadMapSelect */


static void spreadMapFollow (SPREAD spread, double x, double y)
{
  SEG *seg ;

  seg = arrp(spread->mapBoxes, spread->activeMapBox, BOX)->seg ;
  if (seg && seg->parent)
    {
      if (seg->grandParent)
	{
	  if (display (seg->parent, seg->grandParent, 0))
	    displayPreserve() ;   /* automatically preserve the called map */
	}
      else
	display (seg->parent, 0, TREE) ;
    }
  return;
} /* spreadMapFollow */

/*****************************************/

void spreadMapCreate(void)
{
  SPREAD spread = currentSpread("spreadMapCreate") ;
	     
  if (!spreadMapFindLimits(spread))
    { messout("No visible numerical column in this spread sheet") ;
      return ;
    }
  
  if (graphExists(spread->mapGraph))
    { graphActivate(spread->mapGraph) ;
      graphPop() ;
    }
  else 
    { if (!displayCreate(DtMULTIMAP))
	return ;
      spread->mapGraph = graphActive() ;
      spread->zoomAll =  TRUE ;
      
      graphAssociate (&GRAPH2SPREAD_ASSOC, spread) ;
      graphRegister (RESIZE,spreadMapResize) ;
      graphRegister (PICK, spreadMapPick) ;
      graphRegister (MIDDLE_DOWN, spreadMapMiddleDown) ;
      graphRegister (MIDDLE_DRAG,  spreadMapDrag) ;
      graphRegister (MIDDLE_UP,  spreadMapUp) ;
      graphMenu (spreadMapMenu) ;
      if (*spread->titleBuffer)
	graphText(spread->titleBuffer,32, 1.5) ;
      graphRetitle(spread->titleBuffer) ;
    }

  
  spreadMapConvert(spread) ;
  spreadMapWhole () ;	/* sets ->centre, ->mag and calls Draw() */

  return;
} /* spreadMapCreate */


/**************************************************/
/**************** drawing info ********************/

static int segpOrder(const void *a, const void *b)
{ const SEG 
    *seg1 = *(const SEG**)a ,
    *seg2 = *(const SEG**)b ;
  float 
    ya = seg1->ys , yb = seg2->ys ;
  
  if (seg1->c && seg2->c)
    return ya > yb ? 1 : ( ya == yb ? 0 : -1 ) ;
  if (seg1->c) /* !seg2-> implied */
    return -1 ; 
  if (seg2->c)
    return 1 ;
  return seg1 < seg2 ? -1 : 1 ;
} /* segpOrder */

/***********************************************/

static int	nx, ny ;	/* window dimensions */
static float	yCentre=0 ;	/* ny/2 */
static float	yLength=0 ;	/* length of picture */
static float	topMargin = 6 ;	/* space at top for buttons etc */
static float	bottomMargin = 1 ; /* space at bottom */
static float	xCursor = 5 ;	/* midline of mini-chromosome */
static float	xScale = 10 ;	/* scale bar text LHS */

static float	symbolSize = 3.0 ; /* size of little squares */

static BOOL isPapDisp = FALSE ;

#define SPRDMAP2GRAPH(c, x) \
  (yCentre  +  c->mag * ((float)(x) - c->centre))
#define SPRDGRAPH2MAP(c,x) \
  (((x) - yCentre) / c->mag + c->centre)

/********* start off with some utility routines ********/

static void getNxNy(void)
{
  graphFitBounds (&nx, &ny) ;
  yLength = (ny - topMargin - bottomMargin) ;
  yCentre = topMargin + 0.5*yLength ;
}

/***************************************/

static void spreadMapWhole (void)
{ float  l=0 ;
  COL *c ;
  int j, maxCol ;
  SPREAD spread = currentSpread("spreadMapWhole") ; 

  getNxNy() ;
  if (!spread->activeMap || spread->zoomAll)
    { maxCol = arrayMax(spread->colonnes) ;
      for(j = 0 ; j < maxCol ; j++)
	{ c = arrp(spread->colonnes,j, COL) ;
	  if (!j || c->map)
	    { 
	      l = c->max - c->min ;
	      c->centre = c->min + 0.5 * l ;
	      if (!l) l = 1 ;
	      c->mag = .9 * yLength / l ;
	      spread->activeMap = c ;
	    }
	}
    }
  else
    { c = spread->activeMap ;
      if (c && c->mag)  /* will work on nmae collone */
	{ 
	  l = c->max - c->min ;
	  if (!l) l = 1 ;
	  c->centre = c->min + 0.5 * l ;
	  c->mag = .9 * yLength / l ;
	}
    }

  spreadMapDraw (spread, 0) ;

  return;
} /* spreadMapWhole */

static void spreadMapZoomIn (void)
{
  COL *c ;
  int j , maxCol ;
  SPREAD spread = currentSpread("spreadMapZoomIn") ; 

  if (spread->zoomAll)
    { maxCol = arrayMax(spread->colonnes) ;
      for(j = 0 ; j < maxCol ; j++)
	{ c = arrp(spread->colonnes,j, COL) ;
	  if (!j || c->map)
	    c->mag *= 2 ; 
	}
    }
  else
    { c = spread->activeMap ;
      if (c->mag)
	c->mag *= 2 ; 
    }

  spreadMapDraw (spread, 0) ;

  return;
} /* spreadMapZoomIn */


static void spreadMapZoomOut (void)
{
  COL *c ;
  int j , maxCol ;
  SPREAD spread = currentSpread("spreadMapZoomIn") ; 

  if (spread->zoomAll)
    { maxCol = arrayMax(spread->colonnes) ;
      for(j = 0 ; j < maxCol ; j++)
	{ c = arrp(spread->colonnes,j, COL) ;
	  if (!j || c->map)
	    c->mag /= 2 ; 
	}
    }
  else
    { c = spread->activeMap ;
      if (c && c->mag)
	c->mag /= 2 ; 
    }

  spreadMapDraw (spread, 0) ;

  return;
} /* spreadMapZoomOut */


static void spreadMapZoomAll (void)
{ 
  SPREAD spread = currentSpread("spreadMapZoomAll") ; 

  spread->zoomAll =  TRUE ;
  graphBoxDraw(spread->allButton, BLACK, LIGHTBLUE) ;
  graphBoxDraw(spread->allButton + 1, BLACK, WHITE) ;

  return;
} /* spreadMapZoomAll */

static void spreadMapZoomActive (void)
{ 
  SPREAD spread = currentSpread("spreadMapZoomActive") ; 

  spread->zoomAll = FALSE ;
  graphBoxDraw(spread->allButton + 1, BLACK, LIGHTBLUE) ;
  graphBoxDraw(spread->allButton, BLACK, WHITE) ;

  return;
} /* spreadMapZoomActive */

static void spreadMapChangeSymbolSize (void)
{
  if (messPrompt ("Change size of little square symbols to",
		   messprintf ("%f", symbolSize), "fz"))
    freefloat (&symbolSize) ;

  return;
} /* spreadMapChangeSymbolSize */

static void spreadMapFlip (void)
{ 
  COL *c ;
  SEG *seg ;
  int i ;
  float tmp=0 ;
  SPREAD spread = currentSpread("spreadMapFlip") ; 

  c = spread->activeMap ;
  if (!c)
    return ;
  
  c->flip = ! c->flip ;

  i = arrayMax(c->segs) ;
  seg = arrp(c->segs, 0, SEG) - 1 ;
  while(seg++, i--)
    if (seg->c == c)
      seg->y = - seg->y ;
  tmp = c->min ;
  c->min = - c->max ; c->max = - tmp ;
  c->centre = - c->centre ;

  spreadMapDraw (spread, 0) ;

  return;
} /* spreadMapFlip */

static void spreadMapSwitch (void)
{
  COL *c ;
  int j ;
  SPREAD spread = currentSpread("spreadMapFlip") ; 

  c = spread->activeMap ;
  if (!c)
    return ;
  
  j = arrayMax(spread->colonnes) ;
  while(j--)
    if ( c == arrp(spread->colonnes,j ,COL))
      { spread->activeColonne = j ;
	graphActivate(spread->graph) ;
	spreadSwitchColonnes() ;
	break ;
      }

  return;
} /* spreadMapSwitch */

/**************************************************************/
/***************** dragging code - middle button **************/

static double	oldy, oldDy, oldx;
static BOOL	dragFast ;
#define DRAGFASTLIMIT xScale

static void spreadMapDrag (double x, double y) 
{
  if (dragFast)
    { graphXorLine (0, oldy - oldDy, DRAGFASTLIMIT, oldy - oldDy) ;
      graphXorLine (0, oldy + oldDy, DRAGFASTLIMIT, oldy + oldDy) ;
    }
  else
    graphXorLine (DRAGFASTLIMIT, oldy, nx, oldy) ;

  oldy = y ;

  if (dragFast)
    { oldDy *= exp ((x - oldx) / 25.) ;
      oldx = x ;
      graphXorLine (0, y - oldDy, DRAGFASTLIMIT, y - oldDy) ;
      graphXorLine (0, y + oldDy, DRAGFASTLIMIT, y + oldDy) ;
    }
  else
    graphXorLine (DRAGFASTLIMIT, y, nx, y) ;

  return;
} /* spreadMapDrag */


static void spreadMapUp (double x, double y) 
{
  COL *c ;
  float x1=0,x2=0,y1=0,y2=0 ;
  int j, maxCol ;
  SPREAD spread = currentSpread("spreadMapUp") ;

  if (spread->zoomAll)
    { maxCol = arrayMax(spread->colonnes) ;
      for(j = 0 ; j < maxCol ; j++)
	{ c = arrp(spread->colonnes,j, COL) ;
	  if (!j || c->map)
	    {
	      if (dragFast)
		{
		  graphBoxDim (spread->cursorBox, &x1, &y1, &x2, &y2) ;
		  c->mag *= (y2 - y1) / (2. * oldDy) ;
		  c->centre = c->min + 
		    (c->max - c->min) * (y - topMargin) / yLength ;
		}
	      else
		c->centre +=  (y - 0.5 - yCentre) / c->mag ;
	    }
	}
    }
  else
    { c = spread->activeMap ;
      if (dragFast)
	{ graphBoxDim (spread->cursorBox, &x1, &y1, &x2, &y2) ;
	  c->mag *= (y2 - y1) / (2. * oldDy) ;
	  c->centre = c->min + 
	    (c->max - c->min) * (y - topMargin) / yLength ;
	}
      else
	c->centre +=  (y - 0.5 - yCentre) / c->mag ;
    }

  spreadMapDraw (spread, 0) ;
}

static void spreadMapMiddleDown (double x, double y) 
{  
  float x1=0,x2=0,y1=0,y2=0 ;
  SPREAD spread = currentSpread("spreadMapMiddleDown") ;

  getNxNy () ; 

  graphBoxDim (spread->cursorBox, &x1, &y1, &x2, &y2) ;
  oldDy = (y2 - y1) / 2. ;

  dragFast = (x < DRAGFASTLIMIT) ? TRUE : FALSE ;

  if(dragFast)
    { graphXorLine (0, y - oldDy, DRAGFASTLIMIT, y - oldDy) ;
      graphXorLine (0, y + oldDy, DRAGFASTLIMIT, y + oldDy) ;
    }
  else
    graphXorLine (DRAGFASTLIMIT, y, nx, y) ;
   
  oldx = x ;
  oldy = y ;
  graphRegister (MIDDLE_DRAG,  spreadMapDrag) ;
  graphRegister (MIDDLE_UP,  spreadMapUp) ;
}

static void spreadMapDragCursor (float *x, float *y, BOOL isDone)
{
  if (isDone)
    { float x1,y1,x2,y2 ;
      COL *c ;
      int j, maxCol ;
      SPREAD spread = currentSpread("spreadMapDragCursor") ;

      getNxNy() ;
      graphBoxDim (spread->cursorBox, &x1, &y1, &x2, &y2) ;
      if (spread->zoomAll)
	{ maxCol = arrayMax(spread->colonnes) ;
	  for(j = 0 ; j < maxCol ; j++)
	    { c = arrp(spread->colonnes,j, COL) ;
	      if (!j || c->map)
		c->centre = c->min + (c->max - c->min) * 
		  ((*y + 0.5*(y2-y1)) - topMargin) / yLength ;
	    }
	}
      else
	{ c= spread->activeMap ;
	  c->centre = c->min + (c->max - c->min) * 
	    ((*y + 0.5*(y2-y1)) - topMargin) / yLength ;
	}

      spreadMapDraw (spread, 0) ;
    }
  else
    *x = xCursor - 0.5 ;
}

/**************************************************/
/**************************************************/

static void drawScale (COL *c, int xScale1)
{
  float cutoff = 5 / c->mag ;
  float unit = 0.01 ;
  float subunit = 0.001 ;
  float x=0, xx=0, y=0, start=0, end=0 ;

  while (unit < cutoff)
    { unit *= 2 ;
      subunit *= 5 ;
      if (unit >= cutoff)
	break ;
      unit *= 2.5 ;
      if (unit >= cutoff)
	break ;
      unit *= 2 ;
      subunit *= 2 ;
    }

  start = SPRDGRAPH2MAP(c, topMargin) ;
  if (start < c->min)
    start = c->min ;
  end = SPRDGRAPH2MAP(c, ny-bottomMargin) ;
  if (end > c->max)
    end = c->max ;
      
  x = unit * (((start > c->min) ? 1 : 0) + (int)(start/unit)) ;
  while (x <= end)
    { y = SPRDMAP2GRAPH(c, x) ;
	
      xx = c->flip ?  -x :  x ;
      graphLine (xScale1-1.5,y,xScale1-0.5,y) ;
      if(unit >= 1)
	graphText (messprintf ("%-4.0f",xx),xScale1,y-.5) ;
      else if(unit >= .1)
	graphText (messprintf ("%-4.1f",xx),xScale1,y-.5) ;
      else if(unit >= .01)
	graphText (messprintf ("%-4.2f",xx),xScale1,y-.5) ;
      else if(unit >= .001)
	graphText (messprintf ("%-4.3f",xx),xScale1,y-.5) ;
    
      x += unit ;
    }
  
  x = subunit * (((start>=c->min)?1:0) + (int)(start/subunit)) ;
  while (x <= end)
    { y = SPRDMAP2GRAPH(c,x) ;
      graphLine (xScale1-1.0,y,xScale1-0.5,y) ;
      x += subunit ;
    }
 
  graphLine (xScale1-0.5, SPRDMAP2GRAPH(c,start), 
	     xScale1-0.5, SPRDMAP2GRAPH(c,end)) ;
}

#define MAP2CHROM(x) \
  (topMargin + yLength * (x - c->min) / (c->max != c->min ? c->max - c->min : 1 ))

static void drawChromosomeBox (SPREAD spread)
{ COL *c = spread->activeMap ;

  if (!c)
    return ;

  if (spread->chromoBox)
    graphBoxClear(spread->chromoBox) ;
  spread->chromoBox = graphBoxStart() ;
  graphFillRectangle (xCursor - 0.25, topMargin, 
		      xCursor + 0.25, ny - bottomMargin) ;
  spread->cursorBox = graphBoxStart() ;
  arrayp(spread->mapBoxes,spread->cursorBox,BOX)->seg = 0 ;
  if (c->max != c->min)
    graphRectangle (xCursor - 0.5, 
		  MAP2CHROM(SPRDGRAPH2MAP(c, topMargin)), 
		  xCursor + 0.5, 
		  MAP2CHROM(SPRDGRAPH2MAP(c, ny - bottomMargin))) ;
  graphBoxEnd () ;
  graphBoxEnd () ;
  graphBoxDraw (spread->chromoBox,BLACK, WHITE) ;
  graphBoxDraw (spread->cursorBox,DARKGREEN,GREEN) ;
}

static void drawChromosomeLine (SPREAD spread)
{
  int i ;

  spread->chromoBox = 0 ;

/*
  graphColor (DARKGRAY) ;
  graphLine (xCursor, MAP2CHROM(SPRDGRAPH2MAP(c,topMargin)), 
	     xScale-0.5, topMargin ) ;
  graphLine (xCursor, MAP2CHROM(SPRDGRAPH2MAP(c,ny - bottomMargin)),
	     xScale-0.5, ny - bottomMargin) ;
  graphColor (BLACK) ;
*/
  graphBoxStart() ;
  graphTextHeight (0.75) ;
  for (i = 0 ; i <= 10 ; ++i)
    graphText (messprintf ("%3d%%",10*i),
	       1, topMargin + yLength * i / 10.0 - 0.25) ;
  graphTextHeight (0) ;
  graphBoxEnd() ;
}

/***********************************************/

static void spreadMapBoxEnd (int ibox, BOX *box)
{
  graphBoxEnd () ;
  if (box->seg && box->seg->flag & FLAG_HIGHLIGHT)
    graphBoxDraw (ibox, BLACK, MAGENTA) ;
  else
    graphBoxDraw (ibox, BLACK, box->color) ;
}

/***********************************************/

static BOOL spreadMapFindLimits(SPREAD spread)
{ COL *c , *c1 ; int i, j , maxCol ; BOOL result = FALSE ;
  Array t = spread->tableau , lineArray ;
  float y=0 ; SPCELL *v=0 ;
  BOOL first, found ;
  
  if (!arrayExists(t) || !arrayMax(t))
    return FALSE ;
  maxCol = arrayMax(spread->colonnes) ;
  for(j = 0 ; j < maxCol ; j++)
    { c = arrp(spread->colonnes,j, COL) ;
      c->map = ( c->type == 'i' || c->type == 'f' ) &&
	!c->hidden ;
      if (!c->map)
	continue ;
      for (first = TRUE, i = 0 ; i < arrayMax(t) ; i++ )
	{ lineArray = arr(spread->tableau, i, Array) ;
	  v = arrp(lineArray, j,SPCELL) ;
	  y = 0 ;
	  found = FALSE ;
	  if (!v->empty)
	    switch (c->type)
	      {
	      case 'i':
		found = TRUE ;
		y = v->u.i ;
		break ;
	      case 'f':
		found = TRUE ;
		y = v->u.f ;
		break ;
	      }
          if (found)
	    { y = c->flip ? -y : y ;
	    if (first)
	      { first = FALSE ;
	      result = TRUE ;
	      c->min = y ;
	      c->max = y ;
	      }
	    if (y < c->min)
	      c->min = y ;
	    if (y > c->max)
	      c->max = y ;
	    }
	}
    }
  for(j = 0 ; j < maxCol ; j++)
    { c = arrp(spread->colonnes,j, COL) ;
      if(!c->map && !c->hidden)
	{ for(i = maxCol ; i > 0 ; i--)
	    { c1 = arrp(spread->colonnes,i-1, COL) ;
	      if(c1->map)
		{ c->min = c1->min ; c->max = c1->max ; break ; }
	    }
	}
    }
	  
  return result ;
}

/**********/

static MENUOPT buttonOpts[] = {
  {spreadMapWhole, "Whole"}, 
  {spreadMapZoomIn, "Zoom In"},
  {spreadMapZoomOut, "Zoom Out"},
  {spreadMapClear, "Clear"},
  {spreadMapFlip, "Flip"},
  {spreadMapSwitch, "Switch"},
  {0, 0}} ;

/**********/

static void spreadMapConvert (SPREAD spread)
{ SEG *seg, *oldSeg ;
  int  i, j, j1, maxCol, maxLine ;
  float y=0 ;
  
  Array lineArray ;
  
  COL *c ;  SPCELL *v ;

  maxCol = arrayMax(spread->colonnes) ;
  maxLine = arrayExists(spread->tableau) ? arrayMax(spread->tableau) : 0 ; 

  for (j = 0 ; j < maxCol ; j++)
    {
      j1 = arr(spread->pos2col,j, int) ;
      c = arrp(spread->colonnes,j1,COL) ;
      if (!j1 || c->map)
	c->segs = arrayReCreate (c->segs, maxLine, SEG) ;
      else
	arrayDestroy(c->segs) ;
    }
   
  if (arrayExists(spread->tableau))
    { for(i=0 ; i < maxLine ; i++)
	{ lineArray = arr(spread->tableau, i, Array) ;
	  
	  if (arr(spread->flags, i, char))
	    continue ;  /* Useful if some colonnes are hidden */
	  oldSeg = 0 ;
	  for(j = 0 ; j < maxCol ; j++)
	    {
	      j1 = arr(spread->pos2col,j, int) ;
	      c = arrp(spread->colonnes,j1,COL) ;
	      
	      if (!c->map)
		continue ;
	      v = arrp(lineArray, j1,SPCELL) ;
	      if (v->empty)
		continue ;

	      seg = arrayp(c->segs, i, SEG) ;
	      seg->line = i ;
	      seg->c = c ;
	      if (oldSeg)
		{ seg->friends = oldSeg->friends ;
		  oldSeg->friends = seg ;
		}
	      else
		{ seg->friends = seg ;    /* loop ! */
		  oldSeg = seg ;
		}

	      switch (c->type)
		{
		case 'i':
		  y = v->u.i ;
		  break ;
		case 'f':
		  y = v->u.f ;
		  break ;
		}
	      seg->y = c->flip ? -y : y ;
	      seg->parent = v->parent ;
	      seg->grandParent = v->grandParent ;
	    }
                 /* Now the names */
	  c = arrp(spread->colonnes,0 ,COL) ;
	  v = arrp(lineArray, 0,SPCELL) ;
	  seg = arrayp(c->segs, i, SEG) ;
	  
	  seg->key = v->u.k ;
	  seg->parent = v->u.k ;
	  seg->grandParent = 0 ;

	  seg->line = i ;
	  seg->flag = FLAG_NAME ;
	  seg->c = c ;
	  if (oldSeg)
	    { seg->friends = oldSeg->friends ;
	      oldSeg->friends = seg ;
	      seg->y = seg->friends->y ;
	    }
	  else
	    seg->friends = seg ;    /* loop ! */
	}
    }


  for(j = 0 ; j < maxCol ; j++)
    { c = arrp(spread->colonnes,j ,COL) ;
      if (c->segs)
	{ SEG *zs, **zsp ;
	  c->segps = arrayReCreate (c->segps, maxLine, SEG*) ;
	  i = maxLine ;
	  zsp = arrayp(c->segps, i - 1, SEG*) ;
	  zs = arrp(c->segs, i - 1, SEG) ;
	  while (i--)
	    *zsp-- = zs-- ;
	}
    }
}

/**********/

static void spreadMapDrawC(COL *c, Array mapBoxes, float x)
{ int i = arrayMax(c->segs), ibox, oldFormat=0 ;
  SEG *seg = arrp(c->segs, 0, SEG) - 1 ;
  BOX *box ;
  int screenMin = topMargin, screenMax = ny - bottomMargin ;
	
  while(seg++, i--)
    { if (!seg->c) 
	continue ;
      if (seg->ys > screenMin && seg->ys < screenMax + 30)
	{ ibox = graphBoxStart() ;
	  box = arrayp(mapBoxes, ibox, BOX) ;
	  box->seg = seg ;
	  seg->ibox = ibox ;
	  box->color = YELLOW ;
	  
	  if (seg->flag & FLAG_NAME)
	    { box->color = WHITE ;
	      if (iskey(seg->key) == 2)
		oldFormat = graphTextFormat(BOLD) ;
/* 	      graphTextHeight (0.7) ; */
	      graphText (name(seg->key), x + seg->xs, seg->ys) ;
/*	      graphTextHeight (0.0) ; */
	      if (iskey(seg->key) == 2)
		graphTextFormat(oldFormat) ;
	      seg->x = x + seg->xs ;
	    }
	  else
	    { graphLine(x , seg->ys, x + symbolSize, seg->ys ) ;
	      graphLine(x , seg->ys , x , seg->ys + .5) ;
	      
	      seg->x = x + seg->xs ;
	    }
	  spreadMapBoxEnd(ibox, box) ;
	}
      else
	{ seg->x = x ;  /* needed to draw lines */
	  seg->ibox = 0 ;
	}
    }
}

/**********/

static void spreadMapDrawLines(Array segs)
{
  float x=0, y=0 ;
  int i = arrayMax(segs) ;
  int screenMin = topMargin, screenMax = ny - bottomMargin ;
  SEG *friend, *seg = arrp(segs, 0, SEG) - 1 ;
  
  while (seg++, i--)
    { friend = seg ;
      x = seg->x ; y = seg->ys ;
      if (seg->c) /* may be zero if line is hidden */
	while (friend = friend->friends, friend != seg)
	  { if ((y > screenMin && y < screenMax) ||
		(friend->ys > screenMin && friend->ys < screenMax))
	      if (x - friend->x < 26)  /* I hope nearest neighbours */
		graphLine(x, y, friend->x + symbolSize, friend->ys) ;
	    x = friend->x ; y = friend->ys ;
	  }
    }
}

/**********/

static void spreadMapBump (COL *c)
{ int i = arrayMax(c->segs) ;
  SEG **segp, *seg ;
  int screenMin = topMargin ; 
  float h = ace_lower(c->type) == 'k' ? 1 : .01 ;
  segp = arrp(c->segps, 0, SEG*) ;


  while (i--)
    { seg = *segp++ ;
      if (seg->c)
	{ seg->xs = 0 ;
	  if ((seg->flag & FLAG_NAME)
	      && seg->friends)
	    seg->ys = SPRDMAP2GRAPH (seg->friends->c, seg->y) ;
	  else
	    seg->ys = SPRDMAP2GRAPH (seg->c, seg->y) ;
	}
    }
  i = arrayMax(c->segs) ;
  arraySort(c->segps,segpOrder) ;
  
  c->bump = bumpReCreate(c->bump, 1, ny) ; 
  segp = arrp(c->segps, 0, SEG*) ;
  while(i--)
    { seg = *segp++ ;
      if (!(seg->flag & FLAG_NAME) || seg->ys >= screenMin)
      bumpItem(c->bump, 1, h, &seg->xs, &seg->ys) ;
    }
}

/**********/

static void spreadMapDraw (SPREAD spread, KEY curr)
{
  int x, j, maxCol;
  COL *c=0, *lastc=0 ;

  spread->mapBoxes = arrayReCreate (spread->mapBoxes, 64, BOX) ;
  spread->activeMapBox = 0 ;
  if (spread->activeMap && !((COL*)spread->activeMap)->map)
     spread->activeMap = 0 ;

  getNxNy () ;
  if (!isPapDisp)
    { graphClear () ;
      graphColor (BLACK) ;
      if (yLength < 6)
	{ messout ("Sorry, this window is too small for a multi map") ;
	  return ;
	}
       drawChromosomeLine (spread) ;
    }
  
  maxCol = arrayMax(spread->colonnes) ;
  if (!maxCol) return ;

  if (*spread->titleBuffer)
    graphText(spread->titleBuffer,32, 1.5) ;
  
  for(x = xScale, j = 0 ; j < maxCol ; j++)
    { c = arrp(spread->colonnes,arr(spread->pos2col, j, int) ,COL) ;
      if (c->map)
	{ drawScale(c, x) ;
	  lastc = c ;
	  x += 8 ;
	  spreadMapBump(c) ;
	  c->mapColBox = graphBoxStart() ;
	  spreadMapDrawC(c, spread->mapBoxes, x) ;
	  if (*c->subtitleBuffer)
	    { c->subTitleBox = graphBoxStart() ;
	      graphText(c->subtitleBuffer, x, 5) ;
	      graphBoxEnd() ;
	    }
	  graphBoxEnd() ;
	  x += 3 * symbolSize + 2 ;
	}
      if (!spread->activeMap)
	spread->activeMap = c ;
    }

  drawChromosomeBox (spread) ;
  c = lastc ;
  arrp(spread->colonnes, 0, COL)->mag = c->mag ;
  arrp(spread->colonnes, 0, COL)->centre = c->centre ;
  arrp(spread->colonnes, 0, COL)->min = c->min ;
  arrp(spread->colonnes, 0, COL)->max = (c->max > c->min ? c->max : c->max + 1) ;
  c = arrp(spread->colonnes, 0, COL) ;
  x += 4 ;
  spreadMapBump(c) ; 
  spreadMapDrawC(c, spread->mapBoxes, x) ;
  spreadMapDrawLines(c->segs) ;

  spread->allButton = graphButton ("Zoom All", spreadMapZoomAll, 5, 0.5) ; 
  graphButton ("Zoom Active", spreadMapZoomActive, 15, 0.5) ; 
  if (spread->zoomAll)
    { graphBoxDraw(spread->allButton, BLACK, LIGHTBLUE) ;
      graphBoxDraw(spread->allButton + 1, BLACK, WHITE) ;
    }
  else
    { graphBoxDraw(spread->allButton + 1, BLACK, LIGHTBLUE) ;
      graphBoxDraw(spread->allButton, BLACK, WHITE) ;
    }
  graphButtons (buttonOpts, 5, 2.5, nx) ; 
  spread->activeMapBox = 0 ;
  if (spread->activeMap && spread->activeMap->subTitleBox)
    graphBoxDraw(spread->activeMap->subTitleBox, BLACK, LIGHTGREEN) ;

  graphRedraw () ;
}

/*********************************************************/

static void spreadMapHighlight (void) 
{
  messout("Not yet done, sorry") ;
}
    
    
#include "graph_.h"

static void myMapPrint(void)
{ float min=0, max=0 ;
  int  oldh = gActive->h, oldtype = gActive->type ;
  float oldcentre=0 ,  olduh = gActive->uh ;
  COL *c ;
  SPREAD spread = currentSpread("myMapPrintWhole") ; 

  c = spread->activeMap ;
  
  min = c->min ; max = c->max ;
  if (!messPrompt("Please state the zone you wish to print",
	     messprintf("%g   %g", min, max), "ff"))
    return ;

  freefloat(&min) ; freefloat(&max) ;
  if (min >= max)
    return ;
  if (min < c->min)
    min = c->min ;
  if (max > c->max)
    max = c->max ;

  oldcentre = c->centre ;
  c->centre = (max + min) / 2 ;  
 
  gActive->uh = 1.05 *  (max - min) * c->mag + topMargin + 5 ;
  gActive->h = gActive->uh * gActive->yFac ;
    
  spreadMapDraw (spread, 0) ;

  gActive->type = TEXT_SCROLL ;
  graphPrint() ;
  gActive->type = oldtype ;

  gActive->h = oldh ;  gActive->uh = olduh ;
  c->centre = oldcentre ;
  spreadMapDraw (spread, 0) ;
}
