/*  File: htiledisp.c
 *  Author: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the AceView genome database package, written by
 *	Danielle et Jean Thierry-Mieg (NCBI) mieg@ncbi.nlm.nih.gov
 *
 * Description: Display of the genetic map
 * Exported functions:
       htileDisplay
 * HISTORY:

 * Created: July 15 2002 (mieg)
 *-------------------------------------------------------------------
 */


/*
#define ARRAY_CHECK
#define MALLOC_CHECK
*/

#include "acedb.h"
#include "bitset.h"

#include "wfiche/gtitle.h"
#include "graph_.h"
#include "hseqdisp.h"
#include "wiggle.h"
#include "cdna.h"
#include "gmap.h"

  /* so signal of 1024  gives 0 */ 
#define FAKEZERO (log((double)1000+256.0)/log ((double)2.0))
#define VERSION 'G'
static void htileWallBoxDrag (float *x, float *y, BOOL isDone) ;
static void htileRZoneBoxDrag (float *x, float *y, BOOL isDone) ;
static char bufMiddleDrag[12] ;
static AC_DB ac_db = 0 ;

#define NSMAX 100
#define SOLEXAVERSION 2
#define SOLEXAMAX 300
#define SOLEXAMAX2 600
#define slxFormatNint 162
#define slxFormatNfloat SOLEXAMAX2

typedef struct HtileStruct *Htile ;

typedef struct TmapStruct {
  magic_t *magic;        /* == Tmap_MAGIC */
  Htile look ;
  AC_HANDLE h ;
  Graph graph ;

  float mag, leftMargin, offset, dragTopMargin ;
  int graphWidth, graphHeight ;
  int a1, a2, min, max ;
  float lnWidth ;
  COLOUROPT bb [SOLEXAMAX] ;

  int   activeBox, minWallBox, maxWallBox, minRZoneBox, maxRZoneBox, bbBox ;
  float dnaOffset ;
  int starBox ;
  char starBuffer [1000] ; /* a place to write the * when i touch an mrna */

  Array segs ;      	/* array of SEG's from the obj */
  Array zones ;      	/* array of RZS from the obj */
  Array solexa ;
  Array genes ;      	/* array of SEG's from the obj */
  Array introns, intronSupport ;
  Array boxIndex ;    /* box -> int (index into look->segs) */
  int lastTrueSeg ;	/* arrayMax(look->segs) after Convert */
  void (*mapDraw) (void) ;
} Tmap ;

 
typedef struct solexaStruct { int a1, flag, signal[SOLEXAMAX2] ; float gaussSignal[SOLEXAMAX2] ; } SLX ;

typedef struct HtileStruct {
  magic_t *magic;        /* == Htile_MAGIC */
  AC_HANDLE h ;
  AC_DB db ;
  BOOL preserve ;
  Tmap *map ;
  int type ; /* 1: gene, 2:mrna, 3:product */
  BOOL onlyNm, onlyTouched, allProbes, showNoSignalProbes, showOnlyExactProbes, showAmbiguousProbes ;
  KEY intMap, key, from, locus, gene, cosmid, tg, mrna, product ;
  int a1, a2 ; /* int map coords of the cosmid object */
  int min, max ;
  unsigned int flag ;
  int maxIntronClone ;
  Array solexaAll ;
  DICT *tissues ;
  Array menu ;
  int rejectAmbiguous ;
  int unique ; /* 0: non_unique, 1: unique, 2: partials */
  int showDot ;
  BOOL hideHeader, hideScale, hideMusic ;
  BOOL index, zlog, showExtrema;
  BOOL isClosed[NSMAX] ;
  BOOL isSolexaClosed[SOLEXAMAX] ;
  BOOL isSmoothed[NSMAX] ;
  BOOL isSolexaSmoothed[SOLEXAMAX] ;
  int solexaStep ;  /* agnostic, will be set by first acces to a wiggle */
  int mxShowNs ;
  int tmFilter, exonicFilter, smoothing, ratio, ends ;
  int filter99, doMask, isMasked ;
  BOOL showGenome, showGenes, showGeneSignal, showCaptureProbes, showMask, showRZones ;
  int romainSmoothing ;
  float zoom ;
  int gaussWidth ;
  Array wallBufs, mask[256] ;
  char geneSearchBuf[301] ;
  char coordSearchBuf[301] ;
  KEYSET activeGenes ;
} HTS ;

static struct HtileStruct htile0 ;
static Htile look0 = &htile0 ;
static float oldLnWidth = .3 ;
static BOOL oldHideHeader = FALSE ;
static BOOL oldHideScale = FALSE ;
static BOOL oldHideMusic = FALSE ;

typedef enum { Tzero, Tgene, Tmrna, Tprobe
	       , MrnaExon, VeryGoodCodingExon, GoodCodingExon, PoorCodingExon, uORFExon, Utr5, Utr3
	       , MrnaIntron, MrnaIntronFuzzy, GAP
	       , Valid5p, V5Capped, V5SL, V5Aggregate, V5stop
	       , Valid3p, V3aataaa, V3variant
	       , MrnaProbe, GENE
} HSTYPE  ;
char *htileType[] = { "Tzero", "Tgene", "Tmrna", "Tprobe"
		     , "MrnaExon", "VeryGoodCodingExon", "GoodCodingExon", "PoorCodingExon", "uORFExon", "Utr5", "Utr3"
		     , "MrnaIntron", "MrnaIntronFuzzy", "GAP"
		     , "Valid5p", "V5Capped", "V5SL", "V5Aggregate", "V5stop"
		     , "Valid3p", "V3aataaa", "V3variant"
		     , "MrnaProbe", "GENE"
} ;

typedef struct maskStruct { int a1, a2 ; } MASK ;

/* ATTENTION KEEP THE COUNTS CORRECT */
#define segFormatNint 25
#define segFormatNfloat (7+2*NSMAX)
#define MOTIF_FILTER 0x1
#define EXONIC_FILTER 0x2
#define NEW_EXONIC_FILTER 0x4
#define INTRONIC_FILTER 0x8
#define INTERGENIC_FILTER 0x10
#define EXPRESSED_PA 0x20
#define EXPRESSED_NUCPAM 0x40
#define EXPRESSED_TOTAL 0x80
#define MASK_FILTER 0x100
#define CAPTURE_FILTER 0x200
typedef struct rzoneStruct
{
  KEY key, intMap ;
  int type ;
  int a1, a2 ;
  float u1, u2, oldx, oldy, dy, x, y ;
}  RZS ;
typedef struct segStruct
{ 
  KEY key, mrna, probe, map ;
  int col, bcol ;
  int a0, a1, a2, b1, b2, x1, x2 , nClo, ln, bumpy ;
  int tissue ; /* dict  (htile->tissues) */
  HSTYPE type, subType ;
  int widgetBox, tableBox, chip ;
  BOOL exactHit ;
  BOOL ambiguous ;
  unsigned int flag ;
  float tm ;
    float a,c,d,b,e,f ;
  float signal[NSMAX] ;
  float gaussSignal[NSMAX] ;
} SEG ;

#define VSURVZZZZ

static int nsAva, nsAvb, nstm, nsGC, nsAT, nsGCm, nsATm ;

static magic_t Htile_MAGIC = "Htile";
static magic_t Tmap_MAGIC = "Tmap";
static BOOL htileConvert (Htile look, BOOL force, ACEOUT ao) ;
static void htileSolexaConvert (Htile look, BOOL force, ACEOUT ao) ;
static BOOL htileSmoothing (Htile look) ;
static void htileGroupAction (KEY key, int box) ;
static void htileGroupActionButton (void *v) ;
static void htileDrawSolexa (Htile look, float offset, float probeOffset, ACEOUT ao) ;
#ifdef BOU_DEF
 static void htileDrawRZone (Htile look, float offset) ;
 static void htileCreateRZone (void) ;
#endif
static void htileDrawProbes (Htile look, float offset, float probeOffset, ACEOUT ao) ;
static void htileGetMask (Htile look, int maskLength) ;
#define BOU_PNS
#include "../wabi/boubou.h"
 
/************************************************************/
/************************************************************/
/************************************************************/
#define TMAP2GRAPH(_tmap,_x) (((_x) - (_tmap)->offset) * (_tmap)->mag + (_tmap)->leftMargin)
#define TGRAPH2MAP(_tmap,_x) ( (((_x - (_tmap)->leftMargin))/(_tmap)->mag)+  (_tmap)->offset)

/************************************************************/

static AC_OBJ htile_key_obj (Htile look, KEY key, AC_HANDLE h)
{
  return ac_get_obj (look->db, className (key), name(key), h) ;
}

/************************************************************/

static Tmap *tmapCurrent (char *caller)
{
  Tmap *map = 0 ;

  if (!graphAssFind (&Tmap_MAGIC,&map))
    messcrash("%s() could not find Tmap on graph", caller);
  if (!map)
    messcrash("%s() received NULL Htile pointer", caller);
  if (map->magic != &Tmap_MAGIC)
    messcrash("%s() received non-magic Tmap pointer", caller);
  
  return map ;
} /* tmapCurrent */

/************************************************************/

static void tmapResize (void)
{
  Tmap *map = tmapCurrent ("tmapResize") ;
  Htile look = (Htile) map->look ;

  pickRememberDisplaySize (DtTiling) ;

  graphFitBounds (&map->graphWidth,&map->graphHeight) ;
  map->leftMargin = 5 ;
  map->mag = ( (map->graphWidth - map->leftMargin -4))/(map->max - map->min) ;
  map->offset = map->min ;

  if (look->from && ! class (look->from))
    {
      float x2 ;
      map->mag *= 8 ; 
      x2 = TGRAPH2MAP (map, map->graphWidth/2.0);
      map->offset += look->from - x2 ;
    }
  map->mapDraw () ;
} /* tmapResize */

/************************************************************/

static void tmapKbd (int k, GraphKeyModifier mod)
{ 
  Tmap *map = tmapCurrent ("tmapKbd") ;
  float x2, x1 ;

  x2 = x1 =  TGRAPH2MAP (map, map->graphWidth/2.0);
  if (map)
    switch (k)
      {
      case LEFT_KEY :
	if (mod & META_MODKEY)
	  map->offset -= map->graphWidth/(1.0 * map->mag) ;
	else if (mod & (SHIFT_MODKEY & CNTL_MODKEY))
	  map->offset -= map->graphWidth/(2.0 * map->mag) ;
	else
	  map->offset -= map->graphWidth/(4.0 * map->mag) ;
	break ;
      case RIGHT_KEY :
	if (mod & META_MODKEY)
	  map->offset += map->graphWidth/(1.0 * map->mag) ;
	else if (mod & (SHIFT_MODKEY & CNTL_MODKEY))
	  map->offset += map->graphWidth/(2.0 * map->mag) ;
	else
	  map->offset += map->graphWidth/(4.0 * map->mag) ;
	break ;
      case UP_KEY :
	if (mod & SHIFT_MODKEY)
	  {
            map->look->zoom *= 2 ;
            look0->zoom = map->look->zoom ;
	  }
	else if (mod & META_MODKEY)
	  {
	    map->mag *= 4 ;
	    x2 = TGRAPH2MAP (map, map->graphWidth/2.0);
	  }
	else if (mod & CNTL_MODKEY)
	  {
	    map->mag *= 2 ;
	    x2 = TGRAPH2MAP (map, map->graphWidth/2.0);
	  }
	else
	  {
	    map->mag *= 1.414 ;
	    x2 = TGRAPH2MAP (map, map->graphWidth/2.0);
	  }
	break ;
      case DOWN_KEY :
	if (mod & SHIFT_MODKEY)
	  {
	    if (map->look->zoom > 0) map->look->zoom /= 2 ;
	    look0->zoom = map->look->zoom ;
	  }
	else if (mod & META_MODKEY)
	  {
	    map->mag /= 4 ;
	    x2 = TGRAPH2MAP (map, map->graphWidth/2.0);
	  }
	else if (mod & CNTL_MODKEY)
	  {
	    map->mag /= 2 ;
	    x2 = TGRAPH2MAP (map, map->graphWidth/2.0);
	  }
	else
	  {
	    map->mag /= 1.414 ;
	    x2 = TGRAPH2MAP (map, map->graphWidth/2.0);
	  }
	break ;
      case PAGE_UP_KEY :
	map->mag *= 4 ;
	x2 = TGRAPH2MAP (map, map->graphWidth/2.0);
	break ;	
      case PAGE_DOWN_KEY :
	map->mag /= 4 ;
	x2 = TGRAPH2MAP (map, map->graphWidth/2.0);
	break ;	
      case HOME_KEY :
	tmapResize () ;
	return ;
      default:
	return ;
      }
  map->offset -= x2 - x1 ;
  map->look->from = 0 ;  
  map->mapDraw () ;
}

/************************************************************/

static void tmapFollow (Tmap *map, int box)
{
  SEG *seg ;
  int iSeg ;
  if (box > 0 &&
      box < arrayMax(map->boxIndex) &&
      (iSeg = arr (map->boxIndex, box, int)) &&
      iSeg > 0 && 
      (seg = arrp(map->segs, iSeg, SEG)))
    display (seg->key, 0, TREE) ;
} /* tmapFollow */

/************************************************************/

static void tmapSelect (Tmap *map, int box)
{
  SEG *seg = 0 ;
  int iSeg ;
  KEY tg ;
  
  if (box > 0 &&
      box < arrayMax(map->boxIndex) && 
      (tg = arr (map->boxIndex, box, int)) &&
      class(tg)
      )
    {
      display (tg, 0, TREE) ;
      return ;
    }


  if (map->activeBox > 0 &&
      map->segs &&
      map->activeBox < arrayMax(map->boxIndex) && 
      (iSeg = arr (map->boxIndex, map->activeBox, int)) &&
      iSeg > 0 &&
      iSeg < arrayMax (map->segs) &&
      (seg = arrp(map->segs, iSeg, SEG)))
    {
      if (seg->widgetBox > 0 &&
	  seg->widgetBox < arrayMax(map->boxIndex))
	graphBoxDraw (seg->widgetBox, seg->col, seg->bcol) ;
      if (seg->tableBox > 0 &&
	  seg->tableBox < arrayMax(map->boxIndex))
	graphBoxDraw (seg->tableBox, BLACK, WHITE) ;
    }
  
  map->activeBox = 0 ; seg = 0 ;
  if (box > 0 &&
      box < arrayMax(map->boxIndex) &&
      (iSeg = arr (map->boxIndex, box, int)) &&
      iSeg > 0 &&
      (seg = arrp(map->segs, iSeg , SEG)))
    { 
      map->activeBox = box ;
      graphBoxDraw (map->activeBox, BLACK, CYAN) ;
      if (seg->widgetBox > 0 &&
	  seg->widgetBox < arrayMax(map->boxIndex))
	graphBoxDraw (seg->widgetBox, BLACK, RED) ;
      if (seg->tableBox > 0 &&
	  seg->tableBox < arrayMax(map->boxIndex))
	graphBoxDraw (seg->tableBox, BLACK, CYAN) ;
     }

  if (seg && map->starBox)
    {
      if (map->activeBox && map->starBox)
	{
	  int x1, x0 = TMAP2GRAPH (map, seg->a1) ;

	  x1 = map->graphWidth ;
	  if (x0 > x1)
	    x0 = 1 ;
	  if (x1 > 999) x1 = 999 ;
	  memset (map->starBuffer, ' ', x1) ;
	  map->starBuffer [x0 - 1] = '*' ;
	  map->starBuffer [x1] = 0 ;
	}
      else
	map->starBuffer [0] = 0 ;

      graphBoxDraw (map->starBox, CERISE, -1) ;
    }
} /* tmapSelect */

/*************************************/
static int boxMiddleDrag = 0 ;
static Tmap *myMapMiddleDrag = 0 ;
static float xMiddleDrag = 0 ;
static BOOL dnaMiddleDrag = FALSE ;
static void tmapMiddleDrag (int box, double x , double y) 
{
  int a1 ;
 
  if (dnaMiddleDrag)
    {
      float dx = myMapMiddleDrag->mag > 1 ? myMapMiddleDrag->mag : 1.0 ;
      a1 = TGRAPH2MAP (myMapMiddleDrag, myMapMiddleDrag->graphWidth/2.0) + (x - myMapMiddleDrag->graphWidth/2.0)/dx ;
    }
  else
    a1 = TGRAPH2MAP (myMapMiddleDrag, x) ;
  graphXorLine (xMiddleDrag, 10, xMiddleDrag, 1000) ; 

  graphBoxClear (boxMiddleDrag) ;
  xMiddleDrag = x ;
  sprintf (bufMiddleDrag, "%d", a1 + myMapMiddleDrag->a1 - 1) ;
  boxMiddleDrag = graphBoxStart () ;
  graphText (bufMiddleDrag, xMiddleDrag + .3, myMapMiddleDrag->dragTopMargin) ;
  graphBoxEnd () ;
  graphBoxDraw (boxMiddleDrag, BLACK, TRANSPARENT) ;

  graphXorLine (x, 10, x, 1000) ;
  return;
} /* tmapDrag */

/*************************************/

static void tmapMiddleDown (int box,  double x , double y) 
{
  Tmap *map = tmapCurrent ("tmapMiddleDown") ;
  int a1 ;
  
  if (y > map->dnaOffset - 1 && y < map->dnaOffset + 1)
    {
      float dx = map->mag > 1 ? map->mag : 1.0 ;
      a1 = TGRAPH2MAP (map, map->graphWidth/2.0) + (x - map->graphWidth/2.0)/dx ;
      dnaMiddleDrag = TRUE ;
    }
  else
    {  
      a1 = TGRAPH2MAP (map, x) ; 
      dnaMiddleDrag = FALSE ;
    }
  sprintf (bufMiddleDrag, "%d", a1 + map->a1 - 1) ;
  myMapMiddleDrag = map ;
  xMiddleDrag = x ;

  boxMiddleDrag = graphBoxStart () ;
  graphText (bufMiddleDrag, x + .3, map->dragTopMargin) ;
  graphBoxEnd () ;
  graphBoxDraw (boxMiddleDrag, BLACK, TRANSPARENT) ;
  graphXorLine (x, 10, x, 1000) ;
  return;
} /* tmapDown */

/*************************************/

static void tmapMiddleUp (int box,  double x , double y) 
{
  float x1 ;
  Tmap *map = myMapMiddleDrag ;

  graphXorLine (xMiddleDrag, 10, xMiddleDrag, 1000) ;
  graphBoxClear (boxMiddleDrag) ;

  if (dnaMiddleDrag) 
    {
      float dx = map->mag > 1 ? map->mag : 1.0 ;
      x1 = (x - map->graphWidth/2.0)/dx ;
    }
  else
    x1 = TGRAPH2MAP (map, x) - TGRAPH2MAP (map, map->graphWidth/2.0) ;
  map->offset += x1 ;
  
  map->look->from = 0 ;  
  map->mapDraw () ;
  return;
} /* tmapMiddleUp */

/*************************************/
static int htileDraggedBox = 0 ;
static void tmapPick (int box,  double x , double y) 
{
  Tmap *map = tmapCurrent ("tmapPick") ;

  if (!box)
    return ;
  if (0 && box <= map->minWallBox && box < map->maxWallBox)
    {
      htileDraggedBox = box ;
      graphBoxDrag (box, htileWallBoxDrag) ;
    }
  else if (box >= map->minRZoneBox && box < map->maxRZoneBox)
    {
      int iZone =  box - map->minRZoneBox ;
      RZS *rzs ;
      if (iZone >= 0 && iZone < arrayMax (map->zones))
	{
	  htileDraggedBox = box ;
	  rzs = arrayp (map->zones, iZone, RZS) ;
	  /* x, y are coords relative to the rzs box, we want to know if the user moves */
	  if (2*x > rzs->u2-rzs->u1) rzs->x = 1 ; else rzs->x = -1 ;
	  graphBoxDrag (box, htileRZoneBoxDrag) ;
	}
    }
  else if (box < arrayMax (map->boxIndex) && arr(map->boxIndex,box,int)) /* a SEG */
    {
      if (box == map->activeBox) /* a second hit - follow it */
	tmapFollow (map, box) ;
      else
	tmapSelect (map, box) ;
    }

  return;
} /* tmapPick */

/************************************************************/

static Tmap *tmapCreate (void *look, AC_HANDLE handle, void (*mapDraw) (void))
{
  Tmap *map = (Tmap *) halloc (sizeof (struct TmapStruct), handle) ;
  map->magic = &Tmap_MAGIC ;
  map->look = look ;
  map->h = handle ; /* do not destroy, belongs to calling routine */
  map->boxIndex = arrayHandleCreate (64,int, handle) ;
  map->segs = arrayHandleCreate (64, SEG, handle) ;
  map->introns = arrayHandleCreate (64, SEG, handle) ;
  map->intronSupport = arrayHandleCreate (64, SEG, handle) ;
  map->activeBox = 0 ;
  map->lastTrueSeg = arrayMax(map->segs) ;
  map->mapDraw = mapDraw ;

  map->lnWidth = oldLnWidth ;

  return map ;
} /* tmapCreate */

/************************************************************/

static void tmapGraphAssociate (Tmap *map)
{
  graphRegister (KEYBOARD, tmapKbd) ;
  graphRegister (PICK, tmapPick) ;
  graphRegister (MIDDLE_DOWN, tmapMiddleDown) ;
  graphRegister (MIDDLE_DRAG, tmapMiddleDrag) ;
  graphRegister (MIDDLE_UP, tmapMiddleUp) ;
  graphAssRemove (&Tmap_MAGIC) ;
  graphAssociate (&Tmap_MAGIC, map) ; /* attach map to graph */
  graphRegister (RESIZE, tmapResize) ;
  map->graph = graphActive() ;

  return ;
} /* tmapCreate */

/************************************************************/
/************************************************************/

static Htile currentHtile (char *caller)
{
  Htile htile;

  if (!graphAssFind (&Htile_MAGIC,&htile))
    messcrash("%s() could not find Htile on graph", caller);
  if (!htile)
    messcrash("%s() received NULL Htile pointer", caller);
  if (htile->magic != &Htile_MAGIC)
    messcrash("%s() received non-magic Htile pointer", caller);
  
  return htile;
} /* currentHtile */


/************************************************************/

static int  htileSegOrder (const void *a, const void *b)  /* for arraySort() call */ 
{
  const SEG *sa = (const SEG*)a, *sb = (const SEG*) b ;
  int n ;

  n = sa->type - sb->type ;
  if (n) return n ;
  if (sa->mrna != sb->mrna)
    {
      n = sa->bumpy - sb->bumpy ;
      if (n) return n ;
      n = sa->a0 - sb->a0 ; 
      if (n) return n ;
      n = lexstrcmp(name (sa->mrna), name (sb->mrna)) ;
      if (n) return n ;
    }
  n = sa->a1 - sb->a1 ;
  if (n) return n ;
  n = sa->a2 - sb->a2 ;
  if (n) return n ;
  return sa->key - sb->key ;
}

/************************************************************/
/***************** Registered routines *********************/

static void htileDestroy (void)
{
  Htile look = currentHtile("htileDestroy") ;

  if (look)
    {
      ac_free (look->h) ;
      
      graphAssRemove (&Htile_MAGIC) ; /* detach look from graph */
      
      look->magic = 0 ;
      ac_free (look) ;
    }

  return;
} /* htileDestroy */

/*************************************************************/

static void htileForceRecalculate(void)
{
  Htile look = currentHtile("htileForceRecalculate") ;

  if (!messQuery ("Are you sure, this may be very slow, it is better to do it in batch mode"))
    return ;
  if (look)
    {
      htileConvert (look, TRUE, 0) ;
      htileSolexaConvert (look, TRUE, 0) ;
      look->map->mapDraw () ;
    }
} /* htileForceRecalculate  */

/*************************************************************/

static void htileExport (void)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao1, ao2 ;
  Htile look = currentHtile("htileExport") ;
  
  if (look && look->map->solexa)
    {
      ao1 = aceOutCreateToFile (messprintf("%s.solexadata", name(look->key)), "w", h) ;
      htileDrawSolexa (look, 0, 0, ao1) ;
    }
  if (look && look->map->segs)
    {
      ao2 = aceOutCreateToFile (messprintf("%s.probedata", name(look->key)), "w", h) ;
      htileDrawProbes (look, 0, 0, ao2) ;
    }
  ac_free (h) ;
} /* htileExport  */

/*************************************************************/
FREEOPT htileWallTypes[] = {
  { 5, "htileWallTypes"},
  { 1, "Cursor"},
  { 2, "Transition_starts_ES"},
  { 3, "Transition_ends_ES"},
  { 4, "Transition_starts_MS"},
  { 5, "Transition_ends_MS"},
  {0, 0}} ;
static void htileDoAddWall (float x, int oldbox)
{
  AC_HANDLE h = ac_new_handle () ;
  vTXT   bfr = vtxtHandleCreate (h) ;
  Htile look = currentHtile("htileAddWall") ;
  KEY mm ;
  int x1 = TGRAPH2MAP(look->map, x) ;
  char *comment = "" ;

  x1 += look->map->a1 - 1 ;
  if (look && graphSelect (&mm,htileWallTypes ))
    {
      vtxtPrintf (bfr, "Locus %s\n", ac_protect(name (look->locus), h)) ;
      if (oldbox)
	vtxtPrintf (bfr, "-D %s\n", array (look->wallBufs, oldbox - look->map->minWallBox, char*)) ;
      vtxtPrintf (bfr, "Tiling_wall %s %s %d %s\n\n"
		  , htileWallTypes[mm].text
		  , name(look->intMap)
		  , x1, comment) ;
      ac_parse (look->db, vtxtPtr (bfr), 0, 0, h) ;

      look->map->mapDraw () ;
    }

  ac_free (h) ;
} /* htileDoAddWall */

/*************************************************************/

static void htileAddWall (void)
{
  htileDoAddWall (graphEventX, 0) ;
} /* htileAddWall */

/*************************************************************/

static void htileAddWallComment (void)
{
  AC_HANDLE h ;
  vTXT bfr ;
  Htile look = currentHtile("htileRecalculate") ;
  AC_OBJ Region ;
  AC_TABLE tbl ;
  int ir ;
  const char *ccp ;
  char *comment ;
  KEY mm = 0 ;

  if (!messPrompt ("Add a comment and stabilize the tmp walls","","t"))
    return ;

  h = ac_new_handle () ;
  Region = htile_key_obj (look, look->locus, h) ;
  bfr = vtxtHandleCreate (h) ;
  comment = strnew (freeword (),h) ;
  Region = htile_key_obj (look, look->locus, h) ;
  tbl = ac_tag_table (Region, "Tiling_wall", h) ;
  vtxtPrintf (bfr, "Locus %s\n", ac_protect (ac_name(Region), h)) ;
  vtxtPrintf (bfr, "-D Tiling_wall tmp\n") ;
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      ccp = ac_table_printable (tbl, ir, 0, "toto") ;
      if (strcmp (ccp, "tmp"))
	continue ;
      graphSelect (&mm,htileWallTypes ) ;
      vtxtPrintf (bfr, "Tiling_wall %s %s %d %s\n"
		  , htileWallTypes[mm].text
		  , ac_table_printable (tbl, ir, 1, "map")
		  , ac_table_int (tbl, ir, 2, 0) 
		  , comment
		  ) ;
    }
  ac_parse (look->db, vtxtPtr (bfr), 0, 0, h) ;

  ac_free (h) ;
  look->map->mapDraw () ;
} /* htileAddWallComment */

/*************************************************************/
/*
static void htileRecalculate(void)
{
  Htile look = currentHtile("htileRecalculate") ;
  
  if (look)
    {
      look->map->segs = arrayReCreate (look->map->segs,64, SEG) ;
      htileConvert (look, FALSE, 0) ;
      look->map->mapDraw () ;
    } 
}
*/
/*************************************************************/
/************************************************************/

static int htileShowSegs (Array segs)
{
  SEG *ss ;
  int ii;

  if (segs)
    for (ii = 0, ss = arrp (segs, 0, SEG) ; ii < arrayMax (segs) ; ss++, ii++)
      {
	printf ("%3d\t%12s\tbump=%d\ta=\t%5d %5d\tb=\t%5d %5d\tx=\t%5d %5d\tnclo=%5d\t%s\n"
		, ii, htileType[ss->type], ss->bumpy, ss->a1, ss->a2, ss->b1, ss->b2, ss->x1, ss->x2
		, ss->nClo, name(ss->mrna)
		) ;
      }
  return 0 ;
} /* htileShowSegs */

/************************************************************/

static BOOL  htileCosmidConvert (Htile look, KEY cosmid, KEYSET knownProbes)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  BOOL ok = FALSE ;
  AC_OBJ Probe, Cosmid = 0 ;
  AC_TABLE tbl, tbl1, tbl2 ;
  int ir, jj, ns, nns, oldns, a1 = 0, a2 = 0, iSeg = 0 ;
  /* int midPoint, midPoint1, midPoint2  */
  const char *ccp, *dna ;
  SEG *seg ;
  PNS *pns ;
  float gRate, cRate, aRate, tRate ;
  int dna1, dna2, nava, navb ;
  double zlog2 = log ((double)2.0), zava, zavb ; 
  double fakezero ; 
  const char *pNam ;

  fakezero = FAKEZERO ;

  Cosmid = htile_key_obj (look, cosmid, h) ;
  dna = ac_obj_dna (Cosmid, h) ;
  dna2 = dna ? strlen(dna) : 0 ;
  iSeg = arrayMax (look->map->segs) ;
  tbl = ac_tag_table (Cosmid, "IntMap", h) ;
  if (tbl)
    {
      dna1 = a1 = ac_table_int (tbl, 0, 1, 0) ; 
      a2 = ac_table_int (tbl, 0, 2, 0) ; 
    }
  else
    { a1 = dna1 = 1 ; a2 = dna2 ; }

  if (a1 > a2) { int a0 = a1 ; a1 = a2 ; a2 = a0 ; }
  if (!look->a1) look->a1 = a1 ;
  if (look->a2 < a2) look->a2 = a2 ;
 
  tbl = ac_tag_table (Cosmid, "Probe_hit", h) ;
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      KEY probe = ac_table_key (tbl, ir, 0 , 0) ;
      if (0 && ! strstr (ac_table_printable (tbl, ir, 0 , "X"), "BBB_"))
	continue ;
      if (! keySetInsert (knownProbes, probe))
	continue ;
      pNam = name(probe) ;
      if (0 && !strcmp (pNam, "BBC_5843_664_528"))
	invokeDebugger () ;
      h1 = ac_new_handle () ;
      Probe = ac_table_obj (tbl, ir, 0 , h1) ;
      if (look->romainSmoothing == 1)
	tbl1 = ac_tag_table (Probe, "Signal_after_image_smoothing", h1) ;
      else if (look->romainSmoothing==2)
	tbl1 = ac_tag_table (Probe, "Signal_eric", h1) ;
      else
	tbl1 = ac_tag_table (Probe, "Signal", h1) ;
      if (tbl1)
	{
	  seg = arrayp (look->map->segs, iSeg++, SEG) ;
	  seg->key = seg->probe = ac_table_key (tbl, ir, 0 , 0) ; 
	  if (!strncmp(ac_key_name(seg->key),"BBB_",4))
	    seg->chip = 1 ;
	  else if (!strncmp(ac_key_name(seg->key),"BBC_",4))
	    seg->chip = 2 ;
	  seg->type = Tprobe ;

	  if (ac_has_tag (Probe, "Single_genome_hit"))
	    seg->ambiguous = 0 ;
	  else if (ac_has_tag (Probe, "near_by"))
	    seg->ambiguous = 1 ;
	  else if (ac_has_tag (Probe, "on_single_chrom"))
	    seg->ambiguous = 2 ;
	  else 
	    seg->ambiguous = 3 ;

	  seg->tm = ac_tag_float (Probe, "TM", 0) ;
	  seg->flag = 0 ;
	  if (keyFindTag (seg->key, str2tag("Capture")))
	    {
	      seg->flag |= CAPTURE_FILTER ;
	      look0->showCaptureProbes = 
		look->showCaptureProbes = TRUE ;
	    }
	  if (
	      keyFindTag (seg->key, str2tag("Forward_exon_hit")) ||
	      keyFindTag (seg->key, str2tag("Reverse_exon_hit"))
	      )
	    seg->flag |= EXONIC_FILTER ;
	  else if (
		   keyFindTag (seg->key, str2tag("Forward_gene_hit")) ||
		   keyFindTag (seg->key, str2tag("Forward_gene_hit"))
		   )
	    seg->flag |= INTRONIC_FILTER ;
	  else
	    seg->flag |= INTERGENIC_FILTER ;
	  if (keyFindTag (seg->key, str2tag("Expressed_pA")))
	    {
	      seg->flag |= EXPRESSED_PA ;
	      if (!(seg->flag & EXONIC_FILTER))
		seg->flag |= NEW_EXONIC_FILTER ;  
	    }
	  if (keyFindTag (seg->key, str2tag("Expressed_total")))
	    {
	      seg->flag |= EXPRESSED_TOTAL ;
	    }
	  if (keyFindTag (seg->key, str2tag("Expressed_Nuc_pA_moins")))
	    {
	      seg->flag |= EXPRESSED_NUCPAM ;
	    }
	  if (seg->flag & NEW_EXONIC_FILTER)
	    seg->flag &= ~(INTRONIC_FILTER | INTERGENIC_FILTER) ;
	  for (ns = 0 ; ns < NSMAX ; ns++)
	    seg->signal[ns] = -1000 ;
	  tbl2 = ac_tag_table (Probe, "IntMap", h1) ;
	  seg->map = ac_table_key (tbl2, 0, 0, 0) ;
	  seg->a1 = ac_table_int (tbl2, 0, 1, 0) ;
	  seg->a2 = ac_table_int (tbl2, 0, 2, seg->a1+50) ;
	  seg->a1 -= look->map->a1 ; 
	  seg->a2 -= look->map->a1 ; 
	  if (1)
	    {
	      int p, ndna, ng, nc, na, nt, ngc, ncg, ncc, ngg, zdna ;

	      p = ndna = na = nt = ng = nc = ngg = ncc = ngc = ncg = zdna = 0 ;
	      if (0) /* density a/t/g/c around the probe */
		{
		  int idna = ac_table_int (tbl2, 0, 1, 0) - dna1, dna3 ;
		  dna3 = idna + 200 ;
		  idna -= 200 ;
		  if (idna < 0) idna = 0 ;
		  for (ccp = dna + idna ; idna < dna2 && idna < dna3 ; idna++, ccp++)
		    { 
		      ndna++ ; 
		      switch (*ccp)
			{
			case 'g':
			case 'G':
			  ng++ ; 
			  break ;
			case 'c':
			case 'C':
			  nc++ ; 
			  break ;
			case 'a':
			case 'A':
			  na++ ; 
			  break ;
			case 't':
			case 'T':
			  nt++ ; 
			  break ;
			}
		    }
		}
	      else /* denstite a/t/g/c inside the probe */
		{
		  const char *motif = ac_tag_printable (Probe, "Motif", "nnnn") ;
		  int i ;

		  for (ccp = motif ; *ccp ; ccp++)
		    { 
		      ndna++ ; 
		      switch (*ccp)
			{
			case 'g':
			case 'G':
			  ng++ ; 
			  if(ccp > motif && (*(ccp-1)=='c' || *(ccp-1)=='C')) ncg += 2 ; 
			  if(ccp > motif && *(ccp-1)!='g' && *(ccp-1)!='G')
			    {
			      for (i = 1; i<12 ; i++)
				if (*(ccp+i) != 'G' && *(ccp+i) != 'g')
				  break ;
			      if (i > ngg) ngg = i ; 
			    }
			  break ;
			case 'c':
			case 'C':
			  nc++ ; 
			  if(ccp > motif && (*(ccp-1)=='g' || *(ccp-1)=='G')) ngc += 2 ; 
			  if(ccp > motif && *(ccp-1)!='c' && *(ccp-1)!='C')
			    {
			      for (i = 1; i<12 ; i++)
				if (*(ccp+i) != 'C' && *(ccp+i) != 'c')
				  break ;
			      if (i > ncc) ncc = i ; 
			    }
			  break ;
			case 'a':
			case 'A':
			  na++ ; 
			  break ;
			case 't':
			case 'T':
			  nt++ ; 
			  break ;
			}
		      switch (*ccp)
			{
			case 'A': p = 0 ; break ;
			case 'G': if(p>0 && (p & 0x1)) {p++ ; if (p>8) zdna++;} else p = 0 ; break ;
			case 'T': case 'C': if(p>-1 && !(p & 0x1)) {p++ ; if (p>8) zdna++;} else p = 1 ; break ;			}

		    }
		}

	    if (ndna == 0) ndna = 1;
	    gRate = ((float) 100.0 * ng)/ndna ; cRate =  ((float)100.0 * nc)/ndna ;
	    aRate = ((float) 100.0 * na)/ndna ; tRate =  ((float) 100.0 * nt)/ndna ;

	    if (1)  /* 2007_11_30 */
	      {	
		/* filtre 1 */
		if (tRate + gRate <= 18 || tRate + gRate > 84)
		  seg->flag |= MOTIF_FILTER ;  
		if (ndna <= 28)
		  seg->flag |= MOTIF_FILTER ;
		
		/* filtre 2 */
		if (zdna > 7)
		  seg->flag |= MOTIF_FILTER ;
		if (ng - nc < -24)
		  seg->flag |= MOTIF_FILTER ;
		if (ncg > 16 || ncg > 18) 
		  seg->flag |= MOTIF_FILTER ;

		/* filtre 3 */
		if (zdna > 4)
		  seg->flag |= MOTIF_FILTER ;
		if (ng == 0)
		  seg->flag |= MOTIF_FILTER ;
		if (aRate + cRate <= 18 || aRate + cRate > 78)
		  seg->flag |= MOTIF_FILTER ;  
		if (tRate + gRate <= 24)
		  seg->flag |= MOTIF_FILTER ;  
		if (ncg > 14) 
		  seg->flag |= MOTIF_FILTER ;

		/* filtre 4 */
		if (gRate - cRate <= -48)
		  seg->flag |= MOTIF_FILTER ;  
		if (cRate > 54)
		  seg->flag |= MOTIF_FILTER ;  
		if (nc == 0 || nc > 24)
		  seg->flag |= MOTIF_FILTER ;

		/* filtre 5 */
		if (gRate <= 6)
		  seg->flag |= MOTIF_FILTER ;  
		if (aRate + cRate > 72)
		  seg->flag |= MOTIF_FILTER ;  

		/* filtre 6 */
		if (ndna <= 31)
		  seg->flag |= MOTIF_FILTER ;  
		if (seg->tm > 84)
		  seg->flag |= MOTIF_FILTER ;  
		if (aRate + tRate <= 18)
		  seg->flag |= MOTIF_FILTER ;  
		if (aRate + cRate <= 24)
		  seg->flag |= MOTIF_FILTER ;  
		if (tRate + gRate > 76)
		  seg->flag |= MOTIF_FILTER ;  

		/* filtre 7 */
		if (aRate + gRate <= 16)
		  seg->flag |= MOTIF_FILTER ;  
		if (aRate - cRate <= -48)
		  seg->flag |= MOTIF_FILTER ;  
		if (aRate - gRate > 32)
		  seg->flag |= MOTIF_FILTER ;  
		if (gRate + cRate > 78)
		  seg->flag |= MOTIF_FILTER ;  
		if (ngc > 16 || ncg > 12)  
		  seg->flag |= MOTIF_FILTER ;  

		/* filtre 8 */
		if (ng - nc <= -18)
		  seg->flag |= MOTIF_FILTER ;  
		if (na - nt > 25)
		  seg->flag |= MOTIF_FILTER ;  
		if (tRate + cRate > 82)
		  seg->flag |= MOTIF_FILTER ;  
		if (gRate - cRate <= -30)
		  seg->flag |= MOTIF_FILTER ;  
	      }
	    }
	  /*
	    midPoint = (look->map->a2 - look->map->a1)/2.0 ;
	    midPoint1 = midPoint - (look->map->a2 - look->map->a1)/20 ;
	    midPoint2 = midPoint + (look->map->a2 - look->map->a1)/20 ;
	  */
	  nava = navb = 0 ; zava = zavb = 0 ;
	  for (jj = 0, oldns = -1, nns = 0 ; tbl1 && jj < tbl1->rows ; jj++)
	    {
	      char *cp1, mybuf[200] ;
	      
	      ccp = ac_table_printable (tbl1, jj, 0, 0) ;
	      if (!ccp) continue ;

	      /* drop .pair from the name of the signal, so it effectivelly becomes optional */
	      strcpy (mybuf, ccp) ;
	      cp1 = strstr (mybuf, ".pair") ;
	      if (cp1) *cp1 = 0 ;
	      ccp = mybuf ;
	      if (ccp)
		{
		  if (look->romainSmoothing==2)
		    {
		      /* sorry for this hack: 
		       *  there was  a mixup in the tramsmission of the romain-smoothed data files 
		       */
		      if (!strcmp (ccp, "92680_532")) 
			ccp =  "92680_635" ;
		      else if (!strcmp (ccp, "92680_635")) 
			ccp =  "92680_532" ;
		      else if (!strcmp (ccp, "84962_532")) 
			ccp =  "84964_532" ;
		      else if (!strcmp (ccp, "84962_635")) 
			ccp =  "84964_635" ;
		    }
		  for (ns = 0, pns = pnsAll ; pns->p ; ns++, pns++)
		    if (!strcmp (pns->nam, ccp))
		      {
			double z1, z = -1000 ;
			switch (look->romainSmoothing)
			  {
			  case 0: 
			    z = ac_table_float (tbl1, jj, 1, 1) ;
			    /* filter extreme values */
			    if (z > 50000)
			      z = -1000 ;
			    else
			      {
				if (strncmp (pns->p, "N_", 2))
				  {
				    if (z < 0) z = 0 ;
				    z += pns->damper ;
				    if (pns->isRna==9 && z < 1024)
				      z = 1024 ;
				    z = log(z)/zlog2 - fakezero ;  /* so signal of 1024  gives 0 */ 
				  }
				else /* directly show the data as provided */
				  z = z + pns->av - fakezero ;
			      }
			    break ;
			  case 1:
			    z1 = ac_table_float (tbl1, jj, 1, 0) ;
			    if (z1 > 0)
			      z = log(z1/256)/zlog2 ;
			    break ;
			  case 2: 
			    z1 = ac_table_float (tbl1, jj, 1, 0) ;
			    if (z1 > 0) z = z1 ; 
			    break ;
			  }

			if (1) /* filterNonStableProbes */
			  {
			    if (ns == oldns && z > -1000 && 
				(seg->signal[ns] > z + 1 ||  seg->signal[ns] < z - 1)
				)
			      { z = seg->signal[ns] = -1000 ; nns = 0 ; }
			  }
			if (z > -1000)
			  {
			    int uuu = 0;
			    if (ns == oldns)
			      {
				seg->signal[ns] = nns * seg->signal[ns] + z ;
				nns++ ;
				seg->signal[ns] /= nns ;
			      }
			    else
			      {
				nns = 1 ;
				seg->signal[ns] = z ;
			      }
			    oldns = ns ;
			    if (isWriteAccess() &&
				(
				 (
				  ns > 0 && (ns & 0x1) && (uuu=2) &&
				  (pns-1)->isRna >= 1 &&  (pns)->isRna == 0 &&
				  (seg->signal[ns - 1] - (pns-1)->av) -  (z - pns->av) > 2
				  ) ||
				 (
				  (pns+1)->p && !(ns & 0x1) && (uuu=1) &&
				  (pns)->isRna >= 1 &&   (pns+1)->isRna == 0 &&
				  seg->signal[ns + 1] + fakezero > (pns+1)->av - 2 &&
				   (z - pns->av) - (seg->signal[ns + 1] - (pns+1)->av) > 2
				  ) 
				 )
				)
			      {
				char *tissue = 0 ;

				if ((pns+1-uuu)->isRna == 2)
				  seg->flag |= EXPRESSED_PA ;
				else
				  {
				    if (strstr ((pns+1-uuu)->nam, "NucPA-"))
				      seg->flag |= EXPRESSED_NUCPAM ;
				    if (strstr ((pns+1-uuu)->nam, "otal"))
				      seg->flag |= EXPRESSED_TOTAL ;
				  }
				if (!(seg->flag & EXONIC_FILTER))
				  seg->flag |= NEW_EXONIC_FILTER ;

				if (strstr (pns->nam, "ES"))
				  tissue = "ES" ;
				else if (strstr (pns->nam, "MS"))
				  tissue = "MS" ;
				else if (strstr (pns->nam, "ERY"))
				  tissue = "ERY" ;

				if (tissue)
				  {
				    char txt[1000] ;

				    if (seg->flag & EXPRESSED_PA)
				      sprintf (txt, "Probe %s\nExpressed_pA %s \"%s\"\n\n"
					       , ac_protect (name (seg->key), h1)
					       , tissue
					       , pns->nam
					       ) ;
				    if (seg->flag & EXPRESSED_NUCPAM)
				      sprintf (txt, "Probe %s\nExpressed_Nuc_pA_moins %s \"%s\"\n\n"
					       , ac_protect (name (seg->key), h1)
					       , tissue
					       , pns->nam
					       ) ;
				    if (seg->flag & EXPRESSED_TOTAL)
				      sprintf (txt, "Probe %s\nExpressed_total %s \"%s\"\n\n"
					       , ac_protect (name (seg->key), h1)
					       , tissue
					       , pns->nam
					       ) ;
				    ac_parse (look->db, txt, 0, 0, h1) ;
				  }
			      }	 

			    if (pns->nam && !strncmp (pns->nam, "G1", 2))
			      {
				if (strstr (pns->nam, "=1")) { nava++ ; zava += z ; }
				if (strstr (pns->nam, "=2")) { navb++ ; zavb += z ; }
				if (strstr (pns->nam, "=3")) { navb++ ; zavb += z ; }
				if (strstr (pns->nam, "=4")) { navb++ ; zavb += z ; }
				if (strstr (pns->nam, "=5")) { navb++ ; zavb += z ; }
			      }
			  }
			break ;
		      }
		}
	    }
	  if (nava) 
	    seg->signal[nsAva] = zava/nava ;
	  if (navb) 
	    seg->signal[nsAvb] = zavb/navb ;
	  if (seg->chip)
	    {
	      seg->signal[nstm] = seg->tm - 72 ;
	      seg->signal[nstm+1] = 0.001 ;
	      
	      seg->signal[nsGC] = gRate - 25;
	      seg->signal[nsGC+1] = cRate - 25;
	      seg->signal[nsGCm] = gRate - 25;
	      seg->signal[nsGCm+1] = cRate - 25;
	      
	      seg->signal[nsAT] = aRate - 25;
	      seg->signal[nsAT+1] = tRate - 25;
	      seg->signal[nsATm] = aRate - 25;
	      seg->signal[nsATm+1] = tRate - 25;
	    }
	  for (ns = 0 ; ns < NSMAX ; ns++)
	    {
	      if (seg->signal[ns] == 0)
		seg->signal[ns] = -1000 ;
	    }
	}
      ac_free (h1) ;
    }
  ok = TRUE ;

  ac_free (h) ;
  return ok ;
} /* htileCosmidConvert */

/************************************************************/

static void htileMrnaConvert (Htile look)
{
  AC_HANDLE h = ac_new_handle () ;
  Array segs = look->map->segs ;
 

  arraySort (segs,  htileSegOrder) ;
  arrayCompress (segs) ;
  ac_free (h) ;
} /* htileMnraConvert */

/************************************************************/

static void htileInit (void)
{
  const char *ccp ;
  int ns ;
  PNS *pns ;

  ccp = "TM" ;
  for (ns = 0, pns = pnsAll ; pns->p ; ns++, pns++)
    if (!strcmp (pns->p, ccp)) break ;
  nstm = ns ;

  ccp = "G" ;
  for (ns = 0, pns = pnsAll ; pns->p ; ns++, pns++)
    if (!strcmp (pns->p, ccp)) break ;
  nsGC = ns ;
  ccp = "G." ;
  for (ns = 0, pns = pnsAll ; pns->p ; ns++, pns++)
    if (!strcmp (pns->p, ccp)) break ;
  nsGCm = ns ;
  ccp = "A" ;
  for (ns = 0, pns = pnsAll ; pns->p ; ns++, pns++)
    if (!strcmp (pns->p, ccp)) break ;
  nsAT = ns ;
  ccp = "A." ;
  for (ns = 0, pns = pnsAll ; pns->p ; ns++, pns++)
    if (!strcmp (pns->p, ccp)) break ;
  nsATm = ns ;

  ccp = "AvG1A" ;
  for (ns = 0, pns = pnsAll ; pns->p ; ns++, pns++)
    if (!strcmp (pns->p, ccp)) break ;
  nsAva = ns ;

  ccp = "AvG1B" ;
  for (ns = 0, pns = pnsAll ; pns->p ; ns++, pns++)
    if (!strcmp (pns->p, ccp)) break ;
  nsAvb = ns ;

  return ;
} /* htileInit */

/************************************************************/

static void htileRemoveMedian (Htile look)
{
  int iSeg, ns ;
  int iSegMax = arrayMax(look->map->segs) ;
  SEG *seg ;
  double fakezero ;

  fakezero = FAKEZERO ;
  if (iSegMax)
    for (iSeg = 0, seg = arrp (look->map->segs, 0, SEG); iSeg < iSegMax ; iSeg++, seg++)
      {
	for (ns = 0 ; ns < NSMAX ; ns++)
	  {
	    if (seg->signal [ns] != -1000)
	      {
		if (pnsAll[ns].isRna==9)
		  seg->signal [ns] -= 10.0 - fakezero ;
		else
		  seg->signal [ns] -= pnsAll[ns].av - fakezero ;
	      }
	  }
      }
} /* htileRemoveMedian */

/************************************************************/

static PGN *pnsPgn (BOOL ratio, PNS *pns)
{
  PGN *pgn ;
  
  for (pgn = (pns->isRna == 2 ? (ratio ? pgnRatioExonic :  pgnExonic) : (ratio ? pgnRatio : pgnAll)) ; pns->nam && pgn->nam ; pgn++)
    if (!strcmp (pgn->nam, pns->nam))
      return pgn ;
  return 0 ;
} /* pnsPgn */

/************************************************************/

static BOOL htileSmoothing (Htile look)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, k, iSeg, sens ;
  double ss[NSMAX] ;
  int NN = 1 ;
  int width = look->gaussWidth ;
  int lastPos, lastiSeg ; 
  int iSegMax = arrayMax (look->map->segs) ;
  double n[NSMAX] ;
  Array aaa[NSMAX] ;
  PNS *pns ;
  PGN *pgn ;
  SEG *seg, *seg1 ;
  double w, sigma = width ;
  double sigma2 = 2 * sigma * sigma ;
  int solexaStep = look->solexaStep  ; 
  switch (look->smoothing)
    {
    case 1:
      NN = 3 * width/solexaStep  ; 
      if (NN < 1) NN = 1 ;
      break ;
    case 2:
    case 3:
    case 4:
      NN = (width/solexaStep - 1) * 2  ;
      if (NN < 1) NN = 1 ;
      break ;
    case 0:
    case 5:
    case 6:
    case 7:
    default:
      NN = 1 ;
      break ;
    }
  if (NN > 1 && look->smoothing == 4)
    for (k = 0, pns = pnsAll ; pns->p && k < NSMAX ; pns++, k++)
      aaa[k] = NN > 1  ? arrayHandleCreate (2 * NN + 1, double, h) : 0 ;

  lastiSeg = lastPos = -1000000 ;
  if (iSegMax)
    for (iSeg = 0, seg = arrp (look->map->segs, 0, SEG); iSeg < iSegMax ; iSeg++, seg++)
      {
	if (look->smoothing > 1 && seg->a1 < lastPos +  width/10 && look->map->mag * (seg->a1 - lastPos) < 1)
	  { /* no need to show too many data points inside the same smoothing window */
	    for (k = 0 ; k < NSMAX ; k++)
	      seg->gaussSignal[k] = -1000 ;
	    continue ;
	  }
	if (look->smoothing == 1 && iSeg - lastiSeg < NN/4 && look->map->mag * (seg->a1 - lastPos) < 1)
	  { /* no need to show too many data points inside the same smoothing window */
	    for (k = 0 ; k < NSMAX ; k++)
	      seg->gaussSignal[k] = -1000 ;
	    continue ;
	  }
	if (messIsInterruptCalled())
	  goto done ;
	memset (ss, 0, sizeof(ss)) ;
	memset (n, 0, sizeof(n)) ;
	lastiSeg = iSeg ;
	lastPos = seg->a1 ;
	if (
	    (look->rejectAmbiguous < seg->ambiguous) ||
	    (look->tmFilter == 1 && seg->tm > 72.2) ||
	    (look->tmFilter == 2 && seg->tm < 75) ||
	    (look->filter99 == 2 && (seg->flag & MOTIF_FILTER)) ||
	    (look->exonicFilter == 1 && !(seg->flag & EXONIC_FILTER)) ||
	    (look->exonicFilter == 2 && !(seg->flag & NEW_EXONIC_FILTER)) ||
	    (look->exonicFilter == 3 && !(seg->flag & (EXONIC_FILTER | NEW_EXONIC_FILTER))) ||
	    (look->exonicFilter == 4 && !(seg->flag & INTRONIC_FILTER)) ||
	    (look->exonicFilter == 5 && !(seg->flag & INTERGENIC_FILTER))
	    )
	  { 
	    for (k = 0 ; k < NSMAX ; k++)
	      {  n[k] = 0 ; ss[k] = -1000 ; seg->gaussSignal[k] = -1000 ; }
	    continue ;
	  }
	for (k = 0, pns = pnsAll ; pns->p && k < NSMAX ; pns++, k++)
	  {
	    if (1 && look->isSmoothed[k])
	      continue ;
	    if (!pns->p)
	      break ;
	    pgn = pnsPgn (look->ratio, pns) ;
	    if (
		(!look->ratio && look->isClosed [k]) ||
		(look->ratio && ! (k & 1) && look->isClosed [k]) ||
		(look->ratio &&   ( k& 1) && look->isClosed [k-1]) ||
		seg->signal[k] == -1000 ||
		(pgn && look->filter99 == 3 && !look->ratio && seg->signal[k] + 10 > pgn->s99) ||
		(pgn && look->filter99 == 3 && look->ratio &&
		 seg->signal[k] - seg->signal[k+1] > pgn->s99
		 )
		)
	      { n[k] = 0 ; ss[k] = -1000 ; }
	    else if (pns->flag & PGG_NOZOOM)
	      {
		arrayMax(aaa[k]) = 0 ; 
		array(aaa[k],0,double) = seg->signal[k] ;
		n[k] = 1 ;
	      } 
	    else
	      {
		switch (look->smoothing)
		  {
		  case 2:
		    w = NN * NN ;
		    n[k] = w ;
		    ss[k] = seg->signal[k] * w ;
		    break ;
		  case 4: 
		    n[k] = 1 ;
		    if (NN > 1)
		      {
			arrayMax(aaa[k]) = 0 ; 
			array(aaa[k],0,double) = seg->signal[k] ;
		      }
		    else
		       seg->gaussSignal[k] = seg->signal[k] ;
		    break ;
		  default:
		    n[k] = 1 ;
		    ss[k] = seg->signal[k] ;
		    break ;
		  }
	      }
	  }
	for (sens = -1 ; sens < 2 ; sens += 2)
	  for (j = 1, i = iSeg + sens, seg1 = seg + sens ; j < NN && i >= 0 && i < arrayMax(look->map->segs)  ;
	       j++, i += sens, seg1 += sens)
	    {
	      if (look->rejectAmbiguous < seg1->ambiguous)
		continue ; 
	      if (
		  (look->tmFilter == 1 && seg1->tm > 72) ||
		  (look->tmFilter == 2 && seg1->tm < 75) ||
		  (look->filter99 == 2  && (seg1->flag & MOTIF_FILTER)) ||
		  (look->exonicFilter == 1 && !(seg->flag & EXONIC_FILTER)) ||
		  (look->exonicFilter == 2 && !(seg->flag & NEW_EXONIC_FILTER)) ||
		  (look->exonicFilter == 3 && !(seg->flag & (EXONIC_FILTER | NEW_EXONIC_FILTER))) ||
		  (look->exonicFilter == 4 && !(seg->flag & INTRONIC_FILTER)) ||
		  (look->exonicFilter == 5 && !(seg->flag & INTERGENIC_FILTER))
		  )
		continue ;
	      
	      if (seg1->map != seg->map)
		continue ;
	      
	      for (k = 0, pns = pnsAll ; pns->p && k < NSMAX ; pns++, k++)
		{
		  if (
		      (!look->ratio && look->isClosed [k]) ||
		      (look->ratio && ! (k & 1) && look->isClosed [k]) ||
		      (look->ratio &&   ( k& 1) && look->isClosed [k-1]) ||
		      (1 && look->isSmoothed[k])
		      )
		    continue ;
		  if (!n[k] || seg1->signal[k]==-1000)
		    continue ;
		  pgn = pnsPgn (look->ratio, pns) ;
		  if (pns->flag & PGG_NOZOOM)
		    continue ;
		  if (pgn && look->filter99 == 3 && pgn && seg1->signal[k] + (look->ratio ? 0 : 10) > pgn->s99)
		    continue ;
		  switch (look->smoothing)
		    {
		    case 1:
		      w = seg->a1 - seg1->a1 ; 
		      w = w * w / sigma2 ;
		      w = exp(-w) ;
		      n[k] += w ;
		      ss[k] += seg1->signal[k] * w ;
		      break ;
		    case 2:
		      w = i - iSeg ;
		      w = NN * NN - w * w ;
		      n[k] += w ;
		      ss[k] += seg1->signal[k] * w ;
		      break ;
		    case 3:
		      n[k]++ ;
		      ss[k] += seg1->signal[k] ;
		      break ;
		    case 4:
		      array(aaa[k],n[k],double) = seg1->signal[k] ;
		      n[k]++ ;
		      break ;
		    }
	      }
	    }
	for (k = 0 ; k < NSMAX ; k++)
	  {
	    if (
		(!look->ratio && look->isClosed [k]) ||
		(look->ratio && ! (k & 1) && look->isClosed [k]) ||
		(look->ratio &&   ( k& 1) && look->isClosed [k-1]) ||
		(1 && look->isSmoothed[k])
		)
	      continue ;

	    if (pnsAll[k].isRna == 92)
	      seg->gaussSignal[k] = seg->signal[k] ;
	    else if (n[k] > 0)
	      {
		if (look->smoothing == 4)
		  {
		    if (NN > 1)
		      {
			arraySort (aaa[k], doubleOrder) ;
			seg->gaussSignal[k] = array (aaa[k], n[k]/2, double) ;
		      }
		    else
		      seg->gaussSignal[k] = seg->signal[k] ;
		  }
		else 
		  seg->gaussSignal[k] = ss[k]/n[k] ;
	      }
	    else
	      seg->gaussSignal[k] = 0 ;
	  }
      }
 
  if (look->ratio && look->smoothing)
    for (k = 0 ; k < NSMAX ; k += look->ratio ? 2 : 1)
      look->isSmoothed[k+1] = look->isSmoothed[k] = !look->isClosed [k] ;
  else  if (!look->ratio && look->smoothing)
    for (k = 0 ; k < NSMAX ; k++)
      look->isSmoothed[k] = !look->isClosed [k] ;

  /*
    for (k = 0 ; k < NSMAX ; k++)
    look->isSmoothed[k] = FALSE ;
  */
 done:
  ac_free (h) ;

  return TRUE ;
} /* htileSmoothing */

/************************************************************/

static BOOL htileSolexaSmoothing (Htile look)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, k, iSlx, sens ;
  double ss[SOLEXAMAX2] ;
  int width = look->gaussWidth ;
  int NN = 1 ;
  int lastPos, lastiSlx ;
  double n[SOLEXAMAX2] ;
  Array mask = 0, aaa[SOLEXAMAX2] ;
  MASK *mm ;
  PNX *pnx ;
  SLX *slx, *slx1 ;
  double w, sigma = width ;
  double sigma2 = 2 * sigma * sigma ;
  int solexaStep = look->solexaStep  ; 

  mask = 0 ;
  if (look->doMask && ! look->isMasked)
    {
      htileGetMask (look, look->doMask) ;
      mask = look->mask[look->doMask] ;
      if (mask == (Array) 1 || !arrayMax(mask)) 
	mask = 0 ;
    }
  if (mask && ! look->isMasked) /* label the probes that should be masked */
    {
      int iMask = 0, iMaskMax = arrayMax (mask) ;
      
      mm = arrp (mask, iMask, MASK) ;
      for (iSlx = 1, slx = arrp (look->map->solexa, 1, SLX); iSlx < arrayMax(look->map->solexa) ; iSlx++, slx++)
	{
	  for ( ; iMask < iMaskMax && mm->a2 < slx->a1 ; mm++, iMask++) ; /* find a candidate zone */
	  if (iMask >= iMaskMax)                                          /* no more masked zone */
	    break ;
	  if (slx->a1 >= mm->a1) 
	    slx->flag |= MASK_FILTER ;               /* probe inside the masked zone */
	}
      look->isMasked = TRUE ;
    }
    
  switch (look->smoothing)
    {
    case 1:
      NN = 3 * width/solexaStep  ; 
      if (NN < 1) NN = 1 ;
      break ;
    case 2:
    case 3:
    case 4:
      NN = (width/solexaStep - 1) * 2  ;
      if (NN < 1) NN = 1 ;
      break ;
    case 0:
    case 5:
    case 6:
    case 7:
    default:
      NN = 1 ;
      break ;
    }
  if (NN > 1 && look->smoothing == 4)
    for (k = 0, pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p &&  k < SOLEXAMAX ; pnx++, k++)
      aaa[k] = NN > 1  ? arrayHandleCreate (2 * NN + 1, double, h) : 0 ;

  lastiSlx = lastPos = -1000000 ;
  for (iSlx = 1, slx = arrp (look->map->solexa, 1, SLX); iSlx < arrayMax(look->map->solexa) ; iSlx++, slx++)
    {
      if (look->smoothing == 0) continue ;
      if (look->doMask && (slx->flag & MASK_FILTER)) continue ;
      if (look->smoothing > 1 && slx->a1 < lastPos +  width/10 && look->map->mag * (slx->a1 - lastPos) < 1)
	{ /* no need to show too many data points inside the same smoothing window */
	  for (k = 0,  pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p &&  k < SOLEXAMAX ; pnx++, k++)
	    slx->gaussSignal[k] = -1000 ;
	  continue ;
	}
      if (look->smoothing == 1 && iSlx - lastiSlx < NN/4 && look->map->mag * (slx->a1 - lastPos) < 1)
	{ /* no need to show too many data points inside the same smoothing window */
	  for (k = 0,  pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p &&  k < SOLEXAMAX ; pnx++, k++)
	    slx->gaussSignal[k] = -1000 ;
	  continue ;
	}
      if (messIsInterruptCalled())
	goto done ;
      lastiSlx = iSlx ;
      if (lastPos == slx->a1)
	continue ;
      lastPos = slx->a1 ;


      if (
	  (look->exonicFilter == 1 && !(slx->flag & EXONIC_FILTER)) ||
	  (look->exonicFilter == 2 && !(slx->flag & NEW_EXONIC_FILTER)) ||
	  (look->exonicFilter == 3 && !(slx->flag & (EXONIC_FILTER | NEW_EXONIC_FILTER))) ||
	  (look->exonicFilter == 4 && !(slx->flag & INTRONIC_FILTER)) ||
	  (look->exonicFilter == 5 && !(slx->flag & INTERGENIC_FILTER))
	  )
	{ 
	  for (k = 0,  pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p &&  k < SOLEXAMAX ; pnx++, k++)
	    {  n[k] = 0 ; ss[k] = -1000 ; slx->gaussSignal[k] = -1000 ; }
	  continue ;
	}
      for (k = 0,  pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p &&  k < SOLEXAMAX ; pnx++, k++)
	{
	  if (1 && look->isSolexaSmoothed[k])
	    continue ;
	  if (pnx->noSmoothing)
	    continue ;
	  if (
	      (!look->ratio && look->isSolexaClosed [k]) ||
	      (look->ratio && ! (k & 1) && look->isSolexaClosed [k]) ||
	      (look->ratio &&   ( k & 1) && look->isSolexaClosed [k-1]) ||
	      (slx->signal[k] == -1000)
	      )
	    { n[k] = 0 ; ss[k] = -1000 ; }
	  else 
	    {
	      switch (look->smoothing)
		  {
		  case 2:
		    w = NN * NN ;
		    n[k] = w ;
		    ss[k] = slx->signal[k] * w ;
		    break ;
		  case 4:
		    n[k] = 1 ;
		    if (NN > 1)
		      {
			arrayMax(aaa[k]) = 0 ; 
			array(aaa[k],0,double) = slx->signal[k] ;
		      }
		    else
		       slx->gaussSignal[k] = slx->signal[k] ;
		    break ;
		  default:
		    n[k] = 1 ;
		    ss[k] = slx->signal[k] ;
		    break ;
		  }
	    }
	}
      for (sens = -1 ; sens < 2 ; sens += 2)
	for (j = 1, i = iSlx + sens, slx1 = slx + sens ; j < NN && i >= 1 && i < arrayMax(look->map->solexa)  ;
	     j++, i += sens, slx1 += sens)
	  {
	    if (look->doMask && (slx1->flag & MASK_FILTER)) continue ;

	    for (k = 0,  pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p &&  k < SOLEXAMAX ; pnx++, k++)
	      {
		if (
		    (!look->ratio && look->isSolexaClosed [k]) ||
		    (look->ratio && ! (k & 1) && look->isSolexaClosed [k]) ||
		    (look->ratio &&   ( k& 1) && look->isSolexaClosed [k-1]) ||
		    (1 && look->isSolexaSmoothed[k])
		    )
		  continue ;
		if (!n[k] || slx1->signal[k]==-1000)
		  continue ;
		switch (look->smoothing)
		  {
		  case 1:
		    w = slx1->a1 - slx->a1 ; 
		    w = w * w / sigma2 ;
		    w = exp(-w) ;
		    n[k] += w ;
		    ss[k] += slx1->signal[k] * w ;
		    break ;
		  case 2:
		    w = i - iSlx ;
		    w = NN * NN - w * w ;
		    n[k] += w ;
		    ss[k] += slx1->signal[k] * w ;
		    break ;
		  case 3:
		    n[k]++ ;
		    ss[k] += slx1->signal[k] ;
		    break ;
		  case 4:
		    array(aaa[k],n[k],double) = slx1->signal[k] ;
		    n[k]++ ;
		    break ;
		  }
	      }
	  }
    for (k = 0,  pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p &&  k < SOLEXAMAX ; pnx++, k++)
	{
	  if (
	      (!look->ratio && look->isSolexaClosed [k]) ||
	      (look->ratio && ! (k & 1) && look->isSolexaClosed [k]) ||
	      (look->ratio &&   ( k& 1) && look->isSolexaClosed [k-1]) ||
	      (1 && look->isSolexaSmoothed[k])
	      )
	    continue ;

	  if (pnsAll[k].isRna == 92)
	    slx->gaussSignal[k] = slx->signal[k] ;
	  else if (n[k] > 0)
	    {
	      if (look->smoothing == 4)
		{
		  if (NN > 1)
		    {
		      arraySort (aaa[k], doubleOrder) ;
		      slx->gaussSignal[k] = array (aaa[k], n[k]/2, double) ;
		    }
		  else
		    slx->gaussSignal[k] = slx->signal[k] ;
		}
	      else 
		slx->gaussSignal[k] = ss[k]/n[k] ;
	    }
	  else
	    slx->gaussSignal[k] = 0 ;
	}
    }
 
  if (look->ratio && look->smoothing)
    for (k = 0,  pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p &&  k < SOLEXAMAX ; pnx+=2, k+=2)
      look->isSolexaSmoothed[k+1] = look->isSolexaSmoothed[k] = !look->isSolexaClosed [k] ;
  else  if (!look->ratio && look->smoothing)
    for (k = 0,  pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p &&  k < SOLEXAMAX ; pnx++, k++)
      look->isSolexaSmoothed[k] = !look->isSolexaClosed [k] ;

  /*
    for (k = 0, pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p &&  k < SOLEXAMAX ; pnx++, k++)
    look->isSolexaSmoothed[k] = FALSE ;
  */
 done:
  ac_free (h) ;

  return TRUE ;
} /* htileSolexaSmoothing */

/************************************************************/

static BOOL htileGeneSmoothingDraw (Htile look, float offset)
{
  int ns, nGene, nGene1, iSeg ;
  float x1, x2, xmax, y ;
  int b1, b2, state ;
  double signal[NSMAX] ;
  int n[NSMAX] ;
  SEG *seg, *gSeg, *gSeg1 ;
  float oldWidth = graphLinewidth (.6) ;

  for (nGene = 0 ; look->map->genes  && nGene < arrayMax(look->map->genes) ; nGene++)
    {
      gSeg = arrp (look->map->genes, nGene, SEG) ; 
      switch (look->smoothing)
	{
	case 5: /* up strand */
	  if (gSeg->x1 > gSeg->x2) gSeg->key = 0 ;
	  break ;
	case 6: /* down strand */
	  if (gSeg->x1 < gSeg->x2) gSeg->key = 0 ;
	  break ;
	}
    }
  for (nGene = 0 ; look->map->genes  && nGene < arrayMax(look->map->genes) ; nGene++)
    {
      gSeg = arrp (look->map->genes, nGene, SEG) ;
      if (!gSeg->key) continue ;
      for (nGene1 = nGene + 1, gSeg1 = gSeg+1 ; nGene1 < arrayMax(look->map->genes) ; gSeg1++, nGene1++)
	if (gSeg1->key && gSeg1->a1 < gSeg->a2)
	  {
	    if (gSeg1->a2 > gSeg->a2)
	      gSeg->a2 = gSeg1->a2 ;
	    gSeg1->key = 0 ; 
	  }
    }
  xmax = look->map->graphWidth ;
  /* detect first relevant gene */
  b1 = 0 ; state = 0 ; /* default start averaging at far left */
  b2 = TGRAPH2MAP (look->map, 1) ; /* left end of window */
  for (gSeg = 0, nGene = 0 ; look->map->genes && nGene < arrayMax(look->map->genes) ; nGene++)
    {
      gSeg = arrp (look->map->genes, nGene, SEG) ;
      if (!gSeg->key) continue ;
      if (gSeg->a2 > b2)
	{
	  if (gSeg->a1 < b2)
	    {
	      state = 1 ;
	      b1 = gSeg->a1 ;
	    }
	  else 
	    {
	      state = 0 ;
	      b1 = nGene ? (gSeg-1)->a2 : b2 ;
	    }
	  break ;
	}
      gSeg = 0 ;
    }

  /* initialise */
  for (ns = 0 ; ns < NSMAX ; ns++)
    { n[ns] = 0 ; signal[ns] = 0 ; }

  for (iSeg = 0, seg = arrp (look->map->segs, 0, SEG); state >=0 && iSeg < arrayMax(look->map->segs) ; iSeg++, seg++)
    {
      if (seg->a1 < b1)
	continue ;

      switch (state)
	{
	case 0: /* i am outside */
	  if (gSeg && seg->a2 > gSeg->a1) /* enter new gene section */
	    state = 3 ; /* report previous section */
	  break ;
	case 1: /* i am inside */
	  if (gSeg && seg->a1 > gSeg->a2) /* exit the gene */
	    state = 2 ; /* report previous section */
	  break ;
	}
      switch (state)
	{
	case 0:
	case 1: /* accumulate */
	  for (ns = 0 ; ns < NSMAX ; ns++)
	    if (seg->signal[ns] != -1000)
	      { n[ns]++ ; signal[ns] += seg->signal[ns] ;}
	  break ;
	case 2:
	case 3: /* draw */
	  if (state == 2) /* exiting */
	    b2 = gSeg->a2 ;
	  else
	    b2 = gSeg->a1 ;
	  x1 = TMAP2GRAPH (look->map, b1) ;
	  if (x1 < 1) x1 = 1 ;
	  x2 = TMAP2GRAPH (look->map, b2) ;
	  if (x2 > xmax)
	    { x2 = xmax ; state = -1 ; }
	  for (ns = 0 ; ns < NSMAX ; ns += look->ratio ? 2 : 1)
	    if (!look->isClosed[ns] && n[ns])
	      {
		y = offset - look->zoom * signal[ns]/n[ns] ;
		graphColor (pnsAll[ns].col) ;
		graphLine (x1, y, x2, y) ;
	      }
	  state -= 2 ; /* reset */
	  for (ns = 0 ; ns < NSMAX ; ns++)
	    {
	      n[ns] = 0 ; signal[ns] = 0 ;
	      if (seg->signal[ns] != -1000)
		{
		  signal[ns] = seg->signal[ns] ;
		  n[ns]++ ;
		}
	    }
	  if (state == 1 && gSeg->a2 > seg->a1) /* entering a gene */
	    b1 = gSeg->a1 ;  
	  else /* exiting a gene */
	    {
	      b1 = gSeg->a2 ;
	      for (nGene++ ; nGene < arrayMax (look->map->genes) ; nGene++)
		{
		  gSeg = arrp (look->map->genes, nGene, SEG) ;
		  if (gSeg->key) break ;
		  gSeg = 0 ;
		}
	    }
	  if (gSeg && gSeg->a1 < seg->a1 && gSeg->a2 > seg->a1)
	    state = 1 ;
	  else
	    state = 0 ;
	  break ;
	}
    }
  graphLinewidth (oldWidth) ;
  return TRUE ;
} /* htileGeneSmoothingDraw */

/************************************************************/

static void htileExportSmoothedValues (Htile look, ACEOUT ao)
{
  int nn = 0, ii, x = 0, dx = 0, ns ;
  SEG *seg ;
  float s, sss[200] ; 
  int nns[200] ;
  PNS *pns ;
  PGN *pgn ;

  look->ratio = TRUE ;
  dx = look->gaussWidth/5 ;
  for (seg = arrp(look->map->segs, 0, SEG), nn = ii = 0   ; ii < arrayMax (look->map->segs) ; ii++, seg++)
    {
      if (ii == 0 || seg->a1 > x + dx)
	{ 
	  if (ii == 0)
	    {
	      aceOutf (ao, "# Chrom\tPosition") ;
	      for (ns = 0, pns = pnsAll ; pns->p  ; ns += look->ratio ? 2 : 1, pns += look->ratio ? 2 : 1)
		aceOutf (ao, "\t%s", look->ratio ? pns->nam2 : pns->nam) ;
	      aceOutf (ao, "\n") ;
	    }
	  if (nn) /* export */
	    {
	      aceOutf (ao, "%s\t%d", name (seg->map), x) ;
	      for (ns = 0, pns = pnsAll ; pns->p  ; ns += look->ratio ? 2 : 1, pns += look->ratio ? 2 : 1)
		aceOutf (ao, "\t%.3f", nns[ns] > 0 ? 10.0 + sss[ns]/nns[ns] : -1000.0) ; 
	      aceOutf (ao, "\n") ;
	    }
	  nn = 0 ;
	  x = seg->a1 - (seg->a1 % dx) ;
	  memset (nns, 0, sizeof (nns)) ;
	  memset (sss, 0, sizeof (sss)) ;
	}
      if (look->rejectAmbiguous < seg->ambiguous)
	continue ;
      switch (look->tmFilter)
	{
	case 1: 
	  if (seg->tm > 72.2)
	    continue ;
	  break ;
	case 2: 
	  if (seg->tm < 75)
	    continue ;
	  break ;
	}
      switch (look->filter99)
	{
	case 2:
	  if (seg->flag & MOTIF_FILTER)
	    continue ;
	  break ;
	}
      nn++ ;
      for (ns = 0, pns = pnsAll ; pns->p && ns < NSMAX  ; ns += look->ratio ? 2 : 1, pns += look->ratio ? 2 : 1)
	{
	  if (!look->ratio && *pns->nam == '.') continue ;
	  pgn = pnsPgn (look->ratio, pns) ;
	  if (
	      (pgn && look->filter99 == 3 && !look->ratio && seg->signal[ns] + 10 > pgn->s99) ||
	      (pgn && look->filter99 == 3 && look->ratio &&
	       seg->signal[ns] - seg->signal[ns+1] > pgn->s99
	       )
	      )
	    continue ;
	  if (0 && look->ratio && ! pgn->nam) continue ;
	  if (seg->signal[ns]==-1000)
	    continue ;
	  if (look->ratio && *(pns+1)->p != '.' && seg->signal[ns+1]==-1000)
	    continue ; 
	  if (look->smoothing && seg->gaussSignal[ns] == -1000)
	    continue ;
	  if (look->smoothing) 
	    {
	      if (look->ratio && *(pns+1)->p != '.' )
		{
		 if (!strcmp (pns->nam2, "G+C") || !strcmp (pns->nam2, "A+T")) 
		   s = seg->gaussSignal[ns] + seg->gaussSignal[ns+1] ;
		 else
		   s = seg->gaussSignal[ns] - seg->gaussSignal[ns+1] ;
		}
	      else
		s = seg->gaussSignal[ns] ;
	    }
	  else
	    {
	      if (look->ratio && *(pns+1)->p != '.' )
		{ 
		  if (!strcmp (pns->nam2, "G+C") || !strcmp (pns->nam2, "A+T")) 
		    s = seg->signal[ns] + seg->signal[ns+1] ;
		  else
		    s = seg->signal[ns] - seg->signal[ns+1] ;
		}
	      else
		s = seg->signal[ns] ;
	    }
	  if (look->index && pgn)
	    s -= pgn->s20 - 10.0 ;
	  if (s > -900)
	    {
	      nns[ns]++ ;
	      sss[ns] += s ;
	    }
	}
    }  
  if (nn) /* export */
    {
      aceOutf (ao, "%s\t%d", name (seg->map), x) ;
      for (ns = 0 ; ns < NSMAX  ; ns += look->ratio ? 2 : 1)
	if (nns[ns] > 0)
	  aceOutf (ao, "\t%.3f", nns[ns] > 0 ? sss[ns]/nns[ns] : -1000.0) ; 
      aceOutf (ao, "\n") ;
    }
} /* htileExportSmoothedValues */

/************************************************************/
Array uArrayHandleGet(KEY key, int size, char *format, AC_HANDLE handle) ;

static BOOL htileConvert (Htile look, BOOL force, ACEOUT ao)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ Cosmid = 0 ;
  AC_TABLE tbl ;
  int ii, a1 = 0, a2 = 0 ;
  BOOL ok = FALSE ;
  KEYSET ks1 = 0 ;
  KEY cacheKey ;
  char hSegFormat[256] ;
  char suffix[64] ;

  if (class(look->key) == _VLocus)
    {
      look->locus = look->key ; look->type = 1 ;
      ks1 = queryKey (look->key, "{>Sequence} $| {>Sequence parts ; >Parts}; IntMap && Solexa") ;
    }
  else if (class(look->key) == _VmRNA)
    {
      look->locus = look->key ; look->type = 1 ;
      ks1 = queryKey (look->key, "CLASS mRNA") ;
    }
  else if (class(look->key) == _VGene)
    {
      look->locus = look->key ; look->type = 1 ;
      ks1 = queryKey (look->key, ">Genomic_sequence;{>Sequence} $| {>Sequence parts ; >Parts}; IntMap && Solexa") ;
    }
  else if (class(look->key) == _VSequence)
    { 
      ks1 = queryKey (look->key, "{CLASS Sequence} $| {>Parts} ; IntMap && (Solexa || COUNT{>intMap  Wiggle} > 0)  ") ;
      look->locus = look->cosmid = look->key  ; look->type = 1 ; 
    }
  
  if (!ks1 || ! keySetMax (ks1))
    goto abort ;
  
  Cosmid = htile_key_obj (look, look->key, h) ;
  if (!Cosmid)
    goto abort ;

  if (class(look->key) == _VmRNA)
    {
      if ((tbl = ac_tag_table (Cosmid, "DNA", h) ))
	{
	  look->intMap = look->key ;
	  a1 = 0 ;
	  a2 = ac_table_int (tbl, 0, 1, 0) ; 
	}
    }
  else if ((tbl = ac_tag_table (Cosmid, "IntMap", h) ))
    {
      look->intMap = ac_table_key (tbl, 0, 0, 0) ; 
      a1 = ac_table_int (tbl, 0, 1, 0) ; 
      a2 = ac_table_int (tbl, 0, 2, 0) ;
      look->from += a1 ;
    }

  if (a1 > a2) { int a0 = a1 ; a1 = a2 ; a2 = a0 ; }
  look->map->a1 = a1 ;
  look->map->a2 = a2 ; /* backwards on purpose */
  ac_free (Cosmid) ;
  
  /* segFormat */
  {
    int i, j ;
    for (i=0 ; i < segFormatNint ; i++)
      hSegFormat[i] = 'k' ;
    for (j=0 ; j < segFormatNfloat ; j++)
      hSegFormat[i+j] = 'f' ;
    hSegFormat[i+j] = 0 ;
    if (i+j > 255 || i * sizeof(int) + j * sizeof (float) != sizeof (SEG))
      messcrash ("bad declaration of htile SEG structure size") ;
  }
  
  switch (look->romainSmoothing)
    {
    case 0: strcpy (suffix, "hTile") ; break ;
    case 1: strcpy (suffix, "hTileSmooth") ; break ;
    case 2: strcpy (suffix, "hTileEric") ; break ;
    }
  suffix[4]= VERSION ;/* ATTENTION: version */
  
  htileInit () ;
  arrayDestroy (look->map->segs) ;


  if (class(look->key) == _VmRNA)
    {
      look->map->segs = arrayHandleCreate (10000, SEG, look->h) ;
      htileMrnaConvert (look) ;
    }
  else if (! force && (cacheKey = gMapCacheKey(look->locus, suffix)) &&
	   (look->map->segs = uArrayHandleGet(cacheKey, sizeof(SEG), hSegFormat, look->h))) ;
  else
    {
      KEYSET knownProbes = keySetCreate () ; 
      look->map->segs = arrayHandleCreate (10000, SEG, look->h) ;
      for (ii = 0 ; ii < keySetMax (ks1) ; ii++)
	htileCosmidConvert (look, keySet (ks1, ii), knownProbes) ;
      arraySort (look->map->segs,  htileSegOrder) ;
      arrayCompress (look->map->segs) ;
      keySetDestroy (knownProbes) ;

      gMapCacheStore(look->locus, suffix, look->map->segs, hSegFormat) ;
    }
    
  ok = TRUE ;
  look->map->max = look->map->a2 - look->map->a1 ;
  keySetDestroy (ks1) ;

  for (ii = 0 ; ii < NSMAX ; ii++)
    look->isSmoothed[ii] = FALSE ;
  for (ii = 0 ; ii < SOLEXAMAX ; ii++)
    look->isSolexaSmoothed[ii] = FALSE ;
  htileRemoveMedian (look) ;
  switch (look->smoothing)
    {
    case 1: /* gauss smoohing */
    case 2: /* parabolic smoohing */
    case 3: /* square smoohing */
    case 4: /* median smoohing */
      htileSmoothing (look) ;
      if (look->map->solexa)
	htileSolexaSmoothing (look) ;
      break ;
    }

  /* export a value every look->gaussWidth/5 */
  if (ao)
    htileExportSmoothedValues (look, ao) ;

 abort:
  ac_free (h) ;
  return ok ;
} /* htileConvert */

/************************************************************/

static void htileLnWidthPlus (void)
{
  Htile look = currentHtile("htileNoAmbiguousProbes") ;
  oldLnWidth = look->map->lnWidth = oldLnWidth + .1 ;
  look->map->mapDraw () ;
}

/************************************************************/

static void htileLnWidthMinus (void)
{
  Htile look = currentHtile("htileNoAmbiguousProbes") ;
  if (oldLnWidth > .1)
    { oldLnWidth = look->map->lnWidth = oldLnWidth - .1 ; }
  look->map->mapDraw () ;
}

/************************************************************/

static void htileHideHeader (void)
{
  Htile look = currentHtile("htileNoAmbiguousProbes") ;
  oldHideHeader = look->hideHeader = ! look->hideHeader ;
  look->map->mapDraw () ;
}

/************************************************************/

static void htileHideMusic (void)
{
  Htile look = currentHtile("htileNoAmbiguousProbes") ;
  oldHideMusic = look->hideMusic = ! look->hideMusic ;
  look->map->mapDraw () ;
}

/************************************************************/

static void htileHideScale (void)
{
  Htile look = currentHtile("htileNoAmbiguousProbes") ;
  oldHideScale = look->hideScale = ! look->hideScale ;
  look->map->mapDraw () ;
}

/************************************************************/

static void htileHideRZones (void)
{
  Htile look = currentHtile("htileNoAmbiguousProbes") ;
  look0->showRZones = look->showRZones = ! look->showRZones ;
  look->map->mapDraw () ; 
}

/************************************************************/

static void htileDisplayPreserve (void)
{
  Htile look = currentHtile("htileDisplayPreserve") ;

  displayPreserve () ;
  look->preserve = look0->preserve = TRUE ;
}

/************************************************************/

static void htileMenuSelect (Htile look)
{
  int ii = 0 ;
  MENUOPT *mm ;
  
  arrayDestroy (look->menu) ;
  look->menu = arrayHandleCreate (NSMAX, MENUOPT, look->h) ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = graphDestroy ; mm->text = "Quit" ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = help ; mm->text = "Help";

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = graphPrint ; mm->text = "Print" ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = htileDisplayPreserve ; mm->text = "Preserve" ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = htileHideHeader ; mm->text = look->hideHeader ? "Show header" : "Hide header" ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = htileHideScale ; mm->text = look->hideScale ? "Show genome scale" : "Hide genome scale" ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = htileHideMusic ; mm->text = look->hideMusic ? "Show Vertical scale" : "Hide vertical scale" ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = htileLnWidthPlus ; mm->text = "+Line width" ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = htileLnWidthMinus ; mm->text = "-Line width" ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = htileAddWall ; mm->text ="Add Wall" ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = htileHideRZones ; mm->text ="Show/Hide R-zones" ;

  /*
    mm = arrayp (look->menu, ii++, MENUOPT) ;
    mm->f = htileAddWallComment ; mm->text ="Comment new Walls" ;
  */
  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = htileExport ; mm->text ="Export" ;
  
  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = htileForceRecalculate ; mm->text ="Recalculate" ;
  
#ifdef BOU_DEF  
  if (!strcmp (getLogin (TRUE), "mieg"))
    {
      mm = arrayp (look->menu, ii++, MENUOPT) ;
      mm->f = htileCreateRZone ; mm->text ="Create R-zone" ;
    }
#endif


  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = 0 ; mm->text = 0;


} /* htileMenuSelect */

/************************************************************/
typedef struct HTILEBPstruct { KEY mrna ; int bumpy ;} HTILEBP ;
/************************************************************/

static void htileDrawMusicalScale (Htile look, float offset)
{
  float oldw, s, x, y, zoom ;
  int j, j0, ss[] = { 1, 2, 5 } ;
#ifdef BOU_DEF
  int *ep, erics[] = {-20, 8, 36, 64, 92, 120, 0} ;
#endif
  int ds = look->index ? 20 : 10 ; /* hauteur effective du dessin */

  x = TMAP2GRAPH (look->map, look->map->min) ;
  s = 0 ;
  zoom = look->zoom ;
#ifdef BOU_DEF
  if (look->index) zoom = look->zoom * .005 ; /* so that 100 gives 10 char size when ratio zoom == 20 */
#endif
  for (x = 1, j0 = -1 ; s == 0 && x < 1e30 ; x*= 10)
    for (j = 0 ;  s == 0 && j < 3 ; j++)
      {
	s = x * ss[j]/1000.0 ;
	if (s * zoom < ds) s = 0 ;
	else j0 = j ;
      }
  if (s == 0) return ;
  graphColor (GRAY) ;
#ifdef BOU_DEF
  if (look->index)
    {
      for (ep = erics ; *ep ; ep++)
	{
	  y = offset - zoom * (*ep) ; graphLine (3, y, look->map->graphWidth -1, y) ;
	}
    }
  else
#endif
    {
      if (j0 == 1)
	{
	  for (x = -s ; x <=  s ; x += s/4.0)
	    {
	      y = offset + zoom * x ; 
	      if (! look->hideMusic && y > 0 && y < look->map->graphHeight)
		graphLine (3, y, look->map->graphWidth -1, y) ;
	    }
	}
      else
	{
	  for (x = -s ; x <= s ; x += s/5.0)
	    {
	      y = offset + zoom * x ; 
	      if (! look->hideMusic && y > 0 && y < look->map->graphHeight) 
		graphLine (3, y, look->map->graphWidth -1, y) ;
	    }
	}
    }
  
  graphColor (BLACK) ;
  y = offset - zoom * s ; 
  if (y > 0) 
    {
      int box = graphBoxStart () ;
      graphText (messprintf (" %g ", s), 1, y) ;
      graphBoxEnd () ;
      graphBoxDraw (box, BLACK, WHITE) ;
    }
  if (! look->hideMusic)
    graphLine (3, y, look->map->graphWidth -1, y) ;
  if (!look->index)
    {
      y = offset + zoom * s ; 
      if (! look->hideMusic && y < look->map->graphHeight)
	{
	  graphLine (3, y, look->map->graphWidth -1, y) ;
	  graphText (messprintf ("%g", - s), 1, y) ;
	}
    }
  oldw = graphLinewidth (.3) ;
  y = offset ; graphLine (1, y, look->map->graphWidth -1, y) ;
  graphLinewidth (oldw) ;
  return ;
} /* htileDrawMusicalScale */

/************************************************************/

static void htileDrawProbes (Htile look, float offset, float probeOffset, ACEOUT ao)
{
  int ns, iSeg, npp1, npp2, npf, npfm ;
  SEG *seg ;
  PNS *pns ;
  PGN *pgn ;
  float s, x, y, yy[NSMAX], oldx[NSMAX], zoom, oldWidth = 0, oldWidth1 ;
  int i, np1[NSMAX], iao = 0 ;
  int iSegMax = arrayMax (look->map->segs) ;
  float mxMx = -1000, mxMy = -1000, mxDelta = 2 ;
  BOOL mxUp = TRUE ;
  float mxMy1 = probeOffset + 1, mxMy2 = look->map->graphHeight - .5 ;

  for (ns = 0 ; ns < NSMAX ; ns++) 
    { np1[ns] = 0 ; oldx[ns] = -1 ; yy[ns] = offset ; }
  npf = npfm = npp1 = npp2 = 0 ;

  if (!ao)
    oldWidth = graphLinewidth (.6) ;
  if (iSegMax)
    for (iSeg = 0, seg = arrp (look->map->segs, 0, SEG); iSeg < iSegMax ; iSeg++, seg++)
      {
	x = TMAP2GRAPH (look->map, seg->a1) ;
	if (x < 0 || x > look->map-> graphWidth)
	  continue ;
	if (seg->chip == 1) npp1++ ; 
	else if (seg->chip == 2) npp2++ ;
	if (look->rejectAmbiguous < seg->ambiguous)
	  { npf++ ; continue ; }
	switch (look->tmFilter)
	  {
	  case 1: 
	    if (seg->tm > 72.2)
	      { npf++ ; continue ; }
	    break ;
	  case 2: 
	    if (seg->tm < 75)
	      { npf++ ; continue ; }
	    break ;
	  }
	switch (look->filter99)
	  {
	  case 2:
	    if (seg->flag & MOTIF_FILTER)
	      { npfm++ ; continue ; }
	    break ;
	  }
	if (
	    (look->exonicFilter == 1 && !(seg->flag & EXONIC_FILTER)) ||
	    (look->exonicFilter == 2 && !(seg->flag & NEW_EXONIC_FILTER)) ||
	    (look->exonicFilter == 3 && !(seg->flag & (EXONIC_FILTER | NEW_EXONIC_FILTER))) ||
	    (look->exonicFilter == 4 && !(seg->flag & INTRONIC_FILTER)) ||
	    (look->exonicFilter == 5 && !(seg->flag & INTERGENIC_FILTER))
	    )
	  continue ;
	
	if (!ao)
	  {
	    graphColor (BLACK) ;
	    if (! look->smoothing)
	      {
		array (look->map->boxIndex, graphBoxStart (), int) = iSeg ;
		graphText ("*", x-.5, probeOffset) ;
		graphBoxEnd () ;
	      }
	  }
	if (ao)
	  iao = 0 ;
	for (ns = 0, pns = pnsAll ; pns->p && ns < NSMAX  ; ns += look->ratio ? 2 : 1, pns += look->ratio ? 2 : 1)
	  {
	    if (look->isClosed[ns]) continue ;
	    if (!look->ratio && *pns->nam == '.') continue ;
	    pgn = pnsPgn (look->ratio, pns) ;
	    if (look->filter99 == 3)
	      {
		if (pgn && ! look->ratio && 
		    seg->signal[ns] + 10 > pgn->s99
		    )
		  continue ;
		if (pgn && look->ratio && 
		    seg->signal[ns+1] - seg->signal[ns] > pgn->s99
		    )
		  continue ;
	      }
	    zoom = look->zoom ;
	    if (!strcmp (pns->nam, "TM")) zoom = 2 ;
	    else if (!strcmp (pns->nam2, "G+C")) zoom = 1 ;
	    else if (!strcmp (pns->nam, "G")) zoom = 1 ;
	    else if (!strcmp (pns->nam, "C")) zoom = 1 ;
	    else if (!strcmp (pns->nam2, "A+T")) zoom = 1 ;
	    else if (!strcmp (pns->nam, "T")) zoom = 1 ;
	    else if (!strcmp (pns->nam, "A")) zoom = 1 ;
	    else if (!strcmp (pns->nam2, "G-C")) zoom = 1 ;
	    else if (!strcmp (pns->nam, "G.")) zoom = 1 ;
	    else if (!strcmp (pns->nam, "C.")) zoom = 1 ;
	    else if (!strcmp (pns->nam2, "A-T")) zoom = 1 ;
	    else if (!strcmp (pns->nam, "T.")) zoom = 1 ;
	    else if (!strcmp (pns->nam, "A.")) zoom = 1 ;
	    else  if (pgn && look->index)  /* so that if zoom == 20, and s == max zoom*s = 10 chars */
	      {
		float zmax ;
		if (pns->flag & (PGG_NUC | PGG_TOTAL))
		  zmax = pgn->s99 ;
		else 
		  zmax = pgn->s80 ;
		zoom *= 1.0/(zmax - 10.0) ;
	      }
	    if (seg->signal[ns]==-1000)
	      continue ;
	    if (look->ratio && *(pns+1)->p != '.' && seg->signal[ns+1]==-1000)
	      continue ; 
	    if (look->smoothing && seg->gaussSignal[ns] == -1000)
	      continue ;
	    if (look->smoothing) 
	      {
		if (look->ratio && *(pns+1)->p != '.' )
		  {
		    if (!strcmp (pns->nam2, "G+C") || !strcmp (pns->nam2, "A+T")) 
		      s = seg->gaussSignal[ns] + seg->gaussSignal[ns+1] ;
		    else
		      s = seg->gaussSignal[ns] - seg->gaussSignal[ns+1] ;
		  }
		else
		  s = seg->gaussSignal[ns] ;
	      }
	    else
	      {
		if (look->ratio && *(pns+1)->p != '.' )
		  { 
		    if (!strcmp (pns->nam2, "G+C") || !strcmp (pns->nam2, "A+T")) 
		      s = seg->signal[ns] + seg->signal[ns+1] ;
		    else
		      s = seg->signal[ns] - seg->signal[ns+1] ;
		  }
		else
		  s = seg->signal[ns] ;
	      }
	    if (0 && look->index && pgn) /* mieg 2008_04_03 :  0 && :: do NOT modify the manual base line pns->av */
	      s -= pgn->s20 - 10 ; /* 10, since the pgn histos are computed after shifting pns->av to 10 */
	    y = offset - zoom * s ;
	    if (ao)
	      {
		if (!iao++)
		  aceOutf (ao, "\n%d", seg->a1) ;
		aceOutf (ao, "\t%g", -y) ;
	      }
	    else
	      {  
		if (y < 2) y = 2 ;
		
		graphColor (pnsAll[ns].col) ;
		if (1) 
		  {
		    if (look->showDot == 1)
		      graphCircle (x, y, .2) ;
		    else if (look->showDot == 2)
		      graphLine (x, offset, x, y) ;
		    else if (look->showDot == 3)
		      graphLine (x, offset + yy[ns], x, y + yy[ns]) ;
		    else if (oldx[ns] > 0 &&  np1[ns])
		      graphLine (oldx[ns], yy[ns], x, y) ;
		    np1[ns]++ ;
		    oldx[ns] = x ;
		    yy[ns] = (look->showDot == 3 ? y + yy[ns] : y) ;
		  }
		graphColor (BLACK) ;
		
		/* show lines at extremal points */
		if (look->showExtrema && (!strcmp (pns->nam,"S=31ES") || (ns && look->mxShowNs == ns)))
		  {
		    if ((mxUp && mxMy < y) || (!mxUp && mxMy > y)) { mxMx = x ; mxMy = y ; }
		    else if ((mxUp && mxMy > y + mxDelta) ||
			     (!mxUp && mxMy < y - mxDelta)
			     ) /* report previous */
		      {
			mxUp = !mxUp ;
			if (mxMx > 0) 
			  {
			    graphColor (LIGHTMAGENTA) ;
			    oldWidth1 = graphLinewidth (.2) ;
			    graphLine (mxMx, mxMy1, mxMx, mxMy2) ;
			    graphLinewidth (oldWidth1) ;
			  }
			mxMx = x ; mxMy = y ;
		      }
		  } 	    /* end of drawing extremal points */
	      }
	  }
      }

  if (!ao)
    {
      graphLinewidth (oldWidth) ;
      graphColor (BLACK) ;
    }

  for (i = ns = 0 ; ns < NSMAX ; ns++)
    if (np1[ns] > i) i = np1[ns] ;
  if (0)
    graphText (messprintf ("chip1 %d probes, chip2 %d, filtered %d, motif-filtered %d, shown %d points", npp1, npp2, npf, npfm, i)
	     , 2, 2.6) ;
  if (ao)
    aceOutf (ao, "\n") ;
} /* htileDrawProbes */

/************************************************************/
/************************************************************/

static int htileWiggleConvert (Htile look, PNX *pnx, int ns, Array wpArray)
{
  int ii, nn2, a1, nn ;
  Array slxs2 = 0 ;
  WIGGLEPOINT *wp ;
  SLX *slx ;
  int solexaStep = look->solexaStep  ; 
  int step = solexaStep ; /* resolution of this experiment */

  /* in old method we parse at actual positions, so the tads have random positions
   * in new method we systematically create a tag every solexastep
   */
  if (!look->map->solexa)
    look->map->solexa = arrayHandleCreate (100000, SLX, look->h) ;
  slxs2 = look->map->solexa ;
   
  /* parse the data in a new array */
  nn =  (look->map->a2 - look->map->a1 )/step + 1 ;
  slx = arrayp (slxs2, nn - 1, SLX) ; /* expand */
  nn = nn2 = 0 ;
  for (ii = 0, wp = arrp (wpArray, 0, WIGGLEPOINT) ; ii < arrayMax (wpArray); wp++, ii++)
    {
      /* wp->x = position on sequence */
      /* wp->y = number of supporting sequences */
      if (wp->x >= look->map->a1 && wp->x < look->map->a2) 
	{
  	  nn++ ;
	  a1 = (wp->x - look->a1 + step/2 + 1) / step ;
	  if (a1 > nn2) nn2 = a1 ;
	  if (a1 >= 1)
	    {
	      slx = arrayp (slxs2, a1 - 1, SLX) ;
	      slx->signal[ns] += wp->y ;
	    }
	}
    }
  arrayMax (slxs2) = nn2 ;
  nn2 = 0 ;
  /* add in the missing coordinates */

  if (nn)
    for (ii = 0, slx = arrp (slxs2, ii, SLX) ; ii < arrayMax (slxs2) ; slx++, ii++)
      slx->a1 = ii * step ;

  if (look->map->solexa)
    {
      nn = arrayMax (look->map->solexa) ;
      slx = arrayp (look->map->solexa, 0, SLX) ;
      slx->signal[ns] = nn+1 ;
    }
  else
    nn = 1 ;
  pnx->signal = (Array)assVoid(nn+1);

  return nn ;
} /* htileWiggleConvert */

/************************************************************/
/* the endRatio
 * ATTENTION  synchronize with w1/wigglelib.c:sxWiggleEndRatio1
 */

static void htileSolexaEndRatios (Htile look, PNX *pnx0, int ns, int NF)
{
  PNX *pnx ;
  int NA = -999 ;
  unsigned int flag1 = 0, flag2 = 0, flag3 = 0, flag4 = 0 ;
  int i, ns1, ns2, ns3, ns4 ;
  
  ns1 = ns2 = ns3 = ns4 = NA ;
  if ((pnx0->flag &  PGG_endRatioLF) ==  PGG_endRatioLF)
    {  flag1 = PGG_ELF ; flag2 = PGG_ERF ;  flag3 = PGG_ELR ; flag4 = PGG_ERR ;}
  if ((pnx0->flag &  PGG_endRatioRF) ==  PGG_endRatioRF)
    {  flag1 = PGG_ERF ; flag2 = PGG_ELF ;  flag3 = PGG_ERR ; flag4 = PGG_ELR ; }
  if ((pnx0->flag &  PGG_endRatioLR) ==  PGG_endRatioLR)
    {  flag1 = PGG_ELR ; flag2 = PGG_ERR ; flag3 = PGG_ELF ; flag4 = PGG_ERF ; }
  if ((pnx0->flag &  PGG_endRatioRR) ==  PGG_endRatioRR)
    {  flag1 = PGG_ERR ; flag2 = PGG_ELR ; flag3 = PGG_ERF ; flag4 = PGG_ELF ; }

  for (i = 0, pnx = pnx0 ; i < NF && ns -i >= 0 ; pnx--, i++)
    if ((pnx->flag & flag1) == flag1) { pnx0->signal = pnx->signal ; ns1 = ns - i ; break ; }
  for (i = 0, pnx = pnx0 ; i < NF && ns -i >= 0 ; pnx--, i++)
    if ((pnx->flag & flag2) == flag2) { ns2 = ns - i ; break ; }
  for (i = 0, pnx = pnx0 ; i < NF && ns -i >= 0 ; pnx--, i++)
    if ((pnx->flag & flag3) == flag3) { ns3 = ns - i ; break ; }
  for (i = 0, pnx = pnx0 ; i < NF && ns -i >= 0 ; pnx--, i++)
    if ((pnx->flag & flag4) == flag4) { ns4 = ns - i ; break ; }

  if (ns1 != NA && ns2 != NA && look->map->solexa)
    {
      unsigned int iMax =  arrayMax (look->map->solexa) ;
      SLX *slx, *slx1 ;

      /* ATTENTION synchronize with sxWiggleEndRatioOne */
      for (int ii = 6 ; ii < iMax - 6 ; ii++)
	{
	  float x = 0, y = 0, z=0, t=0, u=0 ;
	  int damper = 10 ;
	  float seuil = .80 ;
	  float zoom = 40 ;
	  int nn = 1 ; /* was 5, so 2n+1 = 110 bases */

	  slx = arrayp (look->map->solexa, ii, SLX) ;

	  for (int j = -nn ; j <= nn ; j++)
	    {
	      slx1 = arrayp (look->map->solexa, ii+j, SLX) ;
	      x += slx1->signal[ns1] ;
	      y += slx1->signal[ns2] ;
	      if (ns3 != NA) z += slx1->signal[ns3] ; 
	      if (ns4 != NA) t += slx1->signal[ns4] ; 
	    }
	  /* do not use the other strand values (s,t) this creates ends echos rather than damping them */
	  x /= 2 * nn + 1 ; y /= 2 * nn + 1 ;
	  u =  (x + damper) / (x + y + 2 * damper) - seuil ;
	  if (u < 0) u = 0 ;
	  slx->signal[ns] = zoom * u * x ;
	}
    }
  if (ns1 != NA)   look->isClosed[ns1] = look->isSolexaClosed[ns1] = 2 ;  
  if (ns2 != NA)   look->isClosed[ns2] = look->isSolexaClosed[ns2] = 2 ;  
  if (ns3 != NA)   look->isClosed[ns3] = look->isSolexaClosed[ns3] = 2 ;  
  if (ns4 != NA)   look->isClosed[ns4] = look->isSolexaClosed[ns4] = 2 ;  
  look->isClosed[ns] = look->isSolexaClosed[ns] = 0 ;  
  /*   pnx0->signal = w1 ; */

} /* htileSolexaEndRatios */

/************************************************************/

static void htileSolexaConvert (Htile look, BOOL force, ACEOUT ao)
{
  /* 2010_01_10 store the wiggle info in IntMap->wiggle */
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ Map = 0 ;
  AC_TABLE tbl = 0 ;
  int ns = 0 ;
  PNX *pnx ;
  
  
  if (force)
    for (ns = 0,  pnx = arrp(look->solexaAll, 0, PNX) ; pnx->nam ; ns++, pnx++)
      pnx->signal = 0 ;
  
  arrayDestroy (look->map->solexa) ;
  if (class(look->key) == _VmRNA)
    {
      Array aa = 0 ;
      WIGGLEPOINT *wp ;
      int ir, jr, x ;
      float y ;

      aa = arrayHandleCreate (1000, WIGGLEPOINT, h) ;
      Map = htile_key_obj (look, look->key, h) ;
      tbl = Map ? ac_tag_table (Map, "Wiggle", h) : 0 ;
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  char *manip = hprintf (h, "%s", ac_table_printable (tbl, ir, 0, "x")) ;
	  aa = arrayReCreate (aa, 1000, WIGGLEPOINT) ;
	  for (jr = 0 ; ir + jr < tbl->rows && ! strcmp (manip, ac_table_printable (tbl, ir+jr, 0, "x")) ; jr++)
	    {
	      x = ac_table_int (tbl, ir+jr, 1, 0) ;
	      y = ac_table_float (tbl, ir+jr, 2, 0) ;
	      wp = arrayp (aa, jr, WIGGLEPOINT) ;
	      wp->x = x ; wp->y = y ;
	    }
	  if (aa && arrayMax (aa))
	    htileWiggleConvert (look, 0 /* was manip, should be a pnx */, ns, aa) ;
	  ir += jr - 1 ;
	}
    }
  else
    {
      int j ;
      int NF = 14 ;
      unsigned int flag[14] =  { PGG_ELF, PGG_ERF, PGG_ELR, PGG_ERR, PGG_nuf, PGG_nur, PGG_ppf, PGG_ppr, PGG_uf, PGG_ur, PGG_endRatioLF, PGG_endRatioRF, PGG_endRatioLR, PGG_endRatioRR  } ;
      const char *suffix[14] = { "u.ELF", "u.ERF", "u.ELR", "u.ERR", "nu.f", "nu.r", "pp.f", "pp.r", "u.f", "u.r", 0, 0, 0, 0} ;
       
      if (look->solexaAll)
	for (pnx = arrp (look->solexaAll, 0, PNX), ns = 0 ; pnx->p ; pnx++, ns++)
	  {
	    char *fNam ;

	    int dn = 0 ;

	    for (j = 0 ; j < NF ; j++)
	      if ((pnx->flag & PGG_endRatios) == 0 && (pnx->flag & flag[j])  ==  flag[j])
		{
		  BOOL ok = FALSE ;
		  if (! ok)
		    {
		      fNam = hprintf (h, "TABIX/%s/R.%s.%s.BF.gz"
				      , pnx->p
				      , name(look->intMap)
				      , suffix[j]
				      ) ;
		      if (filCheckName(fNam, 0, "r"))
			ok = TRUE ;
		    }
		  if (! ok)
		    {
		      dn = 2 ;
		      fNam = hprintf (h, "TABIX/%s/R.%s.%s.BF.gz"
				      , pnx->p
				      , name(look->intMap) + dn
				      , suffix[j]
				      ) ;
		      if (filCheckName(fNam, 0, "r"))
			ok = TRUE ;
		    }
		  if (! ok)
		    {
		      fNam = hprintf (h, "TABIX/%s/%s/R.chrom.%s.BF.gz"
				      , pnx->p
				      , name(look->intMap)
				      , suffix[j]
				      ) ;
		      if (filCheckName(fNam, 0, "r"))
			ok = TRUE ;
		    }
		      

		  if (ok)
		    {
		      Array aa = sxGetWiggleZone (0, fNam, "BF", &(look->solexaStep), name(look->intMap), look->a1, look->a2, h) ;
		      if (aa && arrayMax (aa))
			htileWiggleConvert (look, pnx, ns, aa) ;
		    }
		}
	    if (pnx->flag & PGG_endRatios)
	      htileSolexaEndRatios (look, pnx, ns, NF) ;
	  }
    }
  
  if (look->map->solexa)
    htileSolexaSmoothing (look) ;
  
  ac_free (h) ;
} /* htileSolexaConvert */

/************************************************************/

static void htileDrawSolexa (Htile look, float offset, float probeOffset, ACEOUT ao)
{
  int ns, iSlx ;
  SLX *slx = 0 ;
  PNX *pnx ;
  Array yyy = 0 ;
  float s, olds = 0, x, y, yy,dy,  oldx, zoom, oldWidth = 0, oldWidth1 ;
  int i, np1, xabs, oldxabs = -1 ;
  float mxMx = -1000, mxMy = -1000, mxDelta = 2 ;
  BOOL mxUp = TRUE ;
  float mxMy1 = probeOffset + 1, mxMy2 = look->map->graphHeight - .5 ;
  double zoom1, damper = 0 ;
  double zlog2 = log ((double)2.0) ;
  double zlog10 = log ((double)10.0) ;

  if (!ao)
    oldWidth = graphLinewidth (look->map->lnWidth) ; /* was .2 */
  else
    {
      for (ns = i = 0,  pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p && ns < SOLEXAMAX  ; ns++, pnx++)
	{
	  if (!pnx->signal) continue ;
	  if (look->isSolexaClosed[ns]) continue ;
	  aceOutf(ao, "%s%s", i++ ? "\t" : "", pnx->nam ) ;
	}
      aceOutf (ao, "\n") ;
    }


  if (look->showDot == 3)
    yyy = arrayCreate (50000, float) ;
  for (ns = 0,  pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p  ; ns ++, pnx ++)
    {
      np1 = 0 ; yy = offset ; oldx = -1 ;
      if (!pnx->signal) continue ;
      if (look->isSolexaClosed[ns]) continue ;

      if (!ao)
	graphColor (pnx->col) ;
      else
	aceOutf (ao, "\n%d", look->map->a1 + slx->a1 - 1) ;

      if (!pnx->noSmoothing && look->index && pnx->maxIndex > 1000000)
	zoom1 = 3000000.0/((double)pnx->maxIndex) ; 
      else
	zoom1 = 1 ;

      zoom = look->zoom ;
      olds = 0 ;
      for (iSlx = 1, slx = arrp (look->map->solexa, 1, SLX); iSlx < arrayMax(look->map->solexa) ; iSlx++, slx++)
	{
	  x = TMAP2GRAPH (look->map, slx->a1 + look->solexaStep) ;
	  if (x < 0 || x > look->map-> graphWidth)
	    continue ;
	  if (look->doMask && (slx->flag & MASK_FILTER)) continue ;

	  /* normalize all tag counts to around 3M */
	  if (
	      (slx->signal[ns]==-1000) ||
	      (!pnx->noSmoothing && look->smoothing && slx->gaussSignal[ns] == -1000)
	      )
	    {
	      if (ao) aceOutf(ao, "\t-999999") ;
	      continue ;
	    }
 	  if (!pnx->noSmoothing && look->smoothing) 
	    {
	      s = slx->gaussSignal[ns] ;
	    }
	  else
	    {
	      s = slx->signal[ns] ;
	    }
	  if (look->zlog)
	    { 
	      damper = 5 ; s += damper ; 
	      s *= zoom1 ; damper *= zoom1 ; 
	      if (look->index)
		s = s > damper ? log (s/damper)/zlog10 : 0 ;
	      else
		s = s > damper ? log (s/damper)/zlog2 : 0 ;
	    }
	  else
	    {
	      s *= zoom1 ; damper *= zoom1 ;
	      if (look->index)
		{ s /= pnx->s99 ; damper /= pnx->s99 ; }
	    }
	  if (ao)
	    aceOutf (ao, "\t%g", s) ;
	  else
	    {
	      xabs =  uToXabs(x) ;
	      if (s > olds) olds = s ;
	      if (xabs != oldxabs)
		{
		  y = offset - zoom * olds ;
		  if (y < 2) y = 2 ;
		  if (look->showDot == 1)
		    graphCircle (x, y, .2) ;
		  else if (look->showDot == 2)
		    graphLine (x, offset, x, y) ;
		  else if (look->showDot == 3)
		    {
		      if (xabs != oldxabs)
			{
			  dy = array (yyy, xabs, float) ; 
			  if (dy < offset)
			    graphLine (x, offset - dy, x, y - dy) ;
			  array (yyy, xabs, float) += zoom * olds ;
			}
		    } 
		  else if (oldx > 0 &&  np1)
		    {
		      graphLine (oldx, yy, x, y) ;
		    }
		  np1++ ;
		  oldxabs = xabs ;
		  oldx = x ;
		  yy = y ;
		  olds = 0 ;
		}
	    }

	  /* show lines at extremal points */
	  if (!ao && look->showExtrema && (!strcmp (pnx->nam,"S=31ES") || (ns && look->mxShowNs == ns)))
	    {
	      if ((mxUp && mxMy < y) || (!mxUp && mxMy > y)) { mxMx = x ; mxMy = y ; }
	      else if ((mxUp && mxMy > y + mxDelta) ||
		       (!mxUp && mxMy < y - mxDelta)
		       ) /* report previous */
		{
		  mxUp = !mxUp ;
		  if (mxMx > 0) 
		    {
		      oldWidth1 = graphLinewidth (.2) ;
		      graphLine (mxMx, mxMy1, mxMx, mxMy2) ;
		      graphLinewidth (oldWidth1) ;
		    }
		  mxMx = x ; mxMy = y ;
		}
	    } 	    /* end of drawing extremal points */
	}
    }
  if (!ao)
    {
      graphLinewidth (oldWidth) ;
      graphColor (BLACK) ;
    }
  else
    aceOutf (ao, "\n") ;
  arrayDestroy (yyy) ;
} /* htileDrawSolexa */

/************************************************************/
/* open/close signal tracks */
static void htileButtonAction (void *v)
{
  int ns = assInt (v) ;
  BOOL newClosed ;
  Htile look = currentHtile("htileButtonAction") ;

  if (look->isClosed[ns])
    newClosed = look->isClosed[ns] = 0 ;
  else
    newClosed = look->isClosed[ns] = 2 ;
  look0->isClosed[ns] = look->isClosed[ns] ;
  if (look->ratio)
    {
      look->isClosed[ns+1] = look->isClosed[ns] ;
    }
  if (!newClosed && look->smoothing && !(1 && look->isSmoothed[ns]) )
    {
      htileSmoothing  (look) ;
    }
  look->map->mapDraw () ;
} /* htileButtonAction */

/************************************************************/
/* open/close signal tracks */
static void htileSolexaButtonAction (void *v)
{
  int ns = assInt (v) ;
  BOOL newSolexaClosed ;
  Htile look = currentHtile("htileSolexaButtonAction") ;

  if (! look) return ;
  newSolexaClosed = look->isSolexaClosed[ns] = look->isSolexaClosed[ns] ? 0 : 2 ;
  look0->isSolexaClosed[ns] = look->isSolexaClosed[ns] ;
  if (look->ratio)
    {
      look->isSolexaClosed[ns+1] = look->isSolexaClosed[ns] ;
    }
  if (!newSolexaClosed && look->smoothing && !(1 && look->isSolexaSmoothed[ns]) )
    {
      if (look->map->solexa)
	htileSolexaSmoothing (look) ;
    }
  look->map->mapDraw () ;
} /* htileSolexaButtonAction */

/************************************************************/

static void htileDrawTileGroupButtons (Htile look, float *offsetp)
{
  COLOUROPT bb[NSMAX] ;
  int box, ns, ng=0, nm ;
  PGG *pgg ;
  PNS *pns ;
  PNX *pnx ;
  FREEOPT *opt, *opts ;
  static int nmm = 0 ;
  BOOL isClosed ;

  for (ng = 0, pgg = pggGroups ; pgg->flag ; ng++, pgg++)
    {
      isClosed = TRUE ;

      for (ns = 0, pns = pnsAll ; pns->p ; ns += look->ratio ? 2 : 1, pns += look->ratio ? 2 : 1)
	{
	  if (! (pns->flag & pgg->flag))
	    continue ;
	  if (!look->isClosed[ns])
	    isClosed = FALSE ;
	}
      for (ns = 0,  pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p ; ns += look->ratio ? 2 : 1, pnx += look->ratio ? 2 : 1)
	{
	  if (! (pnx->flag & pgg->flag))
	    continue ;
	  if (!look->isSolexaClosed[ns])
	    isClosed = FALSE ;
	}
      
      bb[ng].f = htileGroupActionButton ;
      bb[ng].text = hprintf (look->h, "%s...", pgg->nam) ;
      bb[ng].fg = (int)BLACK ;
      if (isClosed)
	bb[ng].bg = PALEGRAY ;
#ifdef BOUBOU_DEF
      else if (pgg->flag & PGG_PA)
	bb[ng].bg = PALERED ;
      else if (pgg->flag & PGG_NUC)
	bb[ng].bg = PALERED ;
      else if (pgg->flag & PGG_TOTAL)
	bb[ng].bg =  PALERED ;
      else if (pgg->flag & PGG_S)
	bb[ng].bg =  PALERED ;
      else if (pgg->flag & PGG_G1)
	bb[ng].bg = PALERED ;
#endif
      else
	bb[ng].bg = PALEGREEN ;
      bb[ng].arg = assVoid (10000 * ng + (isClosed ? 1 : 2)) ;
      bb[ng].next = 0 ;
    }

  bb[ng].text = 0 ;
  box = graphColouredButtons (bb, 2, offsetp, look->map->graphWidth) ; 

  for (ng = 0, pgg = pggGroups ; ng < 20 && pgg->flag ; ng++, pgg++)
    {
      nm = 0 ;
      for (ns = 0, pns = pnsAll ; pns->p ; ns += look->ratio ? 2 : 1, pns += look->ratio ? 2 : 1)
	{
	  if (! (pns->flag & pgg->flag))
	    continue ;
	  /*
	    if (pgg->flag != PGG_SOLEXA && !look->isClosed[ns]) continue ;
	  */
	  nm++ ;	      
	}
      for (ns = 0,  pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p ; ns += look->ratio ? 2 : 1, pnx += look->ratio ? 2 : 1)
	{
	  if (! (pnx->flag & pgg->flag))
	    continue ;
	  if (!pns->signal)
	    continue ;
	  /*
	    if (pgg->flag == PGG_SOLEXA && !look->isSolexaClosed[ns]) continue ;
	  */
	  nm++ ;	      
	}

      if (!nm)
	continue ;
      opts = halloc ((nm+7) * sizeof (FREEOPT), look->h) ;
      opt = opts ;
      opt->key = /* nm + */  5 ; opt->text = hprintf (look->h, "hfg_%d", nmm++) ; opt++ ;
      opt->key = 10000 * ng + 1 ; opt->text = "Open all" ; opt++ ;
      opt->key = 10000 * ng + 2 ; opt->text = "Close all" ; opt++ ;
      opt->key = 1000000 + 35 ; opt->text = "Mask 35" ; opt++ ;
      opt->key = 1000000 + 25 ; opt->text = "Mask 25" ; opt++ ;
      /*       opt->key = 1000000 + 18 ; opt->text = "Mask 18" ; opt++ ; */
      opt->key = 1000000 + 0  ; opt->text = "UnMask" ; opt++ ;
      graphBoxFreeMenu (box + ng, htileGroupAction, opts) ;
    }
} /* htileDrawTileGroupButtons */

/************************************************************/

static void solexaColourChange (KEY key, int box)
{
  Htile look = currentHtile("solexaColourChange") ;
  int nt, n = box - look->map->bbBox ;

  if (n >= 0 && n < SOLEXAMAX)
    {
      nt = assInt (look->map->bb[n].arg) ;  
      if (nt >= 0 && look->solexaAll && nt < arrayMax (look->solexaAll))
	arrp (look->solexaAll, nt, PNX)->col = key ;
    }
  look->map->mapDraw () ;
}

/************************************************************/

static void htileDrawTileButtons (Htile look, float *offsetp)
{
  COLOUROPT *bb = look->map->bb ;
  int ii, ns, nt1=0, nt ;
  PNS *pns ;
  PNX *pnx ;
  const char *ccp ;

  for (ns = 0, pns = pnsAll ; pns->p && nt1 < SOLEXAMAX ; ns += look->ratio ? 2 : 1, pns += look->ratio ? 2 : 1)
    {
      ccp =  look->ratio ? pns->nam2 : pns->nam ;
      if (*ccp == '.') continue ; 
      if (look->isClosed[ns] == 1) continue ;
      bb[nt1].f = htileButtonAction ; 
      bb[nt1].text = ccp ;
      bb[nt1].fg = (int)BLACK ;
      bb[nt1].bg = look->isClosed[ns] ? PALEGRAY : pnsAll[ns].col ;
      bb[nt1].arg = assVoid (ns) ;
      bb[nt1].next = 0 ;
      nt1++ ;
    }
  for (ccp = 0, pnx = arrp(look->solexaAll, 0, PNX), nt = 0 ; pnx->p && nt1 < SOLEXAMAX ; nt ++, pnx++)
    {
      if (!pnx->signal)
	continue ;
      if (ccp == pnx->nam)
	continue ;
      ccp =  pnx->nam ;
      if (*ccp == '.') continue ;
      if (look->isSolexaClosed[nt] == 1) continue ;
      bb[nt1].f = htileSolexaButtonAction ;  
      bb[nt1].text = ccp ;
      bb[nt1].fg = (int)BLACK ;
      bb[nt1].bg = look->isSolexaClosed[nt] ? PALEGRAY : pnx->col ;
      bb[nt1].arg = assVoid (nt) ;
      bb[nt1].next = 0 ;
      nt1++ ;
    }
  bb[nt1].text = 0 ;
  if (! isGifDisplay)
    {
      look->map->bbBox = graphColouredButtons (bb, 2, offsetp, look->map->graphWidth) ;
      
      for (ii = 0 ; ii < nt1 ; ii++)
	graphBoxFreeMenu(look->map->bbBox + ii, solexaColourChange, graphColors);
    }
  else
    {
      float oldh = graphTextHeight (1.5) ;
      look->map->bbBox = graphColouredButtons (bb, 2, offsetp, look->map->graphWidth) ;
      graphTextHeight (oldh) ;
    }
}

/************************************************************/

static void htileShowSignal (void *v)
{
  int ns, showRatio = assInt (v) ;
  Htile look = currentHtile("htileShowSignal") ;
  BOOL needSmoothing = FALSE ;
  int x ;
  PNS *pns ;
  PNX *pnx ;
  
  look->ratio = showRatio ;
  look0->ratio = showRatio ;
  switch (showRatio) 
    {
    case 0:
      for (ns = 0, pns = pnsAll ; pns->p ; ns += 2, pns += 2)
	look0->isClosed[ns+1] = look->isClosed[ns+1] = look->isClosed[ns] ;
      for (ns = 0, pnx =  arrp(look->solexaAll, 0, PNX) ; pnx->p ; ns += 2, pnx += 2)
	look0->isSolexaClosed[ns+1] = look->isSolexaClosed[ns+1] = look->isSolexaClosed[ns] ;
      break ;
    case 1:
      for (ns = 0, pns = pnsAll ; pns->p ; ns += 2, pns += 2)
	{
	  x = look->isClosed[ns] | look->isClosed[ns+1] ;
	  if (x==3) x = 2 ;
	  look0->isClosed[ns] = look->isClosed[ns] = x ;
	}
      for (ns = 0, pnx =  arrp(look->solexaAll, 0, PNX) ; pnx->p ; ns += 2, pnx += 2)
	{
	  x = look->isSolexaClosed[ns] | look->isSolexaClosed[ns+1] ; 
	  if (x==3) x = 2 ;
	  look0->isSolexaClosed[ns] = look->isSolexaClosed[ns] = x ;
	}
      break ;
    case 2:
      for (ns = 0, pns = pnsAll ; pns->p ; ns ++, pns++)
	{
	  if (
	      ((pns->flag & PGG_ELF) && ((pns+1)->flag & PGG_ERF)) ||
	      ((pns->flag & PGG_ELR) && ((pns+1)->flag & PGG_ERR)) 
	      )
	  look0->isClosed[ns+1] = look->isClosed[ns+1] = look->isClosed[ns] ;
	}
      for (ns = 0, pnx =  arrp(look->solexaAll, 0, PNX) ; pnx->p ; ns ++, pnx++)
	{
	  if (
	      ((pnx->flag & PGG_ELF) && ((pnx+1)->flag & PGG_ERF)) ||
	      ((pnx->flag & PGG_ELR) && ((pnx+1)->flag & PGG_ERR)) 
	      )
	  look0->isSolexaClosed[ns+1] = look->isSolexaClosed[ns+1] = look->isSolexaClosed[ns] ;
	}
      break ;
    }
    
  if (look->smoothing)
    {
      needSmoothing = FALSE ;
      for (ns = 0, pns = pnsAll ; pns->p ; ns++, pns++)
	if (!look->isClosed[ns] && !look->isSmoothed[ns])
	  needSmoothing = TRUE ;
      if (needSmoothing)
	htileConvert (look, FALSE, 0) ;

      needSmoothing = FALSE ;
      for (ns = 0, pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p ; ns++, pnx++)
	if (!look->isSolexaClosed[ns] && !look->isSolexaSmoothed[ns])
	  needSmoothing = TRUE ;
      if (needSmoothing)
	htileSolexaConvert (look, FALSE, 0) ;
    }
  
  look->map->mapDraw () ;
} /*  htileShowSignal */

/************************************************************/

static void htileFilterProbes (void *v)
{
  int ns = assInt (v) ;
  Htile look = currentHtile("htileFilterProbes") ;

  switch (ns)
    {
    case 60:
      look->showExtrema = FALSE ;
      look0->showExtrema = FALSE ;
      break ;
    case 61:
      look->showExtrema = TRUE ;
      look0->showExtrema = TRUE ;
      break ;
    case 90: /* show lines */
      look->showDot = 0 ;
      look0->showDot = 0 ;
      break ;
    case 91: /* show dots */
      look->showDot = 1 ;
      look0->showDot = 1 ;
      break ;
    case 92: /* show bars */
      look->showDot = 2 ;
      look0->showDot = 2 ;
      break ;
    case 93: /* show pileup */
      look->showDot = 3 ;
      look0->showDot = 3 ;
      break ;
    case 200: /* all probes */
     if (look->exonicFilter)
       {
	 look->exonicFilter = 0 ;
	 look0->exonicFilter = 0 ;
	 htileConvert (look, FALSE, 0) ;
       }
      break ;
    case 201: /* EXONIC_PROBE */
     if (look->exonicFilter != 1)
       {
	 look->exonicFilter = 1 ;
	 look0->exonicFilter = 1 ;
	 htileConvert (look, FALSE, 0) ;
       }
      break ;
    case 202: /* NEW_EXONIC */
     if (look->exonicFilter != 2)
       {
	 look->exonicFilter = 2 ;
	 look0->exonicFilter = 2 ;
	 htileConvert (look, FALSE, 0) ;
       }
      break ;
    case 203: /* EXONIC | NEW_EXONIC */
     if (look->exonicFilter != 3)
       {
	 look->exonicFilter = 3 ;
	 look0->exonicFilter = 3 ;
	 htileConvert (look, FALSE, 0) ;
       }
      break ;
    case 204: /* INTRONIC_PROBE */
     if (look->exonicFilter != 4)
       {
	 look->exonicFilter = 4 ;
	 look0->exonicFilter = 4 ;
	 htileConvert (look, FALSE, 0) ;
       }
      break ;
    case 205: /* INTERGENIC_PROBE */
     if (look->exonicFilter != 5)
       {
	 look->exonicFilter = 5 ;
	 look0->exonicFilter = 5 ;
	 htileConvert (look, FALSE, 0) ;
       }
      break ;
    case 300: /* any TM */
     if (look->tmFilter)
       {
	 look->tmFilter = 0 ;
	 look0->tmFilter = 0 ;
	 htileConvert (look, FALSE, 0) ;
       }
      break ;
    case 301: /* Tm < 72 */
     if (look->tmFilter != 1)
       {
	 look->tmFilter = 1 ;
	 look0->tmFilter = 1 ;
	 htileConvert (look, FALSE, 0) ;
       }
     break ;
    case 302: /* Tm > 75 */
      if (look->tmFilter != 2)
	{
	  look->tmFilter = 2 ;
	  look0->tmFilter = 2 ;
	  htileConvert (look, FALSE, 0) ;
	}
      break ;
    case 400: /* multi mapping */
      if (look->unique)
	{
	  look->unique = 0 ;
	  look0->unique = 0 ;
	  htileSolexaConvert (look, FALSE, 0) ;
	}
      break ;
    case 401: /* unique */
       if (look->unique != 1)
	{
	  look->unique = 1 ;
	  look0->unique = 1 ;
	  htileSolexaConvert (look, FALSE, 0) ;
	}
      break ;
     case 402: /* partial */
       if (look->unique != 2)
	{
	  look->unique = 2 ;
	  look0->unique = 2 ;
	  htileSolexaConvert (look, FALSE, 0) ;
	}
      break ;
    case 500: /* unique mapping */
      if (look->rejectAmbiguous)
	{
	  look->rejectAmbiguous = 0 ;
	  look0->rejectAmbiguous = 0 ;
	  htileConvert (look, FALSE, 0) ;
	}
      break ;
    case 501: /* keep local repeats */
       if (look->rejectAmbiguous != 1)
	{
	  look->rejectAmbiguous = 1 ;
	  look0->rejectAmbiguous = 1 ;
	  htileConvert (look, FALSE, 0) ;
	}
      break ;
    case 502: /* single chromo */
       if (look->rejectAmbiguous != 2)
	{
	  look->rejectAmbiguous = 2 ;
	  look0->rejectAmbiguous = 2 ;
	  htileConvert (look, FALSE, 0) ;
	}
      break ;
    case 503: /* any mapping */
       if (look->rejectAmbiguous != 3)
	{
	  look->rejectAmbiguous = 3 ;
	  look0->rejectAmbiguous = 3 ;
	  htileConvert (look, FALSE, 0) ;
	}
      break ;
    case 600: /* any signal or probe motif */
       if (look->filter99 != 0)
	{
	  look->filter99 = 0 ;
	    look0->filter99 = 0 ;
	  htileConvert (look, FALSE, 0) ;
	}
      break ;
    case 601: /* filter probe with non central G1 */
       if (look->filter99 != 1)
	{
	  look->filter99 = 1 ;
	  look0->filter99 = 1 ;
	  htileConvert (look, FALSE, 0) ;
	}
      break ;
    case 602: /* filter probes with suspect motif */
       if (look->filter99 != 2)
	{
	  look->filter99 = 2 ;
	  look0->filter99 = 2 ;
	  htileConvert (look, FALSE, 0) ;
	}
      break ;
    case 603: /* filter signal in 99% quantile */
       if (look->filter99 != 3)
	{
	  look->filter99 = 3 ;
	  look0->filter99 = 3 ;
	  htileConvert (look, FALSE, 0) ;
	}
      break ;
    case 40:
      look->romainSmoothing = 1 ;
      look0->romainSmoothing = 1 ;
      htileConvert (look, FALSE, 0) ;
      break ;
    case 41:
      look->romainSmoothing = 2 ;
      look0->romainSmoothing = 2 ;
      htileConvert (look, FALSE, 0) ;
      break ;
    case 42:
      look->romainSmoothing = 0 ;
      look0->romainSmoothing = 0 ;
      htileConvert (look, FALSE, 0) ;
      break ;
    case 70:
      look->smoothing = 0 ;
      look0->smoothing = 0 ;  
      htileConvert (look, FALSE, 0) ;
      break ;
    case 71:
      look->smoothing = 1 ;
      look0->smoothing = 1 ;
      htileConvert (look, FALSE, 0) ;
      break ;
    case 72: 
      look->smoothing = 2 ;
      look0->smoothing = 2 ;
      htileConvert (look, FALSE, 0) ;
      break ;
    case 73: /* apparently never used */
     look->smoothing = 3 ;
      look0->smoothing = 3 ;
      htileConvert (look, FALSE, 0) ;
      break ;
    case 74:
     look->smoothing = 4 ;
      look0->smoothing = 4 ;
      htileConvert (look, FALSE, 0) ;
      break ;
    case 75:
      look->smoothing = 5 ;
      look0->smoothing = 5 ;
      htileConvert (look, FALSE, 0) ;
      break ;
    case 76:
     look->smoothing = 6 ;
      look0->smoothing = 6 ;
      htileConvert (look, FALSE, 0) ;
      break ;
    case 77:
      look->smoothing = 7 ;
      look0->smoothing = 7 ;
      htileConvert (look, FALSE, 0) ;
      break ;
    case 8:
      look->zoom *= 2 ;
      look0->zoom *= 2 ;
      break ;
    case 9:
      if (look->zoom > 0) look->zoom /= 2 ;
      look0->zoom = look->zoom ;
      break ;
    case 10:
      {
	int w, i, i0 = 1, j, jj[3] = {1,2,5} ;
	w = look->gaussWidth ;
	for (i = 1 ; i0 && i < 1000000 ; i *= 10)
	  for (j = 0 ; i0 && j < 3 ; j++)
	    if (w == i * jj[j]) 
	      { w = j < 2 ? i * jj[j+1] : 10*i ; i0 = 0 ;}
	look->gaussWidth = look0->gaussWidth = w ;
      }
      htileConvert (look, FALSE, 0) ;
      break ;
    case 11:
      {
	int w, i, i0 = 1, j, jj[3] = {1,2,5} ;
	w = look->gaussWidth ;
	for (i = 1 ; i0 && i < 1000000000 ; i *= 10)
	  for (j = 0 ; i0 && j < 3 ; j++)
	    if (w == i * jj[j]) 
	      { w = j > 0 ? i * jj[j-1] : i/2 ; i0 = 0 ; }
	if (w < 10) w = 10 ;
	look->gaussWidth = look0->gaussWidth = w ;
      }
      htileConvert (look, FALSE, 0) ;
      break ;

    case 9000:
      if (look->smoothing)
	{
	  look->smoothing = 0 ;
	  htileConvert (look, FALSE, 0) ;
	}
      else
	return ;
      break ;
    case 9002:
    case 9005:
    case 9011:
    case 9012:
    case 9015:
    case 9021:
    case 9022:
    case 9025:
    case 9031:
    case 9032:
    case 9035:
    case 9041:
    case 9042:
    case 9045:
    case 9051:
    case 9052:
    case 9055:
    case 9061:
    case 9062:
    case 9065:
      {  /* 90...number of zeroes...digit */
	int w = 1 ;
	switch (ns)
	  {
	  case 9002: w = 2 ; break ;
	  case 9005: w = 5 ; break ;
	  case 9011: w = 10 ; break ;
	  case 9012: w = 20 ; break ;
	  case 9015: w = 50 ; break ;
	  case 9021: w = 100 ; break ;
	  case 9022: w = 200 ; break ;
	  case 9025: w = 500 ; break ;
	  case 9031: w = 1000 ; break ;
	  case 9032: w = 2000 ; break ;
	  case 9035: w = 5000 ; break ;
	  case 9041: w = 10000 ; break ;
	  case 9042: w = 20000 ; break ;
	  case 9045: w = 50000 ; break ;
	  case 9051: w = 100000 ; break ;
	  case 9052: w = 200000 ; break ;
	  case 9055: w = 500000 ; break ;
	  case 9061: w = 1000000 ; break ;
	  case 9062: w = 2000000 ; break ;
	  case 9065: w = 5000000 ; break ;
	  }
	if (!look->smoothing) look->smoothing = look0->smoothing ;
	look->gaussWidth = look0->gaussWidth = w ;
      }
      htileConvert (look, FALSE, 0) ;
      break ;
    case 12: /* raw count */
      look->index = FALSE ;
      look0->index = FALSE ;
      look->zlog = FALSE ;
      look0->zlog = FALSE ;
      break ;
    case 13: /* Count per Gb */
      look->index = TRUE ;
      look0->index = TRUE ;
      look->zlog = FALSE ;
      look0->zlog = FALSE ;
      break ;
    case 14:  /* Log2 raw */
      look->index = FALSE ;
      look0->index = FALSE ;
      look->zlog = TRUE ;
      look0->zlog = TRUE ;
      break ;
    case 15:  /* Log10 raw */
      look->index = TRUE ;
      look0->index = TRUE ;
      look->zlog = TRUE ;
      look0->zlog = TRUE ;
      break ;
    }
  look->map->mapDraw () ;
} /*  htileFilterProbes */

/************************************************************/

static FREEOPT htileSmootherMenu[] = {
  { 8, "Smoother menu"},
  {70, "No smoothing"},
  {71, "Gaussian smoothing"},
  {72, "Parabolic smoothing"},
  {73, "Square smoothing"},
  {74, "Median smoothing"},
  {75, "Gene smoothing up"},
  {76, "Gene smoothing down"},
  {77, "Gene smoothing merged"}
} ;

static void htileSmootherAction (KEY key, int box)
{ 
  htileFilterProbes (assVoid(key)) ;
}

/************************************************************/

static FREEOPT htileRatioMenu[] = {
  { 3, "Ratio menu"},
  {1, "Ratio"},
  {2, "Ends"},
  {0, "Signal"}
} ;

static void htileRatioAction (KEY key, int box)
{ 
  htileShowSignal (assVoid(key)) ;
}

/************************************************************/

static FREEOPT htileVppMenu[] = {
  { 19, "Ratio menu"},
  {9000, "No smoothing"},
  /*
    {9002, "  2"},
    {9005, "  5"},
  */
  {9011, " 10"},
  {9012, " 20"},
  {9015, " 50"},
  {9021, "100"},
  {9022, "200"},
  {9025, "500"},
  {9031, "  1 k"},
  {9032, "  2 k"},
  {9035, "  5 k"},
  {9041, " 10 k"},
  {9042, " 20 k"},
  {9045, " 50 k"},
  {9051, "100 k"},
  {9052, "200 k"},
  {9055, "500 k"},
  {9061, "  1 M"},
  {9062, "  2 M"},
  {9065, "  5 M"}
} ;

static void htileVppAction (KEY key, int box)
{ 
  htileFilterProbes (assVoid(key)) ;
}

/************************************************************/

static FREEOPT htileValueMenu[] = {
  { 4, "Value menu"},
  {12, "Raw count"},
  {13, "Per Gb"},
  {14, "Log2 raw"},
  {15, "Log10 raw"}
} ;

static void htileValueAction (KEY key, int box)
{ 
  htileFilterProbes (assVoid(key)) ;
}

/************************************************************/

static FREEOPT htileDotMenu[] = {
  { 4, "Dot menu"},
  {90, "Line"},
  {91, "Dot"},
  {92, "Bar"},
  {93, "Stack"}
} ;

static void htileDotAction (KEY key, int box)
{ 
  htileFilterProbes (assVoid(key)) ;
}

/************************************************************/

static FREEOPT htileExonicMenu[] = {
  { 6, "Exonic menu"},
  {200, "Any Probe"},
  {201, "Exonic"},
  {202, "New exonic"},
  {203, "Exonic + new exonic"},
  {204, "Intronic"},
  {205, "Intergenic"},
} ;

static void htileExonicAction (KEY key, int box)
{ 
  htileFilterProbes (assVoid(key)) ;
}

/************************************************************/

static FREEOPT htileUniqueMenu[] = {
  { 3, "Unique menu"},
  {400, "Show non-unique tags"},
  {401, "Show tags mapping uniquely"},
  {402, "Show tags with partial of ambiguous mapping"}
} ;

static void htileUniqueAction (KEY key, int box)
{ 
  htileFilterProbes (assVoid(key)) ;
}

/************************************************************/

static FREEOPT htileAmbiguousMenu[] = {
  { 4, "Ambiguous menu"},
  {500, "Reject any probe mapping at multiple sites"},
  {501, "Keep probes mapping in close repeats"},
  {502, "Keep probes mapping on a single chromosome"},
  {503, "Keep all probes"}
} ;

static void htileAmbiguousAction (KEY key, int box)
{ 
  htileFilterProbes (assVoid(key)) ;
}

/************************************************************/

static FREEOPT htile99Menu[] = {
  { 4, "99 menu"},
  {600, "Keep all probes"},
  {601, "Reject probes with average G1 signal outside 5-95%"},
  {602, "Reject probes with possibly biased sequence"},
  {603, "Reject top 1% signals"}
} ;

static void htile99Action (KEY key, int box)
{ 
  htileFilterProbes (assVoid(key)) ;
}

/************************************************************/

static FREEOPT htileTmMenu[] = {
  { 3, "TM menu"},
  {300, "Any TM"},
  {301, "TM < 72.2"},
  {302, "TM > 75"}
} ;

static void htileTmAction (KEY key, int box)
{ 
  htileFilterProbes (assVoid(key)) ;
}

/************************************************************/

static FREEOPT htileFilterMenu[] = {
  { 8, "Rnai menu"},
  {200, "Any probe"},
  {201, "Exonic probe"},
  {202, "New exonic"},
  {203, "Intronic"},
  {204, "Intergenic"},
  {207, "Sequence filter"},
  {208, "Just TM < 72.2"},
  {209, "Just TM > 75"}
} ;

static void htileFilterAction (KEY key, int box)
{ 
  htileFilterProbes (assVoid(key)) ;
}

/************************************************************/
static void htileGroupActionDo (KEY key)
{
  int ng = key / 10000 ;
  int ns, nn = (key % 10000) ;
  PGG *pgg ;
  PNS *pns ;
  PNX *pnx ;
  BOOL needSmoothing = FALSE ;
  Htile look = currentHtile("htileGroupAction") ;
 
  pgg = pggGroups + ng ;

  if (key >= 1000000 && key < 1000000 + 256)
    {
      switch (key - 1000000)
	{
	case 18:
	case 25:
	case 35:
	  if (look->doMask != key - 1000000)
	    {
	      look->isMasked = 0 ;
	      look->doMask = key - 1000000 ;
	      if (look->map->solexa)
		htileSolexaSmoothing (look) ;
	    }
	  break ;
	default:
	  look->doMask = 0 ;
	  break ;
	}
      goto done ;
    }

    for (ns = 0, pns = pnsAll ; pns->p ; ns += look->ratio ? 2 : 1, pns += look->ratio ? 2 : 1)
      {
	if (! (pns->flag & pgg->flag))
	  continue ; 
	if (nn == 1) 
	  look0->isClosed[ns] = look->isClosed[ns] = FALSE ;
	if (nn == 2) 
	  look0->isClosed[ns] = look->isClosed[ns] = TRUE ;
      }
    for (ns = 0, pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p ; ns += look->ratio ? 2 : 1, pnx += look->ratio ? 2 : 1)
      {
	if (! (pnx->flag & pgg->flag))
	  continue ; 

	if (nn == 1) 
	  look0->isSolexaClosed[ns] = look->isSolexaClosed[ns] = FALSE ;
	if (nn == 2) 
	  look0->isSolexaClosed[ns] = look->isSolexaClosed[ns] = TRUE ;
      }
  if (look->smoothing)
    {
      needSmoothing = FALSE ;
      for (ns = 0, pns = pnsAll ; pns->p ; ns++, pns++)
	if (!look->isClosed[ns] && !look->isSmoothed[ns])
	  needSmoothing = TRUE ;
      if (needSmoothing)
	htileSmoothing (look) ;
      
      needSmoothing = FALSE ;
      for (ns = 0, pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p ; ns++, pnx++)
	if (!look->isSolexaClosed[ns] && !look->isSolexaSmoothed[ns])
	  needSmoothing = TRUE ;
      if (needSmoothing && look->map->solexa)
	 htileSolexaSmoothing (look) ;
    }
  
 done:
  look->map->mapDraw () ;
  
  return ;
} /* htileGroupActionDo */
static void htileGroupActionButton (void *v)
{ 
  KEY key = assInt (v) ;
  htileGroupActionDo (key) ;
}

static void htileGroupAction (KEY key, int box)
{ 
  htileGroupActionDo (key) ;
} /* htileGroupAction */

/************************************************************/

static void htileDrawActionButtons (Htile look, float *offsetp)
{
  COLOUROPT bb[80] ;
  int box, ns = 0
    , nsSmoother = 0, nsVpp = 0, nsVmm = 0, nsRatio = 0, nsValue = 0
    , nsDot = 0, nsExonic = 0, nsUnique = 0, nsAmbiguous = 0, nsTm = 0
    , nsFilter = 0, ns99 = 0 
    ;
  static char buf1[64] ;
  
  if (0)
    {
      bb[ns].f = htileFilterProbes;
      bb[ns].text = look->showExtrema ? "Show extrema" : "Hide extrema" ;
      bb[ns].fg = BLACK ;
      bb[ns].bg = PALEBLUE ;
      bb[ns].arg =  look->showExtrema ? assVoid(60)  : assVoid(61) ;
      bb[ns].next = 0 ;
      ns++ ;
    }

  bb[ns].f = htileFilterProbes;
  bb[ns].text = "y++" ;
  bb[ns].fg = BLACK ;
  bb[ns].bg = PALEBLUE ;
  bb[ns].arg = assVoid(8) ;
  bb[ns].next = 0 ;
  ns++ ;

  bb[ns].f = htileFilterProbes;
  bb[ns].text = "y--" ;
  bb[ns].fg = BLACK ;
  bb[ns].bg = PALEBLUE ;
  bb[ns].arg = assVoid(9) ;
  bb[ns].next = 0 ;
  ns++ ;

  if (1)
    {
      int w = look->gaussWidth ;
      if (w < 1000)
	sprintf(buf1, "v++:%d bp", look->gaussWidth) ; 
      else if (w < 1000000)
	sprintf(buf1, "v++:%d kb", look->gaussWidth/1000) ; 
      else 
	sprintf(buf1, "v++:%d Mp", look->gaussWidth/1000000) ; 
    }
   
  bb[ns].f = htileFilterProbes;
  bb[ns].text = buf1 ;
  bb[ns].fg = BLACK ;
  bb[ns].bg = PALEBLUE ;
  bb[ns].arg = assVoid(10) ;
  bb[ns].next = 0 ;
  nsVpp = ns ;
  ns++ ;

  bb[ns].f = htileFilterProbes;
  bb[ns].text = "v--" ;
  bb[ns].fg = BLACK ;
  bb[ns].bg = PALEBLUE ;
  bb[ns].arg = assVoid(11) ;
  bb[ns].next = 0 ;
  nsVmm = ns ;
  ns++ ;

  bb[ns].f = htileShowSignal ;
  bb[ns].fg = BLACK ;
  bb[ns].bg = PALEBLUE ;
  switch (look->ratio)
    {
      case 0: 
	bb[ns].text = "Signal" ;
	bb[ns].arg = assVoid(1) ; 
	break ;
      case 1: 
	bb[ns].text = "Ratio" ;
	bb[ns].arg = assVoid(0) ; 
	break ;
    }
  bb[ns].next = 0 ;
  nsRatio = ns ;
  ns++ ;

  bb[ns].f =  htileFilterProbes; ;
  bb[ns].fg = BLACK ;
  bb[ns].bg = PALEBLUE ;
  switch (look->index + 2*look->zlog)
    {
    case 0: bb[ns].arg = assVoid(13) ; bb[ns].text = "Raw count" ; break ;
    case 1: bb[ns].arg = assVoid(14) ; bb[ns].text = "Per Gb" ; break ;
    case 2: bb[ns].arg = assVoid(15) ; bb[ns].text = "Log2 raw" ; break ;
    case 3: bb[ns].arg = assVoid(12) ; bb[ns].text = "Log10 raw" ; break ;
    }
  bb[ns].next = 0 ;
  nsValue = ns ;
  ns++ ;

  bb[ns].f = htileFilterProbes; 
  switch (look->showDot)
    {
    case 0:
      bb[ns].text = "Line" ;
      bb[ns].arg =  assVoid(91) ;
      break ;
    case 1:
      bb[ns].text = "Dot" ;
      bb[ns].arg =  assVoid(92) ;
      break ;
    case 2:
      bb[ns].text = "Bar" ;
      bb[ns].arg =  assVoid(93) ;
      break ;
    case 3:
      bb[ns].text = "Stack" ;
      bb[ns].arg =  assVoid(90) ;
      break ;
    }
  bb[ns].fg = BLACK ;
  bb[ns].bg = PALEBLUE ;
  bb[ns].next = 0 ;
  nsDot = ns ;
  ns++ ;

  bb[ns].f = htileFilterProbes;
  switch (look->exonicFilter)
    {
    case 0: 
      bb[ns].text = "Probes" ;
      bb[ns].arg = assVoid(201) ;
      break ;
    case 1:
      bb[ns].text = "Exonic" ;
      bb[ns].arg = assVoid(202) ;
      break ;
    case 2:
      bb[ns].text = "New exonic" ;
      bb[ns].arg = assVoid(203) ;
      break ;
    case 3:
      bb[ns].text = "All exonic" ;
      bb[ns].arg = assVoid(204) ;
      break ;
    case 4:
      bb[ns].text = "Intronic" ;
      bb[ns].arg = assVoid(205) ;
      break ;
    case 5:
      bb[ns].text = "Intergenic" ;
      bb[ns].arg = assVoid(200) ;
      break ;
    }
  bb[ns].fg = BLACK ;
  bb[ns].bg = PALEBLUE ;  
  bb[ns].next = 0 ;
  nsExonic = ns ;
  ns++ ;


  bb[ns].f = htileFilterProbes;
  switch (look->unique)
    {
    case 0: 
      bb[ns].text = "Multi" ;
      bb[ns].arg = assVoid(401) ;
      break ;
    case 1:
      bb[ns].text = "Unique" ;
      bb[ns].arg = assVoid(402) ;
      break ;
     case 2:
      bb[ns].text = "Partial" ;
      bb[ns].arg = assVoid(400) ;
      break ;
    }
  bb[ns].fg = BLACK ;
  bb[ns].bg = PALEBLUE ;
  bb[ns].next = 0 ;
  nsUnique = ns ;
  ns++ ;

#ifdef BOU_DEF

  bb[ns].f = htileFilterProbes;
  switch (look->rejectAmbiguous)
    {
    case 0: 
      bb[ns].text = "Unique" ;
      bb[ns].arg = assVoid(501) ;
      break ;
    case 1:
      bb[ns].text = "Local repeats" ;
      bb[ns].arg = assVoid(502) ;
      break ;
    case 2:
      bb[ns].text = "Same chrom" ;
      bb[ns].arg = assVoid(503) ;
      break ;
    case 3:
      bb[ns].text = "Multi chrom" ;
      bb[ns].arg = assVoid(500) ;
      break ;
    }
  bb[ns].fg = BLACK ;
  bb[ns].bg = PALEBLUE ;
  bb[ns].next = 0 ;
  nsAmbiguous = ns ;
  ns++ ;

  bb[ns].f = htileFilterProbes;
  switch (look->tmFilter)
    {
    case 0: 
      bb[ns].text = "Any TM" ;
      bb[ns].arg = assVoid(301) ;
      break ;
    case 1:
      bb[ns].text = "TM < 72.2" ;
      bb[ns].arg = assVoid(302) ;
      break ;
    case 2:
      bb[ns].text = "TM > 75" ;
      bb[ns].arg = assVoid(300) ;
      break ;
    }
  bb[ns].fg = BLACK ;
  bb[ns].bg = PALEBLUE ;  
  bb[ns].next = 0 ;
  nsTm = ns ;
  ns++ ;

  bb[ns].f = htileFilterProbes;
  switch (look->filter99)
    {
    case 0:
      bb[ns].text = "No Filter" ;
      bb[ns].arg = assVoid(601) ;
      break ;
    case 1:
      bb[ns].text = "Medium G1" ;
      bb[ns].arg = assVoid(602) ;
      break ;
    case 2:
      bb[ns].text = "Sequence" ;
      bb[ns].arg = assVoid(603) ;
      break ;
    case 3:
      bb[ns].text = "Signal<99%" ;
      bb[ns].arg = assVoid(600) ;
      break ;
    }

  bb[ns].fg = BLACK ;
  bb[ns].bg = PALEBLUE ;
  bb[ns].next = 0 ;
  nsFilter = ns99 = ns ;
  ns++ ;

#endif

  bb[ns].f = htileFilterProbes;
  switch (look->smoothing)
    {
    case 0:
      bb[ns].text = "No smoothing..." ;
      bb[ns].arg = assVoid(71) ;
      break ;
    case 1:
      bb[ns].text = "Gaussian smoothing..." ;
      bb[ns].arg = assVoid(72) ;
      break ;
    case 2:
      bb[ns].text = "Parabolic smoothing..." ;
      bb[ns].arg = assVoid(73) ;
      break ;
    case 3:
      bb[ns].text = "Square smoothing..." ;
      bb[ns].arg = assVoid(74) ;
      break ;
    case 4:
      bb[ns].text = "Median smoothing..." ;
      bb[ns].arg = assVoid(75) ;
      break ;
    case 5:
      bb[ns].text = "Gene smoothing up..." ;
      bb[ns].arg = assVoid(76) ;
      break ;
    case 6:
      bb[ns].text = "Gene smoothing down..." ;
      bb[ns].arg = assVoid(77) ;
      break ;
    case 7:
      bb[ns].text = "Gene smoothing merged..." ;
      bb[ns].arg = assVoid(70) ;
      break ;
    }
  bb[ns].fg = BLACK ;
  bb[ns].bg = PALEBLUE ;
  bb[ns].next = 0 ;
  nsSmoother = ns ;
  ns++ ;

  bb[ns].text = 0 ;
  box = graphColouredButtons (bb, 1, offsetp, look->map->graphWidth) ;
  graphBoxFreeMenu(box + nsVpp, htileVppAction, htileVppMenu) ; 
  graphBoxFreeMenu(box + nsVmm, htileVppAction, htileVppMenu) ; 
  if (nsRatio)
    graphBoxFreeMenu(box + nsRatio, htileRatioAction, htileRatioMenu) ; 
  if (nsValue)
    graphBoxFreeMenu(box + nsValue, htileValueAction, htileValueMenu) ; 
  if (nsDot)
    graphBoxFreeMenu(box + nsDot, htileDotAction, htileDotMenu) ; 
  if (nsExonic)
    graphBoxFreeMenu(box + nsExonic, htileExonicAction, htileExonicMenu) ; 
  if (nsUnique)
    graphBoxFreeMenu(box + nsUnique, htileUniqueAction, htileUniqueMenu) ; 
  if (nsAmbiguous)
    graphBoxFreeMenu(box + nsAmbiguous, htileAmbiguousAction, htileAmbiguousMenu) ; 
  if (nsTm)
    graphBoxFreeMenu(box + nsTm, htileTmAction, htileTmMenu) ; 
  if (nsFilter)
    graphBoxFreeMenu(box + nsFilter, htileFilterAction, htileFilterMenu) ; 
  if (ns99)
    graphBoxFreeMenu(box + ns99, htile99Action, htile99Menu) ; 
  graphBoxFreeMenu(box + nsSmoother, htileSmootherAction, htileSmootherMenu) ; 
} /* htileDrawActionButtons */

/************************************************************/
typedef enum { ZDNA, RGpC, RGmC, RApG, RCpT, RApC, RTpG, RCG, RA, RT, RG, RC, NA, NT, NG, NC, HZERO } HTP ;
typedef struct typeStruct { HTP type ; const char *nam ; int min, max, col ; BOOL open ; float y ;} HTYPE ;
static  HTYPE histoTypes [] = {
    { ZDNA, "ZDNA", 0, 100, RED, FALSE, 0},
    { RA, "rate A", 15, 35, GREEN, FALSE, 0},
    { RT, "rate T", 15, 35, RED, FALSE, 0},
    { RG, "rate G", 15, 35, BLACK, FALSE, 0},
    { RC, "rate C", 15, 35, BLUE, FALSE, 0},
    { RGpC, "rate G+C", 40, 60, DARKBLUE, FALSE, 0 },
    { RGmC, "rate G-C", -10, 10, CYAN, FALSE, 0 },
    { RApG, "rate A+G", 40, 60, RED, TRUE, 0 }, 
    { RCpT, "rate C+T", 40, 60, GREEN, TRUE, 0 }, 
    { RApC, "rate A+C", 40, 60, VIOLET, FALSE, 0 },
    { RTpG, "rate T+G", 40, 60, LIGHTMAGENTA, FALSE, 0 },
    { RCG, "rate CpG", 0, 10, CERISE, FALSE, 0  },
    { HZERO, 0, 0, 0, BLACK, FALSE, 0 }} ;

/************************************************************/

static void htileGenomeAction (void *v)
{
  int ns = assInt (v) ;
  HTYPE *tt ;
  Htile look = currentHtile("htileFilterProbes") ;

  switch (ns)
    {
    case 2000:
      look0->showGenes = look->showGenes = FALSE ;
      break ;
    case 2001:
      look0->showGenes = look->showGenes = TRUE ;
      break ;
    case 4000:
      look0->showCaptureProbes = look->showCaptureProbes = FALSE ;
      break ;
    case 4001:
      look0->showCaptureProbes = look->showCaptureProbes = TRUE ;
      break ;
    case 3000:
      look0->showGeneSignal = look->showGeneSignal = FALSE ;
      break ;
    case 3001:
      look0->showGeneSignal = look->showGeneSignal = TRUE ;
      break ;
    case 1000:
      look0->showGenome = look->showGenome = FALSE ;
      break ;
    case 1001:
      look0->showGenome = look->showGenome = TRUE ;
      break ;
    default:
     for (tt = histoTypes ; tt->nam ; tt++)
       if (tt->type == ns)
	 tt->open = ! tt->open ;
    }
  look->map->mapDraw () ;
} /*  htileGenomeAction */

static void htileDrawGenomeButtons (Htile look, float *offsetp)
{
  COLOUROPT bb[NSMAX] ;
  int ns = 0 ;
  HTYPE *tt ;

  bb[ns].f = htileGenomeAction ;
  bb[ns].text = look->showGenes ? "close gene tracks" : "open gene tracks" ;
  bb[ns].fg = BLACK ;
  bb[ns].bg = PALEGREEN ;
  bb[ns].arg = look->showGenes ? assVoid(2000) : assVoid(2001) ;
  bb[ns].next = 0 ;
  ns++ ;

  bb[ns].f = htileGenomeAction ;
  bb[ns].text = look->showCaptureProbes ? "close capture tracks" : "open capture tracks" ;
  bb[ns].fg = BLACK ;
  bb[ns].bg = PALEGREEN ;
  bb[ns].arg = look->showCaptureProbes ? assVoid(4000) : assVoid(4001) ;
  bb[ns].next = 0 ;
  ns++ ;

  bb[ns].f = htileGenomeAction ;
  bb[ns].text = look->showGenome ? "close" : "open genome tracks" ;
  bb[ns].fg = BLACK ;
  bb[ns].bg = PALEGREEN ;
  bb[ns].arg = look->showGenome ? assVoid(1000) : assVoid(1001) ;
  bb[ns].next = 0 ;
  ns++ ;

  if (look->showGenome)
    for (tt = histoTypes ; tt->nam ; tt++)
      {
	bb[ns].f = htileGenomeAction ;
	bb[ns].text = tt->nam ;
	bb[ns].fg = BLACK ;
	bb[ns].bg = look->showGenome && tt->open ? tt->col : WHITE ;
	bb[ns].arg = assVoid(tt->type) ;
	bb[ns].next = 0 ;
	ns++ ;
      }
	  
  bb[ns].text = 0 ;
  graphColouredButtons (bb, 2, offsetp, look->map->graphWidth) ;
  return ;
} /* htileDrawGenomeButtons */

/************************************************************/

static void htileDrawScale (Htile look, float offset)
{
  float u10, u20, u1, u2, u3, y1, y2 ;
  int sBox, i, i0, i1 = 0, sc ;
  int sc1, sc3, sc4 ;
  int v1 = TGRAPH2MAP (look->map, 2) ;
  int v2 = TGRAPH2MAP (look->map, look->map->graphWidth - 2) ;

  if (v1 < look->map->min) v1 = look->map->min ;
  if (v2 > look->map->max) v2 = look->map->max ;

  v1 += look->map->a1 - 1 ;
  v2 += look->map->a1 - 1 ;
  /* we need to set sc,sc1 twice becuase if v2 - v1 = 8, sc = only 5 the first time */
  sc = utMainRoundPart (v2 - v1) ;
  sc1 =  utMainPart (sc/20) ;
  if (sc1 < 1) sc1 = 1 ;
  v1 -= v1 % sc1 ;
  v2 = v2 - (v2 % sc1) + sc1 ;
  sc = v2 - v1 ;
  sc1 =  utMainPart (sc/20) ;
  if (sc1 < 1) sc1 = 1 ;
  u10 = u1 = TMAP2GRAPH (look->map, v1 - look->map->a1 + 1) ;
  u20 = u2 = TMAP2GRAPH (look->map, v2 - look->map->a1 + 1) ;

  sBox = graphBoxStart () ;

  if (v1 - look->map->a1 > look->map->min) u10 = 0 ;
  if (v2 - look->map->a1 < look->map->max) u20 =  look->map->graphWidth ;
  y2 = offset ; /* look->map->graphHeight - .07 ; */
  graphLine (u10, y2, u20, y2) ;
  y1 = y2 - .4 ;
  graphLine (u1, y1, u1, y2) ; graphLine (u2, y1, u2, y2) ;
  y1 = y2 - .2 ;

  for (i0 = i = 0 ; i * sc1 <= sc ; i++)
    { 
      sc3 = v1 + i*sc1 ;
      u3 = TMAP2GRAPH (look->map, sc3 - look->map->a1 + 1) ;
      if (u3 < 1 || u3 > look->map->graphWidth - 2) continue ;
      if (
	  (
	   (sc > 20 && ( (20*sc1!=sc && (sc3/sc1) % 10 == 0) || (20*sc1==sc && (sc3/sc1) % 5 == 0))) ||
	   (sc >= 10 && sc <= 20 && i %5 == 4) ||
	   (sc >= 5 && sc < 10 && i%2 == 1) ||
	   (sc < 5)
	   )
	  )
	{
	  graphLine (u3, y2 - .5, u3, y2) ; 
	  if (i > 0)
	    {
	      if (i1 == 0) i1 = i ;
	      sc4 = 10*sc1 + look->map->a1 ; ;
	      if (sc4 >= 1000000 && sc3 % 1000000 == 0)
		graphText (messprintf ("%dM", sc3/1000000), u3 - .3 , y1 - 1.4) ;
	      else if (sc4 >= 1000 && sc3 % 1000 == 0)
		graphText (messprintf ("%dk", sc3/1000), u3 - .3 , y1 - 1.4) ;
	      else
		graphText (messprintf ("%d", sc3), u3 - .3 , y1 - 1.4) ;
	    }
	  else if (sc3 == 0)
	    graphText (messprintf ("%d", sc3), u3 - .3 , y1 - 1.4) ;
	  i0 = i ;
	}
      else
	graphLine (u3, y1, u3, y2) ;
    }

  sc3 = v1 + i0*sc1 ;
  u3 = TMAP2GRAPH (look->map, sc3 - look->map->a1 + 1) ;
  sc4 = 10*sc1 + look->map->a1 ; ;
  if (sc4 >= 1000000 && sc3%1000000 == 0)
    graphText (messprintf ("%dM", sc3/1000000), u3 - .3 , y1 - 1.4) ;
  else if (sc4 >= 1000 && sc3%1000 == 0)
    graphText (messprintf ("%dk", sc3/1000), u3 - .3 , y1 - 1.4) ;
  else
    graphText (messprintf ("%d", sc3), u3 - .3 , y1 - 1.4) ;

  graphText (" ",  1, y2 - .6) ;
  graphBoxEnd () ;
  graphBoxDraw (sBox, BLACK, WHITE) ;
  /* draw some dna */
  y2++ ;
} /* htileDrawScale */

/************************************************************/

static void htileDrawDna (Htile look, float offset)
{
  AC_HANDLE h = ac_new_handle () ;

  float y2 ;
  int jj[32], ccol[32], width[32] ;

  int  x0 = look->map->graphWidth/2.0 ;
  int i, j, a1, a2, ai, ir, ln, a0 = TGRAPH2MAP (look->map, x0) ;
  AC_OBJ Cosmid = 0, Region = htile_key_obj (look, look->locus, h) ;
  AC_KEYSET cosmids = ac_objquery_keyset (Region, "{CLASS mRNA} SETOR {CLASS Sequence} SETOR {>sequence;>Parts} $| {>Sequence};CLASS mRNA OR Is_cosmid", h) ;
  char *qq = "select c,m,a1,a2 from c in @, m in c->intmap, a1 in m[1], a2 in a1[1]" ;
  char *dna ;
  char buf[2] ;
  float x, dx = look->map->mag > 1 ? look->map->mag : 1.0 ;
  AC_TABLE cosmidTable = ac_bql_table (look->db, qq, cosmids, "+2+3+4+1", 0, h) ;
  
  buf[1] = 0 ;
  
  y2 = look->map->dnaOffset = offset ;
  for (ir = 0 ; cosmidTable && ir < cosmidTable->rows ; ir++)
    {
      if (ac_table_key (cosmidTable, ir, 2, 0))
	{
	  a1 = ac_table_int (cosmidTable, ir, 2, 0) - look->map->a1 + 1 ;
	  a2 = ac_table_int (cosmidTable, ir, 3, 0) - look->map->a1 + 1 ;
	}
      else
	{
	  a1 = look->map->a1 ; 
	  a2 = look->map->a2 ;
	}
      a1 = a1 - look->map->a1 + 1 ; 
      a2 = a2 - look->map->a1 + 1 ; 
      if (a1 > a0 || a2 < a0) /* area to be shown */
	continue ;
      Cosmid = ac_table_obj (cosmidTable, ir, 0, h) ;
      dna = ac_obj_dna (Cosmid, h) ;
      ln = strlen (dna) ;
      j = 0 ;
      memset (jj, 0, sizeof(jj)) ;
      if (a0 > 100 && a0 + 100 < ln)
	{
	  char c = *(dna + a0 + 100) ;
	  *(dna + a0 + 100) = 0 ;
	  width[j] = 9 ; ccol[j]=RED; jj[j++] = pickMatch (dna + a0 - 100, "*ttc???gaa*") ;
	  width[j] = 9 ; ccol[j]=ORANGE; jj[j++] = pickMatch (dna + a0 - 100, "*ttc???taa*") ;
	  width[j] = 9 ; ccol[j]=ORANGE; jj[j++] = pickMatch (dna + a0 - 100, "*tta???gaa*") ;
	  width[j] = 5 ; ccol[j]=BLUE; jj[j++] = pickMatch (dna + a0 - 100, "*gtttc*") ;
	  width[j] = 5 ; ccol[j]=BLUE; jj[j++] = pickMatch (dna + a0 - 100, "*gaaac*") ;
	  width[j] = 6 ; ccol[j]=GREEN7; jj[j++] = pickMatch (dna + a0 - 100, "*aggaag*") ;
	  width[j] = 6 ; ccol[j]=GREEN7; jj[j++] = pickMatch (dna + a0 - 100, "*cttcct*") ;
	  *(dna + a0 + 100) = c ;
	}
      for (i = a0 - 100 ; i < a0 + 100 ; i++)
	if (i >= a1 && i < a2)
	  {
	    x = look->map->graphWidth/2.0 + (i - a0) * dx - .5 ;
	    ai = i - a1 ;
	    if (ai >= 3 && ai < ln)
	      {
		buf[0] = dna[ai] ;
		if (x > 20 && x < look->map->graphWidth - 2)
		  {
		    int col = 0 ;
		    for (j = 0 ; j < 32 ; j++)
		      if (jj[j] && i >= jj[j] + a0 - 100 && i <= jj[j] + a0 - 101 + width[j])
			{
			  col = ccol[j] ;
			  break ;
			}
		    if (col) graphColor (col) ;
		    graphText (buf, x, y2-.5) ;
		    if (col) 	  graphColor (BLACK) ;
		  }
	      }
	  }
    }
  ac_free (h) ;
} /* htileDrawDna */

/************************************************************/

static void htileDrawGenes (Htile look, float offset, BOOL isUp)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  AC_OBJ Locus = ac_get_obj (look->db, className (look->locus), name(look->locus), h) ;
  AC_ITER tiles = 0 ;
  AC_OBJ Gene = 0, Tg = 0, tile = 0 ;
  AC_TABLE mm, tgs, spls, ss ;
  SEG *gSeg ;
  KEY gene, tg ;
  const char *ccp ;
  float u1, u2, v1, v2, oldv1 = 0, oldv2, y1, y2, signal ;
  int pass, ir, jr,  g1, g2, a1, a2, b1, b2, x1, x2 , geneBox, col = 0, bcol ;
  BitSet bb = bitSetCreate (1000000, h) ;
  vTXT buf = vtxtHandleCreate (h) ;
  KEY _Source_exons = str2tag("Source_exons") ;

  tiles = ac_objquery_iter (Locus, " {CLASS Sequence } SETOR {>Sequence} SETOR {>In_junction} ; {IS * } SETOR {>Parts} ; IntMap && Transcribed_gene", h) ; /* miegTOTO */
  if (isUp) 
    {
      if (arrayExists (look->map->genes))
	arrayDestroy (look->map->genes) ;
      look->map->genes = arrayHandleCreate (300, SEG, look->h) ;
    }
  mm = ac_tag_table (Locus, "IntMap", h) ;
  g1 = ac_table_int (mm, 0, 1, 1) ;
  g2 = ac_table_int (mm, 0, 2, 1) ;
  while (ac_free (tile), tile = ac_iter_obj (tiles))
    {
      mm = ac_tag_table (tile, "IntMap", h) ;
      a1 = ac_table_int (mm, 0, 1, 1) ;
      a2 = ac_table_int (mm, 0, 2, 1) ;
      if (a2 < g1 || a1 > g2)
	continue ;
      	
      for (pass = 0 ; pass < 5 ; pass++) /* first the genes with a geneid */
	{
	  switch (pass)
	    {
	    case 0:
	    case 1:
	      tgs = ac_tag_table (tile, "Transcribed_gene", h) ;
	      y1 = offset + (isUp ? -.6 : 0) ; y2 = y1 + .6 ;   /* Defines thickness of the gene. Was set to .6 originally*/
	      break ;
	    case 2: 
	      {
		AC_KEYSET kpg = ac_objquery_keyset (tile, ">Genes ; {CLASS Gene }  SETOR {>Genefinder;>Source;>genes;}; NOT Transcribed_gene",h) ;
		tgs = ac_keyset_table (kpg, 0, -1, 1, h) ;
	      }
	      y1 = offset + (isUp ? -2 : 2) ; y2 = y1 + .6 ;   /* Defines thickness of the gene. Was set to .6 originally*/
	      break ;
	    case 3:
	      {
		AC_KEYSET kpg = ac_objquery_keyset (tile, ">Genes;{>Transcribed_gene;>Matching_genefinder_gene} SETOR{>Genefinder};>Source;>subsequence;Is_predicted_gene;NOT Method=\"Sasha*\" AND NOT Method=\"Vega*\"",h) ;
		tgs = ac_keyset_table (kpg, 0, -1, 1, h) ;
	      }
	      y1 = offset + (isUp ? -5 : 5) ; y2 = y1 + .6 ;   /* Defines thickness of the gene. Was set to .6 originally*/
	      break ;
	    case 4: 
	      {
		AC_KEYSET kpg = ac_objquery_keyset (tile, ">Subsequence ; Is_predicted_gene",h) ;
		tgs = ac_keyset_table (kpg, 0, -1, 1, h) ;
	      }
	      y1 = offset + (isUp ? -2 : 2) ; y2 = y1 + .6 ;   /* Defines thickness of the gene. Was set to .6 originally*/
	      break ;
	    }
	  for (ir = 0 ; tgs && ir < tgs->rows ; ir++)
	    {
	      tg = ac_table_key (tgs, ir, 0, 0) ;
	      ac_free (h1) ;
	      h1 = ac_new_handle () ;
	      
	      switch (pass)
		{
		case 0:
		  Tg = ac_table_obj (tgs, ir, 0, h1) ;
		  Gene = ac_tag_obj (Tg, "Gene", h) ;
		  gene = Gene ? ac_obj_key (Gene) : 0 ;
		  if (!keyFindTag (gene, _GeneId))
		    continue ;
		  col = BLUE ;
		  break ;
		case 1:
		  Tg = ac_table_obj (tgs, ir, 0, h1) ;
		  Gene = ac_tag_obj (Tg, "Gene", h) ;
		  gene = Gene ? ac_obj_key (Gene) : 0 ;
		  if (keyFindTag (gene, _GeneId))
		    continue ;
		  if (ac_has_tag (Gene, "cloud"))
		    { col = BLUE4 ;  y1 = offset + (isUp ? -2.7 : 2.7) ; y2 = y1 + .4 ; }
		  else
		    col = BLUE7 ;
		  break ;
		case 2:
		  gene = tg ; tg = 0 ; Tg = 0 ;
		  Gene = Tg = ac_table_obj (tgs, ir, 0, h1) ;
		  if (ac_has_tag (Gene, "Non_protein_coding"))
		    col = ORANGE ;
		  else
		    col = GREEN ;
		  break ;
		case 3:
		case 4:
		  if (! keyFindTag (tg, _Source_exons))
		    continue ;
		  col = GRAY ;
		  gene = tg ; 
		  Gene = Tg = ac_table_obj (tgs, ir, 0, h1) ;
		  break ;
		}
	      mm = ac_tag_table (Tg, "IntMap", h1) ;
	      a1 = ac_table_int (mm, 0, 1, 1) ;
	      a2 = ac_table_int (mm, 0, 2, 1) ;
	      if ((a1 >= a2 && isUp) || (a1 <= a2 && !isUp))
		continue ;
	      x1 = a1 - look->map->a1 + 1 ;
	      x2 = a2 - look->map->a1 + 1 ;
	      u1 = TMAP2GRAPH (look->map, x1) ;
	      u2 = TMAP2GRAPH (look->map, x2) ;
	      if (u1 > u2) { float u0 = u1 ; u1 = u2 ; u2 = u0 ; }
	      if (u2 < 0 || u1 > look->map->graphWidth)
		continue ;
	      if (u1 < 0)
		u1 = 0 ;
	      if (u2 >  look->map->graphWidth)
		u2 =  look->map->graphWidth ;
	      if (ac_has_tag (Gene, "MAQC"))
		{
		  gSeg = arrayp (look->map->genes, arrayMax(look->map->genes), SEG) ;
		  gSeg->key = gene ;
		  gSeg->x1 = x1 ;
		  gSeg->x2 = x2 ;
		  gSeg->a1 = x1 < x2 ? x1 : x2 ;
		  gSeg->a2 = x1 < x2 ? x2 : x1 ;
		  signal = -128 ;
		  ss = ac_tag_table (Gene, "Specific_exon_signal", 0) ;
		  if (ss)
		    {   /* average of pA+ signal in exon specific probes */
		      int ir1, nns = 0, n1, nng=0;
		      float s1,signalg = 0 ;
		      const char *ccp ;
		      signal = 0 ;
		      for (ir1 = 0 ; ir1 < ss->rows ; ir1++)
			{
			  s1 = ac_table_float (ss, ir1, 1, -128) ;
			  n1 = ac_table_int (ss, ir1, 2, 0) ;
			  ccp = ac_table_printable (ss, ir1, 0, "toto") ;
			  if (n1 && strstr (ccp, "pA+"))
			    {
			      signal += s1*n1 ; nns += n1 ;
			    }
			  if (n1 && strstr (ccp, "G1"))
			    {
			      signalg += s1*n1 ; nng += n1 ;
			    }
			}
		      ac_free (ss) ;
		      if (nng)  signalg /= nng ;
		      if (nns)  signal /= nns ;
		      if (nns && nng)  signal = signal - signalg + 10 ;
		      else   signal = -128 ; 
		    }
		  if (signal < 1) col = GRAY ; /* not tested */
		  else if (signal < 9.62) col = BLUE1 ;
		  else if (signal < 9.755) col = BLUE2 ;
		  else if (signal < 10.0) col = BLUE4 ;
		  else if (signal < 10.2) col = BLUE5 ;
		  else if (signal < 10.6) col = BLUE6 ;
		  else if (signal < 99) col =  BLUE8 ;
		  gSeg->col = col ;
		}
	      array (look->map->boxIndex, geneBox = graphBoxStart (), int) = gene ;
	      if (u2 > 0 && u1 < look->map->graphWidth)
		{
		  switch (pass)
		    {
		    case 2:
		      v1 = TMAP2GRAPH (look->map, x1 - .5) ;
		      v2 = TMAP2GRAPH (look->map, x2 + .5) ;
		      if (v1 < 0) v1 = 0 ;
		      if (v2 < 0) v2 = 0 ;
		      if (v1 >  look->map->graphWidth)
			v1 =  look->map->graphWidth ;
		      if (v2 >  look->map->graphWidth)
			v2 =  look->map->graphWidth ;
		      graphRectangle (v1, y1, v2, y2) ;
		      break ;
		    case 3:
		    case 4:
		      oldv2 = 0 ;
		      
		      spls = ac_tag_table (Tg, "Source_exons", h1) ;
		      for (jr = 0 ; spls && jr < spls->rows ; jr++)
			{
			  b1 = ac_table_int (spls, jr, 0, 0) ;
			  b2 = ac_table_int (spls, jr, 1, 0) ;
			  if (x1 < x2)
			    {
			      v1 = TMAP2GRAPH (look->map, x1 + b1 - 1 - .5) ;
			      v2 = TMAP2GRAPH (look->map, x1 + b2 - 1 + .5) ;
			    }
			  else
			    {
			      v1 = TMAP2GRAPH (look->map, x1 - b2 + 1 - .5) ;
			      v2 = TMAP2GRAPH (look->map, x1 - b1 + 1 + .5) ;
			    }
			  if (v1 < 0) v1 = 0 ;
			  if (v2 < 0) v2 = 0 ;
			  if (v1 >  look->map->graphWidth)
			    v1 =  look->map->graphWidth ;
			  if (v2 >  look->map->graphWidth)
			    v2 =  look->map->graphWidth ;
			  if (jr > 0)
			    {  
			      if (x1 < x2 && v1 > oldv2 )
				graphLine (oldv2, (y1 + y2)/2, v1, (y1+y2)/2) ;			
			      if (x1 > x2 && v2 < oldv1 )
				graphLine (v2, (y1 + y2)/2, oldv1, (y1+y2)/2) ;			
			    }
			  if (v1 < v2)
			    graphRectangle (v1, y1, v2, y2) ;
			  oldv1 = v1 ;
			  oldv2 = v2 ;
			}
		      break ;
		    case 0:
		    case 1:
		      spls = ac_tag_table (Tg, "Splicing", h1) ;
		      for (jr = 0 ; spls && jr < spls->rows ; jr++)
			{
			  ccp = ac_table_printable (spls, jr, 2, "") ;
			  b1 = ac_table_int (spls, jr, 0, 0) ;
			  b2 = ac_table_int (spls, jr, 1, 0) ;
			  if (x1 < x2)
			    {
			      v1 = TMAP2GRAPH (look->map, x1 + b1 - 1 - .5) ;
			      v2 = TMAP2GRAPH (look->map, x1 + b2 - 1 + .5) ;
			    }
			  else
			    {
			      v1 = TMAP2GRAPH (look->map, x1 - b2 + 1 - .5) ;
			      v2 = TMAP2GRAPH (look->map, x1 - b1 + 1 + .5) ;
			    }
			  if (v1 < 0) v1 = 0 ;
			  if (v2 < 0) v2 = 0 ;
			  if (v1 >  look->map->graphWidth)
			    v1 =  look->map->graphWidth ;
			  if (v2 >  look->map->graphWidth)
			    v2 =  look->map->graphWidth ;
			  
			  if (v1 < v2)
			    {
			      if (strstr (ccp, "xon"))
				graphRectangle (v1, y1, v2, y2) ;
			      else 
				graphLine (v1, (y1 + y2)/2, v2, (y1+y2)/2) ;
			    }
			}
		      break ;
		    }
		}
	      if (isGifDisplay)
		{
		  KEY pg = 0 ;
		  /* mouse over for the web */
		  vtxtClear (buf) ;
		  vtxtPrintf (buf, "%s %s"
			      , ac_key_name(gene)
			      , ac_tag_printable (Gene, "Title", "wormbase WS190 model")
			      ) ;
		  if (pass == 3)
		    pg = keyGetKey (gene, str2tag("Predicted_mRNA")) ;
		  if (pg)
		   graphBubbleInfo (geneBox, pg, "mrna",  vtxtPtr (buf)) ;
		  else
		   graphBubbleInfo (geneBox, gene, "Gene", vtxtPtr (buf)) ;
		}
	      
	      if (gene && col != GRAY)
		{
		  float uu=0, dx ;
		  int dn = strlen (name(gene)) ;
		  char buf[256] ;
		  int ib ;
		  BOOL ok = TRUE ;
		  
		  if (u2 > look->map->graphWidth)
		    u2 = look->map->graphWidth ;
		  if (u1 < 0)
		    u1 = 0 ;
		  dx = u2 - u1 ;
		  if (dx > 0)
		    {
		      memset (buf, 0, sizeof(buf)) ;
		      if (dx < 8) dx = 8 ;
		      if (dx > dn) dx = dn ;
		      if (dx > 255) dx = 255 ;	
		      while (dx > 1)
			{
			  ok = TRUE ;
			  uu = (u1 + u2)/2 - dx/2 ;
			  for (ib = uu ; ib < uu + dx ; ib++)
			    if (ib >= 0 && bit(bb, ib))
			      ok = FALSE ;
			  if (ok) break ;
			  dx-- ;
			}
		      if (dx >= 1 && ok)
			{
			  strncpy (buf, name(gene), dx) ;
			  {
			    float oldh = graphTextHeight (1.5) ;
			    graphText (buf, uu+.1, (isUp ? y1 - 1.0 : y2+.4)) ;
			    graphTextHeight (oldh) ;
			  }
			  
			  for (ib = uu ; ib < uu + dx ; ib++)
			    if(ib >= 0) bitSet (bb, ib) ;
			}
		    }
		}
	      if (look0->activeGenes && keySetFind (look0->activeGenes, gene, 0))
		bcol = YELLOW ;
	      else
		bcol = TRANSPARENT ;
	      
	      graphBoxEnd () ;
	      graphBoxDraw (geneBox, col, bcol) ;
	    }
	}
    }
  graphColor (BLACK) ;
  ac_free (h1) ;
  ac_free (h) ;
  arraySort (look->map->genes, htileSegOrder) ;
} /* htileDrawGenes */

/************************************************************/

static void htileDrawCaptureProbes (Htile look, float offset, BOOL isUp)
{
#ifdef JUNK
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  AC_OBJ Locus = ac_get_obj (look->db, className (look->locus), name(look->locus), h) ;
  AC_ITER tiles = 0 ;
  AC_OBJ Gene = 0, Tg = 0, tile = 0 ;
  AC_TABLE mm, tgs, spls, ss ;
  SEG *gSeg ;
  KEY gene, tg ;
  const char *ccp ;
  float u1, u2, v1, v2, oldv1 = 0, oldv2, y1, y2, signal ;
  int pass, ir, jr,  g1, g2, a1, a2, b1, b2, x1, x2 , geneBox, col = 0, bcol ;
  BitSet bb = bitSetCreate (1000000, h) ;
  vTXT buf = vtxtHandleCreate (h) ;
  KEY _Source_exons = str2tag("Source_exons") ;

  tiles = ac_objquery_iter (Locus, " {CLASS Sequence } SETOR {>Sequence} SETOR {>In_junction} ; {IS * } SETOR {>Parts} ; IntMap && Transcribed_gene", h) ; /* miegTOTO */
  if (isUp) 
    {
      if (arrayExists (look->map->genes))
	arrayDestroy (look->map->genes) ;
      look->map->genes = arrayHandleCreate (300, SEG, look->h) ;
    }
  mm = ac_tag_table (Locus, "IntMap", h) ;
  g1 = ac_table_int (mm, 0, 1, 1) ;
  g2 = ac_table_int (mm, 0, 2, 1) ;
  while (ac_free (tile), tile = ac_iter_obj (tiles))
    {
      mm = ac_tag_table (tile, "IntMap", h) ;
      a1 = ac_table_int (mm, 0, 1, 1) ;
      a2 = ac_table_int (mm, 0, 2, 1) ;
      if (a2 < g1 || a1 > g2)
	continue ;
      	
      for (pass = 0 ; pass < 1 ; pass++)
	{
	  switch (pass)
	    {
	    case 0:
	    case 1:
	      tgs = ac_tag_table (tile, "Transcribed_gene", h) ;
	      y1 = offset + (isUp ? -.6 : 0) ; y2 = y1 + .6 ;   /* Defines thickness of the gene. Was set to .6 originally*/
	      break ;

	  for (ir = 0 ; tgs && ir < tgs->rows ; ir++)
	    {
	      tg = ac_table_key (tgs, ir, 0, 0) ;
	      ac_free (h1) ;
	      h1 = ac_new_handle () ;
	      
	      switch (pass)
		{
		case 0:
		  Tg = ac_table_obj (tgs, ir, 0, h1) ;
		  Gene = ac_tag_obj (Tg, "Gene", h) ;
		  gene = Gene ? ac_obj_key (Gene) : 0 ;
		  if (!keyFindTag (gene, _GeneId))
		    continue ;
		  col = BLUE ;
		  break ;
		}

	      mm = ac_tag_table (Tg, "IntMap", h1) ;
	      a1 = ac_table_int (mm, 0, 1, 1) ;
	      a2 = ac_table_int (mm, 0, 2, 1) ;
	      if ((a1 >= a2 && isUp) || (a1 <= a2 && !isUp))
		continue ;
	      x1 = a1 - look->map->a1 + 1 ;
	      x2 = a2 - look->map->a1 + 1 ;
	      u1 = TMAP2GRAPH (look->map, x1) ;
	      u2 = TMAP2GRAPH (look->map, x2) ;
	      if (u1 > u2) { float u0 = u1 ; u1 = u2 ; u2 = u0 ; }
	      if (u2 < 0 || u1 > look->map->graphWidth)
		continue ;
	      if (u1 < 0)
		u1 = 0 ;
	      if (u2 >  look->map->graphWidth)
		u2 =  look->map->graphWidth ;
	      if (ac_has_tag (Gene, "MAQC"))
		{
		  gSeg = arrayp (look->map->genes, arrayMax(look->map->genes), SEG) ;
		  gSeg->key = gene ;
		  gSeg->x1 = x1 ;
		  gSeg->x2 = x2 ;
		  gSeg->a1 = x1 < x2 ? x1 : x2 ;
		  gSeg->a2 = x1 < x2 ? x2 : x1 ;
		  gSeg->col = RED ;
		}
	      array (look->map->boxIndex, geneBox = graphBoxStart (), int) = gene ;
	      if (u2 > 0 && u1 < look->map->graphWidth)
		{
		  switch (pass)
		    {
		    case 0:
		      v1 = TMAP2GRAPH (look->map, x1 - .5) ;
		      v2 = TMAP2GRAPH (look->map, x2 + .5) ;
		      if (v1 < 0) v1 = 0 ;
		      if (v2 < 0) v2 = 0 ;
		      if (v1 >  look->map->graphWidth)
			v1 =  look->map->graphWidth ;
		      if (v2 >  look->map->graphWidth)
			v2 =  look->map->graphWidth ;
		      graphRectangle (v1, y1, v2, y2) ;
		      break ;
		    }
		}
	      if (isGifDisplay)
		{
		  KEY pg = 0 ;
		  /* mouse over for the web */
		  vtxtClear (buf) ;
		  vtxtPrintf (buf, "%s %s"
			      , ac_key_name(gene)
			      , ac_tag_printable (Gene, "Title", "Capture probe")
			      ) ;
		   graphBubbleInfo (geneBox, gene, "Probe", vtxtPtr (buf)) ;
		}
	      
	      if (gene && col != GRAY)
		{
		  float uu=0, dx ;
		  int dn = strlen (name(gene)) ;
		  char buf[256] ;
		  int ib ;
		  BOOL ok = TRUE ;
		  
		  if (u2 > look->map->graphWidth)
		    u2 = look->map->graphWidth ;
		  if (u1 < 0)
		    u1 = 0 ;
		  dx = u2 - u1 ;
		  if (dx > 0)
		    {
		      memset (buf, 0, sizeof(buf)) ;
		      if (dx < 8) dx = 8 ;
		      if (dx > dn) dx = dn ;
		      if (dx > 255) dx = 255 ;	
		      while (dx > 1)
			{
			  ok = TRUE ;
			  uu = (u1 + u2)/2 - dx/2 ;
			  for (ib = uu ; ib < uu + dx ; ib++)
			    if (ib >= 0 && bit(bb, ib))
			      ok = FALSE ;
			  if (ok) break ;
			  dx-- ;
			}
		      if (dx >= 1 && ok)
			{
			  strncpy (buf, name(gene), dx) ;
			  {
			    float oldh = graphTextHeight (1.5) ;
			    graphText (buf, uu+.1, (isUp ? y1 - 1.0 : y2+.4)) ;
			    graphTextHeight (oldh) ;
			  }
			  
			  for (ib = uu ; ib < uu + dx ; ib++)
			    if(ib >= 0) bitSet (bb, ib) ;
			}
		    }
		}
	      bcol = TRANSPARENT ;
	      
	      graphBoxEnd () ;
	      graphBoxDraw (geneBox, col, bcol) ;
	    }
	}
    }
  graphColor (BLACK) ;
  ac_free (h1) ;
  ac_free (h) ;
  arraySort (look->map->genes, htileSegOrder) ;
#endif
} /* htileDrawCaptureProbes */

/************************************************************/

static BOOL htilePBSCol (KEY type, int *colp, int *dyp)
{
  Tmap *map = tmapCurrent ("tmapResize") ;
  Htile look = (Htile) map->look ;
  int n, col = GRAY ;
  static int n0 = 0 ;
  static DICT *dict = 0 ;
  static Array cols = 0 ;
  PNX *pnx ;
  char *cp, buf[128] ;

  if (!type) return col ;
  if (!dict)
    {
      dict = dictCreate (32) ;
      cols = arrayCreate (32, int) ;
      for (pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p ; pnx++)
	{
	  strncpy (buf, pnx->p, 127) ;
	  for (cp = buf ; *cp; cp++)
	    switch ((int)*cp)
	      {
	      case '/': *cp='_'; break ;
	      case '+': *cp='p'; break ;
	      case '-': *cp='m'; break ;
	      }
	  dictAdd (dict, buf, &n) ;
	  array (cols, n, int) = pnx->col ;
	  n0 = n ;
	}
    }
  if (dictFind (dict, name(type), &n))
    {
      *dyp = n > n0 ? n - n0 : n ;
      *colp = arr (cols, n, int) ;
      return TRUE ;
    }
  return FALSE ;
} /* htilePBSCol */

/************************************************************/

static void htileDrawPBS (Htile look, float offset, BOOL isUp)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  AC_OBJ Locus = ac_get_obj (look->db, className (look->locus), name(look->locus), h) ;
  AC_ITER tiles = 0 ;
  AC_OBJ Pbs = 0, tile = 0 ;
  AC_TABLE mm, pbsT, types ;
  KEY pbs ;
  PNX *pnx ;
  float u1, u2 ;
  int ir, g1, g2, a1, a2, x1, x2 , pbsBox, col, dy, width, ns ;

  /*
  float y1 = offset ; 
  float y2 = y1 + .6 ;   // Defines thickness of the gene. Was set to .6 originally
  */ 
  
  tiles = ac_objquery_iter (Locus, " {CLASS Sequence } SETOR {>Sequence} ; {IS * } SETOR {>Parts} ", h) ; /* miegTOTO */
  mm = ac_tag_table (Locus, "IntMap", h) ;
  g1 = ac_table_int (mm, 0, 1, 1) ;
  g2 = ac_table_int (mm, 0, 2, 1) ;
  while (ac_free (tile), tile = ac_iter_obj (tiles))
    {
      mm = ac_tag_table (tile, "IntMap", h) ;
      a1 = ac_table_int (mm, 0, 1, 1) ;
      a2 = ac_table_int (mm, 0, 2, 1) ;
      if (a2 < g1 || a1 > g2)
	continue ;

      pbsT = ac_tag_table (tile, "PBS", h) ;
      for (ir = 0 ; pbsT && ir < pbsT->rows ; ir++)
	{
	  pbs = ac_table_key (pbsT, ir, 0, 0) ;
	  
	  ac_free (h1) ;
	  h1 = ac_new_handle () ;
	  Pbs = ac_table_obj (pbsT, ir, 0, h1) ;
	  mm = ac_tag_table (Pbs, "IntMap", h1) ;
	  a1 = ac_table_int (mm, 0, 1, 1) ;
	  a2 = ac_table_int (mm, 0, 2, 1) ;
	  if ((a1 > a2 && isUp) || (a1 < a2 && !isUp))
	    continue ;
	  width = ac_tag_int (Pbs, "Width", 0);
	  if (width)
	    {
	      int a0 = (a1 + a2)/2 ;
	      a1 = a0 - width/2 ;
	      a2 = a0 + width/2 ;
	    }
	  x1 = a1 - look->map->a1 ;
	  x2 = a2 - look->map->a1 ;
	  u1 = TMAP2GRAPH (look->map, x1 - .5) ;
	  u2 = TMAP2GRAPH (look->map, x2 + .5) ;
	  if (u1 > u2) { float u0 = u1 ; u1 = u2 ; u2 = u0 ; }
	  if (u2 < 0 || u1 > look->map->graphWidth)
	    continue ;
	  if (ac_has_tag (Pbs, "Peak"))
	    {
	       if (ac_has_tag (Pbs, "IgG"))
		 col = GRAY ;
	       else
		 col = RED ;
	       dy = 0 ; 
	    }
	  else
	    {
	      col = GRAY ; dy = 1 ;
	      types = ac_tag_table (Pbs, "Type", h1) ;
	      for (ns = 0, pnx = arrp(look->solexaAll, 0, PNX) ;pnx->p &&  ns < SOLEXAMAX ; pnx++, ns++)
		for (pnx = arrp(look->solexaAll, 0, PNX) ; pnx->p ; pnx++)
		  if (!strcmp (pnx->p, ac_table_printable (types, 0, 1, "toto")) &&
		      look->isSolexaClosed[ns])
		    { col = -1 ; break ; }
	      if (col == -1) /* do not draw closed lanes */
		continue ; 
	      htilePBSCol (ac_table_key (types, 0, 1, 0), &col, &dy) ;
	    }
	  array (look->map->boxIndex, pbsBox = graphBoxStart (), int) = pbs ;
	  graphRectangle (u1, offset +(dy-1)/3.0, u2, offset +dy/3.0) ;
	  graphBoxEnd () ;
	  graphBoxDraw (pbsBox, BLACK, col) ;
	}
    }
  graphColor (BLACK) ;
  ac_free (h1) ;
  ac_free (h) ;
  arraySort (look->map->genes, htileSegOrder) ;
} /* htilePBS */

/************************************************************/
/************************************************************/
#ifdef BOU_DEF

/* draw replication start site */
FREEOPT htileRZoneMenu[] = {
  { 2, "R-zone menu"},
  { 1, "Hide all R-zones"},
  { 2, "Delete this zone"},
  { 0, 0 } } ;

/* manage signal tracks */
static void htileRZoneAction (KEY k, int box)
{
  Htile look = currentHtile("htileWallAction") ;
  
  switch (k)
    {
    case 1:
      look->showRZones = FALSE ;
      break ;
    case 2:
      {
	AC_HANDLE h = ac_new_handle () ;
	vTXT bfr = vtxtHandleCreate (h) ;

	if (box >= look->map->minRZoneBox && box < look->map->maxRZoneBox)
	  {
	    int iZone =  box - look->map->minRZoneBox ;
	    RZS *rzs ;
	    if (iZone >= 0 && iZone < arrayMax (look->map->zones))
	      {
		rzs = arrayp (look->map->zones, iZone, RZS) ;
		vtxtPrintf (bfr, "-D PBS %s\n\n", ac_protect (ac_key_name(rzs->key), h)) ;
		ac_parse (look->db, vtxtPtr (bfr), 0, 0, h) ;
	      }
	  }
	ac_free (h) ;
	break ;
      }
    }
  look->map->mapDraw () ;
} /* htileWallAction */

/************************************************************/

static void htileDoDrawRZone (Htile look, float offset, AC_OBJ Pbs, 
			      AC_TABLE mm, int type)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE rep = 0 ;
  int iZone, box, col ;
  char *phase ;
  float dy = 0 ;
  RZS *rzs;
  KEY intMap = ac_table_key (mm, 0, 0, 0) ;
  int a1 = ac_table_int (mm, 0, 1, 1) ;
  int a2 = ac_table_int (mm, 0, 2, 1) ;
  int x1 = a1 - look->map->a1 ;
  int x2 = a2 - look->map->a1 ;
  float u1 = TMAP2GRAPH (look->map, x1 - .5) ;
  float u2 = TMAP2GRAPH (look->map, x2 + .5) ;
  if (u2 < 0 || u1 > look->map->graphWidth)
    return ;


  switch (type)
    {
    case 1: 
      phase = "S1" ;
      col = GREEN3 ;
      break ;
    case 2: 
      phase = "S2" ;
      col = GREEN5 ;
      break ;
    case 3: 
      phase = "S3" ;
      col = GREEN7 ;
      break ;
    default:
      goto done ;
    }

  if (!(rep = ac_tag_table (Pbs, phase, h)))
    goto done ;
  dy=ac_table_float (rep, 0, 0, 0.1) ;
  iZone = arrayMax (look->map->zones) ;
  array (look->map->boxIndex, box = graphBoxStart (), int) = iZone ;
  rzs = arrayp (look->map->zones, iZone, RZS) ;
  rzs->key = ac_obj_key (Pbs) ;
  rzs->intMap = intMap ;
  rzs->type = type ;
  rzs->a1 = a1 ;
  rzs->a2 = a2 ;
  rzs->u1 = u1 ;
  rzs->u2 = u2 ;
  rzs->dy = dy ;
  rzs->x = u1 ; rzs->y = offset - look->zoom * dy - .5 ;
  graphRectangle (u1, offset - look->zoom *dy - .5, u2, offset - look->zoom *dy + .5) ;
  graphBoxEnd () ;
  graphBoxDraw (box, BLACK, col) ;
  graphBoxFreeMenu (box, htileRZoneAction, htileRZoneMenu) ;

 done:
  ac_free (h) ;

  return ;
} /* htileDoDrawRZone */

/************************************************************/

static void htileDoCreateRZone (float x)
{
  Htile look = currentHtile("htileCreateRZone") ;
  int a1, a2 ;
  const char *buf ;
  int x1 = TGRAPH2MAP(look->map, x) ;
   
  buf = ac_new_name (look->db, "PBS", "zone_%d", look->h) ;
  if (look)
    {
	AC_HANDLE h = ac_new_handle () ;
	vTXT bfr = vtxtHandleCreate (h) ;

	vtxtPrintf (bfr, "PBS %s\n", buf) ;
	a1 = x1 - 100000 ; a2 = x1 + 100000 ;
	vtxtPrintf (bfr, "S1 .1\nS2 .5\nS3 .8\n") ;
	vtxtPrintf (bfr, "IntMap %s %d %d\n", name(look->intMap), a1 + look->map->a1 - 1, a2 + look->map->a1 - 1) ;
	vtxtPrintf (bfr, "Locus %s\n\n", name(look->key)) ;
	vtxtPrintf (bfr, "Locus %s\n", name(look->key)) ;
	vtxtPrintf (bfr, "PBS %s %d %d\n", buf, a1, a2) ;

	ac_parse (look->db, vtxtPtr (bfr), 0, 0, h) ;
	ac_free (h) ;

	look->map->mapDraw () ;
    }
} /* htileDoCreateRZone */
 
/***************/

static void htileCreateRZone (void)
{
  htileDoCreateRZone (graphEventX) ;
} /* htileCreateRZone */


/************************************************************/

static void htileDrawRZone (Htile look, float offset)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  AC_OBJ Locus = ac_get_obj (look->db, className (look->locus), name(look->locus), h) ;
  AC_OBJ Pbs = 0 ;
  AC_TABLE mm, pbsT ;
  KEY pbs ;
  int ir, jr ;

  mm = ac_tag_table (Locus, "IntMap", h) ;
  /*
    g1 = ac_table_int (mm, 0, 1, 1) ;
    g2 = ac_table_int (mm, 0, 2, 1) ;
  */
  arrayDestroy (look->map->zones) ;
  look->map->minRZoneBox = arrayMax(look->map->boxIndex) ;
  pbsT = ac_tag_table (Locus, "PBS", h) ;
  for (ir = jr = 0 ; pbsT && ir < pbsT->rows ; ir++)
    {
      pbs = ac_table_key (pbsT, ir, 0, 0) ;
      if (! keyFindTag (pbs, str2tag("Replication")))
	continue ;
      if (!jr++)
	look->map->zones = arrayHandleCreate (32, RZS, look->h) ;
      ac_free (h1) ;
      h1 = ac_new_handle () ;
      Pbs = ac_table_obj (pbsT, ir, 0, h1) ;
      mm = ac_tag_table (Pbs, "IntMap", h1) ;
      htileDoDrawRZone (look, offset, Pbs, mm, 1) ;
      htileDoDrawRZone (look, offset, Pbs, mm, 2) ;
      htileDoDrawRZone (look, offset, Pbs, mm, 3) ;
    }

  look->map->maxRZoneBox = arrayMax(look->map->boxIndex) ;
  graphColor (BLACK) ;
  ac_free (h1) ;
  ac_free (h) ;
  arraySort (look->map->genes, htileSegOrder) ;
} /* htileDrawRZone */
#endif
/************************************************************/
/************************************************************/

static int htileDrawGeneSignal (Htile look, float offset)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ Gene = 0 ;
  AC_TABLE tbl ;
  SEG *gSeg ;
  PNS *pns ;
  PGN *pgn ;
  KEY gene = 0 ;
  const char *tag1 ;
  char *cp1, mybuf[200] ;
  double z1, zz, signal ;
  float u1, u2, y1, y2, dy ;
  int ii, ns, n1, nz, jr, x1, x2 , geneBox, col ;

  y1 = offset ;
  y2 = y1 + 1 ;
  for (ii = 0 ; look->map->genes && ii < arrayMax(look->map->genes) ; ii++)
    {
      gSeg = arrp (look->map->genes, ii, SEG) ;
      x1 = gSeg->x1 ;
      x2 = gSeg->x2 ;
      u1 = TMAP2GRAPH (look->map, x1) ;
      u2 = TMAP2GRAPH (look->map, x2) ;
      if (u1 > u2) { float u0 = u1 ; u1 = u2 ; u2 = u0 ; }
      if (u1 < 15) u1 = 15 ;
      if (u2 < u1 || u1 > look->map->graphWidth)
	continue ;
      for (ns = 0, dy = -1 ; ns < NSMAX ; ns += look->ratio ? 2 : 1)
	{
	  pns = pnsAll + ns ;
	  if (!pns->nam)
	    break ;
	  pgn = pnsPgn (look->ratio, pns) ;
	  if (
	      (!look->ratio && look->isClosed [ns]) ||
	      (look->ratio && ! (ns & 1) && look->isClosed [ns]) ||
	      (look->ratio &&   (ns & 1) && look->isClosed [ns-1]) 
	      )
	    continue ;
	  if (pns->flag & PGG_PA) 
	    {
	      tag1 = "Specific_exon_signal" ;
	    }
	  else if (pns->flag & (PGG_NUC | PGG_TOTAL))
	    {
	      tag1 = "Specific_intron_signal" ;
	    }
	  else
	    continue ;
	  dy++ ;
	  Gene =  htile_key_obj (look, gene = gSeg->key, h) ;
	  nz = 0 ; zz = 0 ;
	    {
	      tbl = ac_tag_table (Gene, tag1 , 0) ;
	      for (jr = 0 ; tbl && jr < tbl->rows ; jr++)
		{
		  strcpy (mybuf, ac_table_printable (tbl, jr, 0, "")) ;
		  cp1 = strstr (mybuf, ".pair") ;
		  if (cp1) *cp1 = 0 ;
		  if (!strcmp (pns->nam, mybuf) || !strcmp (pns->p, mybuf)) 
		    {
		      z1 = ac_table_float (tbl, jr, 1, -1) ;
		      n1 = ac_table_int (tbl, jr, 2, 0) ;
		      nz += n1 ; zz += n1 * z1 ;
		    }
		}
	      ac_free (tbl) ;
	    }
	  ac_free (Gene) ;
	  if (!nz)
	    continue ;
	  signal = zz/nz + 10.0 - pns->av ; 
	  col = GRAY ; /* not tested */
	  if (pns->flag & PGG_PA) /* pA+ */
	    {
	      if (signal < 1) col = GRAY ; /* not tested */
	      else if (signal < pgn->s20) col = RED1 ;
	      else if (signal < pgn->s40) col = RED2 ;
	      else if (signal < pgn->s60) col = RED3 ;
	      else if (signal < pgn->s80) col =  RED5 ;
	      else if (1) col = RED6 ;
	      else continue ;
	    }
	  else if (pns->flag & PGG_TOTAL) /* Total */
	    {
	      if (signal < 1) col = GRAY ; /* not tested */
	      else if (signal < pgn->s20) col = BLUE1 ;
	      else if (signal < pgn->s40) col = BLUE2 ;
	      else if (signal < pgn->s60) col = BLUE4 ;
	      else if (signal < pgn->s80) col =  BLUE6 ;
	      else if (1) col = BLUE ;
	      else continue ;
	    }
	  else /* pA- */
	    {
	      if (signal < 1) col = GRAY ; /* not tested */
	      else if (signal < pgn->s20) col = GREEN1 ;
	      else if (signal < pgn->s40) col = GREEN2 ;
	      else if (signal < pgn->s60) col = GREEN3 ;
	      else if (signal < pgn->s80) col =  GREEN5 ;
	      else if (1) col = GREEN6 ;
	      else continue ;
	    }

	  y1 = offset  + dy ; y2 = y1 + .7 ;
	    {
	      array (look->map->boxIndex, geneBox = graphBoxStart (), int) = gene ;
	      graphRectangle (u1, y1, u2, y2) ;
	      graphBoxEnd () ;
	      graphBoxDraw (geneBox, /* gSeg->col */ col, col) ;
	    }
	}
    }
  graphColor (BLACK) ;
  for (ns = 0, dy = -1 ; ns < NSMAX ; ns += look->ratio ? 2 : 1)
    {
      pns = pnsAll + ns ;
      if (!pns->nam)
	break ;
      if (
	  (!look->ratio && look->isClosed [ns]) ||
	  (look->ratio && ! (ns & 1) && look->isClosed [ns]) ||
	  (look->ratio &&   (ns & 1) && look->isClosed [ns-1]) 
	  )
	continue ;
      if (strstr(pns->nam, "pA+")) ;
      else if (pns->isRna) ;
      else continue ;
      dy++ ;
      graphText (look->ratio ? pns->nam2 : pns->nam, 1, offset  + dy) ;
    }      
  ac_free (h) ;
  return dy ;
} /* htileDrawGeneSignal */

/************************************************************/
/************************************************************/

static void htileDrawGenome (Htile look, float offset)
{
  AC_HANDLE h1 = 0, h = 0 ;
  const char *dna, *ccp ;
  AC_KEYSET cosmids = 0 ;
  AC_TABLE cosmidTable = 0 ;
  AC_OBJ Region = 0, Cosmid = 0 ;
  int ir, i, a1, a2, b, b1, b2,  x1, x2 ;
  int w = 1000 ; /* average over 1000 bp */
  int na, nt, nc, ng, ncg, zdna, p ;
  float ra, rt, rc, rg, rcg ;  
  float oldx = 0, x, y ;
  BOOL first = TRUE ;
  HTYPE *tt ;
  const char* qq ;

  h = ac_new_handle () ;
  x1 = look->map->a1 + TGRAPH2MAP (look->map,  1) ;
  x2 = look->map->a1 + TGRAPH2MAP (look->map,  look->map->graphWidth) ;

  /* average over at least 20 bp */
  w = look->gaussWidth ;
  if (w < 20) w = 20 ;

  Region = htile_key_obj (look, look->locus, h) ;
  cosmids = ac_objquery_keyset (Region, "{CLASS mRNA} $| {CLASS Sequence} $| {>Sequence} ; {!Parts} $| {>Parts}", h) ;
  qq = "select c,m,a1,a2 from c in @, m in c->intmap, a1 in m[1], a2 in a1[1]" ;
  cosmidTable = ac_bql_table (look->db, qq, cosmids, "+2+3+4+1", 0, h) ;
  for (ir = 0 ; cosmidTable && ir < cosmidTable->rows ; ir++, ac_free (h1))
    {
      h1 = ac_new_handle () ;
      Cosmid = ac_table_obj (cosmidTable, ir, 0, h1) ;
      if (ac_table_key (cosmidTable, ir, 2, 0))
	{
	  a1 = ac_table_int (cosmidTable, ir, 2, 0) - look->map->a1 + 1 ;
	  a2 = ac_table_int (cosmidTable, ir, 3, 0) - look->map->a1 + 1 ;
	}
      else
	{
	  a1 = look->map->a1 ; 
	  a2 = look->map->a2 ;
	}
      if (a2 < x1 || a1 > x2) /* area to be shown */
	continue ;

      dna = ac_obj_dna (Cosmid, h1) ;
      if (!dna) continue ;
      /* base stretch to be analysed */
      b1 = x1 - a1 + 1 ;
      if (b1 < 0) b1 = 0 ;
      b2 = x2 - a1 + 1 ;
      if (b2 > a2 - a1)
	b2 = a2 - a1 ;
      x1 = b2 + a1 ; /* start of next tile */
      if (0)
	{
	  char *c ;
	  for (b = b1 ; b < b2 ; b++)
	    {
	      if ((c = strstr(dna+b, "gatct")))
		{
		  b = c - dna ; 
		  c+= 6 ;
		  if (*c++=='t' && *c++=='t' && *c++>0 && *c++=='t' && *c++=='t' && *c++=='t' && *c++=='t')
		    {
		      graphColor (RED) ;
		      x = TMAP2GRAPH (look->map,  b + a1 - look->map->a1) ;
		      graphLine (x, 1, x, 100) ;
		      graphColor (BLACK) ;
		    }
		}
	      else if ((c = strstr(dna+b, "aaaa")))
		{
		  b = c - dna ; 
		  c+= 5 ;
		  if (*c++=='a' && *c++=='a' && *c++>0 && *c++=='a' && *c++=='g'
		      && *c++=='a' && *c++=='t' && *c++=='c')
		    {
		      graphColor (RED) ;
		      x = TMAP2GRAPH (look->map,  b + a1 - look->map->a1) ;
		      graphLine (x, 1, x, 100) ;
		      graphColor (BLACK) ;
		    }
		}
	      else
		break ;
	    }
	}
      for (b = b1+w/2 ; b < b2 ; b += w)
	{
	  zdna = p = na = nc = ng = nt = ncg = 0 ;
	  ccp = dna + b + 1 -w/2 ;
	  for (i = 0 ; i < w && b+i < b2 ; ccp++, i++)
	    {
	      switch (*ccp)
		{
		case 'g': ng++ ; break ;
		case 'c': nc++ ; if(*(ccp+1)=='g') ncg++ ; break ;
		case 'a': na++ ; break ;
		case 't': nt++ ; break ;
		}  
	      if (!zdna)
		switch (*ccp)
		  {
		  case 'a': p = 0 ; break ;
		  case 'g': if(p>0 && (p & 0x1)) {p++ ; if (p>18) zdna++;} else p = 0 ; break ;
		  case 't': case 'C': if(p>-1 && !(p & 0x1)) {p++ ; if (p>18) zdna++;} else p = 1 ; break ;
		  }
	    }
	  i = na + nt + nc + ng ;
	  if (i < w/10)  continue ;
	  zdna *= 100 ;
	  ra = (100.0 * na)/i ; rt =  (100.0 * nt)/i ; 
	  rg =  (100.0 * ng)/i ; rc =  (100.0 * nc)/i ; rcg =  (200.0 * ncg)/i ;
	  x = TMAP2GRAPH (look->map,  b + w/2 + a1 - 1 - look->map->a1) ;
	  for (tt = histoTypes ; tt->nam ; tt++)
	    {
	      if (! tt->open)
		continue ;
	      switch (tt->type)
		{
		case RA: y = ra ; break ;
		case RT: y = rt ; break ;
		case RG: y = rg ; break ;
		case RC: y = rc ; break ;
		case RGpC: y = rg + rc ; break ;
		case RApG: y = ra + rg ; break ;
		case RCpT: y = rt + rc ; break ;
		case RGmC: y = rg - rc ; break ;
		case RTpG: y = rt + rg ; break ;
		case RCG: y = rcg ; break ;
		case ZDNA: y = zdna ; break ;
		default: y = 0 ; break ;
		}
	      y -= (tt->min + tt->max)/2.0 ;
	      y = offset -  20.0 * y/ (tt->max - tt->min) ;
	      graphColor (tt->col) ;
	      if (oldx > 10)
		graphLine (oldx, tt->y, x, y) ;
	      if (oldx > 10 && first)
		graphText (tt->nam, 2, tt->y) ; 
	      tt->y = y ;
	    }
	  if (oldx > 10 && first)
	    first = FALSE ;
	  oldx = x ;
	}
  }
  graphColor (BLACK) ;
  ac_free (h1) ;
  ac_free (h) ;
} /* htileDrawGenome */

/************************************************************/
/************************************************************/

static void htileGetMask (Htile look, int maskLength)
{ 
  AC_HANDLE h ;
  MASK *mm ;
  int ii=0, a1, a2 ;
  Array aa, mask ; 
  ACEIN ai ;
  
  mask = look->mask[maskLength] ;
  if (mask) return ;
  h = ac_new_handle () ;
  ai = aceInCreateFromFile (messprintf ("SOLEXA/Mask%d.%s", maskLength, name (look->intMap)), "r", 0, h) ;
  if (ai)
    {
      aa = arrayHandleCreate (10000, MASK, look->h) ;
      while (ai && aceInCard (ai))
	{
	  if (aceInWord (ai) &&
	      aceInInt (ai, &a1) &&
	      aceInInt (ai, &a2)
	      )
	    {
	      mm = arrayp (aa, ii++, MASK) ;
	      mm->a1 = a1 ; 
	      mm->a2 = a2 ; 
	    }
	}
    }
  else
    aa = (Array)1;
  look->mask[maskLength] = aa ;
  ac_free (h) ;
} /* htileGetMask */

/************************************************************/
#ifdef BOU_DEF

static void htileDrawMaskLength (Htile look, float offset, int maskLength)
{
  float dy = 10 ; /* Height of the histograms */
  int ii, a1, a2, b1, b2 ;
  float x1, x2 ;
  MASK *mm ;
  Array mask ;

  htileGetMask (look, maskLength) ;
  mask = look->mask[maskLength] ;
  if (mask == 0 || mask == (Array)1) return ;

  x1 = TGRAPH2MAP (look->map,  1) ;
  x2 = TGRAPH2MAP (look->map,  look->map->graphWidth) ;

  a1 = x1 + look->map->a1 - 1 ; /* int mapp coords of the visible region */
  a2 = x2 + look->map->a1 - 1 ;

  if (look->smoothing == 0)
    {
      for (ii = 0, mm = arrp (mask, 0, MASK) ; ii < arrayMax (mask) ; mm++, ii++)
	{
	  if (mm->a2 < a1) continue ;
	  if (mm->a1 > a2) break ;
	  b1 = mm->a1 > a1 ? mm->a1 : a1 ;
	  b2 = mm->a2 < a2 ? mm->a2 : a2 ; 
	  x1 = TMAP2GRAPH (look->map, b1 - look->map->a1) ;
	  x2 = TMAP2GRAPH (look->map, b2 - look->map->a1) ;
	  
	  if (x1 < x2 - .5)
	    graphFillRectangle (x1, offset+dy/3, x2, offset) ;
	  else
	    graphLine ((x1+x2)/2, offset, (x1+x2)/2, offset - dy) ;
	}
    }
  else
    {
      int jj, kk, c1, c2, step = look->gaussWidth/10 ;
      Array a = arrayCreate (10000, int) ;
      Array aa = arrayCreate (10000, float) ;
      float dy1, dy2, zz[11] = {1,3,6, 8, 9, 10, 9, 8, 6, 3, 1} ;

      for (ii = 0, mm = arrp (mask, 0, MASK) ; ii < arrayMax (mask) ; mm++, ii++)
	{
	  if (mm->a2 < a1) continue ;
	  if (mm->a1 > a2) break ;
	  b1 = mm->a1 > a1 ? mm->a1 - a1 : 0 ;
	  b2 = mm->a2 < a2 ? mm->a2 - a1 : a2 - a1 ; 
	  for (jj = b1/step ; jj < b2/step + 1 ; jj++)
	    {
	      c1 = b1 > jj * step ? b1 : jj * step ;
	      c2 = b2 < (jj+1) * step ? b2 : (jj+1) * step ;
	      if (c1 < c2)
		array (a, jj , int) += c2 - c1 ;
	    }
	}
      /* smooth the curve */
      for (jj = 0 ; jj < arrayMax (a) ; jj++)
	{
	  dy1 = 0 ; dy2 = 0 ;
	  for (kk = -5 ; kk <= 5 ; kk++)
	    if (jj + kk >= 0 && jj + kk < arrayMax(a))
	      {
		dy1 += zz[kk + 5] * ((float)arr (a, jj + kk, int) )/step ;
		dy2 += zz[kk + 5] ;
	      }
	  if (dy2 > 0) dy1 /= dy2 ;
	  array (aa, jj, float) = dy1 ;
	}
      for (x2 = dy2 = 0, jj = 0 ; jj * step < a2 - a1 ; jj++)
	{
	  dy1 = dy2 ;
	  x1 = x2 ;
	  dy2 = 20 * array (aa, jj, float) ;
	  x2 =  TMAP2GRAPH (look->map, jj * step + a1 - look->map->a1) ;
	  if (jj)
	    graphLine (x1, offset - dy1, x2, offset - dy2) ;
	}
      ac_free (a) ;
      ac_free (aa) ;
      x1 =  TMAP2GRAPH (look->map, a1 - look->map->a1) ;
      x2 =  TMAP2GRAPH (look->map, a2 - look->map->a1) ;
      dy1 = 20 ;
      graphLine (x1, offset - dy1, x2, offset - dy1) ; 
    }
} /* htileDrawMaskLength */

/************************************************************/

static void htileDrawMask (Htile look, float offset)
{
  /*
    graphColor (GREEN) ;
    htileDrawMaskLength (look, offset, 18) ;
  */
  graphColor (BLUE) ;
  htileDrawMaskLength (look, offset, 25) ;
  graphColor (RED) ;
  htileDrawMaskLength (look, offset, 35) ;

  graphColor (BLACK) ;
} /* htileDrawMask */
#endif
/************************************************************/
FREEOPT htileWallMenu[] = {
  { 3, "Create zone"},
  { 2, "walls"},
  { 1, "Comment"},
  { 2, "Delete"},
  { 0, 0 } } ;

/* manage signal tracks */
static void htileWallAction (KEY k, int box)
{
  Htile look = currentHtile("htileWallAction") ;
  
  switch (k)
    {
    case 1:
      htileAddWallComment () ;
      break ;
    case 2:
      {
	AC_HANDLE h = ac_new_handle () ;
	vTXT bfr = vtxtHandleCreate (h) ;
	Htile look = currentHtile("htileRecalculate") ;
	char *cp ;

	vtxtPrintf (bfr, "Locus %s\n", ac_protect (name(look->locus), h)) ;
	cp = array (look->wallBufs, box - look->map->minWallBox, char *) ;
	vtxtPrintf (bfr, "-D %s\n\n", cp) ;
	ac_parse (look->db, vtxtPtr (bfr), 0, 0, h) ;
	
	ac_free (h) ;
	}
      look->map->mapDraw () ;
      break ;
    case 3:
      {
	AC_HANDLE h = ac_new_handle () ;
	vTXT bfr = vtxtHandleCreate (h) ;
	Htile look = currentHtile("htileRecalculate") ;
	char *cp ;

	vtxtPrintf (bfr, "PBS %s\n", ac_protect (name(look->locus), h)) ;
	cp = array (look->wallBufs, box - look->map->minWallBox, char *) ;
	vtxtPrintf (bfr, "-D %s\n\n", cp) ;
	ac_parse (look->db, vtxtPtr (bfr), 0, 0, h) ;
	
	ac_free (h) ;
	}
      look->map->mapDraw () ;
      break ;
    }
} /* htileWallAction */

/************************************************************/

static void htileWallBoxDrag (float *x, float *y, BOOL isDone)
{
  static float oldY ;
  static BOOL isDragging = FALSE ;

  if (isDragging)
    *y = oldY ;
  else
    { oldY = *y ;
      isDragging = TRUE ;
    }

  if (isDone)
    { 
      htileDoAddWall (*x, htileDraggedBox) ;
      isDragging = FALSE ;
    }

  return;
} /* htileWallBoxDrag */

/************************************************************/

static void htileRZoneMove (float x, float y, int box)
{
  AC_HANDLE h = ac_new_handle () ;
  Tmap *map = tmapCurrent ("tmapPick") ;
  int iZone ;
  RZS *rzs ;
  float dy, dx ;
  vTXT txt = vtxtHandleCreate (h) ;
  Htile look = currentHtile("htileRZoneMove") ;

  if (box > 0 &&
      box < arrayMax(map->boxIndex) && 
      (iZone = arr (map->boxIndex, box, int)) >= 0 &&
      arrayExists (map->zones) &&
      iZone >= 0 &&
      iZone < arrayMax (map->zones)
      )
    {
      rzs = arrayp (map->zones, iZone, RZS) ;
      dy = rzs->dy - (y - rzs->y)/look->zoom ;
      dx = (x - rzs->u1)/map->mag ;
      if (rzs->x==1)  rzs->a2 += dx ;
      else  rzs->a1 += dx ;
      vtxtPrintf (txt, "PBS %s\n", ac_protect(ac_key_name (rzs->key), h)) ;
      switch (rzs->type)
	{   
	case 1:
	  vtxtPrintf (txt, "S1 %f\n", dy) ;
	  break ;
	case 2:
	  vtxtPrintf (txt, "S2 %f\n", dy) ;
	  break ;
	case 3:
	  vtxtPrintf (txt, "S3 %f\n", dy) ;
	  break ;
	} 
      vtxtPrintf (txt, "IntMap %s %d %d\n", name(rzs->intMap), rzs->a1, rzs->a2) ;
      ac_parse (look->db, vtxtPtr (txt), 0, 0, h) ;
    }
  ac_free (h) ;
  map->mapDraw () ;
} /* htileRZoneMove */

/************************************************************/

static void htileRZoneBoxDrag (float *x, float *y, BOOL isDone)
{
  static float oldY ;
  static BOOL isDragging = FALSE ;

  if (0 && isDragging)
    *y = oldY ;
  else
    { oldY = *y ;
      isDragging = TRUE ;
    }

  if (isDone)
    { 
      htileRZoneMove (*x, *y, htileDraggedBox) ;
      isDragging = FALSE ;
    }

  return;
} /* htileRZoneBoxDrag */

/************************************************************/

static void htileDrawWalls (Htile look, float offset)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ Type, Region = 0 ;
  AC_TABLE tbl = 0 ;
  int box, ir, a1, col ;
  BOOL first = TRUE ;
  float x, y1, y2, width = .2, oldWidth ;
  const char *comment = "toto" ;
  char *cp ;

  y1 = offset ; y2 = look->map->graphHeight - .5 ;

  oldWidth = graphLinewidth (.2) ;	
  Region = htile_key_obj (look, look->locus, h) ;
  tbl = ac_tag_table (Region, "Tiling_wall", h) ;
  look->wallBufs = arrayHandleCreate (12, char*, look->h) ;
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      Type = ac_table_obj (tbl, ir, 0, h) ;
      if (Type)
	{
	  col = ac_tag_key (Type,  "Colour", _RED) - _WHITE + WHITE ;
	  width = ac_tag_float (Type,  "Width", .2) ;
	}
      else
	col = RED ;
      a1 = ac_table_int (tbl, ir, 2, 0) ;
      comment = ac_table_printable (tbl, ir, 3, 0) ;
      if (comment)
	comment = ac_protect(comment, h) ;
      if (a1)
	{
	  a1 -= look->map->a1 - 1 ;
	  x = TMAP2GRAPH (look->map, a1) ;
	  box = graphBoxStart () ;
	  if (first)
	    look->map->minWallBox = box ;
	  look->map->maxWallBox = box + 1 ;
	  first = FALSE ;
	  cp = hprintf (look->h, "Tiling_wall %s %s %d %s", ac_name (Type), name(look->intMap), a1 + look->map->a1 - 1, comment ? comment : "") ;
	  array (look->wallBufs, box - look->map->minWallBox, char *) = cp ;	    
	  graphColor (col) ;
	  graphLinewidth (width) ;	
	  graphLine (x, y1, x, y2) ;
	  graphColor (TRANSPARENT) ;
	  graphLine (x-.7, y2-1, x+.7, y2-1) ;
	  graphBoxEnd () ;
	  graphBoxDraw (box, BLACK, TRANSPARENT) ;
	  graphBoxFreeMenu (box, htileWallAction, htileWallMenu) ;
	}
    }
  graphLinewidth (oldWidth) ;
  graphColor (BLACK) ;
  ac_free (h) ;
  return ;
} /* hTileDrawWalls */

/************************************************************/

static void hTileGeneSearch (char *gene)
{
  AC_HANDLE h = ac_new_handle () ;
  Htile look = currentHtile("htileGeneSearch") ;
  AC_ITER  iter ;
  AC_KEYSET ks ;
  AC_OBJ Gene = 0 ;
  BOOL same = TRUE ;

  if (gene && *gene)
    {
      ks = ac_dbquery_keyset (look->db, messprintf ("Find Locus %s ; >gene ; IS \"%s\"" 
						    , ac_protect (name (look->locus), h)
						    , gene
						    ), h) ;
      if (!ac_keyset_count (ks))
	ks = ac_dbquery_keyset (look->db, messprintf ("Find Locus %s ; >gene ; IS \"%s*\"" 
						    , ac_protect (name (look->locus), h)
						    , gene
						    ), h) ;
      if (!ac_keyset_count (ks))
	ks = ac_dbquery_keyset (look->db, messprintf ("Find Locus %s ; >gene ; IS \"*%s*\"" 
						    , ac_protect (name (look->locus), h)
						    , gene
						    ), h) ;
      if (!ac_keyset_count (ks)) same = FALSE ;
      if (!ac_keyset_count (ks))
	ks = ac_dbquery_keyset (look->db, messprintf ("Find gene \"%s*\"; COUNT {>Locus ; Sequence} > 0"), h) ;
      if (!ac_keyset_count (ks))
	ks = ac_dbquery_keyset (look->db, messprintf ("Find gene \"*%s*\"; COUNT {>Locus ; Sequence} > 0"), h) ;
      if (!ac_keyset_count (ks))
	ks = ac_dbquery_keyset (look->db, messprintf ("Find gene \"%s\"; COUNT {>Locus ; Sequence} > 0"), h) ;

      if (!ac_keyset_count (ks))
	ks = ac_dbquery_keyset (look->db, messprintf ("Find gene \"%s*\"; COUNT {>Locus ; Sequence} > 0"), h) ;
      if (!ac_keyset_count (ks))
	ks = ac_dbquery_keyset (look->db, messprintf ("Find gene \"*%s*\"; COUNT {>Locus ; Sequence} > 0"), h) ;

      if (!ac_keyset_count (ks))
	messout ("Sorry, i do not know gene %s", gene) ;
      else
	{
	  keySetDestroy (look0->activeGenes) ;
	  look0->activeGenes = keySetCreate () ;
	  iter = ac_keyset_iter (ks, FALSE,  h) ;
	  while (ac_free (Gene), Gene = ac_iter_obj (iter))
	    {
	      KEY gene = ac_obj_key (Gene) ;
	      keySetInsert (look0->activeGenes, gene) ;
	    }
	  look->from = 0 ;  
	  if (same)
	    look->map->mapDraw () ;
	}
  }
  ac_free (h) ;
} /* hTileGeneSearch */

/************************************************************/

static void hTileCoordSearch (char *txt)
{
  Htile look = currentHtile("htileCoordSearch") ;
  int a1, x ;
  float a0 ;
  
  if (sscanf(txt, "%d", &x))
    {
      a1 = x - look->map->a1 ;
      a0 = TGRAPH2MAP (look->map, look->map->graphWidth/2.0) ;
      look->map->offset += a1 - a0 ;
      look->from = 0 ;  
      look->map->mapDraw () ;
    }
} /* hTileCoordSearch */

/************************************************************/

static void htileDraw (void)
{
  float xx = 2 ;
  float offset, y0, wiggleOffset ; 
  int a0, pass ;
  Htile look = currentHtile("htileDraw") ;

  if (look->from) /* centre */
    {
      AC_HANDLE h = ac_new_handle () ;
      AC_OBJ Gene = htile_key_obj (look, look->from, h) ;
      AC_TABLE tbl = ac_tag_table (Gene, "IntMap", h) ;
      int a1 = ac_table_int (tbl, 0, 1, 0) ;
      int a2 = ac_table_int (tbl, 0, 2, 0) ;
      int a0, centre ;

      if (! Gene)
	{
	  int da = 3000 ;  /* initial zone-widht of the wiggle */
	  centre = look->from + look->map->min - look->map->a1 ;
	  look->map->mag = (3.0 * look->map->graphWidth) / (7.0 * da) ;
	  /* center the center */
	  a0 = TGRAPH2MAP(look->map, look->map->graphWidth/2.0) ;
	  look->map->offset += centre - a0 ;
	  a0 = TGRAPH2MAP(look->map, look->map->graphWidth/2.0) ;
	}
      else if (!strcasecmp (ac_class (Gene), "Gene") || !strcasecmp (ac_class (Gene), "PBS"))
	{
	  if (a1 > a2) { int a3 = a1 ; a1 = a2 ; a2 = a3 ; }
	  if (a1 && a2 && a2 > look->map->a1 && a1 < look->map->a2)
	    {
	      int da = a2 - a1 ; if (da < 5000) da = 5000 ;
	      centre = (a1 + a2)/2 + look->map->min - look->map->a1 ;
	      look->map->mag = (3.0 * look->map->graphWidth) / (7.0 * da) ;
	      if(!strcasecmp (ac_class (Gene), "PBS"))
		look->map->mag /= 10 ;
	      /* center the center */
	      a0 = TGRAPH2MAP(look->map, look->map->graphWidth/2.0) ;
	      look->map->offset += centre - a0 ;
	      a0 = TGRAPH2MAP(look->map, look->map->graphWidth/2.0) ;
	    }
	  look->showGenes = TRUE ;
	  look->from = 0 ;
	}
      ac_free (h) ;
    }

  htileMenuSelect (look) ;
  graphClear () ;

  for (wiggleOffset = 1, pass = 0 ; pass < 2 ; pass++)
    {
      /* in pass 0 compute wiggleOffset */
      if (pass) /* actually draw on pass 1 */
	{
	  if (look->smoothing >= 5)
	    htileGeneSmoothingDraw (look, wiggleOffset) ;
	  else
	    htileDrawProbes (look, wiggleOffset, y0, 0) ;
	  if (look->map->solexa)
	    htileDrawSolexa (look, wiggleOffset, y0, 0) ;
	  if (look->showGenome)
	    htileDrawGenome (look, wiggleOffset - 20) ;

	  graphMenu (arrp(look->menu, 0, MENUOPT)) ;
	  strncpy (look->geneSearchBuf, name(look->key), 300) ;
	  offset = 1 ;
	  graphText ("Region:", xx, offset) ; xx += 7 ;
	  graphTextScrollEntry (look->geneSearchBuf, 300, 20, xx, offset, hTileGeneSearch) ; xx += 23 ;
	  graphText ("Position:", xx, offset) ; xx += 9 ;
	  a0 = TGRAPH2MAP (look->map, look->map->graphWidth/2.0) + look->map->a1 + .49 ;
	  sprintf (look->coordSearchBuf, "%d", a0) ;
	  graphTextScrollEntry (look->coordSearchBuf, 300, 20, xx, offset, hTileCoordSearch) ; xx += 21 ; 
	  /* draw cosmid */
	  xx = 1 ;
	  
	  arrayMax(look->map->boxIndex) = 0 ;
	  offset = 4 ;
	  if (!look->hideHeader)
	    {
	      htileDrawActionButtons (look, &offset) ;
	      offset += 2.3 ;
	      htileDrawGenomeButtons (look, &offset) ;
	      offset += 2.3 ;
	      htileDrawTileGroupButtons (look, &offset) ;
	      offset += 2.3 ;
	    }
	  htileDrawTileButtons (look, &offset) ;
	  offset += 2.3 ;
	}
      else
	wiggleOffset = 4 + 2.3 * 3 ;
      
      if (look->showCaptureProbes)
	{
	  if (pass)
	    {
	      htileDrawCaptureProbes (look, offset + 6.2, TRUE) ;
	      offset += 8 ;
	    }
	  else
	    wiggleOffset += 8 ;
	}
      if (look->showGenes)
	{
	  if (pass)
	    {
	      htileDrawGenes (look, offset + 6.2, TRUE) ;
	      offset += 8 ;
	    }
	  else
	    wiggleOffset += 8 ;
	}
      if (pass)
	{
	  if (! look->hideScale)
	    htileDrawScale (look, offset + .5) ;
	  if (! look->hideScale && ! isGifDisplay)
	    htileDrawDna (look, offset + 1.5) ;
	  offset += 2.5 ;
	}
      else
	wiggleOffset += 2.5 ;
      
      if (look->showGenes)
	{
	  if (pass)
	    {
	      htileDrawGenes (look, offset, FALSE) ;
	      offset += 6 ;
	    }
	  else
	    wiggleOffset += 6 ;
	}
      if (pass && look->showGenes)
	htileDrawPBS (look, offset - 2, TRUE) ;
      if (look->showGenes)
	{
	  if (pass)
	    {
	      htileDrawGeneSignal (look, offset) ;
	      offset += 8 ;
	    }
	  else
	    wiggleOffset += 8 ;
	}
      
      
      /* offset is the position of the ZERO line */
      y0 = offset + 8 ;
      look->map->dragTopMargin  = y0 ;
      
      if (look->ratio)
	{
	  offset += (look->map->graphHeight  - offset - 6)/2 + offset ; /* was offset += 28 */
	  wiggleOffset += (look->map->graphHeight  - offset - 6)/2 + offset ; /* was offset += 28 */
	}
      else
	 wiggleOffset = offset = look->map->graphHeight - 3 ;

    }

#ifdef BOU_DEF
  if (look->showRZones)
    htileDrawRZone (look, offset) ;
#endif
  if (0)   htileDrawMusicalScale (look, offset) ;
#ifdef BOU_DEF
  if (look->showMask)
    htileDrawMask (look, offset) ;
#endif
  if (look->showExtrema)
    htileDrawWalls (look, y0 + 1) ;
  xx = look->map->graphWidth/2.0 ;
  graphColor (GRAY) ;
  if (look->showGenes && ! look->hideHeader)
    graphLine (xx, y0, xx, look->map->graphHeight) ;
  graphColor (BLACK) ;
  if (1)   htileDrawMusicalScale (look, offset) ;
  graphRedraw () ;
  graphRegister (KEYBOARD, tmapKbd) ;
}

/************************************************************/

static BOOL solexaAllInit (Htile look)
{
  AC_HANDLE h = ac_new_handle () ;
  static Array aa = 0 ;
  BOOL ok = FALSE ;
  AC_TABLE tbl = 0 ;
  AC_OBJ Run = 0, Ali = 0 ;
  AC_ITER iter = 0 ;
  KEY run ;
  KEY _pA = str2tag ("polyA") ;
  KEY _Cap = str2tag ("Cap") ;
  KEY _Total = str2tag ("Total") ;
  int ir, j ;
  PNX *pnx ;

  if (! arrayExists (look->solexaAll))
    {
      int NF = 14 ;
      KEY col ;
      unsigned int flag[14] =  { PGG_ELF, PGG_ERF, PGG_ELR, PGG_ERR, PGG_nuf, PGG_nur, PGG_ppf, PGG_ppr, PGG_uf, PGG_ur, PGG_endRatioLF, PGG_endRatioRF, PGG_endRatioLR, PGG_endRatioRR  } ;
      const char *suffix[14] = {".LF", ".RF", ".LR", ".RR", "nu+", "nu-", "pp+", "pp-", "+", "-", ".eLF", ".eRF", ".eLR", ".eRR" } ;
      const int colors[14] = { BROWN, GREEN, MAGENTA, CYAN, PALEORANGE, YELLOW, BLACK, GRAY, ORANGE, DARKBLUE, PALEGREEN, PALEVIOLET, PALEMAGENTA, PALECYAN  } ;
      aa = look->solexaAll = arrayHandleCreate (128, PNX, look->h) ;
      iter = ac_dbquery_iter (look->db, "find run w_colour || wiggle", h) ;
      ir = -1 ;
      while (ac_free (Run), ir++, Run = ac_next_obj (iter))
	{
	 for (j = 0 ; j < NF ; j++)
	   {
	     pnx = arrayp (aa, NF * ir + j, PNX) ;
	     run = ac_obj_key (Run) ;
	     pnx->p = ac_name (Run) ;    /* this is the name in Map->Wiggle */
	     pnx->nam =  hprintf (look->h, "%s%s", ac_name (Run), suffix[j]) ;  /* this is the nickname for the button */
	     pnx->flag = PGG_rna | flag[j] ;
	     if (keyFindTag (run, _pA)) pnx->flag |= PGG_PA ;
	     if (keyFindTag (run, _RNA)) pnx->flag |= PGG_rna ;
	     if (keyFindTag (run, _Total)) pnx->flag |= PGG_TOTAL ;
	     if (keyFindTag (run, _Cap)) pnx->flag |= PGG_CAP ;
	     
	     pnx->s99 = 1 ;
	     if ( 
		 (Ali = ac_tag_obj (Run, "Ali", h)) &&
		 (tbl = ac_tag_table (Ali, "Accepted", h)) &&
		 tbl->cols >= 5
		  )
	       pnx->s99 = ac_table_float (tbl, 0, 4, 1.0e+9)/1.0e+9 ; /* count in accepted Tb, since z is in kb */
	     col = 0 ;
	     if (j >= 18)
	       {
		 if (flag[j] & PGG_f) col = keyGetKey (run, str2tag ("W_colour_plus")) ; 
		 if (flag[j] & PGG_r) col = keyGetKey (run, str2tag ("W_colour_minus")) ; 
		 if (flag[j] & PGG_ns) col = keyGetKey (run, str2tag ("W_colour_ns")) ; 
	       }
	     pnx->col = col ? (WHITE + col - _WHITE) : colors[j] ;
	     pnx->open = TRUE ;
	     pnx->maxIndex = 1000 ;
	     if (NF * ir + j < NSMAX)
	       look0->isSolexaClosed [NF*ir + j] = ! pnx->open ;
	   }
	}
      pnx = arrayp (aa, NF * ir + 1, PNX) ; /* zero terminate with 2 empty lines */
    }
  ac_free (h) ;
  return ok ;
} /*  solexaAllInit */

/************************************************************/
/*************************************************************/

BOOL htileDisplay (KEY key, KEY from, BOOL isOldGraph)
{
  Htile  oldlook = 0, look ;
  AC_HANDLE handle = 0 ;
  KEY key0 = key ;
  static BOOL firstPass = TRUE ;
  KEYSET ks ;
		/* centre on parent object if poss */
  if (!key)
    return FALSE ;

  cDNAAlignInit () ;
  if (!ac_db)
    ac_db = ac_open_db (0,0) ;

  /* we had in the first part SETELSE {>Locus }  for boubou in 2006-2008 */
  ks = queryKey (key, "{{CLASS sequence OR CLASS Locus} SETELSE {(CLASS Gene || CLASS PBS)  ; >Genomic_sequence}}SETELSE {CLASS mRNA} ") ;
  key = keySet (ks, 0) ;
  if (!key)
    return FALSE ;

  if (
      isOldGraph && 
      graphAssFind (&Htile_MAGIC, &oldlook) &&
      oldlook->magic == &Htile_MAGIC &&
      oldlook->key == key &&
      oldlook->map &&
      graphActivate (oldlook->map->graph)
      )
    { 
      oldlook->from = from + oldlook->map->a1 ;
      tmapResize () ;
      return TRUE ;
    }
  if (oldlook) 
    htileDestroy () ;

  handle = handleCreate () ;
  if (firstPass)
    {
      firstPass = FALSE ;
      look0->smoothing = 0 ;
      look0->ratio = 0 ;
      look0->showDot = 3 ;  /* 0: line , 1: line, 2: bar, 3: stack */
      look0->showExtrema = FALSE ;
      look0->rejectAmbiguous = 0 ;
      look0->romainSmoothing = 0 ;
      if (!strcasecmp(className (from)," PBS"))
	{ look0->gaussWidth = 200 ;  look0->zoom = 5 ; }
      else if (from && class (from) == 0)
	{ look0->gaussWidth = 10 ;  look0->zoom = 5/32.0  ; }
      else if (class (key) == _VSequence)
	{ look0->gaussWidth = 1000 ;  look0->zoom = 5/32.0 ; }
      else if (class (key) == _VmRNA)
	{ look0->gaussWidth = 20 ;  look0->zoom = 5 ; }
      else
	{ look0->gaussWidth = 100000 ; look0->zoom = 40 ; }
      look0->zoom = 5/256.0 ; 
      if (isGifDisplay)
	{ 
	  look0->showRZones = FALSE ; 
	  look0->gaussWidth = 10 ; 
	  look0->zoom = 10.0/128 ; 
	  oldHideHeader = TRUE ; 
	  oldLnWidth = .3 ;
	} 
      look0->unique = 0 ;
      look0->smoothing = 0 ;
      look0->index = FALSE ;
      look0->tmFilter = 0 ;
      look0->exonicFilter = 0 ;
      look0->filter99 = 0 ;
      look0->showMask = FALSE ;
      look0->showGenome = FALSE ;
      look0->showRZones = FALSE ;
      look0->showGenes = TRUE ;
      look0->showCaptureProbes = TRUE ;
      look0->showGeneSignal = TRUE ; 
      {
	int i ;
	switch (look0->ratio)
	  {
	  case 0:
	    for (i = 0 ; i < NSMAX && pnsAll[i].p ; i++)
	      look0->isClosed [i] = ! pnsAll[i].open ;
	    break ;
	  case 1:
	    for (i = 0 ; i < NSMAX && pnsAll[i].p ; i += 2)
	      look0->isClosed [i] = look0->isClosed [i+1] = !(pnsAll[i].open || pnsAll[i+1].open) ;
	    break ;
	  case 2:
	    for (i = 0 ; i < NSMAX && pnsAll[i].p ; i++)
	      look0->isClosed [i] = ! pnsAll[i].open ;
	    break ;
	  }
      }
    }
  
  look = (Htile) halloc (sizeof (struct HtileStruct), 0) ;  
  *look = *look0 ;
  look->magic = &Htile_MAGIC;

  look->db = ac_db ;          /* cosmetic, since we are inside xace */
  look->h = handle ;
  solexaAllInit (look) ;
  look->key = key ;
  look->from = from ;

  look->hideHeader = oldHideHeader ;
  look->hideMusic = oldHideMusic ;
  look->hideScale = oldHideScale ;

  memset (look->geneSearchBuf, 0, 301) ;
  look->map = tmapCreate (look, look->h, htileDraw) ;
  look->map->lnWidth = .3 ;

  if (!htileConvert (look, FALSE, 0) || ! (isOldGraph ||  displayCreate (DtTiling)))
    { 
      handleDestroy (handle) ;
      ac_free (look) ;
      display (key0,from,TREE) ;
      return FALSE ;
    }
    
  htileSolexaConvert (look, FALSE, 0) ;
  if (look->preserve)
    displayPreserve () ;
  graphRegister (DESTROY, htileDestroy) ;

  graphRetitle (name(key)) ;

  graphAssociate (&Htile_MAGIC, look) ; /* attach look to graph */
  tmapGraphAssociate (look->map) ;
  tmapResize () ; /* redraws */

  htileShowSegs(0) ; /* for compiler happiness */

  return TRUE ;
} /* htileDisplay */

/*************************************************************/
/* precompute all gmaps of all regions in the active keyset */
void htilePrecompute (KEYSET ks0)
{
  AC_HANDLE h = ac_new_handle () ;
  KEYSET ks = 0 ;
  int ii, nn = 0 ;
  Htile look = 0 ;
  char *cp ;
  ACEOUT ao = 0 ;

  cp = freeword () ;
  if (cp && !strcmp (cp, "-x") && (cp = freeword ()))
    ao = aceOutCreateToFile (cp, "w", h) ;
  ks = query (ks0, "CLASS Region") ;
  cDNAAlignInit () ;
  if (!ac_db)
    ac_db = ac_open_db (0,0) ;

  look = (Htile) halloc (sizeof (struct HtileStruct), 0) ;  
  look->magic = &Htile_MAGIC;

  look->db = ac_db ;          /* cosmetic, since we are inside xace */
  
  for (ii = nn = 0 ; ks && ii < keySetMax (ks) ; ii++)
    {
      look->key = keySet (ks, ii) ;
      look->h = ac_new_handle () ;
      solexaAllInit (look) ;
      look->map = tmapCreate (look, look->h, 0) ;
      look->smoothing = 0 ;
      look->gaussWidth = 1000 ;
      if (htileConvert (look, TRUE, ao))
	nn++ ;
      ac_free (look->h) ;
    }
  keySetDestroy (ks) ;
  freeOutf ("// precomputed %d regions\n", nn) ;
  ac_free (look) ;
  ac_free (h) ;
}  /* htilePrecompute */

/*************************************************************/
/*************************************************************/
/*************************************************************/
