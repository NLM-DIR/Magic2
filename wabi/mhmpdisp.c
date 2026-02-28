/*  File: mhmpdisp.c
 *  Author: Jean Thierry-Mieg, Michel Potdevin (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@ncbi.nlm.nih.gov
 *
 * Description: Display of the genetic map
 * Exported functions:
       mhmpDisplay
 * HISTORY:

 * Created: July 15 2002 (mieg)
 *-------------------------------------------------------------------
 */
#include "acedb.h"
#include "bump.h"

#include "display.h"
#include "pick.h"
#include "session.h"
#include "systags.h"
#include "tags.h"
#include "classes.h"
#include "sysclass.h"
#include "bs.h"
#include "lex.h"
#include "a.h"
#include "graph.h"
#include "query.h"
#include "cdna.h"

typedef struct MhmStruct {
  magic_t *magic;        /* == MHM_MAGIC */
  AC_HANDLE handle ;
  char title[1024] ;
  char unit[5] ;
  int isIntMap ; /* 1: use the IntMap tag,  2: use the Map tag */
  int type ; /* 1: gene, 2:mrna, 3:product */
  int colour, importBox, glyphBox ;
  int decimate, showNT ;
  BOOL showButtons ;
  int showCaption ;
  KEY glyph ;
  float offset, yoffset, mag ;
  float topMargin, bottomMargin, leftMargin, rightMargin ;
  float xMag, yMag ;
  KEYSET ks ;
  Array mapSet, geneSet ;
  int graphWidth, graphHeight ;
  int maxLen, nActiveMaps ;
  int   activeBox ;
  int	messageBox ;
  char  messageText[128] ;
  unsigned int flag ;
  Array boxIndex ;    /* box -> int (index into look->segs) */
  Array captions ;
  float xCaption, yCaption ;
  int captionBox ;
  Graph graph ;
} *MHM;



typedef struct MhmCaptionStruct {
  KEY glyph ;
  int col, nn ;
  char txt[1024] ;
} MHMCAPTION ;


typedef struct MhmsetStruct {
  KEY map ;
  int min, centre, max, len, nn ;
  float y ;
} MHMSET ;


typedef enum { Tzero, Tgene, Tmrna, GENE } HSTYPE  ;
typedef struct MhgeneStruct {
  KEY gene, map, glyph ;
  int col, bcol ;
  int a1, a2 ;
  int mapIndex, caption ;
  int type ;
  HSTYPE hType ;
  unsigned int flag ;
} MHGENE ;

static void  mhmImport (void) ;
static void mhmDraw (MHM look) ;
static BOOL mhmGlyphBoxDraw (MHM look) ;
static magic_t MHM_MAGIC = "Mhm";
static void  mhmImportTableColumn (void) ;

/************************************************************/
/************************************************************/
/* to be called from the debugger */
static void mhmGeneShow (Array aa)
{
  int i ;
  MHGENE *gg ;

  for (i = 0 ; aa && i < arrayMax(aa) ; i++)
    {
      gg = arrp (aa, i, MHGENE) ;
      printf("%s %d %d\tglyph=%d mapIndex=%d %s\n"
	     ,name(gg->map), gg->a1, gg->a2, gg->glyph, gg->mapIndex, name(gg->gene)
	     );
    }
  if (aa) mhmGeneShow (0) ; /* for compiler happiness */
}

/************************************************************/
/************************************************************/

static MHM currentMhm (char *caller)
{
  void *v = 0 ;
  MHM mhm;

  if (!graphAssFind (&MHM_MAGIC,&v))
    messcrash("%s() could not find Mhm on graph", caller);
  if (!(mhm=v))
    messcrash("%s() received NULL Mhm pointer", caller);
  if (mhm->magic != &MHM_MAGIC)
    messcrash("%s() received non-magic Mhm pointer", caller);
  
  return mhm;
} /* currentMhm */


/************************************************************/
/***************** Registered routines *********************/

static void mhmDestroy (void)
{
  MHM look = currentMhm("mhmDestroy") ;

  if (look)
    {
      look->magic = 0 ;
      graphAssRemove (&MHM_MAGIC) ; /* detach look from graph */
      
      handleDestroy (look->handle) ;
    }

  return;
} /* mhmDestroy */

/*************************************************************/

static void mhmRecalculate(void)
{
  MHM look = currentMhm("mhmRecalculate") ;
  
  if (look)
    mhmDraw (look) ;
  pickRememberDisplaySize ("DtMHM") ;
}

/*************************************/
/*************************************/

static void mhmKbd (int k)
{ 
  switch (k)
    {
    case UP_KEY :
      break ;
    default:
      return ;
    }
}

/*************************************/

static void mhmFollow (MHM look, int box)
{
  int ii ;
  MHGENE *gg ;
 
  if (box > 0 &&
      box < arrayMax(look->boxIndex) &&
      (ii = arr(look->boxIndex, box, int)) &&
      ii >= 0 && ii < arrayMax (look->geneSet) &&
      (gg = arrp (look->geneSet, ii, MHGENE)) &&
      gg->gene)
    display (gg->gene, 0, TREE) ;
} /* mhmFollow */

/*************************************/

static void mhmSelect (MHM look, int box)
{
  int ii ;
  MHGENE *gg = 0 ;
 
  if (look->activeBox > 0 &&
      look->activeBox < arrayMax(look->boxIndex) &&
      (ii = arr(look->boxIndex, look->activeBox, int)) &&
      ii >= 0 && ii < arrayMax (look->geneSet) &&
      (gg = arrp (look->geneSet, ii, MHGENE)))
    graphBoxDraw (look->activeBox, gg->col, gg->bcol) ;
  
  look->activeBox = 0 ;
  if (box > 0 &&
      box < arrayMax(look->boxIndex) &&
       (ii = arr(look->boxIndex, box, int)) &&
      ii >= 0 && ii < arrayMax (look->geneSet) &&
      (gg = arrp (look->geneSet, ii, MHGENE)))
    { 
      look->activeBox = box ;
      graphBoxDraw (look->activeBox, BLUE, YELLOW) ;
    }

  if (look->activeBox && look->messageBox)
    { 
      if (gg->gene) strncpy (look->messageText, name (gg->gene), 100) ;
      strcat (look->messageText, messprintf(" %d  %d", gg->a1, gg->a2)) ;
      
      graphBoxDraw (look->messageBox, BLACK, CYAN) ;
    }
}

/*************************************/

void mhmCaptionDrag (float *x, float *y, BOOL isDone)
{
  if (isDone)
    { 
      MHM look = currentMhm("mhmCaptionDrag") ;
      
      look->xCaption = (*x > 1 ? *x : 1) ;
      look->yCaption = (*y > 1 ? *y : 1) + 1 ; 
      mhmDraw (look) ;
    }
  
  return;
} /* mhmCaptionDrag */

/*************************************/

static void mhmPick (int box,  double x , double y) 
{
  MHM look = currentMhm("mhmPick") ;

  if (!box)
    return ;
  else if (box == look->captionBox)
    graphBoxDrag (box, mhmCaptionDrag) ;
  else if (box < arrayMax (look->boxIndex) && arr(look->boxIndex,box,int)) /* a gene */
    {
      if (box == look->activeBox) /* a second hit - follow it */
	mhmFollow (look, box) ;
      else
	mhmSelect (look, box) ;
    }

  return;
} /* mhmPick */

/*************************************************************/
/************************************************************/

static void mhmToggleButtons (void)
{
  MHM look = currentMhm("mhmToggleButtons") ;

  look->showButtons = !look->showButtons ;
  mhmDraw (look) ;
} /* mhmToggleButtons */

/************************************************************/

static void mhmToggleCaption (void)
{
  MHM look = currentMhm("mhmToggleCaption") ;

  switch (look->showCaption)
    {
    case 2: look->showCaption = 1 ; break ;
    case 1: look->showCaption = 0 ; break ;
    case 0: look->showCaption = 2 ; break ;
    }

  mhmDraw (look) ;
} /* mhmToggleCaption */

/************************************************************/

static void mhmToggleEditCaption (void)
{
  MHM look = currentMhm("mhmToggleeditCaption") ;

  switch (look->showCaption)
    {
    case 2: look->showCaption = 3 ; break ;
    case 3: look->showCaption = 2 ; break ;
    }

  mhmDraw (look) ;
} /* mhmToggleEditCaption */

/************************************************************/

static FREEOPT mhmLegendMenu [] =
{
  { 3, "Legend"},
  { 0, "Hide the legend"},
  { 1, "Show the legend"},
  { 3, "Edit the legend and compress the data"}
} ;

  /* the range 10 000 to 20 000 is forbidden */

static BOOL mhmLegendSelector (KEY k, int box)
{
  MHM look = currentMhm("mhmColourSelector") ;

  look->showCaption = k ;
  mhmDraw (look) ;

  return TRUE ; 
} /* mhmGlyphSelector */

/************************************************************/

static FREEOPT mhmGlyphMenu [] =
{
  { 22, "Glyphs"},
  /* the range 10 000 to 20 000 is forbidden */
  { 21009, "Rectangle Over" }, /* > 20000 clustered in geneOrder and merged */
  { 22009, "Rectangle Under" },
  { 01, "Line Over"},
  { 00, "Line On"},
  { 02, "Line Under"},
  { 11, "V Over"},
  { 10, "V On"},
  { 20, "^ On"},
  { 22, "^ Under"},
  { 41, "Circle-above Over"},
  { 40, "Circle-above On"},
  { 80, "Circle-under On"},
  { 82, "Circle-under Under"},

  { 03, "  Line Over/Line under"},
  { 13, "     V Over/Line under"},
  { 43, "Circle Over/Line under"},

  { 23, "  Line Over/^ under"},
  { 33, "     V Over/^ Under"},
  { 63, "Circle Over/^ Under"},

  { 83, "  Line Over/Circle Under"},
  { 93, "     V Over/Circle Under"},
  {123, "Circle Over/Circle Under"},
} ; /* mhmGlyphMenu */

/************************************************************/

static void mhmGlyphDraw (KEY glyph, float x, float y, float dx)
{
  float y1, y2 ;

  if (glyph / 10 == 2100)
    graphFillRectangle (x, y - .8, x + dx, y ) ;
  else if (glyph / 10 == 2200)
    graphFillRectangle (x, y, x + dx, y + .8) ;
  else
    {
      switch (glyph % 10)
	{
	case 1: y1 = y - .8 ; y2 = y ; graphLine (x, y1, x, y2) ; break ; 
	case 0: y1 = y - .4 ; y2 = y + .4 ; graphLine (x, y1, x, y2) ; break ; 
	case 3: y1 = y - .8 ; y2 = y + .8 ; graphLine (x, y1, x, y2) ; break ; 
	case 2: y1 = y ; y2 = y + .8 ; graphLine (x, y1, x, y2) ; break ; 
	default: y1 = y2 = y ; break ;
	}
      if ((glyph / 10) & 0x1) /* V over */
	{
	  graphLine (x - .3, y1 - .4, x, y1) ;
	  graphLine (x + .3, y1 - .4, x, y1) ; 
	}
      if ((glyph / 10) & 0x2) /* ^ Under */
	{
	  graphLine (x - .3, y2 + .4, x, y2) ;
	  graphLine (x + .3, y2 + .4, x, y2) ; 
	}
      if ((glyph / 10) & 0x4) /* circle over */
	graphCircle (x, y1 - .2, .4) ;
      if ((glyph / 10) & 0x8) /* circle under */
	graphCircle (x, y2 + .2, .4) ;
    }

  return ;
} /* mhmGlyphDraw */

/************************************************************/

static BOOL mhmGlyphSelector (KEY glyph, int box)
{
  MHM look = currentMhm("mhmColourSelector") ;

  look->glyph = glyph ;
  graphBoxClear (look->glyphBox) ;
  mhmGlyphBoxDraw (look) ;

  return TRUE ;
} /* mhmGlyphSelector */

/************************************************************/

#define graphBoxBox(_box) { \
	       float _x1, _y1, _x2, _y2 ; \
	       graphBoxDim (_box, &_x1, &_y1, &_x2, &_y2) ; \
	       graphRectangle (_x1 - 0.4, _y1 - 0.1, _x2 + 0.4, _y2 + 0.1) ; \
		}	   

static BOOL mhmGlyphBoxDraw (MHM look)
{
  float x = look->graphWidth - 20, y = 1.4 ;

  look->glyphBox = graphBoxStart () ;
  graphText ("Glyph", x - 9, y) ;
  mhmGlyphDraw (look->glyph, x - 1, y + .7, 1.0) ;
  graphText ("..", x, y - .2) ;
  graphBoxEnd () ;
  graphBoxBox (look->glyphBox) ;
  graphBoxDraw (look->glyphBox, BLACK, WHITE) ; 
  graphBoxFreeMenu (look->glyphBox, (FreeMenuFunction) mhmGlyphSelector, mhmGlyphMenu) ;

  return TRUE ;
} /* mhmGlyphBoxDraw */

/************************************************************/
/************************************************************/

static BOOL mhmMapSetCreate (MHM look, KEYSET ksAll)
{
  int ii ;
  float a1, a2 ;
  KEYSET ks1 ;
  MHMSET *mm ;
  OBJ Map ;

  arrayDestroy (look->mapSet) ;
  look->mapSet = arrayHandleCreate (32, MHMSET, look->handle) ;
 
  ks1 = query (ksAll, "{>IntMap} $| {>Map}") ;
  arraySort (ks1, keySetAlphaOrder) ;
  
  ii = keySetMax (ks1) ;
  while (ii--)
    {
      mm = arrayp (look->mapSet, ii, MHMSET) ;
      mm->map = keySet (ks1, ii) ;
      mm->min = 1 ;
      mm->centre = 0 ;
      mm->max =   1 ;
      if (keyFindTag (mm->map, str2tag ("Extent")) &&
	  (Map = bsCreate (mm->map)))
	{
	  if (bsGetData (Map, str2tag ("Extent"), _Float, &a1) &&
	      bsGetData (Map, _bsRight, _Float, &a2))
	    { mm->min = a1 ; mm->max = a2 ; }
	  if (bsGetData (Map, str2tag ("Centre"), _Float, &a1))
	    mm->centre = a1 ;
	  bsDestroy (Map) ;
	}
      mm->len = mm->max - mm->min + 1 ;
      if (keyFindTag (mm->map, str2tag ("Centromere_telomeres")) &&
	  (Map = bsCreate (mm->map)))
	{
	  Array aa = arrayCreate (3, BSunit) ;
	  BSunit *uu ;

	  if (bsGetArray (Map, str2tag ("Centromere_telomeres"), aa, 4) &&
	      arrayMax (aa) >= 4
	      )
	    {
	      uu = arrp (aa, 0, BSunit) ;
	      mm->min = uu[0].f ; mm->centre = (uu[1].f + uu[2].f)/2 ; mm->max = uu[3].f ; 
	    }
	  arrayDestroy (aa) ;
	  bsDestroy (Map) ;
	}
    }
  keySetDestroy (ks1) ;

  return TRUE ;
}

/************************************************************/

static int mhgeneOrder (const void *va, const void *vb)
{
  const MHGENE *ma = (const MHGENE*)va,  *mb = (const MHGENE*)vb ;

  if (ma->mapIndex != mb->mapIndex) return ma->mapIndex < mb->mapIndex ? -1 : 1 ;
  if ((int)mb->glyph > ma->glyph + 10000) return 1 ; /* careful , glyph is unsigned */
  if ((int)ma->glyph > mb->glyph + 10000) return -1 ;
  if (ma->glyph > 20000 && ma->glyph != mb->glyph) 
    return ma->glyph < mb->glyph ? -1 : 1 ;
  if (ma->glyph > 20000 && ma->col != mb->col) 
    return ma->col < mb->col ? 1 : -1 ;
  if (ma->a1 != mb->a1) return ma->a1 - mb->a1 ;
  if (ma->a2 != mb->a2) return ma->a2 - mb->a2 ;
  return ma->gene < mb->gene ? -1 : 1 ;
}

/**************/ 

static BOOL mhmGeneSetCreate (MHM look, KEYSET ks0, char *title)
{
  Array aa ;
  int jj, i, ii, i0, index ;
  BSunit *uu ;
  Array units = 0 ;
  MHGENE *gg ;
  MHMSET *mm ;
  MHMCAPTION *gc ;
  OBJ Gene = 0 ;
  char *cp ;
  KEY gene ;
  static char buf[1024] ;
  KEYSET ksNew = look->ks ? keySetMINUS (ks0, look->ks) : ks0 ;

  if (!look->captions)
    {
      look->captions = arrayHandleCreate (32, MHMCAPTION, look->handle) ;
      gc = arrayp (look->captions, 0, MHMCAPTION) ; /* avoid 0 */
      memset (buf, 0, 1024) ;
    }
  if (! messPrompt ("Caption", buf, "w"))
    memset (buf, 0, 1024) ;
  else
    {
      cp = freepos () ;
      strncpy (buf, cp, 1023) ;
    }
  if (title)
    strncpy (title, buf, 1023) ;

  for (ii = index = 0 ; !index && ii < arrayMax (look->captions) ; ii++)
    {
      gc = arrp (look->captions, ii, MHMCAPTION) ;
      if (gc->col == look->colour && gc->glyph == look->glyph)
	index = ii ;
    }
  if (!index)
    {
      index = ii ;
      gc = arrayp (look->captions, index, MHMCAPTION) ;
      gc->col = look->colour ;
      gc->glyph = look->glyph ;
    }
  strncpy (gc->txt, buf, 1023) ;
  
  if (!look->geneSet)
    look->geneSet = arrayHandleCreate (32, MHGENE, look->handle) ;
  aa = look->geneSet ; /* alias to aa for convenience */
  /* recalculate the mapIndex since we changed the mapSet */
  ii = i0 = arrayMax (aa) ;
  while (ii--)
    {
      gg = arrp (aa, ii, MHGENE) ;
      if (keySetFind (ks0, gg->gene, 0))
	{
	  gg->col = look->colour ;
	  gg->glyph = look->glyph ;
	  gg->caption = index ;
	}
      for (jj = 0 ; jj < arrayMax (look->mapSet) ; jj++)
	{
	  mm =  arrayp (look->mapSet, jj, MHMSET) ;
	  if (mm->map == gg->map)
	    {
	      gg->mapIndex = jj ;
	      if (gg->a1 > mm->max) mm->max = gg->a1 ;
	      break ;
	    }
	}	  
    }
  
  /* include the new genes */
  ii = keySetMax (ksNew) ;

  while (bsDestroy (Gene), ii--)
    {
      int a1, a2 ;

      gene = keySet (ksNew, ii) ;
      units = arrayReCreate (units, 300, BSunit) ;

      if (!look->isIntMap)
	{
	  if (bsIsTagInClass (class (gene), _IntMap))
	    {
	      look->isIntMap = 1 ;
	      strcpy (look->unit, "b") ;
	    }
	  else
	    {
	      look->isIntMap = 2 ;
	      strcpy (look->unit, "") ;
	    }
	}
      if ((Gene = bsCreate (gene)))
	{
	  a1 = a2 = 0 ;

	  if (look->isIntMap == 1)
	    bsGetArray (Gene, _IntMap, units, 3) ;
	  else  if (look->isIntMap == 2)
	    {
	      float u1, u2 ;
	      if (bsGetKey (Gene, _Map, &gg->map) &&
		  bsPushObj (Gene))
		{
		  if (bsGetData (Gene, _Position, _Int, &u1))
		    a2 = a1 = u1 ;
		  else if (bsGetData (Gene, _Position, _Float, &u1))
		    bsGetData (Gene, _Right, _Float, &u2) ;
		}
	      a1 = u1 ; a2 = u2 ;
	      uu = arrayp (units, 2, BSunit) - 3 ;
	      uu[0].k = gg->map ;
	      uu[1].i = a1 ;
	      uu[2].i = a2 ;
	    }
	}
      
      for (i = 0 ; i < arrayMax (units) ; i += 3)
	{
	  uu = arrayp (units, i, BSunit) ;
	  
	  gg = arrayp (aa, i0++, MHGENE) ;
	  gg->gene = gene ;
	  gg->col = look->colour ;
	  gg->glyph = look->glyph ;
	  gg->caption = index ;
	  gg->map = uu[0].k ;
	  a1 = uu[1].i ; a2 = uu[2].i ;
	  if (a1 < a2) { gg->a1 = a1 ; gg->a2 = a2 ; }
	  else { gg->a1 = a2 ; gg->a2 = a1 ; }
	  for (jj = 0 ; jj < arrayMax (look->mapSet) ; jj++)
	    {
	      mm =  arrayp (look->mapSet, jj, MHMSET) ;
	      if (mm->map == gg->map)
		{
		  gg->mapIndex = jj ;
		  if (gg->a1 > mm->max) mm->max = gg->a1 ;
		  if (gg->a2 > mm->max) mm->max = gg->a2 ;
		  if (gg->a1 < mm->min) mm->min = gg->a1 ;
		  if (gg->a2 < mm->min) mm->min = gg->a2 ;
		  break ;
		}
	    }	  
	}
    }
  arraySort (aa, mhgeneOrder) ; 
  
  for (jj = 0 ; jj < arrayMax (look->mapSet) ; jj++)
    {
      mm =  arrayp (look->mapSet, jj, MHMSET) ;
      mm->len = mm->max - mm->min + 1 ;
    }

  /* count the members of each caption category */
  for (ii = 0, gc = arrp (look->captions, ii, MHMCAPTION) ; ii < arrayMax (look->captions) ; gc++, ii++)
    gc->nn = 0 ;

  ii = arrayMax (aa) ; 
  gg = arrp (aa, 0, MHGENE)  - 1 ;
  gc = arrp (look->captions, 0, MHMCAPTION) ; 
  while (gg++, ii--)
    ((gc + gg->caption)->nn)++ ;

  if (look->ks) keySetDestroy (ksNew) ;
  arrayDestroy (units) ;
  return TRUE ;
}

/************************************************************/

static BOOL mhmConvert (MHM look, KEYSET ks0, char *title)
{
  KEYSET ksAll = 0, ks1 = query (ks0, "(IntMap AND NEXT AND NEXT) || (Map && NEXT AND NEXT)") ;

  /* restrict to mapable objects */
  if (!keySetMax (ks1))
    {
      messout ("Sorry, there is no object in this display with IntMap or Map coordinates") ;
      keySetDestroy (ks1) ;
      return FALSE ;      
    } 

  if (look->ks)
    ksAll= keySetOR (look->ks, ks1) ;
  else
    ksAll = ks1 ;

  /* create the map set */
  mhmMapSetCreate (look, ksAll) ;
  
  /* create the mapped objects */
  mhmGeneSetCreate (look, ks1, title) ;

  /* remember the total set current genes */
  keySetDestroy (look->ks) ;
  look->ks = arrayHandleCopy (ksAll, look->handle) ;

  if (ksAll != ks1)
    keySetDestroy (ksAll) ;  
  keySetDestroy (ks1) ;

  return TRUE ;
}

/************************************************************/

static MENUOPT mhmMenu[] = {
  { graphDestroy, "Quit" },
  { help, "Help" },
  { graphPrint, "Print" },
  { mhmToggleButtons, "Show/Hide Buttons" },
  { mhmImportTableColumn , "Import from a file" },
  { 0, 0 }
} ;

/************************************************************/

static void mhmDrawMaps (MHM look)
{
  int ii, jj, maxLen = 1 ;
  MHMSET *mm ;
  MHGENE *gg ;
  float y, x0 = look->leftMargin ;
  char *cp ;

  /* count active map members */
  for (ii = 0 ; ii < arrayMax (look->mapSet) ; ii++)
    {
      mm = arrayp (look->mapSet, ii, MHMSET) ;
      mm->nn = 0 ;
    }
  mm = arrp (look->mapSet, 0, MHMSET) ;
  for (ii = 0, gg = arrp (look->geneSet, 0, MHGENE) ; ii < arrayMax (look->geneSet) ; ii++, gg++)
    if (gg->col != WHITE) 
      ((mm + gg->mapIndex)->nn)++ ;

  look->nActiveMaps = 0 ; /* avoid zero */
  for (ii = 0 ; ii < arrayMax (look->mapSet) ; ii++)
    {
      mm = arrayp (look->mapSet, ii, MHMSET) ;
      if (! look->showNT)
	{
	  cp = name(mm->map) ;
	  if (strstr (cp, "NT"))
	    {
	      mm->nn = 0 ;
	      continue ;
	    }
	}
      if (mm->nn)
	{
	  if (mm->len > maxLen)
	    maxLen = mm->len ;
	  look->nActiveMaps++ ;
	}
    }
  if (!look->nActiveMaps)
    { 
      look->nActiveMaps = 1 ; /* avoid zero */
      (mm->nn)++ ;
      maxLen = mm->len ;
    }
  look->maxLen = maxLen ;
  look->xMag = (look->graphWidth - look->leftMargin - look->rightMargin)/maxLen ;
  look->yMag = (look->graphHeight - look->topMargin - look->bottomMargin)/look->nActiveMaps ;

  for (ii = jj = 0 ; ii < arrayMax (look->mapSet) ; ii++)
    {
      mm =  arrayp (look->mapSet, ii, MHMSET) ;
      if (!mm->nn)
	continue ;
      y = look->topMargin + look->yMag * (jj++) ;
      mm->y = y ;
      cp = name(mm->map) ;
      if (0 && !strncasecmp(cp, "CHROMOSOME_", strlen ( "CHROMOSOME_")))
	cp += strlen ( "CHROMOSOME_") ;
      graphText (cp, 1, y - .5) ;
      graphLine (x0, y, x0 + look->xMag * mm->len, y) ;
      if (mm->centre)
	graphFillRectangle (x0 + look->xMag * mm->centre - .5, y - .5,
			    x0 + look->xMag * mm->centre + .5, y + .5) ;
    }
  graphText (" ", look->graphWidth - 1, look->graphHeight - 1) ;
} /* mhmDrawMaps */

/************************************************************/

static void mhmDrawScale (MHM look)
{
  int i, iMax = 0, step, Step, nStep, exp, expScale, trueScale ;
  float x, oldx, DX, y, x0 = look->leftMargin ;
  char unit, *cp ;

  y = 1 + look->topMargin + look->yMag * look->nActiveMaps ;

  /* count the number of needed small tick in the page */
  DX = look->graphWidth -  look->leftMargin - look->rightMargin ;
  nStep = utMainRoundPart (DX/2) ; /* we like a minor tick every 2 char */
  step = utMainRoundPart (DX / (nStep * look->xMag)) ;
  if (step < 1) step = 1 ;

  /* draw the scale line */
  for (i = 0 ; ; i++)
    {
      x = x0 + i * step * look->xMag ;
      if (x > look->graphWidth - 1)
	{
	  iMax = i ;
	  x -= step * look->xMag ;
	  graphLine (x0, y, x, y) ;
	  break ;
	}
    }
  /* draw the small ticks */
  for (i = 0 ; i < iMax ; i++)
    {
      x = x0 + i * step * look->xMag ;
      graphLine (x, y, x, y - .4) ;
    }
  /* draw the big ticks, with numbers */
  Step = step ; exp = 0 ; trueScale = 1 ;
  while ((Step % 10) == 0) { exp++ ; Step /= 10 ; trueScale *= 10 ; }
  switch (exp % 3)
    {
    case 2: expScale = 100 ; break ;
    case 1: expScale = 10 ; break ;
    case 0:
    default: expScale = 1 ; break ;
    }
  switch (exp/3)
    {
    case 1: unit = 'k' ; break ;
    case 2: unit = 'M' ; break ;
    case 3: unit = 'G' ; break ;
    case 4: unit = 'T' ; break ;
    case 0:
    default: unit = ' ' ; break ;
    }
  switch (Step % 10)
    {
    case 1: 
    case 2: Step = 5 * step ; break ;
    case 5: Step = 2 * step ; break ;
    }
  for (i = 0, oldx = -100.0 ; i * Step < iMax * step ; i++)
    {
      x = x0 + i * Step * look->xMag ;
      if (((i*Step)/trueScale) % 10 == 5)
	graphLine (x, y, x, y - .75) ;
      else
	graphLine (x, y, x, y - 1.0) ;
      if (x > oldx + 10)
	{
	  if (i && (x + x - oldx > look->graphWidth))
	    {
	      cp = messprintf ("%d%c%s", expScale * i * (Step/trueScale), unit, look->unit) ;
	      graphText (cp, x - strlen (cp) + 2, y + .3) ;
	    }
	  else
	    {
	      cp = messprintf ("%d", expScale * i * (Step/trueScale)) ;
	      graphText (cp, x - (strlen (cp)/2.0) , y + .3) ;
	    }
	  oldx = x ;
	}
    }
  graphText (" ", look->graphWidth - 1, look->graphHeight - 1) ;
} /* mhmDrawScale */

/************************************************************/

static Array mhmCondense (Array aa, float xMag, int decimate)
{
  Array bb ;
  int ii, jj ;
  MHGENE *gg, *gg2 ;

  if (!aa || !arrayMax (aa))
    return 0 ;
  bb = arrayCreate (arrayMax (aa), MHGENE) ;
  for (ii = jj = 0, gg = arrp (aa, 0, MHGENE), gg2 = arrp (bb, 0, MHGENE) ; 
       ii < arrayMax (aa) ; ii++, gg++)
    {
      if (gg->col == WHITE) continue ;
      if (decimate > 0 && (ii % decimate))
	continue ;
      if (gg->glyph < 20000 || !jj || 
	  gg->glyph != gg2->glyph ||
	  gg->map != gg2->map ||
	  (gg->a1 - gg2->a2) * xMag > .2)	  
	{
	  gg2 = arrayp (bb, jj++, MHGENE) ;
	  *gg2 = *gg ;
	}
      else
	gg2->a2 = gg->a2 ;
    }
  return bb ;
} /* mhmCondense */

/**************/ 

static void mhmDrawGenes (MHM look)
{
  int ii, box ;
  MHGENE *gg ;
  MHMSET *mm ;
  float y, x0 = look->leftMargin ;
  Array geneSet = mhmCondense (look->geneSet, look->xMag, look->decimate) ;
  
  for (ii = 0 ; ii < arrayMax (geneSet) ; ii++)
    {
      gg = arrayp (geneSet, ii, MHGENE) ;
      if (!gg->glyph) continue ;
      mm =  arrayp (look->mapSet, gg->mapIndex, MHMSET) ;
      if (! look->showNT)
	{
	  const char *cp = name(mm->map) ;
	  if (strstr (cp, "NT"))
	    continue ;
	}
      y = mm->y ;
      
      array (look->boxIndex, box = graphBoxStart(), int) = ii ;
      mhmGlyphDraw (gg->glyph, x0 + look->xMag * (gg->a1 - mm->min) - .1, y
		    , look->xMag * (gg->a2 - gg->a1)) ;
      graphBoxEnd () ;
      gg->bcol = TRANSPARENT ;
      graphBoxDraw (box, gg->col, gg->bcol) ;
    }
  arrayDestroy (geneSet) ;
} /* mhmDrawGenes */


/*** polygon
  static Array xy = 0 ;
  
  if (!xy)
    { xy = arrayCreate (10, float) ;
      arrayMax(xy) = 10 ;
    }
      if (0)
	{
	  arr (xy, 0, float) = u1 ; arr (xy, 1, float) = 4 ;
	  arr (xy, 2, float) = u2 - 1.5 ; arr (xy, 3, float) = 4 ;
	  arr (xy, 4, float) = u2  ; arr (xy, 5, float) = 5 ;
	  arr (xy, 6, float) = u2 - 1.5 ; arr (xy, 7, float) = 6 ;
	  arr (xy, 8, float) = u1 ; arr (xy, 9, float) = 6 ;
	  graphPolygon (xy) ;
	}
*****/

/************************************************************/
/************************************************************/

static BOOL mhmColourSelector (KEY col, int box)
{
  MHM look = currentMhm("mhmColourSelector") ;
  
  look->colour = col ;
  graphBoxDraw (look->importBox, BLACK, col) ;
  return TRUE ;
}

/************************************************************/

static void  mhmImport (void)
{
  KEYSET ks = 0, ks1 = 0 ; 
  char title[1024] ;
  MHM look = currentMhm("mhmImport") ;

  if (!keySetActive (&ks, 0) ||
      !keySetMax (ks))
    {
      messout ("First select a keyset containing objects having the IntMap or the Map tag") ;
      goto abort ;
    }
  ks1 = query (ks, "IntMap AND NEXT") ;
  if (!keySetMax (ks1))
    {
      messout ("First select a keyset containing objects having the IntMap or the Map tag") ;
      goto abort ;
    }
  mhmConvert (look, ks1, title) ;
  mhmDraw (look) ;
  keySetNewDisplay (keySetCopy (ks1), title) ;

 abort:
  ks = 0 ; /* do not destroy */
  keySetDestroy (ks1) ;
}

/************************************************************/

static void  mhmImportTableColumn (void)
{
  AC_HANDLE h = ac_new_handle () ;
  KEYSET ks = 0 ;
  int i, column = 0, myClasse = 0, nn = 0, nnr = 0, nnm = 0 ;
  KEY key = 0;
  ACEIN ai = 0 ;
  char cutter, *word = 0 ;
  static char firstPass = 0, fileName[FIL_BUFFER_SIZE] , dirName[DIR_BUFFER_SIZE], suffixName[12], colBuf[32], classBuf[32], title[1024] ;
  MHM look = currentMhm("mhmImportTableColumn") ;

  if (firstPass == 0)
    {
      firstPass = 'a' ;
      memset (fileName, 0, sizeof(fileName)) ;
      memset (colBuf, 0, sizeof(colBuf)) ;
      memset (classBuf, 0, sizeof(classBuf)) ;
      memset (dirName, 0, sizeof(dirName)) ;
      memset (suffixName, 0, sizeof(suffixName)) ;
    }

  /* select a file, column and class name */
  ai = aceInCreateFromChooser ("Please select a tab delimited file containing a column with known mapped object"
			       , dirName, fileName, suffixName, "r", h
			       ) ;
  if (! ai)
    goto abort ;

  if (!messPrompt ("Specify the column containing known mapped objects", colBuf, "iz"))
     goto abort ;
  if (! freeint (&column) || column < 1)
    goto abort ; ;

  if (!messPrompt ("Specify the class (Gene, mRNA ...) of the mapped objects", classBuf, "wz"))
    goto abort ;
  word = freeword () ;
  if (! lexword2key (word, &key, _VMainClasses))
    {
      messout ("Unknown class, sorry") ;
      goto abort ;
    }
  myClasse = KEYKEY(key) ;

  ks = keySetHandleCreate (h) ; nn = 0 ;
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      for (i = 0 ; i < column - 1 ; i++)
	{
	  word = aceInWordCut (ai, "\t", &cutter) ;
	}
       word = aceInWordCut (ai, "\t", &cutter) ;
       if (! word)
	 continue ;
       
       if (! lexword2key (word, &key, myClasse))
	 { nnr++ ; continue ; }
       if (! keyFindTag (key, _IntMap))
	 { nnm++ ; continue ; }

       keySet (ks, nn++) = key ;
    }

  messout ("%d mapped %s, %d names not recognized, %d not mapped", nn, classBuf, nnr, nnm) ;

  keySetSort (ks) ; keySetCompress (ks) ;
  mhmConvert (look, ks, title) ;
  mhmDraw (look) ;
  keySetNewDisplay (keySetCopy (ks), title) ;

 abort:
  ac_free (h) ;
}

/************************************************************/

static void  mhmButtons (MHM look)
{ 
  int box ;
  float line = .5 ;

  if (look->showCaption != 3)
    graphText (look->title, 3, 1) ;
  if (look->showButtons &&
      !isGifDisplay)
    {
      if (look->showCaption < 2)
	{
	  graphButton ("Quit", graphDestroy, look->graphWidth - 7, line) ;
	/* mmhpColorSelector */
	  box = look->importBox = graphButton("Import..", mhmImport, look->graphWidth - 17., line) ;
	  graphBoxDraw (box, BLACK, look->colour) ;
	  graphBoxFreeMenu (box, (FreeMenuFunction) mhmColourSelector, &graphColors[0]) ;
	  mhmGlyphBoxDraw (look) ;
	  box = graphButton("Legend..", mhmToggleCaption, look->graphWidth - 39., line) ;
	  graphBoxFreeMenu (box, (FreeMenuFunction) mhmLegendSelector, mhmLegendMenu) ;
	}
      else if (look->showCaption == 2)
	{
	  graphButton ("Edit",  mhmToggleEditCaption, look->graphWidth - 30, line) ;
	  graphButton ("Back to Map-display",  mhmToggleCaption, look->graphWidth - 23, line) ;
	}
      else if (look->showCaption == 3)
	{
	  graphButton ("Apply",  mhmToggleEditCaption, look->graphWidth - 23, line) ;
	}
    }
}

/************************************************************/

static void mhmEditCaption (MHM look)
{
  MHMCAPTION *gc ;
  int ii ;
  float x, y, dx ;

  x = 3 ; y = 11 ;
  dx = look->graphWidth - x - 4 ;
    
  graphTextScrollEditor ("Title ", look->title, 1024, dx, 1, 2, 0) ;
  graphTextEditor ("Unit ", look->unit, 4, 1, 3.5, 0) ;
  graphText ("appears at the end of the x axis, automatically scaled to k, M, G, Tunit", 20, 3.5) ;
  graphIntEditor ("Decimate ", &look->decimate, 1, 5.2, 0) ;
  graphText ("if set to 10, only one gene in 10 is drawn", 20, 5.2) ;
  graphIntEditor ("Show NT", &look->showNT, 1, 6.9, 0) ;
  graphText ("if set to 0, just show chrom 1, 2 .. X, Y ; if 1 also show NT fragments", 20, 6.9) ;

  for (ii = 1, gc = arrp (look->captions, ii, MHMCAPTION) ;
       ii < arrayMax (look->captions) ; gc++, ii++)
    {
      if (gc->col == WHITE || !gc->nn) continue ;
      graphColor (gc->col) ;
      mhmGlyphDraw (gc->glyph, x, y, 1.0) ;
      graphColor (BLACK) ;
      graphTextScrollEditor (messprintf ("%5d", gc->nn), gc->txt, 1024, dx, x+2, y - .5, 0) ;
      y += 3.0 ;
    }
  graphText (" ", look->graphWidth, y+4) ;
} /* mhmEditCaption */

/************************************************************/

static void mhmDrawCaption (MHM look)
{
  MHMCAPTION *gc ;
  int box, ii ;
  float x, y ;

  if (look->xCaption == 0)
    {
      look->xCaption = look->graphWidth - 15 ;
      look->yCaption = look->graphHeight - 2.8 * arrayMax (look->captions) - 4 ;
    }
  x = look->xCaption ;
  y = look->yCaption ;
  if (look->showCaption == 2)
    {  x = 4 ; y = 6 ; } 
  if (y < 3) y = 3 ;
  if (x < 0) x = 0 ;
  
  box = graphBoxStart () ;
  if (look->showCaption == 1)
    look->captionBox = box ;

  for (ii = 1, gc = arrp (look->captions, ii, MHMCAPTION) ;
       ii < arrayMax (look->captions) ; gc++, ii++)
    {
      if (gc->col == WHITE || !gc->nn) continue ;
      graphColor (gc->col) ;
      mhmGlyphDraw (gc->glyph, x, y, 1.0) ;
      graphColor (BLACK) ;
      graphText (messprintf("%d %s", gc->nn, gc->txt), x + 2, y - .5) ;
      y += 2.8 ;
    }
  graphBoxEnd () ;
  
  graphText (" ", look->graphWidth, y+4) ;
} /* mhmDrawCaption */

/************************************************************/

static void mhmDraw (MHM look)
{
  graphFitBounds (&look->graphWidth,&look->graphHeight) ;
  graphRetitle (look->title) ;
  graphClear () ;
  
  mhmGeneShow (0) ; /* for compiler happiness */
  graphMenu (mhmMenu) ;
  
  arrayMax(look->boxIndex) = 0 ;
  mhmButtons (look) ;

  switch (look->showCaption)
    {
    case 0:
      look->captionBox = 0 ;
      mhmDrawMaps (look) ;
      mhmDrawGenes (look) ;
      mhmDrawScale (look) ;
      break ;
    case 1: 
      mhmDrawMaps (look) ;
      mhmDrawGenes (look) ;
      mhmDrawScale (look) ;
      mhmDrawCaption (look) ;
      break ;
    case 2:
      mhmDrawCaption (look) ;
      break ;
    case 3:
      mhmEditCaption (look) ;
      break ;
    }

  graphRedraw () ;
}

/************************************************************/
/*************************************************************/

BOOL mhmpDisplay (KEYSET ks, int glyph, KEY colour)
{
  MHM look ;
  AC_HANDLE handle = 0 ;

  if (!ks)
    return FALSE ;

  cDNAAlignInit () ;

  handle = handleCreate () ;
  look = (MHM) halloc (sizeof (struct MhmStruct), handle) ; 
  look->magic = &MHM_MAGIC;
  look->handle = handle ;
  look->boxIndex = arrayHandleCreate (64,int, handle) ;
  look->activeBox = 0 ;

  look->topMargin = 6 ;
  look->bottomMargin = 3 ; /* space for the scale */
  look->leftMargin = 12 ;
  look->rightMargin = 2 ;
  look->colour = colour ? colour : RED ;
  look->glyph = glyph ? glyph : 11 ;
  look->showButtons = TRUE ;
  look->showCaption = 1 ;
  strcpy (look->title, "Chromosome map") ;

  if (!mhmConvert (look, ks, 0))
    { 
      handleDestroy (handle) ;
      return FALSE ;
    }
  
  if (!displayCreate ("DtMHMP"))
    return FALSE ;

  graphRetitle (look->title) ;
  graphRegister (RESIZE,mhmRecalculate) ;
  graphRegister (DESTROY, mhmDestroy) ;
  graphRegister (KEYBOARD, mhmKbd) ;
  graphRegister (PICK, mhmPick) ;

  look->graph = graphActive() ;

  graphAssociate (&MHM_MAGIC, look) ; /* attach look to graph */

  mhmDraw (look) ;

  return TRUE ;
} /* mhmpDisplay */

/*************************************************************/

BOOL mhmpDrawGif (KEYSET ks)
{
  BOOL ok = FALSE, noNT = FALSE ;
  KEY colour = 0, nnc ;
  int nn, glyph = 0 ;
  char *word ;
  KEYSET ks1 = 0 ;

  while (freestep ('-'))	/* options */
    if ((word = freeword()))
      {
	if (!strcmp (word, "noNT"))
	  noNT = TRUE ;
	if (!strcmp (word, "glyph") && 
	    freeint (&nn))
	  glyph = nn ;
	if (!strcmp (word, "colour"))
	  {
	    word = freeword() ;
	    if (word)
	      {
		if (*word == '-')
		  freeback () ;
		else if (lexword2key (word, &nnc, 0) &&
			  nnc > _WHITE &&
			  nnc < _WHITE + NUM_TRUECOLORS) 
		  colour= WHITE + nnc - _WHITE ;
	      }
	  }
      }
	
  if (noNT) ks1 = query (ks, "! IntMap = *N* && ! Map = *N*") ;
  else ks1 = query (ks, "IntMap || Map") ;
  ok = mhmpDisplay (ks1, glyph, colour) ; 
  keySetDestroy (ks1) ;
  return ok ;
}


/*************************************************************/
/*************************************************************/
/*************************************************************/
