/*  File: fmaptrace.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
   Display virtual subsequence a la maniere des intervals in gmap
   The virtual parent is a long green box.
   No dna recursion is performed.
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 22 10:47 1998 (fw)
 * Created: Thu Dec  9 00:01:50 1993 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */

/* #define ARRAY_CHECK */

#include "fmap_.h"

#include "display.h"
#include "bs.h"
#include "a.h"
#include "dna.h"
#include "cdna.h"
#include "sysclass.h"
#include "classes.h"
#include "systags.h"
#include "bump.h"
#include "tags.h"
#include "dnaalign.h"
#include "query.h"
#include "dotter.h"
#include "acembly.h" 
#include "pick.h"
#include "session.h"

/************************************************************/


void fMapTraceSuppress (KEY key, KEY tag, BOOL keep) ; /* also used by trace.c */

/************************************************************/

static Array findAssemblyErrors (KEY link, KEY contig, KEY key, Array dna, Array dnaR, int a1, int a2, BOOL reverse) ;
static void showAssemblyErrors (LOOK look, float x, SEG *seg, BOOL reverse) ;
static void showVirtualDna(LOOK look, float x, SEG *seg, Array dna, BOOL rev) ;

static void fMapTraceAllErrors (void) ;
static void contigDoMakePair (KEY contig1, KEY contig2) ;
static void fMapToggleTraceAssembly (void) ;
static void fMapTraceReAssembleContig (LOOK look, KEY contig) ;
static void fMapTraceRemoveContig (LOOK look, KEY contig) ;
static void fMapShowVirtualMenu (KEY key, int box) ;
static void showTaggedBases (LOOK look, float x, SEG *seg, BOOL reverse) ;
static void fMapTraceCleanContig (KEY contig) ;
static void fMapTraceAddContig (KEY link, KEY contigAbove, KEY newContig, int ll) ;

/************************************************************/

static BOOL isReversed = FALSE ;
static LOOK x1Look = 0 ;
static int selectedTraceMagic = 455424;
static int multipletMagic = 25445;

static FREEOPT virtualBoxMenu[] = {
  {4, "Virtual box menu"},
  {'t', "Who am I"},
  {'f', "Looks like (common 15-mers)"},
  {'s', "Move to new contig"},
  {'S', "Discard"}
  /*  'c', "Clip end here"  this is bugged */
} ;

#define BC_TAG 0x80

typedef struct { char sup, sdo, bup, bdo, nup, ndo ; } QUALITY ;
static int newErrorTracking = 0 ;

/***************************************************************************************/

static int x1Order (const void *a, const void *b)
{
  int 
    a1 = arrp(x1Look->segs, *(const int*)a, SEG)->x1,
    a2 = arrp(x1Look->segs, *(const int*)a, SEG)->x2,
    b1 = arrp(x1Look->segs, *(const int*)b, SEG)->x1,
    b2 = arrp(x1Look->segs, *(const int*)b, SEG)->x2,
    au = a1 < a2 ? a1 : a2 ,
    bu = b1 < b2 ? b1 : b2 ;

  return au - bu ;
}

/***************************************************************************************/
#ifdef JUNK
static void fMapUnSelectTrace (void)
{ 
  FMAPLOOKGET("unselectTrace") ;
  
  graphAssRemove (&selectedTraceMagic) ;
  mapDrawColumns (look->map) ;
}
/***************************************************************************************/

static void fMapSelectTrace (void)
{
  KEYSET ks, aa, bb, cc ;
  void *v ;
  int i , j ;
  KEY k1, k2 ;
  KEYSET selectedTraces = 0 ;
  FMAPLOOKGET("fMapselectTrace") ;
  
  if (!keySetActive(&ks, &v))
    { messout ("First Select a keySet containing the Fragments you want to display") ;
      return ;
    }
  if (graphAssFind (&selectedTraceMagic, &v))
    { selectedTraces = (KEYSET) v ;
      keySetDestroy (selectedTraces) ;
    }
  bb = query (ks, "CLASS Sequence AND DNA") ;
  cc = query (ks, "CLASS DNA") ;
  i = arrayMax(cc) ;
  j = 0 ;
  for (i=0 ; i < keySetMax(cc) ; i++)
    { k1 = keySet (cc, i) ;
      if (dnaReClass (k1, &k2))
	keySet(cc, j++) = k2 ;
    }
  keySetMax (cc) = j ;
  
    
  aa = keySetOR (bb, cc) ;
  keySetDestroy (bb) ;
  keySetDestroy (cc) ;
  v = (void*) aa ;
  graphAssociate (&selectedTraceMagic, v) ;
  
  mapDrawColumns (look->map) ;
}
#endif
/***************************************************************************************/

static void arrow (float x, float y1, float y2, float tp)
{ float ymin, ymax ;

  ymin = (y1 < y2 ? y1 : y2 ) ;
  ymax = (y1 > y2 ? y1 : y2 ) ;
	   
  if (ymin < tp) ymin = tp ;
  if (ymax > mapGraphHeight)
    ymax = mapGraphHeight ;

  if (ymax > tp && ymin < mapGraphHeight-1) ;
  else return ;

  /* show begin */
  if (y1  > tp && y1 < mapGraphHeight-1)
    graphLine (x - 0.25, y1, x + 0.25, y1) ; 
  /* show end */
  if (y2  > tp && y2 < mapGraphHeight-1)
    { graphLine (x - 0.25, y2 + (y1 < y2 ? - .5 : .5) , x, y2) ; 
      graphLine (x + 0.25, y2 + (y1 < y2 ? - .5 : .5) , x, y2) ; 
    }
  /* show cuts */
  if (ymin == tp)
    { if (y1 > y2)
	{ graphLine (x - 0.25, ymin + 1, x, ymin + .5) ; 
	  graphLine (x + 0.25, ymin + 1, x, ymin + .5) ; 
	}
      graphCircle (x, ymin, .5) ;
    }
  if (ymax == mapGraphHeight)
    { if (y1 < y2)
	{ graphLine (x - 0.25, ymax - 1, x, ymax - .5) ; 
	  graphLine (x + 0.25, ymax - 1, x, ymax - .5) ; 
	}
      graphCircle (x, ymax, .5) ;
    }
  graphLine (x, ymin, x, ymax) ;
}

/***************************************************************************************/

static void showArrows (LOOK look, Array aa, float *offset, BOOL reverse)
{ int i, j, ix, ibox, tp = topMargin + 2, col ;
  float x, y1, y2, ymin, ymax ;
  SEG   *seg ; BOOL isUp ;
  BUMP  bump = bumpCreate (mapGraphWidth,0) ;
  
  for (j = 0 ; j < arrayMax(aa) ; j++)
    { i = arr(aa, j, int) ;
      seg = arrp(look->segs, i, SEG) ;
      
      if (seg->type & 1) 
	{ y2 = MAP2GRAPH(look->map,seg->x1) ; 
	  y1 = MAP2GRAPH(look->map,seg->x2) ;
	}
      else
	{ y1 = MAP2GRAPH(look->map,seg->x1) ; /* begin */
	  y2 = MAP2GRAPH(look->map,seg->x2) ; /* end */
	}
      ymin = (y1 < y2 ? y1 : y2 ) ;
      ymax = (y1 > y2 ? y1 : y2 ) ;
      
      if (ymax > tp && ymin < mapGraphHeight-1)
	{ array(look->boxIndex,ibox=graphBoxStart(),int) = i ;
	  ix = 0 ; 
	  if (ymin < tp) ymin = tp ;
	  if (ymax > mapGraphHeight) ymax = mapGraphHeight ;
	  bumpItem (bump,1,(ymax - ymin + 1) + 0.2, &ix, &ymin) ;
	  x = *offset + ix ;
	  if (reverse) 
	    arrow (x, ymax, ymin, tp) ;
	  else
	    arrow (x, ymin, ymax, tp) ;
	  
	  if (look->flag &  FLAG_VIRTUAL_ERRORS)
	    { isUp = (seg->type & 1)  /* UP obj */ ? TRUE : FALSE ;
	      showAssemblyErrors (look, x, seg, isUp) ;
	      if (seg->data.i &   ((unsigned int) 1 << 31))
		showTaggedBases (look, x, seg, isUp) ;
	    }
	  graphBoxEnd () ;
	  graphBoxFreeMenu (ibox, fMapShowVirtualMenu, virtualBoxMenu) ;
	  
	  col = TRANSPARENT ; 
	  i = seg->data.i >> 24 ;
	  if (i)
	    { if (i & 1) col = PALEYELLOW ;  /* stolen */
	      if (i & 2) col = GREEN ; /* proposed */
	      if (i & 4) col = CYAN ;  /* vector */
	      if (i & 8) col = LIGHTGREEN ; /* new read */
	      /* 28, 29, 30 for tag vector, 31 for tagged_base */
	    }
	  i = (seg->data.i >> 22) & 3 ;
	  switch (i)
	    { 
	      case 1 : col = PALEYELLOW ;  break ;
	      case 2 : col = PALEBLUE ;  break ;
	      case 3 : col = PALEGREEN ;  break ;
	      /* 22, 23 are for queryColor */
	    }
	  if (col== TRANSPARENT)  
	    { 
	      i = (seg->data.i >> 16) & 0x3f ;  /* registered in fmapcontrol */
	      if (i < TRANSPARENT) col = i ;
	    }

	  /*	  if ((look->flag & FLAG_COLOR_CONTIGS) && */	      
	  graphBoxDraw (ibox, BLACK, col) ;
	}
    }
  *offset += bumpMax (bump) + 1 ;
  bumpDestroy (bump) ;
}

/***************************************************************************************/
/***************************************************************************************/

void fMapVirtualColor (LOOK look)
{ Array dna = look->dna ;
  int i, j ; char *cp, *cq ;
  float y ;
  if (!dna || !look->colors || !arrayMax(dna) ||
      look->min < 0 || look->min > arrayMax(dna) ||
      look->max < 0 || look->max > arrayMax(dna) )
    return ;

  i = look->min ; j = 0 ;
  cp = arrp (dna, i, char) ;
  cq = arrp (look->colors, i, char) ;
  for (; i < look->max ; i++, cp++, cq++)
    if (!*cp) 
      { if (!j) j = i ;
	*cq |= TINT_LIGHTGRAY ;
      }
    else 
      { if (i && j)
	  { y = MAP2GRAPH(look->map, (i + j)/2) ;
	    if (i > j + 20) graphLine(1,y,mapGraphWidth,y) ;
	    j = 0 ;
	  }
      }
}

/***************************************************************************************/
/***************************************************************************************/

void fMapShowPrimers (LOOK look, float *offset)
{ int is, ibox ;
  float x = *offset, y, oldw = graphLinewidth (.2) ;
  SEG *ss ;
  BOOL isUp, found = FALSE ;

  graphColor (BLACK) ;
  
  for (is=0 ; is < arrayMax(look->segs) ; is++)
    { ss = arrp(look->segs, is, SEG) ;
      if (ss->type != PRIMER && ss->type != PRIMER_UP)
	continue ;
      isUp = ss->type & 0x01 ? TRUE : FALSE ;
      y = MAP2GRAPH(look->map, isUp ? ss->x2 - 3 : ss->x1 + 3 ) ;
      if (y > mapGraphHeight || y < topMargin)
	continue ;
      array(look->boxIndex,ibox=graphBoxStart(),int) = is ; 
      if (isUp)
	{ graphLine (x, y, x+1.3, y - .7) ;
	  graphLine (x, y, x+1.3, y) ;
	  graphLine (x+1.3, y, x+1.3, y -.7) ;
	}
      else
	{ graphLine (x, y, x+1.3, y + .7) ;
	  graphLine (x, y, x+1.3, y + .7) ;
	  graphLine (x+1.3, y, x+1.3, y +.7) ;
	}
      graphBoxEnd() ;
      graphBoxDraw (ibox, BLACK, TRANSPARENT) ;
      found = TRUE ;
    }
  graphColor (BLACK) ;
  oldw = graphLinewidth (oldw) ;
  if (found)
    *offset += 3 ;
}

/***************************************************************************************/
/***************************************************************************************/

static char* mainCloneName(KEY *kp)
{ static KEY mainClone = 0 ;
  
  if (!mainClone)
    { KEYSET clones ;
      clones = query (0, ">? Clone Main_Clone") ;
      if (keySetMax(clones) == 1)
	mainClone = keySet (clones, 0) ;
      keySetDestroy (clones) ;
    }
  if (kp) *kp =  mainClone ;
  return mainClone > 1 ? name(mainClone) : "" ;
}


typedef struct keyseg { KEY key, subclone ; int type, i, g, x1, x2 ; } MPLT ;

void fMapShowVirtualMultiplets (LOOK look, float *offset)
{ int j, is, ix, ibox, tp = topMargin + 2 , x1, x2, g ;
  float x, y1, y2, ymin, yymin, ymax ;
  SEG *ss ;
  BUMP  bump ;
  Array mm = look->multiplets ;
  MPLT* mpp ;
  BOOL foundArrow ;
  int insertLength = 0 ;

  if (!insertLength)
    { KEY key, _subclone_length ;
      OBJ obj ;

      lexaddkey ("Subclone_Length", &_subclone_length, 0) ;
      insertLength = 3500 ;
      if (*mainCloneName (&key)  &&
	  (obj = bsCreate(key)))
	{ bsGetData (obj, _subclone_length, _Int, &insertLength) ;
	  bsDestroy (obj) ;
	}
    }
  if (!mm)
    return ;
  bump = bumpCreate (mapGraphWidth,0) ;

  for (is=0 ; is < arrayMax(look->segs) ; is++)
    { ss = arrp(look->segs, is, SEG) ;
      if (ss->type != VIRTUAL_MULTIPLET_TAG &&
	  ss->type != VIRTUAL_MULTIPLET_TAG_UP)
	continue ;
      ymin = MAP2GRAPH(look->map, ss->x1) ;
      ymax = MAP2GRAPH(look->map, ss->x2) ;
      if (ymin > mapGraphHeight || ymax < tp)
	continue ;
      j = ss->data.i & 0xffff ;
      mpp = arrp(mm, j, MPLT) ;
      if (mpp->type != VIRTUAL_SUB_SEQUENCE_TAG_UP) 
	continue ;
      g = mpp->g ; foundArrow = FALSE ;
      for (; mpp->g == g ; mpp++)
	{ x1 = mpp->x1 ; x2 = mpp->x2 ;
	  y1 = MAP2GRAPH(look->map,x1) ;
	  y2 = MAP2GRAPH(look->map,x2) ;
	  if ((y1 > tp && y1 < mapGraphHeight-1) ||
	      (y2 > tp && y2 < mapGraphHeight-1))
	    { foundArrow = TRUE ; break ;}
	}
      if (!foundArrow) {} ; /* continue ; show red line even if full across screen   */
      ymin = ymin > tp ? ymin : tp ;
      ymax = ymax < mapGraphHeight ? ymax : mapGraphHeight ;

      yymin = ymin ; ix = 0 ;
      bumpItem (bump,1,(ymax - ymin + 1) + 0.2, &ix, &yymin) ;
      x = *offset + ix ;
	
      array(look->boxIndex,ibox=graphBoxStart(),int) = is ; 
      graphColor (ss->x2 > ss->x1 + insertLength ? RED : GREEN) ;
      graphLine (x, ymin, x, ymax) ;
      graphColor (BLACK) ;

      j = ss->data.i & 0xffff ;
      mpp = arrp(mm, j, MPLT) ;
      g = mpp->g ;
      for (; mpp->g == g ; mpp++)
	{ 
	  switch (ss->type)
	    { 
	    case VIRTUAL_MULTIPLET_TAG:
	      x1 = ss->x1 + mpp->x1 ; x2 = ss->x1 + mpp->x2 ;
	      break ;
	    case VIRTUAL_MULTIPLET_TAG_UP:
	      x1 = ss->x2 - mpp->x1 ; x2 = ss->x2 - mpp->x2 ;
	      break ;
	    default: break ;
	    }
	  y1 = MAP2GRAPH(look->map,x1) ;
	  y2 = MAP2GRAPH(look->map,x2) ;

	  if ((y1 > tp && y1 < mapGraphHeight-1) ||
	      (y2 > tp && y2 < mapGraphHeight-1))
	      arrow(x, y1, y2, tp) ;
	}
      graphBoxEnd() ;
    }

  *offset += bumpMax (bump) + 1 ;
  bumpDestroy (bump) ;
  fMapVirtualColor(look) ;

/*  fMapShowPrimers (look, offset) ; chain these 2 columns */
}

/***************************************************************************************/

void fMapSelectVirtualMultiplet (LOOK look, SEG *seg) 
{ MPLT* mpp ;
  Array aa ;
  int j, g ;
  Stack s ;

  if (seg->type != VIRTUAL_MULTIPLET_TAG &&
      seg->type != VIRTUAL_MULTIPLET_TAG_UP)      
    return ;
  j = seg->data.i & 0xffff ;
  s = stackCreate(0) ;
  aa = look->multiplets ;
  mpp = arrp(aa, j, MPLT) ;
  g = mpp->g ;
  do
    { catText(s, 
	      messprintf("%s %d:%d ;",
			 name(mpp->key), mpp->x1, mpp->x2)) ;
    } while  ((++mpp)->g == g) ;
  strncpy (look->segTextBuf, stackText(s,0), 125) ;
  stackDestroy (s) ;
}

/***************************************************************************************/

static int mpltOrder (const void *a, const void *b)
{ 
  KEY ka = ((const MPLT*)a)->subclone, kb = ((const MPLT*)b)->subclone ;
  return ka > kb ? 1 : -1 ;
}

void fMapTraceFindMultiplets (LOOK look) 
{ OBJ obj ; KEY subclone ;
  int i, j , i1, i2, i3, g, min, max ;
  Array segs = look->segs ;
  MPLT *mmp, *mmp1, *mmp2 ;
  SEG *seg ;
  
  Array 
    mm2, mm = arrayCreate (32, MPLT) ;
  
  for (j = 0, i = 1 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      
      switch (seg->type)
	{
	case VIRTUAL_SUB_SEQUENCE_TAG:
	case VIRTUAL_SUB_SEQUENCE_TAG_UP:
	  if ((obj = bsCreate (seg->key)))
	    { if (bsGetKey (obj, _Subclone, &subclone))
		{ mmp = arrayp (mm, j++, MPLT) ;
		  mmp->subclone = subclone ;
		  mmp->key = seg->key ;
		  mmp->type = seg->type ;
		  mmp->i = i ;
		  switch (seg->type)
		    {
		    case VIRTUAL_SUB_SEQUENCE_TAG:
		      mmp->x1 = seg->x1 ;
		      mmp->x2 = seg->x2 ;
		      break ;
		    case VIRTUAL_SUB_SEQUENCE_TAG_UP:
		      mmp->x1 = seg->x2 ;
		      mmp->x2 = seg->x1 ;
		      break ;
		    default: break ;
		    }
		}
	      bsDestroy (obj) ;
	    }
	default: break ;
	}
    }

  if (!arrayMax(mm))
    goto abort ;

  arraySort (mm, mpltOrder) ;
  arrayCompress (mm) ;

  mmp1 = 0 ; i2 = 0 ;
  mm2 = arrayCreate (32, MPLT) ;
  for (i = 0 , mmp = arrp(mm,0, MPLT) ; i < arrayMax(mm) ;  mmp++, i++)
    { subclone = mmp->subclone ;
      mmp1 = mmp + 1 ; i1 = i + 1 ;
      while (i1 < arrayMax(mm) &&
	     mmp1->subclone == subclone)
	{ i1++ ; mmp1++ ; }
      seg = arrp(segs, mmp->i, SEG) ;
      if (i1 > i + 1)
	{ mmp->type = seg->type | 1 ;
	  min = seg->x1 ;
	  max = seg->x2 ;
	  g = i2 ;
	  if (min > max) { max = seg->x1 ; min = seg->x2 ;}
	  for (i3 = i ; i3 < i1 ; i3++)
	    { mmp1 = mmp + i3 - i ;
	      seg = arrp(segs, mmp1->i, SEG) ;
	      seg->data.i = (seg->data.i & 0xffff0000) | g ; /* keep flags color | proposed */
	      seg->parent = subclone ;
	      mmp2 = arrayp(mm2, i2++, MPLT) ;
	      *mmp2 = *mmp1 ;
	      mmp2->g = mmp->i ;
	      if (min > seg->x1) min = seg->x1 ;
	      if (min > seg->x2) min = seg->x2 ;
	      if (max < seg->x1) max = seg->x1 ;
	      if (max < seg->x2) max = seg->x2 ;
	    }   /* i make a loop */
	  i2 -= i3 - i ; /* bactrack */
	  for (i3 = i ; i3 < i1 ; i3++)
	    { mmp2 = arrayp(mm2, i2++, MPLT) ;
	      mmp2->x1 -= min ;
	      mmp2->x2 -= min ;
	    }   /* i make a loop */
	  seg = arrayp(segs, arrayMax(segs), SEG) ;
	  seg->type = VIRTUAL_MULTIPLET_TAG ;
	  seg->parent = subclone ;
	  seg->key = subclone ;
	  seg->x1 = min ;
	  seg->x2 = max ;
	  seg->data.i = (seg->data.i & 0xffff0000) | g ; /* keep flags proposed */
	}
      if (i1 > i + 1)
	{ mmp = mmp + i1 - i - 2 ;
	  i = i1 - 2 ; 
	}
    }
    
  look->multiplets = mm2 ;
  
 abort:
  arrayDestroy (mm) ;
}

/***************************************************************************************/
/***************************************************************************************/

static void fMapDoShowVirtual (LOOK look, float *offset, int type)
{
  float y1, y2, ymin, ymax ;
  int   i, ibox, dummy ;
  SEG   *seg = 0 ;
  int  virtualParent = 0 ;
  KEYSET selectedTraces = 0 ;
  Array 
    indices = arrayCreate (32, int) , 
    indicesR = arrayCreate (32, int) ;

  *offset += 1.25 ;
/*
  if (graphAssFind (&selectedTraceMagic, &v))
     { selectedTraces = (KEYSET) v ;
       graphButton ("UnSelect", fMapUnSelectTrace, *offset - 2, topMargin) ;
     }
  else
    graphButton ("Select", fMapSelectTrace, *offset, topMargin) ;
*/
  for (i = 1 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (selectedTraces &&
	  !keySetFind (selectedTraces, seg->key, &dummy))
	continue ;
        
      if (seg->type & 1) 
	{ y2 = MAP2GRAPH(look->map,seg->x1) ; 
	  y1 = MAP2GRAPH(look->map,seg->x2) ;
	}
      else
	{ y1 = MAP2GRAPH(look->map,seg->x1) ; /* begin */
	  y2 = MAP2GRAPH(look->map,seg->x2) ; /* end */
	}
      ymin = (y1 < y2 ? y1 : y2 ) ;
      ymax = (y1 > y2 ? y1 : y2 ) ;
      
      if (ymax > topMargin && ymin < mapGraphHeight-1)
	{
	  if (seg->type == type)
	    array (indices, arrayMax(indices), int) = i ;
	  else if (seg->type == type + 1)
	    array (indicesR, arrayMax(indicesR), int) = i ;
	  else if (seg->type == VIRTUAL_PARENT_SEQUENCE_TAG)
	    virtualParent = i ;
	}
    }
  
  if (seg && type == VIRTUAL_SUB_SEQUENCE_TAG  && virtualParent)
    { array(look->boxIndex,ibox=graphBoxStart(),int) = virtualParent ;
      y1 = MAP2GRAPH(look->map,seg->x1) ;
      y2 = MAP2GRAPH(look->map,seg->x2) ;
      graphRectangle(*offset, y1, *offset + .8 , y2) ;
      graphBoxEnd() ;
      graphBoxDraw(ibox, BLACK , LIGHTGREEN) ;
      *offset += 1.5 ;
    }

  if (!arrayMax(indices) &&
      !arrayMax(indicesR) )
    goto abort ;
    
  /*   look->flag & FLAG_COMPLEMENT ? TRUE : FALSE ; */

  x1Look = look ;
  arraySort (indices, x1Order) ;
  arraySort (indicesR, x1Order) ;

  graphColor (BLACK) ;

  showArrows(look, indicesR, offset, TRUE) ;
  showArrows(look, indices, offset, FALSE) ;
  *offset += 1.5 ;

 abort:
  arrayDestroy (indices) ;
  arrayDestroy (indicesR) ;  
}

extern void fMapToggleDna (void) ;

void fMapShowAlignments (LOOK look, float *offset)
{
  fMapDoShowVirtual (look, offset, VIRTUAL_ALIGNED_TAG) ;
}

static void fMapTraceToggleAssemblyDna (void)
{ FMAPLOOKGET("fMapTraceToggleDna") ;

  mapColSetByName ("Assembly DNA", 2) ;
  fMapDraw(look,0) ; 
}

void fMapShowPreviousContigs (LOOK look, float *offset)
{ fMapDoShowVirtual (look, offset, VIRTUAL_PREVIOUS_CONTIG_TAG) ;
}

static void fMapToggleTraceAssembly (void)
{ FMAPLOOKGET("toggleTraceAssembly") ;

  mapColSetByName("Trace Assembly", 2) ;
  fMapDraw (look, 0) ;
}

void fMapShowTraceAssembly (LOOK look, float *offset)
{ float x = *offset, y = MAP2GRAPH (look->map, look->map->centre) ;
  fMapDoShowVirtual (look, offset, VIRTUAL_SUB_SEQUENCE_TAG) ;
  graphLine (x, y, *offset, y) ;
}


/***************************************************************************************/
/***************************************************************************************/

static void fMapFindVector (SEG *seg)
{ KEY key = seg->key, dnaKey ;
  OBJ obj = 0 ;
  int cTop, cEnd, vTop, vEnd, maxDna ;

  if (seg->data.i &  (1 << 28))  /* no recursion */
    return ;
  seg->data.i |=  (1 << 28) ;
  seg->data.i &= ~(1 << 29) ;
  seg->data.i &= ~(1 << 30) ;
  
  if (!(obj = bsCreate (key)))
    return ;
  if (bsGetData (obj, _Clipping, _Int, &cTop) &&
      bsGetData (obj, _bsRight, _Int, &cEnd) &&
      bsGetData (obj, _Vector_Clipping, _Int, &vTop) &&
      bsGetData (obj, _bsRight, _Int, &vEnd) &&
      bsGetKey (obj, _DNA, &dnaKey) &&
      bsGetData (obj, _bsRight, _Int, &maxDna))
    { if ((cTop - vTop) * (cTop - vTop) < 10 && cTop > 10)
	seg->data.i |= (1 << 29) ;
      if ((cEnd - vEnd) * (cEnd - vEnd) < 10 && cEnd < maxDna - 10)  
	seg->data.i |= (1 << 30) ;
    }
  bsDestroy (obj) ;
}

/*************/

static void showAssemblyErrors (LOOK look, float x, SEG *seg, BOOL reverse)
{ MAP map = look->map ;
  KEY key = seg->key ;
  char * complement = look->flag & FLAG_COMPLEMENT  ? (char*)0 : ((char*)0) + (1 << 23) ;
  char* vp = complement + (int)key ;
  void *wp ;
  A_ERR *ep ;
  float y, yy, z1, z2 ;
  Array a ;
  int i, color = 0, oldColor = -1 ;
  float a1 = 0, a2 = 0 ;
  int x1 = seg->x1, x2 = seg->x2 ;
  Array dna = look->dna, dnaR ;

  if (!look->dnaR)
    { look->dnaR = dnaCopy (look->dna) ;
      reverseComplement (look->dnaR) ;
    }
  dnaR = look->dnaR ;

  if (x1 < x2)
    { z1 = x1 ; z2 = x2 ; }
  else
    { z1 = x2 ; z2 = x1 ; }
	
  if (!assExists (look->virtualErrors))
    look->virtualErrors = assCreate () ;
  if (assFind(look->virtualErrors, vp, &wp))
    a = (Array) wp ;
  else
    { a = findAssemblyErrors (look->seqKey, 0, key, dna, dnaR, x1, x2, reverse) ;
      wp = (void*) a ;
      assInsert(look->virtualErrors, vp, wp) ;
    }

  if (!(seg->data.i &   (1 << 28)))
    fMapFindVector(seg) ;
  if (seg->data.i &   (1 << (reverse ? 30 : 29)))
    { y = MAP2GRAPH(map, seg->x1) ;
      if (y > topMargin && y < mapGraphHeight)
	{ graphColor (ORANGE) ;
	graphCircle (x, y, .6) ;
	graphColor (BLACK) ;
	}
    }
  if (seg->data.i &  (1 << (reverse ? 29 : 30)))
    { y = MAP2GRAPH(map, seg->x2) ;
      if (y > topMargin && y < mapGraphHeight)
	{ graphColor (ORANGE) ;
	graphCircle (x, y, .6) ;
	graphColor (BLACK) ;
	}
    }
    
  if (!arrayExists(a) || !arrayMax(a) || a->size != sizeof(A_ERR))
    return ;
  ep = arrp(a, 0, A_ERR)  - 1 ;
  i = arrayMax(a) ; a2 = 0 ;
  while (ep++, i--)
    {
      if (ep->iLong > z2)
	break ;
      if (ep->iLong < z1)
	continue ;
      switch(ep->type)
	{
	case ERREUR: 
	  color = RED ; 
	  break ;
	case INSERTION: /* _AVANT: */
/*	case INSERTION_APRES: */
	case INSERTION_DOUBLE:  /*_AVANT:  */
/*	case INSERTION_DOUBLE_APRES:  */
	case TROU:
	case TROU_DOUBLE:
	  color = BLUE ; 
	  break ;
	case AMBIGUE:
	default:
	  color = DARKGREEN ; 
	  break ;
	}
      y = ep->iLong - .5 ;
      yy = ep->iLong + .5 ;
	
      y = MAP2GRAPH(map, y) ;
      yy = MAP2GRAPH(map, yy) ;
      if (yy - y < .3) yy = y + .3 ;
      if (y > topMargin + 2 && y < mapGraphHeight )
	{ 
	  if (oldColor == -1) a1 = y ;
	  if (y < a2 - .1) /* && (color == oldColor */
	    { a2 = yy ; /* coalesce the 2 boxes */
	      if (color != oldColor)
		oldColor = DARKRED ;
	    }
	  else if (y < a2 && color == oldColor)
	    { a2 = yy ; /* coalesce the 2 boxes */
	      if (color != oldColor)
		oldColor = DARKRED ;
	    }
	  else
	    { if (oldColor != -1)
		{ graphColor(oldColor) ;
		  graphFillRectangle(x - .25, a1, x + .25, a2) ;
		  graphColor (BLACK) ;
		}
	      /* register */
	      a1 = y ; a2 = yy ; oldColor = color ;
	    }
	}
    }
  
  if (oldColor != -1)
    { graphColor(oldColor) ;
      graphFillRectangle(x - .25, a1, x + .25, a2) ;
      graphColor (BLACK) ;
    }
}

/***************************************************************************************/

static Array findAssemblyErrors (KEY link, KEY contig, KEY key, Array dnaD, Array dnaR, int a1, int a2, 
				 BOOL upSequence)
{ Array a = 0, mydna = 0 ;
  OBJ obj ;
  int x1, x2, x3, b1, b2, b3, u ;
  KEY mydnakey = 0 ;

  if ((obj = bsCreate(key)))
    { if (bsGetKey (obj, _DNA, &mydnakey))
	{ if ((mydna = blyDnaGet(link, contig, mydnakey)))
	    {
	      x1 = 1 ; x2 = arrayMax (mydna) ;
		{
		  /* the window may have only part of the dna */
		  if (a1 < 0) 
		    {
		      if (upSequence) 
			{ x2 += a1 ; a1 -= a1 ; }
		      else
			{ x1 -= a1 ; a1 -= a1 ; }
		    }
		  if ((u = a2 - arrayMax(dnaD)) > 0)
		    {
		      if (upSequence) 
			{ x1 += u ; a2 -= u ; }
		      else
			{ x2 -= u ; a2 -= u ; }
		    }
		  if (x1 < x2 && a1 < a2)
		    { 
		      if (upSequence) 
			{ b1 = a2 ; b2 = a1 ; }
		      else 
			{ b1 = a1 ; b2 = a2 ; }
		      b3 = b2 ; x3 = x2 ; x1-- ;
		      a = baseCallCptErreur (dnaD, dnaR, mydna, upSequence, 
					     b1, b2, b3,
					     &x1, &x2, &x3, newErrorTracking) ;
		    }
		  else
		    a = arrayCreate (10, A_ERR) ;
		}
	    }
	}
      bsDestroy (obj) ;
    }
  return a ;
}

/***************************************************************************************/
/***************************************************************************************/

static Array findTaggedBases (KEY key)
{ Array a = 0 ;
  static Array units = 0 ;
  OBJ obj ;
  char cc, *cp ;
  int i, j = 0, clipTop = 1, x ;
 
  obj = bsCreate (key) ;
  if (!obj)
    return 0 ;
  units = arrayReCreate (units, 12, BSunit) ;
  bsGetData (obj, _Clipping, _Int, &clipTop) ;
  if (bsFindTag (obj, _Significant_bases) &&
      bsFlatten (obj, 2, units) &&
      arrayMax (units))
    { a = arrayCreate (20, KEY) ;
      for (i = 0 ; i < arrayMax (units) ; i += 2)
	{ x = arr (units, i, BSunit).i ;  
	  cp = arr (units, i + 1, BSunit).s ;  
	  cc = dnaEncodeChar [(int)*cp] ;
	  if (x >= clipTop) 
	    array(a, j++, KEY) = (cc << 24) + x - clipTop ;
	}
    }
  bsDestroy (obj) ;
  return a ;
}

/***************************************************************************************/

#define A_COLOR GREEN
#define T_COLOR RED
#define G_COLOR DARKGRAY /* lightgray in base call explicit */
#define C_COLOR CYAN
#define N_COLOR YELLOW

static void showTaggedBases (LOOK look, float x, SEG *seg, BOOL reverse)
{ MAP map = look->map ;
  KEY key = seg->key ;
  char* vp = assVoid(key) ;
  void *wp ;
  KEY *ip ;
  float y, yy ;
  Array a ;
  int i, color = BLACK ;
  int z, z1, z2 ;
  int x1 = seg->x1, x2 = seg->x2 ;

  if (x1 < x2)
    { z1 = x1 ; z2 = x2 ; }
  else
    { z1 = x2 ; z2 = x1 ; }
	
  if (!assExists (look->taggedBases))
    look->taggedBases = assCreate () ;
  if (assFind(look->taggedBases, vp, &wp))
    a = (Array) wp ;
  else
    { a = findTaggedBases (key) ;
      wp = (void*) a ;
      assInsert(look->taggedBases, vp, wp) ;
    }
    
  if (!arrayExists(a) || !arrayMax(a) || a->size != sizeof(KEY))
    return ;
  ip = arrp(a, 0, KEY)  - 1 ;
  i = arrayMax(a) ;
  while (ip++, i--)
    { z = *ip & 0xffffff ;
      if (reverse)
	z = x2 - z ;
      else
	z = x1 + z ;
      if (z > z2)
	break ;
      if (z < z1)
	continue ;
      if (!reverse)
	switch(((*ip) >> 24) & 0x0f)
	  {
	  case A_: color = A_COLOR ; break ;
	  case T_: color = T_COLOR ; break ;
	  case G_: color = G_COLOR ; break ;
	  case C_: color = C_COLOR ; break ;
	  }
      else
	switch(((*ip) >> 24) & 0x0f)
	  {
	  case A_: color = T_COLOR ; break ;
	  case T_: color = A_COLOR ; break ;
	  case G_: color = C_COLOR ; break ;
	  case C_: color = G_COLOR ; break ;
	  }
	
      y = z - .5 ;
      yy = z + .5 ;
	
      y = MAP2GRAPH(map, y) ;
      yy = MAP2GRAPH(map, yy) ;
      if (yy - y < .45) yy = y + .45 ;
      if (y > topMargin + 2 && y < mapGraphHeight )
	{ graphColor(color) ;
	  graphFillRectangle(x - .45, y, x + .45, yy) ;
	  graphColor (BLACK) ;
	}
    }
}

/*****************************************************************************/
/*****************************************************************************/

static void fMapShowVirtualMenu (KEY key, int box) 
{ KEY tag ;
  SEG *seg ;
  FMAPLOOKGET("fMapShowVirtualMenu") ;

  seg = BOXSEG(box) ;
  if (!seg)
    return ;

  switch (key)
    {
    case 't':
      display (seg->key, look->seqKey, TREE) ;
      break ;
    case 'f':
      { KEY dnaKey ; Array segDna ; 

	look->zoneMin = 0 ;
	look->zoneMax = arrayMax(look->dna) ;
	if (!dnaSubClass (seg->key, &dnaKey) ||
	    !(segDna = dnaGet (dnaKey)))
	  return ;
	dnaAlignFindRepeat (segDna, look->dna, look->colors, 50, 400, 15) ;
	mapColSetByName ("Summary bar", TRUE) ;
	fMapDraw (look, 0) ; /* NO recompute */
	arrayDestroy (segDna) ;
	return ;
      }
      break ;
    case 's': case 'S':
      switch (seg->type | 1)
	{
	case VIRTUAL_ALIGNED_TAG_UP:
	  tag = _Aligned_into ; key = 'S' ;
	  break ;
	case VIRTUAL_PREVIOUS_CONTIG_TAG_UP:
	  tag = _Later_part_of ; key = 'S' ;
	  break ;
	case VIRTUAL_SUB_SEQUENCE_TAG_UP:
	  tag = _Assembled_into ; 
	  break ;
	default:
	  tag = _Assembled_into ;
	  break ;
	}
      fMapTraceSuppress (seg->key, tag, key == 's' ? TRUE : FALSE) ;
      traceGraphDestroy () ;
      break ;
    case 'c':
 /* bugged because a1 a2 in contig not updated */
      { float xx, yy ; int z1 ;
	OBJ obj = bsUpdate (seg->key) ;
        int i, c1, c2, c3 ;

	if (!obj) return ;
    
	graphBoxAt (graphEventX, graphEventY, &xx, &yy) ;
	z1 = MAP2GRAPH(look->map, seg->x1) ;
	if (z1 < topMargin) yy += topMargin - z1 ;	  
	z1 = yy / look->map->mag ;

	if (bsGetData (obj, _Clipping, _Int, &c1))
	  bsGetData (obj, _bsRight, _Int, &c2) ;
	if (!bsGetData (obj, _Old_Clipping, _Int, &i))
	  { bsAddData (obj, _Old_Clipping, _Int, &c1) ;
	    bsAddData (obj, _bsRight, _Int, &c2) ;
	  }
	bsAddData (obj, _Clipping, _Int, &c1) ;
	c3 = c1 + z1 ;
	bsAddData (obj, _bsRight, _Int, &c3) ;
	bsSave (obj) ;
	
      }
      traceGraphDestroy () ;
      look->pleaseRecompute = TRUE ;
      fMapDraw (look, 0) ;

      break ;
    }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* attention a1 base incluse, a2 base exclue */
/* Attention les definitions doivent etre dans l'ordre inverse de la peinture
   faite par fMapDrawColor */
#define OT_TAG 1
#define AC_TAG 2
#define BAD_QU 4
#define UN_DO 8
#define UN_UP 16
#define NO_DO 32
#define NO_UP 64
#define TRAC 128
#define FAIT 256

/*****************************************************************************/

static void fMapComputeColor (LOOK look, int a1, int a2)
{ int i, j, x1, x2, y1, y2, n ;
  KEY key ;
  A_ERR *ep, *epmax ;
  QUALITY *qual ;
  SEG *segp ;
  Array inter = arrayCreate (a2 - a1, QUALITY), units ;
  Array errors, color = look->contigCol ;
  char * complement = look->flag & FLAG_COMPLEMENT  ? (char*)0 : ((char*)0) + (1 << 23) ;
  char *ind ;
  int *ip ;
  void *vp ;
  BOOL reverse ;

  if (!assExists (look->virtualErrors))
    goto abort ;
  units = look->segs ;
  i = arrayMax (units) - 1 ;
  if (i < 1)
    goto abort ;
  segp = arrp (units, 0, SEG) ; 
  arrayp (inter, a2 - a1 - 1, QUALITY)->sup = 0 ;
  while (segp++, i--)
    { if ((segp->type != VIRTUAL_SUB_SEQUENCE_TAG) && 
	  (segp->type != VIRTUAL_SUB_SEQUENCE_TAG_UP))
	continue ;
      key = segp->key ;
      x1 = segp->x1 ;
      x2 = segp->x2 ;
      if ((x1 < a1 && x2 < a1) || (x1 > a2 && x2 > a2))
	continue ;
      if (segp->type == VIRTUAL_SUB_SEQUENCE_TAG)
	reverse = FALSE ;
      else
	reverse = TRUE ;
      if (x1 < x2)
	{ y1 = x1 ;
	  y2 = x2 ;
	}
      else 
	{ y1 = x2 ;
	  y2 = x1 ;
	}
      if (!look->dnaR)
	{ look->dnaR = dnaCopy (look->dna) ;
	  reverseComplement (look->dnaR) ;
	}

      ind = complement + (int)key ;
      if (assFind (look->virtualErrors, ind, &vp))
	errors = (Array)vp ;
      else
	{ errors = findAssemblyErrors (look->seqKey, 0, key, look->dna, look->dnaR, x1, x2, reverse) ;
	  vp = (void*)errors ;
	  assInsert(look->virtualErrors, ind, vp) ;
	}
      if (!errors)
	continue ;
      y1 = y1 < a1 ? a1 : y1 ;
      y2 = y2 > a2 ? a2 : y2 ;
      j = y2 - y1 ;
      if (y1 - a1 >= arrayMax (inter)) continue ; /* mieg */
      qual = arrp (inter, y1 - a1, QUALITY) ;
      if (arrayMax (errors))
	{ ep = arrp (errors, 0, A_ERR) ;
	  epmax = arrp (errors, arrayMax (errors) - 1, A_ERR) + 1 ;
	}
      else 
	ep = epmax = 0 ;
      while (ep < epmax)
	{ if (ep->iLong < y1)
	    ep++ ;
	  else break ;
	}
      if (ep < epmax)
	n = ep->iLong ;
      else n = -1 ;
      while (j--)
	{ if (reverse)
	    (qual->sup)++ ;
	  else
	    (qual->sdo)++ ;
	  if (n == y1)
	    { switch (ep->type)
		{
		case ERREUR:
		case INSERTION: case INSERTION_DOUBLE: case INSERTION_TRIPLE:
		case TROU: case TROU_DOUBLE: case TROU_TRIPLE:
		  if (reverse)
		    (qual->bup)++ ;
		  else
		    (qual->bdo)++ ;
		  break ;
		case AMBIGUE:
		  if (reverse)
		    (qual->nup)++ ;
		  else
		    (qual->ndo)++ ;
		  break ;
		case TYPE80:
		  break ;
		}
	      ep++ ;
	      if (ep < epmax)
		n = ep->iLong ;
	      else n = -1 ;
	    }
	  qual++ ;
	  y1++ ;
	}
    }
  ip = arrp (color, a1, int) - 1 ;
  qual = arrp (inter, 0, QUALITY) - 1 ;
  j = a2 - a1 ;
  while (ip++, qual++, j--)
    { *ip = FAIT ;
      if (qual->sup + qual->sdo)
	*ip |= TRAC ;
      else continue ;
      if (qual->sup < 2)
	*ip |= NO_UP ;
      if (qual->sdo < 2)
	*ip |= NO_DO ;
      if ((100 * (qual->bup + qual->bdo + qual->nup + qual->ndo)) / (qual->sup + qual->sdo) > 50)
	*ip |= BAD_QU ;
      if (qual->sup == qual->nup)
	*ip |= UN_UP ;
      if (qual->sdo == qual->ndo)
	*ip |= UN_DO ;
    }
 abort:
  arrayDestroy (inter) ;
}

/*****************************************************************************/

static void fMapComputeTagColor (LOOK look, int a1, int a2)
{ OBJ obj = 0 ;
  SEG *segp ;
  KEY key ;
  int i, j, k, x1, x2, y1, y2, z1, z2, *ip, clip, dummy ;
  Array flags = 0 ;
  char *cp ;
  BSunit *u ;
  BOOL isUp ;

  i = arrayMax (look->segs) - 1 ;
  if (i < 1) return ;
  segp = arrp (look->segs, 0, SEG) ;
  while (segp++, i--)
    { 
      if (segp->type == VIRTUAL_SUB_SEQUENCE_TAG)
	isUp = FALSE ;
      else if (segp->type == VIRTUAL_SUB_SEQUENCE_TAG_UP)
	isUp = TRUE ;
      else
	continue ;
      key = segp->key ;
      x1 = segp->x1 ;
      x2 = segp->x2 ;
      if ((x1 < a1 && x2 < a1) || (x1 > a1 && x2 > a2))
	continue ;
      flags = arrayReCreate (flags, 10, BSunit) ;
      if ((obj = bsCreate (key)) && bsGetArray (obj, _Assembly_tags, flags, 3))
	{ if (!dnaAlignGetClip (look->seqKey, 0, key, &clip, &dummy))
	    clip = 1 ;
	  for (j = 0 ; j < arrayMax (flags) ; j += 3)
	    { u = arrp (flags, j, BSunit) ;
	      cp = u[0].s ;
	      y1 = u[1].i ;
	      y2 = u[2].i ;
	      if (!isUp)
		{ z1 = x1 + y1 - clip ;
		  z2 = x1 + y2 - clip ;
		}
	      else
		{ z1 = x2 - y2 + clip ;
		  z2 = x2 - y1 + clip ;
		}
	      if (z1 > z2) { int zz = z1 ; z1 = z2 ; z2 = zz ; }
	      if (z1 < 0) z1 = 0 ;
	      if (z2 > arrayMax (look->contigCol))
		z2 = arrayMax (look->contigCol) ;
	      k = z2 - z1 ;
	      if (k < 0) k = 0 ; /* to prevent crash */
	      ip = arrp (look->contigCol, z1, int) - 1 ;
	      if (!strcmp (cp, "Compression") || !strcmp (cp, "Ambiguity"))
		while (ip++, k--)
		  *ip |= AC_TAG ;
	      else
		while (ip++, k--)
		  *ip |= OT_TAG ;
	    }
	}
      bsDestroy (obj) ;
    }
  arrayDestroy (flags) ;
}

/*****************************************************************************/

static void fMapDoDrawColor (LOOK look, float offset, int min, 
			     int max, int type)
{ float x1 = 0, x2 = 0, y1 = 0, y2 = 0 ;
  int z1, z2, i, *ip, *ip0, typedr = 0 ;

/* Attention a l'ordre des dessins de couleurs par rapport #define */
  switch (type)
    {
    case OT_TAG:
      x1 = offset - .23 ;
      x2 = offset + .25 ;
      graphColor (LIGHTGREEN) ;
      typedr = OT_TAG ;
      break ;
    case AC_TAG:
      x1 = offset - .23 ;
      x2 = offset + .25 ;
      graphColor (RED) ;
      typedr = OT_TAG | AC_TAG ;
      break ;
    case BAD_QU:
      x1 = offset - .5 ;
      x2 = offset + .5 ;
      graphColor (RED) ;
      typedr = BAD_QU ;
      break ;
    case NO_UP:
      x1 = offset - .5 ;
      x2 = offset ;
      graphColor (LIGHTBLUE) ;
      typedr = NO_UP | UN_UP | BAD_QU ;
      break ;
    case NO_DO:
      x1 = offset ;
      x2 = offset + .5 ;
      graphColor (LIGHTBLUE) ;
      typedr = NO_DO | UN_DO | BAD_QU ;
      break ;
    case UN_UP:
      x1 = offset - .5 ;
      x2 = offset ;
      graphColor (BLACK) ;
      typedr = UN_UP | BAD_QU ;
      break ;
    case UN_DO:
      x1 = offset ;
      x2 = offset + .5 ;
      graphColor (BLACK) ;
      typedr = UN_DO | BAD_QU ;
      break ;
    case TRAC:
      x1 = offset - .5 ;
      x2 = offset + .5 ;
      graphColor (YELLOW) ;
      typedr = type | (type - 1) ;
      break ;
    case FAIT:
      x1 = offset - .5 ;
      x2 = offset + .5 ;
      graphColor (BLACK) ;
      typedr = type | (type - 1) ;
      break ;
    }
  ip0 = arrp (look->contigCol, 0, int) ;
  ip = ip0 + min ;
  i = max - min ;
  while (i)
    { while (i && !(*ip & type))
	{ ip++ ;
	  i-- ;
	}
      if (!i)
	break ;
      z1 = ip - ip0 ;
      while (i && (*ip & typedr))
	{ ip++ ;
	  i-- ;
	}
      z2 = ip - ip0 ;
      y1 = MAP2GRAPH (look->map, z1) ;
      y2 = MAP2GRAPH (look->map, z2) ;
      if (y2 < y1 + .2) 
	y2 = y1 + .2 ;
      graphFillRectangle (x1, y1, x2, y2) ;
    }
}

/*****************************************************************************/
/* Attention : ordre de peinture doit correspondre aux #define */
static void fMapDrawColor (LOOK look, float offset, int min, int max)
{ fMapDoDrawColor (look, offset, min, max, FAIT) ;
  fMapDoDrawColor (look, offset, min, max, TRAC) ;
  fMapDoDrawColor (look, offset, min, max, NO_UP) ;
  fMapDoDrawColor (look, offset, min, max, NO_DO) ;
  fMapDoDrawColor (look, offset, min, max, UN_UP) ;
  fMapDoDrawColor (look, offset, min, max, UN_DO) ;
  fMapDoDrawColor (look, offset, min, max, BAD_QU) ;
  fMapDoDrawColor (look, offset, min, max, AC_TAG) ;
  fMapDoDrawColor (look, offset, min, max, OT_TAG) ;
  graphColor (BLACK) ;
}

/*****************************************************************************/
/* penser a ajouter que si les tableaux d'erreurs sont detruits ou que des
modifications sur les sequences a l'ecran sont faites il faut remettre a 0 le 
tableau des couleurs */
static void fMapTraceColorReport (LOOK look, float offset, int min, int max)
{ int x1, x2, i, *ip, *ip0 ;

  if (!assExists (look->virtualErrors))
    look->virtualErrors = assCreate () ;
  if (min < 0)
    min = 0 ;
  i = arrayMax (look->dna) ;
  if (max > i)
    max = i ;
  if (!look->contigCol)
    { look->contigCol = arrayCreate (i, int) ;
      array (look->contigCol, i - 1, int) = 0 ;
    }
/* Normalement la fonction est appelee sur au plus un contig => noir par defaut
   i e il y a un contig sur la zone min-max  pour le tableau couleur et donc FAIT
   correspond aux zones ou il y a un contig */
  ip0 = arrp (look->contigCol, 0, int) ;
  ip = ip0 + min ;
  i = max - min ;
  while (i > 0)
    { while (i && (*ip & FAIT))
	{ ip++ ;
	  i-- ;
	}
      if (!i)
	break ;
      x1 = ip - ip0 ;
      while (i && !(*ip & FAIT))
	{ ip++ ;
	  i-- ;
	}
      x2 = ip - ip0 ;
      fMapComputeColor (look, x1, x2) ; /* x1 inclue, x2 exclue */
      fMapComputeTagColor (look, x1, x2) ;
    }
  fMapDrawColor (look, offset, min, max) ;
}

/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/

int fMapShowVirtualMultiTracePos ;
BOOL fMapFollowVirtualMultiTrace (int box) 
{
  SEG *seg, *s ;
  KEY kk, from ;
  KEYSET aa = 0 ;
  int i, j, n ; 
  SEQINFO *sinf ;
  FMAPLOOKGET("showVirtualTrace") ;
  
  seg = BOXSEG(box) ;
  if (!seg)
    return FALSE ;

  if (isDisplayBlocked())
    { display (seg->key, 0, 0) ;
    return TRUE ;
    }
  
  from = GRAPH2MAP (look->map, graphEventY + .2) ; /* + look->map->mag * .5) ; */
  
  if (look->seqKey == seg->key)
    { display (seg->key, from, DtMultiTrace) ;
    goto done ;
    }
  
  aa = queryKey (seg->key, "ABI ; >Assembled_into") ;
  if (!aa || keySetMax (aa) == 0)
    { keySetDestroy (aa) ;
    aa = queryKey (seg->key, "SCF_FILE") ;
    }
  
  /* check for whole look */
  if (aa && (i = keySetMax (aa)))
    while (i--)
      { kk = keySet (aa, i) ;
      if (look->seqKey == kk)
	{ display (kk, from, DtMultiTrace) ;
	goto done ;
	}
      }
  /* now check for a subcontig */
  if (aa && (i = keySetMax (aa)))
    while (i--)
      { kk = keySet (aa, i) ;
      s = arrp(look->segs, 0, SEG) - 1 ;
      j = arrayMax(look->segs) ;
      while (s++, j--)
	if (
	    (s->type == SEQUENCE ||
	     s->type == SEQUENCE_UP) &&
	    s->key == kk &&
	    s->data.i &&
	    (sinf = arrayp (look->seqInfo, s->data.i, SEQINFO)) &&
	    sinf->flags & SEQ_VIRTUAL_ERRORS)
	  { 
	    /* n1 = from ;  must cast to int before compare */
	    n = from - s->x1 ;
	    if (look->origin != s->x1)
	      { look->origin = s->x1 ;
		fMapDraw (look, 0) ;
	      }	    
	    display (kk, n, DtMultiTrace) ;
	    goto done ;
	  }
      }
  display (seg->key, look->seqKey,0) ;
 done:
  keySetDestroy (aa) ;
  
  return TRUE ;
}

/*****************************************************************************/

void fMapShowVirtualDna (LOOK look, float *offset)
{ BUMP  bump = bumpCreate (mapGraphWidth,0) ;
  float y1, y2, old ;
  int   i, j, ix, ibox, dummy ;
  SEG   *seg = 0 ;
  int  virtualParent = 0 ;
  KEYSET selectedTraces = 0 ;
  Array indices = arrayCreate (32, int) , indicesR = arrayCreate (32, int) ;
  void *v ;

  if (!(look->flag &  FLAG_VIRTUAL_ERRORS))
    return ;
  if (!mapColSetByName ("DNA Sequence", -1))
    return ;

  if (graphAssFind (&selectedTraceMagic, &v))
    selectedTraces = (KEYSET) v ;

  if (!look->dnaSkip)  /* force open the dna */
    { int box = graphBoxStart() ;
      old = *offset ;
      fMapShowDNA (look, offset) ;    
      graphBoxEnd() ;
      graphBoxClear (box) ;
      *offset = old ;  /* since i boxCleared */
    }

  *offset += 1.25 ;

  for (i = 1 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      switch (seg->type)
	{
	case VIRTUAL_SUB_SEQUENCE_TAG:
	case VIRTUAL_SUB_SEQUENCE_TAG_UP:
	case VIRTUAL_PARENT_SEQUENCE_TAG:
	  break ;
	default:
	  continue ;
	}
      if (selectedTraces &&
	  !keySetFind (selectedTraces, seg->key, &dummy))
	continue ;
      if (seg->type == VIRTUAL_SUB_SEQUENCE_TAG)
	array (indices, arrayMax(indices), int) = i ;
      else if (seg->type == VIRTUAL_SUB_SEQUENCE_TAG_UP)
	array (indicesR, arrayMax(indicesR), int) = i ;
      else if (seg->type == VIRTUAL_PARENT_SEQUENCE_TAG)
	virtualParent = i ;
    }
  
  if (seg && virtualParent)
    { array(look->boxIndex,ibox=graphBoxStart(),int) = virtualParent ;
      y1 = MAP2GRAPH(look->map,seg->x1) ;
      y2 = MAP2GRAPH(look->map,seg->x2) ;
      graphRectangle(*offset, y1, *offset + .8 , y2) ;
      graphBoxEnd() ;
      graphBoxDraw(ibox, BLACK , LIGHTGREEN) ;
      *offset += 1.5 ;
    }

  if (!arrayMax(indices) &&
      !arrayMax(indicesR) )
    goto abort ;
      
  isReversed = FALSE ; /* look->flag & FLAG_REVERSE ? TRUE : FALSE ; */

  x1Look = look ;
  arraySort (indices, x1Order) ;
  arraySort (indicesR, x1Order) ;

  graphColor (BLACK) ;
  graphTextFormat (FIXED_WIDTH) ;

  for (j = 0 ; j < arrayMax(indicesR) ; j++)
    { i = arr(indicesR, j, int) ;
      seg = arrp(look->segs, i, SEG) ;
      y1 = MAP2GRAPH(look->map,seg->x1) ;
      y2 = MAP2GRAPH(look->map,seg->x2) ;
   
      if (y2 > topMargin && y1 < mapGraphHeight-1)
	{ array(look->boxIndex,ibox=graphBoxStart(),int) = i ;
	  ix = 0 ; 
	  if (y1 < topMargin) y1 = topMargin ;
	  if (y2 > mapGraphHeight) y2 = mapGraphHeight ;
	  bumpItem (bump,1,(y2-y1+1)+0.2,&ix,&y1) ;
	  showVirtualDna(look, *offset + look->dnaWidth*ix + ix,
			  seg, look->dna, TRUE) ;
	  graphBoxEnd() ;
	  graphBoxFreeMenu (ibox, fMapShowVirtualMenu, virtualBoxMenu) ;
	}
      
    }

  *offset += 2 + (1 + look->dnaWidth) * bumpMax(bump) ;
  if (*offset > mapGraphWidth)
    goto abort ;
  bump = bumpReCreate (bump, mapGraphWidth, 0) ;
  
  for (j = 0 ; j < arrayMax(indices) ; j++)
    { i = arr(indices, j, int) ;
      seg = arrp(look->segs, i, SEG) ;
      y1 = MAP2GRAPH(look->map,seg->x1) ;
      y2 = MAP2GRAPH(look->map,seg->x2) ;
   
      if (y2 > topMargin && y1 < mapGraphHeight-1)
	{ array(look->boxIndex,ibox=graphBoxStart(),int) = i ;
	  ix = 0 ; 
	  if (y1 < topMargin) y1 = topMargin ;
	  if (y2 > mapGraphHeight) y2 = mapGraphHeight ;
	  bumpItem (bump,1,(y2-y1+1)+0.2,&ix,&y1) ;
	  showVirtualDna(look, *offset + look->dnaWidth*ix + ix,
			 seg, look->dna, FALSE) ;
	  graphBoxEnd() ;
	  graphBoxFreeMenu (ibox,  fMapShowVirtualMenu, virtualBoxMenu) ;
	}
    }

  *offset += 2 + (1 + look->dnaWidth) * bumpMax(bump) ;
  
 abort:
  bumpDestroy (bump) ;
  arrayDestroy(indices) ;
  arrayDestroy(indicesR) ;
  graphTextFormat (PLAIN_FORMAT) ;
}

/***************************************************************************************/

static void showVirtualDna(LOOK look, float x,
			       SEG *seg, Array dna, BOOL upSequence)
{ KEY key = seg->key ; int x1 = seg->x1 , x2 = seg->x2 ;
  char* complement = look->flag & FLAG_COMPLEMENT  ? (char*)0 : ((char*)0) + (1<< 23) ;
  int i, nn , dx ;
  char *vp = complement + (int)key ;
  void *wp ;
  A_ERR *ep ;
  float ddx = 0, y = 0 , y1 = topMargin , y2 = mapGraphHeight ;
  Array a, dnaR ;
  int start, skip, width ;
  MAP map = look->map ;
  char c, *buffer, buf[2] ;
  BOOL backwards ;
  
  if (!vp)
    return ;
  if (x > mapGraphWidth)
    return ;
  buf[1] = 0 ;
  y = MAP2GRAPH(map, x1) - 0.5 ;
  if (y > mapGraphHeight - 1)
    return ;
	
  if (!look->dnaR)
    { look->dnaR = dnaCopy (look->dna) ;
      reverseComplement (look->dnaR) ;
    }
  dnaR = look->dnaR ;
  if (!assExists (look->virtualErrors))
    look->virtualErrors = assCreate () ;
  if (assFind(look->virtualErrors, vp, &wp))
    a = (Array) wp ;
  else
    { a = findAssemblyErrors (look->seqKey, 0, key, dna,dnaR, x1, x2, upSequence) ;
      wp = (void*) a ;
      assInsert(look->virtualErrors, vp, wp) ;
    }
  
  if (!arrayExists(a))
    return ;

  start = look->dnaStart ;
  skip = look->dnaSkip ;
  width = look->dnaWidth ;

  nn = width + 12 ;
  buffer = messalloc(nn) ;
  memset(buffer, '.', nn - 1) ;
  buffer [width + 11] = 0 ;
  nn-- ;
  start += ((x1 - start)/skip)*skip ;
  while (start < x2)
    { y = MAP2GRAPH(look->map, start) - 0.5 ;
      if (y > y1)
	break ;
      start += skip ;
    }
  if (x1 > start)
    { dx = x1 - start ;
      if (dx < width)
	graphText(buffer + nn - (width - dx), x + dx, y + .3) ;
      start += skip ;
    }
  while (start + skip <= x2)
    { y = MAP2GRAPH(look->map, start) - 0.5 ;
      if (y > y2)
	break ;
      graphText(buffer + nn - width, x, y + .3) ;
      start += skip ;
    }
  if (start <= x2)
    { y = MAP2GRAPH(look->map, start) - 0.5 ;
      if (y < y2)
	{ dx = x2 - start + 1 ;
	  if (dx > width)
	    dx = width ;
	  graphText(buffer + nn - dx, x, y + .3) ;
	}
    }
  messfree (buffer) ;
  start = look->dnaStart ;
  if (!arrayExists(a) || !arrayMax(a))
    return ;
    
  backwards = FALSE ;
  if (look->flag & FLAG_COMPLEMENT_SUR_PLACE)
    upSequence = !upSequence ;
  if (look->flag & FLAG_COMPLEMENT)
    { backwards = !backwards ;
    }
  ep = arrp(a, 0, A_ERR)  - 1 ;
  i = arrayMax(a) ;
  while (ep++, i--)
    { c = dnaDecodeChar[upSequence ? (int)complementBase(ep->baseShort) : (int)ep->baseShort] ;
      ddx = 0 ;
      switch(ep->type)
	{
	case ERREUR: 
	  /* color = RED ;  */
	  c = ace_lower(c) ;
	  ddx = 0 ;
	  break ;
	case TROU:
	  /* color = BLUE ;  */
	  c = '*' ;
	  ddx = 0 ;
	  break ;
	case TROU_DOUBLE:
	case TROU_TRIPLE:
	  /* color = BLUE ;  */
	  c = '#' ;
	  ddx = 0 ;
	  break ;
	case INSERTION: /*_AVANT: */
	  /* color = BLUE ; */
	  c = ace_upper(c) ;
	  ddx = -.5 ;
	  break ;
/*
	case INSERTION_APRES:
	// color = BLUE ; 
	  c = ace_upper(c) ;
	  ddx = .5 ;
	  break ;  
*/
	case INSERTION_DOUBLE: /* _AVANT: */
	case INSERTION_TRIPLE: 
	  /* color = BLUE ; */
	  c = 'X' ;
	  ddx = -.5 ;
	  break ;
/*	case INSERTION_DOUBLE_APRES:
	case INSERTION_TRIPLE_APRES:
	// color = BLUE ; 
	  c = 'X' ;
	  ddx = .5 ;
	  break ;  */
	case AMBIGUE: 
	  /* color = GREEN ;  */
	  c = ace_lower(c) ;
	  ddx = 0 ;
	  break ;
	case TYPE80:
	  break ;
	}
      if (backwards)
	ddx *= -1 ;
      nn = ep->iLong ;
      dx = (nn - start) % skip ;
      y = nn + (backwards ? - dx : - dx) ;  /* not sure */
      y = MAP2GRAPH(map, y) - .5 ;
      buf[0] = c ;
      if (nn >= x1 && nn <= x2 && dx < width && y > y1 && y < y2)
	{ /* graphColor(BLACK) ;  color */
	  graphText(buf, x + dx + ddx, y) ;
	  /* graphColor (BLACK) ; */
	}
    }
}

/***************************************************************************************/
/***************************************************************************************/

static void fMapTraceAllErrors (void)
{ FMAPLOOKGET("fMapTraceAllErrors") ;

  look->flag ^= FLAG_VIRTUAL_ERRORS ; 
  fMapDraw(look,0) ; 
}

/***************************************************************************************/
/***************************************************************************************/

static void fMapTraceAddContig (KEY link, KEY contigAbove, KEY newContig, int ll)
{ KEY key ;
  OBJ Link = 0 ;
  int i, x1, x2, max ;
  Array units, order ;

  Link = bsCreate (link) ;

  units = arrayCreate(50, BSunit) ;
  bsGetArray (Link, _Subsequence, units, 3) ;
  bsDestroy (Link) ;
  order = arrayCreate (50, BSunit) ;
  max = 0 ;
  for (i = 0 ; i < arrayMax(units) ; i += 3)
    { key = arr (units, i, BSunit).k ;
      x1 = arr (units, i + 1, BSunit).i ;
      x2 = arr (units, i + 2, BSunit).i ;
      array (order, max++, BSunit).k = key ;
      array (order, max++, BSunit).i = x2 - x1 ;
      if (key == contigAbove)
	{ array (order, max++, BSunit).k = newContig ;
	  array (order, max++, BSunit).i = ll ;
	}
    }
  arrayDestroy (units) ;
  alignToolsAdjustLink (link, 0, order) ;
  arrayDestroy (order) ;
}

/***************************************************************************************/

static void fMapTraceCutContig (LOOK look, KEY contig, int where, char action)
{ KEY key, key1, key2, dnaKey, link = look->seqKey ;
  OBJ Link = 0 ;
  int i, x1, x2, max, max1, max2 ;
  char buff[7] ;
  Array units, order ;

  Link = bsCreate (link) ;
  if (!Link || !bsFindKey (Link, _Subsequence, contig) ||
      !bsGetData (Link, _bsRight, _Int, &x1) ||
      !bsGetData (Link, _bsRight, _Int, &x2))
    goto abort ;

  if (x1 > x2)
    where = x1 - x2 - where ;
  switch (action)
    {
    case 'c':
      sprintf (buff, "cut") ;
      break ;
    case 'u':
      sprintf (buff, "unjoin") ;
      break ;
    default:
      goto abort ;
    }
  if (!dnaAlignCutContig (contig, where, &key1, &key2, &max1, &max2, action))
    goto abort ;

  bsFindTag (Link, _Subsequence) ;
  units = arrayCreate(50, BSunit) ;
  bsFlatten (Link, 3, units) ;
  bsDestroy (Link) ;
  order = arrayCreate (50, BSunit) ;
  max = 0 ;
  for (i = 0 ; i < arrayMax(units) ; i += 3)
    { key = arr (units, i, BSunit).k ;
      x1 = arr (units, i + 1, BSunit).i ;
      x2 = arr (units, i + 2, BSunit).i ;
      if (key == contig)
	{ 
	  if (x1 < x2)
	    { array (order, max++, BSunit).k = key1 ;
	      array (order, max++, BSunit).i = max1 ;
	      array (order, max++, BSunit).k = key2 ;
	      array (order, max++, BSunit).i = max2 ;
	    }
	  else
	    { array (order, max++, BSunit).k = key2 ;
	      array (order, max++, BSunit).i = - max2 ;
	      array (order, max++, BSunit).k = key1 ;
	      array (order, max++, BSunit).i = - max1 ;
	    }
	}
      else
	{ array (order, max++, BSunit).k = key ;
	  array (order, max++, BSunit).i = x2 - x1 ;
	}
    }
  arrayDestroy (units) ;
  alignToolsAdjustLink (link, 0, order) ;
  arrayDestroy (order) ;

  /* after looking for contig in link */
  sessionGainWriteAccess () ;
  if (dnaSubClass (contig, &dnaKey))
    {
      lexAlias (&dnaKey, messprintf("%s.%s", buff, name(contig)), FALSE, FALSE) ;
      lexAlias (&contig, messprintf("%s.%s", buff, name(contig)), FALSE, FALSE) ;
    }

 abort:
  bsDestroy (Link) ;
}

/***************************************************************************************/
 
static KEY movedContig1 = 0 ;
static LOOK movedLook = 0 ;

static void contigFlip (int box)
{ BOOL done = FALSE ;
  SEG *seg ;
  int a1, a2  ;
  OBJ obj ;
  FMAPLOOKGET("contig flip") ;
  
  if ((seg = BOXSEG(box)) &&
       (obj = bsUpdate(look->seqKey)))
    { if (bsFindKey (obj, _Subsequence, seg->key) &&
	  bsGetData (obj, _bsRight, _Int, &a1) &&
	  bsGetData (obj, _bsRight, _Int, &a2))
	{ bsFindKey (obj, _Subsequence, seg->key) ;
	  bsAddData (obj, _bsRight, _Int, &a2) ;
	  bsAddData (obj, _bsRight, _Int, &a1) ;
	  done = TRUE ;
	}
      bsSave (obj) ;
    }
  if (done)
    { if (seg->x1 <= look->map->centre && look->map->centre <= seg->x2)
	look->map->centre = seg->x1 + seg->x2 - look->map->centre ;
      traceGraphDestroy () ;
      look->pleaseRecompute = TRUE ;
      fMapDraw (look, 0) ;
    }
}

static void contigFixExtend (int box, char action)
{ SEG *seg ;
  BSunit *u ;
  OBJ obj = 0 ;
  KEYSET ks = 0 ;
  KEY contig ;
  int i, j, dx, max ;
  Array units = 0, order = 0 ; 
  FMAPLOOKGET("fixConsensus") ;

  seg = BOXSEG(box) ;
  if (!seg ||
      movedContig1
      )
    return ;
  
  contig = seg->key ;
  
  switch (action)
    {
    case 'C': 
      fMapTraceCleanContig (contig) ;
      break ;
    case 'd': /* double */
      fMapTraceCleanContig (contig) ;
      j = 0 ;
      i = abiFixDoubleContig (look->seqKey, contig, &j) ;
      messout (" %d clips moved by an average of %d bases",
	               i, i > 0 ? j/i : 0 ) ;
      break ;
    case 'e': 
      fMapTraceCleanContig (contig) ;
      abiFixExtendContig (look->seqKey, contig) ;
      break ;
    case 'f': 
      dnaAlignFixContig (look->seqKey, contig) ; 
      break ;
    case 'a':
      i = 0 ;
      baseCallPatchContig (contig, &i) ;
      messout ("%d edits", i) ;
      goto done ;
      break ;
    case 'E': case 'G': case 'F':  case 'H': case 'A': 
      baseCallUnclipContig (contig, action, &dx) ; 
      dnaAlignFixContig (look->seqKey, contig) ;
      break ;
    case 'T':
      baseCallTileContig (contig, action, &dx) ; 
      dnaAlignFixContig (look->seqKey, contig) ;
      break ;
    case '6':
      { float xx, yy, z1 ; int yy1, x1, x2 ;
	graphBoxAt (graphEventX, graphEventY, &xx, &yy) ;
	z1 = MAP2GRAPH(look->map, seg->x1) ;
	if (z1 < topMargin) yy += topMargin - z1 ;	  
	yy1 = yy / look->map->mag ; /* coord now in bp in contig */
	yy1 += seg->x1 ;
	x1 = yy1 - 500 ; if (x1 < seg->x1) x1 = seg->x1 ;
	x2 = yy1 + 500 ; if (x2 >= seg->x2) x2 = seg->x2 - 1 ;
	look->zoneMin = 0 ;
	look->zoneMax = arrayMax(look->dna) ;
	dnaAlignFindRepeat (look->dna, look->dna, look->colors, x1, x2, 6) ;
	mapColSetByName ("Summary bar", TRUE) ;
	fMapDraw (look, 0) ; /* NO recompute */
	return ;
      }
    case 't':
      { float xx, yy, z1 ; int yy1, x1, x2 ;
	graphBoxAt (graphEventX, graphEventY, &xx, &yy) ;
	z1 = MAP2GRAPH(look->map, seg->x1) ;
	if (z1 < topMargin) yy += topMargin - z1 ;	  
	yy1 = yy / look->map->mag ; /* coord now in bp in contig */
	yy1 += seg->x1 ;
	x1 = yy1 - 500 ; if (x1 < seg->x1) x1 = seg->x1 ;
	x2 = yy1 + 500 ; if (x2 >= seg->x2) x2 = seg->x2 - 1 ;
	look->zoneMin = 0 ;
	look->zoneMax = arrayMax(look->dna) ;
	dnaAlignFindRepeat (look->dna, look->dna, look->colors, x1, x2, 12) ;
	mapColSetByName ("Summary bar", TRUE) ;
	fMapDraw (look, 0) ; /* NO recompute */
	return ;
      }
    case 'c': 
      { float xx, yy, z1 ; int yy1 ;
	graphBoxAt (graphEventX, graphEventY, &xx, &yy) ;
	z1 = MAP2GRAPH(look->map, seg->x1) ;
	if (z1 < topMargin) yy += topMargin - z1 ;	  
	yy1 = yy / look->map->mag ;
	fMapTraceCleanContig (contig) ;
	fMapTraceCutContig (look, contig, yy1, 'c') ; goto done ;
      }
    case 'u':
      { float xx, yy, z1 ; int yy1 ;
	graphBoxAt (graphEventX, graphEventY, &xx, &yy) ;
	z1 = MAP2GRAPH(look->map, seg->x1) ;
	if (z1 < topMargin) yy += topMargin - z1 ;	  
	yy1 = yy / look->map->mag ;
	fMapTraceCutContig (look, contig, yy1, 'u') ; goto done ;
      }
    case 'r': 
      fMapTraceCleanContig (contig) ;
      fMapTraceReAssembleContig (look, contig) ;
      return ;
    case 'R': 
      fMapTraceRemoveContig (look, contig) ;
      return ;
    case 'p': 
      fMapTraceCleanContig (contig) ;
      ks = alignToolsPurifyContig (look->seqKey, contig) ;
      keySetDestroy (ks) ;
      goto done ;
    default: return ;
    }

  obj = bsCreate (look->seqKey) ;
  if (!obj) goto done ;

  units = arrayCreate (90, BSunit) ; 
  if (bsFindTag (obj, _Subsequence) && 
      bsFlatten (obj, 3, units))
    { order = arrayCreate (90, BSunit) ;
      max = 0 ;
      for (i = 0 ; i < arrayMax(units) ; i += 3)
	{ u = arrp (units, i, BSunit) ;
	  array (order, max++, BSunit).k = u[0].k ;
	  array (order, max++, BSunit).i = u[2].i - u[1].i ;
	}
      alignToolsAdjustLink (look->seqKey, 0, order) ;
      arrayDestroy (order) ;
    }
 done:
  arrayDestroy (units) ;
  bsDestroy (obj) ;
  
  traceGraphDestroy () ;
  look->pleaseRecompute = TRUE ;
  fMapDraw (look, 0) ;
}
  
static void contigExtend (int box)
{ contigFixExtend (box, 'e') ;
}

static void contigDouble (int box)
{ contigFixExtend (box, 'd') ;
}

static void contigFixConsensus (int box)
{ contigFixExtend (box, 'f') ;
}

static void contigAutoEdit (int box)
{ contigFixExtend (box, 'a') ;
}

static void contigUnclipHand (int box)
{ contigFixExtend (box, 'H') ;
}

static void contigUnclipTile (int box)
{ contigFixExtend (box, 'T') ;
}

static void contigUnclipExcellent (int box)
{ contigFixExtend (box, 'E') ;
}

static void contigUnclipGood (int box)
{ contigFixExtend (box, 'G') ;
}

static void contigUnclipFair (int box)
{ contigFixExtend (box, 'F') ;
}

static void contigReAssemble (int box)
{ contigFixExtend (box, 'r') ;
}

static void contigRemove (int box)
{ contigFixExtend (box, 'R') ;
}

static void contigCleanUp (int box)
{ contigFixExtend (box, 'C') ;
}

static void contigFind6Repeat (int box)
{ contigFixExtend (box, '6') ;
}

static void contigFind12Repeat (int box)
{ contigFixExtend (box, 'T') ;
}

static void contigCut (int box)
{ contigFixExtend (box, 'c') ;
}

static void contigUnJoin (int box)
{ contigFixExtend (box, 'u') ;
}

static void  contigMoveCancelled (void)
{ graphUnMessage () ;
  messStatus ("assembling") ;
  movedContig1 = 0 ;
  movedLook = 0 ;
  displayUnBlock () ;
}

static void contigDoMove (KEY key2)
{ OBJ Link ;
  KEY key1 = movedContig1, key, link = movedLook->seqKey ;
  int 
    i, j, x1, x2, p1 = 0, p2 = 0 ;
  Array units = 0, bilan = 0 ; 
  LOOK look = movedLook ;

  graphActivate (look->graph) ;
  
  contigMoveCancelled () ;
  if (!key2 ||
      key2 == key1)
    return ;
	    
  Link = bsCreate (link) ;
  units = arrayCreate (50, BSunit) ;
  if (!Link ||
      !bsFindKey (Link, _Subsequence, key2) ||
      !bsFindTag (Link, _Subsequence) ||
      !bsFlatten (Link, 3, units))
    { bsDestroy (Link) ;
      arrayDestroy (units) ;
      return ;
    }
  bsDestroy (Link) ;
  for (i = 0 ; i < arrayMax(units) ; i += 3)
    { key = arr (units, i, BSunit).k ;
      if (key == key1)
	{ p1 = arr (units, i + 1, BSunit).i ;
	  p2 = arr (units, i + 2, BSunit).i ;
	  break ;
	}
    }
  bilan = arrayCreate (50, BSunit) ;
  j = 0 ;
  if (p1 || p2)
    for (i = 0 ; i < arrayMax(units) ; i += 3)
      {
	key = arr (units, i, BSunit).k ;
	x1 = arr (units, i + 1, BSunit).i ;
	x2 = arr (units, i + 2, BSunit).i ;
	if (key == key1)
	  continue ;
	array (bilan, j++, BSunit).k = key ;
	array (bilan, j++, BSunit).i = x2 - x1 ;
	if (key == key2)
	  { array (bilan, j++, BSunit).k = key1 ;
	  array (bilan, j++, BSunit).i = p2 - p1 ;
	  }
      }
  arrayDestroy (units) ;
  alignToolsAdjustLink (link, 0, bilan) ;
  arrayDestroy (bilan) ;
  traceGraphDestroy () ;
  look->pleaseRecompute = TRUE ;
  fMapDraw (look, 0) ;
}

static void contigMove (int box)
{ SEG *seg ;
  FMAPLOOKGET("contig flip") ;
  
  seg = BOXSEG(box) ;
  if (!seg ||
      movedContig1
      )
    { contigMoveCancelled () ;
      return ;
    }
  movedContig1 = seg->key ;
  movedLook = look ;
  displayBlock (contigDoMove, 
		"Pick in same window a contig under which you want to move") ;
  graphRegister (MESSAGE_DESTROY, contigMoveCancelled) ;
}


static void contigMoveToBottom (int box)
{ SEG *seg ;
  KEY lastKey = 0 ;
  OBJ Link = 0 ;
  FMAPLOOKGET("contig flip") ;
  
  seg = BOXSEG(box) ;
  if (!seg ||
      movedContig1 ||
      !(Link = bsCreate (look->seqKey)))
    return ;
  
  if(bsGetKey (Link, _Subsequence, &lastKey))
    while (bsGetKey (Link, _bsDown, &lastKey)) ;
  bsDestroy (Link) ;
  
  if (!lastKey ||
      lastKey == seg->key)
    return ;
  
  movedContig1 = seg->key ;
  movedLook = look ;
  contigDoMove (lastKey) ;
}


static void contigDoPair (KEY key2)
{ KEY key1 = movedContig1 ;
  
  contigMoveCancelled () ;
  if (!key2 ||
      key2 == key1)
    return ;

  contigDoMakePair (key2, key1) ;
}

static void contigMakePair (int box)
{
  SEG *seg ;
  FMAPLOOKGET("contig makePair") ;
  
  seg = BOXSEG(box) ;
  if (!seg ||
      movedContig1
      )
    { contigMoveCancelled () ;
      return ;
    }
  movedContig1 = seg->key ;
  movedLook = look ;
  displayPreserve() ;
  displayBlock (contigDoPair, 
		"Pick in any window a second contig as your target") ;

}

/*****************************************/

static void contigDoForcePair (LOOK look, KEY target, KEY petit)
{
  KEYSET ks = queryKey (petit, ">Assembled_from") ;

  if (!keySetMax(ks))
    keySet (ks, 0) = petit ;
  fMapTraceCleanContig (target) ;
  if (!dnaAlignContigForcePair (look->seqKey, target, ks))
      messout ("Sorry, no fit") ;
  else
    {  OBJ obj ;
       if ((obj = bsUpdate (petit)))
	 bsKill (obj) ;  

       look->pleaseRecompute = TRUE ;
       traceGraphDestroy () ;
       fMapDraw (look, 0) ;
     }
}

/*****************************************/

static void contigDoForce (KEY key2)
{ KEY key1 = movedContig1 ;
  LOOK look = movedLook ;
  
  contigMoveCancelled () ;
  if (!key2 ||
      key2 == key1)
    return ;

  contigDoForcePair (look, key2, key1) ;
}

static void contigForce (int box)
{
  SEG *seg ;
  FMAPLOOKGET("contig force") ;
  
  displayPreserve () ;
  seg = BOXSEG(box) ;
  if (!seg ||
      movedContig1
      )
    { contigMoveCancelled () ;
      return ;
    }
  movedContig1 = seg->key ;
  movedLook = look ;
  displayBlock (contigDoForce, 
		"Pick a second contig as your target") ;

}

/*****************************************/
#ifdef JUNK
static void contigDoMakeDotter (KEY key2)
{ KEY key1 = movedContig1 ;
  Array dna1, dna2 ;
  KEY d1, d2 ;
  contigMoveCancelled () ;
  if (!key2 ||
      key2 == key1)
    return ;

  if (key1 && dnaSubClass (key1, &d1) &&
      key2 && dnaSubClass (key2, &d2))
    { dna1 = dnaGet (d1) ; dna2 = dnaGet (d2) ;
      if (dna1 && dna2)
	{ dnaDecodeArray (dna1) ;   dnaDecodeArray (dna2) ;
/*
	  dotter ('N', 0, name(key1), 
		  arrp(dna1, 0, char), 0, 
		  name(key2), 
		  arrp (dna2, 0, char), 0,0,0,0,0,0,0,0,0) ;
*/

	}
    }
}

static void contigMakeDotter (int box)
{
  SEG *seg ;
  FMAPLOOKGET("contig flip") ;
  
  seg = BOXSEG(box) ;
  if (!seg ||
      movedContig1
      )
    { contigMoveCancelled () ;
      return ;
    }
  movedContig1 = seg->key ;
  movedLook = look ;
  displayBlock (contigDoMakeDotter, 
		"Pick a second contig") ;

}
#endif
static void contigDoMakeCompare (LOOK look, KEY key1, KEY key2)
{ Array dna1 = 0, dna2 = 0 ;
  KEY dnaKey1 = 0, dnaKey2 = 0 ;
  OBJ obj ;

  obj = bsCreate (key1) ;
  if (obj)
    { bsGetKey (obj, _DNA, &dnaKey1) ;
      bsDestroy (obj) ;
    }

  obj = bsCreate (key2) ;
  if (obj)
    { bsGetKey (obj, _DNA, &dnaKey2) ;
      bsDestroy (obj) ;
    }

  if (dnaKey1)
    dna1 = dnaGet (dnaKey1) ;
  if (dnaKey2)
    dna2 = dnaGet (dnaKey2) ;

  if (!dna2)
    messout("PLease pick a contig") ;
  if (dna1 && dna2)
    { messStatus ("Comparing") ;
      if (!dnaAlignCompare (key1, key2, dna1, dna2))
	messout ("Sorry, i cannot align these 2 contigs ") ;
      else
	{ look->pleaseRecompute = TRUE ;
	  fMapDraw (look, 0) ;
	}
    }
  arrayDestroy(dna1) ;
  arrayDestroy(dna2) ;
}
 
static void contigDoCompare (KEY key2)
{ KEY key1 = movedContig1 ;
  LOOK look = movedLook ;

  graphActivate (look->graph) ;

  contigMoveCancelled () ;
  if (!key2 ||
      key2 == key1)
    return ;
  
  contigDoMakeCompare (look, key1, key2) ;
}

static void contigCompare (int box)
{
  SEG *seg ;
  FMAPLOOKGET("contigCompare") ;
  
  displayPreserve () ;
  seg = BOXSEG(box) ;
  if (!seg ||
      movedContig1
      )
    { contigMoveCancelled () ;
      return ;
    }
  movedContig1 = seg->key ;
  movedLook = look ;
  displayBlock (contigDoCompare, 
		"Pick a second contig as your target") ;

}

static void contigWhoAmI (int box)
{
  SEG *seg ;
  FMAPLOOKGET("contigWhoAmI") ;
  
  seg = BOXSEG(box) ;
  if (seg && seg->key)
    display (seg->key, 0, TREE) ;
}

static MENUOPT contigMenu[] = { 
 { (GraphFunc)contigWhoAmI, "Who am I ?"},
 { menuSpacer,""},
 { (GraphFunc)contigCut, "Cut Contig here"},
 { (GraphFunc)contigUnJoin, "Unjoin Contig here"},
 { (GraphFunc)contigFlip, "Flip"},
 { (GraphFunc)contigMove,"Move under"}, 
 { (GraphFunc)contigMoveToBottom,"Move to Bottom"},
 { (GraphFunc)contigRemove, "Remove Contig"},
 { (GraphFunc)contigCleanUp, "Clean Previous"},
 {  menuSpacer,""},
 { (GraphFunc)contigMakePair, "Join to"},
 { (GraphFunc)contigReAssemble, "Re-assemble"},
 { (GraphFunc)contigCompare, "Compare"},
   { (GraphFunc)contigForce, "Assemble into"},
/*  {(GraphFunc)contigMakeDotter, "Dot Plot with",  works, but very slow on whole contigs}, */
 {  menuSpacer,""},
 { (GraphFunc)contigUnclipExcellent, "Clip on Excellent"},
 { (GraphFunc)contigUnclipGood, "Clip on Good"},
 { (GraphFunc)contigUnclipFair, "Clip on Fair"},
 { (GraphFunc)contigUnclipTile, "Clip to Tile"},
 { (GraphFunc)contigUnclipHand, "Restore Hand Clips"},
   { (GraphFunc)contigExtend, "Extend Extremities"},
   { (GraphFunc)contigDouble, "Consolidate"}, 
 { menuSpacer,""},

/*  (GraphFunc)contigUnclipAll, "Unclip To Max", */
/*  (GraphFunc)contigPurify, "Split", */
 { (GraphFunc)contigFixConsensus, "Fix"},
 { (GraphFunc)contigAutoEdit, "Auto-edit"},
 { (GraphFunc)contigFind6Repeat, "Looks like (common 6-mers)"},
 { (GraphFunc)contigFind12Repeat, "Looks like (common 12-mers)"},
   { 0, 0} 
} ;

static void fMapTraceCleanKeepContig (KEY contig, KEYSET keep)
{ 
  OBJ obj = bsUpdate (contig) ;
  Array aa, aa1 ;
  int i, ii, j, dummy ;
  KEY key, target = _Aligned ;

  if (!obj)
    return ;

  aa = arrayCreate (100, BSunit) ;
  
lao:
  if (bsGetArray (obj, _Aligned, aa, 5))
    { aa1 = arrayCreate (100, BSunit) ;
      bsFindTag (obj, _Aligned) ;
      bsRemove (obj) ;
      ii = 0 ;
      for (i = 0 ; i < arrayMax(aa) ; i += 5)
	{ 
	  key = array (aa, i, BSunit).k ;
	  if (keySetFind (keep, key, &dummy))
	    for (j = 0 ; j < 5 ; j++)
	      array(aa1, ii++, BSunit) = array(aa, i+j, BSunit) ;
	}
      if (ii > 0)
	bsAddArray (obj, _Aligned, aa1, 5) ;
      arrayDestroy (aa1) ;
    }
  if (target == _Aligned) /* fake enumerated loop */
    { target = _Previous_contig ; goto lao ; }

  arrayDestroy (aa) ;
  bsSave (obj) ;
}

static void fMapTraceCleanContig (KEY contig)
{ OBJ obj = bsUpdate (contig) ;

  if (obj)
    {  
      if (bsFindTag (obj, _Aligned))
	bsRemove (obj) ;
      if (bsFindTag (obj, _Previous_contig))
	bsRemove (obj) ;

      bsSave (obj) ;
    } 
}

static void fMapTraceCleanAssembly (KEY aa, KEYSET keep)
{ KEYSET ks = queryKey (aa, "> Subsequence") ;
  int i = keySetMax (ks) ;
  
  while (i--)
    if (keep)
      fMapTraceCleanKeepContig (keySet (ks, i), keep) ;
    else
      fMapTraceCleanContig (keySet (ks, i)) ;
  keySetDestroy (ks) ;
}
   
static void fMapTraceRetitle (LOOK look)
{ KEY dummy, link = look->seqKey ;
  KEYSET contigs = queryKey (link, ">Subsequence") ;
  int nc = keySetMax (contigs) ;
  KEYSET reads = query (contigs, "> Assembled_from") ;
  int nr = keySetMax (reads) ;
  int nn = 0, nbc = 0, i = nc, x ;
  OBJ obj = 0 ;
  
  while (i--)
    if ((obj = bsCreate (keySet (contigs, i))))
      { if (bsGetKey (obj, _DNA, &dummy) &&
	    bsGetData (obj, _bsRight, _Int, &x))
	  { nn += x ;
	    if (x > 3000) nbc++ ;
	  }
	bsDestroy (obj) ;
      }
  keySetDestroy (contigs) ;
  keySetDestroy (reads) ;
  
  graphRetitle 
    (messprintf ("%s %s : %d Contigs (%d > 3 kb), %d reads, %d bases",
		mainCloneName(0), name(link), nc, nbc, nr, nn)) ;
  pickRememberDisplaySize ("FMAP") ;
}

typedef struct { KEY key ; int nr, nb ;} REPORT ;

static int reportOrder (const void *a, const void *b)
{ int x = ((const REPORT*)a)->nb,  y = ((const REPORT*)b)->nb ;
  return y - x ;
}

static void fMapTraceReport (void)
{ KEY key, dummy, link ;
  KEYSET contigs ;
  int line, x, i, nc, nrtotal = 0, a1 = 0, a2 = 0 ; 
  KEYSET reads ;
  OBJ obj = 0 ;
  Array aa ;
  REPORT *r ;
  Graph old = graphActive () ;
  FMAPLOOKGET("fMapTraceReport") ;

  link = look->seqKey ;
  contigs = queryKey (link, ">Subsequence") ;
  i = nc = keySetMax (contigs) ;
  aa = arrayCreate (nc, REPORT) ;

  while (i--)
    if ((obj = bsCreate (key = keySet (contigs, i))))
      { r = arrayp (aa, i, REPORT) ;
	if (bsGetKey (obj, _DNA, &dummy) &&
	    bsGetData (obj, _bsRight, _Int, &x))
	  r->nb = x ;
	bsDestroy (obj) ;
	reads = queryKey (key, "> Assembled_from") ;
	r->nr = keySetMax (reads) ;
	r->key = key ;
	keySetDestroy (reads) ;
      }
  keySetDestroy (contigs) ;
  
  arraySort (aa, reportOrder) ;
  
  graphCreate (TEXT_FULL_SCROLL, 
	       messprintf("%s %s", mainCloneName(0), name(look->seqKey)),
	       .2, .2, .6,.6) ;
  
  graphText ("    Length,     cumul,    # reads,    cumul     Contig", 1, 1) ;
  for (line = 4, i = 0, a1 = 0, a2 = 0 ; i < arrayMax (aa) ; i++, line++)
    { r = arrayp (aa, i, REPORT) ;
      graphText (messprintf ("%8d", r->nb), 2, line) ; a1 += r->nb ;
      graphText (messprintf ("%8d", a1), 12, line) ;
      graphText (messprintf ("%8d", r->nr), 22, line) ; a2 += r->nr ;
      graphText (messprintf ("%8d", a2), 32, line) ;
      graphText (name(r->key), 42, line) ;
      nrtotal += r->nr ;
    }
  graphText (messprintf (" %d", nrtotal), 10, 2) ;
  graphTextBounds (50, line + 3) ;
  graphRedraw () ;
  graphActivate (old) ;
  arrayDestroy (aa) ;
}

static void fMapTraceToggleTrack (void)
{ 
  FMAPLOOKGET("newTracking") ;
  
  switch (newErrorTracking)
    {
    case 0: 
    case 1: newErrorTracking = 2 ; break ;
    case 2: newErrorTracking = 1 ; break ;
    }
  traceGraphDestroy () ;
  look->pleaseRecompute = TRUE ;
  fMapDraw (look, 0) ;
}

#ifdef JUNK
static void fMapToggleKb (void)
{
  FMAPLOOKGET("toggleKb") ;
  
  look->flag ^= FLAG_HIDE_SMALL_CONTIGS ;
  look->pleaseRecompute = TRUE ;
  traceGraphDestroy () ;
  fMapDraw (look, 0) ;
}
#endif

static void fMapTraceAddWall (void)
{
  SEG *seg ;
  int i, centre, x1, x2 ;
  BOOL done = FALSE ;
  OBJ Seq = 0 ;
  FMAPLOOKGET("fMapTraceAddWall") ;

  centre = look->map->centre ;
  cDNAAlignInit() ;

  for (i = 1 ; i < arrayMax(look->segs) ; ++i)
    { 
      seg = arrp(look->segs,i,SEG) ;
    
      if ((seg->type | 0x1) != SEQUENCE_UP)
	continue ;
      if (seg->x1 > centre || seg->x2 < centre)
	continue ;
      if (!keyFindTag (seg->key, _Genomic))
	continue ;
      if (!(Seq = bsUpdate(seg->key)))
	continue ;
      if (seg->type == SEQUENCE)
	{ x1 = centre - seg->x1 + 1 ; x2 = x1 + 10 ; }
      else /* SEQUENCE UP */
	{ x1 = -centre + seg->x2 + 1 ; x2 = x1 - 10 ; }
      if (bsAddData (Seq, _Gene_wall, _Int, &x1))
	bsAddData (Seq, _bsRight, _Int, &x2) ;
      bsSave (Seq) ;
      done = TRUE ;
      break ;
    }
  if (done)
    fMapDraw (look, 0) ;	
}

/***********************************************************************/
static KEYSET selectedReads[3] ;

static void fMapColorQuery (void)
{
  FMAPLOOKGET("fMapColorQuery") ;
  
  look->flag ^= FLAG_COLOR_CONTIGS ;

  fMapDraw (look, 0) ;
}

void fMapColorSeg (SEG *seg)
{
  int col, j ;

  for (col = 1 ; col < 3 ; col++)
    if (keySetExists(selectedReads[col]) && 
      keySetFind(selectedReads[col], seg->key, &j))
      seg->data.i |= (col & 0x3) << 22 ;
}

int fMapQueryColor (KEY key)
{
  int col = 0, i, i1 = 0 , j ;

  for (i = 1 ; i < 3 ; i++)
    if (keySetExists(selectedReads[i]) && 
      keySetFind(selectedReads[i], key, &j))
      i1 = i ;
  switch (i1)
    {
      case 1: col = YELLOW ; break ;
      case 2: col = LIGHTBLUE ; break ;
      case 3: col = LIGHTGREEN ; break ;
    }
  return col ;
}

static void fMapUnQuery (LOOK look)
{
  SEG *seg ;
  int i ;
  
  for (i = 1 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;

      if (seg->type != VIRTUAL_SUB_SEQUENCE_TAG &&
	  seg->type != VIRTUAL_SUB_SEQUENCE_TAG_UP)
	continue ;
      seg->data.i &= 0xff3fffff ; /* zeroes bits 22 and 23 */
    }
  i = 3 ;
  while (i--) keySetDestroy (selectedReads[i]) ;
  look->flag |= FLAG_COLOR_CONTIGS ;
  fMapDraw (look, 0) ;
}

static void fMapQuery (KEY col, int box)
{ 
  char *cp = 0 ;
  SEG *seg ; 
  int i, j ;
  KEYSET ks ;
  FMAPLOOKGET("fMapQuery") ;
  
  if (!col)
    fMapUnQuery (look) ;
  else if (col < 3)
    {
      if (!messPrompt("Give a query to select some sequence","IS *","t"))
	return ;
      cp = strnew (freeword(), 0) ;
      keySetDestroy (selectedReads[col]) ;
      ks = keySetCreate () ; j = 0 ;
      for (i = 1 ; i < arrayMax(look->segs) ; i++)
	{
	  seg = arrp (look->segs, i, SEG) ;
	  keySet(ks, j++) = seg->key ;
	}
      keySetSort (ks) ; keySetCompress (ks) ;
      selectedReads[col] = query (ks, cp) ;
      keySetDestroy (ks) ;
      messfree(cp) ;
      fMapDraw (look, 0) ;
    }
}

static FREEOPT fMapQueryMenu[] = 
{ {3,"Colors"},
  {1, "Yellow"},
  {2, "Blue"},
  {3, "Green"},
  {0, "UnColor"},
} ;

/***********************************************************************/

static Array getPhrapQuality (KEY key)
{ KEY pk ;
  Array aa = 0 ;
  OBJ obj = bsCreate (key) ;

  if (obj &&
      bsGetKey (obj, _Quality, &pk))
    aa = arrayGet (pk, char, "c" ) ;
  bsDestroy (obj) ;
  return aa ;
}

/***********************************************************************/

void fMapShowPhrapQuality (LOOK look, float *offset)
{ 
  Array aa = 0 ;
  float y = 0, y1 = 0, y2 = 0 ;
  int   ii, jj, ibox, x, x2, new, y1col, dx, dx0 ;
  SEG   *seg ;
  BOOL found = FALSE ;
  unsigned  char *ucp ;
  int mycol[] = { BLACK, RED, ORANGE, BLUE, GREEN, YELLOW } ; 
  int ncolors = 6,    /* nb of col in above enum */
    colorBinSize = 100/(ncolors -1) ;
  SEQINFO *sinf ;
  
  if (!(look->flag & FLAG_VIRTUAL_ERRORS))
    return ;
  for (ii = 1 ; ii < arrayMax(look->segs) ; ii++)
    {
      seg = arrp(look->segs,ii,SEG) ;
      
      if ( (!seg->data.i) ||
	   (seg->type != SEQUENCE &&
	    seg->type != SEQUENCE_UP)
	   )
	continue ;
      sinf = arrayp (look->seqInfo, seg->data.i, SEQINFO) ;
      if ( !(sinf->flags & SEQ_VIRTUAL_ERRORS))
	continue ;
      
      y1 = MAP2GRAPH(look->map, seg->x1) ;
      y2 = MAP2GRAPH(look->map, seg->x2 - 1) ;
      if (y1 >= mapGraphHeight || y2 <= topMargin)
	continue ;

      if (!(aa = getPhrapQuality (seg->key) ))
	continue ;
      if (arrayMax(aa) != 1 + seg->x2 - seg->x1)
	{
	  messout ("max = %d  != x2-x1 = %d", arrayMax(aa),seg->x2 - seg->x1) ;
	  arrayDestroy (aa) ;
	  continue;
	}

	/* test, alternating  linear ramps 
      { int i, j ;

	j = arrayMax(aa) ;
	for (i=0; i < j ; i++)
	  array (aa,i,unsigned char) = ii %2 ? (100 * i) / j : (100 * (j - i)) / j ;
      }
      */


	

      found = TRUE ;

      array(look->boxIndex,ibox=graphBoxStart(),int) = ii ;
      if (y1 < topMargin) y1 = topMargin ;
      if (y2 > mapGraphHeight) y2 = mapGraphHeight ;
      graphRectangle(*offset, y1 - .1, *offset + 1.3 , y2 + .1) ;
      x = GRAPH2MAP (look->map, y1) ; new = 0 ; 
      x2 = GRAPH2MAP (look->map, y2) ;
      if (x2 > x)
	{ jj = x2 - x ;
	  ucp = arrp (aa, x - seg->x1, unsigned char) -1 ;
	}
      else 
	{ jj = 0 ;
	  ucp = arrp (aa, 0, unsigned char) -1 ;
	}
      /* the idea is to proviledge the worst color at small zoom */
      /* dx is the number of steps before changing color */
      dx0 =  GRAPH2MAP (look->map, y1 + .2) - GRAPH2MAP (look->map, y1) ;
      if (dx0 < 1) dx0 = 1 ;
        
      dx = dx0 ; y = MAP2GRAPH(look->map, x) ; y1col =  ACEDB_MAXINT ;
      if (jj) while (ucp++, x++, jj--)
	{
	  y = MAP2GRAPH(look->map, x) ;  
	  new = *ucp ? 1 + *ucp / colorBinSize : 0 ;  if (new >= ncolors) new = ncolors - 1 ;  
	  if (!jj || !dx--) 
	    {
	      ibox = graphBoxStart () ;
	      graphRectangle (*offset + .1, y1, *offset + 1.2 , y) ;
	      graphBoxEnd () ;
	      graphBoxDraw (ibox, mycol[y1col], mycol[y1col]) ;
	      y1 = y  ; y1col = new ; dx = dx0 ;
	    }
	  else
	    {
	      if (new < y1col) y1col = new ; /* select worst color */
	    }
	}

      arrayDestroy (aa) ;
      graphBoxEnd() ;

     }
  *offset += found ? 2 : 0 ;
}

/***********************************************************************/

void fMapShowContig (LOOK look, float *offset)
{
  float y1, y2 ;
  int   ii, j, ibox, mxName = 0 , x1, x2 ;
  SEG   *seg ;
  BOOL found = FALSE ;
  char *cp ;
  SEQINFO *sinf ;

  if (!isGifDisplay)
    {
      if (mapColSetByName ("Summary bar", -1))
	graphButton ("Clear", fMapClear, 1, 4.8) ;
      ibox = graphButton ("Show...", fMapColorQuery, 8, 4.8) ;
      graphBoxFreeMenu (ibox,  fMapQuery, fMapQueryMenu) ;
      graphButton ("Arrows", fMapToggleTraceAssembly, 12, 4.8) ;
      /*   graphButton ("< 1kb", fMapToggleKb, *offset + 12, 4.8) ; */
      graphButton ("Quality", fMapTraceAllErrors, 22, 4.8) ;
      graphButton ("Report", fMapTraceReport, 32, 4.8) ;
      if (newErrorTracking)
	{ 
	  if (newErrorTracking == 2)
	    graphButton ("New Track On", fMapTraceToggleTrack, 41, 4.8) ;
	  else
	    graphButton ("New Track Off", fMapTraceToggleTrack, 41, 4.8) ;
	}
      graphButton ("DNA", fMapToggleDna, 53, 4.8) ;
      if (mapColSetByName ("DNA Sequence", -1))
	graphButton ("Variations", fMapTraceToggleAssemblyDna, 58, 4.8) ;
    }

  for (ii = 1 ; ii < arrayMax(look->segs) ; ii++)
    { seg = arrp(look->segs,ii,SEG) ;

      if ( (!seg->data.i) ||
	  (seg->type != SEQUENCE &&
	   seg->type != SEQUENCE_UP)
	  )
	continue ;
      sinf = arrayp (look->seqInfo, seg->data.i, SEQINFO) ;
      if ( !(sinf->flags & SEQ_VIRTUAL_ERRORS))
	continue ;
      found = TRUE ;

      y1 = MAP2GRAPH(look->map, seg->x1) ;
      y2 = MAP2GRAPH(look->map, seg->x2 - 1) ;
      if (y1 >= mapGraphHeight || y2 <= topMargin)
	continue ;
      array(look->boxIndex,ibox=graphBoxStart(),int) = ii ;
      if (y1 < topMargin) y1 = topMargin ;
      if (y2 > mapGraphHeight) y2 = mapGraphHeight ;
      x1 = GRAPH2MAP (look->map, y1) ;
      x2 = GRAPH2MAP (look->map, y2) ;
      graphRectangle(*offset, y1 - .1, *offset + 1.3 , y2 + .1) ;
      if (look->flag & FLAG_VIRTUAL_ERRORS)
	fMapTraceColorReport (look, *offset + .65, x1, x2) ;
      graphBoxEnd() ;
      if (seg->type == SEQUENCE)
	{ graphLine (*offset - .65, y2 - .5, *offset + .65 , y2 + .7) ;
	  graphLine (*offset +1.95, y2 - .5, *offset + .65 , y2 + .7) ;
	}
      else
	{ graphLine (*offset - .65, y1 + .5, *offset + .65 , y1 - .7) ;
	  graphLine (*offset +1.95, y1 + .5, *offset + .65 , y1 - .7) ;
	}
      if (!(look->flag & FLAG_VIRTUAL_ERRORS))
	graphBoxDraw(ibox, BLACK , LIGHTBLUE) ; /* default color */
      graphBoxMenu (ibox, contigMenu) ;
      cp = name(seg->key) ; j = strlen (cp) ; cp += j ; 
      while (*cp != '.' && j--) cp-- ;
      graphText (cp, *offset + 1.6, y1 + .8) ;
      j = strlen (cp) ;
      if (j > mxName) mxName = j ;
    }
  if (found)
    {
      *offset += 2 + mxName ;
      fMapShowPhrapQuality (look, offset) ;
      fMapTraceRetitle (look) ;
    }
}

/***************************************************************************************/
/***************************************************************************************/
 
 /* public, called from multitrace editor */
void fMapTraceForget (KEY key)
{ Array a ;
  char *vp ;
  void *wp ;
  Associator all = 0 ;

  if (!x1Look  || x1Look->magic != &FMAPLOOK_MAGIC)
    return ;
  all = x1Look->virtualErrors ;
  if (!assExists(all))
    return ;

    /* direct trace */
  vp = (char*)0 +  key ;
  if (vp && assFind(all, vp, &wp))
    { a = (Array) wp ;
      if (arrayExists(a))
	  arrayDestroy(a) ;
      assRemove (all, vp) ;
    }
   /* reverse trace */
  vp = (char *)0 + (1 << 23) + key ;
  if (vp && assFind(all, vp, &wp))
    { a = (Array) wp ;
      if (arrayExists(a))
	  arrayDestroy(a) ;
      assRemove (all, vp) ;
    }
}  

/***************************************************************************************/

void fMapTraceForgetContig (KEY contig)
{ KEYSET ks ;
  KEY *keyp ;
  int i ;

  ks = queryKey (contig, ">Assembled_from") ;
  i = keySetMax (ks) ;
  if (i)
    { keyp = arrp (ks, 0, KEY) - 1 ;
      while (keyp++, i--)
	fMapTraceForget (*keyp) ;
    }
  keySetDestroy (ks) ;
}
/***************************************************************************************/
/***************************************************************************************/

 /* public, called from multitrace editor */
void fMapTraceSuppress (KEY key, KEY tag, BOOL keep)
{
  SEG *seg1 ;
  KEY seq = 0, dna ;
  OBJ obj, Seq ;
  KEYSET aa ; int ii, dummy, x1, x2, a1, a2 ;
  FMAPLOOKGET("fMapTraceSuppress") ;
  
  if ((tag == _Assembled_into) && (obj = bsUpdate(key)))
    { KEY _By_hand ;
      lexaddkey ("By_hand", &_By_hand,0) ;
      bsAddTag (obj, _By_hand) ;
      bsSave (obj) ;
    }
  aa = queryKey (key, messprintf ("> %s", name(tag))) ;
  if (!keySetMax (aa))
    { keySetDestroy (aa) ;
      return ;
    }
  ii = keySetMax(look->segs) ;
  seg1 = arrp(look->segs, 0, SEG) - 1 ;
  while (seg1++, ii--)
    { if (!(seg1->data.k &  FLAG_VIRTUAL_ERRORS) ||
	  (seg1->type != SEQUENCE &&
	   seg1->type != SEQUENCE_UP)
	  )
	continue ;
      if (keySetFind (aa, seg1->key, &dummy))
	{ if (tag ==  _Assembled_into &&
	      (Seq = bsCreate (seg1->key)))
	  { seq = seg1->key ;
	    if (bsFindKey(Seq, _Assembled_from, key) &&
		bsGetData (Seq, _bsRight, _Int, &a1) &&
		bsGetData (Seq, _bsRight, _Int, &a2) &&
		bsGetData (Seq, _bsRight, _Int, &x1) &&
		bsGetData (Seq, _bsRight, _Int, &x2))
	    bsDestroy (Seq) ;
	  }
	    

	  if ((obj = bsUpdate(key)))
	    { if (bsFindKey (obj, tag, seg1->key))
	      bsRemove(obj) ;
	      bsSave (obj) ;
	    }
	  if ((tag == _Assembled_into) && dnaSubClass (key, &dna))
	    defCptAddSeqIn (look->seqKey, dna) ;
	  if (keep && seq)
	    { KEY newSeq, link = look->seqKey ;
	      int i, j = 1 ;
	      while (lexword2key (messprintf("%s.%d", name(link),j),
				   &newSeq, _VSequence)) j++ ;
	      lexaddkey (messprintf("%s.%d", name(link),j), &newSeq, _VSequence) ;
	      if ((Seq = bsUpdate (newSeq)))
		{ Array dna, newDna ;
		  KEY dnaKey, newDnaKey ;
		  bsAddKey (Seq, _Assembled_from, key) ;
		  a2 = x2 - x1 + 1 ; a1 = 1 ;
		  bsAddData (Seq, _bsRight, _Int, &a1) ;
		  bsAddData (Seq, _bsRight, _Int, &a2) ;
		  bsAddData (Seq, _bsRight, _Int, &x1) ;
		  bsAddData (Seq, _bsRight, _Int, &x2) ;
		  bsSave (Seq) ;
		  lexaddkey (name(newSeq), &newDnaKey, _VDNA) ;
		  newDna = arrayCreate (x2 - x1 + 2, char) ;
		  dnaSubClass(key, &dnaKey) ;
		  dna = dnaGet (dnaKey) ;

		  i = x2 - x1 + 1 ;
		  while (i--) array(newDna,i, char) = arr (dna, x1 + i, char) ;
		  dnaStoreDestroy (newDnaKey, newDna) ;
		  arrayDestroy (dna) ;
		  bsSave (Seq) ;
		}
		  /* need to put the contig in the assembly */
	      fMapTraceAddContig (look->seqKey, seq, newSeq, x2 - x1 + 1) ; 
	    }

	  goto ok ;
	}
    }
  keySetDestroy (aa) ;
  return ;
 ok:
  keySetDestroy (aa) ;
  look->pleaseRecompute = TRUE ;
  fMapDraw (look, 0) ;
}

/***************************************************************************************/
/***************************************************************************************/
typedef struct { int x1, x2, x3, a1, a2, a3 ; KEY  key, dnaKey ; } READ ;

static void contigDoMakePair (KEY contig1, KEY contig2) 
{ KEYSET ks ;
  int i ;
  OBJ obj ;
  KEY link, join, key ;
  BOOL done = FALSE ;
  FMAPLOOKGET("contigDoMakePair") ;

  fMapTraceCleanContig (contig1) ;
  fMapTraceCleanContig (contig2) ;

  messStatus ("Joining") ;
  join = dnaAlignAsmbPaire21 (contig1, contig2, 0, 0, 7, 0, 0) ;
  if (!join)
    join = dnaAlignAsmbPaire21 (contig1, contig2, 0, 0, 4, 0, 0) ;
  if (!join)
    join = dnaAlignAsmbPaire21 (contig1, contig2, 0, 0, -1, 0, 0) ;
  if (!join)
    { messout ("Sorry, I could not find a join") ;
      return ;
    }
  link = look->seqKey ;
  i = 0 ; ks = keySetCreate() ;
  if ((obj = bsCreate(link)))
    { if (bsGetKey (obj, _Subsequence, &key))
	do 
	  { 
	    if (key != contig1 && key != contig2)
	      keySet (ks, i++) = key ;
	    else if (!done)
	      { keySet (ks, i++) = join ;
		done = TRUE ;
	      }
	  } while (bsGetKey (obj, _bsDown, &key)) ;
      bsDestroy (obj) ;
    }
/* we want to keep the order */
/* ks = queryKey (link, ">Subsequence") ; */
  alignToolsAdjustLink (link, ks, 0);
  traceGraphDestroy () ;
  look->pleaseRecompute = TRUE ;
  fMapDraw (look, join) ;
}

/***************************************************************************************/
/***************************************************************************************/

static void fMapTraceNewName (KEY key, char *bb, char cc)
{ int i ;
  char *cp, *cr;
  char buf[256] ;

  strncpy (buf, name(key), 254) ; buf[255] = 0 ;
  cp = buf + strlen(buf) - 1 ;
  if (cc) { *(++cp) = cc ; *(cp + 1) = 0 ; }
  cr = buf ;
  while (*cr == '_') cr++ ;
  if (!strncmp (buf,"Link.", 5))
    cr += 5 ; /* so not to pass word link to Save command */
  if (*cp == '.') *cp-- = 0 ;
  while (*cp >= '0' && *cp <= '9') *cp-- = 0 ;
  i = 0 ; 
  while (++i)
    if (!lexword2key (messprintf("%s%d", buf, i),
		      &key, _VSequence) &&
	!lexword2key (messprintf("%s%d.", buf, i),
		      &key, _VSequence))
      break ;
  strncpy (bb, messprintf("%s%d", cr, i), 127) ;
}

static BOOL fMapTraceShouldRename(KEY key, char *buf)
{
  BOOL yes = FALSE ;
  char *tt = "This operation will modify the present assembly\n"
             "Press cancel to show the result in the same window\n"
             "Or give here the new name of the resutling assembly" ;
  yes = messPrompt(tt, buf,"wz") ;

  if (yes)
    { char *cp = freeword() ;
      strncpy (buf, cp, 127) ;
    }
  else
     strncpy (buf, name(key), 127) ;
  return yes ;
}
#ifdef JUNK

static void fMapTraceExplodeAll (void)
{  FMAPLOOKGET("fMapTraceExplodeAll") ;

 messout("BOOM !") ;return ;
}

static void fMapTraceRemoveTag (KEY key, KEY tag)
{ OBJ obj = bsUpdate (key) ;
  if (obj)
    { if (bsFindTag ( obj, tag))
	bsRemove (obj) ;
      bsSave (obj) ;
    }
}
#endif

static void fMapTraceReInsertLoners (void)
{ FMAPLOOKGET("fMapTraceReInsertLoners") ;

  messStatus ("Reinsert loners") ;
  if (!messQuery("You may want first to save this assembly,"
		 "do you want to continue"))
    return ;

  fMapTraceCleanAssembly (look->seqKey, 0) ;
  dnaAlignReInsertLoners (look->seqKey, 'r') ;
  traceGraphDestroy () ;
  look->pleaseRecompute = TRUE ;
  fMapDraw (look, 0) ;
}

static void fMapTraceReInsertSubclones (void)
{ FMAPLOOKGET("fMapTraceReInsertSuclones") ;

  fMapTraceCleanAssembly (look->seqKey, 0) ;
  dnaAlignReInsertLoners (look->seqKey, 's') ;
  traceGraphDestroy () ;
  look->pleaseRecompute = TRUE ;
  fMapDraw (look, 0) ;
}

static void fMapTraceGlobalCommand (char *com, BOOL cleanLeft)
{ 
  KEY key ;
  KEYSET ks ;
  FMAPLOOKGET("fMapTraceGlobalAction") ;

  ks = queryKey (look->seqKey, ">Subsequence") ;
  messStatus (com) ;
  key = defCptExecuteCommand (look->seqKey, 0, 0, com, 0) ;
  traceGraphDestroy () ;
  if (key)
    { 
      look->seqKey = key ;
      if (cleanLeft) 
	fMapTraceCleanAssembly (key, ks) ;
      look->pleaseRecompute = TRUE ;
      fMapDraw (look, 0) ;
      sessionDoSave (TRUE) ; /* save */
    }
  keySetDestroy (ks) ;
}

static void fMapTraceAddNewReads (void)
{ 
  char com[300] ;
  char nm [128] ;
  int n1, n2, n3 ; 
  KEYSET ks1 = 0, ks2 = 0, ks3 = 0 ;
  FMAPLOOKGET("fMapTraceAddNewReads") ;

  ks1 = query(0, "Find New_Read") ;
  n1 = keySetMax(ks1) ;
  keySetDestroy (ks1) ;
  if (!n1)
    { messout ("Zero new_read found, sorry") ;
      keySetDestroy (ks1) ;
      return ;
    }
  ks2 = queryKey (look->seqKey, "> subsequence ; >Assembled_from ") ;
  n2 = keySetMax(ks2) ;
  keySetDestroy (ks2) ;

  fMapTraceNewName (look->seqKey, nm, 0) ;
  if (fMapTraceShouldRename (look->seqKey, nm))
    {
      sprintf (com, "%s%s%s%s%s%s",
	       "NewSCF\n",
	       "Add >? New_Read\nGet 12\nSort 100\n",
	       "Assemble 5\nFix\nGet 40\nSort 100\n",
	       "Assemble 8\nJoin\nFix\nSave ", nm, 
	       "\nOrder_by_Size\n") ;
    }
  else
    sprintf (com, "%s%s%s%s%s%s",
	     "NewSCF\n",
	     "Add >? New_Read\nGet 12\nSort 100\n",
	     "Assemble 5\nFix\nGet 40\nSort 100\n",
	     "Assemble 8\nJoin\nFix", 
	     "\nOrder_by_Size\nRename ", name(look->seqKey)) ;

  fMapTraceGlobalCommand (com, TRUE) ;

  ks3 = queryKey (look->seqKey, "> subsequence ; >Assembled_from ") ;
  n3 = keySetMax(ks3) ;
  ks2 = keySetMINUS(ks1, ks3) ;

  keySetDestroy (ks1) ;
  keySetDestroy (ks3) ;

  messout ("Found %d new reads, %d where added, %d did not fit",
	   n1, n3 - n2, n1 -n3 + n2) ;
  if (keySetMax(ks2))
    { Graph g = graphActive () ;
      displayCreate(DtKeySet) ;
      graphRetitle("Exported keyset") ;
      keySetShow (ks2,0) ; ks2 = 0 ;
      keySetSelect () ;
      graphActivate (g) ;
    }
  else
    keySetDestroy (ks2) ;
}

static void fMapTraceAssembleAll    (void)
{ char nm [128], buf [256] ;
  static int taux = 8 ;
  FMAPLOOKGET("fMapTraceAssembleAll") ;

  sprintf (buf, "%d", taux) ;
  messStatus ("Assembling") ;
  if (!messPrompt ("Specify the %% of errors [0, 60]", buf, "iz"))
    return ;
  freeint (&taux) ;
  if (taux < 0 || taux > 60)
    { messout ("Rate out of range [0,60], sorry") ;
      return ;
    }
  fMapTraceCleanAssembly (look->seqKey, 0) ;
  fMapTraceNewName (look->seqKey, nm, 0) ;
  sprintf (buf, 
    "Get 40\nSort 100\nAssemble %d\nFix\nSave %s\n", taux, nm) ;
  fMapTraceGlobalCommand (buf, TRUE) ;
}

static void fMapTraceJoinAll (void)
{ char nm [128], buf [256] ;
  FMAPLOOKGET("fMapTraceJoinAll") ;

  messStatus ("Joining") ;
  fMapTraceCleanAssembly (look->seqKey, 0) ;
  fMapTraceNewName (look->seqKey, nm, 0) ;
  sprintf (buf, 
     "Join\nFix\nSave %s\n", nm) ;
  fMapTraceGlobalCommand (buf, TRUE) ;
}

static void fMapTraceJoinDiagonal (void)
{ FMAPLOOKGET("fMapTraceJoinDiagonal") ;

  messStatus ("Joining") ;
  fMapTraceGlobalCommand ("Join diagonal\n", TRUE) ;
}


void fMapTraceReassembleAll (void)
{ char *nm = 0, *np = 0 ;
  KEYSET ks = 0 ;
  int level ;
  FILE *ff = 0 ;
  KEY last = 0 ;

  messStatus ("Assembling") ;
  if (!messPrompt("To visualize an existing assembly, cancel and double "
		  "click an object\n of the Assembly Class\n\n"
                  "The present command will perform a complete reassembly "
		  "according to the script firstpass.smb\n"
		  "This operation may take several minutes\n\n"
		  "Please cancel or specify what sequences should be assembled\n\n"
		  "-active, assembles the active keyset\n"
		  "-all (the default), assembles all reads","-all","w"))
    return ;
  freenext () ;
  np = strnew (freeprotect(freepos ()), 0) ;
  if (!strncasecmp (np, "-active", 7) && !keySetActive (&ks, 0))
    { messout ("To use -active, please first select a keySet containing sequences") ;
      goto abort ;
    }
  if (!messPrompt("Operations will be echoed in your startup xterm\n\n"
		  "Give a name-prefix for the new assembly","b1","w"))
    goto abort ;
  nm = strnew (freeword (), 0) ;
  ff = filopen ("firstpass","smb", "r") ;
  if (!ff) 
    goto abort ;
  displayPreserve () ;
  level = freesetfile (ff, messprintf (" %s %s ", np, nm)) ;
  defComputeTace (level, ks) ;

  last = (KEY) lastAssembly () ;
  if (last) display (last, 0, 0) ;
 abort:
  messfree (np) ;
  messfree (nm) ;
}

static void fMapTraceFixAll (void)
{ KEYSET ks ;
  int i ;
  KEY key ;
  OBJ obj ;
  FMAPLOOKGET("fMaptraceFixAll") ;
  
  obj = bsCreate (look->seqKey) ;
  if (!obj)
    return ;
  ks = keySetCreate() ; i = 0 ;
  if (bsGetKey (obj, _Subsequence, &key))
    do 
      { keySet (ks, i++) = key ;
       } while (bsGetKey (obj, _bsDown, &key)) ;
  bsDestroy (obj) ;
  
  i = keySetMax (ks) ;
  while (i--)
    dnaAlignFixContig (look->seqKey, keySet (ks, i)) ;
  alignToolsAdjustLink (look->seqKey, ks, 0);
  keySetDestroy (ks) ;
  traceGraphDestroy () ;
  look->pleaseRecompute = TRUE ;
  fMapDraw (look, 0) ;
}

static void fMapTraceAutoEditAll (void)
{ Graph old = graphActive () ;
  messStatus ("Auto Edit") ;

  graphMessage ("This takes a few seconds per read, Type F4 to interupt") ;

  fMapTraceGlobalCommand ("AutoEdit", FALSE) ;
  if (graphActivate (old))
    graphUnMessage () ;
}

static void fMaptraceClipOnTile (void)
{ messStatus ("Clip on Tile") ;

  fMapTraceGlobalCommand ("Clip_on T\n", FALSE) ;
}

static void fMaptraceClipOnExcel (void)
{ messStatus ("Clip on Excellent") ;

  fMapTraceGlobalCommand ("Clip_on E\n", FALSE) ;
}

static void fMaptraceClipOnGood (void)
{ messStatus ("Clip on Good") ;

  fMapTraceGlobalCommand ("Clip_on G\n", FALSE) ;
}

static void fMaptraceClipOnFair (void)
{ messStatus ("Clip on Fair") ;

  fMapTraceGlobalCommand ("Clip_on F\n", FALSE) ;
}

static void fMaptraceClipOnHand (void)
{ messStatus ("Clip on Fair") ;

  fMapTraceGlobalCommand ("Clip_on H\n", FALSE) ;
}

static void fMapTraceCleanUpAssembly (void)
{ 
  FMAPLOOKGET("fMapTraceCleanUpAssembly") ;
  messStatus ("Clean up") ;
  
  fMapTraceCleanAssembly (look->seqKey, 0) ;
  traceGraphDestroy () ;
  look->pleaseRecompute = TRUE ;
  fMapDraw (look, 0) ;
}

static void fMapTraceCopyAssembly (void)
{ char nm[128], buf[256] ;
  KEY key ;
  FMAPLOOKGET("fMapTraceCopyAssembly") ;

  fMapTraceNewName (look->seqKey, nm, 0) ;
  sprintf (buf, "Save_as %s\n", nm) ;
  key = defCptExecuteCommand (look->seqKey, 0, 0, buf, 0) ;
  fMapTraceCleanAssembly (key, 0) ;
  if (key)
    display (key, 0, 0) ;
}

static void fMapTraceOrderContigsBySize (void)
{ char buf[256] ;

  sprintf (buf, "Order_by_Size\n") ;
  fMapTraceGlobalCommand (buf, FALSE) ;
}

static void fMapTraceOrderContigsBySubclones (void)
{ char buf[256] ;

  sprintf (buf, "Order_by_Subclones\n") ;
  fMapTraceGlobalCommand (buf, FALSE) ;
}

static void fMapTraceMeasles (void)
{ FMAPLOOKGET("fMapTraceMeasles") ;
     /* open the assembly window */
  defCptOpen (look->seqKey) ;
}

#ifdef JUNK

static void fMapTraceTrainNeuralNet (void)
{ Graph old = graphActive () ;
  messStatus ("Training") ;

  graphMessage ("This takes a few seconds per read, Type F4 to interupt") ;

  fMapTraceGlobalCommand ("TrainNN", FALSE) ;
  if (graphActivate (old))
    graphUnMessage () ;
}


static void fMapTraceCleanUp (void)
{  
  FMAPLOOKGET("fMapTraceGlobalAction") ;

  messStatus ("Clean up") ;
  traceGraphDestroy () ;

  fMapTraceCleanAssembly (look->seqKey, 0) ;
  look->pleaseRecompute = TRUE ;
  fMapDraw (look, 0) ;
  sessionDoSave (TRUE) ; /* save */
}

static void fMapTraceCompareClip (void)
{ FMAPLOOKGET("fMapTraceCompareClip") ;
  dnaAlignCompareClip (look->seqKey) ;
}

/************/

static void fMapTraceDoAssemblyDotter (KEY key1, KEY key2)
{ Array dna1 = 0, dna2 = 0 ;
  KEY d1, d2 ;
  char *cp1, *cp2 ;

  if (key1 && dnaSubClass (key1, &d1) &&
      key2 && dnaSubClass (key2, &d2))
    { dna1 = dnaGet (d1) ; dna2 = dnaGet (d2) ;
      if (dna1 && dna2)
	{
	  dnaDecodeArray (dna1) ;   dnaDecodeArray (dna2) ;
	  cp1 = messalloc (arrayMax(dna1));
	  cp2 = messalloc (arrayMax(dna2));
	  memcpy(cp1, arrp(dna1,0,char), arrayMax(dna1));
	  memcpy(cp2, arrp(dna2,0,char), arrayMax(dna2));
	  dotter ('N', 0, 
		  name(key1),cp1, 0,
		  name(key2),cp2, 0,
		  0, 0, 0, 0, 0, 
		  0, 0, 0, 0, 0, 0, 0) ;
	}
    }
  arrayDestroy (dna1) ;
  arrayDestroy (dna2) ;
}
#endif

static void fMapTraceDoAssemblyCompare (KEY key)
{ LOOK look = movedLook ;

  graphActivate (look->graph) ;
  contigMoveCancelled () ;
  if (!key || key == look->seqKey)
    return ;
  fMapTraceCleanAssembly (look->seqKey, 0) ;
  if (dnaAlignAssemblyCompare (look->seqKey, key))
    { look->pleaseRecompute = TRUE ;
      fMapDraw (look, 0) ;
    }
  else messout ("Sorry, I cannot compare these objects") ;

}

static void fMapTraceAssemblyCompare (void)
{ FMAPLOOKGET("AssemblyCompare") ;

  displayPreserve () ;
  movedLook = look ;
  displayBlock (fMapTraceDoAssemblyCompare,
		"Pick a contig or an assembly as your target") ;
}

static void fMapTraceDoFuseAssembly (KEY key)
{ char buf[256], nm[128] ;
  KEYSET ks = 0 ;
  KEY link ;
  LOOK look = movedLook ;

  graphActivate (look->graph) ;
  contigMoveCancelled () ;
  if (!key || key == look->seqKey)
    return ;
  if (class (key) != _VSequence && !dnaReClass (key, &key))
    { messout ("Please pick a contig or an assembly as your target") ;
      return ;
    }
  ks = keySetCreate () ;
  keySet (ks, 0) = look->seqKey ;
  keySet (ks, 1) = key ;
  fMapTraceNewName (look->seqKey, nm, 0) ;
  sprintf (buf, "load -active\nget 30\nsort\nassemble 12\nsave %s\n", nm) ;
  traceGraphDestroy () ;
  if ((link = defCptExecuteCommand (0, 0, ks, buf, 0)))
    { look->seqKey = link ;
      look->pleaseRecompute = TRUE ;
      fMapDraw (look, 0) ;
      sessionDoSave (TRUE) ; /* save */
    }
  keySetDestroy (ks) ;
}

static void fMapTraceFuseAssembly (void)
{ FMAPLOOKGET("AssemblyCompare") ;

  displayPreserve () ;
  movedLook = look ;
  displayBlock (fMapTraceDoFuseAssembly,
		"Pick a contig or an assembly as your target") ;
}

/************/

MENUOPT fMapAcemblyOpts[] = {
  {fMapTraceAssembleAll, "Assemble"},
  {fMapTraceJoinAll, "Join"},
  {fMapTraceJoinDiagonal, "Join diagonal"},
  {fMapTraceReassembleAll, "Reassemble from scratch"},
/*   fMapTraceExplodeAll, "Destroy contigs < 1kb", */
  {fMapTraceReInsertLoners, "Re Insert Loners"},
  {fMapTraceReInsertSubclones, "Re Insert Subclones"},
  {fMapTraceAddNewReads, "Add new reads"},
  {fMapTraceFixAll, "Fix"},  
  {fMapTraceAutoEditAll, "AutoEdit"},
  {fMaptraceClipOnTile, "Clip to Tile"},
  {fMaptraceClipOnExcel, "Clip on Excellent"},
  {fMaptraceClipOnGood, "Clip on Good"},
  {fMaptraceClipOnFair, "Clip on Fair"},
  {fMaptraceClipOnHand, "Restore Hand Clips"},
/*   fMapTraceTrainNeuralNet, "Train Neural-Net", */
  {fMapTraceOrderContigsBySize, "Order By Size"},
  {fMapTraceOrderContigsBySubclones, "Order By Subclones"},
  {fMapTraceCleanUpAssembly, "Clean"},
  {fMapTraceCopyAssembly, "Save"},
  {fMapTraceAssemblyCompare, "Compare"},
  {fMapTraceFuseAssembly, "Fuse"},
  {fMapTraceMeasles, "Command window"},
/*  fMapTraceCompareClip, "Clip Comparison", */
   {fMapTraceToggleTrack, "New Error Tracking"},
   {fMapTraceAddWall, "Add Wall at screen center"},
   {0, 0}
} ;

/***************************************************************************************/
/***************************************************************************************/

static void fMapTraceReAssembleContig (LOOK look, KEY contig)
{ KEY link = look->seqKey, subLink, key ;
  int i, j, k, ii ;
  static int taux = 3 ;
  char buf [256], nm[128] ;
  KEYSET reads = 0, newContigs = 0, baddies = 0, ks  = 0 ;
  OBJ obj ;
  Array unit = 0, order = 0 ;
  BSunit *u ;

  sprintf (buf, "%d", taux) ;
  if (!messPrompt ("Specify the %% of errors [0, 60]", buf, "iz"))
    return ;
  freeint (&taux) ;
  if (taux < 0 || taux > 60)
    { messout ("Rate out of range [0,60], sorry") ;
      return ;
    }
  reads = queryKey (contig, ">Assembled_from ; >DNA") ;
  if (keySetMax(reads) < 2)
    { keySetDestroy (reads) ;
      return ;
    }
  displayPreserve () ;
  messStatus ("ReAssembling") ;

     /* extract reads off contig and reassemble it */
  fMapTraceNewName (contig, nm, '_') ;
  sprintf (buf, "Get 20\nSort 100\nAssemble %d\nFix\nSave %s\n", taux, nm) ;
  subLink = defCptExecuteCommand (0, reads, 0, buf, 0) ;  /* destroys reads */
  /* keySetDestroy (reads) ; those are DNA, I want Sequence */
  if (!subLink) return ;
     /* replace it in the look */
     /* in correct order, so queryKey which sorts does not work */
  k = 0 ; 
  unit = arrayCreate (90, BSunit) ;
  if ((obj = bsCreate(link)))
    { if (bsFindTag (obj, _Subsequence))
	bsFlatten (obj, 3, unit) ;
      bsDestroy (obj) ;
    }
  k = 0 ; newContigs = keySetCreate() ;
  if ((obj = bsCreate(subLink)))
    { if (bsGetKey (obj, _Subsequence, &key))
	do 
	  { keySet (newContigs, k++) = key ;
	  } while (bsGetKey (obj, _bsDown, &key)) ;
      bsDestroy (obj) ;
    }
  reads = queryKey (contig, ">Assembled_from") ;
  ks = query (newContigs, ">Assembled_from") ;
  baddies = keySetMINUS (reads, ks) ;
  keySetDestroy (ks) ;
  ks = dnaAlignMakeSubSequence (link, baddies, nm) ;
  order = arrayCreate (90, BSunit) ;
  k = 0 ;
  for (i = 0 ; i < arrayMax (unit) ; i += 3)
    { u = arrp (unit, i, BSunit) ;
      key = u[0].k ;
      if (key != contig)
	{ array (order, k++, BSunit).k = key ;
	  array (order, k++, BSunit).i = u[2].i - u[1].i ;
	  continue ;
	}
      for (j = 0 ; j < keySetMax (newContigs) ; j++)
	{ key = keySet (newContigs, j) ; /* to prevent isolated sequences */
	  if (!keySetFind (baddies, key, &ii))
	    { array (order, k++, BSunit).k = key ;
	      array (order, k++, BSunit).i = 1 ;
	    }
	}
      for (j = 0 ; j < keySetMax (ks) ; j++)
	{ array (order, k++, BSunit).k = keySet (ks, j) ;
	  array (order, k++, BSunit).i = 1 ;
	}
    }
  alignToolsAdjustLink (link, 0, order);

  keySetDestroy (ks) ;
  keySetDestroy (reads) ;
  keySetDestroy (newContigs) ;
  keySetDestroy (baddies) ;
  arrayDestroy (unit) ;
  arrayDestroy (order) ;

  traceGraphDestroy () ;
  look->pleaseRecompute = TRUE ;
  fMapDraw (look, 0) ;
}

/***************************************************************************************/

static void fMapTraceRemoveContig (LOOK look, KEY contig)
{ KEY link = look->seqKey ;
  OBJ Link = 0 ;
  int i, max = 0 ;
  Array unit = 0, order = 0 ;
  BSunit *u ;

  Link = bsCreate (link) ;
  if (!Link)
    return ;
  unit = arrayCreate (90, BSunit) ;
  if (!bsFindKey (Link, _Subsequence, contig) || !bsFindTag (Link, _Subsequence) ||
      !bsFlatten (Link, 3, unit))
    goto abort ;
  order = arrayCreate (90, BSunit) ;
  for (i = 0 ; i < arrayMax (unit) ; i += 3)
    { u = arrp (unit, i, BSunit) ;
      if (u[0].k == contig)
	continue ;
      array (order, max++, BSunit).k = u[0].k ;
      array (order, max++, BSunit).i = u[2].i - u[1].i ;
    }
  alignToolsAdjustLink (link, 0, order) ;
  arrayDestroy (order) ;
 abort:
  bsDestroy (Link) ;
  arrayDestroy (unit) ; 
  traceGraphDestroy () ;
  look->pleaseRecompute = TRUE ;
  fMapDraw (look, 0) ;
}

/***************************************************************************************/
/***************************************************************************************/

void abiFixFinish (void)
{
  KEY link ;
  FMAPLOOKGET("abiFixFinish") ;

  messStatus ("Finish") ;
  link = look->seqKey ;
  abiFixDoFinish (link) ;
}

/***************************************************************************************/
/***************************************************************************************/

void fMapGelDisplay (void)
{
  extern void gelComparativeDisplay (KEY clone, KEY link) ; /* from geldisp.c */
  KEYSET ks = query (0, "FIND Clone") ;
  FMAPLOOKGET("abiGelDisplay") ;

  if (keySetMax (ks) )
    gelComparativeDisplay (keySet(ks, 0), look->seqKey) ;
  keySetDestroy (ks) ;  
}

/***************************************************************************************/
/***************************************************************************************/

void fMapReDrawWindow (void)
{
  FMAPLOOKGET("Redraw") ;
  
  look->pleaseRecompute = TRUE ;
  fMapDraw (look, 0) ;
}

/***************************************************************************************/
/***************************************************************************************/

void fMapTraceDestroy (LOOK look)
{ Array a ;
  Associator 
    all1 = look ? look->taggedBases : 0 ,
    all2 = look ? look->virtualErrors : 0 ;
  void *vp = 0 ;
  
  defCptDestroyLook (look->seqKey) ;
  if (assExists(all1))
    { while (assNext (all1, &vp, &a))
	{ if (arrayExists(a))
	    arrayDestroy (a) ;
	}
      assDestroy (look->taggedBases) ;
    }
  vp = 0 ;
  if (assExists(all2))
    { while (assNext (all2, &vp, &a))
	{ if (arrayExists(a))
	    arrayDestroy (a) ;
	}
      assDestroy (look->virtualErrors) ;
    }
  if (arrayExists (look->contigCol))
    arrayDestroy (look->contigCol) ;
  if (graphAssFind (&multipletMagic, &a))
    { if (arrayExists(a))
	arrayDestroy (a) ;
      graphAssRemove (&multipletMagic) ;
    }
}

/****************************************/

