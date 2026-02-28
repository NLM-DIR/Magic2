/*  File: ksetdisp.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Manipulate Key Sets.
      keyset windows have their own menus and  KEYSET. When they 
      are destroyed they (VERY IMPORTANT) destroy the keyset.
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 16 16:13 1999 (fw)
 *			-	fn() to void(void) to please WIN32 compiler
 * * Jun 10 16:53 1996 (il)
 * * Oct 15 15:15 1992 (mieg): Mailer
 * * Nov  5 18:27 1991 (mieg): I introduced keySetAlpha
 * * Oct 23 18:53 1991 (mieg): changed name dump to a keySetDump
 * * May 1996 (il) : added scroll bar to keyset's and rearranged the menus.
 * *               : ability to select/highlight members of the keyset and 
 * *               : perform the usual operation's on them. 
 * Created: Wed Oct 23 18:52:00 1991 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: ksetdisp.c,v 1.25 2020/05/30 16:50:32 mieg Exp $ */

#include "acedb.h"

#include "display.h"
#include "lex.h"
#include "dump.h"
#include "tags.h"		
#include "systags.h"
#include "dna.h"
#include "a.h"
#include "bs.h"
#include "sysclass.h"
#include "session.h"     
#include "query.h"
#include "pick.h"
#include "biblio.h"
#include "peptide.h"
#include "tree.h"
#include "forest.h"
#include "main.h"
#include "help.h"


	/* keySetShow () displays a KEYSET (a sorted array of KEY's) in 
	   multicolumn format, allowing the user to pick and move around 
	   with arrow keys.  Objects are displayed with default displaytype.

	   To initialise the system, call with handle==0.  It returns the new 
	   handle.  To reuse the graph call first with (0,handle) to
	   kill the alphabetic keyset. The aim of this is to
	   allow the user to change the displayed KEYSET.

	   keySetCurrent() allows the user to set the key directly (a bit
	   like pick).  It takes the key as argument, not the box number.
	   It returns the display type which the keyset window has been set to
           use (zero for default) suitable for use as the third argument of 
            display()

	   Although we can save to permanent objects (by the menu Save option)
	   and recover from them (via display() and keySetDispCreate()),
	   a keyset display window is NEVER attached to an object - it always
	   deals in copies of the permanent keyset.  This allows much greater
	   flexibility in changing the contents of a window, avoids tangles
	   with the cache, and is absolutely necessary.
	*/
int alignDumpKeySet (KEYSET kSet, FILE *fil);
BOOL parseBuffer (char *text, KEYSET ks);
static void queByExam(void);
void qbeCreateFromKeySet (KEYSET ksAlpha, int i0);
void acedbMailer (KEY key, KEYSET ks, Stack sText);
BOOL doMultiMapDisplay (KEYSET ks);

BOOL (*mhmpDisplayp) (KEYSET ks, int glyph, KEY colour) = NULL ;

#define COLWIDTH 16
#define MAXOBS 1000
#define SCROLLBARGAP 4
#define SCROLLBARPOS 1.0
#define MINSCROLLBARLENTH 0.7
#define BOTGAP 1
#define DEFAULT 0

/* check that these numbers are ok in menuOptions */
#define OR 21
#define AND 22
#define XOR 23
#define MINUS 24
#define UNDO 25

#define UP TRUE
#define DOWN FALSE
#define MAX_CLASSES 200

static BOOL isOtherTools = TRUE ; /* extended set of buttons for specialists */
static float TOPGAP ;
static int MAGIC = 187641 ;	/* also use address as graphAss handle */

typedef struct LOOKstruct
{ 
  int   magic ;
  KEYSET kSet , keySetAlpha ,selected, undo;
  BOOL fullAlpha ; /* complete call to keySetAlphaHeap done */
  KEY   curr ;
  Graph graph, owner ;
  BOOL  isBlocking ;
  int   base, top ;		/* top is 1 + index of last displayed key */
  int numperpage, showtype,colWidth[50];
  Array box2key ;
  BOOL showClass ;
  int showQuery; /* 0: no, 1:yes-nofocus, 2:yes-focus */
  float scrolllen;
  char txtBuffer[60], queryBuffer[256] ;
  int txtbox,numselpage,queryBox;
  int scrollBarBox,upBox,downBox,firstBox,lastBox;
  int itemBoxStart,itemBoxEnd;
  int nitem, nselected ; char touched ;
  char *message ; 
  KEY keySetKey ; /* the key if displaying a KEYSET class object */
} *LOOK ;

static int  N3rdLINE, NQueryLine ;
static FREEOPT menuOptions[] = { 
  {25,"Key Set ?"},
  {99,"Quit"},
  {98,"Help"},
  {97, "Print"},
  {0,""},
  {7,"New empty keyset"},
  {10,"Import keyset"},
  {2,"Copy whole keyset"},
  {1,"Save as"},  
  {0,""},
  {6,"Add-Remove toggle"},
  {OR,"OR"},  /* reserver 20-26 for these */
  {AND,"AND"},
  {MINUS,"MINUS"},
  {XOR,"XOR"},
  {UNDO,"UNDO"},
  {0,""},
/* 47,"Show/Hide Classes", */
  {41,"Mail"},  
  {51,"Keyset Summary"},
  {45,"Show as Multimap"},
  {43, "Show as Chrom map"},
  {44, "Show as h map"},
  {49, "Show as locator"},
  {60, "Show as expression"},
  {101, "Show as wiggle"},
  {8,"Show as Text"},
} ;

static void localDestroy (void) ;
static void localPick (int box,double x,double y) ;
static void localKbd (int k) ;
static void ksetlocalMenu (KEY k) ;
static void pageUp(void), pageDown(void) ;
static void* keySetShow2 (KEYSET kSet, LOOK look) ;
static void ksetAddKey (KEY key) ;
static void ksetDeleteKey (KEY key) ;
static void ksetUnBlock (void) ;
static void controlLeftDown (double x, double y);
static void controlLeftUp (double x, double y);
static void controlLeftDrag (double x, double y);
static void selectedToNewKeySet(void);
static BOOL columnWidth(int startpos,int direction,int *width,int *endpos,int numrows);
static void drawMenus(void);
static void drawScrollBar(int max);
static void keySetMix(LOOK qry, int combinetype);
static void keySetSummary(void);
static void changeShowtype(KEY key);
static void nextShowtype(void);
static void scrollMiddleUp (double x, double y);

#define LOOKGET(name)     LOOK look ; \
       			  if (!graphAssFind (&MAGIC, &look)) \
		            messcrash ("graph not found in %s", name) ; \
                          if (look->magic != MAGIC) \
                            messcrash ("%s received a wrong pointer",name)

static LOOK selectedKeySet = 0 ;

static int startbox,lastbox,numboxes,oldcol, oldy ;

static int scrollTriangle(BOOL direction, float x, float y,int colour)
{ int i = 3, box,col;
  float incr;
  Array temp ;
  LOOKGET ("scrollTriangle") ;

  if(direction)
    incr = -1.0;
  else
    incr = 1.0;

  box = graphBoxStart();
  col = graphColor(colour);

  temp = arrayCreate(2*i, float) ;
  array(temp, 0, float) = x;
  array(temp, 1, float) = y;
  array(temp, 2, float) = x+2.0;
  array(temp, 3, float) = y;
  array(temp, 4, float) = x+1.0;
  array(temp, 5, float) = y+incr;
  graphPolygon(temp);
  col = graphColor(col);

  graphLine(x, y, x+2.0, y);
  graphLine(x+2.0, y, x+1.0, y+incr);
  graphLine(x+1.0, y+incr, x, y);
  arrayDestroy(temp);
  
  graphBoxEnd();
  
  return box;
}

/******************************/

static int nItem (LOOK look)
{ 
  int i, nitem = 0 ;

  if (!(look->touched & 2)) return look->nitem ;
  if (look->fullAlpha && keySetExists(look->keySetAlpha))
    return keySetMax(look->keySetAlpha) ;
  
  if (keySetExists(look->kSet))
    { KEY *kp = arrp (look->kSet, 0, KEY) - 1 ;
    
      i = keySetMax(look->kSet) ;
      while (kp++, i--)
	if (lexIsKeyVisible(*kp)) nitem++ ;
    }
  look->nitem = nitem ;
  look->touched &= 0x1 ;  /* not 2  */
  return nitem ;
}
/* mhmp 04.09.98 forcement necessaire comme nItem */
static int nSelected (LOOK look)
{ 
  unsigned int i ; 
  int nselected = 0 ;

  if (!(look->touched & 1)) return look->nselected ;
  if (keySetExists(look->selected) && keySetMax (look->selected))
    { KEY *kp = arrp (look->selected, 0, KEY) - 1 ;
    
      i = keySetMax(look->selected) ;
      while (kp++, i--)
	if (*kp && lexIsKeyVisible(*kp)) nselected++ ;
    }
  look->nselected = nselected ;
  look->touched &= 0x2 ;  /* not 1  */
  return nselected ;
}

static void shownItem (LOOK look)
{

  if (look->touched)
    {
      nItem (look) ;
      nSelected (look) ;
    }
  if(look->numselpage == look->nselected)
    sprintf(look->txtBuffer,"%d items %d selected", look->nitem, look->nselected);
  else
    sprintf(look->txtBuffer,"%d items %d selected (%d off screen)", 
    look->nitem, look->nselected,look->nselected-look->numselpage) ;
}


/******************************/
/** fake biblio functions   **/
/** for model free kernel   **/
/******************************/
static BOOL setPage(float *z)
{
  int w,h ;
  float temp,top,bottom;
  int max;

  LOOKGET ("setPage") ;

  if (look != selectedKeySet)
    keySetSelect () ; 

  look->curr = 0;
  max = nItem (look) ;
  graphFitBounds(&w,&h);

  top  =  (float)TOPGAP;
  bottom = (float)(h - BOTGAP);

  if(*z<=top)
    look->base = 0;
  else if(*z+look->scrolllen >= bottom){    
    look->base = max;
    if (keySetMax(look->keySetAlpha) < keySetMax(look->kSet) &&    /* load extra items needed */
	keySetMax(look->keySetAlpha) < look->base)
      { keySetDestroy(look->keySetAlpha) ;
	look->keySetAlpha 
	  = keySetAlphaHeap(look->kSet, keySetMax(look->kSet)) ;
	look->fullAlpha = TRUE ;
      }
  }
  else{
    if((int)look->scrolllen==MINSCROLLBARLENTH)
      temp = (*z-top)/(bottom-top-look->scrolllen);
    else
      temp = (*z-top)/(bottom-top);
    look->base = (int)(temp * max);
    if(look->base > max){
      look->base=max;
    }
  }
  return TRUE;
}

static void keyBandDrag(float *x, float *y, BOOL isUP)
{
  float top,bottom;
  int nx,ny;
  LOOKGET ("keyBandDrag") ;
  *x =  0.5; /* x is set to the middle */

  top  =  (float)TOPGAP;
  graphFitBounds (&nx, &ny) ;
  bottom = (float)(ny - BOTGAP);

  if(*y > bottom - look->scrolllen)
    *y = bottom - look->scrolllen;
  if(*y < top )
    *y = top;
  if(isUP) /* after scrolling has finished */
    {
      if(setPage(y)){
	if(look->base == keySetMax(look->keySetAlpha)){
	  --look->base;
	  pageUp();
	  return;
	}
	if (look != selectedKeySet)
	  keySetSelect () ; 
	keySetShow2 (look->kSet, look) ;
      }
    }
}

static void pageUp(void)
{
  int nx,ny,colX,startpos,endpos,numrows,cWidth;
  BOOL room= TRUE;
  LOOKGET ("pageUp") ;

  if (look->base)
    {
      graphFitBounds (&nx, &ny) ;
      colX = SCROLLBARGAP;
      startpos = look->base+1; /* il +1 to keep last item still on page */
      numrows = ny-(TOPGAP+BOTGAP);
      while(room){    
	startpos--;
	room = columnWidth(startpos,FALSE,&cWidth,&endpos,numrows);
	/* mhmp comme dans keySetShow2 */
	if((cWidth && (colX + cWidth) < nx) || 
	   colX == SCROLLBARGAP){ /* if it fits or is the first */ 
	  startpos= endpos;
	  look->base = endpos;
	}
	else
	  room = FALSE;
	colX += cWidth +2;
      }
      if(look->base < 0) /* returns -1 if before start of first page */
	look->base = 0;
      look->curr = 0;
      keySetShow2 (look->kSet, look) ;
    }
}

static void pageDown(void)
{
  LOOKGET ("pageDown") ;
  if (look->top+1 < keySetMax(look->keySetAlpha))
    { 
	  look->base = look->top ; /* removed +1 to keep last item on next page */
	  look->curr = 0;
	  keySetShow2 (look->kSet, look) ;
    }
}

static void localResize (void)
{
  LOOKGET("localResize") ;

  look->curr = 0;
  keySetShow2 (look->kSet, look) ;
}

static void none(void)
{
}

static void acedump(void)
{
  /* .ace dump */
  FILE *fil;
  static char fileName[FIL_BUFFER_SIZE] , dirName[DIR_BUFFER_SIZE] ;
  register int i;

  LOOKGET("acedump");
  if (!look->keySetAlpha)
    return ;
  if (!(fil = filqueryopen (dirName, fileName, "ace", "w",
			    "Where do you wish to export the data ?")))
    return ;
  fprintf (fil,"// data dumped from keyset display\n\n") ;
  if (keySetMax(look->kSet) > keySetMax(look->keySetAlpha))
    { keySetDestroy(look->keySetAlpha) ;
      look->keySetAlpha = 
	    keySetAlphaHeap(look->kSet, keySetMax(look->kSet)) ;
      look->fullAlpha = TRUE ;
    }
  for (i = 0 ; i < keySetMax(look->keySetAlpha) ; ++i)
    if (!dumpKey (keySet(look->keySetAlpha,i), fil, 0))
      if (!messQuery("Do you want to proceed ?"))
	break ;
  filclose (fil) ;

}
static void acedumptimestamp(void)
{
  /* .ace dump */
  FILE *fil;
  static char fileName[FIL_BUFFER_SIZE] , dirName[DIR_BUFFER_SIZE] ;
  register int i;
  LOOKGET("acedump+(timestamp)");

  if (!look->keySetAlpha)
	return ;
  if (!(fil = filqueryopen (dirName, fileName, "ace", "w",
			    "Where do you wish to export the data ?")))
    return ;
  fprintf (fil,"// data dumped from keyset display\n\n") ;
  if (keySetMax(look->kSet) > keySetMax(look->keySetAlpha))
    { keySetDestroy(look->keySetAlpha) ;
      look->keySetAlpha = 
	keySetAlphaHeap(look->kSet, keySetMax(look->kSet)) ;
      look->fullAlpha = TRUE ;
    }
  dumpTimeStamps = TRUE ;
  for (i = 0 ; i < keySetMax(look->keySetAlpha) ; ++i)
    if (!dumpKey (keySet(look->keySetAlpha,i), fil, 0))
      if (!messQuery("Do you want to proceed ?"))
	break ;
  filclose (fil) ;
  dumpTimeStamps = FALSE ;
  return;
}

static void nameDump(void)
{	/* name dump */
  FILE *fil;
  static char fileName[FIL_BUFFER_SIZE] , dirName[DIR_BUFFER_SIZE] ;
  LOOKGET("namedump");
  if (!look->keySetAlpha)
    return ;
  if (!(fil = filqueryopen (dirName, fileName, "list", "w",
			    "Where do you wish to export the names ?")))
    return ;
  fprintf(fil,"KeySet : \"%s\"\n",fileName) ;
  if (keySetMax(look->kSet) > keySetMax(look->keySetAlpha))
    { keySetDestroy(look->keySetAlpha) ;
      look->keySetAlpha = 
	keySetAlphaHeap(look->kSet, keySetMax(look->kSet)) ;
      look->fullAlpha = TRUE ;
    }
  keySetDump(fil,0 ,look->keySetAlpha) ;
  filclose (fil) ;
  return ;
}
static void fastaDump(void)
{
     /* FastA format sequence dump */
  FILE *fil;
  LOOKGET("fastaDump");

  if (!look->keySetAlpha)
    return ;
  if (!(fil = dnaFileOpen()))
    return ;
  if (keySetMax(look->kSet) > keySetMax(look->keySetAlpha))
    { keySetDestroy(look->keySetAlpha) ;
      look->keySetAlpha = 
	keySetAlphaHeap(look->kSet, keySetMax(look->kSet)) ;
      look->fullAlpha = TRUE ;
    }
  dnaDumpFastAKeySet (look->keySetAlpha, fil, 0) ;
  filclose (fil) ;
  return ;
}
static void protDump(void)
{
  /* FastA format sequence dump */
  FILE *fil;
  LOOKGET("Prot FASTA Dump");
  
  if (!look->keySetAlpha)
    return ;
  if (!(fil = pepFileOpen()))
    return ;
  if (keySetMax(look->kSet) > keySetMax(look->keySetAlpha))
    { keySetDestroy(look->keySetAlpha) ;
      look->keySetAlpha = 
	keySetAlphaHeap(look->kSet, keySetMax(look->kSet)) ;
      look->fullAlpha = TRUE ;
      }
  pepDumpFastAKeySet (look->keySetAlpha, fil, 0) ;
  filclose (fil) ;
  return ;
}

static void alignDump(void)
{
   /* Alignment dump */
  FILE *fil;
  LOOKGET("alignDump");

  if (!(fil = pepFileOpen()))
    return ;
  alignDumpKeySet (look->kSet, fil) ;
  filclose (fil) ;
  return ;
}

static void selectAllKeySet(void)
{
  LOOKGET("SelectAllKeySet");
  
  keySetDestroy(look->selected);
  look->touched |= 1 ;
  look->selected = keySetCopy(look->kSet);  
  keySetShow2(look->kSet,look) ;
}

static void reverseSelectedKeySet(void)
{
  KEYSET temp;
  LOOKGET("SelectAllKeySet");

  look->curr = 0 ;
  look->touched |= 1 ;
  temp = keySetCopy(look->selected);
  look->selected = keySetXOR(temp,look->kSet);
  keySetDestroy(temp);  
  keySetShow2(look->kSet,look) ;

}

static void removeSelectedKeySet(void)
{
  LOOKGET("removeSelectedKeySet");

  keySetDestroy (look->undo) ;
  look->undo = look->kSet ;

  look->kSet = keySetXOR(look->undo,look->selected);

  look->curr = 0;
  look->touched |= 3 ;
  keySetDestroy(look->keySetAlpha) ;
  keySetReCreate(look->selected);

  keySetShow2(look->kSet,look) ;
}

static void clearSelected(void)
{
  LOOKGET("clearSelected");

  look->curr = 0 ;
  look->touched |= 1 ;
  keySetReCreate(look->selected);
  keySetShow2(look->kSet,look) ;
  
}

static void editSelected(void)
{
  int i, max, max1;
  KEY *kp;
  char *cp ;
  KEYSET ksNew = 0 ;
  static char* buffer ;
  Stack localStack = 0 ;
  LOOKGET("editSelected");

  if (!look->selected || !(max = keySetMax(look->selected)))
    { messout("first select the objects you want to edit") ;
      return ;
    }	
  i = max ;
  kp = arrp (look->selected, 0, KEY) - 1 ; max1 = 0 ;
  while (kp++, i--)
    if (pickType(*kp) == 'B')  /* don t edit array this way */
      max1++ ;
  if (max != max1)
    { messout("You have selected %d  non-text objects, \n"
	      "they cannot be edited in this way", max - max1) ;
    return ;
    }
  if(!checkWriteAccess())
    return;
  cp = messprintf("%s%d%s",
		 "Type a .ace command, it will be applied to the ",
		  max, " currently selected objects") ;
  if (!messPrompt(cp, buffer ? buffer : "", "t"))
    return ;
  cp = freeword() ; if (!cp) return ;
  messfree(buffer) ;
  buffer= strnew(cp, 0) ;

  i = max ;
  localStack = stackCreate (40*i) ;

  kp = arrp (look->selected, 0, KEY) - 1 ;
  while (kp++, i--)
    if (pickType(*kp) == 'B')  /* don t edit array this way */
      { catText (localStack, className (*kp)) ;
      catText (localStack, " ") ;
      catText (localStack, freeprotect(name (*kp))) ;
      catText (localStack, "\n") ;
      catText (localStack, buffer) ;
      catText (localStack, "\n\n") ;
      }

  ksNew = keySetCreate () ;
  parseBuffer (stackText (localStack,0), ksNew) ;
  stackDestroy (localStack) ;  /* it can be quite big */
  i = keySetMax(ksNew) ;
  for (cp = buffer ; *cp ; cp++)
    if (*cp == '%') *cp = '.' ;
  messout (messprintf ("// I updated %d objects with command\n %s\n", i, buffer)) ;
  messdump (messprintf ("// I updated %d objects with command\n %s\n", i, buffer)) ;
  /*  parseKeepGoing = FALSE ; */
  keySetDestroy (ksNew) ;
}
 
static void killSelected(void)
{
  int max;
  register int i;
  KEY key;

  LOOKGET("killSelected");

  if (!look->selected || !(max = keySetMax(look->selected)))
    { messout("first select the objects you want to kill") ;
      return ;
    }
  if(!checkWriteAccess())
    return;
  if (!messQuery (messprintf ("Do you really want to destroy %d items", max)))
    return ;
  look->touched |= 3 ;
  i= 0;
  while(i<max){
    key = keySet(look->selected,i) ;
    switch ( pickType(key))
      {
      case 'A':
	arrayKill (key) ;
	break ;
      case 'B':
	{ OBJ obj ;
	  if ((obj = bsUpdate (key)))
	    bsKill (obj) ;	      
	}
	break ;
      }
    { KEY *kp = arrp(look->kSet, 0, KEY) ;
      int max = keySetMax (look->kSet) ;
      
      while (max-- && *kp != key) kp++ ;
      while (kp++, max--)
	*(kp - 1) = *kp ;
    }
    keySetMax (look->kSet)-- ;
    i++;
  }
  keySetDestroy(look->keySetAlpha) ;
  keySetReCreate(look->selected);
  look->curr = 0; /* else crashes if remove last item */
  keySetShow2(look->kSet,look) ;
}

#ifdef JUNK
static void editAll (void)
{
  selectAllKeySet () ;
  editSelected () ;
}

static void killAll (void)
{
  selectAllKeySet () ;
  killSelected () ;
}
#endif

static void beforecombine( int type)
{
  LOOKGET("before Combine");

  if (selectedKeySet && selectedKeySet->magic == MAGIC)
    { displayPreserve() ; look->owner = -1 ;
      clearSelected() ;
      look->touched  |= 3 ;
      keySetMix(look, type) ;
    }
  else
    messout
      ("%s%s",
       "First select a key Set window by picking ",
       "its background with the left mouse button") ;
}


static void combineXOR(void)
{
  beforecombine(XOR);
}

static void combineAND(void)
{
  beforecombine(AND);
}

static void combineOR(void)
{
  beforecombine(OR);
}

static void combineMINUS(void)
{
  beforecombine(MINUS);
}

static void follow(void)
{
  int i;
  KEY key = 0;

  LOOKGET("follow");

  if (!look->keySetAlpha)
    return ;
  i = look->curr-2+look->base ;
  if (i < 0 || i > keySetMax(look->keySetAlpha))
    i = 0 ;
  if (!keySetMax(look->keySetAlpha))
    return ;
  key = keySet(look->keySetAlpha,i) ;
  if (pickType(key) != 'B')
    return ;
  { int type, tcl ;
    KEY tag ;
    Stack sta = stackCreate (50) ;
    if (treeChooseTagFromModel (&type, &tcl, class(key), &tag, sta, 0))
      { KEYSET tmp = query(look->kSet, messprintf(">%s",stackText (sta, 0))) ;

	keySetDestroy(look->keySetAlpha) ;
	keySetDestroy (look->undo) ;
	look->undo = look->kSet ;

	look->kSet = tmp ;
	look->touched |= 3 ;
	look->curr = 0 ;
	keySetShow2(look->kSet,look) ;
      }
    stackDestroy (sta) ;
  }
}

static void grepset (void)
{
  char *cp ;
  int i ;
  KEYSET ks = 0 ;
  Stack s = 0 ;
  LOOKGET("grepset");

  if (!messPrompt("Type a few words (at least 3 letters long), upper/lower case do not matter","", "t"))
    return ;
  s = stackCreate(50) ;
  catText (s, "*") ;
  while ((cp = freeword()))
    { catText (s, cp) ; catText (s, "*") ; }
  cp = stackText (s, 0) ; i = 0 ;
  while (*cp) { if (*cp != ' ' && *cp != '*') i++ ; cp++ ;}
  if (i >= 3) 
    {
      ks = queryGrep (look->kSet, stackText (s, 0)) ;
      keySetNewDisplay (ks, stackText(s,0)) ; /* will destroy ks */
      keySetSelect () ;
    }
  else
    messout ("word too short") ;
  stackDestroy (s) ;
}

static void relatedBiblio(void)
{
  LOOKGET("relatedBiblio");
  
  biblioKeySet ("Biblio for keyset",look->kSet) ;
}
 
void *keySetNewDisplay (KEYSET ks, const char *title)
{
  LOOK look = 0 ;
  if (keySetExists(ks))
    { 
      Graph g = graphActive () ;
      displayCreate(DtKeySet) ;
      graphRetitle(title) ;
      look = keySetShow (ks,0) ;
      keySetSelect () ;
      displayPreserve() ;
      graphActivate (g) ;
    }
  return look ;
}

static void keySetDoQuery(char *txt)
{
  KEYSET ks = 0 ;
  char *cp ;
  LOOKGET("keySetDoQuery");

  cp = look->queryBuffer ;
  while (*cp == ' ') cp++ ;
  if (!*cp) return ;
  ks = query (look->kSet, cp) ;
  keySetNewDisplay (ks, look->queryBuffer) ; /* will destroy ks */
  keySetSelect () ;
}

static void openQueryBox(void)
{
  LOOKGET("openQueryBox");
  if(look->curr)
    look->curr -= (look->showQuery ? 1 : -1 ) * NQueryLine ; 
  look->showQuery = look->showQuery ? 0 : 2 ; /* toggle */  
  keySetShow2 (look->kSet, look) ;
}

static void findneigh(void)
{ 
  LOOKGET("findneigh");
  if (!look->keySetAlpha)
    return ;
  keySetDestroy (look->undo) ;
  look->undo = look->kSet ;

  look->kSet = keySetNeighbours(look->undo) ; /* will destroy tmp */

  look->touched |= 3 ;
  look->curr = 0 ;
  keySetDestroy(look->keySetAlpha) ;
  keySetShow2(look->kSet,look) ;
}

static void queByExam(void)
{
  int i;
  LOOKGET("queByExam");

  if (!look->keySetAlpha)
    return ;
  i = look->curr-2+look->base ;
  if (i < 0 || i >= keySetMax (look->keySetAlpha))
    i = 0 ;
  qbeCreateFromKeySet (look->keySetAlpha,i - 1) ;
}

static void showClasses(void)
{
  LOOKGET("showClasses");

  look->curr = 0;
  look->showClass = !look->showClass;
  keySetShow2(look->kSet,look) ;
}

static void* keySetShow2 (KEYSET kSet, LOOK look)
{  
  int	i, y, colX, col, nx, ny, numcol;
  int   max = ( kSet ? keySetMax(kSet) : 0 ) ;
  const char	*cp ;
  KEY   key ;
  int   box=0;
  BOOL FIRST = TRUE, room = TRUE ;
  int startpos=0, cWidth=0, endpos=0, numrows=0 ;

  if (!look->isBlocking)	/* so as not to hide message */
    graphPop() ;

   if (look->base >= max)
    look->base = 0 ;

  { KEYSET dummyks ;
    void *dummylook ;
    if (!keySetActive (&dummyks, &dummylook))
      selectedKeySet = look ;
  }

  graphFitBounds (&nx, &ny) ;
  if (look->kSet != kSet)
    {
      look->showQuery = 0 ;
      keySetDestroy (look->undo) ;
      look->undo = keySetCopy (look->kSet) ;
      if (keySetExists(look->keySetAlpha))
	keySetDestroy(look->keySetAlpha) ;
      /* don't destroy look->kSet here, calling routine should */
    }
  look->kSet = kSet ;
  look->touched |= 3 ;
  look->curr = 0 ;

  drawMenus();  /* here, hence after evaluation of biblio possible */
  if (max)
    { col = 1 + nx/2 ;	/* max number possible columns */
      if (!look->keySetAlpha)
        { look->keySetAlpha = keySetAlphaHeap(kSet, look->base + 4*col*ny) ;
	  look->fullAlpha = FALSE ;
        }
      else
	if (keySetMax(look->keySetAlpha) < keySetMax(look->kSet) &&
	    keySetMax(look->keySetAlpha) < look->base + col*ny)
	  { keySetDestroy(look->keySetAlpha) ;
	    look->keySetAlpha 
	      = keySetAlphaHeap(look->kSet, keySetMax(look->kSet)) ;
	    look->fullAlpha = TRUE ;
	  }
    }

  look->box2key = arrayReCreate(look->box2key, 100, KEY) ;
  y = TOPGAP;
  if (look->message)
    { graphText (look->message, 2, TOPGAP -1.98) ; TOPGAP += 1.2 ; }
  if (max){
  look->numselpage = 0;
  colX = SCROLLBARGAP;
  endpos = look->base -1;
  numrows = (ny-TOPGAP)-BOTGAP;
  if (numrows < 0) 
    return look; /* window too small to start with */
  room = TRUE;
  i = numcol = 0;
  while (room)
    {    
      startpos = endpos+1;
      room = columnWidth(startpos,TRUE,&cWidth,&endpos,numrows);
      look->colWidth[i++] = cWidth;
      if (cWidth && ((colX + cWidth) < nx || colX == SCROLLBARGAP))
	{ /* if it fits or is the first */ 
	  look->top = endpos;
	  y = TOPGAP;
	  numcol++;
	  while (startpos <= endpos)
	    {
	      key = arr(look->keySetAlpha,startpos++,KEY) ;
	      cp = 0 ;
	      if (!nextName (key, &cp) || *cp == '\177')
		continue ;

	      box = graphBoxStart() ;
	      array (look->box2key, box, KEY) = key ;
	      if (iskey(key) == 2)
		graphTextFormat (BOLD) ;  /* pickable object */
	      else
		graphTextFormat (PLAIN_FORMAT) ;
	      if (look->showClass)
		graphText (messprintf ("%s:%s", className(key), cp), 
			   colX, y++) ;
	      else
		graphText (cp, colX, y++) ;
	      graphBoxEnd () ;
	      graphBoxFreeMenu (box, (FreeMenuFunction) ksetlocalMenu, menuOptions) ;
	      if (FIRST)
		{ look->firstBox = box;
		  FIRST = FALSE;
		}
	      if (keySetFind (look->selected, key, 0))
		{ graphBoxDraw(box,BLACK,LIGHTGRAY);
		  look->numselpage++;
		}
	    }
	  colX += (cWidth+2);
	}
      else if (cWidth != 0) /* Do not stop on blank keys. For some reason
			       some databases have multiple ""'s */
	room = FALSE;
    }
  look->lastBox = box;

  /*  look->numperpage = numcol*numrows; mhmp */
  look->numperpage = look->lastBox - look->firstBox + 1 ;
  array(look->box2key, arrayMax(look->box2key), KEY) = 0 ; /* terminate */

  graphTextFormat (PLAIN_FORMAT) ;

  if(look->curr)
    {
      if(look->curr >= look->firstBox && look->curr <= look->lastBox) /* current may not be on the page now */
	graphBoxDraw(look->curr,WHITE,BLACK);                         /* due to other tools being expanded  */
      else
	look->curr = 0;
    }

  shownItem (look) ;

  look->txtbox = graphBoxStart();
  graphTextPtr(look->txtBuffer, (nx/2)-15, TOPGAP - 1.98, 60) ; 
  graphBoxEnd();
  /* mhmp 11.09.98 */
  drawScrollBar(look->nitem);
}
  else{
    sprintf(look->txtBuffer,"0 items 0 selected"); 
    look->txtbox = graphBoxStart();
    graphTextPtr(look->txtBuffer, (nx/2)-15,  TOPGAP - 1.98,60) ;
    graphBoxEnd();
  }
  
  graphRedraw();
  messfree (look->message) ;
  return look ;
}

void* keySetShow (KEYSET kSet, void* handle)
{
  return graphActive() ? keySetMessageShow (kSet, handle, 0) : 0 ;
}

void* keySetMessageShow (KEYSET kSet, void* handle, char *message)
{ Graph old = graphActive() ;
  LOOK  look = handle ;

  if (look && look->magic != MAGIC)
    messcrash("keySetShow called with corrupted handle") ;
  if (look && !graphActivate (look->graph))
    messout ("keySetShow lost its graph - taking over this one") ;

  if (!look)
    { look = (LOOK) messalloc (sizeof (struct LOOKstruct)) ;
      look->magic = MAGIC ;
      look->graph = graphActive() ;
      look->owner = old ;
      look->showtype = DEFAULT;
      look->selected = keySetCreate();
      look->showClass = FALSE;
      look->showQuery = 0 ;
      look->curr = 0;
      look->numselpage = 0;
      graphRegister (DESTROY, localDestroy) ;
      graphRegister (MESSAGE_DESTROY, ksetUnBlock) ;
      graphRegister (PICK, localPick) ;
      graphRegister (KEYBOARD, localKbd) ;
      graphRegister (RESIZE, localResize) ;
      graphRegister (MIDDLE_UP, scrollMiddleUp) ;
      graphHelp ("KeySet") ;
      graphFreeMenu ((FreeMenuFunction) ksetlocalMenu, menuOptions) ;

      graphAssociate (&MAGIC, look) ;
      if (!selectedKeySet)
	selectedKeySet = look ;
    }
  else
    { 
      if (look->magic != MAGIC)
	messcrash ("keySetShow called with bad handle") ;
      look->curr = 0 ;
      look->base = 0 ;
      messfree (look->message) ;
      keySetDestroy(look->keySetAlpha) ;
      keySetReCreate(look->selected);
    }
  if (message) 
    look->message = strnew (message, 0) ;
  look->touched |= 3 ;
  look->curr = 0 ;
  return keySetShow2(kSet, look) ;
}

/*****************************/

/* create a new keyset list display */
void newKeySet (char *title)
{
  KEYSET k = keySetCreate () ;
  Graph  g = graphActive () ;
  
  displayCreate(DtKeySet) ;
  graphRetitle(title && *title ? title : "New KeySet") ; 
  keySetShow (k,0) ;
  keySetSelect () ;
  graphActivate (g) ;
}

/*static void newEmptyKeySet(void)
{ static int n = 0 ;
  newKeySet(messprintf("New %d",++n)) ;
}
*/
/***********************************/

BOOL importKeySet (KEYSET kSet, FILE *fil)
{
  int table, nbad = 0 ;
  int level = freesetfile (fil,"") ;
  char *word ;
  KEY key = 0;

  keySetMax(kSet) = 0 ;
  while (freecard (level))	/* closes file itself */
    if ((word = freeword()))
      { if (!(table = pickWord2Class (word)))
	{ ++nbad ; continue ; }
      if (table == _VKeySet && !keySetMax(kSet))
	continue ;
      freestep (':') ;
      if (!(word = freeword()))
	{ messout ("No object name, line %d", 
		   freestreamline (level)) ;
	continue ;
	}
      if (lexword2key (word, &key, table))
	keySet(kSet, keySetMax(kSet)) = key ;
      else
	++nbad ;
      }
  if (nbad)
    messout ("%d objects not recognized", nbad) ;
  keySetSort(kSet) ;
  keySetCompress(kSet) ;
  return TRUE;
}


/***********************************/

BOOL keySetDisplay (KEY key, KEY from, BOOL isOldGraph)
{
  KEYSET kset ;
  LOOK  look ;

  if (!(kset = arrayGet (key,KEY,"k")))
    return FALSE ;

  if (!isOldGraph)
    displayCreate (DtKeySet) ;
  else
    localDestroy () ;
  graphRetitle (name(key)) ;
  displayPreserve () ;
  
  look = keySetShow (kset, 0) ;
  look->keySetKey = key ;
  if (look)
    {
      return TRUE ; 
    }
  return FALSE ;
}

/***********************************/

/* BEWARE: look->curr (a box) runs 2..keySetMax+1, not 0..max-1 */
/* box 0 is the background, 1 is the "selected" message */

char *keySetCurrent (void *handle, KEY k)
{
  LOOK look = (LOOK) handle ;
  /* this is called only once from main pick, i removed the garbage, mieg 8/98*/
  if (look->magic != MAGIC)
    messcrash ("keySetCurrent called with bad look handle") ;
    
  if (look->showtype)
    return name(look->showtype) ;
  else
    return 0;
}

/*************************************/
/*************************************/

static void localDestroy (void)
{
  LOOKGET("keySetDestroy") ;

  look->magic = 0 ;
  if (selectedKeySet == look)
    selectedKeySet = 0 ;
  keySetDestroy (look->kSet) ;
  keySetDestroy (look->undo) ;
  keySetDestroy (look->keySetAlpha) ;
  messfree (look) ;

  graphAssRemove (&MAGIC) ;
}

/*************************************/

       /* selects keyset for mix, query etc */
void keySetSelect (void)
{
  LOOKGET ("keySetSelect") ;

  if (    selectedKeySet 
      &&  selectedKeySet != look
      && selectedKeySet->magic == MAGIC
      && graphExists(selectedKeySet->graph))
    { Graph hold = graphActive () ;
      graphActivate (selectedKeySet->graph) ;
      
      graphBoxDraw (1,BLACK,WHITE) ;    /* box 1 is the message */
      graphActivate (hold) ;
    }
  graphBoxDraw (1,BLACK,PALERED) ;
  selectedKeySet = look ;
}

/***************/

BOOL keySetActive(KEYSET *setp, void** lookp) 
{
  if (selectedKeySet && selectedKeySet->magic == MAGIC
      && graphExists(selectedKeySet->graph)
      && arrayExists(selectedKeySet->kSet))
    { if (setp) *setp = selectedKeySet->kSet ;
      if (lookp) *lookp = selectedKeySet ;
      return TRUE ;
    }
  selectedKeySet = 0 ;
  if (setp) *setp = 0 ; 
  if (lookp) *lookp = 0 ;
  return FALSE ;
}

/***************/

/* BEWARE: look->curr (a box) runs 2..keySetMax+1, not 0..max-1 */
/* box 0 is the background, 1 is the "selected" message */

static void localPick (int box,double x,double y)
{

  LOOKGET("keySetPick");

  if (look != selectedKeySet)
    keySetSelect () ;
  if (look->showQuery)
    {
      if (box == look->queryBox) 
	{
	  graphTextScrollEntry (look->queryBuffer,0,0,0,0,0) ;
	  look->showQuery = 2 ; /* keep focus */
	  return ;
	}
      else
	{
	  graphEntryDisable() ;
	  graphRegister (KEYBOARD, localKbd) ;
	  look->showQuery = 1 ;/* release focus */
	}
    }

  if (!box) 
    return ;
  if(box == look->upBox){
    pageUp();
    return;
  }
  else if(box == look->downBox){
    pageDown();
    return;
    }
  else if(box == look->scrollBarBox){
    graphBoxDrag(box, keyBandDrag);
    return ;
  }
  
  if (box == look->curr)                  /* double click */
    { KEY key = array(look->box2key, look->curr, KEY) ;
      if (key)
	{
	  if (look->isBlocking)
	    ksetDeleteKey (key) ;
	  else 
	    {  
	      graphPostBuffer (name(key)) ;
	      if (class (key) == _VKeySet)
		{
		  Graph g = graphActive () ;
		  KEYSET ks2 ;

		  displayCreate(DtKeySet) ;
		  graphRetitle (name (key)) ;
		  ks2 = arrayGet (key, KEY, "k") ;
		  keySetShow (ks2,0) ;
		  ks2 = 0 ; /* given to new display */
		  displayPreserve() ;
		  graphActivate (g) ;      
		}
	      else if(look->showtype != 0)
		display(key,0,name(look->showtype));
	      else
		display(key,0,0);
	    }
	}
      return;
    }
  else
    { 
      if(look->curr > 1)
	graphBoxDraw(look->curr,BLACK,LIGHTGRAY);
      look->curr=box;
      if(look->curr >1){
	startbox = box;
	controlLeftDown (x,y);
	if (look->curr < arrayMax(look->box2key))
	  graphPostBuffer (name(array(look->box2key, look->curr, KEY))) ;
      }
    }
}

/***********************************/

static void keySetMix(LOOK qry, int combinetype)
{
  KEYSET s=0;
  LOOK  q = selectedKeySet ;
  
  if(combinetype == AND) 
      s = keySetAND(q->kSet,qry->kSet) ;
  else if(combinetype  == OR)
      s = keySetOR(q->kSet,qry->kSet) ;
  else if(combinetype == MINUS)
      s = keySetMINUS(qry->kSet,q->kSet) ;
  else if(combinetype == XOR)
      s = keySetXOR(q->kSet,qry->kSet) ;
  else
    messcrash("keySetMix received wrong operator %d",combinetype) ;

  keySetDestroy (qry->undo) ;
  qry->undo = qry->kSet ;
  qry->kSet = s ;
  qry->touched |= 3 ;
  qry->curr = 0 ;

  keySetDestroy(qry->keySetAlpha) ;
  keySetShow2(qry->kSet,qry ) ;
}

/***********************************/

static void localKbd (int k)
{
  float x0,y0,x1,y1,x2,y2,x3,y3,x4,y4 ;
  int max, before ;
  KEY key;
  LOOKGET("keySetPick") ;

  if (look->curr < 2)
    return ;

  max = keySetMax (look->kSet) ;
  
  key = arr(look->box2key,look->curr, KEY);
  graphBoxDraw(look->curr, BLACK, LIGHTGRAY);

  switch (k)
    {
    case UP_KEY :
      if (arr(look->box2key, look->curr-1, KEY))
	--look->curr ;
      else if (look->base)
	{ pageUp () ;
	  for (look->curr = arrayMax(look->box2key)-1 ;
	       look->curr && !arr(look->box2key, look->curr, KEY) ;
	       --look->curr) ;
	}
      break ;
    case DOWN_KEY :
      if(arr(look->box2key,look->curr,KEY) !=  arrayMax(look->kSet)){ /* do not want to go down if this is the last key */
	if (arr(look->box2key, look->curr+1, KEY))
	  ++look->curr ;
	else if (look->top < max-1)
	  { pageDown () ;
	    for (look->curr = 1 ;
		 look->curr < arrayMax (look->box2key) && 
		 !arr(look->box2key, look->curr, KEY) ;
		 ++look->curr) ;
	  }
      }
      break ;
    case RIGHT_KEY :
      graphBoxDim (look->curr,&x0,&y0,&x2,&y2) ;
      graphBoxDim (look->lastBox,&x3,&y3,&x4,&y4);
      if(x0 != x3 || look->lastBox == look->curr){
	y1 = y0 ; x1 = x0 ;
	while (arr(look->box2key,look->curr+1,KEY) &&
	       (y1 < y0-0.5 || x1 == x0))
	  graphBoxDim (++look->curr,&x1,&y1,&x2,&y2) ;
	if ((y1 < y0 - 0.5 || x1 == x0) && look->top < max-1)
	  { pageDown() ;
	    for (look->curr = 1 ;
		 look->curr < arrayMax (look->box2key) && 
		 !arr(look->box2key, look->curr, KEY) ;
		 ++look->curr) ;
	  }
      }
      break ;
    case LEFT_KEY :
      graphBoxDim (look->curr,&x0,&y0,&x2,&y2) ;
      graphBoxDim (look->firstBox,&x3,&y3,&x4,&y4); 
      if(x3 != x0 || look->curr == look->firstBox){ /* i.e. the first column */
	y1 = y0 ; x1 = x0 ;
	before = look->curr;
	while (look->curr > 2 && (y1 > y0+0.5 || x1 == x0) &&look->curr > look->firstBox)
	  graphBoxDim (--look->curr,&x1,&y1,&x2,&y2) ;
	if (
	    ((y1 > y0+0.5 || x1 == x0) && look->base) ||
	    ( (before == look->firstBox) && 
	      (arr(look->box2key,before,KEY)!= arr(look->keySetAlpha,0,KEY))
	      )
	    )
	  { pageUp () ;
	    for (look->curr = arrayMax(look->box2key)-1 ;
		 look->curr && !arr(look->box2key, look->curr, KEY) ;
		 --look->curr) ;
	  }
	if(look->curr < look->firstBox)
	  look->curr = look->firstBox;
      }
      break ;
    case HOME_KEY :
      if (look->base)
	{ look->base =0;
	  keySetShow2 (look->kSet, look) ;
	}
      for (look->curr = 1 ;
	   look->curr < arrayMax (look->box2key) && 
	   !arr(look->box2key, look->curr, KEY) ;
	   ++look->curr) ;
      break ;
    default:
      graphBoxDraw(look->curr, WHITE, BLACK) ;
      return ;
    }

  key = arr(look->box2key,look->curr, KEY);
  if ( key && !keySetFind (look->selected, key, 0))
    { keySetInsert (look->selected, key) ;
      look->numselpage ++ ;
      look->touched |= 1 ;
    }
  graphBoxDraw (look->curr, WHITE, BLACK) ; 
  shownItem (look) ;

  if (look->txtbox)
    graphBoxDraw(look->txtbox,BLACK,WHITE);
  graphPostBuffer(name(key));
   { KEY key = array(look->box2key, look->curr, KEY) ;
     if(look->showtype != 0)
       display(key,0,name(look->showtype));
     else
       display(key,0,0);
   }
}

/***************************/

static void ksetAddKey (KEY key)
{
  LOOKGET("ksetAddKey") ;

  if (key) keySetInsert(look->kSet,key) ;
  look->touched |= 3 ;
  keySetDestroy(look->keySetAlpha) ;

  keySetShow2(look->kSet,look) ;

  displayRepeatBlock () ;
}

/***************************/

static void ksetDeleteKey (KEY key)
{
  LOOKGET("ksetDeleteKey") ;

  keySetRemove(look->kSet,key) ;
  look->touched |= 3 ;
  keySetRemove(look->selected,key) ; /* mhmp add-remove */
  keySetDestroy(look->keySetAlpha) ;

  look->curr = 0; /* otherwise if it is the last box it crashes */
  keySetShow2(look->kSet,look) ;

  displayRepeatBlock () ;
}

/*****************************/

static void ksetUnBlock (void)
{
  LOOKGET("ksetUnBlock") ;

  if (look->isBlocking)
    { look->isBlocking = FALSE ;
      displayUnBlock() ;
    }
}

/*************************************/

static void selectedToNewKeySet(void)
{
  Graph  g ;
  LOOKGET("selectedToNewKeySet") ;

  g = graphActive () ;
  displayCreate(DtKeySet) ;
  graphRetitle ("Copied selected keyset") ;
  keySetShow (keySetCopy (look->selected),0) ;
  displayPreserve() ;
  graphActivate (g) ;
}

static void ksetlocalMenu (KEY k)
{
  extern void Bdump (KEY key) ;
  Graph  g ;
  KEY	 key = 0 ;
  KEYSET ks2 = 0 ;
  FILE   *fil ;
  char *word = 0 ;
  static char fileName[FIL_BUFFER_SIZE] , dirName[DIR_BUFFER_SIZE] ;
  LOOKGET("localMenu") ;

  switch (k)
    {
    case 99: 
      graphDestroy () ; 
      break ;
    case 98: 
      helpOn("KeySet") ;
      break ;
    case 97: 
      graphPrint() ;
      break ;
    case 51:
      keySetSummary();
      break;
    case 1:	/* save */
      if (!look->keySetAlpha)
	break ;
      if(!checkWriteAccess())
	return ;

      word = look->keySetKey ? name(look->keySetKey) : "" ;
      while (!key && messPrompt ("Give a name",word,"t"))
	if (!(word = freeword()) || !*word ||
	    (!lexaddkey (word,&key,_VKeySet) &&
	     !messQuery ("Overwrite existing keyset ?")))
	  return ;
      
      if (!key)
	break ;
      
      if (!lexlock(key))
	{ messout ("Sorry, key is locked by someone else") ;
	  break ;
	}
      arrayStore (key,look->kSet,"k") ;
      lexunlock (key) ;
      look->keySetKey = key ;
      /*
      if (TRUE || messQuery ("Do you want to save the session now?"))
	sessionClose (TRUE) ;  
	save */
      break ;
    case 2:	/* copy */
      g = graphActive () ;
      displayCreate(DtKeySet) ;
      graphRetitle ("Copied keyset") ;
      keySetShow (keySetCopy (look->kSet),0) ;
      displayPreserve() ;
      graphActivate (g) ;
      break ;
     case 8: 
      if (look->curr > 1 && look->curr < arrayMax(look->box2key))
	{ key = array(look->box2key, look->curr, KEY) ;
	  if (key) 
	    display(key, 0, TREE);
	}
      break;
    case 10:	/* Read from file */
      if (!(fil = filqueryopen (dirName, fileName, "list", "r",
				"File of object names ?")))
	return ;

      g = graphActive () ;
      displayCreate(DtKeySet) ;
      graphRetitle (fileName) ;
      ks2 = keySetCreate () ;
      importKeySet(ks2, fil);
      keySetShow (ks2,0) ;
      ks2 = 0 ; /* givenb to new display */
      displayPreserve() ;
      graphActivate (g) ;
      break ;
    case 6:     /* Add key */
      if (look->isBlocking)
	{
	  ksetUnBlock() ;
	  break ;
	}
      look->isBlocking = TRUE ;
      displayPreserve() ; look->owner = -1 ;
      displayBlock (ksetAddKey,
		    "Double-clicking in a different window will add entries, "
		    "inside the current KeySet window will delete them from this KeySet.  "
		    "This will continue until you remove this message "
		    "or select again this option in the same menu") ;
      break ;
    case 7: 
      newKeySet (0) ;
      break ;
    case 41 :  /* Mailer */
      if (!look->keySetAlpha)
	break ;
      acedbMailer(0, look->kSet, 0) ;
      break ;
    case 42:     /* JTM Hexadecimal dump of disk blocks */
      if (!look->keySetAlpha)
	break ;
      if (look->curr < 2)
	break ;
      Bdump(keySet(look->keySetAlpha,look->curr-2+look->base)) ;
      break ;
    case 43:
      if (!look->keySetAlpha)
	break ;
      if (*mhmpDisplayp)
	(*mhmpDisplayp) (look->kSet, 0, 0) ;
      break ;
    case 44:
      if (look->curr > 1 && look->curr < arrayMax(look->box2key))
	{ key = array(look->box2key, look->curr, KEY) ;
	  if (key) 
	    display(key, 0, DtHSEQ) ;
	}
      break ; 
    case 101:
      if (look->curr > 1 && look->curr < arrayMax(look->box2key))
	{ key = array(look->box2key, look->curr, KEY) ;
	  if (key) 
	    display(key, 0, DtTiling) ;
	}
      break ; 
    case 45:
      if (!look->keySetAlpha)
	break ;
      doMultiMapDisplay (look->kSet) ;
      break ;
    case 46:
     if (!look->keySetAlpha)
       break ;
     forestDisplayKeySet (0, keySetCopy(look->kSet), FALSE) ;
     break ; 
    case 47:
      look->curr = 0;
      look->showClass = !look->showClass;
      keySetShow2(look->kSet,look) ;
      break ; 
    case 48:
     findneigh() ;
     break ;
    case 49:
      if (look->curr > 1 && look->curr < arrayMax(look->box2key))
	{ key = array(look->box2key, look->curr, KEY) ;
	  if (key) 
	    display(key, 0, DtGLOC) ;
	}
      break ; 
    case 60:
      if (look->curr > 1 && look->curr < arrayMax(look->box2key))
	{ key = array(look->box2key, look->curr, KEY) ;
	  if (key) 
	    display(key, 0, DtGeneExp) ;
	}
     break ;
   case OR:  
     beforecombine(OR);
     break ;
   case AND:  
     beforecombine(AND);
     break ;
   case MINUS:  
     beforecombine(MINUS);
     break ;
   case XOR:  
     beforecombine(XOR);
     break ;
    case UNDO:
	if (keySetExists (look->undo))
	  {
	    arrayDestroy(look->kSet) ;
	    look->kSet = look->undo ;	
	    look->undo = 0 ;
	    look->touched |= 3 ;
	    look->curr = 0 ;
	    
	    keySetDestroy(look->keySetAlpha) ;
	    keySetShow2(look->kSet, look) ;
	  }
	else
	  {
	    messout ("Sorry, no previous keyset available") ;
	    return ;
	  }
    default: 
      break ;
    }
}

/*************************************************/
/************* left button for thumb **********/

static void controlLeftDrag (double x, double y)
{
  int i,nx,ny,ix,numrows,incr,iy,col,colour,addbit ;
  int index;
  BOOL inlimits=TRUE;
  float right,ncol;
  KEY key;
  LOOKGET("controlLeftDrag");

  if(oldcol >= 0){
    graphFitBounds (&nx, &ny) ;
    
    iy = (int)y;
    numrows = ny-(TOPGAP+BOTGAP);
    col = 0;
    if(x < SCROLLBARGAP)
      x= (float)SCROLLBARGAP+1.0;

    right = (float)SCROLLBARGAP;
    ncol = ((float)look->numperpage/(float)numrows);
    i=0;
    while(i <= ncol){
      right+=(float)look->colWidth[i++]+2.0; 
    }
    if(x > right)
      x= right;
    ix = SCROLLBARGAP;
    while(ix < x){
      ix +=look->colWidth[col++] + 2;
    }
    
    if(iy  > ny-BOTGAP-1)
      iy = ny-BOTGAP-1;
    if(iy < TOPGAP-1)
      iy = TOPGAP-1;
    incr = iy-oldy;
    oldy=iy;
    incr += (col-oldcol)*numrows;
    oldcol = col;
    numboxes += incr;
    if(incr > 0)
      addbit = -1;
    else
      addbit =-2;
    if(incr != 0){
      while(lastbox != startbox+(numboxes+addbit) && inlimits){
	if((incr>0 && lastbox > startbox) || (incr<0 && lastbox < startbox) || lastbox == startbox)
	  colour = LIGHTGRAY;
	else{
	  key = arr(look->box2key,lastbox,KEY);
	  if(!keySetFind(look->selected,key,&index))
	    colour = WHITE;
	  else
	    colour = LIGHTGRAY;
	}
	if(lastbox >= look->firstBox && lastbox <=  look->lastBox)
	  graphBoxDraw(lastbox,BLACK,colour);
	else
	  inlimits = FALSE;
	if(incr > 0)
	  lastbox++;
	else
	  lastbox--;
      }
    }
  }
  else{
    col = 0;
    ix = SCROLLBARGAP;
    while(ix < x){
      ix +=look->colWidth[col++] + 2;
    }
    oldy = (int)y-1;
    oldcol = col;
  }    
}    

static void controlLeftUp (double x, double y)
{
  KEY key;
  int endbox=0 ;
  int index ;
  LOOKGET("controlLeftUp");

  look->curr = 0;
  if(numboxes > 1){
    numboxes--;
    endbox = startbox+numboxes-1;
  }
  if(numboxes < 1){
    numboxes--;
    endbox = startbox;
    startbox = startbox + numboxes;
  }
  else if(numboxes == 1){
    endbox = startbox;
    look->curr = startbox;
  }
  while(startbox <= endbox && startbox < arrayMax (look->box2key)){
    key = arr(look->box2key,startbox, KEY);
    if(key && !keySetFind(look->selected,key,&index))
      {
	keySetInsert(look->selected,key); 
	look->touched |= 1 ;
	look->numselpage++;
      }
    startbox++;
  }
  graphRegister (LEFT_DRAG, none) ; 
  graphRegister (LEFT_UP, none) ;

  shownItem (look) ;

  graphBoxDraw(look->txtbox,BLACK,WHITE);
  if(numboxes==1 && look->curr > 0){
    graphBoxDraw(look->curr,WHITE,BLACK);
    graphPostBuffer(name(arr(look->box2key,look->curr, KEY)));
  }
}

static void controlLeftDown (double x, double y) 
{ 

  LOOKGET("controlLeftDown");
  numboxes = 1;
  oldcol = -1;
  lastbox = startbox;
  graphBoxDraw(lastbox,BLACK,LIGHTGRAY);
  graphRegister (LEFT_DRAG, controlLeftDrag) ; 
  graphRegister (LEFT_UP, controlLeftUp) ;
}


static void nextShowtype(void)
{
  FREEOPT *menu;
  int i ;
  LOOKGET("changeShowtype");
  
  menu = pickGetDisplayMenu();
  for(i=1;i<=menu[0].key;i++)
    if(look->showtype == menu[i].key)
      { i++ ; if (i > menu[0].key) i = 1 ;
      look->showtype = menu[i].key ;
      keySetShow2(look->kSet,look) ;
      break ;
      }
}

static void changeShowtype(KEY key)
{
  LOOKGET("changeShowtype");
  
  look->showtype = key;
  keySetShow2(look->kSet,look) ;
}

#define MAXLEN 45
static void keySetSummary(void)
{
  int	i, *summary ;
  KEY	*kp ;
  char	*buf ;
  LOOKGET("keySetSummary") ;

  summary = (int *) messalloc ((sizeof(int) * MAX_CLASSES)) ;
  
  if (look->kSet)
    for (i = keySetMax(look->kSet),
	 kp = arrp (look->kSet, 0, KEY) - 1 ; i-- ; kp++)
      if (lexIsKeyVisible(*kp))
	summary[class(*kp)]++;
  
  buf = (char *) messalloc ((MAX_CLASSES+1)*MAXLEN) ;
  sprintf (buf,"KeySet Class Information\n");
  for (i = 0 ; i < MAX_CLASSES ; i++)
    if (summary[i] != 0)
      strcat (buf, messprintf ("%s %d items\n",
			       pickClass2Word(i), summary[i])) ;
  graphMessage (buf) ;
  messfree (buf) ;

  messfree (summary) ;
}

static BOOL columnWidth (int startpos, int direction,
			 int *cWidth, int *endpos, int numrows)
{
  int k=0,i,max,j;
  const char *cp,*cp2,*cp3;
  KEY key;
  LOOKGET("columnWidth");

  if (numrows < 1)
    return FALSE;

  i= startpos;
  *cWidth = 0;
  max = keySetMax(look->keySetAlpha);
  
  while(k<numrows){  
    if(i >= max){               /* check end if moving forward */
      return FALSE;
    }
    if(i < 0)                   /* check for start if going back a page */
      return FALSE; 
    cp = cp2 =0 ;
    key = arr(look->keySetAlpha,i,KEY) ;
    if (!nextName (key, &cp) || *cp == '\177'){
      if(direction) 
	*endpos = i++;
      else
	*endpos = i--;
      continue ;
    }  
    if(look->showClass){
      cp2 = className (key);	    
      cp3 = strnew(messprintf("%s:%s",cp2,cp),0);
      /*sprintf(cp3,"%s:%s\n",cp2,cp);*/
      if ((j = strlen(cp3)) > *cWidth){
	*cWidth = j ;
      }
      messfree(cp3);
    }
    else{
      if ((j = strlen(cp)) > *cWidth){
	*cWidth = j ;
      }
    }
    if(direction)
      *endpos = i++;
    else
      *endpos = i--;
    k++;
  }
  return TRUE;
}

static void drawScrollBar(int max)
{
  int nx, ny ;
  float xpt, top, bottom, boxstart, boxend, len ;
  LOOKGET("drawScrollBar");
  

  graphFitBounds (&nx, &ny) ;
  if(max > look->numperpage && max != 0){
    top  =  (float)TOPGAP;
    bottom = (float)(ny - BOTGAP);
    
    xpt = SCROLLBARPOS;
    graphLinewidth(0.4);
    graphLine(xpt,top,xpt,bottom); 
    graphLinewidth(0); /* mieg 2002 , used to be stupid int k */
    
    if(look->base !=0)
      look->upBox = scrollTriangle(UP,xpt-1.0,top,GREEN);
    else
      look->upBox =  scrollTriangle(UP,xpt-1.0,top,WHITE);
    if(look->top != max-1)
      look->downBox = scrollTriangle(DOWN,xpt-1.0,bottom,GREEN)  ;
    else
      look->downBox = scrollTriangle(DOWN,xpt-1.0,bottom,WHITE)  ; 

    len =  (((float)look->numperpage)/(float)max)*(bottom-top);
    if(len < MINSCROLLBARLENTH)
      len = MINSCROLLBARLENTH;
    look->scrolllen = len;

    if(look->base == 0)
      boxstart = top;
    else{
      if((int)look->scrolllen == MINSCROLLBARLENTH)
	boxstart = (((float)look->base/(float)max)*(bottom-top-look->scrolllen)) + top;
      else
	boxstart = (((float)look->base/(float)max)*(bottom-top)) + top;
    }
    boxend = boxstart + len;

    if(boxend > bottom){
      boxend =  bottom;
      boxstart = bottom-len;
    }

    look->scrollBarBox = graphBoxStart();

    graphRectangle(xpt-0.5,boxstart,xpt+0.5,boxend);
    graphBoxEnd();
    graphBoxDraw(look->scrollBarBox,BLACK,GREEN);
  } 
  else
    look->scrollBarBox=look->upBox=look->downBox =0; 
}

static void beforeForest (void)
{     
  LOOKGET ("beforeForest") ;

  if (!look->keySetAlpha || !look->kSet || !keySetMax(look->kSet))
    return ;
  forestDisplayKeySet (0, keySetCopy(look->kSet), FALSE) ;
}

static void otherTools (void)
{
  LOOKGET("otherTools") ;

/* changing the number of buttons changes the current selected */

  if(look->curr)
    look->curr -= (isOtherTools ? 1 : -1 ) * N3rdLINE ;
  isOtherTools = !isOtherTools ;
  keySetShow2(look->kSet,look) ;
}

static void beforeAddRemove(void)
{
  ksetlocalMenu (6) ;
}

static void beforeSaveAs(void)
{
  checkWriteAccess() ;
  if (!isWriteAccess())
    return ;
  ksetlocalMenu (1) ;
}

static void beforeEmpty(void)
{
  ksetlocalMenu (7) ;
}


static void beforeCopy(void)
{
  ksetlocalMenu (2) ;
}

static void drawMenus(void){
  int i,nx,ny,menuregion;
  float line ;
  static int wp = 0 ;
  FREEOPT *menu;
  KEYSET kSet;

  static MENUOPT dumpAs[] = {
    {menuSpacer,"Mail/Write data as:.."}, 
    {acedump,"  ace file"},
    {acedumptimestamp,"  ace file with timestamps"},
    {nameDump,"  list of names"},
    {fastaDump,"  DNA in FASTA format"},
    {protDump,"  Protein in FASTA format"},
    {alignDump,"  Protein alignment "},
    {0,0}
  };

  static MENUOPT combine[] = {
    {combineOR,"OR"},
    {combineAND,"AND"},
    {combineMINUS,"MINUS"},
    {combineXOR,"XOR"},
    {0,0}
  };

  static MENUOPT query[] = {
    {openQueryBox, "Query"},
    {follow,"Follow"},
    {findneigh,"Find Neighbours"}, /* redundant on morer info -> background ? */
    {queByExam,"Query by Example"},
    {grepset,"Text search"},
    {0,0}
  };
  
  /*  static MENUOPT selectMenu[] = { 
    {clearSelected,"Clear Selection"},
    {selectAllKeySet,"Select All"},
    {reverseSelectedKeySet,"Reverse"},
    {removeSelectedKeySet,"Remove Selected"},
    {selectedToNewKeySet,"Selected To New Keyset"},
    {0,0}
  };
  */

  static MENUOPT modifSelectMenu[] = {
    {beforeAddRemove, "Add/remove items one by one"}, 
    {menuSpacer,""},
    {beforeEmpty,"New empty keyset"},
    {beforeCopy,"Copy whole keyset"},    
    {beforeSaveAs,"Save current keyset as"},
    {menuSpacer,""},
    {selectAllKeySet,"Select All"},
    {clearSelected,"Clear Selection"},
    {reverseSelectedKeySet,"Reverse selection"},
    {removeSelectedKeySet,"Remove selected items"},
    {selectedToNewKeySet,"Copy selected items"},
    {0,0}
  };

  static MENUOPT editMenu[] = {  
    /*   {editAll,"Select and Edit all these objects"}, */
    /* {killAll,"Kill all these objects"}, */  
    {selectAllKeySet,"Select All"},
    {editSelected,"Edit the selected objects"},
    {killSelected,"Kill the selected objects"},
    {0,0}
  } ;

  LOOKGET("drawMenus");

  kSet = look->kSet;
  if (!wp) wp = writeAccessPossible() ? 1 : 2 ;
  graphRegister(LEFT_DOWN, controlLeftDown);
  graphClear () ;
  graphFitBounds (&nx, &ny) ;
  
  menuregion = graphBoxStart();
  
  /****  first line of buttons ****/
  line = .5 ;
  /* i change the box to be the show as, rather than the text 
  char str1[15];
  graphText("Show As:", 1.0 ,line + .1 );
  menu = pickGetDisplayMenu();
  for(i=1;i<=menu[0].key;i++){
    if(look->showtype == menu[i].key){
      sprintf(str1,"%s..",menu[i].text);
      break;
    }
  }
  i = graphButton(str1,none,10.0,line) ;
  graphBoxFreeMenu(i,(FreeMenuFunction)changeShowtype,menu);
  */
  i = graphButton("Show As..:", nextShowtype, 1.0 ,line); 
  menu = pickGetDisplayMenu();
  graphBoxFreeMenu(i,(FreeMenuFunction)changeShowtype,menu);
  for(i=1;i<=menu[0].key;i++)
    if(look->showtype == menu[i].key)
      {
        graphText(menu[i].text, 13, line + .2) ;
	break;
      }
  
  if (biblioPossible(kSet))
    graphButton("Biblio",relatedBiblio,23.5,line);
    
  graphButton("More Info",beforeForest,31.1,line) ;

  graphButton("Quit",graphDestroy,42.0,line);

  /****  second line of buttons ****/
  line = 2.0 ;

  i = graphButton("Query..",openQueryBox,1.0,line);
  graphBoxMenu(i,query);

  i=graphButton("Select/Modif..",clearSelected,9.9,line) ; /* was 0.5); */
  graphBoxMenu(i,modifSelectMenu) ;

  i = graphButton("Export..",none,25.5,line) ; 
  graphBoxMenu(i,dumpAs);

  if (0) graphButton("Other",otherTools,35.1,line) ; /* was 0.5); */

  graphButton("Help",help,42.0,line);

  /****  optional third line of buttons ****/
  N3rdLINE = wp ? 3 : 2 ;  /* number of boxes of third line */
  if (isOtherTools)
    { 
      /* ATTENTION, adjust N3rdLINE if you edit this peice of code
       * ++ it for each additional box */
      line = 3.5 ;
  
      i=graphButton("Combine..",none,1.0,line); 
      graphBoxMenu(i,combine);

      graphButton(look->showClass ? "Hide Class" : "Class Names" ,showClasses,12.0,line) ;

      if (wp)
      i=graphButton("Edit/Kill..",editSelected,25.5,line); 
      graphBoxMenu(i,editMenu);
    } 
  NQueryLine = 2 ; /* number of boxes of query line */
  if (look->showQuery)
    { 
      int nx, ny ; graphFitBounds (&nx, &ny) ;
      line += 1.7 ;
      graphText ("Query: ", 3, line) ;
      look->queryBox = graphTextScrollEntry (look->queryBuffer,254, nx-10, 9, line, keySetDoQuery) ;  
    }
  if (look->showQuery != 2) /* focus desired */
    { graphEntryDisable() ;
    graphRegister (KEYBOARD, localKbd) ;
    }

  graphLine(0.0,0.0,(float)nx+1.0,0.0);
  line += 1.6 ;
  graphLine(0.0,line,(float)nx+1.0,line);
  graphBoxEnd();

  if (look == selectedKeySet)
    graphBoxDraw (menuregion,BLACK,PALERED) ;
  else
    graphBoxDraw (menuregion,BLACK,WHITE) ;

  TOPGAP = line + 2.4 ;
}

static void scrollMiddleUp(double x, double y)
{
  float newy;
  LOOKGET("scrollMiddleUp");
  
  if(look->keySetAlpha && look->scrollBarBox){
    newy = (float)y-(look->scrolllen/2.0);
    setPage(&newy);

    if(look->base == keySetMax(look->keySetAlpha)){ 
      --look->base;
      pageUp();
      return;
    }
    if (look != selectedKeySet)
      keySetSelect () ; 
    keySetShow2 (look->kSet, look) ; 
  }
}
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
