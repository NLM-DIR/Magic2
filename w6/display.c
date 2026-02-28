/*  File: display.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 **       Class sensitive object display                           **
          Shares and reuse display windows.
 * Exported functions:
   display 
     displays a Key according to diplayType 
     if displayType==0, it defaults according to classes.wrm.
     Each display type reuses its own window until that window has 
       been preserved with displayPreserve.

   displayBlock
     puts up a message over the current active window and blocks 
       the display system until a key has been picked somewhere.
     the rest of the program still runs.
     The registered function receives as parameter the picked key 
       and, implicitly, the graph it was picked from, which is 
       active.

   displayUnBlock
     Clears out.

   displayRepeatBlock
     Reactivates the current blocking function.

   isDisplayBlocked
     should be called by clients to prevent access to blocked displays

   displayPreserve
     preserves the object, so that the next object of that class
     will not reuse its window.
 * HISTORY:
 * Last edited: Nov 24 15:21 1998 (fw)
 * Created: Sun Dec  1 17:29:54 1991 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: display.c,v 1.18 2019/02/26 05:01:00 mieg Exp $ */

#include "acedb.h"
#include "aceio.h"
#include "lex.h"
#include "graph.h"
#include "display.h"
#include "sysclass.h"
#include "pick.h"
#include "client.h"
#include "query.h"
#include "freeout.h"
#if defined(applec)
#include "MPWIncludes.h"
#endif

/* KEY isGifDisplay = FALSE ; now in w2/graphAcedbInterface.c,  used to unset useless buttons etc */

#define MAXNDISPLAY 64
static int nDisplays = 0 ;

static Array  previous = 0 ; 
Array displayFuncArray = 0 ;
static Graph blockGraph = 0 ;
static BOOL isRepeatBlock = FALSE ;
static BlockFunc blockFunc ;
BOOL prefValue(char * name);

extern void (*gifEntry)(void*, int, BOOL) ; /* use as flag for giface */
BOOL displayReportGif = TRUE ;

/*********************************************************************/
/* displayName is now the name of an element of the display sysclass */

BOOL display (KEY key, KEY from, char *displayName)
{ KEY displayKey ;
  Graph g = graphActive() ;	/* save */
  BOOL gReuse ;
  DisplayFunc ff ;

  if (!previous)
    previous = arrayCreate (32, Graph) ;

  if (graphActivate (blockGraph))
    { (*blockFunc) (key) ;
      if (isRepeatBlock)
	isRepeatBlock = FALSE ;
      else
	{ 
	  if (graphActivate (blockGraph))
	    graphUnMessage () ;
	  blockGraph = 0 ;
	}
      graphActivate (g) ;	
      return FALSE ;
    }

  if (externalServer) externalServer (key, 0, 0, TRUE) ;/* get object and neighbours */

  if (iskey(key) != 2)
    return FALSE ;

  if (nDisplays >= MAXNDISPLAY)
    { messout ("Too many displayed objects - delete some first") ;

      return FALSE ;
    }
  
  if (!KEYKEY(key) && pickType(key) == 'B')
    displayName = "TREE" ; /* display must be TREE for the  model */
  
  if (!displayName)
    displayKey = pickDisplayKey (key) ;
  else
    if (!lexword2key (displayName, &displayKey, _VDisplay))
      { messerror ("Display type %s unknown, please edit wspec/displays.wrm",
		   displayName ? displayName : "NULL") ;
	return FALSE ;    /* not displayable */
      } 

  if (!displayKey) 
    return FALSE ;    /* not displayable */

  if ((g = array(previous,KEYKEY(displayKey), Graph)))
    { if (graphActivate(g))
	{ gReuse = TRUE ;
	}
      else
	{ array(previous,KEYKEY(displayKey), Graph) = 0 ;
	  gReuse = FALSE ;
	}
    }
  else 
    gReuse = FALSE;

  if ((ff = array(displayFuncArray,KEYKEY(displayKey), DisplayFunc)) &&
      ff (key,from,gReuse))
    { array(previous,KEYKEY(displayKey), Graph) = graphActive() ;
      graphPop () ;
      if (0 && gifEntry && displayReportGif)
	freeOutf ("// %s %s \"%s\"\n", 
		  name(displayKey), className(key), name(key));
      return TRUE ;
    }
  else
    return FALSE ;
}
  
/**********************************/

void displayPreserve (void)
{
  Graph g = graphActive() ;
  int i ;

  if (!previous)
    previous = arrayCreate (32, Graph) ;
  i = arrayMax(previous) ;
  while (i--)
    if (array(previous,i, Graph) == g)
      array(previous,i, Graph) = 0 ;
}

/**************************************/

BOOL displayBlock (BlockFunc func, char *message)
{
  static int n = 0 ;
  messStatus ("Select by double-clicking") ; /* mhmp nov 97 */
  if (graphExists(blockGraph))  /* mieg june 1 */
    return FALSE ;
  blockGraph = graphActive() ;
  isRepeatBlock = FALSE ;
  blockFunc = func ;
/*  if (!getenv("ACEDB_NO_GRAB_WINDOW"))*/
  if(!prefValue("NO_MESSAGE_WHEN_DISPLAY_BLOCK") && n++ < 3)
    {
      graphMessage (messprintf ("%s\n\n%s",
			      "Select an object by double-clicking as "
			      "though you were going to display it.",
			      message ?  message : "")) ;
    }
  return TRUE ;
}

/**************************************/

void displayRepeatBlock (void)
{
  isRepeatBlock = TRUE ;
}

/**************************************/

void displayUnBlock(void)
{
  if (blockGraph)		/* prevent callback recursion */
    { blockGraph = 0 ;
      graphUnMessage () ;
    }
  isRepeatBlock = FALSE ;
  blockFunc = 0 ;
}

/**************************************/

BOOL isDisplayBlocked (void) 
{
  return (blockGraph != 0) ;
}

/******************************************/
/******************************************/

/**** next section supports some simple and generic display types ***/

#include "bs.h"
#include "systags.h"

BOOL displayAs (KEY key, KEY from, BOOL isOldGraph)
{ 
  OBJ obj ;
  KEY tag, newkey ;

  if (!(obj = bsCreate (key)))
    return FALSE ;

  if (bsFindTag (obj, str2tag ("Display_as")) &&
      bsGetKeyTags (obj, _bsRight, &tag) &&
      bsGetKey (obj, _bsRight, &newkey))
    { bsDestroy (obj) ;
      display (newkey, from, 0) ;
      return TRUE ;
    }
  
  bsDestroy (obj) ;
  display (key, from, TREE) ;
  return FALSE ;
}


static int escapeUrl(char *old, char *new)
{ int len = 0;
  char c;
  while ((c = *old++))
    { 
      switch (c) 
	{ 
	case ' ':
	  c = '+'; /* fall through */
	default:
	  if (new) *new++ = c;
	  len++;
	  continue;
	case '?':
	case '=':
	case '+':
	case '/':
	case '#':
	case '&':
	case '%':
	  if (new)
	    { sprintf(new, "%%%2.2X", c);
	      new += 3; 
	    }
	  len += 3;
	  continue;
	}
    }
  if (new) *new = 0;
  return len+1; /* +1 for null terminator */
}
	


BOOL doWWWDisplay (KEY key, KEY from, int recursecount)
{
  OBJ obj, fromObj ;
  KEY template = 0, ref_key;
  char *url = 0, *ref_tag = 0 ;
  AC_HANDLE handle;
  char *str, *new, *cname,
    *urlsave=0;			/* init for compiler happiness */
  BOOL result;
  BSMARK mark = 0;
  Array a;
  int i;

  if (recursecount == 20)
    { messerror("Wild recursion detected in wwwDisplay!");
      return FALSE;
    }

  if (!strncmp(name(key),"http",4))   /* mieg, no nonsense code 
					 directly use name of URL object
				      */
     return graphWebBrowser(name(key)) ;

  if (!(obj = bsCreate (key)))
    return FALSE;

  if (bsGetData (obj, str2tag ("Url"), _Text, &url)) 
    { bsDestroy (obj);
      return graphWebBrowser(url);
    }
  
  if (bsGetKey (obj, str2tag ("Web_location"), &template))
    { bsDestroy(obj);
      return doWWWDisplay(template, key, recursecount+1);
    }
  
  if (from==key || !(fromObj = bsCreate(from)))
    { bsDestroy(obj);
      return FALSE;
    }
  
  url = 0;
  ref_key=0;
  
  if (bsGetData (obj, str2tag ("Reference_tag"), _Text, &ref_tag))
    do { /* if there's a class limitation on this tag, enforce it */
      mark = bsMark(obj, mark);
      if (!bsGetData(obj, _bsRight, _Text, &cname) ||
          strcmp(cname, className(from)) == 0)
        { lexaddkey(ref_tag, &ref_key, 0);
          if (bsGetData(fromObj, ref_key, _Text, &url))
            break;
          else
            ref_key = 0; /* so ref_key is valid when we match */
        } 
      bsGoto(obj, mark);
    }
    while (bsGetData(obj, _bsDown, _Text, &ref_tag));
  
  bsDestroy(fromObj);

  if (!ref_key && bsFindTag(obj, str2tag("Use_name")))
    { 
      if (!bsGetData(obj, _bsRight, _Text, &cname))
        url = name(from); /* no class constraint */
      else
        do { /* search for matching class */
          if (strcmp(cname, className(from)) == 0)
            { url = name(from);
              break;
            }
        } while (bsGetData(obj, _bsDown, _Text, &cname));
    }

  if (!url)
    { bsDestroy(obj);
      return FALSE;
    }
  /* ref_key may be zero here, is we Use_name */

  handle = handleCreate();
  url = strnew(url, handle); /* cause we're going to mung it */
  a = arrayHandleCreate( 10, BSunit, handle);

  if (bsFindTag (obj, str2tag("Rewrite")) && bsFlatten(obj, 4, a))
    for( i=0; i < arrayMax(a); i +=4)
      {  
        if ((str = arr(a, i+2, BSunit).s) &&
            strlen(url) >= strlen(str) &&
            strncmp(url, str, strlen(str)) == 0 )
          { urlsave = url; /* in case we fail on the next test */
            url += strlen(str);
          }
        else
          if (str) continue; /* no string matches */
      
        if ((str = arr(a, i+3, BSunit).s) &&
            (strlen(url) >= strlen(str) &&
             strcmp(url+strlen(url)-strlen(str), str) == 0))
          *(url+strlen(url)-strlen(str)) = 0;
        else
          if (str) /* no string matches */
            { url=urlsave;/* restore */
              continue; 
            }
        
        /* now url has pref/postfixes removed, escape url-dangerous chars */
        /* before adding new ones, only do this when using object name or  */
        /* reference tag, not a hard-specified URL */
        new = (char *)halloc(escapeUrl(url, 0), handle);
        (void)escapeUrl(url, new);
        url = new;
        
        /* now add prefix and postfix */
        new = (char *)halloc(1+strlen(url)+
                             strlen(arr(a, i, BSunit).s)+
                             strlen(arr(a, i+1, BSunit).s),
                             handle);
        strcpy(new, arr(a, i, BSunit).s);
        strcat(new, url);
        strcat(new, arr(a, i+1, BSunit).s);
        result = graphWebBrowser (new) ;
        goto ok;
      }
  
 result = FALSE;
ok:
  bsDestroy(obj);
  handleDestroy(handle);
  return result;
  
}

BOOL wwwDisplay (KEY key, KEY from, BOOL isOldGraph)
{ return doWWWDisplay( key, from, 0); }


/******************************************/
/********* Display Customization  *********/
/******************************************/

/* This list must be constructed by hand and match Graph.h */

static FREEOPT graphTypeOptions[] =
{
  {6, "graph types"},
  {PLAIN, "PLAIN"},
  {TEXT_SCROLL, "TEXT_SCROLL"},
  {TEXT_FIT, "TEXT_FIT"},     
  {MAP_SCROLL, "MAP_SCROLL"},
  {PIXEL_SCROLL, "PIXEL_SCROLL"},
  {TEXT_FULL_SCROLL, "TEXT_FULL_SCROLL"}
} ;

#define MAXDISPMENU 64 
static FREEOPT dispOptMenu[MAXDISPMENU] = {
  {1, ""},
  {0, "Graph"},   /* shorter name, not too inacurate */
  {0, "Text"},
};

static FREEOPT displayOptions[] =
{
  {8, "display options"},
  {'g', "GraphType"},
  {'t', "Title"},
  {'w', "Width"},
  {'h', "Height"},
  {'x', "Xposition"},
  {'y', "Yposition"},
  {'m', "Menu"},
  {'a' , "Help"}
  } ;

/******************************************************************/

static Array dArray = 0 ;

/******************************************************************/

Graph displayCreate(char *displayName)
{
 AcedbDisplay *dp = 0 ;
 Graph g ;
 float x, y ;
 static Array a = 0 ;
 int i, d ;
 KEY key = 0 ;

 if (!displayName)
   return 0 ;

 lexword2key(displayName, &key, _VDisplay) ;
 d = KEYKEY (key) ;
 if (d == ZERO)
   return 0 ;

 if (d > 0 && d < arrayMax(dArray))
   dp = arrp(dArray, d, AcedbDisplay) ;

 if (!dp || dp->width <= 0)
   { messerror ("%s%s%s\n%s\n%s", 
		"The geometry of display ", displayName, 
		" is missing", 
		"Please edit the file  wspec/display.wrm",
		"Then quit and restart acedb") ;
     return 0 ;
   }
 if (!*dp->title)
   strcpy (dp->title, "Sorry, no title in display.wrm") ;
 
 /***** A trick to prevent window superposition ***/
 if (!a)
   a = arrayCreate(20, Graph) ;
 for (i = 0 ; i < arrayMax(a) ; i++)
   if (!graphExists(arr(a,i,Graph)))
     break ;

 x = dp->x + i/30.0 ;
 while (x > 1.0) x -= 1.0 ; /* prevent overflow HJC */

 y = dp->y + i/30.0 ;
 while (y > 1.0) y -= 1.0 ;

  /* This is useful if the window manager does not
   * prompt the user, and is called by displayCreate
   */

 g =  graphCreate (dp->type , dp->title, 
		x, y, dp->width, dp->height) ;

 if (g)
   graphHelp(dp->help) ;
 
 array(a,i,Graph) = g ;

 return g ;
} /* displayCreate */

/**************************************/

void pickDefaultDisplays (void)
{
  KEY key ;
  int   n, d ;
  AcedbDisplay *dp ;
  const char **type, *types[] = {"TREE", "DtKeySet", "DtLongText", "DtBqlDisp", 0} ;
  float xx[]            = {  .001, .1, .12, 0.2, 0} ;
  float yy[]            = {  .001, .1, .12, 0.2, 0} ;
  float ww[]            = {  .55,  .2, .30, 1.9, 0} ;
  float hh[]            = {  .3,   .3, .30, 0.9, 0} ;
  
  dArray = arrayReCreate (dArray, 12, AcedbDisplay) ;

  lexaddkey("ZERO", &key, _VDisplay) ;
  d = KEYKEY (key) ;
  dp = arrayp(dArray, d, AcedbDisplay) ;
  dp->type = 0 ;  

  for (type = types, n = 0 ; *type ; type++, n++)
    {
      lexaddkey(*type, &key, _VDisplay) ;
      d = KEYKEY (key) ;
      dp = arrayp(dArray, d, AcedbDisplay) ;
      dp->type = TEXT_FULL_SCROLL ;
      dp->x = xx[n] ; dp->y = yy[n] ; dp->width = ww[n] ; dp->height = hh[n] ;
      strncpy(dp->title, *type, 31) ;
      strncpy(dp->help, *type, 31) ;
    }
} /* pickDefaultDisplays */


FREEOPT *pickGetDisplayMenu(void)
{ 
  return dispOptMenu;
} /* pickGetDisplayMenu */

void pickGetDisplayTypes(void)
{
  AcedbDisplay *dp ;
  int d ;
  int line = 0 ;
  char *cp ; float x ;
  KEY option, displayKey, treeKey;
  FILE *fil ;
  AC_HANDLE h = ac_new_handle () ;
  /*   ACEIN ai = aceInCreate ("wspec/displays.wrm", FALSE, h) ; */

  dispOptMenu[0].key = 2; /* since this gets called twice... */
  lexaddkey("TREE", &treeKey, _VDisplay);
  dispOptMenu[2].key = treeKey ; /* provide a tree default unless overridden */
    
  fil = filopen("wspec/displays", "wrm", "r") ;
  if (!fil)
    messExit("Cannot find display definition file : wspec/displays.wrm");

  
  /* Defaults */
  freespecial ("\n\t\"\\") ;
  while(line ++ ,freeread(fil))
    /* read _D... */
    if ((cp = freeword()) && *cp++ == '_' && *cp++ == 'D')
      {
	if(!*cp)
	  messExit("Error parsing line %d of wspec/displays.wrm", line);
	
	lexaddkey (cp, &displayKey, _VDisplay) ;
	d = KEYKEY (displayKey) ;
	if (d > arrayMax(dArray))
	  { dp = arrayp(dArray, d, AcedbDisplay) ;
	    dp->type = TEXT_FIT ;
	    dp->x = .001 ; dp->y = 0.001 ; dp->width = .3  ; dp->height = .5 ;
	    strncpy(dp->title, "acedb", 31) ;
	    strncpy(dp->help, "No_help", 31) ;
	  }
	
	dp = arrayp(dArray, d, AcedbDisplay) ;
	
	while (TRUE)
	  {
	    freenext() ;
	    if (!freestep('-')) 
	      {
		if ((cp = freeword()))
		  messcrash ("In displays.wrm line %d, no - at start of option %s", line, freepos()) ;
		else
		  break ;
	      }
	    
	    freenext() ;
	    if (!freekey(&option, displayOptions))
	      messcrash ("In displays.wrm line %d, unknown option %s", line , freepos()) ;
	    
	    switch (option) 
	      {
	      case 'm':
		cp = freeword();
		if (dispOptMenu[0].key < MAXDISPMENU-1)
		  { int where;
		    if (displayKey == treeKey)
		      where = 2;
		    else
		      where = ++dispOptMenu[0].key;
		    dispOptMenu[where].key = displayKey;
		    dispOptMenu[where].text = strnew(cp, 0);
		  }
		break;
	      case 't':
		cp = freeword() ;
		strncpy(dp->title, cp, 41) ;
		    break ;
	      case 'a':
		if ((cp = freeword()))
		  strncpy(dp->help, cp, 31) ;
		break ;
	      case 'g':                       /* GraphType (enum) */
		freenext() ;
		if (!freekey(&option, graphTypeOptions))
		  messcrash("Bad graph type  in line %d of displays at %s ", 
			    line,  freepos()) ;
		dp->type = option ;
		break ;
	      case 'x':
		if (!freefloat(&x))
		  messcrash("Missing float x value in line %d of wspec/displays.wrm at %s",
			    line,  freepos()) ;
		if (x < 0 || x > 1.3)
		      messout ("In wsepc/displays.wrm, line %d, x value %f out of range 0 1.3",
			       line, x ) ;
		else
		  dp->x  = x ;
		break ;
	      case 'w':
		if (!freefloat(&x))
		  messcrash("Missing float width value in line %d of wspec/displays.wrm at %s",
				line,  freepos()) ;
		if (x < 0 || x > 1.99)
		  messout ("In wsepc/displays.wrm, line %d, width value %f out of range .0 1.3",
			   line, x ) ;
		else
		  dp->width  = x ;
		break ;
	      case 'y':
		if (!freefloat(&x))
		  messcrash("Missing float y value in line %d of wspec/displays.wrm at %s",
			    line,  freepos()) ;
		if (x < 0 || x > 1.0)
		  messout ("In wsepc/displays.wrm, line %d, y value %f out of range 0 1.0",
			   line, x ) ;
		else
		  dp->y  = x ;
		break ;
	      case 'h':
		if (!freefloat(&x))
		  messcrash("Missing float height value in line %d of wspec/displays.wrm at %s",
			    line,  freepos()) ;
		if (x < 0 || x > 1.0)
		  messout ("In wsepc/displays.wrm, line %d, height value %f out of range 0 1.0",
			   line, x ) ;
		else
		  dp->height  = x ;
		break ;
	      }
	  }
      }
  filclose ( fil ) ;
  ac_free (h) ;
} /* pickGetDisplayTypes */



static void pickSetGraphType (char *cp, int type)
{ KEY d, displayKey ;
  AcedbDisplay *dp ;

  lexaddkey(cp, &displayKey, _VDisplay) ;
  d = KEYKEY (displayKey) ;
  if (d > arrayMax(dArray))
    { dp = arrayp(dArray, d, AcedbDisplay) ;
      dp->x = .001 ; dp->y = 0.001 ; dp->width = .3  ; dp->height = .5 ;
      strncpy(dp->title, "acedb", 31) ;
      strncpy(dp->help, "No_help", 31) ;
    }
  dp = arrayp(dArray, d, AcedbDisplay) ;
  dp->type = type ;
} /* pickSetGraphType */


void pickSetGraphTypes (void)
/* reset here the type of all graphs, since this is a code issue */
{
  pickSetGraphType ("DtMain", TEXT_FIT) ;
  pickSetGraphType ("DtChrono", TEXT_SCROLL) ;
  pickSetGraphType ("DtFile_Chooser", TEXT_SCROLL) ;
  pickSetGraphType ("DtHelp", TEXT_FULL_SCROLL) ;
  pickSetGraphType ("DtAce_Parser", TEXT_SCROLL) ;
  pickSetGraphType ("DtDump", TEXT_SCROLL) ;
  pickSetGraphType ("DtKeySet", TEXT_FIT) ;
  pickSetGraphType ("DtSession", TEXT_FULL_SCROLL) ;
  pickSetGraphType ("DtStatus", TEXT_SCROLL) ;
  pickSetGraphType ("DtUpdate", TEXT_SCROLL) ;
  pickSetGraphType ("DtLongText", TEXT_FULL_SCROLL) ;
  pickSetGraphType ("DtQuery", TEXT_SCROLL) ;
  pickSetGraphType ("DtBqlDisp", TEXT_FULL_SCROLL) ;
  pickSetGraphType ("DtQueryBuilder", TEXT_FULL_SCROLL) ;
  pickSetGraphType ("DtQueryByExample", TEXT_FULL_SCROLL) ;
  pickSetGraphType ("DtSpreadSheet", TEXT_FULL_SCROLL) ;
  pickSetGraphType ("DtMULTIMAP", TEXT_FIT) ;
  pickSetGraphType ("DtAction", TEXT_FIT) ;
  pickSetGraphType ("TREE", TEXT_SCROLL) ;
	pickSetGraphType ("DtBiblio", TEXT_FIT) ; /*  mhmp 03.12.02 */
  pickSetGraphType ("CMAP", TEXT_FIT) ;
  pickSetGraphType ("FMAP", TEXT_FIT) ;
  pickSetGraphType ("PEPMAP", TEXT_FIT) ;
  pickSetGraphType ("DtColControl", TEXT_FIT) ;
  pickSetGraphType ("DtGel", TEXT_FIT) ;
  pickSetGraphType ("GMAP", TEXT_FIT) ;
  pickSetGraphType ("VMAP", TEXT_FIT) ;
  pickSetGraphType ("GRID", TEXT_FIT) ;
  pickSetGraphType ("PMAP", TEXT_FIT) ;
  pickSetGraphType ("tMULTIMAP", TEXT_FIT) ;
  pickSetGraphType ("DtPmapPadSheet", TEXT_FULL_SCROLL) ;
  pickSetGraphType ("DtPmapFingerprint", TEXT_FULL_SCROLL) ;
  pickSetGraphType ("DtAlign", TEXT_FULL_SCROLL) ;
  pickSetGraphType ("DtCodons", TEXT_FIT) ;
  pickSetGraphType ("DtDnaTool", TEXT_SCROLL) ;
  pickSetGraphType ("DtAlias", TEXT_FIT) ;
  pickSetGraphType ("DtPlotPolygon", TEXT_FIT) ;
  pickSetGraphType ("DtHistogram", TEXT_FIT) ;
  pickSetGraphType ("DtMultiTrace", TEXT_FIT) ;
  pickSetGraphType ("DtHSEQ", TEXT_FIT) ;
  pickSetGraphType ("DtTiling", TEXT_FIT) ;
  pickSetGraphType ("DtGLOC", TEXT_FIT) ;
  pickSetGraphType ("DtGLOCBIG", TEXT_FIT) ;
  pickSetGraphType ("DtFiche", TEXT_SCROLL) ;
  pickSetGraphType ("DtImage", PIXEL_SCROLL) ;
  pickSetGraphType ("DtAlignment", TEXT_FULL_SCROLL) ;
  pickSetGraphType ("DtDendrogram", TEXT_FULL_SCROLL) ;
  pickSetGraphType ("DtGeneExp",  TEXT_FIT) ;

  return;
} /* pickSetGraphTypes */


char *pickMainTitle(void)
{ 
  
  static char *title = 0 ;

  if (!title)
    {
      KEYSET ks = query (0, "Find clone MainTitle") ;
      char *cp = 0 ;

      if (keySetMax (ks))
	{
	  OBJ Clone = 0 ;
	  if ((Clone = bsCreate (keySet(ks,0))))
	    {
	      bsGetText (Clone, str2tag ("MainTitle"), &cp) ;
	      bsDestroy (Clone) ;
	    }
	  keySetDestroy (ks) ;
	}
      if (!cp)
	{ 
	  int d ;
	  KEY key = 0 ;

	  lexword2key ("DtMain", &key, _VDisplay) ;
	  d = KEYKEY (key) ;
	  
	  cp = arrp(dArray, d, AcedbDisplay)->title ;
	}
      
      title = cp ? strnew (cp, 0) : "" ;
    }
  return title ;
} /* pickMainTitle */



/***** following needed for giface *****/
void pickSetDisplaySize (const char *cp, float x, float y, float w, float h)
{ KEY d, displayKey ;
  AcedbDisplay *dp ;

  lexaddkey(cp, &displayKey, _VDisplay) ;
  d = KEYKEY (displayKey) ;
  if (d > arrayMax(dArray))
    { dp = arrayp(dArray, d, AcedbDisplay) ;
      dp->x = .001 ; dp->y = 0.001 ; dp->width = .3  ; dp->height = .5 ;
      strncpy(dp->title, "acedb", 31) ;
      strncpy(dp->help, "No_help", 31) ;
    }
  dp = arrayp(dArray, d, AcedbDisplay) ;
  dp->width = w ;
  dp->height = h ;
  if (x > 0) dp->x = x ;
  if (y > 0) dp->y = y ;

  return;
} /* pickSetDisplaySize */

void pickRememberDisplaySize (const char *display)
{
  float sx, sy, sw, sh ;
  if (graphWindowSize (&sx, &sy, &sw, &sh))
    pickSetDisplaySize (display, sx, sy, sw, sh) ;
} /* pickRememberDisplaySize */


/*************************************************************/
/**************** Interactive class control ****************/
/*************************************************************/

#ifdef CLASSCONTROL

static KEY classChosen = 0 ;
static Array classTree = 0 ;
static Graph classGraph = 0 ;

static void localDestroy(void)
{ 
  classGraph =  0 ;
  arrayDestroy(classTree) ;
}

/************************************************************/  

static void classPick(KEY box)
{ 

/*
static int previousBox = 0 ;
  int i ;
  CT *ct ;

  for (i = 0, ct = arrp(classTree,0,CT) ; i < arrayMax(classTree) ; ct++, i++)
    if (ct->box == previousBox)
      graphBoxDraw(ct->box, BLACK, ct->color) ;
      	
  for (i = 0, ct = arrp(classTree,0,CT) ; i < arrayMax(classTree) ; ct++, i++)
    if (ct->box == box)
      { 
	  classChosen = ct->key ;
	
	graphBoxDraw(ct->box, BLACK, RED) ;
	if (ct->filter)
	  { filterPtr = stackText(filterStack, ct->filter) ;
	    graphBoxDraw(filterBox, -1, -1) ;
	  }
      }
  previousBox = box ;
*/
}

/************************************************************/  

static void classDestroy(void)
{
}

/****************************************************************/

static void classRename(void)
{ 
/*
OBJ Class ;
  char *cp ;

  if (!isWriteAccess())
    { messout("Sorrry, you do not have write access") ;
      return ;
    }
 
  if (iskey(classChosen) != 2)
    { messout ("This class no longer exists, sorry") ;
      return ;
    }
*/

  /* I save once here to let the system a chance to crash if need be */
 
}

/****************************************************************/

static MENUOPT
   classMenu[] = {
     graphDestroy,"Quit",
     help,"Help",
     0,0
     } ;

/****************************************************************/

static void classDrawBranch(CT* ct, Array x)
{ /*
int i, g = ct->generation , n = arrayMax(classTree) ;
  CT *ct1 ;

  if(!array(x, g, int))
    array(x, g, int ) = 4 ;
*/
         /* Draw self */
/*
  ct->box = graphBoxStart() ;
  ct->x = array(x, g, int) ;
  ct->y = 7 + 3*g ;
  graphText(name(ct->key), ct->x, ct->y ) ;
  ct->len =strlen(name(ct->key)) ;
  array(x, g, int) += ct->len + 5 ;

  graphBoxEnd() ;
  if(TRUE)
    ct->color = WHITE ;
  else
    ct->color = YELLOW ;
  if (ct->key == _This_class || ct->key == KEYMAKE(_VClass, 1))
    ct->color = LIGHTBLUE ;

  graphBoxDraw(ct->box, BLACK, ct->color) ;
  for (i = 0, ct1 = arrp(classTree,0,CT) ; i < n ;
       ct1++, i++)
    if (ct->ancester == ct1->key)
      graphLine(ct1->x + ct1->len/2, ct1->y + 1.0, ct->x + ct->len/2 , ct->y - .0) ;
*/  
               /* Recursion */
/*
  for (i = 0, ct1 = arrp(classTree,n - 1,CT) ; i < n ;
       ct1--, i++)
     if (ct1->cta == ct)
       classDrawBranch(ct1, x) ;
*/
}

/****************************************************************/

static void classDraw() 
{ int i, max = 0 , maxGen = 0 , n ;
  CT *ct ;
  Array maxX = arrayCreate(12, int) ;
  
/*
  if(!classTree)
    classTreeConstruct() ;
*/
  
  n = arrayMax(classTree) ;
  graphClear() ;
  graphPop() ;
  for (i = 0, ct = arrp(classTree,n - 1,CT) ; i < n ; ct--, i++)
    { if (ct->generation == 1)
	classDrawBranch(ct, maxX) ;
      if (ct->generation > maxGen)
	maxGen = ct->generation ;
    }

  for (i=0 ; i< arrayMax(maxX); i++)
    if (max < arr(maxX,i, int))
      max = arr(maxX,i,int) ;
  if (!max)
    max = 40 ;
 
  graphTextBounds (max + 3 , 7 + 3 * maxGen) ;
  graphButtons(classMenu, 2., 2., 40) ;
  graphRedraw() ;
  arrayDestroy(maxX) ;
}

/*************************************************************/

void classControl(void)
{ 
 
  if ( graphActivate(classGraph))
    graphPop() ;
  else
    classGraph =  displayCreate(DtClass) ;
  graphMenu(classMenu) ;

  graphRegister(DESTROY,localDestroy) ;
  graphRegister(RESIZE, classDraw) ;
  graphRegister(PICK,  classPick) ;
  classDraw() ;
}

#endif /* CLASSCONTROL */


/************ end of file ****************/
