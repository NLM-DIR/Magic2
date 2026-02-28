/*  file: mainpick.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: May 13 14:11 1999 (fw)
 * * Dec  3 14:42 1998 (edgrif): Change calls to new interface to aceversion.
 * * Oct 22 11:45 1998 (edgrif): Remove dec. of pickDraw, put in acedb.h
 * * Jun 11 1:38 1996 (rbrusk):
 *		-	WIN32 letting double pickPick() do its thing in pickCreate()
 * * Jun  6 17:44 1996 (rd)
 * * Jun 5 1:38 1996 (rbrusk):
 *	-	TEXT_ENTRY_FUNC_CAST added to classDisplay() and grepDisplay()
 *		callback arguments to graphTextScrollEntry()
 * * Jun 5 1:38 1996 (rbrusk):
 *	-	WIN32 port changes to 4.3
 * * May 10 21:13 1996 (rd)
 * * Sep 12 17:45 1992 (rd): in classDisplay() added quick check if
	template is a valid object name (Alan's suggestion)
 * * Dec  6 13:28 1991 (mieg): autocompletion system
 * * Dec  2 20:59 1991 (mieg): if(isdisplayBlock) no autoPick
 because this create a textBound on wrong graphtype somehow ?
 * * Nov  5 18:28 1991 (mieg): call to keySetShow changed
 * Created: Tue Nov  5 18:26:56 1991 (mieg)
 *-------------------------------------------------------------------
 */


/* $Id: mainpick.c,v 1.14 2017/07/23 15:03:18 mieg Exp $ */

/************************************************************/

#include "display.h"
#include "main.h"
#include "bitset.h"
#include "graph.h"
#include "key.h"
#include "keyset.h"
#include "lex.h"
#include "pick.h"
#include "sysclass.h"
#include "../wh/menu.h"
#include "client.h"
#include "pref.h"
#include "freeout.h"
#include "help.h"
#include "chrono.h"
#include "spreaddisp.h"
#include "model.h"

#include "query.h"
#include "querydisp.h"
#include "bqldisp.h"

#include "session.h"
#include "sessiondisp.h"

#include "spread.h"

#include "longtext.h"

#include <ctype.h>		/* for isprint */

/************************************************************/

static BOOL isNewLayout = TRUE ;    /* mieg:  new layout */

static int PICK_MAGIC = 49662 ;    /* also use address as associator handle */

typedef struct PickStruct
  { int    magic ;  /* == PICK_MAGIC */
    int    curr, currBox ;
    KEY    currClasse ;
    char   template[256] ;
    int    templateBox ;
    char*  grepText ;
    int    grepBox , longGrepBox , m1, m2, m3, modeBox1, modeBox2 ;
    BOOL   longGrep ;
    char   whatdoIdoText[256] ;
    int    whatdoIdoBox, classListBox ;
    Graph  mainGraph ;
    Graph  dispGraph ;
    KEYSET oblist ;
    void*  dispHandle ;
  } *Pick ;
#define PICKSIZE sizeof(struct PickStruct)

#define PICKGET(myname)     Pick pick ;\
                                       \
                          if (!graphAssFind (&PICK_MAGIC, &pick)) \
		            messcrash ("%s can't find pointer",myname) ; \
                          if (pick->magic != PICK_MAGIC) \
                            messcrash ("%s received a wrong pointer",myname)

typedef struct {
  KEY key ;			/* in _VClass */
  int x, y ;
} LAYOUT ;

/************************************************************/

Array classLayout = 0 ;		/* NB Global */

/************************************************************/

static int   firstClassBox = 0 ;
static Pick  mainPick = 0 ;
static Array classList = 0 ;
static MENU  mainMenu ;	

static VoidRoutine xacedbQuitFunction = 0;

static void pickDestroy (void) ;
static int pickClassComplete (char *cp, int len) ;
static BOOL classListMenuCall (KEY k, int box) ;
static void pickDoDrawNew(void) ;
static void getLayout (void);

/*************************************************************/

void messStatusDisplayInWindow(const char * text)
/* registered as the messStatusRoutine using
   messStatusRegister() in acedbGraphInit() */
{
#ifndef MACINTOSH
  /* non-MACINTOSH */
  Pick pick = mainPick ;
  static const char *old = "dummy" ; /* some impossible value */
  Graph oldGraph ;

  if (old == text || !pick) 
    return ;

  oldGraph = graphActive() ;
  graphActivate(pick->mainGraph) ;
  if (!text)
    { old = 0 ;
      strcpy(pick->whatdoIdoText, "Ready") ;
      graphBoxDraw(pick->whatdoIdoBox,BLACK,PALEGREEN) ;
    }
  else if (text != old)
    { old = text ;
      strncpy(pick->whatdoIdoText, text, 250) ;
      if (pick->whatdoIdoBox)
	graphBoxDraw(pick->whatdoIdoBox,BLACK,RED) ;
    }
  
  graphActivate(oldGraph) ;
#else
  /* MACINTOSH */
  extern void StartAsyncSpinning (short period) ;
  extern void StopAsyncSpinning (void) ;

  static char *old = "dummy" ; /* some impossible value */

  if (old != text )
  {
    if (!text)
      { old = 0 ;
        StopAsyncSpinning () ;
      }
    else if (text != old)
      { old = text ;
        StartAsyncSpinning (8) ;
      }
  }
#endif /* MACINTOSH */
} /* messStatusDisplayInWindow */

/****************************************************************/

void quitFunctionRegister (VoidRoutine func)
{
  if (func)
    xacedbQuitFunction = func;
} /* quitFunctionRegister */

/****************************************************************/
/* state = 1: template unchanged
           2,  added a star at end
           3, autocomplete
	   */

static void obDisplay (Pick pick, BOOL show1, int state)
{
	/* then switch to the display graph (make if necessary) */
  char *message;

  if (pick->currClasse)
    message = messprintf("%s:%s",name(pick->currClasse),pick->template) ;
  else
    message = messprintf("All classes:%s",pick->template) ;

  switch (state)
    {
    case 1: 
      if (!keySetMax(pick->oblist))
	{
	  if (pick->currClasse)
	    message = messprintf("No %s:%s %s",
				 name(pick->currClasse),
				 pick->template, 
				 "is the class correct ? ") ;
	  else
	    message = messprintf("No object named %s in any class", 
				 pick->template) ;
	}
      break ;

    case 2: 
      if (keySetMax(pick->oblist))
	{
	  if (pick->currClasse)
	    message =  messprintf("%s:%s %s",
				  name(pick->currClasse),
				  pick->template, 
				  "No exact match, * appended") ;
	  else
	    message =  messprintf("%s %s",
				  pick->template, 
				  "No exact match in any class, * appended") ;
	}
      else
	{
	  if (pick->currClasse)
	    message = messprintf("No %s:%s %s",
				 name(pick->currClasse),
				 pick->template, 
				 " is the class correct ? ") ;
	  else
	    message = messprintf("No object named %s in any class",
				 pick->template) ;
	}
      break ;

    case 3:
      if (keySetMax(pick->oblist))
	message = "No possible completion in this class" ;
      break;

        case 4:  /* mhmp 11.07.02 */
			message = "" ;
			messout ("To search text:\n Please specify a search string\n with at least 3 letters or digits.\n Authorized wild chars are\n * (any string), ? (single char).") ;
      break;

    case 5: 
      if (!keySetMax(pick->oblist))
	{ state = 4 ; message = messprintf("%s Not found",pick->grepText) ; }
      else
	message = messprintf("Search %s",pick->grepText) ;
      break ;
    case 6: 
      if (!keySetMax(pick->oblist))
	{
	  state =  4 ;
	  message = messprintf("%s Not found",pick->grepText) ;
	}
      else
	message = messprintf("Long Search %s",pick->grepText) ;
      break ;

    case 7: 
      if (!keySetMax(pick->oblist))
	{ 
	  state = 4 ; 
	  message = messprintf("Text %s Not found in any class",
			       pick->grepText) ; }
      else
	message = messprintf("Search %s in all classes",pick->grepText) ;
      break ;

    case 8: 
      if (!keySetMax(pick->oblist))
	{ 
	  state =  4 ; 
	  message = messprintf("%s Not found in any class",
			       pick->grepText) ;
	}
      else
	message = messprintf("Long Search %s in all classes",
			     pick->grepText) ;
      break ;

    case 9: 
      message = messprintf("%s not found hit in selected class, "
			   "searched all classes (no hit)",
			   pick->grepText) ;
      break ;

    case 10: 
      message = messprintf("%s not found in selected class, "
			   "long-searched all classes (no hit)",
			   pick->grepText) ;
      break ;
    }

  switch (state)
    {
    case 5: 
    case 6: 
    case 7: 
    case 8: /* costly search, reply in new window */
      displayCreate(DtKeySet) ;
      graphRetitle(message) ;
      keySetMessageShow (pick->oblist,0,message) ;
      break ;
    default:       /* fast search, reply in main window */
      if (graphExists (pick->dispGraph))
	 keySetMessageShow (pick->oblist,pick->dispHandle,message) ;
      else 
	{ 
	  pick->dispGraph = displayCreate(DtKeySet) ;
	  graphRetitle("Main KeySet") ;
	  pick->dispHandle = keySetMessageShow (pick->oblist,0,message) ;
	}
      break ;
    }

  if (keySetMax(pick->oblist) == 1 && show1)
    { 
      KEY key = keySet(pick->oblist,0) ;
      
      display(key, 0, keySetCurrent (pick->dispHandle, key));
    }

  pick->oblist = 0 ;		/* since given to keySetShow */

#if !defined(WIN32)  /* The WIN32 user interface is nicer
			if I don't do this */
   graphActivate (pick->mainGraph) ;
#endif
}

/**************** text action **********************/

static BOOL classDisplay ()  /* displays current option list */
{ 
  char *cp ;
  PICKGET("classDisplay") ;
  
  if (!pick->curr)
    { messout("Please pick a class") ;
      return FALSE ;
    }

  cp = pick->template ;
  while (*cp == ' ') ++cp ;
  if (!*cp)
    { strcpy (pick->template, "*") ;
      graphBoxDraw (pick->templateBox, -1, -1) ;
    }

  graphCompletionEntry (pickClassComplete,pick->template,0,0,0,0) ;

  return classListMenuCall
    (arrp(classLayout,pick->curr,LAYOUT)->key, 0) ;
}

/**************** text action **********************/

static BOOL classDisplayNew (int curr)  /* displays current option list */
{ 
  char *cp ;
  PICKGET("classDisplay") ;

  cp = pick->template ;
  while (*cp == ' ') ++cp ;
  if (!*cp)
    { strcpy (pick->template, "*") ;
      graphBoxDraw (pick->templateBox, -1, -1) ;
    }

  graphCompletionEntry (pickClassComplete,pick->template,0,0,0,0) ;

  return classListMenuCall
    (curr ? arrp(classLayout,curr,LAYOUT)->key : 0, 0) ;
}

/*************************************************/

static void grepDisplayNew (int curr)
{ KEYSET ob1 = 0, ob2 = 0, ob3 = 0;
  int n , state ;
  static int confirm = 0 ;
  char *cp,*cq ;
  static char localText[256] , lastText[256] ;
  PICKGET("grepDisplay") ;
  
/* first check that the search text is non-trivial */

  for (n = 0, cp = pick->grepText ; *cp ; ++cp)
    if (isprint((int)*cp) && *cp != '*' && *cp != '?')
      ++n ;
  if (n < 3)
    {
      /*
	messout ("To search text please specify a search string " 
	       "with at least 3 letters or digits."
	       "   Authorized wild chars are\n"
	       "* (any string), ? (single char).") ;
	       */
      
      pick->longGrep = FALSE ;
      graphBoxDraw(pick->longGrepBox, BLACK, WHITE) ;
      pick->oblist = keySetCreate () ;
      obDisplay (pick, FALSE, 4) ;
      return ;
    }

/* make localText with a '*' in front of and behind the grep text */
  
  cp = pick->grepText ;
  if (*cp == '*')
    cq = localText ;
  else
    { *localText = '*' ;
      cq = localText+1 ;
    }
  while ((*cq++ = *cp++)) ;
  if (cq[-2] != '*')
    { cq[-1] = '*' ; *cq = 0 ; }

/* If the text has not changed, wait for confirmation */
  if(!strcmp(localText,lastText))
    if(!confirm++)
      { graphBoxDraw(pick->longGrepBox, BLACK, WHITE) ;
	return ;
      }
  strncpy(lastText,localText,255) ;
  strncpy(pick->grepText,localText,255) ;
  confirm = 0 ;

/* then do the search */
  graphBoxDraw(pick->grepBox,BLACK,RED) ;
  if (pick->longGrep)
    graphBoxDraw(pick->longGrepBox, BLACK, YELLOW) ;
  if (externalServer)
    ob1 = externalServer (0, 0, localText, FALSE) ;
  else 
    ob1 = queryGrep(0, localText) ;
  if (pick->longGrep)
    { graphBoxDraw(pick->grepBox,BLACK,YELLOW) ;
      graphBoxDraw(pick->longGrepBox, BLACK, RED) ;
      ob2 = longGrep(localText) ;
      ob3 = keySetOR(ob1, ob2) ;
    }
  else
    { ob3 = ob1 ; ob1 = 0 ; }

  keySetDestroy(ob1) ;
  keySetDestroy(ob2) ; 

  graphBoxDraw(pick->grepBox,BLACK,YELLOW) ;
  state = pick->longGrep ? 6 : 5 ;
  pick->longGrep = FALSE ;
  graphBoxDraw(pick->longGrepBox, BLACK, WHITE) ;

  if (keySetMax(ob3) && curr) /* filter by single class */
    {
      int n2 = keySetMax(ob3), n3 = 0 ;
      KEY *kp2, *kp3 ;
      KEY cc = 	arrp(classLayout, curr, LAYOUT)->key ;
     
      /* it is better to filter locally than to query, becausee of client/server */
      ob2 = arrayCopy (ob3) ;
      kp2 = arrp(ob2, 0, KEY) ; kp3 = arrp(ob3, 0, KEY) ;
      while (n2--)
	{
	  if (lexIsInClass(*kp2, cc))
	    { *kp3++ = *kp2 ; n3++ ;}
	  kp2++ ;
	}
      keySetMax(ob3) = n3 ;
      if (n3)
	keySetDestroy (ob2) ;
      else
	{
	  state += 4 ; /* report enlarged search */
	  keySetDestroy (ob3) ;
	  ob3 = ob2 ; ob2 = 0 ;
	}
    }
  else
    state += 2 ; /* report general failure */
  pick->oblist = ob3 ;
  obDisplay (pick, FALSE, state) ;
}

static void grepDisplay (void)
{  grepDisplayNew(0) ; }

/**************** pick Action ************/

static int pickClassComplete (char *cp, int len)
{ 
  PICKGET("pickClassComplete") ;

  if (!pick->curr)
    return 0 ;

  return ksetClassComplete (cp, len, 
		arrp(classLayout, pick->curr, LAYOUT)->key) ;
}

static void pickModeChoice1 (void) ;
static void pickPick (int k)
{ 
  PICKGET("pickPick") ;
  
  if (k == pick->templateBox)
    { graphCompletionEntry (pickClassComplete, pick->template,0,0,0,0) ; 
      return ;
    }

  if (!k) 
    {
      if (pick->currBox)
	graphBoxDraw (pick->currBox, BLACK, WHITE) ;
      pick->curr = 0 ; pick->currBox = 0 ;
      return ;
    }

  if (k == pick->grepBox)
    {
      if (pick->currBox)
	graphBoxDraw (pick->currBox, BLACK, WHITE) ;
      pick->curr = 0 ; pick->currBox = 0 ;
      graphTextScrollEntry (pick->grepText,0,0,0,0,0) ;
    }
  else if (k == pick->classListBox)
    { KEY key ;
      if (graphSelect (&key, arrp (classList, 0, FREEOPT)))
	classListMenuCall (key, 0) ;
    }

  if (k > firstClassBox && k < firstClassBox + arrayMax(classLayout))
    { 
      
      pickModeChoice1() ;
      if (k != pick->currBox)
	{ if (pick->currBox)
	    graphBoxDraw (pick->currBox, BLACK, WHITE) ;
	  pick->currBox = k ;
	  if (isNewLayout)graphBoxDraw (pick->currBox, BLACK,PALEBLUE) ;
	  else graphBoxDraw (pick->currBox, WHITE, BLACK) ;
	  pick->curr = k - firstClassBox ;
	  graphCompletionEntry (pickClassComplete,pick->template,0,0,0,0) ;
	}
      else if (!classDisplay () && strcmp (pick->template,"*"))
	{ 		/* inappropriate template - display whole list */
	  strcpy (pick->template,"*") ;
	  graphCompletionEntry (pickClassComplete,pick->template,0,0,0,0) ;
	  classDisplay () ;
	}
    }
}

/****************** destroy routine ************/

static void pickDestroy (void)
{ 
  Graph g = graphActive() ;
  PICKGET("obDestroy") ;

  if (graphActivate(pick->dispGraph))
      graphDestroy () ;

  if (graphActivate(pick->mainGraph))
     graphDestroy () ;

  graphActivate (g) ;
  graphAssRemove (&PICK_MAGIC) ;

  pick->magic = 0 ;
  if (pick->grepText) 
    messfree (pick->grepText);

  messfree (pick) ;
  mainPick = 0 ;
}

/***************** public routine *****************************/

static Graph mainGraph = 0 ;
void pickPopMain(void)
{
  if(graphActivate(mainGraph))
    graphPop() ;
}

void pickReReadClasses (void)	/* called when models change */
{
  arrayDestroy (classLayout) ;
  classLayout = 0;		/* force a re-create */

  getLayout();
} /* pickReReadClasses */

static void getLayout (void)
{ 
  FILE *fil ;
  LAYOUT *z ;
  int ii, n = 0, level, x, y = 0 ;
  KEY key ;
  char *card, *word ;

  if (classLayout)
    return ;

  classListCreate (&classList) ;

  classLayout = arrayCreate (32, LAYOUT) ;

  if (filName ("wspec/layout","wrm","r") &&
      (fil = filopen ("wspec/layout","wrm","r")))
    { level = freesetfile (fil, 0) ;
      freespecial ("\n\t\"") ;
      while ((card = freecard (level)))
	{ while (TRUE)
	    { x = freepos() - card ;
	      if (!(word = freeword()))	/* exit loop */
		break ;
	      if (lexword2key (word, &key, _VClass))
		{ z = arrayp(classLayout, ++n, LAYOUT) ;
		  z->key = key ;
		  z->x = x ;
		  z->y = y ;
		  ii = pickWord2Class (word) ;
		  if (ii > 0) pickList[ii & 255].visible = 'V' ;
		}
/*
	      else
		messout ("%s in wspec/layout.wrm is not a class name", word) ;
*/
	    }
	  ++y ;
	}
      /* filclose (fil) ; already closed, crashes linux, srk */
    }
  else
    { 
      int x, i, ncol ;		/* old code using pickVocList */
      int graphWidth, graphHeight ;
      int maxCol[10] ;
      graphFitBounds (&graphWidth,&graphHeight) ;

      for (i = 1 ; i <= pickVocList[0].key ; ++i)
	if (lexword2key (pickVocList[i].text, &key, _VClass))
	  { z = arrayp(classLayout, ++n, LAYOUT) ;
	    z->key = key ;
	    z->x = strlen(name(key)) + 2 ;
	  }
	else
	  messout ("pickVocList[%d] = %s is not a class name",
		   i, pickVocList[i].text) ;

      for (ncol = 10 ; ncol > 1 ; --ncol)
	{ for (i = 0 ; i < ncol ; ++i)
	    maxCol[i] = 0 ;
	  x = 0 ;
	  for (i = 1 ; i < arrayMax(classLayout) ; ++i)
	    { z = arrp(classLayout, i, LAYOUT) ;
	      if (z->x > maxCol[(i-1)%ncol])
		{ x += z->x - maxCol[(i-1)%ncol] ;
		  if (x > graphWidth)
		    break ;
		  maxCol[(i-1)%ncol] = z->x ;
		}
	    }
	  if (i == arrayMax(classLayout))
	    break ;
	}
      /* use it */
      x = 2 ;
      for (i = 1 ; i < arrayMax(classLayout) ; ++i)
	{ z = arrp(classLayout, i, LAYOUT) ;
	z->x = x ;
	z->y = y ;
	if (i % ncol) x += maxCol[(i-1)%ncol] ;
	else { x = 2 ; ++y ; }
	}

    }
}

static void resizeLayout (void)
{   
  isNewLayout = !prefValue ("OLD_STYLE_MAIN_WINDOW") ;
  if (!isNewLayout) 
    {
      arrayDestroy (classLayout) ;
      classLayout = 0;		/* force a re-create */
    }
  pickDraw() ;
}

/***************** main create routine ************************/

static void pickLongGrep (void)
{    
  PICKGET("pickLong") ;
  
  pick->longGrep = TRUE ;
  grepDisplay() ;   
}

void classListCreate (Array *classListP)
{ int max ;
  KEYSET ks, ksa ;
  
  ks = queryLocalParametrized
    (0, "FIND Class HIDDEN OR VISIBLE", ""); /* not burried */

  max = keySetMax (ks) ;
  ksa = keySetAlphaHeap (ks, max) ;
  *classListP = arrayReCreate (*classListP, max, FREEOPT) ;
  array (*classListP, 0, FREEOPT).key = max ;
  array (*classListP, 0, FREEOPT).text = "" ;
  while (max--)
    { array (*classListP, max + 1, FREEOPT).key = keySet(ksa, max) ;
      array (*classListP, max + 1, FREEOPT).text = name(keySet(ksa, max)) ;
    }
  keySetDestroy (ks) ;
  keySetDestroy (ksa) ;
  
  return;
} /* classListCreate */

static BOOL classListMenuCall (KEY classe, int box)
{ char *cp ; KEY key ; KEYSET ob ;
  BOOL resul ;
  char *nom ;
  int n, state = 1, ii , i, j ;
  PICKGET ("classListMenu") ;
  
  pick->currClasse = classe ;
  pick->oblist = keySetCreate() ; j = 0 ;
lao:
  /* classe is a real Classe-key in _VClass */
  for (ii = 0 ; ii < 256 ; ii++)
    {
      if (!classe && (pickType(ii) == 'X' || pickList[ii].visible == 'h' || !pickList[ii].name))
	continue ;
      if (classe) 
	{ ii = classe ; nom = name(classe) ; }
      else
	nom = className(ii << 24) ;
      
      if (!strcmp("*", pick->template))
	{
	  if (classe && lexword2key (pick->template, &key, classe))
	    { 
	      ob = keySetCreate () ;
	      keySet(ob, 0) = key ;
	    }
	  else
	    { cp = messprintf( "FIND %s", nom) ;
	    if (externalServer) 
	      { 
		ob = externalServer (200, cp, 0, FALSE) ;
	      if (keySetMax(pick->oblist) == 200)
		messout("I limited the importation to 200, type ** to get all %s(s)", name(classe)) ;
	      }
	    else
	      ob = query(0, cp) ;
	    }
	}
      else
        ob = query(0, messprintf( "FIND %s IS %s",nom,
				  freeprotect(pick->template))) ;
      n = keySetMax(ob) ;
      for (i = 0 ; i < n ; i++)
	keySet(pick->oblist,j++) = keySet(ob,i) ;
      if (externalServer && !strcmp("*", pick->template) && j >= 200)
	{ 
	  messout("I limited the importation to 200, type ** to get all %s(s)", nom) ;
	  break ;
	}
      if (classe)
	break ;
    }
 
  if (state == 1 && !keySetMax(pick->oblist))
    { cp = pick->template ; n = strlen (cp) ;
      if (!n || (*(cp + n - 1) != '*' && n < 19))
	{ strcat (pick->template, "*") ;
	  graphBoxDraw (pick->templateBox, -1, -1) ;
	  state = 2 ;
	  goto lao ;
	}
    }
	/* mhmp 23.06.02 */ 
	arraySort (pick->oblist, keySetOrder) ;
	arrayCompress(pick->oblist);
  resul = (keySetMax(pick->oblist) != 0) ;
  obDisplay (pick, TRUE, state) ;
  return resul ;
}



/***********************************/

int getCurrentClass (void)
{
  KEY class = 0;
  unsigned char mask ;
  
  if (mainPick->curr)
    {
      class = pickWord2Class (pickVocList[mainPick->curr].text);
      if (!class)
	{
	  /* deal with subclasses */
	  lexword2key (pickVocList[mainPick->curr].text, &class, _VClass);
	  pickIsA (&class, &mask);
	}
    }
  return class;
}

/***********************************/


void helpOnClass (void)
{
  PICKGET ("aboutClass") ;
  
  if (!pick->curr)
    messout ("Please pick a class") ;
  else
    helpOn (name(arrp(classLayout,pick->curr,LAYOUT)->key)) ;
} /* helpOnClass */


static MENUOPT helpButtonMenu[] = 
{
  {help, "General help"},
  {helpOnClass, "Help on active class"},
  {0, 0}
} ;

#ifdef JUNK
extern void  dnaAnalyse(void),
	     queryCreate(void), 
             qbuildCreate(void), 
             qbeCreate(void), 
             spreadCreate(void), 
             gMapCompute(void) ,
	     blyControl(void),
             parseControl(void),
             dumpAll(void) ,
#ifdef DEBUG
	     acedbtest(void),
#endif /* DEBUG acedbtest */
             fp2 (void),
	     acedbstatus(void),
             chronoShow(void),

             readModels(void),
	     dumpAsnDefinitions(void),
#ifdef WCS
		wcsShow(void),
		wcsAnnotate(void),
	     	wcsSearch(void),
#endif
	     alignMaps(void), 
             updateData(void),
	     addKeys(void),
/*	     metaCheckCode(void),	*/
/* 	     setUpdate (void),  detlef's test code */
             invokeDebugger(void),
  /*             editPreferences(void), */
             help(void);
#else

#ifdef ACEMBLY
 extern void blyControl (void);
#endif /* ACEMBLY */

#endif /* JUNK */


static MENUOPT pickNewMenu[] = 
{
  { help,	"Help"},
  { graphPrint,	"Print"},
#ifdef ACEMBLY
  { menuSpacer,	""},
  { blyControl,	"Acembly"},
#endif
  { graphDestroy,"Exit"},
  {0,0}
} ;

#ifndef ACEMBLY
extern void oxgridCreate(void) ;
extern BOOL oxgridPossible(void) ;
static MENUOPT pickNewOxMenu[] = 
{
  {graphDestroy,"Exit"},
  {help,"Help"},
  {graphPrint,"Print"},
  {menuSpacer,""},
  { oxgridCreate, "Ox Grids"},
  {0,0}
} ;
#endif

static void localCleanUp (void)
{
  graphCleanUp() ;
  /* mainPick->curr = 0; */
  pickDoDrawNew() ;
}

/*************************************/

static void queryButtonFreeAction(KEY key, int box)
{
  switch (key)
    {
    case 0: break ;
    case 'h': helpOn ("Query") ; break ;
    case 'q': queryCreate() ; break ;
    case 'B': bqldCreate() ; break ;
    case 'b': qbeCreate() ; break ;
    case 'u': qbuildCreate() ; break ;
    case 's': spreadDispCreate(FALSE) ; break ;
    case 'd': dnaAnalyse() ; break ;
    }
}

static FREEOPT queryButtonFreeMenu[] = 
{
  {7, "Query menu"},
  {'h', "Help on Query"},
  {0, ""},
  { 'q',      "Query"},   
  { 'B', "BQL Table Query"},
#if !defined(MACINTOSH)
  { 'b',        "Query by Examples"},  
#else  /* always keep same number of lines */
  {0, ""},
#endif
  { 'u',     "Query builder"},   
  { 's',     "Table Maker"},
  { 'd',       "DNA Analysis"},
} ;

static void querySelector(void)
{
  KEY k ;
  
  if (!graphSelect (&k,queryButtonFreeMenu))
    return ;
  queryButtonFreeAction(k,0) ;
}

/***********************************/

static void adminButtonFreeAction(KEY key, int box)
{
  switch (key)
    {
    case 0: break ;
    case 's': acedbstatus(); break ;
    case 'S': sessionControl() ; break ;
    case 'p': editPreferences() ; break ;
    case 'P': paletteDisplay() ; break ;
    case 'd': dumpAll() ; break ;
    case 'c': chronoShow () ; break ;
#ifdef DEBUG
    case 't': acedbtest() ; break ;
#endif /* DEBUG acedbtest */
    }
}

static FREEOPT adminButtonFreeMenu[] = 
{
#if (defined DEBUG)
  {9, "Admin menu"},
#else
  {7, "Admin menu"},
#endif /* DEBUG */
  {'s', "Program Status"},
  /*   {'S', "Session Control"}, suppressed 2016_09_01, never used and i do not want to supprt it */
  {'p', "Preference"},
  {'P', "Palette"},
  { 0,  ""},
  {'d', "Dump All"},

  {0,""},
  {'c', "Chronometer"},

#ifdef DEBUG
  {0,""}
  {'t', "Test subroutine"},
#endif /* DEBUG acedbtest */
} ;


static void adminSelector(void)
{
  KEY k ;
  
  if (!graphSelect (&k,adminButtonFreeMenu))
    return ;
  adminButtonFreeAction(k,0) ;
} /* adminSelector */

/************************************************************/
/* interactive way to check, whether we have write access.
   Used in by functions that need write access, but want to
   enable the user to grab it, if we don't have it already,
   and otherwise (when it returns FALSE) return from that func */
/************************************************************/
BOOL checkWriteAccess(void)
{ 
  if (isWriteAccess())
    return TRUE;

  if (writeAccessPossible())
    { 
      if (messQuery("You do not have Write Access, "
		    "shall I get it for you?"))
	{
	  /* try to get it */
	  if (sessionGainWriteAccess())
	    /* it worked */
	    return TRUE;
	  else
	    /* we couldn't get it */
	    return FALSE;
	}
      else
	/* didn't want to gain write access */
	return FALSE;		
    }

  
  /* impossible to get write access */
  messout("Sorry, you cannot gain Write Access");
  return FALSE;
} /* checkWriteAccess */

/************************************************************/

static void editButtonFreeAction(KEY key, int box)
{
  switch (key)
    {
    case 0: break ;
    case 'r':			/* Read .ace file */
      if (checkWriteAccess())
	parseControl(); 
      break ;

    case 'a':			/* Add/Alias/Rename */
      if (checkWriteAccess())
	addKeys() ;
      break ;

    case 'A':			/* Align maps */
      if (checkWriteAccess())
	alignMaps() ;
      break ;

    case 'D':
      sessionAutoSave (0, 0) ;   /* drop auto save */
      break;

    case 'd':
      sessionAutoSave (0, 0) ;
      sessionReleaseWriteAccess(); /* don't save, drop write access */
      break;

    case 'W': 
      sessionAutoSave (60 * 30, 60 * 30) ;  /*  (60 * 30, 60 * 5)  */

    case 'w': 
      sessionGainWriteAccess();	/* Gain Write Access */
      /* will re-draw if write access status changes */
      break;

    case 'm':
      readModels() ;
      if (mainPick) mainPick->curr = 0 ; /* layout may change */
      break ;
    case 'u':
      updateData() ;
      break ;
    }
  pickDoDrawNew() ;
} /* editButtonFreeAction */

/* this menu is only displayed if write access is possible */

static FREEOPT editButtonFreeMenu1[] = 
{
  {9, "Edit menu"},
  {'u', "Add Update"},
  {'A', "Align maps"},
  {0,""},
  {'r', "Read .ace file"},
  {'a', "Add/Edit/Rename"},
  {0,""},
  {'m', "Read Models"},
  {'W', "Gain write access + Auto save"},
  {'w', "Gain write access"} 

} ;

static FREEOPT editButtonFreeMenu2[] = 
{
  {9, "Edit menu"},
  {'u', "Add Update"},
  {'A', "Align maps"},
  {0,""},
  {'r', "Read .ace file"},
  {'a', "Add/Edit/Rename"},
  {0,""},
  {'m', "Read Models"},
  {'W', "Gain Auto save"},
  {'d', "Drop write access (no save)"}
} ;

static FREEOPT editButtonFreeMenu3[] = 
{
  {9, "Edit menu"},
  {'u', "Add Update"},
  {'A', "Align maps"},
  {0,""},
  {'r', "Read .ace file"},
  {'a', "Add/Edit/Rename"},
  {0,""},
  {'m', "Read Models"},
  {'D', "Drop auto save (you will need to save manually)"},
  {'d', "Drop write access (no save)"}

} ;
static void editSelector1(void)
{
  KEY k ;
  
  /*   graphEvent (RIGHT_DOWN, 3,3) ; */
  if (!graphSelect (&k,editButtonFreeMenu1))
    return ;
  editButtonFreeAction(k,0) ;
}

static void editSelector2(void)
{
  KEY k ;
  
  /*   graphEvent (RIGHT_DOWN, 3,3) ; */
  if (!graphSelect (&k,editButtonFreeMenu2))
    return ;
  editButtonFreeAction(k,0) ;
}

static void editSelector3(void)
{
  KEY k ;
  
  /*   graphEvent (RIGHT_DOWN, 3,3) ; */
  if (!graphSelect (&k,editButtonFreeMenu3))
    return ;
  editButtonFreeAction(k,0) ;
}

/***********************************/

static void saveAndKeepWriteAccess(void)
     /* used by save button freemenu and the save button itself */
{
  if (isWriteAccess())
    {
      messStatus ("Saving") ;
      sessionClose (TRUE);	/* close and save session,
				   nothing will happen, 
				   if we don't have write access */

      if (!sessionGainWriteAccess()) /* re-gain write access */
	messout ("Sorry, it was not possible to re-gain write access "
		 "after saving the last session.");
    }
  else
    {
      messout ("You don't have write access. Nothing to save.");
    }

} /* saveAndKeepWriteAccess */

/***********************************/

static void saveAndLoseWriteAccess(void)
     /* used by save button freemenu */
{
  messStatus ("Saving") ;
  sessionClose (TRUE) ;		/* close and save session,
				   nothing will happen, 
				   if we don't have write access */
  sessionAutoSave (0, 0) ;
} /* saveAndLoseWriteAccess */

/***********************************/

static MENUOPT save[] = {
	{saveAndKeepWriteAccess, "Save session and keep write access"},
	{saveAndLoseWriteAccess, "Save session and lose write access"},
        {0,0}
} ;

/***********************************/

void exitButtonAction (void)
     /* non-static, because extern'd by quovadis.c for now */
{
  if (thisSession.session != 1 &&
      !messQuery ("Do you really want to quit acedb?"))
    return;

  if (!xacedbQuitFunction)
    (*xacedbQuitFunction)();
  else
    {
      graphCleanUp () ;  /* kills all open windows,
			    forces displayed objects to cache */
      
      writeAccessChangeRegister (0); /* avoid unnecessary re-draws
					just before exiting */

      /* release read/write locks and clean up temp-files 
	 is all done by aceQuit() */
      if (isWriteAccess() && 
	  messQuery ("You did not save your work, should I ?"))
	aceQuit (TRUE);
      else
	aceQuit (FALSE);
      
      if (!getenv("ACEDB_NO_BANNER"))
	printf ("\n\nA bientot !\n") ;
    }
  exit (EXIT_SUCCESS);
} /* exitButtonAction */

/***********************************/

static void serverOn(void)
{
  BOOL bb=FALSE;
  if (externalServerState)
    bb = externalServerState (1) ; /* set */
  mainPickReport (bb ? 1 : 2) ;
}

static void serverOff(void)
{
  BOOL bb=FALSE;
  if (externalServerState)
    bb = externalServerState (0) ; /* unset */
  mainPickReport (bb ? 1 : 0) ;
}

void mainPickReport (int n)
{
  BOOL bb = FALSE ;
  float line = .5 ;
  static int boxOn = 0, boxOff ; 
  char *cp = 0 ;
  if (!externalServerState)
    return ;

  switch (n)
    {
    case -1: /* ask and draw */
      bb = externalServerState (-1) ; /* ask */
      graphText ("Server:", .5, line) ;
      boxOn = graphButton ("On", serverOn, 8.5, line) ;
      boxOff = graphButton ("Off", serverOff, 12.0, line) ;
      break ;
    case 1:
      bb = TRUE ;
      break ;
    case 0: case 2:
      bb = FALSE ;
      break ;
    default:
      break ;
    }
  if (!boxOn) return ;

  if (externalServerName) 
    cp = externalServerName () ;
  if (cp)
    graphText(cp, 18, line);
  
  graphBoxDraw(boxOn, BLACK, bb ? PALEBLUE : WHITE) ;
  graphBoxDraw(boxOff, BLACK, !bb ? (n == 2 ? RED : PALEBLUE) : WHITE) ;      
}


static char *pickFirstClass = 0;
static int pickFirstBox = 0;
static char *pickFirstTemplate = 0;

void pickGetArgs(int *argcp, char **argv)
{ 
  int i;
  for (i=0; i < *argcp-1; i++)
    if (strcmp(argv[i], "-initclass") == 0)
      pickFirstClass = argv[i+1];
      
  for (i=0; i < *argcp-1; i++)
    if (strcmp(argv[i], "-inittemplate") == 0)
      pickFirstTemplate = argv[i+1];
}

static void pickDoDrawOld()
{
  int box, i, y ; float line;
  Pick pick = mainPick ;
  LAYOUT *z ;
  extern MENUOPT quovadis[] ;

  if (!graphActivate (mainGraph))
    return ; 
  graphClear() ;

  if (!mainMenu)
    {
      mainMenu = menuInitialise ("Main menu", (MENUSPEC*)quovadis) ;
      
      if (!writeAccessPossible())
	{
	  menuSuppress (mainMenu, "Save") ;
	  menuSuppress (mainMenu, "Write Access") ;
	  menuSuppress (mainMenu, "Add-Alias-Rename") ;
	}
    }

  if (writeAccessPossible())
    {
      if (isWriteAccess())
	{
	  menuSetFlags (menuItem (mainMenu, "Write Access"), MENUFLAG_DISABLED) ;
	  menuUnsetFlags (menuItem (mainMenu, "Save"), MENUFLAG_DISABLED) ;
	}
      else 
	{
	  menuUnsetFlags (menuItem (mainMenu, "Write Access"), MENUFLAG_DISABLED) ;
	  menuSetFlags (menuItem (mainMenu, "Save"), MENUFLAG_DISABLED) ;
	}
    }
  graphNewMenu (mainMenu);


  /* drawing */
  line = 0.5 ;
  if (externalServerState)
    { mainPickReport (-1) ; line += 2 ; }
  graphText ("Search: ", 0.5, line) ;
  box = graphButton("Help..",help,37.,line) ;
  graphBoxMenu (box, helpButtonMenu) ;
  if (pickFirstTemplate)
    strcpy(pick->template, pickFirstTemplate);
  else
    strcpy (pick->template,"*") ;
  pick->templateBox = graphCompletionEntry (pickClassComplete, pick->template, 
					    19, 12.0, line, 
					    TEXT_ENTRY_FUNC_CAST classDisplay) ;
  *pick->whatdoIdoText = 0 ;
 
  line += 1.5 ;
  pick->whatdoIdoBox = graphBoxStart() ;
  graphTextPtr(pick->whatdoIdoText, 18, line, 28) ;
  graphBoxEnd() ;

  line += .5 ;
  graphText ("In Class:", 0.5, line) ;
  pick->classListBox = box = graphMenuTriangle(TRUE, 9.7, line) ;
  graphBoxFreeMenu(box, (FreeMenuFunction) classListMenuCall, 
		   arrp (classList, 0, FREEOPT) ) ;
  line += 1 ;
  y = 0 ;
  firstClassBox = graphBoxStart () ;
  for (i = 1 ; i < arrayMax(classLayout) ; ++i)
    { int b = graphBoxStart () ;
      z = arrp(classLayout, i, LAYOUT) ;
      if (z->y > y) y = z->y ;
      graphText (name(z->key), z->x, z->y + line) ;
      if (pickFirstClass && (strcmp(name(z->key), pickFirstClass) == 0))
	pickFirstBox = b;
      graphBoxEnd () ;
    }
  graphBoxEnd () ;
  graphNewBoxMenu(firstClassBox, mainMenu);

  line += y + 2 ;

  graphBoxStart() ;
  graphText ("Global Search: ", 0.5, line) ;
  graphBoxEnd () ;

  *(pick->grepText) = 0 ;
  pick->grepBox = graphTextScrollEntry (pick->grepText, 255, 16, 14.5, line,
					TEXT_ENTRY_FUNC_CAST grepDisplay) ; 
  if (lexMax(_VLongText) > 2)
    pick->longGrepBox = graphButton ("Long Search", pickLongGrep, 32.,line) ;

  graphRedraw () ;
}

static void pickModeChoice (int m1, int m2, int m3)
{ 
  Pick pick = mainPick ;
  int i ;

  if (m1)
    { 
      if (pick->m1 == m1) { m1-- ; if (!m1) m1 = 1 ; }
      if (!m2)
	{ if (m1 == 1) m2 = 1 ; else m2 = 2 ; }
      i = 3 ; while (i--) 
	graphBoxDraw (pick->modeBox1 + i + 1, BLACK, WHITE) ;    
      pick->m1 = m1 ;
      i = m1 ; while (i--)
	graphBoxDraw (pick->modeBox1 + i + 1, BLACK, PALEBLUE) ;
    }
  if (m2)
    {
      graphBoxDraw (pick->modeBox2 + pick->m2, BLACK, WHITE) ;    
      pick->m2 = m2 ;
      graphBoxDraw (pick->modeBox2 + pick->m2, BLACK, PALEBLUE) ;
      if (m2 == 2)
	pickPick (0) ;
    }
}

static void pickModeChoice10 (void) { pickModeChoice(1,0,1) ; }
static void pickModeChoice20(void) { pickModeChoice(2,0,3) ; }
static void pickModeChoice30(void) { pickModeChoice(3,0,7) ; }
static void pickModeChoice1 (void) { pickModeChoice(0,1,0) ; }
static void pickModeChoice2 (void) { pickModeChoice(0,2,0) ; }
/*
static void pickModeChoice11 (void) { pickModeChoice(1,1,1) ; }
static void pickModeChoice12 (void) { pickModeChoice(1,2,2) ; }
static void pickModeChoice22 (void) { pickModeChoice(2,2,3) ; }
static void pickModeChoice32 (void) { pickModeChoice(3,2,4) ; }
*/
static void modeSearch (int box)
{  
  Pick pick = mainPick ;
  pick->longGrep = FALSE ;
  switch (10 * pick->m1 + pick->m2)
    {
    case 11: /* names, selected class */ 
      if (!pick->curr)
	{
	  messout("Please select a class or search \"in all classes\"") ;
	  return ;
	}
      classDisplayNew (pick->curr) ;
      break ;
    case 12: /* names, all classes */
      classDisplayNew (0) ;
      break ;
    case 21: /* text, selected class */  
      if (!pick->curr)
	{ 
	  messout("Please select a class or search \"in all classes\"") ;
	  return ;
	}
      strcpy (pick->grepText, pick->template) ;
      grepDisplayNew (pick->curr) ;
      break ;
    case 22: /* text, all classes */  
      strcpy (pick->grepText, pick->template) ;
      grepDisplayNew (0) ;
      break ;
    case 31:
    case 32: /* long texts */
      pick->longGrep = TRUE ;
      strcpy (pick->grepText, pick->template) ;
      grepDisplayNew (0) ;
      break ;
    }
  messStatus (0) ;
}

static void classChoice (void) ;

static void otherClassChoice (void)
{
  KEY key ;
  if (graphSelect (&key, arrp (classList, 0, FREEOPT)))
    classListMenuCall (key, 0) ;
}

static void pickDoDrawNew()
{
  int box, i, y ;
  float line;
  Pick pick = mainPick ;
  LAYOUT *z ;
  int graphWidth, graphHeight ;
  Graph old = graphActive() ;

  if (!graphActivate (mainGraph))
    return ; 

#ifndef ACEMBLY  
  if (oxgridPossible())
    graphMenu (pickNewOxMenu) ;
  else
#endif
    /* background menu */
    graphMenu (pickNewMenu) ;

  /**********************/

  graphActivate (mainGraph);
  graphClear() ;		/* do this after oxgridpossible()
				 which will do a messStatus that 
				 wants to draw into the mainwindow */

  graphFitBounds (&graphWidth,&graphHeight) ;

  /* drawing */
  line = 0.5 ;

  pick->currBox = 0 ; 

  /* server/client */
  if (externalServerState)
    { mainPickReport (-1) ; line += 2 ; }

  /* menu bar */
  if (writeAccessPossible())
    {
      if (isWriteAccess())
	{ 
	  if (sessionAutoSave (-2, 0))
	    {
				box = graphButton("Save..", saveAndKeepWriteAccess, 1., line) ;
				graphBoxDraw (box, BLACK, PALEBLUE) ;
				graphBoxMenu(box,save);
	      box = graphButton("Edit..", editSelector3,9.,line) ;
	      graphBoxFreeMenu (box, editButtonFreeAction,
			editButtonFreeMenu3) ;
	    }
	  else
	    {
				box = graphButton("Save..", saveAndKeepWriteAccess, 1., line) ;
				graphBoxMenu(box,save);
	      box = graphButton("Edit..", editSelector2,9.,line) ;
	      graphBoxFreeMenu (box, editButtonFreeAction,
			editButtonFreeMenu2) ;
	    }
	}
      else
	{ 
	  box = graphButton("Edit..", editSelector1,9.,line) ;
	  graphBoxFreeMenu (box, editButtonFreeAction,
			editButtonFreeMenu1) ;
	}
    }

  box = graphButton("Query..", querySelector,17.,line) ;
  graphBoxFreeMenu (box, queryButtonFreeAction, queryButtonFreeMenu) ;

  box = graphButton("Admin..", adminSelector,26.,line) ;
  graphBoxFreeMenu (box, adminButtonFreeAction, adminButtonFreeMenu) ;

  box = graphButton("Clear", localCleanUp,37.,line) ;
  graphBoxMenu(box, 0);

  box = graphButton("Help", help, 44.,line + 1.3) ;
  graphBoxMenu(box, 0);

  box = graphButton("Exit", exitButtonAction, 44.,line) ;
  graphBoxMenu(box, 0);

  line += 1.5 ;

  /* status */
  *pick->whatdoIdoText = 0 ;
  graphText ("Status:", 1, line) ;
  pick->whatdoIdoBox = graphBoxStart() ;
  graphTextPtr(pick->whatdoIdoText, 10, line, 32) ;
  graphBoxEnd() ;
  line += 1.8 ;

  /* search box */


  /*
  strcpy (pick->template,"*") ;
    pick->templateBox = graphCompScrollEntry (pickClassComplete, pick->template, 
					    250, 28,10.0, line, 
					    TEXT_ENTRY_FUNC_CAST modeSearch) ;
					    */
  graphText ("Search", 1, line +.1) ;
  pick->m1 = 1 ; pick->m2 = 1 ; pick->m3 = 1 ; 

  pick->modeBox2 = graphBoxStart () ;  
  box = graphButton ("in the selected class", pickModeChoice1, 8,line) ; 
  graphBoxMenu(box, 0);
  graphText ("or", 31, line +.1) ;
  box = graphButton ("in all classes", pickModeChoice2, 34,line) ;
  graphBoxMenu(box, 0);
  graphBoxEnd () ;
  graphBoxDraw (pick->modeBox2 + pick->m2, BLACK, PALEBLUE) ;

  line += 1.6 ;
  pick->modeBox1 = graphBoxStart () ;
  graphText ("objects", .5, line + .1) ;
  box = graphButton ("named", pickModeChoice10, 8,line) ;
  graphBoxMenu(box, 0);
  graphText ("or", 14.9, line + .1) ;
  box = graphButton ("related to", pickModeChoice20, 17.9,line) ;
  graphBoxMenu(box, 0);
  if (lexMax(_VLongText) > 2)
    {
      /* graphText ("+",31.2, line) ;  */
      box = graphButton ("even in long texts", pickModeChoice30, 30.2,line) ;
      graphBoxMenu(box, 0);
    } 
  else { graphBoxStart () ;  graphBoxEnd () ;}/* empty box needed for counts */

 graphBoxEnd () ;
  i = pick->m1 ;
  while (i--)
    graphBoxDraw (pick->modeBox1 + i + 1, BLACK, PALEBLUE) ;

  line += 1.5 ;
   
  if (pickFirstTemplate)
    strcpy(pick->template, pickFirstTemplate);
  else
    strcpy (pick->template,"*") ;
  pick->templateBox = graphCompScrollEntry (pickClassComplete, pick->template, 
					    250, 30,10.0, line, 
					    TEXT_ENTRY_FUNC_CAST modeSearch) ;
  line += 1.5 ;
  /* class list */
  graphLine (0,line, graphWidth, line) ; line += .7 ;

  graphText ("Select a class:", 0.5, line) ;

  box = graphButton ("Other class...", otherClassChoice, 23, line - .1);
  graphBoxFreeMenu (box,(FreeMenuFunction) classListMenuCall,
		    arrp (classList, 0, FREEOPT)) ;

  box = graphButton ("Selection", classChoice, 39, line - .1) ; 
  graphBoxMenu(box, 0);
  line += 1.3 ;
  y = 0 ;
  firstClassBox = graphBoxStart () ;
  for (i = 1 ; i < arrayMax(classLayout) ; ++i)
    { int b = graphBoxStart () ;
      z = arrp(classLayout, i, LAYOUT) ;
      if (z->y > y) y = z->y ;
      graphText (name(z->key), z->x, z->y + line) ;
      if (pickFirstClass && (strcmp(name(z->key), pickFirstClass) == 0))
	 pickFirstBox = b;
       graphBoxEnd () ;
       if (pick->curr == i)
	 {
	   pick->currBox = b ;
	   graphBoxDraw (b, BLACK,PALEBLUE) ; 
	 }
    }
  graphBoxEnd () ;

  /* the first classbox that contain the class-chooser boxes forms
     part of the background as well, so attach the bg-menu as well */
#ifndef ACEMBLY  
  if (oxgridPossible())
    graphBoxMenu (firstClassBox, pickNewOxMenu) ;
  else
#endif
    graphBoxMenu (firstClassBox, pickNewMenu) ;

  line += y + 2 ;

  *(pick->grepText) = 0 ;

  graphRedraw () ;

  graphActivate (old) ;
  return ;
} /* pickDoDrawNew */

static char* myaceGetVersionString (void)
{
  static char *cp = 0, *cq  ;
  if (!cp)
    {
      cq = strnew (aceGetVersionString(), 0) ;
      if (!strncmp(cq, "ACEDB Version",13))
	cp = strnew (messprintf ("ACEDB%s",cq+13), 0) ;
      else
	cp = strnew (cq, 0) ;
      messfree (cq) ;
    }
  return cp ;
}

void pickDraw (void)
{ 
  Graph oldGraph = graphActive () ;
  BOOL old ;

  if (!graphActivate (mainGraph))
    messcrash ("pickDraw() called before pickCreate() !");

  old = isNewLayout;
  isNewLayout = !prefValue ("OLD_STYLE_MAIN_WINDOW") ;
  if (old != isNewLayout)
    { mainPick->whatdoIdoBox = 0;
      mainPick->currBox = 0 ;
      mainPick->curr = 0 ;
    }

  if (thisSession.subDataRelease)
    {
      if (sessionDbName())
	graphRetitle(messprintf("%s, %s, %s %d-%d",
				myaceGetVersionString(),
				pickMainTitle(), sessionDbName(),
				thisSession.mainDataRelease, 
				thisSession.subDataRelease)) ;
      else
	graphRetitle(messprintf("%s, %s, %d-%d",
				myaceGetVersionString(),
				pickMainTitle(),
				thisSession.mainDataRelease, 
				thisSession.subDataRelease)) ;
    }
  else
    {
      if (sessionDbName())
	graphRetitle(messprintf("%s, %s %s",
				myaceGetVersionString(), pickMainTitle(), sessionDbName())) ;
      else if (pickMainTitle()) 
	graphRetitle(messprintf("%s",
				pickMainTitle())) ;
      else
	graphRetitle(messprintf("%s",
				myaceGetVersionString())) ;
    }

  getLayout () ;

  if (isNewLayout) 
    pickDoDrawNew() ;
  else
    pickDoDrawOld() ;

  graphActivate (oldGraph) ;
} /* pickDraw */

/************************************************************/

Graph pickCreate (void)
{
  BOOL selectFirst = !isNewLayout ;
  Pick pick = 0 ;
  
  writeAccessPossible();	/* establish access status, get
				   possible warning message out of the way */

  pick = (Pick) messalloc (sizeof(struct PickStruct));
  pick->magic = PICK_MAGIC ;
  pick->grepText = (char*) messalloc (256) ; /* free'd in pickDestroy */

  mainPick = pick ;

  pick->mainGraph = mainGraph =
    displayCreate (DtMain) ;
  graphAssociate (&PICK_MAGIC, pick) ;
  graphRegister (PICK, pickPick) ;

  graphRegister (RESIZE, resizeLayout) ;
  graphRegister (DESTROY, pickDestroy) ;
  
  selectFirst = prefValue ("OLD_STYLE_MAIN_WINDOW");

  
  pick->curr = 0 ;
  pickFirstBox = 0; /* gets set if pickFirstClass matches a class */  
  pickDraw() ;
  pickFirstClass = 0;
  pickFirstTemplate = 0; /* pickDraw gets called again sometimes */
  
  if (!pickFirstBox) 
    pickFirstBox  = 1 + firstClassBox ;    /* pick the first voc */
  else 
    selectFirst = TRUE;


  if (selectFirst) /* srk - replace this with select if desired behaviour
		      is no keyset window _except_ when -initclass 
		      flag used */
    {
      pickPick (pickFirstBox) ;
      pickPick (pickFirstBox) ;	/* call twice to achieve selection */
    }

  graphActivate(pick->mainGraph) ;

  return pick->mainGraph ;
} /* pickCreate */

/*************************************/
/**** Auto completion system *********/

void mainKeySetComplete (KEYSET ks, char *text, int len)
{
  int n, i ;
  char *cp, *cq ;
  KEYSET ks2 ;
  
  if (!keySetMax(ks))		/* no match */
    { messbeep() ;
      goto ok ;
    }
  else if (keySetMax(ks) == 1)	/* unique match */
    { strncpy (text, name(keySet(ks,0)), len-1) ;
    goto ok ;
    }
				/* else multiple matches */
  n = strlen(text) ;
  for (i = 0 ; i < n ; ++i)
    text[i] = 0 ;

  ks2 = keySetAlphaHeap(ks,keySetMax(ks)) ;
  cp = name (keySet (ks2, 0)) ;
  cq = name (keySet (ks2, keySetMax(ks2)-1)) ;
  keySetDestroy(ks2) ;

  for (i = 0 ; *cp && *cp == *cq && i < len-1 ; ++i, ++cp, ++cq)
    text[i] = *cp ;
  text[i] = 0 ;

  if (i > n)		/* length increase */
    goto ok ;
 
ok:
  if (mainPick)
    { Graph oldGraph = graphActive() ;
      
      if (ks != mainPick->oblist) /* probably a useless kludge */
	{ keySetDestroy (mainPick->oblist) ;
	  mainPick->oblist = ks ;
	}
      if (graphActivate (mainPick->mainGraph))
	obDisplay (mainPick, FALSE, 3) ; /* no automatic display() if 1 */
    
      graphActivate (oldGraph) ;
    }
}

/****************/

   /* return number of possible completions */
int ksetClassComplete (char *text, int len, int classe)
{
  static KEYSET ks = 0 ;
  static int previousClass = 0 ;
  static int myKeySetId = 0 ;
  static char buffer[80] ;
  KEYSET ks2 ;
  char *cp, *cq ;
  int i, j , n ;

  if (myKeySetId != keySetExists(ks))
    ks = 0 ;
  if (!text)
    { keySetDestroy(ks) ;
      return 0 ;
    }

  if (ks)
    {
      if (previousClass != classe)
	keySetDestroy(ks) ;
      else
	{     /* check if search field just got longer */
	  cp = buffer ; cq = text ; n = 80 ; /* size of buffer */
	  while (--n && *cp && *cq == *cp){cp++; cq++;} 
	  if (!n || *cp)
	    keySetDestroy (ks) ;
	}
    }
  
  strncpy(buffer, text, 79) ;
  if (!ks)
    {
      if (class(classe) == _VClass)
	ks = query (0, messprintf ("FIND %s IS \"%s*\"", 
				   name(classe), text)) ;
      else			/* assume a true classe */
	ks = query (0, messprintf ("FIND %s IS \"%s*\"", 
				   pickClass2Word(classe), text)) ;
    }
  else
    { 
      ks2 = arrayCreate(keySetMax(ks), KEY) ;
      text[len-2] = 0 ;		/* prevents running off end */
      cp = text + strlen(text) ;
      *cp = '*' ; *(cp+1) = 0 ;
      for (i=0, j=0 ; i < keySetMax(ks) ; i++ )
	if (pickMatch (name(keySet(ks,i)), text))
	  keySet(ks2,j++) = keySet(ks,i) ;
      *cp = 0 ;			/* remove * from end */
      keySetDestroy(ks) ;
      ks = ks2 ;
    }

  mainKeySetComplete (ks, text, len) ;
    
  myKeySetId = keySetExists(ks) ;
  previousClass = classe ;
  strncpy(buffer, text, 79) ;
  return keySetMax(ks) ;
}

/****************************************/
/****************************************/

static Graph cGraph = 0;
static int 
 cFirstBox = 0, cFirstBoxEnd = 0 , yLimit, firstXX, firstYY, 
  cSecondBox = 0, cSecondBoxEnd = 0 ,
  cZeroBox = 0, cZeroBoxEnd = 0 ;
static BitSet cBs = 0 ;
static KEYSET cBox2key = 0 ;
static KEY draggedKey = 0 ;
static Array cLayout = 0 ;

static void cDestroy (void)
{
  cGraph = 0 ; cBs = 0 ; cBox2key = 0 ; cLayout = 0 ;
}

static void cDraw (void) ;

static int cLayoutOrder(const void *va, const void *vb)
{
  const LAYOUT *lla = (const LAYOUT*)va,  *llb = (const LAYOUT*)vb ;
  if (lla->key)
    {
      if (llb->key)
	return lla->y != llb->y  ? lla->y - llb->y  : lla->x - llb->x ;
      return -1 ;
    }
  return llb->key ? 1 : 0 ;
}


static void cApply (void)
{
  int i, j ;
  LAYOUT *ll1, *ll2 ;

  for (i = 0 , j = 1 ; i < arrayMax(cLayout) ; i++)
    { 
      ll1 = arrp(cLayout, i, LAYOUT) ;
      if (bitt(cBs,KEYKEY(ll1->key)))
	{ ll2 = arrayp(classLayout, j++, LAYOUT) ;
	*ll2 = *ll1 ;
	}
    }  
  arrayMax(classLayout) = j ;
  graphDestroy () ;
  if (mainPick) mainPick->curr = 0 ; /* layout may change */
  pickDoDrawNew() ;
}

static void cSave (void)
{
  int i, level ;
  FILE *f = 0 ;
  LAYOUT *ll ;  
  static char fileName[FIL_BUFFER_SIZE] , dirName[DIR_BUFFER_SIZE] ;
  static BOOL firstpass = TRUE ;

  arraySort (cLayout, cLayoutOrder) ;

  if (firstpass)
    { firstpass = FALSE ;
      strcpy (fileName,"layout.wrm") ;
      strcpy (dirName,sessionFilName ("", 0, 0)) ;
      strcat (dirName,"/wspec") ;
    }
  f = filqueryopen(dirName, fileName, "", "w",
			    "Where do you write the new layout file ?") ;
  if (!f) return ;
  level = freeOutSetFile(f) ;

  for (i = 0 ; i < arrayMax(cLayout) ; i++)
    { 
      ll = arrp(cLayout, i, LAYOUT) ;
      if (bitt(cBs,KEYKEY(ll->key)))
	{ 
	  freeOutxy (name(ll->key), ll->x, ll->y) ;
	}
    }  
  freeOut ("\n") ;
  freeOutClose (level) ;
  filclose (f) ;
}


static MENUOPT cButtons[] = 
{
  {cApply, "Apply"},
  {cSave, "Save"},
  {graphDestroy,"Cancel"},
  {0,0}
} ;

static void cBoxDrag (float *x, float *y, BOOL isDone)
{
  int i, n, ii ;
  LAYOUT *ll = 0 ;

  if (!isDone)
    return ;
  *x += .2 ; *y -= .2 ;/* before rounding to integer,  easier interface */
  n = arrayMax(cLayout) ;
  if (*y > yLimit)
    { 
      if (class(draggedKey))  /* not on tools */
	bitUnSet (cBs, KEYKEY(draggedKey)) ;  
      else
	{
	  if (draggedKey == 10)  /* suppress the whole column */
	    {
	      for (i=0; i <= n ; i++)
		{
		  ll = arrayp(cLayout, i, LAYOUT) ;
		  if (ll->key && ll->x > 1)
		    ll->x -- ;
		}
	    }
	  if (draggedKey >= 100)  /* suppress the whole line */
	    {
	      ii = draggedKey - 100 - 1 ;
	      for (i=0; i <= n ; i++)
		{
		  ll = arrayp(cLayout, i, LAYOUT) ;
		  if (ll->key && ll->y == ii)
		    bitUnSet (cBs, KEYKEY(ll->key)) ;  
		  if (ii > -1 && ll->key && ll->y > ii)
		    ll->y -- ;
		}
	    }
	}
    }
  else
    {
      if (class(draggedKey))  /* not on tools */
	{
	  ii = *y - firstYY ;
	  if (ii < -1) ii = -1 ;
	  bitSet (cBs, KEYKEY(draggedKey)) ;
	  if (ii<0) /* shift everybody down one line */
	    { ii ++ ;
	      for (i=0; i <= n ; i++)
		{
		  ll = arrayp(cLayout, i, LAYOUT) ;
		  if (ll->key)
		    ll->y += 1 ;
		}
	    }
	  if (*x < firstXX)   /* shift everybody right */
	    {
	      for (i=0; i <= n ; i++)
		{
		  ll = arrayp(cLayout, i, LAYOUT) ;
		  if (ll->key)
		    ll->x += firstXX - *x ;
		}
	      *x = firstXX - 1 ;
	    }
	  for (i=0; i <= n ; i++)
	    {
	      ll = arrayp(cLayout, i, LAYOUT) ;
	      if (!ll->key || ll->key == draggedKey)
		{ 
		  ll->key = draggedKey ;
		  ll->x = *x - firstXX + 1 ;
		  ll->y = ii ;
		  break ;
		}
	    }
	}
      else
	{
	  switch (draggedKey)
	    {
	    case 1:   /* add empty line */
	      ii = *y - firstYY ; 
	      for (i=0; i <= n ; i++)
		{
		  ll = arrayp(cLayout, i, LAYOUT) ;
		  if (ll->key && ll->y >= ii)
		    ll->y++ ;
		}
	      break ;
	    }
	}
    }
  cDraw() ;
}

static void cPick (int box)
{
  int n ;

  n = box - cFirstBox - 1 ;
  if (box > cZeroBox && box < cZeroBoxEnd)
    {
      n = box - cZeroBox - 1 ; 
      if (n<0 || !keySetExists(cBox2key) || n >= keySetMax(cBox2key))
	return ;
      draggedKey = keySet(cBox2key,box) ;
      graphBoxDrag (box, cBoxDrag) ;
    }
  else if (box > cFirstBox && box < cFirstBoxEnd)
    {
      n = box - cFirstBox - 1 ;
      if (n<0 || !keySetExists(cBox2key) || n >= keySetMax(cBox2key))
	return ;
      draggedKey = keySet(cBox2key,box) ;
      graphBoxDrag (box, cBoxDrag) ;
    }
  else if (box > cSecondBox && box < cSecondBoxEnd)
    {
      n = box - cSecondBox - 1 ;
      if (n<0 || !keySetExists(cBox2key) || n >= keySetMax(cBox2key))
	return ;
      draggedKey = keySet(cBox2key,box) ; 
      graphBoxDrag (box, cBoxDrag) ;
    }
}

static void cDraw (void)
{
  int box, i, j, x, xmax, ymax ; float line = 1, y = 0, dy = 0 ;
  KEYSET ks, ksa ;
  KEY key ;
  int graphWidth, graphHeight ;
  LAYOUT *ll = 0 ;

  graphFitBounds (&graphWidth,&graphHeight) ;
  line += 1 ;
  ks = query (0, "FIND Class NOT Buried") ;
  ksa = keySetAlphaHeap(ks, keySetMax(ks)) ;
  keySetDestroy (cBox2key) ;
  cBox2key = arrayHandleCreate(30,KEY,graphHandle()) ;
  y = 0 ;
  graphClear () ; 
  line = 1.2 ;
  graphButtons (cButtons, 1, line, 100) ;

  line += 1.8 ;
  graphText("Drag the lines and classes with the left mouse button", 1,line); 
  line += 2 ;
  graphText("In yellow is your current class list:", 1,line) ; line += 1.5 ;

  arraySort (cLayout, cLayoutOrder) ;

  /* go to next line to prevent overlappings */
  xmax = 0 ; ymax = 0 ; x = 0 ; y = 0 ; dy = 0 ;
  for (i = 0 , j = 0 ; i < arrayMax(cLayout) ; ++i)
    { 
      ll = arrp(cLayout, i, LAYOUT) ;
      if (!ll->key || !bitt(cBs,KEYKEY(ll->key)))
	continue ;
      if (ll->y > y) { x = 0 ; y = ll->y ; }
      if (ll->x < x) { x = 0 ; dy++ ; } 
      if (ll->x + strlen(name(ll->key)) + 2 > xmax)
	xmax = ll->x + strlen(name(ll->key)) + 2 ;  
      if (ll->key && ll->y + 2 > ymax)
	ymax = ll->y + 2 ;
      x = ll->x + strlen(name(ll->key)) + 2 ;
      ll->y += dy ;
    }
  /* search min x, y , max = */
  x = y = 1000000 ; 
  for (i = 0 , j = 0 ; i < arrayMax(cLayout) ; ++i)
    { 
      ll = arrp(cLayout, i, LAYOUT) ;
      if (!ll->key || !bitt(cBs,KEYKEY(ll->key)))
	continue ;
      if (ll->y < y) y = ll->y ; 
      if (ll->x < x) x = ll->x ;  
    }
  /* adjust min x, y , to be 2, 0 */
  x -= 2 ;
  for (i = 0 , j = 0 ; i < arrayMax(cLayout) ; ++i)
    { 
      ll = arrp(cLayout, i, LAYOUT) ;
      if (!ll->key || !bitt(cBs,KEYKEY(ll->key)))
	continue ;
      ll->x -= x ;
      ll->y -= y ;
    }
  /* search xmax ymax */
  xmax = 0 ; ymax = 0 ; 
  for (i = 0 , j = 0 ; i < arrayMax(cLayout) ; ++i)
    { 
      ll = arrp(cLayout, i, LAYOUT) ;
      if (!ll->key || !bitt(cBs,KEYKEY(ll->key)))
	continue ;
      if (ll->x + strlen(name(ll->key)) + 2 > xmax)
	xmax = ll->x + strlen(name(ll->key)) + 2 ;  
      if (ll->y + 2 > ymax)
	ymax = ll->y + 2 ;
    }
  if (!ymax) ymax = 2 ; 
  if (!xmax) xmax = 12 ;

  cZeroBox = graphBoxStart () ;  graphBoxEnd() ; 
  /* vertical line removal, rather, i will always stack left and top 
  keySet(cBox2key,box = graphBoxStart ()) = 10 ;
  for (j = 0 ; j <= ymax ; j++)
    graphText("|", 1, line+j) ;
  graphBoxEnd() ;
  graphBoxDraw(box,BLACK,YELLOW) ;
  */
  for (j = 0 ; j <= ymax ; j++)
    {
      keySet(cBox2key,box = graphBoxStart ()) = j + 100 ; /* allow 100 other tools */
      for (i = 0 ; i < xmax ; i++)
	graphText("_", i+2, line+j) ;
      graphBoxEnd() ;
      graphBoxDraw(box,BLACK,PALEYELLOW) ;
    }
  cZeroBoxEnd = graphBoxStart () ;  graphBoxEnd() ; 

  cFirstBox = graphBoxStart () ;  graphBoxEnd() ;  
  firstXX = 3 ; firstYY = line + 1 ;
  for (i = 0 , j = 0 ; i < arrayMax(cLayout) ; ++i)
    { 
      ll = arrp(cLayout, i, LAYOUT) ;
      if (!ll->key || !bitt(cBs,KEYKEY(ll->key)))
	continue ;
      keySet(cBox2key,box = graphBoxStart ()) = ll->key ;
      graphText(name(ll->key),ll->x + 2, line + ll->y + 1) ;
      graphBoxEnd() ;
    }  
  cFirstBoxEnd = graphBoxStart () ;  graphBoxEnd() ; 
  line += ymax + 2 ; 
  yLimit = line ;
  graphLine (0,line,200,line) ;
  line += 1 ;
  
  graphText("These are the remaining classes:", 1,line) ;
  cSecondBox = graphBoxStart () ;  graphBoxEnd() ; 
  keySet(cBox2key,box = graphBoxStart ()) = 1 ;
  graphText ("___add_line___", 38,line) ;
  graphBoxEnd () ;
  graphBoxDraw (box, BLACK, PALEYELLOW) ;
  line += 1.5 ;
  for (i = 0 , xmax = 1 ; i < keySetMax (ksa) ; i++)
    { 
      key = keySet (ksa, i) ; 
      if (bitt(cBs, KEYKEY(key)))
	continue ;
      x = xmax + strlen(name(key)) ;
      if (x > graphWidth)
	{ xmax = 1 ; line++ ; }
      keySet(cBox2key,box = graphBoxStart ()) = key ;
      graphText (name(key), xmax, line) ;
      graphBoxEnd () ;
      xmax += strlen(name(key)) +1 ;
      xmax = ((xmax + 7)/8) * 8 ;
      if (bitt(cBs, KEYKEY(key)))
	graphBoxDraw (box, BLACK, PALEBLUE) ;
    } 
  cSecondBoxEnd = graphBoxStart () ;  graphBoxEnd() ; 
  graphRedraw () ;
  keySetDestroy (ks) ;
  keySetDestroy (ksa) ;
}

static void classChoice (void) 
{
  Graph old = graphActive () ;
  int i, j ;
  LAYOUT *ll1, *ll2 ;

  if (!graphActivate (cGraph))
    {
      cGraph = graphCreate (TEXT_FIT,"Class Selector", .3, .3, .6, .5) ;
      graphHelp ("Class Selector") ;
      graphRegister (DESTROY, cDestroy) ;
      graphRegister (RESIZE, cDraw) ;
      graphRegister (PICK, cPick) ;
 
      bitSetDestroy (cBs) ;
      cBs = bitSetCreate (80, graphHandle()) ;
      cLayout = arrayHandleCreate (50, LAYOUT,graphHandle()) ;
      for (i = 1 , j = 0 ; i < arrayMax(classLayout) ; ++i)
	{ 
	  ll1 = arrp(classLayout, i, LAYOUT) ;
	  bitSet (cBs, KEYKEY(ll1->key)) ;
	  ll2 = arrayp(cLayout, j++, LAYOUT) ;
	  *ll2 = *ll1 ;
	}  
    }
  cDraw() ;
  graphActivate (old) ;
}

