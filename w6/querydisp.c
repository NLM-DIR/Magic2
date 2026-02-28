/*  File: querydisp.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 **  Constructs general queries.                              **
 ***    Gary Aochi modified version of querydisp.c
 * Exported functions:
 **  queryCreate
 * HISTORY:
 * Last edited: Dec 11 10:20 1998 (fw)
 * * Jun  3 12:19 1996 (rd)
 * * Feb 09 15:30 1993 (aochi): split off qbe to qbedisp.c and query_.h
 *                              and cleaned up some code in querydisp.c
 * * Jun 14 10:14:40 PDT 1993 (aochi) : split querydisp.c into 2 files:
 *                                      querydisp.c (ctrl), querybuild.c
 * * Jan 18 12:55 1993 (mieg): Changed MAGIC and added some static def
 * (original querydisp.c created Mon Dec 30 19:05:47 1991 (mieg))
 * Created: Jun 14 10:14:40 PDT 1993 (aochi)
 *-------------------------------------------------------------------
 * gha 03/02/93  (changed selection list to keyset)
 * gha 03/09/93  (querpick ifs, fixed queryCommandEdit(off1.9)) 
 * gha 04/09/93  added > to indicate last query, hit rtn to go thru queries
 * gha 04/28/93  jean cleaned up stuff, temporarily have hidden stuff,
 *               changed >? to Find, > to Follow, fixed classlist
 *               added code to unhighlight some buttons
 * gha 05/25/93  added jean's new Search stuff and changes he made.
 * gha 06/14/93  reduced new querydisp.c from querydisp.c
 */

 /* $Id: querydisp.c,v 1.4 2009/06/22 19:35:05 mieg Exp $ */

#include <ctype.h>
#include "acedb.h"
#include "lex.h"
#include "bs.h"         /* used for bsTagsInClass() procedure */
#include "query_.h"
#include "query.h"
#include "parse.h"
#include "command.h"
#include <wac/ac.h>
#include <wac/ac_.h>
#include <wac/acinside_.h>

typedef struct QCTRLSTUFF
{ int   magic;        /* == MAGIC*/
  Graph graph ;
  char fileName[FIL_BUFFER_SIZE] ;
  char dirName[DIR_BUFFER_SIZE] ;
  int pgmBox ;
  int curr ;
  int last_find ;
  FILE *fil ;
  BOOL modif ;
  Array pgm ; /* array of Stack, each one a program line */
  char flagmask;       /* control mask, 8 bits of togglable settings */
} *QUER ;

Graph querGraph = 0 ;
KEYSET query_undo_set ;

#define FILENAME_SIZE  24

static int MAGIC = 92443 ;	/* also use address as graphAss handle */

#define NO_AUTOSV     0x01

#define SETFLAG(y)  quer->flagmask = (quer->flagmask | y)
#define TESTFLAG(y) (quer->flagmask & y)
#define RMFLAG(y)   quer->flagmask = (quer->flagmask & ~y)

#include "display.h"	/* must come after QUER definition */
     
static void querDisplay(void) ;
static void querPick (int k) ;
void queryCreate (void) ;
     
#define QUERGET(name)  QUER quer ; \
                       \
                       if(!graphActivate(querGraph)) \
                          messcrash("%s can't find querGraph", name); \
                       if (!graphAssFind (&MAGIC, &quer)) \
                           messcrash ("%s can't find graph", name) ; \
                       if (quer->magic != MAGIC) \
                           messcrash ("%s received a wrong pointer", name) ; \
                       displayPreserve()  
     

/*********** stuff local to the quer routine and initialiser **********/
/************************************************************************/
     
static Array dirMarks = 0 ;
static Stack dirStack = 0 ;

/************************* AUTOSAVE STUFF **************************/

/* gha 11-2-92: will save to date-time-stamp file if no filename specified */

static char auto_filename[FILENAME_SIZE];

static char *autosaveCreateName(void)
{
  char timeBuf[25] ;
  char *tstamp = timeShow (timeNow (), timeBuf, 25) ;

  tstamp[13] = 0 ;		/* truncate at 10's of minutes */
  
  return tstamp ;
}

/* if active, save the query commands to a file when "Search" invoked */
static void autoSave(void)
{
  int i ;
  FILE *fil ;
  QUERGET("autoSave");
  
  if (TESTFLAG(NO_AUTOSV))
    return;

  if (!graphExists (quer->graph) || !quer->pgm 
      || arrayMax(quer->pgm) <= 1)  /* nothing to save */
    return ;
  
  strcpy (auto_filename, autosaveCreateName()) ;

  if (!strcmp (quer->fileName, auto_filename)) 
    { messout ("Warning: autosave-query filename conflict: %s.\n"
	       "Using different autosave name: %sX.", 
	       auto_filename, auto_filename) ;
      strcat (auto_filename, "X");
    }

  fil = filopen (messprintf ("wquery/%s", auto_filename), "qry", "a+") ;
  if (!fil)
    return ;
  
  for(i=0; i < arrayMax(quer->pgm) ; i++)
    fprintf (fil,"%s\n", ARR2STRING(quer->pgm, i));
  
  filclose(fil) ;
}

/*********** action routines (may be externally called) ***********/

/* external load (from file) */
static void querParseFile(void)
{
  Array a ;
  int i, j, box ;
  Array s ;
  Array pgm ;
  char c ;
  FILE *fil ;
  QUERGET("querGetFile");
  

  fil = filqueryopen (quer->dirName, quer->fileName, "qry",
		      "r", "Choose a Query File") ;

  if (!fil || !graphExists (quer->graph))
    return ;

  box = quer->curr;
  pgm = quer->pgm ;
  if(!pgm)
    pgm = quer->pgm = arrayCreate(5,Stack) ;
  if (arrayMax(pgm) > 0)   /* gha 2-24 append, not overwrite */
    i = arrayMax(pgm) - 1;
  else 
    i = 0;

  s = 0 ;
  j = 0 ;
  while (c = fgetc(fil), c != (char)EOF)  
    {
      if (c == '\n' || !s) {
	if (s && j)
	  array(s,j,char) = 0 ;
	s = array(pgm,i,Array) ;
	if (s)
	  arrayMax(s) = 0 ;
	else
	  s = array(pgm,i,Array) = arrayCreate(QBUFF_MULT*BUFFER_SIZE, char) ;
	j = 0 ;
	i++ ;
      }
      if(isspace((int)c)) 
	{
	  if(j) array(s,j++,char) =  ' ' ; 
	}
      else
	if (c!='\n' && isprint((int)c))
	  array(s,j++,char) = c ;
    }
  if (s)
    {
      if (j)
	array(s,j,char) = 0 ;
      else {
	arrayDestroy(s) ;
	array(pgm,--i,Array) = 0 ;
	arrayMax(pgm) = i ;
      }
    }
  filclose(fil) ;
  /* One empty entry for a new instruction */ 
  array (pgm,i,Array) = a = arrayCreate(QBUFF_MULT*BUFFER_SIZE, char) ;
  array(a,0,char) = 0 ;
  quer->last_find = -1 ;
  querDisplay () ;
  querPick(box);
}

/* external save (to file) */
static void querSaveFile(void)
{
  int i ;
  Array s ;
  Array pgm ;
  FILE *fil ;
  QUERGET("querSaveFile");
  
  pgm = quer->pgm ;
  if(!pgm || arrayMax(pgm) <= 1)  {
    messout("No program to save") ;
    return ;
  }
  
  fil = filqueryopen (quer->dirName, quer->fileName, "qry", 
		      "w","Choose a Query File") ; 

  if (!fil || !graphExists (quer->graph))
    return ;
  
  pgm = quer->pgm ;
  
  for(i=0; i + 1 < arrayMax(pgm); i++)
    {
      s = array(pgm,i,Array) ;
      if (s)
	fprintf (fil,"%s\n", arrp(s,0,char)) ;
    }
  
  filclose(fil) ;
}


/**********************************************************/

/* Called by QUIT */
static void querDestroy (void)
{    
  Array pgm ;
  int i ;
  char * cp ; 
  QUER quer ;
  
  if ( graphAssFind (&MAGIC, &quer)
      && (quer->magic == MAGIC))  {
    cp =  quer->fileName ;
    if( quer->pgm 
       && quer->modif
       && *quer->fileName
       && strlen(quer->fileName)
       && messQuery (messprintf("%s not saved, Save ?", cp)) )
      querSaveFile() ;
    
    quer->magic = 0 ;

    if ((pgm = quer->pgm)) {
      i = arrayMax(pgm) ;
      while(i--) 
	arrayDestroy(arr(pgm,i,Array)) ;
      arrayDestroy(pgm) ;
    }
    if (!graphExists(qbuildGraph) && !graphExists(qbeGraph))
      arrayDestroy(query_undo_set);
    messfree (quer) ;
    
    stackDestroy(dirStack) ;
    arrayDestroy(dirMarks) ;
    
    if(graphActivate(querGraph))
      graphDestroy () ;
    
    querGraph = 0 ;
  }
}

/* modified gha 9/4/92
 * Mouse pick function for query graph window...based upon 
 * format of boxes as called in queryCreate()...caveat programmer.
 */
static void querPick (int k)
{    
  int i ;
  QUERGET("querPick");
  
  if (k <= 0)
    return ;

  if (k > quer->pgmBox ) {    /* for query entries */
    i = (k - quer->pgmBox - 1)/2;
    if (i < arrayMax(quer->pgm)) {
      quer->curr = k ;
      graphTextScrollEntry (ARR2STRING(quer->pgm, i), 0, 0, 0, 0, 0) ;
    }
  }
  else
    return;
  graphActivate(quer->graph) ;
}

/*******************  add 03.05.1992       kocab DKFZ  **********/
/****************************************************************/

#ifdef FRED

#include "igdevent.h"
extern BOOL fredExists ;
static void querySend (void)
{
  fredscEvent      event;
  int  i;
  FredId    j;
  Array s ;
  Array pgm ;
  static Stack buf = 0 ;
  QUERGET("querySend");
  
  
  if (!fredExists)
    { messout ("aceServer is not initialised, sorry") ;
      return ;
    }
  
  pgm = quer->pgm ;
  if(!pgm || arrayMax(pgm) <= 1)
    { messout("No program to send") ;
      return ;
    }
  
  buf = stackReCreate(buf, 200) ;
  for(i=0; i < arrayMax(pgm) - 1; i++) 
    {
      s = array(pgm,i,Array) ;
      if (arrayExists(s))
	{ catText(buf, arrp(s,0,char)) ; catText(buf,";") ; }
/*	catText(buf, arrp(s,0,char)) ;   gha debug from jtm 5/25 */
    }
  catText(buf, "\000") ;
  
  event.typeE = ACE_QUERY;
  event.typeD = DATA_IN_MEM_BUF;
  event.replyTo = 0;
  event.buffer.buff.bptr = stackText(buf,0) ;
  event.buffer.buff.blen = strlen(stackText(buf,0)) ;
  j = fredSendEvent(&event) ;
  printf("fredsendevent = %i \n",(int)j   );
  fredKillEvent(&event);
  if (wait_to_events (  &event, DATA_OK | DATA_NOT_OK ) < 0 ) 
    messout(" error in communication with face");
  else 
    if (event.typeE == DATA_NOT_OK ) 
      messout("data not received by face");
  fredKillEvent(&event);
}
#else
static void querySend (void)
{ messout("To send query, please recompile with -D FRED") ;
}
#endif
/*******************  add 03.05.1992       kocab DKFZ  **********/
/****************************************************************/



/*************** alternative aceclient/server external query **************/

/****************************************************************/
/*************** beginning of menu stuff *******************/
static void clearQueryCommands(void);

static void qbeCallup(void) { 
  qbeCreate();
  querDisplay();  /* gha added to unhighlight by example button */
}

static void qbuildCallup(void) { 
  qbuildCreate();
  querDisplay();  /* gha added to unhighlight builder button */
}

static void querButtonUndo(void) {
  queryCommandUndo();
  querDisplay();  /* gha added to unhighlight undo button */
}

static void qNewKeySet (void)
{
  newKeySet ("Query Answer") ;
}

static MENUOPT querButtons[] =
{
  { querDestroy,        "Quit" },
  { help,               "Help" },
  { graphPrint,         "Print" },
  { clearQueryCommands, "Clear" },
  { qNewKeySet,         "New KeySet" },
  { qbeCallup,          "By Example" },
  { qbuildCallup,       "Builder" },
  { querParseFile,      "Load" },
  { querSaveFile,       "Save" },
  { queryCommandEdit,   "Search" },
  { querButtonUndo,     "Undo" },
  { 0, 0 }
} ;

static void querAutoSv(void);

/* Look at following functions to make sure not to disturb [#] indexing */
/* querAutoSv and MENU_AUTO */
#define MENU_AUTO  13
static MENUOPT querMenu[] =
{
  { querDestroy,        "Quit" },
  { help,               "Help" },
  { graphPrint,         "Print Window" },
  { menuSpacer,         "" },
  { querDisplay,        "Refresh" },
  { clearQueryCommands, "Clear Entries" },
  { qNewKeySet,         "New KeySet" },
  { qbeCreate,          "Query By Example" },
  { qbuildCreate,       "Query Builder" },
  { menuSpacer,         "" },
  { querParseFile,      "Load From File" },
  { querSaveFile,       "Save To File" },
  { querySend,          "Send" },
  { querAutoSv,         "Enable Autosave" }, /* MUST BE #13, see querAutoSv(), MENU_AUTO */
  { menuSpacer,         "" },
#ifdef IGD
  { queryCommandEdit,   "Search Locally" },
  { querySend,          "Search IGD" },
  { querySend2,         "Search Server" },
#else
  { queryCommandEdit,   "Search" },
#endif
  { queryCommandUndo,   "Undo Search" },
  { 0, 0 }
} ;

static void querAutoSv(void) 
{
  QUERGET("querAutoSv");

  if (TESTFLAG(NO_AUTOSV))
    {
      RMFLAG(NO_AUTOSV);
      querMenu[MENU_AUTO].text = "Disable AutoSave";
    }
  else
    {
      SETFLAG(NO_AUTOSV);
      querMenu[MENU_AUTO].text = "Enable AutoSave ";
    }
  graphMenu(querMenu) ;  
}


/*************** end of menu stuff *******************/

static void clearQueryCommands(void)
{
  QUERGET("clearQueryCommands");

  if (quer->pgm) {
    strcpy(ARR2STRING(quer->pgm, 0), "");
    arrayMax(quer->pgm) = 1;
  }
  querDisplay();
}


/* open/show query-interface window */ 
#define QC_OFFSET 7  /* should be at least 3 */
static void querDisplay(void)
{
  int i = 0;
  QUERGET("querDisplay");
  
  graphRegister (PICK,(GraphFunc)querPick) ;
  graphActivate(quer->graph) ;
  graphClear() ;
  
  graphButtons (querButtons, 1, 1, 75) ;  /* put buttons on screen */

  graphTextFormat(BOLD);
  quer->pgmBox = graphBoxStart() ;
  
  if(quer->pgm)  {
    graphTextFormat(BOLD);
    graphText("Direct Query Commands:", 1.5 , QC_OFFSET - 2.2);
    graphTextFormat(PLAIN);
    graphText("Click line,   select as first word implies BQL syntax",
	       24.5, QC_OFFSET - 2.5);
    graphText("(Press Return key to execute, F4 to interrupt.)",   
	       24.5, QC_OFFSET - 1.5);
    for(i=0; i<arrayMax(quer->pgm); i++)  {
      if (i == quer->last_find)
	graphText(">", 1, i*1.5 + QC_OFFSET);
      graphText("Query:", 2, i*1.5 + QC_OFFSET);
      quer->curr =
	graphTextScrollEntry(ARR2STRING(quer->pgm, i),
			     QBUFF_MULT*BUFFER_SIZE - 1, 183, 9, 
			     i*1.5 + QC_OFFSET, 
			     TEXT_ENTRY_FUNC_CAST queryCommandEdit) ;
    }
  }
  graphBoxEnd() ;
  
  graphTextBounds (180, (int)(QC_OFFSET + 1.5*i + 2)) ;        /* see queryCreate */
  graphRedraw() ;
}


/*********************************************************************/
/********************  public routines   *****************************/


/* show formed-query in window's query command line */
void queryWindowShowQuery(char *buffer) 
{
  int box;
  QUERGET("queryWindowShowQuery");

  box = quer->curr;
  if (strlen(buffer) >= QBUFF_MULT*BUFFER_SIZE)
    messout("Warning, the query command may have exceeded the buffer capacity: %d.\n",
	    QBUFF_MULT*BUFFER_SIZE-1);

  strncpy(ARR2STRING(quer->pgm, (quer->curr - quer->pgmBox -1)/2), 
	  buffer, QBUFF_MULT*BUFFER_SIZE-1);
  quer->curr = graphTextScrollEntry(ARR2STRING(quer->pgm, 
					       (quer->curr-quer->pgmBox- 1)/2),
				    0, 0, 0, 0, 0);
  querPick(box);
}


/* Calls query() on query syntax and re-displays
 * window with new blank entry line, setting current pointer to that line.
 */
void queryCommandEdit (void)
{ 
  Array a ;
  int i ;
  KEYSET oldSet, newSet ;
  void *look ;
  char *querText, *cp;  
  QUERGET("queryCommandEdit");
  
  quer->modif = TRUE ;
  
  i = (quer->curr - quer->pgmBox - 1) / 2 ;
  quer->last_find = i;

  if(i + 1 == arrayMax(quer->pgm)  &&
     arrayExists(array(quer->pgm,i,Array)))
    { cp = arrp(arr(quer->pgm,i,Array), 0, char) ;
      if (cp)
	while(*cp && *cp == ' ') cp++ ;
      if(cp && *cp) /* but only if line is not empty */
	{
	  /* One empty entry for a new instruction */ 
	  array(quer->pgm,i+1,Array) = 
	    a = arrayCreate(QBUFF_MULT*BUFFER_SIZE, char);
	  array(a,0,char) = 0 ;
	  querDisplay() ;
	}
    }
  
  /* Execute */
  if (arrayExists(array(quer->pgm,i,Array))) {
    querText = arrp(arr(quer->pgm,i,Array), 0, char) ;
    cp = querText ;
    if (cp)
      while(*cp && *cp == ' ') cp++ ;
    if(!*cp)
      return ; /* do nothing on empty lines */
    if(!keySetActive(&oldSet, &look))
      newKeySet("Query Answer") ;
    if(keySetActive(&oldSet, &look))  { 
      messStatus("Searching") ;
      if (!strncasecmp (querText, "select ",7))
	{
	  AC_HANDLE h = ac_new_handle () ;
	  const char *errors = 0 ;
	  AC_DB db = ac_open_db (0, &errors) ;
	  AC_KEYSET aksOld = ac_new_keyset (db , h) ;
	  aksOld->ks = oldSet ;
	  AC_TABLE tbl = ac_bql_table (db, querText, aksOld, 0, &errors, h) ;
	  aksOld->ks = 0 ; /* so oldSet will not be destroyed */
	  AC_KEYSET aks = ac_table_keyset (db, tbl,0, h) ;
	  newSet = keySetCopy (aks->ks) ; 
	  ac_free (h) ;
	}
      else
	newSet = query(oldSet, querText) ;
      if (newSet != oldSet )
	{ if(keySetExists(query_undo_set) &&
	     query_undo_set != oldSet )
	    keySetDestroy(query_undo_set) ;
	  query_undo_set = oldSet ;
	}
      keySetShow(newSet,look) ;
      if (newSet && keySetMax(newSet))  /* update class, taglist */
	{ qbuildNewClass(class(keySet(newSet, 0)));
	  autoSave();
	}
      querDisplay();
      querPick(quer->pgmBox+ 2*(i+1) + 1);  /* go to next line */
    }
    else
      messout("First select a keySet, thank you") ;
  }
}


/* Called by UNDO  */
void queryCommandUndo (void)
{
  KEYSET oldSet, newSet ;
  void   *look ;

  if (keySetExists(query_undo_set) &&
      keySetActive(&oldSet,&look)) {
    keySetDestroy (oldSet) ;
    newSet = keySetCopy(query_undo_set);
    if (newSet)    /* keySet should check for NULL; reset Taglist */
      qbuildNewClass(class(keySet(newSet, 0)));
    keySetShow (newSet, look) ;
  }
}

/*************************************************************/

/* Initializes and creates all data structure necessary for
 * query-window interface.  Also, calls querDisplay() to pop up
 * image.
 */
void queryCreate (void)
{
  Array a ;
  QUER quer ;
  char *cp ;
  
  if(graphActivate(querGraph))
    {
      graphPop() ;
      return ;
    }
  
  querGraph = displayCreate(DtQuery) ;
  if (!querGraph)
    return ;
  
  quer=(QUER)messalloc(sizeof(struct QCTRLSTUFF));
  quer->magic = MAGIC;
  quer->graph = querGraph ; /* provision for multi windows */
  quer->last_find = -1 ;
  
  graphTextBounds (180,40) ;   /* needed for text box sizing */
  graphRegister (DESTROY, querDestroy) ;
  graphRegister (PICK, querPick) ;
  graphRetitle ("Query Commands"); /* acedb 1.9.1 prepends "DtQuery :" */
  
  if (getenv("ACEDB_AUTOSAVE")) {
    RMFLAG(NO_AUTOSV);
    querMenu[MENU_AUTO].text = "Disable AutoSave";
    }
  else {
    SETFLAG(NO_AUTOSV);
    querMenu[MENU_AUTO].text = "Enable AutoSave ";
  }

  graphMenu (querMenu) ;  
  graphAssociate (&MAGIC, quer);
  graphHelp ("Query") ;
  
  if ((cp = filName ("wquery", 0, "r")))
    strcpy (quer->dirName, cp) ;
  strcpy (quer->fileName, "") ;
  
  quer->pgm = arrayCreate(5,Array) ;
  array(quer->pgm,0,Array) = a = arrayCreate(QBUFF_MULT*BUFFER_SIZE,char) ;
  array(a,0,char) = 0 ;
  
  querDisplay();
}

/********** end of file querdisp.c *************/

