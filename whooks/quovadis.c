/*  File: quovadis.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Feb  4 10:52 1999 (fw)
 * * Feb  1 21:41 1996 (rd)
 * Created: Thu Sep  8 15:42:54 1994 (mieg)
 *-------------------------------------------------------------------
 */
/* $Id: quovadis.c,v 1.14 2017/01/17 21:16:43 mieg Exp $ */

/* main menu */

#include "acedb.h"
#include "session.h"
/* #include "../wh/menu.h" do not, it may be obtained with conflict from /usr/include */
#include "systags.h"
#include "sysclass.h"  
#include "classes.h"  
#include "lex.h"
#include "orepeats.h"
#include "hseqdisp.h"
#include "orepeats.h"
#include "dump.h"
#include "bqldisp.h"

extern DumpFunc dumpFunc[] ;
extern ParseFunc parseFunc[] ;
extern KillFunc killFunc[] ;

#ifndef NON_GRAPHIC

#define DEBUG    1  /* mieg , i want acdbtest linked in */

extern Array  displayFuncArray ;
#include "display.h"

/* includes for kernel sub-modules, 
   some containing prototypes for display functions */
#include "session.h"
#include "model.h"		/* for readModels() */
#include "dna.h"		/* for dnaAnalyse() */
#include "query.h"		/* for queryCreate() */

/* includes for the display modules, which declare functions
 * that are the top-level entry functions to the mainMenu items */
#include "sessiondisp.h"	/* for sessionControl() */
#include "spreaddisp.h"		/* for spreadDispCreate() */
#include "dendrogram.h"
#include "hseqdisp.h" 
#include "querydisp.h"		/* for bqlCreate qbeCreate, qbuildCreate */
#include "bqldisp.h"

/* these externs have to be removed and replaced with includes
   that declare the prototypes for these functions */
extern void
	     exitButtonAction(void),
             gMapCompute(void) ,
	     blyControl(void),
             parseControl(void),
             dumpAll(void) ,
             oxgridCreate(void),
#ifdef DEBUG
	     acedbtest(void),
#endif /* DEBUG acedbtest */
             fp2 (void),
	     acedbstatus(void),
             chronoShow(void),
	     dumpAsnDefinitions(void),
	     alignMaps(void), 
             updateData(void),
	     addKeys(void),
/*	     metaCheckCode(void),	*/
/* 	     setUpdate (void),  detlef's test code */
             editPreferences(void);

/************************************************************/

static void menuSaveAction (void)	/* menu function */
{
  if (isWriteAccess())
    { BOOL reGet = messQuery("Shall I retain Write Access after saving?");
      sessionClose (TRUE) ;
      if (reGet)
	sessionGainWriteAccess();
    }
  else
    sessionGainWriteAccess () ; /* try to grab it */
} /* menuSaveAction */


static void menuWriteAccess (void)
{
  /* can't use this function directly in menu definition,
     it returns a boolean success value (ignored here) */
  sessionGainWriteAccess();
} /* menuWriteAccess */

/************************************************************/

static void qspreadDispCreate (void)
{ spreadDispCreate (0) ; }

/************************************************************/

static void menuReadModels (void)
{
  /* can't use this function directly in menu definition,
     it returns a boolean success value (ignored here) */
  readModels();
} /* menuReadModels */

/************************************************************/

MENUOPT quovadis[] = {
            { exitButtonAction, "Quit"},
            { help  ,           "Help"},
	    { graphCleanUp,     "Clean up"},
            { acedbstatus,      "Program status"},
#ifdef ACEMBLY
	    { blyControl, "Acembly"},
/* 	    { fp2, "histo genethon clones"}, */
#else
            { updateData,       "Add Update File"},
	    { editPreferences,  "Preferences"},
	    { sessionControl,   "Session control"},

#endif
            { menuSpacer,            ""},
#if defined(WIN32)
	    { 0,		"Query"},
#endif
            { queryCreate,      "Query"},   
#if !defined(MACINTOSH)
            { qbeCreate,        "Query by Examples"},   
#endif
            { qbuildCreate,     "Query builder"}, 
            { bqldCreate,        "BQL table queries"},   
  
            { qspreadDispCreate, "Table Maker"},
#if defined(WIN32)
	    { 0,  	  0}, /* submenu trailer */
#endif
            { dnaAnalyse,       "DNA Analysis"},
            { menuSpacer,       ""},
#ifndef ACEMBLY
#if defined(WIN32)
	/* (rbrusk) July 17, 1998 - Dendrogram functions */
	{ 0,		     "Phylogenetic Trees"},
#endif
	{ readTaxonomicTree, "Read Taxonomy File"},   
	{ readDNATree,       "Read DNA Tree File"},   
	{ readProteinTree,   "Read Protein Tree File"},   
#if defined(WIN32)
	{ 0,  	  0}, /* submenu trailer */
#endif
 /* option below added by jld */
            { oxgridCreate,    "Comparative Maps"},
#endif /* !ACEMBLY */
            { menuSpacer,       ""},
            { addKeys,          "Add-Alias-Rename"},
            { parseControl,     "Read .ace files"},
            { alignMaps,        "Align maps"},   
            { menuReadModels,   "Read models"},
            { dumpAll,          "Dump"}, 

	    { sessionControl,   "Session control"},
 
	    { menuSpacer,            ""},
            { chronoShow,       "Chronometer"}, 

#ifdef DEBUG
	    { menuSpacer,            ""},
            { acedbtest,        "Test subroutine"},
#endif /* DEBUG acedbtest */

/*	    { metaCheckCode,    "MetaData"}, */
            { menuSpacer  ,          ""},
            { menuWriteAccess,	     "Write Access"} ,
            { menuSpacer  ,          "" },
            { menuSaveAction,        "Save" },
            { 0, 0 }		/* end of menu marker */
};

      /* For each display type, you must provide
	 a static display function and register it
	 by a line of code in the following function,

	 you must also register a special function
	 to dump and parse unprotected A type objects,
	 B types are handled by the acedb kernel.

	 The prototypes are in display.h
      */

/* includes that provide the prototypes for the DisplayFunc's */

#include "tree.h"
#include "forest.h"
#include "vmap.h"
#include "fmap.h"


/* some display modules don't have public includes 
   for the DisplayFunc yet */
extern BOOL
	keySetDisplay (KEY, KEY, BOOL),
        longTextDisplay (KEY, KEY, BOOL),	
	pMapDisplay (KEY, KEY, BOOL),
	gMapDisplay (KEY, KEY, BOOL),
	multiMapDisplay (KEY, KEY, BOOL),
	cFicheDisplay (KEY, KEY, BOOL),
        pepDisplay (KEY, KEY, BOOL),
        alignDisplay (KEY, KEY, BOOL),
	gridDisplay (KEY, KEY, BOOL),
	actionDisplay (KEY, KEY, BOOL),
        fingerPrintDisplay (KEY, KEY, BOOL),
	multiTraceDisplay (KEY, KEY, BOOL),
	metabDisplay (KEY, KEY, BOOL),
        cMapDisplay (KEY, KEY, BOOL),
	displayAs (KEY, KEY, BOOL),
	wwwDisplay (KEY, KEY, BOOL),
#if defined(WIN32)
	scriptDisplay(KEY, KEY, BOOL),
#endif
	dendrogramDisplay(KEY, KEY, BOOL),
        imageDisplay (KEY, KEY, BOOL),
        geneExpDisplay (KEY, KEY, BOOL) ;

/************************************************************/

static void hook(char *displayName, DisplayFunc f)
{ KEY displayKey = 0 ;	
	
  lexaddkey (displayName, &displayKey, _VDisplay) ;
  array (displayFuncArray , KEYKEY (displayKey), DisplayFunc) = f ;
}

BOOL (*mhmpDrawGifp) (KEYSET ks) = 0 ; /* used in gifacemain */
static void displayInit (void)
{ 
#ifdef ACEMBLY
  extern BOOL (*mhmpDisplayp) (KEYSET ks, int glyph, KEY colour) ; /* in ksetdisp.c */
  BOOL mhmpDisplay (KEYSET ks, int glyph, KEY colour) ;
  BOOL mhmpDrawGif (KEYSET ks) ;
#endif
  displayFuncArray = arrayCreate (32,DisplayFunc) ;

                       /*  System hooks */
  hook (TREE, treeDisplay) ;
  hook (DtKeySet, keySetDisplay) ;
  hook ("DtForest", forestDisplay) ;
  hook (DtLongText, longTextDisplay) ;
  /*   hook ("DtTableResultDisplay", spreadTableDisplay) ; ridiculous, wrong prototype */
  hook (DISPLAY_AS, displayAs) ;
  hook (WWW, wwwDisplay) ;

#if defined(WIN32)
  hook (SCRIPT, scriptDisplay) ;
#endif

                      /* Application hooks */
  hook (PMAP, pMapDisplay) ;
  hook (GMAP, gMapDisplay) ;
  hook (VMAP, vMapDisplay) ;

  hook (DtMULTIMAP, multiMapDisplay) ;
  hook (DtPmapFingerprint, fingerPrintDisplay) ;

#ifdef ACEMBLY
  mhmpDisplayp = mhmpDisplay ;
  mhmpDrawGifp = mhmpDrawGif ;
  hook (DtMultiTrace, multiTraceDisplay) ;
  hook (CFICHE, cFicheDisplay) ;
  hook (DtHSEQ, hseqDisplay) ;
  hook (DtTiling, htileDisplay) ;
  hook (DtGLOC, glocDisplay) ;
  hook (DtGLOCBIG, glocBigDisplay) ;
  hook (DtGeneExp, geneExpDisplay) ;
#endif

  hook (DtBqlDisp, bqlDisplay) ;

  hook (FMAP, fMapDisplay) ;
  hook (PEPMAP, pepDisplay) ;
  hook (DtAlignment, alignDisplay) ;
  hook (GRID, gridDisplay) ;
  hook (CMAP, cMapDisplay) ;
  hook (DtImage, imageDisplay) ;

  hook (METAB, metabDisplay) ;
  hook (DtDendrogram, dendrogramDisplay) ;
}

#endif /* NON_GRAPHIC */

/*************** A class parsing functions ****************/

extern BOOL keySetDumpFunc (FILE* f, Stack s, KEY k) ;
extern BOOL keySetParse (int level, KEY) ;

extern BOOL longTextDump (FILE* f, Stack s, KEY k) ;
extern BOOL longTextParse (int level, KEY) ;

extern BOOL tableDump (FILE* f, Stack s, KEY k) ;
extern BOOL tableParse (int level, KEY) ;

extern BOOL dnaDump (FILE* f, Stack s, KEY k) ;
extern BOOL dnaParse (int level, KEY key) ;

extern BOOL peptideDump (FILE* f, Stack s, KEY k) ;
extern BOOL peptideParse (int level, KEY key) ;

extern BOOL baseQualityDump (FILE* f, Stack s, KEY k) ;
extern BOOL baseQualityParse (int level, KEY key) ;

extern BOOL basePositionDump (FILE* f, Stack s, KEY k) ;
extern BOOL basePositionParse (int level, KEY key) ;

extern BOOL binaryDump (FILE* f, Stack s, KEY k) ;
extern BOOL binaryParseParse (int level, KEY key) ;


static void parseArrayInit (void)
{
  parseFunc[_VKeySet]  = keySetParse ;
  dumpFunc[_VKeySet]   = keySetDumpFunc ;

  parseFunc[_VDNA]  = dnaParse ;
  dumpFunc[_VDNA]   = dnaDump ;

  parseFunc[_VPeptide]  = peptideParse ;
  dumpFunc[_VPeptide]   = peptideDump ;

  parseFunc[_VBinary]  = aceBinaryParse ;
  dumpFunc[_VBinary]   = aceBinaryDump ;

  parseFunc[_VLongText]  = longTextParse ;
  dumpFunc[_VLongText]   = longTextDump ;

  parseFunc[_VTableResult]  = tableParse ;
  dumpFunc[_VTableResult]   = tableDump ;

  parseFunc[_VBaseQuality]  = baseQualityParse ;
  dumpFunc[_VBaseQuality]   = baseQualityDump ;

  parseFunc[_VBasePosition]  = basePositionParse ;
  dumpFunc[_VBasePosition]   = basePositionDump ;

  parseFunc[_VOligoRepeat]  = oligoRepeatParse ;
  dumpFunc[_VOligoRepeat]   = oligoRepeatDump ;

  parseFunc[_VFicheView]  = longTextParse ;
  dumpFunc[_VFicheView]   = longTextDump ;
}

/*********************************************************/

void qvInit (void)		/* entry point called from sessionInit */
{
#ifdef ACEMBLY
  extern void  acemblyInit (void) ;
 acemblyInit () ;
#endif

 parseArrayInit () ;

#ifndef NON_GRAPHIC
  displayInit () ;
#endif
}

/******************** end of file ************************/
 
 
 
 
 
 
