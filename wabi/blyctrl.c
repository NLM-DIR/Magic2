
/*  File: blyctrl.c
 *  Author: mieg
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Graphic interface to acembly
 * HISTORY:
 * Last edited: Nov 26 18:42 1998 (fw)
 * * Nov 19 15:00 1998 (edgrif): Fixed callback func. dec. for
 *              graphTextScrollEditor.
 * Created: April 1997 (mieg)
 *-------------------------------------------------------------------
 */

/* @(#)blyctrl.c	1.2 5/23/97 */
#include "acedb.h"
#include "graph.h"
#include "mytime.h"
#include "query.h"
#include "bs.h"
#include "systags.h"
#include "session.h"
#include "display.h"
#include "keyset.h"
#include "pick.h"


/* #efine DEBUG */


enum  MODE { UNDEFINED = 0, SHOTGUN, MUTANT, CDNA, VECTOR} ;

static Graph blyControlGraph = 0 ;

static int state = 0 ;
static int BLYMAGIC = 0 ;

#define NN 30

typedef struct lookstruct 
{
  void *magic ;
  int level ;
  int mode ;
  char *buffer ;
  KEY wildSeq ;
  int wildLength ;
  int state ; char *title ; AC_HANDLE hh ; 
  BOOL isShotgun ; mytime_t date ;  
  KEYSET seqVectors ;
  char vname[NN] ;
  char leftSite[NN], rightSite[NN] ; int nn1 ;
} *LOOK ;
static void blyctrlControlDraw (LOOK look) ;
static LOOK globalLook = 0 ;
static LOOK BLYGET (char *title) ;

/***************************************************************************************/

static void blyctrlControlDestroy (void)
{
  LOOK look = BLYGET ("blyctrlControlDestroy") ;
  look->magic = 0 ;
  handleDestroy (look->hh) ;
  keySetDestroy (look->seqVectors) ;
  blyControlGraph = 0 ;
  globalLook = 0 ;
}


/***************************************************************************************/

static void blyctrlInit (LOOK look)
{
  KEYSET ks ;
  KEY clone = 0, mode = 0  ;
  OBJ obj = 0 ;
  
  look->state = 0 ;
  look->title = "unc-32" ;
  look->date = timeNow() ;
  look->isShotgun = TRUE ;

  ks = query (0, "FIND Clone Main_clone") ;
  if (!keySetMax(ks) || !(clone = keySet(ks, 0)) || !(obj = bsCreate (clone)))
    goto abort ;
  keySetDestroy (ks) ;
  look->seqVectors = queryKey (clone, "FOLLOW Sequencing_Vector") ;
  
  if (bsFindTag (obj, str2tag("Project_type")))
    bsGetKeyTags (obj, _bsRight, &mode) ;
  if (pickMatch ("*mutation*", name(mode)))
    look->mode = CDNA ;
  else if (pickMatch ("*sequence*", name(mode)))
    look->mode = SHOTGUN ;
  else if (pickMatch ("*cdna*", name(mode)))
    look->mode = CDNA ;
  else
    look->mode = UNDEFINED ;

  bsGetData (obj, str2tag("Level"), _Int, &look->level) ;

abort:
  bsDestroy (obj) ;
  keySetDestroy (ks) ;  
}

/***************************************************************************************/
/***************************************************************************************/
/****************** Bly control graph, actual calls to bly*********************************************/
/***************************************************************************************/

static void popBlyControl(void)
{
  if (graphActivate (blyControlGraph))
    graphPop() ;
}

static LOOK BLYGET (char *title)
{
  LOOK look = globalLook ; 

  
  if (!graphActivate (blyControlGraph))
    { blyControlGraph = 0 ; globalLook = 0 ; return 0 ; }
  
  graphCheckEditors(blyControlGraph, TRUE) ;
  return look ;
}

/******************************/

static void blyctrlRedraw (void)
{
  LOOK look = BLYGET ("blyctrlRedraw") ;

  popBlyControl() ;
  blyctrlControlDraw (look) ;
}

/******************************/

static void beforeDefCompute (void)
{
  extern void defCompute(void) ;
  graphDestroy() ;
  defCompute () ;
}

/******************************/

static void beforeDefSimple (void)
{
  extern void fMapTraceReassembleAll (void) ;

  graphDestroy() ;
  fMapTraceReassembleAll () ;
}

/*****************************************************************/

static MENUOPT blyControlMenu[] =
{
  /*
 { graphDestroy, "Quit" },
  { help, "Help" }, 
  */
  { beforeDefSimple,       "Assembly: simple interface"},
  { beforeDefCompute,       "Assembly: full toolbox"}, 
  { 0, 0 }
} ;


/*****************************************************************/
#ifdef JUNK
static BOOL checkTC (int tc)
{
  if (tc < 20 || tc > 99)
    return FALSE ;
  return TRUE ;
}

static void createVector(void)
{
}

static void selectVector(void)
{
}

static void showVectors(void)
{
  KEYSET ks = query (0,"FIND Vector") ;
  keySetNewDisplay (ks,"Vectors") ;
}


static void vectorSelect(LOOK look)
{
  int line = 8 ; look-> nn1 = 40 ;

  memset (look->vname,0, NN) ;
  memset (look->leftSite,0, NN) ;
  memset (look->rightSite,0, NN) ;
  graphButton ("List Known Vectors", showVectors, 4, line) ; line += 2 ;
  graphButton ("Select a Known Vectors", selectVector, 4, line) ; line += 2 ;
  graphButton ("Create a new vector", createVector, 4, line) ; line += 2 ;
  graphWordEditor ("Vector name",look->vname,18, 4, line++,0) ;
  graphWordEditor ("upStream bases",look->leftSite,28, 4, line++,0) ;
  graphWordEditor ("downStream bases",look->rightSite,28, 4, line++,0) ;
  graphIntEditor ("Primer to insert distance", &look->nn1, 4, line++, 0) ; 
}
#endif

/***********************************************************************************/

static void sSh (LOOK look)
{
  look->mode = SHOTGUN ;
  look->level = 0 ;
}

static void scDNA (LOOK look)
{  
  look->mode = SHOTGUN ;
  look->level = 0 ;
}

static void sMuOk (void)
{
  LOOK look = BLYGET ("sMuOk") ;
  
  switch (look->level)
    {
    case 0:
      look->level = 1 ;  /* should check */
      break ;
    } 
}


static void sMuChange (void)
{
  LOOK look = BLYGET ("sMuChange") ;
  
  switch (look->level)
    {
    case 1:
      look->level = 0 ; 
      break ;
    } 
}

static BOOL dEntry (char *unused1, int unused2)
{
  return TRUE ;
}

static void acefileselect (void)
{
  LOOK look = BLYGET ("acefileselect") ;
   blyctrlControlDraw (look) ;
}


static void dnafileselect (void)
{
  LOOK look = BLYGET ("dnafileselect") ;
  blyctrlControlDraw (look) ;
}


static void dbselect (void)
{
  graphClear () ;
  graphText ("hello world", 3, 3) ;

  graphRedraw () ;
}


static void sMu (LOOK look)
{
  int line = 3 ;
  
  switch (look->level)
    {
    case 0:
      if (!look->buffer)
	look->buffer = messalloc (32001) ;
      graphText ("Please enter the wild type sequence", 2, line++) ;
      graphText ("You may:", 2, line++) ;
      graphTextScrollEditor ("enter the dna here", look->buffer, 32000, 32, 4, line, dEntry) ;
      line += 1.3 ; graphText ("or select a file containing plain dna", 4, line) ;
      graphButton (" ", dnafileselect, 60, line) ;
      line += 1.3 ; graphText ("or select a fasta file", 4, line) ;
      graphButton (" ", dnafileselect, 60, line) ;
      line += 1.3 ; graphText ("or select a .ace file", 4, line) ;
      graphButton (" ", acefileselect, 60, line) ;
      line += 1.3 ; graphText ("or select among the sequences known it this database", 4, line) ;
      graphButton (" ", dbselect, 60, line) ;
      line += 2 ; 
      graphButton ("OK", sMuOk, 8, line) ;
      break ;
    case 1:
      graphText 
	( messprintf ("The wild type sequence is called %s", name(look->wildSeq)), 
	  2, line++) ;
      graphText 
	( messprintf ("its length is %d base pairs", look->wildLength),
	  2, line++) ;
      graphButton ("Change", sMuChange, 8, line) ;
      
      /*
      
  wild ? ;
    couper coller ou lire un fichier fasta ou un .ace ;

pcr  || cloned fragments

  pcr:
nomenclature ;

    allele+3.date
  
   f122+3.30av
   e2310-7.4aug
   f111.toto.4avr
   f122+3.1.30avr


cloned:
   vector ?
   primer ?

   f122.clone.4avr 
   

SCF_dir ?


==========


  */
    default:
	break ;
    }
}


/***********************************************************************************/
/* Undefined Project */

#ifdef JUNK

static void ppM (void)
{
  LOOK look = BLYGET ("ppS") ;
  
  look->mode = MUTANT ;
  look->level = 0 ;
  blyctrlControlDraw (look) ;
}

static void ppS (void)
{
  LOOK look = BLYGET ("ppS") ;
  
  look->mode = SHOTGUN ;
  look->level = 0 ;
  blyctrlControlDraw (look) ;
}

static void ppD (void)
{
  LOOK look = BLYGET ("ppS") ;
  
  look->mode = CDNA ;
  look->level = 0 ;
  blyctrlControlDraw (look) ;
}

static MENUOPT pleaseMenu[]  = 
{
  {ppS, "Sequence a new area"},
  {ppD, "Sequence cDNA and align on genomic"},
  {ppM, "Sequence mutants and compare to wild type"},
  {0, 0}
} ;
#endif

/***********************************************************************************/

static void blyctrlControlDraw (LOOK look)
{ 
  int line = 2 ;

  graphClear () ;
  
  switch (look->mode)
    {
    case UNDEFINED:
      graphButtons (blyControlMenu, 2, line, 20) ;
      /*
      graphText ("Please define your project", 2, line++) ;
      graphButtons (pleaseMenu, 2, line, 20) ;
      */
      break ;
    case SHOTGUN:
      sSh (look) ;
      break ;
    case MUTANT:
      sMu (look) ;
      break ;
    case CDNA:
      scDNA(look) ;
      break ;
    }
  graphRedraw () ;
}

void blyControl (void)
{
  Graph old = graphActive () ;
  AC_HANDLE hh ;
  LOOK look ;

  if (graphActivate (blyControlGraph))
    { graphPop () ;
    return ;
    }
  state = 0 ;
  blyControlGraph = graphCreate (TEXT_SCROLL, "Project parameters", 
				 0, 0, 0.6, 0.5) ;
  graphHelp("Acembly.control") ;
  graphRegister (DESTROY, blyctrlControlDestroy) ;
  graphRegister (RESIZE, blyctrlRedraw) ;
  graphMenu (blyControlMenu) ;
  
  hh = handleCreate() ;
  look = (LOOK) halloc(sizeof(struct lookstruct), hh) ;
  globalLook = look ;
  look->magic = &BLYMAGIC ;
  look->hh = hh ;
  blyctrlInit(look) ;
  blyctrlControlDraw (look) ;
  graphActivate (old) ;
}

/***************************************************************************************/
/***************************************************************************************/
