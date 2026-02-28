/*  File: bqldisplay.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: display bql tables using the table.h package
 *    a modernisation of the table system using bql syntax
 *    the visible/hidden box contriols the select bql syntax
 *    mandatory is implied ?
 *    The tab key allows to prompt for the spelling
 *
 * wspec/displays.wrm needs the line
 *   _DDtBql   -g TEXT_FULL_SCROLL -t "BQL table queries" -w .90  -height 0.75 -help Bql
 * Exported functions:
 * HISTORY:
 * Last edited: 
 * * Mar 29 2024 (mieg)
 * Created: Thu Mar 29 2024 (mieg)
 *-------------------------------------------------------------------
 */

#include "acedb.h"
#include "ac.h"
#include "bql.h"

#include "display.h"
#include "bqldisp.h"

magic_t GRAPH2BQLD_ASSOC = "BQLDISPLAY" ;
#define BQLD_MAGIC  716512935 
typedef struct bqldStruct {
  AC_HANDLE h ;
  Graph graph ;
  int magic ;
  int line ;
  BOOL syntaxHelp ;
  Array vars ;
} *BQLD ;

typedef struct bqlVarStruct {
  char buf[128] ; 
  int type ;
  int box ;
} BqlVar ;
static BQLD currentBqlDisp (char *caller) ;
static void bqldDraw (BQLD bqld) ;

#ifdef JUNK

#include "lex.h"
#include "pick.h"
#include "bs.h"
#include "dna.h"
#include "peptide.h"
#include "systags.h"
#include "sysclass.h"
#include "query.h"
#include "session.h"
#include "tree.h"

#include <ctype.h>
static void bqldChooseTag(BQLD bqld, int box, COL *c, COL *fromC, int cnt) ;
static void bqldChooseDna(BQLD bqld, int box, COL *c, COL *fromC, int cnt) ;
static void bqldChoosePep(BQLD bqld, int box, COL *c, COL *fromC, int cnt) ;

#define graphBoxBox(_box) { \
	       float _x1, _y1, _x2, _y2 ; \
	       graphBoxDim (_box, &_x1, &_y1, &_x2, &_y2) ; \
	       graphRectangle (_x1 - 0.4, _y1 - 0.1, _x2 + 0.4, _y2 + 0.1) ; \
		}	   

static FREEOPT *classeListe = 0 ;
static Array classeListeArray = 0 ;
static void classeListeInit (void) ;
static BOOL isCheckEditor = FALSE ; /* mhmp 16.04.98 */

static FREEOPT extendChoice[] =
{ 
  {4, "Define"},
  {'f', "From"},
  {'r', "Right_of"},
  {'c', "Copy"},
  {'p', "Compute"}
} ;

static FREEOPT typeChoice2[] =
{ 
  {12, "Class ..."},
  {'x', "Show Data"},
  {'b', "Show Tag"},
  {'n', "Show Next Tag"},
  {'K', "Show Next Key"},
  {SHOW_MULTI, "Data;Data;Data"},
  {SHOW_MIN, "Compute Min"},
  {SHOW_MAX, "Compute Max"},
  {SHOW_AVG, "Compute Average"},
  {SHOW_VAR, "Compute Variance"},
  {SHOW_SUM, "Compute Sum"}, 
  {SHOW_COMPUTE, "Compute equation"}, 
  {'c', "Count"}
} ;

/*****************************************************/
/******* Colonne Definition Interface ****************/ 

static void bqldDoInitColonne (BQLD bqld, int newColonne)
     /* private within bqldDisp-package */
{
  COL *c ;
  Array t ;

  c = arrayp(bqld->colonnes, newColonne, COL) ;  
  c->tagStack = 0 ;/*mhmp 17.04.98*/
  /*strcpy (c->subtitleBuffer, "") ;mhmp 16.04.98 + 20.11.02*/
  strcpy (c->subtitleBuffer, messprintf("Column #%d", newColonne + 1)) ;
  strcpy (c->legendBuffer, "") ;
  c->type = 0 ; c->realType = 0 ;
  c->colonne = newColonne ;
  c->hidden = FALSE ;
  c->extend = 'f' ;
  c->extendp = freekey2text(c->extend, extendChoice) ;
  c->mandatory = 1 ;
  c->showType = SHOW_ALL ;
  c->nonLocal = FALSE ;
  c->from = 1 ;

  if (arrayExists(t = bqld->tableau))
    {
      int i = arrayMax(t) ;
      int max = arrayMax (bqld->colonnes) ;

      while (i--)
	array(array(t,i,Array), max, SPCELL).u.k = 0 ;
    }
  c->width = 12 ;
  *c->conditionBuffer = 0 ;
  bqld->modified = TRUE ;

  return;
} /* bqldDoInitColonne */


void bqldInitColonne (BQLD bqld)
     /* private within bqldDisp-package */
{
  COL *c ;
  int max ;

  max = arrayMax(bqld->colonnes) ;
  array(bqld->pos2col,max, int) = max ;
  bqld->activeColonne = max ;

  c = arrayp(bqld->colonnes, bqld->activeColonne, COL) ;  
  if (c) {} ; /* for compiler happiness */
  bqldDoInitColonne (bqld, max) ;

  return;
} /* bqldInitColonne */

/*****************************************************/

static COL *bqldDefineActivate(int box)
{
  COL *c ; 
  BQLD bqld = currentBqlDisp("bqldDefineActivate") ; 

  bqld->activeColonne = bqld->defBoxPerCol ? 
      (box - bqld->definitionBox -1) / bqld->defBoxPerCol : 0 ;
  c = arrayp(bqld->colonnes, bqld->activeColonne, COL) ;  
  if (bqld->activeDefBox) 
    graphBoxDraw(bqld->activeDefBox, BLACK, WHITE) ; 
  bqld->activeDefBox = box ;

  return c ;
} /* bqldDefineActivate */

/******** Optional ********/

static FREEOPT optionalChoice[] =
{ 
  {3, "Define"},
  {1, "Optional"},
  {2, "Mandatory"},
   {0, "Null"}
} ;

static void bqldDefineOptional (KEY k, int box)
{
  COL *c = bqldDefineActivate(box) ;
  KEY old = c->mandatory ;
  BQLD bqld = currentBqlDisp("defineoptional") ;

  c->mandatory = k ;
  c->optionalp = freekey2text(k, optionalChoice) ;
  graphBoxDraw(box, BLACK, LIGHTBLUE) ;
  if (k != old)
    {
      bqld->modified = TRUE ;
      bqld->modif = TRUE ; /*30.11.98 */
    }
}

/******** Extend ********/

static void bqldDefineExtend (KEY k, int box)
{
  COL *c, *fromC ;
  KEY old ;
  BQLD bqld = currentBqlDisp("defineextend") ;

  /* mhmp 11.01.98 suite a un crash (pas reproduit) */
  if (box <= 0) 
    {
      messout ("Warning: bqldDefineExtend called with box <= 0") ;
      return ;
    }
  c = bqldDefineActivate(box) ; 
  old = c->extend ;
  c->extend = k ;
  c->extendp = freekey2text(k, extendChoice) ;
  graphBoxDraw(box, BLACK, LIGHTBLUE) ;
  if (k != old)
    {
      bqld->modified = TRUE ;
      bqld->modif = TRUE ; /*30.11.98 */ 
      isCheckEditor = TRUE ;
      
      switch (c->extend)
	{
	case 'f': /* from */ 
	  c->from = 1 ; 
	  c->showType = SHOW_ALL ; 
	  c->nonLocal = FALSE ;
	  c->type = 0 ;
	  c->classe = 0 ;
	  c->realType = 0 ;
	  break ;
	case 'r': /* right_of */
	  c->from = c->colonne ;
	  c->showType = SHOW_ALL ; 
	  c->nonLocal = FALSE ;
	  c->type = 0 ;
	  c->classe = 0 ;
	  c->realType = 0 ;
	  break ;
	case 'c': /* copy */
	  c->from = 1 ; 
	  c->showType = SHOW_COPY ; 
	  c->nonLocal = FALSE ;
	  fromC =  arrp(bqld->colonnes, c->from - 1 , COL) ;
	  c->type = fromC->type ;
	  c->classe = fromC->classe ;
	  c->realType = c->type ; 
	  break ;
	case 'p': /* compute */
	  c->from = 1 ; 
	  c->showType = SHOW_COMPUTE ;
	  c->showtypep = freekey2text(k, typeChoice2) ; 
	  c->nonLocal = TRUE ;  
	  fromC =  arrp(bqld->colonnes, c->from - 1 , COL) ;
	  c->type = 'f' ;
	  c->classe = 0 ;
	  c->realType = 'f' ;
	  break ;
	}      
      stackDestroy (c->tagStack) ;
      c->tagStack = stackCreate (50) ;
    }

  bqldDefineColonne(TRUE) ;
} /* bqldDefineExtend */

/******** Hide ********/

static FREEOPT hideChoice[] = { 
  {2, "Hide"},
  {TRUE, "Hidden"},
  {FALSE, "Visible"}
} ;

static void bqldDefineHide (KEY k, int box)
{
  COL *c = bqldDefineActivate(box) ;
  KEY old = c->hidden ;
  BQLD bqld = currentBqlDisp("definehide") ;

  c->hidden = k ;
  c->hiddenp = freekey2text(k, hideChoice) ;
  graphBoxDraw(box, BLACK, LIGHTBLUE) ;
  if (k != old)
    {
      bqld->modified = TRUE ; 
      bqld->modif = TRUE ; /*30.11.98 */
    }
} /* bqldDefineHide */

/******** Type ********/

static FREEOPT typeChoice[] =
{ 
   {9, "Class ..."},
   {'k', "Class"},
   {'i', "Integer"},
   {'f', "Float"},
   {'d', "Date"},
   {'t', "Text"},
   {'D', "DNA"},
   {'P', "Peptide"},
   {'b', "Show_Tag"},
   {'n', "Next_Tag"},
   {'K', "Next_Key"},
   {'5', "Count"}
} ;

static void bqldDefineType (KEY k, int box)
{
  KEY old ;
  COL *c = bqldDefineActivate(box) ;
  BQLD bqld = currentBqlDisp("defineType") ;

  old = c->type ;
  switch (k)
    {
    case 'b': 
      c->nonLocal = FALSE ;
      break ;
    case 'x': 
      c->nonLocal = FALSE ;
      c->showType = SHOW_ALL ;
      c->type = c->realType ;
      c->showtypep = freekey2text(k, typeChoice2) ;
      goto ok ;
      break ;
    case 'n':
      c->classe = _VSystem ;
      c->nonLocal = FALSE ;
      break ;
    case 'K':
	{ KEY key ;
	classeListeInit () ;
	if (graphSelect (&key, classeListe))
	{ c->classe = key ;
	}
	else
	c->classe = 0 ;
	}
      c->nonLocal = FALSE ;
      break ;	
    case SHOW_MIN:
    case SHOW_MAX:	
    case SHOW_AVG:	
    case SHOW_VAR:	
    case SHOW_SUM:
    case SHOW_COMPUTE:
      switch (c->realType)
	{
	case 'i': case 'f': case 'd':
          c->showType = k ; c->type = c->realType ;
	  bqld->modified = TRUE ;
	  bqld->modif = TRUE ; /*mhmp 30.11.98 */
	  c->nonLocal = TRUE ;  
	  c->showtypep = freekey2text(k, typeChoice2) ;
	  goto ok ;
	  break ;
	default:
	  return ;
	}
    case SHOW_MULTI:
      c->showType = k ; 
      c->type = c->realType ;
      bqld->modified = TRUE ;
      bqld->modif = TRUE ; /*mhmp 30.11.98 */
      c->nonLocal = TRUE ;  
      c->showtypep = freekey2text(k, typeChoice2) ;
      {
	int i1 ;
	COL * c1;
	for (i1 = bqld->activeColonne + 1 ; i1 < arrayMax(bqld->colonnes) ; i1++)
	  { 
	    c1 = arrayp(bqld->colonnes, i1, COL) ;  
	    if (c1->from == bqld->activeColonne + 1)
	      messout ("You should modify colonne %d, since you cannot derive from a multi valued column", i1 + 1) ;
	  }
      }
      goto ok ;
      break ;
    case 'c':
      c->nonLocal = TRUE ;
      break ;
    }
  c->showType = SHOW_ALL ;
  c->type = k ;
  c->typep = freekey2text(c->type, typeChoice2) ;
ok:
  if (old != c->type)
    {
      bqld->modified = TRUE ;
      bqld->modif = TRUE ; /* mhmp 30.11.98 */
    }

  graphBoxDraw(box, BLACK, LIGHTBLUE) ;
  if (1) bqldDefineColonne(TRUE) ;
}

/********* Pick boxes ************/

static void bqldPickDefinition(int box, double x, double y)
{
  COL *c ; int i ;
  BQLD bqld = currentBqlDisp("bqldPickDefinition") ;

  if (!graphCheckEditors (graphActive(), 0))
    return ;

  if (!box)
    return ;
  bqld->modif = TRUE ; /*mhmp 25.11.98 */
  if (box == bqld->titleBox)
    { graphTextEntry(bqld->titleBuffer,0,0,0,0) ;
      return ;
    }
  if (box == bqld->paramBox)
    { graphTextEntry(bqld->paramBuffer,0,0,0,0) ;
      return ;
    }
  if (box == bqld->sortColonneBox) return ;
   
  i = bqld->activeColonne = (box - bqld->definitionBox -1) / bqld->defBoxPerCol ; 
  if (i < 0 || i >= arrayMax(bqld->colonnes))
    return ;
  c = arrayp(bqld->colonnes, bqld->activeColonne, COL) ;  

  if (!( (box - bqld->widthBox) % bqld->defBoxPerCol))
    /*  graphTextEntry(c->widthBuffer,0,0,0,0) ;*/return ;
  else if (!( (box - bqld->fromBox) % bqld->defBoxPerCol) &&
	   bqld->activeColonne)
    /*   graphTextEntry(c->fromBuffer,0,0,0,0) ;*/ return  ;
  else if (!( (box - bqld->conditionBox) % bqld->defBoxPerCol))
    graphTextScrollEntry(c->conditionBuffer,0,0,0,0,0) ;
  else if (!( (box - bqld->subTitleBox) % bqld->defBoxPerCol))
    graphTextEntry(c->subtitleBuffer,0,0,0,0) ;
  else if (!( (box - bqld->legendBox) % bqld->defBoxPerCol))
    graphTextEntry(c->legendBuffer,0,0,0,0) ;
  /* mhmp 12.07.02 peins-le en vert */
  else if (!( (box - bqld->dna1Box) % bqld->defBoxPerCol))
    graphTextScrollEntry(c->dna1Buffer,0,0,0,0,0) ;
  else if (!( (box - bqld->dna2Box) % bqld->defBoxPerCol))
    graphTextScrollEntry(c->dna2Buffer,0,0,0,0,0) ;
  else if (!( (box - bqld->typeBox) % bqld->defBoxPerCol))
    { KEY key = c->type ;
      if (bqld->activeColonne)
 	{ if (graphSelect (&key, typeChoice2))
	    bqldDefineType (key, box) ;
	}
      else
 	{ if (graphSelect (&key, classeListe))
	    { c->type = c->realType = 'k' ;
	      c->classe = key ;
	      c->nonLocal = FALSE ;
	      c->typep = freekey2text(c->classe, classeListe) ;
	    }
	}
      graphBoxDraw(box, BLACK, WHITE) ;
      bqldDefineColonne(TRUE) ;/* added by mhmp 01.07.98, seems useless, mieg 
is very useful mhmp */	
     }
  else if (!( (box - bqld->dnaBox) % bqld->defBoxPerCol))
    bqldChooseDna(bqld, box, c, arrp(bqld->colonnes, c->from - 1 , COL), c->extend) ;
  else if (!( (box - bqld->pepBox) % bqld->defBoxPerCol))
    bqldChoosePep(bqld, box, c, arrp(bqld->colonnes, c->from - 1 , COL), c->extend) ;
  else if (!( (box - bqld->tagBox) % bqld->defBoxPerCol))
    bqldChooseTag(bqld, box, c, arrp(bqld->colonnes, c->from - 1 , COL), c->extend) ;
  else if (!( (box - bqld->tagTextBox) % bqld->defBoxPerCol))
    graphTextScrollEntry(c->tagTextBuffer,0,0,0,0,0) ;
  else if (!( (box - bqld->hideBox) % bqld->defBoxPerCol))
    { 
      c->hidden++ ;
      c->hidden = c->hidden % hideChoice->key ;
      c->hiddenp = freekey2text(c->hidden, hideChoice) ;
      graphBoxDraw(box, BLACK, LIGHTBLUE) ;
    }
  else if (!( (box - bqld->optionalBox) % bqld->defBoxPerCol))
    { 
      c->mandatory++ ;
      c->mandatory = c->mandatory % optionalChoice->key ;
      c->optionalp = freekey2text(c->mandatory, optionalChoice) ;
      bqld->modified = TRUE ;
      graphBoxDraw(box, BLACK, LIGHTBLUE) ;
    }
  else if (!( (box - bqld->extendBox) % bqld->defBoxPerCol))
    { 
      int extend = 'f' ;

      switch (c->extend)
	{
	case 'f': /* from */     extend = 'r' ; break ;
	case 'r': /* right_of */ extend = 'c' ; break ;
	case 'c': /* copy */     extend = 'f' ; break ;
	}
      c->extend = 0 ; /* insures a redraw */
      bqldDefineExtend (extend, bqld->extendBox) ; 
    }
  else
    { if (bqld->activeDefBox) 
	graphBoxDraw(bqld->activeDefBox, BLACK, WHITE) ; 
      bqld->activeDefBox = box ;
    }

  return;
} /* bqldPickDefinition */

/********* DNA/Pep Chooser 1 *******************/

static void bqldChooseDna(BQLD bqld, int box,
			    COL *c, COL *fromC, int continuation)
{
  COL *c1 = c ;

  while (c1->extend != 'f' && fromC->from > 0  && fromC->from < arrayMax(bqld->colonnes))
    { c1 = fromC ;
      fromC = arrayp(bqld->colonnes, fromC->from - 1, COL) ;   
    }
  if (fromC->type == 'K' && !fromC->classe)
    { messout("First define the class of column %d",c->from) ;
      return ;
    }
  if (ace_lower(fromC->type) != 'k' || !fromC->classe)
    { messout("First define the colonne you construct from") ;
      return ;
    }
  if (ace_lower(fromC->showType) == SHOW_MULTI)
    { messout("You cannot construct onto a multi-valued column") ;
      return ;
    }
  if (!dnaInClass (superClass(fromC->classe)))
    { messout("You cannot select DNA from class: %s",
	      name(fromC->classe)) ;
      return ;
    } 
  c->nonLocal = FALSE ;
  c->showType = SHOW_DNA ;
  c->showtypep = freekey2text('D', typeChoice2) ;
  c->type = 'D' ; 
  c->dna1= 1 ; c->dna2 = 20 ;
  memset (c->dna1Buffer, 0, 256) ;
  memset (c->dna2Buffer, 0, 256) ;
  strcpy (c->dna1Buffer, "1") ;
  strcpy (c->dna2Buffer, "20") ;
  bqld->modified = TRUE ;
  bqldDefineColonne(TRUE) ;

  return;
} /* bqldChooseDna */

/********* DNA/Pep Chooser 1 *******************/

static void bqldChoosePep(BQLD bqld, int box,
			    COL *c, COL *fromC, int continuation)
{
  COL *c1 = c ;

  while (c1->extend != 'f' && fromC->from > 0  && fromC->from < arrayMax(bqld->colonnes))
    { c1 = fromC ;
      fromC = arrayp(bqld->colonnes, fromC->from - 1, COL) ;   
    }
  if (fromC->type == 'K' && !fromC->classe)
    { messout("First define the class of column %d",c->from) ;
      return ;
    }
  if (ace_lower(fromC->type) != 'k' || !fromC->classe)
    { messout("First define the colonne you construct from") ;
      return ;
    }
  if (ace_lower(fromC->showType) == SHOW_MULTI)
    { messout("You cannot construct onto a multi-valued column") ;
      return ;
    }
  if (!pepInClass (superClass(fromC->classe)))
    { messout("You cannot select Peptide from class: %s",
	      name(fromC->classe)) ;
      return ;
    } 
  c->nonLocal = FALSE ;
  c->showType = SHOW_PEP ;
  c->showtypep = freekey2text('P', typeChoice2) ;
  c->type = 'P' ; 
  c->dna1= 1 ; c->dna2 = 20 ;
  memset (c->dna1Buffer, 0, 256) ;
  memset (c->dna2Buffer, 0, 256) ;
  strcpy (c->dna1Buffer, "1") ;
  strcpy (c->dna2Buffer, "20") ;
  bqld->modified = TRUE ;
  bqldDefineColonne(TRUE) ;

  return;
} /* bqldChoosePep */

/*********Tag Chooser 1 *******************/

static void bqldChooseTag(BQLD bqld, int box,
			    COL *c, COL *fromC, int continuation)
{
  int classe ;
  COL *c1 = c ;
  Stack s1 = 0 ; KEY tag1 = 0, tag2 ; int type1 ;

  while (c1->extend && c1->extend != 'f' && fromC->from > 0  && fromC->from < arrayMax(bqld->colonnes))
    { c1 = fromC ;
      if (class(fromC->tag) || (fromC->tag && fromC->tag < _Date ))
	tag1 = 0 ;
      else
	{
	  tag1 = fromC->tag ;
	  if (lexword2key (c1->tagp, &tag2, _VSystem))
	    tag1 = tag2 ;
	}
      fromC = arrayp(bqld->colonnes, fromC->from - 1, COL) ; 
    }
  if (fromC->type == 'K' && !fromC->classe)
    { messout("First define the class of column %d",c->from) ;
      return ;
    }
  if (ace_lower(fromC->type) != 'k' || !fromC->classe)
    { messout("First define the colonne you construct from") ;
      return ;
    }
  if (ace_lower(fromC->showType) == SHOW_MULTI)
    { messout("You cannot construct onto a multi-valued column") ;
      return ;
    }
  if (pickList[superClass(fromC->classe)].type != 'B')
    { messout("You cannot select a tag in non B-class: %s",
	      name(fromC->classe)) ;
      return ;
    }
  /* if user cancels, do not touch anything in c->  */
  s1 = stackCreate (32) ; 
  type1 = c->type ; 
  if (!treeChooseTagFromModel(&type1, &classe, fromC->classe, &tag1, s1
			      , continuation == 'r' ? 2 : 1))
    return ;
  c->tag = tag1 ; c->type = type1 ;
  stackDestroy (c->tagStack) ;
  c->tagStack = s1 ;
  c->nonLocal = FALSE ; c->showType = SHOW_ALL ; c->realType = c->type ; 
  lexword2key(pickClass2Word(classe), &c->classe, _VClass) ;
  c->tagp = stackText(c->tagStack, 0) ;
  bqld->modified = TRUE ;
  bqldDefineColonne(FALSE) ;

  return;
} /* bqldChooseTag */

/*********Tag Chooser 2 *******************/

static void bqldEditTagText(void)
{
  COL *c ;
  char *cp, *cq ;
  BQLD bqld = currentBqlDisp("bqldEditTagText") ;

  c = arrayp(bqld->colonnes, bqld->activeColonne, COL) ;  
  cp = c->tagTextBuffer ;
  if (!*cp)
    strcpy(cp, "?")  ;
  c->tagStack = stackReCreate(c->tagStack, 32) ;
  cq = cp + strlen(cp) -1 ;
  while (cq >= cp && *cq == ' ')
    *cq-- = 0 ;
  pushText(c->tagStack, cp) ;
  c->tagp = stackText(c->tagStack, 0) ;
  bqld->modified = TRUE ;
  bqldDefineColonne(FALSE) ;

  return;
} /* bqldEditTagText */

/******** TagChooser 3 *****  Less good than the other method

static void bqldDefineTag (KEY k, int box)
{
  COL *c = bqldDefineActivate(box) ;
  KEY old = c->tag ;
  BQLD bqld = currentBqlDisp("bqldDefineTag") ;

  c->tag = k ;
  c->tagp = freekey2text(k, arrp(c->tagMenu, 0, FREEOPT)) ;
  c->type = 'b' ;
  c->typep =  freekey2text(c->type, typeChoice) ;
  stackDestroy(c->tagStack) ;
  graphBoxDraw(box, BLACK, LIGHTBLUE) ;
  if (k != old)
    bqld->modified = TRUE ;
}



static Array bqldDefineTagMenu (int classe)
{ KEYSET ks ; Array a = 0 ; int i ;
    
  ks = bsTagsInClass(classe) ;
  if (ks)
    { a = arrayCreate(12, FREEOPT) ;
      for (i=0; i<keySetMax(ks); i++)
	{ array(a, i+1, FREEOPT).key = keySet(ks,i) ;
	  array(a, i+1, FREEOPT).text = name(keySet(ks,i)) ;
	}
      
      array(a, 0, FREEOPT).key = keySetMax(ks) ;
      array(a, 0, FREEOPT).text = "?" ;
      keySetDestroy(ks) ;
    }
  return a ;
}
  
*************************************/

/******** Class ********/
static void classeListeInit ()
{
  KEYSET ks ;
  KEY key, model, table ;
  int i, j ;
  unsigned char mask ;
  
  if (classeListe)
    return ;
  
  ks = queryLocalParametrized (0, "FIND Class","") ;
  arraySort (ks, keySetAlphaOrder) ;
  classeListeArray = arrayCreate (keySetMax(ks) + 1 , FREEOPT) ;

  j = 1 ;
  for (i=0 ; i < arrayMax(ks) ; i++)
    { key = table = keySet(ks, i) ;
      pickIsA (&table, &mask) ;
      if (pickList [table].type == 'B' && 
	  (model = pickList [table].model) &&
	  iskey(model) == 2)
	    { arrayp(classeListeArray, j, FREEOPT)->key = key ;
	      arrayp(classeListeArray, j, FREEOPT)->text = name(key) ;
	      j++ ;
	    }
    }

  keySetDestroy (ks) ;
  arrayp(classeListeArray, 0, FREEOPT)->key = j - 1 ;
  arrayp(classeListeArray, 0, FREEOPT)->text = "Classes" ;
  
  classeListe = arrp(classeListeArray, 0, FREEOPT) ;

  return;
} /* classeListeInit */

static void bqldDefineClass (KEY k, int box)
{
  COL *c = bqldDefineActivate(box) ;
  int old = c->classe ; 
  BQLD bqld = currentBqlDisp("bqldDefineClass") ;

  c->type = c->realType = 'k' ;
  c->nonLocal = FALSE ;
  c->classe = k ;
  c->classp = name(k) ;
  
  bqldDefineColonne(FALSE) ;
  if (k != old)
    bqld->modified = TRUE ;

  return;
} /* bqldDefineClass */

/************************************************/

static void bqldDefineFrom (void)
{ 
  int level ;
  COL *col ;
  int old ;
  BQLD bqld = currentBqlDisp("bqldDefineFrom") ;
  
  col = arrayp(bqld->colonnes, bqld->activeColonne, COL) ;  
  old = col->from ;
  if (!bqld->activeColonne)
    return ;
  
  level= freesettext(col->fromBuffer, "") ;
  if (! freecard(level) ||
      !freecheck("i"))
    { messout("Please type an integer") ;
    sprintf(col->fromBuffer,"%d", col->from) ;
    graphTextEntry (col->fromBuffer, 0, 0, 0, 0) ;
    }
  else
    { freeint(&col->from) ;
    if (col->from <= 0 ||
	col->from >= arrayMax(bqld->colonnes))
      { messout("Please build from an existing colonne") ;
      col->from = 0 ;
      }
    bqldDefineColonne(TRUE) ;
    if (col->from != old)
      bqld->modified = TRUE ;
    }
}

static BOOL bqldDefineWidth (int n)
{ 
  return (n > 0) ;
}

static BOOL bqldDefineDnaCoord (char *text, int box)
{ 
	return TRUE ;
}

/*
static BOOL bqldDefinePepCoord (char *text, int box)
{ 
	return TRUE ;
}
*/
/************************************************/
/*
static void bqldGetWidths (BQLD bqld)
{ int i, max = arrayMax(bqld->colonnes), level ;
  COL *c ;
  
  for (i = 0; i < max ; i++)
    { c = arrayp(bqld->colonnes, i, COL) ;  
      level= freesettext(c->widthBuffer, "") ;
      if (! freecard(level) ||
	  ! freecheck("i"))
	{ if(messPrompt(messprintf("Width of colonne %d", i), "12", "i"))
	    freeint(&c->width) ;
	  else
	    c->width = 12 ;
	  sprintf(c->widthBuffer,"%d", c->width) ;
	  graphTextEntry (c->widthBuffer, 0, 0, 0, 0) ;
	}
      else
	freeint(&c->width) ;
    }
}
*/
  
/************************************************/

static void bqldDefineCondition (void)
{
  COL *c ;
  char *cp ;
  BQLD bqld = currentBqlDisp("bqldDefineCondition") ;

  c = arrayp(bqld->colonnes, bqld->activeColonne, COL) ;  
  cp = c->conditionBuffer ;
  while (*cp) cp ++ ;
  while (--cp >= c->conditionBuffer && *cp == ' ')
    *cp = 0 ;
  if (*c->conditionBuffer &&
      !condCheckSyntax(messprintf(" %s", freeprotect(c->conditionBuffer))))
    messout("Please correct this syntax error") ;
  bqld->modified = TRUE ;

  return;
} /* bqldDefineCondition */
  
/************************************************/

static void bqldDefineSubtitle (void)
{
  COL *c ;
  char *cp ;
  BQLD bqld = currentBqlDisp("bqldDefineSubtitle") ;

  c = arrayp(bqld->colonnes, bqld->activeColonne, COL) ;  
  cp = c->subtitleBuffer ;
  while (*cp) cp ++ ;
  while (--cp >= c->subtitleBuffer && *cp == ' ')
    *cp = 0 ;
  bqld->modified = TRUE ;

  return;
} /* bqldDefineSubtitle */
  
/************************************************/

static void bqldDefineLegend (void)
{
  COL *c ;
  char *cp ;
  BQLD bqld = currentBqlDisp("bqldDefineSubtitle") ;

  c = arrayp(bqld->colonnes, bqld->activeColonne, COL) ;  
  cp = c->legendBuffer ;
  while (*cp) cp ++ ;
  while (--cp >= c->legendBuffer && *cp == ' ')
    *cp = 0 ;
  bqld->modified = TRUE ;

  return;
} /* bqldDefineLegend */
  
/************************************************/

static void bqldDefineTitle (void)
{ 
  char *cp ;
  BQLD bqld = currentBqlDisp("bqldDefineTitle") ;

  cp = bqld->titleBuffer ;
  while (*cp) cp ++ ;
  while (--cp >= bqld->titleBuffer && *cp == ' ')
    *cp = 0 ;
  bqld->modified = TRUE ;

  return;
} /* bqldDefineTitle */
  
/************************************************/

static void bqldDefineParam (void)
{ 
  char *cp ;
  BQLD bqld = currentBqlDisp("bqldDefineParam") ;

  cp = bqld->paramBuffer ;
  while (*cp) cp ++ ;
  while (--cp >= bqld->paramBuffer && *cp == ' ')
    *cp = 0 ;
  bqld->modified = TRUE ;

  return;
} /* bqldDefineParam */
  
/************************************************/

void bqldShow(void)
     /* private within bqldDisp-package */
{ 
  BQLD bqld = currentBqlDisp("bqldShow") ;

  if (!graphCheckEditors (graphActive(), 0))
    return ;
  /*  bqldGetWidths (bqld) ;*/
  bqld->modified = TRUE ;	/* force recompute */
  if (bqldDoRecompute (bqld))
    bqldDisplayData (bqld) ;

  return;
} /* bqldShow */

/***********/

static void bqldForceOpenColonne (BOOL force)
{
  BQLD bqld = currentBqlDisp("bqldForceOpenColonne2") ;

  bqldInitColonne(bqld) ;
  bqldDefineColonne(force) ;
  bqld->modified = TRUE ;

  return;
} /* bqldForceOpenColonne */

static void bqldOpenColonne (void)
{
  BQLD bqld = currentBqlDisp("bqldOpenColonne") ;

  bqldInitColonne(bqld) ;
  bqldDefineColonne(FALSE) ;
  bqld->modified = TRUE ;

  return;
} /* bqldOpenColonne */

/***********/

static BOOL bqldShiftParameters (char *buffer, int x0, int dx)
{
  int n ;
  char *cp0, *cp, *cq ;
  
  if (! buffer || ! strstr (buffer, "%"))
    return FALSE ;

  cp = cp0 = strnew (buffer, 0) ;
  buffer[0] = 0 ;
  while (*cp)
    {
      cq = strstr (cp, "%") ;
      if (cq) *cq = 0 ;
      strcat (buffer, cp) ;
      if (cq)
	{
	  cp = cq + 1 ;
	  n = 0 ;
	  while (*cp >= '0' && *cp <= '9')
	    n = 10 * n + (*cp++ - '0') ;
	  strcat (buffer, messprintf ("%%%d", n > x0 ? n + dx : n)) ;
	}
      else 
	break ;
    }

  messfree (cp0) ;
  return TRUE ;
} /* bqldShiftParameters */

/***********/

static void bqldInsertColonne (void)
{
  COL *c ;
  int i, nn ;
  BQLD bqld = currentBqlDisp("bqldInsertColonne") ;

  nn = arrayMax(bqld->colonnes) ;
  if (! nn)
    return ;

  c = arrayp(bqld->colonnes, bqld->activeColonne, COL) ;  

  if (!messQuery(messprintf("Do you really want to insert a colonne after colonne %d",
		 bqld->activeColonne + 1)))
    return ;

  c = arrayp(bqld->colonnes, nn, COL) ;   /* create an extra colonne at the end */
  for (i = nn ; i >  bqld->activeColonne + 1 ; i--)
    {
      c = arrayp(bqld->colonnes, i, COL) ;
      *c = *(c-1) ;
    }
  bqldDoInitColonne (bqld, bqld->activeColonne + 1) ;

  for (i = bqld->activeColonne + 2 ; i < arrayMax(bqld->colonnes) ; i++)
    {
      c = arrayp(bqld->colonnes, i, COL) ;  
      if (c->from > bqld->activeColonne + 1)
	c->from++ ;
      bqldShiftParameters (c->conditionBuffer, bqld->activeColonne + 1, 1) ;
      bqldShiftParameters (c->dna1Buffer, bqld->activeColonne + 1, 1) ;
      bqldShiftParameters (c->dna2Buffer, bqld->activeColonne + 1, 1) ;
       if (stackExists(c->tagStack))
	 {
	   bqldShiftParameters (c->tagTextBuffer, bqld->activeColonne + 1, 1) ;
	   stackClear (c->tagStack) ;
	   pushText (c->tagStack, c->tagTextBuffer) ;
	 }
    }
  isCheckEditor = TRUE ;
  
  bqld->modified = TRUE ;
  bqldDefineColonne(FALSE) ;
}

/***********/

static void bqldRemoveColonne (void)
{ COL *c, *c1 ;
  int i ;
  BQLD bqld = currentBqlDisp("bqldRemoveColonne") ;

  if (!arrayMax(bqld->colonnes))
    return ;
  c = arrayp(bqld->colonnes, bqld->activeColonne, COL) ;  
  if (!bqld->activeColonne)
    { messout ("You cannot remove the first colonne") ;
      return ;
    }
  if (!messQuery(messprintf("Do you really want to remove colonne %d",
		 bqld->activeColonne + 1)))
    return ;
  bqldDestroyCol(c) ;
  for ( i = bqld->activeColonne, c = arrayp(bqld->colonnes, i, COL) ;  
        i + 1 < arrayMax(bqld->colonnes) ; i++)
    { c1 = c++ ;
      memcpy(c1, c, sizeof(COL)) ;
    }

  arrayMax(bqld->colonnes) -- ;

  /*    for (i = bqld->activeColonne + 1 ; i < arrayMax(bqld->colonnes) ; i++)
    { c = arrayp(bqld->colonnes, i, COL) ;  
      if (c->from > bqld->activeColonne)
	c->from -- ;
      else if (c->from == bqld->activeColonne)
 messout ("You should modify colonne %d, which is derived from the colonne you just destroyed", i + 1) ;
    }*/
  for (i = bqld->activeColonne ; i < arrayMax(bqld->colonnes) ; i++)
    {
      c = arrayp(bqld->colonnes, i, COL) ;  
      if (c->from > bqld->activeColonne + 1)
	c->from -- ;
      else if (c->from == bqld->activeColonne + 1)
	messout ("You should modify colonne %d, which is derived from the colonne you just destroyed", i + 2) ;
      bqldShiftParameters (c->conditionBuffer, bqld->activeColonne, -1) ;
      bqldShiftParameters (c->dna1Buffer, bqld->activeColonne, -1) ;
      bqldShiftParameters (c->dna2Buffer, bqld->activeColonne, -1) ;
    }
  isCheckEditor = TRUE ;
  /*     bqld->activeColonne-- ; mhmp 21.04.98  12.05.99*/
  bqld->activeColonne = arrayMax(bqld->colonnes) - 1 ;
  bqld->modified = TRUE ;
  bqldDefineColonne(FALSE) ;
}

/**********************************************************/

static void shouldWeDestroy(void)
{ BQLD bqld = currentBqlDisp("shouldwedestroy") ;

  if (bqld->quitWithoutConfirmation ||
      messQuery("Do you really want to quit the Table_Maker ?"))
    graphDestroy() ;
}

/***********************************************************/
/*********** Input Output of the definitions ***************/ 
/* The actual operation are in a separate non graphic file */


static char dirName[DIR_BUFFER_SIZE], fileName[FIL_BUFFER_SIZE] ;
static BOOL firstDirPass = TRUE ; 

static void bqldDoGetDefinitions (KEY key)
{ BQLD bqld = currentBqlDisp("bqldDoGetDefinitions") ;

  displayUnBlock () ;
  if (class(key) != _VTable)
   messout ("Sorry, you must pick a Table object") ;
  else if (!iskey(key))
    messout ("Sorry, this object is empty.") ;
  else
    {
      bqldDoReadDefinitions (bqld, key, 0, 0,"", TRUE) ; /* do not substitute param */
  /* mhmp 12.05.99*/
      bqld->activeColonne = arrayMax(bqld->colonnes) - 1 ;
      bqldDefineColonne(TRUE) ;  /* To redraw */
      bqld->fileName = 0 ;
    }
}

static void bqldGetDefinitions (void)
{ FILE *f ;
  char *cp ;
  BQLD bqld = currentBqlDisp("bqldGetDefinitions") ;
  Graph old = graphActive () ;
  KEYSET ks = query (0,">?Table") ;
  cp = bqld->fileName ;
  if( cp && *cp && strlen(cp) && bqld->modif
      && messQuery (messprintf("%s not saved, Save ?", cp)) )
    { 
      f = filqueryopen(bqld->dirName, bqld->fileName, "def", "w",
			"Choose a File to store the Table Definition") ;
      if (!f)
	messout("Sorry, not done") ;
      else
	bqldDoSaveDefinitions (bqld, f) ;
    }  
  displayCreate(DtKeySet) ;
  graphRetitle("Table definitions") ;
  keySetShow (ks,0) ;
  keySetSelect () ;
  graphActivate (old) ;

  displayBlock (bqldDoGetDefinitions,
		"Pick a Table object") ;
}

static void bqldReadDefinitions (void)
{ FILE *f ;
  char *cp ;/*mhmp 11.98*/

  BQLD bqld = currentBqlDisp("bqldReadDefinitions") ;
  
  if (firstDirPass)
    { if (filName ("wquery", 0, "r"))
	strncpy (dirName, filName ("wquery", 0, "r"), DIR_BUFFER_SIZE-1) ;
      strcpy (fileName, "table") ;
      firstDirPass = FALSE ;
    }
  cp = bqld->fileName ; /* mhmp 26.11.98 pour saver le fichier en cours */
  if( cp && *cp && strlen(cp) && bqld->modif
      && messQuery (messprintf("%s not saved, Save ?", cp)) )
    { 
      f = filqueryopen(bqld->dirName, bqld->fileName, "def", "w",
			"Choose a File to store the Table Definition") ;
      if (!f)
	messout("Sorry, not done") ;
      else
	bqldDoSaveDefinitions (bqld, f) ;
      bqld->modif = FALSE ; 
      return ;
    }  

  f = filqueryopen (dirName, fileName, "def","r", 
		    "Choose a Table-definition file") ;
  if (!f)
    return ;

  bqld->modif = FALSE ;
  bqld->fileName = fileName ;
  bqld->dirName = dirName ;
  if (!bqldDoReadDefinitions (bqld, 0, f, 0,"", TRUE))
 /*  do not substitute param */
    {
      messout ("Sorry, your definition file is not correct") ;
	bqldDestroy (bqld) ;
	bqldDispCreate (TRUE) ;
	return ;
    }
  /* mhmp 12.05.99*/
  bqld->activeColonne = arrayMax(bqld->colonnes) - 1 ;
  bqldDefineColonne(TRUE) ;  /* To redraw */
  if (!graphCheckEditors (graphActive(), 0))
    {
      messout ("Sorry, your definition file is not correct") ;
      bqldDestroy (bqld) ;
      bqldDispCreate (TRUE) ;
    }
  /* repeat after eventual destroy */
  bqld = currentBqlDisp("bqldReadDefinitions") ;
  bqld->modif = FALSE ;
  bqld->fileName = fileName ;
  bqld->dirName = dirName ;
}

static void bqldWriteDefinitions (void)
{ FILE *f = 0 ;
  COL *c = 0 ;
  int  maxCol ;
  
  BQLD bqld = currentBqlDisp("bqldWriteDefinitions") ;

 
  if (!graphCheckEditors (graphActive(), 0))
      return ; 
  maxCol = arrayMax (bqld->colonnes) ;
  if (maxCol)
    c = arrp(bqld->colonnes, 0 , COL) ;
   /* if(!bqld->modif) return ; mieg: I may want to save an unmodified file */
  if (!c || !c->type)
    return ; 
  if (firstDirPass)
    { if (filName ("wquery", 0, "r"))
	strncpy (dirName, filName ("wquery", 0, "r"), DIR_BUFFER_SIZE-1) ;
      strcpy (fileName, "table") ;
      firstDirPass = FALSE ;
    }
  /*f = filqueryopen(dirName, "", "def", "w", mhmp 27.03.98*/
  /*  cp = fileName + strlen(fileName) - 4 ;
  if (cp > fileName && !strcmp(cp, ".def")) *cp = 0 ;mhmp 19.11.98 */  

  f = filqueryopen(dirName, fileName, "def", "w",
		   "Choose a File to store the Table Definition") ;
  if (!f)
    return ;
  bqld->modif = FALSE ;
  bqldDoSaveDefinitions (bqld, f) ;
}

static void bqldSaveDefinitions (void)
{ static KEY tableKey = 0 ;
  char *nam ;
  BQLD bqld = currentBqlDisp("bqldSaveDefinitions") ;
  if (!graphCheckEditors (graphActive(), 0))
      return ;
  nam = tableKey ? name(tableKey) : "" ;
  if (!messPrompt (
    "Please give a Name to save this table as an acedb object",
    nam ? nam : "", "w")) return ;
  nam = strnew (freeword(), 0) ;
  if (lexword2key (nam, &tableKey, _VTable))
    { if (!messQuery ("This table allready exists, do you want to overwrite it"))
      goto abort ;
    }
  sessionGainWriteAccess() ;
  lexaddkey (nam, &tableKey, _VTable) ;
  bqldDoSaveInObj (bqld, tableKey) ;
abort: 
  messfree (nam) ;
}

/************/

static void bqldExportBql (void)
{
  COL *c = 0 ;
  int  maxCol ;
  BQLD bqld = currentBqlDisp("bqldEportBql") ;

 
  if (!graphCheckEditors (graphActive(), 0))
      return ; 
  maxCol = arrayMax (bqld->colonnes) ;
  if (maxCol)
    c = arrp(bqld->colonnes, 0 , COL) ;
   /* if(!bqld->modif) return ; mieg: I may want to save an unmodified file */
  if (!c || !c->type)
    return ; 
 
  bqldDoExportBql (bqld, 0, TRUE) ;
} /* bqldExportBql */

/************/

static void bqldGraphTop (void)
{
  graphGoto (1,1) ;
}
/************/

static void bqldNewTable (void)
{
  bqldDispCreate (0) ; 
}

/************/

extern void bqldImportKeySet (void) ;
/* static void bqldComment(void) ; */
static void bqldCallImportKeySet (void) ;

static MENUOPT bqldDefMenu[] =
  {
   {shouldWeDestroy, "Quit"},
   {help, "Help"},
   {graphPrint, "Print"},
   {bqldShow, "Search Whole Class"},
   {bqldCallImportKeySet, "Search Active KeySet"},
   {bqldGetDefinitions, "Standard Query"},
   {bqldSaveDefinitions, "Save Standard Query"},/* mhmp 09.07.02  write save*/
   {bqldReadDefinitions, "Read Query from file"},
   {bqldWriteDefinitions, "Write Query to file"},/* mhmp 09.07.02  write save*/
   {bqldExportBql, "Write Bql"},
   {bqldOpenColonne, "Add Column"}, 
   {bqldInsertColonne, "Insert Column"}, 
   {bqldRemoveColonne, "Suppress Column"}, 
   {bqldNewTable, "New table"},
/*    {bqldComment, "Comments"}, Not very useful */
   {0, 0}
   } ;

/**********************************************************/
static void bqldCallImportKeySet (void)
{
  if (!graphCheckEditors (graphActive(), 0))
      return ;
  bqldImportKeySet () ;
} 
/*
static void bqldComment(void)
{  int x = 1 , y = 7, maxx = 30  ;
   Stack s ; char *cp ;
  BQLD bqld = currentBqlDisp("bqldComment") ;

  s = bqld->comments ;

  graphClear() ;
  graphHelp("Table_Maker") ;
  graphMenu(bqldDefMenu) ;
  graphButtons (bqldDefMenu, 1, 1, 68) ;

   graphText 
     ("Comments are those lines in the Command file starting, with #.", 3, 4) ;

   if (stackExists(s))
     {
       stackCursor(s, 0) ;
       while (cp = stackNextText(s))
	 {
	   graphText (cp, x, y) ;
	   x += strlen(cp) + 1 ;
	   if (x > maxx)
	     maxx = x ;
	   x = 1 ;
	   y++;
	 }
     }
   else
     graphText("No comments in this file", 3, 6) ;
  graphTextBounds( maxx + 2 , y+3) ;
  graphRedraw() ;
}
*/

/****************************************/
/*
static void setSortColumn (char *buf)
{
  int i ;
  BQLD bqld = currentBqlDisp("setSortColumn") ;

  freeforcecard (buf) ;
  if (!freecheck ("iz") ||
      !freeint (&i) || 
      i < 0 || 
      i > arrayMax(bqld->colonnes))
    messout ("Sorry, value must be between a valid "
	     "column number or 0, which means left "
	     "to right ordering.") ;
  bqld->sortColonne = i ;
}
*/
static BOOL setSortColumn (int n)
{
  BQLD bqld = currentBqlDisp("setSortColumn") ;
  if (n < 0 || n > arrayMax(bqld->colonnes))
    {
      messout ("Sorry, value must be between a valid "
	       "column number or 0, which means left "
	       "to right ordering.") ;
      return FALSE ;
    }
  return TRUE ;
}

/*********************************************/

void bqldDefineColonne (BOOL force)
{ 
  COL *c ;
  int i, max, box, from ;
  float line, centralLine = 1 ;
  BQLD bqld = currentBqlDisp("bqldDefineColonne") ;

  if (!force && !isCheckEditor && !graphCheckEditors (graphActive(), 0))
      return ;
  isCheckEditor = FALSE ;
  if (!arrayMax(bqld->colonnes))
    { bqldForceOpenColonne(force) ;
      return ;
    }
  classeListeInit () ;
  graphClear() ;
  bqld->activeDefBox = 0 ;
  bqld->defBoxPerCol = 0 ;
  graphRegister (PICK, bqldPickDefinition) ;
  graphHelp("Table_Maker") ;
  graphMenu(bqldDefMenu) ;

  box = graphBoxStart() ;
  graphButtons (bqldDefMenu, 1, 1, 68) ;
  graphBoxEnd() ;
  graphBoxDim (box, 0, 0, 0, &line) ;
  line += 1 ;

  graphText ("Sort column:", 3, line) ;
  strncpy (bqld->sortBuffer, 
	   messprintf("%d", bqld->sortColonne), 6) ;
  bqld->sortColonneBox =
    /*    graphTextEntry (bqld->sortBuffer, 6, 16, line, 
		    setSortColumn) ;*/
	graphIntEditor ("", &bqld->sortColonne, 16, line, setSortColumn) ;
  graphText ("F4 to interrupt", 40, line) ;
  line += 1.5 ;

  graphText ("Title", 3, line) ;
  bqld->titleBox = 
    graphTextScrollEntry (bqld->titleBuffer, 300, 60, 9, line, 
		    TEXT_ENTRY_FUNC_CAST bqldDefineTitle) ;
  line += 1.5 ;

  graphText ("Parameters", 3, line) ;
  bqld->paramBox = 
    graphTextEntry (bqld->paramBuffer, 179, 14, line, 
		    TEXT_ENTRY_FUNC_CAST bqldDefineParam) ;
  line += 1.5 ;

  bqld->definitionBox = graphBoxStart() ;
  graphText ("Column", 1, line) ;
  line += 1.5 ;

  max = arrayMax(bqld->colonnes) ;
  for (i = 0, c = arrayp(bqld->colonnes, 0, COL) ;
       i < max ; i++, c++, line += 2.5 )
    {
      if (i == bqld->activeColonne)
	centralLine = line ;
      graphText(messprintf("%3d", i + 1), 2, line) ;
 
      graphText ("Title", 7, line) ;
      bqld->subTitleBox = 
	graphTextEntry (c->subtitleBuffer, 59, 17, line, 
			TEXT_ENTRY_FUNC_CAST bqldDefineSubtitle) ;

      graphText ("Legend", 7, line+=1.3) ;
      bqld->legendBox = 
	graphTextScrollEntry (c->legendBuffer, 1023, 59, 17, line, 
			TEXT_ENTRY_FUNC_CAST bqldDefineLegend) ;
  
      graphText ("Width", 7, line += 1.3) ;
      sprintf(c->widthBuffer,"%d", c->width) ;
      bqld->widthBox = 
	/*	graphTextEntry (c->widthBuffer, 4, 14, line, 
			TEXT_ENTRY_FUNC_CAST bqldDefineWidth) ;mhmp 30.03*/
	graphIntEditor ("", &c->width, 14, line, bqldDefineWidth) ;
  
      box = bqld->hideBox = graphBoxStart() ;
      graphTextPtrPtr(&c->hiddenp, 24, line, 7) ;
      c->hiddenp = freekey2text(c->hidden, hideChoice) ;
      graphBoxEnd() ;
      graphBoxBox(box) ;
      graphBoxFreeMenu(box, (FreeMenuFunction) bqldDefineHide,
		       hideChoice) ;
      
      box = bqld->optionalBox = graphBoxStart() ;
      graphTextPtrPtr(&c->optionalp, 33, line, 9) ;
      c->optionalp = freekey2text(c->mandatory, optionalChoice) ;
      graphBoxEnd() ;
      if (i) 
	{ graphBoxBox(box) ;
	  graphBoxFreeMenu(box, (FreeMenuFunction) bqldDefineOptional,
			   optionalChoice) ;
	}
      else
	graphBoxMarkAsClear(box) ;
      box = bqld->typeBox = graphBoxStart() ;
      switch (c->showType)
	{
	case SHOW_ALL:
	  if (c->type != 'c')
	    {
	      graphTextPtrPtr(&c->typep, 45, line, 7) ;
	      graphTextPtrPtr(&c->classp, 53, line, 26) ;   
	    }
	  else 
	    {  
	      c->showtypep = "COUNT" ;
	      c->realType = c->type ; /*mhmp 22.02.99 */
	      graphTextPtrPtr(&c->showtypep, 45, line, 16) ;
	    }
	  break ;
	case SHOW_DNA:
	  c->showtypep = "DNA" ;
	  graphTextPtrPtr(&c->showtypep, 45, line, 16) ;
	  break ;
	case SHOW_PEP:
	  c->showtypep = "Peptide" ;
	  graphTextPtrPtr(&c->showtypep, 45, line, 16) ;
	  break ;
	case SHOW_COPY:
	  c->showtypep = "COPY" ;
	  graphTextPtrPtr(&c->showtypep, 45, line, 16) ;
	  break ;
	case SHOW_MIN:
	case SHOW_MAX:
	case SHOW_AVG:
	case SHOW_VAR:
	case SHOW_SUM:
	case SHOW_COMPUTE:
	  c->showtypep = freekey2text (c->showType, typeChoice2) ;
	  c->nonLocal = TRUE ;  
	  graphTextPtrPtr(&c->showtypep, 45, line, 16) ;
	  break ;
	case SHOW_MULTI:
	  c->showtypep = freekey2text (c->showType, typeChoice2) ;
	  c->type = c->realType ;
	  c->nonLocal = TRUE ;  
	  graphTextPtrPtr(&c->showtypep, 45, line, 16) ;
	}
      c->typep = freekey2text(c->type, typeChoice) ;
      c->classp = c->classe && ace_lower(c->type) == 'k' ?
	name(c->classe) : 0 ;
      graphBoxEnd() ;
      graphBoxBox(box) ;
      if (!i)
	graphBoxFreeMenu(box, (FreeMenuFunction) bqldDefineClass,
		       classeListe) ;
      else
	graphBoxFreeMenu(box, (FreeMenuFunction) bqldDefineType,
		       typeChoice2) ;

      box = bqld->extendBox = graphBoxStart() ;
      graphTextPtrPtr(&c->extendp, 9, line += 1.5, 8) ;
      c->extendp = freekey2text(c->extend, extendChoice) ;
      graphBoxEnd() ;
      if (i) graphBoxBox(box) ;
      if (i)
	graphBoxFreeMenu(box, (FreeMenuFunction) bqldDefineExtend,
			 extendChoice) ;
      else
	graphBoxMarkAsClear(box) ;
      if (!c->from) c->from = 1 ;
      sprintf(c->fromBuffer,"%d", from = c->from) ;
      bqld->fromBox = box =
	graphTextEntry (c->fromBuffer, 4, 18, line, 
			TEXT_ENTRY_FUNC_CAST bqldDefineFrom) ;
      /* intEditor does not allow to call a full redraw */

      if (!i)
	graphBoxMarkAsClear(box) ;      
      box = bqld->dnaBox = graphBoxStart() ;
      graphText ("DNA" , 28, line) ;
      graphBoxEnd() ; 
       if (i &&
	  c->extend == 'f' &&
	  c->from > 0 &&
	  arrp(bqld->colonnes, c->from - 1 , COL) &&
	  dnaInClass (superClass( arr(bqld->colonnes, c->from - 1 , COL).classe))
	  )  { graphBoxBox(box) ; }
      else graphBoxMarkAsClear(box) ;
      box = bqld->pepBox = graphBoxStart() ;
      graphText ("PEP" , 28, line) ;
      graphBoxEnd() ; 
       if (i &&
	  c->extend == 'f' &&
	  c->from > 0 &&
	  arrp(bqld->colonnes, c->from - 1 , COL) &&
	  pepInClass (superClass( arr(bqld->colonnes, c->from - 1 , COL).classe))
	  )  { graphBoxBox(box) ; }
      else graphBoxMarkAsClear(box) ;
      box = bqld->tagBox = graphBoxStart() ;
      if (i) graphText ("Tag:" , 33, line) ;
      graphBoxEnd() ;
      if (i)  { graphBoxBox(box) ; } /* {} needed around stupid macro bobox */
      else graphBoxMarkAsClear(box) ;

      /* a left over which seems useless as of  mars 22 2002

      {
      COL * fromC ;
      if (from > 0)
	fromC = arrp(bqld->colonnes, from - 1, COL) ;
      else 
	fromC = 0 ;

	if (!c->extend)
	  if (!i || !fromC || ace_lower(fromC->type) != 'k' || !fromC->classe)
	graphBoxMarkAsClear(box) ;
	}
      */
      if (!stackExists(c->tagStack))
	 strcpy(c->tagTextBuffer, "?") ;
      else
	strncpy(c->tagTextBuffer, stackText(c->tagStack,0), 359) ;
      
      /* always create all boxes for correct count and clear the unneeded */
      if (c->type == 'D' || c->type == 'P')
	bqld->dna1Box = box =
	  graphTextScrollEditor ("from", c->dna1Buffer, 200, 10, 38, line, bqldDefineDnaCoord) ;
      else
	{
	  graphBoxStart () ; graphBoxEnd() ; 
	  graphBoxStart () ; graphBoxEnd() ; 
	}
      if (c->type == 'D' || c->type == 'P')
	bqld->dna2Box = box =
		graphTextScrollEditor ("to", c->dna2Buffer, 200, 10, 54, line, bqldDefineDnaCoord) ; /* mhmp 11.07.02 from --> to */
      else
	{
	  box = graphBoxStart () ; graphBoxEnd() ; 
	  box = graphBoxStart () ; graphBoxEnd() ; 
	}
      box = graphBoxStart () ; graphBoxEnd() ; 
      if (i && c->type != 'D'&& c->type != 'P')
	{
	  bqld->tagTextBox = box =
	    graphTextScrollEntry(c->tagTextBuffer, 359, 
				 300, 38, line, 
				 TEXT_ENTRY_FUNC_CAST bqldEditTagText) ;
	}
      else
	{
	  box = graphBoxStart () ; graphBoxEnd() ; 
	  box = graphBoxStart () ; graphBoxEnd() ; 
	}
      box = graphBoxStart () ; graphBoxEnd() ; 
    
      if (!i)
	graphText ("To restrict the search, use the condition box or search a keyset", 9 , line) ;
      graphText ("Condition", 7, line += 1.3) ;
      bqld->conditionBox = 
	graphTextScrollEntry (c->conditionBuffer, 359, 
			      300, 17, line, 
			TEXT_ENTRY_FUNC_CAST bqldDefineCondition) ;
      if (c->type == 'D' || c->type == 'P')
	graphText ("limit to sequence matching a UNIX \"RegExp\", i.e: atgc  ^g[tc].*ag$  ", 7, line += 1.1) ;

          /* last paragraph of the loop */
      box = graphBoxStart() ;
      graphBoxEnd() ;
      if (!bqld->defBoxPerCol)
	bqld->defBoxPerCol = box - bqld->definitionBox ;
    }
  graphBoxEnd() ;
  graphButton("Add column", bqldOpenColonne, 2, line) ;
  graphButton("Page top", bqldGraphTop, 15, line) ;
  graphButton("Search Whole Class", bqldShow, 26, line) ;
  graphButton("Search Active KeySet", bqldCallImportKeySet, 48, line) ;
  graphTextBounds (300, line + 3) ;

  graphRedraw() ;
  graphGoto (1, centralLine) ;
  return ;
} 
static void bqldExportKeySet (void)
{
}

#endif /* JUNK */

/**********************************************************/
/**************** Select **********************************/

static void bqldSelect (BQLD bqld)
{
} /* bqldSelect */

/**********************************************************/
/****************** From **********************************/

static void bqldNewVar (char *cp)
{
  BQLD bqld = currentBqlDisp("bqldNewvar") ;
  Array vars = bqld->vars ;
  BqlVar *vp = 0 ;
  int ii ;

  for (ii = 0 ; ii < arrayMax (vars) ; ii++)
    {
      vp = arrp (vars, ii, BqlVar) ;
      if (cp == vp->buf)
	break ;
    }
  vp->type = *cp ? 2 : 1 ; /* this var is defined or waiting */
  if (ii == arrayMax (vars) - 1)
    {
      ii++ ;
      vp = arrayp (vars, ii, BqlVar) ;
      vp->type = 1 ;
    }
  bqldDraw (bqld) ; 

  return  ;  
} /* bqldNewVar */
  
/**********************************************************/

static void bqldFromNewClass (void)
{
  BQLD bqld = currentBqlDisp("bqldSyntaxHelpOpen") ;
    
  /* bqld->type = 31 ; */
  bqldDraw (bqld) ; 
} /* bqldFromNewClass */

/**********************************************************/

static void bqldFromNewVar (void)
{
  BQLD bqld = currentBqlDisp("bqldFromNewVar") ;
    
  /* bqld->type = 41 ; */
  bqldDraw (bqld) ; 
} /* bqldFromNewVar */

/**********************************************************/

static void bqldFromNewConstant (void)
{
  BQLD bqld = currentBqlDisp("bqldFromNewConstant") ;
    
  /* bqld->type = 51 ; */
  bqldDraw (bqld) ; 
} /* bqldFromNewConstant */

/**********************************************************/

static void bqldFromNewEquation (void)
{
  BQLD bqld = currentBqlDisp("bqldFromNewEquation") ;
    
  /* bqld->type = 61 ; */
  bqldDraw (bqld) ; 
} /* bqldFromNewEquation */

/**********************************************************/

static void bqldFromNewClear (void)
{
  BQLD bqld = currentBqlDisp("bqldFromNewClear") ;
    
  /* bqld->type = 2 ; */
  bqldDraw (bqld) ; 
} /* bqldFromNewClear */

/**********************************************************/

static void bqldFromNewHelp (void)
{
  BQLD bqld = currentBqlDisp("bqldFromNewHelp") ;
    
  
  bqldDraw (bqld) ; 
} /* bqldFromNewHelp */

/**********************************************************/

static void bqldFrom (BQLD bqld)
{
  int ii, line = bqld->line ;
  BqlVar *vp ;
  Array vars = bqld->vars ;
  
  for (ii = 0, vp = arrp (vars, 0, BqlVar) ; ii < arrayMax (vars) ; ii++, vp++)
    {
      switch (vp->type)
	{
	case 1:
	  graphText ("Name a new variable", 2, line ) ;
	  graphTextEntry(vp->buf,8,24,line,bqldNewVar) ;
	  line++ ;
	  break ;
	case 2:
	  graphText ("Define variable", 4, line) ;
	  graphText (vp->buf, 25, line) ;
	  graphButton ("From class", bqldFromNewClass, 30, line) ;
	  graphButton ("From previous variable", bqldFromNewVar, 30, line+1) ;
	  graphButton ("As a constant", bqldFromNewConstant, 30, line+2) ;	  
	  graphButton ("From equation", bqldFromNewEquation, 30, line+3) ;
	  graphButton ("Clear", bqldFromNewClear, 60, line) ;
	  graphButton ("Help", bqldFromNewHelp, 60, line+1) ;
	  line += 6 ;
	  break ;
	case 31:
	  graphText (messprintf ("For all %s in class ", vp->buf), 4, line) ;
	  vp->box = graphButton ("Select a class...", 0, 40, line - .1);
	  break ;
	case 41:
	  graphText (messprintf ("For all %s in %s->class ", vp->buf), 4, line) ;
	  vp->box = graphButton ("Select a class...", 0, 40, line - .1);
	  break ;
	case 51:
	  graphText (messprintf ("Constant %s in class ", vp->buf), 4, line) ;
	  vp->box = graphButton ("Select a class...", 0, 40, line - .1);
	  break ;
	case 61:
	  graphText (messprintf ("Equation %s in class ", vp->buf), 4, line) ;
	  vp->box = graphButton ("Select a class...", 0, 40, line - .1);
	  break ;
	}
    }


  bqld->line = line ;
} /* bqldFrom */

/**********************************************************/
/**************** Where **********************************/

static void bqldWhere (BQLD bqld)
{
} /* bqldWhere */

/**********************************************************/
/***************** Syntax Help ****************************/

static void bqldSyntaxHelpDisplay (BQLD bqld)
{
  int line = bqld->line ;
  graphText ("hello from Select ... from ... where ..", 3, line++) ;
  bqld->line = line ;
} /* bqldSyntaxHelpDisplay */

/**********************************************************/

static void bqldSyntaxHelpOpen (void)
{
  BQLD bqld = currentBqlDisp("bqldSyntaxHelpOpen") ;
    
  bqld->syntaxHelp = TRUE ;
  bqldDraw (bqld) ; 
} /* bqldSyntaxHelpOpen */

/**********************************************************/

static void bqldSyntaxHelpClose (void)
{
  BQLD bqld = currentBqlDisp("bqldSyntaxHelpClose") ;

  bqld->syntaxHelp = FALSE ;
  bqldDraw (bqld) ; 
} /* bqldSyntaxHelpClose */

/**********************************************************/

static void bqldSyntaxHelp (BQLD bqld)
{
  if (bqld->syntaxHelp)
    {
      graphButton("Close Help", bqldSyntaxHelpClose, 48, 1) ;
      bqld->line += 2 ;
      bqldSyntaxHelpDisplay (bqld) ;
    }
  else
    {
      graphButton("Help", bqldSyntaxHelpOpen, 48, 1) ;
      bqld->line += 2 ;
    }
} /* bqldSyntaxHelp */

/***********************************************************/

  static void bqldSyntax (BQLD bqld)
{
  bqld->line = 1 ;

  graphText ("Select ... from ... where ..", 3, 1) ;
  graphText ("on the Acedb Query Language syntax", 54 + 6 * bqld->syntaxHelp, 1) ;
  
  bqldSyntaxHelp (bqld) ;

  
} /* bqldSyntax */

/***********************************************************/
/***********************************************************/

static void bqldDraw (BQLD bqld)
{
  graphClear () ;

  bqldSyntax (bqld) ;
  bqldSelect (bqld) ;
  bqldFrom (bqld) ;
  bqldWhere (bqld) ;

  graphRedraw () ;
} /* bqldDraw */

/***********************************************************/
/************** Display Graph Creationand Destruction ******/
/***********************************************************/

static void bqldPick (int k)
{
  return ;
}
static BQLD currentBqlDisp (char *caller)
{
  /* find and verify BQLD struct on active graph */
  BQLD bqld = 0 ;

  if (!(graphAssFind(&GRAPH2BQLD_ASSOC, &bqld)))
    messcrash("%s() could not find BQLD on graph", caller);
  if (!bqld)
    messcrash("%s() received NULL BQLD pointer", caller);
  if (bqld->magic != BQLD_MAGIC)
    messcrash("%s() received non-magic BQLD pointer", caller);
  
  return bqld; 
} /* currentBqlDisp */

/*****************************************/

static void bqldDestroy (void)
{    
  Graph bqlG ;
  BQLD bqld = currentBqlDisp("bqldDestroy");
  
  bqlG = bqld->graph ;
  if (graphExists(bqlG))
    {
      graphActivate(bqlG) ;
      graphDestroy() ;
    }
  ac_free (bqld->h) ;
  bqld->magic = 0 ;
  ac_free (bqld) ;

  return;
} /* bqldDestroy */

/*****************************************/
/******** public interface ***************/
void bqldCreate (void)
{
  BQLD bqld = 0 ;

  if (!displayCreate (DtBqlDisp)) 
    return ;

  bqld = (BQLD) halloc (sizeof (struct bqldStruct), 0) ;
  memset (bqld, 0, sizeof (struct bqldStruct)) ;

  bqld->magic = BQLD_MAGIC ;
  bqld->h = ac_new_handle () ;
  bqld->vars = arrayHandleCreate (32, BqlVar, bqld->h) ;
  array (bqld->vars, 0, BqlVar).type = 1 ;  /* make room */
  graphRegister (DESTROY, bqldDestroy) ;
  bqld->graph = graphActive() ;
  /*  if (oldGraph)
    graphAssRemove (&GRAPH2BQLD_ASSOC) ;
  */
  graphAssociate (&GRAPH2BQLD_ASSOC, bqld) ;
  graphRegister (PICK, bqldPick) ;
  graphRegister (RESIZE, bqldDraw) ;
  graphRegister (DESTROY, bqldDestroy) ;
  
  bqldDraw (bqld) ;
  
  return;
} /* bqldDispCreate */

/*************** entry point ****************/
BOOL bqlDisplay (KEY new, KEY old, BOOL isOldGraph) 
{
  bqldCreate () ;
  return 0 ;
} /* bqlDisplay */

/**********************************************************/
/**********************************************************/
/**********************************************************/

 
