/*  File: sprdctrl.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: part of the spreadDisp-package
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 22 17:40 1998 (fw)
 *	-	removed (GraphFunc) typecasts on graphTextEntry() etc. callbacks
 *		replacing them with TEXT_ENTRY_FUNC_CAST on spreadDefineTitle() etc.
 * * Feb 28 01:22 1996 (mieg)
 * Created: Thu Jun 11 22:04:35 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* @(#)sprdctrl.c	1.42 10/29/97 */

#include "acedb.h"

#include "display.h"
#include "spreaddisp_.h"

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

static void spreadChooseTag(SPREAD spread, int box, COL *c, COL *fromC, int cnt) ;
static void spreadChooseDna(SPREAD spread, int box, COL *c, COL *fromC, int cnt) ;
static void spreadChoosePep(SPREAD spread, int box, COL *c, COL *fromC, int cnt) ;

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

static void spreadDoInitColonne (SPREAD spread, int newColonne)
     /* private within spreadDisp-package */
{
  COL *c ;
  Array t ;

  c = arrayp(spread->colonnes, newColonne, COL) ;  
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

  if (arrayExists(t = spread->tableau))
    {
      int i = arrayMax(t) ;
      int max = arrayMax (spread->colonnes) ;

      while (i--)
	array(array(t,i,Array), max, SPCELL).u.k = 0 ;
    }
  c->width = 12 ;
  *c->conditionBuffer = 0 ;
  spread->modified = TRUE ;

  return;
} /* spreadDoInitColonne */


void spreadInitColonne (SPREAD spread)
     /* private within spreadDisp-package */
{
  COL *c ;
  int max ;

  max = arrayMax(spread->colonnes) ;
  array(spread->pos2col,max, int) = max ;
  spread->activeColonne = max ;

  c = arrayp(spread->colonnes, spread->activeColonne, COL) ;  
  if (c) {} ; /* for compiler happiness */
  spreadDoInitColonne (spread, max) ;

  return;
} /* spreadInitColonne */

/*****************************************************/

static COL *spreadDefineActivate(int box)
{
  COL *c ; 
  SPREAD spread = currentSpread("spreadDefineActivate") ; 

  spread->activeColonne = spread->defBoxPerCol ? 
      (box - spread->definitionBox -1) / spread->defBoxPerCol : 0 ;
  c = arrayp(spread->colonnes, spread->activeColonne, COL) ;  
  if (spread->activeDefBox) 
    graphBoxDraw(spread->activeDefBox, BLACK, WHITE) ; 
  spread->activeDefBox = box ;

  return c ;
} /* spreadDefineActivate */

/******** Optional ********/

static FREEOPT optionalChoice[] =
{ 
  {3, "Define"},
  {1, "Optional"},
  {2, "Mandatory"},
   {0, "Null"}
} ;

static void spreadDefineOptional (KEY k, int box)
{
  COL *c = spreadDefineActivate(box) ;
  KEY old = c->mandatory ;
  SPREAD spread = currentSpread("defineoptional") ;

  c->mandatory = k ;
  c->optionalp = freekey2text(k, optionalChoice) ;
  graphBoxDraw(box, BLACK, LIGHTBLUE) ;
  if (k != old)
    {
      spread->modified = TRUE ;
      spread->modif = TRUE ; /*30.11.98 */
    }
}

/******** Extend ********/

static void spreadDefineExtend (KEY k, int box)
{
  COL *c, *fromC ;
  KEY old ;
  SPREAD spread = currentSpread("defineextend") ;

  /* mhmp 11.01.98 suite a un crash (pas reproduit) */
  if (box <= 0) 
    {
      messout ("Warning: spreadDefineExtend called with box <= 0") ;
      return ;
    }
  c = spreadDefineActivate(box) ; 
  old = c->extend ;
  c->extend = k ;
  c->extendp = freekey2text(k, extendChoice) ;
  graphBoxDraw(box, BLACK, LIGHTBLUE) ;
  if (k != old)
    {
      spread->modified = TRUE ;
      spread->modif = TRUE ; /*30.11.98 */ 
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
	  fromC =  arrp(spread->colonnes, c->from - 1 , COL) ;
	  c->type = fromC->type ;
	  c->classe = fromC->classe ;
	  c->realType = c->type ; 
	  break ;
	case 'p': /* compute */
	  c->from = 1 ; 
	  c->showType = SHOW_COMPUTE ;
	  c->showtypep = freekey2text(k, typeChoice2) ; 
	  c->nonLocal = TRUE ;  
	  fromC =  arrp(spread->colonnes, c->from - 1 , COL) ;
	  c->type = 'f' ;
	  c->classe = 0 ;
	  c->realType = 'f' ;
	  break ;
	}      
      stackDestroy (c->tagStack) ;
      c->tagStack = stackCreate (50) ;
    }

  spreadDefineColonne(TRUE) ;
} /* spreadDefineExtend */

/******** Hide ********/

static FREEOPT hideChoice[] = { 
  {2, "Hide"},
  {TRUE, "Hidden"},
  {FALSE, "Visible"}
} ;

static void spreadDefineHide (KEY k, int box)
{
  COL *c = spreadDefineActivate(box) ;
  KEY old = c->hidden ;
  SPREAD spread = currentSpread("definehide") ;

  c->hidden = k ;
  c->hiddenp = freekey2text(k, hideChoice) ;
  graphBoxDraw(box, BLACK, LIGHTBLUE) ;
  if (k != old)
    {
      spread->modified = TRUE ; 
      spread->modif = TRUE ; /*30.11.98 */
    }
} /* spreadDefineHide */

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

static void spreadDefineType (KEY k, int box)
{
  KEY old ;
  COL *c = spreadDefineActivate(box) ;
  SPREAD spread = currentSpread("defineType") ;

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
	  spread->modified = TRUE ;
	  spread->modif = TRUE ; /*mhmp 30.11.98 */
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
      spread->modified = TRUE ;
      spread->modif = TRUE ; /*mhmp 30.11.98 */
      c->nonLocal = TRUE ;  
      c->showtypep = freekey2text(k, typeChoice2) ;
      {
	int i1 ;
	COL * c1;
	for (i1 = spread->activeColonne + 1 ; i1 < arrayMax(spread->colonnes) ; i1++)
	  { 
	    c1 = arrayp(spread->colonnes, i1, COL) ;  
	    if (c1->from == spread->activeColonne + 1)
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
      spread->modified = TRUE ;
      spread->modif = TRUE ; /* mhmp 30.11.98 */
    }

  graphBoxDraw(box, BLACK, LIGHTBLUE) ;
  if (1) spreadDefineColonne(TRUE) ;
}

/********* Pick boxes ************/

static void spreadPickDefinition(int box, double x, double y)
{
  COL *c ; int i ;
  SPREAD spread = currentSpread("spreadPickDefinition") ;

  if (!graphCheckEditors (graphActive(), 0))
    return ;

  if (!box)
    return ;
  spread->modif = TRUE ; /*mhmp 25.11.98 */
  if (box == spread->titleBox)
    { graphTextEntry(spread->titleBuffer,0,0,0,0) ;
      return ;
    }
  if (box == spread->paramBox)
    { graphTextEntry(spread->paramBuffer,0,0,0,0) ;
      return ;
    }
  if (box == spread->sortColonneBox) return ;
   
  i = spread->activeColonne = (box - spread->definitionBox -1) / spread->defBoxPerCol ; 
  if (i < 0 || i >= arrayMax(spread->colonnes))
    return ;
  c = arrayp(spread->colonnes, spread->activeColonne, COL) ;  

  if (!( (box - spread->widthBox) % spread->defBoxPerCol))
    /*  graphTextEntry(c->widthBuffer,0,0,0,0) ;*/return ;
  else if (!( (box - spread->fromBox) % spread->defBoxPerCol) &&
	   spread->activeColonne)
    /*   graphTextEntry(c->fromBuffer,0,0,0,0) ;*/ return  ;
  else if (!( (box - spread->conditionBox) % spread->defBoxPerCol))
    graphTextScrollEntry(c->conditionBuffer,0,0,0,0,0) ;
  else if (!( (box - spread->subTitleBox) % spread->defBoxPerCol))
    graphTextEntry(c->subtitleBuffer,0,0,0,0) ;
  else if (!( (box - spread->legendBox) % spread->defBoxPerCol))
    graphTextEntry(c->legendBuffer,0,0,0,0) ;
  /* mhmp 12.07.02 peins-le en vert */
  else if (!( (box - spread->dna1Box) % spread->defBoxPerCol))
    graphTextScrollEntry(c->dna1Buffer,0,0,0,0,0) ;
  else if (!( (box - spread->dna2Box) % spread->defBoxPerCol))
    graphTextScrollEntry(c->dna2Buffer,0,0,0,0,0) ;
  else if (!( (box - spread->typeBox) % spread->defBoxPerCol))
    { KEY key = c->type ;
      if (spread->activeColonne)
 	{ if (graphSelect (&key, typeChoice2))
	    spreadDefineType (key, box) ;
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
      spreadDefineColonne(TRUE) ;/* added by mhmp 01.07.98, seems useless, mieg 
is very useful mhmp */	
     }
  else if (!( (box - spread->dnaBox) % spread->defBoxPerCol))
    spreadChooseDna(spread, box, c, arrp(spread->colonnes, c->from - 1 , COL), c->extend) ;
  else if (!( (box - spread->pepBox) % spread->defBoxPerCol))
    spreadChoosePep(spread, box, c, arrp(spread->colonnes, c->from - 1 , COL), c->extend) ;
  else if (!( (box - spread->tagBox) % spread->defBoxPerCol))
    spreadChooseTag(spread, box, c, arrp(spread->colonnes, c->from - 1 , COL), c->extend) ;
  else if (!( (box - spread->tagTextBox) % spread->defBoxPerCol))
    graphTextScrollEntry(c->tagTextBuffer,0,0,0,0,0) ;
  else if (!( (box - spread->hideBox) % spread->defBoxPerCol))
    { 
      c->hidden++ ;
      c->hidden = c->hidden % hideChoice->key ;
      c->hiddenp = freekey2text(c->hidden, hideChoice) ;
      graphBoxDraw(box, BLACK, LIGHTBLUE) ;
    }
  else if (!( (box - spread->optionalBox) % spread->defBoxPerCol))
    { 
      c->mandatory++ ;
      c->mandatory = c->mandatory % optionalChoice->key ;
      c->optionalp = freekey2text(c->mandatory, optionalChoice) ;
      spread->modified = TRUE ;
      graphBoxDraw(box, BLACK, LIGHTBLUE) ;
    }
  else if (!( (box - spread->extendBox) % spread->defBoxPerCol))
    { 
      int extend = 'f' ;

      switch (c->extend)
	{
	case 'f': /* from */     extend = 'r' ; break ;
	case 'r': /* right_of */ extend = 'c' ; break ;
	case 'c': /* copy */     extend = 'f' ; break ;
	}
      c->extend = 0 ; /* insures a redraw */
      spreadDefineExtend (extend, spread->extendBox) ; 
    }
  else
    { if (spread->activeDefBox) 
	graphBoxDraw(spread->activeDefBox, BLACK, WHITE) ; 
      spread->activeDefBox = box ;
    }

  return;
} /* spreadPickDefinition */

/********* DNA/Pep Chooser 1 *******************/

static void spreadChooseDna(SPREAD spread, int box,
			    COL *c, COL *fromC, int continuation)
{
  COL *c1 = c ;

  while (c1->extend != 'f' && fromC->from > 0  && fromC->from < arrayMax(spread->colonnes))
    { c1 = fromC ;
      fromC = arrayp(spread->colonnes, fromC->from - 1, COL) ;   
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
  spread->modified = TRUE ;
  spreadDefineColonne(TRUE) ;

  return;
} /* spreadChooseDna */

/********* DNA/Pep Chooser 1 *******************/

static void spreadChoosePep(SPREAD spread, int box,
			    COL *c, COL *fromC, int continuation)
{
  COL *c1 = c ;

  while (c1->extend != 'f' && fromC->from > 0  && fromC->from < arrayMax(spread->colonnes))
    { c1 = fromC ;
      fromC = arrayp(spread->colonnes, fromC->from - 1, COL) ;   
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
  spread->modified = TRUE ;
  spreadDefineColonne(TRUE) ;

  return;
} /* spreadChoosePep */

/*********Tag Chooser 1 *******************/

static void spreadChooseTag(SPREAD spread, int box,
			    COL *c, COL *fromC, int continuation)
{
  int classe ;
  COL *c1 = c ;
  Stack s1 = 0 ; KEY tag1 = 0, tag2 ; int type1 ;

  while (c1->extend && c1->extend != 'f' && fromC->from > 0  && fromC->from < arrayMax(spread->colonnes))
    { c1 = fromC ;
      if (class(fromC->tag) || (fromC->tag && fromC->tag < _Date ))
	tag1 = 0 ;
      else
	{
	  tag1 = fromC->tag ;
	  if (lexword2key (c1->tagp, &tag2, _VSystem))
	    tag1 = tag2 ;
	}
      fromC = arrayp(spread->colonnes, fromC->from - 1, COL) ; 
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
  spread->modified = TRUE ;
  spreadDefineColonne(FALSE) ;

  return;
} /* spreadChooseTag */

/*********Tag Chooser 2 *******************/

static void spreadEditTagText(void)
{
  COL *c ;
  char *cp, *cq ;
  SPREAD spread = currentSpread("spreadEditTagText") ;

  c = arrayp(spread->colonnes, spread->activeColonne, COL) ;  
  cp = c->tagTextBuffer ;
  if (!*cp)
    strcpy(cp, "?")  ;
  c->tagStack = stackReCreate(c->tagStack, 32) ;
  cq = cp + strlen(cp) -1 ;
  while (cq >= cp && *cq == ' ')
    *cq-- = 0 ;
  pushText(c->tagStack, cp) ;
  c->tagp = stackText(c->tagStack, 0) ;
  spread->modified = TRUE ;
  spreadDefineColonne(FALSE) ;

  return;
} /* spreadEditTagText */

/******** TagChooser 3 *****  Less good than the other method

static void spreadDefineTag (KEY k, int box)
{
  COL *c = spreadDefineActivate(box) ;
  KEY old = c->tag ;
  SPREAD spread = currentSpread("spreadDefineTag") ;

  c->tag = k ;
  c->tagp = freekey2text(k, arrp(c->tagMenu, 0, FREEOPT)) ;
  c->type = 'b' ;
  c->typep =  freekey2text(c->type, typeChoice) ;
  stackDestroy(c->tagStack) ;
  graphBoxDraw(box, BLACK, LIGHTBLUE) ;
  if (k != old)
    spread->modified = TRUE ;
}



static Array spreadDefineTagMenu (int classe)
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

static void spreadDefineClass (KEY k, int box)
{
  COL *c = spreadDefineActivate(box) ;
  int old = c->classe ; 
  SPREAD spread = currentSpread("spreadDefineClass") ;

  c->type = c->realType = 'k' ;
  c->nonLocal = FALSE ;
  c->classe = k ;
  c->classp = name(k) ;
  
  spreadDefineColonne(FALSE) ;
  if (k != old)
    spread->modified = TRUE ;

  return;
} /* spreadDefineClass */

/************************************************/

static void spreadDefineFrom (void)
{ 
  int level ;
  COL *col ;
  int old ;
  SPREAD spread = currentSpread("spreadDefineFrom") ;
  
  col = arrayp(spread->colonnes, spread->activeColonne, COL) ;  
  old = col->from ;
  if (!spread->activeColonne)
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
	col->from >= arrayMax(spread->colonnes))
      { messout("Please build from an existing colonne") ;
      col->from = 0 ;
      }
    spreadDefineColonne(TRUE) ;
    if (col->from != old)
      spread->modified = TRUE ;
    }
}

static BOOL spreadDefineWidth (int n)
{ 
  return (n > 0) ;
}

static BOOL spreadDefineDnaCoord (char *text, int box)
{ 
	return TRUE ;
}

/*
static BOOL spreadDefinePepCoord (char *text, int box)
{ 
	return TRUE ;
}
*/
/************************************************/
/*
static void spreadGetWidths (SPREAD spread)
{ int i, max = arrayMax(spread->colonnes), level ;
  COL *c ;
  
  for (i = 0; i < max ; i++)
    { c = arrayp(spread->colonnes, i, COL) ;  
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

static void spreadDefineCondition (void)
{
  COL *c ;
  char *cp ;
  SPREAD spread = currentSpread("spreadDefineCondition") ;

  c = arrayp(spread->colonnes, spread->activeColonne, COL) ;  
  cp = c->conditionBuffer ;
  while (*cp) cp ++ ;
  while (--cp >= c->conditionBuffer && *cp == ' ')
    *cp = 0 ;
  if (*c->conditionBuffer &&
      !condCheckSyntax(messprintf(" %s", freeprotect(c->conditionBuffer))))
    messout("Please correct this syntax error") ;
  spread->modified = TRUE ;

  return;
} /* spreadDefineCondition */
  
/************************************************/

static void spreadDefineSubtitle (void)
{
  COL *c ;
  char *cp ;
  SPREAD spread = currentSpread("spreadDefineSubtitle") ;

  c = arrayp(spread->colonnes, spread->activeColonne, COL) ;  
  cp = c->subtitleBuffer ;
  while (*cp) cp ++ ;
  while (--cp >= c->subtitleBuffer && *cp == ' ')
    *cp = 0 ;
  spread->modified = TRUE ;

  return;
} /* spreadDefineSubtitle */
  
/************************************************/

static void spreadDefineLegend (void)
{
  COL *c ;
  char *cp ;
  SPREAD spread = currentSpread("spreadDefineSubtitle") ;

  c = arrayp(spread->colonnes, spread->activeColonne, COL) ;  
  cp = c->legendBuffer ;
  while (*cp) cp ++ ;
  while (--cp >= c->legendBuffer && *cp == ' ')
    *cp = 0 ;
  spread->modified = TRUE ;

  return;
} /* spreadDefineLegend */
  
/************************************************/

static void spreadDefineTitle (void)
{ 
  char *cp ;
  SPREAD spread = currentSpread("spreadDefineTitle") ;

  cp = spread->titleBuffer ;
  while (*cp) cp ++ ;
  while (--cp >= spread->titleBuffer && *cp == ' ')
    *cp = 0 ;
  spread->modified = TRUE ;

  return;
} /* spreadDefineTitle */
  
/************************************************/

static void spreadDefineParam (void)
{ 
  char *cp ;
  SPREAD spread = currentSpread("spreadDefineParam") ;

  cp = spread->paramBuffer ;
  while (*cp) cp ++ ;
  while (--cp >= spread->paramBuffer && *cp == ' ')
    *cp = 0 ;
  spread->modified = TRUE ;

  return;
} /* spreadDefineParam */
  
/************************************************/

void spreadShow(void)
     /* private within spreadDisp-package */
{ 
  SPREAD spread = currentSpread("spreadShow") ;

  if (!graphCheckEditors (graphActive(), 0))
    return ;
  /*  spreadGetWidths (spread) ;*/
  spread->modified = TRUE ;	/* force recompute */
  if (spreadDoRecompute (spread))
    spreadDisplayData (spread) ;

  return;
} /* spreadShow */

/***********/

static void spreadForceOpenColonne (BOOL force)
{
  SPREAD spread = currentSpread("spreadForceOpenColonne2") ;

  spreadInitColonne(spread) ;
  spreadDefineColonne(force) ;
  spread->modified = TRUE ;

  return;
} /* spreadForceOpenColonne */

static void spreadOpenColonne (void)
{
  SPREAD spread = currentSpread("spreadOpenColonne") ;

  spreadInitColonne(spread) ;
  spreadDefineColonne(FALSE) ;
  spread->modified = TRUE ;

  return;
} /* spreadOpenColonne */

/***********/

static BOOL spreadShiftParameters (char *buffer, int x0, int dx)
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
} /* spreadShiftParameters */

/***********/

static void spreadInsertColonne (void)
{
  COL *c ;
  int i, nn ;
  SPREAD spread = currentSpread("spreadInsertColonne") ;

  nn = arrayMax(spread->colonnes) ;
  if (! nn)
    return ;

  c = arrayp(spread->colonnes, spread->activeColonne, COL) ;  

  if (!messQuery(messprintf("Do you really want to insert a colonne after colonne %d",
		 spread->activeColonne + 1)))
    return ;

  c = arrayp(spread->colonnes, nn, COL) ;   /* create an extra colonne at the end */
  for (i = nn ; i >  spread->activeColonne + 1 ; i--)
    {
      c = arrayp(spread->colonnes, i, COL) ;
      *c = *(c-1) ;
    }
  spreadDoInitColonne (spread, spread->activeColonne + 1) ;

  for (i = spread->activeColonne + 2 ; i < arrayMax(spread->colonnes) ; i++)
    {
      c = arrayp(spread->colonnes, i, COL) ;  
      if (c->from > spread->activeColonne + 1)
	c->from++ ;
      spreadShiftParameters (c->conditionBuffer, spread->activeColonne + 1, 1) ;
      spreadShiftParameters (c->dna1Buffer, spread->activeColonne + 1, 1) ;
      spreadShiftParameters (c->dna2Buffer, spread->activeColonne + 1, 1) ;
       if (stackExists(c->tagStack))
	 {
	   spreadShiftParameters (c->tagTextBuffer, spread->activeColonne + 1, 1) ;
	   stackClear (c->tagStack) ;
	   pushText (c->tagStack, c->tagTextBuffer) ;
	 }
    }
  isCheckEditor = TRUE ;
  
  spread->modified = TRUE ;
  spreadDefineColonne(FALSE) ;
}

/***********/

static void spreadRemoveColonne (void)
{ COL *c, *c1 ;
  int i ;
  SPREAD spread = currentSpread("spreadRemoveColonne") ;

  if (!arrayMax(spread->colonnes))
    return ;
  c = arrayp(spread->colonnes, spread->activeColonne, COL) ;  
  if (!spread->activeColonne)
    { messout ("You cannot remove the first colonne") ;
      return ;
    }
  if (!messQuery(messprintf("Do you really want to remove colonne %d",
		 spread->activeColonne + 1)))
    return ;
  spreadDestroyCol(c) ;
  for ( i = spread->activeColonne, c = arrayp(spread->colonnes, i, COL) ;  
        i + 1 < arrayMax(spread->colonnes) ; i++)
    { c1 = c++ ;
      memcpy(c1, c, sizeof(COL)) ;
    }

  arrayMax(spread->colonnes) -- ;

  /*    for (i = spread->activeColonne + 1 ; i < arrayMax(spread->colonnes) ; i++)
    { c = arrayp(spread->colonnes, i, COL) ;  
      if (c->from > spread->activeColonne)
	c->from -- ;
      else if (c->from == spread->activeColonne)
 messout ("You should modify colonne %d, which is derived from the colonne you just destroyed", i + 1) ;
    }*/
  for (i = spread->activeColonne ; i < arrayMax(spread->colonnes) ; i++)
    {
      c = arrayp(spread->colonnes, i, COL) ;  
      if (c->from > spread->activeColonne + 1)
	c->from -- ;
      else if (c->from == spread->activeColonne + 1)
	messout ("You should modify colonne %d, which is derived from the colonne you just destroyed", i + 2) ;
      spreadShiftParameters (c->conditionBuffer, spread->activeColonne, -1) ;
      spreadShiftParameters (c->dna1Buffer, spread->activeColonne, -1) ;
      spreadShiftParameters (c->dna2Buffer, spread->activeColonne, -1) ;
    }
  isCheckEditor = TRUE ;
  /*     spread->activeColonne-- ; mhmp 21.04.98  12.05.99*/
  spread->activeColonne = arrayMax(spread->colonnes) - 1 ;
  spread->modified = TRUE ;
  spreadDefineColonne(FALSE) ;
}

/**********************************************************/

static void shouldWeDestroy(void)
{ SPREAD spread = currentSpread("shouldwedestroy") ;

  if (spread->quitWithoutConfirmation ||
      messQuery("Do you really want to quit the Table_Maker ?"))
    graphDestroy() ;
}

/***********************************************************/
/*********** Input Output of the definitions ***************/ 
/* The actual operation are in a separate non graphic file */


static char dirName[DIR_BUFFER_SIZE], fileName[FIL_BUFFER_SIZE] ;
static BOOL firstDirPass = TRUE ; 

static void spreadDoGetDefinitions (KEY key)
{ SPREAD spread = currentSpread("spreadDoGetDefinitions") ;

  displayUnBlock () ;
  if (class(key) != _VTable)
   messout ("Sorry, you must pick a Table object") ;
  else if (!iskey(key))
    messout ("Sorry, this object is empty.") ;
  else
    {
      spreadDoReadDefinitions (spread, key, 0, 0,"", TRUE) ; /* do not substitute param */
  /* mhmp 12.05.99*/
      spread->activeColonne = arrayMax(spread->colonnes) - 1 ;
      spreadDefineColonne(TRUE) ;  /* To redraw */
      spread->fileName = 0 ;
    }
}

static void spreadGetDefinitions (void)
{ FILE *f ;
  char *cp ;
  SPREAD spread = currentSpread("spreadGetDefinitions") ;
  Graph old = graphActive () ;
  KEYSET ks = query (0,">?Table") ;
  cp = spread->fileName ;
  if( cp && *cp && strlen(cp) && spread->modif
      && messQuery (messprintf("%s not saved, Save ?", cp)) )
    { 
      f = filqueryopen(spread->dirName, spread->fileName, "def", "w",
			"Choose a File to store the Table Definition") ;
      if (!f)
	messout("Sorry, not done") ;
      else
	spreadDoSaveDefinitions (spread, f) ;
    }  
  displayCreate(DtKeySet) ;
  graphRetitle("Table definitions") ;
  keySetShow (ks,0) ;
  keySetSelect () ;
  graphActivate (old) ;

  displayBlock (spreadDoGetDefinitions,
		"Pick a Table object") ;
}

static void spreadReadDefinitions (void)
{ FILE *f ;
  char *cp ;/*mhmp 11.98*/

  SPREAD spread = currentSpread("spreadReadDefinitions") ;
  
  if (firstDirPass)
    { if (filName ("wquery", 0, "r"))
	strncpy (dirName, filName ("wquery", 0, "r"), DIR_BUFFER_SIZE-1) ;
      strcpy (fileName, "table") ;
      firstDirPass = FALSE ;
    }
  cp = spread->fileName ; /* mhmp 26.11.98 pour saver le fichier en cours */
  if( cp && *cp && strlen(cp) && spread->modif
      && messQuery (messprintf("%s not saved, Save ?", cp)) )
    { 
      f = filqueryopen(spread->dirName, spread->fileName, "def", "w",
			"Choose a File to store the Table Definition") ;
      if (!f)
	messout("Sorry, not done") ;
      else
	spreadDoSaveDefinitions (spread, f) ;
      spread->modif = FALSE ; 
      return ;
    }  

  f = filqueryopen (dirName, fileName, "def","r", 
		    "Choose a Table-definition file") ;
  if (!f)
    return ;

  spread->modif = FALSE ;
  spread->fileName = fileName ;
  spread->dirName = dirName ;
  if (!spreadDoReadDefinitions (spread, 0, f, 0,"", TRUE))
 /*  do not substitute param */
    {
      messout ("Sorry, your definition file is not correct") ;
	spreadDestroy (spread) ;
	spreadDispCreate (TRUE) ;
	return ;
    }
  /* mhmp 12.05.99*/
  spread->activeColonne = arrayMax(spread->colonnes) - 1 ;
  spreadDefineColonne(TRUE) ;  /* To redraw */
  if (!graphCheckEditors (graphActive(), 0))
    {
      messout ("Sorry, your definition file is not correct") ;
      spreadDestroy (spread) ;
      spreadDispCreate (TRUE) ;
    }
  /* repeat after eventual destroy */
  spread = currentSpread("spreadReadDefinitions") ;
  spread->modif = FALSE ;
  spread->fileName = fileName ;
  spread->dirName = dirName ;
}

static void spreadWriteDefinitions (void)
{ FILE *f = 0 ;
  COL *c = 0 ;
  int  maxCol ;
  
  SPREAD spread = currentSpread("spreadWriteDefinitions") ;

 
  if (!graphCheckEditors (graphActive(), 0))
      return ; 
  maxCol = arrayMax (spread->colonnes) ;
  if (maxCol)
    c = arrp(spread->colonnes, 0 , COL) ;
   /* if(!spread->modif) return ; mieg: I may want to save an unmodified file */
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
  spread->modif = FALSE ;
  spreadDoSaveDefinitions (spread, f) ;
}

static void spreadSaveDefinitions (void)
{ static KEY tableKey = 0 ;
  char *nam ;
  SPREAD spread = currentSpread("spreadSaveDefinitions") ;
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
  spreadDoSaveInObj (spread, tableKey) ;
abort: 
  messfree (nam) ;
}

/************/

static void spreadExportBql (void)
{
  COL *c = 0 ;
  int  maxCol ;
  SPREAD spread = currentSpread("spreadEportBql") ;

 
  if (!graphCheckEditors (graphActive(), 0))
      return ; 
  maxCol = arrayMax (spread->colonnes) ;
  if (maxCol)
    c = arrp(spread->colonnes, 0 , COL) ;
   /* if(!spread->modif) return ; mieg: I may want to save an unmodified file */
  if (!c || !c->type)
    return ; 
 
  spreadDoExportBql (spread, 0, TRUE) ;
} /* spreadExportBql */

/************/

static void spreadGraphTop (void)
{
  graphGoto (1,1) ;
}
/************/

static void spreadNewTable (void)
{
  spreadDispCreate (0) ; 
}

/************/

extern void spreadImportKeySet (void) ;
/* static void spreadComment(void) ; */
static void spreadCallImportKeySet (void) ;

static MENUOPT spreadDefMenu[] =
  {
   {shouldWeDestroy, "Quit"},
   {help, "Help"},
   {graphPrint, "Print"},
   {spreadShow, "Search Whole Class"},
   {spreadCallImportKeySet, "Search Active KeySet"},
   {spreadGetDefinitions, "Standard Query"},
   {spreadSaveDefinitions, "Save Standard Query"},/* mhmp 09.07.02  write save*/
   {spreadReadDefinitions, "Read Query from file"},
   {spreadWriteDefinitions, "Write Query to file"},/* mhmp 09.07.02  write save*/
   {spreadExportBql, "Write Bql"},
   {spreadOpenColonne, "Add Column"}, 
   {spreadInsertColonne, "Insert Column"}, 
   {spreadRemoveColonne, "Suppress Column"}, 
   {spreadNewTable, "New table"},
/*    {spreadComment, "Comments"}, Not very useful */
   {0, 0}
   } ;

/**********************************************************/
static void spreadCallImportKeySet (void)
{
  if (!graphCheckEditors (graphActive(), 0))
      return ;
  spreadImportKeySet () ;
} 
/*
static void spreadComment(void)
{  int x = 1 , y = 7, maxx = 30  ;
   Stack s ; char *cp ;
  SPREAD spread = currentSpread("spreadComment") ;

  s = spread->comments ;

  graphClear() ;
  graphHelp("Table_Maker") ;
  graphMenu(spreadDefMenu) ;
  graphButtons (spreadDefMenu, 1, 1, 68) ;

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
  SPREAD spread = currentSpread("setSortColumn") ;

  freeforcecard (buf) ;
  if (!freecheck ("iz") ||
      !freeint (&i) || 
      i < 0 || 
      i > arrayMax(spread->colonnes))
    messout ("Sorry, value must be between a valid "
	     "column number or 0, which means left "
	     "to right ordering.") ;
  spread->sortColonne = i ;
}
*/
static BOOL setSortColumn (int n)
{
  SPREAD spread = currentSpread("setSortColumn") ;
  if (n < 0 || n > arrayMax(spread->colonnes))
    {
      messout ("Sorry, value must be between a valid "
	       "column number or 0, which means left "
	       "to right ordering.") ;
      return FALSE ;
    }
  return TRUE ;
}

/*********************************************/

void spreadDefineColonne (BOOL force)
{ 
  COL *c ;
  int i, max, box, from ;
  float line, centralLine = 1 ;
  SPREAD spread = currentSpread("spreadDefineColonne") ;

  if (!force && !isCheckEditor && !graphCheckEditors (graphActive(), 0))
      return ;
  isCheckEditor = FALSE ;
  if (!arrayMax(spread->colonnes))
    { spreadForceOpenColonne(force) ;
      return ;
    }
  classeListeInit () ;
  graphClear() ;
  spread->activeDefBox = 0 ;
  spread->defBoxPerCol = 0 ;
  graphRegister (PICK, spreadPickDefinition) ;
  graphHelp("Table_Maker") ;
  graphMenu(spreadDefMenu) ;

  box = graphBoxStart() ;
  graphButtons (spreadDefMenu, 1, 1, 68) ;
  graphBoxEnd() ;
  graphBoxDim (box, 0, 0, 0, &line) ;
  line += 1 ;

  graphText ("Sort column:", 3, line) ;
  strncpy (spread->sortBuffer, 
	   messprintf("%d", spread->sortColonne), 6) ;
  spread->sortColonneBox =
    /*    graphTextEntry (spread->sortBuffer, 6, 16, line, 
		    setSortColumn) ;*/
	graphIntEditor ("", &spread->sortColonne, 16, line, setSortColumn) ;
  graphText ("F4 to interrupt", 40, line) ;
  line += 1.5 ;

  graphText ("Title", 3, line) ;
  spread->titleBox = 
    graphTextScrollEntry (spread->titleBuffer, 300, 60, 9, line, 
		    TEXT_ENTRY_FUNC_CAST spreadDefineTitle) ;
  line += 1.5 ;

  graphText ("Parameters", 3, line) ;
  spread->paramBox = 
    graphTextEntry (spread->paramBuffer, 179, 14, line, 
		    TEXT_ENTRY_FUNC_CAST spreadDefineParam) ;
  line += 1.5 ;

  spread->definitionBox = graphBoxStart() ;
  graphText ("Column", 1, line) ;
  line += 1.5 ;

  max = arrayMax(spread->colonnes) ;
  for (i = 0, c = arrayp(spread->colonnes, 0, COL) ;
       i < max ; i++, c++, line += 2.5 )
    {
      if (i == spread->activeColonne)
	centralLine = line ;
      graphText(messprintf("%3d", i + 1), 2, line) ;
 
      graphText ("Title", 7, line) ;
      spread->subTitleBox = 
	graphTextEntry (c->subtitleBuffer, 59, 17, line, 
			TEXT_ENTRY_FUNC_CAST spreadDefineSubtitle) ;

      graphText ("Legend", 7, line+=1.3) ;
      spread->legendBox = 
	graphTextScrollEntry (c->legendBuffer, 1023, 59, 17, line, 
			TEXT_ENTRY_FUNC_CAST spreadDefineLegend) ;
  
      graphText ("Width", 7, line += 1.3) ;
      sprintf(c->widthBuffer,"%d", c->width) ;
      spread->widthBox = 
	/*	graphTextEntry (c->widthBuffer, 4, 14, line, 
			TEXT_ENTRY_FUNC_CAST spreadDefineWidth) ;mhmp 30.03*/
	graphIntEditor ("", &c->width, 14, line, spreadDefineWidth) ;
  
      box = spread->hideBox = graphBoxStart() ;
      graphTextPtrPtr(&c->hiddenp, 24, line, 7) ;
      c->hiddenp = freekey2text(c->hidden, hideChoice) ;
      graphBoxEnd() ;
      graphBoxBox(box) ;
      graphBoxFreeMenu(box, (FreeMenuFunction) spreadDefineHide,
		       hideChoice) ;
      
      box = spread->optionalBox = graphBoxStart() ;
      graphTextPtrPtr(&c->optionalp, 33, line, 9) ;
      c->optionalp = freekey2text(c->mandatory, optionalChoice) ;
      graphBoxEnd() ;
      if (i) 
	{ graphBoxBox(box) ;
	  graphBoxFreeMenu(box, (FreeMenuFunction) spreadDefineOptional,
			   optionalChoice) ;
	}
      else
	graphBoxMarkAsClear(box) ;
      box = spread->typeBox = graphBoxStart() ;
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
	graphBoxFreeMenu(box, (FreeMenuFunction) spreadDefineClass,
		       classeListe) ;
      else
	graphBoxFreeMenu(box, (FreeMenuFunction) spreadDefineType,
		       typeChoice2) ;

      box = spread->extendBox = graphBoxStart() ;
      graphTextPtrPtr(&c->extendp, 9, line += 1.5, 8) ;
      c->extendp = freekey2text(c->extend, extendChoice) ;
      graphBoxEnd() ;
      if (i) graphBoxBox(box) ;
      if (i)
	graphBoxFreeMenu(box, (FreeMenuFunction) spreadDefineExtend,
			 extendChoice) ;
      else
	graphBoxMarkAsClear(box) ;
      if (c->from <= 0) c->from = 1 ;
      sprintf(c->fromBuffer,"%d", from = c->from) ;
      spread->fromBox = box =
	graphTextEntry (c->fromBuffer, 4, 18, line, 
			TEXT_ENTRY_FUNC_CAST spreadDefineFrom) ;
      /* intEditor does not allow to call a full redraw */

      if (!i)
	graphBoxMarkAsClear(box) ;      
      box = spread->dnaBox = graphBoxStart() ;
      graphText ("DNA" , 28, line) ;
      graphBoxEnd() ; 
       if (i &&
	   c->extend == 'f' &&
	   c->from <= arrayMax (spread->colonnes) &&
	   arr(spread->colonnes, c->from - 1 , COL).classe &&
	   dnaInClass (superClass( arr(spread->colonnes, c->from - 1 , COL).classe))
	  )  { graphBoxBox(box) ; }
      else graphBoxMarkAsClear(box) ;
      box = spread->pepBox = graphBoxStart() ;
      graphText ("PEP" , 28, line) ;
      graphBoxEnd() ; 
      if (i &&
	  c->extend == 'f' &&
	  c->from <= arrayMax (spread->colonnes) &&
	  arr(spread->colonnes, c->from - 1 , COL).classe &&
	  pepInClass (superClass( arr(spread->colonnes, c->from - 1 , COL).classe))
	  )  { graphBoxBox(box) ; }
      else graphBoxMarkAsClear(box) ;
      box = spread->tagBox = graphBoxStart() ;
      if (i) graphText ("Tag:" , 33, line) ;
      graphBoxEnd() ;
      if (i)  { graphBoxBox(box) ; } /* {} needed around stupid macro bobox */
      else graphBoxMarkAsClear(box) ;

      /* a left over which seems useless as of  mars 22 2002

      {
      COL * fromC ;
      if (from > 0)
	fromC = arrp(spread->colonnes, from - 1, COL) ;
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
	spread->dna1Box = box =
	  graphTextScrollEditor ("from", c->dna1Buffer, 200, 10, 38, line, spreadDefineDnaCoord) ;
      else
	{
	  graphBoxStart () ; graphBoxEnd() ; 
	  graphBoxStart () ; graphBoxEnd() ; 
	}
      if (c->type == 'D' || c->type == 'P')
	spread->dna2Box = box =
		graphTextScrollEditor ("to", c->dna2Buffer, 200, 10, 54, line, spreadDefineDnaCoord) ; /* mhmp 11.07.02 from --> to */
      else
	{
	  box = graphBoxStart () ; graphBoxEnd() ; 
	  box = graphBoxStart () ; graphBoxEnd() ; 
	}
      box = graphBoxStart () ; graphBoxEnd() ; 
      if (i && c->type != 'D'&& c->type != 'P')
	{
	  spread->tagTextBox = box =
	    graphTextScrollEntry(c->tagTextBuffer, 359, 
				 300, 38, line, 
				 TEXT_ENTRY_FUNC_CAST spreadEditTagText) ;
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
      spread->conditionBox = 
	graphTextScrollEntry (c->conditionBuffer, 359, 
			      300, 17, line, 
			TEXT_ENTRY_FUNC_CAST spreadDefineCondition) ;
      if (c->type == 'D' || c->type == 'P')
	graphText ("limit to sequence matching a UNIX \"RegExp\", i.e: atgc  ^g[tc].*ag$  ", 7, line += 1.1) ;

          /* last paragraph of the loop */
      box = graphBoxStart() ;
      graphBoxEnd() ;
      if (!spread->defBoxPerCol)
	spread->defBoxPerCol = box - spread->definitionBox ;
    }
  graphBoxEnd() ;
  graphButton("Add column", spreadOpenColonne, 2, line) ;
  graphButton("Page top", spreadGraphTop, 15, line) ;
  graphButton("Search Whole Class", spreadShow, 26, line) ;
  graphButton("Search Active KeySet", spreadCallImportKeySet, 48, line) ;
  graphTextBounds (300, line + 3) ;

  graphRedraw() ;
  graphGoto (1, centralLine) ;
  return ;
}

/***********************************************************/
/***********************************************************/
 
 
