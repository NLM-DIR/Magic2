/*  File: sprdop.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 11 23:43 1999 (rd)
 * Created: Fri Dec 11 13:53:37 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: sprdop.c,v 1.54 2020/01/12 12:58:02 mieg Exp $ */


#include "acedb.h"

#include "bitset.h"
#include "spread_.h"
#include "vtxt.h"
#include "freeout.h"
#include "lex.h"
#include "pick.h"
#include "systags.h"
#include "tags.h"
#include "table.h"
#include "sysclass.h"
#include "query.h"
#include "bs.h"
#include "dna.h"
#include "peptide.h"
#include "java.h"
#include "client.h"


/************************************************************/

static Array colonnes = 0 ;
Array pos2col = 0 ;
static void spreadDestroyConditions (SPREAD spread);
static void spreadStorePrecomputation (SPREAD spread)  ;
static void spreadTable2Tableau (SPREAD spread, TABLE *tt) ;
static int scLoopRight (SPREAD spread, SC sc) ;

TABLE* (*externalTableMaker) (char *quer) = 0 ; /* allocation */

/*****************************************************/
/************* Public utilities **********************/
/*****************************************************/

magic_t SPREAD_MAGIC = "SPREAD_MAGIC";

/*****************************************************/

SPREAD spreadCreate (void)
     /* public function */
{
  SPREAD spread ;

  spread = (SPREAD) messalloc (sizeof (struct SpreadStruct));

  spread->magic = &SPREAD_MAGIC; 
  spread->colonnes = arrayCreate (256, COL) ;
  spread->tableau = arrayCreate (10, Array) ;
  spread->pos2col = arrayCreate (8, int) ;
  spread->modif = FALSE ;  /* mhmp24.11.98 */

  return spread ;
} /* spreadCreate */

/*****************************************/

void uSpreadDestroy (SPREAD spread)
     /* public access via spreadDestroy macro */
{    
  Array t , c ;
  int i ;
  char * cp ; 
    
  if (!spread)
    return;

  if (spread->magic != &SPREAD_MAGIC)
    messcrash ("uSpreadDestroy received corrupt spread->magic");

  if (spread->precomputing && !spread->table)
    spreadStorePrecomputation (spread) ;
  
  tableDestroy (spread->table) ;
  spread->magic = 0 ;
  cp =  spread->fileName ;
/* mhmp 25.11.98 */
  if ( cp && *cp && strlen (cp) && spread->modif
      && messQuery (messprintf ("%s not saved, Save ?", cp)) )
    {
      FILE *f = filqueryopen (spread->dirName, cp, "def", "w",
		       "Choose a File to store the Table Definition") ;
      if (!f)
	messout ("Sorry, not done") ;
      else
	spreadDoSaveDefinitions (spread, f) ;
    }
  
  spreadDestroyConditions (spread) ;
  
  if ((t = spread->tableau))
    { i = arrayMax (t) ;
      while (i--) 
	arrayDestroy (arr (t,i,Array)) ;
      arrayDestroy (t) ;
    }
  arrayDestroy (spread->flags) ;
  
  if ((c = spread->colonnes))
    { i = arrayMax (c) ;
      while (i--)
	spreadDestroyCol (arrp (c,i,COL)) ;
      arrayDestroy (spread->colonnes) ;
    }
  arrayDestroy (spread->pos2col) ;
  arrayDestroy (spread->mapBoxes) ;
  arrayDestroy (spread->oldClass) ;
  stackDestroy (spread->comments) ;
  keySetDestroy (spread->exportedKeySet) ;
  messfree (spread->href) ;
  messfree (spread) ;
  spread = 0;
} /* uSpreadDestroy */

/*****************************************/
/* clear the DNA and other long arrays already exported */
void spreadCleanUp (SPREAD spread)
     /* public access via spreadDestroy macro */
{    
  Array cols ;
  COL* c ;
  int i ;
    
  if (!spread)
    return;

  if (spread->magic != &SPREAD_MAGIC)
    messcrash ("spreadCleanUp received corrupt spread->magic");

  
  if ((cols = spread->colonnes))
    { i = arrayMax (cols) ;
      while (i--)
	{
	  c = arrp (cols,i,COL) ;
	  stackDestroy (c->text) ;
	  stackDestroy (c->dnaStack) ;
	  stackDestroy (c->pepStack) ;
	  arrayDestroy (c->dnaD) ;
	  arrayDestroy (c->dnaR) ;
	  arrayDestroy (c->pep) ;
	}
    }
} /* spreadCleanUp */

/*****************************************************/

void spreadDestroyCol (COL* c)
     /* private within SpreadPackage */
{
  stackDestroy (c->text) ;
  stackDestroy (c->dnaStack) ;
  stackDestroy (c->pepStack) ;
  arrayDestroy (c->dnaD) ;
  arrayDestroy (c->dnaR) ;
  arrayDestroy (c->pep) ;
  stackDestroy (c->tagStack) ;
  arrayDestroy (c->tagMenu) ;
  arrayDestroy (c->segs) ;
  arrayDestroy (c->segps) ;
  bumpDestroy (c->bump) ;
  condDestroy (c->conditionCondition) ;  

  return;
} /* spreadDestroyCol */

/*****************************************************/

static void spreadDestroyConditions (SPREAD spread)
{
  COL *colonne = arrp (spread->colonnes, 0, COL) - 1 ;

  int max = arrayMax (spread->colonnes)  ;
  while (colonne++, max--)
    { 
      condDestroy (colonne->tagCondition) ;
      condDestroy (colonne->countCondition) ;
      condDestroy (colonne->conditionCondition) ;
      queryRegExpDestroy (colonne->dnaCondition) ;
    }

  return;
} /* spreadDestroyConditions */

/*****************************************/

static void spreadCleanUpTableau (SPREAD spread)
{
  int i = arrayMax (spread->colonnes) ;
  COL *col ;

  while (i--)
    { col = arrp (spread->colonnes, i, COL) ;
      if (col->text) 
	{ col->text = stackReCreate (col->text, 100) ;
	  pushText (col->text, "toto") ; /* avoid zero */
	}
    }
      
  if (arrayExists (spread->tableau))
    { i = arrayMax (spread->tableau) ;
      while (i--)
	arrayDestroy (arr (spread->tableau, i, Array)) ;
    }
  spread->tableau = arrayReCreate (spread->tableau, 50, Array) ;

  return;
} /* spreadCleanUpTableau */

/*****************************************************/
/* a good function i think
static void spreadSaveNthCol (SPREAD spread, KEYSET ks, int col)
{
  Array t = spread->tableau ;
  int i , j = 0 , max = arrayMax (t) ;
  KEY key, lastKey = 0 ;

  ks = keySetReCreate (ks) ;
  for (i = 0 ; i < max ; i++)
    { key = arr (arr (t,i,Array), col, SPCELL).u.k ;
      if (key && key != lastKey)
	keySet (ks, j++) = lastKey = key ;
    }
  keySetSort (ks) ;
  keySetCompress (ks) ;
}
*/

/**************************************************************/
/**************************************************************/

         /* self managed calloc */

static AC_HANDLE scHandle = 0 ;
static Stack freeScStack = 0 ;
static int nScUsed = 0 , nScAlloc = 0 ;    /* useful for debug */

static SC scAlloc (void)       /* self managed calloc */
{
  static int blocSize = 256 ; /* i need no more SC than columns */
  SC p = 0 ;
  int i ;

  if (!freeScStack)
    {
      scHandle = ac_new_handle () ;
      nScUsed = nScAlloc = 0 ;
      blocSize = 256 ; 
      freeScStack = stackHandleCreate (sizeof (SC)*blocSize, scHandle) ;
    }
  if (stackEmpty (freeScStack))
    { p = (SC) halloc (blocSize * sizeof (struct spread_cell_struct), scHandle) ;
      for (i = blocSize ; i-- ; ++p)
        push (freeScStack,p,SC) ;
      nScAlloc += blocSize ;
      blocSize *= 2 ;
    }
  p = pop (freeScStack, SC) ;
  memset (p, 0, sizeof (struct spread_cell_struct)) ;
  nScUsed ++ ;

  return p ;
} /* scAlloc */

static void scFree (SC sc)
{
  if (sc && freeScStack)
    { 
      if (sc->obj)
	bsDestroy (sc->obj) ;
      if (sc->mark) messfree (sc->mark) ;
      push (freeScStack, sc, SC) ;
      nScUsed-- ;
    }
  else
    invokeDebugger () ;

  return;
} /* scFree */

void scShutDown (void)
{
  ac_free (scHandle) ;
}

void scStatus (int *used, int *alloc)
{ *used = nScUsed ; *alloc = nScAlloc ;
}

/*****************************************************/
/************* Sorting of lines **********************/
/*****************************************************/

static int sortColonne ;

int spreadOrder (const void *x, const void *y)
     /* private within SpreadPackage */
{
  const Array a = * (const Array*)x, b = * (const Array*)y ;
  int i, j, max = arrayMax (colonnes), r ;
  COL *c ;
  SPCELL *su, *sv ;
  BSunit u, v ;


  if (sortColonne >= 0)  /* if -1 sort all columns in natural order */
    {
      j = sortColonne ;
      su = arrp (a,j,SPCELL) ;
      sv = arrp (b,j,SPCELL) ;
      if (su->empty && ! sv->empty) return -1 ;
      if (sv->empty && ! su->empty) return  1 ;
      c = arrp (colonnes,j,COL) ;
      if (c->type /* && !c->hidden */) 	
				/* RD 990311 - comment out && !c->hidden to allow sorting on hidden columns */
				/* I see no reason why not and it is clearly useful and used */
	{ 
	  u = su->u ;
	  v = sv->u ;
	  
	  if (!arr (b,j,SPCELL).empty && !arr (a,j,SPCELL).empty && (u.k != v.k))
	    {
	      if (c->showType == SHOW_MULTI) 
		{
		  r = lexstrcmp
		    (stackText (c->text, u.i),
		     stackText (c->text, v.i)) ;
		  if (r)
		    return r ;
		}
	      else switch (c->type)
		{
		case 'b':
		  return u.k ? -1 : 1 ;
		case 'c':
		  return u.i > v.i ? 1 : -1 ;
		case 'i':
		case 'f':
		  return su->z > sv->z ? 1 : -1 ;
		case 'd':
		  return u.time > v.time ? 1 : -1 ;
		case 'k': case 'K': case 'n':
		  return keySetAlphaOrder (&u.k, &v.k) ;
		case 't':
		  r = lexstrcmp
		    (stackText (c->text, u.i),
		     stackText (c->text, v.i)) ;
		  if (r)
		    return r ;
		  break ;
		case 'D':
		  r = lexstrcmp
		    (stackText (c->dnaStack, u.k),
		     stackText (c->dnaStack, v.k)) ;
		  if (r)
		    return r ;
		  break ;
		case 'P':
		  r = lexstrcmp
		    (stackText (c->pepStack, u.k),
		     stackText (c->pepStack, v.k)) ;
		  if (r)
		    return r ;
		  break ;
		default:
		  messcrash ("Unknown type in spreadOrder") ;
		}
	    }
	}
    }

  for (i = 0 ; i < max ; i++)
    {
      j = arr (pos2col,i, int) ;
      if (j == sortColonne) continue ;
      c = arrp (colonnes,j,COL) ;
      if (c->type && !c->hidden)
	{ 
	  su = arrp (a,j,SPCELL) ;
	  sv = arrp (b,j,SPCELL) ;
	  u = su->u ;
	  v = sv->u ;

	  if (su->empty && sv->empty) continue ;
	  if (su->empty) return -1 ;
	  if (sv->empty) return  1 ;
	  if (u.k != v.k)
	    {
	      if (c->showType == SHOW_MULTI) 
		{
		  r = lexstrcmp
		    (stackText (c->text, u.i),
		     stackText (c->text, v.i)) ;
		  if (r)
		    return r ;
		}
	      else switch (c->type)
		{
		case 'b':
		  return u.k ? -1 : 1 ;
		case 'c':
		  return u.i > v.i ? 1 : -1 ;
		case 'i': 
		case 'f':
		  return su->z > sv->z ? 1 : -1 ;
		case 'd':
		  return u.time > v.time ? 1 : -1 ;
		case 'k': case 'K': case 'n': 
		  return keySetAlphaOrder (&u.k, &v.k) ;
		case 't':
		  r = lexstrcmp
		    (stackText (c->text, u.i),
		     stackText (c->text, v.i)) ;
		  if (r)
		    return r ;
		  break ;
		case 'D':
		  r = lexstrcmp
		    (stackText (c->dnaStack, u.k),
		     stackText (c->dnaStack, v.k)) ;
		  if (r)
		    return r ;
		  break ;
		case 'P':
		  r = lexstrcmp
		    (stackText (c->pepStack, u.k),
		     stackText (c->pepStack, v.k)) ;
		  if (r)
		    return r ;
		  break ;
		default:
		  messcrash ("Unknown type in spreadOrder") ;
		}
	    }
	}
    }
  return 0 ;
} /* spreadOrder */

/*********************/

void spreadReorder (SPREAD spread)
     /* private within SpreadPackage */
{
  int i, j, max = arrayMax (spread->colonnes) ;
  
  colonnes = spread->colonnes ;  /* needed by spreadOrder */
  pos2col = spread->pos2col ;
  sortColonne = spread->sortColonne - 1 ;
  
  /*  this will deconnect sortColonne once one switches them  */
  for (i = 0 ; i < max ; i++)
    { 
      j = arr (pos2col,i, int) ;
      if (i != j)
	sortColonne = - 1 ;
    }
  
  if (spread->tableau)
    arraySort (spread->tableau, spreadOrder) ;
  
  return;
} /* spreadReorder */

/*****************************************************/
/************* Syntax checking  **********************/
/*****************************************************/

static BOOL conditionHasParam (char *cp)
{
  if (cp) while (*cp)
    { if (*cp++ == '%')
	return TRUE ;
    }
  return FALSE ;
} /* conditionHasParam */

/********************************/

static BOOL spreadCheckOneCondition (SPREAD spread, char *buffer, int icol, BOOL isCalcul)
{ /* must check here that substitution is compatible with syntax */
  Stack  s ;
  int i, mine, n1, n2 ;
  char cc, *cp = buffer, *cq ;
  
  cp-- ;
  while (*++cp) 
    if (*cp == '%')
      {
	cc = * (cp + 1) ;
	if (cc == '%')
	  { messout ("Give a value to %% before running in column %d at:\n%s",
		    icol , cp) ;
	  return FALSE ;
	  }
	cq = cp + 1; i = 0 ;
	while (*cq >= '0' && *cq <= '9')
	  { i *= 10 ; i += (*cq - '0') ; cq++ ; }
	if (i <= 0 || i > icol)
	  { messout ("Cannot substitute parameter in column %d at:\n%s",
		    icol, cp) ;
	  return FALSE ;
	  }
      }
  
  s = stackCreate (500) ;
  if (isCalcul)
    catText (s,"[") ;
  catText (s, buffer) ;
  if (isCalcul)
    catText (s,"]") ;
  n1 = stackMark (s) ; pushText (s, " ") ;
  for (i = 0 ; i < icol ; i++)
    { /* use same value as in scMatch */
      switch (arrp (spread->colonnes, i, COL)->type)
	{
	case 'b':  case 'k': case'n': case 'K':
	  catText (s," __toto__ ") ;
	  break ;
	  case't':
	    catText (s," \"__please_protect__ __<text>__ __with_double_quotes__\"") ;
	  break ;
	case 'i': case 'f': case 'c':
	  catText (s," 1 ") ;
	  break ;
	case 'd':
	  catText (s," now ") ;
	  break ;
	case 'D':
	  catText (s," ATGC ") ;
	  break ;
	}
    }  
  mine = freesettext (stackText (s, 0), stackText (s,n1)) ;
  freespecial ("%\\\"") ;    /* No subshells ($ removed) */
  n2 = stackMark (s) ;
  while (freecard (mine))
    pushText (s, freepos ()) ;
  if (isCalcul)
    {
      if (!queryCalcul (stackText (s,n2), &n1, 0))
	{
	  messout (
		  "After test substitution, I cannot parse calcul of col %d :\n%s",
		  icol, stackText (s,n2)) ;
	  stackDestroy (s) ;
	  return FALSE ;
	}
    }
  else
    {
      if (!condCheckSyntax (stackText (s,n2)))
	{
	  messout (
		  "After test substitution, I cannot parse condition of col %d :\n%s",
		  icol, stackText (s,n2)) ;
	  stackDestroy (s) ;
	  return FALSE ;
	}
    }
  stackDestroy (s) ;
  return TRUE ;
}

/************************************/

BOOL spreadCheckConditions (SPREAD spread)
     /* private within SpreadPackage */
{
  COL *colonne = arrp (spread->colonnes, 0, COL) - 1 ;
  COND cond = 0 ;
  int from, icol = 0 , max = arrayMax (spread->colonnes) ;
  char cc ;

  spread->nCol = arrayMax (spread->colonnes) ;
  while (icol++, colonne++, max--)
    { 
      condDestroy (colonne->tagCondition) ;
      condDestroy (colonne->countCondition) ;
      condDestroy (colonne->conditionCondition) ;
       
      from =  colonne->from - 1 ;
      if (icol > 1)
	{
	  if (from + 1 >= icol || from < 0)
	    {  
	      messout ("From field %d of column %d is wrong",
		       from + 1, icol) ;
	      return FALSE ;
	    }
	  cc = arrp (spread->colonnes, from, COL)->type ;
	  if (colonne->extend == 'f' && 
	      cc != 'k' && cc != 'K' && cc != 'n')
	    { 	
	      messout (
		       "Column %d is computed from column %d which is not a list of objects",
		       icol, from + 1) ;
	      return FALSE ;
	    }
	  if (colonne->extend == 'f' && 
	      arrp (spread->colonnes, from, COL)->showType == SHOW_MULTI)
	    { 	
	      messout (
		       "Column %d is computed from column %d which is multi valued",
		       icol, from + 1) ;
	      return FALSE ;
	    }
	}

      if (colonne->type == 't' && !colonne->text)
	{ colonne->text = stackReCreate (colonne->text, 100) ;
	  pushText (colonne->text, "\177NULL") ; /* so 0 means undefined */
	 }

      arrayDestroy (colonne->pep) ;
      arrayDestroy (colonne->dnaD) ;
      arrayDestroy (colonne->dnaR) ;

      if (colonne->type == 'D')  /* DNA */
	{
	  if (!colonne->dnaStack)
	    { colonne->dnaStack = stackReCreate (colonne->dnaStack, 100) ;
	    pushText (colonne->dnaStack, "\177NULL") ; /* so 0 means undefined */
	    }
	  if (!spreadCheckOneCondition (spread, colonne->dna1Buffer, icol, TRUE) ||
	      !spreadCheckOneCondition (spread, colonne->dna2Buffer, icol, TRUE))
	    return FALSE ;
	}
      else if (colonne->type == 'P')  /* Peptide */
	{
	  if (!colonne->pepStack)
	    { colonne->pepStack = stackReCreate (colonne->pepStack, 100) ;
	    pushText (colonne->pepStack, "\177NULL") ; /* so 0 means undefined */
	    }
	  if (!spreadCheckOneCondition (spread, colonne->dna1Buffer, icol, TRUE) ||
	      !spreadCheckOneCondition (spread, colonne->dna2Buffer, icol, TRUE))
	    return FALSE ;
	}
      else if (colonne->showType == SHOW_COMPUTE)  /* compute */
	{
	  if (colonne->tagStack)
	    { 
	      if (!spreadCheckOneCondition (spread, stackText (colonne->tagStack,0), icol, TRUE) )
		return FALSE ;
	    }
	}
      else if (colonne->type == 'c')  /* Count */
	{
	  if (colonne->extend == 'c') ; /* COPY, no check needed on the tag */
	  else if (colonne->tagStack && stackText (colonne->tagStack, 0))
	    {
	      char *cp = strstr (stackText (colonne->tagStack,0), ";") ;

	      if (cp) *cp = 0 ;
	      if (!condConstruct (stackText (colonne->tagStack,0), &cond))
		{ messout ("syntax error in Tag field :%s",
			   stackText (colonne->tagStack, 0)) ;
		return FALSE ;
		}
	      condDestroy (colonne->tagCondition) ;
	      colonne->tagCondition = cond ;
	      if (cp)
		{
		  if (!condConstruct (cp + 1, &cond))
		    {
		      messout ("syntax error in COUNT Tag field :%s",
			       cp+1) ;
		      return FALSE ;
		    }   
		  condDestroy (colonne->countCondition) ;
		  colonne->countCondition = cond ;
		  *cp = ';' ;
		}
	    }
	  else
	    {
	      if (!colonne->tag)
		cond = 0 ;
	      else if (!condConstruct (name (colonne->tag), &cond))
		{ 
		  messout ("syntax error in Tag field :%s",
			   name (colonne->tag)) ;
		  return FALSE ;
		}
	      condDestroy (colonne->tagCondition) ;
	      colonne->tagCondition = cond ;
	    }
	}
      else  /* any other type */
	{
	  if (colonne->extend == 'c') ; /* COPY, no check needed on the tag */
	  else if (colonne->type == 'D' || colonne->type == 'P')  ;
	  else if (colonne->tagStack && stackText (colonne->tagStack, 0))
	    {
	      char *cp = strstr (stackText (colonne->tagStack,0), ";") ;

	      if (cp) *cp = 0 ;
	      if (!condConstruct (stackText (colonne->tagStack,0), &cond))
		{
		  messout ("syntax error in Tag field :%s",
			   stackText (colonne->tagStack, 0)) ;
		  return FALSE ;
		}
	      condDestroy (colonne->tagCondition) ;
	      colonne->tagCondition = cond ;
	      if (cp)
		{
		  if (!condConstruct (cp + 1, &cond))
		    {
		      messout ("syntax error in COUNT Tag field :%s",
			       cp+1) ;
		      return FALSE ;
		    } 
		  condDestroy (colonne->countCondition) ;
		  colonne->countCondition = cond ;
		  *cp = ';' ;
		}
	    }
	  else
	    {
	      if (!colonne->tag)
		cond = 0 ;
	      else if (!condConstruct (name (colonne->tag), &cond))
		{
		  messout ("syntax error in Tag field :%s",
			   name (colonne->tag)) ;
		  return FALSE ;
		}
	      condDestroy (colonne->tagCondition) ;
	      colonne->tagCondition = cond ;
	    }
	}
      
      colonne->condHasParam = conditionHasParam (colonne->conditionBuffer) ;
      condDestroy (colonne->conditionCondition) ;
      if (colonne->type == 'D' || colonne->type == 'P')
	{
	  if (colonne->conditionBuffer[0])
	    {
	      queryRegExpDestroy (colonne->dnaCondition) ;
	      colonne->dnaCondition = queryRegExpCreate (colonne->conditionBuffer, FALSE, 0) ;
	      if (! colonne->dnaCondition)
		return FALSE ;
	    }
	}
      else if (!colonne->condHasParam)
	{
	  if (*colonne->conditionBuffer)
	    {
	      /* allow condition on Text to unify syntax
               if ( ace_lower (colonne->type) != 'k' &&
		   colonne->type != 'n')
		{
		  messout ("Condition %s \n cannot be applied to a non-object",
			   colonne->conditionBuffer) ;
		  return FALSE ;
		}
	      else 
	      */
	      condDestroy (colonne->conditionCondition) ;
	      if (!condConstruct (colonne->conditionBuffer, &cond))
		{ messout ("syntax error in Tag field :%s",
			   colonne->conditionBuffer) ;
		return FALSE ;
		}
	      else
		colonne->conditionCondition = cond ;
	    }
	}
      else
	{
	  if (!spreadCheckOneCondition (spread, colonne->conditionBuffer, icol, FALSE))
	    return FALSE ; 
	  condDestroy (colonne->conditionCondition) ;
	  colonne->conditionCondition = 0 ; /* do not register in this case */
	}
    }

  return TRUE ;
} /* spreadCheckConditions */
  
/**********************************************************************/
/************************ Actual calculation **************************/
/**********************************************************************/

static Stack scStack = 0 ;

static BOOL scMatch (SPREAD spread, SC sc, OBJ *objp, KEY key)
{  
  COL *col = sc->col ; SC cc ;
  int mine, n0, n1, n2 ;
  long int li ;
  BOOL found = FALSE ;
  static Stack s = 0 ; char *cp = 0 ;
  char timeBuf[25] ;

  switch (sc->col->type)
    {
    case 'i': case 'c':
      li = sc->z ;
      cp = strnew (messprintf ("%ld", li), 0) ; key = _Text ;
      break ;
    case 'f':
      cp = strnew (messprintf ("%lg", sc->z), 0)  ; key = _Text ;
      break ;
    case 'd':
      cp = strnew (messprintf ("\"%s\"", timeShow (sc->u.time, timeBuf, 25)), 0) ; key = _Text ;
      break ;
    case 't':
      cp = strnew (stackText (sc->col->text, sc->u.i), 0) ; key = _Text ;
      break ;
    case 'D':
      cp = strnew (stackText (sc->col->dnaStack, sc->u.k), 0) ; key = _Text ;
      break ;
    case 'P':
      cp = strnew (stackText (sc->col->pepStack, sc->u.k), 0) ; key = _Text ;
      break ;
    default:
      break ;
    }
  switch (sc->col->showType)
    {
    case SHOW_MULTI:
      messfree (cp) ;
      cp = strnew (stackText (sc->col->text, sc->u.i), 0) ; key = _Text ;
      break ;
    default:
      break ;
    }

  s = stackReCreate (s, 256) ;
  /* prepare the text of the query with param substituted */
  n0 = stackMark (s) ;
  pushText (s, col->conditionBuffer) ;
  n1 = stackMark (s) ; pushText (s, " ") ; /* push now, so I can then cat */
  
  stackClear (scStack) ;
  cc = sc ;
  do
    { push (scStack, cc, SC) ;
    } while ((cc = cc->up)) ;
  while (!stackEmpty (scStack))
    { 
      catText (s, " ") ;
      cc = pop (scStack, SC) ; 
      if (!cc->iCol)
	{ catText (s, freeprotect (name (cc->u.k))) ;
	continue ; 
	}
      if (cc->empty) /* use same value as in checkSyntax */
	switch (cc->col->type)
	  {
	  case 'i': case 'f': case 'c':
	    catText (s," 1 ") ;
	    break ;
	  case 'd':
	    catText (s," now ") ;
	    break ;
	  case 'D':
	    catText (s," ATGC ") ;
	    break ;
	  case 'b':  case 'k': case'n': case 'K': case't':
	  default:
	    catText (s," __toto__ ") ;
	    break ;
	  }
      else
	switch (cc->col->type)
	  {
	  case 'i': case 'c':
	    li = cc->z ;
	    catText (s,messprintf ("%ld", li) ) ;
	    break ;
	  case 'f':
	    catText (s,messprintf ("%lg", cc->z) ) ;
	    break ;
	  case 'd':
	    catText (s,messprintf ("\"%s\"", timeShow (cc->u.time, timeBuf, 25))) ;
	    break ;
	  case't':
	    catText (s, freeprotect (stackText (cc->col->text, cc->u.i))) ;
	    break ;
	  case'D':
	    catText (s, cc->col->dnaStack ? stackText (cc->col->dnaStack, cc->u.k) : "0") ;
	    break ;
	  case'P':
	    catText (s, cc->col->pepStack ? stackText (cc->col->pepStack, cc->u.k) : "0") ;
	    break ;
	  case 'b':  case 'k': case'n': case 'K':
	  default:
	    catText (s, freeprotect (name (cc->u.k))) ;
	    break ;
	  }
    }
  
  /* rely on freesubs to actually substitute the param in the query */
  mine = freesettext (stackText (s, n0), stackText (s,n1)) ;
  freespecial ("%\\\"") ;    /* No subshells ($ removed) */
  n2 = stackMark (s) ;
  while (freecard (mine))
    pushText (s, freepos ()) ;

   /* evaluate the query */
  switch (sc->col->type)
    {
    case 'D':
    case 'P':
      n2 = queryRegExpFind (cp, sc->col->dnaCondition, FALSE) ;
      if (n2 > 0)
	found = TRUE ;
      else if (n2 == 0) ;
      else if (n2 < 0)
	spread->interrupt = TRUE ;
      break ;
    default:
      if (condCheckSyntax (stackText (s,n2)))
	{
	  if (queryFind21 (objp,key,stackText (s,n2), cp))
	    found = TRUE ;
	}
      else
	{
	  spread->interrupt = TRUE ;
	}
      break ;
    }
  messfree (cp) ;

  return found ;
} /* scMatch */

/**********************************************************************/

static BOOL scDnaBounds (SC sc, char *buffer, int *dnap, float *fdnap)
{  
  SC cc ;
  int mine, n0, n1 ;
  BOOL ok = TRUE ;
  static   Stack s = 0 ;
  char *cp = buffer ;
  s = stackReCreate (s, 256) ;
  /* prepare the text of the query with param substituted */
  n0 = stackMark (s) ; 

  pushText (s, "[") ;
  catText (s, buffer) ;
  catText (s, "]") ; 

  n1 = 0 ;
  while (*cp && *cp != '%')
    cp++ ;
  if (*cp) /* some param should be substituted */
    {
      n1 = stackMark (s) ;  /* push now, so I can then cat */
      pushText (s, " ") ;
      stackClear (scStack) ;
      cc = sc ;
      do
	{ push (scStack, cc, SC) ;
	} while ((cc = cc->up)) ;
      while (!stackEmpty (scStack))
	{ 
	  catText (s, " ") ;
	  cc = pop (scStack, SC) ; 
	  if (!cc->iCol)
	    { catText (s, freeprotect (name (cc->u.k))) ;
	    continue ; 
	    }
	  if (cc->empty) /* use same value as in checkSyntax */
	    {
	      if (1)
		{
		  catText (s,"____MISSING_PARAMETER____") ;
		}
	      else
		{
		  switch (cc->col->type)
		    {
		    case 'i': case 'f': case 'c':
		      catText (s," 1 ") ;
		      break ;
		    case 'd':
		      catText (s," now ") ;
		      break ;
		    case 'D':
		      catText (s," ATGC ") ;
		      break ;
		    case 'b':  case 'k': case'n': case 'K': case't':
		    default:
		      catText (s," __toto__ ") ;
		      break ;
		    }
		}
	    }
	  else /* we only which to substitute numbers, not text */
	    switch (cc->col->type) /* mhmp 12.07.02  %d -> (%d)  et %f -> (%f)*/
	      {
	      case 'c':
		catText (s,messprintf ("(%d)", cc->u.i) ) ;
		break ;
	      case 'i': 
		catText (s,messprintf ("(%lld)", (long long int)cc->z) ) ;
		break ;
	      case 'f':
		catText (s,messprintf ("(%lf)", cc->z) ) ; /* DO NOT use %g, because 1e+06 is parsed as an addition */
		break ;
	      case 'd':
		catText (s,"0") ; /*messprintf ("\"%s\"", timeShow (cc->u.time, timeBuf, 25)) ; */
		break ;
		case't':
		  catText (s,"0") ; /* freeprotect (stackText (cc->col->text, cc->u.i)) ; */
		break ;
		case'D':
		  catText (s,"0") ; /* freeprotect (stackText (cc->col->dnaStack, cc->u.k)) ; */
		break ;
		case'P':
		  catText (s, "0") ; /*freeprotect (stackText (cc->col->pepStack, cc->u.k)) ; */
		break ;
	      case 'b':  case 'k': case'n': case 'K':
	      default:
		catText (s, freeprotect (name (cc->u.k))) ;
		break ;
	      }
	  
	}
      
      /* rely on freesubs to actually substitute the param in the query */
      /* NB: mieg, jan-2004
	 freesettext calls freenewstream which continuously expands parmStack
	 so after a very long time, acedb crashes
	 we should rather substitute the parameters here in line 
      */
      mine = freesettext (stackText (s, n0), stackText (s, n1)) ;
      freespecial ("%\\\"") ;    /* No subshells ($ removed) */
      
      while (freecard (mine))
	{
	  if (strstr (freepos (), "____MISSING_PARAMETER____"))
	    ok = FALSE ;
	  else
	    ok = queryCalcul (freepos (), dnap, fdnap) ;
	}
    }
  else
    ok = queryCalcul (stackText (s, n0), dnap, fdnap) ;

  return ok ;
} /* scDnaBounds */

/*****************************************************/

static BOOL scQuery (SPREAD spread, SC s2, SC sc)
{
  KEY k = 0, direction = _bsRight ;
  OBJ obj = s2->obj ;
  int nn = 0, type = _Int ;
  double xf = 0, xxf = 0, xxf2 = 0  ;
  long long int xi, xxi = 0, xxi2 = 0 ;
  mytime_t xtm , xxtm = 0 ;
  char *xt ;

  if (!sc->iCol)    /* accept column one */
    return TRUE ;
  if (sc->mark) 
    { 
      if (!obj) messcrash ("Error 1 in scQuery, sorry") ;
      bsGoto (obj, sc->mark) ;
      direction = _bsDown ;
    }
  else
    { 
      switch (sc->col->extend)
	{
	case 'r': /* right of: start inside obj */
  	  if (!obj && sc->scFrom->mark) messcrash ("Error 1 in scQuery, sorry") ;
	  if (!obj && ! (obj = s2->obj = bsCreate (sc->scGrandParent->key)))
	    return FALSE ;
	  bsGoto (obj, sc->scFrom->mark) ;
	  break ;
	case 'f': /* from: start at root */
	  if (obj)
	    bsGoto (obj, 0) ;
	  break ;
	case 'c': /* COPY */
	  sc->u = sc->scFrom->u ;
	  sc->key = sc->scFrom->key ;
	  if (obj)
	    bsGoto (obj, sc->scFrom->mark) ;
	  goto done ;
	  break ;
	case 'p': /* COMPUTE */
	  sc->key = 0 ; 
	  sc->u.f = 0 ;
	  {
	    float xxf ;
	    char *cq =  sc->col->tagTextBuffer ;
	    if (! *cq &&  sc->col->tagStack)
	      cq = stackText (sc->col->tagStack, 0) ;
	    if (scDnaBounds (sc, cq, 0, &xxf))
	      sc->u.f = xxf ;
	    else
	      return FALSE ;
	  }
	  goto done ;
	  break ;
	}
      if (sc->col->type == 'D' && s2->key)
	{
	  Array dna = 0, dnaPiece = 0 ;
	  int i, dna1, dna2, d1, d2 ;
	  BOOL reverse = FALSE ;
	  char *cp, *cq ;

	  if (scDnaBounds (sc, sc->col->dna1Buffer, &dna1, 0) &&
	      scDnaBounds (sc, sc->col->dna2Buffer, &dna2, 0) &&
	      ((dna = s2->col->dnaD) || (dna = dnaGet (s2->key))))
	    {
	      s2->col->dnaD = dna ;
	      if (dna1 < dna2)
		{ d1 = dna1 ; d2 = dna2 ; reverse = FALSE ; }
	      else
		{ d1 = dna2 ; d2 = dna1 ; reverse = TRUE ; }
	      if (reverse)
		{
		  if (!s2->col->dnaR)
		    {
		      s2->col->dnaR = dnaCopy (dna) ;
		      reverseComplement (s2->col->dnaR) ;
		    }
		  dna = s2->col->dnaR ;
		}
	      if (d1 > 0 && d1 >  arrayMax (dna)) /* attaention becuase arrayMax is unsigned int */
		d1 = arrayMax (dna) ;
	      if (d2> 0 && d2 >  arrayMax (dna))
		d2 = arrayMax (dna) ;
	      if (d1 < 1) d1 = 1 ;
	      if (d2 < 1) d2 = 1 ;
	      if (d1 < d2 && d1 > 0)
		{
		  if (reverse) 
		    dna1 = arrayMax (dna) - d2 + 1 ; 
		  else
		    dna1 = d1 ;
		  dnaPiece = arrayCreate (d2 - d1 + 2, char) ;
		  array (dnaPiece, d2 - d1 + 1, char) = 0 ; /* zero terminate */
		  arrayMax (dnaPiece) = d2 - d1 + 1 ; /* restore */
		  for (i = 0, cp = arrp (dnaPiece, 0, char), cq = arrp (dna, dna1 - 1, char) ;
		       i < d2 - d1 + 1 ; cp++, cq++, i++)
		    *cp = dnaDecodeChar [ (int)*cq] ;
		  
		  sc->u.k = stackMark (sc->col->dnaStack) ;
		  pushText (sc->col->dnaStack, arrp (dnaPiece, 0, char)) ;
		  arrayDestroy (dnaPiece) ;
		  return TRUE ;
		}

	      return FALSE ;
	    }
	  return FALSE ; /* mhmp 13.07.02 */
	}
      if (sc->col->type == 'P' && s2->key)
	{
	  Array pep = 0, pepPiece = 0 ;
	  int i, pep1, pep2, d1, d2 ;
	  char *cp, *cq ;

	  if (scDnaBounds (sc, sc->col->dna1Buffer, &pep1, 0) &&
	      scDnaBounds (sc, sc->col->dna2Buffer, &pep2, 0) &&
	      ((pep = s2->col->pep) || (pep = peptideGet (s2->key))))
	    {
	      s2->col->pep = pep ;
	      if (pep1 < pep2)
		{ d1 = pep1 ; d2 = pep2 ; }
	      else
		{ return FALSE ; }
	      if (d2 > arrayMax (pep))
		d2 = arrayMax (pep) ;
	     
	      if (d1 < d2 && d1 > 0)
		{
		  pepPiece = arrayCreate (d2 - d1 + 2, char) ;
		  array (pepPiece, d2 - d1 + 1, char) = 0 ; /* zero terminate */
		  arrayMax (pepPiece) = d2 - d1 + 1 ; /* restore */
		  for (i = 0, cp = arrp (pepPiece, 0, char), cq = arrp (pep, pep1 - 1, char) ;
		       i < d2 - d1 + 1 ; cp++, cq++, i++)
		    *cp = pepDecodeChar [ (int)*cq] ;
		  
		  sc->u.k = stackMark (sc->col->pepStack) ;
		  pushText (sc->col->pepStack, arrp (pepPiece, 0, char)) ;
		  arrayDestroy (pepPiece) ;
		  return TRUE ;
		}

	      return FALSE ;
	    }
		return FALSE ; /* mhmp 13.07.02 */
	}
      if (!obj && sc->col->tag && !bIndexFind (s2->key, sc->col->tag))
	return FALSE ;
      if (!queryFindLocalise3 (sc->col->tagCondition, & (s2->obj), sc->scGrandParent->key))
	return FALSE ;  
      obj = s2->obj ; /* in case we just created obj */
      if (!obj && ! (obj = s2->obj = bsCreate (sc->scGrandParent->key)))
	return FALSE ; /* if queryFind3 did not create the object */
      direction = _bsRight ;
    }
     
  /* do not iterate inside a column on tags or on arithmetic evaluators */
  if (sc->mark && 
      (sc->col->type == 'b' || sc->col->showType != SHOW_ALL))
    return FALSE ;

  switch (sc->col->showType)
    {
    case SHOW_MULTI:
      switch (sc->col->type)
	{
	case 'k': case 'K': case 'n':
	  if (!sc->col->text)
	    { sc->col->text = stackCreate (256) ; pushText (sc->col->text, "toto") ; }
	  direction = _bsDown ; nn = 0 ; 
	  if (bsGetKeyTags (obj, _bsRight, &k))
	    do
	      {
		OBJ obj2 = 0 ;
		
		if (class (k) != _VComment &&
		    class (k) != _VUserSession &&
		    (
		     !sc->col->countCondition ||
		     queryFind3 (sc->col->countCondition, &obj2, k)
		     )) 
		  { 
		    if (nn++)
		      { catText (sc->col->text, "; ") ; catText (sc->col->text, name (k)) ; }
		    else
		      {
			sc->u.i = stackMark (sc->col->text) ;
			pushText (sc->col->text, name (k)) ;
		      }
		  }
		bsDestroy (obj2) ;
	      } while (bsGetKeyTags (obj, direction, &k)) ;
	  if (!nn) return FALSE ;
	  break ;
	default:
	  return FALSE ;	  
	}
      goto done ;
    default:
      break ;
    }

  while (TRUE)
    {
      switch (sc->col->type)
	{
	case 'b':
	  if (sc->mark ||
	      !bsGetKeyTags (obj, _bsHere, &k))
	    return nn ? TRUE : FALSE ; 
	  direction = _bsDown ;
	  while (k &&
		 ((class (k) == _VComment) || (class (k) == _VUserSession))
		 )
	    { k = 0 ; bsGetKeyTags (obj, _bsHere, &k) ; }
          if (!k)   
	    return nn ? TRUE : FALSE ; 
	  sc->u.k = sc->key = k ;
	  break ;
	case 'c': /* count */
	  if (sc->mark ||
	      !bsGetKeyTags (obj, direction, &k))
	    return FALSE ;     
	  direction = _bsDown ; nn = 0 ;
	  do
	    {
	      OBJ obj2 = 0 ;

	      if (class (k) != _VComment &&
		  class (k) != _VUserSession &&
		  (
		   !sc->col->countCondition ||
		   queryFind3 (sc->col->countCondition, &obj2, k)
		   ))nn++ ;
	      bsDestroy (obj2) ;
	    } while (bsGetKeyTags (obj, direction, &k)) ;
	  if (!nn) return FALSE ;
	  sc->u.i = nn ;
	  goto done ;
	  break ;
	case 'k': case'n': case 'K': 
	  if (!bsGetKeyTags (obj, direction, &k))
	    return nn ? TRUE : FALSE ;
	  direction = _bsDown ;
	  while (k &&
		 ((class (k) == _VComment) || (class (k) == _VUserSession))
		 )
	    { k = 0 ; bsGetKeyTags (obj, _bsHere, &k) ; }
          if (!k)   
	    return nn ? TRUE : FALSE ; 
	  sc->u.k = sc->key = k ;
	  break ;
	case 'i': /* LongInt */
	  type = bsType (obj, direction) ;
	  if (type == _LongInt)
	    {
	      if (!bsGetData (obj, direction, type, &xt) ||
		  ! sscanf (xt, "%lli", &(xi))
		  )
		return nn ? TRUE : FALSE ;
	    }
	  else
	    {
	      int x ;
	      if (!bsGetData (obj, direction, _Int, &x))
		return nn ? TRUE : FALSE ;
	      xi = x ;
	    }
	  sc->z = xi ;
	  sc->u.i = xi ;
	  break ;
	case 'f':
	  type = bsType (obj, direction) ;
	  if (type == _LongFloat)
	    {
	      if (!bsGetData (obj, direction, type, &xt) ||
		  ! sscanf (xt, "%lg", &(xf))
		  )
		return nn ? TRUE : FALSE ;
	    }
	  else
	    {
	      float x ;
	      if (!bsGetData (obj, direction, _Float, &x))
		return nn ? TRUE : FALSE ;
	      xf = x ;
	    }
	  sc->z = xf ;
	  sc->u.f = xf ;
	  break ;
	case 'd':
	  if (!bsGetData (obj, direction, _DateType, &xtm))
	    return nn ? TRUE : FALSE ;
	  sc->u.time = xtm ;
	  break ;
	case 't':
	  {
	    /* this is a hack, but it seems that a Text is preceded by a cell
	       of type Text but with no text inside it, at least sometimes */
	  encore:
	    if (!bsGetData (obj, direction, _Text, &xt))
	      return nn ? TRUE : FALSE ;
	    if (!xt || !*xt)
	      { direction = _bsDown ; goto encore ; }
	  }
	  sc->u.i = stackMark (sc->col->text) ;
	  pushText (sc->col->text, xt) ;
	  break ;
	}
      /* SHOW_ALL=0, SHOW_MIN, SHOW_MAX, SHOW_AVG */
      switch (sc->col->showType)
	{
	case SHOW_ALL:
	  goto done ;
	case SHOW_MIN:
	  switch (sc->col->type)
	    {
	    case 'i':
	      if (!nn || xxi > xi) xxi = xi ;
	      sc->u.i = xxi ; /* reset */
	      sc->z = xxi ; /* reset */
	      break ;
	    case 'f':
	      if (!nn || xxf > xf) xxf = xf ;
	      sc->u.f = sc->z = xxf ;
	      break ;
	    case 'd':
	      if (!nn || xxtm > xtm) xxtm = xtm ;
	      sc->u.time = xxtm ;
	      break ;
	    default: /* undefined on other types */
	      return FALSE ;
	    }
	  break ;
	case SHOW_MAX:
	  switch (sc->col->type)
	    {
	    case 'i':
	      if (!nn || xxi < xi) xxi = xi ;
	      sc->z =  xxi ;
	      sc->u.i = xxi ; /* reset */
	      break ;
	    case 'f':
	      if (!nn || xxf < xf) xxf = xf ;
	      sc->u.f = sc->z = xxf ;
	      break ;
	    case 'd':
	      if (!nn || xxtm < xtm) xxtm = xtm ;
	      sc->u.time = xxtm ;
	      break ;
	    default: /* undefined on other types */
	      return FALSE ;
	    }
	  break ;
	case SHOW_AVG:
	  switch (sc->col->type)
	    {
	    case 'i':
	      if (!nn) xxi = 0 ;
	      xxi += xi ;
	      sc->z = xxi/ (nn + 1) ; /* reset */
	      sc->u.i = xxi/ (nn + 1) ; /* reset */
	      break ;
	    case 'f':
	      if (!nn) xxf = 0 ;
	      xxf += xf ;
	      sc->u.f = sc->z = xxf/ (nn + 1) ;
	      break ;
	    default: /* undefined on other types */
	      return FALSE ;
	    }
	  break ;
	case SHOW_VAR:
	  switch (sc->col->type)
	    {
	    case 'i':
	      if (!nn) { xxi = 0 ; xxi2 = 0 ; sc->u.i = 0 ; }
	      xxi += xi ; xxi2 += xi * xi ;
	      if (nn > 1)
		sc->u.i = sc->z = (xxi2  - (xxi*xxi/ ((float)nn)))/ (nn - 1.0); /* reset */
	      break ;
	    case 'f':
	      if (!nn) { xxf = 0 ; xxf2 = 0 ; sc->u.f = 0 ; }
	      xxf += xf ; xxf2 += xf * xf ;
	      if (nn > 1)
		sc->u.f = sc->z = (xxf2 - xxf*xxf/nn)/ (nn - 1.0) ;
	      break ;
	    default: /* undefined on other types */
	      return FALSE ;
	    }
	  break ;
	case SHOW_SUM:
	  switch (sc->col->type)
	    {
	    case 'i':
	      if (!nn) xxi = 0 ;
	      xxi += xi ;
	      sc->z = xxi ; /* reset */
	      sc->u.i = xxi ; /* reset */
	      break ;
	    case 'f':
	      if (!nn) xxf = 0 ;
	      xxf += xf ;
	      sc->u.f = sc->z = xxf ;
	      break ;
	    default: /* undefined on other types */
	      return FALSE ;
	    }
	  break ;
	default:
          break ;
	}
      nn++ ;
      direction = _bsDown ;
    }
done:
  if (obj)
    {
      sc->mark = bsMark (obj, sc->mark) ;
      sc->parent = bsParentKey (obj) ;
    }

  return TRUE ;
} /* scQuery */

/*****************************************************/

static BOOL scCondition (SPREAD spread, SC sc)
{
  BOOL found = FALSE;
  BOOL isCond = *sc->col->conditionBuffer ? TRUE : FALSE;
  OBJ obj = 0 ;
  KEY key ;
  char *cp = 0 ;
  char timeBuf[25] ;

  if (!isCond) return TRUE ;
  switch (sc->col->type)
    {
    case 'b': 
      return TRUE ;
      break ;
    case 'n': case 'k': 
      if (!lexIsInClass (sc->key, sc->col->classe))
	return FALSE ;
      key = sc->key ;
      break ;
    case 'K':  /* no class check in case K */
      key = sc->key ;
      break ;
    case 'c':
       cp = strnew (messprintf ("%d", sc->u.i),0 ) ; key = _Text ;
      break ;
    case 'i': 
      cp = strnew (messprintf ("%lld", (long long int)sc->z),0 ) ; key = _Text ;
      break ;
    case 'f':
       cp = strnew (messprintf ("%lf", sc->z),0 ) ; key = _Text ;
      break ;
    case 'd':
      cp = strnew (messprintf ("\"%s\"", timeShow (sc->u.time, timeBuf, 25)),0) ; key = _Text ;
      break ;
    case 't':
      cp = strnew (stackText (sc->col->text, sc->u.i),0) ; key = _Text ;
      break ;
    case 'D': 
      cp = strnew (stackText (sc->col->dnaStack, sc->u.k),0) ; key = _Text ;
      break ;
    case 'P': 
      cp = strnew (stackText (sc->col->pepStack, sc->u.k),0) ; key = _Text ;
      break ;
    default:
      key = _Text ;
      break ;
    }

  switch (sc->col->showType)
    {
    case SHOW_MULTI:
      messfree (cp) ;
      cp = strnew (stackText (sc->col->text, sc->u.i), 0) ; key = _Text ;
      break ;
    default:
      break ;
    }


  if (sc->col->conditionCondition)
    found = sc->col->conditionCondition && 
      queryFind31 (sc->col->conditionCondition, &obj, key, cp) ;
  else if (sc->col->dnaCondition)
    found = queryRegExpFind (cp, sc->col->dnaCondition, FALSE) ;
  else
    found = scMatch (spread, sc, &obj, key) ;

  messfree (cp) ;
  bsDestroy (obj) ;

  return found ;
} /* scCondition */

/*****************************************************/

static int scLoopDown (SPREAD spread, SC s2, SC sc)
{ 
  int n = 0 , mand = sc->col->mandatory ;
  
  if (
      ! (sc->scFrom && sc->scFrom->empty) &&
      ! (sc->scFrom && sc->scFrom->col && sc->scFrom->col->nonLocal))
    while (TRUE) /* loop on the whole column */
    {
      if (messIsInterruptCalled ())
	{ 
	  spread->interrupt = TRUE ;
	  return n ;
	}

      if (!scQuery (spread, s2, sc))
	break ;  /* no more target */
      if (mand == 0) /* NULL required */
	return 1 ; /* do not export */
      if (scCondition (spread, sc))  /* target acceptable */
	n += scLoopRight (spread, sc) ; /* expand and export */
      if (sc->col->showType != SHOW_ALL) 
	break ;
    }
  if (!n && mand != 2) /* null or optional */
   { sc->empty = TRUE ;
     n = scLoopRight (spread, sc) ; /* explore once */
   }
  return n ;  
} /* scLoopDown */
  
/*****************************************************/

static void scExport (SPREAD spread, SC sc)
{ 
  SC cc ;
  int i = 0, nCol = spread->nCol ;
  SPCELL *sp ;
  COL *colonne ;
  Array a = 0 ;

  stackClear (scStack) ;

  
  for (i = sc->iCol + 1, colonne = arrp (spread->colonnes, sc->iCol, COL) + 1 ;
       i < nCol ; i++, colonne++) 
    if (colonne->mandatory == 2) /* mandatory required */
      return ; /* missing non optional column, do not export */

  a = arrayCreate (nCol, SPCELL) ;cc = sc ;
  do
    { push (scStack, cc, SC) ;
    } while ((cc = cc->up)) ;
  for (i = 0 ; i < nCol ; i++)
    { 
      sp = arrayp (a, i, SPCELL) ;
      if (!stackEmpty (scStack))
	{ cc = pop (scStack, SC) ;
	sp->empty = cc->empty ;
	sp->u = cc->u ;sp->z = cc->z ;
	sp->parent = cc->parent ;
	sp->grandParent = cc->scFrom ?  cc->scFrom->key : cc->key ;
	}	
      else
	sp->empty = TRUE ;
    }
  array (spread->tableau, arrayMax (spread->tableau), Array) = a ;

  return;
} /* scExport */

/*****************************************************/

static int scLoopRight (SPREAD spread, SC sc)
{ 
  int n = 0, from ;
  SC sr = 0, s2, srFrom = 0 ;
  
  if (!spread->interrupt && (!sc->iCol || (sc->col &&  sc->col->type)) &&
     sc->iCol + 1 < spread->nCol)
    { 
      sr = scAlloc () ;
      sr->up = sc ; 
      sr->iCol = sc->iCol + 1 ;
      sr->col = arrp (spread->colonnes, sr->iCol, COL) ;
      srFrom = sr ; from = sr->col->from - 1 ;
      while (srFrom->iCol != from && srFrom->up) srFrom = srFrom->up ;
      if (srFrom == sr || srFrom->iCol != from)
	{ messerror ("Confusion 1 in scLoopRight") ; return 1 ; }
      sr->scFrom = srFrom ;
      s2 = sr ;
      while (s2 && s2->col && s2->scFrom && s2->col->extend != 'f') 
	s2 = s2->scFrom ;  
			if (s2 && s2->col && s2->scFrom) /* mhmp 08.07.02  + && s2->scFrom */
				s2 = s2->scFrom ;
      sr->scGrandParent = s2 ;

      if (sr->col && sr->col->type)
	n = scLoopDown (spread, s2, sr) ;

      scFree (sr) ;
    }
  if (!n) /* then export self */
    { n = 1 ; scExport (spread, sc) ;}

  if (sc->obj)
    bsDestroy (sc->obj) ; /* release memory NOW !!!!
			     this is the main purpose this SC method */
  arrayDestroy (sc->col->pep) ;
  arrayDestroy (sc->col->dnaD) ;
  arrayDestroy (sc->col->dnaR) ;
  return n ;
} /* scLoopRight */

/*****************************************************/
/************* Public functions **********************/
/*****************************************************/

KEY spreadRecomputeKeySet (SPREAD spread, KEYSET ks, int last, int minx, int maxx) 
     /* public function */
{
  int i, new = 0, mx = keySetMax (ks) ;
  KEY key ;
  BOOL stop = last && spread->sortColonne == 1 && !spread->precompute ? TRUE : FALSE ; 
  int nn = 200 ;
  SC sc ;

  if (!spreadCheckConditions (spread))
    {
      spreadDestroyConditions (spread) ;
      return 0 ;
    }

  if (externalServer && spread->isActiveKeySet)
    {   /* export active keyset and compute on server side */
      Stack ss = 0 ;
      int n, level ;
      char *cp ;
      TABLE *tt = 0 ;

       /* export active keyset was done by the filter call, compute on server side */
      ss = stackCreate (100) ;
      level = freeOutSetStack (ss) ; 
      spreadDoExportDefinitions (spread) ;
      freeOutClose (level) ;
      n = stackMark (ss) ;
      if (!n) { stackDestroy (ss) ; return 0 ; }

      pushText (ss, "Table -active -a = ") ;
      cp =stackText (ss, 0) ;
      catText (ss, cp) ;
      cp = stackText (ss, n) ;
      tt = externalTableMaker (cp) ;
      stackDestroy (ss) ;
      spreadTable2Tableau (spread,tt) ;
      tableDestroy (tt) ;
      spreadDestroyConditions (spread) ;
      return 0 ;
    }

  if (last == 1 && ! (minx == 1)) last = 0 ;
  scStack = stackReCreate (scStack, 16) ;

  spread->interrupt = FALSE ;
  spreadCleanUpTableau (spread) ; 
  if (maxx <= 0)
    maxx = mx ;

  if (last || spreadCheckConditions (spread))
    for (i = last ; i < maxx ; i++)
      { 
	new = i + 1 ; /* for encore */
	if (minx > 0 && maxx > 0 &&
	    (i < minx - 1 || i > minx + maxx - 2))
	  continue ;
	key = keySet (ks,i) ;
/*	
	if (i + 2 < mx)
	  bsAnticipate (keySet (ks, i + 2)) ;
*/

	sc = scAlloc () ;
	sc->key = key ; sc->u.k = key ;
	sc->parent = key ; 
	sc->iCol = 0 ; sc->up = 0 ;
	sc->col = arrayp (spread->colonnes, sc->iCol, COL) ;
	scLoopRight (spread, sc) ; 
	scFree (sc) ;
	if (spread->interrupt)
	  break ;
	if (stop && --nn <= 0) /* will not stop if starts at 0 */
	  break ;
	new = 0 ; /* remains 0 at end of loop */
      }
  spreadReorder (spread) ;
  spread->modified = FALSE ;
  if (last >= maxx) new = 0 ;  /* destroy */
  if (!stop || !new) spreadDestroyConditions (spread) ;

  return new ;
} /* spreadRecomputeKeySet */

/*************/

  /* return the alphabetically ordered 1st column */
KEYSET spreadGetKeyset (SPREAD spread) 
     /* public function */
{
  KEYSET ksTmp ; 
  KEYSET ks ;
  COL *c ;

  if (!spread->modified)
    return 0 ;
  
  if (!arrayMax (spread->colonnes))
    { messout ("// First define the table") ;
      return 0 ;
    }
  
  c = arrp (spread->colonnes, 0, COL) ;
  if (ace_lower (c->type) != 'k')
    { messout ("Column 1 should be type Object") ;
      return 0 ;
    }
  if (!c->classe)
    { messout ("First choose the class of column 1 ") ;
      return 0 ;
    }
  
  if (*c->conditionBuffer &&
      ! condCheckSyntax (messprintf (" %s", c->conditionBuffer)))
    { 
      messout ("First fix the syntax error in column 1's condition") ;
      return 0 ;
    }

  ksTmp = query ( 0, messprintf (">?%s %s", 
			     name (c->classe), 
			     c->conditionBuffer) ) ;

  ks = keySetAlphaHeap (ksTmp, keySetMax (ksTmp)) ;
  keySetDestroy (ksTmp) ;

  return ks ;
} /* spreadGetKeyset */

/*************/

KEYSET spreadFilterKeySet (SPREAD spread, KEYSET ks)
     /* public function */
{
  KEYSET ks1 = 0, ks2 = keySetCreate () ;
  COL *c ;
 
  if (!spread || !arrayMax (spread->colonnes))
    return ks2 ;

  c = arrp (spread->colonnes, 0, COL) ;
  if (!c->classe)
    return ks2 ;
  
  if (*c->conditionBuffer &&
      ! condCheckSyntax (messprintf (" %s", c->conditionBuffer)))
    return ks2 ;

  ks1 = query (ks, messprintf ("CLASS %s", name (c->classe))) ;

  if (*c->conditionBuffer)
    ks2 = query (ks1, c->conditionBuffer) ;
  else
    { ks2 = ks1 ; ks1 = 0 ; }
  keySetDestroy (ks1) ;

  return ks2 ;
} /* spreadFilterKeySet */


/*************/

static void spreadTable2Tableau (SPREAD spread, TABLE *tt)
{
  int nn, j, j1, tMax = tt ? tableMax (tt) : 0 ;
  COL *c ;  
  int  maxCol = arrayMax (spread->colonnes) ;
  Array a = 0 ;
  SPCELL *su ;
  BSunit *u ;
  char *cp ;


  colonnes = spread->colonnes ;  /* needed by spreadOrder later on */
  pos2col = spread->pos2col ;
  sortColonne = spread->sortColonne - 1 ;
  spreadCleanUpTableau (spread) ;

  if (!tt || !tMax || !maxCol)
    return ;

  for (nn = 0 ; nn < tMax ; nn++)
    {
      a = array (spread->tableau, nn, Array) = arrayCreate (maxCol, SPCELL) ;
      array (a,maxCol - 1,SPCELL).empty = FALSE ; /* make room */
      for (j = 0 ; j < maxCol; j++)
	{
	  j1 = arr (pos2col,j, int)  ;
	  c = arrp (spread->colonnes,j1,COL) ;
	  if (c->hidden || tabEmpty (tt,nn,j1))
	    {
	      array (a,j1,SPCELL).empty = TRUE ;
	      continue ;
	    }
	  su = arrayp (a, j1,SPCELL) ;
	  u = &(su->u) ;
	  if (c->showType == SHOW_MULTI) 
	    {
	      cp =  tabString (tt,nn,j1) ;
	      if (cp && *cp)
		{
		  if (!stackExists (c->text))
		    { c->text = stackCreate (80) ; pushText (c->text, "toto") ; }
		  u->i = stackMark (c->text) ;
		  pushText (c->text, cp) ;
		}
	      else
		array (a,j1,SPCELL).empty = TRUE ;
	    }
	  else switch (c->type)
	    {
	    case 0:
	      break ;
	    case 'c':
	      u->i = tabInt (tt,nn,j1) ;
	      break ;
	    case 'i': 
	      su->z = u->i = tabInt (tt,nn,j1) ;
	      break ;
	    case 'f': 
	      su->z = u->f = tabFloat (tt,nn,j1) ;
	      break ;
	    case 't':
	      cp =  tabString (tt,nn,j1) ;
	      if (cp && *cp)
		{
		  if (!stackExists (c->text))
		    { c->text = stackCreate (80) ; pushText (c->text, "toto") ; }
		  u->i = stackMark (c->text) ;
		  pushText (c->text, cp) ;
		}
	      else
		array (a,j1,SPCELL).empty = TRUE ;
	      break ;
	    case 'd': 
	      u->time = tabDate (tt,nn,j1) ;
	      break ;
	    case 'k': case 'K': case 'n': case 'b':
	      u->k = array (a, j1,SPCELL).parent = tabKey (tt,nn,j1) ;  /* parent has to be guessed */
	      break ;
	    }
	}
    }
}

/*************/

BOOL spreadDoRecompute (SPREAD spread)
     /* private within SpreadPackage */
{ 
  KEYSET ks;
  int level, n ;
  char *cp ;
  TABLE *tt = 0 ;

  if (!spreadCheckConditions (spread))
    return FALSE ;

  if (externalServer) /* compute on server side */
    {
      Stack ss = stackCreate (100) ;
      level = freeOutSetStack (ss) ; 
      spreadDoExportDefinitions (spread) ;
      freeOutClose (level) ;
      n = stackMark (ss) ;
      if (!n) { stackDestroy (ss) ; return FALSE ; }
      pushText (ss, "Table -a = ") ;
      cp =stackText (ss, 0) ;
      catText (ss, cp) ;
      cp = stackText (ss, n) ;
      tt = externalTableMaker (cp) ;
      stackDestroy (ss) ;
      spreadTable2Tableau (spread,tt) ;
      tableDestroy (tt) ;
      return TRUE ;
    }
  ks = spreadGetKeyset (spread) ;
  if (!ks)
    return FALSE ;
  spreadRecomputeKeySet (spread, ks, 0, 0, 0) ;
  keySetDestroy (ks) ;

  return TRUE ;
} /* spreadDoRecompute */

/******************************************************/
/************ Dumper **********************************/
/******************************************************/

BOOL spreadDumpLegend (SPREAD spread, char style)
{
  char *cp, *cq, *cr ;
  int j, iCol, n1 = 0 ;
  COL *c ;
  int  maxCol = arrayMax (spread->colonnes) ;

  if (style == 'x' && spread->showTitle)
    {
      for (j = iCol = 0 ; j < maxCol; j++)
	{ 
	  c = arrp (spread->colonnes, arr (spread->pos2col,j, int) ,COL) ;
	  if (!c->hidden)
	    { 
	      cp = c->subtitleBuffer ;
	      cq = c->legendBuffer ; 
	      if (cp && cq && n1 < strlen (cp))
		n1 = strlen (cp) ;
	    }
	}
      if (n1)
	{
	  cr = "<b>Legend : </b>" ;
	  if (* (spread->titleBuffer))
	    freeOutf ("<h3>%s</h3>\n", spread->titleBuffer) ;

	  freeOutf ("<table border=0>\n") ;
	  for (j = iCol = 0 ; j < maxCol; j++)
	    { 
	      c = arrp (spread->colonnes, arr (spread->pos2col,j, int) ,COL) ;
	      if (!c->hidden)
		{ 
		  cp = c->subtitleBuffer ;
		  cq = c->legendBuffer ;
		  if (cp && cq && *cp && *cq) 
		    {
		      freeOutf ("<tr><td>%s</td><td><b>%s</b></td><td>%s</td></tr>\n", cr, cp, cq) ;
		      cr = "" ;
		    }
		}
	    }
	  freeOut ("</table>\n") ;
	}
    }
  return n1 ? TRUE : FALSE ;
}

  /* case noHide, always put something so that %n counts to correct column */
static int lastPrintedCol = -1 ; /* lastPrintedCol col ever printed */
static void spreadAscii (SPREAD spread, 
			 int sep0, int maxCol, 
			 Array lineArray, int iLine, BOOL noHide, 
			 char style, BOOL title, int nLine, char *href)
{ 
  int j, j1, x=0, colMax = 0 , iCol ;
  COL *c ;
  SPCELL *su ;
  BSunit *u ;
  char sep[2], *cp = 0, oldType ;
  BOOL empty, NEW_JAVA = FALSE ;
  static KEY lastUk, xColor ;
  char bufvoid[1] = {'v'} ;
  char bufinteger[1] = {'i'} ;
  char buffloat[1] = {'f'} ;
  char bufdate[1] = {'d'} ;
  char buftext[1] = {'t'} ;
  char buf0[2] = {0, '\n'} ;
  char bufk[2] = {'k', 0} ;
  char bufK[2] = {'K', 0} ;
  char buftag[1] = {'g'} ;
  char beginline[1] ={'.'};
  char timeBuf[25] ;

  /*   if (getenv ("NEWJAVA") ) NEW_JAVA = TRUE ;  */
  pos2col = spread->pos2col ;

  sep[0] = sep0 ? sep0 : '\t' ; sep [1] = 0 ;
  if (!style) style = 'a' ;

  if (style != 'C' && title)
    {
      if (style == 'x') 
	freeOut ("<tr bgcolor=#afafff>\n") ;
      else
	freeOut ("# ") ;
      for (j = iCol = x = 0 ; j < maxCol; j++)
	{
	  j1 = arr (pos2col,j, int) ;
	  c = arrp (spread->colonnes,j1,COL) ;
	  if (noHide || !c->hidden)
	    { 
	      if (style == 'x') 
		freeOut ("<td>") ;
	      else if (iCol++) 
		{
		  if (*sep == ' ')
		    { 
		      freeOut (" ") ; x++ ;
		      for (; x < c->width ; x++)
			freeOut ( " ") ;
		    }
		  else
		    freeOut ( sep) ;
		}	
	      cp = c->subtitleBuffer ; x = 0 ;
	      if (cp && *cp)
		{ freeOut (cp) ; x = strlen (cp) ; }
	      else if (style == 'x')
		freeOut ("&nbsp;") ;
	      if (style == 'x') 
		freeOut ("</td>\n") ;
	    }
	}
      if (style == 'x') 
	freeOut ("</tr>") ;
      freeOut ("\n") ;
    }
   	
  if (iLine == 0)
    {
      if (style == 'C')
	{
	  beginline[0] = '>' ;
	}
      if (style == 'x')
	{
	  xColor = 1 ;
	}
    }

  if (style == 'x') 
    {
      KEY newUk = 0 ;
      
      u = 0 ;
      for (j = iCol = 0 ; j < maxCol; j++)
	{
	  j1 = arr (pos2col,j, int)  ;
	  c = arrp (spread->colonnes, j1 ,COL) ;
	  if (noHide || !c->hidden)
	    {
	      newUk = (arr (lineArray, j1,SPCELL).u).k ;
	      break ;
	    }
	}
      xColor = (xColor + (!newUk || (newUk != lastUk) ? 1 : 0)) % 2 ;
      if (xColor)
	freeOut ("<tr bgcolor=#efefff>") ;
      else
	freeOut ("<tr>") ;
	
      lastUk = newUk ? newUk : 1 ;
    }

  for (j = iCol = 0 ; j < maxCol; j++)
    {
      j1 = arr (pos2col,j, int) ;
      c = arrp (spread->colonnes, j1 ,COL) ;
      oldType = c->type ;
      if (c->showType == SHOW_MULTI)
	c->type = 't' ;
      if (noHide || !c->hidden)
	{ 
	  iCol++ ;
	  if (style == 'x')
	    freeOut ("<td>") ;
	  else if (style != 'C' && iCol > 1) 
	    {
	      if (*sep == ' ')
		{ 
		  freeOut ( " ") ; x++ ;
		  for (; x < c->width ; x++)
		    freeOut ( " ") ;
		    }
	      else
		freeOut ( sep) ;
	    }	

	  su = arrp (lineArray, j1,SPCELL) ;
	  u = & (su->u) ;
	  empty = arr (lineArray, j1,SPCELL).empty ;
	  if (!NEW_JAVA) lastPrintedCol = -1 ;
	  x = 0 ;
   	  if (style == 'C')
	    {
	      freeOutBinary (beginline, 1);
	      beginline[0]='>';
	      colMax++ ;
	    }

          if (empty)
	    switch (style)
	      {
	      case 'a':
		freeOut ("NULL") ; x = 4 ;
		break ;
	      case 'C':
		freeOutBinary (bufvoid, 1) ;
		break ;
	      case 'x':
		freeOut ("&nbsp;") ; x = 4 ;
	      default:
		break ;
	      }
	  else  switch (c->type)
	    {
	    case 0:
	      break ;
	    case 'k': case 'K': case 'b': case 'n': /* key nextKey tag nextTag */
	      if (iCol == spread->exportedColumn && 
		  spread->exportedKeySet)
		keySet (spread->exportedKeySet, keySetMax (spread->exportedKeySet)) = u->k ;
	      switch (style)
		{
		case 'j': 
		  if (! NEW_JAVA || j > lastPrintedCol)
		    { cp = messprintf ("?%s?%s?", className (u->k), 
				      freejavaprotect (name (u->k))) ;
		    lastPrintedCol = j ;
		    }
		  else
		    cp = messprintf ("#%s?", 
				    freejavaprotect (name (u->k))) ;
		  freeOut (cp) ;
		  x = strlen (cp) ;
		  break ;
		case 'a': 
		  if (arrayExists (spread->oldClass)
		      && j1 < arrayMax (spread->oldClass)
		      && class (u->k) != array (spread->oldClass,j1,int))
		    { array (spread->oldClass,j1,int) = class (u->k) ;
		    cp = messprintf ("%s:%s", class (u->k) ? className (u->k) : "Tag",
				    freeprotect (name (u->k))) ;
		    }
		  else
		    cp = freeprotect (name (u->k)) ;
		  freeOut (cp) ;
		  x = strlen (cp) ;
		  break ;
		case 'C':
		  if (! class (u->k)) /* tag */
		    {
		      freeOutBinary (buftag,1) ;
		      freeOutBinary ((char*)& (u->k), 4) ;
		    }
		  else
		    {
		      if (iskey (u->k) == 2)
			{ bufK[1] = class (u->k) ; freeOutBinary (bufK,2) ; }
		      else
			{ bufk[1] = class (u->k) ; freeOutBinary (bufk,2) ; }
		      cp = name (u->k) ;
		      freeOutBinary (cp, strlen (cp)) ;
		      freeOutBinary (buf0, 2) ;
		    }
		  break ;
		case 'x': 
		  cp = name (u->k) ; 
		  if (href &&
		      (
		       !strcasecmp (className (u->k), "Sequence") ||
		       !strcasecmp (className (u->k), "Element") ||
		       !strcasecmp (className (u->k), "cDNA_clone") ||
		       !strcasecmp (className (u->k), "Gene") ||
		       !strcasecmp (className (u->k), "Transcribed_gene") ||
		       !strcasecmp (className (u->k), "mRNA") ||
		       !strcasecmp (className (u->k), "Procuct") ||
		       !strcasecmp (className (u->k), "GeneId") ||
		       !strcasecmp (className (u->k), "LocusId") 
		       ))
		    freeOutf ("%sc=%s&l=%s\">%s</a>", href, className (u->k), cp, cp) ;
		  else
		    freeOut (cp) ;
		  x = strlen (cp) ;
		  break ;
		default: 
		  cp = freeprotect (name (u->k)) ; 
		  freeOut (cp) ;
		  x = strlen (cp) ;
		  break ;
		}
	      break ;
	    case 'i': case 'c': /* integer count */
	      switch (style)
		{
		case 'j': 
		  if (! NEW_JAVA || j > lastPrintedCol)
		    { 
		      cp = messprintf ("?int?%lli?", (long long int) su->z) ; 
		      lastPrintedCol = j ;
		    }
		  else
		    cp = messprintf ("#%lld?", (long long int) su->z) ; 
		  freeOut (cp) ;
		  x = strlen (cp) ;
		  break ;
		case 'C': 
		   freeOutBinary (bufinteger,1) ;
		   freeOutBinary ((char*)& (u->i), 4) ;
		  break ;
		default:
		  cp = messprintf ("%d", u->i) ; 
		  freeOut (cp) ;
		  x = strlen (cp) ;
		  break ;
		}
	      break ;
	    case 'f':
	      switch (style)
		{
		case 'j': 
		  if (! NEW_JAVA || j > lastPrintedCol)
		    { cp = messprintf ("?float?%lg?", su->z) ; 
		    lastPrintedCol = j ;
		    }
		  else
		    cp = messprintf ("#%lg?", su->z) ; 
		  freeOut (cp) ;
		  x = strlen (cp) ;
		  break ;
		case 'C': 
		  freeOutBinary (buffloat,1) ;
		  freeOutBinary ((char*)& (u->f), 4) ;
		  break ;
		default: 
		  {
		    float zf = u->f > 0 ? u->f : -u->f ;
		    int izf = u->f + .1 ;
		    if (izf > 0 && zf - izf < ACE_FLT_RESOLUTION)
		      cp = messprintf ("%d", u->f > 0 ? izf : -izf) ;
		    else
		      cp = messprintf ("%g", u->f) ;
		  }
		  freeOut (cp) ;
		  x = strlen (cp) ;
		  break ;
		}
	      break ;
	    case 'd':
	      switch (style)
		{
		case 'j': 
		  if (! NEW_JAVA || j > lastPrintedCol)
		    { cp = messprintf ("?date?%s?", timeShowJava (u->time, timeBuf, 25)) ; 
		    lastPrintedCol = j ;
		    }
		  else
		    cp = messprintf ("#%s?", timeShowJava (u->time, timeBuf, 25)) ; 
		  freeOut (cp) ;
		  x = strlen (cp) ;
		  break ;
		case 'C': 
		  freeOutBinary (bufdate,1) ;
		  freeOutBinary ((char*)& (u->time), 4) ;
		  break ;
		case 'x':
		  cp = timeShow (u->time, timeBuf, 25) ;
		  freeOut (cp) ;
		  x = strlen (cp) ;
		  break ;
		default: 
		  cp = messprintf ("\"%s\"", timeShow (u->time, timeBuf, 25)) ;
		  freeOut (cp) ;
		  x = strlen (cp) ;
		  break ;
		}
	      break ;
	    case 't':
	      if (stackExists (c->text) && u->i)
		{
		  cp =stackText (c->text, u->i) ;
		  switch (style)
		    {
		    case 'j': 
		      if (! NEW_JAVA || j > lastPrintedCol)
			{ cp = messprintf ("?txt?%s?", freejavaprotect (cp)) ; 
			lastPrintedCol = j ;
			}
		      else
			cp = messprintf ("#%s?", freejavaprotect (cp)) ; 
		      freeOut (cp) ;
		      x = strlen (cp) ;
		      break ;
		    case 'C': 
		      cp = freeprotect (cp) ;
		      freeOutBinary (buftext,1) ;
		      freeOutBinary (cp, strlen (cp)) ;
		      freeOutBinary (buf0, 2) ;
		      break ;
		    case 'x':
		      freeOut (cp) ;
		      x = strlen (cp) ;
		      break ;
		    default: 
		      cp = freeprotect (cp) ; 
		      freeOut (cp) ;
		      x = strlen (cp) ;
		      break ;
		    }
		}
	    case 'D':
	      if (stackExists (c->dnaStack) && u->k)
		{
		  cp =stackText (c->dnaStack, u->k) ;
		  switch (style)
		    {
		    case 'j': 
		      if (! NEW_JAVA || j > lastPrintedCol)
			{ cp = messprintf ("?txt?%s?", freejavaprotect (cp)) ; 
			lastPrintedCol = j ;
			}
		      else
			cp = messprintf ("#%s?", freejavaprotect (cp)) ; 
		      freeOut (cp) ;
		      x = strlen (cp) ;
		      break ;
		    case 'C': 
		      cp = freeprotect (cp) ;
		      freeOutBinary (buftext,1) ;
		      freeOutBinary (cp, strlen (cp)) ;
		      freeOutBinary (buf0, 2) ;
		      break ;
		    case 'x':
		      freeOut (cp) ;
		      x = strlen (cp) ;
		      break ;
		    default: 
		      freeOut (cp) ;
		      x = strlen (cp) ;
		      break ;
		    }
		}
	      else
		{
		  if (noHide && style != 'C')
		    { 
		      freeOut ("000") ;
		      x = 3 ;
		    }
		  if (lastPrintedCol > j) lastPrintedCol = j ;
		}   
	      break ;
	    case 'P':
	      if (stackExists (c->pepStack) && u->k)
		{
		  cp =stackText (c->pepStack, u->k) ;
		  switch (style)
		    {
		    case 'j': 
		      if (! NEW_JAVA || j > lastPrintedCol)
			{ cp = messprintf ("?txt?%s?", freejavaprotect (cp)) ; 
			lastPrintedCol = j ;
			}
		      else
			cp = messprintf ("#%s?", freejavaprotect (cp)) ; 
		      freeOut (cp) ;
		      x = strlen (cp) ;
		      break ;
		    case 'C': 
		      cp = freeprotect (cp) ;
		      freeOutBinary (buftext,1) ;
		      freeOutBinary (cp, strlen (cp)) ;
		      freeOutBinary (buf0, 2) ;
		      break ;
		    case 'x':
		      freeOut (cp) ;
		      x = strlen (cp) ;
		      break ;
		    default: 
		      freeOut (cp) ;
		      x = strlen (cp) ;
		      break ;
		    }
		}
	      else
		{
		  if (noHide && style != 'C')
		    { 
		      freeOut ("000") ;
		      x = 3 ;
		    }
		  if (lastPrintedCol > j) lastPrintedCol = j ;
		}   
	      break ;
	    default:
	      messcrash ("Unknown type in spreadDump") ;
	    }
	  if (style == 'x')
	    freeOut ("</td>\n") ;
	}
      c->type = oldType ;
    }
  
  if (style == 'C')
    { 
      char bfret[3];
      
      bfret[0] = 'l';
      colMax-- ;
      while (colMax > 255)
	{
	  bfret[1]= (char)255; freeOutBinary (bfret, 2) ;
	  colMax -= 255;
	}
      
      bfret[1] = colMax ; bfret[2] = '\n' ;
      freeOutBinary (bfret, 3) ;
    }
  else if (style == 'x') 
    freeOut ("</tr>") ;
  else
    freeOut ("\n") ;
    
  return;
} /* spreadAscii */

/*****************/

static void spreadDoDumpFormat (SPREAD spread)
{
  int i, j, j1 ;
  int  maxCol = arrayMax (spread->colonnes) ;
  COL *c ;
  BSunit *u ;
  Array lineArray = 0 ;
  KEY key ;

  pos2col = spread->pos2col ;
  spread->oldClass = arrayReCreate (spread->oldClass,maxCol, int) ;

  if (maxCol)
    freeOut ("Format ") ;
  for (j = 0 ; j < maxCol; j++)
    {
      j1 = arr (pos2col,j, int) ;
      c = arrp (spread->colonnes,j1,COL) ;
    if (!c->hidden)
      switch (c->type)
	{
	case 0:
	  break ;
	case 'i': case 'c':
	  freeOut (" Int ") ;
	  break ;
	case 'f': 
	  freeOut (" Float ") ;
	  break ;
	case 't': 
	  freeOut (" Text ") ;
	  break ;
	case 'd': 
	  freeOut (" Date ") ;
	  break ;
	case 'k': case 'K': case 'n': case 'b':
	  for (i=0 ; i < arrayMax (spread->tableau) ; i++)
	    { 
	      lineArray = arr (spread->tableau, i, Array) ;
	      if (lineArray && j1 < arrayMax (lineArray) &&
		  !arr (lineArray, j1,SPCELL).empty)
		{ 
		  u = & (arr (lineArray, j1,SPCELL).u) ;
		  key = u->k ;
		  if (key)
		    { freeOutf (" %s ", class (key) ? className (key) : "Tag") ;
		      array (spread->oldClass,j1,int) = class (key) ;
		      goto ok ;
		    }
		}
	    }
	  array (spread->oldClass,j1,int) = -1 ;
	  freeOut (" Key ") ;
	ok:
	  break ;
	}
    }
  freeOut ("\n") ;

  return;
} /* spreadDoDumpFormat */

/*****************/

int spreadDoDump (SPREAD spread, int sep0, char style, BOOL beginTable)
     /* style 0, like 'a' but no format line */
{ 
  BOOL title = spread->showTitle ;
  int n = 0, i, maxCol = arrayMax (spread->colonnes) ;
  Array lineArray , oldLineArray = 0 ; 
  Stack s = stackCreate (80) ;
  char *href = 0 ;
  static char *href1 = 0 ;

  if (beginTable && (style == 'a' || style == 'h'))
    { spreadDoDumpFormat (spread) ;}

  if (!href1)
    {
      OBJ Clo = 0 ;
      char sp[256], *cp ;
      KEYSET clones = query (0, "Find clone strategy && Species") ;

      strcpy (sp, "human") ;
      if (keySetMax(clones))
	{
	  if ((Clo = bsCreate (keySet (clones, 0))))
	    {
	      if (bsGetData (Clo, str2tag ("Species"), _Text, &cp) && cp && *cp)
		strncpy (sp, cp, 255) ;
	      bsDestroy (Clo) ;
	    }
	}
      keySetDestroy (clones) ;

      cp = sp - 1 ;
      while (*++cp) *cp = ace_lower(*cp) ;
      href1=hprintf(0, "www.aceview.org/av.cgi?db=%s", sp) ;
    }

  if (beginTable && style == 'x' && 
      messPrompt ("Do you wish to call back to", href1,"wz"))
    {
      char *cp = freeword () ;
      messfree (href1) ;
      href1 = strnew (cp, 0) ;
      href = hprintf (0, "<a href=\"http://%s&", cp) ;
    }
  if (style == 'X') 
    {
      style = 'x' ; /* X: x but non interactive */
      if (spread->href)
	href = hprintf (0, "<a href=\"http://%s&", spread->href) ;
    }
    
  lastPrintedCol = -1 ;
  spreadReorder (spread) ;
  if (arrayExists (spread->tableau) && arrayMax (spread->tableau))
    { spread->flags = 
	arrayReCreate (spread->flags, arrayMax (spread->tableau), char) ;
      array (spread->flags, arrayMax (spread->tableau) - 1 , char) = 0 ;
 	
      for (i = 0 ; i < arrayMax (spread->tableau); i++)
	{ lineArray = arr (spread->tableau, i, Array) ;
	  if (oldLineArray && !spreadOrder (&oldLineArray, &lineArray))
	    { arr (spread->flags, i, char) = 1 ; continue ; }
	          /* Useful if some colonnes are hidden */
	  oldLineArray = lineArray ;
	  spreadAscii (spread, sep0, maxCol,lineArray, i, FALSE, style, title, n, href) ;
	  title = 0 ; /* show just once */
	  n++ ;
	}
    }
  if (style == 'C')
    {
    /*
    * We are now at the end of an encore block.  At the beginning of the
    * next encore block, we will move the cursor right one before depositing
    * data, so here we move the cursor left to get it back to column -1.
    */
    freeOutBinary ("l\001.",3);
    }
  stackDestroy (s) ;

  if (spread->exportedKeySet && keySetMax (spread->exportedKeySet) > 1)
    {
      keySetSort (spread->exportedKeySet) ;
      keySetCompress (spread->exportedKeySet) ;
    }
  messfree (href) ;
  return n ;
} /* spreadDoDump */

/******************************************************/
/******************************************************/

static char* spreadTableTypes (SPREAD spread)
{
  int i, j, j1,maxCol = arrayMax (spread->colonnes) ;
  COL *c ;
  char *types = messalloc (maxCol + 1) ;

  for (i = 0, j = 0 ; j < maxCol; j++)
    { 
      j1 = arr (pos2col,j, int) ;
      c = arrp (spread->colonnes,j1 ,COL) ;
    if (!c->hidden)
      switch (c->type)
	{
	case 'c': 
	case 'i':types[i++] = 'i' ; break ;
	case 'f':types[i++] = 'f' ; break ;
	case 'd':types[i++] = 't' ; break ; /* dates */
	case 'k': case 'K': case 'n': case 'b':
	  types[i++] = 'k' ; break ;
	case 't': types[i++] = 's' ; break ; /* text */
	}
    }
  types[i] = 0 ;

  return types ;
} /* spreadTableTypes */


TABLE *spreadToTable (SPREAD spread, AC_HANDLE h)
     /* public function */
{ 
  int ii, j, j1, j2, nLines, maxCol = arrayMax (spread->colonnes) ;
  Array lineArray , oldLineArray = 0 ; 
  TABLE *table ;
  SPCELL * su ;
  BSunit u ;
  BOOL empty ;
  COL *c ;
  char *types  ;

  if (!arrayExists (spread->tableau) || ! (nLines = arrayMax (spread->tableau)))
    return 0 ;
 
  types = spreadTableTypes (spread) ;
  table =  tableHandleCreate (nLines, types, h) ;
  messfree (types) ;

  if (arrayExists (spread->tableau) && arrayMax (spread->tableau))
    { spread->flags = 
	arrayReCreate (spread->flags, arrayMax (spread->tableau), char) ;
      array (spread->flags, arrayMax (spread->tableau) - 1 , char) = 0 ;
 				    
      for (ii = 0 ; ii < arrayMax (spread->tableau); ii++)
	{ lineArray = arr (spread->tableau, ii, Array) ;
	  if (oldLineArray && !spreadOrder (&oldLineArray, &lineArray))
	    { arr (spread->flags, ii, char) = 1 ; continue ; }
	          /* Useful if some colonnes are hidden */
	  oldLineArray = lineArray ;	

	  for (j = 0, j2 = 0 ; j < maxCol; j++)
	    {
	      j1 = arr (pos2col,j, int) ;
	      c = arrp (spread->colonnes,j1 ,COL) ;
	      if (!c->hidden)
		{
		  su = arrp (lineArray, j1,SPCELL) ;
		  u = su->u ;
		  empty = arr (lineArray, j1,SPCELL).empty ;
		  if (!empty && c->showType == SHOW_MULTI)
		    tableSetString (table, ii, j2, stackText (c->text, u.i)) ;
		  else if (!empty) switch (c->type)
		    {
		    case 'c': tableInt (table,ii,j2) = u.i ; break ;
		    case 'i':
		      {
			int xi = su->z ;
			long long int lli = su->z ;
			if (xi == lli)
			  tableInt (table,ii,j2) = xi ;
			else
			  {
			    char buf[128] ;
			    sprintf (buf, "%lld", lli) ;
			    tableSetString (table, ii, j2, buf) ;
			  }
			break ;
		      }
		    case 'f': 
		      {
			int xf = su->z ;
			if (xf == su->z)
			  tableFloat (table,ii,j2) = xf ; 
			else
			  {
			    char buf[128] ;
			    sprintf (buf, "%lg", su->z) ;
			    tableSetString (table, ii, j2, buf) ;
			  }
			break ;
		      }
		    case 'd': tableDate (table,ii,j2) = u.time ; break ;
		    case 'k': case 'K': case 'n': case 'b':
		       tableKey (table,ii,j2) = u.k ; break ;
		    case 't':  /* text */
		      tableSetString (table, ii, j2, stackText (c->text, u.i)) ;
		      break ;
		    }
		  else switch (c->type)
		    {
		    case 'i': case 'c': tableInt (table,ii,j2) = UT_NON_INT ; break ;
		    case 'f': tableFloat (table,ii,j2) = UT_NON_FLOAT ; break ;
		    case 'd': tableDate (table,ii,j2) = u.time ; break ;
		    case 'k': case 'K': case 'n': case 'b': case 't':
		      tableKey (table,ii,j2) = 0 ; break ;
		    }
		  j2++ ;
		}
	    }
	}
    }
  return table ;
}

/*****************/

/*
TABLE *spreadGet (KEY key)
{
  OBJ obj = bsCreate (key) ;
  TABLE *table = 0 ;

  if (obj) 
    {
      tableDestroy (table) ;
      bsDestroy (obj) ;
    }
  return 0 ;
}
*/

static void spreadStorePrecomputation (SPREAD spread) 
{
  OBJ obj = 0 ;
  KEY _Precompute ;

  lexaddkey ("Precompute", &_Precompute, 0) ;
  if (!spread->tmKey || !spread->tKey ||
      ! (obj = bsUpdate (spread->tmKey)))
    return ;
  spread->table = spreadToTable (spread, 0) ;
  if (spread->table) tableStore (spread->tKey, spread->table) ;

  bsAddKey (obj, _Precompute,spread->tKey) ; 
  bsSave (obj) ;
}


KEYSET spreadGetPrecalculatedKeySet (SPREAD spread, KEY key, char *cr)
     /* public function */
{	
  KEY kk1 ; KEYSET ks ; int j ;
  OBJ obj = 0 ;
  TABLE *table ; 
  char *fullname = 0 , *cp, *cq ;
  KEY _Precompute ;

  lexaddkey ("Precompute", &_Precompute, 0) ;

  if (pickType (key) != 'B') return 0 ;
  obj = bsCreate (key) ;
  if (!obj || 
      !bsFindTag (obj, _Precompute))
    {
      spread->precomputing = FALSE ;
      bsDestroy (obj) ;
      return 0 ;
    }

  spread->precomputing = TRUE ;

  fullname = messalloc (strlen (name (key)) + (cr ? strlen (cr) : 0) + 8 ) ;

  cp = fullname ; cq = name (key) ;
  while ((*cp++ = *cq++)) ; 
  if (cr) { cp-- ; *cp++ = '#' ;}
  cq = cr ? cr : "" ; 
  while (*cq == ' ' && *cq) cq++ ;
  while (*cq) {
    if (*cq == ' ') { *cp++ = '#' ; while (*cq == ' ') cq++ ; }
    *cp++ = *cq++ ;
  }
    
  lexaddkey (fullname, &kk1, _VTableResult) ;

  table = tableGet (kk1) ; spread->tKey = kk1 ; spread->tmKey = key ;
  if (!table)
    return 0 ;
  spread->precomputing = FALSE ; /* we have it already */
  spread->table = table ;
  /* if possible return keys */
  for (j=0 ; j < table->ncol ; j++)
    if ((tableTypes (table))[j] == 'k')
      return keySetCopy (table->col[j]) ;
  /* return non empty set */
  ks = keySetCreate () ;
  keySet (ks, 0) = key ;

  return ks ;
} /* spreadGetPrecalculatedKeySet */

/******************** eof **********************************/
 
