/*  File: queryexe.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 *         Used anywhere for non interactive executions.       
 * Exported functions:
         query, queryKey, queryGrep                           
 * HISTORY:
 * Last edited: Apr 12 12:11 2002 (edgrif)
 * * Oct  7 20:06 1993 (mieg): Introduced set operations {} delimited
 * * Oct  5 20:06 1993 (mieg): Introduced arithmetic operations [] delimited
 * * Nov 18 15:17 1992 (mieg): Changed # precedence and pushpop sub-objs
 * * May 26 17:21 1992 (mieg): queryIsA  for dynamic subclass system
 * Created: Tue Nov  5 16:53:55 1991 (mieg)
 *-------------------------------------------------------------------
 */


/* $Id: queryexe.c,v 1.48 2020/05/30 16:50:31 mieg Exp $ */

#define CONDITION_DEFINED
typedef struct condition *COND ;

#include "acedb.h"
#include "client.h"
#include "bindex.h"
#include "bitset.h"

extern KEYSET cacheKeySet (void) ;

 /* Query operators, the reordering loop stops on TAG */
typedef enum { /*  0 */ NULLQOP, SUBQUERY, SUBEXP, NAMEDSET, 
	       /*  4 */ TEXT, CALCUL, NUMBER, TIME, TAG, 
	       /*  9 */ FINALPIPE, PIPEQUERY, WHERE,
	       /* 12 */ FIND, FROM, FOLLOW_TAG, WEBQUERY, NEIGHBOURS, GREP, EXPAND,
	       /* 19 */ OR, XOR, AND, NOT, 
	       /* 23 */ SETELSE, SETDEFINE, SETOR, SETXOR, SETAND, SETMINUS, 
	       /* 29 */ COMPOSITE_TAG, IS, 
	       /* 31 */ LT, LE, GT, GE, NEQ, EQ, LIKE,
	       /* 38 */ PLUS, MINUS, MULTIPLY, DIVIDE, 
	       /* 42 */ COUNT, AVG, SUM, MIN, MAX, CLASS, VECTOR_TAG 
	       /* 49 */  } QOP ;

char *qName[] = {
               /*  0 */ "NULLQOP", "SUBQUERY", "SUBEXP", "NAMEDSET", 
	       /*  4 */ "TEXT", "CALCUL", "NUMBER", "TIME", "TAG",
	       /*  9 */ "FINALPIPE", "PIPEQUERY", "WHERE",
	       /* 12 */ "FIND", "FROM", "FOLLOW", "WEBQUERY", "NEIGHBOURS","GREP", "EXPAND",
	       /* 19 */ "OR", "XOR", "AND", "NOT",
               /* 23 */ "SETELSE", "SETDEFINE", "SETOR", "SETXOR", "SETAND", "SETMINUS", 
	       /* 29 */ "COMPOSITE_TAG", "IS", 
	       /* 31 */ "LT", "LE", "GT", "GE", "NEQ", "EQ", "LIKE",
	       /* 38 */ "PLUS", "MINUS", "MULTIPLY", "DIVIDE",
	       /* 42 */ "COUNT", "AVG", "SUM", "MIN", "MAX", "CLASS", "VECTOR_TAG" 
	       /* 49 */  } ;

/* RD added numbers to help when debugging */

typedef struct conditionNode *CDT;
struct conditionNode
{ QOP qop ;
  CDT up, left, right ; /* tree structure */
  union {  KEY tag ; KEY classe ; float x ; Array a ; } n ;
  double x ;
  int mark ;   /* mark in the associated stack */
} ;


static void* MAGIC ;
struct condition
{ 
  void*   magic ;		/* == &MAGIC */
  Stack stack ;
  CDT cdt ;
} ;

COND queryFindCond = 0 ;

static BOOL doInterrupt ;
static BOOL localiseTag = FALSE ;
static Stack subObjMarks = 0 ;
AC_HANDLE subObjMarksHandle = 0 ;
static BOOL fuzzy ;
static Array subQ = 0 ; /* Array of subQuery results */
static KEYSET queryExecute (KEYSET oldSet, CDT cdt, Stack stack) ;
static BOOL checkParentheses (const char *text) ;
static void querySetTextName (const char *tname)  ;

/* mieg, july 2003,
 * webquery is an addition to the query language registered from acembly
 * with this system, tace compiles ok and the word 'webquery'
 * is tokenised as an operator if and only if the function is registered
 */
static WEBQUERYFUNC webQueryFunc = 0 ;
BOOL queryRegisterWebQuery (WEBQUERYFUNC f)
{ webQueryFunc = f ; return TRUE ; }

/**********/

static BOOL isCalcul (CDT cdt)
{
  while (cdt)
    {
      if (cdt->qop == CALCUL)
	return TRUE ;
      cdt = cdt->up ;
    }
  return FALSE ;
}

/**********/

static BOOL isNumber(QOP qop)
{
  switch (qop)
    {
    case CALCUL:
    case NUMBER:
    case TIME:
    case PLUS:
    case MINUS:
    case MULTIPLY:
    case DIVIDE:
    case COUNT:
    case AVG:
    case SUM:
    case MIN:
    case MAX:
      return TRUE ;
    default: break ;
    }
  return FALSE ;
}

/**********/

static BOOL isLocatorOp (QOP qop)
{
  switch (qop)
    {
    case TAG:
    case COMPOSITE_TAG:   /* example:: Int #type */
      return TRUE ; /* a tag or composite tag is necessarilly unique */
    case VECTOR_TAG:      /* example:: tag:2   */  
      return FALSE ; /* should be true if we could handle ALL matches
                        not just the first one
			if we fix this, this will influence calls to queryExpand
			and the usage of localiseTag
		     */
    default: break ;
    }
  return FALSE ;
}

/**************************************************************/

/* Manage CDT allocation
   hope they speed things up/prevent fragmentation
   in any case, they will help monitoring memory usage
   NB they must be paired.
*/

static AC_HANDLE CDThandle = 0 ;
static Stack freeCDTstack = 0 ;
static int nCDTused = 0, nCDTalloc = 0 ;        /* useful for debug */

static CDT CDTalloc (void)       /* self managed calloc */
{
  static int blocSize = 512 ;
  CDT p ;
  int i ;
 
  nCDTalloc++ ;

  if (!freeCDTstack)
    {
      CDThandle = ac_new_handle () ;
      freeCDTstack = stackHandleCreate (4*blocSize, CDThandle) ;
    }
  if (stackEmpty (freeCDTstack))
    { p = (CDT) halloc (blocSize * sizeof (struct conditionNode), CDThandle) ;
      for (i = blocSize ; i-- ; ++p)
        push (freeCDTstack,p,CDT) ;
/*       printf ("Adding %d to CDT free list\n",blocSize) ; */
      blocSize *= 2 ;
    }
  p = pop (freeCDTstack,CDT) ;
  memset (p, 0, sizeof (struct conditionNode)) ;
  ++nCDTused ;
  return p ;
}

static void CDTfree (CDT cdt)
{
  push (freeCDTstack,cdt,CDT) ;
  --nCDTused ;
}

void CDTstatus (int *used, int *alloc)
{ *used = nCDTused ; *alloc = nCDTalloc ;
}

void CDTShutDown (void)
{ 
  ac_free (CDThandle) ;
}

/**************************************************************/

/* Manage COND allocation
   hope they speed things up/prevent fragmentation
   in any case, they will help monitoring memory usage
   NB they must be paired.
*/

static AC_HANDLE condHandle = 0 ;
static Stack freeCONDstack = 0 ;
static int nCONDused = 0, nCONDalloc = 0 ;        /* useful for debug */

static COND CONDalloc (void)       /* self managed calloc */
{
  static int blocSize = 512 ;
  COND p ;
  int i ;
 
  nCONDalloc++ ;

  if (!freeCONDstack)
    {
      condHandle = ac_new_handle () ;
      freeCONDstack = stackHandleCreate (4*blocSize, condHandle) ;
    }
  if (stackEmpty (freeCONDstack))
    { 
      p = (COND) halloc (blocSize * sizeof (struct condition), condHandle) ;
      for (i = blocSize ; i-- ; ++p)
        push (freeCONDstack,p,COND) ;
/*       printf ("Adding %d to COND free list\n",blocSize) ; */
      blocSize *= 2 ;
    }
  p = pop (freeCONDstack,COND) ;
  memset (p, 0, sizeof (struct condition)) ;
  p->magic = &MAGIC ;
  ++nCONDused ;
  return p ;
}

static void CONDfree (COND cond)
{ cond->magic = 0 ;
  push (freeCONDstack,cond,COND) ;
  --nCONDused ;
}

void condShutDown (void)
{
  ac_free (condHandle) ;
}

void CONDstatus (int *used, int *alloc)
{ *used = nCONDused ; *alloc = nCONDalloc ;
}

/**************************************************************/
       /* Recursively destroys a tree of conditions */

static void cdtDestroy(CDT c)
{
  if (c)
    {
      if (c->left)
	cdtDestroy(c->left) ;
      if (c->right)
	cdtDestroy(c->right) ;
      c->left = c->right = 0 ;  /* stops eventual loop */
      CDTfree(c) ;
    }
}

/**************************************************************/
       /* Call with 0 if you want to rescan the same text */
       /* itemizes a text considering quotes and nested parentheses */
static QOP cdtItem (char *text, char **item, int context)
{
  static char *previousText = 0 , *cp , backup = 0 ;
  int np = 0 , ncb ;
  static  QOP qop = 0 ;
  
  *item = 0 ;
  if (!text || (previousText != text))
    { cp = text ; previousText = text ; qop = 0; }
  else 
    {
      *cp = backup ;
      switch (qop)
	{   /*skip the former end delimiter */
	case SUBEXP: case SUBQUERY: case TEXT: case CALCUL:
	  if (*cp) 
	    cp++ ;
	  break ;
	default:
	  break ;
	}
    }

  qop = 0 ;
  if (!text)
    return 0 ;

  while (TRUE) switch(*cp)
    {
    case 0:
      goto done ;

    case '(':
      if (qop)
	goto done ;
      *item = ++cp ; /* inside the parenthese */
      np = 1 ;
      while (np)
	switch (*cp++)
	  {
          case 0:
	    messout("Unbalanced parenthesis:\n %s",text) ;
	    return NULLQOP ;
	  case '(':
	    np++ ;
	    break ;
	  case ')':
	    np--; 
	    break ;
	  }
      cp-- ; /* The paranthesis will be masked */
      qop = SUBEXP ;
      goto done ;

    case '{':
      if (qop)
	goto done ;
      *item = ++cp ; /* inside the parenthese */
      ncb = 1 ;
      while (ncb)
	switch (*cp++)
	  {
          case 0:
	    messout ("Unbalanced curly bracket:\n %s",text) ;
	    return NULLQOP ;
	  case '{':
	    ncb++ ;
	    break ;
	  case '}':
	    ncb--; 
	    break ;
	  }
      cp-- ; /* The paranthesis will be masked */
      qop = SUBQUERY ;
      goto done ;

    case '[':
      if (qop)
	goto done ;
      *item = ++cp ; /* inside the parenthese */
      ncb = 1 ;
      while (ncb)
	switch (*cp++)
	  {
          case 0:
	    messout ("Unbalanced square bracket:\n %s",text) ;
	    return NULLQOP ;
	  case '[':
	    ncb++ ;
	    break ;
	  case ']':
	    ncb--; 
	    break ;
	  }
      cp-- ; /* The paranthesis will be masked */
      qop = CALCUL ; 
      goto done ;

    case '\"':
      if (qop)
	goto done ;
      *item = ++cp ;
      while (*cp && !(*cp == '\"' && *(cp - 1) != '\\')) cp++ ;
      if (!*cp)
	{ messout("Unbalanced double quote:\n %s",text) ;
	  return NULLQOP ;
	}
           /* The final double quote will be masked */
      qop = TEXT ;
      goto done ;

    case ')':
      messout("Unexpected closing parenthesis:\n %s",text) ;
      return NULLQOP ;

    case '}':
      messout("Unexpected closing curly bracket:\n %s",text) ;
      return NULLQOP ;

    case ']':
      messout("Unexpected closing square bracket:\n %s",text) ;
      return NULLQOP ;

    case ' ': case '\t': case '\n': case ',' :  
      if (qop)
	goto done ;
      cp++ ; /* otherwise ignore */
      break ;
    
    case '!':
      if (qop)
	goto done ;
      cp++ ;
      if (*cp == '=')
	{ cp++ ;
	  qop = NEQ ;
	}
      else
	qop = NOT ;
      goto done ;

    case ';':
      if (qop)
	goto done ;
      cp++ ;
      qop = PIPEQUERY ;
      goto done ;

    case '~':
      if (qop)
	goto done ;
      cp++ ;
      qop = LIKE ;
      goto done ;

    case '^':
      if (qop)
	goto done ;
      cp++ ;
      qop = XOR ;
      goto done ;

    case '<':
      if (qop)
	goto done ;
      cp++ ;
      if (*cp == '=')
	{ cp++ ;
	  qop = LE ;
	}
      else
	qop = LT ;
      goto done ;

    case '>':  /* Multiple meaning, >Tag with nothing on the left
		  means FOLLOW tag */
      if (qop)
	goto done ;
      cp++ ;
      if (*cp == '=')
	{ cp++ ;
	  qop = GE ;
	}
      else if (*cp == '$')
	{ cp++ ;
	  qop = FROM ;
	}
      else if (*cp == '?')
	{ cp++ ;
	  qop = FIND ;
	}
      else 
	qop = GT ;
      goto done ;

    case '$':  /* Set operations */
      if (qop)
	goto done ;
      cp++ ;
      if (*cp == '&')
	{ cp++ ;
	  qop = SETAND ;
	}
      else if (*cp == '|')
	{ cp++ ;
	  qop = SETOR ;
	}
      else if (*cp == '^')
	{ cp++ ;
	  qop = SETXOR ;
	}
      else if (*cp == '-')
	{ cp++ ;
	  qop = SETMINUS ;
	}
      else if (*cp == '/')
	{ cp++ ;
	  qop = SETELSE ;
	}
      else if (*cp == '=')
	{ cp++ ;
	  qop = SETDEFINE ;
	}

      if (qop)   /* recognised set operator */
	goto done ;
                 /* otherwise interpreted as a named set */
      if (!*item)
	{ qop = NAMEDSET ; /* to be reexamined when the word is complete */
	  *item = cp - 1 ;
	}
      break ;
    
    case '=':   /* I treat = and == the same way */
      if (qop)
	goto done ;
      cp++ ;
      if (*cp == '=')
	{ cp++ ;
	}
      qop = EQ;
      goto done ;

    case '&':   /* I treat & and && the same way */
      if (qop)
	goto done ;
      cp++ ;
      if (*cp == '&')
	{ cp++ ;
	}
      qop = AND ;
      goto done ;

    case '|':   /* I treat | and || the same way */
      if (qop)
	goto done ;
      cp++ ;
      if (*cp == '|')
	{ cp++ ;
	}
      qop = OR ;
      goto done ;

    case '#':
      if (qop)
	goto done ;
      cp++ ;
      qop = COMPOSITE_TAG ;
      goto done ;

    case ':':  
      if (qop)
	goto done ;
      cp++ ;
      qop = VECTOR_TAG ;
      goto done ;

    case '+': case '-': case '*': case '/':
      if (context == CALCUL)
	switch (*cp)
	  {
	  case '+':  
	  case '-':  
	    if (qop)
	      goto done ;
	    { 
	      char *cq = cp ;
	      int signe = 1 ;
	      while (*cq)
		switch (*cq++)
		  {
		  case '+':
		    cp = cq ;
		    continue ;
		  case '-':
                    signe = - signe ;
		    cp = cq ;
		    continue ;
		  case ' ':
		    continue ;
		  default:
		    goto signDone ;
		    break ;
		  }
	    signDone:
	      qop = (signe == 1 ? PLUS : MINUS) ;
	    }
	    goto done ;
	    
	  case '*':  
	    if (qop)
	      goto done ;
	    cp++ ;
	    qop = MULTIPLY ;
	    goto done ;

	  case '/':  
	    if (qop)
	      goto done ;
	    cp++ ;
	    qop = DIVIDE ;
	    goto done ;

	  }
     else ;  /* ATTENTION: else fall thru  onto default */

    default:
      if (!*item)
	{ qop = TAG ; /* to be reexamined when the word is complete */
	  *item = cp ;
	}
      cp++ ;
      break ;
    }
 done:
  backup = *cp ;
  *cp = 0 ;
  return qop ;
}

/**************************************************************/
/**************************************************************/
       /* Recursive construction of a tree of conditions */
/**************************************************************/

static KEY goodClass(KEY classe)
{ 
  KEY subClass = classe ; unsigned char c ;
  
  if (pickIsComposite(classe))
    return classe;
  
  if (!pickIsA(&classe, &c))
    return 0 ;
  
  return c ? subClass : classe ;  /* return super class only if mask == 0 */
}

       /* Scans the cdt Tree
	  looking for the occurence of a query operator 
	  in order of increasing priority
	  */
static BOOL cdtReorder(CDT *top, Stack s)
{
  CDT c, cc, new, ca, cb, cr, cl ;
  QOP q ;
  BOOL applyLeft, applyRight, setWhere ;

     /* Progressivelly transform the linear
        expression into a tree by moving every recogised
        expression in the left branch of as SUBEXP
	*/

     /* A leading GT is a follow */
   c = *top ;
   if (!c) 
     return FALSE;
   while (c)
     {
       if (c->qop == GT &&
	   (!c->up ||
	     c->up->qop == PIPEQUERY ||
	     (c->up->qop == SUBQUERY && c == c->up->left)
	   ))
	   c->qop = FOLLOW_TAG ;
       c = c->right ;
     }

  c = *top ;
  while(c)
    {
      if (c->left)
	if (!cdtReorder(&(c->left), s)) 
	  return FALSE ;
      c = c->right ;
    }

     /* Now we reorder the d1 side recursively 
      * We find highest operator and make of it a subexpression 
      * In case of equality I work left to right as usual
      * i.e. the rightest binds more strongly
      * 
      * This is wrong for MINUS which is not associative
      * so i change to leftest binds more strongly for MINUS
      */

  while (TRUE)
    {
      cc = c = *top ;
      q = c->qop ;
      while(c->right)
	{ 
	  c = c->right ;
	  /*  if (c->qop >=q)  rightest of given precedence */
	  if (c->qop > q ||  /* leftest of given precedence */
	      (q != MINUS && q == c->qop)) /* rightest if ! MINUS */
	    { q = c->qop ;
	      cc = c ;
	    }
	}

      if (q <= FOLLOW_TAG)  /* only one of those between any pair of ;; */
	{ int n = 0 ;
	  c = cc ;
	  while (c)
	    { switch (c->qop)
	      {
	      case PIPEQUERY:
		n = 0 ; break ;
	      case FIND: case FOLLOW_TAG: case FROM:
		if (n)
		  { messerror("Two occurences of FIND/FOLLOW/FROM in %s", 
			      stackText(s, (*top)->mark)) ;
		  return FALSE ;
		  }
		break ;
	      default: break ;
	      }
	    c = c->up ;
	    /* c = c->right ;  leftest of given */
	    }
	}

      if (q <= PIPEQUERY ) /* was <= TAG  || ((q <= FOLLOW_TAG) && cc->right))   */
	{ /* now search for PIPEQUERY */ 
	  c = *top ; cc = 0 ;
	  while(c)
	    { 
	      if (c->qop  == PIPEQUERY) /* rightest of given precedence */
		{ cc = c ; q = PIPEQUERY ;}
	      c = c->right ;
	    }
	  if (!cc) break ;
	}

      setWhere = FALSE ;
      applyLeft = FALSE ; applyRight = TRUE ; /* default is unary */
      if (q == NEIGHBOURS)  /* Solitary operator */
	{ cc->qop = NEIGHBOURS ;
	  setWhere = TRUE ;
	  applyRight = FALSE ;
	}
      else if (q == GREP)  /* Unary operator */
	{ char *cp ;
	  if (!cc->right)
	    { messerror("Nothing to the right of %s: %s",
			qName[cc->qop], stackText(s, cc->mark)) ;
	      return FALSE ;
	    }
	  cp = stackText(s,cc->right->mark) ;
	  while (*cp=='*' || *cp == '?') cp++ ;
	  if (strlen(cp) < 3)
	     { messout("The text to the right of GREP should have at least 3 characters %s: %s",
			qName[cc->qop], stackText(s, cc->mark)) ;
	      return FALSE ;
	    }
	}
      else if (q == WEBQUERY)  /* Unary operator */
	{ char *cp ;
	  if (! webQueryFunc)
	    {
	      messout ("WEBQUERY is not available in this version of acedb, sorry") ;
	      return FALSE ;
	    }
	  if (!cc->right)
	    { 
	      messout ("Nothing to the right of %s: %s",
			qName[cc->qop], stackText(s, cc->mark)) ;
	      return FALSE ;
	    }
	  cp = stackText(s,cc->right->mark) ;
	  while (*cp=='*' || *cp == '?') cp++ ;
	  if (strlen(cp) < 3)
	     { 
	       messout("The text to the right of WEBQUERY should have at least 3 characters %s: %s",
		       qName[cc->qop], stackText(s, cc->mark)) ;
	       return FALSE ;
	     }
	}
      else if (q == NOT)  /* Unary operator */
	{ if (!cc->right)
	    { messout("Nothing to the right of %s: %s",
			qName[cc->qop], stackText(s, cc->mark)) ;
	      return FALSE ;
	    }
	}
      else if (q == IS)  /* Unary operator */
	{ if (!cc->right)
	    { messerror("Nothing to the right of %s: %s",
			qName[cc->qop], stackText(s, cc->mark)) ;
	      return FALSE ;
	    }
          if (cc->right->qop != TEXT &&
	     cc->right->qop < LT &&
	     cc->right->qop > LIKE)
	    { messout("IS must be followed by a text or a comparator: %s",
		      stackText(s,cc->right->mark)) ;
	      return FALSE ;
	    }
	}
      else if (q == CLASS)  
	{ if (!cc->right)
	    { messerror("Nothing to the right of %s: %s",
			qName[cc->qop], stackText(s, cc->mark)) ;
	      return FALSE ;
	    }
	  cr = cc->right ;
	  ca = cc->up ; cb = cr->right ;          

	    { /* Autocomplete to a Class */
	      KEY classe = 0;
	      
	      if (lexword2key(stackText(s,cr->mark), &classe, _VClass))
		{
		  cr->n.classe = goodClass(classe) ;
		  cr->qop = CLASS ;
		}
	      else 
		while (lexNext(_VClass,&classe))
		  if (pickMatch
		      (name(classe),
		       stackText(s,cr->mark) ))
		    { 
		      cr->n.classe = goodClass(classe) ;
		      if (strcasecmp(stackText(s,cr->mark),
					name(classe)))
			messout("Autocompleting CLASS %s to CLASS %s",
				stackText(s,cr->mark),
				name(classe)) ;
		      
		      break ;
		    }
	      if (!classe)
		{ messout( "CLASS <class> Missing or unrecognised <class> at %s: %s ",
			   qName[cc->qop], stackText(s,cr->mark)) ;
		  return FALSE ;
		}
	    }
	     /* now transform it into a subexp */
	            
	  if (cr->qop != CLASS)
	    { messout( "Missing or unrecognised Class after %s: %s",
			qName[cc->qop], stackText(s, cc->mark)) ;
	      return FALSE ;
	    }
	}
      else if (q == FIND )  /* FIND Class */
	{ if (!cc->right)
	    { messerror("Nothing to the right of %s: %s",
			qName[cc->qop], stackText(s, cc->mark)) ;
	      return FALSE ;
	    }
          cr = cc->right ;
	  ca = cc->up ; cb = cr->right ;          

	    { /* Autocomplete to a Class */
	      KEY classe = 0;
	      
	      if (lexword2key(stackText(s,cr->mark), &classe, _VClass))
		{ 
		  cr->n.classe = goodClass(classe) ;
		  cr->qop = CLASS ;
		}
	      else 
		while (lexNext(_VClass,&classe))
		  if (pickMatch
		      (name(classe),
		       stackText(s,cr->mark) ))
		    { cr->qop = CLASS ;
		      cr->n.classe = goodClass(classe) ;
		      if (strcasecmp (stackText(s,cr->mark),name(classe)))
			messout("Autocompleting FIND %s to FIND %s",
				stackText(s,cr->mark),
				name(classe)) ;
		      break ;
		    }
	    }
	   
	  if (cr->qop != CLASS)
	    { messout( "Missing or unrecognised Class #%s# after FIND",
		      stackText(s,cr->mark)) ;
	      return FALSE ;
	    }
	  if (cr->left)
	    { messout( "Sorry, parsing error after a FIND in:  %s", 
		       stackText(s, cc->mark)) ;
	      return FALSE ;
	    }
	  setWhere = TRUE ;
	}
      else if (q == EXPAND)  /* EXPAND a keyset, may 15 2003, */
	{ 
	  applyRight = FALSE ;
	}
      else if (q == FOLLOW_TAG)  /* FOLLOW tag */
	{ 
	  if (!cc->right ||
	      (cc->right && cc->right->qop == PIPEQUERY))
	    { 
	      messerror("Nothing to the right of %s: %s",
			qName[cc->qop], stackText(s, cc->mark)) ;
	      return FALSE ;
	    }
	  else
	    {
	      cr = cc->right ;
	      ca = cc->up ; cb = cr->right ;          
	      
	      if (cr->qop == TEXT)
		{ /* Autocomplete to a tag */
		  KEY t = 0;
		  while(lexNext(0,&t))
		    if (pickMatch(name(t),
				  stackText(s,cr->mark)))
		      { cr->qop = TAG ;
		      cr->n.tag = t ;
		      if (strcasecmp (stackText(s,cr->mark),name(t)))
			messout("Autocompleting >%s to >%s",
				stackText(s,cr->mark),
				name(t)) ;
		      break ;
		      }
		}
	      
	      if (!isLocatorOp (cr->qop) && 
		  !(cr->qop == SUBEXP && 
		    cr->left && isLocatorOp (cr->left->qop)))
		{ messout("Sorry, unrecognized expression to the right of %s: %s\n%s%s",
			  qName[cc->qop], stackText(s, cc->mark),
			  "// Should be FOLLOW tag or FOLLOW tag#tag\n",
			  "// Please use table maker for more complex operations\n") ;
		return FALSE ;
		}  
	      setWhere = TRUE ;
	    }
	}
      else if (q == COUNT || q == AVG || q == SUM || q == MIN || q == MAX) 
	/* COUNT <locator> */
	{ cr = cc->right ;
          if (!cr)
	    {
	      messerror("Nothing to the right of %s: %s",
			qName[cc->qop], stackText(s, cc->mark)) ;
	      return FALSE ;
	    }
	  if (cr->qop == TEXT)
	    { /* Autocomplete to a tag */
	      KEY t = 0;
	      while(lexNext(0,&t))
		if (pickMatch(name(t),
			     stackText(s,cr->mark)))
		  { cr->qop = TAG ;
		    cr->n.tag = t ; 
		    if (strcasecmp (stackText(s,cr->mark),name(t)))
		      messout("Autocompleting ... %s to ... %s",
			      stackText(s,cr->mark),
			      name(t)) ;
		    break ;
		  }
	    }

	  if (cr->qop != TAG && cr->qop != SUBEXP &&
	     !(cc->qop == COUNT && cr->qop == SUBQUERY))
	    { messout( "Unrecognised expression after %s: %s",
		        qName[cc->qop], stackText(s, cc->mark)) ;
	      return FALSE ;
	    }  
	}
      else if (q == PIPEQUERY) /* optionally binary operator */
	{  
	  if (cc->left)  /* was buried automatically */
	    {   
	      new = CDTalloc() ;
	      new->left = cc->left ; new->right = cc->right ; new->up = cc ;
	      cc->left = new ; cc->right = 0 ;
	      if (new->left) new->left->up = new ;
	      if (new->right) new->right->up = new ;
	      new->qop = cc->qop ; cc->qop = SUBEXP ;
	      continue ;
	    }
	  if (cc->up && cc->up->qop != PIPEQUERY)
	    applyLeft = TRUE ;
	}

      else  if (cc->qop == PLUS || cc->qop == MINUS) /* optionally unary */
	{ 
	  if (cc->up && 
	      cc->up->qop != PIPEQUERY &&
	      cc->up->qop != PLUS &&
	      cc->up->qop != MINUS &&
	      cc != cc->up->left) /* buried as a ()  or {} or [] */
	    applyLeft = TRUE ;

	  if (!cc->right || cc->right->qop == PIPEQUERY) 
	    { messout("Missing expression to the right of %s: %s", 
		      qName[cc->qop], stackText(s, cc->mark)) ;
	    return FALSE ;
	    }
	}
      else /* binary operator */
	{ 
	  if (cc->up && cc->up->qop != PIPEQUERY &&
	      cc != cc->up->left) /* buried as a ()  or {} or [] */
	    applyLeft = TRUE ;
	  else
	    { messout("Missing expression to the left of %s: %s", 
		      qName[cc->qop], stackText(s, cc->mark)) ;
	    return FALSE ;
	    }

	  if (!cc->right || cc->right->qop == PIPEQUERY) 
	    { messout("Missing expression to the right of %s: %s", 
		      qName[cc->qop], stackText(s, cc->mark)) ;
	    return FALSE ;
	    }
	  if (cc->qop == VECTOR_TAG && /* FOLLOW tag */
	      cc->right->qop != NUMBER)
	    { 
	      messout("Sorry, unrecognized expression to the left of %s: %s\n",
		      qName[cc->qop], stackText(s, cc->mark),
		      "// Should be tag:number#tag\n",
		      "// Please use table maker for more complex operations\n") ;
	      /* doInterrupt = TRUE ; */
	      return FALSE ;
	    }  
	}
    
  
  /* turn it into a SUBEXP */
  
      new = CDTalloc() ;
      
      cl = cr = ca = cb = 0 ; 
      if (applyRight)
	{
	  cr = cc->right ; 
	  if (cr)
	    { 
	      cb = cr ? cr->right : 0 ;  
	      new->right = cr ; cr->up = new ; cr->right = 0 ;
	      cc->right = cb ; if (cb) cb->up = cc ;
	    }
	  
	}
      if (applyLeft)
	{
	  cl = cc->up ;
	  if (*top == cl)
	    *top = cc ;
	  if (cl)
	    { 
	      ca = cl->up ;
	      new->left = cl ; cl->up = new ; cl->right = 0 ;

	      if (ca)  /* possibly alreay points on cc */
		{
		  if (ca->left == cl) 
		    ca->left = cc ;
		  else if (ca->right == cl)  
		    ca->right = cc ;
		}
	      cc->up = ca;
	    }
	  
	}


      new->qop = cc->qop ; 
      new->up = cc ; 
      if (cc->left)
	{ messout("Sorry, kernel error while parsing %s: %s" ,
		  qName[cc->qop], stackText(s, cc->mark)) ;
	  return FALSE ;
	}
      cc->left = new ; 
      cc->qop = setWhere ? PIPEQUERY : SUBEXP ;
    }
  return TRUE ;
}

/***********************/

       /* construct a chain of all operators at a given level */
static CDT cdtConstructLevel (char *text, Stack s, int context)
{
  CDT cdt = 0, top = 0, next ;
  char dummy ;
  char *item, *cp ;
  QOP qop ;
  KEY tag;
  double xx ;

  if (!text || !*text)
    return 0 ;
  cdtItem(0,&item, context) ;
  while ((qop = cdtItem (text, &item, context)))
    {
      next = CDTalloc() ;
      if (!top)
	top = next ;
      if (cdt)
	cdt->right = next ;
      next->up = cdt ;
      cdt = next ;
      cdt->qop = qop ;
      switch(qop)
	{
	case TAG:
	  if (!strcasecmp(item,"FOLLOW"))
	    { cdt->qop = FOLLOW_TAG ;
	      break ;
	    }
	  if (!strcasecmp(item,"EXPAND"))
	    { cdt->qop = EXPAND ;
	      break ;
	    }
	  else if (!strcasecmp(item,"FIND"))
	    { cdt->qop = FIND ;
	      break ;
	    }
	  else if (!strcasecmp(item,"FROM"))
	    { cdt->qop = FROM ;
	      break ;
	    }
	  else if (!strcasecmp(item,"WHERE"))
	    { cdt->qop = PIPEQUERY ; /* synonym of WHERE */
	      break ;
	    }
	  else if (!strcmp(item,"NEQ"))
	    { cdt->qop = NEQ ;
	      break ;
	    }
	  else if (!strcmp(item,"LIKE"))
	    { cdt->qop = LIKE ;
	      break ;
	    }
	  else if (!strcmp(item,"AND"))
	    { cdt->qop = AND ;
	      break ;
	    }
	  else if (!strcmp(item,"OR"))
	    { cdt->qop = OR ;
	      break ;
	    }
	  else if (!strcmp(item,"NOT"))
	    { cdt->qop = NOT ;
	      break ;
	    }
	  else if (!strcmp(item,"CLASS"))
	    { cdt->qop = CLASS ;
	      break ;
	    }
	  else if (!strcmp(item,"IS") || !strcmp(item,"Is") ||
		  !strcmp(item, "SELF") || !strcmp(item, "Self"))
	    { cdt->qop = IS ;
	      break ;
	    }
	  else if (!strcmp(item,"COUNT"))
	    { cdt->qop = COUNT ;
	      break ;
	    }
	  else if (!strcmp(item,"AVG"))
	    { cdt->qop = AVG ;
	      break ;
	    }
	  else if (!strcmp(item,"MIN"))
	    { cdt->qop = MIN ;
	      break ;
	    }
	  else if (!strcmp(item,"MAX"))
	    { cdt->qop = MAX ;
	      break ;
	    }
	  else if (!strcmp(item,"SUM"))
	    { cdt->qop = SUM ;
	      break ;
	    }
	  else if (!strcmp(item,"XOR"))
	    { cdt->qop = XOR ;
	      break ;
	    }
	  else if (!strcmp(item,"SETOR"))
	    { cdt->qop = SETOR ;
	      break ;
	    }
	  else if (!strcmp(item,"SETXOR"))
	    { cdt->qop = SETXOR ;
	      break ;
	    }
	  else if (!strcmp(item,"SETAND"))
	    { cdt->qop = SETAND ;
	      break ;
	    }
	  else if (!strcmp(item,"SETMINUS"))
	    { cdt->qop = SETMINUS ;
	      break ;
	    }
	  else if (!strcmp(item,"SETELSE"))
	    { cdt->qop = SETELSE ;
	      break ;
	    }
	  else if (!strcmp(item,"SETDEFINE"))
	    { cdt->qop = SETDEFINE ;
	      break ;
	    }
	  else if (!strcasecmp(item,"NEIGHBOURS"))
	    { cdt->qop = NEIGHBOURS ;
	      break ;
	    }
	  else if (!strcasecmp(item,"GREP"))
	    { cdt->qop = GREP ;
	      break ;
	    }
	  else if (webQueryFunc && !strcasecmp(item,"WEBQUERY"))
	    { cdt->qop = WEBQUERY ;
	      break ;
	    }
	  else if (!strcmp(item,"GE"))
	    { cdt->qop = GE ;
	      break ;
	    }
	  else if (!strcmp(item,"GT"))
	    { cdt->qop = GT ;
	      break ;
	    }
	  else if (!strcmp(item,"LE"))
	    { cdt->qop = LE ;
	      break ;
	    }
	  else if (!strcmp(item,"LT"))
	    { cdt->qop = LT ;
	      break ;
	    }
	  else if (!strcmp(item,"NEXT"))
	    { cdt->n.tag = _bsRight ;
	      break ;
	    }
	  else if (!strcmp(item,"HERE"))
	    { cdt->n.tag = _bsHere ;
	      break ;
	    }
	  /* NOTE: the order of the next two clauses is important,
	     we need to recognise numbers in preference to tagnames.
	     Tagnames which are valid numbers are illegal and flagged, 
	     but the read-models process on a models.wrm file
	     with them can still put tags with names like "10"
	     into class zero. Hence we don't look for tags
	     until after we've looked from numbers. 
	  */
	  else if (sscanf (item,"%lg%c",&xx,&dummy) == 1)
                  /* dummy trick checks no unmatched characters */
	    { cdt->qop = NUMBER;
	      cdt->x = xx ;
	      /* fall through to match 1967 with a _Text 1967 */
	    }
	  else if (lexword2key(item, &tag, _VSystem))
	    { cdt->n.tag = tag ;
	      /* fall through to match tag with a _Text */
	    }
	  else
	    { cdt->qop = TEXT;
	      /* fall through */
	    }
	case TEXT: 
	  cdt->mark = stackMark(s) ;
	  if ((cp = freeunprotect(item)))
	    pushText(s,cp) ;
	  break ;

	case SUBEXP : case NAMEDSET : case CALCUL: case SUBQUERY: 
	  cdt->mark = stackMark(s) ;
	  pushText(s,item) ;
	  break ;
	default:
	  break ;
	}
    }
  return top ;
}

/***********************/
 /* Expands subexpressions presented in parenthesis to the scanner */

static BOOL cdtExpand(CDT cdt, Stack s, int context) ;

static CDT cdtConstructExpand(char *text, Stack s, int context) 
{
 CDT cdt =  cdtConstructLevel(text,s, context) ;
 cdtExpand(cdt,s, context) ;
 return cdt ;
}

static BOOL cdtExpand(CDT cc, Stack s, int context)
{
 while(cc)
   {
     if (cc->qop == SUBEXP || cc->qop == SUBQUERY || cc->qop == CALCUL)
	 { cc->left = cdtConstructExpand(stackText(s,cc->mark), s, 
				       cc->qop == CALCUL ? CALCUL : 
				         cc->qop == SUBQUERY ? 0 : context) ;
	   if (cc->left) cc->left->up = cc ;
	 }
     cc = cc->right ;
   }
 return TRUE ;
}

/**********/

static CDT cdtRemoveSubExp(Stack s, CDT cc, BOOL *problemp) 
{
  if (cc->left)
    cc->left = cdtRemoveSubExp(s, cc->left, problemp) ;
  if (cc->right)
    cc->right = cdtRemoveSubExp(s, cc->right, problemp) ;

  if (cc->qop == VECTOR_TAG )  /* FOLLOW tag */
    {
      if (!cc->right)
	{ *problemp = TRUE ;
	messerror("Nothing to the right of %s: %s",
		  qName[cc->qop], stackText(s, cc->mark)) ;
	return cc ;
	}
      
      if (!cc->left)
	{ *problemp = TRUE ;
	messerror("Nothing to the left of %s: %s",
		  qName[cc->qop], stackText(s, cc->mark)) ;
	return cc ;
      }
    if (cc->right->qop != NUMBER ||
	!isLocatorOp (cc->left->qop))
      { *problemp = TRUE ;
	messout("Sorry, unrecognized expression around %s: %s\n",
		qName[cc->qop], stackText(s, cc->mark),
		"// Should be tag:number\n",
		"// Please use table maker for more complex operations\n") ;
	return NULL ;
      }  
    }
  
  if (cc->qop != SUBEXP)
    return cc ;
  
  if (cc->right)
    { *problemp = TRUE ;
      messout("Warning: Please reformulate your query, the () are bizare in %s", stackText(s,0)) ;
      return cc ;
    }
  if (cc->left)
    { CDT cl = cc->left ;
      cc->left->up = cc->up ;
      CDTfree(cc) ;
      return cl ;
    }
  return cc ;
}

/**********/
       /* Destruction */
void  condDoDestroy(COND cond) 
{
  if (cond && cond->magic == &MAGIC)
    { stackDestroy (cond->stack) ;
      cdtDestroy(cond->cdt) ;
      CONDfree(cond) ;
      cond = 0 ;
    }
  if (cond && cond->magic != &MAGIC)
    invokeDebugger() ;
}

/**********/
       /* Initialisation part */
BOOL condConstruct (const char *text, COND* cdp)
{
  COND cond = 0 ;
  char *text1 ;
  while(*text == ' ')
    text++ ;
  if (!text || !*text || !checkParentheses (text))
    return FALSE ;
  
  text1 = strnew (text, 0) ;
  cond = CONDalloc() ;

  cond->stack = stackCreate(12*strlen(text)) ; /* wild big guess */
  pushText(cond->stack,text) ;
  cond->cdt = cdtConstructLevel (text1, cond->stack, 0) ;
  messfree (text1) ;
  if (cdtExpand (cond->cdt,cond->stack, 0) 
      &&  cdtReorder (&(cond->cdt),cond->stack) )
    { BOOL problem = FALSE ;
      /*      cdtDump(cond->stack, cond->cdt, 0, FALSE) ; */
      cond->cdt = cdtRemoveSubExp (cond->stack, cond->cdt, &problem) ;
      if (problem)
	{ 
	  condDestroy(cond) ;
	  *cdp = 0 ;
	  return FALSE ;
	}
      *cdp = cond ;
      return TRUE ;
    }
  else
    { condDestroy(cond) ;
      *cdp = 0 ;
      return FALSE ;
    }
}

/**********/
       /* Number comparisons */
static BOOL condCheck (double x, QOP qop, double y)
{
  double delta = x - y , a, b, d ; 
  BOOL egaux = FALSE ;
  float xf = x, yf = y ;

  a = x > 0 ? x : -x ;
  b = y > 0 ? y : -y ;
  d = delta > 0 ? delta : -delta ;
  if (
      (xf == yf) || 
      (x == 0 && y == 0) ||
      (d < a * .25e-12 &&  d < b * .25e-12)
      )
    egaux = TRUE ;

  switch(qop)
    {
    case LIKE:
      if (
	  (x == 0 && y == 0) ||
	  (20*d < a &&  20*d < b)  /* 5% error rate */
	  )
  	return TRUE ;
      return FALSE ;
    case EQ: 
      return egaux ;
    case NEQ:
      return !egaux ;
    case LE:
      return  egaux || delta < 0 ;
    case LT:
      return !egaux && delta < 0 ;
    case GE:
      return  egaux || delta > 0 ;
    case GT:
      return !egaux && delta > 0 ; ;
     default:
      messcrash("Unknown operator %d in condCheck.",qop) ;
    }
  return FALSE ; /* for compiler happiness */
}
/**********/
       /* Date comparisons remplacement de condCheckDate
	par timeComparison mhmp 22.10.98 */

/**********/
       /* Time comparisons
	  suppress that, because nearly anything is a date
          so it creates havoc in untyped string comparisons

static BOOL condCheckTime (char *cp, QOP qop, char *cq)
{ float x, y ;
  if ((x = (float) timeParse(cp)) && 
      (y = (float) timeParse(cq)))
    return condCheck (x, qop, y) ;
  else
    return FALSE ;
}
    */
     /* String comparisons */
static BOOL condCheckText (const char *cp, QOP qop, char *cq)
{
  switch(qop)
    {
    case EQ:
      return pickMatch(cp, cq)  ? TRUE : FALSE ;
    case LIKE:
      return queryRegExpMatch(cp, cq, FALSE) > 0  ? TRUE : FALSE ;
    case NEQ:
      return !pickMatch(cp, cq)  ? TRUE : FALSE ;
    case LE:
      return lexstrcmp(cp,cq) <= 0  ? TRUE : FALSE ;
    case LT:
      return lexstrcmp(cp,cq) < 0  ? TRUE : FALSE ;
    case GE:
      return lexstrcmp(cp,cq) >= 0  ? TRUE : FALSE ;
    case GT:
      return lexstrcmp(cp,cq) > 0  ? TRUE : FALSE ;
    default:
      messcrash("Unknown operator %d in condCheckText.",qop) ;
    }
  return FALSE ; /* for compiler happiness */
}

/*******************************/
static  BOOL  writeMess = TRUE ;
static const char* myName = 0 ;
static void querySetTextName (const char *cp) 
{
  messfree(myName) ;
  if (cp) myName = strnew (cp, 0) ;
}

static BOOL myNextName(KEY key, const char **cpp)
{
  if (key <= _LastC)
    {
      if (*cpp) 
	return FALSE ;
      else 
	{ *cpp = myName ; return myName ? TRUE : FALSE ; }
    }
  return nextName (key, cpp) ;
}

/************************************************************************/

static BindexFindResult condMatch(Stack s, KEY key, OBJ* objp, CDT cdt)
{
  KEY k ; BindexFindResult found ; 
  const char *ccp ;
  char *cp ;  float f ; int i = 0 ; mytime_t tm ;
  /*  QOP oper ; mhmp 22.10.98 */

  if (fuzzy && (*objp)) return BINDEX_TAG_UNCLEAR ; /* obj already open, use true search */
  if (cdt) 
  switch (cdt->qop)
    {
    case FOLLOW_TAG:
      messout("Internal error : condMatch should not parse FOLLOW") ;
      return BINDEX_TAG_ABSENT ;
    case WHERE:
      messout("Internal error : condMatch should not parse WHERE") ;
      return BINDEX_TAG_ABSENT ;
    case FIND:
      messout("Internal error : condMatch should not parse FIND") ;
      return BINDEX_TAG_ABSENT ;
    case FROM:
      messout("Sorry: set operation FROM not yet implemented") ;
      return BINDEX_TAG_ABSENT ;
    case SUBQUERY:
    case NAMEDSET: case SETDEFINE: case NEIGHBOURS:
    case SETOR: case SETAND: case SETXOR: case SETMINUS: case SETELSE:
      messout("Internal error : condMatch should not parse a SET operation") ;
      return BINDEX_TAG_ABSENT ;
    case OR:
      if (!cdt->left) return BINDEX_TAG_ABSENT ;
      if (!cdt->right) return BINDEX_TAG_ABSENT ;
      switch (condMatch(s,key,objp, cdt->left))
	{
	case BINDEX_TAG_ABSENT: 
	  return condMatch(s,key,objp, cdt->right) ;
	case BINDEX_TAG_UNCLEAR:
	  switch (condMatch(s,key,objp, cdt->right))
	    {
	    case BINDEX_TAG_ABSENT: return BINDEX_TAG_UNCLEAR ;
	    case BINDEX_TAG_UNCLEAR: return BINDEX_TAG_UNCLEAR;
	    case BINDEX_TAG_PRESENT: return BINDEX_TAG_PRESENT ;
	    }
	case BINDEX_TAG_PRESENT:
	  return BINDEX_TAG_PRESENT ;
	}
	
    case AND:
      if (!cdt->left) return BINDEX_TAG_ABSENT ;
      if (!cdt->right) return BINDEX_TAG_ABSENT ;
      switch (condMatch(s,key,objp, cdt->left))
	{
	case BINDEX_TAG_ABSENT: 
	  return BINDEX_TAG_ABSENT ;
	case BINDEX_TAG_UNCLEAR:
	  switch (condMatch(s,key,objp, cdt->right))
	    {
	    case BINDEX_TAG_ABSENT: return BINDEX_TAG_ABSENT ;
	    case BINDEX_TAG_UNCLEAR: return BINDEX_TAG_UNCLEAR;
	    case BINDEX_TAG_PRESENT: return BINDEX_TAG_UNCLEAR ;
	    }
	case BINDEX_TAG_PRESENT:
	  return condMatch(s,key,objp, cdt->right) ;
	}
      
    case XOR:
      if (!cdt->left) return BINDEX_TAG_ABSENT ;
      if (!cdt->right) return BINDEX_TAG_ABSENT ;
      switch (condMatch(s,key,objp, cdt->left))
	{
	case BINDEX_TAG_ABSENT: 
	  switch (condMatch(s,key,objp, cdt->right))
	    {
	    case BINDEX_TAG_ABSENT: return BINDEX_TAG_ABSENT ;
	    case BINDEX_TAG_UNCLEAR: return BINDEX_TAG_UNCLEAR;
	    case BINDEX_TAG_PRESENT: return BINDEX_TAG_PRESENT ;
	    }
	case BINDEX_TAG_UNCLEAR:
	  return BINDEX_TAG_UNCLEAR ;
	case BINDEX_TAG_PRESENT:
	  switch (condMatch(s,key,objp, cdt->right))
	    {
	    case BINDEX_TAG_ABSENT: return BINDEX_TAG_PRESENT ;
	    case BINDEX_TAG_UNCLEAR: return BINDEX_TAG_UNCLEAR;
	    case BINDEX_TAG_PRESENT: return BINDEX_TAG_ABSENT ;
	    }
	}

    case NOT:
      if (!cdt->right) 
	{
	  messout("Internal error : condMatch should not parse NOT with nothing to the right") ;
	  return BINDEX_TAG_ABSENT ;
	}
      switch (condMatch(s,key,objp, cdt->right)) 
	{
	case BINDEX_TAG_ABSENT: return BINDEX_TAG_PRESENT ;
	case BINDEX_TAG_UNCLEAR: return BINDEX_TAG_UNCLEAR;
	case BINDEX_TAG_PRESENT: return BINDEX_TAG_ABSENT ;
	}
      
    case SUBEXP:
      if (!cdt->left) 
	{
	  messout("Internal error : condMatch should not parse a subexp with nothing to the left") ;
	  return BINDEX_TAG_ABSENT ;
	}
      if (cdt->right) /* juxtaposition understood as AND */
	{ 
	  messout ("Internal error : missing operator at %s: %s",
		   qName[cdt->qop], stackText(s, cdt->mark)) ;
	  return BINDEX_TAG_ABSENT ; 
	}
      return condMatch(s,key,objp, cdt->left) ;
      
    case CALCUL: 
      if (!cdt->left)
	return BINDEX_TAG_ABSENT ;
      if (cdt->right) /* meaningless juxtaposition */
	return BINDEX_TAG_ABSENT ; 

      if (isNumber(cdt->left->qop))
	  switch (condMatch(s,key,objp, cdt->left))
	    {
	    case BINDEX_TAG_ABSENT: return BINDEX_TAG_ABSENT ;
	    case BINDEX_TAG_UNCLEAR: return BINDEX_TAG_UNCLEAR;
	    case BINDEX_TAG_PRESENT: cdt->x = cdt->left->x  ; return BINDEX_TAG_PRESENT ;
	    }
      else if (isLocatorOp (cdt->left->qop))
	switch (condMatch(s,key,objp, cdt->left))
	  {
	  case BINDEX_TAG_ABSENT: return BINDEX_TAG_ABSENT ;
	  case BINDEX_TAG_UNCLEAR: break ;
	  case BINDEX_TAG_PRESENT: break ; /* we still have to look inside the object */
	  }
      else
	return BINDEX_TAG_ABSENT ; /* not a number, not a locator => meaningless */

      if (fuzzy) return BINDEX_TAG_UNCLEAR ; /* quit fuzzy analysis */
      if (!*objp && !(*objp = bsCreate(key)))
	return BINDEX_TAG_ABSENT ;

      if (condMatch(s,key,objp, cdt->left) &&
	       bsGetKeyTags(*objp, cdt->left->qop == VECTOR_TAG ? _bsHere : _bsRight, &k))
	{ 
	  double x ;
	  char *endptr ;
	  if (k == _Int && bsGetData(*objp, _bsHere, k, &i)) 
	    x = i ;
	  else if (k == _Float && bsGetData(*objp, _bsHere, k, &f))
	    x = f ;
	  else if (k == _DateType && bsGetData(*objp, _bsHere, k, &tm))
	    x = tm ;
	  else if (k <= _LastC && bsGetData(*objp, _bsHere, k, &cp) && 
		   (x = strtod(cp, &endptr))) ;
	  else
	    return BINDEX_TAG_ABSENT ;
	  cdt->x = x ;
	  return BINDEX_TAG_PRESENT ;
	}
      return BINDEX_TAG_ABSENT ;

    case PLUS: case MINUS: case MULTIPLY: case DIVIDE:
      if (!cdt->right) return BINDEX_TAG_ABSENT ;
      if (fuzzy &&
	  isCalcul(cdt) &&
	  (
	   (cdt->left && !isNumber (cdt->left->qop)) ||
	   (cdt->right && !isNumber (cdt->right->qop))
	   )
	  )
	return BINDEX_TAG_UNCLEAR ;

      if (! isCalcul(cdt) &&
	  (
	   (cdt->left && !isNumber (cdt->left->qop)) ||
	   !isNumber (cdt->right->qop))
	  )
	return BINDEX_TAG_ABSENT ;

      if (cdt->left) switch (condMatch(s,key,objp, cdt->left))
	{
	case BINDEX_TAG_ABSENT: return BINDEX_TAG_ABSENT ;
	case BINDEX_TAG_UNCLEAR: return BINDEX_TAG_UNCLEAR ;
	case BINDEX_TAG_PRESENT: break ;
	}

      switch (condMatch(s,key,objp, cdt->right))
	{
	case BINDEX_TAG_ABSENT: return BINDEX_TAG_ABSENT ;
	case BINDEX_TAG_UNCLEAR: return BINDEX_TAG_UNCLEAR ;
	case BINDEX_TAG_PRESENT: break ;
	}

      switch (cdt->qop)
	{
	case PLUS:
	  cdt->x = (cdt->left ? cdt->left->x : 0 ) + cdt->right->x ;
	  break ;
	case MINUS:
	  cdt->x = (cdt->left ? cdt->left->x : 0 ) - cdt->right->x ;
	  break ;
	  /* case MULTIPLY DIVIDE are actually tested to have a cdt->left */
	case MULTIPLY:
	  cdt->x = (cdt->left ? cdt->left->x : 1 ) * cdt->right->x ;
	  break ;
	case DIVIDE:
	  if (! cdt->right->x)
	    return BINDEX_TAG_ABSENT ;
	  cdt->x = (cdt->left ? cdt->left->x : 1 ) / cdt->right->x ;
	  break ;
	default: break ;
	}
      return BINDEX_TAG_PRESENT ;
	      
    case EQ:
    case NEQ: case LT: case LE: case GT: case GE: case LIKE:
      if (!cdt->left) return BINDEX_TAG_ABSENT ;
      if (!cdt->right) return BINDEX_TAG_ABSENT ;
      if (cdt->left->qop == IS )
	{
	  ccp = 0 ;
	  found =  BINDEX_TAG_ABSENT ;
	  while (myNextName(key,&ccp))
	    if (condCheckText (ccp, cdt->qop, stackText(s,cdt->right->mark)))
	      { found = BINDEX_TAG_PRESENT ; break ; }
	  return found ;
	}

      /* = tag means = text_value of tag, so fail is no problem */
      switch (condMatch(s,key,objp, cdt->left))
	{
	case BINDEX_TAG_ABSENT: return BINDEX_TAG_ABSENT ;
	case BINDEX_TAG_UNCLEAR: return BINDEX_TAG_UNCLEAR ;
	case BINDEX_TAG_PRESENT: break ;
	}

      switch (condMatch(s,key,objp, cdt->right))
	{
	case BINDEX_TAG_ABSENT: 
	  if (cdt->right->qop != TAG) return BINDEX_TAG_ABSENT ;
	  else break ;
	case BINDEX_TAG_UNCLEAR: 
	  if (cdt->right->qop != TAG) return BINDEX_TAG_UNCLEAR ;  
	  else break ;
	case BINDEX_TAG_PRESENT: break ;
	}

      if (isNumber (cdt->left->qop) &&
	  isNumber (cdt->right->qop))
	return 
	  condCheck(cdt->left->x,cdt->qop,cdt->right->x) ? BINDEX_TAG_PRESENT : BINDEX_TAG_ABSENT ;

      if (cdt->left->qop == TEXT) /*  || class(key) < _LastN) what was that junk ? mieg feb 10 2004 */
	return 
	  condCheckText(stackText(s,cdt->left->mark), cdt->qop, 
			stackText(s,cdt->right->mark)) ? BINDEX_TAG_PRESENT : BINDEX_TAG_ABSENT ;

       
      if (cdt->left->n.tag != _bsRight &&
	  cdt->left->n.tag != _bsHere &&
	  cdt->left->n.tag < lexMax(0) && 
	  (!*objp && !bIndexFind(key,cdt->left->n.tag))) /* if *objp, i may be in a subobj */
	  return BINDEX_TAG_ABSENT ;

      if (fuzzy) return BINDEX_TAG_UNCLEAR ;
      if (!*objp && !(*objp = bsCreate(key)))
	return BINDEX_TAG_ABSENT ;

	/* RD handle _bsRight separately, else takes two steps right,
	   one from condMatch and one here
	   -- this is an ugly fix: the condMatch calls above should return
	   values, which are compared here.  i.e. the recursive process
	   should deal in values, not BOOLs
        */

      if ((cdt->left->n.tag == _bsRight && bsGetKeyTags(*objp, _bsHere, &k)) ||
	  (cdt->left->n.tag < lexMax(0) && bsGetKeyTags(*objp, cdt->left->n.tag,&k)) )
	do
	  { if (k == _Int)
	      { if (isNumber(cdt->right->qop) &&
		    bsGetData(*objp, _bsHere, k, &i) &&
		    condCheck((float) i,cdt->qop,cdt->right->x))
		  return BINDEX_TAG_PRESENT ;
	      }
	    else if (k == _Float)
	      { if (isNumber(cdt->right->qop) &&
		    bsGetData(*objp, _bsHere, k, &f) &&
		    condCheck(f,cdt->qop,cdt->right->x))
		  return BINDEX_TAG_PRESENT ; 
	      }
	    else if (k == _DateType)
	      { mytime_t x ;  /* et suppression des (float) pour x et tm */
		if ( (x = timeParse(stackText(s,cdt->right->mark))) &&
		    bsGetData(*objp, _bsHere, k, &tm) )
		  switch (cdt->qop) /* mhmp 22.10.98 */
		    {
		    case EQ:
		      return timeComparison (0, tm, x) ? BINDEX_TAG_PRESENT : BINDEX_TAG_ABSENT ;
		    case NEQ:
		      return !timeComparison (0, tm, x) ? BINDEX_TAG_PRESENT : BINDEX_TAG_ABSENT ;
		    case LT:
		      return timeComparison (-1, tm, x) ? BINDEX_TAG_PRESENT : BINDEX_TAG_ABSENT ;
		    case GE:
		      return !timeComparison (-1, tm, x) ? BINDEX_TAG_PRESENT : BINDEX_TAG_ABSENT ;
		    case GT:
		      return timeComparison (1, tm, x)  ? BINDEX_TAG_PRESENT : BINDEX_TAG_ABSENT; 
		    case LE:
		      return !timeComparison (1, tm, x) ? BINDEX_TAG_PRESENT : BINDEX_TAG_ABSENT ;
		    default:
		      if (!x && writeMess) 
			{ messout ("Sorry - not a good date") ; 
			writeMess = FALSE ;
			}
		      return BINDEX_TAG_ABSENT ;
		    }   
	      }
	    else if (k == _LongInt || k == _LongFloat)
	      { if (bsGetData(*objp, _bsHere, k, &cp))
		  { if (isNumber(cdt->right->qop))
		      { double xx ;
			if (sscanf(cp, "%lg",&xx) &&
			    condCheck(xx,cdt->qop,cdt->right->x))
			  return BINDEX_TAG_PRESENT ;
		      }
		    if (condCheckText(cp,cdt->qop,
				      stackText(s,cdt->right->mark)))
		      return BINDEX_TAG_PRESENT ;
		  }
	      }
	    else if (k <= _LastC)
	      { if (bsGetData(*objp, _bsHere, k, &cp))
		  { if (0 && isNumber(cdt->right->qop))
		      { double xx ;
			if (sscanf(cp, "%lf",&xx) &&
			    condCheck((float)xx,cdt->qop,cdt->right->x))
			  return BINDEX_TAG_PRESENT ;
		      }
		    if (condCheckText(cp,cdt->qop,
				      stackText(s,cdt->right->mark)))
		      return BINDEX_TAG_PRESENT ;
		  }
	      }
	    else             /* ordinary Key */
	      { 
		ccp = 0 ;
		while (nextName(k, &ccp))
		  { if (0 && /* 2017_10_12 do not force conversion to a number because of rounding errors */
			isNumber(cdt->right->qop))
		      { double xx ;
			if (sscanf(ccp, "%lf",&xx) &&
			    condCheck((float)xx,cdt->qop,cdt->right->x))
			  return BINDEX_TAG_PRESENT ;
		      }
		    if (condCheckText(ccp,cdt->qop,
				      stackText(s,cdt->right->mark)))
		      return BINDEX_TAG_PRESENT ;
		  }
	      }
	  } while (bsGetKeyTags(*objp,_bsDown, &k)) ; /* do while construct */
      return BINDEX_TAG_ABSENT ;
            
    case COMPOSITE_TAG:
      if (!cdt->left) return BINDEX_TAG_ABSENT ;
      if (!cdt->right) return BINDEX_TAG_ABSENT ;

      switch (condMatch(s,key,objp, cdt->left))
	{
	case BINDEX_TAG_ABSENT: return BINDEX_TAG_ABSENT ;
	case BINDEX_TAG_UNCLEAR: return BINDEX_TAG_UNCLEAR ;
	case BINDEX_TAG_PRESENT: break ;
	}

      if (fuzzy) return BINDEX_TAG_UNCLEAR ;
      if (!*objp && !(*objp = bsCreate(key)))
	return BINDEX_TAG_ABSENT ;

      if (! localiseTag)
	push(subObjMarks,bsHandleMark(*objp,0,subObjMarksHandle), BSMARK) ;
      i = 0 ;
      if (bsPushObj(*objp))
	i = condMatch(s,key,objp, cdt->right) ;
      if (! localiseTag)   /* mieg 2007_01_26, table maker Map # position no longer worked */
	bsGoto(*objp,pop(subObjMarks,BSMARK)) ;
	      
      return i ;

    case VECTOR_TAG:
      if (!cdt->left) return BINDEX_TAG_ABSENT ;
      if (!cdt->right) return BINDEX_TAG_ABSENT ;
       switch (condMatch(s,key,objp, cdt->left))
	{
	case BINDEX_TAG_ABSENT: return BINDEX_TAG_ABSENT ;
	case BINDEX_TAG_UNCLEAR: return BINDEX_TAG_UNCLEAR ;
	case BINDEX_TAG_PRESENT: break ;
	}

      if (fuzzy) return BINDEX_TAG_UNCLEAR ;

      if (!isNumber(cdt->right->qop))
	return BINDEX_TAG_ABSENT ;
      cdt->n.tag = _bsHere ;
      i = cdt->right->x ;
      if (i < 0)
	return BINDEX_TAG_ABSENT ;
      if (fuzzy)
	return BINDEX_TAG_UNCLEAR ;
      if (!*objp && !(*objp = bsCreate(key)))
	return BINDEX_TAG_ABSENT ;
      while (i--)
	if (!bsFindTag(*objp,_bsRight)) 
	  return BINDEX_TAG_ABSENT ;
      return BINDEX_TAG_PRESENT ;

    case TAG:
      /* if !fuzzy, i need to position */
      if (fuzzy)
	{
	  if (isCalcul(cdt))
	    return BINDEX_TAG_UNCLEAR ;
	  return bIndexFind(key,cdt->n.tag) ;
	}

      if (!*objp && !(*objp = bsCreate(key)))
	return BINDEX_TAG_ABSENT ;
      if (!bsFindTag(*objp,cdt->n.tag))
	return BINDEX_TAG_ABSENT ;
      if (cdt->right)  /* big blue, interpreted as big AND blue */
	return condMatch(s,key,objp, cdt->right) ;
      if (cdt->up && cdt->up->qop == COUNT)
	return BINDEX_TAG_PRESENT ;
      if (isCalcul (cdt))
	{
	  double x ; 
	  char *endptr ;

	  if (!bsGetKeyTags(*objp, cdt->n.tag,&k))
	    return BINDEX_TAG_ABSENT ;
	  if (k == _Int && bsGetData(*objp, _bsHere, k, &i)) 
	    { cdt->x = i ; return BINDEX_TAG_PRESENT ; }
	  else if (k == _Float && bsGetData(*objp, _bsHere, k, &f))
	    { cdt->x = f ; return BINDEX_TAG_PRESENT ; }
	  else if (k <= _LastC && bsGetData(*objp, _bsHere, k, &cp) && (x = strtod(cp, &endptr)))
	    { cdt->x = x ; return BINDEX_TAG_PRESENT ; }
	  else
	    return BINDEX_TAG_ABSENT ;
	}
      return BINDEX_TAG_PRESENT ;

    case COUNT:
      if (!cdt->right) 
	return BINDEX_TAG_ABSENT ;
      if (fuzzy && cdt->right && cdt->right->qop == TAG &&
	  bIndexFind(key,cdt->right->n.tag) == BINDEX_TAG_ABSENT)
	return BINDEX_TAG_ABSENT ;	  
      if (fuzzy)
	return BINDEX_TAG_UNCLEAR ;
      if (cdt->right && cdt->right->qop == SUBQUERY)
	{ 
	  Array aa =  queryKey (key, stackText(s,cdt->right->mark)) ;
	  cdt->x = aa ? keySetMax (aa) : 0 ;
	  keySetDestroy (aa) ;
	  return BINDEX_TAG_PRESENT ;
	}
      if (
	 !condMatch(s,key,objp, cdt->right) ||
	 (!*objp && !((*objp = bsCreate(key))))
	 )
	return BINDEX_TAG_ABSENT ;
      if (!bsGetKeyTags(*objp,_bsRight, &k) )
	{ cdt->x = 0 ;
	  return BINDEX_TAG_PRESENT ;
	}
      i = 0 ;
      do
	{
	  if (class(k) != _VComment)  i++ ;
	} while(bsGetKeyTags(*objp,_bsDown, &k)) ;
	
      cdt->x = i ;
      return BINDEX_TAG_PRESENT ;
	
    case AVG: case SUM: case MIN: case MAX:
      if (!cdt->right) return BINDEX_TAG_ABSENT ; 
      if (!*objp && !bIndexFind(key,cdt->right->n.tag))
	return BINDEX_TAG_ABSENT ;
      if (!bsIsTagInObj (*objp, key, cdt->right->n.tag))
	return BINDEX_TAG_ABSENT ;
      if (fuzzy)
	return BINDEX_TAG_UNCLEAR ;
      if (!*objp && !(*objp = bsCreate(key)))
	return BINDEX_TAG_ABSENT ;
      if (!bsGetKeyTags(*objp, cdt->right->n.tag,&k))
	return BINDEX_TAG_ABSENT ;

      { 
	double x, sum = 0, min = 0, max = 0 ; 
	int nn = 0 ;
	do
	{
	  if (k == _Int && bsGetData(*objp, _bsHere, k, &i)) 
	    x = i ;
	  else if (k == _Float && bsGetData(*objp, _bsHere, k, &f))
	    x = f ;
	  else if (k == _DateType && bsGetData(*objp, _bsHere, k, &tm))
	    x = tm ;
	  else if (k <= _LastC && bsGetData(*objp, _bsHere, k, &cp) && (x = (float) atof(cp))) ;
	  else
	    continue ;

	  if (!nn++)
	    min = max = sum = x ;
	  else
	    { if (x < min) min = x ;
	      if (x > max) max = x ;
	      sum += x ;
	    }
	} while(bsGetKeyTags(*objp,_bsDown, &k)) ;

	if (!nn)
	  return BINDEX_TAG_ABSENT ;
	switch(cdt->qop)
	  { 
	  case SUM:
	    cdt->x = sum ;
	    break ;
	  case AVG:
	    cdt->x = sum/nn ;
	    break ;
	  case MIN:
	    cdt->x = min ;
	    break ;
	  case MAX:
	    cdt->x = max ;
	    break ;
	  default: break ;
	  }
	return BINDEX_TAG_PRESENT ;
      }
	
    case TEXT:   
    case NUMBER:
    case TIME:	
      if (cdt->up)
	switch (cdt->up->qop)
	  {
	  case FIND: case FOLLOW_TAG: case FROM: case WHERE:
	  case PIPEQUERY: case SUBQUERY:
	  case AND: case OR: case XOR: case NOT:
	  case TEXT: case NUMBER: case TIME: case TAG:
             /* leftmost text or logical operator or 
	      * FIND Author a* will be recognised as a name selection 
	      * FOLLOW Author a* will be recognisewd as a name selection 
	      */
	    break;
	  default:
	    if (cdt->right)
	      return condMatch(s,key,objp, cdt->right) ; /* check the rest of the expression */
	    return BINDEX_TAG_PRESENT ;
	  }
	  
      ccp = 0 ;
      found = BINDEX_TAG_ABSENT ;
      while (myNextName(key,&ccp))
	if (pickMatchCaseSensitive (ccp, stackText(s,cdt->mark), 
				    pickCaseSensitive (key)))
	  { found = BINDEX_TAG_PRESENT ; break ; }
      if (found && cdt->right)
	{
	  /* messout("Interpolating AND just right of %s: %s",
		  qName[cdt->qop], stackText(s, cdt->mark)) ; */
	  return condMatch(s,key,objp, cdt->right) ; /* check the rest of the expression */
	  }
      return found ;

    case IS:  /* matches to text to its right */
      if (!cdt->right) 
	return BINDEX_TAG_ABSENT ;
      ccp = 0 ; 
      found = BINDEX_TAG_ABSENT ;
      while (myNextName(key,&ccp))
	if (pickMatchCaseSensitive(ccp, stackText(s,cdt->right->mark),
				   pickCaseSensitive (key)))
	  { found = BINDEX_TAG_PRESENT ; break ; }
      return found ;

    case CLASS:
      if (!cdt->right) 
	return BINDEX_TAG_ABSENT ;
      {
	Array a = pickIsComposite(cdt->right->n.classe);
	int j;
	
	if (!a)
	  return 
	    lexIsInClass(key, cdt->right->n.classe) ? 
	    BINDEX_TAG_PRESENT : BINDEX_TAG_ABSENT ;
	
	for (j=0; j <arrayMax(a); j++)
	  if (lexIsInClass(key, array(a, j, KEY)))
	    return BINDEX_TAG_PRESENT;
        
	return BINDEX_TAG_ABSENT ;
      }
      
    case FINALPIPE:
      return BINDEX_TAG_PRESENT ;

    case GREP: /* should not happen inside a query */
      return BINDEX_TAG_ABSENT ;

    case WEBQUERY: /* should not happen inside a query */
      return BINDEX_TAG_ABSENT ;

    case PIPEQUERY:
      return condMatch(s,key,objp, cdt->right) ; /* check the rest of the expression */

    default:
      messout("Confusion in condMatch at: %s: %s",
	       qName[cdt->qop], stackText(s, cdt->mark)) ;
      return BINDEX_TAG_ABSENT ;
    }
  return BINDEX_TAG_ABSENT ; /* for compiler happiness */
}

/***********************************************************/
    /* printf since for debugging */
#ifdef DEBUG
static void cdtDump(Stack s, CDT cdt, int level, BOOL isLeft)
{
  int i = level ;
  if (!cdt)
    { printf("\n cdtDumping NULL \n") ;
      return ;
    }
  while(i--) 
    printf("\t") ;
  printf(isLeft ? "< " : "> ") ;
  printf("%s:", qName[cdt->qop]) ;
  if (cdt->qop == TAG && cdt->n.tag)
    printf(" tag = %s ", name(cdt->n.tag) );
  else  if (cdt->qop == CLASS && cdt->n.classe)
    printf(" class = %s ", name(cdt->n.classe)) ;
  if ( cdt->mark)
      printf(" cp = %s ",  stackText(s, cdt->mark)) ;
  printf("\n") ;
  if (cdt->left)
    cdtDump(s,cdt->left, level + 1, TRUE) ;
  if (cdt->right)
    cdtDump(s,cdt->right, level + 1, FALSE) ;
}
#endif

/**********************************************************/

static BOOL condFMatch(Stack s, KEY key, OBJ* objp, CDT cdt)
{
  if (bIndexVersion(-1) && !(*objp) )
    {
      fuzzy = TRUE ;
      switch (condMatch(s,key,objp,cdt))
	{
	case BINDEX_TAG_ABSENT: fuzzy = FALSE ; return FALSE ;
	case BINDEX_TAG_UNCLEAR: fuzzy = FALSE ; break ;
	case BINDEX_TAG_PRESENT: fuzzy = FALSE ; 
	  if ( ! localiseTag) return TRUE ;
	  else break ; /* because we want to position in the object */
	}
    }

  switch (condMatch(s,key,objp,cdt))
    {
    case BINDEX_TAG_ABSENT: return FALSE ;
    case BINDEX_TAG_UNCLEAR: messcrash ("inconsistency in condFMatch") ;
    case BINDEX_TAG_PRESENT: return TRUE ;
    }
  return FALSE ;
}

/**********************************************************/
   /* Works en place */
static void queryFilter (KEYSET ks, CDT cdt, Stack stack) 
  /* NOTE: oldSet is destroyed here */
{    
  register int i = 0, imax, j = 0 ;
  KEY k ;
  OBJ obj = 0 ;   /* created if needed by condMatch */

  if (!ks)
    messcrash ("query filter received a null ks") ;

  subObjMarks = stackReCreate(subObjMarks, 8) ;
  ac_free (subObjMarksHandle)  ;
  subObjMarksHandle = ac_new_handle () ;
  imax = keySetMax(ks) ;
  /* works on a single object at a time */
  for(i=0; !doInterrupt && i<imax; i++)
    { 
      k = keySet(ks,i) ;
      if (i && (!(i%10)) && messIsInterruptCalled())
	{ 
	  doInterrupt = TRUE ;
	  break ;
	}
   
      obj = 0 ;
      if (condFMatch(stack, k, &obj, cdt))
	keySet(ks,j++) = k ;	/* Keeps set sorted as a subset */
      if (obj)
	bsDestroy(obj) ;
    }

  keySetMax (ks) = j ;
}

/**********************************************************/

static KEYSET  querySetOperation (KEYSET oldSet, CDT cdt, Stack stack) 
  /* NOTE: oldSet is destroyed in here or by functions called in here */
{ 
  KEYSET aa = 0, bb = 0, cc = 0;

  switch (cdt->qop)
    {
    case SUBQUERY:
      if (cdt->left)
	{
#ifdef DEBUG
	  int id = oldSet ? oldSet->id : 0 ;
#endif
	  KEYSET newSet =
	    queryExecute(oldSet, cdt->left, stack); /* destroys oldSet */
#ifdef DEBUG
	  if (oldSet && oldSet->magic == ARRAY_MAGIC && oldSet->magic == id && newSet->id != id)
	    messcrash ("querySetOperation(): case SUBQUERY:\n"
		       " oldSet should have been destroyed!");
#endif
	  return newSet ;
	}
      else
	{
	  keySetDestroy (oldSet);
	  return keySetCreate () ;
	}

    case SETDEFINE:
      keySetDestroy (oldSet);
      return keySetCreate () ;

    case PIPEQUERY: 
    case WHERE:
      {
	aa = 0;
	cc = oldSet ;

	if (cdt->left)
	  {
#ifdef DEBUG
	    int id = cc ? cc->id : 0 ;
#endif
 	    aa = queryExecute(cc, cdt->left, stack) ; /* destroys cc */
#ifdef DEBUG
	  if (cc && cc->magic == ARRAY_MAGIC && cc->id == id && aa->id != id)
	    messcrash ("querySetOperation(): case WHERE left:\n"
		       " cc should have been destroyed!");
#endif
	    cc = aa ;
	  }

	if (cdt->right)
	  { 
#ifdef DEBUG
	    int id = cc ? cc->id : 0 ;
#endif
	    aa = queryExecute(cc, cdt->right, stack) ; /* destroys cc */
#ifdef DEBUG
	    if (cc && cc->magic == ARRAY_MAGIC && cc->id == id && aa->id != id) 
	      messcrash ("querySetOperation(): case WHERE right:\n"
		       " cc should have been destroyed!");
#endif
	    cc = aa ; 
	  }
	return cc ;
      }
      break;
      
    default:
      /* evaluate the left and right query and then combine
	 the two sets according to the set-operator */
      {
	KEYSET oldSet1 = oldSet; /* gets destroyed by first call */
	KEYSET oldSet2 = keySetCopy(oldSet); /* gets destroyed by 2nd call */
#ifdef DEBUG
	int id1 = oldSet1 ? oldSet1->id : 0 ;
	int id2 = oldSet2 ? oldSet2->id : 0 ;
#endif

	aa = queryExecute(oldSet1, cdt->left, stack) ; /* destroys oldSet */
	if (cdt->qop == SETELSE &&
	    keySetMax (aa))
	  {
	    keySetDestroy (oldSet2) ;
	    return aa ;
	  }
	bb = queryExecute(oldSet2, cdt->right, stack) ; /* destroys oldSet's copy */

#ifdef DEBUG
	if (oldSet1 && oldSet1->magic == ARRAY_MAGIC && oldSet1->id == id1 && aa->id != id1)
	  messcrash ("querySetOperation(): case default:\n"
		     " oldSet should have been destroyed!");
	if (oldSet2 && oldSet2->magic == ARRAY_MAGIC && oldSet2->id == id2 && bb->id != id2)
	  messcrash ("querySetOperation(): case default:\n"
		     " oldSet's copy should have been destroyed!");
#endif
      
	switch (cdt->qop)
	  {
	  case SETOR:
	    cc = keySetOR(aa, bb) ;
	    break ;  
	  case SETXOR:
	    cc = keySetXOR(aa, bb) ;
	    break ;  
	  case SETAND:
	    cc = keySetAND(aa, bb) ;
	    break ;  
	  case SETMINUS:
	    cc = keySetMINUS(aa, bb) ;
	    break ;  
	  case SETELSE:
	    cc = keySetOR (aa, bb) ;
	    break ;  
	  default: break ;
	  }
	keySetDestroy(aa) ; 
	keySetDestroy(bb) ; 
	return cc ;
      }
    }
}

/**********************************************************/
   /* Works en place */
static KEYSET queryFrom(CDT cdt, Stack stack) 
{    
  messout("Sorry Named Sets are not yet implemented, please complain") ;
  return keySetCreate() ;
}

/***********************************/

static KEYSET queryExpand(KEYSET oldSet, CDT cdt, Stack stack) 
{
  KEYSET s = keySetCreate() ;
  KEY k , k2 ;
  int  i , imax = oldSet ? keySetMax(oldSet) : 0  , j = 0 ;
  OBJ obj = 0 ;
  QOP qop = cdt->right->qop ;
  KEY tag  = qop == TAG ? cdt->right->n.tag : 0 ;
  subObjMarks = stackReCreate(subObjMarks, 8) ;
  ac_free (subObjMarksHandle)  ;
  subObjMarksHandle = ac_new_handle () ;

  for(i=0;i<imax;i++)
    { 
      k = keySet(oldSet,i) ;
      if (i && (!(i%10)) && messIsInterruptCalled())
	{ 
	  doInterrupt = TRUE ;
	  break ;
	}
      
      if (tag)
	{ 
	  if (bIndexFind(k,tag) &&
	      bsIsTagInObj (0, k, tag) && /* open only if tag in model */
	      (obj = bsCreate(k)))
	    {
	      if (bsGetKey(obj,tag,&k2))
		{ 
		  if (class(k2))
		    keySet(s, j++) = k2 ;
		  while(bsGetKey(obj,_bsDown,&k2))
		    if (class(k2))
		      keySet(s, j++) = k2 ;
		}
	      bsDestroy(obj) ;
	    }
	}
      else /* case tag#stuff: locate, then get right column */
	{
	  BOOL old = localiseTag ;
	  obj = 0 ; 
	  localiseTag = TRUE ; /* force localising to tag in condFMatch */
	  
	  if (condFMatch(stack, k, &obj, cdt->right) &&
	      obj &&
	      bsGetKey(obj,_bsRight,&k2))
	    {
	      do
		{  
		  if (class(k2))
		    keySet(s, j++) = k2 ;
		} while (bsGetKey(obj,_bsDown,&k2)) ; /* do-while */
	    }
	  bsDestroy(obj) ;  
	  localiseTag = old ;
	}
    }
  keySetSort (s) ;
  keySetCompress(s) ;
  if (cdt->left)
    queryFilter(s, cdt->left, stack) ;

  keySetDestroy(oldSet);

  return s ;
} /* queryExpand */

/***********************************/
/* OR the content of each keyset obj into the result */
static KEYSET queryExpandKeysets (KEYSET oldSet, CDT cdt, Stack stack) 
{
  KEYSET s = keySetCreate(), ks = 0 ;
  KEY k ;
  int  i , imax = oldSet ? keySetMax(oldSet) : 0  , i1, j = 0 ;

  for(i=0;i<imax;i++)
    { 
      k = keySet(oldSet,i) ;
      if (class (k) != _VKeySet)
	continue ;
      
      if (i && (!(i%10)) && messIsInterruptCalled())
	{ 
	  doInterrupt = TRUE ;
	  break ;
	}

      ks = arrayGet (k, KEY, "k") ;
      if (ks)
	{
	  for (i1 = 0 ; i1 < keySetMax(ks) ; i1++)
	    keySet(s, j++) = keySet(ks, i1) ;
	}
      keySetDestroy (ks) ;
    }

  keySetSort (s) ;
  keySetCompress(s) ;

  keySetDestroy(oldSet);

  return s ;
} /* queryExpandKeysets  */

/***********************************/

static KEYSET queryClass(CDT cdt, Stack stack) 
{
  KEYSET s = keySetCreate() ;
  char *cp ;
  KEY k = 0 ; /* primes lexanext */
  KEY t =  cdt->right->n.classe ; int i = 0 ;
  Array a = pickIsComposite(t);
  
  /* deal with composite classes */
  if (a)
    {
      int j;
      
      for (j=0; j <arrayMax(a); j++)
	{
	  k = 0;
	  while(lexNext(array(a, j, KEY) ,&k))
	    keySet(s,i++) = k ;
	}

      keySetSort (s) ;
      keySetCompress(s) ;
      
      if (cdt->left)
	queryFilter(s, cdt->left, stack) ;
      
      return s ;
    }
  else  if (t == _VText)
    { 
      CDT c = cdt->up ;
      
      if (!c || !c->right || c->right->qop != TEXT)
	{
	  messout ("FIND Text should be followed by a Text, try to \"double quote it\"") ;
	  return s ;
	}
      cp = stackText(stack, c->right->mark) ;
      if (strlen(cp) < 3 )
	{ 
	  messout ("FIND Text should be followed by at least 3 characters") ;
	  return s ;
	}
      i = stackMark(stack) ;
      pushText(stack, "*") ;
      catText(stack, cp) ;
      catText(stack, "*") ;

      keySetDestroy(s) ;
      return 
	queryGrep(0, stackText(stack,i)) ;
    }

  
  while(lexNext(t,&k))
   keySet(s,i++) = k ;

  /* keySetSort(s) ;  keysetCompress not needed */
  if (cdt->left)
    queryFilter(s, cdt->left, stack) ;

  return s ;
}

/**********************************************************/

static KEYSET queryExecute (KEYSET oldSet, CDT cdt, Stack stack)
  /* NOTE: oldSet is destroyed in here or by functions called in here */
{    
  CDT cc = cdt ;

  while(cc->qop == SUBEXP && cc->right)
    cc = cc->right ;

  switch(cc->qop) /* was cdt->qop */
    {
    case FOLLOW_TAG:
      if (cc->right && isLocatorOp (cc->right->qop))
	return queryExpand(oldSet, cc, stack);	/* destroys oldSet */
      break ;

    case EXPAND:
      return queryExpandKeysets (oldSet, cc, stack);	/* destroys oldSet */
      break ;

    case NEIGHBOURS:
      return keySetNeighbours(oldSet); /* destroys oldSet */
      break ;

    case GREP:
      if (cc->right)
	{
	  char *cp ;
	  int i = stackMark(stack) ;
	  cp = stackText(stack,cc->right->mark) ;
	  if (strlen(cp) >= 3)
	    {
	      KEYSET grepKs ;
	      pushText(stack,"*") ;
	      catText(stack,cp) ;
	      catText(stack,"*") ;
	      grepKs = queryGrep (oldSet, stackText(stack,i)) ;
	      keySetDestroy (oldSet);
	      return grepKs ;
	    }
	}
      keySetDestroy (oldSet);
      return keySetCreate() ;
      break ;

    case WEBQUERY:
      if (cc->right && webQueryFunc)
	{
	  char *cp ;
	  cp = stackText(stack,cc->right->mark) ;
	  if (strlen(cp) >= 3)
	    {
	      KEYSET webKs ;
	      webKs = webQueryFunc (oldSet, cp) ;
	      keySetDestroy (oldSet);
	      return webKs ;
	    }
	}
      keySetDestroy (oldSet);
      return keySetCreate() ;
      break ;

    case FIND:
      keySetDestroy (oldSet);

      if (cc->right && cc->right->qop == CLASS)
	{ /* fishy acec patch J&V4D oct 2001 */
	  /* treat the specaile case :: query find Sequence IS \"T28A11\" */
	  CDT cc1 = 0 ;

	  if (cc->up && cc->up->qop == PIPEQUERY &&
	      cc->up->right &&
	      (
	       (
		cc->up->right->qop == IS &&
		cc->up->right->right && cc->up->right->right->qop == TEXT && 
		cc->up->right->right->mark &&  cc->up->right->right->mark < stackMark (stack) &&
		(cc1 = cc->up->right->right)
		) ||
	       (
		cc->up->right->qop == TEXT &&
		cc->up->right->mark &&  cc->up->right->mark < stackMark (stack) &&
		(cc1 = cc->up->right)
		) 
	       )
	      )
	    {
	      KEY dummy = 0 ;
	      KEYSET ks = keySetCreate () ; 
	      char * cp = stackText(stack,cc1->mark) ;
	      if (lexword2key (cp, &dummy,  cc->right->n.classe)) 
		{
		  keySet (ks, 0) = dummy ;
		  return ks ; /* success */
		}
	      while (*cp && *cp != '*' && *cp != '?') cp++ ;
	      if (!*cp) return ks ; /* no wild char, we are done */
		
	      keySetDestroy (ks) ;
	      /* fall thru on general case */
	    }
	  
	  return queryClass(cc, stack) ; 
	}

      return keySetCreate() ;
      break ;

    case FROM:
      keySetDestroy (oldSet);

      if (cc->right && cc->right->qop == NAMEDSET)
	return queryFrom(cc, stack) ; 
      break ;

    case PIPEQUERY: 
    case SUBQUERY: 
    case WHERE:
    case SETAND: 
    case SETOR: 
    case SETXOR: 
    case SETMINUS: 
    case SETDEFINE:
    case SETELSE:
      return querySetOperation(oldSet, cc, stack) ;

    case NAMEDSET:
      keySetDestroy (oldSet);

      messout("NAMEDSET not yet coded, please complain, sorry") ;
      return keySetCreate () ;

    default:
      if (oldSet)
	{
	  queryFilter(oldSet, cc, stack) ;
	  return oldSet ;
	}
    }
  return keySetCreate() ; 
}

/**********************************************************/
  /* Hopefully clear local error message on easy errors */

static BOOL checkParentheses (const char *text)
{
  static Stack bra = 0 ;
  char c, old ;
  const char *cp = text ;
  char *text1 ;
  BOOL tt = FALSE ;

  bra = stackReCreate (bra, 40) ;
  while ((c = *cp++))
    { 
      if (c == '\\')
	{ if (*cp) cp++ ; /* jump protected char */
	  continue ; 
	}
      if (c == '\"') tt = !tt ;
      if (!tt) 
	switch (c)
	  {
	  case '{': push (bra, '}', char) ; break ;
	  case '[': push (bra, ']', char) ; break ;
	  case '(': push (bra, ')', char) ; break ;
	  case '}': case ']': case ')':
	    if (stackEmpty (bra))
	      { 
		text1 = strnew (text, 0) ;
		text1[cp - text] = 0 ;
		messout ("Unbalanced %c in:\n %s", c, text1) ;
		messfree (text1) ;
		return FALSE ;
	      }
	    old = pop (bra, char) ;
	    if (old != c)
	      { 
		text1 = strnew (text, 0) ;
		text1[cp - text] = 0 ;
		messout ("Mismatched %c in:\n %s", c, text1) ;
		messfree (text1) ;
		return FALSE ;
	      }
	    break ;
	  }
    }
  
  if (!stackEmpty (bra))
    { messout ("Missing final %c in:\n %s", pop(bra, char), text) ;
      return FALSE ;
    }
  if (tt)
    { messout ("Missing final \" in:\n %s", text) ;
      return FALSE ;
    }
  return TRUE ;
}

/**********************************************************/

static KEYSET querySingleCommand (KEYSET oldSet, char *text)
     /* NOTE: this destroys oldSet */
{ COND cond = 0 ;   
  KEYSET newSet = 0 ;

  if (condConstruct (text, &cond))
    { 
#ifdef DEBUG
      int id = oldSet ? oldSet->id : 0 ;
#endif
      newSet = queryExecute (oldSet,   /* NOTE: destroys oldSet */
			     cond->cdt, cond->stack);
#ifdef DEBUG
	  if (oldSet && oldSet->magic == ARRAY_MAGIC && oldSet->id == id && newSet->id != id)
	    messcrash ("querySingleCommand():\n"
		       " oldSet should have been destroyed!");
#endif
      condDestroy (cond) ;
    }
  else
    keySetDestroy(oldSet);

  return newSet ;
}

/**********************************************************/
KEYSET oldSetForServer = 0 ;
  /* Can process several command separated by semi columns */
KEYSET queryParametrized (KEYSET oldSet, const char *text, const char* parms)
{ 
  KEYSET newSet = 0, runningSet = 0 ;
  static Stack myText = 0 , runningText = 0 ;
  int m, mine ;
  char *cp ;

  static int subQueryLevel = 0 ;
  if (externalServer)
    {
      /* ignore parms, i think they are never used */
      oldSetForServer = oldSet ;
      cp = strnew (text,0) ; /* stabilize text, it may be a messprintf */
      newSet = externalServer (-1, cp, 0, FALSE) ; 
      messfree(cp) ;
      return newSet ;
      /* -1 means export oldset to server for evaluation */
    }


  /* oldSet belongs to calling routine */
  runningSet = keySetExists(oldSet) ? keySetCopy(oldSet) : 0 ;

  if (!subQueryLevel++)
    { doInterrupt = FALSE ;
      subQ = arrayReCreate(subQ, 8, Array) ;
      myText = stackReCreate(myText, 80) ;
    }
 
  m = stackMark(myText) ;
  if (text)
    pushText(myText, text) ;
  mine = freesettext (stackText(myText, m), parms) ;
  freespecial ("\n\t\"/%\\@") ;    /* No subshells ($ removed) */
  while (!doInterrupt && freecard(mine))
    { 
      if ((cp = freepos()) && *cp)
	{
#ifdef DEBUG
	  int id = runningSet ? runningSet->id : 0 ;
#endif
	  runningText = stackReCreate(runningText, 80) ;
	  pushText(runningText, cp) ;

	  newSet = 
	    querySingleCommand(runningSet, /* destroys runningSet */
			       stackText(runningText,0)) ;

#ifdef DEBUG
	  if (newSet->id != id &&
	      runningSet && runningSet->magic == ARRAY_MAGIC &&
	      runningSet->id == id)
	    messcrash("queryParametrized() - "
		      "runningSet should have been destroyed by "
		      "querySingleCommand()!\n");
#endif

	  runningSet = newSet ;
	}
    }
  while (freecard (mine)) ; /* if interrupted */
  if (!newSet)  /* happens in case text = 0, 
		   or when, the query was incorrect
		   (e.g. a messerror in cdtReorder) */
    newSet = keySetCreate() ;
  if (!--subQueryLevel)
    { m = arrayMax(subQ) ;
      while(m--)
	arrayDestroy(arr(subQ,m,Array)) ;
      arrayDestroy(subQ) ;
    }
  return newSet ;
}

/***************************************************************************/
   /* filter text applied to a single key */
   /* avoids having to create s each time */
KEYSET queryKeyParametrized (KEY key,  const char *text, const char *parms)
{
  static KEYSET s = 0 ;
  if (!s)
    s = arrayCreate(1,KEY) ;
  keySet(s,0) = key ;
  return queryParametrized (s, text, parms) ;
}

/***************************************************************************/
   /* runs on the client, never on server */

KEYSET queryLocalParametrized (KEYSET ks1, const char *text, const char *parms)
{
  KEYSET ks = 0 ;
  KEYSET (*myExternalServer) (KEY key, char *query, char *grep, BOOL getNeighbours) ;

  myExternalServer = externalServer ;
  if (externalServer)
    externalServer (-2,0,0,0) ;   /* disable */
  ks = queryParametrized (ks1, text, parms) ;
  if (myExternalServer)
     myExternalServer (-3,0,0,0) ; /* reenable */
  return ks ;
}

/***************************************************************/
/* export the keySet containing the transitive closure or key->tag */
KEYSET queryKeyTransitiveClosure (KEY key, KEY tag) 
{
  KEYSET ks = keySetCreate (), ks1 = 0 ;
  int nn = 1, n2 = 0 ;
  char *question = strnew (messprintf (" { IS *} SETOR {> %s} ", name (tag)), 0) ;
  keySet (ks, 0) = key ;
  
  while (nn > n2)
    {
      ks1 = query (ks, question) ;
      n2 = nn ;
      nn = keySetMax (ks1) ;
      keySetDestroy (ks) ;
      ks = ks1 ; ks1 = 0 ;
    }
  return ks ;
} /* queryKeyTransitiveClosure */

/**********************************************************************/

        /* Search all the X classes and all visible names */

KEYSET queryGrep (KEYSET oldSet, const char *text)
{ 
  KEY k, k1 ;
  int t, n ;
  OBJ obj ; 
  const char *cp ;
  char *buffer = 0, *cq, cc ;
  KEYSET resultKs = 0, autoKs = 0, namedKs = 0 ;
  
  /* clean up spaces except if \ protected */
  {
    cq = buffer = messalloc (strlen(text) + 1) ;
    cp = text ; 
    while ((cc = *cp++))
      switch (cc)
	{
	case '\\':
	  if (*cp) /* read the char following the \ */
	    *cq++ = *cp++ ;
	  break ;
	case ' ':
	  *cq++ = '*' ;
	  break ;
	default:
	  *cq++ = cc ;
	  break ;
	}
  }


  /* collect is autoKs the autoXrefed classes */
  autoKs = keySetCreate() ; n = 0 ;
  t = 256 ;
  while(!doInterrupt && t--)
    {
      
      if (pickXref(t))  /* pick all cross referenced objects */
	{ for (k = 0 ; lexNext (t,&k) ;)
	    { cp = 0 ; 
	      cp = name(k) ;

	      if (messIsInterruptCalled())
		{ doInterrupt = TRUE ;
		  break ;
		}

	      if (pickMatch (cp, buffer))
		{ if (bIndexFind(k,_Quoted_in) && (obj = bsCreate(k)))
		    {
		      if (bsGetKey (obj, _Quoted_in, &k1)) 
			do keySet(autoKs,n++) = k1 ;
		      while (bsGetKey (obj, _bsDown, &k1)) ;
		      bsDestroy(obj);
		    }
		}
	    }
	}
    }
  if (n > 1)
    {
      keySetSort(autoKs) ;
      keySetCompress(autoKs) ;
    }

  namedKs = keySetCreate () ; n = 0 ;
  if (oldSet) 
    {
      /* collect the oldSetNames fitting the template */
      for (t = 0 ; !doInterrupt && t < keySetMax(oldSet) ; t++)
	{
	  k = keySet(oldSet, t);  /* pick matching objects in old only */
	  cp = 0 ;
	  while (nextName(k, &cp))
	    if (pickMatch (cp, buffer)) 
	      keySet(namedKs, n++) =  k ;  /* pick Self */
	}
    }
  else
    {
      /* collect any name fitting the template */
      int i ;
      t = 256 ;
      while(!doInterrupt && t--)  /* pick matching objects everywhere */
	if (!pickXref(t))
	  for (k = i = 0 ; !doInterrupt && lexNext (t,&k) ; i++)
	    { 

	      if (i && (!(i%10)) && messIsInterruptCalled())
		{ doInterrupt = TRUE ;
		  break ;
		}

	      cp = 0 ;
	      while (nextName(k, &cp))
		if (pickMatch (cp, buffer))
		   keySet(namedKs, n++) =  k ;  /* pick Self */
	    }
    }

  if (n > 1)
    {
      keySetSort(namedKs) ;
      keySetCompress(namedKs) ;
    }

  if (externalServer)   /* merge with objects from the server */
    {
      KEYSET tmp = autoKs ;
      KEYSET externalKs = externalServer (0, 0, buffer,FALSE) ;
      autoKs = keySetOR (tmp, externalKs) ;

      keySetDestroy (externalKs) ;
      keySetDestroy (tmp) ;
    }

  if (oldSet) /* all that needed to catch the ?Text refered by old */
    { 
      KEYSET tmp = autoKs ;

      autoKs = keySetAND (oldSet, tmp) ;
      keySetDestroy (tmp) ;
    }

  resultKs = keySetOR (namedKs, autoKs) ;
  keySetDestroy (autoKs) ;
  keySetDestroy (namedKs) ;
  /* DO NOT  keySetDestroy (oldSet), it belongs to the calling routine  */
  messfree(buffer);

  return resultKs ;
}

/*********************************************************/
   /* checks the syntax of a query */
BOOL  condCheckSyntax (const char *cp) 
{ COND cond ;
  
  if (cp &&
      condConstruct (cp, &cond))
    { condDestroy(cond) ; return TRUE ; }
  else
    return FALSE ;
}

/*********************************************************/
  /* Checks values of filters on object */
unsigned char  queryIsA(OBJ obj, KEY key, Array a, BitSet bb)
{ COND *condp = a && arrayMax(a) ? arrp(a,0,COND) : 0 ;
  int i = a ? arrayMax(a) : 0 ; 
  unsigned char mask = 0 ; 
  unsigned int m1 ;
  OBJ oldObj = obj ;

  if (!i) return mask ;

  subObjMarks = stackReCreate(subObjMarks, 8) ;
  ac_free (subObjMarksHandle)  ;
  subObjMarksHandle = ac_new_handle () ;

  for (i = 0 ; i < arrayMax(a) ; i++)
    { if (obj)
	bsGoto (obj, 0) ; /*  RD 960913 IMPORTANT: reset model to root */
      if (condFMatch((*condp)->stack, key, &obj, (*condp)->cdt))
	bitSet (bb, i) ;
      condp++ ;
    }
  if (!oldObj && obj) /* locally created */
    bsDestroy (obj) ;
  m1 = bigArrayMax(bb) ? bigArr (bb, 0, unsigned int) : 0 ;
  mask = (unsigned char) (m1 & 0xff) ; /* extract lowest byte */
  return mask ;  
}

/*********************************************************/
  /* Checks validity of a query if !key, or if key belongs to subclass */
BOOL queryIsAInit (Stack s, Array a)
{ COND cond ;
  char *cp ;
  int n = 0 ;

  stackCursor(s, 0) ;
  if (stackNextText(s))
    while ((cp = stackNextText(s)))
      { if (cp &&
	    condConstruct (cp, &cond))
	  array(a, n++, COND) = cond ;
      else
        messcrash ("Bad condition %s in subclass definition ", cp) ;
      }
  return TRUE ;
}

/*********************************************************/
  /* Checks validity of a query if !key, or if key belongs to subclass */
BOOL constraintsInit (Stack s, Array a)
{ COND cond ;
  char *cp ;
  int n = 0 ;

  stackCursor(s, 0) ;
  if (stackNextText(s))
    while ((cp = stackNextText(s)))
      { if (cp &&
	    condConstruct (cp, &cond))
	  array(a, n++, COND) = cond ;
      else
        messcrash ("Bad condition %s in constraint definition ", cp) ;
      }
  return TRUE ;
}

/*********************************************************/

BOOL queryFindLocalise (OBJ obj, KEY key, const char *cp)
{
  static COND cond = 0 ;
  BOOL result ;
  BOOL old = localiseTag ;

    /* if !cp, i use again the same condition */
  if (cp)   /* construct cond */
    { if (cond)
	{ condDestroy(cond) ;
	  cond = 0 ;
	}
      if ( !condConstruct (cp, &cond))
	return FALSE ;
    }
  else
    if (!cond)
      messcrash ("First call to queryFind must give a query text") ;

  queryFindCond = cond ;
  if (!key)
    return TRUE ;
  
  localiseTag = TRUE ;
  result = condFMatch(cond->stack, key, (OBJ*)(&obj), cond->cdt) ;
  localiseTag = old ;
  
  return result ;
}

/*********************************************************/

static BOOL queryFind2 (OBJ *objp, KEY key, const char *cp)
{
  static COND cond = 0 ;
  BOOL result ; 
  BOOL old = localiseTag ;

    /* if !cp, i use again the same condition */
  if (cp)   /* construct cond */
    { if (cond)
	{ condDestroy(cond) ;
	  cond = 0 ;
	}
      if (!condConstruct (cp, &cond)) 
	return FALSE ;
    }
  else
    if (!cond)
      messcrash ("First call to queryFind2 must give a query text") ;

  if (!key)
    return FALSE ;   /* not same as above ! */
    
  localiseTag = FALSE ;
  result = 
    condFMatch(cond->stack, key, objp, cond->cdt) ;
  localiseTag = old ;
  
  return result ;
}

BOOL queryFind21 (OBJ *objp, KEY key, const char *cp, const char *tname)
{
  querySetTextName (tname) ;
  return  queryFind2 (objp,key,cp) ;
}

/*********************************************************/

BOOL queryFindLocalise3 (COND cond, OBJ *objp, KEY key)
{ BOOL result ;
  BOOL old = localiseTag ;

  if (!key || !cond)
    return TRUE ;
  localiseTag = TRUE ;
  result = 
    condFMatch(cond->stack, key, objp, cond->cdt) ;
  localiseTag = old ;
  
  return result ;
}

BOOL queryFind3 (COND cond, OBJ *objp, KEY key)
{ BOOL result ;
  BOOL old = localiseTag ;

  subObjMarks = stackReCreate(subObjMarks, 8) ;
  ac_free (subObjMarksHandle)  ;
  subObjMarksHandle = ac_new_handle () ;

  if (!key || !cond)
    return TRUE ;
  localiseTag = FALSE ;
  result = 
    condFMatch(cond->stack, key, objp, cond->cdt) ;
  localiseTag = old ;
  
  return result ;
}

BOOL queryFind31 (COND cond, OBJ *objp, KEY key, const char *tname)
{
  querySetTextName (tname) ;
  return  queryFind3 (cond, objp, key) ;
}

/*********************************************************/
/*
void cdtTest(void)
{ COND cond ;
  char * cp ;
  while(messPrompt("Test queryScan : ", "", "t") )
   if (cp = freeword() &&
      condConstruct(cp,&cond) )
     {
       cdtDump(cond->stack, cond->cdt, 0, FALSE) ;
       condDestroy(cond) ;
     }
}
*/
/*********************************************************/
/*********************************************************/

BOOL queryCalcul (const char *buf, int *ip, float *fp)
{  
  BOOL ok = FALSE ;
  OBJ obj = 0 ;
  COND cond = 0 ;
  CDT cdt = 0 ;

  if (!ip && !fp)
    messcrash ("queryCalcul received 2 void pointers") ;

  if (condCheckSyntax (buf))
    {
      if (condConstruct (buf, &cond))
	{
	  cdt = cond->cdt ;
	  if (cdt->qop == CALCUL && condMatch(cond->stack,0,&obj, cdt))
	    {
	      if (ip)
		{
		  *ip = cdt->x + (cdt->x > 0 ? .5 : -.5) ; /* rounding */
		}
	      else if (fp)
		{
		  *fp = cdt->x ;
		}
	      ok = TRUE ;
	    }
	}
    }
  condDestroy (cond) ;
  return ok ;
}

/*********************************************************/
/*********************************************************/

BOOL queryCheckConstraints (OBJ obj, KEY key, void *v)
{
  COND cond = (COND) v ;
   
  subObjMarks = stackReCreate(subObjMarks, 8) ;
  ac_free (subObjMarksHandle)  ;
  subObjMarksHandle = ac_new_handle () ;

  return 
    condFMatch(cond->stack, key, &obj, cond->cdt) ;
}

/*********************************************************/

void queryClearConstraints (void)
{ COND cond ;
  int i = 256 ;

  while (i--)
    if ((cond = pickList [i].constraints))
      { condDestroy (cond) ;
	pickList [i].constraints = 0 ;
      }
}

/*********************************************************/

BOOL queryConstraintsInit (const char *text, int classe)
{ COND cond ;

  if (condConstruct (text, &cond))
    { pickList [classe].constraints = cond ;
      return TRUE ;
    }
  return FALSE ;
}

/*********************************************************/
/*********************************************************/

 
 
