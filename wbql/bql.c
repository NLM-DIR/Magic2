
#include "acedb.h"
#include "bql.h"
#include "bs.h"
#include "bitset.h"

typedef enum {Zero = 0
	      , SEMICOLUMN, TITLE, ORDER_BY
	      , FROM, SELECT, COMA,  WHERE
	      , WHERE_AVOID, WHERE_IF    /* reserved words */
	      , IN, SET, SET_EQ
	      , OR, OR2, XOR, XOR2, AND, AND2, NOT, NOT2
	      , LIKE2, LIKEREGEXP, LIKE, EQ, NEQ, LEQ, SEQ, LT, ST, LLK, SLK
	      , ISA
	      , MODULO, PLUS, MINUS, MULT, DIVIDE, POWER
	      , DNA, PEPTIDE, PEP, CODING, DATEDIFF
	      , COUNT, MIN, MAX, SUM, AVERAGE, STDDEV
	      , CLNAM, NAM, TIMESTAMP
	      , TTAG, TAG, HASTAG, MERGETAG
	      , RIGHTOF   /* square brackets x[2], meaning right of */
	      , CLASSE   /* , INCLASS */
	      , NUMBER, DOLLAR, VAR, TAGKEY, ACTIVE
	      , OBJECT
	      , ROUND, CURLY, SQUARE  /* () and {} brackets */
	      , DATE, QUOTE, DOUBLEQUOTE 
	      , RAW
	      , EMPTY
	      , LAST

} BQLTYPE ; 

static const char *bqlName[] = { 
  "Zero"
  , ";", "TITLE", "order_by"
  , "from", "select", ",", "where"
  , "__where_avoid__", "__where_if__"   /* reserved words */
  , "in" , ":=", "=" 
  , "||", "OR", "^^", "XOR", "&&" , "AND", "!", "NOT"
  , "like", "=~", "~", "==", "!=", ">=", "<=", ">", "<", ">~", "<~"
  , "ISA"
  , "modulo", "+", "-", "*", "/", "^"
  , "DNA", "PEPTIDE", "PEP", "CODING", "DATEDIFF"
  , "COUNT", "MIN", "MAX", "SUM", "AVERAGE", "STDDEV"
  , ".class", ".name", ".timestamp"
  , ">>", "->", "#", "=>", ":"
  , "class"  /* , ".class" */
  , "number", "$", "var", "key", "@"
  , "object"
  , "()", "{}", "[]"
  , "`", "QUOTE", "DOUBLEQUOTE" 
  , "RAW"
  , "EMPTY"
  , 0
} ;

static const int bqlSide[] = { 
  /* 0   no check, 
   * 1   down optional
   * 2   down mandatory
   * 4   right optional
   * 8   right mandatory
   */
  /* "Zero" */ 0 
  /* ";", "title", "order_by" */  , 5, 2, 2
  /* "from", "select", ",", "where" */, 10, 8, 6, 9
  /* "__where_avoid__", "__where_if__" */, 9, 9
  /* , "in" , ":=", "="  */  , 10, 10, 10
  /* , "||", "OR", "^^", "XOR", "&&" , "AND", "!", "NOT" */  , 10, 10, 10, 10, 10, 10, 8, 8
  /* , ISA */ , 1
  /* , "like", "=~", "~", "==", "!=", ">=", "<=", ">", "<", ">~", "<~" */ , 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10
  /* , "modulo" , "+", "-", "*", "/", "^" */ , 10 , 10, 9, 10, 10, 10
  /* , "DNA", "PEPTIDE", "PEP", "CODING", "DATEDIFF" */  , 8, 8, 8, 8, 8
  /* , "COUNT", "MIN", "MAX", "SUM", "AVERAGE", "STDDEV" */ , 5, 5, 5, 5, 5, 5
  /* , ".class", ".name", ".timestamp" */ , 0, 0, 0
  /* , ">>", "->", "#", "=>", ":" (: means x[3]) */ , 8, 8, 8, 8, 8
  /* , "class"*/ , 1
  /*, "number", "$", "var", "key", "@" */ , 1, 1, 1, 1, 1
  /* , "object" */ , 1
  /*  "()", "{}", "[]" */ , 1, 1, 1
  /* "`", , "QUOTE", "DOUBLEQUOTE" */ , 1, 1, 1
  /* , "RAW" */ , 6
  /* , "EMPTY" */ , 6
  , 0
} ;

typedef struct bqlNode NODE ;
typedef struct bqlDcl  DCL ;

struct  bqlDcl {
  int var ;
  NODE *dclNode ;
} ;

struct bqlNode {
  NODE *up, *down, *right, *dclNode, *parent, *child ;
  BQLTYPE type ;
  int mark, var ;

  double z ;
  char *txt ;
  KEY key ; 
  KEY timeStamp ;
  BSunit uu ;
  KEY uType ;

  BOOL isBool, isNumber, isText, isDate, isVar, isIter, isSet ;

  BOOL avoid, where, where_if ;

  int depth, row, col, nCol ;
  BOOL myObj, myAa, myBsmark ;
  KEY dnaKey ;
  OBJ obj ;
  Array aa ;
  Array dnaD, pep ;
  Stack dnaStack, pepStack ;
  vTXT vtxt ;
  BSMARK bsMark ;
  RegExp br ;
} ;


struct bqlStruct {
  AC_DB db ;
  KEYSET ksIn ;
  KEYSET ksOut ;
  AC_TABLE results ;

  AC_HANDLE h ;
  BOOL debug ;
  BOOL warning ;
  void *magic ;
  BOOL where, where_if ;
  BOOL openObj ;  /* if 0: do not open obj but set mayBe = TRUE  */
  BOOL mayBe, isSorted, doNotSort ;
  BOOL force_zero ;
  Stack s ;
  NODE *node, *from ;
  DICT *dict ;
  vTXT errTxt ;
  Array titles ;
  Array dcls, froms, wheres ;    
  int maxLine ;
  int dummyVar ;
  int semicol ;
  int counting ;
  int minmaxavstd ;
  int minmaxavstdN ;
  int minmaxavstdX ;
  int minmaxavstdX2 ;
  int tableRow ;
  int inCurly ;
  char *order_by ;
  char *title_by ;
} ;

struct bqlIterStruct {
  AC_HANDLE h ;
  void *magic ;
  BQL *bql ;
  KEYSET ksIn, ksCurrent ;
  int currentRow ;
} ;

typedef struct pStruct {
  AC_HANDLE h ;
} PP ;

typedef struct sortStruct {
  BitSet uses, dcls ;
  NODE *node ;
  int level, type ;
} SS ;


static BOOL bqlExpandComa (BQL *bql, NODE *node) ;
static BOOL bqlExpandCount (BQL *bql, NODE *node, NODE *coma) ;
static BOOL bqlExpandMinMax (BQL *bql, NODE *node, NODE *coma) ;
static BOOL bqlExpandWhere (BQL *bql, NODE *node, NODE *coma) ;
static BOOL bqlExpandFrom (BQL *bql, NODE *node, NODE *coma) ;
static BOOL bqlExpandIn (BQL *bql, NODE *node, NODE *coma) ;
static BOOL bqlParseAcedbQuery (BQL *bql, const char *query) ;
static void bqlCleanDcls (NODE *node) ;
static BOOL bqlExpandCurly (BQL *bql, NODE *node) ;

/*************************************************************************************/
/*************************************************************************************/

static int bqlShowNode (BQL *bql, NODE *node, int x, int y) 
{
  int i,dx = 0, dy = 0 ;
  char *cp ;
  NODE *dcl = node->dclNode  ? node->dclNode : node ;

  if (! node->type)
    return 0 ;

  cp = messprintf( "%s ", bqlName[node->type]) ;
  dx += strlen (cp) ;
  fprintf (stderr, "%s", cp) ;
  if (dcl->mark && (dcl->type == CLASSE || dcl->type >= NUMBER))
    {
      cp = stackText (bql->s, dcl->mark) ;
      dx += strlen (cp) ;
      fprintf (stderr, "%s", cp) ;
    }

  if (0) fprintf (stderr, " up:%s", node->up ? bqlName[node->up->type] : "root") ;
  y++ ;

  if (node->down) 
    {
      fprintf (stderr, "      ") ; dx += 6 ;
      dy += bqlShowNode (bql, node->down, x + dx, y) ;
    }
  if (node->right)
    {
      fprintf (stderr, "\n") ;
      for (i = 0 ; i < x ; i++)
	fprintf (stderr, " ") ;
      dy += bqlShowNode (bql, node->right, x, y) ;
    }
  dy++ ;
  if (y == 1)
    fprintf (stderr, " ") ;
  return dy ;
}

static void showBql (BQL *bql)
{
  bqlShowNode (bql, bql->node, 0, 0) ;
  fprintf (stderr, "\n") ;
}

/*************************************************************************************/

static BOOL bqlVar (BQL *bql, NODE *node)
{
  BOOL ok = TRUE ;
  
  if (node->type == RAW)
    {
      char *cp, *cq, *cp0 ;
      int i = 0 ; 
      
      cp0 = node->mark ? stackText (bql->s, node->mark) : 0 ;
      if (cp0 && strlen (cp0) > 0)
	{
	  /* clean up spaces and tabs */
	  for (i = 0, cp = cp0 ; *cp == ' ' || *cp == '\t' ; i++, cp++)
	    {} ;
	  cq = cp ;
	  if (strlen (cp) > 0)
	    for (cq = cp + strlen (cp) - 1 ; cq > cp && (*cq == ' ' || *cq == '\t') ; i++, cq--)
	      {} ;
	  if (i > 0)
	    {
	      char cc = *++cq ; 
	      *cq = 0 ;
	      if (strlen (cp) > 0)
		{
		  node->mark = stackMark (bql->s) ;
		  pushText (bql->s, cp) ;
		}
	      else if (! node->down && ! node->right)
		{
		  node->type = ROUND ;
		  node->mark = 0 ;
		}
	      else
		{
		  node->type = ROUND ;
		  node->mark = 0 ;
		}
	      *cq = cc ;
	    }
	}
      cp0 = node->mark ? stackText (bql->s, node->mark) : 0 ;
      if (! cp0 || strlen (cp0) == 0)
	{
	  node->type = ROUND ;
	  node->mark = 0 ;
	}
    }
  if (node->type == RAW)
    {
      char *cp0 = stackText (bql->s, node->mark) ;
      node->type = VAR ; dictAdd (bql->dict, cp0, &node->var) ; 
    }
  if (node->right)
  bqlVar (bql, node->right) ;
  if (node->down)
    bqlVar (bql, node->down) ;
  
  return ok ;
} /* bqlVar */

/*************************************************************************************/
/* reorder the node->right->right->... into a tree
 * by localizing the least binding operator 
 */
static BOOL bqlReorder (BQL *bql, NODE *node, BOOL *okp)
{
  BOOL ok = FALSE ;
  NODE *right ;
  BQLTYPE type ;
  NODE *best = node ;

  if (! *okp || ! node)
    return FALSE ;
  right = node->right ; type = node->type ;
  while (right)
    {
      if (right->type < type) 
	{ best = right ; type = best->type ; }
      right = right->right ;
    }
  if (node != best)  /* reorder */
    {
      NODE *up = best->up ;
      if (up && up->right == best) up->right = 0 ;
      if (up && up->down == best) up->down = 0 ;
      ok = TRUE ;
      if (! best->down || best->down->type == EMPTY ||
	       ( best->down->type == ROUND && !  best->down->right && !  best->down->down)
	  )
	best->down = node ; 
      else if (! best->right || best->right->type == EMPTY ||
	       ( best->right->type == ROUND && !  best->right->right && !  best->right->down)
	       )
	{
	  best->right = best->down ;	
	  best->down = node ; 
	}
      else if (node->type < EMPTY)
	{
	  *okp = FALSE ;
	  vtxtPrintf (bql->errTxt, "// ... Bad parsing, sorry\n") ;
	  return FALSE ;
	}
	

      up = node->up ;
      if (up && up->right == node) up->right = best ;
      if (up && up->down == node) up->down = best ;
      best->up = up ;
      node->up = best ;
      if (node->right == best) 
	node->right = 0 ;

      if (bql->node == node)
	bql->node = best ;

      return TRUE ; 
    }
  
  if (node->down) while (*okp && bqlReorder (bql, node->down, okp)) {} ;
  if (node->right) while (*okp && bqlReorder (bql, node->right, okp)) {} ;
  return ok ;
} /* bqlReorder */

/*************************************************************************************/
/* recognize the symbols and reserved words
 *  on spaces
 */
static BOOL bqlGetTypes (BQL *bql, NODE *node, BOOL *okp)
{
  BOOL ok = FALSE ;
  int type ;

  while (*okp && node)
    {
       if (node->type != DOLLAR && node->type != QUOTE)
	{
	  if (node->type == RAW && node->mark)
	    {
	      double z = 0 ;
	      char cc ;
	      const char **nam = bqlName ; 
	      char *cp = stackText (bql->s, node->mark) ;
	      for (type = 0 ; *nam ; nam++, type++)
		if (! strcasecmp (cp, *nam) &&
		    (strcasecmp (cp, "TITLE") || ! strcmp (cp, "TITLE")) &&

		    (strcasecmp (cp, "DNA") || ! strcmp (cp, "DNA")) &&
		    (strcasecmp (cp, "PEP") || ! strcmp (cp, "PEP")) &&
		    (strcasecmp (cp, "PEPTIDE") || ! strcmp (cp, "PEPTIDE")) &&
		    (strcasecmp (cp, "CODING") || ! strcmp (cp, "CODING")) &&
		    (strcasecmp (cp, "DATEDIFF") || ! strcmp (cp, "DATEDIFF")) &&

		    (strcasecmp (cp, "OBJECT") || ! strcmp (cp, "OBJECT")) &&

		    /* (strcasecmp (cp, "COUNT") || ! strcmp (cp, "COUNT")) && */
		    (strcasecmp (cp, "MIN") || ! strcmp (cp, "MIN")) &&
		    (strcasecmp (cp, "MAX") || ! strcmp (cp, "MAX")) &&
		    (strcasecmp (cp, "SUM") || ! strcmp (cp, "SUM")) &&
		    (strcasecmp (cp, "AVERAGE") || ! strcmp (cp, "AVERAGE")) &&
		    (strcasecmp (cp, "STDDEV") || ! strcmp (cp, "STDDEV")) 


		    )


		  { 
		    node->type = type ; ok = TRUE ; 
		    if (node->type == LIKE2) node->type = LIKE ;
		    if (node->type == AND2) node->type = AND ;
		    if (node->type == OR2) node->type = OR ;
		    if (node->type == XOR2) node->type = XOR ;
		    if (node->type == NOT2) node->type = NOT ;
		    if (node->type == PEP) node->type = PEPTIDE ;

		    if (0)
		      { /* 2022-03-04, remove this clause bqltest: see the example : tom where ! p->parent && p->papers
			 * it wrongly parses   "where !a && b" into "where ! (a && b)" rather than the correct "where (!a) && b"
			 */
			if (node->type == NOT && node->up && node->up->type == WHERE)
			  node->type = WHERE_AVOID ;
		      }
		    if (node->type == SET_EQ) 
		      {
			NODE *up, *node2 = node ;
			node->type = SET ; 
			while ((up = node2->up) && up->down == node2 && up->type != COMA)
			  { /* if we are inside a where, we rather want EQ */
			    if (up->type == WHERE || up->type == WHERE_IF)
			      {
				node->type = EQ ; /* if we are inside a where, we rather want EQ */
				break ;
			      }
			    node2 = up ;
			  }
			break ;
		      }
		    break ;
		  }
	      if (sscanf (cp, "%lg%c", &z, &cc) == 1)
		{ node->type = NUMBER ; node->z = z ; }
	    }
	  if (node->type == RAW && ! node->mark)
	    {  
	      if (node->right && ! node->down && node->right && node->up)
		{ /* happens with nested parenthesis */
		  if (node == node->up->right)
		    { node->up->right = node->right ; node->right->up = node->up ; }
		  if (node == node->up->down)
		    { node->up->down = node->right ; node->right->up = node->up ; }
		}
	      else
		node->type = EMPTY ;
	    }
	  if (node->down) ok |= bqlGetTypes (bql, node->down, okp) ;
	}
      node = node->right ;
    }
  return ok ;
} /* bqlGetTypes */

/*************************************************************************************/
/* tokenize the text 
 *  on spaces
 */
static BOOL bqlTokenize (BQL *bql, NODE *node, BOOL *okp)
{ 
  BOOL ok = FALSE ;
  NODE *new ;
  char *cp, *cq, *buf = 0 ;
  
  if (! *okp || ! node)
    return FALSE ;
  ok |= bqlTokenize (bql, node->down, okp) ;
  ok |= bqlTokenize (bql, node->right, okp) ;
  
  if (*okp && node->type == RAW)
    {
      NODE *oldRight = node->right ;
      cp = stackText (bql->s, node->mark) ;
      buf = strnew (cp, bql->h) ; /* buf is needed because the stack may be reallocated */
      cp = buf ;
      while (*cp == ' ') cp++ ;
      cq = strchr (cp, ' ') ;
      if (cq && *cp)
	{
	  *cq = 0 ;
	  node->mark = stackMark (bql->s) ;
	  pushText (bql->s, cp) ;
	  new = (NODE *) halloc (sizeof (NODE), bql->h) ;
	  new->type = RAW ;
	  new->mark = stackMark (bql->s) ;
	  cp = cq + 1 ;
	  while (*cp == ' ') cp++ ;
	  pushText (bql->s, *cp ? cp : " ") ;
	  new->right = node->right ;
	  new->up = node ;
	  node->right = new ;
	  new->right = oldRight ;
	  if (oldRight) oldRight->up = new ;
	  node = new ;
	  bqlTokenize (bql, new, okp) ;
	  ok = TRUE ;
	}
    }
  messfree (buf) ;
  return ok ;
} /* bqlTokenize */

/*************************************************************************************/

static BOOL bqlSpaceProtectSymbols (BQL *bql, NODE *node, BOOL *okp)
{
  char *cp0, *cp ;
  int n ;

  if (! *okp || ! node)
    return FALSE ;
  if (*okp) bqlSpaceProtectSymbols (bql, node->right, okp) ;
  if (*okp) bqlSpaceProtectSymbols (bql, node->down, okp) ;

  if (*okp && node->type == RAW && node->mark)
    {
      cp0 = stackText (bql->s, node->mark) ;
      if (cp0 && (n = strlen (cp0)) > 0)
	{
	  char *cq, buffer[3*n + 1] ;
	  
	  n = 0 ;
	  cq = buffer ;
	  for (cp = cp0 ; *cp ; cp++)
	    {
	      switch ((int) *cp)
		{
		case '\'':
		case '\"':
		  messcrash (" bqlSpaceProtectSymbols found char %c", *cp) ;
		  return FALSE ;
		  break ;
		}
	      switch ((int) *cp)
		{
		case '=':
		  n++ ;
		  *cq++ = ' ';
		  *cq++ = *cp ;
		  if (cp[1] == '=' || cp[1] == '~' || cp[1] == '>')  /* eat up following symbol in case of == and =~ and => */
		    *cq++ = *++cp ;
		  *cq++ = ' ' ;
		  break ;
		case '-': /* minus */
		  n++ ;
		  *cq++ = ' ';
		  *cq++ = *cp ;
		  if (cp[1] == '=' || cp[1] == '>')  /* eat up following symbol in case of -= and -> */
		    *cq++ = *++cp ;
		  *cq++ = ' ' ;
		  break ;
		case '^':
		case '|':
		case '&':
		  n++ ;
		  *cq++ = ' ';
		  *cq++ = *cp ;
		  if (cp[1] == cp[0])  /* eat up following symbol in case of &&  || */
		    *cq++ = *++cp ;
		  *cq++ = ' ' ;
		  break ;
		case '>':
		  n++ ;
		  *cq++ = ' ';
		  *cq++ = *cp ;
		  if (cp[1] == '=' || cp[1] == '>' || cp[1] == '~')  /* eat up following symbol in case of >= >~ and >> */
		    *cq++ = *++cp ;
		  *cq++ = ' ' ;
		  break ;
		case '<':
		  n++ ;
		  *cq++ = ' ';
		  *cq++ = *cp ;
		  if (cp[1] == '=' || cp[1] == '~')  /* eat up following = in case of <= <~ */
		    *cq++ = *++cp ;
		  *cq++ = ' ' ;
		  break ;
		case ':':
		case '!':
		case '+':
		case '*':
		case '/':
		case '%':
		  n++ ;
		  *cq++ = ' ';
		  *cq++ = *cp ;
		  if (cp[1] == '=')  /* eat up following = in case of := <= != += *= /= %= */
		    *cq++ = *++cp ;
		  *cq++ = ' ' ;
		  break ;
		case '#':
		case ',':
		  n++ ;
		  *cq++ = ' ';
		  *cq++ = *cp ;
		  *cq++ = ' ' ;
		  break ;
		case '.':
		  if (
		      (!strncasecmp(cp, ".class ",7) || !strcasecmp(cp, ".class")) ||
		      (!strncasecmp(cp, ".name ",6) || !strcasecmp(cp, ".name")) ||
		      (!strncasecmp(cp, ".timestamp ",11) || !strcasecmp(cp, ".timestamp"))
		      )
		    {
		      n++ ;
		      *cq++ = ' ';
		    }
		  *cq++ = *cp ;
		  break ;
		default:
		  *cq++ = *cp ;
		  break ;
		}
	    }
	  *cq++ = 0 ;
	  if (n > 0) /* text is modified, re-register */
	    {
	      node->mark = stackMark (bql->s) ;
	      pushText (bql->s, buffer) ;
	    }	
	}
    }

  return TRUE ;
} /* bqlSpaceProtectSymbols */
  
/*************************************************************************************/

static BOOL bqlDelay (const char *cp0, struct tm *tsp, double *zp)
{
  int x = 0, n = 0, v = 0, unit ;
  const char *cp = cp0 ;
  BOOL wantMonth, wantDay, wantHours, wantMins, wantSecs ;

  if (! tsp || ! cp0) return FALSE ;

  memset (tsp, 0, sizeof(struct tm)) ;
  if ((v = sscanf (cp, "%d%n", &x, &n)) != 1)
    return FALSE ;
  /* check for a unit */
  cp += n ;
  unit = *cp++ ; 
  switch ((int) unit)
    {
    case 0:   /* not a Delay */
      return FALSE ;
    case '-': /* not a Delay */
      return FALSE ;
    case 'Y':
      tsp->tm_year = x ;
      unit = 'Y' ;
      break ;
    case 'M':
      wantMonth = TRUE; 
      tsp->tm_mon = x ;
      if (tsp->tm_mon > 12 || tsp->tm_mon < 1)
	return 0;
       break ;
    case 'D': 
      wantDay = TRUE; 
      tsp->tm_mday = x ;
      if (tsp->tm_mday > 31)
	return 0;
       break ;
    case 'h':
      wantHours = TRUE;  
      tsp->tm_hour = x ;
      if (tsp->tm_hour > 23)
	return 0;
       break ;
    case 'm':
      wantMins = TRUE; 
      tsp->tm_min = x ;
      if (tsp->tm_min > 59)
	return 0;
       break ;
    case 's':
      wantSecs = TRUE;
      tsp->tm_sec = x ;
      if (tsp->tm_sec > 59)
	return 0;
       break ;
    }

  *zp = aceTime(tsp, wantMonth, wantDay, wantHours, wantMins, wantSecs) ;

  return TRUE ;
} /* bqlDelay */

/*************************************************************************************/
/* 
 * on brackets
 * each item is exported in a node type RAW or symbol
 * the parenthetized block is put under and analyzed recursivelly
 */
static BOOL bqlGetQuotes (BQL *bql, NODE *node, BOOL *okp)
{
  BOOL ok = FALSE ;
  char *cp, *cq ;
  int pos ;
  NODE *new, *new2 ;
  char *cp0, *buf = 0 ;
  
  if (! node)
    return FALSE ;
  while (*okp && bqlGetQuotes (bql, node->down, okp)) ;
  while (*okp && bqlGetQuotes (bql, node->right, okp)) ;

  if (! *okp || ! node->mark)
    return FALSE ;
  
  cp = stackText (bql->s, node->mark) ;
  buf = strnew (cp, bql->h) ; /* buf is needed because the stack may be reallocated */
  cp0 = cp = buf ;
  while (*cp)
    {
      cq = cp ; 
      while (*cp == ' ')  /* gobble leading spaces */
	{ cp++ ; pos++ ; }
      if (*cp == '\'' || *cp =='\"' || *cp =='`') /* gobble texts*/
	{
	  char cc = *cp ;
	  ok = TRUE ;
	  *cp = 0 ;
	  node->mark = stackMark (bql->s) ;
	  pushText (bql->s, cp0) ;
	  cq = strchr (cp+1, cc) ;
	  
	  if (! cq)
	    {
	      *okp = FALSE ;
	      vtxtPrintf (bql->errTxt, "// ... Missing %c at the end of %s\n", cc, cp) ;
	      return FALSE ;
	    }
	  *cq = 0 ; 
	  new = (NODE *) halloc (sizeof (NODE), bql->h) ;
	  new->type = DOLLAR ;
	  new->mark = stackMark (bql->s) ; 
	  pushText (bql->s, cp+1) ;
	  new->up = node ;  
	  if (cc == '`') /* a date */  
	    {
	      struct tm ts ;
	      double z1 = 0 ;
	      new->type = DATE ;
	      new->isDate = TRUE ;
	      if (bqlDelay (0, &ts, &z1)) /* use bqlDelay (cp+1, &z1) to activate this code */
		{
		  new->z = z1 ;
		  new->isNumber = TRUE ;
		  new->isDate = FALSE ;
		}
	      else
		{
		  new->z = timeParse (cp + 1) ; 
		  if (! new->z)
		    {
		      *okp = FALSE ;
		      vtxtPrintf (bql->errTxt, "// ... Cannot parse the date `%s`\n", cp + 1) ;
		      return FALSE ;
		    } 
		}
	    } 
	  if (1)
	    {
	      if (cq[1])
		{
		  cp0 = cq + 1 ;
		  new2 = (NODE *) halloc (sizeof (NODE), bql->h) ;
		  new2->type = RAW ;
		  new2->mark = stackMark (bql->s) ;
		  pushText (bql->s, cp0) ;
		  new->right = new2 ;
		  new2->up = new ;
		  new2->right = node->right ;
		}
	      else
		{
		  new->right = node->right ;
		}
	    }
	  node->right = new ; 
	  messfree (buf) ;
	  return TRUE ;
	}
      cp++ ;
    }
  messfree (buf) ;
  return ok ;
} /* bqlGetQuotes */

/*************************************************************************************/
/* 
 * runs after bqlGetBrackets
 * any remaining closing brackets is an error
 */
static BOOL bqlCheckForExtraBrackets (BQL *bql, NODE *node, BOOL *okp)
{
  char *cp, *cq ;
  int pos ;
  BOOL ok = FALSE ;
  char *buf = 0 ;
  
  if (! *okp || ! node)
    return FALSE ;
  if (node->type == DOLLAR)
    return FALSE ;
  while (*okp && bqlCheckForExtraBrackets (bql, node->down, okp)) ;
  while (*okp && bqlCheckForExtraBrackets (bql, node->right, okp)) ;

  if (! node->mark)
    return FALSE ;
  cp = stackText (bql->s, node->mark) ;
  buf = strnew (cp, bql->h) ; /* buf is needed because the stack may be reallocated */
  cp = buf ;
  while (*cp)
    {
      cq = cp ; pos = 0 ;
      while (*cp == ' ')  /* gobble leading spaces */
	{ cp++ ; pos++ ; }
      if (*cp == '\'' || *cp =='\"') /* gobble texts*/
	{
	  cq = strchr (cp+1, *cp) ;
	  if (! cq)
	    {
	      *okp = FALSE ;
	      vtxtPrintf (bql->errTxt, "// ... Missing %c at the end of %s\n", *cp, cp) ;
	      return FALSE ;
	    }
	  cp = cq + 1 ;
	  continue ;
	}
      if (*cp == ')' || *cp == ']' || *cp == '}')
	{ 
	  *okp = FALSE ;
	  cq = buf ;
	  if (strlen (cq) == 1 && node->up && node->up->mark)
	    cq = stackText (bql->s, node->up->mark) ;
	  else if (strlen (cq) == 1 && node->up && node->up->up && node->up->up->mark)
	    cq = stackText (bql->s, node->up->up->mark) ;

	  vtxtPrintf (bql->errTxt, "// ... Unbalanced %c at the end of :\n\t%s\n", *cp, cq) ;
	  return FALSE ;
	}
      cp++ ;
    }
  messfree (buf) ;
  return ok ;
} /* bqlCheckForExtraBrackets */

/*************************************************************************************/
/* 
 * on brackets
 * each item is exported in a node type RAW or symbol
 * the parenthetized block is put under and analyzed recursivelly
 */
static BOOL bqlGetBrackets (BQL *bql, NODE *node, BOOL *okp)
{
  char *cp, *cq ;
  int pos ;
  BOOL ok = FALSE ;
  char *buf = 0 ;
  
  if (! *okp || ! node)
    return FALSE ;
  if (node->type == DOLLAR)
    return FALSE ;
  while (*okp && bqlGetBrackets (bql, node->down, okp)) ;
  while (*okp && bqlGetBrackets (bql, node->right, okp)) ;

  if (! node->mark)
    return FALSE ;
  cp = stackText (bql->s, node->mark) ;
  buf = strnew (cp, bql->h) ; /* buf is needed because the stack may be reallocated */
  cp = buf ;
  while (*cp)
    {
      cq = cp ; pos = 0 ;
      while (*cp == ' ')  /* gobble leading spaces */
	{ cp++ ; pos++ ; }
      if (*cp == '\'' || *cp =='\"') /* gobble texts*/
	{
	  cq = strchr (cp+1, *cp) ;
	  if (! cq)
	    {
	      *okp = FALSE ;
	      vtxtPrintf (bql->errTxt, "// ... Missing %c at the end of %s\n", *cp, cp) ;
	      return FALSE ;
	    }
	  cp = cq + 1 ;
	}
      if (*cp == '(' || *cp == '[' || *cp == '{')
	{
	  BQLTYPE type = ROUND ;
	  char cc1, cc = *cp ;
	  char *cp1 = cp + 1 ;
	  int level = 1 ;
	  cq = cp + 1 ;
	  
	  if (cc == '(') { cc1 = ')' ; type = ROUND ; }
	  if (cc == '[') { cc1 = ']' ; type = RIGHTOF ; }
	  if (cc == '{') { cc1 = '}' ; type = CURLY ; }
	  *cp = 0 ;
	  if (buf[pos])
	    {
	      node->mark = stackMark (bql->s) ;
	      pushText (bql->s, buf) ;
	    }
	  else
	    node->mark = 0 ;
	  while (level)
	    {
	      if (! *cq)
		{
		  *okp = FALSE ;
		  vtxtPrintf (bql->errTxt, "// ... Missing %c at the end of :\n\t%c %s\n", cc1, cc, cp1) ;
		  return FALSE ;
		}
	      else if (*cq == cc)
		level++ ;
	      else if (*cq == cc1)
		level-- ;
	    
	      if (! level)
		{
		  NODE *new, *newR, *down ;

		  *cq = 0 ;
		  new = (NODE *) halloc (sizeof (NODE), bql->h) ;
		  new->type = type ;
		  down = (NODE *) halloc (sizeof (NODE), bql->h) ;
		  down->type = RAW ;
		  down->up = new ;
		  new->down = down ;
		  down->mark = stackMark (bql->s) ;
		  pushText (bql->s, cp1) ;
		  *cq = ')' ;
		  
		  if (cq[1])
		    {
		      newR =(NODE *) halloc (sizeof (NODE), bql->h) ;
		      newR->type = RAW ;
		      newR->up = new ;
		      new->right = newR ;
		      newR->right = node->right ;
		      if (newR->right)
			newR->right->up = newR ;
		      node->right = new ;
		      new->up = node ;
		      newR->mark = stackMark (bql->s) ;
		      pushText (bql->s, cq + 1) ;
		    }
		  else
		    {
		      new->right = node->right ;
		      if (new->right)
			new->right->up = new ;
		      node->right = new ;
		      new->up = node ;
		    }
		  messfree (buf) ;
		  return TRUE ;
		}
	      cq++ ;
	    }
	}
      cp++ ;
    }
  messfree (buf) ;

  return ok ;
} /* bqlGetBrackets */

/**************************************************************/

static KEY goodClass (KEY classe, unsigned char *classFilterp)
{ 
  if (pickIsComposite(classe))
    return classe;
  
  if (!pickIsA (&classe, classFilterp))
    return 0 ;
  
  return classe ;  
} /* goodClass */

/*************************************************************************************/
/* WHERE_AVOID xxx   become  COUNT { 1 where xxx} == 0 */
static BOOL bqlSetAvoid (BQL *bql, NODE *node)
{
  NODE *up = node->up ;

  if (up && node->type == WHERE_AVOID && up->type != WHERE)
    node->type = ROUND ; /* do nothing */
  if (node->type == WHERE_AVOID && ! node->right)
     node->type = NOT ; /* no neeed for complications */
  if (node->type == WHERE_AVOID && (node->right->type == SET || node->right->type >= ISA))
     node->type = NOT ; /* no neeed for complications */
  if (up && node->type == WHERE_AVOID && up->type == WHERE && up->up)
    {
      NODE *stuffNode = up ;
      NODE *eqNode = (NODE *) halloc (sizeof (NODE), bql->h) ;
      NODE *zeroNode = (NODE *) halloc (sizeof (NODE), bql->h) ;
      NODE *countNode = (NODE *) halloc (sizeof (NODE), bql->h) ;
      NODE *curlyNode = (NODE *) halloc (sizeof (NODE), bql->h) ;
      NODE *fromNode = (NODE *) halloc (sizeof (NODE), bql->h) ;
      NODE *selectNode = (NODE *) halloc (sizeof (NODE), bql->h) ;
      NODE *comaNode = (NODE *) halloc (sizeof (NODE), bql->h) ;
      NODE *oneNode = (NODE *) halloc (sizeof (NODE), bql->h) ;
      NODE *whereNode = node ;  /* WHERE_AVOID becomes the where inside COUNT { 1 WHERE ..} == 0 */

      /* new top WHERE node */
      whereNode->type = WHERE ;  
      whereNode->up = up->up ;
      up->up->down = whereNode ;

      if (whereNode->down == 0 &&
	  up->up && up->up->type == COMA &&
	  up->up->up && up->up->up->type == COMA &&
	  up->up->up->down && up->up->up->down->type != WHERE
	  )
	{
	  whereNode->down =  up->up->up->down ;
	  whereNode->down->up = whereNode ;
	  up->up->up->down = whereNode ;
	  whereNode->up = up->up->up ;
	  up->up->up->right = up->up->right ;
	  if (up->up->right)
	    up->up->right->up = up->up->up ;
	}

      if (0 && stuffNode->down)
	{
	  whereNode->down = stuffNode->down ;
	  stuffNode->down->up = whereNode ;
	  stuffNode->down = 0 ;
	}

      /* disconnect the WHERE_AVOID node */
      if (node->right)
	node->right->up = up ;
      if (node == up->right)
	up->right = node->right ;
      if (node == up->down)
	up->down = node->right ;

   
      eqNode->type = EQ ;
      whereNode->right = eqNode ;
      eqNode->up = whereNode ;
      
      countNode->type = COUNT ;
      eqNode->down = countNode ;
      countNode->up = eqNode ;
      eqNode->right = zeroNode ;
      zeroNode->up = eqNode ;

      curlyNode->type = CURLY ;
      countNode->right = curlyNode ;
      curlyNode->up = countNode ;

      fromNode->type = FROM ;
      curlyNode->down = fromNode ;
      fromNode->up = curlyNode ;

      comaNode->type = COMA ;
      fromNode->right = comaNode ;
      comaNode->up = fromNode ;

      comaNode->down = stuffNode ;
      stuffNode->up = comaNode ;

      selectNode->type = SELECT ;
      fromNode->down = selectNode ;
      selectNode->up = fromNode ;
      selectNode->right = oneNode ;
      oneNode->up = selectNode ;

      oneNode->type = NUMBER ;
      oneNode->isNumber = TRUE ;
      oneNode->uu.i = 0 ;
      oneNode->uType = ac_type_int ;
      oneNode->z = 1 ;
  
      zeroNode->type = NUMBER ;
      zeroNode->isNumber = TRUE ;
      zeroNode->uu.i = 0 ;
      zeroNode->uType = ac_type_int ;
      zeroNode->z = 0 ;
   }

  if (node->right)
    bqlSetAvoid (bql, node->right) ;
  if (node->down)
    bqlSetAvoid (bql, node->down) ;
  
  return TRUE ;
}  /* bqlSetAvoid */

/*************************************************************************************/
/* COUNT { xxx } becomes COUNT {select 1 where xxx } */
static BOOL bqlSetCurly (BQL *bql, NODE *node)
{
  NODE *up = node->up ;

  if (0 &&
      up && node->type == CURLY && node->down && node->down->type != WHERE)
    {
      NODE *stuffNode = node->down ;
      NODE *fromNode = (NODE *) halloc (sizeof (NODE), bql->h) ;
      NODE *selectNode = (NODE *) halloc (sizeof (NODE), bql->h) ;
      NODE *oneNode = (NODE *) halloc (sizeof (NODE), bql->h) ;
      NODE *comaNode = (NODE *) halloc (sizeof (NODE), bql->h) ;
      NODE *whereNode = (NODE *) halloc (sizeof (NODE), bql->h) ;

      fromNode->type = FROM ;
      fromNode->up = node ;
      node->down = fromNode ;

      selectNode->type = SELECT ;
      fromNode->down = selectNode ;
      selectNode->up = fromNode ;
      selectNode->right = oneNode ;
      oneNode->up = selectNode ;

      comaNode->type = COMA ;
      comaNode->up = fromNode ;
      fromNode->right = comaNode ;

      whereNode->type = WHERE ;  
      whereNode->up = comaNode ;
      comaNode->down = whereNode ;
      whereNode->right = stuffNode ;
      stuffNode->up = whereNode ;

      oneNode->type = NUMBER ;
      oneNode->isNumber = TRUE ;
      oneNode->uu.i = 0 ;
      oneNode->uType = ac_type_int ;
      oneNode->z = 1 ;
     }

  if (node->right)
    bqlSetCurly (bql, node->right) ;
  if (node->down)
    bqlSetCurly (bql, node->down) ;
  
  return TRUE ;
}  /* bqlSetCurly */

/*************************************************************************************/

static BOOL bqlCleanParenthesis (BQL *bql, NODE *node)
{
  BOOL ok = FALSE ;
  NODE *up = node->up ;

  /* clean out empty cells */
  if (node->type == Zero)
    {
      if (! node->down && ! node->right)
	{
	  if (node == up->down) up->down = 0 ;
	  else if (node == up->right) up->right = 0 ;
	  return TRUE ;
	}
      else
	node->type = ROUND ;
    }
  
  if (node->type == VAR 
      && node->up 
      && (node->up->type == IN || node->up->type == SELECT)
      && node->mark
      )
    {
      char *cp = stackText (bql->s, node->mark) ;
      if (*cp++ == '?')
	{
	  char cc, *cq = cp ;
					
	  KEY classe = 0 ;
	  while (*++cq && *cq != ' ')
	    {}
	  cc = *cq ; *cq = 0 ;
	  if (*cp && lexword2key(cp, &classe, _VClass))  
	    {
	      unsigned char cc = 0 ;

	      node->type = CLASSE ;
	      node->var = classe ;
	      node->key = goodClass (classe, &cc) ;
	      node->uType = cc ;
	    }
	  else if (! *cp)
	    {
	  vtxtPrintf (bql->errTxt, "// ... SYNTAX ERROR no pattern visible to the right of IN CLASS ~ operator") ;
	  return FALSE ;
	    }
	  else
	    {
	      bql->warning = TRUE ;
	      vtxtPrintf (bql->errTxt, "// ... WARNING: select class %s : unrecognized class"
			  , cq) ;
	    }
	  *cq = cc ;
	}
    }
 
  /* simplify parenthesis */
  if (node == bql->node && 
      node->type == ROUND)
    {
      if (! node->right)
	{
	  node = bql->node = node->down ;
	  if (node->down) node->down->up = 0 ;
	  return TRUE ;
	}
      if (! node->down)
	{
	  bql->node = node->right ;
	  if (node->right) node->right->up = 0 ;
	  return TRUE ;
	}
    }

  if (up && node->type == ROUND)
    {
      if (! node->right)
	{
	  if (up)
	    {
	      if (node == up->down) up->down = node->down ;
	      else if (node == up->right) up->right = node->down ;
	    }
	  node->type = 0 ; node->mark = 0 ; 
	  if (node->down) node->down->up = up ;
	  ok = TRUE ;
	}
      else if (! node->down)
	{
	  if (up)
	    {
	      if (node == up->down) up->down = node->right ;
	      else if (node == up->right) up->right = node->right ;
	    }
	  node->type = 0 ; node->mark = 0 ; 
	  if (node->right) node->right->up = up ;
	  node = up ? up : node->right ; /* backtrack */
	  ok = TRUE ;
	}
      else
	{
	  if (node->right->right && ! node->right->down)
	    { 
	      node->right->down = node->down ;
	      node->down->up = node->right ;
	      if (node == up->down) up->down = node->right ;
	      else if (node == up->right) up->right = node->right ;
	      ok = TRUE ;
	    }
	}
    }
  if (node->right && node->type != QUOTE)
    ok &= bqlCleanParenthesis (bql, node->right) ;
  if (node->down && node->type != QUOTE)
    ok &= bqlCleanParenthesis(bql, node->down) ;

  switch (node->type)
    {
    case ROUND:
    case CURLY:
    case DOLLAR:
    case QUOTE:
      if ( node->right && node->right->down == 0 && node == node->up->right)
	{
	  NODE *right = node->right ;
	  node->right = 0 ;
	  node->up->right = right ;
	  node->up = right ;
	  right->down = node ;      
	}
      break ;
    case COMA:
      switch (up->type)
	{
	case FROM:
	case SELECT:
	case COMA:
	  break ;
	default:
	  if (node == up->right && ! up->up->right)
	    {
	      NODE *right = node->right ;
	      up->up->right = right ;
	      if (right) right->up = up->up ;
	      up->right = node->down ;
	      if (node->down)
		node->down->up = up ;
	    }
	  break ;
	}
      break ;
    default:
      break ;
    }
    

  return ok ;
} /* bqlCleanParenthesis */
  
/*************************************************************************************/

static BOOL bqlCheckVariableDeclarations (BQL *bql, NODE *node, int pass)
{
  BOOL ok = TRUE ;
  
  if (! node)
    {
      NODE *fromNode = bql->node ;

      if (pass == 0)
	bqlCleanDcls (bql->node) ;

      if (fromNode)  /* check that the variables are declared in the proper order */
	{
	  if (fromNode->right && ! bqlCheckVariableDeclarations (bql, fromNode->right, pass))
	    {
	      vtxtPrintf (bql->errTxt, "// ... BQL ERROR 12: Each variable in the from clause must be declared and initialised before it is used, maybe you forgot a coma between definitions\n") ;
	      bqlShowNode (bql, fromNode->right, 0, 0) ;
	      return FALSE ;      
	    }
	  if (fromNode->down && ! bqlCheckVariableDeclarations (bql, fromNode->down, pass))
	    {
	      vtxtPrintf (bql->errTxt, "// ... BQL ERROR 11: The variable in the select clause use variables not declared in the from clause\n") ;
	      bqlShowNode (bql, fromNode->down, 0, 1) ;
	      return FALSE ;      
	    }
	}
      return TRUE ;
    }
  
  if (node->type == CLNAM || node->type == NAM || node->type == TIMESTAMP)
    {
	  node->dclNode = node ;
    }

  if ((node->type == QUOTE || node->type == DOLLAR) && node->up && node->up->type == CLASSE)
    {
      KEY classe ;
      const char *cp0 = node->mark ? stackText (bql->s, node->mark) : 0 ;
      char *cq = ac_unprotect (cp0, bql->h) ;
      if (cq && strlen (cq) > 0 &&
	  lexword2key(cq, &classe, _VClass))  
	{
	  unsigned char cc = 0 ;
	  
	  node->up->var = classe ; 
	  node->up->key = goodClass (classe, &cc) ;
	  node->up->uType = cc ;
	  node->up->mark = node->mark ;
	  node->type = 0;
	}
      else if (! *cq)
	{
	  vtxtPrintf (bql->errTxt, "// ... SYNTAX ERROR no pattern visible to the right of IN CLASS - operator") ;
	  return FALSE ;
	}
      else
	{
	  bql->warning = TRUE ;
	  vtxtPrintf (bql->errTxt, "// ... WARNING: select CLASS %s : unrecognized class"
		      , cq) ;
	}
    }
  if (node->type == LIKEREGEXP)
    {
      NODE *right = node->right ;
      const char *cp0 = right && right->mark ? stackText (bql->s, right->mark) : 0 ;
      if (! cp0 ||  !*cp0)
	{
	  vtxtPrintf (bql->errTxt, "// ... no pattern visible to the right of \'like\' or \'=~\' string matching request.\nPossibly your regexp involves [] or other symbols as in\n      x =~ a[bc].*\nand you should quote it or double quote it as in\n      x =~ \'a[bc].*\'\n") ;
	  return FALSE ;
	}
      node->br = regExpCreate (cp0, bql->h) ; 
      if (! node->br)
	{
	  vtxtPrintf (bql->errTxt, "// ... bad no pattern found to the right of a LIKEREGEXP =~ operator, expecting a Unix regular expression, found %s", cp0) ;
	  return FALSE ;
	}
    }
  if (node->type == LIKE)
    {
      NODE *right = node->right ;
      const char *cp0 = right && right->mark ? stackText (bql->s, right->mark) : 0 ;
      if (! cp0 ||  !*cp0)
	{
	  vtxtPrintf (bql->errTxt, "// ... no pattern visible to the right of a LIKE ~ operator") ;
	  return FALSE ;
	}
    }

  if (node->type == ISA)
    {
      NODE *right = node->right ;
      char *cp0 = right && right->mark ? stackText (bql->s, right->mark) : 0 ;
      if (! cp0 ||  !*cp0)
	{
	  vtxtPrintf (bql->errTxt, "// ... SYNTAX ERROR no pattern visible to the right of ISA ~ operator") ;
	  return FALSE ;
	}
      if (1)
	{
	  char cc, *cq = cp0 ;
					
	  KEY classe = 0 ;
	  while (*++cq && *cq != ' ')
	    {}
	  cc = *cq ; *cq = 0 ;
	  if (*cp0 && lexword2key(cp0 , &classe, _VClass))  
	    {
	      unsigned char cc = 0 ;

	      node->var = classe ;
	      node->key = goodClass (classe, &cc) ;
	      node->uType = cc ;
	    }
	  else
	    {
	      bql->warning = TRUE ;
	      vtxtPrintf (bql->errTxt, "// ... WARNING: ISA %s : unrecognized class"
			  , cp0) ;
	    }
	  *cq = cc ;
	}

    }

  if (node->type == VAR)
    {
      NODE *up = node->up ;
      DCL *dcl ;

      if (! node->var)  
	{
	  char *cp0 = node->mark ? stackText (bql->s, node->mark) : 0 ;
	  if (cp0 && strlen (cp0) > 0)
	    dictAdd (bql->dict, cp0, &node->var) ;
	}

      dcl = arrayp (bql->dcls, node->var, DCL) ;

      if (node == up->down && (up->type == SET || up->type == IN))
	{
	  /* should not already be declared, declare it */
	  if (dcl->dclNode && pass == 0)
	    {
	      vtxtPrintf (bql->errTxt, "// ... variable %s is redeclared\n", stackText (bql->s,node->mark)) ;
	      return FALSE ;
	    }
	  if (! dcl->dclNode)
	    {
	      dcl->dclNode = node ; 
	      dcl->var = node->var ;
	    }
	  node->dclNode = dcl->dclNode ;
	}
      else if (up->type == CLASSE)
	{
	  KEY classe ;
	  const char *cp0 = node->mark ? stackText (bql->s, node->mark) : 0 ;
	  char *cq = ac_unprotect (cp0, bql->h) ;
	  if (cq && strlen (cq) > 0 &&
	      lexword2key(cq, &classe, _VClass))  
	    {
	      unsigned char cc = 0 ;
	      up->var = classe ; 
	      up->key = goodClass (classe, &cc) ;
	      up->uType = cc ;
	  
	      up->mark = node->mark ;
	      node->type = 0 ;
	    }
	  else
	    {
	      vtxtPrintf (bql->errTxt, "// ... class %s unknown\n", cp0 ? cp0 : "NULL") ;
	      return FALSE ;
	    }
	}
      else if ( 
	       ( node == up->right && (up->type == TAG || up->type == TTAG || up->type == HASTAG || up->type == MERGETAG)) ||
	       ( node == up->down && up->up && (up->type == TAG || up->type == TTAG || up->type == HASTAG || up->type == MERGETAG) && (up->up->type == TAG || up->up->type == TTAG || up->up->type == HASTAG || up->up->type == MERGETAG || up->up->type == CLASSE)) ||
	       ( node == up->down && (up->type == RIGHTOF && (up->up->type == TAG || up->type == TTAG || up->up->type == MERGETAG) ))
		)
	{
	  KEY tag ;
	  const char *cp0 = node->mark ? stackText (bql->s, node->mark) : 0 ;
	  if (cp0 && strlen (cp0) > 0 &&
	      lexword2key(cp0, &tag, 0)
	      )
	    { node->key = tag ; node->type = TAGKEY ; }
	  else
	    {
	      bql->warning = TRUE ;
	      vtxtPrintf (bql->errTxt, "// ... WARNING: tag %s unknown\n", cp0 ? cp0 : "NULL") ; 
	      node->key = 0 ; node->type = TAGKEY ;
	    }
	}
      else if ( 
	       ( arrayMax(bql->dcls) < node->var || ! arrayp (bql->dcls, node->var, DCL)->dclNode) &&
	       node == up->right && 
	       node->type == VAR && ! node->down && ! node->right &&
	       (up->type ==  EQ || up->type == NEQ || up->type == LEQ || up->type == SEQ || up->type == LT  || up->type == ST || up->type == LIKE || up->type == LIKEREGEXP || up->type == ISA)
		)
	node->type = DOLLAR ;
      else
	{
	  /* should already be declared */
	  dcl = arrayp (bql->dcls, node->var, DCL) ;
	  if (! dcl || ! dcl->dclNode)
	    {
	      vtxtPrintf (bql->errTxt, "// ... BQL ERROR 1: variable %s is used  before it is declared\n"
			  , stackText (bql->s,node->mark)) ;
	      return FALSE ;
	    }
	  node->dclNode = dcl->dclNode ;
	}
    }

  if (node->right && node->type == FROM)  /* usueful in from clause embedded inside a {} : COUT { select x from ... } */
    ok &= bqlCheckVariableDeclarations (bql, node->right, pass) ;
  if (node->down && node->type != QUOTE)
    ok &= bqlCheckVariableDeclarations (bql, node->down, pass) ;
  if (node->right  && node->type != FROM && node->type != QUOTE)
    ok &= bqlCheckVariableDeclarations (bql, node->right, pass) ;

  return ok ;
} /*  bqlCheckVariableDeclarations */

/*************************************************************************************/
/* translate variable names into column numbers
 *  x+y-z   -> 1+2-3
 */
static BOOL bqlCheckOrderBy (BQL *bql)
{
  char cc, *cp, *cp0 = bql->order_by, *cr, *cs ;
  BOOL ok = TRUE ;
  int var = 0 ;
  char *bb, *buf = 0 ;

  if (! cp0)
    return TRUE ;
  
  cp = cp0 + strlen (cp0) - 1 ;
  while (cp > cp0 && *cp == ' ') *cp-- = 0 ;
  for (cp = cp0 ; ok && *cp ; cp++)
    switch (*cp)
      {
      case '+':
      case '-':
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
	break ;
      default:
	ok = FALSE ;
	break ; 
      }
  if (ok) /* we are already in column number format */
    return TRUE ;
  ok = TRUE ;
  bb = buf = strnew (bql->order_by, bql->h) ;
  for (cp = cp0 ; ok && *cp ; cp++)
    switch (*cp)
      {
      case '+':
      case '-':
	*bb++ = *cp ;
	break ;
      default:
	cr = cp ;
	while (*cr && *cr != '+' && *cr != '-')
	  cr++ ;
	cc = *cr ;
	*cr = 0 ;
	var = 0 ;
	dictFind (bql->dict, cp, &var) ;
	*cr = cc ;
	if (! var)
	  {
	    ok = FALSE ;
	    break ;
	  } 
	cs = hprintf (bql->h, "%d", var) ;
	while (*cs) *bb++ = *cs++ ;
	cp = cr - 1 ; *cr = cc ;
      }
  if (ok) /* all named variables are known */
    {
      *bb++ = 0 ;
      bql->order_by = buf ;
      return TRUE ;
    }

  vtxtPrintf (bql->errTxt, "// BQL ERROR 100\n// Sorry, bql cannot parse :\n//                  order_by %s\n", bql->order_by) ;
  vtxtPrintf (bql->errTxt, "// ...  The expected syntax is:\n// -either          order_by +3+4-2+5\n// -or              order_by -x+y-z+t\n// where the +- signs mean ascending or descending order (the leading + is optional)\n// and apply to exported variables refered by their names or their column numbers\n") ;
  return FALSE ;
} /* bqlCheckOrderBy */

/*************************************************************************************/
/* translate variable names into column numbers
 *  x+y-z   -> 1+2-3
 */
static BOOL bqlCheckTitleBy (BQL *bql)
{
  char *cp, *cp0 = bql->title_by, *cr, *cs, *ct ;
  int var = 0 ;

  if (! cp0)
    return TRUE ;

  cp0 = strnew (cp0, bql->h) ;
  cp = cp0 + strlen (cp0) - 1 ;
  while (cp > cp0 && *cp == ' ') *cp-- = 0 ;
  
  while (cp0)
    {
      /* split on , */
      cp = cp0 ;
      cr = strchr (cp, ',') ;
      if (cr) *cr = 0 ;
    
      /* split on : */
      cs = strchr (cp, ':') ;
      var = 0 ;
      if (cs)
	{
	  *cs = 0 ;
	  while (*cp == ' ') 
	    cp++ ;
	  ct = cs ;
	  while (ct > cp && *ct == ' ') 
	    ct-- ;
	  /* recognize the numbers of the variable names before the : */
	  if (sscanf (cp, "%d", &var) == 1)
	    ;
	  else
	    dictFind (bql->dict, cp, &var) ;

	  if (var && var > 0 && var <= dictMax (bql->dict))
	    {
	      /* store the tite in the titles array */
	      array (bql->titles, var - 1, const char*) = strnew (cs + 1, bql->h) ;
	    }
	  else
	    {
	      vtxtPrintf (bql->errTxt, "// BQL ERROR 110\n// Sorry, the variable name  or column number %s is not recognized in TITLE %s\n"
			  , cs, bql->title_by) ;
	      return FALSE ;
	    } 
	}
      cp0 = cr ? cr + 1 : 0 ;
    }

  return TRUE ;
} /* bqlCheckTitleBy */

/*************************************************************************************/
/*************************************************************************************/
 
static void bqlCleanDcls (NODE *node)
{
  node->dclNode = 0 ;
  if (node->right)
    bqlCleanDcls (node->right) ;
  if (node->down)
    bqlCleanDcls (node->down) ;
  return ;
} /* bqlCleanDcls */

/*************************************************************************************/
/*************************************************************************************/

static BOOL bqlSplitOnSemicolumn (BQL *bql, NODE *node0)
{
  BOOL ok = FALSE ;
  char *cp = 0, *cr = 0 ;

  cp = stackText (bql->s, node0->mark) ;
  if (! node0 ||
      node0->type == QUOTE ||
      node0->type == DOLLAR
      )
    goto done ;

  cr = cp ?  strchr (cp, ';') : 0 ;
 again:
  if (cr)
    {
      /* check if we are inside a ' */
      if (1)
	{
	  char *cs = cp ;
	  int i = 0 ;
	  for (;cs < cr ; cs++)
	    if (*cs == '\'' && cs[-1] != '\\')
	      i++ ;
	  if (i % 2)
	    {
	      cr = strchr (cr + 1, ';') ;
	      goto again ;
	    }
	}
      /* check if we are inside a " */
      if (1)
	{
	  char *cs = cp ;
	  int i = 0 ;
	  for (;cs < cr ; cs++)
	    if (*cs == '"' && cs[-1] != '\\')
	      i++ ;
	  if (i % 2)
	    {
	      cr = strchr (cr + 1, ';') ;
	      goto again ;
	    }
	}
      /* check if we are inside a ( */
      if (1)
	{
	  char *cs = cp ;
	  int i = 0, j = 0, k = 0 ;
	  for (;cs < cr ; cs++)
	    switch ((int) *cs)
	      {
	      case '(': i++ ; break ; 
	      case '[': j++ ; break ; 
	      case '{': k++ ; break ; 
	      case ')': i-- ; break ; 
	      case ']': j-- ; break ; 
	      case '}': k-- ; break ;
	      } 
	  if (i != 0 || j != 0 || k != 0)
	    {
	      cr = strchr (cr + 1, ';') ;
	      goto again ;
	    }
	}
      *cr++ = 0 ;
      while (cp && *cp == ' ') cp++ ;
      while (cr && *cr == ' ') cr++ ;

      if (*cp || *cr)
	{
	  node0->mark = stackMark (bql->s) ;
	  pushText (bql->s, ";") ;
	  node0->type = SEMICOLUMN ;
	}

      if (cp && *cp)
	{
	  NODE *down = (NODE *) halloc (sizeof (NODE), bql->h) ;
	  node0->down = down ;
	  down->up = node0 ;
	  down->type = RAW ;
	  down->mark = stackMark (bql->s) ;
	  pushText (bql->s, cp) ;
	  
	}

      if (cr && *cr)
	{
	  NODE *right = (NODE *) halloc (sizeof (NODE), bql->h) ;
	  right->up = node0 ;
	  right->type = RAW ;
	  node0->right = right ;
	  right->mark = stackMark (bql->s) ;
	  pushText (bql->s, cr) ;
	}
      ok = TRUE ;
    }

 done:
  if (node0->down)
    ok |= bqlSplitOnSemicolumn (bql, node0->down) ;
  if (node0->right)
    ok |= bqlSplitOnSemicolumn (bql, node0->right) ;
  return ok ;
} /* bqlSplitOnSemicolumn */

/*************************************************************************************/

static BOOL bqlFromActive (BQL *bql, NODE *node, NODE *var)
{
  BOOL ok = FALSE ;

  if (node->type == CURLY)
    return ok ;

  if (node->type == ACTIVE)
    {
      ok = TRUE ;
      node->type = VAR ;
      node->var = var->var ;
      node->dclNode = var ;
      node->mark = var->mark ;
    }

  if (node->down)  ok |= bqlFromActive (bql, node->down, var) ;
  if (node->right) ok |= bqlFromActive (bql, node->right, var) ;

  return ok ;
} /* bqlFromActive */

/*************************************************************************************/
/* replace all non atomic selected block by a new variable
 * select XXX from YYY  -> select -x from YYY, _x in XXX
 * i.e.  cases with implicit select from where
 */
static BOOL bqlDoSetComas (BQL *bql, NODE *node0)
{
  if (node0->type == SELECT)
    {
      NODE *select = node0 ;
      NODE *from = 0, *fromComa, *sComa ;
      NODE *new ;

      if (! select->up || select->up->type != FROM)
	{
	  new = (NODE *) halloc (sizeof (NODE), bql->h) ;
	  new->type = FROM ;
	  new->down = select ;
	  new->up = select->up ;
	  if (new->up && new->up->down == select)
	    new->up->down = new ;
	  if (new->up && new->up->right == select)
	    new->up->right = new ;
	  select->up = new ;
	}
      if (select->up &&
	  select->up->type == FROM
	  )
	from = select->up ;
      if (from)
	{
	  fromComa = from ;
	  while (fromComa->right && fromComa->right->type == COMA)
	    fromComa = fromComa->right ;

	  new = (NODE *) halloc (sizeof (NODE), bql->h) ;
	  new->type = COMA ;
	  new->mark = stackMark (bql->s) ;
	  pushText (bql->s, ",") ;
	  new->up = fromComa ;
	  new->down = fromComa->right ;
	  if (new->down) new->down->up = new ;
	  fromComa->right = new ;
	  fromComa = fromComa->right ;
	}

      /* Force select to be followed by a succession of COMA */
      sComa = select  ;
      while (sComa->right && sComa->right->type == COMA)
	sComa = sComa->right ;
      if (sComa->right) /* not a coma */
	{
	  new = (NODE *) halloc (sizeof (NODE), bql->h) ;
	  new->type = COMA ;
	  new->mark = stackMark (bql->s) ;
	  pushText (bql->s, ",") ;
	  new->up = sComa ;
	  new->down = sComa->right ;
	  if (new->down) new->down->up = new ;
	  sComa->right = new ;
	}
      sComa = select  ;
    }

  if (node0->down)
    bqlDoSetComas (bql, node0->down) ;
  if (node0->right)
    bqlDoSetComas (bql, node0->right) ;

  return TRUE ;
} /*  bqlDoSetComas */

/*************************************************************************************/
/*************************************************************************************/

static int ssLevelOrder (const void *va, const void *vb)
{
  const SS *sa = (const SS*)va ;
  const SS *sb = (const SS*)vb ;
  int nn = sa->level - sb->level ; 
  return nn ;
} /* ssLevelOrder */

/*************************************************************************************/

static BOOL bqlSetDcls (BitSet bb, NODE *node)
{
  if (node->type == VAR && node == node->dclNode)
    bitSet (bb, node->var) ;

  if (node->down)
    bqlSetDcls (bb, node->down) ;
  if (node->right)
    bqlSetDcls (bb, node->right) ;
  return TRUE ;
} /* bqlSetDcls */

/*************************************************************************************/

static BOOL bqlSetUses (BitSet bb, NODE *node)
{
  if (node->type == VAR)
    bitSet (bb, node->var) ;

  if (node->down)
    bqlSetUses (bb, node->down) ;
  if (node->right)
    bqlSetUses (bb, node->right) ;
  return TRUE ;
} /* bqlSetFroms */

/*************************************************************************************/
/* atomize the query (between ;, should also be between {} into its components */

static BOOL bqlAtomize (BQL *bql, NODE * node, Array aa, BQLTYPE type)
{
  BOOL ok = FALSE ;

  if (node->type == CURLY) return FALSE ;
  if (node->type == type && node->down && (type == WHERE || type == WHERE_IF))
    {
      NODE *up = node->up ;
      NODE *down = node->down ;
      NODE *right = node->right ;
      SS *ss = arrayp (aa, arrayMax(aa), SS) ;

      ok = TRUE ;
      ss->node = right ; ss->node->up = 0 ; ss->type = type ;
      *node = *down ;
      node->up = up ; ss->node->up = 0 ;
    }

  else if (node->type == type)
    {
      SS *ss ;
      NODE *sComa = node ;
      
      while (sComa)
	{
	  sComa = sComa->right ;
	  if (sComa && sComa->type == COMA && sComa->down)
	    {
	      ok = TRUE ;
	      ss = arrayp (aa, arrayMax(aa), SS) ;
	      ss->node = sComa->down ;
	      ss->node->up = sComa->down = 0 ;
	    }
	}
      node->down = 0 ;
    }
      
  if (node->down)
    ok |= bqlAtomize (bql, node->down, aa, type) ;
  if (node->right)
    ok |= bqlAtomize (bql, node->right, aa, type) ;

  return ok ;
} /* bqlAtomize */

/*************************************************************************************/

static NODE *bqlNewFrom (BQL *bql, Array froms, NODE *node)
{
  char buf[64] ;
  NODE *in = (NODE *) halloc (sizeof (NODE), bql->h) ;
  NODE *var = (NODE *) halloc (sizeof (NODE), bql->h) ;
  NODE *var2 = (NODE *) halloc (sizeof (NODE), bql->h) ;
  NODE *up = node->up ;
  SS *s1 ;
  BOOL isDown = (up && up->down == node) ? TRUE : FALSE ;
  int nsc ;

  in->type = IN ; 
  in->right = node ;
  node->up = in ;
  in->down = var2 ;
  var2->up = in ;
    
  var->type = var2->type = VAR ;
  var->dclNode = var2->dclNode = var2 ;
  var->mark = var2->mark = stackMark (bql->s) ;

  nsc = ++bql->dummyVar ;
  sprintf (buf, "__%d", nsc) ;
  pushText (bql->s, buf) ;
  dictAdd (bql->dict, buf, &(var->var)) ;
  var2->var = var->var ;
  arrayp (bql->dcls, var->var, DCL)->dclNode = var2 ;

  if (up)
    {
      if (isDown)
	up->down = var ;
      else
	up->right = var ;
    }
  var->up = up ;
  s1 = arrayp (froms, arrayMax (froms), SS) ;
  s1->node = in ; in->up = 0 ;
 
  return var ;
} /* bqlNewFrom */

/*************************************************************************************/

static BOOL bqlAtomizeWhere (BQL *bql, Array froms, NODE *node)
{
  BOOL ok = FALSE ;

  switch (node->type)
    {
    case SQUARE:
    case CURLY:
      return ok ;
    case ACTIVE:
      break ;
    case DNA:
    case PEPTIDE:
    case CODING:
    case DATEDIFF:
    case TAG:
      if (node->up)
	switch (node->up->type)
	  {
	  case COUNT:
	  case MIN:
	  case MAX:
	  case SUM:
	  case AVERAGE:
	  case STDDEV:
	    break ;
	  default:
	    if (1)
	      {
		ok = TRUE ;
		bqlNewFrom (bql, froms, node) ;
	      }
	    break ;
	  }
      break ;
    case NUMBER:
    case VAR:
      break ;
    default:
      break ;
    }

  if (node->down)
    ok |= bqlAtomizeWhere (bql, froms, node->down) ;
  if (node->right)
    ok |= bqlAtomizeWhere (bql, froms, node->right) ;
  return ok ;
} /* bqlAtomizeWhere */

/*************************************************************************************/
/* atomize everything except atoms */
static NODE *bqlDoAtomizeFrom (BQL *bql, Array froms, NODE *node)
{
  NODE *var = 0 ;
  if (node)
    switch (node->type)
      {
      case ACTIVE:
      case NUMBER:
      case VAR:
      case DOLLAR:
      case QUOTE:
      case DOUBLEQUOTE:
      case COMA:
      case SEMICOLUMN:
	break ;
      case PLUS:
      case MINUS:
      case MULT:
      case DIVIDE:
      case POWER:
      case MODULO:
      case COUNT:
      case MIN:
      case MAX:
      case SUM:
      case AVERAGE:
      case IN:
      case TAG:
      case TTAG:
      case RIGHTOF:
      case OBJECT:
	var = bqlNewFrom (bql, froms, node) ;
	break ;
      default:
	break ;
      }  
  return var ;
} /* bqlDoAtomizeFrom */

/*************************************************************************************/
static BOOL bqlAtomizeFrom (BQL *bql, Array froms, NODE *node)
{
  BOOL ok = FALSE ;

  switch (node->type)
    {
    case SEMICOLUMN:
    case CURLY:
      return ok ;
    case DNA:
    case PEPTIDE:
    case CODING:
    case DATEDIFF:
      if (node->right && node->right->type == COMA)
	{
	  bqlDoAtomizeFrom (bql, froms, node->right->down) ; 
	  if (node->right->right &&  node->right->right->type == COMA)
	    {
	      bqlDoAtomizeFrom (bql, froms, node->right->right->down) ; 
	      bqlDoAtomizeFrom (bql, froms, node->right->right->right) ; 
	    }
	  else
	    bqlDoAtomizeFrom (bql, froms, node->right->right) ; 
	}
      break ;
    case CLASSE:
    case TAG:
    case TTAG:
      if (node->down)
	{
	  switch (node->down->type)
	    {
	    case TAG:
	    case TTAG:
	    case RIGHTOF:
	    case OBJECT:
	    case CLASSE:
	      bqlNewFrom (bql, froms, node->down) ;
	      break ;
	    default:
	      break ;
	    }
	}
     if (node->right)
	{
	  switch (node->right->type)
	    {
	    case CLASSE:
	    case TAG:
	    case TTAG:
	    case MERGETAG:
	      if (1)
		{
		  NODE *S = node ;
		  NODE *T = node->right ;
		  NODE *a = node->down ;
		  NODE *b = T->down ;
		  NODE *c = T->right ;
		  BQLTYPE t = S->type ;

		  ok = TRUE ;
		  S->type = T->type ;
		  T->type = t ;
		  T->right = b ; if (b) b->up = T ;
		  T->down = a ; if (a) a->up = T ;
		  S->down = T ; if (T) T->up = S ;
		  S->right = c ; if (c) c->up = S ;
		  bqlNewFrom (bql, froms, T) ;
		}

	      break ;
	    default:
	      break ;
	    }
	}
      break ;
    default:
      break ;
    }

  if (node->down)
    ok |= bqlAtomizeFrom (bql, froms, node->down) ;
  if (node->right)
    ok |= bqlAtomizeFrom (bql, froms, node->right) ;
  return ok ;
} /* bqlAtomizeFrom */

/*************************************************************************************/
/* Transfer all the where clauses down of from select
 */
static BOOL bqlDesintegrate (BQL *bql, NODE *node0)
{
  AC_HANDLE h = ac_new_handle () ;
  Array selects = arrayHandleCreate(64, SS, h) ; 
  Array wheres = arrayHandleCreate(64, SS, h) ; 
  Array froms = arrayHandleCreate(64, SS, h) ; 
  int nn = 0, level = 0 ;
  BitSet bb = bitSetCreate (64, h) ;
  BitSet bbw = bitSetCreate (64, h) ;
  BOOL ok = TRUE ;

  /* prepare a declaration variable @@ to replace @ */
  if (1)
    {
      NODE *var = (NODE *) halloc (sizeof (NODE), bql->h) ;
      const char *ccp ;
      var->type = VAR ;
      var->var = arrayMax (bql->dcls) ;
      var->mark = stackMark (bql->s) ;
      bql->dummyVar++ ;
      ccp = hprintf(bql->h,"@@%d", bql->dummyVar) ;
      pushText (bql->s, ccp) ; 
      dictAdd (bql->dict, ccp, (&var->var)) ;
      var->dclNode = var ;
      arrayp (bql->dcls, var->var, DCL)->dclNode = var ;

      if (bqlFromActive (bql, node0, var))
	{
	  SS *s1 = arrayp (froms, 0, SS) ;
	  NODE *in = (NODE *) halloc (sizeof (NODE), bql->h) ;
	  NODE *a = (NODE *) halloc (sizeof (NODE), bql->h) ;

	  s1->node = var->up = a->up = in ; s1->node->up = 0 ;
	  in->down = var ;
	  in->right = a ;
	  a->type = ACTIVE ;
	  in->type = IN ;
	}
    }

  while (ok)
    {
      /* Isolate recursivelly the select from where branches */
      int i ;
      ok = FALSE ;

      ok |= bqlAtomize (bql, node0, selects, SELECT) ;
      ok |= bqlAtomize (bql, node0, froms, FROM) ;
      ok |= bqlAtomize (bql, node0, wheres, WHERE) ;
      ok |= bqlAtomize (bql, node0, wheres, WHERE_IF) ;

      for (i = 0 ; i < arrayMax (selects) ; i++)
	{
	  ok |= bqlAtomize (bql, arrayp (selects, i, SS)->node, selects, SELECT) ;
	  ok |= bqlAtomize (bql, arrayp (selects, i, SS)->node, froms, FROM) ;
	  ok |= bqlAtomize (bql, arrayp (selects, i, SS)->node, wheres, WHERE) ;
	  ok |= bqlAtomize (bql, arrayp (selects, i, SS)->node, wheres, WHERE_IF) ;
	}
      for (i = 0 ; i < arrayMax (froms) ; i++)
	{
	  ok |= bqlAtomize (bql, arrayp (froms, i, SS)->node, selects, SELECT) ;
	  ok |= bqlAtomize (bql, arrayp (froms, i, SS)->node, froms, FROM) ;
	  ok |= bqlAtomize (bql, arrayp (froms, i, SS)->node, wheres, WHERE) ;
	  ok |= bqlAtomize (bql, arrayp (froms, i, SS)->node, wheres, WHERE_IF) ;
	}
      for (i = 0 ; i < arrayMax (wheres) ; i++)
	{
	  ok |= bqlAtomize (bql, arrayp (wheres, i, SS)->node, selects, SELECT) ;
	  ok |= bqlAtomize (bql, arrayp (wheres, i, SS)->node, froms, FROM) ;
	  ok |= bqlAtomize (bql, arrayp (wheres, i, SS)->node, wheres, WHERE) ;
	  ok |= bqlAtomize (bql, arrayp (wheres, i, SS)->node, wheres, WHERE_IF) ;
	}
  
      /* creata a new variable for all the dummy cases inside the selects
       * by transferring the subtree in the froms array
       */

      for (i = 0 ; i < arrayMax (selects) ; i++)
	{
	  SS *s1, *ss = arrayp (selects, i, SS) ; 
	  if (ss->node && 
	      (ss->node->type == IN || ss->node->type == SET) &&
	      ss->node->down && ss->node->down->type == VAR
	      )
	    {
	      NODE *var = (NODE *) halloc (sizeof (NODE), bql->h) ;
	      NODE *in = ss->node ;
              s1 = arrayp (froms, arrayMax(froms), SS) ;
	      s1->node = in ; in->up = 0 ;
	      
	      ok = TRUE ;
	      ss->node = var ;
	      *var = *(in->down) ;
	      var->dclNode = in->down->dclNode = in->down ;
	      var->up = 0 ;
	    }	    
	  else if (ss->node && ss->node->type != VAR)
	    {
	      ok = TRUE ;
	      ss->node = bqlNewFrom (bql, froms, ss->node) ;
	    }
	}
    }
  if (1)
    {
      /* split where on AND 
       */
      int i ;
      BOOL ok = TRUE ; /* trun at least once */
      while (ok) /* as often as something change */
	{
	  for (i = 0, ok = FALSE ; ! ok && i < arrayMax (wheres) ; i++)
	    {  /* break if somethng has changed */
	      SS *ss = arrayp (wheres, i, SS) ; 
	      if (ss->node && ss->node->type == AND)
		{
		  NODE *d = ss->node->down ;
		  NODE *r = ss->node->right ;

		  if (d && r && d->type && r->type)
		    {
		      ok = TRUE ;  /* split where a AND B into where A ; where B */
		      ss->node = d ; d->up = 0 ;
		      SS *s1 = arrayp (wheres, arrayMax (wheres), SS) ;
		      s1->node = r ; r->up = 0 ;
		    }
		}
	    }
	  /* break if nothing has changed */
	}
    }

  if (1)
    {
      /* create a new variable for all pieces inside a where clause 
       */
      int i ;
      for (i = 0 ; i < arrayMax (wheres) ; i++)
	{
	  SS *ss = arrayp (wheres, i, SS) ; 
	  if (ss->node)
	    bqlAtomizeWhere (bql, froms, ss->node) ;
	  ss->type = WHERE ;
	}
    }

  if (1)
    {
      /* create a new variable for all pieces double dereferencing from clause
       */
      int i ;
      for (i = 0 ; i < arrayMax (froms) ; i++)
	{
	  SS *ss = arrayp (froms, i, SS) ; 
	  if (ss->node)
	    bqlAtomizeFrom (bql, froms, ss->node) ;
	}
    }

  if (arrayMax(wheres))
    {         /* bitSet the dependencies */
      int ii = 0 ;
      for (ii = 0 ; ii < arrayMax (wheres) ; ii++)
	{
	  SS *ss = arrayp (wheres, ii, SS) ;
	  ss->dcls = bitSetCreate (24, h) ;
	  ss->uses = bitSetCreate (24, h) ;
	  bqlSetDcls (ss->dcls, ss->node) ;
	  bqlSetUses (ss->uses, ss->node) ;
	  
	  bitSetMINUS (ss->uses, ss->dcls) ; /* each node may freely use its own dcls */
	  if (bitSetCount (ss->uses))
	    ss->level = 1 ;
	}
    }

  if (arrayMax (froms))
    {      /* bitSet the dependencies */
      int ii = 0 ;
      for (ii = 0 ; ii < arrayMax (froms) ; ii++)
	{
	  SS *ss = arrayp (froms, ii, SS) ;
	  ss->dcls = bitSetCreate (24, h) ;
	  ss->uses = bitSetCreate (24, h) ;
	  bqlSetDcls (ss->dcls, ss->node) ;
	  bqlSetUses (ss->uses, ss->node) ;
	  
	  bitSetMINUS (ss->uses, ss->dcls) ;  /* each node may freely use its own dcls */
	  if (bitSetCount (ss->uses))
	    nn = ss->level = 1 ;
	  if (ss->node && ss->node->type == IN && ss->node->right && ss->node->right->type == CLASSE)
	    ss->level = 2 ;
	}
    }

  /* order the froms by level */
  nn = 1 ; level = 0 ;
  while (level++, level < 128 && nn)
    {
      int ii, jj ;
      nn = 0 ;
      for (ii = 0 ; ii < arrayMax (froms) ; ii++)
	{
	  SS *s1 = arrp (froms, ii, SS) ;
	  if (s1->level >= level)
	    {
	      if (! bitSetCount (s1->dcls))
		continue ;
	      bitSetOR (bb, s1->dcls) ;
	      bitSetAND (bb, s1->dcls) ;
	      for (jj = 0 ; jj < arrayMax (froms) ; jj++)
		{
		  SS *s2 = arrp (froms, jj, SS) ;
		  if (ii == jj) continue ;
		  if (bitSetAND (bb, s2->uses))
		    { nn++ ; if (s2->level < s1->level + 1) s2->level = s1->level + 1 ; }
		  bitSetOR (bb, s1->dcls) ;
		}
	    }
	}
    }
  arraySort (froms, ssLevelOrder) ;
  if (arrayMax (froms))
    {  /* by now froms only depend on previous froms, set the level */
      int ii ;
      for (ii = 0 ; ii < arrayMax (froms) ; ii++)
	{
	  SS *s1 = arrp (froms, ii, SS) ;
	  s1->level = ii ;
	}
    }

  if (arrayMax (froms))
    {  /* establish the where dependencies */
      int ii, jj ;
      for (ii = 0 ; ii < arrayMax (froms) ; ii++)
	{
	  SS *s1 = arrp (froms, ii, SS) ;
	  if (!s1->node || ! s1->dcls)
	    continue ;
	  bitSetOR (bb, s1->dcls) ;
	  bitSetAND (bb, s1->dcls) ;
	  for (jj = 0 ; jj < arrayMax (wheres) ; jj++)
	    {
	      SS *s2 = arrp (wheres, jj, SS) ;
	      if (bitSetAND (bb, s2->uses))
		{ 
		  bitSetOR (bbw, s2->uses) ;
		  bitSetMINUS (bbw, s1->dcls) ;
		  
		  if (s2->level < s1->level ||
		      ! bitSetCount (bbw))
		    s2->level = s1->level ;
		}
	      bitSetOR (bb, s1->dcls) ;
	    }
	}
    }
  /* by now each where can be hooked back to its lowest from */
  if (arrayMax (wheres))
    {
      int jj ;

      for (jj = 0 ; jj < arrayMax (wheres) ; jj++)
	{
	  SS *s2 = arrp (wheres, jj, SS) ;
          int level = s2->level < arrayMax (froms) ? s2->level : arrayMax(froms) - 1 ;
	  SS *s1 = arrp (froms, level, SS) ;

	  if (! s1->node)
	    continue ;
	  if (s1->node->type != WHERE && s1->node->type != WHERE_IF)
	    {
	      NODE *where = (NODE *) halloc (sizeof (NODE), bql->h) ;
	      where->type = s2->type ; /* WHERE ; */
	      where->down = s1->node ;
	      where->right = s2->node ;
	      if (where->right)  where->right->up = where ;
	      if (where->down) where->down->up = where ;
	      s1->node = where ; s1->node->up = 0 ;
	    }
	  else
	    {
	      NODE *and = (NODE *) halloc (sizeof (NODE), bql->h) ;

	      and->type = AND ;
	      and->down = s1->node->right ;
	      and->right = s2->node ;
	      if (and->right) and->right->up = and ;
	      if (and->down) and->down->up = and ;
	      s1->node->right = and ; and->up = s1->node ;
	    }
	}
    }
  /* hook back the froms */
  if (1)
    {  /* by now froms only depend on previous froms, set the level and hook it */
      int ii ;
      NODE *from = node0 ;
      
      from->right = 0 ;
      for (ii = arrayMax (froms) ;  ii >= 1 ; ii--)
	{
	  SS *s1 = arrp (froms, ii - 1, SS) ;

	  if (s1->node)
	    {
	      NODE *fComa = (NODE *) halloc (sizeof (NODE), bql->h) ;
	      if (from->right)
		{ 
		  fComa->right = from->right ;
		  fComa->right->up = fComa ;
		}
	      from->right = fComa ;
	      fComa->up = from ;
	      fComa->type = COMA ;	  
	      fComa->down = s1->node ; s1->node->up = fComa ;
	    }
	}
    }
  /* hook bact the selects */
  if (1)
    {  /* by now froms only depend on previous froms, set the level and hook it */
      int ii ;
      NODE *select = (NODE *) halloc (sizeof (NODE), bql->h) ;
      NODE *from = node0 ;

      select->up = from ; from->down = select ;
      select->type = SELECT ;
      for (ii = arrayMax (selects) ;  ii >= 1 ; ii--)
	{
	  SS *s1 = arrp (selects, ii - 1, SS) ;
	  NODE *sComa = (NODE *) halloc (sizeof (NODE), bql->h) ;

	  if (select->right)
	    { 
	      sComa->right = select->right ;
	      sComa->right->up = sComa ;
	    }
	  select->right = sComa ;
	  sComa->up = select ;
	  sComa->type = COMA ;	  
	  sComa->down = s1->node ; s1->node->up = sComa ;
	}
    }
  
  ac_free (h) ;
  return TRUE ;
} /*  bqlDesintegrate */

/****************/

static BOOL bqlDezingue (BQL *bql, NODE* node)
{
  if (
      (node->type == WHERE || node->type == WHERE_IF) &&
      node->right &&
      node->right->type == COUNT
      )
    { /* transform where COUNT{} into where COUNT {} > 0*/
      NODE *new = (NODE *) halloc (sizeof (NODE), bql->h) ;
      NODE *zero = (NODE *) halloc (sizeof (NODE), bql->h) ;
      new->type = LT ;
      new->right = zero ; new->down = node->right ;
      node->right->up = zero->up = new ;
      node->right = new ;
      zero->type = NUMBER ; zero->z = 0 ;
    }
  if (node->down)
    bqlDezingue (bql, node->down) ;
  if (node->right)
    bqlDezingue (bql, node->right) ;
  if (node->up && node->up->type == SEMICOLUMN && node->up->down == node)
    bqlDesintegrate (bql, node) ; 
  if (node->up && node->up->type == CURLY && node->up->down == node)
    bqlDesintegrate (bql, node) ; 
  return TRUE ;
} /* bqlDezingue */

/****************/

static BOOL bqlSetComas (BQL *bql)
{
  NODE *semicol = semicol = bql->node ;
  while (semicol && semicol->type == SEMICOLUMN)
    {
      if (semicol->down)
	{
	  bqlDoSetComas (bql, semicol->down) ;
	  bqlDezingue (bql, semicol->down) ;
	}
      semicol = semicol->right ;
    }

  return TRUE ;
} /*  bqlCreateStComas */

/*************************************************************************************/

static BOOL bqlSetSemiCols  (BQL *bql)
{
  NODE *semicol = 0 ;

  if (! bql->node)
    return FALSE ;
  if (bql->node->type != SEMICOLUMN)
    {
      NODE *new = (NODE *) halloc (sizeof (NODE), bql->h) ;
      new->type = SEMICOLUMN ;
      new->down = bql->node ;
      bql->node->up = new ;
      bql->node = new ;
    }
  semicol = bql->node ;
  while (semicol && semicol->type == SEMICOLUMN)
    semicol = semicol->right ;
  if (semicol && semicol->type != SEMICOLUMN)
    {
      NODE *new = (NODE *) halloc (sizeof (NODE), bql->h) ;
      new->type = SEMICOLUMN ;
      new->down = semicol ;
      new->up = semicol->up ;
      if (new->up->down == semicol)
	new->up->down = new ;
      if (new->up->right == semicol)
	new->up->right = new ;
      semicol->up = new ;
    }
  return TRUE ;
} /*  bqlSetSemiCols */

/*************************************************************************************/

static void bqlCatQuoteWord (BQL *bql, char *cs)
{
  char cc, *cq ;
  if (! cs || ! *cs)
    return ;
  while (*cs == ' ') cs++ ;
  cq = cs + strlen (cs) - 1 ;
  while (*cq == ' ') *cq-- = 0 ;
  cq = cs ;
  if (*cs == '\"' || *cs == '\'')
    {
      catText (bql->s, cs) ;
      return ; 
    }

  while (*cq && *cq != ' ') cq++ ;
  cc = *cq ; *cq = 0 ;
  catText (bql->s, "\"") ;
  catText (bql->s, cs) ;
  catText (bql->s, "\"") ;
  if (! cc)
    return ;

  catText (bql->s, " ") ;
  catText (bql->s, cq + 1) ;

  return ;
} /* bqlCatQuoteWord */

/*************************************************************************************/
/* expand allowed abbreviations
 * i.e.  cases with implicit select from where
 */
static BOOL bqlCreatePhonyVariables  (BQL *bql, NODE *node0)
{
  BOOL debug = bql->debug ;
  char *cp = 0, *cz = 0 ;
  if (node0->down)
    bqlCreatePhonyVariables (bql, node0->down) ;
  if (node0->right)
    bqlCreatePhonyVariables (bql, node0->right) ;

  if ( node0->type == SEMICOLUMN)
    return TRUE ;

  if ( node0->type == CURLY)
    return TRUE ;

  if (node0->up)
    {
      switch (node0->up->type)
	{
	case SEMICOLUMN: break ;
	case CURLY:
	  if (node0->up->down == node0) 
	    break ;
	  /* else fall through */
	default:
	  return TRUE ;
	}
    }
  cp = stackText (bql->s, node0->mark) ;
  while (*cp == ' ') cp++ ;
  if (cp && strncasecmp (cp, "select ", 7))  /* assume select */
    {
      char *cs = strnew (cp, bql->h) ; /* because the bql->s may be reallocated */
      node0->mark = stackMark (bql->s) ;
      pushText (bql->s, "select ") ;
      catText (bql->s, cs) ;
      cp = stackText (bql->s, node0->mark) ;
    }
  cp += 7 ;  
  while (*cp == ' ') cp++ ;
  if (cp && ! strncasecmp (cp, "Find ", 5))  /* assume select */
    {
      memset (cp, ' ', 4) ; cp+=4 ; *cp = '?' ;
    }
  else if (cp && ! strncasecmp (cp, "class ", 6))  /* assume select */
    {
      memset (cp, ' ', 5) ; cp+=5 ; *cp = '?' ;
    }

  else if (cp && ! strncasecmp (cp, "isa ", 4))  /* assume select */
    {
      char *cs = strnew (cp + 4, bql->h) ; /* because the bql->s may be reallocated */
      node0->mark = stackMark (bql->s) ;
      pushText (bql->s, hprintf (bql->h, "select @ where @ ISA %s", cs)) ;
    }

  if (cp && *cp == '>' && (cp[1] == ' ' || cp[1] == '\"' || cp[1] == '=' || cp[1] == '~')) /* select >= x   ->   select @ where @ >= x */
    {
      char *cs = strnew (cp, bql->h) ; /* because the bql->s may be reallocated */
      node0->mark = stackMark (bql->s) ;
      pushText (bql->s, hprintf (bql->h, "select @ where @ %c%c", cs[0], cs[1])) ;
      bqlCatQuoteWord (bql, cs + 2) ;
    }
  else if (cp && *cp == '=' && (cp[1] == '=' || cp[1] == '~')) /* select >= x   ->   select @ where @ >= x */
    {
      char *cs = strnew (cp, bql->h) ; /* because the bql->s may be reallocated */
      node0->mark = stackMark (bql->s) ;
      pushText (bql->s,  hprintf (bql->h, "select @ where @  %c%c", cs[0], cs[1])) ;
      bqlCatQuoteWord (bql, cs+2) ;
    }
  else if (cp && (*cp == '-' || *cp == '>' || *cp == '=' ) && cp[1] == '>') /* select =>tag   ->   select @=>tag */
    {
      char *cs = strnew (cp, bql->h) ; /* because the bql->s may be reallocated */
      node0->mark = stackMark (bql->s) ;
      pushText (bql->s,  hprintf (bql->h, "select @ %c%c", cs[0], cs[1])) ;
      catText (bql->s, cs + 2) ;
    }
  else if (cp && cp[1] == '=' && (cp[0] == '!' || cp[0] == '>' || cp[0] == '<')) /* select tida or ~= or ~> */
    {
      char *cs = strnew (cp, bql->h) ; /* because the bql->s may be reallocated */
      node0->mark = stackMark (bql->s) ;
      pushText (bql->s,  hprintf (bql->h, "select @ where @ %c%c", cs[0], cs[1])) ;
      bqlCatQuoteWord (bql, cs+2) ;
    }
  else if (cp && cp[1] == '~' && (cp[0] == '!' || cp[0] == '>' || cp[0] == '<')) /* select tida or ~= or ~> */
    {
      char *cs = strnew (cp, bql->h) ; /* because the bql->s may be reallocated */
      node0->mark = stackMark (bql->s) ;
      pushText (bql->s,  hprintf (bql->h, "select @ where @ %c%c", cs[0], cs[1])) ;
      bqlCatQuoteWord (bql, cs+2) ;
    }
  else if (cp && *cp == '>') /* select >tag   ->   select @->tag */
    {
      char *cs = strnew (cp, bql->h) ; /* because the bql->s may be reallocated */
      node0->mark = stackMark (bql->s) ;
      pushText (bql->s, "select @ ->") ;
      catText (bql->s, cs + 1) ;
    }
  else if (cp && !strncasecmp (cp, "Follow ", 7))
    {
      char *cs = strnew (cp+7, bql->h) ; /* because the bql->s may be reallocated */
      node0->mark = stackMark (bql->s) ;
      pushText (bql->s, "select @->") ;
      catText (bql->s, cs) ;
    }
  else if (cp && (*cp == '~' || *cp == '<')) /* select tida or ~= or ~> */
    {
      char *cs = strnew (cp, bql->h) ; /* because the bql->s may be reallocated */
      node0->mark = stackMark (bql->s) ;
      pushText (bql->s,  hprintf (bql->h, "select @ where @ %c ", cs[0])) ;
      bqlCatQuoteWord (bql, cs+1) ;
    }
 else if (cp &&  !strncasecmp (cp, "IS ", 3))
    {
      char *cs = strnew (cp + 3, bql->h) ; /* because the bql->s may be reallocated */
      node0->mark = stackMark (bql->s) ;
      pushText (bql->s, "select @ where @ ~ ") ;
      bqlCatQuoteWord (bql, cs) ;
    }
  else if (cp && *cp == '#') /* select #tag   ->   select @ where @#tag */
    {
      char *cs = strnew (cp, bql->h) ; /* because the bql->s may be reallocated */
      node0->mark = stackMark (bql->s) ;
      pushText (bql->s,  hprintf (bql->h, "select @ where @  %c ", cs[0])) ;
      catText (bql->s, cs + 1) ;
    }
  cp = stackText (bql->s, node0->mark) ;


  if (1 && cp && ! strstr (cp, "from") && ! strncasecmp (cp, "select ", 7) /* && selectAll */)
    {
      /* assume select x,x->t,... from x in @ */
      char *cp = stackText (bql->s, node0->mark) + 7 ;
      while (*cp == ' ') cp++ ;
      int nw = 0 ;
      
      if (*cp == '?')  /* this would edit ;; bql ?run into  bql s from s in class run */
	{
	  KEY classe = 0 ;
	  char cc, *cq = cp + 1, *cr,*ccl ;
	  int nsc = ++bql->dummyVar ;

	  while (isalnum(*++cq))
	    {}
	  cc = *cq ; *cq = 0 ;
	  ccl = strnew (cp +  1, bql->h) ;
	  *cq = cc ;
	  if (cp[1] && lexword2key(ccl, &classe, _VClass))  
	    {
	      int ntilde = 0, nT = 0 ;
	      char *ncoma = 0 ;
	      cr = cq ;
	      ncoma = strchr (cq, ',') ;

	      while (*cq == ' ') cq++ ;
	      cr = cq ;
	      while (*cr == ' ') cr++ ;
	      if (*cr == '~') 
		{
		  ntilde  = 1 ; cr++ ;
		  while (*cr == ' ') cr++ ;
		}
	      if (*cr) 
		{ 
		  nw = 1 ;
		  while (1)
		    {
		      while (*cr && (isalnum(*cr) || *cr == '_' || *cr == '*' || *cr == '?')) cr++ ;
		      while (*cr == ' ') cr++ ;
		      if (*cr && 
			  (
			   (cr[0] == '-' && cr[1] =='>') ||
			   (cr[0] == '>' && cr[1] =='>') ||
			   (cr[0] == '#') 
			   )
			  )
			{ cr += 2 ;nT = 1; }
		      else break ;
		    }
		  if (*cr) nw++ ;
		}
	      node0->mark = stackMark (bql->s) ;
	      if (ncoma == 0 && ntilde == 0 && nT == 1)
		{
		  while (*cq == ' ') cq++ ;
		  if (*cr)
		    {
		      char *cr2 = strnew (cr, bql->h) ;
		      *cr = 0 ;
		      cz = messprintf ("select _%d from _%d in  class %s where _%d ~ %s", nsc, nsc, cp+1, nsc, ac_protect (lexcleanup (cr2, bql->h), bql->h)) ;
		    }
		  else
		    cz = messprintf ("select _%d from _%d in  class %s", nsc, nsc, cp+1) ;
		}
	      else if (ncoma == 0 && ntilde == 0 && nw == 1 && !strstr (cq, "where"))
		{
		  KEY tag = 0 ;
		  while (*cq == ' ') cq++ ;
		  if (lexword2key(cq, &tag, 0))
		    cz = messprintf ("select _%d from _%d in class %s where (_%d ~ %s || _%d#%s)", nsc, nsc, ccl, nsc, ac_protect (lexcleanup (cq, bql->h), bql->h), nsc, cq) ;
		  else 
		    cz = messprintf ("select _%d from _%d in class %s where _%d ~ %s", nsc, nsc, ccl, nsc, ac_protect (lexcleanup (cq, bql->h), bql->h)) ;
		}
	      else if (ncoma == 0 && ntilde == 1 && nw == 1 && !strstr (cq, "where"))
		{
		  while (*cq == ' ') cq++ ;
		  if (*cq == '~') cq++ ;
		  while (*cq == ' ') cq++ ;
		  cz = messprintf ("select _%d from _%d in class %s where _%d ~ %s", nsc, nsc, ccl, nsc, ac_protect (lexcleanup (cq, bql->h), bql->h)) ;
		}
	      else if (ncoma == 0 && nw == 0 && node0->right && node0->right->type == DOLLAR)
		cz = messprintf ("select _%d from _%d in class %s where _%d ~ ", nsc, nsc, ccl, nsc) ;
	      else if (ncoma == 0 && !strstr (cq, "where"))
		cz = messprintf ("select _%d from _%d in class %s where _%d %s", nsc, nsc, ccl, nsc, cq) ;
	      else
		{
		  char *cfrom ;
		  char *cw, *cs, *ct, *newDna = "" ;
		  int k = 0;
		  /* after the coma we want to edit ', >tag' into ', ?->tag' */
		  if (0 && ncoma)
		    {
		      cs = cr + 1 ;
		      while (*cs == ' ' ) cs++ ;
		      if (*cs == '>')
			{
			  cr = hprintf (bql->h, ", ?-%s", cs) ;
			}
		      else if (*cs == '#')
			{
			  cr = hprintf (bql->h, ", ?%s", cs) ;
			}
		    }

		  /* we need to edit @ into _%d */
		  ct = cs = halloc (1 + 5 * strlen(cq), bql->h) ; /* room to add the nsc numbers after each occurence of ? */
		  
		  cfrom = (char *) strcasestr (cr, "from") ;
		  cr-- ;
		  while (*++cr)
		    switch ((int)*cr)
		      {
		      case ',':
			k = 0 ; 
			*ct++ = *cr ;
			break ;
		      case '>':
			if (k == 0 && (cr < cfrom || ! cfrom))
			  {
			    *ct++ = '_' ; 
			    {
			      char *cv = messprintf ("%d", nsc) ;
			      while (*cv)
				*ct++ = *cv++ ;
			    } 
			    *ct++ = '-' ;
			    *ct++ = '>' ;
			  }
			else 
			  *ct++ = *cr ;
			k = 1 ;
			break ;
		      case '#':
			if (k == 0 && (cr < cfrom || ! cfrom))
			  {
			    *ct++ = '_' ; 
			    {
			      char *cv = messprintf ("%d", nsc) ;
			      while (*cv)
				*ct++ = *cv++ ;
			    } 
			    *ct++ = '#' ;
			  }
			else 
			  *ct++ = *cr ;
			k = 1 ;
			break ;
		      case '?':
			*ct++ = '_' ; 
			{
			  char *cv = messprintf ("%d", nsc) ;
			  while (*cv)
			    *ct++ = *cv++ ;
			} 
			k = 1 ;
			break ;
		      case 'D':
			if (k == 0 && (cr < cfrom || ! cfrom) && cr[1] == 'N' && cr[2] == 'A' )
			  {
			    char *ct1 = ct ;
			    char *ct2 = cr ;
			    char cc3, *ct3 = strchr (cr, ',') ;
			    
			    if (ct3 && cfrom && ct3 > cfrom)
			      ct3 = cfrom ;
			    if (!ct3 && cfrom)
			      ct3 = cfrom ;
			    if (ct3)
			      { cc3 = *ct3 ; *ct3 = 0 ; }
			    *ct++ = '_' ; 
			    *ct++ = 'd' ;
			    {
			      char *cv = messprintf ("%d", nsc) ;
			      while (*cv)
				*ct++ = *cv++ ;
			    } 
			    newDna = hprintf (bql->h, "%s , %s in %s (_%s)", newDna, ct1, ct2, ct1+2) ;
			    if (ct3)
			      { *ct3 = cc3 ; cr = ct3 ; }
			    else
			      cr += 2 ;
			    k = 1 ;
			    break ;
			  }
			else 
			  *ct++ = *cr ;
			k = 1 ;
			break ;
		      case ' ':
			*ct++ = *cr ;
			break ;
		      default:
			*ct++ = *cr ;
			k = 1 ;
			break ;
		      }
		  *ct = 0 ;
		  cq = cs ; /* the modified string */

		  cw = strstr (cq, " where") ;
		  if (cw) *cw++ = 0 ; else cw = " " ;
		  if (ncoma == 0)
		    cz = messprintf ("select _%d from _%d in class %s %s %s", nsc, nsc, ccl, cq, cw) ;
		  else
		    cz = messprintf ("select _%d %s from _%d in class %s %s %s", nsc, cq, nsc, ccl, cw, newDna) ;
		}
	      pushText (bql->s, cz) ;
	      
	    }
	}
      
      cp = stackText (bql->s, node0->mark) ;
      cp += 7 ;
      nw = 1 ;
      if (node0->right && node0->right->type != SEMICOLUMN && node0->right->type != DOLLAR)
	nw++ ;
      if (nw == 1)
	{
	  char *cr = cp ;
	  if (*cr) 
	    { 
	      if (*cr >= '0' && *cr <= '9') nw++ ; /* a variable cannot start by a number */
	      while (*cr && (isalnum(*cr)  || *cr == '{' || *cr == '[' || *cr == '(' || *cr == '_' || *cr == '*' || *cr == '?')) cr++ ;
	      while (*cr == ' ') cr++ ;
	      if (*cr) nw++ ;
	    }
	}

      if (nw == 1 && ! strstr (cp, "from") && ! strstr (cp, "where"))
	{
	  int nsc = ++bql->dummyVar ;
	  KEY tag = 0 ;
	  char *cs = strnew (cp, bql->h) ; /* because the bql->s may be reallocated */
	  node0->mark = stackMark (bql->s) ;
	  if (! strncmp (cs, "DNA", 3))
	    pushText (bql->s, messprintf ("select _DNA%d  from _%d in @, _DNA%d in DNA(_%d)", nsc, nsc, nsc, nsc)) ;
	  else if (! strncmp (cs, "PEPTIDE", 7))
	    pushText (bql->s, messprintf ("select _PEPTIDE%d  from _%d in @, _PEPTIDE%d in PEPTIDE(_%d)", nsc, nsc, nsc, nsc)) ;
	  else if (! strncmp (cs, "CODING", 6))
	    pushText (bql->s, messprintf ("select _CODING%d  from _%d in @, _CODING%d in CODING(_%d)", nsc, nsc, nsc, nsc)) ;
	  else if (lexword2key(cs, &tag, 0))
	    pushText (bql->s, messprintf ("select _%d  from _%d in @ where (_%d #%s || _%d ~ %s)", nsc, nsc, nsc, cs, nsc, ac_protect (lexcleanup (cs, bql->h), bql->h), cs)) ;
	  else
	    pushText (bql->s, messprintf ("select _%d from _%d in @ where _%d ~ %s"
					  , nsc, nsc, nsc
					  , ac_protect (lexcleanup (cs, bql->h), bql->h)) 
		      ) ;
	  cp = stackText (bql->s, node0->mark) ;
	}
      else if (nw == 2 &&
	       (*cp == '!' || !strncasecmp (cp, "not ", 4)) &&
	       ! strstr (cp, "from") && ! strstr (cp, "where")
	       )
	{
	  int nsc = ++bql->dummyVar ;
	  KEY tag = 0 ;
	  char *cs = strnew (cp, bql->h) ; /* because the bql->s may be reallocated */
	  node0->mark = stackMark (bql->s) ;
	  if (*cs == '!') cs++ ;
	  else cs += 4 ;
	  
	  while (*cs == ' ') cs++ ;
	  if (lexword2key(cs, &tag, 0))
	    pushText (bql->s, messprintf ("select _%d  from _%d in @ where NOT (_%d #%s || _%d ~ %s)", nsc, nsc, nsc, cs, nsc, ac_protect (lexcleanup (cs, bql->h), bql->h), cs)) ;
	  else
	    pushText (bql->s, messprintf ("select _%d from _%d in @ where NOT (_%d ~ %s)"
					  , nsc, nsc, nsc
					  , ac_protect (lexcleanup (cs, bql->h), bql->h))
		      ) ;
	  cp = stackText (bql->s, node0->mark) ;
	}
    }
      
  cp = stackText (bql->s, node0->mark) ;
  if (cp && cp == strstr (cp,"select ") &&  /* 2018_09_03 protect  select ?paper ; "k*" */
      node0->type == DOLLAR &&
      node0->up &&
      node0->up->type == SEMICOLUMN)
    node0->type = RAW ;

 if (debug && cz)
    fprintf (stderr, "....... %s\n", cz) ;
  return TRUE ;
} /* bqlCreatePhonyVariables */

/*************************************************************************************/

static BOOL bqlCheckSides (BQL *bql, NODE *node)
{
  BOOL ok = TRUE ;
  NODE *down = node->down ;
  NODE *right = node->right ;
  int side = bqlSide[node->type] ;

  if (side & 0x2 && ! down)
    {
      vtxtPrintf (bql->errTxt, "// ... Bad parsing, nothing left of %s\n", bqlName[node->type]) ;
      ok = FALSE ;
    }
  if (side & 0x8 && ! right)
    {
      vtxtPrintf (bql->errTxt, "// ... Bad parsing, nothing right of %s\n", bqlName[node->type]) ;
      ok = FALSE ;
    }

  if (side && ok && down) ok = bqlCheckSides (bql, down) ;  
  if (side && ok && right) ok = bqlCheckSides (bql, right) ;  
  return ok ;
}

/*************************************************************************************/

static BOOL bqlCheckSyntax (BQL *bql, NODE *node0, int pass)
{
  BOOL ok = TRUE ;
  NODE *node, *selectNode ;
  
  /* select from are mandatory */
  node = node0 ;
  while (ok && node && node->type == SEMICOLUMN)
    {
      NODE *down = node->down ;
      if (down)
	ok = bqlCheckSyntax (bql, down, pass) ;
      node = node->right ;
    }
  if (! node)
    return ok ;
  if (node->type == SELECT && node->right && node->right->type == FROM)
    {
      NODE *fromNode = node->right ;
      NODE *downNode = fromNode->down ;

      bql->node = fromNode ;
      fromNode->up = 0 ;
      node->right = downNode ;
      node->up = fromNode ;
      fromNode->down = node ;
      if (downNode) downNode->up = node ;
    }

  if (0 && node->type != FROM)
    {
      vtxtPrintf (bql->errTxt, "// ... from clause Missing, all BQL queries should be of the form select ... from ...") ;
      return FALSE ;      
    }
  if (0 && (! node->down || node->down->type != SELECT))
    {
      vtxtPrintf (bql->errTxt, "// ... select clause Missing, all BQL queries should be of the form select ... from ...") ;
      return FALSE ;      
    }
  if (node->type == ORDER_BY)
    node = node->down ;
  if (node->type == FROM && node->down && node->down->type == SELECT)
    selectNode = node->down ;
  else
    selectNode = node ;

  if (node->type == VAR 
      && node->up 
      && node->up->type == CLASSE
      && node->mark
      && ! node->key
      )
    {
      char *cp = stackText (bql->s, node->mark) ;
      if (1)
	{
	  char cc, *cq = cp ;
					
	  KEY classe = 0 ;
	  while (*++cq && *cq != ' ')
	    {}
	  cc = *cq ; *cq = 0 ;
	  if (*cp && lexword2key(cp, &classe, _VClass))  
	    {
	      unsigned char cc = 0 ;

	      node->var = classe ;
	      node->key = goodClass (classe, &cc) ;
	      node->uType = cc ;
	    }
	  else
	    {
	      vtxtPrintf (bql->errTxt, "// ... Class %s : unrecognized class"
			  , cp) ;
	      return FALSE ; 
	    }
	  *cq = cc ;
	}
    }
     

  
  /* select should contain a vertical list of coma, with all leaves being variables */
  
  for (node = selectNode->down ; ok && node ; node = node->down)
    {
      NODE *leaf ;
      if (node->type == COMA)
	leaf = node->right ;
      else
	leaf = node ;
      if (leaf->type == ROUND && ! leaf->down && ! leaf->right)
	continue ;
      else
	{
	  vtxtPrintf (bql->errTxt, "// ... There is something left of the select clause\n") ;
	  bqlShowNode (bql, selectNode, 0, 0) ;
	  return FALSE ;      
	}
      
      for (node = selectNode->right ; ok && node ; node = node->down)
	{
	  NODE *leaf ;
	  if (node->type == COMA)
	    leaf = node->right ;
	  else
	    leaf = node ;
	  if (leaf->type == ROUND && ! leaf->down && ! leaf->right)
	    continue ;
	  if (leaf->right || leaf->down || leaf->type != VAR)
	    {
	      vtxtPrintf (bql->errTxt, "// ... BQL ERROR 4: select clause should be of the form select a, b, c where a, b , are variables defined in the from clause\nThere is an error around \n") ;
	      bqlShowNode (bql, node, 0, 0) ;
	      return FALSE ;      
	    }
	  
	}
    }     
  return ok ;
} /* bqlCheckSyntax */

/*************************************************************************************/
/* foreach node->type==var, 
 * set node->dclNode array depth to prepare for x[n] constructs
 */
static BOOL bqlSetDclNodeDepth (BQL *bql, NODE *node)
{ 
  NODE *right = node ? node->right : 0 ;
  NODE *down = node ? node->down : 0 ;
  BOOL ok = TRUE ;

  if (node && node->type == VAR)
    {   
      node->nCol = 1 ;
      if (node->up 
	  && node->up->type == IN 
	  && node->up->right
	  && node->up->right->type == TAG
	  && node->up->right->right
	  && node->up->right->right->type == RIGHTOF
	  && node->up->right->right->right
	  && node->up->right->right->right->type == NUMBER
	  )
	{
	  node->nCol = node->up->right->right->right->z ;
	  if (node->nCol < 0)
	    node->nCol = node->up->right->right->right->z = 0 ;
	}
	  
      if (! node->dclNode)
	{
	  node->dclNode = arrayp (bql->dcls, node->var, DCL)->dclNode ;
	}
      if (! node->dclNode)
	ok = TRUE ;
    }
  if (node && node->type == RIGHTOF)
    {
      NODE *tagNode = node->down ; /* tagNode is the b of "a in b:2" */
      NODE *depthNode = node->right ;

      if (tagNode && depthNode && depthNode->type == NUMBER)
	{
	  int depth = depthNode->z ;
	  NODE *varNode  = 0 ;

	  if (node->up && 
	      node->up->type == IN && 
	      node->up->down &&
	      node->up->down->type == VAR)
	    varNode = node->up->down ;
	  if (varNode)
	    {   /* replace b in a[1], c in a[2], 
		 *      by b in a[1], c in b[1]
		 */
	      NODE *parent ;

	      if (! varNode->dclNode)
		messcrash ("no dclNode") ;
	      varNode->dclNode->col = depth ;
	      parent = tagNode->dclNode ;
	      if (parent)
		{
		  if (parent->child && parent->child->depth < depth)
		    {
		      node->down->dclNode = parent->child ;
		      depth -= parent->child->col ;
		    }
		  
		  parent->child = varNode->dclNode ;
		  varNode->depth = depth ;  /* varNode is the a of "a in b:2" */
		  while (parent->parent)
		    parent = parent->parent ;
		  varNode->parent = parent ;
		  parent->nCol += depth ; /* over estimating is not a problem */
		}
	    }
	}
    }
  if (ok && down)
    bqlSetDclNodeDepth (bql, down) ;
  if (ok && right)
    bqlSetDclNodeDepth (bql, right) ;

  return ok ;
} /* bqlSetDclNodeDepth */

/*************************************************************************************/
/* in case of error: return false, details in bqlError() 
 * *  Sort the lines of the table according to the spec
	*  Return FALSE in case of errors, detailed in errTxt
	*     spec should look like "-2+4+1" which would sort by
	*     column 2 reversed, then column 4, then column 1.
	*  KEY and objects are sorted by their names
	*  Names and strings are sorted acedb way, understanding numbers
	*     James_Bond7 comes before James_Bond12
	*  Numbers are sorted numerically.
	*  Boolean false comes before True
	*  Dates are sorted chronologically
 */
static BOOL bqlSort (BQL *bql)
{
  BOOL ok = FALSE ;
  if (bql && bql->results)
    ok = ac_table_sort_compress (bql->results, bql->order_by, bql->errTxt) ; 
  bql->isSorted = TRUE ;
  return ok ;
}

/*************************************************************************************/

BOOL bqlExport (BQL *bql, char beauty, BOOL showTitle, int beginLine, int maxLine, ACEOUT ao) 
{
  int col ;
  Array cols = arrayHandleCreate (12, int, bql->h) ;
  const char *more = 0 ;
  vTXT txt = vtxtHandleCreate (0) ;
  const char **titles = bql->titles ? arrayp (bql->titles, 0, const char*) : 0 ;

  if (! bql->isSorted && ! bql->doNotSort)
    bqlSort (bql) ;

  for (col = 0 ; titles[col] && col < arrayMax (bql->titles) ; col++)
    array (cols, col, int) = col + 1 ;
  array (cols,  bql->results->cols, int) = 0 ; /* zero terminate */

  
  if (bql->results->rows && bql->results->cols)
    ac_table_display (txt, bql->results
		      , showTitle ? titles : 0
		      , arrayp (cols, 0, int)
		      , 0
		      , beginLine, maxLine, more
		      , beauty
		      ) ;
  if (vtxtPtr (txt))
    {
      if (ao)
	aceOut (ao, vtxtPtr (txt)) ;
      else  /* freeOut is needed on the server side */
	freeOut (vtxtPtr (txt)) ;
    }
  ac_free (txt) ;
  return TRUE ;
} /* bqlExport */

/*************************************************************************************/

void uBqlDestroy (void *vp)
{
  BQL *bql = (BQL *) vp ;
  if (bql)
    {
      if (bql->magic == uBqlDestroy)
	bql->magic = NULL ;
      else
	messcrash ("bqlDestroy received a non BQL pointer") ;
      ac_free (bql->h) ;
    }
} /* uBqlDestroy */

/*************************************************************************************/

void uBqlIterDestroy (void *vp)
{
  BQL_ITER *iter = (BQL_ITER *) vp ;
  if (iter)
    {
      if (iter->magic == uBqlIterDestroy)
	iter->magic = NULL ;
      else
	messcrash ("bqlDestroy received a non BQL_ITER pointer") ;
      ac_free (iter->h) ;
    }
} /* uBqlIterDestroy */

/*************************************************************************************/
/*************************************************************************************/

static void bqlClean (BQL *bql, NODE *node)
{
  if (bql->where) return ;

  if (node->type == VAR)
    {
      if (node->myObj)
	bsDestroy (node->obj) ;
      node->obj = 0 ;
      node->myObj = FALSE ;

      if (node->myAa)
	arrayDestroy (node->aa) ;
      node->aa = 0 ;
      node->myAa = FALSE ;
      ac_free (node->vtxt) ;
      arrayDestroy (node->pep) ;
      stackDestroy (node->dnaStack) ;
      stackDestroy (node->pepStack) ;
      ac_free (node->br) ;
      node->key = 0 ;
      node->uu.s = 0 ;
      node->timeStamp = 0 ;
      node->uType = 0 ;
      node->z = 0 ;
      node->isNumber = FALSE ; 
      node->isDate = FALSE ; 
      node->parent = 0 ;
      if (node->myBsmark)
	ac_free (node->bsMark) ;
      node->bsMark = 0 ; node->myBsmark = FALSE ;
    }
  if (node->down)
    bqlClean (bql, node->down) ;
  if (node->right)
    bqlClean (bql, node->right) ;
  return ;
} /* bqlClean */

/*************************************************************************************/

static BOOL bqlDoExportLine (BQL *bql, BOOL test)
{
  NODE *select ;
  NODE *var ;
  NODE *coma ;
  BOOL ok = FALSE ;
  BOOL isFirst = TRUE ;
  AC_TABLE table = bql->results ;
  int row = test ? bql->tableRow : bql->tableRow++ ;
  int col = -1 ;
  BOOL debug = bql->debug ;

  select = bql->from->down ;
  if (select && select->type == SELECT)
    {
      var = select->right ; coma = 0 ;
      while (var)
	{
	  if (! test) ok = TRUE ;
	  if (var->type == COMA)
	    {
	      coma = var ; var = coma->down ;
	    }
	   if (! var->dclNode)
	     return FALSE ;

	   if (var->dclNode)
	     { 
	       if (row == 0)
		 {
		   char *cr = array (bql->titles, col + 1, char*) ;
		   if (! cr || ! *cr)
		     array (bql->titles, col + 1, char*) = stackText (bql->s, var->dclNode->mark) ;
		 }
	       if (!test)
		 ac_table_insert_type (table, row, ++col, 0, ac_type_empty) ;
	       if (debug) fprintf (stderr, "\t") ;
	       if (var->dclNode->vtxt)
		 {
		   if (vtxtPtr (var->dclNode->vtxt))
		     {
		       if (test) return TRUE ;
		       ac_table_insert_text (table, row, col, vtxtPtr (var->dclNode->vtxt) + 2) ;
		     }
		 }
	       else if (var->dclNode->isNumber)
		 {
		   int z1 = var->dclNode->z ;
		   float z2 = var->dclNode->z ;
		   long long int z3 = var->dclNode->z ;
		   double z4 = var->dclNode->z ;
		   if (z1 == var->dclNode->z)
		     {
		       if (test) return TRUE ;
		       ac_table_insert_type (table, row, col, &z1, ac_type_int) ;
		     }
		   else if (z2 == var->dclNode->z)
		     {
		       if (test) return TRUE ;
		       ac_table_insert_type (table, row, col, &z2, ac_type_float) ;
		     }
		   else if (z3 == var->dclNode->z)
		     {
		       char buf[256] ;
		       if (test) return TRUE ;
		       sprintf (buf, "%lld",  z3) ;
		       ac_table_insert_text (table, row, col, buf) ;
		     }
		   else if (z4 == var->dclNode->z)
		     {
		       char buf[256] ;
		       if (test) return TRUE ;
		       sprintf (buf, "%lg",  z4) ;
		       ac_table_insert_text (table, row, col, buf) ;
		     }
		   if (debug) fprintf (stderr, "%g", var->dclNode->z) ;
		 }
	       else if (var->dclNode->isDate)
		 {
		   mytime_t z1 = var->dclNode->z ;
		   if (test) return TRUE ;
		   ac_table_insert_type (table, row, col, &z1, ac_type_date) ;
		   if (debug) 
		     {
		       char buf[64] ;
		       buf[0] = 0 ;
		       timeShow (z1, buf, 64) ;
		       fprintf (stderr, "%s", buf) ;
		     }
		 }
	       else if (var->dclNode->key)
		 {
		   if (0)
		     {
		       if (test) return TRUE ;
		       ac_table_insert_type (table, row, col, &var->dclNode->key, ac_type_key) ;
		     }
		   if (debug) fprintf (stderr, "%s", name(var->dclNode->key)) ;
		 }
	       else if (var->dclNode->uType)
		 {
		   BSunit uu = var->dclNode->uu ;
		   switch (var->dclNode->uType)
		     {
		     case 0:
		       break ;
		     case _Int:
		       if (test) return TRUE ;
		       ac_table_insert_type (table, row, col, &(uu.i), ac_type_int) ;
		       if (debug) fprintf (stderr, "%d", uu.i);
		       break ;
		     case _Float:
		       if (test) return TRUE ;
		       ac_table_insert_type (table, row, col, &(uu.f), ac_type_float) ;
		       if (debug) fprintf (stderr, "%g", uu.f);
		       break ;
		     case _LongInt:
		     case _LongFloat:
		     case _Text:
		       if (uu.s)
			 { 
			   if (test) return TRUE ;
			   ac_table_insert_text (table, row, col, uu.s) ;
			   if (debug) fprintf (stderr, "%s", uu.s) ;
			 }
		       break ;
		     case __DNA1:
		       if (var->dclNode->dnaStack)
			 { 
			   char *cp = stackText (var->dclNode->dnaStack, 0) ;
			   if (test) return TRUE ;
			   ac_table_insert_text (table, row, col, cp) ;
			   if (debug) fprintf (stderr, "%s", uu.s) ;
			 }
		       break ;
		     case __Protein1:
		       if (var->dclNode->pepStack)
			 { 
			   char *cp = stackText (var->dclNode->pepStack, 0) ;
			   if (test) return TRUE ;
			   ac_table_insert_text (table, row, col, cp) ;
			   if (debug) fprintf (stderr, "%s", uu.s) ;
			 }
		       break ;
		     case _DateType:
		       {
			 char buf25[25] ;
			 uu.time = (mytime_t)var->dclNode->z ;
			 if (test) return TRUE ;
			 ac_table_insert_type (table, row, col, &(uu.time), ac_type_date) ;
			 timeShow (uu.time, buf25, 25) ;
			 if (debug) fprintf (stderr, "%s", buf25) ;
		       }
		       break ;
		     default:
		       break ;
		     }
		 }
	       if (var->dclNode->key && var->dclNode->key > _LastC)
		 {
		   KEY key = var->dclNode->key ;
		   int i = keySetMax (bql->ksOut) ;
		   KEY old = i ? keySet (bql->ksOut, i-1) : 0 ;
		   if (bql->ksOut && isFirst && key && class (key) && key != old)
		     keySet (bql->ksOut, i) = key ;
		 }
	       if (var->dclNode->key)
		 {
		   if (test) return TRUE ;
		   ac_table_insert_type (table, row, col, &(var->dclNode->key), ac_type_key) ;
		 }
	       if (debug) fprintf (stderr, "%s", name(var->dclNode->key)) ;
	     }
	   else 
	     {
	       if (debug) fprintf (stderr, "\tNULL") ;
	     }
	   if (coma)
	     var = coma = coma->right ;
	   else
	     var = 0 ;
	   isFirst = FALSE ;
	}
    }
  if (debug) fprintf (stderr, "\n") ;

  return ok ;
} /* bsDoExportLine */

static BOOL bqlExportLine (BQL *bql)
{
  if (bqlDoExportLine (bql, TRUE))
    return bqlDoExportLine (bql, FALSE) ;
  return FALSE ;
} /* bqlExportLine */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

static BOOL bqlExpandNext (BQL *bql, NODE *coma)
{
  NODE *right = coma ? coma->right : 0 ;
  BOOL ok = FALSE ;

  if (bql->where) return TRUE ;
  if (right)
    ok = bqlExpandComa (bql, right) ;
  else
    {
      if (bql->counting)
	{
	  ok = TRUE ;
	  if(bql->inCurly) bql->counting++ ;
	}
      else if (bql->mayBe)
	return TRUE ;
      else
	ok = bqlExportLine (bql) ;
    } 
  return ok ;
} /* bqlExpandNext  */

/**************************************************************/

static BOOL bqlExpandClass (BQL *bql, NODE *node, NODE *coma)
{
  BOOL ok = FALSE, ok2 = FALSE ;
  NODE *var = node->down ;
  NODE *right = node->right ;
  int i = 0 ;
  KEYSET ks = 0 ;
  NODE *where = node->up ;
  KEY subClasse = right->var ;
  KEY classe = right->key ;
  KEY classFilter = right->uType ;

  if (where && where->type != WHERE && where->type != WHERE_IF)
    where = 0 ;
  
  for (i = 1 ; i < lexMax (classe) ; i++)
    {
      if (bql->maxLine && bql->results->rows > bql->maxLine)
	break ;
      ok = TRUE ;
      if (var->myObj)
	{
	  bsDestroy (var->obj) ;
	  var->myObj = FALSE ;
	}
      var->obj = 0 ;
      var->parent = 0 ;
      var->key =  KEYMAKE (classe, i) ;
      if (!lexIsKeyVisible(var->key) ||
	  (classFilter && ! lexIsInClass(var->key, subClasse))
	  )
	continue ;
      if (where)
	{
	  bql->where = TRUE ;
	  ok = bqlExpandWhere (bql, where, coma) ;
	  bql->where = FALSE ;
	}
      if (ok)
	{
	  bql->openObj = TRUE ;
	  bql->mayBe = FALSE ;
	  if (0)
	    {
	      bql->openObj = FALSE ;
	      ok = bqlExpandNext (bql, coma) ; 
	      bql->openObj = TRUE ;
	      if (ok && bql->mayBe) /* may be */
		ok = bqlExpandNext (bql, coma) ; 
	    }
	  else
	    ok = bqlExpandNext (bql, coma) ; 
	}
      ok2 |= ok ;
      bqlClean (bql, node) ;
    }
  keySetDestroy (ks) ;
  return ok2 ;
} /* bqlExpandClass */

/*************************************************************************************/

static BOOL bqlNumberCompare (BQLTYPE type, double z1, double z2)
{
  BOOL ok = FALSE ;
  float z1f = z1, z2f = z2 ;
  z1 = z1f ; z2 = z2f ;

  switch (type)
    {
    case EQ:
      ok = (z1 * 1.000000001 >=  z2 && z2 * 1.000000001 >= z1)  ? TRUE : FALSE ;
      break ;
    case NEQ:
      ok = (z1 * 1.000000001 >=  z2 && z2 * 1.000000001 >= z1) ? FALSE : TRUE ;
      break ;
    case LEQ:
    case LLK:
      ok = (z1 >= z2 ? TRUE : FALSE) ;
      break ;
    case SEQ:
    case SLK:
      ok = (z1 <= z2 ? TRUE : FALSE) ;
      break ;
    case LT:
      ok = (z1 > z2 ? TRUE : FALSE) ;
      break ;
    case ST:
      ok = (z1 < z2 ? TRUE : FALSE) ;
      break ;
    default:
      break ;
    }

  return ok ;
} /*  bqlNumberCompare */

/*************************************************************************************/
/* 
 * In < <= == >= > we use number comparison, i.e. round up to most precise  
 * In <~ ~ >~ we round up to less precise so that  
 *          `2012` ~  `2012-3` ---> true 
 *          `2012` == `2012-3` ---> false
 *          `2012` <~ `2012-3` ---> true 
 *          `2012` <= `2012-3` ---> true 
 *          `2012` >~ `2012-3` ---> true 
 *          `2012` >= `2012-3` ---> false
 */
static BOOL bqlDateCompare (BQLTYPE type, double z1, double z2)
{
  mytime_t t1, t2 ;

  if (z1 < 0 || z2 < 0)
    return FALSE ;
  z1 += .000000001 ;
  z2 += .000000001 ;
  t1 = (mytime_t) z1 ;
  t2 = (mytime_t) z2 ;
  switch (type)
    {
    case EQ:
      if (t1 == t2)
       return TRUE ;
     break ;
    case NEQ:
     if (timeComparison (0, t1, t2))
       return FALSE ;
     return TRUE ;
      break ;
    case LT:
      if (t1 > t2)
	return TRUE ;  
      break ;
    case ST:
      if (t1 < t2)
	return TRUE ;
      break ;
    case LEQ:
      if (t1 >= t2)
	return TRUE ;
      break ;
    case SEQ:
      if (t1 <= t2)
	return TRUE ;
      if (z1 == z2)
	return TRUE ;
      break ;
    case LIKE:
    case LIKEREGEXP:
     if (t1 == t2 || timeComparison (0, t1, t2))
       return TRUE ;
     break ;
    case LLK:
      if (t1 >= t2 || timeComparison (0, t1, t2) || timeComparison (+1, t1, t2))
	return TRUE ;  
      break ;
    case SLK:
      if (t1 <= t2 || timeComparison (0, t1, t2) || timeComparison (-1, t1, t2))
	return TRUE ;
      break ;
    default:
      break ;
    }

  return FALSE ;
} /*  bqlDateCompare */

/*************************************************************************************/
/* store the position in var->row, call with 
 *     pass == 0 : new           ...   set row = 0
 *     pass == 1 : do not move   ...   row unchanged
 *     pass == 2 : iterate       ...   row++  (incremented before returning)
 *
 * return 
 *     0 : no value in given row
 *     1 : a number was returned in *zp
 *     2 : a date was returned in *zp
 *     the returned value is a valid char*, the printable value, for example a key name
 */
static const char *bqlExpandEquation (BQL *bql, NODE *node, int pass, double *zp)
{
  const char *ccp = 0 ;
  NODE *left ;
  NODE *right ;

  if (! node)
    return 0 ;
  left = node->down ; right = node->right  ;
  
  switch (node->type)
    {
    case NOT:
      *zp = 0 ;
      if (! right)
	return 0 ;
      ccp = bqlExpandEquation (bql, right, pass, zp) ;
      if (ccp == assVoid (1))  /* a number */
	{
	  if (*zp == 0) /* FALSE */ 
	    *zp = 1 ;   /* return TRUE */
	  else       
	    *zp = 0 ;   /* TRUE */
	  *zp = 1 ;     /* return FALSE */
	  return ccp ;  /* return a number */
	}
      else if (ccp)   /* TRUE, a text or the name of a variable */
	{
	   *zp = 0 ;   /* return FALSE */
	   ccp = assVoid (1) ; /* a number */
	   return ccp ;  /* return a number */
	}
      else /* !ccp || ! *ccp : i.e.  FALSE */
	{
	  *zp = 1 ;   /* return TRUE */
	  ccp = assVoid (1) ; /* a number */
	  return ccp ;  /* return a number */
	}
      break ;
    case ROUND:
    case CURLY:
    case SET:
      if (pass == 2)
	return 0 ;
      if (left)
	return bqlExpandEquation (bql, left, pass, zp) ;
      return 0 ;
      break ;
    case NUMBER:
      if (pass == 2)
	return 0 ;
      if (zp) *zp = node->z ;
      ccp = assVoid (1) ;
      return ccp ;
      break ;
    case DATE:
      if (pass == 2)
	return 0 ;
      if (zp) *zp = node->z ;
      ccp = assVoid (2) ;
      return ccp ;
      break ;
    case DOLLAR:
    case QUOTE:
    case DOUBLEQUOTE:
      if (pass == 2)
	return 0 ;
      ccp = stackText (bql->s, node->mark) ; 
      return ccp ;
      break ;
    case VAR:
      if (pass == 2)
	return 0 ;
      node = node->dclNode ;
      if (node && node->isNumber)
	{
	  if (zp) *zp = node->z ;
	  ccp = assVoid (1) ;
	} 
      else if (node && node->isDate)
	{
	  if (zp) *zp = node->z ;
	  ccp = assVoid (2) ;
	} 
      else if (node && node->dnaStack)
	{
	  ccp = stackText (node->dnaStack, 0) ;
	}
       else if (node && node->pepStack)
	{
	  ccp = stackText (node->pepStack, 0) ;
	}
      else if (node && node->type == VAR && ! node->key)
	ccp = node->dclNode->uu.s ; /* 2020_10_06, was ccp = 0 */
      else if (node)
	{     
	  if (node->dclNode->uType == _Text)
	    ccp = node->dclNode->uu.s ;
	  else if (node->dclNode->key)
	    ccp = name (node->dclNode->key) ;
	  else if (node->dclNode->uType >= _LastC)
	    ccp = node->dclNode->uu.k ? name (node->dclNode->uu.k) : (bql->counting ? 0 : (char *)1) ;
	} 
      return ccp ;
      break ;

    case COUNT:
      {
	int counting = bql->counting ; /* store */
	if (pass == 2)
	  return 0 ;
	NODE *new = 0 ;
	if (! node->down)
	  {
	    new = (NODE *) halloc (sizeof (NODE), bql->h) ;
	    new->up = node ;
	    node->down = new ;
	    new->nCol = 0 ;
	    new->dclNode = new ;
	  }
	bql->counting = 1 ; /* reinitialize */
	if (node->right)
	  {
	    if (node->right->type == CURLY)
	      {
		bqlExpandCurly (bql, node->right) ;
	      }
	    else if (node->right->type == FROM)
	      bqlExpandFrom (bql, node, 0) ;
	    else
	      bqlExpandIn (bql, node, 0) ;
	  }
	if (zp) *zp = bql->counting - 1 ;
	bql->counting = counting ; /* restore */
	
	if (new)
	  {
	    node->down = 0 ;
	    ac_free (new) ;
	  } 
	ccp = assVoid (1) ;  /* return a number */
	return ccp ;
      }
      break ; 

    case MINUS:
      if (! left)
	{
	  double z2 = 0 ;
	  ccp = bqlExpandEquation (bql, right, pass, &z2) ;
	  if (ccp == assVoid (1))
	    {
	      if (zp) *zp = - z2 ;
	      return ccp ;
	    }
	  return 0 ;
	}
      else /* fall thru */
	{} ;
    case PLUS:
    case MULT:
    case DIVIDE:
    case POWER:
    case MODULO:
      {
	const char *ccp1, *ccp2 ;
	double z1 = 0, z2 = 0 ;
	int pass1, pass2 ;

	if (! zp)
	  return 0 ;
	switch (pass)
	  {
	  case 0: /* init */
	    pass1 = pass2 = 0 ;
	    break ;
	  case 1: /* stay on same rows */
	    pass1 = pass2 = 1 ;
	    break ;
	  case 2:  /* iterate right, if impossible iterate left */
	    pass1 = 1 ; pass2 = 2 ;
	    break ;
	  }

      lao:
	ccp1 = bqlExpandEquation (bql, left, pass1, &z1) ;
	if (ccp1 != assVoid (1) &&  /* not a number */
	    ccp1 != assVoid (2)     /* not a date */
	    )
	  return 0 ;
	ccp2 = bqlExpandEquation (bql, right, pass2, &z2) ;
	if (ccp1 == assVoid (2) || ccp2 == assVoid (2)) /* a date */
	  {
	    if (ccp1 == assVoid (2) && ccp2 == assVoid (2)) /* a pair of dates */
	      {
		if (node->type == MINUS && z1 >= 0 && z2 >= 0)
		  {
		    int dt = 0 ;
		    timeDiffSecs (z1, z2, &dt) ;
		    *zp = dt ;
		    return assVoid (2) ;
		  }
	      }
	    if (ccp1 == assVoid (2) && ccp2 == assVoid (1)) /* z1 is a date */
	      {
		if (node->type == PLUS && z1 >= 0)
		  {
		    *zp = z1 ; *zp += z2 ;
		    return assVoid (1) ;
		  }
		if (node->type == MINUS && z1 >= 0)
		  {
		    *zp = z1 ; *zp -= z2 ;
		    if (*zp > 0)
		      return assVoid (1) ;
		  }
	      }
	    if (ccp1 == assVoid (1) && ccp2 == assVoid (2)) /* z2 is a date */
	      {
		if (node->type == PLUS && z2 >= 0)
		  {
		    *zp = z1 ; *zp += z2 ;
		    return assVoid (1) ;
		  }
		if (node->type == MINUS && z2 >= 0)
		  {
		    *zp = z1 ; *zp -= z2 ;
		    if (*zp > 0)
		      return assVoid (1) ;
		  }
	      }


	    if (pass == 2) /* try to iterate */
	      {
		pass1 = 2 ; pass2 = 0 ; goto lao ;
	      }
	    else
	      return 0 ;
	  }
	if (ccp2 != assVoid (1)) /* not a number */
	  {
	    if (pass == 2) /* try to iterate */
	      {
		pass1 = 2 ; pass2 = 0 ; goto lao ;
	      }
	    else
	      return 0 ;
	  }

	if (1)
	  {
	    int k1, k2 ;
	    switch (node->type)
	      {
	      case PLUS:    *zp = z1 + z2 ; break ;
	      case MINUS:    *zp = z1 - z2 ; break ;
	      case MULT:    *zp = z1 * z2 ; break ;
	      case DIVIDE:  if(z2) *zp = z1 / z2 ; else ccp2 = 0 ; break ;
	      case POWER: 
		k2 = z2 ;
		if (z2 == 0) *zp = 1 ;
		else if (z1 == 1) *zp = 1 ;
		else if (k2 == z2) 
		  {
		    int i ;
		    k1 = 1 ;
		    if (k2 < 0) { k2 = -k2 ; k1 = -1 ; }
		    for (*zp = z1, i = 1 ; i < k2 ; i++)
		      *zp *= z1 ;
		  
		    if (k1 == -1)
		      *zp = 1.0/(*zp) ;
		  }
		else if (z1 < 0)  return 0 ; /* ERREUR */
		else 
		  {
		    *zp = log (z1) * z2 ; 
		    *zp = exp (*zp) ;
		  }
		break ;			 
	      case MODULO:  
		k1 = z1 + .5 ; k2 = z2 + .5 ; 
		if (k2 < 0)
		  ccp2 = 0 ;
		else if (k2 == 0)
		  *zp = 0 ;
		else
		  {
		    if (k1 <= 0)
		      {
			int k3 = -k1/k2 ;
			k1 = z1 + .5 + (k3 +1) * k2 ;
		      }
		    *zp = k1 % k2 ;
		  }
		break ;  
	      default:
		break ;
	      }
	    return ccp2 ;
	  }
      }

    case TAG:
      {
	/* on ouvre l'objet et on scan jusqu'a la ligne var->row */
	BOOL ok = FALSE ;
	NODE *var = node->down ;
	NODE *tagNode = node->right ;
	KEY key = var && var->dclNode ? var->dclNode->key : 0 ;
	KEY tag = tagNode ? tagNode->key : 0 ;
	OBJ obj = 0 ;
	int row ;
	int rightOf = 0 ;
	
	if (! tag)
	  return 0 ;
	
	if (tag)
	  {
	    obj = var->dclNode->obj ;
	    if (obj && bsFindTag (obj, tag))
	      {
		ok = TRUE ; 
	      }
	    else if (keyFindTag (key, tag))
	      {
		if (! bql->openObj)
		  {
		    bql->mayBe = TRUE ;
		  }
		else
		  {
		    var->dclNode->myObj = TRUE ;
		    obj = var->dclNode->obj = bsCreate (key) ;
		    ok = TRUE ;
		  }
	      }
	    else
	      ok = FALSE ;
	  } 
	if (ok && obj)
	  {
	    int nCol = var->dclNode->nCol ;
	    Array aa ;
	    
	    if (var->aa && var->key != key)
	      ac_free (var->aa) ; 
	    aa = var->aa ;
	    var->col = 0 ; var->uu.k = 0 ;
	    var->key = key ;
	    var->col = 0 ;
	    if (pass == 0) 
	      var->row = 0 ;
	    if (pass == 2)
	      var->row++ ;
	    if (! aa)
	      {
		if (bsFindTag (obj, tag))
		  {
		    if (var->dclNode->nCol)
		      {
			var->bsMark = bsHandleMark (obj, var->bsMark, bql->h) ; 
			var->timeStamp = bsGetTimeStamp (obj) ;
			var->myBsmark = TRUE ; 
		      }
		    var->parent = var->dclNode ;
		  }
		else
		  return 0 ;
		ok = FALSE ;      
		var->uType = bsType (obj, _bsRight) ;
		
		aa = var->aa = arrayHandleCreate (32, BSunit, 0) ;
		var->myAa = TRUE ;
		bsGetArray (obj, tag, var->aa, nCol) ;
	      }
	    
	    if (!aa)
	      return 0 ;
	    
	    row = var->row ;
	    if (aa && row * nCol < arrayMax (aa))
	      {
		BSunit *up = arrp (aa, row * nCol + rightOf, BSunit) ;
		var->uu = up[0] ;
		if (var->uType == _Int)
		  {
		    var->isNumber = TRUE ;
		    var->isDate = FALSE ;
		    var->z = var->uu.i ;
		    ccp = assVoid(1) ; 
		    if (zp) *zp = var->z ;
		  }
		else if (var->uType == _Float)
		  {
		    var->isNumber = TRUE ;
		    var->isDate = FALSE ;
		    var->z = var->uu.f ; 
		    ccp = assVoid(1) ;  
		    if (zp) *zp = var->z ;
		  }
		else if (var->uType == _LongInt)
		  {
		    long long int lli = 0 ;
		    var->isNumber = TRUE ;
		    var->isDate = FALSE ;
		    var->z = 0 ;
		    if (var->uu.s) sscanf (var->uu.s, "%lli", &(lli)) ;
		    var->z = lli ;
		    ccp = assVoid(1) ;  
		    if (zp) *zp = var->z ;
		  }
		else if (var->uType == _LongFloat)
		  {
		    var->isNumber = TRUE ;
		    var->isDate = FALSE ;
		    var->z = 0 ;
		    if (var->uu.s) sscanf (var->uu.s, "%lf", &(var->z)) ; 
		    ccp = assVoid(1) ;  
		    if (zp) *zp = var->z ;
		  }
		else if (var->uType == _DateType)
		  {
		    var->isNumber = FALSE ;
		    var->isDate = TRUE ;
		    var->z = var->uu.time ;
		    ccp = assVoid(2) ;  
		    if (zp) *zp = var->uu.time ;
		  }
		else if (var->uType > _LastC)
		  {
		    ccp = name (var->uu.k) ;
		  }
		else
		  var->isNumber = var->isDate = FALSE ;
	      }
	  }
	return ccp ;
      }
    case RIGHTOF:
      {
	/* on ouvre l'objet et on scan jusqu'a la ligne var->row */
	
      } 
      break ;
    default:
      break ;
    }

  return ccp ;
} /* bqlExpandEquation */

/*************************************************************************************/

static BOOL bqlExpandObject (BQL *bql, NODE *node, NODE *coma)
{
  BOOL ok = FALSE, ok2 = FALSE ;
  NODE *var = node->down ;
  NODE *obj = node->right ;
  NODE *objComa = obj ? obj->right : 0 ;
  NODE *down3 = objComa ? objComa->right : 0 ;
  NODE *objClass = objComa ? objComa->down : 0 ;

  int i = 0 ;
  KEYSET ks = keySetCreate () ;
  NODE *where = node->up ;

  if (where && where->type != WHERE && where->type != WHERE_IF)
    where = 0 ;
  
  if (objClass && down3)
    {
      KEY subClasse = 0 ;
      KEY classe = 0 ;
      unsigned char classFilter ;
      const char * clNam = 0 ;
      
      if (objClass->dclNode)
	clNam = bqlExpandEquation (bql, objClass, 0, 0) ;
      else if ( objClass->var)
	clNam = dictName (bql->dict, objClass->var) ;
      else if ( objClass->mark)
	clNam = stackText (bql->s, objClass->mark) ;

      if (clNam && lexword2key(clNam, &subClasse, _VClass)) 
	{ 
	  const char *objNam = 0 ;
	  KEY key ;
	  double z = 0 ;
	  int n = 0, ii ;
	  BOOL caseSensitive ;
	  char buff[128] ;
	  
	  if (down3->dclNode)
	    objNam = bqlExpandEquation (bql, down3, 0, &z) ;
	  else if (down3->var)
	     objNam = dictName (bql->dict, down3->var) ;
	  else if (down3->mark)
	     objNam = stackText (bql->s, down3->mark) ;

	  classe = goodClass (subClasse, &classFilter) ;
	  key = KEYMAKE (classe, 1) ;
	  caseSensitive = pickCaseSensitive (key) ;

	  if (objNam == assVoid (1)) /* a number */
	    {
	      sprintf (buff, "%g", z) ;
	      objNam =  buff ;
	    }
	  else if (objNam == assVoid (2)) /* a date */
	    {
	      timeShow (z , buff, 127) ;
	      objNam =  buff ;
	    }
	  if (objNam)
	    {
	      if (lexword2key (objNam, &key, classe))
		{
		  if (lexIsKeyVisible(key) &&
		      (! classFilter ||  lexIsInClass(key, subClasse))
		      )
		    keySet (ks, n++) = key ;
		}
	      else
		{	      
		  for (ii = 1 ; ii < lexMax (classe) ; ii++)
		    {
		      key = KEYMAKE (classe, ii) ;
		      if (pickMatchCaseSensitive (name (key), objNam, caseSensitive) &&
			  lexIsKeyVisible(key))
			keySet (ks, n++) = key ;
		    }
		}
	    }
	  keySetSort (ks) ;
	  keySetCompress (ks) ;
	}
    }
  for (i = 0 ; ks && i < keySetMax (ks) ; i++)
    {
      ok = TRUE ;
      if (var->myObj)
	{
	  bsDestroy (var->obj) ;
	  var->myObj = FALSE ;
	}
      var->obj = 0 ;
      var->parent = 0 ;
      var->key = keySet (ks, i) ;
      if (!lexIsKeyVisible(var->key))
	continue ;
         
      if (where)
	{
	  bql->where = TRUE ;
	  ok = bqlExpandWhere (bql, where, coma) ;
	  bql->where = FALSE ;
	}
      if (ok)
	{
	  bql->openObj = TRUE ;
	  bql->mayBe = FALSE ;
	  if (0)
	    {
	      bql->openObj = FALSE ;
	      ok = bqlExpandNext (bql, coma) ; 
	      bql->openObj = TRUE ;
	      if (ok && bql->mayBe) /* may be */
		ok = bqlExpandNext (bql, coma) ; 
	    }
	  else
	    ok = bqlExpandNext (bql, coma) ; 
	}
      ok2 |= ok ;
      if (coma && coma->right)
	bqlClean (bql, coma->right) ;
      bqlClean (bql, var) ;
    }
  keySetDestroy (ks) ;
  return ok2 ;
} /* bqlExpandObject */

/*************************************************************************************/
/*************************************************************************************/
/* syntax
   select s,pep from s in class sequence, pep in PEPTIDE(s)
   select s,pep from s in class sequence, pep in PEPTIDE(s,x, x+7) or any equation
*/


static BOOL bqlExpandPeptide (BQL *bql, NODE *node, NODE *coma, BOOL isCoding)
{
  BOOL ok = FALSE, ok2 ;
  NODE *var = node->down ;
  NODE *where = node->up ;
  NODE *dnaNode = node->right ;
  NODE *dnaVarComa = dnaNode ? dnaNode->right : 0 ;
  NODE *dnaVar = dnaVarComa ? dnaVarComa->down : 0 ;
  NODE *dna1Coma = dnaVarComa ? dnaVarComa->right : 0 ;
  NODE *dna1Node = dna1Coma ? dna1Coma->down : 0 ;
  NODE *dna2Coma = dna1Coma ? dna1Coma->right : 0 ;
  NODE *dna2Node = dna2Coma ? dna2Coma->down : 0 ;
  int dna1 = -99999999, dna2 =  -99999999 ;

  if (where && where->type != WHERE && where->type != WHERE_IF)
    where = 0 ;

  if (! dnaVar && dnaVarComa->type == VAR)
    { dnaVar = dnaVarComa ; dna1Coma = 0 ; }

  ok = TRUE ;
  /* if we need coordinates and cannot find them, no need to load the DNA */
  if (ok && dna1Coma && dna1Coma->type != COMA)
    { dna1Node = dna1Coma ; dna2Node = dna2Coma = 0 ; }
  if (ok && dna1Node)
    {
      double xx = 0 ;
      if (bqlExpandEquation (bql, dna1Node, 0, &xx) == assVoid(1))
	dna1 = dna1Node->z = xx ;
      else
	ok = FALSE ;
    }
  if (ok && dna2Coma && dna2Coma->type != COMA)
    dna2Node = dna2Coma ;
  if (ok && dna2Node)
    {
      double xx = 0 ;
      if (bqlExpandEquation (bql, dna2Node, 0, &xx) == assVoid (1))
	dna2 = dna2Node->z = xx ;
      else
	ok = FALSE ;
    }
  if (ok && dnaVar->type == VAR && dnaVar->dclNode && dnaVar->dclNode->key)
    {
      var->dnaD = peptideTranslate (dnaVar->dclNode->key, isCoding) ;
    }
  if (ok && var->dnaD)
    {
      int d1, d2 ;
      BOOL reverse = FALSE ;
      Array dna = var->dnaD ;

      if (dna1 ==  -99999999)
	dna1 = 1 ;
      if (dna2 ==  -99999999)
	dna2 = arrayMax (dna) ;

      if (dna1 <= dna2)
	{ d1 = dna1 ; d2 = dna2 ; reverse = FALSE ; }
      else
	{ d1 = dna2 ; d2 = dna1 ; reverse = TRUE ; }
      
      if ( 
	  (
	   d1 > 0 && d1 >  arrayMax (dna) && /* attention because arrayMax is unsigned int */
	   d2 > 0 && d2 >  arrayMax (dna)
	   ) ||
	  ( d1 < 1 && d2 < 1)
	   )
	;
      else
	{	
	  if (d1 > 0 && d1 >  arrayMax (dna)) /* attention because arrayMax is unsigned int */
	    d1 = arrayMax (dna) ;
	  if (d2 > 0 && d2 >  arrayMax (dna))
	    d2 = arrayMax (dna) ;
	  if (d1 < 1) d1 = 1 ;
	  if (d2 < 1) d2 = 1 ;
	  if (d1 <= d2 && d1 > 0 && !reverse)  /* we cannot reverse a peptide */
	    {
	      int i ;
	      unsigned char *cp, *cq ;
	      Array pepPiece ;
	      
	      dna1 = d1 ;
	      pepPiece = arrayCreate (d2 - d1 + 2, char) ;
	      array (pepPiece, d2 - d1 + 1, unsigned char) = 0 ; /* zero terminate */
	      arrayMax (pepPiece) = d2 - d1 + 1 ; /* restore */
	      for (i = 0, cp = arrp (pepPiece, 0, unsigned char), cq = arrp (dna, dna1 - 1, unsigned char) ;
		   i < d2 - d1 + 1 ; cp++, cq++, i++)
		{
		  *cp = pepDecodeChar [ (int)*cq] ;
		  if (isCoding && *cp == '*')
		    { cp[1] = 0 ; break ;
		    }
		}
	      var->pepStack = stackCreate (0) ;
	      pushText (var->pepStack, arrp (pepPiece, 0, char)) ;
	      arrayDestroy (pepPiece) ;
	      ok = TRUE ;
	    }
	}
      arrayDestroy (var->dnaD) ;
      var->uType = __Protein1 ;
      var->dclNode->uType = __Protein1 ;
    }

  ok = TRUE ;
  if (where)
    {
      bql->where = TRUE ;
      ok = bqlExpandWhere (bql, where, coma) ;
      bql->where = FALSE ;
    }
  if (ok)
    ok = bqlExpandNext (bql, coma) ; 
  ok2 = ok ;

  bqlClean (bql, var) ;

  return ok2 ;
} /* bqlExpandPeptide */

/*************************************************************************************/
/* syntax
   select s,dna from s in class sequence, dna in DNA(s)
   select s,dna from s in class sequence, dna in DNA(s,x, x+7) or any equation
*/

static BOOL bqlExpandDNA (BQL *bql, NODE *node, NODE *coma)
{
  BOOL ok = FALSE, ok2 ;
  NODE *var = node->down ;
  NODE *where = node->up ;
  NODE *dnaNode = node->right ;
  NODE *dnaVarComa = dnaNode ? dnaNode->right : 0 ;
  NODE *dnaVar = dnaVarComa ? dnaVarComa->down : 0 ;
  NODE *dna1Coma = dnaVarComa ? dnaVarComa->right : 0 ;
  NODE *dna1Node = dna1Coma ? dna1Coma->down : 0 ;
  NODE *dna2Coma = dna1Coma ? dna1Coma->right : 0 ;
  NODE *dna2Node = dna2Coma ? dna2Coma->down : 0 ;
  int dna1 = -99999999, dna2 =  -99999999 ;

  if (where && where->type != WHERE && where->type != WHERE_IF)
    where = 0 ;

  if (! dnaVar && dnaVarComa && dnaVarComa->type == VAR)
    { dnaVar = dnaVarComa ; dna1Coma = 0 ; }
  
  ok = TRUE ;
  /* if we need coordinates and cannot find them, no need to load the DNA */
  if (ok && dna1Coma && dna1Coma->type != COMA)
    { dna1Node = dna1Coma ; dna2Node = dna2Coma = 0 ; }
  if (ok && dna1Node)
    {
      double xx = 0 ;
      if (bqlExpandEquation (bql, dna1Node, 0, &xx) == assVoid(1))
	dna1 = dna1Node->z = xx ;
      else
	ok = FALSE ;
    }
  if (ok && dna2Coma && dna2Coma->type != COMA)
    dna2Node = dna2Coma ;
  if (ok && dna2Node)
    {
      double xx = 0 ;
      if (bqlExpandEquation (bql, dna2Node, 0, &xx) == assVoid(1))
	dna2 = dna2Node->z = xx ;
      else
	ok = FALSE ;
    }
  if (ok && dnaVar && dnaVar->type == VAR && dnaVar->dclNode && dnaVar->dclNode->key)
    {
      if (! var->dnaD || var->dnaKey != dnaVar->dclNode->key) 
	{
	  arrayDestroy (var->dnaD) ;
	  var->dnaD = dnaHandleGet (dnaVar->dclNode->key, bql->h) ;
	  var->dnaKey = dnaVar->dclNode->key ;
	}
    }
  else
    ok = FALSE ;
  if (ok && var && var->dnaD)
    {
      int d1, d2 ;
      BOOL reverse = FALSE ;
      Array dna = var->dnaD ;
 
      if (dna1 ==  -99999999)
	dna1 = 1 ;
      if (dna2 ==  -99999999)
	dna2 = arrayMax (dna) ;

      if (dna1 <= dna2)
	{ d1 = dna1 ; d2 = dna2 ; reverse = FALSE ; }
      else
	{ d1 = dna2 ; d2 = dna1 ; reverse = TRUE ; }
      if (
	  (
	   d1 > 0 && d1 >  arrayMax (dna) && /* attention because arrayMax is unsigned int */
	   d2 > 0 && d2 >  arrayMax (dna)
	   ) ||
	  ( d1 < 1 && d2 < 1)
	  )
	;
      else
	{	
	  if (d1 > 0 && d1 >  arrayMax (dna)) /* attention because arrayMax is unsigned int */
	    d1 = arrayMax (dna) ;
	  if (d2 > 0 && d2 >  arrayMax (dna))
	    d2 = arrayMax (dna) ;
	  if (d1 < 1) d1 = 1 ;
	  if (d2 < 1) d2 = 1 ;
	  if (d1 <= d2 && d1 > 0)
	    {
	      int i ;
	      unsigned char *cp, *cq ;
	      Array dnaPiece ;
	      
	      dnaPiece = arrayCreate (d2 - d1 + 2, char) ;
	      array (dnaPiece, d2 - d1 + 1, unsigned char) = 0 ; /* zero terminate */
	      arrayMax (dnaPiece) = d2 - d1 + 1 ; /* restore */
	      if (! reverse)
		{
		  dna1 = d1 ;
		  for (i = 0, cp = arrp (dnaPiece, 0, unsigned char), cq = arrp (dna, dna1 - 1, unsigned char) ;
		       i < d2 - d1 + 1 ; cp++, cq++, i++)
		    *cp = dnaDecodeChar [ (int)*cq] ;
		}
	      else
		{
		  dna1 = arrayMax (dna) - d2 + 1 ; 
		  dna1 = d2 ;
		  dna = var->dnaD ;
		  for (i = 0, cp = arrp (dnaPiece, 0, unsigned char), cq = arrp (dna, dna1 - 1, unsigned char) ;
		       i < d2 - d1 + 1 ; cp++, cq--, i++)
		    *cp = dnaDecodeChar [(int)complementBase(*cq)] ;
		}
	      var->dnaStack = stackCreate (0) ;
	      pushText (var->dnaStack, arrp (dnaPiece, 0, char)) ;
	      arrayDestroy (dnaPiece) ;
	    }
	}
      var->uType = __DNA1 ;
      var->dclNode->uType = __DNA1 ;
    }

  ok = TRUE ;
  if (where)
    {
      bql->where = TRUE ;
      ok = bqlExpandWhere (bql, where, coma) ;
      bql->where = FALSE ;
    }
  if (ok)
    ok = bqlExpandNext (bql, coma) ; 
  ok2 = ok ;

  bqlClean (bql, var) ;

  return ok2 ;
} /* bqlExpandDNA */

/*************************************************************************************/
/* syntax
   select dd from dd in DATEDIFF(unit,date_1, date_2)
   unit is one of year, month, day, hour, minute, second
                  y, m, d, h, n, s  
                              n is for miNute as is SQL
select p1, p2, d1, d2 from a in class author where a == tom, p1 in a->papers, p2 in a->papers, d1 in p1->published, d2 in p2->published where DATEDIFF('year',d1,d2)> 0
*/
static BOOL bqlExpandDateDiff (BQL *bql, NODE *node, NODE *coma)
{
  BOOL ok = FALSE, ok2 ;
  NODE *var = node->down ;
  NODE *where = node->up ;
  NODE *dd = node->right ; /* DATEDIFF */
  NODE *ddComa = dd ? dd->right : 0 ;
  /*  NODE *unit = ddComa ? ddComa->down : 0 ;  $year */
  NODE *d1Coma = ddComa ? ddComa->right : 0 ;
  NODE *d1Node = d1Coma ? d1Coma->down : 0 ;   /* date1 */
  NODE *d2Node = d1Coma ? d1Coma->right : 0 ;  /* date2 */

  if (where && where->type != WHERE && where->type != WHERE_IF)
    where = 0 ;

  if (d2Node)
    ok = TRUE ;
  if (ok && ! d1Node->dclNode)
    ok = FALSE ;
  if (ok && ! d2Node->dclNode)
    ok = FALSE ;
  if (ok && ! d1Node->dclNode->isDate)
    ok = FALSE ;
  if (ok && ! d2Node->dclNode->isDate)
    ok = FALSE ;

  ok = TRUE ;
  if (ok)
    {
      var->uType = _Int ;
      var->dclNode->uType = _Int ;
      var->dclNode->isNumber = TRUE ;
      var->z = d1Node->dclNode->uu.time - d2Node->dclNode->uu.time ;
    }
  ok = TRUE ;
  if (where)
    {
      bql->where = TRUE ;
      ok = bqlExpandWhere (bql, where, coma) ;
      bql->where = FALSE ;
    }
  if (ok)
    ok = bqlExpandNext (bql, coma) ; 
  ok2 = ok ;

  bqlClean (bql, var) ;

  return ok2 ;
} /* bqlExpandDateDiff */

/*************************************************************************************/
/* syntax
   select a from a in @, n in a.name, t in a.timetamp
*/
static BOOL bqlExpandName (BQL *bql, NODE *node, NODE *coma)
{
  BOOL ok = FALSE, ok2 = FALSE ;
  NODE *where = node->up ;
  NODE *dd = node->right ; /* DATEDIFF */
  NODE *var = node->down ;
  NODE *dd1 = dd ? dd->down : 0 ;

  if (where && where->type != WHERE && where->type != WHERE_IF)
    where = 0 ;

  if (var &&  dd1->dclNode && dd1->dclNode->key)
    {
      ok = TRUE ;
      var->uType = _Text ;
      var->uu.s = name (dd1->dclNode->key) ;
    }

 if (! ok2)
    {
      if (!bql->counting) ok = TRUE ;
      if (where)
	{
	  bql->where = TRUE ;
	  ok = bqlExpandWhere (bql, where, coma) ;
	  bql->where = FALSE ;
	}
      if (ok)
	ok = bqlExpandNext (bql, coma) ; 
      ok2 |= ok ;
    }

  bqlClean (bql, var) ;

  return ok2 ;
} /* bqlExpandName */

/*************************************************************************************/
/* syntax
   select a from a in @, n in a.name, t in a.timetamp
*/
static BOOL bqlExpandClassName (BQL *bql, NODE *node, NODE *coma)
{
  BOOL ok = FALSE, ok2 = FALSE ;
  NODE *where = node->up ;
  NODE *dd = node->right ; /* DATEDIFF */
  NODE *var = node->down ;
  NODE *dd1 = dd ? dd->down : 0 ;

  if (where && where->type != WHERE && where->type != WHERE_IF)
    where = 0 ;

  if (var &&  dd1->dclNode && dd1->dclNode->key)
    {
      ok = TRUE ;
      var->uType = _Text ;
      var->uu.s = className (dd1->dclNode->key) ;
    }

 if (! ok2)
    {
      if (!bql->counting) ok = TRUE ;
      if (where)
	{
	  bql->where = TRUE ;
	  ok = bqlExpandWhere (bql, where, coma) ;
	  bql->where = FALSE ;
	}
      if (ok)
	ok = bqlExpandNext (bql, coma) ; 
      ok2 |= ok ;
    }

  bqlClean (bql, var) ;

  return ok2 ;
} /* bqlExpandClassName */

/*************************************************************************************/
/* syntax
   select a from a in @, n in a.name, t in a.timetamp
*/
static BOOL bqlExpandTimeStamp (BQL *bql, NODE *node, NODE *coma)
{
  BOOL ok = FALSE, ok2 = FALSE ;
  NODE *where = node->up ;
  NODE *dd = node->right ; /* DATEDIFF */
  NODE *var = node->down ;
  NODE *dd1 = dd ? dd->down : 0 ;

  if (where && where->type != WHERE && where->type != WHERE_IF)
    where = 0 ;

  if (var &&  dd1->dclNode && dd1->dclNode->timeStamp)
    {
      ok = TRUE ;
      var->uType = _Text ;
      var->uu.s = name (dd1->dclNode->timeStamp) ;
    }

 if (! ok2)
    {
      if (!bql->counting) ok = TRUE ;
      if (where)
	{
	  bql->where = TRUE ;
	  ok = bqlExpandWhere (bql, where, coma) ;
	  bql->where = FALSE ;
	}
      if (ok)
	ok = bqlExpandNext (bql, coma) ; 
      ok2 |= ok ;
    }

  bqlClean (bql, var) ;

  return ok2 ;
} /* bqlExpandTimeStamp */

/*************************************************************************************/

static BOOL bqlExpandTag (BQL *bql, NODE *node, NODE *coma)
{
  BOOL ok = FALSE, ok2 = FALSE ;
  NODE *var = node->down ;
  NODE *tagNode = node->right ;
  NODE *tagVar = tagNode ? tagNode->down : 0 ;
  NODE *tagTag = tagNode ? tagNode->right : 0 ;
  KEY key = 0 ;
  KEY tag = tagTag ? tagTag->key : 0 ;
  OBJ obj = 0 ;
  NODE *where = node->up ;
  int rightOf = 0 ;

  if (where && where->type != WHERE && where->type != WHERE_IF)
    where = 0 ;

  if (tagTag && tagTag->type == RIGHTOF && 
      tagTag->down && (tagTag->down->type == VAR || tagTag->down->type == TAGKEY)  && 
      tagTag->right && tagTag->right->type == NUMBER)
    {
      var->depth = tagTag->right->z >= 0 ?  tagTag->right->z  -1 : -1 ;
      tagTag = tagTag->down  ;
  }
  else
    var->depth = 0 ;

  if (! var || ! var->dclNode)
    return FALSE ;

  rightOf = var->col = var->depth ;
  if (var->dclNode->nCol < rightOf + 1 )
    var->dclNode->nCol = rightOf + 1 ;

  if (! tag && tagTag)
    tag = str2tag (stackText (bql->s,tagTag->mark)) ; 
  if (! tag)
    return FALSE ;
  if (! tagVar || ! tagVar->dclNode)
    return FALSE ;

  if (tagNode->type == HASTAG)
     {
       rightOf = var->col = var->depth  = -1 ;
       if (var->dclNode->nCol == 0)
	 {
	   if (tagVar->dclNode && tagVar->dclNode->key)
	     key = tagVar->dclNode->key ;
	   if (keyFindTag (key, tag))
	     {
	       var->key = var->uType = tag ; 
	       ok = TRUE ;
	     }
	   else
	     {
	       var->key = var->uType = 0 ;
	       ok = FALSE ;
	     }
	   if (where)
	     {
	       bql->where = TRUE ;
	       ok = bqlExpandWhere (bql, where, coma) ;
	       bql->where = FALSE ;
	     }
	   if (ok)
	     ok = bqlExpandNext (bql, coma) ; 
	   ok2  |= ok ;
	   return ok2 ;
	 }	   
     }
  if (tag)
    {
      if (tagVar->dclNode && tagVar->dclNode->key)
	key = tagVar->dclNode->key ;
      obj = tagVar->dclNode->obj ;
      if (obj && bsFindTag (obj, tag))
	{
	  ok = TRUE ; 
	}
      else if (keyFindTag (key, tag))
	{
	  if (! bql->openObj)
	    {
	      bql->mayBe = TRUE ;
	    }
	  else
	    {
	      tagVar->dclNode->myObj = TRUE ;
	      obj = tagVar->dclNode->obj = bsCreate (key) ;
	      ok = TRUE ;
	    }
	}
      else
	ok = FALSE ;
    }
  if (ok && obj)
    {
      int nCol = var->dclNode->nCol ;
      Array aa = var->aa ;

      if (! aa)
	{
	  aa = var->aa = arrayHandleCreate (32, BSunit, 0) ;
	  var->myAa = TRUE ;
	}
      arrayMax (aa) = 0 ;
      var->uu.k = 0 ;
      if (bsFindTag (obj, tag))
	{
	  if (var->dclNode->nCol >= 0)
	    {
	      var->bsMark = bsHandleMark (obj, var->bsMark, bql->h) ;  
	      var->timeStamp = bsGetTimeStamp (obj) ;
	      var->myBsmark = TRUE ; 
	    }
	  var->parent = tagVar->dclNode ;
	}
      else
	return FALSE ;

      ok = FALSE ;      
      if (tagNode->type == HASTAG || rightOf == -1)
	{ 
	  var->uType = tag ;
	  var->key = tag ;
	  var->row = 0 ;
	  var->col = var->depth  = -1 ;

	  if (var->dclNode->nCol > 0)
	    bsGetArray (obj, tag, var->aa, var->dclNode->nCol) ;
	  ok = TRUE ;
	  if (where)
	    {
	      bql->where = TRUE ;
	      ok = bqlExpandWhere (bql, where, coma) ;
	      bql->where = FALSE ;
	    }
	  if (ok)
	    ok = bqlExpandNext (bql, coma) ;  
	  ok2 |= ok ;
	  if (coma && coma->right)
	    bqlClean (bql, coma->right) ;
	}
      else if (bql->counting && (node->type == COUNT || node->up->type == COUNT))
	{ 
	  bsGetArray (obj, tag, var->aa, 1) ;
	  bql->counting = arrayMax (var->aa) + 1 ;
	  if (var->myBsmark)
	    ac_free (var->bsMark) ;
	  if (var->myAa)
	    arrayDestroy (var->aa) ;
	  var->myAa = var->myBsmark = FALSE ;
	  var->aa = 0 ;
	  return TRUE ;
	}
      else
	{
	  var->uType = bsType (obj, _bsRight) ;
	  var->key = 0 ;
	  if (rightOf > 0)
	    {
	      int i ;
	      for(i = 0 ; i < rightOf ; i++)
		{
		  if (bsFindTag (obj, _bsRight))
		    var->uType = bsType (obj, _bsRight) ;
		}
	    }
	 
	  if (bsGetArray (obj, tag, var->aa, var->dclNode->nCol))
	    {
	      int row ;
	      int iMax = arrayMax (aa) ;
	      BSunit *old = 0 ;
	      char timeBuf[25] ;

	      var->parent = tagVar->dclNode ;
	      if (tagNode->type == HASTAG)
		{
		  var->isNumber = FALSE ; /* bs->timeStamp */
		  var->isDate = FALSE ;
		  var->key = var->uType = tag ; 
		  var->row = 0 ; 
		  var->col = var->depth = -1 ;
		  ok = TRUE ;
		  if (where)
		    {
		      bql->where = TRUE ;
		      ok = bqlExpandWhere (bql, where, coma) ;
		      bql->where = FALSE ;
		    }
		  if (ok)
		    ok =bqlExpandNext (bql, coma) ; 
		  return ok ;
		}
	      else if (tagNode->type == MERGETAG)
		{
		  KEYSET ks= 0 ;
		  DICT *dict = 0 ;
		  if (rightOf > 0 && iMax > 2 * nCol)
		    {
		      if (var->uType == _Text)
			dict = dictCaseSensitiveHandleCreate (iMax/nCol + 1, 0) ;
		      else
			ks = keySetCreate () ;
		    }
		  var->isNumber = FALSE ;
		  var->isDate = FALSE ;
		  var->vtxt = vtxtHandleCreate (0) ;
		  nCol = var->dclNode->nCol ;
		  for (row = 0 ; row * nCol + rightOf < iMax ; row ++)
		    {
		      BSunit *up = arrp (aa, row * nCol + rightOf, BSunit) ;
		      if (old && old[0].k == up[0].k)
			continue ;
		      if (ks)
			{
			  if (keySetFind (ks,  up[0].k, 0))
			    continue ;
			  keySetInsert (ks, up[0].k) ;
			}
		      if (dict)
			{
			  if (! up[0].s || dictFind (dict,  up[0].s, 0))
			    continue ;
			  dictAdd (dict, up[0].s, 0) ;
			}
		      var->uu = up[0] ;
		      if (var->uType == _Int)
			vtxtPrintf (var->vtxt, "; %d", var->uu.i) ;
		      else if (var->uType == _Float)
			vtxtPrintf (var->vtxt, "; %g", var->uu.f) ;
		      else if (var->uType == _DateType)
			vtxtPrintf (var->vtxt, "; %s"
				    , timeShow (var->uu.time, timeBuf, 25)
				    ) ;
		      else if (var->uType > _LastC)
			{
			  if (var->uu.k)
			    vtxtPrintf (var->vtxt, "; %s", name(var->uu.k)) ;
			}
		      else
			{
			  if (var->uu.s)
			    vtxtPrintf (var->vtxt, "; %s", var->uu.s) ;
			}
		      old = up ;
		    }
		  ac_free (dict) ;
			keySetDestroy (ks) ;
		  ok = bqlExpandNext (bql, coma) ; 
		}
	      else
		{
		  nCol = var->dclNode->nCol ;
		  for (row = 0 ; row * nCol + rightOf < iMax ; row ++)
		    {
		      BSunit *up = arrp (aa, row * nCol + rightOf, BSunit) ;
		      double zzZ = 0 ;
		      int zzN = 0 ;
		      
		      var->uu = up[0] ;
		      if (var->uType == _Int)
			{
			  var->isNumber = TRUE ;
			  var->isDate = FALSE ;
			  zzZ = var->z = var->uu.i ;
			  zzN++ ;
			  var->key = 0 ;
			}
		      else if (var->uType == _Float)
			{
			  var->isNumber = TRUE ;
			  var->isDate = FALSE ;
			  zzZ = var->z = var->uu.f ;
			  zzN++ ;
			  var->key = 0 ;
			}
		      else if (var->uType == _LongFloat)
			{
			  var->isNumber = TRUE ;
			  var->isDate = FALSE ;
			  var->z = 0 ;
			  if (var->uu.s) sscanf (var->uu.s, "%lf", &(var->z)) ; 
			  zzN++ ;
			  zzZ = var->z ;
			  var->key = 0 ;
			  if (bql->minmaxavstd)
			    {
			      bql->minmaxavstdN++ ;
			      bql->minmaxavstdX += var->z ;
			      bql->minmaxavstdX2 += (var->z) * (var->z) ;
			    }
			}
		      else if (var->uType == _LongInt)
			{
			  var->isNumber = TRUE ;
			  var->isDate = FALSE ;
			  var->z = 0 ;
			  if (var->uu.s) 
			    {
			      long long int lli = 0 ;
			      sscanf (var->uu.s, "%lli", &(lli)) ;
			      zzZ = var->z = lli ;
			      zzN++ ;
			      if (bql->minmaxavstd)
				{
				  bql->minmaxavstdN++ ;
				  bql->minmaxavstdX += lli ;
				  bql->minmaxavstdX2 += lli * lli ;
				}
			    }
			  var->key = 0 ;
			}
		      else if (var->uType == _DateType)
			{
			  var->isNumber = FALSE ;
			  var->isDate = TRUE ;
			  var->z =  var->uu.time ;
			  var->key = 0 ;
			}
		      else if (var->uType > _LastC)
			{
			  if (var->key != var->uu.k)
			    {
			      if (var->myObj)
				{
				  bsDestroy (var->obj) ;
				  var->myObj = FALSE ;
				}
			      var->obj = 0 ;
			    }
 			  var->key = var->uu.k ;
			}
		      else
			var->isNumber = var->isDate = FALSE ;
		      if (zzN)
			switch (bql->minmaxavstd)
			  {
			  case MIN:
			    if (bql->minmaxavstdN == 0 || zzZ < bql->minmaxavstdX)
			      {
				bql->minmaxavstdN++ ;
				bql->minmaxavstdX = zzZ ;
			      }
			    break ;
			  case MAX:
			    if (bql->minmaxavstdN == 0 || zzZ > bql->minmaxavstdX)
			      {
				bql->minmaxavstdN++ ;
				bql->minmaxavstdX = zzZ ;
			      }
			    break ;
			  case STDDEV:
			    bql->minmaxavstdX2 += zzZ * zzZ ;
			  case AVERAGE:
			  case SUM:
			    bql->minmaxavstdN++ ;
			    bql->minmaxavstdX += zzZ ;
			    break ;
			  }
		      
		      var->dclNode->row = row ; 
		      if (!old  || old->k != up->k)
			{
			  ok = TRUE ;
			  if (where)
			    {
			      bql->where = TRUE ;
			      ok = bqlExpandWhere (bql, where, coma) ;
			      bql->where = FALSE ;
			    }
			  if (ok)
			    ok = bqlExpandNext (bql, coma) ; 
			  ok2 |= ok ;
			  if (coma && coma->right)
			    bqlClean (bql, coma->right) ;
			}
		      old = up ;
		    }
		}
	    }
	}
    }
  if (! ok2)
    {
      if (!bql->counting) ok = TRUE ;
      if (where)
	{
	  bql->where = TRUE ;
	  ok = bqlExpandWhere (bql, where, coma) ;
	  bql->where = FALSE ;
	}
      if (ok)
	ok = bqlExpandNext (bql, coma) ; 
      ok2 |= ok ;
    }
  bqlClean (bql, var) ;

  return ok2 ;
} /* bqlExpandTag */

/*************************************************************************************/
/* transitive closure */
static BOOL bqlExpandTTag (BQL *bql, NODE *node, NODE *coma)
{
  BOOL ok = FALSE, ok2 = FALSE ;
  NODE *var = node->down ;
  NODE *tagNode = node->right ;
  NODE *tagVar = tagNode ? tagNode->down : 0 ;
  NODE *tagTag = tagNode ? tagNode->right : 0 ;
  KEY key = 0 ;
  KEY tag = tagTag ? tagTag->key : 0 ;
  NODE *where = node->up ;
  int ii3, n3 = 0 ;
  KEYSET ks1 = 0, ks2 = 0, ks3 = 0 ;

  if (where && where->type != WHERE && where->type != WHERE_IF)
    where = 0 ;

  if (! tag && tagTag)
    tag = str2tag (stackText (bql->s,tagTag->mark)) ; 
  if (! tag)
    return FALSE ;
  if (! tagVar || ! tagVar->dclNode)
    return FALSE ;

  if (tagVar->dclNode && tagVar->dclNode->key)
    key = tagVar->dclNode->key ;
  
  ks1 = keySetCreate () ; /* list to be explored */
  ks2 = keySetCreate () ; /* additional keys */
  ks3 = keySetCreate () ; /* final list */
  
  if (key) 
    keySet (ks1, 0) = key ;

  while (keySetMax (ks1))
    {
      KEY key1 ;
      int n2 = 0, i = keySetMax (ks1) ;
      keySetMax (ks2) = 0 ;

      while (i--)
	{
	  key1 = keySet (ks1, i) ;
	  if(key1 != key) 
	    keySet (ks3, n3++) = key1 ;
	  if (keyFindTag (key1, tag))
	    {	
	      if (! bql->openObj)
		{
		  bql->mayBe = TRUE ;
		}
	      else
		{
		  
		  OBJ obj = bsCreate (key1) ;
		  KEY k, t = tag ;
		  int j ;
		  if (obj)
		    {
		      while (bsGetKey (obj, t, &k))
			{
			  t = _bsDown  ;
			  if (k != key && ! keySetFind (ks3,  k, &j))
			    keySet (ks2, n2++) = k ;		      
			} 
		      bsDestroy (obj) ;
		    }
		}
	    }
	}
      keySetSort (ks2) ;
      keySetCompress (ks2) ;
      keySetSort (ks3) ;
      keySetCompress (ks3) ;
      n3 = keySetMax (ks3) ;
      i = keySetMax (ks2) ;
      if (keySetMax (ks1) > i)
	keySetMax (ks1) = i ;
      while (i--)
	keySet (ks1, i) = keySet (ks2, i) ;
    }

  n3 = keySetMax (ks3) ;
  for (ii3 = 0 ; ii3 < n3 ; ii3++)
    {
      key = keySet (ks3, ii3) ;
      if (var->myObj)
	{
	  bsDestroy (var->obj) ;
	  var->myObj = FALSE ;
	}
      var->obj = 0 ;
      if (var->myAa)
	{
	  arrayDestroy (var->aa) ;
	  var->myAa = FALSE ;
	}
      var->aa = 0 ;
      if (var->myBsmark)
	{  
	  ac_free (var->bsMark) ;
	  var->myBsmark = FALSE ;
	}
      var->bsMark = 0 ;

      var->isNumber = var->isDate = FALSE ;
      var->uu.k = var->key = key ;
      var->col = var->row = 0 ;
      var->nCol = 1 ;
      var->parent = 0 ; 

      ok = TRUE ;
      if (where)
	{
	  bql->where = TRUE ;
	  ok = bqlExpandWhere (bql, where, coma) ;
	  bql->where = FALSE ;
	}
      if (ok)
	ok = bqlExpandNext (bql, coma) ;  
      ok2 |= ok ;
      if (coma && coma->right)
	bqlClean (bql, coma->right) ;
    }

  if (! ok2)
    {
      ok = TRUE ;
      if (where)
	{
	  bql->where = TRUE ;
	  ok = bqlExpandWhere (bql, where, coma) ;
	  bql->where = FALSE ;
	}
      if (ok)
	ok = bqlExpandNext (bql, coma) ; 
      ok2 |= ok ;
    }
  bqlClean (bql, var) ;

  keySetDestroy (ks1) ;
  keySetDestroy (ks2) ;
  keySetDestroy (ks3) ;

  return ok2 ;
} /* bqlExpandTTag */

/*************************************************************************************/
/* at this stage the only possibility is x[1]
 * more complex z in x[2]
 * have already been transformed as
 *  _z in x[1] , z in _z[1] 
 * which allows looping correctly on tag 2 construct
 */

static BOOL bqlExpandRightOf (BQL *bql, NODE *node, NODE *coma)
{
  BOOL ok = FALSE, ok2 = FALSE ;
  NODE *var = node->down ;
  NODE *tag = node->right ;
  NODE *tagVar = tag ? tag->down : 0 ;
  NODE *dcl = tagVar ? tagVar->dclNode : 0 ;
  NODE *where = node->up ;
  if (where && where->type != WHERE && where->type != WHERE_IF)
    where = 0 ;

  if (dcl && dcl->aa)
    {
      Array aa = dcl->aa ;
      int iMax = aa ? arrayMax (aa) : 0 ;
      int i, nCol = dcl->nCol ;
      int row, col ;
      BSunit *old = 0, *old0 ;
      
      if (var->depth < 0)
	var->depth = 0 ;
      col = dcl->col + var->depth ;
      var->col = col ; var->aa = dcl->aa ; var->nCol = dcl->nCol ;
      var->bsMark = dcl->bsMark ;
      if (dcl->parent && dcl->parent->obj)
	{
	  var->dclNode->aa = aa ; var->dclNode->nCol = dcl->nCol ;
	  var->dclNode->parent = dcl->parent ;
	  bsGoto (dcl->parent->obj, dcl->bsMark) ;
	  if (var->depth == 0)
	    {
	      var->uType = dcl->uType ;
	      var->key = dcl->key ;
	      var->isNumber = dcl->isNumber ;
	      var->isDate = dcl->isDate ;
	      var->z = dcl->z ;
	      var->uu = dcl->uu ;
	    }
	  else
	    {
	      for (i = 0 ; i < col ; i++)
		bsFindTag (dcl->parent->obj, _bsRight) ;
	      var->uType = bsType (dcl->parent->obj, _bsRight) ;
	      old0 = arrp (aa,dcl->row * nCol + dcl->col, BSunit) ;
	      for (row = dcl->row ; row * nCol + col < iMax ; row ++)
		{
		  BSunit *up0 = arrp (aa, row * nCol + dcl->col, BSunit) ;
		  BSunit *up = arrp (aa, row * nCol + col, BSunit) ;
		  
		  if (1)
		    {
		      BOOL new = FALSE ;
		      for (i = 0 ; !new && i <= + dcl->col ; i++)
			if ((up0 - i)->k != (old0 - i)->k)
			  new = TRUE ; /* we must not change the x[0] tag */
		      if (new)
			break ;
		    }
		  var->uu = up[0] ;
		  var->dclNode->row = row ;
		  if (var->uType == _Int)
		    {
		      var->isNumber = TRUE ;
		      var->isDate = FALSE ;
		      var->z = var->uu.i ;
		      var->key = 0 ;
		    }
		  else if (var->uType == _Float)
		    {
		      var->isNumber = TRUE ;
		      var->isDate = FALSE ;

		      var->z = var->uu.f ; 
		      var->key = 0 ;
		    }
		  else if (var->uType == _LongFloat)
		    {
		      var->isNumber = TRUE ;
		      var->isDate = FALSE ;
		      var->z = 0 ;
		      if (var->uu.s) sscanf (var->uu.s, "%lf", &(var->z)) ; 
		      var->key = 0 ;
		    }
		  else if (var->uType == _LongInt)
		    {
		      var->isNumber = TRUE ;
		      var->isDate = FALSE ;
		      var->z = 0 ;
		      if (var->uu.s) 
			{
			  long long int lli = 0 ;
			  sscanf (var->uu.s, "%lli", &(lli)) ;
			  var->z = lli ;
			}
		      var->key = 0 ;
		    }
		  else if (var->uType == _DateType)
		    {
		      var->isNumber = FALSE ;
		      var->isDate = TRUE ;

		      var->z = var->uu.time ; 
		      var->key = 0 ;
		    }
		  else if (var->uType > _LastC)
		    {
		      if (var->key != var->uu.k)
			{
			  if (var->myObj)
			    {
			      bsDestroy (var->obj) ;
			      var->myObj = FALSE ;
			    }
			  var->obj = 0 ;
			}
		      var->key = var->uu.k ;
		    }
		  else
		    var->isNumber = var->isDate = FALSE ;
		  
		  if (!old || old->k != up->k)
		    { 
		      ok = TRUE ;
		      if (where)
			{
			  bql->where = TRUE ;
			  ok = bqlExpandWhere (bql, where, coma) ;
			  bql->where = FALSE ;
			}
		      if (ok)
			{ 
			  if (where)
			    {
			      bql->where = TRUE ;
			      ok = bqlExpandWhere (bql, where, coma) ;
			      bql->where = FALSE ;
			    }
			  if (ok)
			    ok = bqlExpandNext (bql, coma) ; 
			  ok2 |= ok ;
			}
		    }
		  old = up ;
		}
	    }
	}
      else
	{
	  ok = TRUE ;
	  if (where)
	    {
	      bql->where = TRUE ;
	      ok = bqlExpandWhere (bql, where, coma) ;
	      bql->where = FALSE ;
	    }
	  if (ok)
	    ok = bqlExpandNext (bql, coma) ; 
	}
    }
  if (!ok2)
    {
      ok = TRUE ;
      var->uType = 0 ;
      if (where)
	{
	  bql->where = TRUE ;
	  ok = bqlExpandWhere (bql, where, coma) ;
	  bql->where = FALSE ;
	}
      if (ok)
	ok = bqlExpandNext (bql, coma) ; 
      ok2 |= ok ;
    }

  bqlClean (bql, var) ;

  return ok2 ;
} /* bqlExpandRightOf */

/*************************************************************************************/

static BOOL bqlExpandActive (BQL *bql, NODE *node, NODE *coma)
{
  BOOL ok = FALSE, ok2 = FALSE ;
  NODE *var = node->down ;
  KEYSET ks = 0 ;
  int i ;
  NODE *where = node->up ;
  if (where && where->type != WHERE && where->type != WHERE_IF)
    where = 0 ;
  
  ks = bql->ksIn ;
  for (i = 0 ; ks && i < keySetMax (ks) ; i++)
    {
      if (bql->maxLine && bql->results->rows > bql->maxLine)
	break ;
      ok  = TRUE ;
      if (var->myObj)
	{
	  bsDestroy (var->obj) ;
	  var->myObj = FALSE ;
	}
      var->obj = 0 ;
      var->parent = 0 ;

      var->key = keySet (ks, i) ;
      if (where)
	{
	  bql->where = TRUE ;
	  ok = bqlExpandWhere (bql, where, coma) ;
	  bql->where = FALSE ;
	}
      if (ok)
	{
	  bql->openObj = TRUE ;
	  bql->mayBe = FALSE ;
	  if (0)
	    {
	      bql->openObj = FALSE ;
	      ok = bqlExpandNext (bql, coma) ; 
	      bql->openObj = TRUE ;
	      if (ok && bql->mayBe) /* may be */
		ok = bqlExpandNext (bql, coma) ; 
	    }
	  else
	    ok = bqlExpandNext (bql, coma) ; 
	}
      ok2 |= ok ; 
      bqlClean (bql, node) ;
    }
  bqlClean (bql, var) ;

  return ok2 ;
} /* bqlExpandActive */

/*************************************************************************************/
/* return the number of line of the curly table implicitelly inside bql->counting */
static BOOL bqlExpandCurly (BQL *bql, NODE *node)
{
  NODE *down = node->down ;
  NODE *right = node->right ;

  bql->inCurly++ ;
  if (right)
    {
      bqlExpandFrom (bql, node, 0) ; 
    }
  else if (down)
    {
      node->right = down ;
      node->down = 0 ;
      bqlExpandFrom (bql, node, 0) ; 
      node->down = down ;
      node->right = right ; 
    }
  bql->inCurly-- ;

  return TRUE ;
} /* bqlExpandCurly  */

/*************************************************************************************/

static BOOL bqlExpandSet (BQL *bql, NODE *node, NODE *coma)
{
  NODE *var = node->down ;
  NODE *set = node->right ;
  NODE *dcl = set ? set->dclNode : 0 ;
  BOOL ok = FALSE ;
  double z = 0 ;
  NODE *where = node->up ;

  if (where && where->type != WHERE && where->type != WHERE_IF)
    where = 0 ;
  
  var->isNumber = var->isDate = FALSE ;
  switch (set->type)
    {
    case NUMBER:
    case DATE:
    case PLUS:
    case MINUS:
    case MULT:
    case DIVIDE:
    case POWER:
    case MODULO:
      {
	const char *ccp = bqlExpandEquation (bql, set, 0, &z) ;
	if (ccp == assVoid(1))
	  {
	    var->isNumber = TRUE ;
	    var->isDate = FALSE ;
	    var->z = z ;
	    var->key = 0 ;
	  }
	else if (ccp == assVoid(2))
	  {
	    var->isNumber = FALSE ;
	    var->isDate = TRUE ;
	    var->z = z ;
	    var->key = 0 ;
	  }
      }
      break ;
    case QUOTE:
    case DOUBLEQUOTE:
    case DOLLAR:
    case RAW:
      var->key = 0 ;
      var->uType = _Text ;
      var->uu.s = (char *)bqlExpandEquation (bql, set, 0, 0) ;
      break ;
    case VAR:
      if (dcl)
	{
	  var->key = dcl->key ;
	  var->z = dcl->z ;
	  var->uType = dcl->uType ;
	  var->uu = dcl->uu ;
	  var->isNumber = dcl->isNumber ;
	  var->isDate = dcl->isDate ;
	}
      break ;
    default:
      var->key = var->uType = 0 ;
      break ;
    }
  ok = TRUE ;
  if (where)
    {
      bql->where = TRUE ;
      ok = bqlExpandWhere (bql, where, coma) ;
      bql->where = FALSE ;
    }
  if (ok)
    bqlExpandNext (bql, coma) ; 
  bqlClean (bql, var) ;
  return ok ;
} /* bqlExpandSet */

/*************************************************************************************/
/* expand a subquery  {select .. from ... where }
 * for example  COUNT{} or {} union {}
 */
static BOOL bqlExpandFrom (BQL *bql, NODE *node, NODE *coma)
{
  NODE *right = node->right ;
  BOOL ok = FALSE ;

  /* store bql registers */
  /* bql->counting is handled by the calling routine bqlExpandCount */
  AC_HANDLE h = bql->h ;
  AC_TABLE table = bql->results ;
  int row = bql->tableRow ;
  Array titles = bql->titles ;
  KEYSET ksOut = bql->ksOut ;
  BOOL where = bql->where ;

  /* reinitialise bql registers */
  bql->h = ac_new_handle () ;
  bql->results = ac_db_empty_table (bql->db, 128, 1, bql->h) ;
  bql->tableRow = 0 ;
  
  bql->ksOut = 0 ;
  bql->where = FALSE ;

  /* restructure the subquery
   * should have been done as part of bqlParse
   * so we assume that we have the classi select.. from.. where 
   */
  /* expand and report */
  if (right && right->right)
    {
      NODE *rr = right->right ;
      if (rr->type == IN)
	ok = bqlExpandIn (bql, rr, 0) ;
      else if (rr->type == COMA && rr->down && rr->down->type == IN)
	ok = bqlExpandIn (bql, rr->down, rr) ;
      else if (rr->type == COMA && rr->down && rr->down->type == WHERE && rr->down->down && rr->down->down->type == IN)
	ok = bqlExpandIn (bql, rr->down->down, rr) ;
      else if (rr->type == COMA && rr->down && rr->down->type == WHERE_IF && rr->down->down && rr->down->down->type == IN)
	ok = bqlExpandIn (bql, rr->down->down, rr) ;
    }

  /* restore bql registers */
  bql->where = where ;
  bql->results = table ;
  bql->tableRow = row ;
  bql->titles = titles ;
  bql->ksOut = ksOut ;

  ac_free (bql->h) ;
  bql->h = h ;

  return ok ;
} /* bqlExpandFrom */

/*************************************************************************************/

static BOOL bqlExpandIn (BQL *bql, NODE *node, NODE *coma)
{
  NODE *right = node->right ;
  BOOL ok = FALSE ;

  if (right)
    switch (right->type)
      {
      case ACTIVE:
	ok = bqlExpandActive (bql, node, coma) ;
	break ;
      case CLASSE:
	ok = bqlExpandClass (bql, node, coma) ;
	break ;
      case OBJECT:
	ok = bqlExpandObject (bql, node, coma) ;
	break ;
      case COUNT:
	ok = bqlExpandCount (bql, node, coma) ;
	break ;
      case CLNAM:
	ok = bqlExpandClassName (bql, node, coma) ;
	break ;
       case NAM:
	ok = bqlExpandName (bql, node, coma) ;
	break ;
      case TIMESTAMP:
	ok = bqlExpandTimeStamp (bql, node, coma) ;
	break ;
      case DNA:
	ok = bqlExpandDNA (bql, node, coma) ;
	break ;
      case PEPTIDE:
	ok = bqlExpandPeptide (bql, node, coma, FALSE) ;
	break ;
      case CODING:
	ok = bqlExpandPeptide (bql, node, coma, TRUE) ;
	break ;
      case DATEDIFF:
	ok = bqlExpandDateDiff (bql, node, coma) ;
	break ;
      case MIN:
      case MAX:
      case SUM:
      case AVERAGE:
      case STDDEV:
	ok = bqlExpandMinMax (bql, node, coma) ;
	break ;
      case HASTAG:
      case MERGETAG:
      case TAG:
	ok = bqlExpandTag (bql, node, coma) ;
	break ;
      case TTAG:
	ok = bqlExpandTTag (bql, node, coma) ;
	break ;
      case RIGHTOF:
	ok = bqlExpandRightOf (bql, node, coma) ;
	break ;
      case CURLY:
	ok = bqlExpandCurly (bql, right) ;
	break ;
      default:
	ok = bqlExpandSet (bql, node, coma) ;
	break ;
      }
  return ok ;
} /* bqlExpandIn */

/*************************************************************************************/

static BOOL bqlExpandCount (BQL *bql, NODE *node, NODE *coma)
{
  NODE *var = node->down ;
  NODE *countTag = node->right ;
  BOOL ok = FALSE ;
  NODE *where = node->up ;
  if (where && where->type != WHERE && where->type != WHERE_IF)
    where = 0 ;
 
  if (countTag && countTag->type == COUNT && var && var->type == VAR)
    {
      var->isNumber = TRUE ;
      var->isDate = FALSE ;
      bql->counting = 1 ;
      if ( countTag->right)
	{
	  countTag->down = var ;
	  bqlExpandIn (bql, countTag, 0) ;
	}
      var->isNumber = TRUE ;
      var->isDate = FALSE ;
      var->z = bql->counting - 1 ;
      var->key = 0 ;
      if (var->myObj)
	{
	  bsDestroy (var->obj) ;
	  var->myObj = FALSE ;
	}
      var->obj = 0 ;

      bql->counting = 0 ;
      ok = TRUE ;
    }

  if (where)
    {
      bql->where = TRUE ;
      ok = bqlExpandWhere (bql, where, coma) ;
      bql->where = FALSE ;
    }
  else
    ok = TRUE ;
  if (ok)
    ok = bqlExpandNext (bql, coma) ; 
  bqlClean (bql, var) ;
  return ok ;
} /* bqlExpandCount */

/*************************************************************************************/
/* s v, t from v in @, t in MAX  v->terrain */
static BOOL bqlExpandMinMax (BQL *bql, NODE *node, NODE *coma)
{
  NODE *var = node->down ;
  NODE *countTag = node->right ;
  BOOL ok = FALSE ;
  NODE *where = node->up ;
  if (where && where->type != WHERE && where->type != WHERE_IF)
    where = 0 ;
 
  if (countTag && ! bql->minmaxavstd  &&
      (countTag->type == MIN || countTag->type == MAX || countTag->type == SUM || countTag->type == AVERAGE || countTag->type == STDDEV)
      && var && var->type == VAR)
    {
      var->isNumber = TRUE ;
      var->isDate = FALSE ;
      bql->counting = 1 ; 
      bql->minmaxavstd = countTag->type ;
      bql->minmaxavstdN = 0 ;
      bql->minmaxavstdX = 0;
      bql->minmaxavstdX2 = 0;
      if ( countTag->right)
	{
	  countTag->down = var ;
	  bqlExpandIn (bql, countTag, 0) ;
	}
      var->isNumber = TRUE ;
      var->isDate = FALSE ;
      switch (countTag->type)
	{
	case MIN:
	case MAX:
	  var->z = bql->minmaxavstdX ;
	  break ;
	case SUM:
	  var->z = bql->minmaxavstdX ;
	  break ;
	case AVERAGE:
	  var->z = bql->minmaxavstdN ? bql->minmaxavstdX/bql->minmaxavstdN : 0 ;
	  break ;
	case STDDEV:
	  var->z = bql->minmaxavstdN ? bql->minmaxavstdX2/bql->minmaxavstdN - (bql->minmaxavstdX/bql->minmaxavstdN)*(bql->minmaxavstdX/bql->minmaxavstdN) : 0 ;
	  break  ;
	default:
	  var->z = 0 ;
	  break ;
	}
      bql->counting = 0 ; 
      bql->minmaxavstd = 0 ;
      bql->minmaxavstdN = 0 ;
      bql->minmaxavstdX = 0;
            bql->minmaxavstdX2 = 0;
      ok = TRUE ;
    }

  if (where)
    {
      bql->where = TRUE ;
      ok = bqlExpandWhere (bql, where, coma) ;
      bql->where = FALSE ;
    }
  else
    ok = TRUE ;
  if (ok)
    ok = bqlExpandNext (bql, coma) ; 
  bqlClean (bql, var) ;
  return ok ;
} /* bqlExpandMinMax */

/*************************************************************************************/
/*
  bql select c,g,x from c in class compare, g in c->runs, a in g->ali, x in a->accepted where x == 256376
*/
static BOOL bqlExpandCondition (BQL *bql, NODE *node)
{
  BOOL ok = FALSE ;
  NODE *left = node->down ;
  NODE *right = node->right ;
  double z1 = 0 ;
  double z2 = 0 ;
  
  switch (node->type)
    {
    case ROUND:
      return left ? bqlExpandCondition (bql, left) : FALSE ;
    case WHERE_AVOID:
      return right ? bqlExpandCondition (bql, right) : TRUE ;
    case NUMBER:
      return node->z == 0 ? FALSE : TRUE ;
    case VAR: 
      node = node->dclNode ;
      if ((node->uType && node->uu.i)  || node->key || node->txt || node->isNumber || node->isDate || node->pepStack || node->dnaStack)
	return TRUE ;
      else
	return FALSE ;
    case DATE:
      if ((node->uType && node->uu.i)  || node->key || node->txt || node->isNumber || node->isDate)
	return TRUE ;
      else
	return FALSE ;
    case ISA:
      {
	KEY key = left->dclNode ?  left->dclNode->key : 0 ;
	KEY classFilter = node->uType ;
	KEY subClass = node->var ;

	if (!left || !lexIsKeyVisible(left->dclNode->key) ||
	    ! (class(key) == node->key) ||
	    (classFilter && ! lexIsInClass(key, subClass))
	    )
	  return FALSE ;
	else
	  return TRUE ;
      }
      break ;
    case AND:
      if (left && right && 
	  bqlExpandCondition (bql, left) && 
	  bqlExpandCondition (bql, right)
	  ) return TRUE ;
      return FALSE ;
      break ;
    case NOT:
      return right ? ! bqlExpandCondition (bql, right) : FALSE ;
      break ;
    case OR:
      if (left && bqlExpandCondition (bql, left)) return TRUE ;
      if (right && bqlExpandCondition (bql, right)) return TRUE ;
      return FALSE ;
      break ;
    case XOR:      
      if (left && bqlExpandCondition (bql, left)) ok = TRUE ;
      if (ok && ! right) return TRUE ;
      return ok ^ bqlExpandCondition (bql, right) ;
      break ;
    default:
      break ;
    }

  if (node->type == MERGETAG)
    return FALSE ;
  else if (node->type == HASTAG)
    {
      KEY key = left && left->dclNode ? left->dclNode->key : 0 ;
      KEY tag = right ? right->key : 0 ;
      if (! tag && right->type == VAR)
	right->key = tag = str2tag (stackText (bql->s, right->mark)) ;
      return key && tag && keyFindTag (key, tag) ;
    }
  else if (node->type == TAG)
    {
      KEY key = left && left->dclNode ? left->dclNode->key : 0 ;
      KEY tag = right ? right->key : 0 ;
      if (0) { return key && tag && keyFindTag (key, tag) ;}
      if (key && tag && keyFindTag (key, tag))
 	{
	  return keyGetKey (key,tag) ;
	}
      return FALSE ;

    }
  else 
    {
      int n ;
      const char *ccp1 = 0, *ccp2 = 0 ;
      char buf1[32], buf2[32] ;
      z1 = 0 ;
      int pass1 = 0, pass2 = 0 ;
      while ((ccp1 = bqlExpandEquation (bql, left, pass1, &z1)))
	{
	  pass1 = 2 ; pass2 = 0 ;
	  z2 = 0 ;
	  while ((ccp2 = bqlExpandEquation (bql, right, pass2, &z2))) 
	    {      
	      pass2 = 2 ;
	      if (ccp1 == assVoid(1))
		{
		  if (ccp2 == assVoid(1)) /* compare 2 numbers */
		    {
		      ok = bqlNumberCompare (node->type, z1, z2) ;
		      if (ok)
			return TRUE ;
		      else
			continue ;
		    }
		  else
		    {
		      sprintf (buf1, "%lli", (long long int)z1) ;
		      ccp1 = buf1 ;
		      if (left->mark)
			ccp1 = stackText (bql->s, left->mark) ;
		    }
		}
	      else if (ccp2 == assVoid(1))
		{
		  sprintf (buf2, "%lli", (long long int)z2) ;
		  ccp2 = buf2 ;
		  if (right->mark)
		    ccp2 = stackText (bql->s, right->mark) ;
		}
	      else if (ccp1 == assVoid(2) && ccp2 == assVoid(2))
		{
		  ok = bqlDateCompare (node->type, z1, z2) ;
		  if (ok)
		    return TRUE ;
		  else
		    continue ;
		}
	      if (
		  ! ccp1 || ! ccp2 ||
		  ccp1 == assVoid(1) || ccp2 == assVoid(1) ||
		  ccp1 == assVoid(2) || ccp2 == assVoid(2)
		  )
		  continue ;
	      if (*ccp1 && *ccp2)
		switch (node->type)
		  {
		  case LIKE:
		    n = pickMatch (ccp1, ccp2) ;
		    ok = (n > 0 ? TRUE : FALSE) ;
		    break ;
		  case LIKEREGEXP:
		    n = node->br ? regExpFind (node->br, ccp1) : 0 ;
		    ok = (n > 0 ? TRUE : FALSE) ;
		    break ;
		  case NEQ:
		  case LEQ:
		  case SEQ:
		  case LT:
		  case ST:
		  case EQ:
		    n = lexstrcmp (ccp1, ccp2) ;
		    switch (node->type)
		      {
		      case EQ:
			ok = (n == 0 ? TRUE : FALSE) ;
			break ;
		      case NEQ:
			ok = (n != 0 ? TRUE : FALSE) ;
			break ;
		      case LEQ:
			ok = (n >= 0 ? TRUE : FALSE) ;
			break ;
		      case SEQ:
			ok = (n <= 0 ? TRUE : FALSE) ;
			break ;
		      case LT:
			ok = (n > 0 ? TRUE : FALSE) ;
			break ;
		      case ST:
			ok = (n < 0 ? TRUE : FALSE) ;
			break ;
		      default:
			break ;
		      }
		  default:
		    break ;
		  }
	      if (ok) return TRUE ; /* true is preferred an breaks the loops */
	    }
	}
    }
  return ok ;
} /* bqlExpandCondition */

/*************************************************************************************/

static BOOL bqlExpandWhere (BQL *bql, NODE *node, NODE *coma)
{
  BOOL ok = FALSE ;
  NODE *down = node->down ;

  if (bql->where)
    {
      ok = TRUE ;
      if (node->right)
	ok = bqlExpandCondition (bql, node->right) ; 
    }
  else
    {
      if (down)
	switch (down->type)
	  {
	  case IN:
	    ok = bqlExpandIn (bql, down, coma) ;
	    break ;
	  case SET:
	    ok = bqlExpandSet (bql, down, coma) ;
	    break ;
	  default:
	    break ;
	  }
      if (node->type == WHERE_IF && ! ok)
	ok = bqlExpandNext (bql, coma) ;
    }

  return ok ;
} /* bqlExpandWhere */

/*************************************************************************************/

static BOOL bqlExpandComa (BQL *bql, NODE *node)
{
  NODE *coma = 0 ;
  BOOL ok = FALSE ;
  
  if (node->type == COMA)
    {
      coma = node ;
      node = coma->down ;
    }

  if (node)
   switch (node->type)
    {
    case WHERE:
      if (bql->where)
	 messcrash ("where inside where") ;
      ok = bqlExpandWhere (bql, node, coma) ;
      break ;
    case WHERE_IF:
      if (bql->where)
	 messcrash ("where_if inside where") ;
      ok = bqlExpandWhere (bql, node, coma) ;
      break ;
    case IN:
      ok = bqlExpandIn (bql, node, coma) ;
      break ;
    case SET:
      ok = bqlExpandSet (bql, node, coma) ;
      break ;
    default:
      break ;
    }
  return ok ;
} /* bqlExpandComa */

/*************************************************************************************/
/*
  bql select c,r from c in class compare,r in  c->runs
*/

static BOOL bqlExpand (BQL *bql, NODE *node)
{
  NODE *right = node ? node->right : 0 ;
  BOOL ok = FALSE ;

  if (right)
    ok = bqlExpandComa (bql, right) ;

  return ok ;
} /* bqlExpand */

/*************************************************************************************/
/*************************  Public interface *****************************************/
/*************************************************************************************/

const char *bqlError (BQL *bql)
{
  if (! bql || ! bql->errTxt)
    return 0 ;
  return vtxtPtr (bql->errTxt) ;
} /* bqlError */ 

/*************************************************************************************/

AC_TABLE bqlResults (BQL *bql)
{
  if (bql &&  bql->results)
    {
      if (! bql->isSorted && ! bql->doNotSort)
	bqlSort (bql) ;
      
      return bql->results ;
    }
  return 0 ;
} /* bqlResults */

/*************************************************************************************/
/* Create a a BQL structure */
BQL *bqlCreate (BOOL debug, AC_HANDLE h) 
{
  BQL *bql = (BQL *) halloc (sizeof (BQL), h) ;
  /* As these structure are availablemade  to the caller, 
   * they are allocated on the parent handle 
   */
  bql->results = ac_db_empty_table (bql->db, 128, 1, h) ;
  bql->tableRow = 0 ;
  bql->errTxt = vtxtHandleCreate (h) ;
  
  bql->debug = debug ;
  bql->h = ac_new_handle () ;	
  bql->dict = dictCaseSensitiveHandleCreate (50, bql->h) ;
  bql->dcls = arrayHandleCreate (128, DCL, bql->h) ;
  bql->froms = arrayHandleCreate (128, NODE*, bql->h) ;
  bql->wheres = arrayHandleCreate (128, NODE*, bql->h) ;
  bql->s = stackHandleCreate (1024, bql->h) ;
  pushText (bql->s, "Zero") ;
  
  
  bql->node = (NODE *) halloc (sizeof (NODE), bql->h) ;
  bql->node->mark = stackMark (bql->s) ;
  bql->node->type = RAW ;
  
  /* register destroy */
  bql->magic = uBqlDestroy ;
  blockSetFinalise (bql, uBqlDestroy) ;
  return bql ;
} /* bqlCreate */

/********************************************/

BQL *bqlDbCreate (AC_DB db, BOOL debug, AC_HANDLE h) 
{
  BQL *bql = bqlCreate (debug, h) ;
  bql->db = db ;
  return bql ;
} /* bqlDbCreate  */

/********************************************************************/
/* To manipulate large tables you may optionally
 * process table one KEY at a time
 * 1: iter = bqlIterCreate ("bql query", ksIn, &errors, h) ;
 * 2: while ((mini_table = bqlIterTable (iter)) 
 *        { do_something_with_the_mini_table ; }
 *     The mini table corresponds to the next key producing some data
 *     returns NULL when no more key produces a non empty mini_table
 *     DO NOT free the mini_table, it is freed internally on each iteration
 * 3: Call ac_free (iter) OR ac_free (h) ;
 */
BQL_ITER *bqlIterCreate (const char *bqlQuery, KEYSET ksIn, const char **errors, AC_HANDLE h) 
{
  BOOL debug = FALSE ;
  BQL_ITER *iter = (BQL_ITER *) halloc (sizeof (BQL_ITER), h) ;

  iter->magic = uBqlIterDestroy ;
  blockSetFinalise (iter, uBqlIterDestroy) ;

  iter->h = ac_new_handle () ;
  iter->bql = bqlCreate (debug, iter->h) ;

  if (!bqlParse (iter->bql, bqlQuery, FALSE))
    {
      const char *ccp = bqlError (iter->bql) ;
      if (ccp) *errors = strnew (ccp, h) ;
      else *errors = strnew ("Unknown error in bqlParse", h) ;
      ac_free (iter) ;
    }
  if (! ksIn)
    messcrash ("bqlIterCreate received a NULL ksIn") ;
  
  if (iter)
    {
      iter->currentRow = 0 ;
      iter->ksIn = arrayHandleCopy (ksIn, iter->h) ;/* make a private copy */
      iter->ksCurrent = keySetHandleCreate (iter->h) ;
    }
  return iter ;
} /* bqlIterCreate  */

/*************************************************************************************/

AC_TABLE bqlIterTable (BQL_ITER *iter)
{
  int nn = 0 ;

  if (! iter)
    messcrash ("bqlIterTable received a NULL iter") ;
  if (iter->magic != uBqlIterDestroy)
    messcrash ("bqlIterTable received an invalid iter") ;

  ac_free (iter->bql->results) ;
  while (nn == 0 && iter->currentRow < keySetMax (iter->ksIn))
    {
      ac_free (iter->bql->results) ;
      keySet (iter->ksCurrent, 0) = keySet (iter->ksIn, iter->currentRow++) ;
      bqlRun (iter->bql, iter->ksCurrent, 0) ;
      nn = iter->bql->results ? iter->bql->results->rows : 0 ;
    }

  if (nn == 0)
    {
      ac_free (iter->h) ;
      return 0 ;
    }
  return iter->bql->results ;
} /* bqlIterTable */

/*************************************************************************************/

void bqlMaxLine (BQL *bql, int maxLine)
{
  if (bql && maxLine >= 0)
    bql->maxLine = maxLine ;
} /* bqlMaxLine */

/*************************************************************************************/

BOOL bqlParse (BQL *bql, const char *query, BOOL acedbQuery)
{
  BOOL ok = TRUE ;
  const char *cp = 0 ;
   
  if (! query || ! *query)
    {
      vtxtPrintf (bql->errTxt, "// BQL received an empty query") ;
      return FALSE ;
    }
  
  bql->titles = arrayHandleCreate (12, char *, bql->h) ;

  if (acedbQuery)
    return bqlParseAcedbQuery (bql, query) ;
  

  vtxtClear (bql->errTxt) ;
  vtxtPrintf (bql->errTxt, "// ... BQL query : \n// ...  %s\n\n", query) ;
  
  bql->node->mark = stackMark (bql->s) ;
  bql->node->type = RAW ;
  cp = query ;
  while (*cp == ' ') cp++ ;
  if (strncasecmp (cp, "select ", 7))
    {
      pushText (bql->s, "select ") ;
      catText (bql->s, cp) ;
    }
  else
    {
      pushText (bql->s, cp) ;
    }

  /* split out the terminal FORCE_ZERO */
  if (1)
    {
      char *cp = stackText (bql->s, bql->node->mark) ;
      char *cq = strstr (cp, "FORCE_ZERO") ;
      
      if (cq)
	{
	  memset (cq, ' ', 10) ;
	  bql->force_zero = TRUE ;
	}
    }

  /* back compatibility with aql */
  {
    char *cp = stackText (bql->s, bql->node->mark) ;
    char *cq = cp ;

    cq = cp - 1;
    while ((cq = strstr (cq + 1, "elect ")))
      {
	if (! strncasecmp (cq - 1, "select ", 7))
	  {
	    cq += 6 ;
	    while (*cq == ' ') cq++ ;       /* replace 'select all ' by 'select' */
	    if (! strncmp (cq, "all ", 4))
	      memset (cq, ' ', 4) ;
	  }
      }

    cq = cp ;
    while ((cq = strstr (cq, "@active:1")))
      memset (cq+1, ' ', 8) ; 
    cp = cq = stackText (bql->s, bql->node->mark) ;
    while ((cq = strstr (cq, ".name")))
      memset (cq, ' ', 5) ; 
    cp = cq = stackText (bql->s, bql->node->mark) ;
    while ((cq = strstr (cq, " exists ")))
      memset (cq, ' ', 8) ; 
    cp = cq = stackText (bql->s, bql->node->mark) ;
    while ((cq = strstr (cq, " exists_tag ")))
      memset (cq, ' ', 12) ; 
    cp = cq = stackText (bql->s, bql->node->mark) ;
    while ((cq = strstr (cq, "->DNA")))
      cq[2] = 'd' ;  /* prevent recognizing tag DNA as function DNA */
    cp = cq = stackText (bql->s, bql->node->mark) ;
    while ((cq = strstr (cq, ">TITLE")))
      cq[1] = 't' ;  /* prevent recognizing tag TITLE as function TITLE */
    cp = cq = stackText (bql->s, bql->node->mark) ;
    while ((cq = strstr (cq, "->PEPTIDE")))
      cq[2] = 'p' ;  /* prevent recognizing tag PEPTIDE as function PEPTIDE */
    cp = cq = stackText (bql->s, bql->node->mark) ;
    while ((cq = strstr (cq, "->CODING")))
      cq[2] = 'c' ;  /* prevent recognizing tag CODING as function CODING */
    cp = cq = stackText (bql->s, bql->node->mark) ;
    while ((cq = strstr (cq, "->DATEDIFF")))
      cq[9] = 'f' ;  /* prevent recognizing tag DATEDIFF as function DATEDIFF */
    cp = cq = stackText (bql->s, bql->node->mark) ;
    while ((cq = strstr (cq, "#DNA")))
      cq[1] = 'd' ;  /* prevent recognizing tag DNA as function DNA */
    cp = cq = stackText (bql->s, bql->node->mark) ;
    while ((cq = strstr (cq, "#PEPTIDE")))
      cq[1] = 'p' ;  /* prevent recognizing tag PEPTIDE as function PEPTIDE */
    cp = cq = stackText (bql->s, bql->node->mark) ;
    while ((cq = strstr (cq, "#CODING")))
      cq[1] = 'c' ;  /* prevent recognizing tag CODING as function CODNG */
    cp = cq = stackText (bql->s, bql->node->mark) ;
    while ((cq = strstr (cq, "#DATEDIFF")))
      cq[8] = 'f' ;  /* prevent recognizing tag DATEDIFF as function DATEDIFF */
 }
  
  /* order by == order   by = order_by */
  if (1)
    {
      char *cp = stackText (bql->s, bql->node->mark) ;
      char *cq = strstr (cp, "no_order") ;
      if (cq)
	{
	  bql->doNotSort = TRUE ; 
	  memset (cq, ' ', 8) ;
	}
    }

  {
    char *cp = stackText (bql->s, bql->node->mark) ;
    while ((cp = strstr (cp, "order")))
      {
	char *cr ;
	cp += 5 ;
	if (! strncmp (cp, "_by", 3))
	  continue ;
	cr = cp ;
	while (cr[1] == ' ') cr++ ;
	if (! strncmp (cr, " by", 3))
	  {
	    cr[1] = cr[2] = ' ' ;
	    cp[0] = '_' ; cp[1] = 'b' ; cp[2] = 'y' ;
	  }
      }
  }
  
  /* split out the terminal title command */
  if (1)
    {
      char *cp = stackText (bql->s, bql->node->mark) ;
      char *cq = strstr (cp, "TITLE") ;
      char *cr ;
      
      while (cq && (cr = strstr (cq + 1, "TITLE")))
	cq = cr ;
      if (cq && strchr (cq, ';'))
	cq = 0 ;
      if (cq)
	{
	  cr = cq ;
	  cq += 5 ;
	  while (*cq == ' ') cq++ ;
	  if (cq && *cq) 
	      bql->title_by = strnew (cq, bql->h) ;
	  while (*cr) *cr++ = ' ' ;
	}
    }


  /* split out the terminal order_by command */
  if (1)
    {
      char *cp = stackText (bql->s, bql->node->mark) ;
      char *cq = strstr (cp, "order_by") ;
      char *cr ;
      
      while (cq && (cr = strstr (cq + 1, "order_by")))
	cq = cr ;
      if (cq && strchr (cq, ';'))
	cq = 0 ;
      if (cq)
	{
	  cr = cq ;
	  cq += 9 ;
	  while (*cq == ' ') cq++ ;
	  if (cq && *cq) 
	      bql->order_by = strnew (cq, bql->h) ;
	  while (*cr) *cr++ = ' ' ;
	}
    }
  
 if (1)
    {
      char *cp = stackText (bql->s, bql->node->mark) ;
      char *cq = cp + strlen (cp) - 1 ;
      while (cq > cp && *cq == ' ') *cq-- = 0 ;
    }

  ok = TRUE ;
  if (1) bqlSplitOnSemicolumn (bql, bql->node) ; 
  while (ok && bqlGetBrackets (bql, bql->node, &ok)) {} ;
  while (ok && bqlCheckForExtraBrackets (bql, bql->node, &ok)) {} ;
  while (ok && bqlGetQuotes (bql, bql->node, &ok)) {} ;


  if (1) bqlCreatePhonyVariables (bql, bql->node) ;  /* may create new "" */
  while (ok && bqlGetBrackets (bql, bql->node, &ok)) {} ;
  while (ok && bqlCheckForExtraBrackets (bql, bql->node, &ok)) {} ;
  while (ok && bqlGetQuotes (bql, bql->node, &ok)) {} ;

  if (ok) bqlSpaceProtectSymbols (bql, bql->node, &ok) ;
  while (ok && bqlTokenize (bql, bql->node, &ok)) {} ;
  while (ok && bqlGetTypes (bql, bql->node, &ok)) {} ;
  if (ok) ok = bqlVar (bql, bql->node) ;
  while (ok && bqlReorder (bql, bql->node, &ok)) {} ;
  while (ok && bqlCleanParenthesis (bql, bql->node)) {} ;
  
  if (ok)
    ok = bqlCheckSides (bql, bql->node) ;
  
  if (ok)
    ok = bqlCheckSyntax (bql, bql->node, 0) ;
  if (ok)
    ok =  bqlCheckVariableDeclarations (bql, 0, 0) ;
 
  bqlSetSemiCols (bql) ;
  bqlSetComas (bql) ;

  if (ok && bqlSetCurly (bql, bql->node)) {} ;
  if (ok && bqlSetAvoid (bql, bql->node)) {} ;

  if (ok)
    ok = bqlSetDclNodeDepth(bql, bql->node) ; /* set node->dclNode and add a coma in FROM clause */
 
  if (ok)
    ok = bqlCheckOrderBy (bql) ;
  if (ok)
    ok = bqlCheckTitleBy (bql) ;
  if (bql->debug)
    showBql (bql) ;

  if (! ok)
    {
      if (bql->debug)
	fprintf (stderr, "BQL ERROR :\n%s\n", bqlError (bql));
    }
 
  return ok ;
} /* bqlParse */

/*************************************************************************************/
/* Parse old style queries
 * Find Author S* ; >Paper journal = nature
 */
static char *bqlParseAcedbCondition (BQL *bql, vTXT txt, char *cp, int nQ)
{
  char *cq, *cr ;
  char nam[12] ;

  sprintf(nam, "a%05d", nQ) ;
  cr = strchr (cp, ';') ;
  if (cr)
    *cr = 0 ;
  
  while (*cp == ' ') cp++ ;
  cq = cp + strlen(cp) - 1 ;
  while (cq > cp && *cq == ' ')
    *cq-- = 0 ;

  if (*cp)
    {
      KEY tag = 0 ;
      if (*cp == '\"')
	vtxtPrintf (txt, " where %s ~ %s", nam, cp) ;
      else if (lexword2key (cp, &tag, 0))
	vtxtPrintf (txt, " where %s#%s ", nam, cp) ;
      else if (! strncasecmp (cp, "IS ", 3) ||
	       ! strncasecmp (cp, "IS<", 3) ||
	       ! strncasecmp (cp, "IS>", 3)
	       ) 
	{
	  cp += 2 ; while (*cp == ' ') cp++ ;
	  if (strchr (cp, '>') ||  strchr (cp, '<'))
	    {
	      vtxtPrintf (txt, " where %s %c", nam, *cp++) ;
	      if (*cp == '=')
		vtxtPrintf (txt, "%c", *cp++) ;
	      while (*cp == ' ') cp++ ;
	      vtxtPrintf (txt, " %s"
			  , *cp == '"' ? cp : ac_protect(cp, bql->h)
			  ) ;
	    }
	  else
	    vtxtPrintf (txt, " where %s ~ %s", nam
			, *cp == '"' ? cp : ac_protect(cp, bql->h)
			) ;
	}
      else if (strchr (cp, '>') ||  strchr (cp, '<') || strchr (cp, '='))
	vtxtPrintf (txt, " where %s->%s ", nam, cp) ;
     else
       vtxtPrintf (txt, " where %s ~ %s", nam
		   , *cp == '"' ? cp : ac_protect(cp, bql->h)
		   ) ;
    }
  else
    {
      vtxtPrintf (txt, " where %s ", nam) ;
    }
  if (cr) cp = cr + 1 ;
  else *cp = 0 ;
  while (*cp == ' ') cp++ ;
  return cp ;
} /* bqlParseAcedbCondition  */

/*************************************************************************************/

static BOOL bqlParseAcedbQuery (BQL *bql, const char *query)
{
  BOOL debug = bql->debug ;
  BOOL ok = FALSE ;
  vTXT txt = vtxtHandleCreate (bql->h) ;
  char nam[32] ;
  int nQ = 0 ;
  char cc, *cp, *cq, *cr ;
  cp = strnew (query, bql->h) ;
  while (*cp == ' ') cp++ ;

  if (! strncasecmp (cp, "find ", 5))
    {
      cq = cp + 5 ;
      while (*cq == ' ') cq++ ;
      cr = cq ;
      if (*cr == '\"')
	{
	  while (*cr && *cr != '\"')
	    cr++ ;
	}
      else
	{
	  while (*cr && *cr != ';' && *cr != ' ')
	    cr++ ;
	}
      cc = *cr ;
      *cr = 0 ;
      ok = TRUE ;
      sprintf(nam, "a%05d", ++nQ) ;
      vtxtPrintf (txt, "select %s from %s in class %s", nam, nam, cq) ;
      *cr = cc ;
      ok = TRUE ;
      cp = bqlParseAcedbCondition (bql, txt, cr, nQ) ;   
     }

  while (*cp == '>' || ! strncasecmp (cp, "follow ", 5))
    {
      if (*cp == '>') 
	cq = cp + 1 ;
      else  
	cq = cp + 7 ;

      while (*cq == ' ') cq++ ;
      cr = cq ;
      if (*cr == '\"')
	{
	  while (*cr && *cr != '\"')
	    cr++ ;
	}
      else
	{
	  while (*cr && *cr != ';' && *cr != ' ')
	    cr++ ;
	}
      cc = *cr ;
      *cr = 0 ;
      ok = TRUE ;
      if (nQ == 0)
	{
	  sprintf(nam, "a%05d", ++nQ) ;
	  vtxtPrintf (txt, "select %s from a in @, %s in a->%s", nam, nam, cq) ;
	}
      else
	{
	  sprintf(nam, "a%05d", ++nQ) ;
	  vtxtPrintf (txt, ", %s ", nam) ;
	  sprintf(nam, "a%05d", nQ - 1) ;
	  vtxtPrintf (txt, " in %s->%s", nam, cq) ;
	  sprintf(nam, "a%05d", nQ) ;
	}
      *cr = cc ;
      ok = TRUE ;
      cp = bqlParseAcedbCondition (bql, txt, cr, nQ) ;   
    }
  if (nQ > 1)
    {
      cp = vtxtPtr (txt) ;
      sprintf(nam, "a%05d", 1) ;
      cq = strstr (cp, nam) ;
      sprintf(nam, "a%05d", nQ) ;
      memcpy (cq, nam, 6) ;
    }

  if (ok)
    {
      if (debug) fprintf (stderr, "// query: %s\n", vtxtPtr (txt)) ; 
      ok = bqlParse (bql, vtxtPtr (txt), FALSE) ;
    }
  return ok ;
} /* bqlParseAcedbQuery */

/*************************************************************************************/

int bqlRun (BQL *bql, KEYSET activeKeyset, KEYSET resultKeyset)
{
  int nn = 0 ;
  BOOL ok = TRUE ;
  BOOL debug = bql->debug ;
  NODE *node ;
  NODE *node0 = bql->node ;
  NODE *from0 = bql->from ;

  bql->ksIn =  activeKeyset ;
  bql->ksOut = resultKeyset ?  resultKeyset : keySetHandleCreate (bql->h) ;
  keySetMax (bql->ksOut) = 0 ;

  bql->openObj = TRUE ;

  node = node0 ;
  while (ok && node)
    {
      NODE *down = node ;
      if (node->type == SEMICOLUMN)
	down = node->down ;
      if (! down ||
	  (down->type == SELECT && ! down->down)
	  )
	{ node = node->right ; continue ;  }
      if (down->type == ORDER_BY)
	down = down->down ;
      if (bql->ksIn == bql->ksOut)
	{
	  bql->ksIn = keySetHandleCopy (bql->ksOut, bql->h) ;
	  keySetMax (bql->ksOut) = 0 ;
	}
      bql->node = down ;
      bql->from = down ;
      bql->results = ac_db_empty_table (bql->db, 128, 1, bql->h) ;
      bql->tableRow = 0 ;
      if (bql->ksIn) arraySort(bql->ksIn,keySetAlphaOrder) ;
      ok = bqlExpand (bql, down) ;
      bql->isSorted = FALSE ;
      if (! ok)
	continue ;
      keySetSort (bql->ksOut) ;
      keySetCompress (bql->ksOut) ;
      bql->ksIn = bql->ksOut ;
      node = node->right ;
      nn = bql->results->rows ;
    }

  bql->node = node0 ;
  bql->from = from0 ;

  if (! ok)
    {
      if (debug)
	fprintf (stderr, "BQL ERROR : %s\n", bqlError (bql));
    }
  
  return nn ;
} /* bqlRun */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/* examples
bql select a from a = "aaa"
bql select a from a = "aaa" where a
bql select a from a = "aaa" where ! a
bql select a from a = "aaa  bbb"
bql select a,b from a = 3*8 , b = "bbb"
bql select a,b from a = 3*8 , b = 2 * a - 1 
bql select a from a =  (2+3) * (1+1)
bql select a,b from a = 3*8 where a > 20 , b = 2 * a - 1

bql select x from x = 2
bql select x from x = -2
bql select x from x = (2)
bql select x from x = (-2)
bql select x from x = (2) * (-2)
bql select x from x = 2+2
bql select x from x = (2+2)
bql select x from x = ((2+2))
bql select x from x = ((2+2)
bql select x from x = (2+2))
bql select x from x = (3*(2+2))
bql select x from x = ((2+2)*3)

bql select a,b from a = "aaa", b="bbb"
bql select c,g,r from c in class compare, g in c->runs, r in g->union_of where g == "SubFR"
bql select c,g,r from c in class compare, g in c->runs where g == "SubFR", r in g->union_of
bql select c,g,r from c in class compare, g in c->runs,r = (2+3) * 2
bql select c,g,r from c in class compare, g in c->runs,r = 2*(2+3)
bql select c,g,r from c in class compare, g in c->runs,r = (2+3) * (1+1)

bql select c,g from c in class compare where c == 2, g in c->runs
bql select c,g from c in class compare where c == "cov", g in c->runs

bql select c,g from c in class compare, g in c#runs
bql select c,g from c in class compare, g in c->runs
bql select c,c#runs from c in class compare
bql select c,c->runs from c in class compare
bql select c,n from c in class compare, n in count c->runs where n > 2
bql select c,n from c in class compare, n in count c->runs where n > 5

bql select a from a in class ali where a >= "e" && a < "q"
bql select a,x from a in class ali where a >= "e" && a < "q", x in a->accepted
bql select a,x,y from a in class ali where a >= "e" && a < "q", x in a->accepted, y in x:2

bql select a,b,x from a in @, b in a->accepted, x in b:4 where x > 	12754
bql select a,b,x from a in @, b in a->accepted  where b > 	12754 , x in b:4

bql select a,b,x from a in class ali, b in a->accepted , x in b[4]

bql select a,b,x,z from a in class ali, b in a->gene_expression , x in b[1] where x == "av", z in x[1]

bql select g,g->union_of from g in class run

bql select a,b from   a in class compare , b in a#apply


bql select a,b,c from   a in class compare , b in a#apply, c in b:1
bql select a,b,c from   a in class compare , b in a#apply, c in b[1] where b and c
find ali
bql select a,b,c from a in @ , b in a->accepted, c in b[4]
bql select a,b,c,d from a in @ , b in a#accepted, c in b[3], d in b[5]

bql select a from a in object ("Compare","cov") where a >= "c"
bql select a from a in object ("Run", "P*")

bql select r  from r in class run   where  COUNT r->ali == 1
bql select r , COUNT r->ali  from r in class run 
bql select a,x from a in class ali, x in a->accepted[2]
bql select a,x from a in class ali, x in a->accepted[5] where a ~ "sub*"


// accept the acedb language
bql -q find compare ; > runs ind* ; > union_of titles

bql select a,x,y,z from a in class ali, x in a->quality_profile , y in x[1], z in x[2]
bql select a,x,y,z from a in class ali, x in a->quality_profile , y in x[1], z in y[1]


find ali subfr
bql select a,h,x from a in @, h in a->h_ali where h == "any", x in h[3]



bql select a,e from a in class ali where a == "subfr", e in a=>error_profile
bql select a,e from a in class ali where a == "subfr", e in a=>error_profile[1]
bql select a,e from a in class ali where a == "subfr", e in a=>error_profile[2]
bql select a,e from a in class ali where a == "subfr", e in a=>error_profile[3]
bql select a,e from a in class ali where a == "subfr", e in a=>error_profile[4]
bql select a,e from a in class ali where a == "subfr", e in a=>error_profile[5]


find sequence toto
bql select s,d from s in @, u=3, d in DNA(s,u,2*u)
bql select s,pep from s in @, pep in PEPTIDE(s,1,1)

bql DNA(@)


bql -q find ali IS >= pa
bql -q find ali IS >= "pa"

find ali
bql select a,x from a in @, x in a->quality_profile order_by -2

bql select g,r,a from g in class run, r in g->union_of, a in r->ali
bql select g,r,a from g in class run, r in g->union_of, a in r->ali where r
bql select g,a from g in class run, r in g->union_of, a in r->ali where r

query find ali 
bql @
bql @,@->run
bql select a,x,y from a in @, x in a->quality_profile , y in x[1] order_by -1

# this one works
query find ali ! accepted
bql select a,b,c,d from a in @ , b in a#accepted, c in b[3], d in b[5]
query find ali accepted
bql select a,b,c,d from a in @ , b in a#accepted, c in b[3], d in b[5]
# this one crashes probably because of a mixup: c->parent->obj is a bad address
query find ali IS >= subfr && IS <= Test.sum
bql select a,b,c,d from a in @ , b in a#accepted, c in b[3], d in b[5]

find ali subfr
bql select a,b,x from a in @, b in a->Genes_touched, x in b[1]

bql select a,b,x from a in @, b in a->gene_expression , x in b:1

find ali subfr
bql select a,b from a in @, b in a->gene_expression where b ~ "Mb_in_high_genes"

bql select a,b,x,z from a in @, b in a->gene_expression , x in b:1 where x == "av", z in x[1]

query find run union_of
bql select s from s in @ where COUNT { select r from r in s->union_of where r#paired_end } > 1
bql select s from s in @ where COUNT { select r from r in s->union_of where COUNT r->ali == 1} == 4

bql select c,d,e from c in class compare ,d in c#apply,e in d[1]
find usersession "*.1_*"
bql -a -c 5 select s,d from s in @active:1, d in s->start

find compare
bql select c,a from c in @, a in c->covariance where a ~ "autosort"

query find variant mm
bql -c 1 select s,d,x,y from s in @, d in s->mm, x in d[1], y in x[1]
bql select s,d,x,y from s in @, d in s->mm, x in d[1], y in d[2]

select g in ?run  where g#union_of AND NOT g#compare

#FAUX probleme x est soit Float soit Int selon la ligne, mais bql pense Int
# il faudrait avoir une fonction ac_tag_table_type qui renvoie le array des types
bql select a,b,x,z from a in @, b in a->gene_expression , x in b:1 where x == "av", z in x[1]

# FAUX on m'a pas le non nombre de lignes 
query find ali tutu
bql select X1, X2, X3 from X1 in @ , X2 in X1->Error_profile where X2 == "xx" , X3 in X2[1]
bql select X1, X2, X3 from X1 in @ , X2 in X1->Error_profile, X3 in X2[1]

bql  select X1, X2, X3, X4 from X1 in @ , X2 in X1->Error_profile, X3 in X2[1], X4 in X3[1]
bql  select X1, X2, X3, X4, X5 from X1 in @ , X2 in X1->Error_profile, X3 in X2[1], X4 in X3[1], X5 in X4[1]


bql  select X1, X2, X3, X4, X5, X6 from X1 in @ , X2 in X1->Error_profile, X3 in X2[1], X4 in X3[1], X5 in X4[1], X6 in X5[1]
bql select X1, X2, X3, X4, X5 from X1 in @ , X2 in X1->Error_profile where X2 == "xx" , X3 in X2[1], X4 in X3[1], X5 in X4[1]
bql select X1, X2, X3, X4, X5, X6 from X1 in @ , X2 in X1->Error_profile where X2 == "xx" , X3 in X2[1], X4 in X3[1], X5 in X4[1], X6 in X5[1]

find ali "PairIndel100_180"
bql -c 2 select X1, X2, X3 from X1 in @ , X2 in X1->Error_profile, X3 in X2[1]
find ali "PairIndel100_180"
bql -c 2 select X1, X2, X3, X4 from X1 in @ , X2 in X1->Error_profile, X3 in X2[1], X4 in X3[1]
find ali "PairIndel100_180"
bql -c 2 select X1, X2, X3, X4, X5, X6 from X1 in @ , X2 in X1->Error_profile, X3 in X2[1], X4 in X3[1], X5 in X4[1], X6 in X5[1]


bql -o tutu2 select X1, X2, X3, X4, X5 from X1 in class "Ali", X2 in X1->Error_profile, X3 in X2[1], X4 in X3[1], X5 in X4[1]
bql  select X1, X2, X3, X4, X5 from X1 in class "Ali", X2 in X1->Error_profile, X3 in X2[1], X4 in X3[1], X5 in X4[1]

# like gives empty answer but ~ works, they are not equivalent
bql select r from r in class "run" where r like "Pa*0"
bql select r from r in class "run" where r ~ "Pa*0"

select ?sequence toto
 select @, DNA(@)

# error  y shows x, it should be empty
select ?sequence 1a50
bql select s,x,y from s in @, x in s->dna[2],y in x[1]   // 1A50	43054	NULL
bql select s,x,y,z from s in @, x in s->dna[2],y in s->dna[1],z in y[1] // 1A50	43054	1A50	43054
bql select s,x,y,z from s in @, x in s->dna[1],y in x[1],z in y[1] // 1A50	1A50	43054	NULL
bql select s,x,y,z from s in @, x in s->dna[0],y in x[0],z in y[1] // DNA DNA  1A50
bql select s,x,y,z from s in @, x in s->dna[0],y in x[0],z in y[2] //DNA DNA 43054
bql select s,x,y,z from s in @, x in s->dna[0],y in x[1],z in y[1] // DNA  1A50 43054
bql select s,x,y,z from s in @, x in s->dna[0],y in x[2],z in y[1] // DNA 43054 null
bql select s,x,y,z from s in @, x in s->dna[0],y in x[3],z in y[1] // DNA null null


show Map // Map      CHROMOSOME_I    1   43054
bql select s,m,x,y from s in @, m in s->map, x in m[1], y in x[2]     // CHROMOSOME_I end 1/43054
bql select s,m,x,y from s in @, m in s->map, x in m[1], y in x[1]     // CHROMOSOME_I ends left/right
bql select s,m,x,y from s in @, m in s->map[0], x in m[1], y in x[1]    // Map chr ends
bql select s,m,x,y from s in @, m in s->map[1], x in m[1], y in x[1]    // CHROMOSOME_I 
bql select s,m,x,y from s in @, m in s->map[2], x in m[1], y in x[1]    //  1
bql select s,m,x,y from s in @, m in s->map[3], x in m[1], y in x[1]    //  43054
bql select s,m,x,y from s in @, m in s->map[4], x in m[1], y in x[1]    // null


bql select s,m,x,y from s in @active:1, m in s->map,x in m[2],y in m[3]
bql select s,m,x,y,z from s in @, m in s->map, x in m[1], y in m[2],z in y[1]
bql select s,m,x,y,z from s in @, m in s->map, x in m[1], y in m[2],z in m[3]
bql select s,m,x,y,z from s in @, m in s->map, x in m[1], y in m[2],z in y[1]


// on human 37a_1 this code is much slower than the equivalent table-maker
query find mrna  "UBR4.aAug10" ; >cdna_clone
select c,r,c,m,a1,a2, x1,x2,len,ali,err,prod,c,inv,ano,ref_seq, tiling, hlib, gene from c in @active:1, m in c->in_mrna where exists m, tg in m->from_gene, gene in tg->gene  where gene like "UBR4" , r in c->read, t in r->tissue, a1 in m->constructed_from, a2 in a1[1], r1 in a2[1] where !r1 OR r1 == r, x1 in r1[1], x2 in x1[1], tg1 in r->from_gene where ! tg1 OR  tg1 == tg, len in tg1[1], ali in len[2], err in ali[1], inv in r->inverted[0] , ano in c->anomalous_clone[0], prod in m->product where ! prod OR prod#best_product,  ref_seq in r->ref_seq[0], tiling in r->mRNA_tiling[0], hlib in c->hinv_libs 

no_order


find MRNA G_t_14_0_1477.a
select m,s,e,n from m in @, s in m#splicing, e in s[5], n in s[8] where n < 1000  // OK

*/

#ifdef JUNK

/* example of a code checking the input more carefully
 *  provided by Nicolas.Thierry-Mieg@univ-grenoble-alpes.fr
 * 2017_10_12
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#define BUFFSIZE 64

/* 
   The string str should hold a positive long int value.
   If it does, return the value.
   Otherwise, print an error message to stdout and return -1.
*/
long get_long_from_str(char *str)
{
    char *endptr;
    long val;

    val = strtol(str, &endptr, 10);

    /* Check for various possible errors */
    if ((val == LONG_MAX) || (val < 0)) {
	printf("value is out of range or negative, try again\n");
	return (-1);
    } else if (endptr == str) {
	printf
	    ("Not a number (no leading digits), try again and type just a number\n");
	return (-1);
    } else if (*endptr != '\0') {
	printf
	    ("Not a number (trash after digits), try again and type just a number\n");
	return (-1);
    } else
	// AOK, we got a positive long int and nothing else
	return (val);
}

/*
  read lines from stdin until a line is entered that holds 
  a positive long int and nothing else.
  Complain whenever the line holds something else and loop.
  Return the positive long int value when the input is correct.
*/
long read_long_from_stdin(void)
{
    char buffer[BUFFSIZE];
    long value = -1;

    do {
	char *gets_retval = fgets(buffer, BUFFSIZE, stdin);
	if (gets_retval != buffer)
	    printf("some kind of read error, please try again\n");
	else if (buffer[strlen(buffer) - 1] != '\n') {
	    // flush remaining chars from stdin
	    while (getchar() != '\n') {
	    }
	    printf("too many characters, please try again\n");
	} else {
	    // remove final newline
	    buffer[strlen(buffer) - 1] = '\0';
	    value = get_long_from_str(buffer);
	}
    } while (value == -1);

    return (value);
}

int main(void)
{
    long value;
    printf("please type a positive long int\n");
    value = read_long_from_stdin();
    printf("OK, found a positive long int: %li\n", value);
    return (0);
}

#endif
