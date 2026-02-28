/*  File: table.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: table.c,v 1.19 2020/01/12 12:58:51 mieg Exp $
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 17 15:58 1999 (fw)
 * * Aug 26 18:19 1998 (fw): proper tableSort over multiple columns
 * * Aug 25 14:25 1998 (fw): tableOut now writes NULL for '0' column-values
 * * Aug 14 12:07 1998 (fw): fixed silly bug in tableCreateFromKeySet()
 * * Jul 15 16:48 1998 (fw): fixed silly bug in tableMax()
 * Created: Thu Oct 17 14:29:36 1996 (rd)
 *-------------------------------------------------------------------
 */

#include "table.h"

#include "acedb.h"
#include "bitset.h"
#include "java.h"
#include "bs.h"
#include "a.h"
#include "sysclass.h"
#include "lex.h"
#include "pick.h"
#include "flag.h"
#include "freeout.h"
#include "dict.h"

static magic_t TABLE_MAGIC = "TABLE";

static void tableFinalise (void *vp)
{
  TABLE *t = (TABLE*)vp;

  if (!t) return ;
  if (t->magic != &TABLE_MAGIC)
    messcrash ("tableFinalise() received non-magic TABLE* pointer") ;
  t->magic = 0 ;

  handleDestroy (t->handle) ;

  return;
} /* tableFinalise */


void uTableDestroy (TABLE *t)
     /* just here for completeness, messfree(table) is enough */
{ 
  messfree (t) ; 

  return;
} /* uTableDestroy */


static BOOL tableExists (TABLE *tt)
{
  if (tt && tt->magic == &TABLE_MAGIC && stackExists (tt->s))
    return TRUE;

  return FALSE;
} /* tableExists */


TABLE *tableHandleCreate (int size, const char *type, AC_HANDLE parent)
/* Note: type has to be 0-terminated string, 
   it is duplicated for the table object */
{ 
  TABLE *t ;
  int i, ncol ;
				/* first check type is OK */
  if (!type || !*type)
    messcrash ("tableCreate: received null or empty type string") ;

  for (i = 0 ; i < strlen(type); ++i)
    switch (type[i])
      { 
      case 'k': case 'g': case 'b': case 'i': case 'f': 
      case 's': case 't': case 'a': case '0': case 'D':
	break ;
      default:
	messcrash ("tableCreate: bad type %c for column %d", type[i], i+1) ;
      }

  if (sizeof(float) != sizeof (KEY) ||
      sizeof(int) != sizeof (KEY) ||
      sizeof(TABLETYPE) != sizeof (KEY) ||
      sizeof(mytime_t) != sizeof (KEY))
    messcrash 
      ("%s\n%s", "On this platform different TABLETYPEs have different byte size",
       " Rewrite the table.h pointer macros tabXxxp ! ") ;

  /* make the table and deal with handles */
  t = (TABLE *)halloc(sizeof(TABLE), parent) ;
  blockSetFinalise (t, tableFinalise) ;

				/* initialise TABLE */
  t->magic = &TABLE_MAGIC ;
  t->handle = handleCreate () ;
  t->ncol = ncol = strlen(type) ;
  t->type = strnew(type, t->handle) ;
  t->name = (TABLETYPE*) halloc(ncol*sizeof(int), t->handle) ;
  t->visibility = (BOOL*) halloc(ncol*sizeof(int), t->handle) ;
  t->s = 0 ;
  t->col = (Array*) halloc(ncol*sizeof(Array), t->handle) ;
  t->parents = (KEYSET*) halloc(ncol*sizeof(KEYSET), t->handle);
  t->grandParents = (KEYSET*) halloc(ncol*sizeof(KEYSET), t->handle);
  t->emptyBits = (BitSet*) halloc(ncol*sizeof(BitSet), t->handle) ;
  t->s = stackHandleCreate (4096, t->handle) ;
  t->s->textOnly = TRUE;
  pushText (t->s, type) ; /* ! occupy mark zero in  a most useful way */
  for (i = 0 ; i < ncol ; ++i)
    { t->name[i].i = stackMark (t->s) ; pushText (t->s, messprintf("f%d",i+1)) ; }
  
  if (size < 64)
    size = 64 ;
  /* initialise each column */
  for (i = 0 ; i < ncol ; ++i)
    {
      /* make each character in the type-string lower-case
       * but preserve the visibility info given by upper-case
       * chars in an extra array */
      char c1 = t->type[i], c2 ;
      c2 = ace_lower (c1) ;
      if (c2 != c1)
	{
	  t->type[i] = c2 ;
	  t->visibility[i] = FALSE;
	}
      else
	t->visibility[i] = TRUE;

      t->name[i].i = stackMark (t->s) ;
      pushText (t->s, "") ;

      t->col[i] = arrayHandleCreate (size, TABLETYPE, t->handle) ;
      t->parents[i] = arrayHandleCreate(size, KEY, t->handle);
      t->grandParents[i] = arrayHandleCreate(size, KEY, t->handle);
      t->emptyBits[i] = bitSetCreate (size, t->handle) ;
    }

  if (!tableExists(t))
    messcrash("tableHandleCreate() - just created invalid TABLE\n"
	      "  something seriously screwy - look at tableExists()");

  return t ;
} /* tableHandleCreate */

TABLE *tableChangeHandle (TABLE *t, AC_HANDLE parent)
     /* XXXXXXXXXX probably bugged XXXXXXXXXXXX */
{
  TABLE *tt = (TABLE *)halloc(sizeof(TABLE), parent) ;
  memcpy(tt,t,sizeof(TABLE)) ;
  blockSetFinalise (t,0) ; /* thus t->handle will not be destroyed */

  return tt ;
} /* tableChangeHandle */

/***********************************/

TABLE *tableCreateFromKeySet (KEYSET keySet, AC_HANDLE handle)
/* make a table of KEYs from a given KEYSET (fw-980805) */
{
  TABLE *table;
  int row;

  if (!keySet)
    return 0;

  table = tableHandleCreate (keySetMax(keySet), "k", handle);

  if (keySetMax(keySet))
    {
      tableKey(table, keySetMax(keySet)-1, 0) = 0; /* set last item to extend size */
      for (row = 0; row < keySetMax(keySet); row++)
	tabKey(table, row, 0) = arr(keySet, row, KEY);
    }

  return table;
} /* tableCreateFromKeySet */


KEYSET  keySetCreateFromTable (TABLE *table, AC_HANDLE handle)
/* makes a keyset from the KEY-type table column (fw-980805) */
{
  int c, i, keyCol=-1;
  KEYSET newKeySet;

  for (c = 0; c < table->ncol; ++c)
    if (table->type[c] == 'k') keyCol = c;

  if (keyCol == -1)
    return 0;		/* table contains no column of type KEY */

  newKeySet = keySetHandleCreate(handle);

  if (tableMax(table))
    {      
      array(newKeySet, tabMax(table,keyCol)-1, KEY) = 0; /* set last item to extend size */
      for (i = 0 ; i < tabMax(table,keyCol) ; ++i)
	arr(newKeySet,i,KEY) = tabKey(table,i,keyCol) ;
    }

  keySetSort(newKeySet) ;
  keySetCompress(newKeySet) ;

  return newKeySet;
} /* keySetCreateFromTable */

/************************************************************/

TABLE *tableCopy (TABLE *old, AC_HANDLE handle)
{
  TABLE *new;
  int colNum;

  if (!tableExists (old))
    messcrash("tableCopy() - received invalid TABLE* pointer");

  new = tableHandleCreate(tableMax(old), old->type, handle);

  /* kill the empty arrays that have been initialised and replace
   * them with copies of the data in the old table */

  stackDestroy(new->s);
  new->s = stackCopy(old->s, new->handle);

  for (colNum = 0; colNum < old->ncol; colNum++)
    {
      arrayDestroy(new->col[colNum]);
      new->col[colNum] = arrayHandleCopy(old->col[colNum], new->handle);

      arrayDestroy(new->parents[colNum]);
      new->parents[colNum] = arrayHandleCopy(old->parents[colNum], new->handle);

      arrayDestroy(new->col[colNum]);
      new->col[colNum] = arrayHandleCopy(old->col[colNum], new->handle);

      bitSetDestroy(new->emptyBits[colNum]);
      new->emptyBits[colNum] = bitSetCopy(old->emptyBits[colNum], new->handle);

      new->name[colNum].i = old->name[colNum].i;
      new->visibility[colNum] = old->visibility[colNum];
    }

  return new;
} /* tableCopy */

/***********************************/

const char *tableTypes (TABLE *t)
{
  return t ? t->type : "" ;
}


const char *tableGetColumnName (TABLE *t, int i)
{
  return t->name[i].i ? stackText (t->s, t->name[i].i) : 0 ;
}

const char *tableSetColumnName (TABLE *t, int i, const char *name)
{ 
  const char *old ;

  if (t->magic != &TABLE_MAGIC)
    messcrash ("tableSetColumnName() received non-magic TABLE* pointer");
  if (i < 0 || i >= t->ncol)
    messcrash ("tableSetColumnName() - col %d is out of bounds 0-%d",
	       i, t->ncol-1) ;
  if (!name || !*name)
    messcrash ("tableSetColumnName() - NULL or empty name") ;

  old = t->name[i].i ? stackText (t->s, t->name[i].i) : 0 ;
  t->name[i].i = stackMark(t->s) ;
  pushText (t->s, name) ;

  return old ;		
} /* tableSetColumnName */

/***********************************/

void tableSetString (TABLE *t, int i, int j, char *string)
{ 
  if (t->magic != &TABLE_MAGIC)
    messcrash ("tableSetString() received non-magic TABLE* pointer");
  if (!stackExists(t->s))
    messcrash ("tableSetString() received TABLE* object with invalid Stack ->s pointer") ;
  if (string && *string)
    {
      tableInt(t,i,j) = stackMark(t->s) ;
      pushText (t->s, string) ;
    }
  else 
    tableInt(t,i,j) = 0 ;

  tableUnSetEmpty (t, i, j);	/* mark table field as non-NULL */

  return;
} /* tableSetString */

/***********************************/

int tableMax (TABLE *t)
/* return the largest number of rows over all columns */
{ 
  int colNum, j = 0, max = 0 ;

  if (!t) return -1;

  if (t->magic != &TABLE_MAGIC)
    messcrash ("tableMax() received non-magic TABLE* pointer");

  for (colNum = 0; colNum < t->ncol; colNum++)
    { 
      j = arrayMax(t->col[colNum]) ; 
      if (j > max) 
	max = j ; 
    }

  return max ;
} /* tableMax */

/***********************************/
   /* no check, could possibly be a macro */
BOOL tabValueNull (TABLE *t, int i, int j)
{ 
  if (i > tabMax(t,j))   
    return TRUE ;
  switch (t->type[j])
    {
    case 'k': case 'g': case 's': case 't': case 'D':
      if (!tabKey(t,i,j)) 
	return TRUE ;
    case 'b': 
       if (tabBool(t,i,j) == UT_NON_BOOL)
	 return TRUE ;
      break ;
    case 'a': 
       if (tabFlag(t,i,j) == UT_NON_FLAG)
	 return TRUE ;
      break ;
    case 'i': 
       if (tabInt(t,i,j) == UT_NON_INT)
	 return TRUE ;
      break ;
    case 'f': 
       if (tabFloat(t,i,j) == UT_NON_FLOAT)
	 return TRUE ;
      break ;
    }

  return FALSE ;
} /* tabValueNull */

/*************/

void tableCheckEmpty (TABLE *t, int j)
{
  int i ;

  for (i = 0 ; i < tabMax(t,j) ; ++i)
    if (tabValueNull (t, i, j))
      tableSetEmpty(t,i,j) ;
    else
      tableUnSetEmpty(t,i,j) ;

  return;
} /* tableCheckEmpty */

/**********************************************************************/

void tableClear (TABLE *t, int n)
{
  int i ;

  for (i = 0 ; i < t->ncol ; ++i)
    { t->col[i] = arrayReCreate (t->col[i], n, TABLETYPE) ;
      t->parents[i] = arrayReCreate (t->parents[i], n, KEY) ;
      t->grandParents[i] = arrayReCreate (t->grandParents[i], n, KEY) ;
      t->emptyBits[i] = bitSetReCreate (t->emptyBits[i], n) ;
    }
}

BOOL tableSetTypes (TABLE *t, const char *types)
/* can only re-set '0' typed columns
   returns FALSE if one or more columns already had a type */
{
  int i ;

  for (i = 0 ; i < t->ncol ; ++i) /* check first */
    if (types[i] != '0' && types[i] != t->type[i] && t->type[i] != '0')
      return FALSE ;

  for (i = 0 ; i < t->ncol ; ++i)
    if (types[i] != '0' && types[i] != t->type[i])
      t->type[i] = types[i] ;

  strcpy (stackText(t->s, 0), t->type) ;

  return TRUE ;
} /* tableSetTypes */

/*********************************************************/

static TABLE *sortTable = 0 ;
static Array sortSpecList = 0;

static int tableCompare (const void* va, const void* vb)
/* comparison function used for qsort() in tableSort,
   that also rearranges an index array of row numbers,
   requires static sortSpec and sortTable to be assigned */
{
  register int a = *(int*)va ;
  register int b = *(int*)vb ;
  register int i, sortCol, sortSpec;


  for (i = 0; i < arrayMax(sortSpecList); ++i)
    {
      sortSpec = arr(sortSpecList,i,int);

      /* get the true sort column in C-numbers */
      if (sortSpec < 0)	/* descending */
	sortCol = -sortSpec - 1;
      else /* sortSpec > 0   -> ascending  */
	sortCol = sortSpec - 1;

      /* we make NULL the smallest value ever */
      if (tabEmpty(sortTable, a, sortCol))
	{
	  if (tabEmpty(sortTable, b, sortCol))
	    continue ;
	  return sortSpec>0 ? -1 : 1 ;
	}
      
      if (tabEmpty(sortTable, b, sortCol))
	return sortSpec>0 ? 1 : -1 ;

      if (tabKey(sortTable, a, sortCol) == tabKey(sortTable, b, sortCol))
	continue ;

      switch (sortTable->type[sortCol])
	{
	case 'k': case 'g': 
	  {
	    KEY
	      xk = tabKey(sortTable, a, sortCol),
	      yk = tabKey(sortTable, b, sortCol) ;
	    int r = keySetAlphaOrder(&xk, &yk);
	    if (r != 0)
	      return sortSpec>0 ? r : -r ;
	  }
	case 'b':
	  {
	    BOOL 
	      ba = tabBool(sortTable, a, sortCol),
	      bb = tabBool(sortTable, b, sortCol) ;
	    int r = ba ? (bb ? -1 : 0 ) : ( bb ? 0 : 1 );
	    if (r != 0)
	      return sortSpec>0 ? r : -r ;
	  }
	case 'i':
	  {
	    int 
	      ba = tabInt(sortTable, a, sortCol),
	      bb = tabInt(sortTable, b, sortCol) ;
	    int r = ba - bb;
	    if (r != 0)
	      return sortSpec>0 ? r : -r ;
	  }
	case 'f':
	  {
	    float 
	      ba = tabFloat(sortTable, a, sortCol),
	      bb = tabFloat(sortTable, b, sortCol) ;
	    float r = ba - bb;
	    if (r != 0)
	      return sortSpec>0 ? r : -r ;
	  }
	case 's': case 'D':
	  {
	    char 
	      *ba = tabString(sortTable, a, sortCol),
	      *bb = tabString(sortTable, b, sortCol) ;
	    int r = lexstrcmp(ba,bb) ;
	    if (r != 0)
	      return sortSpec>0 ? r : -r ;
	  }
	case 't':
	  {
	    mytime_t
	      ba = tabDate(sortTable, a, sortCol),
	      bb = tabDate(sortTable, b, sortCol) ;
	    int r = ba - bb;
	    if (r != 0)
	      return sortSpec>0 ? r : -r ;
	  }
	} /* end switch */
    }

  /* default keep original order */
  if (a > b) 
    return +1 ;

  return -1 ;
} /* tableCompare */
/************************************************************/

BOOL tableSort (TABLE *t, const char *spec)
/* sorts elements in the table according to
   the sort specifier :
   spec should look like "-2+4+1" which would specify first priority
   column 2 reversed, then column 4, then column 1.

   returns FALSE if table is wrong, or spec is out of range
*/
{
  AC_HANDLE scopeHandle;
  int i, i1, numRows, specLen, sortCol;
  int *indexArray ;
  const char *ccp ;
  BOOL isPlus = TRUE ;

  if (!t) return FALSE;
  if (t->magic != &TABLE_MAGIC)
    messcrash ("tableSort() received non-magic TABLE* pointer") ;

  if (!spec)
    messcrash ("tableSort() received NULL sort specifier");

  specLen = strlen(spec);

  sortSpecList = arrayReCreate(sortSpecList, (specLen), int);
  for (i = 0, ccp = spec; *ccp ; i++)
    {
      isPlus = TRUE ;
      if (*ccp == '+') ccp++ ;
      else if (*ccp == '-') { isPlus = FALSE ; ccp++ ; } ;
      sortCol = 0 ;
      while (*ccp && *ccp >= '0' && *ccp <= '9')
	{
	  i1 = *ccp++ - '0' ;
	  sortCol = 10 * sortCol + i1 ;
	}
      if (!sortCol) break ;
      if (sortCol < 1 || sortCol > t->ncol) 
	/* column numbers have to range from 1 to colNum */
	return FALSE;
      if (isPlus)
	array(sortSpecList, i, int) = sortCol;
      else 
	array(sortSpecList, i, int) = -sortCol;
    }

  numRows = tableMax(t);
  if (numRows <= 1) return TRUE; /* 1 or no elements are always sorted */

  /***************************************/

  scopeHandle = handleCreate();

  /* create the indices that keep track of which row went where
     during comparative sorting */
  indexArray = (int*)halloc (numRows * sizeof(int), scopeHandle);

  /* fill the indices with row numbers */
  for (i = 0; i < numRows; ++i)
    indexArray[i] = i;

  /* set global table-var used by compare() to be this one */
  sortTable = t;

  /* rearrange lines */
  qsort (indexArray, numRows, sizeof(int), tableCompare) ;
  
  /**** use index to rearrange column by column ****/

#ifdef NEW_RICHARD_CODE
  {
    TABLETYPE tmp ;
    BitSet bs=0 ;
    int k, j;

    for (i = 0 ; i < t->ncol ; ++i)
      { 
	bs = bitSetReCreate (bs, tabMax(t,i));
	for (j = 0 ; j < tabMax(t,i) ; ++j)
	  if (!bit(bs,j))         /* not done yet */
	    if (indexArray[j] != j)
	      { 
		tmp = arr((t)->col[i],j,TABLETYPE) ;
		for (k = j ;
		     !bitSet(bs, indexArray[k]) ;
		     k = indexArray[k])
		  { 
		    arr((t)->col[i],k,TABLETYPE) = 
		      arr((t)->col[i], indexArray[k],TABLETYPE) ;
		    bitSet(bs,k) ;
		  }
		arr((t)->col[i],k,TABLETYPE) = tmp ;
		bitSet(bs,k) ;
	      }
	    else
	      bitSet(bs,j) ;
      }
    
    bitSetDestroy (bs);
  }
#else
 /****** still use previuos code, until we get the new one to work ****/
  {
    Array aa, dummy ;
    BitSet bb, oldbb ;
    int n = numRows, ii, j, j1;
    
    for (ii = 0 ; ii < t->ncol ; ii++)
      {
	dummy = t->col[ii] ;
	aa = t->col[ii] = arrayHandleCreate (n, TABLETYPE, t->handle) ;
	array(aa,n-1,TABLETYPE).k = 0 ; /* make room */
	oldbb = t->emptyBits[ii] ;
	bb =  t->emptyBits[ii] = bitSetCreate (n, t->handle) ;
	for (j=0; j < n ; j++)
	  {
	    j1 = indexArray[j] ;
	    arr(aa,j,TABLETYPE).k = array(dummy,j1,TABLETYPE).k ;
	    if (bit(oldbb,j1)) bitSet(bb,j) ;
	  }
	arrayDestroy(dummy) ;
	bitSetDestroy (oldbb) ;
      }
  }
#endif /* !NEW_RICHARD_CODE */
  handleDestroy (scopeHandle);

  return TRUE;
} /* tableSort */
/************************************************************/

void tableMakeUnique (TABLE *t) 	/* assumes sorted */
{
  int i1, ii, j , mx = tableMax(t), nCol = t->ncol ; ;
  
  if (mx < 2 || !nCol) return ;

  for (i1 = 0, ii=1 ; ii < mx ; ii++ )
    {
      /* dry pass */
      for (j=0; j<nCol ; j++)
	if(tabKey(t,ii,j) != tabKey(t,i1,j)) goto different ;
      continue ; /* line ii == line i1, drop line ii */
    different:
      i1++ ;   /* accept line ii */
      if (i1 < ii)  /* copy if needed */
	{
	  for (j=0; j<nCol ; j++) 
	    {
	      tabKey(t,i1,j) = tabKey(t,ii,j) ;
	      if (tabEmpty(t,ii,j))
		bitSet(t->emptyBits[j],i1) ;
	      else
		bitUnSet(t->emptyBits[j],i1) ;
	    }
	}
    }
  for (j=0; j<nCol ; j++)
    arrayMax(t->col[j]) = i1 + 1 ;

  return;
} /* tableMakeUnique */
/************************************************************/

/************************************************/
/******* table operator utility functions *******/
/************************************************/

/* Function to check whether two tables are compatible, that
   they can be used for the subsequent table operators together. */
/*
   NOTE: Once this check is performed on both tables the short version
   tabXXX() of the operations can be used. In case it is only a one-off
   statement that uses a table operation, the tableXXXX() routines 
   should be used to include the type checking from this function
*/
BOOL tableIsColumnCompatible (TABLE *leftTable, 
			      TABLE *rightTable)
{
  int colCount;

  if (leftTable == NULL || rightTable == NULL)
    return FALSE;

  /* we check that the value-types of the columns match in both tables */
  if (strcmp(leftTable->type, rightTable->type) != 0)
    return FALSE;

  /* check for unevaluated columns, ideally we'd just like to exclude them */
  for (colCount = 0; colCount < leftTable->ncol; colCount++)
    {
      if (leftTable->type[colCount] == '0' || rightTable->type[colCount] == '0')
	return FALSE;
    }

  return TRUE;
} /* tableIsColumnCompatible */

/***********************************************************************************/

/* function to compare a complete row in a table with a row in a second table. */
/*
   Returns FALSE if tables are incompatible, type mismatch or uninitialised
   ONLY returns TRUE if full compatability and if all values in the two rows match
*/
BOOL tableIsRowEqual (int    rowCountLeft, 
		      TABLE *leftTable, 
		      int    rowCountRight, 
		      TABLE *rightTable)
{
  if (!tableIsColumnCompatible(leftTable, rightTable))
    return FALSE;

  return (tabIsRowEqual(rowCountLeft, leftTable, rowCountRight, rightTable));
} /* tableIsRowEqual */

/*************************************************************************************/


/* simple function to compare a complete row in a table with a row in a second table. */
/*
   NOTE: no type checking for column compatability is performed
   this function is to be used for quick access in long loops,
   for sporadic occasional row-comparison that includes all type checking
   please use tableIsRowEqual().
*/
BOOL tabIsRowEqual (int    rowCountLeft, 
		    TABLE *leftTable, 
		    int    rowCountRight, 
		    TABLE *rightTable)
{
  int  colCount;
  BOOL isValueDifferent = FALSE;

  for (colCount = 0; colCount < leftTable->ncol; colCount++)
    {
      switch (leftTable->type[colCount])
	{
	case 'i':
	  if (tabInt(leftTable, rowCountLeft, colCount) != tabInt(rightTable, rowCountRight, colCount))
	    isValueDifferent = TRUE;
	  break;
	case 'f':
	  if (tabFloat(leftTable, rowCountLeft, colCount) != tabFloat(rightTable, rowCountRight, colCount))
	    isValueDifferent = TRUE;
	  break;
	case 'b':
	  if (tabBool(leftTable, rowCountLeft, colCount) != tabBool(rightTable, rowCountRight, colCount))
	    isValueDifferent = TRUE;
	  break;
	case 't':
	  if (tabDate(leftTable, rowCountLeft, colCount) != tabDate(rightTable, rowCountRight, colCount))
	    isValueDifferent = TRUE;
	  break;
	case 'k':		/* keys */
	case 'g':		/* tags */
	  if (tabKey(leftTable, rowCountLeft, colCount) != tabKey(rightTable, rowCountRight, colCount))
	    isValueDifferent = TRUE;
	  break;
	case 's': case 'D':
	  if (strcmp(tabString(leftTable, rowCountLeft, colCount), tabString(rightTable, rowCountRight, colCount)) != 0)
	    isValueDifferent = TRUE;
	  break;
	}
    }
  if (!isValueDifferent)
    return TRUE;

  return FALSE;
} /* tabIsRowEqual */

/**********************************************************************/

AC_HANDLE globalTableCreateHandle;

/* function to copy all values in a row of the srcTable to the end of the destTable */
/*
   The destTable is created if it is not existent.
   The function returns FALSE, if an existing destTable is not type-compatible
   with the source table.
   Otherwise, the values in the specified row are appended to the end
   of the destTable and TRUE is returned.
*/
BOOL tableCopyRow (int    rowCountSrc, 
		   TABLE *srcTable, 
		   TABLE *destTable)
{
  if (!destTable)
    {
      /* we have to initialise the destTable to be of type and value of srcTable */
      destTable = tableHandleCreate (tableMax(srcTable), tableTypes(srcTable), globalTableCreateHandle) ;
    }
  else
    {
      /* a destination table exists, but it is compatible ? */
      if (!tableIsColumnCompatible(srcTable, destTable))
	return FALSE;
    }

  tabCopyRow (rowCountSrc, srcTable, destTable);

  return TRUE;
} /* tableCopyRow */

/**********************************************************************/

/* simple function to copy all values in a row of the srcTable to the end of the destTable */
/*
   NOTE: 
   (1) No type checking for column compatability is performed.
   (2) The destTable is assumed to be initialised with the same types as the srcTable.

   This function is to be used for quick access in long loops,
   for sporadic occasional row-copy that includes all type checking
   and initialisation please use tableCopyRow().
*/
void tabCopyRow (int    rowCountSrc, 
		 TABLE *srcTable, 
		 TABLE *destTable)
{
  int colCount;
  int rowCountDest = tableMax(destTable);

  for (colCount = 0; colCount < srcTable->ncol; colCount++)
    {
      switch (srcTable->type[colCount])
	{
	case 'i':
	  tableInt(destTable, rowCountDest, colCount) = tabInt(srcTable, rowCountSrc, colCount);
	  break;
	case 'f':
	  tableFloat(destTable, rowCountDest, colCount) = tabFloat(srcTable, rowCountSrc, colCount);
	  break;
	case 'b':
	  tableBool(destTable, rowCountDest, colCount) = tabBool(srcTable, rowCountSrc, colCount);
	  break;
	case 't':
	  tableDate(destTable, rowCountDest, colCount) = tabDate(srcTable, rowCountSrc, colCount);
	  break;
	case 'k':
	  tableKey(destTable, rowCountDest, colCount) = tabKey(srcTable, rowCountSrc, colCount);
	  break;
	case 's': case 'D':
	  tableSetString(destTable, rowCountDest, colCount, tabString(srcTable, rowCountSrc, colCount));
	  break;
	case '0':
	  tableSetEmpty(destTable, rowCountDest, colCount);
	  break;

	}
    }
} /* tabCopyRow */
/************************************************************/

/* produce a combined table that contains rows from both
   tables but no duplicates */
/* fails with return code FALSE, if the two tables are incombinable
   or in case an existing table pointer was passed, this has
   already been assigned column types, otherwise a new table
   will be returned in *destTable 
*/
BOOL tableUnion (TABLE *leftTable,
		 TABLE *rightTable,
		 TABLE **destTable)
{
  int rowCountLeft, rowCountRight;
  int lengthLeft, lengthRight;
  BOOL isInLeftTable;

  if (!destTable)	   /* has to be called with a pointer to TABLE* 
			      in order to return results */
    messcrash ("table.c:tableUnion() - destTable pointer to TABLE* is NULL");


  if (!tableIsColumnCompatible(leftTable, rightTable))
    return FALSE;

  lengthLeft = tableMax(leftTable);
  lengthRight = tableMax(rightTable);

  if (!*destTable)		/* table it points to is still NULL */
    {
      *destTable = tableHandleCreate((lengthLeft + lengthRight/2),
				     leftTable->type, globalTableCreateHandle);
    }
  else				/* we got passed a pointer of an existing table */
    {
      if (strlen(leftTable->type) != strlen((*destTable)->type))
	messcrash ("table.c:tableUnion() - previously created destTable has incorrect column number");

      if (!tableSetTypes(*destTable, leftTable->type))
	return FALSE;		/* fails if types were already set */
    }

  /********/
  /* include all rows from leftTable and the ones from rightTable we haven't got already */
  /********/

  /* copy all rows from the left table */
  for (rowCountLeft = 0; rowCountLeft < lengthLeft; rowCountLeft++)
    {
      tabCopyRow(rowCountLeft, leftTable, *destTable);
    }
  
  /* only copy rows from right table that we haven't got already */
  for (rowCountRight = 0; rowCountRight < lengthRight; rowCountRight++)
    {
      /* look if we have this row already */
      isInLeftTable = FALSE;
      for (rowCountLeft = 0; rowCountLeft < lengthLeft; rowCountLeft++)
	if (tabIsRowEqual(rowCountLeft, leftTable,
			  rowCountRight, rightTable))
	  {
	    isInLeftTable = TRUE;
	    break;
	  }
      if (!isInLeftTable)
	tabCopyRow(rowCountRight, rightTable, *destTable);
    }


  return TRUE;
} /* tableUnion */
/************************************************************/



/* produce a table with rows that are in the first but 
   not the second table */
/* fails with return code FALSE, if the two tables are incombinable
   or in case an existing table pointer was passed, this has
   already been assigned column types, otherwise a new table
   will be returned in *destTable 
*/
BOOL tableDiff (TABLE *leftTable,
		TABLE *rightTable,
		TABLE **destTable)
{
  int rowCountLeft, rowCountRight;
  int lengthLeft, lengthRight;
  BOOL isInRightTable;

  if (!destTable)	   /* has to be called with a pointer to TABLE* 
			      in order to return results */
    messcrash ("table.c:tableDiff() - destTable pointer to TABLE* is NULL");


  if (!tableIsColumnCompatible(leftTable, rightTable))
    return FALSE;

  lengthLeft = tableMax(leftTable);
  lengthRight = tableMax(rightTable);

  if (!*destTable)		/* table it points to is still NULL */
    {
      *destTable = tableHandleCreate((lengthLeft - lengthRight/2),
				     leftTable->type, globalTableCreateHandle);
    }
  else				/* we got passed a pointer of an existing table */
    {
      if (strlen(leftTable->type) != strlen((*destTable)->type))
	messcrash ("table.c:tableDiff() - previously created destTable has incorrect column number");

      if (!tableSetTypes(*destTable, leftTable->type))
	return FALSE;		/* fails if types were already set */
    }

  /********/
  /* only include rows from the leftTable that aren't in the rightTable */
  /********/

  for (rowCountLeft = 0; rowCountLeft < lengthLeft; rowCountLeft++)
    {
      /* look for this row in the rightTable */
      isInRightTable = FALSE;
      for (rowCountRight = 0; rowCountRight < lengthRight; rowCountRight++)
	if (tabIsRowEqual(rowCountLeft, leftTable,
			  rowCountRight, rightTable))
	  {
	    isInRightTable = TRUE;
	    break;
	  }
      if (!isInRightTable)
	/* not found -> copy left row into destinationTable */
	tabCopyRow(rowCountLeft, leftTable, *destTable);
    }

  return TRUE;
} /* tableDiff */
/************************************************************/



/* produce a table with rows that are in the first but 
   not the second table */
/* fails with return code FALSE, if the two tables are incombinable
   or in case an existing table pointer was passed, this has
   already been assigned column types, otherwise a new table
   will be returned in *destTable 
*/
BOOL tableIntersect (TABLE *leftTable,
		     TABLE *rightTable,
		     TABLE **destTable)
{
  int rowCountLeft, rowCountRight;
  int lengthLeft, lengthRight;
  BOOL isInRightTable;

  if (!destTable)	   /* has to be called with a pointer to TABLE* 
			      in order to return results */
    messcrash ("table.c:tableIntersect() - destTable pointer to TABLE* is NULL");


  if (!tableIsColumnCompatible(leftTable, rightTable))
    return FALSE;

  lengthLeft = tableMax(leftTable);
  lengthRight = tableMax(rightTable);

  if (!*destTable)		/* table it points to is still NULL */
    {
      *destTable = tableHandleCreate((lengthLeft/2 + lengthRight/2),
				     leftTable->type, globalTableCreateHandle);
    }
  else				/* we got passed a pointer of an existing table */
    {
      if (strlen(leftTable->type) != strlen((*destTable)->type))
	messcrash ("table.c:tableIntersect() - previously created destTable has incorrect column number");

      if (!tableSetTypes(*destTable, leftTable->type))
	return FALSE;		/* fails if types were already set */
    }

  /********/
  /* only include rows that are shared between leftTable and rightTable */
  /********/

  for (rowCountLeft = 0; rowCountLeft < lengthLeft; rowCountLeft++)
    {
      /* look for this row in the rightTable */
      isInRightTable = FALSE;
      for (rowCountRight = 0; rowCountRight < lengthRight; rowCountRight++)
	if (tabIsRowEqual(rowCountLeft, leftTable,
			  rowCountRight, rightTable))
	  {
	    isInRightTable = TRUE;
	    break;
	  }
      if (isInRightTable)
	/* left was found in right -> copy left row into destinationTable */
	tabCopyRow(rowCountLeft, leftTable, *destTable);
    }

  return TRUE;
} /* tableIntersect */
/************************************************************/


/************************************************************/
/********** functions to output/export tables ***************/
/************************************************************/

int tableOut (TABLE *t, char separator, char style)
{ 
  return tableSliceOut (0, tableMax(t), t, separator, style) ;
} /* tableOut */
/************************************************************/

int tableSliceOut (int begin, 
		   int count, 
		   TABLE *t, 
		   char separator, 
		   char style)
{
  int i, n = 0, nMax, end;
  char sepString[2] ;
  BOOL showEmpty = FALSE ;
  int ii = -1 ;
  char timeBuf[25] ;
  char *cp ; int xi ; float xf ; mytime_t xt ; KEY xk ; 
  Array oldClass = 0 ; 
  DICT **colFlagDict = 0 ;
  const char *ccp ;
  const char **colFlagSet = 0 ;
  char bufvoid[1] = {'v'} ;
  char bufinteger[1] = {'i'} ;
  char buffloat[1] = {'f'} ;
  char bufdate[1] = {'d'} ;
  char buftext[1] = {'t'} ;
  char bufk[2] = {'k', 0} ;
  char bufK[2] = {'K', 0} ;
  char buftag[1] = {'g'} ;
  char beginline[1] ={'>'};

  if (style == 'a' || style == 'j' ||  style == 'J' || style == 'C' )
    showEmpty = TRUE ; 

  if (!t) return 0;
  if (t->magic != &TABLE_MAGIC)
    messcrash ("tableSliceOut() received non-magic TABLE* pointer") ;
  if (t->ncol <= 0) return 0 ;

  nMax = tableMax(t);

  /* limit the begin-value to be within bounds */
  begin = (begin < 0) ? 0 : begin;
  begin = (begin > nMax-1) ? nMax-1 : begin;
  /* limit the number of rows to be within bounds */
  count = (count < 0) ? 0 : count;
  /* calc the end of the table */
  end = begin + count;
  end = (end > nMax) ? nMax : end;


  oldClass = arrayCreate (t->ncol, int) ;
  freeOutInit() ;		/* just in case someone forgot */

  sepString[0] = separator ;
  sepString[1] = 0 ;
				/* preparation for flags */
  for (i = 0 ; i < t->ncol ; ++i) 
    /*
      RD 000615 - this looks wrong to me
      should always be -1 so you write class first time 
      so I am commenting it out
      
      mieg oct 2, 2004
      removes the commenting out
      we do not want the class name on table exports
    */
    if ((t->type[i] == 'k' || t->type[i] == 'g') && 
	tabMax(t,i) && !tabEmpty(t,0,i))
      array (oldClass, i, int) = class(tabKey(t,0,i));
    else
      array (oldClass, i, int) = -1 ;
  
  for (i = 0 ; i < t->ncol ; ++i) 
    if (t->type[i] == 'a')
      { if (!colFlagDict)
	  { colFlagDict = (DICT**) messalloc (t->ncol * sizeof(DICT*)) ;
	    colFlagSet = (const char**) messalloc (t->ncol * sizeof(char*)) ;
	  }
	if (!(ccp = tableGetColumnName(t,i)))
	  ccp = "Default" ;
	if (!strncmp (ccp, "Flag", 4) ||
	    !strncmp (ccp, "flag", 4))
	  { ccp += 4 ;
	    if (*ccp == 's') ++ccp ; /* remove "flags" as well as "flag" */
	    if (!*ccp)
	      ccp = "Default" ;
	  }
	colFlagSet[i] = ccp ;
	colFlagDict[i] = flagSetDict (ccp) ;
      }
   	
  if(style=='C')
    {
      xk=tabKey (t,begin,0) ;
      beginline[0] = '>' ;
    }

  ii = -1 ;
  for (n = begin ; n < end ; n++)
    {
      for (i = 0 ; i < t->ncol ; ++i)
	{
   	  if(style == 'C')
	    {
	      freeOutBinary (beginline, 1);
	      beginline[0]='>';
	    }
	  else if (i > 0)
	    freeOut (sepString) ;

	  if (n >= tabMax(t,i) || tabEmpty(t,n,i) || t->type[i] == '0')
	    {
	      if (showEmpty) 
		switch (style)
		  {
		  case 'j':  case 'J':
		    freeOut("???") ;
		    break ;
		  case 'C':
		    freeOutBinary (bufvoid, 1) ;
		    break;
		  default:
		    freeOut("NULL") ;
		    break ;
		  }
	    }
	  else
	    { 
	      switch (t->type[i])
		{ 
		case 'k':
		case 'g':
		  xk = tabKey(t,n,i) ;
		  switch (style)
		    {
		    case 'h': 
		      freeOut (name(xk)) ;
		      break ;
		    case 'C': 
		      {
			switch (t->type[i])
			  {
			  case 'g':
			    freeOutBinary (buftag,1) ;
			    freeOutBinary ((char*)&xk, 4) ;
			    break ;
			  case 'k':
			    if (iskey (xk) == 2)
			      { bufK[1] = class (xk) ; freeOutBinary (bufK,2) ; }
			    else
			      { bufk[1] = class (xk) ; freeOutBinary (bufk,2) ; }
			    cp = name(xk) ;
			    freeOutBinary (cp, strlen(cp) + 1) ;
			    break ;
			  }
			}
		      break ;
		    case 'j': 
		      freeOutf("?%s?%s?", className(xk), 
			       freejavaprotect(name(xk))) ;
		      break ;
		    case 'J': 
		      if (i <= ii && class(xk) == array (oldClass, i, int)) 
			freeOutf("?%s?", 
				 freejavaprotect(name(xk))) ;
		      else
			freeOutf("?%s?%s?", className(xk), 
				 freejavaprotect(name(xk))) ;
		      array (oldClass, i, int) = class(xk) ;
		      break ;
		    case 'a': 
		      if (class(xk) == array (oldClass, i, int)) 
			{
			  freeOut (freeprotect(name(xk))) ;
			  break ;
			}
		      array (oldClass, i, int) = class(xk) ;
		      /* else fall through */			
		    case 'A': 
		      freeOutf("%s:%s", className(xk), 
			       freeprotect(name(xk))) ;
		      break ;
		    default:
		      freeOut (freeprotect(name(xk))) ;
		      break ;
		}
		  break ;
		case 'a':
		  cp = flagNameFromDict(colFlagDict[i], tabFlag(t,n,i)) ;
		  switch (style)
		    {
		    case 'j': 
		      freeOutf ("?flag%s?%s?", colFlagSet[i], cp) ;
		      break ;
		    case 'J': 
		      if (i > ii) freeOutf ("tag") ;
		      freeOutf("?%s?", cp) ;
		      break ;
		    case 'C': 
		      freeOutBinary (buftext,1) ;
		      freeOutBinary (cp, strlen(cp) + 1) ;
		      break ;
		    default:
		      freeOut (cp) ;
		      break ;
		    }
		  break ;
		case 'b':
		  cp = tabBool(t,n,i) ? "TRUE" : "FALSE" ;
		  switch (style)
		    {
		    case 'j': 
		      freeOutf ("?BOOL?%s?", cp) ;
		      break ;
		    case 'J': 
		      if (i > ii) freeOutf ("BOOL") ;
		      freeOutf("?%s?", cp) ;
		      break ;
		    case 'C':
		      freeOutBinary (bufinteger,1) ;
		      if (tabBool(t,n,i))
		        xi = 1 ;
		      else
			xi = 0 ;
		      freeOutBinary ((char*)&xi, 4) ;
		      break ;
		    default:
		      freeOut (cp) ;
		      break ;
		    }
		  break ;
		case 'i':
		  xi = tabInt(t,n,i) ;
		  switch (style)
		    {
		    case 'j': 
		      freeOutf ("?int?%d?", xi) ;
		      break ;
		    case 'J': 
		      if (i > ii) freeOutf ("int") ;
		      freeOutf("?%d?", xi) ;
		      break ; 
		    case 'C':
		      freeOutBinary (bufinteger,1) ;
		      freeOutBinary ((char*)&xi, 4) ;
		      break ;
		    case 'A': 
		      freeOutf("Int:%d", xi) ;
		      break ;
		    default:
		      freeOutf ("%d", xi) ;
		      break ;
		    }
		  break ;
		case 'f':
		  xf = tabFloat (t,n,i) ;
		  switch (style)
		    {
		    case 'j': 
		      freeOutf ("?float?%g?", xf) ;
		      break ;
		    case 'J': 
		      if (i > ii) freeOutf ("float") ;
		      freeOutf("?%g?", xf) ;
		      break ;
		    case 'C':
		      freeOutBinary (buffloat,1) ;
		      freeOutBinary ((char*)&xf, 4) ;
		      break ;
		    case 'A': 
		      freeOutf("Float:%g", xf) ;
		      break ;
		    default:
		      {
			float zf = xf > 0 ? xf : -xf ;
			int izf = xf + .1 ;
			if (izf > 0 && zf - izf < ACE_FLT_RESOLUTION)
			  freeOutf ("%d", xf > 0 ? izf : -izf) ;
			else
			  freeOutf ("%g", xf) ;
		      }
		      break ;
		    }
		  break ;
		case 's': case 'D':
		  cp = tabString(t,n,i) ; 
		  switch (style)
		    {
		    case 'j': 
		      freeOutf ("?txt?%s?", freejavaprotect(cp)) ; 
		      break ;
		    case 'J': 
		      if (i > ii) freeOutf ("txt") ;
		      freeOutf ("?%s?", freejavaprotect(cp)) ; 
		      break ; 
		    case 'h': 
		      freeOut (cp) ;
		      break ;
		    case 'C':
		      freeOutBinary (buftext,1) ;
		      freeOutBinary (cp, strlen(cp) + 1) ;
		      break ;
		    case 'A': 
		      freeOutf ("Text:%s", freeprotect(cp)) ; 
		      break ;
		    default:
		      if (t->type[i] == 'D')
			freeOut (cp) ;
		      else
			freeOut (freeprotect(cp)) ;
		      break ;
		    }
		  break ;   
		case 't':
		  xt = tabDate(t,n,i) ;
		  switch (style)
		    {
		    case 'j': 
		      freeOutf ("?date?%s?", timeShowJava(xt, timeBuf, 25)) ; 
		      break ;
		    case 'J': 
		      freeOutf ("date") ;
		      freeOutf ("?%s?", timeShowJava(xt, timeBuf, 25)) ; 
		      break ; 
		    case 'h': 
		      freeOut (timeShow(xt, timeBuf, 25)) ;
		      break ; 
		    case 'C':
		      freeOutBinary (bufdate,1) ;
		      freeOutBinary ((char*)&xt, 4) ;
		      break ;
		    case 'A': 
		      freeOutf ("Date:%s", timeShow(xt, timeBuf, 25)) ; 
		      break ;
		    default:
		      freeOut (freeprotect(timeShow(xt, timeBuf, 25))) ;
		      break ;
		    }
		  break ;
		default:
		  messcrash ("table.c:tableSliceOut() - invalid column type %c", t->type[i]) ;
		  break ;
		}
	      if (i > ii) ii = i ; /* i's type has now been exported to java */
	    }
	}

      if (style == 'C')
	{ 
	  char bfret[3];
	  
	  bfret[0] = 'l';
	  xi = t->ncol-1 ;
	  while (xi > 0 && xi > 255)
	    {
	      bfret[1]= (char)255; freeOutBinary (bfret, 2) ;
	      xi -= 255;
	    }

	  bfret[1] = xi ; bfret[2] = '\n' ;
	  freeOutBinary (bfret, 3) ;
	  beginline[0] = '.';
	}
      else
	freeOut ("\n") ;
    }

  if (style == 'C')
    {
    /*
    * We are now at the end of an encore block.  At the beginning of the
    * next encore block, we will move the cursor right one before depositing
    * data, so here we move the cursor left to get it back to column -1.
    */
    freeOutBinary("l\001.",3);
    }

  messfree (colFlagSet) ;
  messfree (colFlagDict) ;
  arrayDestroy (oldClass) ;

  return n ;
} /* tableSliceOut */

/************************************************************/

const char *tablePrintable (TABLE *t, int i, int j)
{
  const char *ccp ;
  if (!t) return 0;
  if (t->magic != &TABLE_MAGIC)
    messcrash ("tableSliceOut() received non-magic TABLE* pointer") ;
  if (i < 0 || j < 0 ||
      i > tableMax(t) || j >= t->ncol ||
      tabEmpty(t,i,j)
      )
    return NULL ;

  switch (t->type[j])
    { 
    case 'k':
    case 'g':
      return tabKey(t,i,j) ? name(tabKey(t,i,j)) : 0 ;
    case 'a':
      return NULL ; /* je sais pas ce que c'est ce truc mieg 2007_03_08 */
    case 'b':
      return tabBool(t,i,j) ? "TRUE" : "FALSE" ;
    case 'i':
      sprintf (t->buf, "%d", tabInt(t,i,j)) ;
      return t->buf ;
    case 'f':
      sprintf (t->buf, "%g", tabFloat (t,i,j)) ;
      return t->buf ;
    case 's': case 'D':
      ccp = tabString(t,i,j) ; 
      return ccp && *ccp ? ccp : 0 ;
    case 't':
      timeShow(tabDate(t,i,j), t->buf, 25) ; 
      return t->buf ;
    default:
      messcrash ("table.c:tablePrintable() - invalid column type %c", t->type[i]) ;
      break ;
    }
  return NULL ;
} /* tablePrintable */

/************************************************************/

static BOOL tableDumpFormat (TABLE *table)
{
  int i, ii ;
  char *ct ;
  KEY key ;
  if (!table || ! tableMax(table))
    return FALSE ;
    
  ct = table->type - 1 ; ii = -1 ;
  freeOut("Format ") ;
  while(ii++, *++ct)
    switch(*ct)
      {
      case 'k':
      case 'g':
	for (i = 0 ; i < tabMax(table, ii) ; i++)
	  if (!tabEmpty(table,i,ii) && class(tabKey(table,i,ii)))
	    { key = tabKey(table,i,ii) ;
	      freeOutf(" %s ", class(key) ? className(key) : "Tag") ; 
	      break ;
	    }
	if (i == tabMax(table,ii))
	  freeOut(" Key ") ;
	continue ;
      case 'b':
	freeOut(" Bool ") ;
	continue ;
      case 'i':
	freeOut(" Int ") ;
	continue ;
      case 'f':
	freeOut(" Float ") ;
	continue ;
      case 's':
	freeOut(" Text ") ;
	continue ;
      case 'D':
	freeOut(" DNA ") ;
	continue ;
      case 't':
	freeOut(" Date ") ;
	continue ;
      case 'a':
	freeOut(" Flag ") ;
	continue ;
      }
  freeOut("\n") ;
  return TRUE ;
}

/* global for now, should go by chunk */
int tablePartOut (int *lastp, TABLE *t, char separator, char style)
{ 
  *lastp = 0 ;  /* global dump */
  if (!*lastp && (style == 'a' || style == 'h'))
    tableDumpFormat(t) ;
  return 
    tableOut (t, separator, style) ;
}

BOOL tableDump (FILE* fil, Stack s, KEY k) 
{
  int level = 0 ;
  TABLE *table = tableGet(k) ; 

  if (!table || ! tableMax(table))
    return FALSE ;

  if (fil)
    level = freeOutSetFile (fil) ;
  else if (s)
    level = freeOutSetStack (s) ;
  /* else freeOut */

  if (tableDumpFormat (table))
    tableOut (table, '\t', 'a') ;
  freeOutf ("\n") ;
  tableDestroy (table) ;
  
  if (level)
    freeOutClose(level) ;
  return TRUE ;
}

static TABLE *tableDoDoParse (int level, char *type, Array oldClass, KEY key)
{
  char *cp, *cq, *ctype ;
  TABLE *tt = 0 ;
  int xi ; float xf ; mytime_t xt ; KEY xk ;
  int ii, nn ;
  char *errText = "" ;

  if (!type || !*type)
    return 0 ;
  ctype = type - 1 ; ii = -1 ;
  while (ii++, *++ctype)
    switch (*ctype)
      {
      case 'k': case 'g': case 'b': case 'i': case 'f': case 's': case 't': case 'D':
	continue ;
      case 'a':
	errText = "tableDoParse cannot parse a type: flags" ;
	goto abort ;
      default:
	return 0 ;
      }
  tt = tableHandleCreate (100,type,0) ;
  nn = -1 ; /* line number */
  while (freecard(level))
    {
      nn++ ;
      ctype = type - 1 ; ii = -1 ;
      while (ii++, *++ctype)
	{
	  freenext() ;
	  switch (*ctype)
	    {
	    case 'k':
	    case 'g':
	      cp = freeword() ;
	      if (!cp) 
		{
		  if (ii) 
		    { errText = "missing object name" ; 
		      goto abort ;
		    }
		  else 
		    return tt ;
		}
	      if (!strcasecmp(cp,"NULL"))
		{ tableKey (tt, nn,ii) = 0 ;
		  break ;
		}
	      cq = cp ;
	      while (*cq && *cq != ':') cq++ ;
	      if (*cq == ':')
		{
		  *cq = 0 ;
		  if ((xi = pickWord2Class(cp)) &&
		      !pickList[xi].protected)
		    { *cq = ':' ;
		      array(oldClass, ii, int) = xi ;
		       /* but beware of "" around name */
		      freeback() ;
		      cp = freepos() ;
		      while(*cp != ':') *cp++ = ' ' ; /* garanteed to work */
		      *cp = ' ' ;
		      freenext() ;
		      cp = freeword() ;
		    }
		  else
		    *cq = ':' ;
		}
	      if (array(oldClass,ii,int) == -1)
		{ errText = "missing or unrecognized or protected class name" ; 
		  goto abort ;
		}
	      if (!cp || ! strlen(cq=lexcleanup(cp, 0)))
		{ errText = "missing object name" ; 
		  messfree (cq) ;
		  goto abort ;
		}
	      messfree (cq) ;
	      if (!array(oldClass,ii,int) &&
		  !lexword2key (cp, &xk,array(oldClass,ii,int)))
		{ errText = "Unrecognized Tag name" ; 
		  goto abort ;

		}
	      if (pickList[array(oldClass,ii,int)].known &&
		  !lexword2key (cp, &xk,array(oldClass,ii,int)))
		{ errText = "Unrecognized Tag name" ; 
		  goto abort ;
		}
	      lexaddkey(cp, &xk, array(oldClass,ii,int)) ;
	      tableKey(tt,nn, ii) = xk ;
	      break ;
	    case 'b':
	      cp = freeword() ;
	      if (!cp) 
		{
		  if (ii) 
		    { errText = "missing bool data" ; 
		      goto abort ;
		    }
		  else 
		    return tt ;
		}
	      if (!strcasecmp(cp,"true"))
		tableBool(tt, nn,ii) = TRUE ;
	      else if (!strcasecmp(cp,"false"))
		tableBool(tt, nn,ii) = FALSE ;
	      else if (!strcasecmp(cp,"NULL"))
		tableBool(tt, nn,ii) = UT_NON_BOOL ;
	      else
		{ errText = "Wrong Bool value, should be TRUE or FALSE" ;
                  goto abort ;
		}
	      break ;
	    case 'i':
	      if (freeint(&xi)) 
		tableInt(tt,nn, ii) = xi ;
	      else if ((cp = freeword()) &&  !strcasecmp(cp,"NULL"))
		tableInt(tt,nn, ii) = UT_NON_INT ;
	      else
		if (ii || cp) 
		  { errText = "wrong or missing  Int value" ;
		    goto abort ;
		  }
		else 
		  return tt ;
	      break ;
	    case 'f':
	      if (freefloat(&xf)) 
		tableFloat(tt,nn, ii) = xf ;
	      else if ((cp = freeword()) &&  !strcasecmp(cp,"NULL"))
		tableFloat(tt,nn, ii) = UT_NON_FLOAT ;
	      else
		if (ii || cp) 
		  { errText = "wrong or missing Float value" ;
		    goto abort ;
		  }
		else 
		  return tt ;
	      break ;
	    case 's': case 'D':
	      cp = freeword() ;
	      if (!cp) 
		{
		  if (ii) 
		    { errText = "Missing text" ;
		      goto abort ;
		    }
		  else 
		    return tt ;
		}
	      if (!strcasecmp(cp,"NULL"))
		tableDate(tt,nn,ii) = 0 ;
	      else
		tableSetString (tt, nn, ii, cp) ;
	      break ;
	    case 't':
	      cp = freeword() ;
	      if (!cp) 
		{
		  if (ii) 
		    { errText = "missing date" ;
		      goto abort ;
		    }
		  else 
		    return tt ;
		}
	      if (!strcasecmp(cp,"NULL"))
		tableDate(tt,nn,ii) = 0 ;
	      else if ((xt = timeParse(cp)))
		tableDate (tt,nn, ii) = xt ;
	      else
		{ errText = "wrong Date value" ;
		  goto abort ;
		}
	      break ;
	    }
	}
    }
  return tt ;
abort: 
  messerror (
	     "  TableResult parse error at line %7d column %d in %.25s : %s\n", 
	     freestreamline(level) -1, ii + 1, name(key), errText) ;
  tableDestroy (tt) ;
  return 0 ;
}

static char* tableParseFormat(int level, Array oldClass)
{
 char *cp ;
 static char type[1000] ;
 int nn, i = 1000 ;
 while(i--) type[i] = 0 ;

 if (!freecard(level))
   return 0 ;
 cp = freeword() ;
 if (!cp || strcasecmp("format", cp))
   return 0 ;
 i = -1 ;
 while (i < 1000 && (cp = freeword()))
   { i++ ; array(oldClass,i, int) = 0 ;
     if (!strcasecmp(cp,"bool"))
       type[i] = 'b' ;
     else if (!strcasecmp(cp,"int"))
       type[i] = 'i' ;
     else if (!strcasecmp(cp,"float"))
       type[i] = 'f' ;
     else if (!strcasecmp(cp,"text"))
       type[i] = 's' ;
     else if (!strcasecmp(cp,"DNA"))
       type[i] = 'D' ;
     else if (!strcasecmp(cp,"date"))
       type[i] = 't' ;
     else if (!strcasecmp(cp,"key"))
       { type[i] = 'k' ;
       array(oldClass,i, int) = -1 ;
       }
     else if (!strcasecmp(cp,"tag")) /* tag means _VSystem */
       { type[i] = 'k' ;
       array(oldClass,i, int) = 0 ;
       }
     else if ((nn = pickWord2Class (cp)))
       { type[i] = 'k' ;
         array(oldClass,i, int) = nn ;
       }
     else
       return 0 ;
   }
 return type ;
}
 
TABLE *tableDoParse (int level, KEY key) 
{
  TABLE *tt = 0 ;
  char *cp ;
  Array oldClass = arrayCreate(12, int) ;

  cp = tableParseFormat(level, oldClass) ;
  if (cp && *cp)
    tt = tableDoDoParse (level, cp, oldClass, key) ;
  arrayDestroy (oldClass) ;
  return tt ;
}

BOOL tableParse (int level, KEY key) 
{
  TABLE *tt = tableDoParse (level, key) ;
  if (!tt)
    return FALSE ;

  tableStore (key, tt) ; 
  tableDestroy (tt) ;
  return TRUE ;
}

/*************************************************/
/****************** Get/Store   ******************/
/*************************************************/

BOOL tableStore (KEY tKey,  TABLE *t) 
{
  int i, max ;
  Array aa = 0 ;
  TABLETYPE *uu, *vv; 

  if (t->magic != &TABLE_MAGIC)
    messcrash ("tableStore() received non-magic TABLE* pointer");

  if(t->ncol <= 0)
    messcrash ("tableStore() received TABLE with no columns");

  if (!tKey || pickType(tKey) != 'A')
    messcrash ("tableStore() called with a non-valid key") ;

  max = tableMax(t) ;

  aa = arrayCreate ((max+1) * t->ncol, TABLETYPE) ;
  uu = arrayp (aa, (max+1) * t->ncol - 1, TABLETYPE) ; /* make room */

  memcpy (arrp(aa, 0, TABLETYPE), t->name, t->ncol*sizeof(TABLETYPE)) ;

  if (max)
    for (i = 0 ; i < t->ncol ; i++)
      { uu = arrp (aa, t->ncol + i * max, TABLETYPE) ;
	vv = arrp(t->col[i], 0, TABLETYPE)  ; 
	memcpy (uu, vv, tabMax(t,i) * sizeof(TABLETYPE)) ;
      }

  arrayStackStore (tKey, aa, "a", t->s) ; /* crashes or succeeds */
  arrayDestroy (aa) ;

  return TRUE ;
} /* tableStore */

/*************************************************/

TABLE *tableHandleGet (KEY tKey, AC_HANDLE h)
{
  Stack s = 0 ;
  Array aa = 0 ;
  int i = 0, max, ncol, maxBitArray = 0 ;
  TABLETYPE *uu, *vv ;
  TABLE *t = 0 ;
  char *types ; 
  int  version ;
  char *vs ;

  if (!tKey || pickType(tKey) != 'A')
    messcrash ("tableHandleGet called with a non-valid key") ;

  if (!arrayStackHandleGet (tKey, &aa, "a", &s, h))
    goto abort ;

  if (!aa || !s || !arrayMax(aa) || !stackMark(s))
    messcrash ("missing data in tableHandleGet") ;

  types = stackText (s, 0) ;
  if (!types || !(ncol = strlen(types)))
    goto abort ;
  if (arrayMax(aa) % ncol)
    messcrash ("inconsistent array size in tableHandleGet") ;
  max = arrayMax(aa)/ncol - 1 ; /* 1 for the names */

  version = 0 ;
  if (!stackAtEnd(s) && (vs = stackNextText(s)) &&
      *vs == '\177')
    version = vs[1] ;		/* 256 possible versions */
  if (version > 0)
    { maxBitArray = 1 + (max-1)/(1 + 32) ;
      max -= maxBitArray ;
      if (sizeof(TABLETYPE) != 32)
	messcrash ("sizeof(TABLETYPE) = %d != 32 in tableGet: rethink!",
		   sizeof(TABLETYPE)) ;
    }

  t = tableHandleCreate (max, types, h) ;
  stackDestroy (t->s) ; t->s = s ;
  memcpy (t->name, arrp(aa, 0, TABLETYPE), ncol*sizeof(TABLETYPE)) ;

  if (max)
    for (i = 0 ; i < ncol ; i++)
      { uu = arrp (aa, ncol + i * max, TABLETYPE) ;
	array(t->col[i], max-1, TABLETYPE).i = 0 ; /* make room */
	array(t->parents[i], max-1, KEY) = 0; /* make room */
	array(t->grandParents[i], max-1, KEY) = 0; /* make room */
	vv = arrp(t->col[i], 0, TABLETYPE)  ; 
	memcpy (vv, uu, max*sizeof(TABLETYPE)) ;
      }
  if (version > 0)
    for (i = 0 ; i < ncol ; i++)
      { uu = arrp (aa, ncol*(1+max) + i*maxBitArray , TABLETYPE) ;
	bitUnSet(t->emptyBits[i], max-1) ; /* make room */
	vv = bigArrp(t->emptyBits[i], 0, TABLETYPE)  ; 
	memcpy (vv, uu, maxBitArray*sizeof(TABLETYPE)) ;
      }
  else
    for (i = 0 ; i < ncol ; i++)
      tableCheckEmpty (t,i) ;
  arrayDestroy (aa) ;
  return t ;  
  
 abort:
  stackDestroy (s) ;
  arrayDestroy (aa) ;
  tableDestroy (t) ;
  return 0 ;
}

/*************************************************/
/****************** end of file ******************/
/*************************************************/ 
