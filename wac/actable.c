#include "wh/acedb.h"
#include "wac/ac_.h"
#include "wh/bitset.h"

struct ac_table_internal
{
  int 	magic;
  int	rows;
  /*
   * number of rows in table == public view of arrayMax(lines)
   */

  int	cols;
  /*
   * number of columns in the widest row of the table
   * (tables extracted from objects do not necessarily
   * have the same number of columns in each row )
   */

  /*------------------------------------------------------------*/
  /* before this line, the structure layout must match exactly  */
  /* the layout of struct ac_table in ac.h                      */
  /*------------------------------------------------------------*/
  AC_HANDLE h ; /* everything allocated inside the table */
  Array lines ; /* Array of arrays of cells */

  int tableWidth ; /* a power of 2, defining the fixed geometry */
  Array lineLengths ;  /* the number of columns of each line */ 
  Array hardLines ; /* the table is always seen as a rectangle */

  int	def_cols;
  /*
   * this is the default number of columns to use when
   * creating a new row.
   */
  AC_DB db ;
  DICT *text_dict ; 		/* place to store the texts */
  DICT *tag_dict;	/* place to store the tags */

  AC_TABLE parent ;
  int parentRow, parentCol ;
  char buf25[25] ;
};

typedef struct cell
{
  ac_type	type;
  /*
   * type of data in this cell.  It is
   * ac_type_empty if there is no data in this cell.
   */
  union
  { int	i ;
    float f ;
    KEY   k ;
    mytime_t time ;
    AC_OBJ object;
  } u ;
} AC_CELL ;

static AC_CELL *ac_table_cell (struct ac_table_internal *t, int row, int col);


static void ac_table_finalize (void *tt)
{
  struct ac_table_internal *t;
  int x,y;
  AC_CELL *c;

  t = (struct ac_table_internal *)tt;
  if (!t->magic || t->parent)
    return ;
  for (x = 0; x< t->rows; x++)
    {
      for (y = 0; y< t->cols; y++)
	{
	  c = ac_table_cell(t, x, y);
	  if (c && (c->type == ac_type_obj))
	    ac_free(c->u.object);
	}
    }
  t->magic = 0;
  ac_free (t->h) ;
}


AC_TABLE ac_empty_table (int rows, int cols, AC_HANDLE handle )
{
  struct ac_table_internal *neuf;
  
  neuf = (struct ac_table_internal *) halloc (sizeof(struct ac_table_internal), handle);
  blockSetFinalise (neuf, ac_table_finalize) ;
  neuf->h = ac_new_handle () ;

  /*
   * rows and cols are a hint of the expected buffer size 
   * but we return an empty table, hence size (0,0)
   */
  neuf->rows = 0 ;
  neuf->cols = 0 ;
  neuf->def_cols = cols ;
  
  for (neuf->tableWidth = 1 ; neuf->tableWidth < cols ; neuf->tableWidth <<= 1) ;

  /* We allocate to the required size */
  neuf->lines = arrayHandleCreate (rows, Array, neuf->h) ;
  neuf->lineLengths = arrayHandleCreate (rows, Array, neuf->h) ;
  neuf->hardLines = arrayHandleCreate (rows * neuf->tableWidth, AC_CELL, neuf->h) ;

  /*
   * each cell is initialised at zero by arrayCreate which
   * is the correct "nothing here" value for all fields of
   * struct row.
   */
  /*
   * all that is left is the magic number, and to return it
   */
  neuf->magic = MAGIC_BASE + MAGIC_AC_TABLE;

  return (AC_TABLE) neuf;
}  /* ac_empty_table */


/*
 * PRIVATE
 */
AC_TABLE ac_db_empty_table (AC_DB db, int rows, int cols, AC_HANDLE handle)
{
  AC_TABLE neuf = ac_empty_table (rows, cols, handle) ;
  ((struct ac_table_internal *)neuf)->db = db ;
  return neuf ;
}

/*
 * PRIVATE
 */
AC_DB ac_table_db (AC_TABLE atable) 
{ return  ((struct ac_table_internal *)atable)->db ; }

/*
 * PRIVATE
 */
static AC_CELL *ac_table_new_cell (struct ac_table_internal *t, int row, int col)
{
  Array line = array (t->lines, row, Array) ; /* automatic reallocation */
  AC_CELL *c ;

  if (!line)
    line = array (t->lines, row, Array) = arrayHandleCreate (t->def_cols, AC_CELL, t->h) ;

  c = arrayp (line, col, AC_CELL) ;

  if (arrayMax (line) > t->cols)
    t->cols = arrayMax (line) ;
  if (arrayMax(t->lines) > t->rows)
    t->rows = arrayMax(t->lines) ;
  return c ; /* cannot fail */  
}

/*
 * PRIVATE
 */
static AC_CELL *ac_table_cell (struct ac_table_internal *t, int row, int col)
{
  Array *linep ;
  
  if (! t)
    return 0;
  if (row < 0)
    return 0;
  if (row >= t->rows)
    return 0;
  if (col >= t->cols)
    return 0;
  while (t->parent)
    { 
      row += t->parentRow ; col += t->parentCol ;
      t =  (struct ac_table_internal *) (t->parent) ;
    }
  linep = arrp (t->lines, row, Array);
  if (!linep)
    return 0;
  if (! * linep)
    return 0;
  if (col < 0)
    return 0;
  if (col >=  arrayMax(*linep))
    return 0;

  return arrp(*linep, col, AC_CELL) ;
}

/****************************************************************/
/* PRIVATE: called by ac_parse_buffer
 *  we make no check 
 */
void ac_tableFillUp (AC_TABLE table, int row)
{
  struct ac_table_internal *t =  (struct ac_table_internal *) table ;
  int i, j ;
  AC_CELL *c, *c2 ;
  Array ll, l0 ;
  
  for (i = row + 1 ; i < t->rows ; i++)
    {
      ll = arr (t->lines, i, Array) ;
      l0 = arr (t->lines, i - 1, Array) ;
      for (j = arrayMax (ll) - 1, c = arrp (ll, j, AC_CELL) ; j >= 0 ; c--, j--)
	if (c->type != ac_type_empty) break ; /* find a filled cell */
      for (; j >= 0 ; c--, j--)
	if (c->type == ac_type_empty) /* find an empty cell to its left */
	  for (c2 = arrp  (l0, j, AC_CELL) ; j >= 0 ; c--, c2--, j--)
	    *c = *c2 ; /* copy the line above */
    }
}

void ac_table_copydown(AC_TABLE table, int row, int col)
{
  struct ac_table_internal *t = (struct ac_table_internal *)table;
  int prev, x;
  AC_CELL *above, *below;

  prev = row - 1;
  for (x = 0; x< col; x++)
    {
      above = ac_table_cell(t, prev, x);
      below = ac_table_new_cell(t, row, x);
      *below = *above;
      if (above->type == ac_type_obj)
	below->u.object = ac_copy_obj (above->u.object, t->h);
    }
}

/****************************************************************/
/* PRIVATE: called by ac_tag_table
 *
 * Implement ac_tag_table().  row, col points at the tag that we are looking
 * for.  That is pparent[row,col] = tag.  We want to find the table in col+1
 * from row to row+n for as many rows as contain tag in [row,col]
 */
AC_TABLE ac_subtable (AC_TABLE pparent, int row, int col, AC_HANDLE handle)
{
  AC_TABLE neuf = 0 ;
  struct ac_table_internal *parent =  (struct ac_table_internal *) pparent ;
  int i, j, rows, cols ;
  KEY k ;
  AC_CELL *c, *c2 ;
  ac_type type ;

  /* 
   * Subtables have no data of their own.  It is all in the parent.  If we
   * are looking in a subtable, we need to walk up through the parent pointer
   * to find the table with the real data in it.
   */
  rows = cols = 0 ;
  while (parent->parent)
    { 
      row += parent->parentRow ; 
      col += parent->parentCol ;
      parent =  (struct ac_table_internal *) (parent->parent) ;
    }

  /*
   * type is the type of the cell we are pointing at.  We assume that it
   * is a tag.  k is the identity of the tag (currently an int pointing
   * into a dict, but in principle any unique identifier is good).
   */
  c = ac_table_cell (parent, row, col) ;
  type = c->type ; 
  k = c->u.i ;

  /*
   * walk down the rows, starting at the current row.  Find the last
   * row that still has tag k in column col.  While doing that, find
   * the widest row that we have.
   */
  for (i = 0 ; ; i++)
    {
      /*
       * if the cell in the next row is still the same thing, 
       */
      c2 = ac_table_cell (parent, row + i, col) ;
      if (c2 && c2->type == type && c2->u.i == k)
	{
	  rows++ ;
	  /*
	   * find the max column in that row
	   */
	  for (j = cols ; ; j++) 
	    {
	      c2 = ac_table_cell (parent, row + i, col + j) ;
	      if (!c2 || c2->type == ac_type_empty)
		break ;
	      if (j >= cols) cols = j + 1 ;
	    }
	}
      else
	{
	  /*
	   * stop when we find a different tag in the first column.
	   */
	  break ;
  	}
    }

  /*
   * I don't understand why this is conditional.  It looks like
   * rows and cols are both guaranteed to be at least 1 because of 
   * the code above.
   */
  if (rows && cols)
    {
      neuf = ac_copy_table (pparent, row, col+1, rows, cols, handle);
    }
  return (AC_TABLE) neuf;
}


AC_TABLE ac_copy_table (AC_TABLE xfrom, int row, int col, int copy_rows, int copy_cols , AC_HANDLE handle)
{
struct ac_table_internal *from =  (struct ac_table_internal *) xfrom ;
AC_TABLE neuf;
int i, j;

/*
* copy a square area out of a table into a new table
*/
neuf = ac_empty_table(copy_rows, copy_cols, handle);

for (i = 0; i < copy_rows; i++)
for (j = 0; j < copy_cols; j++)
  {
    /*
      neuf[i,j] = old[i+row, j+col]
    */
    AC_CELL *c;
    c = ac_table_cell(from, i+row, j+col);
    if (! c)	
      continue;
    switch (c->type)
      {
      case ac_type_empty:
	break;
      case ac_type_bool:
      case ac_type_key:
      case ac_type_int:
      case ac_type_float:
      case ac_type_date:
	ac_table_insert_type(neuf , i, j, &c->u, c->type);
	break;
      case ac_type_text:
	{
	  const char *s = dictName (from->text_dict, c->u.i) ;
	  ac_table_insert_type (neuf , i, j, &s, c->type);
	}
	break;
      case ac_type_tag:
	{
	  const char *s = dictName (from->tag_dict, c->u.i) ;
	  ac_table_insert_type(neuf , i, j, &s, c->type);
	}
	break;
      case ac_type_obj:
	{
	  /*
	   * do not need to ac_copy_obj () the c->u.object here because
	   * ac_table_insert_type does that internally.
	   */
	  ac_table_insert_type(neuf , i, j, &c->u.object, c->type);
	}
	break;
      }
  }
return neuf;
}


ac_type ac_table_type ( AC_TABLE tbl, int row, int col )
{
  struct ac_table_internal *t = (struct ac_table_internal *) tbl;
  AC_CELL *c = ac_table_cell (t, row, col) ;

  return c ? c->type : ac_type_empty;
} /* ac_table_type */


/*
 * to get values from a (row,col) of the table, we have a common
 * theme:  if the type of the cell is right, return the value,
 * else return a default.  Note that if the type of the cell is
 * not ac_type_empty, then you know that the cell exists, so you
 * do not need any array checking when fetching the value.
 */

BOOL ac_table_bool ( AC_TABLE tbl, int row, int col, BOOL deflt )
{
  struct ac_table_internal *t = (struct ac_table_internal *) tbl;
  AC_CELL *c = ac_table_cell (t, row, col) ;

  if (c && c->type == ac_type_bool)
    return c->u.k ? TRUE : FALSE ;
  return deflt;
}

KEY ac_table_key (AC_TABLE tbl, int row, int col, KEY deflt )/* removed: breaks the interface */
{
  struct ac_table_internal *t = (struct ac_table_internal *) tbl;
  AC_CELL *c = ac_table_cell (t, row, col) ;

  if (c && c->type == ac_type_key)
    return c->u.k ;
  else if (c && c->type == ac_type_tag)
    return c->u.k ;
  else if (c && c->type == ac_type_obj)
    return ac_obj_key (c->u.object) ;
  return deflt;
}

int ac_table_int (AC_TABLE tbl, int row, int col, int deflt )
{
  struct ac_table_internal *t = (struct ac_table_internal *) tbl;
  AC_CELL *c = ac_table_cell (t, row, col) ;

  if (c && c->type == ac_type_int)
    return c->u.i ;
  return deflt;
}

const char *ac_table_tag (AC_TABLE tbl, int row, int col, const char *deflt )
{
  struct ac_table_internal *t = (struct ac_table_internal *) tbl;
  AC_CELL *c = ac_table_cell (t, row, col) ;

  if (c && c->type == ac_type_tag)
    return dictName (t->tag_dict, c->u.i) ;
  else if (c && c->type == ac_type_key)
    return ac_key_name (c->u.k) ;
  return deflt;
}


float ac_table_float (AC_TABLE tbl, int row, int col, float deflt )
{
  struct ac_table_internal *t = (struct ac_table_internal *) tbl;
  AC_CELL *c = ac_table_cell (t, row, col) ;

  if (c && c->type == ac_type_float)
    return c->u.f ;
  return deflt;
}

mytime_t ac_table_date (AC_TABLE tbl, int row, int col, mytime_t deflt)
{
  struct ac_table_internal *t = (struct ac_table_internal *) tbl;
  AC_CELL *c = ac_table_cell (t, row, col) ;

  if (c && c->type == ac_type_date)
    return c->u.time ;
  return deflt;
}

const char * ac_table_text (AC_TABLE tbl, int row, int col, const char * deflt )
{
  struct ac_table_internal *t = (struct ac_table_internal *) tbl;
  AC_CELL *c = ac_table_cell (t, row, col) ;

  if (c && c->type == ac_type_text )
    return dictName (t->text_dict, c->u.i) ;
  else if (c && c->type == ac_type_obj && !strcmp (ac_class(c->u.object),"Text"))
    return ac_name(c->u.object) ;
  return deflt;
}

AC_OBJ ac_table_obj (AC_TABLE tbl, int row, int col, AC_HANDLE handle)
{
  struct ac_table_internal *t = (struct ac_table_internal *) tbl;
  AC_CELL *c = ac_table_cell (t, row, col) ;

  if (c && c->type == ac_type_obj)
    return ac_copy_obj (c->u.object, handle);
  else if (c && c->type == ac_type_key)
    return ac_key2obj (t->db, c->u.k, TRUE, handle) ;
  return NULL;
}

/*
 * inserting data into the table
 */

int ac_table_insert_type (AC_TABLE tbl, int row, int col, const void *dat, ac_type type)
{
  struct ac_table_internal *t = (struct ac_table_internal *) tbl;
  AC_CELL *c = ac_table_new_cell (t, row, col) ;

  c->type = type;
  if (dat)
    switch (type)
      {
      case ac_type_empty:
	c->type = ac_type_empty ;
	break;
      case ac_type_bool:
	c->u.k = *(BOOL*)dat ? TRUE : FALSE ;
	break;
      case ac_type_key:
	c->u.k = *(KEY *)dat;
	break;
      case ac_type_int:
	c->u.i = *(int *)dat;
	break;
      case ac_type_float:
	c->u.f = *(float *)dat;
	break;
      case ac_type_text:
	{
	  char *cp = *(char**) dat ;
	  if (cp && *cp)
	    {
	      if (! t->text_dict) 
		t->text_dict = dictCaseSensitiveHandleCreate(200, t->h);
	      dictAdd (t->text_dict, cp, &(c->u.i)) ;
	    }
	  else
	    c->type = ac_type_empty ;
	}
	break;
      case ac_type_date:
	c->u.time = *(mytime_t*) dat ;
	break;
      case ac_type_tag:
	{
	  char *cp = *(char**) dat ;
	  if (! t->tag_dict) 
  	    t->tag_dict = dictHandleCreate(200, t->h);
	  if (cp)
	    dictAdd (t->tag_dict, cp, &(c->u.i)) ;
	}
	break;
      case ac_type_obj:
	c->u.object = ac_copy_obj (* (AC_OBJ *)dat,t->h);
	break;
      default:
	messcrash ("bad type provided to ac_table_insert_type") ;
      }
  return TRUE ; /* we garantee success */
} /* ac_table_insert_type */

/*
 * an easier to use insert for text.
 */

int ac_table_insert_text (AC_TABLE tbl, int row, int col, const char *value)
{
  return ac_table_insert_type (tbl, row, col, (void *)&value, ac_type_text );
}

/*
 * a function to help in displaying a table.  I originally considered
 * a complete print_table() with callbacks, but it would be so complicated
 * to use that it is not worth doing.  This converts a single cell
 * to a string.
 */

const char *ac_table_printable (AC_TABLE tbl, int row, int col, const char *deflt)
{
  struct ac_table_internal *t = (struct ac_table_internal *) tbl;
  AC_CELL *c = ac_table_cell (t, row, col) ;

  if (c)
    switch (c->type)
      {
      case ac_type_empty:
	return deflt;
      case ac_type_bool:
	return c->u.k ? "TRUE" : "FALSE" ; 
      case ac_type_int:
	sprintf(t->buf25,"%d", c->u.i);
	return t->buf25;
      case ac_type_float:
	sprintf(t->buf25,"%g", c->u.f);
	return t->buf25;
      case ac_type_text:
	return c->u.i > 0 && c->u.i <= dictMax(t->text_dict) ? dictName (t->text_dict, c->u.i) : "" ;
      case ac_type_date:
	return timeShow (c->u.time, t->buf25, 25) ;
      case ac_type_tag:	
	return dictName (t->tag_dict, c->u.i) ;
      case ac_type_key:	
	return ac_key_name (c->u.k) ;
      case ac_type_obj:
	return ac_name(c->u.object);
      default:
	messcrash ("bad type in ac_table_printable") ;
      }
  return deflt ;
}

/***************************************************************************/

const char *ac_table_class (AC_TABLE tbl, int row, int col, const char *deflt)
{
  struct ac_table_internal *t = (struct ac_table_internal *) tbl;
  AC_CELL *c = ac_table_cell (t, row, col) ;

  if (c)
    switch (c->type)
      {
      case ac_type_obj:
	return ac_class (c->u.object);
      default:
	return deflt;
      }
  return deflt ;
}

/***************************************************************************/
/* 
 * Export a table
 * return the number of lines actually printed
 *
 * -if vtxtMarkup (blkp) was called previously, the function exports
 *  an html table, otherwise, a tab delimited table.
 * -titles: optional, in non-html the title line stars with a #
 * -cols[]: list of columns to be displayed: i.e. {3,1,4, 0}
 *    if cols == 0, all columns are exported in their natural order
 *    else the list must end on 0.
 *    In addition, empty columns are systematically hidden.
 * -colorControlColumn: ignored, defaults to the first displayed column
 *   The background color alternates between withe and pale blue when the
 *   content of this column changes.
 *   Note that colorControlColumn may be hidden (not listed in cols[])
 * -maxLine: maximal number of lines to be exported
 *   If 0, the whole table is exported
 * -more:  i.e. "<a href='wholetable.cgi'>...>/a>"
 *   if non zero and the table exceeds maxLine
 *   *more is printed in each colun of the last line
 *
 * the rows must be sorted as wished before this function is called
 *
 * In html5, the table is downloadable and the rows may be reordered
 * on the client side thanks to a server provided javascript.
 */
int ac_table_display (vTXT blkp 
		      , AC_TABLE t, const char **titles
		      , int *cols, int colorControlColumn
		      , int beginLine, int maxLine, const char *more
		      , char beauty
		      )
{
  struct ac_table_internal *ti = (struct ac_table_internal *) t ;
  int ii, i1, jj, jjj, nLine = 0 ;
  const char *ccp ;
  char previous[256] ;
  int mycols [t->cols + 1] ;
  BOOL markUp = blkp->markUp ;
  char bufvoid[1] = {'v'} ;
  char bufinteger[1] = {'i'} ;
  char buffloat[1] = {'f'} ;
  char bufdate[1] = {'d'} ;
  char buftext[1] = {'t'} ;
  char bufk[2] = {'k', 0} ;
  char bufK[2] = {'K', 0} ;
  char buftag[1] = {'g'} ;
  char beginline[1] ={'>'};
  static int tableId = 0 ;
  
  if (beauty == 'J') beauty = 'j' ;
  if (!cols)
    {
      for (jjj = 0 ; jjj < t->cols ; jjj++)
	mycols[jjj] = jjj ;
      mycols[jjj] = -1 ;
    }
  else   /* verify */
    {
      for (jjj = 0 ; cols &&  cols[jjj] > 0 ; jjj++)
	{
	  jj = cols[jjj] -1  ; /* actual column to be displayed */
	  if (jj >= t->cols)
	    {
	      messerror ("ac_table_display was asked to diplay in column %d the value t[%d], but %d > %d== t->cols"
			 , jjj, jj, jj, t->cols) ;
	      continue ;
	    }
	  mycols[jjj] = cols[jjj] - 1 ;
	}
      mycols[jjj] = -1 ;
    }
  
  if (markUp)
    {
      tableId++ ;
      vtxtPrintf (blkp, "\n<br>\n<button onclick=\"exportToCSV(\'table%d\', \'table%d.csv\')\">Download this table</button>\t\t<i><font color=\"#007f7f\">Click on the headers to reorder the lines</font></i>", tableId, tableId) ;
      vtxtPrintf (blkp, "\n<table id=\"table%d\" width=\"98%%\" border=3>\n", tableId) ;
    }
  
  if (maxLine <= 0)
    maxLine = t->rows ;
  /* eliminates empty columns */
  if (! beauty)
    for (jjj = 0 ; mycols[jjj] >= 0 ; jjj++)
      {
	BOOL empty = TRUE ;
	jj = mycols[jjj]  ; /* actual column to be displayed */
	for (ii = 0 ; empty && ii < t->rows && ii < maxLine ; ii++)
	  {
	    if (ac_table_type (t, ii, jj) !=  ac_type_empty)
	      empty = FALSE ;
	  }
	if (empty)
	  mycols[jjj] = 999999 ;
      }
  /* export the titles */
  if (titles)
    {
      int jCol = 0 ;
      if (markUp)
	vtxtPrint (blkp, "<thead>\n  <tr VALIGN=TOP bgcolor=\"#afafff\">\n") ;
      else
	vtxtPrint (blkp, "#") ;
      for (jjj = 0 ; mycols[jjj] >= 0 ; jjj++)
	{
	  jj = mycols[jjj]  ; /* actual column to be displayed */
	  if (jj ==  999999)
	    continue ;
	  ccp = titles[jj] ;
	  if (!ccp) 
	    ccp = markUp ? "&nbsp;" : "" ;
	  if (markUp)
	    vtxtPrintf (blkp, "    <th onclick=\"sortTable(\'table%d\',%d)\" data-order=\"asc\">%s</th>\n", tableId, jCol++, ccp) ;
	  else
	    vtxtPrintf (blkp, "%s%s", jjj ? "\t" : "", ccp) ;
	}
      vtxtPrint (blkp, "\n") ; nLine++ ;
    }
  if (markUp)
    vtxtPrint (blkp, "  </tr>\n</thead>\n\n<tbody>\n") ;

  /* export the table */
  if (1 || !colorControlColumn)
    colorControlColumn = mycols[0] ;
  else
    colorControlColumn-- ;
  if(beauty == 'C')
    {
      beginline[0] = '>' ;
    }
  if (beginLine < 0) beginLine = 0 ;
  previous[0] = 0 ;
  for (ii = beginLine, i1 = 0 ; ii < t->rows && ii < maxLine + beginLine ; ii++)
    {
      /* alternate colors when the item in first column changes */
      jj = colorControlColumn ; 
      ccp = ac_table_printable (t, ii, jj, "") ;
      if (ii && strncmp (ccp, previous, 255))
	i1++ ;
      strncpy (previous, ccp, 255) ;
     
      if (markUp)
	{
	  if (i1 & 0x1)
	    vtxtPrint (blkp,"  <tr VALIGN=TOP data-parity=\"even\">\n") ; 
	  else
	    vtxtPrint (blkp, "  <tr VALIGN=TOP data-parity=\"odd\">\n") ;
	}
      for (jjj = 0 ; mycols[jjj] >= 0 ; jjj++)
	{
	  if(beauty == 'C')
	    {
	      freeOutBinary (beginline, 1);
	      beginline[0]='>';
	    }
	  jj = mycols[jjj]  ; /* actual column to be displayed */
	  if (jj ==  999999)
	    continue ;
	  ccp = ac_table_printable (t, ii, jj, 0) ;
	  if (ccp && ! strcmp (ccp, "\177\177(NULL KEY)"))
	    ccp = 0 ;
	  if (beauty == 'a')
	    { 
	      AC_CELL *c = ac_table_cell (ti, ii, jj) ;
	      char *ccp1 = 0 ;
	      if (ccp)
		{
		  switch (c->type)
		    {
		    case ac_type_empty:
		    case ac_type_bool:
		    case ac_type_int:
		    case ac_type_float:
		      break ;
		    case ac_type_date:
		    default:
		      ccp = ccp1 = ac_protect (ccp, 0) ;
		      break ;
		    }
		}
	      else
		ccp = "NULL" ;
	      vtxtPrintf (blkp, "%s%s", jjj ? "\t" : "", ccp) ;
	      if (ccp1)
		ac_free (ccp1) ;
	    }
	  else if (beauty == 'C')
	    { 
	      AC_CELL *c = ac_table_cell (ti, ii, jj) ;
	      char *ccp1 = 0 ;
	      KEY xk ;
	      float xf ;
	      int xi ;
	      if (ccp)
		{
		  switch (c->type)
		    {
		    case ac_type_tag:
		      xk = ac_table_key (t, ii, jj, 0) ;
		      freeOutBinary (buftag,2) ; 
		      freeOutBinary ((char*)&xk, 4) ;
		      break ;
		    case ac_type_key:
		      xk = ac_table_key (t, ii, jj, 0) ;
		      if (iskey (xk) == 2)
			{ bufK[1] = class (xk) ; freeOutBinary (bufK,2) ; }
		      else
			{ bufk[1] = class (xk) ; freeOutBinary (bufk,2) ; }
		      freeOutBinary (ccp, strlen(ccp) + 1) ;
		      break ;
		    case ac_type_bool:
		      freeOutBinary (bufinteger,1) ;
		      xi = ac_table_int(t, ii, jj, 0) ;
		      if (xi)
		        xi = 1 ;
		      else
			xi = 0 ;
		      freeOutBinary ((char*)&xi, 4) ;
		      break ;
		    case ac_type_int:
		      xk = ac_table_int (t, ii, jj, 0) ;
		      freeOutBinary (bufinteger,1) ;
		      freeOutBinary ((char*)&xi, 4) ;
		      break ;
		    case ac_type_float:
		      xk = ac_table_float (t, ii, jj, 0) ;
		      freeOutBinary (buffloat,1) ;
		      freeOutBinary ((char*)&xf, 4) ;
		      break ;
		    case ac_type_date:
		      xi = ac_table_date (t, ii, jj, 0) ;
		      freeOutBinary (bufdate,1) ;
		      freeOutBinary ((char*)&xi, 4) ;
		      break ;
		    case ac_type_empty:
		      freeOutBinary (bufvoid, 1) ;
		      break ;
		    case ac_type_text:
		    default:
		      freeOutBinary (buftext, 1) ;
		      freeOutBinary (ccp, strlen(ccp) + 1) ;
		      break ;
		    }
		}
	      else
		freeOutBinary (bufvoid, 1) ;
	      ac_free (ccp1) ;
	    }
	  else if (beauty == 'j')
	    { 
	      AC_CELL *c = ac_table_cell (ti, ii, jj) ;
	      char *ccp1 = 0 ;
	      if (ccp)
		{
		  switch (c->type)
		    {
		    case ac_type_empty:
		      vtxtPrintf (blkp, "%s???", jjj ? "\t" : "") ;
		      break ;
		    case ac_type_bool:
		      vtxtPrintf (blkp, "%s?bool?%s?", jjj ? "\t" : ""
				  , ccp) ;
		      break ;
		    case ac_type_int:
		      vtxtPrintf (blkp, "%s?int?%s?", jjj ? "\t" : ""
				  , ccp) ;
		      break ;
		    case ac_type_float:
		      vtxtPrintf (blkp, "%s?float?%s?", jjj ? "\t" : ""
				  , ccp) ;
		      break ;
		    case ac_type_date:
		      vtxtPrintf (blkp, "%s?date?%s?", jjj ? "\t" : ""
				  , ccp) ;
		      break ;
		    case ac_type_text:
		      vtxtPrintf (blkp, "%s?txt?%s?", jjj ? "\t" : ""
				  , ccp) ;
		      break ;
		    default:
		      vtxtPrintf (blkp, "%s?%s?%s?", jjj ? "\t" : ""
				  , ac_key_class (ac_table_key(t, ii, jj,0))
				  , ccp
				  ) ;
		      continue ;
		      break ;
		    }
		}
	      else
		vtxtPrintf (blkp, "%s???", jjj ? "\t" : "") ;
	      if (ccp1)
		ac_free (ccp1) ;
	    }
	  else  if (markUp)
	    {
	      if (!ccp)
		ccp = "&nbsp;" ;
	      vtxtPrintf (blkp,"    <td>%s</td>\n", ccp) ;
	    }
	  else if (beauty == 'h')
	    {
	      vtxtPrintf (blkp, "%s%s", jjj ? "\t" : "", ccp ? ccp : "NULL") ;
	    }
	  else
	    {
	      vtxtPrintf (blkp, "%s%s", jjj ? "\t" : "", ccp ? ccp : "") ;
	    }
	} 
      if (beauty == 'C')
	{ 
	  char bfret[3];
	  int xi = t->cols - 1 ;
	  
	  bfret[0] = 'l';
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
	{    
	  if (markUp)
	    vtxtPrint (blkp, "</tr>") ;
	  vtxtPrint (blkp, "\n") ; nLine++ ;
	}
    }
  if (t->rows > maxLine && more && *more)
    {
      ccp = more ; 
      if (markUp)
	vtxtPrint (blkp, "  <tr VALIGN=TOP bgcolor=#d0d0ff>\n") ; 
      for (jjj = 0 ; mycols[jjj] >= 0 ; jjj++)
	{
	  jj = mycols[jjj]  ; /* actual column to be displayed */
	  if (jj ==  999999)
	    continue ;
	  if (markUp)
	    vtxtPrintf (blkp,"  <td>%s</td>", ccp) ;
	  else
	    vtxtPrintf (blkp, "%s%s", jjj ? "\t" : "", ccp) ;
	}
      if (1)
	{
	  if (markUp)
	    vtxtPrint (blkp, "</tr>") ;
	  vtxtPrint (blkp, "\n") ; nLine++ ;
	} 
    }
  if (beauty == 'C')
    {
      /*
       * We are now at the end of an encore block.  At the beginning of the
       * next encore block, we will move the cursor right one before depositing
       * data, so here we move the cursor left to get it back to column -1.
       */
      freeOutBinary("l\001.",3);
    }
  else
    {
      if (markUp)   
	vtxtPrint (blkp, "</tbody>\n</table>\n") ;
      vtxtBreak (blkp) ;
    }
  return nLine ;
} /* ac_table_display */

/*
* This function is really just here for use in the debugger.
*/
int ac_print_table(FILE *fp, AC_TABLE t)
{
  int x,y;

  if (!fp)
    fp = stderr;
  for (x = 0; x<t->rows; x++)
    {
      for (y = 0; y< t->cols; y++)
	fprintf(fp,"%d %d -%s-\n",x,y,ac_table_printable(t,x,y,"")) ;
      fprintf(fp,"\n");
    }
  return 0;
}

/*********************************************************/

static struct ac_table_internal *sortTable = 0 ;
static Array sortSpecList = 0;

static int acTableOrder (const void* va, const void* vb)
/* comparison function used for qsort() in tableSort,
   that also rearranges an index array of row numbers,
   requires static sortSpec and sortTable to be assigned */
{
  register int rowA = *(int*)va ;
  register int rowB = *(int*)vb ;
  register int i, col, sortSpec;
  ac_type typeA, typeB ;
  AC_CELL *cellA, *cellB ;
  struct ac_table_internal *t = sortTable ;
  int x ; float z ;

  for (i = 0 ; i < arrayMax (sortSpecList) ; ++i)
    {
      sortSpec = arr(sortSpecList,i,int);
      
      /* get the true sort column in C-numbers */
      if (sortSpec < 0)	/* descending */
	col = - sortSpec - 1;
      else /* sortSpec > 0   -> ascending  */
	col = sortSpec - 1;
      
      /* we make NULL the largest value ever */
      cellA = ac_table_cell (t, rowA, col) ;
      cellB = ac_table_cell (t, rowB, col) ;
      typeA = cellA->type ;
      typeB = cellB->type ;
      
      if (typeA == ac_type_empty)
	{
	  if (typeB == ac_type_empty)
	    continue ;
	  return sortSpec > 0 ? 1 : -1 ;
	}
      
      if (typeB == ac_type_empty)
	return sortSpec>0 ? -1 : 1 ;

      if (cellA->u.k == cellB->u.k)
	continue ;
      
      switch (typeA)
	{
	case ac_type_empty:
	  break ;
	case ac_type_tag:
	case ac_type_key:
	  {
	    KEY keyA = cellA->u.k ;
	    KEY keyB = cellB->u.k ;
	    if (acCaseSensitive (keyA))
	      x = lexstrCasecmp (ac_key_name(keyA), ac_key_name(keyB));
	    else
	      x = lexstrcmp (ac_key_name(keyA), ac_key_name(keyB));
	    if (x)
	      return sortSpec > 0 ? x : -x ;
	  }
	  break ;
	case ac_type_obj:
	  {
	    KEY keyA = ac_obj_key (cellA->u.object) ;
	    KEY keyB = ac_obj_key (cellB->u.object) ;
	    x = lexstrcmp (ac_key_name(keyA), ac_key_name(keyB));
	    if (x)
	      return sortSpec > 0 ? x : -x ;
	  }
	  break ;
	case ac_type_bool:
	  x = cellA->u.k - cellB->u.k ;
	  if (x)
	    return x * sortSpec > 0 ? 1 : -1 ;
	  break ;
	case ac_type_int:
	  x = cellA->u.i - cellB->u.i ;
	  if (x)
	    return x * sortSpec > 0 ? 1 : -1 ;
	  break ;
	case ac_type_float:
	  z = cellA->u.f - cellB->u.f ;
	  if (z != 0)
	    return z * sortSpec > 0 ? 1 :-1 ;
	  break ;
	case ac_type_text:
	  {
	    const char *sA = dictName (t->text_dict, cellA->u.i) ;
	    const char *sB = dictName (t->text_dict, cellB->u.i) ;
	    x = lexstrcmp (sA, sB) ;
	    if (x)
	      return x * sortSpec > 0 ? 1 : -1 ;
	  }
	  break ;
	case ac_type_date:
	  {
	    mytime_t tA = cellA->u.time ;
	    mytime_t tB = cellB->u.time ;
	    x = tA - tB ;
	    if (x)
	      return x * sortSpec > 0 ? 1 : -1 ;
	  }
	  break ;
	} /* end switch */
    }

  /* default keep original order */
  return 0 ;
} /* acTableOrder */

/************************************************************/
/* adapted from tableSort, used in Aql
*/
static BOOL ac_table_do_sort (AC_TABLE table, const char *spec)
{
  AC_HANDLE h ;
  int i, i1, j, rows, specLen, sortCol ; 
  Array indexArray = 0 ;
  const char *ccp ;
  BOOL isPlus = TRUE ;
  struct ac_table_internal *t = (struct ac_table_internal *)table ;
  vTXT vSpec = 0 ;

  if (!t) return FALSE; 
  if (t->magic != MAGIC_BASE + MAGIC_AC_TABLE)
    messcrash ("tableSort() received non-magic TABLE* pointer") ;
  rows = t->rows ;
  if (rows <= 1) return TRUE; /* 1 or no elements are always sorted */

  h = ac_new_handle () ;

  if (!spec)
    {
      vSpec = vtxtHandleCreate (h) ;
      for (i = 1 ; i <= t->cols ; i++)
	vtxtPrintf (vSpec, "+%d", i) ;
      spec = vtxtPtr (vSpec) ;
    }

  specLen = strlen(spec);
  sortSpecList = arrayReCreate(sortSpecList, (specLen), int);

  for (i = j = 0, ccp = spec ; *ccp ; i++)
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
      if (sortCol < 1 || sortCol > t->cols) 
	/* column numbers have to range from 1 to colNum */
	continue ;
      if (isPlus)
	array(sortSpecList, j, int) = sortCol;
      else 
	array(sortSpecList, j, int) = -sortCol;
      j++ ;
    }

  /***************************************/

  /* create the indices that keep track of which row went where
     during comparative sorting */
  indexArray = arrayHandleCreate (rows, int, h);

  /* fill the indices with row numbers */
  for (i = 0; i < rows ; ++i)
    array (indexArray, i, int) = i ;

  /* set global table-var used by compare() to be this one */
  sortTable = t;

  /* rearrange lines */
  arraySort (indexArray, acTableOrder) ;
  
  /**** use index to rearrange column by column ****/

  {
    BitSet bs = 0 ;
    int row, col, k, k2 ;
    AC_CELL tmp ;
    int n = sizeof (AC_CELL) ;

    for (col = 0 ; col < t->cols ; col++)
      { 
	bs = bitSetReCreate (bs, t->rows) ;
	for (row = 0 ; row < t->rows ; row++)
	  if (!bit(bs,row))         /* not done yet */
	    {
	      bitSet(bs,row) ;
	      if (arr (indexArray, row, int) != row)
		{ 
		  memcpy (&tmp, ac_table_cell (t, row, col), n) ;
		  for (k = row, k2 = arr (indexArray, k, int) ;
		       !bit(bs, k2) ;
		       k = k2, k2 = arr (indexArray, k, int)
		       )
		    { 
		      memcpy (ac_table_cell (t, k, col), ac_table_cell (t, k2, col), n) ;
		      bitSet(bs,k) ;
		    }
		  memcpy (ac_table_cell (t, k, col), &tmp, n) ;
		}
	    }
      }
    
    bitSetDestroy (bs);
  }

  ac_free (h);

  return TRUE;
} /* ac_table_do_Sort */

/************************************************************/
/* Sort a table according to a specification
   the sort specifier :
   spec should look like "-2+4+1" which would specify first priority
   column 2 reversed, then column 4, then column 1.

   returns FALSE if table is wrong, or spec is out of range
*/
BOOL ac_table_sort (AC_TABLE table, const char *spec, vTXT errTxt)
{
  return ac_table_do_sort (table, spec) ;
} /* ac_table_sort */

/************************************************************/

AC_TABLE ac_table_paragraph (AC_TABLE table, int *rowp, AC_HANDLE h) 
{
  AC_TABLE tbl = 0 ;
  int ir, ir0, irMax ;
  int col, colMax ;
  struct ac_table_internal *t2, *t = (struct ac_table_internal *)table ;
  AC_CELL *c0, *c ;

  if (! table)
    messcrash ("ac_table_paragraph received a null table") ;
  if (t->magic != MAGIC_BASE + MAGIC_AC_TABLE)
    messcrash ("tableSort() received non-magic TABLE* pointer") ;
  if (! rowp)
    messcrash ("ac_table_paragraph received a null rowp") ;
  if (*rowp < 0)
    messcrash ("tableSort() received negative row value %d", *rowp) ;

  if (*rowp >= table->rows)
    return 0 ;

  irMax = ir0 = *rowp ;
  colMax = table->cols ;
  c0 = ac_table_cell (t, ir0, 0) ;
  while (c0->u.object == ac_table_cell (t, irMax, 0)->u.object)
    irMax++ ;
  
  tbl = ac_db_empty_table (ac_table_db (table), irMax - ir0, table->cols, h) ;
  t2 = (struct ac_table_internal *)tbl ;
  for (ir = ir0 ; ir < irMax ; ir++)
    for (col = 0 ; col < colMax ; col++)
      {
	c = ac_table_new_cell (t2, ir - ir0, col) ;
	c0 = ac_table_cell (t, ir, col) ;
	*c = *c0 ;
      }
  *rowp = irMax ;
  return tbl ;
}

/************************************************************/
/* Sort and compress a table
 * Eliminate succesive identical lines
 * The resulting table will have no duplicates only if it is fully sorted 
 * Hence the spec is completed with all remaing columns in natural order 
 * Example  :
 *   spec =  -2+5  in a table with 7 column is completed as
 *           -2+5+1+3+4+6
 */
BOOL ac_table_sort_compress (AC_TABLE table, const char *spec, vTXT errTxt)
{
  AC_HANDLE h ;
  int i, i1, j, k, n, cols, rows, specLen, sortCol ; 
  const char *ccp ;
  BOOL isPlus = TRUE ;
  struct ac_table_internal *t = (struct ac_table_internal *)table ;
  vTXT vSpec = 0 ;

  if (!t) return FALSE; 
  if (t->magic != MAGIC_BASE + MAGIC_AC_TABLE)
    messcrash ("tableSort() received non-magic TABLE* pointer") ;
  rows = t->rows ;
  if (rows <= 1) return TRUE; /* 1 or no elements are always sorted */

  h = ac_new_handle () ;
  vSpec = vtxtHandleCreate (h) ;

  rows = t->rows ;

  /* set sortTable and spec */
  if (!spec)
    {
      for (i = 1 ; i <= t->cols ; i++)
	vtxtPrintf (vSpec, "+%d", i) ;
    }
  else
    vtxtPrintf (vSpec, "%s", spec) ;
  spec = vtxtPtr (vSpec) ;

  specLen = strlen(spec);
  sortSpecList = arrayReCreate(sortSpecList, (specLen), int);

  for (i = j = 0, ccp = spec ; *ccp ; i++)
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
      if (sortCol < 1 || sortCol > t->cols) 
	/* column numbers have to range from 1 to colNum */
	continue ;
      if (isPlus)
	array(sortSpecList, j, int) = sortCol;
      else 
	array(sortSpecList, j, int) = -sortCol;
      j++ ;
    }

  /* if spec does not contain all columns add them in their natural order */
  for (i = 1 ; i <= t->cols ; i++)
    {
      BOOL found = FALSE ;
      for (j = 0 ; ! found && j < arrayMax (sortSpecList) ; j++)
	{
	  k = arr (sortSpecList, j, int) ;
	  if (i == k || i == -k)
	    found = TRUE ;
	}
      if (! found)
	vtxtPrintf (vSpec, "+%d", i) ;
    }
  spec = vtxtPtr (vSpec) ;

  /* sort the table completely */
  ac_table_do_sort (table, spec) ;

  /* remove duplicate lines */
  /* sortTable and sortSpecList are positionned by ac_table_do_sort */
  rows = t->rows ;
  cols = t->cols ;
  n = sizeof (AC_CELL) ;
  for (i = 0, j = -1 ; i < rows ; i++)
    {
      if (j >= 0 &&
	  ! acTableOrder (&i, &j)
	  )
	continue ;  /* candidate line j is identical to line i */
      j++ ;
      if (i > j) /* transfer *j to *i */
	{
	  for (k = 0 ; k < cols ; k++)
	    memcpy (ac_table_cell (t, j, k), ac_table_cell (t, i, k), n) ;
	}
    }
  t->rows = j + 1 ;

  return TRUE ;
} /* ac_table_sort_compress */

/************************************************************/
/************************************************************/

