
#define DEFINE_OBJ
typedef struct sobj *OBJ ;

#include <wh/acedb.h>
#include <wh/bs_.h>
#include <wh/command.h>
#include <wh/table.h>
#include <wh/dna.h>
#include <wh/peptide.h>
#include <wh/parse.h>
#include <wh/spread_.h>
#include <wh/spread.h>
#include <wh/bitset.h>

#include <wac/ac.h>
#include <wac/ac_.h>
#include <wac/acinside_.h>

#define ACE_NCBI
static void ac_finalise (void *);

BOOL acCaseSensitive (KEY key)
{
  return pickCaseSensitive (key) ;
}

/********************************************************************/
/********************************************************************/
/*
 * connect to the internal existing database 
 */

static int ac_open_counter = 0;
/*
 * counter to know when we should call aceQuit().
 */

int ac_db_nActiveClients (AC_DB db)
{ /* not public on server side */
  return -1 ;
}


static void ac_close_on_exit()
{
  if (ac_open_counter > 0)
    aceQuit(TRUE);
}

AC_DB ac_open_db (const char *database, const char **error)
{
  AC_DB db;

  /*
   * thisSession.session > 0 if the database is already open
   */
  if (thisSession.session != 0)
    {
      /*
       * If it is initialized already, we are embedded in tace or something
       * similiar, and we do not want to initialize the database.
       */
    }
  else
    {
      /*
       * If it is NOT initialized already, we are running in a standalone
       * Ace C program, and we need to initialize the database and see
       * that we de-initialize it later.
       */
      if (ac_open_counter == 0)
	{
	  aceInit (database);
	  atexit (ac_close_on_exit);
	}
      ac_open_counter++;
    }

  db = (AC_DB) halloc (sizeof (struct ac_db), 0) ;
  blockSetFinalise (db, ac_finalise) ;

  if (error) *error = 0 ; /* internal opening cannot fail */
  db->magic = MAGIC_BASE + MAGIC_AC_DB ;
  db->handle = ac_new_handle () ;
  db->command_stack = stackHandleCreate (1000, db->handle) ;
  db->look = aceCommandCreate (255, 0xff) ;

  return db ;
}

/********************************************************************/

static void ac_finalise (void *v)
{
  int magic = v ? *(int *)v : 0 ;
  int type = v && magic ? magic - MAGIC_BASE : 0 ;

  if (v && magic)
    switch (type)
      {
      case MAGIC_AC_DB:
	{
	  AC_DB db = (AC_DB) v ;
	  db->magic = 0 ;
	  aceCommandDestroy (db->look) ;
	  handleDestroy (db->handle) ;
	  if (ac_open_counter > 0)
	    {
	      ac_open_counter--;
	      if (ac_open_counter == 0)
		aceQuit(TRUE);
	    }
	}
	break ;
	
      case MAGIC_AC_KEYSET:
	{
	  AC_KEYSET aks = (AC_KEYSET) v ;
	  aks->magic = 0 ;
	  keySetDestroy (aks->ks) ;
	}
	break ;
	
      case MAGIC_AC_ITERATOR:
	{
	  AC_ITER iter = (AC_ITER)v ;
	  iter->magic = 0 ;
	  keySetDestroy (iter->ks) ;
          /* iter->ac_obj belongs to the client */
	}
	break ;
	
      case MAGIC_AC_OBJECT: /* AC_OBJ */
	{
	  AC_OBJ aobj = (AC_OBJ)v ;
	  aobj->magic = 0 ;
	  bsDestroy (aobj->obj) ;
	}
	break ;      
	
      default:
	messcrash ("Bad magic in ac_free, type = %d", type) ;
      }
  return ;
}

/********************************************************************/
/********************************************************************/
/*
 * raw command
 */

static unsigned char * ac_command_binary (AC_DB db, const char *command, int *response_length, AC_HANDLE handle, BOOL bin)
{
  unsigned char *result = 0 ;
  Stack s = 0 ;

  if (!db || db->magic != MAGIC_BASE + MAGIC_AC_DB)
    messcrash ("ac_command received a null or invalid db handle") ;
  if (!command || ! *command)
    return 0 ;
  s = commandStackExecute (db->look, command) ;

  if (bin)
    {
      char *cp =  stackText (s, 0) ;
      int n = stackMark(s) ;

      cp++; n-- ;
      while (n > 0 && *cp != '\n') {cp++; n-- ;}
      cp++; n-- ;
      if (*cp == 'c') 
	{ 
	  cp++ ;
	  memcpy (&n, cp, 4) ;
	  cp += 5 ;
	  if (n > 0)
	    { 
	      result = halloc (n+1, handle) ;
	      memcpy (result, cp, n) ; 
	      result[n] = 0 ;
	    }
	  else
	    n = 0 ;
	}
       else
	 n = 0 ;
      if (response_length)
	*response_length = n ;
    }
  else
    {
      if (response_length)
	*response_length = stackMark(s) - 1;
      result = (unsigned char *) strnew (stackText (s, 0), handle) ;
    }

  stackDestroy (s) ;
  return result ;
}

unsigned char * ac_command (AC_DB db, const char *command, int *response_length, AC_HANDLE handle)
{ 
  return ac_command_binary (db, command, response_length, handle, FALSE) ;
}

/********************************************************************/
/*
 * execute an ACEDB command and make a keyset of whatever the active
 * keyset is at the end of the command.  If iks is not NULL, it
 * is made the active keyset before executing the command.
 */
AC_KEYSET ac_command_keyset (AC_DB db, const char *command, AC_KEYSET iks,
			     AC_HANDLE handle) 
{
  AC_KEYSET aks = 0 ;
  unsigned char *cp = 0 ;
  
  if (!db || db->magic != MAGIC_BASE + MAGIC_AC_DB)
    messcrash ("ac_command_keyset received a null or invalid db handle") ;
  if (iks)
    aceCommandSwitchActiveSet (db->look, iks->ks, 0) ;

  cp = ac_command (db, command, 0, 0) ;
  messfree (cp) ;

  aks = ac_new_keyset (db, handle) ;
  aceCommandSwitchActiveSet (db->look, 0, aks->ks) ;

  return aks ;
}

/********************************************************************/
/********************* KeySet ***************************************/
/* note that ac_keyset is created on the handle 
   but not ak->ks which is private to this package
   and freed inside ac_finalise
*/

AC_KEYSET ac_new_keyset (AC_DB db, AC_HANDLE handle)
{
  KEYSET ks = keySetCreate () ;
  AC_KEYSET aks = (AC_KEYSET) halloc (sizeof (struct ac_keyset), handle) ;
  blockSetFinalise (aks, ac_finalise) ;

  if (!db || db->magic != MAGIC_BASE + MAGIC_AC_DB)
    messcrash ("ac_new_keyset received a null or invalid db handle") ;
  aks->db = db ;
  aks->magic = MAGIC_BASE + MAGIC_AC_KEYSET ;
  aks->ks = ks ;

  return aks ;
}

/********************************************************************/

AC_KEYSET ac_copy_keyset (AC_KEYSET aksold, AC_HANDLE handle)
{
  KEYSET ks = keySetCopy (aksold->ks) ;
  AC_KEYSET aks = (AC_KEYSET) halloc (sizeof (struct ac_keyset), handle) ;
  blockSetFinalise (aks, ac_finalise) ;

  aks->db = aksold->db ;
  aks->magic = MAGIC_BASE + MAGIC_AC_KEYSET ;
  aks->ks = ks ;

  return aks ;
} /* ac_keyset_copy */

/********************************************************************/
/* performs a query over the initial keyset,  creates a keyset of the result */
static AC_KEYSET ac_doquery_keyset (AC_DB db,  const char *queryText, AC_KEYSET initial_keyset,
			   AC_HANDLE handle )
{
  KEYSET ks = query (initial_keyset ? initial_keyset->ks : 0, queryText) ;
  AC_KEYSET aks = (AC_KEYSET) halloc (sizeof (struct ac_keyset), handle) ;
  blockSetFinalise (aks, ac_finalise) ;

  aks->db = db ;
  aks->magic = MAGIC_BASE + MAGIC_AC_KEYSET ;
  aks->ks = ks ;

  return aks ;
} /* ac_doquery_keyset  */

/* PUBLIC
 * performs a query, creates a keyset of the result
 */
AC_KEYSET ac_dbquery_keyset (AC_DB db, const char *queryText, AC_HANDLE handle)
{
  if (!db || db->magic != MAGIC_BASE + MAGIC_AC_DB)
    messcrash ("ac_dbquery_keyset received a null or invalid db handle") ;
  return ac_doquery_keyset (db, queryText, 0, handle) ;
}

AC_KEYSET ac_ksquery_keyset (AC_KEYSET initial_keyset, const char *queryText, AC_HANDLE handle)
{
  return ac_doquery_keyset (initial_keyset->db, queryText, initial_keyset, handle) ;
}

AC_KEYSET ac_table_keyset (AC_DB db, AC_TABLE atbl, int column,  AC_HANDLE handle)
{
  int ir, jr = 0 ;
  KEYSET ks = keySetCreate () ;
  AC_KEYSET aks = 0 ;
  KEY key ;

  if (atbl && column < atbl->cols)
    for (ir = 0 ; ir < atbl->rows ; ir++)
      {
	key = ac_table_key (atbl, ir, column, 0) ;
	if (key)
	  keySet (ks, jr++) = key ;
      }
  if (jr) 
    { keySetSort (ks) ; keySetCompress (ks) ; }
  
  aks = (AC_KEYSET) halloc (sizeof (struct ac_keyset), handle) ;
  blockSetFinalise (aks, ac_finalise) ;

  aks->db = db ;
  aks->magic = MAGIC_BASE + MAGIC_AC_KEYSET ;
  aks->ks = ks ;

  return aks ;
}

AC_KEYSET ac_objquery_keyset (AC_OBJ aobj, const char *queryText,  AC_HANDLE handle) 
{
  AC_KEYSET aks1 = ac_new_keyset (aobj->db, handle) , aks2 ;

  if (aobj)
    ac_keyset_add (aks1, aobj) ;
  else
    messcrash ("ac_objquery_keyset called with null obj") ;

  if (queryText)
    {
      aks2 = ac_doquery_keyset (aobj->db, queryText, aks1, handle) ;
      ac_free (aks1) ;

      return aks2 ;
    }
  else
    return aks1 ;
}

/********************************************************************/
/*
 * read/write keyset in non-volatile storage.  In the current
 * client implementation, name is a file name on the server's
 * disk.
 * ac_write writes a standard ace file, 
 *   the name of the file becomes the name of the keyset
 * ac_read only imports the recognised keys
 *
 * ac_read_keyset returns NULL if there is no file with that name.
 */

/* returns NULL if there is no file with that name.
 * imports the recognised keys otherwise
 */
AC_KEYSET ac_read_keyset (AC_DB db, const char *fname, AC_HANDLE handle)
{ 
  AC_KEYSET aks = 0  ;
  KEYSET ks = 0 ;
  int i ;

#ifdef ACE_NCBI
  FILE *f = 0 ;
  int level ;

  if (!db || db->magic != MAGIC_BASE + MAGIC_AC_DB)
    messcrash ("ac_read_keyset received a null or invalid db handle") ;

  if (!fname || !*fname ||
      ! (f = filopen (fname, 0, "r")))
    {
      freeOutf ("// Sorry, I could not open file %s\n", fname ? fname : "NULL") ;
    }
  else
    {
      level = freesetfile (f, "") ;
      ks = keySetRead (level) ; /* will close level */
      freeclose (level) ; /* redundant but more elegant */
      
      aks = ac_new_keyset (db, handle) ;
      i = keySetMax (ks) ;
      if (i) while (i--) 
	keySet (aks->ks, i) = keySet (ks, i) ;
      keySetDestroy (ks) ;
    }
#else
  ACEIN fi = 0 ;

  if (!db || db->magic != MAGIC_BASE + MAGIC_AC_DB)
    messcrash ("ac_read_keyset received a null or invalid db handle") ;

  if (!fname || !*fname ||
      ! (fi = aceInCreateFromFile (fname, "r", 0, 0)))
    {
      freeOutf ("// Sorry, I could not open file %s\n", fname ? fname : "NULL") ;
    }
  else
    {
      Stack errors = stackCreate (1000) ;
      ACEOUT fo = aceOutCreateToStack (errors, 0) ;

      ks = keySetRead (fi, fo) ;
      aceOutDestroy (fo) ;
      stackDestroy (errors) ;  /* BUG: i should possibly look at the errors */
      aceInDestroy (fi) ;

      aks = ac_new_keyset (db, handle) ;
      i = keySetMax (ks) ;
      if (i) while (i--) 
	keySet (aks->ks, i) = keySet (ks, i) ;
      keySetDestroy (ks) ;
    }
#endif
    
  return aks ;
}

/********************************************************************/
/* ac_keyset_write writes a standard ace file, 
 *   the name of the file becomes the name of the keyset
 *   returns FALSE if !ks OR if the file cannot be opened
 */

BOOL ac_keyset_write (AC_KEYSET aks, const char *name)
{ 
  KEYSET ks = aks->ks ;

#ifdef  ACE_NCBI
  FILE *f = filopen (name, 0, "w") ;
  
  if (name && *name && aks && 
      (f = filopen (name, 0, "w")))
    {
      char *cp, *name2 ;

      name2 = strnew (name, 0) ;
      cp = name2 + strlen (name2) - 1 ;
      while (cp > name2 && *cp != '/')
	cp-- ;
      if (*cp=='/') cp++ ;
      fprintf(f, "KEYSET %s\n", cp) ;
      if (keySetMax (ks))
	keySetDump (f, 0, ks) ;
      filclose (f) ;
      messfree (name2) ;
      return TRUE ;
    }
#else
  ACEOUT fo = 0 ;
  if (name && *name && aks && 
      (fo = aceOutCreateToFile (name, "w", 0)))
    {
      const char *cp = name + strlen (name) - 1 ;
      while (cp > name && *cp != '/')
	cp-- ;
      if (*cp=='/') cp++ ;
      aceOutf (fo, "KEYSET %s\n", cp) ;
      if (keySetMax (ks))
	keySetDump (fo, 0, ks) ;
      aceOutDestroy (fo) ;
      return TRUE ;
    }
#endif
  return FALSE ;
}

/********************************************************************/

int ac_keyset_and (AC_KEYSET aks1, AC_KEYSET aks2 )
{
  KEYSET ks = keySetAND (aks1->ks, aks2->ks) ;

  keySetDestroy (aks1->ks) ;
  aks1->ks = ks ;

  return keySetMax (ks) ;
} /* ac_keyset_and */

/********************************************************************/

int ac_keyset_or (AC_KEYSET aks1, AC_KEYSET aks2 )
{
  KEYSET ks = keySetOR (aks1->ks, aks2->ks) ;

  keySetDestroy (aks1->ks) ;
  aks1->ks = ks ;

  return keySetMax (ks) ;
} /* ac_keyset_or */

/********************************************************************/

int ac_keyset_xor (AC_KEYSET aks1, AC_KEYSET aks2)
{
  KEYSET ks = keySetXOR (aks1->ks, aks2->ks) ;

  keySetDestroy (aks1->ks) ;
  aks1->ks = ks ;

  return keySetMax (ks) ;
} /* ac_keyset_xor */

/********************************************************************/

int ac_keyset_minus (AC_KEYSET aks1, AC_KEYSET aks2 )
{
  KEYSET ks = keySetMINUS (aks1->ks, aks2->ks) ;

  keySetDestroy (aks1->ks) ;
  aks1->ks = ks ;

  return keySetMax (ks) ;
} /* ac_keyset_minus */

/********************************************************************/

BOOL ac_keyset_add (AC_KEYSET aks, AC_OBJ aobj)
{
  KEY key = aobj->key ;

  return keySetInsert (aks->ks, key) ;
} /* ac_keyset_add */

/********************************************************************/

BOOL ac_keyset_contains (AC_KEYSET aks, AC_OBJ aobj)
{
  KEY key = aobj->key ;

  return keySetFind (aks->ks, key, 0) ;
} /* ac_keyset_contains */

/********************************************************************/

BOOL ac_keyset_remove (AC_KEYSET aks, AC_OBJ aobj)
{
  int i = 0, imax ;
  KEY *kp = 0, key = aobj->key ;
  BOOL found = FALSE ;
  
  imax = keySetMax (aks->ks) ;
  if (imax)
    for (i = 0, kp = arrp (aks->ks, 0, KEY) ; i < imax ; kp++, i++)
      if (key == *kp) { found = TRUE ; break ; }
  if (found)
    {
      imax -= 1 ;
      for (;i < imax ; i++, kp++)
	*kp = *(kp+1) ;
      aks->ks->max --;
    }
  return found ;
} /* ac_keyset_remove */

/********************************************************************/
/*
 * returns a subset strating at x0 of length nx
 * index start at zero and are ordered alphabetically
 * which is usefull to get slices for display
 */
AC_KEYSET ac_subset_keyset (AC_KEYSET aks, int x0, int nx, AC_HANDLE handle) 
{
  AC_KEYSET aks1 = 0  ;
  KEYSET kA ;
  int i, j ;

  aks1 = ac_new_keyset (aks->db, handle) ;
  if (aks->ks)
    {
      if (x0 < 1) 
	x0 = 1 ;
      if (nx > keySetMax (aks->ks) - x0 + 1) 
	nx = keySetMax (aks->ks) - x0 + 1 ; 
      kA = keySetAlphaHeap(aks->ks, x0 + nx - 1) ;

      for (i = x0 - 1, j = 0 ; nx > 0 ; i++, j++, nx--)
	keySet (aks1->ks, j) = keySet (kA, i) ;
      keySetSort (aks1->ks) ;
    }

  return aks1 ;
} /* ac_subset_keyset */

/********************************************************************/

int ac_keyset_count (AC_KEYSET aks)
{
  return  aks && keySetExists (aks->ks) ? keySetMax (aks->ks) : 0 ;
}

/********************************************************************/
/********************************************************************/
/* Keyset Iterator functions */

AC_ITER ac_query_iter (AC_DB db, int fillhint, const char *query, AC_KEYSET initial_keyset, AC_HANDLE handle) 
{ 
  KEYSET ks = 0 ;
  AC_ITER iter = (AC_ITER) halloc (sizeof (struct ac_iter), handle) ;
  blockSetFinalise (iter, ac_finalise) ;

  if (!db || db->magic != MAGIC_BASE + MAGIC_AC_DB)
    messcrash ("ac_query_iter received a null or invalid db handle") ;
  iter->db = db ;
  iter->magic = MAGIC_BASE + MAGIC_AC_ITERATOR ;
  ks = query (initial_keyset ? initial_keyset->ks : 0, query) ; 
  iter->ks = ks && keySetMax (ks) ? keySetAlphaHeap (ks, keySetMax (ks)) : keySetCreate () ;
  iter->curr = 0 ;
  iter->max = keySetMax (iter->ks) ;

  keySetDestroy (ks) ;
  return iter ;
} /* ac_query_iter */

/****************************************************************/
/* PUBLIC ac_dbquery_iter ()
 * make a query, dump resulting keyset directly into an iterator
 * Short hand 1 of previous function, fillhint = true by default
 *	db = database to query
 *	query = query to use
 */
AC_ITER ac_dbquery_iter (AC_DB db, const char *query, AC_HANDLE handle)
{
  return ac_query_iter (db, TRUE, query, 0, handle) ;
} /* ac_dbquery_iter */

/****************************************************************/
/* PUBLIC ac_ksquery_iter ()
 * make a query, dump resulting keyset directly into an iterator
 * Short hand 2 of previous function, fillhint = true by default
 *	ks = keyset to query
 *	query = query to use
 */
AC_ITER ac_ksquery_iter (AC_KEYSET aks, const char *query, AC_HANDLE handle)
{
  if (aks)
    return ac_query_iter (aks->db, TRUE, query, aks, handle) ;
  return 0 ;
} /* ac_ksquery_iter */

/****************************************************************/
/* PUBLIC ac_objquery_iter ()
 * make a query, dump resulting keyset directly into an iterator
 * Short hand 3 of previous function, fillhint = true by default
 *	obj = single object on which we initialise the query
 *	query = query to use
 */
AC_ITER ac_objquery_iter (AC_OBJ aobj, const char *query, AC_HANDLE handle)
{
  AC_ITER iter = 0 ;
  if (aobj)
    {
      AC_KEYSET aks = ac_objquery_keyset (aobj, "IS *",  0) ;
      iter = ac_query_iter (aobj->db, TRUE, query, aks, handle) ;
      ac_free (aks) ;
    }
  return iter ;
} /* ac_objquery_iter */

/********************************************************************/

AC_ITER ac_keyset_iter (AC_KEYSET aks, int fillhint, AC_HANDLE handle) 
{ 
  AC_ITER iter = (AC_ITER) halloc (sizeof (struct ac_iter), handle) ;
  blockSetFinalise (iter, ac_finalise) ;

  iter->db = aks->db ;
  iter->magic = MAGIC_BASE + MAGIC_AC_ITERATOR ;
  iter->ks = aks->ks && keySetMax (aks->ks) ? keySetAlphaHeap (aks->ks, keySetMax (aks->ks)) : keySetCreate () ;
  iter->curr = 0 ;
  iter->max = keySetMax (iter->ks) ;

  return iter ;
} /* ac_keyset_iter */

/********************************************************************/
/* Causes the datase to save its state to disk
 * useful in very long client codes to ensure that the work so far cannot be lost
 */
void ac_db_commit (AC_DB db)
{
  sessionClose (TRUE) ;
} /* ac_db_commit */

/********************************************************************/
/* refresh, is declared by compatibility with acclient */
void ac_db_refresh (AC_DB db) {} ; /* this function is not needed in acinside
				      since bsCreate() always get the current object
				      for the acedb cache 
				   */

AC_OBJ ac_key2obj (AC_DB db, KEY key, BOOL fillIt, AC_HANDLE handle)
{
  AC_OBJ aobj = (AC_OBJ) halloc (sizeof (struct ac_object), handle) ;
  blockSetFinalise (aobj, ac_finalise) ;

  aobj->db = db ;
  aobj->magic = MAGIC_BASE + MAGIC_AC_OBJECT ;
  aobj->key = key ;
  if (fillIt) aobj->obj = bsCreate (key) ;
  aobj->isEmpty = aobj->obj ? FALSE : TRUE ;

  return aobj ;
} /* ac_key2obj */

/********************************************************************/
/* fetch the next object in the keyset being iterated over */

AC_OBJ ac_next_obj (AC_ITER iter)
{
  /* do not ac_free (iter->ac_obj) ; it belongs to the application */
  if (iter->curr < iter->max)
    {
      iter->ac_obj = ac_key2obj (iter->db, keySet (iter->ks, iter->curr), TRUE, 0) ;
      iter->curr++ ;
      return iter->ac_obj ;
    }
  return 0 ;
} /* ac_next */

/********************************************************************/
/* resets the iterator to its top */
BOOL ac_iter_rewind (AC_ITER iter) 
{ 
  if (iter && iter->max)
    {
      ac_free (iter->ac_obj) ;
      iter->curr = 0 ;
      return TRUE ;
    }
  return FALSE ;
} /* ac_rewind  */

/********************************************************************/
/********************************************************************/
/* object related functions */

/********************************************************************/
#include <sysclass.h>
KEY ac_get_key (AC_DB db, const char *classe, const char *nam)
{
  KEY cl = 0, key = 0 ;
  
  if (!db || db->magic != MAGIC_BASE + MAGIC_AC_DB)
    messcrash ("ac_get_key received a null or invalid db handle") ;
    
  if (lexword2key(classe, &cl, _VMainClasses))
    lexword2key (nam, &key, KEYKEY(cl)) ;

  return key ;  
} /* ac_get_key */

/********************************************************************/

AC_OBJ ac_get_obj (AC_DB db, const char *class, const char *nam, AC_HANDLE handle)
{
  AC_OBJ aobj = 0 ;
  char buf[4096] ;
  KEYSET ks = 0 ;

  if (!db || db->magic != MAGIC_BASE + MAGIC_AC_DB)
    messcrash ("ac_get_obj received a null or invalid db handle") ;
  if (strlen(nam) < 3000)
    {
      sprintf (buf, "Find %s IS \"%s\"", class, nam) ;
      ks = query (0, buf) ;
    }
  else /* never use messprintf in this library */
    {
      char *cp = messalloc (30 + strlen(class) + strlen(nam)) ;
      sprintf (cp, "Find %s IS \"%s\"", class, nam) ;
      ks = query (0, cp) ;
      messfree (cp) ;
    }
  if (keySetMax (ks))
    aobj = ac_key2obj (db, keySet (ks, 0), FALSE, handle) ;
  keySetDestroy (ks) ;

  return aobj ;  
} /* ac_get_obj */

/********************************************************************/

AC_OBJ ac_copy_obj (AC_OBJ aobj, AC_HANDLE handle)
{
  return aobj ? ac_key2obj (aobj->db, aobj->key, FALSE, handle) : 0 ;
} /* ac_dup */

/***************************************************************************************/
/* returns the first obj satisfying the query, starting on obj:src
 * remember that acedb offers NO garantee on the order of the objects
 * following a tag, so the result of this query is random
 * if there are severeal solutions to the query
 * yet, this call is very convenient if the query is unambiguous
 * for example if the query 'follows' a UNIQUE tag.
 */
AC_OBJ ac_objquery_obj (AC_OBJ src, const char *question, AC_HANDLE h)
{
  AC_ITER iter = ac_objquery_iter (src, question, h) ;
  AC_OBJ obj = 0, objCopy = 0 ;
  
  if (iter)
    {
      obj = ac_iter_obj (iter) ;
      if (obj)
	{
	  objCopy = ac_copy_obj (obj, h) ;
	  ac_free (obj) ;
	}
      ac_free (iter) ;
    }
  return objCopy ;
} /* ac_objquery_obj */

/********************************************************************/

const char *ac_key_class (KEY key)
{
  return className (key) ;
} /* ac_name */

/********************************************************************/

const char *ac_key_name (KEY key)
{
  return name (key) ;
} /* ac_name */

/********************************************************************/

const char *ac_new_name (AC_DB db, const char *class_nam, const char *format, AC_HANDLE handle)
{
  AC_HANDLE h = ac_new_handle () ;
  unsigned char *buf ;
  char *command, *cp, *cq ;
  const char *ccp = "NULL" ;

  cp = strnew (format, h) ;
  cq = strstr(cp, "%") ;
  if (cq) 
    {
      *cq = 0 ;
      command = hprintf (h, "New %s %s%s%s", class_nam, cp, "\\%", cq+1) ;
    }
  else
    command = hprintf (h, "New %s ", class_nam) ;
  buf = ac_command (db, command, 0, h) ;
  if (
      (cp = strstr ((char*)buf, "\"")) &&
      (cq = strstr (cp+1, "\""))
      )
    {
      *cq = 0 ;
      ccp = strnew (cp+1, handle) ; /* success */
    }
  ac_free (h) ;

  return ccp ;
} /* ac_new_name */
/* new unique objects */
#ifdef JUNK
we used to have a fancy code, verify if it was transfered to command.c
AceKey aceMakeKeyFromFormat(char* classname, char* format)
{
  char *cp = format ;
  int i, classe = -1 ;
  KEY key, kk ;
  static Stack localStack = 0 ;

  classe = aceClassNumber(classname) ;
  if (classe < 0)
    { 
      aceSetErrorMessage(messprintf ("// Bad classname %s", classname)) ;
      return 0 ;
    }
  if (pickList[classe].protected)
    { 
      aceSetErrorMessage(messprintf ("// Class %s is protected", classname)) ;
      return 0 ;
    }

  if (!cp) 
    { 
      aceSetErrorMessage(messprintf ("// Usage New Class Format: Construct a new name"
				     " according to C-language printf(format,n) \n"
				     "// Example: New Oligo xxx.\\%%d\n")) ;
      return 0 ;
    }
  if (!aceGetWriteAccess (FALSE) )
    { 
      aceSetErrorMessage("// Sorry, you do not have write access") ;
      return 0 ;
    }
  localStack = stackReCreate(localStack, 80) ;
  pushText(localStack, cp) ;
  
  i = -2 ; /* important flag */
  cp = stackText(localStack,0) ;
  while (*cp && *cp != '%')
  cp++ ;
  if (*cp++ == '%')
    { 
      while(*cp && isdigit((int)*cp)) cp++ ;  
      if (*cp++ == 'd')
	{ i = 1 ;
	while (*cp && *cp != '%') cp++ ;
	if (!*cp)
	  goto okformat ;
	}
      aceSetErrorMessage ("// Only allowed format should contain zero or once %%[number]d") ;
      return 0 ;
    }
okformat:
  key = kk = 0 ;
  cp = stackText(localStack,0) ;
  while (i < lexMax(classe) + 3 && /* obvious format obvious problem */
	 !(kk = lexaddkey(messprintf(cp,i), &key, classe)))
    if (i++ < 0) break ;  /* so i break at zero f format has no %d */
  if (kk)
    return key ;
  aceSetErrorMessage (messprintf ("Name %s allready exists\n", cp)) ;
  return 0 ;
}
#endif
/********************************************************************/

KEY ac_obj_key (AC_OBJ aobj)
{
  if (aobj)
    return aobj->key ;
  else
    return 0 ;
}

/********************************************************************/

const char *ac_name (AC_OBJ aobj)
{
  return aobj ? name (aobj->key) : "(NULL)" ;
} /* ac_name */

/********************************************************************/
   
const char *ac_class (AC_OBJ aobj)
{
  return aobj ?  className (aobj->key) : "(NULL)" ;
} /* ac_class */

/********************************************************************/
   
BOOL ac_obj_equal (AC_OBJ aobj1, AC_OBJ aobj2)
{
  if (!aobj1  && !aobj2)
    return TRUE ;

  if (
      (!aobj1  && aobj2) ||
      (!aobj2  && aobj1) ||
      (aobj1->key != aobj2->key)
      )
    return FALSE ;

  return TRUE ;
} /* ac_obj_equal */

/********************************************************************/
   
int ac_obj_cmp (AC_OBJ aobj1, AC_OBJ aobj2)
{
  void *a, *b ;
  
  if (!aobj1 || ! aobj2) 
    return 0 ;
  a = assVoid (aobj1->key) ;
  b = assVoid (aobj2->key) ;

  return keySetAlphaOrder(&a, &b) ;
} /* ac_obj_cmp */

/********************************************************************/

BOOL ac_has_tag	(AC_OBJ aobj, const char * tagName)
{
  KEY tag = 0 ;
  if (!lexword2key (tagName, &tag, 0))
    return FALSE ;
  
  if (!aobj)
    return FALSE ;
  else if (aobj->obj)
    return bsFindTag (aobj->obj, tag) ;
  else
    return keyFindTag (aobj->key, tag) ;
}  /* ac_has_tag */

/********************************************************************/

static BOOL ac_getTag (AC_OBJ aobj, const char *tagName, KEY *tagp)
{
  if (!lexword2key (tagName, tagp, 0))
    return FALSE ;
  if (!aobj || (!aobj->obj && !keyFindTag (aobj->key, *tagp)))
    return FALSE ;
  if (!aobj->obj)
    aobj->obj = bsCreate (aobj->key) ;
  if (aobj->obj)
    return bsFindTag (aobj->obj, *tagp) ;
  
  return FALSE ; 
} /* ac_getTag */

/********************************************************************/
/*
 * returns the number of lines just behind a tag or zero
 * if the tag is not there or there is no data item there
 */
int ac_tag_count (AC_OBJ aobj, const char *tagName)
{
  KEY tag = 0 ;
  int x = 0 ;
  KEY k = 0 ;

  if (aobj && ac_getTag (aobj, tagName, &tag))
    if (bsGetKeyTags (aobj->obj, _bsRight, &k))
      do { x++ ; } while (bsGetKeyTags (aobj->obj, _bsDown, &k)) ;
  return x ;
} /* ac_tag_count */

/********************************************************************/

ac_type ac_tag_type (AC_OBJ aobj, const char *tagName)
{
  /*
   * returns the type of the data after the tag, or ac_none
   * if the tag is not there or there is no data item there
   */ 
  KEY tag = 0, key = 0 ;
  
  if (aobj &&  ac_getTag (aobj, tagName, &tag) &&
      (key = bsType (aobj->obj, _bsRight)))
    {
      if (class(key) == _VText) return ac_type_text ;
      if (class(key)) return ac_type_obj ;
      if (key == _Int) return ac_type_int ;
      if (key == _Float) return ac_type_float ;
      if (key == _DateType) return ac_type_date ;
      if (key == _Text) return ac_type_text ;
      if (key) return ac_type_tag ;

    }
  return ac_type_empty ;
} /* ac_tag_type */

/********************************************************************/
/*
 * returns integer after the tag, or dflt if the tag
 * is not there or the data is not an integer
 * (Should it convert float?)
 */
int ac_tag_int (AC_OBJ aobj, const char *tagName, int dflt)
{
  KEY tag = 0 ;
  int x = dflt ;

  if (aobj && ac_getTag (aobj, tagName, &tag))
    bsGetData (aobj->obj, tag, _Int, &x) ;
  return x ;
} /* ac_tag_int */

/********************************************************************/

float ac_tag_float (AC_OBJ aobj, const char *tagName, float  dflt)
{
  KEY tag = 0 ;
  float x = dflt ;

  if (aobj && ac_getTag (aobj, tagName, &tag))
    bsGetData (aobj->obj, tag, _Float, &x) ;
  return x ;
} /* ac_tag_float */

/********************************************************************/

const char *ac_tag_text (AC_OBJ aobj, const char *tagName, const char *dflt)
{
  KEY tag = 0, key = 0 ;
  char *cp ;
  
  if (aobj && ac_getTag (aobj, tagName, &tag))
    {
      if (bsGetKey (aobj->obj, tag, &key) &&
	  class (key) == _VText)
	return name(key) ;
      if (bsGetData (aobj->obj, tag, _Text, &cp))
	return cp ;
    }
  return dflt ;
} /* ac_tag_text */

/********************************************************************/

mytime_t ac_tag_date (AC_OBJ aobj, const char *tagName, mytime_t dflt)
{
  KEY tag = 0 ;
  mytime_t x = dflt ;
  
  if (aobj && ac_getTag (aobj, tagName, &tag))
    bsGetData (aobj->obj, tag, _DateType, &x) ;
  return x ;
} /* ac_tag_float */

/********************************************************************/

KEY ac_tag_key (AC_OBJ aobj, const char *tagName, KEY dflt)
{
  KEY tag = 0, key = 0 ;

  if (aobj && ac_getTag (aobj, tagName, &tag) &&
      bsGetKeyTags (aobj->obj, tag, &key))
    return key ;
  return dflt ;
} /* ac_tag_key */

/********************************************************************/

AC_OBJ ac_tag_obj (AC_OBJ aobj, const char *tagName, AC_HANDLE handle)
{
  KEY tag = 0, key = 0 ;

  if (aobj && ac_getTag (aobj, tagName, &tag) &&
      bsGetKey (aobj->obj, tag, &key))
    return ac_key2obj (aobj->db, key, FALSE, handle) ;
  return 0 ;
} /* ac_tag_obj */

/********************************************************************/

const char * ac_tag_printable( AC_OBJ aobj, const char *tagName, const char *dflt )
{
  /*
   * returns the string value of the data after the tag, or dflt
   * if the tag is not there or there is no data item there
   */ 
  KEY tag = 0, key = 0 ;
  
  if (aobj && ac_getTag (aobj, tagName, &tag) &&
      (key = bsType (aobj->obj, _bsRight)))
    {
      if (class(key))
	{
	  key = 0 ;
	  if (bsGetKey (aobj->obj, tag, &key))
	    return name(key);	/* it is an object */
	}
      else if (key == _Int)
	{ 
	  int i;
	  if (bsGetData (aobj->obj, tag, _Int, &i))
	    {
	      sprintf(aobj->buf25,"%d",i);
	      return aobj->buf25;
	    }
	}
      else if (key == _Float)
	{ 
	  float f;
	  if (bsGetData (aobj->obj, tag, _Float, &f))
	    {
	      sprintf(aobj->buf25,"%g",f);
	      return aobj->buf25;
	    }
	}
      else if (key == _DateType) 
	{ 
	  int i;
	  if (bsGetData (aobj->obj, tag, _DateType, &i))
	    return timeShow(i, aobj->buf25, 25);
	}
      else if (key == _Text) 
	{
	  char *cp;
	  if (bsGetData (aobj->obj, tag, _Text, &cp))
	    return cp ;
	}
      else if (key >= _Date) 
	{
	  if (bsGetKeyTags (aobj->obj, _bsRight, &key))
	    return name(key);	/* it is a tag */
	}

      /* bug: where is type bool ? */

    }
  return dflt;
} /* ac_tag_printable */

/********************************************************************/
/********************************************************************/
/*
 * exports an object in ace format
 * useful, in conjunction with ac_tag_ace to 
 * edit objects and recover all the data 
 * without having to know and control
 * all the details of the schema
 */
char *ac_obj_ace (AC_OBJ aobj, AC_HANDLE handle)
{
  char *cp = 0, *cq ;
  unsigned char *buf = 0 ;
  char *command ;
  const char *nam ;
  AC_HANDLE h ;

  if (aobj && iskeyold(aobj->key))
    {
      h = ac_new_handle () ;
      nam = className (aobj->key) ;

      command = hprintf (h, "query find %s IS \"%s\"", nam, name(aobj->key)) ;
      buf = ac_command (aobj->db, command, 0, h) ;
      buf = ac_command (aobj->db, "show -a", 0, h) ;
      if ((cp = strstr ((char*)buf, nam)) &&
	  (cq = strstr (cp, "\n\n"))) 
	{ /* clean up at end of paragraph */
	  *(cq+2) = 0 ;
	  cp = strnew (cp , handle) ; /* success */
	}
      ac_free (h) ;
    }

  return cp ;
} /* ac_obj_ace */

/********************************************************************/
/*
 * exports everything right of a tag in ace format
   * If the tag can be found in the object
   *  - The returned value is a char*, allocated on handle, 
        containing a valid piece to be used as an ace file
        but notice that the mandatory first line of
	any .ace paragraph (class   object_name) is not included
   * Otherwise:
   *  - returns NULL
 */
char *ac_tag_ace (AC_OBJ aobj, const char *tag, AC_HANDLE handle)
{
  char *cp = 0, *cq, *nam ;
  int nn ;
  KEY kTag = tag && *tag ? str2tag(tag): 0 ;
  unsigned char *buf = 0 ;
  char *command ;
  AC_HANDLE h ;

  if (kTag && aobj && keyFindTag (aobj->key, str2tag (tag))) 
    {
      h = ac_new_handle () ;
      nam = className (aobj->key) ;
      nn = strlen (nam) ;
      command = hprintf (h, "query find %s IS \"%s\"", nam, name(aobj->key)) ;
      buf = ac_command (aobj->db, command, 0, h) ;
      command = hprintf (h, "show -a %s", tag) ; /* just export this (multi-line) tag in .ace format */
      buf = ac_command (aobj->db, command, 0, h) ;
      if ((cp = strtok ((char *)buf, "\n"))) /* jump comments */
	do
	  {
	    if (!strncasecmp (cp, nam, nn) && *(cp + nn) == ' ' &&
		(cp = strtok (0, "")) &&
		(cq = strstr (cp, "\n\n"))) 
	      { /* clean up at end of paragraph */
		*(cq+1) = 0 ;
		cp = strnew (cp , handle) ; /* success */
		break ;
	      }
	  } while ((cp = strtok (0, "\n"))) ;
      ac_free (h) ;
    }

  return cp ;
} /* ac_tag_ace */

/********************************************************************/
/********************************************************************/
/* recursive  dump cell  routine */
static int ac_tagBS (AC_DB db, BS bs, BS bsm, AC_TABLE t, int row, int col)
{
  int oldrow ;

  while (bs)
    {
      oldrow = row ;
      bsModelMatch (bs, &bsm) ;

      /* dump right */
      if (bs->right)
	row = ac_tagBS (db, bs->right, bsModelRight(bsm), t, row, col + 1) ;
      for (; oldrow >= 0 && oldrow <= row ; oldrow++)
	{ /* fill in all the intermediate lines in my column */
	  if (bs->key <= _LastC)
	    ac_table_insert_text (t, oldrow, col, bsText(bs)) ;
	  else if (bs->key == _Int)
	    ac_table_insert_type (t, oldrow, col, &(bs->n), ac_type_int) ;
	  else if (bs->key == _Float)
	    ac_table_insert_type (t, oldrow, col, &(bs->n), ac_type_float) ;
	  else if (bs->key == _DateType)
	    ac_table_insert_type (t, oldrow, col, &(bs->n), ac_type_date) ;
#if 1
	  else
	    ac_table_insert_type (t, oldrow, col, &(bs->key), ac_type_key) ;
#else
	  else if (class(bs->key))
	    {
	      AC_OBJ aobj = ac_key2obj (db, bs->key, 0, 0) ;
	      ac_table_insert_type (t, oldrow, col, &(aobj), ac_type_obj) ;
	    }
	  else if (!class(bs->key)) /* a tag */
	    {
	      const char *cp = name(bs->key) ;
	      ac_table_insert_type (t, oldrow, col, &cp, ac_type_tag) ;
	    }
#endif
	}
      /* move down, dump in the cuurent while loop */
      bs = bs->down ;	;
      if (!bs)
	break ;	
      else
	row++ ;
    }
  return row ;
}


/** fetch a table from a tag in an object   */
AC_TABLE ac_tag_table (AC_OBJ aobj, const char *tagName,  AC_HANDLE handle)
{
  AC_TABLE atable = 0 ;

  if (tagName == NULL)
    {
      BS model = bsModelRoot (aobj->obj);
      atable = ac_db_empty_table (aobj->db, 1, 1, handle) ;
      if (!aobj->obj)
	aobj->obj = bsCreate (aobj->key) ;
      ac_tagBS( aobj->db, aobj->obj->root->right,
		bsModelRight(model), atable, 0, 0);
      return atable;
    }

  if (ac_tag_type (aobj, tagName) == ac_type_empty)
    return 0 ;

  atable = ac_db_empty_table (aobj->db, 1, 1, handle) ;

  if (aobj->obj->curr &&
      aobj->obj->curr->right) 
    ac_tagBS (aobj->db, aobj->obj->curr->right, bsModelRight(aobj->obj->modCurr), atable, 0, 0) ;
  
  return atable ;
} /* ac_tag_table */

/********************************************************************/
#ifdef JUNK
/* performs an AQL query */
static AC_TABLE ac_newAql (AC_DB db, AC_KEYSET aks, const char *myquery,  AC_HANDLE handle, const char* orderBy, const char **error_messages)
{
  AC_TABLE atable = 0 ;
  AQL aql = 0 ;
  TABLE *result = 0 ;
  int row, col, maxrow, maxcol ;
  const char *cp, *types ;

  if (!aks)
    {
      aql = aqlCreate2 (myquery, 0, 0, 0, 0, 0) ;
      result = aqlExecute (aql, 0, 0, 0, 0) ;
    }
  else
    {
      TABLE *activeTable = tableCreateFromKeySet (aks->ks ? aks->ks : keySetCreate(), 0);
      aql = aqlCreate2 (myquery, 0, 0, 0, 0, 1, "@active", 0) ; 
      result = aqlExecute (aql,  0, 0, 0, "TABLE @active", activeTable, 0);
      tableDestroy (activeTable);
    }

  if (result && !(aqlIsError(aql)))
    { 
      atable = ac_db_empty_table (db, maxrow = tableMax (result), maxcol = result->ncol, handle) ;
      types = tableTypes(result) ;
      if (orderBy)
	tableSort (result, orderBy) ;
      for (row = 0 ; row < maxrow ; row++)
	for (col = 0 ; col < maxcol ; col++)
	  {
 	    if (tabEmpty(result,row,col))
	      {
		/* if the cell is empty, we should not insert any data - it will remain ac_type_empty */
	      }
	    else
	      {
		switch (types[col])
		  {
#if 0
		  case 'k': case 'K':
		    {
		      AC_OBJ aobj = ac_key2obj (db, (tabKey (result, row, col)), 0, 0) ;
		      ac_table_insert_type (atable, row, col, &(aobj), ac_type_obj) ;
		      ac_free (aobj) ;
		    }
		    break ;
		  case 'g':
		    cp = name(tabTag (result, row, col)) ;
		    ac_table_insert_type (atable, row, col, &cp, ac_type_tag) ;
		    break ;
#endif
		  case 'k': 
		  case 'K':
		  case 'g':
		    ac_table_insert_type (atable, row, col, &(tabKey (result, row, col)), ac_type_key) ;
		    break ;
		  case 'b':
		    ac_table_insert_type (atable, row, col, &(tabBool (result, row, col)), ac_type_bool) ;
		    break ;
		  case 'i':
		    ac_table_insert_type (atable, row, col, &(tabInt (result, row, col)), ac_type_int) ;
		    break ;
		  case 'f':
		    ac_table_insert_type (atable, row, col, &(tabFloat (result, row, col)), ac_type_float) ;
		    break ;
		  case 't':
		    ac_table_insert_type (atable, row, col, &(tabDate (result, row, col)), ac_type_date) ;
		    break ;
		  case 's':
		    cp = tabString (result, row, col) ;
		    ac_table_insert_text (atable, row, col, cp) ;
		    break ;
		  }
	      }
	  }
    }
  else
    {
      /*
	if (error_messages)
	*error_messages = hprintf (handle, "\n// Aql: error %d : %s\n", 
	aqlGetErrorNumber(aql), aqlGetErrorMessage(aql)) ;
      */
    }
  aqlDestroy (aql) ;
  tableDestroy (result) ;

  return (AC_TABLE)atable ;
} /* ac_newAql */

/********************************************************************/

AC_TABLE ac_aql_table (AC_DB db, const char *myquery, AC_KEYSET initial_keyset, const char* orderBy, const char **error_messages, AC_HANDLE handle)
{
  if (!db || db->magic != MAGIC_BASE + MAGIC_AC_DB)
    messcrash ("ac_aql_table received a null or invalid db handle") ;
  
  return ac_newAql (db, initial_keyset, myquery, handle, orderBy, error_messages) ;
} /* ac_aql_table */
#endif
/********************************************************************/
/********************************************************************/
/* performs a BQL query, please free the hanlde afterwards
 * The results table and error messages and private data  are allocated on the handle
 */
#include "bql.h"
AC_TABLE ac_bql_table (AC_DB db, const char *query, AC_KEYSET initial_keyset, const char* orderBy, const char **error_messages, AC_HANDLE handle)
{
  AC_HANDLE h = ac_new_handle () ;
  BQL *bql ;
  char *command = 0 ;
  AC_TABLE results = 0 ;
  KEYSET resultKs = keySetHandleCreate (h) ;

  if (!db || db->magic != MAGIC_BASE + MAGIC_AC_DB)
    messcrash ("ac_aql_table received a null or invalid db handle") ;
  
  if ( ! query)
    messcrash ("ac_bql_table received a NULL query") ;
  if ( ! *query)
    messcrash ("ac_bql_table received an empty query") ;
   
  bql = bqlDbCreate (db, FALSE, handle) ; /* (BOOL debug, AC_HANDLE h) */
  command = halloc (strlen (query) + 100 + (orderBy ? strlen(orderBy) + 12 : 0), 0) ;
  sprintf (command, "%s %s %s", query, orderBy ? "order_by" : "", orderBy ? orderBy : "") ;
  if (bqlParse (bql, command, FALSE) &&
      bqlRun (bql,  initial_keyset ? initial_keyset->ks : 0, resultKs)
      )
    results = bqlResults (bql) ;
 
  if (error_messages)
    *error_messages = bqlError (bql) ;
  ac_free (h) ;
  return results ;
} /* ac_bql_table */

/********************************************************************/

AC_TABLE ac_obj_bql_table (AC_OBJ aobj, const char *query, const char* orderBy, const char **error_messages, AC_HANDLE handle)
{
  AC_KEYSET aks = 0 ;
  if (aobj == 0)
    messcrash ("Null obj in call to ac_obj_bql_table") ;
  aks = ac_new_keyset (aobj->db, handle) ;
  ac_keyset_add (aks, aobj) ;
  return ac_bql_table (aobj->db, query, aks, orderBy, error_messages, handle) ;
} /* ac_obj_bql_table */

/********************************************************************/
#if 0
const char * ac_keyname (AC_DB db, KEY key)
{
  return name (key) ;
} /* ac_db_keyname */
#endif
/********************************************************************/
/********************************************************************/
/*
 * PRIVATE
 */

static AC_TABLE ac_db_stuffed_table (AC_DB db, TABLE *tt, AC_HANDLE handle )
{
  int rows = tableMax (tt) ;
  int cols = tt->ncol ;
  AC_TABLE tbl = ac_db_empty_table (db, rows, cols, handle) ;
  char type ;
  const char *types = tableTypes (tt) ;
  int i, j, xi ;
  float xf = 0 ;
  mytime_t d ;
  const char *ccp ;
  BOOL b ;
  KEY key ;

  /* make room */
  ac_table_insert_type (tbl, rows -1, cols -1, 0, ac_type_empty) ;
  for (j = 0 ; j < cols ; j++)
    {
      type = types[j] ;
      if (!type || type == 'a')
	continue ; 
      tableCheckEmpty (tt, j) ;
      for (i = 0 ; i < tabMax(tt, j) ; i++)
	{
	  if (tabEmpty (tt, i, j))
	    continue ;
 /* legal types and corresponding characters are:
     k KEY
     g tag
     b bool
     i Int
     f Float
     s Text (string)
     t DateType  // d means DISK, kernel reserved
     a flagset   // flag set name is column name, minus "Flag" prefix if present
*/
	  switch ((int) type) 
	   {
	     case 'k':
	       key = tabKey (tt, i, j) ;
	       ac_table_insert_type (tbl, i, j, &key, ac_type_key) ;
	       break ;
	     case 'g':
	       key = tabTag (tt, i, j) ;
	       ac_table_insert_type (tbl, i, j, &key, ac_type_tag) ;
	       break ;
	     case 'b':
	       b = tabBool (tt, i, j) ;
	       ac_table_insert_type (tbl, i, j,&b, ac_type_bool) ;
	       break ;
	     case 'i':
	       xi = tabInt (tt, i, j) ;
	       ac_table_insert_type (tbl, i, j, &xi, ac_type_int) ;
	       break ;
	     case 'f':
	       xf = tabFloat (tt, i, j) ;
	       ac_table_insert_type (tbl, i, j, &xf, ac_type_float) ;
	       break ;
	     case 's':
	       ccp = tabString (tt, i, j) ;
	       ac_table_insert_type (tbl, i, j, (void**) (&ccp), ac_type_text) ;
	       break ;
	     case 't':
	       d = tabDate (tt, i, j) ;
	       ac_table_insert_type (tbl, i, j, &d, ac_type_date) ;
	       break ;
	   }
	}
    }
  return tbl ;
}  /* ac_stuffed_table */

/********************************************************************/

/* BUG not implemented */
/*
 * Table maker related functions
 *
 * Three constructors:
 *	tablemaker query
 *	aql query
 *	extract data from object
 */
/*
enum ac_tablemaker_where
{ 
ac_tablemaker_db, 
ac_tablemaker_file, 
ac_tablemaker_server_file, 
ac_tablemaker_text 
};
*/
/*
 * AC_TABLEmaker_db => query is name of tablemaker query in database
 * AC_TABLEmaker_file => query is file name known by server
 *		Not a particularly secure feature
 * AC_TABLEmaker_text => query is the actual text of a tablemaker
 *	query that you read into your source code.
 */

AC_TABLE ac_tablemaker_table (AC_DB db, const char *queryNam, AC_KEYSET iks,
	 enum ac_tablemaker_where where,
	 const char *parameters , const char* orderBy, const char **error_messages, AC_HANDLE handle)
{
  char *error_message = 0 ;
  AC_TABLE tbl = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  Stack defStack = 0 ;
  KEY tableDefKey = 0 ;
  SPREAD spread  = 0 ;
  FILE *f = 0 ;
  KEYSET kA = 0, ksTmp = 0 ;
  
  if (iks)
    aceCommandSwitchActiveSet (db->look, iks->ks, 0) ;
  
  if (!queryNam || !strlen(queryNam))
    {
      error_message = hprintf (handle, "// the definition passed to  ac_tablemaker_table is empty\n" ) ;
      goto abort ;
    }
  spread = spreadCreate () ;
  switch (where)
    {
    case ac_tablemaker_file:
      f = filopen (queryNam, 0, "r") ;
      if (!f)
	{
	  error_message = hprintf (handle, "// Sorry, could not find the table definitions file %s\n" , queryNam) ;
	  goto abort ;
	}
      if (! spreadDoReadDefinitions (spread, 0, f, 0, parameters, FALSE)) /* will close f */
	{
	  error_message = hprintf (handle, "// Sorry, there is a syntax error in the table definitions file %s, please test it in graphic acedb mode\n" , queryNam) ;
	  goto abort ;
	  }
      break ;
    case ac_tablemaker_text:
      defStack = stackHandleCreate (strlen(queryNam) + 24, h) ;
      pushText (defStack, queryNam) ;
      if (! spreadDoReadDefinitions (spread, 0, 0, defStack, parameters, FALSE)) /* will close f */
	{
	  error_message = 
	    hprintf (handle, "// Sorry, there is a syntax error in table definition passed to ac_tablemaker_table, please test it in graphic acedb mode\n" ) ;
	  goto abort ;
	}
	break ;
	
    case ac_tablemaker_db:
      if (! lexword2key (queryNam, &tableDefKey, _VTable)) 
	{
	  error_message = 
	    hprintf (handle, "// Sorry, currently is no table definition registered on the server under the name %s\n" , queryNam) ;
	  goto abort ;
	}
      if (! spreadDoReadDefinitions (spread, tableDefKey, 0, 0, parameters, FALSE)) /* will close f */
	{
	  error_message = 
	    hprintf (handle, "// Sorry, there is a syntax error in table definition registered on the server under the name %s, please test it in graphic acedb mode\n" , queryNam) ;
	  goto abort ;
	}
      break ;
    default:
      goto abort ;
    }

  if (iks)
    {	      
      spread->precompute = FALSE ; /* cannot store precomputed 
				    * table, if active keyset is 
				    * used for master column */
      spread->isActiveKeySet = TRUE ;
      ksTmp = spreadFilterKeySet (spread, iks->ks) ;
      kA = keySetAlphaHeap(ksTmp, keySetMax(ksTmp)) ;
      keySetDestroy (ksTmp) ;
    }
  else
    kA = spreadGetKeyset (spread) ;

  spreadRecomputeKeySet (spread, kA, 0, 0, keySetMax (kA)) ;
  /* now we must convert spread->table into AceC */ 
  {
    TABLE *tt = spreadToTable (spread, handle) ;
    if (tt)
      {
	if (orderBy &&
	    !tableSort (tt, orderBy))
	  error_message = 
	    hprintf (handle, "// Sorry, the sort parameter %s is incorrect , i should look like:  \"-2+4+1\" for " 
		     "first priority column 2 reversed, "
		     "then column 4" , orderBy) ;
	
      }
    if (tt) 
      tbl = ac_db_stuffed_table (db, tt, handle) ;
    ac_free (tt) ;
  }
  spreadDestroy (spread) ;
  if (0 && *error_messages) *error_messages = error_message ;

 abort:
  ac_free (h) ;

  return tbl ;
}
/*
 * perform tablemaker query:
 *	db is the database
 * 	query is the query
 * 	where describes how to interpret the value of query
 *	parameters are parameters to substitute in the query
 */

/********************************************************************/
/********************************************************************/
/********************************************************************/


/*
 * PUBLIC ac_keyset_table
 *
 * Create a 1 column table containing all objects in the keyset.
 */

AC_TABLE ac_keyset_table (AC_KEYSET aks, int min, int max, int fillhint, AC_HANDLE handle )
{
  AC_TABLE tbl;
  int mx = 0 ;

  if (aks && keySetExists (aks->ks))
    mx = keySetMax (aks->ks) ;
  if (max < 1)
    max = mx ? mx : 1 ;
  tbl = ac_db_empty_table (aks->db, max, 1, handle);

  if (mx)
    {
      KEYSET ksAlpha = keySetAlphaHeap (aks->ks, mx) ;
      KEY *kp ;
      int i ;
      
      for (i = 0, kp = arrp (ksAlpha, 0, KEY) ; i < mx ; i++, kp++)
	ac_table_insert_type (tbl, i, 0, kp, ac_type_key) ;
      
      keySetDestroy (ksAlpha) ;
    }
  return tbl;
}



/* x1 x2 in biolog coords [1,max] */
char *ac_zone_dna (AC_OBJ aobj, int x1, int x2, AC_HANDLE h)
{
  Stack s = 0 ;
  Array dna = dnaGet (aobj->key) ;
  char* result = 0;

  if (!dna) return 0 ;
  if (!x1 && !x2) { x1 = 1 ; x2 = arrayMax (dna) ; }
  if (x1 > x2)
    {
      reverseComplement (dna) ;
      x1 = arrayMax(dna) - x1 + 1 ;
      x2 = arrayMax(dna) - x2 + 1 ;
    }
  if (dna && x1 >= 1 && x2 >= 1 && x1 <= arrayMax(dna) && x2 <= arrayMax(dna) && x1 != x2)
    {
      dnaDecodeArray (dna) ;
      array (dna, x2, char) = 0 ;
      result = strnew (arrp (dna, x1 - 1, char), h) ;
    }
  arrayDestroy (dna) ;
  stackDestroy (s) ;
  return result ;
}

char *ac_obj_dna (AC_OBJ aobj, AC_HANDLE h)
{
  return ac_zone_dna (aobj, 0, 0, h) ;
}

char *ac_obj_peptide (AC_OBJ aobj, AC_HANDLE h)
{
  Array pep = peptideGet (aobj->key) ;
  int nn = pep ? arrayMax (pep) : 0 ;
  char* result = 0;

  if (nn)
    {
      pepDecodeArray (pep) ;
      result =  strnew (arrp (pep, 0, char), h) ;
    }
  arrayDestroy (pep) ;
  return result ;
}

char *ac_longtext (AC_OBJ aobj, AC_HANDLE h)
{
  Stack s = aobj ? stackGet (aobj->key) : 0 ;
  char *cp, *cq = 0, *cq0 = 0 ;

  if (s)
    {
      cq = cq0 = halloc (stackMark (s) + 10, h) ;
      stackCursor(s,0) ;
      while ((cp = stackNextText(s)))
	 {
	   sprintf(cq, "%s\n", cp) ; /* \n replaces one or more zero, so we do not go over size */
	   cq += strlen (cq) ;
	 }
 
      cp = cq0 - 1 ;
      while ((cp = strchr (cp + 1, '\n'))) *cp= ' ' ;
    }
  return cq0 ;
}

unsigned char *ac_a_data (AC_OBJ aobj, int *len, AC_HANDLE h) 
{
  /*
   * INCOMPLETE unless you add code in dump.c !!
   * nothing will get exported except for DNA Peptide
   *
   * other type A object
   *
   * If the object is of type A :
   *  - The returned value is the an array of bytes, allocated on handle
   *  - If len is not NULL, it is filled in with the number of bytes.
   *  - There is a \0 after the returned string (not included in the length).
   * Otherwise:
   *  - returns NULL and *len is undefined.
   *
   * The content of of the returned array is not specified, in
   * practice it is "whatever show -C does" and usually contains
   * formating characters that should be cleansed in the client.
   *	
   * Properly speaking, this should not be considered as a public function
   * but as a startup point if you want to handle a new type of A data.
   * A data types are declared on the server side in whooks/quovadis.c
   * and require registering dedicated code. In the same way one
   * should register here new AceC code to handle the specifics of that
   * class.
   */
  unsigned char *buf = 0 ;
  char *command ;

  if (aobj && pickType (aobj->key) == 'A')
    {
      command = hprintf (0, "query find %s IS \"%s\"", className (aobj->key), name(aobj->key)) ;
      buf = ac_command (aobj->db, command, len, h) ;
      messfree (buf) ;
      messfree (command) ;
      buf = ac_command_binary (aobj->db, "show -C", len, h, TRUE) ;
    }
  return buf ;
}

BOOL ac_parse (AC_DB db,
	       const char *text,
	       const char **error_messages,
	       AC_KEYSET *nks,
	       AC_HANDLE h)
{
  AC_KEYSET ks;
  BOOL success = TRUE;

  if (!db || db->magic != MAGIC_BASE + MAGIC_AC_DB)
    messcrash ("ac_parse received a null or invalid db handle") ;
  if (! sessionGainWriteAccess())
    {
      if (error_messages)
	*error_messages = strnew ("You do not have write access", 0) ;
      return FALSE ;
    }
  
  ks = ac_new_keyset (db, h);
  success = parseAceCBuffer (text, db->command_stack, 0, ks->ks) ;
  
  if (strstr( stackText(db->command_stack, 0), "error") )
    success = FALSE;
  
  if (error_messages)
    *error_messages = strnew( stackText(db->command_stack, 0), h);
  
  if (nks)
    *nks = ks;
  else
    ac_free(ks);
  
  return success;
}

