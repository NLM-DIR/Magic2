/*  File: filsubs.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 *                   cross platform file system routines
 *              
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  8 16:33 1998 (fw)
 * * Dec  8 10:20 1998 (fw): new function filAge to determine time since
 *              last modification of file
 * * Oct 22 16:17 1998 (edgrif): Replace unsafe strtok with strstr.
 * * Oct 15 11:47 1998 (fw): include messSysErrorText in some messges
 * * Sep 30 09:37 1998 (edgrif): Replaced my strdup with acedb strnew.
 * * Sep  9 14:07 1998 (edgrif): Add filGetFileName routine that will 
 *               return the filename given a pathname 
 *              (NOT the same as the UNIX basename).
 * * DON'T KNOW WHO DID THE BELOW..assume Richard Bruskiewich (edgrif)
 *	-	fix root path detection for default drives (in WIN32)
 * * Oct  8 23:34 1996 (rd)
 *              filDirectory() returns a sorted Array of character 
 *              strings of the names of files, with specified ending 
 *              and spec's, listed in a given directory "dirName";
 *              If !dirName or directory is inaccessible, 
 *              the function returns 0
 * * Jun  6 17:58 1996 (rd)
 * * Mar 24 02:42 1995 (mieg)
 * * Feb 13 16:11 1993 (rd): allow "" endName, and call getwd if !*dname
 * * Sep 14 15:57 1992 (mieg): sorted alphabetically
 * * Sep  4 13:10 1992 (mieg): fixed NULL used improperly when 0 is meant
 * * Jul 20 09:35 1992 (aochi): Add directory names to query file chooser
 * * Jan 11 01:59 1992 (mieg): If file has no ending i suppress the point
 * * Nov 29 19:15 1991 (mieg): If file had no ending, we were losing the
                               last character in dirDraw()
 * Created: Fri Nov 29 19:15:34 1991 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: filsubs.c,v 1.27 2020/05/30 16:50:30 mieg Exp $	 */

#include <stdio.h>
#include "regular.h"
#include "mytime.h"
#include "call.h"		/* for callScript (to mail stuff) */

/********************************************************************/

#include "mydirent.h"

#include <sys/file.h>
#define HOME_DIR_ENVP "HOME"

#define ABSOLUTE_PATH(path) *path == SUBDIR_DELIMITER


/********************************************************************/

static Stack dirPath = 0 ;

void filAddDir (const char *s)	/* add to dirPath */
{
  char *home ;

  if (!dirPath)
    {
      dirPath = stackCreate (128) ;
    }
  /* if the user directory is specified */
  if (*s == '~' &&
		(home = getenv (HOME_DIR_ENVP))) /* substitute */
    {
      pushText (dirPath, home) ;
      catText (dirPath, ++s) ;
    }
  else
    pushText (dirPath, s) ;

  catText (dirPath, SUBDIR_DELIMITER_STR) ;

  return;
} /* filAddDir */

/*********************************************/

void filAddPath (const char *cp0)
{
  char *cp1, *cp, *cq ;
  
  cp1 = cq = cp = strnew (cp0, 0) ;

  while (TRUE)
    { 
      while (*cq && *cq != PATH_DELIMITER)
	++cq ;
      if (*cq == PATH_DELIMITER)
	{
	  *cq = 0 ;
	  filAddDir (cp) ;
	  cp = ++cq ;
	}
      else
	{ 
	  filAddDir (cp) ;
	  break ;
	}
    }

  messfree (cp1) ;
  return;
} /* filAddPath */


/*****************************************************************************/
/* This function returns the filename part of a given path,                  */
/*                                                                           */
/*   Given   "/some/load/of/directories/filename"  returns  "filename"       */
/*                                                                           */
/* The function returns NULL for the following errors:                       */
/*                                                                           */
/* 1) supplying a NULL ptr as the path                                       */
/* 2) supplying "" as the path                                               */
/* 3) supplying a path that ends in "/"                                      */
/*                                                                           */
/* NOTE, this function is _NOT_ the same as the UNIX basename command or the */
/* XPG4_UNIX basename() function which do different things.                  */
/*                                                                           */
/* The function makes a copy of the supplied path on which to work, this     */
/* copy is thrown away each time the function is called.                     */
/*                                                                           */
/*****************************************************************************/

char *filGetFileName (const char *path, AC_HANDLE h)
{
  char *path_copy = NULL ;
  const char *path_delim = SUBDIR_DELIMITER_STR ;
  char *result2, *result = NULL, *tmp ;
    
  if (path != NULL)
    {
      if (strcmp((path + strlen(path) - 1), path_delim) != 0) /* Last char = "/" ?? */
	{
	  if (path_copy != NULL) messfree(path_copy) ;
	  
	  path_copy = strnew(path, 0) ;
	  
	  tmp = path_copy ;
	  while (tmp != NULL)
	    {
	      result = tmp ;
	      
	      tmp = strstr(tmp, path_delim) ;
	      if (tmp != NULL) tmp++ ;
	    }
	}
    }
  result2 = strnew (result, h) ;
  messfree (path_copy) ;
  
  return result2 ;
} /* filGetFileName */


/*****************************************************************************/
/* This function returns the file-extension part of a given path/filename,   */
/*                                                                           */
/*   Given   "/some/load/of/directories/filename.ext"  returns  "ext"        */
/*                                                                           */
/* The function returns NULL for the following errors:                       */
/*                                                                           */
/* 1) supplying a NULL ptr as the path                                       */
/* 2) supplying a path with no filename                                      */
/*                                                                           */
/* The function returns "" for a filename that has no extension              */
/*                                                                           */
/* The function makes a copy of the supplied path on which to work, this     */
/* copy is thrown away each time the function is called.                     */
/*                                                                           */
/*****************************************************************************/
char *filGetExtension(const char *path)
{
  static char *path_copy = NULL ;
  char *extension = NULL, *cp ;
    
  if (path == NULL)
    return NULL;

  if (strlen(path) == 0)
    return NULL;

  if (path_copy != NULL) messfree(path_copy) ;
  path_copy =  (char *) messalloc ((strlen(path)+1) * sizeof(char));
  strcpy (path_copy, path);

  cp = path_copy + (strlen(path_copy) - 1);
  while (cp > path_copy &&
	 *cp != SUBDIR_DELIMITER &&
	 *cp != '.')
    --cp;

  extension = cp+1;
    
  return(extension) ;
} /* filGetExtension */


/**********************************************************************/
/* This function takes a directory name and does the following:
   1. Returns the name if it is "complete" 
      (an absolute path on a given platform)
   2. On WIN32 platforms, for onto rooted paths	lacking a 
      drive specification, returns the directory name prefixed with 
      the default drive letter 
   3. Otherwise, assumes that the directory name resides within the 
      current working directory and thus, returns it prefixes the
      directory name with the working directory path */
/**********************************************************************/
char *filGetFullPath(char *dir, AC_HANDLE handle)
{
  char *pwd ;
  char dirbuf[MAXPATHLEN] ;
  char *ret = 0 ; /* ret = 0  signals error that the path was not found */

  /* Return dir if absolute path already */
  if (ABSOLUTE_PATH(dir))
    { 
      ret =  strnew (dir, handle);
    }


  /* else if I can, then prefix "dir" with working directory path... */
  else if ((pwd = getcwd (dirbuf, MAXPATHLEN)))
    { 
      ret =  (char *) halloc (strlen(pwd) + strlen(dir) + 2, handle) ;

      strcpy (ret, pwd) ;
      strcat (ret, SUBDIR_DELIMITER_STR) ;
      strcat (ret, dir) ;
    }
  
  return ret ;    
} /* filGetFullPath */

/*******************************/

static BOOL filCheck (char *name, const char *spec)
	/* allow 'd' as second value of spec for a directory */
{
  char *cp ;
  BOOL result ;
  struct stat status ;

  if (!spec) /* so filName returns full file name (for error messages) */
    return TRUE ;
				/* directory check */
  if (spec[1] == 'd'  &&
      (stat (name, &status) || !(status.st_mode & S_IFDIR)))
    return FALSE ;

  switch (*spec)
    {
    case 'r':
      return !(access (name, R_OK)) ? TRUE : FALSE  ;
    case 'w':
    case 'a':
      if (!access (name, W_OK))	/* requires file exists */
	return TRUE ;
				/* test directory writable */
      cp = name + strlen (name) ;
      while (cp > name)
	if (*--cp == SUBDIR_DELIMITER) break ;
      if (cp == name)
	return !(access (".", W_OK)) ? TRUE : FALSE  ;
      else
	{ *cp = 0 ;
	  result = !(access (name, W_OK)) ? TRUE : FALSE  ;
	  *cp = SUBDIR_DELIMITER ;
	  return result ;
	}
    case 'x':
      return !(access (name, X_OK)) ? TRUE : FALSE  ;
    default:
      messcrash ("Unknown spec %s passed to filName", spec) ;
    }
  return FALSE ;
}

/************************************************/

static char *filDoName (const char *name, const char *ending, const char *spec, BOOL strict)
{
  static Stack part = 0, full = 0 ;
  char *dir, *result ;

  if (!name)
    messcrash ("filName received a null name") ;

  /*
  * we create these two buffer areas the first time we come through here
  */
  if (!part)
    { 
      part = stackCreate (128) ;
      full = stackCreate (MAXPATHLEN) ;
    /* 
    * add cwd as default to search.  This use of getcwd() is ok because
    * if you look above, you see that the Stack full was created with
    * MAXPATHLEN bytes in it.
    */
      /* filAddDir (".") ; 2025_04_22 now done in freeinit() */
      /* getcwd (stackText (full, 0), MAXPATHLEN )) ; */
    }

  /*
  * now empty out the workspace and put the desired file name in it
  */
  stackClear (part) ;

  catText (part, name) ;

  /*
  * append an extension if they asked for one.
  */
  if (ending && *ending)
    {
      if (*ending != '.' && *ending != '/') catText (part, ".") ;
      catText (part, ending) ;
    }
	/* NB filName is reentrant in the sense that it can be called 
	   on strings it generates, because they first get copied into 
	   part, and then the new name is constructed in full.
	*/
  /*
  * if they gave us a fully qualified file name, we don't need to search
  * the path - just use the name they gave us.
  */
  if (ABSOLUTE_PATH(name))
    {
      stackClear (full) ;
      catText (full, stackText (part, 0)) ;
      result = stackText (full, 0) ;

      if (filCheck (result, spec))
	return result ;
      else
	return 0 ;
    }
  
  /*
  * walk down the directory path, looking for the specified file
  * in each directory
  */
  if (! dirPath && ! strict)
    filAddDir (".") ; 
  stackCursor (dirPath, 0) ;
  dir = stackNextText (dirPath) ;
  if (! dir && ! strict)
    filAddDir (".") ; 
  stackCursor (dirPath, 0) ;
  while ((dir = stackNextText (dirPath)))
    { 
      stackClear (full) ;
      catText (full, dir) ;
      catText (full, stackText (part, 0)) ;
      result = stackText (full, 0) ;
      if (filCheck (result, spec))
	return result ;
      if (strict && strcmp (dir, "./"))
	break ;
    }
  return 0 ;
} /* filDoName */

/************************************************************/

char *filName (const char *name, const char *ending, const char *spec)
{ return filDoName(name, ending, spec, FALSE) ; }

/************************************************************/

BOOL filCheckName(const char *name, const char *ending, const char *spec)
{
  char *cp = filName (name, ending, spec) ;
  BOOL ok = cp ? TRUE : FALSE ;
  
  return ok ;
}

/************************************************************/

char *filStrictName (const char *name, const char *ending, const char *spec)
{ return filDoName(name, ending, spec, TRUE) ; }

/************************************************************/

BOOL filremove (const char *name, const char *ending) 
				/* TRUE if file is deleted. -HJC*/
{
  char *s = filName (name, ending, "r") ;
  if (s)
    return unlink(s) ? FALSE : TRUE ;
  else
    return FALSE ;
} /* filremove */

/************************************************************/

FILE *filopen (const char *name, const char *ending, const char *spec)
{
  char *s = filName (name, ending, spec) ;
  FILE *result = 0 ;
   
  if (!s)
    {
      if (spec[0] == 'r')
	messerror ("Failed to open for reading: %s (%s)",
		   filName (name, ending,0),
		   messSysErrorText()) ;
      else if (spec[0] == 'w')
	messerror ("Failed to open for writing: %s (%s)",
		   filName (name, ending,0),
		   messSysErrorText()) ;
      else if (spec[0] == 'a')
	messerror ("Failed to open for appending: %s (%s)",
		   filName (name, ending,0),
		   messSysErrorText()) ;
      else
	messcrash ("filopen() received invalid filespec %s",
		   spec ? spec : "(null)");
    }
  else 
    {
      int pass = 6 ; /* up to 1 minute wait */
      int delay = 1 ;
      result = 0 ;
      while (! result && pass--)
	{
	  result = fopen (s, spec) ;
	   if (! result)
	     sleep (delay) ;
	   delay *= 2 ;
	}

      if (!result)
	messerror ("Failed to open %s (%s)",
		   s, messSysErrorText()) ;
    }
  return result ;
} /* filopen */

/********************* temporary file stuff *****************/

/*
* create or open a temporary file
*
* for open, it is much like filopen()
*
* for create, it generates a temp file name, creates and opens
* it, and returns the file name created.  
*
*/
static Associator tmpFiles = 0 ;

FILE *filtmpopenOLD (char **nameptr, const char *spec)
{
#if defined(SUN) || defined(SOLARIS)
  /*
  * not sure why we avoid /tmp on SUN, but this is backward compatible
  */
  char *prefix="/var/tmp";
#else
  char *prefix="/tmp";
#endif
  static char *name_prototype = 0;
  char *newname;
  int fd;

  if (!nameptr)
    messcrash ("filtmpopen requires a non-null nameptr") ;

  /*
  * if opening the file for reading, we don't need to do anything
  * special.
  */
  if (!strcmp (spec, "r"))
    return filopen (*nameptr, 0, spec) ;

  /*
  * retaining a degree of compatibility with tempnam - use TMPDIR
  * if set, otherwise use our own default temp directory.
  */
  if (! name_prototype )
    {
    /*
    * pick either TMPDIR or our own prefix, and generate prototype
    * filename for mkstemp to use
    */
    char *s = getenv("TMPDIR");
    if (!s) s = prefix;
    name_prototype= (char *) messalloc(strlen(s)+1+5+8+1);
    strcpy(name_prototype,s);
    strcat(name_prototype,"/ACEDB");
    strcat(name_prototype,"XXXXXXXX");
    }

  /*
  * the newly allocated filename is something you can free() when you
  * are done with it.  (This is compatible with tempnam() which was in
  * the old version.)  You copy it because mkstemp() mangles the string
  * that you pass in.
  */
  newname = strnew (name_prototype, 0);

  /*
  * generate the new name AND open it.  mkstemp is supposed to be
  * immune to races.
  */
  fd = mkstemp(newname);
  if ( fd < 0 )
    { 
      messerror ("failed to create temporary file %s", newname) ;
      messfree (newname) ;
      return 0 ;
    }

  /*
  * make the name available to calling program
  */

  *nameptr = newname ;

  /*
  * store the file name in an associator - we can remove temp files
  * later by iterating over the associator.  We probably don't just
  * open the file and then unlink it immediately because that doesn't 
  * work on Windows.
  */

  if (!tmpFiles)
    tmpFiles = assCreate () ;
  assInsert (tmpFiles, *nameptr, *nameptr) ;

  /*
  * we already have the file open, but we want stdio.
  */
  return fdopen (fd, spec) ;
} /* filtmpopen */

/********************* temporary file stuff *****************/

static FILE *filTmpOpenWithSuffix (char **nameptr, char *suffix, const char *spec)
{
  char *realName;
  int fd = -1 ;

#ifdef __CYGWIN__
  char *tmpenv = getenv("TEMP");
  char tmppath[PATH_MAX];
  
  if (tmpenv)
    cygwin_conv_to_full_posix_path(tmpenv, tmppath);
  else
    cygwin_conv_to_full_posix_path("c:\\Temp", tmppath);
#endif

  if (!nameptr)
    messcrash ("filtmpopen requires a non-null nameptr") ;


  
#if defined(SUN) || defined(SOLARIS)
  realName = strnew (messprintf ("%s/ACEDB.XXXXXXXX", "/var/tmp"), 0) ;
#elif defined(__CYGWIN__)
  realName = strnew (messprintf ("%s/ACEDB.XXXXXXXX", tmppath), 0) ;
#else
realName = strnew (messprintf ("%s/ACEDB.XXXXXXXX", "/tmp"), 0) ;
#endif
  fd = mkstemp (realName) ;
  if (fd < 0)
    { 
      messerror ("failed to create temporary file (%s)",
		 messSysErrorText()) ;
      return 0 ;
    }

    
  if (!tmpFiles)
    tmpFiles = assCreate () ;
  assInsert (tmpFiles, realName, realName) ;
  *nameptr = realName ;
  return fdopen (fd, spec) ;
} /* filTmpopenWithSuffix */

FILE *filtmpopen (char **nameptr, const char *spec)
{
  return filTmpOpenWithSuffix(nameptr, 0, spec);
} /* filtmpopen */


/************************************************************/

BOOL filtmpremove (const char *name)	/* delete and free()  */
{
  BOOL result = filremove (name, 0) ;
  assRemove (tmpFiles, name) ;
  messfree(name);
  return result ;
}

/************************************************************/

void filtmpcleanup (void)
{ 
  const void *name = 0 ;
 
  if (tmpFiles)
    while (assNext (tmpFiles, &name, 0))
      { 
	char *name1 = (char *) name;
	filremove (name1, 0) ;
	messfree (name1) ; /* Since messfree zeros it. */
      }
}

/************* filqueryopen() ****************/

static QueryOpenRoutine queryOpenFunc = 0 ;

QueryOpenRoutine filQueryOpenRegister (QueryOpenRoutine neuf)
{ QueryOpenRoutine old = queryOpenFunc ; queryOpenFunc = neuf ; return old ; }

FILE *filqueryopen (char *dname, char *fname, char *end
				  , const char *spec, const char *title)
{
  Stack s ;
  FILE*	fil = 0 ;
  int i ;
				/* use registered routine if available */
  if (queryOpenFunc)
    return (*queryOpenFunc)(dname, fname, end, spec, title) ;

  /* otherwise do here and use messprompt() */
  s = stackCreate(50);

  if (dname && *dname)
    { pushText(s, dname) ; catText(s,"/") ; }
  if (fname)
    catText(s,fname) ; 
  if (end && *end)
    { catText(s,".") ; catText(s,end) ; }

 lao:
  if (!messPrompt("File name please", stackText(s,0), "w")) 
    { stackDestroy(s) ;
      return 0 ;
    }
  i = stackMark(s) ;
  pushText(s, freepath()) ;	/* freepath needed by WIN32 */
  if (spec[0] == 'w' && 
      (fil = fopen (stackText(s,i), "r")))
    { if ( fil != stdin && fil != stdout && fil != stderr)
	fclose (fil) ; 
      fil = 0 ;
      if (messQuery (messprintf ("Overwrite %s?",
				 stackText(s,i))))
	{ 
	  if ((fil = fopen (stackText(s,i), spec)))
	    goto bravo ;
	  else
	    messout ("Sorry, can't open file %s for writing",
		     stackText (s,i)) ;
	}
      goto lao ;
    }
  else if (!(fil = fopen (stackText(s,i), spec))) 
    messout ("Sorry, can't open file %s",
	     stackText(s,i)) ;
bravo:
  stackDestroy(s) ;
  return fil ;
} /* filqueryopen */

/*********************************************/

Associator mailFile = 0, mailAddress = 0 ;

void filclose (FILE *fil)
{
  void *address ;
  void *filename ;

  if (!fil || fil == stdin || fil == stdout || fil == stderr)
    return ;
  if (mailFile && assFind (mailFile, fil, &filename))
    { if (assFind (mailAddress, fil, &address))
	callScript ("mail", messprintf ("%s < %s", address, filename)) ;
      else
	messerror ("Have lost the address for mailfile %s", filename) ;
      assRemove (mailFile, fil) ;
      assRemove (mailAddress, fil) ;
      unlink ((char *) filename) ;
      messfree (filename) ;
    }
  fclose (fil) ;
} /* filclose */

/***********************************/

FILE *filmail (const char *address)	/* requires filclose() */
{
  char *filename = "fil_tmp_open_failed" ;
  FILE *fil ;

  if (!mailFile)
    { mailFile = assCreate () ;
      mailAddress = assCreate () ;
    }
  if (!(fil = filtmpopen (&filename, "w")))
    { messout ("failed to open temporary mail file %s", filename) ;
      return 0 ;
    }
  assInsert (mailFile, fil, filename) ;
  assInsert (mailAddress, fil, address) ;
  return fil ;
} /* filmail */

/******************* directory stuff *************************/

static int dirOrder(const void *a, const void *b)
{
  const char *cp1 = *(const char **)a, *cp2 = *(const char**)b;
  return strcmp(cp1, cp2) ;
} /* dirOrder */

/* returns an Array of strings representing the filename in the
   given directory according to the spec. "r" will list all files,
   and "rd" will list all directories.
   The behaviour of the "w" spec is undefined.
   The array has to be destroyed using filDirectoryDestroy,
   because the memory of the strings needs to be reclaimed as well. */

Array filDirectoryCreate (char *dirName,
					char *ending, 
					char *spec)
{
  Array a ;

  DIR	*dirp ;
  char	*dName, entryPathName[MAXPATHLEN], *leaf ;
  int	dLen, endLen ;
  MYDIRENT *dent ;

  if (!dirName || !(dirp = opendir (dirName)))
    return 0 ;

  if (ending)
    endLen = strlen (ending) ;
  else
    endLen = 0 ;

  strcpy (entryPathName, dirName) ;
  strcat (entryPathName, "/") ;
  leaf = entryPathName + strlen(dirName) + 1 ;

  a = arrayCreate (16, char*) ;
  while ((dent = readdir (dirp)))           
    { dName = dent->d_name ;
      dLen = strlen (dName) ;
      if (endLen && (dLen <= endLen ||
		     dName[dLen-endLen-1] != '.'  ||
		     strcmp (&dName[dLen-endLen],ending)))
	continue ;

      strcpy (leaf, dName) ;
      if (!filCheck (entryPathName, spec))
	continue ;

      if (ending && dName[dLen - endLen - 1] == '.') /* remove ending */
	dName[dLen - endLen - 1] = 0 ;

      /* the memory of these strings is freed my 
	 the messfree()'s in filDirectoryDestroy() */
      array(a,arrayMax(a),char*) = strnew(dName, 0) ;
    }
  
  closedir (dirp) ;
  
  /************* reorder ********************/
    
  arraySort(a, dirOrder) ;
  return a ;

} /* filDirectoryCreate */

/*************************************************************/

void filDirectoryDestroy (Array filDirArray)
{
  int i;
  char *cp;

  for (i = 0; i < arrayMax(filDirArray); ++i)
    {
      cp = arr(filDirArray, i, char*);

      messfree (cp);
    }
  arrayDestroy (filDirArray);

  return;
} /* filDirectoryDestroy */

/************************************************************/
/* determines the age of a file, according to its last modification date.

   returns TRUE if the age could determined and the int-pointers
   (if non-NULL will be filled with the numbers).

   returns FALSE if the file doesn't exist, is not readable,
   or the age could not be detrmined. */
/************************************************************/
BOOL filAge (const char *name, const char *end,
	     int *diffYears, int *diffMonths, int *diffDays,
	     int *diffHours, int *diffMins, int *diffSecs)
{
  struct stat status;
  mytime_t time_now, time_modified;
  char time_modified_str[25];
	  
  /* get the last-modification time of the file,
     we parse the time into two acedb-style time structs
     in order to compare them using the timediff functions */
  
  if (!(filName (name, end, "r")))
    return FALSE;

  if (stat (filName (name, end, "r"), &status) == -1)
    return FALSE;

  {
    time_t t = status.st_mtime;
    struct tm *ts;

    /* convert the time_t time into a tm-struct time */
    ts = localtime(&t);		

    /* get a string with that time in it */
    strftime (time_modified_str, 25, "%Y-%m-%d_%H:%M:%S", ts) ;

    time_now =      timeNow();
    time_modified = timeParse(time_modified_str);

    if (diffYears)
      timeDiffYears (time_modified, time_now, diffYears);

    if (diffMonths)
      timeDiffMonths (time_modified, time_now, diffMonths);

    if (diffDays)
      timeDiffDays (time_modified, time_now, diffDays);

    if (diffHours)
      timeDiffHours (time_modified, time_now, diffHours);

    if (diffMins)
      timeDiffMins (time_modified, time_now, diffMins);

    if (diffSecs)
      timeDiffSecs (time_modified, time_now, diffSecs);
  }
  return TRUE;
} /* filAge */


/*************************************************************/


/* This function copies a file (surprisingly there is no standard Posix      */
/* function to do this). The file is created with read/write permissions     */
/* for the user. Note that it is permissible for the new file to exist, in   */
/* which case its contents will be OVERWRITTEN. This is to allow the caller  */
/* to create a file with a unique name, close it and then supply it as an    */
/* argument to this call (e.g. caller may use aceTmpCreate() to create the   */
/* file).                                                                    */
/*                                                                           */
/* The function returns FALSE for the following errors:                      */
/*                                                                           */
/* 1) supplying a NULL ptr or empty string for either file name.             */
/* 2) if the file to be copied does not exist/is not readable.               */
/* 3) if the file to be created is not writetable.                           */
/* 4) if the copy fails for some other reason (e.g. read/write failed)       */
/*                                                                           */
/* This code is adapted from ffcopy.c from "POSIX Programmers Guide"         */
/* by Donald Levine, publ. by O'Reilly.                                      */
/*                                                                           */
/* It would be handy to use access() to check whether we can read a file     */
/* but this only uses the real UID, not the effective UID. stat() returns    */
/* info. about the file mode but we would need to get hold of all sorts of   */
/* other information (effective UID, group membership) to use it.            */
/*                                                                           */
BOOL filCopyFile(char *curr_name, char *new_name)
{
  BOOL status = TRUE ;
  struct stat file_stat ;
  size_t buf_size ;
  size_t bytes_left ;
  enum {BUF_MAX_BYTES = 4194304} ;			    /* Buffer can be up to this size. */
  char *buffer = NULL ;
  int curr_file = -1, new_file = -1 ;

  /* File names supplied ?                                                   */
  if (!curr_name || !(*curr_name) || !new_name || !(*new_name))
    status = FALSE ;

  /* Make sure file to be copied exists and can be opened for reading.       */
  if (status)
    {
      if (stat(curr_name, &file_stat) != 0)
       status = FALSE ;
      else
	{
	  bytes_left = buf_size = file_stat.st_size ;
	  if (buf_size > BUF_MAX_BYTES)
	    buf_size = BUF_MAX_BYTES ;
	  
	  if ((curr_file = open(curr_name, O_RDONLY, 0)) == -1)
	    status = FALSE ;
	}
    }

  /* Make sure file to be written can be opened for writing (O_TRUNC means   */
  /* existing contents of the file will be overwritten).                     */
  if (status)
    {
      if ((new_file = open(new_name, O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR)) == -1)
	status = FALSE ;
    }

  
  /* allocate a buffer for the whole file, or the maximum chunk if the file  */
  /* is bigger.                                                              */
  if (status)
    {
      if((buffer = (char *)messalloc(buf_size)) == NULL)
	status = FALSE ;
    }
  

  /* Copy the file.                                                          */
  if (status)
    {
      while (bytes_left > 0 && status == TRUE)
	{
	  if (read(curr_file, buffer, buf_size) != buf_size)
	    {
	      status = FALSE ;
	    }

	  if (status == TRUE)
	    {
	      if(write(new_file, buffer, buf_size) != buf_size)
		  status = FALSE ;
	    }

	  bytes_left -= buf_size;
	  if (bytes_left < buf_size) buf_size = bytes_left ;
	}

    }

  /* Clear up buffer and files.                                              */
  if (buffer != NULL)
    messfree(buffer) ;

  if (curr_file > -1)
    {
      if (close(curr_file) != 0)
	messcrash("close() of file being copied failed.") ;
    }
  if (new_file > -1)
    {
      if (close(new_file) != 0)
	messcrash("close() of file being copied failed.") ;
    }

  return status ;
}



/************************************************************/
/********************** end of file *************************/

