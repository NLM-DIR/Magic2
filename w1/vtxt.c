/*  File: vtxt.c
 *  Author: Jean Thierry-Mieg (mieg@kaa.crbm.cnrs-mop.fr)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
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
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Contains variable length text interface
                to simplify the preparation of text documents
                The code is implemented over the acedb stack functions
 * Exported functions: see wh/vtxt.h
 * HISTORY:
 * Created: Thu May 8 2003: mieg
 *-------------------------------------------------------------------
 */

#include "ac.h"
#include "call.h"

/*********************************************************************/
/* Typical usage
 *
 * we use recursive calls vtxt to fill 
 * chapters, paragraphs, sentences
 * and in case of success, we print the title
 * and return the whole thing


f1 (vTXT v1)
{
  vTXT v2 = vtxtCreate () ;
  if (ptr = f2 (s2,..)) // return 0 or bigStackText (s2, 0)
    {
      vtxtPrintf (v1, "Title\n") ;
      vtxtPrint (v1, ptr) ; // unformated, so we do not interpret % and \ a second time
    }
  vtxtDestroy (v2) ; // allways needed
}
*/

/*********************************************************************/

static void vtxtFinalise (void *v)
{
  vTXT blkp = (vTXT) v ;

  if (blkp)
    {
      if (bigStackExists (blkp->s))
	bigStackDestroy (blkp->s) ;
    }
  return ;
}

/*********************************************************************/

vTXT vtxtHandleCreate (AC_HANDLE h)
{
  vTXT blkp = (vTXT) halloc (sizeof (struct vTXT_struct), h) ;
  
  blkp->s = bigStackHandleCreate (1024, h) ;
  if (! h) blockSetFinalise (blkp, vtxtFinalise) ;
  
  return blkp ;
} /* vtxtCreate */

/*********************************************************************/

vTXT vtxtCreate (void)
{
  return vtxtHandleCreate (0) ;
} /* vtxtCreate */

/*********************************************************************/

BOOL vtxtClear (vTXT s)
{
  bigStackClear (s->s) ;
  return  bigStackExists(s->s)  ? TRUE : FALSE ;
} /* vtxtPtr */

/*********************************************************************/

char *vtxtPtr (vTXT s)
{
  return bigStackExists(s->s) && bigStackMark (s->s)? bigStackText (s->s, 0) : 0 ;
} /* vtxtPtr */

/*********************************************************************/

long int vtxtLen (vTXT s)
{
  return bigStackExists(s->s) ? bigStackMark (s->s) : 0 ;
} /* vtxtLen */

/*********************************************************************/

BOOL vtxtMarkup (vTXT s)
{
  s->markUp = TRUE ;
  return TRUE ;
} /* vtxtMarkup */

/*********************************************************************/
/* Writable position, to be used in vtxtAt */
long int vtxtMark (vTXT s)
{
  long int nn = bigStackMark (s->s) ;
  return nn ;
} /* vtxtMark */

/*********************************************************************/
/* the content strating at pos*/
char *vtxtAt (vTXT s, long int pos) 
{
  char *ptr = vtxtPtr (s) ;
  return (ptr && vtxtLen (s) > pos) ? ptr + pos : 0 ;
} /* vtxtAt */

/*********************************************************************/
/* change a into b in the whole vTXT 
 * possibly realloc (therefore invalidates previous vtxtPtr return values
 *
 * returns the number of replaced strings
 */

int vtxtReplaceString (vTXT vtxt, char *a, char *b)
{
  char *cp, *cq, *buf = 0 ;
  int nn = 0 ;

  cp = vtxt ? vtxtPtr (vtxt) : 0 ;
  cq = cp && a ? strstr (cp, a) : 0 ;

  if (!cq)
    return 0 ;
  buf = cp = strnew (cp, 0) ;
  vtxtClear (vtxt) ;
  while ((cq = strstr (cp, a)))
    {
      nn++ ;
      *cq = 0 ;
      vtxtPrint (vtxt, cp) ;
      if (b && *b) vtxtPrint (vtxt, b) ;
      cp = cq + strlen (a) ;
    }
  vtxtPrint (vtxt, cp) ;
  ac_free (buf) ;

  return nn ;
} /* vtxtReplaceString */

/*********************************************************************/
/* removes all occurences of begin...end,
 *
 * returns the number of removed strings
 */

int vtxtRemoveBetween (char *txt, char *begin, char *end)
{
  char *cp, *cq, *cr ;
  int nn = 0 ;

  cp = txt && *txt && begin && *begin && end && *end 
    ? strstr (txt, begin) : 0 ;

  while (cp && (cq = strstr (cp, end)))
    {
      nn++ ;
      cr = cp ;
      while ((*cp++ = *cr++)) ;
      cp = strstr (cp, end) ;
    }
  return nn ;
} /* vtxtRemoveBetween */

/*********************************************************************/

char *vtextUpperCase (char *blkp) /* acts on a normal char buffer */
{
  char *cp = blkp ;
  while (*cp) { *cp = ace_upper(*cp) ; cp++; }

  return blkp ;
} /* vtextUpperCase */

/*********************************************************************/

char *vtextLowerCase (char *blkp) /* acts on a normal char buffer */
{
  char *cp = blkp ;
  while (*cp) { *cp = ace_lower(*cp) ; cp++; }
  
  return blkp ;
} /* vtextLowerCase */

/*********************************************************************/
/* equivalent to vstrFindReplaceSymbols(sBuf,sBuf,0,"\t \n"," ",0,1,1); */
/* change \n\t into space, collapse multiple spaces into single one
 * removes starting and trailing spaces 
 */
void vtextCollapseSpaces (char *buf) /* acts on a normal char buffer */
{
  char *cp, *cq ;
  BOOL isWhite = TRUE ; /* to remove starting spaces */
  
  cp = cq = buf ; cp-- ;
  if (buf)
    while (*++cp)
      switch (*cp)
	{
	case ' ': /* case '\n': case '\t': */
	  if (!isWhite) *cq++ = ' ' ;
	  isWhite = TRUE ;
	  break ;
	default:
	  if (cp != cq) *cq = *cp ;
	  cq++ ; 
	  isWhite = FALSE ;
	  break ;
	}
  *cq = 0 ;
  while (--cq > buf && *cq == ' ') *cq = 0 ; /* remove trailing space */
} /* vtextCollapseSpaces  */

/*********************************************************************/
/* change a into b in the whole of buf */
void vtextReplaceSymbol (char *buf, char a, char b) /* acts on a normal char buffer */
{
  char *cp, *cq ;
  
  cp = cq = buf ; cp-- ;
  if (cp)
    while (*++cp)
      {
	if (*cp == a)
	  { if (b) *cq++ = b ; }
	else
	  {
	    if (cp != cq) *cq = *cp ;
	    cq++ ;
	  } 
      }
  *cq = 0 ;
} /* vtextReplaceSymbol */

/*********************************************************************/
/* removes html markup, i.e. everything inside <> */
void vtxtCleanHtml (vTXT vtxt)
{
  char *cp, *cq, *cp0 ;
  int inHtml = 0 ;
  int inQuote = 0 ;
  int i ;
  
  cp = cq = cp0 = vtxtPtr (vtxt) ; cp-- ;
  if (cq && *cq)
    while (*++cp)
      {
	if (!inQuote && *cp == '<')
	  inHtml++ ;
	if (*cp == '\"')
	  {
	    for (i = 1 ; cp - i >= cp0 && *(cp-i) == '\\' ; i++) ;
	    if (i%2)  /* distinguish "  \"  \\" \\\" */
	      inQuote = 1 - inQuote ;
	  }
	if (!inHtml)
	  *cq++ = *cp ; /* copy one char */
	if (inHtml && !inQuote && *cp == '>')
	  inHtml-- ;
      }
  while (cp >= cq) *cp-- = 0 ;
} /* vtextCleanHtml */

/*********************************************************************/
/* Unformatted print, appends to what is already there */

long int vtxtPrint (vTXT s, const char *txt)
{
  long int nn, len ;
  char *cp ;
  
  nn = bigStackMark (s->s) ;
  if (txt && (len = strlen (txt)))
    {
      /* programmer's note
       * this function differs from catText or catBinary
       * if the bigStack is empty and we vtxtPrint (s,"aa")
       * as a result we have aa0  (zero terminated)
       * and s->ptr, used for the next vtxtPrint or vtxtPrintf
       * points on that zero
       * whereas in catBinary, we point on the 4th char, after the zero
       * but the next catBinary or catText will move back one char
       * so the 2 types of calls should not be mixed
       * although a VTXT is implemented as a BigStack
       */
      bigStackExtend (s->s, len + 1) ;
      cp = bigStackText (s->s, nn) ;
      
      memcpy (cp, txt, len) ;
      *(cp + len) = 0 ;
      s->s->ptr += len ;
    }

  return nn ;
} /* vtxtPrint */

/*************************************************************************************/
/* export a percentage z, always with at least 4 significant digits, so we can see 99.999993 */
void vtxtPercent (vTXT s, float z)
{
  if (z == 0)
    vtxtPrintf (s, "0") ;
  else if (z == 100)
    vtxtPrintf (s, "100") ;
  else  if (z < 0)
    vtxtPrintf (s, "%.2f", z) ;
  else if (z > 100)
    vtxtPrintf (s, "%.2f", z) ;
  else if (z == 50)
    vtxtPrintf (s, "50") ;
  else if (z < 50 && z > 40)
    {
      float z1  = 50 - z ;
      int n = 0 ;
      char *f ;

      while (z1 < 100) { n++ ; z1 *= 10 ; }
      if (n < 2) n = 2 ;
      f = messprintf ("%%.%df", n) ;
      vtxtPrintf (s, f, z) ;
    }
  else if (z > 50 && z < 60)
    {
      float z1  = z - 50 ;
      int n = 0 ;
      char *f ;

      while (z1 < 100) { n++ ; z1 *= 10 ; }
      if (n < 2) n = 2 ;
      f = messprintf ("%%.%df", n) ;
      vtxtPrintf (s, f, z) ;
    }
  else if (z < 41)
    {
      float z1  = z ;
      int n = 0 ;
      char *f ;

      while (z1 < 100) { n++ ; z1 *= 10 ; } 
      if (n < 2) n = 2 ;
      f = messprintf ("%%.%df", n) ;
      vtxtPrintf (s, f, z) ; 
    }
  else 
    {
      float z1  = 100 - z ;
      int n = 0 ;
      char *f ;

      while (z1 < 100) { n++ ; z1 *= 10 ; }
      if (n < 2) n = 2 ;
      f = messprintf ("%%.%df", n) ;
      vtxtPrintf (s, f, z) ;
    }
} /* vtxtPercent */

/*********************************************************************/
/* General formatted printf, appends to what is already there */
long int vtxtPrintf (vTXT s, const char * format, ...)
{
  long int nn, len ;
  char *cp ;
  
  va_list ap;
  va_start(ap, format);
  len = utPrintfSizeOfArgList (format, ap) ;
  va_end (ap) ;

  bigStackExtend (s->s, len + 1) ;
  nn = bigStackMark (s->s) ;
  cp = bigStackText (s->s, nn) ;

  va_start(ap, format);
  len = vsprintf(cp, format, ap) ;
  va_end (ap) ;

  s->s->ptr += len;
  return nn ;
} /* vtxtPrintf */

char *vtxtPrintWrapped (vTXT s, char *text, int lineLength) 
{
  char cc, *cp, *cq ;
  int ii ;

  cp = text ;
  while (*cp)
    {
      ii = lineLength ; cq = cp ;
      while (ii-- && *cq) cq++ ;
      cc = *cq ; *cq = 0 ;
      vtxtPrintf (s, "%s%s", cp, cc ? "\n" : "") ;
      *cq = cc ;
      cp = cq ;
    }
  return bigStackText (s->s, 0) ;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

BOOL vtxtComma (vTXT blkp)
{
  char 
    *cp = vtxtPtr (blkp),
    *cq = cp ? cp + strlen (cp) - 1 : 0,
    *cr = cq ;

  if (cq) 
    while (cq > cp)
      switch (*cq)
	{
	case '.': case ';':
	  *cq = ',' ;
	  /* fall trhu */
	case ',': 
	  if (cq == cr) vtxtPrint (blkp, " ") ;
	  return  FALSE ;
	case ' ': cq-- ; break ;
	default: 
	  if (cq == cr) { vtxtPrint (blkp, ", ") ; return TRUE ; }
	  if (cq == cr - 1) { *(cq+1) = ',' ; vtxtPrint (blkp, " ") ; return TRUE ; }
	  if (cq <= cr - 1) { *(cq+1) = ',' ; return TRUE ; }
	}
  return FALSE ;
} /* vtxtComma */

/******************************************************************************/

BOOL vtxtDot (vTXT blkp)
{
  char 
    *cp = vtxtPtr (blkp),
    *cq = cp ? cp + strlen (cp) - 1 : 0,
    *cr = cq ;

  if (cq) 
    while (cq > cp)
      switch (*cq)
	{
	case '\n': return  FALSE ;
	case ',': case ';':
	  *cq = '.' ;
	  /* fall trhu */
	case '.': case ':': /* we may pass a line but we keep the : */
	  if (cq == cr) vtxtPrint (blkp, " ") ;
	  return TRUE ; 
	case ' ': cq-- ; break ;
	default: 
	  if (cq == cr) { vtxtPrint (blkp, ". ") ; return TRUE ; }
	  if (cq == cr - 1) { *(cq+1) = '.' ; vtxtPrint (blkp, " ") ; return TRUE ; }
	  if (cq <= cr - 1) { *(cq+1) = '.' ; return TRUE ; }
	}
  return FALSE ;
} /* vtxtDot */

/******************************************************************************/

BOOL vtxtBreak (vTXT blkp)
{
  if (vtxtDot (blkp))
    {
      if (blkp->markUp)
	vtxtPrint (blkp, "<br>\n") ;
      else
	vtxtPrint (blkp, "\n") ;
      return TRUE ;
    }
  return FALSE ;
} /* vtxtBreak */

/******************************************************************************/

BOOL vtxtEmptyLine (vTXT blkp, int count)
{
  int i ;
  vtxtBreak (blkp) ;
  
  for (i = 0 ; i < count ; i++)
    {
      if (blkp->markUp) vtxtPrint (blkp, "<br>\n") ;
      else vtxtPrint (blkp, "\n") ;
    }
  return TRUE ;
} /* vtxtEmptyLine */

/******************************************************************************/

BOOL vtxtHr (vTXT blkp, int above, int below)
{
  vtxtBreak (blkp) ;
  vtxtEmptyLine (blkp, above) ;

  if (blkp->markUp) 
    vtxtPrint (blkp, "<hr></hr>\n") ;
  else
    vtxtPrint (blkp, "\n") ;

   vtxtEmptyLine (blkp, below) ;
  return TRUE ;
} /* vtxtHr */

/******************************************************************************/

BOOL vtxtBold (vTXT blkp, char *text) 
{
  if (blkp->markUp) vtxtPrint (blkp, "<b>") ;
  vtxtPrint (blkp,  text) ;
  if (blkp->markUp) vtxtPrint (blkp, "</b>") ; 
  return TRUE ;
} /* vtxtBold */

/******************************************************************************/

BOOL vtxtItalic (vTXT blkp, char *text) 
{
  if (blkp->markUp) vtxtPrint (blkp, "<i>") ;
  vtxtPrint (blkp,  text) ;
  if (blkp->markUp) vtxtPrint (blkp, "</i>") ; 
  return TRUE ;
} /* vtxtItalicc */

/******************************************************************************/

BOOL vtxtSequence (vTXT blkp, char *text) 
{
  char cc, *cp, *cq ;
  int nn = strlen (text) ;
  int dx = 100 ;

  if (blkp->markUp)
    vtxtPrint (blkp, "<font face='courier'>\n<pre DNA_START>\n") ;
  cp = text ; 
  
  while (nn > dx)
    {
      cq = cp + dx ;
      cc = *cq ;
      *cq = 0 ;
      vtxtPrint (blkp,  cp) ;
      *cq = cc ;
      cp = cq ; nn -= dx ; 
      if (blkp->markUp)
	vtxtPrint (blkp, "<DNA>\n") ;
      else
	vtxtPrint (blkp, "\n") ;
    }
  if (nn)
    vtxtPrint (blkp,  cp) ;
  if (blkp->markUp) 
    vtxtPrint (blkp, "<DNA>\n</pre DNA_END>\n</font>\n<br>\n") ; 
  else
    vtxtPrint (blkp, "\n") ;
  return TRUE ;
} /* vtxtSequence */  

/******************************************************************************/

BOOL vtxtPeptide (vTXT blkp, char *text, BOOL addStop)
{
  BOOL ok ;

  if (addStop) 
    {
      vTXT buf = vtxtCreate () ;
      vtxtPrintf (buf, "%s*", text) ;
      ok  = vtxtSequence (blkp, vtxtPtr (buf)) ;
      vtxtDestroy (buf) ;
    }
  else
    ok = vtxtSequence (blkp, text) ;

  return ok ;
} /*vtxtPeptide  */

/******************************************************************************/
/******************************************************************************/

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  FILE BASIC FUNCTIONS                    _/
_/  copied from Vahan                       _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

#include <wh/vtxt_.h>

vFILE vtxtFileCreate (const char *fileName, int mode)  
{
  int hFile = creat (fileName, mode) ;
  return hFile == -1 ? 0 : hFile ;
} /* vtxtFileCreate */

/******************************************************************************/

vFILE vtxtFileOpen (const char *fileName, int flag, int mode)  
{
  int hFile ;
  flag |= O_BINARY ;
  hFile = open (fileName, flag, mode) ;
  return hFile == -1 ? 0 : hFile ;
} /* vtxtFileOpen */

/******************************************************************************/

int vtxtFileRead (vFILE fileHandle, char *Buffer, int BuffSize)  
{
  return read (fileHandle, Buffer, BuffSize) ;
}

/******************************************************************************/

char *vtxtFileReadUntil (vFILE hFile,int * pPos,char * lookFor, int minOfs, AC_HANDLE h)
{
  char * ptr = 0, bufr[vFILEBLOCKSIZE+1];
  int	len = 1, opos, pos = 0, ln = strlen (lookFor);
  vTXT blkp = vtxtHandleCreate (h) ;
  
  vtxtFileSetPos (hFile,*pPos);
  
  while (!ptr)
    {
      if (!(len = vtxtFileRead (hFile, bufr, sizeof(bufr)-1)))
	break;
      bufr[len]=0;
      
      opos = vtxtPrintf (blkp, "%s", bufr) - ln;
      if (opos < minOfs)
	opos = minOfs;
      
      if ((ptr = strstr (vtxtPtr(blkp) + opos, lookFor)) )
	{ pos = ptr - vtxtPtr(blkp) ; break ; }
      pos += len ;
      minOfs = 0 ;
    }
  
  if (ptr) *ptr = 0 ;
  (*pPos) += pos ;
  return vtxtPtr(blkp) ;
} /* vtxtFileReadUntil */

/******************************************************************************/

int vtxtFileWrite (vFILE fileHandle, const char *Buffer, int BuffSize)
{
  return write (fileHandle, Buffer, BuffSize) ;
}

/******************************************************************************/

void vtxtFileClose (vFILE fileHandle)
{
  if (fileHandle != 0 && fileHandle != -1)
    close (fileHandle) ;
}

/******************************************************************************/

int vtxtFileSetPos (vFILE fileHandle, myoff_t offset)
{
  return offset == ACEDB_MAXINT ?
    lseek (fileHandle, 0, SEEK_END) : 
    lseek (fileHandle, offset, SEEK_SET) ;
}

/******************************************************************************/

int vtxtFileGetPos (vFILE fileHandle)
{
  return lseek (fileHandle, 0, SEEK_CUR) ;
}

/******************************************************************************/

BOOL vtxtFileRemove (char *fileName)
{
  return filremove (fileName, 0) ;
}

/******************************************************************************/

int vtxtFileGetLength (vFILE fileHandle)
{
  int curPos, len ;
  
  curPos = lseek (fileHandle, 0, SEEK_CUR) ;
  len = lseek (fileHandle, 0, SEEK_END) ;
  lseek (fileHandle, curPos, SEEK_SET) ;
  return len ;
}

/******************************************************************************/

BOOL vtxtFileIsEndReached (vFILE fileHandle)
{
  return vtxtFileGetLength (fileHandle) == vtxtFileGetPos (fileHandle) ? TRUE : FALSE  ;
}

/******************************************************************************/
/******************************************************************************/
/*
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  _/                                          _/
  _/  FILENAME BASED FUNCTIONS                _/
  _/                                          _/
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

char *vtxtFileMakePath (const char *dnam,  const char *flnm, const char *ext, AC_HANDLE h)
{
  int len = strlen (dnam) ;
  char *destination = hprintf (h
			       ,"%s%s%s%s%s", dnam, 
			       (len == 0 || dnam[len-1] == vFILEDIRSEPARATOR[0]) ?  "" : vFILEDIRSEPARATOR, 
			       flnm,  (ext[0] == '.' || ext[0] == 0) ? "" : ".", ext
			       ) ;
  return destination ;
}

/******************************************************************************/
/* if file is readeable return it's size else return 0 ; */
int vtxtFileIsReadeable (const char *fileName)
{
  int len = 0 ;
  vFILE hFile = 0 ;
  
  if ((fileName) && 
      (hFile = vtxtFileOpen (fileName, O_RDONLY, vFILEDEFAULT)))
    len = vtxtFileGetLength (hFile) ;
  
  vtxtFileClose (hFile) ;
  return len ;
}

/******************************************************************************/
/* allocates buffer reads the file content and returns it's address */
char *vtxtFileGetContent (const char *fileName, int *lengthp, AC_HANDLE h)
{
  vFILE hFile = 0 ;
  int  fileLength ;
  char *contentBuf = 0 ;
  
  if (lengthp)
    *lengthp = 0 ;
  if ((fileName) &&
      (hFile = vtxtFileOpen (fileName, O_RDONLY, vFILEDEFAULT)) &&
      (fileLength = vtxtFileGetLength (hFile)) &&
      (contentBuf = halloc (fileLength+2, h)) && /* two for double zero at the end */
      (fileLength = vtxtFileRead (hFile, contentBuf, fileLength))
      )
    {
      contentBuf[fileLength] = contentBuf[fileLength+1] = 0 ;
      if (lengthp) *lengthp = fileLength ;
    }
  else
    ac_free (contentBuf) ;
  
  vtxtFileClose (hFile) ;
  return contentBuf ;
}


/******************************************************************************/
/* writes the buffer into a file */
int vtxtFileSetContent (const char *fileName, const char *contentBuf)
{
  int len = 0 ;
  vFILE hFile = 0 ;
  
  if (fileName &&
      contentBuf &&
      (hFile = vtxtFileOpen (fileName, O_TRUNC|O_WRONLY|O_CREAT, vFILEDEFAULT))) 
    len = vtxtFileWrite (hFile, contentBuf, strlen (contentBuf)) ;
  
  vtxtFileClose (hFile) ;
  return len ;
}

/******************************************************************************/
/* writes the buffer into the end of the file */
char *vtxtFileAppend (const char *fileName, char *contentBuf)
{
  int len = 0 ;
  vFILE hFile = 0 ;
  BOOL nonStdout = FALSE ;
  
  if ((fileName) &&
	(nonStdout = strcmp (fileName, "stdout")) && 
	(contentBuf) &&
	(hFile = vtxtFileOpen (fileName, O_WRONLY|O_CREAT, vFILEDEFAULT)) && 
	(vtxtFileSetPos (hFile, ACEDB_MAXINT) != -1))
    len = vtxtFileWrite (hFile, contentBuf, strlen (contentBuf)) ;
  
  if (!nonStdout && contentBuf)
    len = printf ("%s", contentBuf) ;
  else 
    vtxtFileClose (hFile) ;
  
  return len ? contentBuf : 0 ;
}

/******************************************************************************/
/******************************************************************************/
/*
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  _/                                          _/
  _/  FILE EXECUTION FUNCTION                 _/
  _/    obsolete: prefer callPipe             _/
  _/                                          _/
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

/* this method is less efficient 'callPipe' defined in call.h */
/* execute the file with given stdin from inputBuf and put the result into stdoutFile  */
char *vtxtFileExecution (const char *tmplt, const char *commandLine, const char *inputBuf, const char *stdoutFile, AC_HANDLE h)
{
  char *res ;
  vTXT cmd = vtxtCreate () ;
  char fileSTDIN [MAXPATHLEN], fileSTDOUT[MAXPATHLEN] ;
  
  /* prepare stdin file and clean up stdout file */
  if  (inputBuf) 
    {
      sprintf (fileSTDIN, "%s.stdin.tmp", tmplt) ;
      vtxtFileRemove (fileSTDIN) ; 
      vtxtFileSetContent (fileSTDIN, inputBuf) ;  
    }
  
  /* if nothing is passed use default stdout name */
  if (!stdoutFile)
    sprintf (fileSTDOUT, "%s.stdout.tmp", tmplt) ; 
  else 
    strcpy (fileSTDOUT, stdoutFile) ; 
  
  vtxtFileRemove (fileSTDOUT) ;
  
  /* execute the command */
  vtxtPrintf (cmd, " ( %s ) > %s %s %s\n", commandLine, 
	   fileSTDOUT, inputBuf ? "<" : "",  inputBuf ? fileSTDIN :"") ;
  
  if  (callSystem (vtxtPtr (cmd)) != -1)
    {
      if  (!stdoutFile)
	res = vtxtFileGetContent (fileSTDOUT, 0, h) ; 
      else 
	res = strnew (stdoutFile, h) ;
    } 
  else
    res = 0 ;
  
  if  (!stdoutFile)
    vtxtFileRemove (fileSTDOUT) ;
  if  (inputBuf)
    vtxtFileRemove (fileSTDIN) ;  /* clean up stdin,  stdout files */
  
  vtxtDestroy (cmd) ;
  return res ;
}

/******************************************************************************/
/******************************************************************************/
/*
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  _/                                          _/
  _/  DIRECTORY FUNCTIONS                     _/
  _/                                          _/
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/


char *vtxtDirGetCurrent (AC_HANDLE h)
{
  char buf [MAXPATHLEN] ;
  getcwd (buf, MAXPATHLEN) ;
  return  strnew (buf, 0) ;
}

/******************************************************************************/

#ifdef JUNK

DICT *vtxtDirList (const char *dnam, AC_HANDLE h)
{   
  DICT *dict = dictHandleCreate (64, h) ;
  int ifiles = 0 ;
  char *dst ;
  
#ifdef  WIN32
  {
    struct _finddata_t c_file ;
    vFILE    hFile ; 
    int res ;
    char   filepath[MAXPATHLEN] ;
    
    vtxtFileMakePath (filepath, dnam, "*", "*") ;
    
    hFile = _findfirst (filepath,  &c_file ) ;
    if (hFile != -1){
      for (res = 0 ;res == 0 ;)
	{
	  if (c_file.name[0] != '.' ||  (c_file.name[1] != '.' && c_file.name[1] != 0))
	    dictAdd (dict, dt->d_name) ;
	  
	  res = _findnext ( hFile,  &c_file ) ;
	}
      _findclose ( hFile ) ;
    }
    ac_free (filepath) ;
  }
#else
  { 
    DIR * mdr ;
    struct dirent *dt ;
    
    mdr = opendir (dnam) ;
    if (mdr)
      while ((dt = readdir (mdr)))
	if (dt->d_name[0] != '.' ||  (dt->d_name[1] != '.' && dt->d_name[1] != 0))
	  dictAdd (dict, dt->d_name) ;
    closedir (mdr) ;
  }
#endif    
  
  if (! dictMax (dict))
    dictDestroy (dict) ;
  return dict ;
}

#endif

/********************************************************************/
/*************************************************************/
/*************************************************************/
/* extracts the data inside <tag>...data...</tag>
 * from the buffer *xmlp, if successful: moves *xmlp pass the tag
 *   beware and correctly include nested blocks with same tag
 * data is copied and allocated on handle h
 */
char *xmlGetTagContent (char **xmlp, const char *tag, char *maxPos, AC_HANDLE h)
{
  char buf1 [1000], buf2 [1000] ;
  char *cp = 0, *cq = 0, *cp1 = 0, *cp2 = 0, *cq1 = 0, *xml1 ;
  int n2, n3;

  if (tag && xmlp && strlen(tag) < 990)
    {
      /* <Gene-commentary_type value="phenotype">19</Gene-commentary_type> */
      sprintf (buf1, "<%s", tag) ;
      sprintf (buf2, "</%s>", tag) ;
      if ((cp = strstr (buf2, " ")))
	{ *cp++ = '>' ; *cp++ = 0 ; }      

      n2 = strlen (buf2) ;
      cp = strstr (*xmlp, buf1) ;
      if (!cp) return 0 ;
      if (maxPos && cp > maxPos)
	return 0 ;
      cp2 = strstr (cp, ">") ;
      if (!cp2) return 0 ;
      cq = strstr (cp2+1, buf2) ;
      if (!cq) /* the block does not exist */
	return 0 ;
      xml1 = cp2 + 1 ;
      cp1 = xmlGetTagContent (&xml1, tag, cq, h) ; /* look for nested repeat */
      if (cp1 && xml1 < cq) /* jump the inner block */
	{
	  cq1 = strstr (xml1, buf2) ;
	  if (!cq1 || cq1 > cq)   /* unbalanced  <tag>..<tag> .. </tag only once> */
	    return 0 ;
	  cq = cq1 ;
	}
      n3 = cq - cp2 - 1 ;
      cp1 = halloc (n3 + 1, h) ;
      memcpy (cp1, cp2 + 1, n3) ;
      cp1[n3] = 0 ;
      if (cp)
	(*xmlp) = cp2 + 1 + n2 + n3 ;
    }
  return cp1 ;
} /* xmlGetTagContent  */

/*************************************************************/
/* usage: cp = xml ; while (cq = xmlGetNextTagContent (&cp, tag, h)) print (cq) ; */
/* extracts iterativelly the content of all occurences of tag in xml */
char *xmlGetNextTagContent (char **xmlpp, const char *tag, AC_HANDLE h)
{
  char *cp ;
  
  cp = xmlGetTagContent (xmlpp, tag, 0, h) ;
  return cp ;
} /* xmlGetNextTagContent */

/*************************************************************/
/* gets as a vTXT the xml part of the URL page
 * i.e. the data in between <dd> and </dd>
 * in addition it recreates the <> symbols masked inside the URL page
 */
vTXT xmlGetDDContent (char *xml, AC_HANDLE h)
{
  /* we wish to say '*(const char)' but i do not know the correct syntax */
  char *cp = xmlGetTagContent (&xml, "dd", 0, h) ;
  vTXT txt = 0 ;

  if (cp)
    {
      txt = vtxtHandleCreate (h) ;
      vtxtPrint (txt, cp) ;
      vtxtReplaceString (txt, "&lt;", "<") ;
      vtxtReplaceString (txt, "&gt;", ">") ;
    }
  return txt ;
} /* xmlGetDDContent  */

/********************************************************************/
/********************************************************************/
/* split a txt on some characters, typically ",:;"
 * neteurn a new txt with the same tokens but in alphbetic order
 */
 
vTXT vtxtTokenSort (vTXT txt0, char *splitOn, AC_HANDLE h0)
{
  char *cp = txt0 ? vtxtPtr (txt0) : 0 ;
  if ( 0 && cp)
    {
      AC_HANDLE h = ac_new_handle () ;
      KEYSET ks = keySetHandleCreate (h) ;
      int i, k, nn = strlen (cp) ;
      char buf[nn+1] ;
      strncpy (buf, cp, nn) ; buf[nn] = 0 ;
      vTXT txt = vtxtHandleCreate (h0) ;
      char *cq ; /* , *cr ; */
      while (cp)
	{
	  for (i = k = 0, cq = cp ; *cq ; i++, cq++)
	    {
	      keySet (ks, k++) = i ;
	      buf[i] = *cq ;
	      if (strchr (splitOn, *cq))
		{
		  buf[i] = 0 ;
		  break ;
		  cp = cq + 1 ;
		}
	    }
	}
      
      return txt ;
      ac_free (h) ;
    }
  return txt0 ;  /* i did nothing  */
} /* vtxtTokenSort */

/********************************************************************/
