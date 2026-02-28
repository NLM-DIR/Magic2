#ifndef vTXT_H_DEF
#define  vTXT_H_DEF

#include "regular.h"
#include "array.h"
#include "utils.h"
/* It is crucial that a struct which is handled/alloc never
   starts on an other handled memory
   It is mandatory to start on a basic type
*/

typedef struct vTXT_struct {  BOOL markUp ; BigStack s;} *vTXT ;

vTXT vtxtCreate (void) ;
vTXT vtxtHandleCreate (AC_HANDLE h) ;
#define vtxtDestroy(_vtxt) messfree (_vtxt)

BOOL  vtxtClear (vTXT blkp) ; /* clears the content */
char *vtxtPtr (vTXT blkp) ; /* the content */
long int vtxtMark (vTXT blkp) ; /* Writable position, to be used in vtxtAt */
char *vtxtAt (vTXT blkp, long int pos) ; /* the content strating at pos, obtained from vtxtMark */
long int vtxtLen (vTXT blkp) ; /* the current length */
long int vtxtPrint (vTXT s, const char *txt) ; /* unformatted, uninterpreted, appends to what is already there */
void vtxtPercent (vTXT s, float z) ;/* format a float with optimal number of decimals */
long int vtxtPrintf (vTXT s, const char * formatDescription, ...) ; /* here an above, return an int valid for vtxtAt */
char *vtxtPrintWrapped (vTXT s, char *text, int lineLength) ; 
int vtxtReplaceString (vTXT vtxt, char *a, char *b) ; /* changes a into b, returns the number of replaced strings */
int vtxtRemoveBetween (char *txt, char *begin, char *end) ; /* removes all occurences of begin...end, returns the number of removed strings */

vTXT vtxtGetURL (const char *urlToGet, int timeOut) ; /* get the content of the url page, in libtcpfn.a */
vTXT vtxtHandleGetURL (const char *url, int timeOut, AC_HANDLE h) ;

char *vtextUpperCase (char *blkp) ; /* acts on a normal char buffer, uppers all letters */
char *vtextLowerCase (char *blkp) ; /* acts on a normal char buffer, lowers all letters */
void vtextCollapseSpaces (char *buf) ; /* change \n\t into space, collapse multiple spaces into single one */
void vtextReplaceSymbol (char *buf, char a, char b) ; /* acts on a normal char buffer */

void vtxtCleanHtml (vTXT vtxt) ;

/* text formatting */
BOOL vtxtMarkup (vTXT blkp) ; /* from now on, use html markup in following commands */
BOOL vtxtComma (vTXT blkp) ;
BOOL vtxtDot (vTXT blkp) ;
BOOL vtxtBreak (vTXT blkp) ;
BOOL vtxtEmptyLine (vTXT blkp, int count) ;
BOOL vtxtHr (vTXT blkp, int above, int below) ;
BOOL vtxtBold (vTXT blkp, char *text) ; 
BOOL vtxtItalic (vTXT blkp, char *text) ; 
BOOL vtxtSequence (vTXT blkp, char *text) ; 
BOOL vtxtPeptide (vTXT blkp, char *text, BOOL addStop) ;

vTXT vtxtTokenSort (vTXT txt, char *splitOn, AC_HANDLE h) ;

/*************************************************************/
/* Functions to analyse an xml document possibly downloaded from the web */

/* returns points to the first occurence of  '<tag>' in xml */
char *xmlGetTag (const char *xml, const char *tag) ;

/* extracts the data inside <tag>...data...</tag>
 * data is copied and allocated on handle h
 */
char *xmlGetTagContent (char **xmlp, const char *tag, char *maxPos, AC_HANDLE h) ;

/* usage: cp = xml ; while (cq = xmlGetNextTagContent (&cp, tag, h)) print (cq) ; */
/* extracts iterativelly the content of all occurences of tag in xml */
char *xmlGetNextTagContent (char **xmlpp, const char *tag, AC_HANDLE h) ;

/* gets as a vTXT the xml part of the URL page
 * i.e. the data in between <dd> and </dd>
 * in addition it recreates the <> symbols masked inside the URL page
 */
vTXT xmlGetDDContent (char *xml, AC_HANDLE h) ;


#endif
