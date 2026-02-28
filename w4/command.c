/*  File: command.c
 *  Author: Danielle et jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1995
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@frmop11.bitnet
 *
 * Description: Command language for acedb
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 11 13:43 1999 (fw)
 * * Dec  3 14:41 1998 (edgrif): Change calls to new interface to aceversion.
 * * Oct  8 11:44 1998 (fw): changed calls to helpPrintf to
                       helpOn() call. The helpsystem now uses
		       a registered funcion for whatever display
		       funcion we want to get at the help.
		       The default is obviously a textual dump.
 * * Aug 14 12:26 1998 (fw): AQL: context var for active keyset @active
 * * Aug  5 17:18 1998 (fw): added AQL functionality, replaced by BQL in 2017
 * * Nov 18 19:27 1997 (rd)
 *		- Read .zip files Assumes that xyz.zip file (only) 
 *		  contains a xyz.ace file; minor menu changes for "p" option
 *		- Used freepath() (see filsubs.c) instead of freeword() for filename
 *		  parses in commands "write", "Keyset read","parse","table reads" etc. 
 * * Jun 11 15:40 (rbrusk): #include "acedb.h" for name(), className() etc.
 * * Apr 23 16:45 1996 (srk)
 * Created: Sun Nov  5 17:41:29 1995 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: command.c,v 1.152 2020/06/07 15:36:51 mieg Exp $ */


/********************* MAIN LOOP *****************************/
#define CHRONO

#include "ac.h"
#include "call.h"
#include "a.h"
#include "../w4/command7_.h"
#include "command.h"
#include "spread_.h"		/* YUCK!!! HACK HACK HACK */
#include "freeout.h"
#include "java.h"		/* for freejavaprotect() */
#include "session.h"
#include "dna.h"
#include "peptide.h"
#include "dump.h"
#include "lex.h"
#include "block.h" /*for blockavail()*/
#include "../whooks/sysclass.h"
#include "../whooks/systags.h"
#include "../whooks/tags.h"
#include "pick.h"
#include "disk.h"
#include "cache.h"
#include "biblio.h"
#include "dump.h"
#include "table.h"
#include "parse.h"
#include "bql.h"
#include "help.h"
#include "alignment.h"
#include "model.h"
#include "longtext.h"		/* for longGrep */
#include "chrono.h"		/* for longGrep */
#include <ctype.h>		/* for isdigit */
#include "whooks/classes.h"
#include "acembly.h"

/************************************************************/

void (*gifEntry)(KEYSET, int, BOOL) = 0 ;	/* global, reset in giface code */
void (*nqEntry)(void) = 0 ;	/* global, reset in giface code */


static void taceTest (KEYSET ks) ; /* at end of this file */

typedef struct { int choix; KEY key; char *text ;} CHOIX ;

/* choix = set correct bits
   1: universal
   2: non-server
   4: acembly
   8: server
   16: gif
*/

/************************************************************/

static CHOIX choixMenu [] = 
{ 
{ 1,  'q', "Quit : Exits the program. Any command abbreviates, $starts a subshell, untill empty line"},
{ 1,  '?', "? : List of commands"},
{ 1,  'h', "Help : for further help"},
{ 1,  'c', "Classes : Give a list of all the visible class names and how many objects they contain."},
{ 1,  'z', "Model [-jaq] model : Shows the model of a class or subtype, useful before Follow or Query commands"},
{ 1,  'n', "Find [-max n] class [name] : Creates a new list with all [at most n]  objects from class, optionally matching name"},
{ 1,  't', "Follow Tag : i.e. Tag Author gives the authors of the papers of the current list"}, 
{ 1,  'g', "Grep [-active] template : searches for *template* everywhere in the database except LongTexts"},
{ 1,  'G', "LongGrep template : searches for *template* also in LongTexts"},
{ 1,  'l', "List [-h | -a | -p | -j | -J | -C] [-b begin] [-c max] [-f filename] [[-t] template] : \n\tlists for [human | ace | perl | java | C] [in file] current names [matching template]"},
{ 1,  's', "Show [-h | -a | -p | -j | -C] [-T] [-b begin] [-c max] [-f filename]  [-C] [[-t] tag] : \n\tshows for [human | ace | perl | java | C],\n\t-T show time stamps, -b begin at object b, -c show max obj,\n \t-C ignore comments, the last parameter restricts to the [tag-branch of] objects of current list"},
{ 1,  'i', "Is [-max n]  [ < <= > >= ~ ] template : keeps in list only those [max n] objects whose name match the template"},
{ 1,  'r', "Remove template : removes from list those objects whose name match the template"},
{ 1,  'x', "Query query_string : performs a complex query - 'help query_syntax' for further info"},
{ 1,  'x', "Where query_string : same as Query, kept for backwards  compatibility, do NOT use"},
{ 1,  'T', "Table-Maker [-active] [-a | -h | -j | -p | -x | -C] [-b begin] [-c count] [-title] [-href www...] [-i bql_query_file] [-o output_file] [-bqlo file_name: export bql query] [-s name:save definition loaded via -f as an internal object] [-n name |  [-f] table_def_file | = command\ncommand] [params] : \n\t[Search active list] [ace,human,java,perl,C style] table_definition_file [parameters (noted %%d in the def]"},
{ 1,  'b', "Biblio [-abstracts] : shows the associated bibliography [including the abstracts]"},
{ 1,  'd', "Dna [-mismatch] [-p | -C | -B] [file] [-x1 u1 -x2 u2] [[-f] [-noClassName] outputfile]: Fasta [or perl or C style]\n// dump of related dna sequences  [allowing mismatches], \n//the values u1 and u2 of the -x1 -x2 options may be the word \"begin\" OR the word \"end\" OR a positive number between 1 and dna length,\n// if u1 > u2, the dna is reverse complemented"},
{ 1,  'v', "Peptide [file] : Fasta dump of related peptide sequences"},
{ 4,  'V', "View [-v view] [-active | -c class -n name] [-xml | -asn | -tc | -flash] [-preview [-mask int]] [-o outputfile] [-p params to end of line]: View configurated data"},
{ 4,  500, "WebQuery class name : Optimised web query creates an active keyset of genes related to the name"},
{ 4,  600, "Magic [-o file] view run : Export an html page using google graphic API"},
{ 0,  'L', "Alignment [file] : Fasta dump of related sequences"},
{ 0,  'a', "Array class :name format : formatted acedump of A class object"},
{ 1,  'K', "KeySet-Read [filename | = Locus a;b;Sequence c] : \n\tread keyset from file filename or from stdin or from command line"},
{ 1,  '1', "spush : push active keyset onto stack"},
{ 1,  '2', "spop  : remove top keyset from stack and make it the active keyset"},
{ 1,  '3', "swap : swap keyset on top of stack and active keyset"},
{ 1,  '4', "sand  : replace top of stack by intersection of top of stack and active keyset"},
{ 1,  '5', "sor  : replace top of stack by union of top of stack and active keyset"},
{ 1,  '6', "sxor  : replace top of stack by exclusive or of top of stack and active keyset"},
{ 1,  '7', "sminus  : removes members of the active keyset from the top of the stack "},
{ 1,  '9', "subset x0 nx : returns the subset of the alphabetic active keyset from x0>=1 of length nx"},
{ 1,  301, "kstore name  : store active keyset in internal tmp space as name  "},
{ 1,  302, "kget name  : copy named set as active keyset "},
{ 1,  303, "kclear name  : clear named set "},
         /* parsing in server/client is done via -ace_in command */
#if defined(WIN32)
{ 2,  'p', "Parse [-NOXREF] file : parses pkzip (.zip) compressed or native ace files or stdin,  until EOF"},
#else
{ 2,  'p', "Parse [-NOXREF] [file] : parse an ace or ace.Z or ace.gz file or stdin,  until ctrl-D or EOF"},
#endif
{ 2,  'P', "PParse [-NOXREF] [file] : as above, but dont stop on first error"},  
{ 8,  'p', "Parse : [-NOXREF] parse a file defined on the server, \n\trun aceclient -ace_in to parse a file or stdin defined on the client side"},
{ 8,  'P', "PParse [-NOXREF] [file] : as above, but dont stop on first error"},  
{ 8,  '_', "serverparse : parse files in server mode"},  
{ 1,  'w', "Write -f filename : acedump current list to file"},
{ 1,  'e', "Edit ace_file_command : Apply the command to every member of the activelist"},
{ 1,  'E', "EEdit ace_file_command : idem No Check"},
{ 1,  'O', "Comment [comment] : Export a comment to the log file"},
{ 8,  'U', "Shutdown : Closes the server (you need special priviledge)"},
{ 8,  'W', "Who : who is connected now (you need special priviledge)"},
{ 8,  'X', "Data_Version : data served now, as stated in wspec/server.wrm"},
{ 0,  'm', "More : If there are more objects to list"}, 
{ 0,  'M', "More : If there are more objects to show"}, 
{ 0,  'B', "More : If there are more lines in the table"}, 
{ 0,  'D', "More : If there are more dna or peptide to dump"},
{ 2,  'R', "Read-models [-index]: reads models,[force reindex of whole database] (admin only)"},
{ 1,  'k', "Kill : Kill and remove from acedb all objects in the current list"},
{ 2,  'Y', "tace [-f command_file] [parameters] : Execute a file of acedb commands"},
{ 4,  'A', "Acembly [-f command_file] [parameters] : Sequence assembly and automatic edition"},
{ 1,  'Z', "Status {on | off} : toggle memory statistics"},
{ 1,  'J', "Date : prints the date and time"},
{ 1,  'I', "Time_Stamps {on | off} : toggle time stamps creation"},
{ 1,  'H', "Chrono {start | stop | show} built in profiler of chrono aware routines"},
{ 5,  300, "Test : test subroutine, may vary"},
/* { 1,  450, "S  [-a | -j | -J | -C] [-o outfile]  : AQL: Acedb Query Language [silent | ace | java | Java | C  style] <aql query>"}, */
{ 1,  451, "Select [-o outfile] [-s] [-h | -a | -j | -J | -C]  [-b begin] [-c count] [-title] [-q] <query>: Acedb Query Language\n\tif the query ends with a ;  (silent mode) do not export, but modify the active list\n\toutfile: name of the output file, relative to $ACEDB\n\t ahjJC: different output formats adapted to different acedb clients\n\t-b 3 -c 12: data exploration tool, export 12 lines, starting approximately at line 3, the edges depend on the size of the objects\n\t-title export a title line prefixed with #\n\t-q old acedb syntax Find ...\n\t-v verbose (for debugging purpose)\n\tsize of the objects, useful only to explore the datacount: ext[output] [silent] {format:human,ace,java,Java,C style.\n"}, 
{ 1,  452, "AQL : alias of Select, same parameters\n"},
{ 1,  453, "BQL : Example:\t\tSelect a,p,j from a in ?Author, p in a->Paper, j in p->Journal where j == \"Nature\"\n\tThe full syntax is presented in wbql/bql.pdf\n\tThe operators are:\n\t\t ; TITLE order_by\n\t\tfrom select , where\n\t\tin := = \n\t\t|| OR ^^ XOR &&  AND ! NOT\n\t\tlike =~ ~ == != >= <= > < >~ <~\n\t\tISA\n\t\tmodulo + - * / ^\n\t\tDNA PEPTIDE PEP CODING DATEDIFF \n\t\tCOUNT MIN MAX SUM AVERAGE STDDEV\n\t\t.class .name .timestamp\n\t\t>> -> # => \n\t\t@ OBJECT\n\t\t() {} [] \n"},
{ 1,  454, "S : alias of Select, same parameters\n"},
{ 1,  'F', "New Class Format : i.e New Plate ya\\%dx.y creates ya8x.y if ya7x.y exists"},
{ 1,  'N', "Count : number of objects in active keyset"},
{ 1,  'C', "Clear : Clear current keyset"},
{ 1,  'S', "Save : Save current state otherwise the system saves every 600 seconds"},
{  2,  '<', "Read-only : Definitively drop present and future write access"},
{ 1,  'u', "Undo : returns the current list to its previous state"},
{ 1,  '/', "// : comment, do nothing"},
{ 0,  '-', "Depad [seqname] : depads an assembly, or currently padded assembly"},
{ 0,  '*', "Pad seqname : pads an assembly"},
{ 16, '#', "GIF : access to gif drawing submenu"},
{ 8, 400, "wspec : export wspec directory as a tar file"},
{ 1, 401, "layout : export the layout file"},
{ 1, 402, "system_list : exports in C format names of all tags and classes"},
{ 2,  'j', "Dump [-s] [-T] [-C] [directory] : dump database, -s splits by class, -T gives timestamps, -C ignore comments"},
{ 0,0,0 }
 } ;

/* note that a tace specific help may be preferable

static FREEOPT helpmenu[] =	// no \n inside FREEOPT text (else rethink)
{ 6, "helpmenu",
  1, "Tacedb : copyright etc", 
  2, "Query_syntax for the find command",
  3, "Useful_sequences : SL1 etc",
  4, "Clone_types used in the physical map (vector etc)",
  5, "DNA_and_amino_acids_nomenclature : codes for bases/amino acids and translation table",
  6, "index"
} ;
*/

/*************************************************************/
/*************************************************************/

static Array choisirMenu (int choix)
{ Array aMenu = arrayCreate (50, FREEOPT) ;
  int i ;
  CHOIX *fc = &choixMenu[0] ;
  FREEOPT *ff ;

  i = 1 ;
  while (fc->key)
    { if (fc->choix & choix)
	{ ff = arrayp (aMenu, i++, FREEOPT) ;
	  ff->key = fc->key ;
	  ff->text = fc->text ;
	}
      fc++ ;
    }
  ff = arrayp (aMenu, 0, FREEOPT) ;
  ff->key = i - 1 ; ff->text = "acedb" ;
  return aMenu ;
}

/*************************************************************/
/******************* Special routines ************************/
/*************************************************************/

static int getOnOff (void)
{ 
  char *cp ;
  
  freenext () ;
  if ((cp = freeword ()))
    {
      if (!strcasecmp (cp, "off"))
	return (int) 'f' ;
      else if (!strcasecmp (cp, "on"))
	return (int) 't' ;
      else
	return (int) '?' ;
    }
  return 0 ;
}

/*************************************************************/

static void myText
(char *cp, int x, int y)
{ 
  char ww[1000] ;
  int i = 990 ;
  char *cq = ww ;

  if (x < 0) x = 0 ;
  while (x--) *cq++ = ' ' ;
  *cq++ = '/' ;   *cq++ = '/' ; *cq++ = ' ' ;
  while (i-- && *cp) *cq++ = *cp++ ;
  *cq++ = '\n' ; *cq++ = 0 ;
  freeOut (ww) ;
}

/**************************/

void tStatus (void)
{
 unsigned long keynum,keynum2, vocnum,vocspace,hspace,tspace ;
 mysize_t aUsed, aMade, aAlloc, aReal ;
 unsigned long int
   dused, dfree, plus, minus, ndread, ndwrite ;
 int
   n,
   bsused, bsalloc, btused, btalloc, 
   nmessalloc, messalloctot,
   bused, bpinned, bfree, bmodif,bsmemory,
   indexT, indexKb, indexV,
   cacheused, cachel, cachek, cachemod, 
   nc1r,nc1w,ncc1r, ncc1w, naread,nawrite,nc2read,nc2write,nccr,nccw,
   level ;
 char timeBuf[25] ;
 int line ;
 extern int NII, NIC, NWC, NWF, NWG ; /* in lexsubs, hash statistics */
 char buffer [2] ;
 static int nPass = 0 ; 
#ifdef ACEMBLY
 int dnaMade = 0, dnaTotal = 0 ;
 monDnaReport (&dnaMade, &dnaTotal) ;
#endif
 cacheList(0) ;
 ndwrite = 123 ; diskavail (&dfree,&dused, &plus, &minus,&ndread,&ndwrite) ;
 lexavail (&vocnum,&keynum,&keynum2,&vocspace, &hspace,&tspace);
 BSstatus (&bsused, &bsalloc, &n) ; bsmemory = n ;
 BTstatus (&btused, &btalloc, &n) ; bsmemory += n ;
 arrayStatus (&aMade, &aUsed, &aAlloc, &aReal) ;
 aStatus (&naread, &nawrite) ;
 nmessalloc = messAllocStatus (&messalloctot) ;
 blockavail(&bused,&bpinned,&bfree,&bmodif) ; 
 blockStatus (&nc1r, &nc1w, &ncc1r, &ncc1w) ;
 cacheStatus(&cacheused,&cachel,&cachek,&cachemod,&nc2read,&nc2write,&nccr,&nccw) ;
 line = 1 ;

 myText("\n************************************************\n",1,1) ;
 myText (messprintf ("  %s, %s", aceGetVersionString(), aceGetLinkDateString()), 2, line++) ;
 line++ ;
 myText (messprintf ("Data directory %s, release %d-%d",
			sessionFilName ("", 0, 0),
			thisSession.mainDataRelease,
			thisSession.subDataRelease), 2, line++) ;
 line++ ;
 myText(messprintf("Session : %d, User %s,  Write Access %s ",
		      thisSession.session, thisSession.name,
		      (isWriteAccess() ? "Yes" : "No" )  ) ,
                      12,line++);
 myText(messprintf
	   ("Global Address %d", thisSession.gAddress), 12, line++) ;
 level = freesettext(buffer,"") ;
 myText(messprintf("Stream level %d", level),12,line++) ;
 freeclose(level) ;
 line++ ;
 myText
   (messprintf
    ("Disk: %ld blocks used,  %ld blocks available.",
                     dused,dfree),2,line++) ;
 myText
   (messprintf
    ("        : %ld blocks allocated and %ld freed in this session.",
                     plus, minus),2,line++);
 line++ ;
 myText(messprintf(
   "Lexiques: %lu classes, %lu kb allocated", vocnum,tspace),2,line++);
 myText(messprintf(
   " %lu keys, %lu in lex2, %lu kb of names %lu hashing keys",
          keynum,keynum2,vocspace/1024, hspace),3,line++);
 myText(messprintf("Hash info: nii = %d, nic = %d, nwc = %d, nwf = %d, nwg = %d",
		       NII  , NIC  , NWC  , NWF  , NWG) , 2, line++) ;

 line++ ;
 blockavail(&bused,&bpinned,&bfree,&bmodif);

 indexV =  bIndexStatus (&indexT, &indexKb) ;
 if (!indexV)
   myText (messprintf("No Indexing (please read models and avoid using code prior to 4_7)"), 2, line++) ;
 else
   myText (messprintf("Indexing Version %d, %d tables using %d kb", indexV, indexT, indexKb), 2, line++) ;
 myText (messprintf("Cache 1: %d kb", blockMax()), 2, line++) ;
 myText(messprintf(
		   "Current content: %d blocks used, %d modified, %d locked, %d free",
		   bused,bmodif,bpinned,bfree),
	4, line++) ;
 myText(messprintf(
		   "Usage : %d obj served, %d saved, %d read from disk, %d saved to disk",
		   nc1r, nc1w, ncc1r, ncc1w),
	4,line++);

 myText(messprintf("Cache 2: %d slots, %d kbytes allocated", cacheused,bsmemory/1024),2,line++);
 myText(messprintf("Current content:  %d obj locked, %d obj known, %d mod.",
                     cachel,cachek,cachemod), 4, line++);
 myText(messprintf(
		   "Usage : %d obj served, %d saved, %d read from cache1, %d saved to cache1",
		   nc2read, nc2write, nccr, nccw),
	4,line++);
 myText(messprintf("memory cells: %d/%d BS + %d/%d BT used/allocated",
		      bsused,bsalloc,btused,btalloc), 4, line++) ;

 line++ ;
 myText(messprintf
   ("Arrays (Including the lexiques):"),
	   2, line++) ;
 myText(messprintf
   ("%ld made, %ld current, %ld Kb allocated, %ld Kb used (0: unknown)",
		      (long)aMade, (long)aUsed, (long)aAlloc/1024, (long)aReal/1024),
	   4, line++) ;

 line++ ;
#ifdef ACEMBLY
 myText(messprintf
   ("Acembly dna cache:"),
	   2, line++) ;
 myText(messprintf
   (" %d arrays, totalling %d Kb", dnaMade, dnaTotal/1024),
	   4, line++) ;
 line ++ ;
#endif
 if (nPass) arrayReport(nPass) ;
 nPass = arrayReportMark () ;
 myText("Cache statistics: each layer serves the known objects or asks the lower layer", 2, line++) ;
 myText(messprintf ("Array: %d read from cache1, %d saved to cache1", naread, nawrite), 4, line++) ;
 myText(messprintf ("Cache2: %d served, %d saved; %d read from cache1, %d saved to cache1", nc2read, nc2write, nccr, nccw), 4, line++) ;
 myText(messprintf ("Cache1: %d served, %d saved; %d read from disk, %d saved to disk", nc1r, nc1w, ncc1r, ncc1w), 4, line++) ;
 myText(messprintf ("Disk blocks: %d served, %d saved", ndread, ndwrite), 4, line++) ;
 myText ("Overall Memory Usage (caches 1 & 2 + arrays + various)", 2, line++) ;
 if (messalloctot)
   myText(messprintf("Messalloc: %d blocks, %d kb allocated", 
			nmessalloc, messalloctot/1024),
	     4, line++);
 else
   myText(messprintf("Messalloc: %d blocks", nmessalloc), 5, line++);
 myText (timeShow(timeNow(), timeBuf, 25), 1,line++) ; 
 line++ ;
 if (0) disknewStatus(); 
 myText("\n************************************************\n",1,1) ;
}

/*************************************************************/
/*************************************************************/

static void showModel(KEY key)
{
#ifdef JUNK
  /* this methos ios much simpler and does export the system models
     but the layout of the models is a bit off which is deadly
     and the comments are lost
     
     it is equivalent to 
     find model ?mmm ; show 
  */

   dumpKeyBeautifully (key, 'h', 0) ;
#else
  FILE *f ;
  char *card, *word ;
  int level ;



   if (!(f = filopen (sessionFilName ("wspec/models", "wrm", 0), 0, "r"))) 
    { freeOut ("// Sorry, can't find the model file wspec/models.wrm") ;
      return ;
    }
  level = freesetfile (f, "") ;
  freespecial("\n") ;
  while ((card = freecard (level)))  /* look for ?className */
    if ((word = freeword()) && !strcmp (word+1, name(key)+1))
      break ;

  if (!card)
    { freeOutf ("// Sorry, model %s not found in models.wrm\n", name(key)) ;
      goto fin ;
    }

  do freeOutf ("%s\n", card) ;
  while ((card = freecard(level)) && !(*card == '?' || *card == '#')) ;

 fin:
  freeclose(level) ;
#endif
}

extern void showModelJaq (KEY key) ;	/* in model.c */

/*************************************************************/
/********* utility for client server *************************/

static int getDataVersion (void)
{
  static int version = -1 ;

  if (version == -1)
    { 
      FILE *f ;
      char *cp ;

      version = 0 ;
      cp = sessionFilName ("wspec/server", "wrm", "r") ;
      if (cp &&
	  (f = filopen (cp, 0, "r")))
	{
	  int ii, level = freesetfile(f,0) ;
	  while (freecard(level))
	    { cp = freeword() ;
	      if (cp &&
		  !strcmp (cp,"DATA_VERSION") &&
		  freeint(&ii))
		{ version = ii ; break ; }
	    }
	  freeclose (level) ;
	}
    }
  return version ;
}

/*************************************************************/

#ifdef ACEMBLY	
BOOL cFicheView (KEY key, char style, char *view, char *params)
{
  BOOL ok = FALSE ;
  vTXT      blk;
  extern char *swormView (vTXT blkp, char *dbName, char *clsName, char *objName, char style, char *view, char *params) ;

  blk = vtxtCreate () ;
  if (swormView (blk, "local", className (key) , name (key) , style, view, params)  )
    ok = TRUE ;
  else
    freeOutf ("// Error: \n")  ;

  if (vtxtPtr (blk)) freeOut (vtxtPtr (blk)) ;
  vtxtDestroy (blk) ;
  return ok ;
}

static BOOL cFlashImage (KEY gene, char style, Stack s, char *view, char *params)
{
  BOOL ok = FALSE ;
  int height = 170 ; /* not a parameter because we want to can the images */
  int width = 680 ;
  int level ;
  OBJ Gene = 0 ;
  vTXT txt = vtxtCreate () ;
  BOOL isVGene = FALSE ;

  if (! (Gene = bsCreate (gene))
      )
    goto done ;

  if (view && !strcasecmp (view, "DtHSeq"))
    {
      KEYSET mrnas = queryKey (gene, ">transcribed_gene ; >mrna") ;
      int nMrna = keySetMax (mrnas) ;
      height = 150 + 25 * nMrna ;
      keySetDestroy (mrnas) ;
    }
  else if (view && !strcasecmp (view, "DtTiling"))
    {
      width = 1100 ; 
      height = 700 ; /* synchronize with aceviewmain.c */ 
    }

  else if (view && !strcasecmp (view, "vgene"))
    width = 1200 ; /* synchronize with aceviewmain.c */
  else if (!strcmp (view, "DtVGene"))
    { 
      height = 600 ; /* not a parameter because we want to can the images */
      width = 1200 ;
      isVGene = TRUE ;
    }
  else if (!strcmp (view, "av_tg"))
    { 
      height =1200 ; /* not a parameter because we want to can the images */
      width = 1200 ;
    }
  else if (!strcmp (view, "av_tg_whole"))
    { 
      height =1200 ; /* not a parameter because we want to can the images */
      width = 1200 ;
    }
  else if (0 && !strcmp (view, "av_mrna_whole"))
    { 
      height =1200 ; /* not a parameter because we want to can the images */
      width = 600 ;
    }

  vtxtPrint (txt, "gif ; "); 
  vtxtPrintf (txt, " dimensions %d %d ;", width, height) ;
  if (!isVGene)
    {
      vtxtPrintf (txt, " display  -D %s %s %s "
		  , view  /* DtHSeq DtGLoc DtVGene */
		  , className (gene)
		  , freeprotect (name(gene))
		  ) ;
    }
  else
    {
      int  x1 = 0, x01, x2 = 0, dx, sign, seqlen ;
      int details = 0 ;
      KEY tg = 0, cosmid = 0 ;
      char *dummy ;
      OBJ Gene = bsCreate (gene), Tg = 0 ;

      if (Gene)
	{
	  if (bsFindTag (Gene, str2tag("cloud_gene")))
	    details = 2 ;
	  if (bsGetKey (Gene, str2tag("Transcribed_gene"), &tg) &&
	      (Tg = bsCreate (tg))
	      )
	    {
	      if (bsGetKey (Tg, str2tag("Genomic_sequence"), &cosmid) &&
		  bsGetData (Tg, _Covers, _Int, &seqlen) &&
		  bsGetData (Tg, _bsRight, _Text, &dummy) &&
		  bsGetData (Tg, _bsRight, _Int, &x1) &&
		  bsGetData (Tg, _bsRight, _Int, &x2)
		  )
		{
		  sign = x1>x2 ? -1 : 1 ;
		  /* xmid = (x1 + x2)/2 ;  offset = sign * (x2 - x1)/10 ; */
		  dx = x1 < x2 ? x2 - x1 : x1 - x2 ;
		  dx /= 10 ;
		  if (dx < 500) dx = 500 ;
		  if (dx > 8000) dx = 8000 ;
		  x01 = x1 ;
		  if (x1 < x2) { x1 -= dx ; x2 += dx ; }
		  else { x1 += dx ; x2 -= dx ; }
		  seqlen = sign * (x2-x1) ;
		  details = 1 ;
		  vtxtPrintf (txt, " seqget %s -coords %d %d -origin %d -view %s ;"
			      " seqdisplay -view %s ",
			      freeprotect (name(cosmid)) ,
			      x1, x2, x01,
			      details ? "av_tg_whole" : "av_tg" , /* seqget -view */
			      details ? "av_tg_whole" : "av_tg" /* seqdisplay -view */
			      ) ;
		}
	      bsDestroy (Tg) ;
	    }
	  bsDestroy (Gene) ;
	}
      }
  if (style == 'v')  vtxtPrintf (txt, "  ; svgdump -") ;
  else if (style == 'f')  vtxtPrintf (txt, "  ;  swfdump -compile -") ;
  else  vtxtPrintf (txt, "  ;  swfdump  -") ;
  level = freesettext (vtxtPtr (txt), 0) ;

  commandExecute (level, TRUE, FALSE, 0, s, 16, 0) ; 
  freeclose (level) ;
  if (stackMark (s))
    ok = TRUE ;
 done:
  bsDestroy (Gene) ;
  vtxtDestroy (txt) ;
  return ok ;
}  /* cFlashImage */

static BOOL cWebImage (KEY tg, int style, char *view, char *params)
{
  BOOL ok = FALSE ;
  char *cp ;
  int height = 2000 ; /* should be a parameter */
  int x1, x2, xmid, dx, origin, sign, offset ;
  OBJ Tg = 0 ;
  KEY cosmid = 0 ;
  vTXT txt = vtxtCreate () ;
  
  if (! (Tg = bsCreate (tg)) ||
      ! bsGetKey (Tg, str2tag ("Genomic_sequence"), &cosmid) ||
      ! bsGetData (Tg, str2tag("covers"), _Int, &dx)  ||
      ! bsGetData (Tg, _bsRight, _Text, &cp) ||
      ! bsGetData (Tg, _bsRight, _Int, &x1) ||
      ! bsGetData (Tg, _bsRight, _Int, &x2)
      )
    goto done ;

  origin = x1 ;  
  sign = x1>x2 ? -1 : 1 ;
  xmid = (x1 + x2)/2 ; offset = sign * (x2 - x1)/10 ;
  x1 = xmid + 1.414 * (x1 - xmid) - offset ;
  x2 = xmid + 1.414 * (x2 - xmid) - offset ;
  dx = x1 < x2 ? x2 - x1 : x1 - x2 ;

  vtxtPrint (txt, "gif ; "); 
  vtxtPrintf (txt, " dimensions %d %d ;", 1200, height) ;
  vtxtPrintf (txt, " seqget \"%s\" -coords %d %d -origin %d ;"
	     ,  name(cosmid), x1, x2, origin) ;

  vtxtPrint (txt, " seqactions -hide_header ;"
	     " seqcolumns -off Locator ;"
	     " seqcolumns -off Clones ;"
	     " seqcolumns -on \"genefinder\" ;"
	     " seqcolumns -off \"-genefinder\" ;"
	     " seqcolumns -off \"NDB_CDS\" ;"
	     " seqcolumns -off \"-NDB_CDS\" ;"
	     " seqcolumns -off \"3 Frame Translation\" ;"
	     " seqcolumns -on \"Trace Assembly\" ;"
	     " seqcolumns -off \"Assembly DNA\" ;"
	     " seqcolumns -off \"Summary bar\" ;"
	     " seqcolumns -on \"RNAi\" ;"
	     
	     " seqcolumns -off TRANSCRIBEDGENE ;"
	     " seqcolumns -off -TRANSCRIBEDGENE ;"
	     " seqcolumns -off -SPLICED_cDNA ;"
	     " seqcolumns -off -EXON_SUPPORT ;"
	     " seqcolumns -off EXON_SUPPORT ;"
	     " seqcolumns -off SPLICED_cDNA ;"
	     " seqcolumns -off SPLICED_cDNA_DECORATE ;"
	     
	     " seqcolumns -on \"Genes\" ;"
	     " seqcolumns -on \"-Genes\" ;"
	     " seqcolumns -off \"mRNA_class\" ;"
	     " seqcolumns -off \"-mRNA_class\" ;"
	     " seqcolumns -off \"DNA Sequence\" ;"
	     " seqcolumns -off \"Text Features\" ;"
	     " seqcolumns -off \"GENE_NAME\" ;"
	     " seqdisplay ;"
	     ) ;
  vtxtPrint (txt,
	     style ? "gif -nobox -" : "bubble -"
	     ) ;

  /*void commandExecute (int level, BOOL noSubshell, 
		     BOOL localIsInteractive,
		     FILE *fil, Stack s, int choix,
		     KEYSET activeKs)
  */
  commandExecute (freesettext (vtxtPtr (txt), 0), TRUE, FALSE, 0, 0, 16, 0) ; 

 done:
  bsDestroy (Tg) ;
  vtxtDestroy (txt) ;
  return ok ;
}  /* cWebImage */

static BOOL cWebBoxes (KEY key, char style, char *view, char *params)
{
  BOOL ok = FALSE ;
  return ok ;
}  /* cWebBoxes */

static BOOL cFicheLazyView (KEY key, char style, char *view, char *params, int preview, int iKeySet, int mask) 
{
  BOOL lazy = (style == 'x' || style == 'f' || style == 'F' || style == 'v') ? TRUE : FALSE  ;
  char *cp, buf[1000], buffer[6] ;
  KEY lazyKey = 0 ;
  Stack s = 0 ;
  int ok = 0, level ;
  int VERSION = 1, jump6 = 0 ; /* = ficheLogicVersion(); */

  if (0 && preview==2) preview = 0 ;  /* prevents saving */
  if (lazy && iskey(key) == 2)
    {
      sprintf(buf,"%s.%d.%c.%s.%s",name(key), class(key),
	      style ? style : '0', 
	      view ? view : "", 
	      params ? params : "") ;
      cp = buf ; cp-- ;
      while (*++cp) if (*cp == ' ') *cp = '_' ;
      if (style == 'f' || style == 'F')
	lexaddkey (buf, &lazyKey, _VBinary) ;
      else
	lexaddkey (buf, &lazyKey, _VFicheView) ;
	  
      if (!mask)
	{
	  s = stackGet (lazyKey) ;
	  if (!s)  /* try to get the images calculated before they where ace dumpable */
	    {
	      KEY lazyKey2 ;
	      lexaddkey (buf, &lazyKey2, _VFicheView) ;
	      s = stackGet (lazyKey2) ;
	    }
	  if (s && stackMark (s) > 6)
	    {
	      strncpy(buffer, stackText(s,0),6) ;
	      if (atoi(buffer) == VERSION)
		jump6 = 6 ; 
	    }
	  if (s)
	    ok = 2 ;
	}
      if (ok != 2 && stackExists(s)) stackDestroy (s) ;
    }
  if (!ok)
    {
      s = stackCreate (64000) ;
      if (0) freeOutf ("%s \n", name(key)) ;

      level = freeOutSetStack (s) ;
      freeOutf("%06d",VERSION) ;
      jump6 = 6 ;

      if (style == 'f' || style == 'F' || style == 'v')
	{
	  ok = cFlashImage (key, style, s, view, params) ;
	}
      else if (view && !strncasecmp (view, "gifImage", 8))
	{
	  ok = cWebImage (key, 0, view, params) ;
	}
      else if (view && !strncasecmp (view, "gifBoxes", 8))
	{
	  ok = cWebBoxes (key, 1, view, params) ;
	}
      else
	{
	  ok = cFicheView (key, style, view, params) ;
	}
      freeOutClose (level) ;
    }
  
  /* if mask, never save directly */
  if (0 && mask) preview &= ~ 0x1 ;
  if (style != 'v' && (preview & 0x1) && lazyKey && ok == 1 && sessionGainWriteAccess() && s && stackMark(s) > 2000)
    stackStore (lazyKey, s) ;
  
  if (preview == 1) /* simply store in the local database */
    {
      if (!(iKeySet%100))
	freeOutf ("%5d %s previewed %s:%s\n", iKeySet, ok ? (ok == 2 ? "KNOWN" : "NEW") : "MISSSING", className(key), name(key)) ;
    }
  else if (preview & 0x6)  /* export as ace file or export as pipe */
    {
      if (preview & 0x4) /* export an ace file */
	{
	  freeOutf ("%s %s\n", className (lazyKey), freeprotect (name(lazyKey))) ;
	  stackCursor(s, 0) ; /* first six bytes = version */
	}
      else if (preview & 0x2)  /* binary dump */
	stackCursor(s, jump6) ; /* first six bytes = version */
      if (style == 'f' || style == 'F') /* compiled flash or uncompiled .tc flash source code */
	{
	  int nn = stackMark (s) - jump6 ;
	  if (preview & 0x4)  /* ace dump */
	    aceBinaryDoDump ((unsigned char*)stackText(s, jump6), nn) ;
	  else if (preview & 0x2)  /* binary dump */
	    freeOutBinary (stackText(s, jump6), nn) ;
	}
      else 
	{
	  while ((cp = stackNextText(s)))
	    { 
	      freeOut(cp) ;
	      freeOut ("\n") ;
	    } 
	  if (preview & 0x4) /* export an ace file */
	    freeOut ("\n***LongTextEnd***\n\n") ;
	}
    }
  stackDestroy (s) ;
 
  if (preview & 0x1)
    {
      if (iKeySet &&  !(iKeySet%1000)) 
	sessionDoSave (TRUE) ;
    }
 
  return ok > 0 ? TRUE : FALSE ;
} /* cFicheLazyView */

#endif

/*************************************************************/
/*************************************************************/

#define CHECK    if (!ksNew || !keySetMax (ksNew)) \
                  { freeOut ("// Active list is empty\n") ; \
		    break ; \
	          }

#define CHECK_WRITE if (look->noWrite) /* in readonly mode */ \
	             { freeOut ("// Sorry, you do not have Write Access.\n") ; \
	               break ; \
	             } \
	            if (!sessionGainWriteAccess())	/* try to grab it */ \
	            { /* may occur is somebody else grabed it */ \
			freeOut ("// Sorry, you cannot gain Write Access\n") ; \
			messfree (cr) ; \
			break ; \
	            }

#define UNDO  { keySetDestroy(ksOld) ; ksOld =  ksNew ; ksNew = 0 ; } /* save undo info */ 

   /* set level with a freeset file command before call
      server, will toggle server commands
      isInteractive: allow messquery
      fil/s is the output flow
   */

char *user_prompt = "acedb> ";

void set_user_prompt(const char *s0)
{
  char *s = strdup(s0);
  if (s[strlen(s)-1] == '$')
    s[strlen(s)-1] = '\n';

  user_prompt = s;
}

KEY aceCommandDoExecute (AceCommand look, int level, 
			 int localIsInteractive, KEY option, int maxChar)
     /* used also by aceserver, jadeserver */
{ 
  extern BOOL  isInteractive ;  /* an acedb global, in freesubs */
  KEYSET ksNew = 0, ksOld = 0, kA = 0, ksTmp = 0 ;
  int nn = 0 , nns = 0, i, ii = 0, j, kk, exportedColumn = 1 ;
  KEY key = 0, classe = 0, *kp = 0, tag = 0, tableKey = 0 ; 
  FILE *f = 0 ;
  char timeBuf[25] ;
  char *cp = 0, *cq =0, *cr = 0, c, *href = 0 ;
  static Stack localStack = 0 ;
  Stack ksStack ;
  BOOL isLongGrep = FALSE, isPipe = FALSE;
  BOOL oldIsInteractive ;
  FREEOPT *qmenu = arrp (look->aMenu, 0, FREEOPT) ;
  SPREAD spread = 0 ;
  const char *bqlo = 0 ;
  static int isView = 0 ;

  if (look->magic != COMMAND_MAGIC)
    messcrash ("Bad look in commandDoExecute") ;

  ksNew = look->ksNew ;
  ksOld = look->ksOld ;
  ksStack = look->ksStack ;
  nns = look->nns ;
  kA = look->kA ;
  
  oldIsInteractive = isInteractive ;
  isInteractive = FALSE ;  /* prevents parsing ? in freelevelselect */
  if (!option && localIsInteractive > 0 && user_prompt ) freeOut (user_prompt) ; /* i..e only if prompt is needed */

  if (option || freelevelselect (level, &option, qmenu))
    {
    switch (option)
      {
      case 'u':   /* undo */
	if (ksOld)
	  { keySetDestroy (ksNew) ;
	    ksNew = ksOld ;
	    ksOld = 0 ;
	    freeOutf ("\n// Recovered %d objects\n", keySetMax(ksNew)) ;
	  }
	else
	  freeOut ("// No more previous list available\n") ;
	break ;
	
	/******************* Trivialities  *****************************/
	
	
      case '/':   /* ignore  comments */
	break ;
	
      case (KEY)(-1): case 'q':   /* exit */
#ifdef KS_ADDED
	{ extern KEYSET ksAdded ;
	  if (ksAdded)
	    { FILE *f = fopen ("seqused.ace","w") ;
	      keySetDump (f, 0, ksAdded) ;
	    }
	}
#endif
        option = 'q' ;
	if (localIsInteractive) freeOut ("\n\n") ;
	break ;
	
      case 'N':   /* Count, always part of stdout    */
	keySetDestroy(kA) ;
	if (ksNew)   /* only keep visible keys. follow aliases, alphasort, compress */
	  kA = keySetAlphaHeap(ksNew, keySetMax(ksNew)) ;
	freeOutf ("// Active KeySet holds %d object names\n", kA ? keySetMax(kA) : 0); 
	keySetDestroy(kA) ;
	break ;

      case '<':   /* Drop write access permanently */
	{
	  BOOL previousIsWriteAccess = isWriteAccess();
	  
	  sessionReleaseWriteAccess();
	  sessionForbidWriteAccess(); /* can never re-gain */


	/* this used to call sessionDropWriteAccess(),
	   which calls sessionDoClose(FALSE), which in turn
	   will just release write access and return
	   and then it would set the dropAccess flag to true
	   so we can never re-gain write access again, 
	   so logic is preserved */

	  if (isWriteAccess())
	    messerror ("// Write access couldn't be dropped\n");
	  else
	    {

	      if (!previousIsWriteAccess)
		freeOut("// Write access dropped permanently\n") ;
	      else
		freeOut("// Write access dropped permanently, recent modifs will not be saved\n") ;
	    }
	  break ;
	}
      case 'S':   /* Save */
	{
	  BOOL reGainWriteAccess = FALSE ;
	  
	  cp = freeword () ;
	  if (cp && !strcasecmp (cp, "-r"))
	    reGainWriteAccess = TRUE ;
	  
	  CHECK_WRITE ;
	  if (!isWriteAccess())
	    { 
	      option = 0 ;
	      freeOut("// Nothing to save\n") ;
	    }
	  else
	    {
	      freeOut("// Saving\n") ;
	      messdump ("// Saving\n") ;
	      sessionDoSave (reGainWriteAccess) ;
	      freeOut("// Saving done\n") ;
	      messdump ("// Saving done\n") ;
	    }
	  break ;
	}


      case 300:   /* Test */
	taceTest (ksNew) ;
	break ;
	
      case 'Z':   /* Status */
	switch (getOnOff())
	  {
	  case 't':
	    look->showStatus = TRUE ;
	    tStatus () ;
	    freeOut ("// Memory statistics on, some commands will dump statistics periodically\n") ; 
	    break ;
	  case 'f':
	    look->showStatus = FALSE ;
	    freeOut ("// Memory statistics off\n") ;
	    break ;
	  case '?':
	    freeOut ("// Usage: Status {on | off}\n") ;
	    break ;
	  default:
	    tStatus () ;
	    break ;
	  }
	break ;
	
      case 'J':   /* Date */
	freeOutf ("\n// %s \n",timeShow(timeNow(), timeBuf, 25)) ;
	break ;

      case 'I':   /* Time Stamps creation */
	switch (getOnOff())
	  {
	  case 't':
	    isTimeStamps = TRUE ;		/* GLOBAL for using timeStamps */
	    freeOut ("// Time-stamps now ON: Each modified data will be time stamped\n") ;
	    break ;
	  case 'f':
	    isTimeStamps = FALSE ;		/* GLOBAL for using timeStamps */
	    freeOut ("// Time-stamps now OFF\n") ;
	    break ;
	  default:
	    freeOut ("// Usage: Time_Stamps {on | off}\n") ;
	    break ;
	  }
	break ;
	
      case 'H':   /* Chrono */
	cp = freeword () ;
	if (cp && !strcasecmp (cp, "start"))
	  chronoStart () ;
	else if (cp && !strcasecmp (cp, "stop"))
	  chronoStop () ;
	else if (cp && !strcasecmp (cp, "show"))
	  chronoReport () ;
	else
	  freeOut ("// Usage: Chrono { start | stop | show}\n") ;
	break ;

      case '?':
	if ((cp = freeword ()))
	  { 
	    char *cq = messprintf("%s*",cp) ; 
	    freeOutf ("// ACEDB Server list of commands matching %s:\n", cq) ;
	    for (i = 1 ; i <= qmenu[0].key ; i++)
	      if (pickMatch (qmenu[i].text, cq))
		freeOutf ("%s\n",qmenu[i].text) ;
	  }
	else
	  {
	    freeOut ("// ACEDB Server list of commands:\n") ;
	    for (i = 1 ; i <= qmenu[0].key ; i++)
	      freeOutf ("%s\n",qmenu[i].text) ;
	    freeOut ("// ? command for options of just that command\n") ;
	  }
	break ;   
	
      case 'h':  /* help */
	/*	if (freekey (&option, helpmenu))
	  { freeforcecard (helpmenu[option].text) ;
	    cp = freeword() ;
	    helpOn (cp) ;
	  }
	  */
	if ((cp = freeword()))
	  helpOn (cp) ;
	else
	  { freeOut
	      (messprintf 
	       (
		"This program maintains an internal list of current objects. \n "
		"You can\n"
		"-Create and modify the list with the simple commands: \n"
		"       Find, Follow, Grep, LongGrep, Biblio, Is, Remove, Undo \n"
		"-Visualise in several format the content of the list with:\n"
		"       List, Show, Count, Dna, Dump\n"
		"-Perform complex queries and print relational tables with: \n"
		"       Query, Table \n"
		"-Modify the data with:\n"
		"       New, Edit, EEdit, Parse, PParse, Kill\n"
		"-See the schema and the dataload with\n"
		"       Model, Classes, Status, Time_stamps\n"
		"-Manipulate several lists a la polonaise with\n"
		"       clear, spush, spop, swap, sand, sor, sxor, sminus\n"
		" \n"
		"Commands are not case sensitive and can be abbreviated.\n"
		"They can be read from a command file with parameters (noted %%1 %%2...) with\n"
		"   @file_name parm1 parm2 parm3...\n"
		"if ACEDB_SUBSHELLS is set, lines starting with $ run as an interactive subshell.\n"
		"Everything following // on a line is ignored (comment).\n"
		"To escape any of the special characters use a backslash, eg \\@.\n\n"

		"To see the syntax of all commands type: ?\n"
		"For further help type: help topic\n"
		)) ;
	  }
	break ;
	

	/******************* Stand Alone   *****************************/
	
      case 'c':  /* Class */
	freeOut ("These are the known classes and the number of objects in each class \n") ;
	for (i=0 ; i < 256 ; i++)
	  /*
	    if (pickClass2Word(i) && !pickList[i].private && lexMaxVisible(i))
	    freeOutf ("%35s %d\n", pickClass2Word(i), lexMaxVisible(i)) ;
	  */
				/* 2nd call to lexMaxVisible() is cheap */
	  if (pickClass2Word(i) && lexMax(i)) 
	    freeOutf ("%35s %d\n", pickClass2Word(i), lexMax(i)) ;
	break ;
	
      case 'z':   /* Z class-model */
	cp = freeword() ;
	if (!cp) 
	  { freeOut ("// Usage Model [-jaq] Class : Shows the model of given class,\n"
		     "//accepts wildcards\n") ;
	    break ;
	  }
	{ BOOL isJaq = FALSE ;
	  BOOL found = FALSE ;

	  if (!strcmp (cp, "-jaq"))
	    { isJaq = TRUE ;
	      cp = freeword () ;
	      if (!cp) 
		{ freeOut ("// Usage Model [-jaq] Class : Shows the model of given class,\n"
			   "//accepts wildcards\n") ;
		  break ;
		}
	    }

	  i = 0 ;
	  for (key = 0 ; lexNext(_VModel, &key) ; )
	    if (pickMatch (name(key)+1, cp))
	      { if (isJaq)
		  showModelJaq (key) ;
		else
		  showModel (key) ;
		found = TRUE ;
	      }
	
	  if (!isJaq)
	    for (key = 0 ; lexNext(_VClass, &key) ; )
	      if (pickMatch (name(key), cp) )
		{ KEY classe = key ;
		  unsigned char mask ;   
		  
		  pickIsA (&classe, &mask) ;
		  if (mask)
		    { found = TRUE ;
		      freeOutf ("// %s is a sub class of %s, here is its filter\n",
				name(key), pickClass2Word(classe)) ;
		      { OBJ obj = bsCreate (key) ;
			if (obj) 
			  { if (bsGetData (obj, _Filter, _Text, &cp))
			      freeOut ( messprintf ("// %s\n\n", cp)) ;
			    bsDestroy (obj) ;
			  }
		      }
		    }
		}

	  if (!found)
	    freeOutf("// Sorry, I can t find a model for : %s\n", cp) ;
	}
	break ;
	
      case 'R':      /* read models */
	{
	  BOOL reindex = FALSE ;
	  cp = freeword() ;
	  if (cp)
	    {
	      if (!strcmp(cp,"-index"))
		{  reindex = TRUE ; }
	      else
		{ freeOut ("// Usage: Read-models [-index]\n") ; break ; }
	    }
	  CHECK_WRITE;
	  if (reindex)
	    bIndexInit(BINDEX_AFTER_READMODELS) ;
	  else
	    readModels() ;
	  freeOut ("// New models read, save your work before quitting\n") ;
	}
	break ;
	
      case 'F':   /* New Class format */
	cp = freeword() ;
	if (!cp) 
	  { freeOutf ("// Usage New Class Format: Construct a new name"
		      " according to C-language printf(format,n) \n") ;
	    freeOut("// Example: New Oligo xxx.\\%%d\n") ;
	    break ;
	  }
	if (!lexword2key(cp, &key, _VMainClasses))
	  { freeOutf ("//! Unkown class %s\n",cp) ;
	    break ;
	  }
	classe = KEYKEY(key) ;
	if (pickList[classe].protected)
	  { freeOutf ("//! Protected class\n",cp) ;
	    break ;
	  }
	
	cp = freeword() ;
	if (!cp) 
	  {
	    freeOut ("// Usage New Class Format: Construct a new name"
		     " according to C-language printf(format,n) \n") ;
	    freeOut ("// Example: New Oligo xxx.\\%%d\n") ;
	    break ;
	  }
	localStack = stackReCreate(localStack, 80) ;
	pushText(localStack, cp) ;
	
	CHECK_WRITE ;
	i = -2 ; /* important flag */
	cp = stackText(localStack,0) ;
	while (*cp && *cp != '%') cp++ ;
	if (*cp++ == '%')
	  { while(*cp && isdigit((int)*cp)) cp++ ;  
	    if (*cp++ == 'd')
	      { i = 1 ;
		while (*cp && *cp != '%') cp++ ;
		if (!*cp)
		  goto okformat ;
	      }
	    freeOut ("// Only allowed format should contain zero or once %%[number]d \n") ;
	    break ;
	  }
      okformat:
	key = kk = 0 ;
	cp = stackText(localStack,0) ;
	while (i < lexMax(classe) + 3 && /* obvious format obvious problem */
	       !(kk = lexaddkey(messprintf(cp,i), &key, classe)))
	  if (i++ < 0) break ;  /* so i break at zero f format has no %d */
	if (kk)
	  freeOutf ("%s  \"%s\"\n", className(key), name(key)) ;
	else
	  freeOutf ("// Name %s allready exists\n", cp) ;
	break ;
	
	/******************* Parsers      *****************************/

      case 'K':  /* read keyset from stack */
	nn = 0 ;
	cp = freepath() ;
	cr = cp ? strnew(cp, 0) : 0 ;
	if (!cp && !(look->choix & 8)) /* server should never open  stdin */
	  nn = freesetfile(stdin,0) ;
	else if (!strcmp (cp,"="))
	  {
	    freenext() ;
	    cp = freepos() ;
	    if (cp && *cp) 
	      cr = strnew(cp, 0) ;
	    else cr = strnew("   ", 0) ; /* reset active list to empty */
	    nn = freesettext (cr, 0) ;
	    freespecial ("\n\"\\/;") ;
	  }
	else
	  { 	
	    cr =  strnew(cp, 0) ;
	    f = fopen (cp, "r") ;
	    if (!f)
	      freeOutf ("// Sorry, I could not open file %s\n", cr) ;
	    else
	      nn = freesetfile (f, "") ;
	    messfree (cr) ; 
	    UNDO ;
	    ksNew = keySetCreate () ; // mieg 2007_01_29 no file impliies non active object
	  }
	if (nn)
	  {
	    UNDO ;
	    ksNew = keySetRead (nn) ; /* will close level */
	    freeOutf ("// Recognised %d objects\n", keySetMax(ksNew)) ;
	    messfree (cr) ;
	  }
	break ;
	
      case '_':  /* server parser */
	CHECK_WRITE ;
	UNDO ;
	ksNew = arrayHandleCreate (50, KEY, look->h) ;
	parseAceCBuffer (freepos(), 0, 0, ksNew) ;
	parseKeepGoing = TRUE ;	  
	freespecial ("\n/\"\\\t") ;
	parseLevel (level, ksNew) ;
	parseKeepGoing = FALSE ;
	break ;

      case 'P':
	parseKeepGoing = TRUE ;	 
	/* fall thru */
      case 'p':
	kk = 0 ;
	look->noXref = FALSE ;
	while (TRUE)
	  {
	    freenext() ;
	    if(freestep('=')) 
	      { 
		kk = 1 ; freenext() ; 
		cp = freepos() ;  
		cr = cp ? strnew(cp, 0) : 0 ;
		break ;
	      }
	    if (freestep('-'))  /* treat her all -options in random order */
	      { cp = freeword () ;
	        if (cp && !strcmp(cp,"NOXREF"))
		  look->noXref = TRUE ;
		continue ;
	      }
	    cp = freepath () ; /* see filsubs.c */
	    cr = cp ? strnew(cp, 0) : 0 ;
	    break ;
	  }
	CHECK_WRITE ;
	UNDO ;
	ksNew = arrayHandleCreate (50, KEY, look->h) ;

	f = 0 ; isPipe = FALSE ;
	if (cr && !kk)
	{
	  /*
	  * they said a file name
	  */
	  f = fopen(cr, "r") ;
	  if (f)
#if defined(WIN32) /* WIN32 has pkunzip? */
	    /* Assumes that xyz.zip file (only) contains a xyz.ace file */
	    {
	      cq = cr + strlen(cr) - 4;
	      if (cq > cp && !strcmp(cq, ".zip"))
		{
		  filclose (f) ; f = 0 ;
		  cq[0] = '\0' ; /* chop off the .zip in cr, leaving root path */
		  if( callScript( "pkunzip",messprintf( " -o %s >%s.err", cr, cr)) != -1)
		    { 
		      isPipe = TRUE ;
		      f = callScriptPipe ("type", messprintf(" %s.ace", cr) ) ;
		    }
		}
	    }
#else  /* UNIX has zcat and gunzip?	*/
	  {
	    cq = cr + strlen(cr) - 2; 
	    if (cq > cr && !strcmp(cq, ".Z")) /* a compressed file */
	      {
		filclose (f) ; f = 0 ;
		f = callScriptPipe ("zcat",cr) ; isPipe = TRUE ;
	      }
	    else if (cq-- && cq > cr && !strcmp(cq, ".gz"))
	      {
		filclose (f) ; f = 0 ;
		f = callScriptPipe ("gunzip",messprintf(" -c %s", cr)) ; isPipe = TRUE ;
	      }
	    
	  }
#endif
	}
	else if (!cr && !kk)
	  { 
            if (!(look->choix & 8)) /* server should never open  stdin */
	      { f = stdin;
#if defined(WIN32)
	      freeOut ("// Type in your data in .ace format, then CTRL-Z<enter> to finish\n") ;
#else
	      freeOut ("// Type in your data in .ace format, then CTRL-D to finish\n") ;
#endif
	      }
	  }
	if (kk)
	  { 
	    if (cr && *cr)
	      {
		BOOL oldXref = XREF_DISABLED_FOR_UPDATE ;
		int nn = freesettext(cr,"") ;
		freeOut ("// Accepting data on command line\n") ;
		freespecial ("\n\t\"\\;/") ; /* recognize ';' as separator, allow continuation lines  */
		XREF_DISABLED_FOR_UPDATE = look->noXref ;
		parseLevel (nn, ksNew) ;
		XREF_DISABLED_FOR_UPDATE = oldXref ;
	      }
	  }
        else if (f) 
	  { 
	    BOOL oldXref = XREF_DISABLED_FOR_UPDATE ;

	    if (cr)
	      {
		messdump("// Parsing file %s %s\n",
			 look->noXref ? "NOXREF" : "",
			 cr) ;
		freeOutf("// Parsing file %s %s\n",
			 look->noXref ? "NOXREF" : "",
			 cr) ;
	      }
	    else
	      messdump("// Parsing stdin\n") ;
	    XREF_DISABLED_FOR_UPDATE = look->noXref ;
	    if (isPipe)
	      parsePipe(f, 0, ksNew) ; /* will pclose f */
	    else
	      parseNamedFile(cr, f, 0, ksNew) ; /* will pclose f */
	    
	    XREF_DISABLED_FOR_UPDATE = oldXref ;
	  }
	else if (cr)
	  freeOutf ("// Sorry, I could not open file %s\n", cr) ;
	else
	  freeOut ("// Use \"aceclient host -ace_in\" to read data on the client side") ;

	parseKeepGoing = FALSE ;
	messfree(cr) ;
	break ;
	
	/******************* Constructors  *****************************/
	
      case 'C':   /* Clear */
	UNDO ;
	ksNew = arrayHandleCreate (50, KEY, look->h) ;
	break ;
	
      case 'j':  /* Dump */
	{ BOOL isSplit = FALSE ;
	  int z;
	  BOOL doDump = TRUE;
	  char dumpDir[DIR_BUFFER_SIZE] = "" ;

	  while ((cp = freeword()))
	    if (!strcmp (cp, "-s"))
	      isSplit = TRUE ;
	    else if (!strcmp(cp, "-T"))
	      dumpTimeStamps = TRUE ;
	    else if (!strcmp(cp, "-C"))
	      dumpComments = FALSE ;
	    else
              { 
		if (filName(cp, "", "rd"))/* if it is a directory continue*/
		  { strcpy (dumpDir,cp) ;
		    z = strlen( dumpDir) ;
		    if (dumpDir[z-1] != '/')
		      { dumpDir[z] = '/' ;
			dumpDir[z+1] = 0 ;
		      }
		  }
	        else
		  { messerror ("Could not open directory %s", cp) ;
		    doDump = FALSE ;
		  }
	      }
	  if (doDump)
	    dumpAllNonInteractive (dumpDir, isSplit) ;
	}
	break ;

      case 'n':   /* FIND */
	{
	  int clientMax = 0 ;
	  cp = freeword() ;
	  if (cp && !strcmp (cp, "-max"))
	    {
	      if (!freeint (&clientMax) || clientMax <= 0)
		{
		  freeOut ("// Usage Find [-max n] Class [name]: Construct a new current list \n"
			   " // please give a positive integer after \"-max\" or avoid this parameter") ;
		  break ; 
		}
	      cp = freeword() ;
	    }
	  if (!cp) 
	    { freeOut ("// Usage Find [-max n] Class [name]: Construct a new current list \n"
		       " // formed with all objects belonging to class \n"
		       " // Example: New Sequence a\\%%.b   // may give a7.b") ;
	    break ;
	    }
	  localStack = stackReCreate(localStack, 80) ;
	  classe = 0 ;
	  lexword2key (cp, &classe, _VClass) ;
	  pushText(localStack, cp) ;
	  freenext () ;
	  i = 0 ;  /* no name */
	  if ((cp = freepos()) && *cp)
	    { 
	      BOOL isProtected = FALSE ;
	      if (*cp == '"') { isProtected = TRUE ; } /* remove flanking " if there */
	      cq = cp + strlen(cp) - 1 ;
	      while (cq > cp && (*cq == ' ' || *cq == '\t')) *cq-- = 0 ;
	      i = 1 ; /* a name */
	      cq = cp - 1 ;
	      while (*++cq)
		if (*cq == '?' || *cq == '*') i = 2 ; /* a regexp */
	      if (classe && i == 1)
		{
		  stackClear (localStack) ;
		  pushText (localStack, isProtected ? freeunprotect (cp) : cp) ;
		}
	      else
		{
		  catText (localStack, " IS ") ;
		  catText (localStack, isProtected ? cp : freeprotect(cp)) ; 
		}
	    }	    
	  UNDO ;
	  switch (i)
	    {
	    case 0: case 2:
	      ksNew = query(0, messprintf("FIND %s",  stackText(localStack, 0))) ;
	      break ;
	    case 1:  /* exact word */
	      ksNew = arrayHandleCreate (50, KEY, look->h) ;
	      if (lexword2key (stackText(localStack, 0), &key, classe))
		keySet (ksNew, 0) = key ; 
	      break ;
	    }
	  if (clientMax > 0 && clientMax < keySetMax (ksNew))
	    {
	      keySetDestroy(kA) ;
	      if (ksNew) 
		kA = keySetAlphaHeap(ksNew, clientMax) ;
	      keySetDestroy (ksNew) ;
	      ksNew = keySetCopy (kA) ;
	      keySetSort (ksNew) ;
	    }
	  freeOutf ("\n// Found %d objects in this class\n", keySetMax(ksNew)) ;
	}
	break ;
      case 'x':   /* Performs a complex query */
	if (!(cp = freewordcut("",&c)))
	  break ;
	cr = strnew (cp, 0) ;
	UNDO ;
	ksNew = query(ksOld,cr) ;
	freeOutf ("\n// Found %d objects\n", keySetMax(ksNew)) ;
	messfree (cr) ;
	break ;
	
      case 'G':
	isLongGrep = TRUE ;
      case 'g':   /* Match == Grep */
	UNDO ;
	if (!(cp = freeword()))
	  break ;
	if (!strcmp(cp, "-active"))
	  { ksTmp = ksOld ; 
	    freenext() ;
	  }
	else
	  { freeback() ;
	    ksTmp = 0 ;
	  }
	if (!(cp = freewordcut("",&c)))
	  break ;
	cr = strnew (messprintf("*%s*",cp), 0) ;
	freeOutf ("// I search for texts or object names matching %s\n", cr) ;
	ksNew = queryGrep(ksTmp, cr) ; ksTmp = 0 ;
	if (isLongGrep)
	  { KEYSET ks1 = ksNew, ks2 = longGrep(cr) ;
	    ksNew = keySetOR(ks1, ks2) ;
	    keySetDestroy(ks1) ;
	    keySetDestroy(ks2) ; 
	    isLongGrep = FALSE ;
	  }
	freeOutf ("\n// Found %d objects\n", keySetMax(ksNew)) ;
	messfree (cr) ;
	break ;
#ifdef ACEMBLY	
      case 'V':   /* View */
	{
	  char *view = 0 ;
	  char *params = 0 ;
	  int cl = 0, mask = 0 ;
	  char *kName = 0 ;
	  char *fName = 0 ;
	  char style = 0 ;
	  FILE *ff = 0 ;
	  int oLevel = 0 ;
	  BOOL ficheOk = FALSE ;
	  int preview = 0 ;
	  
	  isView = 1 ;

	  freenext() ; key = 0 ;
	  preview = 2 ; /* default behaviour is to dump for the web */
	  while (freestep('-'))
	    { 
	      freenext () ;
	      cp = freeword() ;
	      if (!cp) break ;
	      else if (!strcasecmp(cp, "v"))
		{ 
		  cp = freeword () ;
		  if (!cp || *cp == '-')
		    {
		      freeOutf ("// View -v %s : Missing view-name, sorry", cp ? cp : "") ;
		      key = 0 ;
		      break ;
		    }
		  else
		    view = strnew (cp, 0) ;
		}
	      else if (!strcasecmp(cp, "active"))
		cl = -1 ;
	      else if (!strcasecmp(cp, "preview"))
		preview = 1 ; /* store, do not export */
	      else if (!strcasecmp(cp, "preview1"))
		preview = 1 ; /* store, do not export */
	      else if (!strcasecmp(cp, "preview2"))
		preview = 2 ; /*  binary export */
	      else if (!strcasecmp(cp, "preview3"))
		preview = 3 ; /*  store and binary export */
	      else if (!strcasecmp(cp, "preview4"))
		preview = 4 ; /* do not store, just .ace export */
	      else if (!strcasecmp(cp, "c"))
		{ 
		  cp = freeword () ;
		  if (!cp || *cp == '-')
		    {
		      freeOutf ("// View -c %s : Missing class-name, sorry", cp ? cp : "") ;
		      key = 0 ;
		      break ;
		    }
		  else
		    {
		      KEY classKey = 0 ;

		      lexword2key (cp, &classKey, _VClass) ;
		      cl = superClass (classKey) ;
		      if (!cl)
			{
			  freeOutf ("// View -c %s : Unknown class, sorry", cp) ;
			  key = 0 ;
			  break ;
			}
		    }
		}
	      else if (!strcasecmp(cp, "n"))
		{ 
		  cp = freeword () ;
		  if (!cp || *cp == '-')
		    {
		      freeOutf ("// View -n %s : Missing object-name, sorry", cp ? cp : "") ;
		      key = 0 ;
		      break ;
		    }
		  else
		    kName = strnew (cp, 0) ;
		}
	      else if (!strcasecmp(cp, "o"))
		{ 
		  cp = freeword () ;
		  if (!cp || *cp == '-')
		    {
		      freeOutf ("// View -f %s : Missing file-name, sorry", cp ? cp : "") ;
		      key = 0 ;
		      break ;
		    }
		  else
		    {
		      fName = strnew (cp, 0) ;
		      ff = filopen (fName, 0, "w") ;
		      if (!ff)
			{
			  freeOutf ("// View -f %s : This file cannot be opened for writing, sorry", cp ? cp : "") ;
			  key = 0 ;
			  break ;
			}
		      oLevel = freeOutSetFile (ff) ;
		    }
		}
	      else if (!strcasecmp(cp, "gf"))
		{ 
		  cp = freeword () ;
		  if (!cp || *cp == '-')
		    {
		      freeOutf ("// View -gf %s : Missing file-name, sorry", cp ? cp : "") ;
		      key = 0 ;
		      break ;
		    }
		  else
		    {
		      char *s, buf[4096] ; 

		      while (*cp == '/' || *cp == '\\') cp++ ; /* forbid absolute paths */
		      while (*cp == '.') cp++ ; /* forbid motions */
		      for (s = cp ; *s ; s++)
			{
			  if (*s == '/' && *(s+1)== '/') *(s+1)= '_' ;
			  if (*s == '\\') *(s+1)= '_' ;
		          if (*s == '/' && *(s+1)== '.') *(s+1)= '_' ;
			}
		      fName = strnew (cp, 0) ;
		      ff = filopen (messprintf ("goldFiles/%s", fName), 0, "r") ;
		      if (!ff)
			{
			  freeOutf ("// View -gf %s : ERROR file not found, sorry", cp ? cp : "") ;
			  key = 0 ;
			  break ;
			}
		      /*
		       * if an external gold file is specified, send it
		       */
		      else
			{
			  while ((s = fgets (buf,4095,ff)))
			    freeOut (s) ;
			}
		      filclose (ff) ; ff = 0 ;
		    }
		  goto viewDone ;
		}
	      else if (!strcasecmp(cp, "svg"))
	  	style = 'v' ;
	      else if (!strcasecmp(cp, "flash"))
	  	style = 'f' ;
	      else if (!strcasecmp(cp, "tc"))
	  	style = 'F' ;
	      else if (!strcasecmp(cp, "xml"))
	  	style = 'x' ;
	      else if (!strcasecmp(cp, "asn"))
	  	style = 's' ;
	      else if (!strcasecmp(cp, "refseq"))
	  	style = 'r' ;
	      else if (!strcasecmp(cp, "mask"))
		freeint (&mask) ;
	      else  /* eat everything else */
		{ 
		  freeback () ;
		  cp = freepos () ;
		  params = strnew (cp, 0) ;
		  break ;
		}
	    }
	  if ((cl > 0 && !kName) || (cl <= 0 && kName))
	     {
	       freeOutf ("// View -c class -n name : You must provide both -c and -n or none , sorry", cp ? cp : "") ;
	       key = 0 ;
	     }
	  if (cl > 0 && 
	      !lexword2key (kName, &key, cl) &&
	      iskey (key) != 2)
	    {
	      freeOutf ("// View -c %s -n %s : Unknown object sorry", className (cl << 24), kName) ;
	      key = 0 ;
	    }
	 
	  if (!key && cl != -1)
	    freeOutf("// Usage View [-v view] [-active | -c class -n name] [-preview] [-xml | -asn] [-o outputfile] [-p params to end of line]: View active set\n") ;
	  {
	    int zzz=1 ; /* 1 do all */
	    KEY key2 ;
	    char *view2 ;
	    
	    if (cl > 0 && cl == _VGene && view && !strcasecmp (view, "fiche"))
	      view2 = "fgene" ;
	    else
	      view2 = view ;
	    
	    if (cl == -1)
	      for (i = 0 ; i < keySetMax (ksNew) ; i++)
		{
		  key2 = keySet(ksNew, i) ;
		  if (!strcmp(name(key2), "3L15"))
		    zzz = 1;
		  if (mask && mask != 999 && ((key2 % 64) != (mask % 64))) continue ;
		  if (zzz) ficheOk |= cFicheLazyView (key2, style, view2, params, preview, i, mask) ;
		}
	    else if (key)
	      ficheOk |= cFicheLazyView (key, style, view2, params, preview, 0, mask) ;
	  }
	  freeOut("\n") ;
      viewDone:
	  messfree (view) ;
	  messfree (params) ;
	  messfree (kName) ;
	  messfree (fName) ;
	  if (oLevel)
	    {
	      freeOutClose (oLevel) ;
	      filclose (ff) ; ff = 0 ;
	    }
	  isView = 2 ;

	  break ;
	}

      case 500: /* optimised web query creates an active set of transcribed_genes */
	UNDO ;
	{
	  extern KEYSET cdnaOptimisedWebQuery (KEYSET ksActive) ;
	  ksNew = (KEYSET) cdnaOptimisedWebQuery (ksOld) ;
	}
	break ;

      case 600: /* [-o file] type run :: export an html page using google API */
	UNDO ;
	{
	  extern KEYSET magicProcessAcedbRequest (KEYSET ksOld) ;
	  ksNew = magicProcessAcedbRequest (ksOld) ;
	}
	break ;

#endif

      case 'T':
	kk = 0 ; cp = 0 ; f = 0 ;
	exportedColumn = 1 ;
	tableKey = 0 ;  /* indicates to save the table def as an object in _VTable */
	look->minN = 0 ;
	look->maxN = -1 ;
	look->beginTable = TRUE ;
	look->outfile = 0 ;
	freenext() ;
	while (freestep('-'))
	  { 
	    freenext () ;
	    cp = freeword() ;
	    if (!cp) break ;
	    if (!strcasecmp(cp, "active"))
	      kk |= 1 ;
            else if (!strcasecmp(cp, "title"))
	      kk |= 2 ;
            else if (!strcasecmp(cp, "j"))
	      kk |= 4 ;
            else if (!strcasecmp(cp, "p"))
	      kk |= 8 ;
            else if (!strcasecmp(cp, "h"))
	      kk |= 16 ;           
            else if (!strcasecmp(cp, "a"))
	      kk |= 32 ;           
            else if (!strcasecmp(cp, "n"))
	      kk |= 64 ;           
            else if (!strcasecmp(cp, "x"))
	      kk |= 256 ;           
            else if (!strcmp(cp, "C"))
	      kk |= 128 ;           
	    else if (!strcmp(cp, "c"))
	      { 
		if (freeint (&i) && i >= 0)
		  look->maxN = i ;
	      }
	    else if (!strcasecmp(cp, "b"))
	      { /* b = 1 means start at line 1 = the beginning, forget that */
		if (freeint (&i) && i >= 2)
		  look->minN = i - 1 ;
	      }
	    else if (!strcasecmp(cp, "xcol"))
	      { /* b = 1 means start at line 1 = the beginning, forget that */
		if (freeint (&i) && i >= 1)
		  exportedColumn = i ;
	      }
	    else if (!strcasecmp(cp, "s"))
	      { 
		cp = freeword () ;
		if (cp && *cp == '-')
		  freeback () ;
		else if (cp)
		  lexaddkey (cp, &tableKey, _VTable) ;
	      }
	    else if (!strcasecmp(cp, "f"))
	      { 
		if ((cp = freepath ()) &&
		    (!(f = filopen(cp, 0, "r")))
		    )
		  { 
		    freeOutf ("// -f %s failed\n",cp ? cp : "Missing Filename") ;
		    while (freeword()) ; /* empty the line, this will create another error */
		  }
	      }
	    else if (!strcasecmp(cp, "o"))
	      { 
		if ((cp = freepath ()) &&
		    (!(look->outfile = fopen(cp,"w")))
		    )
		  { 
		    freeOutf ("// out_fileopen %s failed\n",cp ? cp : "Missing Filename") ;
		    while (freeword()) ; /* empty the line, this will create another error */
		  }
	      }
	    else if (!strcasecmp(cp, "bqlo"))
	      { 
		if ((cp = freepath ()))
		  {
		    if (!(look->outfile = fopen(cp,"w")))
		      { 
			freeOutf ("// bqlo fileopen %s failed\n",cp ? cp : "Missing Filename") ;
			while (freeword()) ; /* empty the line, this will create another error */
		      }
		    else
		      bqlo = strnew (cp, 0) ;
		  }
	      }
	    else if (!strcasecmp(cp, "href"))
	      { 
		cp = freeword () ;
		if (cp && *cp == '-')
		  freeback () ;
		else if (cp)
		  href = strnew (cp, 0) ;
	      }
	    
	    freenext () ;
	  }      
        cp = freeword() ;
	if (cp && !strcmp (cp,"="))
	  { 
	    Stack ss1 = 0 ;
	    
	    freenext() ;
	    cp = freepos() ;
	    if (cp && *cp) 
	      cr = strnew(cp, 0) ;
	    else cr = strnew("   ", 0) ; /* reset active list to empty */

	    spread = spreadCreate() ;
	    ss1 = stackCreate(100) ;
	    pushText (ss1,cr) ;

	    spreadDoReadDefinitions (spread, 0, 0, ss1, 0, FALSE) ; 
	    key = 3 ; /* needed on next call to spreadDoReadDefinitions */
	    if (!spreadCheckConditions(spread))
	      {
		spreadDestroy(spread);

		freeOut ("// Sorry, Bad table definition,  please debug it in graphic mode\n") ;
		break ;
	      }
	    UNDO ;
	    if (kk & 1)
	      ksNew = keySetCopy (ksOld) ;
	    else
	      ksNew = spreadGetKeyset (spread) ;

	    if (!ksNew)
	      { 
		freeOut ("// This table is empty in its first column, it may be misconstructed\n") ;
		ksNew = arrayCreate (50, KEY) ;
	      }
	    keySetSort (ksNew) ;
	    stackDestroy (ss1) ;
	    goto readDirectTableDef ;
	  }
	
	if (!f && (!cp || !*cp))
	  { freeOut ("// Table-Maker  [-active] [-a | -h | -j | -p | -C] [-b begin] [-c count] [-title] [-o output_file] [-n name |  [-f] table_def_file | = command\ncommand] [params] : \n\t[Search active list] [ace,human,java,perl,C style] table_defition_file [parametrs (noted %%d in the def]") ;
	  break ;	  
	  }

	if (look->spread)
	  spreadDestroy(look->spread);

	key = 0 ;
	if (kk & 64)
	  {  if (!lexword2key (cp, &key, _VTable) ||
		 iskey (key) != 2)
	    { freeOutf ("Sorry, I could not find table %s\n", cp) ;
	      break ;
	    }
	  }
	else
	{ 
	  if (!f) 
	    { 
	      freeOutf ("Sorry, I could not open file %s\n", cp) ;
	      break ;
	    }
	}
	
	freenext() ;   /* get parameters */
	cp = freepos() ;
	messfree (cr) ;
	if (cp)
	  cr = strnew (cp, 0) ;
	else 
	  cr = 0 ;
	
	spread = spreadCreate() ;
	
	if (tableKey)
	  {
	    CHECK_WRITE ;
	    spreadDoReadDefinitions (spread, key, f, 0, cr, FALSE) ; /* will close f */
	    if (spread)
	      { 
		spreadDoSaveInObj (spread, tableKey) ;

		spreadDestroy (spread);

		freeOut ("\n//Table saved\n") ;
	      }
	    option = 'T' ;
	    break ;
	  }
		
	if (bqlo)
	  {
	    spreadDoReadDefinitions (spread, key, f, 0, cr, FALSE) ; /* will close f */
	    if (spread)
	      { 
		spreadDoExportBql (spread, bqlo, TRUE) ;
		spreadDestroy (spread);
	      }
	    option = 'T' ;
	    break ;
	  }
	
	spread->exportedColumn = exportedColumn ;
	spread->exportedKeySet = keySetCreate () ;

      readDirectTableDef:
	if ( kk & 256) spread->style = 'X' ;  /* html, no questions asked */
	if ( kk & 32) spread->style = 'a' ;
	if ( kk & 16) spread->style = 'h' ;
	if ( kk & 8) spread->style = 'p' ;
	if ( kk & 4) spread->style = 'j' ;
	if ( kk & 128) spread->style = 'C' ;
	if ( kk & 2) spread->showTitle = TRUE ;
	
	spread->href = strnew (href, 0) ;
	messfree (href) ;
	keySetDestroy (kA) ;
	if (kk & 1) /* -active */
	  { 
	    if (key || f || cr) /* f == 0 in case of direct table */
	      spreadDoReadDefinitions (spread, key, f, 0, cr, FALSE) ; /* will close f */
	    spread->precomputing = FALSE ; /* cannot store ! */
	    spread->isActiveKeySet = TRUE ;
	    keySetDestroy (ksTmp) ;
	    ksTmp = spreadFilterKeySet (spread, ksNew) ;
	    kA = keySetAlphaHeap(ksTmp, keySetMax(ksTmp)) ;
	    keySetDestroy (ksTmp) ;
	  }
	else /* ! -active */
	  { 
	    UNDO ;
	    kA = spreadGetPrecalculatedKeySet (spread, key, cr) ;
	    if (!kA)
	      { 
		if (key || f || cr)  /* f == 0 in case of direct table */
		  spreadDoReadDefinitions (spread, key, f, 0, cr, FALSE) ;
		kA = spreadGetKeyset (spread) ;
	      }
	    else
	      filclose (f) ; 
	    f = 0 ;
	    if (!kA)
	      { freeOut ("// This table has a problem, please debug it in graphic mode") ;
		kA = arrayHandleCreate (50, KEY, look->h) ;
	      }
	    ksNew = arrayHandleCopy (kA, look->h) ;
	    keySetSort (ksNew) ;
	  }
	
	look->spread = spread ;
	look->lastSpread = look->minN > 0 ? look->minN  : 1 ;
	if (look->maxN >= 0) look->maxN += look->minN ;
	look->cumulatedTableLength  = 0 ;
	
	messfree (cr) ;
	look->lastCommand = 'B' ;
	/* fall thru to case encore */
	
      case 'B': /* encore: suite du dump des tables */
	{
	  SPREAD spread = look->spread ;
	  int last = look->spread ? look->lastSpread : 0 ;
	  
	  if (look->lastCommand != 'B')
	    break ;
	  
	  if (look->outfile)
	    look->outLevel = freeOutSetFile (look->outfile) ;

	  if (look->beginTable && ace_lower(spread->style) == 'x')
	    {
	      freeOut ("<html>\n<body bgcolor=white>\n<h2>") ;
	      if (*spread->titleBuffer) freeOutf ("<h2>\n%s\n</h2>", spread->titleBuffer) ;
	      spreadDumpLegend (spread, 'x') ;
	      freeOut ("<p>\n<table border=1>\n") ;
	    }

	  if (last && spread && kA)
	    { 
	      if (!spread->table)
		{ 
		  last = spreadRecomputeKeySet (spread, kA, last, look->minN, look->maxN) ;
	          
		  look->cumulatedTableLength += spreadDoDump (spread, '\t', spread->style, 
							      look->beginTable) ;
		  look->beginTable = FALSE ;
		}
	      else
		{ 
	          look->cumulatedTableLength += 
		    tablePartOut (&last, spread->table, '\t', spread->style) ;
		}
	      spread->showTitle = FALSE ;
	    } 
	  if (look->maxN >=0 && !spread->precomputing &&
	       look->cumulatedTableLength >= look->maxN) /* stop dumping */
	    { last = 0 ; spread->precomputing = FALSE ; /* cannot store ! */ }
	  if (!last)
	    { 
	      if (ace_lower(spread->style) == 'x')
		freeOut ("</table>\n</body>\n</html>\n") ;
	      else
		freeOutf ("\n//# %d lines in this table\n", 
			  look->cumulatedTableLength) ;
	      UNDO ;
	      ksNew = spread->exportedKeySet ; spread->exportedKeySet = 0 ;
              if (spread)
		spreadDestroy(spread);
	    }
	  else
	    { /* clear the DNA and other long arrays already exported */
	      spreadCleanUp (spread);
	    }
	  look->spread = spread ;
	  look->lastSpread = last ;
	  option = last ? 'B' : 'T' ; /* correct option for automatic looping */
	  if (look->outfile)
	    {
	      freeOutClose (look->outLevel--) ;
	      if (!last)
		{ filclose (look->outfile) ; look->outfile = 0 ; }
	    }
	}
	break ;
	
	/******************* Server        *****************************/
	
      case 'U':  /* shUtdown [now] */
      case 'W':  /* who */
	break ;
      case 'X':  /* data release */
	freeOutf("version %d = data_version, see wspec/server.wrm\n", 
		 getDataVersion ()) ;
	break ;
	
	/******************* Editors       *****************************/
      case 'O':  /* Comment */
	cp = freepos() ;
	if (cp) messdump("Comment %s", cp) ;
	freeOutf ("// Comment %s\n", cp ? cp : "") ;
        break ;
      case 'E':  /* EEdit */
	CHECK ;
	parseKeepGoing = TRUE ;
      case 'e': /* edit */
	CHECK ;
	freenext () ;
	cp = freepos() ;
	if (!cp || !*cp)
	  { freeOut ("// Usage Edit ace_command e.g. \n"
		     "// Edit Author Tom, adds author tom to the whole active set\n") ;
	    parseKeepGoing = FALSE ;
	    break ;
	  }
	cr = strnew (cp, 0) ;
	CHECK_WRITE ;
	i = freesettext (cr,"") ; tag = 0 ;
	if (freecard(i) && (cp = freeword()) && !strcmp (cp, "-D") &&
	    (cp = freeword()))
	  tag = str2tag (cp) ;
	freeclose (i) ;
	i = keySetMax (ksNew) ;
	localStack = stackReCreate (localStack, 40*i) ;
	kp = arrp (ksNew, 0, KEY) - 1 ;
	while (kp++, i--)
	  if (pickType(*kp) == 'B')  /* don t edit array this way */
	    { 
	      if (tag && !bIndexFind (*kp, tag)) /* no deletion needed */
		  continue ;
	      catText (localStack, className (*kp)) ;
	      catText (localStack, " ") ;
	      catText (localStack, freeprotect(name (*kp))) ;
	      catText (localStack, "\n") ;
	      catText (localStack, cr) ;
	      catText (localStack, "\n\n") ;
	    }
	ksTmp = keySetReCreate (ksTmp) ; /* do not UNDO, keep same active set */
	parseBuffer (stackText (localStack,0), ksTmp) ;
	keySetDestroy (ksTmp) ;
	stackDestroy (localStack) ;  /* it can be quite big */
	j = keySetMax(ksNew) ;
	freeOutf ("// I updated %d objects with command %s\n", j, cr) ;
	parseKeepGoing = FALSE ;
	messfree (cr) ;
	/* absolutely no reason to save here, mieg march 2000 */
	break ;
	
      case 'k':   /* Kill objects */
	CHECK ;
	CHECK_WRITE ;
	if (localIsInteractive &&
	    !messQuery (
			messprintf ("// Do you really want to destroy these %d objects", 
				    keySetMax (ksNew))))
	  break ;
	keySetKill (ksNew) ;
	ksNew = keySetReCreate (ksNew) ; /* no possible undo */
	break ;
	
      case 'i':   /* IS : keeps names matching template */
	{
	  int clientMax = 0 ;
	  
	  /* CHECK ; never CHECK above UNDO, this make sundo command unpredictable in scripts */
	  UNDO ;
	  
	  if ((cp = freeword())) 
	    {
	      if (!strcmp (cp, "-max"))
		{
		  if (!freeint (&clientMax) || clientMax <= 0)
		    {
		      freeOut ("// Usage IS   [-max nn]  [ < <= > >= ~ ] template : keeps names from current list matching template\n"
			       " // please give a positive integer after \"-max\" or avoid this parameter") ;
		      break ; 
		    }
		}
	      else
		freeback() ; 
	    }
	  freenext() ;
	  cp = freepos () ;
	  if (cp && *cp)
	    { 
	      cr = strnew(cp, 0) ;
	      ksNew = query (ksOld, messprintf ("IS %s", cr)) ;
	      messfree (cr) ;
	    }
	  else	if (clientMax)
	    ksNew = keySetCopy (ksOld) ;
	  else
	    freeOut ("// Usage IS   [-max nn]  [ < <= > >= ~ ] template : keeps names from current list matching template") ;
	  if (clientMax > 0 && clientMax < keySetMax (ksNew))
	    {
	      keySetDestroy(kA) ;
	      if (ksNew) 
		kA = keySetAlphaHeap(ksNew, clientMax) ;
	      keySetDestroy (ksNew) ;
	      ksNew = keySetCopy (kA) ;
	      keySetSort (ksNew) ;
	    }
	}
  break ;
	
      case 'r':   /* remove names matching template */
		  /* CHECK ; never CHECK above UNDO, this make sundo command unpredictable in scripts */
	UNDO ;
	if ((cp = freeword()))
	  { cr = strnew(cp, 0) ;
	    keySetDestroy (ksTmp) ;
	    ksTmp = query (ksOld, messprintf ("IS %s", cr)) ;
	    ksNew = keySetMINUS (ksOld, ksTmp) ;
	    keySetDestroy (ksTmp) ;
	    messfree (cr) ;
	  }
	else
	  freeOut ("// Usage Remove  [ < <= > >= ~ ] template : remove names matching template") ;
	break ;

      case 't':   /* Tag follow */
		  /* CHECK ; never CHECK above UNDO, this make sundo command unpredictable in scripts */
	UNDO ;
	tag = 0 ;
	if ((cp = freeword()))
	  { /* Autocomplete to a tag */
	    while (lexNext(0, &tag))
	      if(pickMatch(name(tag), cp))
		{ if (strcasecmp(cp, name(tag)))
		    freeOutf ("// Autocompleting %s to %s",
			      cp, name(tag)) ;
		  goto goodTag ;
		}
	    freeOut ("// Usage Follow tag-name, change to list of objects pointed by tag \n") ;
	    freeOutf ("// Sorry, Tag %s unknown \n", cp) ;
	    freeOut ("// type:  model class to see the model of the class\n\n") ;
	    break ;
	  }
	else
	  { freeOut ("// Follow without tag-name, enters active keySets") ;
	  }
      goodTag:
	ksTmp = keySetReCreate (ksTmp) ;
	for (i = 0, j = 0 ; i < keySetMax(ksOld) ; i++)
	  if (class(keySet(ksOld,i)) == _VKeySet) 
	    { KEY k = keySet (ksOld, i) ;
	      int i1 = 0 ;
	      KEYSET ks = arrayGet (k, KEY, "k") ; 
	      
	      if (ks)
		{
		  for (i1 = 0 ; i1 < keySetMax(ks) ; i1++)
		    keySet(ksTmp, j++) = keySet(ks, i1) ;
		}
	      keySetDestroy (ks) ;
	    }
	ksNew = tag ? query (ksOld, messprintf(">%s", name(tag))) : arrayHandleCreate (50, KEY, look->h) ;
	if (j)
	  { KEYSET ks = keySetOR (ksNew, ksTmp) ;
	    keySetDestroy (ksNew) ;
	    ksNew = ks ;
	    keySetSort (ksNew) ;
	    keySetCompress (ksNew) ;
	  }
	keySetDestroy (ksTmp) ;
	freeOutf ("\n// Found %d objects\n", keySetMax(ksNew)) ;
	break ;
	
      case 'b':
		  /* CHECK ; never CHECK above UNDO, this make sundo command unpredictable in scripts */
	UNDO ;
	i = 0 ;
	if ((cp = freeword()))
	  { if (!strcmp(cp,"-abstracts"))
	      i = 1 ;
	    else
	      freeback() ;
	  }
	kA = biblioFollow (ksOld) ;
	ksNew = keySetCopy (kA) ; keySetSort (ksNew) ; /* ! backwards */
	biblioDump (kA, i, 80) ;
	break ;

/******************* Dumpers       *****************************/

      case 'l':   /* lists names matching template */
	CHECK ;
	if (look->outfile)
	  filclose (look->outfile) ;
	look->outfile = 0 ;
	look->minN = 0 ;
	look->maxN = -1 ;
	look->notListable = 0 ;
	look->beauty = 'h' ;   
	cp = freeword();
	cr = 0 ;		/* used to hold possible name template */
	while(cp){
	  if (!strncmp(cp, "-a", 2))  /* ace format required */
	    look->beauty = 'a' ;
	  else if (!strcmp(cp, "-p"))  /* perl format required */
	    look->beauty = 'p' ;
	  else if (!strcmp(cp, "-perl"))  /* perl format required */
	    look->beauty = 'p' ;
	  else if (!strncmp(cp, "-j", 2))  /* java format required */
	    look->beauty = 'j' ;
	  else if (!strncmp(cp, "-J", 2))  /* java format required */
	    look->beauty = 'J' ;
	  else if (!strncmp(cp, "-C", 2))  /* ace c format required */
	    look->beauty = 'C' ;
	  else if (!strncmp(cp, "-h", 2))  /* human format explicitly required */
	    look->beauty = 'H' ;
	  else if (!strncmp(cp, "-b", 2))  /* minN required */
	    { if (freeint (&i) && i >= 0)
	      look->minN = i ;
	    }
	  else if (!strncmp(cp, "-c", 2))  /* maxN required */
	    { if (freeint (&i) && i >= 0)
	      look->maxN = i ;
	    }
	  else if (!strncmp(cp, "-f", 2))  /* -f filename */
	  {
	    if (!(cp = freepath()))  /* filename */
	      {
		freeOutf ("// fopen %s failed\n",cp ? cp : "Missing Filename") ;
		break ;
	      }
	    else
	      {
		if (look->beauty == 'h')  /* change default */
		  look->beauty = 'a' ;
		if (!strcmp(cp, "-"))
		  goto toutpret ;
		if (!(look->outfile=fopen(cp,"w")))
		  {
		    freeOutf ("// fopen %s failed\n",cp ? cp : "Missing Filename") ;
		    break ;
		  }
	      }
	    look->outLevel = freeOutSetFile (look->outfile) ;
	    maxChar = 0 ; /* mieg sept 2000 do no interrupt */
	  } 
	  else if (!strncmp(cp, "-t", 2)) /* -t [template] */
	    { if ((cp = freeword()))
	        cr = strnew (cp, 0) ;
	    } 
	  else			/* assume a tag */
	    cr = strnew (cp, 0) ;
	  freenext () ;
	  cp = freeword(); 
	}

      toutpret:
	keySetDestroy(kA) ;
	if (cr)
	  { ksTmp = arrayCreate (keySetMax (ksNew), KEY) ;
	    for (i = 0, j = 0 ; i < keySetMax(ksNew) ; i++)
	      if (pickMatch(name(keySet(ksNew,i)), cr))
		keySet(ksTmp, j++) = keySet(ksNew, i) ;
	    kA = keySetAlphaHeap(ksTmp, keySetMax(ksTmp)) ;
	    messfree (cr) ;
	    keySetDestroy (ksTmp) ;
	  }
	else
	  kA = keySetAlphaHeap(ksNew, keySetMax(ksNew)) ;
	look->nextObj = look->minN; 
	if (look->maxN >= 0) look->maxN += look->minN ;
	look->lastCommand = 'm'; 
	
	if (look->beauty == 'C')
	  {
	    char cc = 'c' ;
	    int nn = keySetMax(kA) ;

	    freeOutBinary (&cc, 1) ;
	    freeOutBinary ((char *)&nn, 4) ;
	  }
	else
	  freeOutf ("\nKeySet : Answer_%d\n",++nn) ;

	/* FALL THROUGH to more (list) */
	
      case 'm':
	if (look->lastCommand != 'm')
	  break ;
	/* stolen from keySetDump so can count char */
	{ int oldc = -1, newc, ii = 0 ;

	  i = look->nextObj ; kp  = arrp (kA,i,KEY) ;
	  if (i) oldc = class(*kp) ;
	  ii = freeOutByte() + maxChar ;
	  for (;i < keySetMax(kA) ; i++, kp++) 
	    { if (maxChar && freeOutByte () > ii)
		break;
	      if (look->maxN >=0 && i >= look->maxN)
		break ;
	      if (look->beauty != 'C' && !lexIsKeyVisible (*kp))
		{ look->notListable++ ; continue ; }
	      
	      switch (look->beauty)
		{
		case 'a':
		  cp = messprintf("%s %s\n", className(*kp), 
				  freeprotect(name(*kp))) ;
		  break ;
		case 'p':
		  cp = messprintf("%s : %s\n", className(*kp), 
				  freeprotect(name(*kp))) ;
		  break ;
		case 'C':
		  {
		    char buf0[2] = {0, '\n'} ;
		    char buf [3] ;
		    
		    if (i == look->minN) 
		      buf[0] = '>' ;
		    else 
		      buf[0] = '.' ;
		    if (iskey (*kp) == 2)
		      buf[1] = 'K' ;
		    else
		      buf[1] = 'k' ;
		    buf[2] = class (*kp) ;
		    freeOutBinary (buf, 3) ;
		    cp = name (*kp) ;
		    freeOutBinary (cp, strlen(cp)) ;
		    freeOutBinary (buf0, 2) ;
		    cp = 0 ; /* do not reprint it later ! */
		  }
                  break ;
		case 'j':
		  cp = messprintf("?%s?%s?\n",className(*kp),
				  freejavaprotect(name(*kp))) ;
                  break ;
		case 'J':
		  newc = class (*kp) ;
		  if (oldc != newc)
		    { oldc = newc ;
		      freeOutf ("?%s", className (*kp)) ;
		    }
		  else
		    freeOutf ("#", className (*kp)) ;
		  if (iskey (*kp) == 2)
		    cp = messprintf("?%s?\n",
				    freejavaprotect(name(*kp))) ;
		  else
		    cp = messprintf("!%s?\n",
				    freejavaprotect(name(*kp))) ;
                  break ;
	        default:
		  newc = class (*kp) ;
		  if (oldc != newc)
		    { oldc = newc ;
		      freeOutf ("%s:\n", className (*kp)) ;
		    }
		  cp = messprintf (" %s\n", name (*kp)) ;
		  break ;
		}
	      if (cp) freeOut (cp) ; /* fprintf over interprets % and \ */
	    }  
	  
	  if (i < keySetMax(kA) && (look->maxN <0 || i < look->maxN)) 
	    { look->nextObj = i; /* start here on more */
	      option = 'm' ;
/*	      if (look->outfile)
		freeOutf ("// %d object to listed \n", i);
 */
	    }
	  else
	    { 
	      option = 'l' ;
	    
	      if (look->beauty == 'C')
		freeOutBinary ("#\n", 2) ;
	      if (look->outfile)
		{ freeOut ("\n\n") ;
		  freeOutClose (look->outLevel--) ;
		  filclose (look->outfile) ;
		  look->outfile = 0 ;
		}
	      if (look->notListable)
		freeOutf("// %d object listed, %d not listable\n", 
			 i - look->minN - look->notListable, look->notListable) ;
	      else
		freeOutf("// %d object listed\n", i - look->minN) ;
	    }
	}
	break ;
	
      case 'w': /* wrtie with no argument is synonymous to show -a
		*/
	cp = freepos () ; /* optional filename */
	if (cp && *cp)
	  freeforcecard (messprintf ("-a %s ", cp)) ;
	else
	  freeforcecard ("-a ") ;
	/* fall through setting "-a -f <rest of old line>" */

      case 's':   /* Show [-p | -j | -a | -h | -C] [-f filename] [tag]*/
	/*
	* must look for -ml-test before it breaks out on empty keyset
	*/
	cp = freeword();
	if (cp && (strcmp(cp, "-ml-test") == 0))
	  {
	  freeOut("-ml-present\n");
	  break;
	  }
	CHECK ;
	keySetDestroy(kA) ;
	if (ksNew && keySetMax(ksNew))
	  kA = keySetAlphaHeap(ksNew, keySetMax(ksNew)) ;
	if (!kA)
	  { freeOut ("// Active list is empty\n") ;
	    break ;
	  }
	nn = 0 ;  /* error code in command line */
	look->minN = 0 ;
	look->maxN = -1 ;
	look->beauty = 'h' ;
	look->dumpTimeStamps = FALSE ;
	look->dumpComments = TRUE ;
	if (look->showCond)
	  { condDestroy (look->showCond) ;
	    look->showCond = 0 ;
	  }

	/*
	* first word was loaded into cp above
	*/

	while (cp) 
	  {
	    if (!strcmp (cp, "-p"))
	      look->beauty = 'p' ;
	    else if (!strcmp (cp, "-perl"))
	      look->beauty = 'p' ;
	    else if (!strcmp (cp, "-j"))
	      look->beauty = 'j' ;
	    else if (!strcmp (cp, "-J"))
	      look->beauty = 'J' ;
	    else if (!strcmp (cp, "-a"))
	      look->beauty = 'a' ;
	    else if (!strcmp (cp, "-C"))
	      look->beauty = 'C' ;
	    else if (!strcmp(cp, "-C"))  /* comments */
	      look->dumpComments = FALSE ; 
	    else if (!strcmp (cp, "-h"))
	      look->beauty = 'H' ;
	    else if (!strcmp (cp, "-ml"))
	      look->beauty = 'm';
	    else if (!strcmp(cp, "-T"))  /* time stamps */
	      look->dumpTimeStamps = TRUE ;
	    else if (!strcmp(cp, "-b"))  /* maxN required */
	      { if (freeint (&i) && i >= 0)
		  look->minN = i ;
		else
		  { freeOut ("// Error: integer required after -b\n") ;
		    nn = 1 ; 
		    break ;
		  }
	      }
	    else if (!strcmp(cp, "-c"))  /* maxN required */
	      { if (freeint (&i) && i >= -1)  /* -1 to get the whole list, used by aceperl */
		  look->maxN = i ;
		else
		  { freeOut ("// Error: integer required after -c\n") ;
		    nn = 1 ;
		    break ;
		  }
	      }
	    else if (!strcmp(cp, "-t")) /* -t [tag] */
	      { if ((cp = freeword()))
		  { char *posKeep = strnew (freepos(), 0) ;
	/* must cache rest of line since condConstruct uses freesubs */
		    if (!condConstruct (cp, &look->showCond))
		      { 
			nn = 1 ; 
			freeOutf ("// Unknown tag -t %s \n", cp) ;
			keySetMax(kA) = 0 ;
		      }
		    freeforcecard (posKeep) ;
		    messfree (posKeep) ;
		  }
		else
		  { freeOut ("// Error: tag name required after -t\n") ;
		    nn = 1 ;
		    break ;
		  }
	      } 
	    else if (!strcmp(cp, "-f"))  /* -f filename */
	      {	if (!(cp = freepath ())) /* filename */
		  { freeOut ("// Error: file name required after -f\n") ;
		    nn = 1 ;
		    break ;
		  }
		if (strcmp (cp, "stdout") && !(look->outfile = fopen(cp,"w")))
		  { freeOutf ("// Error: file open failed for %s\n",cp) ;
		    nn = 1 ;
		    break ;
		  }
	      }
	    else              /* assume a tag */
	      { 
		char *posKeep;

		if (cp[0] == '-')
		  { freeOutf ("// Error: invalid option %s \n"
			      " // valid options are: -a -h -p -j -J -ac -C -b -c -t -f\n",cp);
		    nn = 1 ;
		    break;
		  }

		posKeep = strnew (freepos(), 0) ;
		if (!condConstruct (cp, &look->showCond))
		  {
		    freeOutf ("// Error: invalid condition  %s \n", cp) ;
		    keySetMax(kA) = 0 ;
		    nn = 1 ;
		  }
		freeforcecard (posKeep) ;
		messfree (posKeep) ;
	      }

	    cp = freeword(); 
	  }
	if (0) freeOutf ("// Ready to predump %d objects\n", keySetMax(kA)) ;
	if (nn) /* error in command line */
	  break ;
	look->lastCommand = 'M' ;
	look->nextObj = look->minN ;
	if (look->maxN >= 0) look->maxN += look->minN ;
  /* Fall thru to case m */

      case 'M':
	if (look->lastCommand != 'M')
	  break ;
	  
	if (look->outfile)
	  look->outLevel = freeOutSetFile (look->outfile) ;

	freeOut ("\n") ;

	classe = 0 ;
	if (0) freeOutf ("// Ready to dump %d objects\n", keySetMax(kA)) ;
	ii = freeOutByte ()  + maxChar ;
	dumpTimeStamps = look->dumpTimeStamps ;
	dumpComments = look->dumpComments ;
	for (j = 0 , i = look->nextObj ; i < keySetMax(kA) ; ++i) 
	  { if (maxChar && freeOutByte () > ii)
	      break;
	    if (look->maxN >=0 && i >= look->maxN)
	      break ;
	    if (look->showStatus && ! (j%1000))
	      { int ll = freeOutSetFile (stderr) ;
		freeOutf ("// %d: %s\n", i, name(keySet(kA,i))) ;
		tStatus () ;
		freeOutClose (ll) ;
	      }
	    if (0) freeOutf ("// ....  dump %s\n", name (keySet (kA, i))) ;
	    dumpKeyBeautifully (keySet(kA,i), look->beauty, look->showCond) ;
	  }
	dumpTimeStamps = FALSE ;
	dumpComments = TRUE ;
	
	/* mieg sept 2000, close and reopen the file */ 
	if (look->outfile) 
	  {
	    freeOutClose (look->outLevel--) ; 
	  }

	if (i < keySetMax(kA) && (look->maxN <0 || i < look->maxN)) 
	  { option = 'M' ;
	    look->nextObj = i; 
/* 	    if (look->outfile)
	      freeOutf ("// encore %d object to come\n", keySetMax (kA) - i); 
*/
	  }
	else
	  { option = 's' ;
	    if (look->outfile)
	      { 
		filclose (look->outfile) ;
		look->outfile = 0 ;
	      }
	    freeOutf ("// %d object dumped\n", i - look->minN) ;
	  }
	break ;

      case 'd':	   /* DNA filename */
      case 'v':	   /* Peptide filename */
	CHECK ;
	keySetDestroy(kA) ;
	if (ksNew && keySetMax(ksNew))
	  kA = keySetAlphaHeap(ksNew, keySetMax(ksNew)) ;
	if (!kA)
	  { freeOut ("// Active list is empty\n") ;
	    break ;
	  }
	look->nextObj = 0;
	look->beauty = 'a' ;
	look->allowMismatches = FALSE ;
	look->noClassName = FALSE ;
	look->minN = look->maxN = 0 ;

	while ((cp = freeword()))
	  { 
	    if (*cp != '-')
	      { freeback () ; break ; }
	    if (!strcmp (cp, "-p"))
	      look->beauty = 'p' ;
	    else if (!strcmp (cp, "-perl"))
	      look->beauty = 'p' ;
	    else if (!strcmp (cp, "-C"))
	      look->beauty = 'C' ;
	    else if (!strcmp (cp, "-B"))
	      look->beauty = 'B' ;
	    else if (!strcmp (cp, "-mismatch"))
	      look->allowMismatches = TRUE ;
	    else if (!strcmp (cp, "-noClassName"))
	      look->noClassName = TRUE ;
	    else if (!strcmp (cp, "-x1"))
	      {
		if (!freeint (&look->minN))
		  {
		    if ((cp = freeword()))
		      {
			if (!strcasecmp (cp, "begin"))
			  look->minN = 1 ;
			else if (!strcasecmp (cp, "end"))
			  look->minN = -1 ;
			else
			  {
			    freeOutf ("-x1 should be \"begin\" OR \"end\" OR a  positive number") ;
			    look->minN = -2 ;
			    break ;
			  }
		      }
		  }
		else
		  {
		    if (look->minN <= 0)
		      {
			 freeOutf ("-x2 should be \"begin\" OR \"end\" OR a  positive number") ;
			look->minN = -2 ;
			break ;
		      }
		  }
	      }
	    else if (!strcmp (cp, "-x2"))
	      {
		if (!freeint (&look->maxN))
		  {
		    if ((cp = freeword()))
		      {
			if (!strcasecmp (cp, "begin"))
			  look->maxN = 1 ;
			else if (!strcasecmp (cp, "end"))
			  look->maxN = -1 ;
			else
			  {
			    freeOutf ("-x1 should be begin | end | positive number") ;
			    look->maxN = -2 ;
			    break ;
			  }
		      }
		  }
		else
		  {
		    if (look->maxN <= 0)
		       {
			 freeOutf ("-x1 should be begin | end | positive number") ;
			 look->maxN = -2 ;
			 break ;
		       }
		  }
	      }
	    else if (!strcmp (cp, "-f")) {} ; /* -f option is optional */
	      
	  }
	if (look->minN == -2 || look->maxN == -2)
	  break ;
	if (look->minN && look->minN == look->maxN)
	  {
	    freeOutf ("-x1 -x2 are equal, i cannot determine orientation") ;
	    break ;
	  }
	if ((look->minN && !look->maxN) || (!look->minN && look->maxN))
	  {
	    freeOutf ("-x1 -x2 must be both specified or absent");
	    break ;
	  }
	if (look->outfile) 
	  {
	    fclose (look->outfile) ; 
	    look->outfile = 0 ;
	  }

	freenext () ;
	if ((cp = freepath()))  /* -f filename */
	  { 
	    if (!strncmp(cp, "-f", 2))  /* -f (optional) */
	      { freenext () ; cp = freepath() ; }
	    if (cp)
	      {	if (!(look->outfile = fopen(cp,"w")))
		  { freeOutf ("// fopen %s failed\n",cp ? cp : "Missing Filename") ;
		    break ;
		  }
	      }
	  }

	look->fastaCommand = option ;
	look->lastCommand = 'D' ;
  /* Fall thru to case D */
	
      case 'D':
	if (look->lastCommand != 'D')
	  break ;
	
	if (look->outfile)
	  look->outLevel = freeOutSetFile (look->outfile) ;
	ii = freeOutByte ()  + maxChar ;
	for (i = look->nextObj ; i < keySetMax(kA) ; ++i) 
	  {
	    if (maxChar && freeOutByte () > ii)
	      break;
	    
	    key = keySet (kA, i) ;
	  
	    if (look->beauty == 'p')
	      {
		if (look->fastaCommand == 'd')
		  { 
		    Array dna = 0 ;
		    if (!(dna = dnaGet (key)))
		      continue ;
		    dnaDecodeArray (dna) ;
		    array(dna, arrayMax(dna), char) = 0 ; /* safeguard */
		    if (!strlen(arrp(dna, 0, char)))
		      {
			arrayDestroy (dna) ;
			continue ;
		      }
		    freeOutf ("{ty=>dna,va=>%s,DNA=>", name(key)) ;
		    freeOut (arrp(dna, 0, char)) ;
		    freeOut ("}\n") ;
		    arrayDestroy (dna) ;
		  }
		else
		  { 
		    Array pep = 0 ;
		    if (!(pep = peptideGet (key)))
		      continue ;
		    pepDecodeArray (pep) ;
		    array(pep, arrayMax(pep), char) = 0 ; /* safeguard */
		    if (!strlen(arrp(pep, 0, char)))
		      {
			arrayDestroy (pep) ;
			continue ;
		      }
		    freeOutf ("{ty=>pep,va=>%s,PEPTIDE=>", name(key)) ;
		    freeOut (arrp(pep, 0, char)) ;
		    freeOut ("}\n") ;
		    arrayDestroy (pep) ;
		  }
	      }
	    else
	      if (look->fastaCommand == 'd')
		{
		  dnaZoneDumpFastAKey (key, 0, 0, look->beauty, look->minN, look->maxN, look->allowMismatches, look->noClassName) ;
		}
	      else
		pepDumpFastAKey (key, 0, 0, look->beauty) ;
	  }
	if (look->outfile)
	  freeOutClose (look->outLevel--) ;
	if (i < keySetMax(kA))
	  { option = 'D' ;
	    look->nextObj = i; 
	  }
	else
	  { option = look->fastaCommand ;
	    look->fastaCommand = 0 ;
	    if (look->outfile)
	      { 
	        filclose (look->outfile) ;
	        look->outfile = 0 ;
	      }
	    freeOutf ("// %d object dumped", keySetMax (kA));
	  }
	break ;

      case 'a':
	CHECK ;
	/*
	   if (lexClassKey (freeword(), &key)
	   && (fmt = freeword())
	   && (a = uArrayGet(key, sizeOfFormattedArray(fmt), fmt)))
	   { dumpFormattedArray (.. should be freeOut...myFile, ...myStack, a, fmt) ;
	   arrayDestroy (a) ;
	   }
	   */
	break ;
	
      case 'L':		/* dump alignments in fancy format */
	CHECK ;
	cp = freepath() ;
	if (!cp || !*cp || !(f = fopen (cp, "w")))
	  f = stdout ;
	keySetDestroy(kA) ;
	kA = keySetAlphaHeap(ksNew, keySetMax(ksNew)) ;
	alignDumpKeySet (kA, f) ;
	if ((f != stdout))
	  fclose (f) ; 
	f = 0 ;
	break ;

      case '-':			/* depad an assembly */
	cp = freeword() ;
	cr = cp ? strnew(cp, 0) : 0 ;
	CHECK_WRITE ;
	{ extern BOOL depad (char*) ;
	  depad (cr) ;
	}
	break ;

      case '*':			/* pad an assembly */
	cp = freeword() ;
	cr = cp ? strnew(cp, 0) : 0 ;
	CHECK_WRITE ;
	{ extern BOOL pad (char*) ;
	  pad (cr) ;
	}
	break ;

      case 400:			/* export wspec */
        freeOut("\n// wspec directory exported\n") ;

	break ;
      case 401:			/* export wspec */
	if ((f = filopen("wspec/layout","wrm","r")))
	  { 
	    i = freesetfile (f,"") ;
	    while (freecard(i))
	      freeOutf("%s\n", freepos()) ;
	    freeOut("\n// end of layout file\n") ;
	  }
        else
	  freeOut("\n// Sorry, no file wspec/layout.wrm available \n") ;
	break ;
      case 402:			/* system_list */
        {
	  KEY tag, cl ; 
	 
	  char buf0[2] = {0, '\n'} ;
	  char buft [2] ;
	  unsigned int ii ;
	  char *cp ;

	  ii = 0x12345678 ;
	  freeOutBinary ((char*)&ii,4) ;
	  ii = lexMax(0) ;
	  freeOutBinary ((char*)&ii,4) ;

	  buft[0] = '>' ;
	  buft[1] = 'g' ;
	  for (tag = 0 ; tag < lexMax(0) ; tag++)
	    { 
	      freeOutBinary (buft,2) ;
	      cp = name(tag) ;
	      freeOutBinary (cp, strlen(cp)) ;
	      freeOutBinary (buf0,2) ;
	      buft[0] = '.' ;
	    }
	  buft[1] = 'k' ;
	  for (cl = 0 ; cl < 256 ; cl++)
	    {  
	     
	      freeOutBinary (buft,2) ;
	      if (pickList[cl].name)
		{
		  cp = pickClass2Word (cl) ;
		  freeOutBinary (cp, strlen(cp)) ;
		}
	      freeOutBinary (buf0,2) ;
	    } 
	  freeOutBinary ("#\n",2) ;
	}

	break ;
      case '#':			/* GIF drawing access */
	if (gifEntry)
	  { 
	    if (freecheck ("w")) /* commands on one line separated by ';' */
	      { cr = strnew(freepos(), 0) ;
 	        i = freesettext (cr, 0) ;
		freespecial ("\n\"\\/%;") ; /* recognize ';' as separator */
		(*gifEntry)(ksNew, i, FALSE) ;
		messfree (cr) ;
	      }
	    else if (localIsInteractive > 0 &&  /* accept multiple lines */
		     !(look->choix & 8)) /* server should never open  stdin */
	      { i = freesetfile (stdin, 0) ;  
		(*gifEntry)(ksNew, i, TRUE) ;
	      }  
	    else 
	      { freeOut ("// Please provide the full gif command on a single line\n") ;
 	        i = freesettext (" ? ", 0) ;
		(*gifEntry)(ksNew, i, FALSE) ;
	      }
	    freeclose (i) ;
	  }
	break ;
	
	/******************* Stack operations  *****************************/
	
      case '1':  /* stack push */
        keySetDestroy (kA) ;
	kA = keySetCopy(ksNew) ;
	freeOut ("// Pushed the active set on stack\n") ; 
	goto finStack ;
	
      case '2': case '3': case '4': case '5': case '6': case '7': case '8':
	if (!nns)
	  { freeOut ("// Sorry, the stack is empty\n") ;
	    break ;
	  }
	    switch (option)
	      {
	      case '2':  /* stack pop */
		UNDO ;
		ksNew = pop(ksStack,KEYSET) ; nns-- ;
		if (nns) 
		  { kA = pop(ksStack,KEYSET) ; nns-- ;
		  }
		else
		  kA = 0 ;
		break ;

	      case '3':  /* stack swap */
		kA = ksNew ;
		ksNew = pop(ksStack,KEYSET) ; nns-- ;
		break ;

	      case '4': case '5': case '6': case '7':
		ksTmp = pop (ksStack,KEYSET) ; nns-- ;
		switch (option)
		  { 
		  case '4':  /* stack AND */
		    kA = keySetAND (ksTmp, ksNew) ;  
		    break ;
		  case '5':  /* stack OR */
		    kA = keySetOR (ksTmp, ksNew) ;  
		    break ;
		  case '6':  /* stack XOR */
		    kA = keySetXOR (ksTmp, ksNew) ;  
		    break ;
		  case '7':  /* stack MINUS */
		    kA = keySetMINUS (ksTmp, ksNew) ;  
		    break ;
		  }
		keySetDestroy (ksTmp) ; 
		break ;

	      }


	  finStack:
	    if (kA) 
	      { push(ksStack, kA, KEYSET) ; nns++ ;
		freeOutf ("// The stack now holds %d keyset(s), top of stack holds %d objects\n", 
			  nns, kA ? keySetMax (kA) : 0 ) ;
		kA = 0 ;
	      }
	    else
	      freeOut ("// The stack is now empty \n") ;
	break ;
  
      case '9': /* subset */
	{
	  int x0, nx, i, j ;
	  
	  if (!freeint (&x0) ||
	      !freeint (&nx) ||
	      x0 < 1 ||
	      nx < 1)
	    freeOutf ("// Usage: subset x0 nx, : returns the subset of the alphabetic active keyset from x0>=1 of length nx") ;
	  else
	    {	      
	      keySetDestroy(kA) ;
	      if (ksNew) 
		{
		  if (x0 < 1) 
		    x0 = 1 ;
		  if (nx > keySetMax (ksNew) - x0 + 1) 
		    nx = keySetMax (ksNew) - x0 + 1 ; 
		  kA = keySetAlphaHeap(ksNew, x0 + nx - 1) ;
		}
	      else 
		nx = 0 ;
	      UNDO ;
	      ksNew = keySetCreate () ;
	      for (i = x0 - 1, j = 0 ; nx > 0 ; i++, j++, nx--)
		keySet (ksNew, j) = keySet (kA, i) ;
	      keySetSort (ksNew) ;
	      keySetDestroy(kA) ;
	    }	  
	}
	break ;	

	/******************* named Stack operations  *****************************/
	
      case 301:  /* kstore name */
	cp = freeword() ;
	if (!cp)
	  freeOutf ("// usage: kstore name :  store active keyset in internal tmp space as name\n") ;
	else
	  { 
	    KEYSET nKs = arrayHandleCopy (ksNew, look->h) ;

	    if (!look->dict)
	      { 
		look->dict = dictHandleCreate (32, look->h) ;
		dictAdd (look->dict, "___toto", 0) ;
	      }
	    if (!look->namedSets)
	      look->namedSets = arrayHandleCreate (32, Array, look->h) ;
	    if (dictFind (look->dict, cp, &nn) &&
		nn < arrayMax (look->namedSets))
	      keySetDestroy (arr (look->namedSets, nn, Array)) ;
	    else
	      dictAdd (look->dict, cp, &nn) ;
	    array (look->namedSets, nn, Array) = nKs ;

	    freeOutf ("// Stored named set %s\n", cp) ;
	  }
	break ;

      case 302:  /* kget name */
	cp = freeword() ;
	if (!cp)
	  freeOutf ("// usage: kget name  : copy named set as active keyset\n") ;
	else
	  { 
	    int nn = 0 ;
	    KEYSET nKs = 0 ;

	    if (look->dict && dictFind (look->dict, cp, &nn) &&
		look->namedSets && 
		nn < arrayMax (look->namedSets) &&
		(nKs = arr (look->namedSets, nn, Array)) &&
		keySetExists (nKs))
	      {
		UNDO ;
		keySetDestroy (kA) ;
		ksNew = keySetCopy (nKs) ;
		freeOutf ("// Recovered named set %s\n", cp) ;
	      }
	    else
 	      {
	      freeOutf ("// empty set %s\n", cp) ;
	      UNDO ;
	      keySetDestroy (kA) ;
	      kA = keySetCreate();
 	      }
	  }
	break ;

      case 303:  /* kclear name */
	cp = freeword() ;
	if (!cp)
	  freeOutf ("// usage:  clear named set\n") ;
	else
	  {
	    int nn = 0 ;

	    if (look->dict && dictFind (look->dict, cp, &nn) &&
		look->namedSets && 
		nn < arrayMax (look->namedSets) &&
		(ksTmp = arr (look->namedSets, nn, Array)) &&
		keySetExists (ksTmp))
	      {
		keySetDestroy (ksTmp) ;
		arr (look->namedSets, nn, Array) = 0 ;
		freeOutf ("// Cleared named set %s\n", cp) ;
	      }
	    else
	      freeOutf ("// Sorry, unkown set %s\n", cp) ;
	  }
	break ;
	

         /******************* external calls *****************************/
	    
      case 'Y':	
        if (look->choix & 8) /* server should never open  stdin */
	  break ;
	f = 0 ;
	if ((cp = freeword())) /* -f filename */
	  { 
	    if (!strncmp(cp, "-f", 2))
	      {
		if (!(cp = freepath()) || /* filename */
		    !(f = fopen(cp,"r")))
		  { freeOutf ("// Sorry: I cannot find file \"%s\"\n",cp ? cp : "") ;
		    break ;
		  }
		else
		  freeOutf ("// Executing tace -f %s\n",cp) ;
	      }
	    else
	      freeback() ; 
	  }
	else 
	  {
	    freeOut ("// Usage tace -f commandfile parameters\n") ;
	    break ;
	  }
	freenext () ;  /* get parameters */
	cr = strnew(freepos(), 0) ;
	i = isInteractive ;
	commandExecute (j = freesetfile (f, cr), FALSE, FALSE, stdout, 0, 3, look->ksNew) ;
	isInteractive = i ;
	freeclose (j) ;
	messfree (cr) ;
	break ;

#ifdef ACEMBLY
      case 'A':	
        if (look->choix & 8) /* server should never open  stdin */
	  break ;
	f = stdin ;
	if ((cp = freeword())) /* -f filename */
	  { 
	    if (!strncmp(cp, "-f", 2))
	      {
		if (!(cp = freepath()) || /* filename */
		    !(f = fopen(cp,"r")))
		  { freeOutf ("Sorry: I cannot find file \"%s\"\n",cp ? cp : "") ;
		  break ;
		  }
	      }
	    else
	      freeback() ; 
	  }
	freenext () ;  /* get parameters */
	cr = strnew(freepos(), 0) ;
	CHECK_WRITE ;
	{ 
	  BOOL oldIsTS = isTimeStamps ;
	  BOOL oldIsIntr = isInteractive ;

	  isTimeStamps  = FALSE ;   /* to accelerate the code */
	  i = freesetfile (f, cr) ;
	  defComputeTace(i, ksNew) ;
	  freeclose (i) ;
	  isTimeStamps = oldIsTS ;
	  isInteractive = oldIsIntr ;
	}
	messfree (cr) ;
	break ;
#endif
	

	/****************************************************************/
#ifdef JUNK	
      case 450:			/* AQL */
	freenext() ;
	look->beauty = 'a' ;
	look->outfile = 0 ;
	while ((cp = freeword()))   /* get aql options */
	  {
	    if (*cp != '-')
	      { freeback () ; break ; }
	    else if (!strncmp (cp, "-s", 2))
	      { look->beauty = 0 ; }
	    else if (!strncmp (cp, "-a", 2))
	      { look->beauty = 'a' ; }
	    else if (!strncmp (cp, "-h", 2))
	      { look->beauty = 'h' ; }
	    else if (!strncmp (cp, "-C", 2))
	      { look->beauty = 'C' ; }
	    else if (!strncmp (cp, "-j", 2))
	      { look->beauty = 'j' ; }
	    else if (!strncmp (cp, "-J", 2))
	      { look->beauty = 'J' ; }
	    else if (!strncmp(cp, "-o", 2))
	      { 
		cp = freepath () ;
		if (!(look->outfile = filopen(cp,"w")))
		  { freeOutf ("// out_fileopen %s failed\n",cp ? cp : "Missing Filename") ;
		    break ;
		  }
	      }

	  }
	cp = freepos () ;         /* get aql query */
	if (!cp || !strlen(cp))
	  { 
	    freeOut ("// Aql: missing query") ;
	    ksNew = arrayHandleCreate (50, KEY, look->h) ;
	    break ;
	  }
	freeOut ("// New Acedb Query Language : AQL \n") ;  
	
	if (look->outfile)
	  look->outLevel = freeOutSetFile (look->outfile) ;

	switch (option)
	  {
	  case 450:
	    cq = strnew(cp,0) ;
	    break ;
	  case 451:  /* select */
	  case 452:  /* select */
	  case 453:  /* select */
	  case 454:  /* select */
	    while (*cp == ' ') cp++ ;
	    if (!strcasecmp (cp, "select "))
	      cp += 7 ;
	    cq = strnew(messprintf("select %s", cp),0) ;
	    break ;
	  }
	if (1)
	  {
	    AQL aql = aqlCreate2 (cq, 0, 0, 0, 0, 1, "@active", 0) ;
	    TABLE *activeTable = tableCreateFromKeySet(ksNew, 0);
	    TABLE *result = aqlExecute (aql, 0, 0, 0, "TABLE @active", activeTable, 0);
	    
	    tableDestroy (activeTable);
	    {
	      DICT *dict = dictCreate (100) ;
	      int iii ;
	      dictAdd (dict, "zorro", &iii) ;
	      dictFind(dict, "zorro", &iii) ;
	      dictDestroy (dict) ;
	    }
	    if (!(aqlIsError(aql)))
	      {
		if (look->beauty) tableOut (result, '\t', look->beauty) ;
		UNDO ;
		ksNew = keySetCreateFromTable (result, look->h) ;
	      }
	    else
	      {
		freeOutf ("\n// Aql: error %d : %s\n", 
			  aqlGetErrorNumber(aql), aqlGetErrorMessage(aql)) ;
	      } 
	    
	    tableDestroy (result) ;
	    aqlDestroy (aql) ;
	  }
      
	
	if (look->outfile)
	  {
	    freeOutClose (look->outLevel--) ;
	    filclose (look->outfile) ; look->outfile = 0 ;
	  }
	messfree (cq) ;
	break ;  /* AQL */
#endif
	
      case 451:  /* AQL alias transferred to BQL y(select ... ) */
      case 452: /* BQL */
      case 453: /* BQL */
      case 454: /* BQL */
      case 455: /* BQL */
	{
	  AC_HANDLE h = ac_new_handle () ;
	  ACEOUT ao = 0 ;
	  ACEIN ai = 0 ;
	  vTXT txt = 0 ;
	  Stack bqlStack = 0 ;
	  BOOL ok = TRUE ;
	  BOOL setOutFile = FALSE ;
	  BOOL debug = FALSE ;
	  BOOL silent = (option == 455 ? TRUE : FALSE) ;
	  BOOL showTitle = FALSE ;
	  char *fileTitle = 0 ;
	  BOOL acedbQuery = FALSE ;
	  int maxLine = 0 ;
	  int beginLine = 0 ;
	  freenext() ;
	  look->beauty = 'h' ;
	  look->outfile = 0 ;

	  while ((cp = freeword()))   /* get bql options */
	    {
	      if (*cp != '-')
		{ freeback () ; break ; }
	      else if (!strncmp (cp, "--title", 7))
		{ 
		  showTitle = TRUE ; 
		  cp = freeword () ;
		  if (cp)
		    fileTitle = strnew (cp, h) ;
		}
	      else if (!strncmp (cp, "-t", 2))
		{ showTitle = TRUE ; }
	      else if (!strncmp (cp, "-a", 2))
		{ look->beauty = 'a' ; }
	      else if (!strncmp (cp, "-h", 2))
		{ look->beauty = 'h' ; }
	      else if (!strncmp (cp, "-s", 2))
		{ silent = TRUE ; }
	      else if (!strncmp (cp, "-C", 2))
		{ look->beauty = 'C' ; }
	      else if (!strncmp (cp, "-j", 2))
		{ look->beauty = 'j' ; }
	      else if (!strncmp (cp, "-J", 2))
		{ look->beauty = 'J' ; }
	      else if (!strncmp (cp, "-q", 2))
		{ acedbQuery = TRUE ; }
	      else if (!strncmp (cp, "-v", 2))
		{ debug = TRUE ; }
	      else if (!strncmp (cp, "-c", 2))
		{ freeint (&maxLine) ; }
	      else if (!strncmp (cp, "-b", 2))
		{ freeint (&beginLine) ; }
	      else if (!strncmp(cp, "-o", 2))
		{ 
		  setOutFile = TRUE ;
		  ao = aceOutCreate (freepath (), 0, FALSE, h) ; 
		  if (ao && fileTitle)
		    aceOutDate (ao, "###", fileTitle) ;

		  if (! ao)
		    { 
		      freeOutf ("// bql -o %s failed, please edit the file name\n",cp ? cp : "Missing Filename") ;
		      ok = FALSE ;
		      break ;
		    }
		}
	      else if (!strncmp(cp, "-i", 2))
		{ 
		  ai = aceInCreate (freepath (), FALSE, h) ; 
		  if (! ai)
		    { 
		      freeOutf ("// bql -i %s failed, please edit the file name\n",cp ? cp : "Missing Filename") ;
		      ok = FALSE ;
		      break ;
		    }
		}
	      else 
		{ freeback () ; break ; }
	    }
	  if (look->outLevel && !setOutFile)
	    {
	      bqlStack = stackHandleCreate (100000, h) ;
	      ao = aceOutCreateToStack (bqlStack, h) ;     /* export to the freeOutf stack */
	    }
	  if (! ao)
	    ao = aceOutCreate (0, 0, FALSE, h) ;  /* export to stdout */
	  if (ai)
	    {
	      txt = vtxtHandleCreate (h) ;
	      while (aceInCard (ai))
		{
		  cp = aceInPos (ai) ;
		  if (cp)
		    vtxtPrintf (txt, " %s ", cp) ;
		}
	      cq = vtxtPtr (txt) ;
	    }
	  else
	    {
	      cp = freepos () ;         /* get bql query */	 
	      if (cp)
		{
		  while (*cp == ' ') cp++ ;
		  if (!strncasecmp (cp, "select ", 7))
		    cp += 7 ;
		  
		  cq = hprintf (h, "%s%s", acedbQuery ? "" : "select ", cp) ;
		}
	      else
		{
		  freeOutf ("// bql : no query was provided, try bql -help\n") ;
		  ok = FALSE ;
		}
	    }
	  
	  if (ok)
	    {
	      BQL *bql = 0 ;
	      int nn = 0 ;
	      ksOld = look->ksNew ;
	      ksNew = arrayHandleCreate (50, KEY, look->h) ;
	      
	      bql = bqlCreate (debug, h) ;
	      if (beginLine < 1)
		beginLine = 1 ;
	      if (0) /* to trim we need to work off keySetAlpha and have normal sort */
		bqlMaxLine (bql, maxLine) ;
	      if (cq)
		{
		  char *cr = cq + strlen(cq) - 1 ;
		  while (*cr == ' ' && cr > cq)
		    cr-- ;
		  if (*cr == ';')
		    silent = TRUE ;
		}
	      if (!bqlParse (bql, cq, acedbQuery))
		{
		  const char *ccp = bqlError (bql) ;
		  freeOutf ("// bql parse error\n%s\n", ccp ? ccp : "") ;
		  ok = FALSE ;
		}
	      else if (! (nn = bqlRun (bql, ksOld, ksNew)))
		{
		  if (0) freeOutf ("// bql run error\n%s\n", bqlError (bql)) ;
		  keySetDestroy (look->ksOld) ;
		  look->ksOld = look->ksNew ;
		  look->ksNew = keySetCreate () ;
		  if (! silent)
		    {
		      freeOutf ("// bql empty table\n%s\n", bqlError (bql)) ;
		      ok = FALSE ;
		    }
		}
	      else 
		{
		  keySetDestroy (look->ksOld) ;
		  look->ksOld = look->ksNew ;
		  look->ksNew = ksNew ;
		  if (beginLine < 1) beginLine = 1 ;  /* do not Plato, to have the same interface -b -c as list and show */
		  if (! silent) /* in silent mode, the active keyset is modified */
		    bqlExport (bql, look->beauty, showTitle, beginLine - 1, maxLine, ao) ;
		  if (bqlStack)
		    {
		      freeOut (stackText (bqlStack, 0)) ;
		      stackDestroy (bqlStack) ;
		    }
		  freeOutf ("// %d lines\n", nn) ;
		}	
	    }	

	  ac_free (h) ;
	}
	break ;
	
      default:
	freeOut ("// Sorry, option not written yet\n") ;
      }
    }
  else
    if (localIsInteractive)
      freeOut ("// Please type ? for a list of commands. \n") ;
  
  
  look->ksNew = ksNew ? ksNew : arrayHandleCreate (50, KEY, look->h) ; /* mieg jan 2000 */
  look->ksOld = ksOld ;
  look->kA = kA ;
  look->ksStack = ksStack ;
  look->nns = nns ;
  look->lastCommand = option ;
  
  stackDestroy (localStack) ;
  switch (option)
    {
    case 'B': case 'm': case 'M': case 'D': case 'V' :
      break ;
    default:
      if (! isView)
        freeOutf ("// %d Active Objects\n", 
		  keySetCountVisible(look->ksNew)) ;
      else if (isView == 2)
	isView = 0 ;
    } 
  if (0) printf("command ends, before restoring isInteractive=%d\n",  isInteractive) ;
  isInteractive = oldIsInteractive ; 
  if (0) printf("command ends, switched back to isInteractive=%d\n",  isInteractive) ;
  return 
    option ;
#if 0
 abort:
  if (0) printf("command abort, before restoring isInteractive=%d\n",  isInteractive) ;
  isInteractive = oldIsInteractive ; 
  if (0) printf("command abort, switched back to isInteractive=%d\n",  isInteractive) ;
  if (isWriteAccess()) /* so that 'echo "parse toto.ace" | tace .' saves implicitely */
    sessionDoSave (FALSE) ;
  return 'q' ;
#endif
}
  
void aceCommandSwitchActiveSet (AceCommand look, KEYSET ks, KEYSET ks2) 
{
  if (ks)
    {
      keySetDestroy (look->ksNew) ;
      look->ksNew = keySetCopy (ks) ;
    } 
  if (look->ksNew && ks2)
    {
      int i = keySetMax (look->ksNew) ;
      while (i--)
	keySet (ks2, i) = keySet (look->ksNew, i) ;
    }
  return ;
}

AceCommand aceCommandCreate (int choix, int toto)
{
  return aceCommandDoCreate (choix, 0, 0) ;
}
AceCommand aceCommandDoCreate (int choix, FILE *fil, Stack s)
{ 
  AceCommand look;

  look = (AceCommand) messalloc (sizeof (struct _AceCommandStruct)) ;

  look->magic = COMMAND_MAGIC ;
  look->choix = choix ;
  look->aMenu = choisirMenu (choix) ;

  look->h = handleCreate () ;
  look->ksNew = arrayHandleCreate (50, KEY, look->h) ;
  look->ksStack = stackHandleCreate(5, look->h) ;


  freeOutInit () ;
  if (fil)
    look->outLevel = freeOutSetFile (fil) ;
  else if (stackExists(s))
     {
       s->textOnly = TRUE ;
       look->outLevel = freeOutSetStack (s) ;
     }

  return look ;
} /* aceCommandCreate */

void aceCommandNoWrite (AceCommand look)
     /* drop write access for this look */
{ 
  if (look->magic != COMMAND_MAGIC)
    messcrash ("Bad magic in aceCommandPublic") ;
  look->noWrite = TRUE ;
}

void aceCommandDestroy (AceCommand look)
{ 
  KEYSET ks ;

  if (look->magic != COMMAND_MAGIC)
    messcrash ("Bad magic in aceCommandDestroy") ;
  if (look->outLevel)
    freeOutClose (look->outLevel) ;
  look->outfile = 0 ;
  arrayDestroy (look->aMenu) ;
  keySetDestroy (look->ksNew) ;
  keySetDestroy (look->ksOld) ;
  keySetDestroy (look->kA) ;
  while (!stackAtEnd (look->ksStack))
    { ks = pop (look->ksStack,KEYSET) ;
      if (keySetExists (ks))
	keySetDestroy (ks) ;
    }
  stackDestroy (look->ksStack) ;
  if (look->spread) 
    spreadDestroy (look->spread);
  if (look->showCond)
    condDestroy (look->showCond) ;

  handleDestroy (look->h) ;  /* deals with look->namedSets and look->dict */

  look->magic = 0 ;
  messfree (look) ;
}

void commandExecute (int level, BOOL noSubshell, 
		     BOOL localIsInteractive,
		     FILE *fil, Stack s, int choix,
		     KEYSET activeKs)
{ 
  AceCommand look;
  static int niveau = 0 ;
  int maxChar = 0 ;
  KEY option ;

  look = aceCommandDoCreate (choix,fil,s) ;
  if (noSubshell || !getenv("ACEDB_SUBSHELLS"))
    noSubshell = TRUE ;

  niveau++ ;
  if (localIsInteractive && !getenv("ACEDB_NO_BANNER"))
    freeOut ("\n\n// Type ? for a list of options\n\n") ;

  option = 0 ;
  look->ksNew = activeKs ? keySetCopy (activeKs) : keySetCreate() ;
  while (TRUE)
    {
      if (noSubshell) /* reset in loop if changed in the command */
	freespecial ("\n\t\"\\/@%") ;  /* forbid sub shells */
      else  
	freespecial ("\n\t\"\\/@%$") ; /* allow sub shells */
      
  
      option = aceCommandDoExecute (look, level,  
			    localIsInteractive, option, maxChar) ;
      switch (option)
	{
	case 'q':
	  goto fin ;
	case 'B':
        case 'D':
	case 'm':
	case 'M':  /* automatic looping */	 
	  break ; 
	default:
	  option = 0 ;
	  break ;
	}
    }

 fin:
  freeclose (level) ;
  aceCommandDestroy (look) ;
  if (!niveau--)
    sessionDoSave (FALSE) ; 
} /* commandExecute */

/* simply exe the command and return the result on the stack */
Stack commandStackExecute (AceCommand look, const char *command)
{
  int level_in = freesettext (command, 0) ; /* establish tace pseudo-stdin */
  int level_out;
  int option = 0 ;
  int maxChar = 0 ;

  Stack s = stackCreate (3000) ;

  freespecial ("\n\t\"\\/%") ;  /* forbid sub shells and includescommand_ct*/

  level_out = freeOutSetStack(s) ;
  while (TRUE)
    { 
      option = aceCommandDoExecute (look, level_in,
				    FALSE /* no prompt */,
				    option, maxChar) ;
      switch (option)
	{
	case 'q':
	  goto fin ;
	case 'B':
        case 'D':
	case 'm':
	case 'M':  /* automatic looping */	 
	  break ; 
	default:
	  option = 0 ;
	  break ;
	}
    }
 fin:

  freeclose(level_in);
  freeOutClose(level_out);
  
  return s ;
}

/****************************************************************/
/****************************************************************/

#include <wac/ac.h>
#include "aceio.h"
#include "sysclass.h"
static void taceTest (KEYSET ks)
{
  KEY key ;
  lexaddkey ("toto1", &key, _VBinary) ;

  if (!sessionGainWriteAccess()) 
    {
      freeOut ("// Sorry, you cannot gain Write Access\n") ; 
      return ;
    }
  
  if (0)
    {
      FILE *f = filopen("totobin","ace", "w") ;
      Array a = arrayCreate (1000, unsigned char) ;
      char *cp = "Hello worldHello worldHello worldHello worldHello worldHello worldHello worldHello worldHello worldHello worldHello worldHello worldHello worldHello worldNext\nHello worldHello worldHello worldHello worldHello worldHello worldHello worldHello worldHello worldHello worldHello worldHello worldFin" ;
      char *cq ;
      int n = strlen(cp) ; 

      arrayMax (a) = n = strlen(cp) ;
      cq = arrp (a, 0, char) ;
      memcpy (cq, cp, n) ;

      arrayStore (key, a, "c") ;
      dumpKey (key, f, 0) ;
      filclose (f) ;
    }
  if (1)
    {
      Array a = arrayGet (key, unsigned char, "c") ;
      unsigned char *cp ;
      
      if (a)
	{
	  cp = arrayp (a, 0, unsigned char) ;
	  printf ("%s\n", cp) ;
	  arrayDestroy (a) ;
	}
      else
	printf ("cannot find object Binary:toto1\n") ;
    }

  return ;
}


/****************************************************************/
/****************************************************************/
