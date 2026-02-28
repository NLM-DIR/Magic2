/*  File: gifacemain.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * SCCS: $Id: gifacemain.c,v 1.33 2020/05/30 16:50:31 mieg Exp $ 
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Apr 22 18:37 1999 (rd)
 * * Jan  8 10:34 1999 (edgrif): Fixed two bugs: 1) call to graphPS was
 *              completely incorrect, wrong params !!
 *              2) display option incorrectly freed keyset supplied to it.
 * * Dec  2 16:20 1998 (edgrif): Corrected decl. of main, added code
 *              to record build time of this module.
 * * Oct 15 11:30 1998 (edgrif): Add new graph/acedb initialisation call.
 * * Sep 15 10:29 1998 (edgrif): Correct faulty main declaration, add
 *              call to allow inserting program name into crash messages.
 * * Jul 23 13:54 1998 (edgrif): Remove naughty redeclarations of
 *      fmap functions and instead include the fmap.h public header.
 * * Mar  4 22:33 1996 (rd): built from tacemain.c
 * Created: Mon Mar  4 22:32:03 1996 (rd)
 *-------------------------------------------------------------------
 */

/******************************************************/
/********* ACEDB gif generator for WWW system  ********/
/******************************************************/

#include "acedb.h"
#include "fmap.h"		/* for fMapDrawGIF etc. */
#include "peptide.h"
#include "session.h"
#include "graph.h"
#include "client.h"
#include "sysclass.h"
#include "query.h"
#include "command.h"
#include "acedbgraph.h"		/* for acedbAppGraphInit */
#include "lex.h"
#include "banner.h"
#include "version.h"
#include "update.h"
#include "../whooks/sysclass.h"
#include "aceio.h"
#include "call.h"

extern int isInteractive ;	/* in freesubs */


#ifndef GIFACESERVER					    /* flag set in truemake to build gifservercontrol.o */

extern void (*gifEntry)(KEYSET, int,BOOL) ;		    /* entry point in command.c */

#endif   /* !GIFACESERVER */


BOOL graphTimerSet(int msecs,int timerIdentifier,TimerFunc func,void * param) { return TRUE ; }


void gifControl (KEYSET ks, int level, BOOL isI) ;  /* public, called by aceserver */

/*MENU mainMenu = 0 ;*/	     /* global so that session, wcs can change */


/* flag set in truemake to build		 gifservercontrol.o */
#ifndef GIFACESERVER

/* Defines a routine to return the compile date of this file.                */
UT_MAKE_GETCOMPILEDATEROUTINE()

/************************************************************/

int main (int argc, char **argv)
{
  int level ;

  messErrorInit(argv[0]) ;			  /* Record program name for crash messages. */

  setbuf (stdout, NULL) ;
  setbuf (stderr, NULL) ;
  
  gifEntry = gifControl ;

  bannerWrite (bannerMainStrings ("giface", FALSE, FALSE)) ;

  freeinit () ;			/* must come before graphInit */


  /* Initialise the graphics package and its interface for acedb.            */
  acedbAppGraphInit(&argc, argv) ;

  aceInit (argc>1 ? argv[1] : 0) ;

  bannerWrite (bannerDataInfoStrings ()); /* must come after aceInit() */
  
  level = freesetfile (stdin, 0) ;
#ifdef ACEMBLY
  acemblyInit () ;
  commandExecute (level, FALSE, TRUE, stdout, 0, 23, 0) ;
#else
  commandExecute (level, FALSE, TRUE, stdout, 0, 19, 0) ;
#endif
 
  /* finish the session and clean up */
  if (isWriteAccess() && 
      messQuery ("You did not save your work, should I ?"))
    aceQuit (TRUE);
  else
    aceQuit (FALSE);

  printf ("\n// A bientot \n\n") ;

  return(EXIT_SUCCESS) ;
}

#endif   /* !GIFACESERVER */

/******************************************************/
/******************************************************/
/********  specific giface callback controller ********/

#include "freeout.h"
#include "display.h"
#include "pick.h"

static FREEOPT options[] =
{ { 28, "giface" },
  { '?', "? : this list of commands.  Multiple commands ';' separated accepted" },

  { 'd', "dimensions x y : in pixels" }, 
  { 'D', "display [-D displaytype] [class key] : graphic display [for the current key]" },
  { 'B', "psdump filename : save active graph as postscript" },
  { 'A', "gifdump [-nobox] filename : save active graph as a gif" },
  { 'T', "svgdump  filename : save active graph as a svg" },
  { 'V', "swfdump [-nobox] filename : save active graph as a swf" },
  { 'W', "jsdump [-nobox] filename : save active graph as a js" },
  { 'b', "bubbledump filename : dump bubble infoof current graph" },
  { 'M', "mouseclick : simulate a left button click on the active graph" },

  { 'm', "map map [-view view] [-coords x1 x2] [-hideheader] [-whole]: does not write - use gif after" },

  { 'P', "pmap [-clone <name>] | [-contig <name> <coord>] : physical map" }, 

  { 'h', "mhmap [-glyph glyph_number] [-colour colour] [-noNT] : multiple IntMap display" }, 

  { 'S', "seqget [-class <class>] <sequence> [-coords x1 x2] [-origin x0] [-view view]: sets current sequence and  origin of coords  for future ops on same line" },
  { 'X', "seqdisplay [-visible_coords v1 v2] [-view view]: works on current sequence, use svg, gif or ps to actually dump the display" },
  { 'n', "seqdna [-file fname] <-coords x1 x2> : works on current sequence <coords deprecated - use seqget>" },
  { 'f', "seqfeatures [-file fname] [-version 1|2] [-list] [-source source(s)] [-feature feature(s)] <-coords x1 x2> : works on current sequence, source(s)/feature(s) are '|' separated lists, <coords deprecated - use seqget>" },
  { 'g', "seqactions  [-dna] [-gf_features] [-hide_header] [-rev_comp] : works on current sequence" },
  { 'c', "seqcolumns {-on name} {-off name} : works on active sequence" },
		/* NB seqcolumns will work BEFORE seqdisplay() now */
  { 'a', "seqalign  [-file fname] [-coords x1 x2] : must follow a seqdisplay <or sequence, which is deprecated>" },
  { 's', "sequence sequence [-coords x1 x2] : <deprecated - use seqget ; seqdisplay>" },
  { 'e', "EMBL filename : embl dump keyset to filename.embl" },

  {500, "pepget <protein | peptide> [-coords x1 x2] : sets current protein for future ops on same line" },  
  {501, "pepseq [-file fname] <-coords x1 x2> : exports current petide sequence" },
  {502, "pepalign [-file fname] <-coords x1 x2> : exports homols of current peptide" },

  { 'K', "makemaps [-gmap] [-pmap] [-cmap] [-alpha] [-all] [-seqmap file] [-seqclonemap file]: makes cached maps and sorted class keysets" },
  { 'U', "update [-all] : add one official update, or all" },

  { 'Q', "Quit" }
} ;

void gifControl (KEYSET ks, int level, BOOL isInteractive)  /* public, called by aceserver */
{ 
  KEY key ;
  char *word ;
  float xfac, yfac, sw, sh ;
  int pw, ph, xClick, yClick ;
  void *seqLook = 0, *pepLook = 0 ;
  KEY oldisGifDisplay = isGifDisplay, dumpBox = TRUE ; ;
  BOOL compile = FALSE ;

  graphScreenSize (&sw, &sh, 0, 0, &pw, &ph) ;
  xfac = sw/pw ; yfac = sh/ph ;

  isGifDisplay = 1 ; /* may be superseeded by a view value */
  while (TRUE)
    { if (isInteractive) freeOut ("acedb-gif> ") ;
      if (freelevelselect (level, &key, options))
	switch (key)
	  {
	  case '?':		/* help */
	    { int i ;
	      for (i = 1 ; i <= options[0].key ; i++)
		freeOutf ("// %s\n", options[i].text) ;
	    }
	    break ;

/*------------------------------------------------------------*/

	  case 'd':		/* dimensions in pixels */
	    if (freeint (&pw) && freeint (&ph))
              { pickSetDisplaySize (FMAP, 0, 0, xfac*pw, yfac*ph) ;
                pickSetDisplaySize (GMAP, 0, 0, xfac*pw, yfac*ph) ;
                pickSetDisplaySize (PMAP, 0, 0, xfac*pw, yfac*ph) ;
                pickSetDisplaySize (DtHSEQ, 0, 0, xfac*pw, yfac*ph) ;
                pickSetDisplaySize (DtTiling, 0, 0, xfac*pw, yfac*ph) ;
                pickSetDisplaySize (DtGLOC, 0, 0, xfac*pw, yfac*ph) ;
                pickSetDisplaySize (DtGLOCBIG, 0, 0, xfac*pw, yfac*ph) ;
                pickSetDisplaySize (DtGeneExp, 0, 0, xfac*pw, yfac*ph) ;
              }
	    else
	      freeOut ("// Usage: DIMENSIONS nx ny : image size in pixels\n") ;
	    break ;

	  case 'D':   /* display an object */
	    /* only displays single objects, multiple objects is an error. */
	    { KEYSET kA ;
	      KEY classKey = 0, displayKey = 0, view = 0 ;
	      char *cp ;
	      KEY _VView, mm;
	      
	      lexaddkey("View", &mm, _VMainClasses);
	      _VView = KEYKEY(mm);
	      
	      while ((cp = freeword()))
		{
		  if (*cp != '-')
		    { freeback() ; break ; }
		  if (!strcmp(cp, "-D")) /* -D displaytype */
		    { if (!(cp = freeword()))
			{ freeOut ("// Error: -D must be followed by a display type\n") ;
			  break ;
			}
		      if (!lexword2key(cp, &displayKey, _VDisplay) && !lexword2key (cp, &view, _VView))
			{ freeOutf ("// Error: bad display type %s\n", cp) ;
			  break ;
			}
		    }
		  else if (!strcmp(cp, "-view")) /* -D displaytype */
		    { if (!(cp = freeword()))
			{ freeOut ("// Error: -view must be followed by a view name\n") ;
			  break ;
			}
		      if (!lexword2key (cp, &view, _VView)  && !lexword2key(cp, &displayKey, _VDisplay) )
			{ freeOutf ("// Error: unknown view %s\n", cp) ;
			  break ;
			}
		    }
		  if (view)
		    {
		      gMapSelectView (view) ;
		      isGifDisplay = view ;
		    }
		}
	      if (displayKey == 1)
		break ;
	      
	      if (freecheck ("ww")) /* [class name] */
		{ cp = freeword() ;
		  if (!lexword2key (cp, &classKey, _VClass))
		    { freeOutf ("// Error: bad class name %s\n", cp) ;
		      break ;
		    }
		  kA = query(0, messprintf("Find %s \"%s\"", 
					   name(classKey), freeword())) ;
		}
	      else if (!ks)
		{ freeOut ("// Error: no active list\n") ;
		  break ;
		}
	      else
		kA = ks ;
		  
	      if (!keySetMax (kA))
		{ freeOut ("// Error: active list is empty\n") ;
		  break ;
		}
	      if (keySetMax(kA) == 1) 
		display (keySet(kA,0), 0, displayKey ? name(displayKey) : 0) ;
	      else
		freeOut ("// Error: can not display multiple objects\n") ;

	      /* Only destroy the key if we created a new one, DON'T destroy */
	      /* the keyset we were passed in !!                             */
	      if (kA != ks) keySetDestroy (kA) ;
	    }
	    break;

	  case 'B':   /* dump the active graph as PostScript */
	    word = freeword() ;
	    if (!word || !*word || !graphActive())
	      freeOut ("// Usage psdump filename: dump the active graph in file as PostScript\n") ;
	    else
	      {
	      BOOL rotation ;
	      float scale ;
	      int pages ;
	      char *fileName ;

	      fileName = messalloc(strlen(word) + 1) ;
	      strcpy (fileName, word);

	      /* Get defaults for PostScript printing and then print to file.*/
	      graphPSdefaults(&rotation, &scale, &pages) ;
	      graphPS(fileName, NULL, NULL, NULL, TRUE, rotation, scale, pages) ;

	      freeOutf ("// I wrote the active graph to file %s\n", fileName) ;
	      messfree (fileName);
	      }
	    break;
	    
	  case 'A':   /* dump as GIF of the active graph */
	  case 'T':   /* dump as SVG of the active graph */
	  case 'V':   /* dump as SWF of the active graph */
	  case 'W':   /* dump as JS of the active graph */
	    word = freeword() ;
	    dumpBox = TRUE ;
	    if (word && !strcmp (word,"-nobox"))
	      {
		dumpBox = FALSE ;
		word = freeword() ;
	      }
	    if (word && !strcmp (word,"-compile"))
	      {
		compile = TRUE ;
		word = freeword() ;
	      }
	    if (!word || !*word || !graphActive())
	      freeOut ("// Usage gifdump filename: dump the active graph in file\n") ;
	    else if (strcmp(word,"-") == 0)
	      { /* LS 2 Mar 98 dump GIF to standard output */
		AC_HANDLE h = handleCreate () ;
		Stack gifStack = stackHandleCreate (100000, h) ;
		ACEOUT fo ;
		fo = aceOutCreateToStack (gifStack, h) ;
		switch (key)
		  {
		    case 'A': 
		      graphGIF (graphActive(), fo, 0) ; 
		      freeOutBinary (stackText(gifStack, 0), stackMark (gifStack)) ;
		      break ;
		    case 'T': 
		      svgGraphExport (graphActive(), fo, dumpBox) ; 
		      dumpBox = FALSE ; /* done internally */
		      freeOutBinary (stackText(gifStack, 0), stackMark (gifStack)) ;
		      break ;
		    case 'V': 
		      swfGraphExport (graphActive(), fo, dumpBox) ; 
		      dumpBox = FALSE ; /* done internally */
		      if (! compile)
			freeOutBinary (stackText(gifStack, 0), stackMark (gifStack)) ;
		      else
			{
			  int nn = 0 ;
			  Stack swfStack =  stackHandleCreate (100000, h) ;
			  ACEOUT foTmp = aceOutCreateToStack (swfStack, h) ;

			  if (filName ("../bin/swfc",0, "x"))
			    {
			      if (1)
				nn = callPipe ("../bin/swfc -o - -" /* "grep flash" */, stackText(gifStack, 0), stackMark (gifStack), foTmp) ;
			      else
				nn = callPipe ("swfc -o stdout " /* "old 32 bit version grep flash" */, stackText(gifStack, 0), stackMark (gifStack), foTmp) ;
			    }
			  if (nn > 0)
			    freeOutBinary (stackText(swfStack, 0), nn) ;
			  else
			    freeOutf ("ERROR in gifacemain.c while calling the SWF compiler ../bin/swfc, sorry\n") ;
			}
		      break ;
		      /* case 'W': jsGraphExport (graphActive(), fo, 0) ; break ; */
		  }
		handleDestroy (h) ;
		if (dumpBox)
		  graphBoxInfoFile(NULL);
	      }
	    else
	      {
		ACEOUT fo = 0 ;
		BOOL error = FALSE ;
		FILE *fil = 0 ;
		char *end ;

		end = word + strlen(word) - 1 ;
		while (end > word && *end != '.') end-- ;
		switch (key)
		  {
		  case 'A':
		    if (!strcmp (end, ".gif"))
		      end = "" ;
		    else
		      end = ".gif" ;
		    if ((fo = aceOutCreateToFile (messprintf ("%s%s", word, end), "wb", 0))) 
		      {
			error = graphGIF (graphActive(), fo, 0) ;
			aceOutDestroy (fo) ;
		      }
		    break ;
		  case 'T':
		    if (!strcmp (end, ".svg"))
		      end = "" ;
		    else
		      end = ".svg" ;
		    if ((fo = aceOutCreateToFile (messprintf ("%s%s", word, end), "wb", 0))) 
		      {
			    error = svgGraphExport (graphActive(), fo, dumpBox) ;
			aceOutDestroy (fo) ;
		      }
		    break ;
		  case 'V': 
		    if (compile)
		      {
			if (!strcmp (end, ".swf"))
			  end = "" ;
			else
			  end = ".swf" ;
		      }
		    else
		      {
			if (!strcmp (end, ".sc"))
			  end = "" ;
			else
			  end = ".sc" ;
		      }
		    if (!compile)
		      {
			if ((fo = aceOutCreateToFile (messprintf ("%s%s", word,end), "w", 0))) 
			  {
			    error = swfGraphExport (graphActive(), fo, dumpBox) ;
			    aceOutDestroy (fo) ;
			  }
		      }
		    else
		      {
			AC_HANDLE h = handleCreate () ;
			Stack gifStack = stackHandleCreate (100000, h) ;
			ACEOUT fo1 ;  
			char *fNam = filName (word, *end ? end+1 : end, "w") ;
			if (1)
			  {
			    FILE *pipe ;
			    
			    if (1)
			      pipe = popen (messprintf ("swfc -o %s -", fNam), "w") ; /* grep flash */
			    else
			      pipe = popen (messprintf ("swfc -o %s ", fNam), "w") ; /*old 32 bit code  grep flash */
			    
			    fo1 = aceOutCreateToStack (gifStack, h) ;
			    swfGraphExport (graphActive(), fo1, 0) ; 
			    fwrite (stackText(gifStack, 0), 1, stackMark (gifStack), pipe) ; 
			    pclose( pipe );
			  }
		      }
		    break ;
		  }
			
		if (dumpBox &&
		    (fil = filopen ((char*) name, "boxes", "w")))
		  { 
		    graphBoxInfoFile (fil) ;
		    filclose (fil) ;
		  }
		
		if (error)
		  freeOutf ("// I wrote the active graph to file %s\n", word) ;
		else
		  freeOutf ("// Sorry, I could not open file %s\n", word) ;
	      }
	    break;

	  case 'b':   /* dump the bubble info */
	    word = freeword() ;
	    if (!word || !*word || !graphActive())
	      freeOut ("// Usage bubbledump filename: dump the active graph in - or file\n") ;
	    else if (graphDumpBubble (word))
	      freeOutf ("// I wrote the active graph to file %s\n", word) ;
	    else
	      freeOutf ("// Sorry, I could not open file %s\n", word) ;
	    break;

	  case 'M' :  /* mouseclick x y */
	    if (!freeint (&xClick) || !freeint (&yClick) || !graphActive())
	      freeOutf ("// Usage:  mouseclick x y - simulates a left button click in the active graph\n") ;
	    else
	      gifLeftDown (xClick, yClick) ;
	    break;

/*------------------------------------------------------------*/

	  case 'm':		/* map/view draw */
	    { extern void gMapDrawGIF (void) ; gMapDrawGIF () ; }
	    break ;

/*------------------------------------------------------------*/

	  case 'P':		/* pmap draw */
	    { extern void pMapDrawGIF (void) ; pMapDrawGIF () ; }
	    break ;

/*------------------------------------------------------------*/
	  case 'h':		/* pmap draw */
	    { 
	      extern void (*mhmpDrawGifp) (KEYSET ks) ;
	      if (mhmpDrawGifp) (*mhmpDrawGifp)(ks) ;
	      break ;
	    }

/*------------------------------------------------------------*/

	  case 'S':		/* get sequence */
	    seqLook = fMapGifGet (seqLook) ;
	    break ;

	  case 'X':		/* sequence display */
	    if (seqLook) fMapGifDisplay (seqLook) ;
	    else freeOut ("// gif seqdisplay without active sequence\n") ;
	    break ;

	  case 'n':		/* dna dump */
	    if (seqLook) fMapGifDNA (seqLook) ;
	    else freeOut ("// gif seqdna without active sequence\n") ;
	    break ;

	  case 'f':		/* feature dump */
	    if (seqLook) fMapGifFeatures (seqLook) ;
	    else freeOut ("// gif seqfeatures without active sequence\n") ;
	    break ;

	  case 'g':		/* actions */
	    if (seqLook) fMapGifActions (seqLook) ;
	    else freeOut ("// gif seqactions without active sequence\n") ;
	    break ;

	  case 'c':		/* column control */
	    if (seqLook) fMapGifColumns (seqLook) ;
	    else freeOut ("// gif seqcolumns without active sequence\n") ;
	    break ;

	  case 'a':		/* alignment */
	    if (seqLook) fMapGifAlign (seqLook) ;
	    else freeOut ("// gif seqalign without displayed active sequence\n") ;
	    break ;

	  case 's':		/* sequence draw - deprecated */
	    seqLook = fMapGifGet (seqLook) ;
	    if (seqLook) fMapGifDisplay (seqLook) ;
	    else freeOut ("// gif seqdisplay without active sequence\n") ;
	    break ;

/*------------------------------------------------------------*/

	  case 500:		/* get protein */
	    pepLook = pepGifGet (seqLook) ;
	    break ;

	  case 501:		/* dna dump */
	    if (pepLook) pepGifSeq (pepLook) ;
	    else freeOut ("// gif pepseq without active protein\n") ;
	    break ;

	  case 502:		/* feature dump */
	    if (pepLook) pepGifAlign (pepLook) ;
	    else freeOut ("// gif pepAlign without active protein\n") ;
	    break ;

/*------------------------------------------------------------*/

	  case 'e':		/* EMBL dump */
	    { extern void emblDumpKeySetFile (KEYSET kset, char *fname) ;
	      if ((word = freeword()))
		emblDumpKeySetFile (ks, word) ;
	      else
		freeOut ("// Usage: EMBL filename\n") ;
	    }
	    break ;

/*------------------------------------------------------------*/

	  case 'K':		/* make cached maps */
	    { extern void gMapMakeAll(void) ;
	      extern void pMapMakeAll(void) ;
	      extern void cMapMakeAll(void) ;
	      extern void sMapMake (FILE *fil, BOOL isSeq) ;
	      char *filname ;
	      FILE *fil ;
	      char *cardKeep ;

	/* must save rest of free card, else sessionWriteAccess() destroys */
	      cardKeep = strnew (freepos(), 0) ;
	      if (!sessionGainWriteAccess())
		{ freeOut ("// Error: can't get write access\n") ;
		  break ;
		}
	      freeforcecard (cardKeep) ; messfree (cardKeep) ;

	      while (freestep ('-') && (word = freeword()))
		{ 
		  cardKeep = strnew (freepos(), 0) ;

		  if (!strcmp (word, "gmap"))
		    gMapMakeAll () ;
		  else if (!strcmp (word, "pmap"))
		    pMapMakeAll () ;
		  else if (!strcmp (word, "cmap"))
		    cMapMakeAll () ;
		  else if (!strcmp (word, "alpha"))
		    lexAlphaMakeAll () ;
		  else if (!strcmp (word, "all"))
		  { 
		    pMapMakeAll () ;
		    gMapMakeAll () ;
		    cMapMakeAll () ;
		    lexAlphaMakeAll () ;
		  }
		  else if (!strcmp (word, "seqmap"))
		    {
		      if (!(word = freeword()))
			freeOut ("// Error: makemaps -seqmap "
				 "requires a file name\n") ;
		      else if (!(fil = fopen (word, "w")))
			freeOutf ("// Error: makemaps -seqmap "
				  "failed to open file %s\n", word) ;
		      else
			{ 
			  filname = strnew (word, 0) ;
			  sMapMake (fil, TRUE) ;
			  freeOutf ("// Sequence-* map file written to %s\n", 
				    filname) ;
			  messfree (filname) ;
			  fclose (fil) ;
			}
		    }
		  else if (!strcmp (word, "seqclonemap"))
		    {
		      if (!(word = freeword()))
			freeOut ("// Error: makemaps -seqclonemap "
				 "requires a file name\n") ;
		      else if (!(fil = fopen (word, "w")))
			freeOutf ("// Error: makemaps -seqclonemap "
				  "failed to open file %s\n", word) ;
		      else
			{ 
			  filname = strnew (word, 0) ;
			  sMapMake (fil, FALSE) ;
			  freeOutf ("// Sequence-* Clone map file written to %s\n", 
				    filname) ;
			  messfree (filname) ;
			  fclose (fil) ;
			}
		    }
		  else
		    freeOut ("// Usage: makemaps [-gmap] [-pmap] [-cmap] [-alpha] [-all] [-seqmap file] [-seqclonemap file]\n") ;

		  freeforcecard (cardKeep) ; messfree (cardKeep) ;
		}

	      while ((word = freeword()))
		{ 
		  freeOutf ("// unrecognized option %s\n", word) ;
		  freeOut  ("// Usage: makemaps [-gmap] [-pmap] [-cmap] [-alpha] [-all] [-seqmap file] [-seqclonemap file]\n") ;
		}
	    }
	    break ;

	  case 'U':		/* add updates */
	    if ((word = freeword()))
	      {
		if (!strcmp (word, "-all"))
		  updateDoAction (TRUE) ;
		else
		  freeOut ("// Usage: update [-all]\n") ;
	      }
	    else
	      updateDoAction (FALSE) ;
	    break ;

/*------------------------------------------------------------*/

	  case 'Q':
	  case (KEY)(-1):	/* end of file */
	    goto finish ;
	  }
    }

 finish:
  graphCleanUp () ;		/* kill all graphs except active graph */
  graphDestroy () ;		/* kill active graph */
  isGifDisplay = oldisGifDisplay ;
  if (seqLook) fMapGifDestroy (seqLook) ;
  if (pepLook) pepGifDestroy (seqLook) ;
}

/**************************************************/
/**************************************************/
 
 
