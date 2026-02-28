/*  File: fmapcontrol.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: display choice control for fmap
 * Exported functions:
 * HISTORY:
 * Last edited: May 18 23:29 1999 (rd)
 * * Oct 21 11:09 1998 (edgrif): Corrected func. proto for fMapFollowVirtualMultiTrace
 *              a whole load of fmap functions are nulled for ACEMBLY way down the file.
 * * Jul 27 10:36 1998 (edgrif): Add calls to fMapIntialise () for externally
 *      called routines.
 * * Jul 16 10:06 1998 (edgrif): Introduce private header fmap_.h
 * * Jun 24 14:17 1998 (edgrif): Fix initialisation bug to make sure isDone
        is set in fmapInitialise. Also remove methodInitialise from fmapInitialise,
	this is not needed, method does its own intialisation.
 * * Dec 21 09:00 1995 (mieg)
 * * Jun 21 00:18 1995 (rd): fused with Jean - remaining discrepancies in
 	readMethod (s) sections - probably he found a bug here with 1, 2 etc.
	and with HOMOL unshading in fMapSelect - I am quite confident here.
 * Created: Sat Jul 25 20:19:11 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: fmapcontrol.c,v 1.77 2019/09/26 20:38:10 mieg Exp $ */
#include "fmap_.h"
#include "bitset.h"
#include "chrono.h"
#include "dna.h"
#include "bindex.h"
#include "peptide.h"
#include "display.h"
#include "key.h"		/* Keyboard keys */
#include "sysclass.h"
#include "classes.h"
#include "systags.h"
#include "a.h"
#include "../wh/menu.h"
#include "query.h"
#include "freeout.h"
#include "orepeats.h"
#include "pref.h"

/********************** globals ******************************/

	/* addresses of these two are used as unique identifiers */
magic_t GRAPH2FMAPLOOK_ASSOC = "FMAP";	/* find the fMap-LOOK on the active graph */
magic_t FMAPLOOK_MAGIC = "FMAP";	/* verify a LOOK pointer to be an fMap-LOOK */

LOOK fMapGifLook = 0 ;		/* global LOOK when no graph for gif package */

/********************** locals *******************************/

static void fMapSelectBox (LOOK look, int box, double x, double y) ;
static void fMapFollow (LOOK look, double x, double y) ;
static void fMapCentre (LOOK look, KEY key, KEY from) ;
static void fMapPick (int box, double x, double y) ;
static void fMapComplementSurPlace (void) ;
static void fMapComplement (void) ;
static void fMapCompHelp (void) ;
static void fMapReverse (void) ;
static void fMapReverseComplement (void) ;
static void fMapDoRecalculate (LOOK look, int start, int stop) ;
static void fMapRecalculate (void) ;
static void fMapZoneKeySet (void) ;
static void fMapKbd (int k) ;
static BOOL findStartStop (KEY *seq, int *start, int *stop) ;

static LOOK fMapCreateLook (KEY seq, int start, int stop, int origin, BOOL doReverse, KEY view, int my_rView);
int fMap_rView ; /* set by fMapCreateLook, used by fMapRecalculate, but unly if view != 0 */

/************** default column specification ************/

static void fMapDefaultColumns (MAP map)
{
/* all column position 100, -90 ... should be different */
#define STL_STATUS_0 
  mapInsertCol (map, -100.0, TRUE, "Locator", fMapShowMiniSeq) ;
  mapInsertCol (map, -90.0, TRUE, "Clones", fMapShowCanonical) ;
#ifdef STL_STATUS  /* mike holman status column */
  mapInsertCol (map, -89.0, TRUE, "Cosmids by group", fMapShowOrigin) ;
#endif
  mapInsertCol (map, -2.1, FALSE, "Up Gene Translation", fMapShowUpGeneTranslation) ;
  mapInsertCol (map, -1.9,  TRUE, "-Confirmed introns", fMapShowSoloConfirmed) ;
  mapInsertCol (map, -0.8,  TRUE, "-Probes", fMapcDNAShowProbe) ;
  mapInsertCol (map, -0.1,  TRUE, "Restriction map", fMapShowCptSites) ;
  mapInsertCol (map, 0.0,  TRUE, "Summary bar", fMapShowSummary) ;
  mapInsertCol (map, 0.1,  TRUE, "Oligo Repeats", fMapcDNAShowOligoRepeats) ;
  mapInsertCol (map, 0.2,  TRUE, "Scale", fMapShowScale) ;
  mapInsertCol (map, 1.9,  TRUE, "Confirmed introns", fMapShowSoloConfirmed) ;
  mapInsertCol (map, 3.0,  TRUE, "EMBL features", fMapShowEmblFeatures) ;
#ifdef STL_STATUS
  mapInsertCol (map, 3.1,  TRUE, "Cosmid Status", fMapShowStatus) ;
#endif
  mapInsertCol (map, 3.23, FALSE, "CDS Lines", fMapShowCDSLines) ;
  mapInsertCol (map, 3.25, FALSE, "CDS Boxes", fMapShowCDSBoxes) ;
  mapInsertCol (map, 3.3,  TRUE, "Alleles", fMapShowAlleles) ;
  mapInsertCol (map, 3.4, FALSE, "cDNAs", fMapShowcDNAs) ;
  mapInsertCol (map, 3.7,  TRUE, "Assembly Tags", fMapShowAssemblyTags) ;
  mapInsertCol (map, 3.8, FALSE, "Oligos", fMapOspShowOligo) ;
  mapInsertCol (map, 3.82, FALSE, "Oligo_pairs", fMapOspShowOligoPairs) ;
/* isFrame starts if either of next 2 are On */
  mapInsertCol (map, 4.0, FALSE, "3 Frame Translation", fMapShow3FrameTranslation) ;
  mapInsertCol (map, 4.05, FALSE, "ORF's", fMapShowORF) ;
  mapInsertCol (map, 4.1, TRUE, "Coding Frame", fMapShowCoding) ;	/* only shows if isFrame */
  mapInsertCol (map, 4.2, FALSE, "ATG", fMapShowATG) ;
/* frame dependent stuff ends */
  /* mapInsertCol (map, 4.98, FALSE, "Gene Translation", fMapShowGeneTranslation) ; */
  mapInsertCol (map, 4.99, FALSE, "Down Gene Translation", fMapShowGeneTranslation) ;

#ifdef ACEMBLY
  mapInsertCol (map, 0.8,  TRUE, "Probes", fMapcDNAShowProbe) ;
  /*  mapInsertCol (map, 3.3,  TRUE, "Solexa", fMapcDNAShowSolexa) ; */
  mapInsertCol (map, 5.5,  TRUE, "Alignements", fMapShowAlignments) ;
  mapInsertCol (map, 5.6,  TRUE, "Previous Contigs", fMapShowPreviousContigs) ;
  mapInsertCol (map, 5.7,  TRUE, "Contigs", fMapShowContig) ;
  mapInsertCol (map, 5.8,  TRUE, "Trace Assembly", fMapShowTraceAssembly) ;
  mapInsertCol (map, 5.9,  TRUE, "Multiplets", fMapShowVirtualMultiplets) ;
#endif

  mapInsertCol (map, 6.0, FALSE, "Coords", fMapShowCoords) ;
  mapInsertCol (map, 6.1, FALSE, "DNA Sequence", fMapShowDNA) ;
#ifdef ACEMBLY
  mapInsertCol (map, 11.2, TRUE, "cDNA_Transcript_Decorate", 0) ;
  mapInsertCol (map, 6.2, FALSE, "Assembly DNA", fMapShowVirtualDna) ;
  mapInsertCol (map, 6.5, FALSE, "Brief Identifications", fMapShowBriefID) ;
  mapInsertCol (map, 6.5, TRUE, "Titles", fMapShowTitles) ;
  mapInsertCol (map, 9.9, TRUE, "GeneWalls", fMapShowGeneWalls) ;
  mapInsertCol (map, 10.0, FALSE, "Gene Names", fMapShowGeneNames) ;
#else
  mapInsertCol (map, 6.5, FALSE, "Brief Identifications", fMapShowBriefID) ;
#endif
  mapInsertCol (map, 10.2, TRUE, "Text Features", fMapShowText) ;
}

/*************** tags and their initialisation ***************/

static KEY _Feature ;
static KEY _Tm, _Temporary, _AceKogN ;
static KEY _EMBL_feature ;
       KEY _Arg1_URL_prefix ;		/* really ugly to make global */
       KEY _Arg1_URL_suffix ;
       KEY _Arg2_URL_prefix ;
       KEY _Arg2_URL_suffix ;
static KEY _Exon_support, _Spliced_cDNA, _Spliced_cDNA_Decorate,  _Gene_name, _Transcribed_gene, _From_gene ;
static KEY _Probe_hit, _Probe_exact_hit ;
static KEY _VIntron = 0, _VTranscribed_gene = 0, _VTranscript = 0 , _VmRNA = 0,  _VmProduct = 0, _VRNAi = 0, _VIST = 0, _VOST = 0 ;
static KEY _Fmap_Header
, _Fmap_Locator
, _Fmap_Summary_bar
, _Fmap_scale
, _Fmap_scale_coords
, _Fmap_DNA 
, _Fmap_3_translation
, _Fmap_Clones
, _Fmap_Gene, _Fmap_Gene_Name, _Fmap_Cloud_gene
, _Fmap_Gene_Up, _Fmap_Gene_Down
, _Fmap_Pg, _Fmap_Pg_Up, _Fmap_Pg_Down
, _Fmap_Tg, _Fmap_Tg_Up, _Fmap_Tg_Down, _Fmap_Tg_Just_Gold, _Gold
, _Fmap_mRNA, _Fmap_mRNA_Up, _Fmap_mRNA_Down
, _Fmap_mProduct 
, _Fmap_Down_Gene_Translation
, _Fmap_Spliced_cDNA, _Fmap_Spliced_cDNA_Up, _Fmap_Spliced_cDNA_Down 
, _Fmap_CDS_Tiling
, _Fmap_SOLEXA, _Fmap_SOLEXA_Up, _Fmap_SOLEXA_Down
, _Fmap_RNAi, _Fmap_RNAi_Up, _Fmap_RNAi_Down
, _Fmap_OST, _Fmap_OST_Up, _Fmap_OST_Down
, _Fmap_Probe
, _Fmap_Oligo
, _Fmap_Homol, _Fmap_AceKogN
, _Fmap_Solexa
, _Fmap_Allele
, _Fmap_Text_Features
;
/* mhmp 25.09.98 */
static int weight, seqDragBox, nbTrans ;

/* This routine is called whenever an fmap interface routine is called       */
/* to make sure that the fmap package is initialised, initialisation is      */
/* only done on the first call.                                              */

void fMapInitialise (void)
{
  static int isDone = FALSE ;
  KEY key ;

  if (isDone) return ;
  isDone = TRUE ;

  _Feature = str2tag ("Feature") ;
  _EMBL_feature = str2tag ("EMBL_feature") ;
  _AceKogN = str2tag ("AceKogN") ;
  _Arg1_URL_prefix = str2tag ("Arg1_URL_prefix");
  _Arg1_URL_suffix = str2tag ("Arg1_URL_suffix") ;
  _Arg2_URL_prefix = str2tag ("Arg2_URL_prefix") ;
  _Arg2_URL_suffix = str2tag ("Arg2_URL_suffix") ;
  _From_gene  = str2tag ("From_gene") ;
  _Exon_support = str2tag ("Exon_support") ;
  _Spliced_cDNA  = str2tag ("Spliced_cDNA") ;
  _Spliced_cDNA_Decorate  = str2tag ("Spliced_cDNA_Decorate") ;
  _Gene_name = str2tag ("Gene_name") ;
  _Transcribed_gene = str2tag ("Transcribed_gene") ;
  _Tm = str2tag ("Tm") ;
  _Temporary = str2tag ("Temporary") ;
  _Probe_hit = str2tag ("Probe_hit") ;
  _Probe_exact_hit = str2tag ("Probe_exact_hit") ;
  if (lexword2key ("Intron", &key, _VMainClasses))
    _VIntron = KEYKEY (key) ;
  if (lexword2key ("Transcribed_gene", &key, _VMainClasses))
    _VTranscribed_gene = KEYKEY (key) ;
  if (lexword2key ("Transcript", &key, _VMainClasses))
    _VTranscript = KEYKEY (key) ;
  if (lexword2key ("mRNA", &key, _VMainClasses))
    _VmRNA = KEYKEY (key) ;
  if (lexword2key ("Product", &key, _VMainClasses))
    _VmProduct = KEYKEY (key) ;
  if (lexword2key ("RNAi", &key, _VMainClasses))
    _VRNAi = KEYKEY (key) ;
  if (lexword2key ("OST", &key, _VMainClasses))
    _VOST = KEYKEY (key) ;
  if (lexword2key ("IST", &key, _VMainClasses))
    _VIST = KEYKEY (key) ;
  
  /* initialise fmapView stuff. i.e. all tag in the model View->Fmap */
   _Fmap_Header = str2tag ("Fmap_Header") ;
   _Fmap_Locator = str2tag ("Fmap_Locator") ;
   _Fmap_Summary_bar = str2tag ("Fmap_Summary_bar") ;
   _Fmap_scale = str2tag ("Fmap_scale") ;
   _Fmap_scale_coords = str2tag ("Fmap_scale_coords") ;
   _Fmap_DNA  = str2tag ("Fmap_DNA") ;
   _Fmap_3_translation = str2tag ("Fmap_3_translation") ;
   _Fmap_Clones = str2tag ("Fmap_Clones") ;
   _Fmap_Gene = str2tag ("Fmap_Gene") ;
   _Fmap_Gene_Name = str2tag ("Fmap_Gene_Name") ;
   _Fmap_Cloud_gene = str2tag ("Fmap_Cloud_gene") ;
   _Fmap_Gene_Up = str2tag ("Fmap_Gene_Up") ;
   _Fmap_Gene_Down = str2tag ("Fmap_Gene_Down") ;
   _Fmap_Pg = str2tag ("Fmap_Pg") ;
   _Fmap_Pg_Up = str2tag ("Fmap_Pg_Up") ;
   _Fmap_Pg_Down = str2tag ("Fmap_Pg_Down") ;
   _Fmap_Tg = str2tag ("_Fmap_Tg") ;
   _Fmap_Tg_Up = str2tag ("_Fmap_Tg_Up") ;
   _Fmap_Tg_Down = str2tag ("_Fmap_Tg_Down") ;
   _Fmap_Tg_Just_Gold = str2tag ("_Fmap_Tg_Just_Gold") ;
   _Gold = str2tag ("Gold") ;
   _Fmap_mRNA = str2tag ("Fmap_mRNA") ;
   _Fmap_mRNA_Up = str2tag ("Fmap_mRNA_Up") ;
   _Fmap_mRNA_Down = str2tag ("Fmap_mRNA_Down") ;
   _Fmap_mProduct = str2tag ("Fmap_mProduct") ;
   _Fmap_Down_Gene_Translation = str2tag ("Fmap_Down_Gene_Translation") ;
   _Fmap_Spliced_cDNA = str2tag ("Fmap_Spliced_cDNA") ; 
   _Fmap_Spliced_cDNA_Up = str2tag ("Fmap_Spliced_cDNA_Up") ; 
   _Fmap_Spliced_cDNA_Down = str2tag ("Fmap_Spliced_cDNA_Down") ; 
   _Fmap_CDS_Tiling = str2tag ("Fmap_CDS_Tiling") ;
   _Fmap_RNAi = str2tag ("Fmap_RNAi") ;
   _Fmap_RNAi_Up = str2tag ("Fmap_RNAi_Up") ;
   _Fmap_RNAi_Down = str2tag ("Fmap_RNAi_Down") ;
   _Fmap_SOLEXA = str2tag ("Fmap_SOLEXA") ;
   _Fmap_SOLEXA_Up = str2tag ("Fmap_SOLEXA_Up") ;
   _Fmap_SOLEXA_Down = str2tag ("Fmap_SOLEXA_Down") ;
   _Fmap_OST = str2tag ("Fmap_OST") ;
   _Fmap_OST_Up = str2tag ("Fmap_OST_Up") ;
   _Fmap_OST_Down = str2tag ("Fmap_OST_Down") ;
   _Fmap_Probe = str2tag ("Fmap_Probe") ;
   _Fmap_Oligo = str2tag ("Fmap_Oligo") ;
   _Fmap_Homol = str2tag ("Fmap_Homol") ;
   _Fmap_AceKogN = str2tag ("Fmap_AceKogN") ;
   _Fmap_Solexa = str2tag ("Fmap_Solexa") ;
   _Fmap_Allele = str2tag ("Fmap_Allele") ;
   _Fmap_Text_Features = str2tag ("Fmap_Text_Features") ;

#ifdef ACEMBLY
  topMargin = 6.5 ;
#else
  topMargin = 4.5 ;		/* fixes minor bug from Sam */
#endif
}

int fmapView (KEY view, KEY tag)
{
  if (view > 1)
    {
      if (keyFindTag (view, tag))
	return 2 ; /* specified true */
      else
	return 0 ; /* specified false */
    }

  return 1 ; /* default */
}

/***************************************************/
/******** Communication toolKit ********************/

static LOOK selectedfMap = 0 ; 

static void fMapUnselect (void)
{
  if (!selectedfMap)
    return ;

  if (selectedfMap->magic != &FMAPLOOK_MAGIC)
    messcrash ("fMapUnselect received corrupted selectfMap->magic");

  if (selectedfMap->selectBox)
    { Graph hold = graphActive () ;
      graphActivate (selectedfMap->graph) ;
      graphBoxDraw (selectedfMap->selectBox, WHITE, WHITE) ;
      graphActivate (hold) ;
    }
  selectedfMap = 0 ;
}

/***************/

static void fMapSelect (LOOK look)
{
  if (selectedfMap && selectedfMap != look)
    fMapUnselect () ;
  selectedfMap = look ;
  if (look->selectBox)
    graphBoxDraw (look->selectBox,BLACK,LIGHTRED) ;
}


/***************/
/* can be called with NULL adresses for undesired parameters */
BOOL fMapActive (Array *dnap, Array *dnaColp, KEY *seqKeyp, void** lookp)
{
  fMapInitialise () ;

  if (selectedfMap && graphExists (selectedfMap->graph))
    { if (selectedfMap->magic != &FMAPLOOK_MAGIC)
	messcrash ("fMapActive has corrupt selectedfMap->magic");

      if (lookp) *lookp = selectedfMap ;
      if (seqKeyp) *seqKeyp = selectedfMap->seqKey ;
      if (dnap) *dnap = selectedfMap->dna ;
      if (dnaColp) *dnaColp = selectedfMap->colors ;
      return TRUE ;
    }
  else
    { selectedfMap = 0 ;
      if (lookp) *lookp = 0 ; 
      if (seqKeyp) *seqKeyp = 0 ; 
      if (dnap) *dnap = 0 ; 
      if (dnaColp) *dnaColp = 0 ; 
      return FALSE ;
    }
}

Graph fMapActiveGraph (void)
{ void *v ; 
  fMapInitialise () ;
  return (fMapActive (0,0,0,&v) ? ((LOOK) v)->graph : 0) ;
}


BOOL fMapActivateGraph (void)
{ void *v ; 
  return (fMapActive (0,0,0,&v) && graphActivate (( (LOOK) v)->graph)) ;
}

/*****************************************************************/
/*************************** Header ******************************/

void fMapToggleDna (void) /* also called in fmaptrace */
{
  FMAPLOOKGET ("fMapToggleDNA") ;

  mapColSetByName ("DNA Sequence", 2) ; /* toggle */
  mapColSetByName ("Coords",  /* set to same state */
		   mapColSetByName ("DNA Sequence", -1)) ;
  fMapDraw (look,0) ; 
}

static void fMapToggleTranslation (void)
{
  int min, max ;
  Array cds, index ;
  FMAPLOOKGET ("fMapToggleTranslation") ;

  if (fMapGetCDS (look, look->translateGene, &cds, &index))
    {
      min = arr (index,0,int) - look->min ; 
      max = arr (index,arrayMax (index)-1,int) - look->min ; 
      if (min > max) 
	mapColSetByName ("Up Gene Translation", 2) ; 
      else
	mapColSetByName ("Down Gene Translation", 2) ; 
      fMapDraw (look,0) ;
    }
  else
    messout ("Select a gene to translate with the "
	     "pulldown menu on the genes") ;
}

static void zoom1 (void)
{ FMAPLOOKGET ("zoom1") ; look->map->mag = 1 ; fMapDraw (look, 0) ; }

static void zoom10 (void)
{ FMAPLOOKGET ("zoom10") ; look->map->mag = 0.1 ; fMapDraw (look, 0) ; }

static void zoom100 (void)
{ FMAPLOOKGET ("zoom100") ; look->map->mag = 0.01 ; fMapDraw (look, 0) ; }

static void zoom1000 (void)
{ FMAPLOOKGET ("zoom1000") ; look->map->mag = 0.001 ; fMapDraw (look, 0) ; }

static void fMapWhole (void)	/* simply show whole object */
{ FMAPLOOKGET ("fMapWhole") ; display (look->seqKey, 0, 0) ; }

#ifdef ACEMBLY
extern MENUOPT fMapAcemblyOpts[] ;
static void fMapTraceJoin (void)
{ FMAPLOOKGET ("fMapTraceJoin") ;
  messout ("This is a menu.  Please use the right mouse button.") ;
}
#else
static void toggleColumnButtons (void)
{ FMAPLOOKGET ("toggleColumnButtons") ;
  look->flag ^= FLAG_COLUMN_BUTTONS ;
  fMapDraw (look, 0) ;
}
#endif

static MENUOPT buttonOpts[] = {
#ifndef ACEMBLY
  { toggleColumnButtons, "Columns"},	/* NB order important */
#endif
  { mapZoomIn, "Zoom In.."},
  { mapZoomOut, "Zoom Out.."},
#ifdef ACEMBLY
  /* dna... called in fmaptrace */
  { fMapReverseComplement, "Complement"},
  { dnaAnalyse, "Analysis.."},
  { fMapTraceJoin, "Assembly.."},
#else
  { fMapClear, "Clear"},
  { fMapReverseComplement, "Rev-Comp.."},
  { fMapToggleDna, "DNA.."},
  { dnaAnalyse, "Analysis.."},
#endif
  { fMapAddGfSegs, "GeneFind.."},
  {0, 0}} ;

static MENUOPT fMapZoomOpts[] = {
  { zoom1, "1 bp/line"},
  { zoom10, "10 bp/line"},
  { zoom100, "100 bp/line"},
  { zoom1000, "1000 bp/line"},
  { fMapWhole, "Whole"},
  { 0, 0 }
} ;

static void fMapChromosomeOrigin (void)
{ 
  KEY source, seq ;
  int i ;
  FMAPLOOKGET ("fMapChromosomeOrigin") ;
  
  seq = look->seqKey ;
  while ((source = keyGetKey (seq, _Source)))
    seq = source ; 
  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
    if (arrp (look->segs,i,SEG)->key == seq)
      { 
	look->origin = arrp (look->segs,i,SEG)->x1 ;
	fMapDraw (look,0) ; /* since i want to rewrite the coordinate system */
	return ;
      }
  

  messout ("Sorry, object %s is not being displayed in this graph",
	   name (seq)) ;
  return ;
}

static MENUOPT fMapOriginOpts[] = {
  { fMapChromosomeOrigin, "Chromosome"},
  { 0, 0 }
} ;

#ifndef ACEMBLY
static MENUOPT fMapDnaOpts[] = {
  { fMapToggleDna, "DNA On-Off" }, 
  { fMapToggleTranslation, "Translate Gene On-Off" },
  { fMapSetDnaWidth, "Set DNA Width" }, 
  { fMapToggleAutoWidth, "Toggle Auto-DNA-width" },
  { 0, 0 }
} ;
#endif

extern void gelDisplay (void) ;
extern void abiFixFinish (void) ;
extern void dnaCptFindRepeats (void) ;
extern void fMapGelDisplay (void) ;
extern void fMapOspInit (void) ;
extern void fMapBlastInit (void) ;

static MENUOPT fMapAnalysisOpts[] = {
  { menuSpacer, "Built in tools"},
  { dnaAnalyse, "   Restriction Analysis, codon usage etc." }, 
  { fMapGelDisplay, "   Artificial gels" },
  { fMapAddGfSegs, "   Genefinder (Green)" }, 
  { menuSpacer, "External tools"},
  { fMapBlastInit, "   BLAST: search external database (Altschul & al)" }, 
  { fMapOspInit, "   OSP: Oligos Selection Program (Hillier)" }, 
#ifdef ACEMBLY
  { dnaCptFindRepeats, "   Repeats: Printout repeat structure (Parsons)" },
  { abiFixFinish, "   Finish: select new reads, StLouis specific (Marth)" },
#endif
  { 0, 0}
} ;

#ifndef ACEMBLY
static MENUOPT fMapComplementOpts[] = {
  { fMapReverseComplement, "Reverse-Complement" },
  { fMapReverse, "Reverse" },
  { fMapComplement, "Complement" },
  { fMapComplementSurPlace, "Complement in place" },
  { fMapCompHelp, "Help-Explanation" },
  { 0, 0}
} ;
#endif

/****************************/

static void fMapShiftOrigin (void)
{
  KEY key ;
  int i,x0 ;
  FMAPLOOKGET ("fMapShiftOrigin") ;

  fMapSelect (look) ;

  if (lexword2key (look->originBuf, &key, _VSequence))
    { for (i = 1 ; i < arrayMax (look->segs) ; ++i)
	if (arrp (look->segs,i,SEG)->key == key)
	  { look->origin = arrp (look->segs,i,SEG)->x1 ;
	    goto done ;
	  }
      messout ("Sorry, object %s is not being displayed in this graph",
	       name (key)) ;
      return ;
    }
  else if (sscanf (look->originBuf, "%d", &x0))
    { look->origin +=  x0 - 1 ;  /* -1 Plato dixit  */
      strcpy (look->originBuf, "1") ;
    }
  else
    { messout ("Give either a sequence name or a position in the "
	       "current units.") ;
      return ;
    }

 done:
  fMapDraw (look,0) ; /* since i want to rewrite the coordinate system */
}

static BOOL fMapSetOriginFromKey (LOOK look, KEY key)
{
  int i ;

  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
    if (arrp (look->segs,i,SEG)->key == key)
      { if (look->flag & FLAG_REVERSE) /* top of object is base 1 */
	  look->origin = look->length - arrp (look->segs,i,SEG)->x2 - 1 ;
        else
	  look->origin = arrp (look->segs,i,SEG)->x1 ;
	return TRUE ;
      }
  return FALSE ;
}

static void fMapDoPickOrigin (KEY key)
{
  FMAPLOOKGET ("fMapDoPickOrigin") ;

  fMapSelect (look) ;

  if (fMapSetOriginFromKey (look, key))
    fMapDraw (look,0) ; /* since i want to rewrite the coordinate system */
  else
    messout ("Sorry, %s is not being displayed in this window",
	     name (key)) ;
}

static void fMapOriginButton (void)
{
  graphRegister (MESSAGE_DESTROY, displayUnBlock) ; /*mhmp 15.05.98*/
  displayBlock (fMapDoPickOrigin, 
		"Pick an object on this window "
		"to set the first base of its sequence as base 1.\n"
		"Remove this message to cancel.") ;
}

/**********************************/

static BOOL setZone (LOOK look, int newMin, int newMax)
{ 
  BOOL res = FALSE ;

  look->zoneMin = newMin ;
  if (look->zoneMin < 0)
    { look->zoneMin = 0 ; res = TRUE ; }
  look->zoneMax = newMax ;
  if (look->zoneMax > look->length)
    { look->zoneMax = look->length ; res = TRUE ; }
  return res ;
}

static void setZoneUserCoords (LOOK look, int x0, int x1)
{ 
  int tmp ;

  if (look->flag & FLAG_REVERSE)
    { tmp = x0 ;
      x0 = look->length - look->origin - x1 ;
      x1 = look->length - look->origin - tmp + 1 ;
    }
  else
    { x0 = x0 + look->origin - 1 ;
      x1 = x1 + look->origin ;
    }

  if (x0 > x1)
    { tmp = x0 ; x0 = x1 ; x1 = tmp ;}

  if (setZone (look, x0, x1))
    messout ("The requested zone extends beyond the sequence "
	     "and is being clipped") ;

  if (look->flag & FLAG_REVERSE)
    strncpy (look->zoneBuf,
	     messprintf ("%d %d",
			 COORD (look, look->zoneMax)+1,
			 COORD (look, look->zoneMin)),
	     15) ;
  else
    strncpy (look->zoneBuf, 
	     messprintf ("%d %d",
			 COORD (look, look->zoneMin),
			 COORD (look, look->zoneMax)-1),
	     15) ;
  if (look->zoneBox) graphBoxDraw (look->zoneBox, -1, -1) ;

  fMapReDrawDNA (look) ;
}

void fMapSetZone (void)
{
  int x0, x1 ;
  FMAPLOOKGET ("fMapSetZone") ;

  fMapSelect (look) ;
  if (sscanf (look->zoneBuf,"%d %d", &x0, &x1) != 2)
    { messout ("The zone must be a pair of integers") ;
      return ;
    }

  setZoneUserCoords (look, x0, x1) ;

  /* mhmp 02.06.98  & 08.09.98 */
  graphEntryDisable () ;
  graphRegister (KEYBOARD, fMapKbd) ;
}

static void fMapDoSetZone (KEY key)
{
  int i ;
  FMAPLOOKGET ("fMapDoSetZone") ;

  fMapSelect (look) ;

  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
    if (arrp (look->segs,i,SEG)->key == key)
      { setZone (look, arrp (look->segs,i,SEG)->x1, 
		       arrp (look->segs,i,SEG)->x2 + 1) ;
	fMapDraw (look,0) ;
	return ;
      }
  messout ("Sorry, sequence %s is not being displayed in this window",
	   name (key)) ;
}

static void fMapSetZoneButton (void)
{
  graphRegister (MESSAGE_DESTROY, displayUnBlock) ; /*mhmp 15.05.98*/
  displayBlock (fMapDoSetZone, 
		"Pick an object in this window "
		"to set the limits of the active zone.\n"
		"Remove this message to cancel.") ;
}

/**********************************/

BOOL fMapFindZone (void *vv, int *min, int *max, int *origin)
{ LOOK look = (LOOK) vv ;

  fMapInitialise () ;

  if (look->magic != &FMAPLOOK_MAGIC)
    return FALSE ;
  *min = look->zoneMin ;
  *max = look->zoneMax ;
  *origin = look->origin ;
  return TRUE ;
}

/**********************************/

/** RD assumes all forward order - should use fMapFindSpan **/

BOOL fMapFindZoneFather (void *vv, int min, int max, KEY *fp, int *originp)
{ LOOK look = (LOOK) vv ;
  KEY parent = look->seqKey ;
  int i, a1 = min ;
  SEG *seg ;

  if (look->magic != &FMAPLOOK_MAGIC)
    return FALSE ;
  /* backwards, to hit the biggest parent last */
  for (i = arrayMax (look->boxIndex) - 1 ; i >= 0 ; i--) 
    { seg = BOXSEG (i) ;
      if (seg->type != SEQUENCE ||
	  seg->x1 > min || seg->x2 < max)
	continue ;
      a1 =  seg->x1 ; parent = seg->key ; break ;
    }

  *originp = min - a1 ;
  *fp = parent ;
  return TRUE ;
}

/***********************************/

static void drawHeader (LOOK look)
{
  int box ;

  look->selectBox = graphBoxStart () ; 
    graphText ("Selected DNA",1.0,.5) ; 
  graphBoxEnd () ;
  if (look == selectedfMap)
    graphBoxDraw (look->selectBox,BLACK,LIGHTRED) ;
  else
    graphBoxDraw (look->selectBox,WHITE,WHITE) ;

  box = graphButton ("Origin..:", fMapOriginButton, 14.5, .4) ;
  graphBoxMenu (box, fMapOriginOpts) ;
  if (!strlen (look->originBuf))
    strcpy (look->originBuf,"1") ;
  look->originBox =
    graphTextEntry (look->originBuf, 16, 25, .5, 
		    TEXT_ENTRY_FUNC_CAST fMapShiftOrigin) ;

  graphButton ("Active zone:", fMapSetZoneButton, 43, .4) ;
  if (look->flag & FLAG_REVERSE)
    strncpy (look->zoneBuf,
	     messprintf ("%d %d",
			 COORD (look, look->zoneMax)+1,
			 COORD (look, look->zoneMin)),
	     15) ;
  else
    strncpy (look->zoneBuf, 
	     messprintf ("%d %d",
			 COORD (look, look->zoneMin),
			 COORD (look, look->zoneMax)-1),
	     15) ;
  look->zoneBox =
    graphTextScrollEntry (look->zoneBuf, 23, 14, 56.5, .5, 
			  TEXT_ENTRY_FUNC_CAST fMapSetZone) ;

  graphText ("Centre:", 72., .5) ;
  sprintf (look->centreBuf, "%d", (int)COORD(look, look->map->centre)) ;
  look->centreBox = graphBoxStart () ;
  graphTextPtr (look->centreBuf, 80., .5, 9) ;
  graphBoxEnd () ;
  /* mhmp 02.06.98 + 08.09.98 */
  graphEntryDisable () ;
  graphRegister (KEYBOARD, fMapKbd) ;			
  if (look->flag & FLAG_COMPLEMENT)
    graphText ("COMPLEMENT", 90, 0) ;
  if (look->flag & FLAG_COMPLEMENT_SUR_PLACE)
    {
      graphText ("COMP_DNA", 101, 0) ;
      graphText ("IN_PLACE", 101, 1) ;
    }
  if (( (look->flag & FLAG_COMPLEMENT) && 
       ! (look->flag & FLAG_REVERSE)) || 
      (! (look->flag & FLAG_COMPLEMENT) &&
       (look->flag & FLAG_REVERSE)))
    graphText ("REVERSED", 90, 1) ;

  look->segBox = graphBoxStart () ;
  graphTextPtr (look->segTextBuf, 2, 2, 126) ;
  graphBoxEnd () ;
  graphBoxDraw (look->segBox, BLACK, PALEBLUE) ;

  box = graphButtons (buttonOpts, 1, 3.2, mapGraphWidth) ;
  mapColMenu (box) ;
#ifdef ACEMBLY
  box-- ; /* because i suppressed whole */
  graphBoxMenu (box + 1, fMapZoomOpts) ;
  graphBoxMenu (box + 2, fMapZoomOpts) ;
  graphBoxMenu (box + 4, fMapAnalysisOpts) ;
  graphBoxMenu (box + 5, fMapAcemblyOpts) ;
  if (look->isGeneFinder)
    graphBoxMenu (box + 6, fMapGeneOpts) ;
  else
    graphBoxClear (box + 6) ;
  topMargin = 6.5 ;
  box++ ; 
#else
  graphBoxMenu (box + 1, fMapZoomOpts) ;
  graphBoxMenu (box + 2, fMapZoomOpts) ;
  graphBoxMenu (box + 4, fMapComplementOpts) ;
  graphBoxMenu (box + 5, fMapDnaOpts) ;
  graphBoxMenu (box + 6, fMapAnalysisOpts) ;
  graphBoxMenu (box + 7, fMapGeneOpts) ;
  topMargin = 4.5 ;
#endif

  look->minLiveBox = 1 + graphBoxStart () ;
  graphBoxEnd () ;
}

/**********************************************************************/
/******************** main menu and functions for it ******************/

static void dumpSegs (void) ;
static void geneStats (void) ;
extern void emblDump (void) ;
static void hideHeaderToggle (void);
static void dumpSeq (void);
static void readFastaFile (void);

static MENUOPT fMapMenu[] = {
  { graphDestroy, "Quit" },
  { help, "Help" },
  { graphPrint, "Print Screen" },
  { mapPrint, "Print Whole Map" },
  { displayPreserve, "Preserve" },
  { fMapRecalculate, "Recalculate" },
  { hideHeaderToggle, "Hide Header" },
  { menuSpacer, "" },
  { dumpSeq, "Export Sequence" },
  { dumpSegs, "Export Features" },
  { fMapExportTranslations, "Export translations" },
  { fMapZoneKeySet, "Make a key set from the active zone" }, 
  { emblDump, "EMBL dump" },
  { readFastaFile, "Import Sequence" },
  { menuSpacer, "" },
#if defined (ACEMBLY)
  { mapColControl, "Columns"} , 
#else
  { geneStats, "Statistics" },
  { fMapAddSegments, "Read Segments" },
  { fMapClearSegments, "Clear Segments" },
#endif
  { 0, 0 }
  } ;

static void hideHeaderToggle (void)
{ FMAPLOOKGET ("hideHeaderToggle") ;
  look->flag ^= FLAG_HIDE_HEADER ;
  {
    MENUOPT *m = fMapMenu;
    
    /* toggle menu text - pretty ugly this way, but easiest way
     * using the MENUOPT system (new menu system not final yet) */
    while (m->text && m->f)
      {
	if (m->f == hideHeaderToggle)
	  {
	    if (look->flag & FLAG_HIDE_HEADER)
	      m->text = "Show Header";
	    else
	      m->text = "Hide Header";
	  }
	++m;
      }
  }
  fMapDraw (look, 0) ;
}

static void dumpSeq (void)
{ extern void dnacptFastaDump (void);
  extern void dnacptDontUseKeySet (void) ;
  dnaAnalyse () ;
  dnacptDontUseKeySet () ;
  dnacptFastaDump () ;
}

static void readFastaFile (void)
{ 
  KEY key ;
  FILE *fil ;
  char c ;
  int level ;

  if (! (fil = filqueryopen (0, 0, "", "r", 
			    "DNA sequence file in plain or FASTA format")))
    return ;

  c = getc (fil) ; ungetc (c,fil) ; /* peek at first character */
  level = freesetfile (fil, 0) ;
  if (c == '>')
    freecard (level) ;		/* throw out title line */

  lexaddkey ("temp_seq", &key, _VDNA) ;
  if (!dnaParse (level, key))
    { freeclose (level) ;	/* close file */
      return ;
    }
  freeclose (level) ;

  lexaddkey ("temp_seq", &key, _VSequence) ; /* now display it */
  display (key, 0, FMAP) ;
}

/**********************************************************************/

void fMapPleaseRecompute (void *vp)  /* set look->pleaseRecompute */
{ LOOK look = (LOOK) vp ;

  if (!look) return ;

  fMapInitialise () ;

  if (look->magic == &FMAPLOOK_MAGIC &&
      graphExists (look->graph))
    look->pleaseRecompute = TRUE ;
    
  /* trace.c may called using a stale look, 
     i do not want to crash on that, the really correct solution
     would be inside fmapdestroy to deregister this look from
     the trace,c code, 
     i do not think it is worth the effort

    messcrash ("fMapPleaseRecompute received corrupt (LOOK)vp->magic");
     */
}

/**********************************************************************/
static int fromTrace = 0, fromTracePos = 0 ;
  /* specialised entry point */

void fMapDrawFromTrace (LOOK look, KEY from, int xx, int type)
{
  fMapInitialise () ;

  fromTrace = type ;
  fromTracePos = xx ;
  fMapDraw (look, from) ;
}

void fMapDraw (LOOK look, KEY from)
{
  SEG *seg = 0 ;
  char *cp ;
  int i, fromIndex = 0 ;
  float absMag ;

  fMapInitialise () ;

  if (look && look->map && graphActivate (look->graph))
    graphPop () ;
  else
    return ;
  
  absMag = (look->map->mag > 0) ? look->map->mag : -look->map->mag ;

  chrono ("fMapDraw") ;

 restart:
  if (!class (from) && look->activeBox && 
      look->activeBox < arrayMax (look->boxIndex) &&
      arr (look->boxIndex, look->activeBox, int) < arrayMax (look->segs) &&
      BOXSEG (look->activeBox)->type != DNA_SEQ)
    { fromIndex = arr (look->boxIndex, look->activeBox, int) ;
      from = BOXSEG (look->activeBox)->parent ;
    }

  /* clear all DNA highlighting */
  if (arrayExists (look->colors))
    { cp = arrp (look->colors, 0, char);
      for ( i=0; i<arrayMax (look->colors); i++)
	*cp++ &= ~ (TINT_HIGHLIGHT1 |  TINT_HIGHLIGHT2);
    }

  graphClear () ;
  look->summaryBox = 0 ;
  look->mrnaOffset = 0 ;
  graphFitBounds (&mapGraphWidth, &mapGraphHeight) ;
  graphMenu (fMapMenu) ;
  if (!look->isOspSelecting)
    graphRegister (PICK, (GraphFunc) fMapPick) ;

  if (look->view) /* impose a minimal width to the graph */
    graphText (".", 1, 60) ; 

  arrayMax (look->segs) = look->lastTrueSeg ;
  look->boxIndex = arrayReCreate (look->boxIndex,200,int) ;
  look->visibleIndices = arrayReCreate (look->visibleIndices,64,int) ;

  if (!fmapView (look->view, _Fmap_Header) || (look->flag & FLAG_HIDE_HEADER))
    { look->selectBox = look->originBox = 
	look->zoneBox = look->segBox = 0 ;
      look->minLiveBox = 1 ;
      topMargin = 1 ;
    }
  else
    drawHeader (look) ;

  if (fmapView (look->view, _Fmap_Down_Gene_Translation) == 2)
    mapColSetByName ("Down Gene Translation" , 1 ) ;
  halfGraphHeight = 0.5* (mapGraphHeight - topMargin) ;

  if (class (from) == _VCalcul)
    { i = KEYKEY (from) ;
      look->map->centre = i ;
      from = 0 ;
    }    
  else if (fromTrace && arrayExists (look->segs))
    { for (i = 0, seg = arrp (look->segs,0,SEG) ; i < arrayMax (look->segs) ; ++i, ++seg)
        if (from == seg->key)
	  { look->map->centre = seg->x1 + fromTracePos ;
	    break ;
	  }
      if (i == arrayMax (look->segs))
	  look->pleaseRecompute = 1 ;
      else if (fromTrace == 1)
	{ from = 0 ; fromIndex = 0 ; }     /* no highlighting */
    }
  fromTrace = 0 ;

  /* mhmp 05.02.98 */
  if (look->map->centre < look->map->min) 
    look->map->centre = look->map->min ;
  if (look->map->centre > look->map->max) 
    look->map->centre = look->map->max ;

  look->min = look->map->centre - (halfGraphHeight-.3)/absMag ;
  look->max = look->map->centre + (halfGraphHeight-.3)/absMag ;
  if (look->pleaseRecompute != -1)
    { if (look->flag & FLAG_COMPLEMENT)
	{ if (look->pleaseRecompute ||
	      (look->min < 0 && look->stop < look->fullLength-1) ||
	      (look->max > look->length && look->start))
	    { look->pleaseRecompute = -1 ;
	      fMapDoRecalculate (look, look->min, look->max) ;
	      goto restart ;
	    }
	}
    else
      { if (look->pleaseRecompute ||
	    (look->min < 0 && look->start) || 
	    (look->max > look->length && look->stop < look->fullLength-1))
	  { look->pleaseRecompute = -1 ;
	    fMapDoRecalculate (look, look->min, look->max) ;
	    goto restart ;
	  }
      }
    }
  look->pleaseRecompute = FALSE ;

  if (look->min < 0)
    look->min = 0 ;
  while (look->min % 3)		/* establish standard frame */
    ++look->min ;
  if (look->max >= look->length)
    look->max = look->length ;

  look->map->fixFrame = &look->min ;	/* fiddle for column drawing */
  seqDragBox = 0 ;/* mhmp 07.10.98 */
  mapDrawColumns (look->map) ;

  /* Display column control buttons...this will be in a separate window.     */
  if (look->flag & FLAG_COLUMN_BUTTONS)  mapDrawColumnButtons (look->map) ;

  look->activeBox = 0 ;
  if (fromIndex)
    for (i = look->minLiveBox ; i < arrayMax (look->boxIndex) ; ++i)
      { if (arr (look->boxIndex, i, int) == fromIndex)
	  { fMapSelectBox (look, i, 0, 0) ;
	    break ;
	  }
      }
  else if (class (from))
    for (i = look->minLiveBox ; i < arrayMax (look->boxIndex) ; ++i)
      { if (BOXSEG (i)->key == from)
	  { fMapSelectBox (look, i, 0, 0) ;
	    break ;
	  }
      }

  chrono ("f-graphRedraw") ;
  if (1 || !mapColSetByName ("Contigs", -1))
    graphRetitle (name (from ? from: look->seqKey)) ;
  graphRedraw () ;

  fMapSelect (look) ;
  chronoReturn () ;
} /* fmapDraw */

/************************************************************/
/****************** Registered routines *********************/

static void fMapDestroy (void)	/* fMap window destructor function */
{
  FMAPLOOKGET ("fMapDestroy") ;

/* RD 981216 Why?  This is wrong.  We may be destroying an fmap 
   independent of a display block somewhere else on the screen.
   kset and tree only displayUnBlock () when they know they are blocking.
   Anyway, the message will be destroyed when the window is by the graph
   package, and you registered that to unblock already.
  displayUnBlock () ;
*/

  handleDestroy (look->handle) ;
  arrayDestroy (look->dna) ;
  arrayDestroy (look->dnaR) ;

  arrayDestroy (look->sites) ;	/* sites come from dnacpt.c */
  stackDestroy (look->siteNames) ;

  arrayDestroy (look->multiplets) ;
  arrayDestroy (look->oligopairs) ;
  stackDestroy (look->geneNamesStack) ;

  mapDestroy (look->map) ;

  if (graphActivate (look->graph))
    { graphAssRemove (&GRAPH2FMAPLOOK_ASSOC) ;
      graphAssRemove (&MAP2LOOK_ASSOC) ;
    }
   
  fMapTraceDestroy (look) ;
  fMapOspDestroy (look) ;

  if (graphActivate (look->selectGraph))   /* from fmapgene.c */
    graphDestroy () ;

  if (look == selectedfMap)
    selectedfMap = 0 ;
  look->magic = 0 ;

  messfree (look) ;
} /* fMapDestroy */

/*************************************************************/

static void fMapDoRecalculate (LOOK look, int start, int stop)
{
  int diff ;
  BOOL isComplement = (look->flag & FLAG_COMPLEMENT) != 0 ;

  if (isComplement)
    { start = look->stop + 2 - start ;
      stop = look->stop + 2 - stop ;
    }
  else
    { start += look->start ;
      stop += look->start ;
    }

  arrayDestroy (look->dna) ;
  arrayDestroy (look->dnaR) ;
  fMapTraceDestroy (look) ;  
  fMapOspDestroy (look) ;
  if (! (look->dna = fMapFindDNA (look->seqKey, &start, &stop)))
    { messerror ("Can't find DNA in fMapDoRecalculate") ; /* mhmp */
      return ;
    }
  look->colors = arrayReCreate (look->colors, arrayMax (look->dna), char) ;
  arrayMax (look->colors) = arrayMax (look->dna) ;
  look->map->max = arrayMax (look->dna) ;

  if (isComplement)
    diff = (start-1) - look->stop ;
  else
    diff = look->start - (start-1) ;
  look->map->centre += diff ;
  look->zoneMin +=diff ;
  look->zoneMax +=diff ;

  if (look->flag & FLAG_REVERSE)
    look->origin = look->length + 1 - look->origin ;
  look->length = arrayMax (look->dna) ;
  look->origin += diff ;
  if (look->flag & FLAG_REVERSE)
    look->origin = look->length + 1 - look->origin ;

  if (isComplement)
    { look->start = stop-1 ;
      look->stop = start-1 ;
    }
  else
    { look->start = start-1 ;
      look->stop = stop-1 ;
    }
  if (!fMapConvert (look, TRUE))
    messcrash ("Can't convert in fMapDoRecalculate") ;

  look->activeBox = 0 ;	seqDragBox = 0 ;
	/* since segs array reorganised */
} /* fMapDoRecalculate  */

static void fMapRecalculate (void) /* menufunc for fMapMenu */
{
  int start, stop ;
  float absMag ;
  KEY curr ;
  FMAPLOOKGET ("fMapRecalculate") ;
  
  if (look->activeBox)
    { curr = BOXSEG (look->activeBox)->parent ;
      if (!class (curr))
	curr = 0 ;	/* avoids key being a tag */
    }
  else
    curr = 0 ;

  absMag = (look->map->mag > 0) ? look->map->mag : -look->map->mag ;
  start = look->map->centre - (halfGraphHeight-2)/absMag ;
  stop = look->map->centre + (halfGraphHeight-2)/absMag ;

  fMapDoRecalculate (look, start, stop) ;
  fMapDraw (look, curr) ;
}

/***********************************/

static double oldx, oldy, oldx2, oldy2, oldx3, oldx4 ;
static float newx2, newy2 , xdeb, xfin, xlimit ;
static float newx, newy ; /*mhmp 22.04.98 */
static int oldbox ;
static LOOK oldlook ;

static void fMapLeftDrag (double x, double y)
{
  graphXorLine (0, oldy, mapGraphWidth, oldy) ;
  oldy = y ;
  int centre = COORD (oldlook, y) ;
  sprintf (oldlook->centreBuf, "%d", centre) ;
  graphBoxDraw (oldlook->centreBox, BLACK, WHITE) ;
  graphXorLine (0, y, mapGraphWidth, y) ;
}

static void frameDNA (void)
{
  FMAPLOOKGET ("frameDNA") ;

  if (look->flag & FLAG_REVERSE)
    {
      graphXorLine (newx2, newy2, xfin, newy2) ;
      graphXorLine (xfin, newy2, xfin, oldy2+1) ;
      graphXorLine (xfin, oldy2+1, oldx2, oldy2+1) ;
      graphXorLine (oldx2, oldy2+1, oldx2, oldy2) ;
      graphXorLine (oldx2, oldy2, xdeb, oldy2) ;
      graphXorLine (xdeb, oldy2, xdeb, newy2-1) ;  
      graphXorLine (xdeb, newy2-1, newx2, newy2-1) ;
      graphXorLine (newx2, newy2-1, newx2, newy2) ;
    }
  else
    {
      graphXorLine (oldx2, oldy2, xfin, oldy2) ;
      graphXorLine (xfin, oldy2, xfin, newy2-1) ;
      graphXorLine (xfin, newy2-1, newx2, newy2-1) ;
      graphXorLine (newx2, newy2-1, newx2, newy2) ;
      graphXorLine (newx2, newy2, xdeb, newy2) ;
      graphXorLine (xdeb, newy2, xdeb, oldy2+1) ;  
      graphXorLine (xdeb, oldy2+1, oldx2, oldy2+1) ;
      graphXorLine (oldx2, oldy2+1, oldx2, oldy2) ;
    }
}

static void fMapLeftDrag2 (double x, double y)
{
  float newx1, newy1, xx1, xx2, yy1, yy2 ;
  int box ;
  FMAPLOOKGET ("fMapLeftDrag2") ;

  frameDNA () ;
  oldx3 = oldx3 + x - oldx4 ;
  box =graphBoxAt (x, y, &newx1, &newy1) ;

  if  (x < xlimit &&
       (box > oldbox || (box == oldbox && newx1 >= (int)oldx)))
    {
      seqDragBox = box ;
      newx = newx1 ;
      newy = newy1 ;
      graphBoxDim (seqDragBox, &xx1, &yy1, &xx2, &yy2) ;
      newy2 = (int)newy + yy1 + 1;
      newx2 = (int)newx + xx1 + 1;
    }
  oldx4 = x ;
  frameDNA () ;
}
static void fMapLeftDrag3 (double x, double y)
{
  float newx1, newy1, xx1, xx2, yy1, yy2 ;
  int box ;
  FMAPLOOKGET ("fMapLeftDrag3") ;

  frameDNA () ;
  box =graphBoxAt (x, y, &newx1, &newy1) ;
  if  (x < xlimit &&
      (box > oldbox || (box == oldbox && newx1 >= (int)oldx)))
    {
      seqDragBox = box ;
      newx = newx1 ;
      newy = newy1 ;
      graphBoxDim (seqDragBox,&xx1, &yy1, &xx2, &yy2) ;
      /*   xdeb = xx1 ;
      xfin = xx2 ;*/
      newy2 = (int)newy + yy1 + 1;
      newx2 = (int)newx + xx1 + 1;
    }
  frameDNA () ;
}

static void fMapLeftUp (double x, double y)
{
  float pos, pos1, min = 1000000 ;/* mhmp 30.04*/
  int i ;
  SEG *seg ;
  KEY keymin = 0 ;
  FMAPLOOKGET ("fMapLeftUp") ;

  graphXorLine (0, oldy, mapGraphWidth, oldy) ;
  graphRegister (LEFT_DRAG, 0) ;
  graphRegister (LEFT_UP, 0) ;

  if (x < look->map->thumb.x)	/* show closest clone or locus! */
    { pos1 =  WHOLE2MAP (look->map, y) ;
      for (i = 1 ; i < arrayMax (look->segs) ; ++i)
	{ seg = arrp (look->segs,i,SEG) ;
	  if (class (seg->key) == _VClone || class (seg->key) == _VLocus)
	    {
	      pos = (seg->x1 + seg->x2)*0.5 ;
	      if (fabs (pos - pos1) < min)
		{ min = fabs (pos - pos1) ; keymin = seg->key ;}
	    }
	}
      if (keymin)
	{
	  if (class (keymin) == _VClone)
	    display (keymin, 0, PMAP) ;
	  else
	    display (keymin, 0, GMAP) ; 
	}
    }
}

static void fMapLeftUp2  (double x, double y)/*mhmp 22.04.98 */
{ 
  frameDNA () ;
  graphRegister (LEFT_UP, 0) ;
  graphRegister (LEFT_DRAG, 0) ;
  fMapSelectBox (oldlook, oldbox, oldx, oldy) ;
}

/***************************************************************/

static void fMapPick (int box, double x, double y)
{
  float x1,x2,y1,y2 ;
  int xxx, ii ;
  FMAPLOOKGET ("fMapPick") ;
  SEG *seg = BOXSEG (box) ;
  
  if (look != selectedfMap)
    fMapSelect (look) ; 

  if (!box)			/* cursor line */
    { graphBoxDim (box, &x1, &y1, &x2, &y2) ;
      y += y1 ;
      oldy = y ;
      oldlook = look ;
      graphColor (BLACK) ;
      graphXorLine (0, y, mapGraphWidth, y) ;
      graphRegister (LEFT_DRAG, (GraphFunc)fMapLeftDrag) ;
      graphRegister (LEFT_UP, (GraphFunc)fMapLeftUp) ;
    }
  else if (box == look->map->thumb.box)
    {  /* coordinates are relative to thumb.box */
      if (!isGifDisplay)
	graphBoxDrag (look->map->thumb.box, mapThumbDrag) ;
      else
	{
	  look->map->centre = 
	    WHOLE2MAP (look->map, y + MAP2WHOLE (look->map,topMargin + 1)) ;
	  fMapDraw (look,0) ;
	}
    }
  else if (box == look->originBox)
    graphTextEntry (look->originBuf,0,0,0,0) ;
  else if (box == look->zoneBox)
    graphTextEntry (look->zoneBuf,0,0,0,0) ;
  else if (box >= look->minLiveBox &&
	   box != look->summaryBox && 
	   seg)
    {
      if (box == look->activeBox &&
	  BOXSEG (box)->type != DNA_SEQ &&
	  BOXSEG (box)->type != TRANS_SEQ &&
	  BOXSEG (box)->type != TRANS_SEQ_UP &&
	  BOXSEG (box)->type != PEP_SEQ)
	fMapFollow (look, x, y) ;
      else if (BOXSEG (box)->type == DNA_SEQ ||
	       BOXSEG (box)->type == PEP_SEQ)
	{
	  graphBoxDim (box,&x1, &y1, &x2, &y2) ; 
	  graphColor (BLACK) ;
	  oldx = x ;
	  oldy = y;
	  oldlook = look ;
	  oldbox = box;
	  oldx2 = (int)x + x1 + 0.0001; /* mhmp 006.10.98 pb d'arrondi */
	  oldy2 = (int)y + y1;
	  newx2 = oldx2 +1;
	  newy2 = oldy2;
	  xfin = x2 ;
	  xlimit = x2 ;
	  xdeb = x1 ;
	  frameDNA () ;
	  oldx3 = oldx ;
	  oldx4 = oldx2 ;
	  fMapLeftDrag3 (oldx2, oldy2) ;
	  graphRegister (LEFT_DRAG, (GraphFunc)fMapLeftDrag3) ;
	  graphRegister (LEFT_UP, (GraphFunc)fMapLeftUp2) ;
	}
      else if ( BOXSEG (box)->type == TRANS_SEQ ||
		BOXSEG (box)->type == TRANS_SEQ_UP)
	{
	  if (look->flag & FLAG_REVERSE)
	    {
	      xxx = COORD (look, BOXSEG (box)->x1) - 3* (int)x ;
	      if (xxx < array (look->maxDna,arrayMax (look->maxDna) - 1 , int))
		return ;
	      for (ii=1; ii<arrayMax (look->maxDna); ii++)
		if (xxx > array (look->maxDna, ii, int))
		  {
		    if (xxx > array (look->minDna, ii-1, int))
		      return ;
		    break ;
		  }
	    }
	  else
	    {
	      xxx = COORD (look, BOXSEG (box)->x1) + 3* (int)x ;
	      if (xxx > array (look->maxDna,arrayMax (look->maxDna) - 1 , int))
		return ;
	      for (ii=1; ii<arrayMax (look->maxDna); ii++)
		if (xxx < array (look->maxDna, ii, int))
		  {
		    if (xxx < array (look->minDna, ii-1, int))
		      return ;
		    break ;
		  }
	    }
	  graphBoxDim (box,&x1, &y1, &x2, &y2) ;
	  graphColor (BLACK) ;
	  oldx = x ;
	  oldy = y;
	  oldlook = look ;
	  oldbox = box;
	  oldx2 = (int)x + x1 + 0.0001; /* mhmp 006.10.98 pb d'arrondi */
	  oldy2 = (int)y + y1;
	  newx2 = oldx2 +1;
	  newy2 = oldy2;
	  xfin = x2 ;
	  xlimit = x2 ;
	  xdeb = x1 ;
	  frameDNA () ;
	  oldx3 = oldx ;
	  oldx4 = oldx2 ;
	  fMapLeftDrag2 (oldx2, oldy2) ;
	  graphRegister (LEFT_DRAG, (GraphFunc)fMapLeftDrag2) ;
	  graphRegister (LEFT_UP, (GraphFunc)fMapLeftUp2) ;
	}
      else			/* not active and not sequence */
	fMapSelectBox (look, box, x, y) ;
    }
}

/***********************************/

static void fMapKbd (int k)
{
  FMAPLOOKGET ("fMapKbd") ;

   if (look != selectedfMap)
    fMapSelect (look) ; 
   
   switch (k)
     {
    case END_KEY : /* reverse complement */
      fMapReverseComplement () ;
      break ;

     default: /* do something with other keys here */
       mapKbdScroom (k) ;
       break ;
     }

   return ;
}

#ifdef JUNK
 this is totally useless in the context of fmap
 i replace that by mapcontrol.c:mapKbdSroom () 

  int  box ;

  if (!look->activeBox)
    return ;

  box = look->activeBox ;
  switch (k)
    {
    case UP_KEY :
      if (box > look->minLiveBox)
	--box ;
      break ;
    case DOWN_KEY :
      if (box < arrayMax (look->boxIndex) - 1)
	++box ;
      break ;

    }
  fMapSelectBox (look, box, 0, 0) ;
  if (look->activeBox != box)
    fMapFollow (look,0.,0.) ;
#endif

/************** little function to act as a callback ***********/

static void fMapRefresh (void)
{
  FMAPLOOKGET ("fMapRefresh") ;

  fMapDraw (look, 0) ;
}

/*********************************************************/
/*********************************************************/

static LOOK fMapCreateLook (KEY seq, int start, int stop, int origin, BOOL doReverse, KEY view, int my_rView)
{
  int  x1, x2 ;
  int  startOrig = start ;
  LOOK look ;
  KEY  parent = 0, dummy = 0 ; 
  KEY  seqOrig = seq ;
  OBJ  Seq ;
  static Array s_children = 0 ;
  BOOL reverse = FALSE ;
  int originShift = origin ; 

  if (! (Seq = bsCreate (seq)))
    return 0 ;
  fMap_rView = my_rView ; /* nasty glogal sorry */
	/* recurse up - done already if from fMapDisplay, but harmless */
  while (Seq)
    {
      if (bsFindTag (Seq, str2tag ("SMAP")) && /* new tag2 recursive system */
	  !bsFindTag (Seq, str2tag ("Splicing")) &&
	  bsFindTag (Seq, str2tag ("S_Parent")) &&
	  bsGetKeyTags (Seq, _bsRight, &dummy) &&
	  bsGetKey (Seq, _bsRight, &parent))
	{ 
	  int i ;
	  BSunit *uu ;
	  
	  bsDestroy (Seq) ;
	  /* find the child in the parent */
	  s_children = arrayReCreate (s_children, 128, BSunit) ;
	  if ((Seq = bsCreate (parent)) &&
	      bsGetArray (Seq, str2tag ("S_Child"), s_children, 4))
	    {
	      for (uu = 0, i = 0 ; i < arrayMax (s_children) ; i += 4)
		{
		  uu = arrp (s_children, i, BSunit) ;
		  if (uu[1].k == seq) break ;
		  uu = 0 ;
		}
	      if (uu)
		{
		  x1 = uu[2].i ; x2 = uu[3].i ;
		  /* finalize */	
		  if (x1 > x2) 
		    reverse = !reverse ;       /* mieg */
		  if (x1 < x2)
		    if (start == 1 && !stop)  /* mieg */
		      { start = x1 ; stop = x2 ; }
		    else
		      { start += x1-1 ; stop += x1-1 ; }
		  else
		    if (start == 1 && !stop)  /* mieg */
		      { start = x2 ; stop = x1 ; }
		    else
		      { int tmp = start ;
		      start = x1 - stop + 1 ; stop = x1 - tmp + 1 ; 
		      }
		  seq = parent ;
		  continue ;
		}
	    }
	}
      if (bsGetKey (Seq, _Source, &parent))
	{ 
	  bsDestroy (Seq) ; 
	  /* find the child in the parent */
	  if (! (Seq = bsCreate (parent)) ||
	      !bsFindKey (Seq, _Subsequence, seq) ||
	      !bsGetData (Seq, _bsRight, _Int, &x1) ||
	      !bsGetData (Seq, _bsRight, _Int, &x2))
	    break ;/* stop the recursion */
	  /* finalize */
	  if (x1 > x2) 
	    reverse = !reverse ;       /* mieg */
	  if (x1 < x2)
	    if (start == 1 && !stop)   /* mieg */
	      { start = x1 ; stop = x2 ; }
	    else
	      { start += x1-1 ; stop += x1-1 ; }
	  else
	    if (start == 1 && !stop)  /* mieg */
	      { start = x2 ; stop = x1 ; }
	    else
	      { int tmp = start ;
	      start = x1 - stop + 1 ; stop = x1 - tmp + 1 ; 
	      }
	  seq = parent ;
	  continue ;
	}
      else
	break ;
    }

  if (Seq)
    bsDestroy (Seq) ;

  look = (LOOK) messalloc (sizeof (struct LookStruct)) ;
  look->handle = handleCreate () ;
  look->segsHandle = handleHandleCreate (look->handle) ;
  look->magic = &FMAPLOOK_MAGIC ;
  look->seqKey = seq ;
  look->seqOrig = seqOrig ;
  look->view = view ;
  
  x1 = start ; x2 = stop ;
  if (! (look->dna = fMapFindDNA (seq, &start, &stop)))
    goto abort ;
  look->zoneMin = x1 - start ;
  look->zoneMax = look->zoneMin + (x2 - x1) + 1 ;

  if (doReverse && reverse)
    look->origin = look->zoneMax - startOrig + 2 - originShift + 1 ;
  else
    look->origin = look->zoneMin - startOrig + 1 + originShift - 1 ;

  look->colors = arrayHandleCreate (arrayMax (look->dna), char, look->handle) ;
  arrayMax (look->colors) = arrayMax (look->dna) ;
  fMapClearDNA (look) ;		/* clear colors array to _WHITE */

  look->probeAss = assHandleCreate (look->handle) ;
  look->chosen = assHandleCreate (look->handle) ;
  look->antiChosen = assHandleCreate (look->handle) ;
  look->length = arrayMax (look->dna) ;
  look->fullLength = sequenceLength (seq) ;
  look->start = start-1 ;
  look->stop = stop-1 ;
  look->featDict = dictHandleCreate (1000, look->handle) ;
  dictAdd (look->featDict, "null", 0) ; /* so real entries are non-zero */
  dictAdd (look->featDict, "null1", 0) ; /* so real entries are non-one */
  look->homolInfo = arrayHandleCreate (256, HOMOLINFO, look->handle) ;
  look->seqInfo = arrayHandleCreate (32, SEQINFO, look->handle) ;
  look->homolBury = bitSetCreate (256, look->handle) ;
  look->homolFromTable = bitSetCreate (256, look->handle) ;
  look->visibleIndices = arrayHandleCreate (64, int, look->handle) ;

  look->map = mapCreate (fMapRefresh) ; 	/* must create map before convert, because it adds columns */
  fMapDefaultColumns (look->map) ;

  if (!fMapConvert (look, FALSE))
    goto abort ;

  look->boxIndex = arrayHandleCreate (64, int, look->handle) ;
  look->activeBox = 0 ; seqDragBox = 0 ;
  look->flag = 0 ;
  look->flag |= FLAG_AUTO_WIDTH ;  /* I honestly cannot understand why you would take that away */

  if (!isGifDisplay && class (look->seqKey) == _VmRNA)
    look->translateGene = 1 ;

  if (doReverse && reverse)
    fMapRC (look) ;

  return look ;

abort :
  messfree (look->handle) ;
  messfree (look) ;
  return 0 ;
} /* fMapCreateLook */

/************************************************************/

BOOL fMapDisplay (KEY key, KEY from, BOOL isOldGraph)
     /* this is being used a DisplayFunc called by display () */
{
  int start, stop ;
  LOOK look ;
  KEY  seq = 0 ;
  OBJ  obj ;
  BOOL shiftOrigin = FALSE ;

  fMapInitialise () ;
 
  if (class (key) == _VSequence || keyFindTag (key, _DNA) ||
      keyFindTag (key, str2tag ("SMAP")))
    seq = key ;
  else if (class (key) == _VDNA)
    { seq = 0 ; dnaReClass (key, &seq) ;}
  else if ((obj = bsCreate (key)))
    { 
      if ((!bsGetKey (obj, _Sequence, &seq)) && 
	  ( class (key) == _VGene || bsGetKey (obj, str2tag ("Genomic_sequence"), &seq))
	  )
	shiftOrigin = TRUE ; 
      bsDestroy (obj) ;
    }
  if (key != seq)
    { from = key ; key = seq ; }


  if (!seq || iskey (seq) != 2)	/* iskey tests a full object */
    return FALSE ;

  start = 1 ;
  stop = 0 ;			/* go to end of sequence */

  /* recurse up to find topmost link - must do here rather than
     wait for fMapCreateLook () so we can check if we have it already 
  */
  findStartStop (&seq, &start, &stop) ;

  if (isOldGraph && graphAssFind (&GRAPH2FMAPLOOK_ASSOC, &look) &&
      look->seqKey == seq &&
      arrayExists (look->dna) &&
      start-1 >= look->start &&
      stop-1 <= look->stop)
    {      /* the requested part of the fMap is already on display */
      look->seqOrig = key ;
      fMapCentre (look, key, from) ;
      if (shiftOrigin) fMapSetOriginFromKey (look, from) ;
      graphRetitle (name (from ? from : key)) ;
      fMapDraw (look, key) ;
      return TRUE ;
    }

				/* make a new look now */
  mapGifMap = (void*)1 ;
  /* so mapCreate () called by fMapCreateLook () does not attach to graph yet */
  /* attention: in  fMapCreateLook () call,  if ever we change to TRUE, and start > stop, then use -1,TRUE */
  look = fMapCreateLook (seq, start, stop, start < stop ? 1 : -1, FALSE, 0, 0) ; 
  mapGifMap = 0 ;
  if (!look)
    { display (seq, 0, TREE) ;	/* the drawing would be meaningless */
      return FALSE ;
    }

 if (isOldGraph && graphAssFind (&GRAPH2FMAPLOOK_ASSOC, 0)) /* reuse */
    {
      KEY kk = from ? from : key ;
      if (class (kk) &&
	  (class (kk) == _VSequence || class (kk) == _VTranscribed_gene || 
	   class (kk) == _VmRNA || class (kk) == _VTranscript || class (kk) == _VGene) )
	graphRetitle (name (kk)) ;
      fMapDestroy () ;
      graphClear () ;
    }
  else				/* open a new window */ 
    { if (!displayCreate (FMAP))
	return FALSE ;

      graphRetitle (name (from ? from : key)) ;
      graphHelp ("Feature-map") ;
      graphRegister (RESIZE, fMapRefresh) ;
      graphRegister (DESTROY, fMapDestroy) ;
      graphRegister (KEYBOARD, fMapKbd) ;
      graphRegister (MESSAGE_DESTROY, displayUnBlock) ;
    }


  look->seqOrig = key ;
  look->map->min = 0 ;
  look->map->max = look->length ;

  look->graph = graphActive () ;
  graphAssociate (&GRAPH2FMAPLOOK_ASSOC, look) ; /* attach look for fmap package */
  graphAssociate (&MAP2LOOK_ASSOC, look) ;       /* attach look for map package */
  mapAttachToGraph (look->map) ;                 /* attach the map itself */

  fMapcDNAInit () ;

  fMapCentre (look, key, from) ; /* centre on key if possible and choose magnification */
  if (shiftOrigin) fMapSetOriginFromKey (look, from) ;

  if (FALSE && !getenv ("ACEDB_PROJECT") && (look->zoneMax - look->zoneMin < 1000))
    { mapColSetByName ("DNA Sequence", TRUE) ;
      mapColSetByName ("Coords", TRUE) ;
    } 
  if (look->translateGene == 1) /* established for class _VmRNA in fMapCreateLook */
    {  
      KEYSET prods = queryKey (look->seqKey, "{>product ; best_product} $/ {>product}") ;
      
      if (keySetMax (prods))	
	look->translateGene = keySet (prods, 0) ;
      else
	look->translateGene = 0 ;
      /*      mapColSetByName ("DNA Sequence", 1) ; */
      mapColSetByName ("Down Gene Translation", 1) ;
    }
  fMapDraw (look, key) ;
  fMapSelect (look) ; 
  if (!look->dna)
    graphClipDraw (0, 0, 30000, 30000) ; /*mhmp 19.06.97*/

  return TRUE ;
} /* fMapDisplay */

/************************************************/
/************************************************/
    /* also called from fmapgene.c */
int fMapOrder (const void *a, const void *b)     
{
  const SEG *seg1 = (const SEG*)a, *seg2 = (const SEG*)b ;
  int diff ;
  
  diff = seg1->type - seg2->type ;
  if (diff > 0)
    return 1 ;
  else if (diff < 0)
    return -1 ;

  if ((seg1->type | 0x1) == SPLICED_cDNA_UP || (seg1->type | 0x1) == EXON_SUPPORT_UP)
    {
      if (! (seg1->source && seg2->source)) /* The transcribed_gene */
	return seg1->source ? -1 : 1 ;
      if (seg1->source < seg2->source) return -1 ;
      if (seg1->source > seg2->source) return 1 ;

      if (seg1->key < seg2->key) return -1 ;
      if (seg1->key > seg2->key) return 1 ;

      if (seg1->x1 < seg2->x1) return -1 ;
      if (seg1->x1 > seg2->x1) return 1 ;

      return 0 ;
    }
  if ((seg1->type | 0x1) == TRANSCRIBEDGENE_UP ||
      (seg1->type | 0x1) == TRANSCRIPT_UP ||
      (seg1->type | 0x1) == PMRNA_UP ||
      (seg1->type | 0x1) == MRNA_UP ||
      (seg1->type | 0x1) == MPRODUCT_UP)
    {
      if (!seg1->source || !seg2->source)
	return seg1->source ? -1 : 1 ;
      diff = seg1->source - seg2->source ; 
      if (diff)
	return diff ;
      diff = seg1->parent - seg2->parent ;
      if (diff) /* seg1->parent &&  seg2->parent) */
	return keySetAlphaOrder (&seg1->parent, &seg2->parent) ; 
      
      /* else fall thru on usual coord comparison */
      /* and see below the seg->key comparison */
    }
    
  diff = seg1->x1 - seg2->x1 ;
  if (diff > 0)
    return 1 ;
  else if (diff < 0)
    return -1 ;

  if ((seg1->type | 0x1) == MPRODUCT_UP)
    {
      diff = seg1->source - seg2->source ;
      if (diff > 0)
	return 1 ;
      else if (diff < 0)
	return -1 ;
    }

  diff = seg1->x2 - seg2->x2 ;
  if (diff > 0)
    return 1 ;
  else if (diff < 0)
    return -1 ;

  if ((seg1->type | 0x1) == PMRNA_UP ||
      (seg1->type | 0x1) == MRNA_UP)
    {
      diff = seg1->key - seg2->key ; /* so that exon is before stolen_exon */
      if (seg1->key &&  seg2->key && diff)
	return diff ;
    }
  return 0 ;
}

/*****************************************************************/

static void fMapCentre (LOOK look, KEY key, KEY from)
		/* need key and from for Calcul positioning */
{ 
  int i, redo = 0, ny, orig = -999999 ;
  SEG *seg, *seg1 ;
  BOOL fromSeq = (class (from) == _VSequence ||
		  (_VTranscribed_gene &&
		   class (from) == _VTranscribed_gene) ||
		  (_VTranscript &&
		   class (from) == _VTranscript) ||
		  (_VmRNA &&
		   class (from) == _VmRNA) ||
		  (_VGene &&
		   class (from) == _VGene) ||
		  (_VIST &&
		   class (from) == _VIST) ||
		  (_VProtein &&
		  class (from) == _VProtein)) ; /* so that protein homol get centered */

  graphFitBounds (&mapGraphWidth,&mapGraphHeight) ;
  halfGraphHeight = 0.5* (mapGraphHeight - topMargin) ;
 
  /* nx = mapGraphWidth ; */ ny = mapGraphHeight - topMargin - 5 ;

lao:
  seg = arrp (look->segs, 0 , SEG) ;
  if (key)
    for (i = 1 ; i < arrayMax (look->segs) ; ++i)
      { seg1 = arrp (look->segs, i, SEG) ;
	if (fromSeq && seg1->key == key)
	  orig = seg->x1 ;
	if ((!fromSeq && seg1->key == key) || (from && seg1->key == from))
	  { seg = seg1 ;
	    break ;
	  }
      }

  if (!redo &&   /* code to present the object on its own strand */
      (
       (_VTranscribed_gene && class (from) == _VTranscribed_gene) ||
       (_VTranscript && class (from) == _VTranscript) ||
       (_VmRNA && class (from) == _VmRNA) ||
       (_VGene && class (from) == _VGene) ||
       (_VmProduct && class (from) == _VmProduct) ||
       ( from && keyFindTag (from, str2tag ("Gene_model"))) ||
       ( !from && keyFindTag (key, str2tag ("Gene_model")))        
       )  &&
      seg->type & 0x1)
    {
      redo = 1 ;
      fMapRC (look) ;
      goto lao ;
    }

  i = KEYKEY (from)*10 ;
  if (class (from) == _VCalcul  &&  i < seg->x2 - seg->x1 + 1)
    { 
      if (seg->type == SEQUENCE_UP)
	look->map->centre = seg->x2 - i ;
      else
	look->map->centre = seg->x1 + i ;
      /* Go to rather high magnification 
	 look->map->mag = (ny-5)/ 5000.0 ; */
    }
  else
    look->map->centre = 0.5 * (seg->x1 + seg->x2) ;
  
  look->map->mag = ny/ (seg->x2 - seg->x1 + 1.0) ;
  
  look->origin = orig != -999999 ? orig : seg->x1 ;
  setZone (look, seg->x1, seg->x2 + 1) ;
  if (fromSeq)
    look->map->mag /= 1.3 ;
  if (look->map->mag <= 0)  /* safeguard */
    look->map->mag = 0.05 ;
}

/*****************************************************************/

static void processMethod (LOOK look, KEY key, float offset,
		   void (*func) (LOOK look, float *offset))
{ 
  METHOD *meth = method (look->view, key) ;
  float p ;

  if (meth->priority)
    p = meth->priority ;
  else if (meth->flags & METHOD_FRAME_SENSITIVE)
    p = 4 + offset ;
  else
    p = 5 + offset ;
  
  mapInsertCol (look->map, p, TRUE, name (key), func) ;
  
  if (meth->flags & METHOD_SHOW_UP_STRAND)
    mapInsertCol (look->map, -p, TRUE, messprintf ("-%s", name (key)), func) ;
}

void fMapProcessMethods (LOOK look)
{
  int i ;
  KEY key ;
  
  if (!look->map)
    return ;
  
  if (methodInfo)		/* first unset METHOD_DONE */
    for (key = 0 ; lexNext (_VMethod, &key) ;)
      method (look->view, key)->flags &= ~METHOD_DONE ;
  
  for (i = 0 ; i < arrayMax (look->segs) ; ++i)
    {
      key = arrp (look->segs, i, SEG)->key ;
      if (class (key) == _VMethod && methodAdd (look->view, key))
	switch (arrp (look->segs, i, SEG)->type)
	  { 
	  case FEATURE: case FEATURE_UP: 
	    processMethod (look, key, 0.7, fMapShowFeature) ; 
	    break ;
	  case SPLICE5: case SPLICE5_UP:
	  case SPLICE3: case SPLICE3_UP:
	    processMethod (look, key, 0.9, fMapShowSplices) ; 
	    break ;
	  case TRANSCRIBEDGENE: case TRANSCRIBEDGENE_UP:
	    processMethod (look, key, 3.1, fMapcDNAShowTranscribedgene) ; 
	    break ;
	  case TRANSCRIPT: case TRANSCRIPT_UP:
	    processMethod (look, key, 3.2, fMapcDNAShowTranscript) ; 
	    break ;
	  case PMRNA: case PMRNA_UP:
	    processMethod (look, key, 1.9, fMapcDNAShowPMrna) ; 
	    break ;
	  case MRNA: case MRNA_UP:
	    processMethod (look, key, 3.2, fMapcDNAShowMrna) ; 
	    break ;
	  case MSOLEXA: case MSOLEXA_UP:
	    processMethod (look, key, 3.3, fMapcDNAShowSolexa) ; 
	    break ;
	  case MGENES: case MGENES_UP:
	    processMethod (look, key, 3.05, fMapcDNAShowGenes) ; 
	    break ;
	  case MRNAI: case MRNAI_UP:
	    processMethod (look, key, 3.22, fMapcDNAShowRNAi) ; 
	    break ;
	  case PROBE: case PROBE_UP:
	    processMethod (look, key, 3.25, fMapcDNAShowProbe) ; 
	    break ;
	  case MOST: case MOST_UP:
	    processMethod (look, key, 3.27, fMapcDNAShowOST) ; 
	    break ;
	  case MPRODUCT: case MPRODUCT_UP:
	    processMethod (look, key, 3.08, fMapcDNAShowMProduct) ; 
	    break ;
	  case EXON_SUPPORT: case EXON_SUPPORT_UP:
	    processMethod (look, key, 3.28, fMapcDNAShowExonSupport) ; 
	  case SPLICED_cDNA: case SPLICED_cDNA_UP:
	    processMethod (look, key, 3.30, fMapcDNAShowSplicedcDNA) ; 
	  case SPLICED_cDNA_DECORATE: case SPLICED_cDNA_DECORATE_UP:
	    processMethod (look, key, 3.30, fMapcDNADecorateSplicedcDNA) ; 
	    break ;
	  case CDNA_GENE_NAME: case CDNA_GENE_NAME_UP:
	    processMethod (look, key, 7.0, fMapcDNAGeneName) ; 
	    break ;
	  case ATG: case ATG_UP:
	    break ;		/* do nothing yet */
	  default:
	    messerror ("Method key %s for unknown seg type %s", name (key), 
		       fMapSegTypeName[arrp (look->segs, i, SEG)->type]) ;
	  }
    }
  /* mieg mars 2001, must come after otherwise all the SMAP->methods get registered as fMapShowSequence */
  for (i = 0 ; i < arrayMax (look->homolInfo) ; ++i)
    {
      key = arrp (look->homolInfo, i, HOMOLINFO)->method ;
      if (methodAdd (look->view, key))
	processMethod (look, key, 0.5, fMapShowHomol) ;
    }
  
  for (i = 0 ; i < arrayMax (look->seqInfo) ; ++i)
    {
      key = arrp (look->seqInfo, i, SEQINFO)->method ;
      if (key && methodAdd (look->view, key))
	processMethod (look, key, -3.0, fMapShowSequence) ;
    }
  

}

/*****************************/

static BOOL pcDown ;
static int pcStart, pcEnd ;
static Array pcArray = 0 ;

static void posConvertInit (SEG *seg, Array exons, BOOL isDown)
{
  int i ;

  pcDown = isDown ;
  pcStart = seg->x1 - 1 ;
  pcEnd = seg->x2 + 1 ;
  pcArray = arrayReCreate (pcArray,20,int) ;
  for (i = 0 ; i + 2 < arrayMax (exons) ; i+=2)
    { array (pcArray,i,int) = arr (exons,i+1,BSunit).i ;
      array (pcArray,i+1,int) = 
	arr (exons,i+2,BSunit).i - arr (exons,i+1,BSunit).i - 1 ;
    }
}

static void posConvert (SEG *seg, int pos1, int pos2)
{
  int i ;

  for (i = 0 ; i < arrayMax (pcArray) ; i += 2)
    if (pos1 > arr (pcArray, i, int))
      pos1 += arr (pcArray, i+1, int) ;
    else
      break ;

  for (i = 0 ; i < arrayMax (pcArray) ; i += 2)
    if (pos2 > arr (pcArray, i, int))
      pos2 += arr (pcArray, i+1, int) ;
    else
      break ;
	
  if (pcDown)
    { seg->x1 = pos1 + pcStart  ;
      seg->x2 = pos2 + pcStart  ;
    }
  else
    { seg->x1 = pcEnd - pos2  ;
      seg->x2 = pcEnd - pos1  ;
    }
}

static void addOldSegs (LOOK look, Array segs, Array oldSegs)
{			/* selective for CALCULATED segs */
  SEG *seg, *oldMaster, *newMaster ;
  int i, j, length, offset ;

  oldMaster = arrp (oldSegs,0,SEG) ;
  newMaster = arrp (segs,0,SEG) ;
  if (newMaster->key != oldMaster->key)
    return ;

  length = newMaster->x2 - newMaster->x1 + 1 ;
  offset = newMaster->x1 - oldMaster->x1 ;

  j = arrayMax (segs) ;
  for (i = 1 ; i < arrayMax (oldSegs) ; ++i)
    { seg = arrp (oldSegs, i, SEG) ;
      switch (seg->type)
	{
	case FEATURE: case FEATURE_UP:
	case ATG: case ATG_UP:
	case SPLICE5: case SPLICE5_UP:
	case SPLICE3: case SPLICE3_UP:
	  if (class (seg->key) != _VMethod)
	    messcrash ("Non-method key in addOldSegs") ;
	  if (method (look->view, seg->key)->flags & METHOD_CALCULATED)
	    { seg->x1 -= offset ;  
	      seg->x2 -= offset ;
	      if (seg->x1 >= 0 && seg->x2 < length)
		array (segs, j++, SEG) = *seg ;
	    }
	  break ;
	default:		/* needed for picky compiler */
	  break ;
	}
    }
}

static KEY defaultSubsequenceKey (char* name, int colour, 
				  float priority, BOOL isStrand)
{
  KEY key ;
  OBJ obj ;
  
  lexaddkey (name, &key, _VMethod) ;  /* mieg, dec 27 */
  if (iskey (key) != 2) /* i.e. no object */
    if ((obj = bsUpdate (key)))	/* make it */
      { float width = 2 ;
        KEY colourTag = _WHITE + colour ;

        if (bsAddTag (obj, _Colour))
	  { bsPushObj (obj) ;
	    bsAddTag (obj, colourTag) ;
	    bsGoto (obj, 0) ;
	  }
	bsAddData (obj, str2tag ("Width"), _Float, &width) ;
	bsAddData (obj, str2tag ("Right_priority"), _Float, &priority) ;
	if (isStrand) bsAddTag (obj, str2tag ("Show_up_strand")) ;
	bsSave (obj) ;
      }
  return key ;
}

BOOL fMapConvert (LOOK look, BOOL force)
{
  KEY	key, parent, view ;
  Array segs = 0 , oldSegs = 0 ;
  BSunit *u ;
  SEG	*seg, *seg1, *seg2 ;
  OBJ	obj ;
  int	i, nsegs = 0, iseg, iseg1, pos1, pos2, tmp ;
  BOOL	isDown, isExons ;
  AC_HANDLE oldSegsHandle ;
  SEQINFO *sinf ;
  KEY _PolyA_found = str2tag("PolyA_found") ;
  KEY M_Transposon, M_EXON_SUPPORT, M_SPLICED_cDNA, M_SPLICED_cDNA_DECORATE, M_GENE_NAME, M_RNA ;
  KEY M_TRANSCRIBEDGENE, M_TRANSCRIPT, M_mRNA, M_PmRNA, M_mProduct, M_GENES, M_RNAI, M_SOLEXA, M_OST, M_Pseudogene, M_Coding ;
  KEY _Best_product = str2tag ("Best_product") ;
  Array units = 0 ;
  static Array s_children = 0 ;
  BOOL isHideNonGold = look->view ? FALSE : prefValue ("HIDE_NON_GOLD_TG") ; /* view dominates prefs */
  static int nErru = 0 ;

  chrono ("fMapConvert") ;
  
  fMapInitialise () ;
  /* outside of gif display, fMapSelectView may be random */
  view = look->view ;

  M_Coding = defaultSubsequenceKey ("Coding", BLUE, 2.0, TRUE) ;
  M_RNA = defaultSubsequenceKey ("RNA", GREEN, 1.7, TRUE) ;
  M_Pseudogene = defaultSubsequenceKey ("Pseudogene", DARKGRAY, 1.8, TRUE) ;  /* mieg dec 27. was wrong arg order */
  M_Transposon = defaultSubsequenceKey ("Transposon", DARKGRAY, 2.05, FALSE) ;
  M_TRANSCRIBEDGENE = defaultSubsequenceKey ("TRANSCRIBEDGENE", DARKGRAY, 3.1, FALSE) ;
  M_mRNA = defaultSubsequenceKey ("mRNA_class", DARKGRAY, 3.2, TRUE) ;
  M_PmRNA = defaultSubsequenceKey ("PmRNA_class", DARKGRAY, 1.902, TRUE) ;
  M_GENES = defaultSubsequenceKey ("Genes", DARKGRAY, 3.0, TRUE) ;
  M_RNAI = defaultSubsequenceKey ("RNAi", DARKGRAY, 3.02, TRUE) ;
  M_SOLEXA = defaultSubsequenceKey ("SOLEXA", DARKGRAY, 3.3, TRUE) ;
  M_OST = defaultSubsequenceKey ("OST", DARKGRAY, 3.06, TRUE) ;
  M_mProduct = defaultSubsequenceKey ("mProduct", DARKGRAY, 3.3, TRUE) ;
  M_TRANSCRIPT = defaultSubsequenceKey ("TRANSCRIPT", DARKGRAY, 3.24, FALSE) ;
  M_EXON_SUPPORT = defaultSubsequenceKey ("EXON_SUPPORT", DARKGRAY, 3.35, TRUE) ;
  M_SPLICED_cDNA = defaultSubsequenceKey ("SPLICED_cDNA", DARKGRAY, 3.4, FALSE) ;
  M_SPLICED_cDNA_DECORATE = defaultSubsequenceKey ("SPLICED_cDNA_DECORATE", DARKGRAY, 3.4, FALSE) ;
  M_GENE_NAME = defaultSubsequenceKey ("GENE_NAME", DARKGRAY, 10.0, FALSE) ;

/* Need POS_TO_SEG1 since posConvert only valid after posConvertInit */
#define POS_TO_SEG1	if (isDown) \
                          { seg1->x1 = seg->x1 + pos1 - 1 ; \
			    seg1->x2 = seg->x1 + pos2 - 1 ; \
			  } \
                        else \
			  { seg1->x1 = seg->x2 - pos2 + 1 ; \
			    seg1->x2 = seg->x2 - pos1 + 1 ; \
			  }

  if (!look || !look->seqKey)
    messcrash ("No sequence for fMapConvert") ;

  if (look->pleaseRecompute > 0) /* -1 should be untouched */
    look->pleaseRecompute = FALSE ;

  if (! (obj = bsCreate (look->seqKey)))
    { messout ("No data for sequence %s", name (look->seqKey)) ;
      chronoReturn () ;
      return 0 ;
    }

  messStatus ("Recalculating DNA window") ;

  oldSegs = look->segs ; oldSegsHandle = look->segsHandle ;
  look->segsHandle = handleHandleCreate (look->handle) ;
  segs = look->segs = 
    arrayHandleCreate (oldSegs ? arrayMax (oldSegs) : 128, SEG, look->segsHandle) ;
  units = arrayReCreate (units, 256, BSunit) ;
  look->homolInfo = arrayReCreate (look->homolInfo, 256, HOMOLINFO) ;
  look->seqInfo = arrayReCreate (look->seqInfo, 32, SEQINFO) ;
  look->homolBury = bitSetReCreate (look->homolBury, 256) ;
  look->homolFromTable = bitSetReCreate (look->homolFromTable, 256) ;
  assClear (look->probeAss) ;
				/* first seg is master */
  seg1 = arrayp (segs,nsegs++,SEG) ;
  seg1->parent = 0 ;		/* important */
  seg1->source = 0 ;
  seg1->key = look->seqKey ;
  seg1->x1 = look->start ;
  seg1->x2 = look->stop ;
  seg1->type = MASTER ;
/*
  printf ("Converting %s: %d %d\n", name (seg1->key), seg1->x1, seg1->x2) ;
*/
				/* second is self */
  seg2 = arrayp (segs,nsegs++,SEG) ;
  seg2->key = seg2->parent = look->seqKey ;
  seg2->source = 0 ;
  seg2->x1 = -look->start ;
  seg2->x2 = look->fullLength - 1 - look->start ;
  seg2->type = SEQUENCE ;
  seg2->data.i = arrayMax (look->seqInfo) ;
  sinf = arrayp (look->seqInfo, seg2->data.i, SEQINFO) ;
  bsDestroy (obj) ;
  oligoRepeatDestroy (look->oligoRepeats) ;

  if (class (look->seqKey) == _VmRNA && keyFindTag (look->seqKey, str2tag("Solexa")))
    processMethod (look, M_SOLEXA, 3.3, fMapcDNAShowSolexa) ; 

  for (iseg = 1 ; iseg < nsegs ; ++iseg)
/* recurses through all sequences since they get added in turn */
    { seg = arrp (segs,iseg,SEG) ;
      if (seg->type == SEQUENCE || seg->type == TRANSCRIBEDGENE || seg->type == PMRNA || seg->type == MRNA  || 
	  seg->type == MPRODUCT || seg->type == TRANSCRIPT)
	isDown = TRUE ;
      else if (seg->type == SEQUENCE_UP || seg->type == TRANSCRIBEDGENE_UP || 
	       seg->type == MPRODUCT_UP || seg->type == PMRNA_UP || seg->type == MRNA_UP ||
	       seg->type == TRANSCRIPT_UP )
	isDown = FALSE ;
      else
	continue ;		/* only interested in sequences */
      if (! (obj = bsCreate (seg->key)))
	continue ;

      if (graphInterruptCalled ())
	break ;
				/* first create subsequences */
      /* new tag2 system */
      if (bsFindTag (obj, str2tag ("S_Child"))) 
	s_children = arrayReCreate (s_children, 128, BSunit) ;
      if (bsFindTag (obj, str2tag ("S_Child")) &&
	  (iseg < 2 || !bsFindTag (obj, str2tag ("Splicing"))) &&   /* mieg mars 2001, block recursions */
	  bsGetArray (obj, str2tag ("S_Child"), s_children, 4))
	{ 
	  for (i = 0 ; i < arrayMax (s_children) ; i += 4)
	    { u = arrp (s_children,i,BSunit) ;
	      if (!u[1].k)
		continue ;
	      if (look->flag &  FLAG_HIDE_SMALL_CONTIGS)
		{ pos1 = u[2].i ; pos2 = u[3].i ;  /* jump if smaller than 1 kb */
		  if ((pos1 - pos2) * ( pos1 - pos2) < 1000000)
		    continue ;
		}  

	      if (nErru++ < 5 && !u[2].i && !u[3].i)
		{ messerror ("Coords of subsequence %s missing in %s",
			     name (u[1].k), name (seg->key)) ;

		}
	      else 
		{
		  pos1 = u[2].i ; pos2 = u[3].i ;
		  if (class (u[1].k) == _VmRNA && fmapView (view, _Fmap_mRNA))
		    {     
		      KEY tmpType = 0, tg, gene ;
		      
		      iseg1 = nsegs ; 
		      if (keyFindTag (u[1].k, str2tag ("From_prediction")))
			{
			  if (isDown)
			    tmpType = (pos2 > pos1) ? PMRNA : PMRNA_UP ;
			  else
			    tmpType = (pos1 > pos2) ? PMRNA : PMRNA_UP ; 
			}
		      else
			{
			  if (isDown)
			    tmpType = (pos2 > pos1) ? MRNA : MRNA_UP ;
			  else
			    tmpType = (pos1 > pos2) ? MRNA : MRNA_UP ; 
			}
		      if ((tg = keyGetKey (u[1].k, _From_gene)) &&
			  (gene = keyGetKey (tg, _Gene)) &&
			  keyFindTag (gene, str2tag ("Cloud_gene")) &&
			  !fmapView (view, _Fmap_Cloud_gene)) ;
		      else if (fmapView (view, (fMap_rView == (tmpType & 0x1)) ? _Fmap_mRNA_Down : _Fmap_mRNA_Up))
			{
			  seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
			  seg1->key = u[1].k ;
			  seg1->parent = seg1->key ;
			  seg1->type = tmpType ;

			  seg1->data.i = 0 ; /* arrayMax (look->seqInfo) ; */
			  /* sinf = arrayp (look->seqInfo, seg1->data.i, SEQINFO) ; */
			  
			  if (pos1 > pos2)
			    { tmp = pos2 ; pos2 = pos1 ; pos1 = tmp ;} 
			  POS_TO_SEG1 ;   
			  seg1->source = seg1->x1 + 1 ; /* needed for correct order */ 
			  seg1->sourceDx = seg1->x2 - seg1->x1 - 1 ;
			  if (seg1->x2 < 0 || seg1->x1 >= look->length)
			    --nsegs ;	/* ignore this mrna */ 
			  else
			    { 
			      KEYSET tmp = queryKey (seg1->key, ">PolyA_after_base") ; 
			      if (keySetMax (tmp))
				seg1->data.i = 0x1 ;
			      keySetDestroy (tmp) ;
			      tmp = queryKey (seg1->key, "Transpliced_to = SL1") ; 
			      if (keySetMax (tmp))
				seg1->data.i |= 2 ;
			      keySetDestroy (tmp) ;
			      tmp = queryKey (seg1->key, "Transpliced_to = SL* && HERE >= SL2 ") ; 
			      if (keySetMax (tmp))
				seg1->data.i |= 4 ;
			      keySetDestroy (tmp) ;
			      if (seg1->data.i >= 6)
				seg1->data.i ^= 0x0F ; /* leave bit 1, zero bit 2 and 4 , set bit 8 */
			      
			      
			      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
			      *seg1 = * (seg1-1) ;
			      if ((seg1->type | 0x1) == PMRNA_UP)
				{
				  seg1->key = M_PmRNA ;
				  seg1->type = PMRNA ;			  
				}
			      else
				{
				  seg1->key = M_mRNA ;
				  seg1->type = MRNA ;			  
				}
			    }
			}
		    }
		  else
		    {
		      iseg1 = nsegs ;
		      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
		      seg1->key = seg1->parent = u[1].k ;
		      seg1->source = seg->key ;
		      seg1->data.i = arrayMax (look->seqInfo) ;
		      sinf = arrayp (look->seqInfo, seg1->data.i, SEQINFO) ;

		      if (isDown)
			seg1->type = (pos2 > pos1) ? SEQUENCE : SEQUENCE_UP ;
		      else
			seg1->type = (pos1 > pos2) ? SEQUENCE : SEQUENCE_UP ;

		      if (pos1 > pos2)
			{ tmp = pos2 ; pos2 = pos1 ; pos1 = tmp ;}
		      POS_TO_SEG1 ;
		      if (seg1->x2 < 0 || seg1->x1 >= look->length)
			--nsegs ;	/* ignore this subsequence */  
		      else
			{ 
			  if (keyFindTag (seg1->key, _Assembled_from))
			    {
			      seg1->data.i = arrayMax (look->seqInfo) ;
			      sinf = arrayp (look->seqInfo, seg1->data.i, SEQINFO) ;   
			      sinf->flags |=  SEQ_VIRTUAL_ERRORS ;
			    }
			  if (class (seg1->key) == _VmProduct && fmapView (view, _Fmap_mProduct))
			    {
			      if (seg1->type & 0x1)
				seg1->type = MPRODUCT_UP ;
			      else
				seg1->type = MPRODUCT ;
			      
			      seg1->parent = seg1->key ; /* to allow fMapGetCDS to work */
			      seg1->source = seg1->key ;
			      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
			      seg1->key = M_mProduct ;
			      seg1->type = MPRODUCT ;
			      seg1->parent = 0 ;
			      seg1->source = 0 ;
			    }
			  else
			    if (class (seg1->key) == _VRNAi && fmapView (view, _Fmap_RNAi))
			    {
			      if (seg1->type & 0x1)
				seg1->type = MRNAI_UP ;
			      else
				seg1->type = MRNAI ;
			      
			      seg1->parent = seg1->key ; 
			      seg1->source = seg1->key ;
			      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
			      seg1->key = M_RNAI ;
			      seg1->type = MRNAI ;
			      seg1->parent = 0 ;
			      seg1->source = 0 ;
			    }
			  else
			    if (class (seg1->key) == _VOST && fmapView (view, _Fmap_OST))
			    {
			      if (seg1->type & 0x1)
				seg1->type = MOST_UP ;
			      else
				seg1->type = MOST ;
			      
			      seg1->parent = seg1->key ; 
			      seg1->source = seg1->key ;
			      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
			      seg1->key = M_OST ;
			      seg1->type = MOST ;
			      seg1->parent = 0 ;
			      seg1->source = 0 ;
			    }
			  else if (_VGene && class (seg1->key) == _VGene &&
				   fmapView (view, (fMap_rView == (seg1->type & 0x1)) ? _Fmap_Gene_Down : _Fmap_Gene_Up)
				   )
			    {
			      OBJ Gene = 0 ;
			      char *fullName = 0 ;
			      KEY gene = seg1->key ;
			      SEG *seg2 ;

			      if (seg1->type & 0x1)
				seg1->type = MGENES_UP ;
			      else
				seg1->type = MGENES ;
			      
			      seg1->parent = 0 ; /* to allow ordering later of product->constructed_from ESTs */
			      seg1->source = gene ;
			      if ( fmapView (view, _Fmap_Gene_Name) &&
				   (bIndexFind (gene, _Title) || 
				   bIndexFind (gene, str2tag ("NewName")) || 
				   bIndexFind (gene, str2tag ("LocusLink"))  ||
				   bIndexFind (gene, _Locus) ) &&
				  (Gene = bsCreate (gene)) &&
				  (bsGetText (Gene, _Title, &fullName)  ||
				   bsGetText (Gene, str2tag ("NewName"),  &fullName) ||
				   bsGetText (Gene, str2tag ("LocusLink"),  &fullName) ||
				   bsGetText (Gene, _Locus,  &fullName)))
				{
				  seg2 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
				  seg2->key = gene ;                 seg1 = arrp (segs,iseg1,SEG) ;
				  seg2->x1 = seg1->x1 ; 
				  seg2->x2 = seg1->x2 ;  
				  seg2->parent = seg1->key ; seg1->source = seg->source /* parent */ ;
				  seg2->type = (isDown) ? CDNA_GENE_NAME : CDNA_GENE_NAME_UP ;
				  if (!look->geneNamesStack)
				    {
				      look->geneNamesStack = stackCreate (1024) ;
				      pushText (look->geneNamesStack, "Toto") ; /* avoid zero */
				    }
				  seg2->data.i = stackMark (look->geneNamesStack) ;
				  pushText (look->geneNamesStack, fullName) ; 
				  seg2 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
				  seg2->key = M_GENE_NAME ;
				  seg2->type = CDNA_GENE_NAME;
				  seg2->parent = 0 ;
				  seg2->source = 0 ;
				}
			      bsDestroy (Gene) ;
			      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
			      seg1->key = M_GENES ;
			      seg1->type = MGENES ;
			      seg1->parent = 0 ;
			      seg1->source = 0 ;
			    }
			}
		    }
		}
	    }
	}


      if (bsFindTag (obj, _Subsequence) && 
	  bsFlatten (obj, 3, units))
	{ for (i = 0 ; i < arrayMax (units) ; i += 3)
	    { u = arrp (units,i,BSunit) ;
	      if (!u[0].k)
		continue ;
	      if (look->flag &  FLAG_HIDE_SMALL_CONTIGS)
		{ pos1 = u[1].i ; pos2 = u[2].i ;  /* jump if smaller than 1 kb */
		  if ((pos1 - pos2) * ( pos1 - pos2) < 1000000)
		    continue ;
		}
	      if (arrp (segs,iseg,SEG)->key == u[0].k)
		{
		  /* about to add another seg for itself, this would go round in circles
		   * until we run out of memory (fw-990511, see SANgc04317) */
		  messerror ("Error : %s is subsequence of itself", name (u[0].k));
		  continue;
		}

	      if (!u[1].i || !u[2].i)
		{ messerror ("Coords of subsequence %s missing in %s",
			     name (seg1->key), name (seg->key)) ;
		/* --nsegs ; removed, i moved nsegs++ in the else clause */
		}
	      else 
		{ 
		  int tmpType ;
		  
		  pos1 = u[1].i ; pos2 = u[2].i ;
		  if (isDown)
		    tmpType = (pos2 > pos1) ? SEQUENCE : SEQUENCE_UP ;
		  else
		    tmpType = (pos1 > pos2) ? SEQUENCE : SEQUENCE_UP ;
	
		  if (keyFindTag (u[0].k, _Source_Exons) &&
		      (!fmapView (view, (fMap_rView == (tmpType & 0x1)) ? _Fmap_Pg_Down : _Fmap_Pg_Up))
		      );
		  else
		    {
		      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
		      seg1->key = seg1->parent = u[0].k ;
		      seg1->source = seg->key ;
		      seg1->data.i = arrayMax (look->seqInfo) ;
		      sinf = arrayp (look->seqInfo, seg1->data.i, SEQINFO) ;
		      
		      seg1->type = tmpType ;
		      if (pos1 > pos2)
			{ tmp = pos2 ; pos2 = pos1 ; pos1 = tmp ;}
		      POS_TO_SEG1 ;
		      if (seg1->x2 < 0 || seg1->x1 >= look->length)
			--nsegs ;	/* ignore this subsequence */  
		      else
			{ 
			  OBJ obj1 = bsCreate (seg1->key) ; 
			  if (obj1 &&
			      bsFindTag (obj1, _Assembled_from))
			    {
			      seg1->data.i = arrayMax (look->seqInfo) ;
			      sinf = arrayp (look->seqInfo, seg1->data.i, SEQINFO) ;   
			      sinf->flags |=  SEQ_VIRTUAL_ERRORS ;
			    }
			  bsDestroy (obj1) ;
			}
		    }
		}
	    }
	}
      
      parent = seg->key ;
      int myFlag = (seg->type | 0x1) == TRANSCRIBEDGENE_UP && seg->key == seg->parent ? seg->data.i : 0 ; /* restore the SL info of the transcribed gene, this looks fragile. 2025_01_09 */
      seg->data.i = arrayMax (look->seqInfo) ;
      sinf = arrayp (look->seqInfo, seg->data.i, SEQINFO) ;
      sinf->flags = myFlag ;
      
      if (bsFindTag (obj, str2tag ("Genomic_canonical")) ||
	  /* bsFindTag (obj, str2tag ("IS_bac")) || */
	  bsFindTag (obj, str2tag ("Genomic")))
	sinf->flags |= SEQ_CANONICAL ;

				/* default methods */
      if (bsFindTag (obj, _Pseudogene))
	sinf->method = M_Pseudogene ;
      else if (bsFindTag (obj, _Transposon))
	sinf->method = M_Transposon ;
      else if (bsFindTag (obj, _tRNA) || 
	       bsFindTag (obj, _rRNA) || 
	       bsFindTag (obj, _snRNA))
	sinf->method = M_RNA ;
      else if (bsFindTag (obj, _CDS))
	sinf->method = M_Coding ;
      else
	sinf->method = 0 ;
				/* explicit method */
      if (bsGetKey (obj, str2tag ("Method"), &sinf->method))
	{
	  if (bsGetData (obj, _bsRight, _Float, &sinf->score))
	    sinf->flags |= SEQ_SCORE ;
	  sinf->method = lexAliasOf (sinf->method) ;
	}
         /* source exons */
      pos1 = 0 ;
      if ( ! ((seg->type | 0x1)  == MPRODUCT_UP) && 
	  bsFindTag (obj, _Source_Exons) && bsFlatten (obj, 3, units))
	{ 
	  KEY gap ;
	  dnaExonsSort3 (units, 3) ;
	  for (i = 0 ; i < arrayMax (units) ; i += 3)
	    { pos2 = arr (units, i, BSunit).i ;
	      if (pos1)
		{ seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
		  POS_TO_SEG1 ;
		  seg1->key = 0 ;
		  seg1->parent = seg1->source = parent ;
		  seg1->data.i = seg->data.i ; /* index in seqInfo */
		  seg1->x1++ ;  /* Because the bsTree gives the coordinates of the exon */
		  seg1->x2-- ;
		  seg1->type = (isDown) ? INTRON : INTRON_UP ;
		}
	      pos1 = pos2 ;
	      pos2 = arr (units, i+1, BSunit).i ;
	      gap = arr (units, i+2, BSunit).k ; /* mieg */
	      if (pos2)
		{ seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
		  POS_TO_SEG1 ;
		  seg1->key = 0 ;
		  seg1->parent = seg1->source = parent ;
		  seg1->data.i = seg->data.i ; /* index in seqInfo */
		  if (gap && !strcmp ("Gap", name (gap)))
		    seg1->type = (isDown) ? EXON_GAP : EXON_GAP_UP ;
		  else
		    seg1->type = (isDown) ? EXON : EXON_UP ;
		}
	      else
		messerror ("Missing pos2 for exon in %s", name (seg->key)) ;
	      pos1 = pos2 ;
	    }
	  sinf->flags |= SEQ_EXONS ;
	  isExons = TRUE ;
	}
      else
	{ arrayMax (units) = 0 ;
	  array (units, 0, BSunit).i = 1 ;
	  array (units, 1, BSunit).i = seg->x2 - seg->x1 + 1 ;
	  isExons = FALSE ;
	}

      posConvertInit (seg, units, isDown) ;

      if ( ! ((seg->type | 0x1)  == MPRODUCT_UP) && 
	   bsFindTag (obj, _CDS))
	{ seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	  seg1->key = _CDS ;
	  seg1->type = isDown ? CDS : CDS_UP ;
	  seg1->parent = seg1->source = parent ;
	  pos1 = 1 ; 
	  bsGetData (obj, _bsRight, _Int, &pos1) ;
	  if (bsGetData (obj, _bsRight, _Int, &pos2))
	    posConvert (seg1, pos1, pos2) ;
	  else
	    { seg1->x1 = seg->x1 - 1 + pos1 ;
	      seg1->x2 = seg->x2 ;
	    }
	  if (!isExons)		/* must have EXON to make CODING */
	    { seg2 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	      *seg2 = *seg ; seg2->source = parent ;
	      seg2->key = 0 ;
	      seg2->type = isDown ? EXON : EXON_UP ;
	      sinf->flags |= SEQ_EXONS ;
	    }
	}

  	 /* remarks, clones etc */
      if (bsFindTag (obj, _Visible) && 
	  bsFlatten (obj, 2, units))
	for (i = 0 ; i < arrayMax (units) ; i += 2)
	  {			/* RD 980414 checks to prevent NULL KEY in display 
				   moved here from Jean's checks in fMapShowVisible ()
				 */
	    if (!iskey (arr (units,i,BSunit).k) || class (arr (units,i,BSunit).k))
	      continue ;
				/* should be a tag */
	    if (!iskey (arr (units,i+1,BSunit).k))
	      continue ;
				/* there is a problem if the model says: Visible foo Text
				   a hack because iskey () could be true for a text pointer
				   should store these as different type, e.g. VISIBLE_TEXT
				   */
	    seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	    *seg1 = *seg ; seg1->source = parent ;
	    seg1->data.k = arr (units, i, BSunit).k ;
	    seg1->key = arr (units, i+1, BSunit).k ;
	    seg1->type = VISIBLE ;
	  }
      
  	 /* assembly tags */
      if (bsFindTag (obj, _Assembly_tags) && 
	  bsFlatten (obj, 4, units))
	for (i = 0 ; i < arrayMax (units)/4 ; ++i)
	  { int length;
	    char *string;
	    char *type_text, *text_text;
	    seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	    *seg1 = *seg ; seg1->source = parent ;
	    /* get start and end of tag */
	    pos1 = arr (units, 4*i+1, BSunit).i;
	    pos2 = arr (units, 4*i+2, BSunit).i;
	    /* Convert co-ords of pos1, pos2 and stuff into seg as x1 and x2 */
	    posConvert (seg1, pos1, pos2);
	    /* Key is the sort of assemly tag */
	    seg1->type = ASSEMBLY_TAG ;
	    seg1->key = _Assembly_tags;
	    /* Get the text */
	    length = 0;
	    type_text = arr (units, 4*i, BSunit).s;
	    if (type_text)
	      length += strlen (type_text);
	    text_text = arr (units, 4*i+3, BSunit).s;
	    if (text_text)
	      length += strlen (text_text);
	    seg1->data.s = string = halloc (length+3, look->segsHandle);
	    string[0] = '\0';

	    if (type_text)
	      strcpy (string, type_text) ;
	    strcat (string, ": ");
	    if (text_text) 
	      strcat (string, text_text);
	    seg1->parent = 0;	/* assembly tags don't need parents */
	  }
      
  	 /* primers (directed sequencing */ /* mieg */
      if (bsFindTag (obj, _Primer) &&
	  bsFlatten (obj, 1, units))
	for (i = 0 ; i < arrayMax (units) ; ++i)
	  { seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	    *seg1 = *seg ; seg1->source = parent ;
	    seg1->key = arr (units,i, BSunit).k ;
	    seg1->type = PRIMER ;
	  }

  	 /* oligo repeats primers (all repeats in genome of various lengths) */ /* mieg */
      if (bsGetKey (obj, str2tag ("OligoRepeat"), &key))
	look->oligoRepeats = oligoRepeatRegister (look->oligoRepeats, seg->x1, seg->x2, key, look->handle) ; 

  	 /* virtual sequence SCF */ /* mieg */
      if (bsFindTag (obj, _SCF_File))
	  { seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	    *seg1 = *seg ; seg1->source = parent ;
	    seg1->type = VIRTUAL_SUB_SEQUENCE_TAG ;
	    if (!isDown)
	      { if (seg->x1 > seg->x2) /* never happens i suppose */
		  { seg1->x1 = seg->x2 ; seg1->x2 = seg->x1 ;  }
		seg1->type = VIRTUAL_SUB_SEQUENCE_TAG_UP ;
	      }
	  }

         /* Transcribed_gene,  mieg */
      pos1 = 0 ;
      if ((fmapView (view, _Fmap_Tg) || fmapView (view, _Fmap_Spliced_cDNA)) &&
	  bsFindTag (obj, _Transcribed_gene) && 
	  bsFlatten (obj, 3, units))
	{ 
	  for (i = 0 ; i < arrayMax (units) ; i += 3)
	    { 
	      KEY tmpType = 0 ;
	      KEY gene ;

	      u = arrp (units,i,BSunit) ;
	      if (!u[0].k)
		continue ;
	      if (!u[1].i || !u[2].i)
		continue ;
	      if (
		  (fmapView (view, _Fmap_Tg_Just_Gold) == 2 || isHideNonGold) &&
		  !keyFindTag (u[0].k, _Gold))
		continue ;
	      pos1 = u[1].i ; pos2 = u[2].i ;
	      if (isDown)
		tmpType = (pos2 > pos1) ? TRANSCRIBEDGENE : TRANSCRIBEDGENE_UP ;
	      else
		tmpType = (pos2 > pos1) ? TRANSCRIBEDGENE_UP : TRANSCRIBEDGENE ;
	      
	      if (!fmapView (view, (fMap_rView == (tmpType & 0x1)) ? _Fmap_Tg_Down : _Fmap_Tg_Up) &&
		  !fmapView (view, (fMap_rView == (tmpType & 0x1)) ? _Fmap_Spliced_cDNA_Down : _Fmap_Spliced_cDNA_Up) 
		  )
		continue ;
	      
	      if ((gene = keyGetKey (u[0].k, _Gene)) &&
		  keyFindTag (gene, str2tag ("Cloud_gene")) &&
		  !fmapView (view, _Fmap_Cloud_gene))
		continue ;
	      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	      seg1->key = seg1->parent = u[0].k ;
	      seg1->data.i = 0 ;
	      seg1->type = tmpType ;
	      
	      if (pos1 > pos2)
		{ tmp = pos2 ; pos2 = pos1 ; pos1 = tmp ;}
	      POS_TO_SEG1 ;
	      seg1->source = seg1->x1 + 1 ;
	      seg1->sourceDx = seg1->x2 - seg1->x1 - 1 ;
	      if (seg1->x2 < 0 || seg1->x1 >= look->length)
		--nsegs ;	/* ignore this transcribedgene */  
	      else
		{ 
		  KEYSET tmp = queryKey (seg1->key, "Transpliced_to = SL1") ; 
		  if (keySetMax (tmp))
		    seg1->data.i |= 2 ;
		  keySetDestroy (tmp) ;
		  tmp = queryKey (seg1->key, "Transpliced_to = SL2") ; 
		  if (keySetMax (tmp))
		    seg1->data.i |= 4 ;
		  keySetDestroy (tmp) ;
		  if (seg1->data.i >= 6)
		    seg1->data.i ^= 0x0F ; /* leave bit 1, zero bit 2 and 4 , set bit 8 */
		}
	    }
	}
      
      /* transcribed_gene pieces,  mieg */
      pos1 = 0 ;
      if ((seg->type | 0x1) == TRANSCRIBEDGENE_UP &&
	  bsFindTag (obj, str2tag ("Splicing")) &&
	  bsFlatten (obj, 4, units))
	{ 
	  seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	  seg1->key = M_TRANSCRIBEDGENE ;
	  seg1->type = (isDown) ? TRANSCRIBEDGENE : TRANSCRIBEDGENE_UP ;
	  seg1->source = seg->x1 + 1 ; /* needed for correct order, avoid zero */
	  seg1->sourceDx = seg->x2 - seg->x1 - 2 ;
	  for (i = 0 ; i < arrayMax (units) ; i += 4)
	    { 
	      char *cp ;
	      int mycol = 0 ;
	      int cover = 0 ;
	      KEY intr ;
	      KEY chrom ;
	      int tg1 , tg2 ;
	      if (arrayMax (units) < i + 3)
		continue ;
	      pos1 = arr (units, i , BSunit).i ;
	      pos2 = arr (units, i + 1, BSunit).i ;
	      if (bsGetKey (obj, _IntMap, &chrom) &&
		  bsGetData (obj, _bsRight, _Int, &tg1) &&
		  bsGetData (obj, _bsRight, _Int, &tg2)
		  )
		{
		  int gpos1, gpos2 ;
		  if (tg1 < tg2)
		    { gpos1 = tg1 + pos1 - 1 ; gpos2 = tg1 + pos2 - 1 ; }
		  else
		    { gpos1 = tg1 - pos1 + 1 ; gpos2 = tg1 - pos2 + 1 ; }
		  if (_VIntron && chrom && 
		      lexword2key (messprintf ("%s__%d_%d", name(chrom), gpos1, gpos2), &intr, _VIntron)
		      )
		    {
		      OBJ Intr = bsCreate (intr) ;
		      if (0 && /* 2023_04_25  obolete color pattern */
			  keyFindTag (intr, _Homol)
			  )
			mycol = RED ;
		      if (Intr)
			{
			  int x = 0 ;
			  if (bsGetData (Intr, _RNA_Seq, _Int, &x) && x >= 1)
			    {
			      float z = log(1+x) ;
			      float z1 = log(2000.0)/8 ;
			      int iCol = .5 + z/z1 ;
			      if (iCol > 8) 
				iCol = 8 ;
			      if (iCol > 0)
				mycol = iCol ;
			      cover = x ;
			    }
			  bsDestroy (Intr) ;
			}
		    }
		}
	      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	      POS_TO_SEG1 ;
	      seg1->key = arr (units, i + 2, BSunit).k ;
	      seg1->parent = seg->key ; 
	      seg1->type = (isDown) ? TRANSCRIBEDGENE : TRANSCRIBEDGENE_UP ;
	      seg1->source = seg->x1 + 1 ; /* needed for correct order, avoid zero */
	      seg1->sourceDx = seg->x2 - seg->x1 - 2 ;
	      if (seg1->key == str2tag ("Intron") ||
		   seg1->key == str2tag ("Alternative_Intron"))
		{ /* in exon case, that init is an int */
		  cp =  arr (units, i + 3, BSunit).s ;
		  if (cp && *cp)
		    {
		      if (!lexword2key (cp, & (seg1->data.k), 0))
			{
			  if (strstr (cp, "Fuzzy"))
			    lexaddkey ("Fuzzy", &seg1->data.k,0) ;
			  else
			    lexaddkey ("Other", &seg1->data.k,0) ;
			}
		    }
		  seg1->data.k = mycol ? KEYMAKE (mycol, cover) : 0 ;
		}
	    }
	}

      /***************/
         /* Transcript, i.e. transcripts direct from cosmid,  mieg */
      pos1 = 0 ;
      if (bsFindTag (obj, str2tag ("Transcripts")) && 
	  bsFlatten (obj, 3, units))
	{ for (i = 0 ; i < arrayMax (units) ; i += 3)
	    { u = arrp (units,i,BSunit) ;
	      if (!u[0].k)
		continue ;
	      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	      seg1->key = seg1->parent = u[0].k ;
	      seg1->data.i = 0 ;	 
	      if (!u[1].i || !u[2].i)
		{ /* messerror ("Coords of transcript %s missing in %s",
			     name (seg1->key), name (seg->key)) ;
			     */
		  --nsegs ;
		}
	      else 
		{ pos1 = u[1].i ; pos2 = u[2].i ;
		  if (isDown)
		    seg1->type = (pos2 > pos1) ? TRANSCRIPT : TRANSCRIPT_UP ;
		  else
		    seg1->type = (pos1 > pos2) ? TRANSCRIPT_UP : TRANSCRIPT ;
		  if (pos1 > pos2)
		    { tmp = pos2 ; pos2 = pos1 ; pos1 = tmp ;}
		  POS_TO_SEG1 ;  
		  seg1->source = seg1->x1 ; /* needed for correct order */
		  seg1->sourceDx = seg1->x2 - seg1->x1 - 1 ;
		  if (seg1->x2 < 0 || seg1->x1 >= look->length)
		    --nsegs ;	/* ignore this transcript */  
		  else
		    { 
		      KEYSET tmp = queryKey (seg1->key, ">PolyA_after_base") ; 
		      if (keySetMax (tmp))
			seg1->data.i = 0x1 ;
		      keySetDestroy (tmp) ;
		      tmp = queryKey (seg1->key, "Transpliced_to = SL1") ; 
		      if (keySetMax (tmp))
			seg1->data.i |= 2 ;
		      keySetDestroy (tmp) ;
		      tmp = queryKey (seg1->key, "Transpliced_to = SL2") ; 
		      if (keySetMax (tmp))
			seg1->data.i |= 4 ;
		      keySetDestroy (tmp) ;
		      if (seg1->data.i >= 6)
			seg1->data.i ^= 0x0F ; /* leave bit 1, zero bit 2 and 4 , set bit 8 */
		    }
		}
	    }
	}
      
         /* transcript pieces, exon part mieg */
      pos1 = 0 ;
      if ((seg->type | 0x1) == TRANSCRIPT_UP &&
	  bsFindTag (obj, str2tag ("Splicing")) &&
	  bsFlatten (obj, 4, units))
	{ 
	  seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	  seg1->key = M_TRANSCRIPT ;
	  seg1->type = (isDown) ? TRANSCRIPT : TRANSCRIPT_UP ;
	  seg1->parent = seg->key ; 
	  seg1->source = seg->source ; 
	  seg1->sourceDx = seg->x2 - seg->x1 - 1 ;
	  for (i = 0 ; i < arrayMax (units) ; i += 4)
	    { 
	      char *cp ;
	      if (arrayMax (units) < i + 3)
		continue ;
	      
	      pos1 = arr (units, i , BSunit).i ;
	      pos2 = arr (units, i + 1, BSunit).i ;
	      
	      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	      POS_TO_SEG1 ;
	      seg1->key = arr (units, i + 2, BSunit).k ;
	      seg1->parent = seg->key ; 
	      seg1->source = seg->source + 1 ; /* needed for correct order */
	      seg1->sourceDx = seg->sourceDx - 2 ;
	      seg1->type = (isDown) ? TRANSCRIPT : TRANSCRIPT_UP ;

	      if (seg1->key == str2tag ("Intron") ||
		   seg1->key == str2tag ("Alternative_Intron"))
		{ /* in exon case, that init is an int */
		  cp =  arr (units, i + 3, BSunit).s ;
		  if (cp && *cp)
		    {
		      if (!lexword2key (cp, & (seg1->data.k), 0))
			{
			   if (strstr (cp, "Fuzzy"))
			    lexaddkey ("Fuzzy", &seg1->data.k,0) ;
			   else
			     lexaddkey ("Other", &seg1->data.k,0) ;
			}
		    }
		}
	    }
	}

         /* mRNA pieces,  mieg */
      pos1 = 0 ;
      if (( (seg->type | 0x1) == PMRNA_UP || (seg->type | 0x1) == MRNA_UP) &&
	  bsGetArray (obj, str2tag ("Splicing"), units, 6))
	{ 
	  seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	  *seg1 = *seg ;
	  if ((seg->type | 0x1) == PMRNA_UP)
	    {
	      seg1->key = M_PmRNA ;
	      seg1->type = (isDown) ? PMRNA : PMRNA_UP ;
	    }
	  else
	    {
	      seg1->key = M_mRNA ;
	      seg1->type = (isDown) ? MRNA : MRNA_UP ;
	    }
	  seg1->parent = seg->key ; 

	  for (i = 0 ; i < arrayMax (units) ; i += 6)
	    { 
	      char *cp ;
	      if (arrayMax (units) < i + 5)
		continue ;
	      /* exons are treated via the Coding tag */
	      if (arr (units, i + 4, BSunit).k == str2tag ("Exon") ||
		   arr (units, i + 4, BSunit).k == str2tag ("Alternative_Exon"))
		continue ;
	      {
		pos1 = arr (units, i , BSunit).i ;
		pos2 = arr (units, i + 1, BSunit).i ;
		
		seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
		POS_TO_SEG1 ;
		seg1->key = arr (units, i + 4, BSunit).k ;
		seg1->parent = seg->key ; 
		seg1->source = seg->source + 1 ; /* needed for correct order */
		seg1->sourceDx = seg->sourceDx - 2 ;
		if ((seg->type | 0x1) == PMRNA_UP)
		  seg1->type = (isDown) ? PMRNA : PMRNA_UP ;
		else
		  seg1->type = (isDown) ? MRNA : MRNA_UP ;
		
		if (seg1->key == str2tag ("Intron") ||
		    seg1->key == str2tag ("Alternative_Intron"))
		  { /* in exon case, that init is an int */
		    cp =  arr (units, i + 5, BSunit).s ;
		    if (cp && *cp)
		      {
			if (!lexword2key (cp, & (seg1->data.k), 0))
			  {
			    if (strstr (cp, "Fuzzy"))
			      lexaddkey ("Fuzzy", &seg1->data.k,0) ;
			    else
			      lexaddkey ("Other", &seg1->data.k,0) ;
			  }
		      }
		  }
	      }
	    }
	}
         /* mRNA pieces, coding part mieg */
      pos1 = 0 ;
      if (( (seg->type | 0x1) == PMRNA_UP || (seg->type | 0x1) == MRNA_UP) &&
	  bsGetArray (obj, str2tag ("Coding"), units, 5))
	{ 
	  char bestP[1024], *bp, *cp ;
	  int ibp = 0, nbp = 0, isExon, nX = 0 ;
	  KEY prod ;
	  BSunit *uu ;

	  /* identify products to be shown in genome view */
	  memset (bestP, 0, sizeof(bestP)) ;
	  bp = bestP ;
	  if (bsGetKey (obj, _Product, &prod)) 
	    do
	      {
		if (keyFindTag (prod, _Best_product) ||
		    keyFindTag (prod, str2tag("Good_product")) ||
		    keyFindTag (prod, str2tag("Very_good_product")))
		  { *bp++ = 'A' + ibp ;nbp++ ; }
		ibp++ ;
	      } while (nbp < 320 && bsGetKey (obj, _bsDown, &prod)) ;
	  
	  /* identify their exons */
	  for (i = 0 ; i < arrayMax (units) ; i += 5)
	    {
	      uu = arrp (units, i, BSunit) ; 
	      cp = uu[4].s ;
	      isExon = 0 ;
	      if (cp)
		while (*cp)
		  {  /* we expect 5a 5b / a 5b / 3a 5b / 3a b / 3a 3b */
		    if ((*cp == '3' || *cp == '5') && (*++cp)) ;
		    else 
		      {
			int i = 1023 ;
			bp = bestP ;
			while (i-- && *bp) 
			  if (*bp++ == *cp)
			    isExon++ ;
		      }
		    cp++ ;
		  }
	      uu[4].i = isExon ;
	    }
	  	      
	  seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	  *seg1 = *seg ;
	  if ((seg->type | 0x1) == PMRNA_UP)
	    {
	      seg1->key = M_PmRNA ;
	      seg1->type = (isDown) ? PMRNA : PMRNA_UP ;
	    }
	  else
	    {
	      seg1->key = M_mRNA ;
	      seg1->type = (isDown) ? MRNA : MRNA_UP ;
	    }
	  seg1->parent = seg->key ;

	  for (i = 0 ; i < arrayMax (units) ; i += 5)
	    { 
	      if (arrayMax (units) < i + 5)
		continue ;
	      uu = arrp (units, i, BSunit) ; 
	      pos1 = uu[0].i ;
	      /* fuse segments with same number of very_good_products */
	      isExon = uu[4].i ;
	      while (i < arrayMax (units) - 5)
		if (isExon == arr (units, i + 9, BSunit).i &&
		    arr (units, i + 1, BSunit).i + 1 == arr (units, i + 5, BSunit).i)
		  i += 5 ;
		else
		  break ;
	      pos2 = arr (units, i + 1, BSunit).i ;

	      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	      POS_TO_SEG1 ;
	      seg1->key = str2tag ("Exon") ;
	      seg1->parent = seg->key ; 
	      seg1->source = seg->source + 1 ; /* needed for correct order */
	      seg1->sourceDx = seg->sourceDx - 2 ;
	      if ((seg->type | 0x1) == PMRNA_UP)
		seg1->type = (isDown) ? PMRNA : PMRNA_UP ;
	      else
		seg1->type = (isDown) ? MRNA : MRNA_UP ;
	      seg1->data.i = 0 ;
	      if (isExon)
		{		     
		     switch (nX % 7)
		       {
		       case 0: seg1->data.i = 2 ; break ;
		       case 1: seg1->data.i = 3 ; break ;
		       case 2: seg1->data.i = 6 ; break ;
		       case 3: seg1->data.i = 7 ; break ;
		       case 4: seg1->data.i = 10 ; break ;
		       case 5: seg1->data.i = 14 ; break ;
		       case 6: seg1->data.i = 15 ; break ;
		       }
		     nX++ ; /* RED YELLOW */
		}

	    }
	}

         /* mRNA pieces, Valid3p  mieg */
      pos1 = 0 ;
      if (( (seg->type | 0x1) == PMRNA_UP || (seg->type | 0x1) == MRNA_UP))
	{
	  if (bsGetArray (obj, str2tag ("Valid3p"), units, 5))
	    { 
	      BSunit *uu ; char *cp ;
	      
	      for (i = 0 ; i < arrayMax (units) ; i += 5)
		{ 
		  if (arrayMax (units) < i + 1)
		    continue ;
		  uu = arrp (units, i, BSunit) ; 
		  pos1 = uu[0].i ;
		  pos2 = pos1 ;
		  
		  seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
		  POS_TO_SEG1 ;
		  seg1->key = str2tag ("Valid3p") ;
		  seg1->parent = seg->key ; 
		  seg1->source = seg->source + 1 ; /* needed for correct order */
		  seg1->sourceDx = seg->sourceDx - 2 ;
		  if ((seg->type | 0x1) == PMRNA_UP)
		    seg1->type = (isDown) ? PMRNA : PMRNA_UP ;
		  else
		    seg1->type = (isDown) ? MRNA : MRNA_UP ;
		  seg1->data.i = uu[2].i ; /* nb of supporting clones */
		  cp = arr (units, i + 3, BSunit).s ;
		  if (cp && *cp && strcasecmp (cp, "AATAAA"))
		    seg1->data.i = - seg1->data.i ;
		  if (bsFindTag (obj, str2tag ("Aggregated_3p_clones")))
		    seg1->data.i |= (1 << 24) ;
		}
	    }
	  else
	    {
	      static Array pA = 0 ;
	      if (bsFindTag (obj, str2tag ("PolyA_found")))
		{
		  /* last mRNA exons, reclip on the polyA */
		  int dx = 0, maxA = 0 ;
		  pA = arrayReCreate (pA, 40, BSunit) ;
		  bsGetArray (obj, _PolyA_found, pA, 3) ;
		  for (int j = 0 ; j < arrayMax (pA) ; j+= 3)
		    {
		      BSunit *pp = arrp (pA, j, BSunit) ;
		      if (pp[2].i > maxA)
			{ dx = pp[0].i ; maxA = pp[2].i ; }
		    }

		  pos1 = pos2 = seg->x2 - seg->x1 + dx ;
		  seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
		  POS_TO_SEG1 ;
		  seg1->key = str2tag ("Valid3p") ;
		  seg1->parent = seg->key ; 
		  seg1->source = seg->source + 1 ; /* needed for correct order */
		  seg1->sourceDx = seg->sourceDx - 2 ;
		  
		  if ((seg->type | 0x1) == PMRNA_UP)
		    seg1->type = (isDown) ? PMRNA : PMRNA_UP ;
		  else
		    seg1->type = (isDown) ? MRNA : MRNA_UP ;
		  seg1->data.i = 10 ;
		  seg1->data.i = - seg1->data.i ;
		  seg1->data.i |= (1 << 24) ;
		}
	      else if (bsFindTag (obj, str2tag ("Aggregated_3p_clones")))
		{
		  pos1 = pos2 = seg->x2 - seg->x1 ;
		  seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
		  POS_TO_SEG1 ;
		  seg1->key = str2tag ("Valid3p") ;
		  seg1->parent = seg->key ; 
		  seg1->source = seg->source + 1 ; /* needed for correct order */
		  seg1->sourceDx = seg->sourceDx - 2 ;
		  
		  if ((seg->type | 0x1) == PMRNA_UP)
		    seg1->type = (isDown) ? PMRNA : PMRNA_UP ;
		  else
		    seg1->type = (isDown) ? MRNA : MRNA_UP ;
		  seg1->data.i = 10 ;
		  seg1->data.i = - seg1->data.i ;
		  seg1->data.i |= (1 << 24) ;
		}
	    }	    
	}

           /* mRNA pieces, Valid5p  mieg */
      pos1 = 0 ;
      
      if (( (seg->type | 0x1) == PMRNA_UP || (seg->type | 0x1) == MRNA_UP) &&
	  bsGetArray (obj, str2tag ("Valid5p"), units, 3))
	{ 
	  /* disregard a valid5p if protein was constructed from AA==1 
	   * it means that the protein was computed before we decided
	   * that the mrna is 'aggregated' or 'capped'
	   */
	  BSunit *uu ;
	  KEYSET testKs = 0 ;
	  testKs = queryKey (seg->key, ">product ; mRNA_5p_complete && best_product && !at_position_1") ; 
	  if (! keySetMax(testKs)) 
	    for (i = 0 ; i < arrayMax (units) ; i += 3)
	      { 
		if (arrayMax (units) < i + 1)
		  continue ;
		uu = arrp (units, i, BSunit) ; 
		pos1 = uu[0].i ;
		pos2 = pos1 ;
		
		seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
		POS_TO_SEG1 ;
		seg1->key = str2tag ("Valid5p") ;
		seg1->parent = seg->key ; 
		seg1->source = seg->source + 1 ; /* needed for correct order */
		seg1->sourceDx = seg->sourceDx - 2 ;
		if ((seg->type | 0x1) == PMRNA_UP)
		  seg1->type = (isDown) ? PMRNA : PMRNA_UP ;
		else
		  seg1->type = (isDown) ? MRNA : MRNA_UP ;
		seg1->data.i = uu[2].i ; /* nb of supporting clones */
		if (bsFindTag (obj, str2tag ("Transpliced_to")))
		  seg1->data.i |= (2 << 24) ;
		if (bsFindTag (obj, str2tag ("Aggregated_5p_clones")))
		  seg1->data.i |= (1 << 24) ;
	      }
	  keySetDestroy (testKs) ;
	}

       /* mProduct pieces,  mieg */
      pos1 = 0 ;
      if ((seg->type | 0x1) == MPRODUCT_UP &&
	  bsFindTag (obj, _Source_Exons))
	{
	  KEY mrna, ztype, zsubtype ;
	  OBJ Mrna = 0 ;
	  int b1, b2, ln ;
	  int frame = 1 ;
	  int z1, z2 ;
	  BOOL isUorf = bsFindTag (obj, str2tag("uORF")) && !bsFindTag (obj, str2tag("Good_product")) ;

	  if (bsGetKey (obj, _mRNA, &mrna) && (Mrna = bsCreate (mrna)))
	    {
	      if (bsFindKey (Mrna, _Product, seg->parent) &&
		  bsGetData (Mrna, _bsRight, _Int, &b1)  &&
		  bsGetData (Mrna, _bsRight, _Int, &b2) ) ;
	      else b1 = b2 = 0 ;
	      
	      bsGetArray (Mrna, str2tag ("Splicing"), units, 6) ;
	      bsGetData (obj, str2tag ("Frame"), _Int, &frame) ;
	      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	      seg1->key = M_mProduct ;
	      seg1->type = MPRODUCT ;
	      seg1->parent = seg->parent ; 
	      seg1->source = seg->key ;

	      for (i = 0 ; i < arrayMax (units) ; i += 6)
		{ 
		  if (arrayMax (units) < i + 6)
		    continue ;
		  z1 = arr (units, i + 2, BSunit).i  ;
		  z2 = arr (units, i + 3, BSunit).i  ;
		  ztype = arr (units, i + 4, BSunit).k ; 
		  ln = 0 ;
		  if (strstr (name (ztype), "Intron") &&
		      i > 5 && i <  arrayMax (units) - 6)
		    {
		      z1 = pos1 = arr (units, i - 6 + 3, BSunit).i  ;
		      z2 = pos2 = arr (units, i + 6 + 2, BSunit).i  ;
		      ln =  arr (units, i + 1, BSunit).i - arr (units, i + 0, BSunit).i + 1 ;
		      if (ln < 0 || ln >= 1<<20) ln = 0 ;
		    }
		  while (TRUE)
		    {
		      int delta = 0 ;
		      /* split the exon/UTR in mrna coordinates */
		      if (!z1 && !z2)
			break ;
		      ztype = arr (units, i + 4, BSunit).k ;
		      /* delta controls the order in which components are drawn in fMapcDNAShowMProductInFrame */
		      if (strstr (name (ztype), "Intron"))
			{ z1 = z2 = 0 ; delta = 3 ;  }
		      else if (z1 < b1 && z2 < b1)
			{ pos1 = z1 ; pos2 = z2 ; z1 = z2 = 0 ; lexaddkey ("5' UTR", &ztype, 0) ; delta = 1 ; }
		      else if (z1 > b2 && z2 > b2)
			{ pos1 = z1 ; pos2 = z2 ; z1 = z2 = 0 ; lexaddkey ("3' UTR", &ztype, 0) ; delta = 2 ; }
		      else if (z1 < b1 && z2 >= b1)
			{ pos1 = z1 ; pos2 = b1 - 1 ; z1 = b1 ; lexaddkey ("5' UTR", &ztype, 0) ; delta = 1 ; }
		      else if (z1 >= b1 && z2 > b2)
			{ pos1 = z1 ; pos2 = b2 ; z1 = b2 + 1;  delta = 5 ; }
		      else  /* if (z1 >= b1 && z2 <= b2) */
			{ pos1 = z1 ; pos2 = z2 ; z1 = z2 = 0 ; delta = 4 ;  }

		      if (isUorf && ztype == str2tag("Exon"))
			ztype = str2tag("uORF") ;

		      pos1 -= b1 -1 ; pos2 -= b1 - 1 ; /* go to Product coords */
		      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
		      POS_TO_SEG1 ;
		      seg1->key = ztype ; zsubtype = 0 ;
		      if (arr (units, i + 5, BSunit).s)
			lexaddkey (arr (units, i + 5, BSunit).s, &zsubtype, 0) ;
		      if (zsubtype > 0xfff) zsubtype = 0 ;
		      seg1->data.i = (zsubtype & 0xfff) | (ln << 12) ; ;
		      seg1->parent = seg->key ; 
		      seg1->source = 10*frame + delta ; /* needed for correct order */
		      seg1->sourceDx = b2 - b1 + 1 ;
		      seg1->type = (isDown) ?  MPRODUCT : MPRODUCT_UP ;
		    }
		}
	      if (bsFindTag (obj, str2tag ("Best_product")))
		{
		  int  delta = 7 ; char *cp ;
		  KEYSET testKs = 0 ;

		  bsGetData (obj, str2tag ("Frame"), _Int, &frame) ;
		  bsGetArray (Mrna, str2tag ("Valid3p"), units, 5) ;
		  for (i = 0 ; i < arrayMax (units) ; i += 5)
		    {
		      pos1 = pos2 = arr (units, i + 1, BSunit).i ;
		      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
		      seg1->x1 = pos1 - 1 ; seg1->x2 = pos1 - 1 ;
		      seg1->key = str2tag ("Valid3p") ;
		      seg1->data.i = arr (units, i + 2, BSunit).i ;
		      seg1->parent = seg->key ; 
		      seg1->source = 10*frame + delta ; /* needed for correct order */
		      seg1->sourceDx = b2 - b1 + 1 ;
		      seg1->type = (isDown) ?  MPRODUCT : MPRODUCT_UP ;
		      cp = arr (units, i + 3, BSunit).s ;
		      if (cp && *cp && strcasecmp (cp, "AATAAA"))
			seg1->data.i = - seg1->data.i ;
		      if (bsFindTag (Mrna, str2tag ("Aggregated_3p_clones")))
			seg1->data.i |= (1 << 24) ;
		    }
		  /* disregard a valid5p if protein was constructed from AA==1 
		   * it means that the protein was computed before we decided
		   * that the mrna is 'aggregated' or 'capped'
		   */
		  testKs = queryKey (mrna, ">product ; mRNA_5p_complete && best_product && !at_position_1") ; 
		  if (!testKs)
		    {
		      bsGetArray (Mrna, str2tag ("Valid5p"), units, 3) ;
		      for (i = 0 ; i < arrayMax (units) ; i += 3)
			{
			  pos1 = pos2 = arr (units, i + 1, BSunit).i ;
			  seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
			  seg1->x1 = pos1 - 1 ; seg1->x2 = pos1 - 1 ;
			  seg1->key = str2tag ("Valid5p") ;
			  seg1->data.i = arr (units, i + 2, BSunit).i ;
			  if (bsFindTag (Mrna, str2tag ("Transpliced_to")))
			    seg1->data.i |= (2 << 24) ;
			  if (bsFindTag (Mrna, str2tag ("Aggregated_5p_clones")))
			    seg1->data.i |= (1 << 24) ;
			  seg1->parent = seg->key ; 
			  seg1->source = 10*frame + delta ; /* needed for correct order */
			  seg1->sourceDx = b2 - b1 + 1 ;
			  seg1->type = (isDown) ?  MPRODUCT : MPRODUCT_UP ;
			}
		    }
		  keySetDestroy (testKs) ;
		}
	      bsDestroy (Mrna) ;
	    }
	}

      /***************/

         /* spliced cDNA ,  mieg */
      pos1 = 0 ;
      if  ( 
	   fmapView (view, (fMap_rView != isDown) ? _Fmap_Spliced_cDNA_Down : _Fmap_Spliced_cDNA_Up) &&
	   (
	    (   /* case 1: (transcribed genes) normal objects with a single coord system  */
	     bsFindTag (obj, _Assembled_from)
	     ) ||
	    (   /* case 2: (mrnas) Janus objects with a double coordinate system */ 
	     iseg < 2 && bsFindTag (obj, str2tag ("Constructed_from"))
	     )
	    ) && 
	   (/* look->max - look->min < 200000) && */
	   bsFlatten (obj, 5, units)
	   ))
	{ 
	  BOOL foundGeneName = FALSE ;
	  KEY fullRead = 0, oldest = 0 ;
	  int i0 ;
	  KEYSET tilingEst = queryKey (seg->key, " {>Read IS NM*} $| {>mrna ;>Mrna_covered_by} $| {>mrna ;>CDS_covered_by} ;>cdna_clone ; >read ") ;

	  if ((fullRead = keyGetKey (seg->key, str2tag ("Best_clone"))) &&
	      keyFindTag (fullRead, _Full_name)) ;
	  else
	    fullRead = 0 ;

	  seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	  seg1->key = M_SPLICED_cDNA ;
	  seg1->type = (isDown) ? SPLICED_cDNA : SPLICED_cDNA_UP ;
	  seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	  seg1->key = M_EXON_SUPPORT ;
	  seg1->type = (isDown) ? EXON_SUPPORT : EXON_SUPPORT_UP ;
	  for (i = 0 ; !foundGeneName && i < arrayMax (units) ; i += 5)
	    { 
	      KEY est = arr (units, i+2, BSunit).k ;

	      if (fullRead && est != fullRead)
		continue ;
	      if (look->view && fmapView (look->view, str2tag ("Fmap_CDS_Tiling")) &&
		  ! keySetFind (tilingEst, est, 0))
		continue ;
	      if (est != oldest && keyFindTag (est, _Full_name))
		{
		  oldest = est ;
		  if (foundGeneName)
		    {
		      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
		      
		      seg1->key = M_GENE_NAME ;
		      seg1->type = (isDown) ? CDNA_GENE_NAME : CDNA_GENE_NAME_UP ;
		      seg1->x1 = seg->x1 ; seg1->x2 = seg->x2 ;
		    }
		}
	    }
	  i0 = nsegs ;
	  for (i = 0 ; i < arrayMax (units) ; i += 5)
	    { 
	      KEY est ;

	      if (arrayMax (units) < i + 5)
		continue ;
	  
	      est = arr (units, i+2, BSunit).k ;
	      if (look->view && fmapView (look->view, str2tag ("Fmap_CDS_Tiling")) &&
		  ! keySetFind (tilingEst, est, 0))
		continue ;

	      pos1 = arr (units, i + 0, BSunit).i ;
	      pos2 = arr (units, i + 1, BSunit).i ;
	      
	      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	      POS_TO_SEG1 ;
	      seg1->key = est ;
	      seg1->source = seg->source /* to order est gene per gene */ ;
	      seg1->sourceDx = seg->x2 - seg->x1 - 1 ;
	      pos1 = arr (units, i + 3, BSunit).i ;
	      pos2 = arr (units, i + 4, BSunit).i ;
	      /* reserve bit 31 for ambiguous, 30 polyA, unaligned end, 28 unaligned begin */
	      seg1->data.i = ((pos1 & 0x3FFF) << 14)  +  (pos2 & 0x3FFF) ;
	      
	      fMapcDNADoFillData (seg1) ; /* implies objCreate, unfortunate, but needed for ordering */ 
	      /* now seg1->parent == KEYKEY (cdna_clone */
	      seg1->type = (isDown) ? SPLICED_cDNA : SPLICED_cDNA_UP ;
	      if (! strncmp (name(seg1->key), "XH_",3)|| ! strncmp (name(seg1->key), "XG_",3))
		seg1->type = (isDown) ? EXON_SUPPORT : EXON_SUPPORT_UP ;
	      seg1->parent = seg->x1 + 1 ;  /* avoid zero , as happens in mRNA case where parent is main object*/
	    }
	  for (i = i0 ; i < nsegs ; i++) /* get the lowest x1 for each clone */
	    {
	      int j, uumin, nuumin = 0 ;
	      KEY myparent, mysource, mytag ;
	      BOOL myup ;
	      seg = arrp (segs, i , SEG) ;
	      if ((seg->type | 0x1) != SPLICED_cDNA_UP &&
		  (seg->type | 0x1) != EXON_SUPPORT_UP)
		continue ;

	      if (class (seg->key) == _VMethod)
		continue ;
	      mytag = seg->type ;
	      myup = mytag & 0x1 ? TRUE : FALSE ;
	      uumin = myup ? seg->x2 : seg->x1 ;
	      nuumin = i - iseg ;
	      myparent = seg->parent ;
	      mysource = seg->source ;
	      if (!myparent)
		{ seg->parent = (1 << 30) | seg->x1 ; continue ; }
	      if (class (myparent))  /* already treated in the recursion */
		continue ; 
	      for (j = i  ; j < nsegs ; j++)
		{
		  seg1 = arrp (segs, j , SEG) ; 
		  if (seg1->type != mytag)
		    continue ;
		  if (seg1->parent != myparent)
		    continue ;
		  if (seg1->source != mysource)
		    continue ;
		  if (!myup && uumin > seg1->x1)
		    uumin = seg1->x1 ;
		  else if (myup && uumin < seg1->x2)
		    uumin = seg1->x2 ;
		}
	      uumin = 100 * uumin  + nuumin ;
	      if (uumin > (1 << 29)) 
		uumin =  (1 << 29) ;
	      if (myup)
		uumin =  - uumin ;
	      uumin |=  (1 << 30) ; /* sets non zero class to stop recursion */
	      for (j = i ; j < nsegs ; j++)
		{
		  seg1 = arrp (segs, j , SEG) ; 
		  if (seg1->type != mytag)
		    continue ;
		  if (seg1->parent != myparent)
		    continue ;
		  if (seg1->source != mysource)
		    continue ;
		  seg1->parent =  uumin ; /* keep clones unequal */
		}
	    }
	  for (i = i0 ; i < nsegs ; i++) /* get the lowest x1 for each clone */
	    {
	      seg1 = arrp (segs, i , SEG) ;
	      switch (seg1->type)
		{
		  /* cheat so that tilingEst appear globaly left of the non tiling est if fMapOrder */
		case SPLICED_cDNA:
		case EXON_SUPPORT:
		  if (keySetFind (tilingEst, seg1->key, 0))
		    seg1->source-- ; 
		  break ;
		case SPLICED_cDNA_UP:
		case EXON_SUPPORT_UP:
		  if (keySetFind (tilingEst, seg1->key, 0))
		    seg1->source++ ; 
		  break ;
		default:
		  break ;
		}
	    }
	  seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;

	  seg1->key = M_SPLICED_cDNA_DECORATE ;
	  seg1->type = (isDown) ? SPLICED_cDNA_DECORATE : SPLICED_cDNA_DECORATE_UP ;
	  seg1->x1 = seg->x1 ; seg1->x2 = seg->x2 ;
	  keySetDestroy (tilingEst) ;
	}

  	 /* virtual sequence assembly tags, */ /* mieg */
      if ((seg->type | 0x1) == SEQUENCE_UP &&
	  bsFindTag (obj, _Assembled_from) && 
	  bsFlatten (obj, 3, units))
	for (i = 0 ; i < arrayMax (units)/3 ; ++i)
	  { int seg1index = nsegs ;
	    seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	    *seg1 = *seg ; seg1->source = parent ;
	    /* get start and end of tag */
	    pos1 = arr (units, 3*i+1, BSunit).i;
	    pos2 = arr (units, 3*i+2, BSunit).i;
	   
	    /* Convert co-ords of pos1, pos2 and stuff into seg as x1 and x2 */
	    posConvert (seg1, pos1, pos2);
	    if (seg1->type == SEQUENCE)
	      {
		if (seg1->x1 > seg1->x2)
		  seg1->type = VIRTUAL_SUB_SEQUENCE_TAG_UP ;
		else
		  seg1->type = VIRTUAL_SUB_SEQUENCE_TAG ;
	      }
	    else /* sequence_up */
	      {
		if (seg1->x1 < seg1->x2)
		  seg1->type = VIRTUAL_SUB_SEQUENCE_TAG_UP ;
		else
		  seg1->type = VIRTUAL_SUB_SEQUENCE_TAG ;
	      }
	    if (seg1->x1 > seg1->x2)
	      { int tmp = seg1->x1 ; seg1->x1 = seg1->x2 ; seg1->x2 = tmp ; }
	    /* Key is the subsequence */
	    seg1->key = arr (units, 3*i, BSunit).k ;
	    seg1->parent = 0 ;	
	    seg1->data.i = 0 ;

	    { OBJ obj1 ;
	      KEY _Proposed, _Stolen, _New_Read ;

	      lexaddkey ("Proposed", &_Proposed, _VSystem) ;
	      lexaddkey ("Stolen", &_Stolen, _VSystem) ;
	      lexaddkey ("New_Read", &_New_Read, _VSystem) ;
	      if ((obj1 = bsCreate (seg1->key)))
		{ KEY key1, myCol ;
		  if (bsFindTag (obj1, _Stolen))
		    seg1->data.i |= 1 << 24 ;
		  if (bsFindTag (obj1, _Proposed))
		    seg1->data.i |= 1 << 25 ;
		  if (bsFindTag (obj1, _Vector))
		    seg1->data.i |= 1 << 26 ;
		  if (bsFindTag (obj1, _New_Read))
		    seg1->data.i |= 1 << 27 ;
		  if (bsFindTag (obj1, _Significant_bases))
		    seg1->data.i |= ((unsigned int)1) << 31 ;
		  if (bsGetKeyTags (obj1, _Colour, &myCol))
		    seg1->data.i |= ((myCol - _WHITE + WHITE) & 0x3f) << 16 ;
#ifdef ACEMBLY
		  fMapColorSeg (seg1) ;
#endif

		  /* primers (directed sequencing */ /* mieg */
		  if (bsGetKey (obj1, _Primer, &key1))
		      do 
		      { seg2 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
			*seg2 = arr (segs,seg1index,SEG) ;
			seg2->key = key1 ;
			seg2->parent = key1 ; /* so both half co color */
			seg2->data.i = 0 ;
			seg2->type = PRIMER | (seg1->type & 1) ;
		      }  while (bsGetKey (obj, _bsDown, &key1)) ;
		  bsDestroy (obj1) ;
		}
	    }
	  }
  	 /* virtual sequence assembly tags, dummy seq */ /* ulrich */
      if (bsFindTag (obj, _Previous_contig) && 
	  bsFlatten (obj, 3, units))
	for (i = 0 ; i < arrayMax (units)/3 ; ++i)
	  { seg1 = arrayp (segs, nsegs++, SEG) ; seg = arrp (segs, iseg, SEG) ;
	    *seg1 = *seg ; seg1->source = parent ;
	    /* get start and end of tag */
	    pos1 = arr (units, 3*i+1, BSunit).i;
	    pos2 = arr (units, 3*i+2, BSunit).i;
	   
	    /* Convert co-ords of pos1, pos2 and stuff into seg as x1 and x2 */
	    posConvert (seg1, pos1, pos2);
	    if (seg1->type == SEQUENCE)
	      {
		if (seg1->x1 > seg1->x2)
		  seg1->type = VIRTUAL_PREVIOUS_CONTIG_TAG_UP ;
		else
		  seg1->type = VIRTUAL_PREVIOUS_CONTIG_TAG ;
	      }
	    else /* sequence_up */
	      {
		if (seg1->x1 < seg1->x2)
		  seg1->type = VIRTUAL_PREVIOUS_CONTIG_TAG_UP ;
		else
		  seg1->type = VIRTUAL_PREVIOUS_CONTIG_TAG ;
	      }
	    if (seg1->x1 > seg1->x2)
	      { int tmp = seg1->x1 ; seg1->x1 = seg1->x2 ; seg1->x2 = tmp ; }

	    /* Key is the subsequence */
	    seg1->key = arr (units, 3*i, BSunit).k ;
	    seg1->parent = 0 ;	
	    seg1->data.i = 0 ;
	  }

  	 /* virtual alignment */ /* mieg */
      if (bsFindTag (obj, _Aligned) && 
	  bsFlatten (obj, 3, units))
	for (i = 0 ; i < arrayMax (units)/3 ; ++i)
	  { seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	    *seg1 = *seg ; seg1->source = parent ;
	    /* get start and end of tag */
	    pos1 = arr (units, 3*i+1, BSunit).i;
	    pos2 = arr (units, 3*i+2, BSunit).i;
	   
	    /* Convert co-ords of pos1, pos2 and stuff into seg as x1 and x2 */
	    posConvert (seg1, pos1, pos2);
	    if (seg1->type == SEQUENCE)
	      {
		if (seg1->x1 > seg1->x2)
		  seg1->type = VIRTUAL_ALIGNED_TAG_UP ;
		else
		  seg1->type = VIRTUAL_ALIGNED_TAG ;
	      }
	    else /* sequence_up */
	      {
		if (seg1->x1 < seg1->x2)
		  seg1->type = VIRTUAL_ALIGNED_TAG_UP ;
		else
		  seg1->type = VIRTUAL_ALIGNED_TAG ;
	      }
	    if (seg1->x1 > seg1->x2)
	      { int tmp = seg1->x1 ; seg1->x1 = seg1->x2 ; seg1->x2 = tmp ; }
	    /* Key is the subsequence */
	    seg1->key = arr (units, 3*i, BSunit).k ;
	    seg1->parent = 0 ;	
	    seg1->data.i = 0 ;
	  }


        	 /* virtual sequence assembly tags, */ /* mieg */
      if (bsGetKey (obj, _Assembled_into, &key)) do
	{ seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	  *seg1 = *seg ; seg1->source = parent ;
	  seg1->key = key ;
	  seg1->type = VIRTUAL_PARENT_SEQUENCE_TAG ;
	  seg1->parent = seg->key ;	
	  seg1->data.i = 0 ;
	  { OBJ obj1 ;
	    KEY Proposed ;

	    lexaddkey ("Proposed", &Proposed, _VSystem) ;
	    if ((obj1 = bsCreate (seg1->key)))
	      { if (bsFindTag (obj1, Proposed))
		  seg1->data.i = -1 ;
		bsDestroy (obj1) ;
	      }
	    }	
	  seg1->x1 = seg->x1 ;	
	  seg1->x2 = seg->x2 ;	
	}  while (bsGetKey (obj, _bsDown, &key)) ;
      
#define POS_PROCESS(Z,Z_UP) \
    if (isDown)  seg1->type = (pos2 > pos1) ? Z : Z_UP ; \
    else         seg1->type = (pos1 > pos2) ? Z : Z_UP ; \
    if (pos1 > pos2) { tmp = pos2 ; pos2 = pos1 ; pos1 = tmp ;} \
    posConvert (seg1, pos1, pos2) ;

         /* Probes, mieg */ 
      if (fmapView (view, _Fmap_Probe) && bsFindTag (obj, _Probe_hit) && bsFlatten (obj, 3, units))
	{
	  for (i = 0 ; i < arrayMax (units) ; i += 3)
	    {
	      pos1 = arr (units,i+1, BSunit).i ; 
	      pos2 = arr (units,i+2, BSunit).i ;
	      if (class (bsKey(obj)) == _VmRNA &&
		  !fMapcDNAProbePosition (look, obj, seg, arr (units,i, BSunit).k, &pos1, &pos2))
		continue ; /* may be the probe is already displayed */
	      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	      seg1->source = parent ;
	      seg1->key = arr (units,i, BSunit).k ;
	      seg1->parent = seg1->key  ; 
	      seg1->data.i = 0 ;
	      POS_PROCESS (PROBE,PROBE_UP) ;
	      if (bsFindKey (obj, _Probe_exact_hit, seg1->key))
		seg1->data.i |=  0x4000 ; /* exact */
	    }
	}
      
         /* Oligos, mieg */ 
      if (fmapView (view, _Fmap_Oligo) && bsFindTag (obj, _Oligo) && bsFlatten (obj, 3, units))
	for (i = 0 ; i < arrayMax (units) ; i += 3)
	  { OBJ obj1 ;
	    seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	    seg1->source = parent ;
	    seg1->key = arr (units,i, BSunit).k ;
	    seg1->parent = seg1->key  ; 
	    seg1->data.i = 0 ;
	    pos1 = arr (units,i+1, BSunit).i ; 
	    pos2 = arr (units,i+2, BSunit).i ;
	    if (!pos1 || !pos2) 
	      fMapOspPositionOligo (look, seg, seg1->key, &pos1, &pos2) ;
	    if (!pos1 || !pos2) 
	      continue ;
	    POS_PROCESS (OLIGO,OLIGO_UP) ;
	    if ((obj1 = bsCreate (seg1->key)))
	      { float xf ; int xi ;
	      if (bsGetData (obj1, _Tm, _Float, &xf))
		{ xi = 10 * xf ; seg1->data.i = (xi & 0x3FF) << 16 ; }
	      if (bsGetData (obj1, _Score, _Float, &xf))
		{ xi = xf ; seg1->data.i |= (xi & 0xfff) ; }     
	      if (!bsFindTag (obj1, _Temporary)) 
		seg1->data.i |=  0x4000 ; /* old */
	      bsDestroy (obj1) ;
	      }
	  }
      
         /* Alleles */
      if (fmapView (view, _Fmap_Allele) && bsFindTag (obj, _Allele) && bsFlatten (obj, 4, units))
	for (i = 0 ; i < arrayMax (units) ; i += 4)
	  { char *string ;
	    seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	    seg1->source = parent ;
	    seg1->key = arr (units,i, BSunit).k ;
	    seg1->parent = 0 ;
	    seg1->data.i = 0 ;
	    pos1 = arr (units,i+1, BSunit).i ; 
	    pos2 = arr (units,i+2, BSunit).i ;
	    if (pos1 || pos2) 
	      seg1->parent = seg1->key ; /* flag to draw it */
	    else
	      seg1->parent = seg->key ; /* couple with parent */
	    POS_PROCESS (ALLELE,ALLELE_UP) ;
	    if ((string = arr (units,i+3,BSunit).s))
	      seg1->data.s = strnew (string, look->segsHandle) ;
	    else
	      seg1->data.s = 0 ;
	  }
    
      /* EMBL features */
      if (bsFindTag (obj, _EMBL_feature) && bsFlatten (obj, 4, units))
	for (i = 0 ; i < arrayMax (units) ; i += 4)
	  { char *string ;
	    seg1 =  arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	    seg1->source = parent ;
	    seg1->parent = nsegs ;
	    seg1->key = arr (units,i, BSunit).k ;
	    pos1 = arr (units,i+1, BSunit).i ; 
	    pos2 = arr (units,i+2, BSunit).i ; 
	    POS_PROCESS (EMBL_FEATURE, EMBL_FEATURE_UP) ;
			/* 4'th column: can't distinguish KEY from Text securely 
			   assume if iskey () that it is a key,  e.g. ?Text
			*/
	    if (iskey (key = arr (units,i+3,BSunit).k))
	      seg1->data.s = strnew (name (key), look->segsHandle) ;
	    else if ((string = arr (units,i+3,BSunit).s))
	      seg1->data.s = strnew (string, look->segsHandle) ;
	    else
	      seg1->data.s = 0 ;
	  }

	 /* Homologies */
      { KEY _Bury = str2tag ("Bury") ;
	HOMOLINFO *hinf ;
	
	if ((1 || fmapView (view, _Fmap_Homol)) && bsFindTag (obj, _Homol) && bsFlatten (obj, 9, units))
	  for (i = 0 ; i < arrayMax (units) ; i += 9)
	    { 
	      u = arrayp (units, i, BSunit) ; /* skip tag2 */
	      if (! fmapView (view, _Fmap_AceKogN) && u[0].k == _AceKogN)
		continue ;
	      u = arrayp (units, i+1, BSunit) ; /* skip tag2 */
	      if (class (u[1].k) != _VMethod || !u[3].i || !u[4].i)
		continue ;
				/* check for #HomolInfo data */
	      if (i && (u[6].i == u[-3].i) && (u[5].i == u[-4].i) &&
		  (u[4].i == u[-5].i) && (u[3].i == u[-6].i) &&
		  (u[2].f == u[-7].f) && (u[1].k == u[-8].k) &&
		  (u[0].k == u[-9].k) && (u[-1].k == u[-10].k))
		{		/* seg1 is correct - check flags */
		  if (u[7].k == _Bury)
		    bitSet (look->homolBury, seg1->data.i) ;
		  continue ;
		}
	      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	      seg1->source = parent ;
	      seg1->key = seg1->parent = u[0].k ; /* target */
	      seg1->data.i = arrayMax (look->homolInfo) ;
	      hinf = arrayp (look->homolInfo, seg1->data.i, HOMOLINFO) ;
	      hinf->method = u[1].k ;
	      hinf->score = u[2].f ;
	      if (u[5].i > u[6].i)
		{ pos1 = u[4].i ; pos2 = u[3].i ;
		  hinf->x1 = u[6].i ; hinf->x2 = u[5].i ;
		}
	      else
		{ pos1 = u[3].i ; pos2 = u[4].i ;
		  hinf->x1 = u[5].i ; hinf->x2 = u[6].i ;
		}
	      POS_PROCESS (HOMOL,HOMOL_UP) ;
	      if (u[7].k == _Bury)
		bitSet (look->homolBury, seg1->data.i) ;
	    }
      }

      /* necessary to ensure array long enough */
      bitExtend (look->homolBury, arrayMax (look->homolInfo)) ;
      bitExtend (look->homolFromTable, arrayMax (look->homolInfo)) ;
      
 
         /* Features */
      if (bsFindTag (obj, _Feature) && bsFlatten (obj, 5, units))
	for (i = 0 ; i < arrayMax (units) ; i += 5)
	  { u = arrayp (units, i, BSunit) ;
	    if (class (u[0].k) != _VMethod || !u[1].i || !u[2].i)
	      continue ;
	    seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	    seg1->source = parent ;
	    seg1->key = u[0].k ;
	    pos1 = u[1].i ; 
	    pos2 = u[2].i ; 
	    POS_PROCESS (FEATURE,FEATURE_UP) ;
	    seg1->data.f = u[3].f ;
	    if (u[4].s)
	      { dictAdd (look->featDict, u[4].s, &tmp) ;
		seg1->parent = -tmp ;
	      }
	    else
	      seg1->parent = 0 ;
	  }

         /* Splice sites */
      { KEY _Predicted_5 = str2tag ("Predicted_5") ;
	KEY _Predicted_3 = str2tag ("Predicted_3") ;

	if (bsFindTag (obj, str2tag ("Splices")) && bsFlatten (obj, 5, units))
	  for (i = 0 ; i < arrayMax (units) ; i += 5)
	    { u = arrayp (units, i, BSunit) ;
	      if (u[0].k != _Predicted_5 && u[0].k != _Predicted_3)
		continue ;
	      if (class (u[1].k) != _VMethod || !u[2].i || !u[3].i)
		continue ;
	      seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
	      seg1->source = parent ;
	      seg1->key = u[1].k ;
	      pos1 = u[2].i ; 
	      pos2 = u[3].i ;
	      if (u[0].k == _Predicted_5)
		{ POS_PROCESS (SPLICE5,SPLICE5_UP) ; }
	      else
		{ POS_PROCESS (SPLICE3,SPLICE3_UP) ; }
	      seg1->data.f = u[4].f ;
	    }
      }

	/* clone end information */
     { int z ;
	KEY tag ;
	for (z = 0 ; z < 2 ; ++z)
	  { tag = z ? _Clone_left_end : _Clone_right_end ;
	    if (bsFindTag (obj, tag) && bsFlatten (obj,2,units))
	      for (i = 0 ; i < arrayMax (units) ; i += 2)
		{ u = arrayp (units, i, BSunit) ;
		  seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
		  seg1->source = parent ;
		  seg1->key = u[0].k ; 		/* Clone */
		  seg1->type = CLONE_END ;
		  pos1 = u[1].i ; pos2 = u[1].i ; posConvert (seg1, pos1, pos2) ;
		  seg1->data.k = tag ;
		  lexaddkey (name (seg1->key), &key, _VSequence) ;
		  seg1->parent = key ;
		}
	  }
      }

  /* confirmed introns: add all now, and remove duplicates with those
     from genes later in fMapFindCoding ()
  */
    
      { KEY _Confirmed_intron = str2tag ("Confirmed_intron") ;
	KEY _EST = str2tag ("EST") ;
	KEY _cDNA = str2tag ("cDNA") ;
	KEY _Homology = str2tag ("Homology") ;
	KEY _UTR = str2tag ("UTR") ;

	if (bsFindTag (obj, _Confirmed_intron) &&
	    bsFlatten (obj, 3, units))
	  for (i = 0 ; i < arrayMax (units) ; i += 3)
	    { u = arrayp (units, i, BSunit) ;
	      if (!u[0].i || !u[1].i)
		continue ;
	      if (!i || u[-3].i != u[0].i || u[-2].i != u[1].i)
		{ seg1 = arrayp (segs,nsegs++,SEG) ; seg = arrp (segs,iseg,SEG) ;
		  seg1->source = parent ;
		  seg1->key = 0 ;
		  pos1 = u[0].i ; 
		  pos2 = u[1].i ; 
		  seg1->data.i = 0 ;
		  POS_PROCESS (INTRON,INTRON_UP) ;
		}
	      seg1->data.i = arrayMax (look->seqInfo) ;
	      sinf = arrayp (look->seqInfo, seg1->data.i, SEQINFO) ;
	      if (u[2].k == _EST)
		sinf->flags |= SEQ_CONFIRM_EST ; 
	      else if (u[2].k == _cDNA)
		sinf->flags |= SEQ_CONFIRM_CDNA ; 
	      else if (u[2].k == _Homology)
		sinf->flags |= SEQ_CONFIRM_HOMOL ; 
	      else if (u[2].k == _UTR)
		sinf->flags |= SEQ_CONFIRM_UTR ; 
	      else
		sinf->flags |= SEQ_CONFIRM_UNKNOWN ; 
	    }
      }

      bsDestroy (obj) ;
    }

/*-------- end of extraction of SEGs ----------*/

  if (look->flag & FLAG_COMPLEMENT)
    { int top = look->length - 1 ;
      seg = arrp (look->segs, 1, SEG) ;
      for (i = 1 ; i < arrayMax (segs) ; ++i, ++seg)
	{ tmp = seg->x1 ;
	  seg->x1 = top - seg->x2 ;
	  seg->x2 = top - tmp ;
	  if (seg->type >= SEQUENCE && seg->type <= ALLELE_UP)
	    seg->type ^= 1 ;
	  if ((seg->type | 0x1) == TRANSCRIBEDGENE_UP)
	    seg->source = top - seg->source - seg->sourceDx ;
	  if ((seg->type | 0x1) == EXON_SUPPORT_UP)
	    seg->source = top - seg->source  - seg->sourceDx;
	  if ((seg->type | 0x1) == SPLICED_cDNA_UP)
	    seg->source = top - seg->source  - seg->sourceDx;
	  if ((seg->type | 0x1) == TRANSCRIPT_UP)
	    seg->source = top - seg->source  - seg->sourceDx;
	  if ((seg->type | 0x1) == MRNA_UP)
	    seg->source = top - seg->source  - seg->sourceDx;
	  if ((seg->type | 0x1) == PMRNA_UP)
	    seg->source = top - seg->source  - seg->sourceDx;
	  if ((seg->type | 0x1) == MPRODUCT_UP)
	    seg->source = top - seg->source  - seg->sourceDx;
	}
    }

  if (oldSegs)
    { arrayMax (oldSegs) = look->lastTrueSeg ; /* still the old value */
      addOldSegs (look, segs, oldSegs) ;
      handleDestroy (oldSegsHandle) ; /* kills old segs and related stuff */
    }

  arraySort (segs, fMapOrder) ;	/* must sort for FindCoding and RemoveSelfHomol */
  fMapFindCoding (look) ;	/* make CODING segs */
  fMapRemoveSelfHomol (look) ;
  fMapTraceFindMultiplets (look) ;  /* make virtual multiplets */

  fMapProcessMethods (look) ;	/* always do this - to pick up edits */

  arraySort (segs, fMapOrder) ;
  /* these 2 routines no longer work after reverse complement because RC calls sort
   * they should be written differently
   */

  fMapOspFindOligoPairs (look) ;    /* why  AFTER sorting ? not sure, but probably ok */
  look->lastTrueSeg = arrayMax (segs) ;

  chronoReturn () ;
  arrayDestroy (units) ;
  return TRUE ;
}

/**************************/

BOOL fMapFindSpan (LOOK look, KEY *key, int *x, int *y)
{
  SEG *seg, *minSeg = 0 ;
  int i, min = look->fullLength ;
  OBJ obj ;

  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs, i, SEG) ;
      if ((seg->type == SEQUENCE || seg->type == SEQUENCE_UP) &&
	  seg->key != *key &&
	  seg->x1 < *x && seg->x1 < *y &&
	  seg->x2+1 >= *x && seg->x2+1 >= *y &&
	  seg->x2 - seg->x1 + 1 <= min &&
	  (obj = bsCreate (seg->key)))
	{ if (!bsFindTag (obj, _Source_Exons)) /* important! */
	    { min = seg->x2 - seg->x1 + 1 ;
	      minSeg = seg ;
	    }
	  bsDestroy (obj) ;
	}
    }

  if (!minSeg)
    return FALSE ;
  if (minSeg->type == SEQUENCE)
    { *x -= minSeg->x1 ;
      *y -= minSeg->x1 ;
    }
  else
    { *x = minSeg->x2 + 2 - *x ;
      *y = minSeg->x2 + 2 - *y ;
    }
  *key = minSeg->key ;
  return TRUE ;
}

/*****************************************************************/
/*****************************************************************/

char* fMapSegTypeName[] = { 
  "MASTER", "ASSEMBLY_TAG",

  "SEQUENCE", "SEQUENCE_UP",
  "CDS", "CDS_UP",
  "INTRON", "INTRON_UP",
  "EXON", "EXON_UP",
  "EXON_GAP", "EXON_GAP_UP",
  "EMBL_FEATURE", "EMBL_FEATURE_UP", 
  "FEATURE", "FEATURE_UP", 
  "ATG", "ATG_UP",
  "SPLICE3", "SPLICE3_UP",
  "SPLICE5", "SPLICE5_UP",
  "CODING", "CODING_UP",
  "TRANSCRIBEDGENE", "TRANSCRIBEDGENE_UP",
  "EXON_SUPPORT", "EXON_SUPPORT_UP",
  "SPLICED_cDNA", "SPLICED_cDNA_UP",
  "SPLICED_cDNA_DECORATE", "SPLICED_cDNA_DECORATE_UP",
  "GENE_NAME", "GENE_NAME_UP",
  "VIRTUAL_SUB_SEQUENCE_TAG",  "VIRTUAL_SUB_SEQUENCE_TAG_UP",
  "VIRTUAL_MULTIPLET_TAG", "VIRTUAL_MULTIPLET_TAG_UP",
  "VIRTUAL_ALIGNED_TAG",  "VIRTUAL_ALIGNED_TAG_UP",
  "VIRTUAL_PREVIOUS_CONTIG_TAG",  "VIRTUAL_PREVIOUS_CONTIG_TAG_UP",
  "HOMOL", "HOMOL_UP",
  "PRIMER", "PRIMER_UP",
  "OLIGO", "OLIGO_UP",
  "OLIGO_PAIR", "OLIGO_PAIR_UP",
  "TRANS_SEQ", "TRANS_SEQ_UP",
  "ALLELE", "ALLELE_UP",
  "VISIBLE",
  "DNA_SEQ", "PEP_SEQ", "ORF", 
  "VIRTUAL_PARENT_SEQUENCE_TAG",
  "VIRTUAL_CONTIG_TAG",
  "CLONE_END"
  } ;

static char *title (KEY key)
{
  OBJ obj ;
  char *result ;
  KEY textKey ;

  result = name (key) ;

  if ((obj = bsCreate (key)))
    { if (bsGetKey (obj, _Title, &textKey))
	result = name (textKey) ;
      bsDestroy (obj) ;
    }

  return result ;
}

/**********************************************************/

void fMapReportLine (LOOK look, SEG *seg, BOOL isGIF, float x)
{ 
  char *cp ;
  HOMOLINFO *hinf ;
  SEQINFO *sinf ;
  METHOD *meth ;
  int i, i2, c1, ii ;
  SEG *seg2 ;  

  memset (look->segTextBuf, 0, 511) ;  /* mieg */

				/* object name */
  if (isGIF)
    *look->segTextBuf = 0 ;
  else if (iskey (seg->key))
    strncpy (look->segTextBuf, name (seg->key), 40) ; 
  else if (seg->parent && iskey (seg->parent))
    strncpy (look->segTextBuf, name (seg->parent), 40) ;
  else
    *look->segTextBuf = 0 ;
  
				/* coordinates */
  if (look->flag & FLAG_REVERSE)
    strcat (look->segTextBuf, 
	    messprintf ("    %d %d (%d)  ", 
			COORD (look, seg->x2), COORD (look, seg->x1),
			1 - COORD (look, seg->x2) + COORD (look, seg->x1))) ;
  else
    strcat (look->segTextBuf, 
	    messprintf ("    %d %d (%d)  ", 
			COORD (look, seg->x1), COORD (look, seg->x2),
			1 - COORD (look, seg->x1) + COORD (look, seg->x2))) ;
  
  if (0 && seg->key && keyFindTag (seg->key, _Composite))
	{
	  OBJ Obj = bsCreate (seg->key) ;
	  int n = 0 ;
	  if (bsGetData (Obj, _Composite, _Int, &n) && n > 0)
	    strcat (look->segTextBuf, 
		    messprintf ("Seen %d", n)) ;
	  bsDestroy (Obj) ;
	}


				/* description */
  if (seqDragBox)
    seg2 = BOXSEG (seqDragBox) ;
  else
    seg2 = seg ;
  switch (seg->type)
    {
    case SEQUENCE: case SEQUENCE_UP:
    case EXON: case EXON_UP:
    case EXON_GAP: case EXON_GAP_UP:
    case INTRON: case INTRON_UP:
      sinf = arrayp (look->seqInfo, seg->data.i, SEQINFO) ;
      if (sinf->method)
	strncat (look->segTextBuf, name (sinf->method), 40) ;
      if (sinf->flags & SEQ_SCORE)
	strncat (look->segTextBuf, messprintf ("%.1f", sinf->score), 10) ;
      if (sinf->flags & SEQ_CONFIRMED)
	{ strcat (look->segTextBuf, "Confirmed") ;
	  if (sinf->flags & SEQ_CONFIRM_CDNA)
	    strcat (look->segTextBuf, " by_cDNA") ;
	  if (sinf->flags & SEQ_CONFIRM_EST)
	    strcat (look->segTextBuf, " by_EST") ;
	  if (sinf->flags & SEQ_CONFIRM_HOMOL)
	    strcat (look->segTextBuf, " by_homology") ;
	  if (sinf->flags & SEQ_CONFIRM_UTR)
	    strcat (look->segTextBuf, " in_UTR") ;
	}
      break ;
    case DNA_SEQ:
      i = seg->x1 + (int)x ;
      i2 = seg2->x1 + (int)newx ;
      if (look->flag & FLAG_REVERSE)
	strcat (look->segTextBuf, 
		messprintf ("Selection %d ---> %d (%d)", 
			    COORD (look, i2), COORD (look, i),
			  1 - COORD (look, i2)+ COORD (look, i))) ;
      else
	strcat (look->segTextBuf, 
		messprintf ("Selection %d ---> %d (%d)", 
			    COORD (look, i), COORD (look, i2),
			    1 + COORD (look, i2)- COORD (look, i))) ;
      if (look->gf.cum && i >= look->gf.min && i < look->gf.max)
	strncat (look->segTextBuf,
		 messprintf ("  cumulative coding score %.2f",
			     look->gf.cum[i - look->gf.min]), 64) ;
      break ;
    case PEP_SEQ:
      i = seg->x1 + 3* (int)x ;
      i2 = seg2->x1 + 3* (int)newx ;
      if ((cp = pepName[ (int) codon (arrp (look->dna, i, char))]))
	{
	  if (look->flag & FLAG_REVERSE)
	    strncat (look->segTextBuf, 
		     messprintf ("%s at %d to %d Selection weight: %d",cp,
				 COORD (look, i+2), COORD (look,i), weight), 64) ;
	  else
	    strncat (look->segTextBuf, 
		     messprintf ("%s at %d to %d Selection weight: %d",cp,
				 COORD (look, i), COORD (look,i+2), weight), 64) ;
	}
      break ;
    case TRANS_SEQ:
      i = seg->x1 + 3* (int)x ;
      c1 =  COORD (look, i) ;
      if ((cp = pepName[ (int) codon (arrp (look->dna, i, char))]))
	{
	  if (look->flag & FLAG_REVERSE)
	    { 
	      if (!seqDragBox) nbTrans = 1 ;
	      strncat (look->segTextBuf,
		       messprintf ("%s at %d ", cp, c1 -  nbTrans + 1), 64) ;
	      if (seqDragBox && c1 <= array (look->minDna, 0, int))
		{
		  for (ii=1; ii<arrayMax (look->maxDna); ii++)
		    if (c1 > array (look->maxDna, ii, int) &&
			c1 <= array (look->minDna, ii- 1, int))
		      { 
			c1 = (array (look->decalCoord, ii-1, int) - c1)/3 ;
			break ;
		      }
		  strcat (look->segTextBuf, 
			  messprintf (" Selection %d ---> %d (%d) weight: %d", 
				      c1, c1 + nbTrans - 1, nbTrans, weight)) ;
		}
	    }
	  else
	    {
	      strncat (look->segTextBuf, 
		       messprintf ("%s at %d ", cp, c1), 64) ;
	      if (seqDragBox && c1 >= array (look->minDna, 0, int))
		{
		  for (ii=1; ii<arrayMax (look->maxDna); ii++)
		    if (c1 < array (look->maxDna, ii, int) &&
			c1 >= array (look->minDna, ii- 1, int))
		      { 
			c1 = (c1 - array (look->decalCoord, ii-1, int))/3 ;
			break ;
		      }
		  strcat (look->segTextBuf, 
			  messprintf ( " Selection %d ---> %d (%d) weight: %d", 
				       c1, c1 + nbTrans - 1, nbTrans, weight)) ;
		}
	    }
	}
      break ;
    case TRANS_SEQ_UP:
      i = seg->x1 + 3* (int)x ;
      c1 =  COORD (look, i) ;
      if ((cp =  pepName[ (int) reverseCodon (arrp (look->dna, i-2, char))]))
	{
	  if (look->flag & FLAG_REVERSE)
	    {
	      if (!seqDragBox) nbTrans = 1 ;
	      strncat (look->segTextBuf, 
		       messprintf ("up-strand %s at %d ", cp, c1 + nbTrans - 1), 64) ;
	      if (seqDragBox && c1 <= array (look->minDna, 0, int))
		{
		  for (ii=1; ii<arrayMax (look->maxDna); ii++)
		    if (c1 > array (look->maxDna, ii, int) &&
			c1 <= array (look->minDna, ii- 1, int) )
		      { 
			c1 = (c1 - array (look->decalCoord, ii-1, int))/3 ;
			break ;
		      }
		  strcat (look->segTextBuf, 
			  messprintf (" Selection %d ---> %d (%d) weight: %d", 
				      c1, c1 - nbTrans + 1, nbTrans, weight)) ;
		}
	    }	  
	  else
	    {
	      strncat (look->segTextBuf, 
		       messprintf ("up-strand %s at %d ", cp, c1), 64) ;
	      if (seqDragBox && c1 >= array (look->minDna, 0, int))
		{
		  for (ii=1; ii<arrayMax (look->maxDna); ii++)
		    if (c1 < array (look->maxDna, ii, int) &&
			c1 >= array (look->minDna, ii- 1, int))
		      { 
			c1 = (array (look->decalCoord, ii-1, int) - c1)/3 ;
			break ;
		      }
		  strcat (look->segTextBuf, 
			  messprintf ( " Selection %d ---> %d (%d) weight: %d", 
				       c1, c1 - nbTrans + 1, nbTrans, weight)) ;
		}
	    }
	}
      break ;
    case ORF:
      strncat (look->segTextBuf, messprintf ("Frame %d", seg->data.i), 64) ;
      break ;
    case SPLICE3: case SPLICE5: case ATG:
      strncat (look->segTextBuf, fMapSegTypeName[seg->type], 32) ;
      strncat (look->segTextBuf, messprintf (" Score %f", seg->data.f), 32) ;
      break ;
    case ASSEMBLY_TAG:
      strncat (look->segTextBuf, messprintf ("%s", seg->data.s), 64);
      break;
    case ALLELE: case ALLELE_UP: 
    case EMBL_FEATURE: case EMBL_FEATURE_UP:
      if (seg->data.s)
	strncat (look->segTextBuf, messprintf ("%s", seg->data.s), 64) ;
      break ;
    case HOMOL: case HOMOL_UP:
      hinf = arrp (look->homolInfo, seg->data.i, HOMOLINFO) ;
      meth = method (look->view, hinf->method) ;
      if (meth->flags & METHOD_BLASTN)
	strncat (look->segTextBuf, 
		 messprintf ("%s %.1f %d%% (%d - %d) ",
			    name (hinf->method),
			    hinf->score,
			    (int) (100 * (4 + hinf->score/ (seg->x2-seg->x1+1)) / 9.0),
			    hinf->x1, hinf->x2),
		 64) ;
      else
	strncat (look->segTextBuf, 
		 messprintf ("%s %.1f (%d - %d) ",
			    name (hinf->method),
			    hinf->score,
			    hinf->x1, hinf->x2),
		 64) ;
      if (!isGIF)		/* title () requires DB access */
	strncat (look->segTextBuf, title (seg->key), 
		 127 - strlen (look->segTextBuf)) ;
      break ;
    case FEATURE: case FEATURE_UP:
      meth = method (look->view, seg->key) ;
      if (meth->flags & METHOD_PERCENT)
	strcat (look->segTextBuf, messprintf ("%.0f%% ", seg->data.f)) ;
      else if (meth->flags & METHOD_SCORE)
	strcat (look->segTextBuf, messprintf ("%.3g ", seg->data.f)) ;
      if (seg->parent)
	strncat (look->segTextBuf, dictName (look->featDict, -seg->parent), 58) ;
      break ;
    case CODING:
      if (look->gf.cum && 
	  seg->x1 >= look->gf.min && seg->x2 < look->gf.max)
	strncat (look->segTextBuf,
		 messprintf ("Coding score %.2f",
			     look->gf.cum[seg->x2 + 1 - look->gf.min] -
			     look->gf.cum[seg->x1 - look->gf.min]), 64) ;
      break ;
    case EXON_SUPPORT: case EXON_SUPPORT_UP:
    case SPLICED_cDNA: case SPLICED_cDNA_UP:
     	fMapcDNAReportLine (look->segTextBuf,seg,256) ;
      break ;
    case TRANSCRIBEDGENE: case TRANSCRIBEDGENE_UP:
     	fMapcDNAReportTranscribedgene (look->segTextBuf,seg,256) ;
      break ;
    case TRANSCRIPT: case TRANSCRIPT_UP:
     	fMapcDNAReportTranscript (look->segTextBuf,seg,256) ;
      break ;
    case PMRNA: case PMRNA_UP:
     	fMapcDNAReportMrna (look->segTextBuf,seg,256) ;
    case MRNA: case MRNA_UP:
     	fMapcDNAReportMrna (look->segTextBuf,seg,256) ;
    case MGENES: case MGENES_UP:
     	fMapcDNAReportGenes (look->segTextBuf,seg,256) ;
      break ;
    case MRNAI: case MRNAI_UP:
     	fMapcDNAReportRNAi (look->segTextBuf,seg,256) ;
      break ;
    case MOST: case MOST_UP:
     	fMapcDNAReportOST (look->segTextBuf,seg,256) ;
      break ;
    case MPRODUCT: case MPRODUCT_UP:
     	fMapcDNAReportMProduct (look->segTextBuf,seg,256) ;
      break ;
    case VIRTUAL_MULTIPLET_TAG:  
      fMapSelectVirtualMultiplet (look, seg) ;
      break ;
    case CLONE_END:
      strncat (look->segTextBuf, name (seg->data.k), 64) ; /* left/right */
      break ;
    case PROBE: case PROBE_UP:
      fMapcDNAReportProbe (look->segTextBuf,seg,256) ;
      break ;
    case OLIGO: case OLIGO_UP:
      i =  (seg->data.i >> 16) & 0x3ff ;
      if (i > 0)
	strcat (look->segTextBuf, messprintf (" Tm: %3.1f ", (float) (i/10.0))) ;
      i =  seg->data.i & 0xfff ;
      if (i > 0)
	strcat (look->segTextBuf, messprintf (" Score: %d ", i)) ;
      break ;
    case OLIGO_PAIR: case OLIGO_PAIR_UP:
      fMapOspSelectOligoPair (look, seg) ;
      break ;
    case VISIBLE:
      if (iskey (seg->data.k))
	strcat (look->segTextBuf, messprintf (" (%s) ", name (seg->data.k))) ; 
      break ;
    default:		/* needed for picky compiler */
      break ;
    }
}

/**********************************************************/

static void fMapSelectBox (LOOK look, int box, double x, double y)
{ 
  int i, i0, max , imax, w=0 ;
  SEG *seg, *seg2, *seg3;
  KEY parent ;
  unsigned char *cp ;
  char *buf, *cq, *cq1 ;
  HOMOLINFO *hinf ;
  METHOD *meth ;
  static int finDna = 0 ;
  int ii, c1 ;
  BOOL isProt ;
  int nExon=0, tframe, oldTframe=0, diffTframe ;
  BOOL isTiret, calcWeight ;
  BOOL doComplementSurPlace = 
    (look->flag & FLAG_COMPLEMENT_SUR_PLACE) ;

  if (box >= arrayMax (look->boxIndex) || box < look->minLiveBox)
    return ;

  if (look->activeBox)
    { seg = BOXSEG (look->activeBox) ;
      parent = seg->parent ;
      for (i = look->minLiveBox ; i < arrayMax (look->boxIndex) ; ++i)
	{ seg2 = BOXSEG (i) ;
	  if (seg2 == seg || (parent && seg2->parent == parent))
	    switch (seg2->type)
	      {
	      case DNA_SEQ: case PEP_SEQ: case TRANS_SEQ: case TRANS_SEQ_UP:
		break ;
	      case HOMOL: case HOMOL_UP:
		hinf = arrp (look->homolInfo, seg2->data.i, HOMOLINFO) ;
		meth = method (look->view, hinf->method) ;
		graphBoxDraw (i, -1, meth->col) ;
		break ;
	      case FEATURE: case FEATURE_UP:
		meth = method (look->view, seg2->key) ;
		graphBoxDraw (i, -1, meth->col) ; 
		break ;
	      case ASSEMBLY_TAG:
		{
		  int color ;
		  unsigned char *c = (unsigned char*) seg2->data.s;
		  switch (c[0] ^ c[1] ^ c[2] ^ c[3])
		    {
		    case 'a'^'n'^'n'^'o': color = BLACK; break;
		    case 'c'^'o'^'m'^'m': color = BLUE; break;
		    case 'o'^'l'^'i'^'g': color = YELLOW; break;
		    case 'c'^'o'^'m'^'p': color = RED; break;
		    case 's'^'t'^'o'^'p': color = RED; break;
		    case 'r'^'e'^'p'^'e': color = GREEN; break;
		    case 'c'^'o'^'s'^'m': color = LIGHTGRAY; break;
		    case 'A'^'l'^'u'^' ': color = GREEN; break;
		    case 'i'^'g'^'n'^'o': color = LIGHTGRAY; break;
		    default:              color = LIGHTBLUE; break;
		    }
		  graphBoxDraw (i, BLACK, color) ;
		}
		break;
	      case INTRON: case INTRON_UP:
		if (arrp (look->seqInfo,seg2->data.i,SEQINFO)->flags & SEQ_CONFIRMED)
		  graphBoxDraw (i, -1, 
				arrp (look->seqInfo,seg2->data.i,SEQINFO)->flags & SEQ_CONFIRM_UTR ? YELLOW : GREEN) ;
		else
		  graphBoxDraw (i, -1, WHITE) ;
		break ;
	      case OLIGO: case OLIGO_UP:
		if (seg2->data.i & 0x4000)      /* old */
		  graphBoxDraw (i, -1, ORANGE) ;
		else if (seg2->data.i & 0x8000) /* selected */
		  graphBoxDraw (i, -1, GREEN) ;
		else
		  graphBoxDraw (i, -1, WHITE) ;
		break ;
	      default:
		if (assFind (look->chosen, SEG_HASH (seg2), 0))
		  graphBoxDraw (i, -1, GREEN) ;
		else if (assFind (look->antiChosen, SEG_HASH (seg2), 0))
		  graphBoxDraw (i, -1, LIGHTGREEN) ;		    
		else
		  graphBoxDraw (i, -1, WHITE) ;
		break ;
	      }
	}
      i = seg->x1 ;  /* mhmp 11.06.98 + 08.09.98 */   
      if (look->isRetro)
	i -=2 ;
      if (i <=  0) i = 0 ;
      if (!finDna)
	finDna = (seg->x2 >= look->length) ? look->length - 1 : seg->x2 ;
      if (look->colors)
	{
	  cp = arrp (look->colors, i, unsigned char) ;
	  for ( ; i <= finDna && i < arrayMax (look->colors) ; ++i)
	    *cp++ &= ~ (TINT_HIGHLIGHT1 | TINT_HIGHLIGHT2);
	}
      for (i = look->minLiveBox ; i < arrayMax (look->boxIndex) ; ++i)
	{ seg2 = BOXSEG (i) ;
	  if ((seg2->type == DNA_SEQ || seg2->type == PEP_SEQ || 
	       seg2->type == TRANS_SEQ || seg2->type == TRANS_SEQ_UP) && 
	      seg2->x1 <= finDna && seg2->x2 >= seg->x1)
	    graphBoxDraw (i, -1, -1) ;
	}
    }

  look->activeBox = box ;
  /* graphBubbleDisplay (box) ; */

  seg = BOXSEG (look->activeBox) ;
  if (seqDragBox)
    seg3 = BOXSEG (seqDragBox) ;
  else
    seg3 = seg ;
    /* record exact dna postition for gf add splice site functions */
  if (seg->type == DNA_SEQ)
    look->DNAcoord = seg->x1 + (int)x; 

  i = (seg->x1 > 0) ? seg->x1 : 0 ;
  max = (seg->x2 >= look->length) ? look->length - 1 : seg->x2 ;
  finDna = max ;
  if (seg->type == DNA_SEQ && seg3->type == DNA_SEQ)
    { 
      i = i + (int)x ;
      i0 = i ;
      imax = seg3->x1 + (int)newx ;
      if (i0 <= imax) /* mhmp  27.04.98 A VOIR: i0 <--> imax*/
	{
	  finDna = (seg3->x2 >= look->length) ? look->length - 1 : seg3->x2 ;
	  cp = arrp (look->colors, i, unsigned char) ;
	  for ( ; i <= imax ; ++i)
	    *cp++ |= TINT_HIGHLIGHT1;
	  /* mhmp 27.04.98 dna en lignes de 50 */
	  buf = messalloc ((imax - i0 + 2) *51/50) ; cq = buf ;
	  cq1 = arrp (look->dna, i0, char) ;
	  i = imax ; while (i-- >= i0 ) 
	    {
	       *cq++ =
		doComplementSurPlace ?
		 dnaDecodeChar[ (int)complementBase(*cq1++)] :
		dnaDecodeChar[ (int)*cq1++] ;
	      if ( ! ((imax - i) % 50) ) 
		*cq++ = '\n' ;
	    }
	  *cq++ = 0 ;
	  graphPostBuffer (buf) ;
	  messfree (buf) ;
	}
    }
  else  if (seg->type == PEP_SEQ || seg->type == TRANS_SEQ || seg->type == TRANS_SEQ_UP)
    { 
      i = i + 3* (int)x ;
      i0 = i ;
      imax = seg3->x1 + 3* (int)newx ;
      finDna = seg3->x2 + 2  ; /* mhmp  07.10.98 */
      finDna = (finDna >= look->length) ? look->length - 1 : finDna ;
      cp = arrp (look->colors, i, unsigned char) ;
      if (i0 <= imax) /* mhmp  27.04.98 A VOIR: i0 <--> imax*/
	{
	  if (i0 < imax)
	    for ( ; i <= imax + 2 ; ++i)
	      *cp++ |= TINT_HIGHLIGHT1;
	  else 
	    { 
	      cp = arrp (look->colors, seg->x1 + 3* (int)x, unsigned char) ;
	      if (seg->type == TRANS_SEQ_UP) /* mhmp 10.06.98 */
		{
		  *cp-- |= TINT_HIGHLIGHT1; 
		  *cp-- |= TINT_HIGHLIGHT2;
		  *cp-- |= TINT_HIGHLIGHT2;
		  look->isRetro = TRUE ;
		}
	      else
		{
		  *cp++ |= TINT_HIGHLIGHT1; 
		  *cp++ |= TINT_HIGHLIGHT2;
		  *cp++ |= TINT_HIGHLIGHT2;
		  look->isRetro = FALSE ;
		}
	    }
	  /* mhmp 27.04.98 proteines en lignes de 50 */
	  buf = messalloc ((imax - i0 + 2) *51/50) ; cq = buf ;
	  cq1 = arrp (look->dna, i0, char) ;
	  i = i0 ; 
	  c1 = COORD (look, i);
	  if (seg->type != PEP_SEQ)
	    {
	      for (ii=1; ii<arrayMax (look->maxDna); ii++)
		if (look->flag & FLAG_REVERSE)
		  {
		    if (c1 > array (look->maxDna, ii, int) &&
			c1 <= array (look->minDna, ii- 1, int))
		      {
			nExon = ii - 1;
			oldTframe = array (look->tframe, ii- 1, int) ;
			break ;
		      }
		  }
		else
		  {
		    if (c1 < array (look->maxDna, ii, int) &&
			c1 >= array (look->minDna, ii- 1, int))
		      {
			nExon = ii - 1;
			oldTframe = array (look->tframe, ii- 1, int) ;
			break ;
		      }
		  }
	      isTiret = FALSE ;
	      nbTrans = 0 ;
	      calcWeight = TRUE ;
	      weight = 18 ;
	      while (i <= imax ) 
		{
		  isProt = FALSE ;
		  c1 = COORD (look, i);
		  for (ii = nExon; ii<arrayMax (look->maxDna); ii++)
		    if (look->flag & FLAG_REVERSE)
		      {
			if (c1 > array (look->maxDna, ii, int) &&
			    c1 <= array (look->minDna, ii- 1, int))
			  { 
			    isProt = TRUE ;
			    isTiret = FALSE ;
			    break ;
			  }
		      }
		    else
		      if (c1 < array (look->maxDna, ii, int) &&
			  c1 >= array (look->minDna, ii- 1, int))
			{ 
			  isProt = TRUE ;
			  isTiret = FALSE ;
			  break ;
			}
		  if (isProt)
		    {
		      if (seg->type == TRANS_SEQ)
			*cq++ = doComplementSurPlace ? 
			  antiCodon (cq1) : codon (cq1) ;
		      else   /* TRANS_SEQ_UP */
			*cq++ = doComplementSurPlace ? 
			  codon (cq1) : antiCodon (cq1) ;
		      nbTrans++ ;
		      if (calcWeight)
			{
			  w = -1 ;
#ifdef JUNK
i do no longer understand the weights chosen
			  if (seg->type == TRANS_SEQ)
			    w = doComplementSurPlace ? 
			      molecularWeight[ (int) antiCodon (cq1)] :
			       molecularWeight[ (int) codon (cq1)] ;
			  else  /* TRANS_SEQ_UP */
			    w = doComplementSurPlace ? 
			      molecularWeight[ (int) codon (cq1)] :
			       molecularWeight[ (int) antiCodon (cq1)] ;
#endif
			}
		      if (w < 0)
			calcWeight = FALSE ;
		      else if (w == 0)
			{
			  weight = 0 ;
			  calcWeight = FALSE ;
			}
		      else
			weight = weight + w - 18 ;
		    }
		  else
		    if (!isTiret)
		      {
			nExon++ ;
			tframe =  array (look->tframe, nExon, int) ; 
			diffTframe = (tframe - oldTframe)%3 ;
			if (diffTframe < 0)
			  diffTframe +=3 ;
			i += diffTframe ;
			cq1 += diffTframe ;
			oldTframe = tframe ;
			isTiret = TRUE ;
		      }
		  cq1 += 3 ;
		  i += 3  ;	
		  if ( ! (nbTrans % 50) ) 
		    *cq++ = '\n' ;
		}
	      *cq++ = 0 ;
	      graphPostBuffer (buf) ;
	      messfree (buf) ;
	    }
	  else /* PEPTIDE */
	    {
	      i = imax ; 
	      while (i >= i0 ) 
		{
		  *cq++ = doComplementSurPlace ? 
		    antiCodon (cq1) : codon (cq1) ;
		  cq1 += 3 ;
		  i -= 3  ;	
		  if ( ! ((imax - i) % 50) ) 
		    *cq++ = '\n' ;
		}
	      *cq++ = 0 ;
	      graphPostBuffer (buf) ;
	      messfree (buf) ;
	      weight = 18 ;
	      cq1 = arrp (look->dna, i0, char) ;
	      i = imax ; 
	      while (i >= i0 ) 
		{
		  w = -1 ;
#ifdef JUNK
		  w = doComplementSurPlace ? 
		    molecularWeight[ (int) antiCodon (cq1)] :
		    molecularWeight[ (int) codon (cq1)] ;
#endif
		  if (w < 0)
		    break ;
		  else if (w == 0)
		    {
		      weight = 0 ;
		      break ;
		    }
		  else
		    weight = weight + w - 18 ;
		  cq1 += 3 ;
		  i -= 3 ;
		}
	    }
	}/* endif i0 */
    } 	/* endif seg->type */
  else if (arrayExists (look->colors) && i >= 0 && 
	   i < arrayMax (look->colors) && max < arrayMax (look->colors))
    { 
      cp = arrp (look->colors, i, unsigned char) ;
      for ( ; i <= max ; ++i)
	*cp++ |= TINT_HIGHLIGHT1;
    }

	/* write blue report line */
  if (! (look->flag & FLAG_HIDE_HEADER) && fmapView (look->view, _Fmap_Header))
    { fMapReportLine (look, seg, FALSE, x) ;
      graphBoxDraw (look->segBox, BLACK, PALEBLUE) ;
    }

  for (i = look->minLiveBox ; i < arrayMax (look->boxIndex) ; ++i)
    if (i != look->map->thumb.box)
      { seg2 = BOXSEG (i) ;
        if ((seg2->type == DNA_SEQ || seg2->type == PEP_SEQ || 
	     seg2->type == TRANS_SEQ || seg2->type == TRANS_SEQ_UP)  &&
	    seg2->x1 <= finDna && seg2->x2 >= seg->x1)
	  graphBoxDraw (i, -1, -1) ; 
      }

				/* colour siblings */
  if ((parent = seg->parent))
    for (i = look->minLiveBox ; i < arrayMax (look->boxIndex) ; i++)
      { if (i == look->map->thumb.box)
	  continue;
	seg2 = BOXSEG (i) ;
	if (seg2->parent == parent)
	  switch (seg2->type)
	    {
	    case DNA_SEQ: case PEP_SEQ: case TRANS_SEQ: case TRANS_SEQ_UP:
	      break ;
	    case HOMOL: case HOMOL_UP:
	      graphBoxDraw (i, -1, PALERED) ;
	      break ;
	    default:
	      if (assFind (look->chosen, SEG_HASH (seg2), 0))
		graphBoxDraw (i, -1, CYAN) ;
	      else if (assFind (look->antiChosen, SEG_HASH (seg2), 0))
		graphBoxDraw (i, -1, DARKRED) ;
	      else
		graphBoxDraw (i, -1, PALEBLUE) ;
	      break ;
	    }
      }

  switch (seg->type)	/* colour box itself over sibling colour */
    {
    case DNA_SEQ: case PEP_SEQ: case TRANS_SEQ: case TRANS_SEQ_UP:
      break ;
    case   VIRTUAL_SUB_SEQUENCE_TAG:
    case   VIRTUAL_SUB_SEQUENCE_TAG_UP:
      graphBoxDraw (box, -1, PALERED) ;
    default:
      if (assFind (look->chosen, SEG_HASH (seg), 0) ||
	  assFind (look->antiChosen, SEG_HASH (seg), 0))
	graphBoxDraw (box, -1, MAGENTA) ;
      else
	graphBoxDraw (box, -1, PALERED) ;
      break ;
    }
}

/**************************************/

static void fMapFollow (LOOK look, double x, double y)
{
  SEG* seg = BOXSEG (look->activeBox) ;
  KEY key = seg->key ;
  int table = class (key) ;

  if (!table || table == _VText)
    { key = seg->parent ;
      table = class (key) ;
    }

  if ((seg->type | 0x1) == SPLICED_cDNA_UP ||
      (seg->type | 0x1) == EXON_SUPPORT_UP ||
      (seg->type | 0x1) == TRANSCRIBEDGENE_UP )
    { if (fMapcDNAFollow (look->activeBox))
	return ;
    }
  else if ((seg->type | 0x1) == TRANSCRIBEDGENE_UP ||
	   (seg->type | 0x1) == MPRODUCT_UP)
    { display (seg->parent, look->seqKey, TREE) ;
      return ;
    }
  else if ((seg->type | 0x1) == PMRNA_UP || (seg->type | 0x1) == MRNA_UP)
    { 
      displayPreserve () ;
      display (seg->parent, look->seqKey, 0) ;
      return ;
    }
  else if (seg->type == VIRTUAL_SUB_SEQUENCE_TAG ||
	   seg->type == VIRTUAL_SUB_SEQUENCE_TAG_UP ||
	   seg->type == VIRTUAL_PARENT_SEQUENCE_TAG ||
	   seg->type == VIRTUAL_CONTIG_TAG)
    { if (fMapFollowVirtualMultiTrace (look->activeBox))
	return ;
    }
  else if (seg->type == SEQUENCE || seg->type == SEQUENCE_UP )
    {
      SEQINFO *sinf = arrayp (look->seqInfo, seg->data.i, SEQINFO) ;
      if (sinf->flags & SEQ_CANONICAL)
	{
#ifdef ACEMBLY
	  KEY map = keyGetKey (seg->key, str2tag("IntMap")) ;
	  if (map && keyFindTag (map, str2tag("Wiggle")))
	    {
	      float x1, y1, x2, y2 ;
	      KEY from ;

	      graphBoxDim (look->activeBox, &x1, &y1, &x2, &y2) ;
	      if (y2 - y1 > 10)
		{
		  if (look->flag & FLAG_COMPLEMENT)
		    from = - GRAPH2MAP(look->map,y+y1) + seg->x2 ; 
		  else
		    from = GRAPH2MAP(look->map,y+y1) - seg->x1 ; 
		  display (seg->key, from,  DtTiling) ;
		}
	      else
		display (key, look->seqKey, TREE) ;
	      return ;
	    }
#endif
	}
    }

  if (table == _VSequence)
    display (key, look->seqKey, TREE) ;
  else if (table)
    display (key, look->seqKey, 0) ;
}

/*************************************************************/
/************* a couple of general utilities *****************/

BOOL fMapFindSegBounds (LOOK look, SegType type, int *min, int *max)
{
  for (*min = 1 ; *min < arrayMax (look->segs) ; ++*min)
    if (arrp (look->segs, *min, SEG)->type == type)
      break ;
  for (*max = *min ; *max < arrayMax (look->segs) ; ++*max)
    if (arrp (look->segs, *max, SEG)->type != type)
      break ;
  return (*min < arrayMax (look->segs)) ;
}

float fMapSquash (float value, float midpoint)
{
  if (value < 0)
    value = 0 ;
  value /= midpoint ;
  return (1 - 1/ (value + 1)) ;
}

/****************************/

static void fMapCompHelp (void)
{
  graphMessage (
"Reverse-Complement: does what you expect (I hope!)\n"
"Complement: keep coordinates, view opposite strand,\n"
"    which runs bottom to top 5' to 3'.\n"
"Reverse: reverse coordinates, view the same strand,\n"
"    which now runs bottom to top 5' to 3'.\n"
" In both the last two cases the DNA display is still\n"
" read 5' to 3' left to right.  Doing both Complement\n"
" and Reverse is equivalent to Reverse-Complement.\n"
"Complement DNA in place: Don't change coords or\n"
"    direction.  DNA display swaps G-C, A-T.  Now\n"
"    the sequence reads 3' to 5' left to right.") ;
}

static void fMapComplementSurPlace (void)
{
  FMAPLOOKGET ("fMapComplementSurPlace") ;

  look->flag ^= FLAG_COMPLEMENT_SUR_PLACE ;
  fMapDraw (look, 0) ;
}

static void fMapReverse (void)
{
  FMAPLOOKGET ("fMapReverse") ;

  look->map->mag = -look->map->mag ;
  look->flag ^= FLAG_REVERSE ;

  look->origin = look->length + 1 - look->origin ;
  fMapDraw (look, 0) ;
}

void fMapRC (LOOK look)
{
  int top = look->length - 1 ;
  int i, tmp ;
  SEG *seg ;
  char *ci, *cj, ctmp ;
  float *ftmp ;

  if (look->dna && ! look->dnaR)
    {
      look->dnaR = arrayCopy (look->dna) ;
      reverseComplement (look->dnaR) ;
    }
  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs, i, SEG) ;
      tmp = seg->x1 ;
      seg->x1 = top - seg->x2 ;
      seg->x2 = top - tmp ;
      if (seg->type >= SEQUENCE && seg->type <= ALLELE_UP)
	seg->type ^= 1 ;
      if ((seg->type | 0x1) == SPLICED_cDNA_UP ||
	  (seg->type | 0x1) == EXON_SUPPORT_UP ||
	  (seg->type | 0x1) == TRANSCRIBEDGENE_UP ||
	  (seg->type | 0x1) == TRANSCRIPT_UP ||
	  (seg->type | 0x1) == PMRNA_UP ||
	  (seg->type | 0x1) == MRNA_UP ||
	  (seg->type | 0x1) == MPRODUCT_UP
	  )
       seg->source = top -  (seg->source + seg->sourceDx) ;
    }

  tmp = arrayMax (look->segs) ;
  arrayMax (look->segs) = look->lastTrueSeg ;
	/* avoid mixing temporary and permanent segs */
  arraySort (look->segs, fMapOrder) ;
  arrayMax (look->segs) = tmp ;

  if (look->dnaR)
    { Array tmp = look->dna ;
      look->dna = look->dnaR ; look->dnaR = tmp ;
    }
  else
    {
      ci = arrp (look->dna, 0, char) ; 
      cj = arrp (look->dna, top, char) ;
      while (ci <= cj)
	{ ctmp = *ci ;
	  *ci++ = complementBase(*cj) ;
	  *cj-- = complementBase(ctmp) ;
	}
    }

  if (look->gf.cum)
    { tmp = look->gf.max ;
      look->gf.max = top + 1 - look->gf.min ;
      look->gf.min = top + 1 - tmp ;
      ftmp = look->gf.cum ;
      look->gf.cum = look->gf.revCum ; /* already reversed */
      look->gf.revCum = ftmp ;
    }

  look->map->centre = top - look->map->centre ;
  tmp = look->zoneMin ;
  look->zoneMin = look->length - look->zoneMax ;
  look->zoneMax = look->length - tmp ;
  look->origin = look->length + 1 - look->origin ;

  look->flag ^= FLAG_COMPLEMENT ;
  oligoRepeatRC (look->oligoRepeats, look->length) ;
  fMapClearDNA (look) ;
}

static void fMapCompToggleTranslation (void)
{
  if ((mapColSetByName ("Up Gene Translation", -1)) ||
      (mapColSetByName ("Down Gene Translation", -1)))
    {
      mapColSetByName ("Up Gene Translation", 2) ;
      mapColSetByName ("Down Gene Translation", 2) ; 
    }
}

static void fMapComplement (void)
{
  FMAPLOOKGET ("fMapComplement") ;

  fMapCompToggleTranslation () ;
  fMapRC (look) ;
  fMapReverse ()  ;		/* does a lookDraw () */
}

static void fMapReverseComplement (void)
{
  FMAPLOOKGET ("fMapComplement") ;

  fMapCompToggleTranslation () ;
  fMapRC (look) ;
  fMapDraw (look, 0) ;
}

void fMapDoReverseComplement (LOOK look, BOOL complement)
{
  Graph old = graphActive () ;
  BOOL iscpl = look->flag & FLAG_COMPLEMENT ;

  if (complement != iscpl)
    { fMapRC (look) ;
      fMapDraw (look, 0) ;
    }
  graphActivate (old) ;
}

/**************************************************************/
/********* major routine to dump out in GFF format ************/

static char *cleanNewlines (const char *s, AC_HANDLE handle)
{ 
  int n = 0 ;
  const char *ccp ;
  char *cq, *copy ;

  for (ccp = s ; *ccp ; ++ccp)
    if (*ccp == '\n') ++n ;

  if (!n)
    return strnew (s, handle) ;

  copy = halloc (ccp-s+n+1, handle) ;
  for (ccp = s, cq = copy ; *ccp ; ++ccp)
    if (*ccp == '\n') 
      { *cq++ = '\\' ; *cq++ = 'n' ; }
    else
      *cq++ = *ccp ;
  *cq = 0 ;
  return copy ;
}

/***********************************/

static BOOL fMapDumpGFF (LOOK look, int version, 
			 KEY *refSeq, int *offPtr, BOOL isList, BOOL minSpan,
			 DICT *sourceSet, DICT *featSet)
{
  int i ;
  SEG *seg = 0 ;
  HOMOLINFO *hinf ;
  SEQINFO *sinf ;
  METHOD *meth ;
  KEY seqKey, sourceKey ;
  int x, y, type, offset = 0 ;
  int chopped ;			/* only used for version 1 */
  int tmp = 0, reversed = 0 ;                /* used by AcePerl */
  char *featName, *sourceName, *tagText = 0 ;
  float score ;
  char strand, frame ;
  BOOL flipped ;
  BOOL isScore ;
  AC_HANDLE handle = handleCreate () ;
  Associator key2sinf = assHandleCreate (handle) ;
  Associator key2source = assHandleCreate (handle) ;
  Associator key2feat = assHandleCreate (handle) ;
  Array stats = arrayHandleCreate (64, int, handle) ;
  DICT *listSet = 0 ;
  char timeBuf[25] ;

				/* first establish seqKey and offset */
  seqKey = 0 ;
  if (!minSpan && refSeq && !*refSeq && version > 1)
    *refSeq = look->seqOrig ;
  if (refSeq && *refSeq)
    { for (i = 1 ; i < arrayMax (look->segs) ; ++i)
        { seg = arrp (look->segs, i, SEG) ;
          if ((seg->type == SEQUENCE || seg->type == SEQUENCE_UP) &&
	      seg->key == *refSeq)
	    break ;
	}
      if (i < arrayMax (look->segs))
#if 0
	{ if (seg->type == SEQUENCE)
	    { seqKey = *refSeq ;
	      offset = - seg->x1 ;
	    }
	  else if (seg->type == SEQUENCE_UP)
	    { messout ("Can't GFF dump from a reversed reference sequence.  Sorry.") ;
	      return FALSE ;
	    }
	}
#endif
      { 
	if (seg->type == SEQUENCE || seg->type == SEQUENCE_UP) {
	  seqKey = *refSeq ;
	  offset = - seg->x1 ;
	  reversed = seg->type == SEQUENCE_UP ? 1 : 0;
	}
      }
      else
	{ messout ("Can't find reference sequence %s", name (*refSeq)) ;
	return FALSE ;
	}
    }
  if (!seqKey)			/* find minimal spanning sequence */
    { x = look->zoneMin + 1 ;
      y = look->zoneMax ;
      seqKey = 0 ; 
      if (!fMapFindSpan (look, &seqKey, &x, &y))
	seqKey = look->seqKey ;
      if (x > y)			/* what if reversed? */
	{ messout ("Can't GFF dump from reversed sequences for now.  Sorry.") ;
	  return FALSE ;
	}
      offset = x - (look->zoneMin + 1) ; /* x changed in fMapFindSpan */
    }
  if (refSeq)
    *refSeq = seqKey ;
  if (offPtr)
    *offPtr = offset ;

  if (!isList)
    { freeOutf ("##gff-version %d\n", version) ;
      freeOutf ("##date %s\n", timeShow (timeParse ("today"), timeBuf, 25)) ;
      freeOutf ("##sequence-region %s %d %d %s\n", name (seqKey), 
		look->zoneMin+1+offset, look->zoneMax+offset, reversed ? " (reversed)" : "") ;
    }

  if (isList)
    listSet = dictHandleCreate (64, handle) ;

  for (i = 0 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs, i, SEG) ;
      if (seg->x1 > look->zoneMax || seg->x2 < look->zoneMin)
	continue ;
      if (seg->type == MASTER ||
	  seg->type == VISIBLE ||
	  seg->type == CDS || seg->type == CDS_UP ||
	  seg->type == DNA_SEQ || seg->type == PEP_SEQ ||
	  seg->type == ORF || seg->type == TRANS_SEQ || seg->type == TRANS_SEQ_UP)
	continue ;

      if (seg->type & 0x01 &&
	  (seg->type >= SEQUENCE && seg->type <= ALLELE_UP))
	{ type = seg->type - 1 ;
	  x = seg->x2+1 ; y = seg->x1+1 ;
	}
      else
	{ type = seg->type ;
	  x = seg->x1+1 ; y = seg->x2+1 ;
	}

      chopped = 0 ;
      if (x <= y)
	{ flipped = FALSE ;
	  strand = '+' ;
	  if (version == 1 && x <= look->zoneMin)
	    chopped = look->zoneMin + 1 - x ;
	}
      else
	{ int tmp = x ; x = y ; y = tmp ;
	  flipped = TRUE ;
	  strand = '-' ;
	  if (version == 1 && y > look->zoneMax)
	    chopped = y - look->zoneMax ;
	}

      if (version == 1)		/* clip */
	{ if (x <= look->zoneMin) 
	    x = look->zoneMin + 1 ;
	  if (y > look->zoneMax)
	    y = look->zoneMax ;
	}

      x += offset ; y += offset ;

      frame = '.' ;		/* other defaults */
      isScore = (version == 1) ? TRUE : FALSE ;
      score = 0.0 ; 
      sourceName = (version == 1 || version == 2)  /* mieg, why not in version2 ? */
	? fMapSegTypeName[type] : 0 ;
      featName = 0 ;
      sourceKey = 0 ;
      sourceName = "";

      switch (type)
	{
	case SEQUENCE: case EXON: case INTRON:
	  sinf = arrayp (look->seqInfo, seg->data.i, SEQINFO) ;
	  if ( seg->key && (version == 2) ) {
	    featName = className (seg->key);
	  } else 
	    featName = (type == SEQUENCE) ? "sequence" : (type == EXON) ? "exon" : "intron" ;
	  if (sinf->method)
	    { 
	      sourceKey = sinf->method ;
	    }
	  if (type == SEQUENCE)
	    { assInsert (key2sinf, assVoid (seg->key), sinf) ; 
	      if (sinf->flags & SEQ_SCORE)
		{ score = sinf->score ;
		  isScore = TRUE ;
		}
	    }
	  break ;
	case CODING:
	  if (version == 1)
	    featName = "coding_exon" ;
	  else
	    featName = "CDS" ;
	  assFind (key2sinf, assVoid (seg->parent), &sinf) ;
	  if (sinf->method)
	    sourceKey = sinf->method ;
	  switch (seg->data.i%3)
	    { case 0: frame = '0' ; break ;
	      case 1: frame = '2' ; break ;
	      case 2: frame = '1' ; break ;
	    }
	  break ;
	case HOMOL:
	  hinf = arrp (look->homolInfo, seg->data.i, HOMOLINFO) ;
	  score = hinf->score ;
	  sourceKey = hinf->method ;
	  featName = "similarity" ;
	  meth = method (look->view, hinf->method) ;
	  if (! (meth->flags & METHOD_STRAND_SENSITIVE))
	    strand = '.' ;
	  if (meth->flags & METHOD_SCORE)
	    isScore = TRUE ;
	  if (meth->flags & METHOD_FRAME_SENSITIVE)
	    frame = '0' ;
	  break ;
	case FEATURE:
	case ATG:
	case SPLICE5:
	case SPLICE3:
	  sourceKey = seg->key ;
	  if (type == SPLICE3) featName = "splice3" ;
	  else if (type == SPLICE5) featName = "splice5" ;
	  else if (type == ATG) featName = "atg" ;
	  else
	    if (version == 1) featName = "feature" ;
	    else featName = "misc_feature" ;
	  meth = method (look->view, seg->key) ;
	  if (meth->flags & METHOD_SCORE)
	    { score = seg->data.f ;
	      isScore = TRUE ;
	    }
	  if (! (meth->flags & METHOD_STRAND_SENSITIVE))
	    strand = '.' ;
	  if (meth->flags & METHOD_FRAME_SENSITIVE)
	    frame = '0' ;
	  break ;
	case ASSEMBLY_TAG:
	  tagText = strnew (seg->data.s, handle) ;
	  sourceName = "assembly_tag" ;
	  featName = tagText ;
	  { char *cp ;
	    for (cp = tagText ; *cp && *cp != ':' && 
		   *cp != ' ' && *cp != '\t' && *cp != '\n' ; ++cp) ;
	    if (*cp)
	      { *cp++ = 0 ;
		while (*cp == ' ' || *cp == '\t' || *cp == '\n') ++cp ;
	      }
	    tagText = cp ;
	  }
	  break ;
	case ALLELE:
	  if (version == 1) sourceName = "variation" ;
	  if (seg->data.s)
	    { if (*seg->data.s == '-')
		featName = "deletion" ;
	      else
		{ char *cp = seg->data.s ;
		  featName = "variation" ;
		  while (*cp)
		    if (!strchr ("AGCTagct", *cp++))
		      featName = (version == 1) ? "insertion_site" : "insertion" ;
		}
	    }
	  break ;
	case EMBL_FEATURE:
	  featName = name (seg->key) ;
	  break ;
	case CLONE_END:
	  featName = name (seg->data.k) ;
	  strand = '.' ;
	  break ;
	default:
	  ;
	}

				/* finalise source and feat fields */
      if (sourceKey)
	{
	  if (assFind (key2source, assVoid (sourceKey), &sourceName))
	    assFind (key2feat, assVoid (sourceKey), &featName) ;
	  else
	    { 
	      OBJ obj = bsCreate (sourceKey) ;
	      sourceName = name (sourceKey) ;
	      if (obj)
		{
		  bsGetData (obj, str2tag ("GFF_source"), _Text, &sourceName) ;
		  if (bsGetData (obj, str2tag ("GFF_feature"), _Text, &featName))
		    assInsert (key2feat, assVoid (sourceKey), 
			       strnew (featName, handle)) ;
		  bsDestroy (obj) ;
		}
	      assInsert (key2source, assVoid (sourceKey), 
			 strnew (sourceName, handle)) ;
	    }
	}

      if (!sourceName)		/* only possible if version > 1 */
	sourceName = "unknown" ;

      if (!featName)		/* from method or from seg->type above */
	{
	  if (version == 1)
	    featName = sourceName ;
	  else
	    featName = fMapSegTypeName[type] ;
	}

      if (frame >= '0' && frame <= '2' && chopped)
	{ frame =  frame + 3 - (chopped % 3) ;
	  if (frame == '3') frame = '0' ;
	  if (frame == '4') frame = '1' ;
	  if (frame == '5') frame = '2' ;
	}

      if (sourceSet && featSet &&
	  (dictMax (sourceSet) || dictMax (featSet)) &&
	  !dictFind (sourceSet, sourceName, 0) &&
	  !dictFind (featSet, featName, 0))
	continue ;

      if (isList)		/* just accumulate stats */
	{ int k ;
	  dictAdd (listSet, messprintf ("%s\t%s",sourceName,featName), &k) ;
	  ++array (stats, k, int) ;
	  continue ;	/* move on to next line */
	}
				/* !isList from here on */
				/* write the main part of the line */

      /* LS/AcePerl: fixup reversed reference sequence */
      if (reversed) {  
	tmp = look->zoneMax + offset + 1 - y;
	y   = look->zoneMax + offset + 1 - x;
	x   = tmp;
	if ( strand == '+' )
	  strand = '-';
	else if ( strand == '-' )
	  strand = '+';
      }

      if (isScore)
	freeOutf ("%s\t%s\t%s\t%d\t%d\t%g\t%c\t%c",
		  name (seqKey), sourceName, featName, 
		  x, y, score, strand, frame) ;
      else 
	freeOutf ("%s\t%s\t%s\t%d\t%d\t.\t%c\t%c",
		  name (seqKey), sourceName, featName, x, y, strand, frame) ;

				/* extras */
      switch (type)
	{
	case SEQUENCE: case EXON: case CODING:
	  if (seg->parent)
	    {
	      if (version == 1) 
		freeOutf ("\t%s", name (seg->parent)) ;
	      else
		/* LS 6/11/99 this isn't always true in sMap world: 
		   freeOutf ("\tSequence \"%s\"", name (seg->parent)) ; 
		*/
		freeOutf ("\t%s \"%s\"",className (seg->parent) , name (seg->parent)) ;
	    }
	  break ;
	case INTRON:
	  sinf = arrayp (look->seqInfo, seg->data.i, SEQINFO) ;
	  if (version == 1)
	    { if (seg->parent)
	        { freeOutf ("\t%s", name (seg->parent)) ;
	          if (sinf->flags & SEQ_CONFIRMED)
		    freeOutf (" Confirmed") ;
		}
	      else if (sinf->flags & SEQ_CONFIRMED) freeOutf ("\tConfirmed") ;
	      if (sinf->flags & SEQ_CONFIRM_EST) freeOutf (" by_EST") ;
	      if (sinf->flags & SEQ_CONFIRM_CDNA) freeOutf (" by_cDNA") ;
	      if (sinf->flags & SEQ_CONFIRM_HOMOL) freeOutf (" by_homology") ;
	      if (sinf->flags & SEQ_CONFIRM_UTR) freeOutf (" in_UTR") ;
	    }
	  else
	    { BOOL isData = FALSE ; 
	      if (seg->parent)
	        {
		  /* LS 6/11/99 this isn't always true in sMap world: 
		     freeOutf ("\tSequence \"%s\"", name (seg->parent)) ; 
		  */
		  freeOutf ("\t%s \"%s\"",className (seg->parent) , name (seg->parent)) ;
		  isData = TRUE ;
		}
	      if (sinf->flags & SEQ_CONFIRM_EST) 
		{ freeOutf ("%s", isData ? " ; " : "\t") ; isData = TRUE ;
		  freeOutf ("Confirmed_by_EST") ;
		}
	      if (sinf->flags & SEQ_CONFIRM_CDNA)
		{ freeOutf ("%s", isData ? " ; " : "\t") ; isData = TRUE ;
		  freeOutf ("Confirmed_by_cDNA") ;
		}
	      if (sinf->flags & SEQ_CONFIRM_HOMOL)
		{ freeOutf ("%s", isData ? " ; " : "\t") ; isData = TRUE ;
		  freeOutf ("Confirmed_by_homology") ;
		}
	      if (sinf->flags & SEQ_CONFIRM_UTR)
		{ freeOutf ("%s", isData ? " ; " : "\t") ; isData = TRUE ;
		  freeOutf ("Confirmed_in_UTR") ;
		}
	    }
	  break ;
	case HOMOL:
	  if (version == 1)
	    freeOutf ("\t%s:%s", className (seg->key), name (seg->key)) ;
	  else
	    freeOutf ("\tTarget \"%s:%s\"", className (seg->key), name (seg->key)) ;
	  hinf = arrp (look->homolInfo, seg->data.i, HOMOLINFO) ;
	  if (flipped)
	    freeOutf (" %d %d", hinf->x2, hinf->x1) ;
	  else
	    freeOutf (" %d %d", hinf->x1, hinf->x2) ;
	  break ;
	case FEATURE:
	case SPLICE5:
	case SPLICE3:
	case ATG:
	  if (seg->parent)
	    {
	      if (version == 1)
		freeOutf ("\t%s", cleanNewlines (dictName (look->featDict, -seg->parent), handle)) ;
	      else
		freeOutf ("\tNote \"%s\"", cleanNewlines (dictName (look->featDict, -seg->parent), handle)) ;
	    }
	  break ;
	case ASSEMBLY_TAG:
	  if (tagText && *tagText) 
	    {
	      if (version == 1)
		freeOutf ("\t%s", cleanNewlines (tagText, handle)) ;
	      else
		freeOutf ("\tNote \"%s\"", cleanNewlines (tagText, handle)) ;
	    }
	  break ;
	case ALLELE:
	  if (version == 1)
	    freeOutf ("\t%s", name (seg->key)) ;
	  else
	    freeOutf ("\tAllele \"%s\"", name (seg->key)) ;
	  if (seg->data.s && *seg->data.s != '-')
	    {
	      if (version == 1)
		freeOutf ("\t%s", cleanNewlines (seg->data.s, handle)) ;
	      else if (!strcmp (featName, "variation"))
		freeOutf (" ; Variant \"%s\"", cleanNewlines (seg->data.s, handle)) ;
	      else
		freeOutf (" ; Insert \"%s\"", cleanNewlines (seg->data.s, handle)) ;
	    }
	  break ;
	case EMBL_FEATURE:
	  if (seg->data.s)
	    {
	      if (version == 1)
		freeOutf ("\t%s", cleanNewlines (seg->data.s, handle)) ;
	      else
		freeOutf ("\tNote \"%s\"", cleanNewlines (seg->data.s, handle)) ;
	    }
	  break ;
	case CLONE_END:
	  if (version == 1)
	    freeOutf ("\t%s", name (seg->key)) ;
	  else
	    freeOutf ("\tClone \"%s\"", name (seg->key)) ;
	  break ;
	default: ;
	}
      freeOutf ("\n") ;
    }

  if (isList)
    for (i = 0 ; i < arrayMax (stats) ; ++i)
      freeOutf ("%s\t%d\n", dictName (listSet, i), arr (stats,i,int)) ;

  handleDestroy (handle) ;
  return TRUE ;
}

/********* entry point to fMapDumpGFF for fmap menu *********/

static void dumpSegs (void)
{
  FILE *fil ;
  int level ;
  FMAPLOOKGET ("debugPrintSegs") ;

  if ((fil = filqueryopen (0, 0, "gff", "w", 
			   "File to write features into")))
    { level = freeOutSetFile (fil) ;
      fMapDumpGFF (look, 1, 0, 0, 0, 0, 0, 0) ;
      freeOutClose (level) ;
      filclose (fil) ;
    }
}

/********* entry point to fMapDumpGFF from dnacpt.c *********/

static BOOL findStartStop (KEY *seq, int *start, int *stop)
{
  OBJ Seq ;
  int x1, x2 ;
  KEY parent ;
  static Array s_children = 0 ;

  if (!*seq || ! (Seq = bsCreate (*seq)))
    return FALSE ;

  while (Seq)			/* recurse up */
    if (bsFindTag (Seq, str2tag ("S_Parent")) &&	/* tag2 system */
	!bsFindTag (Seq, str2tag ("Splicing")) &&
	bsGetKeyTags (Seq, _bsRight, 0) && /* a tag */
	bsGetKey (Seq, _bsRight, &parent))
      {
	int i ;

	bsDestroy (Seq) ;
				/* find the child in the parent */
	s_children = arrayReCreate (s_children, 128, BSunit) ;
	if (! (Seq = bsCreate (parent)) ||
	    !bsGetArray (Seq, str2tag ("S_Child"), s_children, 4))
	  break ;		/* stop the recursion */
	for (i = 0 ; i < arrayMax (s_children) ; i += 4)
	  if (arrp (s_children, i+1, BSunit)->k == *seq)
	    break ;
	 if (i >= arrayMax (s_children))
	   break ;		/* stop the recursion */

	 x1 = arrp (s_children, i+2, BSunit)->i ; 
	 x2 = arrp (s_children, i+3, BSunit)->i ; 
	 if (x1 < x2) {
	   if (!*stop)
	     { *start = x1 ; *stop = x2 ; }
	   else
	     { *start += x1-1 ; *stop += x1-1 ; }
	 }
	 else {
	   if (!*stop)
	     { *start = x2 ; *stop = x1 ; }
	   else
	     { int tmp = *start ;
	       *start = x1 - *stop + 1 ; *stop = x1 - tmp + 1 ; 
	     }
	 }
	 *seq = parent ;
      }
#ifdef ACEDB4
    else if (bsGetKey (Seq, _Source, &parent))
      { bsDestroy (Seq) ;
				/* find the child in the parent */
        if (! (Seq = bsCreate (parent)) ||
	    !bsFindKey (Seq, _Subsequence, *seq) ||
	    !bsGetData (Seq, _bsRight, _Int, &x1) ||
	    !bsGetData (Seq, _bsRight, _Int, &x2))
	  break ;
	if (x1 < x2)
	  {
	  if (!*stop)
	    { *start = x1 ; *stop = x2 ; }
	  else
	    { *start += x1-1 ; *stop += x1-1 ; }
	  }
	else
	  {
	  if (!*stop)
	    { *start = x2 ; *stop = x1 ; }
	  else
	    { int tmp = *start ;
	      *start = x1 - *stop + 1 ; *stop = x1 - tmp + 1 ; 
	    }
	  }
	*seq = parent ;
      }
#endif /* ACEDB4 */
    else
      break ;

  if (Seq)
    bsDestroy (Seq) ;

  return (*start != 1 || *stop != 0) ;
}

#ifdef OLD_4_7_CODE

static BOOL findStartStop (KEY *seq, int *start, int *stop)
{
  OBJ Seq ;
  int x1, x2 ;
  KEY parent ;
  Array dna = 0 ;

  if (!*seq || ! (Seq = bsCreate (*seq)))
    return FALSE ;

  *start = 1 ;
  *stop = 0 ;			/* go to end of sequence */
  while (bsGetKey (Seq, _Source, &parent)) /* recurse up */
    { bsDestroy (Seq) ;
      if (! (Seq = bsCreate (parent)) ||
	  !bsFindKey (Seq, _Subsequence, *seq) ||
	  !bsGetData (Seq, _bsRight, _Int, &x1) ||
	  !bsGetData (Seq, _bsRight, _Int, &x2))
	break ;
      if (x1 < x2)
	{
	  if (!*stop)
	    { *start = x1 ; *stop = x2 ; }
	  else
	    { *start += x1-1 ; *stop += x1-1 ; }
	}
      else
	{
	  if (!*stop)
	    { *start = x2 ; *stop = x1 ; }
	  else
	    { *start = x1 - *stop + 1 ; *stop = x1 - *start + 1 ; }
	}
      *seq = parent ;
    }
 
  if (*start == 1 && *stop == 0)
    if ((dna = dnaGet (*seq)))
      {	*stop = arrayMax (dna) ;
	arrayDestroy (dna) ;
      }

  if (Seq)
    bsDestroy (Seq) ;

  return (*start != 1 || *stop != 0) ;
}

#endif

void fMapDumpSegsKeySet (KEYSET kSet)
{
  int i ;
  LOOK look ;

  fMapInitialise () ;

  look = (LOOK) messalloc (sizeof (struct LookStruct)) ;
  look->origin = 0 ;
  look->zoneMin = 0 ;
  for (i = 0 ; i < keySetMax (kSet) ; ++i)
    { look->seqKey = keySet (kSet, i) ;
      look->start = 1 ;
      look->stop = 0 ;			/* find whole sequence */
      if (!findStartStop (&look->seqKey, &look->start, &look->stop))
	{ Array dna ;
	  if ((dna = dnaGet (look->seqKey))) /* try DNA */
	    { look->start = 1 ;
	      look->stop = arrayMax (dna) ;
	      arrayDestroy (dna) ;
	    }
	  else
	    continue ;		/* can't get limits */
	}
      --look->start ; --look->stop ;
      look->fullLength = look->length = look->stop - look->start + 1 ;
      look->zoneMax = look->length - 1 ;

      if (fMapConvert (look, TRUE))
	fMapDumpGFF (look, 1, 0, 0, 0, 0, 0, 0) ;
      arrayDestroy (look->segs) ;
    }

  messfree (look) ;
}

/****************/

static void fMapDumpAlign (LOOK look, BOOL isPep)
{
  int i, ii, bufMax ;
  char *cp, *cq ;
  Array dna = 0 ;
  
  SEG *seg ;
  HOMOLINFO *hinf = 0 ;
  KEY seqKey = look->seqOrig ;
  int a1, a2, x1, x2, x, y ;
  char *buf = 0 ;

  BOOL first = TRUE ;
  AC_HANDLE handle = handleCreate () ;

  x = look->zoneMin ;
  y = look->zoneMax ;
  bufMax = y - x + 1 ;
  if (isPep) bufMax /= 3 ;
  messfree (buf) ;
  buf = messalloc (bufMax + 1) ;
  buf[bufMax] = 0 ;
				/* extras */
  for (ii = 0 ; ii < arrayMax (look->segs) ; ++ii)
    { seg = arrp (look->segs, ii, SEG) ;
      if (seg->x1 > look->zoneMax || seg->x2 < look->zoneMin)
	continue ;

      switch (seg->type)
	{
	case HOMOL:
	  if (first)
	    {
	      first = FALSE ;
	      dna = dnaCopy (look->dna) ;
	      x = look->zoneMin ;
	      y = look->zoneMax ;  
	      memset (buf, (int)'.', bufMax) ;
	      freeOutf ("%s/%d-%d ", name (seqKey), 
			isPep ? (x + 1 - look->origin)/3 : x + 1 - look->origin , 
			isPep ? (y - look->origin)/3  : y - look->origin) ;
	      if (x > 0 && x < y && y < arrayMax (dna))
		for (i = x , cp = arrp (dna, i, char), cq = buf ; i < y ; i++, cp++, cq++)
		  {
		    if (isPep)
		      { *cq = codon (cp) ; i +=2 ; cp += 2 ; }
		    else
		      *cq = dnaDecodeChar[ (int)*cp] ;
		  }
	      freeOutf ("%s %s\n", buf, name (seqKey)) ;
	      arrayDestroy (dna) ;
	    }
	  dna = dnaGet (seg->key) ;
	  if (dna && arrayMax (dna))
	    {
	      hinf = arrp (look->homolInfo, seg->data.i, HOMOLINFO) ;
	      if (hinf->x1 > hinf->x2)   /* flipped) */
		{
		}
	      else
		{
		  x1 = hinf->x1 ;  x2 = hinf->x2 ;
		  a1 = seg->x1 ; a2 = seg->x2 ;
		  i = look->zoneMin - a1 ;
		  if (i > 0) { a1 += i ; x1 += i ; }
		  i = seg->x2 - look->zoneMax ;
		  if (i > 0) { a2 -= i ; x2 -= i ; }
		  if (x1 >= 0 && x1 < x2 && x2 < arrayMax (dna))
		    {
		      freeOutf ("%s/", name (seg->key)) ;
		      freeOutf ("%d-%d ", isPep ? x1/3 : x1, isPep ? x2/3 : x2 ) ;
		      
		      memset (buf, (int)'.', bufMax) ;
		      for (i = a1 , cp = arrp (dna, x1 - 1, char) , 
			     cq = (isPep ? buf + (a1 - look->zoneMin)/3 : buf + a1 - look->zoneMin) ; 
			   i <= a2 ; i++, cp++, cq++)
			if (isPep)
			  { *cq = codon (cp) ; i += 2 ; cp += 2 ; }
			else
			  *cq = dnaDecodeChar[ (int)*cp] ;
		      freeOutf ("%s %s\n",buf, name (seg->key)) ;
		    }
		}
	      arrayDestroy (dna) ;
	    }
	  break ;
	default: ;
	}
    }
  messfree (buf) ;
  
  handleDestroy (handle) ;
}

/*************************************/

static void geneStats (void)
{
#if !defined (MACINTOSH)
  int i ;
  SEG *seg ;
  int len = 0 ;
  int nIntron = 0 ;
  int nExon = 0 ;
  int nCDS = 0 ;
  int intronLength = 0 ;
  int exonLength = 0 ;
  int cdsLength = 0 ;
  int nLoci = 0 ;
  int ncDNA = 0 ;
  int nMatch = 0 ;
  char *n ;
  FMAPLOOKGET ("geneStats") ;

  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs, i, SEG) ;
      n = name (seg->parent) ;
      if (n[strlen (n)-1] == 'a')
	continue ;
      switch (seg->type)
	{
	case INTRON: case INTRON_UP:
	  ++nIntron ;
	  intronLength += (seg->x2 - seg->x1 + 1) ;
	  break ;
	case EXON: case EXON_UP:
	  ++nExon ;
	  exonLength += (seg->x2 - seg->x1 + 1) ;
	  break ;
	case CDS: case CDS_UP:
	  ++nCDS ;
	  cdsLength += (seg->x2 - seg->x1 + 1) ;
	  printf ("  %s\n", name (seg->parent)) ;
	  break ;
	case VISIBLE:
	  printf ("    %s: %s %s\n", name (seg->parent), 
		  name (seg->data.k), name (seg->key)) ;
	  if (seg->data.k == _Locus) ++nLoci ;
	  if (seg->data.k == _Matching_cDNA) ++ncDNA ;
	  if (seg->data.k == _Brief_identification) ++nMatch ;
	  break ;
	default:		/* needed for picky compiler */
	  break ;
	}
    }
  for (i = 0, n = arrp (look->dna,0,char) ; i < arrayMax (look->dna) ; ++i)
    if (*n++)
      ++len ;

  printf ("%d bp sequenced DNA from %d in LINK\n", len, look->length) ;
  printf ("%d introns, total %d bp\n", nIntron, intronLength) ;
  printf ("%d exons, total %d bp\n", nExon, exonLength) ;
  printf ("%d CDS (predicted genes), total %d bp\n", nCDS, cdsLength) ;
  printf ("%d database matches\n", nMatch) ;
  printf ("%d cDNAs\n", ncDNA) ;
  printf ("%d Genetic loci\n", nLoci) ;
#else
  messout ("Sorry, this does nothing on a macintosh "
           " (you aren't missing much - it was a quick hack)") ;
#endif
}


#ifndef ACEMBLY			/* trace, assembly stuff */
void fMapTraceDestroy (LOOK look) { ;}
void fMapShowTraceAssembly (LOOK look, float *offset) { ; }
BOOL fMapFollowVirtualMultiTrace (int box) { return FALSE ; }
void abiFixFinish (void) { ; }
void fMapGelDisplay (void) { gelDisplay () ; }
void fMapTraceFindMultiplets (LOOK look) { ; }
void fMapSelectVirtualMultiplet (LOOK look, SEG *seg) { ; } 
#endif

/****************** giface hooks *****************/

void* fMapGifGet (void *opaqueLook)
{
  KEY key, classKey, classSeq, view = 0 ;
  int x1, x2, origin ;
  char *word ;
  LOOK look = (LOOK) opaqueLook ;
        
  fMapInitialise () ;

  if (look)
    { fMapGifLook = look ; mapGifMap = look->map ;
      fMapDestroy () ; /* could check for possible reuse */
      fMapGifLook = 0 ; mapGifMap = 0 ;
    }

  if (! (word = freeword ()))
    goto usage ;

  lexword2key ("Sequence", &classSeq, _VClass);
  if (strcmp (word, "-class")) {
    classKey = classSeq ;
  } else {
    if (! (word = freeword ()) || !lexword2key (word, &classKey, _VClass)
	|| ! (word = freeword ())) 
      goto usage;
  }
  
  if (!classKey) {
      freeOut ("// gif sequence error : class not found\n") ;
      goto usage ;
  }
    
  if (!lexword2key (word, &key, classKey))
    { 
      freeOutf ("// gif sequence error : %s \"%s\" not known\n", name (classKey), word) ;
      goto usage ;
    }

  if (classKey != classSeq && ! keyFindTag (key, _DNA) && !  keyFindTag (key, str2tag ("SMAP"))) {
    freeOutf ("// gif sequence error : %s \"%s\" contains no DNA or SMAP\n", name (classKey), word) ;
    goto usage ;
  }

  x1 = 1 ; x2 = 0 ;		/* default for whole sequence */
  origin = UT_NON_INT ;
  while (freecheck ("w"))
    { 
      word = freeword () ;
      if (!strcmp (word, "-coords"))
	{ 
	  if (!freeint (&x1) || !freeint (&x2) )
	    goto usage ;
	}
      else if (!strcmp (word, "-origin"))
	{ 
	  if (!freeint (&origin) )
	    goto usage ;
	}
      else if (!strcmp (word, "-view"))
	{
	  if (! (word = freeword ()))
	    { 
	      freeOut ("// Error: -view must be followed by a view name\n") ;
	      goto usage ;
	    }
	  if (!lexword2key (word, &view, _VView))
	    {
	      freeOutf ("// Error: unknown view %s\n", word) ;
	      goto usage ;
	    }
	}
      else
	goto usage ;
    }
  if (origin == UT_NON_INT)
    origin = 1 ;

  mapGifMap = (MAP) 1 ;		/* so mapCreate doesn't access the graph */
  if (x1 < x2 || (x1 == 1 && x2 == 0))
    look = fMapCreateLook (key, x1, x2, origin, TRUE, view, 0) ;
  else
    look = fMapCreateLook (key, x2, x1, origin + 2, TRUE, view, 0x1) ;
  mapGifMap = 0 ;
				/* set the zone */
  if (look)
    {
      look->view = view ;
      if (x1 == 1 && x2 == 0)
	{ int i ;
	  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
	    if (arrp (look->segs,i,SEG)->key == key)
	      { setZone (look, arrp (look->segs,i,SEG)->x1, 
			 arrp (look->segs,i,SEG)->x2 + 1) ;
	      break ;
	      }
	}
      else if (x2 < x1)
	fMapRC (look) ;
      return look ;
    }
  else
    {
      freeOutf ("// gif warning, probably this sequence has no associated dna\n") ;
      return 0 ;
    }
 usage:
  freeOutf ("// gif seqget error: usage: seqget [-class class] sequence [-coords x1 x2] [-origin x0] \n") ;
  return 0 ;
}

/**************************************************************/

void fMapGifDestroy (void *opaqueLook)
{
  LOOK look = (LOOK) opaqueLook ;

  if (look)
    { fMapGifLook = look ; mapGifMap = look->map ;
      fMapDestroy () ;
      fMapGifLook = 0 ; mapGifMap = 0 ;
    }
} 

/**************************************************************/

void fMapGifDisplay (void *opaqueLook)
{
  int x1, x2 ;
  char *word ;
  LOOK look = (LOOK) opaqueLook ;
        
  fMapInitialise () ;

  fMapGifLook = look ; mapGifMap = look->map ;

  if (look->flag & FLAG_REVERSE)
    { x1 = COORD (look, look->zoneMax) + 1 ;
      x2 = COORD (look, look->zoneMin) ;
    }
  else
    { x1 = COORD (look, look->zoneMin) ;
      x2 = COORD (look, look->zoneMax) - 1 ;
    }

  while ((word = freeword ()))
    {
      if (*word != '-')
	{ freeback () ; break ; }
      if (!strcmp (word, "-visible_coords"))
	{ 
	  if (!freeint (&x1) || !freeint (&x2) || (x2 == x1))
	    goto usage ;
	  if (x1 > x2)
	    { int t = x1 ; x1 = x2 ; x2 = t ; }
	}
      else if (!strcmp (word, "-view")) /* -D displaytype */
	{
	  KEY view = 0 ;
	  if (! (word = freeword ()))
	    { 
	      freeOut ("// Error: -view must be followed by a view name\n") ;
	      goto usage ;
	    }
	  if (!lexword2key (word, &view, _VView))
	    {
	      freeOutf ("// Error: unknown view %s\n", word) ;
	      goto usage ;
	    }
	  look->view = view ;
	}
    }
  
  if (! (look->graph = displayCreate (FMAP)))
    { fMapGifLook = 0 ; mapGifMap = 0 ; return ; }

  topMargin = 0 ;
  graphFitBounds (&mapGraphWidth,&mapGraphHeight) ;
  halfGraphHeight = 0.5* (mapGraphHeight - topMargin) ;

  look->map->min = 0 ;
  look->map->max = look->length ;
  look->map->centre = look->origin + (x1+x2) / 2.0 ;
  look->map->mag = (mapGraphHeight  - topMargin) / (1.05 * (x2 - x1)) ;
  if (x1 > x2)
    { 
      look->map->mag = -look->map->mag ;
      fMapRC (look) ; 
      /*  look->origin = look->origin - look->length - 1 ; mieg ??? coords are wrong */
    }

  if (FALSE && !getenv ("ACEDB_PROJECT") && (look->zoneMax - look->zoneMin < 1000))
    { 
      mapColSetByName ("DNA Sequence", TRUE) ;
      mapColSetByName ("Coords", TRUE) ;
    }    
  if (look->translateGene == 1) /* established for class _VmRNA in fMapCreateLook */
    {  
      look->translateGene = keyGetKey (look->seqKey, str2tag ("Product")) ;
      /*      mapColSetByName ("DNA Sequence", 1) ; */
      mapColSetByName ("Down Gene Translation", 1) ;
    }

  fMapDraw (look, 0) ;
  fMapSelect (look) ; 

  { float d = (mapGraphHeight - topMargin - 5.0) / (look->map->mag * 1.05) ;
    float xx2 = look->map->centre + 0.5*d ;
    x1 = (int) (xx2 - d + 0.5) ;
    x2 = (int) (xx2 + 0.5) ;
    if (0)  /* this kills the svg exportation */
      {
	freeOutf ("// FMAP_COORDS %s %d %d", 
		  freeprotect (name (look->seqKey)), x1, x2) ;
	freeOutf (" -visible_coords %d %d\n", COORD (look,x1), COORD (look, x2)) ;
      }
  }

  fMapGifLook = 0 ; mapGifMap = 0 ; 
  return ;

usage:
  freeOutf ("// gif sequence error: usage: SEQUENCE sequence [-coords x1 x2]\n") ;
  fMapGifLook = 0 ; mapGifMap = 0 ; 
}

/**************************************************************/

void fMapGifActions (void* opaqueLook)
{
  LOOK look = (LOOK) opaqueLook ;

  fMapGifLook = look ; mapGifMap = look->map ;

  while (freecheck ("w"))
    { char *word = freeword () ;
      if (!strcmp (word, "-dna"))
	fMapToggleDna () ;
      else if (!strcmp (word, "-gf_features"))
	fMapAddGfSegs () ;
      else if (!strcmp (word, "-hide_header"))
	hideHeaderToggle () ;
      else if (!strcmp (word, "-rev_comp"))
	fMapReverseComplement () ;
      else
	goto usage ;
    }
  if (look->map)
    fMapDraw (look, 0) ;

  fMapGifLook = 0 ; mapGifMap = 0 ; 
  return ;

usage:
  freeOutf ("// gif seqactions error: usage: SEQACTIONS [-dna] [-gf_features] [-hide_header] [-rev_comp] [-visible_coords x1 x2] [-view view] \n") ;
  fMapGifLook = 0 ; mapGifMap = 0 ; 
}

/**************************************************************/

void fMapGifColumns (void* opaqueLook)
{
  LOOK look = (LOOK) opaqueLook ;
  int state = -1 ;

  if (!look->map) return ;

  fMapGifLook = look ; mapGifMap = look->map ;

  fMapInitialise () ;

  while (freecheck ("w"))
    { char *word = freeword () ;
      if (!strcmp (word, "-on"))
	state = 1 ;
      else if (!strcmp (word, "-off"))
	state = 0 ;
      else if (state > -1)
	mapColSetByName (word, state) ;
      else
	goto usage ;
    }

  fMapGifLook = 0 ; mapGifMap = 0 ; 
  return ;

usage:
  freeOutf ("// gif seqcolumns error: usage: SEQCOLUMNS  {-on name} {-off name}\n") ;
  fMapGifLook = 0 ; mapGifMap = 0 ; 
}

/**************************************************************/

void fMapGifAlign (void* opaqueLook)
{
  int x1, x2, pepFrame = 0 ;
  int level = 0 ;
  FILE *fil = 0 ;
  char *word ;
  BOOL isPep = FALSE ;
  LOOK look = (LOOK) opaqueLook ;

  fMapInitialise () ;

  fMapGifLook = look ; mapGifMap = look->map ; 

  while (freecheck ("w"))
    { word = freeword () ;
      if (!strcmp (word, "-coords"))
	{ if (!freeint (&x1) || !freeint (&x2) || (x2 == x1))
	    goto usage ;
	  setZoneUserCoords (look, x1, x2) ;
	}
      else if (!strcmp (word, "-file"))
	{ if (!freecheck ("w") || ! (fil = filopen (freeword (), 0, "w")))
	    goto usage ;
	  level = freeOutSetFile (fil) ;
	}
      else if (!strcmp (word, "-peptide"))
	{ isPep = TRUE ;
	}
      else
	goto usage ;
    }

  pepFrame = (pepFrame + 999) % 3 ;
  fMapDumpAlign (look, isPep) ;

  if (level) freeOutClose (level) ;
  if (fil) filclose (fil) ;	/* should be same condition as level */

  fMapGifLook = 0 ; mapGifMap = 0 ; 
  return ;

usage:
  if (level) { freeOutClose (level) ; filclose (fil) ; }
  freeOutf ("// gif seqaligns error: usage: SEQALIGN [-coords x1 x2] [-file fname] \n") ;
  fMapGifLook = 0 ; mapGifMap = 0 ;
}

/**************************************************************/

static void addToSet (DICT *dict)  /* little utility for parsing below */
{ 
  char cutter, *cutset, *word ;

  if (freestep ('"'))
    cutset = "|\"" ;
  else
    cutset = "| \t" ;
  while ((word = freewordcut (cutset, &cutter)))
    { dictAdd (dict, word, 0) ;
      if (cutter != '|') break ;
    }
}

/******************/

void fMapGifFeatures (void* opaqueLook)
{
  KEY refseq = 0 ;
  int x1, x2, offset ;
  int level = 0 ;
  int version = 1 ;
  FILE *fil = 0 ;
  char *word ;
  AC_HANDLE handle = handleCreate () ;
  DICT *sourceSet = dictHandleCreate (16, handle) ;
  DICT *featSet = dictHandleCreate (16, handle) ;
  BOOL isList = FALSE ;
  BOOL minSpan = FALSE ;
  BOOL result ;
  LOOK look = (LOOK) opaqueLook ;

  fMapInitialise () ;

  fMapGifLook = look ; mapGifMap = look->map ; 

  while (freecheck ("w"))
    { word = freeword () ;
      if (!strcmp (word, "-coords"))
	{ if (!freeint (&x1) || !freeint (&x2) || (x2 == x1))
	    goto usage ;
	  setZoneUserCoords (look, x1, x2) ;
	}
      else if (!strcmp (word, "-file"))
	{ if (!freecheck ("w") || ! (fil = filopen (freeword (), 0, "w")))
	    goto usage ;
	  level = freeOutSetFile (fil) ;
	}
      else if (!strcmp (word, "-refseq"))
	{ if (!freecheck ("w") || ! (lexword2key (freeword (), &refseq, _VSequence)))
	    goto usage ;
	}
      else if (!strcmp (word, "-version"))
	{ if (!freeint (&version))
	    goto usage ;
	}
      else if (!strcmp (word, "-list"))
	isList = TRUE ;
      else if (!strcmp (word, "-minspan"))
	minSpan = TRUE ;
      else if (!strcmp (word, "-source") && freecheck ("w"))
	addToSet (sourceSet) ;
      else if (!strcmp (word, "-feature") && freecheck ("w"))
	addToSet (featSet) ;
      else
	goto usage ;
    }

  result = fMapDumpGFF (look, version, &refseq, &offset, 
			isList, minSpan, sourceSet, featSet) ;

  if (level) freeOutClose (level) ;
  if (fil) filclose (fil) ;	/* should be same condition as level */

  if (result) 
    freeOutf ("// FMAP_FEATURES %s %d %d\n", freeprotect (name (refseq)), 
	      look->zoneMin+1+offset, look->zoneMax+offset) ;
  fMapGifLook = 0 ; mapGifMap = 0 ; 
  handleDestroy (handle) ;
  return ;

usage:
  if (level) { freeOutClose (level) ; filclose (fil) ; }
  freeOut ("// gif seqfeatures error: usage: SEQFEATURES [-coords x1 x2] ") ;
  freeOut ("[-file fname] [-refseq sequence] [-version 1|2] ") ;
  freeOut ("[-list] [-minspan] [-source source (s)] [-feature feature (s)]\n") ;
  handleDestroy (handle) ;
  return ;
}

/**************************************************************/

void fMapGifDNA (void* opaqueLook)
{
  KEY key ;
  int x1, x2 ;
  int level = 0 ;
  FILE *fil = 0 ;
  char *word ;
  LOOK look = (LOOK) opaqueLook ;

  fMapInitialise () ;

  fMapGifLook = look ; mapGifMap = look->map ;

  while (freecheck ("w"))
    { word = freeword () ;
      if (!strcmp (word, "-coords"))
	{ if (!freeint (&x1) || !freeint (&x2) || (x2 == x1))
	    goto usage ;
	  setZoneUserCoords (look, x1, x2) ;
	}
      else if (!strcmp (word, "-file"))
	{ if (!freecheck ("w") || ! (fil = filopen (freeword (), 0, "w")))
	    goto usage ;
	  level = freeOutSetFile (fil) ;
	}
      else
	goto usage ;
    }

  key = 0 ; x1 = look->zoneMin+1 ; x2 = look->zoneMax ;
  if (!fMapFindSpan (look, &key, &x1, &x2)) key = look->seqKey ;

  dnaDumpFastA (look->dna, look->zoneMin, look->zoneMax-1,
		messprintf ("%s/%d-%d", name (key), x1, x2), 0, 0) ;
  
  if (level) freeOutClose (level) ; 
  if (fil) filclose (fil) ;
  freeOutf ("// FMAP_DNA %s %d %d\n", freeprotect (name (key)), x1, x2) ;
  fMapGifLook = 0 ; mapGifMap = 0 ;
  return ;

usage:
  if (level) { freeOutClose (level) ; filclose (fil) ; }
  freeOutf ("// gif seqdna error: usage: SEQDNA [-coords x1 x2] [-file fname]\n") ;
  fMapGifLook = 0 ; mapGifMap = 0 ;
}

/**************************************************************/
/**************************************************************/
/**************************************************************/

static void fMapZoneKeySet (void)
{  
  SEG *seg ;
  int i, n = 0 ;
  KEYSET ks = keySetCreate () ;
  FMAPLOOKGET ("zoneKeySet") ;

  if (0) 
    { /* keep compiler happy */
      fMapCompHelp () ;
      fMapComplement () ;
      fMapComplementSurPlace () ;
      fMapToggleTranslation () ;
      geneStats () ; /* keep compiler happy */
    }

  for (i = 0 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs, i, SEG) ;
      if (seg->x1 > look->zoneMax || seg->x2 < look->zoneMin)
	continue ;
      if (class (seg->key))
	keySet (ks, n++) = seg->key ;
    }
  keySetSort (ks) ;
  keySetCompress (ks) ;

  if (!keySetMax (ks))
    { messout ("The active zone does not contain any key") ;
      keySetDestroy (ks) ;
    }
  else 
    keySetNewDisplay (ks, messprintf ("%s %d %d", 
				     name (look->seqKey),look->zoneMin - look->origin + 1, 
				     look->zoneMax - look->origin)) ;
}

/****************** end of file *****************/
