/*  File: fmapcdna.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
   Display cDNA
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 15 15:21 1998 (fw)
 * Created: Thu Dec  9 00:01:50 1993 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */

#include "fmap_.h"

#include "main.h"
#include "display.h"

#include "bs.h"
#include "a.h"
#include "call.h"
#include "systags.h"
#include "classes.h"
#include "display.h"
#include "dna.h"
#include "peptide.h"
#include "dnaalign.h"
#include "makemrna.h"
#include "bump.h"
#include "query.h"
#include "session.h"
#include "cdna.h"
#include "acembly.h"
#include "dotter.h"
#include "parse.h"
#include "annot.h"
#include "pick.h"
#include "bindex.h"

static int numIntronCol = 12 ;
static KEY intronCol[] = { PALEYELLOW , PALEGREEN, PALECYAN,   PALEMAGENTA ,  PALEVIOLET , 
			   RED2 , PALEBLUE, PALEGRAY, GREEN2, RED1, BLUE1, PALEORANGE,   } ;

#ifdef ACEMBLY
static ANNOTATION  notes[] = 
{
  /*   { 1, "Ready for submission",  "Ready_for_submission",  0, 'k', 0, {0}}, */
  { 1, "Fully edited",  "Fully_edited",  0, 'k', 0, {0}},
  { 1, "Remark", "Remark", 0, _Text, 1000, {0}},
  { 1, "Sequencing error in genomic DNA",  "Sequencing_error_in_genomic_dna",  0, _Text, 1000, {0}},
  /*
  { 4,    "No Gap", "No_Gap", 0, 'k', 0, {0}},
  { 4,    "Gap in single open exon", "Single_gap", 0, 'k', 0, {0}},
  { 4,    "Gap", "Gaps", 0, 'k', 0, {0}},
  { 4,    "Single clone", "Single_clone", 0, 'k', 0, {0}},
  */
  { 4,    "5-prime", 0, 0, 'k', 0, {0}},
  /*
  { 8,          "SL1",  "SL1", 0, 'k', 0, {0}},
  { 8,          "SL2",  "SL2", 0, 'k', 0, {0}},
  { 8,          "gccgtgctc",  "gccgtgctc", 0, 'k', 0, {0}},
  */
  { 8,          "Non standard 5p motif",  "5p_motif", 0, _Text, 1000, {0}},
  { 8,          "Clone with 5p inversion",  "5p_inversion", 0, 'k', 0, {0}},
  /*
  { 8,          "Many clones start here",  "Many_clones_start_here", 0, 'i', 0, {0}},
  { 8,          "5p alternatif",  "Alternative_5prime", 0, 'k', 0, {0}},
  */
  { 4,          "Alternative splicing",  "Alternative_splicing", 0, 'k', 0, {0}},
  /*
  { 12,             "Coupled to alt 3p",  "coupled_to_alt_3p", 0, 'i', 0, {0}},
  { 12,             "Alternative exon",  "Alternative_exon", 0, 'k', 0, {0}},
  { 12,             "Overlapping alternative exons",  "Overlapping_alternative_exons", 0, 'k', 0, {0}},
  { 12,             "3_6_9 splice variations",  "3_6_9_splice_variations", 0, 'k', 0, {0}},
  */
{ 8,          "Confirmed non gt_ag",  "Confirmed_non_gt_ag", 0, 'k', 0, {0}},
  { 4,    "3-prime", 0, 0, 'k', 0, {0}},
  /*   { 8,        "Poly A seen in 5p read",  "Poly_A_seen_in_5p_read", 0, 'k', 0, {0}}, */
  { 8,        "Alternative poly A",  "Alternative_poly_A", 0, 'k', 0, {0}},
  { 8,        "No Alternative poly A",  "No_Alternative_poly_A", 0, 'k', 0, {0}},
  { 8,        "Fake_internal poly A",  "Fake_internal_poly_A", 0, 'k', 0, {0}},
  /*
 { 1, "Translation", 0, 0, 'k', 0, {0}},
  { 4,    "Open from first to last exon",  "Open_from_first_to_last_exon", 0, 'k', 0, {0}},
  { 4,    "Nb introns after stop", "Nb_introns_after_stop", 0, _Int, 0, {0}},
  */
  { 1, "Problem",  0, 0, 'k', 0, {0}},
  { 4,    "Jack Pot Bug",  "Jack_Pot_Bug", 0, 'k', 0, {0}},
  { 4,    "Fuse Bug",  "Fuse_Bug", 0, 'k', 0, {0}},
  { 4,    "Other Program Bug",  "Other_Bug", 0, _Text, 1000, {0}},
  { 1, "Genefinder", 0, 0, 'k', 0, {0}},
  { 4,    "Compatible", "Compatible_prediction", 0, 'k', 0, {0}},
  { 4,    "Different", "Different", 0, 'k', 0, {0}},
  /*   { 4,    "Corresponds to one predicted gene", "Matching_genefinder_gene", 0, 'k', 0, {0}}, */
  { 4,    "No predicted gene", "No_prediction", 0, 'k', 0, {0}},
  /*
  { 4,    "cDNA covers several predicted genes", "cDNA_covers_several_predicted_genes", 0, 'k', 0, {0}},
  { 4,    "Predicted covers several genes", "Predicted_covers_several_genes", 0, 'k', 0, {0}},
*/
  { 0, 0, 0, 0, 0, 0, {0}}
} ;
#endif

static void fMapSplicedcDNADecorate (LOOK look, float x, SEG *seg, BOOL isUp, int x1, int x2) ;
static Array findSpliceErrors (KEY link, KEY key, Array dnaD, Array dnaR, int a1, int a2, 
				 int x1, int x2, BOOL upSequence) ;
static Array spliceCptErreur (Array dnaD, Array dnaR, Array dnaCourt, BOOL isUp,
             int x1, int x2, int x3, 
             int *clipTopp, int *cEndp, int *cExp) ;


static void fMapcDNAAction (KEY key, int box) ;
static KEY myTranslatedGene = 0 ;
static FREEOPT fMapcDNAMenu[] = {
#ifdef ACEMBLY
  {  21, "cDNA menu"},
  {'c', "Show cDNA clone"},
  {'t', "Show EST"},
  /* {'d', "Show EST Sequence"}, */
  /* {'s', "Show EST Trace"}, does not work well */
  {'w', "Change Strand"},
  {'g', "Discard to trash"},
  {'e', "Discard from here"},
  {'v', "Ignore this clone"},
  {'u', "Ignore duplicate clone"},
  {'o', "Fuse to clone"},
  {'L', "Split as 1 clone/read"},
  {'R', "Restrict to active zone"},
  {'f', "Suspected internal deletion"},
  {'m', "Suspected mosaic"},
  {'6', "Resequence gap"},
  {'7', "Internal priming"},
  {'8', "Internal capping"},
  {'9', "Anomalous length"},
  {'A', "Possibly unspliced"},
  {'B', "Unspliced non coding"},
  {'5', "Real 5' end"},
  {'3', "Real 3' end"},
  {'C', "Comment"},
  /*   {'p', "Dot plot"} useless */
#else
  {  2, "cDNA menu"},
  {'c', "Show cDNA clone"},
  {'t', "Show EST"}
  /* {'d', "Show EST Sequence"}  not useful */
#endif
} ;

static void fMapProductAction (KEY key, int box) ;
static FREEOPT fMapProductMenu[] = {
  {7, "Product menu"},
  {'w', "Who am I"},
  {'m', "Show mRNA"},
  {'U', "Recompute title"},
  {'t', "Translate"},
  {'T', "3 frame Translate on/off"},
  {'g', "Genome View"},
  {'f', "Fiche-display"},  
} ;

static void fMapTranscribedgeneAction (KEY key, int box) ;
static FREEOPT fMapTranscribedgeneMenu[] = {
  {13, "T_Gene menu"},
  {'w', "Who am I"},
  {'t', "Translate"},
  {'s', "Recompute Splicing"},
  {'S', "Limit this gene to the active zone"},
  {'T', "Shed weakly connected variants"},
  {'F', "Fuse genes present in the active zone"},
  {'E', "Fuse and search repeats in active zone"},
  {'J', "Fuse and search rubber repeats in active zone"},
  /*   
       {'R', "Revert to original base call"},
       {'r', "Rename this gene"},
       {'c', "Sequin"}
  */
  {'k', "Eliminate this gene"},
  {'P', "Change strand"} ,
  {'G', "Create genebox"},
  {'H', "Absorb in genebox"},
  {'U', "Split geneBox by geneID"}
  /*
    {'a', "Annotate this gene"},
    {'b', "New annotator"}
  */
} ;

/* uses same fMapTranscribedgeneAction */
static FREEOPT fMapTranscriptMenu[] = {
  { 3, "T_Gene menu"},
  {'w', "Who am I"},
  /*   {'B', "New annotator"}, */
  {'D', "Show cDNA is new window"},
  {'K', "Elimintate this transcript"},
} ;

static void fMapMrnaAction (KEY key, int box) ;
static FREEOPT fMapMrnaMenu[] = {
  { 4, "mRNA menu"},
  {'w', "Who am I"},
  {'D', "Show cDNA is new window"},
  {'t', "Translate"},
  {'T', "3 frame Translate on/off"},
} ;

static void fMapGenesAction (KEY key, int box) ;
static FREEOPT fMapGenesMenu[] = {
  { 8, "Genes menu"},
  {'w', "Who am I"},
  {'R', "Rename, keep position name"},
  {'r', "Rename"},
  {'k', "Kill"},
  {'f', "RNA editing"},
  {'u', "Use AM"},
  {'s', "Resize as tg"},
  {'t', "Resize as genefinder"}
} ;

static FREEOPT fMapRnaiMenu[] = {
  { 1, "Rnai menu"},
  {'w', "Who am I"},
} ;

/***********************************************************************/

void fMapcDNAInit (void)
{
  static BOOL firstPass = TRUE ;
  
  if (!firstPass)   
    cDNAAlignInit () ;

  firstPass = FALSE ;
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***************************************************************************************/
/***********************************************************************/

static BOOL isTooLong (SEG *seg)
{
  KEY mrna = 0, clone = 0 ;
  OBJ Mrna = 0 ;
  BOOL tooLong = FALSE ;
  int nn = 0 ;

  if (seg->key)
    {
      mrna = keyGetKey (seg->key, _In_mRNA) ;
      clone = keyGetKey (seg->key, _cDNA_clone) ;
    }
  if (mrna)
    Mrna = bsCreate (mrna) ;
  if (Mrna && 
      bsFindKey (Mrna, _cDNA_clone, clone) &&
      bsGetData (Mrna, _bsRight, _Int, 0) &&
      bsGetData (Mrna, _bsRight, _Text, 0) &&
      bsGetData (Mrna, _bsRight, _Int, 0) &&
      bsGetData (Mrna, _bsRight, _Text, 0) &&
      bsGetData (Mrna, _bsRight, _Float, 0) &&
      bsGetData (Mrna, _bsRight, _Text, 0) &&
      bsGetData (Mrna, _bsRight, _Int, &nn))
    {
      if (nn > 30 || nn < -30)
	tooLong = TRUE ;
    }
  bsDestroy (Mrna) ;
  
  return tooLong ;
}

static BOOL isAmbiguous (SEG *seg)
{
  return seg->data.i &  0x80000000 ?  TRUE : FALSE ;
}

/***********************************************************************/
/***********************************************************************/
/* set bit 31 for ambiguous, 30 polyA, unaligned end, 28 unaligned begin */

BOOL fMapcDNADoFillData (SEG *seg)
{
  KEY dummy, cDNA_clone = 0, parentClone = 0 ;
  int x, xx, xx1 = 0, xx2 = 0, v1 ; 
  OBJ obj2 = 0;
  BOOL isComplete = FALSE ;
  KEY _Fully_sequenced = 0 ;
  static BSMARK mark = 0 ;
  static KEY _Parent_clone = 0 ;

  if (! _Parent_clone)
    _Parent_clone = str2tag ("Parent_clone") ;
  lexword2key ("Fully_sequenced", &_Fully_sequenced, 0) ;
  cDNAAlignInit () ;

  seg->data.i &= 0x0FFFFFFF ;
  if (keyFindTag(seg->key, _PolyA_after_base) ||
      keyFindTag(seg->key, _cDNA_clone))
    {
      obj2 = bsCreate (seg->key) ;
      if (obj2)  
	{ 
	  isComplete = FALSE ;
	  if (bsGetKey (obj2, _cDNA_clone, &cDNA_clone))
	    {
	      seg->parent = KEYKEY(cDNA_clone) ; /* needed for bumping */
	      isComplete = keyFindTag (cDNA_clone, _Fully_sequenced) ;
	      if ((parentClone = keyGetKey (cDNA_clone, _Parent_clone)))
		seg->parent = parentClone ;	      
	    }
	  xx1 = seg->data.i & 0x3FFF ; xx2 =  ((seg->data.i >> 14) & 0x3FFF) ; 

	  if (bsGetData (obj2, _PolyA_after_base, _Int, &xx) &&
	      xx == xx1)   /* before eventual swapping */
	    seg->data.i |=  0x40000000 ; 

	  if (xx1 > xx2) { int xx0 = xx2 ; xx2 = xx1 ; xx1 = xx0 ; } /* swap */
	  
	  xx = 1 ;
	  if (bsGetKey (obj2, _Transpliced_to, &dummy))
	    do 
	      {  /* get longest annotated transplicing */
		mark = bsMark(obj2, mark);
		if (bsGetData (obj2, _bsRight, _Int, &x) &&
		    x > xx) xx = x ;
	      	bsGoto(obj2, mark);
	      } while (bsGetKey (obj2, _bsDown, &dummy)) ;

	  if (xx == 1)
	    bsGetData (obj2, _Vector_clipping, _Int, &xx) ;
	  if (
	      ( bsFindTag (obj2, _Forward) && xx + (isGifDisplay ? 30 : 5) < xx1) ||
	      ( bsFindTag (obj2, _Reverse) && xx + (isGifDisplay ? 40 : 10) < xx1) 
	      )
	    seg->data.i |=  0x10000000 ;   /* up incomplete */	   
 
	  if ( !isComplete &&
	      ( 
	       (bsGetData (obj2, _Vector_clipping, _Int, &v1) && 
		bsGetData (obj2, _bsRight, _Int, &xx)) ||
	       (bsGetKey (obj2, _DNA, &dummy) && 
		bsGetData (obj2, _bsRight, _Int, &xx))
	       ) &&
	      (xx - 40 > xx2)
	      )
	    seg->data.i |=  0x20000000 ;    /* down incomplete */	    	    
	  if (bsGetKey (obj2, _From_gene, &dummy) && 
	      bsGetKey (obj2, _bsDown, &dummy))  /* at least 2 genes */
	    seg->data.i |=  0x80000000 ;   
	  bsDestroy (obj2) ;
	}
    }
  return TRUE ;
}

/*********************************************************************/
/* economize function calls in the seg loop */

#define cDNAFillDataZZ(__seg) (( (TRUE || (__seg)->data.i & 0x80000000)) ? TRUE : fMapcDNADoFillData(__seg))
#define cDNAFillData(__seg) (( (__seg)->data.i & 0x80000000) ? TRUE : fMapcDNADoFillData(__seg))

/*********************************************************************/

void fMapcDNAReportLine (char *buffer, SEG *seg1, int maxBuf)
{
  int x1, x2, y1, y2, ln = 0, polyA = 0, composite = 0 ; ;
  SEG *seg ;
  OBJ obj = 0 ;
  float mlength = 0 ;
  KEY dnaKey, cdna_clone = 0  ;
  char *cp = strnew(buffer,0) ;

  x1 = x2 = seg1->data.i >> 14 & 0x3FFF ; /* start equal, since order unknown */
  seg = seg1 ;
  while (seg1->type == seg->type && seg1->key == seg->key)
    {
      y1 = seg->data.i >> 14 & 0x3FFF ;
      y2 = seg->data.i & 0x3FFF ;
      if (x1 > y1) x1 = y1 ;
      if (x2 < y2) x2 = y2 ;
      if (x1 > y2) x1 = y2 ;
      if (x2 < y1) x2 = y1 ;
      seg++ ;
    }
  seg = seg1 - 1 ;
  while (seg1->type == seg->type && seg1->key == seg->key)
    {
      y1 = seg->data.i >> 14 & 0x3FFF ;
      y2 = seg->data.i & 0x3FFF ;
      if (x1 > y1) x1 = y1 ;
      if (x2 < y2) x2 = y2 ;
      if (x1 > y2) x1 = y2 ;
      if (x2 < y1) x2 = y1 ;
      seg-- ;
    }
  if ((obj = bsCreate(seg1->key)))
    {
      int n = 0, n1 = 0 ;
      if (bsGetKey(obj,_DNA,&dnaKey))
	bsGetData (obj, _bsRight, _Int, &ln) ;
      bsGetData (obj, _PolyA_after_base, _Int, &polyA) ; 
      bsGetKey (obj, _cDNA_clone, &cdna_clone) ;

      if (bsGetData (obj, _Composite, _Int, &n) && n > 0)
	{
	  if (bsGetData (obj, _bsRight, _Int, &n1) && n1 > 0 &&
	      bsGetData (obj, _bsRight, _Int, &n1) && n1 > n
	      )
	    n = n1 ;
	}
      composite = n ;
      bsDestroy(obj) ;
    }

  if (cdna_clone && (obj = bsCreate(cdna_clone)))
    {
      if (bsGetData (obj, _PCR_product_size, _Float, &mlength))
	mlength -= .200 ;
      bsDestroy(obj) ;
    }

  if (cdna_clone && strcmp (name(seg1->key), name(cdna_clone) ))
    strncpy (buffer, messprintf ("cDNA Clone %s ", name(cdna_clone)), maxBuf) ;
  if (composite > 0)
    strncat (buffer, messprintf ("#%d ", composite), maxBuf) ;
  if (mlength > 0)
    strncat (buffer, messprintf ("(%g kb) ", mlength), maxBuf) ;
  if (0) strncat (buffer, cp, maxBuf) ;
  strncat (buffer, messprintf ("ln %d, Match %d to %d", ln, x1, x2), maxBuf) ;
  if (polyA)
    strncat (buffer, messprintf ("PolyA at %d", polyA, x2), maxBuf) ; 
  messfree (cp) ;
}

/*********************************************************************/

void fMapShowGeneWalls (LOOK look, float *offset)
{
  SEG *seg ;
  int i, ii, x1 = 0, x2, xg1, xg2 ;
  float y ;
  KEY cosmid ;
  OBJ Cosmid = 0 ;
  BOOL isDownStrand, isDown ;
  BSunit *u ;
  Array aa = 0 ;
  static int xScale = 0 ; /* memorize between up and down 
			     should realy be the offset of the scale
			     */

  if (*look->map->activeColName == '-')
    {
      xScale = *offset + 2 ;
      isDownStrand = FALSE ;
    }
  else
    {
      xScale = xScale ? xScale + 6 : *offset + 2 ;
      isDownStrand = TRUE ;
    }


  graphColor (BLUE) ;
  xScale = 20 ;

  for (ii = 1 ; ii < arrayMax(look->segs) ; ++ii)
    { seg = arrp(look->segs,ii,SEG) ;
      switch (seg->type)
	{
	case SEQUENCE: isDown = TRUE ; break ;
	case SEQUENCE_UP: isDown = FALSE ; break ;
	default: continue ;
	}
      cosmid = seg->key ;
      if (!keyFindTag (cosmid, _Gene_wall))
	continue ;

      Cosmid = bsCreate (cosmid) ;
      if (!Cosmid)
	continue ;
      aa = arrayReCreate (aa, 36, BSunit) ; 

      if (bsGetArray (Cosmid,_Gene_wall, aa, 2))   
	for (i = 0 ; i < arrayMax(aa) ; i += 2)
	  { 
	    u = arrp (aa, i, BSunit) ;
	    if (isDown)
	      { 
		x1 = seg->x1 + u[0].i - 1 ;
		x2 = seg->x1 + u[1].i - 1 ;
	      }
	    else
	      {
		x1 = seg->x2 - u[0].i ;
		x2 = seg->x2 - u[1].i ;
	      }
	    if (isDownStrand)
	      {
		if (x1 > x2)
		  continue ;
		if (x2 >= x1 + 10)
		  graphColor (CYAN) ;
		xg1 = xScale ;
		xg2 = *offset ; /* so i do not over use the empty screen */
	      }
	    else
	      {
		if (x1 < x2)
		  continue ;
		if (x1 >= x2 + 10)
		  graphColor (CYAN) ;
		xg1 = 1 ;
		xg2 = xScale ;
	      }
	    y = MAP2GRAPH(look->map,x1) ;
	    if (y > .2 + topMargin &&
		y < mapGraphHeight - .2)
	      graphLine (xg1, y, xg2, y) ;
	    graphColor (BLUE) ;
	  }
      bsDestroy (Cosmid) ;
  }
  graphColor (BLACK) ;
  arrayDestroy (aa) ;
}

/*********************************************************************/

static void fMapShowSequencingErrors (LOOK look, float *offset)
{
  SEG *seg ;
  int i, ii, xx, oldColor ;
  float y ;
  KEY cosmid ;
  OBJ Cosmid = 0 ;
  static Array e2a = 0 ;
  static KEY _error = 0, _error2 ;
  BOOL isDown ;
  BSunit *uu ;

  if (!_error) 
    {
      _error = str2tag ("Possible_error_at_base") ;    
      _error2 = str2tag ("Possible_genomic_error") ;
    }

  if (*look->map->activeColName == '-')
    return ;

  for (ii = 1 ; ii < arrayMax(look->segs) ; ++ii)
    { seg = arrp(look->segs,ii,SEG) ;
      switch (seg->type)
	{
	case SEQUENCE: isDown = TRUE ; break ;
	case SEQUENCE_UP: isDown = FALSE ; break ;
	default: continue ;
	}
      cosmid = seg->key ;
      if (!keyFindTag (cosmid, _error) && !keyFindTag (cosmid, _error2))
	continue ;

      Cosmid = bsCreate (cosmid) ;
      if (!Cosmid)
	continue ;
      oldColor = graphColor (RED) ;
      if (bsGetData (Cosmid,_error, _Int, &xx))
	do
	  {
	    if (isDown) xx = seg->x1 + xx - 1 ;
	    else xx = seg->x2 - xx ;

	    y = MAP2GRAPH(look->map,xx) ;
	    if (y > .2 + topMargin &&
		y < mapGraphHeight - .2)
	      graphLine (1, y,  30, y) ;
	  } while (bsGetData (Cosmid, _bsDown, _Int, &xx)) ; 

      graphColor (ORANGE) ;
      if (bsFindTag (Cosmid, _error2))
	{
	  if (!e2a) e2a = arrayCreate(12, BSunit) ;
	  if (bsGetArray (Cosmid,_error2, e2a, 2))
	    for (i = 0, uu = arrp(e2a,0,BSunit) ; i < arrayMax(e2a); i+= 2, uu += 2)
	      {
		xx = uu[1].i ;
		if (isDown) xx = seg->x1 + xx - 1 ;
		else xx = seg->x2 - xx ;
		
		y = MAP2GRAPH(look->map,xx) ;
		if (y > .2 + topMargin &&
		    y < mapGraphHeight - .2)
		  graphLine (1, y + .5,  *offset < 30 ? *offset : 20, y + .5) ;
	      }
	}
      bsDestroy (Cosmid) ;
      graphColor (oldColor) ;
  }
  graphColor (BLACK) ;
}

/*********************************************************************/

void fMapcDNAGeneName (LOOK look, float *offset)
{
  BOOL found = FALSE ;
  char *cp = 0, *cq ;
  float y = 0, y1 = 0, y2 = 0 ;
  int 
    box = 0, dy,
    ii, x1, x2, width = 28,   /* width of the text */
    widthUsed = 0 ;
  SEG *seg ;

   if (*look->map->activeColName == '-')
    { 
      return ;
    }

  for (ii = 1 ; ii < arrayMax(look->segs) ; ++ii)
    { 
      seg = arrp(look->segs,ii,SEG) ;
      switch (seg->type)
	{
	case CDNA_GENE_NAME: 
	case CDNA_GENE_NAME_UP: break ;
	default: continue ;
	}
      if (seg->x2 < look->min) 
	continue ;
      if (seg->x1 > look->max) /* do not break here, because we mix gn and gn_UP in same loop */
	continue ;
      if (!seg->data.i)
	continue ;
      x1 = seg->x1 ; x2 = seg->x2 ;
      if (x1 < look->min)
	x1 = look->min ;
      if (x2 > look->max)
	x2 = look->max ;
      if (x1 > x2)
	continue ;
      y1 = MAP2GRAPH(look->map,x1) ;
      y2 = MAP2GRAPH(look->map,x2+1) ;	/* to cover full base */
      
      if (y2 > y1)
	{
	  if (stackExists(look->geneNamesStack) &&
	      seg->data.i > 0 &&
	      (cp = stackText (look->geneNamesStack, seg->data.i)) &&
	      strlen (cp))
	    {
	      box = graphBoxStart() ;
	      dy = uLinesText (cp, width) ;
	      y = (y1 + y2 - dy)/2 ;
	      if (y < y1) y = y1 ;
	      while ((cq = uNextLine(cp)))
		{
		  if (widthUsed < strlen (cq))
		    widthUsed = strlen (cq) ;
		  graphText (cq, *offset + 2, -.35 + y++) ;
		  if (y > y2 + 1) break ;
		}
	      graphBoxEnd () ;
	      graphBoxInfo (box, seg->key, cp) ;
	    }
	  found = TRUE ;
	}
    }
  if (found)
    *offset += widthUsed + 3 ;
}

/*********************************************************************/
/* this function is necessary to be able to toggle decoration */
void fMapcDNADecorateSplicedcDNA (LOOK look, float *offset)
{
  return ;
}

/***************************************************************************************/
/* type 0:exon, 1: intron  , 2: pair , 3: unaligned top */
static void estBubbleInfo (int box, SEG *seg, int origin, int type)
{
  KEY est, clone ;

  est = seg->key ;
  clone = keyGetKey (est, _cDNA_clone) ;

  if (keyFindTag (seg->key, _RefSeqMaker))
    graphBubbleInfo (box, 0, 0, "AceView reference sequence reconstructed from all aligned cDNAs and ESTs with best match to the genome.") ;
  else if (type == 2)
    graphBubbleInfo (box, clone, "cDNA_clone", 
		     messprintf("Paired reads\nsequence %s"
				, name(est))
		     ) ;
  else if (type == 3)
    graphBubbleInfo (box, clone, "cDNA_clone", 
		     messprintf("PolyA visible in\nsequence %s"
				, name(est))
		     ) ;
  else if (type == 4)
    graphBubbleInfo (box, clone, "cDNA_clone", 
		     messprintf("Unaligned bases in\nsequence %s"
				, name(est))
		     ) ;
  else if (type == 0 || type == 1)
    {
      char support [80] ;
      KEY tissue = keyGetKey (est, str2tag ("Tissue")) ;
      BOOL inverted = keyFindTag (est, str2tag ("Inverted")) ;
      BOOL pa = keyFindTag (est, str2tag ("polyA_after_base")) ;
      BOOL pp = keyFindTag (clone, str2tag ("Internal_priming_on_A_rich")) ;
      BOOL pp1 = keyFindTag (clone, str2tag ("Internal_priming")) ;
      BOOL deleted = (
		      (keyFindTag (clone, str2tag ("Suspected_internal_deletion"))  || 
		      keyFindTag (clone, str2tag ("Manual_internal_deletion")) ) &&
		      ! keyFindTag (clone, str2tag ("Manual_no_internal_deletion"))) ;
      KEYSET capping = queryKey (est, "Transpliced_to = SL0") ;
      KEYSET pairs = queryKey (clone, ">read ; from_gene") ;
      BOOL isPair = keySetMax (pairs) > 1 ;
      BOOL isCapped = keySetMax (capping) ;
      BOOL isRefSeq = keyFindTag (est, str2tag("Ref_seq")) ;
      char bf1[80] ;

      sprintf (bf1, "%d bp%s"
	       , seg->x2 - seg->x1 + 1
	       , type == 1 ? " Intron\n" : "\n") ;

      support [0] = 0 ;
      if (est && keyFindTag (est, _Composite))
	{
	  OBJ Est = bsCreate (est) ;
	  int n = 0, n1 = 0 ;
	  if (bsGetData (Est, _Composite, _Int, &n) && n > 0)
	    {
	      if (bsGetData (Est, _bsRight, _Int, &n1) && n1 > 0 &&
		  bsGetData (Est, _bsRight, _Int, &n1) && n1 > n
		  )
		n = n1 ;
	    }
	    sprintf (support, "Seen %d times", n) ;
	  bsDestroy (Est) ;
	}
      graphBubbleInfo (box, clone, "cDNA_clone", 
		       messprintf("%s%s%s%s%s%s%s%s%s%s%s%s"
				  , tissue ? name(tissue) : ""
				  , tissue ? "\n" : ""
				  , bf1
				  , isRefSeq ? "NCBI RefSeq " : "Sequence " 
				  , name(est)
				  , isCapped ? ", capped clone" : ""
				  , isPair ? ", paired reads" : ""
				  , inverted ?  ", submitted on opposite strand" : ""
				  , deleted ? ", insert partially deleted" : ""
				  , pp ? ", primed in A rich region of mRNA or genome" : ""
				  , pp1 ? ", internal priming" : ""
				  , pa ? ", polyA included" : ""
				  , support[0] ? support : ""
				  )) ;
      keySetDestroy (capping) ;
      keySetDestroy (pairs) ;
    }
}


void fMapcDNAShowSplicedcDNAorExon (LOOK look, float *offset, BOOL showExonSupport)
{
  int ii, box = -1, x1, x2; KEY bgCol ;
  SEG *seg ;
  float x = 0, y1, y2, oldy2 = 0, oldWidth=0, xMax = 0 ;
  BOOL is3p, isReadUp ;
  BOOL isDown ; 
  KEY bumpedParent = 0 ;
  BUMP bump = 0 ; 
  float ytop= 0, ybottom, yz ;
  int ix = 0 ;
  BSMARK mark = 0 ;
  int MY_MAGENTA = LIGHTGRAY ;
  KEY _aForward = 0, _aReverse = 0 ;
  BOOL foundCompositeColour = FALSE ;

  if (*offset > mapGraphWidth)
    return ;

  if (*look->map->activeColName == '-')
    { 
      isDown = FALSE ;
      xMax = *offset + 1.7 * 6 ;   /* limit to 6 columns */
    }
  else
    { 
      isDown = TRUE ;
      if (!isGifDisplay) 
	xMax = mapGraphWidth ;
      else
	xMax = 80 ;
    }

  for (ii = 1 ; ii < arrayMax(look->segs) ; ++ii)
    { 
      foundCompositeColour = FALSE ;	       
      seg = arrp(look->segs,ii,SEG) ;
      if (showExonSupport)
	{
	  switch (seg->type)
	    {
	    case EXON_SUPPORT: if (isDown) break ; else continue ;
	    case EXON_SUPPORT_UP: if (!isDown) break ; else continue ;
	    default: continue ;
	    }
	}
      else
	{
	  switch (seg->type)
	    {
	    case SPLICED_cDNA: if (isDown) break ; else continue ;
	    case SPLICED_cDNA_UP: if (!isDown) break ; else continue ;
	    default: continue ;
	    }
	}


     if (class(seg->key) == _VMethod)
       continue ;

      y1 = MAP2GRAPH(look->map,seg->x1 - .5) ;
      y2 = MAP2GRAPH(look->map,seg->x2 + .5) ;	/* to cover full base */
      if (y1 < .2 + topMargin) y1 = .2 + topMargin ;
      if (y2 < .2 + topMargin) y2 = .2 + topMargin ;
      if (y2 > mapGraphHeight - .2)	y2 = mapGraphHeight - .2 ;  

      if (seg->x2 < look->min)
	{
	  if (box > 0) box = -1 ;
	  continue ;
	}
      if (seg->x1 > look->max)
	{
	  if (x < xMax && box > 0 && seg->type == (seg - 1)->type &&
	      oldy2 < mapGraphHeight &&
	      cDNAFillData(seg) && cDNAFillData(seg - 1) && 
	      seg->parent == (seg - 1)->parent) /* link to previous */
	    { 
	      if (seg->key == (seg - 1)->key)
		{
		  float mid = (mapGraphHeight + oldy2)/2.0 ;

		  graphColor (MY_MAGENTA) ;  /* YELLOW */  
		  graphLine (x +.4, oldy2, x + 1.5 , mid) ; /* .4 was .8 */
		  graphLine (x +1.5, mid, x + .4, mapGraphHeight) ;
		}
	      else
		{
		  graphColor (isTooLong(seg) ? RED : GREEN) ;  /* YELLOW */
		}
	      graphColor (BLACK) ;
	    }
	  if (box > 0)
	    box = -1 ;
	  continue ;
	}
 
      x1 = seg->data.i >> 14 & 0x3FFF ;
      x2 = seg->data.i & 0x3FFF ;
      is3p = isReadUp = x1 > x2 ? TRUE : FALSE ;
      if (!isDown) isReadUp = !isReadUp ;
      cDNAFillData (seg) ;
      if (seg->parent == (seg - 1)->parent)
	{
	  oldy2 = MAP2GRAPH (look->map, (seg - 1)->x2 - .5) ;
	  if (oldy2 < .2 + topMargin) oldy2 = .2 + topMargin ;

	  if (y1 > oldy2)
	    {
	      ytop = oldy2 + .1 ;
	      yz = ytop ; ix = 0 ;
	      if ((seg-1)->x2 < look->min  && seg->parent != bumpedParent)
		{
		  SEG *seg2 = seg ;
		  int utop, ubottom ;
		  
		  ubottom = seg->x2 ;
		  utop = GRAPH2MAP (look->map, oldy2) ;
		  bumpedParent = seg->parent ;
		  while (seg2->parent == bumpedParent)
		    {
		      if (seg2->x1 < utop) utop = seg2->x1 ;
		      if (seg2->x2 > ubottom) ubottom = seg2->x2 ;
		      seg2++ ;
		    }
		  ytop = MAP2GRAPH (look->map, utop - .5) + .1 ;
		  if (ytop < .2 + topMargin) ytop = .2 + topMargin ;
		  ybottom = MAP2GRAPH (look->map, ubottom - .5) - .1 ;
		  if (ybottom > mapGraphHeight - .2)	ybottom = mapGraphHeight - .2 ;   

		  if (!bump)
		    bump = bumpCreate (mapGraphWidth, 0) ;
		  if (ybottom > ytop + .2)
		    bumpItem (bump,1, ybottom - ytop,&ix,&yz) ;
		  x = *offset + 1.7*ix ;
		  bumpedParent = seg->parent ;
		}
	      if (x < xMax && y1 > oldy2 && seg->key == (seg-1)->key)
		{
		  float mid = (y1 + oldy2) / 2 ;
		  int bb = graphBoxStart () ;

		  graphColor (MY_MAGENTA) ; /* BLACK MAGENTA  ORANGE) ; */
		  graphLine (x +.4, oldy2, x + 1.5 , mid) ; /* .4 was .8 */
		  graphLine (x +1.5, mid, x + .4, y1) ;
		  graphColor (BLACK); 
		  graphLine (x + .1, y1, x + 1.5, y1) ; 
		  
		  graphBoxEnd () ;
		  {
		    int ux1 = seg->x1, ux2 = seg->x2 ;
		    seg->x1 = (seg-1)->x2 + 1 ;
		    seg->x2 = ux1 - 1 ;
		    fMapBoxInfo (look, bb, seg) ;
		    array(look->boxIndex, bb, int) = ii ;
		    graphBoxFreeMenu(bb, fMapcDNAAction, fMapcDNAMenu) ; 
		    estBubbleInfo (bb, seg, look->origin, 1) ;
		    seg->x1 = ux1 ;
		    seg->x2 = ux2 ;
		  }
		}
	      else /* different component */
		{
		  if (x < xMax && oldy2 < y1)
		    {
		      float mid = (y1 + oldy2)/2.0 ;
		      int bx1 = graphBoxStart () ;
		      graphColor (isTooLong(seg) ? RED : GREEN) ;  /* YELLOW */	
		      graphLine (x +.8, oldy2, x + .8, y1) ;
		      graphColor (TRANSPARENT) ;
		      graphLine (x +.2, mid, x +1.4, mid) ;
		      graphColor (BLACK); 
		      graphBoxEnd () ;
		      estBubbleInfo (bx1, seg, look->origin, 2) ;
		    }  
		}	  
	    }
	}
      else
	{

	  if ((seg->type | 0x1) == ((seg - 1)->type | 0x1) && 
	      cDNAFillData(seg - 1) &&
	      seg->parent == (seg - 1)->parent &&
	      (seg-1)->x2 < look->min) /* was clipped by window top */
	    ytop = .2 + topMargin ;
          else
	    ytop = y1 ;

	  yz = ytop ; ix = 0 ;
	  if (seg->parent != bumpedParent)
	    {
	      SEG *seg2 = seg ;
	      int utop, ubottom ;

	      bumpedParent = seg->parent ;
	      ubottom = seg->x2 ;
	      utop = GRAPH2MAP (look->map, ytop) ;
	      while (seg2->parent == bumpedParent)
		{
		  if (seg2->x1 < utop) utop = seg2->x1 ;
		  if (seg2->x2 > ubottom) ubottom = seg2->x2 ;
		  seg2++ ;
		}
	      ytop = MAP2GRAPH (look->map, utop - .5) + .1 ;
	      if (ytop < .2 + topMargin) ytop = .2 + topMargin ;
	      ybottom = MAP2GRAPH (look->map, ubottom - .5) - .1 ;
	      if (ybottom > mapGraphHeight - .2)	ybottom = mapGraphHeight - .2 ;   

	      if (!bump)
		bump = bumpCreate (mapGraphWidth, 0) ;
	      bumpItem (bump,1,(ybottom - ytop + 1)+0.2,&ix,&yz) ;
	      x = *offset + 1.7*ix ;
	    }

	  if (x < xMax &&
	      seg->type == (seg - 1)->type && 
	      cDNAFillData(seg) && cDNAFillData(seg - 1) &&
	      seg->parent == (seg - 1)->parent) /* was clipped by window top */
	    {
	      if (seg->key == (seg - 1)->key)
		graphColor (MY_MAGENTA) ;  /* YELLOW */
	      else
		graphColor (isTooLong(seg) ? RED : GREEN) ;  /* YELLOW */
	      graphLine (x +.8, .2 + topMargin, x + .8, y1) ;
	      graphColor (BLACK) ;
	    }

	  else if (x < xMax &&
		   look->max - look->min < 30000)
	    {
	      KEY motif = 0 ;
	      int xx = 0 ; char *tt = 0 ;
	      OBJ obj = bsCreate (seg->key) ;
	      KEY bestMotif = 0 ;
	      BOOL bold = FALSE ;

	      if (obj)
		{
		  if (bsFindTag (obj, _Real_5prime))
		    {
		      bestMotif = _Real_5prime ; bold = TRUE ;
		    }
		  if (bsGetKey (obj, _Transpliced_to, &motif))
		    do 
		      {
			mark = bsMark(obj, mark);
			xx = 0 ; tt = 0 ;
			if (!strncmp (name(motif), "SL",2))
			  {
			    if (!bestMotif ||
				!strcmp (name(bestMotif), "SL0") ||
				(strcmp (name(motif), "SL0") && 
				 strcmp (name(motif), name(bestMotif)) == -1))
			      bestMotif = motif ;
			    if (strcmp (name(motif), "SL0") && 
				bsGetData (obj, _bsRight, _Int, &xx) &&
				bsGetData (obj, _bsRight, _Text, &tt) &&
				strstr (tt, "Exact"))
			      bold = TRUE ;
			  }
			bsGoto(obj, mark);
		      } while (bsGetKey (obj, _bsDown, &motif)) ;
		  else if (bsFindTag (obj, _Real_3prime))
		    {
		      bestMotif = _Real_3prime ; bold = TRUE ;
		    }

		  if (bestMotif) 
		    motif = bestMotif ;
		  if (motif)
		    {
		      float old = 0 ; int oldFormat = 0 ;

		      if (!bold) old = graphTextHeight (0.75) ;
		      else oldFormat = graphTextFormat(BOLD) ;

		      if (!strncasecmp(name(motif),"SL",2))
			{
			  int bx1  = graphBoxStart () ;
			  graphText(name(motif) + 2,x + .2, y1 - .95) ;
			  graphBoxEnd () ;
			  if (!strcasecmp(name(motif),"SL0"))
			    graphBubbleInfo (bx1, 0, 0, "Capped clone") ;
			  else
			    graphBubbleInfo (bx1, 0, 0, messprintf ("Transpliced leader %s", name(motif))) ;
			}
		      else if (motif == _Real_5prime)
			{
			  if (isDown)
			    {
			      float y11 = MAP2GRAPH (look->map, seg->x1 - 30) ;
			      graphLine (look->mrnaOffset, y11, x + 5, y11) ;
			    }
			  graphLine (x + .1, y1 - .1, x + 1.5, y1 - .1) ; 
			  graphText("=", x + .3, y1 - .95) ;
			}
		      else if (motif == _Real_3prime)
			{
			  if (isDown)
			    {
			      float y22 = MAP2GRAPH (look->map, seg->x2 + 30) ;
			      graphLine (look->mrnaOffset, y22, x + 5, y22) ;
			    }
			  graphLine (x + .1, y2 + .1, x + 1.5, y2 + .1) ; 
			}
		      else if (strcmp(name(motif),"gccgtgctc") &&
			       strncmp(name(motif),"gagaga",6)) ; /* rien du tout */
		      else if (*name(motif) == 'a')
			graphText("a", x + .2, y1 - .95) ;
		      else
			graphText("X", x + .2, y1 - .95) ;

		      if (!bold) graphTextHeight (old) ;
		      else  graphTextFormat (oldFormat) ;
		    }
		  bsDestroy (obj) ;
		}
		  
	    }
	}

      oldy2 = y2 ;

      if (seg->parent != bumpedParent)
	{
	  SEG *seg2 = seg ;
	  int utop, ubottom ;

	  ubottom = seg->x2 ;
	  utop = GRAPH2MAP (look->map, ytop) ;
	  bumpedParent = seg->parent ;
	  while (seg2->parent == bumpedParent)
	    {
	      if (seg2->x1 < utop) utop = seg2->x1 ;
	      if (seg2->x2 > ubottom) ubottom = seg2->x2 ;
	      seg2++ ;
	    }

	  ytop = MAP2GRAPH (look->map, utop - .5) + .1 ;
	  if (ytop < .2 + topMargin) ytop = .2 + topMargin ;
	  ybottom = MAP2GRAPH (look->map, ubottom - .5) - .1 ;
	  if (ybottom > mapGraphHeight - .2)	ybottom = mapGraphHeight - .2 ;   
	  
	  yz = ytop ; ix = 0 ;
	  if (!bump)
	    bump = bumpCreate (mapGraphWidth, 0) ;
	  bumpItem (bump,1,(ybottom - ytop + 1)+0.2,&ix,&yz) ;
	  x = *offset + 1.7*ix ;
	}

      if (x >= xMax) continue ;
      if (0) printf ("est:%s ix=%d x1=%d x2=%d", name(seg->key), ix, seg->x1, seg->x2) ;

      if (
	  (!strncmp (name(seg->key), "XH_",3) || !strncmp (name(seg->key), "XG_",3))  &&
	  keyFindTag (seg->key, _Is_read) && keyFindTag (seg->key, _Composite))
	{
	  OBJ obj = bsCreate(seg->key) ;
	  Array aa =  arrayCreate (24, BSunit) ;
	  
	  if (bsGetArray (obj, _Composite, aa, 3) &&
	      arrayMax (aa) >= 3)
	    {
	      int r, rMax = arrayMax (aa) ;
	      for (r = 0 ; r < rMax ; r+=3)
		{
		  BSunit *u = arrp (aa, r, BSunit) ;
		  int x1 = u[0].i ;
		  int x2 = u[1].i ;
		  int cover = u[2].i ; 
		  float z = log(1.0+cover) ;
		  float z1 = log(8000.0)/8 ;
		  float yy1, yy2 ;
		  int iCol = .5 + z/z1 ;
		  if (iCol > 8) 
		    iCol = 8 ;
		  if (iCol < 0)
		    iCol = 0 ;
		  bgCol = GREEN1 + iCol - 1 ;

		  if (isDown)
		    {
		      yy1 = MAP2GRAPH(look->map,seg->x1 + x1 - 1 - .5) ;
		      yy2 = MAP2GRAPH(look->map,seg->x1 + x2 - 1 + .5) ;
		    }
		  else
		    {
		      yy1 = MAP2GRAPH(look->map,seg->x2 - x2 + 1 - .5) ;
		      yy2 = MAP2GRAPH(look->map,seg->x2 - x1 + 1 + .5) ;
		    }
		  if (yy1 < .2 + topMargin) yy1 = .2 + topMargin ;
		  if (yy2 > mapGraphHeight - .2)  yy2 = mapGraphHeight - .2 ;  
		  if (yy1 < yy2 - .01)
		    {
		      int box2 = graphBoxStart() ; 
		      int oldColour = graphColor (bgCol) ;
		      graphRectangle (x+.1, yy1, x+1.5, yy2) ;
		      graphBoxEnd () ;
		      graphBoxDraw (box2, WHITE, bgCol) ;
		      graphColor (oldColour) ;
		    }
		  bgCol = TRANSPARENT ;
		  foundCompositeColour = TRUE ;		}
	    }
	  arrayDestroy (aa) ;
	  bsDestroy (obj) ;
	}

      /* draw the beginning of the est line as arrow circle or line*/
      if (seg->key != (seg-1)->key)
	{
	  if (is3p)
	    {
	      if (!isReadUp)
		{
		  if (seg->data.i & 0x40000000)
		    {
		      int bx1 = graphBoxStart () ;
		      graphCircle (x + .9, y1, .6) ;
		      graphBoxEnd () ;
		      estBubbleInfo (bx1, seg, look->origin, 3) ;
		    }
		  else
		    {
		      if (seg->data.i & 0x10000000) 
			{ 
			  float yy = MAP2GRAPH(look->map, seg->x1 + .5 - .5) ; 
			  int bx1 = graphBoxStart () ;
			  
			  graphColor (RED) ;
			  if (yy - y1 < .20) yy = y1 + .20 ;
			  graphFillRectangle (x + .1, y1, x + 1.5, yy) ; 
			  graphColor (TRANSPARENT) ;
			  graphLine (x + .1, yy - 1, x + 1.5, yy) ; 
			  graphColor (BLACK); 
			  graphBoxEnd () ;
			  estBubbleInfo (bx1, seg, look->origin, 4) ;
			}
		      else
			graphLine (x + .1, y1, x + 1.5, y1) ;
		    }
		}
	      else  /*  (is3p && isReadUp) */
		{ 
		  if (seg->data.i & 0x20000000) graphColor (RED) ;
		  graphLine (x + .1 , y1 + .6, x + .8, y1) ;
		  graphLine (x + 1.5, y1 + .6, x + .8, y1) ; 
		  if (seg->data.i & 0x20000000) graphColor (BLACK) ;
		}
	    }
	  else  /* !is3p */
	    {
	      if (isReadUp)
		{
		  if (seg->data.i & 0x40000000)
		    {  
		      int bx1 = graphBoxStart () ;
		      graphCircle (x + .9, y1, .6) ;  
		      graphBoxEnd () ;
		      estBubbleInfo (bx1, seg, look->origin, 3) ;
		    }
		  else
		    { 
		      if (seg->data.i & 0x20000000) 
			graphColor (RED) ;
		      graphLine (x + .1 , y1 + .6, x + .8, y1) ;
		      graphLine (x + 1.5, y1 + .6, x + .8, y1) ; 
		      if (seg->data.i & 0x20000000) graphColor (BLACK) ;
		    }
		}
	      else /* (!is3p && !isReadUp) */
		{  
		  if (seg->data.i & 0x10000000) 
		    { 
		      float yy = MAP2GRAPH(look->map, seg->x1 + .5 - .5) ; 
		       int bx1 = graphBoxStart () ;
		       
		       graphColor (RED) ; 
		      if (yy - y1 < .20) yy = y1 + .20 ;
		      graphFillRectangle (x + .1, y1, x + 1.5, yy) ;  
		      graphColor (TRANSPARENT) ;
		      graphLine (x + .1, yy - 1, x + 1.5, yy) ; 
		      graphColor (BLACK); 
		      graphBoxEnd () ;
		      estBubbleInfo (bx1, seg, look->origin, 4) ;
		    }
		  else
		    graphLine (x + .1, y1, x + 1.5, y1) ;
		}
	    }
	}
      else
	graphLine (x + .1, y1, x + 1.5, y1) ;

      graphColor (BLACK) ;

      box = graphBoxStart() ; 
      graphColor (TRANSPARENT) ;
      graphLine (x + .1, (y1 + y2)/2.0 , x + 1.5,  (y1 + y2)/2.0) ; 
      graphColor (BLACK) ;

      if (!isGifDisplay) oldWidth = graphLinewidth (isAmbiguous(seg) ? .5 : .3) ;
      if (isAmbiguous(seg))
	{
	  graphColor (BLUE) ;
	  graphLine (x + .9, y1, x + .9, y2) ;
	  graphColor (BLACK) ;
	}
      else if (keyFindTag (seg->parent, _Is_AM))
	{
	  graphColor (RED) ;
	  graphLine (x + .9, y1, x + .9, y2) ;
	  graphColor (BLACK) ;
	}
      else
	graphLine (x + .9, y1, x + .9, y2) ;

      if (!isGifDisplay) graphLinewidth (oldWidth) ; 

      if (look->dna &&
	  fmapView (look->view, _Fmap_cDNA_Decorate) &&
	  mapColSetByName ("SPLICED_cDNA_DECORATE", -1) )
	fMapSplicedcDNADecorate (look, x + .8, seg, isReadUp, x1, x2) ;
      graphBoxEnd () ;

      bgCol = 0 ;
      if (bgCol) {} ;
#ifdef ACEMBLY
      bgCol =  fMapQueryColor (seg->key) ;

#endif
      if (foundCompositeColour)
	bgCol = TRANSPARENT ;
      else if (keyFindTag (seg->key, _Colour))
	{
	  OBJ obj = bsCreate(seg->key) ;
	  bsGetKeyTags (obj,_Colour, &bgCol) ;
	  bsDestroy (obj) ;
	  bgCol = WHITE + bgCol - _WHITE ;
	}
      else if (_aForward  &&
	  (
	   (keyFindTag (seg->key, _Forward) && keyFindTag (seg->key, _aReverse)) ||
	   (keyFindTag (seg->key, _aForward) && keyFindTag (seg->key, _Reverse)) 
	   )
	  )
	bgCol = PALEYELLOW ;
      else
	bgCol = TRANSPARENT ;
      graphBoxDraw (box, BLACK, bgCol) ;

      fMapBoxInfo (look, box, seg) ;

      array(look->boxIndex, box, int) = ii ;
      graphBoxFreeMenu(box, fMapcDNAAction, fMapcDNAMenu) ; 
      estBubbleInfo (box, seg, look->origin, 0) ;

      /* draw the end of the est line as arrow circle or line*/
      if (seg->key != (seg+1)->key)
	{
	  if (is3p)
	    {
	      if (isReadUp)
		{
		  if (seg->data.i & 0x40000000)
		    {
		      /* int bx1 = 0 ; obsured the mouse click  graphBoxStart () ; */
		      graphCircle (x + .9, y2, .6) ; 
		      /*
			graphBoxEnd () ;
			estBubbleInfo (bx1, seg, look->origin, 3) ;
		      */
		      y2 += .6 ;
		    }
		  else
		    {
		      if (seg->data.i & 0x10000000) 
			{ 
			  float yy = MAP2GRAPH(look->map, seg->x2 + .5 - .5) ; 
			  int bx1 = graphBoxStart () ;

			  graphColor (RED) ; 
			  if (yy - y2 < .20) yy = y2 + .20 ;
			  graphFillRectangle (x + .1, y2, x + 1.5, yy) ; 
			  graphColor (TRANSPARENT) ;
			  graphLine (x + .1, yy - 1, x + 1.5, yy) ; 
			  graphColor (BLACK); 
			  graphBoxEnd () ;
			  estBubbleInfo (bx1, seg, look->origin, 4) ;
			}
		      else
			graphLine (x + .1, y2, x + 1.5, y2) ;
		    }
		}
	      else /* (is3p && !isReadUp) */
		{  
		  if (seg->data.i & 0x20000000)
		    graphColor (RED) ;
		  graphLine (x + .1 , y2 - .6, x + .8, y2) ;
		  graphLine (x + 1.5, y2 - .6, x + .8, y2) ;  
		  if (seg->data.i & 0x20000000) graphColor (BLACK) ;
		}
	    }
	  else /* !is3p */
	    {
	      if (!isReadUp)
		{
		  if (seg->data.i & 0x40000000)
		    {
		      /* int bx1 = 0 ; obscured the mouse click */
		      
		      graphCircle (x + .9, y2, .6) ; 
		      /*
			graphBoxEnd () ;
			estBubbleInfo (bx1, seg, look->origin, 3) ;
		      */
		      y2 += .6 ;
		      if (0) /* test TRANSPARENT */
			{ 
			  int bx = graphBoxStart () ;
			  graphRectangle (x - 4, y2 - 4 , x + 4, y2 + 4) ;
			  graphBoxEnd() ;
			  graphBoxDraw (bx, RED, TRANSPARENT) ;
			}
		    }
		  else
		    { 
		      if (seg->data.i & 0x20000000) 
			graphColor (RED) ;
		      graphLine (x + .1 , y2 - .6, x + .8, y2) ;
		      graphLine (x + 1.5, y2 - .6, x + .8, y2) ; 
		      if (seg->data.i & 0x20000000) graphColor (BLACK) ;
		    }
		}
	      else /* (!is3p && isReadUp) ; */
		{
		  if (seg->data.i & 0x10000000) 
		    { 
		      float yy = MAP2GRAPH(look->map, seg->x2 + .5 - .5) ; 
		      int bx1 = graphBoxStart () ;
			  
		      graphColor (RED) ; 
		      if (yy - y2 < .20) yy = y2 + .20 ;
		      graphFillRectangle (x + .1, y2, x + 1.5, yy) ; 
		      graphColor (TRANSPARENT) ;
		      graphLine (x + .1, yy - 1, x + 1.5, yy) ; 
		      graphColor (BLACK); 
		      graphBoxEnd () ;
		      estBubbleInfo (bx1, seg, look->origin, 4) ;
		    }
		  else
		    graphLine (x + .1, y2, x + 1.5, y2) ;
		}
	    }
	}
      else
	 graphLine (x + .1, y2, x + 1.5, y2) ;

      if (seg->key != (seg+1)->key && /* isReadUp && */
	  look->max - look->min < 200000) /* was 30000 for worm */
	{ 
	  KEY clone = keyGetKey (seg->key, _cDNA_clone) ;
	  KEY clone2 = keyGetKey ((seg+1)->key, _cDNA_clone) ;
	  KEY lib = keyGetKey (clone, _Library) ;
	  int bb = graphBoxStart () ;

	  if (clone != clone2 && lib && bIndexFind (lib, _Symbol))
	    {
	      OBJ Library = bsCreate (lib) ;
	      char *symbol ;
	      
	      if (Library && bsGetData (Library, _Symbol, _Text, &symbol))
		{
		  float old = graphTextHeight (0.75) ;
		  graphText(symbol, x +.45, y2 + .3) ; 
		  y2 += .8 ; 
		  graphTextHeight (old) ;
		}
	      bsDestroy (Library) ;

	    } 
	  if (clone != clone2  && keyFindTag (clone, _Resequence))
	    { graphText("*",x + .45, y2 + .3) ; y2 += .5 ; }
	  if (clone != clone2 &&
	      keyFindTag (clone, _Anomalous_clone))
	    {
	      BOOL lengthAnomaly =FALSE ;
	      int bb2 = graphBoxStart() ;
	      float pbSize = .3 ;   /* little problem */

	      if (keyFindTag (clone, _Length_anomaly))
		{ lengthAnomaly = TRUE  ; }
	      if (keyFindTag (clone, _Suspected_internal_deletion))
		{ lengthAnomaly = FALSE ; pbSize = .9 ; }                 /* small problem  but ignore length info */
	      if (keyFindTag (clone, _Internal_priming))
		{ lengthAnomaly = FALSE ; }                 /* small problem  but ignore length info */
	      if (keyFindTag (clone, _Mosaic))
		{ lengthAnomaly = FALSE ; pbSize = .9 ; } /* big problem */
	      
	      if (lengthAnomaly) 
		{
		  int oldFormat = graphTextFormat(BOLD) ;
		  graphText ("?", x + .2, y2 + .3) ; y2 += .8 ;
		  graphTextFormat(oldFormat) ;
		  graphBoxEnd () ;
		  graphBoxDraw (bb2, RED, TRANSPARENT) ;
		}
	      else
		{	      
		  graphPointsize (pbSize) ;
		  graphPoint (x + .9, y2 + .3) ; y2 += .8 ;
		  graphBoxEnd () ;
		  graphBoxDraw (bb2, RED, RED) ;
		}
	    }
	  graphBoxEnd () ;
	  graphBoxDraw (bb, BLACK, TRANSPARENT) ;
	}


    }
  
  if (bump) 
    {
      *offset += 1.7*bumpMax (bump) ;
      bumpDestroy (bump) ; 
    }
  if (*offset > xMax) *offset = xMax ;

  bsMarkFree(mark);
  fMapShowSequencingErrors (look, offset) ;
  graphColor (BLACK) ;

} /* fMapcDNAShowSplicedcDNA */

void fMapcDNAShowSplicedcDNA (LOOK look, float *offset)
{
  return fMapcDNAShowSplicedcDNAorExon (look, offset, FALSE) ;
}
void fMapcDNAShowExonSupport (LOOK look, float *offset)
{
  return fMapcDNAShowSplicedcDNAorExon (look, offset, TRUE) ;
}

/*********************************************************************/

static void fMapSplicedcDNADecorate (LOOK look, float x, SEG *seg, BOOL isUp, int x1, int x2)
{
  A_ERR *ep ;
  Array a = 0 ;
  Array dna = look->dna, dnaR ;
  int oldColor, color = 0, i, a1 = seg->x1, a2 = seg->x2 ;
  float y, yy, y1, y2 ;
  int isAm = -1 ;

  if (isUp && !look->dnaR)
    { look->dnaR = dnaCopy (look->dna) ;
      reverseComplement (look->dnaR) ;
    } 
  dnaR = look->dnaR ;

  a = findSpliceErrors (look->seqKey, seg->key, dna, dnaR, a1, a2, x1, x2, isUp) ;
  if (!a || !arrayMax (a))
    goto done ;
  
  i = arrayMax(a) ;
  ep = arrp(a, 0, A_ERR)  - 1 ;
  y1 = y2 = 0 ; oldColor = -1 ;
  while (ep++, i--)
    {
      if (ep->iLong > a2)
	break ;
      if (ep->iLong < a1)
	continue ;
      switch(ep->type)
	{
	case ERREUR: 
	  color = RED ; 
	  break ;
	case INSERTION: /* _AVANT: */
	case INSERTION_DOUBLE:  /*_AVANT:  */
	case TROU:
	case TROU_DOUBLE:
	  color = BLUE ; 
	  break ;
	case AMBIGUE: 
	  if (isAm == -1) /* it always pays to be lazy ! */
	    isAm = keyFindTag (seg->key, _Is_AM) ;
	  color = isAm ? DARKGREEN : -1  ; /* show the N_ just for the AM sequences */
	  break ;
	default:
	  break ;
	}
      y = ep->iLong  ;
      yy = ep->iLong  ;
	
      y = MAP2GRAPH(look->map, y - .5) ;
      yy = MAP2GRAPH(look->map, yy + .5) ;
      if (yy - y < .10) yy = y + .10 ;
      if (y > topMargin + .2 && y < mapGraphHeight - .2)
	{ 
	  if (y < y2 - .1)
	    { y2 = yy ; /* coalesce the 2 boxes */
	      if (color != oldColor)
		oldColor = DARKRED ;
	    }
	  else
	    { if (oldColor != -1)
		{ 
		  graphColor(oldColor) ;
		  if (y2 - y1 < .12)
		    graphLine(x - .35, y1, x + .35, y1) ;
		  else
		    graphFillRectangle(x - .35, y1, x + .35, y2) ;
		  graphColor (BLACK) ;
		}
	      /* register */
	      y1 = y ; y2 = yy ; oldColor = color ;
	    }
	}
    }
  
  if (oldColor != -1)
    { graphColor(oldColor) ;
      graphFillRectangle(x - .35, y1, x + .35, y2) ;
      graphColor (BLACK) ;
    }
  graphColor (BLACK) ;

 done:
  if (arrayExists(a))
    arrayDestroy (a) ; 
  return ;
}  /* fMapSplicedcDNADecorate */

/***************************************************************************************/

static Array findSpliceErrors (KEY link, KEY key, Array dnaD, Array dnaR, int a1, int a2, 
				 int x1, int x2, BOOL upSequence)
{
  Array a = 0, mydna = 0 ;
  OBJ obj ;
  int x3, b1, b2, b3, u ;
  KEY mydnakey = 0 ;

  if ((obj = bsCreate(key)))
    { if (bsGetKey (obj, _DNA, &mydnakey))
	{ 
	  if ((mydna = dnaGet(mydnakey)))
	    {
	      /* the window may have only part of the dna */
	      if (a1 < 0) 
		{
		  if (upSequence) 
		    { x2 -= a1 ; a1 -= a1 ; }
		  else
		    { x1 -= a1 ; a1 -= a1 ; }
		}
	      if ((u = a2 - (long)arrayMax(dnaD)) > 0)
		{
		  if (upSequence) 
		    { x1 -= u ; a2 -= u ; }
		  else
		    { x2 -= u ; a2 -= u ; }
		}
	      if (upSequence) 
		{  
		  b1 = a2 ;
		  b2 = a1 ;
		}
	      else
		{ 
		  b1 = a1 ; b2 = a2 ;
		}
	      if (x1 > x2) 
		{ x3 = x1 ; x1 = x2 ; x2 = x3 ; } 
	      b3 = b2 ; x3 = x2 ; x1-- ;
	      if (a1 < a2 && x1 < x2)
		a = spliceCptErreur (dnaD, dnaR, mydna, upSequence, 
				     b1, b2, b3,
				     &x1, &x2, &x3) ;
	    }
	  arrayDestroy (mydna) ;
	}
      bsDestroy (obj) ;
    }
  if (a && !arrayMax(a))
    arrayDestroy (a) ;
  return a ;
}

/***************************************************************************************/
/***************************************************************************************/
static void cleanErr (Array err, Array dnaLong, Array dnaShort, BOOL isUp)
{ int i, j, n = arrayMax(err), nn, start, stop ;
  A_ERR *eq, *ep, *epMax ;
  char *cl, *cs, cc ;
  int sens = isUp ? -1 : 1 ;
  ALIGN_ERROR_TYPE type ;

  if (n < 2) return ;
  if (isUp)
    { ep = arrp (err, 0, A_ERR) - 1 ;
      n-- ;
      while (ep++, n--)
	if ((ep->type == INSERTION) && ((ep+1)->type == ERREUR) &&
	    (ep->iShort == (ep+1)->iShort + 1))
	  { type = ep->type ; ep->type = (ep+1)->type ; (ep+1)->type = type ;
	    i = ep->iLong ; ep->iLong = (ep+1)->iLong ; (ep+1)->iLong = i + 1;
	  }
    }

  for (n = 0 ; n < arrayMax(err) ; n++)
    { 
      /* we cannot ep++ because err may get extended inside this loop, mieg 2006/07/29 */
      ep = arrp(err, n, A_ERR) ;
      nn = arrayMax (err) - 1 ;
      switch (ep->type)
	{
	case AMBIGUE:
	  if (isUp)
	    { 
	      if (n > 0 &&
		  (ep - 1)->type == INSERTION &&
		  (ep - 1)->iShort < arrayMax(dnaShort) &&
		  (ep - 1)->iLong >= 0 &&
		  (ep - 1)->iLong < arrayMax(dnaLong))
		{ cl = arrp(dnaLong, (ep - 1)->iLong, char) ; 
		  cs = arrp(dnaShort, (ep - 1)->iShort, char) ;
		  i = j = 0 ;
		  cc = *cl ; while (*cl++ == cc) i++ ;
		  if (sens == -1)
		    { cc = complementBase(cc) ;
		      while (*cs-- == cc) j++ ;
		    }
		  else
		    while (*cs++ == cc) j++ ; 
		  if (i >= j && (ep - 1)->iShort == ep->iShort + j)
		    { 
		      eq = ep - 2 ;
		      i = nn - n + 1 ;
		      while (eq++, i-- > 0) *eq = *(eq + 1) ;
		      arrayMax(err) -= 1 ;
		      ep-- ;
		      ep->type = INSERTION ;
		      ep->iLong++ ;
		    }
		}
	      if (n < nn &&
		  (ep + 1)->type == INSERTION &&
		  (ep + 1)->iShort < arrayMax(dnaShort) &&
		  (ep + 1)->iLong >= 0 &&
		  (ep + 1)->iLong < arrayMax(dnaLong))
		{ cl = arrp(dnaLong, (ep + 1)->iLong, char) ; 
		  cs = arrp(dnaShort, (ep + 1)->iShort, char) ;
		  i  = j = 0 ;
		  cc = *cl ; while (*cl++ == cc) i++ ;
		  if (sens == -1)
		    { cc = complementBase(cc) ;
		      while (*cs-- == cc) j++ ;
		    }
		  else
		    while (*cs++ == cc) j++ ; 
		  if (i >= j && (ep + 1)->iShort == ep->iShort + j)
		    { 
		      eq = ep - 1 ;
		      i = nn - n  ;
		      while (eq++, i-- > 0) *eq = *(eq + 1) ;
		      arrayMax(err) -= 1 ;
		      ep-- ;
		      ep->type = INSERTION ;
		      ep->iLong++ ;
		    }
		}
	    }
	  else
	    { if (n < nn &&
		  (ep + 1)->type == INSERTION &&
		  ep->iShort < arrayMax(dnaShort) &&
		  ep->iLong >= 0 &&
		  ep->iLong < arrayMax(dnaLong))
/*		  (ep + 1)->type == INSERTION_DOUBLE) plus difficile */
		{ cl = arrp(dnaLong, ep->iLong, char) ; 
		  cs = arrp(dnaShort, ep->iShort, char) + sens ;
		  i = j = 0 ;
		  cc = *cl ; while (*cl++ == cc) i++ ;
		  if (sens == -1)
		    { cc = complementBase(cc) ;
		      while (*cs-- == cc) j++ ;
		    }
		  else
		    while (*cs++ == cc) j++ ; 
		  if (i == j && (ep + 1)->iShort == ep->iShort + i * sens)
		    { ep->type = INSERTION ;
		      eq = ep ;
		      i = nn - n - 1 ;
		      while (eq++, i-- > 0) *eq = *(eq + 1) ;
		      arrayMax(err) -= 1 ;
		    }
		}
	    }
	  break ;
	case INSERTION: case INSERTION_DOUBLE:
	  start = ep->iShort ;
	  epMax = arrp (err, nn, A_ERR) + 1 ;
	  if (n <= nn && !isUp && ep->iShort >= 0 &&
	      ep->iShort < arrayMax(dnaShort))
	    { cs = arrp(dnaShort, ep->iShort, char) ;
	      cc = *cs ; i = arrayMax(dnaShort) - ep->iShort ;
	      while (i-- > 0 && *(++cs) == cc)
		{ ep->iShort++ ;
		  ep->iLong++ ;
		}
	      if (ep->type == INSERTION_DOUBLE)
		ep->iLong-- ;
	    }
	  else if (n <= nn && isUp && ep->iShort >= 0 &&
		   ep->iShort < arrayMax(dnaShort))
	    { cs = arrp(dnaShort, ep->iShort, char) ;
	      cc = *cs ; i = ep->iShort ;
	      while (i-- > 0 && *(--cs) == cc && ep->iShort > 0)
		{ ep->iShort-- ;
		  ep->iLong++ ;
		}
	    }
	  stop = ep->iShort ;
	  eq = ep ;
	  if (isUp)
	    while (++eq < epMax && eq->iShort < start && 
		   eq->iShort >= stop)
	      eq->iShort++ ;
	  else
	    while (++eq < epMax && eq->iShort > start && 
		   eq->iShort <= stop)
	      eq->iShort-- ;
	  /* condense insertion-error, nov 1 , 2003 */
	  if (n < nn &&
	      ep->type == INSERTION && n &&
	      (ep+1)->type == ERREUR &&
	      ep->iLong == (ep+1)->iLong &&
	      ep->iShort - 1 == (ep+1)->iShort)
	    { 
	      i = nn - n - 2 ;
	      eq = ep ;
	      while (++eq, i-- > 0)
		*eq = *(eq + 1) ;
	      n-- ;
	      arrayMax (err) -= 1 ;
	    }
	  if (ep->type == INSERTION_DOUBLE)
	    { 
	      i = nn - n + 1 ;
	      eq = arrayp (err, arrayMax (err), A_ERR) + 1 ;
	      ep = arrp (err, n, A_ERR) ;
	      while (eq--, i-- > 0)
		*eq = *(eq - 1) ;
	      ep->type = INSERTION ;
	      ep->baseShort = W_ ;
	      if (!isUp)
		ep->iShort-- ;
	      ep++ ;
	      ep->type = INSERTION ;
	      if (isUp)
		ep->iShort++ ;
	      ep->baseShort = W_ ;
	      ep++ ; if (n) n-- ;
	    }
	  break ;
	case TROU: case TROU_DOUBLE:
	  if (ep->iLong >= 0 && ep->iLong < arrayMax(dnaLong))
	    { cs = arrp(dnaLong, ep->iLong, char) ;
	      cc = *cs ; i = arrayMax(dnaLong) - ep->iLong ;
	      while (i-- > 0 && *(++cs) == cc)
		{ ep->iShort += sens ;
		  ep->iLong++ ;
		}
	      if (ep->type == TROU_DOUBLE && !isUp)
		ep->iShort -- ;
	    }
	  break ;
	case ERREUR:
	  break ;
	default:
	  break ;
	}
    }
}

/*********************************************************/

static void fuseErr (Array err1, Array err2)
{ int i, n1 = arrayMax(err1), n2 = arrayMax(err2) ;
  for (i = 0 ; i < n2 ; i++)
    array(err1, n1 + i, A_ERR) = arr(err2, i, A_ERR) ;
}

/*********************************************************/

static Array spliceCptErreur (Array dnaD, Array dnaR, Array dnaCourt, BOOL isUp,
			 int x1, int x2, int x3, 
			 int *clipTopp, int *cEndp, int *cExp)
{ int amax, ctop, c2, c3, i, nn ;
  Array err1 = 0, err2 = 0, dnaLong ;
  A_ERR *ep, *eq, ee ;

  amax = arrayMax (dnaD) ;

  if (isUp)
    { x1 = amax - x1 - 1 ;
      x2 = amax - x2 - 1 ;
      x3 = amax - x3 - 1 ;
      dnaLong = dnaR ;
    }
  else
    dnaLong = dnaD ;

  ctop = *clipTopp ;
  x2-- ;  c2 = *cEndp -1 ;
  err1 = /* makeErr2 (dnaLong, dnaCourt, &x1, &x2, &ctop, &c2) ;*/
    aceDnaTrackErrors (dnaCourt, ctop, &c2, dnaLong, x1, &x2, &nn, err1, 8, 0, FALSE, 0, FALSE) ;
  *cEndp = c2 = c2 + 1 ; x2++ ;  /* may have moved a bit */      
  *clipTopp = ctop ;
  c3 = *cExp - 1 ;  
  if (x2 < x3 && c2 < c3)
    { err2 =  /* makeErr2 (dnaLong, dnaCourt, &x2, &x3, &c2, &c3) ; */
	aceDnaTrackErrors (dnaCourt, c2, &c3, dnaLong, x2, &x3, &nn, err2, 8, 0, FALSE, 0, FALSE) ;
      x3++ ;
      *cExp = c3 + 1 ;
      fuseErr (err1, err2) ;  /* accumulates in err1 */
      arrayDestroy (err2) ;
    }
  if (isUp)
    { x1 = amax - x1 - 1 ;
      x2 = amax - x2 - 1 ;
      x3 = amax - x3 - 1 ;
      i = arrayMax(err1) ;
      if (i)
	{ ep = arrp (err1, 0, A_ERR) - 1 ;
	  while (ep++, i--)
	    { 
	      switch (ep->type)
		{
		case INSERTION:
		case INSERTION_DOUBLE:
		  ep->iLong-- ;
		  break ;
		case TROU:
		case TROU_DOUBLE:
		  ep->iShort-- ;
		  break ;
		case AMBIGUE:
		case ERREUR:
		  break ;
		default:
		  break ;
		}
	      ep->iLong = amax - ep->iLong - 1 ;
	    }
	  i = arrayMax(err1) ;
	  ep = arrp (err1, 0, A_ERR) ;
	  eq = arrp (err1, i - 1, A_ERR) ;
	  while (ep < eq)
	    { ee = *ep ; *ep++ = *eq ; *eq-- = ee ; }
	}
    }
  cleanErr(err1, dnaD, dnaCourt, isUp) ;
  return err1 ;
}

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

void fMapcDNAReportGenes (char *buffer, SEG *seg1, int maxBuf)
{
  char *cp = strnew(buffer,0) ;
  KEY newName = keyGetKey (seg1->key, _NewName) ;

  if (newName && strcmp(name(newName), name(seg1->key)))
    strncpy (buffer, messprintf ("Gene %s(%s), %s", 
				 name(seg1->key), name(newName),
				 cp + strlen(name(seg1->key)) + 1),
	     maxBuf) ;
  else
    strncpy (buffer, messprintf ("Gene %s, %s", 
				 name(seg1->key),
				 cp),
	     maxBuf) ;
  messfree (cp) ;
}

/*********************************************************************/

void fMapcDNAReportRNAi (char *buffer, SEG *seg1, int maxBuf)
{
  char *cp = strnew(buffer,0) ;
  KEY ph = keyGetKey (seg1->key, _Phenotype) ;

  strncpy (buffer, messprintf ("RNAi %s, %s %s", 
			       name(seg1->key),
			       cp, ph ? name(ph) : ""),
	   maxBuf) ;
  messfree (cp) ;
}

/*********************************************************************/

void fMapcDNAReportProbe (char *buffer, SEG *seg1, int maxBuf)
{
  strncpy (buffer, messprintf ("Probe %s, %s"
			       , name(seg1->key)
			       , seg1->data.i & 0x4000 ? "exact" : "approximate"
			       )
	   , maxBuf) ;
}

/*********************************************************************/

void fMapcDNAReportOST (char *buffer, SEG *seg1, int maxBuf)
{
  char *cq, *cp = strnew(buffer,0), *ca, *cb, *ct ;
  OBJ OST = 0 ;

  cq = cp + strlen(cp) - 1 ;
  while (cq > cp && *cq == ' ') *cq-- = 0 ;

  if (keyFindTag (seg1->key, str2tag ("Amplified")))
    ca = "Amplified" ;
  else if (keyFindTag (seg1->key, str2tag ("Not_amplified")))
    ca = "Not amplified" ;
  else
    ca = "" ;

  if (keyFindTag (seg1->key, str2tag ("Ambiguous")))
    cb = "Non unique alignement" ;
  else if (keyFindTag (seg1->key, str2tag ("Amplified_not_ambiguous")))
    cb = "unique alignement" ;
  else
    cb = "" ;

  ct = "" ;
  if (keyFindTag (seg1->key, _Target  ) &&
      (OST = bsCreate (seg1->key)))
    bsGetData (OST, _Target, _Text, &ct) ;

  strncpy (buffer
	   , messprintf ("OST %s, %s %s %s%s" 
			 , cp
			 , ca
			 , cb
			 , *ct ? "target = " : ""
			 , ct
			 ) 
	   ,maxBuf
	   ) ;
  messfree (cp) ;
  bsDestroy (OST) ;
}

/*********************************************************************/

void fMapcDNAReportMrna (char *buffer, SEG *seg1, int maxBuf)
{
  char *cp = strnew(buffer,0) ;
  if (seg1->data.i && !strstr (name(seg1->key), "Exon"))
    strncpy (buffer, messprintf ("Gene %s, %s  %s", 
				 name(seg1->parent), 
				 cp, name(seg1->data.i)),
	     maxBuf) ;
  else
    strncpy (buffer, messprintf ("Gene %s, %s", 
				 name(seg1->parent),
				 cp),
	     maxBuf) ;
  messfree (cp) ;
}

/*********************************************************************/

void fMapcDNAReportMProduct (char *buffer, SEG *seg1, int maxBuf)
{
  char *cp = strnew(buffer,0) ;
  if (seg1->data.i)
    strncpy (buffer, messprintf ("Product %s, %s  %s", 
				 name(seg1->parent), 
				 cp, name(seg1->data.i)),
	     maxBuf) ;
  else
    strncpy (buffer, messprintf ("Product %s, %s", 
				 name(seg1->parent),
				 cp),
	     maxBuf) ;
  messfree (cp) ;
}

/*********************************************************************/

void fMapcDNAReportTranscript (char *buffer, SEG *seg1, int maxBuf)
{
  char *cp = strnew(buffer,0) ;
  if (seg1->data.i)
    strncpy (buffer, messprintf ("Transcript %s, %s  %s", 
				 name(seg1->parent), 
				 cp, name(seg1->data.i)),
	     maxBuf) ;
  else
    strncpy (buffer, messprintf ("Transcript %s, %s", 
				 name(seg1->parent),
				 cp),
	     maxBuf) ;
  messfree (cp) ;
}

/*********************************************************************/

void fMapcDNAReportTranscribedgene (char *buffer, SEG *seg1, int maxBuf)
{
  char *cp = strnew(buffer,0) ;
  strncpy (buffer, messprintf ("Gene %s, %s #%d", name(seg1->parent), 
			       cp, seg1->data.k ? KEYKEY(seg1->data.k) : 0), maxBuf) ;
  messfree (cp) ;
}

/***************************************************************************************/
/***************************************************************************************/

void fMapcDNAShowRNAi (LOOK look, float *offset)
{
  int ii, box = 0, ix = 0 ;
  SEG *seg ;
  float x = 0, y1, y2, ddx ;
  BOOL isDown ; 
  KEY rnai = 0 ;
  BUMP bump = 0 ;
  float yz ;
  int color, color2 ;
  char pheno [10000] ; /* hope large enough */
  KEY _nvp = str2tag ("No_obvious_phenotype") ;

  cDNAAlignInit () ; 

  if (*look->map->activeColName == '-')
    isDown = FALSE ;
  else
    isDown = TRUE ;

  for (ii = 1 ; ii < arrayMax(look->segs) ; ++ii)
    { seg = arrp(look->segs,ii,SEG) ;
      switch (seg->type)
	{
	case MRNAI: if (isDown) break ; else continue ; 
	case MRNAI_UP: if (!isDown) break ; else continue ;
	default: continue ;
	}

      if (class(seg->key) == _VMethod)
	continue ;

      if (!keyFindTag (seg->key, _Phenotype) &&
	  !keyFindTag (seg->key, _nvp))
	continue ;

      y1 = MAP2GRAPH(look->map,seg->x1 - .5) ;
      y2 = MAP2GRAPH(look->map,seg->x2+1 - .5) ;	/* to cover full base */

      if (y2 < .3 + topMargin ||
	  y1 > mapGraphHeight - .3)
	continue ;

      if (!bump)
	bump = bumpCreate (mapGraphWidth, 0) ;
      if (y1 < .2 + topMargin) 
	y1 = .2 + topMargin ;
      if (y2 > mapGraphHeight - .3)
	y2 = mapGraphHeight - .3;  

      yz = y1 + .1 ; ix = 0 ;   /* cheat a little to avoid bump by rounding */
      bumpItem (bump,1,(y2 - y1 - .2) , &ix, &yz) ;  

      x = *offset + .5 + .9 *ix ;

      rnai = seg->key ;
      array(look->boxIndex, box = graphBoxStart (), int) = ii ;
      
      ddx = .2 ;
      sprintf (pheno, "RNAi experiment") ;
      
      if (keyFindTag (rnai, _Phenotype))
	{
	  color = RED ; color2 = LIGHTRED ; 
	  if (y2 - y1 > 1)
	    {
	      KEY ph = keyGetKey (rnai, _Phenotype) ;
	      if (ph && strlen (name(ph)) < 8000) 
		sprintf (pheno, "RNAi experiement, %s", name (ph)) ;
	    }
	}
      else
	{ 
	  color = GRAY ; color2 = LIGHTGRAY ;
	  sprintf (pheno, "RNAi experiment yielding no obvious phenotype") ;
	}
      /* mixture was no good 
	 if (keyGetKey (rnai, _Gene))
	 color = RED ;
	 else
	 {
	 color = BLACK ;
	 if (color2 == LIGHTRED) color2 = RED ;
	 }  mixture was no good */
      
      graphRectangle (x , y1, x + 2 * ddx, y2) ;
  
      graphBoxEnd () ;
      
      graphBoxDraw (box, color, color2) ;

      fMapBoxInfo (look, box, seg) ;
      graphBoxFreeMenu(box, fMapGenesAction, fMapRnaiMenu) ; 

      graphBubbleInfo (box, seg->key, 0, pheno) ;
    }   
      
  *offset += bumpMax(bump) ? .9 * bumpMax (bump) + 1.0 : 0 ;
  bumpDestroy (bump) ;
} /* fMapcDNAShowRNAi */

/***************************************************************************************/
/***************************************************************************************/

void fMapcDNAShowOST (LOOK look, float *offset)
{
  int ii, box = 0, ix = 0 ;
  SEG *seg ;
  float x = 0, y1, y2, dy, ddx ;
  BOOL isDown ; 
  KEY ost = 0 ;
  BUMP bump = 0 ;
  BOOL hasOst = FALSE ;
  float yz ;
  int color, color2 ;
  char pheno [1000] ; /* hope large enough */

  cDNAAlignInit () ; 

  if (*look->map->activeColName == '-')
    isDown = FALSE ;
  else
    isDown = TRUE ;

  for (ii = 1 ; ii < arrayMax(look->segs) ; ++ii)
    { seg = arrp(look->segs,ii,SEG) ;
      switch (seg->type)
	{
	case MOST: if (isDown) break ; else continue ; 
	case MOST_UP: if (!isDown) break ; else continue ;
	default: continue ;
	}

      if (class(seg->key) == _VMethod)
	continue ;

      y1 = MAP2GRAPH(look->map,seg->x1 - .5) ;
      y2 = MAP2GRAPH(look->map,seg->x2+1 - .5) ;	/* to cover full base */

      if (y2 < .3 + topMargin ||
	  y1 > mapGraphHeight - .3)
	continue ;

      if (!bump)
	bump = bumpCreate (mapGraphWidth, 0) ;
      if (y1 < .2 + topMargin) 
	y1 = .2 + topMargin ;
      if (y2 > mapGraphHeight - .2)
	y2 = mapGraphHeight - .2 ;  

      yz = y1 + .1 ; ix = 0 ;   /* cheat a little to avoid bump by rounding */
      ix = 0 ; 
      if (0) /* they take too much space */ 
	bumpItem (bump,1,(y2 - y1 - .2) , &ix, &yz) ;  
      hasOst = TRUE ;
      x = *offset + .5 + .9 *ix ;

      ost = seg->key ;
      array(look->boxIndex, box = graphBoxStart (), int) = ii ;
      
      ddx = .4 ;
      dy = MAP2GRAPH(look->map,seg->x1 + 30)  - MAP2GRAPH(look->map,seg->x1)  ;
      if (dy < .8) dy = .8 ;
      if (y2 - y1 > 2 * dy)
	{
	  graphLine (x + ddx , y1 + dy, x + ddx, y2 - dy) ;
	  graphColor (BLACK) ;
	  graphLine (x + ddx , y1, x + ddx, y1 + dy) ;
	  graphLine (x + ddx , y2, x + ddx, y2 - dy) ;
	  graphLine (x , y1, x + 2 * ddx, y1) ;
	  graphLine (x , y2, x + 2 * ddx, y2) ;
	  graphLine (x, y1 + dy - .5, x+ddx, y1 + dy) ;
	  graphLine (x + 2 * ddx, y1 + dy - .5, x+ddx, y1 + dy) ;
	  graphLine (x, y2 - dy + .5, x+ddx, y2 - dy) ;
	  graphLine (x+ 2 * ddx, y2 - dy + .5, x+ddx, y2 - dy) ;
	}
      else
	{
	  graphColor (BLACK) ;
	  graphLine (x , y1, x + 2 * ddx, y2) ;
	  graphLine (x , y2, x + 2 * ddx, y1) ;
	  graphLine (x , y1, x + 2 * ddx, y1) ;
	  graphLine (x , y2, x + 2 * ddx, y2) ;
	}
      
      if (keyFindTag (ost,  str2tag ("Amplified")))
	{ 
	  color = DARKGREEN ; color2 = TRANSPARENT ;
	  sprintf (pheno, "Successful PCR amplification of cDNA library") ;
	}
      else
	{ 
	  color = DARKGRAY ; color2 =  TRANSPARENT /*LIGHTGRAY*/ ;
	  sprintf (pheno, "Failed PCR amplification of cDNA library") ;
	}

      graphBoxEnd () ;
      
      graphBoxDraw (box, color, color2) ;

      fMapBoxInfo (look, box, seg) ;
      graphBoxFreeMenu(box, fMapGenesAction, fMapRnaiMenu) ; 
      graphBubbleInfo (box, seg->key, "Vidal", pheno) ;
    }   
      
  if (bumpMax(bump))
    *offset += .9 * bumpMax (bump) + 1.0 ;
  else if (hasOst)
    *offset += 1.9 ;
  bumpDestroy (bump) ;
} /* fMapcDNAShowOST */

/***************************************************************************************/

void fMapcDNAShowGenes (LOOK look, float *offset)
{
  int n, ii, box = 0, ix = 0 ;
  SEG *seg ;
  float x = 0, y1, y2, oldTextHeight ;
  BOOL isDown, isArrow ; 
  KEY gene = 0 ;
  BUMP bump = 0 ; 
  float yz ;
  int color, color2 ;
  BOOL isCloud ;
  int GENE_COLOR = CYAN ;       /* RED ; ORANGE DARKGREEN ;  magenta for the web */
  int GENE_COLOR2 =  PALECYAN ; /* PALERED ;  YELLOW LIGHTGREEN ;  palemagenta for the web */
  KEY _Cloud_gene = str2tag ("Cloud_gene") ;
  cDNAAlignInit () ; 

  oldTextHeight = graphTextHeight (1.0) ;
  if (*look->map->activeColName == '-')
    isDown = FALSE ;
  else
    isDown = TRUE ;

  for (ii = 1 ; ii < arrayMax(look->segs) ; ++ii)
    { seg = arrp(look->segs,ii,SEG) ;
      switch (seg->type)
	{
	case MGENES: if (isDown) break ; else continue ; 
	case MGENES_UP: if (!isDown) break ; else continue ;
	default: continue ;
	}

      if (class(seg->key) == _VMethod)
	continue ;

      isCloud = keyFindTag (seg->key, _Cloud_gene) ;
      /*
	if (!keyFindTag (seg->key, _Transcribed_gene) &&
	!keyFindTag (seg->key, _Phenotype))
	continue ;
      */

      y1 = MAP2GRAPH(look->map,seg->x1 - .5) ;
      y2 = MAP2GRAPH(look->map,seg->x2 + .5) ;	/* to cover full base */

      if (y2 < .3 + topMargin ||
	  y1 > mapGraphHeight - .3)
	continue ;

      if (!bump)
	{
	  *offset += .3 ;
	  bump = bumpCreate (mapGraphWidth, 0) ;
	}
      if (y1 < .2 + topMargin) 
	y1 = .2 + topMargin ;
      if (y2 > mapGraphHeight -.2)
	y2 = mapGraphHeight - .2 ;  
      
      yz = y1 + .1 ; ix = 0 ;   /* cheat a little to avoid bump by rounding */
      bumpItem (bump,1,(y2 - y1 - .2) , &ix, &yz) ;  

      x = *offset + 1.7 * ix ;
      gene = seg->key ;
      
      array(look->boxIndex, box = graphBoxStart (), int) = ii ;

      isArrow = FALSE ; 
      if (keyFindTag (gene, _Title) || keyFindTag (gene, _NewName))
	{
	  color = GENE_COLOR ; color2 = GENE_COLOR2 ;
	  if (isCloud) 
	    {
	      color2 = LIGHTGRAY ;
	      if (isDown)
		graphRectangle (x + .3 , y1, x + 1.1, y2 - .2) ;
	      else
		graphRectangle (x + .3, y1 + .2 , x + 1.1, y2) ;
	    }
	  else
	    {
	      if (isDown)
		graphRectangle (x , y1, x + 1.4, y2 - .4) ;
	      else
		graphRectangle (x , y1 + .4 , x + 1.4, y2) ;
	    }
	  isArrow = TRUE ;
	}
      else if (keyFindTag (gene, str2tag("GeneBox")))
	{
	  color = GREEN7 ; color2 = GREEN4 ;
	  graphRectangle (x , y1, x + .4, y2) ;
	}
      else
	{
	  color = GENE_COLOR ; color2 = WHITE ; /* PALEGRAYGRAY ; */
	  graphRectangle (x , y1, x + 1.4, y2) ;
	}
      graphBoxEnd () ;

      if (keyFindTag (gene, _Colour))
	{
	  KEY geneCol = color2 ;
	  OBJ Gene = bsCreate (gene) ;
	  if (Gene && bsGetKeyTags (Gene, _Colour, &(geneCol)))
	    color2 = geneCol - _WHITE + WHITE ;
	  bsDestroy (Gene) ;
	}

      graphBoxDraw (box, color, color2) ;

      n = strlen(name(gene)) ;
      if (y2 - y1 > n + 3)
	{
	  int oldFormat = graphTextFormat(BOLD) ;
	  graphColor (BLACK) ;
	  graphTextUp (name(gene), x + .25, (y1 + y2)/2 + n/2) ; 
	  graphTextFormat(oldFormat) ;
	}

            
      if (isCloud && isArrow)
	{
	  if (isDown)
	    {
	      graphColor (BLUE) ;
	      graphLine (x - .1, y2 - .5, x + .7, y2) ;
	      graphLine (x + .7, y2, x + 1.5, y2 - .5) ;
	      graphColor (BLACK) ;
	    }
	  else
	    {
	      graphColor (BLUE) ;
	      graphLine (x - .1, y1 + .5, x + .7, y1) ;
	      graphLine (x + .7, y1, x + 1.5, y1 + .5) ;
	      graphColor (BLACK) ;
	    }
	}
      else if (! isCloud && isArrow)
	{
	  if (isDown)
	    {
	      graphColor (BLUE) ;
	      graphLine (x - .4, y2 - .7, x + .7, y2) ;
	      graphLine (x + .7, y2, x + 1.8, y2 - .7) ;
	      graphColor (BLACK) ;
	    }
	  else
	    {
	      graphColor (BLUE) ;
	      graphLine (x - .4, y1 + .7, x + .7, y1) ;
	      graphLine (x + .7, y1, x + 1.8, y1 + .7) ;
	      graphColor (BLACK) ;
	    }
	}

      fMapBoxInfo (look, box, seg) ;
      graphBoxFreeMenu(box, fMapGenesAction, fMapGenesMenu) ; 
      if (isGifDisplay)
	{
	  KEY title = keyGetKey (seg->key, _Title) ;
	  graphBubbleInfo (box, seg->key, className(seg->key), 
			   messprintf("%s", title ? name(title) : "")
			   ) ;
	}
    }   
 
  if (bumpMax (bump))
    *offset += 1.7 * bumpMax (bump) + .5 ;
  bumpDestroy (bump) ;
  graphTextHeight (oldTextHeight) ;
  graphColor (BLACK) ;
} /* fMapcDNAShowGenes */

/***************************************************************************************/
/***************************************************************************************/

static void fMapcDNADecorateTranscript (LOOK look, int color, int a1, int a2) 
{
  char *cp = arrp (look->colors, a1, char) ;
  int n = a2 - a1;

  if (a2 > a1)
    {
      n = a2 - a1 + 1 ;
      cp = arrp (look->colors, a1, char) ;
      while (n--) 
	*cp++ |= color ;
    }
  else
    {
      n = a1 - a2 + 1 ;
      cp = arrp (look->colors, a1, char) ;
      while (n--) 
	*cp-- |= color ;
    }
}

/*********/

void fMapcDNAShowTranscript (LOOK look, float *offset)
{
  int ii, box = 0, ix = 0, oldix = 0, color = 0 ;
  SEG *seg ;
  float x = 0, y1, y2 ;
  BOOL isDown ;
  KEY transcript = 0, type, subtype, ts ;
  BUMP bump = 0 ; 
  float yz ;
  int TR_COLOR = DARKGREEN ; /* magenta for the web */
  int TR_COLOR2 =  PALEGREEN ; /* palemagenta for the web */
  cDNAAlignInit () ; 
  



  if (*look->map->activeColName == '-')
    { 
      isDown = FALSE ;
      if (FALSE && mapColSetByName ("cDNA_Transcript_Decorate", -1))
	fMapClearDNA (look) ;
    }
  else
    { isDown = TRUE ;
    }

  for (ii = 1 ; ii < arrayMax(look->segs) ; ++ii)
    { seg = arrp(look->segs,ii,SEG) ;
      switch (seg->type)
	{
	case TRANSCRIPT: if (isDown) break ; else continue ;
	case TRANSCRIPT_UP: if (!isDown) break ; else continue ;
	default: continue ;
	}

      if (class(seg->key) == _VMethod)
	continue ;

      y1 = MAP2GRAPH(look->map,seg->x1 - .5) ;
      y2 = MAP2GRAPH(look->map,seg->x2+1 - .5) ;	/* to cover full base */

      if (y2 < .3 + topMargin ||
	  y1 > mapGraphHeight - .3)
	continue ;

      if (!bump)
	bump = bumpCreate (mapGraphWidth, 0) ;
      if (y1 < .2 + topMargin) 
	y1 = .2 + topMargin ;
      if (y2 > mapGraphHeight - .2)
	y2 = mapGraphHeight - .2 ;  

      type = seg->key ;
      subtype = seg->data.i ; 
      yz = y1 + .1 ; ix = 0 ;   /* cheat a little to avoid bump by rounding */
      bumpItem (bump,1,(y2 - y1 - .2) , &ix, &yz) ;  
      if (transcript == seg->parent && /* do not move inside a given transcript */
	  ix != oldix)
	ix = oldix ;

      oldix = ix ; 
      x = *offset + 1.4 *ix ;

      if (isDown && transcript != seg->parent &&  /* entering a new transcript */
	  (ts = keyGetKey (seg->parent, _Transpliced_to)) &&
	  strcmp (name(ts), "SL0") &&
	  ! strncmp (name(ts), "SL", 2))
	{
	  int bx1 = graphBoxStart () ;
	  graphText (name(ts) + 2, x, y1 - 1) ;
	  graphBoxEnd () ;
	  graphBubbleInfo (bx1, 0, 0, messprintf ("Transpliced leader %s", name(ts))) ;
	}
      
      transcript = seg->parent ;
      
      array(look->boxIndex, box = graphBoxStart (), int) = ii ;
      graphColor  (TR_COLOR) ; 
      
      if (strstr (name(type), "ntron"))
	{
	  float mid = (y1 + y2) / 2 ;
	  if (subtype)
	    {
	      if (subtype == _gt_ag) ;
	      else if (subtype == _gc_ag) ;
	      else
		graphColor (BLUE) ;
	    }
	  /*
	    graphLine (x +.2, y2, x + 1.6 , mid) ;
	    graphLine (x +1.6, mid, x + .2, y1 +.1) ;
	  */
	    graphLine (x + 1.1, y2, x + 2.6 , mid) ;
	    graphLine (x + 2.6, mid, x + 1.1, y1) ;
	    color = TRANSPARENT ; 
	}
      else if (strstr (name(type), "xon"))
	{ 
	  graphRectangle (x +.2, y1, x + 1.4 , y2) ;
	  color = TR_COLOR2 ;
	  if (FALSE && mapColSetByName ("cDNA_Transcript_Decorate", -1))
	    fMapcDNADecorateTranscript (look, TINT_LIGHTGREEN, seg->x1, seg->x2) ;
	}
      else if (type == _UTR_3prime ||
	       type == _UTR_5prime)
	{ 
	  graphRectangle (x +.4, y1, x + 1.2 , y2) ;
	  color = WHITE ;
	}
      else if (type == _Gap)
	{ 
	  int oldColor = graphColor (BLACK) ;
	  graphLine (x +.8, y1, x + .8 , y2) ; 
	  color = oldColor ;
	}
      graphColor (BLACK); 
      graphBoxEnd () ;
      
      graphBoxDraw (box, BLACK, color) ;
      /* avant graphBoxInfo (box, transcript, info) ;*/
      /* mhmp 30.09.99
	 graphBoxInfo (box, transcript, name(seg->key)) ;*/
      fMapBoxInfo (look, box, seg) ;
      graphBoxFreeMenu(box, fMapTranscribedgeneAction, fMapTranscriptMenu) ; 
      graphBubbleInfo (box, seg->parent, "Transcribed_gene", "Transcribed gene") ;
    }      

  *offset += (1.4 + 1.0) * bumpMax (bump) + .5 ;
  bumpDestroy (bump) ;
  graphColor (BLACK) ;
}   /* fMapcDNAShowTranscript */

/*********/

void fMapcDNADoShowMrna (LOOK look, float *offset, BOOL isPredicted)
{
  int ii, box = 0, ix = 0, color = 0 ;
  SEG *seg ;
  float x = 0, y1, y2 ;
  BOOL isDown, useAm = FALSE ; 
  KEY mrna = 0, bumpedParent = 0, type, subtype, ts, pg ;
  BUMP bump = 0, bumpIntrons = 0 ; 
  float yz ;
  KEY _Well_supported ;
  int TR_COLOR = DARKGREEN  ;/*  MAGENTA ;  ORANGE DARKGREEN ;  magenta for the web */
  int TR_COLOR2 = LIGHTGREEN ; /* PALEMAGENTA ;  YELLOW LIGHTGREEN ;  palemagenta for the web */
  int TR_COLOR3 =  PALEBLUE ; /* YELLOW LIGHTGREEN ;  palemagenta for the web */
  int TR_COLOR4 =  PALEORANGE ; /* YELLOW LIGHTGREEN ;  palemagenta for the web */
  KEY _Valid3p = str2tag ("Valid3p") ;
  KEY _Valid5p = str2tag ("Valid5p") ; /* RED YELLOW */

  cDNAAlignInit () ; 

  look->mrnaOffset = *offset ;
  if (0 && look->view && isPredicted)
    return ;

  _Well_supported= str2tag ("Well_supported") ;

  if (*look->map->activeColName == '-')
    { 
      isDown = FALSE ;
      if (FALSE && mapColSetByName ("cDNA_Transcript_Decorate", -1))
	fMapClearDNA (look) ;
    }
  else
    { 
      isDown = TRUE ;
    }

  for (ii = 1 ; ii < arrayMax(look->segs) ; ++ii)
    { seg = arrp(look->segs,ii,SEG) ;
      switch (seg->type)
	{
	case MRNA: if (isDown && !isPredicted) break ; else continue ; 
	case MRNA_UP: if (!isDown && !isPredicted) break ; else continue ;
	case PMRNA: if (isDown && isPredicted) break ; else continue ; 
	case PMRNA_UP: if (!isDown && isPredicted) break ; else continue ;
	default: continue ;
	}

      if (class(seg->key) == _VMethod)
	continue ;

      if (class(seg->key))
	continue ;

      if (!class(seg->parent))
	continue ;

      if (seg->x1 > look->max || seg->x2 < look->min)
	continue ;

      y1 = MAP2GRAPH(look->map,seg->x1 - .5) ;
      y2 = MAP2GRAPH(look->map,seg->x2 + .5) ;	/* to cover full base */
      if (y1 < .2 + topMargin) y1 = .2 + topMargin ;
      if (y2 < .3 + topMargin) continue ;
      if (y2 > mapGraphHeight - .2)	y2 = mapGraphHeight- .2 ;  
      if (y1 > mapGraphHeight - .3) continue ;

      type = seg->key ;
      subtype = seg->data.i ; 

      /* color selection */
      if (isPredicted)
	{
	  if ((pg = keyGetKey (seg->parent, _From_prediction)))
	    {
	      if (!strncmp (name(seg->parent), "hmm", 3))
		{  TR_COLOR = LIGHTGREEN ; TR_COLOR2 = WHITE ; }
	      else if (keyFindTag (pg, str2tag ("has_exact_est"))) 
		{ TR_COLOR = MAGENTA ; TR_COLOR2 = PALEMAGENTA ; }
	      else
		{ TR_COLOR = BLUE ; TR_COLOR2 =  PALEBLUE ; }
	    }
	}
      else
	{
	  if (keyFindTag (seg->parent, _From_AM))
	    { TR_COLOR = RED ; TR_COLOR2 = LIGHTMAGENTA ; } 
	  else if (keyFindTag (seg->parent, _Bad_quality))
	    { TR_COLOR = ORANGE ; TR_COLOR2 = YELLOW ; }
	  else if (keyFindTag (seg->parent, _Colour))
	    { 
	      KEY myCol = 0 ;
	      OBJ obj = bsCreate(seg->parent) ;
	      bsGetKeyTags (obj,_Colour, &myCol) ;
	      bsDestroy (obj) ;
	      TR_COLOR =  TR_COLOR2 = WHITE + myCol - _WHITE ;
	    }
	  else
	    { TR_COLOR = MAGENTA ; TR_COLOR2 = LIGHTMAGENTA ; } 
	}

      if (seg->parent != bumpedParent) /* do not move inside a given mrna */
	{
	  float ytop = 0, ybottom = 0 ;
	  
	  useAm = keyFindTag (seg->parent, _From_AM) ;

	  bumpedParent = seg->parent ;
	  ytop = y1 + .1 ; ix = 0 ;   /* cheat a little to avoid bump by rounding */
	  ybottom = MAP2GRAPH (look->map, seg->x1 + seg->sourceDx - .5) - .1 + 1.2 ;
	  if (ybottom > mapGraphHeight - .2)	ybottom = mapGraphHeight - .2 ;   

	  if (!bump)
	    bump = bumpCreate (mapGraphWidth, 0) ;
	  yz = ytop ;
	  if (ybottom > ytop)
	    bumpItem (bump,1, ybottom - ytop,&ix,&yz) ;
	  /* printf("\nbumping %s %g %g ix=%d\n", name(bumpedParent), ytop, ybottom, ix) ;   */
	}
      x = *offset + 2.4 * ix ;

      if (!isDown && ix > 3)
	continue ;

      if (isDown && mrna != seg->parent &&  /* entering a new mrna */
	  keyFindTag (seg->parent, _Transpliced_to))
	{
	  int i ;
	  static Array units = 0 ;
	  BSunit *uu ;
	  OBJ Parent = bsCreate (seg->parent) ;

	  units = arrayReCreate (units, 12, BSunit) ;
	  if (bsGetArray (Parent, _Transpliced_to, units, 2))
	    for ( i = 0 ; i < arrayMax(units) ; i += 2)
	      {
		uu = arrp (units, i, BSunit) ;
		ts = uu[0].k ;
		if (strcmp (name(ts), "SL0") &&
		    ! strncmp (name(ts), "SL", 2))
		  graphText (name(ts) + 2, x, y1 - 1) ;
	      }
	  bsDestroy (Parent) ;
	  arrayDestroy (units) ;
	}
      if (0 && isDown && mrna != seg->parent)  /* entering a new mrna */
	{
	  KEYSET ks = queryKey (seg->parent, "Follow Product ; Follow kantor ; Blastp_date && (!blastp || Taxblast_date) && Pfam_date && Psort_date && Expasy_date && Kantor_date") ;
	  
	  if (keySetMax(ks))
	    { 
	      int oldFormat = graphTextFormat(BOLD) ;
	      graphText ("f", x, y1 - 2) ;
	      graphTextFormat(oldFormat) ;
	    }
	  keySetDestroy (ks) ;
	}	  
      
      mrna = seg->parent ;
      
      array(look->boxIndex, box = graphBoxStart (), int) = ii ;
      graphColor (TR_COLOR) ; 

      type = seg->key ;
      if (strstr (name(type), "ntron") &&  ! strstr (name(subtype), "Fuzzy"))
	{
	  float mid = (y1 + y2) / 2 ;

	  /* if (type == _Stolen_intron)
	     graphColor (CERISE) ;
	     else */ if (subtype)
	       {
		 if (subtype == _gt_ag) ;
		 else if (subtype == _gc_ag) ;
		 else
		   graphColor (BLUE) ;
		}
	  if (y2 > y1 + .2)
	    {
	      graphLine (x +.4, y2, x + 1.4 , y2) ; /* give width to the intron box */
	      graphLine (x +1.4, y2, x + 2.4 , mid) ;
	      graphLine (x +2.4, mid, x + 1.4, y1) ;
	    }
	  else
	    graphLine (x +1.4, mid, x + 2.4 , mid) ;
	  
	  color = TRANSPARENT ;
	  {
	    int z;
	    void *vp ;
	    
	    z =  ((171 * (seg->x1)) % 30269) ^ ((172 * seg->x2) % 30307) ; /* some large int */
	    if (z == 0) z = 1 ; /* unlucky ! */
	    if (graphAssFind (assVoid (z), &vp))
	      color = intronCol [assInt (vp) % numIntronCol] ;
	    else
	      {
		int intronX = ix ; 
		float intronY = y1 ;
		
		if (!bumpIntrons)
		  bumpIntrons = bumpCreate (300, 0) ;
		bumpItem (bumpIntrons,1,(y2 - y1) - .06, &intronX, &intronY) ;
		color = intronCol [intronX % numIntronCol] ;
		graphAssociate (assVoid(z), assVoid(intronX % numIntronCol)) ;
	      }
	  }
	}
      else if (type == _ORF_Gap || type == _Stolen_Exon)
	{ 
	  graphRectangle (x +.4, y1, x + 1.4 , y2) ;
	  color = PALEMAGENTA ; /* PALE_GRAY */
	}
      else if (type == _Predicted_Exon)
	{ 
	  graphRectangle (x +.4, y1, x + 1.4 , y2) ;
	  color = LIGHTBLUE ; /* PALE_GRAY */
	}
      else if (strstr (name(type), "xon"))
	{ 
	  if (subtype & 0x2)
	    graphRectangle (x +.4, y1, x + 1.4 , y2) ;
	  else
	    graphRectangle (x +.6, y1, x + 1.2 , y2) ;
	  switch (subtype)
	    {
	      /* control the colors of the successive orf of a given frame */
	      /* just now, they are all same color */
	    case 0: color = WHITE ; break ;  /* utr */
	    case 1: color = WHITE ; break ;  /* non best product */
	    case 2:
	      color = TR_COLOR2 ; color = YELLOW ; break ;
	    case 3:
	      color = TR_COLOR2 ; color=RED ; break ; /* best product */
	    case 6: color = TR_COLOR3 ; color=GREEN ; break ;
	    case 7: color = TR_COLOR4 ; color=ORANGE ; break ;
	    case 10: color = TR_COLOR4 ; color=PALECYAN ; break ;
	    case 14: color = TR_COLOR4 ; color=GREEN3 ; break ;
	    case 15: color = TR_COLOR4 ; color=PALEYELLOW ; break ;
	    }
	  if (0 && useAm)
	    color = TR_COLOR ;	      
	}
      else if (subtype == _Fuzzy_gt_ag || subtype == _Fuzzy_gc_ag) 
	{
	  graphLine (x +.6, y2, x + 1.2 , y2) ; 
	  graphLine (x +.9, y1, x + .9 , y2) ; 
	  color = TRANSPARENT ;
	}
      else if (strstr (name(subtype), "Fuzzy"))
	{
	  graphColor (BLUE) ;
	  graphLine (x +.6, y2, x + 1.2 , y2) ; 
	  graphLine (x +.9, y1, x + .9 , y2) ; 
	  color = TRANSPARENT ;
	}
      else if (type == _Gap)
	{ 
	  float mid = (y1 + y2) / 2 ;
	  graphColor (BLACK) ; 

	  graphLine (x +.9, y1, x + .9 , y2) ; 
	  if (1) graphRectangle (x + .6, mid - .2, x + 1.2, mid + .2) ;
	  color = TRANSPARENT ;
	}
      else if (type == _Valid3p)
	{ 
	  int n = seg->data.i ;
	  if (n < 0) { graphColor (BLUE); n = -n ; }
	  else graphColor (BLACK); 
	  graphLine (x + 1.2, y2, x + 2.4 , y2) ; 
	  if (n > 5)
	    graphFillRectangle (x + 2.0, y2 + (isDown ? -.5 : .5), x + 2.4 , y2) ; 
	  else
	    graphRectangle (x + 2.1, y2 + (isDown ? -.4 : .4), x + 2.4 , y2) ; 
	  color = TRANSPARENT ;
	}
      else if (type == _Valid5p)
	{ 
	  if ((seg->data.i >> 24) >= 1)
	    graphColor (BLACK); 
	  else
	    graphColor (BLUE);
	  graphLine (x + 1.2, y1, x + 2.4 , y1) ; 
	  if ((seg->data.i >> 24) >= 2)
	    graphFillRectangle (x + 2.0, y1 + (isDown ? +.5 : -.5), x + 2.4 , y1) ; 
	  else	
	    graphRectangle (x + 2.1, y1 + (isDown ? +.4 : -.4), x + 2.4 , y1) ; 
	  color = TRANSPARENT ;
	}
      graphColor (BLACK); 
      graphBoxEnd () ;
      
      graphBoxDraw (box, BLACK, color) ;
      /* mhmp 30.09.99
	 graphBoxInfo (box, mrna, name(seg->key)) ;*/
      fMapBoxInfo (look, box, seg) ;
      graphBoxFreeMenu(box, fMapMrnaAction, fMapMrnaMenu) ; 
      /* geneGoodName(seg->parent) */
      {
	char *war = 0, *var = name(seg->parent) + strlen(name(seg->parent)) - 2 ;

	if (strlen(name(seg->parent)) > 2 && *var == '.')
	  { var++ ; war = " variant " ; }
	else
	  var = 0 ;
	
	if (strstr (name(type), "xon"))
	  {
	     graphBubbleInfo (box, seg->parent,  className(seg->parent),
			       messprintf("%d bp exon", seg->x2 - seg->x1 + 1)
			      ) ;
	     /*
	       if (subtype)
	       graphBubbleInfo (box, seg->parent,  className(seg->parent), 
	       messprintf("%s%s Exon %d bp", var ? war : "", var ? var : "", seg->x2 - seg->x1 + 1)) ;
	       else
	       graphBubbleInfo (box, seg->parent,  className(seg->parent), 
	       messprintf("%s%s UTR %d bp", var ? war : "", var ? var : "", seg->x2 - seg->x1 + 1)) ;
	     */
	  }
	else if (type == _Valid3p)
	  {
	    int n = seg->data.i ;
	    if (n < 0) n = -n ;
	    graphBubbleInfo (box, 0, 0
			     , messprintf("Validated 3'end, %d accession%s", n, n > 1 ? "s" : "")
			     ) ;
	  }
	else if (type == _Valid5p)
	  {
	    int n =  seg->data.i & 0xffffff ;
	    int nt = seg->data.i >> 24 ;

	     if (nt >= 2)
	       graphBubbleInfo (box, 0, 0
				, messprintf("Capped 5' end, %d accession%s", n, n > 1 ? "s" : "")
				) ;
	     
	     else if (n >= 1)
	       graphBubbleInfo (box, 0, 0
				, messprintf("5' end, %d aggregated accession%s",  n, n > 1 ? "s" : "")
			   ) ;
	  }
	else /* intron */
	  {
	    char *cheat = "" ; /* Stolen_Intron cheat, do not show them on web we have too many in build 31 */
	    if (type && name(type)) cheat = name(type) ;
	    if (strstr(cheat, "tron") )
	      cheat = "Intron" ;
	      
	    graphBubbleInfo (box, seg->parent,  className(seg->parent), 
			     messprintf("%s%s %s %s %d bp", var ? war:"", var ? var : "", cheat, subtype > _Date ? name(subtype) : "", seg->x2 - seg->x1 + 1)) ;
	  }

      }
      if (isDown && (seg->type != (seg+1)->type || mrna != (seg+1)->parent))  /* leaving an mrna i hope */
	{
	  char buf[2], *cp = name (seg->parent), *cq ;
	  int i = strlen (cp) ;
	  BOOL well = keyFindTag (seg->parent, _Well_supported) ;

	  cq = cp + i - 1 ;
	  while (cp < cq && *cq != '.') cq-- ;
	  if (*cq++ == '.' && strlen (cq) < 8)
	    {
	      buf[0] = *cq ; buf[1] = 0 ;
	      graphText (buf, x + .5, y2 + .0) ;
	      if (well)
		{
		  graphLine (x + .5, y2 + 1.0, x + 1.5, y2 + 1.0) ;
		}
	    }
	}      
    }   
  
  if (bumpMax (bump))
    *offset += 2.4 * (!isDown && bumpMax (bump) > 3 ? 3 : bumpMax (bump) ) + 1.5 ;
  bumpDestroy (bump) ;
  bumpDestroy (bumpIntrons) ;
  graphColor (BLACK) ;
} /* fMapcDNADoShowMrna */

void fMapcDNAShowMrna (LOOK look, float *offset)
{
  fMapcDNADoShowMrna (look, offset, FALSE) ;

}

void fMapcDNAShowPMrna (LOOK look, float *offset)
{
  fMapcDNADoShowMrna (look, offset, TRUE) ;

}

/*********/

static BOOL fMapcDNAShowStopsInFrame (LOOK look, float *offset, BOOL isDown, int frame, int color)
{
  float x1, x2, y1, y2 ;
  char *cp, cc ;
  int box, globalBox, ii, i1, i2, di ;
  int iStop, iMet, jj, j1, j2 ;
  int STOP_COLOR = BLACK ;
  int M_COLOR = GREEN ;
  char *translationTable = pepGetTranslationTable (look->seqKey, 0) ; 
  static Array stops = 0, mets = 0 ;

  if (! look->dna)
    return FALSE ;
  if (!stops)
    stops = arrayReCreate (stops, 200, int) ;
  if (!mets)
    mets = arrayReCreate (mets, 200, int) ;
  x1 = *offset ;
  x2 = *offset + .6 ;
  *offset += .8 ;
  y1 = .2 + topMargin ;
  y2 = mapGraphHeight - .2 ; 
  
  i1 = GRAPH2MAP (look->map, y1) ;
  i2 = GRAPH2MAP (look->map, y2) ;
  
  i1 = 3 * ( i1/3) + frame - 1 ;
  i2 = 3 * ( i2/3) + frame - 1 ;
  
  if (i1 < 0) i1 = frame - 1 ;
  if (i1 > arrayMax(look->dna) - 3)
    return FALSE ;
  if (i2 < 0) 
    return FALSE ;
  if (i2 > arrayMax(look->dna) - 2)
    i2 = arrayMax(look->dna) - 2 ;

  y1 = MAP2GRAPH(look->map,i1 - .5) ;
  y2 = MAP2GRAPH(look->map,i2 - .5) ;
  
  /* graphRectangle (x1, y1, x2, y2) ; */

  globalBox = graphBoxStart () ;
  if (isDown)
    {
      di = 3 ; 
      cp = arrp (look->dna, i1, char) ;
    }
  else
    {
      di = - 3 ; 
      cp = arrp (look->dnaR, arrayMax(look->dna) - i1 - 3, char) ;
    }


  for (iStop = iMet = 0, ii = i1; ii < i2 ; ii += 3, cp += di)
    {
      /* do not call codon here, too slow */
      cc = e_codon (cp, translationTable) ;
      
      switch (cc)
	{
	case '*': 
          array (stops, iStop++, int) = ii ;
          break ;
	case 'M': 
	  array (mets, iMet++, int) = ii ;
          break ;
	case 'X': 
          array (stops, iStop++, int) = ii | 0x40000000 ;
          break ;
	}
    }

  /* draw the ORF boxes */
  y1 = MAP2GRAPH (look->map,i1 - .5) ;
  j1 = i1 ;
  for (jj = 0 ; jj <= iStop ; jj++)
    {
      j2 = jj < iStop ? arr (stops, jj, int) : i2 ;
      j2 &= 0x00ffffff ;
      if (j2 - j1 > 300) /* bp */
	{
	  y1 = MAP2GRAPH (look->map,j1 - .5) ;
	  y2 = MAP2GRAPH (look->map,j2 + .5) ;
	  box = graphBoxStart() ;
	  graphRectangle (x1, y1, x2, y2) ;
	  graphBoxEnd () ;
	  graphBoxDraw (box, LIGHTGRAY, WHITE) ;
	  graphBubbleInfo (box, 0, 0, messprintf("ORF %d bp, frame %d, ATG (green) and stops (black)", j2 - j1 - 3, frame)) ;
	}
      j1 = j2 + 3 ;
    }
    
  /* draw the Met boxes */
  for (jj = 0 ; jj < iMet ; jj++)
    {
      j2 = arr (mets, jj, int) ;
      y1 = MAP2GRAPH(look->map,j2 - .5) ;
      y2 = MAP2GRAPH(look->map,j2 + 2.5) ;
      if (color != M_COLOR)
	graphColor (M_COLOR) ;
      color = M_COLOR ;
      box = graphBoxStart() ;
      graphFillRectangle (x1, y1, x2, y2) ; /*  (x1-.5, y1, x2+.3, y2) ; */
      graphBoxEnd () ;
      graphBubbleInfo (box, 0, 0, messprintf("ATG in frame %d at %d", frame, j2)) ;
    }
  
  /* draw the stop boxes after the ATG box */
  for (jj = 0 ; jj < iStop ; jj++)
    {
      j2 = arr (stops, jj, int) ;
      if (j2 & 0x40000000)
	continue ;
      y1 = MAP2GRAPH(look->map,j2 - .5) ;
      y2 = MAP2GRAPH(look->map,j2 + 2.5) ;
      if (color != STOP_COLOR)
	graphColor (STOP_COLOR) ;
      color = STOP_COLOR ;
      box = graphBoxStart() ;
      graphFillRectangle (x1, y1, x2, y2) ;
      graphBoxEnd () ;
      graphBubbleInfo (box, 0, 0, messprintf("Stop in frame %d at %d", frame, j2)) ;
    }

  graphBoxEnd () ;
  graphBubbleInfo (globalBox, 0, 0, "3 frames translation: Met in green, stops in black") ;
  graphColor (BLACK) ;

  return TRUE ;
} /* fMapcDNAShowStopsInFrame */

/************************************************************/
/* x, y position de la pointe */
static void fmapcdnaDrawTriangle (int box, float x, float y)
{
  int dx = 0, dy;
  float size = .8 ; /* size of the triangle */
  float d2x = 0, d2y = 0 ;
  float width = 0, height = 0 ;
  float aspect = 1.625; /* default aspect of the fonts */

  if (graphActive()) 
    {
      graphTextInfo(&dx, &dy, &width, &height); 
      if(dx) aspect =  ((float) dy)/dx;
    }
  /* triangle equilateral de cote size*0.8 hauteur de caractere */
  d2x = size * aspect * 0.866 ;
  d2y = size * 0.5 ;

  graphLine (x - d2x, y + d2y     , x - d2x, y - d2y) ;
  graphLine (x - d2x, y + d2y     , x, y) ;
  graphLine (x - d2x, y - d2y     , x, y) ;
}

/************************************************************/

static BOOL fMapcDNAShowMProductInFrame (LOOK look, float *offset, BOOL isDown, int frame
					 , int TR_COLOR, int TR_COLOR2)
{
  int ii, x1, x2, box = 0, ln, color, bgColor ;
  SEG *seg ;
  OBJ Mrna = 0 ;
  float x = 0, y1, y2 ;
  KEY product = 0, type, subtype, title = 0, mrna = 0 ;
  BOOL foundIt = FALSE ;
  KEY pg ;
  KEY _Valid3p = str2tag ("Valid3p") ;
  KEY _Valid5p = str2tag ("Valid5p") ;
  BOOL isUorf = FALSE ;
  /*
     if (frame)
     fMapcDNAShowStopsInFrame (look, offset, isDown, frame, TR_COLOR) ;
     *offset += 2 ;
     */
  for (ii = 1 ; ii < arrayMax(look->segs) ; ++ii)
    { 
      int il = 0 ;
      isUorf = FALSE ;
      seg = arrp(look->segs,ii,SEG) ;
      switch (seg->type)
	{
	case MPRODUCT: if (isDown) break ; else continue ; 
	case MPRODUCT_UP: if (!isDown) break ; else continue ;
	default: continue ;
	}
      
      if ((mrna  = keyGetKey (seg->parent, _mRNA)) &&
	  (pg = keyGetKey (mrna, _From_prediction)) &&
	  !keyFindTag (pg, str2tag ("has_exact_est")))
	{ TR_COLOR = BLUE ; TR_COLOR2 =  LIGHTBLUE ; }
      else if (keyFindTag (seg->parent, str2tag ("Best_product")) &&
	       keyFindTag (seg->parent, str2tag ("Good_product")))
	{ TR_COLOR = MAGENTA ; TR_COLOR2 = PALEMAGENTA ; }
      else if (keyFindTag (seg->parent, str2tag ("Very_good_product")))
	{ TR_COLOR = MAGENTA ; TR_COLOR2 = PALEVIOLET ; }
      else if (keyFindTag (seg->parent, str2tag ("Good_product")))
	{ 
	  TR_COLOR = MAGENTA ; TR_COLOR2 = PALEVIOLET ;
	  if (isGifDisplay)
	    continue ;
	}
      else if (keyFindTag (seg->parent, str2tag ("uORF_candidate")))
	{ 
	  TR_COLOR = MAGENTA ; TR_COLOR2 = PALEGREEN ;
	  isUorf = TRUE ;
	  if (0 && isGifDisplay)
	    continue ;
	}
      else
	{ 
	  TR_COLOR = MAGENTA ; TR_COLOR2 = PALEYELLOW ;
	  if (0 && isGifDisplay)
	    continue ;
	}
      if (frame == 1 && !title && class(seg->key) == _VmProduct)
	{
	  int oldT = graphTextFormat (BOLD) ;
	  KEY bestProd = seg->key ;
	  
	  if ((mrna = keyGetKey (seg->parent, _mRNA)) &&
	      keyFindTag (mrna, _From_AM) &&
	      (Mrna = bsCreate (mrna)))
	    {
	      int nerr = 0 ;
	      bsGetData (Mrna, _Tiling_error, _Int, &nerr) ;
	      if (nerr)
		graphText ("THIS mRNA DIFFERS FROM THE GENOME", 14, topMargin + il++) ;
	      bsDestroy (Mrna) ; 
	    }
	  if (mrna)
	    {
	      KEYSET prods = queryKey (mrna, ">product ; best_product") ;
	      if (keySetMax(prods))
		bestProd = keySet (prods, 0) ;
	      keySetDestroy (prods) ;
	    } 
	  title = 0 ;
	  if (!isGifDisplay)
	    {
	      title = keyGetKey (bestProd, _Title) ;
	      if (!title)
		title = keyGetKey (bestProd, _Blastp_title) ;
	    }
	  if (title)
	    {
	      float oldH = graphTextHeight (3.0) ;
	      
	      if (5*strlen(name(title)) < 3*mapGraphWidth - 5)
		graphText (name(title), 15, topMargin + 1) ;
	      else
		{
		  char **cpp = uBrokenLines (name(title), (3*mapGraphWidth - 5)/5.0) - 1;
		  while (*++cpp) graphText (*cpp, 14, topMargin + il++) ;
		}
	      
	      graphTextHeight (oldH) ;
	    }
	  graphTextFormat (oldT) ;
	}

      if (class(seg->key))  /* we are only interested for now in exon/intron */
	continue ;

      if (seg->x2 < look->min ||
	  seg->x1 > look->max) 
	continue ;
      

      if (!frame)
	{ foundIt = TRUE ; break ; }
     
      /*** invert the 2 clauses to hide the met/stop bar in non open frames ***/
      if ((int) (seg->source/10) != frame)   /* draw only the oRFS of this frame */
	continue ;
      /*************/

      x1 = seg->x1 > look->min ?  seg->x1 : look->min ;
      x2 = seg->x2 + 2 < look->max ?  seg->x2 + 1 : look->max ;

      type = seg->key ;
      subtype = seg->data.i & 0xfff ; 
      ln =  ((unsigned int)seg->data.i) >> 12 ;
      if (ln < 0) ln = 0 ;
      product = seg->parent ;
      
      if (!foundIt)
	{
	  *offset += 1.4 ;
	  x = *offset ;
	  graphColor (TR_COLOR) ;
	  graphRectangle (x +.5, MAP2GRAPH(look->map, look->min - .5),
			  x + 1.1, MAP2GRAPH(look->map, look->max - .5)) ; 
	  graphColor (BLACK) ;
	}
 
      
      y1 = MAP2GRAPH(look->map, x1 - .5) ;
      y2 = MAP2GRAPH(look->map, x2 - .5) ;	
      array(look->boxIndex, box = graphBoxStart (), int) = ii ;
      color = TR_COLOR ; bgColor = TR_COLOR2 ;
      
      if (strstr (name(type), "Exon"))
	{ 
	  /*if (strstr (name(type), "lternative")) 
	    bgColor = TR_COLOR ; */
	  graphRectangle (x +.2, y1, x + 1.4 , y2) ;
	  ln = seg->x2 - seg->x1 + 1 ;
	  graphBubbleInfo (box, product, 0, messprintf("%s %d bp", isUorf ? "uORF" : name(type), ln)) ;
	}
      else if (strstr (name(type), "uORF"))
	{ 
	  graphRectangle (x +.5, y1, x + 1.1 , y2) ;
	  ln = seg->x2 - seg->x1 + 1 ;
	  graphBubbleInfo (box, product, 0, messprintf("%s %d bp", name(type), ln)) ;
	}
      else if (!foundIt && strstr (name(type), "5' UTR"))
	{  
	  bgColor = WHITE ;
	  graphRectangle (x +.5, y1, x + 1.1 , y2) ;
	  ln = seg->x2 - seg->x1 + 1 ;
	  graphBubbleInfo (box, product, 0, messprintf("%s %d bp", name(type), ln)) ;
	}
      else if (!foundIt && strstr (name(type), "uORF"))
	{  
	  bgColor = PALEGREEN ;
	  graphRectangle (x +.5, y1, x + 1.1 , y2) ;
	  ln = seg->x2 - seg->x1 + 1 ;
	  graphBubbleInfo (box, product, 0, messprintf("%s %d bp", name(type), ln)) ;
	}
      else if (strstr (name(type), "3' UTR"))
	{  
	  bgColor = WHITE ;
	  graphRectangle (x +.5, y1, x + 1.1 , y2) ;
	  ln = seg->x2 - seg->x1 + 1 ;
	  graphBubbleInfo (box, product, 0, messprintf("%s %d bp", name(type), ln)) ;
	}
      else if (strstr (name(type), "Intron"))
	{ 
	   if (subtype == _gt_ag || subtype == _gc_ag ||
	       subtype == _Fuzzy_gt_ag || subtype == _Fuzzy_gc_ag) ;
	   else
	     color = BLUE ; 
	   if (strstr (name(type), "lternative"))
	     bgColor = (color == TR_COLOR) ? TR_COLOR2 : LIGHTBLUE ;
	   else
	     bgColor = WHITE ;

	  fmapcdnaDrawTriangle (box, x + .2, (y1+y2)/2.0) ; 
	  if (ln > 0)
	    graphBubbleInfo (box, product, 0, messprintf("%s scar %s %d bp", name(type), name(subtype), ln)) ;
	  else
	    graphBubbleInfo (box, product, 0, messprintf("%s scar %s", name(type), name(subtype))) ;
	}
      else if (type == _ORF_Gap)
	{ 
	  bgColor = YELLOW ;
	  graphRectangle (x +.2, y1, x + 1.4 , y2) ; 
	  ln = seg->x2 - seg->x1 + 1 ;
	  graphBubbleInfo (box, product, 0, messprintf("%s %d bp", name(type), ln)) ;
	}
      else if (type == _Gap)
	{ 
	  float mid = (y1 + y2) / 2 ;
	  graphColor (BLACK) ;
	  graphLine (x +.8, y1, x + .8 , y2) ; 
	  graphRectangle (x + .5, mid - .2, x + 1.1, mid + .2) ;
	}
      else if (type == _Valid3p)
	{ 
	  int n = seg->data.i ;
	  if (n < 0) { color = BLUE ; n = -n ; }
	  else color = BLACK ;
	  graphLine (x - 1.3, y2, x + .5 , y2) ; 
	  if (n > 5)
	    graphFillRectangle (x -1.3, y2 -.5 , x -.8 , y2) ; 
	  else
	    graphRectangle (x -1.3, y2 -.4 , x -.9 , y2) ; 
	  bgColor = TRANSPARENT ;
	  graphBubbleInfo (box, product, 0
			   , messprintf("3' end with poly-A signal supported by %d clones", seg->data.i)
			   ) ;
	}
      else if (type == _Valid5p)
	{ 
	  if ((seg->data.i >> 24) >= 1)
	    color = BLACK ;
	  else
	    color = BLUE ;
	  graphLine (x - 1.3, y1, x + .5 , y1) ; 
	  if ((seg->data.i >> 24) >= 2)
	    graphFillRectangle (x - 1.3, y1 +.5 , x -.8 , y1) ; 
	  else
	    graphRectangle (x - 1.3, y1 +.4 , x -.9 , y1) ; 
	  bgColor = TRANSPARENT ;
	  graphBubbleInfo (box, product, 0
			   , messprintf("5' end supported by %d clones", seg->data.i & 0xffffff)
			   ) ;
	}
      graphBoxEnd () ;
      graphBoxDraw (box, color, bgColor) ;
      fMapBoxInfo (look, box, seg) ;
      graphBoxFreeMenu(box, fMapProductAction, fMapProductMenu) ; 
      
      foundIt = TRUE ;
    }

  if (frame)
    {	
      graphColor (BLACK) ;
      if (foundIt)
	*offset += 2.4 ;
      else
	*offset += 1.4 ;
    }
      
  return foundIt ;
  graphColor (BLACK) ;
}  /* fMapcDNAShowMProductInFrame */

/*********/

void fMapcDNAShowMProduct (LOOK look, float *offset)
{
  int frame ;
  BOOL isDown ;
  int TR_COLOR =  MAGENTA ; /* ORANGE magenta for the web */
  int TR_COLOR2 =  LIGHTMAGENTA ; /* YELLOW palemagenta for the web */
  
cDNAAlignInit () ; 

  if (*look->map->activeColName == '-')
    isDown = FALSE ;
  else
    isDown = TRUE ;
   
  if (fMapcDNAShowMProductInFrame (look, offset, isDown, 0, TR_COLOR, TR_COLOR2)) /* something to do */
    {
      *offset += 2 ;
      for (frame = 1 ; frame < 4 ; frame++)
	fMapcDNAShowStopsInFrame (look, offset, isDown, frame, TR_COLOR) ;
      for (frame = 1 ; frame < 4 ; frame++)
	fMapcDNAShowMProductInFrame (look, offset, isDown, frame, TR_COLOR, TR_COLOR2) ;
    }
} /* fMapcDNAShowMProduct */

/****************************************************/

BOOL fMapcDNAProbePosition (LOOK look, OBJ Mrna, SEG *seg, KEY probe, int *p1, int *p2)
{
  int ii, a1, x1, x2, pos1 = *p1, pos2 = *p2 ;
  Array units = 0 ;
  BSunit *uu ;
  BOOL ok = TRUE ;
  void *v1 ;
  const void *vv ;
  Array bucket = 0 ;
  int iBucket ;

  if (!pos1 || !pos2) 
    return FALSE ;
  units = arrayCreate (500, BSunit) ;
  if (class(look->seqKey) != _VmRNA)
    {
      ok = FALSE ;
      if (bsGetArray (Mrna, _Splicing, units, 5))
	for (ii = 0 ; !ok && ii < arrayMax (units) ; ii += 5)
	  {
	    uu = arrayp (units, ii, BSunit) ;
	    a1 = uu[0].i ;	/* a2 = uu[1].i ; */
	    x1 = uu[2].i ;	x2 = uu[3].i ;
	    if (pos1 >= x1 && pos1 <= x2)
	      {
		*p1 = a1 + pos1 - x1 ;
		*p2 = a1 + pos2 - x1 ;
		pos1 = *p1 ; pos2 = *p2 ;
		ok = TRUE ;
	      }
	  }
    }
  pos1 = seg->type & 0x1 ? seg->x2 - pos1 + 1 : seg->x1 - pos1 + 1 ;
  v1 = assVoid (pos1) ;
  iBucket = 0 ;
  bucket = arrayCreate (32, const void *) ;
  while (ok && assFindNext (look->probeAss, assVoid(probe), &vv, &bucket, &iBucket))
    if (vv == v1)
      ok = FALSE ;
  if (ok)
    assMultipleInsert (look->probeAss, assVoid(probe), assVoid (pos1)) ;

  arrayDestroy (units) ;
  return ok ;
} /* fMapcDNAProbePosition */

/****************************************************/

void fMapcDNAShowProbe (LOOK look, float *offset) 
{ 
  int n, ii, box = 0, ix = 0 ;
  SEG *seg ;
  float x = 0, y1, y2, oldTextHeight ;
  BOOL isDown ;
  KEY gene = 0 ;
  BUMP bump = 0 ; 
  float yz ;
  int color, color2 ;
  int PROBE_COLOR = ORANGE ;       /* RED ; ORANGE DARKGREEN ;  magenta for the web */
  int PROBE_COLOR2 =  PALEORANGE ; /* PALERED ;  YELLOW LIGHTGREEN ;  palemagenta for the web */
  cDNAAlignInit () ; 

  oldTextHeight = graphTextHeight (1.0) ;
  if (*look->map->activeColName == '-')
    isDown = FALSE ;
  else
    isDown = TRUE ;

  for (ii = 1 ; ii < arrayMax(look->segs) ; ++ii)
    { seg = arrp(look->segs,ii,SEG) ;
      switch (seg->type)
	{
	case PROBE: if (isDown) break ; else continue ; 
	case PROBE_UP: if (!isDown) break ; else continue ;
	default: continue ;
	}

      if (class(seg->key) == _VMethod)
	continue ;

      y1 = MAP2GRAPH(look->map,seg->x1 - .5) ;
      y2 = MAP2GRAPH(look->map,seg->x2+1 - .5) ;	/* to cover full base */

      if (y2 < .3 + topMargin ||
	  y1 > mapGraphHeight - .3)
	continue ;

      if (!bump)
	{
	  *offset += .3 ;
	  bump = bumpCreate (mapGraphWidth, 0) ;
	}

      if (y1 < .2 + topMargin) 
	y1 = .2 + topMargin ;
      if (y2 > mapGraphHeight - .2)
	y2 = mapGraphHeight - .2 ;  
      
      yz = y1 + .1 ; ix = 0 ;   /* cheat a little to avoid bump by rounding */
      bumpItem (bump,1,(y2 - y1 - .2) , &ix, &yz) ;  

      x = *offset + 1.7 * ix ;
      gene = seg->key ;
      
      array(look->boxIndex, box = graphBoxStart (), int) = ii ;

      color = PROBE_COLOR ; 
      color2 = seg->data.i & 0x4000 ? PROBE_COLOR2 : PALEGRAY ; 
      graphRectangle (x , y1, x + 1.4, y2) ;
	
      n = strlen(name(gene)) ;
      if (y2 - y1 > n)
	{
	  int oldFormat = graphTextFormat(BOLD) ;
	  graphColor (BLACK) ;
	  graphTextUp (name(gene), x + .25, (y1 + y2)/2 + n/2) ; 
	  graphTextFormat(oldFormat) ;
	}

      graphBoxEnd () ;
      graphBoxDraw (box, color, color2) ;
            
      fMapBoxInfo (look, box, seg) ;
      if (0) graphBoxFreeMenu(box, fMapGenesAction, fMapGenesMenu) ; 
      if (isGifDisplay)
	{
	  graphBubbleInfo (box, seg->key, className(seg->key) 
			   , seg->data.i & 0x4000 ? "exact" : "approximate"
			   ) ;
	}
    }   
 
  if (bumpMax (bump))
    *offset += 1.7 * bumpMax (bump) + .5 ;
  bumpDestroy (bump) ;
  graphTextHeight (oldTextHeight) ;
  graphColor (BLACK) ;
} /* fMapShowProbe */

/***************************************************************************************/
#include "../wabi/boubou.h"

static int fMapcDNATissueColor (const char *tissue, int nn)
{
  int col[] = {RED7, CYAN, GREEN7, BLUE7, BLACK} ;
  
  /*
    PNX *pnx ;
    for (pnx = solexaAll ; tissue && pnx->p ; pnx++)
    if (!strcmp (tissue, pnx->p))
      return pnx->col ;
  */
  return col[nn % 5] ;
} /* fMapcDNATissueColor */

/*****/

void fMapcDNAShowSolexa (LOOK look, float *offset) 
{ 
  int ii, x, s, ntissue = 0 ;
  int sMax, sMax0 = 50 ; /* limit the horizontal "height" of the wiggle */
  float z0 = 0, z1, y1, y0 ;
  OBJ Mrna ;
  Array aa = 0 ;
  BSunit *uu ;
  const char *tissue = 0 ;
  
  if (0 && isGifDisplay) /* on the web */
    return ;
  cDNAAlignInit () ; 

  if (*look->map->activeColName == '-')
    { 
      return ;
    }
  else
    {
      if (!keyFindTag (look->seqKey, str2tag("Solexa")))
	return ;
    }

  Mrna = bsCreate (look->seqKey) ;
  aa = arrayCreate (32000, BSunit) ;
  bsGetArray (Mrna, str2tag("Solexa"), aa, 3) ;

  if (arrayMax (aa))
    {

      y0 = MAP2GRAPH(look->map, look->map->min) ; 
      y1 = MAP2GRAPH(look->map, look->map->max) ; 
      if (y0 < .2 + topMargin) y0 =  .2 + topMargin ; 
      if (y1 > mapGraphHeight - .2) y1 = mapGraphHeight - .2 ;
      graphLine (*offset, y0, *offset, y1) ;
      for (ii = 0, sMax = sMax0 ; ii < arrayMax(aa) ; ii+= 3)
	{ 
	  uu = arrp (aa, ii, BSunit) ;
	  s = uu[2].i ;
	  if (s > sMax) sMax = s ; 
	}
      for (ii = 0 ; ii < arrayMax(aa) ; ii+= 3)
	{
	  uu = arrp (aa, ii, BSunit) ;
	  if (!uu[0].s) continue ;
	  if (!tissue || strcmp (tissue, uu[0].s))
	    {
	      tissue = uu[0].s ;
	      y0 = 0 ; z0 = *offset ;

	      graphColor (fMapcDNATissueColor (tissue, ntissue)) ;

	      ntissue++ ;
	    }
	  x = uu[1].i ; s = uu[2].i * sMax0/sMax ;
	  z1 = *offset + s ;
	  y1 = MAP2GRAPH(look->map, x) ;  
	  if (y1 < .2 + topMargin) continue ;
	  if (y1 > mapGraphHeight - .2) continue ;
	  if (y0 == 0)
	    y0 = MAP2GRAPH(look->map, x - 1) ;  
	  if (z0 == *offset) y0 = y1 ;
	  if (z1 > *offset)
	    graphLine (z0, y0, z1, y1) ;
	  else if (z0 > *offset)
	    graphLine (z0, y0, *offset, y0) ;
	  y0 = y1 ; z0 = z1 ;    
	}
      graphColor (BLACK) ;
    }
  bsDestroy (Mrna) ;
  arrayDestroy (aa) ;
  return ;
} /* fMapcDNAShowSolexa */

/***************************************************************************************/

void fMapcDNAShowTranscribedgene  (LOOK look, float *offset)
{
  int ii, box = 0 ;
  SEG *seg ;
  float x = 0, y1, y2 ;
  BOOL isDown ;
  BUMP bump = 0, bumpIntrons = 0 ; 
  float yz ;
  int ix = 0, color = 0, lastParent = -1, lastIx = -1 ;
  int TG_COLOR = DARKCYAN ; /* DARKGREEN ;  magenta for the web */
  int TG_COLOR2 =  LIGHTCYAN ; /* LIGHTGREEN ;  palemagenta for the web */
  extern int fMap_rView ;

  cDNAAlignInit () ;

  if (*look->map->activeColName == '-')
    { isDown = FALSE ;
    }
  else
    { isDown = TRUE ;
    }

  if (!fmapView (look->view, isDown != fMap_rView ? str2tag ("Fmap_Tg_Down") : str2tag ("Fmap_Tg_Up")))
    return ;

  for (ii = 1 ; ii < arrayMax(look->segs) ; ++ii)
    { seg = arrp(look->segs,ii,SEG) ;
      switch (seg->type)
	{
	case TRANSCRIBEDGENE: if (isDown) break ; else continue ;
	case TRANSCRIBEDGENE_UP: if (!isDown) break ; else continue ;
	default: continue ;
	}

      if (class(seg->key) == _VMethod)
	continue ;

      if (seg->x1 > look->max || seg->x2 < look->min)
	continue ;

      if (seg->key == seg->parent) /* the overall transcribedgene, do not draw */
	{ int dy = 1 ;
	  if (seg->data.i)
	    {
	      SEQINFO *sinf = arrayp (look->seqInfo, seg->data.i, SEQINFO) ;
	      if (sinf->flags & 0x1)
		{
		  if (seg->type == TRANSCRIBEDGENE)
		    graphText ("*", *offset, MAP2GRAPH(look->map,seg->x1) - 1) ;
		  else if (seg->type == TRANSCRIBEDGENE_UP)
		    graphText ("*", *offset, MAP2GRAPH(look->map,seg->x2) + 1) ; 
		}
	      if (sinf->flags & 0x1) dy = 2 ;
	      if (sinf->flags & 2)
		{
		  if (seg->type == TRANSCRIBEDGENE)
		    graphText ("SL1", *offset, MAP2GRAPH(look->map,seg->x1) - dy) ;
		  else if (seg->type == TRANSCRIBEDGENE_UP)
		    graphText ("SL1", *offset, MAP2GRAPH(look->map,seg->x2) + dy) ;
		}
	      if (sinf->flags & 4)
		{
		  if (seg->type == TRANSCRIBEDGENE)
		    graphText ("SL2", *offset, MAP2GRAPH(look->map,seg->x1) - dy) ;
		  else if (seg->type == TRANSCRIBEDGENE_UP)
		    graphText ("SL2", *offset, MAP2GRAPH(look->map,seg->x2) + dy) ;
		}
	      if (sinf->flags & 8)
		{
		  if (seg->type == TRANSCRIBEDGENE)
		    graphText ("SL1&2", *offset, MAP2GRAPH(look->map,seg->x1) - dy) ;
		  else if (seg->type == TRANSCRIBEDGENE_UP)
		    graphText ("SL1&2", *offset, MAP2GRAPH(look->map,seg->x2) + dy) ;
		}
	    }
	  continue ;
	}

      y1 = MAP2GRAPH(look->map,seg->x1 - .5) ;
      y2 = MAP2GRAPH(look->map,seg->x2+1 - .5) ;	/* to cover full base */
      if (y1 < .2 + topMargin) y1 = .2 + topMargin ;
      if (y2 < .2 + topMargin) y2 = .2 + topMargin ;
      if (y2 > mapGraphHeight - .2)	y2 = mapGraphHeight - .2 ;  

      yz = y1 + .03 ; ix = 0 ;  /* cheat a little to avoid bump by rounding */

      if (!bump)
	bump = bumpCreate (mapGraphWidth, 0) ;
      bumpItem (bump,1,(y2 - y1) - .06, &ix, &yz) ;

      if (seg->parent == lastParent && ix < lastIx)
	ix = lastIx ; 
      if (seg->parent != lastParent)
	{ lastIx = ix ; lastParent = seg->parent ; }
      x = *offset + 1.7*ix ;
      array(look->boxIndex, box = graphBoxStart (), int) = ii ;
      graphColor (TG_COLOR) ; 
      if (strstr (name(seg->key), "ntron"))
	{
	  color = TRANSPARENT ;
	  if (seg->data.k && 
	      (KEYKEY(seg->data.k) == _Fuzzy_gt_ag || seg->data.k == _Fuzzy_gc_ag)
	      )
	    { 
	      if (class(seg->data.k) == 1)
		graphColor (RED) ;
	      graphLine (x +.6, y2, x + 1.2 , y2) ; 
	      graphLine (x +.9, y2, x + .9, y1) ;
	    }
	  else if (seg->data.k && 
		   strstr(name(KEYKEY(seg->data.k)), "Fuzzy"))
	    { 
	      graphColor (BLUE) ;
	      graphLine (x +.6, y2, x + 1.2 , y2) ; 
	      graphLine (x +.9, y2, x + .9, y1) ;
	    }
	  else
	    {
	      float mid = (y1 + y2) / 2 ;
	      graphColor (TG_COLOR) ;
	      if (seg->data.k && 
		  KEYKEY(seg->data.k) != _gt_ag &&
		  KEYKEY(seg->data.k) != _gc_ag)
		graphColor (BLUE) ;
	      color = TRANSPARENT ;
	      if (class(seg->data.k) > 0)
		color = GREEN1 + class(seg->data.k) - 1 ;
	      if (0) /* 2023_04_25 obsolete color pattern */
		graphColor (RED) ;
	      if (color > GREEN6 && color <= GREEN8)
		graphColor (WHITE) ;
	      else
		graphColor (BLACK) ;
	      graphLine (x +.8, y2, x + 1.6 , mid) ;
	      graphLine (x +1.6, mid, x + .8, y1 +.1) ;
	    }

	  if (0)
	    {
	      int z, intronX = 0 ; 
	      float intronY = y1 ;
	      
	      if (!bumpIntrons)
		bumpIntrons = bumpCreate (300, 0) ;
	      bumpItem (bumpIntrons,1,(y2 - y1) - .06, &intronX, &intronY) ;
	      color = intronCol [intronX % numIntronCol] ;
	      z = ((171 * (seg->x1)) % 30269) ^ ((172 * seg->x2) % 30307) ; /* some large int */
	      graphAssociate (assVoid(z), assVoid(intronX % numIntronCol)) ;
	    }
	}
      else if (strstr (name(seg->key), "xon") &&
	       strstr (name(seg->key), "artial"))
	{ 
	  graphRectangle (x +.2, y1, x + 1.4 , y2) ;
	  color = WHITE ;
	}
      else if (strstr (name(seg->key), "xon"))
	{ 
	  char *cp = name(seg->parent) ;
	  color = TG_COLOR2 ;
	  graphColor (TG_COLOR) ; 
	  if (0 && *(cp+1) == '_')
	    switch (*cp)
	      {
	      case 'A': /* AceView */
	      case 'B': /*  AceView */
		graphColor (MAGENTA) ; 
		color = LIGHTMAGENTA ;
		break ;
	      case 'J':	
	      case 'G':	/* JBIRC */
		graphColor (CYAN) ; 
		color = PALECYAN ;
		break ;
	      case 'N':  /* NEDO */
		graphColor (BLUE) ; 
		color = LIGHTBLUE ;
		break ;
	      case 'K':	 /* Kent */
		graphColor (ORANGE) ; 
		color = YELLOW ;
		break ;
	      case 'E':	
	      case 'F':  /* EBI */
		graphColor (VIOLET) ; 
		color = PALEVIOLET ;
		break ;
	      case 'H': /* NCBI */
		break ;
	      case 'S':	  /* Splign */
		graphColor (RED) ; 
		color = LIGHTRED ;
		break ;
	      case 'U':	  /* Unigene */
		graphColor (GRAY) ; 
		color = PALEGRAY ;
		break ;
	      default:
		break ;
	      }
	  graphRectangle (x +.2, y1, x + 1.4 , y2) ;
	}
      else if (seg->key == _Gap)
	 {  
	   float mid = (y1 + y2) / 2 ;
	   graphColor (BLACK) ;
	   graphLine (x +.8, y1, x + .8 , y2) ; 
	   graphRectangle (x + .5, mid - .2, x + 1.1, mid + .2) ;
	   color = WHITE ;
	 }
      else
	{
	  graphColor (RED) ;
	  graphRectangle (x +.1, y1, x + .4 , y2) ;
	}
      graphColor (BLACK); 

      graphBoxEnd () ;
      if (seg->key == _Last_exon)
	graphCircle (x + .8, isDown ? y2 : y1 , .6) ;
      
      array(look->boxIndex, box, int) = ii ;
      graphBoxDraw (box, BLACK, color) ;
      /* mhmp 30.09.99 
      graphBoxInfo (box, seg->key, name(seg->key)) ; */
      fMapBoxInfo (look, box, seg) ; /*  */
      graphBoxFreeMenu(box, fMapTranscribedgeneAction, fMapTranscribedgeneMenu) ;
      graphBubbleInfo (box, seg->parent, className(seg->parent),
		       messprintf("%s %s", seg->key ? name(seg->key) : "", seg->data.k ? name(seg->data.k) : "")) ;
    }  
  
  *offset += 1.7*bumpMax (bump) ;
  bumpDestroy (bump) ;
  bumpDestroy (bumpIntrons) ;
  graphColor (BLACK) ;
  /* fMapcDNAShowTranscribedgene2 (look, offset) ; */
}

/*********************************************************************/
/*********************************************************************/


Array fMapcDNAReferenceDna = 0,  fMapcDNAReferenceHits = 0 ;

#ifdef ACEMBLY

KEY  fMapcDNAReferenceGene = 0 ;
KEY  fMapcDNAReferenceEst = 0 ;
int fMapcDNAReferenceOrigin = 0 ;

static void annotStart (KEY gene)
{ 
  KEY key ;
  OBJ Gene = bsUpdate(gene) ;
  lexaddkey(name(gene), &key, _VAnnotation) ;
  
  sessionAutoSave (60 * 30, 60 * 5) ;   /* save every 30 minutes, or after 5 minutes inactivity */

  bsAddKey (Gene, str2tag ("Annotations"), key) ;
  bsSave (Gene) ;
  autoAnnotateGene (gene, key) ;
  annotate (key, notes) ;
}

/*********************************************************************/

BOOL fMapcDNAPickTrace (KEY gene, KEY est, int from)
{
  KEY key = 0 ; 
  SEG *seg2  ;
  Array aa = 0 ;
  int i, j, origin = 0 ;
  char *cp, *cq ;
  Array hits = 0 ;
  FMAPLOOKGET("fMapcDNAPickTrace") ;

  
  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
      {	
	seg2 = arrp(look->segs,i, SEG) ;  
	if (seg2->key != gene || seg2->key != seg2->parent)  /* I want the overall transcribedgene */
	  continue ;
	if (seg2->type & 0x1) /* wrong strand */
	  continue ;
	origin = seg2->x1 > 5000 ? 5000 : seg2->x1 ;
	j = seg2->x2 + origin + 5000 ; /* take gene plus 5 kb downstream */
	if (j > arrayMax(look->dna))
	  j = arrayMax(look->dna) ;
	if (j > seg2->x1 - origin)
	  j = j -  seg2->x1 + origin ;
	if (j < 0) /* happens if zoom is incorrectly set */
	  continue ; 
	aa = arrayCreate (j, char) ;
	array(aa,j - 1,char) = 0 ; 
	arrayMax(aa)-- ; /* make room, add terminal 0 */
	cp = arrp(aa,0,char) ; cq = arrp(look->dna, seg2->x1 - origin,char) ;
	while(j--) *cp++ = *cq++ ;
	hits = cDNAGetReferenceHits (gene, origin) ;
	if (!hits)
	  { arrayDestroy(hits) ; arrayDestroy(aa) ;
	  continue ;
	  }
	fMapcDNAReferenceOrigin = origin ;
	fMapcDNAReferenceDna = aa ;
	fMapcDNAReferenceHits = hits ;
	fMapcDNAReferenceGene = gene ;
	fMapcDNAReferenceEst = est ;
	key = gene ;
	origin = seg2->x1 ;
      }
  
  if (key)
    {
      if (look->origin != origin)
	{ look->origin = origin ;
	fMapDraw (look, 0) ;
	}	    
      display (key, from, DtMultiTrace) ; 
      /*      annotStart (gene) ; */
      return TRUE ;
    } 
  return FALSE ;
}
#endif

BOOL fMapcDNAFollow (int box)
{
  KEY clone, est ;
  SEG *seg ;  
#ifdef ACEMBLY
  SEG *seg2 ;  
  int from, i ;
  KEY key ;
  KEYSET geneSet = 0 ;
#endif  
  FMAPLOOKGET("fMapcDNAFollow") ;
  
  seg = BOXSEG(box) ;
  if (!seg)
    return FALSE ; 
  
  if (isDisplayBlocked())
    { 
      if ((seg->type | 0x1) == TRANSCRIBEDGENE_UP)
	display (seg->parent, 0, 0) ;
      else
	display (seg->key, 0, 0) ;
      return TRUE ;
    }

  if ((seg->type | 0x1) == TRANSCRIBEDGENE_UP)
    { 
#ifdef ACEMBLYJUNK
      KEY gbAnnotKey, annot = 0 ;

      if ((annot = keyGetKey (seg->parent, str2tag ("Annotations"))) &&
	  (gbAnnotKey =  keyGetKey (annot, str2tag ("Gb_annot"))))
	display (gbAnnotKey, 0, 0) ;
      geneAnnotDisplay (0, seg->parent, 0) ;
#else
      display (seg->parent, look->seqKey, TREE) ;
#endif
      return TRUE ;
    }
  else if ((seg->type | 0x1) != SPLICED_cDNA_UP)
    return FALSE ;

  est = seg->key ;
  fMapcDNAReferenceDna = 0 ;

#ifdef ACEMBLY

  from = GRAPH2MAP (look->map, graphEventY + .2) ; /* + look->map->mag * .5) ; */
  sessionAutoSave (60 * 30, 60 * 5) ;   /* save every 30 minutes, or after 5 minutes inactivity */

  geneSet = queryKey (est, "FOLLOW From_gene") ;
  if (keySetMax(geneSet))
    for (i = key = 0 ; !key && i < arrayMax(look->segs) ; ++i)
      {
	seg2 = arrp(look->segs,i, SEG) ;
	if (seg2->type == TRANSCRIBEDGENE && keySetFind (geneSet, seg2->key, 0) && 
	    seg2->x2 > seg2->x1 + 10 &&
	    from >= seg2->x1 && from < seg2->x2 &&
	    seg2->x1 >=0 && seg2->x2 <= arrayMax(look->dna) )
	  {
	    if (seg2->type & 0x1) /* wrong strand */
	      continue ; 
	    if (seg2->key != seg2->parent)  /* I want the overall transcribedgene */
	      continue ; 
	    if (fMapcDNAPickTrace (seg2->key, est, from - seg2->x1))
	      { keySetDestroy (geneSet) ; return TRUE ; }
	  }
      }
  keySetDestroy (geneSet) ;
#endif

  clone = keyGetKey (est, _cDNA_clone) ;
  display (clone ? clone : est, 0, TREE) ;
  return TRUE ;
}


/***************************************************************************************/
/***********************************************************************/
#ifdef ACEMBLY
static KEY fMapTranscribedGeneFuseToCloneClone = 0 ;
static void fMapTranscribedGeneFuseToClone (KEY clone)
{
  OBJ Clone = 0 ;
  FMAPLOOKGET("fMapTranscribedGeneFuseToClone") ;

  if (!(class(clone) == _VcDNA_clone))
    {
      messout ("Please pick a cDNA clone") ;
      return ;
    }
  if ((Clone = bsUpdate (fMapTranscribedGeneFuseToCloneClone)))
    {
      bsAddKey (Clone, _Fuse_to_clone, clone) ;
      bsSave (Clone) ;
    }
}
#endif
/***************************************************************************************/

static void fMapcDNAAction (KEY key, int box)
{ 
  SEG *seg ;
  KEY clone = 0, dna = 0 ;
#ifdef ACEMBLY
  static mytime_t lastStrandDate = 0 ;
  OBJ obj = 0 ;
  KEY gene = 0;
  KEYSET genes = 0 ;
  int i ;
#endif
  FMAPLOOKGET("fMapcDNAAction") ;

  seg = BOXSEG(box) ;
  if (!seg)
    return ;
  displayPreserve() ;
  switch (key)
    {
    case 'c':
      clone = keyGetKey (seg->key, _cDNA_clone) ;
      display (clone ? clone : seg->key, 0, TREE) ;
      break ;
    case 'd': 
      dna = keyGetKey (seg->key, _DNA) ;
      if (dna) display (dna, 0, 0) ;
      else display (seg->key, 0, TREE) ;
      break ;
    case 't':
      display (seg->key, 0, TREE) ;
      break ;
#ifdef ACEMBLY
    case 's':
      if (keyFindTag (seg->key, _SCF_File))
	display (seg->key, 50, DtMultiTrace) ;
      else display (seg->key, 0, TREE) ;
      break ;
    case 'w':
      if (!sessionGainWriteAccess())
	return ; 
      sessionAutoSave (60 * 30, 60 * 5) ;   /* save every 30 minutes, or after 5 minutes inactivity */

      if (lastStrandDate && 
	  (timeNow () - lastStrandDate > 120) &&
	  !messQuery("Do you really believe that the strand of this est is wrong"))
	return ;
      lastStrandDate = timeNow () ;
      if ((obj = bsUpdate(seg->key)))
	{
	  int i ;
	  OBJ obj1 = 0 ;
	  BOOL isF = FALSE ;
	  KEYSET tg2, tg3, tg1 = queryKey (seg->key, ">from_gene") ;

	  clone = keyGetKey(seg->key, _cDNA_clone) ;
	  if (bsFindTag (obj, _Forward))
	    bsAddTag (obj, _Reverse) ;
	  else if (bsFindTag (obj, _Reverse))
	    { isF = TRUE ; bsAddTag (obj, _Forward) ; }
	  if (bsFindTag (obj, _Intron_boundaries))
	    bsRemove (obj) ;

	  if (isF)
	    {
	      bsAddTag (obj, _mForward) ;
	      bsAddTag (obj, _Colour) ;
	      bsPushObj (obj) ;
	      if (0) bsAddTag (obj, _YELLOW) ;
	    }
	  else
	    {
	      bsAddTag (obj, _mReverse) ;
	      bsAddTag (obj, _Colour) ;
	      bsPushObj (obj) ;
	      if (0) bsAddTag (obj, _ORANGE) ;
	    }
	  bsSave (obj) ;
	  look->pleaseRecompute = TRUE ;
	  {
	    KEYSET ks = query (0, "Find mrna !from_gene && !from_prediction") ;
	    keySetKill (ks) ;
	    keySetDestroy (ks) ;
	  }
	  tg2 = queryKey (seg->key, ">from_gene") ;
	  tg3 = keySetOR (tg1, tg2) ;
	  keySetDestroy (tg1) ;
	  keySetDestroy (tg2) ;
	  tg1 = query (tg3, ">read ; !Ref_seq && !( forward && areverse || reverse && aforward); COLOUR = paleyellow") ;
	  tg2 = query (tg3, ">read ; !Ref_seq && (forward && areverse || reverse && aforward); ! COLOUR = paleyellow") ;
	  if (0)
	    {
	      for (i = 0 ; i < keySetMax (tg1) ; i++)
		{
		  if ((obj1 = bsUpdate (keySet(tg1, i))))
		    {
		      if (bsFindTag (obj1, _Colour))
			bsRemove (obj1) ;
		      bsSave (obj1) ;
		    }
		}
	      for (i = 0 ; i < keySetMax (tg2) ; i++)
		{
		  if ((obj1 = bsUpdate (keySet(tg2, i))))
		    {
		      bsAddKey (obj1, _Colour, _PALEYELLOW) ;
		      bsSave (obj1) ;
		    }
		}
	    }
	  keySetDestroy (tg1) ;
	  keySetDestroy (tg2) ;
	  keySetDestroy (tg3) ;
	  fMapDraw(look,0) ;
	}
      break ;
#ifdef ACEMBLY
    case 'e': case 'f': case 'g': case 'm': case 'o': case 'L': case 'R': case 'u': case 'v': case '5': case '3': case 'C':
      clone = keyGetKey(seg->key, _cDNA_clone) ;
      if (!clone || !sessionGainWriteAccess())
	return ;
      sessionAutoSave (60 * 30, 60 * 5) ;   /* save every 30 minutes, or after 5 minutes inactivity */

      genes = queryKey (seg->key, "> From_gene") ;
      for (i = 0 ; i < keySetMax(genes) ; i++)
	{
	  SEG *s = arrp (look->segs, 0, SEG) - 1 ;
	  int j = arrayMax(look->segs) ;
	  gene = keySet (genes, i) ;
	  while (s++, j--)
	    if (s->key == gene &&
		s->x1 < seg->x2 &&  /* intersect non void */
		s->x2 > seg->x1)
	      break ;
	  if (j >= 0)
	    break ; /* relevant gene was not found */
	  gene = 0 ;
	}
      if (1)
	{
	  switch (key)
	    {
	    case 'e':
	      if (gene)
		{
		  if (!messQuery("Do you really want to remove this cDNA clone from this gene"))
		    return ;
		  cDNARemoveCloneFromGene (clone,gene) ;
		}
	      else
		messout ("please remove from the genove view, they may also belong to another mRNA") ;
	      break ;
	    case 'm':
	      if (!messQuery(
"Tag this clone as Mosaic"))
		return ;
	      {
		OBJ Clone = bsUpdate (clone) ;
		if (Clone)
		  {
		    bsAddTag (Clone, _Mosaic) ;
		    bsSave (Clone) ;
		  }
	      }
	      break ;
	    case 'f':
	      if (!messQuery(
"Tag this read at this position as Suspected_internal_deletion, it will be shown but not use in transcripts"))
		return ;
	      {
		OBJ Clone = bsUpdate (clone) ;
		OBJ Gene = 0 ;
		KEY ff ;
		if (Clone)
		  {
		    int n = 0, nest = 0 ;

		    bsAddKey (Clone, _Suspected_internal_deletion, seg->key) ;
		    if (gene && (Gene = bsCreate(gene)))
		      {
			int i, i0 = 0, j0 = 1, loopMax = 3 ;
			Array units = arrayCreate (60, BSunit) ;
			
			while (loopMax-- > 0 && !nest && i0 < j0)
			  /* i must loop because a single clone may have several errors */
			  {
			    n = 0 ;
			    if (!n && bsGetArray (Gene, _Intron_boundaries, units, 6))
			      for (i = i0 ; i < arrayMax(units) ; i+= 6, i0 += 6)
				{ 
				  ff = arr (units, i + 0, BSunit).k ;
				  if (!ff || ff != _Other)
				    continue ;
				  j0 = arrayMax(units) ;
				  if (arr (units, i + 4, BSunit).k == clone)
				    { n = arr (units, i + 2, BSunit).i ; break ; }
				  if (arr (units, i + 5, BSunit).k == clone)
				    { n = arr (units, i + 3, BSunit).i ; break ; }
				}   
			    else
			      j0 = 100 ; /* break the loop */
			    if (n && /* n is the coordinate in the gene, we have to match it to
				      * a coordinate in the est 
				      */
				bsGetArray (Gene, _Assembled_from, units, 5))
			      for (i = 0 ; i < arrayMax(units) ; i+= 5)
				{
				  if (arr (units, i + 2, BSunit).k == seg->key &&
				      arr (units, i + 1, BSunit).i >= n - 5 &&
				      arr (units, i + 1, BSunit).i <= n + 5)
				    { nest = arr (units, i + 4, BSunit).i ; break ; }
				}
			  }
			arrayDestroy (units) ;
			bsDestroy (Gene) ;
		      }
		    if (nest) 
		      bsAddData (Clone, _bsRight, _Int, &nest) ;
		    bsAddKey (Clone, _Manual_internal_deletion, seg->key) ;
		    if (nest) 
		      bsAddData (Clone, _bsRight, _Int, &nest) ;
		    bsSave (Clone) ;
		  }	
	      }
	      break ;
	    case '5':
	      if (!messQuery("Flag this read as a real 5' end"))
		return ;
	      {
		OBJ Est = bsUpdate (seg->key) ;
		if (Est)
		  {
		    bsAddTag (Est, _Real_5prime) ;
		    bsSave (Est) ;
		  }	
	      }
	      break ;
	    case '3':
	      if (!messQuery("Flag this est as a real 3' end"))
		return ;
	      {
		OBJ Est = bsUpdate (seg->key) ;
		if (Est)
		  {
		    bsAddTag (Est, _Real_3prime) ;
		    bsSave (Est) ;
		  }	
	      }
	      break ;
	    case 'u':
	      if (!messPrompt(
			      "Tag this read as Duplicate_clone, it will be shown but not use in transcripts","","t"))
		return ;
	      {
		OBJ Clone = bsUpdate (clone) ;
		if (Clone)
		  {
		    char *comment = strnew (freeword (),0) ;
		    KEY cc = 0 ;

		    bsAddTag (Clone, _Duplicate_clone) ;
		    if (comment && *comment) 
		      {
			lexaddkey (comment, &cc, _VText) ;
			bsAddKey (Clone, _Comment, cc) ;
		      }
		    bsSave (Clone) ;
		  }	
	      }
	      break ;
	    case 'o': /* fuse to clone */
	      fMapTranscribedGeneFuseToCloneClone = clone ;
	      displayBlock (fMapTranscribedGeneFuseToClone, "Select a clone you want to fuse to") ;
	      break ;
	    case 'R': /* restrict to active zone */
	      if (!messPrompt(
			      "Restrict the EST to be mapped to the active zone","","t"))
		return ;
	      {
		KEY map = 0 ;
		int a1, a2 ;
		OBJ Est = bsUpdate (seg->key) ;
		if (map && Est)
		  {
		    bsAddKey (Est, _IntMap, map) ;
		    bsAddData (Est, _bsRight, _Int, &a1) ;
		    bsAddData (Est, _bsRight, _Int, &a2) ;
		    bsSave (Est) ;
		  }	
	      }
	      break ;
	    case 'L': /* split as 1 clone per read */
	      {
		KEY newClone = 0 ;
		KEYSET reads = queryKey (clone, ">Read Is_read") ;
		int i, n, nn = keySetMax (reads) ;
		char *cp, *cq, *cr, *cg ;

		if (nn < 2)
		  messout ("This clone contains a single read, I cannot split it, sorry") ;
		else if (!checkWriteAccess()) ; 
		else if (! gene)
		  messout ("Sorry, I do not know which gene this clone belongs to") ;
		else if (messQuery ("Create a separate clone for each read contained in this clone"))
		  {
		    vTXT txt = vtxtCreate () ; 
		    cg = ac_protect (name (gene), 0) ;
		    for (n = 1 ; n < nn ; n++)
		      {
			for (i = 1 ; i < 20 ; i++)
			  {
			    cp = messprintf ("%s%c", name(clone), (char)('a'+i)) ;
			    if (! lexword2key (cp, &newClone, _VcDNA_clone))
			      break ;
			  }
			if (i >= 20)
			  break ;
			cq = ac_protect (cp, 0) ;
			cr = ac_protect (name (keySet (reads, n)), 0) ;
			lexaddkey (cp, &newClone, _VcDNA_clone) ;
			vtxtPrintf (txt, "cDNA_clone %s\nFrom_gene %s\nRead %s\n\n", cq, cg, cr) ;
			ac_free (cq) ;
			ac_free (cr) ;
		      }
		    parseBuffer (vtxtPtr (txt), 0) ;
		    ac_free (txt) ;
		    ac_free (cg) ;
		  }
		keySetDestroy (reads) ;
		if (0) messout ("Please recompute your gene %s", name(gene)) ;
	      }
	      break ;
	    case 'C': /* add a comment in the clone */
	      if (messPrompt ("Add a comment to the clone", "", "t"))
		{
		  char *comment = strnew (freeword (),0) ;
		  
		  if (comment && *comment) 
		    {
		      OBJ Clone = bsUpdate (clone) ;
		      if (Clone)
			{
			  KEY cc = 0 ;
			  lexaddkey (comment, &cc, _VText) ;
			  bsAddKey (Clone, _Comment, cc) ;
			}
		      bsSave (Clone) ;
		    }	
		}
	      break ;
	    case 'v':
	      if (!messPrompt(
			      "Ignore this clone, it will be shown but not use in transcripts","","t"))
		return ;
	      {
		OBJ Clone = bsUpdate (clone) ;
		if (Clone)
		  {
		    char *comment = strnew (freeword (),0) ;
		    KEY cc = 0 ;

		    bsAddTag (Clone, _Ignore_this_clone) ;
		    if (comment && *comment) 
		      {
			lexaddkey (comment, &cc, _VText) ;
			bsAddKey (Clone, _Comment, cc) ;
		      }
		    bsSave (Clone) ;
		  }	
	      }
	      break ;
	    case 'g':
	      if (!messQuery("Do you really want to trash this read, it will rm it from all genes"))
		return ;
	      {
		OBJ Est = bsUpdate (seg->key) ;
		KEYSET
		  ks = queryKey (seg->key, 
      "{>DNA} $| {>BaseCall} $| { >SCF_Position } $| {>cDNA_clone ; COUNT read = 1}") ;
		if (Est)
		  {
		    if (bsFindTag (Est, _Is_read))
		      bsRemove (Est) ;
		    if (bsFindTag (Est, _cDNA_clone))
		      bsRemove (Est) ;
		    bsAddTag (Est, _Bad_quality) ;
		    bsAddData (Est, _bsRight, _Text, "Trashed manually") ;
		    bsSave (Est) ;
		  }	
		keySetKill (ks) ;
		keySetDestroy (ks) ;
	      }
	    break ;
	    }
	  
	  look->pleaseRecompute = TRUE ;
	  fMapDraw(look,0) ;
	}
      else
	messout ("removal failed, please club jean") ;
      break ;

    case '6':
      sessionAutoSave (60 * 30, 60 * 5) ;   /* save every 30 minutes, or after 5 minutes inactivity */
      clone = keyGetKey(seg->key, _cDNA_clone) ;
      if (clone && (obj = bsUpdate (clone)))
	{
	  bsAddTag (obj, _r_gap) ;
	  bsSave (obj) ;
	}
      break ;
    case '7':
      sessionAutoSave (60 * 30, 60 * 5) ;   /* save every 30 minutes, or after 5 minutes inactivity */
      clone = keyGetKey(seg->key, _cDNA_clone) ;
      {
	KEYSET ests = queryKey (clone, ">Read") ;
	int i = keySetMax(ests) ;
	KEY est ;

	while (i--)
	  if ((est = keySet (ests, i))&& (obj = bsUpdate (est)))
	    {
	      bsAddTag (obj, _Fake_internal_poly_A) ;
	      bsSave (obj) ;
	    }
	if ((obj = bsUpdate (clone)))
	  {
	    bsAddTag (obj, _Internal_priming) ;
	    bsAddTag (obj, _Internal_priming_manual) ;
	    bsSave (obj) ;
	  }
      }
      break ;
    case '8':
      sessionAutoSave (60 * 30, 60 * 5) ;   /* save every 30 minutes, or after 5 minutes inactivity */
      clone = keyGetKey(seg->key, _cDNA_clone) ;
      {
	KEYSET ests = queryKey (clone, ">Read") ;
	int i = keySetMax(ests) ;
	KEY est ;

	while (FALSE && i--)
	  if ((est = keySet (ests, i))&& (obj = bsUpdate (est)))
	    {
	      if (bsFindTag (obj, _Transpliced_to))
		bsRemove (obj) ;
	      bsSave (obj) ;
	    }	
	if ((obj = bsUpdate (clone)))
	  {
	    bsAddTag (obj, _Internal_capping) ;
	    bsSave (obj) ;
	  }
      }
      break ;
    case '9':
      sessionAutoSave (60 * 30, 60 * 5) ;   /* save every 30 minutes, or after 5 minutes inactivity */
      clone = keyGetKey(seg->key, _cDNA_clone) ;
      {
	if ((obj = bsUpdate (clone)))
	  {
	    bsAddTag (obj, _Length_anomaly) ;
	    bsSave (obj) ;
	  }
      }
      break ;
    case 'A':
      sessionAutoSave (60 * 30, 60 * 5) ;   /* save every 30 minutes, or after 5 minutes inactivity */
      clone = keyGetKey(seg->key, _cDNA_clone) ;
      {
	if ((obj = bsUpdate (clone)))
	  {
	    bsAddTag (obj, str2tag ("Possibly_unspliced")) ;
	    bsSave (obj) ;
	  }
      }
      break ;
    case 'B':
      sessionAutoSave (60 * 30, 60 * 5) ;   /* save every 30 minutes, or after 5 minutes inactivity */
      clone = keyGetKey(seg->key, _cDNA_clone) ;
      {
	if ((obj = bsUpdate (clone)))
	  {
	    bsAddTag (obj, str2tag ("Unspliced_non_coding")) ;
	    bsSave (obj) ;
	  }
      }
      break ;
#endif
    case 'p':   /* dotter */
      {
	Array dnaGene, dnaEst = dnaGet (seg->key) ;
	static char *ss1, *ss2 ;
	KEY est = seg->key ;
	KEYSET geneSet = 0 ;
	int i, origin, j, from ;
	char *cp, *cq ;
	SEG *seg2 ;

	est = seg->key ;
	key = gene = 0 ; fMapcDNAReferenceDna = 0 ;
	
	origin = look->origin  ;
	from = (seg->x1 + seg->x2) / 2 ;
	
	geneSet = queryKey (est, "FOLLOW From_gene") ;
	if (keySetMax(geneSet))
	  for (i = 0 ; !key && i < arrayMax(look->segs) ; ++i)
	    {
	      seg2 = arrp(look->segs,i, SEG) ;
	      if (seg2->type == TRANSCRIBEDGENE && keySetFind (geneSet, seg2->key, 0) && 
		  seg2->x2 > seg2->x1 + 10 &&
		  from >= seg2->x1 && from < seg2->x2 &&
		  seg2->x1 >=0 && seg2->x2 <= arrayMax(look->dna) )
		{
		  if (seg2->type & 0x1) /* wrong strand */
		    continue ; 
		  gene = seg2->key ;
		  origin = seg->x1 > 5000 ? 5000 : seg->x1 ;
		  j = seg2->x2 + origin + 5000 ; /* take gene plus 5 kb downstream */
		  if (j > arrayMax(look->dna))
		    j = arrayMax(look->dna) ;
		  j = j -  seg2->x1 + origin ;
		  dnaGene = arrayCreate (j, char) ;
		  array(dnaGene,j - 1,char) = 0 ; 
		  arrayMax(dnaGene)-- ; /* make room, add terminal 0 */
		  cp = arrp(dnaGene,0,char) ; cq = arrp(look->dna, seg2->x1 - origin,char) ;
		  while(j--) *cp++ = *cq++ ;

		  from = from - seg2->x1 ;
		  origin = seg2->x1 ;

		  dnaEst = dnaGet (est) ;
		  if (!dnaEst || ! arrayMax(dnaEst))
		    goto abort ;
		  dnaDecodeArray (dnaEst) ;
		  dnaDecodeArray (dnaGene) ;
		  ss1 = messalloc (arrayMax(dnaEst) + 1) ;
		  memcpy (ss1, arrp(dnaEst, 0, char), arrayMax(dnaEst)) ;
		  ss2 = messalloc (arrayMax(dnaGene) + 1) ;
		  memcpy (ss2, arrp(dnaGene, 0, char), arrayMax(dnaGene)) ;
		  ss1[arrayMax(dnaEst)] = 0 ;
		  ss2[arrayMax(dnaGene)] = 0 ;
		  
		  dotter ('N', "     ", 
			  
			  name(gene), ss2, origin,  /* name, dna and offset */
			  name(seg->key), ss1, 0,
			  
			  
			  from, 0, /* central coord */
			  NULL, NULL, NULL, NULL, 0.0, /* irrelevant */
			  0, NULL, /*DBhits.next */ 0, /* qoffset */  0, 0) ;
		  
		abort:
		  arrayDestroy (dnaEst) ;
		  arrayDestroy (dnaGene) ;
		  
		  break ;
		}
	    }
	keySetDestroy (geneSet) ;
      }
      break ;
#endif
    }
}

/***********************************************************************/

static void fMapProductAction (KEY key, int box)
{ 
  KEY product = 0, gene = 0, mrna = 0 ;
  SEG *seg ;
  Graph old = graphActive () ;
  FMAPLOOKGET("fMapProductAction") ;

  seg = BOXSEG(box) ;
  if (!seg)
    return ;
  product = seg->parent ;

  switch (key)
    {
    case 'w':
      display (product, 0, TREE) ;
      graphActivate (old) ;
      return ;
    case 'g':
      gene = keyGetKey (product, _From_gene) ;
      if (gene)
	{
	  displayPreserve () ;
	  display (gene, product, 0) ;
	}
      return ;
    case 'm':
      mrna = keyGetKey (product, _mRNA) ;
      if (mrna)
	{
	  display (mrna, product, TREE) ;
	}
      return ;
#ifdef ACEMBLY
    case 'U': /* Recompute title */ 
      if (!checkWriteAccess()) 
	return ;
      mrna = keyGetKey (product, _mRNA) ;
      if (mrna)
	mrnaAddKantorInfo (mrna) ;
      break ;
    case 'f': /* Fiche-display  */
      {
	KEY mRNA = keyGetKey (product, _mRNA) ;
	
	if (mRNA)
	  cFicheDisplay (mRNA, 0, FALSE) ;
      }
      return ;
/*Commented_by_vahan 
    case 'F': 
      ficheGraph(product);
      return ;
*/
#endif           
    case 't':
      {
	int option = 2 ;  /* toggle */
	
	if (look->translateGene != seg->parent)
	  option = 1 ;     /* true */
	
	look->translateGene = seg->parent ;
	if (seg->type & 0x1)
	  {
	    mapColSetByName ("Up Gene Translation", option) ;
	    mapColSetByName ("Down Gene Translation", FALSE) ;
	  }
	else
	  {
	    mapColSetByName ("Down Gene Translation", option) ;
	    mapColSetByName ("Up Gene Translation", FALSE) ;
	  }
      }
      break ;
    case 'T': 
      mapColSetByName ("3 Frame Translation", 2) ;
      break ;
    } 
  look->pleaseRecompute = TRUE ;
  graphActivate (old) ;
  fMapDraw(look,0) ;
}

/***********************************************************************/

static void fMapGenesAction (KEY key, int box)
{ 
  KEY gene = 0 ;
  SEG *seg ;  
  BOOL resizeAsGenefinder = FALSE ;
  BOOL keepPosName = FALSE ;
  Graph old = graphActive () ;
  FMAPLOOKGET("fMapgenesAction") ;

  seg = BOXSEG(box) ;
  if (!seg)
    return ;

  gene = seg->key ;
  switch (key)
    {
    case 'w':
      display (gene, 0, TREE) ;
      graphActivate (old) ;
      return ;
      break ;
    case 'R': 
      keepPosName = TRUE ;
      /* fall thru */
    case 'r': /* rename */
      if (!checkWriteAccess()) 
	break ;
      if (messPrompt ("Please propose a new name:", name(gene), "w"))
	{
	  char *cp = strnew (freeword(), 0) ;
	  KEY dummy = 0, dummy2 ;
	  KEY seq = 0, newname = 0 ;
	  int x1 = 0, x2 = 0 ;
	  BOOL gotSeq = FALSE ;
	  OBJ Seq = 0, Gene = 0 ;

	  if ((seq = keyGetKey (gene, _Genomic_sequence)) &&
	      (Seq = bsCreate (seq)) &&
	      bsFindKey (Seq, _Genes, gene) &&
	      bsGetData (Seq, _bsRight, _Int, &x1) &&
	      bsGetData (Seq, _bsRight, _Int, &x2))
	    gotSeq =TRUE ;
	  bsDestroy (Seq) ;
	  
	  if (gene)
	    newname = keyGetKey (gene, _NewName) ;
	  if (lexword2key (cp, &dummy, _VGene) &&
	      keyGetKey (dummy, _SMAP))
	    messout ("Sorry, this gene is already mapped elsewhere") ;
	  else if (lexword2key (cp, &dummy, _VGene) &&
		   (dummy2=keyGetKey (dummy, _NewName)) &&
		   newname && dummy2 != newname)
	    messout ("Sorry, this gene already has a different new name") ;
	  else
	    {
	      if (gotSeq &&
		  (Seq = bsUpdate (seq)))
		{
		  if (bsFindKey (Seq, _Genes, gene))
		    bsRemove (Seq) ;
		  bsSave (Seq) ;
		}

	      lexAlias (&gene, cp, FALSE, FALSE) ;

	      if ((Gene = bsUpdate (gene)))
		{
		  lexaddkey (cp, &dummy, _VNewName) ;
		  if (!keepPosName)
		    bsAddKey (Gene, _NewName, dummy) ;
		  bsSave (Gene) ;
		}
	      
	      if (gotSeq &&
		  (Seq = bsUpdate (seq)))
		{
		  bsAddKey (Seq, _Genes, gene) ;
		  bsAddData (Seq, _bsRight, _Int, &x1) ;
		  bsAddData (Seq, _bsRight, _Int, &x2) ;
		  bsSave (Seq) ;
		}
	      look->pleaseRecompute = TRUE ;
	    }
	}
      break ;
    case 't': 
      resizeAsGenefinder = TRUE ;
      /* fall thru to case 's' */
    case 's': /* resize */
      if (!checkWriteAccess()) 
	break ;
      {
	KEY tg = 0, map = 0, cosmid = 0 ;
	int i, ix1, ix2, ia1, ia2 ;
	OBJ Tg = 0, Gene = 0, Seq = 0 ;
	KEYSET ks, maps, cosmids ;

	if (resizeAsGenefinder)
	  {
	    ks = queryKey (gene, ">genefinder") ;
	    maps = query (ks, ">IntMap") ;
	    cosmids = query (ks, ">Source") ;
	  }
	else /* as transcribed_gene */
	  {
	    ks = queryKey (gene, ">Transcribed_gene") ;
	    maps = query (ks, ">IntMap") ;
	    cosmids = query (ks, ">Genomic_sequence") ;
	  }

	if (keySetMax (ks) == 0)
	  {	
	    if (resizeAsGenefinder)
	      messout ("This gene is not linked to any genefinder model"
		       "Sorry, i cannot proceed") ;
	    else
	      messout ("This gene is not linked to any transcribed genes"
		       "Sorry, i cannot proceed") ;
	  }
	else if (keySetMax (ks) > 1)
	  {	
	    if (resizeAsGenefinder)
	      messout ("This gene is linked to several gene genefinder models"
		       "Sorry, i cannot proceed") ;
	    else
	      messout ("This gene is linked to several transcribed genes"
		       "Sorry, i cannot proceed") ;
	  }
	else if (keySetMax (maps) != 1)
	  {	
	    if (resizeAsGenefinder)
	      messout ("This gene is linkded to a genefinder which does not yet have an IntMap tag"
		       "Sorry, i cannot proceed") ;
	    else
	      messout ("This gene is linkded to a transcribed_gene which does not yet have an IntMap tag"
		       "Sorry, i cannot proceed") ;
	  }
	else if (keySetMax (cosmids) != 1)
	  {
	    if (resizeAsGenefinder)
	      messout ("This gene is linkded to a genefinder which does not have a Source tag"
		     "Sorry, i cannot proceed") ;
	    else
	      messout ("This gene is linkded to a genefinder which does not have a genomic_sequence tag"
		     "Sorry, i cannot proceed") ;
	  }
	else
	  {
	    for (i = 0 ; i < keySetMax (ks) ; i++)
	      {
		tg = keySet (ks, i) ;
		if ((Tg = bsCreate (tg))) /* ok for tg and for pg */
		  {
		    if (bsGetKey (Tg, _IntMap, &map) &&
			bsGetData (Tg, _bsRight, _Int, &ix1) &&
			bsGetData (Tg, _bsRight, _Int, &ix2))
		      { 
			if (ix1 < ix2)
			  {
			    if (!i || ix1 < ia1) ia1 = ix1 ;
			    if (!i || ix2 > ia2) ia2 = ix2 ;
			  }
			else
			  {
			    if (!i || ix1 > ia1) ia1 = ix1 ;
			    if (!i || ix2 < ia2) ia2 = ix2 ;
			  }
		      }
		    bsDestroy (Tg) ;
		  }
	      }
	    map = keySet (maps, 0) ;
	    if ((Gene = bsUpdate (gene))) /* ok for tg and for pg */
	      {
		if (bsFindTag (Gene, _Genomic_sequence))
		  bsRemove(Gene) ;
		if (bsFindTag (Gene, _IntMap))
		  bsRemove (Gene) ;
		if (bsAddKey (Gene, _IntMap, map) &&
		    bsAddData (Gene, _bsRight, _Int, &ia1))
		  bsAddData (Gene, _bsRight, _Int, &ia2) ;
		bsSave (Gene) ;
	      }
	    
	    cosmid = keySet (cosmids, 0) ;
	    if ((Seq = bsUpdate (cosmid)))
	      { 
		for (i = 0 ; i < keySetMax (ks) ; i++)
		  {
		    tg = keySet (ks, i) ;

		    if (
			(
			 (resizeAsGenefinder &&  bsFindKey (Seq, _Subsequence, tg)) ||
			 (! resizeAsGenefinder &&  bsFindKey (Seq, _Transcribed_gene, tg))
			) &&
			bsGetData (Seq, _bsRight, _Int, &ix1) &&
			bsGetData (Seq, _bsRight, _Int, &ix2))
		      {
			if (ix1 < ix2)
			  {
			    if (!i || ix1 < ia1) ia1 = ix1 ;
			    if (!i || ix2 > ia2) ia2 = ix2 ;
			  }
			else
			  {
			    if (!i || ix1 > ia1) ia1 = ix1 ;
			    if (!i || ix2 < ia2) ia2 = ix2 ;
			  }
		      }
		  }
		if (bsFindKey (Seq, _Genes, gene))
		  bsRemove (Seq) ;
		bsAddKey (Seq, _Genes, gene) ;
		bsAddData (Seq, _bsRight, _Int, &ia1) ;
		bsAddData (Seq, _bsRight, _Int, &ia2) ;
		bsSave (Seq) ;
	      }
	    look->pleaseRecompute = TRUE ;
	  }
	bsSave (Seq) ;
	keySetDestroy (ks) ;
	keySetDestroy (maps) ;
	keySetDestroy (cosmids) ;
      }
      break ;
    case 'k': /* kill */
      if (messQuery ("Do you really wish to kill this gene"))
	{
	  OBJ obj = 0 ;
	  if ((obj = bsUpdate (gene)))
	    bsKill (obj) ;	      

	  look->pleaseRecompute = TRUE ;
	}
      break ;
    case 'u': /* use AM */
	{
	  OBJ obj = 0 ;
	  if ((obj = bsUpdate (gene)))
	    {
	      bsAddTag (obj, _Use_AM) ;
	      bsSave (obj) ;	      
	    }

	  look->pleaseRecompute = TRUE ;
	}
      break ;
    case 'f': /* fuse */ 
      if (messPrompt ("Do you see RNA editing in this gene ?", "", "t"))
	{
	  KEY kk = 0 ;
	  char *cp ;
	  OBJ Gene = bsUpdate (gene) ;

	  if (Gene)
	    {
	      if ((cp = freeword ()) && *cp && lexaddkey (cp, &kk, _VText))
		bsAddKey (Gene, _RNA_editing, kk) ;
	      else
		bsAddTag (Gene, _RNA_editing) ;
	      bsSave (Gene) ;
	    }
	}
      break ;
     default:
      return ;
    } 
  graphActivate (old) ;
  if (look->pleaseRecompute)
    fMapDraw(look,0) ;
}

/***********************************************************************/

static void fMapMrnaAction (KEY key, int box)
{ 
  KEY gene = 0 ;
  SEG *seg ;
  Graph old = graphActive () ;
  FMAPLOOKGET("fMapMrnaAction") ;
  
  seg = BOXSEG(box) ;
  if (!seg)
    return ;
  gene = seg->parent ;
  
  switch (key)
    {
    case 'w':
      display (gene, 0, TREE) ;
      graphActivate (old) ;
      return ;
      break ;
    case 'D':
      displayPreserve () ;
      display (gene, 0, FMAP) ;  /* view this mrna as the main sequence object */
      graphActivate (old) ;
      return ;
      break ;
    case 't':
      look->translateGene = seg->parent ;
      
      if (seg->type & 0x1)
	{
	  mapColSetByName ("Down Gene Translation", FALSE) ;
	  mapColSetByName ("Up Gene Translation", TRUE) ;
	}
      else
	{
	  mapColSetByName ("Down Gene Translation", TRUE) ;
	  mapColSetByName ("Up Gene Translation", FALSE) ;
	}
      break ;
    case 'T': 
      mapColSetByName ("3 Frame Translation", 2) ;
      break ;
    default:
      return ;
    } 
  look->pleaseRecompute = TRUE ;
  graphActivate (old) ;
  fMapDraw(look,0) ;
}

/***********************************************************************/
#ifdef ACEMBLY
static KEY fMapTranscribedGeneAbsorbInGeneBoxTg = 0 ;
static void fMapTranscribedGeneAbsorbInGeneBoxPick (KEY gene)
{
  OBJ Tg = 0 ;
  FMAPLOOKGET("fMapTranscribedGeneAbsorbInGeneBoxPick") ;
  
  
  if (!(class(gene) == _VGene))
    {
      messout ("Please pick a genebox") ;
      return ;
    }
  if ((Tg = bsUpdate (fMapTranscribedGeneAbsorbInGeneBoxTg)))
    {
      bsAddKey (Tg, _Gene, gene) ;
      look->pleaseRecompute = TRUE ;
      bsSave (Tg) ;
    }
  fMapDraw (look, 0) ;
}
#endif
static void fMapTranscribedgeneAction (KEY key, int box)
{ 
  KEY gene = 0 ;
  SEG *seg ;
  Graph old = graphActive () ;
#ifdef ACEMBLY
  OBJ Gene = 0 ;
  int searchRepeats = 0 ;
#endif
  FMAPLOOKGET("fMapTranscribedgeneAction") ;
  
  seg = BOXSEG(box) ;
  if (!seg)
    return ;
  gene = seg->parent ;
#ifdef ACEMBLY
lao:
#endif
  switch (key)
    {
    case 'w':
      display (gene, 0, TREE) ;
      graphActivate (old) ;
      return ;
      break ;
    case 't':
      if (gene == myTranslatedGene)
	myTranslatedGene = 0 ;
      else
	myTranslatedGene = gene ;
      mapColSetByName ("3 Frame Translation", TRUE) ;
      break ;
    case 'r': /* rename */
      if (messPrompt("New name for this gene", name(gene),"w"))
	{
	  char *cp = strnew(freeword(), 0) ;
	  if (cp) 
	    { lexAlias (&gene, cp, TRUE, FALSE) ;
	    messfree (cp) ;
	    }
	}
      break ;
#ifdef ACEMBLY
    case 's': /* Recompute Splicing */
    case 'S': /* Limit this gene to the active zone */
    case 'T': /* Shed weakly connected variants */
    case 'U': /* split genebox by geneId */
     if (gene)
	{
	  KEY cosmid = keyGetKey (gene, _Genomic_sequence) ;
	  int z1= 0, z2 = 0, i = arrayMax (look->segs) ;
	  int myMin, myMax ;
	  SEG *seg1 = 0 ;
	  sessionAutoSave (60 * 30, 60 * 5) ;   /* save every 30 minutes, or after 5 minutes inactivity */
	  while (i--)
	    {
	      seg1 = arrp (look->segs, i, SEG) ;
	      if (seg1->key == cosmid && (seg1->type | 0x1) == SEQUENCE_UP)
		{ z1 = seg1->x1 ; z2 = seg1->x2 ; break ; }
	    }
	  if (z1 >= z2) break ;
	  if ( key == 'S')
	    { 
	      myMin = look->zoneMin ; myMax = look->zoneMax;  
	      if (seg1->type & 0x1) /* up cosmid */
		cDNARealignGene (gene, z2 - myMax, z2 - myMin, 1, FALSE, 1, 0, 2) ;
	      else
		cDNARealignGene (gene, myMin - z1, myMax - z1, 2, FALSE, 1, 0, 2) ;
	    }
	  else if ( key == 'U')
	    { 
	      myMin = look->zoneMin ; myMax = look->zoneMax;  
	      if (seg1->type & 0x1) /* up cosmid */
		cDNARealignGene (gene, z2 - myMax, z2 - myMin, 1, FALSE, 1, 0, 4) ;
	      else
		cDNARealignGene (gene, myMin - z1, myMax - z1, 2, FALSE, 1, 0, 4) ;
	    }
	  else 
	    { 
	      OBJ Cosmid = 0 ;
	      int direction = 2 ;

	      myMin = seg->x1 ; myMax = seg->x2 ; 
	      if ((Cosmid = bsCreate (cosmid)))
		{
		  if (bsFindKey (Cosmid, _Transcribed_gene, gene) &&
		      bsGetData (Cosmid, _bsRight, _Int, &myMin) &&
		      bsGetData (Cosmid, _bsRight, _Int, &myMax))
		    {
		      if (myMin > myMax)
			{ int dummy = myMin ; myMin = myMax ; myMax = dummy ; direction = 1 ; }
		    }
		  bsDestroy (Cosmid) ;
		}
	      if ( key == 's') /* just an indication */
		cDNARealignGene (gene, 0, 0, direction, TRUE, 0, 0, 2) ;
	      else if (key == 'T')
		{
		  /* doShed = 1 so we do not force calcul if there is a single mrna */
		  cDNARealignGene (gene, 0, 0, 0, FALSE, 2, 0, 1) ;
		}
	    }
	}
      sessionAutoSave (60 * 30, 60 * 5) ;   /* save every 30 minutes, or after 5 minutes inactivity */
      break ;
    case 'J':  /* fuse get rubber repeats*/
      searchRepeats = 2 ;
    case 'E':  /* fuse get repeats*/
      if (!searchRepeats)
	searchRepeats = 1 ;
      /* fall thru */
    case 'F':  /* fuse */
      if (!gene)
	break ;
      {
	KEYSET clones = 0, genes = 0, ks1 = 0, cosmids = 0, aa = 0; 
	KEY cosmid ;
	int z1= 0, z2 = 0, zz1 = 0, zz2 = 0, dz, 
	  u1, u2, du, bestDu, bestDz, 
	  a1, a2, i = arrayMax (look->segs), j, jcl = 0, jg = 0, jcs = 0 ;
	SEG *seg1 = 0 ;
	OBJ Gene = 0, Cosmid = 0 ;
	BOOL cosmidUp = FALSE ;
	BOOL doFuse = FALSE ; /* becuase the fusion is done inside this module */

	/* accumulate clones of other local genes */
	clones = keySetCreate () ; 
	genes = keySetCreate () ; 
	cosmids = keySetCreate () ; 
	while (i--)
	  {
	    seg1 = arrp (look->segs, i, SEG) ;
	    if (seg1->key == seg1->parent && seg1->type == TRANSCRIBEDGENE)
	      {
		z1 = seg1->x1 ; z2 = seg1->x2 ;  
		zz1 = z1 ; zz2 = z2 ;
	        cosmid = keyGetKey (seg1->key, _Genomic_sequence) ;
	        if (z2 > look->zoneMin && z1 < look->zoneMax)
		  {
		    keySet (genes, jg++) = seg1->key ;
		    ks1 = queryKey (seg1->key, "FOLLOW cDNA_clone") ;
		    for (j = 0 ; j < keySetMax(ks1) ; j++)
		      keySet(clones, jcl++) = keySet (ks1, j) ;
		    if (keySetMax(ks1))
		      keySet (cosmids, jcs++) = cosmid ; 
		    keySetDestroy (ks1) ;
		  }
	      }
	  }
	keySetSort (genes) ;
	keySetCompress (genes) ;
	keySetSort (clones) ;
	keySetCompress (clones) ;

	keySetSort (cosmids) ;
	keySetCompress (cosmids) ;
	cosmid = keySet(cosmids, 0) ; 
	keySetDestroy (cosmids) ;

	if (keySetMax(genes)) 
	  {
	    aa = keySetAlphaHeap (genes, keySetMax(genes)) ;
	    for (j = 0 ; j < keySetMax(aa) ; j++)
	      cdnaCleanUp (0, keySet(aa, j), 0) ;
	    gene = keySet(aa,0) ;
	    keySetKill (aa) ;
	    Gene = bsUpdate (gene) ;
	    if (Gene)
	      {
		KEYSET reads = query (clones, ">Read") ;
		for (j = 0 ; j < keySetMax(clones) ; j++)
		    bsAddKey (Gene, _cDNA_clone, keySet(clones, j)) ;
		for (j = 0 ; j < keySetMax(reads) ; j++)
		  bsAddKey (Gene, _Read, keySet(reads, j)) ;
		keySetDestroy (reads) ;
	      }
	    bsAddKey (Gene, _Genomic_sequence, cosmid);
	    if (bsFindTag (Gene, str2tag("To_be_fused_with")))
	      bsRemove (Gene) ; /* no more fusions than what i did here */
	    bsSave (Gene) ;
	  }
	else
	  goto fusionfailed ;

	/* get coord of cosmid */

	i = arrayMax (look->segs) ; 
	bestDz = bestDu = 0 ;
	while (i--)
	  {
	    seg1 = arrp (look->segs, i, SEG) ;
	    if ((seg1->type | 0x1) == SEQUENCE_UP &&
		keyFindTag (seg1->key, _Genomic))
	      { 
		/* find shortest genomic seg with longest intersect with zone */
		z1 = seg1->x1 ; z2 = seg1->x2 ;
		dz = z2 - z1 ;

		u1 = look->zoneMin > z1 ? look->zoneMin : z1 ;
		u2 = look->zoneMax < z2 ? look->zoneMax : z2 ;
		du = u2 - u1 ;

		if (du < bestDu) continue ;
		if (bestDz && dz > bestDz) continue ;
		
		bestDz = dz ; bestDu = du ;
		cosmid = seg1->key ; 
		zz1 = z1 ; zz2 = z2 ;
		
		if (seg1->type & 0x1) /* up cosmid */
		  cosmidUp = TRUE ;
		else
		  cosmidUp = FALSE ;
	      }
	  }
	if (cosmidUp) /* up cosmid */
	  { a1 = zz2 - look->zoneMin + 5 ; a2 = zz2 - look->zoneMax - 5 ; }
	else
	  { a1 = look->zoneMin - zz1 - 5 ; a2 = look->zoneMax - zz1 + 5 ; } 
	if ((Cosmid = bsUpdate (cosmid)))
	  {
	    bsAddKey (Cosmid, _Transcribed_gene, gene) ;
	    bsAddData (Cosmid, _bsRight, _Int, &a1) ;
	    bsAddData (Cosmid, _bsRight, _Int, &a2) ;
	    bsSave (Cosmid) ;
	  }
	cDNARealignGene (gene, a1, a2, cosmidUp ? 1 : 2, doFuse, 1, searchRepeats, 2) ;
	/* search and split_cloud */

      fusionfailed:
	keySetDestroy (clones) ;
	keySetDestroy (genes) ; 
	keySetDestroy (aa) ; 
      }
      break ;
    case 'k': /* kill */

      if (messQuery("Do you really want to eliminate this gene and all reads inside it ?"))
	{
	  KEYSET gks = queryKey (gene, "> cDNA_clone") ;
	  int i ;
	  
	  for (i = 0 ; i < keySetMax(gks) ; i++)
	    cDNARemoveCloneFromGene (keySet(gks,i),gene) ;

	  Gene = bsUpdate (gene) ; 
	  if (Gene) bsKill(Gene) ;
	}
      break ;
    case 'K': /* kill transcript */
       if (messQuery("Do you really want to eliminate this transcript ?") &&
	   checkWriteAccess())
	 {
	   Gene = bsUpdate (seg->parent) ; 
	   if (Gene) bsKill(Gene) ;
	 }
      break ;
     
    case 'R': /* revert to abi base call */
      if (messQuery 
   ("You will lose all the editions of this gene\nDo you wish to proceed ?"))
      {
	KEYSET ks = queryKey (gene, "{FOLLOW BaseCall} $| {FOLLOW SCF_Position}") ;
	int i = keySetMax(ks) ;
	while (i--)
	  arrayKill (keySet(ks,i)) ;
	key = 's' ;
	goto lao ;
      }
      else
	return ;
      break ;
    case 'P': /* change strand for whole gene */
      if (messQuery 
   ("Do you want to change the strand of all reads of this gene ?"))
	{ 
	  OBJ Est = 0 ;
	  KEYSET ks = queryKey (gene, ">Read") ;
	  int i = keySetMax(ks) ;
	  while (i--)
	    {
	      Est = bsUpdate (keySet(ks, i)) ;
	      if (Est)
		{
		  if (bsFindTag (Est, _Forward))
		    {
		      bsAddTag (Est, _Reverse) ;
		      bsAddTag (Est, _mReverse) ;
		      bsAddTag (Est, _Colour) ;
		      bsPushObj (Est) ;

		    }
		  else
		    {
		      bsAddTag (Est, _Forward) ; 
		      bsAddTag (Est, _mForward) ;
		      bsAddTag (Est, _Colour) ;
		      bsPushObj (Est) ;
		    }
		  bsSave (Est) ;
		}
	    }
	key = 's' ;
	 {
	    KEYSET ks = query (0, "Find mrna !from_gene && !from_prediction") ;
	    keySetKill (ks) ;
	    keySetDestroy (ks) ;
	 }
	goto lao ;
	}
      else
	return ;
      break ;
    case 'a': /* annotate */
      annotStart (gene) ;
      return ; /* no fmapredraw needed */
    case 'b': /* new annotator */
      geneAnnotDisplay (0, gene, 0) ;
      return ; /* no fmapredraw needed */
    case 'G': /* create gene box */
      {
	KEY gg = 0 ;

	if (lexword2key (name(gene), &gg, _VGene) &&
	    keyFindTag (gg, _SMAP))
	  {
	    messout ("Gene %s already exist, please chose another name", name(gg)) ;
	    return ;
	  }
	else
	  {
	    OBJ GG = 0, Tg = 0 ;
	    KEY map = 0, cosmid = 0 ;
	    char *cp ;
	    int x1, x2 ;
	    
	    lexaddkey (name(gene), &gg, _VGene) ;
	    
	    if ((GG = bsUpdate (gg)) && (Tg = bsCreate (gene)))
	      {
		if (keyFindTag (gene, _mRNA))
		  {
		    int im ;
		    KEYSET mrnas = queryKey (gene, ">mRNA") ;

		    for (im = 0 ; im < keySetMax (mrnas) ; im++)
		      bsAddKey (GG, _mRNA, keySet (mrnas, im)) ;
		  }		      
		bsAddKey (GG, _Transcribed_gene, gene) ;
		
		if (bsGetKey (Tg, _IntMap, &map) &&
		    bsGetData (Tg, _bsRight, _Int, &x1) &&
		    bsGetData (Tg, _bsRight, _Int, &x2))
		  {
		    bsAddKey (GG, _IntMap, map) ;
		    bsAddData (GG, _bsRight, _Int, &x1) ;
		    bsAddData (GG, _bsRight, _Int, &x2) ;
		  }
		bsSave (GG) ;
		if (bsGetKey (Tg, _Genomic_sequence, &cosmid) &&
		    bsGetData (Tg, _Covers, _Int, &x1) &&
		    bsGetData (Tg, _bsRight, _Text, &cp) &&
		    bsGetData (Tg, _bsRight, _Int, &x1) &&
		    bsGetData (Tg, _bsRight, _Int, &x2))
		  {
		    OBJ Cosmid = bsUpdate (cosmid) ;
		    
		    if (Cosmid)
		      {
			if (bsFindKey (Cosmid, _Transcribed_gene, gene))
			  bsAddKey (Cosmid, _Genes, gg) ;
			bsAddData (Cosmid, _bsRight, _Int, &x1) ;
			bsAddData (Cosmid, _bsRight, _Int, &x2) ;
			bsSave (Cosmid) ;
		      }
		  }
	      }
	    bsDestroy (Tg) ;
	    bsDestroy (GG) ;
	  }
      }
      break ; /* fmapredraw needed */
    case 'H': /* absorb in genebox */
      if (!checkWriteAccess ())
	return ; 
      fMapTranscribedGeneAbsorbInGeneBoxTg = gene ;
      displayBlock (fMapTranscribedGeneAbsorbInGeneBoxPick, 0) ;
      return ; /* no fmapredraw needed */
     case 'B': /* new annotator */
      if ((key = keyGetKey (gene, str2tag ("Derived_from_gene"))))
	geneAnnotDisplay (0, key, 0) ;
      return ; /* no fmapredraw needed */
    case 'D': /* show sequence in new window */
      {
	KEY dnaGene ;

	gene = seg->parent ;
	dnaGene = keyGetKey (gene, str2tag ("Spliced_sequence")) ;
	if (!dnaGene)
	  {
	    messout ("Sorry, dna not available ") ;
	    return ;
	  }

	displayPreserve () ;
	{
	  Graph oldg = graphActive() ;
	  void *look2 ;

	  display(dnaGene, 0, 0)  ;
	  mapColSetByName ("Gene Translation", TRUE) ;
	  mapColSetByName ("DNA", FALSE) ; 
	  mapColSetByName ("Gene Names", FALSE) ; 
	  if (fMapActive (0, 0, 0, &look2))
	    {
	      LOOK look3 = (LOOK) look2 ;
	      look3->pleaseRecompute = TRUE ;  
	      fMapDraw(look3,0) ;
	    }
	  graphActivate (oldg) ;
	}
	return ;
      }
      break ;
    case 'c': /* new annotator */
      callScript ("ace_sequin",name(gene)) ;
      return ; /* no fmapredraw needed */
#endif
    default:
      return ;
    } 
#ifdef ACEMBLY
    cDNAEliminateDeadMrnas () ;
#endif

  look->pleaseRecompute = TRUE ;
  graphActivate (old) ;
  fMapDraw(look,0) ;
}

/*******************************************************************/
/*******************************************************************/
