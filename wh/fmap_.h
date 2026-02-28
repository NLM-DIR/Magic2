/*  File: fmap_.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description:
 *      This is the internal header file for the fmap package, this file
 *      _must_only_ be included by fmap source files. The public header for
 *      users of the fmap package is fmap.h.
 *
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 17 13:46 1998 (fw)
 * * Sep  2 10:32 1998 (edgrif): Correct declaration for fMapFollowVirtualMultiTrace,
 *      should be void, not BOOL.
 * * Jul 28 16:26 1998 (edgrif): Add fMapInitialise() to function declarations.
 * * Jul 27 10:48 1998 (edgrif): Finish comments/tidy up.
 * * Jul 23 09:13 1998 (edgrif): Move fMapFindDNA into public header for use
 *      by dna code. (dnacpt.c and others).
 * Created: Thu Jul 16 09:29:49 1998 (edgrif)
 *-------------------------------------------------------------------
 */

#ifndef ACEDB_FMAP_P_H
#define ACEDB_FMAP_P_H
/*             NO CODE BEFORE THIS                          */

#include "ac.h"
#include "fmap.h"		/* Include our public header. */

#include "lex.h"		/* Other needed headers. */
#include "classes.h"
#include "tags.h"
#include "method.h" 
#include "bump.h"
#include "orepeats.h"

/************************************************************/

typedef struct
{ int   min, max ;
  float *cum, *revCum ;
} GFINFO ;

/************** FMAP completion of LOOK structure *******************/
struct LookStruct {
  magic_t *magic;			/* == &FMAPLOOK_MAGIC */
  AC_HANDLE handle, segsHandle ;
  KEY   seqKey, mapKey, dnaKey ;
  int   start, stop ;		/* limits of information in this look */
				/* pythagorean coords - oriented as dna */
  int   min, max ;		/* min and max DNA pos on screen */
  int   length, fullLength ;
  int   origin, zoneMin, zoneMax ;
  int   activeBox, centreBox, minLiveBox, summaryBox ;
  int   selectBox, zoneBox, originBox, segBox ;
  int   lastTrueSeg ;
  int   dnaWidth, dnaStart, dnaSkip ;		  /* set in makeSartSkip, used in virtualDna */
  float mrnaOffset ;
  char  originBuf[24], zoneBuf[24], centreBuf[24], oligoNameBuffer [16] ;
  char  segNameBuf[64], segTextBuf[512] ;
  Graph graph ;
  Array segs ;			/* array of SEG's from the obj */
  Array boxIndex ;		/* SEG if >minLiveBox, <maxLiveBox */
  Array dna, dnaR ;		/* dnaR in acembly */
  Array colors ;
  Array sites ;
  Stack siteNames ;
  Associator chosen, antiChosen ;
  Associator probeAss ;
  Array multiplets, contigCol ;
  Array homolInfo ;
  Array seqInfo ;
  Array minDna, maxDna ;
  Array decalCoord ;
  Array tframe ;
  int   flag ;
  MAP   map ;
  GFINFO gf ;
  KEY translateGene ;
  Graph selectGraph ;		/* graph for showing selections */
  int pleaseRecompute ;		/* force recompute (invoked from dnacpt etc) */
  Associator virtualErrors, taggedBases ;
  int DNAcoord ;		/* valid when BOXSEG(look->activeBox)->type == DNA_SEQ */
  Array oligopairs ;
  BOOL isGeneFinder ;		/* triggers the genefinder menu in acembly */
  BitSet homolBury ;		/* same index as homolInfo */
  BitSet homolFromTable ;	/* MatchTable or Tree? (index as homolInfo) */
  DICT *featDict ;
  DICT *bcDict ;
  Array bcArray ;		/* for local column info */
  Array visibleIndices ;
  BOOL isRetro ; /* needed to correctly paint the triplet corresponding to an amino acid */
  KEY  seqOrig ;		/* the orginal Sequence key that fmapCreateLookfm() was called on */
  BOOL isOspSelecting ; /* set/unset in fmaposp, to prevent deregistering PICK */
  BOOL decorateSplicedcDNA ; /* set to true by  fMapcDNAShowSplicedcDNA */
  Stack geneNamesStack ; /* push genenames there */
  OREPEAT *oligoRepeats ;   /* to show profiles of all short repeats in genome */
  KEY view ;  /* to control the display */
} ;

/************************************************************/

#define FLAG_RNA			0x0001
#define FLAG_COMPLEMENT 		0x0002
#define FLAG_REVERSE			0x0004
#define FLAG_COLOR_HID			0x0008
#define FLAG_HAVE_GF			0x0010
#define FLAG_COMPLEMENT_SUR_PLACE 	0x0020
#define FLAG_HIDE_HEADER 		0x0040
#define FLAG_AUTO_WIDTH 		0x0080
#define FLAG_VIRTUAL_ERRORS 		0x0100
#define FLAG_HIDE_SMALL_CONTIGS 	0x0200
#define FLAG_COLOR_CONTIGS		0x0400
#define FLAG_COLUMN_BUTTONS		0x0800

/************************************************************/

#define COORD(look, x) \
   ((look->flag & FLAG_REVERSE) ? \
   look->length - (x) - look->origin \
   : (x) - look->origin + 1)

#define BOXSEG(i) \
  ((i)<arrayMax(look->boxIndex) && arr(look->boxIndex, (i), int) < arrayMax(look->segs) ? \
  arrp(look->segs, arr(look->boxIndex, (i), int), SEG) : 0)

/************************************************************/

		/* addresses of the following are unique ids */
extern magic_t GRAPH2FMAPLOOK_ASSOC ; /* for graph -> FMAPLOOK  */
extern magic_t FMAPLOOK_MAGIC ;	  /* for extra verification on FMAPLOOK */

#define FMAPLOOKGET(name) \
  LOOK look ; \
  if (!(look = fMapGifLook) && \
      !graphAssFind (&GRAPH2FMAPLOOK_ASSOC, &look)) \
    messcrash ("graph not found in %s", name) ; \
  if (look->magic != &FMAPLOOK_MAGIC) \
    messcrash ("%s received a wrong FMAPLOOK pointer", name)

/***************** segs ***********************************/

/************************
 **** VERY IMPORTANT ****
 ************************
 If you change this enum then you must correspondingly change
 the list of names in static char* fMapSegTypeName[] in fmapcontrol.c
 *************************/

typedef enum { 
  MASTER=0, ASSEMBLY_TAG,

  SEQUENCE, SEQUENCE_UP,	/* UP must be odd */
  CDS, CDS_UP,
  INTRON, INTRON_UP,
  EXON, EXON_UP,		/* draw exons after introns */
  EXON_GAP, EXON_GAP_UP,		/* draw exons after introns */
  EMBL_FEATURE, EMBL_FEATURE_UP, 
  FEATURE, FEATURE_UP,
  ATG, ATG_UP,
  SPLICE3, SPLICE3_UP,
  SPLICE5, SPLICE5_UP,
  CODING, CODING_UP,
  TRANSCRIBEDGENE, TRANSCRIBEDGENE_UP,
  PMRNA, PMRNA_UP,
  MRNA, MRNA_UP,
  MPRODUCT, MPRODUCT_UP,
  TRANSCRIPT, TRANSCRIPT_UP,
  EXON_SUPPORT,  EXON_SUPPORT_UP,
  SPLICED_cDNA, SPLICED_cDNA_UP,
  SPLICED_cDNA_DECORATE, SPLICED_cDNA_DECORATE_UP,
  MGENES, MGENES_UP,
  MSOLEXA, MSOLEXA_UP,
  MRNAI, MRNAI_UP,
  MOST, MOST_UP,
  CDNA_GENE_NAME, CDNA_GENE_NAME_UP,
  VIRTUAL_SUB_SEQUENCE_TAG,	/* a subseq without recursion */ 
  VIRTUAL_SUB_SEQUENCE_TAG_UP,
  VIRTUAL_MULTIPLET_TAG,
  VIRTUAL_MULTIPLET_TAG_UP,
  VIRTUAL_ALIGNED_TAG,
  VIRTUAL_ALIGNED_TAG_UP,
  VIRTUAL_PREVIOUS_CONTIG_TAG,
  VIRTUAL_PREVIOUS_CONTIG_TAG_UP,
  HOMOL, HOMOL_UP,
  PRIMER, PRIMER_UP,  /* useful for Directed sequencing */
  PROBE, PROBE_UP,
  OLIGO, OLIGO_UP,
  OLIGO_PAIR, OLIGO_PAIR_UP,
  TRANS_SEQ, TRANS_SEQ_UP,
  ALLELE, ALLELE_UP,    /* !!!!!!!!!!!!! allele is last upable tag !!!!!!!!!!!!! */

  VISIBLE,
  DNA_SEQ, PEP_SEQ, ORF, 
  VIRTUAL_PARENT_SEQUENCE_TAG, /* mieg, a subseq, but without recursion */ 
  VIRTUAL_CONTIG_TAG,
  CLONE_END
} SegType ;

typedef struct {
  KEY	  key, parent ;
  KEY     source ;		/* object the field comes from */
  SegType type ;
  int	  x1,x2 ;
  int sourceDx ;   /* length of source, needed in fMapRC :reverseComplement */
  union {
    float f ;
    int i ;
    KEY k ; 
    char *s ;
  } data ;		/* utility slot */
} SEG ;

#define segFormat "kkiiif"

#define HASH(x,y)		((char*)((char *)0)  + ((x) << 24) + (y))
#define SEG_HASH(seg)		HASH((seg)->type,(seg)->x2)
#define UNHASH_TYPE(x)		((int)((x) - ((char *)0)) >> 24)
#define UNHASH_X2(x)		((int)((x) - ((char *)0)) & 0xffffff)

/*************** extra sequence info ********************/

#define SEQ_CANONICAL		0x00000001
#define SEQ_VISIBLE		0x00000002
#define SEQ_EXONS		0x00000004
#define SEQ_VIRTUAL_ERRORS	0x00000008
#define SEQ_SCORE		0x00000010

#define SEQ_CONFIRM_UNKNOWN	0x00010000
#define SEQ_CONFIRM_EST		0x00020000
#define SEQ_CONFIRM_HOMOL	0x00040000
#define SEQ_CONFIRM_CDNA	0x00080000
#define SEQ_CONFIRM_UTR		0x00100000
#define SEQ_CONFIRMED		0x001f0000

typedef struct {
  KEY method ;
  float score ;
  KEY flags ;
} SEQINFO ;

typedef struct { KEY key ; float y1, y2 ; int iseg ; } KEYZONE ;

/**************** tags and classes ***********************/

extern KEY _Arg1_URL_prefix ;
extern KEY _Arg1_URL_suffix ;
extern KEY _Arg2_URL_prefix ;
extern KEY _Arg2_URL_suffix ;

/**********************************************************/

/* Private, non-static routines.                                             */
/*                                                                           */
				/* fmapcontrol.c */
void fMapInitialise (void) ;
BOOL fMapFindSegBounds (LOOK look, SegType type, int *min, int *max) ;
void fMapDumpSegs (LOOK look, int version, 
		   KEY *refSeq, int *offPtr, 
	/* refseq is the coordinate basis (take minimal spanning seq if 0) */
		   KEYSET methset) ;
	/* if non-zero, take methset as a set of methods to dump */
void fMapSetZone (void) ;
BOOL fMapFindZoneFather (void *vv, int min, int max, KEY *fp, int *originp) ;
float fMapSquash (float value, float midpoint) ;
BOOL fMapConvert (LOOK look, BOOL force) ;
int fMapOrder (const void *a, const void *b) ;
void fMapRC (LOOK look) ;
extern char* fMapSegTypeName[] ;
void fMapProcessMethods (LOOK look) ;
BOOL fMapActivateGraph (void) ;
BOOL fMapFindSpan (LOOK look, KEY *key, int *x, int *y) ;
void fMapReportLine (LOOK look, SEG *seg, BOOL isGIF, float x) ;

				/* fmapsequence.c */
int sequenceLength (KEY seq) ;
void fMapShowDNA (LOOK look, float *offset) ;
void fMapSetDnaWidth (void) ;
void fMapToggleAutoWidth (void) ;
void fMapShowCoords (LOOK look, float *offset) ;
void fMapShowSummary (LOOK look, float *offset) ;
void fMapClear (void) ;
void fMapColorIntronsExons (int box) ;
void fMapExportTranslation (int box) ;
void fMapExportTranslations (void) ;
void fMapExportcDNA (int box) ;
void fMapFindCoding (LOOK look) ;
void fMapShowATG (LOOK look, float *offset) ;
void fMapShowCoding (LOOK look, float *offset) ;
void fMapShowOrigin (LOOK look, float *offset) ;
void fMapShowStatus (LOOK look, float *offset) ;
BOOL fMapGetCDS (LOOK look, KEY parent, Array *cds, Array *index);

				/* fmapfeatures.c */
void fMapShowMiniSeq (LOOK look, float *offset) ;
void fMapShowScale (LOOK look, float *offset) ;
void fMapShowSequence (LOOK look, float *offset) ;
void fMapShowSoloConfirmed (LOOK look, float *offset) ;
void fMapShowEmblFeatures (LOOK look, float *offset) ;
void fMapShowText (LOOK look, float *offset) ;
void fMapShowFeature (LOOK look, float *offset) ;
void fMapShowHomol (LOOK look, float *offset) ;
void fMapShowBriefID (LOOK look, float *offset) ;
void fMapShowTitles (LOOK look, float *offset) ;
void fMapAddSegments (void) ;
void fMapClearSegments (void) ;
				/* mieg, shift right by score */
void fMapShowGF_seg (LOOK look, float *offset) ;
void fMapShowSplices (LOOK look, float *offset) ;
void fMapShowCptSites (LOOK look, float *offset) ;
void fMapShowAssemblyTags (LOOK look, float *offset) ;
void fMapShowCDSBoxes (LOOK look, float *offset) ;
void fMapShowCanonical (LOOK look, float *offset) ;
void fMapShowGeneNames (LOOK look, float *offset) ;
void fMapShowCDSLines (LOOK look, float *offset) ;
void fMapShowCDNAs (LOOK look, float *offset) ;
void fMapShowText (LOOK look, float *offset) ;
void fMapShowUserSegments (LOOK look, float *offset) ;
void fMapRemoveSelfHomol (LOOK look) ;
void fMapShowAlleles (LOOK look, float *offset) ;
void keyZoneAdd (Array a, KEY key, float y1, float y2, int iseg);
int  keyZoneOrder (const void *va, const void *vb) ;
 
				/* fmapgene.c */
extern FREEOPT fMapChooseMenu[] ;
void fMapChooseMenuFunc (KEY key, int box) ;
extern MENUOPT fMapGeneOpts[] ;
void fMapAddGfSegs (void) ;

                 /*fmaposp.c*/
void fMapOspInit (void) ;
void fMapOspDestroy (LOOK look) ;
void fMapOspShowOligo (LOOK look, float *offset) ;               
BOOL fMapOspPositionOligo (LOOK look, SEG *seg, KEY oligo, int *a1p, int *a2p) ;
void fMapOspFindOligoPairs (LOOK look) ;    /* make virtual multiplets */
void fMapOspShowOligoPairs (LOOK look, float *offset) ;               
void fMapOspSelectOligoPair (LOOK look, SEG *seg) ;
				/* fmapassembly.c */
void fMapShowAlignments (LOOK look, float *offset) ;
void fMapShowPreviousContigs (LOOK look, float *offset) ;
void fMapShowContig (LOOK look, float *offset) ;
void fMapShowTraceAssembly (LOOK look, float *offset) ;
void fMapShowVirtualDna (LOOK look, float *offset) ;
void fMapShowVirtualMultiplets (LOOK look, float *offset) ;
BOOL fMapFollowVirtualMultiTrace (int box) ;
void fMapTraceDestroy (LOOK look) ;
void fMapTraceFindMultiplets (LOOK look) ;  /* make virtual multiplets */
void fMapSelectVirtualMultiplet (LOOK look, SEG *seg) ;
void fMapShowCloneEnds (LOOK look, float *offset) ;
void fMapTraceForget (KEY key) ;
void fMapTraceForgetContig (KEY contig) ;

void fMapcDNAShowSolexa (LOOK look, float *offset) ;               
void fMapcDNAShowProbe (LOOK look, float *offset) ;               
BOOL fMapcDNAProbePosition (LOOK look, OBJ Mrna, SEG *seg, KEY probe, int *p1, int *p2) ;
void fMapcDNADecorateSplicedcDNA (LOOK look, float *offset) ;
void fMapcDNAGeneName (LOOK look, float *offset) ;
void fMapcDNAShowExonSupport (LOOK look, float *offset) ;
void fMapcDNAShowSplicedcDNA (LOOK look, float *offset) ;
void fMapcDNAShowTranscribedgene (LOOK look, float *offset) ;
void fMapcDNAShowTranscript (LOOK look, float *offset) ;
void fMapcDNAShowMrna (LOOK look, float *offset) ;
void fMapcDNAShowPMrna (LOOK look, float *offset) ;
void fMapcDNAShowGenes (LOOK look, float *offset) ;
void fMapcDNAShowRNAi (LOOK look, float *offset) ;
void fMapcDNAShowOST (LOOK look, float *offset) ;
void fMapcDNAShowMProduct (LOOK look, float *offset) ;
void fMapShowcDNAs (LOOK look, float *offset) ;
void fMapcDNAShowOligoRepeats (LOOK look, float *offset) ;
void fMapShowGeneWalls (LOOK look, float *offset) ;
BOOL fMapcDNAFollow (int box) ;
void fMapcDNAInit (void) ;
BOOL fMapcDNADoFillData(SEG *seg) ;
void fMapcDNAReportLine (char *buffer, SEG *seg1, int maxBuf) ;
void fMapcDNAReportTranscribedgene (char *buffer, SEG *seg1, int maxBuf) ;
void fMapcDNAReportTranscript (char *buffer, SEG *seg1, int maxBuf) ;
void fMapcDNAReportMrna (char *buffer, SEG *seg1, int maxBuf) ;
void fMapcDNAReportGenes (char *buffer, SEG *seg1, int maxBuf) ;
void fMapcDNAReportRNAi (char *buffer, SEG *seg1, int maxBuf) ;
void fMapcDNAReportOST (char *buffer, SEG *seg1, int maxBuf) ;
void fMapcDNAReportProbe (char *buffer, SEG *seg1, int maxBuf) ;
void fMapcDNAReportMProduct (char *buffer, SEG *seg1, int maxBuf) ;
void fMapColorSeg (SEG *seg) ;
BOOL dnaParse (int level, KEY key) ;
BOOL hexAddSegs (char *name, int type, KEY key, 
		 int step, float thresh, BOOL isRC,
		 char *dna, int len, float* partial) ;
BOOL geneFinderAce (char *seq, void *gf) ;
void cDNAAlignInit (void) ;

/*                  NO CODE AFTER THIS                                       */
#endif /* ACEDB_FMAP_P_H */

