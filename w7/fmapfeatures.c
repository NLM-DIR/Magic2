/*  file: fmapfeatures.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: feature display for fmap package
 * Exported functions:
 * HISTORY:
 * Last edited: Apr 20 05:02 1999 (rd)
 * * Sep 16 10:47 1998 (edgrif): Insert macro to define bitField for
 *              bitset.h ops.
 * * Jul 16 10:06 1998 (edgrif): Introduce private header fmap_.h
 * Created: Sun Jul 26 13:02:46 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: fmapfeatures.c,v 1.36 2019/03/04 22:44:17 mieg Exp $ */

#include "fmap_.h"
#include "call.h"
#include "bitset.h"

#include "display.h"
#include "forest.h"
#include "query.h"
#include "sysclass.h"
#include "systags.h"
#include "bump.h"
#include "tags.h"
#include "a.h"
#include "peptide.h"
#include "freeout.h"
#include "dna.h"
#include "parse.h"
#include "cdna.h"
#include "main.h"

typedef struct {
  LOOK look ;
  KEY method ;
				/* following under user control */
  BOOL byWidth, byOffset, byHist ;
  BOOL frameSensitive, strandSensitive ;
  BOOL blixemX, blixemN ;
  BOOL isBump, isCluster ;
  BOOL isShowText ;
  float minScore, maxScore ;
  int colour ;
  float histBase ;
  float width ;
  float minMag, maxMag ;
				/* following not under user control */
  float offset ;
  BOOL isDown ;
  BOOL isFrame ;
  BOOL autoWidth ;		/* set FALSE if method or user sets width */
  enum { DEFAULT=0, WIDTH, OFFSET, HIST } mode ;
  float fmax ;
  BUMP bump ;
  Associator cluster ;		/* only non-zero if isCluster */
  int clusterCount ;
} BoxCol ;

static void   bcCheck (BoxCol *bc);   /* checks before using a column */
static BoxCol *bcFromName (LOOK look, char *name) ;
static BOOL   bcTestMag (BoxCol *bc, float mag) ;
static int    bcDrawBox (BoxCol *bc, SEG *seg, float score) ;
extern void mrnaTransferPg2PredictedMrna (KEYSET ks) ;
static void addVisibleInfo (LOOK look, int i) ;

/*************************************/
/* mhmp 01.10.99
static void fMapBoxInfo (LOOK look, int box, SEG *seg)*/
void fMapBoxInfo (LOOK look, int box, void *aseg)
{ 
  extern void (*gifEntry) (KEYSET) ;	/* entry point in command.c */
  SEG *seg = (SEG *) aseg ;

  if (!gifEntry)		/* don't waste time if not in giface */
    return ;

  fMapReportLine (look, seg, TRUE, 0) ;

  graphBoxInfo (box, class (seg->key) ? seg->key : seg->parent, 
		strnew (look->segTextBuf, look->segsHandle)) ;
}

/*************************************/

void fMapShowMiniSeq (LOOK look, float *offset)
{
  int i ;
  float y, y0 = -10, newoff, x0 ;
  SEG *seg ;

  if (!fmapView (look->view, str2tag ("Fmap_Locator")))
    return ;
  newoff = *offset + 5 ;
  mapShowLocator (look, &newoff) ;
  x0 = (look->flag & FLAG_REVERSE) ? look->map->max : 0 ;

     /* Position main Markers = loci just now */
  graphTextHeight (0.75) ;
  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs,i,SEG) ;
      if (class (seg->key) == _VLocus)
	{ y = MAP2WHOLE (look->map, (seg->x1 + seg->x2 - 2*x0)/2.0) ;
	  if (y > y0 + 2)
	    graphText (name (seg->key), *offset+0.5, y-0.4) ;
	  y0 = y ;
	}
    }
  graphTextHeight (0) ;
  *offset = newoff ;
}

/************************************/

void fMapShowScale (LOOK look, float *offset)
{
  float unit, subunit ;
  int iUnit, iSubunit, type, unitType ;
  int x, xf, xr, width = 0 ;
  float y ;
  char *cp, unitName[] = { 0, 'k', 'M', 'G', 'T', 'P' }, buf[2] ;
  BOOL isReverse ;

  xf = look->origin - 1 ;
  xr = look->length - look->origin ;
  isReverse = ((look->flag & FLAG_REVERSE) != 0) ;

  unit = subunit = 1.0 ;
  mapFindScaleUnit (look->map, &unit, &subunit) ;
  iUnit = unit + 0.5 ;
  iSubunit = subunit + 0.5 ;

  if (isReverse)
    x = xr - iUnit * ((xr - look->min)/iUnit) ;
  else
    x = xf + iUnit * ((look->min - xf)/iUnit) ;
  if (x < look->min - subunit) x += unit ;
  
  for (type = 1, unitType = 0 ; unit > 0 && 1000 * type < unit && unitType < 5; 
       unitType++, type *= 1000) ;
  
  for ( ; x < look->max +  subunit ; x += unit)
    { y = MAP2GRAPH (look->map, x) ;
      graphLine (*offset, y, *offset+1, y) ;
      buf[0] = unitName[unitType] ; buf[1] = 0 ;
      if (fmapView (look->view, str2tag ("Fmap_scale_coords")))
	{
	  cp = messprintf ("%d%s", COORD (look,x)/type, buf) ;
	  if (width < strlen (cp) + 2)
	    width = strlen (cp) ;
	  graphText (cp, *offset+1.5, y-0.5) ;
	}
    }

  if (isReverse)
    x = xr - iSubunit * ((xr - look->min)/iSubunit) ;
  else
    x = xf + iSubunit * ((look->min - xf)/iSubunit) ;
  if (x < look->min - subunit) x += subunit ;
  for ( ; x < look->max + subunit; x += subunit)
    { y = MAP2GRAPH (look->map, x) ;
      graphLine (*offset+0.5, y, *offset+1, y) ;
    }

  graphLine (*offset+1, topMargin + .5, *offset+1, mapGraphHeight-0.5) ;
  *offset += width + 2 ;
}

/**********************************************************/
/**********************************************************/

/*********************** geneMenu ***********************/

static void translateGene (int box)
{
  BOOL isNewGene = FALSE;
  SEG *seg ;
  FMAPLOOKGET ("translateGene") ;

  if (box >= arrayMax (look->boxIndex) || ! (seg = BOXSEG (box)))
    { messerror ("problem picking boxes for translation");
      return;
    }

  if (look->translateGene)	/* there was a previous selection */
    if (look->translateGene != seg->parent)
      {
	/* the user has selected a different gene for translation */
	isNewGene = TRUE ;
	look->pleaseRecompute = TRUE ;
      }

  /* select the gene for translation, from the selected SEG */
  look->translateGene = seg->parent ;

  /* mhmp 09.09.98 */
  if (seg->type == EXON_UP)
    {
      mapColSetByName ("Down Gene Translation", FALSE) ;
      mapColSetByName ("Up Gene Translation", isNewGene ? TRUE : 2);
    }
  else
    {
      mapColSetByName ("Down Gene Translation", isNewGene ? TRUE : 2);
      mapColSetByName ("Up Gene Translation", FALSE) ;
    }

  fMapDraw (look, 0) ;

  return;
} /* translateGene */

static void showInPepDisp (int box)
{ 
  SEG *seg ;
  Array pep ;
  FMAPLOOKGET ("translateGene") ;

  if (box < arrayMax (look->boxIndex) && (seg = BOXSEG (box)) &&
      (pep = peptideGet (seg->parent)))
    { arrayDestroy (pep) ;
      display (seg->parent, look->seqKey, PEPMAP) ;
    }
  else
    messerror ("Can't generate protein sequence") ;
}

void fMapShowcDNA (int box)
{
  KEY key ;
  Array cds = 0 ;
  SEG *seg ; 
  Stack s ;
  FMAPLOOKGET ("fMapShowcDNA") ;
    
  if (box >= arrayMax (look->boxIndex) ||
      ! (seg = BOXSEG (box)) ||
      !fMapGetCDS (look, seg->parent, &cds, 0))
    { messerror ("Show cDNA only works on coding genes") ;
      return ;
    }

  s = stackCreate (500) ;
  freeOutSetStack (s) ;
  lexaddkey (messprintf ("%s.cDNA",name (seg->parent)),  &key, _VSequence) ;
  freeOutf ("Sequence %s\ncDNA\nCDS\n\nDNA %s\n",
			 name (key), name (key)) ;
  dnaDecodeArray (cds) ;
  freeOut (arrp (cds,0,char)) ;
  freeOut ("\n\n") ;
  arrayDestroy (cds) ;

  parseBuffer (stackText (s,0), 0) ;
  stackDestroy (s) ;
  displayPreserve () ;
  display (key, 0, 0)  ;
}

#ifdef ACEMBLY
static void fMapCreateGeneBox (int box)
{
  KEY seq = 0, gg = 0 ;
  SEG *seg ; 
  FMAPLOOKGET ("fMapCreategeneBox") ;
    
  if (! _VGene ||
      box >= arrayMax (look->boxIndex) ||
      ! (seg = BOXSEG (box)) ||
      ! (seq = seg->parent))
    { messerror ("Confusion in  fMapCreategeneBox, sorry") ;
      return ;
    }

  if (!checkWriteAccess ())
    return ;
  if (lexword2key (name (seq), &gg, _VGene) &&
      keyFindTag (seq, str2tag ("SMAP")))
    {
      messout ("Gene %s already exist, please chose another name", name (gg)) ;
      return ;
    }
  else
    {
      OBJ GG = 0, Seq = 0 ;
      KEY map = 0, cosmid = 0 ;
      int x1, x2 ;
      
      lexaddkey (name (seq), &gg, _VGene) ;
      
      if ((GG = bsUpdate (gg)) && (Seq = bsCreate (seq)))
	{
	  bsAddKey (GG, str2tag ("Genefinder"), seq) ;
	  if (bsGetKey (Seq, _IntMap, &map) &&
	      bsGetData (Seq, _bsRight, _Int, &x1) &&
	      bsGetData (Seq, _bsRight, _Int, &x2))
	    {
	      bsAddKey (GG, _IntMap, map) ;
	      bsAddData (GG, _bsRight, _Int, &x1) ;
	      bsAddData (GG, _bsRight, _Int, &x2) ;
	    }
	  bsSave (GG) ;
	  bsGetKey (Seq, _Source, &cosmid) ;
	  bsDestroy (Seq) ;
	    
	  if (cosmid)
	    {
	      OBJ Cosmid = bsUpdate (cosmid) ;
	      
	      if (Cosmid)
		{
		  if (bsFindKey (Cosmid, _Subsequence, seq) &&
		      bsGetData (Cosmid, _bsRight, _Int, &x1) &&
		      bsGetData (Cosmid, _bsRight, _Int, &x2))
		    {
		      bsAddKey (Cosmid, str2tag ("Genes"), gg) ;
		      bsAddData (Cosmid, _bsRight, _Int, &x1) ;
		      bsAddData (Cosmid, _bsRight, _Int, &x2) ;
		    }
		  bsSave (Cosmid) ;
		}
	    }
	}
      bsDestroy (Seq) ;
      bsDestroy (GG) ;
    }
  look->pleaseRecompute = TRUE ;
  fMapDraw (look, 0) ;
}

static KEY  fMapAbsorbInGeneBoxSeq = 0 ;
static void fMapAbsorbInGeneBoxPick (KEY gene)
{
  OBJ Pg = 0, Gene = 0 ;
  int i ;
  KEY meth = 0 ;
  KEYSET ks = 0 ;
  FMAPLOOKGET ("fMapCreategeneBoxPick") ;
    

  if (! (class (gene) == _VGene))
    {
      messout ("Please pick a genebox") ;
      return ;
    }
  if ((Pg = bsUpdate (fMapAbsorbInGeneBoxSeq)))
    {
      if (bsFindTag (Pg, _Method))
	{
	  if (strncmp (name (fMapAbsorbInGeneBoxSeq), "hmm", 3))
	    {
	      bsRemove (Pg) ;
	      lexaddkey ("Absorbed", &meth, _VMethod) ;
	      bsAddKey (Pg, _Method, meth) ;
	    }
	}
      bsAddKey (Pg, str2tag ("Model_of_gene"), gene) ;
      look->pleaseRecompute = TRUE ;
      bsSave (Pg) ;
    }
  if ((Gene = bsUpdate (gene)))
    {
      ks = queryKey (fMapAbsorbInGeneBoxSeq, ">NM_id") ;
      for (i = 0 ; i < keySetMax (ks) ; i++)
	bsAddKey (Gene, str2tag ("NM_Id"), keySet (ks,i)) ;
      keySetDestroy (ks) ;
      ks = queryKey (fMapAbsorbInGeneBoxSeq, ">Locusid") ;
      for (i = 0 ; i < keySetMax (ks) ; i++)
	bsAddKey (Gene, str2tag ("LocusId"), keySet (ks,i)) ;
      keySetDestroy (ks) ;

      bsSave (Gene) ;
    }
  fMapDraw (look, 0) ;
}

static void fMapAbsorbInGeneBox (int box)
{
  KEY seq = 0 ;
  SEG *seg ; 
  FMAPLOOKGET ("fMapCreategeneBox") ;
    
  if (! _VGene ||
      box >= arrayMax (look->boxIndex) ||
      ! (seg = BOXSEG (box)) ||
      ! (seq = seg->parent))
    { messerror ("Confusion in  fMapAbsorbInGeneBox, sorry") ;
      return ;
    }

  if (!checkWriteAccess ())
    return ; 
  fMapAbsorbInGeneBoxSeq = seq ;
  displayBlock (fMapAbsorbInGeneBoxPick, 0) ;
}

static void fMapCreatePgMrna (int box)
{
  KEY seq = 0 ;
  SEG *seg ; 
  FMAPLOOKGET ("fMapCreategeneBox") ;
    
  if (! _VGene ||
      box >= arrayMax (look->boxIndex) ||
      ! (seg = BOXSEG (box)) ||
      ! (seq = seg->parent))
    { messerror ("Confusion in  fMapAbsorbInGeneBox, sorry") ;
      return ;
    }

  if (!checkWriteAccess ())
    return ; 

  if (0 && keyGetKey (seq, str2tag ("Predicted_mRNA")))
    messout ("This gene already has a corresponding predicted_mRNA") ;
  else
    {
      KEYSET ks = keySetCreate () ;
      keySet (ks, 0) = seq ;
      mrnaTransferPg2PredictedMrna (ks) ;
      keySetDestroy (ks) ;
      look->pleaseRecompute = TRUE ;
      fMapDraw (look, 0) ;
    }
}
#endif

#if !defined (MACINTOSH)
static void fMaptRNAfold (int box) ;
#endif

static MENUOPT geneMenu[] = {
  { (GraphFunc)translateGene, "Show translation" },
  { (GraphFunc)fMapShowcDNA, "Display cDNA in new window" },
  { (GraphFunc)showInPepDisp, "Display translation in protein window" },
  { (GraphFunc)fMapColorIntronsExons, "Color Exons" },
#ifdef ACEMBLY
  { (GraphFunc)fMapCreateGeneBox, "Create genebox" },
  { (GraphFunc)fMapAbsorbInGeneBox, "Absorb in genebox" },
  { (GraphFunc)fMapCreatePgMrna, "Create predicted mRNA" },
#endif
#ifndef JUNK     /* really this is junk */
  { (GraphFunc)fMapExportTranslation, "Export translation" },
  { (GraphFunc)fMapExportcDNA, "Export cDNA" },
#if !defined (MACINTOSH) 
  { (GraphFunc)fMaptRNAfold, "Fold tRNA" },
#endif
#endif
  { 0, 0 }
} ;

/*********************** intronMenu ***********************/

#ifndef ACEMBLY

static void intronConfirm (int box, int flag)
{
  SEG *seg ;
  OBJ obj ;
  KEY key ;
  int x, y ;
  SEQINFO *oldSinf, *newSinf ;
  FMAPLOOKGET ("intronConfirm") ;
  
  if (box >= arrayMax (look->boxIndex) || ! (seg = BOXSEG (box)))
    { messerror ("Problem picking box in intronConfirm") ; 
      return ; 
    }

  key = 0 ; x = seg->x1+1 ; y = seg->x2+1 ;
  if (!fMapFindSpan (look, &key, &x, &y) || ! (obj = bsUpdate (key)))
    { messerror ("Can't find object to store confirmation info") ;
      return ;
    }

  bsAddData (obj, str2tag ("Confirmed_intron"), _Int, &x) ;
  bsAddData (obj, _bsRight, _Int, &y) ;
  bsPushObj (obj) ;
  switch (flag)
    { 
    case SEQ_CONFIRM_EST: bsAddTag (obj, str2tag ("EST")) ; break ;
    case SEQ_CONFIRM_HOMOL: bsAddTag (obj, str2tag ("Homology")) ; break ;
    case SEQ_CONFIRM_CDNA: bsAddTag (obj, str2tag ("cDNA")) ; break ;
    case SEQ_CONFIRM_UTR: bsAddTag (obj, str2tag ("UTR")) ; break ;
    case 0: bsPrune (obj) ; break ;
    default: messerror ("Unrecognised flag in intronConfirm") ; break ;
    }
  bsSave (obj) ;

  newSinf = arrayp (look->seqInfo, arrayMax (look->seqInfo), SEQINFO) ;
  oldSinf = arrp (look->seqInfo, seg->data.i, SEQINFO) ;	/* must follow */
  seg->data.i = arrayMax (look->seqInfo)-1 ;
  *newSinf = *oldSinf ;
  if (flag)
    { newSinf->flags |= flag ;
      graphBoxDraw (box, -1, newSinf->flags & SEQ_CONFIRM_UTR ? YELLOW : GREEN) ;
    }
  else
    { newSinf->flags &= ~SEQ_CONFIRMED ; /* clear all CONFIRM flags */
      graphBoxDraw (box, -1, WHITE) ;
    }
}

static void intronConfirmEST (int box) 
{ intronConfirm (box, SEQ_CONFIRM_EST) ; }

static void intronConfirmHomol (int box)
{ intronConfirm (box, SEQ_CONFIRM_HOMOL) ; }

static void intronConfirmCDNA (int box)
{ intronConfirm (box, SEQ_CONFIRM_CDNA) ; }

static void intronConfirmUTR (int box)
{ intronConfirm (box, SEQ_CONFIRM_UTR) ; }

static void intronUnconfirm (int box)
{ intronConfirm (box, 0) ; }

static MENUOPT intronMenu[] = {
  { (GraphFunc)intronConfirmEST, "Confirm by EST" },
  { (GraphFunc)intronConfirmHomol, "Confirm by Homol" },
  { (GraphFunc)intronConfirmCDNA, "Confirm by CDNA" },
  { (GraphFunc)intronConfirmUTR, "Confirm in UTR" },
  { (GraphFunc)intronUnconfirm, "Unconfirm" },
  { 0, 0 }
} ;
#endif

/*****************************/

static void fMapFuseCosmids (int box)
{
  SEG *seg ;
  OBJ Cosmid1 = 0, Cosmid2 = 0 ;
  KEY superLink = 0, cosmid1 = 0, cosmid2 = 0, cosmid3 = 0, link12, link23 = 0, link13,  dna1Key = 0, dna2Key = 0 ;
  KEYSET links1 = 0, links2 = 0, links12 = 0, links22 = 0 ;
  int i, overlap1, overlap2, dnaLength1 = 0, dnaLength2 = 0, dnaLengthNew = 0 ;
  Array dna1 = 0, dna2 = 0 ;
  char *cp, *cq ;
  FMAPLOOKGET ("fMapFuseCosmids") ;
  
  if (box >= arrayMax (look->boxIndex) || ! (seg = BOXSEG (box)))
    { messerror ("Problem picking box in fMapFuseCosmids") ; 
      return ; 
    }

  cosmid1 = seg->key ;
  if (!cosmid1)
    { messerror ("Can't find cosmid to be fused") ;
      return ;
    }
  
  cosmid2 = keyGetKey (cosmid1, str2tag ("Overlap_right")) ;
  if (!cosmid2)
    return ;
  cosmid3 = keyGetKey (cosmid2, str2tag ("Overlap_right")) ;
  
  if (!messQuery (messprintf ("Are you sure you want to fuse %s %s", name (cosmid1), name (cosmid2))))
    return ;
  if (!checkWriteAccess ())
    return ;

  if ((Cosmid1 = bsCreate (cosmid1)))
    {
      if (bsGetKey (Cosmid1, _DNA, &dna1Key))
	bsGetData (Cosmid1, _bsRight, _Int, &dnaLength1) ;
      else
	{
	  KEY dummy ; int u1, u2 ;
	  if (bsGetKey (Cosmid1,str2tag ("IntMap"), &dummy) &&
	      bsGetData (Cosmid1, _bsRight, _Int, &u1) &&
	      bsGetData (Cosmid1, _bsRight, _Int, &u2))
	    {
	      if (u1 < u2) dnaLength1 = u2 - u1 + 1 ;
	      else dnaLength1 = u1 - u2 + 1 ;
	    }
	}

      if (bsGetKey (Cosmid1, str2tag ("Overlap_right"), &cosmid2))
	bsGetData (Cosmid1, _bsRight, _Int, &overlap1) ;
      bsGetKey (Cosmid1, _Source, &superLink) ;
      bsDestroy (Cosmid1) ;
    }
  if ((Cosmid2 = bsCreate (cosmid2)))
    {
      if (bsGetKey (Cosmid2, _DNA, &dna2Key))
	bsGetData (Cosmid2, _bsRight, _Int, &dnaLength2) ;
      else
	{
	  KEY dummy ; int u1, u2 ;
	  if (bsGetKey (Cosmid2,str2tag ("IntMap"), &dummy) &&
	      bsGetData (Cosmid2, _bsRight, _Int, &u1) &&
	      bsGetData (Cosmid2, _bsRight, _Int, &u2))
	    {
	      if (u1 < u2) dnaLength2 = u2 - u1 + 1 ;
	      else dnaLength2 = u1 - u2 + 1 ;
	    }
	}

      if (bsGetKey (Cosmid2, str2tag ("Overlap_right"), &cosmid3))
	bsGetData (Cosmid2, _bsRight, _Int, &overlap2) ;
      bsDestroy (Cosmid2) ;
      dnaLengthNew = dnaLength2 + overlap1 - 1 ;
    }
  
  if (!dnaLength1 || !dnaLength2 || dna1Key || dna2Key)
    { /* they must both be virtual and have intmap  or both be real */
      if (!dna1Key || !dna2Key || overlap1 > dnaLength1 + 1)
	{
	  messerror ("cannot fuse accross a dna gap, sorry") ;
	  goto abort ;
	}
      
      dna1 = dnaGet (dna1Key) ;
      dna2 = dnaGet (dna2Key) ;
      if (!dna1 || !dna2 ||
	  arrayMax (dna1) != dnaLength1 ||
	  arrayMax (dna2) != dnaLength2)
	{
	  messerror ("cannot correctly read the dna in fMapFuseCosmids") ;
	  goto abort ;
	}
      array (dna1, dnaLength2 + overlap1 - 1, char) = 0 ; /* make room including terminal zero */
      dnaLengthNew = arrayMax (dna1) = dnaLength2 + overlap1 - 1 ; /* remove terminal zero */
      cp = arrp (dna1, overlap1 - 1, char) ;
      cq = arrp (dna2, 0, char) ;
      for (i = overlap1 ; i <= dnaLength1 ; i++, cp++, cq++)
	if (! (*cp & *cp))
	  {
	    messerror ("missmatch in overlapping cosmid dna, cannot fuse, sorry") ;
	    goto abort ;
	  }
      /* all verifs ok, i fuse the dna */
      cp = arrp (dna1, overlap1 - 1, char) ;
      cq = arrp (dna2, 0, char) ;
      for (i = 0 ; i < dnaLength2 ; i++, cp++, cq++)
	*cp = *cq ;
    }
  /* now i update cosmid1 */
  if ((Cosmid1 = bsUpdate (cosmid1)))
    {
      if (dna1)
	{
	  i = arrayMax (dna1) ;
	  if (bsGetKey (Cosmid1, _DNA, &dna1Key))
	    bsAddData (Cosmid1, _bsRight, _Int, &i) ;
	}
      if (bsFindTag (Cosmid1, str2tag ("Overlap_right")))
	bsRemove (Cosmid1) ;
      if (cosmid3)
	{
	  i = overlap1 + overlap2 - 1 ;
	  bsAddKey (Cosmid1, str2tag ("Overlap_right"), cosmid3) ;
	  bsAddData (Cosmid1, _bsRight, _Int, &i) ;
	}
      bsSave (Cosmid1) ;
    }
  /* change its coordinates in superlink */
  {
    OBJ SuperLink = 0 ;
    int a1, a2 ;

    if ((SuperLink = bsUpdate (superLink)))
      {
	if (bsFindKey (SuperLink, _Subsequence, cosmid1) &&
	    bsGetData (SuperLink, _bsRight, _Int, &a1) &&
	    bsGetData (SuperLink, _bsRight, _Int, &a2))
	  {
	    if (a1 < a2)
	      a2 = a1 + dnaLengthNew - 1 ;
	    else
	      a2 = a1 - dnaLengthNew + 1 ;
	    if (bsFindKey (SuperLink, _Subsequence, cosmid1) &&
		bsGetData (SuperLink, _bsRight, _Int, &a1))
	      bsAddData (SuperLink, _bsRight, _Int, &a2) ;
	  }
	bsSave (SuperLink) ;
      }
  }

  links1 = queryKey (cosmid1, ">In_junction") ;
  links2 = queryKey (cosmid2, ">In_junction") ;
  /* deal with the junctions */
  
  links12 = keySetAND (links1, links2) ;
  links22 = keySetMINUS (links2, links1) ;
  link12 = link23 = 0 ;
  link12 = keySet (links12, 0) ;
  link23 = keySet (links22, 0) ;

  /* remove the genomic tag, hence the parts of links12 */
  if ((Cosmid2 = bsUpdate (link12)))
    {
      if (bsFindTag (Cosmid2, str2tag ("Parts")))
	bsRemove (Cosmid2) ;
      if (bsFindTag (Cosmid2, str2tag ("Junction")))
	bsRemove (Cosmid2) ;
      bsSave (Cosmid2) ;
    }

  if (link23 && (Cosmid2 = bsUpdate (link23)))
    {
      if (bsFindTag (Cosmid2, str2tag ("Parts")))
	bsRemove (Cosmid2) ;
      if (bsFindTag (Cosmid2, str2tag ("Junction")))
	bsRemove (Cosmid2) ;
      bsSave (Cosmid2) ;
    }

  if ((Cosmid2 = bsUpdate (cosmid2)))
    {
      if (bsFindTag (Cosmid2, str2tag ("Genomic")))
	bsRemove (Cosmid2) ;
      bsSave (Cosmid2) ;
    }

  if (cosmid3)
    {
      OBJ Link13 = 0 ;
      OBJ SuperLink = 0 ;
      int u1, u2 ;

      lexaddkey (messprintf ("%s_%s", name (cosmid1), name (cosmid3)), &link13, _VSequence) ;
      if ((Link13 = bsUpdate (link13)))
	{
	  bsAddKey (Link13, str2tag ("Parts"), cosmid1) ;
	  bsAddKey (Link13, str2tag ("Parts"), cosmid3) ;
	  bsSave (Link13) ;
	}
      if ((SuperLink = bsUpdate (superLink)))
	{
	  if (bsFindKey (SuperLink, _Subsequence, cosmid3) &&
	      bsGetData (SuperLink, _bsRight, _Int, &u1) &&
	      bsGetData (SuperLink, _bsRight, _Int, &u2) &&
	      bsFindKey (SuperLink, _Subsequence, cosmid1) &&
	      bsGetData (SuperLink, _bsRight, _Int, &u1))
	    {
	      bsAddKey (SuperLink, _Subsequence, link13) ;
	      bsAddData (SuperLink, _bsRight, _Int, &u1) ;
	      bsAddData (SuperLink, _bsRight, _Int, &u2) ;
	    }
	  bsSave (SuperLink) ;
	}
    }
  /* restore the tgs */
  {
    OBJ Link13 = 0, Link  = 0 ;
    int x ;
    Array units = arrayCreate (30, BSunit) ;
    BSunit *uu ;

    if ((Link = bsCreate (link12)))
      {
	if (bsGetArray (Link, str2tag ("Transcribed_gene"), units, 3))
	  {
	    if ((Link13 = bsUpdate (cosmid1)))
	      {
		for (i = 0 ; i < arrayMax (units) ; i += 3)
		  {
		    uu = arrp (units, i, BSunit) ;
		    bsAddKey (Link13,  str2tag ("Transcribed_gene"), uu[0].k) ;
		    x = uu[1].i ;
		    bsAddData (Link13, _bsRight, _Int, &x) ; 
		    x = uu[2].i ;
		    bsAddData (Link13, _bsRight, _Int, &x) ;
		  }
		bsSave (Link13) ;
	      }
	  }
	bsDestroy (Link) ;
	if ((Link = bsUpdate (link12)))
	  {
	    if (bsFindTag (Link, str2tag ("Transcribed_gene")))
	      bsRemove (Link) ;
	    bsSave (Link) ;
	  }
      }
    if ((Link = bsCreate (cosmid2)))
      {
	if (bsGetArray (Link, str2tag ("Transcribed_gene"), units, 3))
	  {
	    if ((Link13 = bsUpdate (cosmid1)))
	      {
		for (i = 0 ; i < arrayMax (units) ; i += 3)
		  {
		    uu = arrp (units, i, BSunit) ;
		    bsAddKey (Link13,  str2tag ("Transcribed_gene"), uu[0].k) ;
		    x = uu[1].i + overlap1 - 1 ;
		    bsAddData (Link13, _bsRight, _Int, &x) ; 
		    x = uu[2].i + overlap1 - 1;
		    bsAddData (Link13, _bsRight, _Int, &x) ;
		  }
		bsSave (Link13) ;
	      }
	  }
	bsDestroy (Link) ; 
	if ((Link = bsUpdate (cosmid2)))
	  {
	    if (bsFindTag (Link, str2tag ("Transcribed_gene")))
	      bsRemove (Link) ;
	    bsSave (Link) ;
	  }
      }
    if (link23 && (Link = bsCreate (link23)))
      {
	if (bsGetArray (Link, str2tag ("Transcribed_gene"), units, 3))
	  {
	    if ((Link13 = bsUpdate (link13)))
	      {
		for (i = 0 ; i < arrayMax (units) ; i += 3)
		  {
		    uu = arrp (units, i, BSunit) ;
		    bsAddKey (Link13,  str2tag ("Transcribed_gene"), uu[0].k) ;
		    x = uu[1].i + overlap1 - 1 ;
		    bsAddData (Link13, _bsRight, _Int, &x) ; 
		    x = uu[2].i + overlap1 - 1;
		    bsAddData (Link13, _bsRight, _Int, &x) ;
		  }
		bsSave (Link13) ;
	      }
	  }
	bsDestroy (Link) ;
	if ((Link = bsUpdate (link23)))
	  {
	    if (bsFindTag (Link, str2tag ("Transcribed_gene")))
	      bsRemove (Link) ;
	    bsSave (Link) ;
	  }
      }
    arrayDestroy (units) ;
  }

  /* restore the genes */
  {
    OBJ Link13 = 0, Link  = 0 ;
    int x ;
    Array units = arrayCreate (30, BSunit) ;
    BSunit *uu ;

    if ((Link = bsCreate (link12)))
      {
	if (bsGetArray (Link, str2tag ("Genes"), units, 3))
	  {
	    if ((Link13 = bsUpdate (cosmid1)))
	      {
		for (i = 0 ; i < arrayMax (units) ; i += 3)
		  {
		    uu = arrp (units, i, BSunit) ;
		    bsAddKey (Link13,  str2tag ("Genes"), uu[0].k) ;
		    x = uu[1].i ;
		    bsAddData (Link13, _bsRight, _Int, &x) ; 
		    x = uu[2].i ;
		    bsAddData (Link13, _bsRight, _Int, &x) ;
		  }
		bsSave (Link13) ;
	      }
	  }
	bsDestroy (Link) ;
	if ((Link = bsUpdate (link12)))
	  {
	    if (bsFindTag (Link, str2tag ("Genes")))
	      bsRemove (Link) ;
	    bsSave (Link) ;
	  }
      }
    if ((Link = bsCreate (cosmid2)))
      {
	if (bsGetArray (Link, str2tag ("Genes"), units, 3))
	  {
	    if ((Link13 = bsUpdate (cosmid1)))
	      {
		for (i = 0 ; i < arrayMax (units) ; i += 3)
		  {
		    uu = arrp (units, i, BSunit) ;
		    bsAddKey (Link13,  str2tag ("Genes"), uu[0].k) ;
		    x = uu[1].i + overlap1 - 1 ;
		    bsAddData (Link13, _bsRight, _Int, &x) ; 
		    x = uu[2].i + overlap1 - 1;
		    bsAddData (Link13, _bsRight, _Int, &x) ;
		  }
		bsSave (Link13) ;
	      }
	  }
	bsDestroy (Link) ; 
	if ((Link = bsUpdate (cosmid2)))
	  {
	    if (bsFindTag (Link, str2tag ("Genes")))
	      bsRemove (Link) ;
	    bsSave (Link) ;
	  }
      }
    if (link23 && (Link = bsCreate (link23)))
      {
	if (bsGetArray (Link, str2tag ("Genes"), units, 3))
	  {
	    if ((Link13 = bsUpdate (link13)))
	      {
		for (i = 0 ; i < arrayMax (units) ; i += 3)
		  {
		    uu = arrp (units, i, BSunit) ;
		    bsAddKey (Link13,  str2tag ("Genes"), uu[0].k) ;
		    x = uu[1].i + overlap1 - 1 ;
		    bsAddData (Link13, _bsRight, _Int, &x) ; 
		    x = uu[2].i + overlap1 - 1;
		    bsAddData (Link13, _bsRight, _Int, &x) ;
		  }
		bsSave (Link13) ;
	      }
	  }
	bsDestroy (Link) ;
	if ((Link = bsUpdate (link23)))
	  {
	    if (bsFindTag (Link, str2tag ("Genes")))
	      bsRemove (Link) ;
	    bsSave (Link) ;
	  }
      }
    arrayDestroy (units) ;
  }

  /* intmap */
  if ((Cosmid1 = bsUpdate (cosmid1)))
    {
      int a1, a2 ;
      KEY map = 0 ;
      OBJ Link13 = 0 ;

      a1 = a2 = 0 ;
      if (bsGetKey (Cosmid1, str2tag ("IntMap"), &map) &
	  bsGetData (Cosmid1, _bsRight, _Int, &a1) &&
	  bsGetData (Cosmid1, _bsRight, _Int, &a2))
	{
	  if (a1 < a2)
	    {
	      a2 = a1 + dnaLengthNew -  1;
	      if (bsGetKey (Cosmid1, str2tag ("IntMap"), &map) &
		  bsGetData (Cosmid1, _bsRight, _Int, &a1))
		bsAddData (Cosmid1, _bsRight, _Int, &a2) ;
	    }
	  else
	    {
	      a2 = a1 - dnaLengthNew + 1;
	      if (bsGetKey (Cosmid1, str2tag ("IntMap"), &map) &
		  bsGetData (Cosmid1, _bsRight, _Int, &a1))
		bsAddData (Cosmid1, _bsRight, _Int, &a2) ;
	    }
	}
      bsSave (Cosmid1) ;
      if (map && a1 != a2 && link13 && (Link13 = bsUpdate (link23)))
	{
	  bsAddKey (Link13, str2tag ("IntMap"), map) ;
	  bsAddData (Link13, _bsRight, _Int, &a1) ;
	  if (a1 < a2)
	    a2 = a1 + overlap1 + dnaLength2 - 1 ;
	  else
	    a2 = a1 - (overlap1 + dnaLength2 - 1) ;
	  bsAddData (Link13, _bsRight, _Int, &a2) ;
	  bsSave (Link13) ;
	}
    }
  if (dna1) { dnaStoreDestroy (dna1Key, dna1) ; dna1 = 0 ; }
  arrayDestroy (dna2) ;
  fMapPleaseRecompute (look) ;
  fMapDraw (look, TRUE) ;


 abort:
  keySetDestroy (links1) ;
  keySetDestroy (links2) ;
  keySetDestroy (links12) ;
  keySetDestroy (links22) ;
  arrayDestroy (dna1) ;
  arrayDestroy (dna2) ;
}

/*****************************/

static MENUOPT fMapClonesMenu[] = { 
  { (GraphFunc)fMapFuseCosmids, "Fuse to Next Cosmid" },
  { 0, 0 }
} ;

/*****************************/

void fMapShowSequence (LOOK look, float *offset)
{
  void  *xx ;
  int	i, box, x ;
  float y1, y2, yb ;
  SEG	*seg ;
  Associator seq2x = assCreate () ;
  SEQINFO *sinf ;
  BoxCol *bc ;
  BOOL isData = FALSE, firstError = TRUE ;
  
  bc = bcFromName (look, look->map->activeColName) ;
  if (!bc || !bcTestMag (bc, look->map->mag)) return ;
  bc->offset = *offset ;

  if (bc->autoWidth) bc->width = 1.25 ;
  if (!bc->colour) bc->colour = BLUE ;

  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs,i,SEG) ;
      if ((seg->type & 0x01) == bc->isDown)
	continue ;
      if (seg->x1 > look->max || seg->x2 < look->min)
	continue ;
      if (look->view &&
	  keyFindTag (seg->source, str2tag ("Model_of_gene")) &&
	  keyFindTag (seg->source, str2tag ("Predicted_mRNA")))
	continue ;
      y1 = MAP2GRAPH (look->map,seg->x1 - .5) ;
      y2 = MAP2GRAPH (look->map,seg->x2+1 - .5) ;	/* to cover full base */
      if (y1 < 2 + topMargin) y1 = 2 + topMargin ;
      if (y2 < 2 + topMargin) y2 = 2 + topMargin ;
      switch (seg->type)
	{ 
	case SEQUENCE: case SEQUENCE_UP:
	  sinf = arrp (look->seqInfo, seg->data.i, SEQINFO) ;
	  if (!seg->parent || sinf->method != bc->method ||
	      ! (sinf->flags & (SEQ_VISIBLE | SEQ_EXONS)))
	    break ;
	  x = 0 ; yb = seg->x1 ; /* bump in seg coords, since increasing */
	  isData = TRUE ;
	  if (bc->bump) 
	    bumpItem (bc->bump, 1, (seg->x2-seg->x1+1), &x, &yb) ;
	  assInsert (seq2x, assVoid (seg->parent), assVoid (x+1)) ;
	  if (bc->isShowText)
	    addVisibleInfo (look, i) ;
	  if (sinf->flags & SEQ_EXONS) /* will draw exons later */
	    break ;
	  box = graphBoxStart () ;
	  array (look->boxIndex,box,int) = i ;
	  graphRectangle (*offset + (x + 0.25)*bc->width, y1, 
			  *offset + (x + 0.75)*bc->width, y2) ;
	  graphBoxEnd () ;
	  fMapBoxInfo (look, box, seg) ;
	  graphBoxDraw (box, bc->colour, WHITE) ;
	  break ;
	case EXON: case EXON_UP:
	case EXON_GAP: case EXON_GAP_UP:
	  sinf = arrp (look->seqInfo, seg->data.i, SEQINFO) ;
	  if (sinf->method != bc->method) break ;
	  if (assFind (seq2x, assVoid (seg->parent), &xx))
	    x = assInt (xx) - 1 ;
	  else
	    { if (firstError) messerror ("Can't find offset for %s in showGraphSequence",
			 name (seg->parent)) ;
	      firstError = FALSE ;
	      break ;
	    }
	  box = graphBoxStart () ;
	  array (look->boxIndex,box,int) = i ;
	  switch (seg->type)
	    {
	    case EXON: case EXON_UP:
	      graphRectangle (*offset + (x + 0.1)*bc->width, y1, 
			      *offset + (x + 0.9)*bc->width, y2) ;
	      break ;	
	    case EXON_GAP: case EXON_GAP_UP: 
	      graphLine (*offset + (x + 0.5)*bc->width, y1, 
			 *offset + (x + 0.5)*bc->width, y2) ;
	      break ;
	    default:
	      break ;
	    }
	  graphBoxEnd () ;
	  fMapBoxInfo (look, box, seg) ;
	  graphBoxDraw (box, bc->colour, WHITE) ;
	  if (isGifDisplay)
	    {
	      KEY pg = seg->parent ;
	      KEY locus = keyGetKey (pg, str2tag ("GeneId_pg")) ;
	      KEY locusLink = 0 ;
	      KEY xm = keyGetKey (pg, str2tag ("NM_id")) ;
	      char *pgNam = name (pg) ;

	      if (*pgNam == '_') pgNam++ ;
	      if (!locus) locus = keyGetKey (pg, _LocusId) ;
	      locusLink = locus ?  keyGetKey (locus, str2tag ("LocusLink")) : locus ;
	      if (!locusLink) locusLink = locus ;
	      if (xm && !strncmp (name (xm),"NM_", 3))
		graphBubbleInfo (box, locus, "URL_LL", messprintf ("NCBI model %s : %s", name (xm), pgNam)) ;
	      else if (xm && !strncmp (name (xm),"XM_", 3))
		graphBubbleInfo (box, locus, "URL_LL", messprintf ("NCBI predicted  model %s %s", name (xm), pgNam)) ;
	      else if (xm)
		graphBubbleInfo (box, xm, "URL_GB", messprintf ("NCBI model %s", pgNam)) ; 
	      else
		graphBubbleInfo (box, seg->parent, "Predicted_gene", "Model") ;
	    }
	  /***************************************/
	  /***** dynamically change geneMenu *****/
	  geneMenu[0].text = "Show translation";
	  if (look->translateGene == seg->parent)
	    /* we have that gene already selected for translation */
	    {
	      if ((seg->type == EXON_UP && mapColSetByName ("Up Gene Translation", -1))
		  ||
		  (seg->type == EXON && mapColSetByName ("Down Gene Translation", -1)))
		{
		  /* the column to show the translation for the already 
		     selected gene is already switched on, so we change the 
		     label to "Hide translation" */
		  geneMenu[0].text = "Hide translation";
		}
	    }

	  graphBoxMenu (box, geneMenu) ;
	  break ;
	case INTRON: case INTRON_UP:
	  sinf = arrp (look->seqInfo, seg->data.i, SEQINFO) ;
	  if (sinf->method != bc->method) break ;
	  if (assFind (seq2x, assVoid (seg->parent), &xx))
	    x = assInt (xx) - 1 ;
	  else
	    { if (firstError) messerror ("Can't find offset for %s in showGraphSequence",
			 name (seg->parent)) ;  
	      firstError = FALSE ;
	      break ;
	    }
	  box = graphBoxStart () ;
	  array (look->boxIndex,box,int) = i ;
	  y2 = MAP2GRAPH (look->map,seg->x2+1) ;	/* to match existing box */
	  graphLine (*offset + (x + 0.5)*bc->width, y1, 
		     *offset + (x + 0.9)*bc->width, y1) ; /* for clean appearance */
	  graphLine (*offset + (x + 0.5)*bc->width, y1, 
		     *offset + (x + 0.9)*bc->width, 0.5* (y1+y2)) ;
	  graphLine (*offset + (x + 0.5)*bc->width, y2, 
		     *offset + (x + 0.9)*bc->width, 0.5* (y1+y2)) ;
	  graphBoxEnd () ;
	  fMapBoxInfo (look, box, seg) ;
	  if (sinf->flags & SEQ_CONFIRMED)
	    graphBoxDraw (box, bc->colour, sinf->flags & SEQ_CONFIRM_UTR ? YELLOW : GREEN) ;
	  else
	    graphBoxDraw (box, bc->colour, WHITE) ;
#ifndef ACEMBLY
	  if (seg->type == INTRON)
	    graphBoxMenu (box, intronMenu) ;
#endif
	  if (isGifDisplay)
	    {
	      KEY pg = seg->parent ;
	      KEY locus = keyGetKey (pg, _LocusId) ;
	      KEY locusLink = 0 && locus ?  keyGetKey (locus, str2tag ("LocusLink")) : locus ;
	      KEY xm = keyGetKey (pg, str2tag ("NM_id")) ;

	      if (xm && !strncmp (name (xm),"NM_", 3))
		graphBubbleInfo (box, locusLink, "URL_LL", messprintf ("NCBI model %s %s", name (xm), name (pg))) ;
	      else if (xm && !strncmp (name (xm),"XM_", 3))
		graphBubbleInfo (box, locusLink, "URL_LL", messprintf ("NCBI predicted  model %s %s", name (xm), name (pg))) ;
	      else if (xm)
		graphBubbleInfo (box, xm, "URL_GB", messprintf ("NCBI model %s", name (pg))) ; 
	      else
		graphBubbleInfo (box, seg->parent, "Predicted_gene", "Model") ;
	    }
	  break ;
	default:
	  break ;
	}
    }

  if (bc->bump) 
    { *offset += bc->width*bumpMax (bc->bump) ;
      bumpDestroy (bc->bump) ; bc->bump = 0 ; 
    }
  else if (isData)
    *offset += bc->width ;
    
  assDestroy (seq2x) ;
}

void fMapShowSoloConfirmed (LOOK look, float *offset)
{
  int	i, box, x ;
  float y1, y2 ;
  SEG	*seg ;
  SEQINFO *sinf ;
  BUMP bump = bumpCreate (20, 0) ;
  BOOL isDown = (*look->map->activeColName == '-') ? FALSE : TRUE ;

  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs,i,SEG) ;
      if (seg->x1 > look->max || seg->x2 < look->min)
	continue ;
      if (! ((seg->type == INTRON && isDown) ||
	    (seg->type == INTRON_UP && !isDown)))
	continue ;
      sinf = arrp (look->seqInfo, seg->data.i, SEQINFO) ;
      if (seg->parent || ! (sinf->flags & SEQ_CONFIRMED))
	continue ;
      y1 = MAP2GRAPH (look->map,seg->x1) ;
      y2 = MAP2GRAPH (look->map,seg->x2+1) ;	/* to cover full base */
      if (y1 < 2 + topMargin) y1 = 2 + topMargin ;
      if (y2 < 2 + topMargin) y2 = 2 + topMargin ;
      x = 1 ; bumpItem (bump, 1, (y2-y1), &x, &y1) ;
      box = graphBoxStart () ;
      array (look->boxIndex,box,int) = i ;
      graphLine (*offset + x, y1, *offset + x + 0.5, 0.5* (y1+y2)) ;
      graphLine (*offset + x, y2, *offset + x + 0.5, 0.5* (y1+y2)) ;
      graphBoxEnd () ;
      graphBoxDraw (box, BLACK, GREEN) ;
      fMapBoxInfo (look, box, seg) ;
    }
  *offset += bumpMax (bump) ;
  bumpDestroy (bump) ;
}

void fMapShowEmblFeatures (LOOK look, float *offset)
{
  int	i, box, x ;
  float y1, y2 ;
  SEG	*seg ;
  BUMP bump = bumpCreate (20, 0) ;

  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs,i,SEG) ;
      if (seg->x1 > look->max || seg->x2 < look->min)
	continue ;
      y1 = MAP2GRAPH (look->map,seg->x1) ;
      y2 = MAP2GRAPH (look->map,seg->x2+1) ;	/* to cover full base */
      if (y1 < 2 + topMargin) y1 = 2 + topMargin ;
      if (y2 < 2 + topMargin) y2 = 2 + topMargin ;
      switch (seg->type)
	{ 
	case EMBL_FEATURE: case EMBL_FEATURE_UP:
	  x = 0 ; bumpItem (bump, 1, (y2-y1), &x, &y1) ;
	  box = graphBoxStart () ;
	  array (look->boxIndex,box,int) = i ;
	  graphRectangle (*offset + x + 0.25, y1, 
			  *offset + x + 0.75, y2) ;
	  graphBoxEnd () ;
	  fMapBoxInfo (look, box, seg) ;
	  graphBoxDraw (box, DARKRED, WHITE) ;
	  array (look->visibleIndices,arrayMax (look->visibleIndices),int) = i ;
	  break ;
	default:
	  break ;
	}
    }
  *offset += bumpMax (bump) ;
  bumpDestroy (bump) ;
}

/****************************************************************/
/******************** text features *****************************/

static BOOL showVisibleTagNames = FALSE ; /* suppress tag names in "Visible" column */

static int bumpOff ;

static void bumpWrite (BUMP bump, const char *text, int x, float y)
{
  bumpItem (bump, strlen (text)+1, 1, &x, &y) ;
  graphText (text, bumpOff+x, y-0.5) ;
}

static BOOL isReversed ;
static LOOK x1Look ;

static int x1Order (const void *a, const void *b)
{
  int diff ;

  if (isReversed)
    diff = arrp (x1Look->segs, * (const  int*)b, SEG)->x2  -
      arrp (x1Look->segs, * (const int*)a, SEG)->x2  ;
  else
    diff = arrp (x1Look->segs, * (const int*)a, SEG)->x1  -
      arrp (x1Look->segs, * (const int*)b, SEG)->x1  ;

  if (diff > 0) 
    return 1 ;
  else if (diff < 0)
    return -1 ;
  else
    return 0 ;
}

static void addVisibleInfo (LOOK look, int i) /* i is the index of a seg */
{
  KEY key = arrp (look->segs,i,SEG)->key ;
  register int j ;
  register SEG *seg = 0 ;
  static int jmin, jmax ;	/* for efficiency cache start/end of VISIBLE segs */
				/* relies on segs not being reordered during drawing */

  if (!arrayMax (look->visibleIndices))
    fMapFindSegBounds (look, VISIBLE, &jmin, &jmax) ;
    
				/* first add the item itself */
  array (look->visibleIndices,arrayMax (look->visibleIndices),int) = i ;

				/* then related VISIBLEs */
  if (jmin < jmax) seg = arrp (look->segs,jmin,SEG) ;
  for (j = jmin ; j < jmax ; ++j, ++seg)
    if (seg->parent == key)
      array (look->visibleIndices,arrayMax (look->visibleIndices),int) = j ;
}

void fMapShowText (LOOK look, float *offset)
{
  int	i, j, box ;
  Stack textStack = 0 ;
  float y ;
  SEG	*seg ;
  HOMOLINFO *hinf ;
  BUMP  bump ;

  if (!fmapView (look->view, str2tag ("Fmap_Text_Features")))
    return ;

  textStack = stackCreate (100) ;
  bump = bumpCreate (40, 0) ;
  bumpOff = *offset ;

  isReversed = ((look->flag & FLAG_REVERSE) > 0) ;
  x1Look = look ;
  arraySort (look->visibleIndices, x1Order) ;

  graphColor (BLACK) ;

  for (j = 0 ; j < arrayMax (look->visibleIndices) ; ++j)
    { i = arr (look->visibleIndices, j, int) ;
      seg = arrp (look->segs, i, SEG) ;
      if (isReversed)
	y = MAP2GRAPH (look->map,seg->x2) ;
      else
	y = MAP2GRAPH (look->map,seg->x1) ;
      if (y < topMargin + 2)
	y = topMargin + 2 ;
      box = 0 ;
      switch (seg->type)
	{ 
	case SEQUENCE: case SEQUENCE_UP:
	  box = graphBoxStart () ;
	  bumpWrite (bump,name (seg->key), 0, y) ; 
	  graphBubbleInfo (box, seg->parent,  className (seg->parent), "") ;
	  graphBoxEnd () ;
	  break ;

	case HOMOL: case HOMOL_UP: 
	  box = graphBoxStart () ;
	  bumpWrite (bump,name (seg->key), 4, y) ;
	  graphBoxEnd () ; 
	  hinf = arrp (look->homolInfo, seg->data.i, HOMOLINFO) ;
	  graphBubbleInfo (box, seg->parent, className (seg->parent), hinf->method ? name (hinf->method) : "Homol") ;
	  break ;

	case FEATURE: case FEATURE_UP:
	  box = graphBoxStart () ;
	  bumpWrite (bump, dictName (look->featDict, -seg->parent), 4, y) ;
	  graphBoxEnd () ;
	  break ;
	  
	case VISIBLE:
	  box = graphBoxStart () ;
	  stackClear (textStack) ;
	  if (showVisibleTagNames && iskey (seg->data.k))  /* mieg: prevents (NULL KEY)) */
	    { catText (textStack, name (seg->data.k)) ;
	      catText (textStack, ":") ;
	    }
	  /* there is a problem if the model says: Visible foo Text
               this is the origin of the NULL KEY in the displays */
	  if (iskey (seg->key))   /* mieg: hack to prevents (NULL KEY)) */
	    catText (textStack, name (seg->key)) ;
	  bumpWrite (bump,stackText (textStack,0),10,y) ; 

	  graphBoxEnd () ;

	  if (iskey (seg->key)) 
	    {
	      graphBoxInfo (box, seg->key,
			    strnew (messprintf ("%s:%s",
						iskey (seg->data.k) ? name (seg->data.k) : "" , 
						name (seg->key)),
				    look->segsHandle)) ;
	      graphBubbleInfo (box, seg->parent, className (seg->parent), "") ;
	    }

	  break ;
	  
	case ALLELE: case ALLELE_UP:
	case EMBL_FEATURE: case EMBL_FEATURE_UP:
	  box = graphBoxStart () ;
	  stackClear (textStack) ;
	  if (iskey (seg->key)) 
	    catText (textStack, name (seg->key)) ;
	  if (seg->data.s)
	    { catText (textStack, ":") ;
	      catText (textStack, seg->data.s) ;
	    }
	  bumpWrite (bump,stackText (textStack,0),10,y) ; 
	  graphBoxEnd () ;
	  break ;

	default:
	  break ;
	}
      if (box)
	{ array (look->boxIndex,box,int) = i ;
	  fMapBoxInfo (look, box, seg) ;
	}

    }

  *offset += bumpMax (bump) ;
  bumpDestroy (bump) ;
  stackDestroy (textStack) ;
}

/*************************************************************/
/******** external stuff: efetch, blxview, trnafold ***********/

static void fMapFetch (int box)
{
  static KEY _Database = 0 ;
  SEG *seg ;
  KEY lib ;
  OBJ obj = 0, libObj = 0 ;
  char *id, *pref, *suff, *url = 0 ;
  FMAPLOOKGET ("fMapFetch") ;

  if (!_Database)
    lexaddkey ("Database", &_Database, 0) ;

  if (! (seg = BOXSEG (box)))
    return ;

  if ((obj = bsCreate (seg->key)) &&
      bsGetKey (obj, _Database, &lib) &&
      bsGetData (obj, _bsRight, _Text, &id) &&
      (libObj = bsCreate (lib)))
    {
      if (bsGetData (libObj, _Arg1_URL_prefix, _Text, &pref))
	{
	  if (bsGetData (libObj, _Arg1_URL_suffix, _Text, &suff))
	    url = messprintf ("%s%s%s", pref, id, suff) ;
	  else
	    url = messprintf ("%s%s", pref, id) ;
	}
      else if (bsGetData (obj, _bsRight, _Text, &id) &&
	       bsGetData (libObj, _Arg2_URL_prefix, _Text, &pref))
	{
	  if (bsGetData (libObj, _Arg2_URL_suffix, _Text, &suff))
	    url = messprintf ("%s%s%s", pref, id, suff) ;
	  else
	    url = messprintf ("%s%s", pref, id) ;
	}
    }
  bsDestroy (obj) ; bsDestroy (libObj) ; /* OK if 0 */

  if (!url)
    url = messprintf ("http://www.sanger.ac.uk/cgi-bin/seq-query?%s", 
		      name (seg->key)) ;
  
  graphWebBrowser (url) ;
}

#define ACEDB			/* for this include of blxview */
#include "blxview.h"
#include "dna.h"
#include "client.h"

static void  showAlignAnticipate (KEYSET ks)
{
  KEYSET ks2 ;

  if (!externalServer)
    return ;
  keySetSort (ks) ;
  keySetCompress (ks) ;
  oldSetForServer = ks ;
  externalServer (-1, 0, 0, TRUE) ; /* get these sequences */
  ks2 = queryLocalParametrized (ks,"{FOLLOW DNA} $| {FOLLOW Peptide}",0) ;
  oldSetForServer = ks2 ;
  externalServer (-1, 0, 0, TRUE) ;
  oldSetForServer = 0 ;
  keySetDestroy (ks2) ;
}

static void doShowAlign (int x1, int x2, char *opts, KEY meth)
{
  SEG *seg ;
  int nks = 0, i, min, max, first, offset;
  char *seq, *cp, *cq ;  
  MSP msp1, *msp, *msp2, *prevmsp=0 ;
  KEY key ;
  HOMOLINFO *hinf ;
  unsigned int methodFlag ;
  Array a ;
  KEYSET ks = 0 ;
  OBJ obj = 0 ;
  FMAPLOOKGET ("fMapShowAlign") ;

  methodFlag = (*opts == 'X') ? METHOD_BLIXEM_X : METHOD_BLIXEM_N ;

  min = x1 - 20000 ;
  max = x2 + 20000 ;
  if (min < 0) min = 0 ;
  if (max >= look->length) max = look->length - 1 ;

  seq = (char*) messalloc (max-min+2) ;
  cp = seq ; cq = arrp (look->dna, min, char) ;
  for (i = min ; i++ <= max ; )
    *cp++ = dnaDecodeChar[ (int)*cq++] ;
  *cp = 0 ;

  offset = min - look->origin ;
  first = x1 - min - 30 ;
  if (first <= 0) first = 1 ;

  ks = keySetCreate () ; nks = 0 ;
  memset (&msp1, 0, sizeof (msp1)) ; /* zero whole structure */
  msp = &msp1 ;
  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
    { 
      seg = arrp (look->segs, i, SEG) ;

      if (seg->x2 <= min || seg->x1 >= max) continue ;

      if (seg->type & 0x1)	/* UP */
	{ msp->qstart = seg->x2 - min + 1 ;
	  msp->qend = seg->x1 - min + 1 ;
	}
      else			/* DOWN */
	{ msp->qstart = seg->x1 - min + 1 ;
	  msp->qend = seg->x2 - min + 1 ;
	}

      switch (seg->type)
	{
 	case HOMOL: case HOMOL_UP:
	  hinf = arrp (look->homolInfo, seg->data.i, HOMOLINFO) ;
	  if (meth && lexAliasOf (hinf->method) != lexAliasOf (meth))
	    break ;
	  if (! (method (0, hinf->method)->flags &
		methodFlag))
	    break ;
	  msp->key = seg->key ;
	  strncpy (msp->sname, name (seg->key), 19) ;
	  keySet (ks, nks++) = seg->key ;
	  if (seg->type & 0x1)
	    strcpy (msp->frame, messprintf (" (-%d)", 
			    1 + ((max-min+1) - msp->qstart) % 3)) ;
	  else
	    strcpy (msp->frame, messprintf (" (+%d)", 
			    1 + (msp->qstart-1) % 3)) ;
	  msp->score = hinf->score + 0.5 ; /* float to int convert */
	  msp->sstart = hinf->x1 ;
	  msp->send = hinf->x2 ;
	  break ;

	case CODING: case CODING_UP:
	  msp->sstart = (seg->data.i + 3) / 3 ;
	  msp->send = (seg->data.i + seg->x2 - seg->x1) / 3 ;
	  msp->score = -1 ;
	  msp->id = 100;
	  sprintf (msp->sname, "%sx", name (seg->parent)) ;
	  keySet (ks, nks++) = seg->parent ;
	  if (! (seg->type & 0x1))  /* x separates from Wormpep entry */
	    sprintf (msp->frame, " (+%d)", 
		     1 + (msp->qstart - 1 - seg->data.i)%3) ;
	  else
	    sprintf (msp->frame, " (-%d)", 
		     1 + ((max-min+1) - msp->qstart - seg->data.i)%3) ;
	  break ;

	case INTRON: case INTRON_UP:
	  msp->sstart = msp->send = msp->id = 0 ;
	  msp->score = -2 ;
	  msp->id = 100;
	  sprintf (msp->sname, "%si", name (seg->parent)) ;  /* i for intron */
	  keySet (ks, nks++) = seg->parent ;
	  if (! (seg->type & 0x1))
	    strcpy (msp->frame, " (+1)") ;
	  else
	    strcpy (msp->frame, " (-1)") ;
	  break ;
	default:
	  break ;
	}

      if (msp->score)		/* filled it */
	{ msp->next = (MSP*) messalloc (sizeof (MSP)) ;
          prevmsp = msp ;
	  msp = msp->next ;
	}
    }
				/* Free last MSP (empty) */
  if (prevmsp) 
    {
      messfree (prevmsp->next) ;
      prevmsp->next = 0;
    }

  if (externalServer)
    showAlignAnticipate (ks) ;
  keySetDestroy (ks) ;

  if (*opts == 'X')		/* Protein */
    for (msp = &msp1 ; msp ; msp = msp->next)
      { if (msp->qstart < 1)
	  { msp->sstart += (3 - msp->qstart) / 3 ;
	    msp->qstart += 3 * ((3 - msp->qstart) / 3) ;
	  }
	if (msp->qend > max-min+1)
	  { msp->send -= (msp->qend + 2 - (max-min)) / 3 ;
	    msp->qend -= 3 * ((msp->qend + 2 - (max-min)) / 3) ;
	  }
				/* Fetch sequences, desc from ACEDB */
	if (!msp->sseq && lexword2key (msp->sname, &key, _VProtein))
	  { if ((a = peptideGet (key)))
	      { pepDecodeArray (a) ;
		msp->sseq = arrp (a, 0, char) ;
		a->base = 0 ; arrayDestroy (a) ;
	      }
	    else
	      msp->sseq = seq ;	/* empty sseq marker */

	    if ((obj = bsCreate (key)))
	      { if (bsGetKey (obj, _Title, &key))
		  msp->desc = strnew (name (key), 0) ;
		bsDestroy (obj) ;
	      }
				/* Avoid duplication */
	    for (msp2 = msp->next; msp2 ; msp2 = msp2->next)
	      if (!strcmp (msp->sname, msp2->sname)) 
		{ msp2->sseq = msp->sseq ;
		  msp2->desc = msp->desc ;
		}
	  }
      }
  else				/* DNA */
    for (msp = &msp1 ; msp ; msp = msp->next)
      { if (msp->qstart < 1)
	  { msp->sstart += 1 - msp->qstart ;
	    msp->qstart = 1 ;
	  }
	if (msp->qend > max-min+1)
	  { msp->send -= max-min+1 - msp->qend ;
	    msp->qend = max-min+1 ;
	  }
				/* Fetch sequences, desc from ACEDB */
	if (!msp->sseq && lexword2key (msp->sname, &key, _VSequence))
	  { if ((a = dnaGet (key)))
	      { dnaDecodeArray (a) ;
		msp->sseq = arrp (a, 0, char) ;
		a->base = 0 ; arrayDestroy (a) ;
	      }
	    else
	      msp->sseq = seq ;	/* empty sseq marker */

	    if ((obj = bsCreate (key)))
	      { if (bsGetKey (obj, _Title, &key))
		  msp->desc = strnew (name (key), 0) ;
		bsDestroy (obj) ;
	      }
				/* Avoid duplication */
	    for (msp2 = msp->next; msp2 ; msp2 = msp2->next)
	      if (!strcmp (msp->sname, msp2->sname)) 
		{ msp2->sseq = msp->sseq ;
		  msp2->desc = msp->desc ;
		}
	  }
      }
				/* remove empty sequence markers */
  for (msp = &msp1 ; msp ; msp = msp->next)
    if (msp->sseq == seq)
      msp->sseq = 0 ;


  if (msp1.score)
    { msp = (MSP*) messalloc (sizeof (MSP)) ;
      *msp = msp1 ;
    }

  blxview (seq, "", first, offset, msp, opts) ;
}

/**************************************************************/

static void fMapShowPepAlign (int box)
{
  SEG *seg ;
  int x1, x2 ;
  BOOL up = FALSE ;
  FMAPLOOKGET ("fMapShowPepAlign") ;

  if (! (seg = BOXSEG (box)))
    return ;

  if (seg->type == HOMOL_UP) up = TRUE ;

  if (look->flag & FLAG_REVERSE) /* reverse back temporarily */
    { x1 = look->length + 1 - seg->x2 ;
      x2 = look->length + 1 - seg->x1 ;
      fMapRC (look) ;
      look->origin = look->length + 1 - look->origin ;
      doShowAlign (x1, x2, up ? "X+Br" : "X-Br" , 0) ;
      fMapRC (look) ;
      look->origin = look->length + 1 - look->origin ;
    }
  else
    doShowAlign (seg->x1, seg->x2, up ? "X-Br" : "X+Br", 0) ;
}

static void fMapShowPepAlignThis (int box)
{
  SEG *seg ;
  HOMOLINFO *hinf ;
  int x1, x2 ;
  BOOL up = FALSE ;
  FMAPLOOKGET ("fMapShowPepAlignThis") ;

  if (! (seg = BOXSEG (box)))
    return ;
  hinf = arrp (look->homolInfo, seg->data.i, HOMOLINFO) ;
  if (seg->type == HOMOL_UP) up = TRUE ;

  if (look->flag & FLAG_REVERSE) /* reverse back temporarily */
    { x1 = look->length + 1 - seg->x2 ;
      x2 = look->length + 1 - seg->x1 ;
      fMapRC (look) ;
      look->origin = look->length + 1 - look->origin ;
      doShowAlign (x1, x2, up ? "X+Br" : "X-Br", hinf->method) ;
      fMapRC (look) ;
      look->origin = look->length + 1 - look->origin ;
    }
  else
    doShowAlign (seg->x1, seg->x2, up ? "X-Br" : "X+Br", hinf->method) ;
}

static void fMapShowDNAAlign (int box)
{
  SEG *seg ;
  BOOL up = FALSE ;
  FMAPLOOKGET ("fMapShowDNAAlign") ;

  if (! (seg = BOXSEG (box)))
    return ; 
  if (seg->type == HOMOL_UP) up = TRUE ;
  doShowAlign (seg->x1, seg->x2, up ? "N-B" : "N+B", 0) ;
}

static void fMapShowDNAAlignThis (int box)
{
  SEG *seg ;
  HOMOLINFO *hinf ;
  BOOL up = FALSE ;
  FMAPLOOKGET ("fMapShowDNAAlignThis") ;

  if (! (seg = BOXSEG (box)))
    return ; 
  if (seg->type == HOMOL_UP) up = TRUE ;
  hinf = arrp (look->homolInfo, seg->data.i, HOMOLINFO) ;
  doShowAlign (seg->x1, seg->x2, up ? "N-B" : "N+B", hinf->method) ;
}

static void fMapExpandAlign (int box)
{
  SEG *seg ; 
  int i, j , xMin, xMax ;
  KEYSET ks ;
  SegType type ;
  FMAPLOOKGET ("fMapExpandAlign") ;

  if (! (seg = BOXSEG (box)))
    return ; 
  type = seg->type | 0x1 ;

  ks = keySetCreate () ;
  xMin = look->zoneMin ; xMax = look->zoneMax ;
  for (i = 0, j = 0 ; i < arrayMax (look->segs) ; ++i)
    { 
      seg = arrp (look->segs, i, SEG) ;
      if (seg->x1 > xMax || seg->x2 < xMin)
	continue ;
      if ((seg->type | 0x1) != type)
	continue ;
      keySet (ks, j++) = seg->key ;
    }
  
  keySetSort (ks) ;
  keySetCompress (ks) ;
  if (keySetMax (ks))
    forestDisplayKeySet (0, ks, FALSE) ; /* will destroy ks */
  else
    keySetDestroy (ks) ;
}

/***************************************************************/

/* fMaptRNAfold uses nip to fold tRNA.  Currently doesn't combine exons - fix!
*/
static void fMaptRNAfold (int box)
{
  SEG *myseg;
  int i, min, max;
  char *seq, *cp, *cq ;  
  FMAPLOOKGET ("fMapRNAfold") ;

  if (! (myseg = BOXSEG (box)))
    return ;

  min = myseg->x1;
  max = myseg->x2;

  seq = (char*) messalloc (max - min +2);
  for (i = min, cp = seq, cq = arrp (look->dna, min, char) ; i++ <= max ; )
    *cp++ = dnaDecodeChar[ (int)*cq++] ;
  *cp = 0 ;

  externalFileDisplay ("RNA fold", callScriptPipe ("tRNAfold", seq), 0) ;

  messfree (seq) ;
}


static MENUOPT pepHomolMenu[] = {
  { (GraphFunc)fMapShowPepAlign,     "Show multiple protein alignment in Blixem"},
  { (GraphFunc)fMapShowPepAlignThis, "Show multiple protein alignment of just this kind of homologies"},
  { (GraphFunc)fMapExpandAlign,      "Develop a table of these homologies" },  
  { (GraphFunc)fMapFetch,            " (slow) Import this protein via the web" },
  { 0, 0 }
} ;

static MENUOPT dnaHomolMenu[] = {
  { (GraphFunc)fMapShowDNAAlign,     "Show multiple dna alignment"},
  { (GraphFunc)fMapShowDNAAlignThis, "Show multiple dna alignment of just this kind of homologies"},
  { (GraphFunc)fMapExpandAlign,      "Develop a table of these homologies" },  
  { 0, 0 }
} ;

/******************** end of external section ***************/

void fMapRemoveSelfHomol (LOOK look)
{
  int i, max ;
  SEG *seg, *cds, *cdsMin, *cdsMax ;
  char *cp ;

  fMapFindSegBounds (look, CDS, &i, &max) ;
  if (i >= max) 
    return ; /* mieg, oct 2004,  otherwise i crash under arraycheck */
  cdsMin = arrp (look->segs, i, SEG) ;
  cdsMax = arrp (look->segs, max-2, SEG) + 2 ; /* mieg, oct 2004, -2+2 for arrayCheck */
  
  fMapFindSegBounds (look, HOMOL, &i, &max) ;
  for ( ; i < max && i < arrayMax (look->segs) ; ++i) 
    /* arrayMax (look->segs) can change */
    { seg = arrp (look->segs, i, SEG) ;
      for (cp = name (seg->key) ; *cp && *cp != ':' ; ++cp) ;
      if (*cp) 
	++cp ;
      else
	cp = name (seg->key) ;
      for (cds = cdsMin ; cds < cdsMax && cds->x1 < seg->x2 ; ++cds)
	if (cds->x1 <= seg->x1 && cds->x2 >= seg->x2 &&
	    !strcmp (cp, name (cds->parent)))
	  { *seg = array (look->segs,arrayMax (look->segs)-1,SEG) ;
	    --arrayMax (look->segs) ;
	    break ;
	  }
    }

  fMapFindSegBounds (look, CDS_UP, &i, &max) ;
  if (i < max && max > 0)   /* mieg oct 2004 */
    {
      cdsMin = arrp (look->segs, i, SEG) ;
      cdsMax = arrp (look->segs, max-1, SEG) + 1 ; /* -1+1 for arrayCheck */
      
      fMapFindSegBounds (look, HOMOL_UP, &i, &max) ;
      for ( ; i < max && i < arrayMax (look->segs) ; ++i) 
	/* arrayMax (look->segs) can change */
	{ seg = arrp (look->segs, i, SEG) ;
	for (cp = name (seg->key) ; *cp && *cp != ':' ; ++cp) ;
	if (*cp) 
	  ++cp ;
	else
	  cp = name (seg->key) ;
	for (cds = cdsMin ; cds < cdsMax && cds->x1 < seg->x2 ; ++cds)
	  if (cds->x1 <= seg->x1 && cds->x2 >= seg->x2 &&
	      !strcmp (cp, name (cds->parent)))
	    { *seg = array (look->segs,arrayMax (look->segs)-1,SEG) ;
	    --arrayMax (look->segs) ;
	    break ;
	    }
	}
    }
}

/****************************************************/

void fMapShowHomol (LOOK look, float *offset) 
{ 
  int i, box ;
  SEG *seg ;
  int frame = look->min % 3 ;
  HOMOLINFO *hinf ;
  BoxCol *bc ;
  char *bubbleClassName, *bubbleMethod ;
  char *cp ;
 
  bc = bcFromName (look, look->map->activeColName) ;
  if (!bc || !bcTestMag (bc, look->map->mag)) return ;
  bcCheck (bc) ;
  bc->offset = *offset ;

  for (i = 0 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs, i, SEG) ;
      if (seg->x1 > look->max || seg->x2 < look->min)
	continue ;
      if (seg->type == HOMOL)
	{ if (!bc->isDown) continue ; }
      else if (seg->type == HOMOL_UP)
	{ if (bc->isDown && bc->strandSensitive) continue ; }
      else
	continue ;
      if (bit (look->homolBury,seg->data.i)) /* bury */
	continue ;
      hinf = arrp (look->homolInfo, seg->data.i, HOMOLINFO) ;
      if (lexAliasOf (hinf->method) != bc->method)
	continue ;
      if (bc->isFrame && (seg->x1 % 3 != frame)) 
	continue ;
      box = bcDrawBox (bc, seg, hinf->score) ;
      if (box)
	{ array (bc->look->boxIndex, box, int) = i ;
	  if (bc->blixemX)
	    graphBoxMenu (box, pepHomolMenu) ;
	  else if (bc->blixemN)
	    graphBoxMenu (box, dnaHomolMenu) ;
	  if (bc->isShowText)
	    addVisibleInfo (look, i) ;
	  
	  if (hinf->method)
	    bubbleMethod = name (hinf->method) ;
	  else
	    bubbleMethod = "Homol" ;
	  if (class (seg->parent))
	    bubbleClassName = className (seg->parent) ;
	  else
	    bubbleClassName = bubbleMethod ;

	  if (bubbleMethod && !strcasecmp (bubbleMethod, "pfam"))
	    {
	      KEY def = keyGetKey (seg->parent, str2tag("Definition")) ;
	      cp = messprintf ("%s [Pfam]", 
			       def ? name(def) : name(seg->parent)
			       ) ;
	      graphBubbleInfo (box, 1, "mrnaPfam", cp) ;
	    }
	  else if (bubbleMethod && !strcasecmp (bubbleMethod, "psort"))
	    {
	      cp = messprintf ("%s [Psort]",  name(seg->parent)) ;
	      graphBubbleInfo (box, 1, "Psort", cp) ;
	    }
	  else if (bubbleMethod && !strcasecmp (bubbleMethod, "blastp"))
	    {
	      cp = "Protein homologies (BlastP with expect 0.001)" ;
	      graphBubbleInfo (box, 1, "BlastP", cp) ;
	    }

	  else
	    graphBubbleInfo (box, seg->parent, bubbleClassName, 
			     messprintf ("%s homologies", bubbleMethod)) ;


	}
    }

  if (bc->bump) { bumpDestroy (bc->bump) ; bc->bump = 0 ; }
  *offset += bc->fmax ;
}

/****************************************************/

void fMapShowFeature (LOOK look, float *offset) 
{ 
  int i, box ;
  SEG *seg ;
  int frame = look->min % 3 ;
  BoxCol *bc ;

  bc = bcFromName (look, look->map->activeColName) ;
  if (!bc || !bcTestMag (bc, look->map->mag)) return ;
  bc->offset = *offset ;
  for (i = 0 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs, i, SEG) ;
      if (seg->x1 > look->max || seg->x2 < look->min)
	continue ;
      if (! ((bc->isDown && (seg->type == FEATURE ||
			   (seg->type == FEATURE_UP && !bc->strandSensitive)))
	    || 
	    (!bc->isDown && seg->type == FEATURE_UP && bc->strandSensitive)))
	continue ;
      if (lexAliasOf (seg->key) != bc->method)
	continue ;
      if (bc->isFrame && (seg->x1 % 3 != frame)) 
	continue ;
      box = bcDrawBox (bc, seg, seg->data.f) ;
      if (box)
	{ array (bc->look->boxIndex, box, int) = i ;
	  if (bc->isShowText && seg->parent) /* -seg->parent = featDict index */
	    addVisibleInfo (look, i) ;
	}
    }

  if (bc->bump) { bumpDestroy (bc->bump) ; bc->bump = 0 ; }
  *offset += bc->fmax ;
}

/********************************************/

void fMapShowAlleles (LOOK look, float *offset)
{
  char *cp ;
  int i, box ;
  SEG *seg ;
  BOOL isTransposon, isDeletion ;
  BOOL isDraw = FALSE ;
  static char DNAchars[] = "AGCTagct" ;
  static Array xy = 0 ;
  
  if (!xy)
    { xy = arrayCreate (6, float) ;
      arrayMax (xy) = 6 ;
    }
	  
  for (i = 0 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs, i, SEG) ;
      if (seg->type > ALLELE_UP)
	break ;
      if ((seg->type != ALLELE && seg->type != ALLELE_UP) ||
	  (seg->parent != seg->key)) /* no explicit coords given */
	continue ;

      if (look->map->mag * (seg->x2 - seg->x1) < 1) /* point symbol */
	{ float y = MAP2GRAPH (look->map, 0.5* (seg->x1+seg->x2)) ;
	  if (y > mapGraphHeight || y < topMargin)
	    continue ;
	  
	  isDeletion = isTransposon = FALSE ;
	  if ((cp = seg->data.s))
	    {
	      if (*cp == '-')
		isDeletion = TRUE ;
	      else
		{
		  while (*cp)
		    if (!strchr (DNAchars, *cp++))
		      isTransposon = TRUE ;
		}
	    }

	  box = graphBoxStart () ;
	  if (isDeletion)
	    { float y1 = MAP2GRAPH (look->map, seg->x1) ;
	      float y2 = MAP2GRAPH (look->map, seg->x2) ;
	      graphLine (*offset, y1, *offset+1.0, y1) ;
	      graphLine (*offset+1.0, y1, *offset+2.0, y-0.7) ;
	      graphLine (*offset, y2, *offset+1.0, y2) ;
	      graphLine (*offset+1.0, y2, *offset+2.0, y+0.7) ;
	      graphLine (*offset+2.0, y-0.7, *offset+2.0, y+0.7) ;
	    }
	  else if (isTransposon)
	    { arr (xy, 0, float) = *offset ; arr (xy, 1, float) = y ;
	      arr (xy, 2, float) = *offset + 1.0 ; arr (xy, 3, float) = y - 0.7 ;
	      arr (xy, 4, float) = *offset + 1.0 ; arr (xy, 5, float) = y + 0.7 ;
	      graphPolygon (xy) ;
	      graphLine (*offset+1.5, y-0.5, *offset+1.5, y+0.5) ;
	      if (seg->type == ALLELE_UP)
		{ graphLine (*offset+1.1, y, *offset+1.5, y-0.5) ;
		  graphLine (*offset+1.9, y, *offset+1.5, y-0.5) ;
		}
	      else
		{ graphLine (*offset+1.1, y, *offset+1.5, y+0.5) ;
		  graphLine (*offset+1.9, y, *offset+1.5, y+0.5) ;
		}
	    }
	  else
	    { graphLine (*offset, y, *offset+0.5, y+0.4) ;
	      graphLine (*offset, y, *offset+0.5, y-0.4) ;
	      graphLine (*offset, y, *offset+1.1, y+0.2) ;
	      graphLine (*offset+1.1, y+0.2, *offset+0.9, y-0.2) ;
	      graphLine (*offset+2, y, *offset+0.9, y-0.2) ;
	    }
	}
      else			/* bracket to show region */
	{ float y1 = MAP2GRAPH (look->map, seg->x1) ;
	  float y2 = MAP2GRAPH (look->map, seg->x2) ;
	  if (y1 > mapGraphHeight || y2 < topMargin)
	    continue ;

	  box = graphBoxStart () ;
	  if (y1 < topMargin) 
	    graphLine (*offset+2.0, topMargin, *offset+2.0, y2) ;
	  else
	    { graphLine (*offset+1.0, y1, *offset+2.0, y1) ;
	      graphLine (*offset+2.0, y1, *offset+2.0, y2) ;
	    }
	  graphLine (*offset+1.0, y2, *offset+2.0, y2) ;
	}
      graphBoxEnd () ;
      graphBoxDraw (box, BLACK, TRANSPARENT) ;
      isDraw = TRUE ;
      array (look->boxIndex, box, int) = i ;
      fMapBoxInfo (look, box, seg) ;
      array (look->visibleIndices,arrayMax (look->visibleIndices),int) = i ;
    }

  if (isDraw)
    *offset += 2.0 ;
}

/********************************************/

static void showSplices (LOOK look, SegType type, BoxCol *bc, float origin)
{
  char  *v ;
  int i, box, background ;
  SEG *seg ;
  float y, delta=0, x ; /*       delta: mieg: shows frame by altering the drawing of the arrow */
  float fac ;
  
  if (!bc) return ;

  fac = bc->width / (bc->maxScore - bc->minScore) ;
  
  for (i = 0 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs, i, SEG) ;
      if (seg->x1 > look->max || seg->x2 < look->min || 
	  seg->key != bc->method || seg->type != type)
	continue ;
      y = MAP2GRAPH (look->map, seg->x2) ;
      if (seg->data.f <= bc->minScore) 
	x = 0 ;
      else if (seg->data.f >= bc->maxScore) 
	x = bc->width ;
      else 
	x = fac * (seg->data.f - bc->minScore) ;
      box = graphBoxStart () ;
      if (x > origin + 0.5 || x < origin - 0.5) 
	graphLine (bc->offset+origin, y, bc->offset+x, y) ;
      else if (x > origin)
	graphLine (bc->offset+origin-0.5, y, bc->offset+x, y) ;
      else
	graphLine (bc->offset+origin+0.5, y, bc->offset+x, y) ;
      switch (type)
	{
	case SPLICE5:
	  delta = (look->flag & FLAG_REVERSE) ? -0.5 : 0.5 ;
	  break ;
	case SPLICE3:
	  delta = (look->flag & FLAG_REVERSE) ? 0.5 : -0.5 ;
	  break ;
        default:
	  messcrash ("Bad type %d in showSplices", type) ;
	}
      graphLine (bc->offset+x, y, bc->offset+x, y+delta) ;
      graphBoxEnd () ;
      v = SEG_HASH (seg) ;
      if (assFind (look->chosen, v, 0))
	background = GREEN ;
      else if (assFind (look->antiChosen, v, 0))
	background = LIGHTGREEN ;
      else
	background = TRANSPARENT ;
      switch (seg->x2 % 3)
	{
	case 0:
	  graphBoxDraw (box, RED, background) ; break ;
	case 1:
	  graphBoxDraw (box, BLUE, background) ; break ;
	case 2:
	  graphBoxDraw (box, DARKGREEN, background) ; break ;
	}
      array (look->boxIndex, box, int) = i ;
      fMapBoxInfo (look, box, seg) ;
      graphBoxFreeMenu (box, fMapChooseMenuFunc, fMapChooseMenu) ;
    }
}

void fMapShowSplices (LOOK look, float *offset)
{
  float x, y1, y2, origin ;
  BoxCol *bc ;

  bc = bcFromName (look, look->map->activeColName) ;
	/* NB old default diff was 5 - now using bc system it is 1 */
  if (!bc || !bcTestMag (bc, look->map->mag)) return ;
  bc->offset = *offset ;

  y1 = MAP2GRAPH (look->map,look->min) ;
  y2 = MAP2GRAPH (look->map,look->max) ;

  graphColor (LIGHTGRAY) ;
  x = *offset + 0.5 * bc->width ;  graphLine (x, y1, x, y2) ;
  x = *offset + 0.75 * bc->width ;  graphLine (x, y1, x, y2) ;

  graphColor (DARKGRAY) ;
  if (bc->minScore < 0 && 0 < bc->maxScore)
    origin = bc->width * (-bc->minScore / (bc->maxScore - bc->minScore)) ;
  else
    origin = 0 ;
  graphLine (*offset + origin, y1, *offset + origin, y2) ;

  graphColor (BLACK) ;
  showSplices (look, SPLICE5, bc, origin) ;
  showSplices (look, SPLICE3, bc, origin) ;

  *offset += bc->width + 1 ;
}

/*************************************************/

void fMapShowCDSBoxes (LOOK look, float *offset)
{
  int i ;
  SEG *seg ;
  float y1, y2 ;

  graphText ("Genes", *offset-1, topMargin+1.5) ;
  *offset += 2.5 ;
  graphLine (*offset, topMargin+2, *offset, mapGraphHeight-0.5) ;

  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs,i,SEG) ;
      if (seg->x1 > look->max || seg->x2 < look->min)
	continue ;
      y1 = MAP2GRAPH (look->map,seg->x1) ;
      y2 = MAP2GRAPH (look->map,seg->x2+1) ;	/* to cover full base */
      if (y1 < 2 + topMargin) y1 = 2 + topMargin ;
      if (y2 < 2 + topMargin) y2 = 2 + topMargin ;
      if (seg->type == CDS)
	graphFillRectangle (*offset, y1, *offset+1, y2) ;
      if (seg->type == CDS_UP)
	graphFillRectangle (*offset, y1, *offset-1, y2) ;
    }
  *offset += 2.5 ;
}

void fMapShowCDSLines (LOOK look, float *offset)
{
  int i ;
  SEG *seg ;
  float y ;

  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs,i,SEG) ;
      if (seg->x1 > look->max || seg->x2 < look->min)
	continue ;
      y = 0.5 * (MAP2GRAPH (look->map,seg->x1) +
		 MAP2GRAPH (look->map,seg->x2+1)) ;
      if (y < 2 + topMargin) y = 2 + topMargin ;
      if (seg->type == CDS || seg->type == CDS_UP)
	graphLine (*offset+1, y, *offset+3, y) ;
    }
  *offset += 5 ;
}

static void generalBump (LOOK look, float *offset, int width, KEY tag)
{
  int i, box ;
  SEG *seg ;
  float y ;
  BUMP bump = bumpCreate (width, 0) ;

  bumpOff = *offset ;

  graphColor (BLACK) ;
  i = (look->flag & FLAG_REVERSE) ? arrayMax (look->segs)-1 : 1 ;
  while (i && i < arrayMax (look->segs))
    { seg = arrp (look->segs,i,SEG) ;
      if (seg->x1 > look->max || seg->x2 < look->min)
	goto loop ;
      y = 0.5 * (MAP2GRAPH (look->map,seg->x1) +
		 MAP2GRAPH (look->map,seg->x2+1)) ;
      if (y < 2 + topMargin) y = 2 + topMargin ;
      if (seg->type == VISIBLE && 
	  seg->data.k == tag) 
	{ box = graphBoxStart () ;
	  bumpWrite (bump, name (seg->key),0,y) ;
	  graphBoxEnd () ;
	  array (look->boxIndex, box, int) = i ;
	  fMapBoxInfo (look, box, seg) ;
	}
    loop: if (look->flag & FLAG_REVERSE) --i ; else ++i ;
    }

  *offset += bumpMax (bump) ;
  bumpDestroy (bump) ;
}

void fMapShowGeneNames (LOOK look, float *offset)
{ generalBump (look, offset, 6, _Locus) ; }

void fMapShowBriefID (LOOK look, float *offset)
{ generalBump (look, offset, 50, _Brief_identification) ; }

void fMapShowTitles (LOOK look, float *offset)
{ generalBump (look, offset, 50, _Title) ; }

/****** code for fMapShowCanonical ******/

void keyZoneAdd (Array a, KEY key, float y1, float y2, int iseg)
{ 
  int i ;
  KEYZONE *z ;

  for (i = 0 ; i < arrayMax (a) ; ++i)
    if ((z = arrp (a, i, KEYZONE))->key == key)
      return ;

  z = arrayp (a, arrayMax (a), KEYZONE) ;
  z->key = key ; z->y1 = y1 ; z->y2 = y2 ; z->iseg = iseg ;
}

int keyZoneOrder (const void *va, const void *vb)
{
  const KEYZONE *a = (const KEYZONE*) va, *b = (const KEYZONE*) vb ;
  return (a->y1 + a->y2 > b->y1 + b->y2) ? 1 : -1 ;
}

void fMapShowCanonical (LOOK look, float *offset)
{
  int i ;
  SEG *seg, *oldseg = 0 ;
  BOOL isDown ;
  int x, box ;
  float y1, y2, y ;
  BoxCol *bc ;
  BUMP bump ;
  Array a ; 

  if (!fmapView (look->view, str2tag ("Fmap_Clone")))
    return ;

  bc = bcFromName (look, "FMAP_Sequences_and_ends") ;
  bump = bumpCreate (3, 0) ; 
  a = arrayCreate (8, KEYZONE) ;

				/* first the bars */
  ++*offset ;
  i = (look->flag & FLAG_REVERSE) ? arrayMax (look->segs)-1 : 1 ;
  while (i && i < arrayMax (look->segs))
    { seg = arrp (look->segs,i,SEG) ;
      if (seg->x1 > look->max || seg->x2 < look->min)
	goto loop1 ;
      y1 = MAP2GRAPH (look->map,seg->x1) ;
      y2 = MAP2GRAPH (look->map,seg->x2+1) ;	/* to cover full base */
      if (y1 < 2 + topMargin) y1 = 2 + topMargin ;
      if (y1 > mapGraphHeight) y1 = mapGraphHeight ;
      if (y2 < 2 + topMargin) y2 = 2 + topMargin ;
      if (y2 > mapGraphHeight) y2 = mapGraphHeight ;
      if (seg->type == SEQUENCE || seg->type == SEQUENCE_UP)
	{ SEQINFO *sinf = arrp (look->seqInfo, seg->data.i, SEQINFO) ;
	  if (sinf->flags & SEQ_CANONICAL)
	    { x = 0 ;
	      y = seg->x1 ;	/* bump in DNA coords so OK reversed */
	      bumpItem (bump, 1, seg->x2-seg->x1+1, &x, &y) ;
	      box = graphBoxStart () ;
	      graphRectangle (*offset+x, y1, *offset+x+0.75, y2) ;
	      graphBoxEnd () ;
	      graphBoxMenu (box, fMapClonesMenu) ; 

	      array (look->boxIndex, box, int) = i;
	      fMapBoxInfo (look, box, seg) ;
	      if (bc && bc->isShowText)
		addVisibleInfo (look, i) ;
	      keyZoneAdd (a, seg->key, y1, y2, i) ;
	    }
	}
    loop1: if (look->flag & FLAG_REVERSE) --i ; else ++i ;
    }
  *offset += bumpMax (bump) ;
  bumpDestroy (bump) ;

				/* then the ends */
  bump = bumpCreate (3, 0) ;
	/* determine isDown: rely on MASTER running Left->Right */
  seg = arrp (look->segs, 0, SEG) ;
  isDown = (seg->x2 > seg->x1) ;
  if (look->flag & FLAG_REVERSE)
    isDown = !isDown ;
  if (look->flag & FLAG_COMPLEMENT) /* because Clone_*_end tags not flipped by fMapRC () */
    isDown = !isDown ;

  i = (look->flag & FLAG_REVERSE) ? arrayMax (look->segs)-1 : 1 ;
  while (i && i < arrayMax (look->segs))
    { seg = arrp (look->segs,i,SEG) ;
      if (seg->type != CLONE_END ||
	  seg->x1 > look->max || seg->x2 < look->min)
	goto loop2 ;
      if (oldseg && seg->key == oldseg->key && seg->x1 == oldseg->x1)
	goto loop2 ;
      oldseg = seg ;
      y = MAP2GRAPH (look->map,seg->x1) ;
      box = graphBoxStart () ;
      array (look->boxIndex, box, int) = i;
      graphLine (*offset, y, *offset+1, y) ;
      if ((seg->data.k == _Clone_left_end && isDown) ||
	  (seg->data.k == _Clone_right_end && !isDown))
	{ graphLine (*offset+0.5, y, *offset+0.5, y+2) ;
	  graphLine (*offset, y+1, *offset+0.5, y+2) ;
	  graphLine (*offset+1, y+1, *offset+0.5, y+2) ;
	}
      else
	{ graphLine (*offset+0.5, y, *offset+0.5, y-2) ;
	  graphLine (*offset, y-1, *offset+0.5, y-2) ;
	  graphLine (*offset+1, y-1, *offset+0.5, y-2) ;
	}
      graphBoxEnd () ;

      fMapBoxInfo (look, box, seg) ;
      keyZoneAdd (a, seg->parent, y, y, i) ;
    loop2: if (look->flag & FLAG_REVERSE) --i ; else ++i ;
    }
  *offset += bumpMax (bump) + 1 ;
  bumpDestroy (bump) ;

				/* then the text next to the clones */
  bump = bumpCreate (12, 0) ;
  bumpOff = *offset ;
  arraySort (a, keyZoneOrder) ;
  for (i = 0 ; i < arrayMax (a) ; ++i)
    { KEYZONE *z = arrp (a, i, KEYZONE) ;
      box = graphBoxStart () ;
      array (look->boxIndex,box,int) = z->iseg ;
      bumpWrite (bump, name (z->key), 0, 0.5* (z->y1 + z->y2)) ;
      graphBoxEnd () ;
      fMapBoxInfo (look, box, seg) ;
    }
  *offset += bumpMax (bump) ;
  bumpDestroy (bump) ;

  arrayDestroy (a) ;
}

/***********************************************************************/

void fMapShowAssemblyTags (LOOK look, float *offset)
{
  int i;
  SEG *seg ;
  float y1, y2 ;
  BUMP bump = bumpCreate (30, 0) ;

				/* first the bars */
  ++*offset ;
  i = (look->flag & FLAG_REVERSE) ? arrayMax (look->segs)-1 : 1 ;
  while (i && i < arrayMax (look->segs))
    { seg = arrp (look->segs,i,SEG) ;
      if (seg->x1 > look->max || seg->x2 < look->min)
	goto loop ;
      y1 = MAP2GRAPH (look->map,seg->x1) ;
      y2 = MAP2GRAPH (look->map,seg->x2+1) ;	/* to cover full base */
      if (y1 < 2 + topMargin) y1 = 2 + topMargin ;
      if (y2 < 2 + topMargin) y2 = 2 + topMargin ;
      if (seg->type == ASSEMBLY_TAG)
	{ int box = graphBoxStart ();
	  int x = 0 ;
	  int color;
	  unsigned char *c = (unsigned char *)seg->data.s;
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
	  y2 -= y1 ;
	  bumpItem (bump, 1, y2, &x, &y1) ;
	  graphRectangle (*offset+x, y1, *offset+x+0.8, y1+y2) ;
	  graphBoxEnd ();
	  fMapBoxInfo (look, box, seg) ;
	  graphBoxDraw (box, BLACK, color);
	  array (look->boxIndex, box, int) = i;
	}
    loop: if (look->flag & FLAG_REVERSE) --i ; else ++i ;
    }
  *offset += bumpMax (bump) ;
  bumpDestroy (bump) ;
}

/*******************************************/

void fMapShowcDNAs (LOOK look, float *offset)
{
  int i ;
  SEG *seg ;
  int x ;
  float y ;
  BUMP bump = bumpCreate (3, 0) ;

  i = (look->flag & FLAG_REVERSE) ? arrayMax (look->segs)-1 : 1 ;
  while (i && i < arrayMax (look->segs))
    { seg = arrp (look->segs,i,SEG) ;
      if (seg->x1 > look->max || seg->x2 < look->min)
	goto loop ;
      y = 0.5 * (MAP2GRAPH (look->map,seg->x1) +
		 MAP2GRAPH (look->map,seg->x2+1)) ;
      if (y < 2 + topMargin) y = 2 + topMargin ;
      if (seg->type == VISIBLE &&
	  seg->data.k == _Matching_cDNA)
	{ x = 0 ;
	  bumpItem (bump, 1, 0.5, &x, &y) ;
	  graphFillRectangle (*offset+x, y-0.25, *offset+x+0.8, y+0.25) ;
	}
    loop: if (look->flag & FLAG_REVERSE) --i ; else ++i ;
    }

  *offset += bumpMax (bump) ;
  bumpDestroy (bump) ;
}

/**********************************************/

static FREEOPT segOpts[] = {
  { 8, "User controllable segments" },
  { SPLICE3, "SPLICE3" },
  { SPLICE5, "SPLICE5" },
  { FEATURE, "FEATURE" },
  { HOMOL, "HOMOL" },
  { SPLICE3_UP, "SPLICE3_UP" },
  { SPLICE5_UP, "SPLICE5_UP" },
  { FEATURE_UP, "FEATURE_UP" },
  { HOMOL_UP, "HOMOL_UP" }
} ;

void fMapAddSegments (void)
{
  void *v ;
  int i, x1, x2, n, level, seqIndex ;
  float score ;
  SEG *seg, *seqSeg ;
  KEY seq, type ;
  static char dirName[DIR_BUFFER_SIZE], fileName[FIL_BUFFER_SIZE] ;
  static Associator seq2seg = 0 ;
  FILE *fil ;
  char *word ;
  FMAPLOOKGET ("fMapAddSegments") ;

  if (! (fil = filqueryopen (dirName, fileName, "useg","r",
	   "File of new items: sequence x1 x2 type [score|text...]")))
    return ;

	/* first delete temp segs */
  arrayMax (look->segs) = look->lastTrueSeg ;
	/* make seq2seg */
  seq2seg = assReCreate (seq2seg) ;
  for (i = 1, seg = arrp (look->segs, 1, SEG) ; i < arrayMax (look->segs) ; ++i, ++seg)
    if (seg->type == SEQUENCE || seg->type == SEQUENCE_UP)
      assInsert (seq2seg, assVoid (seg->key), assVoid (i)) ;

  n = 0 ;
  level = freesetfile (fil,0) ;
  while (freecard (level))
    { if (messIsInterruptCalled ())
	break ;
      if (! (word = freeword ()))
	continue ;
      if (!lexword2key (word, &seq, _VSequence))
	{ messout ("Can't find sequence %s line %d - F4 to interrupt",
		   word, freestreamline (level)) ;
	  continue ;
	}
      if (!assFind (seq2seg, assVoid (seq), &v))
	continue ;
      seqIndex = assInt (v) ;
      if (!freeint (&x1) || !freeint (&x2) || !freekey (&type, segOpts))
	{ messout ("No coords then segtype line %d - F4 to interrupt", 
		   freestreamline (level)) ;
	  continue ;
	}
      if (x1 > x2)
	{ i = x1 ; x1 = x2 ; x2 = i ; 
	  if (type & 0x01) --type ; else ++type ;
	}
      seg = arrayp (look->segs, arrayMax (look->segs), SEG) ;
      seg->parent = 0 ;
      seqSeg = arrp (look->segs, seqIndex, SEG) ;
      if (seqSeg->type == SEQUENCE)
	{ seg->x1 = seqSeg->x1 + x1 - 1 ;
	  seg->x2 = seqSeg->x1 + x2 - 1 ;
	}
      else			/* SEQUENCE_UP */
	{ seg->x1 = seqSeg->x2 - x2 + 1 ;
	  seg->x2 = seqSeg->x2 - x1 + 1 ;
	  if (type & 0x01) --type ; else ++type ;
	}
      seg->type = type ;
      ++n ;

#define ABORT(z) { --arrayMax (look->segs) ; \
	      messout ("%s line %d - F4 to interrupt", z, freestreamline (level)) ; \
	      continue ; }

      switch (type)
	{
	case SPLICE3: case SPLICE3_UP: 
	case SPLICE5: case SPLICE5_UP:
	case FEATURE: case FEATURE_UP:
	  if (! (word = freeword ()))
	    ABORT ("No method") ;
	  if (lexaddkey (word, &seg->key, _VMethod))
	    messout ("Previously unknown method %s line %d", 
		     word, freestreamline (level)) ;
	  methodSet (word, WHITE, 0, 0, 0.5, 0, 0, 0) ;
		/* default - will read from database entry if exists */
	  seg->data.f = freefloat (&score) ? score : 0.0 ;
	  break ;
	case HOMOL: case HOMOL_UP:
	  { HOMOLINFO *hinf ;
	    seg->data.i = arrayMax (look->homolInfo) ;
	    hinf = arrayp (look->homolInfo, seg->data.i, HOMOLINFO) ;
	    if (! (word = freeword ()))
	      ABORT ("No method for HOMOL") ;
	    if (lexaddkey (word, &hinf->method, _VMethod))
	      messout ("Previously unknown method %s for HOMOL, line %d", 
		       word, freestreamline (level)) ;
	    methodSet (word, WHITE, 0, 0, 0.5, 0, 0, 0) ;
	    if (!freefloat (&hinf->score))
	      ABORT ("No score for HOMOL") ;
	    if (! (word = freeword ()))
	      ABORT ("No sequence for HOMOL") ;
	    lexaddkey (word, &seg->key, _VSequence) ;
	    if (!freeint (&hinf->x1) || !freeint (&hinf->x2))
	      ABORT ("No target coords for HOMOL") ;
	  }
	  break ;
	}
    }
  freeclose (level) ;
  messout ("Added %d new segments", n) ;

  arraySort (look->segs, fMapOrder) ;
  look->lastTrueSeg = arrayMax (look->segs) ;

  fMapProcessMethods (look) ;
  fMapDraw (look, 0) ;
}

void fMapClearSegments (void)
{
  int i ;
  KEY method ;
  SEG *seg ;
  FMAPLOOKGET ("fMapClearSegments") ;

  if (!messPrompt ("Give method of segment to clear", "", "wz") ||
      !lexword2key (freeword (), &method, _VMethod))
    return ;

	/* first delete temp segs */
  arrayMax (look->segs) = look->lastTrueSeg ;

  for (i = 1 ; i < arrayMax (look->segs) ; )
    if ((seg = arrp (look->segs, i, SEG))->key == method)
      *seg = arr (look->segs, --arrayMax (look->segs), SEG) ;
    else
      ++i ;

  arraySort (look->segs, fMapOrder) ;
  look->lastTrueSeg = arrayMax (look->segs) ;

  fMapDraw (look, 0) ;
}

/************************************************************/
/******* generic "box" column - for homols, features ********/

static void bcCheck (BoxCol *bc)	/* checks before using a column */
{
 if (bc)
   {
     bc->isFrame = bc->look->map->isFrame && bc->frameSensitive ;
     
     bc->mode = WIDTH ;
     if (bc->byOffset && !bc->look->map->isFrame)
       bc->mode = OFFSET ;
     else if (bc->byHist)
       bc->mode = HIST ;
     
     if (bc->autoWidth)
       {
	 if (bc->mode == OFFSET)
	   bc->width = 7 ;
	 else
	   bc->width = 2 ;
       }
     else if (bc->mode == WIDTH && bc->byOffset && bc->width > 2)
       bc->width = 2 ;
     
     if (bc->mode == HIST)
       { if (bc->histBase < 0) bc->histBase = 0 ;
       if (bc->histBase > 1) bc->histBase = 1 ;
       bc->fmax = bc->width*bc->histBase + 0.2 ;
       }
     else
       bc->fmax = 0 ;
     /* checks to prevent arithmetic crash */
     if (bc->mode == OFFSET && !bc->minScore)
       bc->minScore = 1 ;
     if (bc->mode != OFFSET && bc->maxScore == bc->minScore)
       bc->maxScore = bc->minScore + 1 ;
     
     if (bc->isBump)
       { if (bc->bump) bumpDestroy (bc->bump) ;
       bc->bump = bumpCreate (1000, 0) ;
       }
     else if (bc->isCluster)
       { 
	 if (bc->cluster) assDestroy (bc->cluster) ;
	 bc->cluster = assHandleCreate (bc->look->handle) ;
	 bc->clusterCount = 0 ;
       }
   }
}

static BoxCol *bcFromName (LOOK look, char *name)
{
  int i ;
  BoxCol *bc ;

  if (!look->bcDict)
    { look->bcDict = dictHandleCreate (64, look->handle) ;
      look->bcArray = arrayHandleCreate (64, BoxCol*, look->handle) ;
    }

  if (dictFind (look->bcDict, name, &i))
    bc = arr (look->bcArray,i,BoxCol*) ;
  else
    {				/* initialise from method */
      METHOD *meth ;
      OBJ obj ;
      
      bc = (BoxCol*) halloc (sizeof (BoxCol), look->handle) ;
      dictAdd (look->bcDict, name, &i) ;
      array (look->bcArray,i,BoxCol*) = bc ;
      
      bc->look = look ;

      if (*name == '-')
	{ bc->isDown = FALSE ;
	  lexaddkey (name+1, &bc->method, _VMethod) ;
	}
      else
	{ bc->isDown = TRUE ;
	  lexaddkey (name, &bc->method, _VMethod) ;
	}
      meth = method (0, bc->method) ;
      
      bc->minScore = meth->min ;
      bc->maxScore = meth->max ;
      bc->colour = meth->col ;
      bc->byWidth = meth->flags & METHOD_SCORE_BY_WIDTH ;
      bc->byOffset = meth->flags & METHOD_SCORE_BY_OFFSET ;
      bc->byHist = meth->flags & METHOD_SCORE_BY_HIST ;
      bc->frameSensitive = meth->flags & METHOD_FRAME_SENSITIVE ;
      bc->strandSensitive = meth->flags & METHOD_STRAND_SENSITIVE ;
      bc->isBump = meth->flags & METHOD_BUMPABLE ;
      bc->blixemX = meth->flags & METHOD_BLIXEM_X ;
      bc->blixemN = meth->flags & METHOD_BLIXEM_N ;
	
      if (meth->width)
	bc->width = meth->width ;
      else
	bc->autoWidth = TRUE ;
	
      if (bc->byHist)
	{ if (bc->minScore == bc->maxScore) bc->maxScore = bc->minScore + 1 ;
	  bc->histBase = (meth->histBase - bc->minScore) / 
	    (bc->maxScore - bc->minScore) ;
	}

      bc->minMag = bc->maxMag = 0 ;
      bc->isCluster = FALSE ;
      bc->isShowText = FALSE ;
      if ((obj = bsCreate (bc->method)))
	{ bsGetData (obj, str2tag ("Min_mag"), _Float, &bc->minMag) ;
	  bsGetData (obj, str2tag ("Max_mag"), _Float, &bc->maxMag) ;
	  if (bsFindTag (obj, str2tag ("Cluster")))
	    bc->isCluster = TRUE ;
	  if (bsFindTag (obj, str2tag ("Show_text")))
	    bc->isShowText = TRUE ;
	  bsDestroy (obj) ;
	}
    }

  if (bc) bcCheck (bc) ;
  return bc ;
}

static BOOL bcTestMag (BoxCol *bc, float mag)
{
  if (mag < 0) mag = -mag ;
  mag = 1/mag ;
  if (bc)
    {
      if (bc->minMag && mag < bc->minMag) return FALSE ;
      if (bc->maxMag && mag > bc->maxMag) return FALSE ;
    }
  return TRUE ;
}

static int bcDrawBox (BoxCol *bc, SEG *seg, float score)
{
  int box, xoff ;
  float y1, y2, dx=0 ;
  double logdeux = log ((double)2.0) ;

  if (!bc) return 0 ;

  if (seg->x1 < bc->look->min)
    y1 = MAP2GRAPH (bc->look->map, bc->look->min) ;
  else
    y1 = MAP2GRAPH (bc->look->map, seg->x1) ;
  if (seg->x2+1 > bc->look->max)
    y2 = MAP2GRAPH (bc->look->map, bc->look->max) ;
  else
    y2 = MAP2GRAPH (bc->look->map, seg->x2+1) ;
  if (y1 >= mapGraphHeight || y2 <= topMargin)
    return 0 ;

  box = graphBoxStart () ;
  switch (bc->mode)
    {
    case OFFSET:
      dx = 1 + ((float)score - bc->minScore) / ((float) (bc->maxScore - bc->minScore)) ;
      if (dx < .8) dx = .8 ; 
      if (dx > 3) dx = 3 ; /* allow some leeway */
      dx = bc->width * log ((double)dx)/logdeux ;
      graphRectangle (bc->offset + dx, y1, 
		      bc->offset + dx + .9, y2) ;
      if (bc->fmax < dx + 1.5)
	bc->fmax = dx + 1.5 ;
      break ;
    case HIST:
      dx = (score - bc->minScore) / (bc->maxScore - bc->minScore) ;
      if (dx < 0) dx = 0 ;
      if (dx > 1) dx = 1 ;
      graphRectangle (bc->offset + bc->width*bc->histBase, y1, 
		      bc->offset + bc->width*dx, y2) ;
      if (bc->fmax < bc->width*dx)
	bc->fmax = bc->width*dx ;
      break ;
    case WIDTH:
      dx = 0.25 + 0.75 * (score - bc->minScore) / 
	 (bc->maxScore - bc->minScore) ;
      if (dx < 0.25) dx = 0.25 ;
      if (dx > 1) dx = 1 ;	/* NB dropthrough */
    case DEFAULT:
      if (bc->mode == DEFAULT) dx = 0.75 ;
      xoff = 1 ;
      if (bc->bump)
	{
	  if (bc->look->map->mag < 0)
	    { y1 = -y1 ; y2 = -y2 ;
	      bumpItem (bc->bump, 2, (y2-y1), &xoff, &y1) ;
	      y1 = -y1 ; y2 = -y2 ;
	    }
	  else
	    bumpItem (bc->bump, 2, (y2-y1), &xoff, &y1) ;
	}
      else if (bc->cluster)	/* one subcolumn per key */
	{ void *v ;
	  if (assFind (bc->cluster, assVoid (seg->key), &v))
	    xoff = assInt (v) ;
	  else
	    { xoff = ++bc->clusterCount ;
	      assInsert (bc->cluster, assVoid (seg->key), assVoid (xoff)) ;
	    }
	}
      graphRectangle (bc->offset + 0.5*bc->width* (xoff - dx), y1, 
		      bc->offset + 0.5*bc->width* (xoff + dx), y2) ;
      if (bc->fmax < 0.5*bc->width* (xoff+1))
	bc->fmax = 0.5*bc->width* (xoff+1) ;
    }
  
  graphBoxEnd () ;
  graphBoxDraw (box, BLACK, bc->colour) ;
  fMapBoxInfo (bc->look, box, seg) ;

  return box ;
}

#ifdef BC_FUNCTIONS_NOT_USED

static BoxCol *bcFromBox (int box)
{
  BoxCol *bc = 0 ;

  return bc ;
}

static void bcSetOverlap (int box)
{
  BoxCol *bc = bcFromBox (box) ;
  FMAPLOOKGET ("bcSetOverlap") ;

  if (bc) bc->isBump = bc->isCluster = FALSE ;
  fMapDraw (look, 0) ;
}

static void bcSetBump (int box)
{
  BoxCol *bc = bcFromBox (box) ;
  FMAPLOOKGET ("bcSetBump") ;

  if (bc) bc->isBump = TRUE ; bc->isCluster = FALSE ;
  fMapDraw (look, 0) ;
}

static void bcSetCluster (int box)
{
  BoxCol *bc = bcFromBox (box) ;
  FMAPLOOKGET ("bcSetCluster") ;

  if (bc) bc->isBump = FALSE ; bc->isCluster = TRUE ;
  fMapDraw (look, 0) ;
}

static void bcSetWidth (int box)
{
  BoxCol *bc = bcFromBox (box) ;
  FMAPLOOKGET ("bcSetWidth") ;

  if (bc) { bc->byWidth = TRUE ; bc->byOffset = FALSE ; bc->byHist = FALSE ; }
  fMapDraw (look, 0) ;
}

static void bcSetOffset (int box)
{
  BoxCol *bc = bcFromBox (box) ;
  FMAPLOOKGET ("bcSetOffset") ;

  if (bc) { bc->byWidth = FALSE ; bc->byOffset = TRUE ; bc->byHist = FALSE ; }
  fMapDraw (look, 0) ;
}

static void bcSetHist (int box)
{ 
  BoxCol *bc = bcFromBox (box) ;
  FMAPLOOKGET ("bcSetHist") ;

  if (bc) { bc->byWidth = FALSE ; bc->byOffset = FALSE ; bc->byHist = TRUE ; }
  fMapDraw (look, 0) ;
}

static MENUOPT bcMenu[] = {
  { menuSpacer, "Handle overlapping entries by:"},
  { (GraphFunc)bcSetOverlap, "  Overlap" }, 
  { (GraphFunc)bcSetBump, "  Bump" },
  { (GraphFunc)bcSetCluster, "  Cluster" },
  { menuSpacer, "Show score by:" },
  { (GraphFunc)bcSetWidth, "  Width" },
  { (GraphFunc)bcSetOffset, "  Offset" },
  { (GraphFunc)bcSetHist, "  Histogram" },
  { 0, 0 }
} ;

#endif /* BC_FUNCTIONS_NOT_USED */

/************************************************************/
/********************** end of file **************************/
