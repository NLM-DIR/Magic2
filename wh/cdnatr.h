/*  File: cdnatr.h
 *  Author: Danielle et Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg. 2001
 *-------------------------------------------------------------------
 * This file is part of the ACEMBLY extension to 
        the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
   Align cDNA
 * Exported functions:
 * HISTORY:
 * Created: Mars 2001 (mieg)
 *-------------------------------------------------------------------
 */

#ifndef CDNATR_H_DEF
#define CDNATR_H_DEF

typedef struct s2mStruct {
  BOOL composite ;
  KEY cosmid, cosmid1 ;
  Array plainHits, geneHits, compositeDesignCovering ;
  Array chain, shadows, contigs, gmrnas, linkPos, dnaD, dnaR ;  
  int bMax, type ; void **magic ; 
  AC_HANDLE h ; 
} S2M ;

typedef struct shadowStruct { 
  KEYSET clones, ignoredClones ;
  BitSet yes, no, ghost, suspect ; 
  int bMax ; 
  BOOL isUp ; 
  int begin3p, end3p ; 
  int a1, a2 ;
} SH ;

typedef struct contigStruct { 
  S2M *s2m ; 
  SH sh ; 
  KEYSET ks ; 
  Array s2mContigs ; 
  BOOL isUp ; 
  KEY gene, cGroup ;
  int a1, a2, d1, d2; 
} SC ;

typedef struct smrnaStruct { 
  KEY chrom, cosmid, gene, cGroup, tr ; 
  int g1, g2, a1, a2, bestDna, nGaps ; 
  int stolen, givenBack  ;
  Array genes, tgs, mrnas, introns, clones, estHits, hits, dnas, mdna, orfs, pA ;
  KEYSET x2a ;
  BOOL composite ;
  BOOL isNH2Complete, isCOOHComplete ;
} SMRNA ;

/**********************************************************************************/

typedef struct orfsStruct {
  SMRNA *smrna ;
  Array dna, hits, introns, pep ;

  int frame, p1, p2, m1, m2 ;
  int openLn, codingLn, nOpen ;
  int upStop, downStop, start, cds ;
  int firstATG, quality, met, leu ;
  int nIntron, nIntronsInsideCDS, nIntronsOutsideCDS ;
  int iDna, topExon, nX, upGap, max ;
  int stolen, givenBack ;
  float weight, pI ;
  BOOL isLeu, hasStart, hasStop, isBest, isSL, isUorf, iGood ;
  char ntgType[8] ;
} ORFS ;

/**********************************************************************************/

#define MINI_LEU2MET 180  /* bp, was 90 until 2006_05_12 */ 
#define MINI_LEU2END 240  /* bp, minimal leu only ORF 2006_10_07 */


KEYSET  mrnaCreate (int type, KEY cosmid, KEY cosmid1, 
		    Array dna, Array dnaR, Array plainHits, Array geneHits, 
		    Array linkPos, KEYSET clipTops, KEYSET clipEnds, 
		    KEYSET clonesWithBadIntron) ;

void mrnaDesignSavePepAndKantor (KEY product, Array pep) ;

extern BOOL isComposite ;
BOOL cdnaReExtendHits (S2M *s2m, SC* sc, SMRNA *gmrna, KEYSET clipTops, KEYSET clipEnds) ;
BOOL mrnaTiling (SC *sc, KEY mrna, BOOL useErrors) ;
void mrnaSaveMrna (S2M *s2m, SC* sc, Array estHits, Array smrnas, SMRNA *smrna, BOOL composite) ;
void mrnaSaveMrnaCoding (SMRNA *smrna) ;
int mrnaPleaseNameBy (void) ;

#endif
