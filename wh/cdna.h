/*  File: cdna.h
 *  Author: mieg
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
 * Description: 
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 15 15:54 1999 (fw)
 * Created: Thu Aug 26 17:55:00 1999 (fw)
 * CVS info:   $Id: cdna.h,v 1.41 2017/08/01 01:17:21 mieg Exp $
 *-------------------------------------------------------------------
 */

#ifndef ACEDB_CDNA_H
#define ACEDB_CDNA_H

#include "acedb.h"
#include "../whooks/sysclass.h"


#define  gI 0x1    /* intron */
#define  gX 0x2    /* exon */
#define  gB 0x4    /* Bis ==alternatif */
#define  gD 0x8    /* debut, sens du gene */
#define  gF 0x10   /* fin, sens du gene */
#define  gS 0x20   /* SL */
#define  gS0 0x40   /* SL0 */
#define  gA 0x80   /* polyA */
#define  gGhost 0x100  /* ghost */
#define  gGap 0x200  /* gap */
#define  gGene 0x400  /* gene */
#define  gDroppedGene 0x800  /* dropped gene */
#define  g3 0x1000
#define  g5 0x2000
#define  gJ 0x4000   /* intron with one error */
#define  gLink 0x8000 /* 3p 5p overlapping read pair, do not reextend */
#define  gCompleteCDS 0x10000 /* complete cds, no jump downstream but */
#define  gOpenGap  0x20000 /* open ORF gap */
#define  gStolen 0x40000   /* exon stolen from another mrna */
#define  gPredicted 0x80000 /* exon stolen from a prediction */
#define  gSuspect 0x100000   /* for suspected problem in exon/intron */
#define  gFuseToGhost 0x200000 /* for fuse_to ghosts */
#define  gReal3p 0x400000
#define  gReal5p 0x800000
#define  gDF 0x1000000 /* debut fuzzy, sens du gene */
#define  gFF 0x2000000 /* fin fuzzy, sens du gene */
#define  gMicro 0x4000000 /* micro intron < minIntronSize */

#define  gDebut (gD | gDF | gS | gS0 | g5 | gReal5p)  /* debut goodies */
#define  gFin (gF | gFF | gA | g3 | gReal3p)  /* fin goodies */
#define  gExact (gD | gF)  /* exact goodies */


extern Array fMapcDNAReferenceDna ;
extern Array fMapcDNAReferenceHits ;
extern KEY fMapcDNAReferenceGene ;
extern KEY fMapcDNAReferenceEst ;
extern int fMapcDNAReferenceOrigin ;
extern int mcAddKeysetAlterSpliceDetails (KEYSET ks, KEY key, BOOL isCoding) ;

void cDNAEliminateDeadMrnas (void) ;
Array cDNAGetReferenceHits (KEY gene, int origin) ;
KEY cDNARealignGene (KEY gene, int z1, int z2, int direction,
		     BOOL doFuse, int doLimit, int searchRepeats, int splitCloud) ;
void fMapcDNADoSelect(KEY k, KEYSET ks, KEYSET taceKs) ;
KEY cDNAMakeDna (KEY gene) ;
void cdnaCleanUp (KEY cosmid, KEY gene, KEYSET genes) ;
void mrnaCleanUp (KEY cosmid, KEY gene, KEYSET genes) ;
int mrnaAddKeysetKantorInfo (KEYSET ks) ;
int mrnaAddKeysetTiling (KEYSET ks)  ;
int cDnaFlipGeneKeySet (KEYSET ks)  ;
int cDnaFlipRenameGenes (KEYSET genes) ;
int cDnaFlagSuspectIntrons (KEYSET ks, BOOL doIgnore) ;
int cDnaExportGeneboxPrimers (KEYSET ks) ;
double intMap2Gmap (int chrom, double X1plusX2half) ;
int cDnaRenameGenes (KEYSET genes, char *fName, char *chromName, 
		     BOOL selectNewName, BOOL renameGene) ;
typedef struct hitStruct { KEY gene, cDNA_clone, est ; BOOL reverse ; 
  int a1, a2, x1, x2, clipTop, clipEnd ; unsigned int type ; int zone, nerr, nerrAll, ea1, ea2, ex1, ex2, maxExact ; } HIT ;

extern  char B2[256],  B2r[256]  ;
extern  KEY _Alternative_exon
, _Alternative_first_exon
, _Alternative_intron
, _Alternative_partial_exon
, _Anomalous_clone
, _Assembled_from
, _Assembled_from_cDNA_clone
, _At_position_1
, _Bad_quality
, _Begin_not_found
, _Best_product
, _Blastp
, _Blastp_title
, _COOH_complete
, _CTF_File
, _Clone_Group
, _Coding_gap
, _Coding_length
, _Complete
, _Composite
, _Confirmed_exon
, _Confirmed_intron
, _Constructed_from
, _Contamination
, _Derived_sequence
, _Discarded_cDNA
, _Discarded_from
, _Duplicate_clone
, _EST_translation
, _Exon
, _Expasy
, _Fake_internal_poly_A
, _First_ATG
, _First_NTG
, _First_Kozak
, _First_exon
, _Fmap_cDNA_Decorate
, _Forward
, _Frame
, _From_AM
, _From_AM_From_prediction
, _From_EST
, _From_gene
, _From_prediction
, _Fuse_to_clone
, _Fuzzy
, _Fuzzy_gt_ag
, _Fuzzy_gc_ag
, _Gap
, _Gene_wall
, _Genes
, _Genomic
, _Genomic_sequence
, _Good_product
, _Hit
, _Hits
, _Ignore_this_clone
, _In_mRNA
, _IntMap
, _Internal_capping
, _Internal_priming
, _Internal_priming_manual
, _Internal_priming_on_A_rich
, _Intron
, _Intron_boundaries
, _Is_AM
, _Is_chain
, _Is_gap
, _Is_read
, _Kantor
, _Last_exon
, _Length_3prime_UTR
, _Length_5prime_UTR
, _Length_anomaly
, _LocusLink
, _Longest_cDNA_clone
, _Longest_CDS
, _Manual_internal_deletion
, _Manual_no_internal_deletion
, _Manual_polyA
, _Mosaic
, _NH2_complete
, _Nb_alternative_exons
, _Nb_confirmed_alternative_introns
, _Nb_confirmed_introns
, _Nb_possible_exons
, _Nb_predicted_exons
, _Nb_stolen_exons
, _NewName
, _No_obvious_phenotype
, _ORF_Gap
, _Open_length
, _Other
, _Overlap_left
, _Overlap_right
, _PCR_product_size
, _Partial_exon
, _Pfam
, _PolyA_after_base
, _Predicted_Exon
, _Predicted_exon
, _Primed_on_polyA
, _Probe_hit
, _Probe_exact_hit
, _Problem
, _Product
, _Protein
, _Psort
, _RNA_editing
, _RT_PCR
, _Read
, _Real_3prime
, _Real_5prime
, _Ref_Seq
, _Ref_mRNA
, _RefSeqMaker
, _Relabelled
, _Representative_product
, _Resequence
, _Reverse
, _Reversed_by
, _SMAP
, _Sage
, _Spliced_sequence
, _Splicing
, _Stolen_Exon
, _Stolen_exon
, _Stolen_intron
, _Suspected_internal_deletion
, _Target
, _Taxblast
, _Tiling_error
, _Title
, _Total_gap_length
, _Total_intron_length
, _Total_length
, _Transcribed_from
, _Transcribed_gene
, _Transcript
, _Transpliced_to
, _UTR_3prime
, _UTR_5prime
, _Use_AM
, _VAnnotation
, _VClone_Group
, _VFiche
, _VIntron
, _VKantor
, _VNewName
, _VRNAi
, _VSage
, _VTranscribed_gene
, _VTranscript
, _VcDNA_clone
, _Vector_clipping
, _Very_good_product
, _VmProduct
, _VmRNA
, _aForward
, _aReverse
, _cDNA_clone
, _ct_ac
, _gc_ag
, _gt_ag
, _mForward
, _mRNA
, _mReverse
, _r_gap
;

void cDNAAlignInit (void) ;
void cDNARemoveCloneFromGene (KEY clone, KEY gene) ; /* cdnaalign.c */
void showHits (Array hits) ;
void cDNASwapA (Array hits) ;
int cDNAOrderByA1 (const void *va, const void *vb) ;
int cDNAOrderGloballyByA1 (const void *va, const void *vb) ;
int cDNAOrderGloballyByA1Errors (const void *va, const void *vb) ;
void cDNASwapX (Array hits) ;
int cDNAOrderByX1 (const void *va, const void *vb) ;
int getPleaseDoRepeats (int set) ;
int cDnaFlagOneSuspectIntrons (KEY tg, BOOL doIgnore, KEYSET result) ;
int giwAddKeysetIntronClass (KEYSET ks, KEY key) ;
int giwAddKeysetProbeWalls (KEYSET ks, KEY key) ;
int giwAddIntronHierarchy (void) ;/* look for cassette introns */

#endif
