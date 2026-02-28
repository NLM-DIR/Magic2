/*  File: cdnainit.h
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
 * Last edited: Aug 26 17:56 1999 (fw)
 * Created: Thu Aug 26 17:55:41 1999 (fw)
 * CVS info:   $Id: cdnainit.c,v 1.11 2009/05/07 02:22:29 mieg Exp $
 *-------------------------------------------------------------------
 */

#include "acedb.h"
#include "lex.h"
#include "dna.h"
#include "cdna.h"

char B2[256] ;
char B2r[256] ;
KEY _Alternative_exon
, _Alternative_first_exon
, _Alternative_intron
, _Alternative_partial_exon
, _Anomalous_clone
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
, _Total_gap_length
, _Total_intron_length
, _Total_length
, _Transcribed_from
, _Transcribed_gene
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
, _mReverse
, _r_gap
 ;


void cDNAAlignInit (void)
{
  int i ; KEY key ;

  if (_Hit)
    return ;

  i = 256 ;  while (i--) B2[i] = 0 ;
  i = 256 ;  while (i--) B2r[i] = 0 ;
  B2[A_] = 0x0 ; B2[T_] = 0x3 ;     /* ATTENTION copied in dnasubs.c:dnagetWordUsage */
  B2[G_] = 0x1 ; B2[C_] = 0x2 ;     /* you must keep the 2 identical */
  B2r[A_] = 0x3 ; B2r[T_] = 0x0 ;     /* ATTENTION copied in dnasubs.c:dnagetWordUsage */
  B2r[G_] = 0x2 ; B2r[C_] = 0x1 ;     /* you must keep the 2 identical */

  lexaddkey ("Annotation", &key, _VMainClasses) ; 
  _VAnnotation = KEYKEY(key) ;
  lexaddkey ("Clone_Group", &key, _VMainClasses) ; 
  _VClone_Group = KEYKEY(key) ;
  lexaddkey ("Fiche", &key, _VMainClasses) ; 
  _VFiche = KEYKEY(key) ;
  lexaddkey ("Intron", &key, _VMainClasses) ; 
   _VIntron = KEYKEY(key) ; 
  lexaddkey ("Kantor", &key, _VMainClasses) ; 
   _VKantor = KEYKEY(key) ; 
  lexaddkey ("Newname", &key, _VMainClasses) ; 
  _VNewName = KEYKEY(key) ;
  lexaddkey ("RNAi", &key, _VMainClasses) ; 
  _VRNAi = KEYKEY(key) ;
  lexaddkey ("Sage", &key, _VMainClasses) ; 
  _VSage = KEYKEY(key) ;

  lexaddkey ("Transcribed_gene", &key, _VMainClasses) ; 
  _VTranscribed_gene = KEYKEY(key) ;
  lexaddkey ("Transcript", &key, _VMainClasses) ; 
  _VTranscript = KEYKEY(key) ;
  lexaddkey ("cDNA_clone", &key, _VMainClasses) ; 
  _VcDNA_clone = KEYKEY(key) ;
  lexaddkey ("Product", &key, _VMainClasses) ; 
  _VmProduct = KEYKEY(key) ;
  lexaddkey ("mRNA", &key, _VMainClasses) ; 
  _VmRNA = KEYKEY(key) ;

 _Alternative_exon = str2tag ("Alternative_exon") ;
 _Alternative_first_exon = str2tag ("Alternative_first_exon") ;
 _Alternative_intron = str2tag ("Alternative_intron") ;
 _Alternative_partial_exon = str2tag ("Alternative_partial_exon") ;
 _Anomalous_clone = str2tag ("Anomalous_clone") ;
  _Assembled_from_cDNA_clone = str2tag ("Assembled_from_cDNA_clone") ;
 _At_position_1 = str2tag ("At_position_1") ;
 _Bad_quality = str2tag ("Bad_quality") ;
 _Begin_not_found = str2tag ("Begin_not_found") ;
 _Best_product = str2tag ("Best_product") ;
 _Blastp = str2tag ("Blastp") ;
 _Blastp_title = str2tag ("Blastp_title") ;
 _COOH_complete = str2tag ("COOH_complete") ;
 _CTF_File = str2tag ("CTF_File") ;
 _Clone_Group = str2tag ("Clone_Group") ;
 _Coding_gap = str2tag ("Coding_gap") ;
 _Coding_length = str2tag ("Coding_length") ;
 _Complete = str2tag ("Complete") ;
 _Confirmed_exon = str2tag ("Confirmed_exon") ;
 _Confirmed_intron = str2tag ("Confirmed_intron") ;
 _Constructed_from = str2tag ("Constructed_from") ;
 _Contamination = str2tag ("Contamination") ;
 _Derived_sequence = str2tag ("Derived_sequence") ;
 _Discarded_cDNA = str2tag ("Discarded_cDNA") ;
 _Discarded_from = str2tag ("Discarded_from") ;
 _Duplicate_clone = str2tag ("Duplicate_clone") ;
 _EST_translation = str2tag ("EST_translation") ;
 _Exon = str2tag ("Exon") ;
 _Expasy = str2tag ("Expasy") ;
 _Fake_internal_poly_A = str2tag ("Fake_internal_poly_A") ;
 _First_ATG = str2tag ("First_ATG") ;
 _First_NTG = str2tag ("First_NTG") ;
 _First_Kozak = str2tag ("First_Kozak") ;
 _First_exon = str2tag ("First_exon") ;
 _Fmap_cDNA_Decorate = str2tag ("Fmap_cDNA_Decorate") ;
 _Forward = str2tag ("Forward") ;
 _Frame = str2tag ("Frame") ;
 _From_AM = str2tag ("From_AM") ;
 _From_AM_From_prediction = str2tag ("From_AM_From_prediction") ;
 _From_EST = str2tag ("From_EST") ;
 _From_gene = str2tag ("From_gene") ;
 _From_prediction = str2tag ("From_prediction") ;
 _Fuse_to_clone = str2tag ("Fuse_to_clone") ;
 _Fuzzy_gt_ag = str2tag ("Fuzzy_gt_ag") ;
 _Fuzzy_gc_ag = str2tag ("Fuzzy_gc_ag") ;
 _Gap = str2tag ("Gap") ;
 _Gene_wall = str2tag ("Gene_wall") ;
 _Genes = str2tag ("Genes") ;
 _Genomic = str2tag ("Genomic") ;
 _Genomic_sequence = str2tag ("Genomic_sequence") ;
 _Good_product = str2tag ("Good_product") ;
 _Hit = str2tag ("Hit") ;
 _Hits = str2tag ("Hits") ;
 _Ignore_this_clone = str2tag ("Ignore_this_clone") ;
 _In_mRNA = str2tag ("In_mRNA") ;
 _IntMap = str2tag ("IntMap") ;
 _Internal_capping = str2tag ("Internal_capping") ;
 _Internal_priming = str2tag ("Internal_priming") ;
 _Internal_priming_manual = str2tag ("Internal_priming_manual") ;
 _Internal_priming_on_A_rich = str2tag ("Internal_priming_on_A_rich") ;
 _Intron = str2tag ("Intron") ;
 _Intron_boundaries = str2tag ("Intron_boundaries") ;
 _Is_AM = str2tag ("Is_AM") ;
 _Is_chain = str2tag ("Is_chain") ;
 _Is_gap = str2tag ("Is_gap") ;
 _Is_read = str2tag ("Is_read") ;
 _Kantor = str2tag ("Kantor") ;
 _Last_exon = str2tag ("Last_exon") ;
 _Length_3prime_UTR = str2tag ("Length_3prime_UTR") ;
 _Length_5prime_UTR = str2tag ("Length_5prime_UTR") ;
 _Length_anomaly = str2tag ("Length_anomaly") ;
 _LocusLink = str2tag ("LocusLink") ;
 _Longest_cDNA_clone = str2tag ("Longest_cDNA_clone") ;
 _Longest_CDS = str2tag ("Longest_CDS") ;
 _Manual_internal_deletion = str2tag ("Manual_internal_deletion") ;
 _Manual_no_internal_deletion = str2tag ("Manual_no_internal_deletion") ;
 _Manual_polyA = str2tag ("Manual_polyA") ;
 _Mosaic = str2tag ("Mosaic") ;
 _NH2_complete = str2tag ("NH2_complete") ;
 _Nb_alternative_exons = str2tag ("Nb_alternative_exons") ;
 _Nb_confirmed_alternative_introns = str2tag ("Nb_confirmed_alternative_introns") ;
 _Nb_confirmed_introns = str2tag ("Nb_confirmed_introns") ;
 _Nb_possible_exons = str2tag ("Nb_possible_exons") ;
 _Nb_predicted_exons = str2tag ("Nb_predicted_exons") ;
 _Nb_stolen_exons = str2tag ("Nb_stolen_exons") ;
 _NewName = str2tag ("NewName") ;
 _No_obvious_phenotype = str2tag ("No_obvious_phenotype") ;
 _ORF_Gap = str2tag ("ORF_Gap") ;
 _Open_length = str2tag ("Open_length") ;
 _Other = str2tag ("Other") ;
 _PCR_product_size = str2tag ("PCR_product_size") ;
 _Partial_exon = str2tag ("Partial_exon") ;
 _Pfam = str2tag ("Pfam") ;
 _PolyA_after_base = str2tag ("PolyA_after_base") ;
 _Predicted_Exon = str2tag ("Predicted_Exon") ;
 _Predicted_exon = str2tag ("Predicted_exon") ;
 _Probe_hit = str2tag ("Probe_hit") ;
 _Probe_exact_hit = str2tag ("Probe_exact_hit") ;
 _Primed_on_polyA = str2tag ("Primed_on_polyA") ;
 _Problem = str2tag ("Problem") ;
 _Protein = str2tag ("Protein") ;
 _Psort = str2tag ("Psort") ;
 _RNA_editing = str2tag ("RNA_editing") ;
 _RT_PCR = str2tag ("RT_PCR") ;
 _Read = str2tag ("Read") ;
 _Real_3prime = str2tag ("Real_3prime") ;
 _Real_5prime = str2tag ("Real_5prime") ;
 _Ref_Seq = str2tag ("Ref_Seq") ;
 _Ref_mRNA = str2tag ("Ref_mRNA") ;
 _RefSeqMaker = str2tag ("RefSeqMaker") ;
 _RefSeqMaker = str2tag ("RefSeqMaker") ;
 _Relabelled = str2tag ("Relabelled") ;
 _Representative_product = str2tag ("Representative_product") ;
 _Resequence = str2tag ("Resequence") ;
 _Reverse = str2tag ("Reverse") ;
 _Reversed_by = str2tag ("Reversed_by") ;
 _Sage = str2tag ("Sage") ;
 _Spliced_sequence = str2tag ("Spliced_sequence") ;
 _Splicing = str2tag ("Splicing") ;
 _Stolen_Exon = str2tag ("Stolen_Exon") ;
 _Stolen_exon = str2tag ("Stolen_exon") ;
 _Stolen_intron = str2tag ("Stolen_intron") ;
 _Suspected_internal_deletion = str2tag ("Suspected_internal_deletion") ;
 _Target = str2tag ("Target") ;
 _Taxblast = str2tag ("Taxblast") ;
 _Tiling_error = str2tag ("Tiling_error") ;
 _Total_gap_length = str2tag ("Total_gap_length") ;
 _Total_intron_length = str2tag ("Total_intron_length") ;
 _Total_length = str2tag ("Total_length") ;
 _Transcribed_from = str2tag ("Transcribed_from") ;
 _Transcribed_gene = str2tag ("Transcribed_gene") ;
 _Transpliced_to = str2tag ("Transpliced_to") ;
 _UTR_3prime = str2tag ("UTR_3prime") ;
 _UTR_5prime = str2tag ("UTR_5prime") ;
 _Use_AM = str2tag ("Use_AM") ;
 _Vector_clipping = str2tag ("Vector_clipping") ;
 _Very_good_product = str2tag ("Very_good_product") ;
 _aForward = str2tag ("aForward") ;
 _aReverse = str2tag ("aReverse") ;
 _cDNA_clone = str2tag ("cDNA_clone") ;
 _ct_ac = str2tag ("ct_ac") ;
 _gc_ag = str2tag ("gc_ag") ;
 _gt_ag = str2tag ("gt_ag") ;
 _mForward = str2tag ("mForward") ;
 _mReverse = str2tag ("mReverse") ;
 _r_gap = str2tag ("r_gap") ;
  /*   _ = str2tag ("") ;  */
}
