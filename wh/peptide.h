/*  File: peptide.h
 *  Author: Richard Durbin (rd@sanger.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
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
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Feb 21 10:38 1997 (il)
 * Created: Wed May 11 01:49:18 1994 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: peptide.h,v 1.8 2015/05/14 01:14:15 mieg Exp $ */ 

#ifndef ACEDB_PEPTIDE_H
#define ACEDB_PEPTIDE_H

extern char pepDecodeChar[] ;	/* [0,21] -> A,C,D... */
extern signed char pepEncodeChar[] ;	/* A,C,D... -> [0,21] : -ve bad */
extern char *pepName[] ; /* A,C,D... -> Alanine, Cysteine... */
extern char *pepShortName[] ;	/* A,C,D... -> Ala, Cys... */
extern int  molecularWeight[] ;  /* A,C,D... -> 89, 121... */

void pepDecodeArray (Array pep) ;
void pepEncodeArray (Array pep) ;


/* k is either a genetic_code object to use for translation. 
 * or a sequence object, in which case we recurse through the
 * parents looking for the genetic_code to be used.
 * 
 * returns a reference to an array for use with e_codon() 
 * and the code that was used
 * for efficiency, the array is NOT a copy but the code itself
 * do NOT free this array, just forget it.
 * returns key of Genetic_code in *geneticCodep if non-NULL. */
char *pepGetTranslationTable (KEY k, KEY *geneticCodep) ;



/*
* codon() and friends 
* assume the standard eucaryot genomic code
* they are kept for backward compatibility
* Prefer e_codon() for any new code.
*/
char codon (const char *cp) ;		/* next 3 bases -> amino acid */
char reverseCodon (const char *cp) ;	/* next 3 bases on opp strand -> amino acid */
char antiCodon (const char* cp) ;

/* extended_codon, needed for non standard genetic codes (e.g. virus, mithocondria)
 * s is dna to translate
 * translationTable is previously obtained from pepGetTranslationTable()
 * returns the (coded) amino acid */
char e_codon (const char *s, const char* translationTable) ; /* next 3 bases -> amino acid */
char e_reverseCodon (const char *s, const char* translationTable) ;
							    /* next 3 bases on opp strand -> amino acid */
char e_antiCodon (const char *s, const char* translationTable) ;

float pepPI (Array pep) ;
int pepWeight (Array pep) ;
				/* main access routines */
Array peptideGet (KEY key) ;
Array e_peptideTranslate(KEY key, BOOL CDS_only, const char *translationTable) ;
Array peptideTranslate(KEY key, BOOL CDS_only) ;
BOOL peptideStore (KEY key, Array pep) ; /* fails if model is missing */
BOOL pepSubClass (KEY protein, KEY *pepKeyp) ; /* gets the parent class of pep, by default a Proteinuence */
BOOL pepReClass (KEY pep, KEY *proteinp) ;  /* gets the parent class of pep, by default a Protein */
BOOL pepInClass (int classe) ; /* does this class supports dna ? */

				/* FastA dumping routines */
BOOL pepDumpFastA (Array a, int from, int to, const char *text, FILE *fil, Stack s) ;
BOOL pepDumpFastAKey (KEY key, FILE *fil, Stack s, char style) ;
int  pepDumpFastAKeySet (KEYSET kSet, FILE *fil, Stack s) ;
				/* if !fil  and !sthese call pepFileOpen */
FILE *pepFileOpen (void) ;
int hashArray(Array a);


/* Used by gifacemain currently.                                             */
void* pepGifGet (void*) ;	/* returns a handle used by other pepGif*() */
void pepGifSeq (void*) ;
void pepGifAlign (void*) ;
void pepGifDestroy (void*) ;


/* THIS SHOULD NOT BE HERE...SIGH...BUT CAN'T PUT IT IN PEPDISP.H BECAUSE    */
/* OF HEADER CLASHES...                                                      */
/* Pep display create data, controls how pep display will be created, passed */
/* in as a void * to pepDisplay() via displayApp()                           */
#ifdef JUNK
typedef struct _PepDisplayData
{
  BOOL CDS_only ;					    /* TRUE: translate only CDS portion of */
							    /* objects DNA.*/
} PepDisplayData ;
#endif


#endif   /* !ACEDB_PEPTIDE_H */
/************** end of file **************/
 
 
