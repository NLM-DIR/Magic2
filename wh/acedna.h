/*  File: dna.h
 *  Author: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1990
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
 * Last edited: feb 2006: mieg , split out from dna.h
 * * May 24 10:26 1999 (edgrif): Add some missing extern func decs.
 * * Jul 23 14:41 1998 (edgrif): Remove redeclarations of fmap functions.
 * * Jul 23 09:06 1998 (edgrif): I have added this header, which was
 *      missing (that's why it shows me as author). I have removed the
 *      declaration of fMapFindDNA which is now in fmaps public header.
 * Created: around 1990, mieg. durbin
 *-------------------------------------------------------------------
 */

/* $Id: acedna.h,v 1.28 2019/11/25 22:39:53 mieg Exp $ */

#ifndef DEF_ACEDNA_H
#define DEF_ACEDNA_H

#include "dna.h"
#include "ac.h"
/*

We use the standard UPAC coding for representing DNA with a single character
per base:

Exactly known bases

A
T  (or U for RNA)
G
C

Double ambiguity

R	AG	~Y		puRine
Y	CT	~R		pYrimidine
M	AC	~K		aMino
K	GT	~M		Keto
S	CG	~S		Strong
W	AT	~w		Weak

Triple ambiguity

H	AGT	~D		not G	
B	CGT	~V		not A
V	ACG	~B		not T
D	AGT	~H		not C

Total ambiguity

N/X	ACGT	~N		unkNown

Possible existence of a base

-       0       padding

Note that:
We do NOT use U internally, but just use it for display if tag RNA is set

NCBI sometimes uses X in place of N, but we only use X for peptides

The - padding character is used when no sequencing is available and
therefore means that the exact length is unknown, 300 - might not mean
exactly 300 unkown bases.

In memory, we use one byte per base, the lower 4 bits correspond to
the 4 bases A_, T_, G_, C_. More tha one bit is set in the ambiguous cases.
The - padding is represented as zero. In some parts of the code, the upper 
4 bits are used as flags.
*/


#define A_ 1
#define T_ 2
#define U_ 2
#define G_ 4
#define C_ 8

#define R_ (A_ | G_)
#define Y_ (T_ | C_)
#define M_ (A_ | C_)
#define K_ (G_ | T_)
#define S_ (C_ | G_)
#define W_ (A_ | T_)

#define H_ (A_ | C_ | T_)
#define B_ (G_ | C_ | T_)
#define V_ (G_ | A_ | C_)
#define D_ (G_ | A_ | T_)

#define N_ (A_ | T_ | G_ | C_)

/* indicator of forward and reverse sequence */
#define FS_ 16
#define RS_ 17

/* indicator for non compatible paired end sequences */
#define NON_COMPATIBLE_PAIR -15

typedef struct  { int dcp, dcq, lng, ok ; } JUMP ;

typedef 
  enum {
    AMBIGUE, INSERTION_TRIPLE, TROU_TRIPLE, INSERTION_DOUBLE, TROU_DOUBLE,
    INSERTION, TROU, ERREUR, TYPE80
       } ALIGN_ERROR_TYPE ;


typedef struct dnaAlignErrStruct 
    { ALIGN_ERROR_TYPE type ;
      int iLong, iShort ;   /* offset in "long" reference dna and "short" tested dna */
      short sens ; /* 1 if short and long are in the same way, otherwise - 1 */
      char baseLong, baseShort ;  /* A 0 for a missing base, or its value */
      }  A_ERR ;

void showDna (Array dna, int n) ; /* for debugging purpose */

extern char dnaDecodeChar[] ;	/* this is the mapping used to decode a single base */
extern char dnaDecodeExtendedChar[] ;	/* this is the mapping used to decode a single base */
extern char dnaEncodeChar[] ;	/* this is the mapping used to encode a single base */
extern char rnaDecodeChar[] ;

#if defined(ARRAY_CHECK)
char checkComplementBase (char cc) ;
#define complementBase(_x) checkComplementBase(_x)
#else
extern char hardComplementBase[] ;	/* complement of single base A_ T_ G_ ... */
#define complementBase(_x) hardComplementBase[(int)(_x)]
#endif

                                /* complement a decoded letter   a t g c */
#define complementLetter(_x) ace_upper(dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)(_x)])]) 
#define complementRnaLetter(_x) ace_upper(rnaDecodeChar[(int)complementBase[(int)dnaEncodeChar[(int)(_x)]]]) 
extern char *aaName[] ;		/* maps single letter code to full name */
char *dnaDecodeString(char *cp) ; /* Decodes up to 255 char in a static buffer */
void dnaEncodeString(char *cp) ; /* reversed encoding, works en place */
void dnaDecodeArray (Array a) ;	/* works in place */
void dnaDecodeExtendedArray (Array a) ;	/* same, but rd added a . special to alignDumpKey */
void rnaDecodeArray (Array a) ;	/* works in place */
void dnaEncodeArray (Array a) ;	/* works in place */
void dnaSolidEncodeArray (Array a, BOOL isDown) ;	/* works in place, transforms ATGC DNA into ATGC transitions, all other bases are coded as repeats of previous base  */
void dnaSolidDecodeArray (Array a, BOOL isDown) ;	/* works in place, transforms back ATGC transitions into standard DNA */

/* Array dnaHandleCopy (dna,h) and dnaCopy(Array dna) ;  // defined as a macro in ac.h->array.h  */

void  reverseComplement(Array dna) ; /* acts in place */

				/* Handles for dnacpt package */
int dnaPickMatch (const unsigned char *cp, int n, const unsigned char *tp, int maxError, int maxN) ;
				/* FastA dumping routines */

BOOL dnaDumpFastA (Array a, int from, int to, char *text, FILE *fil, Stack s) ;
BOOL dnaDoDump (Array dna, int debut, int fin, int offset) ;
Array dnaParseLevel (int level, unsigned char *c1p, char **seqNamep, char **commentp, AC_HANDLE h) ;

float oligoTm (Array dna, int x1, int x2, float *GC_rate) ; /* Maniatis formula */

/* Bp equivalent of the complexity of this oligo, if minEntopy > 0 return 0 if below  */
int oligoEntropy (unsigned const char *dna, int ln, int minEntropy) ;

 /* if TRUE, use the jumper adapted to the sequencing technology 
  * Ilm is the default 
  */
void aceDnaSetIlmJumper (BOOL ok) ;   /* Illumina */
void aceDnaSetSolidJumper (BOOL ok) ; /* LIF/SolID transition csfasta DNA */
void aceDnaSetRocheJumper (BOOL ok) ;
void aceDnaSetPacBioJumper (BOOL ok) ;
void aceDnaSetNanoporeJumper (BOOL ok) ;
void aceDnaSetEditGenomeJumper (BOOL ok) ;

Array aceDnaDoubleTrackErrors (Array  dna1, int *x1p, int *x2p, BOOL isDown,
			       Array dna2, Array dna2R, int *a1p, int *a2p, 
			       int *NNp, Array err, int maxJump, int maxError, BOOL doExtend, int *maxExactp) ;

Array aceDnaTrackErrors (Array  dna1, int pos1, int *pp1, 
			 Array dna2, int pos2, int *pp2, 
			 int *NNp, Array err, 
			 int maxJump, int maxError, BOOL doExtend, int *maxExactp
			 , BOOL isErrClean) ; /* err was just created and is already set t zero */

void newLocalCptErreur(Array longDna, int xl1, int xl2, int pol,
		    Array shortDna, int xs1, int xs2, int pos, int sens,
		    int *NNp, int *startp, int *stopp, int *recouvp, Array errArray) ;
void aceDnaShowErr (Array err) ;

/* the text values are filled in dnasubs.c

G	Guanine
A	Adenine 
T	Thymine
C	Cytosine
R	Purine			A,G
Y	Pyrimidine		T,C
M	Amino			A,C
K	Ketone			G,T
S	Strong interaction	C,G
W	Weak interaction	A,T
H	not-G			A,C,T
B	not-A			G,C,T
V	not-T			A,C,G
D	not-C			A,G,T
N	any                     G,A,T,C

*/

/* BAM/SAM cigar format */
typedef struct cigaretteStruct {
  int a1, a2, x1, x2, type, dx ;
} SAMCIGAR ;

BOOL samParseCigar (char *cigar, Array cigarettes, int a1, int *a2p, int *x1p, int *x2p, int *alip, int *insp, int *delp, int *bigDelp) ;
BOOL samCheckCigar (const char * readName, char *cigar, Array cigarettes, int a1, int ln) ;
int  samScoreCigar (const char * readName, char *cigar, Array cigarettes, int a1, int ln) ;

int samFileCheckCigars (ACEIN ai, const char *target) ; /* ceck all cigars in a sam file */
int samFileExportIntrons (ACEIN ai, const char *target, DICT *targetDict, DICT *intronDict, KEYSET intronSupport) ;
int samFileScore (ACEIN ai, ACEOUT ao, const char *target, DICT *targetDict, Array targetDnas, DICT *intronDict, KEYSET intronSupport) ;

#endif
