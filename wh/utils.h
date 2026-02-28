/*  File: utils.h
 *  Author: Rob Clack (rnc@sanger.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 2003
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
 * Description: header file for utils.c - general utility functions
 * HISTORY:
 * Last edited: Jun 30 11:28 2003 (rnc)
 * Created: Fri May 2nd, 2003 (rnc)
 * CVS info:   $Id: utils.h,v 1.21 2019/01/31 16:01:47 mieg Exp $
 *-------------------------------------------------------------------
 */

#ifndef UTILS_DEF
#define UTILS_DEF

BOOL isExecutableInPath (const char *name) ;

/* case-insensitive version of strstr */
char *strcasestr (const char *str1, const char *str2);
/* Wrapper to allow lexstrcmp to be used with arraySort.y */
int arrstrcmp (const void *s1, const void *s2);
/* lextstrcmp does an intuitive sort of strings that include numbers */
int lexstrcmp (const char *a, const char *b);
int lexstrCasecmp (const char *a, const char *b);

/* utility to find and consume an argument on the unix command line */
/* call getArg  (..,..,"-this") to consume: -this */
/* call getArgV (..,..,"-this") to consume: -this value */

BOOL getArg (int *argcp, const char **argv, const char *arg) ; /* arg exists */
const char* getArgV (int *argcp, const char **argv, const char *arg) ; /* arg and value found */
int pickMatch (const char *cp, const char *tp) ;
int pickMatchCaseSensitive (const  char *cp,const char *tp, BOOL caseSensitive) ;

/* utilities for natural ordering od rrays of numbers */
int intOrder (const void *a, const void *b) ; /* used to sort an int array */
int unsignedIntOrder (const void *a, const void *b) ; /* used to sort an int array */
int unsignedIntReverseOrder (const void *a, const void *b) ; /* used to sort an int array */
int floatOrder (const void *a, const void *b) ; /* used to sort an int array */
int doubleOrder (const void *a, const void *b) ; /* used to sort an int array */

/* statistics */
/* chi-square statistics, computed with the -1/2 corection for small mumbers */
int chi2 (int a1, int a2, int b1, int b2, float *c2p) ;
/* Wicoxon Mann Whitney statistics, i.e. the probablibility that u >= U1 is observed in classes of size n1, n2
 * U1 is the unscaled AUC : u = n1 * n2 * auc
 * copy the exact value in *pValue, up to n1 < N and n2 < N where currently N = 21
 * copy the Gaussin approximation in *pValueG
 * copy the Gaussian  approximation in *pValue if n1 or n2 >= N
 * return TRUE, in case of success
 * return FALSE if the argumant are negative or all zero

 * call wilcoxon(-99, -99, -99, 0, 0) to run a test tabulation  
 */
BOOL wilcoxon (int U1, int n1, int n2, double *pValuep, double *pGaussValuep) ;

/* log(n!) approximated via Striling formula when n > 1000 */
double utLogFac (int n) ;

/* arrondis */
int utArrondi (float x) ;
int utMainPart (float p) ;
int utMainRoundPart (float p) ;
int utUpperRoundPart (float p) ;
double utDoubleMainPart (double p) ;

void bufferToUpper (char *buf) ; 
void bufferToLower (char *buf) ; 
void bufferSwitchCase (char *buf) ; 

/* oneByte data compressor library
 * see utils.c for a longer descrition
 * The purpose of this library is to obtain a very compressed format 
 * for storing int values up to 10^5 in a single byte with 1% relative accuracy
 */
void oneByteInitialize (int showCode) ;
unsigned char oneByteEncode (unsigned int x) ;
unsigned int oneByteDecode (unsigned char c) ;
int oneByteTest (void) ;

/* Smoothed Bayesian percentage histogram 
 * hh should be initialized as an array(101, double)
 * and zeroed before a new distrib is computed
 * the function will accumulate in hh a 'Gaussain' of surface one for each pair
 * m black balls observed in n draws 
 */
int utSmoothHisto (Array hh, int m, int n, int mult) ;
void utSmoothHistoTest (void) ;  /* show some examples */

/* multiplicity in a fastc style sequence name xxxx#421a200b21h100 */
int fastcMultiplicity (const char *ccp, int *mult, int multMax) ;

int fetch_and_add (int *variable, int value) ; /* x86 specific multithread-safe addition */


/* REGEXP package */

typedef struct regExpStruct *RegExp ;

RegExp regExpCreate (const char *pattern, AC_HANDLE h) ;
int regExpFind (RegExp br, const char *data) ;
int regExpMatch (const char *data, const char *pattern) ; 

/* Usage
 * br = regExpCreate (pattern, 0, h) ;
 *   // 0: bad pattern, !0 : reusable query handle 
 * n = regExpFind (br, "data") ;
 *   // 0: not found, n > 0 : found at position, pos = (n-1) if getPos was set
 * n = regExpMatch (data, pattern, 0) ;
 *   // -2: data NULL or empty, -1: bad pattern, 0: not found, >0: found a position
 *   //    all in one call, for a non reusable query
 * ac_free (br) ;
 */

#endif

/*********************** end of file ********************************/
