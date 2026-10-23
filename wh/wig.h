/*
 *  File wig.h  : 2026 new header file for ACEDB wiggles (density profiles over a sequence area)
 *  Author: J Thierry-Mieg, 2010
 *  Copyright (C) J Thierry-Mieg, 2010
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
 * This file is part of the ACEVIEW genome database package, written by
 *	Jean Thierry-Mieg (NCBI) mieg@ncbi.nlm.nih.gov
 *
 *-------------------------------------------------------------------
 */


#ifndef ACEDB_WIG_DEF
#define ACEDB_WIG_DEF

#include "../wac/ac.h"
#include "dict.h"
#include "aceio.h"

typedef struct wigStruct {
  AC_HANDLE h ;
  int magic ;
  int step ;
  int scale ;  /* divide all values by scale if exporting */
  int posMin ; /* coordinate of offSet 0, divisible by step */
  int iMin ; /* first offset with non zero value */
  int iMax ; /* last  offset with non zero value */
  int shift ;  /* used in makeGroupWiggles and in findPeaks, divisible by step */
  Array aa ;   /* uint every step starting */
  long int nBp ;
  unsigned int median, medianNoZero, minCover ;   /*  set by wigFindPeaks */ 
  double chromLength ;
  double bases ;
  Array peaks ;
  const char *runName, *chromName ;
} *WIG ;

WIG wigCreate (const char *runName, const char *chromName, int step, AC_HANDLE h) ;
/* wigParseAZ can be called on several files to cumulate the data */
void wigParseAZ (WIG wig, const char *fName, int shift) ;
/* feed in a WIG an array of WSEG alignment segments */
void wigFeed (WIG wig, BigArray wsegs, int scale) ;
void wigCompress (WIG wig, int step) ;
void wigCountGenes (WIG wig, BigArray geneB, BigArray geneC, long int *cumuls, char wigStrand) ;

void wigExportAZ (WIG wig, const char *fName) ;
void wigExportBF (WIG wig, const char *fNam, BOOL gzo) ;
void wigExportBigWig (WIG wig, const char *fName, const char *indexName) ;

void wigExportPeaks (WIG wig, const char *fName, BOOL gzo) ;
int  wigFindPeaks (WIG wig, int minCover) ;     /* minCover = 0 =>automatic setting */
int wigSetShift (WIG wig, int shift) ; /* returns previous value */

/*
 * Disk compact AW file structures
 * In each block, the values are byte splitted creating long flows odf zeroes that compress better
 */

typedef struct azStruct {
  AC_HANDLE h ;
  int magic ;
  int step ;      /* recommended 10 for large genomes */
  char *fName ;
  char *target ;
  int posMin, posMax ;
  int xMin ;      /* minimal position */
  int xMax ;      /* maximal position */
  int NB ;        /* number of files */
  int bMax ;      /* 0x1 << BMAX values in each file, default BMAX=18 => 2560 kilobases if step==10 */
  int cNB ;       /* currently cached  file */
  long int c1, c2 ;    /* min max coordinates [c1,c2] in current cache */
  Array cache ;   /* values in current cache */
  int chromLength ;
  float bases ;   /* total number of bases in wig */
} AZZ ;


AZZ *wigAzOpen (const char *fName, AC_HANDLE h0) ;
BOOL wigAzZone (AZZ *az, int x1, int x2, int shift, WIG wig, Array wPoints, int *nPosp, long int *nBpp, BOOL cumul) ;
AZZ *wigAzWrite (const char *fName, WIG wig, AC_HANDLE h0) ;


typedef struct cstStruct { int magic ; int step, dxMax, bestShift ; long int nu ; double u1, v1, u2, v2, bestUv, *uu, *uv ; } *CST ;
CST cisTransCreate (int dxMax, AC_HANDLE h) ;
void cisTransCumulate (CST genomeCst, CST chromCst) ;
void cisTransFeed (CST cst, int step, WIG waF, WIG waR) ;
void cisTransCumulateArrays (CST cst, int step, Array aaF, Array aaR) ;
void cisTransNormalize (CST cst) ;
void cisTransExport (CST cst, ACEOUT ao) ;


#endif
