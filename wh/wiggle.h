/*
 *  File wiggle.h  : header file for ACEDB wiggles (density profiles over a sequence area)
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


#ifndef ACEDB_WIGGLE_DEF
#define ACEDB_WIGGLE_DEF

#include "../wac/ac.h"
#include "dict.h"
#include "aceio.h"

typedef enum { BF, BV, BG, AF, AM, AG, AW, AZ, BHIT, COUNT } WFORMAT ;
typedef struct wiggleStruct { 
  AC_HANDLE h ; const char *title ; 
  ACEIN ai ; ACEOUT ao, aoPeaks ;
  Array aaa, aaaCopy ;
  WFORMAT in, out ;

  int in_x1, in_step, in_span, out_step, out_span, delta, gauss ;
  float minCover, maxCover ;
  float wiggleScale ;  /* multiply wp->y by scale when exporting */
  float wiggleScale1 ;  /* multiply wiggle1 wp->y by scale when exporting */
  float wiggleScale2 ;  /* multiply wiggle2 wp->y by scale when exporting and add wiigle 1 and 2 */
  char fileType[16] ;
  const char *nmWiggle, *support, *cosmidWiggle
    , *trackName, *hitFile, *wiggleDir, *wiggleFeature, *ignoredProbeFileName, *mapFileName
    , *sxxWiggleFileName
    , *wiggleFileName1, *wiggleFileName2
    , *swiggleFileName1, *swiggleFileName2   /* controlling stranded wiggles */
    , *selectFileName, *rejectFileName
    , *inFileOfFileList
    , *inFileName, *outFileName, *target_class, *min_target_class, *sxxChromosome
    , *transcriptsEndsFileName ;
  BOOL unique ; /* BHIT format: only consider unique hits */
  BOOL forceUnique ; /* BHIT format: reset multiplicity = 1 */
  BOOL non_unique ; /* BHIT format: only consider unique hits */
  BOOL partial ; /* used in -ventialte option: keep reads with overhand or orphans or incompatible pairs */
  BOOL noPartial ; /* treat partial as if the were complete, useful for nanopore, PacBio */
  BOOL gzi, gzo ; /* gzip the input or the output */
  BOOL noRemap ; /* all coordinates go in aaa[0], target is discarded */
  Array targets ;
  KEYSET map2remap ;
  DICT *dict, *mapDict, *target_mapDict, *remapDict, *selectDict, *rejectDict ;
  BOOL strand, antistrand, stranded, mrnaWiggle ;
  int strandShift_max, minErrRate, minAliRate, minAliLength, pair ;
  int maxErr, maxErrRate ;
  const char *strandShift_f, *strandShift_r ;
  const char *strategy ;
  BOOL RNA_seq ;
  BOOL ventilate, cumul, peaks, hierarchic, lengthCoverage, flagEnds ;
  int wiggleRatio, wiggleRatioDamper ;
  int multiPeaks ;  
  int BF_predictor ; /* degree of the polynome used to compress the BF format */
  int BF_compressor ; /* degree of the polynome used to compress the BF format */
} WIGGLE ;

typedef struct sxwtMapStruct { int target_map, map, x1, x2, remap, a1, a2 ; float y ;
  BOOL isDown ; KEYSET hits, upHits, downHits ;} SXMAP ;

typedef struct wigglePointStruct { int x ; float y ; } WIGGLEPOINT ;

BOOL sxCheckFormat (const char *io, WFORMAT *ip, const char *ccp, char *ftype) ;

void sxWiggleParse (WIGGLE *sx, int z1, int z2) ;
int sxWiggleGauss (WIGGLE *sx) ;
void sxWiggleScale (WIGGLE *sx, float scale) ;
void sxWiggleShift (WIGGLE *sx, float delta) ;
void sxWiggleFloor (WIGGLE *sx, float mini) ;
void sxWiggleEndRatio (WIGGLE *sx) ;
Array sxWiggleStrandedRatio (WIGGLE *sx, Array aaa, Array bbb, int damper, AC_HANDLE h) ;
Array sxWiggleMultiply (WIGGLE *sx, Array aaa, AC_HANDLE h) ;
void sxWiggleMultiplyLocally (WIGGLE *sx, Array bbb, Array aaa1) ;
void sxWiggleCopy (WIGGLE *sx) ;
void sxWiggleExport (WIGGLE *sx) ;

Array sxGetWiggleZone (Array aa, const char *fNam, char *type, int *stepp, const char *chrom, int a1, int a2, AC_HANDLE h) ; /* returns an array of WIGGLEPOINT */


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
} AZZ ;

AZZ *wigAzWrite (const char *fName, const char *target, Array aa, Array wPoints, int step, int posMin, int posMax, int bMax, AC_HANDLE h0) ; /* Write the data */
AZZ *wigAzOpen (const char *fName, AC_HANDLE h) ;    /* Open an AZZ file, return NULL on error. */
BOOL wigAzAt (AZZ *az, int x, unsigned int *value) ; /* Value at position x, FALSE if out of range */
BOOL wigAzZone (AZZ *az, int x1, int x2, Array wPoints, int *nPos, long int *nBpp) ; /* Out of range values are just set to zero */
void azDoClose (AZZ *az) ;     /* free az */
#define azClose(_az) {azDoClose(_az);_az=0}
#endif
