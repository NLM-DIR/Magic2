/*  bgzf.h — Blocked GNU Zip Format (BGZF) reader and GZI index
 *
 *  WHAT IS BGZF:
 *    BGZF is an extension of the gzip format defined in the SAM/BAM
 *    specification (https://samtools.github.io/hts-specs/SAMv1.pdf, Appendix).
 *    A BGZF file is a valid gzip file whose deflate stream is partitioned into
 *    independent compressed blocks of at most 65,535 bytes each.  Because each
 *    block is self-contained, random access is possible: seek to a block
 *    boundary in the compressed file, decompress that block alone, and resume
 *    reading from any uncompressed offset within it.  This makes parallel
 *    processing of large .gz files feasible without decompressing from the
 *    start.  The format is frozen as part of the SAM specification and is
 *    therefore stable indefinitely.
 *
 *  WHAT IS GZI:
 *    A .gzi file is a binary index produced by  bgzip -r file.fasta.gz
 *    It records the compressed and uncompressed byte offset of every BGZF
 *    block boundary, as pairs of little-endian uint64_t values, preceded by
 *    a uint64_t entry count.  Given a target uncompressed offset, a binary
 *    search on the GZI locates the block to seek to.
 *
 *  WHY WE REIMPLEMENT RATHER THAN USE HTSLIB:
 *    This module is a self-contained reimplementation of the BGZF reader and
 *    GZI loader.  We deliberately avoid the htslib dependency because external
 *    libraries introduce version-skew risk: APIs change, bugs are introduced,
 *    and scientific results can silently differ across versions.  The BGZF and
 *    GZI formats are frozen standards; this code will remain correct without
 *    maintenance.  The only external dependency is zlib (RFC 1950/1951), which
 *    is part of the base system on every POSIX platform and has been stable
 *    since 1996.
 *
 *  USAGE SUMMARY:
 *    GZIIndex *idx = bgzfGZILoad (fNam) ;        // load .gzi sidecar
 *    BGZFFile *bgzf = bgzfOpen (fNam) ;           // open .gz file
 *    uint64_t  cOff = bgzfGZIFloor (idx, uOff) ;   // compressed seek point
 *    bgzfSeek (bgzf, cOff) ;                                // seek to block boundary
 *    bgzfRead (bgzf, buf, len) ;                          // decompress into buf
 *    bgzfClose (bgzf) ;
 *    bgzfGZIFree (idx) ;
 *
 *  Author: derived from BGZF spec; integrated into acedb library style.
 *  This file may be copied freely as a self-contained unit.
 */

#ifndef BGZF_H_DEFINED
#define BGZF_H_DEFINED

#include <stdint.h>
#include <stddef.h>
#include <sys/types.h>          /* ssize_t, off_t                           */

/*---------------------------------------------------------------------------
 * GZI index
 *-------------------------------------------------------------------------*/

typedef struct {
    uint64_t cOff ;             /* compressed byte offset of block start    */
    uint64_t uOff ;             /* uncompressed byte offset of block start  */
} GZIEntry ;

typedef struct {
    GZIEntry *entries ;       /* sorted array, entries[0] is always (0,0) */
    size_t    nEntries ;        /* number of entries including (0,0)         */
    uint64_t  uSize ;          /* total uncompressed file size              */
} GZIIndex ;

/* Load the .gzi sidecar file for fNam (i.e. fNam + ".gzi").
 * Returns NULL if the file cannot be opened or is malformed.
 * The implicit entry (cOff=0, uOff=0) is always prepended.
 */
GZIIndex   *bgzfGZILoad  (const char *fNam) ;

/* Free a GZIIndex returned by bgzfGZILoad.                              */
void        bgzfGZIFree  (GZIIndex *idx) ;

/* Return the compressed offset of the BGZF block whose uncompressed
 * start is the largest value <= uTarget.  This is the correct seek
 * point for reading uncompressed data at uTarget.                       */
uint64_t    bgzfGZIFloor (const GZIIndex *idx, uint64_t uTarget) ;

/*---------------------------------------------------------------------------
 * BGZF file handle
 *-------------------------------------------------------------------------*/

#define BGZF_MAX_BLOCK  65536   /* maximum uncompressed block size in bytes */

typedef struct BGZFFile_ BGZFFile ;

/* Open a BGZF-format .gz file for reading.
 * Returns NULL on failure.
 */
BGZFFile   *bgzfOpen     (const char *fNam) ;

/* Close and free a BGZFFile.                                            */
void        bgzfClose    (BGZFFile *bgzf) ;

/* Seek to a compressed byte offset (must be a block boundary from GZI).
 * Returns 0 on success, -1 on error.
 */
int         bgzfSeek     (BGZFFile *bgzf, uint64_t cOffset) ;

/* Decompress and read up to len bytes into buf.
 * Transparently crosses block boundaries.
 * Returns number of bytes read, 0 at EOF, -1 on error.
 */
ssize_t     bgzfRead     (BGZFFile *bgzf, void *buf, size_t len) ;

/* Return current uncompressed position.                                 */
uint64_t    bgzfTell     (const BGZFFile *bgzf) ;

#endif /* BGZF_H_DEFINED */
