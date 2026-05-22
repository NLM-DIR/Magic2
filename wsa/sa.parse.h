#ifndef SA_PARSE_H
#define SA_PARSE_H

/* May 14, 2026
 * New struct to parse the reads
 * the concept is to put all reads of a BB block together in minimla space
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH

 * This code is public.
 */

/* Define 3 large buffers, possibly we need their size to reuse them
 * In SRA streaming mode, they are allocated and recycled in wsra C++ code
 * In local file mode, they are handled in wsa/sa.parse.c 
 *
 * Allocate the buffers as aligned addresses use
 *      posix_memalign ((void**)&cp, 64, bufferSize) ;
 * Start each buffer by 16 0, so the 3 DnaRecord offsets are never 0
 * Allocate the buffers as size 2^n, n >= 24 (4M) to facilitate recycling
 * The idBuffers and qualBuffers are only used sequentially at export
 * The dnaBuffer is used every where and all its adresses must be aligned
 * the dna is filled with ascii dna sequence (indifferently upper or lower case)
 * it will be remapped to binary via a look-up table dnaEncodeChar later.
 * each dna sequence must be followed by at least 2 zeroes (up to 16)
 * Its offset must be divisible by 16: (dna & 0xf == 0x0)
 */

/* usage: DnaRecord *up ;
 *   unsigned char *cp = idBuffer + (up->id >> 1) is the name of the pair 
 *   up->id & 0x1 == 0x0 : is first read in a pair
 *   up->id & 0x1 == 0x1 : is second read in a pair
 *   unsigned char *cp = dnaBuffer + up->dna : is the ATGC sequence
 *   unsigned char *cp = qual : is the quality factors
 */
 
typedef struct dnaRecordStruct {
  unsigned int xId ;   // offset in idBuffer, offset <<1 | (isMate ? 0x1 : 0x0)
  unsigned int xDna ;  // offset in dnaBuffer, divisible by 16
  unsigned int xQual ; // offset in idBuffer
  int dnaLn ;         // length of the dna sequence
} DnaRecord ;

typedef struct saParseStruct {
  unsigned char *idBuffer ;      // all ids, packed, zero terminated, a pair shares a single id
  unsigned char *dnaBuffer ;     // all dnas, aligned, 2 to 16 terminal 0
  unsigned char *qualityBuffer ; // all quals, packed, zero terminated
  DnaRecord *records ;

  unsigned int idBufferSize ;      // sizes may be useful for recycling
  unsigned int dnaBufferSize ;
  unsigned int qualityBufferSize ;
  unsigned int nRecords ;          // usage: for (i=0;i<nRecords;i++)

  unsigned int minDnaLn, maxDnaLn, nBases ;
  unsigned int isPaired ;           // 0: single end reads,  1: paired end reads
} SAPARSE ;

/* notice that read pairs are stored in a single SRAPARSE obj
 * preferably but not necessarily interleaved (memory optimization)
 */

SAPARSE *sraParseGet (const char *srrId, int BMAX /* , format info */ ) ;
void sraParseClose (SAPARSE *saParse) ;  // free or recycle the buffers


#endif
