#ifndef SA_PARSE_H
#define SA_PARSE_H

/* May 14, 2026
 * New struct to parse the reads
 * the concept is to put all reads of a BB block together in minimal contiguous space
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH

 * This code is public.
 */

/* Define 4 large buffers, records, id, dna, quality (see SAPARSE definition below)
 *  possibly we need their size to reuse them
 *  Allocate the 4 buffers as size 2^n, n >= 24 (4M) to facilitate recycling
 *    The records buffer, 4 u_int per record, contains offsets in the other buffers.
 *    The dnaBuffer should hold approximatelly BMAX MegaBases
 *    The idBuffer is probably much smaller, but identifier lengths vary a lot with samples
 *    The quality buffer is optional, skip it fasta mode, fill it in fastq mode.
 *    Stop filling when the dnaBuffer is full, rather than reallocating
 *
 * In C++ SRA case, Allocate the buffers as aligned addresses use
 *      posix_memalign ((void**)&cp, 64, bufferSize) ;
 * The records buffer is homogeneous and contains 4 unsigned int
 *  It must be populated with offsets in the 3 other buffers and dnalength
 *  Start polulating on line zero, report its length in SAPARSE->nRecords
 * Start the idBuffer, dnaBuffer and qualityBuffer by 16 0, so the 3 DnaRecord offsets are never 0
 * The idBuffer is only used in sa.parse.c to create a dictionary
 * The qualityBuffer is only used sequentially at export
 * The dnaBuffer is used every where and all its adresses must be aligned
 *   Each dna line is filled with ascii dna sequence (no check needed, indifferently upper or lower case)
 *   Each dna sequence must be followed by at least 2 zeroes (up to 16)
 *   The dna offset xDna must be divisible by 16: (xDna & 0xf == 0x0)
 *
 * In C local file case:
 *  The idBuffer is not allocated, the bb->dict dictionary is created directly
 *  The dnaBuffer and quality buffer are indentical to the C++ SRA case
 *
 * In both cases:
 *  The record buffer is used as is to flesh the BB->dnaRecords virtual array.
 *  The identifiers are decorated by the target class character, mapped in a dictionary, and dropped. 
 *  A virtual-virtual DNA array is created using only 2 calls to malloc, rather than 2 per dna sequence.
 *    This virtual uses the recod buffer data to point into the dnaBuffer
 *    The bases in dnaBuufer are checked and translated to binary A_,T_,G_,C_ using the dnaEncodeChar[] look-up table.
 *  The optional quality buffer is kept as is, to be used in SAM/BAM exportation.
 *
 * Deallocation:
 *  All buffers are freed when the BB data block processing is completed.
 *  In SRA streaming mode, the buffers are allocated and recycled in wsra C++ code
 *  In local file mode, they are allocated on BB->handle and freed implicitely
 */

/* usage: DnaRecord *r ;
 *   unsigned char *cp = idBuffer + r->xId >> 1 ;   is the dna identifer
 *   unsigned char *cp = dnaBuffer + r->xDna ;      is the ATGC sequence
 *   unsigned char *cp = qualityBuffer + r->xQual ; is the quality factors
 */
 
typedef struct dnaRecordStruct {
  unsigned int xId ;   // offset in idBuffer, offset <<1 | (isMate ? 0x1 : 0x0)
  unsigned int xDna ;  // offset in dnaBuffer, divisible by 16
  unsigned int xQual ; // offset in idBuffer
  unsigned int dnaLn ;         // length of the dna sequence
} DnaRecord ;

typedef struct saParseStruct {
  Array dnaArray ;               // dnaBuffer = arrp (dnaArray, 0, unsigned char) 
  unsigned char *idBuffer ;      // all ids, packed, zero terminated, a pair shares a single id
  unsigned char *dnaBuffer ;     // all dnas, aligned, 2 to 16 terminal 0
  unsigned char *qualityBuffer ; // all quals, packed, zero terminated
  DnaRecord *records ;

  unsigned int nRecords ;          // usage: for (i=0;i<nRecords;i++)
  unsigned int isPaired ;           // 0: single end,  1: paired end, discovered from NCBI/SRA

  unsigned int idBufferSize ;      // optional, only if useful for C++ recycling
  unsigned int dnaBufferSize ;
  unsigned int qualityBufferSize ;


} SAPARSE ;

/* notice that read pairs are stored in a single SRAPARSE obj
 * preferably but not necessarily interleaved (memory optimization)
 */

/* public functions */
SAPARSE *sraParseGet (const char *srrId, int BMAX, int format) ; // format = 0/1  (fasta/fastq)
void sraParseClose (SAPARSE *saParse) ;  // free or recycle the buffers


#endif
