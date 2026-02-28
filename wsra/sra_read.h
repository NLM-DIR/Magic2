/** A library that downoads NGS runs from SRA using NCBI NGS SDK:
 *    https://github.com/ncbi/sra-tools/wiki/09.-Downloading-NGS-SDK
 */

#ifndef _SRA_READ_H
#define _SRA_READ_H

#ifdef __cplusplus
extern "C" {
#endif

/* A type for opaque pointer to a C++ structure that holds read iterator for
   an SRA run and a string that serves as a buffer for downloaded reads in
   FASTA format */
typedef void SRAObj;

/* Create a new SRAObj object for an SRA run */
SRAObj* SraObjNew(const char* accession);

/* Free SRAObj object */
SRAObj* SraObjFree(SRAObj* sra);

/* Download a batch of reads with num_bases bases and return a pointer to
   to reads in FASTA or FASTQ format. The function downloads complete reads,
   so in most cases the number of downloaded bases will be just above
   num_bases. If quality_scores == 0, then the downloaded reads will be
   in the FASTA format, otherwise FASTQ with quality scores. */
const char* SraGetReadBatch(SRAObj* sra, int num_bases, int quality_scores);


#ifdef __cplusplus
}
#endif

#endif
