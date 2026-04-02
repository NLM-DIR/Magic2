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

/* A struct to return SRA data, must be allocated with SRAReadBatchNew() */
typedef struct SRAReadBatch
{
    /* Internal SRA data */
    SRAObj* sra_obj;

    /* Forward reads in FASTA or FASTQ format (do not deallocate the buffer) */
    const char* seq;

    /* Reverese reads in FASTA or FASTQ format (do not deallocate the buffer) */
    const char* seq2;

    /* Number of bases read in the batch */
    long unsigned int num_bases;

    /* 1 is the SRA run is paired, otherwise 0 */
    int is_paired;
} SRAReadBatch;

/* Create a new SRAReadBatch object for an SRA run */
SRAReadBatch* SRAReadBatchNew(const char* accession);

/* Free SRAReadBatch object */
SRAReadBatch* SRAReadBatchFree(SRAReadBatch* sra);


/* Download a batch of reads with num_bases bases and return a pointers to
   to reads in FASTA or FASTQ format. The function downloads complete reads,
   so in most cases the number of downloaded bases will be just above
   num_bases. If quality_scores == 0, then the downloaded reads will be
   in the FASTA format, otherwise FASTQ with quality scores.
   If split_spot == 0 then both mates of paired reads are downloaded to
   sra->seq, otherwise they are split into sra->seq and sra->seq2. */
int SraGetReadBatch(SRAReadBatch* sra, long int num_bases, int quality_scores,
                    int split_spot);


#ifdef __cplusplus
}
#endif

#endif
