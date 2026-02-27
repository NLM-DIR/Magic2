/*
 * sa.sequenceParallelParse.c

 * This module is part of the sortalign package
 * Created February, 2026

 * This code is public.

 * This module reads and decompresses fasta/fast files
 * relying on the io_uring fast library
 */

#include "sa.h"
#include <linux/io_uring.h>

/**************************************************************/
/**************************************************************/

typedef struct runFilesStruct {
    CHAN *outChan;                 // channel of BB*
    const char *filename;          // R1
    const char *filename2;         // R2 (NULL if !pairedEnd)
    int BMAX;                      // target MB per emitted BB
    BOOL gzi;                      // TRUE = gzipped
    BOOL pairedEnd;
    BOOL isFasta;                  // TRUE = FASTA, FALSE = FASTQ
    AC_HANDLE h;
} RF ;

typedef struct fileStreamStruct {
    int fd;
  /*     struct io_uring ring; */
    char *compBuf[2];
    struct inflate_state *inflate;
    char *decompBuf;
    char *remnantBuf;
    size_t remnantSize;
    off_t offset;
    int bufIdx;
    int pending;
    size_t totalRecords;
    size_t totalBytes;
} URFS ;


/**************************************************************/
/************************* utilities **************************/
/**************************************************************/

static inline void *memrchr(const void *s, int c, size_t n)
{
    const unsigned char *p = (const unsigned char *)s + n ;
    while (p > (const unsigned char *)s)
      {
        --p;
        if (*p == (unsigned char)c)
	  return (void *)p ;
      }
    return NULL;
} /* memrchr */

/**************************************************************/
/* Locate last raw record in buffer */
static int getLastRawRecord (char *start, char *end)
{
  for (char *p = end ; p >= start ; p--)
    if (*p == '\n')
      return p + 1 - start ;
  return -1 ; /* failed */
} /* getLastRawRecord */

/**************************************************************/
/* Locate last fasta/fastc record in buffer */
static unsigned char  *getLastFastaRecord (unsigned char *start, unsigned char *end)
{
  unsigned char *p = end ;
  while (1)
    {
      unsigned char *p = memrchr (end, '>', end - start) ;
      if (p &&  (p == start || p[-1] == '\n'))
	break ;
    } ;
  return p ; 
} /* getLastFastaRecord */

/**************************************************************/
/* Locate last SRA pair record in buffer */
static int getLastSraPairRecord (char *start, char *end)
{
  for (char *p = end ; p >= start ; p--)
    if (*p == '>')
      {
	if (p == start)
	  return 0 ;
	char *cp = strchr (p, '\n') ;
	if (cp)
	  {
	    *cp = 0 ;
	    char *cq = strstr (p, ".1\n") ;
	    *cp = '\n' ;
	    if (!cq)  /* second member of the pair, iterate */ 
	      continue ;
	    return p - start ;
	  }
      }
  return -1 ; /* failed */
} /* getLastSraPairRecord */

/**************************************************************/
/* Locate last fastq record in buffer */
static unsigned char *getLastFastqRecord (unsigned char *start, unsigned char *end)
{
  for (unsigned char *p = end ; p >= start ; p--)
    if (*p == '@')
      {
	if (p == start)
	  return 0 ;
	unsigned char *cp = memchr (p, '\n', end - start) ;
	if (cp)
	  { /* check is we see dna */
	    BOOL ok = TRUE ;
	    for (unsigned char *cq = cp + 1 ; ok && *cq && *cq != '\n' ; cq++)
	      if (! dnaEncodeChar[(int)*cq])
		ok = FALSE ;
	    if (! ok) /* probably the  quality line */
	      continue ;
	    return p ;
	  }
      }
  return NULL ; /* failed */
} /* getLastFastqRecord */

/**************************************************************/
/**************************************************************/
/* Submit one read to the io_uring decompressor */
static void submitRead (URFS *fs, size_t want)
{
#ifdef JUNK
  struct io_uring_sqe *sqe = io_uring_get_sqe (&fs->ring) ;
  io_uring_prep_read (sqe, 0
		      , fs->compBuf[fs->bufIdx]
		      , want, fs->offset
		      ) ;
  sqe->user_data = (unsigned long)fs->bufIdx ;
  io_uring_submit (&fs->ring) ;
  fs->offset += want ;
  fs->bufIdx ^= 1 ;
  fs->pending++ ;
#endif
} /* submitRead */

/**************************************************************/
/* Process one completed read for a file */
static void process_file_completion (URFS *fs, BOOL isFasta, size_t targetSize, BB **outBB, size_t *outRecords)
{
#ifdef JUNK
  struct io_uring_cqe *cqe ;
  io_uring_wait_cqe (&fs->ring, &cqe) ;
  size_t nread = cqe->res ;
  char *raw = fs->compBuf[ (int)cqe->user_data] ;
  io_uring_cqe_seen (&fs->ring, cqe) ;
  fs->pending-- ;
  
  if  (nread == 0) return ;
  
  // Decompress or direct copy
  size_t produced = 0 ;
  if  (fs->inflate)
    {
      fs->inflate->next_in   = raw ;
      fs->inflate->avail_in  = nread ;
      fs->inflate->next_out  = fs->decompBuf + fs->remnantSize ;
      fs->inflate->avail_out = /* decompBuf size */ - fs->remnantSize ;
      
      isal_inflate (fs->inflate) ;
      produced = /* decompBuf size */ - fs->remnantSize - fs->inflate->avail_out ;
    }
  else
    {
      memcpy (fs->decompBuf + fs->remnantSize, raw, nread) ;
      produced = nread ;
    }
  
  size_t totalAvail = fs->remnantSize + produced ;
  
  // Find last complete record
  char *lastRecordStart = NULL ;
  size_t newRecords = 0 ;
  if  (isFasta)
    {
      //  lastRecordStart = find_last_fasta_record (fs->decompBuf, fs->decompBuf + totalAvail - 1) ;
       lastRecordStart = getLastFastaRecord ((fs->decompBuf, fs->decompBuf + totalAvail - 1) ;
    }
  else
    {
      lastRecordStart = getLastFastqRecord (fs->decompBuf, fs->decompBuf + totalAvail - 1, &newRecords) ;
    }

  size_t completeLen = lastRecordStart ?  (lastRecordStart - fs->decompBuf) : 0 ;

  // Emit if we have enough data
  if  (completeLen >= targetSize)
    {
      BB *bb = bbNew (fs->h) ;  // your bbNew function
      bb->gzBuffer = halloc (completeLen + 1, bb->h) ;
      memcpy (bb->gzBuffer, fs->decompBuf, completeLen) ;
      bb->gzBuffer[completeLen] = 0 ;
      bb->nSeqs = 100 ;   // your guess
      // set other fields as in original code
      
      *outBB = bb ;
      *outRecords = newRecords ;
    }

  // Save remnant
  fs->remnantSize = totalAvail - completeLen ;
  if  (fs->remnantSize > 0) 
    memmove (fs->remnantBuf, fs->decompBuf + completeLen, fs->remnantSize) ;
  else 
    fs->remnantSize = 0 ;
#endif

} /* process_file_completion */

/**************************************************************/
/* Main agent - fully functional */
void saUringSequenceParser (const PP *pp, RC *rc, TC *tc, BB *bb)
{
  getLastRawRecord (0, 0) ;
  getLastFastaRecord (0, 0) ;
  getLastSraPairRecord (0, 0) ;
  getLastFastqRecord (0, 0) ;
#ifdef JUNK
    AsyncReaderData *rd = arg ;
    AC_HANDLE h = rd->h ;

    // Open files
    int fd1 = open (rd->filename, O_RDONLY | O_DIRECT | O_CLOEXEC) ;
    int fd2 =  (rd->pairedEnd && rd->filename2) ? open (rd->filename2, O_RDONLY | O_DIRECT | O_CLOEXEC) : -1 ;
    if  (fd1 < 0) messcrash ("cannot open %s", rd->filename) ;
    if  (rd->pairedEnd && fd2 < 0) messcrash ("cannot open %s", rd->filename2) ;

    struct stat st1, st2 ;
    fstat (fd1, &st1) ;
    if  (rd->pairedEnd) fstat (fd2, &st2) ;

    // Initialize streams
    URFS fs1, fs2 ;
    init_file_stream (&fs1, fd1, rd->gzi,  (rd->BMAX << 20) +  (8 << 20), h) ;
    if  (rd->pairedEnd) {
        init_file_stream (&fs2, fd2, rd->gzi,  (rd->BMAX << 20) +  (8 << 20), h) ;
    }

    size_t records1 = 0, records2 = 0 ;

    while  ( (fs1.offset < st1.st_size || fs1.pending > 0 || fs1.remnantSize > 0) &&
            (!rd->pairedEnd ||  (fs2.offset < st2.st_size || fs2.pending > 0 || fs2.remnantSize > 0))) {

        // Submit reads
        while  (fs1.pending < 8 && fs1.offset < st1.st_size) {
            size_t want = 8 << 20 ;
            if  (fs1.offset + want > st1.st_size) want = st1.st_size - fs1.offset ;
            submit_read (&fs1, want) ;
        }
        if  (rd->pairedEnd) {
            while  (fs2.pending < 8 && fs2.offset < st2.st_size) {
                size_t want = 8 << 20 ;
                if  (fs2.offset + want > st2.st_size) want = st2.st_size - fs2.offset ;
                submit_read (&fs2, want) ;
            }
        }

        // Process completions
        BB *bbR1 = NULL ; size_t recR1 = 0 ;
        process_file_completion (&fs1, rd->isFasta, rd->BMAX << 20, &bbR1, &recR1) ;

        BB *bbR2 = NULL ; size_t recR2 = 0 ;
        if  (rd->pairedEnd) {
            process_file_completion (&fs2, rd->isFasta, rd->BMAX << 20, &bbR2, &recR2) ;
        }

        // Emit when both sides are ready and synced
        if  (rd->pairedEnd) {
            if  (bbR1 && bbR2 && recR1 == recR2 && recR1 > 0) {
                BB *bb = bbNew (h) ;
                bb->gzBuffer = bbR1->gzBuffer ;
                bb->gzBuffer2 = bbR2->gzBuffer ;
                bb->nSeqs = recR1 ;
                bb->nSeqsR2 = recR2 ;
                channelPut (rd->outChan, &bb, BB*) ;
                records1 += recR1 ;
                records2 += recR2 ;
            }
        } else if  (bbR1 && recR1 > 0) {
            BB *bb = bbNew (h) ;
            bb->gzBuffer = bbR1->gzBuffer ;
            bb->nSeqs = recR1 ;
            channelPut (rd->outChan, &bb, BB*) ;
            records1 += recR1 ;
        }
    }

    channelClose (rd->outChan) ;
    io_uring_queue_exit (&fs1.ring) ;
    if  (rd->pairedEnd) io_uring_queue_exit (&fs2.ring) ;
    close (fd1) ;
    if  (rd->pairedEnd) close (fd2) ;
    ac_free (h) ;
#endif
} /* saUringSequenceParser */

/**************************************************************/
/**************************************************************/
/**************************************************************/
