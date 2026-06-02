#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define BMAX        (3 * 1024 * 1024)   /* 3 MB nominal chunk size          */
#define SCAN_BUF    4096                /* bytes to read when hunting '>'   */

/******************************************************************/
/* BGZF is partitioned gzip with block addresses */
static BOOL isGZ (const char *fNam)
{
  BOOL ok = FALSE ;
  int ln = strlen (fNam) ;
  if (ln > 3 && ! strcmp (fNam + ln - 3, ".gz"))
    ok = TRUE ;
  return ok ;
}  /* isGZ */

/******************************************************************/
/* BGZF is partitioned gzip with block addresses */
static BOOL isBGZF (const char *fNam)
{
  BOOL ok = FALSE ;
  unsigned char h[18] ;
  int fd = open (fNam, O_RDONLY) ;
  if (read(fd, h, 18) == 18
      && h[0] == 0x1f && h[1] == 0x8b    /* gzip magic          */
      && h[3]  & 0x04                      /* FEXTRA flag set     */
      && h[12] == 'B' && h[13] == 'C'      /* BGZF extra tag      */
      && h[14] == 0x02 && h[15] == 0x00  /* extra field size=2  */
      )
    ok = TRUE ;
  
  close(fd);
  return ok ;
}  /* isBGZF */

/******************************************************************/
/* BGZF is partitioned gzip with block addresses */
static BOOL hasPIGZ (void)
{
  BOOL ok = (system("which pigz > /dev/null 2>&1") == 0 ? TRUE : FALSE)  ;
  return ok ;
}  /* hasPIGZ */
    
/******************************************************************/
static void readsParseForwardScan (const PP *pp, FILE *fp)
{
  /* forward scan for pipe input (pigz or gzopen)
   * cannot seek, so dispatch on '>' threshold crossings
   * last chunk flushed after fread exhaustion
   */
  off_t  BMAX  = (((off_t)pp->bMax) << 20) ;  /* megabytes */
  off_t  pos   = 0 ;
  off_t  start = 0 ;

  size_t         capBuf   = (size_t)BMAX + SCAN_BUF ;
  unsigned char *chunkBuf = malloc (capBuf) ;

  unsigned char buf[SCAN_BUF] ;
  size_t nr ;

  while ((nr = fread (buf, 1, sizeof buf, fp)) > 0)
    {
      unsigned char *p   = buf ;
      unsigned char *end = buf + nr ;

      while (p < end)
        {
          unsigned char *hit = memchr (p, '>', (size_t)(end - p)) ;
          if (!hit)
            {
              /* append remainder to chunkBuf, no '>' found in this read */
              size_t len = (size_t)(end - p) ;
              if (pos + (off_t)len - start > (off_t)capBuf)
                {
                  capBuf   *= 2 ;
                  chunkBuf  = realloc (chunkBuf, capBuf) ;
                }
              pos += len ;
              p    = end ;
              break ;
            }

          off_t hitPos = pos + (off_t)(hit - p) ;

          if (hitPos - start >= BMAX)
            {
              /* chunk is large enough — dispatch it */
              RB rb = { .fd = -1, .buf = chunkBuf, .start = start, .stop = hitPos } ;
              chanPut (pp->rbChannel, &rb, RB) ;   /* blocking — backpressure */

              /* new chunk starts at hitPos */
              capBuf   = (size_t)BMAX + SCAN_BUF ;
              chunkBuf = malloc (capBuf) ;
              start    = hitPos ;
            }

          p   = hit + 1 ;
          pos = hitPos + 1 ;
        }
    }

  /* flush final chunk */
  if (pos > start)
    {
      RB rb = { .fd = -1, .buf = chunkBuf, .start = start, .stop = pos } ;
      chanPut (pp->rbChannel, &rb, RB) ;
    }
  else
    free (chunkBuf) ;

  return ;
} /* readsParseForwardScan */

/******************************************************************/

static void readsParseBackwardScan (const PP *pp, int fd)
{
  /* recursivelly cut the file in blocks of size BMAX,
   * each block will start on '>' (fasta) or '@' (fastq)
   * for interleaved file, always start on a pair
   */
  struct stat st ;
  fstat (fd, &st);
  off_t s = st.st_size;
  
  off_t start = 0 ;
  off_t stop = s ;

  RB rb ;
  /* --- chunk count ---------------------------------------------------- */
  off_t BMAX = (((off_t)pp->bMax) << 20) ; /* megabytes */
  long int n = (s + BMAX - 1) / BMAX;    /* ceiling division                 */

  /* --- boundary scan buffer (stack is fine at 4 KB) ------------------- */
  /* pread: thread-safe, no lseek, reads at absolute offset */
  unsigned char buf[SCAN_BUF] ;
  for (long int ii = n - 1; ii >= 0; ii--)
    {
      if (! ii)
	start = 0 ;
      else
	{
	  off_t pos = (off_t) ii * BMAX ;
	  start = pos ;          /* fallback if '>' not found            */
	  
	  while (pos < stop)
	    {
	      ssize_t nr = pread (fd, buf, SCAN_BUF, pos) ;
	      if (nr <= 0) break ;
	      
	      unsigned char *hit = memchr (buf, '>', (size_t)nr);
	      if (hit)
		{
		  start = pos + (off_t)(hit - buf);
		  break;
		}
	      pos += nr;
	    }
	}
      rb.fd = fd ;  rb.start = start ; rb.stop = stop ; 
      chanPut(pp->rbChan, &rb, RB) ;   /* blocking — backpressure here         */
      stop = start ;
    }
  return ;
} /* readsParseBackwardScan */

/******************************************************************/
static void readsParsePlain (const PP *pp, const char *fNam)
{
  int fd = open (fNam, O_RDONLY | O_NOATIME) ;
  if (fd < 0)
    messcrash ("Cannot open file %s\n", fNam) ;
  /* just after open => prefetching */
  posix_fadvise (fd, 0, 0, POSIX_FADV_SEQUENTIAL | POSIX_FADV_WILLNEED) ;
  readsParseBackwardScan (pp, fd) ;
  close (fd) ;
} /* readsParsePlain */

/******************************************************************/
/* To construct a BGZF file :
 *   bgzip file.fasta          # replaces gzip, output is valid .gz
 *   bgzip -r file.fasta.gz    # builds the .gzi index alongside
 *
 * BGZF backward scan strategy:
 *   We work in uncompressed coordinate space using the GZI index.
 *   For each chunk boundary ii * BMAX we find the nearest BGZF block
 *   via bgzfGZIFloor, seek there, decompress forward until we hit '>',
 *   then read the full chunk into a malloc'd RB buffer and dispatch.
 *   posix_fadvise is valid here because we own the fd to the .gz file.
 */
static void readsParseBGZF (const PP *pp, const char *fNam)
{
  /* --- open file and load GZI index ----------------------------------- */
  int fd = open (fNam, O_RDONLY | O_NOATIME) ;
  if (fd < 0)
    messcrash ("Cannot open BGZF file %s\n", fNam) ;

  /* prefetching is valid: we own the fd */
  posix_fadvise (fd, 0, 0, POSIX_FADV_SEQUENTIAL | POSIX_FADV_WILLNEED) ;

  GZIIndex *idx = bgzfGZILoad (fNam) ;
  if (!idx)
    messcrash ("Cannot load GZI index for %s\n"
               "  Run: bgzip -r %s\n", fNam, fNam) ;

  /* --- chunk count in uncompressed coordinate space ------------------- */
  off_t    BMAX  = (((off_t)pp->bMax) << 20) ;   /* megabytes            */
  uint64_t uSize = idx->uSize ;
  long int n     = (long int)((uSize + (uint64_t)BMAX - 1) / (uint64_t)BMAX) ;

  uint64_t stop  = uSize ;

  /* --- backward scan -------------------------------------------------- */
  for (long int ii = n - 1 ; ii >= 0 ; ii--)
    {
      uint64_t start ;

      if (!ii)
        {
          start = 0 ;
        }
      else
        {
          uint64_t probe = (uint64_t)ii * (uint64_t)BMAX ;

          /* find compressed offset of the BGZF block containing probe   */
          uint64_t  cOff = bgzfGZIFloor (idx, probe) ;

          /* open a private BGZFFile for this boundary scan
           * each scan is independent: no shared state                   */
          BGZFFile *bgzf = bgzfOpen (fNam) ;
          if (!bgzf)
            messcrash ("Cannot open BGZF file %s\n", fNam) ;
          bgzfSeek (bgzf, cOff) ;

          /* decompress forward until we find '>' at or after probe      */
          unsigned char scanBuf[SCAN_BUF] ;
          uint64_t      uPos = bgzfTell (bgzf) ;
          start = probe ;           /* fallback if '>' not found          */

          while (uPos < stop)
            {
              ssize_t nr = bgzfRead (bgzf, scanBuf, SCAN_BUF) ;
              if (nr <= 0) break ;
              unsigned char *hit = memchr (scanBuf, '>', (size_t)nr) ;
              if (hit)
                {
                  start = uPos + (uint64_t)(hit - scanBuf) ;
                  break ;
                }
              uPos += (uint64_t)nr ;
            }
          bgzfClose (bgzf) ;
        }

      /* --- read full chunk into RB buffer ----------------------------- */
      size_t         chunkLen = (size_t)(stop - start) ;
      unsigned char *chunkBuf = malloc (chunkLen) ;

      BGZFFile *bgzf = bgzfOpen (fNam) ;
      if (!bgzf)
        messcrash ("Cannot open BGZF file %s\n", fNam) ;

      /* seek to the BGZF block containing start, then discard bytes
       * before start so bgzfTell reaches exactly start                  */
      uint64_t cOff = bgzfGZIFloor (idx, start) ;
      bgzfSeek (bgzf, cOff) ;

      uint64_t uPos = bgzfTell (bgzf) ;
      while (uPos < start)
        {
          unsigned char discard[SCAN_BUF] ;
          size_t   want = (size_t)(start - uPos) ;
          if (want > SCAN_BUF) want = SCAN_BUF ;
          ssize_t nr = bgzfRead (bgzf, discard, want) ;
          if (nr <= 0) break ;
          uPos += (uint64_t)nr ;
        }

      /* read chunkLen bytes into chunkBuf                               */
      size_t done = 0 ;
      while (done < chunkLen)
        {
          ssize_t nr = bgzfRead (bgzf, chunkBuf + done, chunkLen - done) ;
          if (nr <= 0) break ;
          done += (size_t)nr ;
        }
      bgzfClose (bgzf) ;

      RB rb ;
      rb.buf    = chunkBuf ;
      rb.bufLen = done ;
      chanPut (pp->rbChan, &rb, RB) ;   /* blocking — backpressure        */

      stop = start ;
    }

  /* --- cleanup -------------------------------------------------------- */
  bgzfGZIFree (idx) ;
  close (fd) ;
  return ;
} /* readsParseBGZF */
/******************************************************************/
static void readsParseGZ (const PP *pp, const char *fNam)
{
  FILE *fp ;
  char cmd[4096] ;
  
  if (hasPIGZ()) /* multithreaded decompression */
    snprintf (cmd, sizeof cmd, "pigz -dc '%s'", fNam) ;
  else
    snprintf (cmd, sizeof cmd, "gzip -dc '%s'", fNam) ;
  fp = popen (cmd, "r") ;
  if (!fp)
    messcrash ("Cannot decompress file %s\n", fNam) ;

  /* prefetching cannot work on gzipped */
  
  setvbuf(fp, NULL, _IOFBF, 256 * 1024);
  readsParseForwardScan (pp, fp) ;
  pclose (fp) ;
  return ;
} /* readsParseGZ */

/******************************************************************/

void saCompressedSequenceParser (const PP *pp, DnaFormat format, const char *fNam)
{
  if (! isGZ (fNam))         /* plain FASTA */
    readsParsePlain (pp, fNam) ;
  else if (isBGZF (fNam))  /* we may need to check dependencies */
    readsParseBGZF (pp, fNam) ;
  else
    {
      fprintf(stderr,
	      "Warning: %s is plain gzip — parallel chunking disabled.\n"
	      "For full performance: bgzip your input (bgzip file.fasta)\n"
	      "Falling back to pigz streaming...\n", fNam
	      ) ;
      readsParseGZ (pp, fNam) ;
    }
} /* saReadsParseDispatch */

/******************************************************************/
/******************************************************************/

