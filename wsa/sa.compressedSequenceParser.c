/* sa.compressedSequenceParser.c
 *
 * Dispatch and parse compressed or plain FASTA/FASTQ sequence files.
 * Supports single-end, interleaved paired (FASTA2/FASTQ2), and
 * explicit paired (R1+R2) files in plain, BGZF, or gzip formats.
 *
 * Compression auto-detection:
 *   plain file     → backward scan with pread (fastest)
 *   .gz + .gzi     → BGZF backward scan (parallel, full performance)
 *   .gz + pigz     → pigz pipe forward scan (multithreaded decompress)
 *   .gz only       → gzip pipe forward scan (universal fallback)
 *
 * BGZF (Blocked GNU Zip Format) is a gzip extension defined in the
 * SAM/BAM specification.  Each block is independently compressed,
 * enabling random access via a .gzi sidecar index.  We reimplement
 * BGZF/GZI support in bgzf.[ch] to avoid external dependencies.
 * See bgzf.h for full rationale.
 *
 * Paired files:
 *   Explicit pairs (R1+R2): fileName2 != NULL, format FASTA or FASTQ
 *   Interleaved pairs:      fileName2 == NULL, format FASTA2 or FASTQ2
 *   In both cases chunk boundaries are synchronized by read identifier.
 *   The larger file drives the chunking; r1/r2 are swapped transparently.
 *
 * All error conditions call messcrash() which logs to .err and aborts.
 * This code is designed to never produce silent incorrect results.
 */

#include "sa.h"
#include "bgzf.h"

#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define SCAN_BUF  4096   /* boundary search read buffer, stack allocated  */

#ifndef O_NOATIME
#define O_NOATIME 0      /* Linux-specific; zero is safe on other systems  */
#endif

/*===========================================================================
 * Forward declarations
 *=========================================================================*/
static int readScanId      (unsigned char *buf, int *start, int *stop, int *suffix) ;
static int readScanDna     (BB *bb, unsigned char *buf, int *start, int *stop) ;
static int readScanQuality (unsigned char *buf, int *start, int *stop) ;

/*===========================================================================
 * Compression detection helpers
 *=========================================================================*/

/******************************************************************/
static BOOL isGZ (const char *fNam)
{
  int ln = ustrlen (fNam) ;
  return (ln > 3 && ustrcmp (fNam + ln - 3, ".gz") == 0) ? TRUE : FALSE ;
} /* isGZ */

/******************************************************************/
/* Try both standard sidecar names:
 *   foobar.fastq.gz.gzi   (bgzip standard: append ".gzi" to full name)
 * Returns TRUE and fills sidecarPath if found, FALSE otherwise.
 */
static BOOL hasSidecar (const char *fNam, char *sidecarPath, size_t bufLen)
{
  /* try fileName + ".gzi" */
  snprintf (sidecarPath, bufLen, "%s.gzi", fNam) ;
  if (access (sidecarPath, R_OK) == 0)
    return TRUE ;
  sidecarPath[0] = 0 ;
  return FALSE ;
} /* hasSidecar */

/******************************************************************/
/* Paranoia check: verify BGZF magic in file header */
static BOOL isBGZF (const char *fNam)
{
  BOOL ok = FALSE ;
  unsigned char h[18] ;
  int fd = open (fNam, O_RDONLY) ;
  if (fd >= 0)
    {
      if (read (fd, h, 18) == 18
          && h[0]  == 0x1f && h[1] == 0x8b
          && (h[3] &  0x04)
          && h[12] == 'B'  && h[13] == 'C'
          && h[14] == 0x02 && h[15] == 0x00)
        ok = TRUE ;
      close (fd) ;
    }
  return ok ;
} /* isBGZF */

/******************************************************************/
static BOOL hasPIGZ (void)
{
  return (system ("which pigz > /dev/null 2>&1") == 0) ? TRUE : FALSE ;
} /* hasPIGZ */

/******************************************************************/
static BOOL pairedSwitchNeeded (const char *fNam1, const char *fNam2)
{
  struct stat st1, st2 ;
  if (stat (fNam1, &st1) < 0) messcrash ("Cannot stat %s\n", fNam1) ;
  if (stat (fNam2, &st2) < 0) messcrash ("Cannot stat %s\n", fNam2) ;
  return (st2.st_size > st1.st_size) ? TRUE : FALSE ;
} /* pairedSwitchNeeded */

/*===========================================================================
 * Sanitizers — called after file open, on first decompressed bytes.
 * Any violation calls messcrash() with a descriptive message.
 * Performance on the failure path is irrelevant.
 *=========================================================================*/

/******************************************************************/
/* Read one line from buf starting at *pos into line[0..maxLen-1].
 * Advances *pos past the '\n'. Returns line length. */
static int sanitizerGetLine (const unsigned char *buf, int *pos, int bufLen,
                              unsigned char *line, int maxLen)
{
  int i = 0 ;
  while (*pos < bufLen && buf[*pos] != '\n' && i < maxLen - 1)
    line[i++] = buf[(*pos)++] ;
  line[i] = 0 ;
  if (*pos < bufLen && buf[*pos] == '\n') (*pos)++ ;
  return i ;
} /* sanitizerGetLine */

/******************************************************************/
/* Read a small header from a plain file or pipe into a stack buffer.
 * Returns bytes read. */
static int sanitizerReadHead (int fd, FILE *fp, unsigned char *buf, int bufLen)
{
  if (fd >= 0)
    {
      ssize_t nr = pread (fd, buf, bufLen - 1, 0) ;
      if (nr < 0) nr = 0 ;
      buf[nr] = 0 ;
      return (int)nr ;
    }
  else
    {
      size_t nr = fread (buf, 1, bufLen - 1, fp) ;
      buf[nr] = 0 ;
      return (int)nr ;
    }
} /* sanitizerReadHead */

/******************************************************************/
static void sanitizeFasta (int fd, FILE *fp, const char *fNam)
{
  unsigned char buf[1024] ;
  sanitizerReadHead (fd, fp, buf, sizeof buf) ;
  if (!buf[0])
    messcrash ("Empty file: %s\n", fNam) ;

  int pos = 0 ;
  unsigned char line1[512], line2[512] ;
  sanitizerGetLine (buf, &pos, sizeof buf, line1, sizeof line1) ;
  sanitizerGetLine (buf, &pos, sizeof buf, line2, sizeof line2) ;

  if (line1[0] != '>')
    messcrash ("FASTA format error in %s:\n"
               "  line 1 does not start with '>': %s\n", fNam, line1) ;
  if (!line2[0])
    messcrash ("FASTA format error in %s:\n"
               "  sequence line is empty\n", fNam) ;
  /* verify sequence chars */
  for (int i = 0 ; line2[i] ; i++)
    if (!dnaEncodeChar[(unsigned char)line2[i]] && line2[i] != '\r')
      messcrash ("FASTA format error in %s:\n"
                 "  invalid DNA char '%c' (0x%02x) in sequence line: %s\n",
                 fNam, line2[i], (unsigned char)line2[i], line2) ;
} /* sanitizeFasta */

/******************************************************************/
static void sanitizeFastq (int fd, FILE *fp, const char *fNam)
{
  unsigned char buf[2048] ;
  sanitizerReadHead (fd, fp, buf, sizeof buf) ;
  if (!buf[0])
    messcrash ("Empty file: %s\n", fNam) ;

  int pos = 0 ;
  unsigned char line1[512], line2[512], line3[512], line4[512] ;
  sanitizerGetLine (buf, &pos, sizeof buf, line1, sizeof line1) ;
  sanitizerGetLine (buf, &pos, sizeof buf, line2, sizeof line2) ;
  sanitizerGetLine (buf, &pos, sizeof buf, line3, sizeof line3) ;
  sanitizerGetLine (buf, &pos, sizeof buf, line4, sizeof line4) ;

  if (line1[0] != '@')
    messcrash ("FASTQ format error in %s:\n"
               "  line 1 does not start with '@': %s\n", fNam, line1) ;
  if (!line2[0])
    messcrash ("FASTQ format error in %s:\n"
               "  sequence line is empty\n", fNam) ;
  if (line3[0] != '+')
    messcrash ("FASTQ format error in %s:\n"
               "  line 3 does not start with '+': %s\n", fNam, line3) ;
  if (!line4[0])
    messcrash ("FASTQ format error in %s:\n"
               "  quality line is empty\n", fNam) ;
  if (ustrlen (line2) != ustrlen (line4))
    messcrash ("FASTQ format error in %s:\n"
               "  sequence length %d != quality length %d\n",
               fNam, (int)ustrlen (line2), (int)ustrlen (line4)) ;
  for (int i = 0 ; line2[i] ; i++)
    if (!dnaEncodeChar[(unsigned char)line2[i]] && line2[i] != '\r')
      messcrash ("FASTQ format error in %s:\n"
                 "  invalid DNA char '%c' (0x%02x) in sequence: %s\n",
                 fNam, line2[i], (unsigned char)line2[i], line2) ;
} /* sanitizeFastq */

/******************************************************************/
static void sanitizeFasta2 (int fd, FILE *fp, const char *fNam)
{
  unsigned char buf[2048] ;
  sanitizerReadHead (fd, fp, buf, sizeof buf) ;
  if (!buf[0])
    messcrash ("Empty interleaved FASTA file: %s\n", fNam) ;

  int pos = 0 ;
  unsigned char line1[512], line2[512], line3[512], line4[512] ;
  sanitizerGetLine (buf, &pos, sizeof buf, line1, sizeof line1) ;
  sanitizerGetLine (buf, &pos, sizeof buf, line2, sizeof line2) ;
  sanitizerGetLine (buf, &pos, sizeof buf, line3, sizeof line3) ;
  sanitizerGetLine (buf, &pos, sizeof buf, line4, sizeof line4) ;

  if (line1[0] != '>')
    messcrash ("FASTA2 format error in %s:\n"
               "  record 1 does not start with '>': %s\n", fNam, line1) ;
  /* check /1 suffix on record 1 */
  int ln1 = ustrlen (line1) ;
  if (ln1 < 3
      || line1[ln1-1] != '1'
      || (line1[ln1-2] != '/' && line1[ln1-2] != '.' && line1[ln1-2] != ' '))
    messcrash ("FASTA2 format error in %s:\n"
               "  record 1 identifier missing /1 suffix: %s\n", fNam, line1) ;
  if (!line2[0])
    messcrash ("FASTA2 format error in %s:\n"
               "  sequence line 1 is empty\n", fNam) ;
  if (line3[0] != '>')
    messcrash ("FASTA2 format error in %s:\n"
               "  record 2 does not start with '>': %s\n", fNam, line3) ;
  /* check /2 suffix on record 2 */
  int ln3 = ustrlen (line3) ;
  if (ln3 < 3
      || line3[ln3-1] != '2'
      || (line3[ln3-2] != '/' && line3[ln3-2] != '.' && line3[ln3-2] != ' '))
    messcrash ("FASTA2 format error in %s:\n"
               "  record 2 identifier missing /2 suffix: %s\n", fNam, line3) ;
  if (!line4[0])
    messcrash ("FASTA2 format error in %s:\n"
               "  sequence line 2 is empty\n", fNam) ;
} /* sanitizeFasta2 */

/******************************************************************/
static void sanitizeFastq2 (int fd, FILE *fp, const char *fNam)
{
  unsigned char buf[4096] ;
  sanitizerReadHead (fd, fp, buf, sizeof buf) ;
  if (!buf[0])
    messcrash ("Empty interleaved FASTQ file: %s\n", fNam) ;

  int pos = 0 ;
  unsigned char line1[512], line2[512], line3[512], line4[512] ;
  unsigned char line5[512], line6[512], line7[512], line8[512] ;
  sanitizerGetLine (buf, &pos, sizeof buf, line1, sizeof line1) ;
  sanitizerGetLine (buf, &pos, sizeof buf, line2, sizeof line2) ;
  sanitizerGetLine (buf, &pos, sizeof buf, line3, sizeof line3) ;
  sanitizerGetLine (buf, &pos, sizeof buf, line4, sizeof line4) ;
  sanitizerGetLine (buf, &pos, sizeof buf, line5, sizeof line5) ;
  sanitizerGetLine (buf, &pos, sizeof buf, line6, sizeof line6) ;
  sanitizerGetLine (buf, &pos, sizeof buf, line7, sizeof line7) ;
  sanitizerGetLine (buf, &pos, sizeof buf, line8, sizeof line8) ;

  if (line1[0] != '@')
    messcrash ("FASTQ2 format error in %s:\n"
               "  record 1 does not start with '@': %s\n", fNam, line1) ;
  int ln1 = ustrlen (line1) ;
  if (ln1 < 3
      || line1[ln1-1] != '1'
      || (line1[ln1-2] != '/' && line1[ln1-2] != '.' && line1[ln1-2] != ' '))
    messcrash ("FASTQ2 format error in %s:\n"
               "  record 1 identifier missing /1 suffix: %s\n", fNam, line1) ;
  if (line3[0] != '+')
    messcrash ("FASTQ2 format error in %s:\n"
               "  record 1 line 3 does not start with '+': %s\n", fNam, line3) ;
  if (ustrlen (line2) != ustrlen (line4))
    messcrash ("FASTQ2 format error in %s:\n"
               "  record 1 seq length %d != quality length %d\n",
               fNam, (int)ustrlen (line2), (int)ustrlen (line4)) ;
  if (line5[0] != '@')
    messcrash ("FASTQ2 format error in %s:\n"
               "  record 2 does not start with '@': %s\n", fNam, line5) ;
  int ln5 = ustrlen (line5) ;
  if (ln5 < 3
      || line5[ln5-1] != '2'
      || (line5[ln5-2] != '/' && line5[ln5-2] != '.' && line5[ln5-2] != ' '))
    messcrash ("FASTQ2 format error in %s:\n"
               "  record 2 identifier missing /2 suffix: %s\n", fNam, line5) ;
  if (line7[0] != '+')
    messcrash ("FASTQ2 format error in %s:\n"
               "  record 2 line 3 does not start with '+': %s\n", fNam, line7) ;
  if (ustrlen (line6) != ustrlen (line8))
    messcrash ("FASTQ2 format error in %s:\n"
               "  record 2 seq length %d != quality length %d\n",
               fNam, (int)ustrlen (line6), (int)ustrlen (line8)) ;
} /* sanitizeFastq2 */

/******************************************************************/
/* sanitizePaired: verify first identifiers of R1 and R2 match.
 * Called after both files are open, on their first bytes.
 * fd1/fd2 >= 0 for plain files, fp1/fp2 for pipes (pass -1/NULL otherwise).
 */
static void sanitizePaired (int fd1, FILE *fp1, const char *fNam1,
                             int fd2, FILE *fp2, const char *fNam2)
{
  unsigned char buf1[512], buf2[512] ;
  sanitizerReadHead (fd1, fp1, buf1, sizeof buf1) ;
  sanitizerReadHead (fd2, fp2, buf2, sizeof buf2) ;

  if (!buf1[0]) messcrash ("Empty file: %s\n", fNam1) ;
  if (!buf2[0]) messcrash ("Empty file: %s\n", fNam2) ;

  /* extract and normalize first identifier from each file */
  int start1 = 0, stop1 = 0, suffix1 ;
  int start2 = 0, stop2 = 0, suffix2 ;

  if (buf1[0] != '>' && buf1[0] != '@')
    messcrash ("Paired file %s does not start with '>' or '@'\n", fNam1) ;
  if (buf2[0] != '>' && buf2[0] != '@')
    messcrash ("Paired file %s does not start with '>' or '@'\n", fNam2) ;

  readScanId (buf1, &start1, &stop1, &suffix1) ;
  readScanId (buf2, &start2, &stop2, &suffix2) ;

  if (ustrcmp (buf1 + start1, buf2 + start2) != 0)
    messcrash ("Paired files do not match at first record:\n"
               "  %s: '%s'\n"
               "  %s: '%s'\n"
               "Are these truly paired files?\n",
               fNam1, buf1 + start1,
               fNam2, buf2 + start2) ;
} /* sanitizePaired */

/*===========================================================================
 * readScanId / readScanDna / readScanQuality
 * Work on a private BB buffer. Edit in place. Cursor maintained by caller.
 *=========================================================================*/

/******************************************************************/
/* readScanId: called when buf[*start] == '>' or '@'
 * Strips paired suffix (/1 /2  .1 .2  ' 1' ' 2'), writes 0 at id end.
 * On return:
 *   buf[*start] is first char of id  (just after '>'/'@')
 *   buf[*stop]  is 0 just written
 *   *suffix = 0 unpaired, 1 or 2 if paired suffix found
 *   dictAdd (dict, buf + *start, &n) is immediately valid
 *   caller does: start = *stop + 1
 * Returns total chars scanned from entry *start including '\n'.
 */
static int readScanId (unsigned char *buf, int *start, int *stop, int *suffix)
{
  unsigned char *cp = buf + *start ;     /* points to '>' or '@'                  */
  unsigned char *p  = cp + 1 ;           /* skip '>' or '@'                       */
  *suffix  = 0 ;

  /* scan to end of line — do NOT use ustrlen, buf may not be terminated  */
  while (*p && *p != '\n' && *p != '\r')
    p++ ;
  unsigned char *eol = p ;

  /* strip trailing whitespace */
  while (p > cp + 1 && (p[-1] == ' ' || p[-1] == '\t'))
    p-- ;

  /* detect and strip paired suffix: /1 /2  .1 .2  ' 1' ' 2'           */
  if (p - cp >= 3
      && (p[-1] == '1' || p[-1] == '2')
      && (p[-2] == '/' || p[-2] == '.' || p[-2] == ' '))
    {
      *suffix = p[-1] - '0' ;
      p -= 2 ;
    }

  *p     = 0 ;
  *start = (int)(cp + 1 - buf) ;   /* offset of first id char            */
  *stop  = (int)(p      - buf) ;   /* offset of 0                     */

  return (int)(eol - cp) + 1 ;     /* total chars scanned incl. '\n'     */
} /* readScanId */

/******************************************************************/
/* readScanDna: called when buf[*start] is first DNA char.
 * Encodes in place via dnaEncodeChar[], skips '\n' '\r'.
 * Stops on '>' '@' or 0.
 * On return:
 *   buf[*start..*stop-1] is binary encoded sequence
 *   buf[*stop] == 0
 *   virtualArrayCreate (buf, *start, *stop, char) is valid
 *   caller does: start = *stop + 1
 * Returns total raw chars scanned.
 */
static int readScanDna (BB *bb, unsigned char *buf, int *start, int *stop)
{
  unsigned char *cp  = buf + *start ;
  unsigned char *out = cp ;
  int n ;

  unsigned char atgcn[256] ;
  memset (atgcn, 4, sizeof(atgcn)) ;
  atgcn[A_] = 0 ;
  atgcn[T_] = 1 ;
  atgcn[G_] = 2 ;
  atgcn[C_] = 3 ;
  

  while (*cp && *cp != '>' && *cp != '@')
    {
      unsigned char c = *cp++ ;
      if (c == '\n' || c == '\r') continue ;
      *out++ = (unsigned char)dnaEncodeChar[(int)c] ;
    }

  *out  = 0 ;
  *stop = (int)(out - buf) ;
  n = (int)(cp - buf) - *start + (*stop - *start) ;

  if (bb)
    {
      int i, iMax = *stop - *start ;

      for (cp  = buf + *start, i = 0; i < iMax ; cp++, i++)
	bb->runStat.p.NATGC[natgc[(int)*cp]]++ ;

      
      if (buf == bb->r2Buffer)
	{
	  for (cp  = buf + *start, i = 0; i < LETTERMAX && i < iMax ; cp++, i++)	    
	    bb->runStat.p.letterProfile2[5*i + natgc[(int)*cp]]++ ;
	  }
      else
	{
	  for (cp  = buf + *start, i = 0; i < LETTERMAX && i < iMax ; cp++, i++)	    
	    bb->runStat.p.letterProfile1[5*i + natgc[(int)*cp]]++ ;
	  }
    
      
      if (iMax > bb->runStat.p.maxReadLength)
	bb->runStat.p.maxReadLength = iMax ;
      if (! bb->runStat.p.minReadLength || iMax < bb->runStat.p.minReadLength)
	bb->runStat.p.minReadLength = iMax ;
      array (bb->runStat.p.lengthDistribution, iMax, long int)++ ;
    }
  
  return n ;
} /* readScanDna */

/******************************************************************/
/* readScanQuality: called when buf[*start] == '+'
 * Skips '+' line (optional repeated id), copies quality in place,
 * skips '\n' '\r', stops on '@' or 0.
 * On return:
 *   buf[*start..*stop-1] is ASCII Phred+33 quality string
 *   buf[*stop] == 0
 *   virtualArrayCreate (buf, *start, *stop, char) is valid
 *   caller does: start = *stop + 1
 * Returns total raw chars scanned.
 */
static int readScanQuality (unsigned char *buf, int *start, int *stop)
{
  unsigned char *cp = buf + *start ;
  unsigned char *p  = cp ;

  /* skip '+' line */
  while (*p && *p != '\n') p++ ;
  if (*p == '\n') p++ ;

  unsigned char *out = p ;
  *start    = (int)(p - buf) ;

  while (*p && *p != '@')
    {
      unsigned char c = *p++ ;
      if (c == '\n' || c == '\r') continue ;
      *out++ = (char)c ;
    }

  *out  = 0 ;
  *stop = (int)(out - buf) ;
  return (int)(p - buf) - (int)(cp - buf) ;
} /* readScanQuality */

/*===========================================================================
 * BGZF identifier helpers for paired backward scan
 *=========================================================================*/

/******************************************************************/
/* Read and normalize the first identifier from fNam at uncompressed
 * offset uStart. Returns TRUE on success, fills id[0..idLen-1].
 */
static BOOL bgzfReadFirstId (const char *fNam, GZIIndex *idx,
                              uint64_t uStart, char *id, size_t idLen)
{
  uint64_t  cOff = bgzfGZIFloor (idx, uStart) ;
  BGZFFile *bgzf = bgzfOpen (fNam, 0) ;
  if (!bgzf) return FALSE ;
  bgzfSeek (bgzf, cOff) ;

  /* discard bytes before uStart */
  uint64_t uPos = bgzfTell (bgzf) ;
  while (uPos < uStart)
    {
      unsigned char discard[SCAN_BUF] ;
      size_t  want = (size_t)(uStart - uPos) ;
      if (want > SCAN_BUF) want = SCAN_BUF ;
      ssize_t nr = bgzfRead (bgzf, discard, want) ;
      if (nr <= 0) { bgzfClose (bgzf) ; return FALSE ; }
      uPos += (uint64_t)nr ;
    }

  unsigned char scanBuf[SCAN_BUF] ;
  BOOL found = FALSE ;
  while (!found)
    {
      ssize_t nr = bgzfRead (bgzf, scanBuf, SCAN_BUF) ;
      if (nr <= 0) break ;
      unsigned char *hit = memchr (scanBuf, '>', (size_t)nr) ;
      if (!hit)
        hit = memchr (scanBuf, '@', (size_t)nr) ;
      if (hit)
        {
          unsigned char tmp[256] ;
          int  tstart = 0, tstop = 0, tsuffix = 0 ;
          size_t cpLen = (size_t)nr - (size_t)(hit - scanBuf) ;
          if (cpLen >= sizeof tmp) cpLen = sizeof tmp - 1 ;
          memcpy (tmp, hit, cpLen) ;
          tmp[cpLen] = 0 ;
          readScanId (tmp, &tstart, &tstop, &tsuffix) ;
          snprintf (id, idLen, "%s", tmp + tstart) ;
          found = TRUE ;
        }
    }

  bgzfClose (bgzf) ;
  return found ;
} /* bgzfReadFirstId */

/******************************************************************/
/* Search forward in fNam from uStart for a record whose normalized
 * identifier matches matchId. Returns uncompressed offset of matching
 * '>'/'@', or uStart as fallback.
 */
static uint64_t bgzfFindIdentifier (const char *fNam, GZIIndex *idx,
                                     uint64_t uStart, const char *matchId)
{
  uint64_t  cOff = bgzfGZIFloor (idx, uStart) ;
  BGZFFile *bgzf = bgzfOpen (fNam, 0) ;
  if (!bgzf) return uStart ;
  bgzfSeek (bgzf, cOff) ;

  /* discard before uStart */
  uint64_t uPos = bgzfTell (bgzf) ;
  while (uPos < uStart)
    {
      unsigned char discard[SCAN_BUF] ;
      size_t  want = (size_t)(uStart - uPos) ;
      if (want > SCAN_BUF) want = SCAN_BUF ;
      ssize_t nr = bgzfRead (bgzf, discard, want) ;
      if (nr <= 0) { bgzfClose (bgzf) ; return uStart ; }
      uPos += (uint64_t)nr ;
    }

  unsigned char scanBuf[SCAN_BUF] ;
  uint64_t result = uStart ;

  while (TRUE)
    {
      ssize_t nr = bgzfRead (bgzf, scanBuf, SCAN_BUF) ;
      if (nr <= 0) break ;

      unsigned char *p   = scanBuf ;
      unsigned char *end = scanBuf + nr ;

      while (p < end)
        {
          unsigned char *hit = memchr (p, '>', (size_t)(end - p)) ;
          if (!hit) hit = memchr (p, '@', (size_t)(end - p)) ;
          if (!hit) { p = end ; break ; }

          uint64_t hitUPos = uPos + (uint64_t)(hit - scanBuf) ;
          unsigned char     tmp[256] ;
          int      tstart = 0, tstop = 0, tsuffix = 0 ;
          size_t   cpLen  = (size_t)(end - hit) ;
          if (cpLen >= sizeof tmp) cpLen = sizeof tmp - 1 ;
          memcpy (tmp, hit, cpLen) ;
          tmp[cpLen] = 0 ;
          readScanId (tmp, &tstart, &tstop, &tsuffix) ;

          if (ustrcmp (tmp + tstart, matchId) == 0)
            {
              result = hitUPos ;
              bgzfClose (bgzf) ;
              return result ;
            }
          p    = hit + 1 ;
          uPos = hitUPos + 1 ;
        }
      uPos += (uint64_t)(p - scanBuf) ;
    }

  bgzfClose (bgzf) ;
  return result ;
} /* bgzfFindIdentifier */

/*===========================================================================
 * Internal helper: decompress a BGZF range into a halloc'd buffer.
 * Seeks to the BGZF block containing uStart, discards bytes before uStart,
 * reads (uStop - uStart) bytes into b->r1Buffer or b->r2Buffer.
 * target: 1 = r1Buffer, 2 = r2Buffer.
 *=========================================================================*/
static void bgzfFillBuffer (BB *b, int target,
                             const char *fNam, GZIIndex *idx,
                             uint64_t uStart, uint64_t uStop)
{
  size_t   sz   = (size_t)(uStop - uStart) ;
  char    *buf  = halloc (sz + 1, b->h) ;

  BGZFFile *bgzf = bgzfOpen (fNam, 0) ;
  if (!bgzf) messcrash ("Cannot open BGZF file %s\n", fNam) ;

  uint64_t cOff = bgzfGZIFloor (idx, uStart) ;
  bgzfSeek (bgzf, cOff) ;

  uint64_t uPos = bgzfTell (bgzf) ;
  while (uPos < uStart)
    {
      unsigned char discard[SCAN_BUF] ;
      size_t  want = (size_t)(uStart - uPos) ;
      if (want > SCAN_BUF) want = SCAN_BUF ;
      ssize_t nr = bgzfRead (bgzf, discard, want) ;
      if (nr <= 0) break ;
      uPos += (uint64_t)nr ;
    }

  size_t done = 0 ;
  while (done < sz)
    {
      ssize_t nr = bgzfRead (bgzf, buf + done, sz - done) ;
      if (nr <= 0) break ;
      done += (size_t)nr ;
    }
  bgzfClose (bgzf) ;

  buf[done] = 0 ;
  if (target == 1) { b->r1Buffer = (unsigned char *)buf ; b->r1BufferSize = done ; }
  else             { b->r2Buffer = (unsigned char *)buf ; b->r2BufferSize = done ; }
} /* bgzfFillBuffer */

/*===========================================================================
 * Forward scan — pipe input (pigz or gzip), single file
 *=========================================================================*/

/******************************************************************/
static void readsParseForwardScan (const PP *pp, int run, FILE *fp, DnaFormat format)
{
  off_t  BMAX  = (((off_t)pp->BMAX) << 20) ;
  off_t  pos   = 0 ;
  off_t  start = 0 ;

  /* sentinel char for chunk boundary */
  unsigned char sentinel = (format == FASTQ || format == FASTQ2) ? '@' : '>' ;

  size_t         capBuf   = (size_t)BMAX + SCAN_BUF ;
  unsigned char *chunkBuf = malloc (capBuf) ;
  size_t         chunkLen = 0 ;
unsigned char buf[SCAN_BUF] ;
  size_t nr ;

  while ((nr = fread (buf, 1, sizeof buf, fp)) > 0)
    {
      unsigned char *p   = buf ;
      unsigned char *end = buf + nr ;

      while (p < end)
        {
          unsigned char c = *p++ ;
          pos++ ;

          if (c != sentinel)
            {
              /* ordinary byte — append to chunkBuf                        */
              if (chunkLen >= capBuf)
                { capBuf *= 2 ; chunkBuf = realloc (chunkBuf, capBuf) ; }
              chunkBuf[chunkLen++] = c ;
              continue ;
            }

          /* c == sentinel — check interleaved /1 requirement             */
          if (format == FASTA2 || format == FASTQ2)
            {
              /* peek at identifier in remaining buf or a small copy       */
              unsigned char tmp[256] ;
              size_t cpLen = (size_t)(end - (p - 1)) ;
              if (cpLen >= sizeof tmp) cpLen = sizeof tmp - 1 ;
              memcpy (tmp, p - 1, cpLen) ;
              tmp[cpLen] = 0 ;
              int tstart = 0, tstop = 0, tsuffix = 0 ;
              readScanId (tmp, &tstart, &tstop, &tsuffix) ;
              if (tsuffix != 1)
                {
                  /* not /1 — treat sentinel as ordinary byte             */
                  if (chunkLen >= capBuf)
                    { capBuf *= 2 ; chunkBuf = realloc (chunkBuf, capBuf) ; }
                  chunkBuf[chunkLen++] = c ;
                  continue ;
                }
            }

          /* valid chunk boundary at pos-1                                 */
          off_t hitPos = pos - 1 ;

          if (hitPos - start >= BMAX)
            {
              /* dispatch — chunkBuf holds bytes [start .. hitPos-1]       */
              BB b ;
              b.h = ac_new_handle () ;
              b.run = run ;
              b.r1BufferSize = chunkLen ;
              b.r1Buffer = (unsigned char *)halloc (b.r1BufferSize + 1, b.h) ;
              memcpy (b.r1Buffer, chunkBuf, b.r1BufferSize) ;
              b.r1Buffer[b.r1BufferSize] = 0 ;
              b.r2Buffer     = NULL ;
              b.r2BufferSize = 0 ;
              b.rc.format    = format ;
              channelPut (pp->plChan, &b, BB) ;

              chunkLen = 0 ;
              start    = hitPos ;
            }

          /* sentinel always goes into chunkBuf (new or continued chunk)  */
          if (chunkLen >= capBuf)
            { capBuf *= 2 ; chunkBuf = realloc (chunkBuf, capBuf) ; }
          chunkBuf[chunkLen++] = c ;
        }
    }
  /* flush final chunk */
  if (chunkLen > 0)
    {
      BB b ;
      b.h = ac_new_handle () ;
      b.run = run ;      
      b.r1BufferSize = chunkLen ;
      b.r1Buffer = (unsigned char *)halloc (b.r1BufferSize + 1, b.h) ;
      memcpy (b.r1Buffer, chunkBuf, b.r1BufferSize) ;
      b.r1Buffer[b.r1BufferSize] = 0 ;
      b.r2Buffer     = NULL ;
      b.r2BufferSize = 0 ;
      b.rc.format = format ;
      channelPut (pp->plChan, &b, BB) ;
    }

  free (chunkBuf) ;
  return ;
} /* readsParseForwardScan */

/*===========================================================================
 * Backward scan — plain file, single file
 *=========================================================================*/

/******************************************************************/
static void readsParseBackwardScan (const PP *pp, int run, int fd, DnaFormat format)
{
  struct stat st ;
  fstat (fd, &st) ;
  off_t s = st.st_size ;

  off_t    BMAX = (((off_t)pp->BMAX) << 20) ;
  long int n    = (s + BMAX - 1) / BMAX ;
  off_t    stop = s ;
  off_t    start ;

  unsigned char sentinel = (format == FASTQ || format == FASTQ2) ? '@' : '>' ;
  unsigned char buf[SCAN_BUF] ;

  for (long int ii = n - 1 ; ii >= 0 ; ii--)
    {
      if (!ii)
        {
          start = 0 ;
        }
      else
        {
          off_t pos = (off_t)ii * BMAX ;
          start = pos ;

          while (pos < stop)
            {
              ssize_t nr = pread (fd, buf, SCAN_BUF, pos) ;
              if (nr <= 0) break ;
              unsigned char *hit = memchr (buf, sentinel, (size_t)nr) ;

              /* for interleaved: only accept /1 records as boundaries   */
              if (hit && (format == FASTA2 || format == FASTQ2))
                {
                  unsigned char  tmp[256] ;
                  size_t cpLen = (size_t)nr - (size_t)(hit - buf) ;
                  if (cpLen >= sizeof tmp) cpLen = sizeof tmp - 1 ;
                  memcpy (tmp, hit, cpLen) ;
                  tmp[cpLen] = 0 ;
                  int tstart = 0, tstop = 0, tsuffix = 0 ;
                  readScanId (tmp, &tstart, &tstop, &tsuffix) ;
                  if (tsuffix != 1) hit = NULL ;
                }

              if (hit)
                {
                  start = pos + (off_t)(hit - buf) ;
                  break ;
                }
              pos += nr ;
            }
        }

      BB b ;
      b.h = ac_new_handle () ;
      b.run = run ;
      b.r1BufferSize = (size_t)(stop - start) ;
      b.r1Buffer = (unsigned char *)halloc (b.r1BufferSize + 1, b.h) ;
      b.r2Buffer     = NULL ;
      b.r2BufferSize = 0 ;

      size_t done = 0 ;
      while (done < b.r1BufferSize)
        {
          ssize_t nr = pread (fd, b.r1Buffer + done,
                              b.r1BufferSize - done,
                              (off_t)start + (off_t)done) ;
          if (nr <= 0) break ;
          done += (size_t)nr ;
        }
      b.r1BufferSize = done ;
      b.r1Buffer[done] = 0 ;
      b.rc.format = format ;
      channelPut (pp->plChan, &b, BB) ;

      stop = start ;
    }
  return ;
} /* readsParseBackwardScan */

/*===========================================================================
 * Forward scan paired — pipe input for R1 + R2
 *=========================================================================*/

/******************************************************************/
static void readsParseForwardScanPaired (const PP *pp, int run,
                                          FILE *fp1, FILE *fp2,
                                          DnaFormat format,
                                          BOOL switchR1R2)
{
  off_t  BMAX = (((off_t)pp->BMAX) << 20) ;
  unsigned char sentinel = (format == FASTQ) ? '@' : '>' ;

  /* R1 accumulation */
  off_t          pos1      = 0 ;
  off_t          start1    = 0 ;
  size_t         capBuf1   = (size_t)BMAX + SCAN_BUF ;
  unsigned char *chunkBuf1 = malloc (capBuf1) ;
  size_t         chunkLen1 = 0 ;

  /* R2 accumulation */
  size_t         capBuf2   = (size_t)BMAX + SCAN_BUF ;
  unsigned char *chunkBuf2 = malloc (capBuf2) ;
  size_t         chunkLen2 = 0 ;

  unsigned char buf1[SCAN_BUF], buf2[SCAN_BUF] ;
  size_t nr1 ;

  while ((nr1 = fread (buf1, 1, sizeof buf1, fp1)) > 0)
    {
      unsigned char *p   = buf1 ;
      unsigned char *end = buf1 + nr1 ;

      while (p < end)
        {
          unsigned char *hit = memchr (p, sentinel, (size_t)(end - p)) ;
          if (!hit)
            {
              size_t len = (size_t)(end - p) ;
              if (chunkLen1 + len > capBuf1)
                { capBuf1 *= 2 ; chunkBuf1 = realloc (chunkBuf1, capBuf1) ; }
              memcpy (chunkBuf1 + chunkLen1, p, len) ;
              chunkLen1 += len ;
              pos1 += len ;
              p     = end ;
              break ;
            }

          off_t hitPos = pos1 + (off_t)(hit - p) ;

          /* consume one matching record from R2 */
          BOOL r2done = FALSE ;
          while (!r2done)
            {
              size_t nr2 = fread (buf2, 1, sizeof buf2, fp2) ;
              if (!nr2) { r2done = TRUE ; break ; }
              unsigned char *q    = buf2 ;
              unsigned char *end2 = buf2 + nr2 ;
              while (q < end2)
                {
                  unsigned char *hit2 = memchr (q, sentinel, (size_t)(end2 - q)) ;
                  if (!hit2)
                    {
                      size_t len2 = (size_t)(end2 - q) ;
                      if (chunkLen2 + len2 > capBuf2)
                        { capBuf2 *= 2 ; chunkBuf2 = realloc (chunkBuf2, capBuf2) ; }
                      memcpy (chunkBuf2 + chunkLen2, q, len2) ;
                      chunkLen2 += len2 ;
                      q = end2 ;
                    }
                  else
                    {
                      size_t len2 = (size_t)(hit2 - q) ;
                      if (chunkLen2 + len2 > capBuf2)
                        { capBuf2 *= 2 ; chunkBuf2 = realloc (chunkBuf2, capBuf2) ; }
                      memcpy (chunkBuf2 + chunkLen2, q, len2) ;
                      chunkLen2 += len2 ;
                      q = hit2 ;
                      r2done = TRUE ;
                      break ;
                    }
                }
            }

          if (hitPos - start1 >= BMAX)
            {
              BB b ;
              b.h = ac_new_handle () ;
	      b.run = run ;
              b.r1BufferSize = (size_t)(hitPos - start1) ;
              b.r1Buffer = (unsigned char *)halloc (b.r1BufferSize + 1, b.h) ;
              memcpy (b.r1Buffer, chunkBuf1, b.r1BufferSize) ;
              b.r1Buffer[b.r1BufferSize] = 0 ;
              b.r2BufferSize = chunkLen2 ;
              b.r2Buffer = (unsigned char *)halloc (b.r2BufferSize + 1, b.h) ;
              memcpy (b.r2Buffer, chunkBuf2, b.r2BufferSize) ;
              b.r2Buffer[b.r2BufferSize] = 0 ;

              if (switchR1R2)
                {
                  unsigned char *tmp = b.r1Buffer ; size_t tsz = b.r1BufferSize ;
                  b.r1Buffer = b.r2Buffer ; b.r1BufferSize = b.r2BufferSize ;
                  b.r2Buffer = tmp        ; b.r2BufferSize = tsz ;
                }
	      b.rc.format = format ;
              channelPut (pp->plChan, &b, BB) ;

              chunkLen1 = 0 ;
              chunkLen2 = 0 ;
              start1    = hitPos ;
            }

          p    = hit + 1 ;
          pos1 = hitPos + 1 ;
        }
    }

  /* flush final chunk */
  if (chunkLen1 > 0)
    {
      BB b ;
      b.h = ac_new_handle () ;
      b.run = run ;
      b.r1BufferSize = chunkLen1 ;
      b.r1Buffer = (unsigned char *)halloc (b.r1BufferSize + 1, b.h) ;
      memcpy (b.r1Buffer, chunkBuf1, b.r1BufferSize) ;
      b.r1Buffer[b.r1BufferSize] = 0 ;
      b.r2BufferSize = chunkLen2 ;
      b.r2Buffer = (unsigned char *)halloc (b.r2BufferSize + 1, b.h) ;
      memcpy (b.r2Buffer, chunkBuf2, b.r2BufferSize) ;
      b.r2Buffer[b.r2BufferSize] = 0 ;
      if (switchR1R2)
        {
          unsigned char *tmp = b.r1Buffer ; size_t tsz = b.r1BufferSize ;
          b.r1Buffer = b.r2Buffer ; b.r1BufferSize = b.r2BufferSize ;
          b.r2Buffer = tmp        ; b.r2BufferSize = tsz ;
        }
      b.rc.format = format ;
      channelPut (pp->plChan, &b, BB) ;
    }

  free (chunkBuf1) ;
  free (chunkBuf2) ;
  return ;
} /* readsParseForwardScanPaired */

/*===========================================================================
 * Backward scan paired — BGZF, identifier-synchronized
 *=========================================================================*/

/******************************************************************/
static void readsParseBackwardScanPaired (const PP *pp, int run,
                                           const char *fNam1,
                                           const char *fNam2,
                                           DnaFormat format,
                                           BOOL switchR1R2)
{
  int fd1 = open (fNam1, O_RDONLY | O_NOATIME) ;
  int fd2 = open (fNam2, O_RDONLY | O_NOATIME) ;
  if (fd1 < 0) messcrash ("Cannot open %s\n", fNam1) ;
  if (fd2 < 0) messcrash ("Cannot open %s\n", fNam2) ;
  posix_fadvise (fd1, 0, 0, POSIX_FADV_SEQUENTIAL | POSIX_FADV_WILLNEED) ;
  posix_fadvise (fd2, 0, 0, POSIX_FADV_SEQUENTIAL | POSIX_FADV_WILLNEED) ;

  char sidecar1[4096], sidecar2[4096] ;
  if (!hasSidecar (fNam1, sidecar1, sizeof sidecar1))
    messcrash ("Cannot find GZI index for %s\n  Run: bgzip -r %s\n",
               fNam1, fNam1) ;
  if (!hasSidecar (fNam2, sidecar2, sizeof sidecar2))
    messcrash ("Cannot find GZI index for %s\n  Run: bgzip -r %s\n",
               fNam2, fNam2) ;

  GZIIndex *idx1 = bgzfGZILoad (fNam1) ;
  GZIIndex *idx2 = bgzfGZILoad (fNam2) ;
  if (!idx1) messcrash ("Cannot load GZI index for %s\n", fNam1) ;
  if (!idx2) messcrash ("Cannot load GZI index for %s\n", fNam2) ;

  /* --- sanity: first identifiers must match --------------------------- */
  char id1first[256], id2first[256] ;
  if (!bgzfReadFirstId (fNam1, idx1, 0, id1first, sizeof id1first))
    messcrash ("Cannot read first identifier from %s\n", fNam1) ;
  if (!bgzfReadFirstId (fNam2, idx2, 0, id2first, sizeof id2first))
    messcrash ("Cannot read first identifier from %s\n", fNam2) ;
  if (ustrcmp (id1first, id2first) != 0)
    messcrash ("Paired BGZF files do not match at first record:\n"
               "  %s: '%s'\n  %s: '%s'\n",
               fNam1, id1first, fNam2, id2first) ;

  off_t    BMAX   = (((off_t)pp->BMAX) << 20) ;
  uint64_t uSize1 = idx1->uSize ;
  uint64_t uSize2 = idx2->uSize ;
  long int n      = (long int)((uSize1 + (uint64_t)BMAX - 1) / (uint64_t)BMAX) ;
  uint64_t stop1  = uSize1 ;
  uint64_t stop2  = uSize2 ;

  unsigned char sentinel = (format == FASTQ) ? '@' : '>' ;

  for (long int ii = n - 1 ; ii >= 0 ; ii--)
    {
      uint64_t start1, start2 ;

      if (!ii)
        {
          start1 = 0 ;
          start2 = 0 ;
        }
      else
        {
          uint64_t probe = (uint64_t)ii * (uint64_t)BMAX ;
          char     id1[256], id2[256] ;

          /* --- find boundary in R1 ------------------------------------ */
          {
            uint64_t  cOff = bgzfGZIFloor (idx1, probe) ;
            BGZFFile *bgzf = bgzfOpen (fNam1, 0) ;
            if (!bgzf) messcrash ("Cannot open %s\n", fNam1) ;
            bgzfSeek (bgzf, cOff) ;
            unsigned char scanBuf[SCAN_BUF] ;
            uint64_t uPos = bgzfTell (bgzf) ;
            start1 = probe ;
            while (uPos < stop1)
              {
                ssize_t nr = bgzfRead (bgzf, scanBuf, SCAN_BUF) ;
                if (nr <= 0) break ;
                unsigned char *hit = memchr (scanBuf, sentinel, (size_t)nr) ;
                if (hit)
                  {
                    start1 = uPos + (uint64_t)(hit - scanBuf) ;
                    unsigned char  tmp[256] ;
                    size_t cpLen = (size_t)nr - (size_t)(hit - scanBuf) ;
                    if (cpLen >= sizeof tmp) cpLen = sizeof tmp - 1 ;
                    memcpy (tmp, hit, cpLen) ; tmp[cpLen] = 0 ;
                    int ts = 0, te = 0, tx = 0 ;
                    readScanId (tmp, &ts, &te, &tx) ;
                    snprintf (id1, sizeof id1, "%s", tmp + ts) ;
                    break ;
                  }
                uPos += (uint64_t)nr ;
              }
            bgzfClose (bgzf) ;
          }

          /* --- find boundary in R2 at same probe ---------------------- */
          {
            uint64_t  cOff = bgzfGZIFloor (idx2, probe) ;
            BGZFFile *bgzf = bgzfOpen (fNam2, 0) ;
            if (!bgzf) messcrash ("Cannot open %s\n", fNam2) ;
            bgzfSeek (bgzf, cOff) ;
            unsigned char scanBuf[SCAN_BUF] ;
            uint64_t uPos = bgzfTell (bgzf) ;
            start2 = probe ;
            while (uPos < stop2)
              {
                ssize_t nr = bgzfRead (bgzf, scanBuf, SCAN_BUF) ;
                if (nr <= 0) break ;
                unsigned char *hit = memchr (scanBuf, sentinel, (size_t)nr) ;
                if (hit)
                  {
                    start2 = uPos + (uint64_t)(hit - scanBuf) ;
                    unsigned char  tmp[256] ;
                    size_t cpLen = (size_t)nr - (size_t)(hit - scanBuf) ;
                    if (cpLen >= sizeof tmp) cpLen = sizeof tmp - 1 ;
                    memcpy (tmp, hit, cpLen) ; tmp[cpLen] = 0 ;
                    int ts = 0, te = 0, tx = 0 ;
                    readScanId (tmp, &ts, &te, &tx) ;
                    snprintf (id2, sizeof id2, "%s", tmp + ts) ;
                    break ;
                  }
                uPos += (uint64_t)nr ;
              }
            bgzfClose (bgzf) ;
          }

          /* --- synchronize by identifier ------------------------------ */
          int cmp = ustrcmp (id1, id2) ;
          if (cmp < 0)
            start1 = bgzfFindIdentifier (fNam1, idx1, start1, id2) ;
          else if (cmp > 0)
            start2 = bgzfFindIdentifier (fNam2, idx2, start2, id1) ;
          /* cmp == 0: perfect sync                                       */
        }

      /* --- fill BB buffers directly ----------------------------------- */
      BB b ;
      b.h = ac_new_handle () ;
      b.run = run ;
      bgzfFillBuffer (&b, 1, fNam1, idx1, start1, stop1) ;
      bgzfFillBuffer (&b, 2, fNam2, idx2, start2, stop2) ;

      if (switchR1R2)
        {
          unsigned char *tmp = b.r1Buffer ; size_t tsz = b.r1BufferSize ;
          b.r1Buffer = b.r2Buffer ; b.r1BufferSize = b.r2BufferSize ;
          b.r2Buffer = tmp        ; b.r2BufferSize = tsz ;
        }
      b.rc.format = format ;
      channelPut (pp->plChan, &b, BB) ;

      stop1 = start1 ;
      stop2 = start2 ;
    }

  bgzfGZIFree (idx1) ;
  bgzfGZIFree (idx2) ;
  close (fd1) ;
  close (fd2) ;
  return ;
} /* readsParseBackwardScanPaired */

/*===========================================================================
 * Single-file wrappers
 *=========================================================================*/
#ifdef JUNK
/******************************************************************/
static void readsParsePlain (const PP *pp, int run, const char *fNam, DnaFormat format)
{
  int fd = open (fNam, O_RDONLY | O_NOATIME) ;
  if (fd < 0) messcrash ("Cannot open file %s\n", fNam) ;
  switch (format)
  {
  case FASTA:  sanitizeFasta  (fd, NULL, fNam) ; break ;
  case FASTQ:  sanitizeFastq  (fd, NULL, fNam) ; break ;
  case FASTA2: sanitizeFasta2 (fd, NULL, fNam) ; break ;
  case FASTQ2: sanitizeFastq2 (fd, NULL, fNam) ; break ;
  default: break ; 
  }
  posix_fadvise (fd, 0, 0, POSIX_FADV_SEQUENTIAL | POSIX_FADV_WILLNEED) ;
  readsParseBackwardScan (pp, fd, format) ;
  close (fd) ;
} /* readsParsePlain */
#endif
/******************************************************************/
static void readsParseBGZF (const PP *pp, int run, const char *fNam, DnaFormat format)
{
  int fd = open (fNam, O_RDONLY | O_NOATIME) ;
  if (fd < 0) messcrash ("Cannot open BGZF file %s\n", fNam) ;
  posix_fadvise (fd, 0, 0, POSIX_FADV_SEQUENTIAL | POSIX_FADV_WILLNEED) ;

  char sidecar[4096] ;
  if (!hasSidecar (fNam, sidecar, sizeof sidecar))
    messcrash ("Cannot find GZI index for %s\n  Run: bgzip -r %s\n",
               fNam, fNam) ;
  if (!isBGZF (fNam))
    messcrash ("File has .gzi sidecar but is not valid BGZF: %s\n", fNam) ;

  GZIIndex *idx = bgzfGZILoad (fNam) ;
  if (!idx) messcrash ("Cannot load GZI index for %s\n", fNam) ;

  off_t    BMAX  = (((off_t)pp->BMAX) << 20) ;
  uint64_t uSize = idx->uSize ;
  long int n     = (long int)((uSize + (uint64_t)BMAX - 1) / (uint64_t)BMAX) ;
  uint64_t stop  = uSize ;

  unsigned char sentinel = (format == FASTQ || format == FASTQ2) ? '@' : '>' ;

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
          uint64_t cOff  = bgzfGZIFloor (idx, probe) ;
          BGZFFile *bgzf = bgzfOpen (fNam, 0) ;
          if (!bgzf) messcrash ("Cannot open BGZF file %s\n", fNam) ;
          bgzfSeek (bgzf, cOff) ;
          unsigned char scanBuf[SCAN_BUF] ;
          uint64_t uPos = bgzfTell (bgzf) ;
          start = probe ;
          while (uPos < stop)
            {
              ssize_t nr = bgzfRead (bgzf, scanBuf, SCAN_BUF) ;
              if (nr <= 0) break ;
              unsigned char *hit = memchr (scanBuf, sentinel, (size_t)nr) ;

              if (hit && (format == FASTA2 || format == FASTQ2))
                {
                  unsigned char   tmp[256] ;
                  size_t cpLen = (size_t)nr - (size_t)(hit - scanBuf) ;
                  if (cpLen >= sizeof tmp) cpLen = sizeof tmp - 1 ;
                  memcpy (tmp, hit, cpLen) ; tmp[cpLen] = 0 ;
                  int ts = 0, te = 0, tx = 0 ;
                  readScanId (tmp, &ts, &te, &tx) ;
                  if (tx != 1) hit = NULL ;
                }

              if (hit)
                {
                  start = uPos + (uint64_t)(hit - scanBuf) ;
                  break ;
                }
              uPos += (uint64_t)nr ;
            }
          bgzfClose (bgzf) ;
        }

      BB b ;
      b.h = ac_new_handle () ;
      b.run = run ;
      b.r2Buffer     = NULL ;
      b.r2BufferSize = 0 ;
      bgzfFillBuffer (&b, 1, fNam, idx, start, stop) ;
      b.rc.format = format ;
      channelPut (pp->plChan, &b, BB) ;

      stop = start ;
    }

  bgzfGZIFree (idx) ;
  close (fd) ;
  return ;
} /* readsParseBGZF */

/******************************************************************/
#ifdef JUNK
static void readsParseGZ (const PP *pp, int run, const char *fNam, DnaFormat format)
{
  char cmd[4096] ;
  if (hasPIGZ ())
    snprintf (cmd, sizeof cmd, "pigz -dc '%s'", fNam) ;
  else
    snprintf (cmd, sizeof cmd, "gzip -dc '%s'", fNam) ;

  FILE *fp = popen (cmd, "r") ;
  if (!fp) messcrash ("Cannot decompress file %s\n", fNam) ;
  setvbuf (fp, NULL, _IOFBF, 256 * 1024) ;
  readsParseForwardScan (pp, fp, format) ;
  pclose (fp) ;
  return ;
} /* readsParseGZ */
#endif
/*===========================================================================
 * Top-level dispatcher
 *=========================================================================*/
/******************************************************************/
/* sanitizeGZFile: open a separate pipe just for sanitization.
 * The pipe is opened, a few hundred bytes read, then closed.
 * The scan pipe is opened separately and is never touched here.
 * This avoids consuming bytes from the scan pipe before scanning starts.
 */
static void sanitizeGZFile (const char *cmd, DnaFormat format, const char *fNam)
{
  FILE *fp = popen (cmd, "r") ;
  if (!fp) messcrash ("Cannot open %s for format verification\n", fNam) ;
  setvbuf (fp, NULL, _IOFBF, 4096) ;
  switch (format)
    {
    case FASTA:  sanitizeFasta  (-1, fp, fNam) ; break ;
    case FASTQ:  sanitizeFastq  (-1, fp, fNam) ; break ;
    case FASTA2: sanitizeFasta2 (-1, fp, fNam) ; break ;
    case FASTQ2: sanitizeFastq2 (-1, fp, fNam) ; break ;
    default: break ;
    }
  pclose (fp) ;
} /* sanitizeGZFile */

/******************************************************************/
/******************************************************************/
void saCompressedSequenceParser (const PP *pp, RC *rc)
{
  DnaFormat format = rc->format ;
  const char *fNam1 = rc->fileName1 ;
  const char *fNam2 = rc->fileName2 ;
  int run = rc->run ;

  /* --- contradiction check ------------------------------------------- */
  if (fNam2 && (format == FASTA2 || format == FASTQ2))
    messcrash ("Format is interleaved (FASTA2/FASTQ2) but two files given:\n"
               "  %s\n  %s\n"
               "Use a single interleaved file or switch to FASTA/FASTQ.\n",
               fNam1, fNam2) ;

  if (fNam2)
    {
      /* ================================================================
       * EXPLICIT PAIRED MODE  (fileName2 != NULL, format FASTA or FASTQ)
       * ============================================================== */

      /* mixed compression is almost certainly user error — messcrash    */
      if (isGZ (fNam1) != isGZ (fNam2))
        messcrash ("Paired files have different compression formats:\n"
                   "  %s\n  %s\n"
                   "Both files must use the same compression method.\n",
                   fNam1, fNam2) ;

      /* swap so fNam1 is always the larger file                         */
      BOOL switchR1R2 = pairedSwitchNeeded (fNam1, fNam2) ;
      if (switchR1R2)
        { const char *tmp = fNam1 ; fNam1 = fNam2 ; fNam2 = tmp ; }

      if (!isGZ (fNam1))
        {
          /* plain paired — sanitize uses pread so FILE position safe     */
          FILE *fp1 = fopen (fNam1, "r") ;
          FILE *fp2 = fopen (fNam2, "r") ;
          if (!fp1) messcrash ("Cannot open %s\n", fNam1) ;
          if (!fp2) messcrash ("Cannot open %s\n", fNam2) ;
          sanitizePaired (-1, fp1, fNam1, -1, fp2, fNam2) ;
          readsParseForwardScanPaired (pp, run, fp1, fp2, format, switchR1R2) ;
          fclose (fp1) ; fclose (fp2) ;
        }
      else
        {
          char sidecar1[4096], sidecar2[4096] ;
          BOOL bgzf1 = hasSidecar (fNam1, sidecar1, sizeof sidecar1) && isBGZF (fNam1) ;
          BOOL bgzf2 = hasSidecar (fNam2, sidecar2, sizeof sidecar2) && isBGZF (fNam2) ;

          /* mixed BGZF / plain-gz is almost certainly user error        */
          if (bgzf1 != bgzf2)
            messcrash ("Paired files have mixed BGZF/gzip compression:\n"
                       "  %s (%s)\n  %s (%s)\n"
                       "Both files must be bgzipped or both plain gzipped.\n",
                       fNam1, bgzf1 ? "BGZF" : "gzip",
                       fNam2, bgzf2 ? "BGZF" : "gzip") ;

          if (bgzf1)
            {
              /* BGZF paired — full backward scan with identifier sync   */
              readsParseBackwardScanPaired (pp, run, fNam1, fNam2,
                                            format, switchR1R2) ;
            }
          else
            {
              /* plain gzip paired — pigz/gzip pipe forward scan         */
              fprintf (stderr,
                       "Warning: %s is plain gzip — parallel chunking disabled.\n"
                       "For full performance: bgzip your inputs.\n", fNam1) ;
              const char *tool = hasPIGZ () ? "pigz" : "gzip" ;
              char cmd1[4096], cmd2[4096] ;
              snprintf (cmd1, sizeof cmd1, "%s -dc '%s'", tool, fNam1) ;
              snprintf (cmd2, sizeof cmd2, "%s -dc '%s'", tool, fNam2) ;

              /* sanitize via separate pipes — do not touch scan pipes   */
              sanitizeGZFile (cmd1, format, fNam1) ;
              sanitizeGZFile (cmd2, format, fNam2) ;

              /* fresh pipes for the actual scan                         */
              FILE *fp1 = popen (cmd1, "r") ;
              FILE *fp2 = popen (cmd2, "r") ;
              if (!fp1) messcrash ("Cannot decompress %s\n", fNam1) ;
              if (!fp2) messcrash ("Cannot decompress %s\n", fNam2) ;
              setvbuf (fp1, NULL, _IOFBF, 256 * 1024) ;
              setvbuf (fp2, NULL, _IOFBF, 256 * 1024) ;
              readsParseForwardScanPaired (pp, run, fp1, fp2, format, switchR1R2) ;
              pclose (fp1) ; pclose (fp2) ;
            }
        }
    }
  else
    {
      /* ================================================================
       * SINGLE FILE MODE  (fileName2 == NULL)
       * FASTA, FASTQ        : single-end
       * FASTA2, FASTQ2      : interleaved paired, /1 /2 suffixes required
       * ============================================================== */
      if (!isGZ (fNam1))
        {
          /* plain file — sanitize uses pread on fd, FILE position safe  */
          int fd = open (fNam1, O_RDONLY | O_NOATIME) ;
          if (fd < 0) messcrash ("Cannot open %s\n", fNam1) ;
          switch (format)
            {
            case FASTA:  sanitizeFasta  (fd, NULL, fNam1) ; break ;
            case FASTQ:  sanitizeFastq  (fd, NULL, fNam1) ; break ;
            case FASTA2: sanitizeFasta2 (fd, NULL, fNam1) ; break ;
            case FASTQ2: sanitizeFastq2 (fd, NULL, fNam1) ; break ;
            default: break ;
            }
          posix_fadvise (fd, 0, 0, POSIX_FADV_SEQUENTIAL | POSIX_FADV_WILLNEED) ;
          readsParseBackwardScan (pp, run, fd, format) ;
          close (fd) ;
        }
      else
        {
          char sidecar[4096] ;
          BOOL bgzf = hasSidecar (fNam1, sidecar, sizeof sidecar) && isBGZF (fNam1) ;

          if (bgzf)
            {
              /* BGZF — sanitize via a temporary BGZFFile, closed before scan */
              BGZFFile *bgzf_h = bgzfOpen (fNam1, 0) ;
              if (!bgzf_h) messcrash ("Cannot open BGZF file %s\n", fNam1) ;
              char    sbuf[4096] ;
              ssize_t snr = bgzfRead (bgzf_h, (unsigned char *)sbuf,
                                      sizeof sbuf - 1) ;
              if (snr <= 0) messcrash ("Empty or unreadable BGZF file %s\n", fNam1) ;
              sbuf[snr] = '\0' ;
              bgzfClose (bgzf_h) ;       /* fully closed before scan opens  */

              /* sanitize on sbuf directly                                */
              switch (format)
                {
                case FASTA:  sanitizeFasta  (-1, NULL, fNam1) ; break ;
                case FASTQ:  sanitizeFastq  (-1, NULL, fNam1) ; break ;
                case FASTA2: sanitizeFasta2 (-1, NULL, fNam1) ; break ;
                case FASTQ2: sanitizeFastq2 (-1, NULL, fNam1) ; break ;
                default: break ;
                }
              /* readsParseBGZF opens its own fresh BGZFFile handles     */
              readsParseBGZF (pp, run, fNam1, format) ;
            }
          else
            {
              /* plain gzip — sanitize via separate pipe, scan via fresh pipe */
              fprintf (stderr,
                       "Warning: %s is plain gzip — parallel chunking disabled.\n"
                       "For full performance: bgzip your input (bgzip file.fasta)\n"
                       "Falling back to %s streaming...\n",
                       fNam1, hasPIGZ () ? "pigz" : "gzip") ;
              char cmd[4096] ;
              if (hasPIGZ ())
                snprintf (cmd, sizeof cmd, "pigz -dc '%s'", fNam1) ;
              else
                snprintf (cmd, sizeof cmd, "gzip -dc '%s'", fNam1) ;

              /* sanitize via a separate pipe — never touches the scan pipe */
              sanitizeGZFile (cmd, format, fNam1) ;

              /* fresh pipe for the actual scan                           */
              FILE *fp = popen (cmd, "r") ;
              if (!fp) messcrash ("Cannot decompress %s\n", fNam1) ;
              setvbuf (fp, NULL, _IOFBF, 256 * 1024) ;
              readsParseForwardScan (pp, run, fp, format) ;
              pclose (fp) ;
            }
        }
    }
} /* saCompressedSequenceParser */

/**************************************************************/
/* parse the fasta/fastq buffers into an array of DNA and qualities */ 
void saParseR12Buffers (const PP *pp, BB *bb)
{
  DnaFormat format = bb->rc.format ;
  unsigned char *buf1 = bb->r1Buffer ;
  unsigned char *buf2 = bb->r2Buffer ;
  int   start1 = 0, stop1 = 0, start2 = 0, stop2 = 0, suffix1 = 0, suffix2 = 0 ;
  Array dna, dnaR,  qual ;
  unsigned char prefix ;

  if (! bb->r1Buffer) /* no buffer */
    return ;

  bb->errors = arrayHandleCreate (256, int, bb->h) ;
  bb->txt1 = vtxtHandleCreate (bb->h) ;
  bb->txt2 = vtxtHandleCreate (bb->h) ;
  bb->length = 0 ;
  bb->dnas = arrayHandleCreate (bb->nSeqs, BigArray, bb->h) ;
  bb->dnasR = arrayHandleCreate (bb->nSeqs, BigArray, bb->h) ;
  bb->dict = dictHandleCreate (bb->nSeqs, bb->h) ;
  bb->runStat.p.lengthDistribution = arrayHandleCreate (1024, long int, bb->h) ;
  bb->runStat.insertLengthDistribution = arrayHandleCreate (1024, long int, bb->h) ;
  bb->nSeqs = 0 ;
  bb->errDict = dictHandleCreate (100000, bb->h) ;
  bb->cpuStats = arrayHandleCreate (128, CpuSTAT, bb->h) ;
 
      switch (format)
    {
    case FASTA:
    case FASTA2:
      prefix = '>' ;
      break ;
    case FASTQ:
    case FASTQ2:
      prefix = '@' ;
      break ;
    default:
      prefix = '>' ;
      break ;
    }

  while (buf1[start1] == prefix)
    {
      int n = 0, nn1 = 0, nn2 = 0 ;

      switch (format)
	{
	case FASTA:
	case FASTQ:
	  readScanId (buf1, &start1, &stop1, &suffix1) ;
	  dictAdd (bb->dict, (char *)buf1 + start1, &nn1) ;
	  start1 = stop1 + 1 ;

	  if (buf2)
	    {
	      readScanId (buf2, &start2, &stop2, &suffix2) ;
	      dictAdd (bb->dict, (char *)buf2 + start2, &nn2) ;
	      start2 = stop2 + 1 ;
	      
	      if (suffix1 != 1)
		messcrash ("In fasta paired end read 1 identifier should end with 1:  %s\n", dictName (bb->dict, nn1)) ;
	      if (suffix2 != 2)
		messcrash ("In fasta paired end read 2 identifier should end with 2:  %s\n", dictName (bb->dict, nn2)) ;
	      if (nn1 != nn2)
		messcrash ("In fasta paired end read 1 and 2 identifiers do not match: %s <> %s\n", dictName (bb->dict, nn1), dictName (bb->dict, nn2)) ;
	    }

	  readScanDna (bb, buf1, &start1, &stop1) ;
	  n = stop1 - start1 ;

	  dna = arrayHandleCreate (n + 1, unsigned char, bb->h) ;
	  arrayMax (dna) = n ;
	  memcpy (arrp(dna, 0, char), buf1 + start1, n) ;
	  array (bb->dnas, (nn1 << 1), Array) = dna ;
	  if (0)
	    {
	      dnaR = arrayHandleCreate (n + 1, unsigned char, bb->h) ;
	      arrayMax (dnaR) = n ;
	      memcpy (arrp(dnaR, 0, char), buf1 + start1, n) ;
	      reverseComplement (dnaR) ;
	      array (bb->dnasR, (nn1 << 1), Array) = dnaR ;
	    }
	  start1 = stop1 + 1 ;

	  if (format == FASTQ)
	    {
	      readScanQuality (buf1, &start1, &stop1) ;
	      qual = arrayHandleCreate (n + 1, unsigned char, bb->h) ;
	      arrayMax (qual) = n ;
	      memcpy (arrp(qual, 0, char), buf1 + start1, n) ;
	      array (bb->quals, (nn1 << 1), Array) = qual ;
	      start1 = stop1 + 1 ;
	    }

	  bb->nSeqs++ ;
	  bb->length += n ;
	  bb->runStat.p.nReads++ ;
	  bb->runStat.p.nBase1 += n ;
	  
	  if (buf2)
	    {
	      readScanDna (bb, buf2, &start2, &stop2) ;
	      n = stop2 - start2 ;
	      
	      dna = arrayHandleCreate (n + 1, unsigned char, bb->h) ;
	      arrayMax (dna) = n ;
	      memcpy (arrp(dna, 0, char), buf2 + start2, n) ;
	      array (bb->dnas, (nn1 << 1) | 0x1, Array) = dna ;

	      if (0)
		{
		  dnaR = arrayHandleCreate (n + 1, unsigned char, bb->h) ;
		  arrayMax (dnaR) = n ;
		  memcpy (arrp(dnaR, 0, char), buf2 + start2, n) ;
		  reverseComplement (dnaR) ;
		  array (bb->dnasR, (nn1 << 1|0x1), Array) = dnaR ;
		}
	      start2 = stop2 + 1 ;

	      if (format == FASTQ)
		{
		  readScanQuality (buf2, &start2, &stop2) ;
		  qual = arrayHandleCreate (n + 1, unsigned char, bb->h) ;
		  arrayMax (qual) = n ;
		  memcpy (arrp(qual, 0, char), buf2 + start2, n) ;
		  array (bb->quals, (nn1 << 1) | 0x1, Array) = qual ;
		  start2 = stop2 + 1 ;
		}
	  
	      bb->nPairs++ ;
	      bb->nSeqs++ ;
	      bb->length += n ;
	      bb->runStat.p.nPairs++ ;
	      bb->runStat.p.nReads++ ;
	      bb->runStat.p.nBase2 += n ;
	    }
	  
	  break ;


	  
	case FASTA2:
	case FASTQ2:
	  readScanId (buf1, &start1, &stop1, &suffix1) ;
	  dictAdd (bb->dict, (char *)buf1 + start1, &nn1) ;
	  start1 = stop1 + 1 ;
		  
	  readScanDna (bb, buf1, &start1, &stop1) ;
	  n = stop1 - start1 ;

	  dna = arrayHandleCreate (n + 1, unsigned char, bb->h) ;
	  arrayMax (dna) = n ;
	  memcpy (arrp(dna, 0, char), buf1 + start1, n) ;
	  array (bb->dnas, (nn1 << 1), Array) = dna ;
	  if (0)
	    {
	      dnaR = arrayHandleCreate (n + 1, unsigned char, bb->h) ;
	      arrayMax (dnaR) = n ;
	      memcpy (arrp(dnaR, 0, char), buf1 + start1, n) ;
	      reverseComplement (dnaR) ;
	      array (bb->dnasR, (nn1 << 1), Array) = dnaR ;
	    }
	  start1 = stop1 + 1 ;

	  bb->nSeqs++ ;
	  bb->length += n ;
	  bb->runStat.p.nReads++ ;
	  bb->runStat.p.nBase1 += n ;
	  
	  if (format == FASTQ2)
	    {
	      readScanQuality (buf1, &start1, &stop1) ;
	      qual = arrayHandleCreate (n + 1, unsigned char, bb->h) ;
	      arrayMax (qual) = n ;
	      memcpy (arrp(qual, 0, char), buf1 + start1, n) ;
	      array (bb->quals, (nn1 << 1), Array) = qual ;
	      start1 = stop1 + 1 ;
	    }
	  

	  readScanId (buf1, &start1, &stop1, &suffix2) ;
	  dictAdd (bb->dict, (char *)buf1 + start1, &nn2) ;
	  start1 = stop1 + 1 ;
		  
	  if (suffix1 != 1)
	    messcrash ("In fasta interleavedread 1 identifier should end with 1:  %s\n", dictName (bb->dict, nn1)) ;
	  if (suffix2 != 2)
	    messcrash ("In fasta interleaved read 2 identifier should end with 2:  %s\n", dictName (bb->dict, nn2)) ;
	  if (nn1 != nn2)
	    messcrash ("In fasta interleaved read 1 and 2 identifiers do not match: %s <> %s\n", dictName (bb->dict, nn1), dictName (bb->dict, nn2)) ;

	  readScanDna (bb, buf1, &start1, &stop1) ;
	  n = stop1 - start1 ;

	  dna = arrayHandleCreate (n + 1, unsigned char, bb->h) ;
	  arrayMax (dna) = n ;
	  memcpy (arrp(dna, 0, char), buf1 + start1, n) ;
	  array (bb->dnas, (nn1 << 1) | 0x1, Array) = dna ;
	  if (0)
	    {
	      dnaR = arrayHandleCreate (n + 1, unsigned char, bb->h) ;
	      arrayMax (dnaR) = n ;
	      memcpy (arrp(dnaR, 0, char), buf1 + start1, n) ;
	      reverseComplement (dnaR) ;
	      array (bb->dnasR, (nn1 << 1) | 0x1, Array) = dnaR ;
	    }
	  start1 = stop1 + 1 ;

	  bb->nPairs++ ;
	  bb->nSeqs++ ;
	  bb->length += n ;
	  bb->runStat.p.nPairs++ ;
	  bb->runStat.p.nReads++ ;
	  bb->runStat.p.nBase2 += n ;
	      
	  if (format == FASTQ2)
	    {
	      readScanQuality (buf1, &start1, &stop1) ;
	      qual = arrayHandleCreate (n + 1, unsigned char, bb->h) ;
	      arrayMax (qual) = n ;
	      memcpy (arrp(qual, 0, char), buf1 + start1, n) ;
	      array (bb->quals, (nn1 << 1) | 0x1, Array) = qual ;
	      start1 = stop1 + 1 ;
	    }
	  
	  break ;
	default:
	  break ;
	}
    }

  globalDnaCreate (bb) ;
  return ;
} /* saParseR12Buffers */

/******************************************************************/
/******************************************************************/
/******************************************************************/
#ifdef JUNK

Here is a suggested test matrix to be systematic:
Format × compression (12 combinations):
InputFormat  Expected path
  plain.fasta  FASTAbackward scan,
preadplain

  .fastqFASTQbackward scan,
 preadplain.fastaFASTA2 interleavedbackward scan,
 /1 boundaryplain.fastqFASTQ2 interleavedbackward scan,
 /1 boundaryfile.fasta.gz (bgzf)FASTABGZF backward scanfile.fastq.gz (bgzf)

  FASTQBGZF backward scanfile.fasta.gz (plain gz)FASTApigz/gzip forward scanfile.fastq.gz (plain gz)FASTQpigz/gzip forward scanR1.fasta + R2.fastaFASTA pairedforward scan pairedR1.fastq + R2.fastqFASTQ pairedforward scan pairedR1.fasta.gz + R2.fasta.gz (bgzf)FASTA pairedBGZF backward pairedR1.fastq.gz + R2.fastq.gz (plain gz)FASTQ pairedpigz forward paired
Specific things to verify on each:

Correct record count matches grep -c '>' or grep -c '@'
Correct DNA encoding — spot check a known sequence
Correct pairing — R1/R2 identifiers match after normalization
Chunk boundaries land on > or @ — no split records
The sanitizer correctly rejects a deliberately malformed file
Empty file → clean messcrash message

For timing on large files,
 the metrics that matter most at your scale:

Wall time vs CPU time ratio — should be near 1.0 for BGZF (I/O bound),
 higher for plain gz (decompression bottleneck)
MB/s throughput for each path
Memory high-water mark — should stay near BMAX × nAgents

Good luck with the testing — the compiler passed, now biology is the referee.

#endif
