/*  bgzf.c — Blocked GNU Zip Format (BGZF) reader and GZI index
 *
 *  See bgzf.h for a full description of BGZF, GZI, and the rationale
 *  for this self-contained reimplementation.
 *
 *  BGZF block layout (all multi-byte fields are little-endian):
 *
 *    Offset  Size  Field
 *    ------  ----  -----
 *      0      2    gzip magic: 0x1f 0x8b
 *      2      1    CM  = 8 (deflate)
 *      3      1    FLG = 0x04 (FEXTRA set)
 *      4      4    MTIME (ignored)
 *      8      1    XFL  (ignored)
 *      9      1    OS   (ignored)
 *     10      2    XLEN = 6
 *     12      2    extra tag SI1='B' SI2='C'
 *     14      2    extra field length = 2
 *     16      2    BSIZE: total block size - 1  (so block = BSIZE+1 bytes)
 *     18      ?    compressed data (deflate, no header/trailer)
 *    -8       4    CRC32
 *    -4       4    ISIZE: uncompressed size of this block
 *
 *  Dependencies: zlib only (present on all POSIX systems since 1996).
 *  No other external library is used.  See bgzf.h for philosophy.
 *
 *  Author: derived from BGZF spec; integrated into acedb library style.
 *  This file may be copied freely as a self-contained unit.
 */

#include "bgzf.h"
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <zlib.h>               /* inflate, crc32 — stable since 1996      */
static int bgzfFlushBlock (BGZFFile *bgzf) ;

/*---------------------------------------------------------------------------
 * BGZF EOF block — mandatory 28-byte sentinel defined in SAM spec.
 * Any BGZF reader uses this to distinguish clean EOF from truncation.
 *-------------------------------------------------------------------------*/
static const unsigned char bgzfEOFBlock[28] = {
  0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
  0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
  0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00
} ;
/* file open mode — stored in BGZFFile.mode */
#define BGZF_MODE_READ  0
#define BGZF_MODE_WRITE 1

/*---------------------------------------------------------------------------
 * Internal helpers — little-endian reads
 *-------------------------------------------------------------------------*/

static uint16_t leU16 (const unsigned char *p)
{ return (uint16_t)(p[0] | ((uint16_t)p[1] << 8)) ; }

static uint32_t leU32 (const unsigned char *p)
{ return (uint32_t)(p[0] | ((uint32_t)p[1] << 8)
                          | ((uint32_t)p[2] << 16)
                          | ((uint32_t)p[3] << 24)) ; }

static uint64_t leU64 (const unsigned char *p)
{ return (uint64_t)(p[0] | ((uint64_t)p[1] <<  8)
                          | ((uint64_t)p[2] << 16)
                          | ((uint64_t)p[3] << 24)
                          | ((uint64_t)p[4] << 32)
                          | ((uint64_t)p[5] << 40)
                          | ((uint64_t)p[6] << 48)
                          | ((uint64_t)p[7] << 56)) ; }

/*===========================================================================
 * GZI index
 *=========================================================================*/

/*---------------------------------------------------------------------------
 * bgzfGZILoad — load the .gzi sidecar file
 *   The .gzi format is:
 *     uint64_t  nEntries          (little-endian)
 *     nEntries × { uint64_t cOff, uint64_t uOff }  (little-endian pairs)
 *   We prepend the implicit entry (0, 0) so callers never need to special-
 *   case the start of file.
 *-------------------------------------------------------------------------*/
GZIIndex *bgzfGZILoad (const char *fNam)
{
  if (!fNam || !*fNam)
    messcrash ("bgzfGZILoad called with NULL or empty filename") ;
  /* construct sidecar path: fNam + ".gzi" */
  size_t  ln      = strlen (fNam) ;
  char   *gziPath = malloc (ln + 5) ;
  memcpy (gziPath, fNam, ln) ;
  memcpy (gziPath + ln, ".gzi", 5) ;   /* includes NUL */

  FILE *f = fopen (gziPath, "rb") ;
  free (gziPath) ;
  if (!f)
    return NULL ;

  /* read entry count */
  unsigned char hdr[8] ;
  if (fread (hdr, 1, 8, f) != 8) { fclose (f) ; return NULL ; }
  uint64_t nRaw = leU64 (hdr) ;

  /* allocate: nRaw file entries + 1 implicit (0,0) entry */
  GZIIndex *idx = malloc (sizeof *idx) ;
  idx->nEntries = (size_t)nRaw + 1 ;
  idx->entries  = malloc (idx->nEntries * sizeof (GZIEntry)) ;

  /* entry 0 is always the implicit start-of-file block */
  idx->entries[0].cOff = 0 ;
  idx->entries[0].uOff = 0 ;

  /* read remaining entries */
  unsigned char pair[16] ;
  for (size_t i = 1 ; i <= (size_t)nRaw ; i++)
    {
      if (fread (pair, 1, 16, f) != 16)
        {
          /* truncated index — return what we have */
          idx->nEntries = i ;
          break ;
        }
      idx->entries[i].cOff = leU64 (pair)     ;
      idx->entries[i].uOff = leU64 (pair + 8) ;
    }
  fclose (f) ;

  /* uncompressed size: last 4 bytes of the .gz file (gzip ISIZE field)
   * valid for files < 4 GB; for larger files the GZI last entry uOff
   * plus one block size is a safe lower bound — we use it as an
   * estimate only for chunk count computation in the dispatcher        */
  idx->uSize = 0 ;
  FILE *gz = fopen (fNam, "rb") ;
  if (gz)
    {
      if (fseeko (gz, -4, SEEK_END) == 0)
        {
          unsigned char tail[4] ;
          if (fread (tail, 1, 4, gz) == 4)
            idx->uSize = leU32 (tail) ;         /* modulo 2^32 for large files */
        }
      fclose (gz) ;
    }

  /* for files > 4 GB the ISIZE field wraps; use last GZI entry as floor */
  if (idx->nEntries > 1)
    {
      uint64_t lastU = idx->entries[idx->nEntries - 1].uOff ;
      if (idx->uSize < lastU)
        idx->uSize = lastU + BGZF_MAX_BLOCK ;   /* conservative estimate    */
    }

  return idx ;
} /* bgzfGZILoad */

/*---------------------------------------------------------------------------
 * bgzfGZIFree
 *-------------------------------------------------------------------------*/
void bgzfGZIFree (GZIIndex *idx)
{
  if (idx)
    {
      free (idx->entries) ;
      free (idx) ;
    }
} /* bgzfGZIFree */

/*---------------------------------------------------------------------------
 * bgzfGZIFloor — binary search: largest uOff <= uTarget
 * Returns the compressed offset to pass to bgzfSeek.
 *-------------------------------------------------------------------------*/
uint64_t bgzfGZIFloor (const GZIIndex *idx, uint64_t uTarget)
{
  if (!idx || idx->nEntries == 0)
    messcrash ("bgzfGZIFloor called with NULL or empty index") ;

  size_t lo = 0 ;
  size_t hi = idx->nEntries ;       /* half-open [lo, hi) */

  while (lo + 1 < hi)
    {
      size_t mid = (lo + hi) / 2 ;
      if (idx->entries[mid].uOff <= uTarget)
        lo = mid ;
      else
        hi = mid ;
    }
  return idx->entries[lo].cOff ;
} /* bgzfGZIFloor */

/*===========================================================================
 * BGZFFile — block decompressor
 *=========================================================================*/

struct BGZFFile_ {
  void          *magic ;    /* == bgzfDestroy when valid, NULL after destroy */
  AC_HANDLE      h ;        /* internal handle for bgzf's own allocations    */
  int            fd ;                         /* raw file descriptor        */
  int            mode ;                       /* O_RDONLY or O_WRONLY       */
  uint64_t       cPos ;                       /* current compressed offset  */
  uint64_t       uPos ;                       /* current uncompressed offset */

  /* read side — decompressed block cache */
  unsigned char  block[BGZF_MAX_BLOCK] ;
  size_t         blockLen ;
  size_t         blockOff ;

  /* write side — uncompressed accumulation buffer */
  unsigned char  wbuf[BGZF_MAX_BLOCK] ;       /* pending uncompressed data  */
  size_t         wbufLen ;                    /* bytes pending in wbuf      */
} ;

/*---------------------------------------------------------------------------
 * BGZFFile acedb-style lifecycle
 *-------------------------------------------------------------------------*/

static void bgzfDestroy (void *vp)
{
  BGZFFile *bgzf = (BGZFFile *)vp ;
  if (!bgzf)
    return ;
  if (bgzf->magic != bgzfDestroy)
    messcrash ("Bad call to bgzfDestroy") ;
  bgzf->magic = NULL ;              /* guard against double destroy        */

  if (bgzf->mode == BGZF_MODE_WRITE)
    {
      bgzfFlushBlock (bgzf) ;       /* flush final partial block           */
      write (bgzf->fd, bgzfEOFBlock, 28) ;  /* mandatory EOF sentinel      */
    }
  close (bgzf->fd) ;
  bgzf->fd = -1 ;
  ac_free (bgzf->h) ;               /* frees internal allocations          */
} /* bgzfDestroy */

/******************************************************************/
BGZFFile *bgzfOpen (const char *fNam, AC_HANDLE h0)
{
  int fd = open (fNam, O_RDONLY) ;
  if (fd < 0) return NULL ;

  BGZFFile *bgzf  = (BGZFFile *) halloc (sizeof *bgzf, h0) ;
  bgzf->h         = handleCreate () ;   /* internal handle                 */
  bgzf->fd        = fd ;
  bgzf->mode      = BGZF_MODE_READ ;
  bgzf->magic     = bgzfDestroy ;       /* type tag + double-free guard    */
  blockSetFinalise (bgzf, bgzfDestroy) ;
  return bgzf ;
} /* bgzfOpen */

/******************************************************************/
BGZFFile *bgzfOpenWrite (const char *fNam, AC_HANDLE h0)
{
  int fd = open (fNam, O_WRONLY | O_CREAT | O_TRUNC, 0666) ;
  if (fd < 0) return NULL ;

  BGZFFile *bgzf  = (BGZFFile *) halloc (sizeof *bgzf, h0) ;
  bgzf->h         = handleCreate () ;
  bgzf->fd        = fd ;
  bgzf->mode      = BGZF_MODE_WRITE ;
  bgzf->magic     = bgzfDestroy ;
  blockSetFinalise (bgzf, bgzfDestroy) ;
  return bgzf ;
} /* bgzfOpenWrite */

/*---------------------------------------------------------------------------
 * bgzfClose
 *-------------------------------------------------------------------------*/
void bgzfClose (BGZFFile *bgzf)
{
  ac_free (bgzf) ;          /* triggers bgzfDestroy via blockSetFinalise    */
} /* bgzfClose */

void bgzfCloseJunk (BGZFFile *bgzf)
{
  if (!bgzf) return ;

  if (bgzf->mode == BGZF_MODE_WRITE)
    {
      /* flush any remaining data as a final partial block */
      bgzfFlushBlock (bgzf) ;
      /* write mandatory BGZF EOF sentinel */
      write (bgzf->fd, bgzfEOFBlock, 28) ;
    }

  close (bgzf->fd) ;
  free (bgzf) ;
} /* bgzfClose */

/*---------------------------------------------------------------------------
 * bgzfSeek — seek to a compressed block boundary
 *-------------------------------------------------------------------------*/
int bgzfSeek (BGZFFile *bgzf, uint64_t cOffset)
{
  if (!bgzf || bgzf->magic != bgzfDestroy)
    messcrash ("Bad call to bgzfSeek") ;
  if (bgzf->mode != BGZF_MODE_READ)
    messcrash ("bgzfSeek called on write-mode handle: %s\n", "use bgzfWrite") ;
   
  off_t rc = lseek (bgzf->fd, (off_t)cOffset, SEEK_SET) ;
  if (rc < 0) return -1 ;
  bgzf->cPos     = cOffset ;
  bgzf->blockLen = 0 ;      /* invalidate cached block */
  bgzf->blockOff = 0 ;
  return 0 ;
} /* bgzfSeek */

/*---------------------------------------------------------------------------
 * bgzfReadBlock — internal: decompress one BGZF block at current cPos
 * Returns number of uncompressed bytes, 0 at EOF, -1 on error.
 *-------------------------------------------------------------------------*/
static ssize_t bgzfReadBlock (BGZFFile *bgzf)
{
  /* --- read the 18-byte BGZF header ----------------------------------- */
  unsigned char hdr[18] ;
  ssize_t nr = pread (bgzf->fd, hdr, 18, (off_t)bgzf->cPos) ;
  if (nr == 0) return 0 ;           /* EOF */
  if (nr < 18) return -1 ;          /* truncated */

  /* validate gzip magic and BGZF extra field */
  if (hdr[0] != 0x1f || hdr[1] != 0x8b) return -1 ;   /* not gzip        */
  if (!(hdr[3] & 0x04))                  return -1 ;   /* FEXTRA not set  */
  if (hdr[12] != 'B' || hdr[13] != 'C') return -1 ;   /* not BGZF        */

  uint16_t bsize = leU16 (hdr + 16) ;    /* total block size - 1           */
  size_t   blockTotalSize = (size_t)bsize + 1 ;
  size_t   compressedSize = blockTotalSize - 18 - 8 ;  /* minus hdr+trailer */

  /* --- read compressed payload + 8-byte trailer ----------------------- */
  size_t         payloadOff = 18 ;
  unsigned char *payload    = malloc (compressedSize + 8) ;

  nr = pread (bgzf->fd, payload, compressedSize + 8,
              (off_t)(bgzf->cPos + payloadOff)) ;
  if (nr < (ssize_t)(compressedSize + 8))
    { free (payload) ; return -1 ; }

  /* --- read trailer: CRC32 and ISIZE ---------------------------------- */
  unsigned char *trailer = payload + compressedSize ;
  uint32_t crc32Expected = leU32 (trailer) ;
  uint32_t isize         = leU32 (trailer + 4) ;

  if (isize > BGZF_MAX_BLOCK)
    { free (payload) ; return -1 ; }   /* malformed block                  */

  /* --- inflate -------------------------------------------------------- */
  z_stream zs ;
  memset (&zs, 0, sizeof zs) ;
  zs.next_in   = payload ;
  zs.avail_in  = (uInt)compressedSize ;
  zs.next_out  = bgzf->block ;
  zs.avail_out = BGZF_MAX_BLOCK ;

  /* -15: raw deflate (no zlib wrapper), matching BGZF spec */
  if (inflateInit2 (&zs, -15) != Z_OK)
    { free (payload) ; return -1 ; }

  int zrc = inflate (&zs, Z_FINISH) ;
  inflateEnd (&zs) ;
  free (payload) ;

  if (zrc != Z_STREAM_END)              return -1 ;   /* decompression error */
  if (zs.total_out != isize)            return -1 ;   /* size mismatch       */

  /* --- verify CRC32 --------------------------------------------------- */
  uint32_t crc = (uint32_t)crc32 (0L, bgzf->block, isize) ;
  if (crc != crc32Expected)             return -1 ;   /* CRC mismatch        */

  /* --- update state --------------------------------------------------- */
  bgzf->blockLen = isize ;
  bgzf->blockOff = 0 ;
  bgzf->cPos    += blockTotalSize ;     /* advance compressed position      */

  return (ssize_t)isize ;
} /* bgzfReadBlock */

/*---------------------------------------------------------------------------
 * bgzfRead — decompress up to len bytes into buf, crossing block boundaries
 *-------------------------------------------------------------------------*/
ssize_t bgzfRead (BGZFFile *bgzf, void *buf, size_t len)
{
  if (!bgzf || bgzf->magic != bgzfDestroy)
    messcrash ("Bad call to bgzfRead") ;
  if (bgzf->mode != BGZF_MODE_READ)
    messcrash ("bgzfRead called on write-mode handle") ;
  if (!buf)
    messcrash ("bgzfRead called on null buf") ;
  if (len == 0) return 0 ;

  size_t         done = 0 ;
  unsigned char *dst  = (unsigned char *)buf ;

  while (done < len)
    {
      /* if block cache is empty, decompress next block */
      if (bgzf->blockOff >= bgzf->blockLen)
        {
          ssize_t rc = bgzfReadBlock (bgzf) ;
          if (rc == 0) break ;          /* EOF                              */
          if (rc  < 0) return -1 ;      /* error                            */
        }

      /* copy from block cache */
      size_t avail = bgzf->blockLen - bgzf->blockOff ;
      size_t take  = (len - done < avail) ? (len - done) : avail ;
      memcpy (dst + done, bgzf->block + bgzf->blockOff, take) ;
      bgzf->blockOff += take ;
      bgzf->uPos     += take ;
      done           += take ;
    }

  return (ssize_t)done ;
} /* bgzfRead */

/*---------------------------------------------------------------------------
 * bgzfTell — current uncompressed position
 *-------------------------------------------------------------------------*/
uint64_t bgzfTell (const BGZFFile *bgzf)
{
  if (!bgzf || bgzf->magic != bgzfDestroy)
    messcrash ("Bad call to bgzfTell") ;
  return bgzf ? bgzf->uPos : 0 ;
} /* bgzfTell */
/*---------------------------------------------------------------------------
 * Little-endian write helpers — write side only
 *-------------------------------------------------------------------------*/

static void putLeU16 (unsigned char *p, uint16_t v)
{
  p[0] = (unsigned char)(v)      ;
  p[1] = (unsigned char)(v >> 8) ;
} /* putLeU16 */

static void putLeU32 (unsigned char *p, uint32_t v)
{
  p[0] = (unsigned char)(v)       ;
  p[1] = (unsigned char)(v >>  8) ;
  p[2] = (unsigned char)(v >> 16) ;
  p[3] = (unsigned char)(v >> 24) ;
} /* putLeU32 */

/*---------------------------------------------------------------------------
 * bgzfFlushBlock — internal: compress and write one BGZF block.
 * Compresses bgzf->wbuf[0..wbufLen-1] and writes the complete block
 * (header + compressed payload + trailer) to bgzf->fd.
 * Returns 0 on success, -1 on error.
 *-------------------------------------------------------------------------*/
static int bgzfFlushBlock (BGZFFile *bgzf)
{
  if (bgzf->wbufLen == 0)
    return 0 ;                              /* nothing to flush             */

  /* compressed payload buffer — worst case is slightly larger than input  */
  unsigned char cbuf[BGZF_MAX_BLOCK + 64] ;
  uLongf        cLen = sizeof cbuf ;

  /* compress using raw deflate (-15 = no zlib wrapper, matching BGZF)    */
  z_stream zs ;
  memset (&zs, 0, sizeof zs) ;
  zs.next_in   = bgzf->wbuf ;
  zs.avail_in  = (uInt)bgzf->wbufLen ;
  zs.next_out  = cbuf ;
  zs.avail_out = (uInt)sizeof cbuf ;

  if (deflateInit2 (&zs, Z_DEFAULT_COMPRESSION,
                    Z_DEFLATED,
                    -15,                    /* raw deflate, no header       */
                    8,                      /* memory level                 */
                    Z_DEFAULT_STRATEGY) != Z_OK)
    return -1 ;

  if (deflate (&zs, Z_FINISH) != Z_STREAM_END)
    { deflateEnd (&zs) ; return -1 ; }
  deflateEnd (&zs) ;
  cLen = zs.total_out ;

  /* total block size = 18 (header) + cLen + 8 (trailer)                  */
  uint16_t bsize = (uint16_t)(18 + cLen + 8 - 1) ;  /* BSIZE field = size-1 */

  /* --- build 18-byte BGZF header --------------------------------------- */
  unsigned char hdr[18] ;
  hdr[ 0] = 0x1f ;  hdr[ 1] = 0x8b ;   /* gzip magic                     */
  hdr[ 2] = 8 ;                          /* CM = deflate                   */
  hdr[ 3] = 4 ;                          /* FLG = FEXTRA                   */
  hdr[ 4] = hdr[5] = hdr[6] = hdr[7] = 0 ; /* MTIME = 0                   */
  hdr[ 8] = 0 ;                          /* XFL                            */
  hdr[ 9] = 255 ;                        /* OS = unknown                   */
  hdr[10] = 6 ; hdr[11] = 0 ;           /* XLEN = 6                       */
  hdr[12] = 'B' ; hdr[13] = 'C' ;       /* BGZF extra tag                 */
  hdr[14] = 2 ;   hdr[15] = 0 ;         /* extra field length = 2         */
  putLeU16 (hdr + 16, bsize) ;           /* BSIZE                          */

  /* --- build 8-byte trailer -------------------------------------------- */
  unsigned char trailer[8] ;
  uint32_t crc = (uint32_t)crc32 (0L, bgzf->wbuf, (uInt)bgzf->wbufLen) ;
  putLeU32 (trailer,     crc) ;
  putLeU32 (trailer + 4, (uint32_t)bgzf->wbufLen) ;  /* ISIZE             */

  /* --- write header + compressed payload + trailer ---------------------- */
  if (write (bgzf->fd, hdr,     18)   != 18)   return -1 ;
  if (write (bgzf->fd, cbuf,    cLen) != (ssize_t)cLen) return -1 ;
  if (write (bgzf->fd, trailer, 8)    != 8)    return -1 ;

  bgzf->cPos   += 18 + cLen + 8 ;
  bgzf->uPos   += bgzf->wbufLen ;
  bgzf->wbufLen = 0 ;                   /* reset write accumulation buffer */
  return 0 ;
} /* bgzfFlushBlock */

/*---------------------------------------------------------------------------
 * bgzfWrite
 *-------------------------------------------------------------------------*/
ssize_t bgzfWrite (BGZFFile *bgzf, const void *buf, size_t len)
{
  if (!bgzf || bgzf->magic != bgzfDestroy)
    messcrash ("Bad call to bgzfWrite") ;
  if (bgzf->mode != BGZF_MODE_WRITE)
    messcrash ("bgzfWrite called on read-mode handle") ;
  if (!buf) 
    messcrash ("bgzfWrite called on null buf") ;
  if (len == 0) return 0 ;
  
  const unsigned char *src  = (const unsigned char *)buf ;
  size_t               done = 0 ;

  while (done < len)
    {
      /* fill wbuf up to BGZF_MAX_BLOCK */
      size_t avail = BGZF_MAX_BLOCK - bgzf->wbufLen ;
      size_t take  = (len - done < avail) ? (len - done) : avail ;
      memcpy (bgzf->wbuf + bgzf->wbufLen, src + done, take) ;
      bgzf->wbufLen += take ;
      done          += take ;

      /* flush when block is full */
      if (bgzf->wbufLen == BGZF_MAX_BLOCK)
        if (bgzfFlushBlock (bgzf) < 0)
          return -1 ;
    }

  return (ssize_t)done ;
} /* bgzfWrite */

