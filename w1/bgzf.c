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
    return 0 ;

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
  int            fd ;                         /* raw file descriptor        */
  uint64_t       cPos ;                       /* current compressed offset  */
  uint64_t       uPos ;                       /* current uncompressed offset */

  /* decompressed block cache */
  unsigned char  block[BGZF_MAX_BLOCK] ;      /* current decompressed block */
  size_t         blockLen ;                   /* valid bytes in block[]     */
  size_t         blockOff ;                   /* read cursor within block[] */
} ;

/*---------------------------------------------------------------------------
 * bgzfOpen
 *-------------------------------------------------------------------------*/
BGZFFile *bgzfOpen (const char *fNam)
{
  int fd = open (fNam, O_RDONLY) ;
  if (fd < 0)
    return NULL ;

  BGZFFile *bgzf = calloc (1, sizeof *bgzf) ;
  bgzf->fd       = fd ;
  bgzf->cPos     = 0 ;
  bgzf->uPos     = 0 ;
  bgzf->blockLen = 0 ;
  bgzf->blockOff = 0 ;
  return bgzf ;
} /* bgzfOpen */

/*---------------------------------------------------------------------------
 * bgzfClose
 *-------------------------------------------------------------------------*/
void bgzfClose (BGZFFile *bgzf)
{
  if (bgzf)
    {
      close (bgzf->fd) ;
      free (bgzf) ;
    }
} /* bgzfClose */

/*---------------------------------------------------------------------------
 * bgzfSeek — seek to a compressed block boundary
 *-------------------------------------------------------------------------*/
int bgzfSeek (BGZFFile *bgzf, uint64_t cOffset)
{
  if (!bgzf) return -1 ;
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
  if (!bgzf || !buf || len == 0) return 0 ;

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
  return bgzf ? bgzf->uPos : 0 ;
} /* bgzfTell */
