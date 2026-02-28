/*  File: jumpalign.c
 *  Author: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg. 2002
 *-------------------------------------------------------------------
 * This file is part of the ACEVIEW PACKAGE
 *	Jean Thierry-Mieg (NCBI and CNRS France) mieg@ncbi.nlm.nih.gov
 *
 * Description:
 *   Align short sequences (reads/solexa) onto a large target (the genome)
 *
 * As far as I know the algo is completly novel.
 * The idea is to construct all 4^w chains or repeated words of length w.
 *   If the chains of length w are known, we find the chains of length w+1
 * by considering the previous chain and comparing one more letter 
 * on average 4 tests per base.
 *   We initialise the w=0 chains by linking every base to the next one
 * and the last base of each read to the first base of the target.
 *
 * Compile and run as follows
 *
    gcc -O4 -o jali jumpalign.c
    jumpalign  -help

 * Exported functions:
 * HISTORY:
 *   2008_02_19: The code is a spin-off of chromorepeats developped to
 *               cope with the Solexa data
 *   
 * Created: Nov 2002 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */  
#define MAXWORD 50
#define DEBUG 1

/* SPLIT =  16 index : mask 2 letters = 4 bits i.e. 0xf
 * SPLIT = 256 index : mask 4 letters = 8 bits i.e. 0xff
 */
#define SPLIT 16

#include "../wh/ac.h"
#include "../wh/channel.h"
#include "../wh/wego.h"
#include "../wh/bitset.h"
#include <sys/mman.h>

typedef struct jpStruct {  
  AC_HANDLE h ;
  BOOL justStrand, justAntistrand;
  BOOL gzi, gzo ;
  BOOL exportExactHits, exportHits ;
  int importIndex, exportIndex, exportRepeats ;
  int importReads, exportReads ;
  int exportHash ;
  int readLengthMin ;
  int readLengthMax ;
  int wMax, wMax0 ;
  long unsigned int NN, NN0, pNN ;
  int NFILTERS ;
  int max_threads ;
  Array pps ;
  Array chroms, reads ;
  BigArray tDnaArray, pDnaArray, dna2chrom ;
  DICT *targetDict, *readDict ;
  char *dDna, *tDna ;
  long unsigned int *zz, *zz2, *zzInit, *pzz, *pzz2 ;  
  int tDnaFd ;/* file descriptor to mmap zzInit */
  size_t  tDnaSize ;
  long int *nRep, *nComp, *nHooks ;
  BOOL clipPolyA ;
  BitSet flags ;
  CHAN *inChan, *outChan, *ixldChan1, *ixldChan2, *ixldChan3 ;
  const char *readFileNameI ;
  const char *readFileNameX ;
  const char *indexFileNameI ;
  const char *indexFileNameX ;
  const char *outFileName ;
  const char *pFileName ;
  const char *tFileName ;
  ACEOUT ao, aoHits ;
} JP ;

typedef struct ppStruct { 
  AC_HANDLE h ;
  long unsigned int NN, pNN ;
  int nam, iFilter, w, wMax, wMax0 ; 
  long unsigned int iMin, iMax ; /* this thread searches in this zone */
  long unsigned int *zz, *zz2, *pzz, *pzz2 ; 
  char *dDna, *tDna ;
  BigArray tDnaArray, pDnaArray ; 
  BitSet flags ;
  long int *nRep, *nComp, *nHooks ;
  CHAN *zoneChan, *outChan ;
  int NFILTERS ;
} PP ;
 
typedef struct readStruct { 
  int nam ;
  long unsigned int ln, start, end ;
  BOOL isDown ;
} SEQ ;

typedef struct indexLoaderStruct { int fd ; long unsigned int size ; unsigned long int nn ; long unsigned int *zzInit ; int suffix ; const char *fNam ; } IXLD ;

typedef struct zoneStruct { int w, iFilter ; long unsigned int iMin, iMax  ;} ZONE ;

static void usage (char *err) ;
#define ZBAD ((long unsigned int)0xffffffffffffffff) 

/***************************************************************/
/* Order by read length to accelerate 'reportExactReads ' */
static int seqOrder (const void *a, const void *b)
{
  const SEQ *pa = (const SEQ*)a,  *pb = (const SEQ*)b ;
  int n = pa->ln - pb->ln ;
  if (n) return n ;             /* short reads first */
  return pa->nam - pb->nam ;    /* first names first */
} /* seqOrder */

/***************************************************************/
/* calls to  reportExactReads must preferably be done in ever increasing values of w */
static int reportExactReads (JP *jp, int w, int oldW, int *oldiip) 
{
  ACEOUT ao = jp->aoHits ;
  int ii = *oldiip ;
  Array reads = jp->reads ;
  int iiMax = arrayMax (jp->reads) ;
  long unsigned int nn = 0, nn1 = 0 ;
  SEQ *seq, *chrom ;
  DICT *targetDict = jp->targetDict ;
  DICT *readDict = jp->readDict ;
  long unsigned int z1, *zz  = jp->zz ;
  long unsigned int *pzz  = jp->pzz ;
  long int d2cMax = jp->dna2chrom ? bigArrayMax (jp->dna2chrom) : 0 ;
  int chromIndex, *chromIndexp = jp->dna2chrom  ? bigArrayp (jp->dna2chrom, 0, int) : 0 ;
  const char *chromNam ;

  if (! ao)
    ao = jp->aoHits = aceOutCreate (jp->outFileName, ".hits", jp->gzo, jp->h) ;
  if (oldW > w) ii = 0 ; /* restart at top, suboptimal */
  /* start below the previously scanned shorter reads */
  for (seq = arrp (reads, ii, SEQ) ; ii < iiMax ; ii++, seq++)
    {
      if (w < seq->ln)
	break ;
      if (w > seq->ln)
	continue ;
      /* success, we found a read of length w */
      
      z1 = pzz[seq->start] ;
      if (z1 == ZBAD) continue ;
      if (z1) nn1++ ;
      while (z1) /* loop along the successive exact hits */
	{
	  nn++ ;
	  if (z1/1024 < d2cMax)
	    { 
	      chromIndex = chromIndexp[z1/1024] ;
	      chrom = arrp (jp->chroms, chromIndex, SEQ) ;
	      chromNam = dictName (targetDict, chrom->nam) ;
	      if (seq->isDown)
		aceOutf (ao, "%s\t%lu\t%lu", chromNam, z1 - chrom->start + 1 , z1 + w  - chrom->start) ;
	      else /* reversed read */
		aceOutf (ao, "%s\t%lu\t%lu", chromNam, z1 + w  -chrom->start, z1 - chrom->start + 1 ) ;
	      aceOutf (ao, "\t%s\t%d\t%d\t0\n", dictName (readDict, seq->nam), 1, w) ;
	    }
	  else
	    aceOutf (ao, "%s\t%lu\t%lu\t%s\t%d\t%d\t0\n", seq->isDown ? "Down" : "Up", z1  , z1 + w - 1, dictName (readDict, seq->nam), 1, w) ;
	  z1 = zz[z1] ;
	}
    }
  *oldiip = ii ;
  fprintf (stderr, "// === w = %d exported %lu reads with %lu hits\n", w, nn1, nn) ;
  return nn ;
} /* reportExactReads */

/***************************************************************/
/* calls to  reportExactReads must preferably be done in ever increasing values of w */
static void reportHitsReport (const JP *jp, SEQ *seq, long unsigned int z1)
{
  DICT *targetDict = jp->targetDict ;
  DICT *readDict = jp->readDict ;
  long int d2cMax = jp->dna2chrom ? bigArrayMax (jp->dna2chrom) : 0 ;
  int chromIndex, *chromIndexp = jp->dna2chrom  ? bigArrayp (jp->dna2chrom, 0, int) : 0 ;
  int w = jp->wMax ;
  const char *chromNam ;
  SEQ *chrom ;
  ACEOUT ao = jp->aoHits ;
 
  if (z1/1024 < d2cMax)
    { 
      chromIndex = chromIndexp[z1/1024] ;
      chrom = arrp (jp->chroms, chromIndex, SEQ) ;
      chromNam = chrom->nam ? dictName (targetDict, chrom->nam) : "toto" ;
      if (seq->isDown)
	aceOutf (ao, "%s\t%lu\t%lu", chromNam, z1 - chrom->start + 1 , z1 + w  - chrom->start) ;
      else /* reversed read */
	aceOutf (ao, "%s\t%lu\t%lu", chromNam, z1 + w  -chrom->start, z1 - chrom->start + 1 ) ;
      aceOutf (ao, "\t%s\t%d\t%d\t0\n", dictName (readDict, seq->nam), 1, w) ;
    }
  else
    aceOutf (ao, "%s\t%lu\t%lu\t%s\t%d\t%d\t0\n", seq->isDown ? "Down" : "Up", z1  , z1 + w - 1, dictName (readDict, seq->nam), 1, w) ;

  return ;
} /* reportHitsReport */

/***************************************************************/

static long unsigned int reportHitsNavigate (const JP *jp, SEQ *seq, Array z1s)
{
  long unsigned int z1, z2, last = 0, nn = 1 ;
  int i, k  = arrayMax (z1s), w = jp->wMax ;

  while (k > 0)
    {
      arraySort (z1s, unsignedIntReverseOrder) ;
      for (i = k - 1 ; i >= 0 ; i--)
	{
	  z1 = arr (z1s, i, long unsigned int) ;
	  if (z1 == ZBAD) continue ;
	  if (z1)
	    {
	      if (! last || z1 > last + w)
		{ reportHitsReport (jp, seq, z1) ; nn++ ; }
	      last = z1 ;
	      z2 = jp->zz[z1] ;
	      if (z2)
		arr (z1s, i, long unsigned int) = z2 ;
	      else
		k-- ;
	      break ;
	    }
	}      
    }
  return nn ;
} /* reportHitsNavigate  */

/***************************************************************/
/* collect all initial hits in an array, so we can navigate
 * with all pointers in parallel and do less multi search 
*/
static long unsigned int reportHitsCollect (const JP *jp, int ii, Array z1s)
{
  SEQ *seq = arrp (jp->reads, ii, SEQ) ; 
  long unsigned int i, k, z1, *pzz = jp->pzz ;
  long unsigned int iMax = seq->end ;
  
  arrayMax (z1s) = 0 ;
  for (i = seq->start, k = 0 ; i < iMax ; i++)
    {
      z1 = pzz[i] ;
      if (z1 == ZBAD) continue ;
      if (z1) 
	array (z1s, k++, long unsigned int) = z1 ;
    } 
  return k ;
} /* reportHitsCollect */

/***************************************************************/

static int reportHits (const JP *jp, int *nReadsp, int *nHitsp)
{
  AC_HANDLE h = ac_new_handle () ;
  Array z1s = arrayHandleCreate (1000, long unsigned int, h) ; 
  int ii, nn = 0, iMax = arrayMax (jp->reads) ;
  int n, nn2 = 0 ;
  int w = jp->wMax ;
  SEQ *seq ;
  
  for (ii = 0, seq = arrp (jp->reads, 0, SEQ) ; ii < iMax ; ii++, seq++)
    {
      if (w > seq->ln)
	break ;
      n = reportHitsCollect (jp, ii, z1s) ;
      nn += n ;
      if (n)
	nn2 += reportHitsNavigate (jp, seq, z1s) ;
      arrayMax (z1s) = 0 ;
    }
  ac_free (h) ;
  *nReadsp = nn ; *nHitsp = nn2 ;
  return nn ;
} /* reportHits  */

/***************************************************************/
/***************************************************************/
/* 
 * create a hash length jp->exportHash
 * export it in binary format
 */
static int jumpAlignExportHash (const JP *jp)
{
  AC_HANDLE h = ac_new_handle () ;
  Array z1s = arrayHandleCreate (1000, long unsigned int, h) ; 
  int ii, nn = 0, iMax = arrayMax (jp->reads) ;
  int w = jp->wMax ;
  SEQ *seq ;
  
  messcrash ("exportHash is not written") ;
  for (ii = 0, seq = arrp (jp->reads, 0, SEQ) ; ii < iMax ; ii++, seq++)
    {
      if (w > seq->ln)
	break ;
      if (reportHitsCollect (jp, ii, z1s))
	reportHitsNavigate (jp, seq, z1s) ;
    }
  ac_free (h) ;
  return nn ;
} /* jumpAlignExportHash  */

/***************************************************************/
/***************************************************************/
/* 
 * create a hash length jp->exportHash
 * export it in binary format
 */
typedef struct rpStruct { long unsigned int pos, k ; } RP ;
/*************************************************************************************/

static int rpOrder (const void *a, const void *b)
{
  const RP *r1 = (const RP *) a, *r2 = (const RP *) b ;
  return r1->pos < r2->pos ? -1 : 1 ;
}

/*************************************************************************************/

static int jumpAlignExportRepeats (const JP *jp)
{
  RP *rp ;
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (jp->outFileName, ".repeats", jp->gzo, h) ;
  long unsigned int *zp, *zz = jp->zz, nn = 0, pos1, pos2, k ;
  register long unsigned int i, NN = jp->NN - jp->wMax, dw = jp->wMax / 2 ;
  BigArray bb = bigArrayHandleCreate (100000, RP, h) ;
  long long unsigned int ii = 0 ;

  for (i = 0, zp = jp->zz ; i <= NN  ; i++, zp++)  /* scan the genome */
    {
      if (*zp)
	{
	  int k = 1 ;
	  register long unsigned int *zq = zp, *zr ;
	  
	  /* count the repeats, then export them */
	  while (*zq)
	    {
	      k++ ; zq = zz + *zq ;
	    }
	  /* report the numbers */
	  zq = zp ;
	  rp = bigArrayp (bb, ii++, RP) ;
	  rp->pos = i + dw ;
	  rp->k = k ;

	  while (*zq)
	    {
	      rp = bigArrayp (bb, ii++, RP) ;
	      rp->pos = *zq + dw ;
	      rp->k = k ;

	      zr = zq ;
	      zq = zz + *zq ;
	      *zr = 0 ;
	    }
	  nn += k ;
	}
    }
  
  bigArraySort (bb, rpOrder) ;
  ii = bigArrayMax (bb) + 1 ; /* we want to loop once more to export the last segment */
  rp = bigArrp (bb, 0, RP) - 1 ;
  k = 0 ; pos1 = pos2 = 0 ; 
  aceOutDate (ao, "###", "Repeats") ;
  aceOutf (ao, "## Number of exact repeats of length %d  centered at each position\n", jp->wMax) ;
  
  aceOutf (ao, "## Sequence\tStrand\tFrom\tTo\tRepeats\tdebug dx mseq start seaq\n") ;
  while (rp++, ii--)
    {
      if (ii == 0 ||  rp->k != k || rp->pos != pos2 + 1)
	{
	  if (k)
	    {
	      long unsigned int dx = (pos1) >> 10 ;
	      int nTarget = bigArray (jp->dna2chrom, dx, int) ;
	      SEQ *seq = arrp (jp->chroms, nTarget, SEQ) ;

	      if (seq->isDown)
		aceOutf (ao, "%s\t%c\t%lu\t%lu\t%lu\tdx=%lu nTarget=%d start=%lu pos1=%lu pos2=%lu\n"
			 , dictName (jp->targetDict, seq->nam)
			 , seq->isDown ? '+' : '-'
			 , pos1 - seq->start + 1, pos2 - seq->start + 1, k
			 , dx, nTarget, seq->start, pos1, pos2
			 ) ;
	    }
	  pos1 = rp->pos ;
	}
      k = rp->k ;
      pos2 = rp->pos ;
    }

  ac_free (h) ;
  fprintf (stderr, "Found %lu repeats of length %d\n", nn, jp->wMax) ;
  return nn ;
} /* jumpAlignExportRepeats  */

/***************************************************************/
/***************************************************************/

static void ppShow (const JP *jp, char *title)
{ 
  char cc, *tDna ;
  long unsigned int i ;
  long unsigned int NN = jp->NN ;
  long unsigned int *zz  = jp->zz ;

  printf ("ppShow %s\n", title ? title : "") ;
  tDna = jp->tDna ? jp->tDna : (bigArrp (jp->tDnaArray, 0, char)) ;
  for (i = 0 ; i < NN ; i++)
    {
      cc = dnaDecodeChar[(int)tDna[i]] ;
      printf ( "%c\t%lu\t%lu", cc, i, zz[i])  ;
      printf ("\n") ;
    }
}  /* ppShow */

/***************************************************************/
/***************************************************************/

static BOOL fdBigWrite (int fd, void *vp, long int nn, const char *fNam)
{
  long int m,  n = 0 ;
  long unsigned int k, gb ;
  char *cp = (char *) vp ;

  gb = 1 ; gb = gb << 30 ; /* 1 Giga */
  while (n < nn)
    {
      if (nn - n < gb)
	k = nn - n ;
      else
	k = gb ;
      m = write (fd, cp + n, k) ;
      if (m < 0) 
	messcrash ("fdBigWrite failed to write %ld bytes, %ld written, %ld remain, file %s\n"
		   , nn, nn - n, n, fNam) ; 
      n += m ;
    }
      return TRUE ;
} /* fdBigWrite  */

/***************************************************************/

static BOOL fdBigRead (int fd, void *vp, long int nn, const char *fNam)
{
  long int m, n = 0 ;
  long unsigned int k, gb ;
  char *cp = (char *) vp ;

  gb = 1 ; gb = gb << 30 ; /* 1 Giga */
  while (n < nn)
    {
      if (nn - n < gb)
	k = nn - n ;
      else
	k = gb ;
      m = read (fd, cp + n, k) ;
      if (m < 0) 
	messcrash ("fdBigRead failed to write %ld bytes, %ld written, %ld remain, file %s\n"
		   , nn, nn - n, n, fNam) ; 
      n += m ;
    }
  return TRUE ;
} /* fdBigRead  */

static BOOL bigmemset (void *vp, unsigned char cc, long int  nn)
{
  long int n = 0 ;
  long unsigned int k, gb ;
  char *cp = (char *) vp ;
 
  gb = 1 ; gb = gb << 30 ; /* 1 Giga */
  while (n < nn)
    {
      if (nn - n < gb)
	k = nn - n ;
      else
	k = gb ;
      memset (cp + n, cc, k) ;
      n += k ;
    }
  return TRUE ;
} /* bigmemset  */

/***************************************************************/
/***************************************************************/

static void jumpAlignExportIndex (JP *jp, int wMax) 
{
  long unsigned int NN = jp->NN ;
  long unsigned int *zz = jp->zz ;
  int suffix = 0 ;
  BOOL useMemoryMapping = FALSE ;

  if (1) /* export the location of the repeats */
    {
      AC_HANDLE h = ac_new_handle () ;
      ACEOUT ao ;
      DICT *targetDict = jp->targetDict ;
      char *fNam = hprintf (h, "%s.jumpalign.chroms", jp->indexFileNameX) ;
      int i, fd =  open (fNam, O_CREAT | O_RDWR | O_BINARY , 0644 );
      
      fdBigWrite (fd, arrp (jp->chroms, 0, SEQ), arrayMax (jp->chroms) * sizeof (SEQ), fNam) ;
      close (fd) ;
      
      ao =  aceOutCreate (jp->indexFileNameX, ".jumpalign.dict", FALSE, h) ;
      {
	int n1 =  (jp->NN >> 16 ) & 0xffff ;
	int n2 = jp->NN & 0xffff ;
	aceOutf (ao, "%d\t%d\t%d\n", jp->wMax, n1, n2) ;
      }
      for (i = 1 ; i <= dictMax (targetDict) ; i++)
	aceOutf (ao, "%s\n", dictName(targetDict, i)) ;
      ac_free (h) ;   
    }

  if (1) /* export the jump table of the chromosomes */
    {
      char *fNam = hprintf (jp->h, "%s.jumpalign.mask", jp->indexFileNameX) ;
      int fd = open (fNam, O_CREAT | O_RDWR | O_BINARY , 0644 );
      long unsigned int n ;

      n = NN ; n = n * sizeof(long unsigned int) ;
      
      fdBigWrite (fd, jp->zz, n, fNam) ;

      close (fd) ;

      fNam = hprintf (jp->h, "%s.jumpalign.dna", jp->indexFileNameX) ;
      fd =  open (fNam, O_CREAT | O_RDWR | O_BINARY , 0644 );
      n =  bigArrayMax (jp->tDnaArray) ;

      fdBigWrite (fd, bigArrp (jp->tDnaArray, 0, char), n, fNam) ;

      close (fd) ;
    }
  if (jp->zz) free (jp->zz) ; 
  jp->zz = 0 ; /* make space */
  if (jp->zz2) free (jp->zz2) ; 
  jp->zz2 = 0 ; /* make space */

   /* create and export the initialization table */
  for (suffix = 0 ; suffix < SPLIT ; suffix++)
    {   
      int fd ;
      unsigned long int NI, NI1 ;
      long unsigned int missing, ii = 0 ;
      unsigned long int oligo, mask ;
      char *cp ;
      int ok ;
      char *fNam = hprintf (jp->h, "%s.%s%d%s.%d", jp->indexFileNameX, "jumpalign.init_",SPLIT,"_table", suffix) ;
      
      mask = 1 ;  mask <<= 2 * wMax ; 
      NI = mask >>  (SPLIT == 256 ? 8 : 4) ; 
      mask-- ; missing = NI ;
      NI1 = NI * sizeof (long unsigned int) ;
      
      unlink (fNam) ;
      fd = open (fNam, O_CREAT | O_RDWR | O_BINARY , 0644 );
      if (useMemoryMapping)
	zz = mmap(NULL, NI1, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0) ;
      else
	{
	  zz = (long unsigned int *) malloc (NI1) ;
	  bigmemset (zz, 0, NI1) ;
	}
      /* we wish to scan the genome and write the first address */
      for (cp = bigArrp(jp->tDnaArray, 0, char), missing = NI, ii = 0, ok = 0, oligo = 0 ; missing > 0 && ii < NN ; cp++, ii++)
	{
	  /* construire l'oligo, si pas connu on l'ajoute */
	  oligo <<= 2 ; oligo &= mask ;
	  switch ((int) *cp)
	    {
	    case A_: oligo |= 0x0 ; ok++; break ;
	    case T_: oligo |= 0x1 ; ok++; break ;
	    case G_: oligo |= 0x2 ; ok++; break ;
	    case C_: oligo |= 0x3 ; ok++; break ;
	    default: ok = 0 ; break ;
	    }
	  
	  /* use the last 8 bits as a suffix flag and create 16 smaller tables */
	  if (
	      ( SPLIT ==  16 && (oligo &  0xf) == suffix) ||
	      ( SPLIT == 256 && (oligo & 0xff) == suffix) 
	      )
	    {
	      unsigned long int oligo2 = oligo >> (SPLIT == 256 ? 8 : 4) ;
	      if (ok >= wMax && ! zz[oligo2]) /* the oligo exists */
		{
		  missing-- ;
		  zz[oligo2] = ii - wMax + 1 ; /* jp->NN0  avoids zero */
		  if (0 && ii-wMax+1 < 12) fprintf (stderr,"%lu %lu\n", ii - wMax + 1, oligo) ;
		}
	    }   
	  
	}
      if (useMemoryMapping)
	{
	  munmap (zz, NI1) ;
	  zz = 0 ; /* zz is unaccessible, but shuld not be freed */
	}
      else
	{
	  fdBigWrite (fd, zz, NI1, fNam) ; 
	  free (zz) ; zz = 0 ; /* make space */
	}
      close (fd) ;
    }
} /* jumpAlignExportIndex */

/***************************************************************/
/***************************************************************/
/* read the initialisation index from disk
 ^ scan the oligos to find where they jump
 * destroy the mask
 * read the chromosome jump table and paste it in zz
 */  
static void jumpAlignImportIndexDict  (JP *jp)
{
  /* read the initialisation index */
  if (1) /* parse the dict */
    {
      AC_HANDLE h = ac_new_handle () ;
      DICT *dict ;
      ACEIN ai = aceInCreate (hprintf (h, "%s.%s", jp->indexFileNameI, "jumpalign.dict"), FALSE, h) ;
      if (! ai)
	messcrash ("cannot open the jumpalign.dict") ;

      ac_free (jp->targetDict) ;
      dict = jp->targetDict = dictHandleCreate (256, jp->h) ;
      aceInCard (ai) ;
      {
	int n1, n2 ; 
	aceInInt (ai, &n1) ; jp->wMax0 = n1 ;
	aceInStep (ai, '\t') ;	aceInInt (ai, &n1) ; 
	aceInStep (ai, '\t') ; aceInInt (ai, &n2) ; 
	jp->NN = n1 ; jp->NN <<= 16 ; jp->NN += n2 ; 
      }
      while (aceInCard (ai))
	dictAdd (dict, aceInWord (ai), 0) ;
      ac_free (h) ;
    }
} /* jumpAlignImportIndexDict */

/***************************************************************/

static void jumpAlignImportIndexChroms  (JP *jp)
{
  if (1) /* parse the chroms coordinates */
    {
      long unsigned int mx, n, fd = open (messprintf ("%s.jumpalign.chroms", jp->indexFileNameI), O_RDONLY | O_BINARY , 0644 );
      if (! fd)
	messcrash ("cannot open the jumpalign.chroms") ;

      ac_free (jp->chroms) ;
      jp->chroms = arrayHandleCreate (1000, SEQ, jp->h) ;
      n = 1 ; mx = 0 ;
      while (n > 0)
	{
	  arrayp (jp->chroms, mx + 1000, SEQ)->ln = 0 ;
	  n = read (fd, arrp (jp->chroms, mx, SEQ), 1000 * sizeof (SEQ)) ;
	  mx += n /  sizeof (SEQ) ;
	  arrayMax (jp->chroms) = mx ;
	}
      close (fd) ;
    }
  if (jp->NN) /* hook the names */
    {
      long int k, kMax = jp->NN/1024 ;
      int i, *ip, iMax =  arrayMax (jp->chroms) ;
      SEQ *chrom = arrp (jp->chroms,0, SEQ) ;
      BigArray d2c = jp->dna2chrom = bigArrayHandleCreate (kMax, int, jp->h) ;
      bigArray (d2c, kMax - 1, int) = 0 ;
      for (i = 0 ; i < iMax ; i++, chrom++)  /* declare the starts */
	bigArray (d2c, chrom->start/1024, int) = i ;
      /* complete the table */
      for (i = 0, k = 0 , ip = 	bigArrayp (d2c, 0, int); k < kMax ; k++)
	if (*ip) 
	  i = *ip ;
	else
	  *ip = i ;
    }      
  return ;
} /* jumpAlignImportIndexDictChroms */

/***************************************************************/

/***************************************************************/

static void jumpAlignImportIndexDna  (JP *jp)
{
  if (1) /* parse the dna */
    {
      long unsigned int NN = jp->NN ;
      char *fNam = hprintf (jp->h, "%s.jumpalign.dna", jp->indexFileNameI) ;
      int fd = open (fNam, O_RDONLY | O_BINARY , 0644 );
      if (! fd)
	messcrash ("cannot open the jumpalign.dna") ;

      ac_free (jp->tDnaArray) ;

      if (1)
	{
	  jp->tDna = mmap(NULL, NN+1, PROT_READ, MAP_PRIVATE, fd, 0) ;
	  jp->tDnaFd = fd ; jp->tDnaSize = NN + 1 ;
	}
      else
	{
	  jp->tDnaArray = bigArrayHandleCreate (NN + 1, long unsigned int, jp->h) ;
	  bigArray (jp->tDnaArray, NN, char) = 0 ; /* double zero terminate */
	  bigArrayMax (jp->tDnaArray) = NN ;
	  
	  fdBigRead (fd, bigArrp (jp->tDnaArray, 0, long unsigned int), NN, fNam) ;
	  
	  close (fd) ;
	}
    }
  
} /* jumpAlignImportIndexDna */

/***************************************************************/
/* parse the jump table */
static void jumpAlignImportIndexMask  (JP *jp)
{
  long unsigned int NN = jp->NN ;
  long int NN1 ;
  char *fNam = hprintf (jp->h, "%s.jumpalign.mask", jp->indexFileNameI) ;
  int fd = open (fNam, O_RDONLY | O_BINARY , 0644 );
  if (! fd)
    messcrash ("cannot open the jumpalign.mask") ;
  
  ac_free (jp->zz) ;
  NN1 = NN ; NN1 = NN1 * sizeof ( long unsigned int) ;
  jp->zz = (long unsigned int *) malloc (NN1) ;
  
  fdBigRead (fd, jp->zz, NN1, fNam) ;
  
  if (jp->NFILTERS > 1)
    {
      jp->zz2 = (long unsigned int *) malloc (NN1) ;
      memcpy (jp->zz2, jp->zz, NN) ;
    }
  close (fd) ;

  return ;
} /* jumpAlignImportIndexMask */

/***************************************************************/

static void jumpAlignImportIndex  (JP *jp)
{
  BOOL debug = FALSE ;

  jumpAlignImportIndexDict (jp) ;
  jumpAlignImportIndexChroms  (jp) ;
  jumpAlignImportIndexDna (jp) ;
  jumpAlignImportIndexMask (jp) ;

  fprintf (stderr, "// %s: Found from the index %lu %d-mers oligos known in the genome\n", timeShowNow (), jp->NN, jp->wMax0) ;
  if (debug) ppShow (jp, "jumpAlignReadInitFromIndex") ;
  return ;
} /* jumpAlignImportIndex */

/***************************************************************/
/***************************************************************/
/***************************************************************/
/* read the Reads to be aligned in binary format
 */  
static void jumpAlignImportReadsDict  (const void *vp)
{
  const JP *jp = vp ;
  int k ;
  /* read the initialisation index */
  if (1) /* parse the dict */
    {
      AC_HANDLE h = ac_new_handle () ;
      DICT *dict = jp->readDict ;
      ACEIN ai = aceInCreate (hprintf (h, "%s.%s", jp->readFileNameI, "jumpalign.reads_dict"), FALSE, h) ;
      if (! ai)
	messcrash ("cannot open the jumpalign.reads_dict") ;

       aceInCard (ai) ;
      {
	int n1, n2, ns ; 

	aceInStep (ai, '\t') ;	aceInInt (ai, &n1) ; 
	aceInStep (ai, '\t') ; aceInInt (ai, &n2) ; 
	aceInStep (ai, '\t') ; aceInInt (ai, &ns) ; 

	channelPut (jp->outChan, &n1, int) ;
	channelPut (jp->outChan, &n2, int) ;
	channelPut (jp->outChan, &ns, int) ;
      }
      while (aceInCard (ai))
	dictAdd (dict, aceInWord (ai), 0) ;
      ac_free (h) ;
    }  
  k = dictMax (jp->readDict) ;
  channelPut (jp->outChan, &k, int) ; /* signal completion */
} /* jumpAlignImportReadsDict */

/***************************************************************/
/* parse the reads coordinates */
static void jumpAlignImportReadsArray  (JP *jp)
{
  char *fNam = hprintf (jp->h, "%s.jumpalign.reads_array", jp->readFileNameI) ;
  long unsigned int mx, fd = open (fNam, O_RDONLY | O_BINARY , 0644 );
  if (! fd)
    messcrash ("cannot open the jumpalign.reads_array") ;
  
  mx = arrayMax (jp->reads) ; mx *= sizeof (SEQ) ;
  fdBigRead (fd, arrp (jp->reads, 0, SEQ), mx, fNam) ;

  if (1)
    {
      SEQ *seq = arrp (jp->reads, 0, SEQ) ;
      int i ;

      mx =   arrayMax (jp->reads) ;
      for (i = 0 ; i < mx ; i++, seq++)
	{
	  if (seq->ln <  jp->readLengthMin)
	    jp->readLengthMin = seq->ln ;
	  if (seq->ln > jp->readLengthMax)
	    jp->readLengthMax = seq->ln ;
	}
    }
  close (fd) ;
} /* jumpAlignImportReadsArray */

/***************************************************************/
/* parse the dna in binary format */
static void jumpAlignImportReadsDna  (JP *jp)
{
  long unsigned int pNN = jp->pNN ;
  char *fNam = hprintf (jp->h, "%s.jumpalign.reads_dna", jp->readFileNameI) ;
  int fd = open (fNam, O_RDONLY | O_BINARY , 0644 );
  if (! fd)
    messcrash ("cannot open the jumpalign.reads.dna") ;
  
  ac_free (jp->pDnaArray) ;
  jp->pDnaArray = bigArrayHandleCreate (pNN + 1, long unsigned int, jp->h) ;
  bigArray (jp->pDnaArray, pNN, char) = 0 ; /* double zero terminate */
  bigArrayMax (jp->pDnaArray) = pNN ;
  
  fdBigRead (fd, bigArrp (jp->pDnaArray, 0, long unsigned int), pNN, fNam) ;
  
  close (fd) ;
} /* jumpAlignImportIndexDna */

/***************************************************************/

static void jumpAlignImportReads  (JP *jp)
{ 
  int k = 1 ;

  /* parallelise the parsing of the readDict as an exercise in multithreading */
  jp->inChan =  channelCreate (2, int, jp->h) ;
  jp->outChan =  channelCreate (2, int, jp->h) ;
  ac_free (jp->readDict) ;
  jp->readDict = dictHandleCreate (256, jp->h) ; /* establish dict before calling wego_go */
 
  wego_go (jumpAlignImportReadsDict, jp, JP) ;
  channelPut (jp->inChan, &k, int) ; /* start jumpAlignImportReadsDict */
  
  if (1)
    {
      int n1 = channelGive (jp->outChan, &n1, int) ;  /* wait on setting jp->pNN */
      int n2 = channelGive (jp->outChan, &n2, int) ;  /* wait on setting jp->pNN */
      int ns = channelGive (jp->outChan, &ns, int) ;  /* wait on setting jp->pNN */
      
      jp->pNN = n1 ; jp->pNN <<= 16 ; jp->pNN += n2 ; 
      ac_free (jp->reads) ;
      jp->reads = arrayHandleCreate (ns, SEQ, jp->h) ;
      array (jp->reads, ns > 0 ? ns - 1 : ns, SEQ).ln = 0 ;
    }
  
  jumpAlignImportReadsArray  (jp) ;
  jumpAlignImportReadsDna (jp) ;
  
  fprintf (stderr, "// %s: Found from the index %d reads contatining %ld base pairs\n"
	   , timeShowNow ()
	   , dictMax (jp->readDict), bigArrayMax (jp->pDnaArray)
	   ) ;
  
  channelGet (jp->outChan, &k, int) ;  /* wait on readDict initialisation */
  
  return ;
} /* jumpAlignImportReads */

/***************************************************************/
/***************************************************************/
/***************************************************************/
/* Export the Reads to be aligned in binary format
 */  
static void jumpAlignExportReadsDict  (JP *jp)
{
  /* read the initialisation index */
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (jp->readFileNameX, ".jumpalign.reads_dict", FALSE, h) ;
  if (! ao)
    messcrash ("cannot create the jumpalign.reads_dict") ;
  
  if (1)
    {
      int n1 =  (jp->pNN >> 16 ) & 0xffff ;
      int n2 = jp->pNN & 0xffff ;
      aceOutf (ao, "%d\t%d\t%d\n"
	       , n1, n2 
	       , arrayMax (jp->reads)
	       ) ;
    }
  
  if (1)
    {
      DICT *dict = jp->readDict ;
      int i, iMax = dictMax (dict) ;
      for (i = 1, iMax = dictMax (dict) ; i <= iMax ; i++)
	aceOutf (ao, "%s\n", dictName (dict, i)) ;
    }
  ac_free (h) ;
} /* jumpAlignExportReadsDict */

/***************************************************************/
/* export the reads coordinates */
static void jumpAlignExportReadsArray  (JP *jp)
{
  char *fNam = hprintf (jp->h, "%s.jumpalign.reads_array", jp->readFileNameX) ;
  long unsigned int mx, fd =  open (fNam, O_CREAT | O_RDWR | O_BINARY , 0644 );
  if (! fd)
    messcrash ("cannot open the jumpalign.reads_array") ;
  
  mx = arrayMax (jp->reads) ; mx *= sizeof (SEQ) ;
  fdBigWrite (fd, arrp (jp->reads, 0, SEQ), mx, fNam) ;

  close (fd) ;
} /* jumpAlignExportReadsArray */

/***************************************************************/
/* parse the dna in binary format */
static void jumpAlignExportReadsDna  (JP *jp)
{
  long unsigned int pNN = jp->pNN ;
  char *fNam = hprintf (jp->h, "%s.jumpalign.reads_dna", jp->readFileNameX) ;
  int fd = open (fNam, O_CREAT | O_RDWR | O_BINARY , 0644 );
  if (! fd)
    messcrash ("cannot open the jumpalign.reads.dna") ;
  
  fdBigWrite (fd, bigArrp (jp->pDnaArray, 0, long unsigned int), pNN, fNam) ;
  
  close (fd) ;
} /* jumpAlignExportIndexDna */

/***************************************************************/

static void jumpAlignExportReads  (JP *jp)
{
  jumpAlignExportReadsDict (jp) ;
  jumpAlignExportReadsArray  (jp) ;
  jumpAlignExportReadsDna (jp) ;

  fprintf (stderr, "// %s: Exported %d reads contatining %ld base pairs\n"
	   , timeShowNow ()
	   , dictMax (jp->readDict), bigArrayMax (jp->pDnaArray)
	   ) ;

  return ;
} /* jumpAlignExportReads */

/***************************************************************/
/***************************************************************/
/* initialize the jump table zz with bridges to the next ATGC base, i.e. for words of length zero */
static long unsigned int jumpAlignInitTargetJumps (JP *jp)
{
  long unsigned int nRepeats = 0 ;
  long unsigned int i, j, *zp, *zz = jp->zz ;
  long unsigned int NN = jp->NN ;
  char *tDna = 0 ;
  
  tDna = jp->tDna ? jp->tDna : bigArrp (jp->tDnaArray, 0, char) ;

  zz = jp->zz = (long unsigned int *) malloc (NN * sizeof(long unsigned int)) ;
  if (jp->NFILTERS > 1) jp->zz2 = (long unsigned int *) malloc (NN * sizeof(long unsigned int)) ;

  
  bigmemset (zz, 0, NN * sizeof (long unsigned int)) ;

  for (zp = zz, i = 0, j = NN + 1 ; i < NN ; zp++, i++)
    { 
      zz[i] = 0 ;
      switch(tDna[i])
	{ 
	case A_: case T_: case G_: case C_:
  /* initialisation, all 'atgc' words of length zero are repeated at next 'atgc' position (zz) */
  /* all 'atgc' position must be tested for repeats (jump) */
	  if (j < i) { zz[j] = i ; }
	  j = i ;
	  nRepeats++ ;
	  break ;
	default: /* immediately skip ambiguous letters like n or - */
	  break ;
	}
    }
 
  return nRepeats ;
}  /* jumpAlignInitTargetJumps */

/***************************************************************/
/***************************************************************/
/* initialize the jump table with direct hooks into the Target index table */
static unsigned long int jumpAlignInitReadJumpsFromFragmentedOligoTableParttN (const JP *jp, IXLD *ixld)
{
  int suffix = ixld->suffix ;
  int wMax0 = jp->wMax0 ;
  int wMax = jp->wMax ;
  long unsigned int pNN = jp->pNN ;

  register long unsigned int *zzInit = ixld->zzInit ;
  register long unsigned int *zz = jp->pzz ;
  BigArray dna =  jp->pDnaArray ;
  /*   BitSet flags = jp->flags ; */

  register const char *cp ;
  unsigned long int oligo, nn = 0 ;
  unsigned long int mask, mask2 ;
  int ok, pass ;
  long unsigned int ii, iMax,  ln = 1 ;

  long unsigned int shift2 = 2 * (wMax0 - 1) ; /* mask the 2 left bits of the oligo to insert a base from the left */
  unsigned long int A__ = 0 ;
  unsigned long int T__ =  ((unsigned long int) 0x1) << shift2 ;
  unsigned long int G__ =  ((unsigned long int) 0x2) << shift2 ;
  unsigned long int C__ =  ((unsigned long int) 0x3) << shift2 ;


  mask = 1 ;  mask <<= 2 * wMax0 ; mask-- ; /* do it this way to avoid left bit overflow */
  mask2 = mask >> 2 ;
  
   /* compute the oligos, hook into the target DNA coordinates, set the flags */
  iMax = pNN  ;
  zz = jp->pzz  ;

  /* scan backwards so that we can ZBAD the short oligos that will, once extended, overlap with absent values */
  /* ATTENTION ii is unsigned, we replace  ii >= 0 by ii + 1 > 0 */
  for (pass = 0 ; pass < 1 ; pass++)  /* pass in 16 times, this will optimize page breaks when dereferencing in zz */
    for (ii = iMax - 1, cp = bigArrp(dna, ii, char), oligo = 0, ok = 0, ln = 0 ; ii + 1 > 0 ; cp--, ii--)
      {
	ln++ ; /* read length */
	/* construire l'oligo, si pas connu on l'ajoute */
	oligo >>= 2 ; oligo = oligo & mask2 ;  /* mask2 hides the left base */
	switch ((int) *cp)
	  {
	  case A_: oligo |= A__ ; ok++; break ;
	  case T_: oligo |= T__ ; ok++; break ;
	  case G_: oligo |= G__ ; ok++; break ;
	  case C_: oligo |= C__ ; ok++; break ;
	  case 0:  /* end of this sequence */
	    ok = 0 ; 
	    break ;
	    
	  default: ok = 0 ; break ;
	  }
	oligo &= mask ;
	
	if (ok >= wMax0) /* the oligo exists */
	  {
	    register long unsigned int iDebut = ii ;
	    /* use the last 4 bits as a suffix flag and create 16 smaller tables */
	    if (
		(SPLIT ==  16 && (oligo &  0xf) == suffix && zz[iDebut] != ZBAD) ||
		(SPLIT == 256 && (oligo & 0xff) == suffix && zz[iDebut] != ZBAD)
		)
	      {
		register unsigned long int oligo2 = oligo >> (SPLIT == 256 ? 8 : 4) ;
		register long unsigned int z1 = zzInit[oligo2] ;
		if (z1)
		  {
		    nn++ ; 
		    zz[iDebut] = z1 ;
 		    /* bitSet (flags, z1) ; */
		    if (pNN < 220) fprintf (stderr,"%lu %c %lx\t%lu\n", iDebut, dnaDecodeChar[(int)(*cp)], oligo, z1) ;
		  }
		else
		  {  /* kill this oligo and everyone over that will overlap with it at the end of the extension */
		    register long unsigned int  dw, dwMax ;
		    dwMax = wMax > wMax0 ?  wMax - wMax0 + 1 : 1 ;
		    if (dwMax > iDebut) dwMax = iDebut ;
		    for (dw = 0 ; dw < dwMax ; dw++)
		      zz[iDebut - dw] = ZBAD ;
		  }
	      }
	  }
      }
  return nn ;
}  /* jumpAlignInitReadJumpsFromFragmentedOligoTableParttN */

/***************************************************************/
/* open or memorymap one of the index parts */
static void jumpAlignInitIndexer (const void *vp)
{
  const JP *jp = (JP *)vp ;  
  IXLD ixld ; 
  BOOL debug = FALSE ;

  memset (&ixld, 0, sizeof (ixld)) ;
  while (channelMultiGet (jp->ixldChan2, &ixld, 1, IXLD)) /* will break on closeChannel */
    {
      if (debug) fprintf (stderr, "%s : suffix %d indexing\n", timeShowNow (), ixld.suffix) ;
      ixld.nn = jumpAlignInitReadJumpsFromFragmentedOligoTableParttN (jp, &ixld) ;
       if (ixld.fd) 
	 { free (ixld.zzInit) ; close (ixld.fd) ;}
       else
	 munmap (ixld.zzInit, ixld.size) ;
       ixld.zzInit = 0 ;
       ac_free (ixld.fNam) ; ixld.fNam = 0 ;
        if (debug) fprintf (stderr, "%s : ---- suffix %d indexed\n", timeShowNow (), ixld.suffix) ;
	channelPut (jp->ixldChan3, &ixld, IXLD) ;
    }
  return ;
} /* jumpAlignInitIndexer */

/***************************************************************/
/* open or memorymap one of the index parts */
static void jumpAlignInitLoader (const void *vp)
{
  const JP *jp = (JP *)vp ;   
  BOOL useMemoryMapping = FALSE ;
  IXLD ixld ;
  unsigned long int NI = 1 ;
  BOOL debug = FALSE ;

  NI <<= 2*jp->wMax0 ; 
  NI *= sizeof ( long unsigned int) ;
  if (SPLIT == 256)
    NI >>= 8 ; /* since we divide in 256 parts */
  else if (SPLIT == 16)
    NI >>= 4 ; /* since we divide in 16 parts */
  
  memset (&ixld, 0, sizeof (ixld)) ;
  while (channelMultiGet (jp->ixldChan1, &ixld, 1, IXLD)) /* will break on closeChannel */
    {
      if (debug) fprintf (stderr, "%s : suffix %d loading\n", timeShowNow (), ixld.suffix) ;
      ixld.fNam = hprintf (0,  "%s.jumpalign.init_%d_table.%d", jp->indexFileNameI, SPLIT,ixld.suffix) ;
      ixld.fd = open (ixld.fNam, O_RDONLY | O_BINARY , 0644 );
      if (ixld.fd < 0)
	messcrash ("jumpAlignInitLoader suffix==ix %d : fd=%d cannot open the jumpalign.init_table\n%s"
		   , ixld.suffix, ixld.fd, ixld.fNam
		   ) ;
      
      if (useMemoryMapping)
	{
	  ixld.size = NI ;
	  ixld.zzInit =  (long unsigned int *) mmap(NULL, NI, PROT_READ, MAP_PRIVATE,  ixld.fd, 0) ;
	}
      else
	{
	  ixld.zzInit = (long unsigned int *) malloc (NI) ; /* malloc : mo need to zero init */
	  fdBigRead (ixld.fd, ixld.zzInit, NI, ixld.fNam) ;
	}
      /* signal the indexer to work on this block */
       if (debug) fprintf (stderr, "%s : --- suffix %d loaded\n",  timeShowNow (), ixld.suffix) ;
      channelPut (jp->ixldChan2, &ixld, IXLD) ; 
    }
  return ;
} /* jumpAlignInitLoader */

/***************************************************************/

static unsigned long int jumpAlignInitReadJumpsFromOligoTable (JP *jp)
{
  int suffix ;
  unsigned long int nn = 0 ;
  IXLD ixld ;

  memset (&ixld, 0, sizeof (ixld)) ;
  memset (jp->pzz, 0, jp->pNN * sizeof(long unsigned int)) ;

  /* create the communication channels before we wego the subroutines */
  jp->ixldChan1 = channelCreate (256, IXLD, jp->h) ;
  jp->ixldChan2 = channelCreate (1, IXLD, jp->h) ;
  jp->ixldChan3 = channelCreate (256, IXLD, jp->h) ;

  /* launch the codes, they do not yet execute */
  ixld.fNam = hprintf (0,  "%s.jumpalign.init_%d_table.%d", jp->indexFileNameI, SPLIT, ixld.suffix) ;
  ixld.fd = open (ixld.fNam, O_RDONLY | O_BINARY , 0644 );
  if (ixld.fd < 0)
    messcrash ("jumpAlignInitLoader suffix==ix %d : fd=%d cannot open the jumpalign.init_table\n%s"
	       , ixld.suffix, ixld.fd, ixld.fNam
	       ) ;
  else
    {
      messout ("bravo\n") ; close (ixld.fd) ; ixld.fd = 0 ;
    }
  /* we can have 2 init loaders to access the disk */

  wego_go (jumpAlignInitLoader, jp, JP) ;
  wego_go (jumpAlignInitLoader, jp, JP) ;


  /* but we must have a single indexer because we erase above in pzz for other suffixes */
  wego_go (jumpAlignInitIndexer, jp, JP) ;
  wego_go (jumpAlignInitIndexer, jp, JP) ;
  wego_go (jumpAlignInitIndexer, jp, JP) ;
  wego_go (jumpAlignInitIndexer, jp, JP) ;
  

  /* send a collection of signals to the loader */
  for (suffix = 0 ; suffix < SPLIT ; suffix++)
    {
      ixld.suffix = suffix ;
      channelPut (jp->ixldChan1, &ixld, IXLD) ;
    }
  channelClose (jp->ixldChan1) ; /* will stop the child threads */

  /* wait for completion, which will be signaled by the indexer */
  for (suffix = 0 ; suffix < SPLIT ; suffix++)
    {
       memset (&ixld, 0, sizeof (ixld)) ;
       channelMultiGet (jp->ixldChan3, &ixld, 1, IXLD) ;
       /* in reality we do not know which suffix is finished, this is in the struct */
       nn += ixld.nn ;
       fprintf (stderr, "// %s: suffix %d jumpAlignInitReadJumps hooked %lu read-oligos of length %d \n"
		, timeShowNow ()
		, ixld.suffix
		, ixld.nn, jp->wMax0
	   ) ;    
    }

  channelClose (jp->ixldChan2) ;

  /* do not free the channels, the child threads must remain able to read the close signal */

  /*
    for (suffix = 0 ; suffix < SPLIT ; suffix++)
    nn += jumpAlignInitReadJumpsFromFragmentedOligoTableParttN (jp, suffix) ;
  */

  return nn ;
}  /* jumpAlignInitReadJumpsFromOligoTable */

/***************************************************************/
/* hook all valid bases to target base 1 */
static unsigned long int jumpAlignInitReadJumpsTrivially (JP *jp)
{
  unsigned long int nn = 0 ;
  register long unsigned int i = jp->pNN, *zp = jp->pzz - 1 ;
  register const char *cp = bigArrp (jp->pDnaArray, 0, char) - 1 ;
  
  if (i > 0)
    while (cp++, zp++, i--)
      switch ((int) *cp)
	{
	case A_:
	case T_:
	case G_:
	case C_:
	  *zp = 1 ; nn++ ;
	  break ;
	default :
	  *zp = 0 ;
	  break ;
	}

  /* set the flag of the only targetted base  i.e. position 1 */
  bitSet (jp->flags, 1) ;

  return nn ;
} /* jumpAlignInitReadJumpsTrivially */

/***************************************************************/

static void jumpAlignInitReadJumps (JP *jp)
{
  unsigned long int nn = 0 ;
  jp->pzz = (long unsigned int *) malloc (jp->pNN * sizeof(long unsigned int)) ;
  if (jp->NFILTERS > 1)
    jp->pzz2 = (long unsigned int *) malloc (jp->pNN * sizeof(long unsigned int)) ;

  jp->flags  =  bitSetCreate (jp->NN + 1, 0) ;

  if (jp->indexFileNameI)
    nn = jumpAlignInitReadJumpsFromOligoTable (jp) ;
  else
    nn = jumpAlignInitReadJumpsTrivially (jp) ;
  
  fprintf (stderr, "// %s: jumpAlignInitReadJumps hooked %lu read-oligos of length %d into the target\n"
	   , timeShowNow ()
	   , nn, jp->wMax0
	   ) ;
  return ;
} /* jumpAlignInitReadJumps */


/***************************************************************/
/***************************************************************/
/* multithreading in this way is not favorable although it works, 
 * it is slower by a factor 3 or so
 * not sure why
 */
static void jumpAlignDoFilter (PP * pp, ZONE * zone)
{
  int NFILTERS = pp->NFILTERS ;
  int iFilter = zone->iFilter ;
  int w = zone->w ;
  long unsigned int i, iMin, iMax ;
  long unsigned int *zp, *zp2, z1, *zz, *zz2 ;
  BigArray dna = pp->tDnaArray ;
  char cc, ccFilter[2], filters[] = {A_,T_,G_,C_} ;
  const char *tDna ;
  long unsigned int NN = pp->NN ;

  zz = (w & 0x1 ? pp->zz : pp->zz2) ;
  zz2 = (w & 0x1 ? pp->zz2 : pp->zz) ;

  ccFilter[0] = filters[iFilter & 0x3] ;
  ccFilter[1] = filters[(iFilter >> 2) & 0x3] ;
  
  if (! bigArrayExists(dna))
    messcrash (" jumpAlignFilter w=%d, iFilter=%d did not receive a good dna", w, iFilter) ;
  tDna = pp->tDna ? pp->tDna : bigArrp (pp->tDnaArray, 0, char) ;

  
  iMin = zone->iMin ; iMax = zone->iMax ;
  if (iMin > NN - w)
    iMin = NN - w ;
  if (iMax == 0 || iMax > NN - w)
    iMax = NN - w ;
  
  if (1) printf("--filter %d:%c w=%d\n", iFilter, dnaDecodeChar[(int)ccFilter[0]], w) ;
  for (i = iMin, zp = zz + i, zp2 = zz2 + i ; i <= iMax ; i++, zp++, zp2++)
    { 
      if (NFILTERS >= 4)
	{
	  cc = tDna[i] ;
	  if (cc != ccFilter[0])
	    continue ;
	}
      if (NFILTERS >= 16)
	{
	  if (tDna[i+1] != ccFilter[1])
	    continue ;
	}
 
     /* search if  word w[i] of length w starting at i is repeated downwards */
      z1 = *zp ;
      if (!z1)      /* w[i] cannot be repeated if (w-1)[i] is not */
	{
	  continue ;
	}
      if (!zp[1])      /* w[i] cannot be repeated if (w-1)[i+1] is not */
	{
	  continue ;
	}
      
      cc = tDna[i + w - 1 ] ;  /* just compare a single new letter */
      while (1)
	{ 
	  pp->nComp++ ;
	  if (cc == tDna[w - 1 +  z1])    /* found the next repeat */
	    { 
	      pp->nRep++ ; /* count all repeats */ 
	      *zp2 = z1 ;
	      break ;        /* done: found repeat */
	    } 
	  z1 = zz[z1] ;
	  if (!z1) 
	    { 
	      *zp2 = 0 ;
	      break ;
	    }       /* done: no repeat */
	}
    }
  /* request new zone calculation or issue termnination */
  if (1)
    {
      ZONE newZone ;
      BOOL done = TRUE ;

      if (iMax < NN - w) /* new coords, same filter */
	{
	  newZone.w = w ;
	  newZone.iFilter = iFilter ;
	  newZone.iMin = iMax + 1 ;
	  newZone.iMax = iMax + NN/(pp->NFILTERS) ;
	  if (newZone.iMax  > NN - w)
	    newZone.iMax  = NN - w ;
	  channelPut (pp->zoneChan, &newZone, ZONE) ;
	  done = FALSE ;
	}
      if (iMax == 0 && iFilter < pp->NFILTERS)
	{
	  newZone.w = w ;
	  newZone.iFilter = iFilter + 1 ;
	  newZone.iMin = iMin ;
	  newZone.iMax = iMax ;
	  channelPut (pp->zoneChan, &newZone, ZONE) ;
	  done = FALSE ;
	}
      if (1 && done)
	{ /* we cannot close the zoneCham, we need it for the next w */
	  newZone.w = -1 ;
	  channelPut (pp->zoneChan, &newZone, ZONE) ;
	}
      /*  else
	  channelClose (pp->zoneChan) ;
      */
    }
} /* jumpAlignDoFilter */

/***************************************************************/
/***************************************************************/
/* this option is faster than the one above by about 5% */
static void jumpAlignDo (const PP *pp)
{
  long unsigned int i ;
  int w = pp->w ;
  int wMax = pp->wMax ;
  BitSet flags = pp->flags ;
  long unsigned int *zp, z1, *zz = pp->zz, *pzz = pp->pzz ;
  long unsigned int comp = pp->nComp[pp->iFilter] ;
  long unsigned int nHooks = 0 ;
  long unsigned int nRep = 0 ; /* pp->nRep[pp->iFilter] ; */
  char cc, *tDna ;
  BOOL debug = FALSE ;

  if (pp->tDnaArray && ! bigArrayExists(pp->tDnaArray))
    messcrash (" jumpAlignFilter w=%d, iFilter=%d did not receive a good dna", w, 0) ;
  if (pp->pDnaArray && ! bigArrayExists(pp->pDnaArray))
    messcrash (" jumpAlignFilter w=%d, iFilter=%d did not receive a good dna", w, 0) ;
  /* tDNa might be a pointer into a BigArray or into a memory mapped from a file */
  tDna = pp->tDna ? pp->tDna : (pp->tDnaArray ? bigArrp (pp->tDnaArray, 0, char) : 0) ;


  if (debug) printf("-- %s: --start single filter %d w=%d\n", timeShowNow (), 0, w) ;
  if (debug) fprintf (stderr, "// found %lu set flags, should be zero\n", bitSetCount (flags)) ;

  if (0)
    {
      char *pDna = bigArrp (pp->pDnaArray, 0, char) ;  
      for (i = 0 ; i <= pp->pNN  ; i++, zp++)
	fprintf (stderr, "//----- Begin: w=%d %c:%lu->%lu\n"
		 , w, pDna[i], i, pzz[i]) ;
    }

  if (pp->pzz)  /* extend the read part and set the flags */
    {
      int pNN = pp->pNN  ;
      char *pDna = bigArrp (pp->pDnaArray, 0, char) ;
      /* scan backwards so that we can ZBAD the short oligos that will, once extended, overlap with absent values */
      /* ATTENTION ii is unsigned, we replace  i >= 0 by i + 1 > 0 */
      for (i = pNN - 1, zp = pzz + i ; i + 1  > 0  ; i--, zp--)  /* scan the reads backwards */
	{
	  if (! *zp || *zp == ZBAD)
	    continue ;   /* w[i] cannot be repeated if (w-1)[i] is not */
	  if (0 && w > 1 &&  ! *(zp+1)) /*  does not apply because we scan backwards */
	    { *zp = 0 ; continue ;  }  /* w[i] cannot be repeated if (w)[i+1] is not */
	  
	  /* search if  word w[i] of length w starting at i is repeated downwards */
	  z1 = *zp ;     
	  cc = pDna[i + w - 1 ] ;  /* just compare a single new letter */
	  if (!cc) 
	    continue ;
	  while (1)
	    { 
	      comp++ ; /* count all comparisons */
	      if (cc == tDna[w - 1 +  z1])    /* found the next repeat */
		{ 
		  nHooks++ ; /* count all repeats */ 
		  *zp = z1 ;  /* pointer in pZZ */
		  /* bitSet (flags, z1) ; do  this in a specialized loop , it is probably faster */
		  break ;        /* done: found repeat */
		} 
	      z1 = zz[z1] ; /* pointer in target */
	      if (!z1) 
		{ 
		  register long unsigned int  dw, dwMax ;
		  dwMax = wMax > w ?  wMax - w + 1 : 1 ;
		  if (dwMax > i) dwMax = i ;
		  for (dw = 0 ; dw < dwMax ; dw++)
		    *(zp - dw) = ZBAD ; /* pointer in pZZ */
		  break ;
		}       /* done: no repeat */
	    }
	}
      if (debug) printf("-- %s: end read single filter %d w=%d\n", timeShowNow (), 0, w) ;
    }
   
  /* trivial Extend the target jumps as fast as possible */
  if (w == 1)  /* hook each acceptable base to its first repeat */
    {
      register char cc, *cp, *cq ;
      register long unsigned int i, j, *zp, NN = pp->NN - 1  ;

      for (i = 0, zp = zz, cp = tDna ; i <  NN  ; i++, zp++, cp++)  /* scan the reads */
	{
	  if ( *zp)
	    {
	      cc = *cp ;
	      for (j = 1, cq = cp + 1 ; *cq != cc && i+j <= NN ; j++, cq++) {} ; /* find the repeat */
	      *zp = i+j ;
	    }
	}
      *zp = 0 ;
      comp = 0 ;
      nHooks = 4 ; /* there are 4 words of length 1 */
      nRep = NN  ;
      goto done ;
    }

  /* flag the chains of words present in the reads */
  if (pp->pDnaArray)
    {
      register long unsigned int i, NN = pp->NN, *zp ;
      char *pDna = bigArrp (pp->pDnaArray, 0, char) ;

      if (debug) fprintf (stderr, "// %s: start flagging chains of words present in the reads\n", timeShowNow ()) ;
      if (debug) fprintf (stderr, "// found %lu set flags\n", bitSetCount (flags)) ;
      for (i = 0, zp = zz ; i <= NN  ; i++, zp++)  /* scan the reads */
	if (*zp) /*  && bit (flags, i) */
	  { bitSet (flags, *zp) ; }
      if (debug) fprintf (stderr, "// %s: flagged %lu words present in the reads\n", timeShowNow (), bitSetCount (flags)) ;
      if (0)
	for (i = 0, zp = zz ; i <= pp->pNN  ; i++, zp++)
	  fprintf (stderr, "//----- w=%d %c:%lu->%lu\n"
		   , w, pDna[i], i, pzz[i]) ;
    }

  if (tDna)  /* extend the target part */
    { 
      BitSet isRepeat = bitSetCreate (pp->NN + 1, 0) ;
      BitSet flags = pp->pDnaArray ? pp->flags : 0 ;
      register char cc ;
      register long unsigned int i, NN = pp->NN - w, z1, *zp ;
      register BOOL firstPass =  pp->w == pp->wMax0 + 1 ? TRUE : FALSE ;
      register int dw, dwMax = pp->wMax - pp->wMax0 ;
      if (0) fprintf (stderr, "wMax = %d wMax0 = %d dwMax =%d\n", pp->wMax, pp->wMax0, dwMax) ;
      dwMax = 1 ;
      for (i = 0, zp = zz ; i <= NN  ; i++, zp++)  /* scan the genome */
	{
	  if ( ! *zp)
	    continue ;   /* w[i] cannot be repeated if (w-1)[i] is not */
	  if (! *(zp+1)) /*  does not apply if we have cancelled the most frequent hooks */
	    { *zp = 0 ; continue ;  }  /* w[i] cannot be repeated if (w)[i+1] is not */
	  if (flags && (! bit (flags, i) || !  bit (flags, i+1)))
	    { *zp = 0 ; continue ;  }  /* w[i] will not extend to a word present in the reads */
	  if (0) 
	    {
	      for (dw = 0 ; dw < dwMax ; dw++)
		if (! bit (flags, i + dw))
		  dw = dwMax + 2 ;
	      if (dw >= dwMax + 2)
		{ *zp = 0 ; continue ;  }  /* w[i] will not extend to a word present in the reads */
	    }
	  /* search if  word w[i] of length w starting at i is repeated downwards */
	  z1 = *zp ;     
	  cc = tDna[i + w - 1 ] ;  /* just compare a target single new letter */
	  if (!cc) 
	    continue ;
	  while (1)
	    { 
	      comp++ ; /* count all comparisons */
	      if (cc == tDna[w - 1 +  z1])    /* found the next repeat */
		{ 
		  BOOL ok = TRUE ;
		  if (1 && firstPass && flags)
		    {
		      for (dw = 0 ; ok && dw < dwMax ; dw++)
			if (! bit (flags, i + dw))
			  ok = FALSE ;  /* we must not chain towards this repeat */
		    }
		  if (ok)
		    {   
		      if (isRepeat && ! bit(isRepeat, i))
			nHooks++ ;
		      nRep++ ; /* count all repeats */ 
		      *zp = z1 ;
		      if (flags)
			bitSet (flags, z1) ;
		      if (isRepeat)
			bitSet (isRepeat, z1) ;
		      break ;        /* done: found repeat */
		    } 
		}
	      z1 = zz[z1] ;
	      if (!z1) 
		{ 
		  *zp = 0 ;
		  break ;
		}       /* done: no repeat */
	    }
	}
      bitSetDestroy (isRepeat) ;
    }

 done:
  pp->nComp[pp->iFilter] = comp ;
  pp->nHooks[pp->iFilter] = nHooks ;
  pp->nRep[pp->iFilter] = nRep + nHooks ;

  if (debug) printf("-- %s: end target single filter %d w=%d\n", timeShowNow (), 0, w) ;
} /* jumpAlignDo */

/***************************************************************/
/* split the function so the thread does not get killed */ 
static void jumpAlignFilter (const void *vp)
{
  PP *pp = (PP *)vp ;
  ZONE zone ;

  memset (&zone, 0, sizeof(zone)) ;
  while (channelGet (pp->zoneChan, &zone, ZONE))
    {
      if (zone.w < 0)
	break ;
      if (pp->NFILTERS == 1)
	{
	  int w = pp->w = zone.w ;
	  pp->iFilter = zone.iFilter ;
	  jumpAlignDo (pp) ;
	  channelPut (pp->outChan, &w, int) ; /* signal completion */
	}
      else
	{
	  jumpAlignDoFilter (pp, &zone) ;
	  /* outchan is written by the last member of the canon */
	}
    }
  channelPut (pp->outChan, &(zone.w), int) ; /* signal completion */
} /* jumpAlignFilter */

/***************************************************************/
/***************************************************************/

static int jumpAlign (JP *jp)
{
  int NFILTERS = jp->NFILTERS ;
  PP *pp, ppp[NFILTERS] ;
  long unsigned int i ;
  int w = 0, iFilter ;
  long int nRepeats, nNewRepeats, nComparisons, nHooks ;
  char *tDna ;
  int nExactHits = 0 ;
  int wMax = jp->wMax ;
  int readLengthMin = jp->readLengthMin ;
  Array pps = jp->pps ;
  BigArray dna = jp->tDnaArray ;
  long unsigned int NN = jp->NN ;
  BOOL debug = FALSE ;
  CHAN *zoneChan = channelCreate (1, ZONE, jp->h) ;
  int oldw = 0, iiReport = 0 ;
  time_t t2, t1, t0 = time (0)  ;

  pp = arrp (pps, 0, PP) ;
  if (1)
    {
      pp->h = ac_new_handle () ;
      pp->NFILTERS = jp->NFILTERS ;
      pp->outChan = channelCreate (2 + NFILTERS, int, pp->h) ;
    }

  tDna = jp->tDna ? jp->tDna : bigArrp (jp->tDnaArray, 0, char) ;

  t1 = time(0) ;

  /* NN is the length of the dna */
  /* zz[i] is the offset of the first repeat of [dna[i]..dna[i+w-1] */
  /* jump is the next position which may have a repeat */

  nRepeats = nNewRepeats =  nComparisons = 0 ; nHooks = 1 ;

  if (NN < 20) 
    {
      for(i = 0 ; i < NN && i < 100 ; i++)
	{
	  printf("%c\t%lu\t%lu\n", dnaDecodeChar[(int)tDna[i]], i, jp->zz[i]) ;
	}
    }
  
  t2 = time(0) ;
  if (1) fprintf (stderr, "%3d%18lu%18lu%12lu%12lu%8d%8d\n", w, nHooks, nNewRepeats, nRepeats, nComparisons
		 , (int)(t2 - t0), (int)(t2 - t1)) ;
  t1 = t2 ;
  /* first recursion */
  /* nRepeats is the number of repeats of length w */
  /* nNewRepeats is the number of different repeats of length w */
  /* nComparisons is the number of comparison needed during this recursion */

  if (1)
    {
      pp->zz = jp->zz ; pp->zz2 = jp->zz2 ;
      pp->tDnaArray = dna ;
      pp->NN = NN ;
      pp->wMax = jp->wMax ;
      pp->wMax0 = jp->wMax0 ;
      pp->tDna = jp->tDna ;
      pp->pDnaArray = jp->pDnaArray ;
      pp->pNN = jp->pNN ;
      pp->pzz = jp->pzz ; pp->pzz2 = jp->pzz2 ;
      pp->flags = jp->flags ;
      jp->nComp =  pp->nComp =  halloc ( NFILTERS * sizeof (long int), pp->h) ;
      jp->nHooks =  pp->nHooks =  halloc ( NFILTERS * sizeof (long int), pp->h) ;
      jp->nRep = pp->nRep =  halloc ( NFILTERS * sizeof (long int), pp->h) ;
      
      for (iFilter = 0 ; iFilter < NFILTERS ; iFilter++)
	{ 
	  PP *pi = ppp + iFilter ;
	  ppp[iFilter] = *pp ;
	  ppp[iFilter].iFilter = iFilter ;
	  ppp[iFilter].zoneChan = zoneChan ; 
	  ppp[iFilter].nComp = pp->nComp ;
	  ppp[iFilter].nRep = pp->nRep ;
	  wego_go (jumpAlignFilter,  pi, PP) ;
	}
      
      if (1) fprintf (stderr, "// Date\tword length\tDistinct repeated words\tRepeats\tCumul\tComparisons\n") ;
      for (w = jp->wMax0 + 1 ; nHooks + nNewRepeats > 0 && w <= wMax ; w++)
	{
	  pp->w = w ;
	  
	  nNewRepeats = nComparisons = nHooks = 0 ;

	  for (iFilter = 0 ; iFilter < NFILTERS ; iFilter++)
	    {
	      /* we cannot manipulate these like this, since ppp has been copied */
	      pp->nComp = pp->nRep = pp->nHooks = 0 ;
	      ppp[iFilter].nRep = ppp[iFilter].nComp = ppp[iFilter].nHooks = 0  ;
	      ppp[iFilter].zz = jp->zz ; ppp[iFilter].zz2 = jp->zz2 ; 
	    }
	  /* finish all initialisations before channelPut because we do not know which client filter will start first */
	  if (debug)   ppShow (jp, messprintf ("before w=%d", w)) ;

	  if (1) /* actual work */
	    {
	      /* start one execution of the actual work 
	       * the code runs 'en canon'
	       * the end of each zone
	       * crates 2 new channelPut
	       * same zone next fileter plus same fileter next zones
	       * until the word length is studied
	       * and we issue channelPut (pp->outChan)
	       */
	      for (iFilter = 0 ; iFilter < 1 ; iFilter++)
		{
		  ZONE zone ;
		  zone.iFilter = iFilter ;
		  zone.w = w ;
		  zone.iMin = 0 ;
		  zone.iMax = NN/NFILTERS ;
		  channelPut ((ppp+iFilter)->zoneChan, &zone, ZONE) ;
		}
	      /* wait for completion */
	      for (iFilter = 0 ; iFilter < NFILTERS ; iFilter++)
		{
		   int wDone = 0 ;
		   channelGet (pp->outChan, &wDone, int) ;  /* wait on all filters */
		}
	     }

	  /* gather the statistics */
	  nNewRepeats = 0 ;
	  for (iFilter = 0 ; iFilter <  NFILTERS ; iFilter++)
	    {
	      nComparisons += jp->nComp[iFilter] ; jp->nComp[iFilter] = 0 ;
	      nNewRepeats += jp->nRep[iFilter] ; jp->nRep[iFilter] = 0 ;
	      nHooks += jp->nHooks[iFilter] ; jp->nHooks[iFilter] = 0 ;
	    }
	  nRepeats += nNewRepeats ;
	  if (debug)   ppShow (jp, messprintf ("after w=%d", w)) ;

	  t2 = time(0) ; 
	  if (1) fprintf (stderr, "// %s:\tw=%3d\t%18lu\t%18lu\t%12lu\t%12lu\t%8d\t%8d\n"
			  , timeShowNow()
			  , w, nHooks, nNewRepeats, nRepeats, nComparisons
			  , (int)(t2 - t0), (int)(t2 - t1)
			  ) ;
	  t1 = t2 ;
	  if (jp->reads && (jp->exportHits || jp->exportExactHits) && w >= readLengthMin)
	    nExactHits += reportExactReads (jp, w, oldw, &iiReport) ;
	  oldw = w ;
	  if (jp->reads && ! jp->exportIndex && jp->exportHits && (! nHooks || w == jp->readLengthMax))
	    break ;
	}
    }

  if (jp->wMax0 == jp->wMax && jp->reads && jp->exportExactHits)
    reportExactReads (jp, jp->wMax, oldw, &iiReport) ;
  for (iFilter = 0 ; iFilter < NFILTERS ; iFilter++)
    channelClose ((ppp+iFilter)->zoneChan) ;
  
  if (jp->exportIndex)
    jumpAlignExportIndex (jp, jp->wMax) ; /* will free zz */
 
  return nExactHits ;
} /* jumpAlign */

/*************************************************************************************/
/* clip 3' polyA and 5' polyT in excess of proportion  4/5
   by adding a flag in the base byte, this way we
   have c & of = A_ but c != A_
 */
static int jumpAlignSoftClipPolyA (JP *jp, SEQ *seq, int *nbp)
{
  long unsigned int i, j, n, k ; 
  int nClipped = 0, nbClipped = *nbp ;
  register unsigned char *cp ;
  BigArray dna = jp->pDnaArray ;
  
  if (seq->end < seq->start + 8)
    return 0 ;

  for (i = seq->start, n = k = 0, cp = bigArrp (dna, i, unsigned char) ; k < 5 ; k++, cp++)
    if (*cp == T_) n++ ;
  if (5 * n >= 4 * k)
    {  /* found a T rich zone */
      for ( ; *cp ; cp++)
	{
	  k++ ;  if (*cp == T_) n++ ;
	  if (5 * n < 4 * k)
	    break ;
	} 
      /* soft clip the 5' T_ */
      for (i = seq->start, cp = bigArrp (dna, i, unsigned char), j = 0 ; j < k ; j++, cp++)
	*cp |= 0x10 ;
      nbClipped += k ;
      nClipped = 1 ;
    }
    
  for (i = seq->end - 1, n = k = 0, cp = bigArrp (dna, i, unsigned char) ; k < 5 ; k++, cp--)
    if (*cp == A_) n++ ;
  if (5 * n >= 4 * k)
    {  /* found a A rich zone */
      for ( ; *cp ; cp--)
	{
	  k++ ;  if (*cp == A_) n++ ;
	  if (5 * n < 4 * k)
	    break ;
	} 
      /* soft clip the 3' A_ */
      for (i = seq->end - 1, cp = bigArrp (dna, i, unsigned char), j = 0 ; j < k ; j++, cp--)
	*cp |= 0x10 ;
      nbClipped += k ; nClipped = 1 ;
    }     
  *nbp = nbClipped ;
  return nClipped ;
} /*  jumpAlignSoftClipPolyA */

/*************************************************************************************/

static int parseFastaFile (JP *jp, long unsigned int *nbp, BOOL isRead)
{
  AC_HANDLE h = ac_new_handle () ;
  long unsigned int nTarget, nb = 0 ;
  int newNam ;
  SEQ *seq = 0 ;
  BOOL last = FALSE ;
  char *cp, *cr ;
  Array seqs = 0 ;
  BigArray dna =  isRead ? jp->pDnaArray : jp->tDnaArray ;
  BigArray dna2chrom = 0 ;
  ACEIN ai ;
  BOOL justStrand = jp->justStrand ;
  DICT *dict = isRead ? jp->readDict : jp->targetDict ;
  long unsigned int NN = isRead ? jp->pNN : jp->NN ;
  int nClipped = 0, nbClipped = 0 ;
  
  if (isRead)
    {
      dict = jp->readDict = dictHandleCreate (256, jp->h) ;
      seqs = jp->reads  = arrayHandleCreate (100000, SEQ, jp->h) ;
      dna = jp->pDnaArray = bigArrayHandleCreate (10000000, char, jp->h) ;
      ai = aceInCreate (jp->pFileName, jp->gzi , h) ;
      if (!ai)
	usage ("Sorry i cannot open the read fasta file") ;
    }
  else
    {
      dna2chrom = jp->dna2chrom = bigArrayHandleCreate (1000, int, jp->h) ;
      dict = jp->targetDict = dictHandleCreate (256, jp->h) ;
      seqs = jp->chroms = arrayHandleCreate (100, SEQ, jp->h) ;
      dna = jp->tDnaArray = bigArrayHandleCreate (1000, char, jp->h) ;
      if (NN == 0)
	{
	  jp->NN0 = NN = 1 ;  /* leave a 1b gap where we can stick the individual reads */
	  bigArray (dna, NN << 4, char) = 0 ;
	}
      ai = aceInCreate (jp->tFileName, 0, h) ;
      if (!ai)
	usage ("Sorry i cannot open the target fasta file") ;
    }

  nTarget = arrayMax (seqs) ;
  
  while (TRUE)
    {
      if (aceInCard (ai))
	{
	  cp = aceInWord (ai) ;
	  if (!cp || !*cp)
	    continue ;
	  {  /* get rid of fastc format */
	    char *cq = strstr (cp, "><") ;
	    if (cq) *cq = 0 ;
	  }
	}
      else
	last = TRUE ;
      if (last || *cp == '>')
	{
	  /* register the read and
	   * add a null base to propagate correctly bacwards
	   * the clause w[i] cannot be repeated if (w-1)[i+1] is not
	   */
	  if (!last)
	    dictAdd (dict, cp+1, &newNam) ;
	  if (seq)
	    {
	      if (isRead)
		{
		  bigArray (dna, NN++, char) = 0 ; /* 00 terminate each sequence */
		  bigArray (dna, NN++, char) = 0 ;
		}
	      else
		{
		  bigArray (dna, NN++, char) = 0 ; /* 00 terminate each sequence */
		  bigArray (dna, NN++, char) = 0 ;
		  if (NN % 1024) /* start on a new page and associate this page to the current newName */
		    NN = NN - (NN % 1024) + 1024 ;
		  bigArray (dna2chrom, NN >> 10, int) = nTarget ;
		}
	      seq->end = NN - 1 ;
	      nb += seq->ln ;
	    }

	  if (seq && seq->ln && isRead)
	    {
	      if (seq->ln <  jp->readLengthMin)
		jp->readLengthMin = seq->ln ;
	      if (seq->ln > jp->readLengthMax)
		jp->readLengthMax = seq->ln ;
	    }
	  
	  if (seq && isRead && jp->clipPolyA)
	    nClipped += jumpAlignSoftClipPolyA (jp, seq, &nbClipped) ;
	  if (seq && (isRead || jp->exportRepeats)  && !justStrand)
	    { /* add the strand inverted read */
	      char *cr ;
	      long unsigned int j, n0 = NN ;
	      

	      bigArray (dna, NN + seq->ln + 2, char) = 0 ; /* make room */
	      for (cr = bigArrp (dna, seq->start + seq->ln - 1, char), j = seq->ln ; j-- ; cr--)
		bigArray (dna, NN++, char) = complementBase(*cr) ;
	      bigArray (dna, NN++, char) = 0 ;
	      bigArray (dna, NN++, char) = 0 ;
	      if (NN % 1024) /* start on a new page and associate this page to the current newName */
		NN = NN - (NN % 1024) + 1024 ;
	      seq = arrayp (seqs, nTarget++, SEQ) ;
	      seq->nam = (seq-1)->nam ;
	      seq->isDown = FALSE ;
	      seq->ln = (seq-1)->ln ;
	      seq->start = n0 ;
	      seq->end = NN - 1 ;  /* last letter of the read (inclusive ? ) */
	      nb += seq->ln ;
	      bigArray (dna2chrom, NN >> 10, int) = nTarget ;
	    }
	  if (last)
	    break ;
	  seq = arrayp (seqs, nTarget++, SEQ) ;
	  seq->nam = newNam ;
	  seq->isDown = TRUE ;
	  seq->start = NN ;
	}
      else
	{
	  seq->ln += strlen(cp) ;
	  for (cr = cp ; *cr ; cr++)
	    bigArray (dna, NN++, char) = dnaEncodeChar [(int)*cr] ;
	}
    }
  *nbp = nb ;

  if (isRead && jp->reads) 
    arraySort (jp->reads, seqOrder) ; /* sort by read length */

  if (isRead) 
    jp->pNN = NN ;
  else
    {
      long int k, kMax = NN/1024 ;
      int chrom, *ip ;

      jp->NN = NN ;
      
      /* complete the dna2chrom table */
      chrom = bigArray (dna2chrom, kMax - 1, int) ;
      for (k = 0, chrom = 0, ip = bigArrayp (dna2chrom, 0, int) ; k < kMax ; k++, ip++)
	if (*ip) 
	  chrom = *ip ;
	else
	  *ip = chrom ;
    } 
  ac_free (h) ;
  
  fprintf (stderr, "// %s: Parsed %lu %s, %lu bases\n", timeShowNow (), nTarget, isRead ? "fragments" : "targets", nb) ;
  if (jp->clipPolyA)
    fprintf (stderr, "//   clipped %d bases in %d reads\n", nbClipped, nClipped) ;
  return nTarget ;
} /* parseFastaFile  */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *err)
{
  fprintf (stderr, "// Usage: jumpalign  aligner program based on iterative jumps\n") ;
  fprintf (stderr, 
	   "// Example:\n"
           "//     jumpalign  -wMax 14 -t ~/NB/TARGET/CHROMS/hs.chrom_2.fasta.gz -exportIndex toto.mask.chrom.2 \n"
           "//     jumpalign  -wMax 30  -importIndex toto.mask.chrom.2 \n"
	   "// -i : the name of a fasta file containing all the reads to be aligned \n"
	   "//      (tested with 1 Million reads, upper limit depends on harware)\n"
	   "// -t fileName : the name of a target fasta file on which we align the reads\n"
	   "//      alternatively, one can align against a precomputed index using -importIndex option\n"
	   "// -o fileNamePrefix : all ouput files start with this prefix\n"
	   "// -strand : if this option is specified the read is aligned on the top strand of the target,\n"
	   "//            this otion should probably be used if the target is mRNA\n"
	   "//            probably not if th target is a genomic DNA file\n"
	   "// -antistrand : if this option is specified the read is aligned on the down strand of the target,\n"
	   "//            this otion should probably be used if the target is mRNA\n"
	   "//            probably not if th target is a genomic DNA file\n"
	   "// -exportIndex fileNamePrefix : a prefix for the precomputed index tables\n"
	   "//       in this case option -i is ignored, the program stops after indexing the target\n"
	   "//       A good value when computing the index is -wMax 14\n"
	   "//       Higher value may block the computer since the needed memory grows like  2^2.wMax\n"
	   "// -importIndex fileNamePrefix : a prefix for the precomputed index tables\n"
	   "//       since the index contains the target DNA, the option -t is ignored\n"
	   "// -exportReads fileNamePrefix : a prefix for the binary reads tables\n"
	   "// -importReads fileNamePrefix : a prefix for the binary reads tables\n"
	   "//       It is way faster to read a binary sequence rather than parse a fasta file\n"
	   "//       so this option is very useful is the same sequences are aligned several times\n"
	   "//n"
	   "// -wMax <int n> : compute the bridges up to words of length n\n" 
	   "//       use n <= 14 in -exportIndex mode, but later any value can be used (say 150)\n"
	   "// -exportHits  -exportExactHits : export hits of length wMax\n"
	   "// -exportRepeats [-: export repeats of length at least wMax\n"
	   "//       -nFilters  : [n=1|4|16, default 1] compute en canon int n voices\n"
	   "// -exportHash : draft code, do not use\n"
	   "// -clipPolyA : soft clip terminal polyA in the reads\n"
	   "// -max_threads <int> : number of simultaneoulsy running functions\n"
	   "//       the program uses concurrent routines synchronized by GO-like channels\n"
	   "//       and accelerates is several concurrent routines actually run in parallel\n"
	   "// -gzo : gzip the large output files\n"
	   "// -gzi : the input file is gzipped\n"
	   "//        gunzip is applied by default on input files called *.gz\n"
	   "//        but this option is useful to gunzip stdin or a file named otherwise\n"
	   ) ;

  if (err)
    fprintf (stderr, "\n// ERROR: %s\n", err) ;
  exit (1) ;
}

/*************************************************************************************/

int main (int argc, const char **argv)
{
  AC_HANDLE h = ac_new_handle () ;
  long unsigned int np=0, npb=0, ntb=0, nFound = 0 ;
  JP jp ;

  if (argc == 1 ||
      getCmdLineBool (&argc, argv, "-h") ||
      getCmdLineBool (&argc, argv, "-help") ||
      getCmdLineBool (&argc, argv, "--help")
      )
    usage (0) ;

  memset (&jp, 0, sizeof (JP)) ;
  if (sizeof (size_t) < 8) messcrash ("Sorry sizeof (size_t) = %d < 8, we are not on a 64-bits machine", sizeof (size_t)) ;
  jp.readLengthMin = 99000000 ;
  jp.readLengthMax = 0 ;
  jp.NFILTERS = 1 ;

  jp.justStrand = getCmdLineBool (&argc, argv, "-strand") ;
  jp.justAntistrand = getCmdLineBool (&argc, argv, "-antistrand") ;

  jp.gzi = getCmdLineBool (&argc, argv, "-gzi") ;
  jp.gzo = getCmdLineBool (&argc, argv, "-gzo") ;
 
  jp.clipPolyA = getCmdLineBool (&argc, argv, "-clipPolyA") ;

  getCmdLineOption (&argc, argv, "-o", &(jp.outFileName)) ;
  getCmdLineOption (&argc, argv, "-t", &(jp.tFileName)) ;
  getCmdLineOption (&argc, argv, "-i", &(jp.pFileName)) ;

  getCmdLineInt (&argc, argv, "-nFilters", &jp.NFILTERS) ;

  if (jp.NFILTERS >= 16) jp.NFILTERS = 16 ;
  else if (jp.NFILTERS >= 4) jp.NFILTERS = 4 ;
  else jp.NFILTERS  = 1 ;
  jp.max_threads = jp.NFILTERS ; jp.max_threads = 5 ;
  getCmdLineInt (&argc, argv, "-max_threads", &jp.max_threads) ;
  getCmdLineInt (&argc, argv, "-wMax", &jp.wMax) ;
  jp.wMax0 = 0 ;
  if (jp.wMax < 3) jp.wMax = 3 ;

  jp.exportHits = getCmdLineBool (&argc, argv, "-exportHits") ;
  jp.exportRepeats = getCmdLineBool (&argc, argv, "-exportRepeats") ;
  jp.exportExactHits = getCmdLineBool (&argc, argv, "-exportExactHits") ;
  jp.importIndex = getCmdLineOption (&argc, argv, "-importIndex", &jp.indexFileNameI) ;
  jp.exportIndex = getCmdLineOption (&argc, argv, "-exportIndex", &jp.indexFileNameX) ;
  getCmdLineInt (&argc, argv, "-exportHash", &jp.exportHash) ;

  if (! jp.exportRepeats && jp.NFILTERS > 1)
    messcrash ("The code is no longer correct for NFLTERS > 1") ;

  if (jp.exportIndex)
    {
      jp.exportIndex = TRUE ;
      jp.NFILTERS = 1 ;
      if (! jp.indexFileNameX)
	messcrash ("option -exportIndex requires option -mask <mask_file_name>") ;
    }
  if (jp.importIndex)
    {
      if (! jp.indexFileNameI)
	 messcrash ("option -importIndex requires option -mask <mask_file_name>") ;
    } 
  jp.importReads = getCmdLineOption (&argc, argv, "-importReads", &jp.readFileNameI) ;
  jp.exportReads = getCmdLineOption (&argc, argv, "-exportReads", &jp.readFileNameX) ;

 

  jp.h = ac_new_handle () ;
  jp.pps = arrayHandleCreate (100, PP, h) ;
  
  if (argc > 1)
    usage (messprintf ("Unknown parameter %s", argv[1])) ;

  fprintf (stderr, "// %s: jumpalign start\n", timeShowNow ()) ;
 
  /* actual wort */ 
  if (jp.importReads)
    {
      jumpAlignImportReads  (&jp) ;
    }
  else if (jp.pFileName)
    {
      np = parseFastaFile (&jp, &npb, TRUE) ;
      if (np && jp.exportReads)
	jumpAlignExportReads  (&jp) ;
      if (jp.max_threads <  jp.NFILTERS + 4)
	jp.max_threads =  jp.NFILTERS + 4 ;
   }

  wego_max_threads (jp.max_threads) ;

 /* open and parse the reads and the target file */
  if (jp.importIndex)
    jumpAlignImportIndex  (&jp) ;
  else if (jp.tFileName)
    {
      parseFastaFile (&jp, &ntb, FALSE) ;
      jumpAlignInitTargetJumps (&jp) ;
    }
  else
    messcrash ("missing option -t, please try jumpalign --help") ;
 
  /* initialize the read jump table and the flags */
  if (jp.pDnaArray && (jp.exportHits || jp.exportExactHits))
     jumpAlignInitReadJumps (&jp) ;

  fprintf (stderr, "// %s:  jumpalign starts\n"
	   , timeShowNow ()
	   ) ;

  if (jp.zz && 
      (jp.exportHits || jp.exportExactHits || jp.exportIndex || jp.exportHash || jp.exportRepeats)
      )
    {
      int nReads = 0, nHits = 0 ;

      jp.aoHits = aceOutCreate (jp.outFileName, ".hits", jp.gzo, jp.h) ;

      nFound = jumpAlign (&jp) ;  /* construct the jump extensions from wMax0 up to wMax */
      if (jp.reads)
	nFound += reportHits (&jp, &nReads, &nHits) ;
      if (jp.exportHash)
	jumpAlignExportHash (&jp) ;
      if (jp.exportRepeats)
	jumpAlignExportRepeats (&jp) ;
      fprintf (stderr, "// %s: found %d reads alignmed in %d poitions\n ", timeShowNow (), nReads, nHits) ;
    }
  fprintf (stderr, "// %s: jumpalign done, found %lu alignments\n ", timeShowNow (), nFound) ;

  ac_free (h) ;
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

