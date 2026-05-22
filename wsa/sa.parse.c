/*
 * sa.sequenceParse.c

 * This module is part of the sortalign package
 * Created May 15, 2026

 * Authors: Danielle Thierry-Mieg, Jean Thierry-Mieg, Greg Boratyn, NCBI/NLM/NIH

 * This code is public.

 * This module downloads runs from the SRA sequence archives
 * And parses local sequence files in fasta/fastq format
 * into SraParse structures
 */

#include "sa.h"
#include <stdint.h>

/* download SRA runs as FASTA or FASTQ */
#define SRA_FASTA 0
#define SRA_FASTQ 1

/********************************************************************************************************************/
/*************************** encode SAPARSE into standard dna arrays ************************************************/
/********************************************************************************************************************/

/*  scan the saParse buffers, encode the dna, create virtual dna arrays */
static int saScanDnaEncode (BB *bb)
{
  Array dnas = bb->dnas ;
  SAPARSE *sap = bb->saParse ;
  int nn = 0 ;
  int ii, iMax = arrayMax (dnas) ;
  Array *dnap = iMax ? arrp (dnas, 0, Array) : 0 ;

  
  int minDnaLn = (0x1 << 30) ;
  int maxDnaLn = 0 ;
  
  for (ii = 0 ; ii < iMax ; ii++, dnap++)
    {
      if (! *dnap)
	continue ; 
      int n = arrayMax (*dnap) ;
      dnaEncodeArray (*dnap) ;
      int n2 = ustrlen (arrp (*dnap, 0, unsigned char)) ;
      if (n2 < n)
	{
	  messcrash ("Wrong character at position %d in dna sequence\n"
		     , n2
		     ) ;
	}
      minDnaLn = (n < minDnaLn ? n : minDnaLn) ;
      maxDnaLn = (n > maxDnaLn ? n : maxDnaLn) ;
      nn += n ;
    }
  sap->nBases = nn ;
  bb->runStat.p.minReadLength = sap->minDnaLn = minDnaLn ;
  bb->runStat.p.minReadLength = sap->maxDnaLn = maxDnaLn ;
  
  return nn ;
} /* saScanDnaEncode */

/********************************************************************************************************************/
/* create virtual dna arrays pointing inside the dnaBuffer
 * to allow transparent usage of the standard acedb dna array libraries
 */

static Array saScanVirtualDnasArrayCreate (BB *bb)
{
  SAPARSE *saParse = bb->saParse ;
  Array dnas ;
  int ii, iMax = saParse->nRecords ;

  Array dnaStructs  = arrayHandleCreate (iMax, struct ArrayStruct, bb->h) ;
  dnas  = arrayHandleCreate (iMax, Array, bb->h) ;
  for (ii = 0 ; ii < iMax ; ii++)
    {
      DnaRecord *r = arrp (bb->dnaRecords, ii, DnaRecord) ;
      unsigned char *cp = saParse->dnaBuffer + r->xDna ;
      int ln = r->dnaLn ;
      if (ln)
	{
	  if ((uintptr_t)cp & 0xf)
	    messcrash ("dna = %d offset is not aligned on a 16 bytes boundary") ;			 
	  struct ArrayStruct *dnaStruct = arrayp (dnaStructs, ii, struct ArrayStruct) ;
	  virtualArrayStructFill (dnaStruct, cp,  ln, unsigned char) ;
	  array (dnas, ii, Array) = (Array) (dnaStruct) ;
	}
    }
  return dnas ;
} /* saScanVirtualArrayCreate */

/********************************************************************************************************************/
/******************************** sequence file parsing utilities ***************************************************/
/********************************************************************************************************************/
/* readScanId
 * Works on a local copy of a fraction of a user provided private BB buffer. Edit in place.
 * Expects to start on >, fasta, or @, fastc
 * Sets suffix (0, 1, 2) for paired end case
 * On return buf (or buf+1 if prefix) is the clean name for the pair stripping /1
 * Returns number of consumed char, including \n
 */
static int readScanId (unsigned char *buf, char prefix, int *suffixp)
{
  unsigned char *cp = buf, *cq ;    
  int ln = 0 ;
  
  if (prefix)
    {
      if (prefix != *cp)    /* cp should point to prefix '>' or '@' */
	{
	  messerror ("Wrong identifier prefix %s\n", cp) ;
	  return -1 ;
	}
      else
	cp++ ;    /* skip prefix */
    }
  cq = ustrchr (cp, '\n') ;
  if (cq)
    {
      *cq = 0 ;
      ln = cq - buf + 1 ;
    }
  else
    ln = ustrlen (buf) ;

  if (suffixp)
    {
      int n = ustrlen (cp) ;
      *suffixp  = 0 ;
      if (n > 2)
	{
	  cq = cp + n - 2 ;
	  if ((*cq == '/' || *cq == '.' || *cq == ' ')  && (cq[1] == '1' || cq[1] == '2'))
	    {
	      *suffixp = cq[1] - '0' ;
	      *cq = 0 ;
	    }
	}
    }
  return ln ;
} /* readScanId */

/******************************************************************/
/* readScanQual
 * accept a full single line of fastq quality factors
 * Returns number of consumed char, including \n
 */
static int readScanQualId (unsigned char *buf)
{
  int ln = 0 ;
  unsigned char *cq ;

  /* scan the quality identifier and drop it */
  cq = ustrchr (buf, '\n') ;
  if (cq)
    {
      *cq = 0 ;
      ln = cq - buf + 1 ;
    }
  else
    ln = ustrlen (buf) ;
  return ln ;
} /* readScanQualId */

/******************************************************************/
/* readScanQual
 * accept a full single line of fastq quality factors
 * Returns number of consumed char, including \n
 */
static int readScanQual (unsigned char *buf)
{
  int ln1 = 0, ln2 = 0 ;
  unsigned char *cq ;

  ln1 = readScanQualId (buf) ;
  /* scan quality */
  if (ln1)
    {
      cq = ustrchr (buf + ln1, '\n') ;
      if (cq)
	{
	  *cq = 0 ;
	  ln2 = cq - buf - ln1 + 1 ;
	}
      else
	ln2 = ustrlen (buf + ln1) ;
    }
  return ln1 + ln2 ;
} /* readScanQual */

/******************************************************************/
/* readScanDna
 * accept full lines untill terminal 0 or > or @
 * Returns number of consumed char, including \n
 */
static int readScanDna (unsigned char *buf)
{
  int ln = 0 ;
  unsigned char *cp = buf - 1, *cq = buf ;

  while (*++cp)
    {
      ln++ ;
      switch ((int)*cp)
	{
	case 0: 
	case '@':
	case '>': *--cp = 0 ; goto done ; break ;
	case '\r': break ;
	case '\n': break ;
	default: *cq++ = *cp ; break ;
	}
    }
 done:
  *cq = 0 ;
  ln = cp - buf + 1 ;
  return ln ;
} /* readScanDna */

/********************************************************************************************************************/

static int saRegisterId (BB *bb, unsigned char *cp, char prefix, Array idArray, int type)
{
  Array dnaRecords = bb->dnaRecords ; 
  int nR = arrayMax (dnaRecords) ;
  DnaRecord *r  = arrayp (dnaRecords, nR, DnaRecord) ;
  int suffix = 0 ;
  int nn = readScanId (cp, prefix, bb->nPairs ? &suffix : 0) ;
  int k, ln = ustrlen (cp + 1) ;
  if (0)
    {
      r->xId = arrayMax (idArray) + 1 ;
      unsigned char *cq = arrayp (idArray, r->xId + ln + 1, unsigned char) ;  /* make room */
      cq = arrayp (idArray, r->xId, unsigned char) ; 
      memcpy (cq, cp + 1, ln + 1) ;
    }
  dictAdd (bb->dict, (char *)(cp + 1), &k) ;
  if (2 * k + type != nR)
    {
      if (type)
	messcrash ("Non matching identifiers in a pair %s and %s\n"
		   , dictName (bb->dict, k - 1)
		   , dictName (bb->dict, k)
		   ) ;
      else
	messcrash ("synchrony problem parsing identifiers %s\n", cp) ;
    }
  if (bb->nPairs)
    {
      if (type == 0 && suffix != 1)
	messcrash ("first identifier of a pair should end as \' 1\', or \'.1\' or \'/1\' : %s\n", cp) ;
      if (type == 1 && suffix != 2)
	messcrash ("second identifier of a pair should end as \' 2\', or \'.2\' or \'/2\' : %s\n", cp) ;
    }
  return nn ;
}

/********************************************************************************************************************/

static int saRegisterDna (BB *bb, unsigned char *cp, Array dnaArray, int type)
{
  Array dnaRecords = bb->dnaRecords ; 
  int nR = arrayMax (dnaRecords) ;
  DnaRecord *r  = arrayp (dnaRecords, nR - 1, DnaRecord) ;

  int nn = readScanDna (cp) ;
  r->xDna = arrayMax (dnaArray) + 1 ;
  r->xDna += 15 ; r->xDna &= ~0xf ;  /* align on 16 bytes boundary */
  int ln = ustrlen (cp) ;
  r->dnaLn = ln ;
  unsigned char *cq = arrayp (dnaArray, r->xDna + ln + 1, unsigned char) ;  /* make room */
  cq = arrayp (dnaArray, r->xDna, unsigned char) ; 
  memcpy (cq, cp, ln + 1) ;

  bb->nSeqs++ ;
  bb->length += ln ;

  bb->runStat.p.nReads++ ;

  bb->runStat.p.nBase1 += type * ln ;
  bb->runStat.p.nBase2 += (1 - type) * ln ;

  return nn ;
}

/********************************************************************************************************************/

static int saRegisterQuality (BB *bb, unsigned char *cp, Array qualityArray)
{
  Array dnaRecords = bb->dnaRecords ; 
  int nR = arrayMax (dnaRecords) ;
  DnaRecord *r  = arrayp (dnaRecords, nR - 1, DnaRecord) ;

  int nn = readScanQual (cp) ;

  if (qualityArray)  /* we may not need the qualities */
    {
      int ln = ustrlen (cp) ;
      r->xQual = arrayMax (qualityArray) + 1 ;
      unsigned char *cq = arrayp (qualityArray, r->xQual + ln + 1, unsigned char) ;  /* make room */
      cq = arrayp (qualityArray, r->xQual, unsigned char) ; 
      memcpy (cq, cp, ln + 1) ;
    }
  return nn ;
}

/********************************************************************************************************************/
/* parse the fasta/fastq buffers into SAPRSE format
 * the 2 buffers contain complete records and start on > or @
 */

static void saParseR12Buffers (const PP *pp, BB *bb)
{
  DnaFormat format = bb->rc.format ;
  unsigned char *buf1 = bb->r1Buffer ;
  unsigned char *buf2 = bb->r2Buffer ;
  SAPARSE *saParse = bb->saParse = halloc (sizeof(SAPARSE), bb->h) ;
  memset (saParse, 0, sizeof (SAPARSE)) ;
  
  if (! bb->r1Buffer) /* no buffer */
    return ;

  unsigned char prefix = '>' ;
  switch (format)
    {
    case FASTQ2:
    case FASTQ:
      prefix = '@' ;
      break ;
    default:
      break ;
    }

  /* we do not need to create baParse->records */
  int nn = (pp->BMAX)<<14 ;  /* estimate 64 bases per record */
  bb->dnaRecords = arrayHandleCreate (nn, DnaRecord, bb->h) ;
  Array idArray = 0 ; // arrayHandleCreate (32 * nn, unsigned char, bb->h) ;
  Array dnaArray = arrayHandleCreate (256 * nn, unsigned char, bb->h) ;
  Array qualityArray = 	(pp->exportSamQuality ? arrayHandleCreate (128 * nn, unsigned char, bb->h) : 0) ;

  unsigned char *cp1 = buf1 ;
  unsigned char *cp2 = buf2 ;
  if (buf2 || format == FASTA2 || format == FASTQ2)
    bb->nPairs = 1 ;
  array (bb->dnaRecords, 1, DnaRecord).xId = 0 ;
  while (*cp1 == prefix)
    {
      if (buf1)
	{
	  cp1 += saRegisterId (bb, cp1, prefix, idArray, 0) ; /* register id or read 1 */
	  cp1 += saRegisterDna (bb, cp1, dnaArray, 0) ; /* register dna of read 1 */
	  if (format == FASTQ || format == FASTQ2)
	    cp1 += saRegisterQuality (bb, cp1, qualityArray) ;      /* register quality of read 1 */
	}
      
      if (buf2)
	{
	  bb->nPairs++ ;
	  cp2 += saRegisterId (bb, cp2, prefix, idArray, 1) ; /* register id or read 2 */
	  cp2 += saRegisterDna (bb, cp2, dnaArray, 1) ; /* register dna of read 1 */
	  if (format == FASTQ)
	    cp2 += saRegisterQuality (bb, cp2, qualityArray) ;      /* register quality of read 2 */
	}

      if (*cp1 && (format == FASTA2 || format == FASTQ2))
	{
	  bb->nPairs++ ;
	  bb->runStat.p.nPairs++ ;
	  cp1 += saRegisterId (bb, cp1, prefix, idArray, 1) ; /* register id or read 2 */
	  cp1 += saRegisterDna (bb, cp1, dnaArray, 1) ; /* register dna of read 2 */
	  if (format == FASTQ2)
	    cp1 += saRegisterQuality (bb, cp1, qualityArray) ;      /* register quality of read 2 */
	}
      if (! bb->nPairs)
	array (bb->dnaRecords, arrayMax (bb->dnaRecords), DnaRecord).xId = 0 ;
    }
  if (buf2 || format == FASTA2 || format == FASTQ2)
    bb->nPairs -= 1 ;
  saParse->nRecords = arrayMax (bb->dnaRecords) ;
  saParse->idBuffer = idArray ? arrayp (idArray, 0, unsigned char) : 0 ;
  saParse->dnaBuffer = arrayp (dnaArray, 0, unsigned char) ;
  saParse->qualityBuffer = (qualityArray ? arrayp (qualityArray, 0, unsigned char) : 0) ;
} /* saParseR12Buffers */

/********************************************************************************************************************/
/*  scan the saParse buffers, encode the dna, create virtual dna arrays */
void saScan (const PP *pp, BB *bb)
{
  SAPARSE *saParse = bb->saParse ;
    
  bb->errors = arrayHandleCreate (256, int, bb->h) ;
  bb->txt1 = vtxtHandleCreate (bb->h) ;
  bb->txt2 = vtxtHandleCreate (bb->h) ;
  bb->length = 0 ;
  bb->runStat.p.lengthDistribution = arrayHandleCreate (1024, long int, bb->h) ;
  bb->runStat.insertLengthDistribution = arrayHandleCreate (1024, long int, bb->h) ;
  bb->nSeqs = 0 ;
  bb->errDict = dictHandleCreate (100000, bb->h) ;
  bb->cpuStats = arrayHandleCreate (128, CpuSTAT, bb->h) ;
  bb->dict = dictHandleCreate (10000, bb->h) ;
  
  if (saParse) /* created by sra_streaming */
    {
      
      int nn = saParse->nRecords ;
      if (nn)
	bb->dnaRecords = virtualArrayHandleCreate (saParse->records, nn, DnaRecord, bb->h) ;
    }
  else
    saParseR12Buffers (pp, bb) ;

  /* create dna virtual arrays to benefit from the acedb dna array library without calling malloc a zillion times */
  bb->dnas = saScanVirtualDnasArrayCreate (bb) ;
  saScanDnaEncode (bb) ;
  return ;
} /* saScan */

/********************************************************************************************************************/
/******************************** sra streaming and caching  ********************************************************/
/********************************************************************************************************************/
/* may 20 2026
 *  f = gzopen ("gilname.gz", r)
 *  cp = gzread (f, 100Mega, buffer) ;
 *  cq = strrchr (cp, '>')
 *    if (!cq || cq == cp) continue reading intil EOF
 *    else { *(cq-1)=0; pass the cp buffer to a new agent which will decode the dna ;
 *    copy cq ... (n bytes)  to a new clean buffer and gzread in buffer+n
 *  gzclose() 
 * the general idea is that parsing big buffers is fast, while decoding them is slow
 * so in this way, even when facing a single large fasta/fastq file,
 * the pipeline will no longer be hanged on the parser
 */

static void saParseFastac (const PP *pp, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  BB b, *bb = 0 ;
  int BMAX = (pp->BMAX << 20) ;
  unsigned char *buffer = halloc (BMAX, h) ;
  unsigned char *buffer2 = halloc (BMAX, h) ;
  int pos = 0 ;
  BOOL done = FALSE ;
  long int nBytes = 0 ;
  int nPuts = 0 ;
  BOOL pairedEnd = bb ? bb->rc.pairedEnd : FALSE ;
  CHAN *chan = pp->plChan ;
  gzFile file = 0 ;
  BOOL debug = FALSE ;
  
  DnaFormat format = rc->format ;
  const char *fileName1 = rc->fileName1 ;
  char tBuf[25] ;
  clock_t t1, t2 ;
  unsigned char prefix = '>' ;
  switch (format)
    {
    case FASTQ2:
    case FASTQ:
      prefix = '@' ;
      break ;
    default:
      break ;
    }

  t1 = clock () ;
  
  file = gzopen (fileName1, "r");
  if (! file)
    messcrash ("\ncannot gzopen target file %s", fileName1) ;

  while (!done)
    {
      int err = 0 ;                    
      int bytes = gzread (file, buffer + pos, BMAX - pos) ;
      unsigned char *cp, *cq ;

      bytes += pos ;
      if (bytes < BMAX)
	{
	  done = TRUE ;
	  if (! gzeof (file)) 
	    messcrash ("Error %s in gzread %s", gzerror (file, & err), fileName1) ;                
        }
      else
	{
	  /* search for beginning of last probably partial sequence */
	  cp = buffer + bytes ; /* just after last byte read */
	  pos = 0 ;
	moveBack:
	  cp-- ;  
	  while (cp > buffer && *cp != '\n')
	    cp-- ;
	  if (cp == buffer)
	    messcrash ("saParseFastac (%s) found a read > BMAX=%d", fileName1, BMAX) ;
	  if (cp[1] != prefix)
	    goto moveBack ;

	  switch (format)
	    {
	    default: break ;
	    case FASTA2:
	      cq = ustrchr (cp+1, '\n') ;
	      if (!cq) /* not a full ID, cannot decide if this is a read 1. Move back 1 record */
		goto moveBack ;
	      if (cq[-1] != '1') /* not read 1 of a pair */
		goto moveBack ;
	      break ;
	    case FASTQ:
	    case FASTQ2:
	      /* i need to be sure that I am on a new record, not on a quality line */
	      cq = cp - 2 ;
	      while (cq > buffer && *cq != '\n' && dnaDecodeChar[(int)(*cp)])
		cq-- ;
	      if (*cq == '\n') /* we were on a sequence line, we need to move back again */
		goto moveBack ;
	      if (format == FASTQ2)
		{
		  cq = ustrchr (cp+1, '\n') ;
		  if (!cq) /* not a full ID, cannot decide if this is a read 1. Move back 1 record */
		    goto moveBack ;
		  if (cq[-1] != '1') /* not read 1 of a pair */
		    goto moveBack ;
		}
	    }

	  cp++ ;  /* move behind \n */
	  pos = buffer + bytes - cp ; /* number of bytes to preserve */

	  memcpy (buffer2, cp, pos) ; /* copy the remnant */
	  bytes -= pos ;
	  cp[0] = 0 ;
	}
      nBytes += bytes ;
      if (nBytes <= 0)
	messcrash ("No sequence found in file %s\n", fileName1) ;
      

      /* create a data block */
      bb = &b ;
      memset (bb, 0, sizeof (BB)) ;
      bb->h = ac_new_handle () ;
	  
      bb->rc.format = format ;
      bb->rc.pairedEnd = pairedEnd ;
      bb->rc.jump5r1 = rc ? rc->jump5r1 : 0 ;
      bb->rc.jump5r2 = rc ? rc->jump5r2 : 0 ;      

      bb->run = rc ? rc->run : 0 ;
      bb->start = timeNow () ;
      /*
	bb->lane = atomic_fetch_add (rc ? &(rc->lane) : &lane, 1) + 1 ;
      */

      bb->lane = atomic_fetch_add_explicit (arrp (pp->runLanes, bb->run, atomic_int), 1, memory_order_relaxed) + 1 ;

      bb->cpuStats = arrayHandleCreate (128, CpuSTAT, bb->h) ;
      bb->rc.fileName1 = fileName1 ;
      /* copy the buffer */
      bb->r1Buffer = halloc (bytes + 32, bb->h) ;
      memcpy (bb->r1Buffer, buffer, bytes) ;
      if (1)
	bb->r1Buffer[bytes] = '0' ;
      else
	{
	  for (int i = 0 ; i < 32 ; i++)
	    bb->r1Buffer[bytes+i] = 0 ;
	}
      bb->r1Buffer[bytes] = '\n' ;
      /* position the remnant */
      if (! done) memcpy (buffer, buffer2, pos) ;
      bb->nSeqs = 100 ;  /* a guess */
      
      /* export the databalock to the channel */
      nPuts++ ;
      channelPut (chan, bb, BB) ;
    }
  channelPut (pp->npChan, &nPuts, int) ; /* global counting of BB blocks accross all sequenceParser agents */
  
  gzclose (file) ;
  ac_free (h) ;
  
  t2 = clock () ;
  saCpuStatRegister ("2.FastaSequenceParser", pp->agent, bb->cpuStats, t1, t2, nBytes) ;

  if (debug)
    {
     int lane = atomic_fetch_add (arrp (pp->runLanes, bb->run, atomic_int), 0) ;
     printf ("--- %s: Stop FastaSequenceParser %d blocks %ld bytes file %s\n", timeBufShowNow (tBuf), lane, nBytes, fileName1) ;
    }
  
  return ;
} /* saParseFastac */

/********************************************************************************************************************/
/******************************** sra streaming and caching  ********************************************************/
/********************************************************************************************************************/

static void sraCachingOut (ACEOUT *ao1p, ACEOUT *ao2p, const char *seq, const char *seq2, const char *sraID
			   , BOOL fastq, BOOL isPaired, BOOL split_pairs, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  char *fNam = 0, *fNam2 = 0 ;

  if (! isPaired)
    split_pairs = FALSE ;

  if (! isPaired)
    {
      if (fastq)
	fNam = hprintf (h, "SRA/%s.fastq", sraID) ;
      else
	fNam = hprintf (h, "SRA/%s.fasta", sraID) ;
    }
  else if (isPaired && !split_pairs)
    {
      if (fastq)
	fNam = hprintf (h, "SRA/%s.sample_12.fastq", sraID) ;
      else
	fNam = hprintf (h, "SRA/%s.sample_12.fasta", sraID) ;
    }
  else /* isPaired && split_pairs */
    {
      if (fastq)
	{
	  fNam  = hprintf (h, "SRA/%s_R1.fastq", sraID) ;
	  fNam2 = hprintf (h, "SRA/%s_R2.fastq", sraID) ;
	}
      else
	{
	  fNam  = hprintf (h, "SRA/%s_R1.fasta", sraID) ;
	  fNam2 = hprintf (h, "SRA/%s_R2.fasta", sraID) ;
	}
    }

  *ao1p = *ao2p = 0 ;
  if (fNam)
    {
      char *cr = filName (fNam, 0, "r") ;
      if (cr)
	fprintf (stderr, "Found cached file %s\n", fNam) ;
      else
	{
	  cr = filName (fNam, ".gz", "r") ;
	  if (cr)
	    fprintf (stderr, "Found cached file %s.gz\n", fNam) ;
	}
      if (!cr)
	{
	  *ao1p = aceOutCreate (fNam, 0, TRUE, h0) ;
	  if (!ao1p)
	    messcrash ("\nCannot create the SRA cache file %s", fNam) ;
	}
    }
  if (fNam && fNam2)
    {
      char *cr = filName (fNam2, 0, "r") ;
      if (cr)
	fprintf (stderr, "Found cached file %s\n", fNam2) ;
      else
	{
	  cr = filName (fNam2, ".gz", "r") ;
	  if (cr)
	    fprintf (stderr, "Found cached file %s.gz\n", fNam2) ;
	}
      if (! cr)
	{
	  *ao2p = aceOutCreate (fNam2, 0, TRUE, h0) ;
	  if (!*ao2p)
	    messcrash ("\nCannot create the SRA cache file %s", fNam2) ;
	}
      if (*ao2p && ! *ao1p)
	messcrash ("Only file %s exists, please remove it", fNam) ; 
    }
  ac_free (h) ;
  return ;
} /* sraCachingOut */

/**************************************************************/

static void sraCacheDo (ACEOUT ao1, ACEOUT ao2,	const char *seq, const char *seq2)
{
  if (seq && ao1) /* all data goes in same file */
    {
      aceOut (ao1, seq) ;
      if (seq2 && ao2) /* all data goes in same file */
	aceOut (ao2, seq2) ;
    }
} /* sraCacheDo */

/**************************************************************/

int saParseSraDownload (PP *pp, const char *sraID)
{
  AC_HANDLE h = ac_new_handle () ;
  char *fNam = 0 ;
  char *cr = 0 ; 
  ACEOUT ao1 = 0, ao2 = 0 ; 
  char tBuf[25] ;
  BOOL fastq = pp->fastq ;
  int split_pairs = pp->split_pairs ;

  if (!filCreateDir("./SRA"))
    messcrash ("\nCannot create or cannot write in the SRA cache directory ./SRA") ;

  {{  /* check in the cache */
      cr = 0 ;

      if (! fastq)  /* fastq contains fasta, no need to download again */
	{
	  if (! cr)
	    {
	      fNam = hprintf (h, "SRA/%s.fasta", sraID) ;
	      cr = filName (fNam, 0, "r") ;
	    }
	  if (! cr)
	    {
	      fNam = hprintf (h, "SRA/%s.sample_12.fasta", sraID) ;
	      cr = filName (fNam, 0, "r") ;
	    }
	  if (! cr)
	    {
	      fNam = hprintf (h, "SRA/%s_R2.fasta", sraID) ;
	      cr = filName (fNam, 0, "r") ;
	    }
	  if (! cr)
	    {
	      fNam = hprintf (h, "SRA/%s.fasta.gz", sraID) ;
	      cr = filName (fNam, 0, "r") ;
	    }
	  if (! cr)
	    {
	      fNam = hprintf (h, "SRA/%s_R2.fasta.gz", sraID) ;
	      cr = filName (fNam, 0, "r") ;
	    }
	  if (! cr)
	    {
	      fNam = hprintf (h, "SRA/%s.sample_12.fasta.gz", sraID) ;
	      cr = filName (fNam, 0, "r") ;
	    }
	}      

      else   /* fastq requested */
	{
	  if (! cr)
	    {
	      fNam = hprintf (h, "SRA/%s.fastq", sraID) ;
	      cr = filName (fNam, 0, "r") ;
	    }
	  if (! cr)
	    {
	      fNam = hprintf (h, "SRA/%s.sample_12.fastq", sraID) ;
	      cr = filName (fNam, 0, "r") ;
	    }
	  if (! cr)
	    {
	      fNam = hprintf (h, "SRA/%s_R2.fastq", sraID) ;
	      cr = filName (fNam, 0, "r") ;
	    }
	  if (! cr)
	    {
	      fNam = hprintf (h, "SRA/%s.fastq.gz", sraID) ;
	      cr = filName (fNam, 0, "r") ;
	    }
	  if (! cr)
	    {
	      fNam = hprintf (h, "SRA/%s_R2.fastq.gz", sraID) ;
	      cr = filName (fNam, 0, "r") ;
	    }
	  if (! cr)
	    {
	      fNam = hprintf (h, "SRA/%s.sample_12.fastq.gz", sraID) ;
	      cr = filName (fNam, 0, "r") ;
	    }
	}
      
      if (cr)
	{            /* file already in cache */
	  fprintf (stderr, "File %s already in cache\n", cr) ;	  
	  ac_free (h) ;
	  return 0 ;
	}
    }}
  /* download */
  
  SRAReadBatch* sra = SRAReadBatchNew (sraID);
  int nn = 0 ;
  BOOL firstPass = TRUE ;
  float Gb = pp->maxSraGb ;
  long int bMax = Gb * (1000000000L) ;  /* to be in decimal Gigabases */
  int BMAX = 1 << 23 ; /* 8 Mb */
    
  fprintf (stderr, "%s : SRA download %s ", timeBufShowNow(tBuf), sraID) ;
  if (Gb) fprintf (stderr, "(top %.3f GigaBases) ", Gb) ;
  
  while (Gb == 0 || bMax > 0)
    {
      int num_bases = BMAX ;
      if (! Gb && num_bases > bMax) num_bases = bMax ;
      

      SraGetReadBatch (sra, num_bases, fastq, split_pairs) ;
      if (sra->seq)
	{
	  if (Gb > 0) bMax -= sra->num_bases ;
	  nn++ ;
	  if (firstPass)
	    sraCachingOut (&ao1, &ao2, sra->seq, sra->seq2, sraID, fastq, sra->is_paired, split_pairs, h) ;
	  firstPass = FALSE ;
	  fprintf (stderr, ".") ;

	  if (!ao1 && ! ao2)
	    break ;
	  if (ao1 || ao2)
	    sraCacheDo (ao1, ao2, sra->seq, sra->seq2) ;
	}
      else
	break ;
    }
  SRAReadBatchFree(sra);
  if (ao1)
    fprintf (stderr, " %s downloaded %d data blocks infile %s\ndone: %s\n", sraID, nn, aceOutFileName (ao1), timeBufShowNow(tBuf)) ;
  
  ac_free (h) ;
  return 0 ;
} /* saSequenceParseSraDownload */

/********************************************************************************************************************/
/******************************** sequence parsing format choice ***************************************************/
/********************************************************************************************************************/

void saParse (const PP *pp, RC *rc) 
{
  DnaFormat format = rc->format ;

  if (format == SRA ||  pp->sraCaching)
    sraSequenceParser (pp, rc, 0, 0, 0) ;
  else if (!rc->fileName2)
    saParseFastac (pp, rc) ;
  else
    otherSequenceParser (pp, rc, 0, 0, 0) ;
  return ;
} /* saSequenceParse */

/********************************************************************************************************************/
/******************************** sequence file parsing utilities ***************************************************/
/********************************************************************************************************************/


