/*
 * sa.sequenceParse.c

 * This module is part of the sortalign package
 * Created November, 2025

 * This code is public.

 * This module downloads runs from the SRA sequence archives
 * Optionally cache them in local directory .//SRA
 * Pases the sequences which can be in format fasta/fastq/fastc or SRA
 */

#include "sa.h"

#define NAMMAX 1024

/* download SRA runs as FASTA or FASTQ */
#define SRA_FASTA 0
#define SRA_FASTQ 1

/* static atomic_int lane = 0 ; */

/**************************************************************/
/**************************************************************/
/* create the GLOBAL DNA
 * Copy all DNAS of the BB block in a single bigArray globalDna
 * Set the individual dna array to point into the globalDna
 * Collect their coordinates 
 */
static void globalDnaCreate (BB *bb)
{
  long int ln = 0 ;
  int n, ii, iMax = arrayMax (bb->dnas) ;
  unsigned char *cp, *cq ;
  Array dna = 0 ;

  if (! bb->isGenome && (iMax & 0x1)) /* complete the last read pair with a zero */
    {
      array (bb->dnas, iMax, Array) = 0 ;
      iMax++ ;
    }
  
  bb->dnaCoords = bigArrayHandleCreate (2 * (iMax + 1), unsigned int, bb->h) ;

  /* compute the length of the global DNA (avoid reallocations) */
  for (ln = 0, ii = 1 ; ii < iMax ; ii++)
    {
      dna = array (bb->dnas, ii, Array) ;
      if (dna)
	{
	  bigArray (bb->dnaCoords, 2*ii, unsigned int) = ln ;    /* off set of dna ii */
	  n = arrayMax (dna) ;
	  bigArray (bb->dnaCoords, 2*ii + 1, unsigned int) = ln + n ;    /* off set of dna ii */
	  ln += n + 16 ;            /* mininal 16 zeroes protection */
	  n = n % 16 ;
	  if (n)          /* align the data */
	    ln += 16 - n ;
	}
    }
  bigArray (bb->dnaCoords, 2 * ii, unsigned int) = ln ;    /* global end */
  bigArray (bb->dnaCoords, 2 * ii + 1, unsigned int) = ln ;    /* global end */ 
  
  /* construct the global DNA  */
  if (!ln) messcrash ("Could not read any sequence in the fasta file") ;
  bb->globalDna = bigArrayHandleCreate (ln, unsigned char, bb->h) ;
  cp = bigArrayp (bb->globalDna, ln -1, unsigned char) ; /* make room */
  for (ln = 0, ii = 1 ; ii < iMax ; ii++)
    {
      /* transfer the sequence to the globalDna array */
      dna = array (bb->dnas, ii, Array) ;
      if (dna)
	{
	  n = arrayMax (dna) ;            /* mininal 16 n protection */
	  cp = bigArrp (bb->globalDna, ln, unsigned char) ;
	  cq = arrp (dna, 0, unsigned char) ;
	  memcpy (cp, cq, n) ; 
	  messfree (dna->base) ;
	  arrayLock (dna) ;
	  dna->base = (char *) cp ; 
	  ln += n ; cp += n ;
	  /* protect each sequence with a terminal A to allow constructing w1 in codeWords */
	  memset (cp, 0, 16) ; /* *cp = A_ */ ; ln += 16 ; cp += 16 ;
	  n = n % 16 ;
	  if (n)          /* align the data */
	    {  memset (cp, 0, 16 - n) ; ln += 16 - n ; cp += 16 - n ; }
	}
    }

  return ;
} /* globalDnaCreate */
  
/**************************************************************/

static BOOL parseOneSequence (DnaFormat format, char *namBuf, ACEIN ai, Array dna, Array qual, int *linep)
{
  int n, nn ;
  char *cp ;
  
  arrayMax (dna) = 0 ;
  if (qual) arrayMax (qual) = 0 ;
  memset (namBuf, 0, NAMMAX) ;

  if (aceInCard (ai))
    {
      (*linep)++ ;
      cp = aceInWord (ai) ;
      switch (format)
	{
	case RAW:
	  if (!cp) return FALSE ;
	  n = strlen ((char *)cp) ;
	  if (!n) return FALSE ;
 	  sprintf (namBuf, "s.%d\n", *linep) ;
	  array (dna, n+1, unsigned char) = 0 ; /* make room */
	  memcpy (arrp (dna, 0, char), cp, n) ;
	  arrayMax (dna) = n ;
	  break ; /* cp already points to the DNA sequence */
	case FASTA:
	case FASTA2:
	  while (! cp || *cp == '#')
	    {
	      cp = 0 ;
	      if (aceInCard (ai))
		{
		  (*linep)++ ;
		  cp = aceInWord (ai) ;
		}
	      else
		return FALSE ;
	    }
	  if (! cp) /* no > identifier in the ramainder of the fasta file */
	    return FALSE ;
	  if (cp[0] != '>' || cp[1] == 0)
	    messcrash ("\nMissing identifier at line %d of fasta sequence file %s\n", *linep, aceInFileName (ai)) ; 
	  strncpy (namBuf, cp + 1, NAMMAX - 2) ;
	  n = nn = 0 ;
	  while (aceInCard (ai))
	    { /* parse the sequence, it may spread over many lines */
	      cp = aceInWord (ai) ;
	      if (!cp)
		break ;
	      if (*cp == '>')
		{
		  aceInCardBack (ai) ; /* fold back on the indentifier of the previous sequence in the ACEIN data flow */
		  break ;
		}
	      (*linep)++ ;
	      n = strlen ((char *)cp) ;
	      if (!n && !nn) return FALSE ;
	      if (1) { char *cq = strstr (cp, "><") ; if (cq) {*cq = 0 ; n = cq - cp ; }}
	      array (dna, nn + n + 1, unsigned char) = 0 ; /* make room */
	      memcpy (arrp (dna, nn, char), cp, n) ;
	      nn += n ;
	      arrayMax (dna) = nn ;
	      }
	  if (! nn)
	    return FALSE ;
	  break ;
	case FASTQ:
	  if (!cp) return FALSE ;
	  if (cp[0] != '@')
	    messcrash ("\nMissing @ identifier, found %s\nat line %d of fastq sequence file %s\n", cp ? cp : "NULL", *linep, aceInFileName (ai)) ;
	  strncpy (namBuf, cp + 1, NAMMAX - 2) ;
	  n = 0 ;
	  if (aceInCard (ai))
	    { /* parse the sequence, it resides on a single line */
	      (*linep)++ ;
	      cp = aceInWord (ai) ;
	      n = cp ? strlen (cp) : 0 ;
	      if (n)
		{
		  array (dna, n + 1, unsigned char) = 0 ; /* make room */
		  memcpy (arrp (dna, 0, char), cp, n) ;
		  arrayMax (dna) = n ;
		}
	    }
	  if (! n)
	    messcrash ("\nMissing sequence at line %d of fastq sequence file %s\n", *linep, aceInFileName (ai)) ;
	  if (aceInCard (ai))
	    (*linep)++ ;
	  else
	    messcrash ("\nMissing quality identifier at line %d of fastq sequence file %s\n", *linep, aceInFileName (ai)) ;
	  if (aceInCard (ai))
	    {

	      (*linep)++ ;
	      if (qual) 
		{
		  cp = aceInWord (ai) ;
		  n = cp ? strlen (cp) : 0 ;
		  if (n)
		    {
		      array (qual, n + 1, unsigned char) = 0 ; /* make room */
		      memcpy (arrp (qual, 0, char), cp, n) ;
		      arrayMax (qual) = n ;
		    }
		}
	    }
	  else
	    messcrash ("\nMissing quality factors at line %d of fastq sequence file %s\n", *linep, aceInFileName (ai)) ;
	  break ;
	case FASTC:
	  while (! cp || *cp == '#')
	    {
	      cp = 0 ;
	      if (aceInCard (ai))
		{
		  (*linep)++ ;
		  cp = aceInWord (ai) ;
		}		  
	      else
		return FALSE ;
	    }
	  if (! cp) /* no > identifier in the fasta file */
	    return FALSE ;
	  if (cp[0] != '>' || cp[1] == 0)
	    messcrash ("\nMissing identifier at line %d of fastc sequence file %s\n", *linep, aceInFileName (ai)) ; 
	  strncpy (namBuf, cp + 1, NAMMAX - 2) ;
	  n = 0 ;
	  if (aceInCard (ai))
	    { /* parse the sequence, it resides on a single line */
	      (*linep)++ ;
	      cp = aceInWord (ai) ;
	      n = strlen (cp) ;
	      if (n)
		{
		  array (dna, n + 1, unsigned char) = 0 ; /* make room */
		  memcpy (arrp (dna, 0, char), cp, n) ;
		  arrayMax (dna) = n ;
		}
	    }
	  if (! n)
	    messcrash ("\nMissing sequence at line %d of fastc sequence file %s\n", *linep, aceInFileName (ai)) ;
	  break ;
	case SRA:
	  messcrash ("\nParsing SRA is not yet implemented, sorry") ;
	  break ;
	default:
	    messcrash ("\nBad format passed to sequenceParser, please edit the sdource code, sorry") ;
	}
      
      /*
	switch (DnaFormat)
	{
	case RAW:
	  break ;
	case FASTA:
	  break ;
	case FASTQ:
	  break ;
	case FASTC:
	  break ;
	case SRA; 
	  break ;
	}
      */
    }
  return arrayMax (dna) ? TRUE : FALSE ;
} /* parseOneSequence */

/**************************************************************/

static BOOL parseOnePair (DnaFormat format, char *namBuf
			  , ACEIN ai1, Array dna1, Array qual1, int *linep1
			  , ACEIN ai2, Array dna2, Array qual2, int *linep2
			  )
{
  BOOL ok1 = FALSE, ok2 = TRUE ;
  int n ;
  int line1 = *linep1 ;
  int line2 = linep2 ? *linep2 : 0 ;
  
  arrayMax (dna1) = 0 ;
  if (dna2) arrayMax (dna2) = 0 ;
  if (qual1) arrayMax (qual1) = 0 ;
  if (qual2) arrayMax (qual2) = 0 ;
    
  memset (namBuf, 0, NAMMAX) ;
  if (ai1)
    ok1 = parseOneSequence (format, namBuf, ai1, dna1, qual1, linep1) ;
  if (ok1 && ai2)
    {
      char namBuf2 [NAMMAX] ;
      int ln = strlen (namBuf) ;
      ok2 = parseOneSequence (format, namBuf2, ai2, dna2, qual2, linep2) ;
      if (strncmp (namBuf, namBuf2, ln-1))
	messcrash ("\nNon matching pair of sequence identifiers %s <> %s\n line %d of file %s\nlne %d of file %s\n"
		   , namBuf, namBuf2, *linep1, aceInFileName (ai1), *linep2, aceInFileName (ai2)) ;
    }

  if (format == FASTC && ok1)
    {    /* fastc format specifies the pair on a single line with >< separator */
      unsigned char *cp = arrp (dna1, 0, unsigned char) ;
      unsigned char *cq = cp ? (unsigned char *) strstr ((char *)cp, "><") : 0 ;
      if (cq)
	{
	  n = strlen ((char *)cq +2) ;
	  if (n)
	    {
	      array (dna2, n + 2, unsigned char) = 0 ; /* make room */
	      memcpy (arrp (dna2, 0, unsigned char), cq + 2, n) ;
	      arrayMax (dna2) = n ;
	    }
	  *cq = 0 ;
	  n = cq - cp ;
	  arrayMax (dna1) = n ;
	}      
    }
  
  n = arrayMax (dna1) ;
  if (n)
    {
      dnaEncodeArray (dna1) ;
      if (0 &&  /* - (not sequenced) i encoded as 0, i do not dare to change that in w1/dnacode.c  */
	  n != strlen (arrp (dna1, 0, char)))  /* will detect an early  zero */
	messcrash ("\nunrecognized DNA character in sequence %s line %d of file %s\n"
		   , namBuf, line1, aceInFileName (ai1)) ;
      
    }
  n = dna2 ? arrayMax (dna2) : 0 ;
  if (n)
    {
      dnaEncodeArray (dna2) ;
      if (n != strlen (arrp (dna2, 0, char)))  /* will detect an early  zero */
	{
	  if (ai2)
	    messcrash ("\nunrecognized DNA character in sequence %s line %d of file %s\n"
		       , namBuf, line2, aceInFileName (ai2)) ;
	  else  /* may happen in fastc case */
	    messcrash ("\nunrecognized DNA character in sequence %s line %d of file %s\n"
		       , namBuf, line1, aceInFileName (ai1)) ;
	}
    }

  return ok1 && ok2 ;
} /* parseOnepair */

/**************************************************************/
/* parse a fasta buffer into an array of DNA */ 
void saSequenceParseGzBuffer (const PP *pp, BB *bb)
{
  if (bb->gzBuffer) /* fasta buffer */
    {
      AC_HANDLE h = ac_new_handle () ;
      Array qual1 = 0, dna1 = arrayHandleCreate (256, unsigned char, bb->h) ;
      int line1 = 0 ; 
      char namBuf [NAMMAX + 12] ;
      ACEIN ai = aceInCreateFromText (bb->gzBuffer, 0, h) ;
      DnaFormat format = bb->rc.format ;

      unsigned char atgcn[256] ;
      memset (atgcn, 4, sizeof(atgcn)) ;
      atgcn[A_] = 0 ;
      atgcn[T_] = 1 ;
      atgcn[G_] = 2 ;
      atgcn[C_] = 3 ;
      
      memset (namBuf, 0, NAMMAX) ;
      bb->errors = arrayHandleCreate (256, int, bb->h) ;
      bb->txt1 = vtxtHandleCreate (bb->h) ;
      bb->txt2 = vtxtHandleCreate (bb->h) ;
      bb->length = 0 ;
      bb->dnas = arrayHandleCreate (bb->nSeqs, BigArray, bb->h) ;
      bb->dict = dictHandleCreate (bb->nSeqs, bb->h) ;
      bb->runStat.p.lengthDistribution = arrayHandleCreate (1024, long int, bb->h) ;
      bb->runStat.insertLengthDistribution = arrayHandleCreate (1024, long int, bb->h) ;
      bb->nSeqs = 0 ;
      bb->errDict = dictHandleCreate (100000, bb->h) ;
      
      while (parseOnePair (format, namBuf, ai, dna1, qual1, &line1, 0, 0, 0, 0)) 
	{
	  int nn1, n1 = arrayMax (dna1) ;
	  int isRead2 = 0 ;
	  char *cp ;
	  
	  bb->nSeqs++ ;
	  bb->length += n1 ;
	  bb->runStat.p.nReads++ ;

	  switch (bb->rc.format)
	    {
	    case FASTA2:
	      cp = namBuf + strlen (namBuf) - 2 ;
	      if (! strcmp (cp, ".2"))
		{
		  isRead2 = 1 ;
		  *cp = 0 ;
		  bb->runStat.p.nPairs++ ;
		}
	      if (! strcmp (cp, ".1"))
		{
		  isRead2 = 0 ;
		  *cp = 0 ;
		}
	      break ;
	    default:
	      break ;
	    }
	  
	  if (isRead2) bb->runStat.p.nBase2 += n1 ;
	  else  bb->runStat.p.nBase1 += n1 ;

	  dictAdd (bb->dict, namBuf, &nn1) ;
	  nn1 = (nn1 << 1) | isRead2 ;
	  array (bb->dnas, nn1, Array) = dna1 ;
	  
	  if (arrayMax (dna1))
	    {
	      int i, iMax = arrayMax (dna1) ;
	      unsigned char *cp = arrp (dna1, 0, unsigned char) ;
	      for (i = 0 ; i < iMax ; i++, cp++)
		bb->runStat.p.ATGCN[atgcn[(int)*cp]]++ ;
	      
	      if (iMax > bb->runStat.p.maxReadLength)
		bb->runStat.p.maxReadLength = iMax ;
	      if (! bb->runStat.p.minReadLength || iMax < bb->runStat.p.minReadLength)
		bb->runStat.p.minReadLength = iMax ;
	      array (bb->runStat.p.lengthDistribution, iMax, long int)++ ;
	    }
	  if (arrayMax (dna1))
	    {
	      int i, iMax = arrayMax (dna1) ;
	      unsigned char *cp = arrp (dna1, 0, unsigned char) ;
	      for (i = 0 ; i < iMax && i < LETTERMAX ; i++, cp++)
		bb->runStat.p.letterProfile1[5*i + atgcn[(int)*cp]]++ ;
	    }

	  dna1 = arrayHandleCreate (256, unsigned char, bb->h) ;
	}
      ac_free (h) ;
      globalDnaCreate (bb) ;
    }
  
  return ;
} /* saSequenceParseGzBuffer */

/**************************************************************/

static int dnaSequenceOrder (const void *va, const void *vb)
{
  const Array *up = va ;
  const Array *vp = vb ;
  const char *cp = *up ? arrp (*up, 0, char) : 0 ;
  const char *cq = *vp ? arrp (*vp, 0, char) : 0 ;

  if (!cp)
    return cq ? 1 : 0 ;
  if (!cq)
    return cp ? -1 : 0 ;
  return strcmp (cp, cq) ;
} /* dnaSequenceOrder */

/**************************************************************/
/* add the multiplicities in the # filed, creating a fastc format */
void saSequenceDeduplicate (const PP *pp, BB *bb)
{
  AC_HANDLE h = ac_new_handle () ;
  Array dnas = bb->dnas ;
  int i, j, iMax = arrayMax (dnas), k, kk = 1 ;
  DICT *newDict = dictHandleCreate (iMax, bb->h) ;
  Array newDnas = arrayHandleCreate (iMax, Array, bb->h) ;
  char buf[64] = {0} ;
  ACEOUT ao = aceOutCreate (pp->outFileName, ".fastc", pp->gzo, h) ;
  
  arraySort (dnas, dnaSequenceOrder) ;
  for (i = 0 ; i < iMax ; i++)
    {
      int mult = 1 ;
      Array dna1 = arr (dnas, i, Array) ;
      if (dna1)
	{
	  const char *cp = arrayMax (dna1) ? arrp (dna1, 0, char) : 0 ;
	  if (! cp)
	    continue ;
	  for (j = i + 1 ; j < iMax ; j++)
	    {
	      Array dna2 = arr (dnas, j, Array) ;
	      if (dna2)
		{
		  const char *cq = arrp (dna2, 0, char) ;
		  if (strcmp (cp, cq))
		    break ;
		  mult++ ;
		}
	      else
		break ;
	    }

	  sprintf (buf, "s%d#%d", kk++, mult) ;
	  dictAdd (newDict, buf, &k) ;
	  array (newDnas, k, Array) = dna1 ;
	  i += mult - 1 ;
	  aceOutf (ao, ">%s\n%s\n", buf, cp) ;
	}
    }
  ac_free (h) ;
  exit (0) ;
}

/**************************************************************/
/* aug 2
 *  f = gzopen ("gilname.gz", r)
 *  cp = gzread (f, 100Mega, buffer) ;
 *  cq = strrchr (cp, '>')
 *    if (!cq || cq == cp) continue reading intil EOF
 *    else { *(cq-1)=0; pass the cp buffer to a new agent which will decode the dna ;
 *    copy cq ... (n bytes)  to a new clean buffer and gzread in buffer+n
 *  gzclose() 
 * the general idea is that parsing big buffers is fast, while decoding them is slow
 * so in this way, even when facing a single large fastq, the pipeline will no longer be hanged on the parser
 */

static void fastaSequenceParser (const PP *pp, RC *rc, TC *tc, BB *bb, int isGenome)
{
  AC_HANDLE h = ac_new_handle () ;
  BB b ;
  int BMAX = isGenome ? 100000 : (pp->BMAX << 20) ;
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
  const char *fileName1 = rc ? rc->fileName1 : tc->fileName ;
  char tBuf[25] ;
  clock_t t1, t2 ;

  t1 = clock () ;
  
  file = gzopen (fileName1, "r");
  if (! file)
    messcrash ("\ncannot gzopen target file %s", fileName1) ;

  while (!done)
    {
      int err = 0 ;                    
      int bytes = gzread (file, buffer + pos, BMAX - pos) ;
      unsigned char *restrict cp ;

      bytes += pos ;


      if (bytes && format == FASTA2)
	{   /* check for identifiers signalling a paired end read */
	  unsigned char *cq = buffer ;
	  int nDots = 0, k = 0 ;
	  while (k++ < bytes && cq && *cq != '\n')
	    nDots += (*cq++ == '.' ? 1 : 0) ;
	  if (nDots == 3)
	    {
	      format = FASTA2 ;
	      pairedEnd = TRUE ;
	    }
	  else
	    {
	      format = FASTA ;
	      pairedEnd = FALSE ;
	    }
	}
      
      if (bytes < BMAX)
	{
	  done = TRUE ;
	  if (! gzeof (file)) 
	    messcrash ("Error %s in gzread %s", gzerror (file, & err), fileName1) ;                
        }
      else
	{
	  /* search for beginning of last probably partial sequence */
	  cp = buffer + bytes - 1 ; /* last byte read */
	  pos = 0 ;
	  while (cp > buffer && *cp != '>')
	    { pos++ ; cp-- ;}
	  if (*cp == '>' && format == FASTA2 && cp > buffer &&  ! strstr ((char *)cp, ".1\n"))
	    {  /* cp is the second read of the pair, search for the previous read */
	      cp-- ; pos++ ;
	      while (cp > buffer && *cp != '>')
		{ pos++ ; cp-- ;}
	    }
	  if (*cp != '>' || (cp > buffer && cp[-1] != '\n'))
	    messcrash ("gzread found a read > BMAX=%d", BMAX) ;
	  pos++ ;

	  memcpy (buffer2, cp, pos) ; /* copy the remnant */
	  bytes -= pos ;
	  cp[0] = 0 ;
	}
      nBytes += bytes ;
      if (!nBytes)
	messcrash ("No sequence found in file %s\n", fileName1) ;
      

      /* create a data block */
      bb = &b ;
      memset (bb, 0, sizeof (BB)) ;
      bb->h = ac_new_handle () ;
	  
      bb->rc.format = format ;
      bb->rc.pairedEnd = pairedEnd ;
      bb->rc.jump5r1 = rc ? rc->jump5r1 : 0 ;
      bb->rc.jump5r2 = rc ? rc->jump5r2 : 0 ;      

      bb->readerAgent = pp->agent ;
      bb->run = rc ? rc->run : 0 ;
      bb->start = timeNow () ;
      /*
	bb->lane = atomic_fetch_add (rc ? &(rc->lane) : &lane, 1) + 1 ;
      */

      bb->lane = atomic_fetch_add (arrp (pp->runLanes, bb->run, atomic_int), 1) + 1 ;

      bb->cpuStats = arrayHandleCreate (128, CpuSTAT, bb->h) ;
      bb->rc.fileName1 = fileName1 ;
      /* copy the buffer */
      bb->gzBuffer = halloc (bytes + 1, bb->h) ;
      memcpy (bb->gzBuffer, buffer, bytes) ;
      bb->gzBuffer[bytes] = 0 ;
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
} /* fastaSequenceParser */

/**************************************************************/
/**************************************************************/
/* we need to download 4 lines to recognize pairs */
static BOOL sraIsPairedEnd (const char *seq, BOOL *hasQualp)
{
  BOOL isPaired = FALSE ;
  const char *cp = seq ;

  cp = seq ;
  for (int i = 0 ; i < 2; i++)
    if (cp) cp = strchr(cp+1, '\n') ;
  *hasQualp = (cp && cp[1] == '+' ? TRUE : FALSE) ;

  cp = seq ;
  for (int i = 0 ; i < 2 ; i++)
    if (cp) cp = strchr(cp+1, '.') ;
  if (cp) cp = strchr(cp+1, '/') ;
  if (cp && cp[1] == '1')
    {
      cp = seq ;
      for (int i = 0 ; i < (*hasQualp ? 4 : 2) ; i++)
	if (cp) cp = strchr(cp+1, '\n') ;
      for (int i = 0 ; i < 2 ; i++)
	if (cp) cp = strchr(cp+1, '.') ;
      if (cp) cp = strchr(cp+1, '/') ;
      if (cp && cp[1] == '2')
	return TRUE ;
    }      
  return isPaired ;
} /* sraIsPairedEnd */

/**************************************************************/
/* we need to download 4 lines to recognize pairs */
static void sraCachingOut (ACEOUT *ao1p, ACEOUT *ao2p, const char *seq, const char *seq2, const char *sraID
			   , BOOL fastq, BOOL split_pairs, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  char *fNam = 0, *fNam2 = 0 ;
  BOOL hasQual = fastq ;
  BOOL isPaired = seq2 ? TRUE : sraIsPairedEnd (seq, &hasQual) ;

  if (! isPaired)
    split_pairs = FALSE ;
  if (! hasQual)
    fastq = FALSE ; 
  else
    fastq = TRUE ; 

  if (! isPaired)
    {
      if (fastq)
	fNam = hprintf (h, "SRA/%s.fastq", sraID) ;
      else
	fNam = hprintf (h, "SRA/%s.fasta", sraID) ;
    }
  else if (!split_pairs)
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
/* nov 2 2025
 * grab data directly from SRA sequence archives
 * f = gzopen ("gilname.gz", r)
 *  cp = gzread (f, 100Mega, buffer) ;
 *  cq = strrchr (cp, '>')
 *    if (!cq || cq == cp) continue reading intil EOF
 *    else { *(cq-1)=0; pass the cp buffer to a new agent which will decode the dna ;
 *    copy cq ... (n bytes)  to a new clean buffer and gzread in buffer+n
 *  gzclose() 
 * the general idea is that parsing big buffers is fast, while decoding them is slow
 * so in this way, even facinga single large fastq, the pipeline will no longer be hanged on the parser
 */

static void sraSequenceParser (const PP *pp, RC *rc, TC *tc, BB *bb, int isGenome)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao1 = 0, ao2 = 0 ;
  BB b ;
  CHAN *chan = pp->plChan ;
  BOOL debug = FALSE ;
  int BMAX = isGenome ? 100000 : (pp->BMAX << 20) ;
  long int bytes = 0, bytes2 = 0,  nBytes = 0 ;
  int nPuts = 0 ;
  DnaFormat format = rc->format ;
  const char *sraID = rc ? rc->fileName1 : tc->fileName ;
  char tBuf[25] ;
  clock_t t1, t2 ;

  if (isGenome || (!pp->sraCaching && format != SRA && format != FASTA2 && format != FASTQ2))
    messcrash ("Bad internal options in sraSequenceParser, please edit the code, sorry") ;
  
  if (format == SRA || pp->sraCaching)  /* check in the cache */
    {
      char *fNam = 0 ;
      fNam = hprintf (pp->h, "SRA/%s", sraID) ;
      char *cr = 0 ;
      
      if (! cr)
	{
	  cr = filName (fNam, ".fasta", "r") ;
	  if (cr)
	    {
	      format = FASTA ;
	      fprintf (stderr, "Found cached file %s\n", fNam) ;
	    }
	}
      if (! cr)
	{
	  cr = filName (fNam, "fasta.gz", "r") ;
	  if (cr)
	    {
	      format = FASTA ;
	      fprintf (stderr, "Found cached file %s.gz\n", fNam) ;
	    }
	}
      if (! cr)
	{
	  cr = filName (fNam, ".sample_12.fasta.gz", "r") ;
	  if (cr)
	    {
	      format = FASTA2 ;
	      fprintf (stderr, "Found cached file %s.sample_12.fasta.gz\n", fNam) ;
	    }
	}
      if (! cr)
	{
	  cr = filName (fNam, ".sample_12.fasta", "r") ;
	  if (cr)
	    {
	      format = FASTA2 ;
	      fprintf (stderr, "Found cached file %s.sample_12.fasta\n", fNam) ;
	    }
	}
      if (! cr)
	{
	  cr = filName (fNam, ".sample_12.fastq.gz", "r") ;
	  if (cr)
	    {
	      format = FASTQ2 ;
	      fprintf (stderr, "Found cached file %s.sample_12.fastq.gz\n", fNam) ;
	    }
	}
      
      if (! cr)
	{
	  cr = filName (fNam, ".sample_12.fastq", "r") ;
	  if (cr)
	    {
	      format = FASTQ2 ;
	      fprintf (stderr, "Found cached file %s.sample_12.fastq\n", fNam) ;
	    }
	}
      if (! cr)
	{
	  cr = filName (fNam, ".fastq.gz", "r") ;
	  if (cr)
	    {
	      format = FASTQ ;
	      fprintf (stderr, "Found cached file %s.fastq.gz\n", fNam) ;
	    }
	}
      if (! cr)
	{
	  cr = filName (fNam, ".fastq", "r") ;
	  if (cr)
	    {
	      format = FASTQ ;
	      fprintf (stderr, "Found cached file %s.fastq\n", fNam) ;
	    }
	}
      
      if (cr)
	{
	  rc->fileName1 = strnew (cr, pp->h) ;
	  rc->format = format ;
	  ac_free (h) ;
	  return fastaSequenceParser (pp, rc, tc, bb, isGenome) ; 
	}
    }
  
  t1 = clock () ;
  SRAObj* sra = SraObjNew(sraID);
  const char *seq, *seq2 ;
  BOOL firstPass = TRUE ;
  BOOL fastq = pp->fastq ;
  BOOL split_pairs = pp->split_pairs ;
  int num_bases = 1 << 27 ; /* 128 M */
  float Gb = pp->maxSraGb ;
  long unsigned int nMax = Gb * (1000000000L) ;  /* to be in decimal Gigabases */

  format = FASTA ;
  if (! pp->sraCaching) /* for internal use, prefer interleaved fasta */
    fastq = FALSE ;
  split_pairs = FALSE ;

  if (num_bases > nMax) num_bases = nMax ;
  nMax =  (nMax + num_bases - 1) / num_bases ; 
  format = SRA ;

  if (! Gb) num_bases = BMAX ;
  while (!Gb || nMax-- > 0)
    {
      SraGetReadBatch(sra, num_bases, fastq, split_pairs, &seq, &seq2) ;
      if (seq)
	{
	  if (firstPass && pp->sraCaching)
	  sraCachingOut (&ao1, &ao2, seq, seq2, sraID, fastq, split_pairs, h) ;
	  firstPass = FALSE ;
	  if (ao1 || ao2)
	    sraCacheDo (ao1, ao2, seq, seq2) ;
	  
	  bytes = strlen (seq) ;
	  bytes2 = seq2 ? strlen (seq2) : 0 ;
	  nBytes += bytes + bytes2 ; 
	  if (!bytes)
	    messcrash ("No sequence found in SRA %s\n", sraID) ;
	  
	  if (format == SRA)
	    {   /* check for identifiers signaling a paired end read */
	      if (split_pairs)
		format = FASTA ;
	      else
		format = FASTA2 ;
	    }
	  
	  /* create a data block */
	  bb = &b ;
	  memset (bb, 0, sizeof (BB)) ;
	  bb->h = ac_new_handle () ;
	  bb->rc.format = format ;
	  bb->rc.jump5r1 = rc ? rc->jump5r1 : 0 ;
	  bb->rc.jump5r2 = rc ? rc->jump5r2 : 0 ;      
	  bb->readerAgent = pp->agent ;
	  bb->run = rc ? rc->run : 0 ;
	  bb->start = timeNow () ;
	  bb->lane = atomic_fetch_add (arrp (pp->runLanes, bb->run, atomic_int), 1) + 1 ;
	  bb->cpuStats = arrayHandleCreate (128, CpuSTAT, bb->h) ;
	  bb->rc.fileName1 = sraID ;
	  /* copy the buffer */
	  bb->gzBuffer = halloc (bytes + bytes2 + 1, bb->h) ;
	  memcpy (bb->gzBuffer, seq, bytes) ;
	  memcpy (bb->gzBuffer + bytes, seq2, bytes2) ;
	  
	  bb->gzBuffer[bytes + bytes2] = 0 ;
	  bb->nSeqs = 100 ;  /* a guess */
	  
	  /* export the databalock to the channel */
	  nPuts++ ;
	  channelPut (chan, bb, BB) ;
	}
    }
  channelPut (pp->npChan, &nPuts, int) ; /* global counting of BB blocks accross all sequenceParser agents */
  
  ac_free (h) ;
  
  t2 = clock () ;
  saCpuStatRegister ("2.sraSequenceParser", pp->agent, bb->cpuStats, t1, t2, nBytes) ;

  if (debug)
    {
     int lane = atomic_fetch_add (arrp (pp->runLanes, bb->run, atomic_int), 0) ;
     printf ("--- %s: Stop sraSequenceParser %d blocks %ld bytes file %s\n", timeBufShowNow (tBuf), lane, nBytes, sraID) ;
    }
  
  return ;
} /* sraSequenceParser */

/**************************************************************/

static void otherSequenceParser (const PP *pp, RC *rc, TC *tc, BB *bb, int isGenome)
{
  AC_HANDLE h = ac_new_handle () ;
  BB b ;
  int BMAX = isGenome ? 100000 : (pp->BMAX << 20 ) ;
  int NMAX = BMAX / 200 ;
  int nPuts = 0 ;
  int nn = 0, nSeqs = 0 ;
  CHAN *chan = 0 ;
  ACEIN ai1 = 0 ;
  ACEIN ai2 = 0 ;
  int line1 = 0, line2 = 0, line10 = 0 ;
  Array dna1 = 0, dna2 = 0, dnas = 0 ;
  Array qual1 = 0, qual2 = 0 ;
  char targetClass = tc ? tc->targetClass : 0 ;
  char namBufG [NAMMAX+12] ;
  char *namBufX, *namBuf = namBufG + 2 ;
  DnaFormat format = rc ? rc->format : tc->format ;
  const char *fileName1 = rc ? rc->fileName1 : tc->fileName ;
  const char *fileName2 = rc ? rc->fileName2 : 0 ;
  BOOL pairedEnd = rc ? rc->pairedEnd : FALSE ;
  char tBuf[25] ;
  clock_t t1, t2 ;
  
  unsigned char atgcn[256] ;
  memset (atgcn, 4, sizeof(atgcn)) ;
  atgcn[A_] = 0 ;
  atgcn[T_] = 1 ;
  atgcn[G_] = 2 ;
  atgcn[C_] = 3 ;
  
  t1 = clock () ;
  
  ai1 = aceInCreate (fileName1, 0, h) ;
  if (!ai1)
    messcrash ("\ncannot read target file %s", fileName1) ;
  aceInSpecial (ai1, "\n") ;
  if (fileName2)
    {
      ai2 = aceInCreate (fileName2, 0, h) ;
      if (!ai2)
	messcrash ("\ncannot read target file %s", fileName2) ;
      aceInSpecial (ai2, "\n") ;
    }
  
  if (isGenome)
    {
      chan = 0 ;
      namBufG[0] = tc->targetClass ;
      namBufG[1] = '.' ;
      namBufX = namBufG ;
    }
  else
    {
      bb = &b ;
      memset (bb, 0, sizeof (BB)) ;
      bb->readerAgent = pp->agent ;
      bb->start = timeNow () ;
      bb->lane = atomic_fetch_add (arrp (pp->runLanes, bb->run, atomic_int), 1) + 1 ;
      chan = pp->plChan ;
      namBufX = namBuf ;
    }
  
  if (! bb->h)
    {  /* isGenome: keep expanding the same bbG if already initialised */
      DnaFormat format = rc ? rc->format : tc->format ;
      
      bb->h = ac_new_handle () ;
      bb->txt1 = vtxtHandleCreate (bb->h) ;
      bb->txt2 = vtxtHandleCreate (bb->h) ;
      bb->errors = arrayHandleCreate (256, int, bb->h) ;
      bb->cpuStats = arrayHandleCreate (128, CpuSTAT, bb->h) ;
      bb->rc.jump5r1 = rc ? rc->jump5r1 : 0 ;
      bb->rc.jump5r2 = rc ? rc->jump5r2 : 0 ;      
      bb->run = rc ? rc->run : 0 ;
      bb->length = 0 ;
      bb->dnas = dnas = arrayHandleCreate (64, BigArray, bb->h) ;
      bb->runStat.p.lengthDistribution = arrayHandleCreate (1024, long int, bb->h) ;
      bb->runStat.insertLengthDistribution = arrayHandleCreate (1024, long int, bb->h) ;
      if (pp->exportSamQuality && format == FASTQ)
	bb->quals = arrayHandleCreate (64, Array, bb->h) ;
      bb->dict = dictHandleCreate (NMAX, bb->h) ;
      bb->errDict = dictHandleCreate (NMAX, bb->h) ;
      bb->rc.fileName1 = fileName1 ;
    }
  if (pp->debug || isGenome)
    printf ("+++ %s: Start sequence parser %s\n", timeBufShowNow (tBuf), fileName1) ;
  
  line1 = line2 = 0 ; line10 = 1 ;
  dna1 = arrayHandleCreate (isGenome ? 1 << 20 : 256, unsigned char, bb->h) ;
  dna2 = arrayHandleCreate (256, unsigned char, bb->h) ;
  if (pp->exportSamQuality && format == FASTQ)
    {
      qual1 = arrayHandleCreate (isGenome ? 1 << 20 : 256, unsigned char, bb->h) ;
      qual2 = arrayHandleCreate (256, unsigned char, bb->h) ;
    }
  while (parseOnePair (format, namBuf, ai1, dna1, qual1, &line1, ai2, dna2, qual2, &line2)) 
    {
      int mult = 1 ;
      char *cr = namBuf + strlen (namBuf) ;
      if (isGenome && targetClass == 'T')
	{
	  /* only keep transcrips called *.a   ||   .a|GENE=...|...  */
	  char *cp = strchr (namBuf, '|') ;
	  cp = cp ? cp : cr ;
	  if (cp < namBuf + 2 || cp[-2] != '.' || cp[-1] != 'a')
	    continue ;  /* reject other transcripts */
	}
      if (format == FASTC && !isGenome)
	{
	  char *cp = strchr (namBuf, '#') ;
	  if (cp)
	    {
	      int k = atoi (cp+1) ;
	      if (k > 1) mult = k ;
	    }
	}
      for (int iMult = 0 ; iMult < mult ; iMult++)
	{
	  memset (cr, 0, 10) ;
	  if (mult > 1)
	    {
	      sprintf (cr, ".%d", iMult+1) ; 
	    }
	  if (pairedEnd || arrayMax(dna2)) 
	    {
	      int k = strlen (namBuf) ;
	      int nn1, n1 = arrayMax (dna1) ;
	      int nn2, n2 = arrayMax (dna2) ;
	      
	      bb->rc.pairedEnd = pairedEnd = TRUE ;
	      nSeqs += 2 ;
	      bb->nSeqs += 2 ;
	      bb->length += n1 + n2 ;
	      bb->runStat.p.nPairs++ ;
	      bb->runStat.p.nReads += 2 ;
	      bb->runStat.p.nBase1 += n1 ;
	      bb->runStat.p.nBase2 += n2 ;
	      if (namBuf[k-1] == '>') k-- ;
	      /* namBuf[k] = '>' ; namBuf[k+1] = 0 ; */
	      namBuf[k] = 0 ;
	      dictAdd (bb->dict, namBuf, &nn1) ;
	      nn1 <<= 1 ;
	      array (bb->dnas, nn1, Array) = dna1 ;
	      if (bb->quals && qual1)
		array (bb->quals, nn1, Array) = qual1 ;

	      if (arrayMax (dna1))
		{
		  int i, iMax = arrayMax (dna1) ;
		  unsigned char *cp = arrp (dna1, 0, unsigned char) ;
		  for (i = 0 ; i < iMax && i < LETTERMAX ; i++, cp++)
		    bb->runStat.p.letterProfile1[5*i + atgcn[(int)*cp]]++ ;
		}
	      if (arrayMax (dna2))
		{
		  int i, iMax = arrayMax (dna2) ;
		  unsigned char *cp = arrp (dna2, 0, unsigned char) ;
		  for (i = 0 ; i < iMax && i < LETTERMAX ; i++, cp++)
		    bb->runStat.p.letterProfile2[5*i + atgcn[(int)*cp]]++ ;
		}
	      if (arrayMax (dna1))
		{
		  int i, iMax = arrayMax (dna1) ;
		  unsigned char *cp = arrp (dna1, 0, unsigned char) ;
		  for (i = 0 ; i < iMax ; i++, cp++)
		    bb->runStat.p.ATGCN[atgcn[(int)*cp]]++ ;

		  if (iMax > bb->runStat.p.maxReadLength)
		    bb->runStat.p.maxReadLength = iMax ;
		  if (! bb->runStat.p.minReadLength || iMax < bb->runStat.p.minReadLength)
		    bb->runStat.p.minReadLength = iMax ;
		  array (bb->runStat.p.lengthDistribution, iMax, long int)++ ;
		}
	      if (arrayMax (dna2))
		{
		  int i, iMax = arrayMax (dna1) ;
		  unsigned char *cp = arrp (dna2, 0, unsigned char) ;
		  for (i = 0 ; i < iMax ; i++, cp++)
		    bb->runStat.p.ATGCN[atgcn[(int)*cp]]++ ;

		  if (iMax > bb->runStat.p.maxReadLength)
		    bb->runStat.p.maxReadLength = iMax ;
		  if (! bb->runStat.p.minReadLength || iMax < bb->runStat.p.minReadLength)
		    bb->runStat.p.minReadLength = iMax ;
		  array (bb->runStat.p.lengthDistribution, iMax, long int)++ ;
		}

	      if (iMult < mult - 1)
		{
		  dna1 = arrayHandleCopy (dna1, bb->h) ;
		  if (bb->quals)
		    qual1 = arrayHandleCopy (qual1, bb->h) ;
		}
	      
	      else
		{
		  dna1 = arrayHandleCreate (256, unsigned char, bb->h) ;
		  if (bb->quals)
		    qual1 = arrayHandleCreate (256, unsigned char, bb->h) ;
		}
		  
	      /*
		namBuf[k] = '<' ;  namBuf[k+1] = 0 ;
		dictAdd (bb->dict, namBuf, &nn2) ;
	      */
	      nn2 = nn1 | 0x1 ;
	      if (arrayMax (dna2))
		array (bb->dnas, nn2, Array) = dna2 ;
	      if (bb->quals && qual2)
		array (bb->quals, nn2, Array) = qual2 ;
	      if (iMult < mult - 1)
		{
		  dna2 = arrayHandleCopy (dna2, bb->h) ;
		  if (bb->quals)
		    qual2 = arrayHandleCopy (qual2, bb->h) ;
		}
	      else
		{
		  dna2 = arrayHandleCreate (256, unsigned char, bb->h) ;
		  if (bb->quals)
		    qual2 = arrayHandleCreate (256, unsigned char, bb->h) ;
		}
	    }
	  else
	    {
	      int nn1, n1 = arrayMax (dna1) ;
	      int new = FALSE ;
	      
	      nSeqs++ ;
	      bb->nSeqs++ ;
	      bb->length += n1 ;
	      bb->runStat.p.nReads++ ;
	      bb->runStat.p.nBase1 += n1 ;
	      new = dictAdd (bb->dict, namBufX, &nn1) ;
	      nn1 <<= (isGenome ? 0 : 1) ;
	      if (isGenome && ! new)
		messcrash ("\nDuplicate target sequence identifier %s at line %d of file %s\nThe doublon may occur in this file or in a previous target file\n"
			   , namBufX
			   , line10
			   , fileName1
			   ) ;
	      line10 = line1 + 1 ;
	      array (bb->dnas, nn1, Array) = dna1 ;
	      if (arrayMax (dna1))
		{
		  int i, iMax = arrayMax (dna1) ;
		  unsigned char *cp = arrp (dna1, 0, unsigned char) ;
		  for (i = 0 ; i < iMax ; i++, cp++)
		    bb->runStat.p.ATGCN[atgcn[(int)*cp]]++ ;

		  if (iMax > bb->runStat.p.maxReadLength)
		    bb->runStat.p.maxReadLength = iMax ;
		  if (! bb->runStat.p.minReadLength || iMax < bb->runStat.p.minReadLength)
		    bb->runStat.p.minReadLength = iMax ;
		  array (bb->runStat.p.lengthDistribution, iMax, long int)++ ;
		}
	      if (arrayMax (dna1))
		{
		  int i, iMax = arrayMax (dna1) ;
		  unsigned char *cp = arrp (dna1, 0, unsigned char) ;
		  for (i = 0 ; i < iMax && i < LETTERMAX ; i++)
		    bb->runStat.p.letterProfile1[5*i + atgcn[(int)*cp]]++ ;
		}

	      if (bb->quals && qual1)
		array (bb->quals, nn1, Array) = qual1 ;
	      if (iMult < mult - 1)
		{
		  dna1 = arrayHandleCopy (dna1, bb->h) ;
		  if (qual1)
		    qual1 = arrayHandleCopy (qual1, bb->h) ;
		}
	      else
		{
		  dna1 = arrayHandleCreate (256, unsigned char, bb->h) ;
		  if (pp->sam)
		    qual1 = arrayHandleCreate (256, unsigned char, bb->h) ;
		}
	    }
	}
      if (chan && (bb->nSeqs > NMAX || bb->length > BMAX))
	{
	  t2 = clock () ;
	  	  
	  saCpuStatRegister ("2.ReadParser", pp->agent, bb->cpuStats, t1, t2, arrayMax (dnas) - 1) ;
	  
	  globalDnaCreate (bb) ;
	  nPuts++ ;
	  channelPut (chan, bb, BB) ;
	  
	  memset (bb, 0, sizeof (BB)) ;
	  t1 = clock () ;
	  
	  nn = 0 ;
	  bb->readerAgent = pp->agent ;
	  bb->start = timeNow () ;
	  bb->lane = atomic_fetch_add (arrp (pp->runLanes, bb->run, atomic_int), 1) + 1 ;
	  bb->h = ac_new_handle () ;
	  bb->txt1 = vtxtHandleCreate (bb->h) ;
	  bb->txt2 = vtxtHandleCreate (bb->h) ;
	  bb->errors = arrayHandleCreate (256, int, bb->h) ;
	  bb->cpuStats = arrayHandleCreate (128, CpuSTAT, bb->h) ;
	  bb->run = rc ? rc->run : 0 ;
	  bb->isGenome = isGenome ;
	  bb->nSeqs = 0 ;
	  bb->length = 0 ;
	  bb->dnas = dnas = arrayHandleCreate (nn, BigArray, bb->h) ;
	  if (pp->exportSamQuality && format == FASTQ)
	    bb->quals = arrayHandleCreate (64, Array, bb->h) ;
	  bb->runStat.p.lengthDistribution = arrayHandleCreate (1024, long int, bb->h) ;
	  bb->runStat.insertLengthDistribution = arrayHandleCreate (1024, long int, bb->h) ;
	  bb->dict = dictHandleCreate (NMAX, bb->h) ;
	  bb->errDict = dictHandleCreate (NMAX, bb->h) ;
	  dna1 = arrayHandleCreate (256, unsigned char, bb->h) ;
	  dna2 = arrayHandleCreate (256, unsigned char, bb->h) ;
	  if (pp->exportSamQuality && format == FASTQ)
	    {
	      qual1 = arrayHandleCreate (256, unsigned char, bb->h) ;
	      qual2 = arrayHandleCreate (256, unsigned char, bb->h) ;
	    }
	}
    }

  if (! bb->nSeqs)
    messcrash ("\nEmpty sequence file %s\n", fileName1) ;
  

  if (chan)
    {
      globalDnaCreate (bb) ;
      t2 = clock () ;
      
      if (dnas)
	saCpuStatRegister (isGenome ? "1.GParser" : "2.ReadParser", pp->agent, bb->cpuStats, t1, t2, arrayMax (dnas) - 1) ;
      nPuts++  ;
      channelPut (chan, bb, BB) ;
      channelPut (pp->npChan, &nPuts, int) ; /* global counting of BB blocks accross all sequenceParser agents */
    }
  
  /* create the REVERSE COMPLEMENT of the GENOME */
  if (isGenome == 2)
    {                       /* we need the minus strand of the genome */
      int iMax = bb->nSeqs ;
      globalDnaCreate (bb) ;
      bb->globalDnaR = bigArrayHandleCopy (bb->globalDna, bb->h) ;
      
      bb->dnasR = arrayHandleCreate (iMax, Array, bb->h) ;
      unsigned char *cp0 = bigArrayp (bb->globalDnaR, 0, unsigned char) ;
      for (int ii = 1 ; ii <= iMax ; ii++)
	{
	  Array dnaR = arrayHandleCreate (8, unsigned char, bb->h) ;
	  unsigned int x1 = bigArr (bb->dnaCoords, 2*ii, unsigned int) ;      /* offset of this DNA */
	  unsigned int x2 = bigArr (bb->dnaCoords, 2*ii + 1, unsigned int) ;
	  messfree (dnaR->base) ;
	  arrayLock (dnaR) ;
	  dnaR->base = (char *) cp0 + x1 ;
	  dnaR->max = dnaR->dim = x2 - x1 ;
	  reverseComplement (dnaR)  ;             /* complement in place */
	  array (bb->dnasR, ii, Array) = dnaR ;
	}
    }
  
  if (1 || pp->debug) printf ("--- %s: Stop sequence parser %d sequences, %d lines from file %s\n", timeBufShowNow (tBuf), nSeqs, line1, fileName1) ;
  
  ac_free (h) ;
  return ;
} /* otherSequenceParser */

/**************************************************************/

void saSequenceParse (const PP *pp, RC *rc, TC *tc, BB *bb, int isGenome)
{
  DnaFormat format = rc ? rc->format : tc->format ;
  if (0 && ! isGenome && ! rc->pairedEnd)
    return saUringSequenceParser (pp, rc, tc, bb) ;

  if (! isGenome)
    {
      if (format == SRA ||  pp->sraCaching)
	return sraSequenceParser (pp, rc, tc, bb, isGenome) ;
      else if ((format == FASTA && !rc->pairedEnd) ||
	       format == FASTA2       
	       )
	return fastaSequenceParser (pp, rc, tc, bb, isGenome) ;
    }

  return otherSequenceParser (pp, rc, tc, bb, isGenome) ;
} /* saSequenceParse */

/**************************************************************/

int saSequenceParseSraDownload (PP *pp, const char *sraID)
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

      if (! fastq)
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

      if (1) /* fastq contains fasta, no need to download again */
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
      else
      
      if (cr)
	{            /* file already in cache */
	  fprintf (stderr, "File %s already in cache\n", cr) ;	  
	  ac_free (h) ;
	  return 0 ;
	}
    }}
  /* download */
  
  SRAObj* sra = SraObjNew(sraID);
  const char *seq, *seq2 ;
  int nn = 0 ;
  BOOL firstPass = TRUE ;
  int num_bases = 1 << 27 ; /* 128 M */
  float Gb = pp->maxSraGb ;
  long unsigned int nMax = Gb * (1000000000L) ;  /* to be in decimal Gigabases */
  if (num_bases > nMax) num_bases = nMax ;
  nMax =  (nMax + num_bases - 1) / num_bases ; 
  fprintf (stderr, "%s : SRA download %s ", timeBufShowNow(tBuf), sraID) ;
  if (Gb) fprintf (stderr, "(top %.3f GigaBases) ", Gb) ;
  
  while (! Gb || nMax-- > 0)
    {
      SraGetReadBatch(sra, num_bases, fastq, split_pairs, &seq, &seq2) ;
      if (seq)
	{
	  nn++ ;
	  if (firstPass)
	    sraCachingOut (&ao1, &ao2, seq, seq2, sraID, fastq, split_pairs, h) ;
	  firstPass = FALSE ;
	  fprintf (stderr, ".") ;

	  if (!ao1 && ! ao2)
	    break ;
	  if (ao1 || ao2)
	    sraCacheDo (ao1, ao2, seq, seq2) ;
	}
    }
  SraObjFree(sra);
  if (ao1)
    fprintf (stderr, " %s downloaded %d data blocks infile %s\ndone: %s\n", sraID, nn, aceOutFileName (ao1), timeBufShowNow(tBuf)) ;
  
  ac_free (h) ;
  return 0 ;
} /* saSequenceParseSraDownload */

/**************************************************************/
/**************************************************************/
/**************************************************************/
 
/* eof */ 
