/*
 * sa.sequenceParse.c

 * This module is part of the sortalign package
 * Created May 15, 2026

 * This code is public.

 * This module downloads runs from the SRA sequence archives
 * And parses local sequence files in fasta/fastq format
 * into SraParse structures
 */

#include "sa.h"

#define NAMMAX 1024

/* download SRA runs as FASTA or FASTQ */
#define SRA_FASTA 0
#define SRA_FASTQ 1

/********************************************************************************************************************/
/*************************** encode SAPARSE into standard dna arrays ************************************************/
/********************************************************************************************************************/

/*  scan the saParse buffers, encode the dna, create virtual dna arrays */
static int saScanDnaEncode (SAPARSE *sap)
{
  int nn = 0 ;
  int ii, iMax = sap->recordMax ;
  DnaRecord *up ;
  unsigned int maxDnaLn = 0, minDnaLn = 1 << 31 ;
  for (ii = 0, up = sap->records ; ii < iMax ; ii++, up++)
    {
      int n = 0 ;
      unsigned char cc, *cp = sap->dnaBuffer + up->dna ;
      while (cc = *cp, dnaEncodeChar[(int)*cp++])
	n++ ;
      if (n < up->dnaLn)
	{
	  cp = sap->dnaBuffer + up->dna ;
	  messcrash ("Wrong character %c at position %d in dna sequence %s\n"
		     , sap->idBuffer + up->id
		     , cc, n
		     ) ;
	}
      minDnaLn = (n < minDnaLn ? n : minDnaLn) ;
      maxDnaLn = (n > maxDnaLn ? n : maxDnaLn) ;
      nn += n ;
    }
  sap->nBases = nn ;
  sap->minDnaLn = minDnaLn ;
  sap->maxDnaLn = maxDnaLn ;
  
  return nn ;
} /* saScanDnaEncode */

/********************************************************************************************************************/
/* create virtual dna arrays pointing inside the dnaBuffer
 * to allow transparent usage of the standard acedb dna array libraries
 */

static void saScanVirtualArrayCreate (BB *bb, SAPARSE *sap)
{
  Array dnas ;
  DnaRecord *up ;
  int ii, iMax = sap->recordMax ;

  dnas  = bb->dnas = arrayHandleCreate (iMax, Array, bb->h) ;
  for (ii = 0, up = sap->records ; ii < iMax ; ii++, up++)
    array (dnas, ii, Array) = virtualArrayCreate (sap->dnaBuffer + up->dna,  up->dnaLn, unsigned char) ;
  return ;
} /* saScanVirtualArrayCreate */

/********************************************************************************************************************/
/*  scan the saParse buffers, encode the dna, create virtual dna arrays */
int saScan (const PP *pp, BB *bb)
{
  SAPARSE *sap = bb->saParse ;
  int iMax = sap->recordMax ;

  if (! iMax) return 0 ;

  saScanDnaEncode (sap) ;
  saScanVirtualArrayCreate (bb, sap) ;
  return iMax ;
} /* saScan */

/********************************************************************************************************************/
/******************************** sequence file parsing utilities ***************************************************/
/********************************************************************************************************************/
/* readScanId
 * Works on a local copy of a fraction of a user provided private BB buffer. Edit in place.
 * Expects to start on >, fasta, or @, fastc
 * Sets suffix (0, 1, 2) for paired end case
 * On return buf (or bu+1 if prefix) is the clean name for the pair stripping /1
 * Returns number of consumed char, including \n
 */
static int readScanId (unsigned char *buf, char prefix, int *suffix)
{
  unsigned char *cp = buf, *cq ;    
  int ln = 0 ;
  
  if (prefix)
    {
      if (prefix != *cp)    /* cp should point to prefix '>' or '@' */
	{
	  messerror ("Wrong identifier %s\n", cp) ;
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

  if (suffix)
    {
      *suffix  = 0 ;
      if (ln > 2)
	{
	  cq = cp + ln - 2 ;
	  if (*cq == '/' || *cq == '.' || *cq == ' ')
	    {
	      cq++ ;
	      if (*cq == '1') *suffix = 1 ;
	      if (*cq == '2') *suffix = 2 ;
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
static int readScanQual (unsigned char *buf)
{
  int ln = 0 ;
  unsigned char *cq = ustrchr (buf, '\n') ;
  if (cq)
    {
      *cq = 0 ;
      ln = cq - buf + 1 ;
    }
  else
    ln = ustrlen (buf) ;
  return ln ;
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
	case '>': *cq++ = 0 ; goto done ; break ;
	case '\r': break ;
	case '\n': break ;
	default: *cq++ = *cp ; break ;
	}
    }
 done:
  ln = cp - buf + 1 ;
  return ln ;
} /* readScanDna */

/******************************************************************/
/* readScan fasta/fastq file chunk
 * store it in SAPARSE struct
 * Chunk is supposed to contain a full number of read (or pairs) records
 */

SAPARSE *saParseCreateFromFileChunk (BB *bb)
{
  SAPARSE *saParse = halloc (sizeof(SAPARSE), bb->h) ;
  // RC *rc = bb->rc ;
  unsigned char *buf = bb->gzBuffer ;

  memset (saParse, 0, sizeof (SAPARSE)) ;
	  


  return saParse ;
} /* saParseCreateFromFileChunk */


/********************************************************************************************************************/
/******************************** sequence file parsing utilities ***************************************************/
/********************************************************************************************************************/

/********************************************************************************************************************/
/******************************** sequnece file parsing utilities ***************************************************/
/********************************************************************************************************************/


