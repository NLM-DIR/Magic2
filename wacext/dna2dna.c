/*
 * Authors: Jean Thierry-Mieg, NCBI, 
 * Oct 2009
 * Multilingual dna parser
 *   Input: raw, fasta, fastq, solid
 *   Output: raw, fasta, solid
 *   Filters: length, entropy
 *
*/

#define BITSET_MAKE_BITFIELD   
#include <ac.h>
#include <waceview/cgi.h>
#include <bitset.h>

typedef enum { FASTA=0, FASTC, FASTQ, CSFASTA, CSFASTC, CCFA, CCFAR, RAW,  CRAW, TAG_COUNT, CTAG_COUNT, COUNT, TENSORFLOW, LETTER, LETTERPLOT, LETTERSHOW, POLYMER} DNAFORMAT ;

#define MAXMULT 16
#define POLYMER_MAX 1010
typedef struct seqStruct { int dna, count, dnaLn ; vTXT id ;  } SEQ ;
typedef struct shadowStruct { int mrna, gene, x1, x2, target, a1, a2, ID, Parent, type, ncRNA, GeneID, gene_name, gene_title, name, gene_biotype, note, description, Dbxref,protein_id,locus_tag ; BOOL isDown ; BOOL exon, cds ; } SHADOW ;

typedef struct sxStruct {  
  DNAFORMAT in, out ; 
  int minEntropy, minLength, maxLength, minQuality, minMultiplicity, split, splitMb ;
  long int nProcessed, nRejected , fProcessed, fRejected, fmProcessed, fmRejected ;
  int  nSplit, nExported, count, maxCount, minLn, maxLn, splitByPrefix, letterNN, polymer ;
  long int *letterCount, *letterProfile, *letterLength ; 
  long int polymerCount[4], polymerProfile[4*POLYMER_MAX] ;
  long int tagProcessed, tagRejected ;
  long int bpSeqProcessed, bpSeqRejected ;
  long int bpTagsProcessed, bpTagsRejected ;
  int clipN, leftClipAt, rightClipAt, maxLineLn, selectPrefixLn ;
  int subsample ;
  int makeTest, makeTestLength, makeTestPairLength, makeTestReverse, makeTestStep, makeTestShift, makeTestPeriod, makeTestInDel, makeTestAddPolyA, makeTestRandom, makeTestSubRate ;
  BOOL n2a, jumpAnchor, complement, decoy, qualityPlot, upper, makeTestGenome ; 
  ACEIN ai, ai1, ai2, aiId ; 
  ACEOUT ao, aoId, aos[256], aos1[256], aos2[256], aoTm, ao1, ao2 ; 
  ACEOUT aoTestSeqs, aoTestSnps ;
  DICT *dict ; 
  vTXT id, dna, qual, qual2 ; 
  Array tagArray, encodedDna ; 
  const char *inFileName, *inFileName1, *inFileName2, *inIdFileName, *outFileName, *outIdFileName, 
    *prefix, *selectFileName, *rejectFileName, *getList, *title, *runTitle, *gtfRemapPrefix ;
  const char *suffix1, *suffix2 ;
  const char *leftClipOn, *rightClipOn, *leftClipInPhase ;
  char *selectPrefix ;
  const char *fastqSelect, *runQualityFileName ;
  const char *shadowFileName, *spongeFileName, *gtfFileName, *gtfGenomeFastaFileName, *refGeneFileName ;
  char fileType[16] ;
  BOOL selectDictHasSpace ;
  DICT *selectDict, *rejectDict, *shadowDict ;
  Array shadowArray, gnam2genes ;
  Array lnDistrib ;
  KEYSET rnaId2mrna ;
  int xor ; /* compare 2 sequences, parameter -xor */
  BOOL gzi, gzo, doCount, getTm, appendByPrefix, keepName , keepName1 , keepName2, keepNameBam, isGff3, gffTAIR ;
  BOOL fastc_paired_end, gffBacteria ;
  BOOL noNameCheck ;
  BOOL makeFastaUnique ;
  unsigned char **entryAdaptor ;
  unsigned char **inPhaseEntryAdaptor ;
  int dnaLn ; /* length of current dna */
  int rawColumn, ynColumn ; /* DNA column in tab delimited raw -I raw%d format: raw3 means DNA is in column 3 */
  int entryAdaptorSeqFound ;
  int entryAdaptorTagFound ;
  long int entryAdaptorBpFound ;
  long int entryAdaptorBpClipped ;
  unsigned char **exitAdaptor ;
  int exitAdaptorSeqFound ;
  int exitAdaptorTagFound ;
  int UMI ;
  long int exitAdaptorBpFound ;
  long int exitAdaptorBpClipped ;
  AC_HANDLE h ;
} SX ;

/*************************************************************************************/
/* a cfa sequence has an Anchor letter, followed by transition codes
 * the length does not change
 */
static void fa2ccfa (SX *sx)
{
  char cc, old, *cp ;
  
  cp = vtxtPtr (sx->dna) ;
  old = ace_lower(*cp) ; 
  *cp = ace_upper(*cp) ; /* set the anchor in upper case */
  if (! index("atgc",old)) old = 'a' ;

  while (*++cp)
    { 
      cc =  (UT_LOWER[(*cp) & 0x7f]) ;
      /* cc = ace_lower [*cp] ; */
      switch (cc)
	{
	default:
	  cc = 'n' ;
	  continue ; /* do not edit the old letter 
			so next time we export the transition jumping the n
		     */
	case '<':
	  *cp = '<' ;
	  old = 'a' ;
	  continue ;
	  break ;
	case '>':
	  *cp = '>' ;
	  old = 'a' ;
	  continue ;
	  break ;
	case 'a':
	  switch (old)
	    {
	    case 'a': *cp = 'a' ; break ;
	    case 't': *cp = 't' ; break ;
	    case 'g': *cp = 'g' ; break ;
	    case 'c': *cp = 'c' ; break ;
	    }
	  break ;
	case 't':
	  switch (old)
	    {
	    case 'a': *cp = 't' ; break ;
	    case 't': *cp = 'a' ; break ;
	    case 'g': *cp = 'c' ; break ;
	    case 'c': *cp = 'g' ; break ;
	    }
	  break ;
	case 'g':
	  switch (old)
	    {
	    case 'a': *cp = 'g' ; break ;
	    case 't': *cp = 'c' ; break ;
	    case 'g': *cp = 'a' ; break ;
	    case 'c': *cp = 't' ; break ;
	    }
	  break ;
	case 'c':
	  switch (old)
	    {
	    case 'n': *cp = 'a' ; break ;
	    case 'a': *cp = 'c' ; break ;
	    case 't': *cp = 'g' ; break ;
	    case 'g': *cp = 't' ; break ;
	    case 'c': *cp = 'a' ; break ;
	    }
	  break ;
	}
      old = cc ;
    }
  return ;
}  /* fa2ccfa */

/*************************************************************************************/
/* various SOLiD color formats, do not edit the Anchor base */
static void cfa2fa (SX *sx)
{
  char cc, old, *cp ;
  
  cp = vtxtPtr (sx->dna) ;
  old = *cp = (UT_LOWER[(*cp) & 0x7f]) ; /* anchor base */
  while ((cc = *++cp))
    {
      switch ((int)cc)
	{
	default:
	  *cp = 'n' ;
	  continue ; /* do not edit old */
	  break ;
	case '<':
	  *cp = '<' ;
	  old = 'a' ;
	  continue ;
	  break ;
	case '>':
	  *cp = '>' ;
	  old = 'a' ;
	  continue ;
	  break ;
	case '0': 
	  switch ((int)old)
	    {
	    case 'a': case 'A': *cp = 'a' ; break ;
	    case 't': case 'T': *cp = 't' ; break ;
	    case 'g': case 'G': *cp = 'g' ; break ;
	    case 'c': case 'C': *cp = 'c' ; break ;
	    }
	  break ;
	case '1': 
	  switch ((int)old)
	    {
	    case 'a': case 'A': *cp = 'c' ; break ;
	    case 't': case 'T': *cp = 'g' ; break ;
	    case 'g': case 'G': *cp = 't' ; break ;
	    case 'c': case 'C': *cp = 'a' ; break ;
	    }
	  break ;
	case '2': 
	  switch ((int)old)
	    {
	    case 'a': case 'A': *cp = 'g' ; break ;
	    case 't': case 'T': *cp = 'c' ; break ;
	    case 'g': case 'G': *cp = 'a' ; break ;
	    case 'c': case 'C': *cp = 't' ; break ;
	    }
	  break ;
	case '3': 
	  switch ((int)old)
	    {
	    case 'a': case 'A': *cp = 't' ; break ;
	    case 't': case 'T': *cp = 'a' ; break ;
	    case 'g': case 'G': *cp = 'c' ; break ;
	    case 'c': case 'C': *cp = 'g' ; break ;
	    }
	  break ;
	}
      old = *cp ;
    }
  return ;
} /* cfa2fa */

/*************************************************************************************/
/* various SOLiD color formats, do not edit the Anchor base */
static void cfa2cfa (SX *sx)
{
  char cc, *cp = vtxtPtr (sx->dna) ;
  
  switch (sx->in)
    {
    default: /* dna has already been mapped to CCFA */
    case CCFA:
      switch (sx->out)
	{
	case CSFASTA:
	case CSFASTC:
	  while ((cc= *(++cp)))
	    switch (cc)
	      {
	      case 'n' : *cp = '.' ; break ;
	      case '<' : *cp = '<' ; break ;
	      case '>' : *cp = '>' ; break ;
	      case 'a' : *cp = '0' ; break ;
	      case 'c' : *cp = '1' ; break ;
	      case 'g' : *cp = '2' ; break ;
	      case 't' : *cp = '3' ; break ;
	      }
	  break ;
	case CCFAR:
	  while ((cc= *(++cp)))
	    switch (cc)
	      {
	      case 'n' : *cp = 'n' ; break ;
	      case '<' : *cp = '<' ; break ;
	      case '>' : *cp = '>' ; break ;
	      case 'a' : *cp = 't' ; break ;
	      case 'c' : *cp = 'g' ; break ;
	      case 'g' : *cp = 'c' ; break ;
	      case 't' : *cp = 'a' ; break ;
	      }
	  break ;
	default:
	  break ;
	}
      break ;
    case CCFAR:
      switch (sx->out)
	{
	case CSFASTA:
	case CSFASTC:
	  while ((cc= *(++cp)))
	    switch (cc)
	      {
	      case 'n' : *cp = '.' ; break ;
	      case '<' : *cp = '<' ; break ;
	      case '>' : *cp = '>' ; break ;
	      case 'a' : *cp = '3' ; break ;
	      case 'c' : *cp = '2' ; break ;
	      case 'g' : *cp = '1' ; break ;
	      case 't' : *cp = '0' ; break ;
	      }
	  break ;
	case CCFA:
	  while ((cc= *(++cp)))
	    switch (cc)
	      {
	      case 'n' : *cp = 'n' ; break ;
	      case '<' : *cp = '<' ; break ;
	      case '>' : *cp = '>' ; break ;
	      case 'a' : *cp = 't' ; break ;
	      case 'c' : *cp = 'g' ; break ;
	      case 'g' : *cp = 'c' ; break ;
	      case 't' : *cp = 'a' ; break ;
	      }
	  break ;
	default:
	  break ;
	}
      break ;
    case CSFASTA:
    case CSFASTC:
      switch (sx->out)
	{
	case CCFA:
	  while ((cc= *(++cp)))
	    switch (cc)
	      {
	      case '.' : *cp = 'n' ; break ;
	      case '<' : *cp = '<' ; break ;
	      case '>' : *cp = '>' ; break ;
	      case '0' : *cp = 'a' ; break ;
	      case '1' : *cp = 'c' ; break ;
	      case '2' : *cp = 'g' ; break ;
	      case '3' : *cp = 't' ; break ;
	      }
	  break ;
	case CCFAR:
	  while ((cc= *(++cp)))
	    switch (cc)
	      {
	      case '.' : *cp = 'n' ; break ;
	      case '<' : *cp = '<' ; break ;
	      case '>' : *cp = '>' ; break ;
	      case '0' : *cp = 't' ; break ;
	      case '1' : *cp = 'g' ; break ;
	      case '2' : *cp = 'c' ; break ;
	      case '3' : *cp = 'a' ; break ;
	      }
	  break ;
	default:
	  break ;
	}
      break ;
    }
  return ;
} /* cfa2ccfa */

/*************************************************************************************/
/*************************************************************************************/

static int seqEntropy (const char *seq, int minEntropy)
{
  const char *ccp = seq ;
  unsigned char buf[1024] ;
  int i, n = strlen (seq) ;
  BOOL isDown = TRUE ;

  if (n > 1023) n = 1023 ;
  for (i = 0, ccp = seq ; i < n ; ccp++, i++)
    switch ((int)*ccp)
      {
      case '<': isDown = FALSE ; buf[i] = FS_ ; break ;
      case '>': isDown = TRUE ; buf[i] = RS_ ; break ;
      case 'a': case 'A': case '0': buf[i] = isDown ? A_ : T_ ; break ;
      case 't': case 'T': case '1': buf[i] = isDown ? T_ : A_ ;break ;
      case 'g': case 'G': case '2': buf[i] = isDown ? G_ : C_ ;break ;
      case 'c': case 'C': case '3': buf[i] = isDown ? C_ : G_ ;break ;
      default: buf[i] = N_ ;
      }
  buf[i] = 0 ;
  return oligoEntropy (buf, n, minEntropy) ;
} /* seqEntropy */

/*************************************************************************************/

static BOOL entropyReject (SX *sx)
{
  if (sx->minMultiplicity && sx->minMultiplicity > sx->count)
    return TRUE ;
  if (sx->minEntropy)
    {
      int s = seqEntropy (vtxtPtr (sx->dna), sx->minEntropy) ;
      if (s < sx->minEntropy)
	return TRUE ;
    }
  return FALSE ;
} /* entropyReject */

/*************************************************************************************/

static const char* getNextLine (SX *sx, BOOL jumpComment, ACEIN ai)
{
  const char *ccp ;
  
  while (aceInCard (ai))
    {
      if (! (ccp = aceInPos (ai)))
	return NULL ;
      if (jumpComment && *ccp == '#') /* jump lines starting with a # */
	continue ;
      if (! sx->gtfGenomeFastaFileName && /* no full name if we compare to the gff */
	  (sx->in == FASTQ || sx->in == FASTA || sx->in == FASTC)
	  ) /* fastq CASAVA has spaces in its identifiers */
	return ccp ;
      return aceInWord (ai) ;
    }
  return NULL ;
} /* getNextLine  */

/*************************************************************************************/

static BOOL aceInDnaCheck (SX* sx, char *ccp, BOOL ignoreAnchor)
{
  char *ccq ;
  char cc, c2 ;
  ccq = ccp - 1 ;
  if (sx->in == CSFASTA || sx->in == CSFASTC)
    {
      if (!ignoreAnchor)
	{
	  cc = dnaEncodeChar[(int)*++ccq] ;
	  if (! cc)
	    messcrash ("Bad Anchor char in dna sequence %s, line %d", ccp, aceInStreamLine (sx->ai)) ;
	  *ccq = dnaDecodeChar[(int)cc] ;
	}
      while ((cc = *++ccq))
	{
	  switch ((int)cc)
	    {
	    case '0':
	    case '1':
	    case '2':
	    case '3': break ;
	    case '.': break ;
	    case '<': break ;
	    case '>': break ;
	    default:
	    messerror ("Bad char %c in cfa sequence %s, line %d, please consider using the -n2a option", cc, ccp, aceInStreamLine (sx->ai)) ;
	    return FALSE ;
	    }
	}
    }
  else
    while ((cc = *++ccq))
      {
	if (cc == '-') cc = 'n' ;
	if (cc == '<') continue ;
	if (cc == '>') continue ;
	c2 = dnaEncodeChar[(int)cc] ;
	if (c2)
	  *ccq = dnaDecodeChar[(int)c2] ;
	else
	  {
	    if (0 && sx->in != FASTQ)
	      messcrash ("Bad char in dna sequence %s, line %d", ccp, aceInStreamLine (sx->ai)) ;
	    else
	      {
		messerror ("Bad char in dna sequence %s, line %d", ccp, aceInStreamLine (sx->ai)) ;
		return FALSE ;
	      }
	  }	  
      }

  return TRUE ;
} /* aceInDnaCheck */

/*************************************************************************************/

static BOOL getDna (SX *sx, ACEIN ai)
{
  static int line = 0 ;
  char cutter, *ccp = 0 ;
  BOOL ok = FALSE ;

  if (1) aceInSpecial (ai, "\n") ;
  while (1)
    {
      line++ ;
      if (!aceInCard (ai))
	return ok ;
      if (sx->gtfGenomeFastaFileName)
	ccp = aceInWord (ai) ;
      else
	ccp = aceInWordCut (ai, "\t", &cutter) ;
      if (! ccp)
	return ok ;
      if (sx->in == RAW && ! sx->ynColumn && sx->rawColumn > 1)
	{
	  int i = 1 ;
	  if (sx->keepName)
	    {
	      vtxtClear (sx->id) ;
	      vtxtPrint (sx->id, ccp) ;
	    }
	  if (sx->keepNameBam && sx->rawColumn > 2)
	    {
	      i++ ; 
	      if (!aceInWordCut (ai, "\t", &cutter))
		return ok ;
	      if (strstr (ccp, "1"))
		vtxtPrint (sx->id, "/1") ;
	      else if (strstr (ccp, "2"))
		vtxtPrint (sx->id, "/2") ;
	    }
	  for (; i < sx->rawColumn - 1 ; i++)
	    aceInWordCut (ai, "\t", &cutter) ;
	  for (; i < sx->rawColumn ; i++)
	    if (! (ccp = aceInWordCut (ai, "\t", &cutter)))
	      return ok ;
	}
      if (sx->in == RAW && sx->ynColumn && sx->rawColumn)
	{
	  BOOL isY = sx->ynColumn <= sx->rawColumn ? TRUE : FALSE ;
	  int i, col = isY ? sx->ynColumn : sx->rawColumn ;

	  if (sx->keepName)
	    {
	      vtxtClear (sx->id) ;
	      vtxtPrint (sx->id, ccp) ;
	    }
	  for (i = 1 ; i < col ; i++)
	    if (!aceInStep (ai, '\t') || ! (ccp = aceInWord (ai)))
	      continue ;
	  if (isY && *ccp != 'Y')
	    continue ;
	  if (isY)  /* yes/No flag is before the dna, we just go on reading */
	    {
	      if (*ccp != 'Y')
		continue ; 
	      for (i = sx->ynColumn ; i < sx->rawColumn ; i++)
		{
		  if (!aceInStep (ai, '\t'))
		    continue ;
		  ccp = aceInWord (ai) ;
		}
	      if (! ccp)
		continue ;
	    }
	  else /* yes/No flag is behind the dna, we strnew the dna, check Y/N and free the memory */
	    {
	      int n = strlen(ccp) ;
	      char dna[n+1] ;

	      strcpy (dna, ccp) ;

	      for (i = sx->rawColumn ; i < sx->ynColumn ; i++)
		{
		  if (!aceInStep (ai, '\t'))
		    {
		      continue ;
		    }
		  ccp = aceInWord (ai) ;
		}
	      if (*ccp != 'Y')
		continue ;
 
	      if (!aceInDnaCheck (sx, dna, FALSE))
		{
		  fprintf(stderr, "Bad sequence %s\n", *dna ? dna : "NULL") ;
		  return FALSE ; 
		}
	      ok = TRUE ; 
	      vtxtPrint (sx->dna, dna) ;
	      return TRUE ;
	    }	  
	}
      if (*ccp == '>')
	{
	  aceInCardBack (ai) ;
	  return ok ; 
	}
      if (!aceInDnaCheck (sx, ccp, ok))
	{
	  if (sx->in == FASTQ)
	    return FALSE ;
	  else
	    continue ;
	}
      ok = TRUE ;
      vtxtPrint (sx->dna, ccp) ;

      if (sx->in == RAW && sx->keepName  && !sx->rawColumn)
	{
	  ccp = aceInPos(ai) ;
	  vtxtClear (sx->id) ;
	  if (ccp)
	    vtxtPrint (sx->id, ccp) ;
	}

      if (sx->in != FASTA && sx->in != FASTC) /* only fasta has dna on several lines */
	break ;
    }
  return ok ;
} /* getDna  */

/*************************************************************************************/

BOOL isIdKosher (SX *sx, vTXT id)
{
  if (sx->selectDict)
    {
      BOOL ok = FALSE ;
      
      if (dictFind (sx->selectDict, vtxtPtr (id), 0)) 
	ok = TRUE ;
      else
	{
	  char cc, *cp = strstr ( vtxtPtr (id), "|Gene|") ;
	  if (!cp) cp = strstr ( vtxtPtr (id), "|GENE|") ;
	  if (cp)
	    {
	      cc = *cp ; *cp = 0 ;
	      if (dictFind (sx->selectDict, vtxtPtr (id), 0)) 
		ok = TRUE ;
	      *cp = cc ;
	    }
	  if (0 && ! ok && sx->in == FASTQ && (cp = strstr ( vtxtPtr (id), "/")))
	    {
	      cc = *cp ; *cp = 0 ;
	      if (dictFind (sx->selectDict, vtxtPtr (id), 0)) 
		ok = TRUE ;
	      *cp = cc ;
	    }
	  if (! ok && ! sx->selectDictHasSpace)
	    {
	      char *cp = strstr ( vtxtPtr (id), " ") ;
	      if (cp)
		{
		  /* char cc = *cp ;  */
		  *cp = 0 ;
		  if (dictFind (sx->selectDict, vtxtPtr (id), 0)) 
		    ok = TRUE ;
		  /* *cp = cc ; 2019_-6_08 */
		}
	    }
	}
      if (!ok) return FALSE ;
    }
  if (sx->rejectDict)
    {
      BOOL ok = FALSE ;
      
      if (dictFind (sx->rejectDict, vtxtPtr (id), 0)) 
	ok = TRUE ;
      else
	{
	  char cc, *cp = strstr ( vtxtPtr (id), "|Gene|") ;
	  if (cp)
	    {
	      cc = *cp ; *cp = 0 ;
	      if (dictFind (sx->rejectDict, vtxtPtr (id), 0)) 
		ok = TRUE ;
	      *cp = cc ;
	    }	
	  if (! ok && ! sx->selectDictHasSpace)
	    {
	      char cc, *cp = strstr ( vtxtPtr (id), " ") ;
	      if (cp)
		{
		  cc = *cp ; *cp = 0 ;
		  if (dictFind (sx->rejectDict, vtxtPtr (id), 0)) 
		    ok = TRUE ;
		  *cp = cc ;
		}
	    }
	}
      if (ok) return FALSE ;
    }

  return TRUE ;
} /* isIdKosher */

/**************************************/

static BOOL sxChecInputHasPairs (SX *sx)
{
  if (sx->fastc_paired_end) 
    return TRUE ;
  if (! sx->inFileName) 
    return FALSE ;
  else
    {
      AC_HANDLE h = ac_new_handle () ;
      int nn = 100000 ;
      BOOL ok = FALSE ;
      ACEIN ai = aceInCreate (sx->inFileName, sx->gzi, h) ;
      while (nn-- && aceInCard (ai))
	{
	  char *cp = aceInPos (ai) ;
	  if (cp && strstr (cp, "><"))
	    { ok = TRUE ; break ; }
	}
      ac_free (ai) ;
      ac_free (h) ;
      return ok ;
    }
}

/**************************************/

static BOOL getSequencePart (SX *sx, ACEIN ai)
{
  vTXT dna = sx->dna ;
  vTXT qual = sx->qual ;
  vTXT id = sx->id ;
  const char *ccp, *ccq ;

  if (! sx->inFileName)
    sx->inFileName =  sx->inFileName1 ;

  vtxtClear (dna) ;
  if (sx->in == FASTQ) vtxtClear (qual) ;
  vtxtClear (id) ;

  switch (sx->in)
    {
    case TENSORFLOW:
    case COUNT:
    case POLYMER:
    case LETTER:
    case LETTERPLOT:
    case LETTERSHOW:
      return FALSE ;
    case FASTA:
    case FASTC:
    case CSFASTA:
    case CSFASTC:
    case CCFA:
    case CCFAR:
      /* get the sequence identifier */
      while (1)
	{
	  vtxtClear (id) ;
	  vtxtClear (dna) ;
	  if (! (ccp = getNextLine (sx, TRUE, ai)))
	    return FALSE ;
	  if (*ccp != '>')
	    {
	      messerror ("Missing > at beginning of sequence, line %d:\n%s\n", aceInStreamLine (ai), ccp) ;
	      continue ;
	    }
	  sx->count = 1 ;
	  if (sx->in == FASTC || sx->in == CSFASTC)
	    {
	      sx->count = fastcMultiplicity (ccp, 0, 0) ;
	    }

	  vtxtPrint (id, ccp+1) ;
	  /* get the sequence dna */
	  if (! getDna (sx, ai))
	    continue ;
	  if (isIdKosher (sx, id))
	    break ;
	}
      break ;

    case FASTQ:
      /* get the sequence identifier */
      while (1)
	{
	fastqIter: 
	  vtxtClear (id) ;
	  vtxtClear (dna) ;
	  while (1) /* wait for a good line, we may loose a few sequences */
	    {
	      if (! (ccp = getNextLine (sx, FALSE, ai)))
		return FALSE ; 
	      if (*ccp != '@') continue ;
	      ccq = 0 ;
              if (! sx->ai1 &&
		  (ccq = strstr (ccp, " ")) && 
		  sx->fastqSelect && strstr (sx->fastqSelect, ":N:") && 
		  strstr (ccq, ":Y:"))
		{
		  aceInCard (ai) ; aceInCard (ai) ; aceInCard (ai) ; 	  
		  continue ; /* chastitize all CASAVA8 reads */
		}
	      if (! sx->fastqSelect || (ccq && strstr (ccq, sx->fastqSelect)))
		break ;
	    }   
	  vtxtPrint (id, ccp+1) ;
	  {  /* drop ' length=  and replace space by / in sequence identifiers */
	    char *cp ;
	    int na ;

	    cp = strstr(vtxtPtr (id), " length=") ;
	    if (cp) *cp = 0 ;
	    na = strlen(vtxtPtr (id)) ;
	    cp = strstr (vtxtPtr (id), ":Y:") ;
	    if (! sx->ai1 && cp && cp - vtxtPtr (id) > na - 8)
	       {
		 aceInCard (ai) ; aceInCard (ai) ; aceInCard (ai) ; 	  
		 continue ; /* chastitize all CASAVA8 reads */
	       }
	    cp = strstr (vtxtPtr (id), ":N:") ;
	    if (cp && cp - vtxtPtr (id) > na - 8)
	      *cp = 0 ;  /* skip the ending */
	    for (cp = vtxtPtr (id) ; *cp ; cp++)
	      if (*cp == ' ') 
		*cp = '/' ;
	    na = strlen(vtxtPtr (id)) ;
	    if (sx->keepName1)
	      {
		cp = vtxtPtr (id) + na - 2 ;
		if (strncmp (cp, "/1",2))
		  vtxtPrintf(id, "/1") ; /* make room */
	      }
	    if (sx->keepName2)
	      {
		cp = vtxtPtr (id) + na - 2 ;
		if (strncmp (cp, "/2",2))
		  vtxtPrintf(id, "/2") ; /* make room */
	      }
	  }
	  /* get the sequence dna */
	  if (! getDna (sx, ai))
	    {
	      goto fastqIter ;
	      messcrash ("Missing DNA in fastq input file, file %s line %d:\n"
			 , sx->inFileName, aceInStreamLine (ai), ccp) ;
	    }
	  if (! (ccp = getNextLine (sx, FALSE, ai)))
	    messcrash ("Missing quality identifier in fastq file, file %s line %d:\n"
		       , sx->inFileName, aceInStreamLine (ai)) ;
	  if (*ccp != '+')
	    messcrash ("Missing + at beginning of quality identifier, file %s line %d:\n%s\n"
		       , sx->inFileName, aceInStreamLine (ai)
		       , ccp ? ccp : "") ;
	  
	  if (ccp[1])
	    {
	      const char *c1, *c2 = vtxtPtr (id) ;
	      int shift = 1 ;

	      /* 2013_07_28: fix a format bug found in fastq file for Ghs4007 Neurobalstome exome where the qual is name +@7:.. rather than +7:... */ 
	      if (ccp[1] == '@' && *c2 != '@' && ccp[2] == *c2) shift = 2 ; 
	      for (c1 = ccp+shift, c2 = vtxtPtr (id) ; *c1 || *c2 ; c1++, c2++)
		if (*c1 != ' ' && *c2 != *c1)
		  {
		    if (0 || !strncmp(c1,"length=",1) ||  !strncmp(c1,"/length=",8) || 
			(sx->keepName1 && (!strcmp(c2, "1") || !strcmp(c2, "/1"))) ||
			(sx->keepName2 && (!strcmp(c2, "2") || !strcmp(c2, "/2")))
			) break ;
		    else
		      {
			messerror ("Missmatch between the sequence identifier %s and the quality identifier %s, file %s line %d"
				   , vtxtPtr (id), ccp + shift
				   , sx->inFileName, aceInStreamLine (ai)) ; 
			goto fastqIter ; 
		      }
		  }
	    }
	  /* get the quality */
	  if (!(ccp = getNextLine (sx, FALSE, ai)))
	    messcrash ("Missing DNA in fastq input file %s line %d:\n"
		       , sx->inFileName, aceInStreamLine (ai)) ;
	  vtxtPrint (qual, ccp) ;
	  if (sx->minQuality)
	    {
	      char mqc = '@' + sx->minQuality  ;
	      char *vq = vtxtPtr(qual), *vs = vtxtPtr (sx->dna) - 1 ;
	      while (*++vs)
		if (*vq++ <= mqc) *vs = 'n' ;
	    }
	  sx->count = 1 ;
	  if (isIdKosher (sx, id))
	    break ;
	}
      break ;

    case RAW:
    case CRAW:
      /* get the sequence dna */
      if (! getDna(sx, ai))
	return FALSE ;
      /* name the sequence */
      if (sx->keepName)
	{
	    /* sx->rawColumn <= 1 &&    aceInStep (ai,'\t') && (ccp = aceInWord (ai)) ; vtxtPrint (sx->id, ccp) ; */

	  if (sx->keepName1)
	    vtxtPrint (sx->id, "/1") ;
	  if (sx->keepName2)
	    vtxtPrint (sx->id, "/2") ;
	}
      else
	{
	  vtxtClear (sx->id) ;
	  vtxtPrintf (sx->id, "%s.%d", sx->prefix, sx->nProcessed + 1) ;
	}
      sx->count = 1 ;
      break ;
    case TAG_COUNT:
    case CTAG_COUNT:
      /* get the sequence dna */
      if (! getDna (sx, ai))
	return FALSE ;
       aceInStep (ai,'\t') ;
      if (! aceInInt (ai, &(sx->count)))
	messcrash ("Missing count in TAG_COUNT input file %s, line %d:\n"
		   , sx->inFileName, aceInStreamLine (ai)) ;
      /* name the sequence */
      if (sx->keepName &&  aceInStep (ai,'\t') && (ccp = aceInPos (ai)))
	vtxtPrintf (sx->id, "%s.%d#%d\t%s", sx->prefix, sx->nProcessed + 1, sx->count, ccp) ;
      else
	{
	  vtxtPrintf (sx->id, "%s.%d", sx->prefix, sx->nProcessed + 1) ;
	  vtxtPrintf (sx->id, "#%d", sx->count) ;
	}
      break ;
    }
  sx->nProcessed ++ ;
  sx->tagProcessed += sx->count ;
  if (vtxtPtr (sx->dna))
    {
      int dn ;
      long int ln = strlen(vtxtPtr (sx->dna)) ; /* to cast the product sx->count * ln */
      char *cp = vtxtPtr (sx->dna) ;
      for (dn = 0, cp = vtxtPtr (sx->dna) ; *cp ; cp++)
	switch ((int)*cp)
	  {
	  case '>':
	  case '<': dn++ ; 
	    sx->fastc_paired_end = TRUE ;
	    break ;
	  }
      if (dn >= 0) 
	{ 
	  ln -= dn ; 
	  if (dn >= 2) 
	    { 
	      sx->nProcessed++ ; 
	      sx->tagProcessed += sx->count ; 
	    }
	  else if (sx->fastc_paired_end)
	    { 
	      sx->nProcessed++ ; 
	      sx->tagProcessed += sx->count ; 
	      sx->nRejected++ ; 
	      sx->tagRejected += sx->count ; 
	    }
	}
	
      sx->bpSeqProcessed += ln ;
      sx->bpTagsProcessed += sx->count * ln ;
    }
  sx->entryAdaptorBpClipped = 0 ;
  sx->exitAdaptorBpClipped = 0 ;
  if (sx->count > sx->maxCount) 
    sx->maxCount = sx->count ;
  return TRUE ;
} /* getSequencePart */

/******************/ 
/* possibly merge the 2 parts of a paired en read */
static BOOL getSequence (SX *sx)
{
  if (sx->ai1)
    {
      char *dna1 = 0, *qual1 = 0, *nam1 = 0 ;
      char *dna2 = 0, *qual2 = 0, *nam2 = 0 ;
      AC_HANDLE h = 0 ;
      BOOL ok = FALSE ;
      int line1, line2 ;

      /* register the future line number of the identifiers */
      line1 = aceInStreamLine (sx->ai1) ;
      line2 = aceInStreamLine (sx->ai2) ;

      if (getSequencePart (sx, sx->ai2))
	{
	  h = ac_new_handle () ;
	  dna2 =  strnew (vtxtPtr (sx->dna), h) ;
	  if (vtxtPtr(sx->qual))
	    qual2 =  strnew (vtxtPtr (sx->qual), h) ;
	  if (vtxtPtr (sx->id))
	    nam2 = strnew (vtxtPtr (sx->id), h) ;
	}
      else
	return FALSE ;
      if (getSequencePart (sx, sx->ai1))
	{
	  
	  nam1 = sx->id ? vtxtPtr (sx->id) : 0 ;
	  if (nam1 && ! nam2)
	    messcrash ("Missing identifier %s, seen in file %s line %d, not in file %s line %d"
		       , nam1, sx->inFileName1, line1, sx->inFileName2, line2) ;
	  if (nam2 && ! nam1)
	    messcrash ("Missing identifier %s, seen in file %s line %d, not in file %s line %d"
		       , nam2, sx->inFileName2, line2, sx->inFileName1, line1) ;
	  if (nam1 && nam2)
	    {
	      if (sx->noNameCheck)
		nam2 = nam1 ;
	      if (strcmp (nam1, nam2) &&
		  sx->suffix1 && sx->suffix2)
		{
		  int i = strlen (nam1) - strlen(sx->suffix1) ;
		  int j = strlen (nam2) - strlen(sx->suffix2) ;
		  
		  if (i == j && i > 0 &&
		      !strncmp (nam1, nam2, i) &&
		      !strcmp (nam1+i, sx->suffix1) &&
		      !strcmp (nam2+i, sx->suffix2)
		      )
			nam1[i] = nam2[i] = 0 ; 
		}
	      if (strcmp (nam1, nam2)) 
		{
		  int i, j ;
		  char *cp1, *cp2 ;

		  cp1 = strstr (nam1, "_F3") ; 
		  cp2 = strstr (nam2, "_F5-") ;   /* LIF SOLiD convention */
		  if (cp1 && cp2) 
		    {
		      *cp1-- = *cp2-- = 0 ;
		    }
		  else
		    {
		      cp1 = strstr (nam2, "_F3") ; 
		      cp2 = strstr (nam1, "_F5-") ;   /* LIF SOLiD convention */
		      if (cp1 && cp2) 
			{
			  *cp1-- = *cp2-- = 0 ;
			}
		      else
			{
			  i = strlen (nam1) ; j = strlen (nam2) ;
			  cp1 = nam1 + i - 1 ; cp2 = nam2 + j - 1 ;
			  while (cp1 > nam1 && cp2 > nam2 && *cp1 == *cp2 && i > 0 && j > 0)
			    { i--; j-- ; cp1-- ; cp2-- ; }
			  if (i == j && nam1[i-2] == '/' && !strncmp (nam1, nam2, i - 2))
			    nam1[i-2] = 0 ;
			  else
			    {
			      i = 0 ;
			      cp1 = nam1 ; cp2 = nam2 ;
			      while (*cp1 && *cp1 != ' ')
				{
				  if (*cp1 != *cp2)
				    { i = 0 ; break ; }
				  cp1++ ; cp2++ ; i++ ;
				}
			      if (i > 0)
				*cp1 = *cp2 = 0 ;
			      else
				messcrash ("Non matching identifiers, I quit\n\t%s seen in file %s line %d\n\t%s seen in file %s line %d\n"
					   , nam1, sx->inFileName1, line1, nam2, sx->inFileName2, line2) ;
			    }
			}
		    }
		}
	    }
	  qual1 = sx->qual ? vtxtPtr (sx->qual) : 0 ;
	  if (qual1 && qual2)
	    vtxtPrintf (sx->qual, "><%s", qual2) ;

	  dna1 = sx->dna ? vtxtPtr (sx->dna) : 0 ;
	  if (dna1 && dna2)
	    vtxtPrintf (sx->dna, "><%s", dna2) ;
	  else if (dna2)
	    vtxtPrintf (sx->dna, "<%s", dna2) ;
	  ok = TRUE ;
	}
      ac_free (h) ;
      return ok ;
    }
  else
    return getSequencePart (sx, sx->ai) ;
} /* getSequence */

/*************************************************************************************/

static void complementSequence (SX *sx)
{
  char cc, *dna = vtxtPtr (sx->dna), *cp, *cq ;
  int n = strlen (dna) ;
  
  for (cp = dna, cq = dna + n - 1 ; cp <= cq ; cp++, cq--)
    {
      cc = *cp ;
      *cp = complementLetter (*cq) ;
      *cq = complementLetter (cc) ;
    }
} /* complementSequence */

/*************************************************************************************/
/* if the output is SOLiD, rotate the fasta space letters, else complement */
static void decoySequence (SX *sx)
{
  char *dna = vtxtPtr (sx->dna), cc, *cp ;
  int i, n = strlen (dna) ;
  BOOL solid = FALSE ;
  
  switch (sx->out)
    {
    case CSFASTA:
    case CSFASTC:
    case CCFA:
    case CCFAR:
      solid = TRUE ;
      break ;
    default:
      break ;
    }
  
  for (i = 0, cp = dna ; i < n ; i++, cp++)
    {
      cc = dnaEncodeChar[(int)*cp] ;
      if (solid)
	{
	  switch ((int)cc)
	    {
	    case A_: cc = T_ ; break ;
	    case T_: cc = G_ ; break ;
	    case G_: cc = C_ ; break ;
	    case C_: cc = A_ ; break ;
	    default: break ;
	    }
	}
      else
	{ 
	  switch ((int)cc)
	    {
	    case A_: cc = T_ ; break ;
	    case T_: cc = A_ ; break ;
	    case G_: cc = C_ ; break ;
	    case C_: cc = G_ ; break ;
	    default: break ;
	    }
	}
      *cp = dnaDecodeChar[(int)cc] ;
    }
} /* decoySequence */

/*************************************************************************************/

static BOOL encodeSequence (SX *sx)
{
  DNAFORMAT old = sx->out ;

  if (sx->complement)
    complementSequence (sx) ;
  if (sx->decoy)
    decoySequence (sx) ;
  switch (sx->in)
    {
    case COUNT:
    case POLYMER:
    case LETTER:
    case LETTERPLOT:
    case LETTERSHOW:
    case TENSORFLOW:
      return TRUE ;
    case FASTA:
    case FASTC:
    case FASTQ:
    case TAG_COUNT:
    case RAW:
	switch (sx->out)
	  {
	  case FASTA:
	  case FASTC:
	  case FASTQ:
	  case TAG_COUNT:
	  case RAW:
	  case COUNT:
	  case POLYMER:
	  case LETTER:
	  case LETTERPLOT:
	  case LETTERSHOW:
	  case TENSORFLOW:
	    break ;
	  case CSFASTA:
	  case CSFASTC:
	  case CCFA:
	  case CCFAR:
	  case CTAG_COUNT:
	  case CRAW:
	    fa2ccfa (sx) ;
	    cfa2cfa (sx) ;
	    break ;
	  }
	break ;
    case CTAG_COUNT:
    case CRAW:
    case CCFA:
    case CCFAR: 
    case CSFASTA:
    case CSFASTC:
      switch (sx->out)
	{
	case FASTA:
	case FASTC:
	case FASTQ:
	case RAW:
	case TAG_COUNT:
	case TENSORFLOW:
	  sx->out = CSFASTA ; cfa2cfa (sx) ; sx->out = old ;
	  cfa2fa (sx) ;
	  break ;
	case CSFASTA:
	case CSFASTC:
	case CCFA:
	case CCFAR:
	case CTAG_COUNT:
	case CRAW:
	  cfa2cfa (sx) ;
	  break ;
	case COUNT:
	case POLYMER:
	case LETTER:
	case LETTERPLOT:
	case LETTERSHOW:
	break ;
	}
      break ;
    }
  return TRUE ;
} /* encodeSequence */

/*************************************************************************************/

static void selectOutFile (SX *sx, const char *ccp)
{
  char cc ;
  int n = 0, ii ;
  
  for (ii = 0 ; ii < sx->splitByPrefix ; ii++)
    {
      n <<= 2 ;
      cc = ccp[ii] ;
      switch ((int)cc)
	{
	case '1': case 'c':	case 'C': n+=1 ; break ;
	case '2': case 'g':	case 'G': n+=2 ; break ;
	case '3': case 't':	case 'T': n+=3 ; break ;
	}
    }
  if (sx->in == FASTC && (sx->out == FASTA || sx->out == FASTQ))
     {
       sx->ao = 0 ;
       sx->ao1 = sx->aos1[n] ;
       sx->ao2 = sx->aos2[n] ;
     }
   else
     {
       sx->ao = sx->aos[n] ;
       sx->ao1 = 0 ;
       sx->ao2 = 0 ;
     }
  return ;
} /* selectOutFile */

/*************************************************************************************/
/* Cumul the histogram of the polymers */
static void sxPolymerAccumulate (SX *sx)
{
  long int *b = sx->polymerCount ;
  long int *bb = sx->polymerProfile ;
  const char *ccp ;
  int n, cc ;
  BOOL isForward = TRUE ;

  for (ccp = vtxtPtr (sx->dna) ; *ccp ; ccp++)
    {
      cc = *ccp ;
      if (cc == '>') continue ;
      if (cc == '<') { if (0) isForward = FALSE ; continue ;}
      for (n = 0 ; *ccp == cc ; ccp++)
	n++ ;
      ccp-- ;

      /* we have localized n repeats of base cc */
      cc = dnaEncodeChar[cc] ;
      if (!  isForward)
	cc = complementBase(cc) ;
      switch (cc)
	{
	case A_: cc = 0 ; break ;
	case T_: cc = 1 ; break ;
	case G_: cc = 2 ; break ;
	case C_: cc = 3 ; break ;
	default: continue ;
	}
      b[cc]+= n ;
      if (n > sx->polymer)
	n = sx->polymer ;
      bb[4*n + cc]++ ;
    }
  return ;
} /* sxPolymerAccumulate */

/**************************/

static void sxPolymerReport (SX *sx)
{
  long int *b = sx->polymerCount ;
  long int *bb = sx->polymerProfile ;

  if (0)
    {
      AC_HANDLE h = ac_new_handle () ;
      ACEOUT ao = aceOutCreate (sx->outFileName, ".polymer.histo", FALSE, h) ;
      aceOutDate (ao, "###", sx->title) ;
      
      aceOutf (ao, "# Base\tA\tT\tG\tC") ;
      aceOutf (ao, "\n# Cumul") ;
      for (int i = 0 ; i < 4 ; i++)
	aceOutf (ao, "\t%ld", b[i]) ;
      aceOutf (ao, "\n# Length\tpolyA\tpolyT\tpolyG\tpolyC") ;
      for (int n = 1 ; n <= sx->polymer && n < POLYMER_MAX ; n++)
	{
	  aceOutf (ao, "\n%d", n) ;
	  for (int i = 0 ; i < 4 ; i++)
	    aceOutf (ao, "\t%ld", bb[4*n + i]) ;
	}
      aceOut (ao, "\n\n") ;
      ac_free (h) ;
    }
  if (1)
    {
      AC_HANDLE h = ac_new_handle () ;
      ACEOUT ao = aceOutCreate (sx->outFileName, ".polymer.tsf", FALSE, h) ;
      aceOutDate (ao, "###", sx->title) ;
      const char *run = sx->runTitle ? sx->runTitle : "X" ;
      char base[4] = "ATGC" ;
      int nMax = sx->polymer ;
      for (int i = 0 ; i < 4 ; i++)
	for (int n = 1 ; n <= nMax + 4 ; n++)
	  if (bb[4*n+i] >  0)
	    aceOutf (ao, "%s\t%c%d\ti\t%ld\n", run, base[i], n, bb[4*n + i]) ;
      ac_free (h) ;
    }
} /* sxPolymerReport */

/*************************************************************************************/
/* Cumul the bp according to their position */
static void sxLetterAccumulate (SX *sx)
{
  int di, i, i1, cc, count = sx->count ;
  long int *b = sx->letterCount ;
  long int *bb = sx->letterProfile ;
  const char *ccp = vtxtPtr (sx->dna) ;

  for (i = i1 = di = 0, ccp = vtxtPtr (sx->dna) ; *ccp && i < 1 * (sx->letterNN + 4) ; i++, ccp++)
    {
      cc = ace_lower(*ccp) ;
      if (cc == '>') { if (i1) sx->letterLength [i1] += count ; i1 = 0 ; continue ; }
      if (cc == '<') { if (i1) sx->letterLength [i1] += count ; i1 = 0 ; di = sx->letterNN + 2 ; continue ; }
      b[cc]+= count ;
      bb[256*(i1+di) + cc] += count ;
      i1++ ;
    }
  if (i1) sx->letterLength [di + i1] += count ; 
  return ;
} /* sxLetterAccumulate */

/*************************************************************************************/
/* Cumul the quality according to their position */
static void sxQualityAccumulate (SX *sx)
{
  int i, cc, count = sx->count ;
  long int *b = sx->letterCount ;
  long int *bb = sx->letterProfile ;
  const char *ccp = vtxtPtr (sx->dna) ;

  for (i = 0, ccp = vtxtPtr (sx->qual) ; ccp && *ccp && i < sx->letterNN ; i++, ccp++)
    {
      cc = *ccp ; b[cc]+= count ;
      bb[256*i + cc] += count ;
    }
  sx->letterLength [i] += count ;
  return ;
} /* sxQualityAccumulate */

/*************************************************************************************/
/* dedicated code to export a modified genome with periodic SNPs */
static void makeTestGenome (SX *sx)
{	 
  AC_HANDLE h = ac_new_handle () ;
  int nerr = 0, ln = strlen (vtxtPtr (sx->dna)) ;
  char *dna2 = halloc (2*ln + 1000, h) ;   /* make lots of room */
  int mult, sub = 0, period = sx->makeTestPeriod ;
  char *atgc = "ATGC" ;
  char *TrTv[256] ;
  char *seq = vtxtPtr (sx->id) ;
  BOOL isRandom = sx->makeTestRandom ;
  BOOL multiLn = TRUE ;
  int subRate = sx->makeTestSubRate ;
  ACEOUT ao =  sx->aoTestSnps ;
  static int nBB = 0, nMM = 0, nSUB = 0, nINS = 0, nDEL = 0, nPLUS = 0, nMINUS = 0, nOK = 0 ; 
  static double alpha = 1.0 ;
  const int shift = sx->makeTestShift ;

  period = period  - 2 ; /* -2 because we add 2 exact bases between each event */
  if (period < 2) period = 2 ;

  TrTv['A'] = "GGGGCT" ;  /* twice as many transition than transversion */ 
  TrTv['G'] = "AAAACT" ;
  TrTv['T'] = "CCCCAG" ;
  TrTv['C'] = "TTTTAG" ;


  /* edit the original sequence */
  memcpy (dna2, vtxtPtr (sx->dna), ln+1) ;  /* include the termnal zero */

  if (1) /* create a sub every period */
    {
      int ii = 0, jj = 0, dx = 0 ;
      char cc, *cp = vtxtPtr (sx->dna), *cq = dna2 ;
      char *cqMax = cq + 2*ln ;
      int i, t = 0 ;
      float type  ;
      BOOL ok ;
      int nbb = nBB ;
      cp += shift ;
      cq += shift ;
      for (ii = jj = shift ; ii < ln - 10 ; ii++, cp++, cq++)
	{
	  int targetNmm = nbb / sx->makeTestPeriod ;

	  if (cq > cqMax)
	    messcrash ("ii=%d, jj=%d",ii,jj);
	  *cq = *cp ;   /* copy the base */
	  if (ii < 10)
	    continue ;

	  nbb = nBB + ii ;
	  if (nbb % 1000000 == 0)
	    fprintf (stderr, "bb:%d t:%d m:%d s:%d i:%d d:%d p:%d +:%d -:%d ok:%d alpha=%g\n"
		     , nbb, targetNmm, nMM,  nSUB, nINS, nDEL, period
		     , nPLUS, nMINUS, nOK, alpha
		     ) ; 
	  if (isRandom)
	    {   /* fine adust the error rate */
	      float dz = 2.0/sqrt(nbb +10 ) ;

	      if (nMM < targetNmm * (1-dz))
		{
		  if (0.9 * alpha * period * randfloat() > 1.0)
		    continue ;
		  if (nMINUS == 100000)
		    fprintf (stderr, "#### bb:%d t:%d m:%d s:%d i:%d d:%d p:%d +:%d -:%d ok:%d dz=%g\n"
		     , nbb, targetNmm, nMM,  nSUB, nINS, nDEL, period
			     , nPLUS, nMINUS, nOK, dz
			     ) ;		
		  alpha *= 0.999999 ;
		  nMINUS++ ;
		}
	      else if (nMM > targetNmm * (1 + dz))
		{
		  if (1.1 * alpha * period * randfloat() > 1.0)
		    continue ;
		  if (nPLUS == 100000)
		    fprintf (stderr, "#### bb:%d t:%d m:%d s:%d i:%d d:%d p:%d +:%d -:%d ok:%d dz=%g\n"
			     , nbb, targetNmm, nMM,  nSUB, nINS, nDEL, period
			     , nPLUS, nMINUS, nOK, dz
		     ) ; 
		  alpha *= 1.000001 ;
		  nPLUS++ ;
		}
	      else
		{
		  if (alpha * period * randfloat() > 1.0 )
		    continue ;
		  nOK++ ;
		}
	    }
	  else
	    if ((ii - shift) % period)
	      continue ;

	  /* select multiplicity 1,2 or 3 in ratio 85:13:2 */
	  if (multiLn)
	    {
	      type = 100 * randfloat () ;
	      if (type <= 2)
		mult = 3 ;
	      else if (type <= 13)
		mult = 2 ;
	      else
		mult = 1 ;
	    }
	  else
	    mult = 1 ;
	  
	  nMM += mult ;
	    
	  /* select sub insert or delete */
	  type = 100 * randfloat () ;
	  ok = FALSE ;
	  if (type <= subRate)
	    { /* create a set of subs favoring transitions over transversions */
	      char *trtv ;
	      
	      cq-- ; cp-- ; ii-- ; jj-- ;
	      for (i = 0 ; i < mult ; i++)
		{
		  *++cq = *++cp ; ii++ ; jj++ ; /* copy the original base */
		  cc = ace_upper(*cq) ;
		  trtv = TrTv[((int)cc) & 0xff] ;
		  
		  sub = 6 * randfloat () ;
		  *cq = trtv[sub%6] ;
		  nerr++ ; 
		  aceOutf (ao, "%s\t%d\t%c>%c\n", seq, ii + 1 + dx, cc, *cq) ;
		}
	      nSUB += mult ;
	      ok = TRUE ;
	      for (i = 0 ; i < 3 ; i++)  /* minimal distance between events */
		{ *++cq = *++cp ; ii++ ; jj++ ; }
	      continue ;
	    }
	  else if (nDEL < nINS) /* favor delete, because they are rejected if the position is sliding */
	    {
	      t = mult ;
	      if (   /* not sliding or next position is not sliding */
		  (cp[0] != cp[t] &&  cp[-1] != cp[t-1]) ||
		  (cp[1] != cp[t+1] &&  cp[0] != cp[t]) 
		  )
		{ /* do not delete a sliding region */
		  while ((cp[0] == cp[t] || cp[-1] == cp[t-1]) )
		    { *cq++ = *cp++ ; ii++ ; jj++ ; }
		  if (ii >= ln)
		    { *cq = 0 ; break ; }
		  aceOutf (ao, "%s\t%d\tDel", seq, ii + 1 + dx) ;
		  for (i = 0 ; i < t ; i++)   /* forget t non sliding bases */
		    { aceOutf (ao, "%c", ace_upper(*cp++)) ; ii++ ; }
		  aceOutf (ao, "\n") ;
		  if (ii >= ln)
		    { *cq = 0 ; break ; }
		  *cq = *cp ;
		  dx -= mult ;
		  nDEL++ ;
		  ok = TRUE ;
		  for (i = 0 ; i < 1 ; i++)  /* minimal distance between events */
		    { *++cq = *++cp ; ii++ ; jj++ ; }
		}
	    }
	  if (! ok) /* insert */
	    {
	      t = mult ;
	      for (i = 0 ; i < t ; i++)
		{
		  int r = randint () % 4  ;
		  *cq = atgc[r % 4] ; /* insert t bases */
		  cq++ ; jj++ ;
		}
	      while (ace_upper(cq[-t]) == ace_upper(cp[0]) || 
		     ace_upper(cq[-t]) == ace_upper(cq[-2*t])
		     )
		cq[-t] = atgc[randint () % 4] ; /* prevent sliding right */
	      while (ace_upper(cq[-1]) == ace_upper(cp[-t -1]))
		cq[-1] = atgc[randint () % 4] ; /* prevent sliding left */
	      *cq = 0 ;
	      aceOutf (ao, "%s\t%d\tIns%s\n", seq, ii + 1 + dx, cq - t) ;		     
	      *cq = *cp ;
	      nerr++ ;
	      dx += mult ;
	      nINS++ ;
	      for (i = 0 ; i < 2 ; i++)  /* minimal distance between events */
		{ *++cq = *++cp ; ii++ ; jj++ ; }
	    }
	  
	}

      /* copy the remainder of the sequence */
      if (ii < ln)
	while ((*cq++ = *cp++))
	  if (cq > cqMax)
	    messcrash ("ii=%d, jj=%d",ii,jj);
      if (strlen (dna2) > 2* ln)
	messcrash ("ln=%d, ln2=%d",ln,strlen (dna2));
    }

  nBB += ln ; /* static cumul of all processed sequences */

  /* export the fasta file */
  if (1)
    {
      char cc, *cq, *cp = dna2 ;
      int ln = strlen(dna2) ;
      ao = sx->aoTestSeqs ;

      aceOutf (ao, ">%s\n", seq) ;
      while (ln > 0)
	{
	  int max = sx->maxLineLn ;
	  if (max == 0) max = 60 ;
	  if (max > ln) max = ln ;
	  cq = cp + max ;
	  cc = *cq ;
	  *cq = 0 ;
	  aceOutf (ao, "%s\n", cp) ;
	  cp = cq ;
	  *cp = cc ;
	  ln -= max ;
	}
   }

  ac_free (h) ; 
  return ;
} /* makeTestGenome */ 

/*************************************************************************************/
/* dedicated code to export modified n-mers and test the aligner code */
static void makeTestWithSNP (SX *sx, BOOL isExact, BOOL forward)
{	 
  AC_HANDLE h = ac_new_handle () ;
  int LN =  sx->makeTestLength ;
  int LN2 =  sx->makeTestPairLength ?  sx->makeTestPairLength : sx->makeTestLength ;
  int period = sx->makeTestPeriod ;
  char cc = 0, cc2, *cp, *cq, namBuf[4000], eNamBuf[1024], gNamBuf[1024], sNamBuf[1024], *atgc = "ATGC", *namType ;
  char eBuf[LN2+1], qBuf[LN2+1], buf[LN2+1] ;
  char eBuf2[LN2+1], qBuf2[LN2+1], buf2[LN2+1] ;
  static int nSeq[4], firstPass = 0 ;
  int ii, jj, k, gene ;
  int MX = sx->makeTest ;
  int ln, maxStart ;
  int makeTestStep= sx->makeTestStep ;
  static DICT *dict = 0, *dicts[4] ;
  ACEOUT ao  ;
  int type =  (isExact ? 1 : 0) + 2 * (forward ? 1 : 0) ;
  const int shift = sx->makeTestShift ;
  vTXT dna2 = vtxtHandleCreate (h) ;
  char *oformat = "fasta" ;

  switch (sx->out)
    {
    case FASTQ: oformat = "fastq" ; break ;
    case RAW: oformat = "raw" ; break ;
    default : oformat = "fasta" ; break ;
    }
  if (firstPass == 0)
    {
      memset (dicts, 0, sizeof(dicts)) ;
      memset (nSeq, 0, sizeof(nSeq)) ;
      firstPass = 1 ;
    }

  if (nSeq[type] > MX) goto done ;

  dict = dicts[type] ;
  if (! dict)
    dict = dicts[type] = dictCreate (10000) ;

  /* register the gene name, so we export a single transcript per gene */
  if (dict)
    {
      char *cp, *cq ;
      sprintf(eNamBuf, "%s",  vtxtPtr (sx->id)) ;
      cp = strstr (eNamBuf, "|Gene|") ;
      if (!cp) cp = strstr (eNamBuf, "|GENE|") ;
      if (cp)
	cp += 6 ;
      else
	cp = eNamBuf ;
      if (cp && *cp)
	{
	  cq = strstr (cp, "|") ; 
	  if (cq) *cq = 0 ;
	  if (dictFind (dict, cp, &gene))  /* gene already exported */
	    return ;
	  dictAdd (dict, cp, &gene) ;
	}
      sprintf(eNamBuf, "%s",  vtxtPtr (sx->id)) ;
      gNamBuf[0] = 0 ;
      cp = strstr (eNamBuf, "|Gene|") ;
      if (cp) {  memcpy(gNamBuf, cp, sizeof (gNamBuf)) ; *cp = 0 ;}
      memcpy(sNamBuf, eNamBuf, sizeof (sNamBuf) ) ;
    }
  if (isExact && forward)
    ao = aceOutCreate (sx->outFileName, messprintf (".%s.exact.forward.%s", dictName (dict, gene), oformat), sx->gzo, h) ;
  else if (! isExact && forward)
    ao = aceOutCreate (sx->outFileName, messprintf (".%s.variant.forward.%s", dictName (dict, gene), oformat), sx->gzo, h) ;
  else if (isExact && ! forward)
    ao = aceOutCreate (sx->outFileName, messprintf (".%s.exact.reverse.%s", dictName (dict, gene), oformat), sx->gzo, h) ;
  else
    ao = aceOutCreate (sx->outFileName, messprintf (".%s.variant.reverse.%s", dictName (dict, gene), oformat), sx->gzo, h) ;

  if (type == 0 && sx->makeTestAddPolyA)
    for (ii = 0 ; ii < LN - 30 ; ii++) /* always keep at least 30 alignable letters */
      vtxtPrint (sx->dna, "A") ;
  ln = strlen (vtxtPtr (sx->dna)) ;

  maxStart = ln - LN2  ; /* 2011_05_11, was 20.0, but now we want a non random positioning */
  /*
    cp = vtxtPtr (sx->id) ;
    cp += strlen(cp) - 2 ;
    if (*(cp+1) != '1' && *cp == '_')
    return ;
  */
  buf[LN2] = 0 ;  eBuf[LN2] = 0 ;
  memset (qBuf, 'i', LN2) ; qBuf[LN2] = 0 ;
  
  /* edit the original sequence */
  vtxtPrintf (dna2, vtxtPtr (sx->dna)) ;  /* holds the edited DNA */
  if (! sx->makeTestInDel)
    {
      int m ;
      namType = "sub" ;
      for (k = 0, m = shift + period/2 ; m < ln ; k++, m += period)
	{
	  cp = vtxtAt (dna2, m) ;
	  cc = *cp ;
	  cc2 = atgc[k % 4] ;
	  if (ace_upper(cc) == ace_upper(cc2))  cc2 = atgc[(k+1) % 4] ;  
	  *cp = ace_upper(cc2) ;
	}
    }
  else
    {  
      int m ;
      namType = "indel" ;
      for (k = 0, m = shift + period/2 ; m < ln ; k++, m += period)
	{
	  switch (k % 2)
	    {
	    case 0: /* insertion */
	      cp = vtxtAt (dna2, m) ;
	      cq = vtxtAt (dna2, 0) + ln ;
	      while (cq > cp)
		{ *cq = *(cq-1) ; cq-- ; }
	      *cq =  atgc[(m + (m>>1) + (m >> 3) + (m>>4)) % 4] ;
	      break ;
	    case 1: /* deletion */
	      cp = vtxtAt (dna2, m) ;
	      while ((*cp = *(cp+1))) cp++ ;
	      break ;
	    }
	}
    }
     

  /* export n-fold coverage */
  for (ii = jj = shift ; ii < maxStart ; ii+=makeTestStep )
    {
      if (++(nSeq[type]) > MX) 
	goto done ;
      strncpy (eBuf, vtxtAt (sx->dna, ii), LN2) ;
      strncpy (buf, vtxtAt (dna2, ii), LN2) ;
      sprintf (namBuf, "%s@%d%s",  sNamBuf, ii, gNamBuf) ;

      jj++ ;
      if (! forward)
	{
	  int i ;
	  char *cp, *cq ;
	  for (i = 0, cp = buf2 + LN2 - 1, cq = buf ; i < LN2 ; i++, cp--, cq++) 
	    *cp =  complementLetter (*cq) ;
	  buf2[LN2] = 0 ;
	  for (i = 0, cp = eBuf2 + LN2 - 1, cq = eBuf ; i < LN2 ; i++, cp--, cq++) 
	    *cp = complementLetter (*cq) ;
	  eBuf2[LN2] = 0 ;
	  for (i = 0, cp = qBuf2 + LN2 - 1, cq = qBuf ; cq < cp ; i++, cp--, cq++) 
	    *cp = *cq ; 
	  qBuf2[LN2] = 0 ;
	}
      else
	{
	  strcpy (buf2, buf) ;
	  strcpy (eBuf2, eBuf) ;
	  strcpy (qBuf2, qBuf) ;
	}
      if (isExact)
	{
	  if (sx->out == FASTQ)
	    {
	      aceOutf (ao, "@%s.exact  //%d\n%s\n", namBuf, nSeq[type], eBuf2) ;
	      aceOutf (ao, "+%s.exact   //%d\n%s\n", namBuf, nSeq[type], qBuf2) ;
	    }	 
	   else if (sx->out == RAW)
	    {
	      aceOutf (ao, "%s\t%s.exact\n", eBuf2, namBuf, namType) ;
	    }
	  else
	    {
	      aceOutf (ao, ">%s.exact   //%d\n", namBuf, nSeq[type]) ;
	      if (sx->makeTestPairLength) 
		{
		  int i ;
		  char *cp, *cq, cc, bufP[LN+1] ;
		  
		  cp = eBuf2 + LN ; cc = *cp ; *cp = 0 ;
		  aceOutf (ao, "%s><",  eBuf2) ;   /* export the first read */
		  *cp = cc ;

		  for (i = 0, cp = bufP, cq = eBuf2 + LN2 - 1 ; i < LN ; cp++, cq--, i++)
		    *cp = complementLetter (*cq) ;
		  *cp = 0 ;
		  aceOutf (ao, "%s\n",  bufP) ;
		}
	      else
		aceOutf (ao, "%s\n",  eBuf2) ;
	    }
	}	    
      else
	{
	   if (sx->out == FASTQ)
	    {
	      aceOutf (ao, "@%s.%s  //%d\n%s\n", namBuf, namType, nSeq[type], buf2) ;
	      aceOutf (ao, "+%s.%s  //%d\n%s\n", namBuf, namType, nSeq[type], qBuf2) ;
	    }
	   else if (sx->out == RAW)
	    {
	      aceOutf (ao, "%s\t%s.%s\n", buf2, namBuf, namType) ;
	    }
	   else
	     {
	       aceOutf (ao, ">%s.%s  //%d\n", namBuf, namType, nSeq[type]) ;
	       if (sx->makeTestPairLength) 
		 {
		   int i ;
		   char *cp, *cq, cc, bufP[LN+1] ;
		   
		   cp = buf2 + LN ; cc = *cp ; *cp = 0 ;
		   aceOutf (ao, "%s><",  buf2) ;   /* export the first read */
		   *cp = cc ;
		   
		   for (i = 0, cp = bufP, cq = buf2 + LN2 - 1 ; i < LN ; cp++, cq--, i++)
		     *cp = complementLetter (*cq) ;
		   *cp = 0 ;
		   aceOutf (ao, "%s\n",  bufP) ;
		 }
	       else
		 aceOutf (ao, "%s\n", buf2) ;
	     }
 	}
    }
      
 done:
  ac_free (ao) ;
  ac_free (h) ; 
  return ;
} /* makeTestWithSNP */

/*************************************************************************************/
/*************************************************************************************/

static void exportTm (SX *sx)
{
  char *dna = vtxtPtr (sx->dna) ; 
  int entropy, length = dna ? strlen (dna) : 0 ;
  float Tm = 0, GC = 0 ;
  Array a = sx->encodedDna ;

  array (a, length + 1, char) = 0 ; /* make space */
  memcpy (a->base, dna, length) ;
  dnaEncodeArray (a) ;

  entropy = oligoEntropy (arrp(a, 0, unsigned char), length, 0) ;
  Tm = oligoTm (a, 1, length, &GC) ;

  aceOutf (sx->aoTm, "%s\t%d\t%d\t%.2f\t%.2f\n"
	   , vtxtPtr (sx->id)
	   , length
	   , entropy
	   , Tm
	   , GC
	   ) ;

  return ;
} /* exportTm */

/*************************************************************************************/

static void exportSequence (SX *sx)
{
  int n, njump = 0 ;
  ACEOUT ao, ao1, ao2 ;
  SEQ *sq ;
  char *cp, cc = '>' ;
  char runTitle[1024] ;

  runTitle[0] = 0 ;
  if (sx->out == FASTQ && sx->runTitle)
    sprintf(runTitle, "%s.", sx->runTitle) ;
  if (sx->splitByPrefix)
    selectOutFile (sx, vtxtPtr (sx->dna)) ;
  ao = sx->ao ;
  ao1 = sx->ao1 ;
  ao2 = sx->ao2 ;
  
  if (sx->upper)
    {
      char *cp =  vtxtPtr (sx->dna) ;
      if (cp--)
	while (*++cp) *cp = ace_upper(*cp) ;
    }

   if (sx->makeTest)
     {
       makeTestWithSNP (sx, FALSE, TRUE) ;
       makeTestWithSNP (sx, TRUE, TRUE) ;
       if (sx->makeTestReverse)
	 {
	   makeTestWithSNP (sx, FALSE, FALSE) ;
	   makeTestWithSNP (sx, TRUE, FALSE) ;
	 }
     }
   else if (sx->makeTestGenome)
     {
       makeTestGenome (sx) ;
     }
  else switch (sx->out)
    {
    case COUNT:
      break ;
    case CSFASTA:
    case CSFASTC:
    case CCFA:
    case CCFAR:
      if (sx->jumpAnchor) njump = 1 ;
    case FASTA:
    case FASTC:  
    case FASTQ:
      cp = 0 ;
      if (sx->aoId && (cp = strchr(vtxtPtr (sx->id), '\t')))
	{
	  aceOutf (sx->aoId, "%s\n", vtxtPtr (sx->id)) ;
	  *cp = 0 ;
	}
      
      cc = (sx->out == FASTQ ? '@' : '>') ;
      if (ao)
	aceOutf (ao, "%c%s%s\n", cc, runTitle, vtxtPtr (sx->id)) ;
      if (ao1)
	aceOutf (ao1, "%c%s%s\n", cc, runTitle, vtxtPtr (sx->id)) ;
      if (ao2)
	aceOutf (ao2, "%c%s%s\n", cc,  runTitle, vtxtPtr (sx->id)) ;
    
      if (cp)
	*cp = '\t' ;

      if (sx->maxLineLn > 0)
	{
	  char cc = 0, *cp = vtxtPtr (sx->dna) + njump ;
	  int n = strlen (cp) ;

	  while (n > 0)
	    {
	      if (n > sx->maxLineLn)
		{
		  cc = *(cp + sx->maxLineLn) ;
		  *(cp + sx->maxLineLn) = 0 ;
		}
	      aceOutf (ao, "%s\n", cp) ;
	      if (n > sx->maxLineLn)
		{
		  *(cp + sx->maxLineLn) = cc ;
		}
	      cp += sx->maxLineLn ;
	      n -= sx->maxLineLn ;
	    }
	}
      else
	{
	  int n1 = 0, n2 = 0 ;
	  if (ao)
	    {
	      aceOutf (ao, "%s\n", vtxtPtr (sx->dna)) ;
	      if (sx->out == FASTQ)
		 aceOutf (ao, "+%s%s\n%s\n",  runTitle, vtxtPtr (sx->id), vtxtPtr (sx->qual)) ;
	    }
	  if (ao1)
	    {
	      cp = strstr (vtxtPtr (sx->dna), "><") ;
	      if (cp) 
		{
		  *cp = 0 ;
		  aceOutf (ao1, "%s\n", vtxtPtr (sx->dna)) ;
		  n1 = strlen (vtxtPtr (sx->dna)) ;
		  *cp = '>' ;
		}
	      else 
		{
		  cp = strstr (vtxtPtr (sx->dna), "<") ;
		  if (cp) 
		    {
		      aceOutf (ao1, "n\n") ; /* forward read is missing, export something */ 
		      n1 = 1 ;
		    }
		  else
		    {
		      aceOutf (ao1, "%s\n", vtxtPtr (sx->dna)) ;
		      n1 = strlen (vtxtPtr (sx->dna)) ;
		    }
		}
	      if (sx->out == FASTQ)
		{
		  char cg, *cq =  vtxtPtr (sx->qual) + n1 ;
		  if (n1 > strlen (vtxtPtr (sx->qual)))
		    messcrash ("%s sequence length %d longer than quality table", vtxtPtr (sx->id), n1) ;
		  cg = *cq ; *cq = 0 ;
		  aceOutf (ao1, "+%s%s\n%s\n", runTitle, vtxtPtr (sx->id), vtxtPtr (sx->qual)) ;
		  *cq = cg ;
		}
	    }
	  if (ao2)
	    {
	      cp = strstr (vtxtPtr (sx->dna), "<") ;
	      if (cp && *(cp+1)) 
		{
		  aceOutf (ao2, "%s\n", cp + 1) ;
		  n2 = strlen (cp + 1) ;
		}
	      else
		{
		  aceOutf (ao2, "n\n") ; /* reverse read is missing, export something */
		  n2 = 1 ;
		}
	      if (sx->out == FASTQ)
		{
		  char cg, *cq =  vtxtPtr (sx->qual2) + n2 ;
		  if (n2 > strlen (vtxtPtr (sx->qual2)))
		    messcrash ("%s sequence length %d longer than quality table", vtxtPtr (sx->id), n2) ;
		  cg = *cq ; *cq = 0 ;
		  aceOutf (ao2, "+%s%s\n%s\n", runTitle, vtxtPtr (sx->id), vtxtPtr (sx->qual2)) ;
		  *cq = cg ;
		}
	    }
	}
      break ;
    case RAW:
    case CRAW:
      for (n = 0 ; n < sx->count ; n++)
	{
	  int n1 = 0 ; /* sx->splitLongDna ; */

	  if (! n1)
	    {
	      if (sx->keepName)
		aceOutf (ao, "%s\t%s\n", vtxtPtr (sx->dna), vtxtPtr (sx->id)) ;
	      else
		aceOutf (ao, "%s\n", vtxtPtr (sx->dna)) ;
	    }
	  else
	    messcrash ("splitLongDna not programmed sorry") ;
	}
      break ;
    case TENSORFLOW:
      for (n = 0 ; n < sx->count ; n++)
	{
	  int i = sx->rightClipAt ;
	  const char *ccp ;
	  char cc = 0 ;
	  aceOutf (ao, "%s,[", vtxtPtr (sx->id)) ;
	  for (ccp = vtxtPtr (sx->dna) ; i > 0 ; i-- , ccp += (*ccp ? 1 : 0), cc = ',')
	    {
	      if (cc) aceOut (ao, ",") ;
	      switch ((int)(*ccp))
		{
		case 'a': aceOut (ao, "[1,0,0,0]") ; break ;
		case 'c': aceOut (ao, "[0,1,0,0]") ; break ;
		case 'g': aceOut (ao, "[0,0,1,0]") ; break ;
		case 't': aceOut (ao, "[0,0,0,1]") ; break ;
		default:  aceOut (ao, "[0,0,0,0]") ; break ;
		}
	    }
	  aceOut (ao, "]\n") ;
	}
      break ;
    case TAG_COUNT:
    case CTAG_COUNT:
      if (0 && sx->nExported > 7647270) fprintf(stderr, "%d:: %s\t", sx->nExported, vtxtPtr (sx->dna)) ;
      dictAdd (sx->dict, vtxtPtr (sx->dna), &n) ;
      if (0 && sx->nExported > 7647270) fprintf(stderr, "%d:: %s\n", sx->nExported, dictName(sx->dict, n)) ;
      if (0 && strncmp (dictName(sx->dict, n), "CTA",3))
	  messcrash ("Got %s rather than %s line %d", dictName(sx->dict, n), vtxtPtr (sx->dna), sx->nExported) ;
      sq = arrayp (sx->tagArray, n - 1, SEQ) ;
      sq->dna = n ; /* needed since we shall sort */
      sq->count += sx->count ;
      sq->dnaLn = sx->dnaLn ;
      if (sx->keepName)
	{
	  if (! sq->id)
	    {
	      sq->id = vtxtCreate () ; 
	      vtxtPrint (sq->id, vtxtPtr (sx->id)) ;
	    }
	  else
	    vtxtPrintf(sq->id, "\t%s", vtxtPtr (sx->id)) ;
	}
      break ;
    case POLYMER:
      sxPolymerAccumulate (sx) ;
      break ;
    case LETTER:
    case LETTERPLOT:
    case LETTERSHOW:
      if (sx->dna && vtxtPtr (sx->dna)) 
	{
	  if (sx->qualityPlot)
	    sxQualityAccumulate (sx) ;
	  else
	    sxLetterAccumulate (sx) ;
	}
      break ;
    }
  sx->nExported++ ;
} /* exportSequence */

/*************************************************************************************/
/*************************************************************************************/

int shadowOrder (const void *a, const void *b)
{
  const SHADOW *up = (const SHADOW *) a, *vp = (const SHADOW *) b ;
  int n ;
  
  n = up->target - vp->target ; if (n) return n ;
  n = up->gene - vp->gene ; if (n) return n ;
  n = up->Parent - vp->Parent ; if (n) return n ;
  if (up->mrna && !vp->mrna) return 1 ;
  if (!up->mrna && vp->mrna) return -1 ;
  n = up->mrna - vp->mrna ; if (n) return n ;
  n = up->cds - vp->cds ; if (n) return n ;
  n = up->x1 - vp->x1 ; if (n) return n ;
  n = up->x2 - vp->x2 ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  n = up->ID - vp->ID ; if (n) return n ;


  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/

int shadowParentOrder (const void *a, const void *b)
{
  const SHADOW *up = (const SHADOW *) a, *vp = (const SHADOW *) b ;
  int n ;
  
  n = up->Parent - vp->Parent ; if (n) return n ;
  if (up->mrna && !vp->mrna) return 1 ;
  if (!up->mrna && vp->mrna) return -1 ;

  return 0 ;
}

/*************************************************************************************/

static void showShadow (SX *sx, DICT *dict, Array shadows)
{
  int n ;
  SHADOW *up ;

  for (n = 0 ; n < arrayMax (shadows) ; n++)
    {
      up = arrp (shadows, n, SHADOW) ;
      fprintf (stderr, "##### %d id=%s parent=%s", n, up->ID ? dictName(dict, up->ID) : "-", up->Parent ? dictName(dict, up->Parent) : "-") ;
      fprintf (stderr, " target=%s %d %d", up->target ? dictName(sx->selectDict, up->target) : "-", up->a1, up->a2) ;
      fprintf (stderr, " x1/x2 %d %d", up->x1, up->x2) ;
      fprintf (stderr, " gene:%s:%d", up->gene ?  dictName(dict, up->gene) : "-" , up->gene) ;
      fprintf (stderr, " mrna:%s:%d", up->mrna ?  dictName(dict, up->mrna) : "-" , up->mrna) ;
      fprintf (stderr, up->cds ? " cds" : " exon") ;      
      fprintf (stderr, "\n") ;
    }
} /* showShadow */

/*************************************************************************************/
/* gather the geometry off the shadow file */
static int parseShadowFile (SX *sx, int type)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = 0 ;
  int nn = 0 ;
  const char *ccp ;
  int a1, a2, x1, x2, oldX2 = 0, mrna = 0, oldMrna = 0, gene = 0, target = 0 ;
  char cutter ;
  SHADOW *up, *vp ;

  if (type == 1)
    ai = aceInCreate (sx->shadowFileName, FALSE, h) ;
  else
    ai = aceInCreate (sx->spongeFileName, FALSE, h) ;
  if (!ai)
    goto done ;

  sx->selectDict = dictHandleCreate (10000, sx->h) ;
  sx->shadowDict = dictHandleCreate (10000, sx->h) ;
  sx->shadowArray = arrayHandleCreate (10000, SHADOW, sx->h) ;
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      ccp = aceInWordCut (ai, "\t", &cutter) ;
      if (! ccp || *ccp == '#')
	continue ;
      a1 = a2 = x1 = x2 = 0 ;
      dictAdd (sx->shadowDict, ccp, &mrna) ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &x1))
	continue ;
      if (type == 1) { aceInStep (ai, '\t') ; aceInInt (ai, &x2) ;  }
      else { x2 = x1 + 1 ; }
      aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ;
      if (! ccp)
	continue ;
      dictAdd (sx->selectDict, ccp, &target) ;
      if (strchr (ccp, ' '))
	sx->selectDictHasSpace = TRUE ;
      aceInStep (ai, '\t') ; aceInInt (ai, &a1) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &a2) ;

      if (a1 < 1 || a2 < 1 || x1 < 1 || x2 < 1)
	messcrash ("missing or null or negative coordinates at lane %d in shadow file %s", 
		   aceInStreamLine (ai)
		   , type ==1 ? sx->shadowFileName : sx->spongeFileName
		   ) ;
      if (type == 2)
	{ 
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ;
	  if (ccp) dictAdd (sx->shadowDict, ccp, &gene) ;
	  if (mrna != oldMrna)
	    {
	      oldMrna = mrna ; oldX2 = 0 ;
	    }
	  if (1)
	    {
	      int da = a2 > a1 ? a2 - a1 : a1 - a2 ;
	      x1 = oldX2 + 1 ; x2 = x1 + da ;
	      oldX2 = x2 ;
	    }	  
	}
      if (x1 > x2) { int dummy = x1 ; x1 = x2 ; x2 = dummy ; dummy = a1 ; a1 = a2 ; a2 = dummy ; }
      if (target)
	{
	  up = arrayp (sx->shadowArray, nn++, SHADOW) ;
	  up->mrna = mrna  ; up->gene = gene ; up->target = target ;
	  up->a1 = a1 ; up->x1 = x1 ;  
	  up->a2 = a2 ; up->x2 = x2 ;
	}
      gene = target = 0 ;
    }

  arraySort (sx->shadowArray, shadowOrder) ;
  arrayCompress (sx->shadowArray) ;

  int ii, jj ;
  nn = arrayMax (sx->shadowArray) ;
  for (ii = jj = 0, up = vp = arrayp (sx->shadowArray, 0, SHADOW) ; ii < nn ; ii++, up++)
    {
      /* faire un clean up si redondant */
      /* crash si incoherent */
      SHADOW *wp = ii < nn - 1 ? up + 1 : 0 ;
      
      if (! up->mrna) continue ;
      if (wp && up->mrna == wp->mrna && up->x2 >= wp->x1) /* potential overlap */
	{
	  if (
	      (up->target != vp->target) ||
	      (up->a1 < up->a2 && wp->a1 > wp->a2) ||
	      (up->a1 > up->a2 && wp->a1 < wp->a2) ||
	      (up->a1 < up->a2 && up->a1 + wp->x1 - up->x1 != wp->a1) ||
	      (up->a1 > up->a2 && up->a1 - wp->x1 + up->x1 != wp->a1)  ||
	      (up->a1 < up->a2 && up->a2 + wp->x2 - up->x2 != wp->a2) ||
	      (up->a1 > up->a2 && up->a2 - wp->x2 + up->x2 != wp->a2)
	      )
	    messcrash ("Parse error, distinct overlap in %s::%d:%d %s::%d:%d %s::%d:%d %s::%d:%d %s line %d"
		       , dictName (sx->shadowDict, up->mrna)
		       , up->x1, up->x2
		       , dictName (sx->selectDict, up->target)
		       , up->a1, up->a2
		       , dictName (sx->shadowDict, wp->mrna)
		       , wp->x1, wp->x2
		       , dictName (sx->selectDict, wp->target)
		       , wp->a1, wp->a2
		       , aceInFileName (ai) ? aceInFileName (ai) : "-"
		       , aceInStreamLine (ai)
		       ) ;
	  if (up->a1 <= up->a2) { up->x2 = wp->x1 - 1 ; up->a2 = wp->a1 - 1 ; }
	  else                  { up->x2 = wp->x1 - 1 ; up->a2 = wp->a1 + 1 ; }
	}
      if (vp <= up && up->mrna && up->x1 <= up->x2)
	{ *vp = *up ; vp++ ; jj++ ;
	  if (0) fprintf (stderr, "%s %d %d %s %d %d\n"
			  , dictName (sx->shadowDict, up->mrna)
			  , up->x1, up->x2
			  , dictName (sx->selectDict, up->target)
			  , up->a1, up->a2
			  ) ;
	} 
    }
  arrayMax (sx->shadowArray) = jj ;
  arraySort (sx->shadowArray, shadowOrder) ;
 done:
  ac_free (ai) ;
  ac_free (h) ;
  return nn ;
} /* parseShadowFile */

/*************************************************************************************/

static DICT *gtf2ace (BOOL isGff3, DICT *dict, KEYSET ks, char *buf, AC_HANDLE h)
{
  DICT *itemDict = dictHandleCreate (100, h) ;
  char *cp, *cq, *cr ;
  int item , n ;
  
  dictAdd (dict, "ZERO", 0) ;
  cq = buf - 1 ; 
  while (cq && (cp = cq + 1) && *cp)
    {
      cq = strstr (cp, ";") ;
      if (cq)
	*cq = 0 ;
      while (*cp == ' ') cp++ ;
      if (isGff3)
	cr = strstr (cp, "=") ;             /* GFF3 as at NCBI */
      else
	cr = strstr (cp, " ") ;  /* GTF as at EBI */
      if (! cr) continue ;
      *cr = 0 ;
      dictAdd (itemDict, cp, &item) ;
      keySet (ks, item) = 0 ;
      cp = cr + 1 ;
      if (*cp)
	{
	  cr = ac_unprotect (cp, h) ;
	  if (0) url_decode_inplace(cr) ;
	  if (1 && ! strcmp (cr, "protein_coding"))
	    cr = "CDS" ;
	  dictAdd (dict, cr, &n) ;
	  keySet (ks, item) = n ;
	}
    }
  if (isGff3 && dictFind (itemDict, "Dbxref", &item))
    {
      n =  keySet (ks, item) ;
      buf = strnew (dictName (dict, n), 0) ;
      cq = buf - 1 ; 
      while (cq && (cp = cq + 1) && *cp)
	{
	  cq = strstr (cp, ",") ;
	  if (cq)
	    *cq = 0 ;
	  while (*cp == ' ') cp++ ;
	  cr = strstr (cp, ":") ;             /* GFF3 as at NCBI */
	  if (! cr) continue ;
	  *cr = 0 ;
	  dictAdd (itemDict, cp, &item) ;
	  keySet (ks, item) = 0 ;
	  cp = cr + 1 ;
	  if (*cp)
	    {
	      cr = ac_unprotect (cp, h) ;
	      if (0) url_decode_inplace(cr) ;	
	      dictAdd (dict, cr, &n) ;
	      keySet (ks, item) = n ;
	    }
	}
      ac_free (buf) ;
    }
  return itemDict ;
} /* gtf2ace */

/*************************************************************************************/

static int gtfItem (const char *tag, DICT *itemDict, KEYSET ks) 
{
  int item = 0 ;
  if (dictFind (itemDict, tag, &item))
    return keySet (ks, item) ;

  return 0 ;
} /* gtfItem */

/* not used yet
static int gtfItemInt (const char *tag, DICT *itemDict, DICT *dict, KEYSET ks) 
{
  int n = 0, item = 0 ;

  if (dictFind (itemDict, tag, &item))
    n = atoi (dictName(dict, keySet (ks, item))) ;

  return n ;
}  gtfItemInt
*/
/*************************************************************************************/

static int shadowGeneName (SX *sx, DICT *dict, SHADOW *shadow, AC_HANDLE h)
{
  int g, gn ;
  KEYSET ks ;

  if (! shadow->gene) 
    shadow->gene = shadow->ID ;

  g =  shadow->gene ;
  gn =  shadow->gene_name ;
  if (! gn)
    gn = shadow->name ;
  if (! gn)
    return g ;
  if (! g)
    return 0 ;

  ks = array (sx->gnam2genes, gn, KEYSET) ;
  if (ks)
    {     
      int j ;
      char *cp ;
      
      for (j = 0 ; j < keySetMax (ks) ; j++)
	if (keySet(ks, j) == g)
	  break ;
      cp = hprintf (h, "%s.%d", dictName(dict, gn), j+1) ;
      dictAdd (dict, cp, &gn) ;
    }
  return gn ;
} /* shadowGeneName */

/*************************************************************************************/
/* gather the geometry off the shadow file, tested on Encode 37.70 gtf file */

static void parseGtfFeature (SX *sx, const char *featureType, const char *fileSuffix, int showAll, BOOL wantCDS, KEYSET mrna2gene)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (sx->gtfFileName, FALSE, h) ;
  KEYSET ks = keySetHandleCreate (h) ;
  int i, j, target, a1, a2, nShadow = 0, type, subtype ;
  char *cp, cutter ;
  int old = 0, oldParent = 0, x1 = 0, x2 = 0 ;
  BOOL isDown = TRUE ;
  BOOL hasCDS = FALSE ;
  BOOL hasEXON = FALSE ;
  BOOL wantGene = !strcasecmp (featureType, "gene") ;
  BOOL wantMRNA = !strcasecmp (featureType, "mRNA") ;
  BOOL wantrRNA = !strcasecmp (featureType, "rRNA") ;
  BOOL wanttRNA = !strcasecmp (featureType, "tRNA") ;
  BOOL wantncRNA = !strcasecmp (featureType, "ncRNA") ;
  BOOL wantmiRNA = !strcasecmp (featureType, "miRNA") ;
  BOOL wantpre_miRNA = !strcasecmp (featureType, "pre_miRNA") ;
  BOOL wantExon = !strcasecmp (featureType, "exon") ;
  int mRNA_type, tRNA_type, rRNA_type, ncRNA_type, miRNA_type, pre_miRNA_type, Exon_type, CDS_type, gene_type, pseudogene_type ;
  DICT *itemDict, *dict ;
  SHADOW *shadow, *shadow2 ;
  Array shadows ; 
  ACEOUT ao = 0 ;

  dict = sx->shadowDict ;
  if (! sx->shadowDict)
    {
      sx->selectDict = dictHandleCreate (10000, sx->h) ;
      dict = sx->shadowDict = dictHandleCreate (10000, sx->h) ;
    }

  if (wantGene) 
    ao = aceOutCreate (sx->outFileName, fileSuffix, sx->gzo, h) ;

  if (wantGene )
    {
      /* the gene shadows must not be part of sx->shadowArray
       * which is used to export the mRNAs fasta file
       */
      shadows = arrayHandleCreate (10000, SHADOW, sx->h) ;
    }
  else
    {
      if (! sx->shadowArray)
	sx->shadowArray = arrayHandleCreate (10000, SHADOW, sx->h) ;
      shadows = sx->shadowArray ;
    }

  if (0) arrayMax (shadows) = 0 ;
  nShadow = arrayMax (shadows) ;

  dictAdd (dict, "ZERO", 0) ;
  dictAdd (dict, "mRNA", &mRNA_type) ;
  dictAdd (dict, "tRNA", &tRNA_type) ;
  dictAdd (dict, "rRNA", &rRNA_type) ;
  dictAdd (dict, "ncRNA", &ncRNA_type) ;
  dictAdd (dict, "miRNA", &miRNA_type) ;
  dictAdd (dict, "pre_miRNA", &pre_miRNA_type) ;
  dictAdd (dict, "exon", &Exon_type) ;
  dictAdd (dict, "CDS", &CDS_type) ;
  dictAdd (dict, "gene", &gene_type) ;
  dictAdd (dict, "pseudogene", &pseudogene_type) ;

  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      /* chromosome */
      cp = aceInWordCut (ai, "\t", &cutter) ;
      if (! cp || *cp == '#')	continue ;
      dictAdd (sx->selectDict, cp, &target) ;

      /* subtype nonsense_mediated_decay */
      aceInStep (ai, '\t') ;
      cp = aceInWordCut (ai, "\t", &cutter) ;
      if (! cp)	continue ;
      dictAdd (dict, cp, &subtype) ;

      /* type gene exon intron CDS start_codon stop_codon ... */
      aceInStep (ai, '\t') ;
      cp = aceInWordCut (ai, "\t", &cutter) ;
      if (! cp)	continue ;
      if (!strcasecmp (cp, "transcript_region"))
	cp = "mrna" ;
      if (!strcasecmp (cp, "antisense_lncRNA"))
	cp = "ncrna" ;
      if (!strcasecmp (cp, "antisense_RNA"))
	cp = "mrna" ;
      if (!strcasecmp (cp, "lnc_RNA"))
	cp = "ncrna" ;
      if (!strcasecmp (cp, "ncRNA"))
	cp = "ncrna" ;
      if (!strcasecmp (cp, "tRNA"))
	cp = "trna" ;
      if (!strcasecmp (cp, "rRNA"))
	cp = "rrna" ;
      if (!strcasecmp (cp, "miRNA"))
	cp = "mirna" ;
      if (!strcasecmp (cp, "pre_miRNA"))
	cp = "pre_mirna" ;
      if (!strcasecmp (cp, "pseudogenic_transcript"))
	cp = "mrna" ;
      if (!strcasecmp (cp, "exon"))
	cp = "exon" ;
      if (!strcasecmp (cp, "pseudogenic_exon"))
	cp = "exon" ;
      if (!strcasecmp (cp, "pseudogene"))
	cp = "gene" ;
      if ( wantGene)
	{
	  if (strcasecmp (cp, "gene") &&
	      (! sx->gffBacteria || strcasecmp (cp, "mrna"))
	      )
	    continue ;
	}
      else if (wantMRNA)
	{
	  if (strcasecmp (cp, "mRNA"))
	    continue ;
	}
      else if (wantrRNA)
	{
	  if (strcasecmp (cp, "rRNA"))
	    continue ;
	}
      else if (wanttRNA)
	{
	  if (strcasecmp (cp, "tRNA"))
	    continue ;
	}
      else if (wantncRNA)
	{
	  if (strcasecmp (cp, "ncRNA"))
	    continue ;
	}
      else if (wantmiRNA)
	{
	  if (strcasecmp (cp, "miRNA"))
	    continue ;
	}
      else if (wantpre_miRNA)
	{
	  if (strcasecmp (cp, "pre_miRNA"))
	    continue ;
	}
      else if (wantCDS)
	{
	  if (!strcasecmp (cp, "stop_codon"))
	    cp = "cds" ;
	  if (strcasecmp (cp, "cds") &&
	      (! sx->gffBacteria || strcasecmp (cp, "mrna"))
	      )
	    continue ;
	}
      else
	{
	  if (! strcasecmp (cp, "cds")) 
	    continue ;
	  
	  if (sx->isGff3)
	    {
	      char *cq = aceInPos (ai) ;
	      if (strcasecmp (cp, "rRNA") && 
		  strcasecmp (cp, "tRNA") && 
		  strcasecmp (cp, "ncRNA") && /*  in mouse the ncRNA seem to have exons */
		  strcasecmp (cp, "exon") && 
		  (strcasecmp (cp, "mRNA") || ! sx->gffBacteria) && 
		  (! cq || ! strstr (cq, "Parent=rna"))
		  )
		continue ;
	    }
	  else
	    {
	      if (strcasecmp (cp, featureType))
		continue ;
	    }
	}
	
      dictAdd (dict, cp, &type) ;

      /* coordinates on chromosome */
      a1 = a2 = 0 ;
      aceInStep (ai, '\t') ;
      aceInInt (ai, &a1) ;
      aceInStep (ai, '\t') ;
      aceInInt (ai, &a2) ;

      /* silly . discard */
      aceInStep (ai, '\t') ;
      cp = aceInWordCut (ai, "\t", &cutter) ;
      if (! cp)	continue ;

      /* strand */
      aceInStep (ai, '\t') ;
      cp = aceInWordCut (ai, "\t", &cutter) ; 
      if (! cp)	continue ;
      if (*cp == '-')
	{ int a0 = a1 ; a1 = a2 ; a2 = a0 ; isDown = FALSE ; }
      else
	isDown = TRUE ;

      /* silly 0 discard */
      aceInStep (ai, '\t') ;
      cp = aceInWordCut (ai, "\t", &cutter) ;
      if (! cp)	continue ;

      /* acedb like stuff, semi-column delimited in column 9 */
      aceInStep (ai, '\t') ;
      cp = aceInWordCut (ai, "\t", &cutter) ;
      if (! cp)	continue ;
      itemDict = gtf2ace (sx->isGff3, dict, ks, cp, h) ;
      
      /* interpret the data */
      {
	if (showAll == 2)
	  {
	    int p, m ;
	    p = gtfItem ("Parent", itemDict, ks) ;
	    m = gtfItem ("transcript_id", itemDict, ks) ;
	    keySet (sx->rnaId2mrna, p) = m ; /* go from exon parent to mrna name */
	    continue ;
	  }
      
	shadow = arrayp (shadows, nShadow++, SHADOW) ;
	shadow->target = target ;
	shadow->a1 = a1 ; shadow->a2 = a2 ;
	if (sx->isGff3)
	  {
	    shadow->ID = gtfItem ("ID", itemDict, ks) ;
	    shadow->Parent = gtfItem ("Parent", itemDict, ks) ;
	    shadow->GeneID = gtfItem ("GeneID", itemDict, ks) ;
	    shadow->ncRNA = gtfItem ("ncrna_class", itemDict, ks) ;
	    shadow->mrna = gtfItem ("transcript_id", itemDict, ks) ;
	    shadow->type = type ;
	    shadow->locus_tag = gtfItem ("locus_tag", itemDict, ks) ;
	    shadow->gene_name = gtfItem ("gene", itemDict, ks) ;
	    shadow->name = gtfItem ("name", itemDict, ks) ;
	    shadow->gene_title = gtfItem ("product", itemDict, ks) ;
	    shadow->note = gtfItem ("Note", itemDict, ks) ;
	    shadow->Dbxref = gtfItem ("Dbxref", itemDict, ks) ;
	    shadow->protein_id = gtfItem ("protein_id", itemDict, ks) ;
	    shadow->gene_biotype = gtfItem ("gene_biotype", itemDict, ks) ;
	    shadow->isDown = isDown ;

	    if (sx->gffTAIR)
	      {
		if (wantMRNA)
		  {
		    shadow->mrna = shadow->ID ;
		    shadow->gene = shadow->Parent ;
		    keySet (mrna2gene, shadow->mrna) = shadow->gene ;
		  }
		if (wantExon || wantCDS)
		  {
		    shadow->mrna = shadow->Parent ;
		    shadow->gene = keySet (mrna2gene,shadow->mrna) ;
		  }
	      }

	    if (wantCDS) shadow->cds = TRUE ;
	    if (! shadow->gene_biotype)
	      shadow->gene_biotype = shadow->type ;

	    if (sx->gffBacteria)
	      {	 
		if (shadow->GeneID)
		  {
		    if (wantGene)
		      shadow->gene = shadow->GeneID ;
		    else
		      {
			shadow->gene = shadow->GeneID ;
			shadow->mrna = shadow->GeneID ;
		      }
		  }
	      }
	    else
	      {
		if (! sx->gffTAIR && wantCDS &&  shadow->Parent) shadow->ID = shadow->Parent ;
		if (! wantCDS)
		  shadow->gene = gtfItem ("gene", itemDict, ks) ;
		if (! shadow->mrna &&
		    ( type == rRNA_type || type == tRNA_type || type == ncRNA_type)
		    )
		  shadow->mrna = shadow->ID ;
	      }
	    if (shadow->mrna && ! shadow->gene) 
	      shadow->gene = keySet (mrna2gene, shadow->mrna) ;
	    if (shadow->gene_name && ! shadow->gene) 
	      shadow->gene = shadow->gene_name ;
	    if (shadow->mrna && ! shadow->gene) 
	      shadow->gene = shadow->mrna ;
	    if (wantMRNA)
	      {
		shadow->gene = shadow->Parent ;
		keySet (mrna2gene, shadow->mrna) = shadow->gene ;
	      }
	    if (!  shadow->GeneID && shadow->Dbxref)
	      {
		const char *cp = strstr (dictName(dict, shadow->Dbxref), "GeneID:") ;
		if (cp)
		  {
		    char *cq = strnew (cp+7,0) ;
		    char *cr = strstr (cq, ",") ;
		    if (cr) *cr = 0 ;
		    dictAdd (dict, cq, &(shadow->GeneID)) ;
		    messfree (cq) ;
		  }
	      }
	  }
	else
	  {
	    shadow->gene = gtfItem ("gene_id", itemDict, ks) ;
	    shadow->GeneID = gtfItem ("NCBI_GeneID", itemDict, ks) ;
	    shadow->mrna = gtfItem ("transcript_id", itemDict, ks) ;
	    shadow->protein_id = gtfItem ("protein_id", itemDict, ks) ;
	    shadow->gene_title = gtfItem ("gene_name", itemDict, ks) ;
	    shadow->description = gtfItem ("description", itemDict, ks) ;	    
	    shadow->gene_biotype = gtfItem ("gene_biotype", itemDict, ks) ;
	    if (shadow->mrna && ! shadow->gene) shadow->gene = shadow->mrna ;
	    shadow->isDown = isDown ;
	    if (wantCDS) shadow->cds = TRUE ;
	    shadow->exon = TRUE ;
	    shadow->Dbxref = gtfItem ("db_xref", itemDict, ks) ;
	    if (!  shadow->GeneID && shadow->Dbxref)
	      {
		const char *cp = strstr (dictName(dict, shadow->Dbxref), "GeneID:") ;
		if (cp)
		  {
		    char *cq = strnew (cp+7,0) ;
		    char *cr = strstr (cq, ",") ;
		    if (cr) *cr = 0 ;
		    dictAdd (dict, cq, &(shadow->GeneID)) ;
		    messfree (cq) ;
		  }
	      }
	  }
      }
    }
  ac_free (ai) ;
  if (showAll == 2)
    goto done ;

  if (wantGene)
    {
      KEYSET gnam2g = keySetHandleCreate (h) ;
      for (i = 0 ; i < arrayMax (shadows) ; i++)
	{
	  shadow = arrp(shadows, i, SHADOW) ;
	  if (! shadow->gene_name)
	    continue ;
	  if (shadow->gene)
	    {
	      int j, g = keySet (gnam2g, shadow->gene_name) ;
	      if (!g)
		keySet (gnam2g, shadow->gene_name) = shadow->gene ;
	      if (keySet (gnam2g, shadow->gene_name) != shadow->gene)
		{
		  KEYSET ks = array (sx->gnam2genes, shadow->gene_name, KEYSET) ;
		  
		  if (! ks)
		    {
		      ks = array (sx->gnam2genes, shadow->gene_name, KEYSET)
			= keySetHandleCreate (sx->h) ;
		      keySet (ks, 0) = g ;
		    }
		  g = shadow->gene ;
		  for (j = 0 ; j < keySetMax (ks) ; j++)
		    if (keySet(ks, j) == g)
		      break ;
		  keySet(ks, j) = g ;
		}
	    }
	}
    }

  if (! sx->gffTAIR && (wantGene || wantCDS)) /* clean doubles */
    {
      KEYSET badgn = keySetCreate () ;
      KEYSET gn2g = keySetCreate () ;
      int g, gn, gnn ;
      
      for (i = 0 ; i < arrayMax (shadows) ; i++)
	{	  
	  shadow = arrp(shadows, i, SHADOW) ;
	  g = shadow->GeneID ;
	  gn = shadow->gene_name ;
	  gnn = keySet (gn2g, gn) ;
	  if (! gnn)
	    keySet (gn2g, gn) = gnn = g ;
	  if (g != gnn)
	    keySet (badgn, gn) = 1  ;
	}
      
      /* rename the top of the list */
      for (i = 0 ; ! sx->gffTAIR && i < arrayMax (shadows) ; i++)
	{
	  shadow = arrp(shadows, i, SHADOW) ;
	  g = shadow->GeneID ;
	  gn = shadow->gene_name ;
	  if (gn)
	    {
	      gnn = keySet (gn2g, gn) ;
	      if (keySet (badgn, gn) == 1)
		{
		  char *cp = hprintf (h, "%s.%s"
				      , dictName (dict, gn)
				      , dictName (dict, g)
				  ) ;
		  dictAdd (dict, cp, &gn) ;
		}
	      if (shadow->mrna == shadow->gene)
		shadow->mrna = gn ;
	      shadow->gene = gn ;
	    }
	}
      keySetDestroy (badgn) ;
      keySetDestroy (gn2g) ;
    }

  if (wantGene)
    {
      int oldGene = 0, oldID = 0 ;
      
      for (i = 0 ; i < arrayMax (shadows) ; i++)
	{
	  int g ;

	  shadow = arrp(shadows, i, SHADOW) ;
	  g =  shadowGeneName (sx, dict, shadow, h) ;
	  if (! g)
	    continue ;
	  if (oldGene == shadow->gene && oldID == shadow->GeneID)
	    continue ;
	  oldGene = shadow->gene ;
	  oldID = shadow->GeneID ;

	  aceOutf (ao, "Gene \"%s\"\n",dictName (dict, shadow->gene)) ;
	  if (shadow->locus_tag)
	    aceOutf (ao, "Locus \"%s\"\n",dictName (dict, shadow->locus_tag)) ;
	  if (shadow->gene_title)
	    {
	      aceOutf (ao, "Title \"%s\"\n",dictName (dict, shadow->gene_title)) ;
	      if (shadow->note)
		aceOutf (ao, "Concise_description \"%s\"\n",dictName (dict, shadow->note)) ;
	      if (shadow->description)
		aceOutf (ao, "Concise_description \"%s\"\n",dictName (dict, shadow->description)) ;
	    }
	  else
	    {
	      if (shadow->description)
		aceOutf (ao, "Title \"%s\"\n",dictName (dict, shadow->description)) ;
	    }
	  if (shadow->gene_name)
	    aceOutf (ao, "LocusLink \"%s\"\n",dictName (dict, shadow->gene_name)) ;
	  if (0 && shadow->locus_tag)
	    aceOutf (ao, "WbId \"%s\"\n",dictName (dict, shadow->locus_tag)) ;
	  if (shadow->Dbxref)
	    {
	      char *cq, *cp = strnew (dictName (dict, shadow->Dbxref), h) ;
	      while (cp)
		{
		  cq = strchr (cp, ',') ;
		  if (cq) *cq = 0 ;
		  if (! strncasecmp (cp, "PMID:", 5))
		    aceOutf (ao, "Reference pm%s\n", cp+5) ;
		  cp = cq ? cq + 1 : 0 ;
		}
	    }
	  if (shadow->GeneID)
	    aceOutf (ao, "GeneID \"%s\"\n",dictName (dict, shadow->GeneID)) ;
	  if (a1 && a2)
	    aceOutf (ao, "IntMap \"%s\" %d %d\n",  dictName (sx->selectDict, shadow->target), shadow->a1, shadow->a2) ;
	  aceOut (ao, "\n") ;
	}
      goto done ;
    }

  /* duplicate in case of multiple inheritance (TAIR) */
  if (sx->gffTAIR) 
    {
      for (i = 0 ; i < arrayMax (shadows) ; i++)
	{
	  const char *ccp ;
	  shadow = arrp(shadows, i, SHADOW) ;
	  if (! shadow->Parent)
	    continue ;
	  ccp = dictName (dict, shadow->Parent) ;
	  if (strchr (ccp, ','))
	    {
	      char *cr, *cq, *cq0 = strnew (ccp, h) ;
	      cq = cq0 ;
	      while (cq)
		{
		  int p1 ;
		  cr = strchr (cq, ',') ;
		  if (cr) *cr = 0 ;
		  dictAdd (dict, cq, &p1) ;
		  if (cq == cq0)
		    {
		      shadow->Parent = p1 ;		
		      if (wantMRNA)
			{
			  shadow->mrna = shadow->ID ;
			  shadow->gene = shadow->Parent ;
			  keySet (mrna2gene, shadow->mrna) = shadow->gene ;
			}
		      if (wantExon || wantCDS)
			{
			  shadow->mrna = shadow->Parent ;
			  shadow->gene = keySet (mrna2gene,shadow->mrna) ;
			}
		    }

		  else
		    {   /* duplicate the record */
		      shadow2 = arrayp (shadows, arrayMax (shadows), SHADOW) ;
		      shadow = arrp(shadows, i, SHADOW) ; /* we may have relocalized shadows */
		      *shadow2 = *shadow ;
		      shadow2->Parent = p1 ;
		      if (wantMRNA)
			{
			  shadow2->mrna = shadow2->ID ;
			  shadow2->gene = shadow2->Parent ;
			  keySet (mrna2gene, shadow2->mrna) = shadow2->gene ;
			}
		      if (wantExon || wantCDS)
			{
			  shadow2->mrna = shadow2->Parent ;
			  shadow2->gene = keySet (mrna2gene,shadow2->mrna) ;
			}
		    }
		  cq = cr ? cr + 1 : 0 ;
		}
	    }
	}
    }

  /* add the mrna name in the CDS */
  if (! wantCDS && sx->isGff3)
    {
      arraySort (shadows, shadowParentOrder) ;
      if (0) showShadow (sx, dict, shadows) ;
      for (i = 0 ; i < arrayMax (shadows) ; i++)
	{
	  shadow = arrp(shadows, i, SHADOW) ;
	  if (shadow->cds && ! shadow->mrna && shadow->Parent)
	    {
	      for (j = i + 1, shadow2 = shadow + 1 ; j < arrayMax (shadows) && shadow2->Parent == shadow->Parent ; j++, shadow2++)
		{
		  if (shadow2->mrna)
		    { 
		      shadow->gene = shadow2->gene ;
		      shadow->mrna = shadow2->mrna ;
		      break ;
		    }
		}
	    }
	}

      if (arrayMax (shadows)) /* remove the exons and cds that do not belong to an mRNA */
	{
	  BitSet bb = bitSetCreate (dictMax (dict) + 1, h) ;
	  for (i = j = 0, shadow = shadow2 = arrp(shadows, i, SHADOW); i < arrayMax (shadows) ; i++, shadow++)
	    {
	      if (! shadow->mrna)
		shadow->mrna = shadow->Parent ;
	      if (shadow->mrna && shadow->type == Exon_type)
		bitSet (bb, shadow->mrna) ;
	      if (i && shadow->type != Exon_type &&	
		  (shadow-1)->type == Exon_type &&
		  (shadow-1)->mrna == shadow->mrna
		  )
		continue ;
	      if (j < i) { *shadow2 = *shadow ; }
	      j++ ; shadow2++ ;
	    }
	  arrayMax (shadows) = j ;
	  
	  /* remove the type for mrna who have an exon */
	  
	  if (arrayMax (shadows))
	    {
	      for (i = j = 0, shadow = shadow2 = arrp(shadows, i, SHADOW); i < arrayMax (shadows) ; i++, shadow++)
		{
		  if (shadow->mrna && shadow->type != CDS_type && shadow->type != Exon_type && bit (bb, shadow->mrna))
		    continue ;
		  if (j < i) { *shadow2 = *shadow ; }
		  j++ ; shadow2++ ;
		}
	      arrayMax (shadows) = j ;
	    }
	}
    }

  if (0) showShadow (sx, dict, shadows) ;
  /* clean up the duplicate gene names: same gene but different ID */
  arraySort (shadows, shadowOrder) ;
  if (! sx->isGff3)
    for (i = 0 ; i < arrayMax (shadows) ; i++)
      {
	int n, j, gene, oldID ;
	const char *ccp ;
	
	shadow = arrp(shadows, i, SHADOW) ;
	ccp = dictName (dict, shadow->gene)  ;
	
	for (shadow2 = shadow, j = i, n = 0, gene = shadow->gene, oldID = shadow->ID ; 
	     j <  arrayMax (shadows) && shadow2->gene == gene ; j++, shadow2++) 
	  if (shadow2->ID != oldID) 
	    { n++ ; oldID = shadow2->ID ; }
	if (n > 1) /* edit the gene names */
	  for (shadow2 = shadow, j = i, n = 0, gene = shadow->gene, oldID = 0 ;
	       j <  arrayMax (shadows) && shadow2->gene == gene ; j++, shadow2++) 
	    {
	      if (shadow2->ID != oldID) 
		{ n++ ; oldID = shadow2->ID ; }
	      dictAdd (dict, messprintf ("%s.%d", ccp, n), &(shadow2->gene)) ;
	      if (shadow2->mrna == gene)
		shadow2->mrna = shadow2->gene ;
	    }
	i = j - 1 ;
      }
  
  /* regularize the coordinates in order to number the exons */
  if (arrayMax (shadows))
    {
      for (i = 0, shadow = arrp(shadows, i, SHADOW); i < arrayMax (shadows) ; i++, shadow++)
	{
	  if (!wantCDS && shadow->cds) continue ;
	  if (! shadow->isDown)
	    { shadow->a1 = -shadow->a1 ; shadow->a2 = -shadow->a2 ; }
	  shadow->x1 = shadow->a1 ; /* serves as ordinal number for the exons */
	  shadow->x2 = 0 ;
	}
      arraySort (shadows, shadowOrder) ;
      arrayCompress (shadows) ;

      /* restore */
      for (i = 0, shadow = arrp(shadows, i, SHADOW); i < arrayMax (shadows) ; i++, shadow++)
	{
	  if (!wantCDS && shadow->cds) continue ;
	  if (! shadow->isDown)
	    { shadow->a1 = -shadow->a1 ; shadow->a2 = -shadow->a2 ; }
	}
      /* replace the ordinal exon number by the spliced mrna coordinates */
      for (old = oldParent = -1, i = x2 = 0, shadow = arrp(shadows, i, SHADOW); i < arrayMax (shadows) ; i++, shadow++)
	{
	  if (!wantCDS && shadow->cds) continue ;
	  if (shadow->mrna != old ||
	      (sx->isGff3 && shadow->Parent != oldParent)
	      )
	    { x1 = x2 = 0 ; j = 0 ; }
	  
	  if (wantCDS || ! shadow->cds) 
	    {
	      x1 = shadow->x1 = x2 + 1 ;
	      shadow->x2 = shadow->isDown ? shadow->x1 + (shadow->a2 - shadow->a1) : shadow->x1 - (shadow->a2 - shadow->a1) ;
	      x2 = shadow->x2 ; 
	      a1 = shadow->a1 ; a2 = shadow->a2 ;
	    }
	  else
	    {
	      if (shadow->isDown && shadow->a1 <= a2)
		{
		  shadow->x1 = x1 + shadow->a1 - a1 ;
		  shadow->x2 = x1 + shadow->a2 - a1 ;
		}
	      else if (shadow->isDown && shadow->a1 > a2)
		{
		  shadow->x1 = x2 + 1 ;
		  shadow->x2 = shadow->x1 + (shadow->a2 - shadow->a1) ;
		}
	      else if (! shadow->isDown && shadow->a1 >= a2)
		{
		  shadow->x1 = x1 + a1 - shadow->a1 ;
		  shadow->x2 = x1 + a1 - shadow->a2 ;
		}
	      else if (! shadow->isDown && shadow->a1 < a2)
		{
		  shadow->x1 = x2 + 1 ;
		  shadow->x2 = shadow->x1 + (shadow->a1 - shadow->a2) ;
		}
	    }
	  old = shadow->mrna ;
	  oldParent = shadow->Parent ;
	}
    }

  if (! sx->gtfRemapPrefix)
    goto done ;

  ac_free (ao) ;

  if (arrayMax (shadows) &&
      /* sx->gffBacteria && */
      wantCDS
      )
    {
      ao = aceOutCreate (sx->outFileName, ".cds2genename.ace", sx->gzo, h) ;
      for (old = i = 0, shadow = arrp(shadows, i, SHADOW); i < arrayMax (shadows) ; i++, shadow++)
	{
	  if (shadow->gene_name &&
	      shadow->gene_title
	      )
	    aceOutf (ao, "Gene \"%s\"\nTitle \"%s\"\n\n"
		     , dictName (dict, shadow->gene)
		     , dictName (dict, shadow->gene_title)
		   ) ;
	}
      ac_free (ao) ;
      
      ao = aceOutCreate (sx->outFileName, ".goodProduct.ace", sx->gzo, h) ;
      for (old = i = 0, shadow = arrp(shadows, i, SHADOW); i < arrayMax (shadows) ; i++, shadow++)
	{
	  if (shadow->name)
	    {
	      int dx = shadow->a2 - shadow->a1 ;
	      int g1 = 0, g = shadow->gene ;
	      int p = shadow->Parent ;

	      if (! sx->gffBacteria && sx->rnaId2mrna && p && p < keySetMax (sx->rnaId2mrna))
		g1 = keySet (sx->rnaId2mrna, p) ;
	      if (g1) g = g1 ;
	      if (dx < 0) dx = -dx ;
	      dx++ ;
	      if (sx->gffBacteria && !g && p)
		{
		  g = p ;
		  aceOutf (ao, "Sequence \"%s\"\nmRNA \"%s\" %d %d\n\n"
			   , dictName (sx->selectDict, shadow->target)
			   , dictName (dict, g)
			   , shadow->a1, shadow->a2
			   ) ;
		  aceOutf (ao, "mRNA \"%s\"\nSource_exons %d %d\n\n"
			   , dictName (dict, g)
			   , 1, dx
			   ) ;
		}
	      if (g)
		aceOutf (ao, "mRNA \"%s\"\nProduct \"%s\" %d %d\n\n"
			 , dictName (dict, g)
			 , dictName (dict, shadow->name)
			 , 1, dx
			 ) ;
	      if (p)
		aceOutf (ao, "Product \"%s\"\nGood_product\nBest_product\n\n"
			 , dictName (dict, p)
			 ) ;
	      if (p && shadow->gene_title)
		aceOutf (ao, "Product \"%s\"\nTitle \"%s\"\n\n"
			 , dictName (dict, p)
			 , dictName (dict, shadow->gene_title)
		       ) ;
	    }
	}
      ac_free (ao) ;
    }


  ao = aceOutCreate (sx->outFileName, fileSuffix, sx->gzo, h) ;
  if (arrayMax (shadows))
    {
      for (old = i = 0, shadow = arrp(shadows, i, SHADOW); i < arrayMax (shadows) ; i++, shadow++)
	{
	  if (shadow->mrna && (sx->gffBacteria || shadow->cds == wantCDS))
	    aceOutf (ao, "%s\t%s%s\t%d\t%d\t%s\t%d\t%d\t%s%s\n"
		     , sx->gtfRemapPrefix ? sx->gtfRemapPrefix : "-"
		     , dictName(dict, shadow->mrna), showAll ? "" : "_CDS"
		     , shadow->x1, shadow->x2
		     , dictName(sx->selectDict, shadow->target)
		     , shadow->a1, shadow->a2
		     , shadow->gene ? dictName(dict, shadow->gene) : ""
		     , showAll ? "" : "_CDS"
		     ) ;
	}
    }

  ac_free (ao) ;
  if (wantGene || wantMRNA || wantCDS)
    goto done ;
  ao = aceOutCreate (sx->outFileName, ".introns", sx->gzo, h) ;
  aceOutDate (ao, "###", sx->title) ;
  if (! sx->gffTAIR && arrayMax (shadows))
    {
      for (old = i = 0, shadow = arrp(shadows, i, SHADOW); i < arrayMax (shadows) ; i++, shadow++)
	{
	  if (shadow->cds) continue ;
	  isDown = shadow->isDown ;
	  if (old && shadow->mrna == old && (shadow->a1 - a2 > 12 || shadow->a1 - a2 < -12))
	    aceOutf (ao, "%s\t%s\t%s\t%d\t%d\t%s\t%09d\t%09d\t%s\n"
		     , sx->gtfRemapPrefix ? sx->gtfRemapPrefix : "-"
		     , dictName(dict, shadow->gene)
		     , dictName(dict, shadow->mrna)
		     , x2 , shadow->x1
		     , dictName(sx->selectDict, shadow->target)
		     , a2 + (isDown ? 1 : -1) , shadow->a1 - (isDown ? 1 : -1)
		     , isDown ? "+" : "-"
		     ) ;
	  old = shadow->mrna ;
	  x2 = shadow->x2 ; a2 = shadow->a2 ;
	}	 
    }
  ac_free (ao) ;

  if (! strcasecmp (fileSuffix, ".mrnaRemap"))
    {
      int typ ;
      BOOL isDown = FALSE, cds = FALSE ;
      char *suffix = "" ;

      for (typ = 0 ; typ < 4 ; typ++)
	{
	  switch (typ)
	    {
	    case 0:
	      suffix = ".f.sponge" ;
	      isDown = TRUE ; cds = FALSE ;
	      break ;
	    case 1:
	      suffix = ".r.sponge" ;
	      isDown = FALSE ; cds = FALSE ;
	      break ;
	    case 2:
	      suffix = ".f.cds_sponge" ;
	      isDown = TRUE ; cds = TRUE ;
	      break ;
	    case 3:
	      suffix = ".r.cds_sponge" ;
	      isDown = FALSE ; cds = TRUE ;
	      break ;
	    }
	  ac_free (ao) ;
	  ao = aceOutCreate (sx->outFileName, suffix, sx->gzo, h) ;
	  if (arrayMax (shadows))
	    {
	      for (old = i = j = 0, shadow = arrp(shadows, i, SHADOW); i <= arrayMax (shadows) ; i++, shadow++)
		{
		  if (! shadow->mrna || ! shadow->gene || ! shadow->target || shadow->isDown != isDown)
		    continue ;
		  if (! sx->gffBacteria && shadow->cds != cds) 
		    continue ;
		  if (shadow->mrna != old)
		    j = 0 ;
		  j++ ; old = shadow->mrna ;
		  aceOutf (ao, "%s\t%d\t%s\t%d\t%d\t%s\n"
			   , dictName(dict, shadow->mrna)
			   , j 
			   , dictName(sx->selectDict, shadow->target)
			   , shadow->a1
			   , shadow->a2
			   , dictName(dict, shadow->gene)
			   ) ;
		}
	    }
	  ac_free (ao) ;
	}
    }
  
  ac_free (ao) ;
  ao = aceOutCreate (sx->outFileName, ".transcripts.ace", sx->gzo, h) ;

  aceOutDate (ao, "//", sx->title) ;

  if (0 && arrayMax (shadows))
    {
      for (old = target = i = 0, shadow = arrp(shadows, i, SHADOW); i <= arrayMax (shadows) ; i++, shadow++)
	{
	  if (i < arrayMax (shadows) && !  shadow->gene)
	    continue ;

	  if (i == arrayMax (shadows) || old != shadow->gene)
	    { 
	      if (i < arrayMax (shadows))
		{
		  old = shadow->gene ; 
		  if (0 && sx->isGff3 && shadow->GeneID)
		    old =  shadow->GeneID ;
		  target = shadow->target ;
		  a1 = shadow->a1 ; a2 = shadow->a2 ;
		}
	    }
	  if (i < arrayMax (shadows))
	    {
	      if (a1 < a2)
		{
		  if (a1 > shadow->a1) a1 = shadow->a1 ;
		  if (a2 < shadow->a2) a2 = shadow->a2 ;
		}
	      else
		{
		  if (a1 < shadow->a1) a1 = shadow->a1 ;
		  if (a2 > shadow->a2) a2 = shadow->a2 ;
		}
	    }
	}
    }

  if (arrayMax (shadows))
    {
      BOOL ok2 = FALSE ;
      if (0) showShadow (sx, dict, shadows) ;
      
      for (old = target = i = 0, shadow = arrp(shadows, i, SHADOW); i <= arrayMax (shadows) ; i++, shadow++)
	{
	  if (i < arrayMax (shadows) && !  shadow->mrna)
	    continue ;
	  if (i == arrayMax (shadows) || old != shadow->mrna)
	    { 
	      if (old)
		{
		  int dxCDS = 0 ;
		  int gene, title, note, gene_name, locus_tag ;
		  aceOutf (ao, "Sequence %s\nSubsequence %s %d %d\n\n"
			   , dictName(sx->selectDict, target)
			   , dictName(dict, old)
			   , a1, a2
			   ) ;
		  
		  aceOutf (ao, "Sequence %s\nIntMap %s %d %d\n"
			   , dictName(dict, old)
			   , dictName(sx->selectDict, target)
			   , a1, a2 
			   ) ;
		  aceOutf (ao, "Method Genefinder\nIs_predicted_gene\n") ;
		  if (sx->gffBacteria)
		    {
		      int dx = a2 - a1 ;
		      if (dx < 0) dx = -dx ;
		      dx++ ;
		      if (0) aceOutf (ao, "-D Source_exons\n") ;
		      aceOutf (ao, "Source_exons 1 %d CDS\n", dx) ;
		    }

		  for (j = i - 1, shadow2 = shadow - 1 ; j >= 0 &&  (! shadow2->mrna || shadow2->mrna == old) ; j--, shadow2--)  {} ;
		  j++ ; shadow2++ ; 
		  for ( gene = title = note = hasCDS = hasEXON = gene_name = locus_tag = 0 ; j < i ; j++, shadow2++)
		    {
		      if (! shadow2->mrna)
			continue ;
		      if (shadow2->ncRNA)
			aceOutf (ao, "%s\n", dictName(dict, shadow2->ncRNA)) ;
		      if (shadow2->gene && ! gene)
			{
			  int gene2 ;
			  gene2 = gene = shadow2->gene ;
			  aceOutf (ao, "Model_of \"X__%s\"\n", dictName(dict, gene2)) ;
			  aceOutf (ao, "Model_of_gene \"%s\"\n", dictName(dict, gene2)) ;
			  if (shadow2->GeneID) 
			    {
			      aceOutf (ao, "GeneId_pg %s\n", dictName(dict, shadow2->GeneID)) ;
			      aceOutf (ao, "GeneId %s\n", dictName(dict, shadow2->GeneID)) ;
			    }
			}
		      if (! shadow2->cds && shadow2->gene_title && ! title)
			{
			  title = shadow2->gene_title ;
			  aceOutf (ao, "Title \"%s"
				   ,  dictName(dict, title)
				   ) ;
			  {
			    char *coma, *buf = strnew (dictName(dict, title), 0) ;
			    coma = strstr (buf, ", transcript variant") ;
			    if (coma) *coma = 0 ;
			    ac_free (buf) ;
			  }
			  if (shadow2->gene_biotype && shadow2->gene_biotype != Exon_type)
			    aceOutf (ao, ", %s"
				     , dictName(dict, shadow2->gene_biotype)
				     ) ;
			  aceOut (ao, "\"\n") ;
			} 
		      if (! shadow2->cds && shadow2->gene_name && ! gene_name)
			{
			  gene_name = shadow2->gene_name ;
			  aceOutf (ao, "LocusLink \"%s\"\n"
				   ,  dictName(dict, gene_name)
				   ) ;
			}
		      if (1 && shadow2->locus_tag && !locus_tag)
			{
			  locus_tag = shadow2->locus_tag ;
			  aceOutf (ao, "Locus \"%s\"\n",dictName (dict, shadow2->locus_tag)) ;
			}

		      if (sx->gffBacteria && shadow2->note && ! note)
			{
			  note = shadow2->note ;
			  aceOutf (ao, "Concise_description %s\n "
				   ,  ac_protect (dictName(dict, note), h)
				   ) ;
			} 
		      if (shadow2->gene_biotype && shadow2->gene_biotype != Exon_type)
			{
			  aceOutf (ao, "%s\n"
				   , dictName(dict, shadow2->gene_biotype)
				   ) ;
			}
		    }
		  for (j = i - 1, shadow2 = shadow - 1 ; j >= 0 &&  (! shadow2->mrna || shadow2->mrna == old) ; j--, shadow2--)  {} ;
		  j++ ; shadow2++ ; 
		  {
		    int v1 = 0 ; /*  v2 = 0 ; */
		    hasCDS = FALSE ;
		    for ( ; j < i && ! hasCDS ; j++, shadow2++)
		      {
			if (! shadow2->mrna)
			  continue ;
			if (shadow2->exon)
			  hasEXON = TRUE ;
			if (shadow2->cds && ! hasCDS)
			  {
			    hasCDS = TRUE ;
			    aceOutf (ao, "CDS // a=%d %d   x=%d %d\n", shadow2->a1, shadow2->a2, shadow2->x1, shadow2->x2) ; 
			    v1 = shadow2->a1 ;
			    /* v2 = shadow2->a2 ;   v1, v2 first cds */
			  }
		      }
		    if (hasCDS) /* start again on the exon to find the start of the CDS */
		      {
			for (j = i - 1, shadow2 = shadow - 1 ; j >= 0 &&  (! shadow2->mrna || shadow2->mrna == old) ; j--, shadow2--)  {} ;
			j++ ; shadow2++ ; 
			for ( ; j < i ; j++, shadow2++)
			  {
			    if (! shadow2->mrna)
			      continue ;
			    if (shadow2->cds) continue ;
			    if (shadow2->isDown)
			      {
				if (shadow2->a2 < v1)
				  dxCDS += shadow2->x2 - shadow2->x1 + 1 ;
				else
				  {
				    dxCDS += v1 - shadow2->a1 ;
				    break ;
				  }
			      }
			    else
			      {
				if (shadow2->a2 > v1)
				  dxCDS += shadow2->x2 - shadow2->x1 + 1 ;
				else
				  {
				    dxCDS += shadow2->a1 - v1 ;
				    break ;
				  }
			      }
			  }
		      }
		  }
		  
		  for (j = i - 1, shadow2 = shadow - 1 ; j >= 0 &&  (! shadow2->mrna || shadow2->mrna == old) ; j--, shadow2--)  {} ;
		  j++ ; shadow2++ ; 
		  for ( ; j < i ; j++, shadow2++)
		    { int dv= 0 ; /* set to gene start to help debugging i.e. dv = 4195563 ; */
		      int jCds, b1, b2, u1, u2, v1, v2 ;
		      SHADOW *cds ;
		      BOOL ok1 ;
		      
		      if (shadow2->cds)
			continue ;
		      if (! shadow2->mrna)
			continue ;
		      
		      /* if there is a CDS inside this exon it is below, thanks to shadowOrder */
		      ok1 = FALSE ;
		      if (hasCDS)
			{ 
			  for(jCds = j + 1, cds = shadow2 + 1 ; ! ok1 && jCds < i ; cds++, jCds++)
			    {
			      if (! cds->cds)
				continue ;

			      b1 = cds->x1 + dxCDS > shadow2->x1 ? cds->x1 + dxCDS : shadow2->x1 ;
			      b2 = cds->x2 + dxCDS < shadow2->x2 ? cds->x2 + dxCDS : shadow2->x2 ;
			      if (b1 <= b2)
				{
				  ok1 = TRUE ;
				  if (! ok2 && b1 > shadow2->x1)
				    {
				      u1 = shadow2->a1 ; u2 = shadow2->a2 ;
				      v1 = shadow2->isDown ? u1 - a1 + 1 : a1 - u1 + 1 ;
				      v2 = shadow2->isDown ? u2 - a1 + 1 : a1 - u2 + 1 ;
				      v2 = v1 + b1 - shadow2->x1 - 1 ;
				      aceOutf (ao, "Source_exons %d %d %s\n" 
					       , v1 + dv , v2+ dv
					       , "UTR_5prime"
					       ) ;
				    }
				  if (1)
				    {
				      ok1 = ok2 = TRUE ;
				      u1 = shadow2->a1 ; u2 = shadow2->a2 ;
				      v1 = shadow2->isDown ? u1 - a1 + 1 : a1 - u1 + 1 ;
				      v2 = shadow2->isDown ? u2 - a1 + 1 : a1 - u2 + 1 ;
				      v2 = v1 + b2 - shadow2->x1 ;
				      v1 = v1 + b1 - shadow2->x1 ;
				      aceOutf (ao, "Source_exons %d %d %s\n"
					       , v1+ dv, v2+ dv
					       , "CDS"
					       ) ;
				    }
				  if (b2 < shadow2->x2)
				    {
				      u1 = shadow2->a1 ; u2 = shadow2->a2 ;
				      v1 = shadow2->isDown ? u1 - a1 + 1 : a1 - u1 + 1 ;
				      v2 = shadow2->isDown ? u2 - a1 + 1 : a1 - u2 + 1 ;
				      v1 = v1 + b2 - shadow2->x1 + 1 ;
				      aceOutf (ao, "Source_exons %d %d %s\n" 
					       , v1+ dv, v2+ dv
					       , "UTR_3prime"
					   ) ;
				    }
				}
			    }
			}
		      
		      if (! ok1)
			{
			  u1 = shadow2->a1 ; u2 = shadow2->a2 ;
			  v1 = shadow2->isDown ? u1 - a1 + 1 : a1 - u1 + 1 ;
			  v2 = shadow2->isDown ? u2 - a1 + 1 : a1 - u2 + 1 ;
			  aceOutf (ao, "Source_exons %d %d %s\n" 
				   , v1+ dv, v2 + dv
				   , hasCDS ? (ok2 ? "UTR_3prime" : "UTR_5prime") : "Exon"
				   ) ;
			}
		    }
		  
		  aceOutf (ao, "\n") ;
		}
	      ok2 = FALSE ;
 	    }

	  if  (i < arrayMax (shadows) && ! sx->gffBacteria && shadow->cds) continue ;
	  if (i < arrayMax (shadows) && old != shadow->mrna)
	    {
	      old = shadow->mrna ;
	      target = shadow->target ;
	      a1 = shadow->a1 ; a2 = shadow->a2 ;
	    }

	  
	  if (i < arrayMax (shadows))
	    {
	      if (a1 < a2)
		{
		  if (a1 > shadow->a1) a1 = shadow->a1 ;
		  if (a2 < shadow->a2) a2 = shadow->a2 ;
		}
	      else
		{
		  if (a1 < shadow->a1) a1 = shadow->a1 ;
		  if (a2 > shadow->a2) a2 = shadow->a2 ;
		}
	    }
	}

      if (1)
	{
	  int oldGene = 0, oldID = 0 ;
	  
	  for (i = 0, shadow = arrp(shadows, i, SHADOW); i < arrayMax (shadows) ; i++, shadow++)
	    if (shadow->gene && shadow->GeneID)
	      {
		if (oldGene != shadow->gene || oldID != shadow->GeneID)
		  aceOutf (ao, "LocusLink \"%s\"\nGeneId \"%s\"\n\n"
			   , dictName(dict, shadow->gene)
			   , dictName(dict, shadow->GeneID)
			   ) ;
		oldGene = shadow->gene ;
		oldID = shadow->GeneID ;
	      }
	}
      
      
 
      if (1)
	for (old = i = 0, shadow = arrp(shadows, i, SHADOW); i + 1 < arrayMax (shadows) ; i++, shadow++)
	  {
	    shadow2 = shadow + 1 ;
	    a1 = shadow->a2 ;
	    a2 = shadow2->a1 ;
	    isDown = shadow->isDown ;
	    if (0) continue ; /* 2024_07_22 */
	    if (shadow->mrna  && shadow->mrna == shadow2->mrna &&
		! shadow->cds && ! shadow2->cds
		)
	      {
		int b1, b2 ;
		if (shadow->isDown)
		  { b1 = a1 + 1 ; b2 = a2 - 1 ; if (b2 - b1 < 0) continue ; }
		else 
		  { b1 = a1 - 1 ; b2 = a2 + 1 ; if (b2 - b1 > 0) continue ;}
		
		aceOutf (ao, "Intron %s__%d_%d\n" , dictName(sx->selectDict, shadow->target), b1, b2) ;
		aceOutf (ao, "Length %d\n", b2 > b1 ? b2 - b1 + 1 : b1 - b2 + 1) ;
		aceOutf (ao, "IntMap %s %d %d\n" , dictName(sx->selectDict, shadow->target), b1, b2) ;
		aceOutf (ao, "Gene %s\nFrom_genefinder %s\n\n"
			 , ac_protect (dictName(dict, shadow2->gene), h)
			 , ac_protect (dictName(dict, shadow2->mrna),h)
			 ) ;
	      }
	  }
    }
  ac_free (ao) ;
 done:   
  ac_free (ai) ;
  ac_free (h) ;
  return ;
}  /* parseGtfFeature */

/**********************************/

static void parseGtfFile (SX *sx) 
{
  KEYSET mrna2gene = keySetCreate () ;
  /* CDS must be called first to allow the definition of the UTRs by substracting the CDS from the exons */
  sx->gnam2genes = arrayHandleCreate (20000, KEYSET, sx->h) ;
  if (! sx->isGff3) /* i.e. when parsing a gtf file */
    {
      /* must come first because we may edit the type when parsing mrnaRemap */
      if (1) parseGtfFeature (sx, "pre_miRNA", ".pre_mirnaRemap", 1, FALSE,  mrna2gene) ; /* mrna */
      if (1) parseGtfFeature (sx, "miRNA", ".mirnaRemap", 1, FALSE,  mrna2gene) ; /* mrna */
      if (1) parseGtfFeature (sx, "rRNA", ".rRNA", 1, FALSE,  mrna2gene) ; /* mrna */
      if (1) parseGtfFeature (sx, "tRNA", ".tRNA", 1, FALSE,  mrna2gene) ; /* mrna */
      if (1) parseGtfFeature (sx, "ncRNA", ".ncRNA", 1, FALSE,  mrna2gene) ; /* mrna */

      /* edited for mouse/rat NCBI.gff3 2016_09_15 */
      /* edited for human T2T gtf 2025_07_19 */
      if (sx->gffBacteria)
	{
	  parseGtfFeature (sx, "mrna",  ".cdsRemap",  0, TRUE,  mrna2gene) ; /* mrna && CDS */
	  parseGtfFeature (sx, "mrna",  ".mrnaRemap",  1, FALSE,  mrna2gene) ; /* mrna && CDS */
	}
      else
	{
	  if (1) parseGtfFeature (sx, "CDS",  ".cdsRemap",  0, TRUE,  mrna2gene) ; /* CDS */
	  if (1) parseGtfFeature (sx, "exon", ".mrnaRemap", 1, FALSE,  mrna2gene) ; /* mrna */
	}
      if (1) 
	parseGtfFeature (sx, "gene",  ".geneTitle.ace",  0, FALSE,  mrna2gene) ; /* gene -> title */
    }
  else /* parse a gff3 file */ 
    {   /* written for the E.coli K12 gff dump from NCBY june 1 2015 
	 * works for D.melano 2016_08
	 * does not work for E coli K12 oct 2017
	 */
      sx->rnaId2mrna = keySetHandleCreate (sx->h) ;
      if (sx->gffBacteria)
  	{
	  parseGtfFeature (sx, "mrna",  ".cdsRemap",  0, TRUE,  mrna2gene) ; /* mrna && CDS */
	  parseGtfFeature (sx, "mrna",  ".mrnaRemap",  1, FALSE,  mrna2gene) ; /* mrna && CDS */
	}
      else
	{
	  if (sx->gffTAIR)
	    {
	      parseGtfFeature (sx, "mRNA", 0, 0, FALSE,  mrna2gene) ; /* mrna */
	    }
	  parseGtfFeature (sx, "cds",  ".cdsRemap",  0, TRUE,  mrna2gene) ; /* CDS */
	  parseGtfFeature (sx, "exon", ".mrnaRemap", 1, FALSE, mrna2gene) ; /* mrna */
	}
    }
  ac_free (mrna2gene) ;
} /* parseGtfFile */

/*************************************************************************************/
/* gather the geometry off the shadow file */
static int parseRefGeneFile (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (sx->refGeneFileName, FALSE, h) ;
  int i, nn = 0, dummy ;
  const char *ccp ;
  int m1, m2, c1, c2, ax1[10000], ax2[10000],x0,  nexons, mrna = 0, gene = 0, target = 0 ;
  BOOL isDown = TRUE ;
  char cutter ;
  SHADOW *up ;


  /* genepred format, UCSC
genename
gene
chrom
strand +-
trtans strat 
trend
cds starts
cds end
nbexon
exonstrat[]
exonend[]
  */


  sx->selectDict = dictHandleCreate (10000, sx->h) ;
  sx->shadowDict = dictHandleCreate (10000, sx->h) ;
  sx->shadowArray = arrayHandleCreate (10000, SHADOW, sx->h) ;
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      if (aceInInt (ai, &dummy))  /* index nin, drop it */
	aceInStep (ai, '\t') ; 
      ccp = aceInWordCut (ai, "\t", &cutter) ;
      if (! ccp || *ccp == '#')
	continue ;

      dictAdd (sx->shadowDict, ccp, &mrna) ; 
      aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ;
      dictAdd (sx->selectDict, ccp, &target) ;
      aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ;
      if (ccp && *ccp == '+') isDown = TRUE ;
      else if (ccp && *ccp == '-') isDown = FALSE ;
      else continue ;

      aceInStep (ai, '\t') ; aceInInt (ai, &m1) ;  /* mrna */
      aceInStep (ai, '\t') ; aceInInt (ai, &m2) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &c1) ;   /* coding */
      aceInStep (ai, '\t') ; aceInInt (ai, &c2) ;

      aceInStep (ai, '\t') ; aceInInt (ai, &nexons) ;
      aceInStep (ai, '\t') ;
      memset (ax1, 0, nexons * sizeof(int)) ;
      memset (ax2, 0, nexons * sizeof(int)) ;
      for (i = 0 ; i < nexons ; i++)
	{ 
	  ccp = aceInWordCut (ai, ",", 0) ;
	  if (ccp)
	    {
	      ax1[i] = atoi (ccp) ;
	      aceInStep (ai, ',') ;
	    }
	}
      aceInStep (ai, '\t') ;		       
      for (i = 0 ; i < nexons ; i++)
	{ 
	  ccp = aceInWordCut (ai, ",", 0) ;
	  if (ccp)
	    {
	      ax2[i] = atoi (ccp) ;
	      aceInStep (ai, ',') ;
	    }
	}

      aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; /* a score, drop it */  
      aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; /* a score, drop it */  
      dictAdd (sx->shadowDict, ccp, &gene) ; 
      if (gene && target)
	{
	  if (isDown)
	    {
	      for (i = x0 = 0 ; i < nexons ; i++)
		{
		  up = arrayp (sx->shadowArray, nn++, SHADOW) ;
		  up->mrna = mrna ; up->gene = gene ; up->target = target ;
		  up->a1 = ax1[i] + 1 ; up->x1 = ++x0 ;
		  up->a2 = ax2[i] + 0 ; up->x2 = x0 = x0 + ax2[i] - ax1[i]  - 1 ;
		}
	    }
	  else
	    {
	      for (i = nexons - 1, x0 = 0 ; i >= 0 ; i--)
		{
		  up = arrayp (sx->shadowArray, nn++, SHADOW) ;
		  up->mrna = mrna ; up->gene = gene ; up->target = target ;
		  up->a2 = ax1[i] + 1 ; up->x1 = ++x0 ;
		  up->a1 = ax2[i] + 0 ; up->x2 = x0 = x0 + ax2[i] - ax1[i]  - 1 ;
		}
	    }
	  gene = target = 0 ;
	}
    }
  ac_free (ai) ;
  ac_free (h) ;

  arraySort (sx->shadowArray, shadowOrder) ;
  arrayCompress (sx->shadowArray) ;
  return nn ;
} /* parseRefGeneFile */

/*************************************************************************************/
/* export the relevant subsequences of the current sequence */
static int exportShadowSequences (SX *sx, Array aa)
{
  AC_HANDLE h = ac_new_handle () ;
  Array cdna = 0 ;  
  char cc, *cp ;
  int nn = 0, ii, iMax = arrayMax (aa), j = 0, mrna = 0, gene = 0, geneId = 0, target = 0 ;
  static int lastii = 0 ;
  SHADOW *up ;
  vTXT id = sx->id ;
  vTXT id1 = vtxtHandleCreate (h) ;
  vTXT dna = sx->dna ;
  vTXT dna1 = vtxtHandleCreate (h) ;
  int len = vtxtLen (dna) ;

  dictFind (sx->selectDict, vtxtPtr (id), &target) ;
  sx->id = id1 ;
  sx->dna = dna1 ;
  if (lastii >= arrayMax(aa))
    { lastii = arrayMax(aa) ; lastii-- ; }
  if (lastii < 0) 
    lastii = 0 ;
  for (ii = lastii, up = arrp (aa, ii, SHADOW) ; ii >= 0 ; ii--, up--)
    {
      if (target > up->target) break ;
    }
  if (ii < 0) ii = 0 ;
  for (up = arrp (aa, ii, SHADOW) ; ii < iMax ; ii++, up++)
    {
      if (up->cds && ! sx->gffBacteria) continue ;
      if (up->mrna != mrna)
	{
	  if (mrna)
	    {
	      vtxtClear (id1) ;
	      vtxtPrint (id1, dictName (sx->shadowDict, mrna)) ;
	      if (gene)
		vtxtPrintf (id1, "|Gene|%s", dictName (sx->shadowDict, gene)) ;
	      if (geneId)
		vtxtPrintf (id1, "|GeneId|%s", dictName (sx->shadowDict, geneId)) ;
	      if (vtxtPtr (sx->dna))
		exportSequence (sx) ;
	      vtxtClear (dna1) ;
	      if (cdna)
		arrayMax (cdna) = 0 ;
	    }
	  mrna = up->mrna ;
	  gene = up->gene ;
	  geneId = up->GeneID ;
	}
      if (target != up->target) continue ;
      if (target < up->target) break ;
      
      if (up->a1 > len) up->a1 = len ;
      if (up->a2 > len) up->a2 = len ;
      if (up->a1 <= up->a2)
	{
	  cp = vtxtPtr (dna) + up->a2 ;
	  cc = *cp ; *cp = 0 ;
	  vtxtPrint (dna1, vtxtPtr (dna) + up->a1 - 1) ;
	  *cp = cc ; 
	}
      else
	{
	  if (! cdna)
	    cdna = arrayHandleCreate (10000, char, h) ;
	  cp = vtxtPtr (dna) + up->a1 ;
	  cc = *cp ; *cp = 0 ;
	  j = up->a1 - up->a2 + 1 ;
	  array (cdna, j, char) = 0 ; arrayMax (cdna) = j ;
	  memcpy (arrp (cdna, 0, char), vtxtPtr (dna) + up->a2 - 1, j) ;
	  *cp = cc ; 
	  dnaEncodeArray (cdna) ;
	  reverseComplement (cdna) ;
	  dnaDecodeArray (cdna) ;
	  vtxtPrint (dna1, arrp (cdna, 0, char)) ;
	}
    }

  if (mrna)
    {
      vtxtClear (id1) ;
      vtxtPrint (id1, dictName (sx->shadowDict, mrna)) ;
      if (gene)
	vtxtPrintf (id1, "|Gene|%s", dictName (sx->shadowDict, gene)) ;
      if (geneId)
	vtxtPrintf (id1, "|GeneId|%s", dictName (sx->shadowDict, geneId)) ;

      if (vtxtPtr (sx->dna))
	exportSequence (sx) ;
    }
  lastii = ii ;
  sx->id = id ;
  sx->dna = dna ;

  return nn ;
} /* exportShadowSequences */

/*************************************************************************************/
/*************************************************************************************/

int tagCountOrder (const void *a, const void *b)
{
  int x1 = ((const SEQ *) a)->count, x2 = ((const SEQ *) b)->count ;
  return x1 < x2 ? 1 : (x1 > x2 ? -1 : 0) ;
}

/*************************************************************************************/
static DICT *myDict = 0 ;

int atgcOrder (const void *a, const void *b)
{
  int n1 = ((const SEQ *) a)->dna ;
  int n2 = ((const SEQ *) b)->dna ;
  
  return strcmp (dictName (myDict, n1), dictName (myDict, n2)) ;
} /* atgcOrder */

/*************************************************************************************/

static void exportTagCountDo (ACEOUT ao, SX *sx, SEQ *seq)
{
  int ii ;
  static int oldDna = 0 ;
  static char *buf = 0 ;
  char cc, *cp ;
  const char *ccp, *ccq ;
  static int bufLn = 0 ;
  
  /* 2011_06 xor is not documented and was probably a transitory idea */
  if (sx->xor == 1)
    {
      ii = 0 ;
      ccp = dictName (sx->dict, seq->dna) ;
      if (oldDna)
	{
	  ccq = dictName (sx->dict, oldDna) ;
	  while (*ccp++ == *ccq++) ii++ ;
	}  
      if (seq->count)
	aceOutf (ao, "%d\t%s\t%d", ii,dictName (sx->dict, seq->dna) + ii, seq->count) ;
      oldDna = seq->dna ;
    }
  else if (sx->xor == 3)
    {
      ccp = dictName (sx->dict, seq->dna) ;
      if ((ii = strlen (ccp)) > bufLn)
	{
	  cp = buf ;
	  ii <<= 1 ;
	  buf = messalloc(ii) ;
	  bufLn = ii ;
	}
      strcpy (buf, ccp) ;
      if (oldDna)
	{
	  cp =buf ;
	  ccp = dictName (sx->dict, oldDna) ;
	  while (*cp | *ccp)
	    {
	      cc = (*cp - 'a') ^ (*ccp - 'a') ;
	      *cp++ = 'a' + cc ; ccp++ ;
	    }
	  *cp = 0 ;
	}  
      if (seq->count)
	aceOutf (ao, "%s\t%d", buf, seq->count) ;
      oldDna = seq->dna ;
    }
  else
    {
      if (seq->count)
	aceOutf (ao, "%s\t%d", dictName (sx->dict, seq->dna), seq->count) ;
    }
  if (sx->keepName)
    {
      aceOutf (ao, "\t%s", vtxtPtr (seq->id)) ;

    }

  aceOutf (ao, "\n") ;
  return ;
} /* exportTagCountDo */

/***********/

static void exportTagCount (SX *sx)
{
  int mx = arrayMax (sx->tagArray) ;
  int ii, ns, nt ;
  float nMb = 0 ;
  SEQ *seq = 0 ;
  DICT *dict = sx->dict ;

  /* sort in sequence order, this may accelerate the alignements and optimize sipping */
  if (1)
    {
      myDict = sx->dict ;
      if (1) arraySort (sx->tagArray, atgcOrder) ;
    }
  else
    arraySort (sx->tagArray, tagCountOrder) ;
  for (ii = ns = nt = 0, seq = arrp (sx->tagArray, 0, SEQ) ; ii < mx ; ii++, seq++)
    if (seq->count)
      {
	if (sx->splitByPrefix)
	  {
	    selectOutFile (sx, dictName (dict, seq->dna)) ;
	    if (seq->count)
	      exportTagCountDo (sx->ao, sx, seq) ;
	  }
	else
	  {
	    exportTagCountDo (sx->ao, sx, seq) ;
	    nMb += seq->dnaLn/((double)1000000.0) ;
	    ns++ ; nt += seq->count ;
	    if (
		(sx->split && (ns % sx->split) == 0) ||
		(sx->splitMb && (sx->splitMb < nMb))
		)
	      { /* rotate the output file */
		char *ccp ;
		
		nMb = ns = 0 ;
		ac_free (sx->ao) ;
		ccp =  hprintf (sx->h, "%s.%d.tc%s", sx->outFileName, ++(sx->nSplit), sx->gzo ? ".gz" : "") ;
		if (sx->gzo)
		  sx->ao = aceOutCreateToGzippedFile (ccp, sx->h) ;
		else
		  sx->ao = aceOutCreateToFile (ccp, "wb", sx->h) ;

		if (sx->aoId)
		  {
		    ac_free (sx->aoId) ;
		    ccp =  hprintf (sx->h, "%s.%d.id", sx->outFileName, ++(sx->nSplit), sx->gzo ? ".gz" : "") ;
		    if (sx->gzo)
		      sx->aoId = aceOutCreateToGzippedFile (ccp, sx->h) ;
		    else
		      sx->aoId = aceOutCreateToFile (ccp, "wb", sx->h) ;
		  }
	      }
	  }
      }

  return ;
} /* exportTagCount */

/*************************************************************************************/
/*************************************************************************************/
static void sxLetterDistribGooglePlot (SX* sx, long int *bb)
{
  AC_HANDLE h = ac_new_handle () ;
  char *cp ;
  int j ;
  ACEOUT ao = aceOutCreate (sx->outFileName, ".html", 0, h) ;

  aceOutf (ao,
		"<html>\n"
		" <head>\n"
		"    <!--Load the AJAX API-->\n"
		"    <script type=\"text/javascript\" src=\"https://www.google.com/jsapi\"></script>\n"
		"    <script type=\"text/javascript\">\n"
		"    "
		"      // Load the Visualization API and the piechart package.\n"
		"      google.load('visualization', '1.0', {'packages':['corechart']});\n"
		"      "
		"      // Set a callback to run when the Google Visualization API is loaded.\n"
		"      google.setOnLoadCallback(drawChart);\n"
		"      "
		"      // Callback that creates and populates a data table, \n"
		"      // instantiates the pie chart, passes in the data and\n"
		"      // draws it.\n"
		"      function drawChart() {\n"
		" "
		"      // Create the data table.\n"
		"      var data = new google.visualization.DataTable();\n"
		"      data.addColumn('string', 'Topping');\n"
		"      data.addColumn('number', 'Slices');\n"
		"      data.addRows([\n"
		) ;

  if (! sx->qualityPlot)
    {
      for (cp = "atgcn" ; *cp ; cp++)
	{
	  j = *cp ;
	  aceOutf (ao,"['%c',%ld]%s\n", *cp, 25+bb[j], j != 'n' ? "," : "") ;
	}
    }
  
  aceOut (ao, 
		"      ]);\n"
		"      "
		"      // Set chart options\n"
		"      var options = {'title':'All other letters present in these sequences',\n"
		"                     'width':400,\n"
		"                     'height':300};\n"
		"       "
		"      // Instantiate and draw our chart, passing in some options.\n"
		"      var chart = new google.visualization.PieChart(document.getElementById('chart_div'));\n"
		"      chart.draw(data, options);\n"
		"    }\n"
		"    </script>\n"
		"  </head>\n"
		"  <body>\n"
		"    <!--Div that will hold the pie chart-->\n"
		"    <div id=\"chart_div\"></div>\n"
		"  </body>\n"
		"</html>\n"
		) ;
  aceOutf (ao,"\n") ;
  ac_free (ao) ;

  ac_free (h) ; /* will close */
} /* sxLetterDistribGooglePlot */

/*************************************************************************************/
/* report 
 * ns: number of sequences of given length
 * b [c] : number of character c seen at any position
 * bb[256*c + i] : number of character c at position i
 */
static  void sxLetterDistribReport (SX* sx)
{
  long int *ns = sx->letterLength ;
  long int *b = sx->letterCount ;
  long int *bb = sx->letterProfile ;
  int NN = 2 * sx->letterNN + 4 ;
  long int n = 0, nn = 0, nTotal = 0 ;
  char cc, *cp, *titre ;
  int i, j, iMax = 0, di, iMax1 = 0 ;
  ACEOUT ao = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  int jmin, jmax ;
  BOOL showTotal = FALSE ;
	
  titre = hprintf (h, "%s", sx->title ? sx->title : "s") ;
  cp = strchr (titre, ' ') ;
  if (cp) *cp = 0 ;
  
  for (jmin = 1 ; jmin < 256 && b[jmin] == 0 ; jmin++) ;
  for (jmax = 255 ; jmax > 1 && b[jmax] == 0 ; jmax--) ;

  ao = aceOutCreate (sx->outFileName, ".profile.txt", FALSE, h) ;
  aceOutDate (ao, "###", sx->title) ;

  for (i = NN  ; i >= 0 ; i--)
    if (ns[i] > 0)
      { iMax = iMax1 = i ; break ; }
  if (iMax < 30)  /* always analyse at least 30 bases, this makes the report format more regular */ 
    iMax = iMax1 = 30 ;
  if (iMax > sx->letterNN)
    for (i = sx->letterNN  ; i >= 0 ; i--)
      if (ns[i] > 0)
	{ iMax1 = i ; break ; }
  aceOutf (ao,"# Position\tNumber of reads\tRun") ;
  if (! sx->qualityPlot)
    {
      for (cp = "ATGCN" ; *cp ; cp++)
	aceOutf (ao,"\t%c", *cp) ;
      aceOutf (ao,"\tRun") ;
      for (cp = "ATGCN" ; *cp ; cp++)
	aceOutf (ao,"\t%%%c", *cp) ;
    }
  else
    {
      for (j = jmin ; j <= jmax ; j++)
	aceOutf (ao,"\t%d", j) ;   
      aceOutf (ao,"\tRun") ;
      for (j = jmin ; j <= jmax ; j++)
	aceOutf (ao,"\t%%%d", j) ;
    }
  for (i = di = 0 ; i < iMax ; i++)
    {
      if (i >= iMax1 && di == 0)
	{
	  i = sx->letterNN + 1 ;
	  di = sx->letterNN + 2 ;
	  aceOut (ao,"\n\n") ;
	  continue ;
	}
      aceOutf (ao,"\n%d", i + 1 - di) ;
      aceOutf (ao,"\t%d", ns[i+1]) ; nTotal += (long int) ns[i+1] ;

      if (i - di == 0) aceOutf (ao,"\t%s/%d", titre, di == 0 ? 1 : 2) ;
      else if ((i - di) % 50 == 49) aceOutf (ao,"\t%d", i - di + 1) ;
      else aceOutf (ao,"\t") ;

      if (! sx->qualityPlot)
	{
	  for (n = 0, cp = "atgcn" ; *cp ; cp++)
	    {
	      j = *cp ;
	      n += bb[256*i+j] ;
	    }
	  nn += n ;
	  if (n == 0) n = 1 ;
	  for (cp = "atgcn" ; *cp ; cp++)
	    {
	      j = *cp ;
	      aceOutf (ao,"\t%d", bb[256*i+j]) ;
	    }

	  if (i - di == 0) aceOutf (ao,"\t%s/%d", titre, di == 0 ? 1 : 2) ;
	  else if ((i - di) % 50 == 49) aceOutf (ao,"\t%d", i - di + 1) ;
	  else aceOutf (ao,"\t") ;

	  for (cp = "atgcn" ; *cp ; cp++)
	    {
	      j = *cp ;
	      aceOutf (ao,"\t%.1f", (100.0*bb[256*i+j])/n) ;
	    }
	}
      else
	{
	  int cumul = 0 ;
	  int cumul2[jmax] ;	
	  float zcumul = 0 ;

	  memset (cumul2, 0, sizeof(cumul2)) ;
	  n = 0 ;
	  for (j = jmax ; j >= jmin ; j--)
	    {
	      n += bb[256*i+j] ;
	      cumul2[j] = n ;
	    }
	  nn += n ;
	  if (n == 0) n = 1 ;
	  for (j = jmin ; j <= jmax ; j++)
	    {
	      cumul = cumul2[j] ;
	      aceOutf (ao,"\t%d", cumul) ;
	    }

	  if (i - di == 0) aceOutf (ao,"\t%s/%d", titre, di == 0 ? 1 : 2) ;
	  else if ((i - di) % 50 == 49) aceOutf (ao,"\t%d", i - di + 1) ;
	  else aceOutf (ao,"\t") ;

	  zcumul = 0 ;
	  for (j = jmin ; j <= jmax ; j++)
	    {
	      zcumul =  (100.0*cumul2[j])/n ;
	      aceOutf (ao,"\t%.1f", zcumul) ;
	    }
	}
    }
  if (showTotal) aceOutf (ao,"\nTotal\t%ld", nTotal) ;
  if (nn == 0) nn = 1 ;
  if (! sx->qualityPlot)
    {
      if (showTotal)
	{
	  for (cp = "atgcn" ; *cp ; cp++)
	    {
	      j = *cp ;
	      aceOutf (ao,"\t%ld", b[j]) ;
	    }
	  for (cp = "atgcn" ; *cp ; cp++)
	    {
	      j = *cp ;
	      aceOutf (ao,"\t%.1f", (100.0*b[j])/nn) ;
	    }
	}
      aceOutf (ao,"\n\nAll other letters present in these sequences\nASCII code\tChar\tNumber of occurences") ;
      nn = 0 ;
      for (j = 1 ; j<256 ; j++)
	{
	  if (!b[j]) continue ;
	  cc = (char) j ;
	  for (cp = "atgcn" ; *cp ; cp++)
	    if (cc == *cp) break ;
	  if (*cp) continue ;
	  nn += b[j] ;
	  aceOutf (ao,"\n%d\t%c\t%d", j, ace_upper (cc), b[j]) ;
	}
      if (!nn) aceOutf (ao,"\nNothing except ATGCN") ;
    }
  else if (0) /* 2010_06_18 danielle no longer wants the % in the quality profile */
    {
      for (j = jmin ; j <= jmax ; j++)
	{
	  aceOutf (ao,"\t%ld", b[j]) ;
	}
      for (j = jmin ; j <= jmax ; j++)
	{
	  aceOutf (ao,"\t%.1f", (100.0*b[j])/nn) ;
	}
    }
  aceOutf (ao,"\n") ;
  aceOutDestroy (ao) ; /* will close */

  sxLetterDistribGooglePlot (sx, b) ;

  if (iMax && (sx->out == LETTERPLOT || sx->out == LETTERSHOW))
    {
      system (messprintf("\\rm %s %s.tcsh %s.ps", sx->outFileName, sx->outFileName, sx->outFileName)) ;
      ao = aceOutCreate (sx->outFileName, ".txt", FALSE, h) ;
      aceOutDate (ao, "###", sx->title) ;

      if (ao)
	{
	  for (i = di = 0 ; i < iMax ; i++)
	    {
	      if (i >= iMax1 && di == 0)
		{
		  i = sx->letterNN + 1 ;
		  di = sx->letterNN + 2 ;
		  aceOut (ao,"\n\n") ;
		  continue ;
		}
	      if (! sx->qualityPlot)
		{
		  for (n = 0, cp = "atgcn" ; *cp ; cp++)
		    {
		      j = *cp ;
		      n += bb[256*i+j] ;
		    }
		  for (cp = "atgcn" ; *cp ; cp++)
		    {
		      j = *cp ;
		      aceOutf (ao, "\t%d", bb[256*i+j]) ;
		    }
		  for (cp = "atgcn" ; *cp ; cp++)
		    {
		      j = *cp ;
		      aceOutf (ao, "\t%.1f", n > 0 ?  (100.0*bb[256*i+j])/n : 0) ;
		    }
		}
	      else
		{
		  int cumul = 0 ;
		  int cumul2[jmax] ;		  
		  float zcumul = 0 ;

		  memset (cumul2, 0, sizeof(cumul2)) ;

		  n = 0 ;
		  for (j = jmax ; j >= jmin ; j--)
		    {
		      n += bb[256*i+j] ;
		      cumul2[j] = n ;
		    }
		  for (j = jmin ; j <= jmax ; j++)
		    {
		      cumul = cumul2[j] ;
		      aceOutf (ao, "\t%d", cumul) ;
		    }
		  for (j = jmin ; j <= jmax ; j++)
		    {
		      zcumul =  n > 0 ? (100.0*cumul2[j])/n : 0 ;
		      aceOutf (ao, "\t%.1f", zcumul) ;
		    }
		}
	      aceOutf (ao, "\n");
	    }
	  ac_free (ao) ;
	  ao = aceOutCreateToFile (messprintf("%s.tcsh", sx->outFileName), "wx", h) ;
	  aceOutf(ao, "#!/bin/tcsh -f\n") ;
	  aceOutf(ao, "  gnuplot -bg white << EOF\n") ;
	  aceOutf(ao, "    set style data linespoints\n") ;
	  aceOutf(ao, "    set terminal postscript color \n") ;
	  aceOutf(ao, "    set terminal postscript solid \n") ;
	  aceOutf(ao, "    set terminal postscript landscape\n") ;
	  aceOutf(ao, "    set output '%s.ps'\n", sx->outFileName) ;
	  if (sx->title) 
	    aceOutf (ao, "set title '%s'\n", sx->title) ;

	  aceOutf(ao, "    plot  ") ;
	  /* order is important since i do not control the colors */
	  if (sx->qualityPlot)
	    {
	      for (j = jmin ; j <= jmax ; j++)
		aceOutf(ao, "%s  '%s.txt' using %d title '%d'", j>jmin ? "," : "", sx->outFileName, j - jmin+1, j) ;
	    }
	  else
	    {
	      aceOutf(ao, "  '%s.txt' using 2 title 'T'", sx->outFileName) ;
	      aceOutf(ao, ", '%s.txt' using 1 title 'A'", sx->outFileName) ;
	      aceOutf(ao, ", '%s.txt' using 3 title 'G'", sx->outFileName) ;
	      aceOutf(ao, ", '%s.txt' using 4 title 'C'", sx->outFileName) ;
	      aceOutf(ao, ", '%s.txt' using 5 title 'N'", sx->outFileName) ;
	    }
	  aceOutf(ao, "\nEOF") ;
	  aceOutDestroy (ao) ; /* will close */

	  if (sx->out == LETTERPLOT)
	    system (messprintf("(sleep 3 ; tcsh %s.tcsh) &", sx->outFileName)) ;
	  else if (sx->out == LETTERSHOW)
	    system (messprintf("(sleep 3 ; tcsh %s.tcsh ; gv %s.ps ) &", sx->outFileName, sx->outFileName)) ;
	}
    }
  sleep (3) ;
  ac_free (h) ;
  return ;
} /* sxLetterDistribReport */

/*************************************************************************************/
/*************************************************************************************/
/* pp->exitAdaptor[0] is initialized with the command line parameter ATC,GGTTTG,ttt */
static void selectPrefixCheck (SX *sx, const char *ccp)
{
  char *cp  ;
  sx->selectPrefix = strnew (ccp, sx->h) ;
  for (cp = sx->selectPrefix ; *cp ; cp++)
    *cp = dnaEncodeChar [(int)*cp] ;
  sx->selectPrefixLn = strlen (sx->selectPrefix) ;
}
/*************************************************************************************/
/* pp->exitAdaptor[0] is initialized with the command line parameter ATC,GGTTTG,ttt */
static void adaptorInit (SX *sx, unsigned const char *eV0, int isExit)
{
  static unsigned char *eV[64] ;
  static unsigned char *eW[64] ;
  static unsigned char *eZ[64] ;
  unsigned const char *ccp = 0 ;
  unsigned char cc, *cq = 0 ;
  int i, j ;
  
  switch (isExit)
    {
    case 1: memset (eV, 0, sizeof (eV)) ; break ;
    case 0: memset (eW, 0, sizeof (eW)) ; break ;
    case 2: memset (eZ, 0, sizeof (eZ)) ; break ;
    }
      
  for (i = j = 0, ccp = eV0 ; j < 64 && *ccp ; ccp++)
    {
      if (i+j == 0 || *ccp == ',')
	{
	  if (j) *cq++ = 0 ;
	  if (isExit == 1)
	    {  cq = eV[j++] = halloc(128, sx->h) ;}
	  else if (isExit == 0)
	    {  cq = eW[j++] = halloc(128, sx->h) ;}
	  else if (isExit == 2)
	    {  cq = eZ[j++] = halloc(128, sx->h) ;}
	  if (*ccp == ',') continue ;
	}
      cc = dnaEncodeChar [(int)*ccp] ;
      if (!cc &&  isExit == 1)
	messcrash ("Unrecognized char in -exitAdaptor %s", sx->exitAdaptor[0]) ;
      if (!cc && isExit == 0)
	messcrash ("Unrecognized char in -exitAdaptor %s", sx->entryAdaptor[0]) ;
      if (!cc && isExit == 2)
	messcrash ("Unrecognized char in -inPhaseEntryAdaptor %s", sx->inPhaseEntryAdaptor[0]) ;
      if (i < 127)
	*cq++ = cc ;
    }
  if (j && cq) *++cq = 0 ;
  
  if (isExit == 1) sx->exitAdaptor = &eV[0] ;
  else if (isExit == 2) sx->inPhaseEntryAdaptor = &eZ[0] ;
  else 
    {
      unsigned char *dna2, *cp, *cq ;
      sx->entryAdaptor = &eW[0] ;
      /* reverse but do not complement */
      for (i = 0 ; eW[i] ; i++)
	{
	  dna2 = eW[i] ;
	  cp = dna2 - 1 ; cq = dna2 + strlen ((char *)dna2) ;
	  while (++cp < --cq)  /* reverse but do not complement */
	    {
	      cc = *cp ;
	      *cp = *cq ;
	      *cq = cc ;
	    }
	}
    }
	  
  return ;
} /* adaptorInit */

/*************************************************************************************/
/*
  "//   -leftClipAt <int> : removes the first <int> - 1 letters of the input sequences \n"
  "//   -rightClipAt <int> : removes the letters after position <int> in the input sequences\n"
  "//   -leftClipInPhase <ATGC..> : left clip this motif if seens starting at base 1, even with 2 errors\n"
  "//   -leftClipOn <ATGC..> : left clips a given motif, for instance a ligation primer\n"
  "//      Clips the longest motif suffix matching the sequence prefix\n"
  "//      The motif is specified using the 16 letter IUPAC code ATGC RYMKSW HBVD N, upper or lower case\n"
  "//         example: left clip ATATgcctgc on GATAT exports gcctgc\n"
  "//   -rightClipOn <ATGC..> : right clips a given motif, for instance a ligation primer\n"
  "//                  Clips the longest motif prefix matching the sequence suffix\n"
  "//                  example: right clip atataTGCTGC on TGCTGC exports atata\n"
  "//   -clipN <int> : clip on N not followed by at least <int> good letters\n"
  "//   -minEntropy <int> : sequence of complexity lower than <int> are dropped\n"
*/
static BOOL clipSequence (SX *sx, int pass)
{
  unsigned char *dna = (unsigned char *) vtxtPtr (sx->dna), *dna1 = 0, *dna2 = 0 ;
  int n, n0 ;

  if (pass == 1 && ! sx->ai2)  /* do not modify this clause without modifying fRejected below after done: */
    return TRUE ;

  n0 = strlen ((const char*)dna) ;
  if (pass == 0 && sx->ai2)
    {
      char *cp = strstr ((char *)dna, "<") ;
      if (cp) 
	{
	  dna2 = (unsigned char *)strnew (cp, 0) ;
	  if (*(cp-1) == '>') cp-- ;
	  *cp = 0 ;
	}
    }
  if (pass == 1)
    {
      char *cp = strstr ((char *)dna, "<") ;
      if (cp) 
	{
	  if (*(cp-1) == '>') *(cp - 1) = 0 ;
	  *cp = 0 ;
	  dna1 = (unsigned char *)strnew ((char *)dna, 0) ;
	  n = strlen (cp+1) ;
	  memcpy (dna, cp+1, n+1) ;
	}
    }
  n = strlen ((const char*)dna) ;

  if (pass == 0 && sx->selectPrefix)
    {
      int ok ;
      unsigned char *cp, *dna3 = (unsigned char *)strnew ((char *)dna, 0) ;

      cp = dna3 - 1 ;
      while (*++cp)
	*cp = dnaEncodeChar[(int)*cp] ;
      ok = dnaPickMatch (dna3, sx->selectPrefixLn, (unsigned char *)sx->selectPrefix, 0, 0) ;
      ac_free (dna3) ;

      if (ok != 1)	
	{ n = 0 ; goto done ; }
    }

  if (sx->rightClipAt > 0 && n > sx->rightClipAt)
    {
      n = sx->rightClipAt ;
      dna[n] = 0 ;
    }
  if (sx->rightClipOn)
    {
      int iv, xv, xv2 ;
      unsigned char *cp, *dna3 = (unsigned char *)strnew ((char *)dna, 0) ;

      cp = dna3 - 1 ;
      while (*++cp)
	*cp = dnaEncodeChar[(int)*cp] ;
      for (iv = 0 ; iv < 64 && sx->exitAdaptor[iv] && *sx->exitAdaptor[iv]; iv++)
	{
	  int n0 = strlen((const char*)sx->exitAdaptor[iv]) ;
	  int n1 = n0 ;
	  xv = 0 ;
	  if (n1 > 12)
	    {
	      while ((xv2 = dnaPickMatch (dna3 + xv, n - xv, sx->exitAdaptor[iv], 2, 2)))
		xv += xv2 ;
	      while (!xv && n1 > 12) /* try to shorten the adptor and reextend the hit */
		{
		  int n2 = n1 / 2, n3 ; 
		  int e1 = 2, e2 = 2 ;
		  char cc = sx->exitAdaptor[iv][n2] ;
		  sx->exitAdaptor[iv][n2] = 0 ;
		  while ((xv2 = dnaPickMatch (dna3 + xv, n - xv, sx->exitAdaptor[iv], 2, 2)))
		    xv += xv2 ;
		  sx->exitAdaptor[iv][n2] = cc ;
		  if (xv) /* verify that the end of the read matches the vector */
		    {
		      n3 = n - xv + 1;
		      cc = sx->exitAdaptor[iv][n3] ;
		      sx->exitAdaptor[iv][n3] = 0 ;
		      if (n - xv < 12) e1 = e2 = 0 ;
		      xv2 = dnaPickMatch (dna3 + xv - 1, n - xv + 1, sx->exitAdaptor[iv], e1, e2) ;
		      sx->exitAdaptor[iv][n3] = cc ;
		      if (! xv2)
			{ xv = 0 ; n1 = 0 ; } /* reject */
		    }
		  n1 /= 2 ;
		}
	    }
	  else
	    {
	      while ((xv2 = dnaPickMatch (dna3 + xv, n - xv, sx->exitAdaptor[iv], 0, 0)))
		xv += xv2 ;
	    }
	  if (xv && n - xv >= 12)
	    {
	      sx->exitAdaptorSeqFound++ ;
	      sx->exitAdaptorTagFound += sx->count ;
	      sx->exitAdaptorBpFound += sx->count * (n - xv + 1) ;
	      sx->exitAdaptorBpClipped += (n - xv + 1) ;
	      dna[xv-1] = 0 ; n = xv-1 ;
	      break ;
	    }
	}
      ac_free (dna3) ;
    }
   if (sx->rightClipAt < 0 && n > -sx->rightClipAt)
    {
      n += sx->rightClipAt ; /* right clip a fix number of bases */
      dna[n] = 0 ;
    }

  if (sx->leftClipInPhase)
    {
      int iv, xv, xvln ;
      unsigned char *cp, *cq, *dna3 = (unsigned char *)strnew ((char *)dna, 0) ;

      cp = dna3 - 1 ;
      while (*++cp)
	*cp = dnaEncodeChar[(int)*cp] ;
      for (iv = 0 ; iv < 64 && sx->inPhaseEntryAdaptor[iv] && *sx->inPhaseEntryAdaptor[iv]; iv++)
	{
	  xv = 0 ;
	  xvln = strlen ((char *)sx->inPhaseEntryAdaptor[iv]) ;
	  if (xvln > 12)
	    {
	      xv = dnaPickMatch (dna3, n, sx->inPhaseEntryAdaptor[iv], 4, 4) ;
	    }
	  else if (xvln > 6)
	    {
	      xv = dnaPickMatch (dna3, n, sx->inPhaseEntryAdaptor[iv], 2, 2) ;
	    }
	  else
	    {
	      xv = dnaPickMatch (dna3, n, sx->inPhaseEntryAdaptor[iv], 0, 0) ;
	    }
	  if (xv == 1)
	    {
	      sx->entryAdaptorSeqFound++ ;
	      sx->entryAdaptorTagFound += sx->count ;
	      sx->entryAdaptorBpFound += sx->count * xvln ;
	      sx->entryAdaptorBpClipped += xvln ;
	      /* shift the sequence */
	      {
		cp = dna ; cq = dna + xvln ;
		while ((*cp++ = *cq++)) ;
		n -= xvln ;
	      } 
	      break ;
	    }
	}
      ac_free (dna3) ;
    }

  if (sx->leftClipOn)
    {
      int iv, xv, xv2 ;
      unsigned char cc, *cp, *cq, *dna3 = (unsigned char *)strnew ((char *)dna, 0) ;

      cp = dna3 - 1 ; cq = dna3 + strlen ((char *)dna3) ;
      while (++cp < --cq)  /* reverse but do not complement */
	{
	  cc = *cp ;
	  *cp = dnaEncodeChar[(int)*cq] ;
	  *cq = dnaEncodeChar[(int)cc] ;
	}
      for (iv = 0 ; iv < 64 && sx->entryAdaptor[iv] && *sx->entryAdaptor[iv]; iv++)
	{
	  xv = 0 ;
	  if (strlen((const char*)sx->entryAdaptor[iv]) > 12)
	    {
	      while ((xv2 = dnaPickMatch (dna3 + xv, n, sx->entryAdaptor[iv], 2, 2)))
		xv += xv2 ;
	    }
	  else
	    {
	      while ((xv2 = dnaPickMatch (dna3 + xv, n, sx->entryAdaptor[iv], 0, 0)))
		xv += xv2 ;
	    }
	  if (xv)
	    {
	      int i ;

	      sx->entryAdaptorSeqFound++ ;
	      sx->entryAdaptorTagFound += sx->count ;
	      sx->entryAdaptorBpFound += sx->count * (n - xv + 1) ;
	      sx->entryAdaptorBpClipped = (n - xv + 1) ;
	      i = n = xv - 1 ; iv = strlen ((char *)dna) - xv ;
	      for (cp = dna ; i > 0 ; cp++, i--) /* translate the dna */
		*cp = *(cp + iv + 1) ;
	      *cp = 0 ;
	      break ;
	    }
	}
      ac_free (dna3) ;
    }
  if (sx->leftClipAt > 1)
    {
      unsigned char *cp = dna, *cq = dna + sx->leftClipAt - 1  ;

      if (sx->leftClipAt < n)
	{
	  while ((*cp++ = *cq++)) ;
	  n -=  sx->leftClipAt ;
	}
      else
	n = 0 ;
      if (n <= 0)
	{ n = 0 ; goto done ; }
    }
  else if (pass == 0 && sx->UMI > 1)
    { /* clip the 5' read at the UMI, but check that it is followed by at least 12t/15 lettres */
      if (sx->UMI + 15 < n)
	{
	  int i = 15, t = 0 ;
	  unsigned char *cq = dna + sx->UMI, *cr = cq ;
	  while (i--)
	    if (*cr++ == 't')
	      t++ ;
	  if (t > 12)
	    {
	      n = sx->UMI ;
	      cq = dna + n ;
	      cr = (unsigned char *) strstr ((char *)cq, "><") ;
	      if (cr)
		while ((*cq++ = *cr++))
		  ;
	      else
		*cq = 0 ;
	    }
	}
     }
  if (sx->clipN > 0)
    {
      unsigned char *cp ;
      int i, i2, j ;

      /* start at end of sequence and clip the Ns not followed by clipN consecutive good letters */
      for (cp = dna + n - 1, i = i2 = n, j = 0 ; i >= 0 ; i--, cp--)
	{
	  if (*cp == 'n') { j = 0 ; i2 = i ? i - 1 : 0 ; } /* candidate clip position */
	  else j++ ;
	  if (i2 == 0 || j >= sx->clipN) 
	    { 
	      dna[i2] = 0 ; n = i2 ; 
	      break ;
	    }
	}

#ifdef JUNK
      /* 2012/2/25  we no longer want to clip prefixes, because of the Illumina spikes that would becode out of phase */
      /* start at first base of sequence and clip the Ns not preceded by clipN consecutive good letters */
      for (cp = dna, i = i2 = 0, j = 0 ; i < n ; i++, cp++)
	{
	  if (*cp == 'n') { j = 0 ; i2 = i + 1 ; } /* candidate clip position */
	  else j++ ;
	  if (j >= sx->clipN)
	    { /* shift the sequence */
	      if (i2 > 0)
		{
		  cp = dna ; cq = dna + i2 ;  
		  while ((*cp++ = *cq++)) ;
		  n -= i2 ; 
		} 
	      break ;
	    }
	}
#endif
    }
 done: 
  if (dna1)
    {
      int n1, n2 = n ;
      unsigned char *dna3 = 0 ;

      if (n2) dna3 = (unsigned char *)strnew ((char *)dna, 0) ;
      n1 = strlen ((char *)dna1) ; 
      if (n1) memcpy (dna, dna1, n1) ;
      if (n1 && n2) 
	{
	  memcpy (dna + n1, "><", 2) ;
	  memcpy (dna + n1 + 2, dna3, n2+1) ;
	  n = n1 + n2 + 2 ;
	}
      if (n1 && ! n2)
	{
	  memcpy (dna, dna1, n1+1) ;
	  n = n1 ; 
	}
      if (n2 && ! n1)
	{
	  dna[0] = '<' ;
	  memcpy (dna + 1, dna3, n2+1) ;
	  n = n2 + 1 ;
	}
      if (!n1 && !n2)
	{
	  dna[0] = 0 ;
	  n = 0 ;
	}
      ac_free (dna3) ;
    }
  else if (dna2)
    {
      int n2 = strlen ((char *)dna2) ;
      if (n > 0) dna[n++] = '>' ;
      memcpy (dna + n, dna2, n2 + 1) ;
      n2 += n ;
    }
  ac_free (dna1) ;
  ac_free (dna2) ;
  sx->dnaLn = n ;
  if ( 
      (pass == 1 || ! sx->ai2) &&
      (
       n < 1 || n < sx->minLength || (sx->maxLength && n > sx->maxLength) || entropyReject (sx)
       )
      )
    { 
      int dn = 0 ;
      const char *cp ;

      if (vtxtPtr (sx->dna))
	{
	  for (dn = 0, cp = vtxtPtr (sx->dna) ; *cp ; cp++)
	    switch ((int)*cp)
	      {
	      case '>':
	      case '<': dn++ ; break ;
	      }
	}
      sx->fRejected++ ;
      sx->fmRejected += sx->count ;
      sx->nRejected++ ;
      sx->tagRejected += sx->count ;
      if (dn > 1)
	{
	  sx->nRejected++ ;
	  sx->tagRejected += sx->count ;
	}
      if (dn) n0 -= dn ;
      sx->bpSeqRejected += n0 ;
      sx->bpTagsRejected += sx->count * n0 ;
      return FALSE ;
    }

  return TRUE ;
} /* clipSequence */

/*************************************************************************************/

static void dna2dnaRun (SX *sx)
{
  double nMb = 0 ;
  int nExported = 0 ;

  sx->id = vtxtHandleCreate (sx->h) ;
  sx->dna = vtxtHandleCreate (sx->h) ;
  sx->qual = vtxtHandleCreate (sx->h) ;
  
  if (sx->in != FASTQ && ! sx->runQualityFileName)
    sx->qual = vtxtHandleCreate (sx->h) ;

  if (sx->out == TAG_COUNT || sx->out == CTAG_COUNT)
    {
      sx->dict = dictHandleCreate (1000000, sx->h) ;
      sx->tagArray = arrayHandleCreate (20000000, SEQ, sx->h) ;
    }

  while (getSequence (sx))
    if (encodeSequence (sx) && clipSequence (sx, 0) && clipSequence (sx, 1))
      {
	sx->fProcessed++ ;
	if (sx->subsample && (sx->fProcessed % sx->subsample))
	  continue ;
	sx->fmProcessed += sx->count ;
	if (sx->doCount && sx->dna && sx->count)
	  {
	    int n = strlen (vtxtPtr (sx->dna)) ;
	    if (n>0 && strchr(vtxtPtr (sx->dna), '>')) n-- ;
	    if (n>0 && strchr(vtxtPtr (sx->dna), '<')) n-- ;
	    if (sx->minLn == 0 || n < sx->minLn) 
	      sx->minLn = n ;
	    if (n > sx->maxLn) 
	      sx->maxLn = n ;
	    keySet (sx->lnDistrib, n) += sx->count ;
	  }
	if (sx->getTm)
	  exportTm (sx) ;  
	else if (sx->spongeFileName || sx->shadowFileName || sx->gtfFileName || sx->refGeneFileName)
	  { /* export the relevant subsequences of the current sequence */
	    exportShadowSequences (sx, sx->shadowArray) ;
	  }
	else
	  { 
	    exportSequence (sx) ;
	    nMb += sx->dnaLn/((double)1000000.0) ;
	    nExported++ ;
	    if ( ! sx->tagArray &&
		 (
		  (sx->split && (nExported % sx->split) == 0) ||
		  (sx->splitMb && sx->splitMb < nMb)
		  )
		 )
	      { /* rotate the output file */
		char *ccp ;
		
		nMb = nExported = 0 ;

		ccp =  hprintf (sx->h, "%s.%d", sx->outFileName, ++(sx->nSplit)) ;
		if (sx->in == FASTC && (sx->out == FASTA || sx->out == FASTQ) && sx->fastc_paired_end)
		  {
		    ac_free (sx->ao1) ;
		    ac_free (sx->ao2) ;
		    sx->ao1 = aceOutCreate (hprintf(sx->h, "%s_1.", ccp), sx->fileType, sx->gzo,  sx->h) ;
		    sx->ao2 = aceOutCreate (hprintf(sx->h, "%s_2.", ccp), sx->fileType, sx->gzo,  sx->h) ;
		  }
		else
		  {
		    ac_free (sx->ao) ;
		    sx->ao = aceOutCreate (hprintf(sx->h, "%s.", ccp), sx->fileType, sx->gzo, sx->h) ;
		  }

		if (sx->aoId)
		  {
		    ac_free (sx->aoId) ; 
		    sx->aoId = aceOutCreate (ccp, ".id", sx->gzo, sx->h) ;
		  }
	      }
	  }
      }

  if (sx->tagArray && arrayMax (sx->tagArray))
    {
      if (sx->out == TAG_COUNT || sx->out == CTAG_COUNT)
	exportTagCount (sx) ;    
    }

  if (sx->letterCount)
    sxLetterDistribReport (sx) ;
  if (sx->polymer)
    sxPolymerReport (sx) ;

  return ;
} /* dna2dnaRun */

/*************************************************************************************/
/*************************************************************************************/
/* hack: inputa fasta file with repeated id/dna entries, clean it */
static void makeFastaUnique (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  DICT *dict = dictHandleCreate (100000,h) ;
  ACEIN ai = sx->ai ;
  ACEOUT ao = sx->ao ;
  int state = 0 ;
  int nIn = 0, nOut = 0 ;

  while (aceInCard (ai))
    {
      const char *cp = aceInWord (ai) ;
      if (cp)
	switch (state)
	  {
	  case 0:
	    if (cp[0] != '>')
	      continue ;
	    nIn++ ;
	    if(dictAdd (dict, cp+1, 0))
	      {
		state = 1 ;
		nOut++ ;
		aceOut (ao, cp) ;
		aceOut (ao, "\n") ;
	      }
	    break ;
	  case 1:
	    state = 0 ;
	    aceOut (ao, cp) ;
	    aceOut (ao, "\n") ;
	    break ;
	  }
    }
  ac_free (h) ;
  fprintf (stderr, "#makeFastaUnique parsed %d ids and exported %d, %.2f%%\n", nIn, nOut, 100.0*nOut/(nIn+.001)) ; 
  return ;
} /* makeFastaUnique */

/*************************************************************************************/
/*************************************************************************************/
/* quality file provides default quality values in -o FASTQ mode */

static void  sxParseQualityFile (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (sx->runQualityFileName, FALSE, h) ;
  Array aa = arrayHandleCreate (1000, int, h) ;
  int x,  q, xMax[3], part ;
  float z ;

  if (! ai)
    messcrash ("Cannot open -qualityFile %s\n", sx->runQualityFileName) ; 
  sx->qual = vtxtHandleCreate (sx->h) ;
  sx->qual2 = vtxtHandleCreate (sx->h) ;

  memset (xMax, 0, sizeof(xMax)) ;
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      if (aceInInt (ai, &part) && aceInStep(ai, '\t') && aceInInt (ai, &x) && aceInStep(ai, '\t') && aceInFloat (ai, &z) && part >= 1 && part <= 2)
	{
	  q = z ;
	  if (q < 0) q = 0 ;
	  if (q > 40) q = 40 ;
	  array (aa, part+10*x, int) = q ;
	  if (x > xMax[part]) xMax[part] = x ;
	}
    }

  for (part = 1 ; part <= 2 ; part++)
    {
      int s = 0, n = 0 ;
      for (s = n = 0, x = 10 ; x <= xMax[part] && x <20 ; x++)
	{
	  q = array (aa, part+10*x, int) ;
	  s += q ; n++ ;
	}
      if (n) s /= n ;
      for (x = 0 ;  x <= xMax[part] && x < 10 ; x++)
	{
	  q = array (aa, part+10*x, int) ;
	  if (q == 0) 
	    array (aa, part+10*x, int)  = s ;
	}
      for (s = n = 0, x = xMax[part] - 10 ; x >= xMax[part] - 20 && x > 0 ; x--)
	{
	  q = array (aa, part+10*x, int) ;
	  s += q ; n++ ;
	}
      if (n) s /= n ;
      for (  x = xMax[part] ; x >= xMax[part] - 10 && x > 0 ; x--)
	{
	  q = array (aa, part+10*x, int) ;
	  if (q == 0) 
	    array (aa, part+10*x, int)  = s ;
	}
    }
  for (part = 1 ; part <= 2 ; part++)
    {
      char buf[xMax[1] + xMax[2] + 1] ;
      for (x = 0 ; x <= xMax[part] ; x++)
	{
	  q = array (aa, part+10*x, int) ;
	  q += 3 ; /* little arbitrary bonus since we did not deduct the true SNPs from the errors */
	  buf[x] = '!' + q ;   /* convention Illumina 1.8+ , see wikipedia fastq format */
	}
      buf[x] = 0 ;
      vtxtPrintf (part == 1 ? sx->qual : sx->qual2 , "%s", buf) ;
    }

  ac_free (ai) ;
  ac_free (h) ;
} /*  sxParseQualityFile */

/*************************************************************************************/
/*************************************************************************************/

static void getSelectionDo (SX *sx, ACEIN ai, DICT *dict)
{
  while (aceInCard (ai))
    {
      register char *cp, *cq ;
      cp = aceInWord (ai) ;
      if (cp)
	{
	  if (sx->in == FASTC)
	    {
	      cq = cp + strlen(cp) - 1 ; 
	      switch ((int)(*cq))
		{
		case '>':
		case '<':
		  *cq = 0 ;
		}
	    }
	  dictAdd (dict, cp, 0) ;
	  if (strchr (cp, ' '))
	    sx->selectDictHasSpace = TRUE ;
	}
    }
} /* getSelectionDo */

/***********/

static void getSelection (SX *sx)
{
  ACEIN ai ;
  char *cp, *cp0, *cq ;

  if (sx->getList)
    {
      sx->selectDict = dictHandleCreate (1000000, sx->h) ;
      cp = cp0 = strnew (sx->getList, 0) ;
      while (1)
	{
	  cq = strstr (cp, ",") ;
	  if (cq) *cq = 0 ;
	  dictAdd (sx->selectDict, cp, 0) ;
	  if (strchr (cp, ' '))
	    sx->selectDictHasSpace = TRUE ;
	  if (cq)
	    {
	      *cq = ',' ;
	      cp = cq + 1 ;
	     }
	  else
	    break ;
	}
      ac_free (cp0) ;
    }
  else if (sx->selectFileName)
    {
      if (sx->in == RAW || sx->in == CRAW)
	messcrash ("Sorry, option -select is incompatible with -I raw/craw, since identifiers are not available") ;
      sx->selectDict = dictHandleCreate (1000000, sx->h) ;
      ai = aceInCreateFromFile (sx->selectFileName, "r", 0, sx->h) ;
      getSelectionDo (sx, ai, sx->selectDict) ;
      ac_free (ai) ;
    }
  if (sx->rejectFileName)
    {
      if (sx->in == RAW || sx->in == CRAW || sx->in == TAG_COUNT || sx->in == CTAG_COUNT)
	messcrash ("Sorry, option -reject is incompatible with -I raw/craw, since identifiers are not available") ;
      sx->rejectDict = dictHandleCreate (1000000, sx->h) ;
      ai = aceInCreateFromFile (sx->rejectFileName, "r", 0, sx->h) ;
      getSelectionDo (sx, ai, sx->rejectDict) ;
      ac_free (ai) ;
    }
} /* getSelection */

/*************************************************************************************/

static void autotest (const char *prog)
{
  AC_HANDLE h = ac_new_handle () ;
  int err = 0 ;



  if (err == 0) 
    fprintf (stderr, "Successful completion of all internal teets") ;
  ac_free (h) ;
  exit (err) ;
} /* autotest */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// dna2dna: Multilingual DNA parser\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, October 2009, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "// Convert a DNA file between several formats: raw, fasta, fastq, color fasta,\n"
	    "// optionally filtering on length or entropy, splitting in fixed size bins or by prefix\n"
	    "//\n"
	    "// Syntax:\n"
	    "// dna2dna [options]\n"
	    "// The options are all optional and may be specified in any order\n"
	    "// DNA Input\n"
	    "//   -i input_file: [default: stdin] sequence file to analyze\n"
	    "//      If the file is gzipped : -i f.fasta.gz\n, it will be gunzipped"
	    "//      Any pipe command prefixed by < is accepted:  -i \'< obj_get f.fasta.gz | gunzip\'"
	    "//  or\n"
	    "//   -i1 file1 -i2 file2: to merge paired end sequence files with matching identifiers\n"
	    "//      In this case the DNA is reexported as    atgctgt...><ttgatta... \n"
	    "//     -suffix1 <aa> -suffix2 <bb> : special pair naming convention\n"
	    "//      In case we parse 2 files -i1 -i2, we expect the smae read names in both files\n"
	    "//      In some convention, the 2 parts are given a suffix like s123a s123b\n"
	    "//      To recognize them as a pair named s123 use -suffix1 a -suffix2 b\n"
	    "//     -noNameCheck: merge file1 file2 assuming the reads are paired\n"
	    "//   -I input_format: [default: fasta] format of the input file, see formats list below\n"
	    "// DNA Output\n"
	    "//   -o output_file: [default: stdout] processed data\n"
	    "//   -O output_format: [default: fasta] format of the processed file\n"
	    "// Complement\n"
	    "//   -complement: reverse complement the sequence, incompatible with SOLiD input.\n"
	    "//   -decoy: exchange in place A<->T and G<->C, creating a non existing DNA with same statistics\n"
	    "//           if -O = ccfa,ccfar,csfasta, rotate in place the fasta letters A->T->G->C->A\n"
	    "//           in color space this keeps 0,1 and exchanges 2<->3, whereas A<->T and G<->C keeps as is 0,1,2,3\n"
	    "// Get/Select/Reject\n"
	    "//   -dnaGet <ATGC..> : only the sequences matching this prefix are exported.\n"
	    "//   -get comma separated list of sequence identifiers, accepting * and ? wild characters.\n"
	    "//   -select fileName: if specified, only the sequences listed in the first column of that file are considered.\n"
	    "//   -reject fileName: if specified, all sequences listed in the first column of that file are rejected.\n"
	    "//                     This option takes precedence over all other options.\n"
	    "// COUNT\n"
	    "//   -count : count the sequences and their multiplicities, in addition to other options\n"
	    "// Tm\n"
	    "//   -getTm : Report length, entropy and melting temperature of each sequence\n"
	    "//     The output is a 4 columns tab delimited file name, ln, S, Tm\n"
	    "//     The file is named <output_file>.TM.txt[.gz], or defaults to stdout or namesn"
	    "//     S is the the dimer Shannon entropy, equal to the length for a fuly random sequence\n"
	    "//     Tm is computed according to Maniatis\n"
	    "//\n"
	    "// gzip (recommended)\n"
	    "//   -gzi -gzo : gzi, decompress the input file; gzo, compress the output file\n"
	    "// \n"
	    "// Formats\n"
	    "//   raw : one DNA sequence per line in 16 letter IUPAC code, no identifiers\n"
	    "//   raw<n> :  Input only, the DNA is in column Number (example: raw1 == raw, raw2, raw3 ...)\n"
	    "//   raw<n1>yn<n2> :  In addition flag in column n2 must be Y\n"
	    "//      example: raw9yn22  imports DNA from column 9 if column 22 says Y\n"
	    "//      in raw or raw1 mode, the identifier of the sequence is assumed to be in column 2\n"
	    "//      in raw<n> mode, n>1, the identifier of the sequence is assumed to be in column 1\n"
	    "//   tc : Tag Count : DNA sequence in 16 letter IUPAC code 'tab' multiplicity\n"
	    "//   fasta : >identifier 'new line' DNA in 16 letter IUPAC code\n"
	    "//     As input, the DNA can be spread over any number of lines\n"
	    "//     As output, the DNA is exported on a single line except if -maxLineLn <number> is specified \n"
	    "//     If the input does not provide identifiers, e.g. in the raw or tc input formats, the option\n"
	    "//          -prefix <txt> serves to construct the identifiers as 'txt.1 txt.2 txt.3 ...\n"
	    "//   fastc : >identifier#tag_count 'new line' DNA in 16 letter IUPAC code\n"
	    "//     Same as fasta, but the identifier includes the tag count multiplicity after the # symbol\n"
	    "//   fastq : @identifier 'new line'  DNA 'new line' +same_identifier 'new line' quality\n"
	    "//     Qualities are specific of the fastq format and dropped when exporting in any other format\n"
	    "//     Hence -O fastq is only possible if -I fastq, for example to split or filter the input file\n"
	    "//     or if we provide \n"
	    "//   tensorflow -rightClipAt n : export identifier, coma, [[1,0,0,0],...[0,1,0,0]] n bases as one hot python arrays\n"
            "//   -runQuality <qualityFile> : 3 columns tab delimited: r x q\n"
	    "//         r=1 for read 1, r=2 for read 2, exported as $out_[12].fastq[.gz]\n"
	    "//         x position integer=1...n where n is at least as long as the longest read\n"
	    "//         q quality interger in range 0 to 50 (-10 log_10(sequencing error rate))\n"
	    "//   -runTitle <text> : The runTitle is added as a prefix in each fastq seqeunce identifier\n"
	    "//\n"
	    "// SOLiD formats and options\n"
	    "//     Recall that a SOLiD sequence starts with an Anchor letter followed by transitions, a natural\n"
	    "//     code corresponding to their ligation sequencing protocol\n"
	    "//   csfasta : LifeTech SOLiD color fasta format: anchor (usually T) followed by transitions 0123\n"
	    "//   csfastc : idem, but with a #multiplicity factor in the sequence identifiers\n"
	    "//   ccfa : modified LifeTech SOLiD format, with 0123 replaced by acgt\n"
	    "//   ccfar : modified LifeTech SOLiD format with 0123 replaced by tgca\n"
	    "//     The advantage of the ccfa format is that it is accepted by standard DNA programs such as\n"
	    "//     Blast/Blat or Bowtie. For example one may align a cfa file onto the direct strand of the genome\n"
	    "//     by transforming both sequences in ccfa format: this is much preferable to decoding cfa to\n"
	    "//     fasta then aligning to the fasta genome, because in such a naive scheme, any missmatch\n"
	    "//     irreversibly disrupts the homology between the genome and the SOLiD sequence\n"
	    "//     downstream of the mismatch\n"
	    "//     The ccfar format is needed to align sequences to the complementary strand. A possible\n"
	    "//     protocol to align a SOLiD file on the genome would be to reformat the file in ccfa,\n"
	    "//     to reformat the genome in both ccfa and ccfar and to run the aligner twice.\n"
	    "//     Preferably, one may align just once using the aceview 'align' program or any other\n"
	    "//     aligner which understands the native SOLiD cfa format and implements the SOLiD error\n"
	    "//     correction mechanism.\n"  
	    "//   -n2a : 'n', 'N' or '.' characters in the input sequence file are converted to 'a', and become\n"
	    "//     exportable in ccfa format. This option should be used for converting a genome file, which\n"
	    "//     usually contains stretches of N, to SOLiD ccfa and ccfar color formats.\n"
	    "//   -jumpAnchor : removes the SOLiD Anchor letter. The resulting color fasta file is\n"
	    "//     no longer decodable into fasta format, but now matches exactly the SOLiD color fasta genome.\n"
	    "//     If the anchor letter was not removed in this way, it would not match the corresponding position \n"
	    "//     on the ccfa genome, resulting in a spurious mismatch in naive aligner programs, such as \n"
	    "//     BLAST/BLAT/Bowtie which, as of December 2009, do not process the color format. \n"
	    "//\n"
	    "// Filters\n"
	    "//   -subsample <int> : Export only every n-th read\n"
	    "//   -selectPrefix <ATGC..> : Consider only the sequences starting by these letters\n" 
	    "//   -leftClipAt <int> : removes the first <int> - 1 letters of the input sequences \n"
	    "//   -rightClipAt <int> : removes the letters after position <int> in the input sequences\n"
	    "//                 a negative value, say -4 , rigth clips the last 4 bases\n"
	    "//   -leftClipOn <ATGC..> : left clips a given motif, for instance a ligation primer\n"
	    "//      The motif is specified using the 16 letter IUPAC code ATGC RYMKSW HBVD N, upper or lower case\n"
	    "//                  Clip on the most 5' position matching the entire motif\n"
	    "//                  exactly for short motifs or with at most 2 errors if > 16bp\n"
	    "//                  or clip the longest motif suffix matching the beginning of the sequence\n"
	    "//         example: left clip ATATgcctgc on GATAT exports gcctgc\n"
	    "//   -rightClipOn <ATGC..> : right clips a given motif, for instance a ligation primer\n"
	    "//                  Clip on the most 3' position matching the entire motif\n"
	    "//                  exactly for short motifs or with at most 2 errors if > 16bp\n"
	    "//                  or clip the longest motif prefix matching the end of the sequence\n"
	    "//         example: right clip atataTGCTGC on TGCTGC exports atata\n"
	    "//   -clipN <int> : clip on terminal N not followed by at least <int> good letters\n"
	    "//   -minQuality <int> : default 30, using @==0, in -I fastq format, letter with quality <= 30  become N\n"
	    "//   -minLength <int> : sequence of clipped length lower than <int> are dropped\n"
	    "//   -minMultiplicity <int> : sequence with multiplicity lower than <int> are dropped\n"
	    "//   -minEntropy <int> : sequence of complexity lower than <int> are dropped\n"
	    "//     Our simple entropy formula only considers the base composition. It is measured in bp \n"
	    "//     equivalent and defined as entropy = - sum over x=A,T,G,C {n(x) log4(n(x)/sequence length)}\n"
	    "//     this means that AAAAAA has entropy zero, but the entropy of a complex sequence\n"
	    "//     of length n is just below n. Letters other than ATGC do not contribute to the calculation.\n"
	    "//   -fastqSelect <filter> : strand/chastitize CASAVA-8 fastq format\n"
	    "//     Keep only the sequences where the second word on the identifier field matches *filter*\n"
	    "//     Example 1:N: 2:N: (1/2 are the pair members, Y/N means reject/keep\n"
	    "//     Notice that in RNA-seq UDG protocol, 1 is Reverse, 2 is Forward\n"
	    "//   -shadow shadow_file_name : deduce the gene DNA from the shadow coordinates and the input chromosome sequence\n"
	    "//     Expects a 6 column tab delimited file\n"
	    "//          gene x1 x2 chromosome a1 a2\n"
	    "//     all coordinates are 1-based, chromosome should be present in the input file (-i input -I format)\n"
	    "//     a1<->x1, a2<->x2 correspond to each other, \n"
	    "//     if (a2-a1)*(x2-x1)<0, the gene sequence is read on the minus strand of the chromosome\n"
	    "//     To export a spliced gene, give one line per exon\n"
	    "//     The shadow file is sorted internally, so the order of the lines is immaterial\n"
	    "//   -sponge sponge_file_name : deduce the gene DNA from the sponge coordinates and the input chromosome sequence\n"
	    "//     Very similar to -shadow option but expects a 6 column tab delimited file\n"
	    "//          transcript exon_number chromosome a1 a2 gene\n"
	    "//   -gtf gtf_file_name -gtfGenome g.fasta : deduce the DNA from gtf coordinates and the genome fasta sequence\n"
	    "//   -gff3 gff3_file_name: idem, but assume the similar gff3 format\n"
	    "//       -gffTAIR : perfectly clean gff3 as provided by TAIR\n"
	    "//      In both cases you can specify\n"
	    "//       -gtfRemap <text> : export a 7 column WIGGLEREMAP file:  text, 6 col shadow format\n"
	    "//                        : export also a predicted gene set of mode in .transcript.ace\n"
	    "//       -gtfGenome <genome.fasta> : extract from this genome the transcript.fasta file\n"
	    "//   -refGene refGene_file_name: idem, but using the UCSC refGene format\n"
	    "//   -keepName : in raw or tc formats export the sequence name(s) in column 2\n"
	    "//               in fastc export a separate file of names for all the multiple sequences\n"
	    "//   -keepName1 : append /1 to the original identifier\n"
	    "//   -keepName2 : append /2 to the original identifier\n"
	    "//   -keepNameBam : append /1 or /2 to the original identifier if there is a 1 or a 2 in the flag in column 2\n"
	    "// Actions\n"
	    "//   -split <int> : split the output files in chunks of specified number of sequences. We recommend 5000000. \n"
	    "//     The multiple output files are named <output_file>.1 .2 .3 ... .n\n"
	    "//   -splitMb <int> : split the output files in chunks of specified number of Mega bases.  We recommend 250. \n"
	    "//     The multiple output files are named <output_file>.1 .2 .3 ... .n\n"
	    "//   -splitByPrefix n : n=1,2,3 or 4, split the output in 4^n files according to the first n letters\n"
	    "//     Note that to merge input file no option is required, pipe them all to stdin, UNIX way\n"
	    "//   -splitPairs : split the output files, separating _1 and _2 reads\n"
	    "//     This option is only used when transforming file format form fastc to fasta : -I fastc -O fasta\n"
	    "//   -appendByPrefix n : idem, but appending to, rather than clearing a pre-existing output file\n"
	    "//     This is useful to merge and fan out by prefix a set of input files\n"
	    /*
	    "//   -splitLongDna n: NOT IMPLEMENTED split long sequences in pieces\n"
	    "//     This options splits a long sequence in pieces at most n bases long\n"
	    "//     setting maxLineLn splits the DNA over several lines of chosen maximal length\n"
	    */
	    "//   -upper export DNA letters in upper casse (default is lower case)\n"
	    "//   -maxLineLn <int>: [default 10000] limit the length of the DNA lanes\n"
	    "//     By default, the DNA is printed on a single, possibly very long, line\n"
	    "//     setting maxLineLn splits, in fast[ac] format the DNA over several lines of chosen maximal length\n"
	    "//     If the value is <=0, the line is never split\n"
	    "//     Notice that if you export a chromosome on a single fasta line, although this is\n"
	    "//     semantically correct, most programs will crash (buffer size exceeded) when reading it\n"
	    "//     this is why the default is chosen large, but not infinite\n"
	    "//   -makeFastaUnique: hack to analyze Ultimal BAM files\n"
	    "//     BamDeal bam2fasta takes a bam and exorts identifier + sequence in fasta format\n"
	    "//     but the same identifier and DNA maybe exported several times\n"
	    "//     This option exports a clean fasta file, with each identifier presnt only once\n"
	    "//     We used this hack in septembre 2023 to clean our first NASA Ultima runs\n"
	    "//\n"
	    "// Letter/quality profiles\n"
	    "//   [-qualityPlot] : plot the profile of the fastq qualities rather than the letter profile\n"
	    "//   -title <text> : The title is copied at the top of the output file\n"
	    "//   -plot <int> : Plots the piled distribution of the first n letters, useful to show biases in the sequences\n"
	    "//   -gnuplot <int> : requires the option -o file_name and the program \"gnuplot\"\n"
	    "//          idem but exports the distribution in fil_name.txt and a picture in file_name.ps\n"
	    "//   -gvplot <int> : requires -o file_name, and the programs \"gnuplot\" and \"gv\"\n"
	    "//          idem but directly displays the postscript file file_name.ps using gv\n"
	    "//   -polymer <int> : histogram of homopolymers of length 1 to n\n"
	    "// Special applications\n"
	    "//   -editGenome : edit the [-i genome] with subs and non-sliding indels\n"
	    "//          -ctsPeriod p : periodic editions,  every p positions\n"
	    "//          -ctsRandom r : (integer) error rate per thousand\n"
	    "//          -ctsSubRate : default 100, substitutions:indels (60 means 60:40)\n"
	    "//   -createTestSet N -ctsL ln [-ctsStep s] [-ctsShift da] -ctsPeriod p] : Create a perfect heterozyguous test set\n"
	    "//                    [-ctsP pair_length] [-ctsIndel] [-ctsAddPolyA]  [-ctsReverse] [-O fastq | raw] \n"
	    "//     The purpose is to create a benchmark to assess the sensitivity/specifity of an aligner/SNV caller program\n" 
	    "//     Parse the template sequence (stdin or specified using -i mysequence -I input_format)\n"
	    "//     Typically, the template would be a 1Mb region of the genome, or the sequences of a set of genes\n"
	    "//     The program creates 8 sets of N tiling sequences (or pairs), say N=1000000, of length ln, say ln=50\n"
	    "//     starting every s base, by default s=1, starting at da, containing known errors using -cstPeriod -ctsRandom -ctsSubrate\n"
           "//\n"
            "//   dna2dna -i f1.fastq -I fastq -fastqSelect 1:N: -O fasta\n"
            "//     Exports as fasta the chastitized reverse reads in a UDG CASAVA-8 fastq file\n"
 	    "//     if fln > 0, the fragment or pair length, paired end sequences are created in fastc format\n" 
	    "//     The \'exact.forward set\' contains exact subsequence of the input, in fasta format\n"
	    "//     The \'variant.forward set\' is derived by modyfing the input every p base.\n"
	    "//     So each short sequence contains one variation\n"
	    "//     The default is to create substitutions, or alternated inserts and deletes if -ctsIndel is specified\n"
	    "//     The two sets are then complemented and exported as 'exact/variant.reverse sets'\n"
	    "//     The four sets are then reexported in solid cs-fasta format\n"
	    "//     Finally, if -O fastq is specified, a quality factor is added to all sequences\n"
	    "//     The -o option is require, it specifies a common prefix for all output files, e.g. -o mydir/test1\n"
	    "//   -autotest : run autotest, in case an error is reported, please email mieg@ncbi.nlm.nih.gov\n"
	    "// Examples:\n"
	    "//   dna2dna -i f1 -I fasta -O raw -minEntropy 16 | sort -u | wc\n"
	    "//     Counts in a sequence file f1 of any format, here fasta, the number of distinct sequences which are too short\n"
	    "//     or too simple, with entropy below the equivalent of 16 well mixed bases.\n"
	    "//\n"
	    "//   dna2dna -i f1.txt  -I raw -O fasta -prefix ZXW > f1.fasta\n"
	    "//     Export the raw sequences present in f1.txt, in fasta format, and names them ZXW.1 ZXW.2 ...\n"
	    "//\n"
	    "//   dna2dna -i f1.fastq -I fastq -O tag_count -o f2.txt\n"
	    "//     Exports in the simple text file f2.txt, without identifiers or quality scores, the sequences contained in the\n"
	    "//     fastq file f1.fastq, and their multiplicity.\n"
	    "//     This is the most compact human readable format\n"
	    "//\n"
	    "//   dna2dna -i f.fasta -rightClipAt 600 -maxLineLn 60 -o g\n"
	    "//     Exports the first 600 bases of each sequence of f.fasta as 10 lines of length 60 in file g.fasta\n"
	    "// Caveat:\n"
	    "//   Lines starting with '#' are considered as comments and dropped out\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  dna2dna --help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/
/*************************************************************************************/

static BOOL checkFormat (const char *io, DNAFORMAT *ip, const char *ccp, char *ftype)
{
  int i ;
  const char **f ; 
  const char *ff[] = { "fasta", "fastc", "fastq", "csfasta",  "csfastc", "ccfa", "ccfar", "raw", "craw", "tc", "ctc", "count", "tensorflow", 0 } ; 

  for (i = 0 , f = ff ; *f ; i++, f++)
    if (! strcasecmp (*f, ccp))
      { 
	*ip = i ; 
	if (*io == 'O')
	  strcpy (ftype, ff[i]) ;
	return TRUE ;
      }

  usage (messprintf ("Unknown format in option -%s %s", io, ccp)) ;

  return FALSE ;
} /* checkFormat */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  const char *ccp = 0 ;
  SX sx ;
  AC_HANDLE h = 0 ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&sx, 0, sizeof (SX)) ;
  sx.h = h ;

  /* freeOut (0) ;  needed be the linker 2017_02_26 */
  /* optional arguments */

  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  getCmdLineOption (&argc, argv, "-i", &sx.inFileName) ;
  getCmdLineOption (&argc, argv, "-i1", &sx.inFileName1) ;
  getCmdLineOption (&argc, argv, "-i2", &sx.inFileName2) ;
  getCmdLineOption (&argc, argv, "-suffix1", &sx.suffix1) ;
  getCmdLineOption (&argc, argv, "-suffix2", &sx.suffix2) ;
  
  getCmdLineOption (&argc, argv, "-o", &sx.outFileName) ;

  if ( getCmdLineOption (&argc, argv, "-autotest", 0))
    autotest (argv[0]) ;

  sx.n2a = getCmdLineBool (&argc, argv, "-n2a") ;
  getCmdLineInt (&argc, argv, "-xor", &(sx.xor)) ;
  sx.jumpAnchor = getCmdLineBool (&argc, argv, "-jumpAnchor") ;
  getCmdLineInt (&argc, argv, "-subsample", &sx.subsample) ;
  sx.fastc_paired_end = getCmdLineBool (&argc, argv, "-splitPairs") ;

  if (getCmdLineInt (&argc, argv, "-appendByPrefix", &sx.splitByPrefix))
    sx.appendByPrefix = TRUE ;
  sx.complement = getCmdLineOption (&argc, argv, "-complement", 0) ;
  sx.upper = getCmdLineOption (&argc, argv, "-upper", 0) ;
  sx.decoy = getCmdLineOption (&argc, argv, "-decoy", 0) ;
  sx.keepName = getCmdLineOption (&argc, argv, "-keepName", 0) ;
  sx.keepNameBam = getCmdLineOption (&argc, argv, "-keepNameBam", 0) ;
  sx.keepName1 = getCmdLineOption (&argc, argv, "-keepName1", 0) ;
  sx.keepName2 = getCmdLineOption (&argc, argv, "-keepName2", 0) ;
  
  sx.gzi = getCmdLineOption (&argc, argv, "-gzi", 0);
  sx.gzo = getCmdLineOption (&argc, argv, "-gzo", 0);

  getCmdLineOption (&argc, argv, "-title", &sx.title) ;
  getCmdLineOption (&argc, argv, "-runTitle", &sx.runTitle) ;
  if (sx.runTitle && strlen (sx.runTitle) > 1000)
    messcrash ("-runTitle parameter is too long (ln=%d > 1000)\n  -runTitle %s\n",  strlen (sx.runTitle), sx.runTitle) ;
  getCmdLineInt (&argc, argv, "-minEntropy", &sx.minEntropy) ;
  getCmdLineInt (&argc, argv, "-minMultiplicity", &sx.minMultiplicity) ;
  getCmdLineInt (&argc, argv, "-minLength", &sx.minLength) ;
  getCmdLineInt (&argc, argv, "-maxLength", &sx.maxLength) ;
  sx.minQuality = 0 ; /* do not look at the quality (beware it may start at '!'==33 or at '@'=64 */
  getCmdLineInt (&argc, argv, "-minQuality", &sx.minQuality) ;
  getCmdLineInt (&argc, argv, "-split", &sx.split) ;
  getCmdLineInt (&argc, argv, "-splitMb", &sx.splitMb) ;
  getCmdLineInt (&argc, argv, "-splitByPrefix", &sx.splitByPrefix) ;

  sx.makeTestGenome = getCmdLineBool (&argc, argv, "-editGenome") ;
  sx.noNameCheck = getCmdLineBool (&argc, argv, "-noNameCheck") ;

  getCmdLineInt (&argc, argv, "-createTestSet", &sx.makeTest) ;
  if (sx.makeTest < 0)
    messcrash ("-createTestSet %d < 0, expecting the number of test sequence to be created, try -help", sx.makeTest) ;

  sx.makeTestLength = 50 ; sx.makeTestStep = 1 ; /* default */

  getCmdLineInt (&argc, argv, "-ctsP", &sx.makeTestPairLength) ;
  if (sx.makeTestPairLength < 0)
    messcrash ("-ctsP %d < 1, expecting the fragment length (pair length) of the test sequences to be created, try -help", sx.makeTestLength) ;
  getCmdLineInt (&argc, argv, "-ctsL", &sx.makeTestLength) ;
  if (sx.makeTestLength < 0)
    messcrash ("-ctsL %d < 1, expecting the length of the test sequences to be created, try -help", sx.makeTestLength) ;
  if (sx.makeTestPairLength > 0 && sx.makeTestPairLength  < sx.makeTestLength)
    sx.makeTestLength = sx.makeTestPairLength ;
  sx.makeTestPeriod = 0 ;

  sx.makeTestReverse =  getCmdLineBool (&argc, argv, "-ctsReverse") ;
  getCmdLineInt (&argc, argv, "-ctsPeriod", &sx.makeTestPeriod) ;
  sx.makeTestSubRate = 100 ;
  getCmdLineInt (&argc, argv, "-ctsSubRate", &sx.makeTestSubRate) ;
  if (getCmdLineInt (&argc, argv, "-ctsRandom", &sx.makeTestPeriod) &&
      sx.makeTestPeriod > 0)
    { 
      sx.makeTestRandom = TRUE ; 
      sx.makeTestPeriod = 1000/sx.makeTestPeriod ; 
    }
  if (sx.makeTestPeriod <= 0)
    sx.makeTestPeriod = sx.makeTestLength ; /* default */
  getCmdLineInt (&argc, argv, "-ctsShift", &sx.makeTestShift) ;
  getCmdLineInt (&argc, argv, "-ctsStep", &sx.makeTestStep) ;
  if (sx.makeTestStep < 0)
    messcrash ("-ctsStep %d < 1, expecting the step between the test sequences to be created, try -help", sx.makeTestStep) ;
  sx.makeTestInDel = getCmdLineOption (&argc, argv, "-ctsInDel", 0) ;
  sx.makeTestAddPolyA = getCmdLineOption (&argc, argv, "-ctsAddPolyA", 0) ;

  getCmdLineInt (&argc, argv, "-UMI", &sx.UMI) ;
  if (sx.UMI < 0)
    messcrash ("-UMI %d < 1, expecting a clipping length if an N is too close to the edge of the read, try -help", sx.clipN) ;

  getCmdLineInt (&argc, argv, "-clipN", &sx.clipN) ;
  if (sx.clipN < 0)
    messcrash ("-clipN %d < 1, expecting a clipping length if an N is too close to the edge of the read, try -help", sx.clipN) ;
  getCmdLineInt (&argc, argv, "-leftClipAt", &sx.leftClipAt) ;
  if (sx.leftClipAt < 0)
    messcrash ("-leftClipAt %d < 1, expecting a positive left clip position, try -help", sx.leftClipAt) ;
  getCmdLineInt (&argc, argv, "-rightClipAt", &sx.rightClipAt) ;
  if (0 && sx.rightClipAt < 0)
    messcrash ("-rightClipAt %d < 1, expecting a positive right clip position, try -help", sx.rightClipAt) ;
  getCmdLineOption (&argc, argv, "-gtfRemap", &sx.gtfRemapPrefix) ;
  getCmdLineOption (&argc, argv, "-gffRemap", &sx.gtfRemapPrefix) ;

  if (getCmdLineOption (&argc, argv, "-rightClipOn", &sx.rightClipOn))
    adaptorInit (&sx, (unsigned const char*) sx.rightClipOn, TRUE) ;
  if (getCmdLineOption (&argc, argv, "-leftClipOn", &sx.leftClipOn))
    adaptorInit (&sx, (unsigned const char*) sx.leftClipOn, FALSE) ;
  if (getCmdLineOption (&argc, argv, "-leftClipInPhase", &sx.leftClipInPhase))
    adaptorInit (&sx, (unsigned const char*) sx.leftClipInPhase, 2) ;
  if (getCmdLineOption (&argc, argv, "-selectPrefix", &ccp))
    selectPrefixCheck (&sx, ccp) ;
  sx.maxLineLn = 0 ;
  getCmdLineInt (&argc, argv, "-maxLineLn", &sx.maxLineLn) ;
  /*   getCmdLineInt (&argc, argv, "-splitLongDna", &sx.splitLongDna) ; */

  getCmdLineOption (&argc, argv, "-get", &sx.getList) ;
  getCmdLineOption (&argc, argv, "-select", &sx.selectFileName) ;
  getCmdLineOption (&argc, argv, "-reject", &sx.rejectFileName) ;
  getCmdLineOption (&argc, argv, "-shadow", &sx.shadowFileName) ;
  getCmdLineOption (&argc, argv, "-sponge", &sx.spongeFileName) ;
  getCmdLineOption (&argc, argv, "-runQuality", &sx.runQualityFileName) ;
  getCmdLineOption (&argc, argv, "-refGene", &sx.refGeneFileName) ;
  getCmdLineOption (&argc, argv, "-gtf", &sx.gtfFileName) ;
  if (getCmdLineOption (&argc, argv, "-gff", &sx.gtfFileName))
     sx.isGff3 = TRUE ;
  if (getCmdLineOption (&argc, argv, "-gff3", &sx.gtfFileName))
    {
      sx.isGff3 = TRUE ;
      if (sx.gtfFileName && strstr (sx.gtfFileName, "Bacteria"))
	sx.gffBacteria = TRUE ;
    }
  sx.gffTAIR = getCmdLineBool (&argc, argv, "-gffTAIR") ;
  getCmdLineOption (&argc, argv, "-gtfGenome", &sx.gtfGenomeFastaFileName) ;
  getCmdLineOption (&argc, argv, "-gffgenome", &sx.gtfGenomeFastaFileName) ;
  sx.makeFastaUnique = getCmdLineBool (&argc, argv, "-makeFastaUnique") ;

  if (getCmdLineInt (&argc, argv, "-polymer", &sx.polymer))
    sx.out = POLYMER ;

  getCmdLineOption (&argc, argv, "-fastqSelect", &sx.fastqSelect) ;
  sx.letterNN = -9999999 ;
  if (getCmdLineInt (&argc, argv, "-plot", &sx.letterNN))
    sx.out = LETTER ;
  sx.qualityPlot = getCmdLineOption (&argc, argv, "-qualityPlot", 0) ;
  if (getCmdLineInt (&argc, argv, "-gnuplot", &sx.letterNN))
    sx.out = LETTERPLOT ;
  if (getCmdLineInt (&argc, argv, "-gvplot", &sx.letterNN))
    sx.out = LETTERSHOW ;
  if (sx.letterNN > 0)
    { 
      if (sx.letterNN < 30) sx.letterNN = 30 ; /* always analyse at least 30 bases, this makes the report format more regular */ 
      sx.letterLength = halloc ((2 * sx.letterNN + 5) * sizeof (long int), h) ;
      sx.letterCount = halloc (256 * sizeof (long int), h) ;
      sx.letterProfile = halloc (256 * (2*sx.letterNN + 5) * sizeof (long int), h) ;
    }
  if (sx.letterNN > -9999999 && sx.letterNN <= 0)
    usage ("The length of the requested letter profile should be a positive number") ;
  if (sx.polymer < 0 || sx.polymer > POLYMER_MAX - 10)
    usage ("The maximal polymer length should be between 1 and 1000") ;

  sx.prefix = "n" ; /* default */
  getCmdLineOption (&argc, argv, "-prefix", &sx.prefix) ;

  strcpy (sx.fileType, "fasta") ;
  if (getCmdLineOption (&argc, argv, "-I", &ccp))
    {
      int col = 0 ;   /* column containing the DNA */
      int yncol = 0 ; /* Y or N flag column, reject the N */
      char cc ;
      /* recognize raw9 as raw format, split column 9 */
      if (sscanf (ccp, "raw%d%c", &col, &cc) == 1 && col < 100 && col > 0)
	{ ccp = "raw" ; sx.rawColumn = col ; }
      /* recognize raw9 as raw format, split column 9 */
      if (sscanf (ccp, "raw%dyn%d%c", &col, &yncol,&cc) == 2 && col < 100 && col > 0)
	{ ccp = "raw" ; sx.rawColumn = col ; sx.ynColumn = yncol ; }
      checkFormat ("I", &sx.in, ccp, sx.fileType) ;
    }
  sx.doCount = getCmdLineOption (&argc, argv, "-count", 0) ;
  sx.getTm = 
    getCmdLineOption (&argc, argv, "-getTm", 0) ||
    getCmdLineOption (&argc, argv, "-getTM", 0)
    ;
  if (getCmdLineOption (&argc, argv, "-O", &ccp))
    checkFormat ("O", &sx.out, ccp, sx.fileType) ;
  else if (sx.doCount)
    sx.out = COUNT ;
  if (sx.out == COUNT)
    sx.doCount = TRUE ;

  /*
    if (sx.splitLongDna && sx.out != RAW && sx.out != FASTA && sx.out != FASTC)
    usage ("splitLongDna only works for out format raw, fasta, fastc, sorry") ;
  */

  if (sx.complement)
    {
      switch (sx.in)
	{
	case CSFASTA:
	case CSFASTC:
	case CCFA:
	case CCFAR:
	  messcrash ("-complement is incompatible with -I CSFASTA/CCFA/CCFAR formats, sorry") ;
	  break ;
	default:
	  break ;
	}
      switch (sx.out)
	{
	case CSFASTA:
	case CSFASTC:
	case CCFA:
	case CCFAR:
	  messcrash ("-complement is incompatible with -O CSFASTA/CCFA/CCFAR formats, sorry") ;
	  break ;
	default:
	  break ;
	}
    }

  if (sx.splitByPrefix  && (sx.splitByPrefix < 1 || sx.splitByPrefix > 4))
    messcrash ("Bad parameter -splitByPrefix %d, should be 1,2,3 or 4", sx.splitByPrefix) ;
  if (sx.appendByPrefix && sx.gzo)
    messcrash ("Sorry, appendByPrefix && gzo are incompatible, this would creat too many pipes on the computer") ;
      if (sx.splitByPrefix && sx.gzo)
	messcrash ("Sorry, splitByPrefix && gzo are incompatible, this would creat too many pipes on the computer") ;
  
  if (sx.appendByPrefix && ! sx.outFileName)
    messcrash ("With option -appendByPrefix, option -o output_file_name is mandatory, sorry") ;
  if (sx.splitByPrefix && ! sx.outFileName)
    messcrash ("With option -splitByPrefix, option -o output_file_name is mandatory, sorry") ;
  if (sx.in == FASTC && (sx.out == FASTA  || sx.out == FASTQ) && ! sx.outFileName)
    messcrash ("-I FASTC -O FASTA/FASTQ conversion requires option -o output_file_name, sorry") ;
  if (sx.in == FASTC && (sx.out == FASTA  || sx.out == FASTQ) && sx.appendByPrefix)
    messcrash ("-I FASTC -O FASTA/FASTQ conversion is incompatible with -appendByPrefix, sorry") ;
  if (sx.in == FASTC && (sx.out == FASTA  || sx.out == FASTQ) && sx.maxLineLn)
    messcrash ("-I FASTC -O FASTA/FASTQ conversion is incompatible with -maxLineLn, sorry") ;
  if (sx.in != FASTQ && sx.out == FASTQ && ! sx.runQualityFileName)
    messcrash ("-O FASTQ requires -I FASTQ or -runQuality, to provide the quality values") ;
  if (sx.out == FASTQ &&  sx.maxLineLn)
    messcrash ("-O FASTQ export format is incompatible with -maxLineLn, sorry") ;
  if (sx.out == TENSORFLOW &&  ! sx.rightClipAt)
    messcrash ("-O TENSORFLOW export fromat requires -rightClipAt, sorry") ;
  if (sx.in == TENSORFLOW)
    messcrash ("-I TENSORFLOW is not an importation format, sorry") ;
  
  fprintf (stderr, "// start: %s\n", timeShowNow()) ;

  if (sx.runQualityFileName)
    sxParseQualityFile (&sx) ;
  
  if (!sx.splitByPrefix && sx.letterNN <= 0 && sx.polymer <= 0 && ! sx.getTm  && ! sx.makeTest)
    {
      if (sx.outFileName)
	{
	   if (sx.split || sx.splitMb)
	     ccp = hprintf(h, "%s.1",sx.outFileName) ;
	   else
	     ccp = sx.outFileName ;
	   if (sx.in == FASTC && (sx.out == FASTA || sx.out == FASTQ) && sxChecInputHasPairs (&sx))
	     {
	       sx.ao1 = aceOutCreate (hprintf(h, "%s_1.", ccp), sx.fileType, sx.gzo, h) ;
	       sx.ao2 = aceOutCreate (hprintf(h, "%s_2.", ccp), sx.fileType, sx.gzo, h) ;
	     }
	   else
	     sx.ao = aceOutCreate (hprintf(h, "%s.", ccp), sx.fileType, sx.gzo, h) ;
	   
	   if (! sx.ao && (! sx.ao1 || ! sx.ao2))
	     exit (1) ;
	}
      else
	{
	  if (sx.gzo)
	    sx.ao = aceOutCreateToPipe ("gzip ", h) ;
	  else
	    sx.ao = aceOutCreateToStdout (h) ;
	  sx.outFileName = "stdout" ;
	}
    }
  else if (sx.splitByPrefix && sx.letterNN <= 0 && sx.polymer <= 0)
    {
      if (sx.outFileName)
	{
	  unsigned int i, ii, j, n = 1 << (2 * sx.splitByPrefix) ;
	  char buf[sx.splitByPrefix + 1], *xx = "ACGT" ;

	  for (ii = 0 ; ii < n ; ii++)
	    {
	      i = ii ;
	      if (ii >= 256) messcrash ("please raise the declaration of aos[256], we need 2^%d", sx.splitByPrefix) ;
	      for (j = 0, i = ii ; j <  sx.splitByPrefix ; j++, i >>= 2)
		buf[sx.splitByPrefix - j - 1] = xx[i & 0x3] ;
	      buf[sx.splitByPrefix] = 0 ;

	      if (sx.in == FASTC && sx.out == FASTA)
		{ 
		  ccp = hprintf (h, "%s.%s", sx.outFileName, buf) ;
		  sx.aos1[ii] = aceOutCreate (hprintf(h, "%s_1.", ccp), sx.fileType, sx.gzo, h) ;
		  sx.aos2[ii] = aceOutCreate (hprintf(h, "%s_2.", ccp), sx.fileType, sx.gzo, h) ;
		  if (! sx.aos1[ii] || !sx.aos2[ii]) 
		    messcrash ("Cannot create split file : %s", ccp) ;
		}
	      else
		{ 
		  if (sx.appendByPrefix)
		    sx.aos[ii] = aceOutCreateToFile (hprintf (h, "%s.%s.%s.", sx.outFileName, buf, sx.fileType), "ab", h) ;
		  else
		    {
		      ccp = hprintf (h, "%s.%s.", sx.outFileName, buf) ;
		      sx.aos[ii] = aceOutCreate (ccp, sx.fileType, sx.gzo, h) ;
		    }
		}
	      if (! sx.aos[ii])
		    messcrash ("Cannot create split file : %s", sx.outFileName) ;
	    }
	}
    }

  sx.nSplit = 1 ;
  
  if (sx.keepNameBam) sx.keepName = TRUE ;
  if (sx.keepName1) sx.keepName = TRUE ;
  if (sx.keepName2) sx.keepName = TRUE ;

  if (sx.keepName)
    {
      switch (sx.out)
	{
	case CSFASTA:
	case CSFASTC:
	case CCFA:
	case CCFAR:
	case FASTA:
	case FASTC: 
	  if (sx.split || sx.splitMb)
	    ccp = hprintf (h, "%s.1.id%s", sx.outFileName, sx.gzo ? ".gz" : "") ;
	  else
	    ccp = hprintf (h, "%s.id%s", sx.outFileName, sx.gzo ? ".gz" : "") ;
	  
	  if (sx.gzo)
	    sx.aoId = aceOutCreateToGzippedFile (ccp, h) ;
	  else
	    sx.aoId = aceOutCreateToFile (ccp, "wb", h) ;
	  if (! sx.aoId)
	    exit (1) ;
	  break ;
	default:
	  break ;
	}
    }

  if (sx.split && ! sx.outIdFileName)
    sx.outIdFileName = "stdout" ;

  if (sx.gtfGenomeFastaFileName)
    { sx.inFileName = sx.gtfGenomeFastaFileName ; }

  if (sx.inFileName1 || sx.inFileName2)
    {
      if (! sx.inFileName1 || ! sx.inFileName2)
	messcrash ("To parse paired files you must specifiy -i1 file_name_1 and  -i2 file_name_2, \ntry -help\n") ;

      sx.ai1 = aceInCreate (sx.inFileName1, sx.gzi, h) ;
      if (! sx.ai1)
	  usage (messprintf("%s not found", sx.inFileName1)) ;
      
      if (sx.in == RAW || sx.in == CRAW || sx.in == TAG_COUNT || sx.in == CTAG_COUNT)
	aceInSpecial (sx.ai1,"\n") ;
      else
	aceInSpecial (sx.ai1,"\t\n") ;
      
      sx.ai2 = aceInCreate (sx.inFileName2, sx.gzi, h) ;
      if (! sx.ai2)
	  usage (messprintf("%s not found", sx.inFileName2)) ;
      
      if (sx.in == RAW || sx.in == CRAW || sx.in == TAG_COUNT || sx.in == CTAG_COUNT)
	aceInSpecial (sx.ai2,"\n") ;
      else
	aceInSpecial (sx.ai2,"\t\n") ;
      
    }
  else if (sx.inFileName)
    {
      sx.ai = aceInCreate (sx.inFileName, sx.gzi, h) ;
      if (! sx.ai)
	usage (messprintf("%s not found", sx.inFileName)) ;
    
      if (sx.in == RAW || sx.in == CRAW || sx.in == TAG_COUNT || sx.in == CTAG_COUNT)
	aceInSpecial (sx.ai,"\n") ;
      else
	aceInSpecial (sx.ai,"\t\n") ;
    }
  if (argc > 1)
    usage (messprintf("sorry unknown arguments %s", argv[argc-1])) ;

  getSelection (&sx) ;

  if (sx.getTm) 
    {
      sx.encodedDna = arrayHandleCreate (1000, char, h) ;
      sx.aoTm = aceOutCreate (sx.outFileName, ".TM.txt", sx.gzo, h) ;
      aceOutDate (sx.aoTm, "###", sx.title) ;
      aceOutf (sx.aoTm, "## N is the length of the sequence in base, GC the percentage of GC\n") ;
      aceOutf (sx.aoTm, "## The melting temperature is computed according to Maniatis Tm = 62.3 + 0.41*(%G+C) - 500/N\n") ;
      aceOutf (sx.aoTm, "## The dimer entropy is defined S = - Sum[over all ij dimers] (n_ij log16 n_ij/N)\n") ;
      aceOutf (sx.aoTm, "## Sequence\tLength\tEntropy\tTm\t%%GC\n") ;
    }
  if (sx.shadowFileName)
    parseShadowFile (&sx, 1) ;
  else if (sx.spongeFileName)
    parseShadowFile (&sx, 2) ;
  else if (sx.refGeneFileName)
    parseRefGeneFile (&sx) ;
  else if (sx.gtfFileName)
    parseGtfFile (&sx) ;

  if (! sx.gtfFileName && ! sx.inFileName && ! sx.inFileName1)
    sx.ai = aceInCreate (0, sx.gzi, h) ;

  if (sx.makeFastaUnique)
    {
      sx.ao = aceOutCreate (sx.outFileName, ".unique.fasta", sx.gzo, h) ;
      makeFastaUnique (&sx) ;
    }

  if (sx.makeTestGenome)
    {
      sx.aoTestSeqs = sx.ao ; /* aceOutCreate (sx.outFileName, ".modified.fasta", sx.gzo, h) ; */
      sx.aoTestSnps = aceOutCreate (sx.outFileName, ".snps", sx.gzo, h) ;
    }

  if (sx.doCount)
    sx.lnDistrib = keySetHandleCreate (h) ;
  if (! sx.gtfFileName || sx.gtfGenomeFastaFileName)
    dna2dnaRun (&sx) ;

  if (sx.doCount)
    {
      ACEOUT ao = aceOutCreate (sx.outFileName, ".count", FALSE, sx.h) ;
      aceOutDate (ao, "###", sx.title) ;

      aceOutf (ao, "Distinct_fragment_processed\t%ld\nDistinct_fragment_kept\t%ld\nDistinct_fragment_rejected\t%ld\n"
	       , sx.fProcessed + sx.fRejected, sx.fProcessed, sx.fRejected
	       ) ;
      aceOutf (ao, "Fragment_processed\t%ld\nFragment_kept\t%ld\nFragment_rejected\t%ld\n"
	       , sx.fmProcessed + sx.fmRejected, sx.fmProcessed, sx.fmRejected
	       ) ;
      aceOutf (ao, "Sequence_processed\t%ld\nSequence_kept\t%ld\nSequence_rejected\t%ld\n"
	       , sx.nProcessed, sx.nProcessed - sx.nRejected, sx.nRejected
	       ) ;
      aceOutf (ao, "Tags_processed\t%ld\nTags_kept\t%ld\nTags_rejected\t%ld\n"
	       , sx.tagProcessed, sx.tagProcessed - sx.tagRejected, sx.tagRejected
	       ) ;
      aceOutf (ao, "Bases_seq_processed\t%ld\nBases_seq_kept\t%ld\nBases_seq_rejected\t%ld\n"
	       , sx.bpSeqProcessed, sx.bpSeqProcessed - sx.bpSeqRejected, sx.bpSeqRejected
	       ) ;
      aceOutf (ao, "Bases_tags_processed\t%ld\nBases_tags_kept\t%ld\nBases_tags_rejected\t%ld\n"
	       , sx.bpTagsProcessed, sx.bpTagsProcessed - sx.bpTagsRejected, sx.bpTagsRejected
	       ) ;
      if (sx.leftClipOn || sx.leftClipInPhase)
	aceOutf (ao,  "Entry_adaptor %d %d %ld\n", sx.entryAdaptorSeqFound, sx.entryAdaptorTagFound, sx.entryAdaptorBpFound) ;
      if (sx.rightClipOn)
	aceOutf (ao,  "Exit_adaptor %d %d %ld\n", sx.exitAdaptorSeqFound, sx.exitAdaptorTagFound, sx.exitAdaptorBpFound) ;
      aceOutf (ao, "Min_probe_length\t%d\nMax_probe_length\t%d\n", sx.minLn, sx.maxLn
	       ) ; 
      aceOutf (ao, "Max_count\t%d\n", sx.maxCount) ;

      if (sx.lnDistrib)
	{
	  KEYSET ks = sx.lnDistrib ;
	  int ln1 = 0, ln5 = 0, ln50 = 0, ln95 = 0, ln99 = 0 ;
	  int i, iMax ;
	  int mode, maxN = 0 ;
	  long int cumul = 0,  n, nn, average ;

	  /* count the reads */
	  iMax = keySetMax (ks) ;
	  for (nn = 0, i = 0 ; i < iMax ; i++)
	    {
	      n = keySet (ks, i) ;
	      nn +=  n ;
	      cumul += i * n ;
	      if (n > maxN) { maxN  = n ; mode = i ; }
	    }
	  average = nn ? cumul / nn : 0 ;
	  /* find the median */
	  for (n = 0, i = 0 ; i < iMax ; i++)
	    { 
	      n += keySet (ks, i) ;
	      if (ln1 == 0 && 100*n >= nn) ln1 = i ;
	      if (ln5 == 0 && 100*n >= 5*nn) ln5 = i ;
	      if (ln50 == 0 && 100*n >= 50*nn) ln50 = i ;
	      if (ln95 == 0 && 100*n >= 95*nn) ln95 = i ;
	      if (ln99 == 0 && 100*n >= 99*nn) ln99 = i ;
	    }
	  aceOutf (ao, "Length_distribution_1_5_50_95_99_mode_av\t%d\t%d\t%d\t%d\t%d\t%d\t%ld\n", ln1, ln5, ln50, ln95, ln99, mode, average) ;
	}
      ac_free (ao) ;
      
    }
  if (0)
    {
      fprintf (stderr, "Sequence_processed\t%ld\nSequence_kept\t%ld\nSequence_rejected\t%ld\n"
	       , sx.nProcessed, sx.nProcessed - sx.nRejected, sx.nRejected
	       ) ;
      fprintf (stderr, "Tags_processed\t%ld\nTags_kept\t%ld\nTags_rejected\t%ld\n"
	       , sx.tagProcessed, sx.tagProcessed - sx.tagRejected, sx.tagRejected
	       ) ;
      fprintf (stderr, "Bases_seq_processed\t%ld\nBases_seq_kept\t%ld\nBases_seq_rejected\t%ld\n"
	       , sx.bpSeqProcessed, sx.bpSeqProcessed - sx.bpSeqRejected, sx.bpSeqRejected
	       ) ;
      fprintf (stderr, "Bases_tags_processed\t%ld\nBases_tags_kept\t%ld\nBases_tags_rejected\t%ld\n"
	       , sx.bpTagsProcessed, sx.bpTagsProcessed - sx.bpTagsRejected, sx.bpTagsRejected
	       ) ;
      if (sx.leftClipOn || sx.leftClipInPhase)
	fprintf (stderr, "Entry_adaptor %d %d %ld\n", sx.entryAdaptorSeqFound, sx.entryAdaptorTagFound, sx.entryAdaptorBpFound) ;
      if (sx.rightClipOn)
	fprintf (stderr, "Exit_adaptor %d %d %ld\n", sx.exitAdaptorSeqFound, sx.exitAdaptorTagFound, sx.exitAdaptorBpFound) ;
      fprintf (stderr, "Min_probe_length\t%d\nMax_probe_length\t%d\n", sx.minLn, sx.maxLn
	       ) ;
      fprintf (stderr, "Max_count\t%d\n", sx.maxCount) ;
    } 
  fprintf (stderr , "Distinct_processed %ld fragments\n", sx.fProcessed) ;
  fprintf (stderr , "Processed %ld fragments\n", sx.fmProcessed) ;
  fprintf (stderr , "Processed %ld sequences\t%ld bp\n", sx.nProcessed, sx.bpSeqProcessed) ;
  fprintf (stderr , "Processed %ld reads\t%ld bp\n", sx.tagProcessed, sx.bpTagsProcessed) ;
  {
    int mx ;
    messAllocMaxStatus (&mx) ;   
    fprintf (stderr, "// done: %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
  }
  {
    int  i = 256 ;
    while (i--)
      {
	ac_free (sx.aos[i]) ;
	ac_free (sx.aos1[i]) ;
	ac_free (sx.aos2[i]) ;
      }
  }
  ac_free (h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

/*

gunzip -c f.3.id.gz | gawk -F '\t'  '{for(i=2;i<=NF;i++){split($i,aa,":");c[i]=aa[5];x[i]=aa[6];y[i]=aa[7];if (0)print "cxy",i,c[i],x[i],y[i], $1,aa[2];}for (i=2;i<NF;i++){for(j=i+1;j<=NF;j++){if(c[i]==c[j]){dx=x[i]-x[j];dy=y[i]-y[j];if(dx<0)dx=-dx;if(dy<0)dy=-dy;if(dxmax<dx)dxmax=dx;if(dymax<dy)dymax=dy;r=int(.5+sqrt(dx*dx+dy*dy));if(rmax<r)rmax=r;if(r>30)r=30;if(dx>30)dx=30;if(dy>30)dy=30;nnn++;nn[dx,dy]++;nx[dx]++;ny[dy]++;rrt++;rr[r]++;if(0) print i,j,y[i],y[j],dx,dy;}}}}END{printf("Run %s Number of pairs of identical reads at given radial distance, saturated at 30\n",run);for(i=1;i<=30;i++)printf("\t%d",i);printf("\n%s",run);for(i=1;i<=30;i++)printf("\t%d",rr[i]);printf("\t%d\t%g\n\nRun %s Number of pairs of identical reads at given (x,y) distance, saturated at 30, in reality dxmax=%d  dymax=%d\n",rrt-rr[30],100*rrt/(rr[30]+rrt),run,dxmax,dymax);for(i=1;i<=30;i++)printf("\t%d",i);printf("\tTotal");for(i=1;i<=30;i++){printf("\n%d",i);for(j=1;j<=30;j++)printf("\t%d",nn[i,j]);printf("\t%d",nx[i]);}printf("\nTotal");for(j=1;j<=30;j++)printf("\t%d",ny[j]);printf("\n");}' run=Rprim1/f.3 | tee Rprim1.f.3.xy.txt




*/
