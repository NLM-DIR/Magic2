/* seqc.c
 * Authors: Jean Thierry-Mieg, NCBI, Joe Meehan, NCTR
 * Oct 2012
 * SAM search in the output of the CALI suite
*/

/* #define ARRAY_CHECK      */

#include "../wac/ac.h"
#include "keyset.h"
#include "query.h"
#include "dict.h"
#include "aceio.h"
#include <wiggle.h>

static void usage (char *message) ;

/*************************************************************************************/

#define SAM_PAIRED_END 0x1
#define SAM_ALIGNED 0x2
#define SAM_UNALIGNED 0x4
#define SAM_MATE_UNALIGNED 0x8
#define SAM_REVERSE 0x10
#define SAM_MATE_REVERSE 0x20
#define SAM_FIRST 0x40
#define SAM_LAST 0x80
#define SAM_SECONDARY 0x100
#define SAM_REJECTED 0x200
#define SAM_DUPLICATE 0x400

typedef struct samRecordStruct { int flag, target, pos, matePos, mapq, targetLn ; char cigar[32000] ; Array dna ; } SR ;
typedef struct samCountStruct {
  long int  record_count
  , perfect, perfect1, perfect2
    , perfectClipped, perfectClipped1, perfectClipped2
    , reads_unaligned, reads_unaligned1, reads_unaligned2
    , reads_in_file, reads_in_file1, reads_in_file2
    , read_bp, read_bp1, read_bp2
    , bp_aligned, bp_aligned1, bp_aligned2
    , bp_unaligned, bp_unaligned1, bp_unaligned2
    , reads_mapped, reads_mapped1, reads_mapped2
    , reads_mapped_on_plus_strand, reads_mapped_on_minus_strand
    , reads_mapped_on_plus_strand1, reads_mapped_on_minus_strand1
    , reads_mapped_on_plus_strand2, reads_mapped_on_minus_strand2
    , ali, ali1, ali2
    , oneAli, oneAli1, oneAli2
    , pairs_mapped, pairs_compatible, pairs_compatible_topology, pairs_compatible_distance_below_1Mb
    
    , soft_clips, nIntrons
    , mismatches, mismatches1, mismatches2
    , substitutions, substitutions1, substitutions2
    , transitions, transitions1, transitions2
    , transversions, transversions1, transversions2
    , insertions, insertions1, insertions2
    , double_insertions, double_insertions1, double_insertions2
    , triple_insertions, triple_insertions1, triple_insertions2

    , delete_len, delete_len1, delete_len2
    , deletions, deletions1, deletions2
    , double_deletions, double_deletions1, double_deletions2
    , triple_deletions, triple_deletions1, triple_deletions2
    ;
  Array errPos, errPos1, errPos2
    , nPos, nPos1, nPos2
    , errType, errType1, errType2
    , errCount, errCount1, errCount2
    , spikes, spikes1, spikes2, introns
    , aliLn, aliLn1, aliLn2
    , multiAli, multiAli1, multiAli2
    ;
} SC ;

typedef struct intronStruct { 
  int chrom, a1, a2, length ;
  int support ;
  char feet[6] ;
  Array tags ;
} INTRON ;

typedef struct rcStruct { 
  int run ;
  Array tags ;
} RC ;

typedef struct tcStruct { 
  int tag ;
  int cols ;
  long int zz[9] ;
} TC ;

typedef struct samStruct {  
  AC_HANDLE h ; BOOL gzi, gzo ;
  const char *inFileName, *outFileName, *title, *fastaFileName ;
  const char *run ;
  int tagMult, n, NN, count, sam ; 

  SC *sc ;
  BOOL report, merge, aceOut, isSorted, BAM, SAM, isFastc ;
  DICT *targetDict, *runDict, *tagDict ;
  DICT *eeDict, *spikeDict ;

  Array errorsP, sams, dnas, chroms, runs ;
  ACEIN ai ; ACEOUT ao ; 
 } SAM ;

/*************************************************************************************/
typedef struct chromStruct { Array dna ; int nam ; } CHROM ;
 
static void samParseFastaRegister (CHROM *chrom, Stack s, int nbp) 
{
  chrom->dna = arrayCreate (nbp+1, char ) ;
  
  array (chrom->dna, nbp+1,char) = 0 ;
  array (chrom->dna, nbp,char) = 0 ;
  arrayMax (chrom->dna) = nbp ;
  memcpy (arrp(chrom->dna, 0, char), stackText(s,0), nbp) ;
  dnaEncodeArray (chrom->dna) ;
  stackClear (s) ;
} /* parseFastaRegister */

/***************************/

static int samParseFastaFile (SAM *sam)
{
  AC_HANDLE h = ac_new_handle () ;
  char *cp ;
  ACEIN ai ;
  int state = 0, n, nn = 0, line = 0, nChrom = 0, pos = 0, target = 0 ;
  long  nbp = 0 ;
  CHROM *chrom = 0 ;
  Stack s = stackHandleCreate (100000, h) ;

  sam->chroms = arrayHandleCreate (1024, CHROM, sam->h) ;

  ai = aceInCreate (sam->fastaFileName, FALSE, h) ;
  aceInSpecial (ai, "\n") ;

  state = nn = 0 ;
  while (aceInCard (ai)) /* will close pp->targetFile */
    {
      line++ ;
      cp = aceInWord (ai) ;
      if (cp && *cp == '#') continue ;
      if (cp) switch (state)
	{
	case 1:
	case 2:
	  if (*cp != '>')
	    {
	      n = strlen (cp) ;
	      nbp += n ;
	      if (nbp > 400000000)
		messcrash ("Sorry the individual reference sequences are limited to 400M") ;
	      if (pos)
		catText (s, cp) ;
	      else
		pushText (s, cp) ;
	      pos += n ;
	      break ;
	    }
      /* fall thru to new sequence */
	case 0: /* expecting   >target_name */ 
	  if (*cp == '>')
	    {
	      if (chrom && nbp)
		samParseFastaRegister (chrom, s, nbp) ;

	      pos = nbp = 0 ;

	      if (! *(cp+1))
		{
		  fprintf (stderr, "// bad character in file %s line %d, expecting a '>', got %s\n"
			   , sam->fastaFileName, line, cp) ;
		  return FALSE ;
		}

	      dictAdd (sam->targetDict, cp+1, &target) ;
	      chrom = arrayp (sam->chroms, target, CHROM) ;
	      chrom->nam = target ;
	      
	      state = 1 ;
	    }
	  break ;
	}
    }

  if (chrom && nbp)
    samParseFastaRegister (chrom, s, nbp) ;

  ac_free (h) ;
  fprintf (stderr, "// %s samParseFastaFile selected %d sequences %ldbp\n", timeShowNow(), nn, nbp) ; 
  return nChrom ;
} /* samParseFastaFile */

/*************************************************************************************/
/*************************************************************************************/

/* all this SS SSM code is a remnant from snp.c, we may need it so i do not erase it */
#define NN 32
typedef struct ssStruct { int pos, count ; short sam, qual, strand ; } SS ;

/*************************************************************************************/

static DICT *eeDictCreate (SAM *sam)
{
  DICT *eeDict = dictHandleCreate (10000, sam->h) ;
  DICT *spikeDict = dictHandleCreate (10000, sam->h) ;

  sam->eeDict = eeDict ;
  sam->spikeDict = spikeDict ;
  /* thse initialisations will drive the order in which the errors are reported */
  /* force wild type single letters fisrt */
  dictAdd (spikeDict, "toto", 0) ;
  dictAdd (eeDict, "toto", 0) ;
  dictAdd (eeDict, "a", 0) ;
  dictAdd (eeDict, "t", 0) ;
  dictAdd (eeDict, "g", 0) ;
  dictAdd (eeDict, "c", 0) ;
  dictAdd (eeDict, "n", 0) ;
  dictAdd (eeDict, "W", 0) ;
  
  dictAdd (eeDict, "a>t", 0) ;
  dictAdd (eeDict, "a>g", 0) ;
  dictAdd (eeDict, "a>c", 0) ;
  dictAdd (eeDict, "a>n", 0) ;
  
  dictAdd (eeDict, "t>a", 0) ;
  dictAdd (eeDict, "t>g", 0) ;
  dictAdd (eeDict, "t>c", 0) ;
  dictAdd (eeDict, "t>n", 0) ;
  
  dictAdd (eeDict, "g>a", 0) ;
  dictAdd (eeDict, "g>t", 0) ;
  dictAdd (eeDict, "g>c", 0) ;
  dictAdd (eeDict, "g>n", 0) ;
  
  dictAdd (eeDict, "c>a", 0) ;
  dictAdd (eeDict, "c>t", 0) ;
  dictAdd (eeDict, "c>g", 0) ;
  dictAdd (eeDict, "c>n", 0) ;
  
  dictAdd (eeDict, "+a", 0) ;
  dictAdd (eeDict, "+t", 0) ;
  dictAdd (eeDict, "+g", 0) ;
  dictAdd (eeDict, "+c", 0) ;
  
  dictAdd (eeDict, "*+a", 0) ;
  dictAdd (eeDict, "*+t", 0) ;
  dictAdd (eeDict, "*+g", 0) ;
  dictAdd (eeDict, "*+c", 0) ;
  
  dictAdd (eeDict, "-a", 0) ;
  dictAdd (eeDict, "-t", 0) ;
  dictAdd (eeDict, "-g", 0) ;
  dictAdd (eeDict, "-c", 0) ;
  
  dictAdd (eeDict, "*-a", 0) ;
  dictAdd (eeDict, "*-t", 0) ;
  dictAdd (eeDict, "*-g", 0) ;
  dictAdd (eeDict, "*-c", 0) ;

  return eeDict ;
} /* eeDictCreate */

/*************************************************************************************/
/*************************************************************************************/

static int samIntronOrder (const void *va, const void *vb)
{
  const INTRON *a = (const INTRON *)va, *b = (const INTRON *)vb ;
  int n ;

  n = a->chrom - b->chrom ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  return 0 ;
} /* samIntronOrder */

/*************************************************************************************/

static int srsOrder (const void *va, const void *vb)
{
  const SR *a = (const SR *)va, *b = (const SR *)vb ;
  int n1, n2 ;

  n1 = a->flag & SAM_REJECTED ;
  n2 = b->flag & SAM_REJECTED ;
  if (n1 != n2) return n1 - n2 ; /* baddy last */

  n1 = a->flag & SAM_SECONDARY ;
  n2 = b->flag & SAM_SECONDARY ;
  if (n1 != n2) return n1 - n2 ; /* baddy last */

  n1 = a->flag & SAM_DUPLICATE ;
  n2 = b->flag & SAM_DUPLICATE ;
  if (n1 != n2) return n1 - n2 ; /* baddy last */

  n1 = a->flag & SAM_ALIGNED ;
  n2 = b->flag & SAM_ALIGNED ;
  if (n1 != n2) return n2 - n1 ; /* aligned first */

  n1 = a->flag & SAM_FIRST ;
  n2 = b->flag & SAM_FIRST ;
  if (n1 != n2) return n2 - n1 ; /* first first */

  return 0 ;
} /* srsOrder */

/*************************************************************************************/
/*************************************************************************************/

static BOOL samSmokeCigar (SAM *sam, SR *sr, BOOL isFirst, BOOL isFirst1, BOOL isFirst2, int bp)
{
  AC_HANDLE h = ac_new_handle () ;
  SC *sc = sam->sc ;
  Array mapped_ref = arrayHandleCreate (2*bp, char, h) ;
  CHROM *chrom = arrayp (sam->chroms, sr->target, CHROM) ;
  Array refDna = chrom->dna ;
  char *cigar, *cpMapped_ref, *cpM, *cpRefDna, *cpR , *cpRead ;
  int i, mLength = 0, pos, posR = 0, nErr = 0, nClipped = 0 ;
  BOOL isPerfect = TRUE, hasSoftClip = FALSE ;
  BOOL isReverse =  (sr->flag & SAM_REVERSE) ? TRUE : FALSE ;  /* we align on minus strand */
  char buf[16] ;

  array (mapped_ref, 2*bp - 1, char) = 0 ; /* make room*/
  cpM = cpMapped_ref = arrp (mapped_ref, 0, char) ;
  cpR = cpRefDna = arrp (refDna, sr->pos-1, char) ; 
  cpRead = arrp (sr->dna, 0, char) ;
  int readlen = arrayMax (sr->dna); /* jfm 10/22/2012 Read length. Used to calculate positions on the minus strand */

  /* Loop through the CIGAR objects to get counts and build the mapped reference string (mapped_ref) */
 
  for (cigar = sr->cigar ; *cigar ; cigar++)
    {
      int mult ;

      /* grab the number */
      for ( mult = 0 ; *cigar && *cigar >= '0' && *cigar <= '9' ; cigar++)
	mult = 10 * mult + (*cigar - '0') ; /* parse a decimal number */
      if (mult == 0)  /* i imagine there may be cases without multipliers between letters */
	mult = 1 ; 
      
      switch ((int) *cigar)  /* parse the CIGAR types */
	{
	case 0:   /* cigar string is over this should not happen after a number */
	  cigar-- ; /* the cigar++ on loop end will reposition on *cigar==0 */
	  break ;

	case 'S': /* Soft Clip (clipped sequences present in SEQ) */
	  hasSoftClip = TRUE ;
	  sc->soft_clips += mult ;
	  /* blank the read */
	  for (i = 0 ; i < mult ; i++)
	    {
	      /* *cpM++ = '-' ; *cpR++ = '-' ; jfm 10/22/2012 */
              *cpM++ = '-' ;
	    }
	  nClipped += mult ;
	  mLength += mult ;
	  posR += mult ;
	  break ;

	case 'M': /* OUARFFFFF, I can't beleive they are so dumb: as of 2012  M can be a sequence match or mismatch */
	case '=': /* this should be fixed in the future */
	case 'X': /* = and X are recent additions to the SAM format, not used by the current the aligners */
	  mLength += mult ;
	  posR += mult ;
	  while (mult--)
	    *cpM++ = *cpR++ ; /* copy over the correct number of (possibly) matching bases */
	  break ;


	case 'N': /*  Skipped region from the reference (for mRNA to genome alignment, represents an intron) */
	case 'D': /* deletions from the reference */
	  isPerfect = FALSE ;
	  /* move forward on the reference */
	  /*  cpR += mult ; */ /* jfm 10/22/2012 Increment after counting deletion types */

	  /* count deletion events as errors at this position (jfm 10/5/2012) */
	  pos = mLength ; /* we are building the mapped reference, so current position is it's length */
	  if (isReverse)  /* we align on munus strand */
	    pos = readlen - pos - 1 + 1 ; /* jfm 10/22/2012 use readlen rather than mLength */ 


 	  isPerfect = FALSE ;

	  if (mult > 20)
	    {
	      INTRON *iip ;
	      Array aa = sam->sc->introns ;

	      sc->nIntrons += 1 ;
	      iip = arrayp (aa, arrayMax (aa), INTRON) ;
	      iip->chrom = sr->target ;
	      if (! isReverse)
		{
		  iip->a1 = cpR - cpRefDna + sr->pos - 1 ;
		  iip->a2 = iip->a1 + mult + 1 ;
		  iip->length = mult ;
		  iip->feet[0] = dnaDecodeChar[(int)cpR[0]] ;
		  iip->feet[1] = dnaDecodeChar[(int)cpR[1]] ;
		  iip->feet[2] = '_' ;
		  iip->feet[3] = dnaDecodeChar[(int)cpR[mult-2]] ;
		  iip->feet[4] = dnaDecodeChar[(int)cpR[mult-1]] ;
		  iip->feet[5] = 0 ;
		}
	      else if (1)
		{
		  iip->a1 = cpR - cpRefDna + sr->pos - 1 ;
		  iip->a2 = iip->a1 + mult + 1 ;
		  iip->length = mult ;
		  iip->feet[0] = dnaDecodeChar[(int)cpR[0]] ;
		  iip->feet[1] = dnaDecodeChar[(int)cpR[1]] ;
		  iip->feet[2] = '_' ;
		  iip->feet[3] = dnaDecodeChar[(int)cpR[mult-2]] ;
		  iip->feet[4] = dnaDecodeChar[(int)cpR[mult-1]] ;
		  iip->feet[5] = 0 ;
		}

	      iip->length = mult ;
	      iip->support = 1 ;
	    }

	  if (1)
	    {
	      sc->deletions ++ ; /* Count deletion events, doublet deletions and triplet deletions */
	      sc->delete_len += mult ;
	      if (mult == 2)
		sc->double_deletions++ ;
	      else if (mult == 3)
		sc->triple_deletions++ ;
	      array (sc->errPos, pos, int) ++ ;
	    }

	  /* Count deletions on first and second fragments */
	  if (isFirst1)
	    {
	      array (sc->errPos1, pos, int) ++ ;
	      sc->deletions1 ++ ; /* Count deletion events, doublet deletions and triplet deletions */
	      sc->delete_len1 = mult ;
	      if (mult == 2)
		sc->double_deletions1++ ;
	      else if (mult == 3)
		sc->triple_deletions1++ ;
	    }
	  if (isFirst2)
	    {
	      array (sc->errPos2, pos, int) ++ ;
	      sc->deletions2 ++ ; /* Count deletion events, doublet deletions and triplet deletions */
	      sc->delete_len2 = mult ;
	      if (mult == 2)
		sc->double_deletions2++ ;
	      else if (mult == 3)
		sc->triple_deletions2++ ;
	    }

	  /* count deletion types */
	  if (mult < 4)
	    {
	      int i, type = 0 ;

	      nErr++ ; /* was called 	mismatch_rate, but it is not a rate */

	      buf[4+mult] = 0 ;
	      strcpy (buf, "Del:") ;
	      if (isReverse)
		{
		  for (i = 0 ; i < mult ; i++)
		    buf[4+i] = ace_upper(dnaDecodeChar[(int)complementBase(cpR[mult - i - 1])]) ; 
		}
	      else
		{
		  for (i = 0 ; i < mult ; i++)
		    buf[4+i] = ace_upper(dnaDecodeChar[(int)cpR[i]]) ; 
		}
	      dictAdd (sam->eeDict, buf, &type) ;
	      array (sc->errType, type, int)++ ;

	      if (isFirst1)
		array (sc->errType1, type, int)++ ;
	      if (isFirst2)
		array (sc->errType2, type, int)++ ;

	      if (1)
		{
		  char *spikeBuf = messprintf ("%s:%d", buf, pos) ;
		  int spike ;

		  dictAdd (sam->spikeDict, spikeBuf, &spike) ;
		  array (sc->spikes, spike, int)++ ;
		  if (isFirst1)
		    array (sc->spikes1, spike, int)++ ;
		  if (isFirst2)
		    array (sc->spikes2, spike, int)++ ;
		}

	    }

          /* move forward on the reference */
          cpR += mult ;
	  break ;

	case 'I': /* deletions to the reference */
	  isPerfect = FALSE ;
	  /* include padding */
	  for (i = 0 ; i < mult ; i++)
	    *cpM++ = '-' ;

	  /* count insertion events at this position (jfm 10/5/2012)
	   * do this before adding to mapped_ref to preserve correct position (jfm 10/16/2012)
	   */
	  pos = mLength ; /* we are building the mapped reference, so current position is it's length */
	  if (isReverse)  /* we align on munus strand */
	    pos = readlen - pos - 1 - 1 ; /* jfm 10/22/2012 Use readlen rather than mLength */
                                      /* jfm 10/22/2012 changed from readlen + 1 - pos */

	  if (1)
	    {
	      array (sc->errPos, pos, int) ++ ;
	      sc->insertions ++ ; /* Count insertion events, doublet insertions and triplet insertions */
	      sc->delete_len += mult ;
	      if (mult == 2)
		sc->double_insertions++ ;
	      else if (mult == 3)
		sc->triple_insertions++ ;
	    }

	  /* Count insertions on first and second fragments */
	  if (isFirst1)
	    {
	      array (sc->errPos1, pos, int) ++ ;
	      sc->insertions1 ++ ; /* Count insertion events, doublet insertions and triplet insertions */
	      sc->delete_len1 = mult ;
	      if (mult == 2)
		sc->double_insertions1++ ;
	      else if (mult == 3)
		sc->triple_insertions1++ ;
	    }
	  if (isFirst2)
	    {
	      array (sc->errPos2, pos, int) ++ ;
	      sc->insertions2 ++ ; /* Count insertion events, doublet insertions and triplet insertions */
	      sc->delete_len2 = mult ;
	      if (mult == 2)
		sc->double_insertions2++ ;
	      else if (mult == 3)
		sc->triple_insertions2++ ;
	    }

	  /* count insertion types */
	  if (mult < 4)
	    {
	      int i, type = 0 ;

	      nErr++ ; /* was called 	mismatch_rate, but it is not a rate */

	      buf[4+mult] = 0 ;
	      strcpy (buf,"Ins:") ;
	      if (isReverse)
		{
		  for (i = 0 ; i < mult ; i++)
		    buf[4+i] = ace_upper(dnaDecodeChar[(int)complementBase(cpRead[posR + mult - 1 - i])]) ;
		}
	      else
		{
		  for (i = 0 ; i < mult ; i++)
		    buf[4+i] = ace_upper(dnaDecodeChar[(int)cpRead[posR + i]]) ;
		}
	      /* posR += mult ;*/

              if ( 0 && mult == 2 && buf[4] == 'N' && buf[5] == 'G' ) { /* Joe's simple minded debug */
                printf ("-----\n");
                printf ("%s\n", sr->cigar); 
                for(i=0; cpRead[i]; i++)
                  putchar( dnaDecodeChar[(int)cpRead[i]] );
                putchar('\n');
                for(i=0; cpMapped_ref[i]; i++) 
                  putchar( dnaDecodeChar[(int)cpMapped_ref[i]] );
                putchar('\n');
              }


	      dictAdd (sam->eeDict, buf, &type) ;
	      array (sc->errType, type, int)++ ;

	      if (isFirst1)
		array (sc->errType1, type, int)++ ;
	      if (isFirst2)
		array (sc->errType2, type, int)++ ;

	      if (1)
		{
		  char *spikeBuf = messprintf ("%s:%d", buf, pos) ;
		  int spike ;

		  dictAdd (sam->spikeDict, spikeBuf, &spike) ;
		  array (sc->spikes, spike, int)++ ;
		  if (isFirst1)
		    array (sc->spikes1, spike, int)++ ;
		  if (isFirst2)
		    array (sc->spikes2, spike, int)++ ;
		}

	    }

          posR += mult;     /* jfm 11/01/2012 Increment position on read for all inserts */
          mLength += mult ; /* jfm 11/1/2012 Increment mapped length on insert */

	  break ;


	default: /* Something  different */
	  messcrash ("unexpected letter in CIGAR %s at line %d file %s"
		     , sr->cigar, aceInStreamLine (sam->ai), sam->inFileName
		     ) ;
	  break ;
	}
    }

  /* ok, we have reconstructed mapped_ref */
  arrayMax (mapped_ref) = mLength ;  /* in priciple we should have use array() to auto adjust the size of mapped_ref
				      *	but we are in an extremely tight loop running billions of times
				      */

  if (1)
    {
      int k, k1 ;
      
      k1 = arrayMax (sr->dna) ; /* len(a.read) */
      k = k1 - nClipped ; /* aligned_length */

      if (1)
	{
	  sc->bp_aligned += k ;
	  sc->read_bp += k1 ;
	  array (sc->aliLn, k, int)++ ;
	}
      if (isFirst1)
	{
	  sc->bp_aligned1 += k ;
	  sc->read_bp1 += k1 ;
	  array (sc->aliLn1, k, int)++ ;
	}
      if (isFirst2)
	{
	  sc->bp_aligned2 += k ;
	  sc->read_bp2 += k1 ;
	  array (sc->aliLn2, k, int)++ ;
	}
    }

  cpRead = arrp (sr->dna, 0, char) ;
  cpM = cpMapped_ref ;

  if ( 0 ) { /* Joe's simple minded debug */
    printf ("-----\n");
    for(i=0; cpRead[i]; i++)
      putchar( dnaDecodeChar[(int)cpRead[i]] );
    putchar('\n');
    for(i=0; cpM[i]; i++)
      putchar( dnaDecodeChar[(int)cpM[i]] );
    putchar('\n');
  }

  /* scan for mismatches the mapped sequence and the original read in parallel, avoiding paddings */
  for ( i = 0 ; *cpM && *cpRead ; cpM++, cpRead++, i++)
    {
      if (*cpRead == N_)
	{
	  pos = isReverse ? arrayMax(sr->dna) - i - 1 : i;
	  array (sc->nPos, pos, int)++ ;
	  if (isFirst1)
	    array (sc->nPos1, pos, int)++ ;
	  if (isFirst2)
	    array (sc->nPos2, pos, int)++ ;
	}
      /*      else if (*cpM != *cpRead && *cpM != '-' && *cpM != N_) */
      /* jfm 10/31/2012 Only count substitutions when both read and reference are unambiguous bases. (to match the Python version) */
      else if ( (*cpM != *cpRead) && ( *cpM == A_ || *cpM == T_ || *cpM == C_ || *cpM == G_ ) && ( *cpRead == A_ || *cpRead == T_ || *cpRead == C_ || *cpRead == G_ ))
	{
	  nErr++ ;
	  pos = isReverse ? arrayMax(sr->dna) - i - 1 : i;
	  
	  strcpy (buf, "Sub:") ;
	  if (isReverse)
	    {
	      buf[4] = ace_upper(dnaDecodeChar[(int)complementBase(*cpM)]) ;
	      buf[5] = ace_upper(dnaDecodeChar[(int)complementBase(*cpRead)]) ;
	    }
	  else
	    {
	      buf[4] = ace_upper(dnaDecodeChar[(int)*cpM]) ;
	      buf[5] = ace_upper(dnaDecodeChar[(int)*cpRead]) ;
	    }
	  buf[6] = 0 ;
	  
	  if (1)
	    {
	      int type, spike ;
	      BOOL isTs = FALSE ;
	      char *spikeBuf = messprintf ("%s:%d", buf, pos) ;
	      char **t, *transition_list[] = {"AG","GA","CT","TC", 0} ;
	      
	      sc->substitutions++ ;
	      if (isFirst1)
		sc->substitutions1++ ;
	      if (isFirst2)
		sc->substitutions2++ ;
	      
	      for (t = transition_list ; *t ; t++)
		if (! strcmp (*t, buf+4))
		  {
		    isTs = TRUE ;
		    break ;
		  }
	    
	      if (isTs)
		{
		  sc->transitions++ ;
		  if (isFirst1)
		    sc->transitions1++ ;
		  if (isFirst2)
		    sc->transitions2++ ;
		}
	      else
		{
		  sc->transversions++ ;
		  if (isFirst1)
		    sc->transversions1++ ;
		  if (isFirst2)
		    sc->transversions2++ ;
		}
	      
	    dictAdd (sam->eeDict, buf, &type) ;
	    array (sc->errType, type, int)++ ;
	    if (isFirst1)
	      array (sc->errType1, type, int)++ ;
	    if (isFirst2)
	      array (sc->errType2, type, int)++ ;
	    
	    dictAdd (sam->spikeDict, spikeBuf, &spike) ;
	    array (sc->spikes, spike, int)++ ;
	    if (isFirst1)
	      array (sc->spikes1, spike, int)++ ;
	    if (isFirst2)
	      array (sc->spikes2, spike, int)++ ;
	    
	    array (sc->errPos, pos, int)++ ;
	    if (isFirst1)
	      array (sc->errPos1, pos, int)++ ;
	    if (isFirst2)
	      array (sc->errPos2, pos, int)++ ;
	    }
	}
    }
  /* Register the mismatch count histogram */
  if (1)
    {
      if (nErr)
	isPerfect = FALSE ;
      
      sc->mismatches += nErr ;
      array (sc->errCount, nErr, int)++ ;
      if (isFirst1)
	{
	  array (sc->errCount1, nErr, int)++ ;
	  sc->mismatches1 += nErr ; 
	}
      if (isFirst2)
	{
	  array (sc->errCount2, nErr, int)++ ;
	  sc->mismatches2 += nErr ; 
	}
    }
  
  /* register the perfect reads */
  if (! nErr && isPerfect)
    {
      if (1)
	{
	  if (hasSoftClip)
	    sc->perfectClipped++ ;
	  else
	    sc->perfect++ ;
	}
      if (isFirst1)
	{
	  if (hasSoftClip)
	    sc->perfectClipped1++ ;
	  else
	    sc->perfect1++ ;
	}
      if (isFirst2)
	{
	  if (hasSoftClip)
	    sc->perfectClipped2++ ;
	  else
	    sc->perfect2++ ;
	}
    }
  
  
  ac_free (h) ;
  return TRUE ;
}/*  samSmokeCigar */

/*************************************************************************************/

static BOOL samCountOneAlignedRecord (SAM *sam, SR *sr, BOOL isFirst, BOOL isFirst1, BOOL isFirst2)
{
  int bp = sr->dna ? arrayMax (sr->dna) : 0 ;
  SC *sc = sam->sc ;

   sc->ali++ ;

   sc->oneAli++ ;
   if (sr->flag & SAM_PAIRED_END)
     {
       if (sr->flag & SAM_FIRST)
	 {
	   sc->ali1++ ;
	   sc->oneAli1++ ;
	   isFirst2 = FALSE ;
	   if (isFirst1)
	     {
	       sc->reads_in_file++ ;
	       sc->reads_mapped++ ;

	       sc->reads_in_file1++ ;
	       sc->reads_mapped1++ ;
	       
	       if (sr->flag &  SAM_REVERSE)
		 {
		   sc->reads_mapped_on_minus_strand++ ;
		   sc->reads_mapped_on_minus_strand1++ ;
		 }
	       else
		 {
		   sc->reads_mapped_on_plus_strand++ ;
		   sc->reads_mapped_on_plus_strand1++ ;
		 }
	     }
	 }
       else
	 {
	   sc->ali2++ ;
	   sc->oneAli2++ ;
	   isFirst1 = FALSE ;
	   if (isFirst2)
	     {
	       sc->reads_in_file++ ;
	       sc->reads_mapped++ ;

	       sc->reads_in_file2++ ;
	       sc->reads_mapped2++ ;

	       if (sr->flag &  SAM_REVERSE)
		 {
		   sc->reads_mapped_on_minus_strand++ ;
		   sc->reads_mapped_on_minus_strand2++ ;
		 }
	       else
		 {
		   sc->reads_mapped_on_plus_strand++ ;
		   sc->reads_mapped_on_plus_strand2++ ;
		 }
	     }
	 }
     }	
   else 
     {
       isFirst1 = isFirst2 = FALSE ;
       if (isFirst)
	 {
	   sc->reads_in_file++ ;
	   sc->reads_mapped++ ;
	       if (sr->flag &  SAM_REVERSE)
		 {
		   sc->reads_mapped_on_minus_strand++ ;
		 }
	       else
		 {
		   sc->reads_mapped_on_plus_strand++ ;
		 }
	 }
     }

   if (isFirst || isFirst1 || isFirst2)
     {
       if ( !(sr->flag & 0x04) ) /* Fragment is mapped (not unmapped) */
         {
	   if ( !(sr->flag & 0x08) ) /* Mate mapped (not unmapped) */
	     {
	       sc->pairs_mapped++;
	       if ( sr->flag & 0x02 ) /* Proper pair (compatible) */
		 sc->pairs_compatible++;
	     }
         }
       
       if (isFirst1 && (!(sr->flag & 0x04 ) &&  !(sr->flag & 0x08))) /* First fragment, not unmapped, mate not unmapped */
	 {
	   /* check topology and distance to define a compatible alignment */ 
	   if (
	       ((sr->flag & 0x30) == 0x20 && sr->matePos > sr->pos) ||
	       ((sr->flag & 0x30) == 0x10 && sr->matePos < sr->pos)
	       )
	     {
	       sc->pairs_compatible_topology++ ;
	       if (
		   ((sr->flag & 0x30) == 0x20 && sr->matePos > sr->pos &&  sr->matePos < sr->pos + 1000000 ) ||
		   ((sr->flag & 0x30) == 0x10 && sr->matePos < sr->pos &&  sr->matePos > sr->pos - 1000000 )
		   )
		 sc->pairs_compatible_distance_below_1Mb++ ;
	     }
	 }
       
       if (bp)
	 samSmokeCigar (sam, sr, isFirst, isFirst1, isFirst2, bp) ;
     }
   
   return TRUE ; 
} /* samCountOneRecord */

/*************************************************************************************/

static BOOL samCountOneUnalignedRecord (SAM *sam, SR *sr, BOOL isFirst, BOOL isFirst1, BOOL isFirst2)
{
  int bp = sr->dna ? arrayMax (sr->dna) : 0 ;
  SC *sc = sam->sc ;

  if (sr->flag & SAM_PAIRED_END)
    {
      if (isFirst1 && (sr->flag & SAM_FIRST))
	{
	  sc->reads_in_file++ ;
	  sc->reads_unaligned++ ;
	  sc->bp_unaligned += bp ;

	  sc->reads_in_file1++ ;
	  sc->reads_unaligned1++ ;
	  sc->bp_unaligned1 += bp ;
	}
      if (isFirst2 && (sr->flag & SAM_LAST))
	{
	  sc->reads_in_file++ ;
	  sc->reads_unaligned++ ;
	  sc->bp_unaligned += bp ;

	  sc->reads_in_file2++ ;
	  sc->reads_unaligned2++ ;
	  sc->bp_unaligned2 += bp ;
	}
    }	
  else if (isFirst)
    {
      sc->reads_in_file++ ;
      sc->reads_unaligned++ ;
      sc->bp_unaligned += bp ;
    }

  

  return TRUE ;
} /* samCountOneRecord */

/*************************************************************************************/

static void samCountOneFragment (SAM *sam, Array srs) 
{     
  SR *sr ;
  SC *sc = sam->sc ;
  BOOL isFirst, isFirst1, isFirst2 ;
  int i ;

  if (! arrayMax (srs))
    return ;

  arraySort (srs, srsOrder) ;

  sr = arrp (srs, 0, SR) ;
  isFirst = TRUE ;
  if (sr && sr->flag & SAM_PAIRED_END)
    isFirst1 = isFirst2  = TRUE ;
  else
    isFirst1 = isFirst2  = FALSE ;

  for (i = 0 ; i < arrayMax (srs) ; i++)
    {
      sr = arrp (srs, i, SR) ;

      sam->sc->record_count++ ;

      if (sr->flag & SAM_REJECTED)
	continue ;
      /*
	if (sr->flag & SAM_DUPLICATE)
	continue ;
      */

      if (! (sr->flag &  SAM_UNALIGNED)) /* jfm 10/22/2012 Was sr->flag &  SAM_ALIGNED */
	samCountOneAlignedRecord (sam, sr, isFirst, isFirst1, isFirst2) ;
      else
	samCountOneUnalignedRecord (sam, sr, isFirst, isFirst1, isFirst2) ;
      isFirst = FALSE ;
      if (sr->flag & SAM_FIRST)
	isFirst1 = FALSE ;
      if (sr->flag & SAM_LAST)
	isFirst2 = FALSE ;
    } 
  if (sr->flag & SAM_PAIRED_END)
    {
      array(sc->multiAli, sc->oneAli1, int)++ ;
      array(sc->multiAli, sc->oneAli2, int)++ ;
      array(sc->multiAli1, sc->oneAli1, int)++ ;
      array(sc->multiAli2, sc->oneAli2, int)++ ;
    }
  else
    array(sc->multiAli, sc->oneAli, int)++ ;

  sc->oneAli = sc->oneAli1 = sc->oneAli2 = 0 ;

}  /* samCountOneFragment */

/*************************************************************************************/

static int samCount (SAM *sam)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  ACEIN ai = 0 ;
  char cutter ;
  char fragName[1024], mateName[1024] ;
  int ln, nn = 0, nsr = 0 ;
  SR *sr = 0 ;
  Array srs = arrayHandleCreate (128, SR, h) ;
  BOOL isNew ;

  fprintf (stderr, "// sam scan starts: %s\n", timeShowNow()) ;

  if (sam->BAM)
    {
      if ( sam->gzi && sam->inFileName)
	ai = aceInCreateFromPipe (hprintf (h, "gunzip -c %s | samtools view -h | sort -k 1,1 ", sam->inFileName), "r", 0, h) ;
      else if (! sam->gzi && sam->inFileName)
	ai = aceInCreateFromPipe (hprintf (h, "samtools view -h   %s | sort -k 1,1 ", sam->inFileName), "r", 0, h) ;
      else if (sam->gzi)
	ai = aceInCreateFromPipe (hprintf (h, "gunzip -c - | samtools view -h -  | sort -k 1,1 ", sam->inFileName), "r", 0, h) ;
      else
	ai = aceInCreateFromPipe (hprintf (h, "samtools view -h - | sort -k 1,1 ", sam->inFileName), "r", 0, h) ;
    }
  else 
    {
      if (! sam->isSorted && sam->gzi)
	ai = aceInCreateFromPipe (hprintf (h, "gunzip -c %s | sort -k 1,1 ", sam->inFileName), "r", 0, h) ;
      else if (! sam->isSorted && ! sam->gzi && sam->inFileName)
	ai = aceInCreateFromPipe (hprintf (h, "sort -k 1,1  %s", sam->inFileName), "r", 0, h) ;
      else if (! sam->isSorted && ! sam->gzi && ! sam->inFileName)
	ai = aceInCreateFromPipe (hprintf (h, "sort -k 1,1 "), "r", 0, h) ;
      else
	ai = aceInCreate (sam->inFileName, sam->gzi, h) ;
    }
    
  fragName[0] = 0 ;
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      nn++ ;
      char *ccp ;

      /* fragment identifier */
      aceInStep (ai, '\t') ;  ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ;
      if (*ccp == '@')
	continue ;
      isNew = !sr || strcmp (fragName, ccp) ? TRUE : FALSE ;
      strcpy (fragName, ccp) ;
      if (sam->isFastc) 
	{
	  int k = 0 ;
	  char *cr = ccp + strlen (ccp) - 1 ;
	  while (*cr != '#' && cr > ccp)
	    cr-- ;
	  if (*cr == '#')
	    {
	      while (*++cr)
		{ k *= 10 ; k += (*cr - '0') ; }
	    }
	  if (k == 0)
	    k = 1 ;
	  sam->tagMult = k ;
	}
      if (isNew)
	{
	  if (nsr)
	    samCountOneFragment (sam, srs) ;
	  ac_free (h1) ;
	  h1 = ac_new_handle () ;
	  arrayMax (srs) = 0 ;
	  nsr = 0 ;
	}
      sr = arrayp (srs, nsr++, SR) ;
      memset (sr, 0, sizeof(SR)) ;
      /* flag */
      aceInStep (ai, '\t') ; if (! aceInInt(ai, &sr->flag)) continue ;

      /* reference identifier */
      aceInStep (ai, '\t') ;  ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ;
      if (*ccp != '*') /* jfm 10/29/2011 Don't lookup reference if it is an asterisk */
        {
          if (! dictFind (sam->targetDict, ccp, &sr->target))
	    {
	      char *cr = ccp + strlen(ccp) - 1 ;
	      if (*cr == '|')   /* makeBlastDB adds a parasit | at the end of the sequences */
		*cr = 0 ;
	      if (! dictFind (sam->targetDict, ccp, &sr->target))
		messcrash ("Sam file, Line %d, target %s not provided in the target fasta file", nn, ccp) ;
	    }
	}

      /* genome position 1-based */
      aceInStep (ai, '\t') ; if (! aceInInt(ai, &sr->pos)) continue ;
      /* quality */
      aceInStep (ai, '\t') ; if (! aceInInt(ai, &sr->mapq)) continue ;
      /* cigar */
      aceInStep (ai, '\t') ;  ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ;
      strcpy (sr->cigar, ccp) ;

      /* rnext */
      aceInStep (ai, '\t') ;  ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ;
      strcpy (mateName, ccp) ;
      /* pnext */
      aceInStep (ai, '\t') ; if (! aceInInt(ai, &sr->matePos)) continue ;
      /* tlen */
      aceInStep (ai, '\t') ; if (! aceInInt(ai, &sr->targetLn)) continue ;

      /* read dna */
      aceInStep (ai, '\t') ;  ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ;
      if (*ccp != '*')
	{
	  ln = strlen (ccp) ;
	  sr->dna = arrayHandleCreate (1024, char, h1) ;
	  
	  array (sr->dna, ln+1, char) = 0 ; /* make room */
	  memcpy (arrp(sr->dna, 0, char), ccp,ln) ;
	  arrayMax (sr->dna) = ln ;
	  dnaEncodeArray (sr->dna) ;   /* ensures uniform fonts and valid alphabet */
	}
    }
  
  if (nsr)
    samCountOneFragment (sam, srs) ;

  ac_free (h1) ;
  ac_free (h) ;

  fprintf (stderr, "// sam scan done, scanned %d lines: %s\n", nn, timeShowNow()) ;
  return nn ;
} /* samCount */

/*************************************************************************************/

static int samParseDo (SAM *sam, ACEIN ai, int run0)
{
  char *cp ;
  int i, tag, run, nn = 0 ;
  RC *rc = 0 ;
  TC *tc ;
  Array tags ;

  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      char cutter ;
      cp = aceInWordCut (ai, "\t", &cutter) ; if (! cp) continue ;
      dictAdd (sam->tagDict, cp, &tag) ;
      aceInStep (ai, '\t') ; cp = aceInWordCut (ai, "\t", &cutter) ; if (! cp) continue ;
      if (run0)
	run = run0 ;
      else
	dictAdd (sam->runDict, cp, &run) ;
      
      nn++ ;
      rc = arrayp (sam->runs, run, RC) ;
      if (! rc->tags)
	rc->tags = arrayHandleCreate (1024, TC, sam->h) ;
      tags = rc->tags ;
      tc = arrayp (tags, tag, TC) ; tc->tag = tag ;
      aceInStep (ai, '\t') ; aceInInt (ai, &tc->cols) ;
      if (tc->cols >= 8)
	messcrash ("Number of columns %d exceeds 8 at line %d of input file %s\n"
		   , tc->cols
		   , aceInStreamLine (ai) 
		   , aceInFileName (ai)
		   ) ;
      for (i = 0 ; i < tc->cols && i < 8 ; i++)
	{
	  double z = 0 ; 
	  aceInStep (ai, '\t') ; aceInDouble (ai, &z) ; tc->zz[i] += z ;
	}
    }  

  return nn ;
} /* samParseDo */

/*************************************************************************************/

static void samExport (SAM *sam)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (sam->outFileName, ".seqc.tsv", sam->gzo, h) ;
  int run, runMax = dictMax (sam->runDict) ;

  for (run = 1 ; run <= runMax ; run++)
    {
      RC *rc = arrp (sam->runs, run, RC) ;
      Array tags = rc ? rc->tags : 0 ;
      int tag, tagMax = tags ? arrayMax (tags) : 0 ;
      const char *runName = dictName (sam->runDict, run) ;

      for (tag = 1 ; tag <= tagMax ; tag++)
	{
	  TC *tc = arrp (tags, tag, TC) ;
	
	  if (tc->cols)
	    {
	      int j ;

	      aceOutf (ao, "%s\t%s\t%d"
		       , dictName (sam->tagDict, tag) 
		       , runName
		       , tc->cols
		       ) ;
	      for (j = 0 ; j < tc->cols ; j++)
		aceOutf (ao, "\t%ld", tc->zz[j]) ;
	      aceOut (ao, "\n") ;
	    }
	}
    }

  ac_free (h) ;
} /* samExport */

/*************************************************************************************/

static void samExportAce (SAM *sam)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (sam->outFileName, ".ace", sam->gzo, h) ;
  int pos, tag, run, runMax = dictMax (sam->runDict) ;
  char buf[256] ;
  int tagMax = dictMax (sam->tagDict) ;
  KEYSET tagPos = keySetHandleCreate (h) ;


  for (tag = 1 ; tag <= tagMax ; tag++)
    { 
      char *cq ;
      const char *tagName = dictName (sam->tagDict, tag) ;
      strncpy (buf, tagName, 200) ; 
      cq = strchr (buf, ':') ; pos = 0 ;
      if (cq) 
	{ *cq = 0 ; strncpy (buf, cq + 1, 200) ; sscanf (buf, "%d", &pos); }
      keySet (tagPos, tag) = pos ;
    }
  for (run = 1 ; run <= runMax ; run++)
    {
      RC *rc = arrp (sam->runs, run, RC) ;
      Array tags = rc ? rc->tags : 0 ;
      const char *runName = dictName (sam->runDict, run) ;

      tagMax = tags ? arrayMax (tags) : 0 ;
      aceOutf (ao, "Ali %s\n", runName) ;
      for (tag = 1 ; tag <= tagMax ; tag++)
	{
	  TC *tc = arrp (tags, tag, TC) ;
	  const char *tagName = dictName (sam->tagDict, tag) ;
	  char *cq ;

	  strncpy (buf, tagName, 200) ;
	  cq = strchr (buf, ':') ; pos = 0 ;
	  if (cq) 
	    { *cq = 0 ; pos = keySet (tagPos, tag) ; }

	  if (!strcmp (buf, "AliLn"))
	    aceOutf (ao, "Aligned_length %d %ld\n", pos, tc->zz[0]) ;

	  else if (!strcmp (buf, "Perfect_reads"))
	    aceOutf (ao, "Perfect_reads %ld\n",  tc->zz[0]) ;

	  else if (!strcmp (buf, "MismatchPerRead"))
	    aceOutf (ao, "Count_mismatch %d %ld\n", pos,  tc->zz[0]) ;

	  else if (!strcmp (buf, "Err_pos"))
	    aceOutf (ao, "Error_position %d %ld\n", pos,  tc->zz[0]) ;
	  else if (!strcmp (buf, "Ambiguous_pos"))
	    aceOutf (ao, "Ambiguous_position %d %ld\n", pos,  tc->zz[0]) ;

	  else if (!strcmp (buf, "Del"))
	    { int p = cq ? strlen (cq+1) : 0 ;
	      const char *mmm = "-----" ;
	      if (p) 
		aceOutf (ao, "Error_profile f1 %s%s %ld\n", mmm + 5 - p, cq+1,  tc->zz[0]) ;
	    }
	  else if (!strcmp (buf, "Ins"))
	    { int p = cq ? strlen (cq+1) : 0 ;
	      const char *mmm = "+++++" ;
	      if (p) 
		aceOutf (ao, "Error_profile f1 %s%s %ld\n", mmm + 5 - p, cq+1,  tc->zz[0]) ;
	    }
	  else if (!strcmp (buf, "Sub"))
	    { 
	      aceOutf (ao, "Error_profile f1 %c>%c %ld\n", ace_lower (buf[0]), ace_lower(buf[1]), tc->zz[0]) ;
	    }
	  else if (!strcmp (buf, "MultiAliZZZZZZZZZZ"))
	    aceOutf (ao, "Multiply-aligned %d %ld\n", pos,  tc->zz[0]) ;
	}
      aceOut (ao, "\n") ;
    }

  ac_free (h) ;
} /* samExportAce */

/*************************************************************************************/

static int samParse (SAM *sam)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0, run = 0 ;

  sam->tagDict = dictHandleCreate (512, sam->h) ;
  sam->runDict = dictHandleCreate (512, sam->h) ;
  sam->runs = arrayHandleCreate (64, RC, sam->h) ;

  if (sam->run)
    dictAdd (sam->runDict, sam->run, &run) ;

  if (! sam->inFileName)
    {
      ACEIN ai = aceInCreate (0, sam->gzi, h) ;
      nn = samParseDo (sam, ai, run) ;
    }
  else
    {
      char *buf = strnew (sam->inFileName, h) ;
      char *cp = buf, *cq ;

      while (cp)
	{
	  cq = strchr (cp, ',') ;
	  if (cq) *cq = 0 ;
	  if (*cp)
	    {
	      ACEIN ai = aceInCreate (cp, sam->gzi, h) ;
	      if (! ai)
		messcrash ("Cannot open input file %s", cp) ;
	      nn += samParseDo (sam, ai, run) ;
	      ac_free (ai) ;
	    }
	  cp = cq ? cp + 1 : 0 ;
	}
    }	

  ac_free (h) ;
  return nn ;
} /* samParse */

/*************************************************************************************/

static int samReport (SAM *sam, const char *run)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (sam->outFileName, ".seqc.tsv", sam->gzo, h) ;
  SC *sc = sam->sc ;
  int type, n, n1, n2 ;

  if (! run)
    run = "xxx" ;

  if (sam->fastaFileName)
    aceOutf (ao, "Reference_File:%s\t%s\t0\n"
	     , sam->fastaFileName
	     , run
	     ) ;
  aceOutf (ao, "Reads_in_file\t%s\t3\t%ld\t%ld\t%ld\n"
	   , run
	   , sc->reads_in_file, sc->reads_in_file1, sc->reads_in_file2
	   ) ;
  aceOutf (ao, "Alignments\t%s\t3\t%ld\t%ld\t%ld\n"
	   , run
	   , sc->ali, sc->ali1, sc->ali2
	   ) ;
  aceOutf (ao, "Unaligned_reads\t%s\t3\t%ld\t%ld\t%ld\n"
	   , run
	   , sc->reads_unaligned, sc->reads_unaligned1, sc->reads_unaligned2
	   ) ;
  aceOutf (ao, "Aligned_bases\t%s\t3\t%ld\t%ld\t%ld\n"
	   , run
	   , sc->bp_aligned, sc->bp_aligned1, sc->bp_aligned2
	   ) ;
  aceOutf (ao, "Unaligned_bases\t%s\t3\t%ld\t%ld\t%ld\n"
	   , run
	   , sc->bp_unaligned, sc->bp_unaligned1, sc->bp_unaligned2
	   ) ;
  aceOutf (ao, "Cumulated_read_length\t%s\t3\t%ld\t%ld\t%ld\n"
	   , run
	   , sc->read_bp, sc->read_bp1, sc->read_bp2
	   ) ;
  aceOutf (ao, "Reads_Mapped\t%s\t3\t%ld\t%ld\t%ld\n"
	   , run
	   , sc->reads_mapped, sc->reads_mapped1, sc->reads_mapped2
	   ) ;
  aceOutf (ao, "Reads_Mapped_per_strand\t%s\t6\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n"
	   , run
	   , sc->reads_mapped_on_plus_strand, sc->reads_mapped_on_minus_strand
	   , sc->reads_mapped_on_plus_strand1, sc->reads_mapped_on_minus_strand1
	   , sc->reads_mapped_on_plus_strand2, sc->reads_mapped_on_minus_strand2
	   ) ;
  aceOutf (ao, "Both_Reads_Mapped\t%s\t4\t%ld\t%ld\t%ld\t%ld\n"
	   , run
	   , sc->pairs_mapped, sc->pairs_compatible, 2*sc->pairs_compatible_topology, 2*sc->pairs_compatible_distance_below_1Mb
	   ) ;
  if (sc->perfect)
    aceOutf (ao, "Perfect_reads\t%s\t3\t%ld\t%ld\t%ld\n"
	     , run
	     , sc->perfect, sc->perfect1, sc->perfect2
	     ) ;
  if (sc->perfectClipped)
    aceOutf (ao, "Perfect_clipped\t%s\t3\t%ld\t%ld\t%ld\n"
	     , run
	     , sc->perfectClipped, sc->perfectClipped1, sc->perfectClipped2
	     ) ;

  for (type = 1 ; type < arrayMax (sc->aliLn) ; type++)
    {
      n = arr (sc->aliLn, type, int) ;
      n1 = type < arrayMax(sc->aliLn1) ? arr (sc->aliLn1, type, int) : 0 ;
      n2 = type < arrayMax(sc->aliLn2) ? arr (sc->aliLn2, type, int) : 0 ;

      if (n > 0)
	aceOutf (ao, "AliLn:%03d\t%s\t3\t%d\t%d\t%d\n"
		 , type
		 , run
		 , n, n1, n2
		 ) ;
    }

  for (type = 0 ; type < arrayMax (sc->errCount) ; type++)
    {
      n = arr (sc->errCount, type, int) ;
      n1 = type < arrayMax(sc->errCount1) ? arr (sc->errCount1, type, int) : 0 ;
      n2 = type < arrayMax(sc->errCount2) ? arr (sc->errCount2, type, int) : 0 ;

      if (n > 0)
	aceOutf (ao, "MismatchPerRead:%03d\t%s\t3\t%d\t%d\t%d\n"
		 , type
		 , run
		 , n, n1, n2
		 ) ;
    }

  if (sc->mismatches)
    aceOutf (ao, "Mismatches\t%s\t3\t%ld\t%ld\t%ld\n"
	     , run
	     , sc->mismatches, sc->mismatches1, sc->mismatches2
	     ) ;

  if (sc->substitutions)
    aceOutf (ao, "Substitutions\t%s\t3\t%ld\t%ld\t%ld\n"
	     , run
	     , sc->substitutions, sc->substitutions1, sc->substitutions2
	     ) ;

  if (sc->transitions)
    aceOutf (ao, "Transitions\t%s\t3\t%ld\t%ld\t%ld\n"
	     , run
	     , sc->transitions, sc->transitions1, sc->transitions2
	     ) ;

  if (sc->transversions)
    aceOutf (ao, "Transversions\t%s\t3\t%ld\t%ld\t%ld\n"
	     , run
	     , sc->transversions, sc->transversions1, sc->transversions2
	     ) ;


  if (sc->insertions)
    aceOutf (ao, "Insertions\t%s\t3\t%ld\t%ld\t%ld\n"
	     , run
	     , sc->insertions, sc->insertions1, sc->insertions2
	     ) ;
  if (sc->double_insertions)
    aceOutf (ao, "DoubletIns\t%s\t3\t%ld\t%ld\t%ld\n"
	     , run
	     , sc->double_insertions, sc->double_insertions1, sc->double_insertions2
	     ) ;
  if (sc->triple_insertions)
    aceOutf (ao, "TripletIns\t%s\t3\t%ld\t%ld\t%ld\n"
	     , run
	     , sc->triple_insertions, sc->triple_insertions1, sc->triple_insertions2
	     ) ;
  
  if ( sc->deletions)
    aceOutf (ao, "Deletions\t%s\t3\t%ld\t%ld\t%ld\n"
	     , run
	     , sc->deletions, sc->deletions1, sc->deletions2
	   ) ;
  if (sc->double_deletions)
    aceOutf (ao, "DoubletDel\t%s\t3\t%ld\t%ld\t%ld\n"
	     , run
	     , sc->double_deletions, sc->double_deletions1, sc->double_deletions2
	     ) ;
  if (sc->triple_deletions)
    aceOutf (ao, "TripletDel\t%s\t3\t%ld\t%ld\t%ld\n"
	     , run
	     , sc->triple_deletions, sc->triple_deletions1, sc->triple_deletions2
	     ) ;
  
  for (type = 1 ; type <= dictMax(sam->eeDict) ; type++)
    {
      n =  type < arrayMax(sc->errType) ? arr (sc->errType, type, int) : 0 ;
      n1 = type < arrayMax(sc->errType1) ? arr (sc->errType1, type, int) : 0 ;
      n2 = type < arrayMax(sc->errType2) ? arr (sc->errType2, type, int) : 0 ;

      if (n > 0)
	aceOutf (ao, "%s\t%s\t3\t%d\t%d\t%d\n"
		 , dictName (sam->eeDict, type)
		 , run
		 , n, n1, n2
		 ) ;
    }

  for (type = 1 ; type <= dictMax(sam->spikeDict) ; type++)
    {
      n =  type < arrayMax(sc->spikes) ? arr (sc->spikes, type, int) : 0 ;
      n1 = type < arrayMax(sc->spikes1) ? arr (sc->spikes1, type, int) : 0 ;
      n2 = type < arrayMax(sc->spikes2) ? arr (sc->spikes2, type, int) : 0 ;

      if (20*n > sc->mismatches)
	aceOutf (ao, "Spike:%s\t%s\t3\t%d\t%d\t%d\n"
		 , dictName (sam->spikeDict, type)
		 , run
		 , n, n1, n2
		 ) ;
    }

  for (type = 1 ; type < arrayMax (sc->multiAli) ; type++)
    {
      n = arr (sc->multiAli, type, int) ;
      n1 = type < arrayMax(sc->multiAli1) ? arr (sc->multiAli1, type, int) : 0 ;
      n2 = type < arrayMax(sc->multiAli2) ? arr (sc->multiAli2, type, int) : 0 ;

      if (n > 0)
	aceOutf (ao, "MultiAli:%03d\t%s\t3\t%d\t%d\t%d\n"
		 , type
		 , run
		 , n, n1, n2
		 ) ;
    }

  for (type = 0 ; type < arrayMax (sc->errPos) ; type++)
    {
      n = arr (sc->errPos, type, int) ;
      n1 = type < arrayMax(sc->errPos1) ? arr (sc->errPos1, type, int) : 0 ;
      n2 = type < arrayMax(sc->errPos1) ? arr (sc->errPos2, type, int) : 0 ;

      if (n > 0)
	aceOutf (ao, "Err_pos:%03d\t%s\t3\t%d\t%d\t%d\n"
		 , type
		 , run
		 , n, n1, n2
		 ) ;
    }

  for (type = 0 ; type < arrayMax (sc->nPos) ; type++)
    {
      n = arr (sc->nPos, type, int) ;
      n1 = type < arrayMax(sc->nPos1) ?  arr (sc->nPos1, type, int) : 0 ;
      n2 = type < arrayMax(sc->nPos2) ? arr (sc->nPos2, type, int) : 0 ;

      if (n > 0)
	aceOutf (ao, "Ambiguous_pos:%03d\t%s\t3\t%d\t%d\t%d\n"
		 , type
		 , run
		 , n, n1, n2
		 ) ;
    }

  ac_free (h) ;

  return 0 ;
} /* samReport */

/*************************************************************************************/
static int samIntronExport (SAM *sam, const char *run)
{
  int nn = 0, i, j, iMax ;
  Array aa = sam->sc->introns ;
  INTRON *iip, *jjp ;
  iMax = arrayMax (aa) ;
  if (! iMax)
    goto done ;

  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (sam->outFileName, ".introns.tsf", sam->gzo, h) ;

  arraySort (aa, samIntronOrder) ;
  for (i = 0, iip = arrp (aa, 0, INTRON) ; i < iMax ; i++, iip++)
    if (iip->support)
      for (j = i+1, jjp= iip+1 ; j < iMax ; j++, jjp++)
	{
	  if (iip->chrom != jjp->chrom ||
	      iip->a1 != jjp->a1 ||
	      iip->a2 != jjp->a2
	      )
	    break ;
	  iip->support += jjp->support ;
	  jjp->support = 0 ;
	}

  for (i = 0, iip = arrp (aa, 0, INTRON) ; i < iMax ; i++, iip++)
    if (iip->a1 && iip->support)
      aceOutf (ao, "%s__%d_%d\t%s\tiit\t%d\t%d\t%s\n"
	       , dictName (sam->targetDict, iip->chrom)
	       , iip->a1, iip->a2, sam->run, iip->support, iip->length, iip->feet
	       ) ;
  ac_free (h) ;
 done:
  return nn ;
} /* samIntronExport */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  fprintf  (stderr,
	    "// SAM.c:\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, October 2012, mieg@ncbi.nlm.nih.gov\n"
	    "//          Joe Meehan, NCTR\n"
	    "// Purpose\n"
	    "//   to analyze : Gather quality metrics form a SAM file, the reference fasta file is required\n"
	    "//   The results are exported in a tab delimited .seqc.tsv to be analysed by a companion python code seqc.py\n"
	    "// Input/output\n"
	    "//   -i input_file  [--gzi] : default stdin, if the file is called .gz or if -gzi, gunzip\n"
	    "//   -o output_file [--gzo] : default stdout, if -gzo, gzip the output and add a .gz endding\n"
	    "// Dependencies\n"
	    "//      To decompress the gz files, the program needs gunzip\n"
	    "//      To read a BAM file, the program needs samtools on the path\n"
	    "// Actions\n"
	    "// PHASE 1:  seqc --count --BAM --run <> --fasta <> -i <x.bam>  [-o <>]\n"
	    "//   --run run_name --fasta f\n"
	    "//   --count : count the alignments and estimate the quality metrics\n"
	    "//   -i file : input file, for example a .bam or .sam file\n"
	    "//   --fastc : reads are in fastc format, use the tag multiplicity\n"
	    "//   --fasta target.fasta : the target fasta file used to produce the bam file\n"
	    "//   --BAM | -B : the input is a BAM file, decompress it by calling samtools which must be in the path\n"
	    "//   --SAM | -S : the input is a SAM file\n"
	    "//   --isSorted : the input is already sorted on column 1, no need to sort it again. \n"
	    "//       ATTENTION BAM files are ususally sorted by chromosome, we need to sort per fragment (column 1)\n"
	    "//       so usually you do NOT want to say --isSorted\n"
	    "//\n"
	    "// PHASE 2: seqc --merge --run <>  -f <f.seqc.tsv,... [-o <>] --aceOut\n"
	    "//   --merge : all input counts will me merged\n"
	    "//   --run group_name : counts are reexported under that group name\n"
	    "//    -i f1,f2,f3.seqc.tsv : files exported from phase 1\n"
	    "//       another way is : cat *.seqc.tsv | seqc --merge --run_name x\n"
	    "//   --aceOut : export in .ace format, compatible with the MAGIC pipeline\n"
	    "//\n"
	    "//   --gzi : gunzip stdin, by default any input file named *.gz is gunzipped\n"
	    "//   --gzo : gzip the output files\n"
	    "//   --help | -h : export this on line help\n"
	    "// Caveat:\n"
	    "//   Lines starting with '#' in the input files are considered as comments and dropped out\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  SAM sam ;
  SC sc ; 
  AC_HANDLE h = 0 ;
 
  freeinit () ; 
  h = ac_new_handle () ;
  memset (&sam, 0, sizeof (SAM)) ;
  memset (&sc, 0, sizeof (SC)) ;
  sam.h = h ;
  sam.sc = &sc ;
  sam.targetDict = dictHandleCreate (1000, h) ;
  eeDictCreate (&sam) ;
  sc.errType = arrayHandleCreate (256, int, h) ;
  sc.errType1 = arrayHandleCreate (256, int, h) ;
  sc.errType2 = arrayHandleCreate (256, int, h) ;

  sc.errPos = arrayHandleCreate (256, int, h) ;
  sc.errPos1 = arrayHandleCreate (256, int, h) ;
  sc.errPos2 = arrayHandleCreate (256, int, h) ;

  sc.nPos = arrayHandleCreate (256, int, h) ;
  sc.nPos1 = arrayHandleCreate (256, int, h) ;
  sc.nPos2 = arrayHandleCreate (256, int, h) ;

  sc.errCount = arrayHandleCreate (256, int, h) ;
  sc.errCount1 = arrayHandleCreate (256, int, h) ;
  sc.errCount2 = arrayHandleCreate (256, int, h) ;

  sc.spikes = arrayHandleCreate (256, int, h) ;
  sc.spikes1 = arrayHandleCreate (256, int, h) ;
  sc.spikes2 = arrayHandleCreate (256, int, h) ;

  sc.introns = arrayHandleCreate (1024, INTRON, h) ;

  sc.aliLn = arrayHandleCreate (256, int, h) ;
  sc.aliLn1 = arrayHandleCreate (256, int, h) ;
  sc.aliLn2 = arrayHandleCreate (256, int, h) ;

  sc.multiAli = arrayHandleCreate (256, int, h) ;
  sc.multiAli1 = arrayHandleCreate (256, int, h) ;
  sc.multiAli2 = arrayHandleCreate (256, int, h) ;


  /* optional arguments */

  if (getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  getCmdLineOption (&argc, argv, "-i", &sam.inFileName) ;
  getCmdLineOption (&argc, argv, "-o", &sam.outFileName) ;
  sam.gzi = getCmdLineBool (&argc, argv, "-gzi") || getCmdLineBool (&argc, argv, "--gzi") ;
  sam.gzo = getCmdLineBool (&argc, argv, "-gzo") || getCmdLineBool (&argc, argv, "--gzo") ;

  sam.BAM = getCmdLineBool (&argc, argv, "-B") || getCmdLineBool (&argc, argv, "--BAM") ;
  sam.SAM = getCmdLineBool (&argc, argv, "-S") || getCmdLineBool (&argc, argv, "--SAM") ;

  sam.isSorted = getCmdLineBool (&argc, argv, "--isSorted") ;
  getCmdLineOption (&argc, argv, "-title", &sam.title) ;
  if (! getCmdLineOption (&argc, argv, "--run", &sam.run))
    getCmdLineOption (&argc, argv, "-r", &sam.run) ;

  sam.isFastc = getCmdLineBool (&argc, argv, "--fastc") ; /* use tag multiplicity */
  sam.count = getCmdLineBool (&argc, argv, "--count") ;
  sam.merge = getCmdLineBool (&argc, argv, "--merge") ;
  sam.aceOut = getCmdLineBool (&argc, argv, "--aceOut") ;
  sam.report = getCmdLineBool (&argc, argv, "--report") ;

  if (! getCmdLineOption (&argc, argv, "--fasta", &(sam.fastaFileName)))
    getCmdLineOption (&argc, argv, "-f", &(sam.fastaFileName)) ;

  fprintf (stderr, "// start: %s\n", timeShowNow()) ;
  if (argc > 1)
    usage (messprintf("sorry unknown argument %s", argv[1])) ;

  if (sam.BAM && sam.SAM)
    messcrash ("Sorry, you cannot specify at the same time --SAM and --BAM, please try --help") ;

  if (sam.count || sam.merge)
    {
      if (! sam.run)
	messcrash ("-count option requires --run run_name") ;
    }

  if (sam.fastaFileName)
    samParseFastaFile (&sam) ;

  if (sam.count)
    {
      samCount (&sam) ;
      samReport(&sam, sam.run) ;
      samIntronExport (&sam, sam.run) ;
    }
  else if (sam.merge)
    {
      samParse (&sam) ;
      samExport (&sam) ;
    }
  else if (sam.aceOut)
    {
      samParse (&sam) ;
      samExportAce (&sam) ;
    }
      
  fprintf (stderr, "// done: %s\n", timeShowNow()) ;
  {
    int mx ;
    messAllocMaxStatus (&mx) ;   fprintf (stderr, "// max memory %d Mb\n", mx) ;
  }
  
  ac_free (h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

