/*
 * sa.main.c : sortalign RNA aligner

 * Created April 18, 2025
 * In collaboration with Greg Boratyn, NCBI

 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality

 * The algorithm extract words from the union of the targets and from the reads creating 2 long tables
 * Both tables are sorted
 *  sorting the targets is a preprocessing, so the time does not count
 *  sorting the reads is done datablock by datablock
 * The two sorted tables are 'merged' which is local in memory
 *  exporting a sequential table.
 * That table is sorted per read
 *  alignments are extended using the magicblast jumper algorithm
 * Finally auxiliary tables, like introns, or wiggles, are exported on demand

 * The code is heavilly parallelized using the channel paradigm borrowed from the Go language
 */

/*
  #define ARRAY_CHECK
  #define MALLOC_CHECK
*/
/* Examples: try   sortalign -h */


#include "sa.h"

static int NN = 1 ;

static void showHitsDo (HIT *hit, long int iMax) ;


/**************************************************************/
/************** utilities *************************************/
/**************************************************************/

static void showHits (BigArray hits)
{
  long int ii, iMax = hits ? bigArrayMax (hits) : 0 ;
  int n = NTARGETREPEATBITS ;
  int mask = (1 << n) - 1 ;
  for (ii = 0 ; ii < iMax && ii < 2700 ; ii++)
    {
      HIT *hit = bigArrp (hits, ii, HIT) ;
      printf (".. %ld:  r=%u\t%u\t%u\t%u\tmult %d %s\n"
	      , ii, hit->read, hit->chrom, hit->a1, hit->x1 >> (n+3), hit->x1 & mask
	      , (hit->x1 >> n) & 0x7 ? "intron" : ""
	      ) ;
    }
  printf ("......... max %ld\n", iMax) ;
} /* showHits */

/**************************************************************/

static void showHitsDo (HIT *hit, long int iMax)
{
  if (hit)
    for (long int ii = 0 ; ii < iMax && ii < 30 ; ii++, hit++)
      {
	printf (".. r=%u\t%u\t%u\t%u\n"
		, hit->read, hit->chrom, hit->a1, hit->x1
		) ;
      }
  printf (".........\n") ;
} /* showHitsDo */

/**************************************************************/

static void showAli (Array aligns)
{
  if (aligns)
    {
      long int ii, iMax = arrayMax (aligns) ;
      
      for (ii = 0 ; ii < iMax && ii < 100 ; ii++)
	{
	  ALIGN *ali = arrp (aligns, ii, ALIGN) ;
	  if (ali->score)
	    printf (".. chain %d\tchainScore %d\tchainErr %d chainAli %d\tscore %d\tr=%u : %d:%d/%d\tchr=%u : %d:%d/%d\tnErr %d\t %d::p %d n %d\tdonor %d acc %d  a1-x1=%d\n"
		    , ali->chain, ali->chainScore
		    , ali->chainErr, ali->chainAli
		    , ali->score
		    , ali->read, ali->x0, ali->x1, ali->x2
		    , ali->chrom, ali->a0, ali->a1, ali->a2
		    , ali->nErr
		    , ali->id, ali->previous, ali->next
		    , ali->donor, ali->acceptor
		    , ali->a1 - ali->x1
		    ) ;
	}
      printf (".........\n") ;
    }
  return ;
} /* showAli */

/**************************************************************/
/**************************************************************/
/* cumulate the number of BB to be analyzed */
static int NTODO = 0 ; 
static void npCounter (const void *vp)
{
  const PP *pp = vp ;
  int nPuts = 0, nn = 0, n = 0 ;
  char tBuf[25] ;

  int nf = pp->nFiles ;
  
  while (n < nf)
    {
      channelGet (pp->npChan, &nPuts, int) ;
      nn += nPuts ;
      printf ("--- %s: npCounter Parsed %d/%d files\t%d blocs\tcumul %d blocks to be analysed\n",
	      timeBufShowNow (tBuf), ++n, nf, nPuts, nn) ;
    }
  printf ("--- %s: npCounter Processed %d files, closing pcChan at %d blocks\n", timeBufShowNow (tBuf), nf, nn) ;
  channelCloseSource  (pp->plChan) ;
  NTODO = nn ; /* this is global but there is no other place where we wrtie this variable */
  return ;
} /* npCounter */

/**************************************************************/
/**************************************************************/

static void loadRegulator (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  int nn = 0, nMax = pp->nBlocks ;
  
  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->plChan, &bb, BB))
    {
      nn++ ;
      while  (nn > nMax)
	{
	  int nd = channelCount (pp->doneChan) ;
	  nMax = pp->nBlocks - nd ;
	  if (nn <= nMax)
	    break ;
	  sleep (1) ;
	}
      channelPut (pp->lcChan, &bb, BB) ;
    }
  channelClose (pp->lcChan) ; /* there is a single instance of the load regulator */
  return ;
} /* load regulator */

/**************************************************************/
/**************************************************************/
/**************************************************************/

static void readParser (const void *vp)
{
  const PP *pp = vp ;
  RC rc ;
  
  while (channelGet (pp->fpChan, &rc, RC))
    saSequenceParse (pp, &rc, 0, 0, 0) ;
  return ;
} /* Readparser */

/**************************************************************/

#ifdef JUNK
include <sys/resource.h>

void* f(void* arg) {
    struct rusage usage_start, usage_end;
    // Get resource usage at start
    getrusage(RUSAGE_THREAD, &usage_start);
    
    // Simulate memory-intensive work
    int size = 1000000;
    char* buffer = malloc(size);
    for (int i = 0; i < size; i += 4096) { // Touch pages
        buffer[i] = 1;
    }
    
    // Get resource usage at end
    getrusage(RUSAGE_THREAD, &usage_end);
    
    // Calculate page faults
    long minor_faults = usage_end.ru_minflt - usage_start.ru_minflt;
    long major_faults = usage_end.ru_majflt - usage_start.ru_majflt;
    
    printf("Thread %ld: Minor faults: %ld, Major faults: %ld\n", 
           (long)arg, minor_faults, major_faults);
    free(buffer);
    return NULL;
}
#endif

/**************************************************************/

#ifndef YANN
static void codeWords (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  char tBuf[25] ;
  long int nnn = 0 ;
  clock_t  t1, t2 ;
  

  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->lcChan, &bb, BB))
    {
      long int nn= 0 ;

      if (pp->debug) printf ("+++ %s: Start code words\n", timeBufShowNow (tBuf)) ;

      t1 = clock () ;
      saSequenceParseGzBuffer (pp, &bb) ;
      if (pp->deduplicate)
	saSequenceDeduplicate (pp, &bb) ;
      saCodeSequenceSeeds (pp, &bb, pp->iStep, FALSE) ;
      t2 = clock () ;

      for (int i = 0 ; i < NN ; i++)
	nn += bigArrayMax (bb.cwsN[i]) ;
      nnn += nn ;
      saCpuStatRegister ("3.CodeWords", pp->agent, bb.cpuStats, t1, t2, nn) ;
      if (pp->debug) printf ("--- %s: Stop code words %ld\n", timeBufShowNow (tBuf), bigArrayMax (bb.cwsN[0])) ;

      channelPut (pp->csChan, &bb, BB) ;
    }

  int n = channelCount (pp->plChan) ;
  if (pp->debug) printf ("..... close csChan at %d,  coded %ld words\n", n, nnn) ;
  channelCloseSource (pp->csChan) ;

  return ;
} /* codeWords */
#endif

/**************************************************************/
/**************************************************************/

#ifndef YANN
static void sortWords (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  char tBuf[25] ;
  long int nnn = 0 ;
  clock_t  t1, t2 ;
  int k ;
  
  memset (&bb, 0, sizeof (BB)) ;
  while ((k = channelGet (pp->csChan, &bb, BB)))
    {
      if (k < 0)
	sleep (1) ;
      else
	{
	  long int nn = 0 ;
	  if (pp->debug) printf ("+++ %s: Start sort words\n", timeBufShowNow (tBuf)) ;
	  
	  t1 = clock () ;
	  for (int k = 0 ; k < NN ; k++)
	    if (bb.cwsN[k])
	      {
		bb.gpu += saSort (bb.cwsN[k], 1) ; /* cwOrder */
		nn += bigArrayMax (bb.cwsN[k]) ;
	      }
	  t2 = clock () ;
	
	  saCpuStatRegister ("4.SortWords", pp->agent, bb.cpuStats, t1, t2, nn) ;
	  if (pp->debug) printf ("--- %s: Stop sort words %ld\n", timeBufShowNow (tBuf), bigArrayMax (bb.cwsN[0])) ;
	  channelPut (pp->smChan, &bb, BB) ;
	}
    }

  int n = channelCount (pp->plChan) ;
  if (pp->debug) printf ("..... close smChan at %d,  sorted %ld words\n", n, nnn) ;
  
  channelCloseSource (pp->smChan) ;

  return ;
} /* sortWords */
#endif

/**************************************************************/

static long int  matchHitsDo (const PP *pp, BB *bbG, BB *bb)
{
  BOOL debug = FALSE ;
  BOOL useIntronSeeds = ! pp->ignoreIntronSeeds ;
  BigArray hitsArray = bigArrayHandleCreate(64, BigArray, bb->h);
  long int nn = 0, k = 0, kkk = 0 ;
  int nHA = 0, hMax = 100000 ;
  const long unsigned int mask26 = (1L << 26) - 1 ;
  BigArray hits = bigArrayHandleCreate(hMax, HIT, bb->h);
  bigArray (hitsArray, nHA++, BigArray) = hits ;
  const int seedLength = pp->seedLength ;
  const int intronBonus = 1 ;
  int absoluteMax = 0x1 << NTARGETREPEATBITS ;
  int absoluteX1Max = 0x1 << ( 31 - 3 - NTARGETREPEATBITS) ;
  int maxTargetRepeats = pp->maxTargetRepeats  ;
  
  for (int kk = 0; kk < NN ; kk++)
    {
      long int i = 0, iMax = bigArrayMax (bbG->cwsN[kk]);
      long int j = 0, jMax = bigArrayMax (bb->cwsN[kk]);


      if (!iMax || !jMax)
	continue ;
      
      unsigned int mask = step1 - 1 ;
      const CW *restrict cw = bigArrp(bbG->cwsN[kk], 0, CW) ;
      const CW *restrict rw = bigArrp(bb->cwsN[kk], 0, CW) ;
      const CW *restrict cw1;
      const CW *restrict cwMax = cw + iMax ;
      HIT *restrict hit;

#ifdef JUNK      
      if (kk == 6)
	{
	  const CW *restrict zw = bigArrp(bb->cwsN[kk], 0, CW) ;
	  for (long int i = 0 ; i < bigArrayMax (bb->cwsN[kk]) ; i++, zw++)
		 fprintf (stderr, "SEED\t%u\t%u\t%u\t%u\n"
			  , kk
			  , zw->nam
			  , zw->seed
			  , zw->pos
			  ) ;
	}
#endif      

      while  (i < iMax && j < jMax)
	{
	  if (0 && kk == 1 && rw->seed == 185667857)
	    printf("(rw->seed == 185667857)\n") ;
	  if ((i & mask) == 0)
	    {
	      if (rw->seed <= (unsigned int) cw->seed)
		{
		  cw++ ;
		  i++ ;
		  bb->skips0++;
		  continue;
		}
	      else if (rw->seed <= (unsigned int) cw->pos)
		{
		  cw += step1 ;
		  i += step1 ;
		  bb->skips1++;
		  continue;
		}
	      else if (rw->seed <= (unsigned int) cw->nam)
		{
		  cw += step2 ;
		  i += step2 ;
		  bb->skips2++;
		  continue;
		}
	      else if (rw->seed <= (unsigned int) cw->intron)
		{
		  cw += step3 ;
		  i += step3 ;
		  bb->skips3++;
		  continue;
		}
	      else
		{
		  cw += step4 ;
		  i += step4 ;
		  bb->skips4++;
		  continue;
		}
	    }
	  if (cw->seed < rw->seed)
	    /* { int di = ((i & 0xf) == 0 && (cw+16)->seed < rw->seed ? 16 : 1) ; i += di ; cw += di ; bb->skipsNotFound++ ;  } */
 	    { i++; cw++; bb->skipsNotFound++ ; } 
	  else if (cw->seed > rw->seed)
	    /* { int dj = ((j & 0xf) == 0 && (rw+16)->seed < cw->seed ? 16 : 1) ; j += dj ; rw += dj ; } */
	    { j++ ; rw++ ; }
	  else
	    {
	      long int a1, x1, i1 = i ;
	      bb->skipsFound++ ;
	      int nTargetRepeats = cw->intron ;
	      
	      if (0 &&   /* avoid, this kills the intron seeds */
		  nTargetRepeats > maxTargetRepeats)
		{ j++ ; rw++ ; continue ; }
	      if (nTargetRepeats >= absoluteMax)
		nTargetRepeats = absoluteMax - 1 ; /* we will report absoluteMax (i.e. 31) even is value is higher */

	      for (cw1 = cw ; cw1 < cwMax ; i1++, cw1++)
		{
		  if ((i1 & mask) == 0)
		    continue ;
		  if (cw1->seed != rw->seed)
		    break ;
		  /* success, non intron case */
		  if (((cw1->intron >> 31) & 0x1) == 0x0)
		    {
		      BOOL readUp = rw->nam & 0x1 ;
		      BOOL chromUp = cw1->nam & 0x1 ;
		      int nTargetRepeats = cw1->intron ;
		      
		      if (1 && nTargetRepeats > maxTargetRepeats)
			continue ; 
		      if (0 && useIntronSeeds) continue ;
		      nn++ ;
		      hit = bigArrayp (hits, k++, HIT) ;
		      hit->read = rw->nam >> 1 ;
		      a1 = cw1->pos ;
		      x1 = rw->pos ;
		      chromUp ^= readUp ; /* we want to be on strand plus of the read */
		      readUp = 0 ;
#ifdef JUNK      		    
		      if (1)
			{
			  fprintf (stderr, "MATCH\t%d\t%d\t%d\t%d\t%d\t%d\n"
				   , kk
				   , rw->nam
				   , rw->seed
				   , rw->pos
				   , cw->nam
				   , cw->pos
				   ) ;
			}
#endif		      
		      if (! chromUp)  /* plus strand of the genome */
			{
			  hit->a1 =
			    a1            /* position of the first base of the seed */
			    - x1          /* locate the virtual position of base 0 of the read */
			    + 1 ;         /* avoid zero */
			  hit->x1 = x1 << 3 ;  /* reserve 3 bits for the intron seeds */
			  hit->chrom = cw1->nam & 0xfffffffe ; /* to select plus strand, kill the last bit */
			  if (hit->x1 > absoluteX1Max)
			    messcrash ("read coordinate x1=%d too large, please edit the source code", x1) ;
			  hit->x1 = ( hit->x1 << NTARGETREPEATBITS) | nTargetRepeats ;  /* all intron seeds are valuable */			
			}
		      else    /* minus strand of the genome */
			{
			  hit->a1 =
			    a1            /* position of the first base of the seed */
			    + (seedLength - 1)  /* go the last base of the seed */
			    + x1   /* locate the virtual position of base 0 of the read */
			    + 1 ;  /* avoid zero */
			  hit->x1 = x1 << 3 ;  /* reserve 3 bits for the intron seeds */
			  hit->chrom = cw1->nam | 0x1 ; /* to select minus strand, set the last bit */
			  if (hit->x1 > absoluteX1Max)
			    messcrash ("read coordinate x1=%d too large, please edit the source code", x1) ;
			  hit->x1 = ( hit->x1 << NTARGETREPEATBITS) | nTargetRepeats ;  /* all intron seeds are valuable */			
			}
		    }
		  else  if (useIntronSeeds)  /* INTRON */
		    {
		      BOOL readUp = rw->nam & 0x1 ;
		      BOOL chromUp = cw1->nam & 0x1 ; 

		      unsigned int z = cw1->intron ;
		      unsigned int isIntronDown = (z >> 28) & 0x4 ;
		      int da1 =  z & 0xf ; /* nb of letters in first exon : 4....11 */
		      int da  =  ((z >> 4) & mask26) ;  /* intron length < 32Mb */

		      chromUp ^= readUp ; /* we want to be on strand plus of the read */
		      readUp = 0 ;
		      hit = 0 ;
		      
		      if (! chromUp)  /* plus strand on the read and on the genome */
			{
			  a1 = cw1->pos ;       /* first base of intron in the genome, in bio coordinates */
		          x1 = rw->pos + da1 ;  /* matching base on the read */
			  
			  /* Create a hit to the last two bases of the donor exon (x1-2 / a1-2) */
			  nn++ ;
			  hit = bigArrayp (hits, k++, HIT) ;
			  hit->read = rw->nam >> 1 ;
			  hit->chrom = cw1->nam & 0xfffffffe ; /* to select plus strand, kill the last bit */
			  if (0 && isIntronDown == 0) hit->chrom |= 0x1 ; /* to cluster on the correct strand of the chroms */
			  hit->x1 =
			    ((x1 - 2) << 3)   /* reserve 3 bits for the intron seeds */
			    | isIntronDown    /* bit 3 gives the intron strand */
			    | 0x1             /* donor site */
			    ;
			  hit->a1 =
			    (a1 - 2)            /* position of the first base of the seed */
			    - (x1 - 2)          /* locate the virtual position of base 0 of the read */
			    - intronBonus       /* prefer intron match to exon match */
			    + 1 ;               /* avoid zero */
			  if (hit->x1 > absoluteX1Max)
			    messcrash ("read coordinate x1=%d too large, please edit the source code", x1) ;
			  hit->x1 = ( hit->x1 << NTARGETREPEATBITS) | 0x1 ;  /* all intron seeds are valuable */
			  
			  /* Create a hit to the first two bases of the acceptor exon (x1/ a1 + da) */
			  nn++ ;
			  hit = bigArrayp (hits, k++, HIT) ;
			  hit->read = rw->nam >> 1 ;
			  hit->chrom = cw1->nam & 0xfffffffe ; /* to select plus strand, kill the last bit */
			  if (0 && isIntronDown == 0) hit->chrom |= 0x1 ; /* to cluster on the correct strand of the chroms */
			  hit->x1 =
			    ((x1) << 3)   /* reserve 3 bits for the intron seeds */
			    | isIntronDown    /* bit 3 gives the intron strand */
			    | 0x2             /* acceptor site */
			    ;
			  hit->a1 =
			    (a1 + da)            /* position of the first base of the seed */
			    - (x1)          /* locate the virtual position of base 0 of the read */
			    - intronBonus       /* prefer intron match to exon match */
			    + 1 ;               /* avoid zero */
			  if (hit->x1 > absoluteX1Max)
			    messcrash ("read coordinate x1=%d too large, please edit the source code", x1) ;
			  hit->x1 = ( hit->x1 << NTARGETREPEATBITS) | 0x1 ;  /* all intron seeds are valuable */
			}

		     else  /* plus strand on the read and minus strand  on the genome */
			{
			  a1 = cw1->pos ;       /* first base of intron in the genome, in bio coordinates */
		          x1 = rw->pos + (seedLength - da1) - 1 ;  /* matching base on the read */
			  
			  /* Create a hit to the first two bases of the acceptor exon (x1+1,x1+2 / a1-1,a1-2) */
			  nn++ ;
			  hit = bigArrayp (hits, k++, HIT) ;
			  hit->read = rw->nam >> 1 ;
			  hit->chrom = cw1->nam | 0x1 ; /* to select minus strand, set the last bit */
			  if (0 && isIntronDown == 0) hit->chrom ^= 0x1 ; /* to cluster on the correct strand of the chroms */
			  hit->x1 =
			    ((x1 + 1) << 3)   /* reserve 3 bits for the intron seeds */
			    | isIntronDown    /* bit 3 gives th intron strand */
			    | 0x2             /* acceptor site */
			    ;
			  hit->a1 =
			    (a1 - 1)       /* position of the first base of the seed */
			    + (x1 + 1)          /* locate the virtual position of base 0 of the read */
			    - intronBonus       /* prefer intron match to exon match */
			    + 1 ;               /* avoid zero */
			  hit->x1 = ( hit->x1 << NTARGETREPEATBITS) | 0x1 ;  /* all intron seeds are valuable */
			  
			  /* Create a hit to the last two bases of the donor exon (x1-1,x1/ a1 + da+1,a1+da) */
			  nn++ ;
			  hit = bigArrayp (hits, k++, HIT) ;
			  hit->read = rw->nam >> 1 ;
			  hit->chrom = cw1->nam | 0x1 ; /* to select minus strand, set the last bit */
			  if (0 && isIntronDown == 0) hit->chrom ^= 0x1 ; /* to cluster on the correct strand of the chroms */
			  hit->x1 =
			    ((x1 - 1) << 3)   /* reserve 3 bits for the intron seeds */
			    | isIntronDown    /* bit 3 gives th intron strand */
			    | 0x1             /* donor site */
			    ;
			  hit->a1 =
			    (a1 + da + 1)            /* position of the first base of the seed */
			    + (x1 - 1)          /* locate the virtual position of base 0 of the read */
			    - intronBonus       /* prefer intron match to exon match */
			    + 1 ;               /* avoid zero */
			  hit->x1 = ( hit->x1 << NTARGETREPEATBITS) | 0x1 ;  /* all intron seeds are valuable */
			}
		    }
		  else
		    continue ;
		  
		  hit->chrom ^= (hit->read & 0x1) ; /* flip chrom for the second read of a pair */
		  if (0 && rw->seed == 168430082)
		    printf("(rw->seed == 168430082)\n") ;
		
		  if (k >= hMax)
		    {
		      bigArrayMax (hits) = k  ;
		      hits = bigArrayHandleCreate(hMax, HIT, bb->h);
		      bigArray (hitsArray, nHA++, BigArray) = hits ;
		      kkk += k ;
		      k = 0 ;
		    }
		}
	      j++ ; rw++ ;
	    }
	}
      bigArrayDestroy (bb->cwsN[kk]) ; /* the read words are no longer needed */
      bigArrayMax (hits) = k ;
      kkk += k ;
    }
  bb->hits = hitsArray ;

  if (debug)
    fprintf (stderr, "..MatchHitsDo found %ld matches\n", kkk) ;

  return nn ;
}  /* matchHitsDo */

/**************************************************************/
/**************************************************************/

#ifndef YANN
static void matchHits (const void *vp)
{
  const PP *pp = vp ;
  BB bb ;
  BB bbG = pp->bbG;
  char tBuf[25] ;
  long int nnn = 0 ;
  
  clock_t  t1, t2 ;
	    
  memset (&bb, 0, sizeof (BB)) ;
  /* grab and match a block of reads */
  while (channelGet (pp->smChan, &bb, BB))
    {
      if (bb.length)
	{
	  if (pp->debug) printf ("+++ %s: Start match %ld bases against %ld target bases\n", timeBufShowNow (tBuf), bb.length, bbG.length) ;


	  t1 = clock () ;
	  long int nn = matchHitsDo (pp, &bbG, &bb) ;
	  if (pp->debug) printf ("--- %s: Stop match hits constructed %ld arrays\n", timeBufShowNow (tBuf), bigArrayMax (bb.hits)) ;
	  t2 = clock () ;

	  nnn += nn ;
	  saCpuStatRegister ("5.MatchHits", pp->agent, bb.cpuStats, t1, t2, nn) ;
	}
      channelPut (pp->moChan, &bb, BB) ;
    }

  int n = channelCount (pp->plChan) ;
  if (pp->debug) printf ("..... close moChan at %d,  found %ld hits\n", n, nnn) ;
  channelCloseSource (pp->moChan) ;

  return ;
} /* matchHits */
#endif

/**************************************************************/
/**************************************************************/
/**************************************************************/

static void sortHitsFuse (const PP *pp, BB *bb)
{
  BigArray aa = bb->hits ;
  long int iMax = aa ? bigArrayMax (aa) : 0 ;
  char tBuf[25] ;
  
  if (pp->debug) printf ("+++ %s: Start sort hits merging %ld arrays\n", timeBufShowNow (tBuf), iMax) ;


  if (iMax == 0)
    bb->hits = 0 ;
  else if (iMax == 1)
    bb->hits = bigArray (aa, 0, BigArray) ;
  else
    {
      BigArray a ;
      long int n = 0 ;
      HIT *up ;
      for (long int i = 0 ; i < iMax ; i++)
	{
	  a = bigArray (aa, i, BigArray) ;
	  n += bigArrayMax (a) ;
	}
      bb->hits = bigArrayHandleCreate (n, HIT, bb->h) ;
      up = bigArrayp (bb->hits, n - 1, HIT) ; /* make room */
      up = bigArrayp (bb->hits, 0, HIT) ;
      for (int i = 0 ; i < iMax ; i++)
	{
	  a = bigArray (aa, i, BigArray) ;
	  n = bigArrayMax (a) ;
	  if (n > 0)
	    {
	      memcpy (up, bigArrp (a, 0, HIT), n * sizeof (HIT)) ;
	      up += n ;
	    }
	  bigArrayDestroy (a) ;
	}
      bigArrayDestroy (aa) ;
    }
  return ;
} /* sortHitsFuse */

/**************************************************************/
/* sort the hits */
#ifndef YANN
static void sortHits (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  char tBuf[25] ;
  long int nnn = 0 ;
  clock_t  t1, t2 ;

  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->moChan, &bb, BB))
    {
      if (pp->align && bb.hits)
	{
	  t1 = clock () ;
	  sortHitsFuse (pp, &bb) ;
	  if (bb.hits)
	    {
	      bb.gpu += saSort (bb.hits, 3) ; /* hitPairOrder */
	      t2 = clock () ;
	      
	      long int nn = bigArrayMax (bb.hits) ;
	      nnn += nn ;
	      saCpuStatRegister ("6.SortHits", pp->agent, bb.cpuStats, t1, t2, nn) ;
	      if (pp->debug) printf ("--- %s: Stop sort hits %ld\n", timeBufShowNow (tBuf), bigArrayMax (bb.hits)) ;
	    }
	}
      channelPut (pp->oaChan, &bb, BB) ;
    }

  if (1)
    {
      memset (&bb, 0, sizeof (BB)) ;
      bb.isGenome = TRUE ;
      channelPut (pp->doneChan, &bb, BB) ; /* destroy bbG.cws, all Matches are already computed */ 
    }
  channelCloseSource (pp->oaChan) ;
  channelCloseSource (pp->doneChan) ;
  
  return ;
} /* sortHits */
#endif

/**************************************************************/
/**************************************************************/
#include <stdint.h>

/**************************************************************/

static void runErrorsCumulate (int run, Array runErrors, Array bb)
{
  Array aa = array (runErrors, run, Array) ;
  int iMax = arrayMax (bb) ;
  if (aa && iMax)
    for (int i = 0 ; i < iMax ; i++)
      array (aa, i, int) += array (bb, i, int) ;
} /* runErrorsCumulate */

/**************************************************************/

static void reportRunErrors (const PP *pp, Array runStats, Array runErrors)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (pp->outFileName, ".runErrors.tsf", 0, h) ;
  int type, run, runMax = arrayMax (runErrors) ;

  for (run = 1 ; run < runMax ; run++)
    {
      char buf[8] ;
      const char *AGCT = "AGCT" ;
      const char *rNam = dictName (pp->runDict, run) ;      
      RunSTAT *rs = arrayp (runStats, run, RunSTAT) ;
      Array errors = array (runErrors, run, Array) ;

      aceOutDate (ao, "##", "Mismatches") ;
      aceOutf (ao, "%s\tAny\ti\t%d\n"
	       , rNam, rs->nErr
	       ) ;

      if (! errors) continue ;
      /* substitutions */
      for (int longBase = 0 ; longBase < 4 ; longBase++)
	{
	  buf[0] = AGCT[longBase] ;
	  buf[1] = '>' ;
	  for (int shortBase = 0 ; shortBase < 4 ; shortBase++)
	    {
	      type = (longBase << 2) | shortBase ;
	      buf[2] = AGCT[shortBase] ;
	      buf[3] = 0 ;

	      aceOutf (ao, "%s\t%s\ti\t%d\n"
		       , rNam, buf, array (errors, type, int)
		       ) ;
	    }
	}
      /* insertions */
      memcpy (buf, "Ins", 3) ;
      for (int shortBase = 0 ; shortBase < 4 ; shortBase++)
	{
	  type = 0x10 + shortBase ;
	  buf[3] = AGCT[shortBase] ;
	  buf[4] = 0 ;
	  aceOutf (ao, "%s\t%s\ti\t%d\n"
		   , rNam, buf, array (errors, type, int)
		   ) ;
	}
      /* deletions */
      memcpy (buf, "Del", 3) ;
      for (int longBase = 0 ; longBase < 4 ; longBase++)
	{
	  type = 0x20  | (longBase << 2) ;
	  buf[3] = AGCT[longBase] ;
	  buf[4] = 0 ;
	  aceOutf (ao, "%s\t%s\ti\t%d\n"
		   , rNam, buf, array (errors, type, int)
		   ) ;
	}
    }
  ac_free (h) ;
} /* runErrorsCumulate */

/**************************************************************/
/********************************************************************/
/********************************************************************/

static void sortAlignTableCaption (const PP *pp, ACEOUT ao)
{
  char tBuf[25] ;
  aceOutf (ao, "## %s Magic aligner\n", timeBufShowNow(tBuf)) ;
  if (1)
    aceOutf (ao, "## Author: Danielle et Jean Thierry-Mieg, Greg Boratyn, NCBI, mieg@ncbi.nlm.nih.gov\n") ;
  if (0)
    aceOutf (ao, "## Latest version and documentation are  available from https://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/Software\n") ;
  aceOutf (ao, "## Sequences are 1-based\n") ;
  if(0)
    {
      aceOutf (ao, "## 5\'All DNA sequences are read on the strand going away from the edges of the alignment\n") ;
      aceOutf (ao, "## So for donors and acceptors the match part (identical on read and target) and the overhang part are read on the opposite strands,\n"
	       "## The advantage of this unusual convention is that the beginning of all these words is independent of the length of the read, \n"
	       "## independent of the strand of the read or of the target, and the overhang at the donor site is identical to the \n"
	       "## match at the acceptor site and vice versa, so the donor acceptor pairs are easy to recognize.\n"
	       ) ;
    }
  aceOutf (ao, "##Read\tAlignment score\tRead multiplicity") ;
  aceOutf (ao, "\tRead length to be aligned, i.e. removing eventual adaptors, leading Ts or trailing As not matching the target") ;
  aceOutf (ao, "\tRead length aligned") ;
  aceOutf (ao, "\tfrom base") ;
  aceOutf (ao, "\tto base, coordinates of the alignment in the read") ;
  
  
  aceOutf (ao, "\tTarget_class (A_mito ... Z_genome) defines a Bayesian hierarchy of most desirable targets, used by some of the post-processing steps. Theses names are prioritized alphabetically") ;
  aceOutf (ao, "\tGene, as given in the target fasta file as : >target_name|Gene|gene_name, to avoid confusions, RefSeq gene names are prefixed by X__") ;
  aceOutf (ao, "\tUniqueness, number of genes or genomic sites on which the read aligns at best score, if the strand of the best alignements varies, the number is negative, hence -2 likely indicates a single alignment on 2 genes antisense to each other") ;
  aceOutf (ao, "\tTarget, the identifier given in the -t fasta file, for example a chromosome or a trancscritpt") ;
  
  aceOutf (ao, "\tfrom base") ;
  aceOutf (ao, "\tto base, coordinates of the alignment in the target") ;
  
  /*
    if (pp->solid)
    aceOutf (ao, "\tNumber of SOLiD corrected bases") ;
    else
  */
  aceOutf (ao, "\tNumber of N bases") ;
  aceOutf (ao, "\tNumber of mismatches") ;
  
  
  
  aceOutf (ao, "\tMissmatches: read position:type,... single base variation  (base in target > base in read); single base insertion (+ base inserted in read); single base deletion (- target base missing from read), given in the coordinates and the orientation of the read") ;
  aceOutf (ao, "\tMissmatches: target position:type,... single base variation  (base in target > base in read); single base insertion (+ base inserted in read); single base deletion (- target base missing from read), given in the coordinates and the orientation of the target") ;
  
  aceOutf (ao, "\t5' overhang read backwards, i.e complementary to the unaligned 5' part of the read, max %d bp", pp->OVLN) ;
  aceOutf (ao, "\t3' overhang read forward, i.e. identical to the unaligned 3' part of the read, max %d bp", pp->OVLN) ;
  
  aceOutf (ao, "\tChain\tFragment length as measured on target in case of paired end sequencing or 0 if not available at this stage of the program") ;
  aceOutf (ao, "\n") ;
  
  aceOutf (ao, "#Run/Read\tScore\tRead multiplicity\tLength to align\tAligned length\tfrom base\tto base in read") ;
  aceOutf (ao, "\tTarget class\tGene (X__name for RefSeq)\tTarget multiplicity\tTarget name (_name for RefSeq)\tfrom base\tto base on target") ;
  
  aceOutf (ao, "\tNumber of %s bases\tNumber of mismatches\tPosition:type in read\tPosition:type in target\t5'overhang (reverse-complemented)\t3'Overhang\tchain\\tPair length", pp->solid ? "corrected" : "N") ;
  
  aceOutf (ao, "\n") ;
} /* sortAlignTableCaption */

/**************************************************************/

static void exportDo (const PP *pp, BB *bb)
{
  ALIGN *ap = bigArrp (bb->aligns, 0, ALIGN) ;
  long int ii, aMax = bigArrayMax (bb->aligns) ;
  DICT *dict = bb->dict ;
  DICT *errDict = bb->errDict ;
  DICT *dictG = pp->bbG.dict ;
  BOOL pairedEnd = bb->runStat.p.nPairs > 0 ;
  const char *run = dictName (pp->runDict, bb->run) ;
  AC_HANDLE h = ac_new_handle () ;
  char *runNam = hprintf (h, "hits/%s.%d.hits", run, bb->lane) ;
  ACEOUT ao = aceOutCreate (pp->outFileName, runNam, pp->gzo, h) ;
  aceOutDate (ao, "###", "sortaling hits") ;
  
  sortAlignTableCaption (pp, ao) ;

  for (ii = 0 ; ii < aMax ; ii++, ap++)
    {
      int read = ap->read ;
      int chrom = ap->chrom ;
      int x1 = ap->x1 ;
      int x2 = ap->x2 ;
      int a1 = ap->a1 ;
      int a2 = ap->a2 ;
      int nerr = ap->nErr ;
      int nN = ap->nN ;

      if (read)
	{
	  aceOutf (ao, "%s/%s%s", run, dictName (dict, read >> 1), (pairedEnd ? ((read & 0x1 )? "<" : ">") : "")) ; 
	  aceOutf (ao, "\t%d", ap->chainScore) ;
	  aceOutf (ao, "\t%d", 1) ; /* ap->multiplicity */
	  aceOutf (ao, "\t%d", ap->readLength) ;
	  aceOutf (ao, "\t%d\t%d\t%d", ap->chainAli, x1, x2) ;
	  
	  aceOutf (ao, "\t%c\t-", ap->targetClass) ;
	  aceOutf (ao, "\t%d", ap->nChains) ; /* ap->targetMultiplicity */
	  aceOutf (ao, "\t%s", dictName (dictG, chrom >> 1)) ; 
	  aceOutf (ao, "\t%d\t%d", a1, a2) ;
	  aceOutf (ao, "\t%d\t%d", nN, nerr) ;
	  aceOutf (ao, "\t%s", ap->errShort ? dictName (errDict, ap->errShort) : "-") ;
	  aceOutf (ao, "\t%s", ap->errLong ? dictName (errDict, ap->errLong) : "-") ;
	  aceOut  (ao, "\t-\t-") ;  /* prefix suffix in genome */
	  aceOutf (ao, "\t%s", ap->leftOverhang ? dictName (dict, ap->leftOverhang) : "-") ;
	  aceOutf (ao, "\t%s", ap->rightOverhang ? dictName (dict, ap->rightOverhang) : "-") ;
	  aceOutf (ao, "\tchain %d 1", ap->chain) ;
	  aceOutf (ao, "\t%d", ap->chainScore) ;
	  aceOut  (ao, "\n") ;
	}
    }

  ac_free (h) ;
  return ;
} /* exportDo */

/**************************************************************/
/**************************************************************/

static void export (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  clock_t  t1, t2, dt = 50 * CLOCKS_PER_SEC ;
  
  char tBuf[25] ;
  int n, nn = 0 ;
  float nG = 0 ;

  int runMax = dictMax (pp->runDict) ;
  ACEOUT aos[runMax + 1]  ;
  ACEOUT aoes[runMax + 1]  ;
  AC_HANDLE h = ac_new_handle () ;
  
  
  memset (aos, 0, sizeof (aos)) ;
  memset (aoes, 0, sizeof (aoes)) ;

  memset (&bb, 0, sizeof (BB)) ;
  if (pp->sam || pp->bam)
    {
      for (int run = 1 ; run <= runMax ; run++)
	{
	  bb.run = run ;
	  aos[bb.run] = saSamCreateFile (pp, &bb, FALSE, pp->bamHandle) ;
	  aoes[bb.run] = saSamCreateFile (pp, &bb, TRUE, pp->bamHandle) ;
	}
    }

  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->aeChan, &bb, BB))
    {
      t1 = clock () ;
      if (bb.aligns && bigArrayMax (bb.aligns))
	{
	  if (pp->sam || pp->bam || pp->hitsFormat)
	    bigArraySort (bb.aligns, saAlignOrder) ;
	  if (pp->hitsFormat)
	    exportDo (pp, &bb) ;
	  if (pp->sam || pp->bam)
	    {
	      ACEOUT ao = aos[bb.run] ;
	      ACEOUT aoe = aoes[bb.run] ;
	      if (! ao)
		messcrash ("canot open the sam/bam file") ;

	      saSamExport (ao, aoe, pp, &bb) ;
	    }
	  t2 = clock () ;
	  saCpuStatRegister ("8.Export_ali", pp->agent, bb.cpuStats, t1, t2, bb.aligns ? bigArrayMax (bb.aligns) : 0) ;
	  n = channelCount (pp->plChan) ;
	  nG += bb.length / 1.0e+9 ;
	  ++nn ;
	  if (t2 - t1 > dt) 
	    {
	      printf ("--- %s: Exported %d/%d blocks, run %s, %lu reads, %lu ali, %.2f Gbases\n", timeBufShowNow (tBuf), nn, n, dictName (pp->runDict, bb.run), bb.nSeqs, bb.nAli, nG) ;
	      t1 = clock () ;
	    }
	}

      channelPut (pp->doneChan, &bb, BB) ;
    }

  n = channelCount (pp->plChan) ;
  if (0) printf ("--- %s: Export closes doneChan at %d\n", timeBufShowNow (tBuf), n) ;
  channelCloseSource (pp->doneChan) ;
  

  ac_free (h) ;
  return ;
} /* export */

/**************************************************************/

#ifdef YANN
static void wholeWork (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  BB bbG = pp->bbG;
  char tBuf[25] ;
  long int nnn = 0 ;
  clock_t  t1, t2, t01, t02 ;

  t01 = clock () ;
  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->lcChan, &bb, BB))
    {
      long int nn = 0 ;

      t1 = clock () ;
      /* code words */
      saSequenceParseGzBuffer (pp, &bb) ;
      saCodeSequenceSeeds (pp, &bb, pp->iStep, FALSE) ;

      if (1 || pp->debug) printf ("+++ %s: Start wholeWork agent %d, lane %d, %ld bases against %ld target bases\n", timeBufShowNow (tBuf), pp->agent, bb.lane, bb.length, bbG.length) ;

      /* sort words */
      for (int k = 0 ; k < NN ; k++)
	if (bb.cwsN[k])
	  bb.gpu += saSort (bb.cwsN[k], 1) ; /* cwOrder */

      
      /* match hits */
      if (bb.length)
	nn = matchHitsDo (pp, &bbG, &bb) ;
      nnn += nn ;

      /* sorthits */
      if (pp->align && bb.hits)
	{
	  sortHitsFuse (pp, &bb) ;
	  bb.gpu += saSort (bb.hits, 3) ; /* hitPairOrder */
	  saAlignDo (pp, &bb) ;
	}

#ifdef JUNK      
      /* export */
      if (bb.aligns && bigArrayMax (bb.aligns))
	{
	  bigArraySort (bb.aligns, saAlignOrder) ;
	  if (! pp->sam)
	    { if (1) exportDo (pp, &bb) ; }
	  else
	    {
	      AC_HANDLE h = ac_new_handle () ;
	      ACEOUT ao = bb.rc.aoSam  ;
	      char *VERSION = "0.1.1" ;
	      DICT *dictG = pp->bbG.dict ;
	      ao = aceOutCreate (pp->outFileName, hprintf (h, ".%s.%d.sam", dictName (pp->runDict, bb.run), bb.lane), pp->gzo, h) ;
	      aceOutf (ao, "@HD VN:1.5\tSO:queryname\n") ;
	      aceOutf (ao, "@PG ID:1\tPN:Magic\tVN:%s\n", VERSION) ;
	      
	      for (int chrom = 1 ; chrom <= dictMax (dictG) ; chrom++)
		{
		  Array dna = arr (pp->bbG.dnas, chrom, Array) ;
		  int ln = dna ? arrayMax (dna) : 0 ;
		  aceOutf (ao, "@SQ\tSN:%s\tLN:%d\n", dictName (dictG, chrom), ln) ;
		}
	      /* aceOutf (ao, "\tCL:%s", commandBuf) ; */
	      exportSamDo (ao, pp, &bb) ;
	      ac_free (h) ;
	    }
	}
#endif
      
      t2 = clock () ;
      saCpuStatRegister ("5.WholeWork", pp->agent, bb.cpuStats, t1, t2, nn) ;
      channelPut (pp->aeChan, &bb, BB) ;
      t02 = clock () ;
      saCpuStatRegister ("5.WholeWorkE", pp->agent, bb.cpuStats, t01, t02, nn) ;
      t01 = t02 ;
    }

  channelCloseSource (pp->aeChan) ;
  return ;
} /* wholeWork */
#endif

/*************************************************************************************/

static void reportRunStats (PP *pp, Array runStats)
{
  AC_HANDLE h = ac_new_handle () ;
  Array aa = runStats ;
  RunSTAT *s0 = arrayp (aa, 0, RunSTAT) ;

  printf ("\n####### Global Statistics, see details in file %srunStats.tsf", pp->outFileName ? pp->outFileName : "") ;
  if (pp->runName) printf ("\n#:Run\t%s", pp->runName) ;
  if (pp->method) printf ("\n#:Method\t%s", pp->method) ;
  
  printf ("\n#:Pairs\t%ld", s0->p.nPairs) ;
  printf ("\n#:Pairs_aligned\t%ld\t%.3f%%", s0->nPairsAligned, 100.0 * s0->nPairsAligned/ (0.0000001 + s0->p.nPairs)) ;
  printf ("\n#:CompatiblePairs\t%ld\t%.3f%%", s0->nCompatiblePairs, 100.0 * s0->nCompatiblePairs/ (0.000001 + s0->nPairsAligned)) ;

  printf ("\n\n#:Reads\t%ld", s0->p.nReads) ;
  if (pp->nRawReads) printf ("\n#:RawReads\t%ld", pp->nRawReads) ;
  printf ("\n#:Reads_aligned\t%ld\t%.3f%%", s0->nMultiAligned[0], 100.0 * s0->nMultiAligned[0]/ (0.0000001 + s0->p.nReads)) ;
  printf ("\n#:PerfectReads\t%ld\t%.2f%%", s0->nPerfectReads, (100.0 * s0->nPerfectReads)/(s0->p.nReads + .000001)) ;

  printf ("\n#:Reads supporting introns %ld, plus %ld, minus %ld, stranding %.2f%%\n"
	  , s0->nIntronSupportPlus + s0->nIntronSupportMinus 
	  , s0->nIntronSupportPlus
	  , s0->nIntronSupportMinus
	  , s0->intronStranding
	  ) ;
  printf ("\n#:SupportedIntrons\t%ld", s0->nSupportedIntrons) ;


  printf ("\n\n#:Bases\t%ld", s0->p.nBase1 + s0->p.nBase2) ;
  if (pp->nRawBases) printf ("\n#:RawBases\t%ld", pp->nRawBases) ;
  printf ("\n#:Bases_aligned\t%ld\t%.3f%%"
	  , s0->nBaseAligned1 + s0->nBaseAligned2
	  , 100.0*(s0->nBaseAligned1 + s0->nBaseAligned2)/(.0001 + s0->p.nBase1 + s0->p.nBase2)
	  ) ;

  printf ("\n#:Mismatches\t%ld\t%.5f%%"
	  , s0->nErr
	  , 100.0 * s0->nErr /(.0001 + s0->nBaseAligned1 + s0->nBaseAligned2)
	  ) ;
  

  printf ("\n") ;
  ac_free (h) ;
} /* reportRunStats */

/*************************************************************************************/
/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

void saUsage (char *message, int argc, const char **argv)
{
  int i ;

  if (message)
    {
      fprintf (stderr, "########## ERROR : %s\n", message) ;
      fprintf (stderr, "########## Please try : sortalign -h\n\n") ;
    }
  else if (argc == 0)
    {

      fprintf (stderr,
	       "// Usage: sortalign  <parameters> \n"
	       "//      try: -h --help --version\n"
	       "//\n"
	       "// EXAMPLES:\n"
	       "//      sortalign --createIndex XYZ -t target.fasta (needed once)\n"
	       "//      sortalign --index XYZ -i f_1.fastq.gz+f_2.fastq.gz --wiggles -o results/xxx \n"
	       "//      sortalign --index XYZ -i SRR35876976  -o results/xxx \n"
	       "//\n"
	       "// OBJECTIVE:\n"
	       "//    Sort align maps deep-sequencing files and can export at once alignments, coverage plots, introns.\n"
	       "//      On the first pass, sortalign analyses the (-t or -T) target(s) and creates a binary index.\n"
	       "//      On subsequent calls, sortalign memory maps the index and analyses the (-i or -I) sequence files.\n"
	       "//      All output files use the -o parameter (we recommend to give the name of a directory) as prefix.\n"
	       "//    The amount of data to be analysed can be large (say 100 Gbases), as it is processed in batches.\n"
	       "//      Depending on the  hardware, using 10 to 20 CPUs, it may process up to 1 Gigabases of human RNA-seq per minute\n"
	       "//      On human or mouse, the program requires around (30 + nB) gigabytes of RAM (see --nB --nA options).\n"
	       "//    Sequence files can be provided locall, or automatically downloaded from SRA as in the EXAMPLES above.\n"
	       "//\n"
	       "// TARGET CREATION:\n"
	       "// --createIndex <directory_name>\n"
	       "// --seedLength <int> : default 18, min 10, max 19\n"
	       "// --maxTargetRepeats <int> : default[81], do not index highly repeated target words\n"
	       "// -t <fasta_FileName>\n"
	       "//    Single target, for example a fasta file containing all the chromosomes of a given organism\n"
	       "// -T <target_configuration_fileName>\n"
	       "//    A tab delimited file describing a more complex protocol\n"
	       "//    Each line describes a target class, in a class specific format:\n"
	       "//      Column 1: one of G, M, C, R, E, A, I, B, V : target class.\n"
	       "//         G: Genome fasta files, RNA will be aligned on the genome jumping introns\n"
	       "//         M: Mitochondial fasta file\n"
	       "//         C: Chloroplast fasta file\n"
	       "//         R: Ribosomal RNA fasta file\n"
	       "//         E: Control sequences (i.e. ERCC)\n"
	       "//         A: Adaptor fasta file \n"


	       "//         I: Known Introns and Transcripts,\n"
	       "//            The I file must be of type .gtf, .gff, .gff3 [.gz], or .introns.\n"
	       "//            In each case the the chromosome name (column 1) must match the G fasta file(s).\n"
	       "//            In gtf/gff format the code expects (at least) 7  columns chrom/method/tag/a1/a2/./[+|-]\n"
	       "//               tag : [exon] other lines are ignored,\n"
	       "//               a1 a2 : the positions of the first and last base of t5he exon.\n"
	       "//               strand : [+|-] indicates the strand, a1 < a2, and the coordinates are 1-based.\n" 
	       "//            In .introns format the code expects 3 columns chrom/a1/a2\n"
	       "//               a1 a2 : give the positions of the 2 G of the Gt_aG motif,\n"
	       "//               on the top strand a1 < a2, on the bottom strand a1 > a2.\n"
	       "//            Please create the .introns file if the gff/gtf file is not available\n"

	       "//         B: Frequent bacteria contaminants fasta file.\n"
	       "//         V: Frequent virus contaminants fasta file.\n"
	       "//      Column 2: f1[,f2] : one or several coma separated file names.\n"
	       "//                For a given class (G, ...) several files can also be specified on several lines.\n"
	       "//      Column 3: optional, the file format if not implied by the file name\n"
	       "//           fasta[default], fastq, fastc, raw.\n"
	       "//\n"
	       "// SEQUENCE FILES  TO BE ANALYZED:\n"
	       "// -x, --index <directory_name> : a directory of index files created  previously using --createIndex\n" 
	       "// -r, --run <runName>   [default -]  : global run name\n"
	       "// -i <sequence_fileName[s]> \n"
	       "//       example 1:   -i f0.fasta,f1.fasta.gz,f2.fastq,f3.fastq.gz,f4.raw,f5.fasta\n"
	       "//        Align files f0 to f5, possibly using different formats.\n"
	       "//       example 2:   -i f1_1.fastq+f1_2.fastq,f2_1.fastq.gz+f2_2.fastq.gz\n"
	       "//        Align read pairs, again possibly using different formats.\n"
	       "//     The file format is implied by the file name, or may be provided explicitly using\n"
	       "//       --raw   | --fasta  | --fastq | --fastc | -sra\n"
	       "//     for example when using  '-i - ' to pipe sequence files into the pipeline\n"
	       "//       example 3: zcat fx.*.fastq.gz | sortalign -x XYZ -i --fastq -o results/fx\n"
	       "// -I <config_fileName>\n"
	       "//     A more precise definition of a set of sequencing files to be analysed\n"
	       "//     Eech line contains 1 to 3 tab separated columns\n"
	       "//          FileName[s]  RunName Descriptors\n"
	       "//     Example:\n"
	       "//          f1.fasta run1 RNA,Nanopore\n"
	       "//          f.R1.gz+f.R2.gz  run2  DNA,fastq,Illumina\n"
	       "//        A table with a single column of sequence files is acceptable, the format will be deduced from the file names\n"
	       "//     1: file name, mandatory\n"
	       "//        For paired end sequencing, provide, as in the second example (f1.F+f1.R), two file names separated by a plus\n"
	       "//        A file named (.fa, .fasta, .fastq, .fastc) implies that format, the file may be gzipped (.gz)\n"
	       "//        A file named (SRR[0-9]*) implies direct SRA download, with optional caching in the ./SRA directory\n"
	       "//        A file named (.sra.fasta.gz) is a fasta file with read-pairs on successive lines\n"
	       "//     2: Run name, optional,the name must not contain spaces or special characters, as it will be used to name subdirectories\n"
	       "//     3: Descriptors, optional, the options are coma separated, possibilities are\n"
	       "//        DNA/RNA : [default RNA], if DNA sortalign will not clip polyA or jump introns\n"
	       "//        Machine : [default Illumina], one of Illumina, PacBio, Nanopore,..\n"
	       "//        Adaptors=atagg,cctg   run specific adaptor sequences\n"
	       "//        Format : [default fasta], one of sra, raw, fasta, fastq, fastc, only needed if not implied by the file name\n"
	       "//  --sraCaching\n"
	       "//    When the requested sequences are is SRA format they are automatically downloaded from NCBI SRA archive\n"
	       "//    If --sraCaching is set, the (large) downloaded fasta.gz files are copied and saved in the ./SRA local directory\n"
	       "//    This is only useful if you intend to align these sequences several times or reuse them in other programs\n" 
	       "//  --sraDownload SRR123,SRR456,SRR999 \n"
	       "//    Just download a set of SRR entries into the cache directory ./SRA/SRR123.fasta.gz ... \n"
	       "//    This command is optional, since sortalign provides direct SRA access\n"
	       "//    --maxGB <int> : [default no limit]\n"
	       "//      Limit the sraDownload to <int> Gigabases per SRR entry\n" 
	       "//  --gzi\n"
	       "//    Forces decompression of the input files, useful when piping into sortalign\n"
	       "//    All files named *.gz are automatically decompressed\n"
	       "// OUTPUT:\n"
	       "// -o <outFileNamePrefix> [default stdout]\n"
	       "//   All output files will share this prefix, large outputs are split\n"
	       "//   Examples : -o x/y/z  all output files will be named x/y/z.something\n"
       	       "//              -o x/y/   all output files will be named x/y/something\n"
	       "// --gzo : the output files will be gziped\n" 

	       "// ALIGNMENT FORMAT:\n"
	       "//   By default, the aligments are exported in .hits format\n"
	       "//   To request bam/sam or to skip these files use\n"
	       "// --no_ali : do not export the alignements (or the introns unless --introns is set)\n"
	       "// --sam : export the alignment in sam format\n"
	       "// --bam : export the alignment in bam format\n"
	       "// --hits : also export .hits in addition to sam/bam\n"
	       "//    --quality_factors: export in sam/am the quality factors read fro SRA or from fastq files\n"
	       "//      These coefficients are not used by this aligner, just echoed to sam/bam on request\n"
	       "//\n"
	       "// DATA FILTERS:\n"
	       "// --minReadLength  <int> : [default 20 (or max read lenght)]\n"
	       "//     Drop reads shorter than this limit\n"
	       "//     The default 20 is lowered to max-read-length if shorter, simplifying short-RNA analysis\n"
	       "// --minEntropy  <int> : [default 20 (or max read lenght/2)]\n"
	       "//     Drop reads withg base-4 entropy shorter than this limit\n"
	       "//     The default 20 is lowered to max-read-length/2 if shorter, simplifying short-RNA analysis\n"
	       "// REQUEST\n"
	       "// --RNA | --DNA : [optional] the code automatically distinguises RNA from DNA\n"
	       "//     By setting one of these option, you impose the alignment mode, otherwise you let the code choose\n"
	       "//     In RNA mode, the code favors introns and detects polyA tails\n"  
	       "// --align : [default] Extend the seed alignments hopefully to full read length\n"
	       "// --do_not_align : do not align, just test the validity of the index directory\n"
	       "//     The score s is computed as s = #aligned bases - errCost * #errors\n"
	       "//   --minAli <int> : [default 30] minimal length of the alignment\n"
	       "//   --minScore <int> : [default 30] minimal score  of the alignment\n"
	       "//   --errCost <int> : [default 8] cost of substitition, or short indel up to 3 bases\n"
	       "//   --errMax <int> : [default NA] maximal number of mismatches in any (partial) alignment\n"
	       "//   --errRateMax <int> : [default 10] maximal percentage of mismatches in any (partial) alignment\n"
	       "//   --no_splice : only accept comtinuous alignments, [by default search also spliced alignments]\n"
	       "//   --maxIntron [default 1000000] : max intron size\n"
	       "//   --ignoreIntronSeeds [default FALSE] : do not use the known intron provided in class I in the -T config file\n"
	       "//   --adaptor1 atgc : Read 1 exit adaptor, as exported in file <run>.overhang3prime.tsf\n"
	       "//   --adaptor2 ggct : Read 2 exit adaptor, as exported in file <run>.overhang5prime.tsf\n"
	       "// --wiggles  : Report target coverage wiggles in UCSC BF (fixed) format\n"
	       "// --introns  : [default] Report intron support\n"
	       "// (--snp : not yet ready) Report candidate SNP counts (substitutions and short indels)\n"
	       "// STEPPING\n"
	       "// --step <int>, take a seed every <int> base\n"
	       "//   while createIndex, the default is 1 for targets < 1Mb, 3 for larger targets, 5 would accelerate the code\n"
	       "//   while analysing the reads, the default is 2. --step 1 would increase sensitivity\n"
	       "//\n"
	       "// OPTIONAL TECHNICAL PARAMETERS\n"
	       "//   Unless the program refuses to run, it is probably best not to specify these parameters\n"
	       "//   The program sortalign is parallelized and synchronized by GO-language like channels\n"	
	       "//     This is very different from the way most C programs are parallelized\n"
	       "//     All code layers execute at the same time, the work-load among layers is self balancing\n"
	       "//     and each agent and communication channel runs in its own virtual thread.\n"
	       "//   There is no limit to the size of the input, as data continuously flow in and out\n"
	       "// --nA or --nAgents <int> : [default 10] number of agent in each code layer\n"
	       "// --nB or --nBlocks <int> : [default 20] number of simultaneous data blocks circulating in the pipeline\n"
	       "// --bMax <int> : [default 20] (range 1--1024) max number of mega-bases in a data block.\n"
	       "// --NN <int> : [default 16] split the seed files in NN parts, allowed values (1,2,4,8,16,32,64)\n"
	       "//              seedLength > 16 imply NN >= 4^(sedLength - 16)\n"
	       "// --max_threads <int>  : [default 128] maximal number of simultaneous UNIX threads\n"
	       "//    Do not impose a low value, like --max_threads 8 or 16, the code will not run correctly\n"
	       "//    Possible values are \n"
	       "//        --nBlocks 1 : for a small test\n"
	       "//        --nAgents 10 --nBlocks 10 --max_threads 128 : default, uses less than 32 Gbytes of RAM\n"
	       "//        --nAgents 20 --nBlocks 30 --max_threads 512 : on a large machine\n"
	       "//    These parameters do not limit the amount of data to be processed, they affect the speed\n"
	       "//    If no message is emited within the first minute, the computer is overloaded,\n"
	       "//    please try to increase --max_threads to 256 or 512. alternativelly lower --nAgents and --nBlocks\n"
	       "//      The code is using virtual threads, at least one by agent and by communication channels.\n"
	       "// --noJump : Do not insert jumper in the index, valid for future GPU version\n"
	       "// --verbose : all kinds of details, mostly usefull for debugging, are reported\n"
	       "//\n\n"
	       ) ;
    }
  else
    {
      fprintf (stderr,
	       "########## ERROR: I do not understand the argument%s ", argc > 2 ? "s" : "") ;
      for (i = 1 ; i < argc ; i++)
	fprintf (stderr, "%s ", argv[i]) ;
      fprintf (stderr, "\n##########\t It may have been misspelled or repeated\n") ; 
      fprintf (stderr, "##########\t Please try : sortalign -h\n\n") ;
    }
  exit (1) ;
} /* saUsage */

/*************************************************************************************/
/************ Public interface, get params, set channels and agents, launch **********/
/*************************************************************************************/

int main (int argc, const char *argv[])
{
  PP p ;
  BB bb ;
  Array cpuStats = 0 ;
  Array runErrors = 0 ;
  Array inArray = 0 ;
  int maxThreads = 128 ;
  long unsigned int nHits = 0, nSeqs = 0, nBaseAligned = 0 ;
  long unsigned int nerr = 0 ;   /* cumulated number of errors */
  long unsigned int aliDx = 0 ; /* cumulated aligned read length */
  long unsigned int aliDa = 0 ;  /* cumulated genome coverage */
  /*   long unsigned int intronSupports = 0 ;  // cumulated intronSupports */ 
  long unsigned int nAli = 0 ;
  BOOL debug = FALSE ;
  int n = 0 ;
  long int skips0 = 0, skips1 = 0, skips2 = 0, skips3 = 0, skips4 = 0, skipsFound = 0, skipsNotFound = 0 ;
  char tBuf[25] ;
  char tBuf0[25] ;
  const char *species = 0 ;
  AC_HANDLE h ;
  int nAgents = 10 ;
  int channelDepth = 10 ;
  mytime_t t0, t1 ;
  
  freeinit () ; 
  messErrorInit (argv[0]) ;

  memset (&p, 0, sizeof (PP)) ;
  h = ac_new_handle () ;
  p.h = h ;

  /**************************  trivial parameters *********************************/
  
  if (argc < 2)
    saUsage (0, 0, argv) ;
  if (getCmdLineBool (&argc, argv, "-help")  ||
      getCmdLineBool (&argc, argv, "--help") ||
      getCmdLineBool (&argc, argv, "-h")
      )
    saUsage (0, 0, argv) ;

  if (getCmdLineBool (&argc, argv, "--version"))
    {
      fprintf (stderr, "sortalign version 0.1.56, jan 2027") ;
      exit (0) ;
    }     

  p.debug  = getCmdLineBool (&argc, argv, "--debug") ;
  p.debug |= getCmdLineBool (&argc, argv, "--verbose") ;
  
  
  /**************************  SRA downloader *********************************/
  
  {{  /* --sraDownload SRR1,SRR2,,...
       * Just download SRA files from NCBI/SRA into the local SRA caching directory
       * This is not needed by the aligner, but is provided as a service to the user
       * since it costs us nothing as we have a real-time downloader a11vailable
       */
      
      const char *sraID = 0 ;
      if (getCmdLineText (h, &argc, argv, "--sraDownload", &sraID))
	{
	  char *cp = 0 ;
	  const char *outFormat = 0 ;

	  p.sraOutFormatPE = TRUE ; /* default */
	  if (getCmdLineText (h, &argc, argv, "--O", &(outFormat)))
	    {
	      cp = strnew (outFormat, h) ;
	      p.sraOutFormatPE = FALSE ; /* default */
	      while (cp)
		{
		  char *cq = strchr (cp, ',') ;
		  if (cq) *cq = 0 ;
		  if (! strcmp (cp, "PE"))
		    p.sraOutFormatPE = TRUE ;
		  else if (! strcmp (cp, "PEQ"))
		    p.sraOutFormatPEQ = TRUE ;
		  else if (! strcmp (cp, "fasta"))
		    p.sraOutFormatFasta = TRUE ;
		  else if (! strcmp (cp, "fastq"))
		    p.sraOutFormatFastq = TRUE ;
		  else
		    saUsage ("-O parameter should of one or several of  PE,PEQ,fasta,fastq", argc, argv) ;
		  
		  cp = cq ? cq + 1 : 0 ;
		}
	    }
		    
	  getCmdLineInt (&argc, argv, "--maxGB", &(p.maxSraGb)) ;
	  if (p.maxSraGb < 0)
	    saUsage ("--maxGB parameter should be positive", argc, argv) ;

	  cp = strnew (sraID, h) ;
	  while (cp)
	    {
	      char *cq = strchr (cp, ',') ;
	      if (cq) *cq = 0 ;
	      saSequenceParseSraDownload (&p, cp) ;
	      cp = cq ? cq + 1 : 0 ;
	    }
	  exit (0) ;
	}
    }}

  /***************************************************************************************/

  {{
      int n = 0 ;
      if (getCmdLineInt (&argc, argv, "--testRG", &(n)))
	{
	  getCmdLineText (h, &argc, argv, "-o", &(p.outFileName)) ;
	  saCreateRandomGenome (&p, n) ;
	  exit (0) ;
	}
    }}

  /******************** NUMA harware  optimizer ******************************************/
  /******************** pin the threads to the least used node ***************************/
  /* numa (non unifirm mmory access) harware are composed of several nodes, each with many CPUs
   * moving data across nodes is slow and costly
   * This module finds the node with the least number of running threads
   * and relaunches the program pinned to that node
   * adding the command line parameter --numactl to prevent recursion
   *
   * This system could be useful in other C programs using multithreading and large memory
   */   

  if (! p.debug &&
      ! getCmdLineBool (&argc, argv, "--numactl")  &&
      !  (getenv("INVOCATION_NOTIFICATIONS") &&
	  strstr(getenv("INVOCATION_NOTIFICATIONS"), "numactl")) &&
      isExecutableInPath ("numactl")   /* see w1/utils.c */
      )
    {
      vTXT txt = vtxtHandleCreate (h) ;
      
      int nodes = saBestNumactlNode () ;
	
      vtxtPrintf (txt, "/usr/bin/numactl  --cpunodebind=%d --membind=%d ", nodes, nodes) ;
      for (int i = 0 ; i < argc ; i++)
	vtxtPrintf (txt, " %s " , argv[i]) ;
      vtxtPrintf (txt, " --numactl ") ;

      fprintf (stderr, "%s\n", vtxtPtr (txt)) ;

      return system (vtxtPtr (txt)) ;
    }
  
  /**************************  debugging modules, ignore *********************************/

  getCmdLineText (h, &argc, argv, "-o", &(p.outFileName)) ;

  p.gzi = getCmdLineBool (&argc, argv, "--gzi") ;   /* decompress input files (implicit for files named .gz) */
  p.gzo = getCmdLineBool (&argc, argv, "--gzo") ;   /* compress most output files */

  getCmdLineInt (&argc, argv, "--maxGB", &(p.maxSraGb));
  if (p.maxSraGb < 0)
    saUsage ("--maxGB parameter should be positive", argc, argv) ;

  if (getCmdLineText (h, &argc, argv, "--species", &species))
    {
      if (species && ! strcasecmp (species, "worm"))
	p.isWorm = TRUE ;
    }

  if ( getCmdLineBool (&argc, argv, "--checkDict"))
    {                                           /* check the dictionary lib */
      DICT *dict = dictHandleCreate (1024, p.h) ;
      ACEIN ai = aceInCreate ("myErrors1.txt", 0, h) ;
      int n = 0, nnn = 0 ;
      char *cp ;
      
      aceInSpecial (ai, "\n") ;
      while (aceInCard (ai))
	{
	  cp = aceInWord (ai) ;
	  if (! cp || ! *cp)
	    continue ;
	  if (!strncmp (cp, "nnE=", 4))
	    {
	      aceInStep (ai, '\t') ;
	      cp = aceInWord (ai) ;
	      if (! cp || ! *cp)
		continue ;
	    }
	  dictAdd (dict, cp, &n) ;
	  nnn++ ;
	}
      fprintf (stderr, "success: loaded %d words\n", nnn) ;
      
      exit (0) ;
    }

  /**************************  debugging modules, ignore *********************************/

  if ( getCmdLineInt (&argc, argv, "--checkSorter", &n))
    {                                           /* check the sorting algorithm */
      BigArray b, a = bigArrayHandleCreate (n, CW, 0) ;

      
      for (int i = 0 ; i < n ; i++)
	bigArray (a, i, CW).seed = ( !(n &0x1) ? (n - i)%7 :  (unsigned int)randint()) ;  

      b = bigArrayHandleCopy (a, 0) ;
      saSort (b, 1) ; /* hitReadOrder */
      
      for (int i = 0 ; i < 50 && i < n ; i++)
	printf("%d %d\n"
	       , bigArr (a,i, CW).seed
	       , bigArr (b,i, CW).seed
	       ) ;
      for (int i = 0 ; i < n-1 ; i++)
	if (bigArr (b,i, CW).seed > bigArr (b,i+1, CW).seed)
	  messcrash ("\nerror line %d\n", i) ;
      showHits (0) ;
      showHitsDo (0, 0) ;
      exit (0) ;
    }

  /**************************  another debugging module, ignore ***************************/
  
  if (getCmdLineInt (&argc, argv, "--makeTest", &n))
    {
      ACEOUT ao = 0 ;
      BigArray aa2 =0, aa = bigArrayHandleCreate (6, int, h) ;
      for (int i = 0 ; i < 6 ; i++)
	bigArray (aa, i, int) = i+1 ;
      bigArrayMapWrite (aa, "mapTest1") ; 
      aa2 = bigArrayMapRead ("mapTest1", int, TRUE, h) ;
      for (int i = 0 ; i < 6 ; i++)
	printf ("\t%d", bigArray (aa2, i, int)) ;
      printf ("\n") ;
      
      /* create a test fasta file with a known 16 mer every 1000 bases */
      ao = aceOutCreate (hprintf (h, "testR.%d", n),  ".fasta", 0, h) ;
      char dna[1000 *n + n] ;
      char w[27] ;
      int i, j, k ;
      char *bb = "atgc" ;
      for (i = k = 0 ; i < n ; i++)
	{
	  for (j = 0 ; j < 26 ; j++)
	    dna[k++] = w[j] = bb[randint () % 4] ;
	  w[26] = 0 ; w[15] = dna[k-11] = 'c' ;
	  aceOutf (ao, ">s.%d\n%s\n", i, w) ;
	  for (j = 0 ; j < 974 ; j++)
	    dna[k++] = bb[randint () % 3] ;
	}
      dna[k++] = 0 ;
      if (n == 1) dna[27] = 0 ; /* expect 10 genome words */
      ao = aceOutCreate (hprintf (h, "testG.%d", n),  ".fasta", 0, h) ;
      aceOutf (ao, ">g.%d\n%s\n", n, dna) ;
      ac_free (h) ;
      exit (0) ;
    }

  /***************** creation of a new target index  *************/
  /* The target index must exist or be created
   *
   * The seed length is only used when creating a new target index
   * otherwise it is inherited from the index
   * The maximal number of repeats of a seed in the targets is used at creation
   * but one is allowed  provide a smaller number when aligning
   * which would drop the index seeds with more repetitions 
   */
  
  p.seedLength = 0 ; /* default */
  if (getCmdLineInt (&argc, argv, "--seedLength", &(p.seedLength)))
    {
      if (p.seedLength < 12)
	messcrash ("\n parameter --seedlength %d must be at least 12, sorry\n", p.seedLength) ;
      else if (p.seedLength > 19)
	messcrash ("\n parameter --seedlength %d cannot excered 19, sorry\n", p.seedLength) ;
    }

  p.maxTargetRepeats = 31 ;  /* was 81  31 12 */
  getCmdLineInt (&argc, argv, "--maxTargetRepeats", &p.maxTargetRepeats) ;


  getCmdLineText (h, &argc, argv, "-x", &(p.indexName)) ;
  getCmdLineText (h, &argc, argv, "--index", &(p.indexName)) ;
  getCmdLineText (h, &argc, argv, "-t", &(p.tFileName)) ;
  getCmdLineText (h, &argc, argv, "-T", &(p.tConfigFileName)) ;
  p.createIndex = getCmdLineText (h, &argc, argv, "--createIndex", &(p.indexName)) ;
  p.noJump = getCmdLineBool (&argc, argv, "--noJump") ;
  
  if (p.createIndex)
    {
      if (! p.tFileName && ! p.tConfigFileName)
	saUsage ("--createIndex requires providing a target parameter -t or -T", argc, argv) ;
      if ( p.tFileName &&  p.tConfigFileName)
	saUsage ("conflicting parameters -t and -T, both define the targets", argc, argv) ;
    }
  
  if (! p.indexName)
    saUsage ("missing parameter -x <indexName> or --createIndex <indexName> or missing indexName", argc, argv) ;

  /******************* sequence files to be analysed , not used if --createIndex *************/

  
  getCmdLineText (h, &argc, argv, "-i", &(p.inFileName)) ;
  getCmdLineText (h, &argc, argv, "-I", &(p.inConfigFileName)) ;

  /* check that input files were provided */
  if (! p.createIndex &&
      ! p.inFileName && ! p.inConfigFileName)
    saUsage ("missing parameter -i or -I inputFileName(s)", argc, argv) ;
  if (p.inFileName && p.inConfigFileName)
    saUsage ("conflicting parameter -i and -I, both define the input files", argc, argv) ;

  /***************** requested analyses, not used if --createIndex *************/

  p.deduplicate = getCmdLineBool (&argc, argv, "--deduplicate") ;
  p.isDna = getCmdLineBool (&argc, argv, "--DNA") ;
  p.isRna = getCmdLineBool (&argc, argv, "--RNA") ;
  p.justStats  = getCmdLineText (h, &argc, argv, "--justStats", 0) ;
  p.align = getCmdLineBool (&argc, argv, "--align") ;          /* left for back compatibility, over-ridden by --do_not_align */
  p.align = ! getCmdLineBool (&argc, argv, "--do_not_align") ; /* default is to align */
  p.ignoreIntronSeeds = getCmdLineBool (&argc, argv, "--ignoreIntronSeeds") ;

  p.hitsFormat = TRUE ; /* default */
  p.introns = TRUE ; /* default */
  p.sam = getCmdLineBool (&argc, argv, "--sam") ;
  p.bam = getCmdLineBool (&argc, argv, "--bam") ;
  p.qualityFactors = getCmdLineBool (&argc, argv, "--quality_factors") ;
  if (p.sam || p.bam)
    {
      p.hitsFormat = FALSE ;
      p.exportSamSequence = TRUE ;
      if (p.qualityFactors)
	p.exportSamQuality = TRUE ;
    }
  /* we may wish both sam/bam and hits */
  /* we may cancel all the alignment files */
  if (getCmdLineBool (&argc, argv, "--no_ali"))
    p.sam = p.bam = p.hitsFormat = p.exportSamQuality = p.exportSamSequence = p.introns = FALSE ;
  if (getCmdLineBool (&argc, argv, "--hits"))
    p.hitsFormat = TRUE ;  /* stay on previous value is not set */
  if (getCmdLineBool (&argc, argv, "--introns"))
    p.introns = TRUE ;  /* stay on previous value is not set */
  
  p.sraCaching = getCmdLineBool (&argc, argv, "--sraCaching");   /* cache the files downlaoded from NCBI/SRA */

  /* action options ***/
  p.wiggle = getCmdLineBool (&argc, argv, "--wiggles") ;
  p.wiggleEnds = getCmdLineBool (&argc, argv, "--wiggleEnds") ;
  p.snps = getCmdLineBool (&argc, argv, "--snp") ;

  getCmdLineText (h, &argc, argv, "--adaptor1", &(p.rawAdaptor1R)) ;
  getCmdLineText (h, &argc, argv, "--adaptor2", &(p.rawAdaptor2R)) ;

  p.wiggle_step = 0 ;  /* examples s=10, 5, 1, if not set by user the default is set in saConfigCheckTargetIndex  */
  getCmdLineInt (&argc, argv, "--wiggleStep", &(p.wiggle_step)) ;

  if (0) { p.sam = TRUE ; p.wiggle = TRUE ; p.wiggleEnds = FALSE ;}
  /*****************  sequence file names and their formats  ************************/
  
  p.runName = 0  ; /* default */

  /* force the sequence file formats,
   * in -I case it is better to provide the format run file by file in the third column of the iConfig file
   * but the best way is not to provide this parameter and let the the format be implied by the file names
          .fasta[.gz]
	  .fastq[.gz]
	  -fna[.gz]
	  .fa[.gz]
	  .fastc[.gz]
	  SRR*
   */
  p.fasta = getCmdLineBool (&argc, argv, "--fasta");
  p.fastq = getCmdLineBool (&argc, argv, "--fastq");
  p.fastc = getCmdLineBool (&argc, argv, "--fastc");
  p.raw = getCmdLineBool (&argc, argv, "--raw");
  p.sra = getCmdLineBool (&argc, argv, "--sra");
  
  /***************** run name  ***********/

  /* force the run name
   * in -I case it is better to provide the run name file by file in the second column of the tConfig file
   * several files can be associated to the same run
   */
  p.run = 0 ;
  getCmdLineText (h, &argc, argv, "-r", &(p.runName)) ;
  getCmdLineText (h, &argc, argv, "--run", &(p.runName)) ;

  /***************** method name, only used in some output files, convenient when optimizing parameters */
  
  getCmdLineText (h, &argc, argv, "--method", &(p.method)) ;
  
  /*****************  Optional technical parameters ************************/
  /*****************  number of parallel target indexes, only used when createIndex  ************************/

  NN = 16 ; /* default */
  getCmdLineInt (&argc, argv, "--NN", &(NN));
  
  if (NN != 1 && NN != 2 && NN != 4 && NN != 8 && NN != 16 && NN != 32 && NN != 64)
    messcrash ("The number of indexes NN=%d must be  a power of 2, say 1 2 4 8 ...", NN) ;
  p.nIndex = NN ;
  
  /* harware check : the code can oly run on 64 bits machnes */
  if (sizeof (long unsigned int) != 8)
    messcrash ("The source code assumes that long unsigned ints use 64 bits not %d, sorry", 8 * sizeof (long unsigned int)) ;

  /***************** amount or parallelization **************************/
  nAgents = NAGENTS ;
  if (! getCmdLineInt (&argc, argv, "--nAgents", &(nAgents)))
    getCmdLineInt (&argc, argv, "--nA", &(nAgents)) ;

  p.nBlocks = NBLOCKS ;  /* max number of BB blocks processed in parallel */
  if (p.debug)      /* under debugger, it is more convenient to run with a single agent and a single block */
    p.nBlocks = 1 ; /* but we wish to be able to reset nBlocks even in debug mode */
  if (!getCmdLineInt (&argc, argv, "--nBlocks", &(p.nBlocks)))
    getCmdLineInt (&argc, argv, "--nB", &(p.nBlocks));
  if (p.nBlocks == 1)
    { nAgents = 1 ; }
  channelDepth = 1 + p.nBlocks/2 ;
  maxThreads = 128 ;  /* UNIX  max on lmem12 machine */
  getCmdLineInt (&argc, argv, "--max_threads", &maxThreads) ;
  if (maxThreads < 24)
    maxThreads = 24 ;
  if (p.nBlocks == 1)
    maxThreads = 16 ;  /* if maxThreads == 8 and the single fasta file is split in 3 parts, the code stalls */
  if (p.createIndex)
    { nAgents = 1 ; maxThreads = 1 ; p.nBlocks = 1 ; }


  /****************** stepping when constructing sequence seeds ********************/
  p.iStep = 0 ;   /* read default = tStep/2 set in  saConfigCheckTargetIndex */
  p.tStep = 0 ;   /* default 2 or 4 set in createIndex or read in existing index */
  if (p.createIndex)
    {
      getCmdLineInt (&argc, argv, "--step", &(p.tStep)) ;
      if (p.tStep & 0x1 ) p.tStep++ ; /* impose an even number so that default iStep=tStep/2 divides tStep */
    }
  else
    getCmdLineInt (&argc, argv, "--step", &(p.iStep)) ;

  getCmdLineLong (&argc, argv, "--nRawReads", &(p.nRawReads)) ;
  getCmdLineLong (&argc, argv, "--nRawBases", &(p.nRawBases)) ;

  /*****************  Aligner filters  ************************/

  p.minLength = -20 ;  /* use a negative default to differentiate from a user's values */
  p.minEntropy = -20 ;
  getCmdLineInt (&argc, argv, "--minLength", &(p.minLength)) ;
  getCmdLineInt (&argc, argv, "--minEntropy", &(p.minEntropy)) ;

  p.minAli = p.minAliPerCent = p.minScore = -1 ;
  p.errMax = 1000000 ; /* on negative values, extendHits stops on first error */
  p.errRateMax = 10 ;
  p.OVLN = 30 ;
  p.splice = TRUE ;
  if (getCmdLineBool (&argc, argv, "--no_splice"))
    p.splice = FALSE ;
  p.errCost = 4 ; /* was 8 */
  getCmdLineInt (&argc, argv, "--errCost", &(p.errCost)) ;
  getCmdLineInt (&argc, argv, "--errMax", &(p.errMax)) ;
  getCmdLineInt (&argc, argv, "--minScore", &(p.minScore)) ;
  getCmdLineInt (&argc, argv, "--minAli", &(p.minAli)) ;
  getCmdLineInt (&argc, argv, "--minAliPerCent", &(p.minAliPerCent)) ;
  getCmdLineInt (&argc, argv, "--errRatMax", &(p.errRateMax)) ;
  p.maxIntron = 1000000 ;
  getCmdLineInt (&argc, argv, "--maxIntron", &(p.maxIntron)) ;

  if (p.minScore < 0) p.minScore = 0 ;
  if (p.minAli < 0) p.minAli = 30 ;
  if (p.minAliPerCent < 0) p.minAliPerCent = 0 ;
  if (p.minAli < p.minScore) p.minAli = p.minScore ;

  p.BMAX = 20 ;
  getCmdLineInt (&argc, argv, "--bMax", &(p.BMAX)) ;
  if (p.BMAX < 1) p.BMAX = 1 ;
  if (p.BMAX > 1024) p.BMAX = 1024 ;

 
  /*****************  Check tat all parameters have been parsed *******************/
  if (argc > 1)
    saUsage (0, argc, argv) ;

  showAli (0) ; /* for compiler happiness */

  /****************** Check the existence of all file names ****************************************/ 

  /* check that the index directory is accessible */
  saConfigIsIndexAccessible (&p) ;

  /* check that the output directoy is writable */
  saConfigIsOutDirWritable (&p) ;   
  oligoEntropy (0, 0, 0) ; /* inialize the logs */
  
  /*****************  Start working ***********************************************/
  t0 = timeNow () ;
  printf ("%s: Start\n", timeBufShowNow (tBuf0)) ;

  cpuStats = arrayHandleCreate (1024, CpuSTAT, h) ;
  p.runStats = arrayHandleCreate (1024, RunSTAT, h) ;
  runErrors = arrayHandleCreate (1024, Array, h) ;
  
  p.runDict = dictHandleCreate (16, p.h) ;
  p.targetClassDict = dictHandleCreate (16, p.h) ;

  dictAdd (p.targetClassDict, "rRNA", 0) ;
  dictAdd (p.targetClassDict, "Mito", 0) ;
  dictAdd (p.targetClassDict, "Genome", 0) ;
  dictAdd (p.targetClassDict, "Bacteria", 0) ;
  dictAdd (p.targetClassDict, "Virus", 0) ;

  if (p.sam || p.bam) p.bamHandle  = ac_new_handle () ;
  
  /*******************  create the index ********************************************/

  p.tArray = saTargetParseConfig (&p) ;
  if (p.createIndex)
    { /* The human genome index consumes around 18 Gigabytes of RAM */
      saTargetIndexCreate (&p) ;
      goto done ;
    }

  /*******************  otherwise verify the existence of the indexes ********************/
  
  NN = saConfigCheckTargetIndex (&p) ;

  /* check the existence of the input sequence files */
  inArray = saConfigGetRuns (&p, p.runStats) ;
  n = dictMax (p.runDict) + 1 ;
  p.runLanes = arrayHandleCreate (n, atomic_int, p.h) ;
  p.runLanesDone = arrayHandleCreate (n, atomic_int, p.h) ;
  array (p.runLanes, n - 1, atomic_int) = 0 ;
  array (p.runLanesDone, n - 1, atomic_int) = 0 ;
  
  /* Set bonus for over represented sequences */
  p.bonus['M'] = 4 ; /* mitochondria */
  p.bonus['C'] = 4 ; /* chloroplast */
  p.bonus['R'] = 8 ; /* rRNA */
  p.bonus['E'] = 8 ; /* ERCC spikeIns */
  p.bonus['A'] = 8 ; /* Adaptors */
  p.bonus['B'] = -12 ; /* Bacteria */
  p.bonus['V'] = -12 ; /* Virus */


  /******************** launch the multiprocessing ***************************************/

  wego_max_threads (maxThreads) ;
  if (0)   channelDepth = 10 ;

  if (1)
    {
      /* Create the communication channels */
      p.fpChan = channelCreate (4096, RC, p.h) ; /* this chan nust be deep enough to accept all file names at once */
      channelDebug (p.fpChan, debug, "fpChan") ;
      p.npChan = channelCreate (6, int, p.h) ; /* count BB per inFile */
      channelDebug (p.npChan, debug, "npChan") ;
      p.gmChan = channelCreate (1, BB, p.h) ;
      channelDebug (p.gmChan, debug, "gmChan") ;
      p.plChan = channelCreate (1, BB, p.h) ;
      channelDebug (p.plChan, debug, "plChan") ;
      p.lcChan = channelCreate (channelDepth, BB, p.h) ;
      channelDebug (p.lcChan, debug, "lcChan") ;
#ifndef YANN
      p.csChan = channelCreate (channelDepth, BB, p.h) ;
      channelDebug (p.csChan, debug, "csChan") ;
      p.smChan = channelCreate (channelDepth, BB, p.h) ;
      channelDebug (p.smChan, debug, "smChan") ;
      p.moChan = channelCreate (channelDepth, BB, p.h) ;
      channelDebug (p.moChan, debug, "moChan") ;
      p.oaChan = channelCreate (channelDepth, BB, p.h) ;
      channelDebug (p.oaChan, debug, "oaChan") ;
#endif
      p.aeChan = channelCreate (channelDepth, BB, p.h) ;
      channelDebug (p.aeChan, debug, "aeChan") ;
      p.doneChan = channelCreate (channelDepth, BB, p.h) ;
      channelDebug (p.doneChan, debug, "doneChan") ;
      
      
      /* create the agents using  "wego_go()"
       *   Their execution is  triggered by their in-channels
       *   They export processed BB data-blocks on their out-channel
       *   They exit when their in-channel closes
       *
       * a non-writable COPY of p is passed to each agent
       * this allows to pass agent specific parameters
       * like p.agent : the identifier of an instance of the agent
       */
      
      /* The genome parsers start immediately */
      p.agent = nAgents ;
      wego_go (saTargetIndexGenomeParser, &p, PP) ;
      /* The load regulator maintains --nB data-blocks in the pipeline */
      wego_go (loadRegulator, &p, PP) ;
      
      /* Start the processing of the sequence files */
      p.nFiles = arrayMax (inArray) ; /* number of files to be processed */
      /* n Files must be processed
       * The readParsers will emit for each file, in parallel, n BB blocks
       * The number is passed to the npChan which feeds the npCounter
       * When n numbers have been received, the npCounter
       * knows the cumulated number N of BB blocks to be analyzed
       * and tells the plChan, hence the load regulator,
       * and recurssibvely all program layers
       * to close after having processed N BB blocks
       */
      wego_go (npCounter, &p, PP) ; channelAddSources (p.plChan, 1) ;
      /* Read preprocessing agents, they do not require the genome */
      for (int pass = 0 ; pass < 2 ; pass++)
	for (int i = 0 ; i < p.nFiles && i < nAgents && i < 10 ; i++)
	  {
	    if (pass) fprintf (stderr, "Launch readParser %d\n", i) ;
	    p.agent = i ;
	    
	    if (pass) wego_go (readParser, &p, PP) ; else channelAddSources (p.npChan, 1) ;
	  }
      for (int pass = 0 ; pass < 2 ; pass++)
	for (int i = 0 ; i < nAgents && i < p.nBlocks ; i++)
	  {
	    p.agent = i ;
	    
#ifndef YANN
	    if (pass) wego_go (codeWords, &p, PP) ; else channelAddSources (p.csChan, 1) ;
	    if (pass) wego_go (sortWords, &p, PP) ; else channelAddSources (p.smChan, 1) ;
#endif
	  }

      
      /* Load the fpChan, triggering the readParsers
       * and recursively the whole pipeline
       */
      for (int k = 0 ; k < p.nFiles ; k++)
	{
	  RC *rc = arrayp (inArray, k, RC) ;
	  array (runErrors, rc->run, Array) = arrayHandleCreate (32, int, p.h) ;
	  
	  channelPut (p.fpChan, rc, RC) ;
	}
      channelClose (p.fpChan) ;
      
      /* wait untill the genome is ready */
      channelGet (p.gmChan, &p.bbG, BB) ;
      p.genomeLength = p.bbG.genomeLength ;
      
      if (! p.bbG.cwsN[0])
	messcrash ("matchHits received no target words") ;
      saGffBinaryParser (&p) ;
      
      /* map the reads to the genome in parallel */

      /* first we register the topology of our workflow (pass == 0)
       * then we launch (pass == 1)
       * otherwise the first agent which finishes could close the out channel
       */
       
      for (int pass = 0 ; pass < 2 ; pass++)
	for (int i = 0 ; i < nAgents && i < p.nBlocks ; i++)
	  {
	    p.agent = i ;
#ifdef YANN
	    if (pass)
	      wego_go (wholeWork, &p, PP) ;
	    else
	      channelAddSources (p.aeChan, 1) ;
#else

	  if (pass)
	    wego_go (matchHits, &p, PP) ;
	  else
	    channelAddSources (p.moChan, 1) ;
		  
	  if (pass)
	    wego_go (sortHits, &p, PP) ;
	  else
	    {
	      channelAddSources (p.oaChan, 1) ;
	      channelAddSources (p.doneChan, 1) ;
	    }
	      
	  if (!i || p.align) /* at least 1 agent */
	    {
	      if (pass)
		{
		  p.agent = 3*i ;
		  wego_go (saAlign, &p, PP) ;
		  p.agent = 3*i + 1 ;
		  wego_go (saAlign, &p, PP) ;
		  p.agent = 3*i + 2 ;
		  wego_go (saAlign, &p, PP) ;
		  p.agent = i ;
		}
	      else
		channelAddSources (p.aeChan, 3) ;
	    }
#endif
	  if ((! p.sam && ! p.bam)  || !i) /* only 1 export agent in sam case */
	    {
	      if (pass)
		wego_go (export, &p, PP) ;
	      else
		channelAddSources (p.doneChan, 1) ;
	    }
	}
    }

  p.confirmedSLs  = arrayHandleCreate (1000, SLS, p.h) ;
  p.confirmedPolyAs  = arrayHandleCreate (1000, POLYA, p.h) ;
  p.confirmedIntrons  = arrayHandleCreate (1000, INTRON, p.h) ;
  p.doubleIntrons  = arrayHandleCreate (1000, DOUBLEINTRON, p.h) ;
  int chromMax = dictMax (p.bbG.dict) + 1 ;
  int runMax = dictMax (p.runDict) + 2 ;
  p.wiggles = 0 ;
  p.wigglesL = 0 ;
  p.wigglesR = 0 ;
  p.wigglesP = 0 ;
  p.wigglesNU = 0 ;
  if (p.wiggle)
    {
      p.wiggles = arrayHandleCreate (2 * chromMax * runMax, Array, h) ;
      array (p.wiggles, 2 * chromMax * runMax - 1, Array) = 0 ; /* initialize */
      p.wigglesP = arrayHandleCreate (2 * chromMax * runMax, Array, h) ;
      array (p.wigglesP, 2 * chromMax * runMax - 1, Array) = 0 ; /* initialize */
      p.wigglesNU = arrayHandleCreate (2 * chromMax * runMax, Array, h) ;
      array (p.wigglesNU, 2 * chromMax * runMax - 1, Array) = 0 ; /* initialize */
    }
  if (p.wiggleEnds)
    {
      p.wigglesL = arrayHandleCreate (2 * chromMax * runMax, Array, h) ;
      p.wigglesR = arrayHandleCreate (2 * chromMax * runMax, Array, h) ;
      array (p.wigglesR, 2 * chromMax * runMax - 1, Array) = 0 ; /* initialize */
      array (p.wigglesR, 2 * chromMax * runMax - 1, Array) = 0 ; /* initialize */
    }

  int nDone = 0 ;
  while (channelGet (p.doneChan, &bb, BB))
    {
      long int n = (bb.hits ? bigArrayMax (bb.hits) : 0) ;
      int laneToDo ;
      int laneDone ;
      if (bb.isGenome && p.bbG.cwsN)
	{
	  saCpuStatCumulate (cpuStats, p.bbG.cpuStats) ;
	  for (int k = 0 ; k < NN ; k++)
	    {
	      if (p.bbG.cwsN) bigArrayDestroy (p.bbG.cwsN[k]) ;
	    }
	  ac_free (p.bbG.cwsN) ;
	  continue ;
	}
      if (p.debug) printf ("%s:Block done\n", timeBufShowNow (tBuf)) ;
      if (p.debug) printf ("Found %ld hits\n", n) ; 
      p.gpu += bb.gpu ;
      skips0 += bb.skips0 ;
      skips1 += bb.skips1 ;
      skips2 += bb.skips2 ;
      skips3 += bb.skips3 ;
      skips4 += bb.skips4 ;
      skipsFound += bb.skipsFound ;
      skipsNotFound += bb.skipsNotFound ;

      nHits += n ;
      nSeqs += bb.nSeqs ;
      nBaseAligned += bb.length ;
      nerr += bb.nerr ;
      nAli += bb.nAli ;
      aliDa += bb.aliDa ;
      aliDx += bb.aliDx ;
      if (bb.cpuStats)
	saCpuStatCumulate (cpuStats, bb.cpuStats) ;
      bb.runStat.nIntronSupportPlus = bb.nIntronSupportPlus ;
      bb.runStat.nIntronSupportMinus = bb.nIntronSupportMinus ;
	
      if (bb.run)
	{
	  saPolyAsCumulate (&p, &bb) ;
	  saSLsCumulate (&p, &bb) ;
	  saIntronsCumulate (&p, &bb) ;
	  saDoubleIntronsCumulate (&p, &bb) ;
	  saRunStatsCumulate (0, &p, &bb) ;
	  saRunStatsCumulate (bb.run, &p, &bb) ;

	  runErrorsCumulate (0, runErrors, bb.errors) ;
	  runErrorsCumulate (bb.run, runErrors, bb.errors) ;

	  if (p.wiggle) saWiggleCumulate (&p, &bb) ;
	}
      
      laneToDo = atomic_fetch_add (arrp (p.runLanes, bb.run, atomic_int), 0) + 0 ;
      laneDone = atomic_fetch_add (arrp (p.runLanesDone, bb.run, atomic_int), 1) + 1 ;

      if (laneToDo == laneDone)
	{
	  fprintf (stderr, "run %s done\n", dictName (p.runDict, bb.run)) ;
	}

      /* release block  memory */
      if (bb.dnas)
	{
	  Array dna = 0 ;
	  int iMax = arrayMax (bb.dnas) ;
	  for (int i = 1 ; i < iMax ; i++)
	    {  /* the dna array points into the memory globalDna array */
	      dna = arr (bb.dnas, i, Array) ;
	      if (dna && dna->lock)
		{
		  dna->base = 0 ;
		  arrayUnlock (dna) ;
		  arr (bb.dnas, i, Array) = 0 ;
		}
	    }
	}
      if (bb.quals)
	{
	  Array qual = 0 ;
	  int iMax = arrayMax (bb.quals) ;
	  for (int i = 1 ; i < iMax ; i++)
	    { 
	      qual = arr (bb.quals, i, Array) ;
	      if (arrayExists (qual))
		ac_free (qual) ;
	      arr (bb.quals, i, Array) = 0 ;
	    }
	}
      if (bb.run)
	{
	  int ns = 0 ;
	  char tBuf[25], tBuf2[25] ;
	  bb.stop = timeNow () ;
	  timeDiffSecs (bb.start, bb.stop, &ns) ;
	  fprintf (stderr, "%s: run %d / slice %d done (%d/%d)  start %s elapsed %d s, nSeqs %ld nBases %.1g strategy %d\n"
		   ,  timeBufShowNow (tBuf), bb.run, bb.lane, ++nDone, NTODO, timeShow (bb.start, tBuf2, 25), ns, bb.nSeqs, (double)bb.length, bb.isRna) ; 
	}
      ac_free (bb.h) ;
    }
  {{
      int n = sizeof(float) * (1+dictMax (p.runDict)) ;
      p.runStranding = halloc (n, p.h) ;
      memset (p.runStranding, 0, n) ;
    }}

  saIntronStranding (&p, p.confirmedIntrons) ;
  if (p.introns)
    {
      saIntronsExport (&p, p.confirmedIntrons) ; /* before wiggleExport to restrand the gene expression */ 
      saDoubleIntronsExport (&p, p.doubleIntrons) ;
    }
  if (p.wiggle)
    saWiggleExport (&p, nAgents) ;
  saCpuStatExport (&p, cpuStats) ;
  saPolyAsExport (&p, p.confirmedPolyAs) ;
  saSLsExport (&p, p.confirmedSLs) ;
  saRunStatExport (&p, p.runStats) ; /* must come afer PolyAsExport and IntronsExport */
  
  wego_log ("Done") ;
  wego_flush () ; /* flush the wego logs to stderr */

 done:
  if (1)
    {
      int ns = 0 ;
      printf ("%s: Program start\n",  tBuf0) ;
      t1 = timeNow () ;

      timeDiffSecs (t0, t1, &ns) ;
      printf ("%s: Program end elapsed %d s\n",  timeBufShowNow (tBuf), ns) ; 
    }

  printf ("\tTarget %d sequences %ld bases\n", p.bbG.dnas ? arrayMax (p.bbG.dnas) - 1 : 0, p.bbG.length) ;
  if (1 || p.debug) printf ("Skips: 0=%ld, %d=%ld, %d=%ld, %d=%ld, %d=%ld, found=%ld, notFound=%ld\n",
			    skips0, step1, skips1, step2, skips2, step3, skips3, step4, skips4, skipsFound, skipsNotFound);
  if (1 || p.debug) printf ("SeedLength %d, tStep=%d, iStep=%d, maxTargetRepeats read/target=%d/%d, nAgents=%d nBlocks=%d NN=%d\n"
			    , p.seedLength, p.tStep, p.iStep
			    , p.maxTargetRepeats, p.tMaxTargetRepeats
			    , nAgents, p.nBlocks, NN
			    ) ;
  if (p.gpu) printf ("GPU called %d times\n", p.gpu) ;
  if (arrayMax (p.runStats))
    reportRunStats (&p, p.runStats) ;
  if (p.align)
    reportRunErrors (&p, p.runStats, runErrors) ;
  if (p.bam) ac_free (p.bamHandle) ;
  /* release memory */
  if (p.bbG.dnas)
    {
      Array dna = 0 ;
      int iMax = arrayMax (p.bbG.dnas) ;
      for (int i = 1 ; i < iMax ; i++)
	{  /* the dna array points into the memory mapped globalTargetDna array */
	  dna = arr (p.bbG.dnas, i, Array) ;
	  if (dna->lock)
	    {
	      dna->base = 0 ;
	      arrayUnlock (dna) ;
	    }
	} 
      if (1) /* use valgrind if double freed errors occur */
	ac_free (p.bbG.h) ;
    }
  /* wego_log is the thread-safe way to pass messages to stderr */
  if (p.justStats && p.outFileName)  system (hprintf (h, "touch %s/toto.BF.gz ; \\rm %s/*.BF.gz %s/*.hits &", p.outFileName  , p.outFileName)) ;  
  if (0)   ac_free (h) ; /* blocks on channel cond destroy */
  return 0 ;
}

/**************************************************************/
/**************************************************************/
/**************************************************************/
