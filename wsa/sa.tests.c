/*
 * sa.tests.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module implements several tests
*/

#include "sa.h"



/* create a random genome of so many Megabases
 * add a gt_ag introns at locations n*10kb + 101-102/198-199 // exon 200 300 // ... until 1000
 * export th corresponding gtf and corresponding mrna.fasta
 */

void saCreateRandomGenome (PP *pp, int nMb)
{
  AC_HANDLE h = ac_new_handle() ;
  int iMax = 10000* nMb ;
  int ii, jj ;
  char *atgc = "atgc" ;
  ACEOUT aog = aceOutCreate (pp->outFileName, "/random.genome.fasta", 0, h) ;
  ACEOUT aom = aceOutCreate (pp->outFileName, "/random.mrna.fasta", 0, h) ;
  ACEOUT aof = aceOutCreate (pp->outFileName, "/random.gtf", 0, h) ;
  ACEOUT aoc = aceOutCreate (pp->outFileName, "/random.tConfig", 0, h) ;
  
  if (iMax < 0) messcrash ("nMb = %d too large in saCreateRandomIndex", nMb) ;

  char buf[iMax+10] ;
  char mrna1[1200] = {0} ;
  char mrna2[1200] = {0} ;

  aceOutf (aoc, "G\t%s/random.genome.fasta\n", pp->outFileName) ;
  aceOutf (aoc, "I\t%s/random.gtf\n", pp->outFileName) ;
  
  for (ii = 0 ; ii < iMax ; ii++)
    buf[ii] = atgc[randint() % 4] ;
  buf[iMax] = 0 ;

  for (ii = 0 ; ii < nMb  ; ii++)
    {
      int a1, a2, da = 0 , dda ;
      for (jj = 0 ; jj < 5 ; jj++)
	{
	  a1 = 10000 * ii + 200 * jj + 1 ;
	  dda = (jj == 2 ? 100 : 100) ;
	  a2 = a1 + dda - 1 ;
	  memcpy (mrna1+da, buf + a1 -1, dda) ;
	  da += dda ;
	  if (a1 > 2) buf[a1-3] = 'a' ;
	  if (a1 > 2) buf[a1-2] = 'g' ;
	  buf[a2] = 'g' ;
	  buf[a2+1] = 't' ;
	  aceOutf (aof, "chr1\trandom\texon\t%d\t%d\t.\t+\t0\tgene_id \"gene_%d ; transcript_id mrna_%d;\"\n", a1, a2, ii, ii) ;
	  if (jj > 0 && jj < 4)
	    aceOutf (aof, "chr1\trandom\tCDS\t%d\t%d\t.\t+\t0\tgene_id \"gene_%d ; transcript_id mrna_%d;\"\n", a1, a2, ii, ii) ;
	}
      /* add a polyA tail */
      jj = 0 ;
      if (0) for (jj = 0 ; jj < 12 ; jj++)
	mrna1[da + jj] = 'a' ;
      mrna1[da + jj] = 0 ;
      /* duplicate the reads with independant errors */
      memcpy (mrna2, mrna1, da+jj+1) ;
      int kk = 5, i, k, dx1, dx2 ;
      /* introduce kk substitutions */
      if (1 )
	for (k = 0 ; k < kk ; k++)
	  {
	    int x1 = randint() % da ;
	    mrna1[x1] = 'c' ;
	    int x2 = randint() % da ;
	    mrna2[x2] = 'g' ;
	  }
      /* introduce kk deletions */
      if (0)
	for (i = dx1 = dx2 = 0 ; i <= da ; i++)  /* include the therminal zero: i <= da */
	  {
	    int x1 = randint() % 100 ;
	    if (x1 ==0)  dx1++ ;
	    mrna1[i] = mrna1[i + dx1] ;
	    int x2 = randint() % 100 ;
	    if (x2 ==0)  dx2++ ;
	    mrna2[i] = mrna2[i + dx2] ;
	  }
      /* introduce kk insertions */
      if (0 )
	for (i = dx1 = dx2 = 0 ; i <= da ; i++)  /* include the therminal zero: i <= da */
	  {
	    int x1 = randint() % 100 ;
	    if (x1 ==0)  dx1++ ;
	    mrna1[i+dx1] = mrna1[i] ;
	    int x2 = randint() % 100 ;
	    if (x2 ==0)  dx2++ ;
	    mrna2[i+dx2] = mrna2[i] ;
	  }
      
      aceOutf (aom, ">s%d\n%s\n", ii, mrna1) ;
      aceOutf (aom, ">t%d\n%s\n", ii, mrna2) ;
    }
  aceOutf (aog, ">chr1\n%s\n", buf) ; /* export the genome after tweaking the gt_ag */
  ac_free (h) ;
  exit (0) ;
}
