/*  File:  dnabloom.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2015
 * -------------------------------------------------------------------
 * AceView is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the AceView project developped by
 *	 D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *
 *  This should be linked with the acedb libraries
 *  available from http://www.acedb.org
 *
 *  Please direct all question about this code to
 *  mieg@ncbi.nlm.nih.gov
 *
 * Purpose of this code
 *  Phase 1:
 *   Take a genome, split it into N sections, or read a table of coordinates defining the sections
 *   so that the section contain an approximately equal number of genes and do not split a gene.
 *   For each section, create a bloom filter of words of length L.
 *   Save these bitsets in binary format.
 *   To allow both strands, we store the dna if the centtal base is A or G, otherwise the complement.
 *
 * Phase 2:
 *   Stream a collection of sequences in fasta or fastc format
 *   For each sequence, count the number of Bloom hits to the N sections.
 *   Assign the sequence to the best or to several best sections.
 *   Ventilate out the sequences in N fasta/fastc files.
 *   The ambiguous sequences are exported in several files and marked
 *   so that they can later be collected and treated and merged using bestali.
 
 * One idea is to write a parallel program for phase 2
 * Another idea is to optimize the hash function
 * If we use 13-mers and 256 sections, each section would be length 10M, 4^13 = 2^26 = 64M
 * but since we select the strand we economize a factore 2->32M words
 * we can stote these a a biSet of size 4Mbytes = 32 Mbits and use the 13mer as an address

 * pseudo code
 parse
 new class -> new section (so that mito, ercc, ... are separated)
 new sequence, add it il old length + new length < 10M, otherwise start new section
 if (sequence > 15M split at 10M and iterate
 create a bitset of size 32Mbits and scan the section in a parallel thread, then export the bitset
 
 ventilate
 grab K sequences and create sets of N counters for each sequence analyse the set in a new thread
 and export that set via a new pipe
 
 usage ...
 main ...

 * in 2108, we tested the code on the NASA data looking for telomeres
 * the -telomere code paralelizes well: on lmem12 n_threads=32, compression 2800%
*/

/* #define ARRAY_CHECK  */

/* 13:5 means take a 13 mer every 6 bases so a 18 mer counts 1, a 19mer counts 2 */
#define WORDLENGTH 13
#define WORDJUMP 5

#define LATESORT 0
#include "ac.h"
#include "bitset.h"
#include "channel.h"
#include "wego.h"

#define IBMAX 10
/* analyse the telomeres inside pices of average length SEGLENGTH */
#define SEGLENGTH 1000   
#define SEGLENGTH2 2000  /* 2*SEGLENGTH */

typedef struct pairStruct { char *buf ; } SHORT_READ ; /* nam:forward seq:reverse seq */
typedef struct seqStruct { 
  int nam ; /* the name of the sequence: defined in the zone dict */
  int ln, x1, x2 ; /* total length of the sequence, the cordicnates of the fragment contained in the zone */
  int pos ; /* offset of the first base of this fragement in the zone stack */
} SEQ ;

typedef struct zoneStruct { 
  DICT *dict ; /* dictionary of the anmes of the sequences containes in this zone */
  Stack dna ;  /* stack containing the DNA of these fragments */
  Array seqs ; /* Array of SEQ belonging to this fragment */
  BitSet bb ;  /* the bitset of all WORDLENGTH  mers present in this fragment */
  CHAN *cin ;
}  ZONE ;

typedef struct teloExportStruct { int ib, nn, mult ; char dna[SEGLENGTH2] ; } TELOEXPORT ; 
typedef enum { FASTA=0, FASTC, FASTQ, CSFASTA, CSFASTC, RAW, GLOBAL } DNAFORMAT ;

typedef struct  pStruct {
  AC_HANDLE h ;
  ACEOUT ao ;
  AC_HANDLE seqH ; /* then names and sequence stack of the zone are allocated on seqH, so they can be released first */
  const char *inFileName, *outFileName, *targetFileName ;
  const char *run ;
  char *titles[IBMAX] ;
  BOOL gzi, gzo, silent, noPolyA ;
  BOOL telomere ;
  int max_threads, ibMax, block ;
  int zoneSize ; /* default 10M */
  int maxZone ;  /* total number of zones */
  Array zones ;
  DNAFORMAT dnaFormat ;
  CHAN *genomeChan ;
  CHAN *teloExportChan ;
  CHAN *teloExportDoneChan ;
  int wordJumps[IBMAX] ;
} PP ;

/*************************************************************************************/
/* analyse a sequence
 *   pass 0: creation : fill a Bloom BitSet, return the number words analysed
 *   pass 1 :counting : return the number of hits to the BitSet
 */
typedef struct int4Struct { int xx[IBMAX] ; } I4 ;

#define SHIFT  (2 * (WORDLENGTH - 1))
#define MASK  ((1 << (2 * WORDLENGTH)) - 1)
#define CENTRALMASK (1 << WORDLENGTH)
#define RIGHTMASK (CENTRALMASK - 1)
#define LEFTMASK ((RIGHTMASK << WORDLENGTH) & MASK)

static I4 dnaBloomOneSequence (int phase, const BitSet blooms[], int wordJumps[], const char *dna, int *nks, int *polyAp, int ibMax)
{
  I4 mynks ;
  register unsigned int f = 0 ; /* the word buffer in the forward direction */
  register unsigned int r = 0 ; /* its reversecomplement r */
  register unsigned int g ;     /* the masked compressed selected word */
  register unsigned int x ;  
  register int n = 0, polyA = 0, ib ;
  register const char *cp = dna - 1 ;
  int jump[IBMAX], jumpA = 0 ;

  memset (&mynks, 0, sizeof(mynks)) ;
  memset (jump, 0, ibMax * sizeof(int)) ;

  while (1)
    {
      cp++ ;
      n++ ;
      switch (*cp)
	{
	case 'a': case 'A': x = 0x0 ; break ;
	case 't': case 'T': x = 0x3 ; break ;
	case 'g': case 'G': x = 0x1 ; break ;
	case 'c': case 'C': x = 0x2 ; break ;
	case ' ' : continue ;
	case '\n' : continue ;  /* merge bases across line breaks */
	case 0:  n = 0 ; x = 0x8 ; break ;
        default: x = n = 0 ; break ;
	}
      f = (x | (f << 2)) & MASK ;
      r = (( (x ^ 0x3) << SHIFT) | (r >> 2)) & MASK ;
      g = (f & CENTRALMASK ? r : f) ;                 /* select the strand */
      g = (g & RIGHTMASK) | ((g >> 1) & LEFTMASK) ;   /* gobble the half base used to select the strand */
      
      if (n >= WORDLENGTH)
	{
	  if (jumpA-- <= 0 && g == 0) 
	    { polyA ++ ; jumpA = WORDJUMP ;}
	  
	  for (ib = 0 ; ib < ibMax ; ib++) 
	    {  
	      BitSet bloom = blooms[ib] ;  
	      switch (phase) 
		{
		case 0:
		  bitSet (bloom, g) ;
		  break ;
		case 1:
		  if (jump[ib]-- <= 0 && bit (bloom, g)) 
		    { 
		      mynks.xx[ib]++ ; jump[ib] = wordJumps[ib] ;
		    } 
		  break ;
		}
	    }
	}
      if (x == 0x8) break ;
    }  
  
  if (polyAp)
    *polyAp = polyA ;
  return mynks ;
} /* dnaBloomOneSequence */

/*************************************************************************************/

static void dblExportTelomere (const void *vp)
{
  const PP *pp = vp ;
  AC_HANDLE h = ac_new_handle() ;
  TELOEXPORT t ;
  int n = 0 ;
  ACEOUT ao = aceOutCreate (pp->outFileName, ".telomere_hits.txt", pp->gzo, h) ; 
 
  aceOutDate (ao, "###", "exact and approximate telomere motifs") ;
  aceOutf (ao, "#Run\tType\tNb motifs\tRead multiplicity\tSequence\n") ;
  while (channelGet (pp->teloExportChan, &t, TELOEXPORT))
    aceOutf (ao, "%s\t%s\t%d\t%d\t%s\n", pp->run, pp->titles[t.ib], t.nn, t.mult, t.dna) ;
  channelPut (pp->teloExportDoneChan, &n, int) ;
  ac_free (h) ;
}

/*************************************************************************************/

typedef struct telomereStruct { 
  int thId ; 
  int nFullSeqs, nFullTags ;
  int nPieceSeqs, nPieceTags ;
  int ibMax ;
  char *buf ; 
  PP *pp ;
  CHAN *chan, *done ; 
  KEYSET kmerFullHistos[IBMAX], kmerPieceHistos[IBMAX], polyaPieceHisto, polyaFullHisto ; 
  BitSet bb, blooms[IBMAX] ; } TELO ;

static void dblDoCountTelomere (const void *vp)
{
  const TELO *tp = (TELO *)vp ;
  int k, nLines = 0, mult = 0, polyA = 0 ;
  int ib, ibMax = tp->ibMax ;
  int nks[ibMax] ;
  char *cp, *cq ;
  PP *pp = tp->pp ;
  BOOL debug = FALSE ;
  DNAFORMAT dnaFormat = tp->pp->dnaFormat ;
  memset (nks, 0, sizeof(nks)) ;
  
  while (1)
    {
       if (tp->pp->max_threads)
	 {
	   if (debug) fprintf (stderr, "dblDoCountTelomere waits on channelGive (chan-%d) nFullTags=%d nPieceTags=%d\n", tp->thId, tp->nFullTags, tp->nPieceTags) ;
	   k = channelGive (tp->chan, &k, int) ; /* wait for data to be ready */
	   if (k == -1)    /* work is over */
	     {
	       if (debug) 
		 fprintf (stderr, "dblDoCountTelomere channel %d got k == -1 and calls channelPut(done)\n", tp->thId) ;
	       k = - tp->thId - 1 ;
	       channelPut (tp->done, &k, int) ; 
	       return ;
	     }
	 }
      cp = cq = tp->buf ;
      mult = 1 ;
      while (cp && *cp)
	{
	  nLines++ ;
	  switch (dnaFormat)
	    {
	    case GLOBAL:
	      cq = 0 ;  /* accept the whole buffer at once */ 
	      break ; 
	    case FASTQ:
	      cq = strchr (cp, '\n') ;
	      if (cq) *cq++ = 0 ;
	      if ((nLines % 4) != 1) /* only keep the DNA lanes */
		cp = 0 ; 
	      break ;
	    case RAW:
	      cq = strchr (cp, '\n') ;
	      if (cq) *cq++ = 0 ;
	      break ;
	    case FASTC:
	    case FASTA:
	    default:
	      switch ((int)*cp)
		{
		case 0: 
		case '#': 
		  cq = strchr (cp, '\n') ;
		  cp = 0 ; /* jump this line */
		  break ;
		case '>':
		  cq = strchr (cp, '\n') ;
		  if (cq) *cq = 0 ;
		  mult = fastcMultiplicity (cp+1, 0, 0) ;
		  cp = 0 ; /* jump this line */
		  break ;
		default:
		  if (dnaFormat == FASTA) /* merge successive lines */
		    cq = strstr (cp, "\n>") ;
		  else
		    cq = strchr (cp, '\n') ;
		  break ;
		} 
	    }
	  if (cq) { *cq++ = 0 ; }
	  
	  if (cp)
	    {
	      I4 nks ; 
	      int i, ln = strlen(cp), iMax = (ln + (SEGLENGTH/2 - 1))/SEGLENGTH ;
	      char cc, *cr = cp, *cs ;
              int delta ;

	      if (iMax < 1)
		iMax = 1 ;
	      delta = ln/iMax ;
	      for (cr = cp, i = 0 ; i < iMax ; cr += delta, i++)
		{
		  cs = 0 ;
		  if (i < iMax - 1)
		    {
		      cs = cr + delta ;
		      cc = *cs ;
		      *cs = 0 ;
		    }
		
		  nks = dnaBloomOneSequence (1, tp->blooms, tp->pp->wordJumps, cr, 0, &polyA, ibMax) ;
		  for (ib = 0 ; ib < ibMax ; ib++)
		    {
		      if (iMax == 1)
			{
			  keySet (tp->kmerFullHistos[ib], nks.xx[ib]) += mult ;
			}
		      keySet (tp->kmerPieceHistos[ib], nks.xx[ib]) += mult ;
		      if (nks.xx[ib] >= 3)
			{
			  TELOEXPORT t ;

			  t.ib = ib ;
			  t.nn = nks.xx[ib] ;
			  t.mult = mult ;
			  memcpy (t.dna, cr, strlen(cr) + 1) ;
			  channelPut (pp->teloExportChan, &t, TELOEXPORT) ;
			}
		    }
		  ((TELO*)tp)->nPieceSeqs++ ; /* analyzing */
		  ((TELO*)tp)->nPieceTags += mult ; 

		  if (iMax == 1 && ! tp->pp->noPolyA)
		    keySet (tp->polyaFullHisto, polyA) += mult ;
		  if (! tp->pp->noPolyA)
		    keySet (tp->polyaPieceHisto, polyA) += mult ;
		  if (cs)
		    *cs = cc ;
		}
	      if (iMax > 1)
		{  /* search again on the full sequence */
		  nks = dnaBloomOneSequence (1, tp->blooms, tp->pp->wordJumps, cp, 0, &polyA, ibMax) ;
		  for (ib = 0 ; ib < ibMax ; ib++)
		    keySet (tp->kmerFullHistos[ib], nks.xx[ib]) += mult ;
		  if (! tp->pp->noPolyA)
		    keySet (tp->polyaFullHisto, polyA) += mult ;
		}
	      ((TELO *)tp)->nFullSeqs++ ; /* analizing */
	      ((TELO *)tp)->nFullTags += mult ; 
	    }
	  cp = cq ; /* start of next line */
	}
      k = tp->thId ;
      if (tp->pp->max_threads)
	{
	  if (debug) fprintf (stderr, "dblDoCountTelomere channel %d and calls channelPut(done)\n", k) ;
	  if (debug) fprintf (stderr, "dblDoCountTelomere offset %d , run %s parse %d lines, analyzed %d fragments, %d distinct %d pieces, %d distinct pieces, %d polyA : %s\n"
			      , tp->thId,  tp->pp->run, nLines, tp->nFullTags, tp->nFullSeqs, tp->nPieceTags, tp->nPieceSeqs, polyA, timeShowNow()) ;
	  channelPut (tp->done, &k, int) ;  /* this thread is ready for a new block */
	}
      else
	break ;
    }
} /* dblDoCountTelomere */

/*************************************************************************************/
/* read a big number of lines, and regularizes to the start of a sequence uing the carryOver */
static int dblFillBuffer (int fd, char *buf, int size, char *carryOverBuf, int carryOver)
{
  int n, n1 = 0, nn = size ;
  char *cp ;
  BOOL debug = FALSE ;

   if (debug) fprintf (stderr, "dblFillBuffercalls read\n") ;
  if (carryOver)
    {
      memcpy (buf, carryOverBuf, carryOver) ;
      nn -= carryOver ;
    } 
   if (debug) fprintf (stderr, "dblFillBuffercalls read\n") ;
  while (nn > 0)   /* read from the file up to size bytes */
    {
      n = read (fd, buf + carryOver + n1, nn) ;
      if (n < 0)
	messcrash ("Error in dblFillBuf:read") ;
      if (debug)  fprintf (stderr, "-----dblFillBuffer requested %d bytes got %d bytes cumul %d bytes\n", nn, n, n1);
      nn -= n ; n1 += n ;
      if (n == 0)
	break ;
    }
   if (debug) fprintf (stderr, "dblFillBuffer got %d bytes\n", n1) ;
  cp = buf + size - nn ;
  cp[0] = 0 ; cp[1] = 0 ;  /* double zero terminate */
  if (nn == 0)
    {
      char *cp0 = cp ;
      while (cp > buf)  /* check for a reminder */
	if (! strncmp (--cp, "\n>", 2))
	  break ;
      if (cp <= buf)   /* this file does not contain a '>' fasta delimiter */
	cp = cp0 ; 
      *cp++ = 0 ;  /* split out the remainder */
      carryOver = strlen (cp) ;
      if (carryOver)
	memcpy (carryOverBuf, cp, carryOver) ;
    }
  else
    {
      cp[0] = 0 ; cp[1] = 0 ; 
      carryOver = 0 ;
    }
  return carryOver ;
}  /* dblFillBuffer */

/*************************************************************************************/
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

static void dblParseFastaFile (PP *pp, CHAN *done, int maxTh, int size, Array telos)
{
  char *carryOverBuf = halloc (1 << 21, 0) ; /* 2 Mbytes */
  int k = 0, n, fd, carryOver = 0, nn = 0 ;
  TELO *tt ;
  BOOL debug = FALSE ;

  if (pp->inFileName)
    fd = open (pp->inFileName, O_RDONLY) ;
  else
    fd = 0 ;
  if (fd < 0)
    messcrash ("dblParseFastaFile cannot open file %s", pp->inFileName) ;
  while (1)
    {
      /* listen on the communication to find an available thread */
      if (pp->max_threads)
	{
	  if (nn < maxTh)
	    k = nn++ ;
	  else
	    k = channelGive (done, &k, int) ; /* the number of the channel which is done */
	}
      tt = arrayp (telos, k, TELO) ; 
      if (debug) fprintf (stderr, "dblParseFastaFile calls dblFillBuffer k=%d\n", k) ;
      carryOver = dblFillBuffer (fd, tt->buf, size -  (1 << 21),  carryOverBuf, carryOver) ;
      if (tt->buf[0] == 0)
	break ;  
      n = 1 ;
      if (pp->max_threads) /* signal that data is ready for thread tt */
	{
	  if (debug) fprintf (stderr, "dblParseFastaFile prepared data and calls channelPut(chan-%d)",k) ;
	  channelPut (tt->chan, &n, int) ; /* data ready */
	}
      else  /* direct procedural call */
	dblDoCountTelomere (tt) ;
    }

  if (pp->max_threads) /* synchronize */
    {
      for (k = 0 ; k < maxTh ; k++)
	{ 
	  tt = arrayp (telos, k, TELO) ;
	  if (debug) fprintf (stderr, "dblParseFastaFile is over %d tags %d pieces in chan-%d\n", tt->nFullTags, tt->nPieceTags, k) ;
	  n = -1 ;
	  if (debug) fprintf (stderr, "dblParseFastaFile is over and calls channelPut(chan-%d)",k) ;
	  channelPut (tt->chan, &n, int) ; /* close the thread */
	}
      n = 0 ; k = maxTh ;
      while (k > 0)
	{
	  if (debug) fprintf (stderr, "dblParseFastaFile waits on channelGive(done)") ;
	  n = channelGive (done, &n, int) ;
	  if (n < 0)  /* count the work done signals */
	    {
	      k-- ;    
	      if (debug) fprintf (stderr, "thread %d done %s\n", -n, timeShowNow()) ;
	    }
	}
    }

  ac_free (carryOverBuf) ;
  return ;
} /* dblParseFastaFile */

/*************************************************************************************/

static char *getMito (PP *pp)
{
  int nn, nn1 = 1<<27 ;
  char *buf = halloc (nn1, pp->h) ;
  int fd = open (pp->targetFileName, O_RDONLY) ;

  nn = read (fd, buf, nn1 - 2) ;
  buf[nn] = 0 ;

  close (fd) ;
  return buf ;
}  /* getMito */

/*************************************************************************************/

static int dblCountTelomeres (PP *pp)
{
  AC_HANDLE h  = ac_new_handle () ;
  int ib, ibMax = pp->ibMax ;
  int ii, j, nn = 0, maxTh ;
  int nFullSeqs = 0, nFullTags = 0 ;
  int nPieceSeqs = 0, nPieceTags = 0 ;
  CHAN *done  = channelCreate (24, int, h) ;
  KEYSET kmerFullHistos[ibMax], kmerPieceHistos[ibMax], polyaFullHisto , polyaPieceHisto ;
  TELO *tp ;
  Array telos = arrayHandleCreate (pp->max_threads + 1, TELO, h) ;
  BitSet blooms[ibMax] ;
  const char *telomereT = "ttagggttagggttagggttagggttagggttagggttagggttagggttagggttaggg" ;
  /* all doublets tNagggtNaggg */
  const char *telomereC = "ttagggttagggttagggtcagggttagggttagggttagggtcagggttagggtaagggttagggtgagggttagggtcagggttagggttagggtcagggtaagggtcagggtgagggtcagggtcagggttagggttagggttagggttagggtgagggtaagggtgagggtgagggtgagggtcagggttagggttagggttagggttagggtaagggtaagggtaagggtgagggtaagggtcagggttagggttaggg" ;
  const char *imagT = "AATCCCAATCCCAATCCCAATCCCAATCCCAATCCCAATCCCAATCCCAATCCCAATCCC" ;
  const char *imagC = "AATCCCAATCCCAATCCCAGTCCCAATCCCAATCCCAATCCCAGTCCCAATCCCATTCCCAATCCCACTCCCAATCCCAGTCCCAATCCCAATCCCAGTCCCATTCCCAGTCCCACTCCCAGTCCCAGTCCCAATCCCAATCCCAATCCCAATCCCACTCCCATTCCCACTCCCACTCCCACTCCCAGTCCCAATCCCAATCCCAATCCCAATCCCATTCCCATTCCCATTCCCACTCCCATTCCCAGTCCCAATCCCAATCCC" ;

  int size = (1 << 22) ; /* 4M */ 
  BOOL debug = FALSE ;
 
  if (pp->block * (1 << 20) > size)
    size = pp->block * (1 << 20) ;
  if (ibMax > IBMAX)
    messcrash ("ibMax = %d > IBMAX = %d", ibMax, IBMAX) ;
  if (debug) channelDebug (done, TRUE, "chan-done") ;
  pp->teloExportChan = channelCreate (256, TELOEXPORT, h) ;	
  
  polyaPieceHisto = keySetHandleCreate (h) ; /* global histogram */
  polyaFullHisto = keySetHandleCreate (h) ; /* global histogram */
     
  for (ib = 0 ; ib < ibMax ; ib++)
    { 
      const char *t = 0 ;
      kmerFullHistos[ib] = keySetHandleCreate (h) ; /* global histogram */
      kmerPieceHistos[ib] = keySetHandleCreate (h) ; /* global histogram */
      blooms[ib] = bitSetCreate (1 << (2*WORDLENGTH - 1), h) ;
      switch(ib)
	{
	case 0: t = telomereC ; pp->titles[ib] = "Tel_C_6" ; pp->wordJumps[ib] = 5 ;  break ; 
	case 1: t = telomereT ; pp->titles[ib] = "Telomeric_6" ; pp->wordJumps[ib] = 5 ; break ; 
	case 2: t = imagC ; pp->titles[ib] = "imagC_6" ; pp->wordJumps[ib] = 5 ; break ; 
	case 3: t = imagT ; pp->titles[ib] = "imagT_6" ; pp->wordJumps[ib] = 5 ; break ; 
	default: 
	  if (pp->targetFileName)
	    { t = getMito(pp) ;  pp->titles[ib] = "extern" ; }
	  break ;
	} 

      dnaBloomOneSequence (0, &blooms[ib], 0, t, 0, 0, 1) ;
	fprintf (stderr , "## ib = %d registered %lu bits int the bloom\n", ib, bitSetCount (blooms[ib])) ;
    }
  fprintf (stderr, "--- data preparation done : %s\n", timeShowNow ()) ;
  wego_go (dblExportTelomere, pp, PP) ; 
      

  maxTh = pp->max_threads ;
  if (maxTh < 1) maxTh = 1 ;
  for (ii = 0 ; ii < maxTh ; ii++)
    {  
      /* allocate all elements */
      tp = arrayp (telos, ii, TELO) ;
      for (ib = 0 ; ib < ibMax ; ib++)
	{
	  tp->kmerFullHistos[ib] = keySetHandleCreate (h) ;
	  tp->kmerPieceHistos[ib] = keySetHandleCreate (h) ;
	  tp->blooms[ib] = blooms[ib] ;  /* the blooms are read only once and common to all threads */
	}
      
      tp->bb =  bitSetCreate (1 << (2 * WORDLENGTH - 1), h) ;
      tp->polyaPieceHisto = keySetHandleCreate (h) ;
      tp->polyaFullHisto = keySetHandleCreate (h) ;
      tp->buf = halloc (size, h) ;

      if (pp->max_threads)
	{
	  tp->chan = channelCreate (24, int, h) ;
	  if (debug) channelDebug (tp->chan, TRUE, messprintf ("chan-%d", ii)) ;
	}
      /* select the slicing */
      tp->thId = ii ; 
      tp->done = done ;
      tp->pp = pp ;
      tp->ibMax = ibMax ;
      
      /* start the thread */
      if (pp->max_threads)
	wego_go (dblDoCountTelomere, tp, TELO) ;
    }
  
  dblParseFastaFile (pp, done, maxTh, size, telos) ;

  channelClose (pp->teloExportChan) ;

  /* cumulate the individual histos */
  for (ii = 0 ; ii < maxTh ; ii++)
    { 
      tp = arrayp (telos, ii, TELO) ;
      for (ib = 0 ; ib < ibMax ; ib++)
	{
	  for (j = keySetMax(tp->kmerFullHistos[ib]) - 1 ; j >= 0 ; j--)
	    keySet (kmerFullHistos[ib], j) +=  keySet(tp->kmerFullHistos[ib],j) ;
	  for (j = keySetMax(tp->kmerPieceHistos[ib]) - 1 ; j >= 0 ; j--)
	    keySet (kmerPieceHistos[ib], j) +=  keySet(tp->kmerPieceHistos[ib],j) ;
	}
      for (j = keySetMax(tp->polyaPieceHisto) - 1 ; j >= 0 ; j--)
	keySet (polyaPieceHisto, j) +=  keySet(tp->polyaPieceHisto,j) ;
      for (j = keySetMax(tp->polyaFullHisto) - 1 ; j >= 0 ; j--)
	keySet (polyaFullHisto, j) +=  keySet(tp->polyaFullHisto,j) ;
      
      nFullTags += tp->nFullTags ;
      nFullSeqs += tp->nFullSeqs ;
      nPieceTags += tp->nPieceTags ;
      nPieceSeqs += tp->nPieceSeqs ;
    }
  fprintf (stderr, "--- histo cumulated %s\n", timeShowNow()) ;
  /* export */
  if (1)
    {
      KEYSET ks ;
      int n, cumulFullK[ibMax], cumulPieceK[ibMax], cumulFullpA = 0,  cumulPiecepA = 0, jMax, jGlobalMax = 0 ;
      ACEOUT ao = aceOutCreate (pp->outFileName, ".telomere_distribution.txt", pp->gzo, h) ;

      aceOutDate (ao, "###", "Histogram of the number of fragments with n 13-mers matching the canonical telomere motif ") ;
      aceOutf (ao, "## Found %d fragments, %d distinct fragments\n", nFullTags, nFullSeqs) ;
      aceOutf (ao, "## Found %d pieces, %d distinct pieces\n", nPieceTags, nPieceSeqs) ;
      aceOutf (ao, "#Number of motifs\tNumber of fragments with x TERRA 13mers TTAGGGTTAGGGT\tNumber of fragments with x homopolymers of 13 consecutive A\tCumul TERRA\tCumul pA\n") ;

      jGlobalMax = jMax = keySetMax(polyaFullHisto) ; 
      for (n = j = 0 ; j < jMax ; j++) 
	n += keySet (polyaFullHisto, j) ;
      cumulFullpA  = n ;

      jGlobalMax = jMax = keySetMax(polyaPieceHisto) ; 
      for (n = j = 0 ; j < jMax ; j++) 
	n += keySet (polyaPieceHisto, j) ;
      cumulPiecepA  = n ;
      
      for (ib = 0 ; ib < ibMax ; ib++)
	{
	  ks = kmerFullHistos[ib] ;
	  jMax = keySetMax(ks) ;
	  if (jGlobalMax < jMax) 
	    jGlobalMax = jMax ;
	  for (n = j = 0 ; j < jMax ; j++) 
	    n += keySet (ks, j) ;
	  cumulFullK[ib] = n ;
	}
      for (ib = 0 ; ib < ibMax ; ib++)
	{
	  ks = kmerPieceHistos[ib] ;
	  jMax = keySetMax(ks) ;
	  if (jGlobalMax < jMax) 
	    jGlobalMax = jMax ;
	  for (n = j = 0 ; j < jMax ; j++) 
	    n += keySet (ks, j) ;
	  cumulPieceK[ib] = n ;
	}

      aceOutf (ao, "\n%s\tA%d", pp->run, WORDLENGTH) ; for (j = 0 ; j < jGlobalMax ; j++) { n =  keySet (polyaFullHisto, j) ;aceOutf (ao, "\t%d", n) ; }
      aceOutf (ao, "\n%s\tA%d-per_read_or_read_pair-cumul", pp->run, WORDLENGTH ) ; for (n = j = 0 ; j < jGlobalMax ; j++) { aceOutf (ao, "\t%d",  cumulFullpA - n) ; n += keySet (polyaFullHisto, j) ; } 

      for (ib = 0 ; ib < ibMax ; ib++)
	{
	  ks = kmerFullHistos[ib] ;
	  aceOutf (ao, "\n%s\t%s-per_read_or_read_pair", pp->run, pp->titles[ib]) ; for (n = j = 0 ; j < jGlobalMax ; j++) { n = keySet (ks, j) ; aceOutf (ao, "\t%d", n) ; }
	  aceOutf (ao, "\n%s\t%s-per_read_or_read_pair-cumul", pp->run, pp->titles[ib]) ; for (n = j = 0 ; j < jGlobalMax ; j++) { aceOutf (ao, "\t%d", cumulFullK[ib] - n) ; n += keySet (ks, j) ; }
	}
      
      aceOutf (ao, "\n") ;
 
      aceOutf (ao, "\n%s\tA%d", pp->run, WORDLENGTH) ; for (j = 0 ; j < jGlobalMax ; j++) { n =  keySet (polyaPieceHisto, j) ;aceOutf (ao, "\t%d", n) ; }
      aceOutf (ao, "\n%s\tA%d-per_1_kb_fragment-cumul", pp->run, WORDLENGTH ) ; for (n = j = 0 ; j < jGlobalMax ; j++) { aceOutf (ao, "\t%d",  cumulPiecepA - n) ; n += keySet (polyaPieceHisto, j) ; } 

      for (ib = 0 ; ib < ibMax ; ib++)
	{
	  ks = kmerPieceHistos[ib] ;
	  aceOutf (ao, "\n%s\t%s-per_1_kb_fragment", pp->run, pp->titles[ib]) ; for (n = j = 0 ; j < jGlobalMax ; j++) { n = keySet (ks, j) ; aceOutf (ao, "\t%d", n) ; }
	  aceOutf (ao, "\n%s\t%s-per_1_kb_fragment-cumul", pp->run, pp->titles[ib]) ; for (n = j = 0 ; j < jGlobalMax ; j++) { aceOutf (ao, "\t%d", cumulPieceK[ib] - n) ; n += keySet (ks, j) ; }
	}
      
      aceOutf (ao, "\n") ;

    }
  ac_free (h) ;	
  return nn ;
} /* dblCountTelomeres */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: dnabloom  -telomere [-run run_name] [-help]... \n"
	   "//      try: -h -help --help \n"
	   "// Example:  dnabloom -i f.fastq -I fastq -run xx -o out -telomere\n"
	   "// -run runName : the name of the run to be exported in the results\n"
	   "// -telomere : count the polyA and the TERRA motifs unsing appropriate k-mers\n"
	   "// -i fileName : the name of a sequence file, possibly .gz, to be analyzed \n"
	   "// -I input_file_format : one of\n"
	   "//    -I fasta: [default] the input is a fasta file\n"
	   "//    -I fastq: the input is a standard fastq file\n"
	   "//    -I fastc: fastc file (contains #multipliers for repeated sequences)\n"
	   "//    -I csfasta: SOLiD color coded transition fasta file\n"
	   "//    -I csfastc: idem (contains #multipliers for repeated sequences)\n"
	   "//    -I raw: the input is a raw file, one sequenc per line, no identifers\n"
	   "// -silent : suppress stderr status report\n"
	   "// -gzi : force gunzipping the input (useful when piping from stdin)\n"
	   "// -gzo : the output file will be gziped\n"
	   "// -max_threads <integer>: [default 1] try 4, 8..., machine dependant\n"
	   ) ;

  
  if (argc > 1)
    {
      fprintf (stderr,
	       "//\n// You said: %s\n", commandBuf) ;
      fprintf (stderr,
	       "// ########## ERROR: I do not understand the argument%s ", argc > 2 ? "s" : "") ;
      for (i = 1 ; i < argc ; i++)
	fprintf (stderr, "%s ", argv[i]) ;
      fprintf (stderr, "\n") ;
    }
  exit (1) ;
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  PP p ;
  AC_HANDLE h = 0 ;
  char commandBuf [1000] ;

  freeinit () ; 
  messErrorInit (argv[0]) ;
  
  h = ac_new_handle () ;
  memset (&p, 0, sizeof (PP)) ;
 
  p.h = h ;
  
  if (argc < 2)
    usage (commandBuf, argc, argv) ;
  if (getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)||
      getCmdLineOption (&argc, argv, "-h", 0)
      )
    usage (commandBuf, 1, argv) ;
  
  /* parse the arguments */
  getCmdLineOption (&argc, argv, "-o", &(p.outFileName)) ;
  getCmdLineOption (&argc, argv, "-i", &(p.inFileName)) ;
  /* the target file opton is not used and not advertised in the online help */
  getCmdLineOption (&argc, argv, "-targetFile", &(p.targetFileName)) ;

  if (1)
    {
      p.dnaFormat = FASTC ; /* default */
      const char *ccp = 0 ;
      if (getCmdLineOption (&argc, argv, "-I", &ccp))
	{
	  if (!strcasecmp (ccp, "global"))
	    p.dnaFormat = GLOBAL ;
	  else if (!strcasecmp (ccp, "raw"))
	    p.dnaFormat = RAW ;
	  else if (!strcasecmp (ccp, "FASTC"))
	    p.dnaFormat = FASTC ;
	  else if (!strcasecmp (ccp, "FASTA"))
	    p.dnaFormat = FASTA ;
	  else if (!strcasecmp (ccp, "FASTQ"))
	    p.dnaFormat = FASTQ ;
	  else
	    messcrash ("Sorry, unknown format in option -I %s, try dnabloom -h", ccp) ;
	}
    }

  p.run = "-"  ; /* default */
  getCmdLineOption (&argc, argv, "-run", &(p.run)) ;
  
  p.gzi = getCmdLineOption (&argc, argv, "-gzi", 0);
  p.gzo = getCmdLineOption (&argc, argv, "-gzo", 0);

  p.noPolyA = getCmdLineOption (&argc, argv, "-noPolyA", 0);
 
  p.silent = getCmdLineOption (&argc, argv, "-silent", 0);
  p.telomere = getCmdLineOption (&argc, argv, "-telomere", 0);

  p.max_threads = 0 ;
  getCmdLineInt (&argc, argv, "-max_threads", &p.max_threads);
  p.ibMax = 4 ;
  getCmdLineInt (&argc, argv, "-ibMax", &p.ibMax);
  getCmdLineInt (&argc, argv, "-block", &p.block) ;
 
  if (p.ibMax > IBMAX)
    messcrash ("Please -ibMax %d cannot exceed IBMAX = %d hard defined in the source code", p.ibMax, IBMAX) ;
  p.zoneSize = 10 ; /* default 10 Mb */
  getCmdLineInt (&argc, argv, "-zoneSize", &p.zoneSize) ;
  p.maxZone = 3500/p.zoneSize ; /* lerger than the genome im MB, so we wont reallocate the arrays */
  p.zoneSize *= 1000000 ;

  if (argc != 1)
    {
      fprintf (stderr, "unknown argument, sorry\n") ;
      usage (commandBuf, argc, argv) ;
    }
  /* check the absolute args */
  if ((WORDLENGTH & 0x1) == 0)
    messcrash ("WORDLENGTH == %d should be odd",  WORDLENGTH) ;
  if (2 *  WORDLENGTH + 1 > 8 * sizeof (int))
    messcrash ("WORDLENGTH = %d should such that, in bits, 2 *  WORDLENGTH + 1 < sizeof(int) = %d"
	       , WORDLENGTH,  8 * sizeof (int)
	       ) ;

  fprintf (stderr, "// %s start\n", timeShowNow()) ;
  aceInWaitAndRetry (0) ; /* does nothink, but triggers the linker on LINUX */
  if (p.max_threads)
    {
      wego_max_threads (p.max_threads + 2) ;
    }

  /* create the zones */
  if (0 && p.targetFileName)
    {
      AC_HANDLE h1 = ac_new_handle () ;
      int iZone ;
      int nMb = 0 ;
      ZONE *zp ;
      ACEIN ai = aceInCreate (p.targetFileName, p.gzi, h) ;
      Stack dna  ;
      const char *ccp ;
      
      if (!ai) messcrash ("Cannot open targetFile %s", p.targetFileName) ;
      
      p.zones = arrayHandleCreate (IBMAX, ZONE, h) ;
      for (iZone = 0 ; iZone < IBMAX ; iZone++)
	{
	  zp = arrayp (p.zones, iZone, ZONE) ;
	  zp->dna = stackHandleCreate (p.zoneSize * 1.1, h) ;
	}
      
      iZone = 0 ; zp = arrayp (p.zones, iZone, ZONE) ;
      dna = zp->dna ;
      while (aceInCard (ai))
	{
	  ccp = aceInWord (ai) ;
	  if (! ccp || ! *ccp || *ccp == '>' || *ccp == '#')
	    continue ;
	  catText (dna, ccp) ;
	  nMb += strlen (ccp) ;
	  
	  if (stackPos(dna) > p.zoneSize)
	    {
	      if (++iZone == p.ibMax)
		break ; 
	       zp = arrayp (p.zones, iZone, ZONE) ;
	       dna = zp->dna ;
	    }
	}
      
      ac_free (h1) ;
      fprintf (stderr, "TargetFile: parsed %d kb : %s\n", nMb, timeShowNow()) ;
    }
  if (p.telomere)
    {
      int n = 0 ;
      p.teloExportDoneChan = channelCreate (8, int, h) ;
      dblCountTelomeres (&p) ;
      channelGet (p.teloExportDoneChan, &n, int) ;
      goto done ;
    }

#ifdef JUNK
  /* create the zones */
  {
    int i ;
    ZONE *zp ;
    AC_HANDLE seqH ;

    seqH = p.seqH = acHandleHandleCreate (h) ;
    p.zones = arrayHandleCreate (p.maxZone, ZONE, h) ;
    for (i = 0 ; i < p.maxZone ; i++)
      {
	zp = arrayp (p.zones, i, ZONE) ;
	zp->dict = dictHandleCreate (10000, seqH) ;
	zp->dna = stackHandleCreate (p.zoneSize * 1.1, seqH) ;
	zp->bb = bitSethandleCreate (1 << (2 * WORDLENGTH - 1), h) ;
      }
  }

  /* create the genome communication channel */

 if (p.max_threads)
   {
     p.genomeChan = channelCreate (10, SHORT_READ, h) ;
     for (i = 0 ; i < NN ; i++)
       p.kmerChannel[k] =  channelCreate (10, BB2, h) ;
   }

  S2K s2kin[4], s2kOut[4] ;
  memset (s2kIn, 0, sizeof(s2kIn)) ;
  memset (s2kOut, 0, sizeof(s2kOut)) ;

  for (ii = 0 ; ii < 4 ; ii++)
    {
      s2k[ii].n = ii ;
      s2k[ii].pp = &p ;
      if (p.max_threads)
	wego_run (p.group, seq2kmer, s2kIn[ii], sizeof(S2K), s2kOut[ii], sizeof(S2K)) ;
    } 

 /* actual work, the top call to parseSeq will start filling the channels, and finish with an empty sequence */
  parseSequence () ;
  synchronize (p.group) ;

#endif

 done:
  
  ac_free (p.h) ;
  if (!p.silent)
    { 
      int mx ;
      messAllocMaxStatus (&mx) ; 
      fprintf (stderr, "// %s done: , max memory %d Mb\n", timeShowNow(),  mx) ;
     }

  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr and to ensure all pipes are closed*/
  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

