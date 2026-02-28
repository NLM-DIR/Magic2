/*  File: sam2gold.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2005
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

 * 2018_04_01
 * strategic idea

 * Count the alignment statistics in a sam file
 
 * Compare a set of sam file contructed by different methods to a gold standard
 * The gold standard describes how some artificial reads were extracted from the target
 * Each method is hoping to rediscover this hidden truth
 * The code was developped on the Baruzzo artifical human and malaria RNA-seq 

 */

#define VERSION "1.1"

/* #define ARRAY_CHECK  */
#include "ac.h"

#ifdef RESTRICT 
#define Restrict restrict
#else
#define Restrict 
#endif

/* snp chaining distance */

typedef struct s2g_struct {
  AC_HANDLE h ;
  const char *outFileName ;
  const char *goldFileName ;
  const char *snpFileName ;
  const char *inFileList ;
  const char *inFileName ;
  BOOL gzo, debug ;
  BOOL addReadPairSuffix ;
  BOOL addReadPairSuffix2 ;
  BOOL addReadPairSuffixForce ;
  BOOL exportIntronSupport ;
  BOOL anne ;
  BOOL exportAliLn ;
  BOOL samStats ;
  BOOL unique ;
  ACEIN ai ; 
  ACEOUT ao ;
  DICT *methodDict ;
  DICT *targetDict ;
  DICT *seqDict ;
  Array snpDicts ;
  DICT *scoreDict ;
  DICT *intronDict ;
  DICT *snpTypeDict ;
  KEYSET mIntron, aliLn ;
  Array i2m ;
  BigArray intronSupport ;
  BigArray hits ;
  BOOL merge, view ;
  int goldMethod ;
  const char *run ;
  const char *method ;
  long int nRawReads, nRawBases ;
  long int nMultiAli[12] ;
  int nIns[32] ;
  int nInsOk[32] ;
  int nGoldIns[32] ;
  int nDel[32] ;
  int nDelOk[32] ;
  int nGoldDel[32] ;
  int nSub[16], nSubOk[16], nGoldSub[16] ;
} S2G ;

typedef struct hit_struct {
  int method, flag, seq, target, strand, type, a1, a2, x1, x2, score, ali, gold, nerr  ;
  BOOL isComplete ;
 } HIT ;

typedef struct i2m_struct {
  int intron, type, flag, support[128] ;
 } I2M ;

#define INTRONTPMAX 101
/*************************************************************************************/

static int hitOrder (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va ;
  const HIT *vp = (const HIT *)vb ;
  int n ;

  n = up->seq - vp->seq ; if (n) return n ;
  n = up->method - vp->method ; if (n) return n ;
  n = up->score - vp->score ; if (n) return n ;
  n = up->target - vp->target ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  n = up->strand - vp->strand ; if (n) return n ;

  return 0 ;
} /* hitOrder */

/*************************************************************************************/

static int intronOrder (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va ;
  const HIT *vp = (const HIT *)vb ;
  int n ;

  n = up->type - vp->type ; if (n) return n ;
  n = up->target - vp->target ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  n = up->strand - vp->strand ; if (n) return n ;
  n = up->method - vp->method ; if (n) return n ;
  n = up->score - vp->score ; if (n) return n ;
  n = up->seq - vp->seq ; if (n) return n ;

  return 0 ;
} /* intronOrder */

/*************************************************************************************/
/*************************************************************************************/

static void s2gInit (S2G *s2g)
{
  s2g->targetDict = dictHandleCreate (256, s2g->h) ;
  s2g->seqDict = dictHandleCreate (10000, s2g->h) ;
  s2g->scoreDict = dictHandleCreate (24, s2g->h) ;
  s2g->methodDict = dictHandleCreate (64, s2g->h) ;
  s2g->intronDict = dictHandleCreate (20000, s2g->h) ;
  s2g->snpDicts = arrayHandleCreate (20000, DICT *, s2g->h) ;
  s2g->snpTypeDict = dictHandleCreate (24, s2g->h) ;
  s2g->hits = bigArrayHandleCreate (100000, HIT, s2g->h) ;
  s2g->i2m = arrayHandleCreate (100000, I2M, s2g->h) ;
  s2g->intronSupport = bigArrayHandleCreate (100000, HIT, s2g->h) ;
  s2g->mIntron = keySetHandleCreate (s2g->h) ;
  s2g->aliLn = keySetHandleCreate (s2g->h) ;

 /* register the goldMethod , it MUST be the first method */
  if (s2g->goldFileName)
    {
      char *buf = strnew (s2g->goldFileName, s2g->h) ;
      char *fNam = strchr (buf, ':') ;
      if (fNam) 
	*fNam++ = 0 ;
      dictAdd (s2g->methodDict, buf, &(s2g->goldMethod)) ;
    }

   return  ;
} /* s2gInit */

/*************************************************************************************/
/*************************************************************************************/

static BOOL s2gRegisterIntron (S2G *s2g, Array cigarettes, int flag, int method, int goldMethod, int target, int seq, int strand, const char *dna)
{
  SAMCIGAR *cgr = arrp (cigarettes, 0, SAMCIGAR) ;
  int x, i, iMax = arrayMax (cigarettes) ;
  int intron0 = 0 ;
  BigArray intronSupport = s2g->intronSupport ;
  long int intron = bigArrayMax (intronSupport) ;
  char buf[256] ;
  int  ali = 0 ;
  int lnDna = dna ? strlen(dna) : 0 ;
  BOOL isComplete = TRUE ;

  if (method != goldMethod  && 
      ! (flag & 0x100)   /* secondary mappings have flag 0x100  except first one */
      )
    {
      for (i = 0 ; i < iMax ; i++, cgr++)
	{
	  if (cgr->type == 'M' || cgr->type == '=' || cgr->type == 'X')
	    ali += cgr->dx ;
	}
      keySet (s2g->aliLn, ali) ++ ;
    }
  
  cgr = arrp (cigarettes, 0, SAMCIGAR) ;
  for (i = 0 ; i < iMax ; i++, cgr++)
    if (cgr->dx > 0)
      {
	isComplete = FALSE ;
	switch (cgr->type)
	  {
	  case 'n':
	  case 'N':
	    {
	      I2M *i2m ;
	      int type = 'N' ;
	      int a1 = (strand == 1 ? cgr->a1 : cgr->a2) ; 
	      int a2 = (strand == 1 ? cgr->a2 : cgr->a1) ; 
	      sprintf (buf, "%s:%c:%d-%d", dictName(s2g->targetDict, target), type, a1, a2) ;
	      x = 0 ;
	      dictAdd (s2g->intronDict, buf, &x) ;
	      i2m = arrayp (s2g->i2m, x, I2M) ;
	      i2m->intron = x ;	
	      i2m->type = cgr->type ;
	      if (method > 128)
		messcrash ("I2M can only store 128 methods") ;
	      i2m->support[method]++ ;
	      if (i2m->support[goldMethod] == 0)
		x = 0 ;
	      if (s2g->exportIntronSupport)
		{
		  HIT *ip = bigArrayp (intronSupport, intron++, HIT) ;
		  if (! intron0) intron0 = intron ;
		  ip->seq = seq ;
		  ip->method = method ;
		  if (method == 0)
		    messcrash ("method= 0 in s2gRegisterIntron") ;
		  ip->target = target ;
		  ip->type = type ;
		  ip->flag = flag ;
		  ip->a1 = a1 ;
		  ip->a2 = a2 ;
		  ip->score = 1 ;
		  ip->strand = strand ;
		    /* not defined in BAM format */	
		  ip->gold = x ;
		}
	      break ;
	    }
	  case 'X':
	    {
	      DICT *dict = 0;
	      int type = 0, type2 = 0 ;
	      char buf[32] ;
	      int x1 = strand == 1 ? cgr->x1 : lnDna -  cgr->x1 + 1 ;

	      sprintf (buf, "n>n") ;
	      dictAdd (s2g->snpTypeDict, buf, &type) ;
	      sprintf (buf, "n>n") ;
	      dictAdd (s2g->snpTypeDict, buf, &type2) ;
	      if (type2 > 30) type2 = 30 ;

	      s2g->nSub[0]++ ;

	      sprintf (buf, "%d_%d", x1, type) ;
	      if ((dict = array (s2g->snpDicts, seq, DICT *)) &&
		  dictFind (dict, buf, 0)
		  )
		{ s2g->nSubOk[0]++ ; if (0) s2g->nSubOk[type2]++ ; }
	      break ;
	    }
	  case 'I':
	    {
	      DICT *dict = 0;
	      int type = 0 ;
	      int ln = cgr->x2 - cgr->x1 + 1 ;
	      char buf[32] ;
	      int x1 = strand == 1 ? cgr->x1 : lnDna -  (cgr->x1 + ln - 1)  + 1 ;

	      if (ln<0 || ln > 30) ln = 30 ;
	      sprintf (buf, "Ins%d", ln) ;
	      dictAdd (s2g->snpTypeDict, buf, &type) ;
	      s2g->nIns[0]++ ;
	      s2g->nIns[ln]++ ;

	      sprintf (buf, "%d_%d", x1, type) ;
	      if ((dict = array (s2g->snpDicts, seq, DICT *)) &&
		  dictFind (dict, buf, 0)
		  )
		{ s2g->nInsOk[0]++ ; s2g->nInsOk[ln]++ ; }
	      break ;
	    }
	  case 'D':
	    {
	      DICT *dict = 0;
	      int type = 0 ;
	      int ln = cgr->a2 - cgr->a1 + 1 ;
	      char buf[32] ;
	      int x1 = strand == 1 ? cgr->x1 : lnDna - cgr->x1 + 2 ;

	      if (ln<0 || ln > 8) 
		break ;
	      sprintf (buf, "Del%d", ln) ;
	      dictAdd (s2g->snpTypeDict, buf, &type) ;
	      s2g->nDel[0]++ ;
	      s2g->nDel[ln]++ ;

	      sprintf (buf, "%d_%d", x1, type) ;
	      if ((dict = array (s2g->snpDicts, seq, DICT *)) &&
		  dictFind (dict, buf, 0)
		  )
		{ s2g->nDelOk[0]++ ; s2g->nDelOk[ln]++ ; }
	      break ;
	    }
	  default:
	    break ;
	  }
      }
  return isComplete ;
} /* s2gRegisterIntron */

/*************************************************************************************/

static void s2gParseGoldFile (S2G *s2g, const char *fNam, int method)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (fNam, 0, h) ;
  char *ccp, *cigar ;
  Array cigarettes = arrayHandleCreate (128, SAMCIGAR, h) ; 
  BigArray hits = s2g->hits ;
  long int nn , nn0 = bigArrayMax (hits) ;
  int seq, target, x1, x2, a1, a2, score, ali, ins, del, bigDel, strand ;
  HIT *hit ;

  nn = nn0 ;
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai)) 
    { /* parse a sam file
       * Expected format is tab delimited
       *  seq.1a	chr1	95672060	95672159	100M	95672060-95672159	+	TCCCTTTCACGCCTCTTCTGATTCATCTTTTTACAATCTACTCCATGATCTCATTCCTTGATGACCAGTGATTTCATCAAGCCCAGTTGGCCTTTACCTT
       */
      BOOL isComplete = FALSE ;

      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#' || *ccp == '/')
	continue ;

      dictAdd (s2g->seqDict, ccp, &seq) ;
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ;
      if (! ccp)
	continue ;
      dictAdd (s2g->targetDict, ccp, &target) ;
      a1 = a2 = 0 ;
      aceInStep (ai, '\t') ; aceInInt (ai, &a1) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &a2) ;
      aceInStep (ai, '\t') ; cigar = aceInWord (ai) ;
      if (! cigar)
	continue ;
      ali = ins = del = bigDel = 0 ;
      samParseCigar (cigar, cigarettes, a1, &a2, &x1, &x2, &ali, &ins, &del, &bigDel) ; 
      score = 1 ;
      if (! strcmp (cigar, "100M") || ! strcmp (cigar, "100="))
	score = 100 ;
      aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", 0) ; /* coords again ignore */
      if (! ccp)
	continue ;
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ;
      if (! ccp)
	continue ;
      if (! strcmp (ccp, "+"))
	strand = 1 ;
      else if (! strcmp (ccp, "-"))
	strand = -1 ;
      else
	continue ;
      if (arrayMax (cigarettes))
	isComplete = s2gRegisterIntron (s2g, cigarettes, 0, method, method, target, seq, strand, 0) ;
      
      hit = bigArrayp (hits, nn++, HIT) ;
      hit->method = method ;
      hit->flag = 0 ;
      hit->seq = seq ; 
      hit->target = target ; 
      hit->a1 = a1 ; 
      hit->a2 = a2 ; 
      hit->strand = strand ;
      hit->score = score ; 
      hit->nerr = 0 ;
      hit->isComplete = isComplete ;
    }
  fprintf (stderr, "//%s : Parsed %ld hits in the gold file %s\n"
	   , timeShowNow()
	   , nn - nn0
	   , fNam
	   ) ;
  ac_free (h) ;
 
  return  ;
} /* s2gParseGoldFile */

/*************************************************************************************/

static void s2gParseOneSamFile (S2G *s2g, const char *fNam, int method, int goldMethod)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (fNam, 0, h) ;

  DICT *unalignedDict = dictHandleCreate (100000, h) ;
  char *ccp, *dna ;
  char *cigar, seqBuf[128] ;
  BigArray hits = s2g->hits ;
  long int nn, nn0 = bigArrayMax (hits) ;
  int dnaLn ;
  int seq, seq2, flag, target, a1, a2, x1, x2, score, nerr, ali, ins, del, bigDel, strand ;
  HIT *hit ;
  Array cigarettes = arrayHandleCreate (128, SAMCIGAR, h) ; 
  KEYSET ksu = keySetHandleCreate (h) ;
  long int nAlignedBases = 0 ;
  long int nErrors = 0, nMID = 0 ;
  long int nBases = 0 ;
  long int nUnalignedReads = 0 ;
  long int nAlignedReads = 0 ;
  long int nPerfectReads = 0 ;
  BOOL isNewPerfectCandidate = FALSE ;
  ACEOUT aoPerfect = aceOutCreate (s2g->outFileName, ".perfect.list", 0, h) ;  
  BOOL nmDoesNotCountInDels = FALSE ;
  nn = nn0 ;

  if (s2g->method && strstr (s2g->method, "STAR")) nmDoesNotCountInDels =  TRUE ;
  if (s2g->method && strstr (s2g->method, "11_MagicBLAST_2018")) nmDoesNotCountInDels =  TRUE ;
  
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai)) 
    { /* parse a sam file
       * Expected format is tab delimited
       *  seq.1a	83	chr2	145633	255	98M2S	=	145500	-231	TCTTGTTAACAAATC
       */
      BOOL isComplete = FALSE ;
      int isMulti ;

      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#' || *ccp == '/' || *ccp == '@')
	continue ;

      /* seqBuf : adjust the suffix indicating paired ends */
      strncpy (seqBuf, ccp, 127) ;

      /* SAM flag
       * 0x4  0x8 unaligned 
       * 0x10 0x20  target strand
       * 0x80   second read of a pair
       * 0x100  secondary mappings
       * 1   	000000000001 	0x1   template having multiple templates in sequencing (read is paired)
       * 2 	000000000010 	0x2   each segment properly aligned according to the aligner (read mapped in proper pair)
       * 4 	000000000100 	0x4   this segment is unmapped (read1 unmapped)
       * 8 	000000001000 	0x8   next segment in the template unmapped (read2 unmapped)
       * 16 	000000010000 	0x10  this read is minus strand: SEQ being reverse complemented (read1 reverse complemented)
       * 32 	000000100000 	0x20  minus strand read2 : SEQ of the next segment in the template being reverse complemented (read2 reverse complemented)
       * 64 	000001000000 	0x40  the first segment in the template (is read1)
       * 128 	000010000000 	0x80  the last segment in the template (is read2)
       * 256 	000100000000 	0x100 not primary alignment
       * 512 	001000000000 	0x200 alignment fails quality checks
       * 1024 	010000000000 	0x400 PCR or optical duplicate
       * 2048 	100000000000 	0x800 supplementary alignment (e.g. aligner specific, could be a portion of a split read or a tied region)
       */
      aceInStep (ai, '\t') ; aceInInt (ai, &flag) ;

      if (flag & 512) /* bad quality */
	continue ;
      if (flag & 1) /* read pairs */
	{
	  char *cr = seqBuf + strlen(seqBuf) ;
	  if (flag & 64)
	    *cr++ = '>' ;
	  else /* 128 */
	    *cr++ = '<' ;
	  *cr = 0 ;
	}

      seq = seq2 = 0 ;
      if (flag & 4) /* unaligned */
	{
	  if (dictAdd (unalignedDict, seqBuf, &seq))
	    {
	      nUnalignedReads++ ;
	      s2g->nMultiAli[0]++ ;
	    }
	  else
	    messcrash ("Read %s unaligned appears twice", seqBuf) ;
	  if (dictFind (s2g->seqDict, seqBuf, &seq2))
	    messcrash ("Read %s is both aligned and unaligned", seqBuf) ;
	  continue ;
	}
      if (dictAdd (s2g->seqDict, seqBuf, &seq))
	{ /* new read */
	  nAlignedReads++ ;
	  isNewPerfectCandidate = TRUE ;
	  if (0) printf ("XXXX\t%s\n", seqBuf) ;
	}
      if (dictFind (unalignedDict, seqBuf, &seq2))
	messcrash ("Read %s is both aligned and unaligned", seqBuf) ;

      /* count the number of alignments of each read */
      isMulti = keySet (ksu, seq) + 1 ;
      keySet (ksu, seq) = isMulti  ;
      s2g->nMultiAli[(isMulti > 10 ? 10 : isMulti)]++ ;

      if (flag & 256)
	{ /* secondary mapping */
	  if (0 && isMulti == 1)
	    messcrash ("Read %s first mapping is called secondary (flag 0x100)", seqBuf) ;
	}
      else
	{
	  if (0 && isMulti > 1)
	    messcrash ("Read %s mapping %d is called primary (NOT flag 0x100)", seqBuf, isMulti) ;
	}

      strand = flag & 16 ? -1 : 1 ;
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '*') /* unaligned */
	messcrash ("why was this unaligned not flagged ") ;
      dictAdd (s2g->targetDict, ccp, &target) ;
      a1 = a2 = 0 ;
      aceInStep (ai, '\t') ; aceInInt (ai, &a1) ; a2 = a1 ;
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ;  /* some score, ignore */
      if (! ccp)
	continue ;
      aceInStep (ai, '\t') ; cigar = aceInWord (ai) ;  /* cigar */
      if (! cigar)
	continue ;
      score = 1 ;
      /*
	if (! strcmp (cigar, "100M") || ! strcmp (cigar, "100="))
	score = 100 ;
      */
      
      arrayMax (cigarettes) = 0 ;
      ali = ins = del = bigDel = 0 ;
      samParseCigar (cigar, cigarettes, a1, &a2, &x1, &x2, &ali, &ins, &del, &bigDel) ; 
      aceInStep (ai, '\t') ; aceInWord (ai) ;  /* mate target */
      aceInStep (ai, '\t') ; aceInWord (ai) ;  /* mate a1 coordinate */
      aceInStep (ai, '\t') ; aceInWord (ai) ;  /* distance to mate */
      aceInStep (ai, '\t') ; dna = aceInWord (ai) ;
      dnaLn = dna ? strlen(dna) : 0 ;
      if (isMulti == 1) 
	nBases += strlen (dna) ;
      if (arrayMax (cigarettes))
	isComplete = s2gRegisterIntron (s2g, cigarettes, flag, method, goldMethod, target, seq, strand, dna) ;
      nerr = -1 ;
      while (1)
	{
	  int n_ ;
	  char c_ ;
	  aceInStep (ai, '\t') ; ccp = aceInWord (ai) ;  /* scan the optional fields */
	  if (! ccp)
	    break ;
	  if (sscanf (ccp, "NM:i:%d%c", &n_, &c_) == 1)
	    nerr = n_ ;	      
	  else if (sscanf (ccp, "nM:i:%d%c", &n_, &c_) == 1)
	    nerr = n_ ;	      
	}

      if (nmDoesNotCountInDels)
	nerr += del + ins ;
      else
	nerr -= bigDel ;
      if (nerr < -1)
	messcrash ("nerr < 0 file %s line %d\n", aceInFileName (ai), aceInStreamLine (ai)) ;
      /* if (nerr >= del) nerr -= del ; */
	if (! (flag & (256 + 2048)))  /* primary alignment */
	  {
	    nAlignedBases += x2 - x1 + 1 - ins ;
	    nErrors += (nerr > 0 ? nerr : 0) ;
	  }

      if (isNewPerfectCandidate)
	{
	  if (ali == dnaLn && nerr == 0)
	    {
	      nPerfectReads++ ;
	      aceOutf (aoPerfect, "%s\n", dictName (s2g->seqDict, seq)) ;
	      isNewPerfectCandidate = FALSE ;
	    }
	}
      hit = bigArrayp (hits, nn++, HIT) ; nn-- ;
      hit->method = method ;
      hit->flag = flag ;
      hit->seq = seq ; 
      hit->target = target ; 
      hit->a1 = a1 ; 
      hit->a2 = a2 ; 
      hit->x1 = x1 ;
      hit->x2 = x2 ;
      hit->strand = strand ;
      hit->score = score ;
      hit->ali = ali ;
      hit->nerr = nerr ;
      hit->isComplete = isComplete ;
   }

  fprintf (stderr, "//%s : Parsed %ld lines in the sam file %s\n"
	   , timeShowNow()
	   , nn - nn0
	   , fNam
	   ) ;
  /* because read appearing 7 times previously occured 6 5 4 3 2 1 , we differentiate */
  for (int i = 1 ; i < 9 ; i++)
    s2g->nMultiAli[i] -= s2g->nMultiAli[i+1] ;

  if (1)
    {
      int nReads = dictMax (s2g->seqDict) ;
      ACEOUT ao = aceOutCreate (s2g->outFileName, ".samStats", 0, h) ;

      nUnalignedReads = s2g->nRawReads - nAlignedReads ;
      
      if (s2g->nRawReads)
	aceOutf (ao, "%s\t%s\tnRawReads\t%ld\n"
		 , s2g->run, s2g->method
		 , s2g->nRawReads
		 ) ;

      if (s2g->nRawBases)
	aceOutf (ao, "%s\t%s\tnRawBases\t%ld\n"
		 , s2g->run, s2g->method
		 , s2g->nRawBases
		 ) ;
      aceOutf (ao, "%s\t%s\tnReads\t%d\n"
	       , s2g->run, s2g->method
	       , nReads
	       ) ;
      aceOutf (ao, "%s\t%s\tnAlignedBases\t%ld\t%.2f%%\n"
	       , s2g->run, s2g->method
	       , nAlignedBases
	       , (100.0 * nAlignedBases) / (s2g->nRawBases + .0000001)
	       ) ;

      aceOutf (ao, "%s\t%s\tnMismaches_and_InDels\t%ld\t%.4f%%\n"
	       , s2g->run, s2g->method
	       , nErrors
	       , (100.0 * nErrors) / (nAlignedBases + .0000001)
	       ) ;

      aceOutf (ao, "%s\t%s\tnAlignedReads\t%d\t%.2f%%\n"
	       , s2g->run, s2g->method
	       , nAlignedReads
	       , (100.0 * nAlignedReads) / (s2g->nRawReads + .0000001)
	       ) ;

      aceOutf (ao, "%s\t%s\tnPerfectReads\t%d\t%.2f%%\n"
	       , s2g->run, s2g->method
	       , nPerfectReads
	       , (100.0 * nPerfectReads) / (s2g->nRawReads + .0000001)
	       ) ;

      for (int i = 0 ; i < 1 ; i++)
	aceOutf (ao, "%s\t%s\tnMultiAligned %ld times\t%ld\t%.2f%%\n"
		 , s2g->run, s2g->method
		 , i
		 ,  nUnalignedReads
		 , (100.0 * nUnalignedReads) / (s2g->nRawReads + .0000001)
		 ) ;
      for (int i = 1 ; i <= 10 ; i++)
	aceOutf (ao, "%s\t%s\tnMultiAligned %ld times\t%d\t%.2f%%\n"
		 , s2g->run, s2g->method
		 , i
		 ,  s2g->nMultiAli[i]
		 , (100.0 * s2g->nMultiAli[i]) / (s2g->nRawReads + .0000001)
		 ) ;
    }
      
  ac_free (h) ;
 
  return  ;
} /* s2gParseOneSamFile */

/*************************************************************************************/
/* create one dist per sequence
 * so we do not repeat the sequence names many times 
 */
static void s2gParseSnpFile (S2G *s2g, const char *fNam, int method)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (fNam, 0, h) ;
  char *ccp, buf[128] ;
  int seq, pos, type ;
  DICT *dict ;
  
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai)) 
    {
      ccp = aceInWord (ai) ; /* sequence name */
      if (!ccp || ! *ccp)
	continue ;
      dictAdd (s2g->seqDict, ccp, &seq) ;
      dict = array (s2g->snpDicts, seq, DICT *) ;
      if (! dict)
	dict = array (s2g->snpDicts, seq, DICT *) = dictHandleCreate (128, s2g->h) ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &pos))
	continue ;
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* type */
      if (!ccp || ! *ccp)
	continue ;
      if (!strncasecmp (ccp, "Ins", 3))
	{
	  int ln = strlen (ccp+3) ;
	  if (ln > 30) ln = 30 ;
	  sprintf (buf, "Ins%d",ln) ;
	  dictAdd (s2g->snpTypeDict, buf, &type) ;
	  (s2g->nGoldIns[0])++ ;
	  (s2g->nGoldIns[ln])++ ;
	}
      else if (!strncasecmp (ccp, "Del", 3))
	{
	  int ln = strlen (ccp+3) ;
	  if (ln > 30) ln = 30 ;
	  sprintf (buf, "Del%d", ln) ;
	  dictAdd (s2g->snpTypeDict, buf, &type) ;
	  (s2g->nGoldDel[0])++ ;
	  (s2g->nGoldDel[ln])++ ;
	}
      else if (ccp[1] == '>')
	{
	  sprintf (buf, "n>n") ;
	  dictAdd (s2g->snpTypeDict, buf, &type) ;
	  (s2g->nGoldSub[0])++ ;
	}
      sprintf (buf, "%d_%d", pos,type) ;
      dictAdd (dict, buf, 0) ;
    }
  return  ;
} /* s2gParseSnpFile */

/*************************************************************************************/

static void s2gParse (S2G *s2g)
{
  AC_HANDLE h = ac_new_handle () ;
  char *cp, *cq, *list ;
  int method, methodFile ;
  ACEIN ai = 0 ;
  KEYSET ks = keySetHandleCreate (h) ;
  DICT *dict = dictHandleCreate (64, h) ;

 
  /* check that all files exist */
  if (s2g->inFileList)
    {
      list = strnew (s2g->inFileList, h) ;
      cp = cq = list ;
      while (cp)
	{
	  if (s2g->samStats)
	    { /* just check the existence of the sam files */
	      if (! (ai = aceInCreate (cp, 0, h)))
		messcrash ("Bad parameter -i, cannot open file %s\n", cp) ;
	      ac_free (ai) ;

	      dictAdd (s2g->methodDict, s2g->method, &method) ;
	      dictAdd (dict, cp, &methodFile) ;
	      keySet (ks, method) = methodFile ;
	      
	      break ;
	    }
	  cq = strchr (cp, ':') ;
	  if (! cq)
	    messcrash ("Bad parameter -i, expected m1:f1,m2:f2...\nMissing : in the inFileList %s\n at %s"
		       ,  s2g->inFileList, cp
		       ) ;
	  *cq = 0 ; method = 0 ;
	  dictAdd (s2g->methodDict, cp, &method) ;
	  cp = cq = cq + 1 ;
	  cq = strchr (cp, ',') ;
	  if (cq) *cq++ = 0 ;
	  dictAdd (dict, cp, &methodFile) ;
	  keySet (ks, method) = methodFile ;
	  if (! (ai = aceInCreate (cp, 0, h)))
	    messcrash ("Bad parameter -i, cannot open file %s\n", cp) ;
	  ac_free (ai) ;
	  cp = cq ;
	}
    }

  if (s2g->goldFileName)
    {
      /* 
	 char *buf = strnew (s2g->goldFileName, h) ;
      const char *fNam = strchr (buf, ':') ;
      if (fNam) 
	fNam++ ;
      else
	fNam = s2g->goldFileName ;
      */
      s2gParseGoldFile (s2g, s2g->goldFileName, s2g->goldMethod) ;
    }

   if (s2g->snpFileName)
    {
      /*  automatically analyzed by aceInCreate (aceInLabel) 
	char *buf = strnew (s2g->snpFileName, h) ;
	const char *fNam = strchr (buf, ':') ;
	if (fNam) 
	fNam++ ;
	else
	fNam = s2g->snpFileName ;
      */
      s2gParseSnpFile (s2g, s2g->snpFileName, s2g->goldMethod) ;
    }

   if (s2g->inFileList)
     {
       for (method = 0 ; method < keySetMax (ks) ; method++)
	 if (keySet (ks, method))
	   s2gParseOneSamFile (s2g, dictName (dict, keySet(ks, method)), method, s2g->goldMethod) ;
    }

  ac_free (h) ;

  return  ;
} /* s2gParse */

/*************************************************************************************/
/*************************************************************************************/

static void s2gAnalyze (S2G *s2g)
{ 
  AC_HANDLE h = ac_new_handle () ;
  BigArray hits = s2g->hits ;
  long int ii, iiMax = bigArrayMax (hits) ;
  int m ;
  int mMax = dictMax(s2g->methodDict) + 1 ;
  int kkk [mMax] ;
  HIT *up ;

  KEYSET nExactRead = keySetHandleCreate (s2g->h) ;
  KEYSET nPartialRead = keySetHandleCreate (s2g->h) ;
  KEYSET nOverlappingRead = keySetHandleCreate (s2g->h) ;
  KEYSET nWildRead = keySetHandleCreate (s2g->h) ;
  KEYSET nMappedOnce = keySetHandleCreate (s2g->h) ;
  KEYSET nMultimapped = keySetHandleCreate (s2g->h) ;
  KEYSET nUnaligned = keySetHandleCreate (s2g->h) ;

  KEYSET nExactAli = keySetHandleCreate (s2g->h) ;
  KEYSET nPartialAli = keySetHandleCreate (s2g->h) ;
  KEYSET nOverlappingAli = keySetHandleCreate (s2g->h) ;
  KEYSET nWildAli = keySetHandleCreate (s2g->h) ;

  KEYSET nFully_alignedExact  = keySetHandleCreate (s2g->h) ;
  KEYSET nFully_alignedWithError = keySetHandleCreate (s2g->h) ;
  KEYSET nPartialExact = keySetHandleCreate (s2g->h) ;
  KEYSET nPartialWithError = keySetHandleCreate (s2g->h) ;

  
  ACEOUT ao ;
  ACEOUT aotsv ;

  bigArraySort (hits, hitOrder) ;
  
  for (up = bigArrp (hits, 0, HIT), ii = 0 ; ii < iiMax ; ii++, up++)
    {
      long int jj ;
      int seq, a1, a2, strand, target ;
      HIT *vp ;

      if (up->method == 1)
	{
	  keySet (nExactRead, 1)++ ;
	  keySet (nExactAli, 1)++ ;

	  seq = up->seq ;
	  target = up->target ;
	  a1 = up->a1 ;
	  a2 = up->a2 ;
	  strand = up->strand ;
	}
      else
	continue ;

      /* non-hiearchic count , i.e. per ali */
      memset (kkk, 0, sizeof(kkk)) ;
      kkk[1]++ ;
      for (jj = ii + 1, vp = up + 1 ; jj < iiMax && vp->seq == seq ; jj++, vp++)
	{ 
	  int m = vp->method ;
	  if (vp->target == target && vp->a1 == a1 && vp->a2 == a2 && vp->strand == strand)
	    { keySet (nExactAli, m)++ ; }
	  else if (vp->target == target && vp->a1 >= a1-10 && vp->a2 <= a2+10 && vp->strand == strand)
	    { keySet (nPartialAli, m)++ ; }
	  else if (vp->target == target && vp->a1 < a2 && vp->a2 > a1 && vp->strand == strand)
	    { keySet (nOverlappingAli, m)++ ; }
	  else
	    keySet (nWildAli, m)++ ;
	  kkk[m]++ ;
	}

      for (m = 1 ; m < mMax ; m++)
	switch (kkk[m])
	  {
	  case 0:
	    keySet (nUnaligned, m)++ ;
	    break ;
	  case 1:
	    keySet (nMappedOnce, m)++ ;
	    break ;
	  default:
	    keySet (nMultimapped, m)++ ;
	    break ;
	  }

      /* hierarchic count , i.e. per read */
      for (m = 1 ; m < mMax ; m++)
	{
	  int ok1 = 0 ;
	  int ok2 = 0 ;
	  int ok3 = 0 ;
	  int ok4 = 0 ;
	  int ok11 = 0 ;
	  int ok12 = 0 ;
	  int ok13 = 0 ;
	  int ok14 = 0 ;

	  /* search for best case */
	  for (jj = ii + 1, vp = up + 1 ; jj < iiMax && vp->seq == seq ; jj++, vp++)
	    {
	      if (vp->method != m)
		continue ;
	      if (vp->target == target && vp->a1 == a1 && vp->a2 == a2 && vp->strand == strand)
		ok1++ ;
	      else if (vp->target == target && vp->a1 >=a1-10  && vp->a2 <= a2+10 && vp->strand == strand)
		ok2++ ;
	      else if (vp->target == target && vp->a1 < a2 && vp->a2 > a1 && vp->strand == strand)
		ok3++ ;
	      else
		ok4++ ;

	      if (up->isComplete)
		{
		  if (! up->nerr)
		    ok11++ ;
		  else
		    ok12++ ;
		}
	      else
		{
		  if (! up->nerr)
		    ok13++ ;
		  else
		    ok14++ ;
		}
	    }
	  if (ok1)
	    keySet (nExactRead, m)++ ;
	  else if (ok2)
	    keySet (nPartialRead, m)++ ;
	  else if (ok3)
	    keySet (nOverlappingRead, m)++ ;
	  else if (ok4)
	    keySet (nWildRead, m)++ ;

	  if (ok11)
	    keySet (nFully_alignedExact, m)++ ;
	  else if (ok12)
	    keySet (nFully_alignedWithError, m)++ ;
	  else if (ok13)
	    keySet (nPartialExact, m)++ ;
	  else if (ok14)
	    keySet (nPartialWithError, m)++ ;
	}
      
      /* proceed to next read */
      while (ii < iiMax - 1 && up[1].seq == seq)
	{ ii++ ; up++ ; }
    }

  ao = aceOutCreate (s2g->outFileName, ".qc", 0, h) ;
  aotsv = aceOutCreate (s2g->outFileName, ".introns.tsv", 0, h) ;

  aceOutDate (ao, "###", "sam to gold") ;
  aceOutf (ao, "# Method\tReads in run\tExact alignments\tPartial alignments\tOverlapping alignments\tWild alignments\tExact reads\tPartial reads\tOverlapping reads\tWild reads\tUniquely mapped reads\tMulti mapped reads\tUnaligned reads\tFully_alignedExact\tFully_alignedWithError\tPartialExact\tPartialWithError\n") ;
  for (m = 1 ; m <= dictMax(s2g->methodDict) ; m++)
    {
      aceOutf (ao, "%s", dictName (s2g->methodDict, m)) ;
      aceOutf (ao, "\t%d", keySet (nExactRead, 1)) ;
      aceOutf (ao, "\t%d", keySet (nExactAli, m)) ;
      aceOutf (ao, "\t%d", keySet (nPartialAli, m)) ;
      aceOutf (ao, "\t%d", keySet (nOverlappingAli, m)) ;
      aceOutf (ao, "\t%d", keySet (nWildAli, m)) ;
      aceOutf (ao, "\t%d", keySet (nExactRead, m)) ;
      aceOutf (ao, "\t%d", keySet (nPartialRead, m)) ;
      aceOutf (ao, "\t%d", keySet (nOverlappingRead, m)) ;
      aceOutf (ao, "\t%d", keySet (nWildRead, m)) ;
      aceOutf (ao, "\t%d", keySet (nMappedOnce, m)) ;
      aceOutf (ao, "\t%d", keySet (nMultimapped, m)) ;
      aceOutf (ao, "\t%d", keySet (nUnaligned, m)) ;
      aceOutf (ao, "\t%d", keySet (nFully_alignedExact, m)) ;
      aceOutf (ao, "\t%d", keySet (nFully_alignedWithError, m)) ;
      aceOutf (ao, "\t%d", keySet (nPartialExact, m)) ;
      aceOutf (ao, "\t%d", keySet (nPartialWithError, m)) ;
      aceOutf (ao, "\n") ;

      aceOutf (aotsv, "GoldMap\t%s", dictName (s2g->methodDict, m)) ;
      aceOutf (aotsv, "\t16") ;
      aceOutf (aotsv, "\t%d", keySet (nExactRead, 1)) ;
      aceOutf (aotsv, "\t%d", keySet (nExactAli, m)) ;
      aceOutf (aotsv, "\t%d", keySet (nPartialAli, m)) ;
      aceOutf (aotsv, "\t%d", keySet (nOverlappingAli, m)) ;
      aceOutf (aotsv, "\t%d", keySet (nWildAli, m)) ;
      aceOutf (aotsv, "\t%d", keySet (nExactRead, m)) ;
      aceOutf (aotsv, "\t%d", keySet (nPartialRead, m)) ;
      aceOutf (aotsv, "\t%d", keySet (nOverlappingRead, m)) ;
      aceOutf (aotsv, "\t%d", keySet (nWildRead, m)) ;
      aceOutf (aotsv, "\t%d", keySet (nMappedOnce, m)) ;
      aceOutf (aotsv, "\t%d", keySet (nMultimapped, m)) ;
      aceOutf (aotsv, "\t%d", keySet (nUnaligned, m)) ;
      aceOutf (aotsv, "\t%d", keySet (nFully_alignedExact, m)) ;
      aceOutf (aotsv, "\t%d", keySet (nFully_alignedWithError, m)) ;
      aceOutf (aotsv, "\t%d", keySet (nPartialExact, m)) ;
      aceOutf (aotsv, "\t%d", keySet (nPartialWithError, m)) ;
      aceOutf (aotsv, "\n") ;

    } 
  
  ac_free (h) ;
  return  ;
} /* s2gAnalyze */   

/*************************************************************************************/
/*************************************************************************************/

 static void s2gExportAliLn (S2G *s2g)
 {
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (s2g->outFileName, ".aliLn", s2g->gzo, h) ;
  KEYSET aliLn = s2g->aliLn ;
  int jj ;

  for (jj = 0 ; jj < keySetMax (aliLn) ; jj++)
    aceOutf (ao, "%d\t%d\n", jj, keySet (aliLn, jj)) ;

  ac_free (h) ;
  return  ;
} /* s2gExportAliLn */   

/*************************************************************************************/

static void s2gExportIntronSupport (S2G *s2g)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = 0 ;
  BigArray intronSupport = s2g->intronSupport ;
  long int ii, jj, iMax = bigArrayMax (intronSupport) ;
  HIT *up, *vp ;
  int goldMethod = s2g->goldMethod ;
  int type ; 
  const char *types = "NID" ;
  const char *titles[] = { ".intron_support", ".insertion_support", ".deletion_support" } ;
  bigArraySort (intronSupport, intronOrder) ;

  /* export */
  for (type = 0 ; type < 3 ; type++)
    {
      ao = aceOutCreate (s2g->outFileName, titles[type], s2g->gzo, h) ;
      aceOutDate (ao, "###", hprintf (h, "%s", titles[type])) ;
      for (ii = 0, up = bigArrp (intronSupport, ii, HIT) ; ii < iMax ; ii++, up++) 
	if (up->type == types[type])
	  aceOutf (ao, "%s\t%s:%d-%d\t%s\t%s\t%s\n"
		   , dictName (s2g->seqDict, up->seq)
		   , dictName (s2g->targetDict, up->target) 
		   , up->a1, up->a2
		   , dictName (s2g->methodDict, up->method) 
		   , up->gold ? "IN" : "OUT_OF" 
		   , goldMethod ? dictName (s2g->methodDict, goldMethod) : ""
		   ) ;
      ac_free (ao) ;
    }

  /* compress the info */
  for (ii = 0, up = bigArrp (intronSupport, ii, HIT) ; ii < iMax ; ii++, up++) 
    {
      int target = up->target ;
      int strand = up->strand ;
      int a1 = up->a1 ;
      int a2 = up->a2 ;
      int method = up->method ;
      int type = up->type ;

      if (method)
	for (jj = ii, vp = up ; jj < iMax ; jj++, vp++)
	  {
	    if (ii < jj)
	      {
		if (vp->target == target &&
		    vp->type == type &&
		    vp->a1 == a1 &&
		    vp->a2 == a2 &&
		    vp->strand == strand &&
		    vp->method == method
		    )
		  { up->score += vp->score ; vp->method = 0 ;}
		else
		  break ; 
	      }
	  }
    }
  
  /* keep happy few */
  for (ii = jj = 0, up = vp =bigArrp (intronSupport, ii, HIT) ; ii < iMax ; ii++, up++) 
    {
      if (up->method)
	{
	  if (up != vp)
	    *vp = *up ;
	  vp++ ; jj++ ;
	}
	}
  iMax = bigArrayMax (intronSupport) = jj ;

  ac_free (h) ;
  return  ;
} /* s2gExportIntronSupport */

/*************************************************************************************/

static void s2gExportIntrons (S2G *s2g)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = 0 ;
  Array I2m = s2g->i2m ;
  KEYSET mIntron = s2g->mIntron ;
  int ii, iMax = arrayMax (I2m) ;
  int m, mMax = dictMax(s2g->methodDict) + 1 ;
  I2M *i2m ;
  int type ;
  const char *types = "NID" ;
  const char *titles[] = { "Introns", "Insertions", "Deletions" } ;
  
  /* export */
  for (type = 0 ; type < 3 ; type++)
    {
      ao = aceOutCreate (s2g->outFileName, hprintf (h,".%s", titles[type]), s2g->gzo, h) ;
      aceOutDate (ao, "###", hprintf (h, "%s support in various methods", titles[type])) ; 
      aceOutf (ao, "# %s", titles[type]) ;
      for (m = 1 ; m < mMax ; m++)
	aceOutf (ao, "\t%s", dictName (s2g->methodDict, m)) ;
      for (ii = 0, i2m = arrp (I2m, ii, I2M) ; ii < iMax ; ii++, i2m++) 
	{
	  if (! i2m->intron || i2m->type != types[type])
	    continue ;
	  aceOutf (ao, "\n%s", dictName (s2g->intronDict, i2m->intron)) ;
	  
	  for (m = 1 ; m < mMax ; m++)
	    {
	      int i ;
	      int *kk = i2m->support ;
	      for (i = 1 ; i < INTRONTPMAX ; i++)
		{
		  if (kk[m] >= i)
		    keySet (mIntron, 128*4 * INTRONTPMAX*type + 4 * INTRONTPMAX * m + 4*i + 0)++ ; /* intron seen */
		  if (kk[1] >= 1 && kk[m] >= i)
		    keySet (mIntron, 128*4 * INTRONTPMAX*type + 4 * INTRONTPMAX * m + 4*i + 1)++ ; /* True positive*/
		  if (kk[1] < 1 && kk[m] >= i)
		    keySet (mIntron, 128*4 * INTRONTPMAX*type + 4 * INTRONTPMAX * m + 4*i + 2)++ ; /* False positive*/
		  if (kk[1] >= 1 && kk[m] < i)
		    keySet (mIntron, 128*4 * INTRONTPMAX*type + 4 * INTRONTPMAX * m + 4*i + 3)++ ; /* False negative*/
		}
	      
	      aceOutf (ao, "\t%d", kk[m]) ;
	    }
	}
      aceOutf (ao, "\n") ;
      ac_free (ao) ;
    }

  if (1)
    {
      ACEOUT ao = aceOutCreate (s2g->outFileName, ".delins.qc", 0, h) ;
      ACEOUT aotsv = aceOutCreate (s2g->outFileName, ".delins.tsv", 0, h) ;
      for (type = 0 ; type < 3 ; type++)
	{
	  const char *titles[] = { "Intron", "Insertion", "Deletion" } ;
	  
	  aceOutf (ao, "\n\n# Method\tMinimal support\t%s seen\tTruth\tTrue positive\tFalse positive\tFalse negative\tPrecision\tRecall\tF score\n", titles[type]) ;
	  aceOutf (aotsv, "# Method\t8\tMinimal support\tTruth\t\t%s seen\tTrue positive\tFalse positive\tFalse negative\tPrecision\tRecall\tF score\n", titles[type]) ;
	  for (ii = 1 ; ii < INTRONTPMAX ; ii++)
	    for (m = 1 ; m < mMax ; m++)
	      {
		KEYSET mIntron = s2g->mIntron ;
		int s0, s, tp, fp, fn ;
		float p, r, f ;
		
		s0 = keySet (mIntron, 128*4 * INTRONTPMAX*type + 4 * INTRONTPMAX * s2g->goldMethod + 4 * 1 + 0) ;
		s = keySet (mIntron, 128*4 * INTRONTPMAX*type + 4 * INTRONTPMAX* m + 4*ii + 0) ;
		tp = keySet (mIntron, 128*4 * INTRONTPMAX*type + 4 * INTRONTPMAX* m + 4*ii + 1) ;
		fp = keySet (mIntron, 128*4 * INTRONTPMAX*type + 4 * INTRONTPMAX* m + 4*ii + 2) ;
		fn = keySet (mIntron, 128*4 * INTRONTPMAX*type + 4 * INTRONTPMAX* m + 4*ii + 3) ;
		
		p = tp ; p = tp ? p/(tp + fp) : 0 ;
		r = tp ; r = tp ? r/ (tp + fn) : 0 ;
		f = p+r > 0 ? 2.0 * p * r / (p + r) : 0 ;
		
		if (! s)
		  continue ;
		
		aceOutf (aotsv,"%s_%d\t%s", titles[type], ii, dictName (s2g->methodDict, m)) ;
		aceOutf (aotsv,"\t8\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\n", s0, s,tp,fp,fn,p,r,f) ;
		
		aceOutf (ao, "%s", dictName (s2g->methodDict, m)) ;
		aceOutf (ao, "\t%d", ii) ;
		aceOutf (ao, "\t%d", s0) ;
		aceOutf (ao, "\t%d", s) ;
		aceOutf (ao, "\t%d", tp) ;
		aceOutf (ao, "\t%d", fp) ;
		aceOutf (ao, "\t%d", fn) ;
		aceOutf (ao, "\t%.4f", p) ;
		aceOutf (ao, "\t%.4f", r) ;
		aceOutf (ao, "\t%.4f", f) ;
		aceOutf (ao, "\n") ; 
	      } 
	  aceOutf (ao, "\n") ;
	}
      
      ac_free (aotsv) ;
      ac_free (ao) ;
    }
  ac_free (h) ;
  return  ;
} /* s2gExportIntrons */

/*************************************************************************************/
/*************************************************************************************/
typedef struct tt_struct {
  int tag, n, run, z[128] ;
 
} TT ;


static int ttOrder (const void *va, const void *vb)
{
  const TT *up = (const TT *)va ;
  const TT *vp = (const TT *)vb ;
  int n ;

  n = up->tag - vp->tag ; if (n) return n ;
  n = up->run - vp->run ; if (n) return n ;

  return 0 ;
} /* ttOrder */

static int ttViewOrder (const void *va, const void *vb)
{
  const TT *up = (const TT *)va ;
  const TT *vp = (const TT *)vb ;
  int n ;

  n = up->tag - vp->tag ; if (n) return n ;
  n = up->run - vp->run ; if (n) return n ;

  return 0 ;
} /* ttOrder */

/*************************************************************************************/
/* copied from aliqc : add up the counts in a tsv file */
static int s2gParseTsv (S2G *s2g, DICT *tagDict, DICT *runDict, Array tags)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = 0 ;
  int nFile = 0, i, iMax = 0, n, tag, run, x ;
  const char *ccp ;
  TT *tt ;

  while ((ai = aceInCreateFromList (ai, &nFile, s2g->inFileList, FALSE, h)))
    {
      aceInSpecial (ai, "\n") ;
      while (aceInCard (ai))
	{
	  ccp = aceInWord (ai) ;
	  if (! ccp || *ccp == '#')
	    continue ;
	  dictAdd (tagDict, ccp, &tag) ;
	  aceInStep (ai, '\t') ;
	  ccp = aceInWord (ai) ;
	  if (! ccp || *ccp == '#')
	    continue ;
	  if (runDict)  /* run_name , ignore in merge case */
	    dictAdd (runDict, ccp, &run) ;
	  aceInStep (ai, '\t') ;
	  if (! aceInInt (ai, &n))
	    continue ;
	  tt = arrayp (tags, iMax++, TT) ;
	  tt->tag = tag ;
	  tt->run = run ;
	  
	  if (n < 0 || n > 128) 
	    messcrash ("Overflow, number of fields = %d > 128 at line %d of file %s.\nWe do not expect more than 128 value following a tag, please edit the code or reformat  the data"
		       , aceInStreamLine (ai)
		       , aceInFileName (ai)
		       , n) ;
	  for (i = 0 ; i < n && i < 128 ; i++)
	    {
	      if (! aceInStep (ai, '\t') || ! aceInInt (ai, &x))
		break ;
	      tt->n = i + 1 ;
	      tt->z[i] = x ;
	    }
	}
    }
  ac_free (h) ;

  return arrayMax (tags) ;
} /*  s2gParseTsv */

/*************************************************************************************/

static int s2gCumulTsv (Array tags)
{
  int i, ii, jj, iMax = arrayMax (tags) ;
  TT *tt, *tt2 ;

  /* sort and cumulate the values */
  arraySort (tags, ttOrder) ;

  for (ii = 0, tt = arrp (tags, ii, TT) ; ii < iMax ; ii++, tt++)
    {
      if ( !tt->tag)
	continue ;

      for (jj = ii + 1, tt2 = tt+1 ; jj < iMax && tt2->tag == tt->tag ; jj++, tt2++)
	{
	  if (tt->n < tt2->n)
	    tt->n = tt2->n ;
	for (i = 0 ; i < tt->n  ; i++)
	    tt->z[i] += tt2->z[i] ;
	  tt2->tag = 0 ;
	}
    } 
  /* clean up */
  for (ii = jj = 0, tt2 = tt = arrp (tags, ii, TT) ; ii < iMax ; ii++, tt++)
    {
      if ( !tt->tag)
	continue ;
      if (tt != tt2)
	*tt2 = *tt ;
      jj++ ; tt2++ ; 
    } 
  arrayMax (tags) = jj ;
  return arrayMax (tags) ;
} /* s2gCumulTsv */

/*************************************************************************************/
/* copied from aliqc : add up the counts in a tsv file */
static void s2gExportMergedTsv (S2G *s2g, DICT *tagDict, const char *runName, Array tags)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (s2g->outFileName, ".tsv", s2g->gzo, h) ;
  int i, ii, iMax = arrayMax (tags) ;
  TT *tt ;
  
  for (ii = 0, tt = arrp (tags, ii, TT) ; ii < iMax ; ii++, tt++)
    {
      aceOutf (ao, "%s", dictName (tagDict, tt->tag)) ;
      aceOutf (ao, "\t%s", runName) ;
      aceOutf (ao, "\t%d", tt->n) ; 
      for (i = 0 ; i < tt->n ; i++)
	aceOutf (ao, "\t%d", tt->z[i]) ;
      aceOutf (ao, "\n") ;
    }
  ac_free (h) ;
  return ;
} /* s2gExportMergedTsv */

/*************************************************************************************/
/* copied from aliqc : add up the counts in a tsv file */
/* Interpret and display a table looking like 
GoldMap	GOLD	8	20000000	20000000	0	214663	0	0	20000000	0	0
GoldMap	HG19t1r1..g_MB4_w60	8	20000000	16140345	54636	0	17877429	1865491	19912192	87808	0
# Method	8	Minimal support	Intron seen	True positive	False positive	False negative	Precision	Recall	F score
Introns_1	GOLD	7	137919	137919	0	0	1.0000	1.0000	1.0000
Introns_1	HG19t1r1..g_MB4_w60	7	135315	129552	5763	8367	0.9574	0.9393	0.9483
Introns_2	GOLD	7	131218	131218	0	6701	1.0000	0.9514	0.9751
Introns_2	HG19t1r1..g_MB4_w60	7	124160	120187	3973	17732	0.9680	0.8714	0.9172
Introns_3	GOLD	7	124935	124935	0	12984	1.0000	0.9059	0.9506
*/
static void s2gViewTsv (S2G *s2g, DICT *tagDict, DICT *runDict, Array tags)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (s2g->outFileName, ".intron_qc.txt", s2g->gzo, h) ;
  int i, ii, iMax = arrayMax (tags) ;
  TT *tt ;
  int gold, itr[21] ;

  arraySort (tags, ttViewOrder) ;

  dictAdd (runDict, "GOLD", &gold) ;
  for (i = 1 ; i < 21 ; i++)
    dictAdd (tagDict, hprintf (h, "Introns_%d", i), &(itr[i])) ;
  aceOutDate (ao, "###", "Comparison of detected introns to the Gold standard") ;
  aceOutf (ao, "Run\tMethod\tMinimal intron support\tGOLD set\tDetected\tTrue Positives\tFalsePositives\tFalse Negatives\tPrecision\tRecall\tF score\n") ;

  for (ii = 0, tt = arrp (tags, ii, TT) ; ii < iMax ; ii++, tt++)
    {
      int run = tt->run ;
      
      aceOutf (ao, "%s", dictName (tagDict, tt->tag)) ;
      aceOutf (ao, "%s", dictName (runDict, run)) ;
      aceOutf (ao, "\t%d", tt->n) ; 
      for (i = 0 ; i < tt->n ; i++)
	aceOutf (ao, "\t%d", tt->z[i]) ;
      aceOutf (ao, "\n") ;
    }
  ac_free (h) ;
  return ;
} /* s2gViewTsv */

/*************************************************************************************/
/*************************************************************************************/
/* ~/SEQC2_RNA_2018 dedicated code, compute fractions and linear genes */

static Array *s2gSeqcAcedbParse (S2G *s2g, const char *types, Array *aa)
{
   AC_HANDLE h = ac_new_handle () ;
   ACEIN ai = aceInCreate (s2g->inFileList, 0, h) ;
   int maxType  = strlen (types) ;
   int t, i, n, p ;
   DICT *dict ;
   long int cc = 0 , cumul[maxType] ;

   memset (cumul, 0, sizeof(cumul)) ;
   aceInSpecial (ai, "\n") ;
   dict = s2g->intronDict ;
   while (aceInCard (ai))
     {  /* expected format   ACEDB gene_name read_count */
       const char *cq, *cp = aceInWord (ai) ;

       if (!cp || cp[1]) continue ;
       cq = strchr (types, *cp) ;
       if (!cq) continue ;
       t = cq - types ;

       aceInStep (ai, '\t') ;
       cp = aceInWord (ai) ;
       if (!cp)
	 continue ;
       dictAdd (dict, cp, &n) ;
       
       aceInStep (ai, '\t') ;
       aceInInt (ai, &p) ;
       array (aa[t], n, float)+= p ;
       cumul[t]+= p ; cc += p ;
     }

   /* normalize */
   for(t = 0 ; t < maxType ; t++)
     {
       int iMax = dictMax (dict) + 1 ;
       Array a = aa[t] ;
       float z = array (a, iMax - 1, float) ; /* make room */
       if (cumul[t])
	 {
	   float * Restrict xp = arrp (a, 0, float) ;
	   z = cc / (5.0 * cumul[t]) ;

	   for (i = 0 ; i < iMax ; i++)
	     xp[i] *= z ;
	 }
     }

   ac_free (h) ;
   return aa ;
} /* s2gSeqcAcedbParse */

/*************************************************************************************/
/* optmize for f the equation f*x + (1-f)*y - z == 0 */
static float s2gSeqcAcedbGetF (S2G *s2g, Array *aa, int t, float f, float df)
{
  AC_HANDLE h = ac_new_handle () ; 
  float * Restrict x = arrp (aa[0], 0, float) ;
  float * Restrict y = arrp (aa[1], 0, float) ;
  float * Restrict z = arrp (aa[t], 0, float) ;
  int k ;
  int ndd, i, iMax = arrayMax (aa[0]) ;
  double bestDd = 0, ddd, dd[iMax] ;
  float bestf = 0 ;
  Array cumul = arrayHandleCreate (iMax, float, h) ; 
  float * Restrict w = arrp (cumul, 0, float) ;
  float junkOver = 0 ; /* kunk the top 5% expressed object */
  float junkBelow = 0 ; /* kunk the top 5% expressed object */
  int minX = 1000 ;

  memset (dd, 0, sizeof(dd)) ;
  fprintf (stderr, "%f\t%f", f, df) ;
  
  array (cumul, iMax -1, float) = 0 ;
  for (i = 0 ; i < iMax ; i++)
    w[i] = x[i] + y[i] + z[i] ;
  arraySort (cumul, floatOrder) ;
  i = iMax - 50 ;
  junkOver = arr (cumul, i, float) ;

  i = iMax/20 ; i = iMax - 1000  ; 
  junkBelow = arr (cumul, i, float) ;

  for (k = 0 ; k < 10 ; k++)
    {
      float u = f + (k - 5) * df ;
      float v = 1 - u ;

      ndd = 0 ;

      memset (dd, 0, sizeof(dd)) ;
      for (i = 0 ; i < iMax ; i++)
	{
	  if (
	      (x[i] > minX && y[i] > minX && z[i] > minX) &&
	      (x[i] + y[i] + z[i] <  junkOver) && 
	      (x[i] + y[i] + z[i] >  junkBelow) 
	      )
	    {
	      dd[i] = log((x[i] * u + y[i] * v)/z[i]) ;
	      dd[i] = dd[i] * dd[i] ;
	      ndd++ ;
	    }
	}
      for (ddd = 0, i = 0 ; i < iMax ; i++)
	ddd += dd[i] ;
      ddd /= ndd ;

      if (k == 0 || ddd < bestDd)
	{ bestDd = ddd ; bestf = u ; }
      fprintf (stderr, "\t%g", ddd) ;
    } 

  f = bestf ;
  fprintf (stderr, "\t%d\t%f\n", ndd, f) ; 
  ac_free (h) ;
  return f ;
} /* s2gSeqcAcedbGetF */

/*************************************************************************************/

static void s2gSeqcAcedb (S2G *s2g)
{
  AC_HANDLE h = ac_new_handle () ; 
  const char types[] = "ABCDE" ;
  int maxType  = strlen (types) ;
  int t ;
  Array *aa ;
  float ff[] = {1, 0, 0, 0, 0} ; 
  int NN = 50000 ;
  
  aa = (Array *) halloc (maxType * sizeof(Array), h) ;
  for(t = 0 ; t < maxType ; t++)
     aa[t] = arrayHandleCreate (NN, float, h) ;

  s2g->intronDict = dictHandleCreate (NN, h) ;
  aa = s2gSeqcAcedbParse (s2g, types, aa) ;

  for (t = 2 ; t < maxType ; t++)
    {
      int k ;
      float f = 0.5 , df = 0.1 ;

      for (k = 0 ; k < 8 ; k++)
	{
	  f = s2gSeqcAcedbGetF (s2g, aa, t, f, df) ;
	  df /= 5 ;
	  if (f < 0)
	    f = 5*df ;
	}
      ff[t] = f ;
    }

  if (1)
    {
      float x, bestx, e, beste ;

      printf ("Mix                      \tA\tE\tC\tD\tB\n") ;
      printf ("Observed                 \t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", 
	      ff[0], ff[4], ff[2], ff[3], ff[1]) ;
      
      beste = 1000 ;
      bestx = 1 ;
      for (x = 0.5 ; x < 1.5 ; x += .01)
	{
	  e = (4/(4+x) - ff[4])*(4/(4+x) - ff[4]) + (1/(1+x) - ff[2])*(1/(1+x) - ff[2]) + (1/(1+4*x) - ff[3])*(1/(1+4*x) - ff[3]) ;
	  if (e < beste)
	    { beste = e ; bestx = x ; }
	}
      x = bestx ;
      printf ("Estimate based on C value\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t x= %.2f, beste=%f\n", 
	      ff[0], 4/(4+x), 1/(1+x), 1/(1+4*x), ff[1], x, beste
	      ) ;
    }

  ac_free (h) ;
} /* s2gSeqcAcedb */


/*************************************************************************************/
/*************************************************************************************/

/* copied from aliqc : add up the counts in a tsv file */
static void s2gTsvAction (S2G *s2g, int action)
{
  AC_HANDLE h = ac_new_handle () ;
  DICT *tagDict = dictHandleCreate (256, h) ;
  DICT *runDict = action == 2 ? dictHandleCreate (256, h) : 0 ;
  Array tags = arrayHandleCreate (64, TT, h) ; 
 
  /* parse the tag file */
  s2gParseTsv (s2g, tagDict, runDict,  tags) ;
  /* export */
  switch (action)
    {
    case 1: /* merge */
      s2gCumulTsv (tags) ;
      s2gExportMergedTsv (s2g, tagDict, s2g->run,  tags) ;
      break ;
    case 2: /* view */
      s2gViewTsv (s2g, tagDict, runDict,  tags) ;
      break ;
    }
  ac_free (h) ;
  return ;
} /* s2gTsvAction  */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/
/* run -run NA12878MOD -laneList tmp/TSNP/NA12878MOD/LaneList -t 20 -target_fasta TARGET/CHROMS/hs.chrom_20.fasta.gz -t1 25000001 -t2 30010000 -minSnpFrequency 18 -minSnpCover 10 -minSnpCount 4 -target_class Z_genome -o tata -maxLanes 4
 */

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: sam2gold -g goldFile -i file_list \n"
	   "//      try: -h --help -v --version\n"
	   "// Parameters\n"
	   "//   --run runName         [short form : -r runName]\n"
	   "//   --method alignerName  [short form : -m alignerName]\n"
	   "//   --samStats : count aligned reads and bases and missmatches\n"
	   "//   -g mi:fi..gold:gold.cig\n"
	   "//     where mi is a method name and -i the corresponding sam file\n"
	   "//     The gold file  defines the gold standard\n"
	   "//     The expected format looks like\n"
	   "//        seq.1a	chr1	95672060	95672159	100M	95672060-95672159	+	TCCCTTTCACGCCTCTTCT\n"
	   "//   -s snp_file \n"
	   "//     Expect a snp file as exported by dna2dna -makeTest\n"
           "//        seq.1a   353 A>G\n"
           "//        seq.1a   402 InsATG\n"
           "//        seq.1a   430 DelCT\n"
	   "//     substitutions are evaluated if there is an MD:Z:error_list column\n"
	   "//     Otherwise only the indels are evaluated\n"
	   "//   -i file_list\n"
	   "//     example m1:f.sam,m2:g.sam,m3:h.sam\n"
	   "//     All sam files are analysed for each method\n"
	   "//   -u : unique\n"
	   "//     Restrict the analysis to uniquely aligned reads\n"
	   "//   -exportIntronSupport\n"
	   "//   -o prefix\n"
	   "//     all output files will be called prefix.*\n"
	   "//     otherwise they are exported to stdout\n"
  	   "//   -gzo : the output files should be gzipped\n"
	   "//   --anne: SEQC acedb titration (expects a single -i file)\n"
	   "//   --aliLn: export the aligned legth histogram as prefix.aliLn\n"
	   "// MERGE\n"
	   "//   -merge -run run_name\n"
	   "//     Merge tsv files, as exported by this program or by aliqc\n"
	   "//     The tsv format is very simple\n"
	   "//        tag run_name n x1 x2 x3 ... xn\n"
	   "//     Where n is the number of expected fields and all xi are integers\n"
	   "//     The merge functionality is to cumul the xi and attribute them\n"
	   "//     to the run_name provided on the command line\n"
	   "//     This method can collate counts, in case data were processed in parallel\n"
	   "// Optional debugging level\n"
	   "//   -exportIntronSupport  : show details\n"
	   "//   -debug : more diagnostics\n"
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
  S2G s2g ;
  AC_HANDLE h = 0 ;
  char commandBuf [4000] ;

  freeinit () ; 
  messErrorInit (argv[0]) ;

  h = ac_new_handle () ;
  memset (&s2g, 0, sizeof (S2G)) ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
  s2g.h = h ;
 
  if (argc < 2)
    usage (commandBuf, argc, argv) ;
  if (getCmdLineBool (&argc, argv, "-help") ||
      getCmdLineBool (&argc, argv, "--help")||
      getCmdLineBool (&argc, argv, "-h")
      )
    usage (commandBuf, 1, argv) ;
  if (getCmdLineBool (&argc, argv, "-v") ||
      getCmdLineBool (&argc, argv, "--version")
      )
    {
      fprintf (stderr, "sam2gold: %s\n", VERSION) ;
      exit (0) ;
    }
  s2g.gzo = getCmdLineBool (&argc, argv, "-gzo") ;
  s2g.unique = getCmdLineBool (&argc, argv, "-u") ;
  s2g.debug = getCmdLineBool (&argc, argv, "-debug") ;
  s2g.samStats = getCmdLineBool (&argc, argv, "--samStats") ;
  s2g.exportIntronSupport = getCmdLineBool (&argc, argv, "-exportIntronSupport") ;
  getCmdLineOption (&argc, argv, "-o", &(s2g.outFileName)) ;
  getCmdLineOption (&argc, argv, "-s", &(s2g.snpFileName)) ;
  getCmdLineOption (&argc, argv, "-g", &(s2g.goldFileName)) ;
  getCmdLineOption (&argc, argv, "-i", &(s2g.inFileList)) ;
  s2g.addReadPairSuffixForce = getCmdLineBool (&argc, argv, "-addReadPairSuffixForce") ;
  s2g.addReadPairSuffix2 = getCmdLineBool (&argc, argv, "-addReadPairSuffix2") ;
  s2g.addReadPairSuffix = getCmdLineBool (&argc, argv, "-addReadPairSuffix") ;
  s2g.merge = getCmdLineBool (&argc, argv, "-m") || getCmdLineBool (&argc, argv, "--merge") ;
  s2g.view = getCmdLineBool (&argc, argv, "-v") || getCmdLineBool (&argc, argv, "--view") ;
  s2g.exportAliLn = getCmdLineBool (&argc, argv, "-aliLn") || getCmdLineBool (&argc, argv, "--aliLn") ;
  getCmdLineOption (&argc, argv, "--run", &(s2g.run)) ;
  getCmdLineOption (&argc, argv, "-r", &(s2g.run)) ;
  getCmdLineOption (&argc, argv, "--method", &(s2g.method)) ;
  getCmdLineLong (&argc, argv, "--nRawReads", &(s2g.nRawReads)) ;
  getCmdLineLong (&argc, argv, "--nRawBases", &(s2g.nRawBases)) ;
  s2g.anne = getCmdLineBool (&argc, argv, "--anne") ;

  if (s2g.anne)
    {
      s2gSeqcAcedb (&s2g) ;
      goto done ;
    }
  else if (s2g.merge)
    {
      if (! s2g.run)
	messcrash ("Missing paramater: --merge needs --run <run_name>, please try --help") ;
      s2gTsvAction (&s2g, 1) ;
    }
  else if (s2g.view)
    {
      if (! s2g.run)
	messcrash ("Missing paramater: --merge needs --run <run_name>, please try --help") ;
      s2gTsvAction (&s2g, 2) ;
    }
  else if (s2g.samStats)
    {
      if (! s2g.run)
	messcrash ("Missing paramater: --samStats needs --run <run_name>, please try --help") ;
      if (! s2g.method)
	messcrash ("Missing paramater: --samStats needs --method <method_name>, please try --help") ;
      if (! s2g.inFileList)
	messcrash ("Missing paramater: --i <sam_file of file_list_of_sam_files> , please try --help") ;
      fprintf (stderr, "//%s : Start \n"
	       , timeShowNow()
	       ) ;

      s2gInit (&s2g) ;
      s2gParse (&s2g) ;
    }

  else if (s2g.goldFileName)
    {
      if (! s2g.inFileList && ! s2g.goldFileName)
	messcrash ("missing -i and -g arguments, please try sam2gold --help\n") ;
      
      if (s2g.goldFileName && ! strchr(s2g.goldFileName, ':'))
	messcrash ("The -g %s argument does not contain a : separator, it should be of the form t1r1..gold:file_name, please try --help", s2g.goldFileName) ;

      fprintf (stderr, "//%s : Start \n"
	       , timeShowNow()
	       ) ;

      s2gInit (&s2g) ;
      s2gParse (&s2g) ;
      /* do NOT sort intronSupport before s2gAnalyze */
      s2gAnalyze (&s2g) ;

      if (s2g.exportIntronSupport)
	s2gExportIntronSupport (&s2g) ;
      if (s2g.exportAliLn)
	s2gExportAliLn (&s2g) ;
    }

  s2gExportIntrons (&s2g) ;
  aceDnaSetIlmJumper (1) ;  /* to please the optimized linker */


  if (s2g.methodDict)
    { 
      int mx ;
      messAllocMaxStatus (&mx) ; 
      fprintf (stderr, "// %s done, %d files %d reads analyzed, max memory %d Mb\n"
	       , timeShowNow()
	       , dictMax (s2g.methodDict), dictMax (s2g.seqDict)
	       , mx) ;
     }

  if (s2g.snpFileName)
    {
      ACEOUT ao = aceOutCreate (s2g.outFileName, ".snp", FALSE, s2g.h) ;
      const char *mm = s2g.method ? s2g.method : "Aligner" ;
      const char *run = s2g.run ? s2g.run : "Run" ;
      int i, *g, *s, *t ;
      float p[4], r[4], f[4] ;

      aceOutDate (ao, "###", "Precision, recall and F-score of substitutions, insertions and deletions calling in the edited RefSeq benchmark") ;
      aceOutf (ao, "\n# Run\tMethod") ;
      aceOutf (ao, "\tSubstitutions in Truth\tSubstitutions called in method") ;
      aceOutf (ao, "\tFP: Substitutions false positives\tTP: Substitutions true positives\tFN: Substitutions false negatives") ;
      aceOutf (ao, "\tSubstitutions precision p=TP/(TP+FP)\tSubstitutions recall r=TP/(TP+FN)\tSubstitutions F-score f=2pr/(p+r)") ;
 
      aceOutf (ao, "\t\tRun\tMethod") ;
      aceOutf (ao, "\tInsertions in Truth\tInsertions called in method") ;
      aceOutf (ao, "\tFP: Insertions false positives\tTP: Insertions true positives\tFN: Insertions false negatives\tInsertions precision p=TP/(TP+FP)\tInsertions recall r=TP/(TP+FN)\tInsertions F-score f=2pr/(p+r)") ;
      aceOutf (ao, "\tSingle insertions precision\tDouble insertions precision\tTriple insertions precision") ;
      aceOutf (ao, "\tSingle insertions recall\tDouble insertions recall\tTriple insertions recall") ;
      aceOutf (ao, "\tSingle insertions F-score\tDouble insertions F-score\tTriple insertions F-score") ;
 
      aceOutf (ao, "\t\tRun\tMethod") ;
      aceOutf (ao, "\tDeletions in Truth\tDeletions called in method") ;
      aceOutf (ao, "\tFP: Deletions false positives\tTP: Deletions true positives\tFN: Deletions false negatives\tDeletions precision p=TP/(TP+FP)\tDeletions recall r=TP/(TP+FN)\tDeletions F-score f=2pr/(p+r)") ;
      aceOutf (ao, "\tSingle deletions precision\tDouble deletions precision\tTriple deletions precision") ;
      aceOutf (ao, "\tSingle deletions recall\tDouble deletions recall\tTriple deletions recall") ;
      aceOutf (ao, "\tSingle deletions F-score\tDouble deletions F-score\tTriple deletions F-score") ;
 
      aceOutf (ao, "\n%s\t%s", run, mm) ;
      g = s2g.nGoldSub ; s = s2g.nSub ; t = s2g.nSubOk ;
      for (i = 0 ; i < 4 ; i++)
	{
	  p[i] = s[i] ? (100.0 * t[i]) / s[i] : 0 ;
	  r[i] = g[i] ? (100.0 * t[i]) / g[i] : 0 ;
	  f[i] = (p[i] + r[i] > 0) ? 2 * p[i] * r[i] / (p[i] + r[i]) : 0 ;
	}
      aceOutf (ao, "\t%d\t%d", g[0], s[0]) ;
      aceOutf (ao, "\t%d\t%d\t%d", s[0] - t[0], t[0], g[0] - t[0]) ;
      aceOutf (ao, "\t%.2f\t%.2f\t%.2f", p[0], r[0], f[0]) ;
 
      aceOutf (ao, "\t\t%s\t%s", run, mm) ;
      g = s2g.nGoldIns ; s = s2g.nIns ; t = s2g.nInsOk ;
      for (i = 0 ; i < 4 ; i++)
	{
	  p[i] = s[i] ? (100.0 * t[i]) / s[i] : 0 ;
	  r[i] = g[i] ? (100.0 * t[i]) / g[i] : 0 ;
	  f[i] = (p[i] + r[i] > 0) ? 2 * p[i] * r[i] / (p[i] + r[i]) : 0 ;
	}
      aceOutf (ao, "\t%d\t%d", g[0], s[0]) ;
      aceOutf (ao, "\t%d\t%d\t%d", s[0] - t[0], t[0], g[0] - t[0]) ;
      aceOutf (ao, "\t%.2f\t%.2f\t%.2f", p[0], r[0], f[0]) ;
      aceOutf (ao, "\t%.2f\t%.2f\t%.2f", p[1], p[2], p[3]) ;
      aceOutf (ao, "\t%.2f\t%.2f\t%.2f", r[1], r[2], r[3]) ;
      aceOutf (ao, "\t%.2f\t%.2f\t%.2f", f[1], f[2], f[3]) ;

      aceOutf (ao, "\t\t%s\t%s", run, mm) ;
      g = s2g.nGoldDel ; s = s2g.nDel ; t = s2g.nDelOk ;
      for (i = 0 ; i < 4 ; i++)
	{
	  p[i] = s[i] ? (100.0 * t[i]) / s[i] : 0 ;
	  r[i] = g[i] ? (100.0 * t[i]) / g[i] : 0 ;
	  f[i] = (p[i] + r[i] > 0) ? 2 * p[i] * r[i] / (p[i] + r[i]) : 0 ;
	}
      aceOutf (ao, "\t%d\t%d", g[0], s[0]) ;
      aceOutf (ao, "\t%d\t%d\t%d", s[0] - t[0], t[0], g[0] - t[0]) ;
      aceOutf (ao, "\t%.2f\t%.2f\t%.2f", p[0], r[0], f[0]) ;
      aceOutf (ao, "\t%.2f\t%.2f\t%.2f", p[1], p[2], p[3]) ;
      aceOutf (ao, "\t%.2f\t%.2f\t%.2f", r[1], r[2], r[3]) ;
      aceOutf (ao, "\t%.2f\t%.2f\t%.2f", f[1], f[2], f[3]) ;

      aceOutf (ao, "\n") ;


    }
 done:
  ac_free (s2g.h) ;

  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderrn and to ensure all pipes are closed*/
  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
