/*  File: sambam.c
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

 * created 2019_05

 * A technical library to
 * manipulate the sam/bam file format
 * in particular, roll and smoke cigars
 */

#define VERSION "1.1"
#include "ac.h"
#include "acedna.h"

/****************/
/****************/
/* BAM/SAM cigar format, a string describing an alignemnt as a sequence of length-tag words
 * example 23S40M1X3M2D20M1200N35M6S

 ** full Cigar spec as of 2019
 * M  match or mismatch
 * I insertion to the reference
 * D deletion from the reference (at least in magicBlast, del > 12 bases  should be removed from the error count NM:i:1234)
 * N skipped region from the reference (intron)
 * S soft clipping (move in read)
 * H hard clipping (not reported in the read-dna)
 * P padding (silent deletion from padded reference) 
 * = sequence match
 * sequence mismatch
 */
static int samDoParseCigar (char *cigar, Array cigarettes, int a1, int *a2p, int *x1p, int *x2p, int *alip, int *insp, int *delp, int *bigDelp, int *clipp) 
{
  int ii = 0, x1 = 1, x2 = 0 ;
  int a2 = a1, ali = 0, clip1 = 0, clip2 = 0, hClip = 0 ;
  int score = 0 ;
  char cc, *cp = cigar ;
  SAMCIGAR *cgr ;
  
  while (*cp)
    {
      char c = *cp ;
      int dx = 0 ;
      while (c >= '0' && c <= '9')
	{
	  dx = 10 * dx + (c - '0') ;
	  c = *++cp ;
	}
      cc = c ; cp++ ;
      if (dx <= 0)
	dx =  0 ;
      switch ((int)cc)
	{
	case 'H': /* hard clip */
	  hClip += dx ;
	case 'S': /* soft clip */
	  a2 = a1 - 1 ;
	  x1 += dx ;
	  x2 = x1 - 1 ;
	  if (ali == 0)
	    clip1 += dx ;
	  else
	    clip2 += dx ;
	  break ;
	case '=':
	  score += dx ;
	case 'M':
	case 'X':
	  a2 = a1 + dx - 1 ;
	  x2 = x1 + dx - 1 ;
	  /*
	    if (! *x1p)
	    *x1p = x1 ;
	  if (x2p)
	    *x2p = x2 ;
	    */
	  ali += dx ;
	  break ;
	case 'I':
	  a2 = a1 - 1 ; /* so a1 will not change */
	  x2 = x1 + dx - 1 ;
	  if (insp) *insp += dx ;
	  score -= 4 * (dx > 3 ? 3 : dx) ;
	  break ;
	case 'D':
	  x2 = x1 - 1 ; /* so x1 will not change */
	  a2 = a1 + dx - 1 ;
	  if (delp) *delp += (dx < 12 ? dx : 0) ;
	  if (bigDelp) *bigDelp += (dx >= 12 ? dx : 0) ;
	  score -= 4 * (dx > 3 ? 3 : dx) ;
	  break ;
	case 'N': /* gt_ag intron */
	case 'n': /* gc_ag intron */
	  x2 = x1 - 1 ; /* so x1 will not change */
	  a2 = a1 + dx - 1 ;
	  break ;
	case 'P': /* padding, ignore */
	  continue ;
	default:
	  messcrash ("Unknown character %c in cigar string %s", cc, cigar) ;
	}

      cgr = arrayp (cigarettes, ii++, SAMCIGAR) ;
      cgr->x1 = x1 ;
      cgr->x2 = x2 ;
      cgr->a1 = a1 ;
      cgr->a2 = a2 ;
      cgr->type = cc ;
      cgr->dx = dx ;
      a1 = a2 + 1 ;
      x1 = x2 + 1 ;
    }
  arrayMax (cigarettes) = ii ;

  if (x1p)  *x1p = 1 + clip1 ;
  if (x2p) *x2p = x2 - clip2 ;
  if (a2p) *a2p = a2 ;
  if (alip) *alip = ali ;
  if (clipp) *clipp = clip1 + clip2 - hClip ;

  return score ;
} /*  samDoParseCigar */

/****************/

BOOL samParseCigar (char *cigar, Array cigarettes, int a1, int *a2p, int *x1p, int *x2p, int *alip, int *insp, int *delp, int *bigDelp) 
{
  int clip = 0 ;

  if (! arrayExists (cigarettes))
    messcrash ("dnaParseCigar passed a null cigarettes array") ;
  if (cigarettes->size != sizeof (SAMCIGAR))
    messcrash ("dnaParseCigar passed a cigarettes array which is not of type SAMCIGAR") ;
  return samDoParseCigar (cigar, cigarettes, a1, a2p, x1p, x2p, alip, insp, delp, bigDelp, &clip)  ;
} /*  samParseCigar */

/****************/

int samScoreCigar (const char *readName, char *cigar, Array cigarettes, int a1, int ln)
{
  int score, x1 = 1, x2 = 0, a2 = 0, ali = 0, clip = 0, ins = 0, del = 0, bigDel = 0 ;
  
  score = samDoParseCigar (cigar, cigarettes, a1, &a2, &x1, &x2, &ali, &ins, &del, &bigDel, &clip) ;
  
  if (ln != clip + x2 - x1 + 1)
    score = -1 ;
  return score ; 
} /* samScoreCigar */

/****************/

BOOL samCheckCigar (const char *readName, char *cigar, Array cigarettes, int a1, int ln)
{
  int x1 = 1, x2 = 0, a2 = 0, ali = 0, clip = 0, ins = 0, del = 0, bigDel = 0 ;
  samDoParseCigar (cigar, cigarettes, a1, &a2, &x1, &x2, &ali, &ins, &del, &bigDel, &clip) ;
  
  if (ln != clip + x2 - x1 + 1)
    {
      messerror ("samCheckCigar failed on %s ln=%d dx=%d cigar=%s"
	       , readName ? readName : "NULL"
	       , ln, clip + x2 - x1 + 1, cigar) ;
      return FALSE ;
    }
  return TRUE ; 
} /* samCheckCigar */

/*************************************************************************************/
/* Parse one sam line */
#ifdef JUNK
static BOOL samParseLine (ACEIN ai, SAMHIT *up, DICT *dict, int target_class, int nn)
{ 
  char *ccp ;
  char buf[1000] ;
  int Z_genome = 0 ;
  HIT2 *up2 ;
  HIT3 *up3 ;
  int flag ;
  char *cigar ;



  memset (up, 0, sizeof (SAMHIT)) ;

  ccp = aceInWord (ai) ;
  if (! ccp || ! *ccp || *ccp == '#' || *ccp == '@')
    goto done  ;

  ccp = aceInWord (ai) ;
  dictAdd (dict,ccp, &(up->tag)) ;
  up->mult = 1 ; /* default */
  
  if ((ccp = strchr (ccp, '#')))
    { /* multiplicity */
      int m = 1 ;
      char cc ;

      if (sscanf (ccp+1, "%d%c", &m, &cc) == 1 && m > 0)
	up->mult = m ;
    }
  aceInStep (ai, '\t') ;   aceInInt (ai, &flag) ;
 
  /*
    int isFirstFragment = 1 ;
    if (flag & 128) isFirstFragment = -1 ;
  */

  aceInStep (ai, '\t') ; ccp = aceInWord (ai) ;
  if (ccp && !strncmp (ccp,"MRNA:",5))ccp += 5 ;
  if (! ccp || ! *ccp)
    goto done  ;
  dictAdd (dict, ccp, &(up->target)) ;
  
  up->a1 = 0 ;
  aceInStep (ai, '\t') ; aceInInt (ai, &(up->a1)) ;
  if (up->a1)
    goto done  ;
 
  aceInStep (ai, '\t') ; aceInInt (ai, &(up2->a1)) ;
  aceInStep (ai, '\t') ;  cigar = aceInWord (ai) ;  /* quality, discard */
  samParseCigar (cigar, cigarettes, up2->a1, &(up2->a2), &(up->x1), &(up->x2), 0, 0, 0) ;

  up->target_class = target_class ;
  aceInStep (ai, '\t') ;  aceInWord (ai) ; /* drop */
  aceInStep (ai, '\t') ;  aceInWord (ai) ; /* drop */
  aceInStep (ai, '\t') ;  aceInInt (ai, &(up->dPair)) ;
  aceInStep (ai, '\t') ;  ccp = aceInWord (ai) ;  /* the actual sequence */
  if (!ccp)
    goto done ;
  up->ln = strlen (ccp) ;

  while (ccp)
    {
      aceInStep (ai, '\t') ;  ccp = aceInWord (ai) ;  /* the actual sequence */
      if (!strncmp (ccp, "AS:", 3))
	sscanf (ccp+5, "%d", &(up->score)) ;
    }

 done:
  ac_free (h) ;  
  return TRUE ; /* true i read a line, bud bad lineif sh->score == 0 */
} /* samParseLine */
#endif

/*************************************************************************************/
/* check all cigars in a sam file without storing or copying the data
 */
static int samFileAnalyze (ACEIN ai, ACEOUT ao
			   , const char *target0
			   , DICT *targetDict, Array targetDnas
			   , DICT *intronDict, KEYSET intronSupport
			   )
{ 
  int nSuccessfullChecks = 0 ;
  int iMax ;
  AC_HANDLE h = ac_new_handle () ;
  Array cigarettes = arrayHandleCreate (1024, SAMCIGAR, h) ;
  Array dnaShort = arrayHandleCreate (1024, char, h) ;
  const char *label = aceInLabel (ai) ; 
  vTXT newCigar = vtxtHandleCreate (h) ;
  vTXT  txt1 = ao ? vtxtHandleCreate (h) : 0 ;
  vTXT  txt2 = ao ? vtxtHandleCreate (h) : 0 ;
  char *fNam = 0 ;
  /* the target_class is supposed to be stored in the file label 
   * but it can be overridden by a column in the bam file
   */
  if (!ai)
    goto done ;
  
  
  if (! label) 
    label = "sam" ;
  aceInSpecial (ai, "\n") ;
  fNam =strnew (aceInFileName (ai), h) ;
  while (aceInCard (ai))
    {
      /* non reentrant: split the fields in aceInPos with copyiing the data */
      int n = 0 ;
      char *cp, *cq ;
      char 
	*readName = 0, 
	*target = 0, 
	*cigar = 0,
	*dna = 0
	; 
      int 
	flag = 0,
	a1 = 0,
	mult = 1
	;

      if (txt1)
	{ 
	  vtxtClear (txt1) ;
	  vtxtClear (txt2) ;
	}
      /* ATTENTION, in sam/bam, never use aceOutf (formatted)
       * because sam files are full of non standard characters
       * in particular they may contain % and quotes 
       */
      cp = aceInPos (ai) ;
      if (!cp || *cp == '#' || *cp == '@' || *cp == '/')
	{
	  if (ao && cp && *cp)
	    {
	      aceOut (ao, cp) ;
	      aceOut (ao, "\n") ;
	    }
	  continue ;
	}
      while (cp && *cp)
	{
	  cq = strchr (cp, '\t') ;
	  if (cq)
	    *cq++ = 0 ;
	  switch (++n)
	    {
	    case 1:
	      readName = cp ;
	      if (txt1)
		vtxtPrint (txt1, cp) ;
	      if ((cp = strchr (cp, '#')))
		{ /* multiplicity */
		  int m = 1 ;
		  char cc ;
		  
		  if (sscanf (cp+1, "%d%c", &m, &cc) == 1 && m > 0)
		    mult = m ;
		}
	      break ;
	    case 2:
	      if (txt1)
		{
		  vtxtPrint (txt1, "\t") ;
		  vtxtPrint (txt1, cp) ;
		}
	      flag = atoi (cp) ;
	      break ;
	    case 3:
	      if (txt1)
		{
		  vtxtPrint (txt1, "\t") ;
		  vtxtPrint (txt1, cp) ;
		}
	      target = cp ;
	      if (target0 && strcmp (target, target0))
		continue ;
	      break ;
	    case 4:
	      if (txt1)
		{
		  vtxtPrint (txt1, "\t") ;
		  vtxtPrint (txt1, cp) ;
		}
	      a1 = atoi (cp) ;
	      break ;
	    case 5:
	      /* do not export the mapping score yet */
	      break ;
	    case 6:
	      cigar = cp ;  /* do not export cigar before it is modified */
	      break ;
	    case 7:
	      if (txt2)
		{
		  vtxtPrint (txt2, "\t") ;
		  vtxtPrint (txt2, cp) ;
		}
	      break ;
	    case 8:
	      if (txt2)
		{
		  vtxtPrint (txt2, "\t") ;
		  vtxtPrint (txt2, cp) ;
		}
	      break ;
	    case 9:
	      if (txt2)
		{
		  vtxtPrint (txt2, "\t") ;
		  vtxtPrint (txt2, cp) ;
		}
	      break ;
	    case 10:
	      if (txt2)
		{
		  vtxtPrint (txt2, "\t") ;
		  vtxtPrint (txt2, cp) ;
		}
	      dna = cp ;
	      break ;
	    case 11:
	      if (txt2)
		{
		  vtxtPrint (txt2, "\t") ;
		  vtxtPrint (txt2, cp) ;
		}
	      break ;
	    default:
	      if (txt2)
		{
		  vtxtPrint (txt2, "\t") ;
		  vtxtPrint (txt2, cp) ;
		}
	      /*
		{
		  int nh, as, nm ;
		  char *xs ;
		  if (! strncmp (cp, "NH:i:", 5))
		    nh = atoi (cp + 5) ;
		  else if (! strncmp (cp, "AS:i:", 5))
		    as = atoi (cp + 5) ;
		  else if (! strncmp (cp, "NM:i:", 5))
		    nm = atoi (cp + 5) ;
		  else if (! strncmp (cp, "XS:A:", 5))
		    xs = cp + 5 ;
		}
	      */
	      break ;
	    }
	  cp = cq ;
	}
      if (!a1 || !cigar[1] || !dna[1])
	continue ;
      if (! samCheckCigar (readName, cigar, cigarettes, a1, strlen (dna)))
	messcrash ("Bad cigar %s %s\n", readName, cigar) ;
      nSuccessfullChecks++ ;

      /* register intron support */
      iMax = arrayMax (cigarettes) ;
      if (iMax && target && intronDict && strlen (target) < 856)
	{ /* export introns */
	  int i, k, t ;
	  SAMCIGAR *cgr ;

	  for (i = 0, cgr = arrp (cigarettes, 0, SAMCIGAR) ; i < iMax ; cgr++, i++)
	    {
	      if (cgr->type == 'N')
		{
		  char buf[1024] ;
		  dictAdd (targetDict, target, &t) ;
		  sprintf (buf, "%d\t%d\t%d", t, cgr->a1, cgr->a2) ;
		  dictAdd (intronDict, buf, &k) ;
		  keySet (intronSupport, k) += mult ;
		}
	    }
	}
      
      /* score the substitutions */
      vtxtClear (newCigar) ;
      if (iMax && target && targetDnas &&  targetDict)
	{ /* score the substitutions and indels */
	  int ii, i, dx, t, score = 0, ln, ali = 0 ;
	  int y1, y2, b1, b2 ;
	  SAMCIGAR *cgr ;
	  Array dnaLong = 0 ;
	  int errCost = 4 ;
	  char *cp, *cq ;

	  if (dictFind (targetDict, target, &t) &&
	      t < arrayMax (targetDnas)
	      )
	    dnaLong = arr (targetDnas, t, Array) ;
	  if (dnaLong)
	    {
	      ln = strlen (dna) ;
	      array (dnaShort, ln , char) = 0 ;
	      arrayMax (dnaShort) = ln ;
	      memcpy (arrp (dnaShort, 0, char), dna, ln) ;
	      dnaEncodeArray (dnaShort) ;
	      for (ii = 0, cgr = arrp (cigarettes, 0, SAMCIGAR) ; ii < iMax ; cgr++, ii++)
		{
		  dx = cgr->dx ;
		  if (ii == 0)
		    { b1 = cgr->a1 ; y1 = cgr->x1 ; }
		  if (1)
		    { b2 = cgr->x2 ; y2 = cgr->x2 ; }
		  
		  switch (cgr->type)
		    {
		    case 'I':
		    case 'D':
		      score -= (dx < 3 ? dx : 3) * errCost ;
		      vtxtPrintf (newCigar, "%d%c", dx, cgr->type) ;
		      break ;
		    case '=':
		      ali += dx ;
		      vtxtPrintf (newCigar, "%d%c", dx, cgr->type) ;
		      score += dx ;
		      break ;
		    case 'X':
		      ali += dx ;
		      vtxtPrintf (newCigar, "%d%c", dx, cgr->type) ;
		      score -= dx * errCost ;
		      break ;
		    case 'M':  /* here is the rub */
		      ali += dx ;
		      cp = arrp (dnaShort, cgr->x1 - 1, char) ;
		      cq = arrp (dnaLong, cgr->a1 - 1, char) ;
		      {
			int m = 0, x = 0 ;
			for (i = 0 ; i < dx ; i++, cp++, cq++)
			  {
			    if ((*cp) & (*cq))
			      {
				if (x)
				  {
				    vtxtPrintf (newCigar, "%dX", x) ;
				    x = 0 ; m = 1 ;
				  }
				else
				  m++ ;
				score++ ;
			      }
			    else
			      {
				if (m)
				  {
				    vtxtPrintf (newCigar, "%d=", m) ;
				    m = 0 ; x = 1 ;
				  }
				else
				  x++ ;
				score -= errCost ;
			      }
			  }
			if (m)
			  vtxtPrintf (newCigar, "%d=", m) ;
			if (x)
			  vtxtPrintf (newCigar, "%dX", x) ;
		      }				  
		      break ;
		    case 'n':
		    case 'N':
		    case 'S':
		    case 'H':
		      vtxtPrintf (newCigar, "%d%c", dx, cgr->type) ;
		      break ;
		    }
		}
	      if (! samCheckCigar (readName, vtxtPtr (newCigar), cigarettes, a1, strlen (dna)))
		messcrash ("Bad newCigar %s %s\n", readName, vtxtPtr (newCigar)) ;
	      if (ao)
		{  
		  aceOut (ao, vtxtPtr (txt1)) ;
		  aceOutf (ao, "\t%d\t%s", score, vtxtPtr (newCigar)) ;
		  aceOut (ao, vtxtPtr (txt2)) ;
		  aceOutf (ao, "\toc:A:%s", cigar) ;
		  aceOut (ao, "\n") ;
		  a1 = a1 + flag + b1 + b2 + y1 + y2 ; /* to please the compiler */
		}
	    }
	}      
    }
  
 done:

  fprintf (stderr, "# Successfully checked %d cigars in file %s\n"
	   , nSuccessfullChecks, fNam
	   ) ;
  
  ac_free (h) ;	

  return nSuccessfullChecks ;
} /* samFileAnalyze */

/*************************************************************************************/
/* check all cigars in a sam file without storing or copying the data
 */
int samFileCheckCigars (ACEIN ai, const char *target) 
{
  ACEOUT ao = 0 ;
  return samFileAnalyze (ai, ao, target, 0, 0, 0, 0) ;
} /* samFileCheckCigars */

/*************************************************************************************/
/* check and export the list of introns
 */
int samFileExportIntrons (ACEIN ai, const char *target, DICT *targetDict, DICT *intronDict, KEYSET intronSupport)
{
  ACEOUT ao = 0 ;
  return samFileAnalyze (ai, ao, target, targetDict, 0, intronDict, intronSupport) ;
} /* samFileExportIntrons */

/*************************************************************************************/
/* Export sam file with scores, requires the genome DNA
 */
int samFileScore (ACEIN ai, ACEOUT ao, const char *target, DICT *targetDict, Array targetDnas, DICT *intronDict, KEYSET intronSupport)
{
  return samFileAnalyze (ai, ao, target, targetDict, targetDnas, intronDict, intronSupport) ;
} /* samFileScoreScoredSamFile  */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
