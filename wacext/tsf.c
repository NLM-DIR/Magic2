/*  File: tsf.c
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
 *
 * 2018_04_01
 * tsf   'table serialised format' 
 * strategic idea
 * Define a format semantically equivalent to a 2 dimensional table format
 * and interchangeable with a table format such that each line in the file
 * is providing the value of a single cell, or a small set of contiguous cells.
 * The order of the lines is immaterial, so they can be produced by the source
 * program emitting the data in random order.
 * When 2 values are associated to the same cell, the program applies
 * an action: update (take last value), add (sum the values), 
 * juxtapose (a then b produces a:b), merge (a then a produces a)
 * compute (a then b produces 7*sin(a) + cos(b)), only in case of 
 * parsing 2 files
 *
 * The format has 3 + n fields
 *     1: T tag  (understood as an object or the nmae of the line of a table)
 *     2: S sample  (understood as a city, where the tag was evaluated, particular to the caption of a column lime)
 *     3: F format (descibes the content of a cell or cell group)
 *   ...: v values conforming to the format
 *
 * There is a similarity with the vcf format (variant cll fromat) used in genetics
 * but the objective and the data are manipulated are different.
 *
 * For a human readble presentation, use tsf -I tsf -O table, to export as table
 */

#define VERSION "1.2"
#include "ac.h"
#include "topology.h"

/* snp chaining distance */

typedef struct tsf_struct {
  AC_HANDLE h ;
  
  const char *inFileList ;
  const char *outFileName ;
  const char *inFileOfFileList ;
  const char *separator ;
  const char *title ;
  const char *corner_title ;
  const char *tag ;
  const char *sample ;
  const char *sortCommand ;
  const char *sampleSelectList ;
  const char *sampleRenameList ;
  const char *tagSelectList ;
  const char *compute ;

  const char *valueSeparator ;  /* default , */
  const char *fieldSeparator ;  /* default \t */
  const char *NA ;              /* default "" */

  
  DICT *sampleDict ;
  DICT *tagDict ;
  DICT *sampleOldNameDict ;
  DICT *sampleNewNameDict ;
  DICT *sampleSelectDict ;
  DICT *tagSelectDict ;
  DICT *tableCaptionDict ;
  DICT *cellCaptionDict ;
  DICT *columnCaptionDict ;
  DICT *valDict ;
  KEYSET sampleRenameKs ;
  KEYSET skipKs ;
  KEYSET sampleDepth ;
  KEYSET col2sample ;
  AC_TABLE results ;

  BOOL gzi, gzo, transpose, debug ;
  BOOL isTableIn, isTableOut ;
  BOOL merge, noMerge ;
  BOOL sum, min, max, replace, sumAll ; 
  BOOL sampleSelectOnly ;
  BOOL tagSelectOnly ;
  BOOL chronoOrder ;
  int skip1, skip2  ;
  int tagDepth ;

  ACEIN ai ; 
  ACEOUT ao ;
  BigArray hits ;
  int  nInputFiles ;
} TSF ;

#define MAXSAMPLEDEPTH 32
#define MAXCOL 64 

typedef struct hit_struct {
  int sample, tag, n ;
  unsigned int flag ;  /* data existence */
  char types[MAXSAMPLEDEPTH] ;
  long int x[MAXSAMPLEDEPTH] ;
  double z[MAXSAMPLEDEPTH] ;
} HIT ;


#ifdef RESTRICT 
#define Restrict restrict
#else
#define Restrict 
#endif
static void tsfParseMetaData (TSF *tsf, ACEIN ai) ;
static BOOL tsfMergeOne (HIT *hit, HIT *hit2, int action) ;

/*************************************************************************************/
static DICT *myDict ;
static int tsfTagOrder (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va ;
  const HIT *vp = (const HIT *)vb ;
  int n, s1, s2 ; ;

  s1 = up->tag ;
  s2 = vp->tag ;
  n = lexstrcmp (dictName (myDict, s1), dictName(myDict,s2)) ; if (n) return n ;
  n = up->sample - vp->sample ; if (n) return n ;

  return 0 ;
} /* tsfTagOrder */

/*************************************************************************************/

static int tsfSampleOrder (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va ;
  const HIT *vp = (const HIT *)vb ;
  int n ;

  n = up->sample - vp->sample ; if (n) return n ;
  n = up->tag - vp->tag ; if (n) return n ;

  return 0 ;
} /* tsfSampleOrder */

/*************************************************************************************/
/*************************************************************************************/

static void tsfTagSelectInit (TSF *tsf)
{
  AC_HANDLE h = ac_new_handle() ;
  char *cq, *cp = strnew (tsf->tagSelectList, h) ;
 
  tsf->tagSelectOnly = TRUE ;
  tsf->tagSelectDict = dictHandleCreate (256, tsf->h) ;
  while (cp)
    {
      cq = strchr (cp, ',') ;
      if (cq) *cq++ = 0 ;
      if (!cp[0])
	{
	  fprintf (stderr, "parameter --tagSelect cannot be parsed, try --help,  expecting comma delimited list a,b,c and so on, received %s"
		 , tsf->tagSelectList) ;
	  exit (1) ;
	}

      if (strchr (cp, '*'))
	tsf->tagSelectOnly = FALSE ;
      else
	{
	  dictAdd (tsf->tagSelectDict, cp, 0) ;
	  dictAdd (tsf->tagDict, cp, 0) ;
	}
      cp = cq ;
    }
  ac_free (h) ;

} /* tsfTagSelectInit */

/*************************************************************************************/
/* create a new sample and associated data */
static int sampleAdd (TSF *tsf, const char *ccp)
{
  int k ;

  if (tsf->sampleSelectOnly && tsf->sampleSelectDict 
      && ! dictFind (tsf->sampleSelectDict, ccp, 0)
      )
    return 0 ;

  dictAdd (tsf->sampleDict, ccp, &k) ;
  return k ;
}  /* sampleAdd */

/*************************************************************************************/

static void tsfSampleSelectInit (TSF *tsf)
{
  AC_HANDLE h = ac_new_handle() ;
  char *cq, *cp = strnew (tsf->sampleSelectList, h) ;
  BOOL sampleSelectOnly = TRUE ;
  tsf->sampleSelectOnly = FALSE ;
  
  tsf->sampleSelectDict = dictHandleCreate (256, tsf->h) ;
  while (cp)
    {
      cq = strchr (cp, ',') ;
      if (cq) *cq++ = 0 ;
      if (!cp[0])
	{
	  fprintf (stderr, "// --sampleSelect ERROR, try --help,  expecting comma delimited list a,b,c and so on, received %s"
		 , tsf->sampleSelectList) ;
	  exit (1) ;
	}

      if (strchr (cp, '*'))
	sampleSelectOnly = FALSE ;
      else
	{
	  dictAdd (tsf->sampleSelectDict, cp, 0) ;
	  sampleAdd (tsf, cp) ;
	}
      cp = cq ;
    }
  tsf->sampleSelectOnly = sampleSelectOnly ;
  ac_free (h) ;

} /* tsfSampleSelectInit */

/*************************************************************************************/

static void tsfSampleRenameInit (TSF *tsf)
{
  AC_HANDLE h = ac_new_handle() ;
  char *cr, *cq, *cp = strnew (tsf->sampleRenameList, h) ;
 
  tsf->sampleOldNameDict = dictHandleCreate (256, tsf->h) ;
  tsf->sampleNewNameDict = dictHandleCreate (256, tsf->h) ;
  tsf->sampleRenameKs = keySetHandleCreate (tsf->h) ;
  while (cp)
    {
      int old, new ;

      cq = strchr (cp, ',') ;
      if (cq) *cq++ = 0 ;
      
      cr = strchr (cp, ':') ;
      if (!cr || ! cr[1] || cp[0] == ':')
	{
	  fprintf (stderr, "// --sampleRename ERROR, try --help,  expecting a1:b1,a2:b2 found %s in string %s"
		 , cp, tsf->sampleRenameList) ;
	  exit (1) ;
	}
      *cr++ = 0 ;
      if (! strncmp (cr, "_P_", 3)) /* do not parse percensamplees */
	{
	  fprintf (stderr, "// --sampleRename ERROR, _P_ is reserved, you cannot rename %s to %s", cp, cr)  ;
	  exit (1) ;
	}

      dictAdd (tsf->sampleOldNameDict, cp, &old) ;
      dictAdd (tsf->sampleNewNameDict, cr, &new) ;
      keySet (tsf->sampleRenameKs, old) = new ;

      cp = cq ;
    }
  ac_free (h) ;
} /* tsfSampleSelectInit */

/*************************************************************************************/

static KEYSET tsfMyCols = 0 ;


static void tsfInit (TSF *tsf)
{
  tsfMyCols = keySetCreate () ;
  tsf->sampleDict = dictHandleCreate (256, tsf->h) ;
  tsf->tagDict = dictHandleCreate (10000, tsf->h) ;
  tsf->hits = bigArrayHandleCreate (100000, HIT, tsf->h) ;
  tsf->tableCaptionDict = dictHandleCreate (256, tsf->h) ;
  tsf->cellCaptionDict = dictHandleCreate (256, tsf->h) ;
  tsf->columnCaptionDict = dictHandleCreate (256, tsf->h) ;
  tsf->valDict = dictHandleCreate (256, tsf->h) ;
  tsf->sampleDepth = keySetHandleCreate (tsf->h) ;
  tsf->col2sample = keySetHandleCreate (tsf->h) ;
  if (tsf->sample)
    dictAdd (tsf->sampleDict, tsf->sample, 0) ;
  if (tsf->tag)
    dictAdd (tsf->tagDict, tsf->tag, 0) ;
  if (tsf->sampleRenameList)
    tsfSampleRenameInit (tsf)  ;
  if (tsf->sampleSelectList)
    tsfSampleSelectInit (tsf)  ;
  if (tsf->tagSelectList)
    tsfTagSelectInit (tsf)  ;
  return  ;
} /* tsfInit */

/*************************************************************************************/

static long int tsfParseTable (TSF *tsf, ACEIN ai)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, kk, tag = 0 ;
  const char *ccp, *ccq ;
  char cutter ;
  HIT *hit ;
  /*
  KEYSET sampleDepth = tsf->sampleDepth ; 
  DICT *sampleDict = tsf->sampleDict ;
  */
  DICT *tagDict = tsf->tagDict ;
  BigArray hits = tsf->hits ;
  long int iMax = bigArrayMax (hits) ;
  int tagDepth = 1 ;
  vTXT tagTxt = vtxtHandleCreate (h) ;
  KEYSET samplesKs = keySetHandleCreate (h) ;
  KEYSET nKs = keySetHandleCreate (h) ;
  KEYSET posKs = keySetHandleCreate (h) ;


  if (tsf->tag)
    {
      ccp = tsf->tag ;
      ccq = strstr (ccp, "..") ;
      tagDepth = 1 ;
      while (ccq)
	{
	  tagDepth++ ;
	  ccq = strstr (ccp + 2, "..") ;
	}
      tsf->tagDepth = tagDepth ;	
      dictAdd (tagDict, ccp, &tag) ;
    }


  while (aceInCard (ai))
    {
      ccp = aceInWordCut (ai, "\t", &cutter) ;
      if (! ccp)
	continue ;
      if (!strncmp (ccp, "#", 1))
	{
	  aceInCardBack (ai) ;
	  tsfParseMetaData (tsf, ai) ;
	  continue ;
	}


      /* parse the tags */ 
      vtxtClear (tagTxt) ;
      for (i = 0 ; i < tagDepth ; i++)
	{
	  if (i>0)
	    vtxtPrintf (tagTxt, "..") ;
	  vtxtPrintf (tagTxt, "%s", ccp) ;
	  ccp = aceInWordCut (ai, "\t", &cutter) ;
	}
      
      if (tsf->tagSelectOnly && tsf->tagSelectDict 
	  && ! dictFind (tsf->tagSelectDict, ccp, 0)
	  )
	continue ;
      dictAdd (tagDict, vtxtPtr (tagTxt), &tag) ;
      
      kk = 0 ;
      while (aceInStep (ai, '\t')) 
	{
	  kk++ ;
	  ccp = aceInWordCut (ai, "\t", &cutter) ;
	  if (ccp)
	    {
	      char cc ;
	      long int x ;
	      if (sscanf (ccp, "%ld%c", &x, &cc) == 1)
		{
		  int pos = keySet (posKs, kk) ;
		  hit = bigArrayp (hits, iMax++, HIT) ;
		  hit->sample = keySet (samplesKs, kk) ;
		  hit->tag = tag ;
		  hit->n = keySet (nKs, kk) ;
		  hit->x[pos] = x ;
		}
	    }
	}
    }
  ac_free (h) ;
  return 0 ;
}/*  tsfParseTable */

/*************************************************************************************/

static long int tsfParseTsf (TSF *tsf, ACEIN ai)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, sample, tag = 0, line = 0 ;
  int action = 0 ;
  long int x ;
  double z ;
  const char *ccp, *ccq ;
  HIT *hit ;
  KEYSET sampleDepth = tsf->sampleDepth ; 
  DICT *valDict = tsf->valDict ;
  DICT *tagDict = tsf->tagDict ;
  BigArray hits = tsf->hits ;
  long int iMax = bigArrayMax (hits) ;
  int tagDepth, nFields ;
  char typeBuf[MAXSAMPLEDEPTH+1] ;

  if (tsf->tag)
    {
      ccp = tsf->tag ;
      ccq = strstr (ccp, "..") ;
      tagDepth = 1 ;
      while (ccq)
	{
	  tagDepth++ ;
	  ccq = strstr (ccp + 2, "..") ;
	}
      tsf->tagDepth = tagDepth ;	
      dictAdd (tagDict, ccp, &tag) ;
    }

  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      char cc ;
      ccp = aceInWordCut (ai, "\t", &cc) ;
      line++ ;
      if (! ccp)
	continue ;
      if (!strncmp (ccp, "#", 1))
	{
	  aceInCardBack (ai) ;
	  tsfParseMetaData (tsf, ai) ;
	  continue ;
	}


      /* parse the line/tag */
     if (tsf->tagSelectOnly && tsf->tagSelectDict 
	 && ! dictFind (tsf->tagSelectDict, ccp, 0)
	 )
       continue ;
     if (! tsf->tag)
       {
	 ccq = strstr (ccp, "..") ;
	 tagDepth = 1 ;
	  while (ccq)
	    {
	      tagDepth++ ;
	      ccq = strstr (ccq + 2, "..") ;
	    }
	  if (tagDepth > tsf->tagDepth)
	    tsf->tagDepth = tagDepth ;
	  if (tagDict)  /* tag_name , ignore in merge case */
	    {
	      static int nnn = 1 ;
	      if (tsf->isTableOut && tsf->noMerge)
		ccp = hprintf (h, "%s#_#_#%d", ccp, nnn++) ;
	      dictAdd (tagDict, ccp, &tag) ;
	    }
       }
     else
       dictAdd (tagDict, tsf->tag, &tag) ;

     /* parse the column/sample */
     aceInStep (ai, '\t') ;
     ccp = aceInWordCut (ai, "\t", &cc) ;
     if (! ccp || *ccp == '#')
       continue ;
      if (! strncmp (ccp, "_P_", 3)) /* do not parse percensamplees */
	continue ; 
      if (! tsf->sample)
	{
	  if (tsf->sampleOldNameDict)
	    {
	      int old, new ;
	      if (dictFind (tsf->sampleOldNameDict, ccp, &old))
		{
		  new = keySet (tsf->sampleRenameKs, old) ;
		  ccp = dictName (tsf->sampleNewNameDict, new) ;
		}
	    }
	  if (tsf->sampleSelectOnly && tsf->sampleSelectDict 
	      && ! dictFind (tsf->sampleSelectDict, ccp, 0)
	      )
	    continue ;
	  sample = sampleAdd (tsf, ccp) ;
	  if (! sample)  /* not selected */
	    continue ;
	}
      else
	sample = sampleAdd (tsf, tsf->sample) ;
     aceInStep (ai, '\t') ;


     /* parse the cell-format */
     aceInStep (ai, '\t') ;
     ccp = aceInWord (ai) ;
     nFields = 0 ;
     if (ccp)
       {
	 char *cq = typeBuf ;
	 const  char *cp = ccp - 1 ;
	 int kk = MAXSAMPLEDEPTH ;
	 
	 while (*++cp)
	   {
	     int k = 0 ;
	     char cc ;

	     while (*cp >= '0' && *cp <= '9')
	       {  k = 10*k + (*cp - '0') ; cp++ ; }
	     if (k == 0) k = 1 ;
	     kk -= k ;
	     nFields += k ;
	     if (kk <=0)
	       messcrash ("Overflow, number of fields = %d > %d at line %d of file %s.\nWe do not expect more than %d value following a sample, please edit the code or reformat  the data"
			  , nFields
			  , MAXSAMPLEDEPTH
			  , aceInStreamLine (ai)
			  , aceInFileName (ai)
			  , MAXSAMPLEDEPTH
			  ) ;
	       
	     cc = *cp ;
	     if (! cc)
	       {
		 cc = 'i' ; /* defaul;ts to Integer */
		 cp-- ;     /* avoid looping */ 
	       }
	     while (k--)
	       *cq++ = cc ;
	   }
	 *cq = 0 ;
       }
     else
       continue ;
     hit = bigArrayp (hits, iMax++, HIT) ;
     hit->sample = sample ;
     hit->tag = tag ;
     strncpy (hit->types, typeBuf, MAXSAMPLEDEPTH) ;
     
     for (i = 0 ; i < nFields ; i++)
       {
	 char cc ;
	 if (! aceInStep (ai, '\t'))
	   break ;
	  hit->n = i + 1 ;
	  switch (hit->types[i])
	    {
	    case 'i':
	    case 'c':
	      ccp = aceInWord (ai) ;
	      if (ccp && sscanf (ccp, "%ld%c", &x, &cc) == 1)
		{
		  hit->x[i] = x ;
		  hit->flag |= (1 << i) ;
		}
	      else
		messcrash ("Non integer number %s in column %d line %d file %s\n"
			   , ccp ?ccp : "NULL"
			   , 4 + i
			   , aceInStreamLine (ai)
			   , aceInFileName (ai)
			   ) ;
	      break ;
	    case 'f':
	      if (aceInDouble (ai, &z))
		{
		  hit->z[i] = z ;
		  hit->flag |= (1 << i) ;
		}
	      else
		messcrash ("Non float number in column %d line %d file %s\n"
			   , 4 + i
			   , aceInStreamLine (ai)
			   , aceInFileName (ai)
			   ) ;
	      break ;
	    default :
	      ccp = aceInWord (ai) ;
	      if (ccp)
		{
		  int k = 0 ;
		  dictAdd (valDict, ccp, &k) ; 
		  hit->x[i] = k ; 
		  hit->flag |= (1 << i) ;
		}
	      else
		messcrash ("Non word in column %d line %d file %s\n"
			   , 4 + i
			   , aceInStreamLine (ai)
			   , aceInFileName (ai)
			   ) ;
	      break ;
	    }
       }
     if (i > keySet (sampleDepth, sample))
       keySet (sampleDepth, sample) = i ;
     if (! tsf->noMerge && iMax > 1 && tsfMergeOne (hit -1 , hit, action))
       {
	 iMax-- ;
	 hit-- ;
       }
    }
  ac_free (h) ;

  return bigArrayMax (hits) ;
} /*  tsfParseTsf */

/*************************************************************************************/

static void tsfParseMetaData (TSF *tsf, ACEIN ai)
{
  char cutter, *ccp ;

  while (aceInCard (ai))
    {
      int type = 0 ;
      ccp = aceInPos (ai) ;
      if (ccp[0] != '#')
	{
	  aceInCardBack (ai) ;
	  return ;
	}
      if (! ccp)
	continue ;

      if (!strncmp (ccp, "####", 4))    /* continue */
	continue ;

      if (!strncmp (ccp, "###", 3))    /* get sticky title */
	{
	  ccp = aceInWordCut (ai, "\t", &cutter) ;
	  if (tsf->title)              /* not if already set */
	    continue ;
	  ccp += 3 ;
	  while (*ccp == ' ')          /* gobble spaces */
	    ccp++ ;
	  if (! ccp[0])
	    continue ;
	  tsf->title = strnew (ccp, tsf->h) ;
	  continue ;
	}

      if (!strncmp (ccp, "##", 2))    /* get captions */
	{
	  ccp += 2 ;
	  while (*ccp == ' ')          /* gobble spaces */
	    ccp++ ;
	  if (! ccp[0])
	    continue ;
	  if (! strncmp (ccp, "CLEAN", 5))
	    tsf->tableCaptionDict = dictHandleCreate (256, tsf->h) ;
	  else
	    dictAdd (tsf->tableCaptionDict, ccp, 0) ; /* gobble repeats */ 
	  continue ;
	}

      type = 0 ;
      if (tsf->isTableIn && !strncmp (ccp, "#@", 2))  type = 1 ;    /* get cell captions */
      else if (! tsf->isTableIn && !strncmp (ccp, "#@", 2))  type = 2 ;  /* get sample ordering */
      else if (tsf->isTableIn && !strncmp (ccp, "#", 1))  type = 2 ;  /* get sample ordering */
      else if (! tsf->isTableIn && !strncmp (ccp, "#", 1)) type = 1 ;    /* get cell captions */
     
      if (type == 2) /* establish sample names ordering */
	  {
	  int col = 2 ;
	  

	  ccp = aceInWordCut (ai, "\t", &cutter) ;
	  aceInStep (ai, '\t') ;
	  ccp = aceInWord (ai)  ;
	  if (tsf->isTableIn) /* jump Line */
	    {
	      aceInStep (ai, '\t') ;
	      ccp = aceInWord (ai)  ;
	    }
	  if (ccp)
	    tsf->corner_title = strnew (ccp, tsf->h) ;
	  if  (! tsf->sample)
	    while (1)
	      {                          /* grab a sticky order for the samples */
		aceInStep (ai, '\t') ;
		ccp = aceInWord (ai)  ;
		if (ccp)            
		  {
		    int sample = sampleAdd (tsf, ccp) ; /* get sorted samples */
		    if (! sample)  /* not selected */
		      continue ;
		    keySet (tsf->col2sample, col++) = sample ;
		  }
		else
		  break ;
	      }
	  else
	    {
	      int sample = sampleAdd (tsf, tsf->sample) ;
	      keySet (tsf->col2sample, col++) = sample ;
	    }
	  keySetMax (tsf->col2sample) = col ;   /* max for this file, if it is a table since all columns must have a sample */
	  continue ;
	}

      if (type == 1 )    /* get cell-captions */
	{
	  ccp += 2 ;
	  while (*ccp == ' ')          /* gobble spaces */
	    ccp++ ;
	  if (! ccp[0])
	    continue ;
	  if (! strncmp (ccp, "CLEAN", 5))
	    tsf->cellCaptionDict = dictHandleCreate (256, tsf->h) ;
	  else
	    dictAdd (tsf->cellCaptionDict, ccp, 0) ; /* gobble repeats */ 
	  continue ;
	}

    }

} /* tsfParseMetaData */

/*************************************************************************************/

static void tsfParse (TSF *tsf)
{
  AC_HANDLE h = ac_new_handle () ;
  
  if (tsf->inFileOfFileList)
    {
      vTXT txt = 0 ;
      char *sep = "" ;
      ACEIN ai = aceInCreate (tsf->inFileOfFileList, 0, h) ;
      if (!ai)
	messcrash ("Cannot open -file_of_file_list %s\n", tsf->inFileOfFileList) ;

      txt = vtxtHandleCreate (h) ;
      while (aceInCard (ai))
	{
	  char cutter, *cp ;
	  while ((cp = aceInWordCut (ai, " ,\t", &cutter)))
	    {
	      if (cp && *cp && *cp != '#')
		{
		  vtxtPrintf (txt, "%s%s",sep, cp) ;
		  sep = "," ;
		}
	    }
	}
      if (*sep)
	tsf->inFileList = vtxtPtr (txt) ;
    }

  if (tsf->inFileList)
    {
      ACEIN ai = 0 ;
      int nn = 0 ;
      const char *error = 0 ;

      /* check that all files exist */
      if (! aceInCheckList (tsf->inFileList, &error, h))
	messcrash ("Bad parameter -i, %s"
		       ,  error
		       ) ;	
 
      /* parse a validated list of files */
      while ((ai = aceInCreateFromList (ai, &nn, tsf->inFileList, tsf->gzi, h)))
	{
	  aceInSpecial (ai, "\n") ;
	  tsfParseMetaData (tsf, ai) ;
	  
	  if (tsf->isTableIn)
	    tsfParseTable (tsf, ai) ;
	  else
	    tsfParseTsf (tsf, ai) ;
	}
     }
  else
    {
      ACEIN ai = aceInCreate (0, tsf->gzi, h) ;
      if (tsf->isTableIn)
	tsfParseTable (tsf, ai) ;
      else
	tsfParseTsf (tsf, ai) ;
    }

  ac_free (h) ;

  return  ;
} /* tsfParse */

/*************************************************************************************/

#ifdef JUNK

static KEYSET tsfMyLines = 0 ;

static int tsfMyOrder (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va ;
  const HIT *vp = (const HIT *)vb ;
  int n ;

  n = keySet (tsfMyLines, up->sample -1) - keySet (tsfMyLines, vp->sample -1) ; if (n) return n ;
  n = keySet (tsfMyCols, up->tag -1) - keySet (tsfMyCols, vp->tag -1) ; if (n) return n ;

  return 0 ;
} /* tsfSampleOrder */

/*************************************************************************************/

static void tsfChronoOrder (TSF *tsf)
{
  AC_HANDLE h = ac_new_handle () ;
  HIT *hit ;
  int ii, iiMax = bigArrayMax (tsf->hits) ;
  int i, tMax = dictMax (tsf->sampleDict) ;
  int j, sMax = dictMax (tsf->tagDict) ;
  KEYSET cols = keySetHandleCreate (h) ;
  KEYSET lines = keySetHandleCreate (h) ;
  Array bb = 0 ;
  BOOL debug = FALSE ;

  if (0) 
    tsfChronoOrder (tsf) ; /* to please the compiler */

  for (i = 0 ; i < tMax ; i++)
    keySet (lines, i) = i ;
  for (j = 0 ; j < sMax ; j++)
    keySet (cols, j) = j ;
  
  if (tsf->chronoOrder)
    {
      bb = arrayHandleCreate (tMax * sMax, double, h) ; 
      array (bb, tMax * sMax - 1, double) = 0 ;
      for (ii = 0, hit = bigArrp (tsf->hits, 0, HIT) ; ii < iiMax ; hit++, ii++)
	{
	  double z  ;
	  
	  i = hit->sample - 1 ;
	  j = hit->tag - 1 ;
	  z = hit->z[0] + hit->x[0] ;
	  array (bb, i * sMax + j, double) = z ;
	}
      topoChronoOrder (bb, lines, cols) ;
    }

  tsfMyLines = arrayHandleCopy (lines, 0) ;
  tsfMyCols = arrayHandleCopy (cols, 0) ;

  for (i = 0 ; i < tMax ; i++)
    keySet (tsfMyLines, keySet (lines, i)) = i ;
  for (j = 0 ; j < sMax ; j++)
    keySet (tsfMyCols, keySet (cols, j)) = j ;

  if (debug)
    for (ii = 0, hit = bigArrp (tsf->hits, 0, HIT) ; ii < iiMax ; hit++, ii++)
      {
	double z  ;
	
	i = hit->sample ;
	j = hit->tag ;
	z = hit->z[0] + hit->x[0] ;
	printf("%d\t%d\t%.1g\n",  i, j, z) ; 
      }

  if (tsf->chronoOrder)
    bigArraySort (tsf->hits, tsfMyOrder) ;

  if (debug)
    for (ii = 0, hit = bigArrp (tsf->hits, 0, HIT) ; ii < iiMax ; hit++, ii++)
      {
	double z  ;
	
	i = hit->sample ;
	j = hit->tag ;
	z = hit->z[0] + hit->x[0] ;
	printf("%d\t%d\t%.1g\n",  i, j, z) ; 
      }



  ac_free (h) ;
  return ;
}
#endif

/*************************************************************************************/
/*************************************************************************************/

static BOOL tsfMergeOne (HIT *hit, HIT *hit2, int action)
{
  int n, n1, n2, i ;
  BOOL ok = FALSE ;

  if ( !hit2->sample)
    return ok ;
  if ( hit->sample != hit2->sample)
    return ok ;
  if ( hit->tag != hit2->tag)
    return ok ;
  n1 = hit->n ;
  n2 = hit2->n ;
  n = n1 < n2 ? n1 : n2 ;
  if (strncmp (hit->types, hit2->types, n))
    return ok ;
  
  hit2->sample = 0 ;
  ok = TRUE ;
  n = n1 > n2 ? n1 : n2 ;
  if (n1 < n2) 
    memcpy (hit->types, hit2->types, n2) ;
  hit->flag |= hit2->flag ;
  
  hit->n = n ;
  for (i = 0 ; i < n ; i++)
    if (i < n2)
      switch ((int)hit->types[i])
	{
	case 'i':   /* add */
	  switch (action)
	    {
	    case 0: 
	    default: /* add */
	      hit->x[i] += hit2->x[i] ;
	      break ;
	    case 1:  /* min */
	      if (i > n1 || hit->x[i] > hit2->x[i])
		hit->x[i] = hit2->x[i] ;
	      break ;
	    case 2:  /* max */
	      if (i > n1 || hit->x[i] < hit2->x[i])
		hit->x[i] = hit2->x[i] ;
	      break ;
	    }
	  break ;
	case 'f':   /* add */
	  switch (action)
	    {
	    case 0: 
	    default: /* add */
	      hit->z[i] += hit2->z[i] ;
	      break ;
	    case 1:  /* min */
	      if (i > n1 || hit->z[i] > hit2->z[i])
		hit->z[i] = hit2->z[i] ;
	      break ;
	    case 2:  /* max */
	      if (i > n1 || hit->z[i] < hit2->z[i])
		hit->z[i] = hit2->z[i] ;
	      break ;
	    }
	  break ;
	case 'c':
	case 't':   /* replace */
	default:
	  switch ((int)hit->types[i])
	    {
	    case 'i':   /* add */
	      switch (action)
		{
		case 0: 
		default: /* replace */
		  hit->x[i]  = hit2->x[i] ;
		  break ;
		case 1:  /* min */
		  if (i > n1 || hit->x[i] > hit2->x[i])
		    hit->x[i] = hit2->x[i] ;
		  break ;
		case 2:  /* max */
		  if (i > n1 || hit->x[i] < hit2->x[i])
		    hit->x[i] = hit2->x[i] ;
		  break ;
		}
	      break ;
	      break ;
	    }
	}

  return ok ;
}  /* tsfMergeOne */

/*************************************************************************************/

static long int tsfMerge (TSF *tsf)
{
  BigArray hits = tsf->hits ;
  long int ii, jj = 0, iMax = bigArrayMax (hits) ;
  HIT *hit, *hit2 ;
  int action = 0 ;

  /* sort and cumulate the values */
  if (iMax)
    {
      bigArraySort (hits, tsfSampleOrder) ;
      for (ii = 0, hit = bigArrp (hits, 0, HIT) ; ii < iMax ; ii++, hit++)
	{
	  if ( !hit->sample)
	    continue ;
	  
	  for (jj = ii + 1, hit2 = hit+1 ; jj < iMax && (! hit2->sample || hit2->sample == hit->sample) && hit2->tag == hit->tag ; jj++, hit2++)
	    tsfMergeOne (hit, hit2, action) ;
	}
      /* clean up */
      for (ii = jj = 0, hit2 = hit = bigArrp (hits, ii, HIT) ; ii < iMax ; ii++, hit++)
	{
	  if ( !hit->sample)
	    continue ;
	  if (hit != hit2)
	    *hit2 = *hit ;
	  jj++ ; hit2++ ; 
	}
    }
  bigArrayMax (hits) = jj ;
  return bigArrayMax (hits) ;
} /* tsfMerge */

/*************************************************************************************/
/*************************************************************************************/

static void tsfExportCaption (TSF *tsf, ACEOUT ao)
{
  int i, iMax ;
  int pass = 0 ;

  aceOutDate (ao, "###"
	      ,  tsf->title ? tsf->title : "No title"
	      ) ;
  
  iMax = dictMax (tsf->tableCaptionDict) ;
  for (i = 1 ; i <= iMax ; i++)
    aceOutf (ao, "## %s\n", dictName (tsf->tableCaptionDict, i)) ;

  for (pass = 0 ; pass < 2 ; pass++)
    {
      if (
	  (pass == 0 && tsf->isTableOut) ||
	  (pass == 1 && ! tsf->isTableOut)
	  )
	{
	  iMax = dictMax (tsf->cellCaptionDict) ;
	  for (i = 1 ; i <= iMax ; i++)
	    aceOutf (ao, "%s %s\n"
		     , tsf->isTableOut ? "#@" : "#"
		     , dictName (tsf->cellCaptionDict, i)) ;
	}
      if (
	  (pass == 1 && tsf->isTableOut) ||
	  (pass == 0 && ! tsf->isTableOut)
	  )
	{
	  iMax = dictMax (tsf->sampleDict) ;
	  if (iMax)
	    {
	      const char *ccp = tsf->corner_title ;
	      if (tsf->isTableOut)
		aceOut (ao, "# Line") ;
	      else
		aceOut (ao, "#@") ;
	      aceOutf (ao, "\t%s", ccp ? ccp : "") ;
	      for (i = 1 ; i <= iMax ; i++)
		aceOutf (ao, "\t%s", dictName (tsf->sampleDict, i)) ;
	      aceOut (ao, "\n") ;
	    }
	}
    }

  return ;
} /* tsfExportCaption */

/*************************************************************************************/

static void tsfExportTsf (TSF *tsf, ACEOUT ao)
{
  long int ii ;
  const long iMax = bigArrayMax (tsf->hits) ;
  const HIT *hit = iMax ? bigArrp (tsf->hits, 0, HIT) : 0 ;
  const DICT * valDict = tsf->valDict ;
  const DICT * sampleDict = tsf->sampleDict ;
  const DICT * tagDict = tsf->tagDict ;
  int k, hasP = 0 ;
  int oldTag = 0, oldSample = 0 ;

  for (ii = 0 ; ii < iMax ; ii++)
    {
      int i, n ;
      HIT hh = hit[ii] ;
      const char *sampleName = dictName (sampleDict, hh.sample) ;

      if (! hh.tag && ! hh.sample)
	continue ;

      if (ii > 0 && tsf->skip1 > 0 && hh.tag != oldTag) 
	{ 
	  for (i = 0 ; i < tsf->skip1 ; i++)
	    aceOut (ao, "\n") ;
	}
      else 
	if (ii > 0 && tsf->skip2 > 0 && hh.sample != oldSample)
	  { 
	    for (i = 0 ; i < tsf->skip2 ; i++)
	      aceOut (ao, "SS\n") ;
	  }

      oldSample = hh.sample ;
      oldTag = hh.tag ;

      aceOutf (ao, "%s", dictName (tagDict, hh.tag)) ;
      aceOutf (ao, "\t%s", sampleName) ;
      n = hh.n ;
      aceOutf (ao, "\t%s", hh.types) ; 
      for (i = 0 ; i < n ; i++)
	if (hh.flag & (1 << i))
	  switch (hh.types[i])
	    {
	    case 0:
	      break ;
	    case 'i':
	    case 'c':
	      aceOutf (ao, "\t%ld", hh.x[i]) ;
	      break ; 
	    case 'f':
	      aceOutf (ao, "\t%g", hh.z[i]) ;
	      break ; 
	    default:
	      k = hh.x[i] ;
	      aceOutf (ao, "\t%s", dictName (valDict, k)) ;
	      break ;
	    }
	else
	  break ;
      aceOutf (ao, "\n") ;

      if (strncmp (sampleName, "P_", 2))
	hasP = ii + 1 ;
    }
  
  /* export the percensamplees */
  for (ii = 0 ; ii < hasP ; ii++)
    {
      int i, n ;
      double z ;
      HIT hh = hit[ii] ;
      const char *sampleName = dictName (sampleDict, hh.sample) ;

      if (strncmp (sampleName, "P_", 2))
	continue ;

      /* Export the percensamplees */
      aceOutf (ao, "%s", dictName (tagDict, hh.tag)) ;
      aceOutf (ao, "\t%s", sampleName) ;
      n = hh.n ;
      aceOutf (ao, "\t%s\t100%%", hh.types) ; 
      z = hh.x[0] ;
      if (z == 0) z = 100 ;
      z = 100.0/z ;
      for (i = 1 ; i < n ; i++)
	if (hh.flag & (1 << i))
	  switch (hh.types[i])
	    {
	    case 'i':
	      aceOut (ao, "\t") ;
	      aceOutPercent (ao, z * hh.x[i]) ;
	      aceOut (ao, "%") ;
	      break ; 
	    case 'f':
	      aceOut (ao, "\t") ;
	      aceOutPercent (ao, z * hh.z[i]) ;
	      aceOut (ao, "%") ;
	    break ; 
	    default:
	      k = hh.x[i] ;
	      aceOutf (ao, "\t%s", dictName (valDict, k)) ;
	      break ;
	    }
	else
	  break ;
      aceOutf (ao, "\n") ;	
    }

  return ;
} /* tsfExportTsf */

/*************************************************************************************/
/* check if the kk first fields are identical with previous line */
#ifdef JUNK
static int tsfSame (int same, int kk, vTXT txt, vTXT oldTxt)
{
  if (same + kk < 0 ) /* we need to study the line */
    {
      const char *cp = vtxtPtr (txt) ;
      const char *cq = vtxtPtr (oldTxt) ;
      
      if (cp && cq)
	{
	  int n = strlen (cp) ;
	  if (strncmp (cp, cq, n))
	    return kk ;
	}
    }
  return same ;
} /* tsfSame */
#endif
/*************************************************************************************/

static void tsfExportTable (TSF *tsf, ACEOUT ao)
{
  AC_HANDLE h = ac_new_handle () ;
  long int ii, jj ;
  long jjMax, iiMax = bigArrayMax (tsf->hits) ;
  const DICT * valDict = tsf->valDict ;
  const DICT * sampleDict = tsf->sampleDict ;
  const DICT * tagDict = tsf->tagDict ;
  const char *NA = tsf->NA ;
  int k ;
  int sampleMax = dictMax (sampleDict) ;
  
  for (ii = 0 ; ii < iiMax ; ii++)
    {
      HIT *vp, *up = up = bigArrp (tsf->hits, ii, HIT) ;
      int tag = up->tag ;
      int n, sample ;
	
      aceOutf (ao, "%d", ii+1) ;  /* line number */
      if (tsf->noMerge)
	{
	  char *cr, *ccp = (char *) dictName (tagDict, tag) ; /* hack unprotect the dictName */
	  cr = strstr (ccp, "#_#_#") ;
	  *cr = 0 ;
	  aceOutf (ao, "\t%s", ccp) ;
	  *cr = '#' ;
	}
      else
	aceOutf (ao, "\t%s", dictName (tagDict, tag)) ;

      /* locate all lines corresponding to this tag */
      for (jj = ii, vp = up ; jj < iiMax && vp->tag == tag ; vp++, jj++)
	;
      jjMax = jj ;
      jj = ii ; vp = up ;

      for (sample = 1 ; sample <= sampleMax ; sample++)
	{
	  char *sep = "" ;

	  /* locate the relevant sample, it may be absent */
	  /* it is probably the next sample */
	  if (0)
	    for ( ; jj < jjMax ; vp++, jj++)
	      if (vp->sample == sample)
		break ;

	  if (jj >= jjMax || vp->tag != tag || vp->sample != sample) /* tough luck, rescan the group */
	    for (jj = ii, vp = up ; jj < jjMax ; vp++, jj++)
	      if (vp->sample == sample)
		break ;

	  aceOut (ao, "\t") ;    /* always create a cell */
	  if (vp->sample != sample || vp->tag != tag || ( (vp->flag & 1) == 0))
	    {
	      if (NA) 
		aceOut (ao, NA) ;
	      continue ;
	    }

	  /* export the cell */
	  for (n = 0 ; n < vp->n ; n++)   /* subcell counter */
	    {
	      if (vp->flag & (1 << n))  /* available data */
		{
		  aceOut (ao, sep) ;
		  sep = "," ;
		  switch (vp->types[n])
		    {
		    case 'i':
		    case 'c':
		    case 'I':
		    case 'C':
		      aceOutf (ao, "%ld", vp->x[n]) ;
		      break ; 
		    case 'f':
		      aceOutf (ao, "%g", vp->z[n]) ;
		      break ; 
		    default:
		      k = vp->x[n] ;
		      aceOutf (ao, "%s", dictName (valDict, k)) ;
		      break ;
		    }
		  /* export only once */
		  vp->flag ^= (1 << n) ;
		}
	    }
	}
      ii = jjMax - 1   ; 
      aceOut (ao, "\n") ;
    }
# ifdef JUNK
  KEYSET skip = keySetHandleCreate (h) ;
  
  for (i = 0 ; i < tsf->skip1 ; i++)
    keySet (skip, i)++ ;
  for (i = 0 ; i < tsf->skip2 ; i++)
    keySet (skip, i) += 2 ;
  for (i = 0 ; i < tsf->skip3 ; i++)
    keySet (skip, i) += 3 ;


	sampleName = dictName (sampleDict, sample) ;
	if (pass == 0 && ! strncmp (sampleName, "P_", 2))
	  continue ;
	if (pass == 1 && strncmp (sampleName, "P_", 2))
	  continue ;

  /* to be used for pass == 1 */
  /* export the sample percensamplees */
  for (i = 1 ; i < hasP ; i++)
    {  
      int jMax = keySet (tsf->sampleDepth, i) ;
      const char *sampleName = dictName (sampleDict, i) ;
      
      if (strncmp (sampleName, "P_", 2))
	continue ;
      
      for (j = 0 ; j < jMax ; j++)
	{
	  kMax++ ;
	  aceOutf (ao, "\t%%%s", sampleName + 2) ;
	  if (jMax > 1)
	    aceOutf (ao, ":%d", j+1) ;
	}
    }

  // la suite je sais pas a quoi cela sert (2020-08-07 */

      
      /* export the percensamplees */
      for (i = 1 ; i < hasP ; i++)
	{  
	  double z ;
	  int jMax = keySet (tsf->sampleDepth, i) ;
	  const char *sampleName = dictName (sampleDict, i) ;
	  
	  if (strncmp (sampleName, "P_", 2))
	    continue ;
	  
	  kk++ ; 
	  vtxtPrintf (txt, "\t100%%") ; 
	  z = (hit->types[0] == 'i' ? hit->x[0] : hit->z[0]) ;
	  if (z == 0) z = 100 ;
	  z = 100.0/z ;
	  same = tsfSame (same, kk, txt, oldTxt) ;
	  
	  for (j = 1 ; j < jMax ; j++)
	    {
	      double zj = (hit->types[i] == 'i' ? hit->x[i] : hit->z[i]) ;
	      kk++ ; 
	      vtxtPrintf (txt, "\t") ;
	      vtxtPercent (txt, z * zj) ;
	      vtxtPrintf (txt, "%%") ;
	      same = tsfSame (same, kk, txt, oldTxt) ;
	    }
	}
      
      if (vtxtPtr (txt))
	aceOut (ao, vtxtPtr (txt)) ;
    }
#endif

  aceOut (ao, "\n") ;
  ac_free (h) ;

  return ;
} /* tsfExportTable */

/*************************************************************************************/

static void tsfExportDo (TSF *tsf, ACEOUT ao)
{
  if (tsf->isTableOut)
    {
      myDict = tsf->tagDict ;
      bigArraySort (tsf->hits, tsfTagOrder) ; /* always needed, to sort the samples */
      tsfExportTable (tsf, ao) ;
    }
  else
    {
      if (! tsf->sortCommand)
	{
	  myDict = tsf->tagDict ;
	  bigArraySort (tsf->hits, tsfTagOrder) ;
	}
      tsfExportTsf (tsf, ao) ;
    }
} /* tsfExportDo */

/*************************************************************************************/

static void tsfTranspose (TSF *tsf)
{
  DICT *dict ;
  HIT *hit = bigArrp (tsf->hits, 0, HIT) ; ;
  long int i, iMax = bigArrayMax (tsf->hits) ;

  dict = tsf->tagDict ; tsf->tagDict = tsf->sampleDict ; tsf->sampleDict = dict ;

  for (i = 0 ; i < iMax ; hit++, i++)
    {
      int x = hit->tag ; hit->tag = hit->sample ; hit->sample = x ;
    }
} /* tsfTranspose */

/*************************************************************************************/

static void tsfExport (TSF *tsf)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = 0 ;
  
  char *fName = tsf->outFileName 
    ? hprintf (h, "%s%s"
	       , tsf->outFileName
	       , tsf->isTableOut ? ".txt" : ".tsf"
	       ) 
    : 0 
    ;
  
  ao = aceOutCreate (fName, 0, tsf->gzo, h) ;

  if (tsf->transpose)
    tsfTranspose (tsf) ;

  tsfExportCaption (tsf, ao) ;
  
  if (tsf->sortCommand)
    {
      const char* command ;

      if (tsf->outFileName)
	{
	  ac_free (ao) ;
	  command = hprintf (h, "sort %s %s  >>  %s%s"
			     , tsf->sortCommand 
			     , tsf->gzo ? " | gzip " : ""
			     , fName
			     , tsf->gzo ? ".gz" : ""
			     ) ;
	  ao = aceOutCreateToPipe (command, h) ;
	}
      else if (! tsf->gzo)
	{
	  command = hprintf (h, "sort %s "
			     , tsf->sortCommand 
			     ) ;
	  ao = aceOutCreateToPipe (command, h) ;
	}  
      else
	{
	  command = hprintf (h, "sort %s %s "
			     , tsf->sortCommand 
			     , tsf->gzo ? " | gzip " : ""
			     ) ;
	  ao = aceOutCreateToPipe (command, h) ;
	}  

    }
  tsfExportDo (tsf, ao) ;
      

  ac_free (h) ;
} /* tsfExport */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/
/* tag -tag NA12878MOD -laneList tmp/TSNP/NA12878MOD/LaneList -t 20 -target_fasta TARGET/CHROMS/hs.chrom_20.fasta.gz -t1 25000001 -t2 30010000 -minSnpFrequency 18 -minSnpCover 10 -minSnpCount 4 -target_class Z_genome -o tata -maxLanes 4
 */

static void usage (FILE *out, const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (out,
	   "// tsf: A generic tool to manipulate numeric data and tables\n"
	   "//   Contact: Jean  Thierry-Mieg, NCBI  mieg@ncbi.nlm.nih.gov\n"
	   "//   This code is public and distributed without restrictions or guarantees\n"
	   "// Usage: tsf [parameters] [-i <file_list>] [-o <output>]\n"
	   "//      try: -h --help -v --version\n"
	   "// Examples:\n"
	   "//     tsf -i f1.tsf,f2.tsf,f3.tsf\n"
	   "//       Combine three .tsf files, into a single one while\n"
	   "//       adding up the numerical values with corresponding captions\n"
	   "//       Many other operations are possible (see below).\n"
	   "//     tsf -i f1.tsf,f2.tsf --setSample newSampleName\n"
	   "//       For each tag, combine the counts of all samples\n"
	   "//     tsf -i f1.tsf -O tabular -o my_table\n"
	   "//       Import a tsf file, export a table called my_table.txt\n"
	   "//     tsf -I tabular -i f1.txt,f2.txt\n"
	   "//       Merge several tables into a probably larger table\n"
	   "//       adding up the values with corresponding captions\n"
	   "//     tsf -I tabular --transpose --lineSelect a,b,c --colSelect A,B\n"
	   "//       Manipulate a generic tab delimited table\n"
	   "//     tsf -i f1.tsf,f2.tsf,f3.tsf  --no_merge\n"
	   "//       Combine and sort the files, maintain the original lines\n"
	   "//\n"
	   "// INPUT:\n"
	   "//   Default: read from stdin.\n"
	   "//   -i, --in_file_list <file_list>\n"
	   "//      A comma delimited list of input files, for example\n"
	   "//       -i f1.tsf,f2.tsf,f3.tsf\n"
	   "//   -f <file_name> : file of file names\n"
	   "//      File contains a list of files, comma or space or tab or line separated\n"
	   "//      ATTENTION all file names must be local or fully qualified starting at /\n"
	   "//   --gzi\n"
	   "//      gunzip the input\n"
	   "//      This is automatic for all input files called .gz\n"
	   "//      The option is useful when reading stdin, or a file not named .gz\n"
	   "// INPUT FORMAT:\n"
	   "//   -I, --In_format [tsf|tabular] : default tsf\n"
	   "//      The input files are provided in tsf or in tabular format (defined below)\n"
	   "//\n"
	   "// OUTPUT:\n"
	   "//   Default: export to stdout\n"
	   "//   -o, --out <file_name>\n"
  	   "//   --gzo : the output file will be gzipped\n"
	   "// OUTPUT FORMAT:\n"
	   "//   -O [tsf|table]: short form of\n"
	   "//   --Out_format [tsf|table] : default tsf\n"
	   "//     The output files are exported in tsf or in tabular format (defined below)\n"
	   "//     Examples:\n"
	   "//       tsf -o foobar                      exports foobar.tsf\n"
	   "//       tsf -o foobar -gzo                 exports foobar.tsf.gz\n"
	   "//       tsf -o foobar -gzo -O tabular      exports foobar.txt.gz\n"
	   "//       tsf -I tsf -O tabular              reformat a tsf file as a table\n"
	   "//\n"
	   "// TITLE AND CAPTIONS :   common to tsf and tabular formats\n"
	   "//    The program exports at the top of each file several lines of metadata\n"
	   "//      ### Title data file-name\n"
	   "//      ##  Captions when available\n"
	   "//      #   Ordered list or samples, which serve as column titles in tabular mode\n"
	   "//    When reading data,\n"
	   "//      The title is recovered from the first available ### line\n"
	   "//      The captions are recovered from the ## lines, dropping the duplicates\n"
	   "//      The sample list accumulates the samples foung in the # lines, droping doubles\n"
	   "//      they are reexported as is, in the same order, ignoring duplicated lines.\n"
	   "//    The title can be created or reset via\n"
	   "//   --title <title> : set a descriptive title for the table\n"
	   "//   --corner_title <corner_title> : set a short title in the left upper corner of the table\n"
	   "//     The captions can be eliminated with the option\n"
	   "//   --no_caption  : eliminate all caption lines\n"
	   "//     or cleaned up then reset by including in the input stream a line\n"
	   "//   ## CLEAN \n"
	   "//     optionally followed by the desired new ## caption lines\n"
	   "//\n"
	   "// TSF FORMAT: (default input and output format)\n"
	   "//     The tsf files are multivalued, self-descriptive and line oriented.\n"
	   "//     Each data line of the file (see examples below)  has the structure:\n"
	   "//       sample tag cell-format values\n"
	   "//     and can be thought of as a 2 dimensional table with\n"
	   "//       line and column names given in the first 2 columns\n"
	   "//       followed by a cell-format (defined below)\n"  
	   "//       followed by the data for the corresponding cell.\n"
	   "//   #  sample sample sample : in TSF format, this line is optional.\n"
	   "//      It is generated automatically and can be manually altered\n"
	   "//      to optimize the ordering of the samples in the output file\n"
	   "//      for example before generating a table (-O tabular)\n"
	   "//      This line is overridden by --sampleSelect (see below)\n"
	   "// TABULAR FORMAT: (alternative format specified as -I tabular or -O tabular)\n"
	   "//   #  sample sample sample : this line, with a single #, must appear before the data\n"
	   "//     Each sample is equivalent to the second column of the tsf format defined above\n"
	   "//     The most recent line starting with single # is treated as the current\n"
	   "//     list of samples used to interpret the data lines. This provides an easy \n"
	   "//     way to merge tables with different sets of columns, into a larger table\n"
	   "//   # FORMAT format   : optional line defining the cell-format as defined below\n"
	   "//   --format <format> : force a cell-format  on the table\n"
	   "//     If neither is provided, data manipulations or tsf exportations are forbidden,\n"
	   "//     the only allowed operation are: sampleSelect, sampleRename, transpose, sort\n"
	   "//\n"
	   "// CELL FORMAT:  \n"
	   "//     The cell-format indicates the number of data values and their types.\n"
	   "//     Types: i (integer), f (float), c (coordinate), t (text).\n"
	   "//            IFCT upper-case types allow comma separated multivalued fields.\n"
	   "//             All numbers are treated using 64bits double precision.\n"
	   "//       A number can be used to repeat a field: 5t3i == tttttiii\n"
	   "//     Example of a 3 lines tsf file:\n"
	   "//     Input:\n"
	   "//        Toyota\tCars\tii\t200\t13\n"
	   "//        Nissan\tCars\tii\t30\t70\n"
	   "//        Toyota\tCars\tii\t60\t41\n"
	   "//        Nissan\tTrucks\tii\t3\t7\n"
	   "//     Command: tsf --merge # cumulate corresponding numeric cells\n"
	   "//     Output:\n"
	   "//        Toyota\tCars\tii\t260\t54\n"
	   "//        Nissan\tCars\tii\t30\t70\n"
	   "//        Nissan\tTrucks\tii\t3\t7\n"
	   "//     Mote that the number of sample tag pairs, i.e. the size of the table, is\n"
	   "//     not limited, however the maximal number of values in a given cell is 24.\n"
	   "//\n"     
	   "// GROUPING TAGS:\n"
	   "//   -t, --setTag <tag>\n"
	   "//     If a \'tag\' is provided, it is forced on each cell, allowing\n"
	   "//     to cumulate the counts of all tags for each given sample.\n"
	   "//   Example using the same input as above:\n"
	   "//     Command: tsf --setTag All_brands\n"
	   "//     Output:\n"
	   "//        All_brands\tCars\tii\t290\t124\n"
	   "//        All_brands\tTrucks\tii\t3\t7\n"
	   "// GROUPING SAMPLES:\n"
	   "//   -s, --setSample <sample>\n"
	   "//     If a \'sample\' is provided, it is forced on each cell, allowing\n"
	   "//     to cumulate the counts of all samples for each given tag.\n"
	   "//   Example using the same input as above:\n"
	   "//     Command: tsf --setSample Vehicles\n"
	   "//     Output:\n"
	   "//        Toyota\tVehicles\tii\t260\t54\n"
	   "//        Nissan\tVehicles\tii\t33\t77\n"

	   "//\n"
	   "// SAMPLE SELECTION and RENAMING:\n"
	   "//   --sampleRename <list>: example  --sampleRename 'c1:d1,c2:d2,c3:d3'\n"
	   "//     Renames sample c1 as d1, sample c2 as d2 and so on. This applies both to\n"
	   "//     the tsf and the table format.\n"
	   "//     Sample renaming is always applied before all other operations.\n"
	   "//   --sampleSelect <list> : \n"
	   "//     example --sampleSelect 't1,t2,t3'\n"
	   "//       export only samples/columns t1,t2,t3 in that order.\n"
	   "//     example --sampleSelect 't1,t2,t3,*'\n"
	   "//       export samples/columns t1,t2,t3 then all other columns.\n"
	   "//     Otherwise, the samples/columns are sorted in order of first occurence\n"
	   "// TAG SELECTION and ORDERING\n"
	   "//   --tagSelect <list> : \n"
	   "//     example --tagSelect 's1,s2,s3'\n"
	   "//       export only lines/tags s1,s2,s3 in that order.\n"
	   "//     Otherwise, all tags are exported in alphanumeric order\n"
	   "//   --sort \"parameters\" (specific of the table format)\n"
	   "//     Calls UNIX sort, to sort the lines of a table.\n"
	   "//     example   --sort \"-k 1,1 -k 2,2 -k 5,5n\"\n"
	   "//     will sort the lines on column 1 (the tag), then 2 (the sample),\n"
	   "//     then numerically on column 5, following the convention of the unix sort.\n"
	   "//     This supersedes the order implied by --sampleSelect or --tagSelect\n"
	   "//\n"
	   "// MERGE ACTIONS\n"
	   "//     --noMerge : sort the input while maintaining the original lines\n"  
	   "//     --merge [default action] : combines corresponding data cells and subcells\n"
	   "//       For any given sample, cell-formats in different data-lines nust have compatible formats,\n"
	   "//       for example: tii and ti are compatible, but ti and tti are not.\n"
	   "//\n"
	   "//     When several lines refer to the same sample and same tag,\n"
	   "//     or to the same type and the --merge parameter forces the same tag name,\n"
	   "//     one must decide how to merge the values in the corresponding cells.\n"
	   "//     This is controlled by the parameters\n"
	   "//     --merge    no_cell_format : perform default actions on all cells and subcells\n"
	   "//       --sum      sample_list       : default action on (i,f) subcells\n"
	   "//       --replace  sample_list       : default action on t subcells\n"
	   "//       // --append   sample_list       : default action on T subcells\n"
	   "//       // --min      sample_list\n"
	   "//       // --max      sample_list\n"
	   "//     The sample_list is a comma delimited list of samples, or the word 'any', for example\n"
	   "//       --min sample1 --max sample2,sample3\n"
	   "//       the min and max will be applied on these samples, the default action on all others\n"
	   "//     For multivalued samples, the subcells are accessed using a square bracket like sample[3]\n"

	   "//     Results of the different actions:\n"
	   "//       Sum : add numerical values, default merge action for numerical values (i,f)\n"
	   "//           7  12 => 19       (missing values default to zero)\n"
	   "//       Replace: last value wins, default merge action for all other types of data (c,t,I,F,C,T)\n"
	   "//           John Paul => Paul  (missing values are ignored)\n"
	   "//       // Min: keep the smallest value (single valued types icft),\n"
	   "//           John Paul => John  (missing values are ignored)\n"
	   "//       // Max: keep the largest value (single value types icft)\n"
	   "//           John Paul => Paul  (missing values are ignored)\n"
	   "//       // Append: union of the 2 lists, ignoring duplicates (multivalued types IFCT)\n"
	   "//           4,5  4 => 4,5\n"
	   "//           4,7  5 => 4,5,7\n"
	   "//\n"
	   "// COMPUTE   (not yet programmed)\n"
	   "//     In addition, one may set missing values or compute certain columns\n"
	   "//       // --missing  any=value | sample=value,sample=value,...\n"
	   "//       // --compute  [ignore]  sample=equation,sample=equation\n"
	   "//       --file_compute  sample=equation,sample=equation\n"
	   "//     For single valued {i,f} numbers, when merging multiple files,\n"
	   "//     one can require a calculation combining the files as follows:\n"
	   "//       // --compute sample=\"sample1 + 3 * sample2 - sample3/100\"\n"
	   "//            will compute a new value for sample according to the vaules of sample,sample2,sample3\n"
	   "//       --file_compute sample=\"a + 3*b - c + 7\" \n"
	   "//     for example, if the same cell in the 3 files has the values:\n"
	   "//            4 9 5,  the exported value is:  4 + 3*9 - 5 + 7 = 33\n"
	   "//   The --file_compute is only valid when parsing a list of files and the\n"
	   "//   variables a,b,c ... in the equation refer to file 1,2,3 ...\n"
	   "//   The calculation is applied just to the given sample\n"
	   "//     Missing values must be given a default value before they are used in an equation\n"
	   "//     unless the [ignore] is specified, and any missing value in the equation returns NULL\n"
	   "//   All other samples not listed as an equation are merged\n"
	   "//\n"
	   "//   It is recommended to run oly one complex operation at a time to avoid side effects\n"
	   "//\n"
	   "// PERCENSAMPLEES\n"
	   "//   P_ and _P_ percensamplee samples\n"
	   "//     Samples called P_ have a special meaning\n"
	   "//     They follow the standard format and are cumulated as usual\n"
	   "//     But the first field is expected to be the denominator of all the other\n"
	   "//     numerical columns and whenever a sample P_ is encountered, the program also\n"
	   "//     exports a _P_ line with the corresponding percensamplees.\n"
	   "//     Please notice that on input the _P_ samples are ignored,\n"
	   "//     and recomputed from the final values in the corresponding P_ lines\n"
	   "// SEPARATORS: used only when parsing externally produced input file\n"
	   "//  this whole section must be revised, it is completely obscure\n"
	   "//  --FS <sep>: default \\t (tab character)\n"
	   "//      Field separator: separates columns in input tables\n"
	   "//  --VS <sep>: default , (comma)\n"
	   "//      Value separator: separates values in multivalued cells\n"
	   "//      Example Ford,Mercedes,Toyota  red,blue,yellow\n"
	   "//  --NA --non_available <na> [default empty]\n"
	   "//     In tabular mode, export <na> in blank cells\n"  
	   "//\n"
	   "// PRESENTATION: TRANSPOSITION and  line skipping\n"
	   "//   --transpose\n"
	   "//       Transpose a table, exchanging lines and columns, i.e. samples and tags.\n"
	   "//   --skip[12]  <int>\n"
	   "//       tsf format : skip n=1,2,3 lines when there is a change in a given column\n"
	   "//         example:  --skip1 4 --skip2 1\n"
	   "//     export 4 blank lines if the data in a column 1 vary, 1 blank line for column 2\n"
	   "//\n"
	   "// HELP\n"
	   "//    -h : short form of \n"
	   "//    --help : this online help\n"
	   ) ;

  
  if (argc > 1)
    {
      fprintf (out,
	       "//\n// You said: %s\n", commandBuf) ;
      fprintf (out,
	       "// ########## ERROR: Sorry, I do not understand the argument%s ", argc > 2 ? "s" : "") ;
      for (i = 1 ; i < argc ; i++)
	fprintf (out, "%s ", argv[i]) ;
      fprintf (out, "// Please try : tsf -h or tsf --help\n") ;
      fprintf (out, "\n") ;
    }
  exit (1) ;
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  TSF tsf ;
  AC_HANDLE h = 0 ;
  char commandBuf [4000] ;
  const char *ccp ;

  freeinit () ; 
  messErrorInit (argv[0]) ;

  h = ac_new_handle () ;
  memset (&tsf, 0, sizeof (TSF)) ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
  tsf.h = h ;
 
  if (argc < 2)
    usage (stderr, commandBuf, argc, argv) ;
  if (getCmdLineBool (&argc, argv, "-help") ||
      getCmdLineBool (&argc, argv, "--help")||
      getCmdLineBool (&argc, argv, "-h")
      )
    usage (stdout, commandBuf, 1, argv) ;
  if (getCmdLineBool (&argc, argv, "-v") ||
      getCmdLineBool (&argc, argv, "--version")
      )
    {
      fprintf (stderr, "tsf: %s\n", VERSION) ;
      exit (0) ;
    }
  tsf.gzi = getCmdLineBool (&argc, argv, "--gzi") ;
  tsf.gzo = getCmdLineBool (&argc, argv, "--gzo") ;
  tsf.transpose = getCmdLineBool (&argc, argv, "--transpose") ;
  tsf.debug = getCmdLineBool (&argc, argv, "--debug") ;

  getCmdLineOption (&argc, argv, "-i", &(tsf.inFileList)) ;
  getCmdLineOption (&argc, argv, "--in_file_list", &(tsf.inFileList)) ;
  getCmdLineOption (&argc, argv, "-f", &(tsf.inFileOfFileList)) ;
  getCmdLineOption (&argc, argv, "-o", &(tsf.outFileName)) ;
  getCmdLineOption (&argc, argv, "--out", &(tsf.outFileName)) ;
  if (getCmdLineOption (&argc, argv, "-O", &(ccp)) || getCmdLineOption (&argc, argv, "--Out_format", &(ccp)))
    {
      if (! strcmp (ccp, "tsf")) ;
      else if (! strcmp (ccp, "tabular"))
	tsf.isTableOut = TRUE ;
      else if (! strcmp (ccp, "table"))
	tsf.isTableOut = TRUE ;
      else
	{
	  fprintf (stderr, "Bad parameter -O %s, please try tsf --help\n", ccp) ;
	  exit (1) ;
	}
    }
  if (getCmdLineOption (&argc, argv, "-I", &(ccp)) || getCmdLineOption (&argc, argv, "--In_format", &(ccp)))
    {
      if (! strcmp (ccp, "tsf")) ;
      else if (! strcmp (ccp, "tabular"))
	tsf.isTableIn = TRUE ;
      else
	{
	  fprintf (stderr, "Bad parameter -I %s, please try tsf --help\n", ccp) ;
	  exit (1) ;
	}
    }

  getCmdLineOption (&argc, argv, "--sampleRename", &(tsf.sampleRenameList)) ;
  getCmdLineOption (&argc, argv, "--sampleSelect", &(tsf.sampleSelectList)) ;
  getCmdLineOption (&argc, argv, "--tagSelect", &(tsf.tagSelectList)) ;
  getCmdLineOption (&argc, argv, "--sort", &(tsf.sortCommand)) ;

  tsf.merge = getCmdLineBool (&argc, argv, "-m") || getCmdLineBool (&argc, argv, "--merge") ;
  tsf.noMerge = getCmdLineBool (&argc, argv, "--no_merge") ;
  tsf.replace = getCmdLineBool (&argc, argv, "--replace") ;
  tsf.min = getCmdLineBool (&argc, argv, "--min") ;
  tsf.max = getCmdLineBool (&argc, argv, "--max") ;
  tsf.sum = getCmdLineBool (&argc, argv, "--sum") ;
  tsf.sumAll = getCmdLineBool (&argc, argv, "--sumAll") ;
  tsf.chronoOrder = getCmdLineBool (&argc, argv, "--chronoOrder") ;
  getCmdLineOption (&argc, argv, "--compute", &(tsf.compute)) ;

  getCmdLineOption (&argc, argv, "--setTag", &(tsf.tag)) ;
  getCmdLineOption (&argc, argv, "--setSample", &(tsf.sample)) ;
  getCmdLineOption (&argc, argv, "-t", &(tsf.tag)) ;
  getCmdLineOption (&argc, argv, "-s", &(tsf.sample)) ;
  getCmdLineOption (&argc, argv, "--title", &(tsf.title)) ;
  getCmdLineOption (&argc, argv, "--corner_title", &(tsf.corner_title)) ;

  tsf.valueSeparator  = "," ;
  tsf.fieldSeparator  = "\t" ;
  getCmdLineOption (&argc, argv, "--VS", &(tsf.valueSeparator)) ;
  getCmdLineOption (&argc, argv, "--FS", &(tsf.fieldSeparator)) ;
  getCmdLineOption (&argc, argv, "--NA", &(tsf.NA)) ;
  getCmdLineOption (&argc, argv, "--non_available", &(tsf.NA)) ;

  getCmdLineInt (&argc, argv, "--skip1", &tsf.skip1) ;
  getCmdLineInt (&argc, argv, "--skip2", &tsf.skip2) ;

  if (tsf.noMerge && tsf.transpose)
    {
      fprintf (stderr, "Sorry, arguments no_merge and transpose are incompatible\n") ;
      exit (1) ;
    }
    

  if (argc > 1)
    {
      fprintf (stderr, "Unknown parameter %s, please try tsf --help\n", argv[1]) ;
      exit (1) ;
    }
  fprintf (stderr, "//%s : Start \n"
	   , timeShowNow()
	   ) ;
  
  tsfInit (&tsf) ;
  tsfParse (&tsf) ;
  if (! tsf.noMerge)
    tsfMerge (&tsf) ;
  tsfExport (&tsf) ;

  if (tsf.sampleDict)
    { 
      int mx ;
      messAllocMaxStatus (&mx) ; 
      fprintf (stderr, "## %s done, %d files %d samples %d tag analyzed, max memory %d Mb\n"
	       , timeShowNow()
	       , tsf.nInputFiles
	       , dictMax (tsf.sampleDict), dictMax (tsf.tagDict)
	       , mx) ;
     }
  ac_free (tsf.h) ;

  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderrn and to ensure all pipes are closed*/
  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
