/*  File: sra_metadata.c
 *  Author: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg. 2002
 *-------------------------------------------------------------------
 * This file is part of the ACEVIEW PACKAGE
 *	Jean Thierry-Mieg (NCBI) mieg@ncbi.nlm.nih.gov
 *
 * Description:
 *   The set of scripts :
 *             SRX_import.tcsh,  which calls SRX_import.[1-5].awk
 *   creates and acedb database summarizing the non redundant
 *   information conatianed in SRA
 *     run-info (SRE), GEO, Biosample (SAM), Project (SRP) and Expreriments (SRX)
 *   the data is transfered in the acedb classes SRE SRP SRX Biosample and GEO
 *   The schema is in waligner/metaData/wspec.SRX 
 *
 *   The idea of the present code is to synthetize the information 
 *   into synthesis tags, in class SRE, SRX, SRP, Biosample etc called
 *       Magic_author, Magic_sample etc.
 *   These tags, plus the file info etc. will later be transferred to the
 *   RumMasterDb database to contrl the MAGIC pipeline
 *   
 *  We may also create a web interface to those data, which would be
 *  managed via sra.cgi and tgifacemblyserver
 *
 * Exported functions:
 * HISTORY:
 *   2015_10_09
 *   
 * Created: Oct 2015 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */  
/* debugging flags -> slower code  */
#define ARRAY_CHECK
#define MALLOC_CHECK

#define DEBUG 1

#include "../wh/ac.h"

typedef struct sraStruct { 
  AC_HANDLE h ;
  const char *dbNam ;
  AC_DB db ;
  AC_OBJ srr, srp, srx, geo, sample ;
  DICT *dict ;
  Array baddies, tokens ;
  int max ;
  const char *template ;
  const char *project ;
  BOOL showSample, showAuthor, dbEdit, old ;
} SRA ;

typedef struct tokenStruct { 
  BOOL masked ; int x1, x2, a1, a2 ;
} TOK ;

/*************************************************************************************/
/*************************************************************************************/
/* Un capitalize name likes UNIVERSITY OF TORONTO into University of Toronto */
static char *sraUncap (const char *t0, AC_HANDLE h)
{
  char *cp, cc, *cr, **cwp ;
  int pass = 0, n, nUp = 0, nLow = 0, nWords = 0, inWord, inShortParenthese ;
  char *t, *tt ;
  char *smallWords[] = { "and ", "of ", "for ", "to ", "the ", 0 } ;
  char *bigWords[] = { "RNA", "DNA", "MRC", "NCBI", "CNRS", "NIH", "UCSC", "USA", "BGI", "JGI", "UCLA", "FDA", 0 } ;

  cp = t = tt = strnew (t0, h) ;
  return tt ;

  cc = 0 ;
  while (1)
    {
      if (cc) { *cr = cc ; t = cr + 1 ; }
      else { if (pass) break ;}
      pass++ ;
      nUp = nLow = nWords = 0 ;
      inWord = inShortParenthese = FALSE ;
      for (cr = t ; *cr && (*cr == ' ' || isalnum (*cr)) ; cr++)
	cr++ ;
      cc = *cr ; *cr = 0 ;
      cp = t - 1 ;
      while (*++cp)
	{
	  if (isalpha (*cp))
	    {
	      if (*cp == ace_upper (*cp)) nUp++ ;
	      else nLow++ ;
	    }
	  if (*cp == ' ') nWords++ ;
	}
      
      if (nLow < 10 * nUp  && nWords > 0)
	{
	  inWord = 0 ;
	  for (cp = t ; *cp ; cp++)
	    {
	      if (cp > t)
		for (cwp = smallWords ; *cwp ; cwp++)
		  {
		    n = strlen (*cwp) ;
		    if (!strncasecmp (cp, *cwp, n))
		      {
			memcpy (cp, *cwp, n) ;
			cp += n - 1 ; 
			goto laba ;
		      }
		  }
	      if (1)
		for (cwp = bigWords ; *cwp ; cwp++)
		  {
		    n = strlen (*cwp) ;
		    if (!strncasecmp (cp, *cwp, n))
		      {
			memcpy (cp, *cwp, n) ;
			cp += n - 1 ; 
			goto laba ;
		      }
		  }
	      if (inWord && ! inShortParenthese) 
		*cp = ace_lower (*cp) ;
	    laba: 
	      inWord = isalnum(*cp) ;
	      if (*cp == '(' &&
		  strchr (cp, ')') &&
		  strchr (cp, ')') < cp + 12
		  )
		inShortParenthese = 1 ;
	      if (*cp == ')')
		inShortParenthese = 0 ;
	    }
	}
    }
  return tt ;
} /* sraUncap */

/*************************************************************************************/

static int sraPrint (vTXT txt, DICT *dict, int pt, int sep)
{
  char *seps[] = {"", " ", ", "} ;

  if (sep < 0 || sep > 2) sep = 0 ;
  if (pt) { 
    const char *ccq = dictName (dict, pt) ;
    while (*ccq == ' ') ccq++ ;
    if (ispunct(*ccq) && sep == 2) sep = 1 ;
    vtxtPrintf (txt, "%s%s", seps[sep],  ccq) ;
    ccq = vtxtAt(txt, vtxtMark (txt) - 1) ;
    sep = ispunct(*ccq) ? 1 : 2 ;
  }
  else
    sep = 0 ;

  return sep ;
} /* sraPrint  */

/*************************************************************************************/

static const char *sraRemoveBaddy (const char *ccp0, Array baddies, AC_HANDLE h)
{
  int ir ;
  char *ccp = strnew (ccp0, h) ;
  if (ccp)
    {
      if (!strncmp (ccp, "MIGS ", 5))
	return (char *)1 ;
      if (!strcasecmp (ccp, "N/A"))
	return 0 ;
      if (!strcasecmp (ccp, "To be published"))
	return 0 ;
      if (!strcasecmp (ccp, "not collected"))
	return 0 ;
      if (!strcasecmp (ccp, "public"))
	return 0 ;
      ccp = strnew (ccp, h) ;
      {
	char *cp = ccp ;
	for (cp = ccp ; *cp ; cp++)
	  if (*cp == '_') *cp = ' ' ;
      } 
      while (*ccp == ' ' || *ccp == ',' || *ccp == ';' || *ccp == '.') ccp++ ;
      for (ir = 0 ; ir < arrayMax (baddies) ; ir++)
	{
	  char *baddy = arr (baddies, ir, char*) ;
	  if (ccp && baddy && ! strncmp (ccp, baddy, strlen (baddy)))
	    {
	      char cc = ccp[strlen (baddy)] ;
	      if (cc == 0 || cc == ' ' || ispunct(cc))
		ccp += strlen (baddy) ;
	    }
	}
      for (ir = 0 ; ir < arrayMax (baddies) ; ir++)
	{
	  char *cp, *cr, *cs, *ct, *baddy = arr (baddies, ir, char*) ;
	  if (ccp && baddy && (cr = strstr (ccp, baddy)))
	    {
	      char c1 = cr > ccp ? *(cr -1) : 0 ;
	      char c2 = cr[strlen(baddy)] ;
	      if (
		  (c1 == 0 || ispunct(c1) || c1 == ' ') &&
		  (c2 == 0 || ispunct(c2) || c2 == ' ') 
		  )
		{
		  cp = strnew (ccp, h) ;
		  cs = cp + (cr - ccp) ;
		  ct = cr + strlen(baddy) ;
		  if (*ct == ' ') ct++ ;
		  while ((*cs++ = *ct++)) ;
		  ccp = cp ;
		}
	    }
	}
    }
  if (ccp && ! strncasecmp (ccp, "E-MTAB-",7))
    {
      char *cq = strchr (ccp,':') ;
      if (cq && cq - ccp < 15)
	ccp = cq+1 ;
    }
  return ccp && *ccp ? ccp : 0 ;
} /* sraRemoveBaddy */

/*************************************************************************************/
/*************************************************************************************/

static int sraSetBaddies (SRA *sra, char cc)
{
  int nBaddy = 0 ;
  Array baddies = sra->baddies ;

  if (1) /* generic baddies */
    {
      array (baddies, nBaddy++, char*) =  "Keywords:" ;
      array (baddies, nBaddy++, char*) =  "Missing" ;
      array (baddies, nBaddy++, char*) =  "Not applicable" ;
      array (baddies, nBaddy++, char*) =  "Not connected" ;
      array (baddies, nBaddy++, char*) =  "Null" ;
      if (1) array (baddies, nBaddy++, char*) =  "N/A" ;
    }

  switch (cc)  /* class specific baddies */
    {
    case 'a': /* authors */
      array (baddies, nBaddy++, char*) = "NCBI (GEO)" ;
      array (baddies, nBaddy++, char*) = "NCBI" ;
      array (baddies, nBaddy++, char*) = "GEO" ;
      break ;
    case 's': /* sample */ 
      array (baddies, nBaddy++, char*) =  "Escherichia coli" ;
      array (baddies, nBaddy++, char*) =  "E. coli" ; 

      if (1)
	{
	  array (baddies, nBaddy++, char*) =  "C_elegans" ;
	  array (baddies, nBaddy++, char*) =  "C elegans" ;
	  array (baddies, nBaddy++, char*) =  "C.elegans" ;
	  array (baddies, nBaddy++, char*) =  "C. elegans" ;
	  array (baddies, nBaddy++, char*) =  "Caenorhabditis elegans" ;
	  array (baddies, nBaddy++, char*) =  "Model organism or animal sample from" ;
	  array (baddies, nBaddy++, char*) =  "Model organism or animal sample" ;
	  array (baddies, nBaddy++, char*) =  "Model organism or animal" ;
	  array (baddies, nBaddy++, char*) =  "Model organism" ;

	  array (baddies, nBaddy++, char*) =  "of Drosophila melanogaster" ;
	  array (baddies, nBaddy++, char*) =  "Drosophila melanogaster" ;
	  array (baddies, nBaddy++, char*) =  "Droshophila melanogaster" ;
	  array (baddies, nBaddy++, char*) =  "Drosophila" ;
	  array (baddies, nBaddy++, char*) =  "D. melanogaster" ; 
	  array (baddies, nBaddy++, char*) =  "D.melanogaster" ; 
	  array (baddies, nBaddy++, char*) =  "D melanogaster" ; 
	  array (baddies, nBaddy++, char*) =  "melanogaster" ; 
	  array (baddies, nBaddy++, char*) =  "E. coli" ; 
	  array (baddies, nBaddy++, char*) =  "E coli" ; 
	  array (baddies, nBaddy++, char*) =  "E.coli" ; 
	  array (baddies, nBaddy++, char*) =  "Homo sapiens" ; 
	  array (baddies, nBaddy++, char*) =  "Human" ; 
	  array (baddies, nBaddy++, char*) =  "Mus musculus" ; 
	  array (baddies, nBaddy++, char*) =  "M. musculus" ; 
	  array (baddies, nBaddy++, char*) =  "M.musculus" ; 
	  array (baddies, nBaddy++, char*) =  "Mouse" ; 
	  array (baddies, nBaddy++, char*) =  "musculus" ; 
	  array (baddies, nBaddy++, char*) =  "Rattus norvegicus" ;
	  array (baddies, nBaddy++, char*) =  "R. norvegicus" ;
	  array (baddies, nBaddy++, char*) =  "R.norvegicus" ;
	  array (baddies, nBaddy++, char*) =  "norvegicus" ;
	  array (baddies, nBaddy++, char*) =  "Rattus" ;
	  array (baddies, nBaddy++, char*) =  "E-MTAB-2544" ;
	}

      array (baddies, nBaddy++, char*) = "str. K-12 substr." ;
      array (baddies, nBaddy++, char*) = "str. K-12" ;
      array (baddies, nBaddy++, char*) = "strain K-12;" ;
      array (baddies, nBaddy++, char*) =  "K12" ;  
      array (baddies, nBaddy++, char*) =  "K-12" ;
      array (baddies, nBaddy++, char*) =  "MG1655" ;
      array (baddies, nBaddy++, char*) =  "sample from" ;
 
      break ;
    default:
      break ;
    }    
 
  arrayMax (baddies) = nBaddy ;
  return nBaddy ;
} /* sraSetBaddies */

/*************************************************************************************/
/* elements of this tag (i.e. species in samples, will be 
 * masked in the following data gathering
 * since we do not want them repeated elsewhere
 *
 * technically they are simply added to the baddy list
 */

static int sraMaskTagContent (SRA *sra, AC_OBJ obj, const char *tag, AC_HANDLE h)
{
  int nNew = 0 ;
  if (obj && ac_has_tag (obj, tag))
    {
      Array baddies = sra->baddies ;
      int nBaddy =  arrayMax (sra->baddies) ;
      const char *ccp = ac_tag_printable (obj, tag, 0) ;

      if (ccp)
	{
	  char *cp = strnew (ccp, h) ;  /* allocate on parent handle, since it is returned */
	  char *cq = strstr (cp, "(taxid:") ;
	  if (cq)
	    {
	      char *cr = strchr (cq, ')') ;
	      *cq = 0 ; if (cr) *cr = 0 ;
	      if (cq[6])
		array (baddies, nBaddy + nNew++, char*) = cq + 6 ;
	    }
	  array (baddies, nBaddy + nNew++, char*) = cp ;
	}      
    }
  return nNew ;
} /* sraMaskTagContent */

/*************************************************************************************/

static BOOL sraTokenRegister (SRA *sra, char *text)
{
  DICT *dict = sra->dict ;
  int i, iMax = dictMax (dict) ;
  char *cp, *cq ;

  for (cp = text ; (ispunct (*cp) || *cp == ' ') && *cp != '+' && *cp != '-' && !(*cp == '.' && isdigit(cp[1])) ; cp++) ;    /* gobble initial spaces */
  for (cq = cp + strlen(cp) - 1 ; cq > cp && (*cq == ' ' || (0 && ispunct (*cq) && *cp != ')')) ; cq--) ;
  cq[1] = 0 ; /* gobble terminal spaces */

  text = cp ;
  for (cq = cp = text ; *cp ; cp++) /* gobble \ */
    if (*cp != '\\') *cq++ = *cp ;
  *cq = 0 ;

  if (text && *text)
    {
      for (i = 1 ; i <= iMax ; i++)
	{
	  if (!strcasecmp (text, "male"))
	    { /* do not absorb male in female */
	      if (dictFind (dict, text, 0))
		return FALSE ;
	    }
	  else if (strstr (dictName(dict, i), text)) 
	    return FALSE ; /* do not double insert a text */
	}
      if (dictAdd (dict, text, &i))
	{
	  array (sra->tokens, arrayMax(sra->tokens), int) = i ;
	  return TRUE ;
	}
    }
  return FALSE ;
} /* sraTokenRegister */

/*************************************************************************************/
/* tokenise a phrase on , ; [] ()
 * add the token to the dict and to tt
 */
static void sraTextTokenize (SRA *sra, const char *text)
{
  AC_HANDLE h = ac_new_handle () ;
  char cc, *cp, *cp2, *cq ;
  char buf[1024] ;
  int nWords = 0 ;

  if (! text || ! *text)
    return ;
  if (dictFind (sra->dict, text, 0))
    return ; /* we do not need to store twce the same token */

  strncpy (buf, text, 1022) ; buf[1023] = 0 ;  /* make it editable */

  for (cp = buf ; *cp == ' ' ; cp++) ;  /* gobble initial spaces */
  for (cq = cp + strlen(cp) - 1 ; cq > cp && *cq == ' ' ; cq--) ;
  cq[1] = 0 ; /* gobble terminal spaces */

  for (cq = cp2 = cp, nWords = 0 ; *cq ; cq++)
    {
      switch ((int) *cq)
	{
	case ' ':
	  nWords++ ; cp2 = cq ;
	  break ;

	case '.': 
	  if (cq < cp2 + 4) /* ignore very short words. they are probably abbreviations */
	    break ;
	case ',':
	  if (isdigit (cq[1])) /* ignore decimal dots in 3.1416 */
	    break ; 
	case ';':
	    *cq = 0 ; sraTokenRegister (sra, cp) ; cp = cp2 = cq + 1 ; nWords = 0 ; 
	  break ;
	  
	case '(':
	case '[':
	  cc = *cq == '(' ? ')' : ']' ;
	  switch ((int)cq[1])  /* do not tokenize (+) (-) (none) (3.123) */
	    {
	    case '+' : cc = 0 ; break ;
	    case '-' : cc = 0 ; break ;
	    case '.' : if (isdigit (cq[2])) cc = 0 ; break ;
	    default : if (isdigit (cq[1])) cc = 0 ; break ;
	    }
	  if (cc == 0) break ;
	  *cq = 0 ;
	  sraTokenRegister (sra, cp) ;
	  cp = cq + 1 ;
	  cq = strchr (cp, cc) ;
	  if (cq) *cq = 0 ;
	  else cq = cp + strlen (cp) - 1 ;
	  sraTextTokenize (sra, cp) ;
	  cp = cp2 = cq + 1 ;
	  break ;
	}
    }

  if (cp && *cp && cp < cq)
    sraTokenRegister (sra, cp) ;

  ac_free (h) ;
  return ;
} /* sraTextTokenize */

/*************************************************************************************/
/* column 0 is the tag itself
 *        1 the direct tag content, equivalent to ac_tag_printable ()
 *        2 the tag_next content etc.
 */
static int sraGetData (SRA *sra, AC_OBJ obj, const char *tag, int column, int maxLength, AC_HANDLE h)
{
  int isMigs = 0 ;
  const char *tag2 = 0, **tag3 ;
  const char *inAuthor[] = {"collected by", "biomaterial provider", "BioSourceProvider", "identified by","INSDC center alias","INSDC center name", "collection date", "worm collection date", 0} ;

  if (column == -1)
    { 
      inAuthor[0] = 0 ;
      column = 2 ;
      tag2 = tag ;
      tag = "Biosample_attribute" ;
    } 
  if (strcmp ("Biosample_attribute", tag))
    {
      inAuthor[0] = 0 ;
    }
  if (obj && ac_has_tag (obj, tag))
    {
      AC_HANDLE h = ac_new_handle () ;
      AC_TABLE table = ac_tag_table (obj, tag, h) ;
      
      if (table && table->cols >= column)
	{
	  int ir ;

	  for (ir = 0 ; ir < table->rows ; ir++)
	    {
	      BOOL ok = TRUE ;
	      BOOL isPerfect = FALSE ;
	      const char *ccp = 0 ;

	      ccp = ac_table_printable (table, ir, 0, "___") ;
	      if (! strcmp (tag, "Description") &&
		  ! strncmp (ccp, "Protocols:", 10)
		  )
		continue ;
	      if (strcasestr (ccp, "External Id"))
		continue ;
	      if (strcasestr (ccp, "Broker name"))
		continue ;
	      if (strcasestr (ccp, "GEO Accession"))
		continue ;
	      if (strcasestr (ccp, "first-public"))
		continue ;
	      if (strcasestr (ccp, "first public"))
		continue ;
	      if (strcasestr (ccp, "last update"))
		continue ;
	      if (strcasestr (ccp, "last-update"))
		continue ;
	      if (tag2)
		{
		  if (strcmp (ccp, tag2))
		    continue ;
		}
	      for (tag3 = inAuthor ; ok && *tag3 ; tag3++) 
		{
		  if (! strcmp (ccp, *tag3))
		    ok = FALSE ;
		}
	      if (!ok)
		continue ;
		
	      if (column == 2)
		{
		  ccp = ac_table_printable (table, ir, column - 2, 0) ;
		  if (! strcasecmp (ccp, "genotype"))
		    isPerfect = TRUE ;
		}
	      ccp = ac_table_printable (table, ir, column - 1, 0) ;
	      if (! strcmp (tag, "Library_name"))
		{
		  if (! strncmp (ccp, "GSM", 3))
		    {
		      const char *cr ;
		      for (cr = ccp + 3 ; *cr ; cr++)
			{
			  if (*cr >= '0' && *cr <= 9)
			    continue ;
			  if (*cr == ' ' || *cr == ':' || *cr == ',')
			    continue ;
			}
		      ccp = *cr ? cr : 0 ;
		    }
		}
	      if ( !ccp)
		continue ;

	      if (isPerfect)
		{
		  int i ;
		  char buf[strlen(ccp)+1], *cq ;

		  for (cq = buf ; *ccp ; ccp++) /* gobble \ */
		    if (*ccp != '\\') *cq++ = *ccp ;
		  *cq = 0 ;
		  
		  dictAdd (sra->dict, buf, &i) ;
		  array (sra->tokens, arrayMax(sra->tokens), int) = i ;
		}
	      else
		{
		  ccp = sraRemoveBaddy (ccp, sra->baddies, h) ;
		  if ( !ccp) continue ;
		  if (ccp == (char *)1)
		    { isMigs = TRUE ; continue ; }
		  
		  if ( maxLength && strlen (ccp) >= maxLength) 
		    continue ;
		  sraTextTokenize (sra, ccp) ; /* tokens will be added to the sra->dict */
		}
	    } 
	}
      ac_free (h) ;
    }
  return isMigs ;
} /*  sraGetData */

/*************************************************************************************/
/* use the second entry to fix the capitatlization of the first entry */
static void sraFixUpper (SRA *sra, int t1, int t2)
{
  char buf1[301] ;
  char buf2[301] ;
  char *cp ;
  int i1, i2 ;
  int nUp1 = 0, nLow1 = 0 ;
  int nUp2 = 0, nLow2 = 0 ;
  
  i1 = array (sra->tokens, t1, int) ;
  i2 = array (sra->tokens, t2, int) ;

  strncpy (buf1, dictName (sra->dict, i1), 299) ;
  strncpy (buf2, dictName (sra->dict, i2), 299) ;
  buf1[300] = 0 ;
  buf2[300] = 0 ;
  
  /* count upper/lower case letters */
  cp = buf1 - 1 ;
  while (*++cp)
    {
      if (*cp == ace_upper (*cp)) nUp1++ ;
      else nLow1++ ;
    }
  cp = buf2 - 1 ;
  while (*++cp)
    {
      if (*cp == ace_upper (*cp)) nUp2++ ;
      else nLow2++ ;
    }
  
   if (nUp2 && nLow2 && nLow1 == 0) /* prefer lower case */
     {
       cp = strcasestr (buf1, buf2) ;
       if (cp)
	 {
	   memcpy (buf1 + (cp - buf1), buf2, strlen(buf2)) ;
	   dictAdd (sra->dict, buf1, &i1) ;
	   array (sra->tokens, t1, int) = i1 ;
	 }
     }
} /* sraFixUpper */

/*************************************************************************************/
/* if a token includes another one, remove the shorter one
 * but be careful to select the best level of capitalization
 */
static void sraFilterTokens (SRA *sra)
{
  DICT *dict = sra->dict ;
  Array tokens = sra->tokens ;
  int t1, t2, i1, i2, tMax = arrayMax (tokens) ;
  const char *ccp1, *ccp2 ;
  BOOL found ;

  /* work backwords
   * if token 7 is included in any token above, kill it
   * if token 21 replaces token 7, kill any token 8 to 20 matching 21, then kill 21 
   */

  for (t2 = 1 ; t2 < tMax ; t2++)
    {
      i2 = array (tokens, t2, int) ;
      if (! i2) continue ;
      ccp2 = dictName (dict, i2) ;
  
      for (t1 = 0 ; i2 && t1 < t2 ; t1++)
	{     /* check if ccp2 is included in previous tokens */
	  i1 = array (tokens, t1, int) ;
	  if (! i1) continue ;
	  ccp1 = dictName (dict, i1) ;
	  
	  if (strcasestr (ccp1, ccp2) &&
	      strcasecmp (ccp2, "male")
	      ) /* drop t2 */
	    {
	      sraFixUpper (sra, t1, t2) ;
	      i2 = array (tokens, t2, int) = 0 ;
	    }
	}
      
      for (t1 = 0, found = FALSE ; i2 && t1 < t2 ; t1++)
	{     /* check if ccp2 is includes previous tokens */
	  i1 = array (tokens, t1, int) ;
	  if (! i1) continue ;
	  ccp1 = dictName (dict, i1) ;
	  
	  if (strcasestr (ccp2, ccp1) &&
	      strcasecmp (ccp1, "male")) /* improve t1 */
	    {
	      if (! found)
		{
		  found = TRUE ;
		  sraFixUpper (sra, t2, t1) ;
		  array (tokens, t1, int)  = i2 ;
		}
	      else
		array (tokens, t1, int) = 0 ;
	    }
	}
      if (found)
	array (tokens, t2, int) = 0 ;
    }
} /* sraFilterTokens */

/*************************************************************************************/

static char *sraMergeTokens (SRA *sra, vTXT txt)
{
  Array tt = sra->tokens ;
  int t, sep = 0 ;
  
  for (t = 0 ; t < arrayMax (tt) ; t++)
    {
      int t1 = arr (tt, t, int) ;
      if (t1)
	sep = sraPrint (txt, sra->dict, t1, sep) ;
    }
  return vtxtPtr (txt) ;
} /* sraMergeTokens */

/*************************************************************************************/

static void sraExportTokens (SRA *sra, AC_OBJ obj, const char *tag)
{
  AC_HANDLE h = ac_new_handle () ; 
  vTXT txt = vtxtHandleCreate (h);
     
  sraFilterTokens (sra) ; 
  sraMergeTokens (sra, txt) ;
  if (obj)
    {
      if (sra->dbEdit)
	{
	  vTXT txt1 = vtxtHandleCreate (h);
	  const char *errors = 0 ;

	  vtxtPrintf (txt1, "%s ", ac_class (obj)) ;
	  vtxtPrint (txt1, ac_protect (ac_name(obj), h)) ;
	  if (vtxtPtr (txt))
	    {
	      int delta = 0 ;
	      vtxtPrintf (txt1, "\n%s ", tag) ;
	      {
		const char *ccp = vtxtPtr (txt) ;
		if (! strncmp (ccp, "GSM", 3))
		  {
		    ccp += 3 ;
		    while (*ccp >= '0' && *ccp <= '9') ccp++ ;
		    if (*ccp == ':')
		      delta = ccp - vtxtPtr (txt) + 1 ;
		  }
		vtxtPrint (txt1, ac_protect (vtxtPtr (txt) + delta, h)) ;
	      }
	    }
	  else
	    vtxtPrintf (txt1, "\n-D %s ", tag) ;
	  vtxtPrint (txt1, "\n\n") ;
	  ac_parse (sra->db, vtxtPtr (txt1), &errors, 0, h) ;
	  if (* errors)
	    fprintf (stderr, "%s\n", errors) ;
	}
      else
	{
	  printf ("%s\t%s\t%s\n", tag, ac_name(obj), vtxtPtr (txt) ? vtxtPtr (txt)  : "-") ;
	}
    }

  arrayMax (sra->baddies) = 0 ;
  arrayMax (sra->tokens) = 0 ;
  ac_free (sra->dict) ;
  sra->dict = dictHandleCreate (128, sra->h) ;
  sra->dict = dictCaseSensitiveHandleCreate (128, sra->h) ;
  ac_free (h) ;
} /* sraExportTokens  */

/*************************************************************************************/
/*************************************************************************************/

static void sraGetPaper5Authors (SRA *sra)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE papers, aut ;
  int ir, jr ;

  papers = sra->srp ? ac_tag_table (sra->srp, "Reference", h) : 0 ;
  for (ir = 0 ; ir < (papers ? papers->rows : 0) ; ir++)
    {
      AC_OBJ pap = ac_table_obj (papers, ir, 0, h) ;
      aut =  ac_tag_table (pap, "Author", h) ;
      if (! aut || aut->rows < 1)
	continue ;
      for (jr = 0 ; jr < aut->rows ; jr++)
	{
	  const char *ccp ;
	  if (aut->rows > 5 && jr == 3)
	    {
	      jr = aut->rows - 2 ;
	      sraTextTokenize (sra, "---") ;
	    }
	  ccp = ac_table_printable (aut, jr, 0, 0) ;
	  if (ccp)
	    {
	      ccp = sraRemoveBaddy (ccp, sra->baddies, h) ; 
	      if (ccp > (char *)1)
		sraTextTokenize (sra, ccp) ;
	    }
	}
      if (aut->rows > 5)
	sraTextTokenize (sra, "et al.") ;
    }
  ac_free (h) ;
} /* sraGetPaper5Authors */

/*************************************************************************************/

static void sraAuthor2 (SRA *sra)
{
  AC_HANDLE h = ac_new_handle () ;
  
  sraSetBaddies (sra, 'a') ;
  
  /* notice the hierarchy g s c p x : geo is consider best, then sample etc 
   * the tokenization is hierarchic
   */
  sraGetData (sra, sra->geo , "Author", 1, 0, h) ;
  sraGetPaper5Authors (sra) ;
  sraGetData (sra, sra->srp , "Identifier", 1, 0, h) ;
  sraGetData (sra, sra->sample , "collected by", -1, 0, h) ;
  sraGetData (sra, sra->sample , "biomaterial provider", -1, 0, h) ;
  sraGetData (sra, sra->sample , "BioSourceProvider", -1, 0, h) ;
  sraGetData (sra, sra->sample , "identified by", -1, 0, h) ;
  sraGetData (sra, sra->sample , "INSDC center alias", -1, 0, h) ;
  sraGetData (sra, sra->sample , "INSDC center name", -1, 0, h) ;
  sraGetData (sra, sra->srr , "Center_name", 1, 0, h) ;

  sraGetData (sra, sra->sample , "Submission", 1, 0, h) ;
  sraGetData (sra, sra->srx , "Submitted_by", 1, 0, h) ;

  sraExportTokens (sra, sra->srr, "Magic_author2") ;

  ac_free (h) ;
} /* sraAuthor2  */

/*************************************************************************************/

static void sraSample2 (SRA *sra)
{
  AC_HANDLE h = ac_new_handle () ;
  int isMigs = 0 ;
  sraSetBaddies (sra, 's') ;
  sraMaskTagContent (sra, sra->srr, "Species", h) ;
  
  isMigs += sraGetData (sra, sra->sample , "Title", 1, 100, h) ;
  isMigs += sraGetData (sra, sra->srr , "Library_name", 1, 200, h) ;
  isMigs += sraGetData (sra, sra->sample , "Biosample_attribute", 2, 100, h) ;
  isMigs += sraGetData (sra, sra->sample , "Identifier", 2, 100, h) ;
  isMigs += sraGetData (sra, sra->sample , "Description", 1, 300, h) ;
  
  sraExportTokens (sra, sra->srr, "Magic_sample2") ;
  
  ac_free (h) ;
} /* sraSample2  */

/*************************************************************************************/
/*************************************************************************************/

static void sraAnalayze (SRA *sra)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
  const char *query ;
  AC_OBJ srr = 0 ;
  AC_ITER iter = 0 ;
  
  if (sra->project)
    query = hprintf (h, "Find project %s ; >SRR ; SPR || SRX || Biosample || GEO ", sra->project) ;
  else
    query = hprintf (h, "Find SRR SPR || SRX || Biosample || GEO ; %s", sra->template ? sra->template : "") ;
  iter = ac_query_iter (sra->db, TRUE, query, 0, h) ;

  if (0 && sra->showAuthor && ! sra->dbEdit)
    printf ("# SRR\tCenter_name\tSRP Identifier\tSRX submitted_by\tSample submission\tGeo author\t\tProposed title\n") ;

  while (ac_free (srr), (srr = ac_iter_obj (iter)))
    {
      AC_HANDLE h1 = ac_new_handle () ;
      nn++ ;
      sra->srr = srr ;
      sra->srp = ac_tag_obj (srr, "SRP", h1) ;
      sra->srx = ac_tag_obj (srr, "SRX", h1) ;
      sra->geo =sra->srp ?  ac_tag_obj (sra->srp, "GEO", h1) : 0 ;
      sra->sample = ac_tag_obj (srr, "Biosample", h1) ;

      if (1)
	{
	  if (sra->showAuthor) sraAuthor2 (sra) ;
	  if (sra->showSample) sraSample2 (sra) ;
	}
      ac_free (h1) ;
      if (sra->max && nn >= sra->max) break ;
    }

  ac_free (h) ;
} /* sraAnalyze  */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *err)
{
  fprintf (stderr, "// Usage: sra_metadata : semantic synthetizer of SRA metadata\n") ;
  fprintf (stderr, 
	   "// Example:\n"
           "//    sra_metadata -db <ACEDB> -q 'SRR4907*' -max 7 -a -s [-dbEdit]\n"
	   "// Options:\n"
	   "//    -db <ACEDB> : mandatory\n"
	   "//      an acedb database with classes SRR, SRP, SRX...\n"
	   "//      as defined by the schema in waligner/metaData/wspec.SRA\n"
	   "//      and populated by the scripts waligner/scripts/SRX_import*\n"
	   "//    -a : show condensed author\n"
	   "//    -s : show condensed sample\n"
	   "//    -p <MAGIC> : limit to SRR in given project\n"
	   "//    -q query : limit to SRR matching the query\n"
	   "//    -max <int>  : optional limit while debugging, default unlimited\n"
	   "//       analyse at most max SRR objects maximal number\n"
	   "//     -dbEdit : edit the database (rather than export the titles)\n"
	   ) ;

  if (err)
    fprintf (stderr, "\n// ERROR: %s\n", err) ;
  exit (1) ;
}

/*************************************************************************************/

int main (int argc, const char **argv)
{
  AC_HANDLE h = ac_new_handle () ;
  SRA sra ;
  const char *errors = 0 ;
  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  memset (&sra, 0, sizeof (SRA)) ;
  if (! getCmdLineOption (&argc, argv, "-db", &(sra.dbNam)))
    usage ("Missing -argument -db") ;

  getCmdLineInt (&argc, argv, "-max", &(sra.max)) ;
  getCmdLineOption (&argc, argv, "-q", &(sra.template)) ;
  getCmdLineOption (&argc, argv, "-p", &(sra.project)) ;

  sra.showAuthor = getCmdLineBool (&argc, argv, "-a") ;
  sra.showSample = getCmdLineBool (&argc, argv, "-s") ;
  sra.dbEdit = getCmdLineBool (&argc, argv, "-dbEdit") ;
  sra.old = getCmdLineBool (&argc, argv, "-old") ;
  /*
    getCmdLineOption (&argc, argv, "-o", &(jp.outFileName)) ;
  getCmdLineOption (&argc, argv, "-t", &(jp.tFileName)) ;
  getCmdLineOption (&argc, argv, "-i", &(jp.pFileName)) ;
  getCmdLineInt (&argc, argv, "-nFilters", &jp.NFILTERS) ;
  */

  if (argc > 1)
    usage (messprintf ("Unknown parameter %s", argv[1])) ;

  sra.h = ac_new_handle () ;
  sra.dict = dictCaseSensitiveHandleCreate (1000, h) ;
  /* sra.dict = dictHandleCreate (1000, h) ; */
  sra.baddies = arrayHandleCreate (32, char *, h) ;
  sra.tokens = arrayHandleCreate (32, int, h) ;
  fprintf (stderr, "// %s: sra_metadata starts\n", timeShowNow ()) ;
 
  sraUncap ("test", h) ; /* to please the compiler */

  /* actual wort */ 
  sra.db = ac_open_db (sra.dbNam, &errors) ;
  if (errors)
    messcrash ("Sorry, i cannot open the database %s : ERROR %s\n", sra.dbNam, errors) ;
      
  sraAnalayze (&sra) ;
  ac_db_close (sra.db) ;

  /* done */
  fprintf (stderr, "// %s: sra_metadata done\n ", timeShowNow ()) ;
 
  ac_free (h) ;
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

