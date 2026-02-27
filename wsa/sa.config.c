/*
 * sa.config.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module parses the various configuration files
 * which define the sequence files to be aligned, the target genome
 * and its annotations (gens, exons, introns ...)
 */

#include "sa.h"

static void saSlInit (PP *pp)
{
  /* if you edit the SL table, also edit the corresponding table in solexa.c */
  char *slu[] = { "GGTTTAATTACCCAAGTTTGAG" ,    /* SL1 */
	 	  "GGTTTTAACCCAGTTACTCAAG" ,    /* SL2 */
	 	  "GGTTTTAACCCAGTTAACCAAG" ,    /* SL3 */
		  "GGTTTTAACCCAGTTTAACCAAG" ,    /* SL4 */
     		  "GGTTTTAACCCAGTTACCAAG" ,    /* SL5 */
 	 	  "GGTTTAAAACCCAGTTACCAAG" ,    /* SL6 */
 	 	  "GGTTTTAACCCAGTTAATTGAG" ,    /* SL7 */
 		  "GGTTTTTACCCAGTTAACCAAG" ,    /* SL8 */
		  "GGTTTATACCCAGTTAACCAAG" ,    /* SL9 */
		  "GGTTTTAACCCAAGTTAACCAAG" ,    /* SL10 */
 		  "GGTTTTAACCAGTTAACTAAG" ,    /* SL11 */
		  "GGTTTTAACCCATATAACCAAG" ,    /* SL12 */
		  /* existe pas
		   * "GTTTTTAACCCAGTTACTCAAG" ,   SL13 
 	 	   * "GGTTTTTAACCCAGTTACTCAAG" ,   SL14 
		   */
		  0 } ;

  pp->SLs = arrayHandleCreate (SLMAX, char*, pp->h) ;  

  for (int i = 0 ; i < SLMAX && slu[i] ; i++)
    {
      char *cs = slu[i] ;
      if (! cs || ! *cs) break ;
      int k = strlen (cs) ;
      char *buf = halloc (k+1, pp->h) ;
      for (int j = 0 ; j < k ; cs++, j++)
	{
	  char c = dnaEncodeChar [(int)*cs] ;
	  buf[j] = c ;
	}
      buf[k] = 0 ;
      array (pp->SLs, i, char *) = buf ;
    }
  return ;
} /* saSlInit */

/*********************************************************************/

void saConfigIsIndexAccessible (PP *pp)
{
  if (! pp->indexName)
    saUsage ("Missing parameter --index", 0, 0) ;
  else
    {
      AC_HANDLE h = ac_new_handle () ;
      const char *cp = pp->indexName ;
      char *cq ;
      if (pp->createIndex)
	{
	  cq = hprintf (h, "mkdir %s ; \\rm %s/* ; echo test > %s/sortalign.testFile", cp, cp, cp) ;
	  system (cq) ;
	}
      cq = hprintf (h, "%s/sortalign.testFile", pp->indexName) ;
      ACEIN ai = aceInCreate (cq, FALSE, h) ;
      if (! ai)
	saUsage (hprintf (h, "Cannot create and access the index directory: %s", pp->indexName), 0, 0) ;
      ac_free (h) ;
    }
} /* saConfigIsIndexAccessible */

/*********************************************************************/

void saConfigIsOutDirWritable (PP *pp)
{
  if (pp->outFileName)
    {
      AC_HANDLE h = ac_new_handle () ;

      char *buf = strnew (pp->outFileName, h) ;
      char *cp = buf + strlen(buf) - 1 ;

      while (cp >= buf && *cp == '/') *cp-- = 0 ;
      
      if (buf[0])
	{
	  char *cmd ;
	  if (pp->wiggle)
	    cmd = hprintf (h, "mkdir -p %s  %s/hits %s/wiggles", buf, buf, buf) ;
	  else
	    cmd = hprintf (h, "mkdir -p %s %s/hits", buf, buf) ;
	  system (cmd) ;
	}
      pp->outFileName = hprintf (pp->h, "%s/", buf) ;
      
      ACEOUT ao = aceOutCreate (buf, "/test", 0, h) ;
      if (ao)
	{
	  int n = randint (), n2 = 0 ;
	  char * nam = strnew (aceOutFileName (ao), h) ;
	  
	  aceOutf (ao, "%d\n", n ) ;
	  ac_free (ao) ;
	  
	  ACEIN ai = aceInCreate (nam, 0, h) ;
	  if (ai)
	    {
	      if (aceInCard (ai))
		aceInInt (ai, &n2) ;
	      if (n != n2)
		saUsage (hprintf (h, "Unable to open a test file -o %s.test", pp->outFileName), 0, 0) ;
	      ac_free (ai) ;
	    }
	}
      ac_free (h) ;
    }
  
  return ;
} /* saConfigIsOutFDirWritable */

/*********************************************************************/

/* check the existence of the sequence files to be aligned
 * identify their absolute file names
 * associate each file to a run name and to its optional parameters
 */
Array saConfigGetRuns (PP *pp, Array runStats)
{
  Array rcs = arrayHandleCreate (64, RC, pp->h) ;
  RC *rc = 0 ;
  int nn = 0 ;
  BOOL isSet = FALSE ;
  
  if(pp->isWorm)
    saSlInit (pp) ;

  if (pp->rawAdaptor2R)  /* we expect this adaptor to be the minus strand exit adaptor of read 2, as exported by our adaptor profile */
    {
      int i, iMax = strlen (pp->rawAdaptor2R) ;
      isSet = TRUE ;
      for (i = 0 ; i < iMax && i < 30; i++)
	pp->adaptors.a1L[i] = dnaEncodeChar[(int)pp->rawAdaptor2R[i]] ;
    }
  
  if (pp->rawAdaptor1R) /* we expect this adaptor to be the top strand exit adaptor of read 1, as exported by our adaptor profile */
    {
      int i, iMax = strlen (pp->rawAdaptor1R) ;
      isSet = TRUE ;
      for (i = 0 ; i < iMax && i < 30; i++)
	pp->adaptors.a2L[i] = dnaEncodeChar[(int)pp->rawAdaptor1R[i]] ;
    }
  
  if (pp->rawAdaptor1R)  /* we expect this adaptor to be the top strand exit adaptor of read 1, as exported by our adaptor profile */
    {
      int i, iMax = strlen (pp->rawAdaptor1R) ;
      isSet = TRUE ;
      for (i = 0 ; i < iMax && i < 30; i++)
	pp->adaptors.a1R[i] = dnaEncodeChar[(int)pp->rawAdaptor1R[i]] ;
    }
  
  if (pp->rawAdaptor2R)  /* we expect this adaptor to be the minus strand exit adaptor of read 2, as exported by our adaptor profile */
    {
      int i, iMax = strlen (pp->rawAdaptor2R) ;
      isSet = TRUE ;
      for (i = 0 ; i < iMax && i < 30; i++)
	pp->adaptors.a2R[i] = dnaEncodeChar[(int)pp->rawAdaptor2R[i]] ;
    }
  
  if (isSet)
    saSetGetAdaptors (2, 0, &(pp->adaptors), 0) ;  /* run=0: valid for all runs */
  
  if (pp->inFileName)
    {   /* Split the individual file names, they are coma separated 
	 * Make the names absolute, and load them in the file-parser fp channel
	 *  all names must be loaded at once, otherwise we may have a deadlock
	 *  because the agents needing genome parser completion are not yet launched 
	 */
      char *buf = strnew (pp->inFileName, pp->h) ;
      char *cp, *cq, *cr, *cs ;
      int run = 0 ;
      char *filName2 = 0 ;
      if (! pp->runName)
	pp->runName = "runX" ;
      dictAdd (pp->runDict, pp->runName, &run) ;
      cp = buf ; *(cp + strlen (cp)) = ',' ;
      while ((cq = strchr (cp, ',')))
	{
	  cs = strchr (cp, ':') ;
	  if (cs)
	    { /* switch runName */
	      *cs = 0 ;
	      dictAdd (pp->runDict, cp, &run) ;
	      cp = cs + 1 ;
	    }
	  *cq = 0 ;
	  if (! *cp)
	    continue ;
	  cs = strchr (cp, '+') ;
	  if (cs)
	    {
	      *cs = 0 ;
	      cr = filName (cs+1, 0, "r") ;
	      if (! cr)
		saUsage (hprintf (pp->h, "Cannot open input file %s\n", cs+1), 0, 0) ;
	      filName2 = strnew (cr, pp->h) ;
	    }
	  

	  rc = arrayp (rcs, nn++, RC) ; 
	  rc->fileName1 = strnew (cp, pp->h) ;
	  rc->fileName2 = filName2 ;
	  rc->run = run ;

	  rc->format = FASTA ; /* default */
	  int ln = strlen (rc->fileName1) ;
	  if (strstr (rc->fileName1, ".sra.fasta.gz") == rc->fileName1 + ln - 13) rc->format = SRACACHE ;
	  else if (strstr (rc->fileName1, ".sra.fasta") == rc->fileName1 + ln - 10) rc->format = SRACACHE ;
	  else if (strstr (rc->fileName1, ".fasta.gz") == rc->fileName1 + ln - 9) rc->format = FASTA ;
	  else if (strstr (rc->fileName1, ".fastq.gz") == rc->fileName1 + ln - 9) rc->format = FASTQ ;
	  else if (strstr (rc->fileName1, ".fasta") == rc->fileName1 + ln - 6) rc->format = FASTA ;
	  else if (strstr (rc->fileName1, ".fastq") == rc->fileName1 + ln - 6) rc->format = FASTQ ;
	  else if (strstr (rc->fileName1, ".fna.gz") == rc->fileName1 + ln - 7) rc->format = FASTA ;
	  else if (strstr (rc->fileName1, ".fa.gz") == rc->fileName1 + ln - 6) rc->format = FASTA ;
	  else if (strstr (rc->fileName1, ".fastc.gz") == rc->fileName1 + ln - 9) rc->format = FASTC ;
	  else if (strstr (rc->fileName1, ".fastc") == rc->fileName1 + ln - 6) rc->format = FASTC ;
	  else if (strstr (rc->fileName1, ".fna") == rc->fileName1 + ln - 4) rc->format = FASTA ;
	  else if (strstr (rc->fileName1, ".fa") == rc->fileName1 + ln - 3) rc->format = FASTA ;
	  else if (! strncmp (rc->fileName1, "SRR", 3)) rc->format = SRA ;
	  else if (! strncmp (rc->fileName1, "DRR", 3)) rc->format = SRA ;
	  else if (! strncmp (rc->fileName1, "ERR", 3)) rc->format = SRA ;
	  
	  /* user can override the defaults */
	  if (pp->raw) rc->format = RAW ;
	  if (pp->fasta) rc->format = FASTA ;
	  if (pp->fastq) rc->format = FASTQ ;
	  if (pp->fastc) rc->format = FASTC ;
	  if (pp->sra) rc->format = SRA ;

	  if (pp->isRna) /* user imposed */
	    {
	      int isRna = 1 ;
	      saSetGetAdaptors (2, &isRna, 0, rc->run) ;
	    }
	  if (pp->isDna) /* user imposed */
	    {
	      int isRna = -1 ;
	      saSetGetAdaptors (2, &isRna, 0, rc->run) ;
	    }

	  if (rc->format == SRA)  /* check in the cache */
	    {
	    }
	  else
	    {
	      cr = filName (rc->fileName1, 0, "r") ;
	      if (! cr)
		saUsage (hprintf (pp->h, "Cannot open input file %s\n", cp), 0, 0) ;
	      rc->fileName1 = strnew (cr, pp->h) ;
	    }
	  
	  RunSTAT *rs = arrayp (runStats, run, RunSTAT) ;
	  rs->p.nFiles++ ;
	  if (rc->fileName2)
	    {
	      rc->pairedEnd = TRUE ;
	      rs->p.nFiles++ ;
	    }
	  if (rc->format == FASTC)
	      rc->pairedEnd = TRUE ;

	  cp = cq + 1 ;
	}
    }
  else if (pp->inConfigFileName)
    {
      AC_HANDLE h = ac_new_handle () ;
      DICT *fDict = dictHandleCreate (64, h) ;
      ACEIN ai = aceInCreate (pp->inConfigFileName, 0, h) ;
      int nRuns = 0, run = 0 ;
      int line = 0 ;

      if (! ai)
	saUsage (hprintf (h, "Cannot find file -I %s specified on the command line", pp->inConfigFileName), 0, 0) ;
      while (aceInCard (ai))
	{
	  char *cq, *cr, *cp = aceInWord (ai) ;
	  
	  line++ ;
	  if (! cp || ! *cp || *cp == '#')
	    continue ;
	  /* file names */
	  rc = arrayp (rcs, nn++, RC) ; 
	  nRuns++ ;
	      /* file pairs */
	  cq = strchr (cp, '+') ;
	  if (cq) 
	    {
	      *cq++ = 0 ;
	      if (*cq)
		{
		  cr = filName (cq, 0, "r") ;
		  if (! cr)
		    saUsage (hprintf (h, "Cannot open input file %s\n", cq), 0, 0) ;
		  rc->pairedEnd = TRUE ;
		  rc->fileName2 = strnew (cr, pp->h) ;
		}
	    }
	  rc->fileName1 = strnew (cp, pp->h) ;
	  if (! dictAdd (fDict, cp, 0))
	    saUsage (hprintf (h, "Duplicate file name %s\n at line %d of file -T %s"
			      , cr
			      , line
			      , pp->inConfigFileName
			      ), 0, 0) ;

	  /* run name */
	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;
	  if (cp && *cp && *cp != '#')
	    dictAdd (pp->runDict, cp, &run) ;
	  else
	    dictAdd (pp->runDict, hprintf (h, "r.%d", nRuns), &run) ;
	  rc->run = run ;

	  /*
	    memcpy (rc->adaptor1L, pp->adaptor1L , 30) ;
	    memcpy (rc->adaptor2L, pp->adaptor2L , 30) ;
	    memcpy (rc->adaptor1R, pp->adaptor1R , 30) ;
	    memcpy (rc->adaptor2R, pp->adaptor2R , 30) ;
	  */
	  /* options */
	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;
	  while (cp)
	    {
	      int k = 0 ;
	      
	      if (*cp == '#')
		break ;
	      cq = strchr (cp, ',') ;
	      if (cq)
		*cq++ = 0 ;
	      if (! strcasecmp (cp, "fasta")) rc->format = FASTA ;
	      else if (! strcasecmp (cp, "fastq")) rc->format = FASTQ ;
	      else if (! strcasecmp (cp, "fastc")) rc->format = FASTC ;
	      else if (! strcasecmp (cp, "raw")) rc->format = RAW ;
	      else if (! strcasecmp (cp, "SRA")) rc->format = SRA ;

	      else if (! strcasecmp (cp, "rna")) /* user imposed */
		{
		  int isRna = 1 ;
		  saSetGetAdaptors (2, &isRna, 0, rc->run) ;
		}

	      else if (! strcasecmp (cp, "dna")) /* user imposed */
		{
		  int isRna = -1 ;
		  saSetGetAdaptors (2, &isRna, 0, rc->run) ;
		}

	      else if (sscanf (cp, "jump1=%d", &k) && k >0)
		rc->jump5r1 = k ;
	      else if (sscanf (cp, "jump2=%d", &k) && k >0)
		rc->jump5r2 = k ;
		
		  
	      cp = cq ;
	    }

	  if (! rc->format)
	    {
	      int ln = strlen (rc->fileName1) ;
	      rc->format = FASTA ; /* default */
	      if (strstr (rc->fileName1, ".sra.fasta.gz") == rc->fileName1 + ln - 13) rc->format = SRACACHE ;
	      else if (strstr (rc->fileName1, ".sra.fasta") == rc->fileName1 + ln - 10) rc->format = SRACACHE ;
	      else if (strstr (rc->fileName1, ".fasta.gz") == rc->fileName1 + ln - 9) rc->format = FASTA ;
	      else if (strstr (rc->fileName1, ".fastq.gz") == rc->fileName1 + ln - 9) rc->format = FASTQ ;
	      else if (strstr (rc->fileName1, ".fasta") == rc->fileName1 + ln - 6) rc->format = FASTA ;
	      else if (strstr (rc->fileName1, ".fastq") == rc->fileName1 + ln - 6) rc->format = FASTQ ;
	      else if (strstr (rc->fileName1, ".fna.gz") == rc->fileName1 + ln - 7) rc->format = FASTA ;
	      else if (strstr (rc->fileName1, ".fa.gz") == rc->fileName1 + ln - 6) rc->format = FASTA ;
	      else if (strstr (rc->fileName1, ".fastc.gz") == rc->fileName1 + ln - 9) rc->format = FASTC ;
	      else if (strstr (rc->fileName1, ".fastc") == rc->fileName1 + ln - 6) rc->format = FASTC ;
	      else if (strstr (rc->fileName1, ".fna") == rc->fileName1 + ln - 4) rc->format = FASTA ;
	      else if (strstr (rc->fileName1, ".fa") == rc->fileName1 + ln - 3) rc->format = FASTA ;
	      else if (! strncmp (rc->fileName1, "SRR", 3)) rc->format = SRA ;
	      else if (! strncmp (rc->fileName1, "DRR", 3)) rc->format = SRA ;
	      else if (! strncmp (rc->fileName1, "ERR", 3)) rc->format = SRA ;
	    }

	  if (rc->format == SRA)  /* check in the cache */
	    {
	    }
	  else
	    {
	      cr = filName (rc->fileName1, 0, "r") ;
	      if (! cr)
		saUsage (hprintf (pp->h, "Cannot open sequence file %s listed in configuration file -I %s \tDid you mean -i %s\n"
				  , cp ? cp : "NULL", pp->inConfigFileName, pp->inConfigFileName
				  ), 0, 0) ;
	      rc->fileName1 = strnew (cr, pp->h) ;
	    }
	  
	  RunSTAT *rs = arrayp (runStats, run, RunSTAT) ;
	  rs->p.nFiles++ ;
	  if (rc->fileName2)
	    rs->p.nFiles++ ;
	  if (rc->fileName2 || rc->format == FASTC)
	      rc->pairedEnd = TRUE ;
	} 
      
      ac_free (h) ;
    }
  pp->nFiles = arrayMax (rcs) ;
  printf ("Found %d sequence file%s\n", arrayMax (rcs), arrayMax (rcs) > 1 ? "s" : "") ;
  return rcs ;
} /* saConfigGetRuns */

/*************************************************************************************/
/*************************************************************************************/
/**
 * Thread-safe accessor for two fixed-length adaptor sequences (31 chars + null).
 * If set == 2, stores copies of aa->adaptor1 and aa->adaptor2. Non re-writable
 * If set == TRUE, stores copies of aa->adaptor1 and aa->adaptor2.
 * If set == FALSE, copies stored values back to *aa (if previously set).
 * Returns TRUE if adaptors have been set at least once, FALSE otherwise.
 */

static pthread_mutex_t adaptor_mutex = PTHREAD_MUTEX_INITIALIZER;

BOOL saSetGetAdaptors (int set, int *isRnap, ADAPTORS *aa, int run)
{
  static AC_HANDLE h = 0 ;
  static Array aaa = 0 ;
  static Array isSets = 0 ;
  static Array isRnas = 0 ;
  int isSet = TRUE ;
  
  pthread_mutex_lock(&adaptor_mutex);

  if (aaa == 0)
    {
      h = ac_new_handle () ;
      aaa = arrayHandleCreate (16, ADAPTORS, h) ;
      isSets = arrayHandleCreate (16, int, h) ;
      isRnas = arrayHandleCreate (16, int, h) ;
    }
  if (run < 0)
    messcrash ("\nBad call to saSetGetAdaptors run = %d < 0", run) ;
  
  ADAPTORS *aR = arrayp (aaa, run, ADAPTORS) ; /* must come first */
  ADAPTORS *a0 = arrayp (aaa, 0, ADAPTORS) ;  /* run may reallocate aa */
  int *isSetRp = arrayp (isSets, run, int) ;
  int *isSet0p = arrayp (isSets, 0, int) ;
  int *isRnaRp = arrayp (isRnas, run, int) ;
  int *isRna0p = arrayp (isRnas, 0, int) ;

  if (set == 0)
    {  /* get */
      int iss = 0 ;
      if (aa) 
	{
	  memset (aa, 0, sizeof (ADAPTORS)) ;
	  if (*isSet0p) /* user imposed gloablly */
	    { iss++ ; memcpy (aa, a0, sizeof (ADAPTORS)) ; }
	  else if (*isSetRp) /* run specific */
	    { iss++ ; memcpy (aa, aR, sizeof (ADAPTORS)) ; }
	}
      if (isRnap)
	{
	  *isRnap = 0 ; /* not set */
	  if (*isRna0p) /* user imposed globally */
	    { iss++ ; *isRnap = *isRna0p ; }
	  if (*isRnaRp) /* run specific */
	    { iss++ ; *isRnap = *isRnaRp ; }
	}
      /* client code wish to block if we return FALSE */
      isSet = (iss > 0 ? TRUE : FALSE) ;
    }

  else if (set == 2)
    { /* user defined */
      if (isRnap)
	*isRnaRp = (*isRnap > 0 ? 9999999 : -9999999) ; 
      if (aa)
	{
	  *isSetRp = 2 ;
	  memcpy (aR, aa, sizeof (ADAPTORS)) ;
	}
      isSet = TRUE ;
    }

  else if (set == 1)
    { /* guess work do not modify previous user set values */
      
      if (! *isRnaRp)
	*isRnaRp = 1 ; /* odd, all future modifs are even, so we allways stay set (non zero) */
      if (isRnap && *isRnap)
	*isRnaRp += (*isRnap > 0 ? 2 : -2) ;
      if (isRnap)
	*isRnap = *isRnaRp ;
      if (aa)
	{
	  switch (*isSetRp)
	    {
	    case 2:
	      break ;
	    case 0:
	    case 1:
	      *isSetRp = 1 ;
	      memcpy (aR, aa, sizeof (ADAPTORS)) ;
	      break ;
	    }
	}
      isSet = TRUE ;
    }
  else
    messcrash ("Bad call set=%d to saSetGetAdaptors", set) ;


  pthread_mutex_unlock (&adaptor_mutex) ;
  return isSet ;
} /* saSetGetAdaptors */

/*************************************************************************************/
/*************************************************************************************/
/*********************************************************************/
/*********************************************************************/

/* check the completeness of the target index directory
 */
int saConfigCheckTargetIndex (PP *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL ok = TRUE ;
  int k = 0 ;
  char *cp ;
  int NN = 0 ;
  while (1)
    {
      char *fNam = hprintf (h, "/cws.sortali.%d", k) ;
      k++ ;
      cp = filName (pp->indexName, fNam, "rb") ;
      if (cp)
	{
	  NN++;
	  if (NN == 1)
	    {
	      char *cq = cp + strlen(cp) - 2 ;
	      *cq = 0 ;
	      pp->tFileBinaryCwsName = strnew (cp, pp->h) ;
	    }
	}
      else
	break ;
    }
  if (NN != 1 && NN != 2 && NN != 4 && NN != 8 && NN != 16 && NN != 32 && NN != 64)
    messcrash ("The number of indexes NN=%d must be  a power of 2, say 1 2 4 8 ...", NN) ;

  cp = filName (pp->indexName, "/dna.sortali", "rb") ;
  if (cp)
    pp->tFileBinaryDnaName = strnew (cp, pp->h) ;
  else
    ok = FALSE ;
  cp = filName (pp->indexName, "/dnaR.sortali", "rb") ;
  if (cp)
    pp->tFileBinaryDnaRName = strnew (cp, pp->h) ;
  else
    ok = FALSE ;
  cp = filName (pp->indexName, "/ids.sortali", "rb") ;
  if (cp)
    pp->tFileBinaryIdsName = strnew (cp, pp->h) ;
  else
    ok = FALSE ;
  cp = filName (pp->indexName, "/coords.sortali", "rb") ;
  if (cp)
    pp->tFileBinaryCoordsName = strnew (cp, pp->h) ;
  else
    ok = FALSE ;
  if (!ok)
    saUsage ("Some of the target binary index files are missing, please run sortalign --createIndex <indexName> ", 0, 0) ;
  

  if (1)
    {
      char *fNam = filName (pp->indexName, "/seedLength", "r") ;
      ACEIN ai = aceInCreate (fNam , 0, h) ;
      BOOL ok = FALSE ;

      if (ai && aceInCard (ai))
	{
	  char *cp = aceInWord (ai) ;
	  if (! cp || strcmp (cp, INDEXVERSION))
	    messcrash ("\nThis index is from a previous version of the progran, please run again --createIndex", cp ? cp : "NULL") ;
	  aceInStep (ai, '\t') ;
	  ok = aceInInt (ai, &(pp->seedLength)) ;
	  aceInStep (ai, '\t') ;
	  ok = aceInInt (ai, &(pp->tStep)) ;
	  aceInStep (ai, '\t') ;
	  ok = aceInInt (ai, &(pp->maxTargetRepeats)) ;
	  if (! pp->wiggle_step)
	    {
	      aceInStep (ai, '\t') ;
	      aceInInt (ai, &(pp->wiggle_step)) ;
	      if (! pp->wiggle_step)
		pp->wiggle_step = 10 ;
	    }
	}
      if (! ok)
	messcrash ("Cannot read the index file %s", fNam) ;
    }
  
  if (!pp->iStep)   /* pp->tStep is always even */
    pp->iStep = (pp->tStep > 1 ? pp->tStep / 2 : 1) ;

  ac_free (h) ;
  return NN ;
} /* saConfigCheckTargetIndex */

/*********************************************************************/
/*********************************************************************/
