/*
 * authors: Danielle and Jean Thierry-Mieg, NCBI, 
 * 15 Sept 2014
 * snpsummary.c
 *   Input: MetaDB->runs and ali
 *   Output: a wholesale SNP table with all kinds of results
 *      all : all columns
 *        r  readFate
 *        p  pairFate
 *        s  sponge
 */

/* #define MALLOC_CHECK    */
/* #define ARRAY_CHECK  */

#include "../wac/ac.h"
#include "../wh/bitset.h"

typedef struct tsnpStruct {
  const char *project ;
  const char *dbName ;
  const char *export ;
  const char *select ;
  int snpType ; /* 0: any, 1=sublib; 2=runs ; 3=group */
  int minSnpFrequency, minSnpCover, minSnpCount ;
  int allC, allM ;
  int chapter ;
  BOOL tsfOut ;
  ACEOUT ao, ao204, aoTsf ;
  KEYSET runs ;
  KEYSET groups ;
  DICT *bloomDict ;
  DICT *runDict ;
  BitSet snpDone ;
  vTXT detectionTxt2 ;
  vTXT titrationTxt5 ;
  vTXT titrationTxt3 ;
  AC_DB db ;
  AC_HANDLE h ;
  AC_TABLE snps, titrationTable, countLibsTable, groupTitles, runTitles ;
  const char *Etargets[4] ;
  const char *EtargetsBeau[4] ;
  const char *capture ;
  const char *orderBy ;
  const char *outFileName ;
  BOOL doub, histo, count_libs, titration, unique, justDetected ;
  BOOL gzi, gzo ;
  Array rrs ;
  Array histos, detectLibs ;
  BOOL Wtrue, Wfalse, DanLi, monomodal ;    /* associated to the current snp */
  KEYSET covers, mutant ;  /* associated to the current snp */
  KEYSET doubleDetect ;
  ACEOUT aoTitration, aoDT ;
  Stack  sorting_titles ;
} TSNP ;

typedef enum { F_ZERO = 0, F_Hide, F_p10, F_perCentRead, T_FMAX } T_FORMAT ;
typedef struct snpStruct {
  AC_HANDLE h ;
  AC_OBJ Snp ;
  BOOL is20_4 ;
} SNP ;

typedef struct runStruct {
  BOOL isGroup ;
  BOOL ignoreLowQ, isLowQ ;
  int sorting_title ;
  int other_title ;
  KEYSET r2g ; /* list of groups g of which r is a member */
  KEYSET g2r ; /* list of runs in this group */
} RR ;

typedef void (*TSNPFunc)(TSNP *tsnp, SNP *snp) ;
typedef struct mmStruct { char cc ; TSNPFunc f ;} MM ;
typedef struct tagtitleStruct { const char *tag, *title ; int col ; int dummy ;  T_FORMAT format ;} TT ;
static void snpBrsHistos (TSNP *tsnp, SNP *snp) ;
#define EMPTY ""

/*************************************************************************************/
/*************************** actual work *********************************************/
/*************************************************************************************/

static void snpChapterCaption (TSNP *tsnp, TT *tts, const char *caption) 
{
  TT *ti ;
  int n = 0 ;
  char buf[8128] ;
  
  memset (buf, 0, sizeof (buf)) ;

  if (caption)
    aceOutf (tsnp->ao, "\t\t%s", caption) ;

  if (tts)
    {
      for (ti = tts ; ti->tag && n < 8127 ; ti++)
	{
	  const char *cp = ti->title ;
	  
	  buf[n++] = '\t' ; /* count the paragraphs */
	  if (cp)
	    {
	      cp-- ;
	      while ((cp = strchr (cp + 1, '\t')))
		buf[n++] = '\t' ;  /* check for multi columns */
	    }
	}  
      if (n < 2) n = 2 ;
      buf[n-2] = 0 ; /* reserve 2 columns */
      aceOutf (tsnp->ao, "%s", buf) ;
    }
  return ;
} /*  snpChapterCaption */

/*************************************************************************************/
static void snpShowTag (TSNP *tsnp, SNP *snp, TT *ti)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  const char *ccp = EMPTY ;
  float z = 0 ;

  if (strcmp (ti->tag, "Spacer"))
    {
      ccp = EMPTY ; tt = ac_tag_table (snp->Snp, ti->tag, h) ;
      if (tt && tt->rows && tt->cols >= ti->col)
	{
	  ccp = ac_table_printable (tt, 0, ti->col, EMPTY) ;
	  z = ac_table_float (tt, 0, ti->col, 0) ;
	}
    }


  switch (ti->format)
    {
    default:
      aceOutf (tsnp->ao, "\t%s", ccp) ;
      break ;
    case F_p10:
      aceOutf (tsnp->ao, "\t%.1f", z/10) ;
      break ;
    case F_Hide:
      break ;
    }
} /* snpShowTag */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

static void snpVCF (TSNP *tsnp, SNP *snp)
{
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    /* { "Line", "Line Number", 0, 0, 0} , */
    { "VCF", "VCF GRCh37 Chromosome\tVCF GRCh37 position\tVCF GRCh37 ID\tVCF Reference\tVCF Variant", 1, 0, 0} ,
    { "Magic", "Magic Identifier", 1, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 

  static int nChrom = 0 ;
  static int nRefA = 0 ;
  static int nRefT = 0 ;
  static int nRefG = 0 ;
  static int nRefC = 0 ;
  static int nSub = 0 ;
  static int nDel = 0 ;
  static int nIns = 0 ;
  static int nA2G = 0 ;
  static int nT2C = 0 ;
  
  static int chapter = 0 ;
  static int line = 0 ;
  const char *caption =
    "VCF GRCh37 Identifiers"
    ;
  if (snp == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;

  AC_HANDLE h = ac_new_handle () ;
  for (ti = tts ; ti->tag ; ti++)
    {
      if (snp == 0)
	{
	  chapter = ++tsnp->chapter ;
	  aceOutf (tsnp->ao, "\t%s", ti->title) ;
	  continue ;
	}
      else if (snp == (void *)2)
	{
	  if (tsnp->aoTsf)
	    {
	      aceOutf (tsnp->aoTsf, "%s\t_%d_10__VCF_Chromosome\ti\t%d\n"
		       , ".0_100__Number"
		       , chapter
		       , nChrom
		       ) ;
	      aceOutf (tsnp->aoTsf, "%s\t_%d_13__VCF_Reference\ti\t%d\n"
		       , ".0_201__Ref_A"
		       , chapter
		       , nRefA
		       ) ;
	      aceOutf (tsnp->aoTsf, "%s\t_%d_13__VCF_Reference\ti\t%d\n"
		       , ".0_202__Ref_T"
		       , chapter
		       , nRefT
		       ) ;
	      aceOutf (tsnp->aoTsf, "%s\t_%d_13__VCF_Reference\ti\t%d\n"
		       , ".0_203__Ref_G"
		       , chapter
		       , nRefG
		       ) ;
	      aceOutf (tsnp->aoTsf, "%s\t_%d_13__VCF_Reference\ti\t%d\n"
		       , ".0_204__Ref_C"
		       , chapter
		       , nRefC
		       ) ;
	      aceOutf (tsnp->aoTsf, "%s\t_%d_14__VCF_Variant\ti\t%d\n"
		       , ".0_220__A>G"
		       , chapter
		       , nA2G
		       ) ;
	      aceOutf (tsnp->aoTsf, "%s\t_%d_14__VCF_Variant\ti\t%d\n"
		       , ".0_221__T>C"
		       , chapter
		       , nT2C
		       ) ;
	      aceOutf (tsnp->aoTsf, "%s\t_%d_14__VCF_Variant\ti\t%d\n"
		       , ".0_210__Sub"
		       , chapter
		       , nSub
		       ) ;
	      aceOutf (tsnp->aoTsf, "%s\t_%d_14__VCF_Variant\ti\t%d\n"
		       , ".0_211__Del"
		       , chapter
		       , nDel
		       ) ;
	      aceOutf (tsnp->aoTsf, "%s\t_%d_14__VCF_Variant\ti\t%d\n"
		       , ".0_212__Ins"
		       , chapter
		       , nIns
		       ) ;
	    }
	  ac_free (h) ;
	  return ;
	}
      else if (! strcmp (ti->tag, "Line"))
	{
	  aceOutf (tsnp->ao, "\t%d", ++line) ;
	}
      else if (! strcmp (ti->tag, "VCF"))
	{
	  AC_TABLE tt = ac_tag_table (snp->Snp, "VCF", h) ;
	  if (tt && tt->cols > 2)
	    aceOutf (tsnp->ao
		     , "\t%s\t%d\t.\t%s\t%s"
		     , ac_table_printable (tt, 0, 0, "-")
		     , ac_table_int (tt, 0, 1, 0)
		     , ac_table_printable (tt, 0, 2, "-")
		     , ac_table_printable (tt, 0, 3, "-")
		     ) ;
	  else
	    aceOut (tsnp->ao, "\t\t\t\t\t") ;
	  if (tsnp->aoTsf)
	    {
	      aceOutf (tsnp->aoTsf, "%s\t_%d_10__VCF_Chromosome\tt\t%s\n"
		       , ac_name (snp->Snp)
		       , chapter
		       , ac_table_printable (tt, 0, 0, "-")
		       ) ;
	      if (ac_table_printable (tt, 0, 0, 0))
		nChrom++ ;
	      aceOutf (tsnp->aoTsf, "%s\t_%d_11__VCF_position\tt\t%d\n"
		       , ac_name (snp->Snp)
		       , chapter
		       , ac_table_int (tt, 0, 1, 0)
		       ) ;
	      aceOutf (tsnp->aoTsf, "%s\t_%d_12__VCF_ID\tt\t.\n"
		       , ac_name (snp->Snp)
		       , chapter
		       ) ;
	      aceOutf (tsnp->aoTsf, "%s\t_%d_13__VCF_Reference\tt\t%s\n"
		       , ac_name (snp->Snp)
		       , chapter
		       , ac_table_printable (tt, 0, 2, "-")
		       ) ;
	      switch ((int) *ac_table_printable (tt, 0, 2, "-"))
		{
		case 'a': case 'A': nRefA++ ; break ;
		case 't': case 'T': nRefT++ ; break ;
		case 'g': case 'G': nRefG++ ; break ;
		case 'c': case 'C': nRefC++ ; break ;
		}
	      aceOutf (tsnp->aoTsf, "%s\t_%d_14__VCF_Variant\tt\t%s\n"
		       , ac_name (snp->Snp)
		       , chapter
		       , ac_table_printable (tt, 0, 3, "-")
		       ) ;
	      if (ac_has_tag (snp->Snp, "Substitution"))
		nSub++ ;
	      else if (ac_has_tag (snp->Snp, "Deletion"))
		nDel++ ;
	      else if (ac_has_tag (snp->Snp, "Insertion"))
		nIns++ ;
	      if (ac_has_tag (snp->Snp, "A2G"))
		nA2G++ ;
	      if (ac_has_tag (snp->Snp, "T2C"))
		nT2C++ ;
	    }
	}
      else if (! strcmp (ti->tag, "Magic"))
	{
	  aceOutf (tsnp->ao, "\t%s", ac_name (snp->Snp)) ;
	  
	  if (tsnp->aoTsf)
	    {
	      aceOutf (tsnp->aoTsf, "%s\t_%d_20__Magic_Identifier\tt\t%s\n"
		       , ac_name (snp->Snp)
		       , chapter
		       , ac_name (snp->Snp)
		       ) ;
	    }
	}
      else  /* gName rName pName Dan_Li Typ */
	snpShowTag (tsnp, snp, ti) ;
    }
  ac_free (h) ;
  return;
}  /* snpVCF */

/*************************************************************************************/

static void snpIdentifiers (TSNP *tsnp, SNP *snp)
{
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    {  "gName", "genome Name", 0, 0, 0} ,
    {  "rName", "RNA Name", 0, 0, 0} ,  
    {  "pName", "Protein Name", 0, 0, 0} ,  
    {  "VCF_hg38", "VCF hg38 name", 0, 0, 0} ,  
    {  "Dan_Li", "Dan Li VCF 38 name", 0, 0, 0} ,  
    {  "Typ", "Type", 0, 0, 0} ,  
    { "Coding", "Coding\tProtein type", 1, 0, 0} ,
    { "Seq_Var", "RNA variation", 10, 0, 0} ,
    { "Seq_Var", "Protein variation", 20, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 


  static int nDanLi = 0 ;
  static int nCoding = 0 ;
  const char *caption =
    "Variant Identifiers"
    ;

  static int chapter = 0 ;
  if (snp == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;

  AC_HANDLE h = ac_new_handle () ;

  for (ti = tts ; ti->tag ; ti++)
    {
      if (snp == 0)
	{
	  chapter = ++tsnp->chapter ;
	  aceOutf (tsnp->ao, "\t%s", ti->title) ;
	  continue ;
	}
      else if (snp == (void *) 2)
	{
	  aceOutf (tsnp->aoTsf, "%s\t_%d_11__Dan_Li\ti\t%d\n"
		   , ".0_300__DanLi"
		   , chapter
		   , nDanLi
		   ) ;
	  aceOutf (tsnp->aoTsf, "%s\t_%d_14__Coding\ti\t%d\n"
		   , ".0_301__Coding"
		   , chapter
		   , nCoding
		   ) ;

	  ac_free (h) ;
	  return ;
	}
      else if (! strcmp (ti->tag, "VCF"))
	{
	  AC_TABLE tt = ac_tag_table (snp->Snp, "VCF", h) ;
	  if (tt && tt->cols > 2)
	    aceOutf (tsnp->ao
		     , "\t%s\t%d\t.\t%s\t%s"
		     , ac_table_printable (tt, 0, 0, "-")
		     , ac_table_int (tt, 0, 1, 0)
		     , ac_table_printable (tt, 0, 2, "-")
		     , ac_table_printable (tt, 0, 3, "-")
		     ) ;
	  else
	    aceOut (tsnp->ao, "\t\t\t\t\t") ;
	}
      else if (! strcmp (ti->tag, "Coding"))
	{
	  const char *ccp = "", *ccq = "" ;
	  aceOut (tsnp->ao, "\t") ;
	  nCoding++ ;


	  if (ac_has_tag (snp->Snp, "Intergenic"))
	    ccp = "Intergenic" ;
	  else if (ac_has_tag (snp->Snp, "Intronic"))
	    ccp = "Intronic" ;
	  else if (ac_has_tag (snp->Snp, "Non_coding_transcript"))
	    ccp = "Non_coding_transcript" ;
	  else if (ac_has_tag (snp->Snp, "UTR_5prime"))
	    ccp = "UTR_5prime" ;
	  else if (ac_has_tag (snp->Snp, "UTR_3prime"))
	    ccp = "UTR_3prime" ;
	  else if (ac_has_tag (snp->Snp, "Synonymous"))
	    { ccp = "Coding Synonymous" ; ccq = ac_tag_text (snp->Snp, "Synonymous", 0) ; }
	  else if (ac_has_tag (snp->Snp, "AA_substitution"))
	    { ccp = "Coding substitution" ; ccq = ac_tag_text (snp->Snp, "AA_substitution", 0) ;  }
	  else if (ac_has_tag (snp->Snp, "Length_variation"))
	    {
	      AC_TABLE tt = ac_tag_table (snp->Snp, "Length_variation", h) ;
	      ccp = hprintf (h, "%s"
			     , ac_table_printable (tt, 0, 0, "")
			     ) ;
	      ccq = hprintf (h, "%s %s"
			     , ac_table_printable (tt, 0, 1, "")
			     , ac_table_printable (tt, 0, 2, "")
			     ) ;
	    }
	  
	  aceOutf (tsnp->ao, "%s\t%s", ccp, ccq) ;
	      
	  if (tsnp->aoTsf)
	    {
	      if (*ccp)
		aceOutf (tsnp->aoTsf, "%s\t_%d_14__Coding\tt\t%s\n"
			 , ac_name (snp->Snp)
			 , chapter
			 , ccp
			 ) ;
	      if (*ccq)
		aceOutf (tsnp->aoTsf, "%s\t_%d_15__Protein_type\tt\t%s\n"
			 , ac_name (snp->Snp)
			 , chapter
			 , ccq
			 ) ;
	    }
	}
      else if (! strcmp (ti->tag, "Seq_Var"))
	{
	  int ic ;
	  const char *tag1, *tag2 ;
	  AC_TABLE tt1 = 0, tt2 = 0 ;

	  switch (ti->col)
	    {
	    case 10: 
	      tag1 = "Reference_RNAexon_sequence" ;
	      tag2 = "Observed__RNAexon_sequence" ;
	      ic = 0 ;
	      break ;
	    case 20:
	      tag1 = "Reference_protein_sequence" ;
	      tag2 = "Observed__protein_sequence" ;
	      ic = 1 ;
	    }
	  tt1 = ac_tag_table (snp->Snp, tag1, h) ;
	  tt2 = ac_tag_table (snp->Snp, tag2, h) ;
	  if (tt1 && tt2 && tt1->cols > ic && tt2->cols > ic)
	    aceOutf (tsnp->ao, "\t%s > %s"
		     , ac_table_printable (tt1, 0, ic, "")
		     , ac_table_printable (tt2, 0, ic, "")
		     ) ;
	  else
	    aceOut (tsnp->ao, "\t") ;
	}
      else  /* gName rName pName Dan_Li Typ */
	{
	  nDanLi++ ;
	  snpShowTag (tsnp, snp, ti) ;
	  if (tsnp->aoTsf)
	    {
	      aceOutf (tsnp->aoTsf, "%s\t_%d_11__Dan_Li\tt\t%s\n"
		       , ac_name (snp->Snp)
		       , chapter
		       , ac_tag_printable (snp->Snp, ti->tag, "") 
		       ) ;
	    }
	}
    }
  ac_free (h) ;
  return;
}  /* snpIdentifiers */

/*************************************************************************************/

static void snpDanLiCounts (TSNP *tsnp, SNP *snp)
{
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "Compute", "Dan Li variant reads in AGLR1_A, sum of 4 libraries", 10, 0, 0} ,
    { "Compute", "Dan Li coverage in AGLR1_A, sum of 4 libraries", 11, 0, 0} ,
    { "Compute", "Dan Li variant reads in AGLR2_A, sum of 4 libraries", 20, 0, 0} ,
    { "Compute", "Dan Li coverage in AGLR2_A, sum of 4 libraries", 21, 0, 0} ,
    { "Compute", "Dan Li variant reads in ROCR1_A, sum of 4 libraries", 30, 0, 0} ,
    { "Compute", "Dan Li coverage in ROCR1_A, sum of 4 libraries",31, 0, 0} ,
    { "Compute", "Dan Li variant reads in ROCR2_A, sum of 4 libraries", 40, 0, 0} ,
    { "Compute", "Dan Li coverage in ROCR2_A, sum of 4 libraries", 41, 0, 0} ,
    /*
      { "Sum", "Dan Li variant reads in sample A, sum of 16 libraries", 0, 0, 0} ,
      { "Sum", "Dan Li coverage in sample A, sum of 16 libraries", 1, 0, 0} ,
    */
    {  0, 0, 0, 0, 0}
  } ; 
  int mm = 0, cc = 0, nn = 0 ;

  const char *caption =
    "Dan Li counts in sample A sum of 4 libraries"
    ;
  if (snp == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;

  AC_HANDLE h = ac_new_handle () ;
  for (ti = tts ; ti->tag ; ti++)
    {
      if (snp == 0)
	{
	  aceOutf (tsnp->ao, "\t%s", ti->title) ;
	  continue ;
	}
      else if (snp == (void *)2)
	{
	  ac_free (h) ;
	  return ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  const char *tag ;
	  AC_TABLE tt = 0 ;
	  char *runs[] = {"toto", "AGLR1", "AGLR2", "ROCR1", "ROCR2"} ;
	  char *run = runs[ti->col/10] ;

	  tag = "DanLi_counts" ;
	  tt = ac_tag_table (snp->Snp, tag, h) ;

	  aceOut (tsnp->ao, "\t") ;
	  if (tt)
	    {
	      int ir ;
	      for (ir = 0 ; ir < tt->rows ; ir++)
		if (!strcmp (ac_table_printable (tt, ir, 0, "toto"), run))
		  {
		    int m = 4*ac_table_int (tt, ir, 1, -1) ;
		    int c = 4*ac_table_int (tt, ir, 3, -1) ;
		    if (m >= 0 && c >= 0)
		      {
			aceOutf (tsnp->ao, "%d", ti->col %2 ? c : m) ;
			mm += m ; cc += c ;
			nn++ ;
		      }
		  }
	    }
	}
      else if (! strcmp (ti->tag, "Sum"))
	aceOutf (tsnp->ao, "\t%d", ti->col %2 ? cc : mm) ;
      else
	snpShowTag (tsnp, snp, ti) ;
    }
  ac_free (h) ;
  return;
}  /* snpDanLiCounts */

/*************************************************************************************/

static void snpDanLiFrequency (TSNP *tsnp, SNP *snp)
{
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "Compute", "Dan Li allele frequency in sample A, sum of 4 libraries AGLR1 v/v+r", 1, 0, 0} ,
    { "Compute", "Dan Li allele frequency in sample A, sum of 4 libraries AGLR2 v/v+r", 2, 0, 0} ,
    { "Compute", "Dan Li allele frequency in sample A, sum of 4 libraries ROCR1 v/v+r", 3, 0, 0} ,
    { "Compute", "Dan Li allele frequency in sample A, sum of 4 libraries ROCR2 v/v+r", 4, 0, 0} ,
    { "Sum", "Dan Li allele frequency in sample A, sum of 16 libraries v/v+r", 10, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 
  int mm = 0, cc = 0, nn = 0 ;

  const char *caption =
    "Dan Li allele frequency in sample A sum of 4 libraries"
    ;
  if (snp == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;

  AC_HANDLE h = ac_new_handle () ;
  for (ti = tts ; ti->tag ; ti++)
    {
      if (snp == 0)
	{
	  aceOutf (tsnp->ao, "\t%s", ti->title) ;
	  continue ;
	}
      else if (snp == (void *)2)
	{
	  ac_free (h) ;
	  return ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  const char *tag ;
	  AC_TABLE tt = 0 ;
	  char *runs[] = {"toto", "AGLR1", "AGLR2", "ROCR1", "ROCR2"} ;
	  char *run = runs[ti->col] ;

	  tag = "DanLi_counts" ;
	  tt = ac_tag_table (snp->Snp, tag, h) ;

	  aceOut (tsnp->ao, "\t") ;
	  if (tt)
	    {
	      int ir ;
	      for (ir = 0 ; ir < tt->rows ; ir++)
		if (!strcmp (ac_table_printable (tt, ir, 0, "toto"), run))
		  {
		    int m = ac_table_int (tt, ir, 1, -1) ;
		    int c = ac_table_int (tt, ir, 3, -1) ;
		    if  (c == 0) c = 1 ;

		    if (m >= 0 && c >= 0)
		      {
			aceOutf (tsnp->ao, "%.2f", 100.0 * m/c) ;
			mm += m ; cc += c ; nn ++ ;
		      }
		  }
	    }
	}
      else if (! strcmp (ti->tag, "Sum"))
	{
	  if  (cc == 0) cc = 1 ;
	  if (mm >= 0 && cc >= 0)
	    aceOutf (tsnp->ao, "\t%.2f", 100.0 * mm/cc) ;
	  else
	    aceOut (tsnp->ao, "\t") ;
	}
      else
	snpShowTag (tsnp, snp, ti) ;
    }
  ac_free (h) ;
  return;
}  /* snpDanLiFrequency */

/*************************************************************************************/

static void snpBrsTitrationDoOne (TSNP *tsnp, SNP *snp, AC_TABLE tt, int ir0, int ir1)
{
  int NN = ir1 - ir0 ;
  if (NN == 3 || NN == 5)
    {
      AC_HANDLE h = ac_new_handle () ;
      int ir, ii ;
      BOOL ok = TRUE ;
      int cov[5], mut[5], ff[5] ;

      for (ii = 0,  ir = ir0 ; ok && ir < ir1 ; ii++, ir++)
	{
	  KEY c = ac_table_key (tt, ir, 0, 0) ;
	  int n = ac_table_int (tt, ir, 1, -1) ;
	  int r = 0 ;
	  
	  if (! dictFind (tsnp->runDict, ac_table_printable (tt, ir, 2, 0), &r))
	    { ok = FALSE ; continue ; }
	  
	  if (!c || !r || n != ii + 1)
	    { ok = FALSE ; continue ; }
	  cov[ii] = keySet (tsnp->covers, r) ;
	  mut[ii] = keySet (tsnp->mutant, r) ;
	  if (cov[ii] < 20)
	    ok = FALSE ; 
	  ff[ii] = .49 + 10000.0 * mut[ii]/cov[ii] ;
	}
      if (ok)
	{  /* correctly found 3 or 5  values, check order */
	  BOOL ok3 = FALSE ;
	  BOOL ok5 = FALSE ;
	  if (NN == 5 && ff[0] > ff[4] + 200 && ff[0] >= 500 && ff[2] > 200)
	    {
	      ok3 = FALSE ; ok5 = TRUE ;
	      for (ii = 1; ii < 5 ; ii++)
		if (ff[ii] > ff[ii-1])
		  ok5 = FALSE ;
	    }
	  else if (NN == 5 && ff[4] > ff[0] + 200 && ff[4] >= 500 && ff[2] > 200)
	    {
	      ok3 = FALSE ; ok5 = TRUE ;
	      for (ii = 1; ii < 5 ; ii++)
		if (ff[ii] < ff[ii-1])
		  ok5 = FALSE ;
	    }
	  else if (NN == 3 && ff[0] > ff[4] - 200 && ff[0] >= 500 && ff[2] > 200)
	    {
	      ok3 = TRUE ; ok5 = FALSE ;
	      for (ii = 2 ; ii < 5 ; ii+=2)
		if (ff[ii] > ff[ii-2])
		  ok3 = FALSE ;
	    }
	  else if (NN == 3 && ff[4] > ff[0] - 200 && ff[4] >= 500 && ff[2] > 200)
	    {
	      ok3 = TRUE ; ok5 = FALSE ;
	      for (ii = 2 ; ii < 5 ; ii+=2)
		if (ff[ii] < ff[ii-2])
		  ok3 = FALSE ;
	    }
	  if (ok5)
	    {
	      const char *ccp = ac_table_printable (tt, ir0, 3, 0)  ;
	      if (! ccp) ccp = ac_table_printable (tt, ir0, 0, "toto")  ;
	      aceOutf (tsnp->aoTitration, "%s\t%s__5\t5f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n"
			  , ac_name (snp->Snp) 
			  , ccp
			  , ff[0]/100.0
			  , ff[1]/100.0
			  , ff[2]/100.0
			  , ff[3]/100.0
			  , ff[4]/100.0
			  ) ;			  
	      vtxtPrintf (tsnp->titrationTxt5, "%s,", ccp) ;
	    }
	  if (ok3)
	    {
	      const char *ccp = ac_table_printable (tt, ir0, 3, 0)  ;
	      if (! ccp) ccp = ac_table_printable (tt, ir0, 0, "toto")  ;

	      aceOutf (tsnp->aoTitration, "%s\t%s__3\t5f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n"
		       , ac_name (snp->Snp) 
		       , ccp
		       , ff[0]/100.0
		       , -10
		       , ff[2]/100.0
		       , -10
		       , ff[4]/100.0
			   ) ;			  
	      vtxtPrintf (tsnp->titrationTxt3, "%s,", ccp) ;
	    }
	}

      ac_free (h) ;
    }
  return;
} /*snpBrsTitrationDoOne */

/*************************************************************************************/

static void snpBrsTitrationDo (TSNP *tsnp, SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tt = tsnp->titrationTable ;
  int ir, irMax = tt ? tt->rows : 0, irOld = -1 ;
  KEY oldC = 0 ;
  
  for (ir = 0 ; ir < irMax ; ir++)
    {
      KEY c = ac_table_key (tt, ir, 0, 0) ;
      if (c == oldC)
	continue ;
      oldC = c ;
      if (irOld >= 0)
	snpBrsTitrationDoOne (tsnp, snp, tt, irOld, ir) ;
      irOld= ir ;
    }
  
  if (ir > irOld + 1)
    snpBrsTitrationDoOne (tsnp, snp, tt, irOld, ir) ;
  
  ac_free (h) ;
} /* snpBrsTitrationDo */

/*************************************************************************************/
/* in each group count how many libs agree of disagree with the average */

static BOOL snpBrsCountLibsDoOne (TSNP *tsnp, SNP *snp, int gg1, int gg2, int gg3)
{
  int cc = tsnp->allC ;   
  int mm = tsnp->allM ;   
  float ff = 100.0 * mm/cc ;
  ACEOUT aoDT = tsnp->aoDT ;
  int type = 0 ;
  Array aaDetect = tsnp->detectLibs ;
  BOOL detected = FALSE ;
  
  if (ff >= 2) type = 1 ;
  if (ff >= 5) type = 2 ;
  if (ff >= 20) type = 3 ;
  if (ff >= 80) type = 4 ;
  if (ff >= 95) type = 5 ;
  
  if (gg1 == gg2)
    {
      int gg = gg1 ;
      RR *rr = arrayp (tsnp->rrs, gg, RR) ;
      KEYSET g2r = rr->g2r ;
      if (g2r)
	{
	  int c = keySet (tsnp->covers, gg) ;
	  int m = keySet (tsnp->mutant, gg) ;
	  float f, df ;
	  int nLib = 0 ;

	  if (c >= 20)
	    {
	      f = 100.0 * m/c ;
	      df = 30 ;
	      if (c < 44)
		df = 200 / sqrt (c) ;
	      
	      for (int j = 0 ; j < keySetMax (g2r) && nLib < 2 ; j++)
		{
		  int r = keySet (g2r, j) ;
		  int c2 = keySet (tsnp->covers, r) ;
		  int m2 = keySet (tsnp->mutant, r) ;
		  if (m2 >= 4 && c2 >= 20 && 100 * m2 >= 2 * c2)
		    nLib++ ;		  
		}
	      if (nLib >= 2)
		{
		  detected = TRUE ;
		  vtxtPrintf (tsnp->detectionTxt2, "%s,", stackText (tsnp->sorting_titles, rr->other_title)) ;
		  for (int type2 = type ; type2 >= 0 ; type2--)
		    {
		      aceOutf (aoDT, "%s\t%s\tftttt\t%.2f", ac_name (snp->Snp), dictName (tsnp->runDict, gg1), f) ;
		      if (f - df < ff && f + df > ff)
			{
			  array (aaDetect, 24 * gg + 6 * 0 + type2, int)++ ;
			  aceOut (aoDT, "\tCompatible") ;
			}
		      else
			{
			  array (aaDetect, 24 * gg + 6 * 1 + type2, int)++ ;
			  aceOut (aoDT, "\tContradition") ;
			}
		      if (tsnp->Wtrue)
			{
			  array (aaDetect, 24 * gg + 6 * 2 + type2, int)++ ;
			  aceOut (aoDT, "\tWtrue") ;
			}
		      if (tsnp->Wfalse)
			{
			  array (aaDetect, 24 * gg + 6 * 3 + type2, int)++ ;
			  aceOut (aoDT, "\tWfalse") ;
			}
		      if (tsnp->DanLi)
			aceOut (aoDT, "\tDanLi") ;
		      aceOut (aoDT, "\n") ;
		    }
		} 
	    }
	}
    }
  return detected ;
} /* snpBrsCountLibsDoOne */

/*************************************************************************************/

static BOOL snpBrsCountLibsDo (TSNP *tsnp, SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tt = tsnp->countLibsTable ;
  int irMax = tt ? tt->rows : 0 ;
  KEYSET done = keySetHandleCreate (h) ;
  int N = dictMax (tsnp->runDict) + 1 ;
  BOOL detected = FALSE ;

  for (int ir = 0 ; ir < irMax ; ir++)
    {
      KEY c1 = ac_table_key (tt, ir, 0, 0) ;
      int n1 = ac_table_int (tt, ir, 1, 0) ;
      int r1 = 0 ;
      if (! dictFind (tsnp->runDict, ac_table_printable (tt, ir, 2, 0), &r1))
	continue ;

      for (int jr = ir ; jr < irMax ; jr++)
	{
	  KEY c2 = ac_table_key (tt, jr, 0, 0) ;
	  int n2 = ac_table_int (tt, jr, 1, 0) ;
	  int r2 = 0, r3 ;
	  if (! dictFind (tsnp->runDict, ac_table_printable (tt, jr, 2, 0), &r2))
	    continue ;
	  r3 = (r1 < r2 ? N * r1 + r2 : N * r2 + r1 ) ;
	  if (keySet (done, r3) == 1)
	    continue ;
	  keySet (done, r3) = 1 ;
	  if (c1 != c2 || (n1/100) != (n2/100))
	    break ;
	  detected |= snpBrsCountLibsDoOne (tsnp, snp, r1, r2, r3) ;
	}
    }

    ac_free (h) ;
    return detected ;
} /* snpBrsCountLibsDo */

/*************************************************************************************/

static BOOL snpBrsParseCounts (TSNP *tsnp, SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
  BOOL detected = FALSE ;
  int minF = tsnp->minSnpFrequency ;
  int minM = tsnp->minSnpCount ;
  int minC = tsnp->minSnpCover > 0 ? tsnp->minSnpCover : 1 ;
  snp->is20_4 = FALSE ;

  if (minM <= 0) minM = 4 ;
  for (int ir = 0 ; ir < keySetMax (tsnp->runs) ; ir++)
    {
      int r = keySet (tsnp->runs, ir) ;
      RR *rr = arrayp (tsnp->rrs, r, RR) ;
      KEYSET r2g = rr->r2g ;
      int c, m ;
      
      if (rr->isGroup)  /* avoid double counting */
	continue ;
      c = keySet (tsnp->covers, r) ;
      m = keySet (tsnp->mutant, r) ;
      if (100 * m >= minF * c && c >= minC && m >= minM)
	snp->is20_4 = TRUE ;
      if (r2g && keySetMax (r2g))
	{
	  for (int i = 0 ; i < keySetMax (r2g) ; i++)
	    {
	      int g = keySet (r2g, i) ;
	      RR *rrg = arrayp (tsnp->rrs, g, RR) ;
	      if (! rr->isLowQ || ! rrg->ignoreLowQ)
		{
		  keySet (tsnp->covers, g) += c ;
		  keySet (tsnp->mutant, g) += m ;
		}
	    }
	}	      
      }

  for (int ir = 0 ; ir < keySetMax (tsnp->runs) ; ir++)
    {
      int r = keySet (tsnp->runs, ir) ;
      RR *rr = arrayp (tsnp->rrs, r, RR) ;
      if (rr->isGroup)
	{
	  int c = keySet (tsnp->covers, r) ;
	  int m = keySet (tsnp->mutant, r) ;
	  
	  if (100 * m >= minF * c && c >= minC && m >= 4)
	    snp->is20_4 = TRUE ;
	}
    }
  if (snp->is20_4)
    {  
      const char *cp, *cq ;

      for (int ir = 0 ; ir < keySetMax (tsnp->runs) ; ir++)
	{
	  int r = keySet (tsnp->runs, ir) ;
	  RR *rr = arrayp (tsnp->rrs, r, RR) ;
	  if (! rr->isGroup)  /* avoid double counting */
	    {
	      int c = keySet (tsnp->covers, r) ;
	      int m = keySet (tsnp->mutant, r) ;
	      if (! rr->isLowQ)
		{
		  tsnp->allC += c ;
		  tsnp->allM += m ;
		}
	    }
	}
      if (tsnp->allC >= 20 && tsnp->count_libs)
	detected = snpBrsCountLibsDo (tsnp, snp) ;
      if (tsnp->justDetected && ! detected)
	return FALSE ;
      cp = "X" ;
      if (ac_has_tag (snp->Snp, "Wtrue")) cp = "Wtrue" ;
      else if (ac_has_tag (snp->Snp, "Wtrue2")) cp = "Wtrue2" ;
      else if (ac_has_tag (snp->Snp, "Wfalse")) cp = "Wfalse" ;
      else if (ac_has_tag (snp->Snp, "Wfalse2")) cp = "Wfalse2" ;

      cq = ac_tag_printable (snp->Snp, "Monomodal", "toto") ;
      if (! cq)  cq = "X" ;
      if (tsnp->ao204)
	aceOutf (tsnp->ao204, "%s\t%s\tift\t%d\t%.2f\t%s\t%s\n", ac_name (snp->Snp), tsnp->project, nn, 100.0 * tsnp->allM/tsnp->allC, cp, cq) ;
      if (tsnp->histo)
	snpBrsHistos (tsnp, snp) ;
      if (tsnp->titration)
	snpBrsTitrationDo (tsnp, snp) ;
    }
  ac_free (h) ;
  return TRUE ;
} /*  snpBrsParseCounts */

/*************************************************************************************/
/* parse the actual count of a snp, and cumulate in tsnp->cover/count for equivalent snp */

static int snpBrsParseTrueCounts (TSNP *tsnp, SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tt = 0 ;
  const char *tag = "BRS_counts" ;
  int r, allC = 0 ;

  tt = ac_tag_table (snp->Snp, tag, h) ;
  for (int ir = 0 ; tt && ir < tt->rows ; ir++)
    if (dictFind (tsnp->runDict, ac_table_printable (tt, ir, 0, "toto"), &r))
      {
	RR *rr = arrayp (tsnp->rrs, r, RR) ;
	
	if (! rr->isGroup)  /* avoid double counting */
	  {
	    int c = ac_table_int (tt, ir, 1, 0) ;
	    int m = ac_table_int (tt, ir, 2, 0) ;
	    
	    keySet (tsnp->covers, r) += c ;
	    keySet (tsnp->mutant, r) += m ;
	    allC += c ;
	    /*
	      if (! strcmp (ac_name (snp->Snp), "NIN.aAug10:4892:Del:tA:t") && ! strncmp ( ac_table_printable (tt, ir, 0, "toto"), "SRR", 3))
	      invokeDebugger () ;
	    */
	  }
      }

  ac_free (h) ;
  return allC ;
}  /* snpBrsParseTrueCounts */

/*************************************************************************************/
/* populate snp0 with the netx [unique] snp */
#define KK(kk)  ((KEY)( ((KEY) (kk)) & ((KEY) 0xffffffL) ))
static BOOL snpBrsNextSnp (TSNP *tsnp, SNP *snp, int iSnp)
{
  /* clean up */
  ac_free (snp->h) ; 
  snp->h = ac_new_handle () ;
  snp->Snp = 0 ;
  tsnp->covers = arrayHandleCreate (dictMax (tsnp->runDict) + 1, KEY, snp->h) ;
  tsnp->mutant = arrayHandleCreate (dictMax (tsnp->runDict) + 1, KEY, snp->h) ;
  tsnp->allC = tsnp->allM = 0 ;
  tsnp->Wtrue= tsnp->Wfalse = tsnp->DanLi = tsnp->monomodal = FALSE ;
  vtxtClear (tsnp->detectionTxt2) ;
  vtxtClear (tsnp->titrationTxt5) ;
  vtxtClear (tsnp->titrationTxt3) ;
  
  if (! tsnp->unique)
    {
      snp->Snp = ac_table_obj (tsnp->snps, iSnp, 0, snp->h) ;
      if (snp->Snp)
	{
	  snpBrsParseTrueCounts (tsnp, snp) ;
	  tsnp->Wtrue  = ac_has_tag (snp->Snp, "Wtrue") ||  ac_has_tag (snp->Snp, "Wtrue2") ;
	  tsnp->Wfalse = ac_has_tag (snp->Snp, "Wfalse") ||  ac_has_tag (snp->Snp, "Wfalse2") ;
	  tsnp->DanLi = ac_has_tag (snp->Snp, "Dan_Li") ;
	  if (! snpBrsParseCounts (tsnp, snp))
	    { ac_free (snp->h) ; return FALSE ; }
	}
    }
  else
    {
      KEY key = ac_table_key (tsnp->snps, iSnp, 0, 0) ;

      if (! bitt (tsnp->snpDone, KK(key)))
	{
	  AC_HANDLE h = ac_new_handle () ;
	  AC_OBJ Snp = ac_table_obj (tsnp->snps, iSnp, 0, h) ;
	  bitSet (tsnp->snpDone, KK(key)) ;
	  
	  if (Snp)
	    {      
	      const char *errors = 0 ;
	      AC_TABLE vv = 0 ;
	      AC_KEYSET aks = ac_new_keyset (tsnp->db, h) ;
	      ac_keyset_add (aks, Snp) ;
	      char *qq = "select t,x1,x2 from s in @, c1 in s->vcf, x1 in c1[1], a1 in c1[2], b1 in c1[3] , g in s->gene, t in g->variant, c2 in t->vcf, x2 in c2[1], a2 in c2[2], b2 in c2[3] where c1 == c2 && x1 - x2 == 0 && a1 == a2 && b1 == b2" ;
	      
	      vv= ac_bql_table (tsnp->db
				, qq
				, aks
				, 0
				, &errors
				, tsnp->h
				) ;
	      if (vv)
		{ /* cumulate the counts in tsnp->covers and select best representative */
		  int bestC = -1, bestIv = -1 ;
		  for (int iv = 0 ; iv < vv->rows ; iv++)
		    {
		      snp->Snp = ac_table_obj (vv, iv, 0, h) ;
		      int c = snpBrsParseTrueCounts (tsnp, snp) ;
		      if (c > bestC) 
			{ bestC = c ; bestIv = iv ;	}
		      key = ac_obj_key (snp->Snp) ;
		      tsnp->monomodal |= ac_has_tag (snp->Snp, "Monomodal") ;
		      bitSet (tsnp->snpDone, KK(key)) ;
		    }
		  snp->Snp = 0 ;
		  if (bestIv >= 0)
		    {
		      snp->Snp = ac_table_obj (vv, bestIv, 0, snp->h) ; 
		      tsnp->Wtrue  = ac_has_tag (snp->Snp, "Wtrue") ||  ac_has_tag (snp->Snp, "Wtrue2") ;
		      tsnp->Wfalse = ac_has_tag (snp->Snp, "Wfalse") ||  ac_has_tag (snp->Snp, "Wfalse2") ;
		      tsnp->monomodal |= ac_has_tag (snp->Snp, "Monomodal") ;
		      if (! snpBrsParseCounts (tsnp, snp))
			{ ac_free (h) ; ac_free (snp->h) ; return FALSE ; }
		    }

		}
	    }
	  ac_free (h) ;
	}
    }
  if (snp->Snp)
    return TRUE ;
  return FALSE ;
}  /* snpBrsNextSnp */

/*************************************************************************************/

static void snpBrsFrequencyCounts (TSNP *tsnp, SNP *snp, int type)
{
  BOOL isFrequency ;
  BOOL isGroup ;
  const char *caption ;
  const char *subcaption, *subcaption2 ;
  KEYSET rg ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "Groups", "", 10, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 

  static KEYSET nSub, nIns, nDel, nA2G, nT2C ;

  static int chapter = 0 ;

  switch (type)
    {
    case 0: 
      isFrequency = TRUE ; isGroup = FALSE ; 
      rg = tsnp->runs ;
      caption = "Allele frequency in individual libraries if more than 20 covering reads, otherwise -20" ;
      subcaption = "" ;
      break ;
    case 1: 
      isFrequency = FALSE ; isGroup = FALSE ; 
      rg = tsnp->runs ;
      caption = "Number of variant reads and total SNV site coverage in individual libraries" ;
      subcaption = "Variant reads in" ;
      subcaption2 = "coverage in" ;
      break ;
    case 2: 
      isFrequency = TRUE ; isGroup = TRUE ; 
      rg = tsnp->groups ;
      caption = "Allele Frequency in groups if more than 20 covering reads, otherwise -20" ;
      subcaption = "Allele Frequency in" ;
      break ;
    case 3: 
      isFrequency = FALSE ; isGroup = TRUE ; 
      rg = tsnp->groups ;
      caption = "Number of variant reads and total SNV site coverage in groups" ;
      subcaption = "Variant reads in" ;
      subcaption2 = "coverage in" ;
      break ;
    default: 
      messcrash ("Bad call to snpBrsFrequencyCounts type=%d", type) ;
      break ;
    }
  if (snp == (void *) 1)
    {
      aceOutf (tsnp->ao, "\t\t%s", caption) ;
      if (! isFrequency) 
	aceOut (tsnp->ao, "\t") ;
      for (int ir = 1 ; ir < keySetMax (rg); ir++)
       	aceOut (tsnp->ao,  isFrequency ? "\t" : "\t\t") ;

      return  ;
    }

  AC_HANDLE h = ac_new_handle () ;
  for (ti = tts ; ti->tag ; ti++)
    {
      if (snp == 0)
	{
	  chapter = ++tsnp->chapter ;
	  if (ti->col == 0)
	    aceOut (tsnp->ao, "\t") ;
	  else
	    {
	      for (int ir = 0 ; ir < keySetMax (rg); ir++)
		{
		  int r = keySet (rg, ir) ;
		  aceOutf (tsnp->ao, "\t%s %s", subcaption, dictName (tsnp->runDict, r)) ;
		  if (! isFrequency)
		    aceOutf (tsnp->ao, "\t%s %s", subcaption2, dictName (tsnp->runDict, r)) ;
		}
	    }
	  nSub = keySetHandleCreate (tsnp->h) ;
	  nDel = keySetHandleCreate (tsnp->h) ;
	  nIns = keySetHandleCreate (tsnp->h) ;
	  nA2G = keySetHandleCreate (tsnp->h) ;
	  nT2C = keySetHandleCreate (tsnp->h) ;
	  
	  continue ;
	}
      else if (snp == (void *)2)
	{
	  if (tsnp->aoTsf && isGroup && isFrequency)
	    {
	      for (int ir = 0 ; ir < keySetMax (rg); ir++)
		{
		  int r = keySet (rg, ir) ;
		  
		  aceOutf (tsnp->aoTsf, "%s\t_%d__%s\ti\t%d\n"
			   , ".0_220__A>G"
			   , chapter
			   , dictName (tsnp->runDict, r)
			   , keySet (nA2G, r) 
			   ) ;
		  aceOutf (tsnp->aoTsf, "%s\t_%d__%s\ti\t%d\n"
			   , ".0_221__T>C"
			   , chapter
			   , dictName (tsnp->runDict, r)
			   , keySet (nT2C, r) 
			   ) ;
		  aceOutf (tsnp->aoTsf, "%s\t_%d__%s\ti\t%d\n"
			   , ".0_210__Sub"
			   , chapter
			   , dictName (tsnp->runDict, r)
			   , keySet (nSub, r) 
			   ) ;
		  aceOutf (tsnp->aoTsf, "%s\t_%d__%s\ti\t%d\n"
			   , ".0_211__Del"
			   , chapter
			   , dictName (tsnp->runDict, r)
			   , keySet (nDel, r) 
			   ) ;
		  aceOutf (tsnp->aoTsf, "%s\t_%d__%s\ti\t%d\n"
			   , ".0_212__Ins"
			   , chapter
			   , dictName (tsnp->runDict, r)
			   , keySet (nIns, r) 
			   ) ;
		}
	    }
	  ac_free (h) ;
	  return ;
	}
      if (ti->col == 10 && rg)
	{
	  for (int ir = 0 ; ir < keySetMax (rg); ir++)
	    {
	      int r = keySet (rg, ir) ;
	      int c = keySet (tsnp->covers, r) ;
	      int m = keySet (tsnp->mutant, r) ;
	      
	      if (!r)
		continue ;
	      if (isFrequency)
		{
		  if (c >= 20)
		    aceOutf (tsnp->ao, "\t%.2f", 100.0*m/c) ;
		  else
		    aceOutf (tsnp->ao, "\t-20") ;
		}
	      else
		aceOutf (tsnp->ao, "\t%d\t%d", m, c) ;
	      if (c >= 10 && 100 * m >= 2 * c && m >= 4)
		{
		  if (ac_has_tag (snp->Snp, "Substitution"))
		    keySet (nSub, r)++ ;
		  else if (ac_has_tag (snp->Snp, "Deletion"))
		    keySet (nDel, r)++ ;
		  else if (ac_has_tag (snp->Snp, "Insertion"))
		    keySet (nIns, r)++ ;
		  if (ac_has_tag (snp->Snp, "A2G"))
		    keySet (nA2G, r)++ ;
		  if (ac_has_tag (snp->Snp, "T2C"))
		    keySet (nT2C, r)++ ;
		}
	    }
	}
      else
	snpShowTag (tsnp, snp, ti) ;
    }
  ac_free (h) ;
  return;
}  /* snpBrsFrequencyCounts */

/*************************************************************************************/

static void snpBrsHistos (TSNP *tsnp, SNP *snp)
{
  KEYSET rg = tsnp->groups ;
  
  if (rg)
    {
      for (int ir = 0 ; ir < keySetMax (rg); ir++)
	{
	  int gg = keySet (rg, ir) ;
	  int cc = keySet (tsnp->covers, gg) ;
	  int mm = keySet (tsnp->mutant, gg) ;
	  float ff = -10 ;
	  Array histo = 0 ;
	  
	  histo = array (tsnp->histos, gg, Array) ;
	  if (! histo)
	    histo = array (tsnp->histos, gg, Array) = arrayHandleCreate (110, int, tsnp->h) ;
	  
	  if (cc >= 20)
	    {
	      ff = 100.0 * mm/cc ;
	      array (histo, ff + .49, int)++ ;
	    } 
	}
    }

  return ;
}  /* snpBrsHistos */

/*************************************************************************************/

static void snpBrsGroup2345 (TSNP *tsnp, SNP *snp)
{
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "Groups", "", 1, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 

  const char *caption = "Number of libraries with coverage >= 20 allowing SNV calling\tCompatibility: number of compatible libraries with coverage >=20 and Delta(AF) < 10%\tInconcistency across libraries: number of libraries with coverage >= 20 and Delta(AF) > 30%\tNumber of libraries with insufficient coverage < 20" ;
    ;
  const char *subcaptions[] =
    { "Libraries with coverage >= 20 in", "Compatible librarie in", "Incompatible libraries in",  "Insufficient coverage < 20 in"} ;
    ;

  if (snp == (void *) 1)
    {
      KEYSET rg = tsnp->groups ;

      aceOutf (tsnp->ao, "\t\t%s", caption) ;
      for (int ir = 1 ; ir < 4*keySetMax (rg); ir++)
	aceOut (tsnp->ao, "\t") ;
      return  ;
    }

  AC_HANDLE h = ac_new_handle () ;
  for (ti = tts ; ti->tag ; ti++)
    {
      KEYSET rg = tsnp->groups ;
      
      if (snp == 0)
	{
	  if (ti->col == 0)
	    aceOut (tsnp->ao, "\t") ;
	  else
	    {
	      for (int ir = 0 ; ir < keySetMax (rg); ir++)
		{
		  int r = keySet (rg, ir) ;
		  for (int type = 0 ; type < 4 ; type++)
		    aceOutf (tsnp->ao, "\t%s %s", subcaptions[type], dictName (tsnp->runDict, r)) ;
		}
	    }
	  continue ;
	}
      else if (snp == (void *)2)
	{
	  ac_free (h) ;
	  return ;
	}
      else if (! strcmp (ti->tag, "Groups"))
	{
	  rg = tsnp->groups ;
	  
	  if (rg)
	    {
	      for (int ir = 0 ; ir < keySetMax (rg); ir++)
		{
		  BOOL ok = TRUE ;
		  
		  int gg = keySet (rg, ir) ;
		  int cc = keySet (tsnp->covers, gg) ;
		  int mm = keySet (tsnp->mutant, gg) ;
		  float ff = -10 ;
		  RR *rr = arrayp (tsnp->rrs, gg, RR) ;
		  int nLow = 0, nHigh = 0, nOk = 0, nContra = 0 ;
		  Array histo = 0 ;

		  if (tsnp->histo)
		    {
		      histo = array (tsnp->histos, gg, Array) ;
		      if (! histo)
			histo = array (tsnp->histos, gg, Array) = arrayHandleCreate (110, int, tsnp->h) ;
		    }


		  if (! rr->g2r)
		    ok = FALSE ;
		  if (cc >= 20)
		    {
		      ff = 100.0 * mm/cc ;
		      if (histo)
			array (histo, ff + .49, int)++ ;
		    } 
		  else 
		    ok = FALSE ;
		  
		  
		  for (int ii = 0 ; ok && ii < keySetMax (rr->g2r) ; ii++)
		    {
		      int r = keySet (rr->g2r, ii) ;
		      int c = keySet (tsnp->covers, r) ;
		      int m = keySet (tsnp->mutant, r) ;
		      
		      if (!r)
			continue ;
		      
		      if (c < 20)
			nLow++ ;
		      else if (c >= 20)
			{
			  float f = 100 * m / c ;
			  nHigh ++ ;
			  /* float df = 200 * sqrt(c)/c ; */
			  if (f > ff - 10 && f < ff+10)
			    nOk++ ;
			  else if (f > ff + 30 || f < ff -30)
			    nContra++ ;
			}
		    }
		  if (ok)
		    {
		      aceOutf (tsnp->ao, "\t%d", nHigh) ;
		      aceOutf (tsnp->ao, "\t%d", nOk) ;
		      aceOutf (tsnp->ao, "\t%d", nContra) ;
		      aceOutf (tsnp->ao, "\t%d", nLow) ;
		    }
		  else
		    aceOut (tsnp->ao, "\t\t\t\t") ;
		}
	    }
	}
      else
	snpShowTag (tsnp, snp, ti) ;
    }
  ac_free (h) ;
  return;
}  /* snpBrsGroup2345set toto=RESULTS/SNV/SnpA2R2.DanLi.SNP_summary.june9.txt */


static void snpBrsRunFrequency (TSNP *tsnp, SNP *snp)
{
  return snpBrsFrequencyCounts (tsnp, snp, 0) ;
}
static void snpBrsRunCounts (TSNP *tsnp, SNP *snp)
{
  return snpBrsFrequencyCounts (tsnp, snp, 1) ;
}
static void snpBrsGroupFrequency (TSNP *tsnp, SNP *snp)
{
  return snpBrsFrequencyCounts (tsnp, snp, 2) ;
}
static void snpBrsGroupCounts (TSNP *tsnp, SNP *snp)
{
  return snpBrsFrequencyCounts (tsnp, snp, 3) ;
}


/*************************************************************************************/

static void snpBrsQC(TSNP *tsnp, SNP *snp)
{
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "Capture", "Targeted by", 0, 0, 0} ,  
    { "Monomodal", "Monomodal", 0, 0, 0} ,  
    { "Wendell", "Wendell Jones trusted site genomic in sample A", 14, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 

  const char *caption =
    "SNP quality control"
    ;

  if (snp == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;

  AC_HANDLE h = ac_new_handle () ;
  for (ti = tts ; ti->tag ; ti++)
    {
      if (snp == 0)
	{
	  aceOutf (tsnp->ao, "\t%s", ti->title) ;
	  continue ;
	}
      else if (snp == (void *)2)
	{
	  ac_free (h) ;
	  return ;
	}
      else if (! strcmp (ti->tag, "Wendell"))
	{
	  AC_TABLE tt = 0 ;

	  tt = ac_tag_table (snp->Snp, "Wtrue", h) ;
	  if (!  tt) 
	    tt = ac_tag_table (snp->Snp, "Wtrue2", h) ;
	  if (tt && tt->cols >=1)
	    aceOut (tsnp->ao, "\tVar") ;
	  else if (ac_has_tag (snp->Snp, "Wfalse") || ac_has_tag (snp->Snp, "Wfalse2"))
	    aceOut (tsnp->ao, "\tRef") ;
	  else
	    aceOut (tsnp->ao, "\t") ;
	}
      else if (! strcmp (ti->tag, "Capture"))
	{
	  AC_TABLE tt = 0 ;
	  AC_OBJ Gene = ac_tag_obj (snp->Snp, "Gene", h) ;
	  char *sep = "" ;
	  const char *errors = 0 ;
	  char *qq = Gene ? hprintf (h, "select c from g in ?Gene where g == \"%s\", c in g->Capture", ac_name(Gene)) : 0 ;

	  aceOut (tsnp->ao, "\t") ;
	  tt = Gene ? ac_obj_bql_table (Gene, qq, "+1", &errors, h) : 0 ;
	  for (int ir = 0 ; tt && ir < tt->rows ; ir++)
	    {
	      const char *ccp = ac_table_printable (tt, ir, 0, 0) ;
	      if (ccp && ccp[2]==0)
		{
		  aceOutf (tsnp->ao, "%s%s", sep, ccp) ;
		  sep = "," ;
		}
	    }
	}
      else if (! strcmp (ti->tag, "Monomodal"))
	{
	  AC_TABLE tt = ac_tag_table (snp->Snp, "Monomodal", h) ;
	  int ir ;
	  char sep = '\t' ;
	  
	  if (tt && ! ac_has_tag  (snp->Snp, "Manual_non_monomodal"))
	    {
	      BOOL ok = FALSE ;
	      for (ir = 0 ; ir < tt->rows ; ir++)
		{
		  const char *cp = ac_table_printable (tt, ir, 0, 0) ;
		  if (cp)
		    {
		      aceOutf (tsnp->ao
			       , "%c%s"
			       , sep
			       , cp
			       ) ;
		      sep = ',' ;
		      ok = TRUE ;
		    }
		}
	      if (! ok && tsnp->monomodal)
		{
		  aceOut (tsnp->ao, "\tAlias") ;
		  ok = TRUE ;
		}
	      if (! ok)
		aceOut (tsnp->ao, "\t") ;
	    }
	  else
	    aceOut (tsnp->ao, "\t") ;
	}
      else
	snpShowTag (tsnp, snp, ti) ;
    }
  ac_free (h) ;
  return;
} /* snpBrsQC */

/*************************************************************************************/

static void snpBrsAllRuns (TSNP *tsnp, SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;

  char *sumTit = hprintf (h, "Variant count in project %s\tCoverage count in project %s", tsnp->project, tsnp->project) ;
  char *avTit = hprintf (h, "Average allele frequency in project %s", tsnp->project) ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "Sum", sumTit, 10, 0, 0} ,
    { "Detection", "Detected in 2 libraries coverage >= 20, variant >= 4, AF >= 2%", 40, 0, 0} ,
    { "Nlib", "Measurable in n libraries\tCalled in n libraries\tNot Measured\tAllele frequency contradictions", 20, 0, 0} ,
    { "Freq", avTit, 30, 0, 0} ,
    { "Wendell", "Genomic allele Frequency in sample A (Jones 2021)", 50, 0, 0} ,
    { "RNA_A", "RNA allele Frequency in sample A", 60, 0, 0} ,
    { "RNA_B", "RNA allele Frequency in sample B", 70, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 

  const char *caption =
    "Allele counts and frequency when summing all runs in this table"
    ;

  float ff = -10 ;

  if (snp == (void *) 1)
    {
      ac_free (h) ;
      return  snpChapterCaption (tsnp, tts, caption) ;
    }

  for (ti = tts ; ti->tag ; ti++)
    {
      if (snp == 0)
	{
	  aceOutf (tsnp->ao, "\t%s", ti->title) ;
	  continue ;
	}
      else if (snp == (void *)2)
	{
      	  ac_free (h) ;
	  return ;
	}
      if (ti->col == 10)
	{	 
	  int c = tsnp->allC ;
	  int m = tsnp->allM ;
	  if (c >= 10)
	    ff = 100.0*m/c ;

	  aceOutf (tsnp->ao, "\t%d\t%d", m, c) ;
	}
      else if (ti->col == 30)
	{	 
	  int c = tsnp->allC ;
	  int m = tsnp->allM ;

	  if (c >= 10)
	    ff = 100.0*m/c ;
	  if (ff >= 0)
	    aceOutf (tsnp->ao, "\t%.2f", ff) ; 
	  else
	    aceOutf (tsnp->ao, "\t-10") ;
	}
      else if (ti->col == 40)
	{	 
	  vTXT txt2 = vtxtTokenSort (tsnp->detectionTxt2, ",", h) ;
	  const char *ccp = vtxtPtr (txt2) ;
	  if (ccp)
	    aceOutf (tsnp->ao, "\t%s", ccp) ;
	  else
	    aceOut (tsnp->ao, "\t") ;
	}
      else if (ti->col == 20)
	{
	  if (ff < 0)
	    aceOutf (tsnp->ao, "\t\t\t\t") ;
	  else
	    {
	      int nLow = 0, nHigh = 0, nOk = 0, nContra = 0 ;
	      for (int r = 1 ; r < keySetMax (tsnp->covers) ; r++)
		{
		  int c = keySet (tsnp->covers, r) ;
		  int m = keySet (tsnp->mutant, r) ;
		  
		  RR *rr = arrayp (tsnp->rrs, r, RR) ;
		  if (rr->g2r)
		    continue ;

		  if (c < 10)
		    nLow++ ;
		  else
		    {
		      float f = 100 * m / c ;
		      float df = 200 * sqrt(c)/c ;
		      if (f < ff + 30 && f > ff - 30)
			nOk++ ;
		      else if (f + df > ff && f - df < ff)
			nOk++ ;
		      else
			nContra++ ;
		      nHigh++ ;		      
		    }
		}
	      aceOutf (tsnp->ao, "\t%d\t%d\t%d\t%d", nHigh, nOk, nLow, nContra) ;
	    }
	}
      else if (ti->col == 50) /* Wendell */
	{
	  AC_TABLE tt = 0 ;

	  tt = ac_tag_table (snp->Snp, "Wtrue", h) ;
	  if (!  tt) 
	    tt = ac_tag_table (snp->Snp, "Wtrue2", h) ;
	  if (tt && tt->cols >=1)
	      aceOutf (tsnp->ao, "\t%.2f", ac_table_float (tt, 0, 0, 0)) ;
	  else
	    aceOut (tsnp->ao, "\t") ;
	}

      else if (ti->col == 50) /* Wendell */
	{
	  AC_TABLE tt = 0 ;

	  tt = ac_tag_table (snp->Snp, "Wtrue", h) ;
	  if (!  tt) 
	    tt = ac_tag_table (snp->Snp, "Wtrue2", h) ;
	  if (tt && tt->cols >=1)
	      aceOutf (tsnp->ao, "\t%.2f", ac_table_float (tt, 0, 0, 0)) ;
	  else
	    aceOut (tsnp->ao, "\t") ;
	}

      else if (ti->col == 60) /* RNA_A */
	{
	  const char *run = "Sample_A_all.70_37sC8-8sR-25lC3" ;
	  int r = 0 ;
	  
	  aceOut ( tsnp->ao, "\t") ;
	  if (dictFind (tsnp->runDict, run, &r))
	    {
	      int c = keySet (tsnp->covers, r) ;
	      int m = keySet (tsnp->mutant, r) ;
	      if (c >= 20)
		aceOutf ( tsnp->ao, "%.2f", 100.0 * m /c) ;
	      else
		aceOut ( tsnp->ao, "-20") ;
	    }
	}

      else if (ti->col == 70) /* RNA_B */
	{
	  const char *run = "Sample_B_all.55_35sC8-8sR-12lC2" ;
	  int r = 0 ;

	  aceOut ( tsnp->ao, "\t") ;
	  if (dictFind (tsnp->runDict, run, &r))
	    {
	      int c = keySet (tsnp->covers, r) ;
	      int m = keySet (tsnp->mutant, r) ;
	      if (c >= 20)
		aceOutf ( tsnp->ao, "%.2f", 100.0 * m /c) ;
	      else
		aceOut ( tsnp->ao, "-20") ;
	    }
	}
      else
	snpShowTag (tsnp, snp, ti) ;
    }
  ac_free (h) ;
  return;
} /* snpBrsAllRuns */

/*************************************************************************************/

static void snpBrsTitration (TSNP *tsnp, SNP *snp)
{
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "AECDB", "ACB titration", 13, 0, 0} ,
    { "AECDB", "AECDB titration", 15, 0, 0} ,
    { "CL", "Significant AF variations among cell lines 1 to 10", 80, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 

  const char *caption =
    "Titration"
    ;

  if (snp == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;

  AC_HANDLE h = ac_new_handle () ;
  for (ti = tts ; ti->tag ; ti++)
    {
      if (snp == 0)
	{
	  aceOutf (tsnp->ao, "\t%s", ti->title) ;
	  continue ;
	}
      else if (snp == (void *)2)
	{
	  ac_free (h) ;
	  return ;
	}
      if (ti->col == 13)
	{	 
	  char *ccp = tsnp->titrationTxt3 ? vtxtPtr (tsnp->titrationTxt3) : "" ;
	  aceOutf (tsnp->ao, "\t%s", ccp ? ccp : "" ) ;
	}
      else if (ti->col == 15)
	{	 
	  char *ccp = tsnp->titrationTxt5 ? vtxtPtr (tsnp->titrationTxt5) : "" ;
	  aceOutf (tsnp->ao, "\t%s", ccp ? ccp : "" ) ;
	}
      else if (ti->col == 80) /* Clones */
	{
	  int gg = 0 ;
	  const char *group = "CL1toCL10_60_20sA1_40lR3" ;

	  aceOut (tsnp->ao, "\t") ;
	  if (dictFind (tsnp->runDict, group, &gg))
	    {
	      RR *rr = arrayp (tsnp->rrs, gg, RR) ;
	      KEYSET g2r = rr ? rr->g2r : 0 ;
	      KEYSET afs = keySetHandleCreate (h) ;
	      int jj = 0, mm = 0, cc = 0 ;
	      int minF = 9999, maxF = -20 ;

	      for (int ii = 0 ; g2r && ii < keySetMax (g2r) ; ii++)
		{
		  int r = keySet (g2r, ii) ;
		  int c = keySet (tsnp->covers, r) ;
		  int m = keySet (tsnp->mutant, r) ;
		  if (c >= 20)
		    {
		      int f = 100.0 * m /c ;
		      if (f < minF) minF = f ;
		      if (f > maxF) maxF = f ;
		      keySet (afs, jj++) = f ; 
		      mm += m ; cc += c ;
		    }
		}
	      if (0 && jj > 1)
		{
		  int kk = 0 ;
		  int ff = 100 * mm / cc ;

		  for (int j = 0 ; j < jj ; j++)
		    {
		      int f = keySet (afs, j) ;
		      if (f >= ff + 30 || f < ff - 30)
			kk++ ;
		    }
		  if (kk)
		    aceOutf (tsnp->ao, "%d", kk) ;
		}
	      if (maxF > minF + 30 && jj > 1)
		{
		  for (int ii = 0 ; g2r && ii < keySetMax (g2r) ; ii++)
		    {
		      int r = keySet (g2r, ii) ;
		      int c = keySet (tsnp->covers, r) ;
		      int m = keySet (tsnp->mutant, r) ;
		      if (c >= 20)
			{
			  int f = 100.0 * m /c ;
			  if (f > minF + 30)
			    aceOutf (tsnp->ao, "%s,", dictName (tsnp->runDict, r)) ;
			}
		    }
		}
	    }
	}
      else
	snpShowTag (tsnp, snp, ti) ;
    }
  ac_free (h) ;
  return;
} /* SnpBrsTitration */

/*************************************************************************************/

static BOOL snpFilter (TSNP *tsnp, SNP *snp)
{
  BOOL ok = snp->is20_4 ;
  int minF = tsnp->minSnpFrequency ;
  int minM = tsnp->minSnpCount ;
  int minC = tsnp->minSnpCover > 0 ? tsnp->minSnpCover : 1 ;

  return ok || tsnp->snpType >= 3 ;
  
  if (tsnp->snpType < 3 && minF + minC > 0)
    {
      AC_HANDLE h = ac_new_handle () ;
      AC_TABLE tt = ac_tag_table (snp->Snp, "BRS_counts", h) ;
      int ir ;

      ok = FALSE ;
      if (tt && tt->rows)
	for (ir = 0 ; ir < tt->rows && !ok ; ir++)
	  {
	    const char *runNam = ac_table_printable (tt, ir, 0, 0) ;
	    int c = ac_table_int (tt, ir, 1, 0) ;
	    int m = ac_table_int (tt, ir, 2, 0) ;

	    if (runNam && dictFind (tsnp->runDict, runNam, 0) 
		&& c >= minC && 100 * m >= minF * c && m >= minM
		)
		  ok = TRUE ;
	  }

      ac_free (h) ;
    }
  return ok ;
} /* snpFilter */

/*************************************************************************************/

static const char *allMethods = "VIQSgGrR" ;
/* VIQSgd   */
/* VIQSDGR     */
/* VIQSG2   */
static MM methods [] = {
  {'V', &snpVCF} ,
  {'I', &snpIdentifiers} ,
  {'Q', &snpBrsQC} ,
  {'S', &snpBrsAllRuns} ,
  {'g', &snpBrsGroupFrequency} , 
  {'d', &snpDanLiFrequency} ,
  {'r', &snpBrsRunFrequency} , 

  {'D', &snpDanLiCounts} ,

  {'R', &snpBrsRunCounts} ,

  {'G', &snpBrsGroupCounts} ,
  {'2', &snpBrsGroup2345} ,  /* meeasured in n runs of this group */
    /* not measured in n runs of this group */
    /* seen in in n runs of this group */
    /* contradiction  in n runs of this group */

  {'T', &snpBrsTitration} ,
{ 0, 0 }
} ;

/*************************************************************************************/

static BOOL snpExportSnp (TSNP *tsnp, int iSnp, int type)
{       
  const char *ccp, *lineName = "" ;
  MM *mm ;
  SNP *snp = 0 ;
  SNP snp0 ;
  
  memset (&snp0, 0, sizeof(SNP)) ;
  
  switch (type)
    {
    case 0: 
      lineName =  "### Caption" ; 
      snp = (void *) 1 ;
      break ;
    case 1: 
      lineName =  "### " ;
      snp = (void *) 0 ;
      break ;
    case 2: 
      snp = &snp0 ;
      if (! snpBrsNextSnp (tsnp, snp, iSnp))
	return FALSE ;
      break ; 
    case 3: 
      snp = (void *) 2 ;
      break ;
    }
  
  if (type != 2 || (type == 2 && snpFilter (tsnp, snp)))
    {
      static int line = 0 ;
      if (type < 2)
	aceOutf (tsnp->ao, "%s", lineName);
      else if (type == 2)
	{ 
	  if (0) aceOutf (tsnp->ao, "%d", ++line) ;
	}
      ccp = tsnp->export - 1 ;
      while (*++ccp)
	{
	  mm = methods - 1 ;
	  while (mm++, mm->cc)
	    if (mm->cc == *ccp)
	      mm->f(tsnp, snp) ;
	}
      if (type < 3)
	aceOutf (tsnp->ao, "\n" ) ;
    }
  ac_free (snp0.h) ;
  
  return TRUE ;
} /* snpExportSnp */
  
/*************************************************************************************/
/*************************************************************************************/

static void snpExportDoubleDetect (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (tsnp->outFileName, ".double_counts.txt", tsnp->gzo, h) ;
  aceOutDate (ao, "##", tsnp->dbName) ;
  aceOutDate (ao, "##", tsnp->project) ;

  int N = dictMax (tsnp->runDict) + 1 ;
  for (int r1 = 1 ; r1 < N ; r1++)
    {
      aceOut (ao, "\n") ;
      for (int r2 = r1 ; r2 <  N ; r2++)
	aceOutf (ao, "%s\t%s\t%d\n"
		 , dictName (tsnp->runDict, r1)
		 , dictName (tsnp->runDict, r2)
		 , keySet (tsnp->doubleDetect, N *r1 + r2) 
		 ) ;
    }
  ac_free (h) ;
} /* snpExportDoubleDetect */

/*************************************************************************************/
static void snpExportHistos (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (tsnp->outFileName, ".group_histos.tsf", tsnp->gzo, h) ;
  aceOutDate (ao, "##", tsnp->dbName) ;
  aceOutDate (ao, "##", tsnp->project) ;

  int jj ;
  aceOut (ao, "# histogram of allele frequencies in groups minCount == 20") ;
  for (int ii = 0 ; ii < keySetMax (tsnp->groups) ; ii++)
    {
      int gg = keySet (tsnp->groups, ii) ;
      Array histo = gg ? array (tsnp->histos, gg, Array) : 0 ;
      if (histo)
	{
	  RR *rr = arrayp (tsnp->rrs, gg, RR) ;
	  for (jj = 0 ; jj < 101 ; jj++)
	    aceOutf (ao, "\n%s___%s\t%d\ti\t%d"
		     , rr->sorting_title ? stackText (tsnp->sorting_titles, rr->sorting_title) : ""
		     , dictName (tsnp->runDict, gg)
		     , jj
		     , array (histo, jj, int)
		     ) ;
	}
      aceOut (ao, "\n") ;
    }
  ac_free (h) ;
} /* snpExportHistos */

/*************************************************************************************/

static void snpExportDetectLibs (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (tsnp->outFileName, ".detect.tsf", tsnp->gzo, h) ;

  Array aa = tsnp->detectLibs ;
  int kMax = arrayMax (aa) ;
  const char *types[] = {"0_Ref","1_Low2","2_Low5","3_Mid","4_High","5_Pure"
			 ,"6_Ref_contradiction","7_Low2_contradiction","8_Low5_contradiction","9_Mid_contradiction","10_High_contradiction","11_Pure_contradiction"
			 ,"12_Ref_True","13_Low2_True","14_Low5_True","15_Mid_True","16_High_True","17_Pure_TRue"
			 ,"18_Ref_False","19_Low2_False","20_Low5_False","21_Mid_False","22_High_False","23_Pure_False",0} ;


  for (int ii = 0 ; ii < keySetMax (tsnp->groups); ii++)
    {
      int gg = keySet (tsnp->groups, ii) ;
      RR *rr = arrayp (tsnp->rrs, gg, RR) ;
      KEYSET g2r = rr->g2r ;
      if (! g2r)
	continue ;
      for (int k = 0 ; k < 24 ; k++) 
	{
	  int kk = 24 * gg + k ;
	  if (kk < kMax)
	    {
	      int n = arr (aa, kk, int) ;
	      if (n >= 0)
		aceOutf (ao, "%s___%s\t%s\ti\t%d\n"
			 , rr->sorting_title ? stackText (tsnp->sorting_titles, rr->sorting_title) : ""
			 , dictName (tsnp->runDict, gg)
			 , types[k]
			 , n
			 ) ;
	    }
	}
    }
  ac_free (h) ;
} /* snpExportDetectLibs */

/*************************************************************************************/

static int snpGetSnpsRunsGroups (TSNP *tsnp)
{
  int ns = 0 ;
  const char *errors = 0 ;
  char *qq ;
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tt = 0 ;
  KEYSET ks ;
  DICT *dict ;

  tsnp->runs = keySetHandleCreate (tsnp->h) ;
  tsnp->groups = keySetHandleCreate (tsnp->h) ;
  tsnp->histos = arrayHandleCreate (128, Array, tsnp->h) ;
  tsnp->detectionTxt2 = vtxtHandleCreate (tsnp->h) ;
  tsnp->titrationTxt3 = vtxtHandleCreate (tsnp->h) ;
  tsnp->titrationTxt5 = vtxtHandleCreate (tsnp->h) ;
  tsnp->snpDone = bitSetCreate (16000, tsnp->h) ;
  tsnp->sorting_titles = stackHandleCreate (200, tsnp->h) ;
  tsnp->sorting_titles->textOnly = TRUE ;
  pushText (tsnp->sorting_titles, "toto") ; /* avoid zero */

  dict = tsnp->runDict = dictHandleCreate (256, tsnp->h) ;

  if (tsnp->project)
    {
      int r, g, nn = 0 ;

      tsnp->rrs = arrayHandleCreate (200, RR, tsnp->h) ;

      qq = hprintf (h, "select r, t from p in ?project where p == \"%s\" , r in p->run where r ISA runs, t in r->%s", tsnp->project, tsnp->orderBy) ;
      tt = ac_bql_table (tsnp->db
			 , qq
			 , 0
			 , "+2"
			 , &errors
			 , h
			 ) ;

      ks = tsnp->runs ; nn = 0 ;
      if (tt)
	for (int ir = 0 ; ir < tt->rows ; ir++)
	  {
	    const char *ccp = ac_table_printable (tt, ir, 0, "toto") ;
	    dictAdd (dict, ccp, &r) ;
	    keySet (ks, nn++) = r ;
	  }

      qq = hprintf (h, "select r, t, s from p in ?project where p == \"%s\", r in p->run where r ISA Groups, t in r->%s, s in r->Other_title", tsnp->project, tsnp->orderBy) ;
      tsnp->groupTitles = tt = ac_bql_table (tsnp->db
					   , qq
					   , 0
					   , "+2"
					   , &errors
					   , tsnp->h
					   ) ;
      ks = tsnp->groups ; nn = 0 ;
      if (tt)
	for (int ir = 0 ; ir < tt->rows ; ir++)
	  {
	    const char *ccp = ac_table_printable (tt, ir, 0, "toto") ;
	    if (dictAdd (dict, ccp, &g))
	      keySet (ks, nn++) = g ;

	    RR *rr = arrayp (tsnp->rrs, g, RR) ;
	    
	    rr->sorting_title = stackMark (tsnp->sorting_titles) ;
	    pushText (tsnp->sorting_titles, ac_table_printable (tt, ir, 1, "_") ) ;
	    rr->other_title = stackMark (tsnp->sorting_titles) ;
	    pushText (tsnp->sorting_titles, ac_table_printable (tt, ir, 2, ccp ) ) ;
	    rr->isGroup = TRUE ;
	  }

      qq = hprintf (h, "select r from p in ?project where p == \"%s\" , r in p->run where r#Ignore_lowQ_runs_in_group_allele_frequency", tsnp->project, 0) ;
      tsnp->runTitles = tt = ac_bql_table (tsnp->db
					   , qq
					   , 0
					   , 0
					   , &errors
					   , tsnp->h
					   ) ;

      if (tt)
	for (int ir = 0 ; ir < tt->rows ; ir++)
	  {
	    const char *ccp = ac_table_printable (tt, ir, 0, "toto") ;
	    if (dictFind (dict, ccp, &r))
	      {
		RR *rr = arrayp (tsnp->rrs, r, RR) ;
		rr->ignoreLowQ = TRUE ;
	      }
	  }

      qq = hprintf (h, "select r from p in ?project where p == \"%s\" , r in p->run where r#Low_sequence_quality", tsnp->project, 0) ;
      tt = ac_bql_table (tsnp->db
			 , qq
			 , 0
			 , 0
			 , &errors
			 , h
			 ) ;

      if (tt)
	for (int ir = 0 ; ir < tt->rows ; ir++)
	  {
	    const char *ccp = ac_table_printable (tt, ir, 0, "toto") ;
	    if (dictFind (dict, ccp, &r))
	      {
		RR *rr = arrayp (tsnp->rrs, r, RR) ;
		rr->isLowQ = TRUE ;
	      }
	  }

      qq = hprintf (h, "select g, r from p in ?project where p == \"%s\", g in p->run where g#union_of && g#snp_add_counts, r in g>>union_of where ! r#union_of, p2 in r->project where p2 == p", tsnp->project) ;
      tt = ac_bql_table (tsnp->db
			 , qq
			 , 0
			 , 0
			 , &errors
			 , h
			 ) ;
      
      /* rr->r2g: list of groups g of which r is a member */
      if (tt)
	for (int ir = 0 ; ir < tt->rows ; ir++)
	  {
	    const char *gNam = ac_table_printable (tt, ir, 0, "toto") ;
	    int g = 0 ;
	    if (dictFind (tsnp->runDict, gNam, &g))
	      {
		int r = 0 ;
		const char *rNam = ac_table_printable (tt, ir, 1, "toto") ;

		if ( dictFind (tsnp->runDict, rNam, &r))
		  {
		    RR *rr = arrayp (tsnp->rrs, r, RR) ;
		    KEYSET ks = rr->r2g ;
		    if (!ks) 
		      ks = rr->r2g = keySetCreate () ;
		    keySet (ks, keySetMax (ks)) = g ;
		    
		    rr = arrayp (tsnp->rrs, g, RR) ;
		    ks = rr->g2r ;
		    if (!ks) 
		      ks = rr->g2r = keySetCreate () ;
		    keySet (ks, keySetMax (ks)) = r ;
		  }
	      }
	  }
    }
  if (! keySetMax (tsnp->runs))
    messcrash ("No run in project %s\n", tsnp->project) ;
  if (! keySetMax (tsnp->groups))
    messout ("No group in project %s\n", tsnp->project) ;

  switch (tsnp->snpType)
    {
    case 0: 
      qq =" where s#VCF  " ;
      break ;
    case 1: 
      qq =" where  s#VCF and (s#Manual_non_monomodal || !  s->monomodal)  " ;
      break ;
    case 3: 
      qq = " where s#VCF and  (s#wtrue || s#wtrue2)  and ! s#Wfalse and ! s#Wfalse2 " ;
      break ;
    case 4: 
      qq = " where s#VCF and  (s#Wfalse || s#Wfalse2) and (s#danli_counts || s#Manual_non_monomodal || !  s->monomodal)" ;
      break ;
    case 5: 
      qq = " where s#VCF and  s#danli_counts and !s#wtrue and ! s#wtrue2  and ! s#Wfalse and ! s#Wfalse2 " ;
      break ;
    case 6: 
      qq = " where s#VCF and  !s#DanLi_counts && ! s#Wendell_counts and (s#Manual_non_monomodal || !  s->monomodal)" ;
      break ;
    case 7: 
      qq = " s#danli_counts and (s#Manual_non_monomodal || !  s->monomodal), m in s->danli_counts where m, x in m[1], y in m[3] where 100 * x > 98 * y" ;
      break ;
    case 20: 
      qq = " where s#VCF and  s#coding  and (s#Manual_non_monomodal || !  s->monomodal)" ;
      break ;
    case 21: 
      qq = " where s#VCF and  s#UTR_3prime  and (s#Manual_non_monomodal || !  s->monomodal)" ;
      break ;
    case 22: 
      qq = " where s#VCF and  s#A2G  and (s#Manual_non_monomodal || !  s->monomodal)" ;
      break ;
    case 23: 
      qq = " where s#VCF and  s#G2A  and (s#Manual_non_monomodal || !  s->monomodal)" ;
      break ;
    }
  char *qqC =  tsnp->capture ? hprintf (h, "g in ?gene where g->capture == \"%s\", s in g->variant ", tsnp->capture) : "s in ?variant" ;
  char *qq2 = hprintf (h, "select s from %s %s", qqC, qq)  ;

  tsnp->snps = ac_bql_table (tsnp->db
			     , qq2
			     , 0
			     , 0
			     , &errors
			     , tsnp->h
			     ) ;
  
  ns = tsnp->snps ? tsnp->snps->rows : 0  ;
  if (! ns)
    messcrash ("No snp in project query %s\n", qq2) ;
  
  fprintf (stderr, " snpGetSnps got %d snps, %d runs , %d groups in project %s\n"
	   , ns
	   , tsnp->runs ? keySetMax (tsnp->runs) : 0
	   , tsnp->groups ? keySetMax (tsnp->groups) : 0 
	   , tsnp->project
	   ) ;
  
  if (tsnp->doub)
    tsnp->doubleDetect = keySetHandleCreate (tsnp->h) ;

  ac_free (h) ;
  return ns ;
} /* snpGetSnpsRunsGroups */

/*************************************************************************************/

static void snpExportHeader (TSNP *tsnp)
{
  aceOutf (tsnp->ao, "## Project %s, quality control report, File %s : %s\n"
	   , tsnp->project
	   , aceOutFileName (tsnp->ao)
	   , timeShowNow () 
	   ) ;

  return ;
} /* snpExportHeader */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// snpsummary: complete SNP summary for the Magic pipeline\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, October 2022, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "// Given a project name\n"
	    "//    import all data from TSNP_DB, export a single SNP table\n"
	    "//    one line per snp, with snp-metadata, then groups or runs counts in comulm chapters\n"
	    "//    a smaller table can be obtained using different options\n"
	    "//\n"
	    "// Syntax:\n"
	    "// snpsummary -db <dir> [options]\n"
	    "//   -db <dir>: snp-database directory (usually: tmp/TSNP_DB/$zone)\n"
	    "// All other parameters are optional and may be specified in any order\n"
	    "//   -p\n" 
	    "//   --project <projectName> export counts to runs and groups belonging to the project (usually $MAGIC)\n" 
	    "//      default : only export the identifiers\n"
	    "//   -o output_file_prefix\n"
	    "//      all exported files will be named output_file_prefix.action\n" 
            "//   --gzo : gzip all output files\n"
	    "//   -e\n" 
	    "//   --export [TPbafpmg...] : only export some groups of columns, in the requested order\n"
	    "//      default: if -export is not specified, export all columns in default order\n"
	    "//            I: snp identifiers\n"
	    "//            G: groups counts\n"
	    "//            g: groups allele frequencies\n"
	    "//            R: runs counts\n"
	    "//            r: runs allele frequencies\n"
	    "//            D: Dan Li counts\n"
	    "//            d: Dan Li allele frequencies\n"
	    "//   --orderBy <tag> : sort the runs by this tag\n"
	    "//      default: sorting_title\n"
	    "//      Insert a blank line if the sorting value looks like A__B and A changes\n"
	    "//   --minSnpCover <int> : coverage default 20\n"
	    "//   --minSnpCount <int> : coverage default 4\n"
	    "//   --minSnpFrequency <int> : per cent mutant / coverage default 10\n"	    
	    "//   --snpType [0,1,2,3] \n"
	    "//       Only export  SNPs seen at least 4 times, cover >=20 m/c>2/100 in a run of this project \n"
	    "//       0 : any_snp\n"	    
	    "//       1 : skip the monomodal snps\n"
	    "//       2 : only the rejected snps\n"
	    "//       3 : only the DanLi snps\n"
	    "//       5 : only the Wendell true positives snps, single transcript per locus\n"
	    "//       6 : only the Wendell true negative snps, single transcript per locus\n"
	    "//       7 : only the Wendell true positives snps, any transcript\n"
	    "//       8 : only the Wendell true negative snps, any transcript\n"
	    "//       10 : gene captured by A1 A2 R1 R2 R3 I1 I2 I3  and snp is Wtrue || Dan_Li\n"
	    "//       9 : no SNP, only the title lines\n"
	    "//    --doubleDetect\n"
	    "//       Count  snp detected at in at least 2 libraries at 10:4 and at least once at 20:8\n"
	    "// Caveat:\n"
	    "//   Caption lines at the top start with a #\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  snpsummary -help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  TSNP tsnp ;
  AC_HANDLE h = 0 ;
  int ns = 0 ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&tsnp, 0, sizeof (TSNP)) ;
  tsnp.h = h ;
  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;
  tsnp.dbName = "NULL" ;
  getCmdLineOption (&argc, argv, "-db", &tsnp.dbName) ; 
  getCmdLineOption (&argc, argv, "--db", &tsnp.dbName) ; 

  tsnp.export =  allMethods ;
  tsnp.project = "NULL" ;
  getCmdLineOption (&argc, argv, "-e", &tsnp.export) ;
  getCmdLineOption (&argc, argv, "--export", &tsnp.export) ;
  getCmdLineOption (&argc, argv, "-p", &tsnp.project) ;
  getCmdLineOption (&argc, argv, "--project", &tsnp.project) ;
  getCmdLineOption (&argc, argv, "--capture", &tsnp.capture) ;
  getCmdLineOption (&argc, argv, "--select", &tsnp.select) ;
  getCmdLineInt (&argc, argv, "--snpType", &tsnp.snpType) ;

  tsnp.histo =   getCmdLineBool (&argc, argv, "--histos") ;
  tsnp.doub =   getCmdLineBool (&argc, argv, "--doubleDetect") ;
  tsnp.titration =   getCmdLineBool (&argc, argv, "--titration") ;
  tsnp.justDetected =   getCmdLineBool (&argc, argv, "--justDetected") ;
  tsnp.count_libs =   getCmdLineBool (&argc, argv, "--countLibs") ;
  tsnp.unique =   getCmdLineBool (&argc, argv, "--unique") ;
  tsnp.tsfOut =   getCmdLineBool (&argc, argv, "--tsf") ;
  tsnp.gzo =   getCmdLineBool (&argc, argv, "--gzo") ;

  tsnp.minSnpCount = 4 ;
  tsnp.minSnpFrequency = 10 ;
  tsnp.minSnpCover = 20 ;
  getCmdLineInt (&argc, argv, "--minSnpFrequency", &tsnp.minSnpFrequency) ;
  getCmdLineInt (&argc, argv, "--minSnpCover", &tsnp.minSnpCover) ;
  getCmdLineInt (&argc, argv, "--minSnpCount", &tsnp.minSnpCount) ;
  if (! tsnp.dbName)
    {
      fprintf (stderr, "Sorry, missing parameter -db, please try snpsummary -help\n") ;
      exit (1) ;
    }
  if (tsnp.dbName)
    {
      const char *errors ;

      tsnp.db = ac_open_db (tsnp.dbName, &errors);
      if (! tsnp.db)
	{
	  fprintf (stderr, "Failed to open db %s, error %s", tsnp.dbName, errors) ;
	  exit (1) ;
	}
    }
  tsnp.orderBy = "Sorting_title_2" ;
  getCmdLineOption (&argc, argv, "--orderBy", &tsnp.orderBy) ;
  
  getCmdLineOption (&argc, argv, "-o", &tsnp.outFileName) ; 
  tsnp.ao = aceOutCreate (tsnp.outFileName, ".SNP_summary.txt", tsnp.gzo, h) ;
  if (tsnp.tsfOut)
    tsnp.aoTsf = aceOutCreate (tsnp.outFileName, ".SNP_summary.tsf", tsnp.gzo, h) ;
  aceOutDate (tsnp.ao, "##", tsnp.dbName) ;

  if (! tsnp.ao)
    messcrash ("cannot create output file %s\n", tsnp.outFileName ? tsnp.outFileName : EMPTY ) ;


  if (argc > 1) usage (messprintf ("Unknown parameters %s", argv[1])) ;

  if (tsnp.doub)
    {
      tsnp.ao204 = aceOutCreate (tsnp.outFileName, ".nDetectingLibsPerSnp.tsf", tsnp.gzo, h) ;
      aceOutDate (tsnp.ao204, "##", messprintf ("Number of libraries in the project measuring the SNP at coverage >= 20, frequency >=2, variant count >= 4 in  at least 2 libraries of project %s", tsnp.project)) ; 
      aceOut (tsnp.ao204, "#Variant\tProject\tFormat\tN libraries\tTrue/False\tMonomodal\n") ;
    }

  if (tsnp.count_libs  && tsnp.project)
    {
      const char *errors = 0 ;
      const char *qq = hprintf (tsnp.h, "select c,n,r from p in ?project where p == \"%s\", c in p->Compare where c#Detection, r in c->runs , n in r[1] where n>0", tsnp.project) ;
      tsnp.countLibsTable = ac_bql_table (tsnp.db, qq, 0, 0, &errors, tsnp.h)  ;
      if (! tsnp.countLibsTable)
	tsnp.count_libs = FALSE ;
      else
	{
	  tsnp.detectLibs = arrayHandleCreate (1024, int, tsnp.h) ;
	  tsnp.aoDT = aceOutCreate (tsnp.outFileName, ".detected_snps.tsf", tsnp.gzo, tsnp.h) ;
	}
    }
  if (tsnp.titration && tsnp.project)
    {
      const char *errors = 0 ;
      const char *qq = hprintf (tsnp.h, "select c,n,r,t from p in ?project where p == \"%s\", c in p->Compare where c#Profile, r in c->runs , n in r[1] where n>0, t in c->sorting_title", tsnp.project) ;
      tsnp.titrationTable = ac_bql_table (tsnp.db, qq, 0, "+t+c+n+r", &errors, tsnp.h)  ;
      if (! tsnp.titrationTable)
	tsnp.titration = FALSE ;
      tsnp.aoTitration = aceOutCreate (tsnp.outFileName, ".titration.tsf", tsnp.gzo, tsnp.h) ;
    }


  snpExportHeader (&tsnp) ;

  ns = snpGetSnpsRunsGroups (&tsnp) ;
  snpExportSnp (&tsnp, 0, 0) ; /* chapter captions */
  snpExportSnp (&tsnp, 0, 1) ; /* column titles */

  if (ns)
    {
      int iSnp ;

      for (iSnp = 0 ; iSnp < ns ; iSnp++)
	snpExportSnp (&tsnp, iSnp, 2) ;
    }
  
  if (tsnp.tsfOut)
    snpExportSnp (&tsnp, 0, 3) ; /* stats */

  if (tsnp.doub)
    snpExportDoubleDetect (&tsnp) ;
  if (tsnp.count_libs)
    {
      snpExportDetectLibs (&tsnp) ;
    }
  
  if (tsnp.histo)
    snpExportHistos (&tsnp) ;

  ac_free (h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
} /* main */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
