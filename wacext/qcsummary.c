/*
 * authors: Danielle and Jean Thierry-Mieg, NCBI, 
 * 15 Sept 2014
 * qcsummary.c
 *   Input: MetaDB->runs and ali
 *   Output: a wholesale QC table with all kinds of results
 *      all : all columns
 *        r  readFate
 *        p  pairFate
 *        s  sponge
 */

#include "../wac/ac.h"

typedef struct qcStruct {
  const char *project ;
  const char *dbName ;
  const char *export ;
  int runType ; /* 0: any, 1=sublib; 2=runs ; 3=group */
  ACEOUT ao ;
  Array runs ;
  DICT *bloomDict ;
  AC_DB db ;
  AC_HANDLE h ;
  const char *Etargets[4] ;
  const char *EtargetsBeau[4] ;
  const char *orderBy ;
  const char *outFileName ;
  BOOL gzi, gzo ;
  BOOL hasSL ;
  BOOL noIntron ;
  BOOL hasCapture ;
  BOOL TT ; /* limit counts to targeted genes */
  const char *externalFiles ;
  Array eTitles, eAaa, eCaptions ;
} QC ;

typedef enum { T_ZERO = 0, T_Raw, T_RawA, T_RawKb,  T_Length, T_Seq, T_Read, T_kb, T_Rejected, T_Unaligned, T_aliFrag, T_cFrag, T_uFrag, T_A, T_T, T_G, T_C, T_MAX } T_TYPE ;
typedef enum { F_ZERO = 0, F_Hide, F_p10, F_perCentRead, T_FMAX } T_FORMAT ;
typedef struct rcStruct {
  AC_HANDLE h ;
  const char *nam ;
  AC_OBJ run, ali, srr, srx, srp, sample ;
  BOOL Paired_end ;
  float var[T_MAX], geneSigExpressed ;
} RC ;

typedef void (*QCFunc)(QC *qc, RC *rc) ;
typedef struct mmStruct { char cc ; QCFunc f ;} MM ;
typedef struct tagtitleStruct { const char *tag, *title ; int col ; T_TYPE setVar ; T_FORMAT format ;} TT ;

#define EMPTY ""

/*************************************************************************************/
/*************************** actual work *********************************************/
/*************************************************************************************/

static void qcChapterCaption (QC *qc, TT *tts, const char *caption) 
{
  TT *ti ;
  int n = 0 ;
  char buf[8128] ;
  
  memset (buf, 0, sizeof (buf)) ;

  if (caption)
    aceOutf (qc->ao, "\t\t%s", caption) ;

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
      aceOutf (qc->ao, "%s", buf) ;
    }
  return ;
} /*  qcChapterCaption */

/*************************************************************************************/

static int qcShowMultiTag (QC *qc, AC_TABLE tt)
{
  AC_HANDLE h = ac_new_handle () ;
  int n, nn = 0, ir, jr ;
  const char *ccp, *ccq ;
  vTXT txt = vtxtCreate () ;
  Stack s = stackHandleCreate (50, h) ;


  vtxtPrintf (txt, "toto") ;
  /* enter the data backwards, do not duplicate up except for column 1 which we always keep */
  for (ir = tt->rows - 1 ; ir >= 0 ; ir--)
    for (jr = tt->cols - 1 ; jr >= 0 ; jr--)
      {
	ccp = ac_table_printable (tt, ir, jr, 0) ;
	if (ccp)
	  {
	    ccq = hprintf(h, "#%s#", ccp) ;
	    if (
		(jr == 0 && !  strstr (vtxtPtr (txt), ccq) && ( ir == tt->rows - 1 ||  strcmp (ccp, ac_table_printable (tt, ir+1, 0, 0)))) ||
		(jr > 0 && ! strstr (vtxtPtr (txt), ccp))
		)
	      {
		vtxtPrintf (txt, "%s", ccq) ;
		pushText (s, ccp) ;
		nn++ ;
	      }
	  }
      }
  /* export */
  n = nn ;
  while (n--)
    {
      ccp = popText(s) ;
      if (! strcmp (ccp, "polyA")) ccp = "polyA selected" ;
      aceOutf (qc->ao, "%s%s"
	       , n == nn - 1 ? "" : ", "
	       , ccp
	       ) ;
      ac_free (txt) ;
    }
  return nn ;
} /* qcShowMultiTag */

/*************************************************************************************/

static int qcShowTable (QC *qc, AC_TABLE tt, int nw)
{
  int ir, jr ;
  
  if (! tt || ! tt->rows)
    return 0 ;
  
  for (ir = 0 ; ir < tt->rows ; ir++)
    for (jr = 0 ; jr < tt->cols ; jr++)
      {
	const char *ccp = ac_table_printable (tt, ir, jr, 0) ;
	if (ccp)
	  {
	    if (!strncasecmp (ccp, "sra", 3)) ccp += 3 ;
	    if (!strcasecmp (ccp, "sraUnspecified_RNA"))
	      continue ;
	    if (! jr)
	      { /* reformat the stage to allow ordering in excell, Danielle, 2025_01_18 */
		int n = 0 ;
	        if (!strcmp (ccp, "EE")) n = 10 ;
	        if (!strcmp (ccp, "L1")) n = 20 ;
	        if (!strcmp (ccp, "L2")) n = 30 ;
	        if (!strcmp (ccp, "L3")) n = 40 ;
	        if (!strcmp (ccp, "L4")) n = 50 ;
	        if (!strcmp (ccp, "YA")) n = 60 ;
	        if (!strcmp (ccp, "OA")) n = 70 ;
	        if (!strcmp (ccp, "Adult")) n = 80 ;
	        if (!strcmp (ccp, "Mixed_stages")) n = 90 ;
	        if (!strcmp (ccp, "Mixed_embryo")) n = 13 ;
		if (!strcmp (ccp, "Embryo")) n = 12 ;
	        if (!strcmp (ccp, "Molt")) n = 29 ;
	        if (!strcmp (ccp, "Dauer")) n = 35 ;
		if (n) aceOutf (qc->ao, "%d ", n) ;
	      }
	    switch (jr)
	      {
	      case 0:
		if (nw) aceOut (qc->ao, "; ") ;
		break ;
	      case 1:
		aceOut (qc->ao, ": ") ;
		break ;
	      default:
		aceOut (qc->ao, ": ") ;
		break ;
	      }		
	    aceOut (qc->ao, ccp) ;
	    nw++ ;
	  }
      }
  return nw ;
} /* qcShowTable */

/*************************************************************************************/

static void qcShowMillions (QC *qc, float z)
{
  if (z == 0)
    aceOutf (qc->ao, "\t0") ;
  else if (z > 100000 || z < -100000)
    aceOutf (qc->ao, "\t%.3f",  z/1000000) ;
  else
    aceOutf (qc->ao, "\t%.6f",  z/1000000) ;
} /* qcShowMillions */

/*************************************************************************************/

static void qcShowPercent (QC *qc, RC *rc, long int z, BOOL p)
{
  if (p)
    {
      float z1 =  rc->var[T_Read] ? 100 * z / rc->var[T_Read] : 0 ;
      
      if (rc->var[T_Read] > 0) 
	{
	  aceOutf (qc->ao, "\t") ;
	  if (0)
	    aceOutf (qc->ao, "%f", z1) ;
	  else
	    aceOutPercent (qc->ao, z1) ;
	}
      else 
	aceOutf (qc->ao, "\t") ;
    }
  else
    {
      if (z == 0)
	aceOutf (qc->ao, "\t0") ;
      else
	aceOutf (qc->ao, "\t%ld", z) ;
    }
} /*  qcShowPercent */

/*************************************************************************************/
static void qcShowTag (QC *qc, RC *rc, TT *ti)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  const char *ccp = EMPTY ;
  float z = 0 ;

  if (strcmp (ti->tag, "Spacer"))
    {
      ccp = EMPTY ; tt = ac_tag_table (rc->ali, ti->tag, h) ;
      if (!tt) tt = ac_tag_table (rc->run, ti->tag, h) ;
      if (tt && tt->rows && tt->cols >= ti->col)
	{
	  ccp = ac_table_printable (tt, 0, ti->col, EMPTY) ;
	  z = ac_table_float (tt, 0, ti->col, 0) ;
	  if (ti->setVar)
	    rc->var[ti->setVar] = z ;
	}
    }

  switch (ti->format)
    {
    default:
      aceOutf (qc->ao, "\t%s", ccp) ;
      break ;
    case F_p10:
      aceOutf (qc->ao, "\t%.1f", z/10) ;
      break ;
    case F_perCentRead:
      aceOutf (qc->ao, "\t%.4f", rc->var[T_Read] ? 100.0 * z/rc->var[T_Read] : 0) ;
      break ;
    case F_Hide:
      break ;
    }
} /* qcShowTag */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

static void qcSetAli (QC *qc, RC *rc)
{
  TT *ti, tts[] = {
    { "Raw_data", "Number of raw reads", 2, T_Raw, F_Hide} ,
    { "Raw_data", "Number of raw reads", 2, T_RawA, F_Hide} ,
    { "Raw_data", "Number of raw bases", 4, T_RawKb, F_Hide} ,
    { "Accepted", "Number of distinct sequences", 0, T_Seq, F_Hide} ,
    { "Accepted", "Number of reads", 2, T_Read, F_Hide} ,
    { "Accepted", "Total sequenced (kb)", 4, T_kb, F_Hide} , 
    { "Rejected", "Number of unaligned reads", 2, T_Rejected, F_Hide} , 
    { "Unaligned", "Number of unaligned reads", 2, T_Unaligned, F_Hide} , 
    { 0, 0, 0, 0, 0 }
  } ;
  rc->ali = ac_tag_obj (rc->run, "Ali", rc->h) ;
  if (rc && rc->run && ! rc->srr) 
    rc->srr = ac_tag_obj (rc->run, "SRR", rc->h) ;
  if (rc && rc->run && ! rc->srx) 
    rc->srx = ac_tag_obj (rc->run, "SRX", rc->h) ;
  if (rc && rc->run && ! rc->srp) 
    rc->srp = ac_tag_obj (rc->run, "SRP", rc->h) ;
  if (rc && rc->run && ! rc->sample) 
    rc->sample = ac_tag_obj (rc->run, "Sample", rc->h) ;
  if (rc && !rc->srr  && rc->run && ac_has_tag (rc->run, "Sublibraries"))
    {
      AC_OBJ subrun = ac_tag_obj (rc->run, "Sublibraries", rc->h) ;
      rc->srr = subrun ? ac_tag_obj (subrun, "SRR", rc->h) : 0 ;
      if (!rc->srp)
	rc->srp = subrun ? ac_tag_obj (subrun, "SRP", rc->h) : 0 ;
      if (!rc->sample)
	rc->sample = subrun ? ac_tag_obj (subrun, "Sample", rc->h) : 0 ;
      ac_free (subrun) ;
    }

  for (ti = tts ; ti->tag ; ti++)
    qcShowTag (qc, rc, ti) ;
  rc->Paired_end = ac_has_tag (rc->run, "Paired_end") ;
  if (rc->Paired_end)
    {
      int i, dx = 0 ;
      i = rc->var[T_Raw] ; if (i > 0 && i%2 == 1 &&  rc->var[T_Raw] > 0)  { rc->var[T_Raw]-- ; dx--; }
      i = rc->var[T_Read] ; if (i > 0 && i%2 == 1 && rc->var[T_Read] > 0)  { rc->var[T_Read]-- ; dx++ ; }
      rc->var[T_Rejected] += dx ;
    }
  /* the raw length is not precomputed, compute it rounding at  .1 base  */
  {
    int ln = 5 + 10 * (1000 *  rc->var[T_kb] / ( rc->var[T_Read] ?  rc->var[T_Read] : 1)) ;
    rc->var[T_Length] = ln/10.0 ;
  }
  if (! rc->var[T_Length]) rc->var[T_Length] = 1 ; /* prevent zero divide */
  if (rc->var[T_Raw] == 0) rc->var[T_Raw] = rc->var[T_Read] +  rc->var[T_Rejected] ;
  if (rc->var[T_RawA] == 0) rc->var[T_RawA] = rc->var[T_Read] ;
  if (rc->var[T_RawKb] == 0) rc->var[T_RawKb] = rc->var[T_kb] ;
} /* qcSetAli */

/*************************************************************************************/
/* Data characteristics before alignment */
static void qcBeforeAli (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  const char *ccp ;
  int ir ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute", "Average read length (nt)", 1, 0, 0} ,
    { "Compute", "Average fragment multiplicity\tMaximal fragment multiplicity\tMillion distinct sequences\tMillion raw reads\tMegabases sequenced", 2, 0, 0} ,
    { "Compute", "%A\t%T\t%G\t%C\t%N\t%GC", 3, 0, 0} ,
    {  0, 0, 0, 0, 0}
  }; 
  const char *caption =
    "Insert sizes measured in pairs"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	aceOutf (qc->ao, "\t%s", ti->title) ;
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  switch (ti->col)
	    {
	    case 1: /* read length */
	      {
		double mbp = 0, mr = 0, mbp1 = 0, mr1 = 0 ; ;

		aceOutf (qc->ao, "\t") ;
		tt = ac_tag_table (rc->ali, "Letter_profile", h) ;
		if (tt)
		  {
		    int ir, nf = 0, n1 ; float n2 ;
		    char buf[64] ; buf[0] = 0 ;
		    for (ir = 0 ; ir < tt->rows  ; ir++)
		      {
			ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
			if (ir == 0) 
			  strncpy (buf, ccp, 64) ;
			if (strncmp (buf, ccp, 64))
			  {
			    mbp1 += mbp ; mr1 += mr ;
			    if (mr)
			      {
				int uk =  .005 + mbp/mr ;
				float uf =  mbp/mr ;
				int u1 = 100 * uf + 0.5 ;
				if (100 * uk == u1)
				  aceOutf (qc->ao, "%s%s:%d", nf ? ", " : EMPTY, buf, uk) ;
				else
				  aceOutf (qc->ao, "%s%s:%.2f", nf ? ", " : EMPTY, buf, uf) ;
			      }
			    else
			      aceOutf (qc->ao, "%s%s:0",  nf ? ", " : EMPTY, buf) ;
			    nf++ ; mbp = mr = 0 ;
			    strncpy (buf, ccp, 64) ;
			  }
			n1 = ac_table_int (tt, ir, 1, 0) ;
			n2 = ac_table_float (tt, ir, 7, 0) ; 
			mr += n2 ;
			mbp += n1 * n2 ;
			if (0) fprintf(stderr, "n1=%d  n2=%f mr=%f\n",n1,n2, mr); 
		      }
		    mbp1 += mbp ; mr1 += mr ;
		    if (tt && tt->rows > 999)
		      {
			mr = rc->var[T_Read] ;
			mbp = rc->var[T_Length] ;
		      }
		    if (mr)
		      {
			int uk =  .005 + mbp/mr ;
			float uf =  mbp/mr ;
			int u1 = 100 * uf + 0.5 ;
			if (100 * uk == u1)
			  aceOutf (qc->ao, "%s%s:%d", nf ? ", " : EMPTY, buf, uk) ;
			else
			  aceOutf (qc->ao, "%s%s:%.2f", nf ? ", " : EMPTY, buf, uf) ;
		      }
		    else
		      aceOutf (qc->ao, "%s%s:0",  nf ? ", " : EMPTY, buf) ;
		  }
		else
		  aceOutf (qc->ao, "%d", (int)(rc->var[T_Length] + .5)) ;
		break ;
	      }
	    case 2:
	      aceOutf (qc->ao, "\t%.2f", rc->var[T_Seq] ? rc->var[T_Read]/rc->var[T_Seq] : 0) ;
	      aceOutf (qc->ao, "\t%d", ac_tag_int (rc->ali, "Maximal_read_multiplicity", 0)) ;
	      if (0) /* round values */
		{
		  qcShowMillions (qc, rc->var[T_Seq]) ; 
		  qcShowMillions (qc, rc->var[T_Read]) ;
		  aceOutf (qc->ao, "\t%.3f", rc->var[T_kb]/1000) ;   
		}
	      else  /* give the decimals */
		{
		  aceOutf (qc->ao, "\t%f", rc->var[T_Seq]/1000000.0) ;
		  aceOutf (qc->ao, "\t%f", rc->var[T_Read]/1000000.0) ;
		  aceOutf (qc->ao, "\t%f", rc->var[T_kb]/1000) ;   
		}
	      break ;
	    case 3: /* ATGC */
	      {
		float z, za, zt, zg, zc, zn ;
		tt = ac_tag_table (rc->ali, "ATGC_kb", h) ;
		z = za = zt = zg = zc = zn = 0 ;
		for (ir = 0 ; tt && ir < tt->rows && ir < 1 ; ir++)
		  {
		    za = ac_table_float (tt, ir, 5, 0) ; 
		    zt = ac_table_float (tt, ir, 6, 0) ; 
		    zg = ac_table_float (tt, ir, 7, 0) ; 
		    zc = ac_table_float (tt, ir, 8, 0) ; 
		    zn = ac_table_float (tt, ir, 9, 0) ; 
		  } 
		z = za + zt + zg + zc + zn ; if (z == 0) z = 1 ;
		aceOutf (qc->ao, "\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f"
			 , 100 * za/z, 100 *zt/z, 100 *zg/z, 100 * zc/z, 100 *zn/z
			 , 100 * (zg + zc) /z
			 ) ;
	      }
	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  ac_free (h) ;
  return;
}  /* qcBeforeAli */

/*************************************************************************************/
/* Data characteristics before alignment */
static void qcAli (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  const char *ccp ;
  float z, zb, zc ;
  int ir ;
  TT *ti, tts[] = { 
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute", "Million raw reads", 20, 0, 0} ,  /* copied from qcBeforeAli */
    { "Compute", "Million reads aligned on any target by Magic\tAverage length aligned per read (nt)\tMb aligned on any target\t% Mb aligned on any target before clipping\t% length aligned on average before clipping", 1, 0, 0} , 
    { "Compute", "% length aligned after clipping adaptors and barcodes (nt)", 2, 0, 0 } ,
    { "Compute", "% Reads aligned on any target", 3, 0, 0 } ,
 
    {  0, 0, 0, 0, 0}
  }; 
  const char *caption =
    "Global alignment statistics"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	aceOutf (qc->ao, "\t%s", ti->title) ;
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  switch (ti->col)
	    {
	    case 999:
	      z =  rc->var[T_G] +rc->var[T_C] ;
	      if (z > 1000) z = 1000 ;
	      aceOutf (qc->ao, "\t%.1f", z/10) ;
	      break ;
	    case 20: /* Million raw reads */
	      qcShowMillions (qc, rc->var[T_Read]) ;
	      break ;
	    case 1:
	      ccp = EMPTY ; tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
 	      z = zb = zc = 0 ;
 	      for (ir = 0 ; tt && ir < tt->rows ; ir++)
 		{
 		  ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
 		  if (! strcasecmp (ccp, "any"))
 		    {
 		      z = ac_table_float (tt, ir, 3, 0) ;
 		      zb = ac_table_float (tt, ir, 5, 0) ;
 		      zc = ac_table_float (tt, ir, 7, 0) ;
 		    }
 		} 

	      aceOutf (qc->ao, "\t%.3f", z/1000000) ;
	      aceOutf (qc->ao, "\t%.2f", zc) ;  /* average length aligned */
	      aceOutf (qc->ao, "\t%.3f", zb/1000) ; /* Mb aligned on any target */

	      aceOutf (qc->ao, "\t%.2f", rc->var[T_kb] ? 100 *zb / rc->var[T_kb] : 0) ; /* % Mb aligned on any target before clipping */
	      {
		float z1 = rc->var[T_kb], z2 = 0 ;
		AC_TABLE tt2 = ac_tag_table (rc->ali, "Unaligned", h) ;
		z2 = tt2 ? ac_table_float (tt2, 0, 4,0) : 0 ; /* kb Unaligned */
		aceOutf (qc->ao, "\t%.2f ", 100*zb/(z1 - z2)) ;  /* % length aligned on average before clipping */
	      }
	   
	 
	      break ;
	    case 2:   /*  Average clipped length */
	      {
		int ir ;
		float z = -1, zClipped = -1 ;
		AC_TABLE tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
		
		for (ir = 0 ; tt &&  ir < tt->rows; ir++)
		  {
		    if (! strcasecmp (ac_table_printable (tt, ir, 0, ""), "any"))
		      {
			z =  ac_table_float (tt, ir, 7, -1) ; 
			zClipped = ac_table_float (tt, ir, 11, -1) ;
			break ;
		      }
		  }
		aceOutf (qc->ao, "\t") ;
		if (z > -1)
		  {
		    aceOutPercent (qc->ao, 100.00 * z/zClipped) ;
		  } 
		  
	      }
	      break ;

	    case 3:
	      ccp = EMPTY ; tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
 	      z = zb = zc = 0 ;
 	      for (ir = 0 ; tt && ir < tt->rows ; ir++)
 		{
 		  ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
 		  if (! strcasecmp (ccp, "any"))
 		    {
 		      z = ac_table_float (tt, ir, 3, 0) ;
 		    }
 		} 

	      aceOutf (qc->ao, "\t%.2f", rc->var[T_Read] ? 100 *z / rc->var[T_Read] : 0) ;
	 
	      break ;
	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }

  ac_free (h) ;
  return;
}  /* qcAli */

/*************************************************************************************/

static void qcAvLengthAli (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  int ir, irAny = -1 ;
  float z, avClipped = 1 ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Jump5", "Bases clipped 5prime in read 1, so that a maximal number of reads alig from base 1", 0, 0, 0} ,
    { "Jump5", "Bases clipped 5prime in read 2, so that a maximal number of reads alig from base 1", 1, 0, 0} ,
    { "Compute", "Average length after clipping adaptors and barcodes (nt)", 1, 0, 0 } ,
    
    { "Compute", "any:Average length aligned (nt)", 2, 0, 0 } ,
    { "Compute", "any1:Average length aligned in read 1 (nt)", 2, 0, 0 } ,
    { "Compute", "any2:Average length aligned in read 2 (paired end)", 2, 0, 0 } ,

    /*  { "Compute", "RT_:Average length aligned on reference transcriptome (nt)", 2, 0, 0 } , :"reference transcriptome used in project like RefSeq or Flybase or Encode (37.104 models) */
    { "Compute", "KT_RefSeq:Average length aligned on RefSeq (nt)", 2, 0, 0 } , /*:"Encode.37.70 Jan2013" */
    { "Compute", "MT_EBI:Average length aligned on EBI (nt)", 2, 0, 0 } , /*:"Encode.37.70 Jan2013" */
    { "Compute", "ET_av:Average length aligned on AceView (nt)", 2, 0, 0 } , /*: AceView" */
    { "Compute", "LT_magic:Average length aligned on magic (nt)", 2, 0, 0 } , /*: Future AceView" */
    { "Compute", "NT_:Average length aligned on Other (nt)", 2, 0, 0 } , /*: other unnoficial  annotation */

    { "Compute", "Z_genome:Average length aligned on Genome (nt)", 2, 0, 0 } ,
    { "Compute", "A_mito:Average length aligned on Mitochondria (nt)", 2, 0, 0 } ,
    { "Compute", "B_rRNA:Average length aligned on rRNA (nt)", 2, 0, 0 } ,

    { "Compute", "C_chloro:Average length aligned on Chloroplast (nt)", 2, 0, 0 } ,
    { "Compute", "QT_smallRNA:Average length aligned on small RNAs (nt)", 2, 0, 0 } ,
    { "Compute", "1_DNASpikeIn:Average length aligned on DNA spikeIn (nt)", 2, 0, 0 } ,
    { "Compute", "0_SpikeIn:Average length aligned on RNA spikeIn (nt)", 2, 0, 0 } ,
    { "Compute", "z_gdecoy:Average length aligned on Imaginary genome specificity control (nt)", 2, 0, 0 } ,
    

    { "Spacer", "", 0, 0, 0} ,

    { "Compute", "any:Average % length aligned (nt)", 2, 0, F_p10 } ,
    { "Compute", "any1:Average % length aligned in read 1 (nt)", 2, 0, F_p10 } ,
    { "Compute", "any2:Average % length aligned in read 2 (paired end)", 2, 0, F_p10 } ,

    { "Compute", "KT_RefSeq:Average % length aligned on RefSeq", 2, 0, F_p10} , /*:"Encode.37.70 Jan2013" */
    { "Compute", "MT_EBI:Average % length aligned on EBI", 2, 0, F_p10} , /*:"Encode.37.70 Jan2013" */
    { "Compute", "ET_av:Average % length aligned on AceView", 2, 0, F_p10} , /*: AceView" */
    { "Compute", "LT_magic:Average % length aligned on magic", 2, 0, F_p10} , /*: Future AceView" */
    { "Compute", "NT_:Average % length aligned on Other", 2, 0, F_p10} , /*: other unnoficial  annotation */

    { "Compute", "Z_genome:Average % length aligned on Genome", 2, 0, F_p10} ,
    { "Compute", "A_mito:Average % length aligned on Mitochondria", 2, 0, F_p10} ,
    { "Compute", "B_rRNA:Average % length aligned on rRNA", 2, 0, F_p10} ,

    { "Compute", "C_chloro:Average % length aligned on Chloroplast", 2, 0, F_p10} ,
    { "Compute", "QT_smallRNA:Average % length aligned on small RNAs", 2, 0, F_p10} ,
    { "Compute", "1_DNASpikeIn:Average % length aligned on DNA spikeIn", 2, 0, F_p10} ,
    { "Compute", "0_SpikeIn:Average % length aligned on RNA spikeIn", 2, 0, F_p10} ,
    { "Compute", "z_gdecoy:Average % length aligned on Imaginary genome specificity control", 2, 0, F_p10} ,
   
    {  0, 0, 0, 0, 0}
  } ; 
  
  const char *caption =
    hprintf (h, "Average length (nt) aligned per target in  %s"
	     "Average length aligned per target in nucleotide (table to the left) and in percentage of the length of the read, after clipping adaptors and barcodes (table to the right)." 
	     "Average lengths are reported only if more than 1000 reads are aligned. "
	     "This is seldom the case for the imaginary decoy genome, used as a mapping specificity control, yet the histogram of length aligned on this target is used to set the thresholds for filtering alignments of low quality." 
	     "Ribosomal RNA (rRNA) is encoded in the genome in many copies with slight sequence variations, hence when we align the archetypic rRNA precursor, usually from RefSeq or GenBank, the quality of the alignments is not expected to be perfect."  
	     "Similarly, the small RNA targets, including tRNA, miRNA and other small RNAs often delineates only the mature product, and the actual genes are larger than the gene models. Hence usually only a fraction of the read lengths align on these targets. This is also true for the imperfect but useful human RNA genes from UCSC."
	     
	     , qc->project
	     ) ;
  if (rc == (void *) 1)
    {
      h = ac_new_handle () ;
      qcChapterCaption (qc, tts, caption) ;
      ac_free (h) ;
      return ;
    }
  
  for (ti = tts ; ti->tag ; ti++)
    {
      int nw = 0 ;
      char buf[256] ;

      if (rc == 0)
	{
	  const char *ccp = strchr (ti->title, ':') ;
	  if (!ccp || ti->col != 2)
	    ccp = ti->title ;
	  else 
	    ccp++ ;
	  aceOutf (qc->ao, "\t%s", ccp) ;
	  continue ; 
	}
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
 
      tt = ac_tag_table (rc->ali, "nh_Ali", h) ;

      if (! strcmp (ti->tag, "Compute"))
	{
	  aceOutf (qc->ao, "\t") ;
	  switch (ti->col)
	    {
	    case 1:   /*  Average clipped length */
	      irAny = ir ;
	      for (ir = 0 ; tt &&  ir < tt->rows; ir++)
		{
		  if (! strcasecmp (ac_table_printable (tt, ir, 0, ""), "any"))
		    irAny = ir ;
		}

	      avClipped = z = irAny > 0 ? ac_table_float (tt, irAny, 11, -1) : -1 ; 
	      if (z > -1)
		{
		  nw++ ;
		  if (ti->format == F_p10)
		    aceOutPercent (qc->ao, 100.00 * z/avClipped) ;
		  else
		    aceOutf (qc->ao, "%.2f", z) ;
		}
	      break ;
	      
	    case 2:   /*  Average aligned length */
	      strncpy (buf, ti->title, 255) ;
	      { char *cp = strchr (buf, ':') ;
		if (cp) *cp = 0 ;
	      } 
	      for (ir = 0 ; tt && ir < tt->rows; ir++)
		{
		  if (! strcasecmp (ac_table_printable (tt, ir, 0, ""), buf))
		    {
		      z =  ac_table_float (tt, ir, 3, -1) ; 
		      avClipped = ac_table_float (tt, ir, 11, -1) ;
		      if (z > 1000) /* at least 1000 tags aligned */
			{
			  z = irAny > 0 ? ac_table_float (tt, ir, 7, -1) : -1 ; 
			  if (z > -1)
			    {
			      nw++ ;
			      if (ti->format == F_p10)
				aceOutPercent (qc->ao, 100.00 * z/avClipped) ;
			      else
				aceOutf (qc->ao, "%.2f", z) ;
			    }
			}
		      break ;
		    }
		}
	      break ;
	    }
	  if (! nw) aceOut (qc->ao, EMPTY) ;
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  
  ac_free (h) ;
  return;
}  /* qcAvLengthAli */
  
/*************************************************************************************/

static void  qcStrandedness (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  int ir ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,

    { "Compute", "A_mito:Reads mapping on plus strand of", 1,0,0},
    { "Compute", "A_mito:Reads mapping on minus strand of", 2,0,0},
    { "Compute", "A_mito:Reads mapping on both strands of", 3,0,0},
    { "Compute", "A_mito:% reads mapping on plus strand of", 1,0,F_p10},
    { "Compute", "A_mito:% reads mapping on minus strand of", 2,0,F_p10},
    { "Compute", "A_mito:% reads mapping on both strand of", 3,0,F_p10},

    { "Compute", "B_rRNA:Reads mapping on plus strand of", 1,0,0},
    { "Compute", "B_rRNA:Reads mapping on minus strand of", 2,0,0},
    { "Compute", "B_rRNA:Reads mapping on both strands of", 3,0,0},
    { "Compute", "B_rRNA:% reads mapping on plus strand of", 1,0,F_p10},
    { "Compute", "B_rRNA:% reads mapping on minus strand of", 2,0,F_p10},
    { "Compute", "B_rRNA:% reads mapping on both strand of", 3,0,F_p10},

    { "Compute", "Reads mapping on plus strand of", 101,0,0},
    { "Compute", "Reads mapping on minus strand of", 201,0,0},
    { "Compute", "Reads mapping on both strands of", 301,0,0},
    { "Compute", "% reads mapping on plus strand of", 101,0,F_p10},
    { "Compute", "% reads mapping on minus strand of", 201,0,F_p10},
    { "Compute", "% reads mapping on both strand of", 301,0,F_p10},

    { "Compute", "Reads mapping on plus strand of", 102,0,0},
    { "Compute", "Reads mapping on minus strand of", 202,0,0},
    { "Compute", "Reads mapping on both strands of", 302,0,0},
    { "Compute", "% reads mapping on plus strand of", 102,0,F_p10},
    { "Compute", "% reads mapping on minus strand of", 202,0,F_p10},
    { "Compute", "% reads mapping on both strand of", 302,0,F_p10},
#ifdef JUNK
    /* 2026_01_13 i added the mito target so i remove annotation 3 so that the number of columns is not modified */
    { "Compute", "Reads mapping on plus strand of", 103,0,0},
    { "Compute", "Reads mapping on minus strand of", 203,0,0},
    { "Compute", "Reads mapping on both strands of", 303,0,0},
    { "Compute", "% reads mapping on plus strand of", 103,0,F_p10},
    { "Compute", "% reads mapping on minus strand of", 203,0,F_p10},
    { "Compute", "% reads mapping on both strand of", 303,0,F_p10},
#endif
    { "Compute", "Z_genome:Reads mapping on plus strand of", 1,0,0},
    { "Compute", "Z_genome:Reads mapping on minus strand of", 2,0,0},
    { "Compute", "Z_genome:Reads mapping on both strands of", 3,0,0},
    { "Compute", "Z_genome:% reads mapping on plus strand of", 1,0,F_p10},
    { "Compute", "Z_genome:% reads mapping on minus strand of", 2,0,F_p10},
    { "Compute", "Z_genome:% reads mapping on both strand of", 3,0,F_p10},

    { "Compute", "0_SpikeIn:Reads mapping on plus strand of", 1,0,0},
    { "Compute", "0_SpikeIn:Reads mapping on minus strand of", 2,0,0},
    { "Compute", "0_SpikeIn:Reads mapping on both strands of", 3,0,0},
    { "Compute", "0_SpikeIn:% reads mapping on plus strand of", 1,0,F_p10},
    { "Compute", "0_SpikeIn:% reads mapping on minus strand of", 2,0,F_p10},
    { "Compute", "0_SpikeIn:% reads mapping on both strand of", 3,0,F_p10},

    {  0, 0, 0, 0, 0}
  } ; 
  
  const char *caption =
    "Strandedness, number of fragments aligning  per strand of the indicated target."
    "The estimation is naive and does not use any a priori knowledge derived from the protocol."
    "The expected values depend on the target. In any experiment, the fragment strandedness measured on the "
    "genome should be close to 50%, but some RNA-Seq experiments use a strand labelling or strand selection. For "
    "instance, the SOLiD RNA-Seq protocol is very accurately and systematically stranded: 99.99% of the fragments"
    "align on the plus strand of transcripts. In the Illumina UDG protocol, between 88% and 99% of the fragments"
    "map on the reverse strand."
    "The most accurate measure is provided by alignements to ribosomal genes which are never ambiguous, but"
    "this method is not usable if ribozero treatment was applied, because in that case the remaining rRNA fragments "
    "are few and strongly rearranged. In the other targets, some genes overlap in trans in the genome, therefore "
    "the fragments aligning on the overlap, where the 2 genes are antisense to each other, can be interpreted "
    "as forward relative to one gene, and reverse relative to the other. These fragments are counted as "
    "ambiguous."
    "Note that the strandedness is only computed when more than 1000 fragments are aligned on a given"
    "target. Note also that the plus percentage is computed only on the informative fragments,"
    "plus + minus, so the sum of the 3 percentages plus+minus+ambiguous may exceed 100%."
    ;

  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {
      int nw = 0 ;
      int readLn = 1 ;
      char buf[256] ;
      char buf2[256] ;

      if (rc == 0)
	{
	  char *cp ;
	  strncpy (buf, ti->title, 255) ;
	  cp = strchr (buf, ':') ;
	  if (ti->col > 100)
	    {
	      int c = (((ti->col-1) % 100)%3) & 0x3 ;
	      const char *E  = qc->Etargets[c] ;
	      if (qc->noIntron)
		E = (c == 0 ? "P_genome" : qc->Etargets[c+1]) ;
	      if (E)
		aceOutf (qc->ao, "\t%s %s", ti->title, E) ;
	    }
	  else if (!cp)
	    aceOutf (qc->ao, "\t%s", buf) ;
	  else /* mamipulate the title */
	    {
	      *cp++ = 0 ;                      /* now buf == A_rrna */
	      aceOutf (qc->ao, "\t%s", cp) ;
	      cp = strchr (buf, '_') ;
	      if (cp)
		aceOutf (qc->ao, " %s", cp+1) ;  /* rrna */
	    }
	  continue ; 
	}
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}

      tt = ac_tag_table (rc->ali, "stranding", h) ;

      if (! strcmp (ti->tag, "Compute"))
	{
	  float up = 0, um = 0, ua = 0, uz = -1 ;

	  int col = ti->col ;
	  if (col > 100) col /= 100 ;
	  aceOutf (qc->ao, "\t") ;
	  switch (col)
	    {
	    case 1:   /*  plus strand */
	    case 2:   /*  minus strand */
	    case 3:   /*  both strands */
	      strncpy (buf, ti->title, 255) ;
	      { char *cp = strchr (buf, ':') ;
		if (cp) *cp = 0 ;
	      }

	      if (ti->col > 100)
		{
		  int c = (((ti->col-1) % 100)%3) & 0x3 ;
		  const char *E  = qc->Etargets[c] ;
		  if (qc->noIntron)
		    {
		      if (c == 0)
			readLn = 5 + 10 * (1000 *  rc->var[T_kb] / ( rc->var[T_Read] ?  rc->var[T_Read] : 1)) ;
		      E = (c == 0 ? "P_genome" : qc->Etargets[c+1]) ;
		    }
		  if (E)
		    strncpy (buf, E, 250) ;
		  else
		    break ;
		}
	      for (ir = 0 ; tt && ir < tt->rows; ir++)
		if (! strcmp (ac_table_printable (tt, ir, 0, ""), "I_introns"))
		  strcpy (buf, "introns") ;  /* specific of Magic2 */
	      for (ir = 0 ; tt && ir < tt->rows; ir++)		
		{
		  strncpy (buf2, ac_table_printable (tt, ir, 0, ""), 255) ;
		   { char *cp = strchr (buf2, '.') ;
		     if (cp) *cp = 0 ;
		   } 
		   if (! strcasecmp (buf, buf2) || ! strcasecmp (buf, buf2 + 3))
		    {
		      up += ac_table_float (tt, ir, 2, 0) ;  /* LOOP because in groups we may have .f and .f2 cases */
		      um += ac_table_float (tt, ir, 4, 0) ;
		      ua += ac_table_float (tt, ir, 6, 0) ;
		      uz += ac_table_float (tt, ir, 2*col, 0) ; 
		    }
		}
	      up = up + um + ua ;
	      if (up > 1000) /* at least 1000 tags aligned */
		{
		  if (uz > -1)
		    {
		      nw++ ;
		      if (ti->format == F_p10)
			aceOutPercent (qc->ao, 100.00 * uz/up) ;
		      else
			aceOutf (qc->ao, "%.2f", uz/readLn) ;
		    }
		}

	      break ;
	    }
	  if (! nw) aceOut (qc->ao, EMPTY) ;
	}
      else
	qcShowTag (qc, rc, ti) ;
    }

  ac_free (h) ;
  return;
} /* qcStrandedness */

/*************************************************************************************/
/*************************************************************************************/

static void qcMismatchTypes (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  int ir ;
  float zAny = 0, denominator = 0 ;
  float zTransition = 0, zTransversion = 0 ;
  float zSlidingInsertion = 0, zSlidingDeletion = 0 ;
  float zInsertion = 0, zDeletion = 0 ;
  float zInsertion1 = 0, zInsertion2 = 0, zInsertion3 = 0 ;
  float zInsertionA = 0, zInsertionT = 0, zInsertionG = 0, zInsertionC = 0 ;
  float zDeletion1 = 0, zDeletion2 = 0, zDeletion3 = 0 ; 
  float zDeletionA = 0, zDeletionT = 0, zDeletionG = 0, zDeletionC= 0 ; 
  float zSub[12] ;

  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
 
    { "Compute", "Megabases uniquely aligned and used for mismatch counts", 1, 0, 0} ,
    { "Compute", "Total number of mismatches", 2, 0, 0} ,
    { "Compute", "Mismatches per kb aligned", 21, 0, 0} ,

    { "Compute", "Transitions", 21, 0, 0} ,
    { "Compute", "Transversions", 21, 0, 0} ,

    { "Compute", "1, 2 or 3 bp insertions not in polymers", 12, 0, 0} ,
    { "Compute", "1, 2 or 3 bp deletions not in polymers", 13, 0, 0} ,
    { "Compute", "1, 2 or 3 bp insertions in polymers", 14, 0, 0} ,
    { "Compute", "1, 2 or 3 bp deletions in polymers", 15, 0, 0} ,

    { "Compute", "A>G", 21, 0, 0} , /* Transitions */
    { "Compute", "T>C", 21, 0, 0} ,
    { "Compute", "G>A", 21, 0, 0} ,
    { "Compute", "C>T", 21, 0, 0} ,

    { "Compute", "A>T", 21, 0, 0} ,  /* Transversions */
    { "Compute", "T>A", 21, 0, 0} ,
    { "Compute", "G>C", 21, 0, 0} ,
    { "Compute", "C>G", 21, 0, 0} ,
    { "Compute", "A>C", 21, 0, 0} ,
    { "Compute", "T>G", 21, 0, 0} ,
    { "Compute", "G>T", 21, 0, 0} ,
    { "Compute", "C>A", 21, 0, 0} ,

    { "Compute", "Insert A", 21, 0, 0} ,
    { "Compute", "Insert T", 21, 0, 0} ,
    { "Compute", "Insert G", 21, 0, 0} ,
    { "Compute", "Insert C", 21, 0, 0} ,
    { "Compute", "Delete A", 21, 0, 0} ,
    { "Compute", "Delete T", 21, 0, 0} ,
    { "Compute", "Delete G", 21, 0, 0} ,
    { "Compute", "Delete C", 21, 0, 0} ,

    { "Compute", "Single Insertions", 20, 0, 0} ,
    { "Compute", "Single Deletions", 21, 0, 0} ,
    { "Compute", "Double insertions", 22, 0, 0} ,
    { "Compute", "Double deletions", 23, 0, 0} ,
    { "Compute", "Triple insertions", 22, 0, 0} ,
    { "Compute", "Triple deletions", 25, 0, 0} ,
 
    { "Compute", "Transitions per kb", 21, 0, 0} ,
    { "Compute", "Transversions per kb", 21, 0, 0} ,
    { "Compute", "Insertions per kb", 21, 0, 0} ,
    { "Compute", "Deletions per kb", 21, 0, 0} ,

    { "Compute", "% transitions", 21, 0, 0} ,
    { "Compute", "% transversions", 21, 0, 0} ,
    { "Compute", "% indels", 21, 0, 0} ,

    { "Compute", "% insertions not in polymers", 21, 0, 0} ,
    { "Compute", "% deletions not in polymers", 71, 0, 0} ,
    { "Compute", "% insertions in polymers", 21, 0, 0} ,
    { "Compute", "% deletions in polymers", 21, 0, 0} ,

    {  0, 0, 0, 0, 0}
  } ; 
  
  const char *caption =
    "Number of mismatches, per type" ;
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;


  memset (zSub, 0, sizeof (zSub)) ;
  /*
Distribution of mismatches in best unique alignments, absolute, observed, counts
*/

  for (ti = tts ; ti->tag ; ti++)
    {
      char buf[256] ;

      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title) ;
	  continue ; 
	}
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}

      tt = ac_tag_table (rc->ali, "Error_profile", h) ;

      if (! strcmp (ti->tag, "Compute"))
	{
	  float zz = 0 ;
	  switch (ti->col)
	    {
	    case 1:   /* Megabases uniquely aligned and used for mismatch counts */
	      for (ir = 0 ; tt &&  ir < tt->rows ; ir++)
		{
		  if (ir == 0)
		    {
		      strcpy (buf, ac_table_printable (tt, ir, 0, "xxx")) ;
		      zz = ac_table_float (tt, ir, 3, 0) ;
		    }
		  else if (strcmp (buf, ac_table_printable (tt, ir, 0, "xxx")))
		    {
		      strcpy (buf, ac_table_printable (tt, ir, 0, "xxx")) ;
		      zz +=  ac_table_float (tt, ir, 3, 0) ;
		    }
		}
	      aceOutf (qc->ao, "\t%.0f", zz) ;
	      denominator = zz ;
	      break ;
	      
	    case 2:   /* Total number of mismatches */
	      for (ir = 0 ; tt && ir < tt->rows ; ir++)
		{
		  const char *ccp = ac_table_printable (tt, ir, 1, "xxx") ;
		  zz =  ac_table_float (tt, ir, 2, 0) ;
		  if (! strcasecmp (ccp, "Any"))
		    continue ;
		  zAny += zz ;
		  if (ccp[1] == '>' && ccp[3] == 0) /* substitution */
		    {
		      int i = 0 ;
		      const char *cq, *subTypes = "a>g t>c g>a c>t a>t t>a g>c c>g a>c t>g g>t c>a" ;

		      cq = strstr (subTypes, ccp) ;
		      i = cq ? (cq - subTypes) / 4 : -1 ;
		      if (i >= 4) zTransversion += zz ;
		      else if (i >= 0) zTransition += zz ;

		      if (i >= 0 && i < 12) zSub[i] += zz ;
		    }
		  else if (! strncmp (ccp, "+", 1))
		    {
		      zInsertion += zz ;
		      if (! strncmp (ccp, "+a", 2))
			zInsertionA += zz ;
		      if (! strncmp (ccp, "+t", 2))
			zInsertionT += zz ;
		      if (! strncmp (ccp, "+g", 2))
			zInsertionG += zz ;
		      if (! strncmp (ccp, "+c", 2))
			zInsertionC += zz ;
		      if (! strncmp (ccp, "+++", 3))
			zInsertion3 += zz ;
		      else if (! strncmp (ccp, "++", 2))
			zInsertion2 += zz ;
		      else 
			zInsertion1 += zz ;
		    }
		  else if (! strncmp (ccp, "-", 1))
		    {
		      zDeletion += zz ;
		      if (! strncmp (ccp, "-a", 2))
			zDeletionA += zz ;
		      if (! strncmp (ccp, "-t", 2))
			zDeletionT += zz ;
		      if (! strncmp (ccp, "-g", 2))
			zDeletionG += zz ;
		      if (! strncmp (ccp, "-c", 2))
			zDeletionC += zz ;
		      if (! strncmp (ccp, "---", 3))
			zDeletion3 += zz ;
		      else if (! strncmp (ccp, "--", 2))
			zDeletion2 += zz ;
		      else 
			zDeletion1 += zz ;
		    }
		  else if (! strncmp (ccp, "*+", 2))
		    {
		      zSlidingInsertion += zz ;
		      if (! strncmp (ccp + 1, "+a", 2))
			zInsertionA += zz ;
		      if (! strncmp (ccp + 1, "+t", 2))
			zInsertionT += zz ;
		      if (! strncmp (ccp + 1, "+g", 2))
			zInsertionG += zz ;
		      if (! strncmp (ccp + 1, "+c", 2))
			zInsertionC += zz ;
		      if (! strncmp (ccp + 1, "+++", 3))
			zInsertion3 += zz ;
		      else if (! strncmp (ccp + 1, "++", 2))
			zInsertion2 += zz ;
		      else 
			zInsertion1 += zz ;
		    }
		  else if (! strncmp (ccp, "*-", 2))
		    {
		      zSlidingDeletion += zz ;   
		      if (! strncmp (ccp + 1, "-a", 2))
			zDeletionA += zz ;
		      if (! strncmp (ccp + 1, "-t", 2))
			zDeletionT += zz ;
		      if (! strncmp (ccp + 1, "-g", 2))
			zDeletionG += zz ;
		      if (! strncmp (ccp + 1, "-c", 2))
			zDeletionC += zz ;
		      if (! strncmp (ccp + 1, "---", 3))
			zDeletion3 += zz ;
		      else if (! strncmp (ccp + 1, "--", 2))
			zDeletion2 += zz ;
		      else 
			zDeletion1 += zz ;
		    }
		}

	      /* report all at once */
	      aceOutf (qc->ao, "\t%.0f", zAny) ;
	      aceOutf (qc->ao, "\t%.5f", denominator > 0 ? zAny / (1000 * denominator) : 0) ;
	      aceOutf (qc->ao, "\t%.0f\t%.0f", zTransition, zTransversion) ;

	      aceOutf (qc->ao, "\t%.0f\t%.0f", zInsertion, zDeletion ) ;
	      aceOutf (qc->ao, "\t%.0f\t%.0f", zSlidingInsertion, zSlidingDeletion) ;

	      { 
		int i ; 
		for (i = 0 ; i < 12 ; i++)
		  aceOutf (qc->ao, "\t%.0f", zSub[i]) ;
	      }

	      aceOutf (qc->ao, "\t%.0f", zInsertionA) ;
	      aceOutf (qc->ao, "\t%.0f", zInsertionT) ;
	      aceOutf (qc->ao, "\t%.0f", zInsertionG) ;
	      aceOutf (qc->ao, "\t%.0f", zInsertionC) ;
	      aceOutf (qc->ao, "\t%.0f", zDeletionA) ;
	      aceOutf (qc->ao, "\t%.0f", zDeletionT) ;
	      aceOutf (qc->ao, "\t%.0f", zDeletionG) ;
	      aceOutf (qc->ao, "\t%.0f", zDeletionC) ;

	      aceOutf (qc->ao, "\t%.0f\t%.0f", zInsertion1, zDeletion1) ;
	      aceOutf (qc->ao, "\t%.0f\t%.0f", zInsertion2, zDeletion2) ;
	      aceOutf (qc->ao, "\t%.0f\t%.0f", zInsertion3, zDeletion3) ;

	      aceOutf (qc->ao, "\t%.5f", denominator > 0 ? zTransition / (1000 * denominator) : 0) ;
	      aceOutf (qc->ao, "\t%.5f", denominator > 0 ? zTransversion / (1000 * denominator) : 0) ;
	      aceOutf (qc->ao, "\t%.5f", denominator > 0 ? (zInsertion + zSlidingInsertion)/ (1000 * denominator) : 0) ;
	      aceOutf (qc->ao, "\t%.5f", denominator > 0 ? (zDeletion + zSlidingDeletion)/ (1000 * denominator) : 0) ;

	      if (zAny)
		{
		  aceOutf (qc->ao, "\t%.2f", 100 * zTransition / zAny) ;
		  aceOutf (qc->ao, "\t%.2f", 100 * zTransversion / zAny) ;
		  aceOutf (qc->ao, "\t%.2f", 100 * (zInsertion + zDeletion + zSlidingInsertion + zSlidingDeletion) / zAny) ;

		  aceOutf (qc->ao, "\t%.2f", 100 * zInsertion / zAny) ;
		  aceOutf (qc->ao, "\t%.2f", 100 * zDeletion / zAny) ;
		  aceOutf (qc->ao, "\t%.2f", 100 * zSlidingInsertion / zAny) ;
		  aceOutf (qc->ao, "\t%.2f", 100 * zSlidingDeletion / zAny) ; 
		}
	      else
		aceOut (qc->ao, "\t\t\t\t\t\t\t") ;
	      break ;
	    default: /* alread  treated in case 2 */
	      break ;

	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  
  ac_free (h) ;
  return;
} /* qcMismatchTypes */

/*************************************************************************************/

static void qcSnpTypesDo (QC *qc, RC *rc, BOOL isRejected)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  int ir ;
  float zAny = 0, denominator = 0 ;
  float zTransition = 0, zTransversion = 0 ;
  float zSlidingInsertion = 0, zSlidingDeletion = 0 ;
  float zInsertion = 0, zDeletion = 0 ;
  float zInsertion1 = 0, zInsertion2 = 0, zInsertion3 = 0 ;
  float zInsertionA = 0, zInsertionT = 0, zInsertionG = 0, zInsertionC = 0 ;
  float zDeletion1 = 0, zDeletion2 = 0, zDeletion3 = 0 ; 
  float zDeletionA = 0, zDeletionT = 0, zDeletionG = 0, zDeletionC= 0 ; 
  float zSub[12] ;

  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,

    { "Compute", "Total number of variant alleles", 2, 0, 0} ,
    { "Compute", "SNP support per kb aligned", 21, 0, 0} ,

    { "Compute", "Transitions", 21, 0, 0} ,
    { "Compute", "Transversions", 21, 0, 0} ,

    { "Compute", "1, 2 or 3 bp insertions not in polymers", 12, 0, 0} ,
    { "Compute", "1, 2 or 3 bp deletions not in polymers", 13, 0, 0} ,
    { "Compute", "1, 2 or 3 bp insertions in polymers", 14, 0, 0} ,
    { "Compute", "1, 2 or 3 bp deletions in polymers", 15, 0, 0} ,

    { "Compute", "A>G", 21, 0, 0} , /* Transitions */
    { "Compute", "T>C", 21, 0, 0} ,
    { "Compute", "G>A", 21, 0, 0} ,
    { "Compute", "C>T", 21, 0, 0} ,

    { "Compute", "A>T", 21, 0, 0} ,  /* Transversions */
    { "Compute", "T>A", 21, 0, 0} ,
    { "Compute", "G>C", 21, 0, 0} ,
    { "Compute", "C>G", 21, 0, 0} ,
    { "Compute", "A>C", 21, 0, 0} ,
    { "Compute", "T>G", 21, 0, 0} ,
    { "Compute", "G>T", 21, 0, 0} ,
    { "Compute", "C>A", 21, 0, 0} ,

    { "Compute", "Insert A", 21, 0, 0} ,
    { "Compute", "Insert T", 21, 0, 0} ,
    { "Compute", "Insert G", 21, 0, 0} ,
    { "Compute", "Insert C", 21, 0, 0} ,
    { "Compute", "Delete A", 21, 0, 0} ,
    { "Compute", "Delete T", 21, 0, 0} ,
    { "Compute", "Delete G", 21, 0, 0} ,
    { "Compute", "Delete C", 21, 0, 0} ,

    { "Compute", "Single Insertions", 20, 0, 0} ,
    { "Compute", "Single Deletions", 21, 0, 0} ,
    { "Compute", "Double insertions", 22, 0, 0} ,
    { "Compute", "Double deletions", 23, 0, 0} ,
    { "Compute", "Triple insertions", 22, 0, 0} ,
    { "Compute", "Triple deletions", 25, 0, 0} ,
 
    { "Compute", "% transitions", 21, 0, 0} ,
    { "Compute", "% transversions", 21, 0, 0} ,
    { "Compute", "% indel", 21, 0, 0} ,

    { "Compute", "% insertions not in polymers", 21, 0, 0} ,
    { "Compute", "% deletions not in polymers", 71, 0, 0} ,
    { "Compute", "% insertions in polymers", 21, 0, 0} ,
    { "Compute", "% deletions in polymers", 21, 0, 0} ,
    { "Compute", "A to G RNA edition: ratio of A>G relative to (A>G + T>C) (1/2 means no visible edition)", 21, 0, 0} ,

    {  0, 0, 0, 0, 0}
  } ; 
  
  const char *caption =
    isRejected ?
    "Number of rejected variant alleles, per type. Designed for clinical cohorts with sequenced diploid tissues. SNPs which are monomodal across all samples, usually with MAF in the 1-30% range, come most probably from systematic noise and are rejected."
    :
    "Number of variant alleles, per type."
    ;
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;


  memset (zSub, 0, sizeof (zSub)) ;
  /*
Distribution of mismatches in best unique alignments, absolute, observed, counts
*/

  for (ti = tts ; ti->tag ; ti++)
    {
      char buf[256] ;

      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s %s", isRejected ? "Rejected " : "",  ti->title) ;
	  continue ; 
	}
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}


      if (! strcmp (ti->tag, "Compute"))
	{
	  float zz = 0 ;
	  int colN = 2 ;


	  tt = ac_tag_table (rc->ali, "SNP_profile", h) ;
	  switch (ti->col)
	    {
	    case 1:   /* Megabases uniquely aligned and used for mismatch counts */
	      for (ir = 0 ; tt &&  ir < tt->rows ; ir++)
		{
		  if (ir == 0)
		    {
		      strcpy (buf, ac_table_printable (tt, ir, 0, "xxx")) ;
		      zz = ac_table_float (tt, ir, colN, 0) ;
		    }
		  else if (strcmp (buf, ac_table_printable (tt, ir, 0, "xxx")))
		    {
		      strcpy (buf, ac_table_printable (tt, ir, 0, "xxx")) ;
		      zz +=  ac_table_float (tt, ir, colN, 0) ;
		    }
		}
	      aceOutf (qc->ao, "\t%.0f", zz) ;
	      denominator = zz ;
	      break ;
	      
	    case 2:   /* Total number of mismatches */
	      for (ir = 0 ; tt && ir < tt->rows ; ir++)
		{
		  const char *ccp = ac_table_printable (tt, ir, 1, "xxx") ;
		  zz =  ac_table_float (tt, ir, isRejected ? 3 : 2, 0) ;
		  if (! strcasecmp (ccp, "Any"))
		    continue ;
		  zAny += zz ;
		  if (ccp[2] == '_') ccp += 3 ;
		  if (ccp[1] == '>' && ccp[3] == 0) /* substitution */
		    {
		      int i = 0 ;
		      const char *cq, *subTypes = "a>g t>c g>a c>t a>t t>a g>c c>g a>c t>g g>t c>a" ;

		      cq = strstr (subTypes, ccp) ;
		      i = cq ? (cq - subTypes) / 4 : -1 ;
		      if (i >= 4) zTransversion += zz ;
		      else if (i >= 0) zTransition += zz ;

		      if (i >= 0 && i < 12) zSub[i] += zz ;
		    }
		  else if (! strncmp (ccp, "+", 1))
		    {
		      zInsertion += zz ;
		      if (! strncmp (ccp, "+a", 2))
			zInsertionA += zz ;
		      if (! strncmp (ccp, "+t", 2))
			zInsertionT += zz ;
		      if (! strncmp (ccp, "+g", 2))
			zInsertionG += zz ;
		      if (! strncmp (ccp, "+c", 2))
			zInsertionC += zz ;
		      if (! strncmp (ccp, "+++", 3))
			zInsertion3 += zz ;
		      else if (! strncmp (ccp, "++", 2))
			zInsertion2 += zz ;
		      else 
			zInsertion1 += zz ;
		    }
		  else if (! strncmp (ccp, "-", 1))
		    {
		      zDeletion += zz ;
		      if (! strncmp (ccp, "-a", 2))
			zDeletionA += zz ;
		      if (! strncmp (ccp, "-t", 2))
			zDeletionT += zz ;
		      if (! strncmp (ccp, "-g", 2))
			zDeletionG += zz ;
		      if (! strncmp (ccp, "-c", 2))
			zDeletionC += zz ;
		      if (! strncmp (ccp, "---", 3))
			zDeletion3 += zz ;
		      else if (! strncmp (ccp, "--", 2))
			zDeletion2 += zz ;
		      else 
			zDeletion1 += zz ;
		    }
		  else if (! strncmp (ccp, "*+", 2))
		    {
		      zSlidingInsertion += zz ;
    		      if (! strncmp (ccp + 1, "+a", 2))
			zInsertionA += zz ;
		      if (! strncmp (ccp + 1, "+t", 2))
			zInsertionT += zz ;
		      if (! strncmp (ccp + 1, "+g", 2))
			zInsertionG += zz ;
		      if (! strncmp (ccp + 1, "+c", 2))
			zInsertionC += zz ;
		      if (! strncmp (ccp + 1, "+++", 3))
			zInsertion3 += zz ;
		      else if (! strncmp (ccp + 1, "++", 2))
			zInsertion2 += zz ;
		      else 
			zInsertion1 += zz ;
		    }
		  else if (! strncmp (ccp, "*-", 2))
		    {
		      zSlidingDeletion += zz ;   
		      if (! strncmp (ccp + 1, "-a", 2))
			zDeletionA += zz ;
		      if (! strncmp (ccp + 1, "-t", 2))
			zDeletionT += zz ;
		      if (! strncmp (ccp + 1, "-g", 2))
			zDeletionG += zz ;
		      if (! strncmp (ccp + 1, "-c", 2))
			zDeletionC += zz ;
		      if (! strncmp (ccp + 1, "---", 3))
			zDeletion3 += zz ;
		      else if (! strncmp (ccp + 1, "--", 2))
			zDeletion2 += zz ;
		      else 
			zDeletion1 += zz ;
		    }
		}

	      /* report all at once */
	      aceOutf (qc->ao, "\t%.0f", zAny) ;
	      aceOutf (qc->ao, "\t%.5f", denominator > 0 ? zAny / (1000 * denominator) : 0) ;
	      aceOutf (qc->ao, "\t%.0f\t%.0f", zTransition, zTransversion) ;

	      aceOutf (qc->ao, "\t%.0f\t%.0f", zInsertion, zDeletion ) ;
	      aceOutf (qc->ao, "\t%.0f\t%.0f", zSlidingInsertion, zSlidingDeletion) ;

	      { 
		int i ; 
		for (i = 0 ; i < 12 ; i++)
		  aceOutf (qc->ao, "\t%.0f", zSub[i]) ;
	      }

	      aceOutf (qc->ao, "\t%.0f", zInsertionA) ;
	      aceOutf (qc->ao, "\t%.0f", zInsertionT) ;
	      aceOutf (qc->ao, "\t%.0f", zInsertionG) ;
	      aceOutf (qc->ao, "\t%.0f", zInsertionC) ;
	      aceOutf (qc->ao, "\t%.0f", zDeletionA) ;
	      aceOutf (qc->ao, "\t%.0f", zDeletionT) ;
	      aceOutf (qc->ao, "\t%.0f", zDeletionG) ;
	      aceOutf (qc->ao, "\t%.0f", zDeletionC) ;

	      aceOutf (qc->ao, "\t%.0f\t%.0f", zInsertion1, zDeletion1) ;
	      aceOutf (qc->ao, "\t%.0f\t%.0f", zInsertion2, zDeletion2) ;
	      aceOutf (qc->ao, "\t%.0f\t%.0f", zInsertion3, zDeletion3) ;


	      if (zAny)
		{
		  aceOutf (qc->ao, "\t%.2f", 100 * zTransition / zAny) ;
		  aceOutf (qc->ao, "\t%.2f", 100 * zTransversion / zAny) ;
		  aceOutf (qc->ao, "\t%.2f", 100 * (zInsertion + zDeletion + zSlidingInsertion + zSlidingDeletion) / zAny) ;
		  aceOutf (qc->ao, "\t%.2f", 100 * zInsertion / zAny) ;
		  aceOutf (qc->ao, "\t%.2f", 100 * zDeletion / zAny) ;
		  aceOutf (qc->ao, "\t%.2f", 100 * zSlidingInsertion / zAny) ;
		  aceOutf (qc->ao, "\t%.2f", 100 * zSlidingDeletion / zAny) ; 
		  if ( zSub[0]+ zSub[1] >= 10)
		    aceOutf (qc->ao, "\t%.2f", zSub[0]/(zSub[0] + zSub[1])) ;
		  else
		    aceOut (qc->ao, "\t") ;
		}
	      else
		aceOut (qc->ao, "\t\t\t\t\t\t\t\t") ;
	      break ;
	    default: /* already  treated in case 2 */
	      break ;

	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  
  ac_free (h) ;
  return;
} /* qcSnpTypesDo */

/***********/

static void qcSnpTypes (QC *qc, RC *rc)
{ 
  qcSnpTypesDo (qc, rc, FALSE) ;
}  /* qcSnpTypes */

/***********/

static void qcSnpRejectedTypes (QC *qc, RC *rc)
{ 
  qcSnpTypesDo (qc, rc, TRUE) ;
}  /* qcSnpRejectedTypes */

/*************************************************************************************/

static void qcSnpCoding (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt, ttS, ttG, ttP ;

  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute", "Number of SNV sites tested\t% of tested SNV sites covered at least 10 times\tMeasured SNV sites\t% rejected monomodals, likely mapping or sequencing errors or RNA-edited sites", 1, 0, 0 } ,
    { "Compute", "Exonic SNV sites\tPure reference SNV (< 5% variant)\tLow frequency SNV (5-20%)\tMid frequency SNV (20-80%)\tHigh frequency SNV (80-95%)\tPure variant SNV (95-100%)", 2, 0, 0} ,
    { "Compute", "Protein changing SNV sites\tPure reference, no change in protein (< 5% variant)\tProtein changing SNV, intermediate (5-95%)\tProtein changing SNV, pure variant (95-100%)", 3, 0, 0} ,
    { "Compute", "Heterozygosity index: heterozygous/homozygous ratio variants", 4, 0, 0 } ,
    /*
    { "Compute", "% Pure reference SNV (< 5% variant)\t% Low frequency SNV (5-20%)\t% Mid frequency SNV (20-80%)\t% High frequency SNV (80-95%)\t% Pure variant SNV (95-100%)", 5, 0, 0} ,
    { "Compute", "% Protein changing SNV sites\t% Pure reference, no change in protein (< 5% variant)\t% Protein changing SNV, intermediate frequency (5-95%)\t% Protein changing SNV, pure variant (95-100%)", 6, 0, 0} ,
    */
    {  0, 0, 0, 0, 0}
  } ; 
  
  const char *caption =
    "SNV"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;
  if (rc)
    {
      ttS = ac_tag_table (rc->ali, "SNP", h) ;
      ttG = ac_tag_table (rc->ali, "Genomic", h) ;
      ttP = ac_tag_table (rc->ali, "Protein_changing", h) ;
    }

  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title) ;
	  continue ; 
	}
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
            
      tt = ttS ;
      
      if (! strcmp (ti->tag, "Compute"))
	{
	  float z1 = 0, z2 = 0 ;
	  int nT, nC, nR, nnM ;
	  switch (ti->col)
	    {
	    case 1:
	      nT = ac_tag_int (rc->ali, "Tested_sites", 0) ;
	      nnM = ac_tag_int (rc->ali, "Not_measurable_sites", 0) ; 
	      nR = ac_tag_int (rc->ali, "Rejected_sites", 0) ;
	      nT += nnM ;
	      if (nT > 0)
		{
		  nC = nT - nnM ;
		  z1 = nT ; z1 = 100.0 * nC / z1 ;
		  z2 = nC ; z2 = nC > 0 ? (100.0 * nR / z2) : 0 ;
		  aceOutf (qc->ao, "\t%d\t%.2f\t%d\t%.2f", nT, z1, nC, z2) ; 
		}
	      else
		aceOut (qc->ao, "\t\t\t\t") ;
	      break ;
	    case 2:
	      tt = ttG ;
	      if (tt)
		{
		  aceOutf (qc->ao, "\t%d\t%d\t%d\t%d\t%d\t%d"
			   , ac_table_int (tt, 0, 0, 0)
			   , ac_table_int (tt, 0, 2, 0)
			   , ac_table_int (tt, 0, 4, 0)
			   , ac_table_int (tt, 0, 6, 0) 
			   , ac_table_int (tt, 0, 8, 0)
			   , ac_table_int (tt, 0, 10, 0)			   
			   ) ;
		}
	      else
		aceOut (qc->ao, "\t\t\t\t\t\t") ;
	     
	      break ;
	    case 3: 
	      tt = ttP ;
	      if (tt)
		{
		  aceOutf (qc->ao, "\t%d\t%d\t%d\t%d"
			   , ac_table_int (tt, 0, 0, 0)
			   , ac_table_int (tt, 0, 2, 0)
			   , ac_table_int (tt, 0, 4, 0) +  ac_table_int (tt, 0, 6, 0) + ac_table_int (tt, 0, 8, 0)
			   , ac_table_int (tt, 0, 10, 0)			   
			   ) ;
		}
	      else
		aceOut (qc->ao, "\t\t\t\t") ;
	      break ;

	    case 4:  
	      nC = nT = 0 ;
	      tt = ac_tag_table (rc->ali, "Genomic", h) ;
	      if (tt)
		{
		  nC = ac_table_int (tt, 0, 6, 0) +  ac_table_int (tt, 0, 8, 0) ;
		  nT = ac_table_int (tt, 0, 10, 0) ;
		}
	      aceOut (qc->ao, "\t") ;  
	      if (nC + nT >= 10)
		{
		  z1 = (nC < 100 * nT ?  nC/(1.0 * nT) : 1000) ;
		  aceOutf (qc->ao, "%.2f", z1) ;
		}
		
	      break ;
#ifdef USELESS
	    case 5:
	      tt = ttG ;
	      nT = tt ? ac_table_int (tt, 0, 0, 0) : 0 ;
	      if (nT > 10)
		{
		  float z = 100.0/nT ;
		  aceOutf (qc->ao, "\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f"
			   , z * ac_table_int (tt, 0, 2, 0)
			   , z * ac_table_int (tt, 0, 4, 0)
			   , z * ac_table_int (tt, 0, 6, 0) 
			   , z * ac_table_int (tt, 0, 8, 0)
			   , z * ac_table_int (tt, 0, 10, 0)			   
			   ) ;
		}
	      else
		aceOut (qc->ao, "\t\t\t\t\t") ;
	      break ;
	    case 6:
	      tt = ttG ;
	      nT = tt ? ac_table_int (tt, 0, 0, 0) : 0 ;
	      tt = ttP ;
	      if (nT > 10)
		{
		  float z = 100.0/nT ;
		  aceOutf (qc->ao, "\t%.2f\t%.2f\t%.2f\t%.2f"
			   , z * ac_table_int (tt, 0, 0, 0)
			   , z * ac_table_int (tt, 0, 2, 0)
			   , z * (ac_table_int (tt, 0, 4, 0) +  ac_table_int (tt, 0, 6, 0) + ac_table_int (tt, 0, 8, 0))
			   , z * ac_table_int (tt, 0, 10, 0)
			   ) ;
		}
	      else
		aceOut (qc->ao, "\t\t\t\t") ;
	      break ;
#endif
	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  
  ac_free (h) ;
  return;
} /*  qcSnpCoding */

/*************************************************************************************/
/*************************************************************************************/

static void qcDoPair (QC *qc, RC *rc, BOOL doPercent)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  const char *ccp ;
  int ir ;
  float z ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    /*     { "Aligned_fragments", "Number of aligned_fragments", 0, 0, 0} , */
    { "Compute", "%s aligned pairs\t%s compatible pairs\t%s incompatible pairs\t%s paired end fragments with only one read aligned\t%s links gene to distant genome\t%s links 2 genes\t%s incompatible topology\t%s links to rRNA\t%s links to mito", 1, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 
  const char *caption =
    "Paired reads compatibility:"
    "  In paired end sequencing, the 2 ends of each fragment are sequenced, yielding two reads. "
    "Considering all reads aligning at their best score, the paired-end  information is used to "
    "retain preferentially the target positions where both ends map in a consistent way, "
    "therefore decreasing the number of ambiguous mappings. It is also used to phase SNVs. "
    "The mapped fragments can then be partitioned into compatible pairs, non-compatible pairs  "
    "and pairs with a single read aligned. In a compatible pair, the 2 reads face each other and span  "
    "a segment no longer than 3 times the median insert length. For RNA-seq fragments best mapping  "
    "on the genome, the limit for compatible pairs is arbitrarily set at 1 Mb, to allow for  "
    "occasional long introns. The compatible and incompatible pairs are then partitioned in subcategories, "
    "according to the targets and the types of incompatibility. Please see the details in the file "
    "RESULTS/Mapping/$MAGIC.pair_stats.txt. Compatible gene extensions refer to fragment for which one end "
    "maps into an annotated gene and the other in a compatible way on the nearby genome."
  ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	{
	  const char *ccp ;
	  ccp = doPercent ? "%" : "Number of" ;
	  aceOut (qc->ao, "\t") ;
	  aceOutf (qc->ao, (char *)ti->title, ccp, ccp, ccp, ccp, ccp, ccp, ccp, ccp) ;
	}
     else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
     else if (! strcmp (ti->tag, "Compute"))
	{
	  switch (ti->col)
	    {
	    case 1:
	      if (! ac_has_tag (rc->run, "Single_end") && ac_has_tag (rc->ali, "Pair_fate"))
		{
		  float za, zb, zc, zd, ze ;

		  ccp = EMPTY ; tt = ac_tag_table (rc->ali, "Orphans", h) ;
		  z = 0 ;
		  for (ir = 0 ; tt && ir < tt->rows ; ir++)
		    {
		      ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		      if (! strcasecmp (ccp, "Any"))
			{
			  z = ac_table_float (tt, ir, 1, 0) ;
			}
		    } 
		  za = ac_tag_float (rc->ali, "Links_gene_to_distant_genome", 0) ;
		  zb = ac_tag_float (rc->ali, "Links_2_genes", 0) ;
		  zc = ac_tag_float (rc->ali, "Incompatible_topology", 0) ;
		  zd = ac_tag_float (rc->ali, "Links_to_rRNA", 0) ;
		  ze = ac_tag_float (rc->ali, "Links_to_mito", 0) ;
		  rc->var[T_cFrag] = ac_tag_float (rc->ali, "Compatible_pairs", 0) ;
		  rc->var[T_uFrag] = ac_tag_float (rc->ali, "Non_compatible_pairs", 0) ;
		  rc->var[T_aliFrag] = rc->var[T_cFrag] + rc->var[T_uFrag] ; /* aligned pairs, not counting orphans */
		  if (! doPercent)
		    {
		      if ( rc->var[T_aliFrag]) aceOutf (qc->ao, "\t%.2f", rc->var[T_aliFrag]) ; else aceOutf (qc->ao, "\t0") ;
		      if ( rc->var[T_cFrag]) aceOutf (qc->ao, "\t%.2f", rc->var[T_cFrag]) ; else aceOutf (qc->ao, "\t0") ;
		      if ( rc->var[T_uFrag]) aceOutf (qc->ao, "\t%.2f", rc->var[T_uFrag]) ; else aceOutf (qc->ao, "\t0") ;
		      aceOutf (qc->ao, "\t%.0f", z) ;
		      aceOutf (qc->ao, "\t%.0f", za) ;
		      aceOutf (qc->ao, "\t%.0f", zb) ;
		      aceOutf (qc->ao, "\t%.0f", zc) ;
		      aceOutf (qc->ao, "\t%.0f", zd) ;
		      aceOutf (qc->ao, "\t%.0f", ze) ;

		    }
		  else
		    {
		      aceOutf (qc->ao, "\t%.2f", rc->var[T_aliFrag] ? 200 * rc->var[T_aliFrag]/ rc->var[T_Read] : 0) ;	  
		      aceOutf (qc->ao, "\t%.2f", rc->var[T_aliFrag] ? 100 * rc->var[T_cFrag]/ rc->var[T_aliFrag] : 0) ;	  
		      aceOutf (qc->ao, "\t%.2f", rc->var[T_aliFrag] ? 100 * rc->var[T_uFrag] / rc->var[T_aliFrag] : 0) ;	  
		      aceOutf (qc->ao, "\t%.2f", rc->var[T_aliFrag] ? 100 * z / rc->var[T_aliFrag] : 0) ;	  
		      aceOutf (qc->ao, "\t%.2f", rc->var[T_aliFrag] ? 100 * za / rc->var[T_aliFrag] : 0) ;	  
		      aceOutf (qc->ao, "\t%.2f", rc->var[T_aliFrag] ? 100 * zb / rc->var[T_aliFrag] : 0) ;	  
		      aceOutf (qc->ao, "\t%.2f", rc->var[T_aliFrag] ? 100 * zc / rc->var[T_aliFrag] : 0) ;	  
		      aceOutf (qc->ao, "\t%.2f", rc->var[T_aliFrag] ? 100 * zd / rc->var[T_aliFrag] : 0) ;	  
		      aceOutf (qc->ao, "\t%.2f", rc->var[T_aliFrag] ? 100 * ze / rc->var[T_aliFrag] : 0) ;	  
		    }
		  break ;
		}
	      else
		aceOutf (qc->ao, "\t\t\t\t\t\t\t\t") ;
	    }
	}
      else
	{
	  qcShowTag (qc, rc, ti) ;
	}
    }
  ac_free (h) ;
  return;
}  /* qcDoPair */

static void qcPairNb  (QC *qc, RC *rc)
{
  return qcDoPair (qc, rc, FALSE) ;
}
static void qcPairPercent  (QC *qc, RC *rc)
{
  return qcDoPair (qc, rc, TRUE) ;
}

/*************************************************************************************/
/* template for a future chapter */
static void qcInsertSize (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Fragment_length_1", "1% fragments in library are shorter than (nt)", 0, 0, 0} ,
    { "Fragment_length_5", "5% fragments in library are shorter than (nt)", 0, 0, 0} ,
    { "Fragment_length_mode", "mode of fragment lengths", 0, 0, 0} ,
    { "Fragment_length_median", "median fragment length", 0, 0, 0} ,
    { "Fragment_length_average", "average of fragment lengths", 0, 0, 0} ,
    { "Fragment_length_95", "5% fragments in library are longer than (nt)", 0, 0, 0} ,
    { "Fragment_length_99", "1% fragments in library are longer than (nt)", 0, 0, 0} ,
    {  0, 0, 0, 0, 0}
  }; 

  const char *caption =
    "Insert size distribution measured from paired alignments, or read lengths in single end protocols"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;
  
  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	aceOutf (qc->ao, "\t%s", ti->title) ;
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (ti->tag[0] == 'F')
	{
	  if (ac_tag_int (rc->ali, "Fragment_length_average", 0))
	    qcShowTag (qc, rc, ti) ;
	  else if (ac_has_tag (rc->ali, "Length_distribution_1_5_50_95_99_mode_av"))
	    {
	      AC_TABLE tt = ac_tag_table (rc->ali, "Length_distribution_1_5_50_95_99_mode_av", h) ;
	      aceOutf (qc->ao, "\t%d", ac_table_int (tt, 0,0,0)) ;
	      aceOutf (qc->ao, "\t%d", ac_table_int (tt, 0,1,0)) ;
	      aceOutf (qc->ao, "\t%d", ac_table_int (tt, 0,5,0)) ;
	      aceOutf (qc->ao, "\t%d", ac_table_int (tt, 0,2,0)) ;
	      aceOutf (qc->ao, "\t%d", ac_table_int (tt, 0,6,0)) ;
	      aceOutf (qc->ao, "\t%d", ac_table_int (tt, 0,3,0)) ;
	      aceOutf (qc->ao, "\t%d", ac_table_int (tt, 0,4,0)) ;
	      break ;
	    }
	  else
	    aceOutf (qc->ao, "\t") ;
	}
      else
	aceOutf (qc->ao, "\t") ;
    }
   ac_free (h) ;
  return;
}  /* qcInsertSize */

/*************************************************************************************/
/* template for a future chapter */
static void qcInsertSizeHisto (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute","", 1, 0, 0} ,
    {  0, 0, 0, 0, 0}
  }; 

  const char *caption =
    "Histogram of insert length (bp) measured as the distance between the first aligned base of the /1 and /2 reads in  compatible pairs"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;
  
  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	aceOutf (qc->ao, "\t%s", ti->title) ;
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
       else  if (! strcmp (ti->tag, "Compute"))
	{
	  switch (ti->col)
	    {
	    default:
	      break ;
	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }

   ac_free (h) ;
  return;
}  /* qcInsertSizeHisto */

/*************************************************************************************/

static float qcShowGeneExp (QC *qc, AC_TABLE tt, const char *target, const char *tag)
{
  int ir ;
  const char *target2 = target ;
  if (qc->TT)
    target2 = hprintf (qc->h, "Captured_%s", target) ;

  if (tt)
    for (ir = 0 ; ir < tt->rows ; ir++)
      if (! strcasecmp (tag   , ac_table_printable (tt, ir, 0, "xxx")) &&
	  (
	   ! strcasecmp (tag, "Targeted_genes")  ||
	   ! strcasecmp (target2, ac_table_printable (tt, ir, 1, "xxx"))
	   )
	  )
	{
	  float z = ac_table_float (tt, ir, 2, 0) + ac_table_int (tt, ir, 2, 0) ;
	  if (z) /* a float */
	    aceOutf (qc->ao, "\t%.1f", z) ;
	  else
	    aceOutf (qc->ao, "\t%s", ac_table_printable (tt, ir, 2, "")) ;
	  return z ;
	}
  aceOut (qc->ao, "\t") ;
  return 0 ;
} /* qcShowGeneExp */

/*************************************************************************************/

static void qcProject (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  const char *ccp ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "Compute", "Project", 40, 0, 0} ,
    { "Compute", "Project abstract", 41, 0, 0} ,
    { "Compute", "Project description", 42, 0, 0} ,
    { "Author", "Laboratory, Author when available", 0, 0, 0} ,
    { "Compute", "Reference", 6, 0, 0} ,    
    { "Compute", "Design", 12, 0, 0} ,
    { "Compute", "Sequencing protocol", 13, 0, 0} ,
    { "Paired_end", "Single or Paired end", 0, 0, 0} ,
    { "Compute", "SRA date", 11, 0, 0} ,
    
   {  0, 0, 0, 0, 0}
  } ; 

  const char *caption =
    "Project"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  if (rc && rc->run && ! rc->srr) 
    rc->srr = ac_tag_obj (rc->run, "SRR", rc->h) ;
  if (rc && rc->run && ! rc->srx) 
    rc->srx = ac_tag_obj (rc->run, "SRX", rc->h) ;
  if (rc && rc->run && ! rc->srp) 
    rc->srp = ac_tag_obj (rc->run, "SRP", rc->h) ;
  if (rc && rc->run && ! rc->sample) 
    rc->sample = ac_tag_obj (rc->run, "Sample", rc->h) ;

  for (ti = tts ; ti->tag ; ti++)
    {
      int nw = 0 ;
      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title) ;
	}
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  DICT *papDict = 0 ;

	  aceOutf (qc->ao, "\t") ;
	  switch (ti->col)
	    {
	    case 1000:  /* RunId */
	      if (ac_has_tag (rc->run, "RunId"))
		tt = ac_tag_table (rc->run, "RunId", h) ;
	      else
		tt = ac_keyset_table (ac_objquery_keyset (rc->run, ">Union_of ; >RunId", h), 0, 0, 0, h) ;
	      if (tt && tt->rows) 
		{
		  if (tt->rows > 30)
		    aceOutf (qc->ao, "%d > 300 runs, list skipped", tt->rows) ;
		  else
		    qcShowMultiTag (qc, tt) ; 
		}
	      else
		aceOutf (qc->ao, "%s", ac_name (rc->run)) ;
	      break ;
	    case 1:  /* SRR */
	      tt = ac_tag_table (rc->run, "SRR", h) ;
	      if (! tt)
		tt = ac_tag_table (rc->run, "SRX", h) ;
	      if (! tt)
		tt = ac_tag_table (rc->run, "sublibraries", h) ;
	      if (! tt)
		aceOutf (qc->ao, EMPTY) ;
	      else 
		qcShowMultiTag (qc, tt) ;
	      break ;
	    case 2:  /* Author */
	      if (rc->srr)
		{
		  ccp = rc->run ? ac_tag_printable (rc->run, "Author", 0) : 0 ;
		  if (! ccp) 
		    ccp = rc->srr ? ac_tag_printable (rc->srr, "Magic_Author2", 0) : 0 ;
		  if (! ccp) 
		    ccp = rc->srr ? ac_tag_printable (rc->srr, "Author", EMPTY) : EMPTY ;
		  aceOutf (qc->ao, "%s", ccp) ;
		}
	      else
		{
		  int nAuthor = 0 ;
		  vTXT txt = vtxtHandleCreate (h) ;

		  ccp = rc->run ? ac_tag_printable (rc->run, "Author", 0) : 0 ;
		  if (ccp)
		    {
		      vtxtPrintf (txt, "%s", ccp) ;
		      nAuthor++ ;
		    }
		   ccp = rc->run ? ac_tag_printable (rc->run, "Sequencing_laboratory", 0) : 0 ;
		   if (ccp && (!nAuthor ||  ! strstr (vtxtPtr (txt), ccp)))
		    {
		      vtxtPrintf (txt, "%s%s"
				  , nAuthor++ ? ", " : ""
				  , ccp
				  ) ;
		    }
		   ccp = rc->run ? ac_tag_printable (rc->run, "Nucleic_prep_author", 0) : 0 ;
		   if (ccp && (!nAuthor ||  ! strstr (vtxtPtr (txt), ccp)))
		    {
		      vtxtPrintf (txt, "%s%s"
				  , nAuthor++ ? ", " : ""
				  , ccp
				  ) ;
		    }
		   if (nAuthor)
		     aceOutf (qc->ao, "%s", vtxtPtr (txt)) ;

		   txt = vtxtHandleCreate (h) ;
		   nAuthor = 0 ;
		   ccp = ac_tag_printable (rc->run, "Submission_date", 0) ;
		   if (ccp && (!nAuthor ||  ! strstr (vtxtPtr (txt), ccp)))
		     {
		       vtxtPrintf (txt, "%s%s"
				   , nAuthor++ ? ", " : ""
				   , ccp
				   ) ;
		     }		    
		   ccp = ac_tag_printable (rc->run, "Library_date", 0) ;
		   if (ccp && (!nAuthor ||  ! strstr (vtxtPtr (txt), ccp)))
		     {
		       vtxtPrintf (txt, "%s library prep %s"
				   , nAuthor++ ? ", " : ""
				   , ccp
				   ) ;
		     }		    
		   ccp = ac_tag_printable (rc->run, "Sequencing_date", 0) ;
		   if (ccp && (!nAuthor ||  ! strstr (vtxtPtr (txt), ccp)))
		     {
		       vtxtPrintf (txt, "%s Sequencing %s"
				   , nAuthor++ ? ", " : ""
				   , ccp
				   ) ;
		     }		    
		   if (nAuthor)
		     aceOutf (qc->ao, " [%s]", vtxtPtr (txt)) ;
	
		}	      
	      break ;
	    case 11:  /* Date */
	      {{
		  mytime_t d1 = ac_tag_date (rc->run, "Submission_date", 0 ) ;
		  mytime_t d2 = ac_tag_date (rc->run, "Release_date", 0) ;
		  char buf25[25] ;
		  if (d2 && d2 < d1) d1 = d2 ;
		  if (d1)
		    aceOutf (qc->ao, timeShow (d1, buf25, 25)) ;
		}}
	      break ;
	    case 3:  /* Species */
	      ccp = rc->run ? ac_tag_printable (rc->run, "Species", EMPTY) : EMPTY ;
	      aceOutf (qc->ao, "%s", ccp) ;
	      break ;
	    case 13:  /* Sequencing protocol */
	      tt = rc->srr ? ac_tag_table (rc->srr, "Nucleic_Acid_Extraction", h) : 0 ;
	      if (! tt)
		 tt = rc->srr ? ac_tag_table (rc->srr, "sraNucleic_Acid_Extraction", h) : 0 ;
	      if (tt)
		{
		  DICT *dict = dictHandleCreate (128,h) ;
		  int ir, jr ;
		  for (ir = 0 ; ir < tt->rows ; ir++)
		    for (jr = 0 ; jr < tt->cols ; jr++)
		      {
			ccp = ac_table_printable (tt, ir, jr, 0) ;
			if (ccp && ! dictFind (dict, ccp, 0))
			  {
			    dictAdd (dict, ccp, 0) ;
			    if (!strncasecmp (ccp, "sra", 3)) ccp += 3 ;
			    if (!strcasecmp (ccp, "sraUnspecified_RNA"))
			      continue ;
			    aceOutf (qc->ao, "%s%s", nw++ ? ", " : EMPTY, ccp) ;			  
			  }
		      }
		  ac_free (dict) ;
		}
	      ccp = rc->srr ? ac_tag_printable (rc->srr, "Annotated_strandedness", 0) : 0 ;
	      if (ccp) 
		aceOutf (qc->ao, ", %s", ccp) ;

	      break ;

	    case 12:  /* Design */
	      ccp = rc->srx ?  ac_tag_printable (rc->srx, "Design", 0) : 0 ;
	      if (ccp) 
		{ aceOutf (qc->ao, "%s ", ccp) ; nw++ ; }
	      ccp = rc->srx ?  ac_tag_printable (rc->srx, "Construction_protocol", 0) : 0 ;
	      if (ccp) 
		{ aceOut (qc->ao, ccp) ; nw++ ; }
	      if (! nw) 
		aceOut (qc->ao, EMPTY) ;
	      break ;
	      
	    case 40:  /* SRP-Project */
	      if (rc->srp) aceOutf (qc->ao, "%s", ac_name (rc->srp)) ;
	      break ;
	    case 41:  /* SRP-Project title */
	      ccp = rc->srp ? ac_tag_printable (rc->srp, "Title", EMPTY) : EMPTY ;
	      aceOutf (qc->ao, "%s", ccp) ;
	      break ;
	    case 42:  /* SRP-Project description */
	      ccp = rc->srp ? ac_tag_printable (rc->srp, "Abstract", EMPTY) : EMPTY ;
	      aceOutf (qc->ao, "%s", ccp) ;
	      break ;
	    case 6:  /* Reference */
	      tt = rc->run ? ac_tag_table (rc->run, "Reference", h) : 0 ;
	      if (!tt) tt = rc->srp ? ac_tag_table (rc->srp, "Reference", h) : 0 ; 
	      if (tt)
		{
		  int ir ;
		  for (ir = 0 ; ir < tt->rows ; ir++)
		    {
		      AC_OBJ pap = ac_table_obj (tt, ir, 0, h) ;
		      if (! papDict)
			papDict = dictHandleCreate (32, h) ;
		      if (! dictAdd (papDict, ac_name(pap), 0))
			continue ;			
		      aceOutf (qc->ao, "%sPMID: %s", nw ? "; " : EMPTY, ac_name(pap)+(!strncmp( ac_name(pap),"pm",2) ? 2 : 0)) ;
		      ccp = pap ? ac_tag_printable (pap, "Title", 0) : 0 ; 
		      if (ccp) 
			{ aceOutf (qc->ao, " : %s", ccp) ; nw++ ; }
		      ccp = pap ? ac_tag_printable (pap, "Citation", 0) : 0 ;
		      if (ccp) aceOutf (qc->ao, "%s%s", nw ? ", " : ": ", ccp) ;
		      ccp = pap ? ac_tag_printable (pap, "Author", 0) : 0 ;
		      if (ccp) 
			{
			  AC_TABLE tt1 = ac_tag_table (pap, "Author", h) ;
			  aceOutf (qc->ao, ", %s", ccp) ;
			  if (tt1 && tt1->rows > 1) aceOutf (qc->ao, " et al.") ;
			}
		    }
		}
	      tt = rc->srp ? ac_tag_table (rc->srp, "Reference", h) : 0 ; 
	      if (tt)
		{
		  int ir ;
		  for (ir = 0 ; ir < tt->rows ; ir++)
		    {
		      AC_OBJ pap = ac_table_obj (tt, ir, 0, h) ;
		      if (! papDict)
			papDict = dictHandleCreate (32, h) ;
		      if (! dictAdd (papDict, ac_name(pap), 0))
			continue ;			
		      aceOutf (qc->ao, "%sPMID: %s", nw ? "; " : EMPTY, ac_name(pap)+(!strncmp( ac_name(pap),"pm",2) ? 2 : 0)) ;
		      ccp = pap ? ac_tag_printable (pap, "Title", 0) : 0 ; 
		      if (ccp) 
			{ aceOutf (qc->ao, " : %s", ccp) ; nw++ ; }
		      ccp = pap ? ac_tag_printable (pap, "Citation", 0) : 0 ;
		      if (ccp) aceOutf (qc->ao, "%s%s", nw ? ", " : ": ", ccp) ;
		      ccp = pap ? ac_tag_printable (pap, "Author", 0) : 0 ;
		      if (ccp) 
			{
			  AC_TABLE tt1 = ac_tag_table (pap, "Author", h) ;
			  aceOutf (qc->ao, ", %s", ccp) ;
			  if (tt1 && tt1->rows > 1) aceOutf (qc->ao, " et al.") ;
			}
		    }
		}
	      tt = rc->srx ? ac_tag_table (rc->srx, "Reference", h) : 0 ; 
	      if (tt)
		{
		  int ir ;
		  for (ir = 0 ; ir < tt->rows ; ir++)
		    {
		      AC_OBJ pap = ac_table_obj (tt, ir, 0, h) ;
		      if (! papDict)
			papDict = dictHandleCreate (32, h) ;
		      if (! dictAdd (papDict, ac_name(pap), 0))
			continue ;			
		      
		      aceOutf (qc->ao, "%sPMID: %s", nw ? "; " : EMPTY, ac_name(pap)+(!strncmp( ac_name(pap),"pm",2) ? 2 : 0)) ;
		      ccp = pap ? ac_tag_printable (pap, "Title", 0) : 0 ; 
		      if (ccp) 
			{ aceOutf (qc->ao, " : %s", ccp) ; nw++ ; }
		      ccp = pap ? ac_tag_printable (pap, "Citation", 0) : 0 ;
		      if (ccp) aceOutf (qc->ao, "%s%s", nw ? ", " : ": ", ccp) ;
		      ccp = pap ? ac_tag_printable (pap, "Author", 0) : 0 ;
		      if (ccp) 
			{
			  AC_TABLE tt1 = ac_tag_table (pap, "Author", h) ;
			  aceOutf (qc->ao, ", %s", ccp) ;
			  if (tt1 && tt1->rows > 1) aceOutf (qc->ao, " et al.") ;
			}
		    }
		}
	      tt = rc->srr ? ac_tag_table (rc->srr, "Reference", h) : 0 ; 
	      if (tt)
		{
		  int ir ;
		  for (ir = 0 ; ir < tt->rows ; ir++)
		    {
		      AC_OBJ pap = ac_table_obj (tt, ir, 0, h) ;
		      if (! papDict)
			papDict = dictHandleCreate (32, h) ;
		      if (! dictAdd (papDict, ac_name(pap), 0))
			continue ;			
		      
		      aceOutf (qc->ao, "%sPMID: %s", nw ? "; " : EMPTY, ac_name(pap)+(!strncmp( ac_name(pap),"pm",2) ? 2 : 0)) ;
		      ccp = pap ? ac_tag_printable (pap, "Title", 0) : 0 ; 
		      if (ccp) 
			{ aceOutf (qc->ao, " : %s", ccp) ; nw++ ; }
		      ccp = pap ? ac_tag_printable (pap, "Citation", 0) : 0 ;
		      if (ccp) aceOutf (qc->ao, "%s%s", nw ? ", " : ": ", ccp) ;
		      ccp = pap ? ac_tag_printable (pap, "Author", 0) : 0 ;
		      if (ccp) 
			{
			  AC_TABLE tt1 = ac_tag_table (pap, "Author", h) ;
			  aceOutf (qc->ao, ", %s", ccp) ;
			  if (tt1 && tt1->rows > 1) aceOutf (qc->ao, " et al.") ;
			}
		    }
		}
	      if (! nw)
		aceOutf (qc->ao, EMPTY) ;
	      ac_free (papDict) ;
	      break ;

	    case 20:   /* sample */
	      if (rc->sample)
		{
		  ccp = ac_tag_printable (rc->sample, "Title", 0) ;
		  if (! ccp)
		    ccp = ac_name (rc->sample) ;
		}
	      else
		{
		  ccp = ac_tag_printable (rc->run, "Sample", 0) ;
		}
	      if (! ccp) ccp = EMPTY ;
	      aceOutf (qc->ao, "%s", ccp) ;
	      break ;
	    case 7:   /* first line */
	      ccp = ac_tag_printable (rc->run, "First_line", 0) ;
	      if (! ccp) ccp = EMPTY ;
	      aceOutf (qc->ao, "%s", ccp) ;
	      break ;
	    case 21:   /* stage */
	      tt = rc->sample? ac_tag_table (rc->sample, "Stage", h) : 0 ;
	      nw = qcShowTable (qc, tt, 0) ;
	      if (! nw)
		aceOutf (qc->ao, EMPTY) ;
	      break ;
	    case 22:   /* system */
	      tt = rc->sample? ac_tag_table (rc->sample, "Tissue", h) : 0 ;
	      nw = qcShowTable (qc, tt, 0) ;
	      tt = rc->sample? ac_tag_table (rc->sample, "Systm", h) : 0 ;
	      nw += qcShowTable (qc, tt, nw) ;
	      if (! nw)
		aceOutf (qc->ao, EMPTY) ;
              break ;
	    case 23: /* bio */
	      tt = rc->run ? ac_tag_table (rc->run, "CT", h) : 0 ;
	      nw = qcShowTable (qc, tt, nw) ;
	      tt = rc->run ? ac_tag_table (rc->run, "Bio", h) : 0 ;
	      nw += qcShowTable (qc, tt, nw) ;
	      if (! nw)
		aceOutf (qc->ao, EMPTY) ;
	      break ;
	      if (! nw)
		aceOutf (qc->ao, EMPTY) ;
	      break ;
	    case 24: /* sex */
	      tt = rc->sample? ac_tag_table (rc->sample, "Sex", h) : 0 ;
	      nw = qcShowTable (qc, tt, nw) ;
	      if (! nw)
		aceOutf (qc->ao, EMPTY) ;
	      break ;
	    }
	}
      else if (! strcmp (ti->tag, "Paired_end"))
	{ 
	  aceOutf (qc->ao, "\t") ;
	  if (ac_has_tag (rc->run, "Is_run"))
	    aceOutf (qc->ao, rc->Paired_end ? "Paired_end" : "Single-end") ;
	  else
	    {
	      int n1, n2 ;
	      n1 = ac_keyset_count (ac_objquery_keyset (rc->run, ">Union_of ;>Runs ;   Paired_end", h)) ;
	      n2 = ac_keyset_count (ac_objquery_keyset (rc->run, ">Union_of ; ! Paired_end", h)) ;
	      if (n1 > 0 && n2 > 0)
		aceOutf (qc->ao, "Mix") ;
	      else if (n1 > 0 && n2 == 0)
		aceOutf (qc->ao, "Paired_end") ;
	      else if (n1 == 0 && n2 > 0)
		aceOutf (qc->ao, "Single-end") ;
	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  ac_free (h) ;
  return;
}  /* qcProject */

/*************************************************************************************/

static void qcMainResults1 (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  const char *ccp ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,  
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute", "Million raw reads", 20, 0, 0} ,  /* copied from qcBeforeAli */
    { "Compute", !qc->TT ? "Million reads aligned on any target by Magic" : "Million reads aligned on targeted genes" , 30, 0, 0 } , /* copied from qcAvLengthAli */
    { "Compute", "Megabases sequenced", 57, 0, 0} ,
    { "Compute", "Mb aligned on any target", 52, 0, 0 } , /* copied from qcAvLengthAli */
    { "Compute", "Average fragment multiplicity, high numbers reflect PCR amplification or very short reads", 40, 0, 0} ,  /* copied from qcBeforeAli */
    { "Compute", "Average read length after clipping adaptors and barcodes (nt)", 50, 0, 0 } , /* copied from qcAvLengthAli */
    { "Compute", "Average length aligned per read (nt)", 60, 0, 0 } , /* copied from qcAvLengthAli */
    { "Compute", "Average length aligned per fragment (nt)", 61, 0, 0 } , /* copied from qcAvLengthAli */
    {  0, 0, 0, 0, 0}
  } ; 
  
  const char *caption =
    "Reads, Megabases, multiplicity and length"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;
  
  char *myAny = qc->TT ? "Captured_av" : "any" ;
  float zr = 0, zs = 0 ;

  for (ti = tts ; ti->tag ; ti++)
    {
      const char *tag ;
      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
       else if (! strcmp (ti->tag, "Compute"))
	{
	  int ir ;
	  float z, zb, zc ;
	  
	  switch (ti->col)
	    {
	    case 20: /* Million raw reads */
	      qcShowMillions (qc, rc->var[T_Read]) ;
	      break ;
	    case 30: /* Million reads aligned on any target by magic */
	      ccp = EMPTY ;
	      tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
 	      z = zb = zc = zr = zs = 0 ;
 	      if (tt)
		for (ir = 0 ; tt && ir < tt->rows ; ir++)
		  {
		    ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		    if (! strcasecmp (ccp, myAny))
		      {
			zs = ac_table_float (tt, ir, 1, 0) ;
			zr = ac_table_float (tt, ir, 3, 0) ;
			zb = ac_table_float (tt, ir, 5, 0) ;
			zc = ac_table_float (tt, ir, 7, 0) ;
		      }
		  } 
	      aceOutf (qc->ao, "\t%.3f", zr/1000000) ;
	      ac_free (tt) ;
	      break ;
	      
	    case 40: /* Average fragment multiplicity */
	      if (qc->TT)
		{
		  aceOutf (qc->ao, "\t%.2f", zs ? zr/zs : 0) ;
		}
	      else
		aceOutf (qc->ao, "\t%.2f", rc->var[T_Seq] ? rc->var[T_Read]/rc->var[T_Seq] : 0) ;
	      break ;
	    case 50:
	      tag = "nh_Ali" ; 
	      tt = ac_tag_table (rc->ali, tag, h) ;
	      z = 0 ;
	      aceOut (qc->ao, "\t") ;
	      if (tt)
		for (ir = 0 ; tt &&  ir < tt->rows; ir++)
		  {
		    if (! strcasecmp (ac_table_printable (tt, ir, 0, ""), qc->TT ? "targeted_genes" : myAny))
		      {
			z = ac_table_float (tt, ir, 11, 0) ;
			aceOutf (qc->ao, "%.2f", z) ;
		      }
		   }
	      ac_free (tt) ;
	      break ;
	      
	    case 52: /* Million reads aligned on any target by magic */
	      ccp = EMPTY ;
	      tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
 	      z = zb = zc = 0 ;
 	      if (tt)
		for (ir = 0 ; tt && ir < tt->rows ; ir++)
		  {
		    ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		    if (! strcasecmp (ccp, "any"))
		      {
			z = ac_table_float (tt, ir, 3, 0) ;
			zb = ac_table_float (tt, ir, 5, 0) ;
			zc = ac_table_float (tt, ir, 7, 0) ;
		      }
		  } 
	      aceOutf (qc->ao, "\t%.3f", zb/1000) ; /* Mb aligned on any target */
	      ac_free (tt) ;
	      break ;

	    case 57:
	      aceOutf (qc->ao, "\t%.3f", rc->var[T_kb]/1000) ;   
	      break ;

	    case 60: /* Average length aligned per read */
	      ccp = EMPTY ;
	      tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
 	      z = zb = zc = 0 ;
 	      if (tt)
		for (ir = 0 ; tt && ir < tt->rows ; ir++)
		  {
		    ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		    if (! strcasecmp (ccp, myAny))
		      {
			z = ac_table_float (tt, ir, 3, 0) ;
			zb = ac_table_float (tt, ir, 5, 0) ;
			zc = ac_table_float (tt, ir, 7, 0) ;
		      }
		  } 
	      aceOutf (qc->ao, "\t%.2f", zc) ;  /* average length aligned */
	      ac_free (tt) ;
	      break ;

	    case 61: /* Average length aligned per fragment */
	      ccp = EMPTY ;
	      tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
 	      z = zb = zc = 0 ;
 	      if (tt)
		for (ir = 0 ; tt && ir < tt->rows ; ir++)
		  {
		    ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		    if (! strcasecmp (ccp, "any1") || ! strcasecmp (ccp, "any2"))
		      {
			zb += ac_table_float (tt, ir, 5, 0) ;
		      }
		  }
	      zc = ac_tag_float (rc->ali,  "Aligned_fragments", 0) ;
	      if (zc > 0) z = 1000.0 * zb / zc ;
	      aceOutf (qc->ao, "\t%.2f", z) ;  /* average length aligned per fragment */
	      ac_free (tt) ;
	      break ;
	    }     
	}
    }
  ac_free (h) ;
  return;
}  /* qcMainResults1 */

/*************************************************************************************/

static void qcMainResults2 (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  const char *ccp ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,  
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute", "% Mb aligned on any target", 54, 0, 0} , 
    { "Compute", "% length aligned after clipping adaptors and barcodes", 55, 0, 0 } ,
    { "Compute", "% Reads aligned on any target", 56, 0, 0 } ,

     {  0, 0, 0, 0, 0}
  } ; 

  const char *caption =
    "Alignment statistics"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;


  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title) ;
	}
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  int ir ;
	  float z, zb, zc ;

	  switch (ti->col)
	    {
	    case 54:
	      ccp = EMPTY ; tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
 	      z = zb = zc = 0 ;
 	      for (ir = 0 ; tt && ir < tt->rows ; ir++)
 		{
 		  ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
 		  if (! strcasecmp (ccp, "any"))
 		    {
 		      z = ac_table_float (tt, ir, 3, 0) ;
 		      zb = ac_table_float (tt, ir, 5, 0) ;
 		      zc = ac_table_float (tt, ir, 7, 0) ;
 		    }
 		} 

	      aceOutf (qc->ao, "\t%.2f", rc->var[T_kb] ? 100 *zb / rc->var[T_kb] : 0) ; /* % Mb aligned on any target after clipping */
	      ac_free (tt) ;
	      break ;

	    case 55:   /*  Average clipped length */
	      {
		int ir ;
		float z = -1, zClipped = -1 ;
		AC_TABLE tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
		
		for (ir = 0 ; tt &&  ir < tt->rows; ir++)
		  {
		    if (! strcasecmp (ac_table_printable (tt, ir, 0, ""), "any"))
		      {
			z =  ac_table_float (tt, ir, 7, -1) ; 
			zClipped = ac_table_float (tt, ir, 11, -1) ;
			break ;
		      }
		  }
		aceOutf (qc->ao, "\t") ;
		if (z > -1)
		  {
		    aceOutPercent (qc->ao, 100.00 * z/zClipped) ;
		  } 
		
		ac_free (tt) ;
	      }
	      break ;

	    case 56:
	      ccp = EMPTY ; 
	      tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
 	      z = zb = zc = 0 ;
 	      for (ir = 0 ; tt && ir < tt->rows ; ir++)
 		{
 		  ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
 		  if (! strcasecmp (ccp, "any"))
 		    {
 		      z = ac_table_float (tt, ir, 3, 0) ;
 		    }
 		} 

	      aceOutf (qc->ao, "\t%.2f", rc->var[T_Read] ? 100 *z / rc->var[T_Read] : 0) ;
	      ac_free (tt) ;
	      break ;
	    }
	}
      else  
	qcShowTag (qc, rc, ti) ;
    }
  ac_free (h) ;
  return;
}  /* qcMainResults2 */

/*************************************************************************************/

static void qcMainResults3 (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,  /* Quality features */
    { "TITLE", "Title", 10, 0, 0} ,

    { "Compute", "% fragments on strand plus of genes", 90, 0, 0} ,  /* copied from qcGeneExpression */
    { "Compute", "Mismatches per kb aligned", 100, 0, 0} , /* copied from qcMismatchTypes */
    { "Accessible_length_8kb", "Homogeneous coverage of transcripts longer than 8kb, max 6kb, limited by the 3' bias in poly-A selected experiments", 0, 0, 0}, /* copied from qc3pBias */
    { "Accessible_length_5kb", "Homogeneous coverage of transcripts longer than 5kb, max 4kb, limited by the 3' bias in poly-A selected experiments", 0, 0, 0}, /* copied from qc3pBias */

    {  0, 0, 0, 0, 0}
  } ; 

  const char *caption =
    "Quality features: strandedness, base quality, 3' bias"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;


  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title) ;
	}
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  int ir ;
	  float z, z1 ;
	  char buf[64] ;

	  switch (ti->col)
	    {

	    case 90:  /* Observed strandedness */
	      tt = ac_tag_table (rc->run, "Observed_strandedness_in_ns_mapping", h) ; 
	      z = 50 ;
	      if (tt)
		for (ir = 0 ; ir < tt->rows ; ir++)
		  {
		    z1 = ac_table_float (tt, 0,1,50) ;
		    if ((z-50)*(z-50) < (z1-50)*(z1-50))
		      z = z1 ;
		  }
	      aceOutf (qc->ao, "\t%.2f", z) ;
	      ac_free (tt) ;
	      break ;
	    case 100:	
	      tt = ac_tag_table (rc->ali, "Error_profile", h) ;	
	      z = z1 = 0 ; buf[0] = 0 ;
	      for (ir = 0 ; tt &&  ir < tt->rows ; ir++)
		{
		  if (ir == 0)
		    {
		      strcpy (buf, ac_table_printable (tt, ir, 0, "xxx")) ;
		      z1 = ac_table_float (tt, ir, 3, 0) ;
		    }
		  else if (strcmp (buf, ac_table_printable (tt, ir, 0, "xxx")))
		    {
		      strcpy (buf, ac_table_printable (tt, ir, 0, "xxx")) ;
		      z1 +=  ac_table_float (tt, ir, 3, 0) ;
		    }
		  if (!strcasecmp ("Any", ac_table_printable (tt, ir, 1, "xxx")))
		    z += ac_table_float (tt, ir, 2, 0) ;
		}
	      aceOutf (qc->ao, "\t%.5f", z1 > 0 ? z / (1000 * z1) : 0) ;
	      ac_free (tt) ;
	      break ;
	    }
	}
      else  /* 110: accessible length 8k/5k */
	qcShowTag (qc, rc, ti) ;
    }
  ac_free (h) ;
  return;
}  /* qcMainResults3 */

/*************************************************************************************/

static void qcMainResults4 (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  const char *ccp ;
  float zz ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,  
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute", "Megabases sequenced", 57, 0, 0} ,
    { "Compute", "Mb aligned on any target", 52, 0, 0 } , /* copied from qcAvLengthAli */
    { "Compute", "Megabases uniquely aligned and used for mismatch counts", 53, 0, 0} ,
    { "Compute", "Mb uniquely and fully aligned in genes", 16, 0, 0 },
    { "Compute", "Mb uniquely aligned in main protein coding genes minus very highly expressed genes", 80, 0, 0 }, /* copied from qcGeneExpression */

    {  0, 0, 0, 0, 0}
  } ; 

  const char *caption =
    "Megabases aligned"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;


  for (ti = tts ; ti->tag ; ti++)
    {
      const char *tag, *target = qc->Etargets[0] ? qc->Etargets[0] : "xx" ;
      if (rc == 0)
	{
	  switch (ti->col)
	     {
	     case 1:
	       if (! strcasecmp (target, "av")) target = "AceView" ;
	       aceOutf (qc->ao, "\t%s %s", target, ti->title) ;
	       break ;
	     default:
	       aceOutf (qc->ao, "\t%s", ti->title) ;
	       break ;
	     }
	}
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  int ir ;
	  float z, zb, zc ;
	  char buf[64] ;

	  switch (ti->col)
	    {
	    case 20: /* Million raw reads */
	      qcShowMillions (qc, rc->var[T_Read]) ;
	      break ;
	    case 30: /* Average fragment multiplicity */
	      aceOutf (qc->ao, "\t%.2f", rc->var[T_Seq] ? rc->var[T_Read]/rc->var[T_Seq] : 0) ;
	      break ;
	    case 40:
	       tag = "nh_Ali" ; 
	       tt = ac_tag_table (rc->ali, tag, h) ;
	       z = 0 ;
	       aceOut (qc->ao, "\t") ;
	       if (tt)
		 for (ir = 0 ; tt &&  ir < tt->rows; ir++)
		   {
		     if (! strcasecmp (ac_table_printable (tt, ir, 0, ""), "any"))
		       {
			 z = ac_table_float (tt, ir, 11, 0) ;
			 aceOutf (qc->ao, "%.2f", z) ;
		       }
		   }
	       ac_free (tt) ;
	       break ;

	    case 16:
	      { 
		const char *target = qc->Etargets[0] ? qc->Etargets[0] : "xx" ;
		tag = "Mb_in_genes" ;
		tt = ac_tag_table (rc->ali, tag, h) ;
		if (tt)
		  {
		    int ii, ir ;
		    for (ii = 0 ; ii < (ti->col > 10 ? 1 : 3) ; ii++)
		      {
			for (ir = 0 ; ir < tt->rows ; ir++)
			  {
			    target = qc->Etargets [ii] ;
			    if (target && ! strcasecmp (target, ac_table_printable (tt, ir, 0, "x")))
			      break ;
			  }
			if (ir < tt->rows)
			  {
			    switch (ti->col)	
			      {
			      case 16:
				rc->geneSigExpressed = ac_table_float (tt, ir, 1, 0) ;
				aceOutf (qc->ao, "\t%.3f", ac_table_float (tt, ir, 1, 0)) ;
				break ;
			      }
			  }
			else
			  aceOutf (qc->ao, "\t") ;
		      }
		  }
		else
		  aceOutf (qc->ao, "\t") ;
		break ;
	      }
	    case 50: /* Million reads aligned on any target by magic */
	      ccp = EMPTY ;
	      tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
 	      z = zb = zc = 0 ;
 	      if (tt)
		for (ir = 0 ; tt && ir < tt->rows ; ir++)
		  {
		    ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		    if (! strcasecmp (ccp, "any"))
		      {
			z = ac_table_float (tt, ir, 3, 0) ;
			zb = ac_table_float (tt, ir, 5, 0) ;
			zc = ac_table_float (tt, ir, 7, 0) ;
		      }
		  } 
	      aceOutf (qc->ao, "\t%.3f", z/1000000) ;
	      ac_free (tt) ;
	      break ;

	    case 51: /* Million reads aligned on any target by magic */
	      ccp = EMPTY ;
	      tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
 	      z = zb = zc = 0 ;
 	      if (tt)
		for (ir = 0 ; tt && ir < tt->rows ; ir++)
		  {
		    ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		    if (! strcasecmp (ccp, "any"))
		      {
			z = ac_table_float (tt, ir, 3, 0) ;
			zb = ac_table_float (tt, ir, 5, 0) ;
			zc = ac_table_float (tt, ir, 7, 0) ;
		      }
		  } 
	      aceOutf (qc->ao, "\t%.2f", zc) ;  /* average length aligned */
	      ac_free (tt) ;
	      break ;

	    case 52: /* Million reads aligned on any target by magic */
	      ccp = EMPTY ;
	      tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
 	      z = zb = zc = 0 ;
 	      if (tt)
		for (ir = 0 ; tt && ir < tt->rows ; ir++)
		  {
		    ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		    if (! strcasecmp (ccp, "any"))
		      {
			z = ac_table_float (tt, ir, 3, 0) ;
			zb = ac_table_float (tt, ir, 5, 0) ;
			zc = ac_table_float (tt, ir, 7, 0) ;
		      }
		  } 
	      aceOutf (qc->ao, "\t%.3f", zb/1000) ; /* Mb aligned on any target */
	      ac_free (tt) ;
	      break ;

	    case 53:   /* Megabases uniquely aligned and used for mismatch counts */
	      tt = ac_tag_table (rc->ali, "Error_profile", h) ;
	      for (ir = 0 ; tt &&  ir < tt->rows ; ir++)
		{
		  if (ir == 0)
		    {
		      strcpy (buf, ac_table_printable (tt, ir, 0, "xxx")) ;
		      zz = ac_table_float (tt, ir, 3, 0) ;
		    }
		  else if (strcmp (buf, ac_table_printable (tt, ir, 0, "xxx")))
		    {
		      strcpy (buf, ac_table_printable (tt, ir, 0, "xxx")) ;
		      zz +=  ac_table_float (tt, ir, 3, 0) ;
		    }
		}
	      aceOutf (qc->ao, "\t%.0f", zz) ;
	      ac_free (tt) ;
	      break ;

	    case 57:
	      aceOutf (qc->ao, "\t%.3f", rc->var[T_kb]/1000) ;   
	      break ;

	    case 80:
	      tag = "Mb_in_genes_with_GeneId_minus_high_genes" ;
	      aceOut (qc->ao, "\t") ;
	      tt = ac_tag_table (rc->ali, tag, h) ;
	      if (tt)
		{
		  for (ir = 0 ; ir < tt->rows ; ir++)
		    {
		      if (target && ! strcasecmp (target, ac_table_printable (tt, ir, 0, "x")))
			break ;
		    }
		  if (ir < tt->rows)
		    {
		      switch (ti->col)	
			{
			case 80:   
			  aceOutf (qc->ao, "%.03f", ac_table_float (tt, ir, 1, 0)) ;
			  break ;
			}
		    }
		}
	      ac_free (tt) ;
	      break ;
	    }
	}
      else 
	qcShowTag (qc, rc, ti) ;
    }
  ac_free (h) ;
  return;
}  /* qcMainResults4 */

/*************************************************************************************/

static void qcMainResults5 (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,  /* Quality features */
    { "TITLE", "Title", 10, 0, 0} ,
    { "Ends", "Supported known starts", 41, 0, 0} , 
    { "Ends", "Candidate starts", 42, 0, 0} ,
    { "Ends", "Supported known ends", 43, 0, 0} ,
    { "Ends", "Candidate ends", 44, 0, 0} ,
    { "Compute", "Supported known exon-exon junctions", 10, 0, 0} ,
    { "Compute", "Candidate new exon-exon junctions", 20, 0, 0} ,
    { "Compute", "Supports of known exon-exon junctions", 30, 0, 0} ,
    /*
    { "Compute2", "Supports of known exon-exon junctions per Megabase uniquely aligned in genes", 40, 0, 0} ,
    { "Compute2", "AceView genes significantly expressed per Megabase uniquely aligned in main protein coding genes", 50, 0, 0} ,
    */
    { "Compute2", "Supports of known exon-exon junction per kilobase uniquely aligned in genes", 60, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 

  const char *caption =
    "Transcripts starts, transcripts ends, and exon junctions"
    ;
  float z10 = 0, z30 = 0 ;
  float zg = -1 ;

  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;


  for (ti = tts ; ti->tag ; ti++)
    {
      const char *tag, *target = qc->Etargets[0] ? qc->Etargets[0] : "xx", *target2 = 0 ;
      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title) ;
	}
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  int ic, ir ;
	  float z = 0 ;
	  target2 = target ;
	  switch (ti->col)
	    {
	    case 10: tag = "Candidate_introns" ;  ic = 3 ; target2 = "AceView" ;  break ;
	    case 20: tag = "Candidate_introns" ;  ic = 5 ; target2 = "AceView" ;  break ;
	    case 30: tag = "Candidate_introns" ;  ic = 13 ; target2 = "AceView" ;  break ;
	    }
	  tt = ac_tag_table (rc->ali, tag, h) ;
	  z = 0 ;
	  if (tt)
	    for (ir = 0 ; ir < tt->rows ; ir++)
	      {
		if (!ir || ! strcasecmp (ac_table_printable (tt, ir, 0, "xxx"), target2))
		  z = ac_table_int (tt, ir, ic, 0) + ac_table_float (tt, ir, ic, 0) ;
	      }
	  aceOutf (qc->ao, "\t%.0f", z) ;
	  if (ti->col == 10)
	    z10 = z ;
	  if (ti->col == 30)
	    z30 = z ;
	}
      else if (! strcmp (ti->tag, "Ends"))
	{
	  int k = 0 ;
	  switch (ti->col)
	    {
	    case 41: tag = "Confirmed_starts" ; break ;
	    case 42: tag = "Candidate_starts" ; break ;
	    case 43: tag = "Confirmed_ends" ; break ;
	    case 44: tag = "Candidate_ends" ; break ;
	    }
	  k = ac_tag_int (rc->ali, tag, 0) ;
	  aceOutf (qc->ao, "\t%d", k) ;
	}
      else if (! strcmp (ti->tag, "Compute2"))
	{
	  float z = 0 ;
	  if (zg == -1)
	    {
	      const char *tag = "Mb_in_genes_with_GeneId_minus_high_genes" ;
	      tt = ac_tag_table (rc->ali, tag, h) ;
	      if (tt)
		{
		  int ir ;
		  const char *target = qc->Etargets [0] ; 
		  for (ir = 0 ; ir < tt->rows ; ir++)
		    {
		      if (target && ! strcasecmp (target, ac_table_printable (tt, ir, 0, "x")))
			break ;
		    }
		  if (ir < tt->rows)
		    {
		      zg = ac_table_float (tt, ir, 1, 0) ;
		    }
		}
	      ac_free (tt) ;
	    }
	  switch (ti->col)
	    {
	    case 40: z = zg > 0 ? z10/zg : 0 ; break ;
	    case 50: z = zg > 0 ? rc->geneSigExpressed/zg : 0 ; break ;
	    case 60: z = 0.001 * (zg > 0 ? z30/zg : 0) ; break ;
	    }
	  aceOutf (qc->ao, "\t%.2f", z) ;
	}
      else  
	qcShowTag (qc, rc, ti) ;
    }
  ac_free (h) ;
  return;
}  /* qcMainResults5 */

/*************************************************************************************/

static void qcSortingTitles (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    {  "Sorting_title", "Sorting_title 1", 0, 0, 0} ,
    {  "Sorting_Title_2", "Sorting title 2", 0, 0, 0} ,  
    {  0, 0, 0, 0, 0}
  } ; 

  const char *caption =
    "Sorting titles"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;
  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title) ;
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  ac_free (h) ;
  return;
}  /* qcSortingTitles */

/*************************************************************************************/

static void qcTitles (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  const char *ccp ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute", "RunId ", 1000, 0, 0} ,
    { "Compute", "Run accession in SRA", 1, 0, 0} ,
    { "Compute", "First line of the fastq file", 7, 0, 0} ,

    {  "Other_title", "Other title", 0, 0, 0} ,
    {  "Sorting_title", "Sorting_title 1", 0, 0, 0} ,
    {  "Sorting_Title_2", "Sorting title 2", 0, 0, 0} ,  
    {  "ERROR", "Download error", 0, 0, 0} ,
    /*    { "Compute", "Sample (summarized from biosample or manual)", 20, 0, 0} , */
    { "Compute", "Stage", 21, 0, 0} , 
    { "Compute", "Tissue, System or cell type", 22, 0, 0} ,
    { "Compute", "Bio: control or treatment, genotype, number of replicates", 23, 0, 0} ,
    { "Compute", "Sex", 24, 0, 0} ,
    { "Compute", "Species", 3, 0, 0} ,
   {  0, 0, 0, 0, 0}
  } ; 

  const char *caption =
    "Metadata for sample, project and laboratory"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {
      int nw = 0 ;
      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title) ;
	}
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  DICT *papDict = 0 ;

	  aceOutf (qc->ao, "\t") ;
	  switch (ti->col)
	    {
	    case 1000:  /* RunId */
	      if (ac_has_tag (rc->run, "RunId"))
		tt = ac_tag_table (rc->run, "RunId", h) ;
	      else
		tt = ac_keyset_table (ac_objquery_keyset (rc->run, ">Union_of ; >RunId", h), 0, 0, 0, h) ;
	      if (tt && tt->rows) 
		{
		  if (tt->rows > 30)
		    aceOutf (qc->ao, "%d > 300 runs, list skipped", tt->rows) ;
		  else
		    qcShowMultiTag (qc, tt) ; 
		}
	      else
		aceOutf (qc->ao, "%s", ac_name (rc->run)) ;
	      break ;
	    case 1:  /* SRR */
	      tt = ac_tag_table (rc->run, "SRR", h) ;
	      if (! tt)
		tt = ac_tag_table (rc->run, "SRX", h) ;
	      if (! tt)
		tt = ac_tag_table (rc->run, "sublibraries", h) ;
	      if (! tt)
		aceOutf (qc->ao, EMPTY) ;
	      else 
		qcShowMultiTag (qc, tt) ;
	      break ;
	    case 2:  /* Author */
	      if (rc->srr)
		{
		  ccp = rc->run ? ac_tag_printable (rc->run, "Author", 0) : 0 ;
		  if (! ccp) 
		    ccp = rc->srr ? ac_tag_printable (rc->srr, "Magic_Author2", 0) : 0 ;
		  if (! ccp) 
		    ccp = rc->srr ? ac_tag_printable (rc->srr, "Author", EMPTY) : EMPTY ;
		  aceOutf (qc->ao, "%s", ccp) ;
		}
	      else
		{
		  int nAuthor = 0 ;
		  vTXT txt = vtxtHandleCreate (h) ;

		  ccp = rc->run ? ac_tag_printable (rc->run, "Author", 0) : 0 ;
		  if (ccp)
		    {
		      vtxtPrintf (txt, "%s", ccp) ;
		      nAuthor++ ;
		    }
		   ccp = rc->run ? ac_tag_printable (rc->run, "Sequencing_laboratory", 0) : 0 ;
		   if (ccp && (!nAuthor ||  ! strstr (vtxtPtr (txt), ccp)))
		    {
		      vtxtPrintf (txt, "%s%s"
				  , nAuthor++ ? ", " : ""
				  , ccp
				  ) ;
		    }
		   ccp = rc->run ? ac_tag_printable (rc->run, "Nucleic_prep_author", 0) : 0 ;
		   if (ccp && (!nAuthor ||  ! strstr (vtxtPtr (txt), ccp)))
		    {
		      vtxtPrintf (txt, "%s%s"
				  , nAuthor++ ? ", " : ""
				  , ccp
				  ) ;
		    }
		   if (nAuthor)
		     aceOutf (qc->ao, "%s", vtxtPtr (txt)) ;

		   txt = vtxtHandleCreate (h) ;
		   nAuthor = 0 ;
		   ccp = ac_tag_printable (rc->run, "Submission_date", 0) ;
		   if (ccp && (!nAuthor ||  ! strstr (vtxtPtr (txt), ccp)))
		     {
		       vtxtPrintf (txt, "%s%s"
				   , nAuthor++ ? ", " : ""
				   , ccp
				   ) ;
		     }		    
		   ccp = ac_tag_printable (rc->run, "Library_date", 0) ;
		   if (ccp && (!nAuthor ||  ! strstr (vtxtPtr (txt), ccp)))
		     {
		       vtxtPrintf (txt, "%s library prep %s"
				   , nAuthor++ ? ", " : ""
				   , ccp
				   ) ;
		     }		    
		   ccp = ac_tag_printable (rc->run, "Sequencing_date", 0) ;
		   if (ccp && (!nAuthor ||  ! strstr (vtxtPtr (txt), ccp)))
		     {
		       vtxtPrintf (txt, "%s Sequencing %s"
				   , nAuthor++ ? ", " : ""
				   , ccp
				   ) ;
		     }		    
		   if (nAuthor)
		     aceOutf (qc->ao, " [%s]", vtxtPtr (txt)) ;
	
		}	      
	      break ;
	    case 11:  /* Date */
	      {{
		  mytime_t d1 = ac_tag_date (rc->run, "Submission_date", 0 ) ;
		  mytime_t d2 = ac_tag_date (rc->run, "Release_date", 0) ;
		  char buf25[25] ;
		  if (d2 && d2 < d1) d1 = d2 ;
		  if (d1)
		    aceOutf (qc->ao, timeShow (d1, buf25, 25)) ;
		}}
	      break ;
	    case 3:  /* Species */
	      ccp = rc->run ? ac_tag_printable (rc->run, "Species", EMPTY) : EMPTY ;
	      aceOutf (qc->ao, "%s", ccp) ;
	      break ;
	    case 40:  /* SRP-Project */
	      if (rc->srp) aceOutf (qc->ao, "%s", ac_name (rc->srp)) ;
	      break ;
	    case 41:  /* SRP-Project title */
	      ccp = rc->srp ? ac_tag_printable (rc->srp, "Title", EMPTY) : EMPTY ;
	      aceOutf (qc->ao, "%s", ccp) ;
	      break ;
	    case 42:  /* SRP-Project description */
	      ccp = rc->srp ? ac_tag_printable (rc->srp, "Abstract", EMPTY) : EMPTY ;
	      aceOutf (qc->ao, "%s", ccp) ;
	      break ;
	    case 6:  /* Reference */
	      tt = rc->run ? ac_tag_table (rc->run, "Reference", h) : 0 ; 
	      if (tt)
		{
		  int ir ;
		  for (ir = 0 ; ir < tt->rows ; ir++)
		    {
		      AC_OBJ pap = ac_table_obj (tt, ir, 0, h) ;
		      if (! papDict)
			papDict = dictHandleCreate (32, h) ;
		      if (! dictAdd (papDict, ac_name(pap), 0))
			continue ;			
		      aceOutf (qc->ao, "%sPMID: %s", nw ? "; " : EMPTY, ac_name(pap)+(!strncmp( ac_name(pap),"pm",2) ? 2 : 0)) ;
		      ccp = pap ? ac_tag_printable (pap, "Title", 0) : 0 ; 
		      if (ccp) 
			{ aceOutf (qc->ao, " : %s", ccp) ; nw++ ; }
		      ccp = pap ? ac_tag_printable (pap, "Citation", 0) : 0 ;
		      if (ccp) aceOutf (qc->ao, "%s%s", nw ? ", " : ": ", ccp) ;
		      ccp = pap ? ac_tag_printable (pap, "Author", 0) : 0 ;
		      if (ccp) 
			{
			  AC_TABLE tt1 = ac_tag_table (pap, "Author", h) ;
			  aceOutf (qc->ao, ", %s", ccp) ;
			  if (tt1 && tt1->rows > 1) aceOutf (qc->ao, " et al.") ;
			}
		    }
		}
	      tt = rc->srp ? ac_tag_table (rc->srp, "Reference", h) : 0 ; 
	      if (tt)
		{
		  int ir ;
		  for (ir = 0 ; ir < tt->rows ; ir++)
		    {
		      AC_OBJ pap = ac_table_obj (tt, ir, 0, h) ;
		      if (! papDict)
			papDict = dictHandleCreate (32, h) ;
		      if (! dictAdd (papDict, ac_name(pap), 0))
			continue ;			
		      aceOutf (qc->ao, "%sPMID: %s", nw ? "; " : EMPTY, ac_name(pap)+(!strncmp( ac_name(pap),"pm",2) ? 2 : 0)) ;
		      ccp = pap ? ac_tag_printable (pap, "Title", 0) : 0 ; 
		      if (ccp) 
			{ aceOutf (qc->ao, " : %s", ccp) ; nw++ ; }
		      ccp = pap ? ac_tag_printable (pap, "Citation", 0) : 0 ;
		      if (ccp) aceOutf (qc->ao, "%s%s", nw ? ", " : ": ", ccp) ;
		      ccp = pap ? ac_tag_printable (pap, "Author", 0) : 0 ;
		      if (ccp) 
			{
			  AC_TABLE tt1 = ac_tag_table (pap, "Author", h) ;
			  aceOutf (qc->ao, ", %s", ccp) ;
			  if (tt1 && tt1->rows > 1) aceOutf (qc->ao, " et al.") ;
			}
		    }
		}
	      tt = rc->srx ? ac_tag_table (rc->srx, "Reference", h) : 0 ; 
	      if (tt)
		{
		  int ir ;
		  for (ir = 0 ; ir < tt->rows ; ir++)
		    {
		      AC_OBJ pap = ac_table_obj (tt, ir, 0, h) ;
		      if (! papDict)
			papDict = dictHandleCreate (32, h) ;
		      if (! dictAdd (papDict, ac_name(pap), 0))
			continue ;			
		      
		      aceOutf (qc->ao, "%sPMID: %s", nw ? "; " : EMPTY, ac_name(pap)+(!strncmp( ac_name(pap),"pm",2) ? 2 : 0)) ;
		      ccp = pap ? ac_tag_printable (pap, "Title", 0) : 0 ; 
		      if (ccp) 
			{ aceOutf (qc->ao, " : %s", ccp) ; nw++ ; }
		      ccp = pap ? ac_tag_printable (pap, "Citation", 0) : 0 ;
		      if (ccp) aceOutf (qc->ao, "%s%s", nw ? ", " : ": ", ccp) ;
		      ccp = pap ? ac_tag_printable (pap, "Author", 0) : 0 ;
		      if (ccp) 
			{
			  AC_TABLE tt1 = ac_tag_table (pap, "Author", h) ;
			  aceOutf (qc->ao, ", %s", ccp) ;
			  if (tt1 && tt1->rows > 1) aceOutf (qc->ao, " et al.") ;
			}
		    }
		}
	      tt = rc->srr ? ac_tag_table (rc->srr, "Reference", h) : 0 ; 
	      if (tt)
		{
		  int ir ;
		  for (ir = 0 ; ir < tt->rows ; ir++)
		    {
		      AC_OBJ pap = ac_table_obj (tt, ir, 0, h) ;
		      if (! papDict)
			papDict = dictHandleCreate (32, h) ;
		      if (! dictAdd (papDict, ac_name(pap), 0))
			continue ;			
		      
		      aceOutf (qc->ao, "%sPMID: %s", nw ? "; " : EMPTY, ac_name(pap)+(!strncmp( ac_name(pap),"pm",2) ? 2 : 0)) ;
		      ccp = pap ? ac_tag_printable (pap, "Title", 0) : 0 ; 
		      if (ccp) 
			{ aceOutf (qc->ao, " : %s", ccp) ; nw++ ; }
		      ccp = pap ? ac_tag_printable (pap, "Citation", 0) : 0 ;
		      if (ccp) aceOutf (qc->ao, "%s%s", nw ? ", " : ": ", ccp) ;
		      ccp = pap ? ac_tag_printable (pap, "Author", 0) : 0 ;
		      if (ccp) 
			{
			  AC_TABLE tt1 = ac_tag_table (pap, "Author", h) ;
			  aceOutf (qc->ao, ", %s", ccp) ;
			  if (tt1 && tt1->rows > 1) aceOutf (qc->ao, " et al.") ;
			}
		    }
		}
	      if (! nw)
		aceOutf (qc->ao, EMPTY) ;
	      ac_free (papDict) ;
	      break ;

	    case 20:   /* sample */
	      if (rc->sample)
		{
		  ccp = ac_tag_printable (rc->sample, "Title", 0) ;
		  if (! ccp)
		    ccp = ac_name (rc->sample) ;
		}
	      else
		{
		  ccp = ac_tag_printable (rc->run, "Sample", 0) ;
		}
	      if (! ccp) ccp = EMPTY ;
	      aceOutf (qc->ao, "%s", ccp) ;
	      break ;
	    case 7:   /* first line */
	      ccp = ac_tag_printable (rc->run, "First_line", 0) ;
	      if (! ccp) ccp = EMPTY ;
	      aceOutf (qc->ao, "%s", ccp) ;
	      break ;
	    case 21:   /* stage */
	      tt = rc->sample? ac_tag_table (rc->sample, "Stage", h) : 0 ;
	      nw = qcShowTable (qc, tt, 0) ;
	      if (! nw)
		aceOutf (qc->ao, EMPTY) ;
	      break ;
	    case 22:   /* system */
	      tt = rc->sample? ac_tag_table (rc->sample, "Tissue", h) : 0 ;
	      nw = qcShowTable (qc, tt, 0) ;
	      tt = rc->sample? ac_tag_table (rc->sample, "Systm", h) : 0 ;
	      nw += qcShowTable (qc, tt, nw) ;
	      if (! nw)
		aceOutf (qc->ao, EMPTY) ;
              break ;
	    case 23: /* bio */
	      tt = rc->run ? ac_tag_table (rc->run, "CT", h) : 0 ;
	      nw = qcShowTable (qc, tt, nw) ;
	      tt = rc->run ? ac_tag_table (rc->run, "Bio", h) : 0 ;
	      nw += qcShowTable (qc, tt, nw) ;
	      if (! nw)
		aceOutf (qc->ao, EMPTY) ;
	      break ;
	      if (! nw)
		aceOutf (qc->ao, EMPTY) ;
	      break ;
	    case 24: /* sex */
	      tt = rc->sample? ac_tag_table (rc->sample, "Sex", h) : 0 ;
	      nw = qcShowTable (qc, tt, nw) ;
	      if (! nw)
		aceOutf (qc->ao, EMPTY) ;
	      break ;
	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  ac_free (h) ;
  return;
}  /* qcTitles */

/*************************************************************************************/

static void qcProtocol (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  AC_KEYSET aks ;
  const char *ccp ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,

    { "Compute", "Sequencing platform and Machine type", 12, 0, 0} ,
    {  "Compute", "Experiment type", 1, 0, 0} ,
    {  "Compute", "Sequencing_date" , 13, 0, 0} ,
    {  "Compute", "Submission_date" , 14, 0, 0} , /* or received date if earlier */
    {  "Compute", "Library preparation", 2, 0, 0} , /* mv to qcTitle */
    {  "RNA_Integrity_Number" , "RNA Integrity Number", 0, 0, 0} ,
    {  "Details", "Details", 0, 0, 0} ,
    {  "MicroArray", "Corresponding microarray", 0, 0, 0} ,
    { "Compute", "Spots, bases, average read length (in SRA)", 30, 0, 0} ,

    /* {  "Machine", "Machine", 0, 0, 0} , */
   /*   {  "Compute", "Strand", 12, 0, 0} , */
    {  0, 0, 0, 0, 0}
  } ; 

  const char *caption =
    "Experimental protocol"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {
      int nw = 0 ;
      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title) ;
	}
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  aceOutf (qc->ao, "\t") ;
	  switch (ti->col)
	    {
	    case 1:  /* RNA selection */
	      tt = 0 ;
              if (! tt) tt = rc->srr ? ac_tag_table (rc->srr, "Nucleic_Acid_Extraction", h) : 0 ;
              if (! tt) tt = rc->srr ? ac_tag_table (rc->srr, "sraNucleic_Acid_Extraction", h)  : 0 ;
              if (! tt) tt = rc->run ? ac_tag_table (rc->run, "Nucleic_Acid_Extraction", h) : 0 ;
	      if (tt) 
		qcShowMultiTag (qc, tt) ; 
	      else
		aceOutf (qc->ao, EMPTY) ;
	      break ;
	      
	    case 2:  /* Library preparation*/
	      tt = 0 ; /* was before SRX_DB: ac_tag_table (rc->run, "Origin",h) ; */
	      if (tt) 
		qcShowMultiTag (qc, tt) ;
	      else
		aceOutf (qc->ao, EMPTY) ;
	      break ;

	    case 1000:  /* RunId */
	      if (ac_has_tag (rc->run, "RunId"))
		tt = ac_tag_table (rc->run, "RunId", h) ;
	      else
		tt = ac_keyset_table (ac_objquery_keyset (rc->run, ">Union_of ; >RunId", h), 0, 0, 0, h) ;
	      if (tt) 
		{
		  if (tt->rows > 30)
		    aceOutf (qc->ao, "%d > 30 runs, list skipped", tt->rows) ;
		  else
		    qcShowMultiTag (qc, tt) ; 
		}
	      else
		 aceOutf (qc->ao, EMPTY) ;
	      break ;
	    case 12:  /* Machine */
	      aks = ac_objquery_keyset (rc->run, "{Machine} SETOR {>sublibraries; Machine}", h) ;
	      tt = aks ? ac_bql_table (qc->db, "select m1,m2 from r in @, m1 in r->machine, m2 in m1[1]", aks, 0, 0, h) : 0 ;
	      if (tt && tt->rows > 1)
		tt = aks ? ac_bql_table (qc->db, "select m1 from r in @, m1 in r->machine", aks, 0, 0, h) : 0 ;
	      if (tt && tt->rows > 1)
		tt = 0 ;
	      if (! tt)
		aceOutf (qc->ao, "-") ;
	      else 
		{
		  aceOutf (qc->ao, "%s", ac_table_printable (tt, 0, 0, EMPTY)) ; 
		  if (tt->cols > 1)
		    aceOutf (qc->ao, ", %s", ac_table_printable (tt, 0, 1, EMPTY)) ; 
		}
	   
	      if (ac_has_tag (rc->run, "Paired_end") ||
		  ac_keyset_count (ac_objquery_keyset (rc->run, ">sublibraries; Paired_end", h))
		  )
		aceOutf (qc->ao, ", %s", "Paired end") ;
  	      else if (nw)
		aceOutf (qc->ao, ", %s", "Single end") ;

	      break ; 

	    case 13:  /* Sequencing_date */
	      {
		ccp = ac_tag_printable (rc->run, "Sequencing_date", 0) ;
		if (! ccp)
		  {
		    AC_OBJ subrun = ac_objquery_obj (rc->run, "{Sequencing_date} SETOR {>sublibraries; Sequencing_date}", h) ;
		    if (subrun)
		      ccp = ac_tag_printable (subrun, "Sequencing_date", 0) ;
		  }
		if (ccp)
		  aceOutf (qc->ao, "%s", ccp) ;
		else
		  aceOutf (qc->ao, "") ;
		}
	      break ;
	    case 14:  /* Submission_date */
	      {
		ccp = ac_tag_printable (rc->run, "Submission_date", 0) ;
		if (! ccp)
		  {
		    AC_OBJ subrun = ac_objquery_obj (rc->run, "{Submission_date} SETOR {>sublibraries; Submission_date}", h) ;
		    if (subrun)
		      ccp = ac_tag_printable (subrun, "Submission_date", 0) ;
		  }
		if (ccp)
		  aceOutf (qc->ao, "%s", ccp) ;
		else
		  aceOutf (qc->ao, "") ;
		}
	      break ;
	    case 30:  /* Spots */
	      tt = rc->srr ? ac_tag_table (rc->run, "Spots", h) : 0 ;
	      if (! tt)
		tt = rc->srr ? ac_tag_table (rc->srr, "Spots", h) : 0 ;
	      if (! tt)
		tt = rc->srx ? ac_tag_table (rc->ali, "Spots", h) : 0 ;
	      if (tt)
		{
		  ccp = ac_table_printable (tt, 0, 0, 0) ;
		  if (ccp)
		    { aceOutf (qc->ao, "%s spots", ccp) ; nw++ ;}
		  ccp = ac_table_printable (tt, 0, 2, 0) ;
		  if (ccp)
		    { aceOutf (qc->ao, ", %s bases", ccp) ; nw++ ;}
		  ccp = ac_table_printable (tt, 0, 4, 0) ; 
		  if (ccp)
		    { aceOutf (qc->ao, ", average read length %s bases", ccp) ; nw++ ;}
		}
	      else
		{ 
		  int row ;
		  double n1 = 0, n2 = 0, n3 ;
		  AC_KEYSET aks = ac_objquery_keyset (rc->run, ">sublibraries; spots", h) ;
		  tt = ac_bql_table (qc->db, "select r, s, b from r in @, s in r->spots, b in s[2]", aks, 0, 0, h) ;
		  if (tt)
		    {
		      for (row = 0 ; row < tt->rows ; row++)
			{
			  n1 += ac_table_int (tt, row, 1, 0) ;
			  ccp = ac_table_printable (tt, row, 2, 0) ;
			  n3 = 0 ;
			  if (ccp)
			    sscanf (ccp, "%lf", &n3) ;
			  n2 += n3 ;
			}
		      if (n1 > 0)
			aceOutf (qc->ao, "%g spots, %g bases, average read length %.0f bases", n1, n2, n2/n1) ;
		    }		  
		}
	      break ;
	    }
	}
      else if (! strcmp (ti->tag, "Machine"))
	{  
	  aceOutf (qc->ao, "\t") ;
	  tt = ac_tag_table (rc->run, "Machine", h) ;
	  if (tt)
	    qcShowMultiTag (qc, tt) ;
	  else
	    aceOutf (qc->ao, EMPTY) ;
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  ac_free (h) ;
  return;
}  /* qcProtocol */

/*************************************************************************************/

static void qcReadFate1 (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  int ir, i, pass ;
  long int z, zb, zc, zd, zm, z31, zdelta ;
  const char *ccp ;
  AC_TABLE tt = 0 ;

  TT *ti, tts[] = {
     { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
     { "Compute", 
       "% Well mapped reads"
       "\t% Reads mapping uniquely  to a single locus (genomic site) and maximum 1 gene"
       "\t% Reads mapping uniquely to a single locus, but to 2 genes in antisense "
       "\t% Reads mapping uniquely to a single locus"
       "\t% Reads mapping to 2 to 9 sites"
       "\t% Reads between 18 and 31 bases, before mapping"  
       "\t% Reads mapping to microbiome"  
       "\t% reads rejected because they map to 10 sites or more"
       "\t% reads rejected by alignment quality filter"
       "\t% reads rejected before alignment because of low entropy (<= 16 bp)"
       "\t% reads with no insert (adaptor only)"
       "\t% reads unaligned from half aligned fragments"
       "\t% reads of good quality from unaligned fragments (missing genome, microbiome or mapping difficulty)"
       "\t% reads with at least 15 bases 5p overhang"
       "\t% reads with at least 15 bases 3p overhang"

       , 1, 0, 0} ,
    
    {  0, 0, 0, 0, 0}
  }; 

  const char *caption =
    "Reads mapped, unmapped or rejected, accounting for all reads in the file"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	aceOutf (qc->ao, "\t%s", ti->title) ;
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  switch (ti->col)
	    {
	    case 1:
	      for (pass = 1 ; pass < 2 ; pass++)
		{
		  /* raw reads */					
		  zdelta = rc->var[T_Read] ;
		  if (! pass) qcShowPercent (qc, rc, zdelta, pass) ;
		  /*   "\t% reads aligned on any target" */

		  z = z31 = 0 ;
		  tt = ac_tag_table (rc->ali, "Preclipped_length", h) ;
		  for (ir = 0 ; tt && ir < tt->rows ; ir++)	
		    {
		      int j =  ac_table_int (tt, ir, 0, 0) ;
		      z = ac_table_float (tt, ir, 1, 0) ;
		      if (j >= 0 && j <= 31) { z31 += z ; }
		    }

		  tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
		  z = zm = 0 ;
		  for (ir = 0 ; tt && ir < tt->rows ; ir++)
		    {
		      ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		      if (! strcasecmp (ccp, "any"))
			{
			  z = ac_table_float (tt, ir, 3, 0) ;  /* tag number */
			}
		      else if (! strcasecmp (ccp, "b_bacteria") || ! strcasecmp (ccp, "v_virus"))
			{
			  zm += ac_table_float (tt, ir, 3, 0) ;
 			}
		    } 

		  qcShowPercent (qc, rc, z, pass) ;
		  
		  /*
		   * "\t% reads mapping uniquely  to a single locus (genomic site) and maximum 1 gene"
		   * "\t% reads mapping uniquely to a single locus, but to 2 genes in antisense "
		   * "\t% reads mapping uniquely to a single locus"
		   * "\t% reads mapping to 2 to 9 sites"
		   */
		  
		  tt = ac_tag_table (rc->ali, "Unicity", h) ;
		  z = zb = zc = 0 ;
		  for (ir = 0 ; tt && ir < tt->rows ; ir++)
		    {
		      ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		      if (! strcasecmp (ccp, "any"))
			{
			  z += ac_table_float (tt, ir, 1, 0) ;
			  for (i  = 2 ; i <= 10 ; i++)
			    zb += ac_table_float (tt, ir, i, 0) ;
			  zc +=  ac_table_float (tt, ir, 11, 0) ; /* multiplicity -2 */
			}
		      else if (! strcasecmp (ccp, "b_bacteria") || ! strcasecmp (ccp, "v_virus"))
			{
			  z -= ac_table_float (tt, ir, 1, 0) ;
			  for (i  = 2 ; i <= 10 ; i++)
			    zb -= ac_table_float (tt, ir, i, 0) ;
			  zc -=  ac_table_float (tt, ir, 11, 0) ; /* multiplicity -2 */
			}
		    }  
		  zd = z + zc ;
		  qcShowPercent (qc, rc, z, pass) ; /* single locus */
		  qcShowPercent (qc, rc, zc, pass) ;  /* 2 genes in antisense */
		  qcShowPercent (qc, rc, zd, pass) ; /* uniquely to single locus */
		  qcShowPercent (qc, rc, zb, pass) ; /* 2 to 9 sites */
		  qcShowPercent (qc, rc, z31, pass) ; /* microbiome */
		  qcShowPercent (qc, rc, zm, pass) ; /* microbiome */
		  zdelta -= (z + zb + zc + zm) ;
		  /* "\t% reads rejected because they map to 10 sites or more" */
		  tt = ac_tag_table (rc->ali, "At_least_10_sites", h) ;
		  z = zb = zc = 0 ;
		  for (ir = 0 ; tt && ir < tt->rows ; ir++)
		    {
		      if (ir == 0)
			{
			  z = ac_table_float (tt, ir, 2, 0) ;
			}
		    } 
		  qcShowPercent (qc, rc, z, pass) ;
		  zdelta -= z ;

		  /*    "\t% reads rejected by alignment quality filter" */
		  tt = ac_tag_table (rc->ali, "Low_quality_mapping", h) ;
		  z = tt ? ac_table_float (tt, 0, 2, 0) : 0 ;
		  qcShowPercent (qc, rc, z, pass) ;
		  zdelta -= z ;
		  

		  /*  "\t% reads rejected before alignment because of low entropy (<= 16 bp)" */
		  tt = ac_tag_table (rc->ali, "Rejected", h) ;
		  z = tt ? ac_table_float (tt, 0, 2, 0) : 0 ;
		  qcShowPercent (qc, rc, z, pass) ;
		  zdelta -= z ;

		   /*   "\t% reads with no insert (adaptor only)" */
		  tt = ac_tag_table (rc->ali, "No_insert", h) ;
		  z = tt ? ac_table_float (tt, 0, 2, 0) : 0 ;
		  qcShowPercent (qc, rc, z, pass) ;
		  zdelta -= z ;

		   /*   "\t% reads unaligned from half aligned fragments" */
		  tt = ac_tag_table (rc->ali, "Orphans", h) ;
		  z = tt ? ac_table_float (tt, 0, 2, 0) : 0 ;
		  for (ir = 0 ; tt && ir < tt->rows ; ir++)
		    {
		      ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		      if (! strcasecmp (ccp, "Any"))
			{
			  z = ac_table_float (tt, ir, 1, 0) ;
			}
		    }  
		  qcShowPercent (qc, rc, z, pass) ;
		  zdelta -= z ;

		  /* remainder = Unaligned reads (missing genome or microbiome or mapping problem)" */
		  qcShowPercent (qc, rc, zdelta, pass) ;

		  tt = ac_tag_table (rc->ali, "Partial_5p", h) ;
		  z = tt ? ac_table_float (tt, 0, 0, 0) : 0 ;
		  qcShowPercent (qc, rc, z, pass) ;

		  tt = ac_tag_table (rc->ali, "Partial_3p", h) ;
		  z = tt ? ac_table_float (tt, 0, 0, 0) : 0 ;
		  qcShowPercent (qc, rc, z, pass) ;
		}
	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }

   ac_free (h) ;
  return;
}  /* qcReadFate1 */

/*************************************************************************************/

static void qcReadFate2 (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  int ir, i, pass ;
  long int z, zb, zc, zd, zm, z31, zdelta ;
  const char *ccp ;
  AC_TABLE tt = 0 ;

  TT *ti, tts[] = {
     { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
     { "Compute", 
       "Raw reads"
       "\tWell mapped reads"
       "\tReads mapping uniquely  to a single locus (genomic site) and maximum 1 gene"
       "\tReads mapping uniquely to a single locus, but to 2 genes in antisense " 
       "\tReads mapping uniquely to a single locus" 
       "\tReads mapping to 2 to 9 sites"    
       "\tReads between 18 and 31 bases, before mapping"  
       "\tReads mapping to microbiome"  
       "\tReads rejected because they map to 10 sites or more" 
       "\tReads rejected by alignment quality filter"
       "\tReads rejected before alignment because of low entropy (<= 16 bp)"
       "\tReads with no insert (adaptor only)"  
       "\tReads unaligned from half aligned fragments"
       "\tReads of good quality from unaligned fragments (missing genome, microbiome or mapping difficulty)"
       "\tReads partially mapped with at least 15 bases 5p overhang"
       "\tReads partially mapped with at least 15 bases 3p overhang"

       , 1, 0, 0} ,
    
    {  0, 0, 0, 0, 0}
  }; 

  const char *caption =
    "Number of reads mapped, unmapped or rejected, accounting for all reads in the file"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	aceOutf (qc->ao, "\t%s", ti->title) ;
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  switch (ti->col)
	    {
	    case 1:
	      for (pass = 0 ; pass < 1 ; pass++)
		{
		  /* raw reads */					
		  zdelta = rc->var[T_Read] ;
		  if (! pass) qcShowPercent (qc, rc, zdelta, pass) ;
		  /*   "\t% reads aligned on any target" */

		  z = z31 = 0 ;
		  tt = ac_tag_table (rc->ali, "Preclipped_length", h) ;
		  for (ir = 0 ; tt && ir < tt->rows ; ir++)	
		    {
		      int j =  ac_table_int (tt, ir, 0, 0) ;
		      z = ac_table_float (tt, ir, 1, 0) ;
		      if (j >= 0 && j <= 31) { z31 += z ; }
		    }

		  tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
		  z = zm = 0 ;
		  for (ir = 0 ; tt && ir < tt->rows ; ir++)
		    {
		      ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		      if (! strcasecmp (ccp, "any"))
			{
			  z = ac_table_float (tt, ir, 3, 0) ;  /* tag number */
			}
		      else if (! strcasecmp (ccp, "b_bacteria") || ! strcasecmp (ccp, "v_virus"))
			{
			  zm += ac_table_float (tt, ir, 3, 0) ;
 			}
		    } 

		  qcShowPercent (qc, rc, z, pass) ;
		  
		  /*
		   * "\t% reads mapping uniquely  to a single locus (genomic site) and maximum 1 gene"
		   * "\t% reads mapping uniquely to a single locus, but to 2 genes in antisense "
		   * "\t% reads mapping uniquely to a single locus"
		   * "\t% reads mapping to 2 to 9 sites"
		   */
		  
		  tt = ac_tag_table (rc->ali, "Unicity", h) ;
		  z = zb = zc = 0 ;
		  for (ir = 0 ; tt && ir < tt->rows ; ir++)
		    {
		      ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		      if (! strcasecmp (ccp, "any"))
			{
			  z += ac_table_float (tt, ir, 1, 0) ;
			  for (i  = 2 ; i <= 10 ; i++)
			    zb += ac_table_float (tt, ir, i, 0) ;
			  zc +=  ac_table_float (tt, ir, 11, 0) ; /* multiplicity -2 */
			}
		      else if (! strcasecmp (ccp, "b_bacteria") || ! strcasecmp (ccp, "v_virus"))
			{
			  z -= ac_table_float (tt, ir, 1, 0) ;
			  for (i  = 2 ; i <= 10 ; i++)
			    zb -= ac_table_float (tt, ir, i, 0) ;
			  zc -=  ac_table_float (tt, ir, 11, 0) ; /* multiplicity -2 */
			}
		    }  
		  zd = z + zc ;
		  qcShowPercent (qc, rc, z, pass) ; /* single locus */
		  qcShowPercent (qc, rc, zc, pass) ;  /* 2 genes in antisense */
		  qcShowPercent (qc, rc, zd, pass) ; /* uniquely to single locus */
		  qcShowPercent (qc, rc, zb, pass) ; /* 2 to 9 sites */
		  qcShowPercent (qc, rc, z31, pass) ; /* microbiome */
		  qcShowPercent (qc, rc, zm, pass) ; /* microbiome */
		  zdelta -= (z + zb + zc + zm) ;
		  /* "\t% reads rejected because they map to 10 sites or more" */
		  tt = ac_tag_table (rc->ali, "At_least_10_sites", h) ;
		  z = zb = zc = 0 ;
		  for (ir = 0 ; tt && ir < tt->rows ; ir++)
		    {
		      if (ir == 0)
			{
			  z = ac_table_float (tt, ir, 2, 0) ;
			}
		    } 
		  qcShowPercent (qc, rc, z, pass) ;
		  zdelta -= z ;

		  /*    "\t% reads rejected by alignment quality filter" */
		  tt = ac_tag_table (rc->ali, "Low_quality_mapping", h) ;
		  z = tt ? ac_table_float (tt, 0, 2, 0) : 0 ;
		  qcShowPercent (qc, rc, z, pass) ;
		  zdelta -= z ;
		  

		  /*  "\t% reads rejected before alignment because of low entropy (<= 16 bp)" */
		  tt = ac_tag_table (rc->ali, "Rejected", h) ;
		  z = tt ? ac_table_float (tt, 0, 2, 0) : 0 ;
		  qcShowPercent (qc, rc, z, pass) ;
		  zdelta -= z ;

		   /*   "\t% reads with no insert (adaptor only)" */
		  tt = ac_tag_table (rc->ali, "No_insert", h) ;
		  z = tt ? ac_table_float (tt, 0, 2, 0) : 0 ;
		  qcShowPercent (qc, rc, z, pass) ;
		  zdelta -= z ;

		   /*   "\t% reads unaligned from half aligned fragments" */
		  tt = ac_tag_table (rc->ali, "Orphans", h) ;
		  z = tt ? ac_table_float (tt, 0, 2, 0) : 0 ;
		  for (ir = 0 ; tt && ir < tt->rows ; ir++)
		    {
		      ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		      if (! strcasecmp (ccp, "Any"))
			{
			  z = ac_table_float (tt, ir, 1, 0) ;
			}
		    }  
		  qcShowPercent (qc, rc, z, pass) ;
		  zdelta -= z ;

		  /* remainder = Unaligned reads (missing genome or microbiome or mapping problem)" */
		  qcShowPercent (qc, rc, zdelta, pass) ;

		  tt = ac_tag_table (rc->ali, "Partial_5p", h) ;
		  z = tt ? ac_table_float (tt, 0, 0, 0) : 0 ;
		  qcShowPercent (qc, rc, z, pass) ;

		  tt = ac_tag_table (rc->ali, "Partial_3p", h) ;
		  z = tt ? ac_table_float (tt, 0, 0, 0) : 0 ;
		  qcShowPercent (qc, rc, z, pass) ;
		}
	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }

   ac_free (h) ;
  return;
}  /* qcReadFate2 */

/*************************************************************************************/
/* telomere */
static void qcBloom (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  int ir ;
  AC_TABLE tt ;
  
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Adaptor", "Exit adaptor of read 1", 1, 0, 0} ,
    { "Adaptor", "Exit adaptor of read 2", 2, 0, 0} ,
    { "Compute", "PolyA index : fragments per million with at least 3 A13 motifs shifted by at least 6 bases", 1, 0, 0} ,
    { "Compute", "Telomere index : fragments per million with at least 3 TTAGGG/CCCTAA motifs", 2, 0, 0} ,
    { "Compute", "Telomere index allowing G/C substitution", 3, 0, 0} ,
    { "Compute", "Noise in telomere index : FPM with at least 3 AATCCC/GGGATT motifs", 4,0,0} ,
    {  0, 0, 0, 0, 0}
  }; 
  const char *caption =
    "Motif search before alignments"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;
  
  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	aceOutf (qc->ao, "\t%s", ti->title) ;
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Adaptor"))
	{
	  const char *ccp ;
	  const char *tag ;
	  switch (ti->col)
	    {
	    case 1: 
	      tag = "Adaptor1" ;
	      break ;
	    case 2: 
	      tag = "Adaptor2" ;
	      break ;
	    } 
	  ccp = ac_tag_printable (rc->run, tag, 0) ;
	  if (! ccp)
	    ccp = ac_tag_printable (rc->ali, tag, 0) ;
	  if (ccp)
	    aceOutf (qc->ao, "\t%s", ccp) ;
	  else
	    aceOutf (qc->ao, "\t%s", EMPTY) ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  float z = -1 ;
	  const char *tag ;
	  /*    tags available in MetaDB 
		"any.A13-per_1_kb_fragment-cumul"
		"any.Tel_C_6-per_1_kb_fragment-cumul" 
		"any.Telomeric_6-per_1_kb_fragment-cumul"
		"any.imagC_6-per_1_kb_fragmen-cumul"
		"any.imagT_6-per_1_kb_fragment-cumul"
		"unaligned.A13-per_1_kb_fragmen-cumul"
		"unaligned.Tel_C_6-per_1_kb_fragmen-cumul"
		"unaligned.Telomeric_6-per_1_kb_fragmen-cumul"
		"unaligned.imagC_6-per_1_kb_fragmen-cumul"
		"unaligned.imagT_6-per_1_kb_fragmen-cumul"
	  */

	  switch (ti->col)
	    {
	    case 1: /* PolyA index */
	      tag = "any.A13-per_1_kb_fragment-cumul" ;
	      break ;
	    case 2: /* Telomere index */
	      tag = "any.Telomeric_6-per_1_kb_fragment-cumul" ;
	      break ;
	    case 3: /* Telomere index with G/C substitution allowed */
	      tag = "any.Tel_C_6-per_1_kb_fragment-cumul" ;
	      break ;
	    case 4: /* Control using the complementary motif */
	      tag = "any.imagT_6-per_1_kb_fragment-cumul" ;
	      break ;
	    } 
	  tt = ac_tag_table (rc->ali, "Bloom", h) ;  
	  for (ir = 0 ; tt && ir < tt->rows ; ir++)
	    {
	      const char *ccp = ac_table_printable (tt, ir, 0, 0) ;
	      if (ccp && ! strcasecmp (tag, ccp))
		{
		  float z2, z5 ;
		  z2 = ac_table_float (tt, ir, 2, -1) ;
		  z5 = ac_table_float (tt, ir, 5, -1) ;
		  if (z2 >= 0)
		    z = 1000000.0 * z5 / z2 ;
		}
	    }
	  if (z >= 0)
	    aceOutf (qc->ao, "\t%.2f", z) ;
	  else
	    aceOutf (qc->ao, "\t%s", EMPTY) ;
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  ac_free (h) ;
  return;
}  /* qcBloom */

/*************************************************************************************/
/* Transpliced leader SL1 SL12 */
static void qcSLsDo (QC *qc, RC *rc, BOOL isSL)
{
  AC_HANDLE h = ac_new_handle () ;
  TT *ti, *tti ;
  TT ttA[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "pA", "polyA site", 1, 0, 0} ,
    { "pA", "polyA support", 2, 0, 0} ,
    { "pA", "polyA support per site", 3, 0, 0}  ,
    {  0, 0, 0, 0, 0}
  } ;
  TT ttS[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "SL1", "SL1 site", 1, 0, 0 } ,
    { "SL2", "SL2 site", 1, 0, 0 } ,
    { "SL3", "SL3 site", 1, 0, 0 } ,
    { "SL4", "SL4 site", 1, 0, 0 } ,
    { "SL5", "SL5 site", 1, 0, 0 } ,
    { "SL6", "SL6 site", 1, 0, 0 } ,
    { "SL7", "SL7 site", 1, 0, 0 } ,
    { "SL8", "SL8 site", 1, 0, 0 } ,
    { "SL9", "SL9 site", 1, 0, 0 } ,
    { "SL10", "SL10 site", 1, 0, 0 } ,
    { "SL11", "SL11 site", 1, 0, 0 } ,
    { "SL12", "SL12 site", 1, 0, 0 } ,
    { "SL1", "SL1 supports", 2, 0, 0 } ,
    { "SL2", "SL2 supports", 2, 0, 0 } ,
    { "SL3", "SL3 supports", 2, 0, 0 } ,
    { "SL4", "SL4 supports", 2, 0, 0 } ,
    { "SL5", "SL5 supports", 2, 0, 0 } ,
    { "SL6", "SL6 supports", 2, 0, 0 } ,
    { "SL7", "SL7 supports", 2, 0, 0 } ,
    { "SL8", "SL8 supports", 2, 0, 0 } ,
    { "SL9", "SL9 supports", 2, 0, 0 } ,
    { "SL10", "SL10 supports", 2, 0, 0 } ,
    { "SL11", "SL11 supports", 2, 0, 0 } ,
    { "SL12", "SL12 supports", 2, 0, 0 } ,
    {  0, 0, 0, 0, 0} 
  } ; 

  const char *caption ;
  if (! isSL)
    {
      caption =  "PolyA sites" ;
      tti = ttA ;
    }
  else 
    {
      caption = "Transpliced leaders sites" ;
      tti = ttS ;
    }
  
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tti, caption) ;
  
  AC_TABLE tbl = 0 ;
  for (ti = tti ; ti->tag ; ti++)
    {
      if (rc == 0)
	aceOutf (qc->ao, "\t%s", ti->title) ;
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
       else  if (*ti->tag)
	{
	  BOOL ok = FALSE ;
	  int  ir ;
	  if (! tbl)
	    tbl = ac_tag_table (rc->ali, "SLs", h) ; 

	  if (tbl)
	    for (ir = 0 ; ir <= tbl->rows ; ir++)
	      {
		if (! strcmp (ac_table_printable (tbl, ir, 0, "toto"), ti->tag))
		  {
		    int np = ac_table_int (tbl, ir, 1, 0) ;
		    int ns = ac_table_int (tbl, ir, 3, 0) ;
		    ok = TRUE ;
		    switch (ti->col)
		      {
		      case 1:
			aceOutf (qc->ao, "\t%d", np) ;
			break ;
		      case 2:
			aceOutf (qc->ao, "\t%d", ns) ;
			break ;
		      case 3:
			if (np == 0) 
			  np = 1 ;
			aceOutf (qc->ao, "\t%.2f", ns/(float)np) ;
			break ;
		      }
		  }
	      }
	  if (! ok)
	    aceOutf (qc->ao, "\t") ;
	}
    }

   ac_free (h) ;
  return;
}  /* qcSLsDo */

static void qcSLs (QC *qc, RC *rc)
{
  qcSLsDo (qc, rc, FALSE) ; /* polyA */
  if (qc->hasSL) 
    qcSLsDo (qc, rc, TRUE)  ;  /* SLs */
} /* qcSLs */

/*************************************************************************************/
/* template for a future chapter */
static void qcOther (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    {  0, 0, 0, 0, 0}
  }; 

  const char *caption =
    "Statistics before alignment"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;
  
  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	aceOutf (qc->ao, "\t%s", ti->title) ;
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else  if (! strcmp (ti->tag, "Compute"))
	{
	  switch (ti->col)
	    {
	    default:
	      break ;
	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }

   ac_free (h) ;
  return;
}  /* qcOther */

/*************************************************************************************/
/* captured genes */
static void qcCapture (QC *qc, RC *rc)
{
  if (! qc->hasCapture) return ;
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tt1 = 0, ttAli = 0 ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute11", "Megabases aligned to captured ILMR1 transcripts\tMegabases aligned in ILMR1 gene boxes",1,0,0},
    { "Compute12", "Megabases aligned to captured ILMR2 transcripts\tMegabases aligned in ILMR2 gene boxes",1,0,0},
    { "Compute13", "Megabases aligned to captured ILMR3 transcripts\tMegabases aligned in ILMR3 gene boxes",1,0,0},
    { "Compute21", "Megabases aligned to captured AGLR1 transcripts\tMegabases aligned in AGLR1 gene boxes",1,0,0},
    { "Compute22", "Megabases aligned to captured AGLR2 transcripts\tMegabases aligned in AGLR2 gene boxes",1,0,0},
    { "Compute31", "Megabases aligned to captured ROCR1 transcripts\tMegabases aligned in ROCR1 gene boxes",1,0,0},
    { "Compute32", "Megabases aligned to captured ROCR2 transcripts\tMegabases aligned in ROCR2 gene boxes",1,0,0},
    { "Compute33", "Megabases aligned to captured ROCR3 transcripts\tMegabases aligned in ROCR3 gene boxes",1,0,0},
    { "Compute34", "Megabases aligned to captured A2I3R2 transcripts\tMegabases aligned in A2I3R2 gene boxes",1,0,0},
    { "Compute11p", "% bases captured by ILMR1 transcripts\t% bases captured by ILMR1 gene boxes",1,0,0},
    { "Compute12p", "% bases captured by ILMR2 transcripts\t% bases captured by ILMR2 gene boxes",1,0,0},
    { "Compute13p", "% bases captured by ILMR3 transcripts\t% bases captured by ILMR3 gene boxes",1,0,0},
    { "Compute21p", "% bases captured by AGLR1 transcripts\t% bases captured by AGLR1 gene boxes",1,0,0},
    { "Compute22p", "% bases captured by AGLR2 transcripts\t% bases captured by AGLR2 gene boxes",1,0,0},
    { "Compute31p", "% bases captured by ROCR1 transcripts\t% bases captured by ROCR1 gene boxes",1,0,0},
    { "Compute32p", "% bases captured by ROCR2 transcripts\t% bases captured by ROCR2 gene boxes",1,0,0},
    { "Compute33p", "% bases captured by ROCR3 transcripts\t% bases captured by ROCR3 gene boxes",1,0,0},
    { "Compute34p", "% bases captured by A2I3R2 transcripts\t% bases captured by A2I3R2 gene boxes",1,0,0},

    {  0, 0, 0, 0, 0}
  }; 

  const char *caption =
    "Capture"
    ;
  if (! getenv ("CAPTURES"))
    return ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;
  
  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	aceOutf (qc->ao, "\t%s", ti->title) ;
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else  if (! strncmp (ti->tag, "Compute", 7))
	{  
	  int ir ;
	  float anyAli = 0 ;
	  float raw = rc->var[T_kb] ;
	  const char *capture = 0 ;
	  char capAli[256] ;
	  char capGbx[256] ;
	  float zba, zbg ;
	  BOOL isP = FALSE ;


	  if (! anyAli) 
	    {
	      ttAli = ac_tag_table (rc->ali, "h_Ali" , h) ;
	      if (ttAli)
		for (ir = 0 ; ir < ttAli->rows ; ir++)
		  {
		    const char *ccp = ac_table_printable (ttAli, ir, 0, 0) ;
		    if (! ccp)
		      continue ;
		    if (! strcasecmp (ccp, "0_SpikeIn") ||
			! strcasecmp (ccp, "1_DNASpikeIn")
			)
		      anyAli -= ac_table_float (ttAli, ir, 5, 0) ;
		    else if (! strcasecmp (ccp, "any"))
		      anyAli += ac_table_float (ttAli, ir, 5, 0) ;
		  }
	      if (anyAli <= 0)
		anyAli = raw ? raw : 1 ;
	    }

	  if (! strcmp (ti->tag, "Compute11"))
	    capture = "I1" ;
	  else if (! strcmp (ti->tag, "Compute12"))
	    capture = "I2" ;
	  else if (! strcmp (ti->tag, "Compute13"))
	    capture = "I3" ;
	  else if (! strcmp (ti->tag, "Compute21"))
	    capture = "A1" ;
	  else if (! strcmp (ti->tag, "Compute22"))
	    capture = "A2" ;
	  else if (! strcmp (ti->tag, "Compute23"))
	    capture = "A3" ;
	  else if (! strcmp (ti->tag, "Compute31"))
	    capture = "R1" ;
	  else if (! strcmp (ti->tag, "Compute32"))
	    capture = "R2" ;
	  else if (! strcmp (ti->tag, "Compute33"))
	    capture = "R3" ;
	  else if (! strcmp (ti->tag, "Compute34"))
	    capture = "A2I3R2" ;

	  else if (! strcmp (ti->tag, "Compute11p"))
	    { capture = "I1" ; isP = TRUE ; } 
	  else if (! strcmp (ti->tag, "Compute12p"))
	    { capture = "I2" ; isP = TRUE ; } 
	  else if (! strcmp (ti->tag, "Compute13p"))
	    { capture = "I3" ; isP = TRUE ; } 
	  else if (! strcmp (ti->tag, "Compute21p"))
	    { capture = "A1" ; isP = TRUE ; } 
	  else if (! strcmp (ti->tag, "Compute22p"))
	    { capture = "A2" ; isP = TRUE ; } 
	  else if (! strcmp (ti->tag, "Compute23p"))
	    { capture = "A3" ; isP = TRUE ; } 
	  else if (! strcmp (ti->tag, "Compute31p"))
	    { capture = "R1" ; isP = TRUE ; } 
	  else if (! strcmp (ti->tag, "Compute32p"))
	    { capture = "R2" ; isP = TRUE ; } 
	  else if (! strcmp (ti->tag, "Compute33p"))
	    { capture = "R3" ; isP = TRUE ; } 
	  else if (! strcmp (ti->tag, "Compute34p"))
	    { capture = "A2I3R2" ; isP = TRUE ; } 
	  
	  if (! capture || strlen (capture) >  210)
	    messcrash ("capture not well defined in qcCapture: capture=#%s#  ti->tag=#%s#", capture ? capture : "NULL", ti->tag) ;
	  sprintf (capAli, "%s_ali", capture) ;
	  sprintf (capGbx, "%s.genebox.ns", capture) ;
	  
	  zba = zbg = 0 ;
	  if (! tt1) tt1 = ac_tag_table (rc->ali, "Capture" , h) ;
	  if (tt1)
	    for (ir = 0 ; ir < tt1->rows ; ir++)
	      {
		if (! strcmp (ac_table_printable (tt1, ir, 0, "x"), capAli))
		  zba = ac_table_float (tt1, ir, 1, 0) ;
		if (! strcmp (ac_table_printable (tt1, ir, 0, "x"), capGbx))
		  zbg = ac_table_float (tt1, ir, 1, 0) ;
	      }
	  if (! isP)
	    {
	      aceOutf (qc->ao, "\t%.3f", zba/1000) ; /* Mb aligned on any target */
	      aceOutf (qc->ao, "\t%.3f", zbg/1000) ; /* Mb aligned on any target */
	    }
	  else 
	    {
	      aceOut (qc->ao, "\t") ; aceOutPercent (qc->ao, 100.0*zba/anyAli) ;
	      aceOut (qc->ao, "\t") ; aceOutPercent (qc->ao, 100.0*zbg/anyAli) ;
	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }

   ac_free (h) ;
  return;
}  /* qcCapture */

/*************************************************************************************/
/* template for a future chapter */
static void qcOneExternalFile (QC *qc, RC *rc, const char *caption, Array titles, Array aaa)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii ;

  if (rc == (void *) 1)
    {
       aceOutf (qc->ao, "\t\t%s", caption) ;
       for (ii = 0 ; ii < arrayMax (titles) ; ii++)
	 aceOut (qc->ao, "\t") ;
       return ;
    }
  
  for (ii = 0 ; ii < arrayMax (titles) ; ii++)
    {
      if (rc == 0)
	{
	  char *title = array (titles, ii, char*) ;
	  if (title)
	    aceOutf (qc->ao, "\t%s", title) ;
	  else
	    aceOut (qc->ao, "\t") ;
	}
      else if (rc && aaa) 
	{
	  int irc = rc - arrp (qc->runs, 0, RC) ;
	  
	  if (irc >= 0 && irc < arrayMax (aaa))
	    {
	      Array aa = arr (aaa, irc, Array) ;
	      char *cp = array (aa, ii, char*) ;
	      if (cp)
		aceOutf (qc->ao, "\t%s", cp) ;
	      else
		aceOut (qc->ao, "\t") ;
	    }
	  else
	    aceOut (qc->ao, "\t") ;
	}
      else
	 aceOut (qc->ao, "\t") ;
    }
  
   ac_free (h) ;
  return;
}  /* qcOneExternalFile */

/*************************************************************************************/

static void qcExternalFiles (QC *qc, RC *rc)
{
  int iiMax  = 0, ii ;
  Array aaa, aa ;
  Array eTitles ;
  
  if (! qc->externalFiles)
    return ;
  if (! qc->eAaa)
    {
      char cutter,  *cp, *cq, *cr, *buf = strnew (qc->externalFiles, qc->h) ; 
      aaa = qc->eAaa = arrayHandleCreate (12, Array, qc->h) ;
      eTitles = qc->eTitles = arrayHandleCreate (12, Array, qc->h) ;
      qc->eCaptions = arrayHandleCreate (12, char *, qc->h) ;
      ii = 0 ;
      cq = buf ;
      while (*cq)
	{
	  AC_HANDLE h = ac_new_handle () ;
	  ACEIN ai = 0 ;
	  int jMax ;
	  cp = cq ; 
	  cq = strchr (cp, ',') ;
	  if (cq) 
	    *cq++ = 0 ;
	  cr = strchr (cp, ':') ;
	  if (cr && cr > cp && cr[1])
	    {
	      *cr++ = 0 ;
	      array (qc->eCaptions, ii, char *) = strnew (cp, qc->h) ;
	      ai = aceInCreate (cr, FALSE, h) ;
	      if (ai) 
		{ 
		  int line = 0, jj ;
		  aceInSpecial (ai, "\n") ;
		  while (aceInCard (ai))
		    {
		      cp = aceInWordCut (ai, "\t", &cutter) ;
		      if (! cp)
			continue ;
		      if (*cp == '#' && cp[1] == '#')
			continue ;
		      if (*cp == '#' && line)
			continue ;
		      if (*cp == '#' && ! line++)
			{
			  jj = 0 ;
			  Array titles = array (eTitles, ii, Array) ;
			  
			  cp++ ; while (*cp == ' ') cp++ ;
			  cq = cp ;
			  while (cq)
			    {
			      cp = cq ;
			      array (titles, jj++, char *) = strnew (cp, qc->h) ;
			      aceInStep (ai, '\t') ;
			      cq = cp = aceInWordCut (ai, "\t", &cutter) ;
			    }
			  jMax = jj ;
			}  
		      if (! cp || ! *cp)
			continue ;
		      {
			int irc ;
			int nr = arrayMax (qc->runs) ;
			for (irc = 0, rc = arrp (qc->runs, 0, RC) ; irc < nr ; irc++, rc++)
			  {
			    if (! strcasecmp (cp, ac_name (rc->run)))
			      {
				aa = array (aaa, irc, Array) =
				  arrayHandleCreate (jMax, char *, qc->h)  ;
				for (jj = 0 ; jj < jMax ; jj++)
				  {
				    cp = aceInWordCut (ai, "\t", &cutter) ;
				    if (cp)
				      array (aa, jj, char *) = strnew (cp, qc->h) ;
				  }
			      }
			  }
		      }
		    }
		}
	      else
		messcrash ("cannot open exteral file %s\n", cr) ;
	    }
	  ac_free (h) ;
	  ii++ ;
	}
    }
  iiMax = arrayMax (qc->eAaa) ;
  for (ii = 0 ; ii < iiMax ; ii++)
    {
      char *caption = array (qc->eCaptions, ii, char*) ;
      Array titles = array (qc->eTitles, ii, Array) ;
      Array aaa = array (qc->eAaa, ii, Array) ;
      if (titles && aaa && arrayMax (aaa))
	qcOneExternalFile (qc, rc, caption, titles, aaa)  ;
    }
} /* qcExternalFiles */

/*************************************************************************************/
/*************************************************************************************/

static void qcCPU (QC *qc, RC *rc)
{
  AC_HANDLE h =  0 ;
  AC_TABLE  tt ;

  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute", "Million reads aligned on all targets per CPU hour", 1, 0, 0 } ,
    { "Compute", "Maximum RAM used (GB)",2, 0, 0} ,
    { "Number_of_lanes", "Number of blocks", 0, 0, 0} ,
    {  0, 0, 0, 0, 0}
  }; 
  
    const char *caption =
    "Aligner performance"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {	
      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title) ;
	}
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
       else if (! strcmp (ti->tag, "Compute"))
	{   
	  int ir, j, k = 0 ;
	  float z = 0 ;
	  switch (ti->col)
	    {
	    case 1: /* million reads per CPU hour  */
	      tt = ac_tag_table (rc->ali, "Elapsed", h) ;
	      if (! tt) tt = ac_tag_table (rc->ali, "CPU", h) ;
 	      for (ir = 0 ; tt && ir < tt->rows ; ir++)	
		k += ac_table_int (tt, ir, 1, 0) ;
	      
	      tt = ac_tag_table (rc->ali, "Accepted", h) ;
	      z = tt ? ac_table_float (tt, 0, 2, 0) : 0 ;

	      z = k ? 3600.0 * z / k : 0 ;
	      aceOutf (qc->ao, "\t%.2f", z/1000000.0) ;
	      break ;
	      
	    case 2:  /* max RAM */ 
	      tt = ac_tag_table (rc->ali, "Max_memory", h) ;
 	      for (ir = 0 ; tt && ir < tt->rows ; ir++)	
		{
		  j = ac_table_int (tt, ir, 1, 0) ;
		  if (k < j) k = j ;
		}
	      aceOutf (qc->ao, "\t%.1f", k/1024.0) ;
	      break ;
	      
	    }
	  ac_free (tt) ;
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  
  ac_free (h) ;
  return;
} /* qcCPU */

/*************************************************************************************/

static void qcMicroRNA (QC *qc, RC *rc)
{
  AC_HANDLE h =  0 ;
  AC_TABLE  tt ;

  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute", "% reads with clipped length 18-31 nt\t18\t19\t20\t21\t22\t23\t24\t25\t26\t27\t28\t29\t30\t31", 1, 0, 0} ,
    { "Compute", "Clipped multiplicity 10\t100\t1000\t10k\t100k\t1M", 2, 0, 0} ,
    { "Compute", "High small",3, 0, 0} ,
    { "Compute", "Mapped on 70k-small_ref\t% Mapped on 70k-small_ref\t% strand plus on 70k small ref", 4, 0, 0} ,
    { "Compute", "Distinct tags clipped  to [18,35] at frequency >= 10^-4 and seen at least 10 times\tTags clipped  to [18,35] at frequency >= 10^-4 and seen at least 10 times\tDistinct tags clipped  to [18,35] at frequency >= 10^-5 and seen at least 10 times\tTags clipped  to [18,35] at frequency >= 10^-5 and seen at least 10 times\tDistinct tags clipped  to [18,35] seen at least 10 times\tTags clipped  to [18,35] seen at least 10 times", 5,0,0} ,
    { "Compute", "Recognized miR [18,35] at frequency >= 10^-4\tSupport\t\t\t\t", 6,0,0} ,
    {  0, 0, 0, 0, 0}
  }; 
  
    const char *caption =
    "Small RNA"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {	
      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title) ;
	}
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
       else if (! strcmp (ti->tag, "Compute"))
	{   
	  BOOL ok ;
	  int ir, i, j ;
	  long int zz[32] ;
	  float z = 0, zzz = 0 ;
	  float up = 0, um = 0 ;

	  memset (zz, 0, sizeof(zz)) ;
	  switch (ti->col)
	    {
	    case 1:  /* Most seen length */
	      tt = ac_tag_table (rc->ali, "Preclipped_length", h) ;
	      for (ir = 0 ; tt && ir < tt->rows ; ir++)	
		{
		  j =  ac_table_int (tt, ir, 0, 0) ;
		  z = ac_table_float (tt, ir, 1, 0) ;
		  if (j >= 0 && j <= 31) { zz[j] = z ; zzz += z ; }
		}
	      qcShowPercent (qc, rc, zzz, TRUE) ;
	      for (j = 18 ; j <= 31 ; j++)
		{ 
		  if (tt)
		    qcShowPercent (qc, rc, zz[j], TRUE) ;
		  else
		    aceOut (qc->ao, "\t") ;
		}
	      break ;
	    case 2:  /* Clipped multiplicity */
	      tt = ac_tag_table (rc->ali, "Clipped_multiplicity", h) ;

	      for (ir = 0 ; tt && ir < tt->cols ; ir+=2)	
		{
		  z = ac_table_int (tt, 0, ir, 0) ;
		  if (ir <= 60) zz[ir/2] = z ;
		}
	      for (j = 0 ; j < 6 ; j++)
		{
		  if (tt)
		    qcShowPercent (qc, rc, zz[j] - zz[j+1], TRUE) ;
		  else
		    aceOut (qc->ao, "\t") ;
		}
	      break ;
	    case 3:  /* high small */
	      tt = ac_tag_table (rc->ali, "High_short", h) ;
 	      aceOut (qc->ao, "\t") ;
	      z = rc->var[T_Read] ;
	      if (z > 0) z = 100/(double)z ;
	         for (ir = 0 ; tt && ir < tt->rows && ir <= 5 ; ir++)	
		{
		  const char *oligo ;
		  j = ac_table_int (tt, ir, 1, 0) ;
		  oligo = ac_table_printable (tt, ir, 0, "") ;
		  if (oligo && j > 0)
		    aceOutf (qc->ao, "%s%s(%dbp;%d;%.2f%%)"
			     , ir > 0 ? ", " : ""
			     , oligo, strlen (oligo) 
			     , j, z*j
			     ) ;
		}
	      break ;
	    case 4:
	      ok = FALSE ;
	      tt = ac_tag_table (rc->ali, "stranding", h) ;
	      for (ir = 0 ; tt && ir < tt->rows; ir++)
		if (! strncmp (ac_table_printable (tt, ir, 0, ""), "70k_small", 9))
		  {		    
		    ok = TRUE ;
		    up += ac_table_float (tt, ir, 2, 0) ;  /* Loop because in groups we may have .f and .f2 cases */
		    um += ac_table_float (tt, ir, 4, 0) ;
		    aceOutf (qc->ao, "\t%.0f\t%.2f\t", up+um,rc->var[T_Read] ? 100 * (up + um) / rc->var[T_Read] : 0) ;
		    if (up + um > 0) aceOutPercent (qc->ao, 100.00 * up/(.0001 + up + um) ) ;
		  }
	      if (! ok)
		aceOutf (qc->ao, "\t\t\t") ;
	      break ;
	    case 5:  /* N_18_35\tKnown miR */
	    case 6:  /* known miR support */
	      tt = ac_tag_table (rc->ali, ti->col == 5 ? "N_18_35" : "Known_miR", h) ;
	      for (i = -4 ; i >= -10 ; i--)
		{
		  int jr = 0 ;
		  int jMax = tt ? tt->rows : 0 ;
		  int t = 0, s = 0 ;
		  for (jr = 0 ; jr < jMax ; jr++)
		    if (i == ac_table_int (tt,jr,0,0))
		      {
			s = ac_table_int (tt,jr,1,0) ;
			t = ac_table_int (tt,jr,2,0) ;
		      }
		  aceOutf (qc->ao, "\t%d\t%d", s, t) ;
		  if (ti->col == 6) break ;
		  if (i == -5) i = -9 ; /* export 10^-4, 10^-5, 10^-10==0 */
		}
	      ac_free (tt) ;
	      tt = ac_tag_table (rc->ali, "Known_miR", h) ;
	      aceOutf (qc->ao, "\t%d", ac_table_int (tt,0,0,0)) ;
	      aceOutf (qc->ao, "\t%d", ac_table_int (tt,0,2,0)) ;
	      break ;
	    }
	  ac_free (tt) ;
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  
  ac_free (h) ;
  return;
} /* qcMicroRNA */

/*************************************************************************************/

static void qcHighVirusBacteria (QC *qc, RC *rc)
{
 AC_HANDLE h = ac_new_handle () ;
 DICT *myDict = 0 ;

  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "HighGenes", "Top endogenous virus and transposable or repeated elements", 1, 0, 0 },
    { "HighGenes", "Top bacteria, microbes, symbionts", 2, 0, 0 },
    { "HighGenes", "Top viruses", 3, 0, 0 },
    {  0, 0, 0, 0, 0}
  }; 
  
  const char *caption =
    "Transposons, bacteria and viruses"
    ;
  const char *targets[3] = {"transposon", "microbe" , "virus"} ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {	
      AC_TABLE tt = 0 ;

      if (rc == 0)
	aceOutf (qc->ao, "\t%s", ti->title) ;
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (ti->col == 0)
	aceOutf (qc->ao, "\t") ;
      else if (ti->col > 0 && ti->col < 4)
	{  
	  int ir, k = 0 ;
	 
	  const char *target = targets[ti->col - 1] ;	  
	  aceOut (qc->ao, "\t") ;
	  tt = ac_tag_table (rc->ali, "High_genes", h) ;
	  for (ir = 0 ; tt && ir < tt->rows ; ir++)
	    {
	      if (! strcasecmp (ac_table_printable (tt, ir, 0, "x"), target))
		{
		  const char *ccp = ac_table_printable (tt, ir, 1, 0) ;
				  
		  if (ccp)
		    {
		      char *cp, *cr ;
		      const char *title = ac_table_printable (tt, ir, 3, 0) ;
		      float x =  ac_table_float (tt, ir, 2, 0) ;
		      if (title) 
			ccp = title ;
		      cp = strnew (ccp, h) ;
		      cr = strstr (cp, ", complete genome") ;
		      if (cr) *cr = 0 ;
		      cr = strstr (cp, "complete genome") ;
		      if (cr) *cr = 0 ;
		      if (! myDict)
			myDict = dictHandleCreate (12, h) ;
		      if (x > 0 && dictAdd (myDict, cp, 0))
			aceOutf (qc->ao, "%s%.0f %s"
				 , (k++ ? ", " : "")
				 , x, cp
				 ) ;
		    }
		}
	    }
	}
    }

   ac_free (h) ;
   return;
} /* qcHighVirusBacteria  */

/*************************************************************************************/

static void qcGeneExpression1 (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;

  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
 
    { "Compute1", "Number of %s genes touched\t%s genes significantly expressed\t%s genes with expression index >= 10 (sFPKM >= 1)\t%s genes with index >= 12  (sFPKM >= 4)\t%s genes with index >= 15 (sFPKM >= 32)\t%s genes with index >= 18 (sFPKM >= 256)\t%s genes with index >= 21 (sFPKM >= 2k)",  100, 0, 0 },
    /* %s genes with index >= 24 (sFPKM >= 16k)\t%s genes with index >= 27 (sFPKM >= 128k)\t%s genes with index >= 30 (sFPKM >= 1M) */
    { 0, 0, 0, 0, 0}
  }; 

  const char *caption =
    "Gene expression"
    ;

  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {	
      int ii ;
      char *cp ;
      const char *target ;
      AC_TABLE tt = 0 ;

      if (rc == 0)
	{
	  if (! strcmp (ti->tag, "Compute"))
	    {
	      char *buf = strnew (ti->title, h) ;
	      
	      cp = strstr (buf, "\t\t") ;
	      if (cp) *cp = 0 ;
	      for (ii = 0 ; ii < 3 ; ii++)
		{
		  target = qc->Etargets [ii] ;
		  aceOutf (qc->ao, "\t") ;
		  if (target && ! strcmp (target, "av")) target = "AceView" ;
		  aceOutf (qc->ao, buf, target ? target : "NA") ;
		}		
	    }
	  else if (! strcmp (ti->tag, "Compute1"))
	    {
	      const char *target2 ;
	      char *buf = strnew (ti->title, h) ;
	      target = qc->Etargets [0] ;
	      if (target && ! strcmp (target, "av")) target = "AceView" ;
	      aceOut (qc->ao, "\t") ;
	      target2 = target ;
	      if (qc->TT)
		target2 = hprintf (h, "Captured_%s", target) ;
	      aceOutf (qc->ao, buf
		       , target2
		       , target2
		       , target2
		       , target2
		       , target2
		       , target2
		       , target2
		       , target2
		       ) ;
	    }
	  else
	    aceOutf (qc->ao, "\t%s", ti->title) ;
	}
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute1"))
	{ 	  
	  const char *tag, *target = qc->Etargets [0] ;

	  tt = ac_tag_table (rc->ali, "Gene_expression" , h) ;
	  /*	  tag = "Targeted_genes" ; qcShowGeneExp (qc, tt, target, tag) ; */
	  tag = "Genes_touched" ; qcShowGeneExp (qc, tt, target, tag) ;
	  tag = "Genes_with_index" ; qcShowGeneExp (qc, tt, target, tag) ;
	  tag = "Genes_with_index_over_10" ; qcShowGeneExp (qc, tt, target, tag) ;
	  tag = "Genes_with_index_over_12" ; qcShowGeneExp (qc, tt, target, tag) ;
	  tag = "Genes_with_index_over_15" ; qcShowGeneExp (qc, tt, target, tag) ;
	  tag = "Genes_with_index_over_18" ; qcShowGeneExp (qc, tt, target, tag) ;
	  tag = "Genes_with_index_over_21" ; qcShowGeneExp (qc, tt, target, tag) ;
	  /*
	    tag = "Genes_with_index_over_24" ; qcShowGeneExp (qc, tt, target, tag) ;
	    tag = "Genes_with_index_over_27" ; qcShowGeneExp (qc, tt, target, tag) ;
	    tag = "Genes_with_index_over_30" ; qcShowGeneExp (qc, tt, target, tag) ;
	  
	    suppressed 6 columns 2025_02_06
	    \tReads in %s genes\tReads in pairs in %s genes\tOrphan reads in %s genes\tMulti-aligned reads in %s genes\tMismatches in %s genes\tA>G in %s genes

	    tag = "Reads_in_genes" ; qcShowGeneExp (qc, tt, target, tag) ;
	    tag = "Reads_in_pairs_in_genes" ; qcShowGeneExp (qc, tt, target, tag) ;
	    tag = "Orphan_reads_in_genes" ; qcShowGeneExp (qc, tt, target, tag) ;
	    tag = "Multi_aligned_reads_in_genes" ; qcShowGeneExp (qc, tt, target, tag) ;
	    tag = "Mismatches_in_genes" ; qcShowGeneExp (qc, tt, target, tag) ;
	    tag = "A2G_in_genes" ; qcShowGeneExp (qc, tt, target, tag) ;
	  */
	}
      else
	qcShowTag (qc, rc, ti) ;
    }

   ac_free (h) ;
   return;
} /* qcGeneExpression1 */

/*************************************************************************************/

static void qcGeneExpression2 (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;

  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,

    { "Compute1", "Zero index in %s\tLow index in %s\tCross over index in %s\tGenes expressed in at least one run of a group",  100, 0, 0 },

    { "Genes", "Genes significantly expressed in %s\t\t", 1, 0, 0 },
    { "Genes", "Transcripts significantly expressed in %s\t\t", 2, 0, 0} ,
    {  0, 0, 0, 0, 0}
  }; 

  const char *caption =
    "Gene and transcripts expression in various annotations"
    ;

  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {	
      int ii ;
      char *cp ;
      const char *target, *tag ;
      AC_TABLE tt = 0 ;

      if (rc == 0)
	{
	  if (! strcmp (ti->tag, "Genes"))
	    {
	      char *buf = strnew (ti->title, h) ;
	      
	      cp = strstr (buf, "\t\t") ;
	      if (cp) *cp = 0 ;
	      for (ii = 0 ; ii < 3 ; ii++)
		{
		  target = qc->Etargets [ii] ;
		  aceOutf (qc->ao, "\t") ;
		  if (target && ! strcmp (target, "av")) target = "AceView" ;
		  aceOutf (qc->ao, buf, target ? target : "NA") ;
		}		
	    }
 	  else if (! strcmp (ti->tag, "Compute1"))
	    {
	      char *buf = strnew (ti->title, h) ;
	      target = qc->Etargets [0] ;
	      if (target && ! strcmp (target, "av")) target = "AceView" ;
	      aceOut (qc->ao, "\t") ;
	      aceOutf (qc->ao, buf
		       , target
		       , target
		       , target
		       ) ;
	    }
	  else
	    aceOutf (qc->ao, "\t%s", ti->title) ;
	}
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute1"))
	{ 	  
	  const char *tag, *target = qc->Etargets [0] ;

	  tt = ac_tag_table (rc->ali, "Gene_expression" , h) ;
	  tag = "Zero_index" ; qcShowGeneExp (qc, tt, target, tag) ;
	  tag = "Low_index" ; qcShowGeneExp (qc, tt, target, tag) ;
	  tag = "Cross_over_index" ; qcShowGeneExp (qc, tt, target, tag) ;
	  tag = "Genes_expressed_in_at_least_one_run" ;
	  qcShowGeneExp (qc, tt, target, tag) ;
	}
      else if (! strcmp (ti->tag, "Genes"))
	{  
	  switch (ti->col)
	    {
	    case 1: 
	      tag = "Genes_with_index" ;
	      break ;
	    case 2: 
	      tag = "mRNA_with_index" ;
	      break ;
	    }
	  tt = ac_tag_table (rc->ali, tag, h) ;
	  if (tt)
	    {
	      int ii, ir ;
	      for (ii = 0 ; ii < 3 ; ii++)
		{
		  aceOut (qc->ao, "\t") ;
		  for (ir = 0 ; ir < tt->rows ; ir++)
		    {
		      target = qc->Etargets [ii] ;
		      if (target && ! strcasecmp (target, ac_table_printable (tt, ir, 0, "x")))
			{
			  aceOutf (qc->ao, "%d", ac_table_int (tt, ir, 1, 0)) ;
			  break ;
			}
		    }
		}
	    }
	  else
	    {
	      aceOut (qc->ao, "\t\t\t") ;
	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }

   ac_free (h) ;
   return;
} /* qcGeneExpression2 */

/*************************************************************************************/

static void qcHighGenes (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;

  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    
    { "Compute", "Mb in very highly expressed genes", 18, 0, 0 },
    { "HighGenes", "Very highly expressed genes, collecting over 2% of the reads", 10, 0, 0 },
    
    {  0, 0, 0, 0, 0}
  }; 
  
  const char *caption =
    "Highly expressed genes"
    ;
  
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;
  
  for (ti = tts ; ti->tag ; ti++)
    {	
      const char *target, *tag ;
      AC_TABLE tt = 0 ;
      
      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title) ;
	  continue ;
	}
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{  
	  switch (ti->col)
	    {
	    case 16:
	      tag = "Mb_in_genes" ;
	      break ;
	    case 18:
	      tag = "Mb_in_high_genes" ;
	      break ;
	    }
	  tt = ac_tag_table (rc->ali, tag, h) ;
	  if (tt)
	    {
	      int ii, ir ;
	      for (ii = 0 ; ii < (ti->col > 10 ? 1 : 3) ; ii++)
		{
		  const char *target2 ;
		  target = qc->Etargets [ii] ;
		  if (qc->TT)
		    target2 = hprintf (h,"Captured_%s", target) ;
		  else
		    target2 = target ;
		  for (ir = 0 ; ir < tt->rows ; ir++)
		    {
		      if (target2 && ! strcasecmp (target2, ac_table_printable (tt, ir, 0, "x")))
			break ;
		    }
		  if (ir < tt->rows)
		    {
		      switch (ti->col)	
			{
			case 16:
			case 18:
			  aceOutf (qc->ao, "\t%.3f", ac_table_float (tt, ir, 1, 0)) ;
			  break ;
			}
		    }
		  else
		    aceOutf (qc->ao, "\t") ;
		}
	    }
	  else
	    {
	      aceOut (qc->ao, "\t") ;
	    }
	}
      else if (! strcmp (ti->tag, "HighGenes"))
	{  
	  int ir, k = 0 ;
	  const char *target2 ;
	  target = qc->Etargets [0] ;
	  if (qc->TT)
	    target2 = hprintf (h,"Captured_%s", target) ;
	  else
	    target2 = target ;
	  aceOut (qc->ao, "\t") ;
	  tt = ac_tag_table (rc->ali, "High_genes", h) ;
	  for (ir = 0 ; tt && ir < tt->rows ; ir++)
	    {
	      if (! strcmp (ac_table_printable (tt, ir, 0, "x"), target2))
		{
		  AC_OBJ Gene = ac_table_obj (tt, ir, 1, h) ;
		  const char *ccp = Gene ? ac_tag_printable (Gene, "Title", ac_name(Gene)) : 0 ;
		  
		  if (ccp)
		    {
		      aceOutf (qc->ao, "%s%s"
			     , (k++ ? ", " : "")
			     , ccp
			     ) ;
		      ccp = ac_table_printable (tt, ir, 2, 0) ;
		      if (ccp)
			aceOutf (qc->ao, " (%sMb)", ccp) ;
		    }
		  ac_free (Gene) ;
		}
	    }
	}		  
      else
	qcShowTag (qc, rc, ti) ;
    }

   ac_free (h) ;
   return;
} /* qcHighGenes */

/*************************************************************************************/

static void qcRunSelection (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;

  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "SRP", "Project", 10, 0, 0} ,
    { "Ends", "Supported known starts", 41, 0, 0} , 
    { "Ends", "Candidate starts", 42, 0, 0} ,
    { "Ends", "Supported known ends", 43, 0, 0} ,
    { "Ends", "Candidate ends", 44, 0, 0} ,
    { "Introns", "Supported known exon-exon junctions", 10, 0, 0} ,
    { "PlusStrand", "% fragments on strand plus of genes", 90, 0, 0} ,  /* copied from qcGeneExpression */
    { "Accessible_length_8kb", "Homogeneous coverage of transcripts longer than 8kb, max 6kb, limited by the 3' bias in poly-A selected experiments", 0, 0, 0}, /* copied from qc3pBias */
    { "Compute", "% bases mapped in intergenic regions (suspected genomic contamination)\t% bases mapped to intronic regions (incompletely spliced transcripts)\t% bases mapped to UTR or non coding regions\t% bases mapped to protein coding regions (CDS)\t% Coding in mRNA: CDS/(CDS+UTR)\t% mapped raw bases in sponge", 1, 0, 0} ,

    { "AnyTarget", "% Mb aligned on any target", 54, 0, 0} , 
    { "Genes", "Genes significantly expressed in %s\t\t", 1, 0, 0 },
    { "Genes", "Genes with index >= 15 in %s", 15, 0, 0 },
    { "Genes", "Genes with Index >= 18 in %s", 18, 0, 0 },
    { "Paired_end", "Single or Paired end", 0, 0, 0} ,
    { "Fragment_length_median", "median insert size or average aligned length", 0, 0, 0} ,
    { "Pairs", "% non compatible pairs", 1, 0, 0} ,
    { "FragAli", "Average length aligned per fragment (nt)", 61, 0, 0 } , /* copied from qcAvLengthAli */
    { "SNV", "Exonic SNV sites\tPure variant SNV (95-100%)", 2, 0, 0} ,
    { "MbAli", "Mb uniquely aligned in main protein coding genes minus very highly expressed genes", 80, 0, 0 }, /* copied from qcGeneExpression */
    { "Ribo", "% Reads aligned in ribosomal RNA", 80, 0, 0 }, /* copied from qcGeneExpression */
    { "Any", "% Reads aligned in any target", 80, 0, 0 }, /* copied from qcGeneExpression */
    
    {  0, 0, 0, 0, 0}
  }; 
  
  const char *caption =
    "Run selection"
    ;
  
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;
  
  const char *target = qc->Etargets[0] ? qc->Etargets[0] : "xx" ;
  for (ti = tts ; ti->tag ; ti++)
    {	
      const char *tag ;
      AC_TABLE tt = 0 ;
      
      if (rc == 0)
	{
	  if (! strcmp (ti->tag, "Genes"))
	    {
	      char *buf = strnew (ti->title, h) ;
	      
	      char *cp = strstr (buf, "\t\t") ;
	      if (cp) *cp = 0 ;
	      for (int ii = 0 ; ii < (ti->col == 1 ? 3 : 1) ; ii++)
		{
		  target = qc->Etargets [ii] ;
		  aceOutf (qc->ao, "\t") ;
		  if (target && ! strcmp (target, "av")) target = "AceView" ;
		  aceOutf (qc->ao, buf, target ? target : "NA") ;
		}		
	    }
	  else
	    aceOutf (qc->ao, "\t%s", ti->title) ;
	  continue ;
	}
       else if (! strcmp (ti->tag, "SRP"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "SRP", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Ends"))
	{
	  int k = 0 ;
	  switch (ti->col)
	    {
	    case 41: tag = "Confirmed_starts" ; break ;
	    case 42: tag = "Candidate_starts" ; break ;
	    case 43: tag = "Confirmed_ends" ; break ;
	    case 44: tag = "Candidate_ends" ; break ;
	    }
	  k = ac_tag_int (rc->ali, tag, 0) ;
	  aceOutf (qc->ao, "\t%d", k) ;
	}
      else if (! strcmp (ti->tag, "Introns"))
	{
	  int ic, ir ;
	  float z = 0 ;
	  tag = "Candidate_introns" ;  ic = 3 ; 
	  tt = ac_tag_table (rc->ali, tag, h) ;
	  z = 0 ;
	  if (tt)
	    for (ir = 0 ; ir < tt->rows ; ir++)
	      {
		if (!ir || ! strcasecmp (ac_table_printable (tt, ir, 0, "xxx"), target))
		  z = ac_table_int (tt, ir, ic, 0) + ac_table_float (tt, ir, ic, 0) ;
	      }
	  aceOutf (qc->ao, "\t%.0f", z) ;
	}
      else if (! strcmp (ti->tag, "PlusStrand"))
	{  /* Observed strandedness */
	  float z, z1 ;
	  tt = ac_tag_table (rc->run, "Observed_strandedness_in_ns_mapping", h) ; 
	  z = 50 ;
	  if (tt)
	    for (int ir = 0 ; ir < tt->rows ; ir++)
	      {
		z1 = ac_table_float (tt, 0,1,50) ;
		if ((z-50)*(z-50) < (z1-50)*(z1-50))
		  z = z1 ;
	      }
	  aceOutf (qc->ao, "\t%.2f", z) ;
	  ac_free (tt) ;
	}
      else if (! strcmp (ti->tag, "AnyTarget"))
	{
	  const char *ccp = EMPTY ; tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
	  float zb = 0 ;
	  for (int ir = 0 ; tt && ir < tt->rows ; ir++)
	    {
	      ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
	      if (! strcasecmp (ccp, "any"))
		{
		  /* z = ac_table_float (tt, ir, 3, 0) ; */
		  zb = ac_table_float (tt, ir, 5, 0) ;
		  /* zc = ac_table_float (tt, ir, 7, 0) ; */
		}
	    } 
	  
	  aceOutf (qc->ao, "\t%.2f", rc->var[T_kb] ? 100 *zb / rc->var[T_kb] : 0) ; /* % Mb aligned on any target after clipping */
	  ac_free (tt) ;
	}
      else if (! strcmp (ti->tag, "Ribo"))
	{
	  const char *ccp = EMPTY ; tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
	  float z = 0, zb = 0 ;
	  for (int ir = 0 ; tt && ir < tt->rows ; ir++)
	    {
	      ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
	      if (! strcasecmp (ccp, "B_rrna"))
		{
		  if (1)  z = ac_table_float (tt, ir, 3, 0) ; 
		  if (0) zb = ac_table_float (tt, ir, 5, 0) ;
		  /* zc = ac_table_float (tt, ir, 7, 0) ; */
		}
	    } 
	  
	  if (0)  aceOutf (qc->ao, "\t%.2f", rc->var[T_kb] ? 100 *zb / rc->var[T_kb] : 0) ; /* % Mb aligned in rRNA */
	  if (1)  aceOutf (qc->ao, "\t%.2f", rc->var[T_Read] ? 100 *z / rc->var[T_Read] : 0) ; /* % Read aligned in rRNA */
	  ac_free (tt) ;
	}
      else if (! strcmp (ti->tag, "Any"))
	{
	  const char *ccp = EMPTY ; tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
	  float z = 0, zb = 0 ;
	  for (int ir = 0 ; tt && ir < tt->rows ; ir++)
	    {
	      ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
	      if (! strcasecmp (ccp, "any"))
		{
		  if (1)  z = ac_table_float (tt, ir, 3, 0) ; 
		  if (0) zb = ac_table_float (tt, ir, 5, 0) ;
		  /* zc = ac_table_float (tt, ir, 7, 0) ; */
		}
	    } 
	  
	  if (0)  aceOutf (qc->ao, "\t%.2f", rc->var[T_kb] ? 100 *zb / rc->var[T_kb] : 0) ; /* % Mb aligned in any */
	  if (1)  aceOutf (qc->ao, "\t%.2f", rc->var[T_Read] ? 100 *z / rc->var[T_Read] : 0) ; /* % Read aligned in any */
	  ac_free (tt) ;
	}
      else if (! strcmp (ti->tag, "Genes"))
	{
	  switch (ti->col)
	    {
	    case 1: 
	      tag = "Genes_with_index" ;
	      break ;
	    case 15: 
	      tag = "Genes_with_index_over_15" ;
	      break ;
	    case 18: 
	      tag = "Genes_with_index_over_18" ;
	      break ;
	    }
	  tt = ac_tag_table (rc->ali, tag, h) ;
	  if (1)
	    {
	      int ii, ir ;
	      for (ii = 0 ; ii < (ti->col == 1 ? 3 : 1) ; ii++)
		{
		  aceOutf (qc->ao, "\t") ;
		  if (tt)
		    {
		      for (ir = 0 ; ir < tt->rows ; ir++)
			{
			  target = qc->Etargets [ii] ;
			  if (target && ! strcasecmp (target, ac_table_printable (tt, ir, 0, "x")))
			    break ;
			}
		      if (ir < tt->rows)
			aceOutf (qc->ao, "%d", ac_table_int (tt, ir, 1, 0)) ;
		    }
		}
	    }
	}
      else if (! strcmp (ti->tag, "Paired_end"))
	{ 
	  aceOutf (qc->ao, "\t") ;
	  if (ac_has_tag (rc->run, "Is_run"))
	    aceOutf (qc->ao, rc->Paired_end ? "Paired_end" : "Single-end") ;
	  else
	    {
	      int n1, n2 ;
	      n1 = ac_keyset_count (ac_objquery_keyset (rc->run, ">Union_of ;>Runs ;   Paired_end", h)) ;
	      n2 = ac_keyset_count (ac_objquery_keyset (rc->run, ">Union_of ; ! Paired_end", h)) ;
	      if (n1 > 0 && n2 > 0)
		aceOutf (qc->ao, "Mix") ;
	      else if (n1 > 0 && n2 == 0)
		aceOutf (qc->ao, "Paired_end") ;
	      else if (n1 == 0 && n2 > 0)
		aceOutf (qc->ao, "Single-end") ;
	    }
	}
      else if (! strcmp (ti->tag, "Pairs"))
	{ /* percent compatible pairs */
	  float z1 = ac_tag_float (rc->ali, "Aligned_fragments", 0) ;
	  float z2 = ac_tag_float (rc->ali, "Non_compatible_pairs", 0) ;
	  aceOutf (qc->ao, "\t%.2f", 100.0 * z2/(z1 + .1)) ;
	}
      else if (! strcmp (ti->tag, "FragAli"))
	{
	  const char *ccp = EMPTY ;
	  tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
	  float z = 0, zb = 0, zc = 0 ;
	  if (tt)
	    for (int ir = 0 ; tt && ir < tt->rows ; ir++)
	      {
		ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		if (! strcasecmp (ccp, "any1") || ! strcasecmp (ccp, "any2"))
		  {
		    zb += ac_table_float (tt, ir, 5, 0) ;
		  }
	      }
	  zc = ac_tag_float (rc->ali,  "Aligned_fragments", 0) ;
	  if (zc > 0) z = 1000.0 * zb / zc ;
	  aceOutf (qc->ao, "\t%.2f", z) ;  /* average length aligned per fragment */
	  ac_free (tt) ;
	}
      else if (! strcmp (ti->tag, "SNV"))
	{
	  tt = ac_tag_table (rc->ali, "Genomic", h) ;
	  if (tt)
	    {
	      int n1 = ac_table_int (tt, 0, 0, 0) ;
	      int n2 = ac_table_int (tt, 0, 10, 0) ;
	      aceOutf (qc->ao, "\t%d\t%d", n1, n2) ;
	    }
	  else
	    aceOut (qc->ao, "\t\t") ;
	}
      else if (! strcmp (ti->tag, "MbAli"))
	{
	  tag = "Mb_in_genes_with_GeneId_minus_high_genes" ;
	  aceOut (qc->ao, "\t") ;
	  tt = ac_tag_table (rc->ali, tag, h) ;
	  if (tt)
	    {
	      int ir ;
	      for (ir = 0 ; ir < tt->rows ; ir++)
		{
		  if (target && ! strcasecmp (target, ac_table_printable (tt, ir, 0, "x")))
		    break ;
		}
	      if (ir < tt->rows)
		{
		  switch (ti->col)	
		    {
		    case 80:   
		      aceOutf (qc->ao, "%.03f", ac_table_float (tt, ir, 1, 0)) ;
		      break ;
		    }
		}
	    }
	  ac_free (tt) ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  float z, zG = 0, zI = 0, zU = 0, zC = 0, zRaw = 0 ;
	  switch (ti->col)
	    {
	    case 1:
	      zG = ac_tag_float (rc->ali, "S_1_intergenic", 0) ;
	      zI = ac_tag_float (rc->ali, "S_1_intronic", 0) ;
	      zU = ac_tag_float (rc->ali, "S_1_UTR", 0) ;
	      zC = ac_tag_float (rc->ali, "S_1_CDS", 0) ;
	      z = zG + zI + zU + zC ;
  	      zRaw = rc->var[T_kb] ? rc->var[T_kb] : z/2 ;
	      if (z > 0)
		aceOutf (qc->ao, "\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f"
			 , 100.0 * zG/z
			 , 100.0 * zI/z
			 , 100.0 * zU/z
			 , 100.0 * zC/z
			 , 100.0 * (zC/(zC+zU+.01))
			 , 100000.0 * z/zRaw
			 ) ;
	      else
		aceOutf (qc->ao, "\t\t\t\t\t\t") ;
	      break ;
	    }
	}
      else if (! strcmp (ti->tag, "Fragment_length_median"))
	{
	  if (ac_has_tag (rc->ali, ti->tag))
	    qcShowTag (qc, rc, ti) ;
	  else
	    {
	      float z = -1 ;
	      aceOut (qc->ao, "\t") ;
	      tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
	      for (int ir = 0 ; tt && ir < tt->rows ; ir++)
		{
		  if (! strcmp ("any", ac_table_printable (tt, ir, 0, "toto")))
		    {
		      z = ac_table_float (tt, ir, 11, 0) ;
		      if (z >= 0)
			aceOutf (qc->ao, "%.0f", z) ;
		      break ;
		    }
		}
	    }
	}
      else /* 110: accessible length 8k/5k, median fragment length */
	qcShowTag (qc, rc, ti) ;
    }

   ac_free (h) ;
   return;
} /* qcRunSelection */

/*************************************************************************************/

static void qcCandidateIntrons (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;

  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Candidate__introns", "Supported %s exon-exon junctions\t\t", 3, 0, 0 },
    { "Candidate__introns", "%% %s exon-exon junctions with support\t\t", 9, 0, 0 },
    { "Candidate__introns", "Number of support for %s exon-exon junctions\t\t", 13, 0, 0 },
    { "Candidate__introns", "Average support of %s exon-exon junctions\t\t", 1000, 0, 0 },
    { "Candidate_introns", "Candidate new exon-exon junctions", 5, 0, 0 },
    { "Candidate_introns", "Minimal support for new exon-exon junctions", 17, 0, 0 },
    {  0, 0, 0, 0, 0}
  }; 
  
  const char *caption =
    "Candidate introns"
    ;
  AC_TABLE ttCI = 0 ;

  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;
  
  for (ti = tts ; ti->tag ; ti++)
    {	
      int ii ;
      char *cp ;
      const char *target ;
      
      if (rc == 0)
	{
	  if (! strcmp (ti->tag, "Candidate__introns"))
	    {
	      char *buf = strnew (ti->title, h) ;
	      
	      cp = strstr (buf, "\t\t") ;
	      if (cp) *cp = 0 ;
	      for (ii = 0 ; ii < 3 ; ii++)
		{
		  target = qc->Etargets [ii] ;
		  aceOutf (qc->ao, "\t") ;
		  if (target && ! strcmp (target, "av")) target = "AceView" ;
		  aceOutf (qc->ao, buf, target ? target : "NA") ;
		}		
	    }
	  else
	    aceOutf (qc->ao, "\t%s", ti->title) ;
	}
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Candidate__introns"))
	{  
	  int ir ;
	  if (! ttCI)
	    ttCI = ac_tag_table (rc->ali, "Candidate_introns", h) ;
	  for (ii = 0 ; ii < 3 ; ii++)
	    {
	      target = qc->Etargets [ii] ;
	      aceOutf (qc->ao, "\t") ;
	      if (! strcmp (target, "av")) target = "AceView" ;
	      if (ttCI && target && *target)
		{
		  for (ir = 0 ; ir < ttCI->rows ; ir++)
		    if (! strcasecmp (target, ac_table_printable (ttCI, ir, 0, "xxx"))) 
		      {
			if (ti->col == 1000)
			  {
			    float u13 = ac_table_float (ttCI, ir, 13, 0) ;
			    int u3 = ac_table_int (ttCI, ir,  3, 1) ;
			    if (u3 == 0) u3 = 1 ;
			    aceOutf (qc->ao, "%.1f", u13/u3) ;
			  }
			else
			  aceOut (qc->ao, ac_table_printable (ttCI, ir, ti->col, EMPTY)) ;
		      }
		}
	    }	
	}
      else if (! strcmp (ti->tag, "Candidate_introns"))
	{  
	  int ir ;
	  if (! ttCI)
	    ttCI = ac_tag_table (rc->ali, "Candidate_introns", h) ;
	  for (ii = 0 ; ii < 1 ; ii++)
	    {
	      target = "any" ;
	      aceOutf (qc->ao, "\t") ;
	      
	      if (ttCI && target && *target)
		{
		  for (ir = 0 ; ir < ttCI->rows ; ir++)
		    if (! strcmp (target, ac_table_printable (ttCI, ir, 0, "xxx"))) 
		      aceOut (qc->ao, ac_table_printable (ttCI, ir, ti->col, EMPTY)) ;
		}
	    }	
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  
  ac_free (h) ;
  return;
} /* qcCandidateIntrons */

/*************************************************************************************/

static void qcSexTissueSignatures (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;

  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
 
    { "Compute1", "Differential sex index in %s\tSex signature in %s",  100, 0, 0 },

    {  0, 0, 0, 0, 0}
  }; 

  const char *caption =
    "Signatures"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {	
      const char *target, *tag ;
      AC_TABLE tt = 0 ;

      if (rc == 0)
	{
	  if (! strcmp (ti->tag, "Compute1"))
	    {
	      char *buf = strnew (ti->title, h) ;
	      target = qc->Etargets [0] ;
	      if (target && ! strcmp (target, "av")) target = "AceView" ;
	      aceOut (qc->ao, "\t") ;
	      aceOutf (qc->ao, buf
		       , target
		       , target
		       ) ;
	    }
	  else
	    aceOutf (qc->ao, "\t%s", ti->title) ;
	}
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute1"))
	{ 	  
	  float z = 0 ;
	  target = qc->Etargets [0] ;

	  tt = ac_tag_table (rc->ali, "Gene_expression" , h) ;
	  tag = "Sex_ratio" ; z = qcShowGeneExp (qc, tt, target, tag) ;
	  if (z < -1) aceOut (qc->ao, "\tFemale") ;
	  else if (z > 4) aceOut (qc->ao, "\tMale") ;
	  else  aceOut (qc->ao, "\t") ;
	}
      else
	qcShowTag (qc, rc, ti) ;
    }

   ac_free (h) ;
   return;
} /* qcSexTissueSignatures */

/*************************************************************************************/

static void qcMappingPerTargetType (QC *qc, RC *rc, int type)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt, ttk ;  
  const char *ccp ;
  float z ;
  long int zhPrevious ;
  long int znhA, znhR, znhE, znhG, znhg, zhBacteria, zhVirus, zSpliced ;
  long int zhe, zhPhiX, zhr, zhm, zhS, zhT, zhA, zhR, zhE, zhG, zhg, zhAny, zhTotal ;
  int ir ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} , 
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute", 
      "%s aligned in microbiome or symbiome"
      "\t%s aligned in viruses"

      "\t%s aligned in PhiX DNA spike-in"
      "\t%s aligned in RNA spike-in (ERCC, Sequin, or yeast enolase)"
      "\t%s aligned in ribosomal RNA"
      "\t%s aligned in mitochondria"
      "\t%s aligned in other small genes"
      "\t%s aligned in transposons"

      "\t%s aligned in %s"
      "\t%s aligned in %s"
      "\t%s aligned in %s"
      "\t%s aligned in genome"

      "\t%s best aligned in any previously annotated transcripts"
      "\t%s inside introns, new exons and intergenic"

      "\t%s aligned on the imaginary genome (mapping specificity control)"
      /* "\t%s supporting known exon junctions" */
      , 1, 0, 0} ,
    {  0, 0, 0, 0, 0}
  }; 

  const char *caption = 
    "Mapping: "
    "In Magic, all reads are mapped in parallel to all targets, nuclear ribosomal (18S, 28S, 5S, 5.8S),  "
    "mitochondrial,  various alternate transcriptomes (e.g AceView public, AceView next, RefSeq, UCSC, Ensembl...),  "
    "the  reference genome, the imaginary genome as a mapping specificity control, the RNA spike in controls  "
    "(such as ERCC) and the DNA Spike in controls (such as PhiX in Illumina experiments).  "
    "  All reads mapping at identical best score in multiple targets are counted and reported in the best target tables.  "
    "  In the hierarchical columns, each read counts only in the first encountered column where it aligns at best score.  "
    "For instance, the counts reported in the 'genome' column map better in the genome than in any transcriptome or  "
    "organelle, and hence corresponds to new exonic intergenic and intronic hits. "
    ;

  switch (type)
    {
    case 1: 
      caption = "Reads mapping in genes, genome and other targets" ;
      break ;
    case 2: 
      caption = "% Reads mapping in genes, genome and other targets" ;
      break ;
    case 3: 
      caption = "Megabases aligned in genes, genome and other targets" ;
      break ;
    }   
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	{
	  char *tit = hprintf (h, "\t%s", ti->title) ;
	  switch (type)
	    {
	    case 1:
	      aceOutf (qc->ao, tit
		       , "Reads"
		       , "Reads"
		       , "Reads"
		       , "Reads"
		       , "Reads"
		       , "Reads"
		       , "Reads"
		       , "Reads"
		       , "Reads", qc->EtargetsBeau[0]
		       , "Reads", qc->EtargetsBeau[1]
		       , "Reads", qc->EtargetsBeau[2]
		       , "Reads"
		       , "Reads"
		       , "Reads"
		       , "Reads"
		       , "Reads"
		       ) ;
	      break ;
	    case 2:
	      aceOutf (qc->ao,tit
		       , "% Reads"
		       , "% Reads"
		       , "% Reads"
		       , "% Reads"
		       , "% Reads"
		       , "% Reads"
		       , "% Reads"
		       , "% Reads"
		       , "% Reads", qc->EtargetsBeau[0]
		       , "% Reads", qc->EtargetsBeau[1]
		       , "% Reads", qc->EtargetsBeau[2]
		       , "% Reads"
		       , "% Reads"
		       , "% Reads"
		       , "% Reads"
		       , "% Reads"
		       ) ;
	      break ;
	    case 3: 
	      aceOutf (qc->ao, tit
		       , "Mb"
		       , "Mb"
		       , "Mb"
		       , "Mb"
		       , "Mb"
		       , "Mb"
		       , "Mb"
		       , "Mb"
		       , "Mb", qc->EtargetsBeau[0]
		       , "Mb", qc->EtargetsBeau[1]
		       , "Mb", qc->EtargetsBeau[2]
		       , "Mb"
		       , "Mb"
		       , "Mb"
		       , "Mb"
		       , "Mb"
		       ) ;
	      break ;
	    }
	}
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      
      else if (! strcmp (ti->tag, "Compute"))
	{
	  int col = (type == 3 ? 5 : 3) ;
	  switch (ti->col)
	    {
	    case 1:
	      ccp = EMPTY ;
	      zhe = zhr = zhm = zhA = zhR = zhE = zhG = zhg = zhAny = zhS = zhT = zhPhiX = zSpliced = zhBacteria = zhVirus = zhPrevious = zhTotal = 0 ;
	      ttk = ac_tag_table (rc->ali, "Known_introns", h) ;
	      for (ir = 0 ; ttk && ir < ttk->rows ; ir++)
		{
		  if (! strcmp ( ac_table_printable (ttk, ir, 0, "x"), qc->Etargets[0]))
		    zSpliced = ttk ? ac_table_float (ttk, ir, 3, 0) : 0 ;
		}
	      tt = ac_tag_table (rc->ali, "h_Ali", h) ;
	      for (ir = 0 ; tt && ir < tt->rows ; ir++)
		{
		  ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		  if (! strcasecmp (ccp, "1_DNASpikeIn"))
		    zhPhiX = ac_table_float (tt, ir, col, 0) ;
		  if (! strcasecmp (ccp, "0_SpikeIn"))
		    zhe = ac_table_float (tt, ir, col, 0) ;
		  if (! strcasecmp (ccp, "B_rRNA"))
		    zhr = ac_table_float (tt, ir, col, 0) ;
		  if (! strcasecmp (ccp, "A_mito"))
		    zhm = ac_table_float (tt, ir, col, 0) ;
		  if (! strcasecmp (ccp, "D_transposon"))
		    zhT = ac_table_float (tt, ir, col, 0) ;

		  if (! strcasecmp (ccp, "QT_smallRNA"))
		    zhS = ac_table_float (tt, ir, col, 0) ;
		  if (! strcasecmp (ccp + 3, qc->Etargets[0]))
		    zhA = ac_table_float (tt, ir, col, 0) ;
		  if (! strcasecmp (ccp + 3, qc->Etargets[1]))
		    zhR = ac_table_float (tt, ir, col, 0) ;
		  if (! strcasecmp (ccp + 3, qc->Etargets[2]))
		    zhE = ac_table_float (tt, ir, col, 0) ;
		  if (ccp[1] == 'T' && ccp[2] == '_')   /* any xT_ is an annotated transcript */
		    zhPrevious += ac_table_float (tt, ir, col, 0) ;

		  if (! strcasecmp (ccp, "v_virus"))
		    zhVirus = ac_table_float (tt, ir, col, 0) ;
		  if (! strcasecmp (ccp, "b_bacteria"))
		    zhBacteria = ac_table_float (tt, ir, col, 0) ;

		  if (! strcasecmp (ccp, "Z_genome"))
		    zhG = ac_table_float (tt, ir, col, 0) ;
		  if (! strcasecmp (ccp, "z_gdecoy"))
		    zhg = ac_table_float (tt, ir, col, 0) ;
		  if (! strcasecmp (ccp, "Any"))
		    zhAny = ac_table_float (tt, ir, col, 0) ;
		}
	      zhTotal = zhPhiX + zhe + zhr + zhm + zhT + zhPrevious + zhVirus + zhBacteria + zhG + zhg ;
	      zhPrevious += zhe + zhr + zhm + zhT  ; /* previously annotated transcripts */ 

	      ccp = EMPTY ; tt = ac_tag_table (rc->ali, "nh_Ali", h) ;
	      znhA = znhR = znhE = znhG = znhg = 0 ;
	      for (ir = 0 ; tt && ir < tt->rows ; ir++)
		{
		  ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
		  if (! strcasecmp (ccp + 3, qc->Etargets[0]))
		    znhA = ac_table_float (tt, ir, col, 0) ;
		  if (! strcasecmp (ccp + 3, qc->Etargets[1]))
		    znhR = ac_table_float (tt, ir, col, 0) ;
		  if (! strcasecmp (ccp + 3, qc->Etargets[2]))
		    znhE = ac_table_float (tt, ir, col, 0) ;
		  if (! strcasecmp (ccp, "Z_genome"))
		    znhG = ac_table_float (tt, ir, col, 0) ;
		  if (! strcasecmp (ccp, "z_gdecoy"))
		    znhg = ac_table_float (tt, ir, col, 0) ;
		} 

	      switch (type)
		{
		case 1:
		  aceOutf (qc->ao, "\t%ld", zhBacteria) ;	 
		  aceOutf (qc->ao, "\t%ld", zhVirus) ;	 
		  aceOutf (qc->ao, "\t%ld", zhPhiX) ;	      
		  aceOutf (qc->ao, "\t%ld", zhe) ;
		  aceOutf (qc->ao, "\t%ld", zhr) ;
		  aceOutf (qc->ao, "\t%ld", zhm) ;
		  aceOutf (qc->ao, "\t%ld", zhS) ;
		  aceOutf (qc->ao, "\t%ld", zhT) ;
		  
		  aceOutf (qc->ao, "\t%ld", znhA) ;
		  aceOutf (qc->ao, "\t%ld", znhR) ;
		  aceOutf (qc->ao, "\t%ld", znhE) ;

		  aceOutf (qc->ao, "\t%ld", znhG) ;
		  
		  aceOutf (qc->ao, "\t%ld", zhPrevious) ;
		  aceOutf (qc->ao, "\t%ld", zhAny - zhPrevious - zhPhiX - zhBacteria -zhVirus - zhg) ; /* intronic, new exon and intergenic */
		  
		  
		  aceOutf (qc->ao, "\t%ld", znhg) ;	 /* decoy */ 
		  
		  /* aceOutf (qc->ao, "\t%ld", zSpliced) ; */
		  break ;
		case 3:
		  aceOutf (qc->ao, "\t%.2f", zhBacteria/1000.0) ;	 
		  aceOutf (qc->ao, "\t%.2f", zhVirus/1000.0) ;	 
		  aceOutf (qc->ao, "\t%.2f", zhPhiX/1000.0) ;	      
		  aceOutf (qc->ao, "\t%.2f", zhe/1000.0) ;
		  aceOutf (qc->ao, "\t%.2f", zhr/1000.0) ;
		  aceOutf (qc->ao, "\t%.2f", zhm/1000.0) ;
		  aceOutf (qc->ao, "\t%.2f", zhS/1000.0) ;
		  aceOutf (qc->ao, "\t%.2f", zhT/1000.0) ;
		  
		  aceOutf (qc->ao, "\t%.2f", znhA/1000.0) ;
		  aceOutf (qc->ao, "\t%.2f", znhR/1000.0) ;
		  aceOutf (qc->ao, "\t%.2f", znhE/1000.0) ;

		  aceOutf (qc->ao, "\t%.2f", znhG/1000.0) ;
		  aceOutf (qc->ao, "\t%.2f", zhPrevious/1000.0) ; 
		  z = zhAny - zhPrevious - zhPhiX - zhBacteria -zhVirus - zhg  ;
		  aceOutf (qc->ao, "\t%.2f", z/1000.0) ; /* intronic, new exon and intergenic */
		  
		  
		  aceOutf (qc->ao, "\t%.2f", znhg/1000.0) ;	 /* decoy */ 
		  
		  /* aceOutf (qc->ao, "\t") ; */
		  break ;
		case 2: 
		  z = zhBacteria ; aceOutf (qc->ao, "\t") ; z = 100 * z/zhTotal ; aceOutPercent (qc->ao, z) ;
		  z = zhVirus ; aceOutf (qc->ao, "\t") ; z = 100 * z/zhTotal ; aceOutPercent (qc->ao, z) ;
		  z = zhPhiX ; aceOutf (qc->ao, "\t") ; z = 100 * z/zhTotal ; aceOutPercent (qc->ao, z) ;
		  z = zhe ; aceOutf (qc->ao, "\t") ; z = 100 * z/zhTotal ; aceOutPercent (qc->ao, z) ;
		  z = zhr ; aceOutf (qc->ao, "\t") ; z = 100 * z/zhTotal ; aceOutPercent (qc->ao, z) ;
		  z = zhm ; aceOutf (qc->ao, "\t") ; z = 100 * z/zhTotal ; aceOutPercent (qc->ao, z) ;
		  z = zhS ; aceOutf (qc->ao, "\t") ; z = 100 * z/zhTotal ; aceOutPercent (qc->ao, z) ;
		  z = zhT ; aceOutf (qc->ao, "\t") ; z = 100 * z/zhTotal ; aceOutPercent (qc->ao, z) ;

		  z = znhA ; aceOutf (qc->ao, "\t") ; z = 100 * z/zhTotal ; aceOutPercent (qc->ao, z) ;
		  z = znhR ; aceOutf (qc->ao, "\t") ; z = 100 * z/zhTotal ; aceOutPercent (qc->ao, z) ;
		  z = znhE ; aceOutf (qc->ao, "\t") ; z = 100 * z/zhTotal ; aceOutPercent (qc->ao, z) ;
		  z = znhG ; aceOutf (qc->ao, "\t") ; z = 100 * z/zhTotal ; aceOutPercent (qc->ao, z) ;

		  z = zhPrevious ; aceOutf (qc->ao, "\t") ; z = 100 * z/zhTotal ; aceOutPercent (qc->ao, z) ;
		  z = zhTotal - zhPrevious - zhPhiX - zhBacteria -zhVirus - zhg  ; aceOutf (qc->ao, "\t") ; z = 100 * z/zhTotal ; aceOutPercent (qc->ao, z) ; /* intronic, new exon and intergenic */
	
		  z = znhg ; aceOutf (qc->ao, "\t") ; z = 100 * z/zhTotal ; aceOutPercent (qc->ao, z) ; /* decoy */

		  /*   z = zSpliced ; aceOutf (qc->ao, "\t") ; z = 100 * z/zhTotal ; aceOutPercent (qc->ao, z) ; */
		  break ;
		}
	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  ac_free (h) ;
  return;
}  /* qcMappingPerTargetType */

/******************/

static void qcMappingPerTarget (QC *qc, RC *rc)
{
  qcMappingPerTargetType (qc, rc, 1) ; /* number of reads */
  if (0)  qcMappingPerTargetType (qc, rc, 2) ; /* percent number of reads */
  qcMappingPerTargetType (qc, rc, 3) ; /* number of kb */

  return ;
} /* qcMappingPerTarget */

/******************/

static void qcMappingPerCentTarget (QC *qc, RC *rc)
{
  if (0) qcMappingPerTargetType (qc, rc, 1) ; /* number of reads */
  if (0) qcHighVirusBacteria (qc, rc) ;
  if (1)  qcMappingPerTargetType (qc, rc, 2) ; /* percent number of reads */
  if (0) qcMappingPerTargetType (qc, rc, 3) ; /* number of kb */

  return ;
} /* qcMappingPerCentPerTarget */

/*************************************************************************************/

static void qc3pBias (QC *qc, RC *rc)
{
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Accessible_length_8kb", "Homogeneous coverage of transcripts longer than 8kb, max 6kb, limited by the 3' bias in poly-A selected experiments", 0, 0, 0}, /* copied from qc3pBias */
    { "Accessible_length", "Number of well expressed transcripts longer than 8kb", 1, 0, 0} ,
    { "Accessible_length_8kb", "Average coverage cumulated over these transcripts", 3, 0, 0} ,
    {  0, 0, 0, 0, 0}
  }; 
  const char *caption =
    "3' bias"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title ? ti->title : "") ;
	  if (0 && rc->ali && ti->title && ! strcmp (ti->title, "Number of well expressed transcripts longer than 8kb"))
	    {
	      AC_HANDLE h = ac_new_handle () ;
	      AC_TABLE tt = ac_tag_table (rc->ali, ti->tag, h) ;
	      if (tt && tt->cols >= 3)
		{
		  int n = ac_table_int (tt, 0, 3, 0) ;
		  if (n > 0)
		    aceOutf (qc->ao, ", out of %d candidates", n) ;
		}
	      ac_free (h) ;
	    }
	}
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  return;
}  /* qc3pBias */

/*************************************************************************************/

static void qcInterGenic (QC *qc, RC *rc)
{
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute", "Intergenic kb\tAverage coverage of intergenic region (suspected genomic contamination)", 1, 0, 0} ,
    {  0, 0, 0, 0, 0}
  }; 
  const char *caption =
    "Intergenic, used to estimate background noise"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title ? ti->title : "") ;
	}
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  float zd = 0, zkb = 0 ;
	  switch (ti->col)
	    {
	    case 1:
	      zd = ac_tag_float (rc->ali, "Intergenic_density", 0) ;
	      zkb = ac_tag_float (rc->ali, "Intergenic", 0) ;
	      if (zd > 0) 
		{
		  float z1 = zd ;
		  int n = 0 ;
		  char *f ;
		  
		  while (z1 < 100) { n++ ; z1 *= 10 ; }
		  if (n < 2) n = 2 ;
		  f = messprintf ("\t%%.0f\t%%.%df", n) ;
		  aceOutf (qc->ao, f, zkb, zd) ;
		}
	      else
		aceOutf (qc->ao, "\t\t") ;
	      break ;
	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  return;
}  /* qcInterGenic */

/*************************************************************************************/

static void qcIntergenic2 (QC *qc, RC *rc)
{
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    /*     { "TITLE", "Title", 10, 0, 0} , */
    { "Compute", "% bases mapped in intergenic regions (suspected genomic contamination)\t% bases mapped to intronic regions (incompletely spliced transcripts)\t% bases mapped to UTR or non coding regions\t% bases mapped to protein coding regions (CDS)\t% Coding in mRNA: CDS/(CDS+UTR)\t% mapped raw bases in sponge (no rRNA)", 1, 0, 0} ,
    {  0, 0, 0, 0, 0}
  }; 
  const char *caption =
    "Intergenic/Intronic/UTR or non-coding/Protein coding region"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title ? ti->title : "") ;
	}
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  float z, zG = 0, zI = 0, zU = 0, zC = 0, zRaw = 0 ;
	  switch (ti->col)
	    {
	    case 1:
	      zG = ac_tag_float (rc->ali, "S_1_intergenic", 0) ;
	      zI = ac_tag_float (rc->ali, "S_1_intronic", 0) ;
	      zU = ac_tag_float (rc->ali, "S_1_UTR", 0) ;
	      zC = ac_tag_float (rc->ali, "S_1_CDS", 0) ;
	      z = zG + zI + zU + zC ;
  	      zRaw = rc->var[T_kb] ? rc->var[T_kb] : z/2 ;
	      if (z > 0)
		aceOutf (qc->ao, "\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f"
			 , 100.0 * zG/z
			 , 100.0 * zI/z
			 , 100.0 * zU/z
			 , 100.0 * zC/z
			 , 100.0 * (zC/(zC+zU+.01))
			 , 100000.0 * z/zRaw
			 ) ;
	      else
		aceOutf (qc->ao, "\t\t\t\t\t\t") ;
	      break ;
	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }
  return;
}  /* qcIntergenic2 */

/*************************************************************************************/
#ifdef JUNK

static void qcDrosoZhenXia (QC *qc, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt, tta, ttb ;
  float zRaw = 0 ;

  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute", "Total reads\tZ Total reads\t % Total reads\t % Z Total reads\t% difference", 1, 0, 0} ,
    { "Compute", "Mapped\tZ Mapped\t% Mapped\t% Z Mapped\t% difference", 2, 0, 0} ,
    { "Compute", "Mapped uniquely\tZ Mapped uniquely\t% Mapped uniquely\t% Z Mapped uniquely\t% difference", 3, 0, 0} ,
    { "Compute", "Mapped non uniquely\tZ Mapped non_unique\t% Mapped non uniquely\t% Z non_uniquely\t% difference", 30, 0, 0} ,
    { "Compute", "Unmapped\tZ unmapped_reads\t% Unmapped\t% Z unmapped_reads\t% difference", 31, 0, 0} ,
    { "Compute", "Intron support: a read supporting 2 introns counts 2\tZ spliced read\t% Intron supports: a read supporting 2 intron counts 2\t% Z spliced read\t% difference", 4, 0, 0} ,
    { "Compute", "Reads mapped in transcripts\tZ Read fragments mapped in genes, a read touching 2 exons counts 2\t% Reads mapped in transcripts\t% Z Read fragments mapped in genes, a read touching 2 exons counts 2\t% difference", 5, 0, 0} , 
    { "Trait", "Stranding\tZ Stranding\t% difference", 1, 0, 0} ,
    { "Trait", "Significantly expressed genes\tTouched Genes\tZ expressed genes : is it touched or significant ?\tDifference", 2, 0, 0} ,
    { "Trait", "Sex\tZ sex", 10, 0, 0} ,

    { "Trait", "Stage\tZ dev_stage", 20, 0, 0} ,

    { "Trait", "Tissue-System", 30, 0, 0} ,

    { "Z_tissue", "Z tissue", 0, 0, 0} ,
    { "Z_cell_type", "Z cell_type", 0, 0, 0} ,
    { "Z_sample_type", "Z sample_type", 0, 0, 0} ,

    { "Trait", "Stage summary", 40, 0, 0} ,
    { "Z_genotype", "Z genotype", 0, 0, 0} ,
    { "Title", "Magic sample title", 0, 0, 0} ,
    { "Sequencing_length", "Sequencing length", 0, 0, 0} ,

    { "N_Stage", "NLM Stage", 0, 0, 0} ,
    { "N_CellType_AnatStage", "NLM CellType", 0, 0, 0} ,
    { "N_Strain_Treatment_Note", "NLM Strain_Treatment_Note", 0, 0, 0} ,


    {  0, 0, 0, 0, 0}
  }; 
  const char *caption =
    "Comparison with Zhen Xia"
    ;
  if (rc == (void *) 1)
    return  qcChapterCaption (qc, tts, caption) ;

  tta = rc ? ac_tag_table (rc->ali, "nh_Ali", h) : 0 ;
  ttb = rc ? ac_tag_table (rc->ali, "h_Ali", h) : 0 ;
  if (rc)
    {
      zRaw = ac_tag_float (rc->ali, "raw_data", 0) ;
      if (! zRaw)
	zRaw = ac_tag_float (rc->run, "Z_Total_reads", 0) ;
    }	      
  for (ti = tts ; ti->tag ; ti++)
    {	
      if (rc == 0)
	{
	  aceOutf (qc->ao, "\t%s", ti->title) ;
	}
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (qc->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{  
	  int ir ;
	  float z1 = 0, z2 = 0 ;
	  switch (ti->col)
	    {
	    case 1:  /* Total */ 
              z1 = ac_tag_float (rc->ali, "raw_data", 0) ;
	      z2 = ac_tag_float (rc->run, "Z_Total_reads", 0) ;	      
	      break ;
	    case 2:
	      for (ir = 0 ; tta ? ir < tta->rows : 0 ; ir++)
		{
		  if (! strcmp (ac_table_printable (tta, ir, 0, ""), "any"))
		    z1 = ac_table_float (tta, ir, 3, 0) ;
		}
	      z2 = ac_tag_float (rc->run, "Z_unique", 0) ;	      
	      z2 += ac_tag_float (rc->run, "Z_non_unique", 0) ;	      
	      break ;
	    case 3:
	      tt = ac_tag_table (rc->ali, "Unicity", h) ;
 	      for (ir = 0 ; tt ? ir < tt->rows : 0 ; ir++)
		{
		  if (! strcmp (ac_table_printable (tt, ir, 0, ""), "any"))
		    z1 = ac_table_float (tt, ir, 1, 0) ;
		}
	      z2 = ac_tag_float (rc->run, "Z_unique", 0) ;	      
	      break ;
	    case 30:
	      for (ir = 0 ; tta ? ir < tta->rows : 0 ; ir++)
		{
		  if (! strcmp (ac_table_printable (tta, ir, 0, ""), "any"))
		    z1 = ac_table_float (tta, ir, 3, 0) ;
		}
	      tt = ac_tag_table (rc->ali, "Unicity", h) ;
 	      for (ir = 0 ; tt ? ir < tt->rows : 0 ; ir++)
		{
		  if (! strcmp (ac_table_printable (tt, ir, 0, ""), "any"))
		    z2 = ac_table_float (tt, ir, 1, 0) ;
		}
	      z1 = z1 - z2 ;   /* mapped - unique */
	      z2 = ac_tag_float (rc->run, "Z_non_unique", 0) ;	      
	      break ;
	    case 31:
	      z2 = ac_tag_float (rc->ali, "raw_data", 0) ;
	      for (ir = 0 ; tta ? ir < tta->rows : 0 ; ir++)
		{
		  if (! strcmp (ac_table_printable (tta, ir, 0, ""), "any"))
		    z1 = ac_table_float (tta, ir, 3, 0) ;
		}
	  
	      z1 = z2 - z1 ;
	      z2 = ac_tag_float (rc->run, "Z_unmapped_reads", 0) ;	      
	      break ;
	    case 4:
	      tt = ac_tag_table (rc->ali, "Known_introns", h) ;
              z1 = tt ? ac_table_float (tt, 0, 2, 0) : 0 ;
	      z2 = ac_tag_float (rc->run, "Z_spliced", 0) ;	      
	      break ;
	    case 5:
	      for (ir = 0 ; ttb ? ir < ttb->rows : 0 ; ir++)
		{	
		  const char **tpp, *targets[] = { "0_SpikeIn", "1_DNASpikeIn", "A_mito", "B_rrna", "MT_EBI", "KT_RefSeq", "ET_av", 0 } ;

		  for (tpp = targets ; *tpp ; tpp++)
		    if (! strcmp (ac_table_printable (ttb, ir, 0, ""), *tpp)) 
		      { z1 += ac_table_float (ttb, ir, 3, 0) ; break ; }
		}
	      tt = ac_tag_table (rc->run, "Z_sponge", h) ;
	      for (ir = 0 ; tt ? ir < tt->rows : 0 ; ir++)
		{
		  if (! strcmp (ac_table_printable (tt, ir, 0, ""), "genes"))
		    z2 = ac_table_float (tt, ir, 1, 0) ;
		}
	      break ;
	    case 6:
              z1 = ac_tag_int (rc->ali, "Genes_with_index", 0) ;
	      z2 = ac_tag_int (rc->run, "Z_expressed_genes", 0) ;	      
	      
	      break ;
	    case 7:
	      break ;

	    }
	  aceOutf (qc->ao, "\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f"
		   , z1, z2
		   , 100.0 * z1/zRaw, 100.0*z2 / zRaw
		   , 100.0 * (z2 - z1) / zRaw
		   ) ;
	}
      else if (! strcmp (ti->tag, "Trait"))
	{  
	  int ir, jr, k ;
	  float z1, z2 ;
	  int g1, gt, g2 ;
	  switch (ti->col)
	    {
	    case 1:
	      tt = ac_tag_table (rc->run, "Observed_strandedness_in_ns_mapping", h) ;
              z1 = tt ? ac_table_float (tt, 0, 1, 0) : 0 ;
	      z2 = ac_tag_float (rc->run, "Z_strand", 0) ;
	      aceOutf (qc->ao, "\t%.2f\t%.2f\t%.2f", z1, z2, z1 - z2) ;
	      break ;
	    case 2:
              g1 = ac_tag_int (rc->ali, "Genes_with_index", 0) ;
              gt = ac_tag_int (rc->ali, "Genes_touched", 0) ;
	      g2 = ac_tag_int (rc->run, "Z_expressed_genes", 0) ;
	      aceOutf (qc->ao, "\t%d\t%d\t%d\t%d", g1, gt, g2, gt - g2) ;
	      break ;
	    case 10:  /* Sex */
	      aceOutf (qc->ao, "\t%s", ac_tag_printable (rc->sample, "Sex", "")) ;
	      aceOutf (qc->ao, "\t%s", ac_tag_printable (rc->run, "Z_sex", "")) ;
	      break ;
	    case 20:  /* dev stage */
	      aceOutf (qc->ao, "\t%s", ac_tag_printable (rc->sample, "Stage", "")) ;
	      aceOutf (qc->ao, "\t%s", ac_tag_printable (rc->run, "Z_dev_stage", "")) ;
	      break ;
	    case 30:  /* Tissue */ 
	      k = 0 ;
	      aceOutf (qc->ao, "\t") ;
	      tt = ac_tag_table (rc->sample, "Tissue", h) ;
	      for (ir = 0 ; tt && ir < tt->rows ; ir++)
		for (jr = 0 ; jr < tt->cols ; jr++)
		  {
		    const char *ccp = ac_table_printable (tt, ir, jr, 0) ;
		    if (ccp)
		      aceOutf (qc->ao, "%s%s"
			       , k++ ? "; " : ""
			       , ccp
			       ) ;
		  }
	      tt = ac_tag_table (rc->sample, "System", h) ;
	      for (ir = 0 ; tt && ir < tt->rows ; ir++)
		for (jr = 0 ; jr < tt->cols ; jr++)
		  {
		    const char *ccp = ac_table_printable (tt, ir, jr, 0) ;
		    if (ccp)
		      aceOutf (qc->ao, "%s%s"
			       , k++ ? "; " : ""
			       , ccp
			       ) ;
		  }

	      break ;
	    case 40:  /* stage summary */
	      k = 0 ;
	      aceOutf (qc->ao, "\t") ;
	      tt = ac_tag_table (rc->sample, "Stage", h) ;
	      for (ir = 0 ; tt && ir < tt->rows ; ir++)
		for (jr = 1 ; jr < tt->cols ; jr++)
		  {
		    const char *ccp = ac_table_printable (tt, ir, jr, 0) ;
		    if (ccp)
		      aceOutf (qc->ao, "%s%s"
			       , k++ ? "; " : ""
			       , ccp
			       ) ;
		  }
	      break ;

	    }
	}
      else
	qcShowTag (qc, rc, ti) ;
    }

   ac_free (h) ;
   return;
}  /* qcDrosoZhenXia */
#endif /* JUNK */

/*************************************************************************************/
/* i: qcMicroRNA not activated */
static const char *allMethods = "R0S1234HE5UWV9zKBFTPbCifmpAdMsr7DgIX8Lv" ;
static MM methods [] = {
  {'R', &qcRunSelection} ,
  {'0', &qcProject} ,
  {'S', &qcSortingTitles} ,
  {'1', &qcMainResults1} ,
  {'2', &qcMainResults2} ,
  {'3', &qcMainResults3} ,
  {'4', &qcMainResults4} ,
  {'H', &qcHighGenes} ,
  {'E', &qcGeneExpression1} ,
  {'5', &qcMainResults5} ,
  {'U', &qcIntergenic2} ,
  {'W', &qcMappingPerCentTarget} ,
  {'V', &qcHighVirusBacteria} ,
  {'9', &qcPairPercent } ,
  {'z', &qcInsertSize} ,
  {'K', &qcSnpCoding } ,
  {'B', &qcBloom} ,  /* adaptors  telomeres polyA */
  {'F', &qcReadFate1 } ,

  {'T', &qcTitles} ,
  {'P', &qcProtocol } ,
  {'b', &qcBeforeAli } ,


  {'C', &qcCPU} ,
  {'i', &qcMicroRNA} ,
  {'f', &qcReadFate2 } ,
  {'m', &qcMappingPerTarget } ,

  {'p', &qcPairNb } ,
  {'Z', &qcInsertSizeHisto} ,
  {'A', &qcAvLengthAli } ,
  {'d', &qcStrandedness } ,
  {'M', &qcMismatchTypes } ,

  {'s', &qcSnpTypes } ,
  {'r', &qcSnpRejectedTypes } ,
  {'7', &qc3pBias} ,
  {'D', &qcInterGenic} ,
  {'g', &qcGeneExpression2} ,
  {'I', &qcCandidateIntrons} ,
  {'X', &qcSexTissueSignatures} ,
  /* {'_', &qcDrosoZhenXia} , */
  {'8', &qcExternalFiles} ,
  {'L', &qcSLs} ,
  {'v', &qcCapture} ,
  {'o', &qcOther} ,
  {'a', &qcAli } ,   /* now included in R */ 
{ 0, 0 }
} ;

/*************************************************************************************/

static BOOL qcExportRun (QC *qc, int irc, int type)
{       
  const char *ccp, *lineName, *sTitle = 0 ;
  static char oldTitle[64] ; 
  MM *mm ;
  RC *rc = 0 ;
  char *cp ;
  KEY srp1, srp2 ;
  
  switch (type)
    {
    case 0: 
      oldTitle[0] = 0 ;
      lineName =  "### Caption" ; 
      memset (oldTitle, 0, sizeof (oldTitle)) ;
      rc = (void *) 1 ;
      break ;
    case 1: 
      lineName =  "### Run" ;
      rc = (void *) 0 ;
      break ;
    case 2: 
      rc = arrp (qc->runs, irc, RC) ;
      if (irc && rc->run && rc[-1].run &&
	  (srp1 = ac_tag_key (rc->run, "SRP", 0)) &&
	  (srp2 = ac_tag_key (rc[-1].run, "SRP", 0)) &&
	  srp1 != srp2
	  )
	aceOut (qc->ao, "\n") ;
      lineName = ac_name (rc->run) ;
      qcSetAli (qc, rc) ;

      sTitle = ac_tag_printable (rc->run, "Sorting_title", 0) ;
      cp = sTitle ? strstr (sTitle, "__") : 0 ;
      if (cp)
	{
	  char cc = *cp ;
	  *cp = 0 ;
	  if (*oldTitle && strcmp(sTitle, oldTitle))
	    aceOutf (qc->ao, "\n" ) ;
	  strncpy (oldTitle, sTitle, 64) ;
	  *cp = cc ;
	}
      break ; 
    }
   
  aceOutf (qc->ao, "%s", lineName) ;

  ccp = qc->export - 1 ;
  while (*++ccp)
    {
      mm = methods - 1 ;
      while (mm++, mm->cc)
	if (mm->cc == *ccp)
	  mm->f(qc, rc) ;
    }

  if (type == 2)
    ac_free (rc->h) ;

  aceOutf (qc->ao, "\n" ) ;
  return TRUE ;
} /* qcExportRun */

/*************************************************************************************/
/*************************************************************************************/

static int qcGetEtargets (QC *qc)
{
  AC_HANDLE h = ac_new_handle () ; 
  const char *ccp, *ccq ;
  int n = 0 ;
  
  memset (qc->Etargets, 0, sizeof (qc->Etargets)) ;
  ccp = getenv ("NO_INTRON") ;
  if (ccp && *ccp == '1') qc->noIntron = TRUE ;
  qc->noIntron = FALSE ;
  ccp = getenv ("Etargets") ;
  fprintf (stderr, "etargets = %s", ccp ? ccp : "NA") ; 
  if (! ccp || ! *ccp) ccp = "av" ;
  if (ccp && *ccp)
    {
      ACEIN ai = aceInCreateFromText (ccp, 0, h) ;
      if (aceInCard (ai))
	while ((ccq = aceInWord (ai)) && n < 3)
	  {
	    qc->Etargets[n] = strnew (ccq, qc->h) ;
	    if (! strcmp (ccq, "av"))
	      ccq = "AceView" ;
	    if (! strcmp (ccq, "transposon"))
	      ccq = "Repbase transposon" ;
	    qc->EtargetsBeau[n] = strnew (ccq, qc->h) ;
	    n++ ;
	  }
      qc->Etargets[n] = 0 ;
    }
  for ( ; n < 3 ; n++)
    {
      qc->Etargets[n] = hprintf (qc->h, "annotation %d", n + 1) ;
      qc->EtargetsBeau[n] = hprintf (qc->h, "annotation %d", n + 1) ;
    }
  ac_free (h) ;
  return n ;
}  /* qcGetEtargets */

/*************************************************************************************/

static int qcGetRuns (QC *qc)
{
  AC_OBJ run = 0 ;
  int nr = 0 ;
  
  qc->runs = arrayHandleCreate (256, RC, qc->h) ;
  
  if (qc->project)
    {
      int pass ;

      for (pass = 1 ; pass < 4 ; pass++)
	{
	  char buf[256], *cprj ;
	  strncpy (buf, qc->project, 255) ;
	  cprj = buf ;
	  if (qc->runType && pass != qc->runType)
	    continue ;

	  while (cprj && cprj[0])
	    {
	      const char *errors = 0 ;
	      int ir ;
	      AC_TABLE tt = 0 ;
	      AC_HANDLE h = ac_new_handle () ; 
	      char *qq, *qq2, *cqrj = strchr (cprj, ',') ;
	      
	      switch (pass)
		{
		case 1: /* sub lib */
		  qq2 = "r#Sublibrary_of" ;
		  break ;
		case 2: /* Runs */
		  qq2 = "(NOT r#Union_of AND NOT r#Sublibrary_of)" ;
		  break ;
		case 3: /* groups */
		  qq2 = "r#Union_of" ;
		  break ;
		}
		  
	      if (cqrj) {*cqrj = 0 ; cqrj++ ; }
	      qq = hprintf (h," select t,r  from p in class project where p == \"%s\", r in p->run where (NOT r#Private AND (r#Is_group OR r#IS_run OR r#Sublibraries) AND %s), t in r->%s", cprj, qq2, qc->orderBy) ;
	      
	      tt = ac_bql_table (qc->db, qq, 0, 0, &errors, h) ;
	      for (ir = 0 ; tt && ir < (tt->rows) ; ir++)
		{
		  AC_HANDLE h1 = ac_new_handle () ;
		  AC_OBJ run2 = 0 ;
		  run = ac_table_obj (tt, ir, 1, h1) ;
		 
		  if (run)
		    {
		      RC *rc = arrayp (qc->runs, nr++, RC) ;
		      rc->run = run ;
		      rc->h = h1 ;
		      run2 = ac_objquery_obj (run, "COUNT {>ali; SLS == SL*} > 0", h1) ;
		      if (run2 )
			qc->hasSL = TRUE ;		      
		      run2 = ac_objquery_obj (run, "COUNT {>ali; Capture} > 0", h1) ;
		      if (run2 )
			qc->hasCapture = TRUE ;		      
		    }
		  else
		    ac_free (h1) ;
		}
	      ac_free (h) ;
	      cprj = cqrj ;
	    }
	}
    }
 
  fprintf (stderr, " qcGetRuns got %d runs in project %s\n", nr, qc->project) ;
     
  return nr ;
} /* qcGetRuns */

/*************************************************************************************/

static void qcExportHeader (QC *qc)
{
  aceOutf (qc->ao, "## Project %s, quality control report, File %s : %s\n"
	   , qc->project
	   , aceOutFileName (qc->ao)
	   , timeShowNow () 
	   ) ;

  return ;
} /* qcExportHeader */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// qcsummary: complete Quality Control summary for the Magic pipeline\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, October 2014, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "// Given a project name\n"
	    "//    import all data from MetaDB, export a single QC table\n"
	    "//    one line per run, groups of columns for each export option\n"
	    "//    a smaller table can be obtained using -O and -run options\n"
	    "//\n"
	    "// Syntax:\n"
	    "// qcsummary -project <project> [options]\n"
	    "//   projectName MANDATORY export qc for all runs in this project (usually $MAGIC)\n" 
	    "// All other parameters are optional and may be specified in any order\n"
	    "//   -db <dir>: MetaDB directory (default MetaDB)\n" 
	    "//   -o output_file_prefix\n"
	    "//      all exported files will be named output_file_prefix.action\n" 
            "//   -gzo : gzip all output files\n"
	    "//   -TT  : limit the stats to the (run dependent) targeted genes\n"
	    "//   -export [TPbafpmg...] : only export some groups of columns, in the requested order\n"
	    "//      default: if -export is not specified, export all columns in default order\n"
	    "//            R:  Run selection\n"
	    "//            0:  Project Info\n"
	    "//            S:  Sorting titles\n"
	    "//            1:  Alignments read counts\n"
	    "//            2:  Alignments percentages\n"
	    "//            3:  Alignments per strand in genes\n"
	    "//            4:  Alignments statistics in Megabases\n"
	    "//            H:  High genes\n"
	    "//            E:  GeneExpression\n"
	    "//            5:  main Results5\n"
	    "//            U:  Intergenic2\n"
	    "//            W:  MappingPerTarget\n"
	    "//            V:  High Virus and bacteria\n"
	    "//            9:  Pair percent\n"
	    "//            z:  Insert size\n"
	    "//            K:  SNP coding\n"
	    "//            B:  Adators, telomeres, poly-A motifs\n"
	    "//            F:  ReadFate\n" 
	    "//            T:  Titles\n"
	    "//            P:  Protocol\n"
	    "//            b:  BeforeAli\n"
	    "//            C: CPU and RAM usage\n"
	    "//            i: Micro-RNA\n"
	    "//            f:  ReadFate2\n" 
	    "//            m:  Mapping per target\n"
	    "//            a:  Ali\n"
	    "//            p:  Pair\n" 
	    "//            A:  Average aligned length\n"
	    "//            d : strandedness\n"
	    "//            M: Mismatch types\n"
	    "//            s:  SNP types\n"
	    "//            r:  Rejected SNP types\n"
	    "//            7: 3p bias\n"
	    "//            D:  Intergenic\n"
	    "//            g: Gene expression 2\n"
	    "//            I:  Candidate introns\n"
	    "//            X:  Sex and tissue signatures\n"
	    "//            U:  Intergenic2\n"
	    "//            L:  SL1-12\n"
	    "//            Z: Zhenxi\n"
	    "//            v: capture\n"
	    "//            o:  Other\n"
	    "//   -orderBy <tag> : sort the lines by this tag\n"
	    "//      default: sorting_title\n"
	    "//      Only export in that order the runs belonging to the project, present in this list\n"
	    "//      Insert a blank line if the sorting value looks like A__B and A changes\n"
	    "//   -runType [0,1,2,3] : which runs\n"
	    "//       0 [default]: any run belonging to the project\n"
	    "//       1 : just the sublibraries\n"
	    "//       2 : just the runs (excluding type 1 and 3 )\n"
	    "//       3: just the groups\n"
	    "// Caveat:\n"
	    "//   Caption lines at the top start with a #\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  qcsummary -help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  QC qc ;
  AC_HANDLE h = 0 ;
  int nr ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&qc, 0, sizeof (QC)) ;
  qc.h = h ;
  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;
  qc.dbName = "MetaDB" ;
  getCmdLineOption (&argc, argv, "-db", &qc.dbName) ; 
  qc.TT = getCmdLineBool (&argc, argv, "-TT") ;
  qc.gzo = getCmdLineBool (&argc, argv, "-gzo") ;

  qc.export =  allMethods ;
  getCmdLineOption (&argc, argv, "-e", &qc.export) ;
  getCmdLineOption (&argc, argv, "-export", &qc.export) ;
  getCmdLineOption (&argc, argv, "--export", &qc.export) ;
  getCmdLineOption (&argc, argv, "-p", &qc.project) ;
  getCmdLineOption (&argc, argv, "-project", &qc.project) ;
  getCmdLineOption (&argc, argv, "--project", &qc.project) ;
  getCmdLineInt (&argc, argv, "-t", &qc.runType) ;
  getCmdLineInt (&argc, argv, "-runType", &qc.runType) ;
  getCmdLineOption (&argc, argv, "-f", &qc.externalFiles) ;
  if (! qc.project)
    {
      fprintf (stderr, "Sorry, missing parameter -project, please try qcsummary -help\n") ;
      exit (1) ;
    }
  if (! qc.dbName)
    {
      fprintf (stderr, "Sorry, missing parameter -db, please try qcsummary -help\n") ;
      exit (1) ;
    }
  if (qc.dbName)
    {
      const char *errors ;

      qc.db = ac_open_db (qc.dbName, &errors);
      if (! qc.db)
	{
	  fprintf (stderr, "Failed to open db %s, error %s", qc.dbName, errors) ;
	  exit (1) ;
	}
    }
  getCmdLineOption (&argc, argv, "-project", &qc.project) ;
  qc.orderBy = "Sorting_title" ;
  getCmdLineOption (&argc, argv, "-orderBy", &qc.orderBy) ;
  
  getCmdLineOption (&argc, argv, "-o", &qc.outFileName) ; 
  qc.ao = aceOutCreate (qc.outFileName, ".QC.txt", qc.gzo, h) ;
  if (! qc.ao)
    messcrash ("cannot create output file %s\n", qc.outFileName ? qc.outFileName : EMPTY ) ;


  if (argc > 1) usage (messprintf ("Unknown parameters %s", argv[1])) ;

  qcGetEtargets (&qc) ;

  qcExportHeader (&qc) ;

  nr = qcGetRuns (&qc) ;
  qcExportRun (&qc, 0, 0) ; /* chapter captions */
  qcExportRun (&qc, 0, 1) ; /* column titles */

  if (nr)
    {
      int irc ;
      RC *rc ;

      for (irc = 0, rc = arrp (qc.runs, 0, RC) ; irc < nr ; irc++, rc++)
	qcExportRun (&qc, irc, 2) ;
    }
  
  ac_free (h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
} /* main */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

